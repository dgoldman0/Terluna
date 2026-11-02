"""Realise a phrase with whole-body dynamics: the score's units become an optimal control problem in Crocoddyl.

    python -m humanities.dance.motion.realise PHRASE [--units i:j] [--span t0:t1] [--windows] [--finish S]
    # needs the dance environment

The grammar fixes what happens and when: which feet bear weight in each phase, how long each flight lasts, which
shapes the body passes through, how fast it travels and turns. This module turns that into a multiple-shooting
problem over the dancer's full rigid-body dynamics at lunar gravity (body.json through pinmodel.py):

- a foot flat on the floor or on demi-pointe holds it in 6D (the flat sole, or the toes' sole), with a wrench-cone
  cost that keeps every contact force inside the friction pyramid, its centre of pressure inside the foot and its yaw
  moment within what friction allows;
- a foot pivoting in a turn holds the floor at one point under the ball, which cannot slip but can swivel, and its two
  metatarsal heads press on the floor without friction, so the foot cannot roll and the floor's push can move across
  the ball for balance as it does on demi-pointe. The floor's turning friction on the ball (the envelope's spin friction
  over its 2 cm effective radius, at the dancer's full lunar weight) acts as a moment against the swivel, so a whirl
  slows unless the other foot pushes;
- a kneeling knee holds the floor at a point on the front of the shank;
- each flight node is free fall of the articulated body, so angular momentum is conserved and the centre of mass
  follows its ballistic arc exactly;
- each landing (a foot from the air, a heel coming down from demi-pointe, a pivot planting, a knee touching down)
  passes an impulse model that stops the new contact, as a rigid contact does;
- joint torques stay inside the body's strength (box-constrained FDDP), joint ranges inside barrier costs, every point
  of the soles and knees that bears no weight at or above the floor;
- the shapes enter as posture targets interpolated through each phase and the turns as the pelvis's yaw. The travel
  enters as a centre-of-mass reference: a fall's follows the envelope's toppling pendulum and a weight shift's moves
  across to the standing foot. The footholds come from the plan: a foot landing from the air goes under the centre of
  mass at mid-stance, a step or a kneel where its shape puts it, a turn's preparation in fourth.

The solution is a torque, force and motion history that satisfies the equations of motion at every node; the replay
(replay.py) runs those torques through MuJoCo's independent contact dynamics as a check. Air drag on the body is left
out of the optimisation: the envelope shows it changes a jump's height by under 4% even with large sleeves. A long
phrase is solved in overlapping windows (solve_windows).
"""
from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pinocchio as pin
import crocoddyl

from humanities.dance import envelope as env
from humanities.dance import grammar as gr
from humanities.dance.motion.pinmodel import build, neutral
from humanities.dance.motion_io import body_spec, save_motion

HERE = Path(__file__).resolve().parent
RUNS = Path('/media/projectspace/terluna-research/humanities/dance-runs/motion')
DT = 0.025
MU = gr.MU
FOOT_HALF = dict(flat=(0.085, 0.04), ball=(0.035, 0.035))     # centre-of-pressure box: half length, half width (m)
STANCE_WIDTH = 0.10
SOLE_POINTS = ('heel_m', 'heel_l', 'ball_m', 'ball_l', 'toe')
FLOOR_POINTS = SOLE_POINTS + ('kneecap',)
FLOOR_MARGIN_M = 0.002
MAX_FOOT_FORCE_N = 900.0     # a foot's peak push, about 10 lunar body weights (1.7 of Earth's): a jump's take-off peak
# Contact kinds: the frame that holds the floor, how many directions it holds (6: no slip and no swivel; 3: a point
# that cannot slip but may swivel), and the body's points that rest on the floor with it.
CONTACTS = {
    'flat': ('sole', 6, ('heel_m', 'heel_l', 'ball_m', 'ball_l')),
    'ball': ('toesole', 6, ('ball_m', 'ball_l', 'toe')),
    'pivot': ('toesole', 3, ('ball_m', 'ball_l', 'toe')),
    'knee': ('kneecap', 3, ('kneecap',)),
}
PIVOT_MOMENT_NM = env.SPIN_FRICTION * env.MASS_KG * gr.G * env.PIVOT_RADIUS_M   # the ball's turning friction at full weight
# The swivel rate over which that friction builds to its full value (a smoothed Coulomb law). Below it the law acts as
# a damper of slope M / w0 on the leg's own twist, whose inertia is a few hundredths of a kg m^2; the explicit DT step
# stays stable only while DT M / (w0 I) stays under about 1, so w0 is 1 rad/s, where tanh still gives 99% of the full
# friction at the whirl's 2.8 rad/s.
PIVOT_SLIP_RAD_S = 1.0
PUSH_REV = 1.0 / 6           # a whirl's push lasts a sixth of a revolution
FOURTH_M = (-0.30, 0.25)     # the free foot behind and beside the turning foot in a turn's preparation: 0.39 m apart
RELEVE = {'ankle': 45.0, 'mtp': 65.0}   # the foot on demi-pointe (poses.json, releve)
PIVOT_SUPPORTS = ('ball_m', 'ball_l')   # a pivot's metatarsal heads: they press on the floor without friction
SPEED_CHANGE_M_S2 = 0.5 * MU * gr.G     # how fast a unit without flights changes its travel: half the floor's grip


# --- shapes to joint targets ---------------------------------------------------------------------------------------

def quat_axis_angle(axis, angle):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    return pin.Quaternion(math.cos(angle / 2), *(math.sin(angle / 2) * a))


def rz(yaw):
    return pin.utils.rpyToMatrix(0, 0, yaw)


def ball_target(spec, side_sign):
    """A ball joint's rotation from elev/az/twist: swing the long axis (down at rest) toward the direction at
    elevation `elev` and azimuth `az` (outward is the limb's own side), then twist about the long axis."""
    elev = math.radians(spec.get('elev', 0.0))
    az = math.radians(spec.get('az', 0.0))
    twist = math.radians(spec.get('twist', 0.0))
    d = np.array([math.cos(az), side_sign * math.sin(az), 0.0])        # horizontal direction of the swing
    axis = np.cross([0.0, 0.0, -1.0], d)
    q_swing = quat_axis_angle(axis, elev) if np.linalg.norm(axis) > 1e-9 else pin.Quaternion.Identity()
    q_twist = quat_axis_angle([0.0, 0.0, side_sign], twist)
    return q_swing * q_twist


SIDED_SIGN = {'lumbar_bend', 'lumbar_rot', 'thoracic_bend', 'thoracic_rot', 'neck_bend', 'neck_rot'}


def shape_targets(shape, side):
    """Joint targets {joint: angle or quaternion} for one key shape (a dict of region -> pose name) with gesture side
    `side`. Mirrored shapes swap G and S."""
    G = side
    S = gr.other(side)
    if shape.get('mirror'):
        G, S = S, G
    out = {}
    for region in ('legs', 'arms', 'torso', 'head', 'feet'):
        name = shape.get(region)
        if not name:
            continue
        joints = gr.POSES['regions'][region][name]['joints']
        for jname, val in joints.items():
            absolute = jname.replace('G_', f'{G}_').replace('S_', f'{S}_')
            limb_side = absolute[0] if absolute[:2] in ('L_', 'R_') else None
            if isinstance(val, dict):
                out[absolute] = ('ball', ball_target(val, 1 if limb_side == 'L' else -1))
            else:
                ang = math.radians(val)
                if jname in SIDED_SIGN and G == 'R':
                    ang = -ang        # positive bends and turns toward G; the model's positive sense is to the left
                out[absolute] = ('hinge', ang)
    return out


CARRIAGE = {   # the dancer's carriage between shapes: lifted, arms low and rounded, slight turnout
    'legs': 'turned_out', 'arms': 'bras_bas', 'torso': 'upright', 'head': 'forward'}
TWO_FOOT_LEGS = {'parallel', 'turned_out', 'demi_plie', 'grand_plie_prep', 'releve', 'lunge_forward', 'kneel_one'}


def place_feet(model, data, q, targets, iters=80, damping=1e-6, com=None):
    """Inverse kinematics that puts each contact frame at its target, moving the root and the joints by the least
    change: a placement (pin.SE3) for a foot held in 6D, a point for a pivot or a knee; with com, the centre of mass
    also goes over that point of the floor."""
    q = q.copy()
    for _ in range(iters if targets else 0):
        pin.framesForwardKinematics(model, data, q)
        errs, Js = [], []
        for fid, M in targets.items():
            if isinstance(M, pin.SE3):
                errs.append(pin.log6(data.oMf[fid].actInv(M)).vector)
                Js.append(pin.computeFrameJacobian(model, data, q, fid, pin.LOCAL))
            else:
                errs.append(np.asarray(M, float) - data.oMf[fid].translation)
                Js.append(pin.computeFrameJacobian(model, data, q, fid, pin.LOCAL_WORLD_ALIGNED)[:3])
        if com is not None:
            Jc = pin.jacobianCenterOfMass(model, data, q)
            errs.append(np.asarray(com, float)[:2] - data.com[0][:2])
            Js.append(Jc[:2])
        e = np.concatenate(errs)
        if np.linalg.norm(e) < 1e-7:
            break
        J = np.vstack(Js)
        dq = J.T @ np.linalg.solve(J @ J.T + damping * np.eye(J.shape[0]), e)
        q = pin.integrate(model, q, dq)
    return q


# --- the floor's turning friction on a pivot -------------------------------------------------------------------------

class PivotFriction(crocoddyl.ActuationModelAbstract):
    """The floating base's actuation plus the floor's turning friction on each pivoting foot: a moment about the
    vertical against the foot's swivel, -M tanh(w / w0), with M the envelope's spin friction at the dancer's full lunar
    weight. The moment acts through the foot's Jacobian on every joint up to the root. Its change with the foot's
    attitude is left out of the derivatives (second order for a level foot); its dependence on the swivel rate is
    exact."""

    def __init__(self, state, frames, moment=PIVOT_MOMENT_NM, w0=PIVOT_SLIP_RAD_S):
        crocoddyl.ActuationModelAbstract.__init__(self, state, state.nv - 6)
        self.frames, self.moment, self.w0 = list(frames), moment, w0
        self.pm = state.pinocchio
        self.pd = self.pm.createData()

    def friction(self, x):
        """The generalised force of the turning friction and its derivative with respect to the velocities."""
        nq, nv = self.pm.nq, self.state.nv
        q, v = x[:nq], x[nq:]
        tau = np.zeros(nv)
        dv = np.zeros((nv, nv))
        for fid in self.frames:
            jz = pin.computeFrameJacobian(self.pm, self.pd, q, fid, pin.LOCAL_WORLD_ALIGNED)[5]
            th = math.tanh(float(jz @ v) / self.w0)
            tau += jz * (-self.moment * th)
            dv += np.outer(jz, jz) * (-self.moment * (1.0 - th * th) / self.w0)
        return tau, dv

    def calc(self, data, x, u):
        tau, _ = self.friction(x)
        tau[6:] += u
        data.tau[:] = tau

    def calcDiff(self, data, x, u):
        _, dv = self.friction(x)
        nv = self.state.nv
        data.dtau_dx[:, :] = 0.0
        data.dtau_dx[:, nv:] = dv
        data.dtau_du[:, :] = 0.0
        data.dtau_du[6:, :] = np.eye(nv - 6)

    def commands(self, data, x, tau):
        f, _ = self.friction(x)
        data.u[:] = (tau - f)[6:]

    def torqueTransform(self, data, x, u):
        data.Mtau[:, :] = 0.0
        data.Mtau[:, 6:] = np.eye(self.state.nv - 6)

    def createData(self):
        return crocoddyl.ActuationDataAbstract(self)


# --- the timeline ----------------------------------------------------------------------------------------------------

@dataclass
class Node:
    t: float
    contacts: dict                 # side -> contact kind (CONTACTS)
    landing: list = field(default_factory=list)   # sides whose contact starts here (an impulse stops it)
    targets: dict = field(default_factory=dict)
    yaw: float = 0.0
    yaw_rate: float = 0.0
    com: np.ndarray | None = None
    com_vel: np.ndarray | None = None
    phase: str = ''
    unit: str = ''
    family: str = ''
    side: str = 'R'                # the unit's gesture side
    phase_id: int = 0


def contact_kind(phase, unit):
    shapes = phase.get('shapes') or []
    legs = [s.get('legs') for s in shapes if isinstance(s, dict)]
    if phase['kind'] == 'turn' or 'releve' in legs or unit in ('skim',):
        return 'ball'
    return 'flat'


def phase_contacts(pl, ph, frac, prev, ball_frac=0.0, revs=0.0):
    """Which sides hold the floor, and how, at fraction frac of a phase. A fall rolls onto the ball once the centre of
    mass passes over it; a turn on one foot pivots on its ball, the other foot touching down for its pushes; a kneel
    rests the gesture foot flat and the supporting knee on the floor. A phase may name its contacts ('contact')."""
    G, S = pl.side, gr.other(pl.side)
    sup = ph.get('support', 'same')
    named = {({'S': S, 'G': G}.get(k, k)): v for k, v in (ph.get('contact') or {}).items()}
    if ph['kind'] == 'flight' or sup == 'none':
        return {}
    if sup == 'kneel':
        return {G: 'flat', S: 'knee'}
    if sup == 'same':
        return dict(prev) if prev else {'L': 'flat', 'R': 'flat'}
    sides = ['L', 'R'] if sup in ('both', 'floor') else [sup]
    if ph['kind'] == 'fall':
        kind = 'flat' if frac < ball_frac else 'ball'
        return {s: named.get(s, kind) for s in sides}
    if ph['kind'] == 'turn' and len(sides) == 1:
        out = {sides[0]: 'pivot'}
        if ph.get('pushes'):
            r = frac * revs
            if any(k + 0.5 <= r < k + 0.5 + PUSH_REV for k in range(int(revs) + 1)):
                out[gr.other(sides[0])] = 'ball'
        return out
    kind = contact_kind(ph, pl.unit)
    return {s: named.get(s, kind) for s in sides}


def fall_lean(a0_deg, a1_deg):
    """A fall's lean (rad) against time from rest at a0, by the envelope's inverted pendulum about the ankle."""
    L = gr._com_height()
    h = 1e-3
    th, w, ts, ths = math.radians(a0_deg), 0.0, [0.0], [math.radians(a0_deg)]
    while th < math.radians(a1_deg):
        w += (gr.G / L) * math.sin(th) * h
        th += w * h
        ts.append(ts[-1] + h)
        ths.append(th)
    return lambda t: float(np.interp(t, ts, ths))


def fall_range(pl):
    """A fall unit's lean from and to (deg), read from its duration, topple_s(a0, a1)."""
    expr = next(p['duration'] for p in gr.UNITS[pl.unit]['phases'] if p['kind'] == 'fall')
    a0, a1 = expr.strip()[len('topple_s('):-1].split(',')
    return gr.evaluate(a0, pl.params), gr.evaluate(a1, pl.params)


def ball_ahead_m():
    """How far the ball of the foot (the metatarsal heads) lies ahead of the ankle."""
    return body_spec()['contacts']['L_ball_m']['pos'][0]


def timeline(placed, dt=DT, finish_s=0.0):
    """Nodes every dt from the placed units: contacts from each phase's support, shape targets interpolated between
    the phase's key shapes (held between phases), the centre of mass's planned path, yaw for turns. finish_s holds the
    last shape at rest at the end, so the phrase closes still."""
    phases = [(pl, ph) for pl in placed for ph in pl.phases]
    t_end = placed[-1].exit.t
    times = np.arange(placed[0].t0, t_end - 1e-9, dt)
    pos = np.array([placed[0].entry.x, placed[0].entry.y], float)
    nodes, last_shape, prev = [], {}, None
    yaw_turn, turn_sign = 0.0, 1.0
    k, cur, pid = 0, -1, -1
    L = gr._com_height()
    for t in times:
        while k < len(phases) - 1 and t >= phases[k][1]['t0'] + phases[k][1]['duration'] - 1e-9:
            k += 1
        pl, ph = phases[k]
        u = gr.UNITS[pl.unit]
        if k != cur:
            # A new phase: set up what its course needs.
            cur, pid = k, pid + 1
            heading0 = math.radians(pl.entry.heading)
            revs = 0.0
            if ph.get('turn_rev'):
                revs = gr.evaluate(ph['turn_rev'], pl.params) if isinstance(ph['turn_rev'], str) else float(ph['turn_rev'])
                if pl.entry.spin < 0.05:      # a new turn goes toward the side the unit names; a spin carries on
                    toward = pl.side if ph.get('toward', 'G') == 'G' else gr.other(pl.side)
                    turn_sign = 1.0 if toward == 'L' else -1.0
            ball_frac = 0.0
            if ph['kind'] == 'fall':
                a0, a1 = fall_range(pl)
                lean = fall_lean(a0, a1)
                origin, along = pos.copy(), np.array([math.cos(heading0), math.sin(heading0)])
                th_ball = math.degrees(math.asin(min(ball_ahead_m() / L, 1.0)))
                ball_frac = gr.topple_s(a0, th_ball) / ph['duration'] if th_ball > a0 else 0.0
            if ph['kind'] == 'shift':
                origin = pos.copy()
                sgn = 1.0 if gr.other(pl.side) == 'L' else -1.0
                across = sgn * pl.params.get('across_m', STANCE_WIDTH) * np.array([-math.sin(heading0), math.cos(heading0)])
            flights = any(p['kind'] == 'flight' for p in pl.phases)
        frac = (t - ph['t0']) / max(ph['duration'], 1e-9)
        ufrac = (t - pl.t0) / max(pl.duration, 1e-9)
        dh = ((pl.exit.heading - pl.entry.heading + 180) % 360) - 180
        heading = math.radians(pl.entry.heading + dh * ufrac)
        before = pos.copy()
        # The centre of mass's path.
        if ph['kind'] == 'fall':
            pos = origin + along * L * (math.sin(lean(t - ph['t0'])) - math.sin(lean(0.0)))
        elif ph['kind'] == 'shift':
            s = min(max(frac, 0.0), 1.0)
            pos = origin + across * s * s * (3 - 2 * s)
        else:
            if flights:
                speed = pl.entry.speed + (pl.exit.speed - pl.entry.speed) * ufrac
            else:
                dv = pl.exit.speed - pl.entry.speed
                speed = pl.entry.speed + math.copysign(min(abs(dv), SPEED_CHANGE_M_S2 * (t - pl.t0)), dv)
            pos = pos + dt * speed * np.array([math.cos(heading), math.sin(heading)])
        yaw_rate = 0.0
        if ph.get('turn_rev'):
            yaw_rate = turn_sign * 2 * math.pi * revs / ph['duration']
            yaw_turn += yaw_rate * dt
        contacts = phase_contacts(pl, ph, frac, prev, ball_frac, revs)
        landing = []
        if prev is not None:
            for s, kind in contacts.items():
                if s not in prev:
                    landing.append(s)                      # from the air
                elif prev[s] != kind and (kind in ('flat', 'knee') or (prev[s] == 'pivot' and kind == 'ball')):
                    landing.append(s)                      # a heel or a knee coming down, a pivot planting
        # Shape targets: the last key shape reached, blended toward the next within the phase.
        shapes = ph.get('shapes') if isinstance(ph.get('shapes'), list) else []
        target = dict(last_shape)
        if len(contacts) == 2 and target.get('legs') not in TWO_FOOT_LEGS:
            target.pop('legs', None)                       # both feet down: a flight's or gesture's legs end
        for s in shapes:
            if s['at'] <= frac + 1e-9:
                target = s
        nxt = next((s for s in shapes if s['at'] > frac + 1e-9), None)
        tg = shape_targets({**CARRIAGE, **{k2: v for k2, v in target.items() if k2 != 'at'}}, pl.side)
        if nxt is not None:
            prev_at = max([s['at'] for s in shapes if s['at'] <= frac + 1e-9], default=0.0)
            w = (frac - prev_at) / max(nxt['at'] - prev_at, 1e-9)
            w = w * w * (3 - 2 * w)
            tn = shape_targets({**CARRIAGE, **{k2: v for k2, v in nxt.items() if k2 != 'at'}}, pl.side)
            for j, (kind, val) in tn.items():
                if j in tg and kind == 'hinge':
                    tg[j] = ('hinge', (1 - w) * tg[j][1] + w * val)
                elif j in tg:
                    tg[j] = ('ball', tg[j][1].slerp(w, val))
                else:
                    tg[j] = (kind, val)
        for s, kind in contacts.items():
            if kind == 'pivot':                            # the turning foot stays on demi-pointe
                for j, deg in RELEVE.items():
                    tg[f'{s}_{j}'] = ('hinge', math.radians(deg))
        nodes.append(Node(t=float(t), contacts=contacts, landing=landing, targets=tg, yaw=heading + yaw_turn,
                          yaw_rate=yaw_rate, com=pos.copy(), com_vel=(pos - before) / dt, phase=ph['kind'],
                          unit=pl.unit, family=u['family'], side=pl.side, phase_id=pid))
        if shapes:
            reached = [s for s in shapes if s['at'] <= frac + 1e-9]
            if reached:
                last_shape = {k2: v for k2, v in reached[-1].items() if k2 != 'at'}
        prev = contacts
    last = nodes[-1]
    for i in range(int(round(finish_s / dt))):
        nodes.append(Node(t=last.t + (i + 1) * dt, contacts=dict(last.contacts), targets=dict(last.targets),
                          yaw=last.yaw, com=last.com.copy(), com_vel=np.zeros(2), phase='hold', unit=last.unit,
                          family=last.family, side=last.side, phase_id=pid + 1))
    return nodes


# --- the optimal control problem -----------------------------------------------------------------------------------------

class Problem:
    def __init__(self, nodes, gravity=gr.G):
        self.model, self.bounds = build(gravity)
        self.data = self.model.createData()
        self.state = crocoddyl.StateMultibody(self.model)
        self.actuation = crocoddyl.ActuationModelFloatingBase(self.state)
        self.nu = self.actuation.nu
        self.nodes = nodes
        self.fid = {(s, k): self.model.getFrameId(f'{s}_{v[0]}') for s in 'LR' for k, v in CONTACTS.items()}
        self.q0 = neutral(self.model)
        self.x_ref_cache = {}
        self.end_at_rest = bool(nodes) and nodes[-1].phase == 'hold'
        d = self.model.createData()
        pin.framesForwardKinematics(self.model, d, self.q0)
        # The toes' sole centre seen from the flat sole's centre at rest (the foot's frame): demi-pointe in place.
        self.toes_ahead = {s: d.oMf[self.fid[(s, 'flat')]].actInv(d.oMf[self.fid[(s, 'ball')]]).translation.copy()
                           for s in 'LR'}
        # How long each landing stays down, and where each phase ends, in this problem's nodes.
        self.phase_last = {}
        for i, n in enumerate(nodes):
            self.phase_last[n.phase_id] = i
        self.stay = []
        for i, n in enumerate(nodes):
            st = {}
            for s in (n.contacts if i == 0 else n.landing):
                j = i
                while j + 1 < len(nodes) and s in nodes[j + 1].contacts and s not in nodes[j + 1].landing:
                    j += 1
                st[s] = j
            self.stay.append(st)

    # A reference state for a node: the carriage with the node's targets applied, the root at the planned travel.
    def reference_q(self, node):
        q = self.q0.copy()
        for jname, (kind, val) in node.targets.items():
            if not self.model.existJointName(jname):
                continue
            jid = self.model.getJointId(jname)
            iq = self.model.idx_qs[jid]
            if kind == 'ball':
                q[iq:iq + 4] = val.coeffs()
            else:
                lo, hi = self.model.lowerPositionLimit[iq], self.model.upperPositionLimit[iq]
                q[iq] = float(np.clip(val, lo, hi))
        q[0:2] = node.com
        q[3:7] = pin.Quaternion(rz(node.yaw)).coeffs()
        return q

    def state_weights(self, node):
        nv = self.model.nv
        w = np.ones(2 * nv)
        w[:6] = 0.0                              # the root is free: dynamics and contacts place it
        w[6:nv] = 2.0
        w[nv:nv + 6] = 0.5
        w[nv + 6:] = 0.2
        return w

    def costs(self, node, landing_targets=None, terminal=False):
        cm = crocoddyl.CostModelSum(self.state, self.nu)
        x_ref = np.concatenate([self.reference_q(node), np.zeros(self.model.nv)])
        act = crocoddyl.ActivationModelWeightedQuad(self.state_weights(node) ** 2)
        cm.addCost('posture', crocoddyl.CostModelResidual(self.state, act, crocoddyl.ResidualModelState(self.state, x_ref, self.nu)), 1.0)
        if not terminal:
            cm.addCost('effort', crocoddyl.CostModelResidual(self.state, crocoddyl.ResidualModelControl(self.state, self.nu)), 1e-3)
            # Smoothness: joint accelerations kept small spread pushes and landings over their stances, as muscles do.
            cm.addCost('smooth', crocoddyl.CostModelResidual(self.state, crocoddyl.ResidualModelJointAcceleration(self.state, self.nu)), 2e-3)
        # Joint ranges as barriers on the hinge coordinates.
        bounds = crocoddyl.ActivationBounds(self._xlb(), self._xub(), 1.0)
        cm.addCost('ranges', crocoddyl.CostModelResidual(self.state, crocoddyl.ActivationModelQuadraticBarrier(bounds),
                                                         crocoddyl.ResidualModelState(self.state, self.nu)), 1e3)
        # Facing: the pelvis stays near upright and its yaw follows the travel. In a turn the yaw's rate follows the
        # plan instead of its angle, so a turn that lags never meets the half-revolution where an angle error wraps.
        pelvis = self.model.getFrameId('pelvis')
        turning = abs(node.yaw_rate) > 1e-9
        cm.addCost('facing', crocoddyl.CostModelResidual(
            self.state, crocoddyl.ActivationModelWeightedQuad(np.array([1.0, 1.0, 0.0 if turning else 4.0]) ** 2),
            crocoddyl.ResidualModelFrameRotation(self.state, pelvis, rz(node.yaw), self.nu)), 10.0)
        if turning:
            cm.addCost('spin', crocoddyl.CostModelResidual(
                self.state, crocoddyl.ActivationModelWeightedQuad(np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0])),
                crocoddyl.ResidualModelFrameVelocity(self.state, pelvis, pin.Motion(np.zeros(3), np.array([0.0, 0.0, node.yaw_rate])),
                                                     pin.LOCAL_WORLD_ALIGNED, self.nu)), 20.0)
        # Travel: the centre of mass follows the planned path in the horizontal while the body holds the floor.
        if node.contacts and node.com is not None:
            com_ref = np.array([node.com[0], node.com[1], 0.0])
            cm.addCost('travel', crocoddyl.CostModelResidual(
                self.state, crocoddyl.ActivationModelWeightedQuad(np.array([1.0, 1.0, 0.0])),
                crocoddyl.ResidualModelCoMPosition(self.state, com_ref, self.nu)), 20.0)
        for side, kind in node.contacts.items():
            fid = self.fid[(side, kind)]
            if CONTACTS[kind][1] == 6:
                # Forces inside the friction pyramid, centre of pressure inside the foot, yaw moment within friction.
                cone = crocoddyl.WrenchCone(np.eye(3), MU, np.array(FOOT_HALF[kind]) * 2, 4, True, 1.0, MAX_FOOT_FORCE_N)
                res = crocoddyl.ResidualModelContactWrenchCone(self.state, fid, cone, self.nu, True)
                b = crocoddyl.ActivationBounds(cone.lb, cone.ub)
                cm.addCost(f'{side}_cone', crocoddyl.CostModelResidual(self.state, crocoddyl.ActivationModelQuadraticBarrier(b), res), 1e2)
                # Light cost on sideways forces and moments: the feet do not squeeze against each other for nothing.
                cm.addCost(f'{side}_shear', crocoddyl.CostModelResidual(
                    self.state, crocoddyl.ActivationModelWeightedQuad(np.array([1.0, 1.0, 0.0, 30.0, 30.0, 30.0]) ** 2),
                    crocoddyl.ResidualModelContactForce(self.state, fid, pin.Force.Zero(), 6, self.nu, True)), 1e-3)
            else:
                # A point: its force inside the friction pyramid, pressing on the floor.
                cone = crocoddyl.FrictionCone(np.eye(3), MU, 4, True, 1.0, MAX_FOOT_FORCE_N)
                res = crocoddyl.ResidualModelContactFrictionCone(self.state, fid, cone, self.nu, True)
                b = crocoddyl.ActivationBounds(cone.lb, cone.ub)
                cm.addCost(f'{side}_cone', crocoddyl.CostModelResidual(self.state, crocoddyl.ActivationModelQuadraticBarrier(b), res), 1e2)
            if kind == 'pivot':
                # The metatarsal heads can only press on the floor.
                for pt in PIVOT_SUPPORTS:
                    res = crocoddyl.ResidualModelContactForce(self.state, self.model.getFrameId(f'{side}_{pt}'),
                                                              pin.Force.Zero(), 1, self.nu, True)
                    b = crocoddyl.ActivationBounds(np.array([1.0]), np.array([MAX_FOOT_FORCE_N]))
                    cm.addCost(f'{side}_{pt}_press', crocoddyl.CostModelResidual(
                        self.state, crocoddyl.ActivationModelQuadraticBarrier(b), res), 1e2)
        # The floor: every sole or knee point that bears no weight at this node stays at or above it (a 2 mm margin),
        # so swing feet, toes under a flat foot, heels over demi-pointe and knees never pass through the floor.
        for side in ('L', 'R'):
            kind = node.contacts.get(side)
            resting = CONTACTS[kind][2] if kind else ()
            for pt in FLOOR_POINTS:
                if pt in resting:
                    continue
                fid = self.model.getFrameId(f'{side}_{pt}')
                floor = crocoddyl.ActivationBounds(np.array([-1e6, -1e6, FLOOR_MARGIN_M]), np.array([1e6, 1e6, 1e6]))
                cm.addCost(f'{side}_{pt}_floor', crocoddyl.CostModelResidual(
                    self.state, crocoddyl.ActivationModelQuadraticBarrier(floor),
                    crocoddyl.ResidualModelFrameTranslation(self.state, fid, np.zeros(3), self.nu)), 1e4)
        # Landings: the arriving foot or knee reaches its planned foothold.
        for side, (fid, M) in (landing_targets or {}).items():
            if isinstance(M, pin.SE3):
                res = crocoddyl.ResidualModelFramePlacement(self.state, fid, M, self.nu)
            else:
                res = crocoddyl.ResidualModelFrameTranslation(self.state, fid, np.asarray(M, float), self.nu)
            cm.addCost(f'{side}_land', crocoddyl.CostModelResidual(self.state, res), 1e4)
        if terminal and self.end_at_rest:
            # The phrase ends at rest.
            w = np.zeros(2 * self.model.nv)
            w[self.model.nv:] = 1.0
            cm.addCost('rest', crocoddyl.CostModelResidual(
                self.state, crocoddyl.ActivationModelWeightedQuad(w),
                crocoddyl.ResidualModelState(self.state, np.concatenate([self.q0, np.zeros(self.model.nv)]), self.nu)), 1e2)
        return cm

    def _xlb(self):
        if 'xlb' not in self.x_ref_cache:
            nv = self.model.nv
            lb = np.full(2 * nv, -1e6)
            ub = np.full(2 * nv, 1e6)
            for jid in range(2, self.model.njoints):
                if self.model.nqs[jid] == 1:
                    iv = self.model.idx_vs[jid]
                    iq = self.model.idx_qs[jid]
                    lb[iv] = self.model.lowerPositionLimit[iq]
                    ub[iv] = self.model.upperPositionLimit[iq]
            self.x_ref_cache['xlb'], self.x_ref_cache['xub'] = lb, ub
        return self.x_ref_cache['xlb']

    def _xub(self):
        self._xlb()
        return self.x_ref_cache['xub']

    def actuation_for(self, node):
        """The floating base's actuation, with the floor's turning friction on any pivoting foot."""
        pivots = [self.fid[(s, k)] for s, k in node.contacts.items() if k == 'pivot']
        return PivotFriction(self.state, pivots) if pivots else self.actuation

    def contact_model(self, node, placements):
        cmult = crocoddyl.ContactModelMultiple(self.state, self.nu)
        gains = np.array([0.0, 30.0])
        for side, kind in node.contacts.items():
            fid = self.fid[(side, kind)]
            M = placements[(side, kind)]
            if CONTACTS[kind][1] == 6:
                c = crocoddyl.ContactModel6D(self.state, fid, M, pin.LOCAL_WORLD_ALIGNED, self.nu, gains)
            else:
                c = crocoddyl.ContactModel3D(self.state, fid, M.translation, pin.LOCAL_WORLD_ALIGNED, self.nu, gains)
            cmult.addContact(f'{side}_{kind}', c)
            if kind == 'pivot':
                for pt in PIVOT_SUPPORTS:
                    cmult.addContact(f'{side}_{pt}', crocoddyl.ContactModel1D(
                        self.state, self.model.getFrameId(f'{side}_{pt}'), 0.0, pin.LOCAL_WORLD_ALIGNED, np.eye(3), self.nu, gains))
        return cmult

    # Footholds -----------------------------------------------------------------------------------------------------

    def _shape_point(self, node, fid):
        """Where the node's shape puts a frame, the pelvis on the planned path, laid on the floor."""
        q = self.reference_q(node)
        q[2] = self.q0[2]
        pin.framesForwardKinematics(self.model, self.data, q)
        self._shape_q = q
        M = self.data.oMf[fid]
        p = M.translation.copy()
        p[2] = 0.0
        return pin.SE3(rz(math.atan2(M.rotation[1, 0], M.rotation[0, 0])), p)

    def _mid_stance(self, i, side):
        m = (i + self.stay[i][side]) // 2
        n = self.nodes[m]
        lateral = (1 if side == 'L' else -1) * 0.06
        p = n.com + rz(n.yaw)[:2, :2] @ np.array([0.0, lateral])
        return n, np.array([p[0], p[1], 0.0])

    def _next_single(self, i):
        """The side a single-support phase stands on, if one follows the phase of node i."""
        j = self.phase_last[self.nodes[i].phase_id] + 1
        if j < len(self.nodes) and len(self.nodes[j].contacts) == 1:
            return next(iter(self.nodes[j].contacts))
        return None

    def foothold(self, i, side, kind, frame_at):
        """Where a landing foot or knee takes the floor (an SE3 of the kind's frame on the floor)."""
        n = self.nodes[i]
        fid = self.fid[(side, kind)]
        before = self.nodes[i - 1].contacts if i > 0 else {}
        down = [s for s in n.contacts if s != side and s in before]           # already on the floor (a pivot may plant)
        together = [s for s in n.landing if s != side and s not in before]    # arriving from the air with this one
        if kind == 'knee' or (down and n.family not in ('travel', 'flight') and n.phase != 'turn'):
            # A step or a kneel inside a shape: where the shape this phase ends in puts it, measured from the foot
            # already down; the centre of mass moves over the phase to where that shape carries it.
            end = self.nodes[self.phase_last[n.phase_id]]
            M = self._shape_point(end, fid)
            if not down:
                return M
            a = down[0]
            A_ref = self.data.oMf[self.fid[(a, n.contacts[a])]].translation.copy()
            com_ref = pin.centerOfMass(self.model, self.data, self._shape_q).copy()
            A = frame_at(a, n.contacts[a]).translation
            p = A + (M.translation - A_ref)
            self.com_targets[n.phase_id] = (i, A[:2] + (com_ref[:2] - A_ref[:2]))
            return pin.SE3(M.rotation, np.array([p[0], p[1], 0.0]))
        if down and n.phase == 'turn':
            # A whirl's push: the free foot beside and behind the pivot.
            pivot = next(frame_at(s, k) for s, k in n.contacts.items() if s in down)
            off = np.array([FOURTH_M[0], (1 if side == 'L' else -1) * FOURTH_M[1], 0.0])
            p = pivot.translation + rz(n.yaw) @ off
            return pin.SE3(rz(n.yaw), np.array([p[0], p[1], 0.0]))
        if together:
            single = self._next_single(i)
            if single is not None:
                # Both feet land into a turn's preparation: the ball of the turning foot under the centre of mass
                # while it turns, the other foot's ball behind and beside it in fourth.
                nm, p = self._mid_stance(i, single)
                p[:2] = nm.com
                if side != single:
                    off = np.array([FOURTH_M[0], (1 if side == 'L' else -1) * FOURTH_M[1], 0.0])
                    p = p + rz(n.yaw) @ off
                    p[2] = 0.0
                M = pin.SE3(rz(n.yaw), p)
                return M * pin.SE3(np.eye(3), -self.toes_ahead[side]) if kind == 'flat' else M
            lateral = (1 if side == 'L' else -1) * STANCE_WIDTH
            p = np.array([n.com[0] - math.sin(n.yaw) * lateral, n.com[1] + math.cos(n.yaw) * lateral, 0.0])
            return pin.SE3(rz(n.yaw), p)
        # One foot from the air: under the centre of mass at mid-stance, a little to its own side.
        nm, p = self._mid_stance(i, side)
        return pin.SE3(rz(nm.yaw), p)

    def plan(self):
        """Footholds at every node from the plan, the contact targets they give, and the starting state. Cheap: no
        solver models are made, so a whole phrase can be planned at once and built window by window."""
        if getattr(self, 'x0', None) is not None:
            return self.x0
        self.com_targets = {}
        data = self.model.createData()
        q = self.reference_q(self.nodes[0])
        q[2] = self.q0[2]
        pin.framesForwardKinematics(self.model, data, q)
        foot, knee = {}, {}            # each foot's flat-sole placement and each knee's point, as they stand

        def frame_at(side, kind):
            if kind == 'knee':
                return knee[side]
            M = foot[side]
            return M if kind == 'flat' else M * pin.SE3(np.eye(3), self.toes_ahead[side])

        def set_frame(side, kind, M):
            if kind == 'knee':
                knee[side] = M
            else:
                foot[side] = M if kind == 'flat' else M * pin.SE3(np.eye(3), -self.toes_ahead[side])

        placed_at = []
        for i, node in enumerate(self.nodes):
            prev = self.nodes[i - 1].contacts if i > 0 else {}
            for side, kind in node.contacts.items():
                if i == 0:
                    # The first footholds stay where the starting shape puts the feet, laid flat on the floor.
                    M = data.oMf[self.fid[(side, kind)]]
                    p = M.translation.copy()
                    p[2] = 0.0
                    set_frame(side, kind, pin.SE3(rz(math.atan2(M.rotation[1, 0], M.rotation[0, 0])), p))
                elif side in node.landing and side in prev and prev[side] == 'pivot' and kind == 'ball':
                    # A pivot plants: the foot keeps its point and takes the dancer's present facing.
                    set_frame(side, 'ball', pin.SE3(rz(node.yaw), frame_at(side, 'pivot').translation))
                elif side in node.landing and not (side in prev and kind == 'flat'):
                    set_frame(side, kind, self.foothold(i, side, kind, frame_at))
                # (A heel coming down from demi-pointe leaves the foot where it stands.)
            placed_at.append({(s, k): frame_at(s, k) for s, k in node.contacts.items()})
        self.placed_at = placed_at
        # Where a step or kneel carries the centre of mass, its planned path moves there over the phase and stays
        # there until the next travelling unit.
        shift, k = np.zeros(2), 0
        while k < len(self.nodes):
            node = self.nodes[k]
            if k > 0 and node.family in ('travel', 'flight') and node.unit != self.nodes[k - 1].unit:
                shift = np.zeros(2)
            target = self.com_targets.get(node.phase_id)
            if target is not None and target[0] == k:
                last = self.phase_last[node.phase_id]
                start, end = node.com + shift, self.nodes[last].com.copy()
                for m in range(k, last + 1):
                    s = (m - k) / max(last - k, 1)
                    self.nodes[m].com = start + (target[1] - start) * s * s * (3 - 2 * s)
                shift = target[1] - end
                k = last + 1
                continue
            node.com = node.com + shift
            k += 1
        self.node_targets = [{self.fid[(s, k)]: (M if CONTACTS[k][1] == 6 else M.translation) for (s, k), M in pa.items()}
                             for pa in placed_at]
        q0 = place_feet(self.model, data, q, self.node_targets[0])
        # A free foot in the starting pose is lifted clear of the floor in the same inverse kinematics.
        pin.framesForwardKinematics(self.model, data, q0)
        targets = dict(self.node_targets[0])
        for side in ('L', 'R'):
            if side in self.nodes[0].contacts:
                continue
            low = min(data.oMf[self.model.getFrameId(f'{side}_{pt}')].translation[2] for pt in FLOOR_POINTS)
            if low < 0.02:
                fid = self.fid[(side, 'flat')]
                M = data.oMf[fid].copy()
                M.translation = M.translation + np.array([0.0, 0.0, 0.02 - low])
                targets[fid] = M
        if len(targets) > len(self.node_targets[0]):
            q0 = place_feet(self.model, data, q0, targets)
        self.x0 = np.concatenate([q0, np.zeros(self.model.nv)])
        return self.x0

    def state_at(self, a):
        """A state at node a that agrees with the plan, for a span cut from the phrase: the node's reference pose with
        its contacts at their planned footholds, moving with the plan's travel and turn as far as its contacts allow."""
        self.plan()
        if a == 0:
            return self.x0
        data = self.model.createData()
        n = self.nodes[a]
        q = self.reference_q(n)
        q[2] = self.q0[2]
        q = place_feet(self.model, data, q, self.node_targets[a], com=n.com)
        v = np.zeros(self.model.nv)
        R = pin.Quaternion(q[6], q[3], q[4], q[5]).toRotationMatrix()
        v[:3] = R.T @ np.array([n.com_vel[0], n.com_vel[1], 0.0])
        v[3:6] = R.T @ np.array([0.0, 0.0, n.yaw_rate])
        if n.contacts:
            pin.computeJointJacobians(self.model, data, q)
            pin.framesForwardKinematics(self.model, data, q)
            rows = [pin.getFrameJacobian(self.model, data, self.fid[(s, k)], pin.LOCAL_WORLD_ALIGNED)[:CONTACTS[k][1]]
                    for s, k in n.contacts.items()]
            rows += [pin.getFrameJacobian(self.model, data, self.model.getFrameId(f'{s}_{pt}'), pin.LOCAL_WORLD_ALIGNED)[2:3]
                     for s, k in n.contacts.items() if k == 'pivot' for pt in PIVOT_SUPPORTS]
            J = np.vstack(rows)
            v = v - J.T @ np.linalg.solve(J @ J.T + 1e-9 * np.eye(J.shape[0]), J @ v)
        return np.concatenate([q, v])

    def _dam(self, node, i, terminal=False, landing_targets=None):
        costs = self.costs(node, landing_targets, terminal=terminal)
        act = self.actuation_for(node)
        if node.contacts:
            return crocoddyl.DifferentialActionModelContactFwdDynamics(
                self.state, act, self.contact_model(node, self.placed_at[i]), costs, 0.0, True)
        return crocoddyl.DifferentialActionModelFreeFwdDynamics(self.state, act, costs)

    def running(self, a, b):
        """The running models of nodes a..b-1 as (node index, model): an impulse before each landing, then the node's
        integrated dynamics over DT."""
        out = []
        for i in range(a, b):
            node = self.nodes[i]
            if node.landing and i > 0:
                imp = crocoddyl.ImpulseModelMultiple(self.state)
                for side in node.landing:
                    kind = node.contacts[side]
                    impulse = crocoddyl.ImpulseModel6D if CONTACTS[kind][1] == 6 else crocoddyl.ImpulseModel3D
                    imp.addImpulse(f'{side}_impulse', impulse(self.state, self.fid[(side, kind)], pin.LOCAL_WORLD_ALIGNED))
                icosts = crocoddyl.CostModelSum(self.state, 0)
                xr = np.concatenate([self.reference_q(node), np.zeros(self.model.nv)])
                icosts.addCost('posture', crocoddyl.CostModelResidual(
                    self.state, crocoddyl.ActivationModelWeightedQuad(self.state_weights(node) ** 2),
                    crocoddyl.ResidualModelState(self.state, xr, 0)), 1.0)
                out.append((i, crocoddyl.ActionModelImpulseFwdDynamics(self.state, imp, icosts)))
            landing_targets = {}
            if i + 1 < len(self.nodes):
                nxt = self.nodes[i + 1]
                for side in nxt.landing:
                    kind = nxt.contacts[side]
                    if side in node.contacts and kind != 'knee':
                        continue                            # a heel or a pivot on a foot already down
                    M = self.placed_at[i + 1][(side, kind)]
                    landing_targets[side] = (self.fid[(side, kind)], M if CONTACTS[kind][1] == 6 else M.translation)
            iam = crocoddyl.IntegratedActionModelEuler(self._dam(node, i, landing_targets=landing_targets), DT)
            iam.u_lb = self.bounds['effort_lower'][6:]
            iam.u_ub = self.bounds['effort_upper'][6:]
            out.append((i, iam))
        return out

    def terminal(self, b):
        """The terminal model of a window ending before node b (the phrase's own end when b is past the last node)."""
        if b < len(self.nodes):
            return crocoddyl.IntegratedActionModelEuler(self._dam(self.nodes[b], b), 0.0)
        return crocoddyl.IntegratedActionModelEuler(self._dam(self.nodes[-1], len(self.nodes) - 1, terminal=True), 0.0)

    def reference_states(self, ms, x_start):
        """Warm-start states for a list of (node index, model): each node's reference pose with its contacts placed
        by inverse kinematics, the height carried along; an impulse keeps the state it receives."""
        data = self.model.createData()
        xs = [x_start]
        for k, (i, m) in enumerate(ms):
            if isinstance(m, crocoddyl.ActionModelImpulseFwdDynamics):
                xs.append(xs[-1].copy())
                continue
            j = min(i + 1, len(self.nodes) - 1)
            qn = self.reference_q(self.nodes[j])
            qn[2] = xs[-1][2]
            if self.nodes[j].contacts:
                qn = place_feet(self.model, data, qn, self.node_targets[j], iters=30, com=self.nodes[j].com)
            xs.append(np.concatenate([qn, np.zeros(self.model.nv)]))
        return xs

    def build(self):
        """The whole shooting problem at once (a short phrase or a span; about 12 MB a node)."""
        self.plan()
        self.ms = self.running(0, len(self.nodes))
        self.problem = crocoddyl.ShootingProblem(self.x0, [m for _, m in self.ms], self.terminal(len(self.nodes)))
        return self.problem

    def warm_start(self):
        xs = self.reference_states(self.ms, self.x0)
        return xs, self.problem.quasiStatic(xs[:-1])

    def solve(self, maxiter=200, verbose=False):
        self.build()
        solver = crocoddyl.SolverBoxFDDP(self.problem)
        if verbose:
            solver.setCallbacks([crocoddyl.CallbackVerbose()])
        xs, us = self.warm_start()
        t0 = time.time()
        solver.solve(xs, us, maxiter, False, 1e-3)
        self.seconds = time.time() - t0
        return Solution(self, self.ms, self.problem, solver, len(self.ms), solver.cost, solver.iter, solver.ffeas,
                        solver.stop)


def contact_forces(model, data):
    """The floor's force on each contact of a running node, by contact name ({} in flight or at an impulse)."""
    if isinstance(model, crocoddyl.ActionModelImpulseFwdDynamics):
        return {}
    try:
        contacts = data.differential.multibody.contacts.contacts.todict()
    except AttributeError:
        return {}
    return {name: np.r_[c.f.linear, c.f.angular] for name, c in contacts.items()}


class Solution:
    """A solved motion as records, one per running model kept: its node, whether it is an impulse, the state it
    starts from, its torques, its feedback gains and the floor's forces. Kept instead of the solver's models so that a
    whole phrase solved window by window never holds every node's model at once."""

    def __init__(self, prob, ms, problem, solver, keep, cost, iters, gap, stop=0.0, records=None, x_end=None):
        self.cost, self.iter, self.ffeas, self.stop = cost, iters, gap, stop
        if records is not None:
            self.records, self.x_end = records, x_end
            return
        self.records = [dict(node=i, impulse=isinstance(m, crocoddyl.ActionModelImpulseFwdDynamics),
                             x=np.array(solver.xs[k]), u=np.array(solver.us[k]), K=np.array(solver.K[k]),
                             forces=contact_forces(m, problem.runningDatas[k]))
                        for k, (i, m) in enumerate(ms[:keep])]
        self.x_end = np.array(solver.xs[keep])

    @property
    def xs(self):
        return [r['x'] for r in self.records] + [self.x_end]

    @property
    def us(self):
        return [r['u'] for r in self.records]

    @property
    def K(self):
        return [r['K'] for r in self.records]


def solve_span(prob, a, b, maxiter=150, verbose=False):
    """Solve nodes a..b-1 of the planned phrase on their own, from the plan's state at node a (a test of a piece)."""
    x = prob.state_at(a)
    ms = prob.running(a, b)
    if ms and isinstance(ms[0][1], crocoddyl.ActionModelImpulseFwdDynamics):
        ms = ms[1:]                                         # the span starts already standing in its contacts
    sub = crocoddyl.ShootingProblem(x, [m for _, m in ms], prob.terminal(b))
    xs = prob.reference_states(ms, x)
    us = list(sub.quasiStatic(xs[:-1]))
    solver = crocoddyl.SolverBoxFDDP(sub)
    if verbose:
        solver.setCallbacks([crocoddyl.CallbackVerbose()])
    t0 = time.time()
    solver.solve(xs, us, maxiter, False, 1e-3)
    prob.seconds = time.time() - t0
    return Solution(prob, ms, sub, solver, len(ms), solver.cost, solver.iter, solver.ffeas, solver.stop)


def solve_windows(prob, window_s=3.5, step_s=1.0, maxiter=150, verbose=False, log=None):
    """Solve the phrase window by window: each window of `window_s` starts from the state the previous one reached
    after `step_s`, so every window is short enough to converge and long enough to see a whole flight from its
    push to its landing. Only one window's models exist at a time; the kept pieces join into one trajectory that
    satisfies the dynamics throughout. `log`, a callable, hears each window's result as it finishes."""
    x = prob.plan()
    N = len(prob.nodes)
    W, S = int(round(window_s / DT)), int(round(step_s / DT))
    records, tail = [], None
    a, total_cost, iters, worst_gap = 0, 0.0, 0, 0.0
    t0 = time.time()
    while a < N:
        b = min(N, a + W)
        ms = prob.running(a, b)
        sub = crocoddyl.ShootingProblem(x, [m for _, m in ms], prob.terminal(b))
        warm_x = prob.reference_states(ms, x)
        warm_u = list(sub.quasiStatic(warm_x[:-1]))
        if tail is not None:
            # Where this window overlaps the last, start from the last window's solution.
            n = min(len(tail[1]), len(ms))
            warm_x[1:n + 1] = tail[0][1:n + 1]
            warm_u[:n] = tail[1][:n]
        solver = crocoddyl.SolverBoxFDDP(sub)
        if verbose:
            solver.setCallbacks([crocoddyl.CallbackVerbose()])
        solver.solve(warm_x, warm_u, maxiter, False, 1e-3)
        keep_nodes = b - a if b >= N else min(S, b - a)
        keep = sum(1 for i, _ in ms if i < a + keep_nodes)
        part = Solution(prob, ms, sub, solver, keep, solver.cost, solver.iter, solver.ffeas)
        records.extend(part.records)
        tail = (list(solver.xs)[keep:], list(solver.us)[keep:])
        total_cost += solver.cost
        iters += solver.iter
        worst_gap = max(worst_gap, solver.ffeas)
        if log:
            log(dict(nodes=[a, b], t=[round(prob.nodes[a].t, 2), round(prob.nodes[b - 1].t + DT, 2)], kept=keep_nodes,
                     iterations=solver.iter, cost=round(solver.cost, 3), gap=float(solver.ffeas),
                     seconds=round(time.time() - t0, 1)))
        x = part.x_end
        a += keep_nodes
        del sub, solver, ms
    prob.seconds = time.time() - t0
    return Solution(prob, None, None, None, 0, total_cost, iters, worst_gap, records=records, x_end=x)


# --- outputs -----------------------------------------------------------------------------------------------------------

def segment_frames(prob, sol, fps=60):
    """Segment world poses at the motion file's rate, from the solution's states (impulse records dropped)."""
    model, data = prob.model, prob.model.createData()
    names = [s['name'] for s in body_spec()['segments']]
    fids = [model.getFrameId(n) for n in names]
    states = [r['x'] for r in sol.records if not r['impulse']] + [sol.x_end]
    t_nodes = np.array([prob.nodes[r['node']].t for r in sol.records if not r['impulse']] +
                       [prob.nodes[sol.records[-1]['node']].t + DT])
    times = np.arange(t_nodes[0], t_nodes[-1] + 1e-9, 1 / fps)
    P, Q, C, QS = [], [], [], []
    for t in times:
        i = min(int(np.searchsorted(t_nodes, t, side='right')) - 1, len(states) - 2)
        a = (t - t_nodes[i]) / (t_nodes[i + 1] - t_nodes[i])
        q = pin.interpolate(model, states[i][:model.nq], states[i + 1][:model.nq], float(np.clip(a, 0, 1)))
        pin.framesForwardKinematics(model, data, q)
        P.append([data.oMf[f].translation.copy() for f in fids])
        Q.append([np.r_[pin.Quaternion(data.oMf[f].rotation).w, pin.Quaternion(data.oMf[f].rotation).vec()] for f in fids])
        C.append(pin.centerOfMass(model, data, q).copy())
        QS.append(q.copy())
    return times - times[0], np.array(P), np.array(Q), np.array(C), np.array(QS)


def contact_record(prob, sol, times):
    """Which floor points bear weight at each output frame and the floor's force on them: a foot's or knee's contact
    force shared equally among the points it rests on (flat: heels and metatarsal heads; demi-pointe and pivot:
    metatarsal heads and toe; knee: the kneecap), for overlays and analysis."""
    names = list(body_spec()['contacts'])
    per_node = []
    for r in sol.records:
        if r['impulse']:
            continue
        node = prob.nodes[r['node']]
        on = np.zeros(len(names), bool)
        force = np.zeros((len(names), 3))
        forces = {}
        for name, f in r['forces'].items():
            forces[name.split('_')[0]] = forces.get(name.split('_')[0], np.zeros(3)) + f[:3]
        for side, kind in node.contacts.items():
            pts = CONTACTS[kind][2]
            for pt in pts:
                i = names.index(f'{side}_{pt}')
                on[i] = True
                force[i] = forces.get(side, np.zeros(3)) / len(pts)
        per_node.append((node.t, on, force))
    t_nodes = np.array([n[0] for n in per_node]) - per_node[0][0]
    idx = np.clip(np.searchsorted(t_nodes, times, side='right') - 1, 0, len(per_node) - 1)
    return np.array([per_node[i][1] for i in idx]), np.array([per_node[i][2] for i in idx])


def lowest_sole(prob, qs):
    """The lowest any sole point goes over the motion (negative: through the floor)."""
    d = prob.model.createData()
    fids = [prob.model.getFrameId(f'{s}_{p}') for s in 'LR' for p in SOLE_POINTS]
    low = 1e9
    for q in qs:
        pin.framesForwardKinematics(prob.model, d, q)
        low = min(low, min(d.oMf[f].translation[2] for f in fids))
    return low


def contact_report(prob, sol):
    """Forces at each node: friction used (tangential over normal), the share of joint torque bounds used, and the
    spin (revolutions a second about the vertical)."""
    rows = []
    bound = np.maximum(np.abs(prob.bounds['effort_lower'][6:]), prob.bounds['effort_upper'][6:])
    nq = prob.model.nq
    for r in sol.records:
        if r['impulse']:
            continue
        node, x, u = prob.nodes[r['node']], r['x'], r['u']
        row = dict(t=round(node.t, 3), phase=node.phase, unit=node.unit, contacts={s: k for s, k in sorted(node.contacts.items())})
        fr, fn = 0.0, 0.0
        for name, f in r['forces'].items():
            fn += max(f[2], 0.0)
            if f[2] > 5.0:                                   # friction as a share of a real press on the floor
                fr = max(fr, math.hypot(f[0], f[1]) / f[2])
        row['normal_N'] = round(fn, 1)
        row['least_press_N'] = round(min((f[2] for f in r['forces'].values()), default=0.0), 1)
        row['friction_used'] = round(fr, 3)
        row['torque_share_max'] = round(float(np.max(np.abs(u) / bound)), 3)
        R = pin.Quaternion(x[6], x[3], x[4], x[5]).toRotationMatrix()
        row['spin_rev_s'] = round(float((R @ x[nq + 3:nq + 6])[2]) / (2 * math.pi), 3)
        rows.append(row)
    return rows


def realise(phrase, units=None, out_dir=RUNS, maxiter=200, tag=None, verbose=False, windows=False, span=None,
            finish_s=0.0, window_log=None):
    """Realise the phrase (or units[i:j] of it, or the nodes in span=(t0, t1) seconds) and write its motion file and
    report. finish_s holds the last shape at rest at the end."""
    placed = phrase.placed if units is None else phrase.placed[units[0]:units[1]]
    nodes = timeline(placed, finish_s=finish_s)
    prob = Problem(nodes)
    if span is not None:
        a = next(i for i, n in enumerate(nodes) if n.t >= span[0] - 1e-9)
        b = next((i for i, n in enumerate(nodes) if n.t >= span[1] - 1e-9), len(nodes))
        solver = solve_span(prob, a, b, maxiter=maxiter, verbose=verbose)
        nodes = nodes[a:b]
    elif windows:
        solver = solve_windows(prob, maxiter=maxiter, verbose=verbose, log=window_log)
    else:
        solver = prob.solve(maxiter=maxiter, verbose=verbose)
    times, P, Q, C, QS = segment_frames(prob, solver)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = tag or phrase.name.split(' (')[0].lower().replace(' ', '_')
    report = contact_report(prob, solver)
    meta = dict(kind='dynamic', gravity=gr.G, air_density=None, phrase=phrase.name,
                units=sorted({n.unit for n in nodes}, key=[n.unit for n in nodes].index), dt=DT,
                span_s=[round(nodes[0].t, 3), round(nodes[-1].t + DT, 3)], finish_s=finish_s,
                solver='crocoddyl.SolverBoxFDDP' + (' in overlapping windows' if windows else ''), iterations=solver.iter,
                cost=solver.cost, stop=solver.stop, dynamic_gap=float(solver.ffeas), seconds=round(prob.seconds, 1),
                evidence='Whole-body optimal control at lunar gravity: rigid-body dynamics with 6D foot contacts, point '
                         'contacts for pivots (with the floor\'s turning friction) and knees, wrench and friction-cone '
                         'costs, landing impulses and torque bounds; air drag left out.',
                producer='humanities/dance/motion/realise.py')
    on, force = contact_record(prob, solver, times)
    meta['floor'] = dict(lowest_sole_point_m=round(float(lowest_sole(prob, QS)), 4))
    save_motion(out_dir / f'{tag}.npz', times, P, Q, C, QS, meta, contact_on=on, contact_force=force)
    (out_dir / f'{tag}_report.json').write_text(json.dumps(dict(meta=meta, nodes=report), indent=1))
    return prob, solver, report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('phrase', choices=['canopy', 'suspension', 'test'])
    ap.add_argument('--units', default=None, help='i:j, a slice of the placed units')
    ap.add_argument('--span', default=None, help='t0:t1, seconds of the phrase to solve on their own')
    ap.add_argument('--maxiter', type=int, default=200)
    ap.add_argument('--windows', action='store_true', help='solve in overlapping windows (a whole phrase)')
    ap.add_argument('--finish', type=float, default=0.0, help='seconds to hold the last shape at rest')
    ap.add_argument('--tag', default=None)
    ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args()
    ph = {'canopy': gr.first_phrase, 'suspension': gr.suspension_study}.get(a.phrase, gr.first_phrase)()
    units = tuple(int(x) for x in a.units.split(':')) if a.units else None
    span = tuple(float(x) for x in a.span.split(':')) if a.span else None
    log = (lambda r: print(json.dumps(r), flush=True)) if a.windows else None
    prob, solver, report = realise(ph, units, maxiter=a.maxiter, verbose=a.verbose, windows=a.windows, span=span,
                                   finish_s=a.finish, tag=a.tag, window_log=log)
    print(f'{len(prob.nodes)} nodes, {solver.iter} iterations, cost {solver.cost:.4g}, gap {solver.ffeas:.2e}, '
          f'{prob.seconds:.0f} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
