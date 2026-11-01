"""Realise a phrase with whole-body dynamics: the score's units become an optimal control problem in Crocoddyl.

    python -m humanities.dance.motion.realise PHRASE [--units i:j] [--out DIR]   # needs the dance environment

The grammar fixes what happens and when: which feet bear weight in each phase, how long each flight lasts, which
shapes the body passes through, how fast it travels and turns. This module turns that into a multiple-shooting
problem over the dancer's full rigid-body dynamics at lunar gravity (body.json through pinmodel.py):

- each stance node carries 6D contacts for the feet on the floor (the flat foot's sole, or the toes' sole on
  demi-pointe) with a wrench-cone cost that keeps every contact force inside the friction pyramid, its centre of
  pressure inside the foot and its yaw moment within what friction allows;
- each flight node is free fall of the articulated body, so angular momentum is conserved and the centre of mass
  follows its ballistic arc exactly;
- each landing passes an impulse model that stops the landing foot, as a rigid contact does;
- joint torques stay inside the body's strength (box-constrained FDDP), joint ranges inside barrier costs;
- the shapes enter as posture targets interpolated through each phase, the travel as a centre-of-mass reference, the
  footsteps as foot placements, turns as the pelvis's yaw.

The solution is a torque, force and motion history that satisfies the equations of motion at every node; the replay
(replay.py) runs those torques through MuJoCo's independent contact dynamics as a check. Air drag on the body is left
out of the optimisation: the envelope shows it changes a jump's height by under 4% even with large sleeves.
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
FLOOR_MARGIN_M = 0.002
MAX_FOOT_FORCE_N = 900.0     # a foot's peak push, about 10 lunar body weights (1.7 of Earth's): a jump's take-off peak


# --- shapes to joint targets ---------------------------------------------------------------------------------------

def quat_axis_angle(axis, angle):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    return pin.Quaternion(math.cos(angle / 2), *(math.sin(angle / 2) * a))


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


def place_feet(model, data, q, targets, iters=80, damping=1e-6):
    """Inverse kinematics that puts each contact frame at its target placement (6D), moving the root and the joints by
    the least change: the soles flat on their footsteps."""
    q = q.copy()
    for _ in range(iters):
        pin.framesForwardKinematics(model, data, q)
        errs, Js = [], []
        for fid, M in targets.items():
            errs.append(pin.log6(data.oMf[fid].actInv(M)).vector)
            Js.append(pin.computeFrameJacobian(model, data, q, fid, pin.LOCAL))
        e = np.concatenate(errs)
        if np.linalg.norm(e) < 1e-7:
            break
        J = np.vstack(Js)
        dq = J.T @ np.linalg.solve(J @ J.T + damping * np.eye(J.shape[0]), e)
        q = pin.integrate(model, q, dq)
    return q


# --- the timeline ----------------------------------------------------------------------------------------------------

@dataclass
class Node:
    t: float
    contacts: dict                 # side -> 'flat' or 'ball'
    landing: list = field(default_factory=list)
    targets: dict = field(default_factory=dict)
    weights: dict = field(default_factory=dict)
    yaw: float = 0.0
    com: np.ndarray | None = None
    phase: str = ''
    unit: str = ''


def contact_kind(phase, unit):
    shapes = phase.get('shapes') or []
    legs = [s.get('legs') for s in shapes if isinstance(s, dict)]
    if phase['kind'] == 'turn' or 'releve' in legs or unit in ('skim',):
        return 'ball'
    return 'flat'


def timeline(placed, dt=DT):
    """Nodes every dt from the placed units: contacts from each phase's support, shape targets interpolated between
    the phase's key shapes (held between phases), travel integrated from the grammar's speeds and headings, yaw for
    turns."""
    nodes = []
    t_end = placed[-1].exit.t
    start = placed[0].t0
    times = np.arange(start, t_end - 1e-9, dt)
    phases = []
    for pl in placed:
        for ph in pl.phases:
            phases.append((pl, ph))
    # Travel: speed linear across each unit, heading interpolated; flights keep their take-off velocity.
    pos = np.array([placed[0].entry.x, placed[0].entry.y])
    yaw_turn = 0.0
    last_shape = {}
    prev_contacts = None
    k = 0
    for t in times:
        while k < len(phases) - 1 and t >= phases[k][1]['t0'] + phases[k][1]['duration'] - 1e-9:
            k += 1
        pl, ph = phases[k]
        frac = (t - ph['t0']) / max(ph['duration'], 1e-9)
        ufrac = (t - pl.t0) / max(pl.duration, 1e-9)
        speed = pl.entry.speed + (pl.exit.speed - pl.entry.speed) * ufrac
        dh = ((pl.exit.heading - pl.entry.heading + 180) % 360) - 180
        heading = math.radians(pl.entry.heading + dh * ufrac)
        pos = pos + dt * speed * np.array([math.cos(heading), math.sin(heading)])
        if ph.get('turn_rev'):
            revs = gr.evaluate(ph['turn_rev'], pl.params) if isinstance(ph['turn_rev'], str) else float(ph['turn_rev'])
            yaw_turn += 2 * math.pi * revs * dt / ph['duration'] * (1 if pl.side == 'L' else -1)
        sup = ph.get('support', 'same')
        if sup in ('L', 'R'):
            contacts = {sup: contact_kind(ph, pl.unit)}
        elif sup == 'both':
            kind = contact_kind(ph, pl.unit)
            contacts = {'L': kind, 'R': kind}
        elif sup in ('none',) or ph['kind'] == 'flight':
            contacts = {}
        else:
            contacts = dict(prev_contacts or {'L': 'flat', 'R': 'flat'})
        landing = [s for s in contacts if prev_contacts is not None and s not in prev_contacts]
        # Shape targets: the last key shape reached, blended toward the next within the phase.
        shapes = ph.get('shapes') if isinstance(ph.get('shapes'), list) else []
        target = dict(last_shape)
        for s in shapes:
            if s['at'] <= frac + 1e-9:
                target = s
        nxt = next((s for s in shapes if s['at'] > frac + 1e-9), None)
        node = Node(t=float(t), contacts=contacts, landing=landing, yaw=heading + yaw_turn, phase=ph['kind'], unit=pl.unit)
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
        node.targets = tg
        node.com = pos.copy()
        nodes.append(node)
        if shapes:
            reached = [s for s in shapes if s['at'] <= frac + 1e-9]
            if reached:
                last_shape = {k2: v for k2, v in reached[-1].items() if k2 != 'at'}
        prev_contacts = contacts
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
        self.frame = {s: {'flat': self.model.getFrameId(f'{s}_sole'), 'ball': self.model.getFrameId(f'{s}_toesole')}
                      for s in ('L', 'R')}
        self.q0 = neutral(self.model)
        self.x_ref_cache = {}

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
        q[3:7] = pin.Quaternion(pin.utils.rpyToMatrix(0, 0, node.yaw)).coeffs()
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
        lo = np.concatenate([self.model.lowerPositionLimit, -1e3 * np.ones(self.model.nv)])
        hi = np.concatenate([self.model.upperPositionLimit, 1e3 * np.ones(self.model.nv)])
        xlo, xhi = self.state.diff(np.zeros(1), np.zeros(1)) if False else (None, None)
        bounds = crocoddyl.ActivationBounds(self._xlb(), self._xub(), 1.0)
        cm.addCost('ranges', crocoddyl.CostModelResidual(self.state, crocoddyl.ActivationModelQuadraticBarrier(bounds),
                                                         crocoddyl.ResidualModelState(self.state, self.nu)), 1e3)
        # Facing: the pelvis's yaw follows the travel or the turn; it stays near upright otherwise.
        R = pin.utils.rpyToMatrix(0, 0, node.yaw)
        cm.addCost('facing', crocoddyl.CostModelResidual(
            self.state, crocoddyl.ActivationModelWeightedQuad(np.array([1.0, 1.0, 4.0]) ** 2),
            crocoddyl.ResidualModelFrameRotation(self.state, self.model.getFrameId('pelvis'), R, self.nu)), 10.0)
        # Travel: the centre of mass follows the planned path in the horizontal while the feet are down.
        if node.contacts and node.com is not None:
            com_ref = np.array([node.com[0], node.com[1], 0.0])
            cm.addCost('travel', crocoddyl.CostModelResidual(
                self.state, crocoddyl.ActivationModelWeightedQuad(np.array([1.0, 1.0, 0.0])),
                crocoddyl.ResidualModelCoMPosition(self.state, com_ref, self.nu)), 20.0)
        # Contacts: forces inside the friction pyramid, centre of pressure inside the foot, yaw moment within friction.
        for side, kind in node.contacts.items():
            fid = self.frame[side][kind]
            cone = crocoddyl.WrenchCone(np.eye(3), MU, np.array(FOOT_HALF[kind]) * 2, 4, True, 1.0, MAX_FOOT_FORCE_N)
            res = crocoddyl.ResidualModelContactWrenchCone(self.state, fid, cone, self.nu, True)
            b = crocoddyl.ActivationBounds(cone.lb, cone.ub)
            cm.addCost(f'{side}_cone', crocoddyl.CostModelResidual(self.state, crocoddyl.ActivationModelQuadraticBarrier(b), res), 1e2)
            # Light cost on sideways forces and moments: the feet do not squeeze against each other for nothing.
            cm.addCost(f'{side}_shear', crocoddyl.CostModelResidual(
                self.state, crocoddyl.ActivationModelWeightedQuad(np.array([1.0, 1.0, 0.0, 30.0, 30.0, 30.0]) ** 2),
                crocoddyl.ResidualModelContactForce(self.state, fid, pin.Force.Zero(), 6, self.nu, True)), 1e-3)
        # The floor: every sole point that bears no weight at this node stays at or above it (a 2 mm margin), so swing
        # feet, toes under a flat foot and heels over demi-pointe never pass through the floor.
        for side in ('L', 'R'):
            kind = node.contacts.get(side)
            for pt in SOLE_POINTS:
                if kind == 'flat' and pt != 'toe':
                    continue                       # held on the floor by the flat foot's contact
                if kind == 'ball' and pt in ('ball_m', 'ball_l', 'toe'):
                    continue                       # on the floor with the toes
                fid = self.model.getFrameId(f'{side}_{pt}')
                floor = crocoddyl.ActivationBounds(np.array([-1e6, -1e6, FLOOR_MARGIN_M]), np.array([1e6, 1e6, 1e6]))
                cm.addCost(f'{side}_{pt}_floor', crocoddyl.CostModelResidual(
                    self.state, crocoddyl.ActivationModelQuadraticBarrier(floor),
                    crocoddyl.ResidualModelFrameTranslation(self.state, fid, np.zeros(3), self.nu)), 1e4)
        # Landings: the landing foot arrives at its planned placement.
        for side, M in (landing_targets or {}).items():
            fid = self.frame[side]['flat']
            cm.addCost(f'{side}_land', crocoddyl.CostModelResidual(
                self.state, crocoddyl.ResidualModelFramePlacement(self.state, fid, M, self.nu)), 1e4)
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

    def contact_model(self, node, placements):
        cmult = crocoddyl.ContactModelMultiple(self.state, self.nu)
        for side, kind in node.contacts.items():
            fid = self.frame[side][kind]
            cmult.addContact(f'{side}_{kind}', crocoddyl.ContactModel6D(
                self.state, fid, placements[(side, kind)], pin.LOCAL_WORLD_ALIGNED, self.nu, np.array([0.0, 30.0])))
        return cmult

    def build(self):
        """The shooting problem: footsteps placed from the plan, a running model per node, impulses at landings."""
        data = self.model.createData()
        q = self.reference_q(self.nodes[0])
        q[2] = self.q0[2]
        pin.framesForwardKinematics(self.model, data, q)
        placements = {}
        running = []
        self.node_targets = []
        for i, node in enumerate(self.nodes):
            # Footsteps: a foot that lands is placed beside the planned path at its landing time.
            for side, kind in node.contacts.items():
                key = (side, kind)
                if key not in placements or side in node.landing:
                    lateral = (1 if side == 'L' else -1) * (STANCE_WIDTH if len(node.contacts) == 2 else 0.06)
                    yaw = node.yaw
                    p = np.array([node.com[0] - math.sin(yaw) * lateral, node.com[1] + math.cos(yaw) * lateral, 0.0])
                    if i == 0:
                        # The first footsteps stay where the starting shape puts the feet, laid flat on the floor.
                        p = data.oMf[self.frame[side][kind]].translation.copy()
                        p[2] = 0.0
                    placements[key] = pin.SE3(pin.utils.rpyToMatrix(0, 0, yaw), p)
                    other_kind = 'ball' if kind == 'flat' else 'flat'
                    # The other frame of the same foot moves with it (a flat foot rises onto its toes in place).
                    off = np.array([0.113 if kind == 'flat' else -0.113, 0.0, 0.0])
                    placements[(side, other_kind)] = pin.SE3(placements[key].rotation, p + placements[key].rotation @ off)
            self.node_targets.append({self.frame[s][k]: placements[(s, k)] for s, k in node.contacts.items()})
            landing_targets = {}
            if i + 1 < len(self.nodes):
                nxt = self.nodes[i + 1]
                for side in nxt.landing:
                    yaw = nxt.yaw
                    lateral = (1 if side == 'L' else -1) * (STANCE_WIDTH if len(nxt.contacts) == 2 else 0.06)
                    p = np.array([nxt.com[0] - math.sin(yaw) * lateral, nxt.com[1] + math.cos(yaw) * lateral, 0.0])
                    landing_targets[side] = pin.SE3(pin.utils.rpyToMatrix(0, 0, yaw), p)
            costs = self.costs(node, landing_targets)
            if node.contacts:
                dam = crocoddyl.DifferentialActionModelContactFwdDynamics(
                    self.state, self.actuation, self.contact_model(node, placements), costs, 0.0, True)
            else:
                dam = crocoddyl.DifferentialActionModelFreeFwdDynamics(self.state, self.actuation, costs)
            iam = crocoddyl.IntegratedActionModelEuler(dam, DT)
            iam.u_lb = self.bounds['effort_lower'][6:]
            iam.u_ub = self.bounds['effort_upper'][6:]
            if node.landing and i > 0:
                imp = crocoddyl.ImpulseModelMultiple(self.state)
                for side in node.landing:
                    kind = node.contacts[side]
                    imp.addImpulse(f'{side}_impulse', crocoddyl.ImpulseModel6D(self.state, self.frame[side][kind], pin.LOCAL_WORLD_ALIGNED))
                icosts = crocoddyl.CostModelSum(self.state, 0)
                xr = np.concatenate([self.reference_q(node), np.zeros(self.model.nv)])
                icosts.addCost('posture', crocoddyl.CostModelResidual(
                    self.state, crocoddyl.ActivationModelWeightedQuad(self.state_weights(node) ** 2),
                    crocoddyl.ResidualModelState(self.state, xr, 0)), 1.0)
                running.append(crocoddyl.ActionModelImpulseFwdDynamics(self.state, imp, icosts))
            running.append(iam)
        last = self.nodes[-1]
        tcost = self.costs(last, terminal=True)
        if last.contacts:
            tdam = crocoddyl.DifferentialActionModelContactFwdDynamics(
                self.state, self.actuation, self.contact_model(last, placements), tcost, 0.0, True)
        else:
            tdam = crocoddyl.DifferentialActionModelFreeFwdDynamics(self.state, self.actuation, tcost)
        terminal = crocoddyl.IntegratedActionModelEuler(tdam, 0.0)
        q0 = place_feet(self.model, data, q, self.node_targets[0])
        # A free foot in the starting pose is lifted clear of the floor in the same inverse kinematics.
        pin.framesForwardKinematics(self.model, data, q0)
        targets = dict(self.node_targets[0])
        for side in ('L', 'R'):
            if side in self.nodes[0].contacts:
                continue
            low = min(data.oMf[self.model.getFrameId(f'{side}_{pt}')].translation[2] for pt in SOLE_POINTS)
            if low < 0.02:
                fid = self.frame[side]['flat']
                M = data.oMf[fid].copy()
                M.translation = M.translation + np.array([0.0, 0.0, 0.02 - low])
                targets[fid] = M
        if len(targets) > len(self.node_targets[0]):
            q0 = place_feet(self.model, data, q0, targets)
        x0 = np.concatenate([q0, np.zeros(self.model.nv)])
        self.problem = crocoddyl.ShootingProblem(x0, running, terminal)
        self.x0 = x0
        return self.problem

    def warm_start(self):
        xs = []
        for i, m in enumerate(self.problem.runningModels):
            pass
        data = self.model.createData()
        xs = [self.x0]
        k = 0
        for m in self.problem.runningModels:
            if isinstance(m, crocoddyl.ActionModelImpulseFwdDynamics):
                xs.append(xs[-1].copy())
                continue
            k += 1
            node = self.nodes[min(k, len(self.nodes) - 1)]
            qn = self.reference_q(node)
            qn[2] = xs[-1][2]
            if node.contacts:
                qn = place_feet(self.model, data, qn, self.node_targets[min(k, len(self.nodes) - 1)], iters=30)
            xs.append(np.concatenate([qn, np.zeros(self.model.nv)]))
        us = self.problem.quasiStatic(xs[:-1])
        return xs, us

    def solve(self, maxiter=200, verbose=False):
        self.build()
        solver = crocoddyl.SolverBoxFDDP(self.problem)
        if verbose:
            solver.setCallbacks([crocoddyl.CallbackVerbose()])
        xs, us = self.warm_start()
        t0 = time.time()
        solver.solve(xs, us, maxiter, False, 1e-3)
        self.solver = solver
        self.seconds = time.time() - t0
        return solver


class WindowedSolution:
    """The concatenated result of solving a long phrase in overlapping windows (xs, us, K and the models used)."""

    def __init__(self, xs, us, K, models, cost, iters, gap):
        self.xs, self.us, self.K, self.models = xs, us, K, models
        self.cost, self.iter, self.ffeas, self.stop = cost, iters, gap, 0.0
        self.problem = type('P', (), {'runningDatas': [m.createData() for m in models]})()
        for m, d, x, u in zip(models, self.problem.runningDatas, xs, us):
            m.calc(d, x, u)


def solve_windows(prob, window_s=3.5, step_s=1.0, maxiter=150, verbose=False):
    """Solve the phrase window by window: each window of `window_s` starts from the state the previous one reached
    after `step_s`, so every window is short enough to converge and long enough to see a whole flight from its
    push to its landing. The kept pieces join into one trajectory that satisfies the dynamics throughout."""
    prob.build()
    models = list(prob.problem.runningModels)
    W, S = int(round(window_s / DT)), int(round(step_s / DT))
    x = prob.x0
    xs_ref, us_ref = prob.warm_start()
    xs_out, us_out, K_out, used = [x], [], [], []
    i, total_cost, iters, worst_gap = 0, 0.0, 0, 0.0
    t0 = time.time()
    while i < len(models):
        j = min(len(models), i + W)
        while j < len(models) and isinstance(models[j], crocoddyl.ActionModelImpulseFwdDynamics):
            j += 1
        if j < len(models):
            term = crocoddyl.IntegratedActionModelEuler(models[j].differential, 0.0)
        else:
            term = prob.problem.terminalModel
        sub = crocoddyl.ShootingProblem(x, models[i:j], term)
        warm_x = [x] + [xs_ref[k + 1] for k in range(i, j)]
        warm_u = [us_ref[k] for k in range(i, j)]
        solver = crocoddyl.SolverBoxFDDP(sub)
        if verbose:
            solver.setCallbacks([crocoddyl.CallbackVerbose()])
        solver.solve(warm_x, warm_u, maxiter, False, 1e-3)
        keep = (j - i) if j >= len(models) else min(S, j - i)
        sx, su, sK = list(solver.xs), list(solver.us), list(solver.K)
        xs_out.extend(sx[1:keep + 1])
        us_out.extend(su[:keep])
        K_out.extend(sK[:keep])
        used.extend(models[i:i + keep])
        # The next window's warm start continues from this window's solution.
        for k in range(keep, j - i):
            xs_ref[i + k + 1] = sx[k + 1]
            us_ref[i + k] = su[k]
        total_cost += solver.cost
        iters += solver.iter
        worst_gap = max(worst_gap, solver.ffeas)
        x = sx[keep]
        i += keep
    prob.seconds = time.time() - t0
    return WindowedSolution(xs_out, us_out, K_out, used, total_cost, iters, worst_gap)


# --- outputs -----------------------------------------------------------------------------------------------------------

def segment_frames(prob, xs, fps=60, models=None):
    """Segment world poses at the motion file's rate, from the solution's states (impulse nodes dropped)."""
    model, data = prob.model, prob.model.createData()
    names = [s['name'] for s in body_spec()['segments']]
    fids = [model.getFrameId(n) for n in names]
    xs = list(xs)
    models = models if models is not None else list(prob.problem.runningModels)
    states = [x for x, m in zip(xs[:-1], models) if not isinstance(m, crocoddyl.ActionModelImpulseFwdDynamics)]
    states.append(xs[-1])
    t_nodes = np.array([n.t for n in prob.nodes] + [prob.nodes[-1].t + DT])[:len(states)]
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


def contact_record(prob, solver, times):
    """Which sole points bear weight at each output frame and the floor's force on them: a foot's 6D contact force
    shared equally among the points it stands on (flat: heels and metatarsal heads; demi-pointe: metatarsal heads and
    toe), for overlays and analysis."""
    names = list(body_spec()['contacts'])
    models = getattr(solver, 'models', None) or list(prob.problem.runningModels)
    datas = solver.problem.runningDatas
    per_node = []
    k = 0
    for m, d in zip(models, datas):
        if isinstance(m, crocoddyl.ActionModelImpulseFwdDynamics):
            continue
        node = prob.nodes[min(k, len(prob.nodes) - 1)]
        k += 1
        on = np.zeros(len(names), bool)
        force = np.zeros((len(names), 3))
        if node.contacts:
            forces = {name.split('_')[0]: c.f.linear.copy() for name, c in d.differential.multibody.contacts.contacts.todict().items()}
            for side, kind in node.contacts.items():
                pts = ('heel_m', 'heel_l', 'ball_m', 'ball_l') if kind == 'flat' else ('ball_m', 'ball_l', 'toe')
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


def contact_report(prob, solver):
    """Forces at each stance node: friction used (tangential over normal), the share of joint torque bounds used."""
    rows = []
    k = 0
    models = getattr(solver, 'models', None) or list(prob.problem.runningModels)
    for m, d, u in zip(models, solver.problem.runningDatas, solver.us):
        if isinstance(m, crocoddyl.ActionModelImpulseFwdDynamics):
            continue
        node = prob.nodes[min(k, len(prob.nodes) - 1)]
        k += 1
        row = dict(t=round(node.t, 3), phase=node.phase, unit=node.unit, contacts=sorted(node.contacts))
        fr = 0.0
        fn = 0.0
        if node.contacts:
            for name, c in d.differential.multibody.contacts.contacts.todict().items():
                f = c.f.linear
                fn += max(f[2], 0.0)
                fr = max(fr, math.hypot(f[0], f[1]) / max(f[2], 1e-6))
        row['normal_N'] = round(fn, 1)
        row['friction_used'] = round(fr, 3)
        bound = np.maximum(np.abs(prob.bounds['effort_lower'][6:]), prob.bounds['effort_upper'][6:])
        row['torque_share_max'] = round(float(np.max(np.abs(u) / bound)), 3)
        rows.append(row)
    return rows


def realise(phrase, units=None, out_dir=RUNS, maxiter=200, tag=None, verbose=False, windows=False):
    placed = phrase.placed if units is None else phrase.placed[units[0]:units[1]]
    nodes = timeline(placed)
    prob = Problem(nodes)
    if windows:
        solver = solve_windows(prob, maxiter=maxiter, verbose=verbose)
    else:
        solver = prob.solve(maxiter=maxiter, verbose=verbose)
    times, P, Q, C, QS = segment_frames(prob, solver.xs, models=getattr(solver, 'models', None))
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = tag or phrase.name.split(' (')[0].lower().replace(' ', '_')
    report = contact_report(prob, solver)
    meta = dict(kind='dynamic', gravity=gr.G, air_density=None, phrase=phrase.name,
                units=[pl.unit for pl in placed], dt=DT,
                solver='crocoddyl.SolverBoxFDDP' + (' in overlapping windows' if windows else ''), iterations=solver.iter,
                cost=solver.cost, stop=solver.stop, dynamic_gap=float(solver.ffeas), seconds=round(prob.seconds, 1),
                evidence='Whole-body optimal control at lunar gravity: rigid-body dynamics with 6D foot contacts, '
                         'wrench-cone costs, landing impulses and torque bounds; air drag left out.',
                producer='humanities/dance/motion/realise.py')
    on, force = contact_record(prob, solver, times)
    meta['floor'] = dict(lowest_sole_point_m=round(float(lowest_sole(prob, QS)), 4))
    save_motion(out_dir / f'{tag}.npz', times, P, Q, C, QS, meta, contact_on=on, contact_force=force)
    (out_dir / f'{tag}_report.json').write_text(json.dumps(dict(meta=meta, nodes=report), indent=1))
    return prob, solver, report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('phrase', choices=['canopy', 'suspension', 'test'])
    ap.add_argument('--units', default=None)
    ap.add_argument('--maxiter', type=int, default=200)
    ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args()
    ph = {'canopy': gr.first_phrase, 'suspension': gr.suspension_study}.get(a.phrase, gr.first_phrase)()
    units = tuple(int(x) for x in a.units.split(':')) if a.units else None
    prob, solver, report = realise(ph, units, maxiter=a.maxiter, verbose=a.verbose)
    print(f'{len(prob.nodes)} nodes, {solver.iter} iterations, cost {solver.cost:.4g}, gap {solver.ffeas:.2e}, '
          f'{prob.seconds:.0f} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
