"""The dancer as a Pinocchio model, built from body.json, with converters to and from MuJoCo's coordinates.

The optimiser (Crocoddyl) and the replay (MuJoCo) must agree on one body. Both are built from body.json: this module
adds every joint and inertia in the spec's order, places frames at the sole and palm points and at each foot's sole
centre for flat and demi-pointe contacts, and maps generalised coordinates between the two programs by joint name,
since their joint orders and quaternion conventions differ (Pinocchio stores x, y, z, w; MuJoCo w, x, y, z).
"""
from __future__ import annotations

import numpy as np
import pinocchio as pin

from humanities.dance.motion_io import body_spec
from shared.constants import MOON_SURFACE_GRAVITY

SOLE_FRAMES = {   # 6D contact frames on the floor at rest: the flat foot's centre and the toes' centre
    'foot': lambda spec_seg: (0.5 * (-0.045 + 0.129), 0.0),
    'toes': lambda spec_seg: (0.030, 0.0),
}


def build(gravity: float = MOON_SURFACE_GRAVITY):
    sp = body_spec()
    model = pin.Model()
    model.name = 'terluna_dancer'
    last_joint = {}
    seg_frame = {}
    lower, upper, effort_lo, effort_hi = {}, {}, {}, {}
    # Depth-first order: Pinocchio's algorithms take each subtree's velocities as one contiguous block.
    children = {}
    for seg in sp['segments']:
        children.setdefault(seg['parent'], []).append(seg)
    ordered = []

    def visit(seg):
        ordered.append(seg)
        for ch in children.get(seg['name'], []):
            visit(ch)

    for root in children[None]:
        visit(root)
    for seg in ordered:
        if seg['parent'] is None:
            # The root's coordinates are the pelvis's world pose, as MuJoCo's free joint stores them.
            jid = model.addJoint(0, pin.JointModelFreeFlyer(), pin.SE3.Identity(), 'root_joint')   # the name Crocoddyl looks for
        else:
            jid = last_joint[seg['parent']]
            placement = pin.SE3(np.eye(3), np.array(seg['offset'], float))
            for k, j in enumerate(seg['joints']):
                if j['kind'] == 'hinge':
                    jm = pin.JointModelRevoluteUnaligned(np.array(j['axis'], float))
                else:
                    jm = pin.JointModelSpherical()
                jid = model.addJoint(jid, jm, placement if k == 0 else pin.SE3.Identity(), j['name'])
                if j['kind'] == 'hinge':
                    lower[j['name']] = np.radians(j['range_deg'][0])
                    upper[j['name']] = np.radians(j['range_deg'][1])
                    effort_lo[j['name']] = [-j['torque_nm'][0]]
                    effort_hi[j['name']] = [j['torque_nm'][1]]
                else:
                    effort_lo[j['name']] = [-t for t in j['torque_nm']]
                    effort_hi[j['name']] = list(j['torque_nm'])
        I = np.diag(seg['inertia_kg_m2'])
        model.appendBodyToJoint(jid, pin.Inertia(seg['mass_kg'], np.array(seg['com'], float), I), pin.SE3.Identity())
        parent_frame = seg_frame.get(seg['parent'], 0)
        seg_frame[seg['name']] = model.addBodyFrame(seg['name'], jid, pin.SE3.Identity(), parent_frame)
        last_joint[seg['name']] = jid
    for name, c in sp['contacts'].items():
        jid = last_joint[c['segment']]
        model.addFrame(pin.Frame(name, jid, seg_frame[c['segment']], pin.SE3(np.eye(3), np.array(c['pos'], float)),
                                 pin.FrameType.OP_FRAME))
    by = {s['name']: s for s in sp['segments']}
    for side in ('L', 'R'):
        foot, toes = by[f'{side}_foot'], by[f'{side}_toes']
        # The flat foot's sole centre (between the heel and the metatarsal heads) and the toes' sole centre, on the floor.
        heel_x = sp['contacts'][f'{side}_heel_m']['pos'][0]
        ball_x = sp['contacts'][f'{side}_ball_m']['pos'][0]
        pf = np.array([0.5 * (heel_x + ball_x), 0.0, -foot['origin'][2]])
        model.addFrame(pin.Frame(f'{side}_sole', last_joint[f'{side}_foot'], seg_frame[f'{side}_foot'], pin.SE3(np.eye(3), pf),
                                 pin.FrameType.OP_FRAME))
        pt = np.array([0.028, 0.0, -toes['origin'][2]])
        model.addFrame(pin.Frame(f'{side}_toesole', last_joint[f'{side}_toes'], seg_frame[f'{side}_toes'], pin.SE3(np.eye(3), pt),
                                 pin.FrameType.OP_FRAME))
    model.gravity.linear = np.array([0.0, 0.0, -gravity])
    lo, hi = model.lowerPositionLimit.copy(), model.upperPositionLimit.copy()
    elo, ehi = -model.effortLimit.copy(), model.effortLimit.copy()
    for jid in range(1, model.njoints):
        name = model.names[jid]
        iq, iv = model.idx_qs[jid], model.idx_vs[jid]
        if name in lower:
            lo[iq], hi[iq] = lower[name], upper[name]
        if name in effort_lo:
            nv = model.nvs[jid]
            elo[iv:iv + nv] = effort_lo[name]
            ehi[iv:iv + nv] = effort_hi[name]
    model.lowerPositionLimit, model.upperPositionLimit = lo, hi
    model.effortLimit = np.maximum(np.abs(elo), np.abs(ehi))
    return model, dict(effort_lower=elo, effort_upper=ehi)


def neutral(model):
    """The rest pose standing on the floor: the pelvis at its rest height."""
    q = pin.neutral(model)
    q[2] = body_spec()['segments'][0]['origin'][2]
    return q


class MujocoMap:
    """Converts configurations and velocities between Pinocchio (free flyer x y z, qx qy qz qw; balls qx qy qz qw;
    velocities in local frames) and MuJoCo (free joint x y z, qw qx qy qz with linear velocity in the world frame and
    angular in the local frame; balls qw qx qy qz)."""

    def __init__(self, model, mj_model):
        import mujoco
        self.mj = mj_model
        self.pairs = []
        for jid in range(2, model.njoints):
            name = model.names[jid]
            mid = mujoco.mj_name2id(mj_model, mujoco.mjtObj.mjOBJ_JOINT, name)
            assert mid >= 0, name
            self.pairs.append((model.idx_qs[jid], model.nqs[jid], model.idx_vs[jid], model.nvs[jid],
                               mj_model.jnt_qposadr[mid], mj_model.jnt_dofadr[mid]))

    @staticmethod
    def _xyzw_to_wxyz(q):
        return np.r_[q[3], q[:3]]

    @staticmethod
    def _wxyz_to_xyzw(q):
        return np.r_[q[1:], q[0]]

    def to_mujoco(self, q, v=None):
        qpos = np.zeros(self.mj.nq)
        qpos[:3] = q[:3]
        qpos[3:7] = self._xyzw_to_wxyz(q[3:7])
        for iq, nq, iv, nv, aq, av in self.pairs:
            qpos[aq:aq + nq] = self._xyzw_to_wxyz(q[iq:iq + 4]) if nq == 4 else q[iq:iq + nq]
        if v is None:
            return qpos
        qvel = np.zeros(self.mj.nv)
        R = pin.Quaternion(q[6], q[3], q[4], q[5]).toRotationMatrix()
        qvel[:3] = R @ v[:3]                       # Pinocchio's free-flyer linear velocity is in the local frame
        qvel[3:6] = v[3:6]                         # both keep the angular velocity in the local frame
        for iq, nq, iv, nv, aq, av in self.pairs:
            qvel[av:av + nv] = v[iv:iv + nv]
        return qpos, qvel

    def from_mujoco(self, qpos, qvel=None):
        q = np.zeros(sum(n for _, n, *_ in self.pairs) + 7)
        q[:3] = qpos[:3]
        q[3:7] = self._wxyz_to_xyzw(qpos[3:7])
        for iq, nq, iv, nv, aq, av in self.pairs:
            q[iq:iq + nq] = self._wxyz_to_xyzw(qpos[aq:aq + 4]) if nq == 4 else qpos[aq:aq + nq]
        if qvel is None:
            return q
        v = np.zeros(sum(n for *_, n, _, _ in self.pairs) + 6)
        R = pin.Quaternion(q[6], q[3], q[4], q[5]).toRotationMatrix()
        v[:3] = R.T @ qvel[:3]
        v[3:6] = qvel[3:6]
        for iq, nq, iv, nv, aq, av in self.pairs:
            v[iv:iv + nv] = qvel[av:av + nv]
        return q, v
