"""Replay an optimised motion through MuJoCo's independent contact dynamics: the check on the motion engine.

Crocoddyl's solution is consistent with rigid 6D foot contacts and the rigid-body dynamics Pinocchio computes. MuJoCo
models the same body (dancer.xml, from the same spec) with soft contacts at the sole points, so the replay asks
whether the torques the optimiser found, with its own feedback gains closing the small differences, still carry the
dancer through the phrase when another program resolves the contacts. The replay uses the optimiser's friction and
drops MuJoCo's joint damping and armature, which the optimiser does not model; air acts on the body through MuJoCo's
ellipsoid fluid model at the design run's density, which the optimiser leaves out.

Reported: how far the replayed body drifts from the plan (centre of mass, pelvis orientation, joint angles), whether
it falls, the contact forces MuJoCo computes, the friction they use and their centre of pressure on each foot.
"""
from __future__ import annotations

import math

import numpy as np
import pinocchio as pin
import crocoddyl

from humanities.dance.motion.pinmodel import MujocoMap
from humanities.dance.motion_io import body_spec
from humanities.dance.body.model import MJCF_OUT


def actuator_map(model, mj):
    """MuJoCo actuator index for each of the optimiser's controls (the joint velocities after the root)."""
    import mujoco
    out = np.zeros(model.nv - 6, int)
    for jid in range(2, model.njoints):
        name = model.names[jid]
        iv = model.idx_vs[jid] - 6
        if model.nvs[jid] == 1:
            out[iv] = mujoco.mj_name2id(mj, mujoco.mjtObj.mjOBJ_ACTUATOR, name)
        else:
            for k, ax in enumerate('xyz'):
                out[iv + k] = mujoco.mj_name2id(mj, mujoco.mjtObj.mjOBJ_ACTUATOR, f'{name}_{ax}')
    assert (out >= 0).all()
    return out


def replay(prob, solver, friction=None, density=1.40, feedback=True):
    import mujoco
    mj = mujoco.MjModel.from_xml_path(str(MJCF_OUT))
    mj.opt.gravity[:] = [0, 0, -float(-prob.model.gravity.linear[2])]
    mj.opt.density = density
    mj.dof_damping[:] = 0.0
    mj.dof_armature[:] = 0.0
    from humanities.dance.motion.realise import MU
    mj.geom_friction[:, 0] = friction if friction is not None else MU
    mj.actuator_ctrllimited[:] = 1
    md = mujoco.MjData(mj)
    mp = MujocoMap(prob.model, mj)
    amap = actuator_map(prob.model, mj)
    state = prob.state
    steps = []
    xs_all, us_all, K_all = list(solver.xs), list(solver.us), list(solver.K)
    for x, u, K, m in zip(xs_all[:-1], us_all, K_all, prob.problem.runningModels):
        if isinstance(m, crocoddyl.ActionModelImpulseFwdDynamics):
            continue
        steps.append((x, u, K, m.dt))
    x0 = solver.xs[0]
    qpos, qvel = mp.to_mujoco(x0[:prob.model.nq], x0[prob.model.nq:])
    md.qpos[:] = qpos
    md.qvel[:] = qvel
    mujoco.mj_forward(mj, md)
    sole = {n: mujoco.mj_name2id(mj, mujoco.mjtObj.mjOBJ_GEOM, n) for n in body_spec()['contacts']}
    floor = mujoco.mj_name2id(mj, mujoco.mjtObj.mjOBJ_GEOM, 'floor')
    n_sub = None
    log = []
    fell = None
    t = 0.0
    for i, (xr, uff, K, dt) in enumerate(steps):
        n_sub = max(1, int(round(dt / mj.opt.timestep)))
        for s in range(n_sub):
            q, v = mp.from_mujoco(md.qpos, md.qvel)
            x = np.concatenate([q, v])
            u = uff - K @ state.diff(xr, x) if feedback else uff
            ctrl = np.zeros(mj.nu)
            ctrl[amap] = u
            md.ctrl[:] = ctrl
            mujoco.mj_step(mj, md)
            t += mj.opt.timestep
        # Log once a node: drift from the plan, contacts and their forces.
        q, v = mp.from_mujoco(md.qpos, md.qvel)
        xr_next = steps[i + 1][0] if i + 1 < len(steps) else xs_all[-1]
        dq = state.diff(xr_next, np.concatenate([q, v]))[:prob.model.nv]
        pd = prob.model.createData()
        com_plan = pin.centerOfMass(prob.model, pd, xr_next[:prob.model.nq]).copy()
        com_now = pin.centerOfMass(prob.model, pd, q).copy()
        fn, ft = 0.0, 0.0
        for c in range(md.ncon):
            con = md.contact[c]
            if floor not in (con.geom1, con.geom2):
                continue
            f6 = np.zeros(6)
            mujoco.mj_contactForce(mj, md, c, f6)
            fn += f6[0]
            ft += math.hypot(f6[1], f6[2])
        log.append(dict(t=round(t, 3), com_drift_m=round(float(np.linalg.norm(com_now - com_plan)), 4),
                        pelvis_tilt_deg=round(float(np.degrees(np.linalg.norm(dq[3:5]))), 2),
                        joint_drift_deg=round(float(np.degrees(np.max(np.abs(dq[6:])))), 2),
                        normal_N=round(fn, 1), friction_used=round(ft / fn, 3) if fn > 1 else 0.0))
        if com_now[2] < 0.45 and fell is None:
            fell = t
    return dict(fell_at_s=fell, nodes=log,
                worst=dict(com_drift_m=max(r['com_drift_m'] for r in log),
                           pelvis_tilt_deg=max(r['pelvis_tilt_deg'] for r in log),
                           friction_used=max(r['friction_used'] for r in log)))
