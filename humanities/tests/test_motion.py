"""The motion engine's body agrees with MuJoCo's, its shapes point limbs where their names say, and its timeline
carries the grammar's contact schedule. Runs in the dance environment (Pinocchio, Crocoddyl, MuJoCo); elsewhere it
skips."""
import math

import numpy as np
import pytest

pin = pytest.importorskip('pinocchio')
pytest.importorskip('crocoddyl')
mujoco = pytest.importorskip('mujoco')

from humanities.dance import grammar as gr
from humanities.dance.body.model import MJCF_OUT
from humanities.dance.motion import realise as rz
from humanities.dance.motion.pinmodel import MujocoMap, build, neutral
from humanities.dance.motion.replay import actuator_map


@pytest.fixture(scope='module')
def models():
    m, b = build()
    mj = mujoco.MjModel.from_xml_path(str(MJCF_OUT))
    mj.dof_armature[:] = 0.0
    mj.dof_damping[:] = 0.0
    mj.jnt_limited[:] = 0
    mj.opt.density = 0.0          # no air: the optimiser's dynamics carry none
    mj.opt.viscosity = 0.0
    mj.opt.disableflags |= mujoco.mjtDisableBit.mjDSBL_CONTACT    # random poses can touch limb to limb
    return m, mj


def random_state(m, rng):
    q = neutral(m)
    q[:3] = rng.normal(size=3)
    q[2] = 5.0                    # well above the floor, so no contact enters the comparison
    for j in range(1, m.njoints):
        iq = m.idx_qs[j]
        if m.nqs[j] >= 4:
            o = 3 if m.nqs[j] == 7 else 0
            qq = rng.normal(size=4)
            q[iq + o:iq + o + 4] = qq / np.linalg.norm(qq)
        else:
            q[iq] = rng.uniform(-0.8, 0.8)
    return q, rng.normal(size=m.nv)


def test_pinocchio_and_mujoco_agree_on_kinematics_and_dynamics(models):
    m, mj = models
    d, md = m.createData(), mujoco.MjData(mj)
    mp, amap = MujocoMap(m, mj), actuator_map(m, mj)
    rng = np.random.default_rng(7)
    for _ in range(5):
        q, v = random_state(m, rng)
        pin.forwardKinematics(m, d, q, v)
        pin.updateFramePlacements(m, d)
        qpos, qvel = mp.to_mujoco(q, v)
        md.qpos[:], md.qvel[:] = qpos, qvel
        mujoco.mj_forward(mj, md)
        for name in ('pelvis', 'L_hand', 'R_toes', 'head'):
            bid = mujoco.mj_name2id(mj, mujoco.mjtObj.mjOBJ_XBODY, name)
            assert np.allclose(d.oMf[m.getFrameId(name)].translation, md.xpos[bid], atol=1e-9)
        ke = pin.computeKineticEnergy(m, d, q, v)
        M = np.zeros((mj.nv, mj.nv))
        mujoco.mj_fullM(mj, md, M)
        assert 0.5 * qvel @ M @ qvel == pytest.approx(ke, rel=1e-6)
        # The same torques give the same joint accelerations in free fall without gravity.
        m.gravity.linear[:] = 0.0
        mj.opt.gravity[:] = 0.0
        tau = np.zeros(m.nv)
        tau[6:] = rng.normal(size=m.nv - 6)
        a = pin.aba(m, d, q, v, tau)
        md.ctrl[:] = 0.0
        md.ctrl[amap] = tau[6:]
        mujoco.mj_forward(mj, md)
        _, a_mj = mp.from_mujoco(md.qpos, md.qacc)
        assert np.allclose(a[6:], a_mj[6:], rtol=1e-6, atol=1e-6 * np.abs(a[6:]).max())


def test_ball_targets_point_limbs_where_named():
    down = np.array([0.0, 0.0, -1.0])
    fwd = rz.ball_target(dict(elev=90, az=0), 1).toRotationMatrix() @ down
    assert np.allclose(fwd, [1, 0, 0], atol=1e-9)
    out_left = rz.ball_target(dict(elev=90, az=90), 1).toRotationMatrix() @ down
    out_right = rz.ball_target(dict(elev=90, az=90), -1).toRotationMatrix() @ down
    assert np.allclose(out_left, [0, 1, 0], atol=1e-9) and np.allclose(out_right, [0, -1, 0], atol=1e-9)
    up = rz.ball_target(dict(elev=180, az=0), 1).toRotationMatrix() @ down
    assert np.allclose(up, [0, 0, 1], atol=1e-9)
    # Turnout turns the left foot's toes outward, to the left.
    toes = rz.ball_target(dict(twist=45), 1).toRotationMatrix() @ np.array([1.0, 0, 0])
    assert toes[1] > 0.7


def test_shapes_mirror_between_sides():
    right = rz.shape_targets({'legs': 'arabesque_90'}, 'R')
    left = rz.shape_targets({'legs': 'arabesque_90'}, 'L')
    assert 'R_hip' in right and 'L_hip' in left
    r = right['R_hip'][1].toRotationMatrix() @ np.array([0, 0, -1.0])
    l_ = left['L_hip'][1].toRotationMatrix() @ np.array([0, 0, -1.0])
    assert r[0] < -0.9 and l_[0] < -0.9 and r[1] == pytest.approx(-l_[1], abs=1e-9)


def test_timeline_carries_the_contact_schedule():
    ph = gr.build('t', gr.motif(('plie_rise', 'R', {'depth_m': 0.12}), ('saute', 'R', {'rise_m': 0.5})))
    nodes = rz.timeline(ph.placed)
    flight = [n for n in nodes if not n.contacts]
    assert flight and len(flight) * rz.DT == pytest.approx(gr.flight_s(0.5), abs=2 * rz.DT)
    landings = [n for n in nodes if n.landing]
    assert len(landings) == 1 and set(landings[0].landing) == {'L', 'R'}
