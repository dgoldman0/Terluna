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
import crocoddyl


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
    # One landing from the air, on both feet; the heels coming down from the rise are landings of their own.
    from_air = [n for k, n in enumerate(nodes) if n.landing and k > 0 and not nodes[k - 1].contacts]
    assert len(from_air) == 1 and set(from_air[0].landing) == {'L', 'R'}
    heels = [n for k, n in enumerate(nodes) if n.landing and k > 0 and nodes[k - 1].contacts]
    assert all(set(nodes[nodes.index(n) - 1].contacts.values()) == {'ball'} for n in heels)


def test_the_whole_phrase_is_scheduled_on_contacts_that_can_do_it():
    ph = gr.first_phrase()
    nodes = rz.timeline(ph.placed, finish_s=1.0)
    for n in nodes:
        if n.phase == 'turn':
            # A turn stands on one pivoting foot; the other foot touches down only to push.
            assert 'pivot' in n.contacts.values(), n.t
    # Two feet held without swivel never turn the body more than the hips can twist between them.
    run, yaw0 = False, 0.0
    for n in nodes:
        held = len(n.contacts) == 2 and all(k in ('flat', 'ball') for k in n.contacts.values())
        if held and not run:
            yaw0 = n.yaw
        if held:
            assert abs(n.yaw - yaw0) < math.radians(70), n.t
        run = held
    fall = [n for n in nodes if n.unit == 'long_fall' and n.phase == 'fall']
    assert fall[0].contacts == {'L': 'flat'} and fall[-1].contacts == {'L': 'ball'}
    whirl = [n for n in nodes if n.unit == 'whirl' and n.phase == 'turn']
    assert sum(1 for n in whirl if n.landing) == 3          # a push each revolution
    last = nodes[-1]
    assert last.phase == 'hold' and set(last.contacts.values()) == {'flat', 'knee'} and not last.com_vel.any()


def test_footholds_are_planned_for_the_whole_phrase():
    ph = gr.first_phrase()
    nodes = rz.timeline(ph.placed, finish_s=1.0)
    prob = rz.Problem(nodes)
    prob.plan()
    for i, n in enumerate(nodes):
        for fid, M in prob.node_targets[i].items():
            p = M.translation if isinstance(M, pin.SE3) else M
            assert abs(p[2]) < 1e-9 and np.linalg.norm(p[:2] - n.com) < 0.8, (n.t, prob.model.frames[fid].name)
    # The whirl's pivot sits under the centre of mass it turns about.
    i = next(k for k, n in enumerate(nodes) if n.unit == 'whirl' and n.phase == 'turn')
    pivot = next(M for f, M in prob.node_targets[i].items())
    assert np.linalg.norm(pivot[:2] - nodes[i].com) < 0.02


def test_pivot_friction_opposes_the_swivel_with_exact_velocity_derivatives():
    m, _ = build()
    state = crocoddyl.StateMultibody(m)
    fid = m.getFrameId('R_toesole')
    act = rz.PivotFriction(state, [fid])
    rng = np.random.default_rng(5)
    q, v = neutral(m), rng.normal(size=m.nv) * 0.5
    x, u = np.concatenate([q, v]), np.zeros(m.nv - 6)
    d = act.createData()
    act.calc(d, x, u)
    act.calcDiff(d, x, u)
    jz = pin.computeFrameJacobian(m, m.createData(), q, fid, pin.LOCAL_WORLD_ALIGNED)[5]
    w = jz @ v
    assert np.sign(jz @ d.tau) == -np.sign(w)                # the floor's moment works against the swivel
    eps, num = 1e-6, np.zeros((m.nv, m.nv))
    for k in range(m.nv):
        v2 = v.copy()
        v2[k] += eps
        d2 = act.createData()
        act.calc(d2, np.concatenate([q, v2]), u)
        num[:, k] = (d2.tau - d.tau) / eps
    assert np.allclose(d.dtau_dx[:, m.nv:], num, atol=1e-4)


def test_the_posed_body_places_its_joints_as_the_engine_does():
    from humanities.dance.body import pose
    m, _ = build()
    d = m.createData()
    shape = dict(legs='star_air', arms='high_v', torso='spiral', head='lifted')
    node = rz.Node(t=0.0, contacts={}, targets=rz.shape_targets({**rz.CARRIAGE, **shape}, 'R'), com=np.zeros(2))
    q = rz.Problem.__new__(rz.Problem)
    q.model, q.q0 = m, neutral(m)
    qq = rz.Problem.reference_q(q, node)
    qq[:3] = 0.0
    pin.framesForwardKinematics(m, d, qq)
    fr = pose.frames(pose.joint_values(shape, 'R'))
    for seg in ('L_hand', 'R_toes', 'head', 'L_shank'):
        R, p = fr[seg]
        assert np.allclose(d.oMf[m.getFrameId(seg)].translation, p, atol=1e-9), seg
        assert np.allclose(d.oMf[m.getFrameId(seg)].rotation, R, atol=1e-9), seg


def _aerial(sequence):
    nodes = rz.timeline(gr.build('figure', gr.motif(*sequence)).placed, finish_s=1.0)
    return nodes, rz.Problem(nodes)


def test_aerial_figures_plan_one_whole_turn_on_the_momentum_their_shapes_need():
    # A backward hanging salto: the pelvis turns once backward about its lateral axis and lands upright; the planned
    # angular momentum is the whole turn over the integral of dt / I through the figure's shapes.
    nodes, prob = _aerial([('hanging_salto', 'R', {'rise_m': 1.6})])
    (plan,) = prob.air_plans
    flight = [n for n in nodes if n.R_air is not None]
    assert plan['air']['somersault'] == -1.0
    first, last = flight[0].R_air, flight[-1].R_air
    assert np.allclose(first, rz.rz(flight[0].yaw), atol=1e-9)
    assert math.degrees(np.linalg.norm(pin.log3(last.T @ rz.rz(flight[-1].yaw)))) < 10.0
    mid = flight[len(flight) // 2].R_air
    assert mid[2, 2] < -0.9                                     # upside down at the middle of the flight
    data = prob.model.createData()
    integral = 0.0
    for n in flight:
        q = prob.reference_q(n)
        q[:3], q[3:7] = 0.0, [0.0, 0.0, 0.0, 1.0]
        pin.ccrba(prob.model, data, q, np.zeros(prob.model.nv))
        integral += rz.DT / data.Ig.inertia[1, 1]
    L = np.array(plan['angular_momentum_N_m_s'])
    assert L[1] == pytest.approx(-2 * math.pi / integral, abs=0.01) and abs(L[0]) < 1e-9 and abs(L[2]) < 1e-9
    # Within what the envelope says a lunar take-off can give (37 N m s), and the take-off may lean into it.
    assert 5.0 < abs(L[1]) < 37.0
    k = nodes.index(flight[0])
    assert all(n.lean_free for n in nodes[k - int(round(rz.LEAN_FREE_S / rz.DT)):k])


def test_the_corkscrew_twists_once_and_the_flat_spin_lies_flat_about_the_vertical():
    nodes, prob = _aerial([('corkscrew', 'R', {'rise_m': 1.6, 'twists': 1})])
    flight = [n for n in nodes if n.R_air is not None]
    assert prob.air_plans[0]['air']['twist'] == -1.0
    assert math.degrees(np.linalg.norm(pin.log3(flight[-1].R_air.T @ rz.rz(flight[-1].yaw)))) < 10.0
    nodes, prob = _aerial([('flat_spin', 'R', {'rise_m': 2.3, 'speed_m_s': 0.5})])
    (plan,) = prob.air_plans
    L = np.array(plan['angular_momentum_N_m_s'])
    assert abs(L[0]) < 1e-9 and abs(L[1]) < 1e-9 and 10.0 < abs(L[2]) < 25.0     # the grip's 25 N m s at most
    flight = [n for n in nodes if n.R_air is not None]
    mid = flight[len(flight) // 2]
    assert math.degrees(math.acos(mid.R_air[2, 2])) == pytest.approx(rz.FLAT_PITCH_DEG, abs=0.5)
    assert flight[0].R_air[2, 2] > 0.99 and flight[-1].R_air[2, 2] > 0.99


def test_phrases_without_aerial_figures_are_planned_as_before():
    nodes = rz.timeline(gr.first_phrase().placed, finish_s=1.0)
    prob = rz.Problem(nodes)
    assert prob.air_plans == [] and not any(n.lean_free or n.z_rel is not None for n in nodes)


def test_self_distance_measures_the_capsules_it_keeps_apart():
    """The pairs kept apart (the legs, the arms, each arm and the head) are measured as the distance between their
    capsules' surfaces, here at rest against the capsules placed from body.json directly."""
    from humanities.dance.motion.apart import SelfDistance, capsules, closest
    from humanities.dance.motion_io import body_spec
    m, _ = build()
    sd = SelfDistance(m)
    d = sd.distances(neutral(m))[0]
    spec = body_spec()
    origin = {s['name']: np.array(s['origin'], float) for s in spec['segments']}
    caps = capsules(spec)
    for (g, s1, s2), dk in zip(sd.prs, d):
        (a1, b1, r1), (a2, b2, r2) = caps[s1], caps[s2]
        P1, Q1, P2, Q2 = (np.array([v]) for v in (origin[s1] + a1, origin[s1] + b1, origin[s2] + a2, origin[s2] + b2))
        s, t = closest(P1, Q1, P2, Q2)
        ref = np.linalg.norm(P1 + s[:, None] * (Q1 - P1) - P2 - t[:, None] * (Q2 - P2)) - r1 - r2
        assert dk == pytest.approx(ref, abs=1e-9)
    assert len(sd.prs) == 24 and d.min() > 0.0                       # at rest no pair touches


def test_self_distance_derivatives_match_finite_differences():
    """In a random pose the analytic derivative of every pair's distance (each closest point's velocity Jacobian
    along the unit vector between the points) matches finite differences in the configuration's tangent space."""
    from humanities.dance.motion.apart import SelfDistance
    m, _ = build()
    sd = SelfDistance(m)
    rng = np.random.default_rng(11)
    q = pin.integrate(m, neutral(m), rng.normal(size=m.nv) * 0.4)
    d, J = sd.jacobian(q)
    eps, num = 1e-7, np.zeros_like(J)
    for k in range(m.nv):
        e = np.zeros(m.nv)
        e[k] = eps
        num[:, k] = (sd.distances(pin.integrate(m, q, e))[0] - sd.distances(pin.integrate(m, q, -e))[0]) / (2 * eps)
    assert np.allclose(J, num, atol=1e-5)


def test_the_barrier_follows_the_body_between_nodes():
    """With an interval T between nodes the residual also measures each distance on the path from the node before,
    which the integrator's step and the motion file's interpolation both take along the node's velocity: q integrated
    back along v for T/4, T/2 and 3T/4 (a fast shin can slide through a standing one between two nodes). Those are the
    configurations pin.interpolate gives between the node before and this one, and the derivatives in the
    configuration and the velocity match finite differences."""
    from humanities.dance.motion.apart import ResidualApart, SelfDistance
    m, _ = build()
    state = crocoddyl.StateMultibody(m)
    T = 0.025
    res = ResidualApart(state, m.nv - 6, interval=T)
    rng = np.random.default_rng(3)
    q = pin.integrate(m, neutral(m), rng.normal(size=m.nv) * 0.4)
    v = rng.normal(size=m.nv)
    x = np.concatenate([q, v])

    def r(x):
        d = res.createData(crocoddyl.DataCollectorMultibody(m.createData()))
        res.calc(d, x)
        return d.r.copy(), d

    r0, data = r(x)
    res.calcDiff(data, x)
    sd, n = SelfDistance(m), res.n
    q_before = pin.integrate(m, q, -T * v)                       # the node before, which the step took here along v
    assert np.allclose(r0[:n], sd.distances(q)[0])
    assert np.allclose(r0[2 * n:3 * n], sd.distances(pin.interpolate(m, q_before, q, 0.5))[0], atol=1e-9)
    eps, num = 1e-6, np.zeros((len(r0), 2 * m.nv))
    for k in range(2 * m.nv):
        dx = np.zeros(2 * m.nv)
        dx[k] = eps
        num[:, k] = (r(state.integrate(x, dx))[0] - r(state.integrate(x, -dx))[0]) / (2 * eps)
    assert np.allclose(data.Rx, num, atol=1e-5)

