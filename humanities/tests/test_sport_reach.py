"""The reach engine follows the envelope's closed forms and the contract's laws for both bodies on both worlds, over the
grips where grip matters, and its single-point, earliest-time and vectorised answers agree."""
import math
import time

import numpy as np
import pytest

from humanities.sport import envelope as env
from humanities.sport import reach
from humanities.sport.state import (BODIES, G_E, G_M, GRIPS, SPRINT_F0_N_KG, SPRINT_V0_M_S, UNITS, V_MAX_M_S, WORLDS,
                                    State, muscle_push)

CASES = [(b, w, g) for b in BODIES.values() for w, g in WORLDS]
ENV = {b.key: env.bodies()[b.name] for b in BODIES.values()}


def sprint(t, mu, g, va=0.0, ns=0.0):
    """The moving push in closed form: x'' = min(mu g, F0 (1 - (va + x')/V0)) from rest at ns, the speed along u held
    at or below V_MAX_M_S: the grip's parabola until the muscles take over, then the sprint profile's exponential."""
    F0, V0 = SPRINT_F0_N_KG, SPRINT_V0_M_S
    r, a, cap = F0 / V0, mu * g, V_MAX_M_S - va
    s = np.maximum(np.asarray(t, float) - ns, 0.0)
    u1 = min(max(V0 * (1 - a / F0) - va, 0.0), cap)
    t1, S0 = u1 / a, va + u1

    def muscles(s2):
        e = np.exp(-r * s2)
        return 0.5 * a * t1 ** 2 + (V0 - va) * s2 - (V0 - S0) / r * (1 - e), V0 - (V0 - S0) * e - va

    if u1 >= cap:
        tc, xc = cap / a, 0.5 * cap * cap / a
    else:
        s2c = math.log((V0 - S0) / (V0 - V_MAX_M_S)) / r
        tc, xc = t1 + s2c, muscles(s2c)[0]
    xm, um = muscles(np.maximum(s - t1, 0.0))
    x, u = np.where(s <= t1, 0.5 * a * s * s, xm), np.where(s <= t1, a * s, um)
    return np.where(s > tc, xc + cap * (s - tc), x), np.where(s > tc, cap, u)


def test_the_support_reaches_behind_the_centre_of_mass_by_posture():
    for b in BODIES.values():
        assert reach.lead0(b, 'upright', 0.0, 1.0) == pytest.approx(b.side_reach_m)
        assert reach.lead0(b, 'upright', -1.0, 0.0) == pytest.approx(b.fore_aft_reach_m)
        assert reach.lead0(b, 'ready', 0.0, -1.0) == pytest.approx(b.ready_half_stance_m)
        assert reach.lead0(b, 'split', 1.0, 0.0) == pytest.approx(b.split_reach_m)
        assert reach.lead0(b, 'split', -1.0, 0.0) == pytest.approx(b.fore_aft_reach_m)
        assert reach.lead0(b, 'split', 0.0, 1.0, facing=math.pi / 2) == pytest.approx(b.split_reach_m)


def test_the_legs_push_is_the_contracts_sprint_profile():
    for s in (-2.0, 0.0, 3.0, 7.0, 9.0):
        assert reach.muscle(s) == pytest.approx(muscle_push(s))


@pytest.mark.parametrize('grip', GRIPS)
def test_an_upright_lean_sideways_follows_lean_start_until_the_grip_binds(grip):
    t = np.arange(0.0, 4.0, 0.001)
    for b, w, g in CASES:
        h, L = b.com_height_m, b.side_reach_m
        _, x, v = reach.profile(b, 'upright', math.pi / 2, grip, g, t=t, steps=False)
        for d in (0.1, 0.3):
            assert L + d < grip * h                      # still leaning at these distances on every grip
            ref = env.lean_start(d, g, h, L)
            assert np.interp(d, x, t) == pytest.approx(ref['seconds'], abs=1.5e-3)
            assert np.interp(d, x, v) == pytest.approx(ref['speed_m_s'], abs=6e-3)
        # past the lean, the push is the grip's or the muscles' (the speed's own law); the drop step only adds
        a = np.gradient(v, t)
        after = (x > grip * h - L + 0.05) & (v < V_MAX_M_S - 0.05) & (t < t[-2])
        assert after.any()
        assert np.allclose(a[after], np.minimum(grip * g, reach.muscle(v[after])), rtol=0.01, atol=1e-3)
        _, xs, _ = reach.profile(b, 'upright', math.pi / 2, grip, g, t=t)
        assert np.all(xs >= x - 1e-9)
    # on the Moon the drop step reaches the first 0.1 m sideways before the lean alone; on Earth it does not
    b = BODIES['dancer']
    for g, quicker in ((G_M, True), (G_E, False)):
        _, lean, _ = reach.profile(b, 'upright', math.pi / 2, 0.8, g, t=t, steps=False)
        _, best, _ = reach.profile(b, 'upright', math.pi / 2, 0.8, g, t=t)
        assert (np.interp(0.1, best, t) < np.interp(0.1, lean, t) - 0.1) == quicker


@pytest.mark.parametrize('grip', GRIPS)
def test_a_moving_player_pushes_at_the_grip_or_the_muscles_up_to_top_speed(grip):
    t = np.arange(0.0, 12.0, 0.005)
    for b, w, g in CASES:
        for va in (0.0, 3.0):
            for ns in (0.0, 0.2):
                _, x, v = reach.profile(b, 'moving', 0.3, grip, g, v_along=va, next_stance_s=ns, t=t)
                xr, vr = sprint(t, grip, g, va, ns)
                # the speed is read linearly between the profile's steps, rounding the cap's corner by mm/s
                assert np.allclose(x, xr, atol=2e-4) and np.allclose(v, vr, atol=5e-3)
                assert v.max() <= V_MAX_M_S - va + 1e-9 and np.all(x[t <= ns] == 0)
                assert v[-1] == pytest.approx(V_MAX_M_S - va)
        # the Moon at low speed: the grip binds, x = mu g t^2 / 2; on Earth the muscles bind from the start at 0.8 up
        _, x, _ = reach.profile(b, 'moving', 0.0, grip, g, t=t)
        one = np.interp(1.0, t, x)
        if g == G_M:
            assert one == pytest.approx(0.5 * grip * g, rel=1e-6)
        elif grip >= 0.8:
            assert one < 0.75 * 0.5 * grip * g


@pytest.mark.parametrize('grip', GRIPS)
def test_the_moving_reach_set_is_a_circle_about_the_drift_point(grip):
    for b in BODIES.values():
        st = State(p=(2.0, -1.0), v=(1.0, 0.5), posture='moving')
        for t in (1.0, 2.0):                             # the Moon's grip binds in every direction below 6 m/s
            ring = reach.reach_set(b, st, t, grip, G_M)
            r = np.hypot(ring[:, 0] - 2.0 - 1.0 * t, ring[:, 1] + 1.0 - 0.5 * t)
            assert np.allclose(r, 0.5 * grip * G_M * t * t, rtol=1e-9)
        for t in (0.5, 1.0):                             # Earth without drift: the muscles' circle
            ring = reach.reach_set(b, State(p=(2.0, -1.0), posture='moving'), t, grip, G_E)
            assert np.allclose(np.hypot(ring[:, 0] - 2.0, ring[:, 1] + 1.0), sprint(t, grip, G_E)[0], rtol=1e-6)
        # the cap, and on Earth the muscles' speed, shorten the reach along the drift and lengthen it against
        for t, g in ((10.0, G_M), (1.0, G_E)):
            ring = reach.reach_set(b, State(v=(3.0, 0.0), posture='moving'), t, grip, g, n=4)
            r = np.hypot(ring[:, 0] - 3.0 * t, ring[:, 1])
            assert r[0] < r[1] < r[2]


def test_a_jump_takes_off_on_the_contracts_timing_and_rides_its_arc():
    for b, w, g in CASES:
        full = b.take_off_m_s(g)
        assert full == pytest.approx(env.take_off(ENV[b.key], g))
        z_to = reach.take_off_height(b)
        assert z_to == pytest.approx(b.com_height_m + b.push_off_m - b.dip_m)
        for posture in ('upright', 'ready', 'split', 'moving'):
            for share in (0.3, 0.6, 0.85, 1.0):
                vz = share * full
                dip, dip_s, push_s = b.jump_timing(posture, vz, g)
                assert reach.prep_s(b, posture, vz, g) == pytest.approx(dip_s + push_s, rel=1e-9)
                j = reach.jumps(b, posture, g, vz=vz)
                assert j['take_off_s'] == pytest.approx(dip_s + push_s) and j['dip_m'] == pytest.approx(dip)
                # back at the take-off height after 2 vz / g, then down to the standing height: State.landing_s
                drop_s = (math.sqrt(vz * vz + 2 * g * (z_to - b.com_height_m)) - vz) / g
                assert j['airtime_s'] == pytest.approx(2 * vz / g + drop_s)
                assert State(posture='aloft', z=z_to, vz=vz).landing_s(b, g) == pytest.approx(j['airtime_s'])
                assert j['hands_top_m'] == pytest.approx(b.hand_band(z_to + vz * vz / (2 * g))[1])
        # the full jump from upright is the envelope's countermovement: the full dip, then 2 push_off / vz
        j = reach.jumps(b, 'upright', g)
        assert j['take_off_s'] == pytest.approx(env.dip_s(b.dip_m, g) + 2 * b.push_off_m / full)
        assert reach.jumps(b, 'moving', g, next_stance_s=0.3)['take_off_s'] == pytest.approx(0.3 + env.STRIDE_S)
        # the hands at the apex of the full jump: a point 2 cm under their top is in reach by 2 cm, 3 cm over is not
        t = j['apex_s']
        r = reach.reach_at(b, State(p=(1.0, 1.0), posture='upright'), (1.0, 1.0, j['hands_top_m'] - 0.02), t, 0.8, g)
        assert r.ok and r.mode == 'jump' and r.margin_m == pytest.approx(0.02, abs=1e-6)
        assert r.take_off_s == pytest.approx(j['take_off_s'], abs=1e-6) and r.moves == ('standing_jump',)
        assert r.committed_s == pytest.approx(full / g) and r.reads_s == pytest.approx(full / g)
        r = reach.reach_at(b, State(p=(1.0, 1.0), posture='upright'), (1.0, 1.0, j['hands_top_m'] + 0.03), t, 0.8, g)
        assert not r.ok and r.margin_m == pytest.approx(-0.03, abs=1e-6)
    # a small jump tells far less than the full one: a 1 m lunar rise leaves upright in about 0.4 s, against 0.87 s
    b = BODIES['dancer']
    small = reach.jumps(b, 'upright', G_M, vz=math.sqrt(2 * G_M * 1.0))
    assert 0.3 < small['take_off_s'] < 0.45 and small['dip_m'] == pytest.approx(0.05)
    assert reach.jumps(b, 'upright', G_M)['take_off_s'] == pytest.approx(0.865, abs=2e-3)
    # the engine finds a way overhead to the small jump's hand height no later than that jump's apex
    r = reach.reach_time(b, State(posture='upright'), (0.0, 0.0, small['hands_top_m']), 0.8, G_M)
    assert r.ok and r.mode == 'jump' and r.seconds <= small['apex_s'] + 0.01 and r.moves == ('standing_jump',)
    assert r.take_off_s + r.committed_s == pytest.approx(r.seconds)


@pytest.mark.parametrize('grip', GRIPS)
def test_a_side_hop_pushes_sideways_at_the_stances_ratio(grip):
    for b in BODIES.values():
        k = float(reach.push_ratio(b, 'ready', 0.0, 1.0, grip))
        assert k == pytest.approx(min(grip, b.ready_half_stance_m / b.ready_com_m))
        assert float(reach.push_ratio(b, 'moving', 0.0, 1.0, grip)) == grip
        hop = env.jumps(ENV[b.key])['lateral_hop']
        if grip == env.MU:
            assert k == pytest.approx(hop['push_ratio'], abs=1e-3)
            for row in hop['hops']:
                vz = row['vertical_m_s']
                for w, g in WORLDS:
                    assert k * vz * 2 * vz / g == pytest.approx(row[w]['covered_m'], abs=2e-3)
        for w, g in WORLDS:
            # the crouch's best hop without a dip, read at its apex: the engine reaches at least that far sideways
            vz = b.take_off_m_s(g, b.crouch_push_m)
            j = reach.jumps(b, 'ready', g, vz=vz)
            assert j['dip_m'] == pytest.approx(0.0, abs=1e-12)
            side, z = k * vz * vz / g + b.reach_out_m - 0.01, reach.take_off_height(b) + j['rise_m']
            r = reach.reach_at(b, State(posture='ready'), (0.0, side, z), j['apex_s'], grip, g)
            assert r.ok and r.margin_m >= 0.01 - 1e-6 and set(r.moves) <= set(UNITS)
            if z > b.reach_up_m:                         # the Moon's hop rises above the standing reach
                assert r.mode == 'jump' and r.moves[-1] in ('side_hop', 'crouch_jump')


@pytest.mark.parametrize('grip', GRIPS)
def test_a_landing_changes_the_travel_by_its_grip_on_the_vertical_impulse(grip):
    for b, w, g in CASES:
        drop = 2.5 ** 2 / (2 * g)                        # lands at 2.5 m/s down, travelling at 4 m/s
        st = State(v=(4.0, 0.0), posture='aloft', z=b.com_height_m + drop, vz=0.0)
        L = reach.landing(b, st, grip, g)
        ref = env.braking_landing(4.0, 2.5, grip, g)
        assert L['v_down_m_s'] == pytest.approx(2.5) and L['stance_s'] == pytest.approx(ref['stance_s'], abs=6e-4)
        assert L['change_m_s'] == pytest.approx(grip * (2.5 + g * L['stance_s']))
        assert min(4.0, L['change_m_s']) == pytest.approx(ref['shed_in_landing_m_s'], abs=6e-4)
        T = L['landing_s'] + L['stance_s']
        t = np.arange(0.0, T + 6.0, 0.001)
        _, x, v = reach.profile(b, 'aloft', math.pi, grip, g, v_along=-4.0, z=st.z, vz=0.0, t=t)
        assert np.all(x[t <= L['landing_s']] == 0)
        _, xT, vT = reach.profile(b, 'aloft', math.pi, grip, g, v_along=-4.0, z=st.z, vz=0.0, t=np.array([T]))
        assert vT[0] == pytest.approx(L['change_m_s'], rel=1e-9)
        assert xT[0] == pytest.approx(0.5 * L['change_m_s'] * L['stance_s'], rel=1e-9)
        if g == G_M:                                     # the grip brakes the rest: the envelope's stop
            stop = np.interp(4.0, v, t)
            assert 4.0 * (stop - L['landing_s']) - np.interp(stop, t, x) == pytest.approx(ref['total_m'], abs=2e-3)


def test_a_fallen_player_gets_up_before_anything_else():
    for b, w, g in CASES:
        up = reach.get_up_s(b, g)
        assert up == pytest.approx(env.floor_and_falls(ENV[b.key])['stand_up_and_stop_s'][w]['0.45_m'], abs=6e-4)
        point = (1.5, 0.5, 1.2)
        down = reach.reach_time(b, State(posture='down'), point, 0.8, g)
        ready = reach.reach_time(b, State(posture='ready'), point, 0.8, g)
        assert down.seconds == pytest.approx(up + ready.seconds, abs=3e-3)
        assert down.moves[0] == 'fall_and_rise' and down.mode == ready.mode
        # before it is up, the hands reach only from the floor to the standing band's low end
        assert reach.reach_at(b, State(posture='down'), (0.3, 0.0, 0.2), 0.5 * up, 0.8, g).ok
        assert not reach.reach_at(b, State(posture='down'), (0.3, 0.0, 1.0), 0.5 * up, 0.8, g).ok


def test_aloft_the_hands_ride_the_arc_until_the_landing():
    for b, w, g in CASES:
        st = State(p=(1.0, 2.0), v=(2.0, -1.0), posture='aloft', z=b.com_height_m + 0.5, vz=1.0)
        T = st.landing_s(b, g)
        for s in (0.25, 0.5, 0.9):
            t = s * T
            zc = st.z + st.vz * t - 0.5 * g * t * t
            cx, cy = 1.0 + 2.0 * t, 2.0 - 1.0 * t
            r = reach.reach_at(b, st, (cx + 0.5 * b.reach_out_m, cy, b.hand_band(zc)[1] - 0.01), t, 0.8, g)
            assert r.ok and r.mode == 'aloft' and r.margin_m == pytest.approx(0.01, abs=1e-6)
            assert r.committed_s == pytest.approx(t) and r.reads_s == pytest.approx(t)
            r = reach.reach_at(b, st, (cx, cy + b.reach_out_m + 0.05, zc), t, 0.8, g)
            assert not r.ok and r.margin_m == pytest.approx(-0.05, abs=1e-6)
        # after the landing, the floor again
        r = reach.reach_at(b, st, (1.0 + 2.0 * (T + 1.0), 2.0 - (T + 1.0), 1.2), T + 1.0, 0.8, g)
        assert r.ok and r.mode == 'landing' and r.moves[0] == 'braking_landing'


def test_reach_at_reach_along_and_reach_time_agree():
    rng = np.random.default_rng(7)
    states = (State(posture='ready', facing=0.4), State(p=(1.0, -1.0), v=(1.5, 0.5), posture='moving',
                                                         next_stance_s=0.1),
              State(posture='upright'), State(posture='split', facing=-1.0))
    for b, w, g in CASES:
        for st in states:
            pts = np.c_[rng.uniform(-6, 6, (12, 2)), rng.uniform(0.0, 4.5, 12)]
            ts = rng.uniform(0.0, 4.0, 12)
            ok, m, c = reach.reach_along(b, st, ts, pts, 0.8, g)
            for i in range(12):
                r = reach.reach_at(b, st, pts[i], ts[i], 0.8, g)
                assert r.ok == ok[i] and r.mode == reach.MODES[c[i]] and r.margin_m == pytest.approx(m[i], abs=1e-9)
                assert set(r.moves) <= set(UNITS)
            for pt in pts[:3]:
                r = reach.reach_time(b, st, pt, 0.8, g)
                if not r.ok:
                    continue
                assert reach.reach_at(b, st, pt, r.seconds, 0.8, g).ok
                early = np.arange(0.0, r.seconds - 1e-3, 0.002)
                assert not reach.reach_along(b, st, early, np.tile(pt, (early.size, 1)), 0.8, g)[0].any()


def test_reach_along_answers_ten_thousand_points_in_well_under_a_second():
    rng = np.random.default_rng(3)
    pts = np.c_[rng.uniform(-8, 8, (10000, 2)), rng.uniform(0.0, 4.0, 10000)]
    ts = rng.uniform(0.0, 6.0, 10000)
    for st in (State(posture='ready', facing=0.7, v=(0.3, 0.1)), State(posture='moving', v=(2.5, -0.5))):
        start = time.perf_counter()                      # a fresh state: the table is built inside the timing
        ok, m, c = reach.reach_along(BODIES['larger'], st, ts, pts, 0.8, G_M)
        seconds = time.perf_counter() - start
        assert seconds < 1.0, seconds
        assert ok.shape == m.shape == c.shape == (10000,) and ok.dtype == bool and set(np.unique(c)) <= {0, 1}
        assert 0.1 < ok.mean() < 0.9
