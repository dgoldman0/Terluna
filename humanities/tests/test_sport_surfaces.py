"""The surface moves keep their physics: pushes follow the force-velocity profile without gravity, every contact stays
inside its friction cone and balances its momentum, the arms hold what they pull, gravity's clock sets the Moon-to-Earth
ratios, closed forms match direct integration, and floor travel is the reach engine's."""
import math

import numpy as np
import pytest

from humanities.dance import envelope as dance
from humanities.sport import envelope as env
from humanities.sport import reach
from humanities.sport import surfaces as S
from humanities.sport.state import BODIES, G_E, G_M, GRIPS, WORLDS, State

CASES = [(b, w, g) for b in BODIES.values() for w, g in WORLDS]


def test_a_push_off_a_surface_follows_the_profile_without_gravity():
    for b in BODIES.values():
        v = S.push_speed(b)
        assert v == pytest.approx(b.take_off_m_s(0.0))
        assert v * v == pytest.approx(2 * b.push_off_m * dance.F0_N_KG * (1 - v / (2 * dance.V0_M_S)))
        assert all(v > b.take_off_m_s(g) for _, g in WORLDS)
        assert S.push_from(b, 0.0, b.push_off_m) == pytest.approx(v)
        for u in (0.5, 1.5, 2.5):
            d = S.push_from(b, u, 0.3)
            assert (u + d) ** 2 - u * u == pytest.approx(2 * 0.3 * dance.F0_N_KG * (1 - (u + d / 2) / dance.V0_M_S))
        assert S.push_from(b, dance.V0_M_S, 0.3) == 0.0
        # the push's mean acceleration carries the body to v over the very stroke the profile needs for v
        T, _ = S.contact_s(b, 0.0, v)
        assert 0.5 * v * T == pytest.approx(b.push_off_m)


@pytest.mark.parametrize('wall_grip', GRIPS)
def test_a_wall_push_stays_in_its_friction_cone_and_balances_momentum(wall_grip):
    for b, w, g in CASES:
        for speed in (2.0, 3.0, 4.0):
            for angle in (30.0, 45.0, 60.0):
                t = S.tic_tac(b, speed, angle, g, wall_grip)
                vn, vt = speed * math.cos(math.radians(angle)), speed * math.sin(math.radians(angle))
                J = vn + t['push_out_m_s']
                assert t['normal_impulse_m_s'] == pytest.approx(J)
                # the friction along the wall and the hold against gravity share the cone mu J; where the cone cannot
                # hold the height through the contact, the body sinks by the rest and nothing is left along the wall
                assert t['holds_height_m_s'] == pytest.approx(g * t['contact_s'])
                if t['sinks_m_s'] > 0:
                    assert t['along_wall_m_s'] == 0.0
                    assert t['holds_height_m_s'] - t['sinks_m_s'] == pytest.approx(wall_grip * J)
                else:
                    assert math.hypot(t['along_wall_m_s'], t['holds_height_m_s']) == pytest.approx(wall_grip * J)
                # the impulse is the change of velocity: in (vt, -vn) toward the wall, out (vt - side, +push)
                out = np.array([vt - t['along_wall_m_s'], t['push_out_m_s']])
                dv = out - np.array([vt, -vn])
                assert np.hypot(*dv) == pytest.approx(t['most']['velocity_change_m_s'])
                assert t['most']['turn_deg'] >= t['kept']['turn_deg'] - 1e-9 and t['yield_ok']
                # the contact stepped in time: the normal speed taken up, then pushed out, at the two accelerations
                a_push = dance.F0_N_KG * (1 - t['push_out_m_s'] / (2 * dance.V0_M_S))
                u, s, dt = -vn, 0.0, 1e-5
                while u < t['push_out_m_s']:
                    u += (t['yield_m_s2'] if u < 0 else a_push) * dt
                    s += dt
                assert s == pytest.approx(t['contact_s'], abs=1e-3)
    # the key contrast: on the Moon a wall's push changes the velocity about twelve times what a floor stance can,
    # and the floor turns the same way only through a bounce of seconds; on Earth the floor is a match for it
    b = BODIES['dancer']
    moon, earth = S.tic_tac(b, 3.0, 45.0, G_M), S.tic_tac(b, 3.0, 45.0, G_E)
    assert moon['kept']['velocity_change_m_s'] > 10 * moon['floor_cut']['one_stance_m_s']
    assert moon['floor_cut']['total_s'] > 10 * moon['contact_s']
    assert earth['kept']['velocity_change_m_s'] < 2.5 * earth['floor_cut']['one_stance_m_s']


def test_the_floor_cut_is_the_envelopes():
    for _, g in WORLDS:
        for speed in (2.0, 3.0, 4.0):
            for turn in (30, 60, 90, 120):
                ref = env.cut_bounce(speed, turn, 0.8, g)
                out = S.floor_cut(speed, turn, 0.8, g)
                assert out['bounce_m_s'] == pytest.approx(ref['vertical_speed_each_way_m_s'], abs=6e-4)
                assert out['velocity_change_m_s'] == pytest.approx(ref['velocity_change_m_s'], abs=6e-4)
                if out['bounce_m_s'] > 0:
                    assert 0.8 * (2 * out['bounce_m_s'] + g * env.STRIDE_S) == pytest.approx(out['velocity_change_m_s'])


def test_a_wall_run_keeps_each_step_in_its_cone_and_the_wall_within_reach():
    for b, w, g in CASES:
        reach_m = S.WALL_REACH_M * b.stature_m / env.STATURE_M
        for approach in (3.0, 4.0, 5.0):
            runs = [S.wall_run(b, approach, g, steps=k) for k in (1, 2, 3)]
            for k, run in zip((1, 2, 3), runs):
                d = run['push_each_m_s']
                Js = [approach + d] + [d] * (k - 1)
                assert run['lifts_m_s'] == pytest.approx([0.8 * J for J in Js])
                for i in range(1, k):                 # each later step reaches the wall and the legs can push
                    off = d * S.STEP_S * i * (i + 1) / 2
                    assert off <= reach_m + 1e-9
                    assert d <= S.push_from(b, i * d, min(b.push_off_m, reach_m - off)) + 1e-9
                # stepped in time: the lifts at each step's start, gravity throughout
                z = vz = 0.0
                for i, lift in enumerate(run['lifts_m_s']):
                    vz += lift
                    span = max(S.STEP_S, S.contact_s(b, approach, d)[0]) if i == 0 else S.STEP_S
                    for _ in range(2000):
                        z += vz * span / 2000 - 0.5 * g * (span / 2000) ** 2
                        vz -= g * span / 2000
                assert z == pytest.approx(run['rise_on_wall_m'], abs=1e-6)
            # a flat wall throws the runner off: more steps keep it closer and lift it less
            assert runs[0]['rise_m'] > runs[1]['rise_m'] > runs[2]['rise_m']
            assert runs[0]['outward_off_m_s'] > runs[2]['outward_off_m_s']
    d = BODIES['dancer']
    assert S.wall_run(d, 4.0, G_M)['rise_m'] > 5 * S.wall_run(d, 4.0, G_E)['rise_m']
    # a curved wall holds a runner by v^2/R; Minetti et al. (2024) give 3.5 m/s for their 4.7 m radius on the Moon
    v = S.wall_of_death_m_s(G_M, 4.7, 0.62)
    assert v == pytest.approx(3.5, abs=0.05) and 0.62 * v * v / 4.7 == pytest.approx(G_M)
    assert S.wall_of_death_m_s(G_E, 4.7) > 7.0 > S.wall_of_death_m_s(G_M, 4.7)


def test_hands_on_an_edge_hold_what_the_arms_pull():
    for b, w, g in CASES:
        for bw in S.ARM_PULL_BRACKET + (S.ARM_PULL_BW,):
            c = S.cat_leap(b, g, arm_bw=bw)
            if c['fastest_catch_m_s'] > 0:
                at = S.cat_leap(b, g, v_down=c['fastest_catch_m_s'], arm_bw=bw)
                assert at['arms_need_bw'] == pytest.approx(bw) and at['holds']
                assert not S.cat_leap(b, g, v_down=1.01 * c['fastest_catch_m_s'], arm_bw=bw)['holds']
            up = S.climb_up(b, g, bw)
            a = bw * G_E - g
            if up['ok'] and up['how'] == 'arms':
                assert up['seconds'] == pytest.approx(env.rise_and_stop_s(up['rise_m'], a, g))
                assert up['rise_m'] * g / (a + g) <= 2 * S.ARM_STROKE_M * b.stature_m / env.STATURE_M + 1e-9
            pd = S.pull_down(b, 0.5, g, bw)
            # stepped in time: falling at g plus the pull until the legs' brake just stops the body at the depth
            acc, z, v, t, dt = g + bw * G_E, 0.0, 0.0, 0.0, 1e-6
            while v * v / (2 * env.BRAKE_M_S2) < 0.5 - z:
                v += acc * dt
                z += v * dt
                t += dt
            assert t + v / env.BRAKE_M_S2 == pytest.approx(pd['seconds'], rel=2e-3)
            assert pd['seconds'] < pd['gravity_alone_s']
    assert not S.climb_up(BODIES['dancer'], G_E, 1.0)['ok'] and S.climb_up(BODIES['dancer'], G_M, 1.0)['ok']
    d = BODIES['dancer']
    assert S.pull_down(d, 0.5, G_M)['saves_s'] > 4 * S.pull_down(d, 0.5, G_E)['saves_s']


def test_a_vault_clears_its_box_and_keeps_the_run():
    for b, w, g in CASES:
        for height in (1.0, 1.2):
            for speed in (2.0, 3.0, 4.0):
                for kind in ('vault', 'kong_vault'):
                    out = S.vault(b, height, 1.5, speed, g, kind)
                    top = height + S.VAULT_CLEAR_M[kind] * b.stature_m / env.STATURE_M
                    hold = out['hands_hold_share'] * g
                    z, vz, x = reach.take_off_height(b), out['take_off_m_s'], -out['take_off_before_m']
                    t, dt = 0.0, 1e-4
                    lowest = math.inf
                    while z >= b.com_height_m or t == 0.0:    # stepped in time through the three phases
                        over = 0.0 <= x <= out['length_m']
                        vz -= (g - hold if over else g) * dt
                        z += vz * dt
                        x += speed * dt
                        t += dt
                        if over:
                            lowest = min(lowest, z)
                    assert lowest >= top - 2e-3
                    assert t == pytest.approx(out['flight_s'], abs=2e-3)
                    assert out['speed_kept_m_s'] == pytest.approx(speed + out['hands_forward_m_s'])
                    assert out['ok'] == (out['take_off_m_s'] <= b.take_off_m_s(g))
                    assert out['take_off_m_s'] <= out['without_hands_m_s'] + 1e-12
        # a 1.0 m box is vaulted on both worlds; with these clearances a 1.2 m box only on the Moon
        assert S.vault(b, 1.0, 1.5, 3.0, g, 'vault')['ok']
        assert S.vault(b, 1.2, 1.5, 3.0, g, 'vault')['ok'] == (g == G_M)


def test_a_drop_lands_at_its_fall_speed_and_a_roll_spreads_it():
    for b, w, g in CASES:
        for level in (1.0, 3.0, 6.0):
            out = S.drop_roll(b, level, g, speed=2.0)
            assert out['landing_m_s'] == pytest.approx(math.sqrt(2 * g * level))
            assert out['fall_s'] * g == pytest.approx(out['landing_m_s'])
            assert out['legs_push_earth_bw'] == pytest.approx((2 * g * level / (2 * env.STROKE_M) + g) / G_E)
            assert out['roll_push_earth_bw'] < out['legs_push_earth_bw']
            ref = env.braking_landing(2.0, out['landing_m_s'], 0.8, g)
            assert out['braking_landing']['shed_m_s'] == pytest.approx(ref['shed_in_landing_m_s'], abs=6e-4)
    # gravity's clock: a 6 m drop on the Moon lands as fast as a 1 m drop on Earth, and its stroke brakes as hard;
    # only the body's own weight on top differs
    d = BODIES['dancer']
    moon, earth = S.drop_roll(d, 6.0, G_M), S.drop_roll(d, 6.0 * G_M / G_E, G_E)
    assert moon['landing_m_s'] == pytest.approx(earth['landing_m_s'])
    assert moon['legs_push_earth_bw'] * G_E - G_M == pytest.approx(earth['legs_push_earth_bw'] * G_E - G_E)


def test_a_precision_landing_stops_within_its_lean():
    for b, w, g in CASES:
        for gap in (0.5, 1.0, 2.0):
            out = S.precision_jump(b, gap, g)
            assert out['limited_by'] == 'balance'
            if not out['ok']:
                continue
            T, vz = out['flight_s'], out['take_off_m_s']
            z_to, z_land = reach.take_off_height(b), b.com_height_m
            assert vz * T - 0.5 * g * T * T == pytest.approx(z_land - z_to, abs=1e-9)
            assert out['horizontal_m_s'] * T == pytest.approx(gap)
            room = S.EDGE_W_M / 2 + S.LAND_LEAN_M * b.stature_m / env.STATURE_M
            assert out['sheds_m_s'] >= out['horizontal_m_s'] and out['horizontal_m_s'] * out['stance_s'] / 2 <= room
            assert out['sheds_m_s'] == pytest.approx(room / b.com_height_m * (out['landing_m_s'] + g * out['stance_s']))
    d = BODIES['dancer']
    assert S.precision_jump(d, 1.0, G_E)['ok'] and S.precision_jump(d, 3.0, G_M)['ok']
    assert S.precision_jump(d, 1.0, G_M)['widest_gap_m'] > 3 * S.precision_jump(d, 1.0, G_E)['widest_gap_m']


def test_a_slide_sinks_at_g_and_slides_at_the_clothings_grip():
    for b, w, g in CASES:
        for speed in (2.0, 4.0):
            out = S.slide(b, speed, g)
            sink = b.com_height_m - S.SLIDE_COM_M * b.stature_m / env.STATURE_M
            assert out['sink_s'] == pytest.approx(math.sqrt(2 * sink / g))
            v1 = speed - S.SLIDE_GRIP * math.sqrt(2 * g * sink)
            # stepped in time: sliding at the clothing's grip from the touchdown's speed
            x, v, dt = 0.0, v1, 1e-5
            while v > 0:
                v -= S.SLIDE_GRIP * g * dt
                x += v * dt
            assert x == pytest.approx(out['slide_m'], rel=1e-3, abs=1e-4)
            assert out['slide']['low_s'] < out['run']['low_s']
    d = BODIES['dancer']
    moon = S.slide(d, 3.0, G_M)
    assert moon['run']['stopped_s'] < moon['slide']['stopped_s'] and moon['run']['low_s'] - moon['slide']['low_s'] > 3.0


def test_a_lache_swings_as_a_physical_pendulum():
    for b, w, g in CASES:
        out = S.lache(b, g, amplitude_deg=5.0, pull=False)
        r = S.HANG_BELOW_M * b.stature_m / env.STATURE_M
        I = b.inertia_kg_m2['layout'] + b.mass_kg * r * r
        k = b.mass_kg * g * r / I
        th, w0, t, dt, crossings = math.radians(5.0), 0.0, 0.0, 1e-4, []
        while len(crossings) < 2:                    # stepped in time: theta'' = -(m g r / I) sin theta
            w0 -= k * math.sin(th) * dt
            prev, th = th, th + w0 * dt
            t += dt
            if prev > 0 >= th or prev < 0 <= th:
                crossings.append(t)
        assert 2 * (crossings[1] - crossings[0]) == pytest.approx(out['period_s'], rel=2e-3)
        for amp in (45.0, 90.0):
            free = S.lache(b, g, amplitude_deg=amp, pull=False)
            th_r = math.radians(free['release_deg'])
            # energy: the release speed from the amplitude's height
            wr2 = 2 * k * (math.cos(th_r) - math.cos(math.radians(amp)))
            assert free['release_m_s'] == pytest.approx(math.sqrt(wr2) * r)
            pulled = S.lache(b, g, amplitude_deg=amp)
            if pulled['pulled_in']:                  # the pull-in keeps the angular momentum
                r1 = r - S.ARM_STROKE_M * b.stature_m / env.STATURE_M
                I1 = b.inertia_kg_m2['layout'] + b.mass_kg * r1 * r1
                assert (pulled['after_pull_m_s'] / r1) * I1 == pytest.approx((pulled['bottom_m_s'] / r) * I)
            T = pulled['flight_s']
            vz = pulled['release_m_s'] * math.sin(th_r)
            assert pulled['release_height_m'] + vz * T - 0.5 * g * T * T == pytest.approx(b.com_height_m, abs=1e-9)
    d = BODIES['dancer']
    assert S.lache(d, G_M, amplitude_deg=5)['period_s'] / S.lache(d, G_E, amplitude_deg=5)['period_s'] == pytest.approx(
        math.sqrt(G_E / G_M))


def test_floor_travel_agrees_with_the_reach_engine():
    d = BODIES['dancer']
    for w, g in WORLDS:
        t = np.arange(0.0, 6.0, 0.01)
        for va in (0.0, 2.0):                        # the closed form is the law reach.py integrates
            x, xd = S.moving_x(t, va, 0.8, g)
            _, xr, xdr = reach.profile(d, 'moving', 0.0, 0.8, g, v_along=va, t=t)
            assert np.allclose(x, xr, atol=1e-3) and np.allclose(xd, xdr, atol=6e-3)
        cases = [(State(posture='ready'), (3.0, 0.0)), (State(posture='upright', facing=0.5), (1.0, 2.0)),
                 (State(p=(1.0, 1.0), v=(1.5, 0.5), posture='moving'), (5.0, -1.0)),
                 (State(p=(0.0, 0.0), v=(-1.0, 0.0), posture='moving'), (2.0, 2.0))]
        for st, target in cases:
            low = d.hand_band(d.com_m(st.posture))[0]
            z = 0.5 * (low + d.reach_up_m)
            ref = reach.reach_time(d, st, (target[0], target[1], z), 0.8, g)
            if ref.mode != 'ground':
                continue
            leg, _ = S.floor_leg(d, st, target, 0.8, g, arrive='hands')
            assert leg == pytest.approx(ref.seconds, abs=3e-3)
            route = S.route_time(d, st, [dict(unit='floor', to=target, arrive='hands')], S.Arena(), 0.8, g)
            assert route['ok'] and route['seconds'] == pytest.approx(ref.seconds, abs=3e-3) and route['units'] == []


def test_floor_paths_go_round_the_boxes():
    path = S.floor_path(S.EXAMPLE, S.EXAMPLE_START, S.EXAMPLE_GOAL)
    pts = [S.EXAMPLE_START] + path
    assert all(S.clear(S.EXAMPLE, p, q) for p, q in zip(pts, pts[1:]))
    assert sum(math.dist(p, q) for p, q in zip(pts, pts[1:])) > math.dist(S.EXAMPLE_START, S.EXAMPLE_GOAL)
    assert S.floor_path(S.Arena(), (0, 0), (5, 5)) == [(5.0, 5.0)]
    assert not S.clear(S.EXAMPLE, S.EXAMPLE_START, S.EXAMPLE_GOAL)


def test_the_quickest_route_over_a_barrier_names_its_moves_and_repeats():
    d = BODIES['dancer']
    for w, g in WORLDS:
        q = S.quickest_route(d, State(p=S.BARRIER_START, posture='ready'), S.BARRIER_GOAL, S.BARRIER, 0.8, g)
        best, floor = q['best'], q['floor_only']
        assert best['ok'] and best['seconds'] <= floor['seconds'] and floor['seconds'] > q['open_floor']['seconds']
        assert best['units'] and set(best['units']) <= set(S.PARKOUR_UNITS)
        again = S.route_time(d, State(p=S.BARRIER_START, posture='ready'), best['moves'], S.BARRIER, 0.8, g)
        assert again['seconds'] == pytest.approx(best['seconds'])
        assert set(S.PARKOUR_UNITS) == {'vault', 'kong_vault', 'tic_tac', 'wall_run', 'cat_leap', 'climb_up',
                                        'precision_jump', 'drop_roll', 'lache', 'slide', 'pull_down'}
