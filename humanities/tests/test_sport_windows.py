"""The ball windows follow the vacuum's closed form, balance drag against weight, match the envelope's best ranges and
throws, hold their limits at both levers, and aim their passes and lobs where they say."""
import math

import numpy as np
import pytest

from humanities.sport import envelope as env
from humanities.sport import reach
from humanities.sport import windows as W
from humanities.sport.state import (BODIES, CATCH_S, G_E, G_M, GRIPS, RELEASE_MAX_M_S, RHO_EARTH, RHO_MOON, WORLDS,
                                    State)

BALL = env.THROW_BALL
AIRS = tuple((w, g, rho) for (w, g), rho in zip(WORLDS, (RHO_MOON, RHO_EARTH)))
ENV = {b.key: env.bodies()[b.name] for b in BODIES.values()}
CHEST = env.AIR_LEVERS_M['chest pass, near the centre of mass']


def test_flights_without_air_follow_the_vacuums_closed_form():
    v0 = np.array([[3.0, 1.0, 4.0], [10.0, 0.0, 10.0], [0.0, -2.0, 6.0], [5.0, 5.0, -1.0]])
    p0 = np.array([1.0, -1.0, 2.0])
    for w, g in WORLDS:
        f = W.fly(BALL, p0, v0, g, 0.0)
        t = f['t']
        for i, v in enumerate(v0):
            T = (v[2] + math.sqrt(v[2] ** 2 + 2 * g * p0[2])) / g
            assert f['land_s'][i] == pytest.approx(T, abs=1e-5)
            live = t < T - 1e-9
            ref = p0 + v * t[:, None] - 0.5 * g * np.outer(t * t, [0.0, 0.0, 1.0])
            assert np.allclose(f['p'][i, live], ref[live], atol=1e-9)
            assert np.all(np.isnan(f['p'][i, t > T + 0.004]))
            assert f['land_p'][i, :2] == pytest.approx(p0[:2] + v[:2] * T, abs=1e-4)
            assert f['land_m_s'][i] == pytest.approx(math.sqrt(v @ v + 2 * g * p0[2]), abs=1e-4)


def test_a_falling_ball_approaches_the_speed_where_drag_balances_its_weight():
    for name, (m, _) in env.BALLS.items():
        for w, g, rho in AIRS:
            k = W.drag_k(name, rho)
            vt = math.sqrt(g / k)
            assert 0.5 * rho * env.drag_area(name) * vt * vt == pytest.approx(m * g)
            f = W.fly(name, (0.0, 0.0, 5000.0), [[0.0, 0.0, 0.0]], g, rho, t_max=20.0)
            assert np.allclose(-f['v'][0, :, 2], vt * np.tanh(g * f['t'] / vt), rtol=1e-6, atol=1e-6)
            assert -f['v'][0, -1, 2] == pytest.approx(vt, rel=2e-3 if g == G_M else 1e-6)


@pytest.mark.parametrize('v0', env.RELEASE_M_S)
def test_best_ranges_match_the_envelopes_ball_flights(v0):
    ang = np.radians(np.arange(5.0, 85.5, 1.0))
    for w, g, rho in AIRS:
        ref = env.ball_flights(BALL, v0, g, rho)
        f = W.fly(BALL, (0.0, 0.0, env.RELEASE_HEIGHT_M), np.c_[v0 * np.cos(ang), 0 * ang, v0 * np.sin(ang)], g, rho)
        i = int(np.argmax(f['land_p'][:, 0]))
        assert f['land_p'][i, 0] == pytest.approx(ref['range_m'], rel=0.01, abs=0.1)
        assert abs(math.degrees(ang[i]) - ref['best_angle_deg']) <= 2.0


@pytest.mark.parametrize('grip', GRIPS)
def test_release_limits_hold_the_envelopes_moment_balance_at_both_levers(grip):
    m = env.BALLS[BALL][0]
    labels = dict(feet_under='feet under the body', split='long split stance', landing='in a landing, feet under')
    for b in BODIES.values():
        held = env.throws(ENV[b.key], RHO_MOON)['fastest_release_held_m_s']
        for w, g in WORLDS:
            for stance, label in labels.items():
                over = W.release_limit(b, stance, BALL, g, grip)
                N, reach_m = over['floor_push_N'], over['reach_m']
                assert over['turning_m_s'] == pytest.approx(held[label][w], abs=0.051)
                assert over['sliding_m_s'] == pytest.approx(grip * N * env.THROW_S / m)
                assert over['m_s'] == min(over['turning_m_s'], over['sliding_m_s'], RELEASE_MAX_M_S)
                # the hand force at the limit is what the balance names: F lever = N reach where turning binds
                F = over['hand_force_N']
                assert F == pytest.approx(m * over['m_s'] / env.THROW_S)
                assert F * b.release_above_com_m <= N * reach_m + 1e-9 and F <= grip * N + 1e-9
                # a chest pass's short lever: the turning limit grows by the levers' ratio, the grip's stays
                chest = W.release_limit(b, stance, BALL, g, grip, lever=CHEST)
                assert chest['turning_m_s'] == pytest.approx(over['turning_m_s'] * b.release_above_com_m / CHEST)
                assert chest['sliding_m_s'] == pytest.approx(over['sliding_m_s'])
                assert chest['m_s'] == min(chest['turning_m_s'], chest['sliding_m_s'], RELEASE_MAX_M_S)
    # the dance body with its feet under it on the Moon: about 3.0 m/s overhand, 9.9 m/s from the chest
    d = BODIES['dancer']
    assert W.release_limit(d, 'feet_under', BALL, G_M, 0.8)['m_s'] == pytest.approx(2.98, abs=0.01)
    chest = W.release_limit(d, 'feet_under', BALL, G_M, 0.8, lever=CHEST)
    assert chest['m_s'] == pytest.approx(9.92, abs=0.01) and chest['limited_by'] == 'turning'
    assert W.release_limit(d, 'feet_under', BALL, G_E, 0.8, lever=CHEST)['limited_by'] == 'arm'
    # the split stance's reach is half its lever, so on the dance floor's 0.5 the turning and sliding limits meet
    split = W.release_limit(d, 'split', BALL, G_M, 0.5)
    assert split['turning_m_s'] == pytest.approx(split['sliding_m_s'])
    assert W.release_limit(d, 'split', BALL, G_M, 0.5, lever=CHEST)['limited_by'] == 'sliding'


def test_a_catch_mirrors_the_throw_and_a_longer_give_softens_it():
    for b in BODIES.values():
        for w, g in WORLDS:
            for stance in ('feet_under', 'split', 'landing'):
                for lever in (None, CHEST):
                    r = W.release_limit(b, stance, BALL, g, 0.8, lever=lever)
                    c = W.catch_limit(b, stance, BALL, g, 0.8, lever=lever)
                    assert c['m_s'] == pytest.approx(r['m_s']) and c['limited_by'] == r['limited_by']
                    soft = W.catch_limit(b, stance, BALL, g, 0.8, give_s=2 * CATCH_S, lever=lever)
                    assert soft['turning_m_s'] == pytest.approx(2 * c['turning_m_s'])
                    assert soft['sliding_m_s'] == pytest.approx(2 * c['sliding_m_s'])


def test_a_throw_in_the_air_is_held_by_the_arm_alone_and_turns_the_body_by_the_landing():
    m = env.BALLS[BALL][0]
    for b in BODIES.values():
        for w, g in WORLDS:
            for lever in (b.release_above_com_m, CHEST):
                a = W.release_limit(b, 'aloft', BALL, g, 0.8, lever=lever)
                L = m * RELEASE_MAX_M_S * lever
                assert a['m_s'] == RELEASE_MAX_M_S and a['angular_momentum_N_m_s'] == pytest.approx(L)
                for shape in ('layout', 'tuck'):
                    assert a['turn_by_landing_deg'][shape] == pytest.approx(
                        env.attitude_error_deg(L, b.inertia_kg_m2[shape], g))
                one = W.release_limit(b, 'aloft', BALL, g, 0.8, lever=lever, descent_s=1.0)
                turn = math.degrees(L / b.inertia_kg_m2['layout'])
                assert one['turn_by_landing_deg']['layout'] == pytest.approx(turn)
    # the envelope's table: 10 m/s overhead at 0.8 m turns the dance body this far laid out; the arm's 20 m/s doubles it
    d = BODIES['dancer']
    row = next(r for r in env.throws(ENV['dancer'], RHO_MOON)['in_the_air']
               if r['release_m_s'] == 10.0 and r['throw'].startswith('overhead'))
    a = W.release_limit(d, 'aloft', BALL, G_M, 0.8, lever=0.8)
    assert a['turn_by_landing_deg']['layout'] == pytest.approx(2 * row['moon']['layout_deg'], abs=1.0)


def test_a_pass_fan_follows_the_vacuum_and_reports_a_moving_throwers_drift():
    b = BODIES['dancer']
    for w, g in WORLDS:
        for lever in (None, CHEST):
            thrower = State(p=(1.0, 2.0), posture='split')
            fan = W.passes(b, thrower, 'split', BALL, (16.0, 2.0), g, 0.0, 0.8, lever=lever)
            z0 = b.com_m('split') + (b.release_above_com_m if lever is None else lever)
            assert fan['release_m'][2] == pytest.approx(z0) and fan['distance_m'] == pytest.approx(15.0)
            assert fan['limit_m_s'] == pytest.approx(W.release_limit(b, 'split', BALL, g, 0.8, lever=lever)['m_s'])
            assert np.all(fan['held']) and fan['speeds_m_s'].max() == pytest.approx(fan['limit_m_s'])
            S, E = np.meshgrid(fan['speeds_m_s'], np.radians(fan['elevations_deg']), indexing='ij')
            T = 15.0 / (S * np.cos(E))
            z = z0 + S * np.sin(E) * T - 0.5 * g * T * T
            ok = np.isfinite(fan['pass_s'])
            assert ok.any() and np.array_equal(ok, z > 0)
            assert np.allclose(fan['pass_s'][ok], T[ok], atol=1e-6)
            assert np.allclose(fan['height_m'][ok], z[ok], atol=1e-4)
            assert np.allclose(fan['drift_m'][ok], 0.0, atol=1e-9)
            assert fan['path'].shape[:2] == S.shape and fan['path'].shape[2] == fan['t'].size
        mover = State(p=(1.0, 2.0), v=(0.0, 1.0), posture='moving')
        fan = W.passes(b, mover, 'split', BALL, (16.0, 2.0), g, 0.0, 0.8, speeds=[15.0], elevations=[20.0])
        assert fan['drift_m'][0, 0] == pytest.approx(1.0 * fan['pass_s'][0, 0], abs=1e-6)
    # in the Moon's air a pass at the split stance's limit reaches 15 m at catchable heights within seconds
    fan = W.passes(b, State(posture='split'), 'split', BALL, (15.0, 0.0), G_M, RHO_MOON, 0.8)
    lo, hi = W.catch_band(b, 'ready', G_M)['standing_m']
    fast = fan['pass_s'][-1][(fan['height_m'][-1] >= lo) & (fan['height_m'][-1] <= hi)]
    assert fast.size and 1.0 < fast.min() < 2.0


def test_a_lob_reaches_its_point_on_time_through_the_air():
    for w, g, rho in AIRS:
        for T in (1.5, 3.0):
            out = W.lob(BALL, (0.0, 0.0, 2.0), (15.0, 3.0, 1.5), T, g, rho)
            assert out['miss_m'] < 1e-4 and out['lowest_m'] > 0
            f = W.fly(BALL, (0.0, 0.0, 2.0), [out['v0']], g, rho, dt=0.001)
            assert np.allclose(f['p'][0, int(round(T / 0.001))], (15.0, 3.0, 1.5), atol=1e-3)
            assert out['arrives_m_s'] == pytest.approx(np.linalg.norm(f['v'][0, int(round(T / 0.001))]), abs=1e-3)
        vac = W.lob(BALL, (0.0, 0.0, 2.0), (15.0, 3.0, 1.5), 2.0, g, 0.0)
        assert vac['v0'] == pytest.approx(np.array([7.5, 1.5, -0.25 + g]), abs=1e-6)


def test_the_catch_band_spans_the_standing_reach_and_the_best_jump():
    for b in BODIES.values():
        for w, g in WORLDS:
            for posture in ('upright', 'ready', 'split', 'moving'):
                c = W.catch_band(b, posture, g)
                assert c['standing_m'] == pytest.approx((b.reach_low_m, b.reach_up_m))
                assert c['in_posture_m'] == pytest.approx(b.hand_band(b.com_m(posture)))
                vz = b.take_off_m_s(g)
                assert c['jump_top_m'] == pytest.approx(reach.take_off_height(b) + vz * vz / (2 * g) + b.reach_up_m
                                                        - b.com_height_m)
        assert W.catch_band(b, 'ready', G_M)['jump_top_m'] > 1.9 * W.catch_band(b, 'ready', G_E)['jump_top_m']
