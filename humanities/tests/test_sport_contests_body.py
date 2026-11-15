"""The contests between bodies keep their physics: equal chasers never close a gap on an open floor, every step of a
chase keeps the floor law and every flight its arc, a cut turns within its bounce's grip, the shoves' forces are the
vocabulary's and balance each body's momentum, and the time-stepped measures converge as the step halves."""
import math

import pytest

from humanities.sport import contest_chase as chase
from humanities.sport import contest_shove as shove
from humanities.sport import envelope as env
from humanities.sport import vocabulary
from humanities.sport.state import BODIES, G_E, G_M, GRIPS, V_MAX_M_S, WORLDS, muscle_push

UNITS = {u['id']: u for u in vocabulary.load()['units']}


@pytest.mark.parametrize('grip', GRIPS)
def test_equal_players_never_close_a_gap_on_an_open_floor(grip):
    for w, g in WORLDS:
        for plan in ('evade', 'circuit'):
            for pursuer in ('dancer', 'larger'):             # the same push per kilogram, whatever the mass
                out = chase.contest(g, grip, None, 'run', 'extrapolating', pursuer, 'dancer', t_max=20.0, plan=plan)
                assert out['tag_s'] is None and out['closest_m'] == pytest.approx(chase.GAP_M, abs=1e-9)


def test_every_step_of_a_chase_keeps_the_floor_law_and_every_flight_its_arc():
    for w, g in WORLDS:
        for tactic in ('hop', 'cut'):
            out = chase.contest(g, 0.8, 10.0, tactic, 'extrapolating', trace=True, t_max=20.0)
            tr = out['trace']
            dt = chase.DT
            for a, b in zip(tr, tr[1:]):
                # the pursuer: never more than its push, never faster than V_MAX_M_S (a wall may only take speed away)
                dv = math.hypot(b[3] - a[3], b[4] - a[4]) / dt
                assert dv <= b[5] + 1e-6 or math.hypot(b[3], b[4]) < math.hypot(a[3], a[4]) + 1e-12
                assert math.hypot(b[3], b[4]) <= V_MAX_M_S + 1e-9
                assert b[5] <= min(0.8 * g, muscle_push(0.0) * 1.5) + 1e-9
                # the evader aloft: the arc is fixed, so its horizontal velocity holds from step to step
                if a[11] == 'air' and b[11] == 'air':
                    assert (b[8], b[9]) == pytest.approx((a[8], a[9]), abs=1e-12)
            assert out['tactics_used'] >= 1


def test_a_cut_turns_within_the_grip_on_its_bounces_vertical_impulse():
    for w, g in WORLDS:
        for speed in (2.0, 3.0, 4.0):
            for angle in (30, 60, 90):
                for grip in GRIPS:
                    wv = chase.cut_bounce_m_s(speed, angle, grip, g)
                    dv = 2 * speed * math.sin(math.radians(angle) / 2)
                    if wv > 0:
                        assert grip * (2 * wv + g * env.STRIDE_S) == pytest.approx(dv)
                    else:
                        assert grip * g * env.STRIDE_S >= dv
                    if grip == env.MU:
                        ref = env.cut_bounce(speed, angle, grip, g)['vertical_speed_each_way_m_s']
                        assert wv == pytest.approx(ref, abs=6e-4)


def test_a_court_chase_converges_as_the_step_halves():
    for plan, court in (('circuit', 10.0), ('evade', 10.0)):
        coarse = chase.contest(G_M, 0.8, court, plan=plan, t_max=30.0)
        fine = chase.contest(G_M, 0.8, court, plan=plan, t_max=30.0, dt=chase.DT / 2)
        assert coarse['tag_s'] is not None and fine['tag_s'] is not None
        assert fine['tag_s'] == pytest.approx(coarse['tag_s'], rel=0.03)


@pytest.mark.parametrize('grip', GRIPS)
def test_the_shoves_steady_push_is_the_vocabularys(grip):
    for b in BODIES.values():
        for w, g in WORLDS:
            figs = vocabulary.times(UNITS['shove'], b, g, grip)['figures']
            brace = vocabulary.times(UNITS['brace'], b, g, grip)['figures']
            split = shove.steady(b, 'split, facing the push', g, grip)
            assert split['push_N'] == pytest.approx(figs['push_N'], rel=1e-9)
            assert split['grip_N'] == pytest.approx(figs['slip_N'])
            assert split['lean_N'] == pytest.approx(figs['topple_N'])
            squared = shove.steady(b, 'wide crouch, squared up', g, grip)
            assert squared['push_N'] == pytest.approx(brace['squared_up_holds_N'])
    d = BODIES['dancer']
    assert shove.steady(d, 'split, facing the push', G_M, grip)['push_N'] == pytest.approx(38.8, abs=0.05)
    assert shove.steady(d, 'split, facing the push', G_E, grip)['push_N'] == pytest.approx(234.5, abs=0.1)
    assert shove.steady(d, 'wide crouch, squared up', G_M, grip)['share_of_weight'] == pytest.approx(0.087, abs=1e-3)


@pytest.mark.parametrize('grip', GRIPS)
def test_moving_an_opponent_balances_each_bodys_momentum(grip):
    for w, g in WORLDS:
        for a_key, b_key in (('dancer', 'dancer'), ('larger', 'dancer'), ('dancer', 'larger')):
            A, B = BODIES[a_key], BODIES[b_key]
            for stance in ('split, facing the push', 'wide crouch, squared up'):
                for deeper in (0.0, 0.1, 0.2):
                    m = shove.move(A, B, g, grip, stance, deeper)
                    if not m['moves']:
                        assert m['pusher_push_N'] <= m['holder_hold_N'] + 1e-9
                        continue
                    a = m['acceleration_m_s2']
                    assert m['contact_N'] - m['holder_hold_N'] == pytest.approx(B.mass_kg * a)
                    assert m['pusher_floor_N'] - m['contact_N'] == pytest.approx(A.mass_kg * a)
                    lead = A.split_reach_m + deeper
                    turn = m['contact_N'] * shove.hands_m(A, 'split') + A.mass_kg * a * A.com_m('split')
                    assert turn <= A.mass_kg * g * lead + 1e-9 and m['pusher_floor_N'] <= grip * A.mass_kg * g + 1e-9
                    assert (turn == pytest.approx(A.mass_kg * g * lead) if m['limited_by'] == 'lean'
                            else m['pusher_floor_N'] == pytest.approx(grip * A.mass_kg * g))
                    # the pair stepped in time covers the distance in the closed form's time
                    for dt in (1e-3, 5e-4):
                        x = v = t = 0.0
                        while x < shove.DISTANCE_M:
                            v += a * dt
                            x += v * dt
                            t += dt
                        assert t == pytest.approx(m['seconds'], abs=2 * dt)
        # equal players braced in the split stance hold each other; a deeper lean breaks the stalemate
        d = BODIES['dancer']
        assert not shove.move(d, d, g, grip)['moves']
        if grip > d.split_reach_m / shove.hands_m(d, 'split'):
            assert shove.move(d, d, g, grip, deeper_m=0.1)['moves']


def test_a_lean_builds_on_gravitys_clock():
    d = BODIES['dancer']
    for deeper in shove.DEEPER_M:
        for held in (0.0, 0.5):
            moon, earth = shove.lean_build_s(d, G_M, deeper, held), shove.lean_build_s(d, G_E, deeper, held)
            assert moon / earth == pytest.approx(math.sqrt(G_E / G_M))
            h, x0 = d.com_m('split'), d.split_reach_m
            for dt in (1e-4, 5e-5):                          # x'' = (g/h)(x - x_eq) stepped in time
                x, v, t = x0, 0.0, 0.0
                while x < x0 + deeper:
                    v += G_M / h * (x - held * x0) * dt
                    x += v * dt
                    t += dt
                assert t == pytest.approx(moon, abs=2e-3)


def test_a_yield_and_its_recovery_converge_as_the_step_halves():
    for b in BODIES.values():
        for w, g in WORLDS:
            for deeper in (0.0, 0.1, 0.2):
                coarse = shove.yielded(b, g, 0.8, deeper, dt=1e-3)
                fine = shove.yielded(b, g, 0.8, deeper, dt=5e-4)
                assert fine['step_needed'] == coarse['step_needed']
                assert fine['travel_m'] == pytest.approx(coarse['travel_m'], rel=0.01, abs=2e-3)
                assert fine['recover_s'] == pytest.approx(coarse['recover_s'], rel=0.01, abs=5e-3)
    # the lurch is gentler on the Moon but its recovery slower; a deeper lean costs more on either world
    d = BODIES['dancer']
    moon, earth = shove.yielded(d, G_M, 0.8), shove.yielded(d, G_E, 0.8)
    assert moon['travel_m'] < earth['travel_m'] and moon['recover_s'] > earth['recover_s']
    assert shove.yielded(d, G_M, 0.8, 0.2)['travel_m'] > moon['travel_m']


def test_a_shove_aloft_conserves_momentum_and_turns_by_its_angular_impulse():
    for a_key in ('dancer', 'larger'):
        A, T = BODIES[a_key], BODIES['dancer']
        for w, g in WORLDS:
            out = shove.jumper(A, T, g, 0.8)
            air = out['aloft']
            F, t = out['shove_N'], air['contact_s']
            assert T.mass_kg * air['speed_change_m_s'] == pytest.approx(F * t)
            left = air['descent_s'] - t / 2
            assert air['landing_moved_m'] == pytest.approx(air['speed_change_m_s'] * left)
            for shape in ('layout', 'tuck'):
                assert air['turn_by_landing_deg'][shape] == pytest.approx(
                    math.degrees(F * t * shove.HANDS_ABOVE_COM_M / T.inertia_kg_m2[shape] * left))
            assert out['dip_fall_s'] < env.dip_s(T.dip_m, g)
    # the take-off's push holds a shove on the Moon, where its floor force is many lunar weights; not on Earth
    assert shove.jumper(BODIES['dancer'], BODIES['dancer'], G_M, 0.8)['push_absorbs']
    assert not shove.jumper(BODIES['dancer'], BODIES['dancer'], G_E, 0.8)['push_absorbs']
