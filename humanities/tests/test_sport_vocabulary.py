"""The athletic vocabulary's product regenerates and binds its inputs, its units are the contract's and join into one
vocabulary, and each unit's times follow the envelope's mechanics, checked by an independent route where one exists."""
import hashlib
import json
import math
from pathlib import Path

import pytest

from humanities.sport import envelope as env
from humanities.sport import state
from humanities.sport import vocabulary as V
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(V.OUT.read_text())
UNITS = {u['id']: u for u in V.load()['units']}
ENVELOPE = json.loads(env.OUT.read_text())
DANCER = state.BODIES['dancer']


def times(uid, body=DANCER, g=state.G_M, mu=state.MU, **given):
    return V.times(UNITS[uid], body, g, mu, given)


def test_the_vocabulary_regenerates_and_binds_its_code_constants_and_inputs():
    assert PRODUCT == json.loads(json.dumps(V.results()))
    assert PRODUCT['schema'] == V.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16] == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_the_units_are_the_contracts_and_every_one_evaluates_on_both_worlds_at_every_grip():
    assert list(UNITS) == list(state.UNITS)
    for u in UNITS.values():
        assert u['lunar'].split(':')[0] in ('native', 'adapted', 'native in its use', 'native in its limit'), u['id']
        for body in state.BODIES.values():
            for _, g in state.WORLDS:
                for mu in state.GRIPS:
                    t = V.times(u, body, g, mu)
                    values = [t['duration_s'], t['tell_s'], t['commits_s'], t['recover_s'], *t['figures'].values()]
                    assert all(math.isfinite(x) and x >= 0 for x in values), (u['id'], body.key, g, mu)
                    overlapping = any(p.get('overlaps') for p in u['phases'])
                    assert t['duration_s'] > 0 or overlapping, u['id']


def test_every_unit_can_be_reached_from_the_ready_crouch_and_lead_back_to_it():
    graph = V.transitions(list(UNITS.values()))
    assert V.reachable(graph, 'ready_crouch') == set(UNITS)
    back = {k: [a for a in graph if k in graph[a]] for k in graph}
    assert V.reachable(back, 'ready_crouch') == set(UNITS)
    assert PRODUCT['transitions']['reach_all_from_ready_crouch']
    # a hand action keeps the posture it found; a landing unit only follows an arc or a stride
    assert 'air_pass' in graph['standing_jump'] and 'air_pass' not in graph['ready_crouch']
    assert 'braking_landing' in graph['side_hop'] and 'braking_landing' not in graph['pivot']


def test_the_highest_jump_shows_itself_as_the_envelopes_countermovement():
    for key, body in state.BODIES.items():
        name = body.name
        for world, g in state.WORLDS:
            t = times('standing_jump', body, g, rise_m=10.0)
            vz = body.take_off_m_s(g)
            assert t['figures']['take_off_m_s'] == pytest.approx(vz)
            assert t['figures']['dip_m'] == pytest.approx(body.dip_m, rel=1e-6)
            assert t['commits_s'] == pytest.approx(2 * vz / g)
            # the full dip and push are the envelope's countermovement (dip at most g, push over push_off at vz/2)
            expect = env.dip_s(body.dip_m, g) + 2 * body.push_off_m / vz
            assert t['tell_s'] == pytest.approx(expect, rel=1e-6)
            if key == 'dancer':
                assert t['tell_s'] == pytest.approx(ENVELOPE['by_body'][name]['jumps']['countermovement']['dip_and_push_s'][world], abs=2e-3)


def test_a_lower_jump_sinks_only_as_deep_as_its_push_needs():
    for _, g in state.WORLDS:
        for vz in (1.0, 1.8, 2.3):
            if DANCER.push_for(vz, g) > DANCER.push_off_m:
                continue
            s = DANCER.push_for(vz, g)
            # the stroke reaches vz exactly on the force-velocity profile, and less stroke would not
            assert DANCER.take_off_m_s(g, s) == pytest.approx(vz, rel=1e-9)
            dip, dip_s, push = DANCER.jump_timing('upright', vz, g)
            lift = DANCER.push_off_m - DANCER.dip_m
            assert dip == pytest.approx(max(state.MIN_DIP_M, s - lift))
            assert dip_s == pytest.approx(env.dip_s(dip, g))
            assert push == pytest.approx((dip + lift + s) / vz)
    # on the Moon a 1 m jump needs only the least dip and leaves within 0.4 s; from the crouch with none, in 0.17 s
    moon = times('standing_jump', rise_m=1.0)
    assert moon['figures']['dip_m'] == pytest.approx(state.MIN_DIP_M)
    assert moon['tell_s'] < 0.4 and times('crouch_jump', rise_m=1.0)['tell_s'] < 0.18
    # beyond the crouch's best take-off the crouch jump caps rather than sinking
    assert times('crouch_jump', rise_m=3.0)['figures']['rise_m'] == pytest.approx(
        DANCER.take_off_m_s(state.G_M, DANCER.crouch_push_m) ** 2 / (2 * state.G_M))


def test_starts_follow_the_lean_then_the_grip_then_the_legs():
    for _, g in state.WORLDS:
        ns = V.physics(DANCER, g, state.MU)
        # inside the lean, the envelope's inverted pendulum
        lean = env.lean_start(0.1, g, DANCER.com_height_m, DANCER.side_reach_m)['seconds']
        assert ns['start_s'](0.1, DANCER.side_reach_m, DANCER.com_height_m) == pytest.approx(lean, abs=2e-3)
    # on the Moon from a lead beyond mu h: the grip alone, x = mu g t^2 / 2
    ns = V.physics(DANCER, state.G_M, state.MU)
    t = ns['start_s'](3.0, 1.0, 0.75)
    assert t == pytest.approx(math.sqrt(2 * 3.0 / (state.MU * state.G_M)), rel=2e-3)
    # on Earth with a high grip and a deep lead: the legs alone, x = V0 (t - V0/F0 (1 - exp(-F0 t / V0)))
    ns = V.physics(DANCER, state.G_E, 1.2)
    F0, V0 = state.SPRINT_F0_N_KG, state.SPRINT_V0_M_S
    t = ns['start_s'](5.0, 1.0, 0.75)
    assert V0 * (t - V0 / F0 * (1 - math.exp(-F0 * t / V0))) == pytest.approx(5.0, rel=3e-3)
    # the drop step's quarter second buys back more than it costs on the Moon
    assert times('drop_step')['figures']['first_metre_s'] < times('lean_start')['figures']['first_metre_s'] - 0.5


def test_a_cut_and_a_braking_landing_are_the_envelopes():
    for _, g in state.WORLDS:
        for mu in state.GRIPS:
            for speed, angle in ((3.0, 45), (2.0, 90), (4.0, 30)):
                t = times('bounce_cut', g=g, mu=mu, speed_m_s=speed, angle_deg=angle)
                assert t['figures']['bounce_rise_m'] == pytest.approx(
                    env.cut_bounce(speed, angle, mu, g)['bounce_rise_m'], abs=6e-4)
            b = times('braking_landing', g=g, mu=mu, speed_m_s=3.0, v_down_m_s=2.0)
            e = env.braking_landing(3.0, 2.0, mu, g)
            assert b['figures']['shed_m_s'] == pytest.approx(e['shed_in_landing_m_s'], abs=6e-4)
            assert b['figures']['stop_m'] == pytest.approx(e['total_m'], abs=2e-3)


def test_throws_and_catches_are_the_envelopes_moment_balance():
    for name, body in ((DANCER.name, DANCER), (state.BODIES['larger'].name, state.BODIES['larger'])):
        held = ENVELOPE['by_body'][name]['throws']['fastest_release_held_m_s']
        for world, g in state.WORLDS:
            cap = state.RELEASE_MAX_M_S
            pairs = (('chest_pass', 'overhand_held_m_s', 'feet under the body'),
                     ('split_throw', 'release_held_m_s', 'long split stance'),
                     ('landing_throw', 'release_held_m_s', 'in a landing, feet under'))
            for uid, key, stance in pairs:
                got = times(uid, body, g)['figures'][key]
                # the envelope's moment balance, then the grip and the arm, which the envelope left out
                assert got == pytest.approx(min(held[stance][world], cap), abs=0.06)
        # a catch held at the chest mirrors the chest pass
        assert times('catch', body)['figures']['chest_holds_m_s'] == pytest.approx(
            times('chest_pass', body)['figures']['release_held_m_s'])


def test_the_pivot_turns_on_the_feet_friction_couple():
    for _, g in state.WORLDS:
        t = times('pivot', g=g, angle_deg=90)['commits_s']
        tau = state.MU * DANCER.mass_kg * g * DANCER.ready_half_stance_m
        alpha = tau / DANCER.inertia_kg_m2['twist']
        # speed the turn for half the time and stop it for the other half
        assert 2 * 0.5 * alpha * (t / 2) ** 2 == pytest.approx(math.pi / 2, rel=1e-9)
    assert times('pivot')['commits_s'] / times('pivot', g=state.G_E)['commits_s'] == pytest.approx(
        math.sqrt(state.G_E / state.G_M), rel=1e-9)


def test_the_boost_lifts_only_on_the_moon_as_the_envelopes_toss():
    toss = ENVELOPE['contact']['toss']['flyers']
    for key, body in state.BODIES.items():
        moon, earth = times('boost', body), times('boost', body, g=state.G_E)
        assert moon['figures']['lifts'] == 1 and earth['figures']['lifts'] == 0
        assert moon['figures']['rise_m'] == pytest.approx(toss[body.name]['moon']['rise_m'], abs=6e-4)
        assert toss[body.name]['earth'] is None


def test_the_vocabularys_starts_agree_with_the_reach_engines_profiles():
    import numpy as np
    from humanities.sport import reach
    for _, g in state.WORLDS:
        for mu in state.GRIPS:
            ns = V.physics(DANCER, g, mu)
            for posture, lead, h, u in (('upright', DANCER.side_reach_m, DANCER.com_height_m, (0.0, 1.0)),
                                        ('ready', DANCER.ready_half_stance_m, DANCER.ready_com_m, (0.0, 1.0))):
                t, x, _ = reach.profile(DANCER, posture, u, mu, g, steps=False)
                for d in (0.3, 1.0, 2.0):
                    # the time the engine's profile first covers d, against the vocabulary's own integration
                    t_engine = float(np.interp(d, x, t))
                    assert ns['start_s'](d, lead, h) == pytest.approx(t_engine, abs=0.02), (posture, g, mu, d)
