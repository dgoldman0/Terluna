"""The sport study's products regenerate and bind their inputs, and each figure follows from its stated mechanics,
checked by an independent route where one exists."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from humanities.sport import envelope as env
from humanities.sport import whole_body as wb
from shared.constants import MOON_SURFACE_GRAVITY as GM, STANDARD_GRAVITY as GE
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(env.OUT.read_text())
DANCER = 'dancer, 1.64 m and 55 kg'
LARGER = 'larger athlete, 1.80 m and 80 kg (scaled)'


def test_envelope_regenerates_and_binds_its_code_constants_and_inputs():
    assert PRODUCT == json.loads(json.dumps(env.results()))
    assert PRODUCT['schema'] == env.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16] == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_the_larger_body_scales_by_mass_and_length_squared():
    small, large = PRODUCT['bodies'][DANCER], PRODUCT['bodies'][LARGER]
    sl, sm = large['stature_m'] / small['stature_m'], large['mass_kg'] / env.MASS_KG
    for k, v in small['inertia_kg_m2'].items():
        assert large['inertia_kg_m2'][k] == pytest.approx(v * sm * sl ** 2, rel=2e-3)
    assert large['com_height_m'] == pytest.approx(small['com_height_m'] * sl, rel=2e-3)
    assert 0.53 < small['com_height_m'] / small['stature_m'] < 0.57
    # the tuck turns faster than the layout, as Yeadon's shares have it
    assert 0.2 < small['inertia_kg_m2']['tuck'] / small['inertia_kg_m2']['layout'] < 0.35


def test_a_lean_follows_the_inverted_pendulum_and_its_clock():
    for g in (GM, GE):
        out = env.lean_start(0.1, g, 0.9, 0.13)
        tau = math.sqrt(0.9 / g)
        assert 0.13 * (math.cosh(out['seconds'] / tau) - 1) == pytest.approx(0.1, rel=0.03)
    moon, earth = env.lean_start(0.3, GM, 0.9, 0.13), env.lean_start(0.3, GE, 0.9, 0.13)
    assert moon['seconds'] / earth['seconds'] == pytest.approx(math.sqrt(GE / GM), rel=0.01)
    # integrated directly: x'' = (g/h)(x + reach) from rest
    x, v, t, dt = 0.0, 0.0, 0.0, 1e-5
    while x < 0.3:
        v += GM / 0.9 * (x + 0.13) * dt
        x += v * dt
        t += dt
    assert t == pytest.approx(moon['seconds'], abs=0.01)


def test_a_cut_needs_the_bounce_whose_vertical_impulse_grip_turns_into_the_change():
    for row in PRODUCT['by_body'][DANCER]['ground']['cuts_through_one_stance']:
        for world, g in (('moon', GM), ('earth', GE)):
            c = row[world]
            dv = 2 * row['speed_m_s'] * math.sin(math.radians(row['angle_deg']) / 2)
            assert c['velocity_change_m_s'] == pytest.approx(dv, abs=1e-3)
            w = c['vertical_speed_each_way_m_s']
            if w > 0:
                assert env.MU * (2 * w + g * env.STRIDE_S) == pytest.approx(dv, abs=2e-3)
                assert c['bounce_rise_m'] == pytest.approx(w * w / (2 * g), abs=2e-3)
    cut = {(c['speed_m_s'], c['angle_deg']): c for c in PRODUCT['by_body'][DANCER]['ground']['cuts_through_one_stance']}
    assert cut[(3.0, 90)]['moon']['bounce_rise_m'] > 1.0 > 0.1 > cut[(3.0, 90)]['earth']['bounce_rise_m']


def test_a_braking_landing_sheds_its_grip_on_the_vertical_impulse_and_slides_off_the_rest():
    out = env.braking_landing(4.0, 2.5, 0.8, GM, stroke=0.4)
    T = 2 * 0.4 / 2.5
    assert out['stance_s'] == pytest.approx(T, abs=1e-3)
    assert out['shed_in_landing_m_s'] == pytest.approx(0.8 * (2.5 + GM * T), abs=1e-3)
    left = 4.0 - out['shed_in_landing_m_s']
    assert out['slide_after_m'] == pytest.approx(left ** 2 / (2 * 0.8 * GM), abs=2e-3)
    # a landing cannot reverse the travel: what it sheds is capped at the speed itself
    assert env.braking_landing(1.0, 2.5, 0.8, GE)['left_m_s'] == 0.0


def test_getting_down_and_standing_up_obey_gravity_at_their_slow_end():
    for g in (GM, GE):
        free_fall = math.sqrt(2 * 0.3 / g)
        assert env.dip_s(0.3, g) > free_fall
        assert env.dip_s(0.3, g, brake=1e9) == pytest.approx(free_fall, rel=1e-3)
        assert env.rise_and_stop_s(0.3, 1e9, g) == pytest.approx(free_fall, rel=1e-3)
        # integrated directly: fall at g until the legs' braking just stops the body at the depth, then brake
        z, v, t, dt = 0.0, 0.0, 0.0, 1e-6
        while v * v / (2 * env.BRAKE_M_S2) < 0.3 - z:
            v += g * dt
            z += v * dt
            t += dt
        t += v / env.BRAKE_M_S2
        assert t == pytest.approx(env.dip_s(0.3, g), rel=2e-3)
    assert env.dip_s(0.3, GM) / env.dip_s(0.3, GE) > 1.9
    cm = PRODUCT['by_body'][DANCER]['jumps']['countermovement']
    assert cm['dip_s']['moon'] > 2 * cm['dip_s']['earth'] - 0.01


def test_ball_flights_reduce_to_the_vacuum_and_terminal_speeds_balance_weight():
    for v0 in (5.0, 10.0):
        out = env.ball_flights('handball (size 3)', v0, GM, 0.0, dt=1e-3)
        th = np.radians(np.arange(5.0, 85.5, 0.25))
        h0 = env.RELEASE_HEIGHT_M
        R = v0 * np.cos(th) * (v0 * np.sin(th) + np.sqrt((v0 * np.sin(th)) ** 2 + 2 * GM * h0)) / GM
        assert out['range_m'] == pytest.approx(R.max(), rel=0.01)
    rho = PRODUCT['air_density_kg_m3']['moon_sea_level']
    assert 1.38 < rho < 1.42
    for name, b in PRODUCT['balls'].items():
        A = b['drag_area_m2']
        assert 0.5 * rho * A * b['terminal_m_s']['moon'] ** 2 == pytest.approx(b['mass_kg'] * GM, rel=0.03), name
        for t in b['throws']:
            assert t['moon']['range_m'] < t['moon_without_air']['range_m']
            assert t['moon']['arrives_m_s'] < t['release_m_s'] + 0.5
    assert PRODUCT['balls']['shuttlecock']['terminal_m_s']['earth'] == pytest.approx(env.SHUTTLE_EARTH_TERMINAL_M_S, rel=0.02)


def test_a_throw_is_one_moment_balance_between_the_hand_and_the_floor():
    body = PRODUCT['bodies'][DANCER]
    t = PRODUCT['by_body'][DANCER]['throws']
    m_ball = t['ball_kg']
    for row in t['standing']:
        F = m_ball * row['release_m_s'] / env.THROW_S
        for world, g in (('moon', GM), ('earth', GE)):
            N = body['mass_kg'] * g
            assert row[world]['reach_needed_recoiling_m'] == pytest.approx(F * body['release_above_com_m'] / N, abs=6e-3)
            assert row[world]['held_feet_under'] == (F * body['release_above_com_m'] <= N * body['fore_aft_reach_m'])
    held = t['fastest_release_held_m_s']
    v = held['feet under the body']['moon']
    F = m_ball * v / env.THROW_S
    assert F * body['release_above_com_m'] == pytest.approx(body['mass_kg'] * GM * body['fore_aft_reach_m'], rel=0.02)
    assert held['feet under the body']['earth'] / v == pytest.approx(GE / GM, rel=0.02)
    assert held['in a landing, feet under']['moon'] > held['long split stance']['moon'] > 4 * v


def test_a_bump_or_a_throw_in_the_air_turns_the_body_by_its_momentum_over_its_inertia():
    body = PRODUCT['bodies'][DANCER]
    rows = PRODUCT['by_body'][DANCER]['landings']['attitude_after_a_bump_at_the_top']['rows']
    for row in rows:
        for world, g in (('moon', GM), ('earth', GE)):
            expect = math.degrees(row['angular_momentum_N_m_s'] / body['inertia_kg_m2']['layout'] * env.APEX_TAKE_OFF_M_S / g)
            assert row[world]['layout_deg'] == pytest.approx(expect, abs=0.6)
        assert row['moon']['layout_deg'] / max(row['earth']['layout_deg'], 1) > 4
    for row in PRODUCT['by_body'][DANCER]['throws']['in_the_air']:
        arm = env.AIR_LEVERS_M[row['throw']]
        assert row['angular_momentum_N_m_s'] == pytest.approx(env.BALLS[env.THROW_BALL][0] * row['release_m_s'] * arm, abs=6e-3)


def test_collisions_and_tosses_keep_momentum_and_the_force_balance():
    for row in PRODUCT['contact']['ground_collisions']:
        ma, mb = (PRODUCT['bodies'][n]['mass_kg'] for n in row['bodies'])
        dva, dvb = row['speed_change_m_s']
        assert ma * dva == pytest.approx(mb * dvb, rel=0.02)
        assert dva + dvb == pytest.approx(row['closing_m_s'], abs=0.02)
        assert row['moon']['slide_m'][0] == pytest.approx(dva ** 2 / (2 * env.MU * GM), abs=0.01)
    toss = PRODUCT['contact']['toss']['flyers'][DANCER]
    a = env.TOSS_FORCE_N / env.MASS_KG - GM
    assert toss['moon']['rise_m'] == pytest.approx(a * env.TOSS_STROKE_M / GM, rel=0.01)
    assert toss['earth'] is None                                  # 350 N does not lift 55 kg on Earth


@pytest.mark.skipif(not wb.MOTION.is_file(), reason='the solved phrase lives on the research drive')
def test_the_whole_body_record_regenerates_and_its_flights_fall_freely():
    product = json.loads(wb.OUT.read_text())
    assert product == json.loads(json.dumps(wb.results()))
    for path, h in product['producer']['inputs'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == h
    for path, h in product['producer']['files'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16] == h
    shares = [f['share_of_gravity'] for f in product['flights']]
    assert all(0.88 < s < 1.06 for s in shares) and abs(np.median(shares) - 1) < 0.02
    for s in product['stances_between_flights']:
        assert s['most_friction_used'] <= wb.FLOOR_FRICTION + 0.01      # the engine's own friction cone
        assert s['travel_change_m_s'] <= s['friction_bound_on_travel_change_m_s'] + 0.05
    # the hard landings: the touchdown impulse carries most of the vertical speed in one instant
    carve = next(s for s in product['stances_between_flights'] if s['unit'].startswith('carve_landing'))
    assert carve['landing_impulse_m_s']['vertical'] > 0.7 * -carve['before']['vertical_m_s']
    assert product['salto_landing']['landing_impulse_m_s']['vertical'] > 0.9 * -product['salto_landing']['before']['vertical_m_s']
    take_off = product['salto_standing_take_off']
    assert take_off['descent_s'] > 3 * take_off['dip_to_lift_off_s']
