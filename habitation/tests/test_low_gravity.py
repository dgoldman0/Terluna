"""The low-gravity design note's figures follow from the shared gravity, regenerate exactly, and keep their design
choices labelled."""
import hashlib
import json
import math
from pathlib import Path

import pytest

from habitation.low_gravity import movement as lg
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(lg.OUT.read_text())
DAILY_LIFE = ROOT / 'research/studies/daily_life/results/daily_life.json'


def test_product_regenerates_and_binds_its_code_and_constants():
    assert PRODUCT == json.loads(json.dumps(lg.results()))
    assert PRODUCT['schema'] == lg.SCHEMA
    for path, h in PRODUCT['producer']['files'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16] == h, path
    assert set(PRODUCT['producer']['constants']) == {'MOON_SURFACE_GRAVITY', 'STANDARD_GRAVITY'}
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_gravity_is_the_shared_one():
    assert math.isclose(lg.RATIO, STANDARD_GRAVITY / MOON_SURFACE_GRAVITY)
    assert PRODUCT['gravity']['ratio'] == pytest.approx(6.038, abs=1e-3)


def test_the_same_push_rises_six_times_as_high_and_lasts_six_times_as_long():
    for hop in PRODUCT['hops']:
        v = lg.take_off(hop['earth_rise_m'])
        assert math.isclose(lg.rise(v, STANDARD_GRAVITY), hop['earth_rise_m'])
        assert math.isclose(lg.rise(v) / hop['earth_rise_m'], lg.RATIO)
        assert math.isclose(lg.airtime(v) / lg.airtime(v, STANDARD_GRAVITY), lg.RATIO)
        assert hop['moon_rise_m'] == pytest.approx(lg.RATIO * hop['earth_rise_m'], abs=0.006)


def test_falls_land_like_a_sixth_of_the_drop_on_earth():
    for f in PRODUCT['falls']:
        assert f['like_earth_m'] == pytest.approx(f['drop_m'] / lg.RATIO, abs=6e-4)
        assert f['impact_m_s'] == pytest.approx(math.sqrt(2 * STANDARD_GRAVITY * f['like_earth_m']), abs=0.01)
    for t in PRODUCT['unguarded_drops'].values():
        assert t['moon_drop_same_impact_m'] == pytest.approx(t['earth_drop_m'] * lg.RATIO, abs=0.006)
        # At the same impact speed v the fall lasts v / g, so it lasts RATIO times as long on the Moon.
        assert t['moon_seconds'] / t['earth_seconds'] == pytest.approx(lg.RATIO, rel=0.03)


def test_short_falls_agree_with_the_daily_life_study():
    if not DAILY_LIFE.is_file():
        pytest.skip('daily-life product absent')
    bodies = json.loads(DAILY_LIFE.read_text())['bodies']
    assert bodies['jump_height_factor'] == pytest.approx(lg.RATIO, abs=0.005)
    ours = {f['drop_m']: f['impact_m_s'] for f in PRODUCT['falls']}
    for theirs in bodies['falls']:
        if theirs['height_m'] in ours and theirs['height_m'] <= 3.0:     # drag lowers theirs by under 1% here
            assert ours[theirs['height_m']] == pytest.approx(theirs['impact_m_s'], rel=0.01)
    assert bodies['traction']['friction'] == lg.FRICTION
    assert lg.stop_m(5.0) == pytest.approx(bodies['traction']['stop_from_5_m_s_m']['moon'], abs=0.05)


def test_the_guard_rule_starts_at_earths_guard_and_gives_the_summit_drawings_2_2_m():
    assert lg.guard_height(0.0) == lg.EARTH_GUARD_M
    assert lg.guard_height(lg.GUARD_HOP_EARTH_M) == pytest.approx(2.2)
    g = PRODUCT['guards']
    assert g['summit_drawings']['covers_earth_hop_m'] == pytest.approx(0.218, abs=1e-3)
    assert g['by_hop'][0]['uk_based_m'] == 1.1 and g['by_hop'][0]['ibc_based_m'] == 1.07
    heights = [r['uk_based_m'] for r in g['by_hop']]
    assert heights == sorted(heights)
    # A fit adult's full jump clears the 2.2 m guard, so its top has to give no grip.
    cmj = next(h for h in PRODUCT['hops'] if h['name'] == 'countermovement jump')
    assert cmj['guard_above_its_top_m'] > 2.2


def test_guard_classes_follow_the_earth_codes_carried_to_the_same_impact():
    classes = PRODUCT['guards']['classes']
    assert [c['kind'] for c in classes] == ['design choice'] * 3
    assert classes[0]['drops_m'][1] == pytest.approx(0.38 * lg.RATIO, abs=0.006)
    assert classes[1]['drops_m'] == [classes[0]['drops_m'][1], pytest.approx(0.762 * lg.RATIO, abs=0.006)]
    assert classes[2]['drops_m'][0] == classes[1]['drops_m'][1] and classes[2]['drops_m'][1] is None
    assert [c['guard_m'] for c in classes] == [None, lg.EARTH_GUARD_M, 2.2]


def test_pushes_fall_with_gravity_and_impacts_do_not():
    p = PRODUCT['guards']['push_and_impact']
    assert p['standing_push_n']['earth'] / p['standing_push_n']['moon'] == pytest.approx(lg.RATIO, rel=0.01)
    assert p['impact_energy_j']['2 m/s'] == round(0.5 * lg.BODY_KG * 4)


def test_stairs():
    s = PRODUCT['stairs']
    for kind in ('everyday', 'bounding'):
        st = s[kind]
        assert st['kind'] == 'design choice'
        assert st['work_per_riser_j_kg'] == pytest.approx(MOON_SURFACE_GRAVITY * st['rise_m'], abs=1e-3)
        assert st['step_down']['landing_m_s'] == pytest.approx(math.sqrt(2 * MOON_SURFACE_GRAVITY * st['rise_m']),
                                                               abs=0.006)
        # A step down lands as softly as an Earth step a sixth as high.
        assert st['step_down']['like_earth_step_m'] == pytest.approx(st['rise_m'] / lg.RATIO, abs=1e-4)
        for b in st['bounds']:
            # The bound's top sits at the landing tread: its rise is the risers' height.
            assert lg.rise(b['take_off_up_m_s']) == pytest.approx(b['rise_m'], rel=0.02)
            assert b['like_earth_hop_m'] == pytest.approx(b['rise_m'] / lg.RATIO, abs=1e-3)
    everyday = s['everyday']
    lo, hi = lg.EARTH_STAIR['general_access_rise_m']
    assert lo <= everyday['rise_m'] <= hi
    assert everyday['two_r_plus_g_m'] <= lg.EARTH_STAIR['two_r_plus_g_m'][1] + 1e-9
    assert s['riser_landing_like_earth_170mm_m'] == pytest.approx(0.17 * lg.RATIO, abs=0.006)
    assert s['flight_limits']['UK general access, 12 risers of 170 mm']['moon_rise_same_tumble_m'] == pytest.approx(
        12 * 0.17 * lg.RATIO, abs=0.06)
    assert s['flight_rise_choice_m'] <= s['flight_limits']['UK general access, 12 risers of 170 mm'][
        'moon_rise_same_tumble_m']
    assert s['landing_need_from_2_m_s_m'] <= s['landing_m']
    assert s['headroom_m']['choice_m'] == pytest.approx(s['headroom_m']['UK escape routes'], abs=0.01)


def test_stops_and_turns_take_six_times_the_room():
    for row in PRODUCT['stops'].values():
        assert row['moon_m'] / row['earth_m'] == pytest.approx(lg.RATIO, rel=0.07)          # rounding at 1 m/s
        assert row['moon_turn_radius_m'] / row['earth_turn_radius_m'] == pytest.approx(lg.RATIO, rel=0.04)


def test_handholds_hold_six_times_their_grip_in_body_mass():
    h = PRODUCT['handholds']
    for g in h['grips'].values():
        assert g['lunar_body_kg_one_hand'] == pytest.approx(g['kgf'] * lg.RATIO, abs=0.6)
        # On Earth no grip in the table stops a fall of more than a few centimetres with both hands.
        assert g['arrest_m']['earth_two_hands'] < 0.15
        assert g['arrest_m']['moon_two_hands'] > g['arrest_m']['moon_one_hand']
    assert h['grips']["women's weak grip (T-score -2.5)"]['lunar_body_kg_one_hand'] > lg.BODY_KG


def test_ramp_run_out_does_not_depend_on_gravity():
    for s in (1 / 12, 1 / 15, 1 / 20):
        theta = math.atan(s)
        for g in (MOON_SURFACE_GRAVITY, STANDARD_GRAVITY):
            run = 10.0
            v2 = 2 * g * math.sin(theta) * run
            assert v2 / (2 * lg.FRICTION * g) / run == pytest.approx(lg.runaway_run_out(s))
    r = PRODUCT['ramps']
    for row in r['by_slope'].values():
        # A 3 m landing stops a runaway from any run that rises up to friction * landing, whatever its slope.
        assert row['rise_of_that_run_m'] == pytest.approx(r['run_rise_for_the_landing_m'], abs=0.01)
        assert row['stop_on_the_ramp_from_walking_m']['moon'] / row['stop_on_the_ramp_from_walking_m']['earth'] == \
            pytest.approx(lg.RATIO, rel=0.05)
    assert r['run_rise_for_the_landing_m'] > lg.EARTH_RAMP['ada_max_rise_m']
    assert r['lunar_slope_with_the_climb_of_earth_1_12_deg'] == pytest.approx(30.1, abs=0.1)
    # The slope that climbs like Earth's 1:12 leaves almost no grip for stopping on the way down.
    steep = r['steep_paths']
    assert steep['30 degrees']['share_of_grip'] > 0.95
    shares = [row['share_of_grip'] for row in steep.values()]
    stops = [row['stop_from_walking_m']['moon'] for row in steep.values()]
    assert shares == sorted(shares) and stops == sorted(stops)


def test_banking_and_halls():
    assert lg.bank_deg(1.0, 10.0) == 0.0                   # grip alone holds a gentle turn
    beta = math.radians(lg.bank_deg(3.0, 3.0))
    a = MOON_SURFACE_GRAVITY * (math.sin(beta) + lg.FRICTION * math.cos(beta)) / (
        math.cos(beta) - lg.FRICTION * math.sin(beta))
    assert a == pytest.approx(9.0 / 3.0, rel=1e-6)
    m = PRODUCT['movement_halls']
    assert m['kind'] == 'design choice'
    assert m['clear_height_m'] >= m['needed_m'] > m['head_top_m']
    assert m['running_airtime_ratio'] == pytest.approx(math.sqrt(lg.RATIO), abs=0.005)
    assert m['running_hop']['length_m'] / m['running_hop']['earth_length_m'] == pytest.approx(lg.RATIO, rel=0.04)
    st = PRODUCT['stepping_stones']
    assert st['moon_length_m'] / st['earth_length_m'] == pytest.approx(lg.RATIO, rel=0.04)
