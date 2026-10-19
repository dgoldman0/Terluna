"""The low-gravity design note's figures follow from the shared gravity and the cited values, regenerate exactly, and
keep their design choices labelled. Each derivation is checked from its inputs, independently of the module's own
arithmetic where it can be."""
import hashlib
import json
import math
from pathlib import Path

import pytest

from habitation.low_gravity import movement as lg
from shared.constants import MOON_SURFACE_GRAVITY as GM, STANDARD_GRAVITY as GE
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(lg.OUT.read_text())
SOURCES = (ROOT / 'habitation/low_gravity/sources.json').read_text()
DAILY_LIFE = ROOT / 'research/studies/daily_life/results/daily_life.json'
R = GE / GM


def test_product_regenerates_and_binds_its_code_and_constants():
    assert PRODUCT == json.loads(json.dumps(lg.results()))
    assert PRODUCT['schema'] == lg.SCHEMA
    for path, h in PRODUCT['producer']['files'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16] == h, path
    assert set(PRODUCT['producer']['constants']) == {'MOON_SURFACE_GRAVITY', 'STANDARD_GRAVITY'}
    assert not constants_changed(PRODUCT['producer']['constants'])
    assert PRODUCT['gravity']['ratio'] == pytest.approx(6.038, abs=1e-3)


def test_cited_values_match_their_sources():
    cited = [  # (as sources.json gives it, whether the module holds that value)
        ('guarding 900 mm along flights and 1100 mm elsewhere', lg.EARTH_GUARD_M == 1.1),
        ('guards of at least 1,067 mm (42 in)', lg.EARTH_GUARD_IBC_M == 1.067),
        ('two risers or 380 mm',
         lg.UNGUARDED_DROP_M['UK, buildings other than dwellings (two risers, or 380 mm)'] == 0.38),
        ('a fall exceeds 600 mm in dwellings', lg.UNGUARDED_DROP_M['UK, dwellings'] == 0.60),
        ('more than 762 mm (30 in)', lg.UNGUARDED_DROP_M['IBC 1015.2 (30 inches)'] == 0.762),
        ('rise 150-170 mm and go 250-400 mm', lg.EARTH_STAIR['general_access_rise_m'] == (0.150, 0.170)
         and lg.EARTH_STAIR['general_access_going_m'] == (0.250, 0.400)),
        ('2R + G between 550 and 700 mm', lg.EARTH_STAIR['two_r_plus_g_m'] == (0.550, 0.700)),
        ('at most 12 risers a flight', lg.EARTH_FLIGHT_RISE_M['UK general access, 12 risers of 170 mm'] == 12 * 0.170),
        ('at most 3,658 mm (12 ft)', lg.EARTH_FLIGHT_RISE_M['IBC 1011.8 (12 feet)'] == 3.658),
        ('2 m clear headroom', lg.EARTH_HEADROOM_M['UK escape routes'] == 2.0),
        ('2,032 mm (80 in) of headroom', lg.EARTH_HEADROOM_M['IBC (80 inches)'] == 2.032),
        ('handrails 900-1000 mm', lg.EARTH_HANDRAIL_M['UK (Approved Document K)'] == (0.900, 1.000)),
        ('handrails 865-965 mm', lg.EARTH_HANDRAIL_M['ADA 505.4 (34-38 inches)'] == (0.865, 0.965)),
        ('no steeper than 1:12, each rising at most 760 mm', lg.EARTH_RAMP['ada_max_rise_m'] == 0.760),
        ('at least 1,525 mm (60 in) long', lg.EARTH_RAMP['ada_landing_m'] == 1.525),
        ('1:12 for flights to 2 m, 1:15 to 5 m and 1:20 to 10 m',
         lg.EARTH_RAMP['adk_max_going_m'] == {'1:12': 2.0, '1:15': 5.0, '1:20': 10.0}),
        ('35.3 +/- 4.5 cm', lg.MARKOVIC_M['countermovement_after_take_off'] == 0.353),
        ('46.2 +/- 4.8 cm', lg.MARKOVIC_M['abalakow_without_arm_swing'] == 0.462),
        ('52.2 +/- 6.2 cm', lg.MARKOVIC_M['abalakow_with_arm_swing'] == 0.522),
        ('comes at 1.42 m/s', lg.WALK_RUN_MOON_M_S == 1.42),
        ('51 kg in men', lg.GRIP_KGF["men's peak median"] == 51.0),
        ('31 kg in women', lg.GRIP_KGF["women's peak median"] == 31.0),
        ('27 kg in men and 16 kg in women', lg.GRIP_KGF["men's weak grip (T-score -2.5)"] == 27.0
         and lg.GRIP_KGF["women's weak grip (T-score -2.5)"] == 16.0),
        ('ending at 17.9 m/s', lg.TERMINAL_M_S == 17.9),
    ]
    for text, matches in cited:
        assert text in SOURCES, text
        assert matches, text


def test_the_same_take_off_speed_rises_six_times_as_high_and_lasts_six_times_as_long():
    for hop in PRODUCT['hops']:
        v = math.sqrt(2 * GE * hop['earth_rise_m'])
        assert hop['take_off_m_s'] == pytest.approx(v, abs=0.006)
        assert v * v / (2 * GM) == pytest.approx(R * hop['earth_rise_m'])
        assert hop['moon_airtime_s'] / hop['earth_airtime_s'] == pytest.approx(R, rel=0.03)


def test_the_belt_jump_keeps_its_push_off_extension_unscaled():
    # Markovic: the countermovement jump by flight time on a mat, the Abalakow jumps by a waist belt from standing.
    extension = 0.462 - 0.353
    flight = 0.522 - extension
    arm = PRODUCT['jump_bounds']['jump with an arm swing']
    assert arm['rise_measured_from'] == 'standing'
    assert arm['earth_take_off_m_s'] == pytest.approx(math.sqrt(2 * GE * flight), abs=0.006)
    assert arm['same_take_off_speed']['moon_rise_m'] == pytest.approx(extension + R * flight, abs=0.006)
    assert arm['same_take_off_speed']['moon_rise_m'] == pytest.approx(2.60, abs=0.01)


def test_the_same_leg_work_bounds_the_jump_from_above():
    for b in PRODUCT['jump_bounds'].values():
        v_e = b['earth_take_off_m_s']
        for row in b['same_leg_work']:
            d, v_m = row['push_off_m'], row['take_off_m_s']
            # The legs' work per kilogram pays the take-off speed and the lift through the push-off on both worlds.
            assert v_m ** 2 / 2 + GM * d == pytest.approx(v_e ** 2 / 2 + GE * d, rel=0.01)
            assert row['moon_rise_m'] > b['same_take_off_speed']['moon_rise_m']
    cmj = PRODUCT['jump_bounds']['countermovement jump']['same_leg_work']
    assert [r['moon_rise_m'] for r in cmj] == pytest.approx([3.89, 4.40], abs=0.01)


def _impact(h, vt=17.9, steps=20000):
    """A fall with quadratic drag, stepped forward in time: speed at the ground and the time taken."""
    v = y = t = 0.0
    dt = math.sqrt(2 * h / GM) / steps
    while y < h:
        a = GM * (1 - (v / vt) ** 2)
        v, y, t = v + a * dt, y + v * dt + 0.5 * a * dt * dt, t + dt
    return v, t


def test_falls_with_drag():
    for f in PRODUCT['falls']:
        v, t = _impact(f['drop_m'])
        assert f['impact_m_s'] == pytest.approx(v, abs=0.01)
        assert f['seconds'] == pytest.approx(t, abs=0.01)
        assert f['like_earth_m'] == pytest.approx(v ** 2 / (2 * GE), abs=2e-3)
    for h in (0.5, 2.0, 7.0):
        assert lg.drop_for_impact(lg.impact_speed(h)) == pytest.approx(h)
    for t in PRODUCT['unguarded_drops'].values():
        v, _ = _impact(t['moon_drop_same_impact_m'])
        assert v == pytest.approx(math.sqrt(2 * GE * t['earth_drop_m']), abs=0.01)


def test_falls_agree_with_the_daily_life_study():
    if not DAILY_LIFE.is_file():
        pytest.skip('daily-life product absent')
    bodies = json.loads(DAILY_LIFE.read_text())['bodies']
    assert bodies['terminal_m_s'] == lg.TERMINAL_M_S
    assert bodies['jump_height_factor'] == pytest.approx(R, abs=0.005)
    ours = {f['drop_m']: f for f in PRODUCT['falls']}
    shared = [f for f in bodies['falls'] if f['height_m'] in ours]
    assert len(shared) >= 4
    for theirs in shared:
        assert ours[theirs['height_m']]['impact_m_s'] == pytest.approx(theirs['impact_m_s'], abs=0.01)
    assert bodies['traction']['friction'] == lg.FRICTION
    assert lg.stop_m(5.0) == pytest.approx(bodies['traction']['stop_from_5_m_s_m']['moon'], abs=0.05)


def test_the_hop_rule_keeps_earths_margin_at_the_same_hops_top():
    j = (2.2 - 1.1) / (R - 1)
    g = PRODUCT['guards']
    assert g['sized_hop']['earth_hop_m'] == pytest.approx(j, abs=1e-4)
    for row in g['by_hop']:
        h_moon = row['uk_based_m']
        for com in (0.6, 0.96, 1.05):
            earth_margin = 1.1 - (com + row['earth_hop_m'])
            moon_margin = h_moon - (com + R * row['earth_hop_m'])
            assert moon_margin == pytest.approx(earth_margin, abs=0.006)
    over = g['sized_hop']['example_adult']['com_over_top_at_the_hops_top_m']
    assert over['earth'] == pytest.approx(0.55 * 1.75 + j - 1.1, abs=1e-3)
    assert over['moon'] == over['earth'] and over['earth'] == pytest.approx(0.08, abs=0.005)


def test_the_energy_rule_costs_earths_lift():
    g = PRODUCT['guards']
    for row in g['by_stature']:
        c = 0.55 * row['stature_m']
        assert row['com_m'] == pytest.approx(c, abs=1e-3)
        assert GM * (row['moon_guard_same_lift_m'] - c) == pytest.approx(GE * (1.1 - c), abs=0.05)
        assert row['earth_lift_speed_m_s'] == pytest.approx(math.sqrt(2 * GE * (1.1 - c)), abs=0.006)
        assert row['moon_lift_speed_over_2_2_m_m_s'] == pytest.approx(math.sqrt(2 * GM * (2.2 - c)), abs=0.006)
    # 2.2 m costs Earth's lift down to the centre of mass where g_M (2.2 - c) = g_E (1.1 - c).
    c_kept = (GE * 1.1 - GM * 2.2) / (GE - GM)
    kept = g['energy_kept']
    assert kept['lowest_com_m'] == pytest.approx(c_kept, abs=1e-3)
    for row in g['by_stature']:
        keeps = row['moon_lift_speed_over_2_2_m_m_s'] >= row['earth_lift_speed_m_s']
        assert keeps == (row['com_m'] >= c_kept)
    assert kept['children_need_m'] == [2.76, 3.87] and kept['for_statures_m'] == [1.0, 1.4]


def test_guard_bands_follow_the_codes_carried_to_the_same_impact():
    bands = PRODUCT['guards']['bands']
    assert bands['kind'].startswith('design choice')
    rows = bands['rows']
    thresholds = {round(t['moon_drop_same_impact_m'], 2) for t in PRODUCT['unguarded_drops'].values()}
    assert rows[0]['drops_m'][1] == min(thresholds) and rows[1]['drops_m'] == [min(thresholds), max(thresholds)]
    assert rows[1]['options_m'] == [1.1, 2.2] and rows[2]['guard_m'] == 2.2
    assert rows[2]['drops_m'][0] == max(thresholds) == pytest.approx(4.71, abs=0.01)


def test_footholds_stand_beyond_the_reach_of_the_sized_hop():
    fh = PRODUCT['guards']['footholds']
    v = fh['take_off_m_s']

    def highest_at(x, base):
        # The highest a hop at speed v launched from height base reaches at distance x, searched over angles.
        best = -1e9
        for k in range(1, 900):
            th = math.radians(k / 10)
            vx, vy = v * math.cos(th), v * math.sin(th)
            t = x / vx
            best = max(best, base + vy * t - 0.5 * GM * t * t)
        return best

    floor_apex = v * v / (2 * GM)
    for spec in fh['reach_m'].values():
        s, d = spec['height_m'], spec['reach_m']
        assert highest_at(d * 0.95, s) > floor_apex > highest_at(d * 1.05, s)


def test_pushes_fall_with_gravity_and_impacts_do_not():
    p = PRODUCT['guards']['push_and_impact']
    assert p['standing_push_n']['earth'] / p['standing_push_n']['moon'] == pytest.approx(R, rel=0.01)
    assert p['impact_energy_j']['2 m/s'] == round(0.5 * lg.BODY_KG * 4)


def test_stairs():
    s = PRODUCT['stairs']
    for kind in ('everyday', 'bounding'):
        st = s[kind]
        assert st['kind'] == 'design choice'
        assert st['work_per_riser_j_kg'] == pytest.approx(GM * st['rise_m'], abs=1e-3)
        assert st['step_down']['landing_m_s'] == pytest.approx(math.sqrt(2 * GM * st['rise_m']), abs=0.006)
        assert st['step_down']['like_earth_step_m'] == pytest.approx(st['rise_m'] / R, abs=1e-4)
        for b in st['bounds']:
            assert b['take_off_up_m_s'] ** 2 / (2 * GM) == pytest.approx(b['rise_m'], rel=0.02)
    lo, hi = lg.EARTH_STAIR['general_access_rise_m']
    assert lo <= s['everyday']['rise_m'] <= hi
    assert s['everyday']['two_r_plus_g_m'] <= lg.EARTH_STAIR['two_r_plus_g_m'][1] + 1e-9
    assert s['flight_limits']['UK general access, 12 risers of 170 mm']['moon_rise_same_tumble_m'] == pytest.approx(
        12 * 0.17 * R, abs=0.06)
    assert s['flight_rise_choice_m'] <= 12 * 0.17 * R
    for kind in ('everyday', 'bounding'):
        landing = s['landings'][kind]
        need = landing['arrival_m_s'] ** 2 / (2 * 0.3 * GM) + 0.3      # half the grip, and a foot
        assert landing['needed_m'] == pytest.approx(need, abs=0.01) and need <= landing['depth_m']
    assert s['headroom_m']['choice_m'] == pytest.approx(2.0 + (R - 1) * 0.10, abs=0.01)


def test_stops_and_turns_take_six_times_the_room_and_twice_that_at_half_grip():
    for row in PRODUCT['stops'].values():
        assert row['moon_m'] / row['earth_m'] == pytest.approx(R, rel=0.07)          # rounding at 1 m/s
        assert row['moon_at_half_grip_m'] == pytest.approx(2 * row['moon_m'], abs=0.011)
        assert row['moon_turn_radius_m'] / row['earth_turn_radius_m'] == pytest.approx(R, rel=0.04)
        assert row['moon_turn_radius_at_half_grip_m'] == pytest.approx(2 * row['moon_turn_radius_m'], abs=0.011)


def test_handholds_hold_six_times_their_grip_in_body_mass():
    h = PRODUCT['handholds']
    for g in h['grips'].values():
        assert g['lunar_body_kg_one_hand'] == pytest.approx(g['kgf'] * R, abs=0.6)
        n = g['kgf'] * GE
        assert g['arrest_m']['moon_two_hands'] == pytest.approx(0.3 * (2 * n / (75 * GM) - 1), abs=0.006)
        assert g['arrest_m']['earth_two_hands'] < 0.15
    assert h['grips']["women's weak grip (T-score -2.5)"]['lunar_body_kg_one_hand'] > lg.BODY_KG


def test_ramps():
    for s in (1 / 12, 1 / 15, 1 / 20):
        theta = math.atan(s)
        for g in (GM, GE):
            run = 10.0
            v2 = 2 * g * math.sin(theta) * run
            assert v2 / (2 * lg.FRICTION * g) / run == pytest.approx(lg.runaway_run_out(s))
    r = PRODUCT['ramps']
    assert r['run_rise_for_the_landing_m'] == pytest.approx(0.5 * 0.6 * 3.0)
    for name, row in r['by_slope'].items():
        s = {'1:12': 1 / 12, '1:15': 1 / 15, '1:20': 1 / 20}[name]
        assert row['run_for_the_landing_m'] * math.sin(math.atan(s)) == pytest.approx(0.9, abs=0.01)
        # The stated departure: each run is longer than Approved Document K's going at its gradient.
        assert row['run_for_the_landing_m'] > row['adk_max_going_m']
        assert row['stop_on_the_ramp_from_walking_m']['moon'] / row['stop_on_the_ramp_from_walking_m']['earth'] == \
            pytest.approx(R, rel=0.05)
    assert r['lunar_slope_with_the_climb_of_earth_1_12_deg'] == pytest.approx(30.1, abs=0.1)
    for name, row in r['steep_paths'].items():
        theta = math.radians(float(name.split()[0]))
        assert row['share_of_grip_up_or_down'] == pytest.approx(math.tan(theta) / 0.6, abs=0.006)
    assert r['steepest_within_half_grip_deg'] == pytest.approx(math.degrees(math.atan(0.3)), abs=0.05)


def test_banks_hold_turns_at_half_the_grip():
    for key, deg in PRODUCT['movement_halls']['banking_deg'].items():
        v, r = (float(x.split()[0]) for x in key.split(' at '))
        b, k = math.radians(deg), v * v / (GM * r)
        # Friction needed on a floor banked at b for this turn, over the normal force.
        needed = (k * math.cos(b) - math.sin(b)) / (math.cos(b) + k * math.sin(b))
        assert needed == pytest.approx(0.3, abs=0.005)
    for key, radius in PRODUCT['movement_halls']['level_turn_radius_m'].items():
        v = float(key.split()[0])
        assert radius == pytest.approx(v * v / (0.3 * GM), abs=0.006)


def test_halls_are_provisional_and_sized_from_the_belt_jump():
    m = PRODUCT['movement_halls']
    assert m['kind'] == 'design choice, provisional'
    rise = (0.462 - 0.353) + R * (0.522 - (0.462 - 0.353))
    assert m['needed_m'] == pytest.approx(1.9 + rise + 1.0, abs=0.006)
    assert m['clear_height_m'] >= m['needed_m']
    assert min(m['needed_at_the_same_leg_work_m']) > m['clear_height_m']
    assert m['running_airtime_ratio'] == pytest.approx(math.sqrt(R), abs=0.005)
    assert m['running_hop']['length_m'] / m['running_hop']['earth_length_m'] == pytest.approx(R, rel=0.04)
    st = PRODUCT['stepping_stones']
    assert st['moon_length_m'] / st['earth_length_m'] == pytest.approx(R, rel=0.04)
