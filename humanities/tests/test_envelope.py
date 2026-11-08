"""The dance envelope regenerates, binds its inputs, and each figure follows from its stated physics, checked by an
independent route where one exists."""
import hashlib
import json
import math
from pathlib import Path

import pytest

from humanities.dance import envelope as env
from shared.constants import MOON_SURFACE_GRAVITY as GM, STANDARD_GRAVITY as GE
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(env.OUT.read_text())


def test_product_regenerates_and_binds_its_code_constants_and_air():
    assert PRODUCT == json.loads(json.dumps(env.results()))
    assert PRODUCT['schema'] == env.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16] == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_samozino_take_off_solves_its_own_equation_and_sits_between_the_brackets():
    for g in (GM, GE):
        v = env.samozino_take_off(g)
        assert v * v == pytest.approx(2 * env.H_PO_M * (env.F0_N_KG * (1 - v / (2 * env.V0_M_S)) - g), rel=1e-9)
    j = PRODUCT['jumps']
    assert 0.25 < j['earth']['rise_m'] < 0.40                       # a trained woman's countermovement jump
    assert (j['moon_earth_take_off']['take_off_m_s'] < j['moon_force_velocity']['take_off_m_s']
            < j['moon_same_work']['take_off_m_s'])
    assert j['moon_earth_take_off']['rise_m'] / j['earth']['rise_m'] == pytest.approx(GE / GM, rel=1e-2)
    # The exact integration of the same push gives a lower take-off on both worlds and a larger lunar gain.
    ratio_closed = j['moon_force_velocity']['take_off_m_s'] / j['earth']['take_off_m_s']
    ratio_ode = j['ode_check']['moon_m_s'] / j['ode_check']['earth_m_s']
    assert j['ode_check']['earth_m_s'] < j['earth']['take_off_m_s'] and ratio_ode > ratio_closed > 1.2


def test_inertia_sums_the_body_and_orders_the_poses():
    M, com, I = env.whole_body_inertia()
    assert M == pytest.approx(55.0, abs=0.02)
    assert 0.54 < com[2] / 1.64 < 0.57                                # standing centre of mass over stature
    I_up = env.whole_body_inertia(arms_up=True)[2]
    I_side = env.whole_body_inertia(arms_side=True)[2]
    assert I_up[1, 1] > I[1, 1] and I_side[2, 2] > 2 * I[2, 2]
    r = PRODUCT['rotations']['inertia_kg_m2']
    assert 7 < r['somersault_layout_arms_up'] < 14 and 0.4 < r['twist_arms_in'] < 1.2


def test_traction_and_balance_follow_their_closed_forms():
    t = PRODUCT['traction']['speeds'][2]
    assert t['build_or_stop_m']['moon'] == pytest.approx(2.0 ** 2 / (2 * env.FRICTION * GM), abs=0.01)
    assert t['tightest_turn_radius_m']['earth'] == pytest.approx(2.0 ** 2 / (env.FRICTION * GE), abs=0.01)
    b = PRODUCT['balance']
    for where, g in (('moon', GM), ('earth', GE)):
        tau = math.sqrt(b['com_height_m'] / g)
        small_angle = tau * math.acosh(15 / 2)                        # theta0 cosh(t/tau), linearised
        assert b['from_2_to_15_deg_s'][where] == pytest.approx(small_angle, rel=0.03)
    assert b['from_2_to_15_deg_s']['moon'] / b['from_2_to_15_deg_s']['earth'] == pytest.approx(math.sqrt(GE / GM), rel=0.01)


def test_drag_flight_reduces_to_the_vacuum_and_terminal_speeds_balance_weight():
    v = 2.0
    nodrag = env.drag_flight(v, GM, 1.4, 0.0)
    assert nodrag['apex_m'] == pytest.approx(v * v / (2 * GM), abs=2e-3)
    assert nodrag['airtime_s'] == pytest.approx(2 * v / GM, abs=2e-3)
    for row in PRODUCT['air']['jumps_with_drag']:
        assert row['moon']['apex_m'] <= row['moon_vacuum']['rise_m'] and row['moon']['landing_m_s'] < v + 0.3
    rho = PRODUCT['air']['density_kg_m3']['moon_sea_level']
    assert 1.38 < rho < 1.42 and 1.15 < PRODUCT['air']['density_kg_m3']['moon_summit'] < 1.19
    vt = PRODUCT['air']['terminal_speed_m_s']['large sleeves or a cloak spread']['moon']
    assert 0.5 * rho * vt ** 2 * 2.0 == pytest.approx(55.0 * GM, rel=0.02)


def test_fabric_falls_and_skirt_flare_balance_their_forces():
    rho = PRODUCT['air']['density_kg_m3']['moon_sea_level']
    for row in PRODUCT['fabric']['fabrics']:
        vm = row['face_on_fall_m_s']['moon']
        assert 0.5 * rho * env.CN_FLAT * vm ** 2 == pytest.approx(row['areal_density_kg_m2'] * GM, rel=0.01)
        assert row['air_over_weight_at_1_m_s']['moon'] / row['air_over_weight_at_1_m_s']['earth'] == pytest.approx(
            (rho / env.RHO_EARTH) * (GE / GM), rel=0.02)
    sk = PRODUCT['fabric']['skirt']
    r, l = sk['waist_radius_m'], sk['length_m']
    for row in sk['rigid_strip_flare']:
        th = math.radians(row['flare_deg'])
        w = row['rpm']['moon'] * 2 * math.pi / 60
        # A uniform strip hinged at the waist: centrifugal moment about the hinge equals gravity's.
        centrifugal = w * w * math.cos(th) * (r * l * l / 2 + math.sin(th) * l ** 3 / 3)
        gravity = GM * math.sin(th) * l * l / 2
        assert centrifugal == pytest.approx(gravity, rel=0.02)


def test_partnering_and_the_leap_stage():
    p = PRODUCT['partnering']
    assert p['throw_earth'] is None and p['one_arm_overhead_hold'] == {'moon': True, 'earth': False}
    a = env.THROW_FORCE_N / 55.0 - GM
    assert p['throw_moon']['take_off_m_s'] == pytest.approx(math.sqrt(2 * a * env.THROW_STROKE_M), abs=1e-3)
    for row in PRODUCT['stage_for_a_leap']['rows']:
        assert row['total_m'] == pytest.approx(row['approach_m'] + row['flight_m'] + row['run_out_m'], abs=0.02)


def test_aerial_figures_follow_their_mechanics():
    a = PRODUCT['aerials']
    I = a['inertia_kg_m2']
    # The corkscrew's tilt: twist rate L sin(tilt) / I_twist equal to the somersault's rate.
    lay = I['layout, arms gathered']
    assert a['corkscrew']['tilt_deg'] == pytest.approx(math.degrees(math.asin(lay['twist'] / lay['somersault'])), abs=0.1)
    # The hanging salto's programme turns the body exactly once.
    Io, It = I['layout, arms overhead']['somersault'], I['tuck']['somersault']
    for where in ('moon', 'earth'):
        h = a['hanging_salto'][where]
        T = a['flights_s'][where]
        turned = h['angular_momentum_N_m_s'] * T * (0.6 / Io + 0.4 / It)
        assert turned == pytest.approx(2 * math.pi, rel=0.01)
        assert h['hang_s'] == pytest.approx(0.4 * T, abs=0.01)
    # The Moon's take-off covers each figure (the flat spin from the high jump); Earth's standing take-off covers none
    # of the corkscrew, star wheel and flat spin.
    cap = a['take_off']
    for fig in ('corkscrew', 'star_wheel'):
        assert a[fig]['angular_momentum_N_m_s']['moon'] < cap['moon']['somersault_or_cartwheel_N_m_s']
        assert a[fig]['angular_momentum_N_m_s']['earth'] > cap['earth']['somersault_or_cartwheel_N_m_s']
    flat = a['flat_spin']['angular_momentum_N_m_s']
    assert flat['moon_high'] < cap['moon_high']['about_the_vertical_N_m_s'] and flat['moon'] < cap['moon']['about_the_vertical_N_m_s']
    assert flat['earth'] > cap['earth']['about_the_vertical_N_m_s']
    assert a['flat_spin']['spread_is_the_largest_axis']
    # The take-off's turning follows from its mean push.
    v = PRODUCT['jumps']['choreographic'][1]['take_off_m_s']
    t = 2 * env.H_PO_M / v
    assert cap['moon']['mean_push_N'] == pytest.approx(env.MASS_KG * (v / t + GM), abs=0.1)
    # Sleeves stand above the falling hands only on the Moon, and never further than their length.
    s = a['sleeve_flight']
    assert 0.5 < s['moon']['sleeves_above_hands_at_landing_m'] <= env.SLEEVE_LENGTH_M
    assert s['earth']['sleeves_above_hands_at_landing_m'] < 0.2 and s['moon']['settle_after_landing_s'] > 10 * s['earth']['settle_after_landing_s']
