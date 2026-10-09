"""Updraft, rain, mixing, lightning and fire balances for the storms component."""
import json
import math

import pytest

from research.studies.floater_viability import storms as s


def winds():
    return json.loads((s.ROOT/s.INPUT_FILES[1]).read_text())


def test_full_body_starts_full_and_lower_fill_rises_higher():
    w = winds()
    h, p, t = w['height_km'], w['pressure_pa'], w['temperature_k']
    assert s.rise_before_full(1., 10., h, p, t) == pytest.approx(10.)
    a = s.rise_before_full(.8, 10., h, p, t)
    b = s.rise_before_full(.5, 10., h, p, t)
    assert 10 < a < b
    # Expansion at the found height closes on 1/f.
    t0, p0 = s._interp(10., h, t), s._interp(10., h, p)
    assert s._interp(a, h, t)/t0*p0/s._interp(a, h, p) == pytest.approx(1/.8, rel=1e-3)


def test_updraft_hold_and_rain_brake_are_consistent():
    hold = s.hold_against_updraft(5., 1.2, 1.6)
    assert hold == pytest.approx(.5*1.2*.47*25/1.6)
    minutes = s.rain_brake_minutes(5., 20., 1.2, 1.6)
    assert minutes*20/60 == pytest.approx(hold)
    assert s.rain_brake_minutes(5., 0., 1.2, 1.6) is None
    assert s.vent_share_to_hold(47., 47.) == pytest.approx(1-math.exp(-1))


def test_rain_surge_and_sink_force_balance():
    held = s.rain_surge(60., 30., 20.)
    assert held == pytest.approx(20.)
    v = s.residual_descent(held, 1.2, 1.6)
    assert .5*1.2*.47*v*v == pytest.approx(held*1.6)
    assert s.rain_surge(60., 30., 70.) == 0.


def test_mixing_limits_and_flammable_reach():
    y = s.mole_to_mass_fraction(.04)
    assert y == pytest.approx(.04*.002016/(.04*.002016+.96*.02896))
    exit_y = s.mole_to_mass_fraction(.98)
    a = s.flammable_reach(.01, .106, 1.204, exit_y)
    b = s.flammable_reach(.1, .106, 1.204, exit_y)
    assert b['reach_m'] == pytest.approx(10*a['reach_m'])
    # Back-substitute: the centreline mass fraction at the reach is the 4% limit.
    d_ef = .01*math.sqrt(.106/1.204)
    assert exit_y*d_ef/(.208*a['reach_m']) == pytest.approx(y)
    assert s.upper_limit_in_air(.2095) == pytest.approx(.761, abs=2e-3)
    assert s.upper_limit_in_air(.175) < s.upper_limit_in_air(.2095)


def test_strikes_and_fire_dose():
    r = s.strike_rate(.001, 100., 100.)
    assert r['capture_area_km2'] == pytest.approx(math.pi*.04)
    assert r['strikes_per_year'] == pytest.approx(.001*math.pi*.04)
    a = s.neighbour_dose(1000., 30., .1, 50.)
    b = s.neighbour_dose(1000., 120., .1, 50.)
    assert a['dose_j_m2'] == pytest.approx(b['dose_j_m2'])
    assert a['dose_j_m2'] == pytest.approx(.1*1000*120e6/(4*math.pi*2500))
    assert s.wet_skin_drying_energy(1.) == 2.45e6


def test_oxygen_ingress_and_scrub():
    a = s.oxygen_ingress_years(24.6, 3.6, 1e-14, 17500., 41.9)
    b = s.oxygen_ingress_years(49.2, 3.6, 1e-14, 17500., 41.9)
    assert b == pytest.approx(2*a)
    assert s.scrub_cost_share(.25, 1., .175, .98) == pytest.approx(2*.25*.175/.98)


def test_evaluate_uses_the_lightning_products_numbers():
    out = s.evaluate()
    assert json.loads(json.dumps(out, allow_nan=False)) == out
    assert out['lightning']['ground_share'] == pytest.approx(.2)
    assert out['lightning']['night_flashes'] == 0
    ring = json.loads((s.ROOT/s.INPUT_FILES[0]).read_text())
    assert out['storm_counts']['rain_max_mm_h'] == ring['storms']['rain_max_mm_h']


def test_venting_follows_the_air_density_and_the_gas_adiabat():
    # An isothermal exponential atmosphere: density falls as exp(-z/H).
    h = [0., 10., 20., 30., 40.]
    p = [1e5*math.exp(-z/50.) for z in h]
    t = [280.]*5
    assert s.vent_share_profile(10., 5., h, p, t) == pytest.approx(1-math.exp(-5/50.), rel=1e-9)
    adiabat = s.vent_share_profile(10., 5., h, p, t, kappa=s.GAS_KAPPA)
    assert adiabat == pytest.approx(1-math.exp(-5/50.)**(1-s.GAS_KAPPA), rel=1e-9)
    assert adiabat < s.vent_share_profile(10., 5., h, p, t)


def test_gas_on_its_adiabat_fills_higher():
    w = winds()
    h, p, t = w['height_km'], w['pressure_pa'], w['temperature_k']
    for f in (.5, .64, .8):
        cool = s.rise_before_full_adiabatic(f, 10., h, p)
        assert cool > s.rise_before_full(f, 10., h, p, t)
        p0 = s._log_pressure_at(10., h, p)
        assert s._log_pressure_at(cool, h, p) == pytest.approx(p0*f**(1/(1-s.GAS_KAPPA)), rel=1e-9)


def test_momentum_length_scales_as_speed_times_root_diameter():
    a = s.momentum_length(.01, 100., .106, 1.204, 1.606)
    b = s.momentum_length(1., 100., .106, 1.204, 1.606)
    c = s.momentum_length(.01, 400., .106, 1.204, 1.606)
    assert b['momentum_length_m'] == pytest.approx(10*a['momentum_length_m'])
    assert c['momentum_length_m'] == pytest.approx(2*a['momentum_length_m'])
    assert a['exit_speed_m_s'] == pytest.approx(math.sqrt(200/.106))
