"""Force balances, scaling laws and product bindings for the structure component."""
import json
import math

import pytest

from research.studies.floater_viability import structure as s


def test_madison_curve_values_and_cap():
    # Wood (1951): the hyperbola passes 100% near 7.5 minutes and 62% near 10 years.
    assert s.duration_of_load(7.5*60/s.YEAR_S) == pytest.approx(1., abs=1e-3)
    assert s.duration_of_load(1/s.HOURS_PER_YEAR/60) == 1.
    ten = s.duration_of_load(10.)
    assert ten == pytest.approx((108.4/(10*s.YEAR_S)**.04635+18.3)/100)
    assert ten == pytest.approx(.621, abs=2e-3)
    assert s.duration_of_load(300.) < s.duration_of_load(100.) < ten < 1
    assert s.duration_of_load(1e12) > .183


def test_assembled_allowable_is_the_product_of_its_chain():
    a = s.assembled_allowable(20e6, .7, 10., 1.5)
    assert a['allowable_pa'] == pytest.approx(20e6*.7*s.duration_of_load(10.)/1.5)
    assert a['allowable_pa'] == pytest.approx(5.79e6, rel=2e-3)
    with pytest.raises(ValueError):
        s.assembled_allowable(20e6, 1.2, 10., 1.5)
    with pytest.raises(ValueError):
        s.assembled_allowable(20e6, .7, 10., .9)


def test_equatorial_resultant_matches_hemisphere_integral():
    r, d, g, base = 100., 1.098, 1.606, 20.
    n = 400
    total = 0.
    for i in range(n):  # rings of the upper hemisphere's horizontal projection
        rho = (i+.5)*r/n
        z = math.sqrt(r*r-rho*rho)  # height above the equator
        p = base + d*g*(r+z)
        total += p*2*math.pi*rho*r/n
    assert s.equatorial_resultant(r, base, d, g) == pytest.approx(total, rel=1e-4)


def test_tendon_share_grows_with_radius_over_breaking_length():
    d, g, sig, rho = 1.098, 1.606, 28e6, 1500.
    a = s.round_tendon_share(100., d, g, sig, rho)['tendon_share_of_lift']
    b = s.round_tendon_share(200., d, g, sig, rho)['tendon_share_of_lift']
    assert b == pytest.approx(2*a)
    assert a == pytest.approx(5*math.pi/4*rho*g*100./sig)
    radius = s.round_ceiling_radius(.25, d, g, sig, rho, 145.)
    assert s.round_tendon_share(radius, d, g, sig, rho, 145.)['tendon_share_of_lift'] == pytest.approx(.25)


def test_end_fittings_and_suspension_are_small_shares():
    assert s.end_fitting_share(.05) == pytest.approx(.1/math.pi)
    assert s.suspension_share(100., 28e6, 1500., 1.606) == pytest.approx(100*1500*1.606/28e6)
    assert s.seam_mass_share(.05, 2.) == .025


def test_module_optimum_is_the_minimum_of_skin_per_volume():
    args = (.1875, 5.8e6, 1500., 1.098, 1.606, 20.)
    best = s.module_optimum(*args)
    r = best['optimum_radius_m']
    f = lambda x: s.skin_per_volume(x, *args)
    assert f(r) == pytest.approx(best['skin_kg_per_m3'])
    assert f(r) < f(.9*r) and f(r) < f(1.1*r)


def test_hex_raft_and_colony_ties_do_not_depend_on_colony_width():
    raft = s.hex_raft(25., 1.098)
    assert raft['cell_over_disc'] == pytest.approx(2*math.sqrt(3)/math.pi)
    assert raft['gross_lift_kg_per_projected_m2'] == pytest.approx(1.098*4/3*math.pi*25/(2*math.sqrt(3)))
    tie = s.colony_network(125., 50., 28e6, 1500.)
    assert tie == pytest.approx(1500*125*50/28e6)
    assert s.colony_network(250., 50., 28e6, 1500.) == pytest.approx(2*tie)


def test_quilt_parts_add_up():
    q = s.quilt_structure(40., 30., 200., 5.8e6, 28e6, 1500., .15, 25e-6)
    assert q['total_kg_m2'] == pytest.approx(q['faces_kg_m2']+q['ties_kg_m2']+q['wall_film_kg_m2'])
    assert q['wall_area_per_projected_m2'] == pytest.approx(80/30)


def test_relaxation_time_halves_the_relative_wind():
    r, cd, u0 = 50., .47, 14.
    tau = s.relaxation_time(r, cd, u0)
    k = 3*cd/(8*1.5*r)
    u, t, dt = u0, 0., tau/20000
    while t < tau:
        u -= k*u*u*dt
        t += dt
    assert u == pytest.approx(u0/2, rel=1e-3)


def test_design_gusts_read_the_ring_and_apply_the_factor():
    winds = json.loads((s.ROOT/s.INPUT_FILES[0]).read_text())
    ring = json.loads((s.ROOT/s.INPUT_FILES[1]).read_text())
    rows = {r['height_km']: r for r in s.design_gusts(winds, ring)}
    ten = ring['wind_aloft_m_s']['10_km']
    assert rows[10]['design_horizontal_m_s'] == pytest.approx(1.4*(ten['max']-ten['median']))
    assert rows[10]['design_vertical_m_s'] == pytest.approx(1.4*ring['vertical_wind_m_s']['10_km']['up_max'])
    assert rows[10]['design_pressure_pa'] == pytest.approx(.5*rows[10]['air_density_kg_m3']*rows[10]['design_gust_m_s']**2)


def test_gusts_dominate_small_bodies_and_the_head_dominates_giants():
    small = s.gust_stress_ratio(10., 1.098, 1.606, 125.)
    giant = s.gust_stress_ratio(1000., 1.098, 1.606, 125.)
    assert small['gust_share'] > .6 > .05 > giant['gust_share']
    assert small['gust_share']+small['steady_share'] == pytest.approx(1)


def test_tear_flow_scaling_and_compartment_rule():
    a = s.tear_outflow(2., .05, 100., .1)
    b = s.tear_outflow(4., .05, 400., .1)
    assert b['volume_m3_s'] == pytest.approx(8*a['volume_m3_s'])
    c = s.compartment_loss(20., 4., 4)
    assert not c['stays_aloft'] and c['minimum_compartments'] == 5
    assert s.compartment_loss(20., 4., 5)['stays_aloft']
    assert s.partition_walls(10., 8, .2) == pytest.approx(.8)


def test_evaluate_is_plain_json_with_the_reference_chain():
    out = s.evaluate()
    assert json.loads(json.dumps(out, allow_nan=False)) == out
    assert out['reference_film']['allowable_pa'] == pytest.approx(5.79e6, rel=2e-3)
    assert out['collagen_tendon']['allowable_pa'] == pytest.approx(10e6*.8/1.5)
    assert all(r['design_gust_m_s'] > r['ring_max_m_s']-r['ring_median_m_s'] for r in out['design_gusts'])
