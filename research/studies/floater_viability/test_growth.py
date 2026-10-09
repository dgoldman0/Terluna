"""Closed-form ages, host supply and phosphorus arithmetic for the growth component."""
import csv
import json
import math

import pytest

from research.studies.floater_viability import growth as g


def test_hydrogen_costs_in_carbon():
    assert g.glucose_per_hydrogen(4.) == pytest.approx(.180156/(.002016*4))
    assert g.fermentation_carbon_per_hydrogen(2.1) == pytest.approx(.4*.180156/(.002016*2.1))
    photo = g.photo_carbon_per_hydrogen(2., 200., .01)
    assert photo == pytest.approx(2*120e6/(.01*200*g.YEAR_S))
    assert g.photo_carbon_per_hydrogen(5., 200., .01) == pytest.approx(2.5*photo)


def test_construction_parts_sum():
    c = g.construction_per_area(1.5, 3., 1.1, .05, 4.)
    assert c['total_c_kg_m2'] == pytest.approx(6+3*.556+1.1*.556+.05)
    assert 0 < c['gas_share'] < 1


def test_round_ages_match_the_closed_form_for_constant_costs():
    # C and S constant per area: total = c pi R², rate = S pi R², t = (2c/S) ln(R/R0).
    c, srate = 10., 1.
    rows = [dict(radius_m=r, construction_c_kg_m2=c, surplus_c_kg_m2_year=srate)
            for r in [20+i*.5 for i in range(161)]]
    out = g.round_body_ages(rows, 20., (40., 100.))
    assert out['years_to_radius']['40.0'] == pytest.approx(2*c/srate*math.log(2), rel=1e-5)
    assert out['years_to_radius']['100.0'] == pytest.approx(2*c/srate*math.log(5), rel=1e-5)
    assert out['growth_stops_at_radius_m'] is None


def test_round_growth_stops_where_surplus_ends():
    rows = [dict(radius_m=r, construction_c_kg_m2=10., surplus_c_kg_m2_year=1-r/100) for r in (20., 50., 80., 120.)]
    out = g.round_body_ages(rows, 20., (60., 110.))
    assert out['years_to_radius']['110.0'] is None
    assert out['growth_stops_at_radius_m'] == pytest.approx(100.)


def test_flat_colonies_double_at_a_fixed_rate():
    out = g.flat_body_ages(10., 1., 100., (200., 800.))
    assert out['doubling_years'] == pytest.approx(10*math.log(2))
    assert out['years_to_area']['200.0'] == pytest.approx(out['doubling_years'])
    assert out['years_to_area']['800.0'] == pytest.approx(3*out['doubling_years'])
    assert g.flat_body_ages(10., 1., 100., (200.,), .5)['doubling_years'] == pytest.approx(2*out['doubling_years'])
    assert g.flat_body_ages(10., 0., 100., (200.,))['doubling_years'] is None


def test_survival_and_recruitment():
    assert g.survival(100., .01) == pytest.approx(math.exp(-1))
    assert g.recruitment_need(.005, .1) == pytest.approx(.05)


def test_host_crown_matches_the_megaforest_allometry():
    crown = g.host_crown(300.)
    assert crown['crown_radius_m'] == pytest.approx(54.)
    n = 2000
    area = sum(2*54*.3*math.sqrt(max(0., 4*s*(1-s)))*(.55*300)/n for s in [(i+.5)/n for i in range(n)])
    assert crown['crown_drag_area_m2'] == pytest.approx(.8*area, rel=1e-4)


def test_wind_moment_ratio_and_threshold():
    row = g.host_wind_moment_ratio(20., 300.)
    crown = g.host_crown(300.)
    assert row['moment_ratio'] == pytest.approx(.47*math.pi*400*300/(crown['crown_drag_area_m2']*crown['crown_centroid_m']))
    assert row['threshold_factor'] == pytest.approx(1/math.sqrt(1+row['moment_ratio']))


def test_nursery_fill_conserves_hydrogen():
    fill = g.nursery_fill_years(2000., 2.1, 5000., 300.)
    rate = 5000/g.fermentation_carbon_per_hydrogen(2.1) + 300
    assert fill['years']*rate == pytest.approx(2000.)
    assert g.nursery_fill_years(2000., 2.1, 0., 0.)['years'] is None


def test_phosphorus_arithmetic():
    need = g.phosphorus_net_need(1.1, .5, 2., .5, 3., .1, .3)
    assert need['net_need_g_m2_year'] == pytest.approx(1.1*.5*2*.5+3*.1*.3)
    assert g.rain_phosphorus(1., .02) == pytest.approx(1*g.JULIAN_YEAR_DAYS*.02/1000)
    assert g.capture_phosphorus(10., .005, 1., .5) == pytest.approx(10e-6*.005*.5*g.YEAR_S)


def test_evaluate_reads_the_canopy_and_reference_tree():
    out = g.evaluate()
    assert json.loads(json.dumps(out, allow_nan=False)) == out
    plant = json.loads((g.ROOT/g.INPUT_FILES[0]).read_text())
    eq = plant['results']['near']['0']['earth_like']
    assert out['canopy']['0']['npp_c_kg_m2_year'] == pytest.approx(eq['growth']*g.CYCLES_PER_YEAR/1000)
    with open(g.ROOT/g.INPUT_FILES[2], newline='') as f:
        row = next(r for r in csv.DictReader(f) if r['case_id'] == '300m-rigid-uniform-reference-30mm')
    assert out['reference_tree']['first_threshold_m_s'] == float(row['first_threshold_m_s'])
    assert out['fruit_share_precedent'] == [.3, .5, .6, .7, .8]


def test_nursery_fill_pays_the_juveniles_gas_loss():
    plain = g.nursery_fill_years(2000., 2.1, 3000., 300.)
    lossy = g.nursery_fill_years(2000., 2.1, 3000., 300., 25.)
    made = lossy['years']*(lossy['hydrogen_from_host_kg_year']+300.-25.)
    assert made == pytest.approx(2000.)
    assert lossy['years'] > plain['years']
    assert g.nursery_fill_years(2000., 2.1, 0., 10., 25.)['years'] is None
