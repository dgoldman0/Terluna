"""The carbon-cycle box model: its binding and regeneration, conservation, convergence and the README's sums."""
import hashlib
import json
import math

import pytest

from research.studies.carbon_cycle import run as C
from shared.provenance import constants_changed


@pytest.fixture(scope='module')
def product():
    return json.loads(C.OUT.read_text())


def close(a, b, rel=2e-3, abs_tol=1e-6):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_tol)


def test_product_binds_its_producer_and_inputs(product):
    assert product['schema'] == C.SCHEMA
    for path, h in product['producer']['files'].items():
        assert hashlib.sha256((C.ROOT / path).read_bytes()).hexdigest() == h, path
    assert set(product['producer']['inputs']) == set(C.INPUTS)
    for path, h in product['producer']['inputs'].items():
        assert hashlib.sha256((C.ROOT / path).read_bytes()).hexdigest() == h, path
    assert not constants_changed(product['producer']['constants'])


def test_a_fresh_run_matches_the_product(product):
    assert json.loads(C.dumps(C.results())) == product


def test_carbon_and_the_cations_are_conserved(product):
    for name, share in product['checks']['conservation_worst_share'].items():
        assert share < 1e-11, name


def test_halving_the_step_changes_little(product):
    c = product['checks']['convergence']
    assert c['held_worst_relative_change'] < 0.01
    assert all(v < 5.0 for v in c['free_threshold_years_change'].values())
    assert abs(c['grown_year']['step'] - c['grown_year']['half_step']) <= 2.0


def test_the_waters_equilibrate_fast_against_the_model(product):
    for box in product['checks']['exchange'].values():
        assert box['years'] < 100.0


def test_the_weathering_law_is_the_co2_studys(product):
    assert product['checks']['law']['cases'] == 12
    assert product['checks']['law']['worst_relative_difference'] < 5e-3


def test_carbonate_constants_at_25_c():
    f = C.fresh_constants(298.15)
    assert close(-math.log10(f['k1']), 6.352, 1e-3)
    assert close(-math.log10(f['k2']), 10.329, 1e-3)
    assert close(math.log10(f['kh']), -1.468, 1e-3)
    assert close(-math.log10(f['ksp']), 8.48, 1e-3)
    assert close(-math.log10(f['kw']), 13.995, 1e-3)
    s = C.sea_constants(298.15, 35.0)
    assert close(s['kh'], 0.0284, 2e-3)
    assert close(-math.log10(s['k1']), 5.847, 1e-3)
    assert close(-math.log10(s['k2']), 8.966, 1e-3)
    assert close(-math.log10(s['kb']), 8.5975, 1e-3)
    assert close(s['ksp'], 4.27e-7, 3e-3)


def test_fresh_water_starts_near_ph_5_6():
    c = C.fresh_constants(C.LAND_AIR_K)
    h, hco3, co3 = C.carbonate(0.0, c['kh'] * 400e-6 * 121590 / 101325, c, 1.0, 1.0)
    assert 5.5 < -math.log10(h) < 5.65
    assert close(hco3, h, 1e-3)                                   # pure water and CO2: bicarbonate balances H+


def test_held_runs_balance_the_air(product):
    """The import is the net sink less the kilns and the dredged return; the net sink is the sum of the sinks."""
    for name, h in product['held'].items():
        for year, a in h['at'].items():
            assert close(a['import_'], a['net_sink'] - a['kilns'] - a['dredged'], 1e-3, 2e-4), (name, year)
            assert close(a['net_sink'], sum(a['sinks'].values()), 1e-3, 2e-4), (name, year)


def test_kiln_power_follows_the_flux(product):
    for h in product['held'].values():
        for a in h['at'].values():
            co2 = a['kiln_gt_co2_yr']
            gw = (co2['calcite'] * C.CALCITE_GJ_T_CO2 + co2['magnesite'] * C.MAGNESITE_GJ_T_CO2) * 1e18 / C.YEAR_S / 1e9
            assert close(a['kiln_gw_minimum'], gw, 2e-3, 1e-3)
            best = C.BEST_KILN_GJ_T_LIME[1] / C.BEST_KILN_GJ_T_LIME[0]
            assert close(a['kiln_gw_best_kilns'][1], best * a['kiln_gw_minimum'], 2e-3, 1e-3)
            assert close(co2['calcite'] + co2['magnesite'], a['kilns'] * C.M_CO2 / C.M_C, 2e-3, 1e-4)


def test_the_builds_carbon_adds_up(product):
    for name, h in product['held'].items():
        b = h['build']
        assert close(b['carbon_on_the_moon'], sum(b['where'].values()), 1e-3), name
        assert close(b['carbon_on_the_moon'], b['at_seeding'] + b['imported_after_seeding'], 1e-3), name
        assert close(b['carbon_on_the_moon'], b['stocks_filled'] + b['taken_by_the_sinks'], 1e-3), name
        assert close(b['where']['air_co2'], product['setting']['co2_gt_c_at_400_ppm'], 1e-3)


def test_methane_holds_at_its_level(product):
    for name, run in product['runs'].items():
        series = run['series']
        late = [ppm for y, ppm in zip(series['year'], series['ch4_ppm']) if y >= 1000]
        if '/held/' in name and not name.endswith('/free'):
            assert all(close(x, C.HELD_METHANE_PPM, 1e-3) for x in late), name
        if '/soils/' in name and not name.endswith('/free'):
            assert 21.0 < late[-1] < 30.0, name                    # the soils' central steady state


def test_co2_stays_at_400_in_held_runs_and_falls_in_free_runs(product):
    for name, run in product['runs'].items():
        ppm = run['series']['co2_ppm']
        if name.endswith('/free'):
            assert ppm[-1] < ppm[0] and ppm[0] == pytest.approx(400.0)
        else:
            assert all(x == pytest.approx(400.0) for x in ppm)


def test_the_free_runs_order_by_corner(product):
    years = {c: product['questions']['holds'][c]['years_to'] for c in C.CORNERS}
    for key in C.THRESHOLDS_PA:
        assert years['high'][key] < years['central'][key] < years['low'][key], key
        assert all(key in years[c] for c in C.CORNERS)
    order = list(C.THRESHOLDS_PA)
    for c in C.CORNERS:
        assert [years[c][k] for k in order] == sorted(years[c][k] for k in order)


def test_the_magnesite_route_lowers_the_import(product):
    for c in C.CORNERS:
        natural = product['held'][f'{c}/held/natural']['at']['10000']
        steered = product['held'][f'{c}/held/magnesite']['at']['10000']
        assert steered['import_'] < natural['import_']
        assert steered['kilns'] > natural['kilns']


def test_salted_seas_arrive_with_earths_carbonate_system(product):
    series = product['runs']['central/held/natural/salted']['series']
    seas = series['waters']['seas']
    assert 7.9 < seas['ph'][0] < 8.1
    assert 3.5 < seas['omega_calcite'][0] < 5.5
    assert 1.9 < seas['dic_mmol_kg'][0] < 2.2


def test_oxygen_is_burial_less_the_rock(product):
    o2 = product['setting']['o2_kg']
    for c in C.CORNERS:
        for year, o in product['questions']['oxygen'][c].items():
            assert close(o['burial_less_rock_share'], (o['burial_kg'] + o['rock_kg']) / o2, 2e-3), (c, year)
            assert o['burial_kg'] > 0 > o['rock_kg']


def test_the_corners_land_is_the_dry_land(product):
    ledger = json.loads((C.ROOT / C.LEDGERS).read_text())['land']
    for c in C.CORNERS:
        land = product['corners'][c]['land']
        assert close(sum(land['classes_km2'].values()), ledger['dry_land_km2'], 1e-3)
        assert close(sum(land['regions_km2'].values()), ledger['dry_land_km2'], 1e-3)


def test_lime_and_the_import_rate_follow_their_fluxes(product):
    for h in product['held'].values():
        for a in h['at'].values():
            co2 = a['kiln_gt_co2_yr']
            assert close(a['lime_gt_yr']['cao'], co2['calcite'] * C.M_CAO / C.M_CO2, 2e-3, 1e-4)
            assert close(a['lime_gt_yr']['mgo'], co2['magnesite'] * C.M_MGO / C.M_CO2, 2e-3, 1e-4)
            assert close(a['import_kg_c_s'], a['import_'] * 1e12 / C.YEAR_S, 2e-3, 1.0)


def test_the_methane_hold_passes_its_share_of_the_air(product):
    per_ppm = product['setting']['tg_ch4_per_ppm']
    for name, h in product['held'].items():
        m = h['methane_after_grown']
        flows = m['cut_at_the_sources'] + m['soils'] + m['removed'] + m['stock_growth']
        assert close(m['made'], flows, 2e-3, 0.01), name
        if '/held/' in name:
            assert close(m['air_share_passed_a_year'], m['removed'] / (per_ppm * C.HELD_METHANE_PPM), 2e-3)
            assert 0.009 < m['air_share_passed_a_year'] < 0.016             # the working case's 1.0-1.4%
