"""The seas' make-up through time: the product's binding and regeneration, its balances, and the states it holds."""
import hashlib
import json
import math

import pytest

from research.studies.sea_chemistry import run as S


@pytest.fixture(scope='module')
def product():
    return json.loads(S.OUT.read_text())


def close(a, b, rel=1e-3, abs_tol=1e-9):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_tol)


def find(product, delivery, sinks, corner='central', sediment='kept_open'):
    for r in product['runs']:
        if (r['delivery'], r['sinks'], r['corner'], r['sediment']) == (delivery, sinks, corner, sediment):
            return r
    raise KeyError((delivery, sinks, corner, sediment))


def test_product_binds_its_producer_and_inputs(product):
    assert product['schema'] == S.SCHEMA
    for path, h in {**product['producer']['files'], **product['producer']['inputs']}.items():
        assert hashlib.sha256((S.ROOT / path).read_bytes()).hexdigest() == h, path
    assert set(product['producer']['inputs']) == set(S.INPUTS)
    assert len(product['runs']) == 3 * len(S.SINKS) * 3 * len(S.SEDIMENT)


def test_a_fresh_run_matches_the_product(product):
    assert json.loads(S.dumps(S.results())) == product


def test_charge_balances_and_chloride_has_no_sink(product):
    """Every operation keeps the alkalinity equal to the conservative ions' charge, less the delivered seawater's own
    offset (its minor ions and borate), and chloride stays as delivered while the seas are open."""
    offset = dict(fresh=0.0, brackish=S.BRACKISH_SHARE, salted=1.0)
    sw = S.SEAWATER
    sw_offset = 1e3 * (2 * sw['ca'] + 2 * sw['mg'] + sw['na'] + sw['k'] - sw['cl'] - 2 * sw['so4'] - sw['alkalinity'])
    for r in product['runs']:
        if r['sediment'] != 'kept_open':
            continue
        for row in r['rows']:
            charge = 2 * row['ca'] + 2 * row['mg'] + row['na_k'] - row['cl'] - 2 * row['so4']
            scale = max(row['na_k'], row['alkalinity_meq_kg'], 1.0)
            assert abs(charge - row['alkalinity_meq_kg'] - offset[r['delivery']] * sw_offset) < 2e-4 * scale, \
                (r['delivery'], r['sinks'], row['year'])
            assert close(row['cl'], 1e3 * offset[r['delivery']] * sw['cl'], rel=1e-5, abs_tol=1e-9)


def test_the_first_thousands_of_years_follow_the_carbon_cycle(product):
    """The fresh seas reach three times calcite saturation when the carbon cycle's held runs find it, both taking the
    weathering that reaches the waters and the floors' uptake."""
    carbon = json.loads((S.ROOT / S.CARBON).read_text())
    for corner in ('low', 'central', 'high'):
        got = find(product, 'fresh', 'natural', corner)['milestones_yr']['calcite_in_the_seas']
        year = carbon['held'][f'{corner}/held/natural']['calcite_from_year']['sea']
        assert close(got, year, rel=0.03), corner
    rows = {row['year']: row for row in find(product, 'fresh', 'natural')['rows']}
    assert 5.5 < rows[0]['ph'] < 5.6 and 7.4 < rows[1000]['ph'] < 7.9


def test_held_seas_stay_where_the_delivery_set_them(product):
    for delivery, salinity, ph in (('fresh', 0.142, 8.32), ('brackish', 3.65, 8.29), ('salted', 35.1, 7.98)):
        rows = {row['year']: row for row in find(product, delivery, 'managed')['rows']}
        for year in (1e5, 1e7, 1e9):
            assert close(rows[year]['salinity_g_kg'], salinity, rel=0.01), (delivery, year)
            assert abs(rows[year]['ph'] - ph) < 0.02, (delivery, year)
        assert close(rows[1e8]['salinity_g_kg'], rows[1e9]['salinity_g_kg'], rel=1e-6)
        assert 'beyond_the_constants' not in find(product, delivery, 'managed')['milestones_yr']


def test_left_alone_the_seas_turn_alkaline_faster_without_clays(product):
    clays = find(product, 'fresh', 'natural')['milestones_yr']
    kept = find(product, 'fresh', 'natural_without_clays')['milestones_yr']
    assert clays['ph_9'] == kept['ph_9'] < 1e5
    assert kept['fresh_limit'] < clays['fresh_limit']
    assert kept['earth_salinity'] < clays['earth_salinity'] < clays['brine_ceiling']
    salted = find(product, 'salted', 'natural')['milestones_yr']
    assert salted['ph_9'] > clays['ph_9']
    for r in product['runs']:
        if r['sediment'] == 'kept_open' and r['sinks'] != 'managed':
            assert r['rows'][-1]['salinity_g_kg'] <= S.BRINE_CEILING_G_KG * 1.001


def test_sediment_fills_the_lakes_before_the_seas(product):
    sed = product['sediment']
    for e in (0, 1):
        assert sed['lakes_full_yr'][e] < sed['seas_full_yr'][e]
    for e in (0, 1):
        assert sed['seas_full_yr'][e] < sed['seas_full_without_deposits_yr'][e]
    for r in product['runs']:
        if r['sediment'] == 'settling':
            assert r['lakes_full_yr'] < r['seas_full_yr'] <= sed['seas_full_yr'][1]
            assert r['rows'][-1]['year'] == r['seas_full_yr']


def test_the_works_follow_weathering_and_phosphorus_follows_the_lakes(product):
    w = product['works_at_the_catchment_rate']
    sh = w['shares']
    assert close(w['minerals_mt_yr']['calcite'], 0.5 * w['weathering_eq_yr'] * sh['ca'] * 100.087 / 1e12, rel=1e-4)
    assert close(w['carbonate_burial_gt_c_yr'] + w['magnesite_burial_gt_c_yr'],
                 0.5 * w['weathering_eq_yr'] * (sh['ca'] + sh['mg']) * 12.011 / 1e15, rel=1e-4)
    years = product['report_years']
    for residence, cases in product['phosphorus_mg_m3'].items():
        for end in ('low', 'high'):
            fresh, salted = cases['fresh'][end], cases['salted'][end]
            assert fresh[-1] > fresh[years.index(1e5)]                      # the lakes fill and the seas get more
            assert all(s >= f for s, f in zip(salted, fresh))               # salt water frees the particles' share
            added = cases['fresh_phosphorus'][end]
            assert added[0] == S.EARTH_P_MG_M3 and close(added[-1], fresh[-1], rel=1e-3)
