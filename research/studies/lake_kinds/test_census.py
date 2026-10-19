"""The lake census: its binding, its agreement with geography's products, and its arithmetic."""
import hashlib
import json

import numpy as np
import pytest

from research.studies.lake_kinds import census as CE


@pytest.fixture(scope='module')
def product():
    return json.loads(CE.OUT.read_text())


def present(path):
    return (CE.ROOT / path).is_file()


def digest(path):
    return hashlib.sha256((CE.ROOT / path).read_bytes()).hexdigest()


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == CE.SCHEMA
    for path, h in product['producer']['files'].items():
        assert digest(path) == h, path
    for path, h in product['producer']['inputs'].items():
        if (path.startswith(('research/runs/', 'geography/inputs/', 'geography/products/', 'climate/gcm/products/'))
                and not present(path)):
            continue                                     # LOLA, the grids and the climatology are fetched or rebuilt
        assert digest(path) == h, path


def test_the_census_is_the_drainage_summarys_lake_set(product):
    t = product['census']['total']
    s = t['drainage_summary']
    assert t['bodies'] == s['bodies'] == 27865
    # the summary's closed bodies are those with any cell in a closed depression; the census's closed lakes are the
    # mostly closed ones, and the rest are open lakes with nested pits
    assert t['bodies_with_any_closed_cell'] == s['closed_bodies'] == 2495
    assert t['closed_bodies'] + t['open_bodies_with_nested_closed_pits']['bodies'] == s['closed_bodies']
    assert t['wholly_closed_bodies'] <= t['closed_bodies']
    assert t['open_bodies_with_nested_closed_pits']['pits_km2'] < 0.1 * t['open_bodies_with_nested_closed_pits'][
        'area_km2']
    assert t['volume_km3'] == pytest.approx(s['volume_km3'], rel=1e-3)
    assert t['median_local_depth_m'] == pytest.approx(s['median_depth_m'], abs=1.0)


def test_every_body_falls_in_one_class(product):
    c = product['census']
    t = c['total']
    rows = c['by_deepest']
    assert sum(r['bodies'] for r in rows) == t['bodies']
    assert sum(r['area_km2'] for r in rows) == pytest.approx(t['area_km2'], rel=1e-3)
    assert sum(r['closed_bodies'] for r in rows) == t['closed_bodies']
    for r in rows:
        assert sum(x['bodies'] for x in r['by_band'].values()) == r['bodies']
        assert sum(x['bodies'] for x in r['by_size'].values()) == r['bodies']
        assert sum(r['floor_by_local_depth_km2']) == pytest.approx(r['area_km2'], rel=2e-3, abs=0.5)
        for depth, floor in r['floor_within_km2'].items():
            assert floor <= r['area_km2'] * (1.0 + 1e-3) + 0.5, (r['deepest_m'], depth)
        if r['deepest_m'][1] <= 12.0:                    # a lake no deeper than 12 m has all its bed within 12 m
            assert r['floor_within_km2']['12'] == pytest.approx(r['area_km2'], rel=1e-3, abs=0.5)
    assert sum(v['bodies'] for v in c['by_size'].values()) == t['bodies']
    assert sum(v['bodies'] for v in c['by_band'].values()) == t['bodies']


def test_the_cumulative_counts_rise_to_the_whole(product):
    c = product['census']
    for band, rows in c['cumulative_bodies_and_area_by_deepest'].items():
        values = list(rows.values())
        assert all(a[0] <= b[0] and a[1] <= b[1] + 1e-6 for a, b in zip(values, values[1:])), band
    whole = list(c['cumulative_bodies_and_area_by_deepest']['all'].values())[-1]
    largest = max(r['deepest_m'][1] for r in c['by_deepest'] if r['bodies'])
    assert whole[0] == c['total']['bodies'] - sum(r['bodies'] for r in c['by_deepest']
                                                  if r['deepest_m'][0] >= CE.DEEPEST_EDGES_M[-2])
    assert largest > 1000.0


def test_groundwater_reproduces_stage_2(product):
    g = product['census']['groundwater']
    assert g['labels_as_saved']
    assert g['into_lakes_km3_yr'] == pytest.approx(g['product_into_lakes_km3_yr'], rel=1e-4)
    assert g['gaining'] == g['product_gaining'] and g['losing'] == g['product_losing']
    springs = json.loads(CE.SPRINGS.read_text())
    gross = springs['cases']['central']['states']['mature']['variants']['central']['moon']['into_lakes_km3_yr']
    assert g['gross_into_lakes_km3_yr'] == pytest.approx(gross, rel=1e-3)
    # the bodies on the 15 km grid take their lakes' gross inflow, a little under the grid's own, which counts the
    # lake cells the 16 pixel/degree bodies leave out
    assert 0.5 * g['gross_into_lakes_km3_yr'] < g['inflow_km3_yr']['resolved'] < 1.05 * g['gross_into_lakes_km3_yr']


def test_closed_parts_follow_the_sea_chemistrys_scaling(product):
    c = product['census']['closed']
    assert c['parts'] == product['census']['total']['closed_bodies']
    turnover, milestones, path = CE.soda_paths()
    assert c['turnover_of_the_waters_yr'] == pytest.approx(turnover, rel=1e-3)
    assert turnover == pytest.approx(1908.0, rel=0.01)                 # 13.8 million km3 against 7,240 a year
    factor = 10.0 / 1.5 / turnover                                    # the sea chemistry's 10 m lake at 1.5 m a year
    one = sorted(m['one_g_kg'] * factor for m in milestones.values())
    assert one[0] == pytest.approx(2500.0, rel=0.03) and one[-1] == pytest.approx(4400.0, rel=0.03)
    assert CE.years_to(path, 1.0) <= milestones['central']['one_g_kg'] * 1.0001
    # deeper closed parts take longer, and the balance with the land water's load is the slower of the two here
    rows = [r for r in c['by_deepest'] if r['parts'] >= 20]
    medians = [r['by_scaling']['years_to_1_g_kg']['central'] for r in rows]
    assert medians[-1] > medians[0]
    for r in rows:
        y = r['by_scaling']['years_to_1_g_kg']
        assert y['high'] <= y['central'] <= y['low']
        chem = r['by_scaling']['years_to_hold_a_chemocline']
        assert chem['0.05/0.003'] <= chem['0.01/0.03']


def test_the_springs_are_denser_than_the_lakes(product):
    for band, v in product['census']['spring_water'].items():
        assert v['salinity_g_kg'][0] > 2.0 * v['lake_salinity_g_kg']
        for t, x in v['at'].items():
            lo, hi = x['density_excess_kg_m3']
            assert 0.03 < lo < hi < 0.2, (band, t)
            k_lo, k_hi = x['warmth_that_cancels_it_k']
            assert 0.1 < k_lo < k_hi < 0.9, (band, t)


def test_renewal_below_a_floor(product):
    r = product['census']['groundwater']['renewal_below']
    for depth, rows in r.items():
        for row in rows:
            assert row['deepest_m'][1] > float(depth)
            q = row['renewed_share_a_year']
            assert 0.0 <= q[0] <= q[1] <= q[2]
            slow, fast = row['years_to_margin']['0.03'], row['years_to_margin']['0.003']
            assert fast['within_10_yr'] >= slow['within_10_yr']


def test_kinds_of_a_threshold_are_read_from_the_classes(product):
    rows = product['census']['by_deepest']
    within = sum(r['bodies'] for r in rows if r['deepest_m'][1] <= 12.0)
    assert within == product['census']['cumulative_bodies_and_area_by_deepest']['all']['12'][0]


def test_weighted_median():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert CE.weighted_median(x, np.array([1.0, 1.0, 1.0, 10.0])) == 4.0
    assert CE.weighted_median(x, np.ones(4)) == 2.0
