"""The domain synthesis's ledgers: their binding, their regeneration and the arithmetic they rest on."""
import hashlib
import json
import math

import pytest

from research.studies.domain_synthesis import ledgers as L


@pytest.fixture(scope='module')
def product():
    return json.loads(L.OUT.read_text())


def close(a, b, rel=1e-2):
    return math.isclose(a, b, rel_tol=rel, abs_tol=1e-9)


def test_product_binds_its_producer_and_inputs(product):
    assert product['schema'] == L.SCHEMA
    for path, h in product['producer']['files'].items():
        assert hashlib.sha256((L.ROOT / path).read_bytes()).hexdigest() == h, path
    assert set(product['producer']['inputs']) == set(L.INPUTS)
    for path, h in product['producer']['inputs'].items():
        assert hashlib.sha256((L.ROOT / path).read_bytes()).hexdigest() == h, path


def test_a_fresh_run_matches_the_product(product):
    assert json.loads(L.dumps(L.results())) == product


def test_the_blends_food_replaces_the_power_heavy_mix(product):
    blend = json.loads((L.ROOT / 'research/studies/food_systems/results/food_blend.json').read_text())
    case = blend['blend'][str(L.MOON_BILLION)][L.RULE]
    for supply in ('fusion', 'array'):
        h = product['heat'][supply]
        assert h['heat_by_component_tw']['food'] == case['heat'][supply]['heat_tw']
        for e in (0, 1):
            assert close(h['heat_tw'][e], sum(v[e] for v in h['heat_by_component_tw'].values()))
            assert close(h['added_dimming_percent'][e], h['heat_tw'][e] * product['heat']['factors']
                         ['dimming_percent_per_tw'][e])
            assert h['heat_tw'][e] > h['as_placed_with_power_heavy']['heat_tw'][e]


def test_the_land_adds_up(product):
    land = product['land']
    assert close(land['kept_wild_km2'] + land['parks_km2'] + land['outside_kept_and_parks_km2'], land['dry_land_km2'])
    for case in land['with_heritage'].values():
        for e in (0, 1):
            used = land['settled_km2'][e] + land['blend_fields_km2'][e] + case['heritage_pasture_km2']
            assert close(case['used_km2'][e], used)
            assert close(case['wild_share_outside_kept'][e], 1.0 - used / land['outside_kept_and_parks_km2'])
    assert land['middle_wet_fields_km2'] == [0.0, 0.0]


def test_the_gases_follow_their_sources_and_sinks(product):
    m = product['gases']['methane']
    for case in [m['low'], m['high']] + m['central']:
        s = case['sources_tg_yr']
        assert close(s['total'], s['wild'] + s['farms'] + s['heritage_cattle'])
        assert close(case['steady_ppm'], s['total'] / case['soil_sink_tg_yr_per_ppm'])
        assert close(case['hold_at_reference_tg_yr'], s['total'] - case['soil_sink_tg_yr_per_ppm'] * m['reference_ppm'])
    assert m['low']['steady_ppm'] < m['central'][0]['steady_ppm'] < m['central'][1]['steady_ppm'] < m['high']['steady_ppm']
    n = product['gases']['n2o']
    for case in (n['low'], n['high']):
        assert case['hold_tg_n_yr'] == case['sources_tg_n_yr']['total']           # N2O has no sink to share it
        assert close(case['rise_ppb_per_year'], case['sources_tg_n_yr']['total'] / n['tg_n_per_ppb'])


def test_people_losses_are_dredged_back(product):
    p = product['phosphorus']
    for land in p['by_land_loss'].values():
        high = land['by_food_end'][1]
        assert close(high['to_lakes'], land['resources_corner']['to_lakes'])
        assert close(high['dredged_share_of_lake_input'], land['resources_corner']['dredged_share_of_lake_input'])
        for row in land['by_food_end']:
            assert row['dredged_to_farms'] == row['people_losses']


def test_the_sky_layers_stack(product):
    h = product['heights']
    bands = list(h['bands'].values())
    assert all(a['top_km'] == b['bottom_km'] for a, b in zip(bands, bands[1:]))
    assert h['storms']['cloud_top_km']['median'] < 30.0 < h['storms']['cloud_top_km']['p90']
    assert 35.0 < h['summit_crown_km'] < 36.0
