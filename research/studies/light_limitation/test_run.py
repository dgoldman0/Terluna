"""The light limitation screen: its binding, its regeneration and the arithmetic it rests on."""
import hashlib
import json
import math

import pytest

from research.studies.light_limitation import run as L


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
    if not (L.ROOT / L.CLIMATOLOGY).is_file():
        pytest.skip('the GCM climatology has not been generated')
    assert json.loads(L.dumps(L.results())) == product


def test_the_rain_holds_the_land_below_its_light_in_every_band(product):
    """Rain-fed growth is the rain's share of the canopy's growth, below the light-limited potential everywhere, and
    the wet margins grow at the potential."""
    bands = product['land']['bands']
    for key, b in bands.items():
        assert max(b['rain_supports']) < 1.0, key
        for e in (0, 1):
            assert b['rain_fed_gt_c'][e] < b['potential_gt_c'][e], key
            share = min(b['light_potential_share_of_equator'], b['sunlight_cap_share_of_equator'][e])
            g = product['land']['totals']['equator_growth_g_c_m2_yr'][e]
            assert close(b['wet_margin_gt_c'][e], b['wet_margin_km2'][e] * 1e6 * g * share / 1e15), key
    shares = [b['rain_supports'][0] for b in bands.values()]
    assert all(a > b for a, b in zip(shares, shares[1:]))                 # drier poleward
    t = product['land']['totals']
    assert close(sum(t['wet_margin_km2']), sum(sum(b['wet_margin_km2']) for b in bands.values()))


def test_watered_land_is_held_by_the_light_or_by_the_sunlight_it_can_evaporate(product):
    for b in product['land']['bands'].values():
        for e in (0, 1):
            binds = b['sunlight_cap_share_of_equator'][e] < b['light_potential_share_of_equator']
            assert b['watered_limit'][e] == ('energy for transpiration' if binds else 'light')
            if binds:
                assert b['watered_elasticity'][e] == 1.0


def test_the_ledger_adds_up_and_weighs_the_light(product):
    led = product['ledger']
    places = led['places']
    for e in (0, 1):
        assert close(led['total_gt_c'][e], sum(p['gt_c'][e] for p in places.values()))
        assert close(sum(s[e] for s in led['share_of_production'].values()), 1.0)
    lo, hi = led['weighted_elasticity']
    assert 0.0 < lo < hi < 1.0
    assert places['rain_fed_land']['limit'] == 'water' and places['rain_fed_land']['elasticity'] == [0.0, 0.0]
    assert places['open_sea']['limit'] == 'nutrients'
    assert max(places['open_sea']['elasticity']) < min(places['lakes']['elasticity'])
    d = product['dimming']
    for e in (0, 1):
        assert close(d['total_percent'][e], d['gases_percent'] + d['people_heat_percent'][e] + d['glow_percent'])
        assert close(product['cost_through_light']['at_the_planned_dimming_percent'][e],
                     led['weighted_elasticity'][e] * d['total_percent'][e])


def test_the_rain_pathway_takes_earths_balanced_dimming_on_the_warming_in_the_air(product):
    r = product['rain_pathway']
    mean = sum(L.G1_SOLAR_REDUCTION) / len(L.G1_SOLAR_REDUCTION)
    assert close(r['g1_mean_solar_reduction_percent'], mean)
    assert close(r['precipitation_cut_per_percent_dimming'][0], L.G1_PRECIPITATION[0] / mean)
    d = product['dimming']
    assert close(r['dimming_answering_warming_in_the_air_percent'][1], d['gases_percent'] + d['glow_percent'])
    for e in (0, 1):
        assert close(r['growth_cut_at_the_planned_dimming_percent'][e],
                     r['growth_cut_per_percent_dimming_if_the_moon_matched'][e]
                     * r['dimming_answering_warming_in_the_air_percent'][e])
