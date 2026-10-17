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


def test_the_held_levels_take_the_deep_columns_pair_and_their_removal_adds_up(product):
    """The held pair's warming is the deep column's, below Earth's 2019 air and the central steady state; methane's
    removal is what the sources biology leaves exceed the soils' uptake at the held level, N2O's all that biology
    leaves, and each passes the air that carries it at that level."""
    h = product['gases']['held']
    forcing = json.loads((L.ROOT / 'atmosphere/radiative_convective/results/reduced_gases.json').read_text())
    pair = forcing['moon_at_design']['both'][f'{L.HELD_CH4_PPM:g}+{L.HELD_N2O_PPM:g}']
    assert close(h['forcing_w_m2']['allsky'], pair['allsky_w_m2'])
    assert h['forcing_w_m2']['allsky'] < h['beside']['earth_2019']['allsky_w_m2']
    assert h['beside']['earth_2019']['allsky_w_m2'] < min(h['beside']['central_steady_state']['allsky_w_m2'])
    m = product['gases']['methane']
    for cut, c in h['methane_by_source_cut'].items():
        for e in (0, 1):
            left = (1.0 - float(cut)) * m['central'][e]['sources_tg_yr']['total']
            k = m['central'][e]['soil_sink_tg_yr_per_ppm']
            assert close(c['sources_left_tg_yr'][e], left)
            assert close(c['removal_tg_yr'][e], max(0.0, left - k * L.HELD_CH4_PPM))
            assert close(c['atmosphere_share_per_year'][e], c['removal_tg_yr'][e] / (L.HELD_CH4_PPM * m['tg_per_ppm']))
    n = product['gases']['n2o']
    for cut, c in h['n2o_by_source_cut'].items():
        for e, case in enumerate((n['low'], n['high'])):
            left = (1.0 - float(cut)) * case['sources_tg_n_yr']['total']
            assert close(c['sources_left_tg_n_yr'][e], left)
            assert close(c['atmosphere_share_per_year'][e], left / (L.HELD_N2O_PPM * n['tg_n_per_ppb'] * 1e3))


def test_the_weighed_pairs_trade_dimming_for_removal_and_the_held_one_is_the_middle(product):
    """The dimming answers the warming with sunlight the Moon would have absorbed; from the low pair to Earth's 2019
    air the dimming grows while N2O's decomposition shrinks and starts later, and the held pair is the middle one."""
    g = product['gases']
    h, w = g['held'], g['pairs_weighed']
    pop = json.loads((L.ROOT / 'provisioning/results/population.json').read_text())
    absorbed = (1.0 - L.DESIGN_PLANETARY_ALBEDO) * pop['heat']['design_sunlight_w_m2']
    assert close(h['dimming_basis']['absorbed_sunlight_w_m2'], absorbed)
    assert close(h['added_dimming_percent'], h['forcing_w_m2']['allsky'] / absorbed * 100.0)
    assert (h['methane_ppm'], h['n2o_ppm']) == (w['middle']['methane_ppm'], w['middle']['n2o_ppm'])
    assert close(w['middle']['added_dimming_percent'], h['added_dimming_percent'])
    low, mid, earth = w['low'], w['middle'], w['earth_2019']
    assert low['added_dimming_percent'] < mid['added_dimming_percent'] < earth['added_dimming_percent']
    for cut in ('0.95', '0.5'):
        for e in (0, 1):
            assert low['n2o'][cut]['decomposition_tw'][e] > mid['n2o'][cut]['decomposition_tw'][e] > \
                earth['n2o'][cut]['decomposition_tw'][e]
            assert low['n2o'][cut]['years_to_reach_from_none'][e] < mid['n2o'][cut]['years_to_reach_from_none'][e] < \
                earth['n2o'][cut]['years_to_reach_from_none'][e]
    for cut in ('0.9', '0.5'):
        assert mid['methane'][cut] == earth['methane'][cut]                # the same methane
        assert all(a > b for a, b in zip(low['methane'][cut]['atmosphere_share_per_year'],
                                         mid['methane'][cut]['atmosphere_share_per_year']))


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


def test_the_glow_is_dimmed_as_provisioning_gives_it(product):
    population = json.loads((L.ROOT / 'provisioning/results/population.json').read_text())['heat']
    glow = product['heat']['fleet_glow']
    for e in (0, 1):
        assert close(glow['added_dimming_percent'][e], glow['w_m2'][e] / population['design_sunlight_w_m2'] * 100)
        for supply in ('fusion', 'array'):
            total = product['heat'][f'{supply}_with_glow']['added_dimming_percent'][e]
            people = product['heat'][supply]['with_close_computing']['added_dimming_percent'][e]
            assert close(total, people + glow['added_dimming_percent'][e])


def test_towns_meet_reefs_by_the_distributions_rule(product):
    towns = product['heights']['sky_towns']
    assert close(towns['towns'] * towns['people_per_town'], L.ALOFT_BILLION * 1e9)
    for cover, reefs in towns['reef_passes_per_town_day'].items():
        small = towns['passes_per_town_day_if_100_m_bodies_made_the_cover'][cover]
        assert L.reef_passes(float(cover), 1.0, towns['relative_drift_m_s'], 1.0, widths=(0.1,))[1] == \
            pytest.approx(small[1], rel=1e-2)
        assert reefs[1] < small[0]


def test_array_power_divides_the_array_heat(product):
    a, h = product['array_power'], product['heat']['array']
    for e in (0, 1):
        total = h['with_close_computing']['heat_tw'][e]
        assert close(a['sky_towns']['intake_tw'][e] + a['ground']['intake_tw'][e], total)
        assert close(a['sky_towns']['w_m2'][e], a['sky_towns']['intake_tw'][e] * 1e12 / (a['sky_towns']['tops_km2'] * 1e6))
    assert a['beams']['5.8_ghz']['receiver_m_from_a_tile'][1] < a['beams']['2.45_ghz']['receiver_m_from_a_tile'][0]
