"""What could take pressure off the shield's dimming: the product's binding and regeneration, and its arithmetic."""
import hashlib
import json
import math

import pytest

from research.studies.heat_levers import run as R


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


def extract_present():
    return (R.ROOT / R.EXTRACT).is_file()


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == R.SCHEMA
    for path, h in product['producer']['files'].items():
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == h, path
    assert set(product['producer']['inputs']) == set(R.INPUTS)
    for path, h in product['producer']['inputs'].items():
        if path == R.EXTRACT and not extract_present():
            continue
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == h, path
    for name, value in product['producer']['constants'].items():
        assert getattr(R.K, name) == value, name


@pytest.mark.skipif(not extract_present(), reason='the GCM extract (research/runs, ignored) is not here')
def test_a_fresh_run_matches_the_product(product):
    assert json.loads(json.dumps(R.results())) == product


def test_points_are_shares_of_the_hot_worlds_absorbed_sunlight(product):
    absorbed = product['basis']['absorbed_sunlight_w_m2']
    assert round(absorbed, 1) == 225.4
    g = product['gases']
    assert math.isclose(g['held_pair_points'], 100.0 * g['held_pair_allsky_w_m2'] / absorbed, rel_tol=1e-12)
    t = product['top_up']
    for w, p in zip(t['first_millennium_w_m2'], t['first_millennium_points']):
        assert math.isclose(p, 100.0 * w / absorbed, rel_tol=1e-12)


def test_lower_gas_pairs_save_in_order_and_cost_more_air(product):
    g = product['gases']
    o = g['options']
    saved = {k: v['saved_w_m2'][1] for k, v in o.items()}
    assert (saved['lowest grid pair, 0.1 and 0.01 ppm'] > saved['low pair, 0.3 and 0.03 ppm']
            > saved['methane low, 0.3 and 0.1 ppm'] > saved['N2O low, 1.87 and 0.03 ppm'] > 0.0)
    assert math.isclose(saved['low pair, 0.3 and 0.03 ppm'], 3.58 - 1.75, abs_tol=0.01)
    assert 0.9 < g['forcing_at_the_hot_worlds_ground_factor'] < 1.0
    held, low = g['costs']['held'], g['costs']['low_pair']
    key = 'methane cut 0.9, N2O cut 0.95'
    assert all(lo > h for lo, h in zip(low[key]['methane_air_share_per_year'], held[key]['methane_air_share_per_year']))
    # what the removal heats is small against what the lower levels save
    assert low[key]['removal_heat_w_m2'][1] < 0.05 * o['low pair, 0.3 and 0.03 ppm']['saved_w_m2'][0]
    w = g['warm_waters']
    assert 3.5 < w['source_factor_lock_to_hot'] < 4.2 and 1.0 < w['source_factor_lock_to_resort'] < w['source_factor_lock_to_hot']
    assert all(h > 5 * lo for h, lo in zip(w['held_methane_removal_tg_yr']['hot'], w['held_methane_removal_tg_yr']['lock']))


def test_co2_saves_little_and_costs_growth_far_beyond_dimming(product):
    c = product['co2']
    lv = c['levels']['280 ppm']
    assert round(lv['clear_sky_w_m2'], 2) == 1.99 and round(c['per_doubling_clear_w_m2'], 2) == 3.98
    assert lv['saved_w_m2'][0] < lv['saved_w_m2'][1] <= lv['clear_sky_w_m2']
    assert round(lv['partial_pressure_pa'], 1) == 34.0
    for t in ('25C', '30C', '35C'):
        assert 0.6 < lv['leaf_rate_over_400_ppm'][t] < 0.75
    assert lv['leaf_rate_over_400_ppm']['35C'] < lv['leaf_rate_over_400_ppm']['25C']   # warmth deepens the loss
    dimming = product['shield']
    worst_dimming = dimming['growth_through_light_percent_per_point'][1] + dimming['land_rain_percent_per_point'][1]
    assert lv['growth_lost_per_point_percent'][0] > 10 * worst_dimming


def test_pressure_helps_only_where_the_moon_is_hot(product):
    p = product['pressure']
    t = p['by_target']
    assert t['lock']['with_the_plants_co2_kept_w_m2'] < 0.0 < t['resort']['with_the_plants_co2_kept_w_m2'] \
        < t['hot world']['with_the_plants_co2_kept_w_m2']
    assert 0.9 < p['co2_kept_w_m2'] < 1.1
    assert 0.19 < p['n2_saved_share'] < 0.21


def test_cirrus_thinning_is_bounded_by_the_refits_high_cloud(product):
    c = product['cirrus']
    assert math.isclose(c['refit_high_cloud_net_w_m2'], c['hot_world_cloud_net_w_m2'] - c['design_case_cloud_net_w_m2'])
    assert 5.0 < c['refit_high_cloud_net_w_m2'] < 6.5
    assert c['thinning_w_m2'][0] == 0.0
    assert math.isclose(c['thinning_w_m2'][1], 0.4 * c['refit_high_cloud_net_w_m2'], rel_tol=1e-12)
    assert c['seeding']['per_year_kg'][0] < c['seeding']['per_year_kg'][1] < 1e9   # under a megatonne a year


def test_the_top_ups_fall_splits_at_the_exobase(product):
    t = product['top_up']
    assert math.isclose(t['exobase_floor_w_m2'][1] + t['above_the_exobase_w_m2'][1], t['first_millennium_w_m2'][1],
                        rel_tol=1e-12)
    assert math.isclose(t['exobase_floor_w_m2'][0], t['first_millennium_w_m2'][0], rel_tol=1e-12)
    assert 1.7 < t['exobase_radius_moon_radii'] < 1.8
    assert min(t['reflecting_sail_lunar_discs']) > 100        # photons cannot take the inflow's momentum
    assert min(t['elevator_shortfall']) > 1e9


def test_beams_food_is_a_part_of_the_beams(product):
    b = product['beams']
    for f, a in zip(b['food_w_m2'], b['all_beams_w_m2']):
        assert 0.0 < f < a
    assert max(b['food_in_orbit_w_m2']) < 0.01


def test_the_glow_matches_climate_targets(product):
    g = product['glow']
    assert round(g['recommended_w_m2'], 3) == 4.237 and round(g['recommended_points'], 2) == 1.88
    assert g['the_roll_took_w_m2'] > g['recommended_w_m2'] > 0.0
    for v in g['set_aside_on_13_october'].values():
        assert 0.0 < v['saved_w_m2'] < g['recommended_w_m2']
    assert g['at_the_corrected_radiations_efficacy_0_59']['change_points'] < 0.0


def test_the_deep_air_halves_the_surface_albedo_lever(product):
    a = product['albedo']
    assert 0.3 < a['clear_sky_kernel'] < 0.4          # against Earth's 0.73
    assert a['all_sky_kernel'][0] < a['all_sky_kernel'][1] < a['clear_sky_kernel']
    assert a['cities_w_m2'][1] < a['parks_w_m2'][1]
    d = a['deserts']
    assert d['lock']['desert_share_of_moon'] > d['hot world']['desert_share_of_moon']
    assert d['lock']['w_m2'][1] > d['hot world']['w_m2'][1]
    assert abs(d['hot world']['desert_share_of_moon'] - 0.72 * 0.189) < 0.03
