"""What could take pressure off the shield's dimming: the product's binding and regeneration, its arithmetic, and the
README's figures against the product."""
import hashlib
import json
import math
from pathlib import Path

import pytest

from research.studies.heat_levers import run as R


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


def extract_present():
    if R.PINNED:
        return Path(R.PINS[R.EXTRACT]['read_from']).is_file()
    return (R.ROOT / R.EXTRACT).is_file()


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == R.SCHEMA
    for path, h in product['producer']['files'].items():
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == h, path
    assert set(product['producer']['inputs']) == set(R.INPUTS)
    pinned = product['producer']['pinned_inputs']
    assert set(pinned) == (set(R.PINS) if R.PINNED else set())
    for path, pin in pinned.items():
        assert pin['sha256'] == R.PINS[path]['sha256'] == product['producer']['inputs'][path], path
        assert pin['commit'] == R.PINS[path]['commit'] and pin['branch'] == 'domain/climate', path
    for path, h in product['producer']['inputs'].items():
        if path == R.EXTRACT and not extract_present():
            continue
        assert hashlib.sha256(R.raw(path)).hexdigest() == h, path
    for name, value in product['producer']['constants'].items():
        assert getattr(R.K, name) == value, name


@pytest.mark.skipif(not extract_present(), reason='the GCM extract (on the research drive, ignored) is not here')
def test_a_fresh_run_matches_the_product(product):
    assert json.loads(json.dumps(R.results())) == product


def test_the_targets_are_the_corrected_ones_with_the_resort_at_301_k(product):
    b = product['basis']
    assert R.INPUTS[R.TARGETS].endswith('/5') and R.INPUTS[R.EXTRACT].endswith('/5')
    assert b['resort_k'] == 301.0 and 0.5 < b['resort_share_of_the_way'] < 0.55
    d = b['targets_dimming_in_all_percent']
    assert d['hot world'][1] < d['resort, 301 K'][0] and d['resort, 301 K'][1] < d['lock, 295 K'][1]
    seas = product['gases']['warm_waters']['sea_median_c']
    assert round(seas['lock'], 2) == 24.12 and round(seas['hot'], 2) == 34.45
    assert seas['lock'] < seas['resort'] < seas['hot']


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
    # each saving carries its own trend to the hot world's ground, which falls less than the held pair's forcing
    low = o['low pair, 0.3 and 0.03 ppm']
    assert 0.97 < low['saving_at_the_hot_worlds_ground_factor'] < 0.99
    assert g['forcing_at_the_hot_worlds_ground_factor'] < low['saving_at_the_hot_worlds_ground_factor']
    assert math.isclose(low['saved_w_m2'][0], low['saved_w_m2'][1] * low['saving_at_the_hot_worlds_ground_factor'])
    held, lp = g['costs']['held'], g['costs']['low_pair']
    key = 'methane cut 0.9, N2O cut 0.95'
    assert all(a > h for a, h in zip(lp[key]['methane_air_share_per_year'], held[key]['methane_air_share_per_year']))
    # what the removal heats is small against what the lower levels save
    assert lp[key]['removal_heat_w_m2'][1] < 0.05 * low['saved_w_m2'][0]


def test_warm_waters_scale_the_wild_sources_alone(product):
    w = product['gases']['warm_waters']
    assert all(math.isclose(x, 5.89) for x in w['wild_sources_left_tg_yr'])      # 10% of the central 58.9 Tg
    for name in ('resort', 'hot'):
        v = w[name]
        f = v['source_factor_from_the_lock']
        for wild, added, rem, lock in zip(w['wild_sources_left_tg_yr'], v['wild_added_tg_yr'], v['removal_tg_yr'],
                                          w['lock']['removal_tg_yr']):
            assert math.isclose(added, wild * (f - 1.0)) and math.isclose(rem, lock + added)
    assert 1.8 < w['resort']['source_factor_from_the_lock'] < 2.0 < 3.4 < w['hot']['source_factor_from_the_lock'] < 3.6
    assert 15.0 < w['hot']['removal_tg_yr'][0] < w['hot']['removal_tg_yr'][1] < 18.5


def test_co2_saves_little_and_costs_growth_far_beyond_dimming(product):
    c = product['co2']
    lv = c['levels']['280 ppm']
    assert c['column_ground_k'] == 290.0
    assert round(lv['clear_sky_infrared_w_m2'], 2) == 1.99 and round(c['per_doubling_clear_w_m2'], 2) == 3.98
    # less CO2 absorbs less sunlight, which adds to the infrared it lets out
    assert c['sunlight_change_280_w_m2'] > 0.0
    assert math.isclose(lv['clear_sky_w_m2'], lv['clear_sky_infrared_w_m2'] + c['sunlight_change_280_w_m2'])
    assert lv['saved_w_m2'][0] < lv['saved_w_m2'][1] == lv['clear_sky_w_m2']
    assert round(lv['points'][0], 2) == 0.86
    assert round(lv['partial_pressure_pa'], 1) == 34.0
    for t in ('25C', '30C', '35C'):
        assert 0.6 < lv['leaf_rate_over_400_ppm'][t] < 0.75
    assert lv['leaf_rate_over_400_ppm']['35C'] < lv['leaf_rate_over_400_ppm']['25C']   # warmth deepens the loss
    dimming = product['shield']
    worst_dimming = dimming['growth_through_light_percent_per_point'][1] + dimming['land_rain_percent_per_point'][1]
    assert lv['growth_lost_per_point_percent'][0] > 10 * worst_dimming


def test_pressure_helps_only_where_the_moon_is_hot_and_less_under_cloud(product):
    p = product['pressure']
    t = p['by_target']
    assert t['lock']['with_the_plants_co2_kept_w_m2'] < 0.0 < t['resort']['with_the_plants_co2_kept_w_m2'] \
        < t['hot world']['with_the_plants_co2_kept_w_m2']
    assert t['resort']['ground_k'] == 301.0
    for name in ('resort', 'hot world', 'hot world settled'):
        assert max(t[name]['with_earths_cloud_masking_w_m2']) < t[name]['with_the_plants_co2_kept_w_m2']
    assert max(t['lock']['with_earths_cloud_masking_w_m2']) < 0.0
    assert 0.4 < t['hot world']['with_earths_cloud_masking_points'][0] < 0.5
    assert 0.9 < p['co2_kept_w_m2'] < 1.1
    assert 0.19 < p['n2_saved_share'] < 0.21


def test_cirrus_thinning_is_bounded_by_the_refits_high_cloud(product):
    c = product['cirrus']
    assert math.isclose(c['refit_high_cloud_net_w_m2'], c['hot_world_cloud_net_w_m2'] - c['design_case_cloud_net_w_m2'])
    assert 5.0 < c['refit_high_cloud_net_w_m2'] < 6.5
    assert c['thinning_w_m2'][0] == 0.0
    assert math.isclose(c['thinning_w_m2'][1], 0.4 * c['refit_high_cloud_net_w_m2'], rel_tol=1e-12)
    assert c['seeding']['per_year_kg'][0] < c['seeding']['per_year_kg'][1] < 1e9   # under a megatonne a year


def test_the_top_ups_fall_splits_at_the_exobase_across_its_range(product):
    t = product['top_up']
    floor = t['exobase_floor_w_m2']
    assert math.isclose(floor['low_flow'], t['first_millennium_w_m2'][0], rel_tol=1e-12)
    for f, a in zip(sorted(floor['high_flow']), t['above_the_exobase_at_the_high_end_w_m2'][::-1]):
        assert math.isclose(f + a, t['first_millennium_w_m2'][1], rel_tol=1e-12)
    for r, s in zip(t['exobase_radius_moon_radii'], t['first_screen_exobase_radius']):
        assert abs(r - s) < 0.01
    # braking to rest at the exobase removes the fall's own speed, whose energy is what lies above the exobase
    above = sorted(t['far_terminal_fall_mj_per_kg'] - x for x in t['exobase_fall_mj_per_kg'])
    for v, e in zip(t['arrival_speed_m_s'], above):
        assert math.isclose(0.5 * v ** 2, e * 1e6)
    assert 1.4e3 < t['arrival_speed_m_s'][0] < t['arrival_speed_m_s'][1] < 1.8e3
    assert min(t['reflecting_sail_lunar_discs']) > 100        # photons cannot take the inflow's momentum
    assert min(t['elevator_shortfall']) > 1e9


def test_beams_food_is_a_part_of_the_beams_and_its_carbon_outgrows_the_import(product):
    b = product['beams']
    for f, a in zip(b['food_w_m2'], b['all_beams_w_m2']):
        assert 0.0 < f < a
    assert max(b['food_in_orbit_w_m2']) < 0.01
    c = b['food_carbon']
    assert all(0.6 < s < 0.7 for s in c['moved_diet_share'])
    assert c['arriving_gt_c_yr'][0] > 0.44 and c['arriving_gt_c_yr'][0] > c['managed_import_gt_c_yr_at_10000'][1]
    assert math.isclose(c['times_the_managed_import'][1],
                        c['arriving_gt_c_yr'][1] / c['managed_import_gt_c_yr_at_10000'][0])


def test_the_glow_at_both_efficacies(product):
    g = product['glow']
    assert [round(e, 3) for e in g['efficacy']] == [0.645, 0.711]
    assert round(g['recommended_w_m2'][1], 3) == 4.237 and round(g['recommended_points'][1], 2) == 1.88
    assert round(g['the_lower_efficacy_takes_points'], 2) == 0.16
    for took, left in zip(g['the_roll_took_w_m2'], g['recommended_w_m2']):
        assert took > left > 0.0
    for v in g['set_aside_on_13_october'].values():
        for saved, left in zip(v['saved_w_m2'], g['recommended_w_m2']):
            assert 0.0 < saved < left
        assert v['points'][0] < v['points'][1]


def test_the_deep_air_halves_the_surface_albedo_lever_and_more_toward_the_poles(product):
    a = product['albedo']
    assert 0.3 < a['clear_sky_kernel'] < 0.4          # against Earth's 0.73
    assert a['all_sky_kernel'][0] < a['all_sky_kernel'][1] < a['clear_sky_kernel']
    k = sorted((x['lat_deg'], x['kernel']) for x in a['kernels_by_column'])
    assert k[0][1] == a['clear_sky_kernel']
    assert max(x for lat, x in k if lat > 50.0) < min(x for lat, x in k if 30.0 < lat < 40.0)
    assert a['cities_w_m2'][1] < a['parks_w_m2'][1]
    d = a['deserts']
    for v in d.values():
        for ratio in v['kernel_over_the_equators'].values():
            assert 0.75 < ratio[0] < ratio[1] < 0.9           # the equator's kernel reads the deserts 10-25% high
        assert v['w_m2'][1] < v['with_the_equators_kernel_w_m2'][1]
    assert min(d['lock']['desert_share_of_moon'].values()) > max(d['hot world']['desert_share_of_moon'].values())
    assert d['lock']['w_m2'][1] > d['hot world']['w_m2'][1]
    assert a['lock_deserts_cloud_cover'] < 0.75 * R.CLOUD_COVER['lock']


def pair_text(v, n=2):
    return f'{v[0]:.{n}f}–{v[1]:.{n}f}'


def test_the_readme_gives_the_products_figures(product):
    text = R.README.read_text()
    rows = {int(line.split('|')[1]): line for line in text.splitlines() if line.startswith('| ') and
            line.split('|')[1].strip().isdigit()}
    r = product['ranking']
    for lever in r['levers']:
        assert lever['rank'] in rows, lever['lever']
    figures = {1: [pair_text(r['levers'][0]['points'])], 2: [pair_text(r['levers'][1]['points'])],
               3: [pair_text(r['levers'][2]['points'])], 4: [pair_text(r['levers'][3]['points'])],
               5: [f"0 to {r['levers'][4]['points'][1]:.1f}"],
               6: [f"0 to {pair_text(r['details']['top_up_high_end_across_the_exobase_points'], 1)}"],
               7: [pair_text(r['details']['pressure_hot_world_clear_sky_points']),
                   pair_text(r['details']['pressure_hot_world_masked_points']),
                   f"−{abs(r['details']['pressure_lock_points'][0]):.2f}"],
               8: [pair_text(r['levers'][7]['points'])], 9: [pair_text(r['levers'][8]['points'])]}
    for rank, wanted in figures.items():
        for figure in wanted:
            assert figure in rows[rank], (rank, figure)
    four = r['first_four_points']
    assert f'The first four together save {pair_text(four, 1)} points' in text
    for target, share in r['first_four_share_of_dimming_in_all'].items():
        assert f'{100 * share[0]:.0f}–{100 * share[1]:.0f}%' in text, target
    for target, d in product['basis']['targets_dimming_in_all_percent'].items():
        assert f'{pair_text(d, 1)}%' in text, target
