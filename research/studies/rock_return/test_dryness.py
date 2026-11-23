"""Tests of the first look at the carriers' two ranges: the helpers' arithmetic, the stored product and the page's figures."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from research.studies.rock_return import carriers as CA
from research.studies.rock_return import dryness as DR
from research.studies.rock_return import evolve as EV
from research.studies.rock_return import run as RU

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


@pytest.fixture(scope='module')
def product():
    return json.loads(DR.OUT.read_text())


def close(a, b, rel=1e-9):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(close(a[k], b[k], rel) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(close(x, y, rel) for x, y in zip(a, b))
    return math.isclose(a, b, rel_tol=rel)


def test_growth_by_rain_rises_and_levels_off():
    assert float(DR.miami(0.0)) == 0.0
    assert math.isclose(float(DR.miami(1000.0)), 3000.0 * (1.0 - math.exp(-0.664)))
    g = DR.miami(np.array([50.0, 500.0, 5000.0]))
    assert g[0] < g[1] < g[2] < 3000.0 and g[1] / g[0] > 5.0 and g[2] / g[1] < 4.0


def test_a_shape_spread_over_the_land_keeps_the_mean():
    shape, area = np.array([0.0, 1.0, 3.0]), np.array([2.0, 1.0, 1.0])
    over = DR.over_the_mean(shape, area)
    assert math.isclose(float((over * area).sum() / area.sum()), 1.0)
    assert over[0] == 0.0 and math.isclose(over[2] / over[1], 3.0)


def test_food_is_a_larger_share_where_plants_grow_less():
    assert DR.food_multiple(146.0, 146.0, 1.0) == 1.0
    assert math.isclose(DR.food_multiple(73.0, 146.0, 0.05), 10.0)


def test_a_longer_haul_costs_hours_and_then_carriers():
    masses = dict(small=1.0, great=1000.0)
    rows = DR.haul_rows(3.0e12, 2.0e13, [1400.0, 1900.0], [15.0, 25.0], [0.2, 0.25], masses, dict(small=1.0e6, great=1.0e8))
    keys = list(rows)
    assert keys == ['9.4 km', '12.5 km', '31 km', '100 km']
    first, far = rows[keys[0]], rows[keys[-1]]
    assert math.isclose(first['hours_a_year'], CA.hours_a_year(1000.0, DR.HAULS_M[0], RU.CRUISE_M_S, RU.HANDLING_S))
    assert math.isclose(first['trips_in_the_first_hauls_hours'], 1000.0)
    assert close(first['by_size']['small']['at_the_same_trips'], first['by_size']['small']['in_the_same_hours'])
    assert far['hours_a_year'] > 5.0 * first['hours_a_year']
    assert all(x > y for x, y in zip(far['kj_per_kg_delivered'], first['kj_per_kg_delivered']))
    # a glide from the height gained covers the straight way home and not a hundred kilometres
    assert math.isclose(first['kj_per_kg_delivered'][0], CA.delivered_cost_j_per_kg(9400.0, 1000.0, 25.0, 0.25, 0.6, 1.0) / 1e3)
    assert math.isclose(far['kj_per_kg_delivered'][0], CA.delivered_cost_j_per_kg(1.0e5, 1000.0, 25.0, 0.25, 0.6, 0.25) / 1e3)
    for name, mass in masses.items():
        same, hours = far['by_size'][name]['at_the_same_trips'], far['by_size'][name]['in_the_same_hours']
        assert math.isclose(same['carriers'] * 1000.0 * DR.LOAD * mass, 3.0e12)
        assert math.isclose(hours['carriers'] * far['trips_in_the_first_hauls_hours'] * DR.LOAD * mass, 3.0e12)
        assert hours['carriers'] > same['carriers'] and hours['dry_food_mt_yr'][0] > same['dry_food_mt_yr'][0]


def test_the_dryness_product_is_whole(product):
    assert product['schema'] == DR.SCHEMA
    assert product['producer']['code_sha256'] == {p: EV.sha256(HERE / p) for p in ('dryness.py', 'rock.py', 'carriers.py', 'run.py')}
    for path, digest in product['inputs'].items():
        assert EV.sha256(ROOT / path) == digest, path
    by = product['by_dryness']
    rock = json.loads((HERE / 'results' / 'rock.json').read_text())
    assert math.isclose(product['rock_kg_yr'], rock['rock']['kg_yr_central'], rel_tol=1e-9)
    assert product['the_two_studies_water_agrees'] == 1.0
    assert math.isclose(sum(r['land_share'] for r in by.values()), 1.0, rel_tol=1e-9)
    assert math.isclose(sum(r['rock_shed_gt_yr'] for r in by.values()) * 1e12, product['rock_kg_yr'], rel_tol=1e-9)
    assert math.isclose(sum(sum(r['rock_settling_here_gt_yr'].values()) for r in by.values()) * 1e12, product['rock_kg_yr'], rel_tol=1e-9)
    assert math.isclose(np.sum(product['rock_from_class_to_class_gt_yr']) * 1e12, product['rock_kg_yr'], rel_tol=1e-9)
    for key in ('by_cover', 'by_rain'):                 # either shape spreads the land's mean and keeps it
        assert math.isclose(sum(r['land_share'] * r['plant_growth_over_the_lands_mean'][key] for r in by.values()), 1.0, rel_tol=1e-6)
    seas = json.loads(DR.SAND_SEAS.read_text())         # the ends' sand and the wind's sand are the sand seas study's
    for kind, km3 in seas['supply']['km3_yr'].items():
        assert math.isclose(sum(r['lasting_sand_at_ends_km3_yr'][kind] for r in by.values()), km3, abs_tol=1e-4)
    assert abs(sum(r['winds_sand_km2'] for r in by.values()) - product['sand']['winds_sand_km2']) <= 2
    for sea in product['sand']['made_seas'].values():
        assert math.isclose(sum(sea['at_ends_within_the_radius_km3_yr'].values()), sea['need_km3_yr'], rel_tol=0.02)
    for region in product['regions'].values():
        assert math.isclose(sum(region['classes'].values()), 1.0, rel_tol=1e-9)


def test_the_stored_hauls_are_what_the_code_gives(product):
    sky = RU.load(RU.SKY)
    sizes = dict(RU.SIZES_KG, albatross_scale=RU.flyer_of(sky)['mass_kg'])
    upkeep = {name: RU.FW.bird_fmr_j_day(mass) for name, mass in sizes.items()}
    rows = DR.haul_rows(product['rock_kg_yr'], product['land_m2'], product['plant_growth_land_mean_g_dry_m2_yr'],
                        sky['food_web']['ranges']['glide_ratio'], sky['food_web']['ranges']['muscle_efficiency'], sizes, upkeep)
    assert close(rows, product['hauls'], rel=1e-9)
    # at the straight way back the food is the sizing's, within the difference between one typical trip and the rock's own ways
    sizing = json.loads(RU.OUT.read_text())['carriers']['food_over_plant_growth']
    crow = rows['9.4 km']['by_size']['crow']['at_the_same_trips']['of_the_lands_plant_growth']
    assert 0.003 < crow[0] < crow[1] < 0.006 and isinstance(sizing, (dict, list))


def test_what_the_write_up_says_of_the_two_ranges(product):
    by, hauls, sand, regions = product['by_dryness'], product['hauls'], product['sand'], product['regions']
    assert [round(100 * by[k]['land_share'], 1) for k in ('hyper-arid', 'arid', 'semi-arid')] == [1.6, 20.5, 25.7]
    assert round(100 * sum(by[k]['land_share'] for k in ('dry sub-humid', 'humid', 'wet')), 1) == 52.2
    assert [round(by[k]['rock_shed_g_m2_yr']) for k in ('hyper-arid', 'arid', 'semi-arid')] == [68, 77, 114]
    assert [round(by[k]['rock_shed_g_m2_yr']) for k in ('dry sub-humid', 'humid', 'wet')] == [158, 184, 210]
    # the plants' growth over the land's mean, and the return's food against its share on the land's mean
    growth = {k: sorted(by[k]['plant_growth_over_the_lands_mean'].values()) for k in by}
    assert [round(x, 2) for x in growth['hyper-arid']] == [0.05, 0.09] and [round(x, 2) for x in growth['arid']] == [0.13, 0.34]
    assert [round(x, 1) for x in growth['semi-arid']] == [0.5, 1.0]
    food = {k: sorted(by[k]['food_over_the_lands_mean_share'].values()) for k in by}
    assert 3.8 < food['hyper-arid'][0] < 3.9 and 8.6 < food['hyper-arid'][-1] < 8.7
    assert 1.2 < food['arid'][0] < 1.3 and 4.1 < food['arid'][-1] < 4.2
    assert 0.7 < food['semi-arid'][0] and food['semi-arid'][-1] < 1.6
    assert all(0.65 < min(food[k]) and max(food[k]) < 1.25 for k in ('dry sub-humid', 'humid', 'wet'))
    # what that makes of the plants' growth, by the sizing's shares for the greatest and the smallest carrier
    sized = {name: hauls['9.4 km']['by_size'][name]['at_the_same_trips']['of_the_lands_plant_growth'] for name in ('albatross_scale', 'hornero')}
    taken = {(k, name): [100 * sized[name][0] * food[k][0], 100 * sized[name][1] * food[k][-1]] for k in ('hyper-arid', 'arid') for name in sized}
    assert [round(x, 1) for x in taken['hyper-arid', 'albatross_scale']] == [0.2, 0.6] and [round(x, 1) for x in taken['hyper-arid', 'hornero']] == [3.5, 10.3]
    assert [round(x, 2) for x in taken['arid', 'albatross_scale']] == [0.06, 0.28] and [round(x, 1) for x in taken['arid', 'hornero']] == [1.1, 4.9]
    same = [hauls[k]['by_size']['crow']['at_the_same_trips']['dry_food_mt_yr'] for k in ('9.4 km', '31 km')]
    assert 1.02 < same[1][0] / same[0][0] < same[1][1] / same[0][1] < 1.07
    great = [hauls[k]['by_size']['albatross_scale']['at_the_same_trips']['dry_food_mt_yr'] for k in ('9.4 km', '31 km')]
    assert 1.19 < great[1][0] / great[0][0] < great[1][1] / great[0][1] < 1.5
    # rock stays in its class, and in the driest land it settles on dry hollow floors
    assert all(0.93 <= round(by[k]['of_the_rock_settling_here_shed_in_this_class'], 2) <= 0.99 for k in by)
    driest = [by[k]['rock_settling_here_gt_yr'] for k in ('hyper-arid', 'arid')]
    assert round(sum(r['dry_hollow'] for r in driest) / sum(sum(r.values()) for r in driest), 2) == 0.94
    # the sand
    assert [by[k]['winds_sand_km2'] for k in ('hyper-arid', 'arid', 'semi-arid')] == [61954, 184687, 378] and sand['winds_sand_km2'] == 247018
    dry_ends = sum(by[k]['lasting_sand_at_ends_km3_yr'][kind] for k in ('hyper-arid', 'arid') for kind in ('lake deltas', 'closed hollows'))
    assert round(dry_ends, 3) == 0.037 and [round(x, 4) for x in sand['winds_take_km3_yr']] == [0.0017, 0.0027]
    assert round(100 * sand['made_seas_rock_over_the_return'], 1) == 4.9
    share = {name: (sea['at_ends_within_the_radius_km3_yr']['hyper-arid'] + sea['at_ends_within_the_radius_km3_yr']['arid'])
             / sum(sea['at_ends_within_the_radius_km3_yr'].values()) for name, sea in sand['made_seas'].items()}
    assert round(share['plain'], 2) == 0.11 and round(share['walled'], 2) == 0.51
    # the hauls
    assert [round(hauls[k]['hours_a_year']) for k in hauls] == [681, 796, 1481, 4037]
    assert [round(100 * hauls[k]['share_of_the_year']) for k in hauls] == [8, 9, 17, 46]
    crow = [hauls[k]['by_size']['crow']['in_the_same_hours']['of_the_lands_plant_growth'] for k in hauls]
    assert [[round(100 * x, 1) for x in pair] for pair in crow] == [[0.4, 0.5], [0.4, 0.6], [0.8, 1.1], [2.2, 3.0]]
    assert [round(x) for x in hauls['31 km']['kj_per_kg_delivered']] == [41, 81]
    # the regions: the run with dry hollows stands in semi-arid country, and Korolev's in wet
    south = regions['southern_highlands']
    assert round(south['classes']['semi-arid'], 2) == 0.56 and round(south['classes']['dry sub-humid'], 2) == 0.31 and round(south['classes']['arid'], 2) == 0.01
    assert round(south['plant_cover'][1], 2) == 0.84 and regions['korolev']['classes']['wet'] == 1.0


def test_what_the_write_up_says_of_a_lapse():
    regions = json.loads((HERE / 'results' / 'evolution.json').read_text())['regions']
    lake = regions['korolev']['wet']['none']['years']
    dry = regions['southern_highlands']['dry']['none']['years']
    years = ('1e+06', '3e+06', '1e+07', '3e+07')
    assert [round(lake[y]['hollow_volume_over_start'], 2) for y in years] == [0.5, 0.12, 0.0, 0.0]
    assert [round(lake[y]['hollow_area_over_start'], 2) for y in years] == [0.77, 0.53, 0.01, 0.0]
    assert [round(lake[y]['forms_kept'], 2) for y in years] == [0.96, 0.79, 0.46, 0.13]
    assert [round(dry[y]['hollow_volume_over_start'], 2) for y in years] == [0.94, 0.83, 0.44, 0.0]
    assert [round(dry[y]['forms_kept'], 2) for y in years] == [0.96, 0.85, 0.46, 0.09]
    assert [round(dry[y]['relief_over_start'], 2) for y in years] == [0.88, 0.78, 0.53, 0.36]
