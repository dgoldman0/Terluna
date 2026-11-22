"""Tests of the living return of rock: the carriers' arithmetic, the loop, the rock module's parts and the products."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from research.studies.rock_return import carriers as CA
from research.studies.rock_return import rock as RK
from research.studies.rock_return import run as RN
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(scope='module')
def rock():
    return json.loads(RK.OUT.read_text())


@pytest.fixture(scope='module')
def product():
    return json.loads(RN.OUT.read_text())


# ---- the carriers' arithmetic ------------------------------------------------------------------------------------
def test_flight_costs_weight_over_glide_and_the_climb_in_full():
    assert math.isclose(CA.payload_cost_j_per_kg(10000.0, 0.0, 20.0, 0.25, 1.6), 1.6 * 10000.0 / 5.0)
    assert math.isclose(CA.payload_cost_j_per_kg(0.0, 1000.0, 20.0, 0.25, 1.6), 1.6 * 1000.0 / 0.25)
    moon = CA.payload_cost_j_per_kg(9000.0, 1000.0, 15.0, 0.2)
    earth = CA.payload_cost_j_per_kg(9000.0, 1000.0, 15.0, 0.2, STANDARD_GRAVITY)
    assert math.isclose(earth / moon, STANDARD_GRAVITY / MOON_SURFACE_GRAVITY)
    with pytest.raises(ValueError):
        CA.payload_cost_j_per_kg(1.0, 1.0, 0.0, 0.2)


def test_a_heavier_load_and_a_glide_home_both_cheapen_each_kilogram():
    alone = CA.payload_cost_j_per_kg(9000.0, 1000.0, 20.0, 0.25)
    light = CA.delivered_cost_j_per_kg(9000.0, 1000.0, 20.0, 0.25, 0.25)
    heavy = CA.delivered_cost_j_per_kg(9000.0, 1000.0, 20.0, 0.25, 1.0)
    assert math.isclose(light, alone * 5.0) and math.isclose(heavy, alone * 2.0)
    powered = CA.delivered_cost_j_per_kg(9000.0, 1000.0, 20.0, 0.25, 1.0, glides_home=0.0)
    assert powered > heavy and math.isclose(powered - heavy, MOON_SURFACE_GRAVITY * 9000.0 / 5.0)


def test_grains_weigh_less_in_water_and_pressure_is_a_sixth():
    dry = MOON_SURFACE_GRAVITY * 1000.0 / 0.25
    assert math.isclose(CA.lift_in_water_j_per_kg(1000.0, 0.25), dry * (1.0 - 1000.0 / 2600.0))
    deep = CA.depth_at_earths_pressure(200.0)
    assert math.isclose(deep, 200.0 * STANDARD_GRAVITY / MOON_SURFACE_GRAVITY)
    assert math.isclose(CA.pressure_atm(deep), CA.pressure_atm(200.0, gravity=STANDARD_GRAVITY))


def test_reach_shares_add_to_the_whole():
    r = CA.reach(0.4, 0.6, 0.1, 0.5, 0.25)
    assert math.isclose(r['flyers_and_waders'], 0.4 + 0.06)
    assert math.isclose(r['with_divers_at_the_delta_fronts'], 0.4 + 0.3)
    assert math.isclose(r['with_divers_on_the_floors_within_a_dive'] + r['left_on_the_deeper_floors'], 1.0)
    with pytest.raises(ValueError):
        CA.reach(0.4, 0.6, 0.6, 0.5, 0.2)


def test_numbers_of_carriers_and_their_hours():
    assert math.isclose(CA.carriers_for(3.0e12, 1000.0, 0.5, 1000), 6.0e6)
    assert math.isclose(CA.hours_a_year(1000, 9000.0, 15.0, 1200.0), 1000 * (1200.0 + 1200.0) / 3600.0)


def test_creep_sets_how_far_a_line_of_rock_can_feed():
    slope = CA.crest_slope(5.0e-5, 100.0, 0.01)
    assert math.isclose(slope, 0.5)
    assert math.isclose(CA.reach_for_slope(slope, 5.0e-5, 0.01), 100.0)


def test_urea_gives_one_calcite_for_two_nitrogens():
    assert math.isclose(CA.calcite_by_urea_kg(2.0 * 14.007), 100.087)


# ---- the loop ----------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('growth,death,half_years', [(0.05, 0.02, 100.0), (0.1, 0.03, 30.0), (0.3, 0.1, 10.0), (1.0, 0.01, 1.0)])
def test_the_loop_rests_where_the_return_equals_the_supply(growth, death, half_years):
    supply, duty = 3.0e12, 600.0
    half = half_years * supply
    s, b = CA.loop_equilibrium(supply, duty, growth, death, half)
    take = s / (s + half)
    assert math.isclose(duty * b * take, supply, rel_tol=1e-12)
    assert math.isclose(growth * take, death, rel_tol=1e-12)
    assert math.isclose(s / supply, half_years * death / (growth - death))
    eig = CA.loop_eigenvalues(supply, duty, growth, death, half)
    assert (eig.real < 0).all()


def test_the_resting_deposit_does_not_depend_on_the_supply_and_the_carriers_follow_it():
    one = CA.loop_equilibrium(1.0e12, 600.0, 0.1, 0.03, 3.0e13)
    two = CA.loop_equilibrium(2.0e12, 600.0, 0.1, 0.03, 3.0e13)
    assert math.isclose(one[0], two[0]) and math.isclose(two[1], 2.0 * one[1])
    with pytest.raises(ValueError):
        CA.loop_equilibrium(1.0, 1.0, 0.02, 0.05, 1.0)


def test_the_loop_recovers_from_losing_half_its_carriers():
    supply, duty, growth, death = 3.0e12, 600.0, 0.3, 0.1
    half = 10.0 * supply
    s, b = CA.loop_equilibrium(supply, duty, growth, death, half)
    t, dep, car, ret = CA.loop_run(supply, duty, growth, death, half, s, 0.5 * b, 400.0, dt=0.1)
    assert math.isclose(ret[0], 0.5 * supply, rel_tol=1e-9)
    assert math.isclose(ret[-1], supply, rel_tol=1e-3) and math.isclose(car[-1], b, rel_tol=1e-3) and math.isclose(dep[-1], s, rel_tol=1e-3)
    assert dep.max() > s
    assert 0.0 < CA.years_until_within(t, ret, supply) < 100.0
    held = CA.loop_run(supply, duty, growth, death, half, s, b, 50.0, dt=0.1)
    assert np.allclose(held[1], s, rtol=1e-9) and np.allclose(held[3], supply, rtol=1e-9)


def test_lakes_stay_open_longer_as_more_is_returned():
    fill = {'0.1': 1.0e4, '0.5': 1.0e6, '0.9': 1.0e8}
    assert math.isclose(CA.open_after(1.0e6, fill, 0.0), 0.5)
    assert math.isclose(CA.open_after(1.0e7, fill, 0.9), 0.5)
    assert CA.open_after(1.0e6, fill, 0.5) > CA.open_after(1.0e6, fill, 0.0)
    assert CA.open_after(1.0e2, fill, 0.0) == 1.0 and CA.open_after(1.0e6, fill, 1.0) == 1.0


# ---- the rock module's parts -------------------------------------------------------------------------------------
def test_band_weights_keep_each_bands_mean_lowering():
    relief = np.array([100.0, 300.0, 50.0, 150.0, 999.0])
    area = np.array([1.0, 1.0, 2.0, 2.0, 5.0])
    eroding = np.array([True, True, True, True, False])
    abs_lat = np.array([5.0, 5.0, 70.0, 70.0, 5.0])
    w = RK.band_weights(relief, area, eroding, abs_lat)
    assert math.isclose((w[:2] * area[:2]).sum() / area[:2].sum(), 1.0)
    assert math.isclose((w[2:4] * area[2:4]).sum() / area[2:4].sum(), 1.0)
    assert w[4] == 0.0 and math.isclose(w[1] / w[0], 3.0)


def test_rock_stops_at_the_first_place_it_can_settle():
    # 0 -> 1 -> 2 (a lake) -> 3 -> 4 (the sea); 5 -> 3
    receiver = np.array([1, 2, 3, 4, 4, 3])
    donor = np.array([True, True, True, True, False, True])
    from research.studies.sea_nutrients import loads as LO
    fronts = LO.kahn_fronts(receiver, donor)
    settles = np.array([False, False, True, False, True, False])
    first = RK.first_settling(receiver, fronts, settles)
    assert first.tolist() == [2, 2, 2, 4, 4, 4]


def test_distances_on_the_grid():
    cols = 5760
    one_cell = MOON_RADIUS * math.radians(1.0 / 16.0)
    a = np.array([1000 * cols + 10])
    assert math.isclose(RK.arc_m(a, a + cols, cols)[0], one_cell, rel_tol=1e-6)
    assert RK.arc_m(a, a, cols)[0] == 0.0
    equator = np.array([1440 * cols + 10])
    assert math.isclose(RK.arc_m(equator, equator + 1, cols)[0], one_cell * math.cos(math.radians(1.0 / 32.0)), rel_tol=1e-4)


def test_weighted_means_floors_and_fill_times():
    w = RK.weighted(np.array([1.0, 3.0]), np.array([1.0, 3.0]))
    assert math.isclose(w['mean'], 2.5)
    assert RK.weighted(np.array([1.0]), np.array([0.0]))['mean'] is None
    floor = RK.floor_by_depth(np.array([1.0, 10.0, 100.0, 1000.0]), np.array([1.0, 1.0, 1.0, 1.0]))
    assert floor == {'under 2 m': 0.25, 'under 20 m': 0.5, 'under 200 m': 0.75}
    volume, kept, area, ok = np.array([1.0e6, 1.0e9, 5.0]), np.array([1.0e6, 1.0e6, 0.0]), np.array([1.0, 1.0, 2.0]), np.array([True, True, True])
    still = RK.still_open(volume, kept, area, ok)['1000 kg/m3']
    assert math.isclose(still['1000'], 0.75) and math.isclose(still['1e+08'], 0.5)       # the lake that keeps nothing never fills
    q = RK.fill_quantiles(volume, kept, area, ok)
    assert math.isclose(q['0.2'], 1.0e3) and math.isclose(q['0.4'], 1.0e6) and q['0.9'] is None


# ---- the products ------------------------------------------------------------------------------------------------
def test_the_rock_product_is_whole(rock):
    assert rock['schema'] == RK.SCHEMA
    assert rock['producer']['code_sha256'] == sha256(HERE / 'rock.py')
    for rel, digest in rock['inputs'].items():
        assert sha256(ROOT / rel) == digest, rel
    back = rock['the_way_back']
    assert math.isclose(sum(back[p]['rock_share'] for p in RK.PLACES), 1.0, rel_tol=1e-9)
    for key, kk in rock['kinds'].items():
        assert math.isclose(sum(kk[n]['rock_share'] for n in ('wild', 'hybrid', 'peoples', 'open_beyond_need', 'not_classed')), 1.0, rel_tol=1e-9)
        assert kk['not_classed']['land_km2'] == 0.0
    for case in rock['carried_on'].values():
        assert math.isclose(case['kept_by_lakes_that_spill'] + case['ends_where_water_ends'] + case['reaches_the_seas'], 1.0, rel_tol=1e-6)
    closed = rock['closed_waters']
    # what ends on dry ground is what first settles there and the little that lakes upstream pass on to it
    assert 0.0 <= closed['rock_share_on_dry_ground'] - back['dry_hollow']['rock_share'] < 1e-3


def test_the_way_back_on_the_design_surface_is_the_terrain_studys(rock):
    """Two surfaces and two routings, one rule for the rock: the mean way back agrees within a tenth."""
    was = rock['setting']['terrain_study']
    assert abs(rock['the_way_back']['all']['lift_m']['mean'] / was['mean_lift_m'] - 1.0) < 0.1
    assert abs(rock['the_way_back']['all']['haul_m']['mean'] / 1e3 / was['mean_distance_km'] - 1.0) < 0.1
    assert abs(rock['the_way_back']['dry_hollow']['rock_share'] - was['where_the_rock_settles']['dry_hollows']) < 0.02
    assert 40.0 < rock['rock']['mean_lowering_m_per_myr'] < 70.0


def test_the_stored_sizing_is_what_the_code_gives(product):
    fresh = json.loads(json.dumps(RN.build()))

    def same(a, b, path=''):
        if isinstance(a, dict):
            assert a.keys() == b.keys(), path
            for k in a:
                same(a[k], b[k], f'{path}/{k}')
        elif isinstance(a, list):
            assert len(a) == len(b), path
            for i, (x, y) in enumerate(zip(a, b)):
                same(x, y, f'{path}[{i}]')
        elif isinstance(a, float) and isinstance(b, float):
            assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12), path
        else:
            assert a == b, path
    same(product, fresh)
    for rel, digest in product['inputs'].items():
        assert sha256(ROOT / rel) == digest, rel


def test_what_the_write_up_says_of_the_carriers(product):
    cost = product['carriers']['cost_kj_per_kg_delivered']
    assert max(cost['sediments_energy_over_the_flight']['earth']) < 1.0 < max(cost['sediments_energy_over_the_flight']['moon'])
    middle = product['carriers']['numbers']['1000 trips, load 0.6 of body']
    assert 2.5e6 < middle['flyers_for_all_the_land_central'] < 4.0e6 and middle['share_of_the_year'] < 0.1
    assert max(product['carriers']['food_over_plant_growth']) < 0.001
    sizes = product['carriers']['by_size_at_1000_trips_and_0p6']
    assert all(math.isclose(v['live_mass_mt'], middle['live_mass_mt'], rel_tol=1e-9) for v in sizes.values())
    assert sizes['crow']['metres_between_carriers'] < 100.0 < 2000.0 < sizes['albatross_scale']['metres_between_carriers']
    assert max(sizes['hornero']['food_over_plant_growth']) < 0.02


def test_what_the_write_up_says_of_reach_cement_and_the_loop(product):
    reach = product['reach']['all_the_land']['reach']
    assert 0.4 < reach['deltas of clay']['flyers_and_waders'] < reach['deltas of sand']['flyers_and_waders'] < 0.7
    assert reach['deltas of sand']['left_on_the_deeper_floors'] < 0.05 < 0.2 < reach['deltas of clay']['left_on_the_deeper_floors']
    assert product['reach']['moons_depth_at_the_pressure_of_earths_deepest_dives_m']['seal'] > 8000.0
    put = product['staying_put']
    assert max(put['that_calcite_over_the_rock_returned']) < 0.001 < 0.03 < put['lime_over_the_particles_by_mass']['fast']
    feeds = put['how_far_each_side_a_line_of_returned_rock_can_feed_m']
    assert 25.0 < feeds['creep 0.003 m2/yr']['central'] < feeds['creep 0.04 m2/yr']['central'] < 450.0
    for name in RN.LOOPS:
        loop = product['loop'][name]
        assert all(re < 0 for re, _ in loop['eigenvalues_per_yr'])
        assert loop['after_half_the_carriers_are_lost']['years_until_the_return_is_within_5_percent'] < 600.0
        assert math.isclose(loop['after_the_erosion_doubles']['carriers_end_at_times_the_start'], 2.0, rel_tol=1e-3)
    assert product['loop']['the_lakes_for_scale']['years_for_a_hundredth_of_their_area_to_fill'] > 1.0e4
    lakes = product['waters']['lake_area_still_open']
    assert lakes['nothing returned']['1e+08'] < 0.15 < 0.6 < lakes['deltas of sand, with_divers_at_the_delta_fronts']['1e+08']


def test_the_sediments_organic_carbon_is_the_carbon_cycle_studys():
    cc = pytest.importorskip('research.studies.carbon_cycle.run')
    assert RN.DELTA_ORGANIC_C == (cc.DELTA_OC[0], cc.DELTA_OC[-1])


@pytest.mark.skipif(os.environ.get('TERLUNA_ROCK_RETURN_REGENERATE') != '1',
                    reason='reads the research drive and takes about four minutes; set TERLUNA_ROCK_RETURN_REGENERATE=1')
def test_the_rock_product_regenerates(rock):
    fresh = RK.build(log=lambda *a: None)
    assert math.isclose(fresh['rock']['kg_yr_central'], rock['rock']['kg_yr_central'], rel_tol=1e-9)
    for place in RK.PLACES:
        assert math.isclose(fresh['the_way_back'][place]['rock_share'], rock['the_way_back'][place]['rock_share'], rel_tol=1e-9)
    assert fresh['lakes']['keep_half'] == rock['lakes']['keep_half']
    assert math.isclose(fresh['carried_on']['central']['reaches_the_seas'], rock['carried_on']['central']['reaches_the_seas'], rel_tol=1e-9)
