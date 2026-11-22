"""Checks of the rivers' oxygen study: its physics by hand, its products against the design routing and against
themselves, and its regeneration where the research drive holds the columns."""
import json
import os
from pathlib import Path

import numpy as np
import pytest

from research.studies.river_oxygen import run as R
from research.studies.river_oxygen import waters as W
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ROUTING = ROOT / 'research' / 'studies' / 'atlas_round' / 'results' / 'design_routing.json'


@pytest.fixture(scope='module')
def mouths():
    return json.loads(R.MOUTHS.read_text())


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


# --- the physics, by hand

def test_the_channel_follows_the_deltas_studys_mature_rivers():
    depth, width = R.channel(1000.0)
    assert depth == pytest.approx(0.153 * 1000.0 ** 0.535) and 5.5 < depth < 7.0
    assert width == pytest.approx(3.875 * 1000.0 ** 0.465)
    assert depth * width * R.VELOCITY_M_S == pytest.approx(1000.0, rel=0.02)     # the fits close on the flow
    assert R.channel(1000.0, 2.0)[0] == pytest.approx(2.0 * depth)


def test_a_river_keeps_less_of_its_sources_temperature_the_farther_and_shallower_it_runs():
    near, _ = R.kept_share(10.0, 20.0, 0.0, 24.0)
    far, travel = R.kept_share(300.0, 20.0, 0.0, 24.0)
    shallow, _ = R.kept_share(300.0, 2.0, 0.0, 24.0)
    through_lakes, _ = R.kept_share(300.0, 20.0, 1.0e5, 24.0)
    assert 0.0 < shallow < far < near < 1.0 and through_lakes < far
    assert travel == pytest.approx(300e3 / 1.69)
    tau = 4.18e6 * 20.0 / 24.0
    assert far == pytest.approx(np.exp(-travel / tau))
    assert through_lakes == pytest.approx(far / (1.0 + 24.0 * 1.0e5 / 4.18e6))


def test_the_fall_warms_the_water_by_its_work():
    assert R.fall_warming_k(1000.0) == pytest.approx(MOON_SURFACE_GRAVITY * 1000.0 / 4180.0)
    assert R.fall_warming_k(1000.0) == pytest.approx(0.389, abs=0.002)
    assert R.fall_warming_k(1000.0, STANDARD_GRAVITY) == pytest.approx(2.346, abs=0.002)     # Joule's waterfall


def test_the_river_draws_its_lakes_water_to_its_own_depth():
    profile = np.array([[30.0, 28.0, 26.0, 26.0]])
    top, bottom = np.array([0.0, 1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0, 4.0])
    assert R.layer_mean(profile, top, bottom, 1.0) == pytest.approx([30.0])
    assert R.layer_mean(profile, top, bottom, 0.5) == pytest.approx([30.0])
    assert R.layer_mean(profile, top, bottom, 2.5) == pytest.approx([(30.0 + 28.0 + 0.5 * 26.0) / 2.5])
    assert R.layer_mean(profile, top, bottom, 99.0) == pytest.approx([27.5])


def test_the_rivers_temperature_is_the_lakes_water_carried_down():
    lake = np.array([25.0, 26.0, 27.0, 28.0])
    kept_all = R.river_temperature(lake, 1, 1.0, 1000.0, 1.0e-3)
    assert kept_all == pytest.approx(np.array([25.0, 25.0, 26.0, 27.0]) + R.fall_warming_k(1000.0))
    kept_none = R.river_temperature(lake, 0, 0.0, 1000.0, 1.0e-3)
    assert kept_none == pytest.approx(lake + 1.0 + R.fall_warming_k(1000.0))


def test_water_goes_deeper_the_colder_it_is_and_mixing_only_lifts_it():
    top, bottom = np.arange(0.0, 40.0, 1.0), np.arange(1.0, 41.0, 1.0)
    one = np.where(top < 10.0, 30.0, np.where(top < 20.0, 29.0, 28.0))            # a warm layer, a step, the nights' water
    sea = np.tile(one, (5, 1))
    river = np.array([31.0, 29.5, 28.5, 27.9, 20.0])
    cls, mixed, depth, surface = R.classes(river, sea, top, bottom, 27.5, 0.0)
    assert list(cls) == [0, 1, 1, 2, 3] and mixed == pytest.approx(river) and surface == pytest.approx(np.full(5, 30.0))
    assert depth[:3] == pytest.approx([0.0, 10.0, 20.0]) and np.isnan(depth[3:]).all()
    for entrained in (0.5, 2.0):
        lifted, _, lifted_depth, _ = R.classes(river, sea, top, bottom, 27.5, entrained)
        assert (lifted <= cls).all() and lifted[0] == 0
        assert (np.nan_to_num(lifted_depth, nan=99.0) <= np.nan_to_num(depth, nan=99.0)).all()
    half, mixed, _, _ = R.classes(river, sea, top, bottom, 27.5, 0.5)
    assert mixed[-1] == pytest.approx((20.0 + 0.5 * 30.0) / 1.5)
    assert list(half) == [0, 1, 1, 1, 3]


def test_the_plunge_lies_deeper_under_the_moons_gravity():
    moon = R.plunge_depth(1000.0, 100.0, 27.0, 29.0)
    earth = R.plunge_depth(1000.0, 100.0, 27.0, 29.0, STANDARD_GRAVITY)
    assert moon / earth == pytest.approx((STANDARD_GRAVITY / MOON_SURFACE_GRAVITY) ** (1.0 / 3.0))
    excess = (R.L.density(27.0) - R.L.density(29.0)) / R.L.density(29.0)
    assert moon == pytest.approx((10.0 ** 2 / (MOON_SURFACE_GRAVITY * excess * 0.25)) ** (1.0 / 3.0))


def test_the_demand_is_the_sea_oxygen_studys():
    setting = json.loads(R.SEA_OXYGEN.read_text())['setting']
    demand = R.demand_g_m2_yr(setting)
    assert demand['oligotrophic'] == pytest.approx(0.32 * 365.25)
    assert demand['eutrophic'] == pytest.approx(0.90 * 365.25)                    # Mueller et al.'s mean
    assert demand['oligotrophic'] < demand['mesotrophic'] < demand['eutrophic'] < demand['new_reservoir']


def test_a_column_looks_for_a_cell_of_its_own_kind():
    lat, lon = np.array([10.0, 0.0, -10.0]), np.array([0.0, 10.0, 20.0])
    water = np.array([[True, True, True], [False, False, True], [False, False, False]])
    assert W.cell_for(0.0, 9.0, lat, lon, water, False) == (1, 1, True)
    assert W.cell_for(0.0, 15.0, lat, lon, water, True) == (1, 2, True)            # water 5 degrees off, past nearer land
    assert W.cell_for(0.0, 9.0, lat, lon, water, True) == (1, 1, False)            # none within reach: the nearest cell
    assert W.cell_for(-10.0, 0.0, lat, lon, water, True) == (2, 0, False)


# --- the products

def test_the_seas_receive_what_the_design_routing_sends_them(mouths):
    summary = json.loads(ROUTING.read_text())['summary']['seas']
    assert set(mouths['seas']) == set(summary)
    for name, sea in mouths['seas'].items():
        assert sea['rivers_m3s'] == pytest.approx(summary[name]['receives_m3s'], abs=1.0)
        assert sea['area_km2'] == pytest.approx(summary[name]['area_km2'], rel=1e-4)
        assert sea['volume_km3'] == pytest.approx(sea['area_km2'] * sea['mean_depth_m'] / 1e3, rel=1e-6)
        areas = [b['area_km2'] for b in sea['below']]
        assert areas == sorted(areas, reverse=True) and areas[0] <= sea['area_km2'] * (1 + 1e-9)
    assert mouths['seas']['south_pole_aitken']['salt'] and 'south_pole_aitken' not in mouths['rivers']


def test_the_rivers_kept_carry_most_of_each_seas_water_and_each_leaves_a_lake(mouths):
    for name, group in mouths['rivers'].items():
        total = mouths['seas'][name]['rivers_m3s']
        assert group['kept_m3s'] + group['small_m3s'] == pytest.approx(total, rel=1e-9)
        assert group['kept_m3s'] >= 0.9 * total and 1 <= len(group['kept']) <= 20
        flows = [m['flow_m3s'] for m in group['kept']]
        assert flows == sorted(flows, reverse=True)
        for m in group['kept']:
            src = m['course']['source']
            assert src['kind'] == 'lake' and src['area_km2'] * 1e6 >= 1.7e5 * src['flow_through_m3s'] * (1 - 1e-6)
            assert m['course']['km'] < 500.0 and src['level_m'] > mouths['seas'][name]['level_m'] - 1.0
            assert m['floor']['lowest_depth_m'] >= m['floor']['first_lowest_depth_m'] - 1e-6


def test_the_shares_add_up_and_the_supply_stays_under_its_bound(product):
    classes = product['classes']
    for label, case in product['cases'].items():
        for name, sea in case.items():
            if name in ('settings', 'all_fresh_seas'):
                continue
            assert sum(sea['shares_of_river_water'][c] for c in classes) == pytest.approx(1.0)
            assert sea['oxygen_to_the_floor_t_yr'] <= sea['oxygen_below_the_nights_mixing_t_yr'] + 1e-6
            for state, met in sea['share_of_demand_met'].items():
                bound = sea['if_every_river_sank']
                limit = bound['share_of_demand_met'][state] * (1.0 + case['settings']['entrained']) * 1.1
                assert 0.0 <= met <= limit
                assert bound['with_two_volumes_of_surface_water'][state] == pytest.approx(3.0 * bound['share_of_demand_met'][state])
    central = product['cases']['central']
    none, two = product['cases']['no surface water taken in'], product['cases']['two volumes of surface water taken in']
    for name in ('near_side_main', 'moscoviense'):
        deep = lambda case: sum(case[name]['shares_of_river_water'][c] for c in classes[2:])      # noqa: E731
        assert deep(none) >= deep(central) - 1e-12 >= deep(two) - 2e-12
        span = product['range_over_the_cases'][name]['share_of_river_water_below_the_nights_mixing']
        assert span['least'] <= deep(central) <= span['most']
        for river in central[name]['rivers']:
            assert sum(river['shares'].values()) == pytest.approx(1.0)
            assert river['sea_foot_c']['coldest'] <= river['sea_foot_c']['mean'] <= river['sea_surface_layer_c']['highest']


def test_the_products_bind_their_inputs(product):
    for rel, digest in product['inputs'].items():
        assert R.sha256(ROOT / rel) == digest, rel
    waters = json.loads(R.WATERS.read_text())
    assert waters['columns_not_yet_run'] == []
    assert waters['inputs'][str(R.MOUTHS.relative_to(ROOT))] == R.sha256(R.MOUTHS)


def test_the_product_regenerates_from_the_columns_on_the_drive(product):
    waters = json.loads(R.WATERS.read_text())
    if not all((W.COLUMNS / spec['file']).exists() for spec in waters['columns'].values()):
        pytest.skip('the water columns are on the research drive (river-oxygen/columns/)')
    again = json.loads(json.dumps(R.build()['cases']))
    assert again == product['cases']


@pytest.mark.skipif(not os.environ.get('TERLUNA_RIVER_OXYGEN_REGENERATE'),
                    reason='reads the design routing on the research drive; set TERLUNA_RIVER_OXYGEN_REGENERATE=1')
def test_the_mouths_regenerate(mouths):
    from research.studies.river_oxygen import mouths as M
    again = M.build(log=lambda *a: None)
    assert json.loads(json.dumps(again['seas'])) == mouths['seas']
    assert json.loads(json.dumps(again['rivers'])) == mouths['rivers']
