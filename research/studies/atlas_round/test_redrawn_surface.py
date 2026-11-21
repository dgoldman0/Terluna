"""Checks of the atlas round's redrawn surface: the creep's change to the load Love numbers, the water over each cell,
the rims, and the product once it has been computed."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import redrawn_surface as R

TABLE = dict(times_years=np.array([0.0, 1.0, 10.0, 100.0, 1000.0]), degrees=np.array([2, 3]),
             central__h_load=np.array([[-0.03, -0.05, -0.06, -0.08, -0.12], [-0.02, -0.03, -0.04, -0.05, -0.07]]),
             central__k_load=np.array([[-0.01, -0.02, -0.03, -0.04, -0.06], [-0.01, -0.01, -0.02, -0.02, -0.03]]))


def test_the_creep_adds_nothing_when_the_water_arrives_and_follows_the_table_after():
    dh, dk = R.creep_increment(TABLE, 'central', 0.0)
    assert not dh.any() and not dk.any()
    dh, dk = R.creep_increment(TABLE, 'central', 100.0)
    assert np.allclose(dh, [-0.05, -0.03]) and np.allclose(dk, [-0.03, -0.01])
    dh, _ = R.creep_increment(TABLE, 'central', 10 ** 1.5)          # halfway in log time between 10 and 100 years
    assert np.allclose(dh, [-0.04, -0.025])
    with pytest.raises(ValueError):
        R.creep_increment(TABLE, 'central', 1.0e4)


def test_the_creep_changes_only_the_degrees_the_table_holds():
    elastic = (np.full(6, -0.06), np.full(6, -0.02))
    h, k = R.love_numbers(elastic, TABLE, 'central', 1000.0, 5)
    assert np.allclose(h[[0, 1, 4, 5]], -0.06) and np.allclose(k[[0, 1, 4, 5]], -0.02)
    assert np.isclose(h[2], -0.15) and np.isclose(h[3], -0.11) and np.isclose(k[2], -0.07)
    h2, _ = R.love_numbers(elastic, TABLE, 'central', 1000.0, 2)
    assert h2.size == 3 and np.isclose(h2[2], -0.15)


def test_an_upper_sea_stands_at_its_rim_over_its_own_land_and_the_rest_at_the_common_level():
    eff = np.array([[-10.0, -2.0, 5.0, 20.0], [-1.0, 3.0, 12.0, 40.0]])
    region = np.array([[False, False, True, True], [False, False, True, False]])
    depth = R.water_depth(eff, 0.0, [(region, 25.0)])
    assert np.allclose(depth, [[10.0, 2.0, 20.0, 5.0], [1.0, 0.0, 13.0, 0.0]])


def test_a_rim_is_the_lowest_pass_to_any_other_sea_and_names_every_sea_joined_there():
    joined = {'a': np.array([0.0, 7.0, 5.0, 9.0, 5.4])}
    seeds = {'a': 0, 'b': 1, 'c': 2, 'd': 3, 'e': 4}
    assert R.rims_from(joined, seeds, ['a', 'b', 'c', 'd', 'e'], ['a']) == {'a': (5.0, ['c', 'e'])}


needs_product = pytest.mark.skipif(not R.OUT.exists(), reason='the redrawn surface has not been computed')


@needs_product
def test_the_atlas_case_reproduces_the_drainage_networks_rims_and_geographys_sills():
    p = json.loads(R.OUT.read_text())
    atlas = p['cases']['atlas']
    net = json.loads((R.ROOT / 'research/studies/drainage_network/results/network.json').read_text())
    seas = next(c for c in net['cases'] if c['rain'] == 1.0)['seas']
    for k in R.UPPER:
        rim = atlas['seas'][k]['surface_above_the_common_level_m'] + R.NW.BELOW_RIM_M
        assert abs(rim - seas[k]['rim_m']) < 0.1
        assert seas[k]['overflow_to'] in atlas['seas'][k]['joins_at_its_rim']
    loading = json.loads((R.ROOT / 'geography/results/crustal_loading.json').read_text())
    for name, sill in loading['base']['sills'].items():
        assert abs(atlas['sills'][name] - sill['sill_depth_m']) < 0.01
    assert atlas['level_change_m'] == 0.0


@needs_product
def test_every_case_holds_28_percent_at_its_common_level_and_the_design_settled():
    p = json.loads(R.OUT.read_text())
    for name, case in p['cases'].items():
        assert abs(case['common_level_share'] - R.SHARE) < 2e-4, name
        assert case['water_share'] > case['common_level_share'], name
        if case['settled'] and case['rounds']:
            assert case['rounds'][-1]['rims_moved_m'] < R.RIM_TOLERANCE_M, name
        if not case['settled']:
            assert all(lo <= hi for lo, hi in case['unsettled_rims_m'].values()), name
    assert p['cases'][R.DESIGN]['settled'] and p['cases']['elastic']['settled']
    assert p['cases']['atlas']['rounds'] == [] and p['cases']['atlas']['settled']


@needs_product
def test_the_sag_and_the_waters_gravity_deepen_the_south_sea_and_its_sill():
    p = json.loads(R.OUT.read_text())
    sill = 'Tenebrarum sub-basin to the South Pole-Aitken Sea'
    atlas = p['cases']['atlas']
    for name in ('elastic', R.DESIGN):
        case = p['cases'][name]
        assert case['sills'][sill] > atlas['sills'][sill]
        assert case['seas']['south_pole_aitken']['greatest_depth_m'] > atlas['seas']['south_pole_aitken']['greatest_depth_m']
        assert case['water_beyond_the_atlas_km3'] > 0
