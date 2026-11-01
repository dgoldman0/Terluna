"""The GCM's rain evaporation: the pooling of runs, the fit to each CM1 box's share at its place and the product's two
measures."""
import json

import numpy as np
import pytest

from climate.gcm import rain_evaporation as re_


def test_runs_at_one_gamma_pool_their_sums():
    a = dict(run='a', years=[1, 2], gamma=0.01, sums=dict(global_=[2.0, 1.0, 1.0]), residual=[1.0, 1.0])
    b = dict(run='b', years=[0], gamma=0.01, sums=dict(global_=[2.0, 1.6, 1.0]), residual=[0.4, 1.0])
    c = dict(run='c', years=[0], gamma=0.02, sums=dict(global_=[2.0, 0.8, 1.0]), residual=[1.2, 1.0])
    for m in (a, b, c):
        m['sums'] = {'global': m['sums'].pop('global_')}
    pooled = re_.pool([c, a, b])
    assert [p['gamma'] for p in pooled] == [0.01, 0.02] and pooled[0]['years'] == 3
    assert pooled[0]['shares']['global']['evaporated_share'] == pytest.approx(1.0 - 2.6 / 4.0)
    assert pooled[0]['shares']['global']['evaporated_share_by_residual'] == pytest.approx((1.4 / 2.0) / (4.0 / 2.0))


def test_the_fit_interpolates_gamma_in_its_log():
    pooled = [dict(gamma=g, shares=dict(near_side_sea=dict(evaporated_share=s))) for g, s in ((0.01, 0.3), (0.04, 0.6))]
    assert re_.fit(pooled, 'near_side_sea', 0.45) == pytest.approx(0.02, rel=1e-9)
    assert re_.fit(pooled, 'near_side_sea', 0.9) is None                 # outside the runs given
    assert re_.fit(pooled, 'near_side_sea', None) is None                # a box that has not run
    flat = [dict(gamma=g, shares=dict(cm1_box=dict(evaporated_share=s))) for g, s in ((0.01, 0.4), (0.04, 0.4))]
    assert re_.fit(flat, 'cm1_box', 0.4) is None                         # no single gamma


def test_each_region_is_its_place_and_surface():
    lat, lon = np.array([5.0, 0.0, -5.0, -10.0]), np.arange(0.0, 360.0, 10.0)
    land = np.zeros((4, 36), bool)
    land[:, :2] = True
    height = np.zeros((4, 36))
    sea = re_.region_mask(lat, lon, land, height, re_.REGIONS['near_side_sea'])
    assert set(np.nonzero(sea)[1]) == {30, 31, 32, 33} and not (sea & land).any()       # 295-335 E
    assert re_.region_mask(lat, lon, land, height, re_.REGIONS["cm1_box"]).sum() == 8         # the land at 0 and 10 E


def test_the_product_sets_each_box_beside_its_place_and_its_two_measures_agree():
    if not re_.OUT.exists():
        pytest.skip('the rain evaporation has not been measured')
    out = json.loads(re_.OUT.read_text())
    assert out['schema'] == re_.SCHEMA
    small = out['cm1']['box_0e_small_lunar_rain_day29']
    assert small['share_of_falling_water_evaporated'] == pytest.approx(0.48, abs=0.005) and small['region'] == 'cm1_box'
    assert out['cm1_limits_control'] == pytest.approx(0.857, abs=0.001)
    gammas = [p['gamma'] for p in out['pooled']]
    assert gammas == sorted(gammas) and gammas[0] < small['gamma'] < gammas[-1]
    for p in out['pooled']:
        g = p['shares']['global']
        assert abs(g['evaporated_share'] - g['evaporated_share_by_residual']) < 0.02
    for case, row in out['cm1'].items():
        assert row['region'] in re_.REGIONS and (row['gamma'] is None) == (row['share_of_falling_water_evaporated'] is None
                                                                         or row['gamma'] is None)
