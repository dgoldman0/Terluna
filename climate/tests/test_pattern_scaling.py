"""Pattern scaling's responses per kelvin, the rule that leaves a disputed response unapplied, and the scaled fields."""
import math

import numpy as np

from climate.gcm import pattern_scaling as ps

LAT = np.array([45.0, 0.0, -45.0])
LON = np.array([0.0, 120.0, 240.0])
LSM = np.array([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 0.0]])


def window(ts, rain_land, rain_water, cover=0.5):
    land = LSM > 0.5
    return dict(lat=LAT, lon=LON, lsm=LSM, ts_k=np.full((3, 3), float(ts)),
                pr_mm_day=np.where(land, rain_land, rain_water), clt_mean=np.full((3, 3), cover),
                cl_layer_mean=np.array([0.1, 0.2]))


def test_a_pair_gives_a_fraction_per_kelvin_over_land_and_water():
    warm, cool = window(302.0, 2.2, 3.3), window(300.0, 2.0, 3.0)
    r = ps.pair_response(warm, cool, {'pr_mm_day': 'ratio', 'ts_k': 'add', 'clt_mean': 'cover', 'cl_layer_mean': 'cover'})
    assert math.isclose(r['delta_K'], 2.0)
    assert math.isclose(r['fields']['pr_mm_day']['land']['per_K'], 0.05)      # 10% over 2 K
    assert math.isclose(r['fields']['pr_mm_day']['water']['per_K'], 0.05)
    assert math.isclose(r['fields']['ts_k']['land']['per_K'], 1.0)
    assert set(r['fields']['cl_layer_mean']) == {'all'}                       # a profile has one response


def test_pairs_that_disagree_in_sign_leave_the_field_unchanged():
    a = {'delta_K': 1.0, 'fields': {'x': {'land': {'per_K': 0.02}, 'water': {'per_K': -0.01}}}}
    b = {'delta_K': 1.0, 'fields': {'x': {'land': {'per_K': 0.04}, 'water': {'per_K': 0.01}}}}
    out = ps.combine([a, b])
    assert math.isclose(out['x']['land']['per_K'], 0.03) and out['x']['land']['resolved']
    assert out['x']['water']['per_K'] == 0.0 and not out['x']['water']['resolved']


def test_the_surface_temperature_moves_the_global_mean_by_the_step_exactly():
    base = window(301.0, 2.0, 3.0)
    applied = {'ts_k': {'land': {'per_K': 1.2, 'resolved': True}, 'water': {'per_K': 0.9, 'resolved': True}}}
    factor = ps.normalise_temperature(applied, base)
    moved = ps.apply(base, {'ts_k': 'add'}, applied, -0.75)
    assert math.isclose(ps.global_mean(moved['ts_k'], LAT, LON) - ps.global_mean(base['ts_k'], LAT, LON), -0.75)
    assert math.isclose(applied['ts_k']['land']['per_K'] / applied['ts_k']['water']['per_K'], 1.2 / 0.9)
    assert factor > 0


def test_apply_scales_amounts_and_keeps_covers_within_bounds():
    base = window(301.0, 2.0, 3.0, cover=0.99)
    applied = {'pr_mm_day': {'land': {'per_K': 0.02}, 'water': {'per_K': 0.01}},
               'clt_mean': {'land': {'per_K': -0.05}, 'water': {'per_K': -0.05}}}
    moved = ps.apply(base, {'pr_mm_day': 'ratio', 'clt_mean': 'cover'}, applied, -1.0)
    land = LSM > 0.5
    assert np.allclose(moved['pr_mm_day'][land], 2.0 * 0.98) and np.allclose(moved['pr_mm_day'][~land], 3.0 * 0.99)
    assert np.all(moved['clt_mean'] <= 1.0) and np.allclose(moved['clt_mean'], 1.0)
    assert moved['lsm'] is base['lsm']                                         # the mask is carried unchanged
