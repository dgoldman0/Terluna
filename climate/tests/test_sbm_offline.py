"""The offline test of the simplified Betts-Miller scheme, its pure pieces: a CM1 column laid on PlaSim's levels keeps
each layer's vapour and heat, the three mappings agree on a smooth column, hour angles fall in their quarters of the
lunar day, the scheme's rain is the Python copy's, and the pass rule takes CM1's rain within half a millimetre a day
or half of it."""
import numpy as np
import pytest

from climate.gcm import sbm
from climate.gcm import sbm_offline as so

SIGMA = np.array([0.038, 0.119, 0.211, 0.317, 0.437, 0.567, 0.699, 0.823, 0.924, 0.983])
EDGES = np.array([0.0, 0.076, 0.162, 0.26, 0.374, 0.5, 0.634, 0.764, 0.882, 0.966, 1.0])
G, R = 1.6242, 287.0


def smooth_column(ps=121500.0):
    """A smooth column on the boxes' levels (cm1_run.vertical_grid, 111 levels to 150 km), the ground up: temperature
    falling 1.6 K/km to 160 K, vapour falling with a 6-km scale, the pressure in hydrostatic balance at lunar
    gravity."""
    from climate.crm.cm1_run import vertical_grid
    zw = np.asarray(vertical_grid(dz0=100.0, z_fine=1000.0, stretch=1.08, dz_max=2000.0, ztop=150000.0), float)
    zh, top_m = 0.5 * (zw[1:] + zw[:-1]), zw[-1]
    temperature = lambda z: np.maximum(300.0 - 1.6e-3 * z, 160.0)
    fine = np.linspace(0.0, top_m, 20001)
    lnp = np.log(ps) - np.concatenate([[0.0], np.cumsum(np.diff(fine) * G / (R * temperature(fine[1:])))])
    p = np.exp(np.interp(zh, fine, lnp))
    return dict(t=temperature(zh), q=0.02 * np.exp(-zh / 6000.0), p=p, ps=ps), zh, zw


def test_layer_means_keep_each_layers_integral_and_a_uniform_profile():
    col, zh, zw = smooth_column()
    pw = so.w_level_pressures(col['p'], col['ps'], zh, zw)
    assert pw[0] == col['ps'] and np.all(np.diff(pw) < 0.0)
    target = EDGES * col['ps']
    q = so.layer_means(col['q'], pw, target)
    covered = np.diff(target)
    covered[0] = target[1] - pw[-1]                                  # the top layer reaches above CM1's top
    assert (q * covered).sum() == pytest.approx((col['q'] * -np.diff(pw)).sum(), rel=1e-12)
    assert np.allclose(so.layer_means(np.full(zh.size, 7.0), pw, target), 7.0)


def test_the_three_mappings_agree_on_a_smooth_column():
    col, zh, zw = smooth_column()
    t_level, q_level, s_level, d_level = so.on_levels(col, 'level', SIGMA, EDGES, zh, zw)
    t_layer, q_layer, _, _ = so.on_levels(col, 'layer', SIGMA, EDGES, zh, zw)
    t_cm1, q_cm1, s_cm1, d_cm1 = so.on_levels(col, 'cm1', SIGMA, EDGES, zh, zw)
    assert np.array_equal(s_level, SIGMA) and np.allclose(d_level, np.diff(EDGES))
    assert np.all(np.abs(t_level[2:] - t_layer[2:]) < 0.5)            # above, CM1's top and the 160 K floor
    assert np.all(q_layer[1:] >= q_level[1:])                          # a layer's mean of falling vapour lies above
    assert t_cm1.size == zh.size and np.all(np.diff(s_cm1) > 0.0)       # its midpoint's value; top first
    assert d_cm1.sum() == pytest.approx(1.0 - so.w_level_pressures(col['p'], col['ps'], zh, zw)[-1] / col['ps'])
    with pytest.raises(ValueError):
        so.on_levels(col, 'spline', SIGMA, EDGES, zh, zw)


def test_hour_angles_fall_in_their_quarters():
    assert so.quarter_of([-90.0, -0.1, 0.0, 89.9, 90.0, 179.9, 180.0, -180.0, -90.1]).tolist() == [0, 0, 1, 1, 2, 2,
                                                                                                    3, 3, 3]
    q = so.quarter_of([-45.0, -45.0, 45.0, 135.0, -135.0])
    assert so.by_quarter([1.0, 3.0, 5.0, 7.0, 9.0], q) == [2.0, 5.0, 7.0, 9.0]
    assert so.by_quarter([1.0, 3.0, 5.0, 7.0, 9.0], q, [3.0, 1.0, 1.0, 1.0, 1.0]) == [1.5, 5.0, 7.0, 9.0]


def test_the_scheme_rains_what_the_copy_removes():
    t = np.array([165.97, 194.74, 221.45, 243.81, 260.84, 272.61, 281.35, 289.48, 295.52, 299.16])
    q = 1.15 * np.array([1e-8, 2.7e-6, 6.02e-5, 4.808e-4, 1.3599e-3, 2.9967e-3, 4.7614e-3, 6.7683e-3, 1.13541e-2,
                         1.32844e-2])
    c = dict(sbm.PLASIM, ga=G)
    got = so.scheme(t, q, SIGMA, np.diff(EDGES), 121490.0, 0.5, 0.0, 0.0, c)
    out = sbm.column(t, q, SIGMA, np.diff(EDGES), 121490.0, tau=so.TAU_S, rh=0.5, step=so.STEP_S, akap=so.KAPPA, c=c)
    assert got['rain'] == pytest.approx(out['rain'].sum() * so.DAY) and got['rain'] > 0.0
    assert got['deep'] and got['buoyant'] and got['cin'] == pytest.approx(out['cin'])


def test_the_pass_rule_takes_cm1s_rain_within_half_a_millimetre_or_half():
    cm1 = dict(sea=[0.0, 1.2, 2.0, 0.2], land=[2.9, 8.6, 1.5, 0.1])
    assert so.passes(cm1, cm1)['passes']
    near = dict(sea=[0.4, 1.7, 2.9, 0.6], land=[50.0, 12.8, 1.9, 9.0])  # the land's morning and night are not tested
    assert so.passes(near, cm1)['passes']
    far = dict(sea=[0.0, 1.2, 2.0, 0.8], land=[2.9, 8.6, 1.5, 0.1])    # the sea's night 0.6 above CM1's 0.2
    result = so.passes(far, cm1)
    assert not result['passes'] and [r['near'] for r in result['quarters']].count(False) == 1
