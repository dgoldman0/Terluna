"""PlaSim's vertical diffusion as pbl_mixing recomputes it, and the module's depth rule and fluxes
(climate/gcm/pbl_mixing.py)."""
import math

import numpy as np

from climate.gcm import levels as lv
from climate.gcm import pbl_mixing as pm

G, R, RV, CP = 1.6242, 287.05, 461.51, 1004.64
AKAP = R / CP


def column(theta_top_first, speed_top_first=None):
    half = lv.half_levels(lv.DEFAULT)
    sigma, edges = lv.full_levels(half), np.r_[0.0, half]
    t = np.asarray(theta_top_first, float) * sigma ** AKAP
    speed = np.zeros(sigma.size) if speed_top_first is None else np.asarray(speed_top_first, float)
    return t, np.zeros(sigma.size), speed, sigma, edges


def coefficients(t, q, speed, sigma, edges, lamm=pm.LAMBDA_M):
    return pm.vdiff_coefficients(t, q, speed, sigma, edges, G, R, R / RV, AKAP, lamm)


def test_in_neutral_air_the_diffusivity_grows_as_the_square_of_the_mixing_length():
    t, q, speed, sigma, edges = column(np.full(10, 290.0), np.linspace(12.0, 3.0, 10))
    short, long = coefficients(t, q, speed, sigma, edges), coefficients(t, q, speed, sigma, edges, 6.0 * pm.LAMBDA_M)
    assert np.allclose(short['ri'], 0.0) and np.allclose(short['f'], 1.0)
    assert np.allclose(long['k_m2_s'] / short['k_m2_s'], (long['l_m'] / short['l_m']) ** 2)
    zz = short['z_vdiff_m']
    lamh = pm.LAMBDA_M * math.sqrt(7.5)
    assert np.allclose(short['l_m'], lamh * pm.KARMAN * zz / (lamh + pm.KARMAN * zz))
    assert np.allclose(zz, -R * pm.TSCAL * np.log(edges[1:-1]) / G)


def test_unstable_air_takes_vdiffs_own_formula_and_still_scales_with_the_mixing_length():
    t, q, speed, sigma, edges = column(320.0 - 2.0 * np.arange(10)[::-1])     # cooler aloft: unstable
    short, long = coefficients(t, q, speed, sigma, edges), coefficients(t, q, speed, sigma, edges, 6.0 * pm.LAMBDA_M)
    k = 4                                                      # the half level near 31 km, where l is near lambda
    sh = edges[k + 1]
    zz = -R * 250.0 / G * math.log(sh)
    dz = R * 250.0 / G * math.log(sigma[k + 1] / sigma[k])
    lamh = 160.0 * math.sqrt(7.5)
    mixh = lamh * 0.4 * zz / (lamh + 0.4 * zz)
    zkfac = math.sqrt(((((zz + dz) / zz) ** (1.0 / 3.0) - 1.0) / dz) ** 3 / zz)
    rdsig = 1.0 / (sigma[k + 1] - sigma[k])
    dthds = (t[k] / sigma[k] ** AKAP - t[k + 1] / sigma[k + 1] ** AKAP) * rdsig
    ri = sh ** AKAP * R * dthds / (sh * rdsig ** 2)                             # the 1 m/s floor on the shear
    f = 1.0 - 15.0 * ri / (1.0 + 75.0 * math.sqrt(-ri) * zkfac * mixh * mixh)
    assert ri < -50.0 and math.isclose(short['ri'][k], ri) and math.isclose(short['f'][k], f)
    rthl = (np.diff(edges)[k + 1] + np.diff(edges)[k]) / (t[k] * np.diff(edges)[k + 1] + t[k + 1] * np.diff(edges)[k])
    assert math.isclose(short['k_m2_s'][k], mixh ** 2 * f * rdsig * G * sh * rthl / R)
    # at lunar heights the length's own term in the denominator stays small, so K keeps most of the l-squared gain
    gain = long['k_m2_s'][k] / short['k_m2_s'][k]
    assert (long['l_m'][k] / short['l_m'][k]) ** 2 > 20.0 and gain > 10.0


def test_stable_air_suppresses_the_diffusion_whatever_the_mixing_length():
    t, q, speed, sigma, edges = column(290.0 + 3.0 * np.arange(10)[::-1])
    long = coefficients(t, q, speed, sigma, edges, 6.0 * pm.LAMBDA_M)
    assert np.all(long['ri'] > 1.0) and np.all(long['f'] < 0.05)


def test_the_diffusion_carries_heat_and_vapour_down_their_gradients():
    t, q, speed, sigma, edges = column(np.linspace(300.0, 290.0, 10))      # warmer aloft: stable, flux downward
    q = np.linspace(0.002, 0.012, 10)                                      # moister below: flux upward
    k_sigma = np.full(9, 1.0e-6)
    heat, vapour = pm.diffusion_fluxes(t, q, 1.2e5, sigma, edges, k_sigma, G, CP, AKAP)
    assert np.all(heat < 0.0) and np.all(vapour > 0.0)
    heat0, vapour0 = pm.diffusion_fluxes(*column(np.full(10, 295.0))[:1], np.full(10, 0.01), 1.2e5, sigma, edges,
                                         k_sigma, G, CP, AKAP)
    assert np.allclose(heat0, 0.0) and np.allclose(vapour0, 0.0)


def test_heights_of_an_isothermal_column_follow_the_scale_height():
    half = lv.half_levels(lv.DEFAULT)
    sigma, edges = lv.full_levels(half), np.r_[0.0, half]
    t, q = np.full(10, 280.0), np.zeros(10)
    zh = pm.half_level_heights(t, q, 1.2e5, edges, G, R, R / RV)
    assert np.allclose(zh, R * 280.0 / G * np.log(1.0 / edges[1:-1]))
    zf = pm.full_level_heights(t, q, 1.2e5, sigma, edges, G, R, R / RV)
    assert np.allclose(zf, R * 280.0 / G * np.log(1.0 / sigma))


def test_the_mixed_layer_ends_where_virtual_potential_temperature_has_risen_half_a_kelvin():
    z = np.array([0.0, 1000.0, 2000.0, 3000.0])
    assert math.isclose(pm.rise_height(z, [300.0, 300.1, 300.3, 301.3], 0.0), 2200.0)
    assert math.isclose(pm.rise_height(z, [300.0, 300.1, 300.3, 301.3], 500.0), 2250.0)
    assert pm.rise_height(z, [300.0, 300.1, 300.2, 300.3], 0.0) is None


def test_the_diagnostics_snapshots_fall_half_a_mean_after_its_midpoint():
    means = (np.arange(-160.0, 560.0, 36.0) + 180.0) % 360.0 - 180.0         # ten 3-day means a lunar day, wrapping
    snaps = pm.snapshot_hours(means)
    assert np.allclose((snaps - means + 180.0) % 360.0 - 180.0, 18.0)
    assert np.all(snaps >= -180.0) and np.all(snaps < 180.0)
    assert pm.quarter_name(snaps[np.argmin(abs(means - (-100.0)))]) == 'sunrise-noon'      # -100 moves to -82
