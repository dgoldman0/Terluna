"""Tests for the polar-light study: product binding, the refraction integral, the ephemeris's geometry, the photon
conversion, the two-state clouds and the terrain's horizons, each against an independent expectation."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from research.studies.polar_light import run
from shared.constants import MOON_ORBIT_INCLINATION_DEG, MOON_RADIUS
from shared.provenance import constants_changed

PRODUCT = run.OUT


@pytest.fixture(scope='module')
def product():
    if not PRODUCT.exists():
        pytest.skip('product not built')
    return json.loads(PRODUCT.read_text())


@pytest.fixture(scope='module')
def per():
    return run.periods()


# ---- binding ---------------------------------------------------------------------------------------------------------
def test_product_schema_and_files(product):
    assert product['schema'] == run.SCHEMA
    for path, digest in product['producer']['files'].items():
        assert run.sha(run.ROOT / path) == digest, path
    for path, digest in product['producer']['inputs'].items():
        assert run.sha(run.ROOT / path) == digest, path


def test_constants_current(product):
    assert constants_changed(product['producer']['constants']) == {}


def test_product_checks(product):
    c = product['checks']
    assert c['solution_is_atlas_source'] and c['atlas_is_solved_sky_atlas']
    d = c['direct_normal_lux_overhead']
    assert abs(d['channels'] / d['solved_sky'] - 1) < 0.005
    f = c['diffuse_lux_at_minus_10']
    assert abs(f['channels'] / f['calendar'] - 1) < 0.03
    assert abs(c['equinox_hours_half_lunar_day'] - 708.73 / 2) < 3.0
    assert 0.55 < c['refraction_earth_control_deg'] < 0.60          # Earth's measured horizontal refraction, about 34'
    p = c['plant_model_twilight_par_at_minus_25']
    assert abs(p['this_study'] / p['plant_model'] - 1) < 0.15


# ---- refraction ------------------------------------------------------------------------------------------------------
def test_refraction_exponential_atmosphere():
    """An isothermal exponential atmosphere: horizontal refraction close to (n0 - 1) sqrt(pi r / 2H)."""
    H, r0 = 50e3, MOON_RADIUS
    z = np.linspace(0.0, 900e3, 2000)
    n_s = run.STANDARD_AIR[0] / (run.BOLTZMANN * run.STANDARD_AIR[1])
    N = 1.2 * n_s * np.exp(-z / H)
    r, n0 = run.refraction(z, N, r0, [0.0, 10.0])
    expect = math.degrees(n0 * math.sqrt(math.pi * r0 / (2 * H)))
    assert abs(r[0] / expect - 1) < 0.03
    # above the horizon: the first-order straight-ray integral through the same spherical atmosphere,
    # (n0 - 1)/H int exp(-s/H) (b/r) / sqrt(1 - (b/r)^2) ds with b = r0 cos(h); the plane-parallel (n0 - 1) cot(h)
    # fails here, since H/r0 is 3%
    s = np.linspace(0.0, 900e3, 400001)
    b = r0 * math.cos(math.radians(10.0))
    q = b / (r0 + s)
    first = n0 / H * np.trapezoid(np.exp(-s / H) * q / np.sqrt(1 - q ** 2), s)
    assert abs(r[1] / math.degrees(first) - 1) < 0.02
    assert r[1] < 0.8 * math.degrees(n0 / math.tan(math.radians(10.0)))


def test_refraction_moon_against_earth():
    ref = run.refraction_summary()
    assert ref['moon']['refraction_deg'][0] < 0.25 * ref['earth']['refraction_deg'][0]
    assert 40.0 < ref['moon']['density_scale_height_km'] < 80.0
    assert np.all(np.diff(ref['moon']['refraction_deg']) < 0)


# ---- geometry --------------------------------------------------------------------------------------------------------
def test_pole_elevation_is_declination(per):
    s = run.lunar_days([90.0], 1.613, per)
    h = run.sun_elevation(90.0, 0.0, s['dec'], s['lam_ss'])
    assert np.allclose(h, s['dec'], atol=1e-9)
    assert s['dec'].max() <= 1.613 + 1e-9


def test_eclipse_year(per):
    assert abs(per['eclipse_year_days'] - 346.62) < 0.1
    assert per['run_tilt_deg'] == pytest.approx(1.613)


def test_equinox_day_is_half(per):
    s = run.lunar_days([0.0], 1.613, per)
    h = run.sun_elevation(75.0, 0.0, s['dec'], s['lam_ss'])
    hours = (h > 0).sum() * run.STEP_H
    assert abs(hours - per['synodic_hours'] / 2) < 3.0


def test_full_earth_stands_high_in_winter(per):
    """At full Earth the Earth stands near the antisolar point plus the orbit's tilt: latitude -(1 + i/eps) dec."""
    eps = 1.613
    s = run.lunar_days(run.SEASONS, eps, per, lon=0.0)
    k = np.argmin(np.abs(s['H']))                           # local noon on the sub-Earth meridian is new Earth;
    mid = np.argmin(np.abs(np.abs(s['H']) - 180.0))          # local midnight there is full Earth
    expect = -(1 + MOON_ORBIT_INCLINATION_DEG / eps) * s['dec'][:, mid]
    assert np.allclose(s['beta_e'][:, mid], expect, atol=0.15)
    assert k != mid


def test_distance_factor_time_mean(per):
    """The time mean of (1 AU / r)^2 over the orbit is 1 / sqrt(1 - e^2)."""
    u = np.arange(0.0, 360.0, 0.1)
    for psi in (0.0, 90.0, 270.0):
        w = run.season_weights(u, psi, dict(per, node_days=1e30))
        mean = np.sum(w * run.distance_factor(u, psi, per['eccentricity']))
        assert abs(mean - 1 / math.sqrt(1 - per['eccentricity'] ** 2)) < 1e-5


# ---- light -----------------------------------------------------------------------------------------------------------
def test_photon_ratios():
    ph = run.photon_ratios()
    assert ph['par_channels'] > 0
    q = np.interp(np.arange(-30.0, 90.1, 5.0), ph['diffuse_elevation_deg'], ph['diffuse_umol_per_lux'])
    assert np.all((q > 0.015) & (q < 0.02))
    assert np.interp(0.0, ph['direct_elevation_deg'], ph['direct_umol_per_lux']) > \
        np.interp(90.0, ph['direct_elevation_deg'], ph['direct_umol_per_lux'])      # the low Sun's beam is red


def test_two_state_cloud_limits():
    t = np.linspace(0, 2 * math.pi, 200, endpoint=False)
    series = (10 ** (2 * np.cos(t)))[:, None] * np.ones((1, 3))
    clear = run.hours_below(series, 1.0, 1.0)
    assert np.allclose(run.hours_below_cloud(series, 1.0, 1.0, np.zeros_like(series), 0.25), clear)
    over = run.hours_below(series * 0.25, 1.0, 1.0)
    assert np.allclose(run.hours_below_cloud(series, 1.0, 1.0, np.ones_like(series), 0.25), over)
    assert np.all(over > clear)


# ---- terrain ---------------------------------------------------------------------------------------------------------
@pytest.fixture(scope='module')
def small_grid(monkeypatch_module):
    monkeypatch_module.setattr(run, 'CAP_DEG', 0.2)
    reach = run.REACH_PX[-1] * run.DX_KM / (MOON_RADIUS / 1000.0) * 180.0 / math.pi + 0.5
    return run.stereo_grid('north', reach)


@pytest.fixture(scope='module')
def monkeypatch_module():
    mp = pytest.MonkeyPatch()
    yield mp
    mp.undo()


def test_horizon_flat_sphere_and_wall(small_grid):
    g = small_grid
    R = MOON_RADIUS / 1000.0
    surf = np.zeros(g['lat'].shape)
    mask = np.zeros(g['lat'].shape, bool)
    mask[g['n'], g['n']] = True                             # the pole
    hz, n0, north, east = run.horizons(surf, g, mask, workers=1)
    expect = math.degrees(math.atan(-math.tan(run.DX_KM / R / 2)))   # the nearest sample, one pixel away
    assert np.allclose(hz, expect, atol=0.01)
    rho = np.hypot(*np.meshgrid(g['x'], g['x']))
    wall = surf.copy()
    wall[(rho > 18.5) & (rho < 23.0)] = 1.0                 # a 1 km ring wall
    hz, *_ = run.horizons(wall, g, mask, workers=1)
    d = 10 * run.DX_KM                                      # the wall's nearest sample on the grid-aligned rays
    a = d / R
    expect = math.degrees(math.atan(((R + 1.0) * math.cos(a) - R) / ((R + 1.0) * math.sin(a))))
    axes = [0, run.N_AZ // 4, run.N_AZ // 2, 3 * run.N_AZ // 4]
    assert np.allclose(hz[0, axes], expect, atol=0.02)
    nine = math.degrees(math.atan(((R + 1.0) * math.cos(0.9 * a) - R) / ((R + 1.0) * math.sin(0.9 * a))))
    assert np.all((hz > 2.3) & (hz < nine + 0.01))          # off-axis rays meet the wall at their rounded samples


def test_sunlit_hours_open_ground(small_grid):
    g = small_grid
    idx = np.flatnonzero((np.hypot(*np.meshgrid(g['x'], g['x'])) < 10.0).ravel())
    n0 = g['nvec'].reshape(-1, 3)[idx]
    north, east = run.local_axes(n0)
    hz = np.full((idx.size, run.N_AZ), -0.5, np.float32)
    lit, flat = run.sunlit_hours(hz, n0, north, east, [1.0], steps=72)
    assert np.all(lit >= flat)
    hz[:] = 30.0
    lit, flat = run.sunlit_hours(hz, n0, north, east, [1.0], steps=72)
    assert np.all(lit == 0) and np.all(flat > 0)


def test_blocked_fraction():
    tab = run.sky_table()
    z = np.zeros((3, run.N_AZ), np.float32)
    z[1] = 10.0
    z[2] = 30.0
    b = run.blocked_fraction(tab, z, -10.0, np.zeros(3))
    assert b[0] == 0.0
    assert abs(b[1] / math.sin(math.radians(10.0)) ** 2 - 1) < 0.3
    assert b[1] < b[2] < 0.5
    half = np.zeros((2, run.N_AZ), np.float32)
    az = np.arange(run.N_AZ) * 360.0 / run.N_AZ
    half[0, (az < 90) | (az > 270)] = 15.0                  # a wall toward the Sun (azimuth 0)
    half[1, (az > 90) & (az < 270)] = 15.0                  # a wall away from it
    b = run.blocked_fraction(tab, half, -10.0, np.zeros(2))
    assert b[0] > b[1]


def test_terrain_section(product):
    for cap in ('north', 'south'):
        t = product['terrain'][cap]
        assert abs(t['land_share'] + t['sea_share'] + t['lake_share'] - 1) < 1e-3     # shares rounded to 4 digits
        y = t['no_direct_sun_all_year']
        e = t['direct_sun']['equinox']['no_direct_sun_share']
        assert y['tilt_2500']['60-90'] <= y['tilt_1613']['60-90'] <= e['60-90']
        assert y['tilt_1613_refracted']['60-90'] <= y['tilt_1613']['60-90']
