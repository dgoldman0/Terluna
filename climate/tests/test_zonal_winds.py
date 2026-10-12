"""The drifting body's solar day, the parcel routes and the stored wind-vector product."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from climate.gcm import zonal_winds as zw
from shared.constants import JULIAN_DAY, MOON_RADIUS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT_PATH = zw.RESULTS / 'zonal_winds_A28_dim5_moon.json'


def product():
    if not PRODUCT_PATH.is_file():
        pytest.skip('the wind-vector product has not been generated')
    return json.loads(PRODUCT_PATH.read_text())


def arr(x):
    return np.array(x, dtype=float)


def test_sun_speed_is_the_subsolar_point_s_analytic_speed():
    expected = 2.0 * math.pi * (MOON_RADIUS + 10e3) / (SYNODIC_MONTH_DAYS * JULIAN_DAY)
    assert zw.sun_speed(0.0, 10e3) == pytest.approx(expected, rel=1e-12)
    assert zw.sun_speed(0.0, 10e3) == pytest.approx(4.30, abs=0.005)
    assert zw.sun_speed(60.0, 10e3) == pytest.approx(0.5 * expected, rel=1e-12)
    assert zw.sun_speed(90.0, 0.0) == pytest.approx(0.0, abs=1e-12)


def test_day_length_limits():
    vs = zw.sun_speed(20.0, 5e3)
    assert zw.solar_day(0.0, vs) == pytest.approx(SYNODIC_MONTH_DAYS)              # at rest: the lunar day
    assert np.isinf(zw.solar_day(-vs, vs)) and zw.day_rate(-vs, vs) == 0.0           # u = -v_sun holds the hour
    assert zw.solar_day(vs, vs) == pytest.approx(SYNODIC_MONTH_DAYS / 2)              # eastward wind shortens it
    assert zw.solar_day(-0.5 * vs, vs) == pytest.approx(2 * SYNODIC_MONTH_DAYS)      # slower westward lengthens it
    assert zw.day_rate(-0.5 * vs, vs) > 0
    for f, length in ((-1.5, 2 * SYNODIC_MONTH_DAYS), (-2.0, SYNODIC_MONTH_DAYS), (-3.0, SYNODIC_MONTH_DAYS / 2)):
        assert zw.day_rate(f * vs, vs) < 0                                            # faster westward: Sun turns back
        assert zw.solar_day(f * vs, vs) == pytest.approx(length)


def test_hour_angles_follow_the_climatology_convention():
    lon = np.arange(64) * 5.625
    hour, k = zw.hour_angles(np.array([90.0, 0.0]), lon)
    assert hour[0, 16] == 0.0 and hour[0, 32] == 90.0 and hour[0, 0] == -90.0       # east of the Sun is afternoon
    assert hour[1, 63] == pytest.approx(-5.625)
    centres = zw.hour_centres()
    assert np.all(np.abs(hour - centres[k]) <= 5.0 + 1e-9)
    assert centres[0] == -175.0 and centres.size == zw.HOUR_BINS


def test_sun_track_recovers_a_line_rounded_to_the_grid():
    step = 36.31827
    true = (154.98 - step * np.arange(1200)) % 360.0
    rounded = (np.round(true / 5.625) * 5.625) % 360.0
    fitted = zw.sun_track(rounded, step)
    assert np.max(np.abs((fitted - true + 180.0) % 360.0 - 180.0)) < 0.05


def at_rest(radius, factor=0.0):
    """Sun-relative velocity of a body in an eastward wind u = factor * v_sun."""
    return lambda t, la, ho: ((1.0 + factor) * zw.SUN_RATE * radius * np.cos(np.radians(la)), np.zeros_like(la))


def test_integrator_gives_the_exact_period_of_solid_rotation():
    lat0 = np.array([0.0, 30.0, -60.0, 80.0])
    hour0 = np.array([0.0, 45.0, -120.0, 170.0])
    radius = MOON_RADIUS + 20e3 + np.zeros(4)
    res = zw.integrate(at_rest(radius), lat0, hour0, radius, 35 * JULIAN_DAY)
    assert res['first_pass_days'] == pytest.approx(np.full(4, SYNODIC_MONTH_DAYS), rel=1e-6)
    assert np.all(res['first_pass_sign'] == 1) and res['end_lat'] == pytest.approx(lat0, abs=1e-6)
    assert np.all(res['sun_up_share'] > 0.4) and np.all(res['sun_up_share'] < 0.6)
    held = zw.integrate(at_rest(radius, -1.0), lat0, hour0, radius, 35 * JULIAN_DAY)
    assert np.all(np.isnan(held['first_pass_days'])) and held['end_hour'] == pytest.approx(hour0, abs=1e-6)
    back = zw.integrate(at_rest(radius, -2.0), lat0, hour0, radius, 35 * JULIAN_DAY)
    assert back['first_pass_days'] == pytest.approx(np.full(4, SYNODIC_MONTH_DAYS), rel=1e-6)
    assert np.all(back['first_pass_sign'] == -1)
    # 35 days move the hour 427 degrees, forward at rest and backward turned back: count the passes of +90 degrees
    moved = 360.0 * 35 / SYNODIC_MONTH_DAYS
    passes = lambda end: np.abs(np.floor((end - 90.0) / 360.0) - np.floor((hour0 - 90.0) / 360.0))
    assert np.array_equal(res['sunsets'], passes(hour0 + moved)) and np.all(res['western_sunrises'] == 0)
    assert np.array_equal(back['western_sunrises'], passes(hour0 - moved)) and np.all(back['sunsets'] == 0)
    assert res['sunsets'].tolist() == [1, 2, 1, 1]
    fixed = np.cos(np.radians(lat0)) / np.pi                                   # a fixed place's mean light
    assert res['fixed_light'] == pytest.approx(fixed, rel=1e-6)
    two = zw.integrate(at_rest(radius), lat0, hour0, radius, 2 * SYNODIC_MONTH_DAYS * JULIAN_DAY)
    assert two['light'] == pytest.approx(fixed, rel=2e-3)                      # two whole lunar days at rest


def test_gridded_solid_rotation_keeps_its_period():
    lat = np.linspace(85.76, -85.76, 32)                       # north first, as in the model
    hours = zw.hour_centres()
    radius = MOON_RADIUS + 10e3
    spin = 0.5 * zw.SUN_RATE                                   # an eastward solid rotation of half the Sun's rate
    u = np.broadcast_to((spin * radius * np.cos(np.radians(lat)))[None, :, None], (1, 32, 36))
    xyz, rows = zw.cartesian(u, np.zeros_like(u), lat, hours)
    lat0, hour0 = np.array([0.0, 40.0, -70.0]), np.array([10.0, -100.0, 160.0])
    layer = np.zeros(3, dtype=int)

    def velocity(t, la, ho):
        ue, vn = zw.interpolate(xyz, rows, hours[0], 10.0, la, ho, lead=(layer,))
        return ue + zw.SUN_RATE * radius * np.cos(np.radians(la)), vn
    res = zw.integrate(velocity, lat0, hour0, np.full(3, radius), 30 * JULIAN_DAY)
    assert res['first_pass_days'] == pytest.approx(np.full(3, SYNODIC_MONTH_DAYS / 1.5), rel=5e-3)
    assert res['end_lat'] == pytest.approx(lat0, abs=0.05)


def test_interpolation_returns_the_grid_s_own_winds_at_its_nodes():
    lat = np.array([60.0, 20.0, -20.0, -60.0])
    lon = np.arange(8) * 45.0
    rng = np.random.default_rng(1)
    u, v = rng.normal(size=(4, 8)), rng.normal(size=(4, 8))
    xyz, rows = zw.cartesian(u, v, lat, lon)
    la, lo = np.repeat(lat, 8), np.tile(lon, 4)
    ue, vn = zw.interpolate(xyz, rows, 0.0, 45.0, la, lo)
    assert ue == pytest.approx(u.ravel(), abs=1e-12) and vn == pytest.approx(v.ravel(), abs=1e-12)


def test_circling_the_subsolar_point_never_passes_through_all_hours():
    radius = MOON_RADIUS + 5e3 + np.zeros(2)
    turn = 2 * np.pi / (10 * JULIAN_DAY)

    def velocity(t, la, ho):
        w = turn * np.cross(np.array([1.0, 0.0, 0.0]), zw.position(la, ho)) * radius[:, None]
        east, north = zw.east_north(la, ho)
        return (w * east.T).sum(axis=1), (w * north.T).sum(axis=1)
    res = zw.integrate(velocity, np.array([10.0, 30.0]), np.array([0.0, 20.0]), radius, 30 * JULIAN_DAY)
    assert np.all(np.isnan(res['first_pass_days']))
    assert np.all(res['hours'].max(axis=0) < 90.0) and np.all(res['hours'].min(axis=0) > -90.0)


def test_fixed_points_and_their_kinds():
    lat = np.linspace(-30.0, 30.0, 13)
    hours = zw.hour_centres()
    radius = MOON_RADIUS
    x = radius * np.radians((hours[None, :] - 5.0)) * np.cos(np.radians(lat[:, None]))   # metres east of 5 deg
    y = radius * np.radians(lat[:, None] - 2.5) + 0 * x                                  # metres north of 2.5 deg
    a = 1e-6
    for su, sv, kind in ((-a * x, -a * y, 'sink'), (a * x, a * y, 'source'), (a * x, -a * y, 'saddle')):
        found = [f for f in zw.fixed_points(su, sv, lat, hours, radius) if abs(f['hour_deg'] - 5.0) < 30]
        assert len(found) == 1
        f = found[0]
        assert f['kind'] == kind and f['hour_deg'] == pytest.approx(5.0, abs=0.1)
        assert f['lat_deg'] == pytest.approx(2.5, abs=0.1)
        assert f['e_folding_days'] == pytest.approx(1.0 / (a * JULIAN_DAY), rel=0.05)


def test_longest_dwell_counts_the_held_stretch():
    t = np.arange(0.0, 200.0, 0.5)                             # days, records half a day apart
    held = np.where(t < 40.0, 10.0 * np.sin(t), (t - 40.0) * 12.2)   # within 20 degrees for 40 days, then moving
    moving = t * 12.2                                          # the lunar day's pace: 30 degrees in 2.46 days
    dwell = zw.longest_dwell(np.stack([held, moving]).T, 0.5)
    assert dwell[0] == 30 and dwell[1] == 2                    # the longest listed lengths that fit


def test_gathering_measures_nearest_neighbours():
    layer, lat0, hour0 = zw.releases(2)
    hours = np.stack([hour0, hour0 + np.where(layer == 1, 0.0, 0.0)])     # two records, parcels not moving
    lats = np.stack([lat0, np.where(layer == 1, 0.0, lat0)])              # second record: height 1 on the equator
    res = dict(hours=hours.astype(np.float32), lats=lats.astype(np.float32))
    groups = np.stack([layer, np.zeros_like(layer)], axis=1)
    g = zw.gathering(res, groups, np.full(layer.size, MOON_RADIUS), 170.0, 1.0, days=(0, 1))
    assert g['days'] == [0, 1]
    q = zw.position(lat0[layer == 0], hour0[layer == 0])
    angle = np.arccos(np.clip(q @ q.T, -1, 1))
    np.fill_diagonal(angle, np.inf)
    expected = np.median(angle.min(axis=1)) * MOON_RADIUS / 1e3
    assert g['nearest_km_p50'][0] == pytest.approx([expected, expected], abs=0.1)
    assert g['nearest_km_p50'][1][1] == 0.0 and g['within_cell_share'][1][1] == 1.0   # stacked parcels coincide


def test_share_explained_by_cell_means():
    rng = np.random.default_rng(2)
    cells = np.repeat(np.arange(10), 200)
    means = rng.normal(scale=3.0, size=10)
    x = means[cells] + rng.normal(scale=0.1, size=cells.size)
    n = np.bincount(cells).astype(float)[None]
    s = np.bincount(cells, x)[None]
    _, b, _ = zw.share_explained(n, s, np.array([x.size], float), np.array([x.sum()]), np.array([(x * x).sum()]))
    sst = (x * x).sum() - x.sum() ** 2 / x.size
    assert b[0] / sst == pytest.approx(1.0, abs=0.01)
    noise = rng.normal(size=cells.size)
    s = np.bincount(cells, noise)[None]
    _, b, _ = zw.share_explained(n, s, np.array([noise.size], float), np.array([noise.sum()]),
                                 np.array([(noise * noise).sum()]))
    assert abs(b[0] / ((noise * noise).sum() - noise.sum() ** 2 / noise.size)) < 0.01


def test_product_binds_producer_inputs_and_constants():
    d = product()
    assert d['schema'] == zw.SCHEMA
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    for path, h in d['producer']['files'].items():
        assert digest(ROOT / path) == h, path
    assert not constants_changed(d['producer']['constants'])
    assert {'MOON_RADIUS', 'SYNODIC_MONTH_DAYS'} <= set(d['producer']['constants'])
    runs = [ROOT / p for p in d['producer']['inputs']]
    if not all(p.is_file() for p in runs):
        pytest.skip('the raw design run is not on this machine')
    for path, h in d['producer']['inputs'].items():
        assert digest(ROOT / path) == h, path


def test_product_numbers_hold_together():
    d = product()
    lat, z = arr(d['axes']['lat_deg']), arr(d['axes']['height_km'])
    assert arr(d['sun_speed_m_s']) == pytest.approx(zw.sun_speed(lat[None, :], z[:, None] * 1e3), abs=1e-3)
    g, h = d['ground_frame'], d['holding']
    for c in ('u', 'v'):
        p = np.stack([arr(g[f'{c}_p{q}_m_s']) for q in zw.PERCENTILES])
        ok = np.all(np.isfinite(p), axis=0)
        assert np.all(np.diff(p[:, ok], axis=0) >= -1e-9)
    cov = arr(g['coverage'])
    assert np.all((cov >= 0) & (cov <= 1)) and np.all(cov[z >= 10.0] > 0.99)
    hold5, hold1, west, back = (arr(h[k]) for k in ('hold_0_5', 'hold_1', 'westward', 'reversed'))
    ok = np.isfinite(hold1)
    assert np.all(hold5[ok] <= hold1[ok] + 1e-9) and np.all(back[ok] <= west[ok] + 1e-9)
    s = d['solar_day']
    normal, reverse = arr(s['normal']), arr(s['reversed'])
    total = np.nansum(normal, axis=0) + np.nansum(reverse, axis=0)
    assert total[ok] == pytest.approx(np.ones(ok.sum()), abs=0.01)
    assert np.nansum(reverse, axis=0)[ok] == pytest.approx(back[ok], abs=0.005)
    days = np.stack([arr(s[f'p{q}_days']) for q in zw.DAY_PERCENTILES])
    assert np.all(np.diff(days[:, ok], axis=0) >= -1e-6)
    shear = d['shear']
    for key in ('speed_p10_m_s', 'speed_p50_m_s', 'speed_p90_m_s'):
        assert np.all(arr(shear[key])[np.isfinite(arr(shear[key]))] >= 0)
    for band in d['bands']:
        assert band['hold_0_5'] <= band['hold_1'] <= 1 and band['reversed'] <= band['westward'] + 1e-9
        assert sum(band['normal_classes']) + sum(band['reversed_classes']) == pytest.approx(1.0, abs=0.01)


def test_sun_frame_means_average_to_the_ground_frame_means():
    d = product()
    npz = PRODUCT_PATH.with_suffix('.npz')
    if not npz.is_file():
        pytest.skip('the supporting arrays have not been generated')
    with np.load(npz) as a:
        n = a['sun_samples'].astype(float)
    for c in ('u', 'v'):
        sun = arr(d['sun_frame'][f'{c}_mean_m_s'])
        mean = np.nansum(np.nan_to_num(sun) * n, axis=-1) / np.maximum(n.sum(axis=-1), 1)
        ground = arr(d['ground_frame'][f'{c}_mean_m_s'])
        ok = np.isfinite(ground)
        assert mean[ok] == pytest.approx(ground[ok], abs=0.01)


def test_stored_routes_hold_together():
    d = product()
    t = d['trajectories']
    for run in (t['steady'], t['sequence']['long_run']):
        sun_up, light = arr(run['sun_up_share']), arr(run['light'])
        assert np.all((sun_up >= 0) & (sun_up <= 1)) and np.all(light <= sun_up + 1e-9)
        first = arr(run['first_pass_days'])
        assert np.all(first[np.isfinite(first)] > 0)
        assert np.all(arr(run['western_sunrises']) >= 0) and np.all(arr(run['sunsets']) >= 0)
        ratio = arr(run['light_ratio'])
        assert np.all((ratio > 0) & (ratio <= np.pi + 1e-6))                 # pi: held at noon
        g = run['gathering']
        near, cell = arr(g['nearest_km_p50']), arr(g['within_cell_share'])
        assert np.all(near >= 0) and np.all((cell >= 0) & (cell <= 1)) and g['days'][0] == 0
        occ = arr(run['occupancy'])
        assert occ.sum(axis=(1, 2)) == pytest.approx(np.ones(len(t["heights_km"])), abs=0.011)   # 216 rounded cells
    assert sum(t['area_share']) * zw.OCCUPANCY_HOURS == pytest.approx(1.0, abs=0.011)   # 18 rounded shares
    for points in t['fixed_points'].values():
        assert all(p['kind'] in ('sink', 'source', 'saddle') for p in points)


def test_read_year_matches_the_global_wind_speeds():
    pytest.importorskip('netCDF4')
    from climate.gcm import global_winds as gw
    from climate.gcm.exoplasim_run import RUNS
    from atmosphere.radiative_convective.thermodynamics import MOON, earthlike_air
    path = RUNS / 'A28_dim5_moon' / 'model' / 'MOST.00029.nc'
    if not path.is_file():
        pytest.skip('the raw design run is not on this machine')
    air = earthlike_air(json.loads((RUNS / 'A28_dim5_moon' / 'progress.json').read_text())
                        ['configuration']['pressure_pa'], 400.0)
    y = zw.read_year(path, air.gas_constant)
    speed = gw.read_year(path, air.gas_constant, MOON.gm_m3_s2, MOON.radius_m)[0]
    mine = np.hypot(y['u'], y['v']).reshape(zw.HEIGHTS_KM.size, -1)
    ok = np.isfinite(speed)
    assert np.array_equal(np.isfinite(mine), ok)               # the same ground leaves out the same samples
    # interpolating the vector gives at most the interpolated speed, and the same below the lowest layer
    assert np.all(mine[ok] <= speed[ok] + 1e-9)
    assert np.mean(speed[ok] - mine[ok]) < 0.02 * np.mean(speed[ok])
    assert mine[0][ok[0]] == pytest.approx(speed[0][ok[0]], abs=1e-9)
