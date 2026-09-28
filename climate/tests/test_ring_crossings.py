"""Checks of the ring-crossing comparison's pure pieces: where great circles cross, patches, local time and the
statistics."""
import numpy as np
import pytest

from climate.crm import cm1_run as c
from climate.crm import ring_analysis as ra
from climate.crm import ring_crossings as rc

A, A_PRIME = dict(tilt_deg=45.0, node_deg=0.0), dict(tilt_deg=45.0, node_deg=180.0)
B, EQUATOR = dict(tilt_deg=45.0, node_deg=90.0), dict(tilt_deg=0.0, node_deg=0.0)


def test_a_and_its_mirror_cross_on_the_equator_at_0_and_180():
    at_0, at_180 = rc.crossings(A, A_PRIME)
    assert (at_0['lat_deg'], at_0['lon_deg']) == (0.0, 0.0) and (at_180['lat_deg'], at_180['lon_deg']) == (0.0, 180.0)
    assert at_0['along'] == pytest.approx([0.0, 0.5], abs=1e-9) and at_180['along'] == pytest.approx([0.5, 0.0], abs=1e-9)
    assert rc.heading_deg(A, 0.0) == pytest.approx(45.0) and rc.heading_deg(A_PRIME, 0.5) == pytest.approx(135.0)


def test_rings_of_different_pairs_cross_at_35_degrees_where_the_paths_say():
    points = rc.crossings(A, B)
    assert [p['lon_deg'] for p in points] == pytest.approx([135.0, 315.0])
    assert [p['lat_deg'] for p in points] == pytest.approx([35.264, -35.264], abs=1e-3)
    for p in points:
        for path, along in zip((A, B), p['along']):
            s = 2 * np.pi * along
            tilt = np.radians(path['tilt_deg'])                           # the runner's path at that point
            assert np.degrees(np.arcsin(np.sin(tilt) * np.sin(s))) == pytest.approx(p['lat_deg'], abs=1e-9)
            lon = (path['node_deg'] + np.degrees(np.arctan2(np.cos(tilt) * np.sin(s), np.cos(s)))) % 360.0
            assert lon == pytest.approx(p['lon_deg'], abs=1e-9)


def test_the_equatorial_ring_is_a_great_circle_of_no_tilt():
    assert rc.ring_path({'case': 'ring', 'latitude_deg': 0.0}) == EQUATOR
    assert rc.ring_path({'path': {'tilt_deg': 45, 'node_deg': 90, 'lat': []}}) == B
    with pytest.raises(ValueError):
        rc.ring_path({'case': 'ring_80n', 'latitude_deg': 80.0})
    at_0, at_180 = rc.crossings(EQUATOR, A)
    assert (at_0['lon_deg'], at_180['lon_deg']) == (0.0, 180.0)
    assert at_0['along'] == pytest.approx([0.0, 0.0], abs=1e-9) and at_180['along'] == pytest.approx([0.5, 0.5], abs=1e-9)
    assert rc.heading_deg(EQUATOR, 0.3) == pytest.approx(90.0)
    with pytest.raises(ValueError):
        rc.crossings(A, A)


def test_patches_wrap_round_the_ring():
    assert rc.patch_columns(100, 0.0, 2.5, 1.0).tolist() == [0, 1, 2, 97, 98, 99]
    assert rc.patch_columns(100, 0.5, 1.0, 1.0).tolist() == [49, 50]


def test_local_time_at_a_crossing_is_what_each_ring_gives_its_column():
    nx, day = 360, 100.0
    for path in (A, A_PRIME):
        track = c.tilted_path(nx, path['tilt_deg'], path['node_deg'])
        geo = dict(x=(np.arange(nx) + 0.5) * 1000.0, sun=dict(h0=-90.0 + path['node_deg'], length=nx * 1000.0,
                                                              dlon=track['dlon']))
        for t in (0.0, 37.0):
            h = ra.hour_angle(geo, t, day)
            expected = [rc.local_hour_angle(-90.0, lon, t, day) for lon in track['lon']]
            assert np.allclose((h - expected + 180.0) % 360.0 - 180.0, 0.0, atol=1e-6)
    assert rc.local_hour_angle(-90.0, 90.0, 0.0, day) == pytest.approx(0.0)       # noon at 90 E at the start
    assert rc.local_hour_angle(-90.0, 0.0, 25.0, day) == pytest.approx(0.0)


def test_dewpoint_inverts_saturation_and_comfort_takes_both_bounds():
    assert rc.dewpoint_c(ra.saturation_pa(293.15)) == pytest.approx(20.0, abs=1e-9)
    assert rc.dewpoint_c(611.2) == pytest.approx(0.0, abs=1e-12)
    air, dew = np.array([18.0, 26.0, 26.1, 22.0, 22.0]), np.array([10.0, 15.0, 10.0, 15.1, 16.0])
    assert rc.comfortable(air, dew, 'strict').tolist() == [True, True, False, False, False]
    assert rc.comfortable(air, dew, 'loose').tolist() == [True, True, True, True, True]


def test_block_bootstrap_keeps_the_mean_and_is_repeatable():
    assert rc.block_bootstrap(np.full(40, 2.0)) == dict(mean=2.0, p5=2.0, p95=2.0)
    assert rc.block_bootstrap([]) is None
    noisy = np.random.default_rng(1).normal(0.0, 1.0, 400)
    first, again = rc.block_bootstrap(noisy), rc.block_bootstrap(noisy)
    assert first == again and first['p5'] < first['mean'] < first['p95']
    assert first['p95'] - first['p5'] == pytest.approx(2 * 1.645 / np.sqrt(400), rel=0.35)
    persistent = np.repeat(noisy[:50], 8)                                # the same draws, each lasting a day
    wide = rc.block_bootstrap(persistent)
    assert wide['p95'] - wide['p5'] > 2.0 * (first['p95'] - first['p5'])


def test_compare_gives_each_ring_and_their_difference_by_local_time():
    hour = np.linspace(-180.0, 180.0, 240, endpoint=False)
    base = 20.0 + 5.0 * np.cos(np.radians(hour))
    a = {key: np.zeros(hour.size) for key in rc.VARIABLES}
    a['air_c'], a['comfortable_strict'] = base, (base > 24.0).astype(float)
    b = {key: v.copy() for key, v in a.items()}
    b['air_c'] = base - 1.0
    out = rc.compare(hour, {'first': a, 'second': b}, lunar_day_h=708.0)
    air = out['air_c']
    assert air['difference']['all'] == dict(mean=pytest.approx(1.0), p5=pytest.approx(1.0), p95=pytest.approx(1.0))
    assert air['difference']['rms_by_local_time'] == pytest.approx(1.0)
    assert air['mean']['first']['day'] > 23.0 > air['mean']['first']['night']
    assert np.argmax(air['by_local_time']['first']) in (4, 5) and out['hour_angle_deg'][4:6] == [-18.0, 18.0]
    hours = out['comfortable_hours_per_lunar_day']['strict']['first']
    assert hours == pytest.approx(708.0 * np.mean(base > 24.0))
    single = rc.compare(hour, {'first': a}, lunar_day_h=708.0)
    assert 'difference' not in single['air_c'] and list(single['air_c']['mean']) == ['first']


def test_values_at_a_sigma_level_interpolate_in_log_pressure():
    prs = np.array([[1000.0, 1200.0], [900.0, 1100.0], [800.0, 1000.0]]) * 100.0
    field = np.array([[0.0, 5.0], [1.0, 6.0], [2.0, 7.0]])
    got = rc.at_sigma(field, prs, np.array([1000.0e2, 1200.0e2]), 0.95)
    w0 = np.log(1000.0 / 950.0) / np.log(1000.0 / 900.0)
    w1 = np.log(1200.0 / 1140.0) / np.log(1200.0 / 1100.0)
    assert got == pytest.approx([w0, 5.0 + w1])
    assert rc.at_sigma(field, prs, np.array([1000.0e2, 1200.0e2]), 0.8)[0] == pytest.approx(2.0)


def test_a_point_on_a_ring_is_found_along_it():
    assert rc.along_ring(A, 0.0, 0.0) == pytest.approx(0.0, abs=1e-9)
    assert rc.along_ring(A, 0.0, 180.0) == pytest.approx(0.5, abs=1e-9)
    assert rc.along_ring(A, 45.0, 90.0) == pytest.approx(0.25, abs=1e-9)
    assert rc.along_ring(A_PRIME, 0.0, 0.0) == pytest.approx(0.5, abs=1e-9)


def test_harmonic_fit_keeps_the_day_and_drops_short_wiggles():
    centres = np.arange(36) * 10.0 - 175.0
    a = np.radians(centres)
    day = 0.3 - 0.8 * np.cos(a) + 0.2 * np.sin(2 * a)
    wiggle = 0.1 * np.cos(9 * a)
    both = np.column_stack([day + wiggle, 2 * day])                           # two levels
    fit = rc.harmonic_fit(both, centres, harmonics=4)
    assert np.allclose(fit[:, 0], day, atol=1e-9) and np.allclose(fit[:, 1], 2 * day, atol=1e-9)


def test_planet_scale_advection_carries_the_long_waves_along_and_up():
    nx, dx, zh = 64, 1000.0, np.array([50.0, 150.0, 250.0])
    x = (np.arange(nx) + 0.5) * dx
    phase = 2 * np.pi * x / (nx * dx)
    wiggle = np.sin(10 * phase)                                           # a storm-sized wave, left out
    d = dict(uinterp=np.full((3, nx), 2.0), winterp=np.broadcast_to(0.1 * np.sin(phase) + 0.05 * wiggle, (3, nx)).copy(),
             th=300.0 + 0.01 * zh[:, None] + np.sin(phase) + 0.5 * wiggle,
             qv=0.01 - 1e-6 * zh[:, None] + 1e-3 * np.cos(phase) + 1e-4 * wiggle)
    i = 7
    dth, dqv = rc.planet_scale_advection(d, [i], 5, dx, zh)
    kx = 2 * np.pi / (nx * dx)
    expect_th = -(2.0 * kx * np.cos(phase[i]) + 0.1 * np.sin(phase[i]) * 0.01)
    expect_qv = -(2.0 * -1e-3 * kx * np.sin(phase[i]) + 0.1 * np.sin(phase[i]) * -1e-6)
    assert dth == pytest.approx(np.full(3, expect_th), rel=1e-9) and dqv == pytest.approx(np.full(3, expect_qv), rel=1e-9)
