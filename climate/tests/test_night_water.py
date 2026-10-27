import math

import numpy as np

from climate.gcm import night_water as nw


def test_quarters_by_hour_angle():
    hours = np.array([-180.0, -135.0, -90.0, -45.0, 0.0, 45.0, 90.0, 135.0, 179.9, 180.0, 270.0])
    names = [nw.QUARTERS[k][0] for k in nw.quarter_index(hours)]
    assert names == ['midnight-sunrise', 'midnight-sunrise', 'sunrise-noon', 'sunrise-noon', 'noon-sunset',
                     'noon-sunset', 'sunset-midnight', 'sunset-midnight', 'sunset-midnight', 'midnight-sunrise',
                     'sunrise-noon']


def test_overlap_keeps_totals():
    src = np.array([0.0, 600.0, 1800.0, 3200.0, 5200.0])
    dst = np.array([0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5200.0])
    share = nw.overlap_matrix(src, dst)
    assert np.allclose(share.sum(axis=1), 1.0)
    values = np.array([3.5, -1.2, 0.4, 0.9])
    assert math.isclose(float((values @ share).sum()), float(values.sum()), rel_tol=1e-12)
    assert np.allclose(nw.overlap_matrix(src, src), np.eye(4))


def test_layers_on_bands_keep_the_column():
    z_half = np.array([9000.0, 5000.0, 2000.0])               # top first, the inner half levels of four layers
    values = np.array([0.3, -1.0, 2.0, 4.0])                   # top first
    edges = np.array([0.0, 1000.0, 3000.0, 6000.0, 9000.0, 13000.0])
    banded = nw.on_bands(values, z_half, edges)
    assert math.isclose(float(banded.sum()), float(values.sum()), rel_tol=1e-12)
    assert math.isclose(float(banded[0]), 4.0 / 2.0, rel_tol=1e-12)   # half the lowest layer lies below 1 km


def test_storage_from_window_means():
    """A series a + b cos(h) + c sin(2h), given as means over windows of 36 degrees, gives back each quarter's change
    over its length once its harmonics are undamped."""
    day_s = 29.53 * 86400.0
    window = 36.0
    a, b, c = 5.0, 2.0, -0.7
    centre = np.arange(-180.0, 180.0, 7.3)
    damp = lambda n: math.sin(math.radians(n * window / 2.0)) / math.radians(n * window / 2.0)
    h = np.radians(centre)
    means = a + b * damp(1) * np.cos(h) + c * damp(2) * np.sin(2.0 * h)
    got = nw.quarter_storage(centre, means[:, None], window, day_s)[:, 0]
    exact = lambda x: a + b * math.cos(math.radians(x)) + c * math.sin(math.radians(2.0 * x))
    for (name, lo, hi), value in zip(nw.QUARTERS, got):
        want = (exact(hi) - exact(lo)) / ((hi - lo) / 360.0 * day_s / 86400.0)
        assert math.isclose(value, want, rel_tol=1e-9, abs_tol=1e-12), name


def test_saturation_matches_the_column_copy():
    from climate.gcm import sbm
    c = dict(sbm.PLASIM)
    t, p = np.array([250.0, 290.0, 303.0]), np.array([60000.0, 100000.0, 121500.0])
    got = nw.saturation(t, p, c)
    for k in range(3):
        assert math.isclose(float(got[k]), sbm.saturation(float(t[k]), float(p[k]), c), rel_tol=1e-12)


def test_group_sums_cover_their_bands():
    edges = np.arange(0.0, 10001.0, 1000.0)
    quarter = dict(vapour_mm_day={'x': np.arange(10.0)}, heat_w_m2={'y': np.ones(10)})
    out = nw.group_sums(quarter, edges)
    assert list(out) == ['0.0-4.0 km', '4.0-7.0 km', '7.0-9.0 km']
    assert out['0.0-4.0 km']['vapour_mm_day']['x'] == 0.0 + 1.0 + 2.0 + 3.0
    assert out['7.0-9.0 km']['heat_w_m2']['y'] == 2.0


def test_storage_from_single_steps_needs_no_undamping():
    """Values taken at single steps (the hour angles themselves, no window) give back each quarter's change as they
    stand."""
    day_s = 29.53 * 86400.0
    hours = np.arange(-180.0, 180.0, 9.1)
    values = 3.0 - 1.5 * np.sin(np.radians(hours)) + 0.4 * np.cos(np.radians(3.0 * hours))
    got = nw.quarter_storage(hours, values[:, None], 36.0, day_s, undamp=False)[:, 0]
    exact = lambda x: 3.0 - 1.5 * math.sin(math.radians(x)) + 0.4 * math.cos(math.radians(3.0 * x))
    for (name, lo, hi), value in zip(nw.QUARTERS, got):
        assert math.isclose(value, (exact(hi) - exact(lo)) / ((hi - lo) / 360.0 * day_s / 86400.0), rel_tol=1e-9,
                            abs_tol=1e-12), name
