"""The turned climatology carried back onto today's grid: exact without a turn, faithful for a smooth field under the
turn's size of rotation, and its area mean kept."""
import numpy as np

from climate.gcm import turned
from climate.gcm import turned_climatology as tc

LAT = np.degrees(np.arcsin(np.polynomial.legendre.leggauss(32)[0]))[::-1]     # T21's Gaussian rows, from the north
LON = np.arange(64) * 5.625


def rotation(angle_deg, axis):
    """A proper rotation by angle_deg about axis (Rodrigues)."""
    k = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    a = np.radians(angle_deg)
    return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K


def climatology(f):
    """A climatology holding f(lat, lon) in every latitude-longitude field."""
    glat, glon = np.meshgrid(LAT, LON, indexing='ij')
    field = f(glat, glon)
    return dict(lat=LAT, lon=LON, **{k: field for k in tc.LATLON})


def test_without_a_turn_the_fields_come_back_unchanged():
    rng = np.random.default_rng(1)
    clim = dict(lat=LAT, lon=LON, **{k: rng.random((32, 64)) for k in tc.LATLON})
    out = tc.to_today(clim, np.eye(3))
    for k in tc.LATLON:
        assert np.allclose(out[k], clim[k], atol=1e-6)


def test_a_smooth_field_returns_to_its_place_today():
    R = rotation(8.05, [0.4, -0.7, 0.3])
    a = np.array([0.3, -0.5, 0.8])
    # The field belongs to the ground: in today's frame it is a . n; the turned run holds it at the turned places.
    today = lambda lat, lon: turned.unit(lat, lon) @ a                                      # noqa: E731
    turned_field = lambda lat, lon: turned.unit(*turned.to_today(lat, lon, R)) @ a          # noqa: E731
    out = tc.to_today(climatology(turned_field), R)
    glat, glon = np.meshgrid(LAT, LON, indexing='ij')
    error = np.abs(out['pr_mm_day'] - today(glat, glon))
    # Bilinear on the T21 grid keeps the field within 0.25% of its range away from the poles; in the two outermost rows,
    # where turned places fall beyond the last Gaussian row and take its values, within 2%.
    assert error[2:-2].max() < 2.5e-3 and error.max() < 2e-2


def test_the_area_mean_is_kept():
    R = rotation(8.05, [0.1, 0.9, -0.2])
    f = lambda lat, lon: 2.0 + np.cos(np.radians(lat)) ** 2 * (1 + 0.3 * np.sin(np.radians(2 * lon)))  # noqa: E731
    clim = climatology(f)
    out = tc.to_today(clim, R)
    w = np.cos(np.radians(LAT))[:, None] * np.ones((32, 64))
    mean = lambda g: (g * w).sum() / w.sum()                                                # noqa: E731
    assert abs(mean(out['pr_mm_day']) / mean(clim['pr_mm_day']) - 1) < 5e-3
