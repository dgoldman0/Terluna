"""Checks of the rough-water reflection: optics, Smith's shadowing, the flat limit and energy."""
import numpy as np
import pytest

from illumination.water_surface import reflection as R

N = 1.34


def unit(elevation_deg, azimuth_deg):
    e, a = np.radians(elevation_deg), np.radians(azimuth_deg)
    return np.array([np.cos(e) * np.sin(a), np.cos(e) * np.cos(a), np.sin(e)])


def test_seawater_index_reproduces_pure_water_at_the_sodium_line():
    assert R.refractive_index(589.3, 20.0, 0.0) == pytest.approx(1.3330, abs=3e-4)
    assert R.refractive_index(589.3, 20.0, 35.0) - R.refractive_index(589.3, 20.0, 0.0) == pytest.approx(0.0064, abs=3e-4)
    assert R.refractive_index(400.0) > R.refractive_index(700.0)


def test_fresnel_limits():
    assert R.fresnel(1.0, N) == pytest.approx(((N - 1) / (N + 1)) ** 2)
    assert R.fresnel(0.0, N) == pytest.approx(1.0)
    brewster = np.cos(np.arctan(N))
    rs = ((brewster - np.sqrt(N * N - 1 + brewster ** 2)) / (brewster + np.sqrt(N * N - 1 + brewster ** 2))) ** 2
    assert R.fresnel(brewster, N) == pytest.approx(rs / 2)


def test_smith_lambda_known_value_and_limits():
    cov = np.eye(2) * 0.5      # sigma = sqrt(0.5) in every azimuth, so a = cot(theta)
    direction = unit(45.0, 30.0)       # cot 45 = 1
    expected = (np.exp(-1) - np.sqrt(np.pi) * R.erfc(1.0)) / (2 * np.sqrt(np.pi))
    assert R.smith_lambda(direction, cov) == pytest.approx(expected)
    assert R.smith_lambda(unit(90.0, 0.0), cov) == 0.0
    assert R.smith_lambda(unit(89.0, 0.0), np.eye(2) * 1e-4) == pytest.approx(0.0, abs=1e-12)


def test_visible_facet_weights_sum_to_one():
    cov = np.array([[0.02, 0.004], [0.004, 0.012]])
    for elevation in (60.0, 30.0, 10.0):
        _, info = R.sky_reflection(unit(elevation, 200.0), cov, N, lambda d: np.ones((len(d), 1)), order=60)
        assert info["weight_sum"] == pytest.approx(1.0, abs=2e-3)


def test_a_nearly_flat_sea_mirrors_the_sky():
    cov = np.eye(2) * 1e-6
    sky = lambda d: np.c_[1 + d[:, 2], 2 + d[:, 0]]
    view = unit(25.0, 120.0)
    mirror = view * [-1, -1, 1]
    reflected, _ = R.sky_reflection(view, cov, N, sky)
    assert reflected == pytest.approx(R.fresnel(view[2], N) * sky(mirror[None])[0], rel=1e-3)


def test_glitter_reflects_the_fresnel_share_of_a_high_sun():
    """Integrated over every view, the glitter of an overhead source returns about Fresnel's 2%."""
    cov = np.eye(2) * 0.004
    e = np.radians(np.linspace(30, 90, 400))
    a = np.radians(np.linspace(0, 360, 361))
    ee, aa = np.meshgrid(e, a, indexing="ij")
    views = np.stack([np.cos(ee) * np.sin(aa), np.cos(ee) * np.cos(aa), np.sin(ee)], -1)
    k = R.glint(np.array([0, 0, 1.0]), views, cov, N)
    flux = np.trapezoid(np.trapezoid(k * np.sin(ee) * np.cos(ee), a, axis=1), e)
    assert flux == pytest.approx(R.fresnel(1.0, N), rel=0.02)


def test_glitter_matches_facet_sampling_for_a_wide_source():
    """Sample facets as the observer sees them and count reflections into a 6-degree disk."""
    rng = np.random.default_rng(3)
    cov = np.array([[0.03, 0.0], [0.0, 0.02]])
    source, view = unit(35.0, 0.0), unit(30.0, 180.0)
    radius = np.radians(6.0)
    g = rng.multivariate_normal([0, 0], cov, size=4_000_000)
    normal = np.c_[-g, np.ones(len(g))]
    normal /= np.linalg.norm(normal, axis=1, keepdims=True)
    cos_v = normal @ view
    weight = np.maximum(cos_v, 0) / (view[2] * normal[:, 2])
    reflected = 2 * cos_v[:, None] * normal - view
    inside = reflected @ source > np.cos(radius)
    # A disk of radiance L and solid angle omega gives E_n = L omega; no shadowing at these angles.
    solid = 2 * np.pi * (1 - np.cos(radius))
    sampled = np.sum(weight[inside] * R.fresnel(cos_v[inside], N)) / len(g) / solid
    predicted = R.disk_glint(source, radius, view, cov, N, points=401)
    lam = R.smith_lambda(view, cov) + R.smith_lambda(source, cov)
    assert lam < 0.01
    assert predicted * (1 + lam) == pytest.approx(sampled, rel=0.03)


def test_the_batch_agrees_with_single_views_and_converges_in_order():
    cov = np.array([[0.018, 0.003], [0.003, 0.011]])
    sky = lambda d: np.c_[1 + 3 * d[:, 2] ** 0.5 + 0.5 * d[:, 0], 2 - d[:, 1]]
    views = np.array([unit(e, a) for e, a in ((40.0, 10.0), (12.0, 200.0), (3.0, 95.0))])
    single = np.array([R.sky_reflection(v, cov, N, sky, order=24)[0] for v in views])
    batch, square = R.sky_reflection_batch(views, cov, N, sky, order=24, block=2)
    assert batch == pytest.approx(single, rel=1e-6)
    assert np.all(square >= batch ** 2 * (1 - 1e-9))
    fine, _ = R.sky_reflection_batch(views, cov, N, sky, order=48)
    assert batch == pytest.approx(fine, rel=1e-2)       # a 3-degree grazing view of a sky with a kink at the horizon
