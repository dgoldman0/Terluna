"""Checks of the sea-scene kernels against the reflection model they share with the sea-appearance study.

    <numba environment>/bin/python -m pytest visualization/reference-renderer/seas/test_kernels.py
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from numba import njit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))
import kernels as k  # noqa: E402
from illumination.water_surface import reflection  # noqa: E402

N_WATER = 1.34


def analytic_sky(d):
    d = np.atleast_2d(d)
    z = np.clip(d[:, 2], 0, 1)
    return np.stack([1 + 2 * z + 0.5 * d[:, 0], 1.5 + 3 * z, 2 + z - 0.3 * d[:, 1]], -1)


def sky_map():
    el = np.r_[np.arange(-1.0, 1.0, 0.05), np.arange(1.0, 90.01, 0.25)]
    az = np.arange(0.0, 360.0, 0.25)
    e, a = np.meshgrid(np.radians(el), np.radians(az), indexing="ij")
    d = np.stack([np.cos(e) * np.sin(a), np.cos(e) * np.cos(a), np.sin(e)], -1)
    return analytic_sky(d.reshape(-1, 3)).reshape(len(el), len(az), 3), el, 0.25


def scene_arrays(cov, earth=None, water=(0.0, 0.0, 0.0)):
    sky, sky_el, step = sky_map()
    params = np.zeros(24)
    params[k.P_R] = 1737.4e3
    params[k.P_SEA_DX] = 1.0
    params[k.P_SEA_HMIN], params[k.P_SEA_HMAX] = -0.01, 0.01
    params[k.P_N] = N_WATER
    params[k.P_F0] = reflection.fresnel(1.0, N_WATER)
    params[k.P_SKY_AZ_STEP] = step
    params[k.P_LUT_AZ_STEP] = 5.0
    params[k.P_PIXEL] = 1e-3
    params[k.P_INTERREFLECT] = 0.0
    params[k.P_SEA_BLOCK], params[k.P_TER_BLOCK] = 4, 4
    params[k.P_SUN_TEST] = 0.2
    lean = np.zeros(5, np.float32)
    lean_offsets = np.zeros((1, 1), np.int64)
    lean_sizes = np.ones((1, 1), np.int64)
    residual = np.array([cov[0, 0], cov[1, 1], cov[0, 1]])
    sea_h = np.zeros((8, 8), np.float32)
    sea_blocks = np.zeros((2, 2), np.float32)
    ter = np.zeros((2, 2), np.float32)
    tex = np.zeros((2, 2, 1, 3))
    lut_s = np.zeros((2, 72, 2, 3))
    lut_t = np.ones((2, 2, 3))
    lut_tb = np.ones((2, 2, 1))
    lut_el, lut_d = np.array([-90.0, 90.0]), np.array([0.0, 1e6])
    sun = np.zeros(10)
    sun[2] = -1.0
    earth = np.zeros(12) if earth is None else earth
    return (params, sea_h, sea_blocks, lean, lean_offsets, lean_sizes, np.array([1.0]), np.array([1], np.int64),
            residual, ter, ter, tex, sky, sky_el, lut_s, lut_tb, lut_el, lut_d, sun, earth, np.array(water, float))


@njit(cache=True)
def average_shade(v_cam, samples, seed, params, sea_h, sea_blocks, lean, lean_offsets, lean_sizes, lean_spacing,
                  lean_levels, residual, ter_h, ter_blocks, tex, sky, sky_el, lut_s, lut_tb, lut_el, lut_d, sun, earth,
                  water):
    np.random.seed(seed)
    p = np.zeros(3)
    d_in = -v_cam
    acc = np.zeros(3)
    for _ in range(samples):
        u = np.random.random(4)
        acc += k.shade_sea(p, d_in, 1.0, 1.0, 0, u[0], u[1], u[2], u[3], params, sea_h, sea_blocks, lean,
                           lean_offsets, lean_sizes, lean_spacing, lean_levels, residual, ter_h, ter_blocks, tex, sky,
                           sky_el, lut_s, lut_tb, lut_el, lut_d, sun, earth, water)
    return acc / samples


def view(elevation_deg, azimuth_deg):
    e, a = math.radians(elevation_deg), math.radians(azimuth_deg)
    return np.array([math.cos(e) * math.sin(a), math.cos(e) * math.cos(a), math.sin(e)])


@pytest.mark.parametrize("elevation,azimuth", [(30.0, 0.0), (8.0, 70.0), (2.0, 200.0), (0.5, 300.0)])
def test_visible_facet_sampling_reflects_the_sky_as_the_study_does(elevation, azimuth):
    cov = np.array([[0.020, 0.004], [0.004, 0.011]])
    args = scene_arrays(cov)
    v = view(elevation, azimuth)
    rendered = average_shade(v, 400000, 7, *args)
    expected, _ = reflection.sky_reflection(v, cov, N_WATER, analytic_sky, order=48)
    np.testing.assert_allclose(rendered, expected, rtol=0.01)


def test_tilted_mean_samples_the_visible_facets():
    """With a mean slope the sampled facets follow p(s) (v.n) G1 over the sheared distribution."""
    mean, cov = np.array([0.05, -0.03]), np.array([0.012, 0.007, 0.002])
    v = view(6.0, 30.0)

    @njit(cache=True)
    def draw(v, mean, cov, n, seed):
        np.random.seed(seed)
        out = np.zeros((n, 2))
        for i in range(n):
            out[i] = k.sample_visible_slope(v, mean, cov, np.random.random(), np.random.random())
        return out

    s = draw(v, mean, cov, 400000, 3)
    c2 = np.array([[cov[0], cov[2]], [cov[2], cov[1]]])
    nodes, weights = reflection.facet_nodes(c2, 60)
    nodes = nodes + mean
    w = weights * np.maximum(v[2] - nodes @ v[:2], 0.0)
    w /= w.sum()
    for f in (lambda g: g[:, 0], lambda g: g[:, 1], lambda g: g[:, 0] ** 2, lambda g: g[:, 0] * g[:, 1]):
        assert np.mean(f(s)) == pytest.approx(np.sum(w * f(nodes)), rel=0.02, abs=2e-4)


@njit(cache=True)
def earth_estimators(v, mean, cov3, earth, n, samples, seed):
    """The Earth's single-reflection glitter by sampling facets and by sampling the disk, separately."""
    np.random.seed(seed)
    nee, brdf = np.zeros(3), np.zeros(3)
    lam_v = k.smith_lambda(v, mean, cov3)
    for _ in range(samples):
        u = np.random.random(4)
        gx, gy = k.sample_visible_slope(v, mean, cov3, u[0], u[1])
        m = k.norm3(np.array([-gx, -gy, 1.0]))
        cvm = max(k.dot3(v, m), 1e-9)
        r = np.array([2 * cvm * m[0] - v[0], 2 * cvm * m[1] - v[1], 2 * cvm * m[2] - v[2]])
        if r[2] > 0:
            escape = (1 + lam_v) / (1 + lam_v + k.smith_lambda(r, mean, cov3))
            brdf += k.fresnel(cvm, n) * escape * k.earth_disk(r, earth)
        z = 1 - u[2] * (1 - math.cos(earth[3]))
        phi = 2 * math.pi * u[3]
        sz = math.sqrt(1 - z * z)
        e = np.array([earth[0], earth[1], earth[2]])
        a = k.norm3(np.array([-e[1], e[0], 0.0]))
        b = np.array([e[1] * a[2] - e[2] * a[1], e[2] * a[0] - e[0] * a[2], e[0] * a[1] - e[1] * a[0]])
        l = z * e + sz * math.cos(phi) * a + sz * math.sin(phi) * b
        h = k.norm3(v + l)
        d_h = k.slope_pdf(-h[0] / h[2], -h[1] / h[2], mean, cov3) / h[2] ** 4
        g2 = 1 / (1 + lam_v + k.smith_lambda(l, mean, cov3))
        nee += k.fresnel(k.dot3(v, h), n) * d_h * g2 / (4 * v[2]) * earth[10] * k.earth_disk(l, earth)
    return nee / samples, brdf / samples


def test_earth_glitter_by_either_sampling_matches_the_disk_formula():
    """Both samplings that the shader combines are unbiased for the Earth's single-reflection glitter; the
    shader adds the Earth seen after a second reflection off a level surface, as the study's model does for
    the sky."""
    cov = np.array([[0.015, 0.0], [0.0, 0.009]])
    e = view(4.0, 180.0)
    radius = math.radians(0.95)
    omega = 2 * math.pi * (1 - math.cos(radius))
    radiance = np.array([0.8, 1.0, 1.3])
    earth = np.zeros(12)
    earth[0:3], earth[3] = e, radius
    earth[4:7] = -e                                    # the Sun behind the viewer: a full Earth
    earth[7:10], earth[10], earth[11] = radiance, omega, 1.0
    for v in (view(3.0, 0.0), view(10.0, 5.0)):
        expected = reflection.disk_glint(e, radius, v[None], cov, N_WATER, points=301)[0] * radiance * omega
        nee, brdf = earth_estimators(v, np.zeros(2), np.array([0.015, 0.009, 0.0]), earth, N_WATER, 1000000, 3)
        np.testing.assert_allclose(nee, expected, rtol=0.01)
        if v[2] < 0.1:                                 # facet sampling reaches the disk often only near specular
            np.testing.assert_allclose(brdf, expected, rtol=0.03)
        shaded = (average_shade(v, 400000, 11, *scene_arrays(cov, earth))
                  - average_shade(v, 400000, 11, *scene_arrays(cov)))
        assert np.all(shaded > 0.98 * expected) and np.all(shaded < 1.2 * expected)


@njit(cache=True)
def scan(o, d, h, inv2r, t_max, step):
    t = 0.0
    while t < t_max:
        if k.gap(o, d, t, h, 0.0, 0.0, 1.0, True, 0.0, inv2r) <= 0:
            return t
        t += step
    return -1.0


def test_heightfield_tracing_finds_the_first_crossing():
    rng = np.random.default_rng(5)
    h = rng.normal(0, 0.4, (64, 64)).astype(np.float32)
    block = 8
    blocks = h.reshape(8, 8, 8, 8).max(axis=(1, 3))
    # A block's cells reach the first nodes of the next blocks: include them.
    blocks = np.maximum(blocks, np.roll(blocks, -1, 0))
    blocks = np.maximum(blocks, np.roll(blocks, -1, 1))
    blocks = np.maximum(blocks, np.roll(np.roll(h.reshape(8, 8, 8, 8).max(axis=(1, 3)), -1, 0), -1, 1))
    blocks = blocks.astype(np.float32)
    inv2r = 0.5 / 1737.4e3
    for _ in range(300):
        o = np.array([rng.uniform(0, 64), rng.uniform(0, 64), rng.uniform(1.6, 3.0)])
        el = math.radians(-rng.uniform(0.5, 40))
        az = rng.uniform(0, 2 * math.pi)
        d = np.array([math.cos(el) * math.sin(az), math.cos(el) * math.cos(az), math.sin(el)])
        t = k.trace_heightfield(o, d, 0.0, 2000.0, h, blocks, block, 0.0, 0.0, 1.0, True, 0.0, inv2r)
        first = scan(o, d, h, inv2r, 2000.0, 0.002)
        assert first > 0 and t == pytest.approx(first, abs=0.01)
