"""Checks of the crust's response to the water: the Love numbers against closed forms and limits,
the harmonic filter, the sea-level equation on a synthetic basin and the sill finder."""
import math

import numpy as np
import pytest

from geography import crustal_loading as cl, love_numbers as ln, topography as tp
from shared.constants import GRAVITATIONAL_CONSTANT

A = 1737.4e3
RHO = 3344.0


def homogeneous(mu_ratio, lame=None):
    g = 4 / 3 * math.pi * GRAVITATIONAL_CONSTANT * RHO * A
    return ln.Interior([ln.Layer(A, RHO, mu_ratio * RHO * g * A, lame)]), g


@pytest.mark.parametrize("n", [2, 3, 10, 40])
def test_homogeneous_incompressible_sphere_matches_closed_forms(n):
    interior, _ = homogeneous(2.0)
    got = ln.elastic_love_numbers(interior, n)
    ref = ln.homogeneous_reference(n, 2.0)
    for key in ("h", "k", "h_load", "k_load"):
        assert got[key] == pytest.approx(ref[key], rel=1e-6)


def test_short_loads_follow_boussinesq():
    """At high degree a compressible sphere answers like a half-space: h' -> -(2/3) g a rho (1 - nu) / mu."""
    interior, g = homogeneous(6.0)
    mu = 6.0 * RHO * g * A
    soft = ln.Interior([ln.Layer(A, RHO, mu, mu)])            # lambda = mu: Poisson's ratio 1/4
    h = ln.elastic_love_numbers(soft, 2000)["h_load"]
    assert h == pytest.approx(-(2 / 3) * g * A * RHO * 0.75 / mu, rel=2e-3)


def test_degree_one_is_rigid_for_a_stiff_moon():
    interior, _ = homogeneous(1e4)
    row = ln.elastic_love_numbers(interior, 1)
    assert abs(row["h_load"]) < 1e-4 and abs(row["gauss_residual"]) < 1e-4


def test_moon_tidal_love_numbers_near_measured():
    """The published interior, read as elastic with a fluid core, against GRAIL's k2 and LOLA's h2."""
    row = ln.elastic_love_numbers(ln.Interior(ln.interior_model(compressible=False)), 2)
    assert row["k"] == pytest.approx(0.02405, rel=0.04)
    assert row["h"] == pytest.approx(0.0371, rel=0.15)


def test_relaxed_fluid_sphere_compensates_fully():
    layers = [ln.Layer(A - 1.0, RHO, 0.0), ln.Layer(A, RHO, 0.0)]
    h, k = ln.relaxed_love_numbers(layers, 6, 0.5, 1e11, 0.25, moho_index=0)
    for n in range(2, 7):
        assert h[n] == pytest.approx(-(2 * n + 1) / 3, rel=1e-5)
        assert k[n] == pytest.approx(-1.0, rel=1e-5)


@pytest.mark.parametrize("n", [3, 5, 10])
def test_thin_shell_agrees_with_thick_shell(n):
    """The thin shell against the same lithosphere integrated as a thick elastic layer (75 km) over the
    same relaxed interior, less the thick layer's own compression (its high-degree response)."""
    layers = ln.interior_model()
    thin, _ = ln.relaxed_love_numbers(layers, n, 75e3, 1e11, 0.25)
    thick, _ = ln.relaxed_shell_ode(layers, n, 75e3, 1e11, 0.25)
    assert thick == pytest.approx(thin[n], rel=0.06)


def test_relaxation_time_grows_with_viscosity():
    t1 = ln.relaxation_time_years(10, 1e21, 40e3, 1e11, 0.25, 3380.0, 1.62)
    t2 = ln.relaxation_time_years(10, 1e22, 40e3, 1e11, 0.25, 3380.0, 1.62)
    assert t2 == pytest.approx(10 * t1) and 1e3 < t1 < 1e6


def test_filter_reproduces_and_selects_degrees():
    filt = cl.Filter(90, 180, 40)
    lat = np.pi / 2 - (np.arange(90) + 0.5) * np.pi / 90
    lon = (np.arange(180) + 0.5) * 2 * np.pi / 180
    p3 = tp.legendre_column(2, 40, np.sin(lat), np.cos(lat))[3 - 2]
    p20 = tp.legendre_column(5, 40, np.sin(lat), np.cos(lat))[20 - 5]
    field = np.outer(p3, np.cos(2 * lon)) + 0.5 * np.outer(p20, np.sin(5 * lon))
    ones = np.ones(41)
    only3 = np.zeros(41)
    only3[3] = 1.0
    same, low = filt.apply(field, np.vstack([ones, only3]))
    assert np.sqrt(np.mean((same - field)**2)) < 2e-3 * np.sqrt(np.mean(field**2))
    assert np.max(np.abs(low - np.outer(p3, np.cos(2 * lon)))) < 5e-3


def test_block_mean_and_upsample_keep_a_smooth_field():
    lat = 90 - (np.arange(64) + 0.5) * 180 / 64
    fine = np.cos(np.radians(lat))[:, None] * np.ones((64, 128))
    coarse = cl.block_mean(fine, 4)
    assert coarse.shape == (16, 32)
    back = cl.upsample(coarse, 4)
    assert np.max(np.abs(back[8:-8] - fine[8:-8])) < 0.01


def test_sea_level_equation_conserves_water_and_sinks_like_airy():
    """A wide flat-floored basin with steep walls, local (Airy-like) compensation: with the volume
    held, the water keeps its depth and the floor sinks by rho_w D / rho_m against the far side
    (degree 0 is a uniform shift, so the land rises by the load's mean over the sphere)."""
    nmax, factor = 47, cl.TRANSFORM_FACTOR
    filt = cl.Filter(48, 96, nmax)
    n = np.arange(nmax + 1)
    rho_m = 3300.0
    h = -(2 * n + 1) * cl.MEAN_DENSITY / (3 * rho_m)
    k = np.full(nmax + 1, -1.0)
    response = cl.Response(h, k, filt)
    nlat, nlon = 48 * factor, 96 * factor
    lat = 90 - (np.arange(nlat) + 0.5) * 180 / nlat
    lon = (np.arange(nlon) + 0.5) * 360 / nlon
    cosd = np.outer(np.sin(np.radians(lat)), np.ones(nlon))          # angular distance from the north pole
    height = np.where(cosd > np.cos(np.radians(50.0)), -2000.0, 3000.0)
    row_area = np.cos(np.radians(lat)) * (np.pi / nlat) * (2 * np.pi / nlon) * cl.MOON_RADIUS**2
    level0 = 0.0
    volume = float(row_area @ np.maximum(level0 - height, 0.0).sum(axis=1))
    sol = cl.solve(height, row_area, response, None, level0, volume=volume)
    eff = height + cl.upsample(sol["dev"], factor)
    assert float(row_area @ np.maximum(sol["level"] - eff, 0.0).sum(axis=1)) == pytest.approx(volume, rel=1e-6)
    U = cl.upsample(sol["U"], factor)
    depth = sol["level"] - eff
    basin = height < level0
    mean_depth = float(row_area @ np.where(basin, depth, 0.0).sum(axis=1)) / float(row_area @ basin.sum(axis=1))
    assert mean_depth == pytest.approx(2000.0, rel=1e-3)
    sag = U[nlat - factor * 8:].mean() - U[:factor * 4].mean()
    assert sag == pytest.approx(cl.WATER_DENSITY * depth[:factor * 4].mean() / rho_m, rel=0.02)


def test_pass_level_finds_a_saddle():
    height = np.full((40, 60), 500.0)
    height[10:30, 5:25] = -300.0
    height[10:30, 35:55] = -250.0
    height[19:21, 25:35] = -120.0                     # a channel whose floor sets the sill
    sill = cl.pass_level(height, [(20, 10), (20, 50)], -300.0, 600.0)
    assert sill == pytest.approx(-120.0, abs=1e-3)
