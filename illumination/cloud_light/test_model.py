"""Independent analytic limits for the spherical molecular beam calculation."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from illumination.cloud_light.model import (
    MolecularColumn, atlas_proxy, cloud_view, photopic_weights, solar_strips,
    solved_moon, sunset_depression_deg,
)
from illumination.sky.atmospheres import MOON_ZERO, rayleigh
from shared.constants import MOON_RADIUS


@pytest.fixture
def shell():
    return MolecularColumn(100., np.array([0., 10., 30.]), np.array([[.1], [.02]]), np.array([550.]))


def test_vertical_column(shell):
    lengths, blocked = shell.path(0., 1.)
    np.testing.assert_allclose(lengths, [10, 20])
    assert not blocked
    assert shell.transmission(0., 1.)[0] == pytest.approx(math.exp(-1.4))


def test_tangent_chords(shell):
    lengths, blocked = shell.path(0., 0.)
    inner, outer = math.sqrt(110**2 - 100**2), math.sqrt(130**2 - 100**2)
    np.testing.assert_allclose(lengths, [inner, outer - inner])
    assert not blocked


def test_occultation_and_finite_segment(shell):
    assert shell.path(5., -1.)[1]
    assert not shell.path(5., -1., 4.)[1]
    assert shell.transmission(5., -1.)[0] == 0
    assert shell.transmission(5., -1., 4.)[0] == pytest.approx(math.exp(-.4))


def test_segment_reciprocity_and_composition(shell):
    a, b = np.array([101., 0]), np.array([111., 14.])
    mid = (a + b) / 2

    def tr(p, q):
        d = q - p
        return shell.transmission(np.linalg.norm(p) - 100, p @ d / np.linalg.norm(p) / np.linalg.norm(d),
                                  np.linalg.norm(d))[0]

    assert tr(a, b) == pytest.approx(tr(b, a), rel=1e-12)
    assert tr(a, b) == pytest.approx(tr(a, mid) * tr(mid, b), rel=1e-12)


def test_solar_disk_vacuum_half_visible():
    vacuum = MolecularColumn(MOON_RADIUS, [0., 100000.], [[0.]], [550.])
    assert vacuum.sun_transmission(0, 0)[0] == pytest.approx(.5)
    assert vacuum.sun_transmission(0, 1)[0] == pytest.approx(1.)
    assert vacuum.sun_transmission(0, -1)[0] == 0
    assert solar_strips(32)[1].sum() == pytest.approx(1.)


def test_geometric_sunset_and_view():
    h = 40000.
    d = sunset_depression_deg(MOON_RADIUS, h)
    assert (MOON_RADIUS + h) * math.cos(math.radians(d)) == pytest.approx(MOON_RADIUS)
    assert sunset_depression_deg(MOON_RADIUS, h, 100000.) - d == pytest.approx(
        math.degrees(100000. / MOON_RADIUS))
    assert cloud_view(MOON_RADIUS, h, 0)['apparent_elevation_deg'] == 90
    view = cloud_view(MOON_RADIUS, h, 400000.)
    assert atlas_proxy().path(1.6, view['mu'], view['ray_length_m'])[1]


@pytest.mark.parametrize('height,mu', [(0., 1.), (40000., 0.), (40000., -.1)])
def test_proxy_against_independent_exponential_quadrature(height, mu):
    model = atlas_proxy(250.)
    r, top = model.radius_m + height, model.radius_m + model.height_m[-1]
    end = -r * mu + math.sqrt(r*r*mu*mu + top*top-r*r)
    column = quad(lambda s: math.exp(-(math.sqrt(r*r + 2*r*mu*s + s*s) - model.radius_m)
                                    / MOON_ZERO.scale_height_m), 0, end, epsrel=1e-10)[0]
    actual = -np.log(model.transmission(height, mu))
    np.testing.assert_allclose(actual, column * rayleigh(model.wavelength_nm, MOON_ZERO), rtol=2e-4)


def test_solved_column_preserves_committed_rayleigh():
    column, facts = solved_moon()
    at550 = list(column.wavelength_nm).index(550)
    tau = float(np.diff(column.height_m) @ column.extinction_m1[:, at550])
    assert tau == pytest.approx(facts['column_facts']['rayleigh_optical_depth_550nm'], abs=1e-11)
    assert tau == pytest.approx(.75674557892616, rel=1e-9)
    assert 0 < photopic_weights(True).sum() < photopic_weights(False).sum()


@pytest.mark.parametrize('height,mu,length', [(-1, 0, 1), (31, 0, 1), (0, 1.1, 1), (0, 0, -1), (0, 0, math.nan)])
def test_invalid_path(shell, height, mu, length):
    with pytest.raises(ValueError):
        shell.path(height, mu, length)
