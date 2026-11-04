import json
import math
from pathlib import Path

import numpy as np
import pytest

from geology import thermal_state as ts
from geology import viscoelastic as ve
from shared.constants import MOON_LOVE_H2, MOON_LOVE_K2, MOON_RADIUS

REF = json.loads((Path(ve.__file__).parent / "interior_sources.json").read_text())["crustal_loading_reference"]


def sphere(eta=1e21, mu=65.6e9, rho=3344.0):
    return ve.Moon([ve.Shell(MOON_RADIUS, rho, mu, eta)])


@pytest.mark.parametrize("n", [2, 5])
def test_homogeneous_sphere_limits(n):
    moon = sphere()
    closed = ve.homogeneous_maxwell(n, moon.mu[0], 1e21 / 65.6e9 / ve.YEAR_S)
    v = moon.love(n, np.array([np.inf, 1e-14]))
    for kind in ("k_load", "h_load", "k", "h"):
        assert v[kind][0].real == pytest.approx(closed[kind][0], rel=1e-4)
        assert v[kind][1].real == pytest.approx(closed[kind][1], rel=1e-4)


def test_homogeneous_sphere_relaxation_time():
    """The fitted spectrum and Gaver-Stehfest both recover the single mode of a homogeneous Maxwell sphere."""
    moon = sphere()
    closed = ve.homogeneous_maxwell(2, moon.mu[0], 1e21 / 65.6e9 / ve.YEAR_S)
    tau = closed["tau_years"]
    elastic, fluid = closed["k_load"]
    t = np.array([0.3, 1.0, 3.0]) * tau
    exact = fluid + (elastic - fluid) * np.exp(-t / tau)
    fitted = ve.response(moon, 2, t, ("k_load",))["k_load"]["series"]
    assert np.max(np.abs(fitted - exact)) < 3e-3
    stehfest = [ve.stehfest_step(moon.love(2, ve.stehfest_points(x))["k_load"].real) for x in t]
    assert np.max(np.abs(np.array(stehfest) - exact)) < 3e-3


def test_elastic_limit_reproduces_crustal_loading_study():
    stack = ts.layers()
    compressible = ve.Moon([ve.Shell(t, d, mu, math.inf, None if mu == 0 else 5 * mu / 3) for t, d, mu in stack])
    for i, n in enumerate(REF["elastic_load"]["n"]):
        v = compressible.love(n, np.array([np.inf]))
        tol = 2e-3 if n <= 20 else 1e-2                     # the study interpolates above degree 40
        assert v["k_load"][0].real == pytest.approx(REF["elastic_load"]["k_load"][i], rel=tol)
        assert v["h_load"][0].real == pytest.approx(REF["elastic_load"]["h_load"][i], rel=tol)
    v = compressible.love(2, np.array([np.inf]))
    assert v["k"][0].real == pytest.approx(REF["tidal"]["compressible"]["k2"], rel=1e-3)
    assert v["h"][0].real == pytest.approx(REF["tidal"]["compressible"]["h2"], rel=1e-3)
    incompressible = ve.Moon([ve.Shell(t, d, mu) for t, d, mu in stack])
    v = incompressible.love(2, np.array([np.inf]))
    assert v["k"][0].real == pytest.approx(REF["tidal"]["incompressible"]["k2"], rel=1e-3)
    assert v["h"][0].real == pytest.approx(REF["tidal"]["incompressible"]["h2"], rel=1e-3)


def test_monthly_tide_near_measured():
    """With the low-viscosity zone at 2e16 Pa s and an inviscid core the monthly k2 and h2 sit a few per cent above the
    measured values: the interior is 1% softer than the Moon when elastic, and the zone's Maxwell relaxation within a
    month adds about 4% (README)."""
    moon = ve.Moon([ve.Shell(*s) for s in ts.shells(ts.CENTRAL)])
    v = moon.love(2, np.array([1j * 2 * math.pi * 365.25 / 27.321661]))
    k2, h2 = abs(complex(v["k"][0])), abs(complex(v["h"][0]))
    assert MOON_LOVE_K2 < k2 < 1.07 * MOON_LOVE_K2
    assert MOON_LOVE_H2 < h2 < 1.2 * MOON_LOVE_H2


def test_step_and_fluid_threshold_converge():
    moon = ve.Moon([ve.Shell(*s) for s in ts.shells(ts.CENTRAL)])
    s = np.array([np.inf, 1e-9, 1e-5, 1e-1])
    for n in (2, 10, 20):
        coarse = moon.love(n, s)["k_load"].real
        fine = moon.love(n, s, step=0.0025)["k_load"].real
        assert np.max(np.abs(coarse - fine) / np.abs(fine)) < 1e-3
    saved = ve.FLUID_FRACTION
    try:
        ve.FLUID_FRACTION = 1e-12
        strict = moon.love(2, s)["k_load"].real
    finally:
        ve.FLUID_FRACTION = saved
    assert np.max(np.abs(strict - moon.love(2, s)["k_load"].real)) < 1e-6


def test_elastic_lid_over_fluid_interior():
    """A stiff lid over an interior that has flowed approaches the thin-shell relaxed values of the crustal-loading
    study; an incompressible lid of E/3 differs from its thin shell's nu = 0.25 by tens of per cent at degree 2."""
    from geology.interior_response import lid_moon
    thin = REF["relaxed_thin_shell"]
    for te in (40, 150):
        moon = lid_moon(te * 1e3)
        for i, n in enumerate(thin["n"][:3]):
            k = moon.love(n, np.array([1e-13]))["k_load"][0].real
            assert k == pytest.approx(thin[f"te_{te}_km"]["k_load"][i], rel=0.35)
    weak, stiff = lid_moon(40e3).love(2, np.array([1e-13]))["k_load"][0].real, lid_moon(400e3).love(2, np.array([1e-13]))["k_load"][0].real
    assert weak < stiff < 0
