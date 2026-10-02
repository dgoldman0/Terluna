"""Checks of the collision and charging integrals over exponential size distributions."""
from math import gamma

import numpy as np
import pytest

from atmosphere.electricity import charging as ch


def test_the_collision_rate_matches_the_closed_form_when_only_graupel_falls():
    n_g, lam_g, n_x, lam_x, c, b = 1.0e3, 1.0 / 2.0e-3, 1.0e4, 1.0 / 0.3e-3, 10.0, 0.37
    rate = ch.collision_rate(n_g, lam_g, n_x, lam_x, lambda d: c * d ** b, lambda d: 0.0 * d)
    # N0 = N lam; expand (D_g + D_x)^2 and integrate each term against both exponentials
    n0g, n0x = n_g * lam_g, n_x * lam_x
    m = lambda k, lam: gamma(k + 1.0) / lam ** (k + 1.0)
    exact = 0.25 * np.pi * c * n0g * n0x * (m(2 + b, lam_g) * m(0, lam_x) + 2 * m(1 + b, lam_g) * m(1, lam_x)
                                            + m(b, lam_g) * m(2, lam_x))
    assert rate == pytest.approx(exact, rel=1e-4)               # D^0.37 is not smooth at 0; 24 nodes give 2e-5


def test_the_charging_rate_follows_a_brute_force_double_integral_with_both_falling():
    n_g, lam_g, n_x, lam_x = 500.0, 1.0 / 1.5e-3, 3.0e4, 1.0 / 0.2e-3
    vg, vx = (lambda d: 8.5 * d ** 0.37), (lambda d: 210.0 * d)
    law = ch.power_law(1.0e-15, 100.0e-6, 3.0, 2.0, 2.5)
    rate = ch.pair_rate(n_g, lam_g, n_x, lam_x, vg, vx, law, 0.9)
    dg = np.linspace(1e-7, 25.0 / lam_g, 3001)[:, None]
    dx = np.linspace(1e-7, 25.0 / lam_x, 3001)[None, :]
    dv = np.abs(vg(dg) - vx(dx))
    f = 0.25 * np.pi * (dg + dx) ** 2 * dv * 0.9 * law(dg, dx, dv) * n_g * lam_g * np.exp(-lam_g * dg) \
        * n_x * lam_x * np.exp(-lam_x * dx)
    brute = np.trapezoid(np.trapezoid(f, dx[0], axis=1), dg[:, 0])
    assert rate == pytest.approx(brute, rel=0.02)


def test_cells_broadcast_and_a_slower_fall_cuts_charging_by_the_speed_power():
    lam_g, lam_x = np.array([1.0 / 1.5e-3, 1.0 / 3.0e-3]), np.array([1.0 / 0.2e-3, 1.0 / 0.3e-3])
    vg = lambda d: 19.3 * d ** 0.37
    law = ch.power_law(1.0e-15, 100.0e-6, 3.0, 2.0, 2.5)
    earth = ch.pair_rate([1e3, 2e3], lam_g, [1e4, 1e4], lam_x, vg, lambda d: 0.0 * d, law)
    moon = ch.pair_rate([1e3, 2e3], lam_g, [1e4, 1e4], lam_x, lambda d: 0.44 * vg(d), lambda d: 0.0 * d, law)
    assert earth.shape == (2,) and moon == pytest.approx(earth * 0.44 ** 3.5)
