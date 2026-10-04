"""Short wind waves at any gravity, from the unified spectrum of Elfouhaily et al. (1997).

Elfouhaily, T., Chapron, B., Katsaros, K. and Vandemark, D. (1997), A unified
directional spectrum for long and short wind-driven waves, J. Geophys. Res.
102(C7), 15781-15796. Equation numbers below are the paper's.

The spectrum holds gravity and surface tension explicitly. Its dispersion is
omega² = g k (1 + (k/k_m)²) (24), with the gravity-capillary wavenumber
k_m = sqrt(rho g / T) and the minimum phase speed c_m = sqrt(2 g / k_m). The
short waves' level follows u*/c_m (44) and the long waves' the inverse wave age
Omega = U10/c_p (31-34); the directional spread follows both (57, 59). The paper
fitted these laws on Earth. At lunar gravity they hold if the same dimensionless
controls set the waves there, and the lunar seas' u*/c_m stays inside the range
the fit spans. The paper's k_m = 370 rad/m at Earth's gravity fixes the water's
surface tension over density.
"""
from __future__ import annotations

import numpy as np

from shared.constants import STANDARD_GRAVITY

EARTH_CAPILLARY_WAVENUMBER = 370.0
TENSION_OVER_DENSITY = STANDARD_GRAVITY / EARTH_CAPILLARY_WAVENUMBER ** 2
FULLY_DEVELOPED = 0.84


def capillary_scales(gravity):
    """The gravity-capillary wavenumber k_m (rad/m) and the minimum phase speed c_m (m/s)."""
    k_m = np.sqrt(gravity / TENSION_OVER_DENSITY)
    return float(k_m), float(np.sqrt(2 * gravity / k_m))


def phase_speed(k, gravity):
    """c(k) from the dispersion relation (24), for deep water."""
    k_m, _ = capillary_scales(gravity)
    k = np.asarray(k, float)
    return np.sqrt(gravity / k * (1 + (k / k_m) ** 2))


def short_wave_level(u_star, gravity):
    """alpha_m (44). The fit reaches zero at u*/c_m = 1/e; lighter stress leaves no short waves."""
    _, c_m = capillary_scales(gravity)
    ratio = np.asarray(u_star, float) / c_m
    with np.errstate(divide="ignore"):
        level = 1e-2 * np.where(ratio < 1, 1 + np.log(ratio), 1 + 3 * np.log(ratio))
    return np.maximum(level, 0.0)


def _peak(u10, omega, gravity):
    if not FULLY_DEVELOPED <= omega <= 5:
        raise ValueError("The spectrum covers inverse wave ages from 0.84 to 5")
    if u10 <= 0:
        raise ValueError("A positive 10 m wind is required")
    k_p = gravity * omega ** 2 / u10 ** 2
    return k_p, phase_speed(k_p, gravity)


def curvature(k, u10, u_star, omega, gravity):
    """Long- and short-wave curvature spectra B_l and B_h (31, 40); B = k³ S(k), S the elevation spectrum."""
    k = np.asarray(k, float)
    k_m, c_m = capillary_scales(gravity)
    k_p, c_p = _peak(u10, omega, gravity)
    c = phase_speed(k, gravity)
    gamma = 1.7 if omega <= 1 else 1.7 + 6 * np.log10(omega)
    sigma = 0.08 * (1 + 4 * omega ** -3)
    shape = np.exp(-(np.sqrt(k / k_p) - 1) ** 2 / (2 * sigma ** 2))
    side = np.exp(-1.25 * (k_p / k) ** 2) * gamma ** shape * np.exp(-omega / np.sqrt(10) * (np.sqrt(k / k_p) - 1))
    long_waves = 0.5 * 6e-3 * np.sqrt(omega) * c_p / c * side
    short_waves = 0.5 * short_wave_level(u_star, gravity) * c_m / c * np.exp(-0.25 * (k / k_m - 1) ** 2)
    return long_waves, short_waves


def spreading(k, u10, u_star, omega, gravity):
    """The upwind-crosswind ratio Delta(k) (57, 59); the spread is [1 + Delta cos 2phi] / 2pi (49)."""
    _, c_m = capillary_scales(gravity)
    _, c_p = _peak(u10, omega, gravity)
    c = phase_speed(k, gravity)
    return np.tanh(np.log(2) / 4 + 4 * (c / c_p) ** 2.5 + 0.13 * u_star / c_m * (c_m / c) ** 2.5)


def slope_variances(u10, u_star, omega, gravity, k_from=None, k_to=None, points=4000):
    """Slope variances along and across the wind from the waves between k_from and k_to (rad/m).

    Each is the integral of k² S(k) (1/2 ± Delta/4) over k (58), taken in ln k
    by the trapezoidal rule. The defaults run from a tenth of the peak wavenumber
    to 15 k_m, where the capillary roll-off has fallen by e^-49.
    """
    k_m, _ = capillary_scales(gravity)
    k_p, _ = _peak(u10, omega, gravity)
    k = np.geomspace(k_p / 10 if k_from is None else k_from, 15 * k_m if k_to is None else k_to, points)
    long_waves, short_waves = curvature(k, u10, u_star, omega, gravity)
    delta = spreading(k, u10, u_star, omega, gravity)
    log_k = np.log(k)
    b = long_waves + short_waves
    total = float(np.trapezoid(b, log_k))
    difference = float(0.5 * np.trapezoid(delta * b, log_k))
    return dict(along=(total + difference) / 2, across=(total - difference) / 2, total=total,
                long_wave_part=float(np.trapezoid(long_waves, log_k)),
                short_wave_part=float(np.trapezoid(short_waves, log_k)))


def short_wave_slopes(u_star, gravity, k_from, peak_speed=np.inf, points=4000):
    """Slope variances along and across the wind from the short-wave regime alone (40, 41), above k_from.

    For use above the highest wavenumber of a wave model that resolves the longer
    waves itself: the long-wave regime (31), which in light winds carries more
    slope than Cox and Munk measured, is left out. The resolved sea's peak phase
    speed enters only the spread (57).
    """
    k_m, c_m = capillary_scales(gravity)
    k = np.geomspace(k_from, 15 * k_m, points)
    c = phase_speed(k, gravity)
    level = float(short_wave_level(u_star, gravity))
    b = 0.5 * level * c_m / c * np.exp(-0.25 * (k / k_m - 1) ** 2)
    delta = np.tanh(np.log(2) / 4 + 4 * (c / peak_speed) ** 2.5 + 0.13 * u_star / c_m * (c_m / c) ** 2.5)
    log_k = np.log(k)
    total = float(np.trapezoid(b, log_k))
    difference = float(0.5 * np.trapezoid(delta * b, log_k))
    return dict(along=(total + difference) / 2, across=(total - difference) / 2, total=total, level=level)


def neutral_wind(u_star, gravity, charnock, roughness_floor, von_karman, height=10.0):
    """Wind at a height over the sea in a neutral surface layer with Charnock's roughness."""
    u_star = np.asarray(u_star, float)
    roughness = np.maximum(charnock * u_star ** 2 / gravity, roughness_floor)
    return u_star / von_karman * np.log(height / roughness)
