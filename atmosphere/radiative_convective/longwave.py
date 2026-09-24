"""Thermal-infrared fluxes by line-by-line integration of the Schwarzschild equation.

Non-scattering, plane-parallel layers with a source function linear in optical
depth; hemispheric fluxes from Gauss-Legendre angular quadrature (a single
diffusivity angle is available for comparison with two-stream codes).
Wavenumbers are cm^-1, radiances W m^-2 sr^-1 (cm^-1)^-1, fluxes W m^-2.
"""
from __future__ import annotations
import math
import numpy as np

from shared.constants import PLANCK, SPEED_OF_LIGHT
from atmosphere.radiative_convective.spectroscopy import C2

_B_FACTOR = 2 * PLANCK * SPEED_OF_LIGHT**2 * 1e8    # 2hc^2 with nu in cm^-1, per cm^-1


def planck(nu, t):
    """Planck radiance per unit wavenumber, W m^-2 sr^-1 (cm^-1)^-1."""
    nu = np.asarray(nu, dtype=float)
    t = np.asarray(t, dtype=float)
    return _B_FACTOR * nu**3 / np.expm1(C2 * nu / t)


def angles(n=4, diffusivity=None):
    """Nodes mu and weights w with flux = 2 pi sum(w mu I); isotropic I gives pi I."""
    if diffusivity is not None:
        mu = np.array([1.0 / diffusivity])
        return mu, np.array([0.5 / mu[0]])
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1), 0.5 * w


def _linear_source_terms(t):
    """Transmission e^-t and the linear-source weight (1 - e^-t - t e^-t)/t."""
    e = np.exp(-t)
    small = t < 1e-4
    safe = np.where(small, 1.0, t)
    c = np.where(small, t / 2 - t * t / 3, (-np.expm1(-t) - t * e) / safe)
    return e, -np.expm1(-t), c


def fluxes(nu, tau, t_levels, t_surface, emissivity=1.0, n_angles=4, diffusivity=None):
    """Upward and downward spectral fluxes at every level.

    tau: vertical optical depth of each layer, shape (layers, nu), surface layer
    first. t_levels: level temperatures, surface level first. Returns spectral
    up and down fluxes (levels, nu).
    """
    nlay = tau.shape[0]
    if len(t_levels) != nlay + 1:
        raise ValueError('Need one more level than layers')
    b = planck(nu[None, :], np.asarray(t_levels, dtype=float)[:, None])
    b_surf = planck(nu, t_surface)
    mu, w = angles(n_angles, diffusivity)
    up = np.zeros((nlay + 1, nu.size))
    down = np.zeros((nlay + 1, nu.size))
    for m, wt in zip(mu, w):
        weight = 2 * math.pi * wt * m
        i_down = np.zeros(nu.size)
        down_levels = [i_down]
        terms = []
        for k in range(nlay - 1, -1, -1):
            e, a, c = _linear_source_terms(tau[k] / m)
            terms.append((e, a, c))
            # Downward ray enters at the top of layer k (level k+1), leaves at level k.
            i_down = i_down * e + b[k] * a - (b[k] - b[k + 1]) * c
            down_levels.append(i_down)
        terms.reverse()
        down_levels.reverse()
        i_up = emissivity * b_surf + (1 - emissivity) * down_levels[0]
        up[0] += weight * i_up
        for k in range(nlay):
            e, a, c = terms[k]
            i_up = i_up * e + b[k + 1] * a - (b[k + 1] - b[k]) * c
            up[k + 1] += weight * i_up
        for k in range(nlay + 1):
            down[k] += weight * down_levels[k]
    return up, down


def integrate(nu, spectral):
    """Integrate spectral fluxes over a uniform grid (trapezoid)."""
    step = nu[1] - nu[0]
    return step * (np.sum(spectral, axis=-1) - 0.5 * (spectral[..., 0] + spectral[..., -1]))
