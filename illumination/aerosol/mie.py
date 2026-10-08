"""Mie scattering by homogeneous spheres: extinction and scattering efficiencies and the asymmetry parameter.

The series follows Bohren and Huffman's BHMIE (Absorption and Scattering of Light by Small Particles, 1983, appendix
A): Riccati-Bessel functions by upward recurrence, the logarithmic derivative D_n(mx) by downward recurrence from
n_max + 15, and n_max = x + 4 x^(1/3) + 2 terms.
"""
from __future__ import annotations
import numpy as np


def efficiencies(x: float, m: complex) -> tuple[float, float, float]:
    """(Q_ext, Q_sca, g) for size parameter x = 2 pi r / lambda and relative refractive index m = n + i k."""
    if x <= 0:
        return 0.0, 0.0, 0.0
    nstop = int(x + 4.0 * x ** (1.0 / 3.0) + 2.0)
    mx = m * x
    nmx = int(max(nstop, abs(mx)) + 15)
    d = np.zeros(nmx + 1, dtype=complex)
    for n in range(nmx, 0, -1):
        d[n - 1] = n / mx - 1.0 / (d[n] + n / mx)
    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = -np.sin(x), np.cos(x)
    xi1 = complex(psi1, -chi1)
    qsca = qext = gsum = 0.0
    an_prev = bn_prev = None
    for n in range(1, nstop + 1):
        fn = (2.0 * n + 1.0) / (n * (n + 1.0))
        psi = (2.0 * n - 1.0) * psi1 / x - psi0
        chi = (2.0 * n - 1.0) * chi1 / x - chi0
        xi = complex(psi, -chi)
        an = ((d[n] / m + n / x) * psi - psi1) / ((d[n] / m + n / x) * xi - xi1)
        bn = ((m * d[n] + n / x) * psi - psi1) / ((m * d[n] + n / x) * xi - xi1)
        qsca += (2.0 * n + 1.0) * (abs(an) ** 2 + abs(bn) ** 2)
        qext += (2.0 * n + 1.0) * (an.real + bn.real)
        gsum += fn * (an * bn.conjugate()).real
        if an_prev is not None:
            gsum += (n - 1.0) * (n + 1.0) / n * (an_prev * an.conjugate() + bn_prev * bn.conjugate()).real
        an_prev, bn_prev = an, bn
        psi0, psi1, chi0, chi1 = psi1, psi, chi1, chi
        xi1 = complex(psi1, -chi1)
    qsca *= 2.0 / x ** 2
    qext *= 2.0 / x ** 2
    g = 4.0 * gsum / (qsca * x ** 2) if qsca > 0 else 0.0
    return float(qext), float(qsca), float(g)


def lognormal_mode(number_cm3: float, median_d_um: float, sigma: float, m: complex, wavelength_um: float,
                   bins: int = 120) -> dict:
    """Extinction and scattering coefficients (1/km) and asymmetry of a lognormal number distribution."""
    lnd = np.linspace(np.log(median_d_um) - 4.0 * np.log(sigma), np.log(median_d_um) + 4.0 * np.log(sigma), bins)
    dlnd = lnd[1] - lnd[0]
    diam = np.exp(lnd)
    weight = np.exp(-0.5 * ((lnd - np.log(median_d_um)) / np.log(sigma)) ** 2) / (np.sqrt(2 * np.pi) * np.log(sigma))
    n_m3 = number_cm3 * 1e6 * weight * dlnd
    area_m2 = np.pi * (diam * 1e-6 / 2.0) ** 2
    ext = sca = gsca = 0.0
    for nd, a, dd in zip(n_m3, area_m2, diam):
        qe, qs, g = efficiencies(np.pi * dd / wavelength_um, m)
        ext += nd * a * qe
        sca += nd * a * qs
        gsca += nd * a * qs * g
    return dict(extinction_per_km=ext * 1e3, scattering_per_km=sca * 1e3, asymmetry=gsca / sca if sca > 0 else 0.0)
