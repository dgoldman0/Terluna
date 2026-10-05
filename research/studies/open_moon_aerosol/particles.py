"""Particle populations for the aerosol budget: lognormal modes of dry particles, each with its hygroscopicity, and what
they do to the air's small ions and to clouds.

A mode is a number concentration (per m3), a dry median diameter (m), a geometric standard deviation and κ. Its mass
follows from the volume of the lognormal, its ion sink from Hõrrak et al.'s (2008) attachment coefficient integrated
over the grown sizes (atmosphere/electricity/conductivity.py), and its cloud condensation nuclei from κ-Köhler theory
(Petters and Kreidenweis 2007, Eq. 10): a particle activates at supersaturation s if its dry diameter exceeds
(4 A^3 / (27 κ ln^2(1 + s)))^(1/3).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

from atmosphere.electricity import conductivity as cd

DENSITY_KG_M3 = 1700.0                 # dry particle density for mass and number (a mixed organic, sulfate and dust value)
NODES = 64                             # points of the size integral, in units of the mode's log-width


@dataclass(frozen=True)
class Mode:
    number_m3: float
    d_m: float                         # dry median diameter
    sigma: float                       # geometric standard deviation
    kappa: float
    name: str = ''

    def mass_kg_m3(self, density: float = DENSITY_KG_M3) -> float:
        return self.number_m3 * density * math.pi / 6.0 * self.d_m ** 3 * math.exp(4.5 * math.log(self.sigma) ** 2)

    def scaled(self, factor: float) -> 'Mode':
        return replace(self, number_m3=self.number_m3 * factor)


def number_for_mass(mass_kg_m3: float, d_m: float, sigma: float, density: float = DENSITY_KG_M3) -> float:
    """The number (per m3) of a lognormal mode holding the given dry mass."""
    return mass_kg_m3 / (density * math.pi / 6.0 * d_m ** 3 * math.exp(4.5 * math.log(sigma) ** 2))


def _sizes(mode: Mode):
    """Quadrature nodes and weights over the mode's dry diameters (Gauss-Hermite in the log of the diameter)."""
    x, w = np.polynomial.hermite_e.hermegauss(NODES)
    d = mode.d_m * np.exp(x * math.log(mode.sigma))
    return d, w / w.sum()


def ion_sink(modes, p_pa: float, t_k: float, rh: float) -> float:
    """The small ions' attachment rate (1/s) to all modes, each grown at the relative humidity by its κ."""
    total = 0.0
    for m in modes:
        if m.number_m3 <= 0.0:
            continue
        d, w = _sizes(m)
        wet = d * cd.growth_factor(rh, m.kappa, d, t_k)
        total += m.number_m3 * float(np.sum(w * cd.attachment_coefficient(wet, p_pa, t_k)))
    return total


def critical_diameter(kappa: float, supersaturation: float, t_k: float) -> float:
    """The dry diameter (m) above which a particle of hygroscopicity κ activates at the supersaturation (a fraction)."""
    a = 4.0 * cd.SURFACE_TENSION * cd.WATER_MOLAR_MASS / (cd.GAS_CONSTANT * t_k * cd.WATER_DENSITY)
    return (4.0 * a ** 3 / (27.0 * kappa * math.log(1.0 + supersaturation) ** 2)) ** (1.0 / 3.0)


def ccn(modes, supersaturation: float, t_k: float) -> float:
    """Cloud condensation nuclei (per m3) at the supersaturation."""
    total = 0.0
    for m in modes:
        if m.number_m3 <= 0.0 or m.kappa <= 0.0:
            continue
        dc = critical_diameter(m.kappa, supersaturation, t_k)
        z = math.log(dc / m.d_m) / (math.sqrt(2.0) * math.log(m.sigma))
        total += m.number_m3 * 0.5 * math.erfc(z)
    return total


def conductivity(modes, q_m3_s: float, p_pa: float, t_k: float, rh: float) -> float:
    """The air's conductivity (S/m) with these particles at this humidity, for ion production q."""
    return float(cd.conductivity_of(q_m3_s, p_pa, t_k, ion_sink(modes, p_pa, t_k, rh))[0])


def night_layer(c0: float, flux: float, depth_m: float, loss_s: float, hours: float) -> float:
    """The mean over the night of a concentration in a surface layer of the given depth that starts at c0, takes a
    surface flux (per m2 per s) and loses at loss_s (1/s): C(t) = C0 e^(-kt) + F/(h k) (1 - e^(-kt))."""
    t = hours * 3600.0
    if loss_s <= 0.0:
        return c0 + flux / depth_m * t / 2.0
    steady = flux / (depth_m * loss_s)
    decay = (1.0 - math.exp(-loss_s * t)) / (loss_s * t)
    return steady + (c0 - steady) * decay
