"""C3 leaf photosynthesis (Farquhar, von Caemmerer and Berry 1980; von Caemmerer 2000).

Gross assimilation is the lesser of the Rubisco-limited rate Vcmax (ci - G*)/(ci + Kc(1 + O/Ko))
and the electron-transport-limited rate J (ci - G*)/(4 ci + 8 G*). J follows the absorbed
photosynthetic photons I through theta J^2 - (I2 + Jmax) J + I2 Jmax = 0, with I2 = I (1 - f)/2
and f = 0.15. Rubisco kinetics are Bernacchi et al. (2001), measured near 100 kPa and applied
as partial pressures; capacities follow Arrhenius functions with the mean activation energies
of Medlyn et al. (2002), and day respiration Bernacchi's. The 25 C capacities are those of
research/studies/atmospheric_co2, whose light-saturated leaf this reproduces.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

R_KJ = 0.008314
MEASUREMENT_PA = 100.0e3
KINETICS = {'gamma_star': (19.02, 37.83, 1e-6), 'kc': (38.05, 79.43, 1e-6), 'ko': (20.30, 36.38, 1e-3)}
EA_VCMAX, EA_JMAX, EA_RD = 65.0, 50.0, 46.39       # kJ/mol


@dataclass(frozen=True)
class Traits:
    vcmax25: float = 60.0            # umol m-2 s-1, at the top of the canopy
    jmax_ratio: float = 100.0 / 60.0
    rd_fraction: float = 0.015       # day respiration as a share of Vcmax
    curvature: float = 0.7           # theta
    unused_share: float = 0.15       # f: absorbed photons not reaching the reaction centres
    ci_ca: float = 0.7               # internal to ambient CO2, a typical C3 value


def kinetics(t_c):
    """(Gamma*, Kc, Ko) in Pa at leaf temperature t_c (C)."""
    t = t_c + 273.15
    return tuple(np.exp(c - dh / (R_KJ * t)) * unit * MEASUREMENT_PA for c, dh, unit in KINETICS.values())


def arrhenius(ea, t_c):
    return np.exp(ea / R_KJ * (1 / 298.15 - 1 / (t_c + 273.15)))


def gross(absorbed, vcmax25, t_c, ca_pa, o2_pa, traits: Traits = Traits()):
    """Gross assimilation (umol CO2 m-2 s-1) for absorbed photosynthetic photons (umol m-2 s-1)."""
    gamma, kc, ko = kinetics(t_c)
    ci = traits.ci_ca * ca_pa
    vcmax = vcmax25 * arrhenius(EA_VCMAX, t_c)
    jmax = traits.jmax_ratio * vcmax25 * arrhenius(EA_JMAX, t_c)
    i2 = np.asarray(absorbed, dtype=float) * (1 - traits.unused_share) / 2
    s = i2 + jmax
    j = (s - np.sqrt(np.maximum(s * s - 4 * traits.curvature * i2 * jmax, 0.0))) / (2 * traits.curvature)
    rubisco = vcmax * (ci - gamma) / (ci + kc * (1 + o2_pa / ko))
    light = j * (ci - gamma) / (4 * ci + 8 * gamma)
    return np.minimum(rubisco, light)


def respiration(vcmax25, t_c, traits: Traits = Traits()):
    """Leaf respiration (umol CO2 m-2 s-1)."""
    return traits.rd_fraction * vcmax25 * arrhenius(EA_RD, t_c)
