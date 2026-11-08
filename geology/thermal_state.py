"""The present Moon's temperature with depth and the viscosity it gives an olivine mantle.

Three published present-day profiles span the Moon's thermal state (interior_sources.json): Kuskov & Kronrod (2009)
from seismic velocities (cool), Karato (2013) from electrical conductivity and tidal Q for dry olivine (central) and
Khan et al. (2014) from electromagnetic sounding, mass, moment of inertia and k2 (warm). They differ by about 400 K
at 200 km depth and by about 200 K at 1,000 km. The near side's mantle is 100-200 K warmer than the far side's (Park
et al. 2025); the near and far profiles add and subtract half that difference.

Laboratory flow laws for olivine (Hirth & Kohlstedt 2003) turn temperature and pressure into viscosity: diffusion
creep, which depends on grain size, and dislocation creep, which depends on stress, acting in parallel. Grain size,
water content and the creep stress are assumptions that come with the laws, stated in the sources file; the laws
themselves are measured on Earth's olivine. Uniaxial laboratory rates convert to viscosity as sigma / (3 strain rate).
The low-viscosity zone above the core keeps its published viscosity, and the cold crust is elastic. The rigid outer
layer follows from the profile as the depth above which the Maxwell time exceeds the time asked about.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import AVOGADRO, BOLTZMANN, GRAVITATIONAL_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS

HERE = Path(__file__).resolve().parent
SOURCES = HERE / "interior_sources.json"
SPEC = json.loads(SOURCES.read_text())
GAS_CONSTANT = BOLTZMANN * AVOGADRO
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
KELVIN = 273.15
LVZ_TOP_M = 500.0e3
LVZ_VISCOSITY_PA_S = (2.0e16, 4.5e16)          # Harada et al. 2014; Briaud et al. 2023
GRAIN_SIZES_M = (0.001, 0.003, 0.01)          # assumption: 3 mm central
STRESSES_MPA = (0.3, 1.0, 3.0)                # assumption: 1 MPa central
WET_C_OH = 100.0 * 16.0                       # 0.01 wt% H2O (Karato 2013's hydrous case), H/1e6 Si
NEAR_FAR_OFFSET_K = (50.0, 75.0, 100.0)       # half of Park et al.'s 100-200 K near-far difference
_LAWS = next(s for s in SPEC["sources"] if s["key"] == "hirth_kohlstedt2003")["values_used"]


@dataclass(frozen=True)
class Rheology:
    """The assumptions that fix a viscosity profile."""
    profile: str = "central"          # 'cool', 'central' or 'warm'
    offset_k: float = 0.0             # added to the profile (near side positive)
    grain_m: float = 0.003
    stress_mpa: float = 1.0
    wet: bool = False
    lvz_pa_s: float = LVZ_VISCOSITY_PA_S[0]

    def label(self):
        side = "" if self.offset_k == 0 else (f", {self.offset_k:+.0f} K")
        return (f"{self.profile}{side}, {'wet' if self.wet else 'dry'}, {1e3 * self.grain_m:g} mm, "
                f"{self.stress_mpa:g} MPa")


CENTRAL = Rheology()
SOFT = Rheology(profile="warm", grain_m=0.001, stress_mpa=3.0, wet=True)
STIFF = Rheology(profile="cool", grain_m=0.01, stress_mpa=0.3, wet=False)


def temperature_k(depth_m, profile="central", offset_k=0.0):
    """Temperature (K) at depth, linear between the digitised nodes and extended at the end gradients."""
    spec = SPEC["temperature_profiles"][profile]
    z = np.asarray(spec["depth_km"], float) * 1e3
    t = np.asarray(spec["temperature_c"], float) + KELVIN
    depth = np.asarray(depth_m, float)
    inside = np.interp(depth, z, t)
    above = t[0] + (depth - z[0]) * (t[1] - t[0]) / (z[1] - z[0])
    below = t[-1] + (depth - z[-1]) * (t[-1] - t[-2]) / (z[-1] - z[-2])
    return np.where(depth < z[0], above, np.where(depth > z[-1], below, inside)) + offset_k


def layers():
    """(top radius m, density, rigidity) of the interior's homogeneous layers, centre outward."""
    return [(float(row["top_km"]) * 1e3, float(row["density_kg_m3"]), float(row["rigidity_pa"]))
            for row in SPEC["interior_model"]["layers"]]


def pressure_pa(depth_m, samples=4000):
    """Hydrostatic pressure of the layered interior at depth."""
    stack = layers()
    r = np.linspace(MOON_RADIUS, 1.0, samples)
    rho = np.select([r <= top for top, _, _ in stack], [d for _, d, _ in stack])
    shell_mass = 4 * math.pi / 3 * np.array([d * (top**3 - (stack[i - 1][0] if i else 0.0)**3)
                                             for i, (top, d, _) in enumerate(stack)])
    below = np.r_[0.0, np.cumsum(shell_mass)[:-1]]
    k = np.searchsorted([top for top, _, _ in stack], r)
    bottoms = np.array([stack[i - 1][0] if i else 0.0 for i in range(len(stack))])
    mass = below[k] + 4 * math.pi / 3 * rho * (r**3 - bottoms[k]**3)
    g = GRAVITATIONAL_CONSTANT * mass / r**2
    p = np.r_[0.0, np.cumsum(0.5 * (rho[1:] * g[1:] + rho[:-1] * g[:-1]) * -np.diff(r))]
    return np.interp(np.asarray(depth_m, float), MOON_RADIUS - r, p)


def creep_viscosity(temperature, pressure, grain_m, stress_mpa, wet):
    """Diffusion, dislocation and composite viscosity (Pa s) of olivine; log-safe for cold rock."""
    out = {}
    for name in ("diffusion", "dislocation"):
        law = _LAWS[name]
        a = law["a_wet"] * WET_C_OH ** law["r"] if wet else law["a_dry"]
        log_rate = (math.log(a) + (law["n"] - 1) * math.log(stress_mpa) - law["p"] * math.log(grain_m * 1e6)
                    - (law["e_j_mol"] + pressure * law["v_m3_mol"]) / (GAS_CONSTANT * temperature))
        # strain rate per MPa of stress, in 1/s; viscosity = sigma / (3 rate) in Pa s
        out[name] = np.exp(np.log(1e6 / 3.0) - log_rate)
    out["composite"] = 1.0 / (1.0 / out["diffusion"] + 1.0 / out["dislocation"])
    return out


def viscosity_pa_s(depth_m, rheology=CENTRAL):
    """Effective viscosity (Pa s) at depth: composite creep in the mantle, the published low-viscosity zone below
    it and an elastic crust (infinite)."""
    depth = np.asarray(depth_m, float)
    t = temperature_k(depth, rheology.profile, rheology.offset_k)
    eta = creep_viscosity(t, pressure_pa(depth), rheology.grain_m, rheology.stress_mpa, rheology.wet)["composite"]
    crust_base = MOON_RADIUS - layers()[-2][0]
    eta = np.where(depth < crust_base, np.inf, eta)
    return np.where(depth > MOON_RADIUS - LVZ_TOP_M, rheology.lvz_pa_s, eta)


def lid_depth_m(years, rheology=CENTRAL):
    """Depth above which the mantle's Maxwell time (viscosity / rigidity) exceeds `years`: the layer that stays
    rigid over that time."""
    mantle_mu = layers()[-2][2]
    z = np.linspace(MOON_RADIUS - layers()[-2][0], MOON_RADIUS - LVZ_TOP_M, 2000)
    tau = viscosity_pa_s(z, rheology) / mantle_mu / YEAR_S
    weak = np.flatnonzero(tau < years)
    return float(z[weak[0]]) if len(weak) else float(z[-1])


def shells(rheology=CENTRAL, step_m=10.0e3):
    """The interior as homogeneous Maxwell shells for viscoelastic.py: (top m, density, rigidity, viscosity),
    the mantle cut into sublayers of about `step_m`, each at its mid-depth viscosity."""
    core, lvz, mantle, crust = layers()
    out = [(core[0], core[1], 0.0, 0.0), (lvz[0], lvz[1], lvz[2], rheology.lvz_pa_s)]
    count = max(1, int(round((mantle[0] - lvz[0]) / step_m)))
    edges = np.linspace(lvz[0], mantle[0], count + 1)
    mids = MOON_RADIUS - 0.5 * (edges[1:] + edges[:-1])
    for top, eta in zip(edges[1:], viscosity_pa_s(mids, rheology)):
        out.append((float(top), mantle[1], mantle[2], float(eta)))
    out.append((crust[0], crust[1], crust[2], math.inf))
    return out


def near_far(rheology=CENTRAL, offset_k=NEAR_FAR_OFFSET_K[1]):
    """The near side's and far side's rheology: the profile plus and minus half the near-far difference."""
    return replace(rheology, offset_k=offset_k), replace(rheology, offset_k=-offset_k)
