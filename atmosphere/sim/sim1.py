"""
Lunar Atmosphere Retention Simulator — v0.5

Adds the two most important upgrades:

1. Upper-atmosphere energy balance
   - Solves upper-atmosphere temperature from a simplified heat budget.
   - Includes residual high-energy solar heating, radiative cooling, optional
     lower-atmosphere thermal coupling, and escape cooling.

2. Transitional outflow / non-Jeans escape check
   - Computes Jeans escape and an isothermal Parker-wind-style outflow estimate.
   - Blends between Jeans and Parker outflow with a tunable coupling parameter.
   - This brackets the uncertain regime where N2 Jeans lambda is below ~10.

Still included:
- Spherical lunar hydrostatic structure.
- Diffusive separation above a homopause.
- Species-specific exobases for N2, O2, and Ar.
- Parametric H/N/O photochemical escape budget.

This remains a reduced screening model, not a production aeronomy code.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import exp, isfinite, log, pi, sqrt
from typing import Dict, Iterable, Tuple

import numpy as np

K_B = 1.380649e-23
AMU = 1.66053906660e-27
G = 6.67430e-11
YEAR = 365.25 * 24 * 3600

R_MOON = 1_737_400.0
M_MOON = 7.342e22
MU_MOON = G * M_MOON
G_SURF_MOON = MU_MOON / R_MOON**2
M_H = 1.00784 * AMU


@dataclass(frozen=True)
class Species:
    name: str
    molar_mass_amu: float
    mole_fraction: float
    collision_diameter_m: float

    @property
    def molecular_mass_kg(self) -> float:
        return self.molar_mass_amu * AMU

    @property
    def collision_cross_section_m2(self) -> float:
        return pi * self.collision_diameter_m**2


DEFAULT_BULK_SPECIES = (
    Species("N2", 28.0134, 0.78084, 3.70e-10),
    Species("O2", 31.9988, 0.20946, 3.46e-10),
    Species("Ar", 39.9480, 0.00934, 3.40e-10),
)


@dataclass
class EnergyBalanceParams:
    """
    Simplified global-mean thermosphere energy budget.

    high_energy_flux_w_m2:
        Unshielded thermosphere-relevant EUV/FUV/soft-X-ray flux normal to the
        Sun. The model converts it to global mean using global_avg_factor.

    residual_heating_fraction:
        Stored in AtmosphereScenario. 1e-3 means 0.1% effective heating leakage.

    lower_coupling_coeff_w_m2_k:
        Effective conductive/eddy coupling to the lower atmosphere. This is a
        major uncertainty. Stronger coupling can warm a cold thermosphere but can
        also power more escape.

    radiative_floor_k:
        Temperature where the simplified incremental radiative cooling term goes
        to zero. This acts as a cold IR-equilibrium floor.
    """
    high_energy_flux_w_m2: float = 0.006
    global_avg_factor: float = 0.25
    heating_efficiency: float = 0.3
    unshielded_upper_temperature_k: float = 1000.0
    radiative_floor_k: float = 240.0
    cooling_exponent: float = 2.0
    lower_coupling_coeff_w_m2_k: float = 0.0
    lower_reservoir_temperature_k: float = 288.0


@dataclass
class OutflowParams:
    """
    Transitional escape model.

    coupling = 0:
        Pure Jeans escape.
    coupling = 1:
        Full isothermal Parker-wind-style flux where that exceeds Jeans.
    coupling = 0.03-0.3:
        Transitional cases. This is the main uncertainty once N2 lambda is below
        ~10.
    """
    enabled: bool = True
    coupling: float = 0.10
    transition_lambda: float = 15.0
    specific_energy_factor: float = 1.0


@dataclass
class PhotochemistryParams:
    enabled: bool = True
    water_vapor_ppm: float = 10.0
    h_source_saturation_leakage: float = 1.0e-3
    h_diffusion_limit_coeff_atoms_m2_s: float = 2.5e17
    n_hot_escape_unshielded_ref_kg_s: float = 3_000.0
    o_hot_escape_unshielded_ref_kg_s: float = 1_500.0


@dataclass
class AtmosphereScenario:
    name: str = "baseline"
    surface_pressure_pa: float = 1.2 * 101_325.0
    surface_temperature_k: float = 288.0
    residual_heating_fraction: float = 1.0e-3

    transition_altitude_m: float = 120_000.0
    transition_width_m: float = 60_000.0
    homopause_altitude_m: float = 100_000.0

    species: Tuple[Species, ...] = DEFAULT_BULK_SPECIES
    energy: EnergyBalanceParams = field(default_factory=EnergyBalanceParams)
    outflow: OutflowParams = field(default_factory=OutflowParams)
    photochemistry: PhotochemistryParams = field(default_factory=PhotochemistryParams)

    z_max_m: float = 120_000_000.0
    dz_m: float = 20_000.0


@dataclass
class SpeciesResult:
    name: str
    exobase_altitude_m: float
    exobase_temperature_k: float
    exobase_number_density_m3: float
    exobase_escape_speed_m_s: float
    jeans_lambda: float
    jeans_escape_kg_s: float
    parker_escape_kg_s: float
    hybrid_escape_kg_s: float
    exobase_found_inside_grid: bool


@dataclass
class PhotochemistryResult:
    h_escape_kg_s: float
    n_hot_escape_kg_s: float
    o_hot_escape_kg_s: float
    total_photochemical_escape_kg_s: float
    water_equivalent_loss_kg_s: float


@dataclass
class EnergyBalanceResult:
    upper_temperature_k: float
    q_high_energy_w_m2: float
    q_lower_coupling_w_m2: float
    q_radiative_cooling_w_m2: float
    q_escape_cooling_w_m2: float
    net_flux_w_m2: float
    floor_limited: bool


@dataclass
class SimulationResult:
    scenario: AtmosphereScenario
    total_atmosphere_mass_kg: float
    upper_temperature_k: float
    energy_balance: EnergyBalanceResult
    escape_kg_s: float
    photochemical_escape_kg_s: float
    total_escape_kg_s: float
    total_cycling_time_years: float
    species_results: Dict[str, SpeciesResult]
    photochemistry_result: PhotochemistryResult
    warnings: Tuple[str, ...]


def gravity_at_radius(r_m):
    return MU_MOON / np.asarray(r_m) ** 2


def escape_speed_at_radius(r_m):
    return np.sqrt(2.0 * MU_MOON / np.asarray(r_m))


def mean_molecular_mass_kg(species: Iterable[Species]) -> float:
    total_x = sum(s.mole_fraction for s in species)
    return sum(s.mole_fraction * s.molecular_mass_kg for s in species) / total_x


def total_atmosphere_mass(surface_pressure_pa: float) -> float:
    return surface_pressure_pa * 4.0 * pi * R_MOON**2 / G_SURF_MOON


def cumulative_integral(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def temperature_profile(z_m: np.ndarray, upper_temperature_k: float, scenario: AtmosphereScenario) -> np.ndarray:
    x = (z_m - scenario.transition_altitude_m) / scenario.transition_width_m
    blend = 0.5 * (1.0 + np.tanh(x))
    return scenario.surface_temperature_k * (1.0 - blend) + upper_temperature_k * blend


def build_bulk_profiles(upper_temperature_k: float, scenario: AtmosphereScenario):
    z = np.arange(0.0, scenario.z_max_m + scenario.dz_m, scenario.dz_m)
    T = temperature_profile(z, upper_temperature_k, scenario)
    r = R_MOON + z
    mbar = mean_molecular_mass_kg(scenario.species)
    integrand = mbar * gravity_at_radius(r) / (K_B * T)
    cumulative = cumulative_integral(integrand, z)
    P = scenario.surface_pressure_pa * np.exp(-cumulative)
    n = P / (K_B * T)
    return z, T, P, n


def build_species_number_density_profile(z_m, T_k, n_bulk_m3, scenario, species):
    total_x = sum(s.mole_fraction for s in scenario.species)
    x_i = species.mole_fraction / total_x
    n_i = np.zeros_like(z_m)

    h_idx = int(np.searchsorted(z_m, scenario.homopause_altitude_m))
    h_idx = min(max(h_idx, 0), len(z_m) - 1)
    n_i[: h_idx + 1] = x_i * n_bulk_m3[: h_idx + 1]

    if h_idx < len(z_m) - 1:
        r = R_MOON + z_m[h_idx:]
        integrand = species.molecular_mass_kg * gravity_at_radius(r) / (K_B * T_k[h_idx:])
        cumulative = cumulative_integral(integrand, z_m[h_idx:])
        n_i[h_idx:] = n_i[h_idx] * np.exp(-cumulative)
    return n_i


def local_density_scale_height_m(z_m, T_k, species_mass_kg):
    return K_B * T_k / (species_mass_kg * gravity_at_radius(R_MOON + z_m))


def mean_free_path_m(n_m3, cross_section_m2):
    return 1.0 / (sqrt(2.0) * np.maximum(n_m3, 1e-300) * cross_section_m2)


def find_exobase(z_m, T_k, n_species_m3, species):
    H = local_density_scale_height_m(z_m, T_k, species.molecular_mass_kg)
    mfp = mean_free_path_m(n_species_m3, species.collision_cross_section_m2)
    ratio = mfp / H
    crossing_idxs = np.where(ratio >= 1.0)[0]
    found = len(crossing_idxs) > 0
    if not found:
        idx = len(z_m) - 1
        return float(z_m[idx]), float(T_k[idx]), float(n_species_m3[idx]), False
    idx = int(crossing_idxs[0])
    if idx == 0:
        return float(z_m[0]), float(T_k[0]), float(n_species_m3[0]), True

    z0, z1 = z_m[idx - 1], z_m[idx]
    y0, y1 = np.log(ratio[idx - 1]), np.log(ratio[idx])
    frac = 0.0 if y1 == y0 else (0.0 - y0) / (y1 - y0)
    z_exo = z0 + frac * (z1 - z0)
    T_exo = np.interp(z_exo, z_m, T_k)
    n_exo = np.exp(np.interp(z_exo, z_m, np.log(np.maximum(n_species_m3, 1e-300))))
    return float(z_exo), float(T_exo), float(n_exo), True


def jeans_particle_flux_m2_s(n_exo_m3, T_exo_k, species_mass_kg, r_exo_m):
    v_th = sqrt(2.0 * K_B * T_exo_k / species_mass_kg)
    lam = MU_MOON * species_mass_kg / (r_exo_m * K_B * T_exo_k)
    flux = n_exo_m3 * v_th / (2.0 * sqrt(pi)) * (1.0 + lam) * exp(-lam)
    return flux, lam


def parker_particle_flux_m2_s(n_exo_m3, T_exo_k, species_mass_kg, r_exo_m):
    """
    Isothermal Parker-wind-style estimate at the exobase.

    This is not assumed to be fully realized. The model blends toward this value
    with OutflowParams.coupling.
    """
    c_s = sqrt(K_B * T_exo_k / species_mass_kg)
    lam = MU_MOON * species_mass_kg / (r_exo_m * K_B * T_exo_k)
    if lam <= 2.0:
        u = c_s
    else:
        # Parker equation on the subsonic branch:
        # y - ln(y) = 4 ln(r/r_s) + 4 r_s/r - 3, y=(u/c_s)^2, r_s/r=lambda/2.
        A = 4.0 * log(2.0 / lam) + 2.0 * lam - 3.0
        if A <= 1.0:
            y = 1.0
        elif A > 50.0:
            y = exp(-A)
        else:
            lo, hi = 1e-300, 1.0
            for _ in range(80):
                mid = sqrt(lo * hi)
                val = mid - log(mid)
                if val > A:
                    lo = mid
                else:
                    hi = mid
            y = sqrt(lo * hi)
        u = c_s * sqrt(y)
    return n_exo_m3 * u, lam


def escape_for_upper_temperature(upper_temperature_k: float, scenario: AtmosphereScenario):
    z, T, _P, n_bulk = build_bulk_profiles(upper_temperature_k, scenario)
    total_escape = 0.0
    details: Dict[str, SpeciesResult] = {}

    for sp in scenario.species:
        n_sp = build_species_number_density_profile(z, T, n_bulk, scenario, sp)
        z_exo, T_exo, n_exo, found = find_exobase(z, T, n_sp, sp)
        r_exo = R_MOON + z_exo
        area = 4.0 * pi * r_exo**2

        jeans_flux, lam = jeans_particle_flux_m2_s(n_exo, T_exo, sp.molecular_mass_kg, r_exo)
        parker_flux, _ = parker_particle_flux_m2_s(n_exo, T_exo, sp.molecular_mass_kg, r_exo)

        jeans_kg_s = jeans_flux * area * sp.molecular_mass_kg
        parker_kg_s = parker_flux * area * sp.molecular_mass_kg
        hybrid_kg_s = jeans_kg_s

        if scenario.outflow.enabled and parker_kg_s > jeans_kg_s:
            # Smoothly blend below the transition lambda.
            weight = scenario.outflow.coupling / (1.0 + exp((lam - scenario.outflow.transition_lambda) / 2.0))
            hybrid_kg_s = jeans_kg_s + weight * (parker_kg_s - jeans_kg_s)

        total_escape += hybrid_kg_s
        details[sp.name] = SpeciesResult(
            name=sp.name,
            exobase_altitude_m=z_exo,
            exobase_temperature_k=T_exo,
            exobase_number_density_m3=n_exo,
            exobase_escape_speed_m_s=float(escape_speed_at_radius(r_exo)),
            jeans_lambda=lam,
            jeans_escape_kg_s=jeans_kg_s,
            parker_escape_kg_s=parker_kg_s,
            hybrid_escape_kg_s=hybrid_kg_s,
            exobase_found_inside_grid=found,
        )
    return total_escape, details


def estimate_photochemical_escape(scenario: AtmosphereScenario) -> PhotochemistryResult:
    p = scenario.photochemistry
    if not p.enabled:
        return PhotochemistryResult(0.0, 0.0, 0.0, 0.0, 0.0)

    fuv = max(scenario.residual_heating_fraction, 0.0)
    source_availability = min(1.0, fuv / p.h_source_saturation_leakage) if p.h_source_saturation_leakage > 0 else 1.0
    f_h_atoms = 2.0 * p.water_vapor_ppm * 1.0e-6
    area_homopause = 4.0 * pi * (R_MOON + scenario.homopause_altitude_m) ** 2
    h_atom_flux = p.h_diffusion_limit_coeff_atoms_m2_s * f_h_atoms * source_availability
    h_escape = h_atom_flux * area_homopause * M_H

    n_hot = p.n_hot_escape_unshielded_ref_kg_s * fuv
    o_hot = p.o_hot_escape_unshielded_ref_kg_s * fuv
    total = h_escape + n_hot + o_hot
    return PhotochemistryResult(h_escape, n_hot, o_hot, total, h_escape * 9.0)


def energy_fluxes_for_temperature(upper_temperature_k: float, scenario: AtmosphereScenario):
    e = scenario.energy
    q_high = (
        e.high_energy_flux_w_m2
        * e.global_avg_factor
        * e.heating_efficiency
        * scenario.residual_heating_fraction
    )
    q_high_unshielded = e.high_energy_flux_w_m2 * e.global_avg_factor * e.heating_efficiency
    q_cond_hot = e.lower_coupling_coeff_w_m2_k * (
        e.lower_reservoir_temperature_k - e.unshielded_upper_temperature_k
    )

    # Calibrate radiative coefficient so the unshielded state can sit near T_hot
    # in the absence of escape cooling.
    beta = e.cooling_exponent
    denom = max(1e-12, e.unshielded_upper_temperature_k**beta - e.radiative_floor_k**beta)
    q_rad_ref = max(1e-12, q_high_unshielded + q_cond_hot)
    cooling_coeff = q_rad_ref / denom
    q_rad = cooling_coeff * max(0.0, upper_temperature_k**beta - e.radiative_floor_k**beta)
    q_cond = e.lower_coupling_coeff_w_m2_k * (e.lower_reservoir_temperature_k - upper_temperature_k)

    escape_kg_s, details = escape_for_upper_temperature(upper_temperature_k, scenario)
    area_surface = 4.0 * pi * R_MOON**2
    escape_power = 0.0
    for sp in scenario.species:
        d = details[sp.name]
        r_exo = R_MOON + d.exobase_altitude_m
        specific_energy_j_kg = MU_MOON / r_exo + 2.5 * K_B * d.exobase_temperature_k / sp.molecular_mass_kg
        escape_power += d.hybrid_escape_kg_s * specific_energy_j_kg * scenario.outflow.specific_energy_factor
    q_escape = escape_power / area_surface

    net = q_high + q_cond - q_rad - q_escape
    return q_high, q_cond, q_rad, q_escape, net


def solve_upper_temperature(scenario: AtmosphereScenario) -> EnergyBalanceResult:
    e = scenario.energy
    lo = max(80.0, e.radiative_floor_k)
    hi = e.unshielded_upper_temperature_k * 1.5
    qh, qc, qr, qe, flo = energy_fluxes_for_temperature(lo, scenario)
    qh_hi, qc_hi, qr_hi, qe_hi, fhi = energy_fluxes_for_temperature(hi, scenario)

    floor_limited = False
    if flo <= 0.0:
        # The prescribed floor is warmer than the formal energy equilibrium.
        # Interpret as floor-limited: the escape rate at this temperature is an
        # upper-bound unless lower-atmosphere coupling supplies the difference.
        floor_limited = True
        return EnergyBalanceResult(lo, qh, qc, qr, qe, flo, floor_limited)

    if fhi >= 0.0:
        # No stable upper bracket. This is a likely runaway/outflow regime.
        return EnergyBalanceResult(hi, qh_hi, qc_hi, qr_hi, qe_hi, fhi, False)

    for _ in range(70):
        mid = 0.5 * (lo + hi)
        qh_m, qc_m, qr_m, qe_m, fm = energy_fluxes_for_temperature(mid, scenario)
        if fm > 0.0:
            lo = mid
        else:
            hi = mid
    T = 0.5 * (lo + hi)
    qh, qc, qr, qe, net = energy_fluxes_for_temperature(T, scenario)
    return EnergyBalanceResult(T, qh, qc, qr, qe, net, False)


def run_scenario(scenario: AtmosphereScenario) -> SimulationResult:
    energy = solve_upper_temperature(scenario)
    escape_kg_s, species_results = escape_for_upper_temperature(energy.upper_temperature_k, scenario)
    photo = estimate_photochemical_escape(scenario)
    total_escape = escape_kg_s + photo.total_photochemical_escape_kg_s
    mass = total_atmosphere_mass(scenario.surface_pressure_pa)
    cycle = mass / total_escape / YEAR if total_escape > 0 else float("inf")

    warnings = []
    if energy.floor_limited:
        warnings.append("Energy solution is floor-limited; escape/outflow at the imposed floor may be an upper-bound or requires extra lower-atmosphere heat transport.")
    for r in species_results.values():
        if not r.exobase_found_inside_grid:
            warnings.append(f"{r.name}: exobase not found inside grid.")
        if r.jeans_lambda < 3.0:
            warnings.append(f"{r.name}: lambda < 3; hydrodynamic/blowoff regime likely.")
        elif r.name == "N2" and r.jeans_lambda < 10.0:
            warnings.append("N2 lambda < 10; transitional non-Jeans escape is important.")

    return SimulationResult(
        scenario=scenario,
        total_atmosphere_mass_kg=mass,
        upper_temperature_k=energy.upper_temperature_k,
        energy_balance=energy,
        escape_kg_s=escape_kg_s,
        photochemical_escape_kg_s=photo.total_photochemical_escape_kg_s,
        total_escape_kg_s=total_escape,
        total_cycling_time_years=cycle,
        species_results=species_results,
        photochemistry_result=photo,
        warnings=tuple(sorted(set(warnings))),
    )


def format_years(years: float) -> str:
    if not isfinite(years):
        return "infinite"
    if years >= 1e9:
        return f"{years / 1e9:.2g} Gyr"
    if years >= 1e6:
        return f"{years / 1e6:.2g} Myr"
    if years >= 1e3:
        return f"{years / 1e3:.2g} kyr"
    return f"{years:.2g} yr"


def make_case(floor=240.0, leakage=1e-3, outflow=0.10, coupling=0.0, name=None):
    return AtmosphereScenario(
        name=name or f"floor={floor}, leakage={leakage}, outflow={outflow}, coupling={coupling}",
        residual_heating_fraction=leakage,
        energy=EnergyBalanceParams(
            radiative_floor_k=floor,
            lower_coupling_coeff_w_m2_k=coupling,
        ),
        outflow=OutflowParams(coupling=outflow),
    )


def print_summary_table(cases):
    print("case | T_upper | escape kg/s | photo kg/s | total kg/s | cycle | N2 lambda | N2 z km")
    print("-" * 108)
    for c in cases:
        r = run_scenario(c)
        n2 = r.species_results["N2"]
        print(
            f"{c.name[:34]:34s} | {r.upper_temperature_k:7.1f} | {r.escape_kg_s:11.3g} | "
            f"{r.photochemical_escape_kg_s:10.3g} | {r.total_escape_kg_s:10.3g} | "
            f"{format_years(r.total_cycling_time_years):>8} | {n2.jeans_lambda:9.2f} | {n2.exobase_altitude_m/1000:7.0f}"
        )


def run_default_suite():
    cases = []
    for coupling in [0.0, 0.03, 0.10, 0.30, 1.0]:
        cases.append(make_case(floor=250, leakage=1e-3, outflow=coupling, coupling=0.0, name=f"250K floor, outflow={coupling}"))
    for floor in [230, 240, 250]:
        cases.append(make_case(floor=floor, leakage=1e-3, outflow=0.10, coupling=0.0, name=f"floor {floor}K, no lower coupling"))
    for floor in [230, 240, 250]:
        cases.append(make_case(floor=floor, leakage=1e-3, outflow=0.10, coupling=1e-7, name=f"floor {floor}K, lower coupling"))
    for leakage in [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]:
        cases.append(make_case(floor=240, leakage=leakage, outflow=0.10, coupling=1e-7, name=f"240K floor, leak={leakage:.0e}"))
    print_summary_table(cases)


if __name__ == "__main__":
    run_default_suite()
