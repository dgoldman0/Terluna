"""Measured-film hydrogen containment screens, not a living-envelope design.

Permeability is mol m/(m² s Pa); its driving pressure is the species partial
pressure difference, not the mechanical gauge pressure. See envelope.md and
envelope_sources.json for the narrow evidence attached to each comparator.
"""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from shared.constants import AVOGADRO, ELEMENTARY_CHARGE, GAS_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS

# Conventional permeability unit, rounded; chemical/reference properties are
# declared inputs, not replacements for the shared physical constants.
BARRER_SI = 3.35e-16
HYDROGEN_MOLAR_MASS_KG = 0.002016
WATER_MOLAR_MASS_KG = 0.01801528
OXYGEN_MOLAR_MASS_KG = 0.031998
WATER_SPLITTING_GIBBS_J_MOL = 237_141.0  # JANAF liquid water, 298.15 K, 1 bar
WATER_SPLITTING_ENTHALPY_J_MOL = 285_830.0


def _number(name: str, value: float, *, positive: bool = False) -> None:
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}")


def _fraction(name: str, value: float, *, positive: bool = False) -> None:
    _number(name, value, positive=positive)
    if value > 1:
        raise ValueError(f"{name} must not exceed one")


def permeation_flux(
    permeability_mol_m_m2_s_pa: float,
    thickness_m: float,
    inside_partial_pressure_pa: float,
    outside_partial_pressure_pa: float = 0.0,
) -> float:
    """Signed outward species flux, mol/(m² s), through an intact uniform film.

    Constant permeability, planar thin-film approximation, no seam/pore leaks.
    A positive total internal pressure does not prevent another species entering.
    """
    _number("permeability", permeability_mol_m_m2_s_pa)
    _number("thickness_m", thickness_m, positive=True)
    _number("inside_partial_pressure_pa", inside_partial_pressure_pa)
    _number("outside_partial_pressure_pa", outside_partial_pressure_pa)
    return permeability_mol_m_m2_s_pa / thickness_m * (
        inside_partial_pressure_pa - outside_partial_pressure_pa
    )


def layered_permeance(layers: Sequence[tuple[float, float]]) -> float:
    """Ideal series permeance mol/(m² s Pa); each tuple is (thickness_m, P_SI).

    Layers share area, have perfect bonding, and retain their measured properties.
    Real cracks, edge bypasses and hydration changes invalidate this ideal bound.
    """
    if not layers:
        raise ValueError("at least one layer is required")
    resistances = []
    for thickness, permeability in layers:
        _number("thickness_m", thickness, positive=True)
        _number("permeability", permeability)
        resistances.append(math.inf if permeability == 0 else thickness / permeability)
    return 1.0 / sum(resistances)


def crossover_current_to_hydrogen_flux(current_ma_cm2: float) -> float:
    """Convert measured H₂ oxidation crossover current to mol/(m² s), using 2F.

    Assumes the reported current has been attributed to hydrogen crossover.
    It is not an intrinsic permeability without the test pressure and thickness.
    """
    _number("current_ma_cm2", current_ma_cm2)
    faraday_c_mol = AVOGADRO * ELEMENTARY_CHARGE
    return current_ma_cm2 * 10.0 / (2.0 * faraday_c_mol)


def measured_flux_budget(
    *, hydrogen_flux_mol_m2_s: float, gas_wetted_area_m2: float,
    projected_area_m2: float, gibbs_efficiency: float = 1.0,
    gibbs_j_mol: float = WATER_SPLITTING_GIBBS_J_MOL,
) -> dict[str, float]:
    """Replacement requirements for an outward measured/predicted H₂ flux.

    Efficiency means reference Gibbs energy divided by supplied useful work.
    Do not divide by this efficiency again in a coupled energy calculation.
    The heat needed for reversible splitting is separate from this work bound.
    """
    _number("hydrogen_flux_mol_m2_s", hydrogen_flux_mol_m2_s)
    _number("gas_wetted_area_m2", gas_wetted_area_m2, positive=True)
    _number("projected_area_m2", projected_area_m2, positive=True)
    _number("gibbs_j_mol", gibbs_j_mol, positive=True)
    _fraction("gibbs_efficiency", gibbs_efficiency, positive=True)
    rate = hydrogen_flux_mol_m2_s * gas_wetted_area_m2
    minimum = rate * gibbs_j_mol
    input_power = minimum / gibbs_efficiency
    return {
        "hydrogen_flux_mol_m2_s": hydrogen_flux_mol_m2_s,
        "hydrogen_loss_mol_s": rate,
        "hydrogen_loss_kg_day": rate * HYDROGEN_MOLAR_MASS_KG * JULIAN_DAY,
        "water_feed_kg_day": rate * WATER_MOLAR_MASS_KG * JULIAN_DAY,
        "oxygen_coproduct_kg_day": rate * 0.5 * OXYGEN_MOLAR_MASS_KG * JULIAN_DAY,
        "minimum_replacement_work_w": minimum,
        "replacement_input_power_w": input_power,
        "replacement_input_w_m2_gas_wetted": input_power / gas_wetted_area_m2,
        "replacement_input_w_m2_projected": input_power / projected_area_m2,
    }


def hydrogen_budget(
    *, gas_volume_m3: float, gas_wetted_area_m2: float, projected_area_m2: float,
    pressure_pa: float, temperature_k: float, hydrogen_mole_fraction: float,
    permeability_mol_m_m2_s_pa: float, thickness_m: float,
    external_hydrogen_partial_pressure_pa: float = 0.0,
    gibbs_efficiency: float = 1.0,
    inflation_time_s: float = JULIAN_DAY,
    gibbs_j_mol: float = WATER_SPLITTING_GIBBS_J_MOL,
) -> dict[str, float | None]:
    """Initial inventory plus maintained-composition leakage and inflation work.

    gas_volume_m3 is actual lifting MIXTURE volume, not necessarily envelope
    volume. H₂ inventory is xPV/RT. Pass root geometry's actual barrier area;
    an internal ballonet does not magically remove its gas-contacting partition.
    initial_turnover_days is n/initial loss, NOT time until an unmaintained
    flexible gasbag empties. Production, purge and lift control are not solved.
    """
    _number("gas_volume_m3", gas_volume_m3, positive=True)
    _number("pressure_pa", pressure_pa, positive=True)
    _number("temperature_k", temperature_k, positive=True)
    _number("inflation_time_s", inflation_time_s, positive=True)
    _fraction("hydrogen_mole_fraction", hydrogen_mole_fraction)
    partial_pressure = pressure_pa * hydrogen_mole_fraction
    flux = permeation_flux(permeability_mol_m_m2_s_pa, thickness_m,
                           partial_pressure, external_hydrogen_partial_pressure_pa)
    if flux < 0:
        raise ValueError("hydrogen_budget models outward H2 loss; use permeation_flux for ingress")
    budget = measured_flux_budget(
        hydrogen_flux_mol_m2_s=flux, gas_wetted_area_m2=gas_wetted_area_m2,
        projected_area_m2=projected_area_m2, gibbs_efficiency=gibbs_efficiency,
        gibbs_j_mol=gibbs_j_mol,
    )
    inventory = partial_pressure * gas_volume_m3 / (GAS_CONSTANT * temperature_k)
    initial_work = inventory * gibbs_j_mol
    initial_input = initial_work / gibbs_efficiency
    loss = budget["hydrogen_loss_mol_s"]
    return {
        **budget,
        "hydrogen_partial_pressure_difference_pa": partial_pressure - external_hydrogen_partial_pressure_pa,
        "hydrogen_inventory_mol": inventory,
        "hydrogen_inventory_kg": inventory * HYDROGEN_MOLAR_MASS_KG,
        "initial_turnover_days": inventory / loss / JULIAN_DAY if loss else None,
        "initial_water_feed_kg": inventory * WATER_MOLAR_MASS_KG,
        "initial_oxygen_coproduct_kg": inventory * 0.5 * OXYGEN_MOLAR_MASS_KG,
        "minimum_inflation_work_j": initial_work,
        "inflation_input_energy_j": initial_input,
        "inflation_input_j_m2_projected": initial_input / projected_area_m2,
        "inflation_input_power_w": initial_input / inflation_time_s,
        "inflation_input_w_m2_projected": initial_input / inflation_time_s / projected_area_m2,
    }


def fixed_volume_exchange(
    *, initial_partial_pressures_pa: Mapping[str, float],
    external_partial_pressures_pa: Mapping[str, float],
    permeabilities_mol_m_m2_s_pa: Mapping[str, float],
    volume_m3: float, gas_wetted_area_m2: float, thickness_m: float,
    temperature_k: float, elapsed_s: float,
) -> dict[str, object]:
    """Exact isothermal, rigid-volume, no-reaction/no-production species exchange.

    p_i(t)=p_external+(p_i(0)-p_external)exp[-P_i A RT t/(t_film V)].
    This diagnostic allows O₂ ingress against positive total overpressure. It
    does NOT predict flexible-envelope buoyancy, collapse, ignition or purge.
    Never borrow oxygen permeability as a hydrogen coefficient.
    """
    species = set(initial_partial_pressures_pa)
    if not species or species != set(external_partial_pressures_pa) or species != set(permeabilities_mol_m_m2_s_pa):
        raise ValueError("all three species mappings must have the same nonempty keys")
    for name, value in (("volume_m3", volume_m3), ("gas_wetted_area_m2", gas_wetted_area_m2),
                        ("thickness_m", thickness_m), ("temperature_k", temperature_k)):
        _number(name, value, positive=True)
    _number("elapsed_s", elapsed_s)
    partial = {}
    for name in sorted(species):
        initial = initial_partial_pressures_pa[name]
        external = external_partial_pressures_pa[name]
        permeability = permeabilities_mol_m_m2_s_pa[name]
        for label, value in (("initial partial pressure", initial), ("external partial pressure", external),
                             ("permeability", permeability)):
            _number(label, value)
        exponent = permeability * gas_wetted_area_m2 * GAS_CONSTANT * temperature_k * elapsed_s / (thickness_m * volume_m3)
        partial[name] = external + (initial - external) * math.exp(-exponent)
    total = sum(partial.values())
    return {
        "partial_pressures_pa": partial,
        "total_pressure_pa": total,
        "mole_fractions": {name: p / total if total else 0.0 for name, p in partial.items()},
    }


def reference_thickness_scan(
    *, pressure_pa: float = 100_000.0, temperature_k: float = 298.15,
    radius_m: float = 100.0, hydrogen_mole_fraction: float = 0.98,
    gibbs_efficiency: float = 0.5,
) -> list[dict[str, object]]:
    """Intact full spheres; material conditions are NOT transported to lunar air.

    1500 kg/m³ is a hypothetical film density for mass sensitivity. Thickness
    changes hold measured P constant; this is a mathematical scan, not evidence
    that a one-micron defect-free wet film can be made. No gas purity control.
    """
    _number("radius_m", radius_m, positive=True)
    area = 4.0 * math.pi * radius_m**2
    cases = []
    comparators = (
        ("kurata_2022_cellophane", 5.78e-17),
        ("kurata_2022_polypropylene", 2.24e-15),
        ("wu_2014_chitin_incomplete_conditions", 0.024 * BARRER_SI),
        ("mfc_2025_incomplete_conditions", 10.7 * BARRER_SI),
        ("lmfc_2025_incomplete_conditions", 12.1 * BARRER_SI),
        ("sensitivity_100_barrer", 100.0 * BARRER_SI),
    )
    for source, permeability in comparators:
        for microns in (1, 10, 25, 100, 1000):
            thickness = microns * 1e-6
            cases.append({
                "comparator": source, "thickness_um": microns,
                "permeability_mol_m_m2_s_pa": permeability,
                "hypothetical_film_areal_mass_kg_m2": 1500.0 * thickness,
                **hydrogen_budget(
                    gas_volume_m3=4.0 * math.pi * radius_m**3 / 3.0,
                    gas_wetted_area_m2=area, projected_area_m2=area / 4.0,
                    pressure_pa=pressure_pa, temperature_k=temperature_k,
                    hydrogen_mole_fraction=hydrogen_mole_fraction,
                    permeability_mol_m_m2_s_pa=permeability, thickness_m=thickness,
                    gibbs_efficiency=gibbs_efficiency, inflation_time_s=JULIAN_YEAR_DAYS * JULIAN_DAY,
                ),
            })
    return cases
