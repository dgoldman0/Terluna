"""Dimensional, species-gradient, energy and geometry checks for containment."""
import math

import pytest

from shared.constants import AVOGADRO, ELEMENTARY_CHARGE, GAS_CONSTANT, JULIAN_DAY
from research.studies.aerophytes.envelope import (
    BARRER_SI, WATER_SPLITTING_GIBBS_J_MOL, crossover_current_to_hydrogen_flux,
    fixed_volume_exchange, hydrogen_budget, layered_permeance, measured_flux_budget,
    permeation_flux, reference_thickness_scan,
)


def budget(**updates):
    kwargs = dict(gas_volume_m3=1000.0, gas_wetted_area_m2=100.0,
                  projected_area_m2=25.0, pressure_pa=100_100.0, temperature_k=298.15,
                  hydrogen_mole_fraction=0.98, permeability_mol_m_m2_s_pa=10 * BARRER_SI,
                  thickness_m=100e-6, gibbs_efficiency=0.5, inflation_time_s=JULIAN_DAY)
    return hydrogen_budget(**(kwargs | updates))


def test_species_gradient_is_not_mechanical_overpressure():
    # 100 Pa shell gauge pressure over 1 bar air does not give 100 Pa H2 gradient.
    result = budget()
    assert result["hydrogen_partial_pressure_difference_pa"] == pytest.approx(98_098)
    assert result["hydrogen_flux_mol_m2_s"] == pytest.approx(3.286283e-6)
    assert result["hydrogen_flux_mol_m2_s"] / permeation_flux(10 * BARRER_SI, 1e-4, 100) == pytest.approx(980.98)
    assert permeation_flux(1e-15, 1e-4, 0, 21_000) < 0  # O2 enters.


def test_analytic_fick_example_and_projection():
    flux = permeation_flux(10 * BARRER_SI, 100e-6, 100_000)
    assert flux == pytest.approx(3.35e-6)
    result = measured_flux_budget(hydrogen_flux_mol_m2_s=flux, gas_wetted_area_m2=4, projected_area_m2=1)
    assert result["replacement_input_w_m2_projected"] == pytest.approx(3.1776894)


def test_barrier_series_resistance_and_zero_permeability():
    assert layered_permeance([(1e-4, 1e-15), (1e-4, 1e-15)]) == pytest.approx(5e-12)
    assert layered_permeance([(1e-4, 0), (1e-4, 1e-15)]) == 0
    with pytest.raises(ValueError):
        layered_permeance([])


def test_thickness_and_gas_volume_are_independent_inputs():
    base = budget()
    thick = budget(thickness_m=1e-3)
    trimmed = budget(gas_volume_m3=500)
    assert thick["hydrogen_loss_kg_day"] == pytest.approx(base["hydrogen_loss_kg_day"] / 10)
    assert thick["hydrogen_inventory_kg"] == base["hydrogen_inventory_kg"]
    assert trimmed["hydrogen_inventory_kg"] == pytest.approx(base["hydrogen_inventory_kg"] / 2)
    assert trimmed["hydrogen_loss_kg_day"] == base["hydrogen_loss_kg_day"]
    assert trimmed["initial_turnover_days"] == pytest.approx(base["initial_turnover_days"] / 2)


def test_inventory_energy_efficiency_and_water_stoichiometry():
    base = budget()
    expected_moles = 0.98 * 100_100 * 1000 / (GAS_CONSTANT * 298.15)
    assert base["hydrogen_inventory_mol"] == pytest.approx(expected_moles)
    assert base["minimum_inflation_work_j"] == pytest.approx(expected_moles * WATER_SPLITTING_GIBBS_J_MOL)
    assert base["inflation_input_energy_j"] == pytest.approx(2 * base["minimum_inflation_work_j"])
    assert base["inflation_input_power_w"] * JULIAN_DAY == pytest.approx(base["inflation_input_energy_j"])
    # Rounded independent chemical molar masses conserve mass to their precision.
    assert base["water_feed_kg_day"] == pytest.approx(base["hydrogen_loss_kg_day"] + base["oxygen_coproduct_kg_day"], rel=2e-5)


def test_no_leakage_has_no_replacement_but_initial_inflation_cost():
    result = budget(permeability_mol_m_m2_s_pa=0)
    assert result["replacement_input_power_w"] == 0
    assert result["initial_turnover_days"] is None
    assert result["inflation_input_energy_j"] > 0


def test_measured_fuel_cell_crossover_uses_two_electrons():
    flux = crossover_current_to_hydrogen_flux(0.03)
    assert flux * (2 * AVOGADRO * ELEMENTARY_CHARGE) == pytest.approx(0.3)
    assert flux == pytest.approx(1.55464045e-6)
    result = measured_flux_budget(hydrogen_flux_mol_m2_s=flux, gas_wetted_area_m2=4, projected_area_m2=1)
    assert result["replacement_input_w_m2_projected"] == pytest.approx(1.474676)


def exchange(elapsed_s):
    return fixed_volume_exchange(
        initial_partial_pressures_pa={"H2": 100_100, "O2": 0, "N2": 0},
        external_partial_pressures_pa={"H2": 0, "O2": 21_000, "N2": 79_000},
        permeabilities_mol_m_m2_s_pa={"H2": 1e-12, "O2": 2e-12, "N2": 1e-12},
        volume_m3=1, gas_wetted_area_m2=1, thickness_m=1e-4,
        temperature_k=300, elapsed_s=elapsed_s,
    )


def test_species_evolution_initial_equilibrium_and_ingress():
    initial = exchange(0)
    assert initial["partial_pressures_pa"] == {"H2": 100_100, "O2": 0, "N2": 0}
    later = exchange(1)
    assert later["partial_pressures_pa"]["O2"] > 0
    assert later["total_pressure_pa"] > 100_000
    equilibrium = exchange(1e9)
    assert equilibrium["partial_pressures_pa"] == pytest.approx({"H2": 0, "O2": 21_000, "N2": 79_000})
    assert sum(equilibrium["mole_fractions"].values()) == pytest.approx(1)


def test_exponential_half_life_for_fixed_volume():
    half_life = math.log(2) * 1e-4 / (1e-12 * GAS_CONSTANT * 300)
    assert exchange(half_life)["partial_pressures_pa"]["H2"] == pytest.approx(100_100 / 2)


@pytest.mark.parametrize("updates", [
    {"thickness_m": 0}, {"gas_volume_m3": -1}, {"hydrogen_mole_fraction": 1.01},
    {"gibbs_efficiency": 0}, {"gibbs_efficiency": 1.1}, {"temperature_k": math.nan},
    {"external_hydrogen_partial_pressure_pa": 200_000},
])
def test_invalid_or_inward_budget_rejected(updates):
    with pytest.raises(ValueError):
        budget(**updates)


def test_scan_contains_wide_mass_energy_tradeoff_and_geometry_scaling():
    rows = reference_thickness_scan()
    assert len(rows) == 30
    thin, thick = rows[0], rows[4]
    assert thick["hypothetical_film_areal_mass_kg_m2"] / thin["hypothetical_film_areal_mass_kg_m2"] == pytest.approx(1000)
    assert thin["replacement_input_power_w"] / thick["replacement_input_power_w"] == pytest.approx(1000)
    larger = reference_thickness_scan(radius_m=200)[0]
    assert larger["replacement_input_w_m2_projected"] == pytest.approx(thin["replacement_input_w_m2_projected"])
    assert larger["inflation_input_j_m2_projected"] == pytest.approx(2 * thin["inflation_input_j_m2_projected"])
