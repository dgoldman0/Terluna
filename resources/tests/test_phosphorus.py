"""Conservation, feasibility and unit checks for the phosphorus return screen."""
import json

import numpy as np
import pytest

from resources import phosphorus as p
from shared.provenance import constants_changed

PRODUCT = json.loads(p.OUT.read_text())


def test_product_binds_inputs_and_constants():
    assert PRODUCT["schema"] == p.SCHEMA
    for path, sha in {**PRODUCT["producer"]["files"], **PRODUCT["producer"]["inputs"]}.items():
        assert p.digest(p.ROOT / path) == sha, path
    assert not constants_changed(PRODUCT["producer"]["constants"])
    assert PRODUCT == p.run()


def test_surface_area_and_stock_flow_units():
    result = p.stock_and_flux()
    assert np.isclose(result["occupied_surface_footprint_m2"], 3.79323281e11)
    assert np.isclose(result["carbon_kg_yr"], 3.79323281e11)
    assert np.isclose(result["dry_kg_yr"], 8.42940624e11)
    assert np.isclose(result["throughput_p_kg_yr"], 9.22830118e9)
    assert np.allclose(result["standing_p_kg"], [4.15273553e9, 4.15273553e10])
    assert np.allclose(result["stock_over_throughput_yr"], [.45, 4.5])
    # Stoichiometry changes P requirements, never area or carbon production.
    lean = p.stock_and_flux(cp_molar=500)
    assert lean["carbon_kg_yr"] == result["carbon_kg_yr"]
    assert np.isclose(lean["throughput_p_kg_yr"] / result["throughput_p_kg_yr"], 106 / 500)


def test_harvest_export_cannot_be_repaired_by_natural_recycling_alone():
    f = PRODUCT["central"]["throughput_p_kg_yr"]
    for upward in PRODUCT["upward_p_kg_yr"]:
        required = p.required_natural_recycling(f, .05, 0, upward)
        assert not required["feasible_at_or_below_unity"]
        assert p.budget(f, .05, 1, 0, upward)["sky_net_p_kg_yr"] < 0
    assert p.required_natural_recycling(f, .05, .999, PRODUCT["upward_p_kg_yr"][0])["feasible_at_or_below_unity"]


def test_algebra_conserves_recycling_and_solves_balance():
    f = PRODUCT["central"]["throughput_p_kg_yr"]
    for row in PRODUCT["recycling_requirements"]:
        if row["feasible_at_or_below_unity"]:
            b = p.budget(f, .05, row["required_fraction"], row["food_return"], row["upward_p_kg_yr"])
            assert abs(b["sky_net_p_kg_yr"]) < 1e-5
    assert p.budget(f, 0, 1, 0, 0)["sky_export_after_immediate_returns_p_kg_yr"] == 0
    assert p.budget(f, 1, 0, 1, 0)["sky_export_after_immediate_returns_p_kg_yr"] == 0


def test_finite_boxes_conserve_mass_including_unrecoverable_and_remain_nonnegative():
    for row in PRODUCT["finite_reservoir_cases"]:
        total = sum(row["initial_p_kg"].values())
        for state in row["states"]:
            assert min(state["p_kg"].values()) >= 0
            # Long-horizon matrix exponentials accumulate sub-ppm roundoff.
            assert abs(sum(state["p_kg"].values()) - total) / total < 1e-6
    initial, matrix, _ = p.transition_model(1e9, 2e9, 1e6, sediment_return_years=1e4,
                                           irreversible_fraction=.001)
    assert np.max(np.abs(matrix.sum(axis=0))) < 1e-15
    diagonal_free = matrix - np.diag(matrix.diagonal())
    assert np.min(diagonal_free) >= 0
    assert np.allclose(p.evolve(initial, matrix, 1000), p.evolve(p.evolve(initial, matrix, 500), matrix, 500), rtol=1e-10)


def test_sediment_is_a_finite_reservoir_not_an_import():
    rows = PRODUCT["finite_reservoir_cases"]
    for row in rows:
        final = row["states"][-1]
        if row["case"] == "no_sediment_return":
            assert final["sky_relative_to_initial"] < 1e-9
            assert final["p_kg"]["unrecoverable"] == 0
            assert final["p_kg"]["sediment"] > 0
        elif row["case"] == "sediment_return_no_irreversible_loss":
            assert final["sky_relative_to_initial"] > 0
            assert final["p_kg"]["unrecoverable"] == 0
        else:
            assert final["sky_relative_to_initial"] < 1e-9
            assert final["p_kg"]["unrecoverable"] > 0


def test_transport_counts_water_carrier_and_invalid_inputs():
    dense = p.lifting(1e9, 1e4, .1)
    dilute = p.lifting(1e9, 1e4, 1e-5)
    assert dilute["power_w"] == pytest.approx(dense["power_w"] * 10000)
    assert p.lifting(1e9, 0, .1)["power_w"] == 0
    for args in ((1e9, -1, .1), (1e9, 1e4, 0), (1e9, 1e4, .1, 0)):
        with pytest.raises(ValueError):
            p.lifting(*args)
    with pytest.raises(ValueError):
        p.budget(1, .05, 1.01, .99, 0)


def test_population_cases_do_not_reimpose_historical_allocation():
    cases = json.loads(p.POPULATION.read_text())["cases"]
    result = PRODUCT["population"]["cases"]
    assert set(result) == set(cases)
    for name, case in cases.items():
        assert result[name]["population_billion"] == case
        orbital = case["orbital"] * 1e9
        assert result[name]["orbital_shield_proxy_kg"] == orbital * 990e3
        assert result[name]["orbital_structure_proxy_kg"] == orbital * 15e3
    assert result["total25"]["lunar_generic_material_proxy_kg"] == [430e12, 3350e12]
