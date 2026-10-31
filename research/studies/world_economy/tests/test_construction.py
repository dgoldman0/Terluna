"""Economic/physical boundary checks for the conditional construction screen."""
import copy
import json
import math
from pathlib import Path

import pytest

from research.studies.world_economy import construction


@pytest.fixture
def inputs():
    return json.loads(Path(construction.__file__).with_name("construction_inputs.json").read_text())


def test_joint_feed_is_not_sum_of_coproduct_feeds(inputs):
    row = construction.import_account(inputs, inputs["cost_cases"]["advanced_industrial"], 0)
    assert row["source_feed_kg"] == max(row["water_feed_requirement_kg"], row["nitrogen_feed_requirement_kg"])
    assert row["source_feed_kg"] < row["water_feed_requirement_kg"] + row["nitrogen_feed_requirement_kg"]


def test_water_dominated_source_case_switches_feed_constraint(inputs):
    row = construction.import_account(inputs, inputs["cost_cases"]["patient_automation"], 5.93e19)
    assert row["source_feed_kg"] == row["water_feed_requirement_kg"]


def test_dilute_nitrogen_costs_more_without_shipping_more_mass(inputs):
    case = copy.deepcopy(inputs["cost_cases"]["advanced_industrial"])
    a = construction.import_account(inputs, case, 0)
    case["nitrogen_mass_fraction_in_source"] /= 2
    b = construction.import_account(inputs, case, 0)
    assert b["delivered_mass_kg"] == a["delivered_mass_kg"]
    assert b["source_feed_kg"] == 2 * a["source_feed_kg"]
    assert b["total_U"] > a["total_U"]


def test_zero_energy_price_keeps_capital_service_and_no_braking_purchase(inputs):
    case = copy.deepcopy(inputs["cost_cases"]["advanced_industrial"])
    case["energy_U_per_MWh"] = 0
    row = construction.import_account(inputs, case, 0)
    assert row["total_U"] > 0
    assert all(v["purchased_energy_U"] == 0 for v in row["components"].values())
    assert row["incoming_braking_energy_charged_as_purchased_MJ"] == 0
    assert row["energy_recovery_revenue_U"] == 0
    assert row["generation_capital_ceiling_U_per_continuous_W"] == 0


def test_scope_and_shielding_are_not_duplicated(inputs):
    assets = construction.installed_assets(inputs, inputs["cost_cases"]["advanced_industrial"])
    by_name = {x["element"]: x for x in assets}
    shield = by_name["optical_shield_fleet"]
    assert shield["beneficiary"] == "moon" and shield["asset_territory"] == "array"
    assert by_name["array_bulk_habitat_shield"]["quantity"] == inputs["array_population"] * 990000
    assert "array_habitat_structure" not in by_name
    assert by_name["ordinary_array_capital"]["unit_cost_U"] == by_name["ordinary_lunar_capital"]["unit_cost_U"]
    assert not any("source_plant" in a["element"] or "import_network" in a["element"] for a in assets)
    assert by_name["aerial_carriers_and_control"]["quantity"] == 1.4e14


def test_all_nine_cases_conserve_cost_and_scale_income(inputs):
    result = construction.calculate(inputs)
    assert len(result["cases"]) == 9
    for row in result["cases"]:
        imports = row["imports"]
        assert math.isclose(imports["total_U"], sum(v["total_U"] for v in imports["components"].values()) + imports["contingency_U"])
        assert math.isclose(row["combined_total_U"], imports["total_U"] + sum(a["total_U"] for a in row["installed_assets"]))
        assert math.isclose(row["combined_total_U"], row["moon_territorial_capital_U"] + row["array_territorial_capital_U"])
        assert row["included_array_structure_comparison"]["added_to_total"] is False
        proxy = row["included_ordinary_physical_proxy"]
        assert proxy["added_to_total"] is False
        assert proxy["total_allowance_U_per_resident"] >= proxy["value_U_per_resident"]
        lo, mid, hi = [row["earth_comparisons"][str(n)] for n in (50000, 100000, 250000)]
        assert math.isclose(lo["combined_in_earth_output_years"], 2 * mid["combined_in_earth_output_years"])
        assert math.isclose(hi["combined_in_earth_output_years"], .4 * mid["combined_in_earth_output_years"])


def test_shorter_schedule_changes_power_and_annual_bill_not_unit_cost(inputs):
    a = construction.calculate(inputs)["cases"][0]
    inputs["build_years"] = 250
    b = construction.calculate(inputs)["cases"][0]
    assert a["combined_total_U"] == b["combined_total_U"]
    assert b["mean_annual_build_U"] == 2 * a["mean_annual_build_U"]
    assert b["imports"]["purchased_energy_mean_PW"] == 2 * a["imports"]["purchased_energy_mean_PW"]


def test_ordinary_total_cannot_omit_the_included_physical_proxy(inputs):
    inputs["cost_cases"]["advanced_industrial"]["ordinary_installed_capital_U_per_resident"] = 100000
    with pytest.raises(ValueError, match="physical proxy"):
        construction.calculate(inputs)


@pytest.mark.parametrize("grade", [0, -1, 1.01, float("nan")])
def test_invalid_recoverable_nitrogen_grade_rejected(inputs, grade):
    inputs["cost_cases"]["advanced_industrial"]["nitrogen_mass_fraction_in_source"] = grade
    with pytest.raises(ValueError):
        construction.calculate(inputs)


def test_bound_inputs_and_committed_construction_reproduce():
    assert construction.load_inputs()["earth_population"] == 9e9
    expected = json.loads((construction.INPUT_PATH.parent / "results/construction.json").read_text())
    assert construction.results() == expected
