"""Accounting identities and scenario boundaries for the productivity screen."""

import copy
import json
from pathlib import Path

import pytest

from research.studies.world_economy import productivity
from research.studies.world_economy.productivity import evaluate_productivity, load_inputs, validate_inputs


@pytest.fixture
def inputs():
    return load_inputs(Path(productivity.__file__).with_name("productivity_inputs.json"))


def run(inputs, **kwargs):
    arguments = dict(earth_population=9e9, moon_population=8e9, array_population=5e9,
                     earth_gdp_per_resident=100_000, scenario_id="close_parity")
    arguments.update(kwargs)
    return evaluate_productivity(inputs, **arguments)


def test_equal_realised_per_resident_gva_returns_population_ratios(inputs):
    scenario = inputs["scenarios"]["close_parity"]
    scenario["moon_gva_per_resident_ratios"] = {sector["id"]: 1.0 for sector in inputs["sectors"]}
    scenario["array_gva_per_resident_ratio"] = 1.0
    result = run(inputs)
    assert result["gross_product_ratios_to_earth"] == pytest.approx({
        "moon": 8 / 9, "array": 5 / 9, "moon_plus_array": 13 / 9,
    })
    other_populations = run(inputs, earth_population=30, moon_population=12, array_population=6)
    assert other_populations["gross_product_ratios_to_earth"] == pytest.approx({
        "moon": 0.4, "array": 0.2, "moon_plus_array": 0.6,
    })


@pytest.mark.parametrize("scenario,expected_resident,expected_moon,expected_combined", [
    ("constrained", 0.8275, 0.7355555555555555, 1.1522222222222223),
    ("close_parity", 1.0490, 0.9324444444444444, 1.488),
    ("realised_specialisation", 1.2570, 1.1173333333333333, 1.9506666666666668),
])
def test_declared_scenarios_have_reviewed_totals(inputs, scenario, expected_resident, expected_moon, expected_combined):
    result = run(inputs, scenario_id=scenario)
    assert result["assumed_gva_per_resident_ratios_to_earth"]["moon"] == pytest.approx(expected_resident)
    assert result["gross_product_ratios_to_earth"]["moon"] == pytest.approx(expected_moon)
    assert result["gross_product_ratios_to_earth"]["moon_plus_array"] == pytest.approx(expected_combined)


@pytest.mark.parametrize("weights,scenario,expected", [
    ("material_weighted", "constrained", 0.7155555555555555),
    ("material_weighted", "close_parity", 0.9204444444444445),
    ("material_weighted", "realised_specialisation", 1.1106666666666665),
    ("service_weighted", "constrained", 0.7493333333333332),
    ("service_weighted", "close_parity", 0.9377777777777779),
    ("service_weighted", "realised_specialisation", 1.1173333333333333),
])
def test_alternative_weights_are_genuine_sensitivity(inputs, weights, scenario, expected):
    assert run(inputs, weight_case=weights, scenario_id=scenario)["gross_product_ratios_to_earth"]["moon"] == pytest.approx(expected)


def test_price_scale_equivariance_and_complete_sector_account(inputs):
    low = run(inputs, earth_gdp_per_resident=50_000)
    high = run(inputs, earth_gdp_per_resident=250_000)
    assert low["gross_product_ratios_to_earth"] == pytest.approx(high["gross_product_ratios_to_earth"])
    for world in ("earth", "moon", "array"):
        a, b = low["worlds"][world], high["worlds"][world]
        assert b["gross_product"] == pytest.approx(5 * a["gross_product"])
        assert sum(a["sector_gva"].values()) == pytest.approx(a["gva"])
        assert sum(a["sector_gross_product"].values()) == pytest.approx(a["gross_product"])
    assert low["worlds"]["earth"]["gross_product"] == pytest.approx(450e12)


def test_tax_bridge_uses_earth_market_gdp_to_define_basic_price_gva(inputs):
    result = run(inputs, earth_population=10, moon_population=10, array_population=10,
                 earth_gdp_per_resident=100, net_product_tax_shares={"earth": .2, "moon": .1, "array": 0})
    earth, moon, array = (result["worlds"][key] for key in ("earth", "moon", "array"))
    assert earth["gross_product"] == pytest.approx(1000)
    assert earth["gva"] == pytest.approx(800)
    assert moon["gva"] == pytest.approx(839.2)
    assert moon["gross_product"] == pytest.approx(839.2 / .9)
    assert array["gross_product"] == pytest.approx(800)
    assert result["gross_product_ratios_to_earth"]["moon"] == pytest.approx(1.049 * .8 / .9)
    for world in (earth, moon, array):
        assert world["gva"] + world["net_product_taxes"] == pytest.approx(world["gross_product"])
        assert sum(world["sector_gross_product"].values()) == pytest.approx(world["gross_product"])


def test_common_tax_share_preserves_relative_gross_product(inputs):
    baseline = run(inputs)
    taxed = run(inputs, net_product_tax_shares={"earth": .12, "moon": .12, "array": .12})
    assert baseline["gross_product_ratios_to_earth"] == pytest.approx(taxed["gross_product_ratios_to_earth"])
    assert taxed["worlds"]["earth"]["gross_product"] == pytest.approx(baseline["worlds"]["earth"]["gross_product"])


def test_array_is_realised_aggregate_without_export_bonus(inputs):
    result = run(inputs, scenario_id="realised_specialisation")
    assert set(result["worlds"]) == {"earth", "moon", "array"}
    assert result["comparison_scope"] == "earth_moon_array_are_three_peer_territorial_economies"
    assert result["moon_plus_array"]["scope"] == "derived_aggregate_not_an_additional_world"
    array = result["worlds"]["array"]
    assert array["gross_product"] == pytest.approx(5e9 * 100_000 * 1.5)
    assert array["sector_allocation"] == "reference_weights_accounting_only_not_array_industry_estimate"
    assert result["evidence"]["probability_interpretation"] == "none"
    assert result["moon_plus_array"]["gross_product"] == pytest.approx(
        result["worlds"]["moon"]["gross_product"] + array["gross_product"])
    inputs["scenarios"]["realised_specialisation"]["export_bonus"] = 1e12
    with pytest.raises(ValueError, match="exactly"):
        run(inputs, scenario_id="realised_specialisation")


def test_inputs_and_returned_metadata_do_not_share_mutable_state(inputs):
    before = copy.deepcopy(inputs)
    result = run(inputs)
    assert inputs == before
    result["evidence"]["statement"] = "changed"
    assert inputs == before


def test_zero_population_has_zero_production_and_no_per_resident_statistic(inputs):
    result = run(inputs, moon_population=0, array_population=0)
    for key in ("moon", "array"):
        assert result["worlds"][key]["gross_product"] == 0
        assert result["worlds"][key]["gva_per_resident"] is None


@pytest.mark.parametrize("invalid", [float("nan"), float("inf"), -1, True, "100"])
def test_invalid_scale_and_population_are_rejected(inputs, invalid):
    with pytest.raises(ValueError):
        run(inputs, earth_gdp_per_resident=invalid)
    with pytest.raises(ValueError):
        run(inputs, moon_population=invalid)


@pytest.mark.parametrize("field", ["weight", "ratio", "array_ratio"])
def test_nonfinite_judgment_inputs_are_rejected(inputs, field):
    if field == "weight":
        inputs["weight_cases"]["central"]["food_biomass"] = float("nan")
    elif field == "ratio":
        inputs["scenarios"]["close_parity"]["moon_gva_per_resident_ratios"]["food_biomass"] = float("inf")
    else:
        inputs["scenarios"]["close_parity"]["array_gva_per_resident_ratio"] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        validate_inputs(inputs)


def test_weight_sum_and_unique_sector_checks(inputs):
    inputs["weight_cases"]["central"]["food_biomass"] += .01
    with pytest.raises(ValueError, match="sum to one"):
        validate_inputs(inputs)
    inputs["weight_cases"]["central"]["food_biomass"] -= .01
    inputs["sectors"].append(copy.deepcopy(inputs["sectors"][0]))
    with pytest.raises(ValueError, match="unique"):
        validate_inputs(inputs)


@pytest.mark.parametrize("taxes", [{"earth": 1}, {"moon": float("nan")}, {"unknown_world": 0}])
def test_invalid_tax_configuration_is_rejected(inputs, taxes):
    with pytest.raises(ValueError):
        run(inputs, net_product_tax_shares=taxes)


def test_reader_checks_duplicate_json_keys(tmp_path, inputs):
    path = tmp_path / "duplicate.json"
    text = json.dumps(inputs)
    path.write_text(text[:-1] + ', "scenarios": {}}', encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON key"):
        load_inputs(path)
