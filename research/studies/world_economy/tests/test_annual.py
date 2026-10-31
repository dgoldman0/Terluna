"""Controls for matched territorial costs, debt ownership and required trade."""

import copy
import json
from pathlib import Path

import pytest

from research.studies.world_economy import annual


def _inputs():
    return json.loads(Path(annual.__file__).with_name("annual_inputs.json").read_text())


def _construction():
    rows = [
        {"element": "ordinary_lunar_capital", "asset_territory": "moon", "beneficiary": "moon", "total_U": 50000},
        {"element": "ordinary_array_capital", "asset_territory": "array", "beneficiary": "array", "total_U": 30000},
    ]
    for name in _inputs()["incremental_asset_lives_years"]:
        world = "moon" if name in {
            "aerial_carriers_and_control", "ground_receivers", "critical_energy_reserves",
            "environment_establishment_allowance",
        } else "array"
        rows.append({"element": name, "asset_territory": world, "beneficiary": "moon", "total_U": 1000})
    return {"cost_case": "unit_test", "water_case": "unit_test", "installed_assets": rows,
            "moon_territorial_capital_U": 100000, "array_territorial_capital_U": 34000}


def _production():
    return {
        name: {"sector_gva": {"goods": gdp * 0.3, "services": gdp * 0.6},
               "net_product_taxes": gdp * 0.1, "gross_product": gdp}
        for name, gdp in {"earth": 10000, "moon": 8000, "array": 5000}.items()
    }


def _run(**kwargs):
    return annual.calculate(_construction(), _production(), _inputs(), **kwargs)


def test_tax_and_disjoint_cfc_identity_assign_assets_to_their_producer():
    result = _run()
    worlds = result["accounts"]["worlds"]
    assert worlds["moon"]["gross_world_product"] == 8000
    assert worlds["moon"]["ordinary_cfc"] == 1200
    assert worlds["moon"]["engineered_cfc"] == pytest.approx(1000 / 40 + 1000 / 40 + 1000 / 20 + 1000 / 100)
    assert worlds["array"]["engineered_cfc"] == pytest.approx(1000 / 15 + 1000 / 40 + 1000 / 100 + 1000 / 500)
    assert worlds["earth"]["engineered_cfc"] == 0
    assert worlds["moon"]["net_world_product"] == pytest.approx(8000 - 1200 - 110)
    shield = next(a for a in result["incremental_cfc_by_asset"] if a["element"] == "optical_shield_fleet")
    assert shield["asset_territory"] == "array" and shield["beneficiary"] == "moon"


def test_import_bill_is_financed_but_not_automatically_depreciated():
    baseline = _run(finance_case="patient")
    construction = _construction()
    construction["moon_territorial_capital_U"] *= 3
    changed = annual.calculate(construction, _production(), _inputs(), finance_case="patient")
    for world in annual.WORLD_IDS:
        assert changed["accounts"]["worlds"][world]["net_world_product"] == baseline["accounts"]["worlds"][world]["net_world_product"]
    assert changed["finance_history"]["moon"]["debt_at_full_population"] == pytest.approx(
        baseline["finance_history"]["moon"]["debt_at_full_population"] * 3
    )


def test_operation_allocation_changes_only_optional_residual():
    baseline, doubled = _run(), _run(operation_scale=2)
    for world in annual.WORLD_IDS:
        before, after = baseline["accounts"]["worlds"][world], doubled["accounts"]["worlds"][world]
        for field in ("gross_world_product", "net_world_product", "net_resident_income", "cash_resources_after_principal"):
            assert before[field] == after[field]
        assert before["resources_after_operations_and_principal"] - after["resources_after_operations_and_principal"] == pytest.approx(before["essential_final_operations"])


def test_debt_is_not_a_product_deduction_and_earth_receives_counterparts():
    domestic = _run(finance_case="patient", creditor_share=0)
    foreign = _run(finance_case="patient", creditor_share=1)
    for world in annual.WORLD_IDS:
        for field in ("gross_world_product", "net_world_product"):
            assert domestic["accounts"]["worlds"][world][field] == foreign["accounts"]["worlds"][world][field]
    receipts = sum(foreign["finance_history"][w]["payment_first_mature_year"] for w in ("moon", "array"))
    earth = foreign["accounts"]["worlds"]["earth"]
    assert earth["cash_resources_after_principal"] == pytest.approx(earth["net_world_product"] + receipts)
    assert foreign["accounts"]["combined"]["cash_resources_after_principal"] == pytest.approx(
        domestic["accounts"]["combined"]["cash_resources_after_principal"]
    )
    assert foreign["required_net_exports_sum_U"] == pytest.approx(0, abs=1e-10)
    for world, row in foreign["trade_requirements"].items():
        assert row["resource_pool_closure_residual_U"] == pytest.approx(0, abs=1e-10)
        if world == "earth":
            assert row["required_net_exports_U"] < 0
        else:
            assert row["required_net_exports_U"] > 0


def test_staged_loans_retire_by_year_750_without_erasing_prior_payments():
    active, retired = _run(finance_case="staged"), _run(finance_case="staged_retired")
    for world in ("moon", "array"):
        assert active["finance_history"][world]["debt_at_full_population"] > 0
        assert retired["finance_history"][world]["debt_at_full_population"] == 0
        assert retired["finance_history"][world]["payment_first_mature_year"] == 0
        assert retired["finance_history"][world]["cash_paid_before_maturity"] > 0
    assert retired["accounts"]["transfers"] == []


def test_current_interest_is_paid_during_construction_not_reborrowed_silently():
    result = _run(finance_case="current_interest")
    moon = result["finance_history"]["moon"]
    assert moon["debt_at_full_population"] == pytest.approx(50000)
    assert moon["total_capitalized_interest"] == 0
    assert moon["total_current_interest"] > 0
    assert moon["cash_paid_before_maturity"] == pytest.approx(50000 + moon["total_current_interest"])


def test_earth_legacy_debt_has_array_creditor_and_does_not_change_production():
    baseline, debt = _run(), _run(earth_legacy_debt_multiple=5)
    assert debt["finance_history"]["earth"]["debt_at_maturity"] == 50000
    payment = debt["finance_history"]["earth"]["payment_first_mature_year"]
    assert debt["accounts"]["worlds"]["earth"]["cash_resources_after_principal"] == pytest.approx(
        baseline["accounts"]["worlds"]["earth"]["cash_resources_after_principal"] - payment
    )
    assert debt["accounts"]["worlds"]["array"]["cash_resources_after_principal"] == pytest.approx(
        baseline["accounts"]["worlds"]["array"]["cash_resources_after_principal"] + payment
    )
    assert debt["required_net_exports_sum_U"] == pytest.approx(0, abs=1e-10)


def test_depletion_reserve_is_visible_and_does_not_invent_a_trade_payment():
    baseline = _run(finance_case="patient")
    changed = _run(finance_case="patient", natural_depletion_shares={"earth": 0.01, "moon": 0.02, "array": 0.005})
    for world in annual.WORLD_IDS:
        reserve = changed["accounts"]["worlds"][world]["natural_depletion"]
        assert baseline["accounts"]["worlds"][world]["net_world_product"] - changed["accounts"]["worlds"][world]["net_world_product"] == pytest.approx(reserve)
        assert baseline["trade_requirements"][world]["required_net_exports_U"] == changed["trade_requirements"][world]["required_net_exports_U"]
        assert changed["trade_requirements"][world]["resource_pool_closure_residual_U"] == pytest.approx(0, abs=1e-10)
    assert "Unpriced" in changed["natural_depletion_status"]


def test_stress_preserves_negative_resources_and_flags_export_requirements():
    construction = _construction()
    construction["moon_territorial_capital_U"] = 1e7
    result = annual.calculate(construction, _production(), _inputs(), "deferred_stress", creditor_share=1)
    moon = result["accounts"]["worlds"]["moon"]
    assert moon["net_world_product"] > 0
    assert moon["cash_resources_after_principal"] < 0
    assert all(result["trade_requirements"]["moon"]["net_exports_exceed_assumed_capacity_screen"].values())
    assert result["ratios_to_earth"]["cash_resources_after_principal"]["ratios"]["moon"] < 0


def test_longer_asset_lives_only_reduce_incremental_cfc_and_inputs_are_unchanged():
    construction, production, inputs = _construction(), _production(), _inputs()
    before = copy.deepcopy((construction, production, inputs))
    baseline = annual.calculate(construction, production, inputs)
    longer = annual.calculate(construction, production, inputs, asset_life_scale=2)
    assert (construction, production, inputs) == before
    for world in annual.WORLD_IDS:
        b, a = baseline["accounts"]["worlds"][world], longer["accounts"]["worlds"][world]
        assert a["ordinary_cfc"] == b["ordinary_cfc"]
        assert a["engineered_cfc"] == b["engineered_cfc"] / 2


@pytest.mark.parametrize("case", list(_inputs()["finance_cases"]))
def test_all_recorded_finance_cases_produce_finite_additive_ledgers(case):
    result = _run(finance_case=case)
    assert json.dumps(result, allow_nan=False)
    for residual in result["accounts"]["closure_residuals"].values():
        assert residual == pytest.approx(0, abs=1e-9)


def test_duplicate_assets_missing_lives_and_inconsistent_production_fail_explicitly():
    construction = _construction()
    construction["installed_assets"].append(copy.deepcopy(construction["installed_assets"][-1]))
    with pytest.raises(ValueError, match="Duplicate"):
        annual.calculate(construction, _production(), _inputs())
    production = _production()
    production["moon"]["gross_product"] += 1
    with pytest.raises(ValueError, match="does not equal"):
        annual.calculate(_construction(), production, _inputs())
