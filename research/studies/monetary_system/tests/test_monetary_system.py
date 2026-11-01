"""Accounting controls distinguish cash, settlement, solvency and external claims."""
import copy
import json

import pytest

from research.studies.monetary_system import ledger as m
from research.studies.monetary_system import run as r


@pytest.fixture
def spec():
    return r.load_inputs()


@pytest.fixture
def book(spec):
    return r.opening(spec)


def all_snapshots(value):
    if isinstance(value, dict):
        if "balance_sheets" in value:
            yield value
        else:
            for item in value.values():
                yield from all_snapshots(item)
    elif isinstance(value, list):
        for item in value:
            yield from all_snapshots(item)


def test_opening_balance_sheets_are_explicit_and_globally_closed(book):
    snap = book.snapshot()
    assert snap["balance_sheets"]["moon_bank_a"]["financial_assets_U"] == 80
    assert snap["balance_sheets"]["moon_bank_a"]["debt_liabilities_U"] == 70
    assert snap["balance_sheets"]["moon_bank_a"]["net_worth_U"] == 10
    assert snap["aggregate_real_assets_U"] == 365
    assert snap["aggregate_net_worth_U"] == 365
    assert snap["net_foreign_assets_U"] == {"Moon": 25, "Earth": -125, "Array": 100}


def test_loan_creates_deposit_without_creating_reserves_and_principal_reverses_it(book):
    before = book.snapshot()
    m.bank_loan(book, "moon_bank_a", "moon_firms", 20)
    assert r.monetary_totals(book)["resident_L_deposits"] == 140
    assert r.monetary_totals(book)["bank_L_reserves"] == 20
    assert book.net_worth() == {a: x["net_worth_U"] for a, x in before["balance_sheets"].items()}
    m.loan_repayment(book, "moon_bank_a", "moon_firms", 20)
    assert book.snapshot()["claims"] == before["claims"]


def test_failed_payment_is_atomic_and_liquidity_is_not_capital(book):
    m.bank_loan(book, "moon_bank_a", "moon_firms", 20)
    before = copy.deepcopy(book.snapshot())
    with pytest.raises(m.LiquidityError):
        m.domestic_payment(book, "moon_firms", "moon_bank_a", "moon_households", "moon_bank_b", 20)
    assert book.snapshot() == before
    m.central_bank_loan(book, "moon_bank_a", 12, 60)
    m.domestic_payment(book, "moon_firms", "moon_bank_a", "moon_households", "moon_bank_b", 20)
    assert book.amount(m.claim("moon_cb", "moon_bank_a", "reserve")) == 0
    assert book.amount(m.claim("moon_cb", "moon_bank_b", "reserve")) == 32
    assert book.net_worth()["moon_bank_a"] == 10
    assert r.monetary_totals(book)["resident_L_deposits"] == 140


def test_same_bank_payment_uses_no_reserves(book):
    m.domestic_payment(book, "moon_firms", "moon_bank_a", "moon_households", "moon_bank_a", 20)
    assert r.monetary_totals(book)["resident_L_deposits"] == 120
    assert r.monetary_totals(book)["bank_L_reserves"] == 20


def test_transfer_between_own_accounts_preserves_net_worth(book):
    before = book.net_worth()
    m.domestic_payment(book, "moon_firms", "moon_bank_a", "moon_firms", "moon_bank_b", 5)
    assert book.net_worth() == before
    assert book.journal[-1]["flow_class"] == "financial_transaction"


def test_cancellation_preserves_deposits_and_liquidity_cannot_repair_loss(spec):
    out = r.domestic_scenario(spec)
    stages = out["stages"]
    assert out["settlement_shortfall_L"] == 12
    for name in ("after_agreed_debt_cancellation", "counterfactual_liquidity_does_not_restore_capital"):
        assert stages[name]["monetary_totals"]["resident_L_deposits"] == 140
        assert stages[name]["monetary_totals"]["bank_a_capital_U"] == -5
    assert stages["after_physical_capital_damage"]["aggregate_real_assets_U"] == 350
    assert stages["after_physical_capital_damage"]["monetary_totals"]["bank_a_capital_U"] == 10
    assert stages["after_agreed_debt_cancellation"]["aggregate_real_assets_U"] == 350
    assert out["ordinary_liquidity_request_rejected_after_loss"]
    assert stages["counterfactual_liquidity_does_not_restore_capital"]["monetary_totals"]["bank_L_reserves"] == 37


def test_collateral_limit_and_insolvency_policy_fail_before_posting(book):
    before = book.snapshot()
    with pytest.raises(ValueError, match="collateral"):
        m.central_bank_loan(book, "moon_bank_a", 11, 10)
    assert book.snapshot() == before
    m.cancel_loan(book, "moon_bank_a", "moon_firms", 15)
    before = book.snapshot()
    with pytest.raises(ValueError, match="Insolvent"):
        m.central_bank_loan(book, "moon_bank_a", 5, 60)
    assert book.snapshot() == before


def test_nonbank_and_bank_primary_bond_purchases_have_different_deposit_paths(spec):
    alternatives = r.fiscal_scenarios(spec)["alternatives"]
    nb = alternatives["nonbank_buyer"]["stages"]
    bank = alternatives["commercial_bank_buyer"]["stages"]
    assert nb["bond_issued"]["monetary_totals"]["resident_L_deposits"] == 110
    assert nb["treasury_spent"]["monetary_totals"]["resident_L_deposits"] == 120
    assert bank["bond_issued"]["monetary_totals"]["resident_L_deposits"] == 120
    assert bank["treasury_spent"]["monetary_totals"]["resident_L_deposits"] == 130
    for case in alternatives.values():
        assert case["stages"]["bond_issued"]["monetary_totals"]["bank_L_reserves"] == 10
        assert case["stages"]["treasury_spent"]["monetary_totals"]["bank_L_reserves"] == 20


def test_secondary_purchases_change_liability_form_without_new_treasury_proceeds(spec):
    for case in r.fiscal_scenarios(spec)["alternatives"].values():
        before = case["stages"]["treasury_spent"]
        after = case["stages"]["central_bank_secondary_purchase"]
        assert before["monetary_totals"]["treasury_bonds_L"] == 30
        assert after["monetary_totals"]["treasury_bonds_L"] == 30
        assert after["monetary_totals"]["resident_L_deposits"] == 130
        assert after["monetary_totals"]["bank_L_reserves"] == 30
        assert before["balance_sheets"]["moon_treasury"] == after["balance_sheets"]["moon_treasury"]
        assert before["aggregate_net_worth_U"] == after["aggregate_net_worth_U"]


def test_external_currency_gap_survives_extra_local_reserves(spec):
    out = r.foreign_scenario(spec)
    assert out["initial_foreign_funding_gap_S"] == 3
    assert out["gap_after_additional_local_credit_S"] == 3
    after = out["stages"]["additional_L_liquidity"]["monetary_totals"]
    assert after["bank_L_reserves"] == 70
    assert after["moon_total_S_settlement_assets"] == 0


def test_export_receipt_is_not_automatically_a_central_bank_reserve(spec):
    stages = r.foreign_scenario(spec)["stages"]
    receipt = stages["external_receipt_owned_by_bank"]["monetary_totals"]
    assert receipt["cb_S_settlement_assets"] == 0
    assert receipt["bank_a_S_settlement_assets"] == 6
    converted = stages["cb_purchases_external_receipt"]["monetary_totals"]
    assert converted["cb_S_settlement_assets"] == 6
    assert converted["bank_a_S_settlement_assets"] == 0
    assert converted["moon_total_S_settlement_assets"] == receipt["moon_total_S_settlement_assets"]


def test_three_world_current_accounts_explain_changes_in_external_claims(spec):
    out = r.foreign_scenario(spec)
    assert out["current_account_U"] == {"Moon": -22, "Earth": 28, "Array": -6}
    assert sum(out["current_account_U"].values()) == 0
    assert all(abs(x) < 1e-9 for x in out["NFA_transaction_reconciliation_residual_U"].values())
    assert out["import_paid_S"] == 28
    assert out["remaining_unpaid_import_S"] == 0
    assert out["stages"]["remaining_import_settled"]["monetary_totals"]["moon_total_S_settlement_assets"] == 3


def test_foreign_residence_does_not_fix_debt_denomination(spec):
    cases = r.external_debt_scenarios(spec)["alternatives"]
    lunar, foreign = cases["L"], cases["S"]
    assert lunar["debt_value_U_after"] == 7
    assert lunar["debt_value_L_after"] == 10
    assert foreign["debt_value_U_after"] == 10
    assert foreign["debt_value_L_after"] == pytest.approx(100 / 7)
    assert lunar["NFA_valuation_U"] == {"Moon": 3, "Earth": -3, "Array": 0}
    assert foreign["NFA_valuation_U"] == {"Moon": 0, "Earth": 0, "Array": 0}
    for case in cases.values():
        assert case["stages"]["opening"]["net_foreign_assets_U"] == case["stages"]["funded_and_spent"]["net_foreign_assets_U"]
        assert sum(case["current_account_U"].values()) == 0


def test_all_reference_snapshots_close_without_negative_claims(spec):
    cases = [r.domestic_scenario(spec), r.fiscal_scenarios(spec),
             r.foreign_scenario(spec), r.external_debt_scenarios(spec)]
    snapshots = list(all_snapshots(cases))
    assert len(snapshots) >= 20
    for snap in snapshots:
        assert snap["global_financial_residual_U"] == pytest.approx(0, abs=1e-9)
        assert snap["net_worth_minus_real_assets_U"] == pytest.approx(0, abs=1e-9)
        for row in snap["balance_sheets"].values():
            assert row["financial_assets_U"] + row["real_assets_U"] - row["debt_liabilities_U"] == pytest.approx(row["net_worth_U"], abs=1e-9)
        assert all(x["amount"] >= 0 for x in snap["claims"])


def test_missing_income_counterpart_fails_atomically(book):
    before = book.snapshot()
    with pytest.raises(ValueError, match="Unexplained"):
        book.post("Unfunded gift", {m.claim("moon_bank_a", "moon_firms", "deposit"): 1})
    assert book.snapshot() == before


@pytest.mark.parametrize("bad", [-1, float("inf"), float("nan"), True])
def test_invalid_quantities_are_rejected(book, bad):
    before = book.snapshot()
    with pytest.raises(ValueError):
        m.bank_loan(book, "moon_bank_a", "moon_firms", bad)
    assert book.snapshot() == before


def test_currency_prices_and_counterparties_are_checked(book):
    for value in (0, -1, float("nan")):
        with pytest.raises(ValueError):
            book.revalue("L", value)
    with pytest.raises(ValueError, match="counterparty"):
        book.post("Unknown issuer", {m.claim("missing", "moon_firms", "loan"): 1})
    with pytest.raises(ValueError, match="currency"):
        book.post("Unknown currency", {m.claim("moon_firms", "moon_bank_a", "loan", "X"): 1})


def test_result_identity_and_reproduction():
    path = r.HERE / "results/reference.json"
    product = r.results()
    assert product == json.loads(path.read_text())
    assert product["producer"]["inputs_sha256"] == r.digest(r.HERE / "inputs.json")
    assert product["producer"]["sources_sha256"] == r.digest(r.HERE / "sources.json")


def test_import_moves_produced_asset_across_territory_and_current_transfer_does_not(book):
    before = book.snapshot()
    with pytest.raises(ValueError, match="territories"):
        m.import_asset(book, "moon_firms", "moon_bank_a", "moon_households", 1)
    assert book.snapshot() == before
    m.buy_foreign_reserves(book, "moon_bank_a", 5)
    m.import_asset(book, "moon_firms", "moon_bank_a", "earth_counterparty", 5)
    assert book.real_assets_U["moon_firms"] == 125
    assert book.real_assets_U["earth_counterparty"] == 220
    assert sum(book.real_assets_U.values()) == 365
    assert book.journal[-1]["flow_class"] == "goods_trade_existing_produced_equipment"
    assert book.snapshot()["net_foreign_assets_U"]["Moon"] == 20


def test_damage_forgiveness_current_transfer_and_revaluation_keep_distinct_flow_classes(spec):
    domestic = r.domestic_scenario(spec)
    classes = {event["event"]: event["flow_class"] for event in domestic["journal"]}
    assert classes["Assumed destruction of produced capital"] == "other_volume_change"
    assert classes["Agreed cancellation of borrower debt"] == "agreed_debt_forgiveness_capital_transfer"
    assert classes["Domestic transfer and reserve settlement"] == "current_transfer"
    debt = r.external_debt_scenarios(spec)["alternatives"]["L"]
    assert debt["journal"][-1]["flow_class"] == "revaluation"


def test_context_hash_change_requires_review(tmp_path):
    path = tmp_path / "context.json"
    path.write_text("{}")
    record = {"local_context": [{"path": "context.json", "sha256": r.digest(path)}]}
    r.verify_context_identities(record, tmp_path)
    path.write_text("{\"changed\": true}")
    with pytest.raises(ValueError, match="context changed"):
        r.verify_context_identities(record, tmp_path)
