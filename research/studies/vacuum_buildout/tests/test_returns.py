"""Cash conservation and incidence controls for the prescribed completion path."""
from copy import deepcopy
import json
import math

import pytest

from research.studies.vacuum_buildout import returns as m


@pytest.fixture
def data():
    cfg = dict(rates=[0., .03], local_U_rate=.025,
               external_and_supplier_rate=.03, outside_supplier_U_fraction=.5)
    q = {}
    eng = dict(
        incremental_modular_cost={"test": {
            "added_steady_resource_provision_USD2041_year": 3.}},
        added_initial_cost={"inherited_installed_quote_USD2041": 100.,
                            "installed_with_preservation_and_contingency_USD2041": 137.5})
    sq = dict(strategies={"test": {"profile": "test", "premium_share": 0.}},
              premium={"staff_per_founding_scale": 0.},
              logistics={"support_staff_per_founding_scale": 0.})
    profile = {
        "domestic": {"USD2041": {"resident_consumption": 1e7,
                                   "domestic_capital_formation": 0.}},
        "commissioning": {"additional_resource_provision_USD2041": 10.},
        "external": {"USD2041": {"resource_provision": 20.}},
        "repricing": {"propellant_kg_year": 0.}}
    d = {"profiles": {"test": profile}}
    b = dict(bootstrap_years=2, commissioning_years=2, seed_scale=.05,
             bootstrap_carrying_fraction=.25,
             common_capital_fraction_of_establishment=.2,
             liquidity_buffer_years=0., seed_no_receipt_runway_years=0.,
             cultural_platform={"staff": 2}, service_kit={"staff": 1},
             service_staff_USD2041_year=10.)
    t = {"source_structure_localization": {"space_chemical": {"cost": {
        "capital_with_contingency_USD2041": 1250.,
        "capital_direct_USD2041": 1000., "preservation_USD2041": 80.,
        "initial_stock_restoration_USD2041": 40.}}}}
    cc = dict(
        underwriting=dict(min_DSCR=1.5, pledgeable_operating_cash_fraction=.7,
                          max_loan_fraction_of_budget=.75,
                          max_debt_fraction_of_seed=.25, principal_grace_years=2),
        monetary=dict(initial_reserve_fraction_of_seed=.05,
                      reserve_contribution_fraction=.2,
                      base_reserve_years_core_imports=1.,
                      U_credit_conversion_buffer_fraction=.25,
                      U_credit_loss_capital_fraction=.15,
                      domestic_holding_days=90., customer_holding_days=7.,
                      settlement_initial_share=.25, settlement_final_share=.85,
                      settlement_transition_years=10, supply_elasticity=1.,
                      incremental_distribution_spending_fraction=.01,
                      cumulative_relative_price_limit=.05,
                      annual_relative_price_limit=.005))
    old = dict(resident_work_fraction=.6, alternative_USD2041_worker_year=100.)
    return cfg, q, eng, sq, d, b, {}, t, cc, old


def annual(year, *, available=200., capital=0., appraisal=100., core_imports=0.):
    gross = max(1000., available+800.)
    return dict(
        year=year, population=100., scale=1., productive_scale=1., service_kits=0,
        starts=.1 if capital else 0., capital=capital,
        industrial_capital=capital, service_capital=0.,
        external_capital_floor=.6*capital, local_capital=.4*capital,
        capital_resource_provision=3. if capital else 0., appraisal=appraisal,
        operation=dict(available_USD2041=available, gross_USD2041=gross,
                       total_cost_USD2041=gross-available,
                       core_imports_USD2041=core_imports,
                       resource_provision_USD2041=10., added_technical_payroll=5.,
                       revenue_by_sector_USD2041={"test": gross}))


def path(rows, *, working=0.):
    return dict(strategy="test", market="test", seed_hardware=10000.,
                seed_carry=20., initial_working=working,
                seed=10020.+working, annual=rows)


def pair(p, data):
    control = m.finance(p, data, "equity_completion")
    credit = m.finance(p, data, "latest_U_completion", control)
    return control, credit


def test_negative_operations_and_construction_receive_dated_equity(data):
    p = path([annual(2, available=-100., capital=200.),
              annual(3, available=-80.)])
    control, credit = pair(p, data)
    assert credit["annual"] == control["annual"]
    first, second = credit["annual"]
    assert first["pre_topup_cash"] == -300.
    assert first["contribution"] == 300.
    assert first["topup_by_purpose"] == dict(operating=100., debt=0.,
                                             capital=200., liquidity=0.)
    assert second["contribution"] == 80.
    assert all(r["draw"] == r["debt"] == r["dividend"] == 0 for r in credit["annual"])
    result = m.metrics(p, credit, data, 3)
    assert result["first_negative_operating_year"] == 2
    assert result["first_external_topup_year"] == 2
    assert result["peak_annual_external_equity"] == 300.
    assert result["additional_external_equity"] == 380.


def test_cash_reserves_principal_and_grace_conserve_funding(data):
    p = path([annual(y, capital=1000. if y == 2 else 0.) for y in range(2, 26)])
    _, sim = pair(p, data)
    assert sim["annual"][0]["draw"] == pytest.approx(750.)
    opening = p["initial_working"]
    debt = 0.
    reserve = sim["initial_reserve"]
    for row in sim["annual"]:
        assert row["cash"] == pytest.approx(
            opening+row["operating_available"]+row["draw"]+row["contribution"]
            -row["interest"]-row["principal"]-row["reserve_contribution"]
            -row["capital"]-row["dividend"])
        debt += row["draw"]-row["principal"]
        reserve += row["reserve_contribution"]
        assert row["debt"] == pytest.approx(debt, abs=1e-10)
        assert row["reserves"] == pytest.approx(reserve)
        assert sum(row["topup_by_purpose"].values()) == pytest.approx(row["contribution"])
        assert row["draw"] == pytest.approx(row["external_draw"]+row["supplier_U_draw"]+row["local_U_draw"])
        assert row["debt"] <= .25*p["seed"]+1e-9
        assert .4*row["U_debt"] <= row["reserves"]+1e-9
        assert abs(row["cash_residual"]) < 1e-9
        opening = row["cash"]
    by_year = {r["year"]: r for r in sim["annual"]}
    assert by_year[2]["interest"] == 0
    assert by_year[3]["principal"] == by_year[4]["principal"] == 0
    assert by_year[5]["principal"] == pytest.approx(750./20)
    assert by_year[24]["debt"] == pytest.approx(0., abs=1e-9)
    assert by_year[25]["interest"] == pytest.approx(0., abs=1e-9)
    assert sum(r["principal"] for r in sim["annual"]) == pytest.approx(750.)


def test_local_relief_preserves_supplier_and_outside_contracts(data):
    p = path([annual(y, capital=1000. if y == 2 else 0.) for y in range(2, 6)])
    _, relieved = pair(p, data)
    neutral_data = deepcopy(data)
    neutral_data[0]["local_U_rate"] = .03
    _, neutral = pair(p, neutral_data)
    issued = relieved["annual"][0]
    assert issued["local_U_draw"] == pytest.approx(300.)
    assert issued["supplier_U_draw"] == pytest.approx(225.)
    assert issued["external_draw"] == pytest.approx(225.)
    for reduced, unchanged in zip(relieved["annual"][1:], neutral["annual"][1:]):
        assert reduced["supplier_interest"] == unchanged["supplier_interest"]
        assert reduced["external_interest"] == unchanged["external_interest"]
        assert reduced["local_U_interest"]/.025 == pytest.approx(unchanged["local_U_interest"]/.03)
        assert unchanged["interest"]-reduced["interest"] == pytest.approx(
            unchanged["local_U_interest"]-reduced["local_U_interest"])


def test_coverage_lock_retains_cash_and_blocks_new_credit(data):
    p = path([annual(2, capital=1000.), annual(3, available=40.), annual(4)])
    _, sim = pair(p, data)
    row = sim["annual"][1]
    assert row["gross_DSCR"] > 1.5
    assert row["DSCR"] == pytest.approx(.7*row["gross_DSCR"])
    assert row["DSCR"] < 1.5
    assert row["covenant_lock"]
    assert row["draw"] == row["dividend"] == 0.
    assert row["cash"] > 0.


@pytest.mark.parametrize("appraisal", [0., -1.])
def test_nonpositive_marginal_value_denies_debt_without_cancelling_build(data, appraisal):
    p = path([annual(2, capital=1000., appraisal=appraisal)])
    control, sim = pair(p, data)
    assert sim["annual"] == control["annual"]
    assert sim["annual"][0]["capital"] == 1000.
    assert sim["annual"][0]["draw"] == 0.
    assert sim["annual"][0]["contribution"] == 800.


def test_price_rejected_local_credit_becomes_equity_not_another_loan(data):
    constrained = deepcopy(data)
    constrained[4]["profiles"]["test"]["domestic"]["USD2041"]["resident_consumption"] = 100.
    p = path([annual(2, capital=1000.)])
    _, easy = pair(p, data)
    _, tight = pair(p, constrained)
    wide, narrow = easy["annual"][0], tight["annual"][0]
    assert narrow["local_U_draw"] < wide["local_U_draw"]
    assert narrow["supplier_U_draw"] == wide["supplier_U_draw"]
    assert narrow["external_draw"] == wide["external_draw"]
    assert narrow["contribution"]-wide["contribution"] == pytest.approx(
        wide["local_U_draw"]-narrow["local_U_draw"])
    assert narrow["relative_price"] <= 1.005+1e-10
    assert not narrow["price_breach"]


def test_zero_reserve_headroom_denies_U_creation(data):
    data[8]["monetary"]["initial_reserve_fraction_of_seed"] = 0.
    p = path([annual(2, capital=1000., core_imports=10000.)])
    _, sim = pair(p, data)
    row = sim["annual"][0]
    assert row["local_U_draw"] == row["supplier_U_draw"] == row["U_debt"] == 0.
    assert row["external_draw"] > 0.
    assert row["relative_price"] == 1.


def test_rescue_equity_does_not_enlarge_the_original_seed_debt_ceiling(data):
    p = path([annual(2, available=-1e6),
              annual(3, available=20000., capital=50000.)])
    _, sim = pair(p, data)
    assert sim["annual"][0]["contribution"] == 1e6
    assert sim["annual"][0]["draw"] == 0.
    assert sim["annual"][1]["debt"] == pytest.approx(.25*p["seed"])
    assert sim["annual"][1]["contribution"] > 0.


def test_custody_deferral_and_later_release_conserve_owner_claims(data):
    data[4]["profiles"]["test"]["domestic"]["USD2041"]["resident_consumption"] = 100.
    p = path([annual(2), annual(3)])
    # Isolate the payout controller with fixed reference prices; no credit is drawn.
    reference = {"annual": [{"year": 2, "price": 1.}, {"year": 3, "price": 2.}]}
    sim = m.finance(p, data, "latest_U_completion", reference)
    first, second = sim["annual"]
    assert first["draw"] == second["draw"] == 0.
    assert first["owner_custody"] > 0.
    assert first["owner_payout"] < first["dividend"]
    assert first["relative_price"] <= 1.005+1e-10
    assert second["owner_custody"] == pytest.approx(0.)
    assert second["owner_payout"] > second["dividend"]
    custody = 0.
    for row in sim["annual"]:
        assert row["owner_payout"]+row["owner_custody"]-custody == pytest.approx(row["dividend"])
        custody = row["owner_custody"]
    assert sum(r["owner_payout"] for r in sim["annual"]) == pytest.approx(sum(r["dividend"] for r in sim["annual"]))
    result = m.metrics(p, sim, data, 2)["discounted"]["0.0"]
    assert result["terminal_net_book_claim_PV"] == pytest.approx(
        first["cash"]+first["reserves"]+first["owner_custody"]-first["debt"])


def test_dividend_retains_next_capital_buffer_increase_and_reserve_contribution(data):
    data[5]["liquidity_buffer_years"] = 1.
    p = path([annual(2, available=2000.),
              annual(3, available=100., capital=100., core_imports=10000.)])
    sim = m.finance(p, data, "equity_completion")
    first, second = sim["annual"]
    # Next year needs a 900 buffer, 100 capital and 20 reserve contribution,
    # offset by 100 of operating cash: retain 920 before distributing 1080.
    assert first["dividend"] == pytest.approx(1080.)
    assert first["cash"] == pytest.approx(920.)
    assert second["reserve_contribution"] == pytest.approx(20.)
    assert second["contribution"] == 0.
    assert second["cash"] >= 900.


def test_initial_allowances_are_separate_from_real_procurement(data):
    _, _, eng, _, _, b, _, t, _, _ = data
    assert m.initial_provisions(b, t, eng) == pytest.approx((100., 50., 12.5))


def test_metrics_include_exact_seed_dates_topups_and_terminal_claims(data):
    p = path([annual(2, available=-100., capital=200.), annual(3, available=-80.)], working=40.)
    _, sim = pair(p, data)
    result = m.metrics(p, sim, data, 3)
    seed_annual = (p["seed_hardware"]+p["seed_carry"])/2
    for rate in data[0]["rates"]:
        expected = -seed_annual-sim["initial_reserve"]-(seed_annual+40.)/(1+rate)
        expected += sum((r["owner_payout"]-r["contribution"])/(1+rate)**r["year"] for r in sim["annual"])
        values = result["discounted"][str(rate)]
        assert values["sponsor_cash_NPV"] == pytest.approx(expected)
        terminal = sim["annual"][-1]
        claim = (terminal["cash"]+terminal["reserves"]+terminal["owner_custody"]-terminal["debt"])/(1+rate)**3
        assert values["terminal_net_book_claim_PV"] == pytest.approx(claim)
        assert values["sponsor_with_terminal_book_claim_NPV"] == pytest.approx(expected+claim)
    assert result["nonnegative_single_root_IRR"] is None


def test_real_accounts_ignore_financing_and_remove_known_transfers_once(data):
    p = path([annual(2, capital=1000.), annual(3)])
    control, sim = pair(p, data)
    base = m.metrics(p, control, data, 3)
    financed = m.metrics(p, sim, data, 3)
    for rate in data[0]["rates"]:
        one, two = base["discounted"][str(rate)], financed["discounted"][str(rate)]
        for field in ("unlevered_operating_less_capital_NPV", "full_price_socioeconomic_NPV",
                      "half_price_socioeconomic_NPV", "discounted_economic_cost",
                      "discounted_exports", "discounted_resident_value"):
            assert one[field] == two[field]
    # Seed removes commissioning, original/stock/new machinery provisions and carrying provision.
    seed_provision = .05*10.+(.2+.05*.8)*100.+.05*(50.+12.5)
    carry_provision = .05*(20.+3.)*.25*2
    initial_cost = 10020.-seed_provision-carry_provision
    # Each occupied year replaces 20 payroll + 5 extra payroll with 60 workers' opportunity cost.
    cost = initial_cost + (800.+1000.-10.-3.-25.+6000.) + (800.-10.-25.+6000.)
    result = financed["discounted"]["0.0"]
    assert result["discounted_economic_cost"] == pytest.approx(cost)
    assert result["unlevered_operating_less_capital_NPV"] == -10020.+200.-1000.+200.
    assert result["full_price_socioeconomic_NPV"] == pytest.approx(2000.+2e7-cost)
    assert result["half_price_socioeconomic_NPV"] == pytest.approx(.5*(2000.+2e7)-cost)


def test_saved_product_reproduction_and_matched_real_accounts():
    saved = json.loads((m.HERE/"results/returns.json").read_text())
    result = m.product()
    assert result == saved
    assert result["producer_sha256"] == m.m.digest(m.__file__)
    for path_name, digest in result["dependencies"].items():
        assert m.m.digest(m.HERE/path_name) == digest
    grouped = {}
    for row in result["horizons"]:
        grouped.setdefault((row["strategy"], row["market"], row["horizon"]), {})[row["mode"]] = row
        assert sum(row["external_equity_by_purpose"].values()) == pytest.approx(row["additional_external_equity"])
        assert not row["price_breach_years"]
        assert row["population"] > 0
    assert len(grouped) == 9*3*2
    for cases in grouped.values():
        control, credit = cases["equity_completion"], cases["latest_U_completion"]
        for rate, baseline in control["discounted"].items():
            other = credit["discounted"][rate]
            for field in ("unlevered_operating_less_capital_NPV", "full_price_socioeconomic_NPV",
                          "half_price_socioeconomic_NPV", "discounted_economic_cost"):
                assert other[field] == baseline[field]
    for path_rows in result["financing_paths"]:
        for previous, following in zip(path_rows["annual"], path_rows["annual"][1:]):
            assert not (previous["dividend"] > 1. and following["contribution"] > 1.), (
                path_rows["strategy"], path_rows["market"], path_rows["mode"],
                previous["year"], previous["dividend"], following["contribution"])
        for row in path_rows["annual"]:
            assert abs(row["cash_residual"]) <= 1.
            assert math.isclose(row["draw"], row["local_U_draw"]+row["supplier_U_draw"]+row["external_draw"], rel_tol=1e-12, abs_tol=1e-5)
            if row["covenant_lock"]:
                assert row["draw"] == row["dividend"] == 0.
    for physical in result["physical_and_operating_paths"]:
        assert all(row["resource"]["feasible"] for row in physical["annual"])
        at_cap = next(row for row in physical["annual"] if row["year"] == 200)
        assert at_cap["population"] == pytest.approx(6468289.238000949)
