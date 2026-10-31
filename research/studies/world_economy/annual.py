"""Three-peer-world annual accounts from supplied cost and production scenarios.

This module performs no file I/O. Pass the separately recorded annual inputs,
one construction case, and productivity's ``worlds`` mapping. Prices, financing
and lifetime assumptions stay in that supplied input record. Operations remain
a separate allocation; neither operations nor debt reduce standard NWP.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping
from numbers import Integral, Real

from . import accounts, finance


WORLD_IDS = ("earth", "moon", "array")
INPUT_SCHEMA = "terluna.world-economy.annual-inputs/1"
OUTPUT_SCHEMA = "terluna.world-economy.annual/1"
ORDINARY_ASSETS = {"ordinary_lunar_capital", "ordinary_array_capital"}
CAPITAL_KEYS = {"moon": "moon_territorial_capital_U", "array": "array_territorial_capital_U"}


def _number(value, name, *, positive=False, fraction=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a number")
    result = float(value)
    if not math.isfinite(result) or result < 0 or (positive and result == 0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}")
    if fraction and result > 1:
        raise ValueError(f"{name} must not exceed one")
    return result


def _year(value, name, *, zero=False):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < (0 if zero else 1):
        raise ValueError(f"{name} must be an integer in the allowed range")
    return int(value)


def _world_shares(values, name):
    if not isinstance(values, Mapping) or set(values) != set(WORLD_IDS):
        raise ValueError(f"{name} must name Earth, Moon and array exactly")
    return {world: _number(value, name + "." + world, fraction=True) for world, value in values.items()}


def _finite_output(value):
    if isinstance(value, Mapping):
        for child in value.values():
            _finite_output(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            _finite_output(child)
    elif isinstance(value, Real) and not isinstance(value, bool) and not math.isfinite(value):
        raise ValueError("Annual calculation exceeded finite numerical range")


def _finance_history(cost, case, build_years):
    mode = case["mode"]
    full_year = _year(case["full_population_year"], "full_population_year")
    if full_year < build_years:
        raise ValueError("full_population_year precedes construction completion")
    gap = full_year - build_years
    f = _number(case["financed_fraction"], "financed_fraction", fraction=True)
    r = _number(case["annual_rate"], "annual_rate")
    cap = _number(case["capitalized_interest_fraction"], "capitalized_interest_fraction", fraction=True)
    term = _year(case["repayment_years"], "repayment_years")
    if mode == "retired":
        return {"mode": mode, "debt_at_full_population": 0.0,
                "interest_first_mature_year": 0.0, "principal_repaid_first_mature_year": 0.0,
                "payment_first_mature_year": 0.0,
                "history_status": "Legacy liabilities assumed retired; past payments are not reconstructed."}
    if mode == "staged":
        history = finance.staged_construction_finance(
            cost, build_years, _year(case["tranche_years"], "tranche_years"), r, term, f, cap, gap
        )
        history.update({"mode": mode, "debt_at_full_population": history["debt_at_maturity"],
                        "gap_treatment": "Tranches continue amortizing during the occupancy gap."})
        return history
    if mode != "completion_loan":
        raise ValueError(f"Unknown finance mode: {mode}")
    history = finance.construction_finance(cost, build_years, r, f, cap, gap)
    service = finance.debt_service_first_year(history["debt_at_full_occupancy"], r, term)
    history.update({"mode": mode, "debt_at_full_population": history["debt_at_full_occupancy"],
                    "interest_first_mature_year": service["interest"],
                    "principal_repaid_first_mature_year": service["principal_repaid"],
                    "payment_first_mature_year": service["payment"],
                    "cash_paid_before_maturity": history["unfinanced_construction"] + history["total_current_interest"],
                    "gap_treatment": "Amortization deferred; interest continues under the construction capitalization policy."})
    return history


def calculate(
    construction_case, productivity_worlds, inputs, finance_case="legacy_retired", *,
    creditor_share=None, earth_legacy_debt_multiple=None, ordinary_cfc_share=None,
    natural_depletion_shares=None, operation_scale=1.0, asset_life_scale=1.0,
    additional_transfers=(),
):
    """Calculate one matched annual scenario, with explicit optional sensitivities.

    Each productivity world supplies ``sector_gva``, ``net_product_taxes`` and
    ``gross_product``. Extra CFC follows installed-asset territory, not its
    beneficiary. The initial atmosphere/water import bill creates no automatic
    depreciation. Future Earth uses the same ordinary CFC share convention.

    Required trade is a financing implication, not an export forecast. The
    resource pool after principal retains the NWP depletion reserve. The
    operations residual is also reported, but the operating allocation changes
    no GDP, NWP, income or pre-allocation cash result.
    """
    inputs = copy.deepcopy(inputs)
    if inputs.get("schema") != INPUT_SCHEMA or set(inputs.get("worlds", [])) != set(WORLD_IDS):
        raise ValueError("Unsupported annual input schema or world boundary")
    if set(productivity_worlds) != set(WORLD_IDS):
        raise ValueError("Production must specify exactly three peer worlds")
    build = _year(inputs["build_years"], "build_years")
    case = inputs["finance_cases"][finance_case]
    shares = _world_shares(inputs["ordinary_cfc_shares"], "ordinary_cfc_shares")
    if ordinary_cfc_share is not None:
        value = _number(ordinary_cfc_share, "ordinary_cfc_share", fraction=True)
        shares = dict.fromkeys(WORLD_IDS, value)
    depletion = _world_shares(
        inputs["natural_depletion_shares"] if natural_depletion_shares is None else natural_depletion_shares,
        "natural_depletion_shares",
    )
    ownership = _number(inputs["earth_creditor_share"] if creditor_share is None else creditor_share,
                        "creditor_share", fraction=True)
    operation_scale = _number(operation_scale, "operation_scale")
    asset_life_scale = _number(asset_life_scale, "asset_life_scale", positive=True)
    life_inputs = inputs["incremental_asset_lives_years"]
    lives = {name: _number(life, name + " life", positive=True) * asset_life_scale
             for name, life in life_inputs.items()}
    operations = inputs["essential_final_operations"]
    base_operation = _number(operations["baseline_share_of_gwp"], "operations baseline", fraction=True)
    extra_operation = _world_shares(operations["extra_shares_of_gwp"], "operations extras")

    asset_rows = []
    seen = set()
    for asset in construction_case["installed_assets"]:
        name, world = asset["element"], asset["asset_territory"]
        if name in seen:
            raise ValueError(f"Duplicate installed asset: {name}")
        seen.add(name)
        if world not in WORLD_IDS:
            raise ValueError(f"Unknown asset territory: {world}")
        if name in ORDINARY_ASSETS:
            continue
        if name not in lives:
            raise ValueError(f"No explicit incremental life for asset: {name}")
        cost = _number(asset["total_U"], name + " installed cost")
        asset_rows.append({"element": name, "asset_territory": world,
                           "beneficiary": asset.get("beneficiary"), "replacement_cost_U": cost,
                           "life_years": lives[name], "annual_cfc_U": cost / lives[name]})
    if set(lives) - seen:
        raise ValueError("Construction case is missing an expected incremental asset")

    world_inputs = {}
    for world, production in productivity_worlds.items():
        gdp = _number(production["gross_product"], world + " gross product", positive=True)
        sectors = dict(production["sector_gva"])
        if "net_product_taxes_adjustment" in sectors:
            raise ValueError("Sector names collide with the separate product-tax adjustment")
        taxes = float(production["net_product_taxes"])
        if not math.isfinite(taxes):
            raise ValueError("Net product taxes must be finite")
        sectors["net_product_taxes_adjustment"] = taxes
        if not math.isclose(math.fsum(sectors.values()), gdp, rel_tol=1e-10, abs_tol=1e-8):
            raise ValueError("Sector GVA plus net product taxes does not equal gross product")
        world_inputs[world] = {
            "sector_gva": sectors, "ordinary_cfc": shares[world] * gdp,
            "engineered_cfc": math.fsum(row["annual_cfc_U"] for row in asset_rows if row["asset_territory"] == world),
            "natural_depletion": depletion[world] * gdp,
            "essential_final_operations": operation_scale * (base_operation + extra_operation[world]) * gdp,
        }

    histories = {}
    flows = list(additional_transfers)
    for world, cost_key in CAPITAL_KEYS.items():
        cost = _number(construction_case[cost_key], cost_key)
        history = _finance_history(cost, case, build)
        history["initial_construction_cost_U"] = cost
        history["earth_creditor_share"] = ownership
        history["limitations"] = inputs["construction_history_limits"]
        histories[world] = history
        flows.extend(accounts.debt_transfers(
            world, "earth", history["interest_first_mature_year"],
            history["principal_repaid_first_mature_year"], ownership,
        ))
    earth_debt = inputs["earth_legacy_debt"]
    multiple = _number(earth_debt["multiple_of_earth_gwp"] if earth_legacy_debt_multiple is None
                       else earth_legacy_debt_multiple, "earth_legacy_debt_multiple")
    if earth_debt["creditor"] != "array":
        raise ValueError("Earth legacy sensitivity names the array as its explicit creditor")
    earth_service = accounts.mature_service_transfers(
        "earth", "array", multiple * productivity_worlds["earth"]["gross_product"],
        _number(earth_debt["annual_rate"], "Earth debt rate"),
        _year(earth_debt["repayment_years"], "Earth debt term"),
    )
    flows.extend(earth_service.pop("transfers"))
    histories["earth"] = {**earth_service, "legacy_multiple_of_gwp": multiple,
                           "status": earth_debt["statement"]}
    ledger = accounts.consolidated_accounts(world_inputs, flows)

    capacity_screens = [_number(x, "net-export screen", fraction=True)
                        for x in inputs["net_export_capacity_share_screens"]]
    trade = {}
    for world, row in ledger["worlds"].items():
        nx = -row["net_primary_income_from_other_worlds"] - row["net_principal_receipts"]
        absorption = row["gross_world_product"] - nx
        consumption_before_depletion = absorption - row["total_cfc"]
        consumption_after_reserve = consumption_before_depletion - row["natural_depletion"]
        trade[world] = {
            "required_net_exports_U": nx, "required_net_exports_share_of_gwp": nx / row["gross_world_product"],
            "gross_domestic_absorption_U": absorption,
            "consumption_capacity_after_cfc_before_depletion_reserve_U": consumption_before_depletion,
            "depletion_reserve_U": row["natural_depletion"],
            "resource_pool_after_capital_maintenance_reserves_U": consumption_after_reserve,
            "resource_pool_closure_residual_U": consumption_after_reserve - row["cash_resources_after_principal"],
            "net_exports_exceed_assumed_capacity_screen": {
                str(screen): nx > screen * row["gross_world_product"] for screen in capacity_screens
            },
            "screen_status": "Assumed sensitivity thresholds only; no trade matrix, export capacity, demand or price equilibrium has been verified.",
        }
    measures = ("gross_world_product", "net_world_product", "net_resident_income",
                "cash_resources_after_principal", "resources_after_operations_and_principal")
    ratios = {}
    for measure in measures:
        denominator = ledger["worlds"]["earth"][measure]
        ratios[measure] = {
            "earth_denominator_U": denominator,
            "denominator_status": "positive" if denominator > 0 else "nonpositive; ratio not meaningful",
            "ratios": {world: ledger["worlds"][world][measure] / denominator if denominator > 0 else None
                       for world in ("moon", "array")},
        }
    result = {
        "schema": OUTPUT_SCHEMA, "evidence": inputs["evidence"],
        "numeraire": "U: the caller's common-price real unit; all monetary outputs are annual U unless identified as a stock or construction history.",
        "world_boundary": "earth_moon_array_are_three_peer_economies",
        "cost_case": construction_case.get("cost_case"), "water_case": construction_case.get("water_case"),
        "finance_case": finance_case, "finance_assumptions": copy.deepcopy(case),
        "applied_parameters": {"ordinary_cfc_shares": shares, "natural_depletion_shares": depletion,
                               "earth_creditor_share": ownership, "earth_legacy_debt_multiple": multiple,
                               "operation_scale": operation_scale, "asset_life_scale": asset_life_scale},
        "incremental_cfc_by_asset": asset_rows, "accounts": ledger,
        "finance_history": histories, "trade_requirements": trade,
        "required_net_exports_sum_U": math.fsum(x["required_net_exports_U"] for x in trade.values()),
        "ratios_to_earth": ratios,
        "ordinary_cfc_scope": inputs["ordinary_cfc_scope"],
        "incremental_cfc_convention": inputs["incremental_cfc_convention"],
        "endowment_convention": inputs["imported_endowment_convention"],
        "natural_depletion_status": inputs["natural_depletion_status"],
        "operations_status": operations["classification"], "trade_closure_status": inputs["trade_closure"],
        "construction_renewal_memo": {
            "rule": "0.5 times mature CFC times build years, assuming a linear capital buildout; gross replacement spending, not additional debt by default. This is not a modeled build schedule or a claim that the full lunar economy operates early.",
            "gross_replacement_spending_U": {
                world: 0.5 * ledger["worlds"][world]["total_cfc"] * build for world in ("moon", "array")
            },
        },
        "inputs_sha256": hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest(),
    }
    _finite_output(result)
    return result


def evaluate(construction_case, production_product, inputs, finance_case="legacy_retired", **overrides):
    """Convenience entry point accepting the complete productivity product."""
    result = calculate(construction_case, production_product["worlds"], inputs,
                       finance_case=finance_case, **overrides)
    result["productivity_scenario_id"] = production_product.get("scenario_id")
    result["productivity_weight_case"] = production_product.get("weight_case")
    result["earth_gdp_per_resident"] = production_product.get("earth_gdp_per_resident")
    return result
