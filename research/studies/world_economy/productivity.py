"""Conditional mature-world value-added scenarios, without a trade forecast.

All judgement inputs live in productivity_inputs.json. The calculation takes
populations and the annual Earth GDP per resident from its caller. Its price
basis is a common reference numeraire, not an observed future PPP estimate.
"""

from __future__ import annotations

import copy
import json
import math
from collections.abc import Mapping
from numbers import Real
from pathlib import Path
from typing import Any


INPUT_SCHEMA = "terluna.world-economy.productivity-inputs/1"
OUTPUT_SCHEMA = "terluna.world-economy.productivity/1"
ARRAY_ALLOCATION = "reference_weights_accounting_only_not_array_industry_estimate"
WORLD_IDS = ("earth", "moon", "array")


def _number(value: Any, name: str, *, minimum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return result


def _mapping(value: Any, name: str) -> Mapping:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _keys(value: Mapping, expected: set[str], name: str) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} must contain exactly {sorted(expected)}")


def _label(value: Any, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")


def _tax_share(value: Any, name: str) -> float:
    result = _number(value, name)
    if result >= 1:
        raise ValueError(f"{name} must be less than one")
    # Negative net taxes represent net product subsidies and are permissible.
    return result


def validate_inputs(inputs: Mapping[str, Any]) -> None:
    """Validate a complete scenario product without mutating it.

    Unknown fields are rejected so that, for example, an extra export bonus
    cannot be silently ignored or mistaken for part of the aggregate scenario.
    """
    inputs = _mapping(inputs, "inputs")
    _keys(inputs, {
        "schema", "evidence", "array_sector_allocation", "net_product_taxes",
        "sectors", "weight_cases", "scenarios",
    }, "inputs")
    if inputs["schema"] != INPUT_SCHEMA:
        raise ValueError("unsupported productivity input schema")
    evidence = _mapping(inputs["evidence"], "evidence")
    _keys(evidence, {"type", "probability_interpretation", "statement"}, "evidence")
    if evidence["type"] != "judgment_scenarios" or evidence["probability_interpretation"] != "none":
        raise ValueError("inputs must identify nonprobabilistic judgement scenarios")
    _label(evidence["statement"], "evidence.statement")
    if inputs["array_sector_allocation"] != ARRAY_ALLOCATION:
        raise ValueError("array sectors require the explicit accounting-only allocation tag")
    taxes = _mapping(inputs["net_product_taxes"], "net_product_taxes")
    _keys(taxes, {"default_share_of_gdp", "statement"}, "net_product_taxes")
    _tax_share(taxes["default_share_of_gdp"], "default tax share")
    _label(taxes["statement"], "net_product_taxes.statement")

    sectors = inputs["sectors"]
    if not isinstance(sectors, list) or not sectors:
        raise ValueError("sectors must be a nonempty list")
    sector_ids = []
    for sector in sectors:
        sector = _mapping(sector, "sector")
        _keys(sector, {"id", "label"}, "sector")
        _label(sector["id"], "sector.id")
        _label(sector["label"], "sector.label")
        sector_ids.append(sector["id"])
    if len(set(sector_ids)) != len(sector_ids):
        raise ValueError("sector identifiers must be unique")
    sector_set = set(sector_ids)

    cases = _mapping(inputs["weight_cases"], "weight_cases")
    if not cases:
        raise ValueError("weight_cases must not be empty")
    for name, weights in cases.items():
        _label(name, "weight case identifier")
        weights = _mapping(weights, f"weights.{name}")
        _keys(weights, sector_set, f"weights.{name}")
        total = math.fsum(_number(v, f"weights.{name}.{s}", minimum=0) for s, v in weights.items())
        if not math.isclose(total, 1.0, rel_tol=0, abs_tol=1e-10):
            raise ValueError(f"weights.{name} must sum to one")

    scenarios = _mapping(inputs["scenarios"], "scenarios")
    if not scenarios:
        raise ValueError("scenarios must not be empty")
    for name, scenario in scenarios.items():
        _label(name, "scenario identifier")
        scenario = _mapping(scenario, f"scenario.{name}")
        _keys(scenario, {"label", "moon_gva_per_resident_ratios", "array_gva_per_resident_ratio"}, f"scenario.{name}")
        _label(scenario["label"], f"scenario.{name}.label")
        ratios = _mapping(scenario["moon_gva_per_resident_ratios"], f"ratios.{name}")
        _keys(ratios, sector_set, f"ratios.{name}")
        for sector_id, ratio in ratios.items():
            _number(ratio, f"ratios.{name}.{sector_id}", minimum=0)
        _number(scenario["array_gva_per_resident_ratio"], f"array ratio.{name}", minimum=0)


def _unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_inputs(path: str | Path) -> dict[str, Any]:
    """Read the specified file; there is no implicit scenario or population."""
    with Path(path).open(encoding="utf-8") as stream:
        inputs = json.load(stream, object_pairs_hook=_unique_json_object)
    validate_inputs(inputs)
    return inputs


def _world(population: float, sectors: dict[str, float], tax_share: float,
           allocation: str) -> dict[str, Any]:
    gva = math.fsum(sectors.values())
    gross_product = gva / (1 - tax_share)
    sector_gross = {key: value / (1 - tax_share) for key, value in sectors.items()}
    return {
        "population": population,
        "sector_gva": sectors,
        "sector_gross_product": sector_gross,
        "sector_allocation": allocation,
        "sector_gross_product_rule": "net_product_taxes_allocated_proportionally_to_gva",
        "gva": gva,
        "net_product_tax_share_of_gdp": tax_share,
        "net_product_taxes": gross_product - gva,
        "gross_product": gross_product,
        "gva_per_resident": gva / population if population else None,
        "gross_product_per_resident": gross_product / population if population else None,
    }


def evaluate_productivity(
    inputs: Mapping[str, Any], *, earth_population: float, moon_population: float,
    array_population: float, earth_gdp_per_resident: float, scenario_id: str,
    weight_case: str = "central",
    net_product_tax_shares: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    """Return annual production in the caller's common reference currency.

    Earth GDP per resident is at market prices. Sector assumptions concern GVA
    at basic prices: subtract Earth's net product taxes before applying them,
    then bridge each territory back to GDP using its own tax share. The array's
    sector allocation is bookkeeping; only its aggregate relative GVA is a
    scenario assumption. There is no additional export term.
    """
    validate_inputs(inputs)
    populations = {
        "earth": _number(earth_population, "earth_population", minimum=0),
        "moon": _number(moon_population, "moon_population", minimum=0),
        "array": _number(array_population, "array_population", minimum=0),
    }
    earth_gdp_per_resident = _number(earth_gdp_per_resident, "earth_gdp_per_resident", minimum=0)
    if not populations["earth"] or not earth_gdp_per_resident:
        raise ValueError("Earth population and GDP per resident must be positive")
    if scenario_id not in inputs["scenarios"]:
        raise ValueError(f"unknown scenario: {scenario_id}")
    if weight_case not in inputs["weight_cases"]:
        raise ValueError(f"unknown weight case: {weight_case}")
    default_tax = inputs["net_product_taxes"]["default_share_of_gdp"]
    supplied_taxes = {} if net_product_tax_shares is None else _mapping(net_product_tax_shares, "net_product_tax_shares")
    if set(supplied_taxes) - set(WORLD_IDS):
        raise ValueError("tax shares contain an unknown territory")
    taxes = {world: _tax_share(supplied_taxes.get(world, default_tax), f"tax share.{world}") for world in WORLD_IDS}
    weights = inputs["weight_cases"][weight_case]
    scenario = inputs["scenarios"][scenario_id]
    moon_ratios = scenario["moon_gva_per_resident_ratios"]
    array_ratio = float(scenario["array_gva_per_resident_ratio"])
    earth_gva_per_resident = earth_gdp_per_resident * (1 - taxes["earth"])
    worlds = {}
    for world in WORLD_IDS:
        per_sector = {}
        for sector_id, weight in weights.items():
            ratio = 1.0 if world == "earth" else (moon_ratios[sector_id] if world == "moon" else array_ratio)
            per_sector[sector_id] = populations[world] * earth_gva_per_resident * weight * ratio
        allocation = ARRAY_ALLOCATION if world == "array" else "sector_judgment_scenario"
        worlds[world] = _world(populations[world], per_sector, taxes[world], allocation)

    earth = worlds["earth"]
    combined_gva = worlds["moon"]["gva"] + worlds["array"]["gva"]
    combined_gross = worlds["moon"]["gross_product"] + worlds["array"]["gross_product"]
    report = {
        "schema": OUTPUT_SCHEMA,
        "evidence": copy.deepcopy(inputs["evidence"]),
        "scenario_id": scenario_id,
        "scenario_label": scenario["label"],
        "weight_case": weight_case,
        "reference_price_basis": "caller_supplied_common_numeraire_not_measured_future_ppp",
        "comparison_scope": "earth_moon_array_are_three_peer_territorial_economies",
        "earth_gdp_per_resident": earth_gdp_per_resident,
        "earth_gva_per_resident": earth_gva_per_resident,
        "export_treatment": "realised_aggregate_output_no_additive_export_bonus_no_trade_equilibrium",
        "trade_flow_treatment": (
            "Interworld power, computing and other purchases belong in a separate flow ledger. "
            "Imported intermediate inputs are already deducted in realised GVA; neither their "
            "imports nor provider export sales are added to or subtracted from these totals again."
        ),
        "worlds": worlds,
        "moon_plus_array": {
            "scope": "derived_aggregate_not_an_additional_world",
            "population": populations["moon"] + populations["array"],
            "gva": combined_gva,
            "gross_product": combined_gross,
        },
        "population_ratios_to_earth": {
            "moon": populations["moon"] / populations["earth"],
            "array": populations["array"] / populations["earth"],
            "moon_plus_array": (populations["moon"] + populations["array"]) / populations["earth"],
        },
        "gva_ratios_to_earth": {
            "moon": worlds["moon"]["gva"] / earth["gva"],
            "array": worlds["array"]["gva"] / earth["gva"],
            "moon_plus_array": combined_gva / earth["gva"],
        },
        "gross_product_ratios_to_earth": {
            "moon": worlds["moon"]["gross_product"] / earth["gross_product"],
            "array": worlds["array"]["gross_product"] / earth["gross_product"],
            "moon_plus_array": combined_gross / earth["gross_product"],
        },
        "assumed_gva_per_resident_ratios_to_earth": {
            "moon": math.fsum(weights[sector] * moon_ratios[sector] for sector in weights),
            "array": array_ratio,
        },
    }
    _finite_output(report)
    return report


def _finite_output(value: Any) -> None:
    if isinstance(value, Mapping):
        for item in value.values():
            _finite_output(item)
    elif isinstance(value, Real) and not isinstance(value, bool):
        if not math.isfinite(value):
            raise ValueError("production calculation overflowed; reduce the input scale")
