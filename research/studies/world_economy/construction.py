"""Conditional construction expenditure in common real units; not a GDP model.

Imported commodities include upstream capital service. Separately priced assets
are final installed assets. A Moon beneficiary does not imply lunar territory,
ownership or CFC. Debt, learning and post-completion depreciation are elsewhere.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS


from . import inventory
from shared.provenance import constants_used

INPUT_PATH = Path(__file__).with_name("construction_inputs.json")



def load_inputs():
    data = json.loads(INPUT_PATH.read_text())
    if data.get("schema") != "terluna.research.world-economy-construction-inputs/1":
        raise ValueError("Unexpected construction input schema")
    _, selected = inventory.load_inputs()
    p = data["physical"]
    first = selected["resources/results/first_screen.json"]["values"]
    exact = {"air_kg":first["inventories.air_kg"],
             "nitrogen_kg":first["inventories.n2_kg"],
             "oxygen_kg":first["inventories.o2_kg"],
             "sea_and_lake_water_kg":first["inventories.water"]["total_kg"][0]}
    for key, expected in exact.items():
        if p[key] != expected:
            raise ValueError(f"Construction physical input differs from bound product: {key}")
    pop = selected["research/studies/distributions/results/distributions.json"]["values"]["distributions.C'.people_billion"]
    for field, world in (("earth_population","earth"),("lunar_population","moon"),
                         ("aerial_population","moon_aloft"),("array_population","array")):
        if data[field] != pop[world]*1e9:
            raise ValueError(f"Population differs from accepted C-prime: {field}")
    fleet = selected["research/studies/solar_shield_array/results/array_heat.json"]["values"]["fleet.optical_mass_Gt"]
    if not fleet[0]*1e12 <= p["shield_fleet_kg"] <= fleet[1]*1e12:
        raise ValueError("Rounded optical fleet is outside source range")
    if data["crust_water_kg"]["deep_crust"] != first["inventories.water"]["crust_kg"][1]:
        raise ValueError("Deep crust quantity differs from bound source")
    for source in data["physical_sources"]:
        if source["sha256"] != selected[source["path"]]["sha256"]:
            raise ValueError("Construction source identity differs from inventory")
    return data

def _positive(value, name, *, zero=False):
    value = float(value)
    if not math.isfinite(value) or value < 0 or (value == 0 and not zero):
        raise ValueError(f"{name} must be finite and {'nonnegative' if zero else 'positive'}")
    return value


def _fraction(value, name):
    value = _positive(value, name)
    if value > 1:
        raise ValueError(f"{name} must not exceed one")
    return value


def import_account(inputs, case, crust_water_kg):
    """Joint ore accounting and the delivered-import bill before financing."""
    p = inputs["physical"]
    w = _positive(p["sea_and_lake_water_kg"], "water") + _positive(crust_water_kg, "crust", zero=True)
    air = _positive(p["air_kg"], "air")
    n = _positive(p["nitrogen_kg"], "nitrogen")
    o = _positive(p["oxygen_kg"], "oxygen")
    if n + o > air * (1 + 1e-12):
        raise ValueError("Constituent gases exceed total air")
    mass = w + air
    water_feed = w + o * _positive(p["water_feed_kg_per_kg_oxygen"], "oxygen water feed")
    water_feed /= _fraction(p["source_water_mass_fraction"], "source water fraction")
    nitrogen_feed = n / _fraction(case["nitrogen_mass_fraction_in_source"], "nitrogen grade")
    feed = max(water_feed, nitrogen_feed)
    price = _positive(case["energy_U_per_MWh"], "energy price", zero=True) / 3600.0
    energy_mj = {
        "source_ore": feed * _positive(case["ore_processing_MJ_per_kg_feed"], "ore energy", zero=True),
        "source_chemistry": n * _positive(case["nitrogen_chemistry_MJ_per_kg_N"], "N energy", zero=True)
        + o * _positive(case["oxygen_chemistry_MJ_per_kg_O"], "O energy", zero=True),
        "departure_propulsion": mass * _positive(case["departure_propulsion_MJ_per_kg_delivered"], "propulsion energy", zero=True),
        "capture_controls": mass * _positive(case["capture_controls_MJ_per_kg_delivered"], "capture energy", zero=True),
        "placement": mass * _positive(case["placement_MJ_per_kg_delivered"], "placement energy", zero=True),
    }
    service = {
        "source_ore": feed * _positive(case["ore_hardware_service_U_per_kg_feed"], "ore service", zero=True),
        "source_chemistry": air * _positive(case["chemical_hardware_service_U_per_kg_air"], "chemistry service", zero=True),
        "departure_propulsion": mass * _positive(case["departure_hardware_service_U_per_kg_delivered"], "departure service", zero=True),
        "capture_controls": mass * _positive(case["capture_hardware_service_U_per_kg_delivered"], "capture service", zero=True),
        "placement": mass * _positive(case["placement_hardware_service_U_per_kg_delivered"], "placement service", zero=True),
    }
    costs = {k: {"purchased_energy_U": v * price, "hardware_service_U": service[k],
                 "total_U": v * price + service[k]} for k, v in energy_mj.items()}
    base = math.fsum(v["total_U"] for v in costs.values())
    contingency = base * _positive(case["import_contingency_fraction"], "contingency", zero=True)
    years = _positive(inputs["build_years"], "build years")
    seconds = years * JULIAN_DAY * JULIAN_YEAR_DAYS
    diagnostic = inputs["generation_cost_diagnostic"]
    life = _positive(diagnostic["life_years"], "generation diagnostic life")
    rate = _positive(diagnostic["annual_real_capital_rate"], "generation diagnostic rate", zero=True)
    recovery_factor = 1 / life if rate == 0 else rate / -math.expm1(-life * math.log1p(rate))
    mwh_per_continuous_w_year = JULIAN_DAY * JULIAN_YEAR_DAYS / 3.6e9
    generation_ceiling = case["energy_U_per_MWh"] * mwh_per_continuous_w_year / recovery_factor
    power_pw = math.fsum(energy_mj.values()) * 1e6 / seconds / 1e15
    return {
        "beneficiary": "moon", "asset_territory": "planetary_environment_with_upstream_services",
        "delivered_water_kg": w, "delivered_air_kg": air, "delivered_mass_kg": mass,
        "source_feed_kg": feed, "water_feed_requirement_kg": water_feed,
        "nitrogen_feed_requirement_kg": nitrogen_feed,
        "components": costs, "purchased_energy_MJ": energy_mj,
        "purchased_energy_mean_PW": power_pw,
        "construction_power_over_mature_moon_transmitting_power": power_pw * 1000 / _positive(case["orbital_power_capacity_TW"], "mature power"),
        "generation_capital_ceiling_U_per_continuous_W": generation_ceiling,
        "generation_diagnostic": {**diagnostic,
            "assumption": "100% availability, zero operation and maintenance cost, no fuel, all electricity revenue available for capital recovery; a ceiling, not an extra construction charge"},
        "incoming_braking_energy_charged_as_purchased_MJ": 0.0,
        "energy_recovery_revenue_U": 0.0, "contingency_U": contingency,
        "total_U": base + contingency, "delivered_average_U_per_kg": (base + contingency) / mass,
    }


def installed_assets(inputs, case):
    """Disjoint installed-stock screen; no upstream delivery plants added."""
    p = inputs["physical"]
    moon = _positive(inputs["lunar_population"], "Moon population")
    sky = _positive(inputs["aerial_population"], "aerial population", zero=True)
    array = _positive(inputs["array_population"], "array population", zero=True)
    rows = []

    def add(name, quantity, price, unit, territory, beneficiary, status="scenario_price_on_physical_screen"):
        quantity = _positive(quantity, name + " quantity", zero=True)
        price = _positive(price, name + " price", zero=True)
        rows.append(dict(element=name, quantity=quantity, quantity_unit=unit, unit_cost_U=price,
                         total_U=quantity * price, asset_territory=territory,
                         beneficiary=beneficiary, evidence=status))

    physical_proxy = (case["ordinary_stock_kg_per_lunar_resident"]
                      * case["ordinary_stock_U_per_kg"])
    ordinary_per_person = _positive(case["ordinary_installed_capital_U_per_resident"], "ordinary capital allowance")
    if ordinary_per_person < physical_proxy:
        raise ValueError("Ordinary capital allowance is smaller than its included physical proxy")
    add("ordinary_lunar_capital", moon, ordinary_per_person, "resident", "moon", "moon",
        "common_per_resident_total_capital_allowance_excluding_the_separate_rows")
    add("ordinary_array_capital", array, ordinary_per_person, "resident", "array", "array",
        "same_ordinary_capital_valuation_as_Earth_and_Moon_includes_habitat_structure_and_normal_computing")
    add("aerial_carriers_and_control", sky * p["aerial_carrier_allowance_kg_per_person"],
        case["aerial_carrier_U_per_kg"], "kg", "moon", "moon")
    add("optical_shield_fleet", p["shield_fleet_kg"], case["shield_fabrication_U_per_kg"], "kg", "array", "moon")
    add("orbital_power_and_beaming", case["orbital_power_capacity_TW"] * 1e12,
        case["orbital_power_U_per_W"], "W", "array", "moon")
    add("ground_receivers", case["ground_receiver_km2"] * 1e6,
        case["ground_receiver_U_per_m2"], "m2", "moon", "moon")
    add("critical_energy_reserves", case["reserve_TWh"] * 1e12,
        case["reserve_U_per_Wh"], "Wh", "moon", "moon")
    add("environment_establishment_allowance", 1, case["environment_establishment_allowance_U"],
        "allowance", "moon", "moon", "unquantified_allowance")
    add("magnetic_protection_allowance", 1, case["magnetic_protection_allowance_U"],
        "allowance", "array", "moon", "unquantified_allowance")
    add("array_bulk_habitat_shield", array * p["array_shield_kg_per_person"],
        case["array_bulk_shield_U_per_kg"], "kg", "array", "array")
    return rows


def calculate(inputs=None):
    """Return 3 cost cases x 3 water inventories and independent income scales."""
    if inputs is None:
        inputs = load_inputs()
    inputs = copy.deepcopy(inputs)
    years = _positive(inputs["build_years"], "build years")
    earth = _positive(inputs["earth_population"], "Earth population")
    rows = []
    for name, case in inputs["cost_cases"].items():
        assets = installed_assets(inputs, case)
        moon_assets = math.fsum(x["total_U"] for x in assets if x["asset_territory"] == "moon")
        array_assets = math.fsum(x["total_U"] for x in assets if x["asset_territory"] == "array")
        ordinary_per_person = case["ordinary_installed_capital_U_per_resident"]
        for water_case, crust in inputs["crust_water_kg"].items():
            imports = import_account(inputs, case, crust)
            moon_total = imports["total_U"] + moon_assets
            combined = moon_total + array_assets
            ratios = {}
            for income in inputs["earth_value_added_per_person"]:
                output = earth * _positive(income, "Earth value added per person")
                ratios[str(income)] = {
                    "earth_annual_output_U": output,
                    "imports_in_earth_output_years": imports["total_U"] / output,
                    "moon_capital_in_earth_output_years": moon_total / output,
                    "array_capital_in_earth_output_years": array_assets / output,
                    "earth_ordinary_capital_in_earth_output_years": earth * ordinary_per_person / output,
                    "combined_in_earth_output_years": combined / output,
                    "annual_build_share_of_earth_output": combined / years / output,
                }
            rows.append(dict(cost_case=name, water_case=water_case, imports=imports,
                             installed_assets=assets, moon_territorial_capital_U=moon_total,
                             array_territorial_capital_U=array_assets, combined_total_U=combined,
                             earth_ordinary_capital_benchmark_U=earth * ordinary_per_person,
                             included_ordinary_physical_proxy=dict(
                                 mass_kg_per_resident=case["ordinary_stock_kg_per_lunar_resident"],
                                 price_U_per_kg=case["ordinary_stock_U_per_kg"],
                                 value_U_per_resident=case["ordinary_stock_kg_per_lunar_resident"] * case["ordinary_stock_U_per_kg"],
                                 total_allowance_U_per_resident=ordinary_per_person,
                                 added_to_total=False),
                             included_array_structure_comparison=dict(
                                 mass_kg=inputs["array_population"] * inputs["physical"]["array_structure_kg_per_person"],
                                 separately_valued_U=inputs["array_population"] * inputs["physical"]["array_structure_kg_per_person"] * case["array_structure_U_per_kg"],
                                 added_to_total=False),
                             beneficiary_memo_U={
                                 "moon": imports["total_U"] + math.fsum(a["total_U"] for a in assets if a["beneficiary"] == "moon"),
                                 "array": math.fsum(a["total_U"] for a in assets if a["beneficiary"] == "array")},
                             mean_annual_build_U=combined / years, earth_comparisons=ratios))
    return {"schema": "terluna.research.world-economy-construction/1",
            "evidence": inputs["evidence"], "numeraire": inputs["numeraire"],
            "build_years": years, "rules": inputs["rules"], "cases": rows,
            "inputs_sha256": hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()}



def results():
    product = calculate()
    paths = [Path(__file__), INPUT_PATH, inventory.HERE / "inventory.py", inventory.HERE / "inputs.json"]
    product["producer"] = dict(study="world_economy", files={
        str(p.relative_to(inventory.ROOT)):inventory.digest(p) for p in paths},
        constants=constants_used([__file__]),
        inputs=json.loads((inventory.HERE / "inputs.json").read_text())["products"])
    product["reading_rule"] = "U is a common transferable real unit. All costs are total construction expenditures before purchaser financing, not annual GDP or CFC. Territories own their installed assets; beneficiaries are a separate memo. Earth capital is a comparator excluded from the Moon-plus-array total. Prices and source-site engineering remain assumptions."
    return product


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=inventory.ROOT / "research/runs/world_economy/construction.json")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2)+"\n")
    print(args.output)

if __name__ == "__main__":
    main()
