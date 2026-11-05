"""Scenario accounts for eel seed supply, feed, cooked yield, oxygen and salts.

Run from the repository root: python -m research.studies.food_systems.eel_budget
These accounts are independent of the historical global food-system products.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path

from shared.constants import JULIAN_YEAR_DAYS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "results" / "eel_budget.json"
SCHEMA = "terluna.research.eel-budget/1"


def positive(value: float, name: str, *, zero: bool = False) -> float:
    if not math.isfinite(value) or value < 0 or (value == 0 and not zero):
        raise ValueError(f"{name} must be finite and {'nonnegative' if zero else 'positive'}")
    return value


def fraction(value: float, name: str, *, zero: bool = False) -> float:
    positive(value, name, zero=zero)
    if value > 1:
        raise ValueError(f"{name} must not exceed 1")
    return value


def cohort(p: dict) -> dict:
    for k, v in p.items():
        positive(v, k)
    for k in ("growout_survival", "nursery_survival", "fertilization_fraction",
              "normal_hatch_fraction_of_fertilized", "larva_to_glass_survival",
              "cut_yield_of_live", "cook_yield_of_cut", "feed_crude_protein_fraction_dry",
              "whole_fish_protein_fraction_wet", "feed_phosphorus_fraction_dry",
              "whole_fish_phosphorus_fraction_wet"):
        fraction(p[k], k)
    if not p["glass_eel_individual_kg"] < p["nursery_exit_individual_kg"] < p["harvest_individual_kg"]:
        raise ValueError("Individual mass must increase from glass eel through harvest")
    harvest_n = p["harvest_live_kg"] / p["harvest_individual_kg"]
    nursery_n = harvest_n / p["growout_survival"]
    glass_n = nursery_n / p["nursery_survival"]
    larvae_n = glass_n / p["larva_to_glass_survival"]
    fertilized_n = larvae_n / p["normal_hatch_fraction_of_fertilized"]
    eggs_n = fertilized_n / p["fertilization_fraction"]
    glass_mass = glass_n * p["glass_eel_individual_kg"]
    nursery_mass = nursery_n * p["nursery_exit_individual_kg"]
    nursery_gain = nursery_mass - glass_mass
    growout_gain = p["harvest_live_kg"] - nursery_mass
    positive(nursery_gain, "nursery net biomass gain")
    positive(growout_gain, "growout net biomass gain")
    nursery_feed = nursery_gain * p["nursery_economic_dry_fcr"]
    growout_feed = growout_gain * p["growout_economic_dry_fcr"]
    feed = nursery_feed + growout_feed
    cut = p["harvest_live_kg"] * p["cut_yield_of_live"]
    cooked = cut * p["cook_yield_of_cut"]
    n_in_feed = feed * p["feed_crude_protein_fraction_dry"] / p["protein_to_nitrogen_factor"]
    n_in_seed = glass_mass * p["whole_fish_protein_fraction_wet"] / p["protein_to_nitrogen_factor"]
    n_in_harvest = p["harvest_live_kg"] * p["whole_fish_protein_fraction_wet"] / p["protein_to_nitrogen_factor"]
    p_in_feed = feed * p["feed_phosphorus_fraction_dry"]
    p_in_seed = glass_mass * p["whole_fish_phosphorus_fraction_wet"]
    p_in_harvest = p["harvest_live_kg"] * p["whole_fish_phosphorus_fraction_wet"]
    n_residual = n_in_feed + n_in_seed - n_in_harvest
    p_residual = p_in_feed + p_in_seed - p_in_harvest
    positive(n_residual, "nonretained nitrogen", zero=True)
    positive(p_residual, "nonretained phosphorus", zero=True)
    return {
        "harvest_live_kg": p["harvest_live_kg"],
        "harvest_individuals": harvest_n,
        "nursery_exit_individuals": nursery_n,
        "glass_eels_required": glass_n,
        "normal_hatchlings_required": larvae_n,
        "fertilized_eggs_required": fertilized_n,
        "eggs_before_fertilization_required": eggs_n,
        "incoming_glass_eel_kg": glass_mass,
        "nursery_exit_live_kg": nursery_mass,
        "nursery_dry_feed_kg": nursery_feed,
        "growout_dry_feed_kg": growout_feed,
        "post_glass_dry_feed_kg": feed,
        "prepared_raw_cut_kg": cut,
        "processing_coproducts_and_trim_kg": p["harvest_live_kg"] - cut,
        "cooking_loss_kg": cut - cooked,
        "cooked_fish_kg": cooked,
        "cooked_portions": cooked / p["portion_cooked_fish_kg"],
        "post_glass_dry_feed_per_kg_cooked": feed / cooked,
        "nutrients_kg": {
            "feed_n": n_in_feed, "seed_n": n_in_seed, "harvest_n": n_in_harvest,
            "nonretained_n": n_residual,
            "feed_p": p_in_feed, "seed_p": p_in_seed, "harvest_p": p_in_harvest,
            "nonretained_p": p_residual,
        },
    }


def oxygen(p: dict) -> dict:
    for k, v in p.items():
        positive(v, k, zero=k in {"feed_dry_kg_per_day", "tan_n_kg_per_kg_dry_feed",
                                 "nitrified_fraction_of_tan", "other_biological_o2_kg_per_kg_dry_feed"})
    fraction(p["nitrified_fraction_of_tan"], "nitrified fraction", zero=True)
    if p["capacity_factor"] < 1:
        raise ValueError("Capacity factor must be at least 1")
    delta = p["inlet_do_mg_per_l"] - p["outlet_floor_do_mg_per_l"]
    positive(delta, "inlet minus outlet dissolved oxygen")
    # 24 h per day; 10^6 mg per kg and 1000 L per cubic metre are unit conversions.
    fish = p["standing_biomass_kg"] * p["fish_respiration_mg_o2_per_kg_per_hour"] * 24 / 1e6
    tan = p["feed_dry_kg_per_day"] * p["tan_n_kg_per_kg_dry_feed"]
    nitrified = tan * p["nitrified_fraction_of_tan"]
    nitrification = nitrified * p["nitrification_o2_kg_per_kg_n"]
    other = p["feed_dry_kg_per_day"] * p["other_biological_o2_kg_per_kg_dry_feed"]
    total = fish + nitrification + other
    design = total * p["capacity_factor"]
    reserve_kg = p["connected_water_m3"] * delta / 1000
    return {
        "standing_biomass_kg": p["standing_biomass_kg"],
        "fish_o2_kg_per_day": fish,
        "tan_n_kg_per_day": tan,
        "nitrification_o2_kg_per_day": nitrification,
        "other_biological_o2_kg_per_day": other,
        "total_o2_kg_per_day": total,
        "capacity_o2_kg_per_day": design,
        "capacity_o2_kg_per_hour": design / 24,
        "nitrification_alkalinity_kg_caco3_per_day": nitrified * p["nitrification_alkalinity_kg_caco3_per_kg_n"],
        "equivalent_flow_without_reoxygenation_m3_per_hour": (design / 24) * 1000 / delta,
        "water_o2_inventory_above_floor_kg": reserve_kg,
        "ideal_inventory_minutes_at_capacity_demand": reserve_kg / (design / 24) * 60,
    }


def hatchery(p: dict, c: dict) -> dict:
    for k, v in p.items():
        positive(v, k, zero=k in {"turnaround_days", "daily_liquid_discharge_m3", "discharged_salt_recovery_fraction"})
    fraction(p["successful_batch_fraction"], "successful batch fraction")
    fraction(p["discharged_salt_recovery_fraction"], "salt recovery", zero=True)
    glass_per_year = c["glass_eels_required"] * p["annual_live_harvest_kg"] / c["harvest_live_kg"]
    duration = p["larval_cycle_days"] + p["turnaround_days"]
    yield_per_year = p["glass_eels_per_successful_batch"] * p["successful_batch_fraction"] * JULIAN_YEAR_DAYS / duration
    return {
        "annual_live_harvest_kg": p["annual_live_harvest_kg"],
        "glass_eels_required_per_year": glass_per_year,
        "expected_glass_eels_per_occupied_batch_slot_year": yield_per_year,
        "batch_slots_continuous_equivalent": glass_per_year / yield_per_year,
        "batch_slots_rounded_up": math.ceil(glass_per_year / yield_per_year),
        "note": "Mean throughput screen; integer rounding is not a reliability or redundancy design.",
    }


def salt_inventory(concentration_g_per_kg: float, p: dict) -> dict:
    positive(concentration_g_per_kg, "salt concentration", zero=True)
    for k in ("water_inventory_m3", "water_density_kg_per_m3"):
        positive(p[k], k)
    positive(p["daily_liquid_discharge_m3"], "liquid discharge", zero=True)
    fraction(p["discharged_salt_recovery_fraction"], "salt recovery", zero=True)
    kg_per_m3 = concentration_g_per_kg * p["water_density_kg_per_m3"] / 1000
    return {
        "salt_concentration_g_per_kg": concentration_g_per_kg,
        "salt_inventory_kg": p["water_inventory_m3"] * kg_per_m3,
        "makeup_salt_after_recovery_kg_per_day": p["daily_liquid_discharge_m3"] * kg_per_m3 * (1 - p["discharged_salt_recovery_fraction"]),
        "note": "Physical salt accounting only; concentration does not specify ion composition or biological suitability.",
    }


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def results() -> dict:
    p = json.loads((HERE / "eel_parameters.json").read_text())
    if p["schema"] != "terluna.research.eel-parameters/1":
        raise ValueError("Unexpected eel parameter schema")
    base = cohort(p["cohort"])
    larval_cases = []
    for survival in p["sensitivity"]["larva_to_glass_survival"]:
        q = dict(p["cohort"], larva_to_glass_survival=survival)
        larval_cases.append({"larva_to_glass_survival": survival, "budget": cohort(q)})
    processing_cases = []
    for case in p["sensitivity"]["processing_and_feed_cases"]:
        overrides = {k: v for k, v in case.items() if k != "name"}
        processing_cases.append({"name": case["name"], "budget": cohort(dict(p["cohort"], **overrides))})
    oxygen_cases = []
    for case in p["sensitivity"]["oxygen_cases"]:
        overrides = {k: v for k, v in case.items() if k != "name"}
        oxygen_cases.append({"name": case["name"], "budget": oxygen(dict(p["oxygen"], **overrides))})
    return {
        "schema": SCHEMA,
        "evidence": p["evidence"],
        "producer": {
            "module": "research.studies.food_systems.eel_budget",
            "files": {str((HERE / name).relative_to(ROOT)): digest(HERE / name)
                      for name in ("eel_budget.py", "eel_parameters.json", "eel_sources.json")},
            "environment_manifest": {"path": "research/studies/food_systems/aquatic_basis.json", "sha256": digest(HERE / "aquatic_basis.json"), "use": "Context only; no environmental carrying capacity is computed."},
            "constants": constants_used([Path(__file__)]),
        },
        "assumptions": copy.deepcopy(p),
        "reference_cohort": base,
        "reference_oxygen": oxygen(p["oxygen"]),
        "reference_hatchery_throughput": hatchery(p["hatchery"], base),
        "larval_survival_sensitivity": larval_cases,
        "processing_and_feed_sensitivity": processing_cases,
        "oxygen_sensitivity": oxygen_cases,
        "salt_inventory_sensitivity": [salt_inventory(s, p["hatchery"]) for s in p["sensitivity"]["salt_concentration_g_per_kg"]],
        "unquantified": ["Broodstock and pre-glass feed and wastes", "Hatchery oxygen and energy", "Grow-out duration and biomass trajectory", "Low-gravity oxygen transfer and electrical energy", "Local prey production and habitat carrying capacity", "Processing energy and coproduct recovery", "Survival or yield of locally spawning or fully freshwater breeding lines"],
    }


def dumps(value: dict) -> str:
    return json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n"


if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(results()))
    print(OUT.relative_to(ROOT))
