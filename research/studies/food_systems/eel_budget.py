"""Scenario accounts for eel seed supply, feed, cooked yield, oxygen and salts; per-person cases for fed eel against
heritage meat and the waters' phosphorus; and the phosphorus a wild spawning run carries downhill.

Run from the repository root: python -m research.studies.food_systems.eel_budget
The per-person cases read the older food study's product, results/food_systems.json, for its fish, gas, crop and
phosphorus factors and pin it by sha256; the older study's code and coefficients are unchanged.
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
OLDER = HERE / "results" / "food_systems.json"
SCHEMA = "terluna.research.eel-budget/2"


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


def older_factors(path: Path = OLDER) -> dict:
    """The older food study's factors the per-person cases reuse, read from its product's basis: N2O per unit of feed
    nitrogen and a recirculating farm's electricity and water per kg of live fish (its aquaponics block), the feed's
    crop energy and the diet's crop yields, N2O's Tg N per ppb, and the phosphorus the waters' life can spare. Its
    feed conversion and feed phosphorus are carried for comparison; the eel cases keep the eel's own."""
    basis = json.loads(path.read_text())["basis"]
    fish, gases, waters, diet = basis["aquaponics"], basis["gases"], basis["seas_phosphorus"], basis["diet"]
    return {
        "n2o_n_share_of_feed_n": list(fish["n2o_share_of_n"]),
        "closed_system_kwh_per_kg_live": fish["kwh_per_kg_fish"],
        "closed_system_water_l_per_kg_live": fish["water_l_per_kg_fish"],
        "older_fcr": list(fish["fcr"]),
        "older_feed_phosphorus_share": list(fish["feed_p_share"]),
        "feed_kcal_per_g": fish["feed_kcal_per_g"],
        "crop_kcal_m2_day": list(diet["crop_kcal_m2_day"]),
        "variety_factor": list(diet["variety_factor"]),
        "n2o_tg_n_per_ppb": gases["n2o_tg_n_per_ppb"],
        "waters_spare_phosphorus_kt_year": list(waters["all_waters_ceiling_kt_year"]),
        "common_level_spare_phosphorus_kt_year": list(waters["common_level_ceiling_kt_year"]),
    }


def pond_methane_kg_per_t(pond: dict) -> float:
    """Methane per tonne of fish from an open pond: its measured flux in mg CH4 per m2 per hour, held for a year, over
    the same ponds' yield in t per ha."""
    flux = positive(pond["ch4_mg_m2_h"], "pond methane flux", zero=True)
    harvest = positive(pond["yield_t_ha"], "pond yield")
    return flux * 24 * JULIAN_YEAR_DAYS * 1e4 / 1e6 / harvest     # mg per m2 a year as kg per ha, over t per ha


def per_person(people_billion: float, cooked_kg: float, c: dict, f: dict, ponds: dict, heritage: list) -> dict:
    """Fed eel for people_billion residents who each eat cooked_kg of cooked eel a year, scaled from cohort c's
    per-tonne accounts: the harvest, the feed and its crop land, the feed's nitrogen and phosphorus, the phosphorus
    the fish leave unretained against what the waters can spare, N2O, a closed system's electricity and water, and
    open ponds' methane for comparison."""
    positive(people_billion, "people (billions)")
    positive(cooked_kg, "cooked eel per person")
    live = people_billion * 1e9 * cooked_kg * c["harvest_live_kg"] / c["cooked_fish_kg"]      # kg live a year
    scale = live / c["harvest_live_kg"]
    n = c["nutrients_kg"]
    feed = c["post_glass_dry_feed_kg"] * scale
    feed_n, feed_p, lost_p = n["feed_n"] * scale, n["feed_p"] * scale, n["nonretained_p"] * scale
    # Crop land for the feed at the older study's crop energy, paired as its aquaponics block pairs them: the best
    # yield with the least variety, then the poorest yield with the most.
    land = [feed * f["feed_kcal_per_g"] * 1e3 * v / (y * JULIAN_YEAR_DAYS) / 1e6
            for v, y in ((min(f["variety_factor"]), max(f["crop_kcal_m2_day"])),
                         (max(f["variety_factor"]), min(f["crop_kcal_m2_day"])))]
    n2o = [feed_n * share / 1e6 for share in f["n2o_n_share_of_feed_n"]]        # kt N a year
    spare = f["waters_spare_phosphorus_kt_year"]
    kwh = live * f["closed_system_kwh_per_kg_live"]
    return {
        "people_billion": people_billion,
        "cooked_eel_kg_per_person_year": cooked_kg,
        "share_of_heritage_meat_allowance": [cooked_kg / max(heritage), cooked_kg / min(heritage)],
        "live_harvest_mt_year": live / 1e9,
        "post_glass_dry_feed_mt_year": feed / 1e9,
        "feed_crop_land_km2": land,
        "feed_nitrogen_kt_year": feed_n / 1e6,
        "feed_phosphorus_kt_year": feed_p / 1e6,
        "harvest_phosphorus_kt_year": n["harvest_p"] * scale / 1e6,
        "nonretained_nitrogen_kt_year": n["nonretained_n"] * scale / 1e6,
        "nonretained_phosphorus_kt_year": lost_p / 1e6,
        "nonretained_phosphorus_over_waters_spare": [lost_p / 1e6 / max(spare), lost_p / 1e6 / min(spare)],
        "n2o_kt_n_year": n2o,
        "n2o_ppb_per_year": [x / (f["n2o_tg_n_per_ppb"] * 1e3) for x in n2o],
        "closed_system_electricity_twh_year": kwh / 1e9,
        "closed_system_mean_power_gw": kwh / (24 * JULIAN_YEAR_DAYS) / 1e6,
        "closed_system_water_km3_year": live * f["closed_system_water_l_per_kg_live"] / 1e12,
        "open_pond_methane_mt_year": {name: live / 1e3 * k / 1e9 for name, k in ponds.items()},
    }


def spawning_run(p: dict, wild_survivals: list) -> dict:
    """The phosphorus a spawning run carries down per tonne of spawners, and the share the returning glass eels carry
    back for each survival from glass eel to spawner, with spawners at the cohort's harvest size and spawners and glass
    eels of the same whole-body phosphorus."""
    spawner, glass = p["harvest_individual_kg"], p["glass_eel_individual_kg"]
    cases = []
    for name, s in [("culture", p["growout_survival"] * p["nursery_survival"])] + [("wild", s) for s in wild_survivals]:
        fraction(s, "glass-eel-to-spawner survival")
        cases.append({"case": name, "glass_eel_to_spawner_survival": s,
                      "phosphorus_share_carried_back": glass / (s * spawner)})
    return {
        "spawner_individual_kg": spawner,
        "glass_eel_individual_kg": glass,
        "phosphorus_carried_down_kg_per_t_spawners": 1000 * p["whole_fish_phosphorus_fraction_wet"],
        "cases": cases,
    }


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def results() -> dict:
    p = json.loads((HERE / "eel_parameters.json").read_text())
    if p["schema"] != "terluna.research.eel-parameters/2":
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
    f = older_factors()
    pp = p["per_person"]
    ponds = {pond["name"]: pond_methane_kg_per_t(pond) for pond in pp["open_ponds"]}
    person_cases = []
    for case in processing_cases:
        for people in pp["people_billion"]:
            for cooked in pp["cooked_eel_kg_per_person_year"]:
                person_cases.append({"processing": case["name"], **per_person(
                    people, cooked, case["budget"], f, ponds, pp["heritage_meat_kg_per_person_year"])})
    return {
        "schema": SCHEMA,
        "evidence": p["evidence"],
        "producer": {
            "module": "research.studies.food_systems.eel_budget",
            "files": {str((HERE / name).relative_to(ROOT)): digest(HERE / name)
                      for name in ("eel_budget.py", "eel_parameters.json", "eel_sources.json")},
            "inputs": {str(OLDER.relative_to(ROOT)): digest(OLDER)},
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
        "older_study_factors": f,
        "open_pond_methane_kg_per_t_live": ponds,
        "per_person_cases": person_cases,
        "spawning_run_phosphorus": spawning_run(p["cohort"], p["spawning_run"]["wild_glass_eel_to_spawner_survival"]),
        "unquantified": ["Broodstock and pre-glass feed and wastes", "Hatchery oxygen and energy", "Grow-out duration and biomass trajectory", "Low-gravity oxygen transfer and electrical energy", "Local prey production and habitat carrying capacity", "Processing energy and coproduct recovery", "Survival or yield of locally spawning or fully freshwater breeding lines", "Methane from closed systems' sludge, which depends on its treatment", "The feed crops' own N2O and the land of feed made from power", "The seawater lake's temperature by depth and the eggs' exposure", "A wild population's size, and so the spawning run's phosphorus in tonnes a year"],
    }


def dumps(value: dict) -> str:
    return json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n"


if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(results()))
    print(OUT.relative_to(ROOT))
