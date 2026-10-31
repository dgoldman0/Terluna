"""Annual biomass comparison: arithmetic on the existing screen and explicit scenarios.

Run python -m research.studies.light_limitation.biomass_comparison.
This does not rerun the GCM, the canopy, or aquatic ecology.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from shared.constants import EARTH_RADIUS, JULIAN_YEAR_DAYS, MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "results/biomass_comparison.json"
INPUTS = {
    "research/studies/light_limitation/results/light_limitation.json": "terluna.research.light-limitation/1",
    "research/studies/light_limitation/biomass_comparison_inputs.json": "terluna.research.biomass-comparison-inputs/1",
}
CONTEXT = (
    "shared/scenarios/water.json",
    "illumination/surface_light/results/surface_light.json",
    "biosphere/ecology/lake.py",
    "research/studies/lake_oxygen/sources.json",
    "research/studies/lake_oxygen/run.py",
    "research/studies/resort_301k/README.md",
    "research/decisions.md",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_inputs():
    loaded = []
    for path, schema in INPUTS.items():
        data = json.loads((ROOT / path).read_text())
        if data.get("schema") != schema:
            raise ValueError(f"Unexpected schema in {path}")
        loaded.append(data)
    return loaded


def aquatic(sea_km2, lake_km2, case):
    """Partition one sea once; convert g C/m2/year times km2 into Gt C/year."""
    share = case["coast_share"]
    efficiency = case["aquatic_npp_over_gpp"]
    if not 0 <= share <= 1 or not 0 <= efficiency <= 1:
        raise ValueError("Area share and NPP/GPP must be fractions")
    coast = sea_km2 * share
    open_sea = sea_km2 - coast
    gross = case["oligotrophic_gpp_g_c_m2_day"] * JULIAN_YEAR_DAYS
    net = gross * efficiency
    parts = {"lakes": lake_km2 * net / 1e9,
             "open_sea": open_sea * net / 1e9,
             "coasts": coast * case["coast_npp_g_c_m2_year"] / 1e9}
    return {"area_km2": {"lakes": lake_km2, "open_sea": open_sea, "coasts": coast},
            "oligotrophic_gpp_g_c_m2_year": gross, "oligotrophic_npp_g_c_m2_year": net,
            "parts_gt_c_year": parts, "total_gt_c_year": sum(parts.values())}


def results():
    old, settings = load_inputs()
    moon_area = 4 * math.pi * (MOON_RADIUS / 1000) ** 2
    # Use the exact water share carried by the original food product's scenario
    # through shared/scenarios/water.json, rather than adding separately rounded
    # open-sea/coast endpoint areas from the old screen.
    water = json.loads((ROOT / "shared/scenarios/water.json").read_text())
    if water.get("schema") != "terluna.shared.water-scenario/1":
        raise ValueError("Unexpected shared water scenario")
    sea_share = water["surface_water_share"]
    sea_area = moon_area * sea_share
    lake_area = old["waters"]["lakes"]["km2"]
    land = old["land"]["totals"]
    land_pair = [land["rain_fed_gt_c"][i] + land["wet_margin_gt_c"][i] for i in (0, 1)]
    cases = {}
    for name, case in settings["cases"].items():
        waters = aquatic(sea_area, lake_area, case)
        t = case["land_endpoint"]
        terrestrial = land_pair[0] * (1 - t) + land_pair[1] * t
        total = terrestrial + waters["total_gt_c_year"]
        cases[name] = {"land_gt_c_year": terrestrial, "waters": waters,
                       "surface_total_gt_c_year": total,
                       "earth_reference_fraction": total / settings["earth_reference_gt_c_year"]}
    sunlight = settings["sunlight"]
    scales = {"completed_301k": [sunlight["completed_301k_scale"]],
              "provisional_revised_301k": [sunlight["completed_301k_scale"] * (1 + f)
                                          for f in sunlight["provisional_increase_over_completed"]]}
    light = {}
    for name, values in scales.items():
        dimming = [1 - s / sunlight["biomass_weather_scale"] for s in values]
        loss = [d * e for d in dimming for e in old["ledger"]["weighted_elasticity"]]
        light[name] = {"sunlight_scales": values,
                       "fraction_less_than_biomass_weather": sorted(dimming),
                       "fraction_npp_loss_if_only_light_changes": [min(loss), max(loss)]}
    plan = settings["planning_gt_c_year"]
    earth = settings["earth_reference_gt_c_year"]
    area_ratio = (MOON_RADIUS / EARTH_RADIUS) ** 2
    return {
        "schema": "terluna.research.biomass-comparison/1",
        "evidence": settings["evidence"],
        "reading_rule": "Gt C/year means new biomass carbon after producer respiration, over a Julian year. Case endpoints are chosen scenarios, not confidence limits. Surface totals exclude an unresolved net sky contribution. Light effects are sensitivities and are not applied as an updated climate forecast.",
        "producer": {"files": {str(Path(__file__).relative_to(ROOT)): digest(Path(__file__))},
                     "inputs": {p: digest(ROOT / p) for p in (*INPUTS, *CONTEXT)},
                     "constants": {"MOON_RADIUS": MOON_RADIUS, "EARTH_RADIUS": EARTH_RADIUS,
                                   "JULIAN_YEAR_DAYS": JULIAN_YEAR_DAYS}},
        "area": {"moon_km2": moon_area, "moon_over_earth": area_ratio,
                 "common_level_seas_km2": sea_area, "additional_rainfed_lakes_km2": lake_area,
                 "dry_land_from_rounded_bands_km2": sum(b["dry_land_km2"] for b in old["land"]["bands"].values())},
        "historical_total_gt_c_year": old["ledger"]["total_gt_c"],
        "cases": cases,
        "sky_host_gt_c_year_before_displaced_surface_production": old["sky"]["gt_c"],
        "light_only_sensitivity": light,
        "all_land_watered_surface_gt_c_year": [land["potential_gt_c"][i] + cases[n]["waters"]["total_gt_c_year"]
                                                for i, n in enumerate(("low", "high"))],
        "planning": {"gt_c_year": plan, "earth_reference_fraction": [v / earth for v in plan],
                     "mean_surface_productivity_over_earth": [v / earth / area_ratio for v in plan],
                     "dry_biomass_gt_year_at_assumed_carbon_fraction": [v / settings["dry_carbon_mass_fraction"] for v in plan],
                     "status": settings["assumptions"]["planning_range"]},
    }


if __name__ == "__main__":
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=2) + "\n")
    print(OUT.relative_to(ROOT))
