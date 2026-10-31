"""Aerial biomass production with one shared photon column and explicit resource demands.

python -m research.studies.aerial_productivity.run
Reads products only. No climate, ray-tracing or organism simulation is launched.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

from shared.constants import JULIAN_YEAR_DAYS, MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "results/aerial_productivity.json"
INPUTS = {
    "research/studies/sky_ecology/results/sky_ecology.json": "terluna.research.sky-ecology/1",
    "research/studies/aerophytes/results/aerophytes.json": "terluna.research.aerophytes/1",
    "research/studies/light_limitation/results/light_limitation.json": "terluna.research.light-limitation/1",
    "research/studies/light_limitation/results/biomass_comparison.json": "terluna.research.biomass-comparison/1",
    "research/studies/light_limitation/biomass_comparison_inputs.json": "terluna.research.biomass-comparison-inputs/1",
    "research/studies/aerial_ecology/results/aerial_ecology.json": "terluna.research.aerial-ecology/1",
    "research/studies/aerial_productivity/parameters.json": "terluna.research.aerial-productivity-parameters/1",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load():
    out = []
    for path, schema in INPUTS.items():
        data = json.loads((ROOT / path).read_text())
        if data.get("schema") != schema:
            raise ValueError(f"Unexpected schema in {path}: {data.get('schema')}")
        out.append(data)
    return out


def photon_column(loading, habitat_fraction=1.0):
    if loading < 0 or not 0 < habitat_fraction <= 1:
        raise ValueError("Nonnegative area loading and a habitat fraction in (0,1] required")
    local = loading / habitat_fraction
    intercepted = -habitat_fraction * math.expm1(-local)
    return {"area_loading": loading, "habitat_fraction": habitat_fraction,
            "intercepted_fraction": intercepted,
            "local_surface_transmission": math.exp(-local),
            "mean_organism_light_fraction": intercepted / loading if loading else 1.0,
            "bottom_plane_mean_light_fraction": math.exp(-local)}


def body_budget(row, light_fraction, hazard, extra_fraction, yield_c):
    """Keep maintenance and replacement paid as incident light changes."""
    if not (0 <= extra_fraction <= 1 and 0 <= hazard <= 1 and light_fraction >= 0):
        raise ValueError("Invalid budget fraction")
    gross = row["gpp_c"] * light_fraction
    fixed = (row["photolytic_gas_not_fixed_c"] + row["maintenance_c"]
             + row["night_gas_substrate_c"] + row["renewal_biomass_c"]
             + row["renewal_respiration_c"] + row["food_to_tenants_c"])
    replacement = hazard * row["construction_c_kg_m2"]
    left = gross - fixed - replacement
    h2_feasible = gross >= row["photolytic_gas_not_fixed_c"]
    feasible = left >= -1e-12 and h2_feasible
    minimum = row["renewal_biomass_c"] + hazard * row["standing_host_c_kg_m2"] + row["food_to_tenants_c"]
    extra = extra_fraction * max(0.0, left) * yield_c
    return {"cls": row["cls"], "diameter_m": row["diameter_m"],
            "gpp_c": gross, "fixed_assimilate_cost_c": fixed, "replacement_cost_c": replacement,
            "remaining_assimilate_c": left, "feasible_at_mean_light": feasible,
            "photolytic_h2_area_time_feasible": h2_feasible,
            "accounted_organic_production_c": minimum if feasible else None,
            "extra_tissue_c": extra if feasible else None,
            "organic_production_with_extra_c": minimum + extra if feasible else None,
            "minimum_light_fraction": (fixed + replacement) / row["gpp_c"],
            "assimilate_balance_residual_c": gross - fixed - replacement - left}


def class_light(column, shares, order, arrangement):
    """The same incoming photon stream is assigned once across all class layers."""
    if set(order) != set(shares) or not math.isclose(sum(shares.values()), 1):
        raise ValueError("One normalized share and one vertical place per class required")
    if arrangement == "mixed":
        return {cls: {"mean": column["mean_organism_light_fraction"],
                      "bottom": column["bottom_plane_mean_light_fraction"]} for cls in shares}
    if arrangement != "stratified":
        raise ValueError("Unknown arrangement")
    h, c = column["habitat_fraction"], column["area_loading"]
    above, out = 0.0, {}
    for cls in order:
        own = shares[cls] * c
        # expm1 avoids cancellation at sparse cover.
        mean = math.exp(-above / h) * (-h * math.expm1(-own / h) / own) if own else 1.0
        out[cls] = {"mean": mean, "bottom": math.exp(-(above + own) / h)}
        above += own
    return out


def surface_loss(budgets, potential_rainfed, column, elasticities):
    """Two alternative surface responses; a water limit cannot exceed light potential."""
    h, t = column["habitat_fraction"], column["local_surface_transmission"]
    e = elasticities
    rain_loss = max(0.0, budgets["rainfed"] - potential_rainfed * t ** e["rainfed_light_potential"])
    local = (rain_loss + budgets["wet_margins"] * (1 - t ** e["wet_margins"])
             + budgets["coasts"] * (1 - t ** e["coasts"])
             + budgets["oligotrophic_water"] * (1 - t ** e["oligotrophic_water"]))
    return {"homogenized_gt_c_year": h * local,
            "persistent_opaque_gt_c_year": column["intercepted_fraction"] * sum(budgets.values()),
            "rainfed_loss_gt_c_year": h * rain_loss}


def microbial_screen(area_m2, carbon_fraction, settings):
    """A fixed optical depth fixes column biomass regardless of how deep it is spread."""
    s = settings
    synthesis_rate = s["favourable_cycle_fraction"] * math.log(2) / s["doubling_days"]
    optimistic = s["optimistic_engineered_case"]
    optimistic_rate = optimistic["favourable_cycle_fraction"] * math.log(2) / optimistic["doubling_days"]
    rows = []
    for diameter in s["particle_diameter_m"]:
        extinction = 3 * s["extinction_efficiency"] / (2 * s["particle_dry_density_kg_m3"] * diameter)
        for tau in s["optical_depths"]:
            column_dry = tau / extinction
            stock_c = column_dry * carbon_fraction * area_m2 / 1e12
            rows.append({"diameter_m": diameter, "optical_depth": tau,
                         "dry_column_kg_m2": column_dry,
                         "dry_concentration_mg_m3": column_dry / s["column_depth_m"] * 1e6,
                         "stock_gt_c": stock_c,
                         "cell_synthesis_gt_c_year": stock_c * synthesis_rate * JULIAN_YEAR_DAYS,
                         "optimistic_cell_synthesis_gt_c_year": stock_c * optimistic_rate * JULIAN_YEAR_DAYS,
                         "growth_minus_removal_per_day": {str(r): synthesis_rate - 1 / r - s["other_mortality_per_day"]
                                                           for r in s["removal_days"]},
                         "light_removed_fraction_geometric": -math.expm1(-tau)})
    return rows


def results():
    sky, aero, old, surface, previous, aerial, p = load()
    area = 4 * math.pi * MOON_RADIUS ** 2
    shares = sky["phosphorus"]["design_guesses"]["class_shares_of_cover"]
    representatives = []
    for cls, diameter in p["representatives"].items():
        representatives.append(next(r for r in sky["food_web"]["production"]
                                    if r["cls"] == cls and r["diameter_m"] == diameter
                                    and r["gas_route"] == p["gas_route"]
                                    and r["gpp_c_kg_m2_year"] == p["reference_gpp_kg_c_m2_year"]))
    if not all(r["all_gates"] for r in representatives) or not math.isclose(sum(shares.values()), 1):
        raise ValueError("Reference organisms must pass their stored gates and class shares must sum to one")
    light_scale = sum(surface["light_only_sensitivity"]["provisional_revised_301k"]["sunlight_scales"]) / 2
    light_ratio = light_scale / previous["sunlight"]["biomass_weather_scale"]
    base = surface["cases"]["illustrative_middle"]
    wet = sum(old["land"]["totals"]["wet_margin_gt_c"]) / 2
    water = base["waters"]["parts_gt_c_year"]
    budgets = {"rainfed": base["land_gt_c_year"] - wet, "wet_margins": wet,
               "coasts": water["coasts"], "oligotrophic_water": water["lakes"] + water["open_sea"]}
    potential = sum(old["land"]["totals"]["potential_gt_c"]) / 2 - wet
    body_keys = {"sky_reef": "sky_reef_60m", "round_giant": "round_giant_1km", "lesser": "lesser_100m"}
    p_need_g_m2 = []
    for edge in ("low", "high"):
        weighted = 0.0
        for cls, body in body_keys.items():
            candidates = [r for r in sky["phosphorus"]["loops"]
                          if r["state"] == "mature" and r["retention_case"] == "earth_forest"
                          and r["edge"] == edge and r["body"] == body]
            choose = min if edge == "low" else max
            weighted += shares[cls] * choose(r["need_g_m2_year"] for r in candidates)
        p_need_g_m2.append(weighted)
    water_traits = {r["name"]: r for r in aero["water"]["needs"]}
    rows = []
    for arrangement, habitat in itertools.product(p["arrangements"], p["habitat_fractions"]):
        for loading in p["area_loadings"]:
            column = photon_column(loading, habitat)
            doses = class_light(column, shares, p["vertical_order"], arrangement)
            losses = surface_loss(budgets, potential, column, p["surface_elasticities"])
            bodies = [body_budget(r, light_ratio * doses[r["cls"]]["mean"],
                                  p["adult_hazard_per_year"], max(p["extra_assimilate_to_tissue_fractions"]),
                                  p["new_tissue_carbon_yield"]) for r in representatives]
            bottom = [body_budget(r, light_ratio * doses[r["cls"]]["bottom"],
                                  p["adult_hazard_per_year"], 0, p["new_tissue_carbon_yield"])
                      for r in representatives]
            feasible = all(r["feasible_at_mean_light"] for r in bodies)
            bottom_feasible = all(r["feasible_at_mean_light"] for r in bottom)
            creditable = feasible and bottom_feasible
            factor = area * loading / 1e12
            mean_production = sum(shares[r["cls"]] * r["accounted_organic_production_c"] for r in bodies) * factor if feasible else None
            production = mean_production if creditable else None
            extra = sum(shares[r["cls"]] * r["extra_tissue_c"] for r in bodies) * factor if creditable else None
            water_volume = {name: water_traits[name]["water_kg_m2_year"] * area * loading / 1e12
                            for name in ("reference", "concentrating_cool_leaves", "concentrating_humid_layer")}
            rows.append({"arrangement": arrangement, "column": column, "class_light_fractions": doses,
                         "projected_area_km2": area * loading / 1e6,
                         "bodies": bodies, "bottom_plane_budgets": bottom,
                         "whole_mixture_feasible_at_mean_light": feasible,
                         "all_classes_feasible_at_bottom_plane_mean": bottom_feasible,
                         "passes_mean_and_bottom_light_budget_checks": creditable,
                         "mean_budget_only_diagnostic_gt_c_year": mean_production,
                         "sky_accounted_gt_c_year": production,
                         "additional_tissue_allocation_sensitivity_gt_c_year": extra,
                         "surface_loss": losses,
                         "additional_world_gt_c_year": {k: production - v if production is not None else None
                                                         for k, v in losses.items() if k != "rainfed_loss_gt_c_year"},
                         "mature_p_resupply_before_fallout_t_year": [n * area * loading / 1e6 for n in p_need_g_m2],
                         "ideal_interception_of_shared_p_fallout_t_year": [v * column["intercepted_fraction"] / 1000
                                                                          for v in aerial["baseline"]["p_fallout_analogue_kg_year"]],
                         "water_throughput_gt_year": water_volume})
    return {"schema": "terluna.research.aerial-productivity/1",
            "evidence": "Conditional column, carbon and resource arithmetic on inherited engineered-organism products. Mutual shading is vertical-beam only: inherited ground-canopy diffuse light is included, but angular transport through the aerial community is unresolved. Occupancy is prescribed, not predicted; no stable rich sky or new 301 K environment has been established.",
            "reading_rule": "Annual production in Gt C per Julian year. Accounted sky organic production is renewed/replacement host tissue plus sugar export, not a complete NPP account. Whole-world increment subtracts one of two alternative surface responses. Credited production requires BOTH class-mean and class-bottom light-budget checks; failure returns null. A positive mean-budget-only diagnostic is not credited. These checks do not establish individual viability under intermittent opaque shadows. Water is throughput in Gt/year; P is return needed, not proven supply.",
            "producer": {"files": {str(Path(__file__).relative_to(ROOT)): sha(Path(__file__))},
                         "inputs": {**{path: sha(ROOT / path) for path in INPUTS},
                                    "research/studies/aerial_productivity/sources.json": sha(HERE / "sources.json")},
                         "constants": {"MOON_RADIUS": MOON_RADIUS, "JULIAN_YEAR_DAYS": JULIAN_YEAR_DAYS}},
            "surface_reference_gt_c_year": sum(budgets.values()), "surface_components_gt_c_year": budgets,
            "rainfed_light_potential_gt_c_year": potential,
            "relative_provisional_sunlight": light_ratio,
            "class_shares": shares, "scenarios": rows,
            "microbial_optical_scenarios": microbial_screen(area, previous["dry_carbon_mass_fraction"], p["microbes"]),
            "vertical_habitat": {"source": "older aerophyte environment; new 301 K fields pending",
                                 "air": [{k: r[k] for k in ("height_km", "temperature_c", "pressure_atm")} for r in aero["environment"]["air"]],
                                 "collecting_shell_area_ratios": {str(z): ((MOON_RADIUS + z * 1000) / MOON_RADIUS) ** 2 for z in (10, 20, 30)}},
            "unresolved": ["Attainable coverage, demographic equilibrium and shadow duration along routes",
                           "Nitrogen fixation energy and trace-element return; additional nutrients for new tissue allocation",
                           "Water collection and storage through dry gaps, especially small reef modules",
                           "Dense-biome clouds, albedo, rain and thermal feedback on the final 301 K climate",
                           "Measured combined physiological, gas-retention and structural traits"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results(), indent=2) + "\n")
    print(args.output)
