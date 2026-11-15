"""Bounded engineering screens for the enclosed off-Earth programme.

Read the admitted staged product; never rerun a founding fleet per resident.
New equations test individual requirements, not integrated mission feasibility.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

from scipy.optimize import brentq
from scipy.special import j0, j1

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
METRICS = (
    "ground_firm_W", "lunar_process_W", "source_process_W", "source_heat_W",
    "depot_W", "source_feed_kg_year", "source_water_kg_year",
    "source_nitrogen_kg_year", "source_cargo_kg_year", "source_tug_count",
    "source_retained_tails_kg_year", "source_tug_component_import_kg_year",
    "LEO_kg_year", "lunar_raw_kg_year", "lunar_bulk_kg_year",
    "lunar_recycled_bulk_kg_year", "lunar_metal_kg_year",
    "lunar_virgin_metal_kg_year",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_pinned(root: Path, entry: dict) -> dict:
    path = root / entry["path"]
    actual = sha256(path)
    if actual != entry["sha256"]:
        raise ValueError(f"Input hash mismatch: {entry['path']}")
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as stream:
            return json.load(stream)
    return json.loads(path.read_text())


def encircled_power(x: float) -> float:
    """Uniform circular aperture, ideal focus; x=pi D r/(lambda R)."""
    if x < 0:
        raise ValueError("Dimensionless radius must be nonnegative")
    return float(1 - j0(x) ** 2 - j1(x) ** 2)


def capture_radius(capture: float, wavelength: float, distance: float, aperture: float) -> float:
    if not 0 < capture < 1 or min(wavelength, distance, aperture) <= 0:
        raise ValueError("Invalid focused-aperture inputs")
    # Monotonic encircled power. Expand bracket for high requested capture.
    hi = 4.0
    while encircled_power(hi) < capture:
        hi *= 2
    x = brentq(lambda z: encircled_power(z) - capture, 0, hi)
    return x * wavelength * distance / (math.pi * aperture)


def radiator_area(heat: float, temperature: float, emissivity: float,
                  sigma: float, margin: float, effective_fraction: float = 1) -> float:
    """Emitting area against cold sky; effective_fraction is an assumed derating."""
    if heat < 0 or temperature <= 0 or not 0 < emissivity <= 1 or not 0 < effective_fraction <= 1:
        raise ValueError("Invalid radiator inputs")
    return heat * margin / (emissivity * sigma * temperature ** 4 * effective_fraction)


def loss_probability(trips: float, probability: float) -> float:
    """Independent identical losses, conditional illustration only."""
    if trips < 0 or not 0 <= probability < 1:
        raise ValueError("Invalid traffic inputs")
    return -math.expm1(trips * math.log1p(-probability))


def refurbishment_purchase_ratio(q: float, parts: float, restored_life: float) -> float:
    """Renewal-reward mass purchase per service-year / all-new reference.

    At every retirement, probability q rebuilds the unit using `parts` of its
    mass and earns restored_life of a new lifetime. Otherwise a new unit is
    purchased. Downtime, workshops, labour and correlated failure are excluded.
    """
    if not 0 <= q < 1 or not 0 <= parts <= 1 or restored_life < 0:
        raise ValueError("Invalid refurbishment inputs")
    return (1 - q + q * parts) / (1 - q + q * restored_life)


def extract_requirements(product: dict, checkpoints: list[int]) -> list[dict]:
    cases = []
    for case in product["cases"]:
        physical = case["physical_case"]
        rows = physical["annual"]
        points = []
        for row in rows:
            if row["year"] in checkpoints:
                point = {"year": row["year"], "population": row["population"],
                         "feasible": row["feasible"],
                         "resident_workforce_ratio": row["capacity_summary"]["resident_workforce"],
                         "source_departures_year": row["physical"]["source"]["trips"]}
                point.update({metric: row["physical"][metric] for metric in METRICS})
                points.append(point)
        peaks = {}
        for metric in METRICS:
            peak = max(rows, key=lambda r: r["physical"][metric])
            peaks[metric] = {"year": peak["year"], "value": peak["physical"][metric]}
        workforce = max(rows, key=lambda r: r["capacity_summary"]["resident_workforce"])
        cases.append({"name": case["name"], "variant": case["variant"],
                      "admitted_completed_path": case["admitted_completed_path"],
                      "first_capacity_shortfall": physical["first_capacity_shortfall"],
                      "checkpoints": points, "operating_peaks": peaks,
                      "workforce_peak": {"year": workforce["year"],
                                         "ratio": workforce["capacity_summary"]["resident_workforce"]},
                      "seed": physical["seed"]})
    return cases


def calculate(config: dict, data: dict) -> dict:
    e, c = data["engineering"], data["constants"]
    p, t, z = (e["industrial"][key] for key in ("power", "transport", "source"))
    sensitivity = config["sensitivity"]
    seconds_day = c["time"]["julian_day_s"]
    days_year = c["time"]["julian_year_days"]
    hours_year = days_year * seconds_day / 3600
    sigma = c["physics"]["stefan_boltzmann_W_m2_K4"]
    margin = e["industrial"]["economics"]["capacity_margin"]
    cases = extract_requirements(data["industry"], config["checkpoints"])
    reference = next(case for case in cases if case["name"] == config["reference_case"])
    terminal = next(row for row in reference["checkpoints"] if row["year"] == 200)
    wavelength = c["physics"]["speed_of_light_m_s"] / p["microwave_Hz"]
    beam = []
    for distance in sensitivity["slant_range_m"]:
        radius = capture_radius(p["beam_capture"], wavelength, distance, p["aperture_m"])
        area = math.pi * radius ** 2
        sites = []
        for mw in sensitivity["beam_delivered_site_MW"]:
            # The inherited 200 W/m2 is delivered electric flux, as defined in
            # initial_industrial_system/run.py and its engineering write-up.
            flux_area = mw * 1e6 / p["receiver_mean_flux_W_m2"]
            x = math.pi * p["aperture_m"] * math.sqrt(flux_area / math.pi) / (wavelength * distance)
            sites.append({"delivered_MW": mw, "flux_sized_area_m2": flux_area,
                          "capture_at_flux_sized_radius": encircled_power(x),
                          "mean_delivered_electric_W_m2_at_80pct_capture": mw * 1e6 / area,
                          "mean_received_RF_W_m2_at_80pct_capture": mw * 1e6 / p["receiver_efficiency"] / area})
        beam.append({"distance_m": distance, "capture_radius_m": radius,
                     "capture_area_m2": area,
                     "delivered_power_at_reference_flux_W": area * p["receiver_mean_flux_W_m2"],
                     "site_comparisons": sites})
    radius_moon = c["moon"]["radius_m"]
    orbit_radius = radius_moon + p["orbital_altitude_m"]
    period_s = 2 * math.pi * math.sqrt(orbit_radius ** 3 / c["moon"]["gm_m3_s2"])
    umbra_s = period_s * math.asin(radius_moon / orbit_radius) / math.pi
    terminal_heat = terminal["source_heat_W"]
    radiators = [{"temperature_K": temp, "effective_fraction": frac,
                  "emitting_area_m2": radiator_area(terminal_heat, temp, p["emissivity"], sigma,
                                                   p["source_capacity_margin"], frac)}
                 for frac in sensitivity["effective_radiator_fraction"]
                 for temp in sensitivity["radiator_temperature_K"]]
    # Recover the inherited storage definition, including discharge efficiency
    # and the protected margin; do not replace it with a new usable fraction.
    operation = data["crater"]["operation"]
    dispatch_fraction = operation["storage_usable_fraction"] * operation["storage_discharge_efficiency"] / (1 + operation["storage_reserve_fraction"])
    city_utilities = data["city"]["utilities"]
    critical_W_person = 1000 * (city_utilities["critical_kW_person"] + city_utilities["food_preservation_critical_kW_person"])
    # This tests energy autonomy, not charging, discharge rating or reliability.
    screens = []
    for row in reference["checkpoints"]:
        reserve_Wh = e["demand"]["reserve_Wh"] * row["population"] / e["normalization"]["population"]
        usable = reserve_Wh * dispatch_fraction
        trips = row["source_departures_year"]
        launches = row["LEO_kg_year"] / t["launch_payload_kg"]
        source_solar_DC = row["source_process_W"] * p["source_capacity_margin"] / p["end_of_life_factor"] / p["source_solar_flux_fraction"]
        beam_DC = row["ground_firm_W"] / (p["RF_efficiency"] * p["beam_capture"] * p["receiver_efficiency"]) * p["stations"] / 2 / p["end_of_life_factor"]
        station_DC_at_end_of_life = beam_DC / p["stations"] * p["end_of_life_factor"]
        station_RF = station_DC_at_end_of_life * p["RF_efficiency"]
        station_conversion_heat = station_DC_at_end_of_life - station_RF
        screens.append({"year": row["year"], "reserve_Wh": reserve_Wh,
                        "dispatchable_energy_Wh_after_margin": usable,
                        "inherited_critical_load_W": row["population"] * critical_W_person,
                        "inherited_critical_load_hours": usable / (row["population"] * critical_W_person),
                        "usable_firm_load_hours": usable / row["ground_firm_W"],
                        "critical_load_hours": {str(f): usable / (row["ground_firm_W"] * f)
                                                for f in sensitivity["critical_load_fractions"]},
                        "beam_installed_DC_W": beam_DC,
                        "beam_hardware_kg": beam_DC / p["space_system_W_kg"],
                        "beam_ideal_solar_area_m2": beam_DC / (c["sun"]["solar_constant_W_m2"] * p["PV_efficiency"]),
                        "beam_station_RF_W_at_end_of_life": station_RF,
                        "beam_station_aperture_average_RF_W_m2": station_RF / (math.pi * (p["aperture_m"] / 2) ** 2),
                        "beam_station_RF_conversion_heat_W": station_conversion_heat,
                        "beam_station_conversion_radiator_m2": radiator_area(station_conversion_heat, p["radiator_K"], p["emissivity"], sigma, p["source_capacity_margin"]),
                        "source_installed_DC_W": source_solar_DC,
                        "source_solar_hardware_kg": source_solar_DC / p["space_system_W_kg"],
                        "source_mine_module_equivalents": row["source_feed_kg_year"] * margin / (z["mine_module_kg_hr"] * hours_year * z["mine_availability"]),
                        "LEO_100t_launches_year": launches,
                        "pad_equivalents": launches / t["launches_per_pad_year"],
                        "source_departures_day": trips / days_year,
                        "source_required_fleet_equivalents": trips * margin * t["nea_cycle_years"] / t["fleet_availability"],
                        "gross_metabolic_O2_kg_year": row["population"] * sensitivity["metabolic_O2_kg_person_day"] * days_year,
                        "one_way_cargo_inventory_kg": {str(years): row["source_cargo_kg_year"] * years
                                                      for years in sensitivity["one_way_transit_years"]}})
    trips = terminal["source_departures_year"]
    max_p = -math.expm1(math.log1p(-sensitivity["illustrative_annual_loss_limit"]) / trips)
    trials = math.ceil(math.log1p(-sensitivity["demonstration_confidence"]) / math.log1p(-max_p))
    traffic = {"trips_year": trips, "loss_cases": [
        {"probability_per_trip": prob, "probability_at_least_one_year": loss_probability(trips, prob)}
        for prob in sensitivity["cargo_loss_probability_per_trip"]],
        "per_trip_limit_for_illustrative_annual_limit": max_p,
        "zero_failure_trials_for_confidence": trials,
        "cargo_kinetic_energy_J": {str(v): 0.5 * t["nea_payload_kg"] * v ** 2 for v in sensitivity["cargo_speed_m_s"]},
        "meaning": "Illustrative independent cargo losses, not Earth-impact risk, qualification or an adopted reliability requirement."}
    pressure = c["physics"]["standard_atmosphere_pa"] * data["crater"]["enclosure"]["reference_pressure_atm"]
    g = c["moon"]["gm_m3_s2"] / radius_moon ** 2
    city = data["city"]["geometry"]
    hall_radius = city["pressure_hall_diameter_m"] / 2
    hall = {"diameter_m": hall_radius * 2, "hoop_force_N_m": pressure * hall_radius,
            "end_force_N": pressure * math.pi * hall_radius ** 2,
            "pressure_only_membrane_thickness_m": pressure * hall_radius / city["effective_allowable_Pa"],
            "detail_mass_factor": city["pressure_detail_factor"]}
    dome = []
    enc = data["crater"]["enclosure"]
    for span in sensitivity["dome_spans_m"]:
        a, h = span / 2, span * enc["cap_height_over_diameter"]
        area = math.pi * (a ** 2 + h ** 2)
        projected = math.pi * a ** 2
        radius = (a ** 2 + h ** 2) / (2 * h)
        shield_pressure = enc["shield_density_kg_m3"] * 3 * g * area / projected
        dome.append({"span_m": span, "cap_radius_m": radius,
                     "uplift_N": pressure * projected,
                     "shield_3m_weight_N": shield_pressure * projected,
                     "shield_3m_uplift_fraction": shield_pressure / pressure,
                     "shield_3m_balance_pressure_Pa": shield_pressure,
                     "pressure_only_thickness_m": pressure * radius / (2 * enc["material_effective_allowable_Pa"])})
    repair = data["reclamation"]["refurbishment"]
    q = repair["collection"] * repair["quality"]
    refurb = [{"restored_life_fraction": alpha, "purchase_ratio": refurbishment_purchase_ratio(q, repair["repair_part_mass_fraction"], alpha)}
              for alpha in sensitivity["refurbished_life_fractions"]]
    primaries = [case for case in cases if case["variant"] == "combined"]
    seed_average_feed = reference["seed"]["source_feed_kg"] / reference["seed"]["bootstrap_years"]
    portfolio_ranges = []
    for portfolio in sorted({case["name"].split(":")[0] for case in primaries}):
        rows = [next(r for r in case["checkpoints"] if r["year"] == 200)
                for case in primaries if case["name"].split(":")[0] == portfolio]
        portfolio_ranges.append({"portfolio": portfolio,
                                 "range_type": "Min/max across scenario combinations, not a probability interval",
                                 **{metric: {"min": min(r[metric] for r in rows),
                                             "max": max(r[metric] for r in rows)}
                                    for metric in ("source_feed_kg_year", "ground_firm_W", "source_process_W", "LEO_kg_year")}})
    return {"schema": "terluna.research.enclosed-feasibility/1",
            "evidence": config["reading_rule"], "boundary": config["boundary"],
            "reference_case": config["reference_case"],
            "requirements": cases, "portfolio_ranges_year200": portfolio_ranges,
            "admission": {variant: {"total": sum(ca["variant"] == variant for ca in cases),
                                     "passed": sum(ca["variant"] == variant and ca["admitted_completed_path"] for ca in cases)}
                          for variant in sorted({ca["variant"] for ca in cases})},
            "beam": {"wavelength_m": wavelength, "ideal_focus_cases": beam,
                     "far_field_boundary_m": 2 * p["aperture_m"] ** 2 / wavelength,
                     "circular_orbit_period_hours": period_s / 3600,
                     "point_sun_central_lunar_shadow_hours": umbra_s / 3600,
                     "qualification": "Ideal uniform focused aperture; no sidelobe, scan, pointing, finite receiver efficiency map or constellation coverage solution."},
            "checkpoint_screens": screens, "source_radiators_year200": radiators,
            "storage_basis": {"usable_fraction": operation["storage_usable_fraction"],
                              "discharge_efficiency": operation["storage_discharge_efficiency"],
                              "reserve_fraction": operation["storage_reserve_fraction"],
                              "dispatch_fraction_after_margin": dispatch_fraction,
                              "critical_W_person": critical_W_person,
                              "target_hours": city_utilities["critical_reserve_hours"]},
            "source_solar_relative_to_flux_assumption": [{"distance_AU": a, "required_area_multiplier": p["source_solar_flux_fraction"] * a ** 2}
                                                        for a in sensitivity["source_distances_AU"]],
            "resource_yield_sensitivity": [{"recovered_yield": y, "feed_for_same_water_kg_year": terminal["source_water_kg_year"] / y}
                                           for y in sensitivity["volatiles_recovered_yield_cases"]],
            "traffic": traffic, "pressure_hall": hall, "crater_domes": dome,
            "refurbishment": {"accepted_fraction": q, "new_parts_fraction": repair["repair_part_mass_fraction"], "cases": refurb},
            "prestart": {"lead_years": config["lead_years"], "programme_start_year": config["programme_start_year"],
                         "opening_year": config["programme_start_year"] + 8,
                         "seed_average_feed_kg_year": seed_average_feed,
                         "seed_average_mine_module_equivalents": seed_average_feed * margin / (z["mine_module_kg_hr"] * hours_year * z["mine_availability"]),
                         "seed_feed_to_opening_annual_feed_ratio": seed_average_feed / reference["checkpoints"][0]["source_feed_kg_year"],
                         "common_programme_USD2026": e["common_programme_USD2041"] / e["normalization"]["USD2026_to_USD2041"],
                         "flat_annual_common_programme_USD2026": e["common_programme_USD2041"] / e["normalization"]["USD2026_to_USD2041"] / config["lead_years"],
                         "budget_status": "Arithmetic distribution of the inherited once-only allowance; not a costed R&D programme or additional charge."}}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE / "results")
    args = parser.parse_args()
    cfg_path = HERE / "feasibility_inputs.json"
    config = json.loads(cfg_path.read_text())
    data = {name: read_pinned(ROOT, entry) for name, entry in config["inputs"].items()}
    result = calculate(config, data)
    result["producer"] = {"path": str(Path(__file__).resolve().relative_to(ROOT)), "sha256": sha256(Path(__file__))}
    result["dependencies"] = config["inputs"] | {"study_inputs": {"path": str(cfg_path.relative_to(ROOT)), "sha256": sha256(cfg_path)}}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "feasibility.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"admission": result["admission"], "reference_case": config["reference_case"],
                      "beam_radius_m": result["beam"]["ideal_focus_cases"][0]["capture_radius_m"],
                      "year200": result["checkpoint_screens"][-2]}, indent=2))


if __name__ == "__main__":
    main()
