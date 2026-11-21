"""Assess dated mining feed against finite hypothetical access allowances.

The resource consumer never constructs a fleet, invents a geological reserve or
rescales feed from population. Physical trajectories arrive as a pinned product.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import gzip
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def nonnegative(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be numeric")
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return float(value)


def fraction(value, name):
    value = nonnegative(value, name)
    if value > 1:
        raise ValueError(f"{name} must be at most one")
    return value


def footprint(feed_kg, parameters):
    feed = nonnegative(feed_kg, "lunar feed")
    density = nonnegative(parameters["bulk_density_kg_m3"], "bulk density")
    depth = nonnegative(parameters["worked_depth_m"], "worked depth")
    access = fraction(parameters["accessible_area_fraction"], "accessible area fraction")
    if min(density, depth, access) <= 0:
        raise ValueError("Density, depth and accessible fraction must be positive")
    excavated = feed/density/depth/1e6
    return {"excavated_area_km2": excavated,
            "prospective_area_km2": excavated/access}


def catalogue_tiers(inventory, scenario):
    groups = inventory["inventory"]["groups"]
    narrow = nonnegative(groups["dark_geometry_excluding_named_science_targets"]
                         ["mass_at_assumed_density_kg"], "narrow catalogue mass")
    broad_group = groups[scenario.get("broader_group", "dark_with_diameter")]
    broad = nonnegative(broad_group["mass_at_assumed_density_kg"], "broad catalogue mass")
    exclusions = []
    if scenario.get("exclude_largest_broad_object", False):
        largest = max(broad_group["largest_objects"], key=lambda r:r["assumed_mass_kg"])
        # Its mass alone exceeds the entire narrow set: it cannot be a member.
        if largest["assumed_mass_kg"] <= narrow:
            raise ValueError("Excluded object's absence from the narrow set is unestablished")
        broad -= largest["assumed_mass_kg"]
        exclusions.append(deepcopy(largest))
    if narrow > broad:
        raise ValueError("Nested catalogue screen exceeds its parent")
    f = fraction(scenario["narrow_extractable_fraction"], "narrow extractable fraction")
    g = fraction(scenario["broader_remainder_extractable_fraction"], "broader extractable fraction")
    return [
        {"id": "narrow_geometry", "catalogue_mass_kg": narrow,
         "extractable_fraction": f, "feed_allowance_kg": narrow*f,
         "available_year": 0, "annual_feed_limit_kg": None},
        {"id": "broader_remainder", "catalogue_mass_kg": broad-narrow,
         "extractable_fraction": g, "feed_allowance_kg": (broad-narrow)*g,
         "catalogue_group": scenario.get("broader_group", "dark_with_diameter"),
         "excluded_objects": exclusions,
         "available_year": 0, "annual_feed_limit_kg": None},
    ]


def allocate_feed(annual, tiers):
    """Allocate raw feed chronologically; deficits do not become hidden backlog."""
    remaining = {}
    for tier in tiers:
        name = tier["id"]
        if name in remaining:
            raise ValueError("Repeated access-tier identity")
        remaining[name] = nonnegative(tier["feed_allowance_kg"], "feed allowance")
        available = tier.get("available_year", 0)
        if isinstance(available, bool) or not isinstance(available, int):
            raise ValueError("Tier availability year must be an integer")
        if tier.get("annual_feed_limit_kg") is not None:
            nonnegative(tier["annual_feed_limit_kg"], "annual feed limit")
    rows = []
    requested_total = supplied_total = 0.
    previous = None
    for record in annual:
        year = record["year"]
        if isinstance(year, bool) or not isinstance(year, int):
            raise ValueError("Physical years must be integers")
        if previous is not None and year <= previous:
            raise ValueError("Physical years must be unique and increasing")
        previous = year
        request = nonnegative(record["source_feed_kg"], "annual source feed")
        shortage = request
        draws = {}
        for tier in tiers:
            name = tier["id"]
            limit = tier.get("annual_feed_limit_kg")
            available = remaining[name] if year >= tier.get("available_year", 0) else 0.
            if limit is not None:
                available = min(available, limit)
            draw = min(shortage, available)
            draws[name] = draw
            remaining[name] -= draw
            shortage -= draw
        supplied = request-shortage
        requested_total += request
        supplied_total += supplied
        rows.append({"year": year, "requested_feed_kg": request,
                     "supplied_under_scenario_kg": supplied,
                     "unallocated_feed_kg": shortage, "tier_draw_kg": draws,
                     "remaining_by_tier_kg": dict(remaining),
                     "cumulative_requested_feed_kg": requested_total,
                     "cumulative_supplied_feed_kg": supplied_total,
                     "cumulative_unallocated_feed_kg": requested_total-supplied_total})
    return {"tiers": deepcopy(tiers), "annual": rows,
            "requested_feed_kg": requested_total,
            "supplied_under_scenario_kg": supplied_total,
            "unallocated_feed_kg": requested_total-supplied_total,
            "first_unallocated_year": next((r["year"] for r in rows if r["unallocated_feed_kg"] > 0), None),
            "within_hypothetical_mass_allowance": all(r["unallocated_feed_kg"] == 0 for r in rows),
            "qualified_reserve_sufficiency": None}


def assess_case(annual, inventory, config):
    """Consume normalized mass amounts; each record covers its declared year."""
    lunar = source = 0.
    rows = []
    for record in annual:
        lunar_feed = nonnegative(record["lunar_virgin_feed_kg"], "annual lunar feed")
        source_feed = nonnegative(record["source_feed_kg"], "annual source feed")
        lunar += lunar_feed
        source += source_feed
        rows.append({"year": record["year"],
                     "period_kind": record.get("period_kind", "operation"),
                     "lunar_virgin_feed_kg": lunar_feed,
                     "source_feed_kg": source_feed,
                     "cumulative_lunar_virgin_feed_kg": lunar,
                     "cumulative_source_feed_kg": source,
                     **footprint(lunar, config["lunar_footprint"])})
    access = {name: allocate_feed(rows, catalogue_tiers(inventory, assumptions))
              for name, assumptions in config["access_scenarios"].items()}
    checkpoints = [r for r in rows if r["year"] in config["checkpoints"]]
    access_checkpoints = []
    for record in checkpoints:
        comparisons = {}
        for name, result in access.items():
            state = next(r for r in result["annual"] if r["year"] == record["year"])
            comparisons[name] = {
                "cumulative_requested_feed_kg": state["cumulative_requested_feed_kg"],
                "cumulative_unallocated_feed_kg": state["cumulative_unallocated_feed_kg"],
                "within_hypothetical_mass_allowance_to_checkpoint": state["cumulative_unallocated_feed_kg"] == 0.,
                "qualified_resource_and_delivery_schedule": None}
        access_checkpoints.append({"year": record["year"], "comparisons": comparisons})
    return {"annual": rows, "source_access": access,
            "checkpoints": checkpoints, "source_access_checkpoints": access_checkpoints,
            "resource_condition": "Financial and physical completion remain conditional on qualifying composition, extraction, permission and dated delivery. Stock allowance comparisons do not admit a mine or validate the inherited chemical route.",
            "cumulative_lunar_virgin_feed_kg": lunar,
            "cumulative_source_feed_kg": source}


def read_field(record, path):
    """Read the admitted schema's field path without importing its producer."""
    for key in ([path] if isinstance(path, str) else path):
        record = record[key]
    return record


def normalize_case(case, fields):
    """Include seed feed once, followed by the producer's actual annual feed.

    The producer supplies an average bootstrap screen and aggregate seed feed.
    Even allocation before the first operating year preserves those totals;
    it does not establish a preproduction deployment or delivery schedule.
    """
    if not case["annual"]:
        raise ValueError("At least one physical operating year is required")
    first = read_field(case["annual"][0], fields["year"])
    if isinstance(first, bool) or not isinstance(first, int):
        raise ValueError("Physical years must be integers")
    seed = case["seed"]
    duration = seed["bootstrap_years"]
    if isinstance(duration, bool) or not isinstance(duration, int) or duration <= 0:
        raise ValueError("Bootstrap duration must be a positive integer")
    seed_lunar = nonnegative(seed["lunar_raw_kg"], "seed lunar feed")
    seed_source = nonnegative(seed["source_feed_kg"], "seed source feed")
    annual = [{"year": y, "period_kind": "bootstrap_average",
               "lunar_virgin_feed_kg": seed_lunar/duration,
               "source_feed_kg": seed_source/duration}
              for y in range(first-duration, first)]
    annual.extend({"year": read_field(r, fields["year"]), "period_kind": "operation",
                   "lunar_virgin_feed_kg": read_field(r, fields["lunar_virgin_feed_kg"]),
                   "source_feed_kg": read_field(r, fields["source_feed_kg"])}
                  for r in case["annual"])
    return annual


def load_checked(spec, base=HERE):
    path = base/spec["path"]
    if spec.get("manifest_path"):
        manifest_path = base/spec["manifest_path"]
        if digest(manifest_path) != spec["manifest_sha256"]:
            raise ValueError("Unreviewed source manifest: "+str(manifest_path))
        manifest = json.loads(manifest_path.read_text())
        if (manifest["sha256"] != spec["sha256"] or manifest["decompressed_sha256"] != spec["content_sha256"]
                or manifest["product_schema"] != spec["schema"] or manifest["path"] != path.name):
            raise ValueError("Source manifest and input identity disagree")
    if digest(path) != spec["sha256"]:
        raise ValueError("Unreviewed source product: "+str(path))
    content = path.read_bytes()
    if spec.get("compression") == "gzip":
        content = gzip.decompress(content)
        if hashlib.sha256(content).hexdigest() != spec["content_sha256"]:
            raise ValueError("Unreviewed decompressed source content: "+str(path))
    elif spec.get("compression") is not None:
        raise ValueError("Unsupported source compression: "+str(spec["compression"]))
    product = json.loads(content)
    if product["schema"] != spec["schema"]:
        raise ValueError("Unexpected source schema: "+str(path))
    return product


def inventory_inputs(config):
    inventory = load_checked(config["inventory_product"])
    geometry = load_checked(config["geometry_product"])
    for name, group in geometry["groups"].items():
        if name in inventory["inventory"]["groups"]:
            raise ValueError("Repeated catalogue group")
        inventory["inventory"]["groups"][name] = group
    return inventory


def read_inputs():
    config = json.loads((HERE/"inputs.json").read_text())
    inventory = inventory_inputs(config)
    if not config["physical_products"]:
        raise ValueError("The physical producer product has not yet been admitted; no feed trajectory is substituted")
    return config, inventory


def assess_physical_cases(physical, inventory, config):
    cases = {}
    for case in physical["cases"]:
        key = case["physical_scenario"]+":"+case["name"]
        if key in cases:
            raise ValueError("Repeated physical case/scenario identity: "+key)
        fields = config["annual_fields"]
        annual = normalize_case(case, fields)
        result = assess_case(annual, inventory, config)
        result["physical_case_name"] = case["name"]
        result["physical_scenario"] = case["physical_scenario"]
        result["physical_convergence"] = deepcopy(case["convergence"])
        result["first_physical_capacity_shortfall"] = case["first_capacity_shortfall"]
        result["seed_evidence"] = case["seed"]["evidence"]
        result["supply_driver_checkpoints"] = []
        for row in case["annual"]:
            if row["year"] not in config["checkpoints"] or "physical" not in row:
                continue
            p = row["physical"]
            result["supply_driver_checkpoints"].append({
                "year": row["year"], "source_feed_kg_year": p["source_feed_kg_year"],
                "source_water_kg_year": p["source_water_kg_year"],
                "source_nitrogen_kg_year": p["source_nitrogen_kg_year"],
                "source_virgin_metal_kg_year": p["source_virgin_metal_kg_year"],
                "depot_water_kg_year": p["depot_water_kg_year"],
                "source_ferry_fuel_kg_year": p["source"]["burned"],
                "lunar_route_propellant_supply_kg_year": p["Moon"]["space_supply"],
                "Earth_route_propellant_supply_kg_year": p["Earth"]["space_supply"],
                "resident_guest_and_initial_water_kg_year": p["source_water_kg_year"]-p["depot_water_kg_year"]-p["source"]["burned"],
                "source_cargo_kg_year": p["source_cargo_kg_year"],
                "source_tug_count": p["source_tug_count"],
                "source_tug_component_import_kg_year": p["source_tug_component_import_kg_year"],
                "source_feed_process_W": p["source_feed_process_W"],
                "source_process_W": p["source_process_W"], "source_heat_W": p["source_heat_W"]})
        cases[key] = result
    return cases


def product():
    config, inventory = read_inputs()
    cases = {}
    identities = []
    dependencies = {"inputs.json": digest(HERE/"inputs.json"),
                    "sources.json": digest(HERE/"sources.json"),
                    config["inventory_product"]["path"]: config["inventory_product"]["sha256"],
                    config["geometry_product"]["path"]: config["geometry_product"]["sha256"]}
    for spec in config["physical_products"]:
        physical = load_checked(spec)
        new = assess_physical_cases(physical, inventory, config)
        if set(cases) & set(new):
            raise ValueError("Repeated physical cases across admitted products")
        cases.update(new)
        dependencies[spec["path"]] = spec["sha256"]
        if spec.get("manifest_path"):
            dependencies[spec["manifest_path"]] = spec["manifest_sha256"]
        identities.append({**deepcopy(spec), "producer": physical["producer"],
                           "evidence": physical["evidence"], "scenarios": physical["scenarios"]})
    return {"schema": "terluna.resources.enclosed-mining/1",
            "producer": {"domain": "resources", "path": "resources/enclosed_mining/run.py",
                         "sha256": digest(__file__)},
            "evidence": config["evidence"],
            "units": {"feed": "kg", "area": "km2", "time": "elapsed programme year"},
            "dependencies": dependencies,
            "physical_products": identities,
            "seed_timing_rule": "Aggregate seed feed is allocated evenly over bootstrap_years immediately before the first operating year. This preserves the producer average screen; it is not a qualified deployment schedule.",
            "reading_rule": "A mass-allowance pass is conditional on the selected fractions and instantaneous access assumption. It proves no geological, operational, permitted or economic reserve. Cost and production responses remain with the physical producer.",
            "rules": config["rules"], "access_rule": config["access_rule"],
            "footprint_assumptions": config["lunar_footprint"], "cases": cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE/"results/assessment.json.gz")
    args = parser.parse_args()
    result = product()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(result, separators=(",", ":"), allow_nan=False)+"\n").encode()
    args.output.write_bytes(gzip.compress(raw, mtime=0) if args.output.suffix == ".gz" else raw)
    print(args.output)


if __name__ == "__main__":
    main()
