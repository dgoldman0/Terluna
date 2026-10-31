"""Schema- and byte-bound inventory; no cross-domain implementation imports."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = "terluna.research.world-economy-inventory/1"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(root=ROOT, manifest_path=HERE / "inputs.json"):
    manifest = json.loads(Path(manifest_path).read_text())
    if manifest.get("schema") != "terluna.research.world-economy-inputs/1":
        raise ValueError("Unexpected inventory manifest schema")
    selected = {}
    for name, spec in manifest["products"].items():
        path = Path(root) / name
        if digest(path) != spec["sha256"]:
            raise ValueError(f"Input changed; review and rebind explicitly: {name}")
        data = json.loads(path.read_text())
        if data.get("schema") != spec["schema"]:
            raise ValueError(f"Unexpected input schema: {name}")
        values = {}
        for field in spec["fields"]:
            value = data
            for key in field.split("."):
                value = value[key]
            values[field] = value
        selected[name] = dict(schema=data["schema"], sha256=spec["sha256"],
                              evidence=data["evidence"], values=values)
    return manifest, selected


def results():
    manifest, sources = load_inputs()
    def value(path, field):
        return sources[path]["values"][field]
    population = value("research/studies/distributions/results/distributions.json",
                       "distributions.C'.people_billion")
    moon, earth, array = [population[k] for k in ("moon", "earth", "array")]
    shares = value("research/studies/inhabited_volume/results/inhabited_volume.json",
                   "assumptions.lift_allocation_shares")
    districts = value("research/studies/inhabited_volume/results/inhabited_volume.json",
                      "assumptions.districts")
    district_gt = population["moon_aloft"] * districts["central50"]["total_t_per_person"]
    gross_lift_gt = district_gt / shares["district_payload"]
    h2_gt = population["moon_aloft"] * value(
        "research/studies/ways_of_living/results/ways_of_living.json",
        "aerial_reference.lifting_species_t_per_resident")
    array_gt = array * value("research/studies/distributions/results/distributions.json",
                            "transcribed.array_habitat_gt_per_billion.value")
    stock = value("resources/results/industry.json", "people.stock_t_per_person")
    return dict(schema=SCHEMA,
        producer=dict(study="world_economy", files={
            str(p.relative_to(ROOT)): digest(p) for p in (HERE / "inventory.py", HERE / "inputs.json")}),
        evidence="Reproducible inventory of committed engineering screens and accepted population; no construction quotation, market estimate, validated capacity or GDP forecast.",
        reading_rule="Copied quantities retain source units and assumptions. Gt=1e12 kg. Population is billion residents. Stocks, annual throughputs, losses and alternative estimates are separate; their rows are not additive. See README for accounting boundaries.",
        baseline_commit=manifest["baseline_commit"], population_billion=population,
        equal_per_resident_benchmark=dict(moon_over_earth=moon/earth, array_over_earth=array/earth,
            moon_and_array_over_earth=(moon+array)/earth,
            moon_per_resident_premium_for_equal_total=earth/moon-1),
        derived=dict(moon_ordinary_material_stock_gt=[moon*x for x in stock],
            aerial_district_gt=district_gt,
            aerial_envelope_flight_allowance_gt=gross_lift_gt*shares["envelope_flight_structure"],
            aerial_control_allowance_gt=gross_lift_gt*shares["control_flight_services"],
            aerial_unmanufactured_reserve_gt=gross_lift_gt*shares["reserve"],
            aerial_hydrogen_gt=h2_gt, array_habitat_gt=array_gt),
        sources=sources,
        exclusions=["Array habitat construction streams are growth, not mature renewal.",
            "Material-stock proxies overlap explicit buildings and infrastructure.",
            "Total delivered electricity already includes food and close computing.",
            "Lift reserve is capacity, not manufactured mass.",
            "Array habitat shielding geometries are alternatives, not additions.",
            "Braking energy is managed energy, not automatically purchased electricity."])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "research/runs/world_economy/inventory.json")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results(), indent=2)+"\n")
    print(args.output)

if __name__ == "__main__":
    main()
