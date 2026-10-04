"""Algebraic authority and inventory scales for a proposed cycling shield fleet.

This does not generate or validate an orbit, handover, or usable-area fraction.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from shared import constants as K
from shared.provenance import constants_changed, constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = HERE / "results/holding.json"
OUTPUT = HERE / "results/tacking_screen.json"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def screen():
    source = json.loads(SOURCE.read_text())
    producer = source["producer"]
    if constants_changed(producer["constants"]):
        raise ValueError("The holding product's shared constants have changed")
    for name, expected in producer["source_hashes"].items():
        if digest(ROOT / name) != expected:
            raise ValueError(f"Holding source has changed: {name}")
    if digest(ROOT / producer["snapshot_exporter"]) != producer["exporter_sha256"]:
        raise ValueError("The holding snapshot exporter has changed")
    case = source["aperture_quadrature"]["variable_distance_finer"]
    area, mass = case["aperture_area_m2"], case["base_optical_mass_kg"]
    sigma = mass / area
    fraction = source["core_diverted_fraction"]
    annual_fuel = case["propellant_kg_s"] * K.JULIAN_DAY * K.JULIAN_YEAR_DAYS
    photon_acceleration = K.SOLAR_CONSTANT / (K.SPEED_OF_LIGHT * sigma)
    inventory = []
    for usable in (0.5, 0.1, 0.01):
        inventory.append(dict(
            assumed_usable_projected_area_fraction=usable,
            deployed_area_multiplier=1 / usable,
            base_optical_inventory_kg=mass / usable,
            inventory_in_reference_propellant_years=mass / usable / annual_fuel,
        ))
    return dict(
        schema="terluna.research.solar-shield-tacking-screen/1",
        producer=dict(lane="research", runner=str(Path(__file__).relative_to(ROOT)),
                      source_hashes={str(Path(__file__).relative_to(ROOT)): digest(__file__)},
                      constants=constants_used([Path(__file__)]),
                      inputs={str(SOURCE.relative_to(ROOT)): digest(SOURCE)}),
        evidence="Algebraic scales from the existing 50 g/m2 holding product. "
                 "No cycling orbit, usable fraction or coverage history has been solved.",
        reading_rule="Acceleration bounds assume an unloaded Sun-facing patch at 1 AU. "
                     "Usable fractions are illustrative assumptions. Inventory numbers "
                     "are necessary area scales under perfect placement, not sufficient "
                     "coverage or mission cost estimates; power hardware, control, "
                     "habitats, replacement, reserves and deployment are excluded. "
                     "Reference propellant years compare masses only.",
        reference=dict(case="aperture_quadrature/variable_distance_finer",
                       base_areal_mass_kg_m2=sigma, aperture_area_m2=area,
                       base_optical_mass_kg=mass,
                       propellant_kg_s=case["propellant_kg_s"],
                       expended_propellant_per_julian_year_kg=annual_fuel,
                       core_diverted_fraction=fraction),
        optical_bounds=dict(
            full_spectrum_max_anti_sunward_mm_s2=2 * photon_acceleration * 1000,
            full_spectrum_max_transverse_relaxed_mm_s2=photon_acceleration * 1000,
            core_max_anti_sunward_mm_s2=2 * fraction * photon_acceleration * 1000,
            core_max_transverse_relaxed_mm_s2=fraction * photon_acceleration * 1000,
            max_sunward_mm_s2=0,
            hardware_and_payload_included=False,
            physical_filter_force_validated=False,
        ),
        illustrative_inventory=inventory,
        units=dict(area="m2", mass="kg", acceleration="mm/s2",
                   time="Julian year", fractions="dimensionless"),
    )


if __name__ == "__main__":
    product = screen()
    OUTPUT.write_text(json.dumps(product, indent=2) + "\n")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}; no orbit feasibility is implied")
