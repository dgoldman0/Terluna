#!/usr/bin/env python3
"""Reproduce the recovered Open Moon industrial scale calculations.

These are transparent scale calculations and do not replace the Terluna domain models.
All SI unless explicitly labelled.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
YEAR = 365.25 * 86400.0
SIGMA = 5.670374419e-8
EPS = 0.9
EV = 1.602176634e-19
AMU = 1.66053906660e-27
SOLAR_LUMINOSITY = 3.828e26
TNT_MT_J = 4.184e15


def kardashev_power(k: float) -> float:
    return 10.0 ** (10.0 * k + 6.0)


def dt_specific_energy() -> float:
    # 17.6 MeV per reaction; reacting mass is D + T.
    e = 17.6e6 * EV
    m = (2.01410177812 + 3.01604928199) * AMU
    return e / m


def fission_specific_energy() -> float:
    # Simple scale: 200 MeV per fission of 235-u nucleus.
    return 200e6 * EV / (235.0 * AMU)


def radiator_area(power_w: float, temp_k: float, emissivity: float = EPS) -> float:
    return power_w / (emissivity * SIGMA * temp_k**4)


def solar_area(power_w: float, net_w_m2: float = 300.0) -> float:
    return power_w / net_w_m2


def kinetic_energy(mass_kg: float, speed_m_s: float) -> float:
    return 0.5 * mass_kg * speed_m_s**2


def main():
    k = 1.17
    P = kardashev_power(k)
    dt_e = dt_specific_energy()
    fission_e = fission_specific_energy()

    summary = {
        "sagan_k": k,
        "civilization_power_W": P,
        "civilization_power_PW": P / 1e15,
        "share_of_solar_luminosity": P / SOLAR_LUMINOSITY,
        "dt_specific_energy_J_kg": dt_e,
        "dt_fuel_kg_per_year_ideal": P * YEAR / dt_e,
        "dt_fuel_kg_per_year_at_50_percent_net": P * YEAR / (0.5 * dt_e),
        "fission_specific_energy_J_kg_200MeV_235u": fission_e,
        "fissioned_mass_kg_per_year": P * YEAR / fission_e,
        "solar_area_100PW_1AU_km2_at_300Wm2": solar_area(100e15) / 1e6,
        "solar_area_500PW_1AU_km2_at_300Wm2": solar_area(500e15) / 1e6,
        "solar_area_500PW_0p1AU_ideal_km2": solar_area(500e15, 300.0 * 100.0) / 1e6,
    }
    summary["radiator_area_km2_for_K1p17"] = {
        str(T): radiator_area(P, T) / 1e6 for T in (350, 600, 1000, 1500)
    }

    # Historical packet arithmetic used 2.7e8 kg/s as a representative bulk flow.
    flow = 2.7e8
    packet_rows = []
    for mass in (1e6, 1e8, 1e9, 1e12, 1e13, 1e14, 1e15):
        interval = mass / flow
        packet_rows.append({
            "packet_mass_kg": mass,
            "packets_per_s": flow / mass,
            "interval_s": interval,
            "interval_minutes": interval / 60.0,
            "interval_hours": interval / 3600.0,
            "interval_days": interval / 86400.0,
        })
    summary["bulk_flow_reference_kg_s"] = flow
    summary["six_year_pipeline_mass_kg"] = flow * 6.0 * YEAR

    kinetic_rows = []
    for mass in (1e6, 1e8, 1e9, 1e12):
        E = kinetic_energy(mass, 10_000.0)
        kinetic_rows.append({
            "packet_mass_kg": mass,
            "speed_km_s": 10.0,
            "kinetic_energy_J": E,
            "TNT_equivalent_megatons": E / TNT_MT_J,
        })

    # Baseline local-arrival heat figures use a ~300-year N2 rate around 2.67e8 kg/s.
    n2_flow = 267_095_731.49278048
    moon_area = 4 * math.pi * (1_737_400.0**2)
    arrival_rows = []
    for speed in (3_000.0, 10_000.0):
        power = kinetic_energy(n2_flow, speed)  # because 0.5 * mdot * v^2 is W
        arrival_rows.append({
            "mass_flow_kg_s": n2_flow,
            "speed_km_s": speed / 1000.0,
            "dissipation_power_W": power,
            "dissipation_power_PW": power / 1e15,
            "global_mean_W_m2": power / moon_area,
        })

    (OUT / "derived_industry_calculations.json").write_text(json.dumps({
        "summary": summary,
        "packet_rates": packet_rows,
        "packet_kinetic_energy": kinetic_rows,
        "arrival_heat": arrival_rows,
        "notes": [
            "K=1.17 is a civilization-scale scenario descriptor, not an asserted Moon-project load.",
            "0.1 AU solar-area scaling assumes inverse-square sunlight and unchanged conversion fraction; thermal/material closure is omitted.",
            "Radiator areas are ideal one-sided deep-space eps=0.9 areas.",
            "Large packet masses are throughput illustrations; recovered design discussion favored smaller 1e6-1e8 kg packets for bounded consequences.",
        ],
    }, indent=2) + "\n")

    for name, rows in [
        ("packet_rates.csv", packet_rows),
        ("packet_kinetic_energy.csv", kinetic_rows),
        ("arrival_heat.csv", arrival_rows),
    ]:
        with (OUT / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)

    radiator_rows = [{"temperature_K": T, "area_km2": radiator_area(P, T) / 1e6}
                     for T in (350, 600, 1000, 1500)]
    with (OUT / "radiator_area.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(radiator_rows[0])); w.writeheader(); w.writerows(radiator_rows)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
