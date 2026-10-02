"""Run the initial wind-wave comparison and numerical checks.

python -m climate.waves.run --build-root /path/to/swan-build

Raw runs remain under ignored research/runs/waves. This exports the selected
numbers and their evidence/provenance to climate/waves/results.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np

from climate.waves.model import (Case, ROOT, RUNS, STANDARD_GRAVITY,
                                 main_cases, read_profile, run_case, sha256, wind_inputs)

HERE = Path(__file__).resolve().parent
EVIDENCE = (
    "SWAN numerical experiments with uniform persistent winds and flat 1000 m water, "
    "initially calm, no incoming swell or currents. Wind magnitudes come from corrected "
    "CM1 water-column statistics; direction, fetch and persistence are imposed. "
    "The air density is held at SWAN's 1.28 kg/m3 and water at 1025 kg/m3. "
    "These are conditional responses, not wave climatology or empirically validated lunar predictions."
)


def sampled_rows(case: Case, cube: np.ndarray, variant: str) -> list[dict]:
    rows = []
    for hours in (1, 3, 6, 12, 24, 48):
        for fetch_km in (10, 25, 50, 100):
            ti = int(round(hours * 3600 / (case.step_s * case.output_every_steps)))
            xi = int(round(fetch_km * 1000 * case.cells / case.length_m))
            row = cube[ti, xi]
            if not np.isclose(row[0], hours * 3600) or not np.isclose(row[1], fetch_km * 1000):
                raise RuntimeError("Requested summary point is not on the output grid")
            rows.append(dict(build=variant, case=case.name, wind_m_s=case.wind_m_s,
                             gravity_m_s2=case.gravity_m_s2, duration_h=hours,
                             fetch_km=fetch_km, hs_m=float(row[2]),
                             peak_period_s=float(row[3]), mean_period_s=float(row[4])))
    return rows


def relative_comparison(reference: np.ndarray, candidate: np.ndarray) -> dict:
    """Compare corresponding samples; avoid relative errors in almost calm cells."""
    mask = reference[..., 2] >= 0.01
    if not mask.any():
        raise ValueError("No finite-amplitude samples for comparison")
    result = {}
    for index, name in ((2, "hs"), (3, "peak_period"), (4, "mean_period")):
        difference = np.abs(candidate[..., index][mask] / reference[..., index][mask] - 1)
        result[name + "_relative_max"] = float(difference.max())
        result[name + "_relative_rms"] = float(np.sqrt(np.mean(difference**2)))
        result[name + "_final_endpoint_relative"] = float(
            candidate[-1, -1, index] / reference[-1, -1, index] - 1)
    result["samples"] = int(mask.sum())
    return result


def experiment(build_root: Path, run_root: Path = RUNS) -> dict:
    build = json.loads((build_root / "build.json").read_text())
    executables = {"stock": build_root / "stock/swan.exe",
                   "patched": build_root / "patched_agrow/swan.exe"}
    for name, executable in executables.items():
        key = "stock_executable" if name == "stock" else "patched_executable"
        if sha256(executable) != build[key]["sha256"]:
            raise RuntimeError(f"Executable does not match build manifest: {name}")
    records, cubes = {}, {}

    def execute(case, variant):
        folder = run_root / variant / case.name
        record = run_case(case, executables[variant], folder)
        key = variant + "/" + case.name
        records[key] = record
        cubes[key] = read_profile(folder / "profile.tbl", case)
        print(f"{key}: {record['elapsed_wall_s']:.1f} s; "
              f"final Hs {cubes[key][-1, -1, 2]:.4g} m", flush=True)
        return cubes[key]

    # Exact sixfold numerical scaling avoids roundoff of the integer SWAN clock.
    # gEarth/6 is a test gravity, distinct from the actual GM/R² used below.
    winds, forcing = wind_inputs()
    similarity_high = Case("similarity_earth", STANDARD_GRAVITY, winds["median"],
                           length_m=20_000, depth_m=200, cells=80,
                           step_s=60, steps=360, output_every_steps=60)
    similarity_low = replace(similarity_high, name="similarity_sixth_g",
                             gravity_m_s2=STANDARD_GRAVITY / 6,
                             length_m=120_000, depth_m=1200, step_s=360)
    similarity = {}
    repeatability = {}
    for variant in executables:
        high = execute(similarity_high, variant)
        repeated = execute(replace(similarity_high, name="similarity_earth_repeat"), variant)
        repeatability[variant] = bool(np.array_equal(high, repeated))
        if not repeatability[variant]:
            raise RuntimeError(f"Identical SWAN inputs were not repeatable: {variant}")
        low = execute(similarity_low, variant).copy()
        low[..., :] /= 6  # time, distance, height, both periods and depth all scale alike.
        similarity[variant] = relative_comparison(high, low)
    similarity["definition"] = dict(length_depth_time_height_period_ratio=6,
                                     low_gravity_m_s2=STANDARD_GRAVITY / 6,
                                     note="Exact numerical scaling control; actual Moon runs use GM/R²")
    similarity["patched_pass"] = bool(similarity["patched"]["hs_relative_max"] < 0.01
                                       and similarity["patched"]["mean_period_relative_max"] < 0.01)
    if not similarity["patched_pass"]:
        raise RuntimeError("Patched gravity similarity failed; inspect raw runs before proceeding")

    rows = []
    cases = main_cases()
    for variant in executables:
        for case in cases:
            rows.extend(sampled_rows(case, execute(case, variant), variant))

    earth_preserved = all(np.array_equal(cubes["stock/" + case.name], cubes["patched/" + case.name])
                          for case in cases if case.name.startswith("earth_"))
    if not earth_preserved:
        raise RuntimeError("AGROW patch changed Earth-reference output; investigate before export")

    baseline = next(case for case in cases if case.name == "moon_p90")
    base = cubes["patched/" + baseline.name]
    refinements = (
        replace(baseline, name="moon_p90_half_dt", step_s=150, steps=1152, output_every_steps=24),
        replace(baseline, name="moon_p90_half_dx", cells=200),
        replace(baseline, name="moon_p90_more_directions", direction_bins=72),
        replace(baseline, name="moon_p90_wider_band", earth_low_hz=0.015,
                earth_high_hz=6.0, frequency_intervals=62),
    )
    checks = {}
    for case in refinements:
        values = execute(case, "patched")
        if case.cells == 2 * baseline.cells:
            values = values[:, ::2]
        checks[case.name] = relative_comparison(base[1:, 10:], values[1:, 10:])
    # A practical first-screen criterion on final Hs/Tm at 100 km, with complete
    # profiles retained so larger early-time/local differences remain visible.
    endpoint_pass = all(abs(item["hs_final_endpoint_relative"]) < 0.05
                        and abs(item["mean_period_final_endpoint_relative"]) < 0.05
                        for item in checks.values())
    edges_pass = all(not record["spectral_edges"]["peak_at_band_edge"]
                     and record["spectral_edges"]["low_two_interval_fraction"] < 0.001
                     and record["spectral_edges"]["high_two_interval_fraction"] < 0.01
                     for record in records.values())
    breaking_inactive = all(record["maximum_breaking_fraction"] == 0 for record in records.values())
    producer_files = [HERE / name for name in ("model.py", "run.py", "build.py")]
    return dict(schema="terluna.climate.wind-waves/1", evidence=EVIDENCE,
                reading_rule="Compare matching wind, fetch, duration and build. Never read the wind labels as wave percentiles. "
                             "Only the tested p90 lunar endpoint has resolution-sensitivity checks; consult their full errors. "
                             "Early growth is sensitive to timestep and AGROW choices. These Earth controls use the same "
                             "lunar-derived winds and do not represent Earth's actual wave climate.",
                units=dict(hs_m="m", peak_period_s="s", mean_period_s="s", wind_m_s="m/s",
                           gravity_m_s2="m/s²", fetch_km="km", duration_h="hours"),
                producer=dict(domain="climate", model="climate/waves/model.py",
                              files={str(path.relative_to(ROOT)): sha256(path) for path in producer_files}),
                inputs=dict(wind=forcing, shared_constants_sha256=sha256(ROOT / "shared/constants.json")),
                build=build, cases=[asdict(case) for case in cases],
                checks=dict(gravity_similarity=similarity, earth_output_preserved=earth_preserved,
                            identical_input_repeatability=repeatability,
                            depth_breaking_inactive=breaking_inactive,
                            refinements=checks, endpoint_resolution_pass=endpoint_pass,
                            spectral_edges_pass=edges_pass,
                            criteria=dict(similarity_hs_tm_relative_max=0.01,
                                          comparison_minimum_reference_hs_m=0.01,
                                          endpoint_hs_tm_absolute_relative_change=0.05,
                                          final_low_edge_fraction_max=0.001,
                                          final_high_edge_fraction_max=0.01,
                                          frequency_bin_convergence_tested=False),
                            empirical_validation="Not established"),
                rows=rows, runs=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, default=RUNS)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    result = experiment(args.build_root.resolve(), args.run_root.resolve())
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "waves.json").write_text(json.dumps(result, indent=2) + "\n")
    with (args.output / "waves.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result["rows"][0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(result["rows"])
    print(json.dumps(result["checks"], indent=2))
    if not (result["checks"]["endpoint_resolution_pass"] and result["checks"]["spectral_edges_pass"]
            and result["checks"]["depth_breaking_inactive"]):
        raise SystemExit("A numerical check needs review; exported results retain the failed check.")


if __name__ == "__main__":
    main()
