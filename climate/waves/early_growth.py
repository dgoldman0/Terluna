"""Separate time-step and frequency-spacing checks for the first six hours.

These runs hold air density at the original experiment's fixed value to
isolate numerical refinement.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from climate.waves.model import Case, ROOT, read_profile, run_case, sha256, wind_inputs
from shared.constants import MOON_SURFACE_GRAVITY


def compare(coarse, fine):
    """Report matched points from 20 minutes onward and at least 1 cm high."""
    mask = (fine[..., 0] >= 1200) & (fine[..., 1] >= 10_000) & (fine[..., 2] >= .01)
    result = {}
    for col, name in ((2, "hs"), (4, "mean_period")):
        error = np.abs(coarse[..., col][mask] / fine[..., col][mask] - 1)
        result[name + "_relative_max"] = float(error.max())
        result[name + "_relative_rms"] = float(np.sqrt(np.mean(error**2)))
        result[name + "_absolute_max"] = float(np.max(np.abs(coarse[..., col][mask] - fine[..., col][mask])))
    result["samples"] = int(mask.sum())
    result["within_five_percent"] = bool(max(result["hs_relative_max"], result["mean_period_relative_max"]) < .05)
    return result


def experiment(build_root, run_root):
    executable = build_root / "patched_agrow/swan.exe"
    build = json.loads((build_root / "build.json").read_text())
    if sha256(executable) != build["patched_executable"]["sha256"]:
        raise ValueError("SWAN executable does not match the build manifest")
    winds, forcing = wind_inputs()
    cubes, records, cases = {}, {}, []
    settings = ([(dt, 48, 21600) for dt in (300, 150, 75, 30, 15)]
                + [(15, n, 21600) for n in (96, 192)] + [(dt, 48, 7200) for dt in (5, 2, 1)])
    for dt, nf, duration in settings:
        case = Case(f"early_dt{dt}_nf{nf}", MOON_SURFACE_GRAVITY, winds["p90"],
                    step_s=dt, steps=duration // dt, output_every_steps=300 // dt,
                    frequency_intervals=nf)
        folder = run_root / case.name
        records[case.name] = run_case(case, executable, folder, timeout_s=600)
        cubes[(dt, nf)] = read_profile(folder / "profile.tbl", case)
        cases.append(asdict(case))
        print(case.name, f"Hs(1h,{duration//3600}h)", cubes[(dt, nf)][12, -1, 2], cubes[(dt, nf)][-1, -1, 2], flush=True)
    checks = {f"dt_{a}_to_{b}": compare(cubes[(a, 48)], cubes[(b, 48)])
              for a, b in zip((300, 150, 75, 30), (150, 75, 30, 15))}
    checks.update({f"frequency_{a}_to_{b}": compare(cubes[(15, a)], cubes[(15, b)])
                   for a, b in ((48, 96), (96, 192))})
    checks.update({f"first_two_hours_dt_{a}_to_{b}": compare(cubes[(a, 48)][:25], cubes[(b, 48)][:25])
                   for a, b in ((15, 5), (5, 2), (2, 1))})
    checks["hours_one_to_two_dt_2_to_1"] = compare(cubes[(2, 48)][12:], cubes[(1, 48)][12:])
    checks["hours_one_to_two_dt_2_to_1"]["scope_note"] = (
        "Follow-up diagnostic for 1–2 hours, examined after the full 20-minute–2-hour check failed for mean period. Earlier times remain unresolved.")
    rows = []
    for (dt, nf), cube in cubes.items():
        for minute in (20, 60, 120, 180, 360):
            if minute // 5 >= len(cube):
                continue
            for x in (10, 25, 50, 100):
                v = cube[minute // 5, x]
                rows.append(dict(step_s=dt, frequency_intervals=nf, duration_min=minute,
                                 fetch_km=x, hs_m=float(v[2]), peak_period_s=float(v[3]), mean_period_s=float(v[4])))
    return dict(schema="terluna.climate.wave-early-growth/1",
                evidence="Numerical refinement of a six-hour, uniform p90-wind episode on the original deep strip. Empirical validation remains open.",
                reading_rule="Comparisons use matched time, fetch and wind, with reference Hs at least 1 cm. Five percent is the chosen numerical tolerance for this wind case; physical accuracy requires separate validation.",
                units=dict(time="s or labelled minutes", length="m or labelled km", period="s", frequency="Hz"),
                producer=dict(domain="climate", files={str(p.relative_to(ROOT)): sha256(p) for p in
                              (Path(__file__), Path(__file__).with_name("model.py"))}),
                inputs=dict(wind=forcing, build_sha256=sha256(build_root / "build.json")),
                checks=checks, cases=cases, rows=rows, runs=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, default=ROOT / "research/runs/waves/early_growth")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("results") / "early_growth.json")
    args = parser.parse_args()
    result = experiment(args.build_root, args.run_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["checks"], indent=2))


if __name__ == "__main__":
    main()
