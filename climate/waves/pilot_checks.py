"""Audit spectra, solver stopping and resolution in completed pilot outputs."""
from __future__ import annotations

import argparse
from datetime import timedelta
import json
import math
from pathlib import Path
import re

import numpy as np

from climate.waves.model import ORIGIN, ROOT, sha256


def spectrum_diagnostics(path, expected_date):
    lines = path.read_text().splitlines()
    dates = [line.split()[0] for line in lines if re.match(r"^\d{8}\.\d{6}\s", line)]
    if dates != [expected_date]:
        raise ValueError("Offshore spectrum is missing its expected final time")
    index = next(i for i, s in enumerate(lines) if s.split() and s.split()[0] == "AFREQ")
    n = int(lines[index + 1].split()[0])
    frequency = np.array([float(s.split()[0]) for s in lines[index + 2:index + 2 + n]])
    locations = [i for i, s in enumerate(lines) if s.split() and s.split()[0] == "LOCATION"]
    if len(locations) != 1 or n < 3 or len(frequency) != n or not np.isfinite(frequency).all() or np.any(np.diff(frequency) <= 0) or frequency[0] <= 0:
        raise ValueError("Invalid offshore frequency grid or location count")
    index = next(i for i, s in enumerate(lines) if s.startswith("VaDens"))
    exception = float(lines[index + 2].split()[0])
    values = np.array([[float(x) for x in s.split()] for s in lines[locations[0] + 1:locations[0] + 1 + n]])
    if values.shape != (n, 3) or not np.isfinite(values).all():
        raise ValueError("Incomplete offshore spectral quantities")
    density = values[:, 0].copy()
    omitted = density == exception
    if np.any((density < 0) & ~omitted) or np.any(values[~omitted, 2] < 0):
        raise ValueError("Invalid variance or directional spread")
    density[omitted] = 0
    m0 = float(np.trapezoid(density, frequency))
    if m0 <= 0:
        raise ValueError("Empty offshore spectrum")
    low = float(np.trapezoid(density[:3], frequency[:3]) / m0)
    high = float(np.trapezoid(density[-3:], frequency[-3:]) / m0)
    peak_at_edge = int(np.argmax(density)) in (0, n - 1)
    return dict(frequencies=n, resolved_hs_m=4 * math.sqrt(m0),
                resolved_tm01_s=float(m0 / np.trapezoid(frequency * density, frequency)),
                low_two_interval_variance_fraction=low, high_two_interval_variance_fraction=high,
                peak_at_edge=peak_at_edge,
                edges_pass=bool(low < .001 and high < .01 and not peak_at_edge),
                omitted_variance_upper_bound_m2=float(np.trapezoid(omitted * (2 * math.pi * 1e-12), frequency)),
                note="Integrals cover the resolved printed spectrum. SWAN's Hs and Tm01 also include a diagnostic high-frequency tail. SPEC1D retains variance and mean direction/spread by frequency; transferring multiple directional peaks requires a full two-dimensional spectrum.")


def audit(product_path, run_root):
    product = json.loads(product_path.read_text())
    if product["schema"] != "terluna.climate.wave-basin-pilot/1":
        raise ValueError("Unsupported pilot product")
    spectra = {}
    for case in product["basin_cases"]:
        expected = (ORIGIN + timedelta(hours=case["hours"])).strftime("%Y%m%d.%H%M%S")
        path = run_root / case["name"] / "offshore.spc"
        if sha256(path) != product["runs"][case["name"]]["output_sha256"]["offshore.spc"]:
            raise ValueError("Offshore spectrum differs from the recorded pilot output")
        spectra[case["name"]] = spectrum_diagnostics(path, expected)
    stationary = {}
    slopes = {s["name"]: s for s in product["slopes"]}
    for name in slopes:
        path = run_root / name / "PRINT"
        if sha256(path) != product["runs"][name]["output_sha256"]["PRINT"]:
            raise ValueError("Stationary solver log differs from the recorded pilot output")
        printed = path.read_text()
        matches = re.findall(r"accuracy OK in\s+([\d.]+) % of wet grid points \(\s*([\d.]+) % required\)", printed)
        if not matches:
            raise ValueError("Missing stationary solver convergence diagnostic")
        actual, required = map(float, matches[-1])
        stationary[name] = dict(final_percent=actual, required_percent=required, passed=actual >= required)
    coarse = np.asarray(slopes["slope_1in50_n400_g73"]["rows"])
    fine = np.asarray(slopes["slope_1in50_n800_g73"]["rows"])[::2]
    if coarse.shape != fine.shape or not np.allclose(coarse[:, 4], fine[:, 4], atol=1e-5):
        raise ValueError("Slope refinement depths are not matched")
    mask = fine[:, 1] >= .1
    slope_error = float(np.max(np.abs(coarse[mask, 1] / fine[mask, 1] - 1)))
    slope_check = dict(hs_relative_max=slope_error, hs_absolute_max_m=float(np.max(np.abs(coarse[:, 1] - fine[:, 1]))),
                       passed=slope_error < .05, comparison_min_hs_m=.1)
    early_path = product_path.with_name("early_growth.json")
    early = json.loads(early_path.read_text())
    short = early["checks"]["first_two_hours_dt_2_to_1"]
    resolution = {name: dict(**check, final_hs_within_five_percent=check["hs_final_relative_max"] < .05,
                            final_hs_tm01_within_five_percent=max(check["hs_final_relative_max"], check["mean_period_final_relative_max"]) < .05)
                  for name, check in product["checks"].items()}
    breaking = {}
    for name, row in slopes.items():
        values = np.asarray(row["rows"])
        indices = np.flatnonzero(values[:, 6] >= .01)
        if not len(indices):
            raise ValueError("No one-percent-breaking region in a stated surf experiment")
        first = values[indices[0]]
        breaking[name] = dict(first_one_percent_breaking_depth_m=float(first[4]),
                              first_one_percent_breaking_hs_m=float(first[1]),
                              maximum_hs_where_at_least_one_percent_break_m=float(values[indices, 1].max()),
                              note="Qb >= 0.01 marks the chosen onset-of-breaking diagnostic; Hs describes the wave spectrum.")
    return dict(schema="terluna.climate.wave-pilot-checks/1",
                evidence="Numerical diagnostics for the wave pilot. Empirical lunar validation and manuscript/source review remain open.",
                reading_rule="Failed local basin resolution bounds restrict geographic interpretation. A converged ideal slope is still conditional on its offshore spectrum, assumed slope and Earth-derived breaking parameters. Early-time dt and frequency checks cover separate runs and only the 4.001 m/s wind.",
                producer=dict(domain="climate", path="climate/waves/pilot_checks.py", sha256=sha256(Path(__file__))),
                inputs=dict(pilot_sha256=sha256(product_path), early_growth_sha256=sha256(early_path)),
                units=dict(relative_errors="dimensionless", heights_and_depths="m", periods="s"),
                spectra=spectra, stationary_stopping=stationary, slope_grid_refinement=slope_check,
                breaking_diagnostics=breaking,
                basin_resolution=resolution, early_first_two_hours_2_to_1_s=short,
                early_hours_one_to_two_2_to_1_s=early["checks"]["hours_one_to_two_dt_2_to_1"],
                frequency_spacing=early["checks"]["frequency_96_to_192"],
                spectral_edges_pass=all(s["edges_pass"] for s in spectra.values()),
                stationary_stopping_pass=all(s["passed"] for s in stationary.values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--product", type=Path, default=ROOT / "climate/waves/results/pilot.json")
    parser.add_argument("--run-root", type=Path, default=ROOT / "research/runs/waves/pilot")
    parser.add_argument("--output", type=Path, default=ROOT / "climate/waves/results/pilot_checks.json")
    args = parser.parse_args()
    result = audit(args.product, args.run_root)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("spectra", "stationary_stopping")}, indent=2))


if __name__ == "__main__":
    main()
