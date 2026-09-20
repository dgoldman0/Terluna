#!/usr/bin/env python3
"""Check imported bytes and existing numerical baselines, with explicit blocked inputs."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import runpy
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_value(actual, expected, location="root"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            raise ValueError(f"Fields differ at {location}")
        for key in expected:
            same_value(actual[key], expected[key], f"{location}.{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise ValueError(f"List differs at {location}")
        for i, (left, right) in enumerate(zip(actual, expected)):
            same_value(left, right, f"{location}[{i}]")
    elif isinstance(expected, bool) or expected is None or isinstance(expected, str):
        if actual != expected:
            raise ValueError(f"Value differs at {location}: {actual!r} vs {expected!r}")
    elif not math.isclose(float(actual), float(expected), rel_tol=1e-7, abs_tol=1e-12):
        raise ValueError(f"Numerical mismatch at {location}: {actual} vs {expected}")


def check() -> dict:
    if sys.flags.optimize:
        raise ValueError("Run without -O: the historical model includes assert-based checks.")
    provenance = json.loads((ROOT / "research/provenance.json").read_text())
    retained = provenance["retained"]
    for item in retained:
        path = (ROOT / item["path"]).resolve()
        if not path.is_relative_to(ROOT) or digest(path) != item["sha256"]:
            raise ValueError(f"Imported bytes differ: {item['path']}")
    result = {"checked_at_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "artifact integrity and numerical implementation only",
              "scientific_validation": False, "manuscript_clearance": False,
              "retained_files": len(retained), "integrity": "PASS",
              "python": platform.python_version(),
              "dependencies": {n: importlib.metadata.version(n) for n in ("numpy", "scipy", "matplotlib")},
              "runner_sha256": digest(Path(__file__)), "checks": {}}
    geometry = ROOT / "illumination/geometry.py"
    if geometry.read_bytes() != (ROOT / "ensemble/planning/twilight_diagnostic.py").read_bytes():
        raise ValueError("Illumination diagnostic differs from its imported source")
    actual = runpy.run_path(str(geometry))["diagnostic"]()
    expected = json.loads((ROOT / "ensemble/planning/twilight_diagnostic.json").read_text())
    same_value(actual, expected, "illumination")
    result["checks"]["illumination"] = {"status": "PASS", "scope": "angular arithmetic only", "output": actual}
    model = ROOT / "research/baselines/feasibility/model.py"
    reference = ROOT / "research/baselines/feasibility/reference.json"
    with tempfile.TemporaryDirectory(prefix="terluna-check-") as directory:
        work = Path(directory)
        out = work / "feasibility"
        completed = subprocess.run([sys.executable, str(model), "--out", str(out)],
                                   capture_output=True, text=True, timeout=180)
        if completed.returncode:
            raise ValueError(f"Feasibility run failed: {completed.stderr[-4000:]}")
        actual = json.loads((out / "results.json").read_text())
        expected = json.loads(reference.read_text())
        same_value({k: actual[k] for k in expected}, expected, "feasibility")
        tables = [i for i in retained if i.get("archive") == "Lunar_Terraforming_Model.zip"
                  and i["member"].endswith(".csv")]
        counts = {}
        for item in tables:
            generated = out / Path(item["member"]).name
            with generated.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            with (ROOT / item["path"]).open(newline="") as stream:
                old = list(csv.DictReader(stream))
            if len(rows) != len(old):
                raise ValueError(f"Row count changed: {item['path']}")
            for i, (left, right) in enumerate(zip(rows, old)):
                if left.keys() != right.keys():
                    raise ValueError(f"CSV fields changed: {item['path']}")
                for key in left:
                    try:
                        a, b = float(left[key]), float(right[key])
                    except (ValueError, TypeError):
                        a, b = left[key], right[key]
                    same_value(a, b, f"{item['path']}[{i}].{key}")
            counts[Path(item["path"]).name] = len(rows)
        result["checks"]["feasibility"] = {"status": "PASS", "model_sha256": digest(model),
                  "reference_sha256": digest(reference), "tables_checked": counts,
                  "numerical_tests": actual["verification"], "baseline": actual["baseline"],
                  "generated_plots": "four historical plots generated in a temporary directory; no visual review"}
        protection = ROOT / "protection"
        inputs = json.loads((protection / "inputs.json").read_text())
        missing = []
        for item in inputs["files"]:
            path = protection / inputs["input_directory"] / item["name"]
            if not path.is_file() or digest(path) != item["sha256"]:
                missing.append(item["name"])
        if missing:
            result["checks"]["protection"] = {"status": "BLOCKED", "reason": "external inputs absent or changed",
                                               "inputs": missing, "verifier_executed": False}
        else:
            target = work / "protection"
            target.mkdir()
            for name in ("model.py", "verify.py"):
                shutil.copy2(protection / name, target / name)
            (target / "results").mkdir()
            shutil.copy2(protection / "results/results.json", target / "results/results.json")
            shutil.copytree(protection / "sources", target / "sources")
            completed = subprocess.run([sys.executable, str(target / "verify.py")],
                                       cwd=target, capture_output=True, text=True, timeout=180)
            if completed.returncode:
                raise ValueError(f"Protection verifier failed: {completed.stderr[-4000:]}")
            result["checks"]["protection"] = {"status": "PASS", "verifier_executed": True,
                "model_sha256": digest(protection / "model.py"),
                "verifier_sha256": digest(protection / "verify.py"),
                "reference_sha256": digest(protection / "results/results.json"),
                "output": json.loads((target / "results/validation.json").read_text()),
                "full_optical_optimization_rerun": False}
    result["not_executed"] = ["April simulator", "regional climate or hydrology", "biological or habitation simulations",
                              "global plasma or dose transport", "full ephemeris optimization"]
    result["result"] = "PASS_WITH_BLOCKED_INPUTS" if any(v["status"] == "BLOCKED" for v in result["checks"].values()) else "PASS"
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "research/runs/checks.json")
    parser.add_argument("--require-inputs", action="store_true", help="Return exit 2 if protection inputs are missing")
    args = parser.parse_args()
    try:
        result = check()
    except (OSError, ValueError, KeyError, subprocess.TimeoutExpired) as exc:
        print(f"CHECK FAILED: {exc}", file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 2 if args.require_inputs and result["result"] != "PASS" else 0


if __name__ == "__main__":
    raise SystemExit(main())
