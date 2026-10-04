"""Run or reuse the bounded cycling experiment matrix, with at most two workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
import subprocess
import sys

from .cycling_analysis import read_raw
from .cycling_search import HERE, ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, choices=(1, 2), default=1)
    parser.add_argument("--case", help="Run just this named case")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    cases = json.loads((HERE/"cycling_scenario.json").read_text())["cases"]
    if args.case:
        cases = [c for c in cases if c["name"] == args.case]
        if not cases:
            raise ValueError("Unknown cycling case")
    def run(case):
        if not args.force:
            try:
                _, metadata = read_raw(case["name"])
                if all(metadata["config"].get(k) == v for k, v in case.items() if k != "name"):
                    return {"name": case["name"], "reused_verified_raw": True}
            except (FileNotFoundError, ValueError):
                pass
        cmd = [sys.executable, "-u", "-m", "research.studies.solar_shield_array.cycling_feedback",
               "--out-name", case["name"]]
        for key, value in case.items():
            if key == "name":
                continue
            cmd.append("--"+key.replace("_", "-"))
            if value is not True:
                cmd.append(str(value))
        env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
        result = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True, timeout=300)
        if result.returncode:
            raise RuntimeError(f"{case['name']}: {result.stderr}")
        return {"name": case["name"], "output": result.stdout}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for result in as_completed([pool.submit(run, case) for case in cases]):
            print(json.dumps(result.result()), flush=True)


if __name__ == "__main__":
    main()
