"""Evolving eastern Smythii spectra, retaining the full first-cycle spin-up."""
from __future__ import annotations

import argparse
from datetime import timedelta
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from climate.waves.coastal import load_terrain
from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.pilot import array_text, execute, load_basin, BasinCase, read_basin
from climate.waves.shore import REFERENCE, RUNS as SHORE_RUNS
from climate.waves.shore_region import perimeter

RUNS = ROOT / "research/runs/waves/shore_history"
END_HOUR = 1419


def parent_case(executable, run_root=RUNS):
    """Replay the established sea and retain hourly directional source spectra."""
    record = json.loads((REFERENCE / "run.json").read_text())
    if sha256(executable) != record["identity"]["executable_sha256"]:
        raise ValueError("Continuous replay requires the recorded executable")
    files = {}
    for name, digest in record["identity"]["inputs"].items():
        path = REFERENCE / name
        if sha256(path) != digest:
            raise ValueError(f"Reference input changed: {name}")
        files[name] = path
    terrain = load_terrain(SHORE_RUNS / "regional_terrain.npz")
    xy, wet = perimeter(terrain, 16)
    basin = load_basin()
    scale = np.array([np.cos(np.deg2rad(np.mean(xy[:,1]))), 1.])
    _, indices = cKDTree(basin["points"]*scale).query(xy[wet]*scale, k=4)
    selected = np.unique(indices)
    files["field.xy"] = array_text(basin["points"][selected])
    start = ORIGIN.strftime("%Y%m%d.%H%M%S")
    extra = ("POINTS 'field' FILE 'field.xy'\n"
             f"SPEC 'field' SPEC2D ABS 'field.spc' OUTPUT {start} 1 HR\n")
    files["INPUT"] = (REFERENCE / "INPUT").read_text().replace("COMPUTE NONSTAT", extra+"COMPUTE NONSTAT")
    directory = Path(run_root) / "parent"
    run = execute(executable, directory, files,
                  ["basin.tbl", "wind.tbl", "offshore.spc", "diagnostics.spc", "field.spc"],
                  timeout_s=7200, threads=2)
    case = BasinCase("history_parent", stride=4, step_s=150, frequency_intervals=48, hours=END_HOUR)
    actual = read_basin(directory/"basin.tbl", case, basin)
    expected = np.asarray(json.loads((REFERENCE/"product.json").read_text())["values"])
    np.testing.assert_array_equal(actual, expected)
    product_path=directory/'product.json'
    if product_path.exists():
        previous=json.loads(product_path.read_text())
        if (previous['schema']!='terluna.climate.shore-history-parent/1' or
                previous['run']!=run or previous['reference_sha256']!=sha256(REFERENCE/'product.json') or
                sha256(directory/'producer_source.py')!=previous['producer_sha256']):
            raise ValueError('Completed basin replay metadata differs from its source record')
        return previous
    result = dict(schema="terluna.climate.shore-history-parent/1", run=run,
                  evidence="Continuous basin replay with hourly directional spectra for the regional boundary.",
                  reading_rule="The first lunar solar cycle supplies spin-up; occurrence statistics use the second cycle.",
                  hours=END_HOUR, spectrum_interval_hours=1, source_nodes=len(selected),
                  reference_sha256=sha256(REFERENCE/"product.json"),
                  reference_values_compared=int(actual.size), reference_values_exact=True,
                  producer_sha256=sha256(Path(__file__)))
    (directory/"producer_source.py").write_bytes(Path(__file__).read_bytes())
    product_path.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result),flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("parent",))
    parser.add_argument("--executable", type=Path, required=True)
    args = parser.parse_args()
    parent_case(args.executable)


if __name__ == "__main__":
    main()
