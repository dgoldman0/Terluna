"""Nest a coastal sector in the twelve-hour Smythii–Marginis wind episode."""
from __future__ import annotations

import argparse
from datetime import timedelta
import json
from pathlib import Path

import numpy as np

from climate.waves.model import ORIGIN, ROOT, sha256, wind_inputs
from climate.waves.pilot import BasinCase, array_text, basin_files, execute, load_basin, read_basin
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
RUNS = ROOT / "research/runs/waves/coastal"
TERRAIN = RUNS / "terrain.npz"


def load_terrain(path=TERRAIN):
    with np.load(path, allow_pickle=False) as f:
        data = {k: f[k].copy() for k in f.files}
    data["metadata"] = json.loads(data["metadata"].item())
    if data["metadata"]["schema"] != "terluna.geography.coastal-grid/1":
        raise ValueError("Unsupported coastal terrain")
    x, y, depth = data["longitude_deg"], data["latitude_deg"], data["depth_m"]
    if depth.shape != (len(y), len(x)) or not np.isfinite(depth).all() or np.any(depth[data["wet"]] <= 0):
        raise ValueError("Invalid coastal grid")
    return data


def parent_files(terrain, boundary_minutes=15):
    case = BasinCase("parent", stride=2, frequency_intervals=48)
    basin = load_basin(stride=2)
    files = basin_files(case, basin, wind_inputs()[0]["p90"])
    x0, y0, width, height = terrain["metadata"]["bounds_deg"]
    start = ORIGIN.strftime("%Y%m%d.%H%M%S")
    nesting = (f"NGRID 'COAST' {x0} {y0} 0 {width} {height} {round(width*2)} {round(height*2)}\n"
               f"NESTOUT 'COAST' 'coast.nest' OUTPUT {start} {boundary_minutes} MIN\n")
    files["INPUT"] = files["INPUT"].replace("COMPUTE NONSTAT", nesting + "COMPUTE NONSTAT")
    return case, basin, files


def nested_files(terrain, boundary, stride=1, directions=36, step_s=75):
    if stride not in (1, 2, 4, 8) or directions not in (36, 72) or step_s <= 0 or 3600 % step_s:
        raise ValueError("Unsupported coastal refinement")
    x, y = terrain["longitude_deg"][::stride], terrain["latitude_deg"][::stride]
    bottom = terrain["depth_m"]
    xi, yi = terrain["longitude_deg"], terrain["latitude_deg"]
    width, height = x[-1]-x[0], y[-1]-y[0]
    start = ORIGIN.strftime("%Y%m%d.%H%M%S")
    end = (ORIGIN + timedelta(hours=12)).strftime("%Y%m%d.%H%M%S")
    ratio = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    wind = wind_inputs()[0]["p90"]
    text = f"""PROJECT 'Smythii coast' '005'
MODE NONSTAT TWODIMENSIONAL
COORDINATES SPHERICAL {MOON_RADIUS:.12g} CCM
SET GRAV {MOON_SURFACE_GRAVITY:.12g} RHO 1025
OUTPUT OPTIONS TABLE 16
CGRID REG {x[0]:.12g} {y[0]:.12g} 0 {width:.12g} {height:.12g} {len(x)-1} {len(y)-1} &
 CIRCLE {directions} {0.03*ratio:.12g} {3*ratio:.12g} 48
INPGRID BOTTOM {xi[0]:.12g} {yi[0]:.12g} 0 {len(xi)-1} {len(yi)-1} {xi[1]-xi[0]:.12g} {yi[1]-yi[0]:.12g}
READINP BOTTOM 1 'depth.bot' 3 0 FREE
BOUNDNEST1 NEST 'boundary.nest' CLOSED
WIND {wind:.12g} 0
GEN3 KOMEN DRAG WU AGROW
BREAKING CONSTANT 1.0 0.73
PROP BSBT
INIT ZERO
QUANTITY TSEC REF {start}
TABLE 'COMPGRID' NOHEAD 'coastal.tbl' TSEC XP YP HS RTP TM01 DEP DIR QB &
 OUTPUT {start} 1 HR
COMPUTE NONSTAT {start} {step_s} SEC {end}
STOP
"""
    return {"INPUT": text, "depth.bot": array_text(bottom), "boundary.nest": boundary}


def read_coastal(path, terrain, stride):
    x, y = terrain["longitude_deg"][::stride], terrain["latitude_deg"][::stride]
    data = np.loadtxt(path, ndmin=2)
    if data.shape != (13 * len(x) * len(y), 9) or not np.isfinite(data).all():
        raise ValueError("Incomplete coastal output")
    cube = data.reshape(13, len(y), len(x), 9)
    wet = terrain["wet"][::stride, ::stride]
    xx, yy = np.meshgrid(x, y)
    points = np.column_stack((xx[wet], yy[wet]))
    values = cube[:, wet]
    if not np.allclose(values[..., :3], np.stack(np.broadcast_arrays(np.arange(13)[:, None]*3600,
                                    points[None, :, 0], points[None, :, 1]), axis=-1), atol=1e-4, rtol=1e-6):
        raise ValueError("Coastal locations or times differ")
    if not np.allclose(values[..., 6], terrain["depth_m"][::stride, ::stride][wet], atol=.01, rtol=6e-5):
        raise ValueError("Coastal wet-node depths differ")
    if np.any(values[..., 3] < 0) or np.any(values[1:, :, 4:6] <= 0) or np.any((values[..., 8] < 0) | (values[..., 8] > 1)):
        raise ValueError("Invalid coastal wave diagnostics")
    return cube, wet


def compare(coarse, fine, coarse_wet, fine_wet, factor=2):
    matched = fine[:, ::factor, ::factor]
    if coarse.shape != matched.shape:
        raise ValueError("Coastal grids must share nodes")
    wet = coarse_wet & fine_wet[::factor, ::factor]
    # Keep a two-coarse-cell buffer from the prescribed nesting boundary.
    interior = np.zeros(wet.shape, bool)
    interior[2:-2, 2:-2] = True
    mask = wet & interior & (coarse[-1, ..., 3] >= .1) & (matched[-1, ..., 3] >= .1)
    if not mask.any():
        raise ValueError("Empty interior comparison")
    result = dict(matched_interior_nodes=int(mask.sum()), boundary_buffer_coarse_cells=2,
                  minimum_hs_m=.1, relative_denominator="refined run")
    for col, name in ((3, "hs"), (5, "tm01")):
        error = np.abs(coarse[-1, ..., col][mask] / matched[-1, ..., col][mask] - 1)
        result[name+"_relative_max"] = float(error.max())
        result[name+"_relative_p95"] = float(np.quantile(error, .95))
        result[name+"_absolute_max"] = float(np.abs(coarse[-1, ..., col][mask]-matched[-1, ..., col][mask]).max())
    result["within_five_percent"] = max(result["hs_relative_max"], result["tm01_relative_max"]) < .05
    return result


def experiment(build_root, run_root=RUNS):
    terrain = load_terrain()
    executable = build_root / "coupled_air/swan.exe"
    manifest = json.loads((build_root / "coupled_air.json").read_text())
    if sha256(executable) != manifest["build"]["executable_sha256"]:
        raise ValueError("Executable differs from the coupled-air build")
    case, basin, files = parent_files(terrain)
    parent = run_root / "parent"
    records = {"parent": execute(executable, parent, files, ["basin.tbl", "offshore.spc", "coast.nest"])}
    read_basin(parent / "basin.tbl", case, basin)
    print("parent complete", records["parent"]["elapsed_wall_s"], flush=True)
    boundary = (parent / "coast.nest").read_text()
    cases, cubes, masks = [], {}, {}
    for stride, directions in ((4, 36), (2, 36), (1, 36), (2, 72)):
        name = f"coast_s{stride}_d{directions}"
        files = nested_files(terrain, boundary, stride, directions)
        directory = run_root / name
        records[name] = execute(executable, directory, files, ["coastal.tbl"], timeout_s=1500 if stride == 1 else 900)
        cube, wet = read_coastal(directory / "coastal.tbl", terrain, stride)
        cubes[name], masks[name] = cube, wet
        cases.append(dict(name=name, stride=stride, direction_bins=directions,
                          spacing_deg=stride/16, wet_nodes=int(wet.sum()),
                          maximum_final_hs_m=float(cube[-1, ..., 3][wet].max()),
                          maximum_breaking_fraction=float(cube[1:, ..., 8][:, wet].max()),
                          longitude_deg=terrain["longitude_deg"][::stride].tolist(),
                          latitude_deg=terrain["latitude_deg"][::stride].tolist(),
                          wet=wet.tolist(), final_values=cube[-1].tolist()))
        print(name, "complete", records[name]["elapsed_wall_s"], flush=True)
    checks = {}
    for a, b, factor in (("coast_s4_d36", "coast_s2_d36", 2), ("coast_s2_d36", "coast_s1_d36", 2),
                          ("coast_s2_d36", "coast_s2_d72", 1)):
        checks[a+"_to_"+b] = compare(cubes[a], cubes[b], masks[a], masks[b], factor)
    return dict(schema="terluna.climate.coastal-resolution/1",
                evidence="Twelve-hour, uniform eastward wind episode nested from the 0.5-degree basin into flooded eastern Smythii terrain. Lunar wind input and breaking calibration remain open.",
                reading_rule="All nested cases use the same 16-pixel/degree terrain and 15-minute full directional boundary spectra. Computational grids sample nested subsets of native terrain nodes. Directional refinement changes the child grid; the parent retains 36 directions. Comparison errors cover final common wet nodes at least two coarse cells inside the boundary.",
                units=dict(columns=["time_s", "longitude_deg", "latitude_deg", "hs_m", "peak_period_s", "mean_period_s", "depth_m", "direction_deg", "breaking_fraction"]),
                producer=dict(domain="climate", files={str(p.relative_to(ROOT)): sha256(p) for p in
                    (Path(__file__), HERE/"pilot.py", HERE/"model.py", ROOT/"shared/constants.json")}),
                inputs=dict(terrain_sha256=sha256(TERRAIN), terrain_metadata=terrain["metadata"],
                            build_manifest_sha256=sha256(build_root/"coupled_air.json"), wind=wind_inputs()[1]),
                checks=checks, cases=cases, runs=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, default=RUNS)
    parser.add_argument("--output", type=Path, default=HERE/"results/coastal.json")
    args = parser.parse_args()
    product = experiment(args.build_root, args.run_root)
    args.output.write_text(json.dumps(product, separators=(",", ":")) + "\n")
    print(json.dumps(product["checks"], indent=2))


if __name__ == "__main__":
    main()
