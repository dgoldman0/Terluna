"""A bounded Smythii–Marginis basin and idealised coastal-slope experiments.

The atlas supplies flooded lunar ground. The experiments use uniform wind
episodes, unreworked shores, omitted bottom friction and assumed breaking
parameters to calculate basin spectra and stationary slope profiles.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
from datetime import timedelta
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import resource
import subprocess
import time

import numpy as np

from climate.waves.build import build_coupled
from climate.waves.model import ORIGIN, ROOT, sha256, wind_inputs
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ATLAS = ROOT / "geography/products/atlas_28pct_4ppd.npz"


def array_text(values):
    stream = io.StringIO()
    np.savetxt(stream, values, fmt="%.10g")
    return stream.getvalue()


def load_basin(path=ATLAS, stride=4, body=874):
    """Sample nested subsets of atlas nodes, retaining their depths and land mask."""
    if stride not in (1, 2, 4):
        raise ValueError("Use a nested atlas stride of 1, 2 or 4")
    with np.load(path, allow_pickle=False) as source:
        metadata = json.loads(source["metadata"].item())
        if metadata["schema"] != "terluna.geography.atlas-grid/1":
            raise ValueError("Unsupported atlas grid")
        mask = source["water_label"][::stride, ::stride] == body
        rows, cols = np.where(mask)
        if not len(rows):
            raise ValueError("Basin is absent")
        # One dry margin beyond every side, retained on the global nested grid.
        rs = slice(max(0, rows.min() - 1) * stride, (rows.max() + 2) * stride, stride)
        cs = slice(max(0, cols.min() - 1) * stride, (cols.max() + 2) * stride, stride)
        lon, lat = source["lon_deg"][cs], source["lat_deg"][rs][::-1]
        wet = (source["water_label"][rs, cs] == body)[::-1]
        height = source["height_m"][rs, cs][::-1]
        physical_depth = metadata["sea_level_m"] - height
        # Finite dry-node elevations avoid exception-value interpolation
        # removing wet nodes through roundoff at a land/water boundary.
        # Other disconnected water bodies in the crop are held dry explicitly.
        depth = np.where(wet, physical_depth, np.minimum(physical_depth, -1.))
    if np.any(depth[wet] <= 0) or any(wet[side].any() for side in (0, -1)) or wet[:, (0, -1)].any():
        raise ValueError("Invalid basin depth or open crop boundary")
    xx, yy = np.meshgrid(lon, lat)
    points = np.column_stack((xx[wet], yy[wet]))
    # Keep the same offshore node in every nested run.
    choice = np.argmin((points[:, 0] - 89.125)**2 + (points[:, 1] + 1.125)**2)
    return dict(lon=lon, lat=lat, wet=wet, depth=depth, points=points, offshore=points[choice],
                offshore_index=int(choice), metadata=metadata, source_sha256=sha256(path))


def execute(executable, directory, files, outputs, timeout_s=900):
    """Hash all run inputs and required outputs, with bounded serial execution."""
    identity = dict(executable_sha256=sha256(executable),
                    inputs={name: hashlib.sha256(data.encode()).hexdigest() for name, data in files.items()})
    record_path = directory / "run.json"
    if directory.exists():
        if not record_path.exists():
            raise FileExistsError(f"Incomplete run retained for inspection: {directory}")
        record = json.loads(record_path.read_text())
        if record["identity"] != identity or set(record["output_sha256"]) != set(outputs) | {"PRINT", "norm_end"}:
            raise ValueError("Cached run identity differs")
        for name, expected in {**identity["inputs"], **record["output_sha256"]}.items():
            if sha256(directory / name) != expected:
                raise ValueError(f"Cached file hash differs: {name}")
        return record
    directory.mkdir(parents=True)
    for name, data in files.items():
        (directory / name).write_text(data)
    def limits():
        resource.setrlimit(resource.RLIMIT_CPU, (timeout_s, timeout_s))
        resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
    start = time.monotonic()
    with (directory / "console.log").open("w") as log:
        result = subprocess.run([str(executable.resolve())], cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                                env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"),
                                timeout=timeout_s, preexec_fn=limits)
    printed = (directory / "PRINT").read_text() if (directory / "PRINT").exists() else ""
    if result.returncode or not (directory / "norm_end").exists() or re.search(r"\*\*\s+(Error|Severe|Fatal)", printed, re.I):
        raise RuntimeError(f"SWAN failed; inspect {directory}")
    record = dict(schema="terluna.climate.swan-pilot-run/1", identity=identity,
                  producer_sha256=sha256(Path(__file__)), elapsed_wall_s=time.monotonic() - start,
                  threads=1, address_space_limit_bytes=2 * 1024**3, cpu_limit_s=timeout_s,
                  warnings=[line.strip() for line in printed.splitlines() if "warning" in line.lower()],
                  output_sha256={name: sha256(directory / name) for name in (*outputs, "PRINT", "norm_end")})
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    return record


@dataclass(frozen=True)
class BasinCase:
    name: str
    stride: int = 4
    step_s: int = 75
    frequency_intervals: int = 96
    hours: int = 12
    turn: bool = False

    def __post_init__(self):
        if not re.fullmatch(r"[a-z0-9_]+", self.name) or self.stride not in (1, 2, 4):
            raise ValueError("Invalid name or atlas stride")
        if self.step_s <= 0 or 3600 % self.step_s or self.hours < 1 or self.frequency_intervals < 12:
            raise ValueError("Invalid time or frequency grid")
        if self.turn and self.hours < 7:
            raise ValueError("Turning episode needs at least seven hours")


def basin_files(case, basin, wind):
    lon, lat = basin["lon"], basin["lat"]
    x0, y0, width, height = lon[0], lat[0], lon[-1] - lon[0], lat[-1] - lat[0]
    nx, ny, delta = len(lon) - 1, len(lat) - 1, lon[1] - lon[0]
    start, end = ORIGIN.strftime("%Y%m%d.%H%M%S"), (ORIGIN + timedelta(hours=case.hours)).strftime("%Y%m%d.%H%M%S")
    ratio = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    forcing = []
    for hour in range(case.hours + 1):
        u, v = (0., wind) if case.turn and hour >= 6 else (wind, 0.)
        forcing.extend((array_text(np.full((2, 2), u)), array_text(np.full((2, 2), v))))
    x, y = basin["offshore"]
    text = f"""PROJECT 'Lunar basin' '002'
MODE NONSTAT TWODIMENSIONAL
COORDINATES SPHERICAL {MOON_RADIUS:.12g} CCM
SET GRAV {MOON_SURFACE_GRAVITY:.12g} RHO 1025
OUTPUT OPTIONS TABLE 16
CGRID REG {x0:.12g} {y0:.12g} 0 {width:.12g} {height:.12g} {nx} {ny} &
 CIRCLE 36 {0.03*ratio:.12g} {3*ratio:.12g} {case.frequency_intervals}
INPGRID BOTTOM {x0:.12g} {y0:.12g} 0 {nx} {ny} {delta:.12g} {delta:.12g}
READINP BOTTOM 1 'depth.bot' 3 0 FREE
INPGRID WIND {x0:.12g} {y0:.12g} 0 1 1 {width:.12g} {height:.12g} &
 NONSTAT {start} 1 HR {end}
READINP WIND 1 'wind.dat' 3 0 0 0 FREE
GEN3 KOMEN DRAG WU AGROW
BREAKING CONSTANT 1.0 0.73
PROP BSBT
INIT ZERO
QUANTITY TSEC REF {start}
TABLE 'COMPGRID' NOHEAD 'basin.tbl' TSEC XP YP HS RTP TM01 DEP DIR QB &
 OUTPUT {start} 1 HR
POINTS 'offshore' {x:.12g} {y:.12g}
SPEC 'offshore' SPEC1D ABS 'offshore.spc' OUTPUT {end} 1 HR
COMPUTE NONSTAT {start} {case.step_s} SEC {end}
STOP
"""
    return {"INPUT": text, "depth.bot": array_text(basin["depth"]), "points.xy": array_text(basin["points"]),
            "wind.dat": "".join(forcing)}


def read_basin(path, case, basin):
    data = np.loadtxt(path, ndmin=2)
    ny, nx = basin["depth"].shape
    n = nx * ny
    if data.shape != ((case.hours + 1) * n, 9) or not np.isfinite(data).all():
        raise ValueError("Incomplete basin output")
    full = data.reshape(case.hours + 1, ny, nx, 9)
    cube = full[:, basin["wet"]]
    if not np.allclose(cube[..., 0], np.arange(case.hours + 1)[:, None] * 3600, atol=.1, rtol=6e-5):
        raise ValueError("Unexpected basin times")
    if not np.allclose(cube[..., 1:3], basin["points"][None, ...], atol=.0001, rtol=6e-5):
        raise ValueError("Unexpected basin locations")
    if not np.allclose(cube[..., 6], basin["depth"][basin["wet"]][None, :], rtol=6e-5, atol=.01):
        raise ValueError("Basin bottom read incorrectly")
    if np.any(cube[..., 3] < 0) or np.any((cube[1:, :, 3:4] > 0) & (cube[1:, :, 4:6] <= 0)):
        raise ValueError("Invalid basin wave diagnostics")
    if np.any(cube[1:, :, 8] < 0) or np.any(cube[1:, :, 8] > 1):
        raise ValueError("Invalid basin breaker fraction")
    return cube


def stationary_boundary(path):
    """Retain the final frequency spectrum, directions and spreads at their supplied values."""
    lines = path.read_text().splitlines()
    dates = [i for i, s in enumerate(lines) if re.match(r"^\d{8}\.\d{6}\s", s)]
    if len(dates) != 1:
        raise ValueError("Expected exactly one offshore spectrum")
    tm = next(i for i, s in enumerate(lines) if s.split() and s.split()[0] == "TIME")
    loc = next(i for i, s in enumerate(lines) if s.split() and s.split()[0] in ("LOCATIONS", "LONLAT"))
    if int(lines[loc + 1].split()[0]) != 1:
        raise ValueError("Expected one offshore location")
    lines[loc] = "LOCATIONS"
    lines[loc + 2] = "0 0"
    lines.pop(dates[0])
    del lines[tm:tm + 2]
    return "\n".join(lines) + "\n"


def slope_files(boundary, slope, cells=400, gamma=.73, depth=50.):
    if not 0 < slope < 1 or cells < 20 or not 0 < gamma < 2:
        raise ValueError("Invalid slope, grid or breaking index")
    length = (depth - .25) / slope
    bottom = np.linspace(depth, .25, cells + 1)
    ratio = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    text = f"""PROJECT 'Lunar slope' '003'
MODE STATIONARY ONEDIMENSIONAL
COORDINATES CARTESIAN
SET GRAV {MOON_SURFACE_GRAVITY:.12g} RHO 1025
OUTPUT OPTIONS TABLE 16
CGRID REG 0 0 0 {length:.12g} 0 {cells} 0 &
 CIRCLE 36 {0.03*ratio:.12g} {3*ratio:.12g} 96
INPGRID BOTTOM 0 0 0 {cells} 0 {length/cells:.12g} 1
READINP BOTTOM 1 'depth.bot' 1 0 FREE
BOUNDSPEC SIDE WEST CONSTANT FILE 'boundary.spc' 1
GEN3 KOMEN DRAG WU
OFF QUAD
OFF WCAP
BREAKING CONSTANT 1.0 {gamma:.12g}
PROP BSBT
NUM ACCUR 0.005 0.005 0.005 99 STAT 100
CURVE 'line' 0 0 {cells} {length:.12g} 0
TABLE 'line' NOHEAD 'slope.tbl' XP HS RTP TM01 DEP DIR QB
COMPUTE
STOP
"""
    return {"INPUT": text, "depth.bot": array_text(bottom[None, :]), "boundary.spc": boundary}


def read_slope(path, files, cells):
    data = np.loadtxt(path, ndmin=2)
    depth = np.loadtxt(io.StringIO(files["depth.bot"]))
    if data.shape != (cells + 1, 7) or not np.isfinite(data).all():
        raise ValueError("Incomplete slope output")
    if not np.allclose(data[:, 4], depth, rtol=6e-5, atol=.0001) or np.any(data[:, 1:5] <= 0):
        raise ValueError("Invalid slope wave/depth output")
    if np.any(data[:, 6] < 0) or np.any(data[:, 6] > 1) or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError("Invalid slope positions/breaking")
    return data


def experiment(build_root, run_root):
    forcing_path = HERE / "results/forcing.json"
    forcing = json.loads(forcing_path.read_text())
    rho = forcing["smythii_marginis_equatorial_density_statistics"]["mean"]
    build = build_coupled(build_root, rho)
    executable = build_root / "coupled_air/swan.exe"
    winds, wind_input = wind_inputs()
    cases = [BasinCase("smythii_east"), BasinCase("smythii_turn", turn=True),
             BasinCase("smythii_east_dt30", step_s=30),
             BasinCase("smythii_east_frequency48", frequency_intervals=48),
             BasinCase("smythii_east_half_grid", stride=2, frequency_intervals=48)]
    basins, cubes, records, summaries = {}, {}, {}, []
    for case in cases:
        basin = load_basin(stride=case.stride)
        files = basin_files(case, basin, winds["p90"])
        directory = run_root / case.name
        records[case.name] = execute(executable, directory, files, ["basin.tbl", "offshore.spc"])
        cube = read_basin(directory / "basin.tbl", case, basin)
        basins[case.name], cubes[case.name] = basin, cube
        oi = basin["offshore_index"]
        summaries.append(dict(**asdict(case), wet_nodes=len(basin["points"]), offshore_lon_lat=basin["offshore"].tolist(),
                              offshore_final=cube[-1, oi].tolist(), final_maximum_hs_m=float(cube[-1, :, 3].max()),
                              final_node_median_hs_m=float(np.median(cube[-1, :, 3])),
                              maximum_breaking_fraction=float(cube[1:, :, 8].max())))
        print(case.name, "max Hs", summaries[-1]["final_maximum_hs_m"], "elapsed", records[case.name]["elapsed_wall_s"], flush=True)
    checks = {}
    for name in ("smythii_east_dt30", "smythii_east_frequency48", "smythii_east_half_grid"):
        reference = "smythii_east_frequency48" if name == "smythii_east_half_grid" else "smythii_east"
        base = cubes[reference]
        candidate = cubes[name]
        lookup = {tuple(p): i for i, p in enumerate(basins[name]["points"])}
        idx = [lookup[tuple(p)] for p in basins[reference]["points"]]
        matched = candidate[:, idx]
        valid = (base[1:, :, 3] >= .1) & (matched[1:, :, 3] >= .1)
        error = np.abs(base[1:, :, 3][valid] / matched[1:, :, 3][valid] - 1)
        final_valid = valid[-1]
        oi = basins[reference]["offshore_index"]
        checks[name] = dict(reference_case=reference, hs_relative_max=float(error.max()),
                            hs_final_relative_max=float(np.max(np.abs(base[-1, final_valid, 3] / matched[-1, final_valid, 3] - 1))),
                            hs_final_absolute_max_m=float(np.max(np.abs(base[-1, :, 3] - matched[-1, :, 3]))),
                            mean_period_final_relative_max=float(np.max(np.abs(base[-1, final_valid, 5] / matched[-1, final_valid, 5] - 1))),
                            offshore_final_hs_relative_change=float(matched[-1, oi, 3] / base[-1, oi, 3] - 1),
                            offshore_final_tm01_relative_change=float(matched[-1, oi, 5] / base[-1, oi, 5] - 1))
    # Feed the calculated SPEC1D frequency spectrum and direction/spread moments to
    # imposed smooth slopes to test their response to that boundary spectrum.
    boundary = stationary_boundary(run_root / "smythii_east_dt30/offshore.spc")
    slopes = []
    for denominator, cells, gamma in ((20, 400, .73), (50, 400, .73), (100, 400, .73), (50, 800, .73), (50, 400, .6), (50, 400, .9)):
        name = f"slope_1in{denominator}_n{cells}_g{round(100*gamma)}"
        files = slope_files(boundary, 1 / denominator, cells=cells, gamma=gamma)
        folder = run_root / name
        records[name] = execute(executable, folder, files, ["slope.tbl"])
        data = read_slope(folder / "slope.tbl", files, cells)
        onset = np.flatnonzero(data[:, 6] >= .01)
        peak = int(np.argmax(data[:, 1]))
        slopes.append(dict(name=name, slope=1 / denominator, cells=cells, breaker_index=gamma,
                           maximum_hs_m=float(data[peak, 1]), depth_at_maximum_hs_m=float(data[peak, 4]),
                           one_percent_breaking_depth_m=float(data[onset[0], 4]) if len(onset) else None,
                           columns=["x_m", "hs_m", "peak_period_s", "mean_period_s", "depth_m", "direction_deg", "breaking_fraction"],
                           rows=data.tolist()))
        print(name, "max Hs", slopes[-1]["maximum_hs_m"], flush=True)
    return dict(schema="terluna.climate.wave-basin-pilot/1",
                evidence="Conditional SWAN basin episodes over atlas bathymetry and idealised stationary slopes driven by one calculated offshore spectrum. Empirical lunar validation, surf climatology, individual breaking waves, run-up, currents and sediment evolution remain open.",
                reading_rule="Wind blows uniformly eastward at the corrected CM1 water p90 magnitude. The turn case interpolates vector components from east at hour 5 to north at hour 6, then holds north; both speed and direction change. Occurrence probabilities require wind histories. The resolution comparisons and pilot_checks.json audit describe geographic uncertainty. Slopes isolate shoaling/refraction and depth breaking, with whitecapping, quadruplets, bottom friction and triads omitted. SPEC1D retains frequency variance and directional moments; SWAN reconstructs the directional distributions from those moments.",
                units=dict(basin_columns=["time_s", "longitude_deg", "latitude_deg", "hs_m", "peak_period_s", "mean_period_s", "depth_m", "direction_deg", "breaking_fraction"]),
                producer=dict(domain="climate", files={str(p.relative_to(ROOT)): sha256(p) for p in (Path(__file__), HERE / "build.py")}),
                inputs=dict(atlas_sha256=sha256(ATLAS), atlas_metadata=basins[cases[0].name]["metadata"], wind=wind_input,
                            forcing_record_sha256=sha256(forcing_path), shared_constants_sha256=sha256(ROOT / "shared/constants.json")),
                air_density_kg_m3=rho, moon_radius_m=MOON_RADIUS, build=build, basin_cases=summaries,
                checks=checks, slopes=slopes, runs=records,
                basin_fields={name: dict(points=basins[name]["points"].tolist(), hourly_values=cube.tolist())
                              for name, cube in cubes.items() if name in ("smythii_east", "smythii_turn")})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, default=ROOT / "research/runs/waves/pilot")
    parser.add_argument("--output", type=Path, default=HERE / "results/pilot.json")
    args = parser.parse_args()
    result = experiment(args.build_root, args.run_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result["checks"], indent=2))


if __name__ == "__main__":
    main()
