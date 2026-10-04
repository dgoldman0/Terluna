"""Waves of the nearside sea through two lunar cycles of recovered weather.

The nearside sea (atlas body 433) holds half of the Moon's standing water and
its Earth-facing coasts. The runner drives it with the physics and numerics of
the Smythii–Marginis cycle (cycle.md): SWAN 41.51 in the coupled-air build,
Komen growth and whitecapping, the GCM's surface stress inverted through
SWAN's Wu drag, the gravity-scaled AGROW knee, depth breaking with index 0.73,
36 directions and 48 frequency intervals from 0.00497 to 0.497 Hz, on the
atlas's 1-degree nodes in spherical coordinates with the lunar radius. The
record starts at the first atmospheric snapshot and runs to hour 1,419: the
first lunar cycle spins the sea up and the second is analysed.

The sea spans the 0-degree meridian, so its grid runs from 85.875 W to
73.125 E and the GCM's periodic longitudes are extended across the seam.
The run is split into segments, each ending in a SWAN hotfile that starts the
next, so a stop costs at most one segment. Completed segments are reused after
their input, executable and output hashes match; a segment interrupted part
way is moved aside and run again.

    python -m climate.waves.nearside audit --executable <coupled-air swan.exe>
    python -m climate.waves.nearside probe --executable ... --step 150 --hours 6
    python -m climate.waves.nearside run --executable ... --step 150
"""
from __future__ import annotations

import argparse
from datetime import timedelta
import json
import math
from pathlib import Path
import shutil
import time

import numpy as np

from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.pilot import array_text, execute
from climate.waves.weather import apply_weights, drag_coefficient, ocean_weights, stress_equivalent_wind
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ATLAS = ROOT / "geography/products/atlas_28pct_4ppd.npz"
ATMOSPHERE = ROOT / "research/runs/waves/cycle/atmosphere.npz"
RUNS = ROOT / "research/runs/waves/nearside"
WATER = ROOT / "shared/scenarios/waves.json"
BODY = 433
END_HOUR = 1419
SEGMENT_HOURS = 48
# Reference points: the deep-water node nearest each IAU feature centre.
FEATURES = ("Oceanus Procellarum", "Mare Imbrium", "Mare Frigoris", "Mare Serenitatis",
            "Mare Tranquillitatis", "Mare Fecunditatis", "Mare Nubium", "Mare Humorum")
REFERENCE_DEPTH_M = 300.0
OUTPUTS = ("sea.tbl", "wind.tbl", "refs1d.spc", "refs2d.spc", "end.hot")
ADDRESS_SPACE = 6 * 1024**3


def stamp(hour):
    return (ORIGIN + timedelta(hours=hour)).strftime("%Y%m%d.%H%M%S")


def circular_crop(occupied):
    """Columns of the narrowest circular window holding every occupied column, one margin each side."""
    occupied = np.asarray(occupied, bool)
    n = len(occupied)
    index = np.flatnonzero(occupied)
    if not len(index) or len(index) == n:
        raise ValueError("A sea must leave some longitudes dry")
    gaps = np.diff(np.r_[index, index[0] + n])
    k = int(np.argmax(gaps))
    first, last = index[(k + 1) % len(index)], index[k]
    width = (last - first) % n + 3
    if width > n:
        raise ValueError("The sea leaves no dry margin in longitude")
    return (first - 1 + np.arange(width)) % n


def great_circle(lat1, lon1, lat2, lon2):
    a, b = np.radians(lat1), np.radians(lat2)
    return np.arccos(np.clip(np.sin(a) * np.sin(b) + np.cos(a) * np.cos(b) * np.cos(np.radians(lon2 - lon1)), -1, 1))


def load_sea(path=ATLAS, body=BODY, stride=4, features=FEATURES, reference_depth=REFERENCE_DEPTH_M):
    """The sea's atlas nodes on a grid continuous across the 0-degree meridian."""
    with np.load(path, allow_pickle=False) as source:
        metadata = json.loads(source["metadata"].item())
        if metadata["schema"] != "terluna.geography.atlas-grid/1":
            raise ValueError("Unsupported atlas grid")
        labels = source["water_label"][::stride, ::stride]
        height = source["height_m"][::stride, ::stride]
        lon_all, lat_all = source["lon_deg"][::stride], source["lat_deg"][::stride]
    rows, cols = np.where(labels == body)
    if not len(rows):
        raise ValueError("The sea is absent from the atlas")
    columns = circular_crop(np.isin(np.arange(len(lon_all)), cols))
    r0, r1 = rows.min() - 1, rows.max() + 1
    if r0 < 0 or r1 >= len(lat_all):
        raise ValueError("The sea reaches the grid's polar rows")
    step = float(lon_all[1] - lon_all[0])
    first = float(lon_all[columns[0]])
    lon = (first - 360.0 if first > 180.0 else first) + step * np.arange(len(columns))
    lat = lat_all[r0:r1 + 1][::-1]
    sub = lambda a: a[r0:r1 + 1][:, columns][::-1]
    wet = sub(labels) == body
    physical_depth = metadata["sea_level_m"] - sub(height)
    # Other waters inside the crop are held dry with finite depths, as in the pilot.
    depth = np.where(wet, physical_depth, np.minimum(physical_depth, -1.0))
    if np.any(depth[wet] <= 0) or wet[(0, -1), :].any() or wet[:, (0, -1)].any():
        raise ValueError("Invalid sea depth or open crop boundary")
    xx, yy = np.meshgrid(lon, lat)
    points = np.column_stack((xx[wet], yy[wet]))
    from geography import nomenclature
    named = nomenclature.by_name(nomenclature.features())
    deep = depth[wet] >= reference_depth
    references = []
    for name in features:
        feature = named[name]
        distance = great_circle(feature["lat"], feature["lon"], points[:, 1], points[:, 0])
        distance[~deep] = np.inf
        index = int(np.argmin(distance))
        references.append(dict(name=name, feature_lon_lat=[feature["lon"], feature["lat"]],
                               node_index=index, lon_lat=points[index].tolist(),
                               depth_m=float(depth[wet][index]),
                               offset_km=float(distance[index] * MOON_RADIUS / 1000)))
    return dict(body=body, stride=stride, lon=lon, lat=lat, wet=wet, depth=depth, points=points,
                references=references, metadata=metadata, source_sha256=sha256(path))


def periodic_weights(lon, lat, water, target_lon, target_lat):
    """ocean_weights on a longitude grid extended one column across the seam."""
    lon_ext = np.r_[lon[-1] - 360.0, lon, lon[0] + 360.0]
    water_ext = np.concatenate((water[:, -1:], water, water[:, :1]), axis=1)
    x = np.mod(np.asarray(target_lon, float) - lon_ext[0], 360.0) + lon_ext[0]
    indices, weights, fallback = ocean_weights(lon_ext, lat, water_ext, x, target_lat)
    rows, cols = np.divmod(indices, len(lon_ext))
    return rows * len(lon) + (cols - 1) % len(lon), weights, fallback


def stress_history(sea, atmosphere=ATMOSPHERE):
    """GCM surface stress at every sea node and archive time, with its interpolation record."""
    with np.load(atmosphere, allow_pickle=False) as d:
        meta = json.loads(d["metadata"].item())
        if meta["schema"] != "terluna.climate.gcm-wave-snapshots/1":
            raise ValueError("Unsupported atmospheric product")
        times = d["time_s"]
        lat, lon = d["latitude_deg"][::-1], d["longitude_deg"]
        mask = np.all((d["land_fraction"][:, ::-1] < .5) & (d["sea_ice_fraction"][:, ::-1] < .01), axis=0)
        xx, yy = np.meshgrid(sea["lon"], sea["lat"])
        indices, weights, fallback = periodic_weights(lon, lat, mask, xx, yy)
        fields = {name: apply_weights(d[name][:, ::-1], indices, weights)
                  for name in ("stress_eastward_Pa", "stress_northward_Pa", "stress_closure_density_kg_m3")}
    if len(times) < 2 or times[0] != 0 or not np.allclose(np.diff(times), 10800):
        raise ValueError("A complete three-hourly record is required")
    wet = sea["wet"]
    lon2, lat2 = np.meshgrid(lon, lat)
    distance = great_circle(yy[..., None], xx[..., None], lat2.ravel()[indices], lon2.ravel()[indices]) * meta["radius_m"]
    used = (weights > 0) & wet[..., None]
    record = dict(atmosphere_sha256=sha256(atmosphere), atmosphere_metadata=meta,
                  coupling="rho_swan * Cd_Wu(Ueq) * |Ueq| * (ueq,veq) = interpolated GCM surface stress vector",
                  spatial_interpolation="Bilinear ocean-corner weights renormalized over cells ice-free throughout the record, on GCM longitudes extended periodically across the 0-degree meridian; nearest open-water cell fills an empty stencil",
                  time_interpolation="Linear surface stress at 15-minute knots, then SWAN linear vector interpolation between knots",
                  regional_open_water_source_cells=int(len(np.unique(indices[used]))),
                  nearest_ocean_fallback_wet_nodes=int(fallback[wet].sum()),
                  farthest_contributing_ocean_cell_km=float(distance[used].max() / 1000))
    return dict(times=times, tx=fields["stress_eastward_Pa"], ty=fields["stress_northward_Pa"],
                rho=fields["stress_closure_density_kg_m3"], record=record)


def window_winds(history, start_h, end_h, density):
    """Stress-equivalent SWAN winds at 15-minute knots over [start_h, end_h]."""
    times = history["times"]
    if start_h < 0 or end_h <= start_h or end_h * 3600 > times[-1]:
        raise ValueError("The window must lie inside the atmospheric record")
    knots = np.arange(start_h * 3600, end_h * 3600 + 1, 900.0)
    position = knots / 10800.0
    left = np.minimum(position.astype(int), len(times) - 2)
    alpha = (position - left)[:, None, None]
    tx = (1 - alpha) * history["tx"][left] + alpha * history["tx"][left + 1]
    ty = (1 - alpha) * history["ty"][left] + alpha * history["ty"][left + 1]
    u, v = stress_equivalent_wind(tx, ty, density)
    # SWAN interpolates the winds linearly between knots; measure the stress that implies.
    um, vm = .5 * (u[:-1] + u[1:]), .5 * (v[:-1] + v[1:])
    speed = np.hypot(um, vm)
    implied = np.hypot(density * drag_coefficient(speed) * speed * um, density * drag_coefficient(speed) * speed * vm)
    wanted = np.hypot(.5 * (tx[:-1] + tx[1:]), .5 * (ty[:-1] + ty[1:]))
    return knots, u, v, float(np.abs(implied - wanted).max())


def segment_input(sea, start_h, end_h, *, step_s, frequency_intervals, water_density, hotstart):
    lon, lat = sea["lon"], sea["lat"]
    x0, y0 = lon[0], lat[0]
    nx, ny, delta = len(lon) - 1, len(lat) - 1, lon[1] - lon[0]
    ratio = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    first_hourly = start_h + (1 if hotstart else 0)
    first_3h = start_h + (3 if hotstart else 0)
    init = "INIT HOTSTART SINGLE 'initial.hot' UNFORMATTED" if hotstart else "INIT ZERO"
    return f"""PROJECT 'Nearside sea' '{sea['body']}'
MODE NONSTAT TWODIMENSIONAL
COORDINATES SPHERICAL {MOON_RADIUS:.12g} CCM
SET GRAV {MOON_SURFACE_GRAVITY:.12g} RHO {water_density:.12g}
OUTPUT OPTIONS TABLE 16
CGRID REG {x0:.12g} {y0:.12g} 0 {lon[-1]-x0:.12g} {lat[-1]-y0:.12g} {nx} {ny} &
 CIRCLE 36 {0.03*ratio:.12g} {3*ratio:.12g} {frequency_intervals}
INPGRID BOTTOM {x0:.12g} {y0:.12g} 0 {nx} {ny} {delta:.12g} {delta:.12g}
READINP BOTTOM 1 'depth.bot' 3 0 FREE
INPGRID WIND {x0:.12g} {y0:.12g} 0 {nx} {ny} {delta:.12g} {delta:.12g} &
 NONSTAT {stamp(start_h)} 15 MIN {stamp(end_h)}
READINP WIND 1 'wind.dat' 3 0 0 0 FREE
GEN3 KOMEN DRAG WU AGROW
BREAKING CONSTANT 1.0 0.73
PROP BSBT
{init}
QUANTITY TSEC REF {stamp(0)}
TABLE 'COMPGRID' NOHEAD 'sea.tbl' TSEC XP YP HS RTP TM01 DEP DIR QB TRANSP &
 OUTPUT {stamp(first_hourly)} 1 HR
TABLE 'COMPGRID' NOHEAD 'wind.tbl' TSEC XP YP WIND UFRI CDRAG &
 OUTPUT {stamp(first_3h)} 3 HR
POINTS 'refs' FILE 'refs.xy'
SPEC 'refs' SPEC1D ABS 'refs1d.spc' OUTPUT {stamp(first_3h)} 3 HR
SPEC 'refs' SPEC2D ABS 'refs2d.spc' OUTPUT {stamp(first_3h)} 3 HR
COMPUTE NONSTAT {stamp(start_h)} {step_s} SEC {stamp(end_h)}
HOTFILE 'end.hot' UNFORMATTED
STOP
"""


def segment_bounds(start_h=0, end_h=END_HOUR, length=SEGMENT_HOURS):
    edges = list(range(start_h, end_h, length)) + [end_h]
    return list(zip(edges[:-1], edges[1:]))


def build_record(executable):
    build = json.loads((executable.parent.parent / "coupled_air.json").read_text())
    if sha256(executable) != build["build"]["executable_sha256"]:
        raise ValueError("SWAN executable differs from its density/build record")
    if not build["build"].get("openmp", False):
        raise ValueError("The nearside sea needs the recorded OpenMP build")
    return build


def run_segment(executable, sea, history, directory, start_h, end_h, *, step_s, frequency_intervals,
                hotstart=None, threads=8, timeout_s=7200, address_space=ADDRESS_SPACE):
    build = build_record(executable)
    water_density = json.loads(WATER.read_text())["water_density_kg_m3"]
    knots, u, v, midpoint_error = window_winds(history, start_h, end_h, build["air_density_kg_m3"])
    files = {"INPUT": segment_input(sea, start_h, end_h, step_s=step_s, frequency_intervals=frequency_intervals,
                                    water_density=water_density, hotstart=hotstart is not None),
             "depth.bot": array_text(sea["depth"]),
             "refs.xy": array_text(np.array([r["lon_lat"] for r in sea["references"]])),
             "wind.dat": "".join(array_text(field) for pair in zip(u, v) for field in pair)}
    if hotstart is not None:
        files["initial.hot"] = hotstart
    if directory.exists() and not (directory / "run.json").exists():
        # An interrupted segment: keep its files for inspection and run it again.
        aside = directory.parent / "interrupted" / f"{directory.name}-{time.strftime('%Y%m%dT%H%M%S')}"
        aside.parent.mkdir(parents=True, exist_ok=True)
        directory.rename(aside)
    record = execute(executable, directory, files, list(OUTPUTS), timeout_s=timeout_s,
                     threads=threads, address_space_bytes=address_space)
    return dict(start_hour=start_h, end_hour=end_h, directory=directory.name,
                wind_midpoint_stress_error_max_Pa=midpoint_error, run=record)


def run(executable, name, *, step_s=150, frequency_intervals=48, start_h=0, end_h=END_HOUR,
        segment_hours=SEGMENT_HOURS, threads=8, run_root=RUNS):
    """Run (or resume) every segment in order, writing a progress record after each."""
    sea = load_sea()
    history = stress_history(sea)
    root = run_root / name
    root.mkdir(parents=True, exist_ok=True)
    progress = root / "progress.json"
    segments, hotstart = [], None
    for a, b in segment_bounds(start_h, end_h, segment_hours):
        directory = root / f"s{a:04d}_{b:04d}"
        began = time.monotonic()
        segment = run_segment(executable, sea, history, directory, a, b, step_s=step_s,
                              frequency_intervals=frequency_intervals, hotstart=hotstart, threads=threads)
        segment["session_wall_s"] = time.monotonic() - began
        segments.append(segment)
        hotstart = directory / "end.hot"
        progress.write_text(json.dumps(dict(
            schema="terluna.climate.nearside-run-progress/1", name=name, step_s=step_s,
            frequency_intervals=frequency_intervals, start_hour=start_h, end_hour=end_h,
            completed_through_hour=b, segments=segments, forcing=history["record"],
            sea=dict(body=sea["body"], stride=sea["stride"], wet_nodes=int(sea["wet"].sum()),
                     grid=[len(sea["lon"]), len(sea["lat"])], references=sea["references"],
                     atlas_sha256=sea["source_sha256"]),
            producer_sha256=sha256(Path(__file__))), indent=1) + "\n")
    return progress


def audit(executable):
    """Print the settings a run would apply, from the actual inputs, before any long run."""
    sea = load_sea()
    history = stress_history(sea)
    build = build_record(executable)
    wet = sea["wet"]
    knots, u, v, midpoint_error = window_winds(history, 0, 48, build["air_density_kg_m3"])
    speed = np.hypot(u, v)[:, wet]
    text = segment_input(sea, 0, 48, step_s=150, frequency_intervals=48,
                         water_density=json.loads(WATER.read_text())["water_density_kg_m3"], hotstart=False)
    # Stress steps between neighbouring columns: across the 0-degree meridian against everywhere else.
    seam = int(np.argmin(np.abs(sea["lon"] - 0.125)))
    stress = np.hypot(history["tx"], history["ty"])
    both = wet[:, 1:] & wet[:, :-1]
    steps = np.abs(np.diff(stress, axis=2))[:, both]
    seam_rows = both[:, seam - 1]
    seam_steps = np.abs(stress[:, seam_rows, seam] - stress[:, seam_rows, seam - 1])
    return dict(
        grid=dict(lon_first=float(sea["lon"][0]), lon_last=float(sea["lon"][-1]), lat_first=float(sea["lat"][0]),
                  lat_last=float(sea["lat"][-1]), shape=list(wet.shape), wet_nodes=int(wet.sum()),
                  depth_wet_min_m=float(sea["depth"][wet].min()), depth_wet_median_m=float(np.median(sea["depth"][wet])),
                  depth_wet_max_m=float(sea["depth"][wet].max()), nodes_shallower_than_10m=int((sea["depth"][wet] < 10).sum())),
        references=sea["references"],
        forcing=dict(history["record"], first_48h_speed_percentiles_m_s=np.percentile(speed, [50, 90, 99, 100]).round(3).tolist(),
                     midpoint_stress_error_max_Pa=midpoint_error,
                     seam_columns_lon=[float(sea["lon"][seam - 1]), float(sea["lon"][seam])],
                     seam_stress_step_p99_Pa=float(np.percentile(seam_steps, 99)),
                     all_column_stress_step_p99_Pa=float(np.percentile(steps, 99))),
        build=dict(air_density_kg_m3=build["air_density_kg_m3"], executable_sha256=build["build"]["executable_sha256"],
                   openmp=build["build"]["openmp"], native_optimization=build["build"]["native_optimization"]),
        input_text=text)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("action", choices=("audit", "probe", "run"))
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--name", default=None)
    parser.add_argument("--step", type=int, default=150)
    parser.add_argument("--hours", type=int, default=6, help="Probe length")
    parser.add_argument("--segment", type=int, default=SEGMENT_HOURS)
    parser.add_argument("--threads", type=int, choices=(1, 2, 4, 8), default=8)
    parser.add_argument("--end", type=int, default=END_HOUR)
    args = parser.parse_args()
    if args.action == "audit":
        print(json.dumps(audit(args.executable), indent=1))
    elif args.action == "probe":
        name = args.name or f"probe_dt{args.step}_h{args.hours}_t{args.threads}"
        print(run(args.executable, name, step_s=args.step, end_h=args.hours, segment_hours=args.segment,
                  threads=args.threads))
    else:
        name = args.name or f"dt{args.step}"
        print(run(args.executable, name, step_s=args.step, end_h=args.end, segment_hours=args.segment,
                  threads=args.threads))


if __name__ == "__main__":
    main()
