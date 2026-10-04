"""The nearside sea's waves through the second lunar cycle, beside Smythii–Marginis.

    python -m climate.waves.nearside_analysis assemble    # hourly fields into one cube
    python -m climate.waves.nearside_analysis check --executable <swan.exe>   # 150 s window
    python -m climate.waves.nearside_analysis summarize   # results/nearside.json

`assemble` reads the segmented run of nearside.py, checks every segment's
recorded hashes, times, coordinates and depths, and keeps hourly Hs, periods,
direction, breaking fraction and energy transport at the wet nodes in
research/runs/waves/nearside/<run>/sea.npz. `summarize` measures the second
lunar cycle as cycle.py does for Smythii–Marginis: heights, time above
thresholds, persistence and direction at every node, the reference points'
episodes and spectra, swell and wind sea at the reference points, and the
energy that reaches each coast. `check` repeats the stormiest three days of the
second cycle at Smythii's 150 s timestep from the main run's hotfile.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from climate.waves import nearside
from climate.waves.coastal_checks import nesting_spectra
from climate.waves.cycle import CYCLE_HOURS, exceedance, window
from climate.waves.cycle_checks import read_spectral_series
from climate.waves.model import ROOT, sha256
from shared.constants import MOON_SURFACE_GRAVITY

HERE = Path(__file__).resolve().parent
RESULT = HERE / "results/nearside.json"
SMYTHII = HERE / "results/cycle.json"
FIELDS = ("hs", "tp", "tm01", "dir", "qb", "tx", "ty")
THRESHOLDS = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0)
WAVE_AGE = 1.2 * 28        # wind sea where 28 u* cos(angle) / c exceeds 1/1.2, the WAM partition


def run_root(name):
    return nearside.RUNS / name


def read_table(path, nodes):
    data = np.fromstring(path.read_text(), sep=" ")
    if data.size % (nodes * 11):
        raise ValueError(f"Incomplete table: {path}")
    return data.reshape(-1, nodes, 11)


def assemble(name="dt300", *, end_hour=nearside.END_HOUR, path=None):
    root = run_root(name)
    progress = json.loads((root / "progress.json").read_text())
    if progress["completed_through_hour"] < end_hour:
        raise ValueError("The run has not reached the requested hour")
    progress["segments"] = [s for s in progress["segments"] if s["end_hour"] <= end_hour]
    sea = nearside.load_sea()
    wet = sea["wet"].ravel()
    nodes = wet.size
    hours, cube, active = [], [], None
    for segment in progress["segments"]:
        directory = root / segment["directory"]
        expected = segment["run"]["output_sha256"]["sea.tbl"]
        if sha256(directory / "sea.tbl") != expected:
            raise ValueError(f"Changed output: {directory}")
        data = read_table(directory / "sea.tbl", nodes)
        t = data[:, 0, 0] / 3600
        first = segment["start_hour"] + (0 if segment["start_hour"] == 0 else 1)
        if not np.allclose(t, np.arange(first, segment["end_hour"] + 1), atol=1e-6):
            raise ValueError(f"Unexpected times in {directory}")
        w = data[:, wet]
        if not np.allclose(w[..., 1:3], sea["points"][None], atol=1e-4):
            raise ValueError("Unexpected node coordinates")
        # Nodes shallower than SWAN's minimum depth stay dry (Hs -9, depth -99) throughout.
        dry = np.all(w[..., 3] < 0, axis=0)
        if active is None:
            active = ~dry
        elif not np.array_equal(active, ~dry):
            raise ValueError("SWAN's dry nodes changed between segments")
        if not np.allclose(w[:, active, 6], sea["depth"][sea["wet"]][None, active], rtol=6e-5, atol=0.01):
            raise ValueError("Bottom read incorrectly")
        w[:, ~active, 3:] = 0.0
        hours.append(t)
        cube.append(w[..., [3, 4, 5, 7, 8, 9, 10]].astype(np.float32))
    hours = np.concatenate(hours)
    cube = np.concatenate(cube)
    if not np.array_equal(hours, np.arange(end_hour + 1)):
        raise ValueError("The hourly record has gaps")
    if np.any(cube[:, active, 0] < 0) or np.any((cube[1:, active, 0] > 0) & (cube[1:, active, 2] <= 0)):
        raise ValueError("Invalid wave diagnostics")
    path = path or root / "sea.npz"
    np.savez_compressed(path, hours=hours, lon=sea["points"][:, 0], lat=sea["points"][:, 1],
                        depth=sea["depth"][sea["wet"]], active=active, **{f: cube[..., i] for i, f in enumerate(FIELDS)},
                        metadata=np.array(json.dumps(dict(
                            schema="terluna.climate.nearside-waves-cube/1", run=name,
                            segments={s["directory"]: s["run"]["output_sha256"]["sea.tbl"] for s in progress["segments"]},
                            units=dict(hs="m", tp="s", tm01="s", dir="degrees counterclockwise from east, toward propagation",
                                       qb="fraction", tx="m3/s per m, eastward", ty="m3/s per m, northward")))))
    return path


def load_cube(name="dt300", path=None):
    with np.load(path or run_root(name) / "sea.npz", allow_pickle=False) as d:
        return {k: d[k] for k in d.files}


def time_above(t, values, threshold):
    """Hours at or above a threshold on piecewise-linear histories, for every column."""
    a, b = values[:-1], values[1:]
    dt = np.diff(t)[:, None]
    both = (a >= threshold) & (b >= threshold)
    cross = (a >= threshold) != (b >= threshold)
    with np.errstate(divide="ignore", invalid="ignore"):
        part = (np.maximum(a, b) - threshold) / np.abs(b - a)
    return np.sum(dt * np.where(both, 1.0, np.where(cross, part, 0.0)), axis=0)


def node_statistics(t, hs, direction):
    dt = np.diff(t)[:, None]
    span = t[-1] - t[0]
    mean = np.sum(dt * (hs[:-1] + hs[1:]) * .5, axis=0) / span
    square = np.sum(dt * (hs[:-1]**2 + hs[:-1] * hs[1:] + hs[1:]**2) / 3, axis=0) / span
    vector = hs**2 * np.exp(1j * np.radians(direction))
    integrated = np.sum(dt * (vector[:-1] + vector[1:]) * .5, axis=0)
    energy = np.sum(dt * (hs[:-1]**2 + hs[1:]**2) * .5, axis=0)
    return dict(max_hs=hs.max(axis=0), mean_hs=mean, rms_hs=np.sqrt(square),
                direction=np.degrees(np.angle(integrated)) % 360,
                resultant=np.divide(np.abs(integrated), energy, out=np.zeros_like(energy), where=energy > 0),
                hours_above={c: time_above(t, hs, c) for c in THRESHOLDS})


def coastal_nodes(sea):
    """Wet nodes beside dry cells, with the unit vector toward the dry side (east, north)."""
    wet = sea["wet"]
    index = -np.ones(wet.shape, int)
    index[wet] = np.arange(wet.sum())
    rows, cols = np.nonzero(wet)
    normal = np.zeros((len(rows), 2))
    for dr, dc, vec in ((0, 1, (1, 0)), (0, -1, (-1, 0)), (1, 0, (0, 1)), (-1, 0, (0, -1))):
        r, c = rows + dr, cols + dc
        dry = ~wet[np.clip(r, 0, wet.shape[0] - 1), np.clip(c, 0, wet.shape[1] - 1)]
        normal[dry] += vec
    length = np.hypot(*normal.T)
    coast = length > 0
    normal[coast] /= length[coast, None]
    return np.flatnonzero(coast), normal[coast]


def swell_partition(name, sea, start_hour, end_hour):
    """Swell and wind sea at the reference points from their three-hourly directional spectra."""
    root = run_root(name)
    progress = json.loads((root / "progress.json").read_text())
    nodes = sea["wet"].size
    refs = sea["references"]
    flat_index = [np.ravel_multi_index(np.argwhere(sea["wet"])[r["node_index"]], sea["wet"].shape) for r in refs]
    times, spectra, ustar, wind_dir = [], [], [], []
    for segment in progress["segments"]:
        directory = root / segment["directory"]
        record = nesting_spectra(directory / "refs2d.spc")
        hours = np.array([(d - nearside.ORIGIN).total_seconds() / 3600 for d in record["dates"]])
        wind = np.fromstring((directory / "wind.tbl").read_text(), sep=" ").reshape(-1, nodes, 7)
        wind_hours = wind[:, 0, 0] / 3600
        keep = (hours >= start_hour) & (hours <= end_hour)
        if not keep.any():
            continue
        if not np.allclose(wind_hours, hours):
            raise ValueError("Wind and spectral times differ")
        times.append(hours[keep])
        spectra.append(record["spectra"][keep])
        w = wind[keep][:, flat_index]
        ustar.append(w[..., 5])
        wind_dir.append(np.degrees(np.arctan2(w[..., 4], w[..., 3])))
        frequency, direction = record["frequency"], record["direction"]
    times, spectra = np.concatenate(times), np.concatenate(spectra)       # (time, ref, f, theta) m2/Hz/deg
    ustar, wind_dir = np.concatenate(ustar), np.concatenate(wind_dir)
    speed = MOON_SURFACE_GRAVITY / (2 * math.pi * frequency)              # deep-water phase speed
    angle = np.radians(direction[None, None, None, :] - wind_dir[:, :, None, None])
    windsea = WAVE_AGE * ustar[:, :, None, None] * np.cos(angle) > speed[None, None, :, None]
    dtheta, df = 360 / len(direction), np.gradient(frequency)
    m0 = lambda e: np.sum(e * dtheta * df[None, None, :, None], axis=(2, 3))
    m1 = lambda e: np.sum(e * dtheta * (df * frequency)[None, None, :, None], axis=(2, 3))
    total, swell = m0(spectra), m0(np.where(windsea, 0, spectra))
    swell_m1 = m1(np.where(windsea, 0, spectra))
    out = []
    for k, r in enumerate(refs):
        ok = total[:, k] > (0.1 / 4)**2
        share = swell[ok, k] / total[ok, k]
        period = np.divide(swell[ok, k], swell_m1[ok, k], out=np.zeros(ok.sum()), where=swell_m1[ok, k] > 0)
        dominant = share > 0.5
        out.append(dict(name=r["name"], samples=int(ok.sum()),
                        swell_energy_share_mean=float(share.mean()),
                        time_swell_dominates=float(dominant.mean()),
                        swell_hs_mean_m=float((4 * np.sqrt(swell[ok, k])).mean()),
                        swell_hs_max_m=float((4 * np.sqrt(swell[:, k])).max()),
                        swell_mean_period_median_s=float(np.median(period[swell[ok, k] > 1e-4])) if np.any(swell[ok, k] > 1e-4) else None))
    return dict(rule=f"Wind sea where {WAVE_AGE:g} u* cos(theta - wind) exceeds the deep-water phase speed g/(2 pi f); the rest is swell. Spectra with resolved Hs of at least 0.1 m.",
                interval_hours=[float(times[0]), float(times[-1])], points=out)


def reference_spectral_checks(name, start_hour, end_hour):
    root = run_root(name)
    progress = json.loads((root / "progress.json").read_text())
    times, density = [], []
    for segment in progress["segments"]:
        t, xy, frequency, d = read_spectral_series(root / segment["directory"] / "refs1d.spc")
        times.append(t)
        density.append(d)
    t, density = np.concatenate(times), np.concatenate(density)
    keep = (t >= start_hour) & (t <= end_hour)
    t, density = t[keep], density[keep]
    m0 = np.trapezoid(density, frequency, axis=-1)
    low = np.divide(np.trapezoid(density[..., :3], frequency[:3], axis=-1), m0, out=np.zeros_like(m0), where=m0 > 0)
    high = np.divide(np.trapezoid(density[..., -3:], frequency[-3:], axis=-1), m0, out=np.zeros_like(m0), where=m0 > 0)
    peak = np.argmax(density, axis=-1)
    edge = (peak == 0) | (peak == len(frequency) - 1)
    included = 4 * np.sqrt(m0) >= .1
    bad = included & ((low >= .001) | (high >= .01) | edge)
    return dict(included_spectra=int(included.sum()), failing_spectra=int(bad.sum()),
                failing_resolved_hs_max_m=float((4 * np.sqrt(m0))[bad].max()) if bad.any() else None,
                high_edge_fraction_max=float(high[included].max()),
                rule="Three-hourly spectra at the reference points with resolved Hs of at least 0.1 m: under 0.1% of variance in the lowest two intervals, under 1% in the highest two, peak inside the band.")


def summarize(name="dt300", *, start=CYCLE_HOURS, end=2 * CYCLE_HOURS, cube_path=None):
    cube_path = cube_path or run_root(name) / "sea.npz"
    cube = load_cube(name, cube_path)
    sea = nearside.load_sea()
    hours = cube["hours"].astype(float)
    t, hs = window(hours, cube["hs"].astype(float), start, end)
    _, tm01 = window(hours, cube["tm01"].astype(float), start, end)
    _, tp = window(hours, cube["tp"].astype(float), start, end)
    unwrapped = np.degrees(np.unwrap(np.radians(cube["dir"].astype(float)), axis=0))
    _, direction = window(hours, unwrapped, start, end)
    direction %= 360
    _, tx = window(hours, cube["tx"].astype(float), start, end)
    _, ty = window(hours, cube["ty"].astype(float), start, end)
    active = cube["active"]
    weights = np.cos(np.radians(cube["lat"])) * active
    weights /= weights.sum()
    stats = node_statistics(t, hs, direction)
    peak_time, peak_node = np.unravel_index(np.argmax(hs), hs.shape)
    # Energy reaching each coast: shoreward component of the transport, times rho g.
    water = json.loads(nearside.WATER.read_text())["water_density_kg_m3"]
    coast, normal = coastal_nodes(sea)
    keep = active[coast]
    coast, normal = coast[keep], normal[keep]
    power = water * MOON_SURFACE_GRAVITY * np.maximum(tx[:, coast] * normal[:, 0] + ty[:, coast] * normal[:, 1], 0)
    dt = np.diff(t)[:, None]
    mean_power = np.sum(dt * (power[:-1] + power[1:]) * .5, axis=0) / (t[-1] - t[0])
    from geography import nomenclature
    features = [f for f in nomenclature.features() if f["type"].split(",")[0] in
                ("Mare", "Oceanus", "Sinus", "Lacus", "Palus", "Promontorium", "Mons", "Montes")]
    def feature_near(lon, lat):
        v = np.array([math.cos(math.radians(lat)) * math.cos(math.radians(lon)),
                      math.cos(math.radians(lat)) * math.sin(math.radians(lon)), math.sin(math.radians(lat))])
        best = max(features, key=lambda f: math.cos(math.radians(f["lat"])) * math.cos(math.radians(f["lon"])) * v[0]
                   + math.cos(math.radians(f["lat"])) * math.sin(math.radians(f["lon"])) * v[1] + math.sin(math.radians(f["lat"])) * v[2])
        return best["name"]
    # The strongest coasts, each at least 300 km from a stronger one already listed.
    coast_lon, coast_lat = np.radians(cube["lon"][coast]), np.radians(cube["lat"][coast])
    unit = np.stack((np.cos(coast_lat) * np.cos(coast_lon), np.cos(coast_lat) * np.sin(coast_lon), np.sin(coast_lat)), axis=-1)
    ranked = []
    for k in np.argsort(mean_power)[::-1]:
        if all(math.acos(min(1.0, float(unit[k] @ unit[j]))) * nearside.MOON_RADIUS > 300e3 for j in ranked):
            ranked.append(k)
        if len(ranked) == 10:
            break
    references = []
    for r in sea["references"]:
        k = r["node_index"]
        references.append(dict(r, max_hs_m=float(hs[:, k].max()), mean_hs_m=float(stats["mean_hs"][k]),
                               tm01_range_s=[float(tm01[hs[:, k] >= .1, k].min()), float(tm01[hs[:, k] >= .1, k].max())],
                               direction_deg=float(stats["direction"][k]), direction_resultant=float(stats["resultant"][k]),
                               thresholds={f"{c:g}": exceedance(t, hs[:, k], c) for c in (0.5, 1.0, 2.0, 3.0)},
                               series=dict(time_hours=t[::3].round(3).tolist(), hs_m=hs[::3, k].round(3).tolist(),
                                           tm01_s=tm01[::3, k].round(2).tolist())))
    smythii = json.loads(SMYTHII.read_text())
    smythii_thresholds = {f"{x['hs_m']:g}": x["area_time_fraction"] for x in smythii["statistics"]["thresholds"]}
    rows, cols = np.nonzero(sea["wet"])
    return dict(
        schema="terluna.climate.nearside-waves/1",
        evidence=("SWAN 41.51 over the nearside sea of the 28% atlas, 1-degree nodes, driven by 60 days of GCM surface "
                  "stress with the Smythii cycle's physics and a 300 s timestep; the first lunar cycle spins the sea up "
                  "and the second is measured."),
        reading_rule=("Statistics use the second solar cycle, interpolated hourly heights and exact cycle endpoints; "
                      "area fractions weight nodes by the cosine of latitude. Coastal power is the time-mean shoreward "
                      "component of SWAN's energy transport at 1-degree nodes beside land, per metre of coast."),
        producer=dict(domain="climate", files={str(p.relative_to(ROOT)): sha256(p) for p in
                      (Path(__file__), HERE / "nearside.py", HERE / "cycle.py", ROOT / "shared/constants.json")}),
        inputs=dict(cube=str(cube_path.relative_to(ROOT)) if cube_path.is_relative_to(ROOT) else str(cube_path),
                    cube_sha256=sha256(cube_path),
                    progress_sha256=sha256(run_root(name) / "progress.json"), atlas_sha256=sea["source_sha256"]),
        run=json.loads((run_root(name) / "progress.json").read_text())["segments"][0]["run"]["identity"]["executable_sha256"],
        clock=dict(report_interval_hours=[start, end], reported_hours=end - start),
        sea=dict(wet_nodes=int(sea["wet"].sum()), swan_dry_nodes=int((~active).sum()), area_km2=5288444,
                 references=len(sea["references"])),
        statistics=dict(
            max_hs_m=float(hs.max()), peak_lon_lat=[float(cube["lon"][peak_node]), float(cube["lat"][peak_node])],
            peak_day_in_reporting_cycle=float((t[peak_time] - start) / 24),
            peak_tm01_s=float(tm01[peak_time, peak_node]), peak_tp_s=float(tp[peak_time, peak_node]),
            peak_depth_m=float(cube["depth"][peak_node]),
            area_time_mean_hs_m=float(np.dot(weights, stats["mean_hs"])),
            area_time_fraction={f"{c:g}": float(np.dot(weights, stats["hours_above"][c]) / (end - start)) for c in THRESHOLDS},
            largest_area_fraction={f"{c:g}": float(((hs >= c) * weights).sum(axis=1).max()) for c in THRESHOLDS},
            node_hs_percentiles_of_max={p: float(np.percentile(stats["max_hs"], int(p[1:]))) for p in ("p50", "p90", "p99")}),
        smythii_marginis=dict(max_hs_m=smythii["statistics"]["max_hs_m"],
                              area_time_mean_hs_m=smythii["statistics"]["area_time_mean_hs_m"],
                              area_time_fraction=smythii_thresholds),
        series=dict(time_hours=t[::3].round(3).tolist(), max_hs_m=hs.max(axis=1)[::3].round(3).tolist(),
                    area_mean_hs_m=(hs * weights).sum(axis=1)[::3].round(4).tolist()),
        references=references,
        reference_spectra=reference_spectral_checks(name, start, end),
        swell=swell_partition(name, sea, start, end),
        coast=dict(nodes=int(len(coast)), mean_power_W_m=dict(median=float(np.median(mean_power)),
                                                               p90=float(np.percentile(mean_power, 90)),
                                                               max=float(mean_power.max())),
                   strongest=[dict(lon_lat=[float(cube["lon"][coast[k]]), float(cube["lat"][coast[k]])],
                                   mean_power_W_m=float(mean_power[k]), near=feature_near(cube["lon"][coast[k]], cube["lat"][coast[k]]))
                              for k in ranked[:10]]),
        maps=dict(active=active.astype(int).tolist(), rows=rows.tolist(), cols=cols.tolist(), lon=cube["lon"].round(3).tolist(), lat=cube["lat"].round(3).tolist(),
                  max_hs_m=stats["max_hs"].round(3).tolist(), mean_hs_m=stats["mean_hs"].round(3).tolist(),
                  direction_deg=stats["direction"].round(1).tolist(), direction_resultant=stats["resultant"].round(3).tolist(),
                  hours_above_1m=stats["hours_above"][1.0].round(1).tolist(), hours_above_2m=stats["hours_above"][2.0].round(1).tolist(),
                  coastal_nodes=coast.tolist(), coastal_mean_power_W_m=mean_power.round(1).tolist()))


def timestep_check(executable, name="dt300", step_s=150):
    """Repeat the stormiest days of the second cycle at a finer timestep from the main run's hotfile."""
    sea = nearside.load_sea()
    history = nearside.stress_history(sea)
    wet = sea["wet"]
    weights = np.cos(np.radians(sea["lat"]))[:, None] * np.ones_like(wet, float)
    stress = np.hypot(history["tx"], history["ty"])[:, wet] @ weights[wet] / weights[wet].sum()
    hours = history["times"] / 3600
    second = (hours >= CYCLE_HOURS) & (hours <= 2 * CYCLE_HOURS)
    peak = float(hours[second][np.argmax(stress[second])])
    start = int(48 * math.floor((peak - 24) / 48))
    start = min(max(start, 48 * math.ceil(CYCLE_HOURS / 48)), 48 * ((nearside.END_HOUR - 96) // 48))
    hotstart = run_root(name) / f"s{start - 48:04d}_{start:04d}" / "end.hot"
    check = run_root(f"check_dt{step_s}_from_{start}")
    check.mkdir(parents=True, exist_ok=True)
    segments = []
    for a, b in ((start, start + 48), (start + 48, start + 96)):
        segments.append(nearside.run_segment(executable, sea, history, check / f"s{a:04d}_{b:04d}", a, b,
                                             step_s=step_s, frequency_intervals=48, hotstart=hotstart))
        hotstart = check / f"s{a:04d}_{b:04d}" / "end.hot"
    nodes = wet.size
    main = np.concatenate([read_table(run_root(name) / f"s{a:04d}_{a + 48:04d}" / "sea.tbl", nodes) for a in (start, start + 48)])
    fine = np.concatenate([read_table(check / f"s{a:04d}_{a + 48:04d}" / "sea.tbl", nodes) for a in (start, start + 48)])
    keep = main[:, 0, 0] / 3600 >= start + 24
    m, f = main[keep][:, wet.ravel()], fine[keep][:, wet.ravel()]
    out = dict(peak_basin_mean_stress_hour=peak, peak_basin_mean_stress_Pa=float(stress[second].max()),
               window_hours=[start + 24, start + 96], restart_hour=start, step_s=step_s, segments=segments)
    for label, floor in (("hs_at_least_0p1m", .1), ("hs_at_least_1m", 1.0)):
        both = (m[..., 3] >= floor) & (f[..., 3] >= floor)
        dh = np.abs(m[..., 3] - f[..., 3])[both] / f[..., 3][both]
        dt = np.abs(m[..., 5] - f[..., 5])[both] / f[..., 5][both]
        out[label] = dict(pairs=int(both.sum()), hs_max=float(dh.max()), hs_p95=float(np.percentile(dh, 95)),
                          tm01_max=float(dt.max()), tm01_p95=float(np.percentile(dt, 95)))
    t = main[keep][:, 0, 0] / 3600
    above_main = time_above(t, m[..., 3], 1.0)
    above_fine = time_above(t, f[..., 3], 1.0)
    out["hours_above_1m_change_max"] = float(np.abs(above_main - above_fine).max())
    (check / "comparison.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("action", choices=("assemble", "check", "summarize"))
    parser.add_argument("--name", default="dt300")
    parser.add_argument("--executable", type=Path)
    args = parser.parse_args()
    if args.action == "assemble":
        print(assemble(args.name))
    elif args.action == "check":
        print(json.dumps({k: v for k, v in timestep_check(args.executable, args.name).items() if k != "segments"}, indent=1))
    else:
        record = summarize(args.name)
        check = sorted(nearside.RUNS.glob("check_dt150_from_*/comparison.json"))
        if check:
            comparison = json.loads(check[-1].read_text())
            comparison.pop("segments")
            record["timestep_check"] = comparison
        RESULT.write_text(json.dumps(record, separators=(",", ":")) + "\n")
        print(RESULT)


if __name__ == "__main__":
    main()
