"""Recover regional wind context from the saved GCM means at their model height."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024**2), b""):
            h.update(block)
    return h.hexdigest()


def mean_output_settings(text):
    """Distinguish PlaSim's accumulated winds from its instantaneous output."""
    values = {}
    for key in ("NLOWIO", "NSNAPSHOT", "NHCADENCE"):
        match = re.search(r"^\s*"+key+r"\s*=\s*(\d+)\b", text, re.MULTILINE | re.IGNORECASE)
        if match is None:
            raise ValueError(f"Missing GCM output setting {key}")
        values[key] = int(match.group(1))
    if values["NLOWIO"] != 1:
        raise ValueError("Regional mean exporter requires NLOWIO=1")
    return values


def export(run_dir, years, output):
    import netCDF4
    progress_path = run_dir/"progress.json"
    namelist_path = run_dir/"model/plasim_namelist"
    settings = mean_output_settings(namelist_path.read_text())
    config = json.loads(progress_path.read_text())["configuration"]
    m = config["model"]
    step_s = m["timestep_min"]*60
    steps_per_year = m["steps_per_lunar_day"]*m["lunar_days_per_year"]
    window_steps = m["steps_per_lunar_day"]//m["writes_per_lunar_day"]
    records, arrays = [], {k: [] for k in ("time_s", "eastward_m_s", "northward_m_s", "water")}
    for year in years:
        path = run_dir/"model"/f"MOST.{year:05d}.nc"
        with netCDF4.Dataset(path) as d:
            lat_all, lon_all = np.asarray(d["lat"][:]), np.asarray(d["lon"][:])
            rows = np.flatnonzero((lat_all >= -15) & (lat_all <= 30))
            cols = np.flatnonzero((lon_all >= 70) & (lon_all <= 110))
            lat, lon = lat_all[rows], lon_all[cols]
            sigma = float(d["lev"][-1])
            times = np.asarray(d["time"][:])
            if (len(times) != steps_per_year//window_steps or
                    not np.allclose(times, np.arange(window_steps-1, steps_per_year, window_steps)) or
                    d["time"].units != "timesteps"):
                raise ValueError("Unexpected GCM output intervals")
            if d["ua"].units != "m s-1" or d["va"].units != "m s-1":
                raise ValueError("Unexpected GCM wind units")
            for variable, key in (("ua", "eastward_m_s"), ("va", "northward_m_s")):
                values = d[variable][:, -1, rows, cols]
                if (values.shape != (len(times), len(lat), len(lon)) or
                        np.ma.getmaskarray(values).any() or not np.isfinite(values).all()):
                    raise ValueError("Incomplete GCM wind field")
                arrays[key].append(np.asarray(values))
            water = np.asarray(d["lsm"][:, rows, cols]) < .5
            arrays["water"].append(water)
            arrays["time_s"].append((year*steps_per_year+times)*step_s)
            if records and (sigma != records[0]["sigma"] or not np.array_equal(lat, prev_lat) or not np.array_equal(lon, prev_lon)):
                raise ValueError("GCM coordinates changed between years")
            records.append(dict(path=str(path.relative_to(ROOT)), sha256=digest(path), year=year, sigma=sigma,
                                near_surface_variable_names=[v for v in ("uas", "vas") if v in d.variables]))
            prev_lat, prev_lon = lat, lon
    data = {k: np.concatenate(v) for k, v in arrays.items()}
    metadata = dict(schema="terluna.climate.gcm-regional-wind-context/1",
                    evidence="Saved means from the settled A28_dim5_moon GCM at its lowest model level. These fields describe circulation across the Smythii–Marginis region.",
                    reading_rule=f"Each vector averages the preceding {window_steps} model steps ({window_steps*step_s/3600:g} Earth hours), on the native {m['resolution']} Gaussian grid, at sigma {sigma:.4f}. Heights follow terrain and atmospheric thickness. Direct SWAN forcing requires 10 m vectors with finer temporal sampling. Timestamp values preserve the model's saved step convention; years index the archived files.",
                    units=dict(time_s="model seconds", longitude_deg="degrees east", latitude_deg="degrees north", eastward_m_s="m/s", northward_m_s="m/s"),
                    producer=dict(domain="climate", files={str(p.relative_to(ROOT)): digest(p) for p in
                                  (Path(__file__), ROOT/"climate/gcm/exoplasim_run.py")}),
                    inputs=dict(progress_sha256=digest(progress_path), namelist_sha256=digest(namelist_path), files=records),
                    output_settings=settings,
                    run=run_dir.name, years=list(years), snapshots=len(data["time_s"]),
                    averaging_interval_s=window_steps*step_s, sigma=sigma, grid_shape=[len(lat),len(lon)],
                    swan_surface_forcing_ready=False,
                    requirements=["Instantaneous 10 m eastward and northward vectors", "Hourly or three-hourly sampling", "Atmospheric surface-layer diagnostic checked at lunar gravity"],
                    source_configuration=config,
                    snapshot_files_available=len(list((run_dir/"model/snapshots").glob("*"))))
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, metadata=json.dumps(metadata), latitude_deg=lat, longitude_deg=lon, **data)
    report = dict(metadata, product=dict(path=str(output.relative_to(ROOT)), sha256=digest(output)))
    speed = np.hypot(data["eastward_m_s"],data["northward_m_s"])
    report["water_vector_mean_speed_percentiles_m_s"] = dict(zip(("median","p90","p99"),np.quantile(speed[data["water"]],[.5,.9,.99]).tolist()))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir",type=Path,default=ROOT/"climate/gcm/runs/A28_dim5_moon")
    p.add_argument("--output",type=Path,default=ROOT/"research/runs/waves/coverage/gcm_region.npz")
    p.add_argument("--record",type=Path,default=ROOT/"climate/waves/results/gcm_coverage.json")
    args = p.parse_args()
    record = export(args.run_dir,range(20,30),args.output)
    args.record.write_text(json.dumps(record,indent=2)+"\n")
    print({k:record[k] for k in ("snapshots","grid_shape","averaging_interval_s","sigma","water_vector_mean_speed_percentiles_m_s")})


if __name__ == "__main__":
    main()
