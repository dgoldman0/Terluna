"""Recover wind vectors and moist surface density from existing CM1 snapshots.

The equatorial ring supplies eastward u and northward v. The exporter reads
selected surface byte ranges and preserves the ring's sampling geometry.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from climate.crm.ring_analysis import case_geometry, grads_layout
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT, CM1_VAPOUR_GAS_CONSTANT

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "terluna.climate.cm1-wave-forcing/1"
FIELDS = ("u10", "v10", "s10", "psfc", "t2", "q2")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def moist_density(pressure_pa, temperature_k, mixing_ratio):
    p, t, r = np.broadcast_arrays(pressure_pa, temperature_k, mixing_ratio)
    if not all(np.isfinite(a).all() for a in (p, t, r)) or np.any(p <= 0) or np.any(t <= 0) or np.any(r < 0):
        raise ValueError("Invalid surface thermodynamics")
    return p * (1 + r) / (t * (CM1_DRY_AIR_GAS_CONSTANT + r * CM1_VAPOUR_GAS_CONSTANT))


def read_surface(path, nx, layout):
    offsets, cursor = {}, 0
    for key, levels in layout:
        if key in FIELDS:
            if levels != 1:
                raise ValueError("Expected a surface variable: " + key)
            offsets[key] = cursor
        cursor += levels * nx * 4
    if set(offsets) != set(FIELDS) or path.stat().st_size != cursor:
        raise ValueError("Missing fields or snapshot size inconsistent with control file")
    fields, h = {}, hashlib.sha256()
    with path.open("rb") as stream:
        for key in FIELDS:
            stream.seek(offsets[key])
            raw = stream.read(nx * 4)
            h.update(raw)
            fields[key] = np.frombuffer(raw, dtype="<f4").copy()
    if not all(np.isfinite(v).all() and np.all(v > -99999) for v in fields.values()):
        raise ValueError("Missing/nonfinite CM1 surface values")
    return fields, dict(file=path.name, size_bytes=cursor, selected_bytes_sha256=h.hexdigest(),
                       byte_offsets=offsets, bytes_per_field=nx * 4)


def export(case_dir, summary_path, output):
    summary = json.loads(summary_path.read_text())
    if summary.get("schema") != "terluna.climate.crm-ring/1":
        raise ValueError("Unsupported ring summary")
    geo = case_geometry(case_dir)
    record = geo["record"]
    track = record["path"]
    if not track or not np.allclose(track["lat"], 0) or track["node_deg"] != 0:
        raise ValueError("This exporter requires the equatorial ring with node zero")
    nx, layout = grads_layout(case_dir / "cm1out_s.ctl")
    if nx != record["grid"]["nx"]:
        raise ValueError("Ring control and case grids disagree")
    meta_nx, meta_layout = grads_layout(case_dir / "cm1out_metadata.ctl")
    if meta_nx != 1 or any(n != 1 for _, n in meta_layout):
        raise ValueError("Unsupported CM1 metadata layout")
    table = np.fromfile(case_dir / "cm1out_metadata.dat", dtype="<f4").reshape(-1, len(meta_layout))
    names = [n for n, _ in meta_layout]
    times, numbers = table[:, names.index("mtime")], table[:, names.index("nwrite")]
    start, end = np.asarray(summary["span_days"]) * 86400
    selected = (times >= start) & (times <= end)
    times, numbers = times[selected], numbers[selected]
    if len(times) != summary["snapshots"] or len(times) < 2 or not np.all(np.diff(times) == record["configuration"]["output_s"]):
        raise ValueError("Incomplete or irregular snapshot times")
    if not np.array_equal(numbers, numbers.astype(int)) or np.any(np.diff(numbers) != 1):
        raise ValueError("Invalid snapshot numbers")
    values, sources = {k: [] for k in FIELDS}, []
    for n in numbers.astype(int):
        fields, source = read_surface(case_dir / f"cm1out_t{n:06d}_s.dat", nx, layout)
        for k in FIELDS:
            values[k].append(fields[k])
        sources.append(source)
    values = {k: np.asarray(v) for k, v in values.items()}
    density = moist_density(values["psfc"], values["t2"], values["q2"])
    water = ~geo["land"]
    lon = np.asarray(track["lon"]) % 360
    local = water & (lon >= 75) & (lon <= 103)
    stats = lambda x: dict(minimum=float(x.min()), median=float(np.median(x)), mean=float(x.mean()),
                          p90=float(np.quantile(x, .9)), p99=float(np.quantile(x, .99)), maximum=float(x.max()))
    speed = np.hypot(values["u10"], values["v10"])
    if not np.allclose(speed, values["s10"], rtol=1e-4, atol=1e-4):
        raise ValueError("CM1 wind-vector magnitudes disagree with s10")
    metadata = dict(schema=SCHEMA, evidence=summary["evidence"],
                    reading_rule="Instantaneous vectors every three hours along one equatorial line, at the original sample locations and times. Density is total moist air at 2 m using CM1's ideal-gas closure, with condensate loading omitted. Selected-byte hashes cover the six listed surface fields in FIELDS order.",
                    units=dict(time_s="model seconds", longitude_deg="east degrees", latitude_deg="north degrees",
                               u10="eastward m/s", v10="northward m/s", s10="m/s", psfc="Pa", t2="K", q2="kg vapour / kg dry air", density_kg_m3="kg/m3"),
                    producer=dict(domain="climate", files={str(p.relative_to(ROOT)): digest(p) for p in
                                  (Path(__file__), ROOT / "climate/crm/ring_analysis.py", ROOT / "shared/constants.json")}),
                    inputs={p.name: digest(p) for p in (summary_path, case_dir / "case.json", case_dir / "cm1out_s.ctl",
                            case_dir / "cm1out_metadata.ctl", case_dir / "cm1out_metadata.dat", case_dir / "terluna_surface.txt")},
                    source_fields=list(FIELDS), sources=sources, span_days=summary["span_days"], snapshots=len(times),
                    water_density_statistics=stats(density[:, water]), water_wind_statistics=stats(speed[:, water]),
                    smythii_marginis_equatorial_water_columns=int(local.sum()),
                    smythii_marginis_equatorial_density_statistics=stats(density[:, local]))
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, metadata=json.dumps(metadata), time_s=times, longitude_deg=lon,
                        latitude_deg=np.asarray(track["lat"]), water=water,
                        density_kg_m3=density.astype(np.float32), **values)
    return dict(schema=SCHEMA, evidence=metadata["evidence"], reading_rule=metadata["reading_rule"],
                units=metadata["units"], producer=metadata["producer"], inputs=metadata["inputs"],
                product=dict(path=str(output.relative_to(ROOT)), sha256=digest(output)),
                **{k: metadata[k] for k in ("span_days", "snapshots", "water_density_statistics", "water_wind_statistics",
                                            "smythii_marginis_equatorial_water_columns", "smythii_marginis_equatorial_density_statistics")})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, default=ROOT / "climate/crm/runs/ring_equator")
    parser.add_argument("--summary", type=Path, default=ROOT / "climate/results/crm/ring_ring_equator.json")
    parser.add_argument("--output", type=Path, default=ROOT / "research/runs/waves/inputs/ring_equator.npz")
    parser.add_argument("--record", type=Path, default=ROOT / "climate/waves/results/forcing.json")
    args = parser.parse_args()
    report = export(args.case_dir, args.summary, args.output)
    args.record.parent.mkdir(parents=True, exist_ok=True)
    args.record.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k.endswith("statistics")}, indent=2))


if __name__ == "__main__":
    main()
