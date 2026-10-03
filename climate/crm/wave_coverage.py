"""Export geographic surface-wind vectors from the three corrected CM1 rings."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from climate.crm.wave_forcing import ROOT, FIELDS, digest, moist_density, read_series

CASES = ("ring_equator", "ring_70_45e", "ring_70_135e")


def geographic_vectors(u, v, track):
    """Rotate along-ring x and left-of-track y into local east and north."""
    lat, lon = np.asarray(track["lat"]), np.asarray(track["lon"])
    n = len(lat)
    s = 2*np.pi*(np.arange(n)+.5)/n
    tilt = np.radians(track["tilt_deg"])
    expected_lat = np.degrees(np.arcsin(np.sin(tilt)*np.sin(s)))
    expected_lon = (track["node_deg"] + np.degrees(np.arctan2(np.cos(tilt)*np.sin(s), np.cos(s)))) % 360
    if not np.allclose(lat, expected_lat, atol=6e-5) or not np.allclose((lon-expected_lon+180)%360-180, 0, atol=6e-5):
        raise ValueError("Ring positions differ from the declared great circle")
    cos_lat = np.cos(np.radians(expected_lat))
    east = np.cos(tilt)/cos_lat
    north = np.sin(tilt)*np.cos(s)/cos_lat
    if not np.allclose(east**2+north**2, 1, atol=1e-12):
        raise ValueError("Invalid geographic wind basis")
    return u*east-v*north, u*north+v*east


def export(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for name in CASES:
        case_dir = ROOT/"climate/crm/runs"/name
        summary_path = ROOT/"climate/results/crm"/f"ring_{name}.json"
        summary, geo, times, values, sources = read_series(case_dir, summary_path)
        track = geo["record"]["path"]
        east, north = geographic_vectors(values["u10"], values["v10"], track)
        if not np.allclose(np.hypot(east, north), values["s10"], atol=1e-4, rtol=1e-4):
            raise ValueError("Geographic wind speed differs from CM1 s10")
        density = moist_density(values["psfc"], values["t2"], values["q2"])
        lon, lat = np.asarray(track["lon"]), np.asarray(track["lat"])
        water = ~geo["land"]
        box = (lon >= 76) & (lon <= 102) & (lat >= -9) & (lat <= 22)
        metadata = dict(schema="terluna.climate.cm1-geographic-winds/1", case=name,
                        evidence=summary["evidence"],
                        reading_rule="Instantaneous three-hourly values at the original ring columns. East/north winds rotate CM1 along-track and left-of-track components using the declared great circle. Separate rings retain their own dynamical realisations. Surface density uses the CM1 moist-gas closure. Selected-byte hashes cover the six named surface fields in source_fields order.",
                        units=dict(time_s="model seconds", longitude_deg="degrees east", latitude_deg="degrees north",
                                   eastward_10m_m_s="m/s", northward_10m_m_s="m/s", density_kg_m3="kg/m3"),
                        producer=dict(domain="climate", files={str(p.relative_to(ROOT)): digest(p) for p in
                            (Path(__file__), ROOT/"climate/crm/wave_forcing.py", ROOT/"climate/crm/ring_analysis.py", ROOT/"shared/constants.json")}),
                        inputs={p.name: digest(p) for p in (summary_path, case_dir/"case.json", case_dir/"cm1out_s.ctl",
                                case_dir/"cm1out_metadata.ctl", case_dir/"cm1out_metadata.dat", case_dir/"terluna_surface.txt")},
                        source_fields=list(FIELDS), sources=sources, tilt_deg=track["tilt_deg"], node_deg=track["node_deg"],
                        snapshots=len(times), columns=len(lon), span_days=summary["span_days"],
                        basin_rectangle_columns=int(box.sum()), basin_rectangle_water_columns=int((box & water).sum()))
        output = output_dir/f"{name}.npz"
        np.savez_compressed(output, metadata=json.dumps(metadata), time_s=times, longitude_deg=lon, latitude_deg=lat,
                            water=water, eastward_10m_m_s=east.astype(np.float32), northward_10m_m_s=north.astype(np.float32),
                            density_kg_m3=density.astype(np.float32))
        records.append({k: v for k, v in metadata.items() if k not in ("sources", "source_fields")})
        records[-1]["product"] = dict(path=str(output.relative_to(ROOT)), sha256=digest(output))
        records[-1]["water_wind_percentiles_m_s"] = dict(zip(("median", "p90", "p99"),
                                                          np.quantile(np.hypot(east, north)[:, water], [.5,.9,.99]).tolist()))
    return dict(schema="terluna.climate.cm1-wind-coverage/1", cases=records,
                evidence="Three corrected CM1 ring realisations with geographic 10 m wind vectors.",
                reading_rule="The tilted rings expand latitude coverage; the equatorial ring supplies the columns crossing the Smythii–Marginis rectangle. Basin-wide interpolation across independent rings would require a separate atmospheric assumption.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, default=ROOT/"research/runs/waves/coverage/rings")
    p.add_argument("--record", type=Path, default=ROOT/"climate/waves/results/ring_coverage.json")
    args = p.parse_args()
    record = export(args.output_dir)
    args.record.write_text(json.dumps(record, indent=2)+"\n")
    print([(c["case"], c["snapshots"], c["basin_rectangle_water_columns"]) for c in record["cases"]])


if __name__ == "__main__":
    main()
