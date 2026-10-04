"""A month at four of the nearside sea's shores: its waves and its tide on one clock.

    python -m climate.waves.shore_month                  # results/shore_month.json
    <python with netCDF4> -m climate.waves.shore_month check-sun

The waves are the nearside run's second lunar cycle (nearside.md) at four
coastal nodes: the most exposed coast, in western Oceanus Procellarum; southern
Mare Imbrium at the eastern end of Montes Carpatus; the largest tide among the
most exposed tenth of the coast, in southern Mare Nubium; and a coast of Mare
Nectaris near the median exposure, also with a large tide. The tide comes from the geography
product (geography/tides.py) at the nearest quarter-degree node of the same sea.

The GCM's month has no calendar date. Its Sun turns evenly once per synodic
month (the synodic clock of climate/gcm/exoplasim_run.py), so every hour of the
wave record has a subsolar longitude. In every month of 2026-2045 the real Sun
stands once over the longitude where the GCM's Sun stands at the middle of the
wave month. Of these months, the one whose range stays closest to typical at
all four shores gives the tide, matched at that moment and run on the wave
month's clock. No month is typical at every shore: across the matched months the
ranges at southern Imbrium and southern Nubium correlate at -0.67.

`check-sun` compares the clock with the sunlight of the GCM's own last year
before the wave run (czen in its output, read with netCDF4) and writes the
comparison beside the nearside run; the build embeds it.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timedelta
import hashlib
import json
import math
from pathlib import Path
import struct

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / "results/shore_month.json"
NEARSIDE = HERE / "results/nearside.json"
ATMOSPHERE = ROOT / "research/runs/waves/cycle/atmosphere.npz"
SNAPSHOTS = ROOT / "research/runs/waves/cycle/snapshots"
SUN_CHECK = ROOT / "research/runs/waves/nearside/sun_clock_check.json"
TIDES = ROOT / "geography/products/tides_28pct_4ppd.npz"
BODY = 433
# Coastal nodes of the nearside run (lon, lat), with why each is shown.
SHORES = (
    ("Western Oceanus Procellarum, by Russell", "W Procellarum", -73.875, 28.875, "the most exposed coast"),
    ("Southern Mare Imbrium, at the eastern end of Montes Carpatus", "S Imbrium", -17.875, 15.875,
     "an exposed shore of the central sea"),
    ("Southern Mare Nubium, by Pitatus", "S Nubium", -11.875, -28.125, "the largest tide on the most exposed tenth of the coast"),
    ("Mare Nectaris, by Fracastorius", "Nectaris", 33.125, -22.125, "a coast near the median exposure, with a large tide"),
)


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def wrap(degrees):
    return (np.asarray(degrees) + 180.0) % 360.0 - 180.0


# ------------------------------------------------------------------ the GCM's Sun

def plasim_value(path, name):
    """The raw bytes of a named variable in a PlaSim restart or status file (Fortran sequential records)."""
    data = Path(path).read_bytes()
    pos = 0
    while pos + 4 <= len(data):
        size = struct.unpack_from("<i", data, pos)[0]
        record = data[pos + 4:pos + 4 + size]
        pos += size + 8
        if size == 16 and record.decode("ascii", "replace").strip() == name:
            size = struct.unpack_from("<i", data, pos)[0]
            return data[pos + 4:pos + 4 + size]
    raise KeyError(f"{name} not in {path}")


def gcm_clock(snapshots=SNAPSHOTS, atmosphere=ATMOSPHERE):
    """The step count and solar rate behind the wave archive, traced through the run's recorded hashes."""
    with np.load(atmosphere, allow_pickle=False) as archive:
        meta = json.loads(archive["metadata"].item())
    run = json.loads((snapshots / "run.json").read_text())
    if sha256(snapshots / "run.json") != meta["run_record_sha256"]:
        raise ValueError("The GCM run record differs from the one the wave archive names")
    if sha256(snapshots / "plasim_status") != run["output_sha256"]["plasim_status"]:
        raise ValueError("The GCM's final state differs from its run record")
    config = meta["source_configuration"]
    if config["model"]["sun_clock"] != "synodic":
        raise ValueError("The alignment needs the GCM's synodic Sun clock")
    end_step = struct.unpack("<i", plasim_value(snapshots / "plasim_status", "nstep"))[0]
    restart_step = end_step - run["output_settings"]["N_RUN_STEPS"]
    planet = config["planet"]
    return dict(restart_step=restart_step, first_snapshot_step=restart_step + meta["first_snapshot_restart_step"],
                step_s=run["timestep_s"], snapshot_interval_s=meta["sample_interval_s"],
                synodic_days=1 / (1 / planet["rotationperiod"] - 1 / planet["year"]),
                sidereal_day_s=planet["rotationperiod"] * 86400.0, sidereal_year_s=planet["year"] * 86400.0)


SUN_RULE = ("radmod.f90's synodic clock sets the hour angle at grid column jlon (0 at 0 E) to "
            "2 pi mod(nstep dt (1/sidereal day - 1/sidereal year), 1) + 2 pi jlon/nlon - pi, so the Sun stands "
            "over 180 - 360 mod(nstep dt (1/sidereal day - 1/sidereal year), 1) degrees east.")


def subsolar_at_step(steps, clock):
    """Longitude (degrees east) under the GCM's Sun at absolute model steps (SUN_RULE)."""
    rate = 1 / clock["sidereal_day_s"] - 1 / clock["sidereal_year_s"]
    turns = np.mod(np.asarray(steps, float) * clock["step_s"] * rate, 1.0)
    return np.mod(180.0 - 360.0 * turns, 360.0)


def subsolar_at_hour(hours, clock):
    """The same at hours since the first atmospheric snapshot, the wave runs' clock."""
    return subsolar_at_step(clock["first_snapshot_step"] + np.asarray(hours, float) * 3600 / clock["step_s"], clock)


def check_sun(snapshots=SNAPSHOTS, atmosphere=ATMOSPHERE, path=SUN_CHECK):
    """The clock against the sunlight of the GCM year that ends at the wave run's restart."""
    import netCDF4
    clock = gcm_clock(snapshots, atmosphere)
    with np.load(atmosphere, allow_pickle=False) as archive:
        model = json.loads(archive["metadata"].item())["source_configuration"]["model"]
    run = json.loads((snapshots / "run.json").read_text())
    restart = Path(run["source_restart"])
    output = restart.with_name(f"MOST.{restart.suffix[1:]}.nc")
    if not model["whole_lunar_days"]:
        raise ValueError("The year's length in steps is known only for whole lunar days")
    year_steps = model["steps_per_lunar_day"] * model["lunar_days_per_year"]
    write_steps = model["steps_per_lunar_day"] // model["writes_per_lunar_day"]
    with netCDF4.Dataset(output) as data:
        czen = np.asarray(data["czen"][:], float)
        lat, lon = np.asarray(data["lat"][:], float), np.asarray(data["lon"][:], float)
        last = np.asarray(data["time"][:], float)            # the last step of each mean, counted within the year
    rows = np.abs(lat) < 20.0
    band = czen[:, rows, :].mean(axis=1)
    centroid = np.degrees(np.angle(band @ np.exp(1j * np.radians(lon)))) % 360.0
    middle = clock["restart_step"] - year_steps + last - (write_steps - 1) / 2
    difference = wrap(centroid - subsolar_at_step(middle, clock))
    record = dict(schema="terluna.climate.gcm-sun-clock-check/1", output=str(output), output_sha256=sha256(output),
                  means=int(len(last)), restart_step=clock["restart_step"],
                  difference_deg=dict(mean=float(difference.mean()), max_abs=float(np.abs(difference).max())),
                  rule=("Each output mean's sunlight-weighted longitude of czen between 20 S and 20 N, against the "
                        "clock's subsolar longitude at the mean's middle step."))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(record, indent=1) + "\n")
    return record


# ------------------------------------------------------------------ the month

def calendar(jd):
    return (datetime(2000, 1, 1, 12) + timedelta(days=float(jd) - 2451545.0)).strftime("%Y-%m-%d %H:%M")


def matching_times(jd, sun_lon, target):
    """Times at which the real subsolar longitude, which falls about 12 degrees a day, passes `target`."""
    d = wrap(np.asarray(sun_lon) - target)
    k = np.flatnonzero((d[:-1] > 0) & (d[1:] <= 0) & (d[:-1] - d[1:] < 90))
    return jd[k] + (jd[k + 1] - jd[k]) * d[k] / (d[k] - d[k + 1])


def build():
    from climate.waves import nearside, nearside_analysis as analysis
    from climate.waves.cycle import CYCLE_HOURS
    from geography import lunar_ephemeris as ephemeris, tides
    from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY
    product = json.loads(NEARSIDE.read_text())
    cube_path = ROOT / product["inputs"]["cube"]
    if sha256(cube_path) != product["inputs"]["cube_sha256"]:
        raise ValueError("The wave cube differs from the one the nearside product measured")
    cube = analysis.load_cube(path=cube_path)
    hours = cube["hours"].astype(float)
    keep = (hours >= CYCLE_HOURS) & (hours <= 2 * CYCLE_HOURS)
    t = hours[keep]
    clock = gcm_clock()
    sun = subsolar_at_hour(t, clock)
    sea = nearside.load_sea()
    coast, normal = analysis.coastal_nodes(sea)
    water = json.loads(nearside.WATER.read_text())["water_density_kg_m3"]
    with np.load(TIDES, allow_pickle=False) as z:
        if json.loads(z["metadata"].item())["schema"] != "terluna.geography.tides-grid/1":
            raise ValueError("Unsupported tide product")
        own = z["labels"] == BODY
        tide_lat, tide_lon = z["lat_deg"][z["rows"][own]], z["lon_deg"][z["cols"][own]]
        basis, typical = z["basis"][:, own], z["monthly_range_median_m"][own]
    unit = lambda lo, la: np.stack((np.cos(np.radians(la)) * np.cos(np.radians(lo)),
                                    np.cos(np.radians(la)) * np.sin(np.radians(lo)), np.sin(np.radians(la))), axis=-1)
    tide_points = unit(tide_lon, tide_lat)
    shores, pixels = [], []
    for name, short, lon, lat, reason in SHORES:
        node = int(np.flatnonzero((np.abs(cube["lon"] - lon) < 1e-3) & (np.abs(cube["lat"] - lat) < 1e-3))[0])
        k = int(np.flatnonzero(coast == node)[0])
        pixel = int(np.argmax(tide_points @ unit(lon, lat)))
        pixels.append(pixel)
        n = normal[k]
        flux = cube["tx"][keep, node] * n[0] + cube["ty"][keep, node] * n[1]
        offset = math.acos(min(1.0, float(tide_points[pixel] @ unit(lon, lat)))) * MOON_RADIUS / 1000
        shores.append(dict(
            name=name, short=short, reason=reason, lon_lat=[lon, lat], depth_m=float(cube["depth"][node]),
            land_toward_deg=float(np.degrees(np.arctan2(n[1], n[0])) % 360),
            mean_power_W_m=float(np.array(product["maps"]["coastal_mean_power_W_m"])[
                np.flatnonzero(np.array(product["maps"]["coastal_nodes"]) == node)[0]]),
            tide_node_lon_lat=[float(tide_lon[pixel]), float(tide_lat[pixel])], tide_node_offset_km=offset,
            typical_monthly_tide_range_m=float(typical[pixel]),
            hs_m=cube["hs"][keep, node].round(3).tolist(), tm01_s=cube["tm01"][keep, node].round(2).tolist(),
            tp_s=cube["tp"][keep, node].round(2).tolist(), dir_deg=cube["dir"][keep, node].round(1).tolist(),
            power_W_m=(water * MOON_SURFACE_GRAVITY * np.maximum(flux, 0)).round(1).tolist()))
    # The tide's mean over the product's record, which the atlas's static level already holds.
    days = np.arange(0.0, tides.YEARS * 365.25, tides.STEP_DAYS)
    record_jd = tides.JD_START + days
    mean = tides.forcing_coefficients(record_jd).mean(axis=0)
    real_sun = ephemeris.earth_and_sun(record_jd)["sun_lon_lat_deg"][0]
    middle_hour = 1.5 * CYCLE_HOURS
    target = float(subsolar_at_hour(middle_hour, clock))
    offsets = (t - middle_hour) / 24
    candidates = matching_times(record_jd, real_sun, target)
    candidates = candidates[(candidates + offsets[0] >= record_jd[0]) & (candidates + offsets[-1] <= record_jd[-1])]
    when = candidates[:, None] + offsets[None, :]
    levels = ((tides.forcing_coefficients(when.ravel()) - mean) @ basis[:, pixels]).reshape(len(candidates), len(t), -1)
    ranges = np.ptp(levels, axis=1)                                  # (candidate, shore)
    typical_ranges = np.median(ranges, axis=0)
    departure = np.max(np.abs(np.log(ranges / typical_ranges)), axis=1)
    pick = int(np.argmin(departure))
    chosen = when[pick]
    mismatch = wrap(ephemeris.earth_and_sun(chosen)["sun_lon_lat_deg"][0] - sun)
    for k, (shore, level) in enumerate(zip(shores, levels[pick].T)):
        shore["tide_m"] = level.round(3).tolist()
        shore["tide_range_this_month_m"] = float(ranges[pick, k])
        shore["tide_range_matched_months_m"] = dict(p5=float(np.percentile(ranges[:, k], 5)), median=float(typical_ranges[k]),
                                                   p95=float(np.percentile(ranges[:, k], 95)))
    producer_files = ("climate/waves/shore_month.py", "climate/waves/nearside_analysis.py", "geography/tides.py",
                      "geography/lunar_ephemeris.py", "shared/constants.json")
    return dict(
        schema="terluna.climate.shore-month/1",
        evidence=("Hourly sea states at four coastal 1-degree nodes from the nearside run's second lunar cycle, and the "
                  "equilibrium tide of the geography product at the nearest quarter-degree node for a real month of "
                  "2026-2045 matched to the GCM's Sun."),
        reading_rule=("time_hours count from the first atmospheric snapshot, as in the wave runs. Heights are "
                      "significant wave heights; tm01_s and tp_s the mean and peak periods; dir_deg the mean "
                      "direction of travel, counterclockwise from east; power_W_m the shoreward component of SWAN's "
                      "energy transport per metre of coast. tide_m is the water level relative to its 2026-2045 mean. "
                      "Air pressure moves the level by under 0.2 m and wind setup by about 1 cm (review.md); neither "
                      "is included. subsolar_longitude_deg is the GCM's Sun; a place is in daylight where its "
                      "longitude lies within 90 degrees of it."),
        producer=dict(domain="climate", files={f: sha256(ROOT / f) for f in producer_files}),
        inputs=dict(nearside_sha256=sha256(NEARSIDE), cube_sha256=product["inputs"]["cube_sha256"],
                    tides_sha256=sha256(TIDES), gcm_run_record_sha256=sha256(SNAPSHOTS / "run.json")),
        gcm_sun=dict(clock, rule=SUN_RULE, check=json.loads(SUN_CHECK.read_text()) if SUN_CHECK.exists() else None),
        tide_month=dict(
            rule=("Among the times in 2026-2045 when the real Sun stands over the GCM's subsolar longitude at the "
                  "middle of the wave month, the one whose tidal range over the month departs least from typical at "
                  "its most atypical shore (each range against that shore's median over the matched months); the "
                  "tide then runs on the wave month's clock from that moment."),
            candidates=int(len(candidates)), start_tt=calendar(chosen[0]), end_tt=calendar(chosen[-1]),
            middle_jd_tt=float(candidates[pick]),
            largest_departure_from_typical_range=float(np.exp(departure[pick]) - 1),
            range_correlation_between_shores=np.corrcoef(ranges.T).round(3).tolist(),
            sun_mismatch_max_deg=float(np.abs(mismatch).max()),
            sun_mismatch_max_hours=float(np.abs(mismatch).max() / 360 * clock["synodic_days"] * 24)),
        time_hours=t.round(3).tolist(), subsolar_longitude_deg=sun.round(2).tolist(), shores=shores)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("action", nargs="?", default="build", choices=("build", "check-sun"))
    args = parser.parse_args()
    if args.action == "check-sun":
        print(json.dumps(check_sun(), indent=1))
        return
    record = build()
    RESULT.write_text(json.dumps(record, separators=(",", ":")) + "\n")
    print(RESULT)


if __name__ == "__main__":
    main()
