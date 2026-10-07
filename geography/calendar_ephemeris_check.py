"""Fetch and check a modest JPL Horizons fixture for the lighting calendar.

Run from the repository root (one BLAS thread)::

    OPENBLAS_NUM_THREADS=1 python -m geography.calendar_ephemeris_check --fetch
    OPENBLAS_NUM_THREADS=1 python -m geography.calendar_ephemeris_check

The default is an offline comparison with byte-pinned reference rows.  Fetching
is explicit, sequential, and fails without replacing the fixture if any request
fails.  No ephemeris coefficients are fitted to these checks.

All epochs are TT, including future dates.  Horizons converts TT to its internal
TDB (the periodic difference is below 2 ms).  The observer reference quantities
14/15 are *apparent* sub-observer/sub-solar coordinates on MOON_ME, seen from
Earth's center.  Their light-time conventions therefore differ slightly from
the instantaneous geometry returned by lunar_ephemeris.  Reported differences
include this effect as well as truncated terms and omitted physical libration.

Horizons MOON_ME is the DE421 mean-Earth/polar-axis standard (a fixed rotation
from the modern principal-axis solution), not IAU_MOON or the PA frame.  A fetch
rejects a response that substitutes another lunar orientation model.  Lunar
radius/oblateness do not enter the direction comparison.  Horizon-event checks
use a level spherical surface and the source-center direction at the Moon's
center; they exclude terrain, refraction, source radius and surface parallax.

Primary references:
https://ssd.jpl.nasa.gov/horizons/manual.html
https://ssd-api.jpl.nasa.gov/doc/horizons.html
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from datetime import datetime, timezone, timedelta
from pathlib import Path
import urllib.parse
import urllib.request

import numpy as np

from geography import lunar_ephemeris as ephemeris

HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "calendar_ephemeris_check.csv"
INPUT_PATH = HERE / "calendar_ephemeris_inputs.json"
RESULT_PATH = HERE / "calendar_ephemeris_check.json"
PHOTOMETRY_PATH = HERE.parent / "research/studies/sea_appearance/results/lighting_calendar.json"
API = "https://ssd.jpl.nasa.gov/api/horizons.api"
START = "2000-01-01"
STOP = "2501-01-01"
SPARSE_STEP_DAYS = 73
TWILIGHT_WINDOW = dict(sample="2030_twilight_month", start="2030-03-01",
                       stop="2030-04-03", step="30m")
# The first half-hour run found near-grazing 0.1-lux crossings here.  These
# predeclared follow-up windows keep the refinement reproducible and compact.
THRESHOLD_REFINEMENTS = [
    dict(sample="2030_procellarum_0p1_darkening", site="W Procellarum", threshold_lux=0.1,
         event="darkening", start="2030-03-10 14:00", stop="2030-03-10 16:00", step="1m"),
    dict(sample="2030_procellarum_0p1_brightening", site="W Procellarum", threshold_lux=0.1,
         event="brightening", start="2030-03-10 17:30", stop="2030-03-10 19:30", step="1m"),
]
# Geographic diversity, near/far sides and both hemispheres across five centuries.
EVENT_SITES = [
    (2000, 0.0, 0.0),
    (2050, 2.34, 93.4),
    (2100, 45.0, 0.0),
    (2200, 80.0, 0.0),
    (2300, -45.0, 180.0),
    (2400, -80.0, 45.0),
    (2500, 0.0, 180.0),
]
FIELDS = ["sample", "jd_tt", "sub_earth_lon_deg", "sub_earth_lat_deg",
          "sub_sun_lon_deg", "sub_sun_lat_deg", "earth_distance_km",
          "sun_distance_km"]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def calendar_jd(text):
    """Gregorian label expressed on TT; datetime is arithmetic only, not UTC."""
    d = datetime.fromisoformat(text)
    return 2451544.5 + (d - datetime(2000, 1, 1)).total_seconds() / 86400


def calendar_label(jd):
    seconds = round((float(jd) - 2451544.5) * 86400)
    return (datetime(2000, 1, 1) + timedelta(seconds=seconds)).isoformat(timespec="seconds")


def horizon_elevation(vectors, latitude_deg, longitude_deg):
    n = ephemeris.unit(math.radians(longitude_deg), math.radians(latitude_deg))
    return np.degrees(np.arcsin(np.clip(np.asarray(vectors) @ n, -1, 1)))


def bisect_event(start, stop, latitude, longitude, target_elevation=0.0):
    def f(jd):
        return float(horizon_elevation(ephemeris.earth_and_sun(jd)["sun"], latitude, longitude)) - target_elevation
    lo, hi, flo = start, stop, f(start)
    if flo * f(stop) > 0:
        raise ValueError("Event bracket does not straddle the requested elevation")
    for _ in range(40):
        mid = (lo + hi) / 2
        if np.sign(f(mid)) == np.sign(flo):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def event_windows():
    """Locate one rising and one setting center-horizon crossing per site/year."""
    windows = []
    for year, latitude, longitude in EVENT_SITES:
        start = calendar_jd(f"{year}-03-01")
        times = start + np.arange(33 * 4 + 1) / 4
        elev = horizon_elevation(ephemeris.earth_and_sun(times)["sun"], latitude, longitude)
        found = set()
        for i in np.flatnonzero(elev[:-1] * elev[1:] < 0):
            event = "sunrise" if elev[i + 1] > elev[i] else "sunset"
            if event in found:
                continue
            center = bisect_event(float(times[i]), float(times[i + 1]), latitude, longitude)
            # Align to an exact half-hour to avoid float/printed-epoch ambiguity.
            rounded = round((center - 0.5) * 48) / 48 + 0.5
            windows.append(dict(sample=f"{year}_{event}", latitude_deg=latitude,
                                longitude_deg=longitude, event=event,
                                model_event_jd_tt=center,
                                start_jd_tt=rounded - 0.5,
                                stop_jd_tt=rounded + 0.5,
                                step_minutes=30))
            found.add(event)
        if len(found) != 2:
            raise ValueError(f"Did not find both crossings for {year}, {latitude}, {longitude}")
    return windows


def request_horizons(parameters, responses):
    params = {"format": "json", "MAKE_EPHEM": "YES", "CSV_FORMAT": "YES",
              "TIME_TYPE": "TT", **parameters}
    url = API + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as response:
        raw = response.read()
    payload = json.loads(raw)
    result = payload.get("result", "")
    if payload.get("error") or "$$SOE" not in result or "$$EOE" not in result:
        raise RuntimeError(f"Horizons failed: {payload.get('error', result[:2000])}")
    if params["EPHEM_TYPE"] == "OBSERVER":
        header = next((s for s in result.splitlines() if s.startswith("Target pole/equ")), "")
        if "MOON_ME" not in header or "IAU_MOON" in header:
            raise ValueError(f"Expected MOON_ME throughout; got {header!r}")
    headers = [s for s in result.splitlines() if s.startswith(
        ("Target body name:", "Center body name:", "Target pole/equ", "Center-site name:",
         "Start time", "Stop  time", "Step-size", "Time format", "Reference frame"))]
    responses.append(dict(parameters=params, response_sha256=sha256(raw),
                          response_bytes=len(raw), signature=payload.get("signature"),
                          retrieved_utc=datetime.now(timezone.utc).isoformat(), headers=headers))
    table = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0]
    return [[v.strip() for v in row] for row in csv.reader(io.StringIO(table)) if row]


def observer_request(start, stop, step):
    return dict(COMMAND="'301'", CENTER="'500@399'", EPHEM_TYPE="OBSERVER",
                QUANTITIES="'14,15'", CAL_FORMAT="'JD'", START_TIME=f"'{start}'",
                STOP_TIME=f"'{stop}'", STEP_SIZE=f"'{step}'")


def parse_observer(rows, sample):
    return [dict(sample=sample, jd_tt=float(row[0]), sub_earth_lon_deg=float(row[3]),
                 sub_earth_lat_deg=float(row[4]), sub_sun_lon_deg=float(row[5]),
                 sub_sun_lat_deg=float(row[6]), earth_distance_km="", sun_distance_km="")
            for row in rows]


def fetch_fixture():
    responses = []
    records = parse_observer(request_horizons(
        observer_request(START, STOP, f"{SPARSE_STEP_DAYS}d"), responses), "73d")
    # A stop epoch is not necessarily included by a fixed-step query.
    end_jd = calendar_jd(STOP) - 1 / 48  # last half hour of 2500
    end_params = observer_request("2500-12-31 23:30", "2501-01-01", "30m")
    final_rows = parse_observer(request_horizons(end_params, responses), "endpoint")
    records.extend(r for r in final_rows if r["jd_tt"] <= end_jd + 1e-8)
    print(f"Fetched {len(records)} sparse/endpoint direction references", flush=True)
    # Independent geometric ranges, same TT epochs; no light-time range mismatch.
    sparse_times = np.array([r["jd_tt"] for r in records if r["sample"] == "73d"])
    for body, center, name in (("301", "500@399", "earth"), ("10", "500@301", "sun")):
        parameters = dict(COMMAND=f"'{body}'", CENTER=f"'{center}'", EPHEM_TYPE="VECTORS",
                          VEC_TABLE="'3'", VEC_CORR="'NONE'", OUT_UNITS="'KM-S'",
                          START_TIME=f"'{START}'", STOP_TIME=f"'{STOP}'",
                          STEP_SIZE=f"'{SPARSE_STEP_DAYS}d'")
        vector_rows = request_horizons(parameters, responses)
        epochs = np.array([float(r[0]) for r in vector_rows])
        if len(epochs) != len(sparse_times) or not np.allclose(epochs, sparse_times, atol=1e-9, rtol=0):
            raise ValueError(f"{name} vector epochs do not match observer TT epochs")
        for record, row in zip(records, vector_rows):
            # VEC_TABLE=3: JD, calendar, x/y/z, vx/vy/vz, light-time, range, range-rate.
            record[f"{name}_distance_km"] = float(row[9])
    windows = event_windows()
    for window in windows:
        start = f"JD{window['start_jd_tt']:.12f}"
        stop = f"JD{window['stop_jd_tt']:.12f}"
        rows = parse_observer(request_horizons(observer_request(start, stop, "30m"), responses),
                              window["sample"])
        if len(rows) != 49:
            raise ValueError(f"Expected 49 half-hour rows for {window['sample']}, got {len(rows)}")
        records.extend(rows)
        print(f"Fetched {window['sample']} ({len(rows)} rows)", flush=True)
    window = TWILIGHT_WINDOW
    rows = parse_observer(request_horizons(observer_request(
        window["start"], window["stop"], window["step"]), responses), window["sample"])
    records.extend(rows)
    print(f"Fetched {window['sample']} ({len(rows)} rows)", flush=True)
    for window in THRESHOLD_REFINEMENTS:
        rows = parse_observer(request_horizons(observer_request(
            window["start"], window["stop"], window["step"]), responses), window["sample"])
        records.extend(rows)
        print(f"Fetched {window['sample']} ({len(rows)} rows)", flush=True)
    save_fixture(records, responses, windows)


def save_fixture(records, responses, windows):
    records.sort(key=lambda r: (r["jd_tt"], r["sample"]))
    output = io.StringIO(newline="")
    output.write("# NASA/JPL Horizons; reference directions are apparent, MOON_ME, east-positive.\n")
    output.write("# All Julian dates are TT. Blank distance cells were not independently sampled.\n")
    output.write("# See calendar_ephemeris_inputs.json for queries, response hashes and event windows.\n")
    writer = csv.DictWriter(output, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    for record in records:
        writer.writerow({k: (f"{v:.12f}" if k == "jd_tt" else f"{v:.9f}")
                         if isinstance(v, (int, float)) else v for k, v in record.items()})
    data = output.getvalue().encode("ascii")
    provenance = dict(schema="terluna.calendar-ephemeris-reference.v1", source=API,
                      source_documentation=["https://ssd.jpl.nasa.gov/horizons/manual.html",
                                            "https://ssd-api.jpl.nasa.gov/doc/horizons.html"],
                      credit="NASA/JPL Solar System Dynamics, Horizons; JPL ephemerides are public domain.",
                      fixture_sha256=sha256(data), model_sha256_at_window_selection=sha256(
                          Path(ephemeris.__file__).read_bytes()),
                      frame="MOON_ME (DE421 mean-Earth/polar-axis standard)",
                      time_scale="TT; Horizons internally converts TT to TDB",
                      sparse_step_days=SPARSE_STEP_DAYS, event_windows=windows,
                      twilight_window=TWILIGHT_WINDOW,
                      threshold_refinements=THRESHOLD_REFINEMENTS, responses=responses)
    CSV_PATH.write_bytes(data)
    INPUT_PATH.write_text(json.dumps(provenance, indent=2) + "\n")


def angular_separation_deg(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), np.sum(a * b, axis=-1)))


def errors_summary(error, epochs):
    err = np.asarray(error, dtype=float)
    index = int(np.argmax(np.abs(err)))
    return dict(max_abs=float(np.abs(err[index])), rms=float(np.sqrt(np.mean(err**2))),
                p95_abs=float(np.percentile(np.abs(err), 95)),
                maximum_at_jd_tt=float(epochs[index]), maximum_at_tt=calendar_label(epochs[index]))


def linear_crossing(epochs, elevations):
    indices = np.flatnonzero(elevations[:-1] * elevations[1:] <= 0)
    if len(indices) != 1:
        raise ValueError(f"Expected one horizon crossing; found {len(indices)}")
    i = int(indices[0])
    return float(epochs[i] - elevations[i] * (epochs[i + 1] - epochs[i]) /
                 (elevations[i + 1] - elevations[i]))


def photometry_comparison(records, epochs, model_vectors, reference_vectors, path=PHOTOMETRY_PATH):
    """Propagate only the direction residual through the committed clear-sky curve."""
    raw = Path(path).read_bytes()
    product = json.loads(raw)
    curve = product["sky"]["ground_lux_by_sun_elevation"]
    angles, lux = np.array(curve["elevation_deg"]), np.array(curve["lux"])
    if not np.all(np.diff(angles) > 0) or not np.all(np.diff(lux) > 0):
        raise ValueError("Photometry check requires a strictly increasing elevation/lux curve")
    sites = [dict(name=s["short"], latitude_deg=s["lon_lat"][1], longitude_deg=s["lon_lat"][0])
             for s in product["sites"]]
    sites += [dict(name=f"{lat:+g}N {lon:g}E check", latitude_deg=lat, longitude_deg=lon)
              for lat, lon in ((0, 0), (45, 0), (60, 0), (80, 0), (-80, 45))]
    bands = [("daylight", 0, 90.000001), ("horizon_twilight", -6, 0),
             ("bright_twilight", -30, -6), ("deep_twilight", -60, -30),
             ("deep_night", -90.000001, -60)]
    month_mask = np.array([r["sample"] == TWILIGHT_WINDOW["sample"] for r in records])
    month_times = epochs[month_mask]
    summaries, threshold_events, event_count_differences = [], [], []
    for site in sites:
        lat, lon = site["latitude_deg"], site["longitude_deg"]
        predicted_el = horizon_elevation(model_vectors, lat, lon)
        reference_el = horizon_elevation(reference_vectors, lat, lon)
        predicted_lux = np.interp(predicted_el, angles, lux)
        reference_lux = np.interp(reference_el, angles, lux)
        delta = predicted_lux - reference_lux
        rel = delta / reference_lux
        band_rows = []
        for name, lo, hi in bands:
            selected = np.flatnonzero((reference_el >= lo) & (reference_el < hi))
            if not len(selected):
                continue
            worst = int(selected[np.argmax(np.abs(rel[selected]))])
            band_rows.append(dict(band=name, samples=len(selected),
                                  max_abs_error_lux=float(np.max(np.abs(delta[selected]))),
                                  max_relative_error_percent=float(np.max(np.abs(rel[selected])) * 100),
                                  rms_relative_error_percent=float(np.sqrt(np.mean(rel[selected]**2)) * 100),
                                  relative_worst=dict(jd_tt=float(epochs[worst]), tt=calendar_label(epochs[worst]),
                                                      reference_elevation_deg=float(reference_el[worst]),
                                                      model_elevation_deg=float(predicted_el[worst]),
                                                      reference_lux=float(reference_lux[worst]),
                                                      model_lux=float(predicted_lux[worst]))))
        summaries.append({**site, "bands": band_rows})
        if not len(month_times):
            continue
        month_reference_el = reference_el[month_mask]
        month_model_el = predicted_el[month_mask]
        for threshold in (100.0, 10.0, 1.0, 0.1):
            target = float(np.interp(threshold, lux, angles))
            shifted = month_reference_el - target
            indices = np.flatnonzero(shifted[:-1] * shifted[1:] <= 0)
            model_shifted = month_model_el - target
            model_indices = np.flatnonzero(model_shifted[:-1] * model_shifted[1:] <= 0)
            if len(indices) != len(model_indices):
                event_count_differences.append(dict(site=site["name"], threshold_lux=threshold,
                                                   reference_crossings=len(indices), model_crossings=len(model_indices)))
            for index in indices:
                i = int(index)
                reference_jd = linear_crossing(month_times[i:i+2], shifted[i:i+2])
                rising = shifted[i + 1] > shifted[i]
                event = "brightening" if rising else "darkening"
                step_minutes = 30
                # Use a one-minute follow-up when the original threshold grazes
                # a brightness extremum, where coarse timing is ill-conditioned.
                refinement = next((w for w in THRESHOLD_REFINEMENTS if w["site"] == site["name"]
                                   and w["threshold_lux"] == threshold and w["event"] == event), None)
                fine_mask = np.array([r["sample"] == refinement["sample"] for r in records]) if refinement else None
                if fine_mask is not None and fine_mask.any():
                    fine_times = epochs[fine_mask]
                    fine_shifted = reference_el[fine_mask] - target
                    reference_jd = linear_crossing(fine_times, fine_shifted)
                    coarse_jd = linear_crossing(fine_times[::2], fine_shifted[::2])
                    step_minutes = 1
                else:
                    coarse_index = min(i // 2 * 2, len(month_times) - 3)
                    coarse_jd = linear_crossing(month_times[coarse_index:coarse_index+3:2],
                                               shifted[coarse_index:coarse_index+3:2])
                candidates = [int(j) for j in model_indices if (model_shifted[j + 1] > model_shifted[j]) == rising]
                if candidates:
                    j = min(candidates, key=lambda j: abs(month_times[j] - reference_jd))
                    # Match within the same event, not one lunar month away.
                    model_jd = bisect_event(float(month_times[j]), float(month_times[j+1]), lat, lon, target) if abs(month_times[j]-reference_jd)<1 else None
                else:
                    model_jd = None
                threshold_events.append(dict(site=site["name"], threshold_lux=threshold,
                                             threshold_elevation_deg=target,
                                             event=event,
                                             reference_jd_tt=reference_jd, reference_tt=calendar_label(reference_jd),
                                             model_jd_tt=model_jd, model_tt=calendar_label(model_jd) if model_jd else None,
                                             model_minus_reference_minutes=(model_jd-reference_jd)*1440 if model_jd else None,
                                             reference_step_minutes=step_minutes,
                                             coarser_reference_step_minutes=step_minutes*2,
                                             reference_interpolation_change_seconds=abs(reference_jd-coarse_jd)*86400))
    return dict(source=str(Path(path).relative_to(HERE.parent)), sha256=sha256(raw),
                source_evidence=product["evidence"], interpolation="linear in lux versus elevation, numpy.interp",
                reading_rule="Direction-error propagation within the fixed clear-sky solar-only photometry curve; this does not validate atmospheric radiative transfer. Earthlight, distance scaling, clouds, terrain, refraction, and source size are held out. Sampling statistics are not climatological averages.",
                sites=summaries, threshold_events=threshold_events, event_count_differences=event_count_differences,
                threshold_timing_span_calendar_tt=[TWILIGHT_WINDOW["start"], TWILIGHT_WINDOW["stop"]],
                threshold_event_abs_error_max_minutes=max((abs(e["model_minus_reference_minutes"]) for e in threshold_events if e["model_minus_reference_minutes"] is not None), default=None),
                threshold_interpolation_sensitivity_max_seconds=max((e["reference_interpolation_change_seconds"] for e in threshold_events), default=None))


def compare_fixture(csv_path=CSV_PATH, input_path=INPUT_PATH):
    data = Path(csv_path).read_bytes()
    metadata = json.loads(Path(input_path).read_text())
    if sha256(data) != metadata["fixture_sha256"]:
        raise ValueError("Pinned Horizons fixture hash mismatch")
    records = list(csv.DictReader(line for line in data.decode("ascii").splitlines()
                                  if not line.startswith("#")))
    epochs = np.array([float(r["jd_tt"]) for r in records])
    model = ephemeris.earth_and_sun(epochs)
    summaries, refs = {}, {}
    for body, prefix in (("earth", "sub_earth"), ("sun", "sub_sun")):
        ref_lon = np.array([float(r[f"{prefix}_lon_deg"]) for r in records])
        ref_lat = np.array([float(r[f"{prefix}_lat_deg"]) for r in records])
        refs[body] = ephemeris.unit(np.radians(ref_lon), np.radians(ref_lat))
        lon, lat = model[f"{body}_lon_lat_deg"]
        lon_error = (lon - ref_lon + 180) % 360 - 180
        summaries[f"{body}_longitude_deg"] = errors_summary(lon_error, epochs)
        summaries[f"{body}_latitude_deg"] = errors_summary(lat - ref_lat, epochs)
        summaries[f"{body}_direction_deg"] = errors_summary(
            angular_separation_deg(model[body], refs[body]), epochs)
        mask = np.array([bool(r[f"{body}_distance_km"]) for r in records])
        ranges = np.array([float(r[f"{body}_distance_km"]) for r in records if r[f"{body}_distance_km"]])
        range_error = model[f"{body}_distance_m"][mask] / 1000 - ranges
        summaries[f"{body}_distance_km"] = dict(samples=int(mask.sum()), **errors_summary(range_error, epochs[mask]))
        summaries[f"{body}_inverse_square_flux_relative"] = errors_summary(
            (ranges / (model[f"{body}_distance_m"][mask] / 1000))**2 - 1, epochs[mask])
    events = []
    for window in metadata["event_windows"]:
        mask = np.array([r["sample"] == window["sample"] for r in records])
        times = epochs[mask]
        lat, lon = window["latitude_deg"], window["longitude_deg"]
        reference_elev = horizon_elevation(refs["sun"][mask], lat, lon)
        predicted_elev = horizon_elevation(model["sun"][mask], lat, lon)
        reference_time = linear_crossing(times, reference_elev)
        # Recompute exact model crossing after code updates; do not reuse fetch-time result.
        predicted_time = bisect_event(float(times.min()), float(times.max()), lat, lon)
        # Removing alternate reference samples estimates event interpolation sensitivity.
        coarse_time = linear_crossing(times[::2], reference_elev[::2])
        events.append(dict(sample=window["sample"], latitude_deg=lat, longitude_deg=lon,
                           reference_event_jd_tt=reference_time, reference_event_tt=calendar_label(reference_time),
                           model_event_jd_tt=predicted_time, model_event_tt=calendar_label(predicted_time),
                           model_minus_reference_minutes=(predicted_time - reference_time) * 1440,
                           interpolation_30_vs_60_minutes_seconds=abs(reference_time - coarse_time) * 86400,
                           elevation_error_max_deg=float(np.abs(predicted_elev - reference_elev).max())))
    decades = []
    for lo, hi in ((2000, 2100), (2100, 2200), (2200, 2300), (2300, 2400), (2400, 2501)):
        mask = (epochs >= calendar_jd(f"{lo}-01-01")) & (epochs < calendar_jd(f"{hi}-01-01"))
        decades.append(dict(start_year=lo, stop_year_exclusive=hi, samples=int(mask.sum()),
                            sun_direction_max_deg=float(angular_separation_deg(model["sun"][mask], refs["sun"][mask]).max()),
                            earth_direction_max_deg=float(angular_separation_deg(model["earth"][mask], refs["earth"][mask]).max())))
    return dict(schema="terluna.calendar-ephemeris-check.v1", evidence="sampled independent numerical comparison",
                producer="geography/calendar_ephemeris_check.py", producer_sha256=sha256(Path(__file__).read_bytes()),
                model="geography/lunar_ephemeris.py", model_sha256=sha256(Path(ephemeris.__file__).read_bytes()),
                fixture_sha256=sha256(data), inputs_sha256=sha256(Path(input_path).read_bytes()),
                samples=len(records), independent_range_samples=summaries["earth_distance_km"]["samples"],
                unique_epochs=int(len(np.unique(epochs))), span_jd_tt=[float(epochs.min()), float(epochs.max())],
                span_calendar_tt=[calendar_label(epochs.min()), calendar_label(epochs.max())],
                frame=metadata["frame"], time_scale=metadata["time_scale"],
                errors=summaries, century_groups=decades, horizon_events=events,
                horizon_event_abs_error_max_minutes=max(abs(e["model_minus_reference_minutes"]) for e in events),
                event_interpolation_sensitivity_max_seconds=max(e["interpolation_30_vs_60_minutes_seconds"] for e in events),
                photometry=photometry_comparison(records, epochs, model["sun"], refs["sun"]),
                reading_rule="Observed maxima at the saved samples, not rigorous all-time error bounds. Direction discrepancies include geometric/apparent convention differences. Ground brightness is not validated by this check.",
                calendar_rule="Epochs are Gregorian labels on TT. For UTC input, add 32.184 s plus the date's TAI-UTC. Future leap seconds are unspecified; a chosen future UTC offset is a calendar convention, not a validated prediction.",
                horizon_rule="Sun-center crossing of an unobstructed level horizon, Moon-center direction. Terrain, refraction, solar radius and surface parallax are excluded. Near-polar grazing events outside these windows can have larger timing errors.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true", help="Refresh pinned public JPL reference rows (network required)")
    args = parser.parse_args()
    if args.fetch:
        fetch_fixture()
    result = compare_fixture()
    RESULT_PATH.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("samples", "span_calendar_tt", "errors",
                     "horizon_event_abs_error_max_minutes", "event_interpolation_sensitivity_max_seconds")}, indent=2))


if __name__ == "__main__":
    main()
