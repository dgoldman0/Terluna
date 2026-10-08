"""Audit the exported finite-source calendar response against pinned JPL inputs.

Run after ``illumination.calendar.build --export``::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m illumination.calendar.validate

No atmospheric simulation or network request is made.  The independent Python
reader below uses disk.samples and linear interpolation of the exported tables.
It emits portable golden cases for the JavaScript reader, a source-disk
quadrature refinement check and sampled geometry-to-brightness discrepancies.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from geography import lunar_ephemeris
from illumination.calendar import disk
from illumination.calendar.geometry import geometry, geometry_from_vectors

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TRANSFER = HERE / "results/transfer.json"
REFERENCE = ROOT / "geography/calendar_ephemeris_check.csv"
REFERENCE_INPUTS = ROOT / "geography/calendar_ephemeris_inputs.json"
HISTORICAL_SOLVED = ROOT / "illumination/sky/results/solved_sky.json"
HISTORICAL_CURVE = ROOT / "research/studies/sea_appearance/results/lighting_calendar.json"
OUTPUT = HERE / "results/validation.json"
SCHEMA = "terluna.illumination.calendar-validation/1"
COMPONENTS = ("solar_direct", "solar_diffuse", "solar", "earth_direct", "earth_diffuse", "earth", "total")
SITES = [
    dict(name="Western Procellarum", longitude=-73.875, latitude=28.875),
    dict(name="Southern Imbrium", longitude=-17.875, latitude=15.875),
    dict(name="Southern Nubium", longitude=-11.875, latitude=-28.125),
    dict(name="Nectaris", longitude=33.125, latitude=-22.125),
    dict(name="Smythii headland", longitude=93.4, latitude=2.34),
    dict(name="Ingenii coast", longitude=162.875, latitude=-37.125),
    dict(name="Equatorial near side", longitude=0.0, latitude=0.0),
    dict(name="Mean eastern limb", longitude=90.0, latitude=0.0),
    dict(name="Equatorial far side", longitude=180.0, latitude=0.0),
    dict(name="Northern high latitude", longitude=0.0, latitude=80.0),
    dict(name="Southern high latitude", longitude=45.0, latitude=-80.0),
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Response:
    """Python implementation independent of the JavaScript interpolation code."""

    def __init__(self, data):
        if data.get("schema") != "terluna.illumination.calendar-transfer/1":
            raise ValueError("Unsupported lighting data")
        self.constants = data["constants"]
        self.direct_angles = np.asarray(data["direct_elevation_deg"], dtype=float)
        self.diffuse_angles = np.asarray(data["diffuse_elevation_deg"], dtype=float)
        self.phase_angles = np.asarray(data["phase_deg"], dtype=float)
        for angles in (self.direct_angles, self.diffuse_angles, self.phase_angles):
            if len(angles) < 2 or not np.all(np.isfinite(angles)) or not np.all(np.diff(angles) > 0):
                raise ValueError("Response nodes must be finite and strictly increasing")
        for angles in (self.direct_angles, self.diffuse_angles):
            if angles[0] > -90 or angles[-1] < 90:
                raise ValueError("The calendar response must cover source elevations -90 through +90 degrees")
        if self.phase_angles[0] > 0 or self.phase_angles[-1] < 180:
            raise ValueError("The Earth response must cover phases 0 through 180 degrees")
        self.solar_direct = np.asarray(data["solar"]["direct_lux"], dtype=float)
        self.solar_diffuse = np.asarray(data["solar"]["diffuse_lux"], dtype=float)
        self.earth_direct = np.asarray(data["earth"]["direct_lux"], dtype=float)
        self.earth_diffuse = np.asarray(data["earth"]["diffuse_lux"], dtype=float)
        for values, shape in ((self.solar_direct, self.direct_angles.shape),
                              (self.solar_diffuse, self.diffuse_angles.shape),
                              (self.earth_direct, (len(self.phase_angles), len(self.direct_angles))),
                              (self.earth_diffuse, (len(self.phase_angles), len(self.diffuse_angles)))):
            if values.shape != shape or not np.all(np.isfinite(values)) or np.any(values < 0):
                raise ValueError("Response arrays must have the documented shape and nonnegative finite lux")

    def light(self, g, include_earth=True, order=12):
        rays, weights = disk.samples(g["sun_elevation_deg"], g["sun_radius_deg"], order=order)
        solar_direct = float(weights @ np.interp(rays, self.direct_angles, self.solar_direct)) * g["sun_distance_factor"]
        solar_diffuse = float(weights @ np.interp(rays, self.diffuse_angles, self.solar_diffuse)) * g["sun_distance_factor"]
        earth_direct, earth_diffuse = 0.0, 0.0
        phase = g["earth_phase_deg"]
        if include_earth and phase < 180:
            p = int(np.clip(np.searchsorted(self.phase_angles, phase, side="right")-1, 0, len(self.phase_angles)-2))
            fraction = float(np.clip((phase-self.phase_angles[p])/(self.phase_angles[p+1]-self.phase_angles[p]), 0, 1))
            rays, weights = disk.samples(g["earth_elevation_deg"], g["earth_radius_deg"], phase,
                                         g["earth_bright_limb_rad"], order)
            direct = ((1-fraction)*np.interp(rays, self.direct_angles, self.earth_direct[p]) +
                      fraction*np.interp(rays, self.direct_angles, self.earth_direct[p+1]))
            diffuse = ((1-fraction)*np.interp(rays, self.diffuse_angles, self.earth_diffuse[p]) +
                       fraction*np.interp(rays, self.diffuse_angles, self.earth_diffuse[p+1]))
            radius = self.constants["earth_radius_m"]
            def solid_angle(distance):
                q = (radius/distance)**2
                return 2*math.pi*q/(1+math.sqrt(1-q))
            factor = (solid_angle(g["earth_distance_m"])/solid_angle(self.constants["earth_reference_distance_m"]) *
                      g["earth_sunlight_factor"])
            earth_direct = float(weights @ direct) * factor
            earth_diffuse = float(weights @ diffuse) * factor
        solar, earth = solar_direct+solar_diffuse, earth_direct+earth_diffuse
        return dict(solar_direct=float(solar_direct), solar_diffuse=float(solar_diffuse), solar=float(solar),
                    earth_direct=float(earth_direct), earth_diffuse=float(earth_diffuse), earth=float(earth),
                    total=float(solar+earth), eclipse=bool(g["source_separation_deg"] < g["earth_radius_deg"]+g["sun_radius_deg"]),
                    earthIncluded=bool(include_earth))


def scalar_geometry(g, index=None):
    return {k:float(v if index is None else v[index]) for k,v in g.items()}


def read_reference():
    metadata = json.loads(REFERENCE_INPUTS.read_text())
    if digest(REFERENCE) != metadata["fixture_sha256"]:
        raise ValueError("Pinned JPL fixture hash mismatch")
    rows = list(csv.DictReader(line for line in REFERENCE.read_text().splitlines() if not line.startswith("#")))
    # The horizon and finer threshold windows intentionally overlap the month.
    unique = {}
    for row in rows:
        key = float(row["jd_tt"])
        if key not in unique or row["earth_distance_km"]:
            unique[key] = row
    records = [unique[key] for key in sorted(unique)]
    epochs = np.array([float(r["jd_tt"]) for r in records])
    native = lunar_ephemeris.earth_and_sun(epochs)
    direction_only = {k:np.array(v, copy=True) for k,v in native.items() if k in ("earth", "sun", "earth_distance_m", "sun_distance_m")}
    for body in ("earth", "sun"):
        direction_only[body] = lunar_ephemeris.unit(
            np.radians([float(r[f"sub_{body}_lon_deg"]) for r in records]),
            np.radians([float(r[f"sub_{body}_lat_deg"]) for r in records]))
    ranged = {k:np.array(v, copy=True) for k,v in direction_only.items()}
    range_mask = np.array([bool(r["earth_distance_km"] and r["sun_distance_km"]) for r in records])
    for body in ("earth", "sun"):
        ranged[f"{body}_distance_m"][range_mask] = np.array([float(r[f"{body}_distance_km"])*1000 for r in records if r[f"{body}_distance_km"]])
    return records, epochs, native, direction_only, ranged, range_mask


def metric(actual, expected, descriptions):
    a, b = np.asarray(actual), np.asarray(expected)
    delta = a-b
    absolute_index = int(np.argmax(np.abs(delta)))
    result = dict(samples=len(delta), max_abs_lux=float(np.max(np.abs(delta))),
                  rms_lux=float(np.sqrt(np.mean(delta*delta))),
                  absolute_worst={**descriptions[absolute_index], "native_lux":float(a[absolute_index]),
                                  "reference_lux":float(b[absolute_index])})
    # Retain weak-source relative errors and separately expose meaningful-light
    # cases. A finite denominator floor is selection, never added illumination.
    result["relative"] = []
    for floor in (1e-8, 0.001, 0.01, 0.1):
        selected = np.flatnonzero(b >= floor)
        if len(selected):
            relative = delta[selected]/b[selected]
            index = int(selected[np.argmax(np.abs(relative))])
            result["relative"].append(dict(reference_at_least_lux=floor, samples=len(selected),
                                           max_abs_percent=float(np.max(np.abs(relative))*100),
                                           rms_percent=float(np.sqrt(np.mean(relative*relative))*100),
                                           worst={**descriptions[index], "native_lux":float(a[index]),
                                                  "reference_lux":float(b[index])}))
    return result


def describe(site, jd, native_g, reference_g):
    return dict(site=site["name"], longitude=site["longitude"], latitude=site["latitude"], jd_tt=float(jd),
                native_sun_elevation_deg=native_g["sun_elevation_deg"], reference_sun_elevation_deg=reference_g["sun_elevation_deg"],
                native_earth_elevation_deg=native_g["earth_elevation_deg"], reference_earth_elevation_deg=reference_g["earth_elevation_deg"],
                native_earth_phase_deg=native_g["earth_phase_deg"], reference_earth_phase_deg=reference_g["earth_phase_deg"])


def geometry_audit(response, epochs, native, directions, ranged, range_mask, order):
    groups = {name:dict(actual=[], expected=[], descriptions=[], eclipses_excluded=0)
              for name in ("directions_only", "directions_and_geometric_ranges")}
    golden_candidates = []
    for site in SITES:
        ng = geometry_from_vectors(native, site["longitude"], site["latitude"])
        dg = geometry_from_vectors(directions, site["longitude"], site["latitude"])
        rg = geometry_from_vectors(ranged, site["longitude"], site["latitude"])
        near_solar = set(np.argsort(np.abs(ng["sun_elevation_deg"]))[:2].tolist())
        near_earth = set(np.argsort(np.abs(ng["earth_elevation_deg"]))[:2].tolist())
        for i, jd in enumerate(epochs):
            g = scalar_geometry(ng, i)
            expected_g = scalar_geometry(dg, i)
            actual = response.light(g, order=order)
            expected = response.light(expected_g, order=order)
            label = describe(site, jd, g, expected_g)
            target = groups["directions_only"]
            if actual["eclipse"] or expected["eclipse"]:
                target["eclipses_excluded"] += 1
            else:
                target["actual"].append([actual[k] for k in COMPONENTS])
                target["expected"].append([expected[k] for k in COMPONENTS])
                target["descriptions"].append(label)
            if range_mask[i]:
                ranged_g = scalar_geometry(rg, i)
                ranged_light = response.light(ranged_g, order=order)
                target = groups["directions_and_geometric_ranges"]
                if actual["eclipse"] or ranged_light["eclipse"]:
                    target["eclipses_excluded"] += 1
                else:
                    target["actual"].append([actual[k] for k in COMPONENTS])
                    target["expected"].append([ranged_light[k] for k in COMPONENTS])
                    target["descriptions"].append(describe(site, jd, g, ranged_g))
            # Selected physically coupled cases for the portable reader, including
            # its Earth-off switch and actual near-horizon phase orientations.
            if (i in near_solar or i in near_earth) and not actual["eclipse"]:
                golden_candidates.append(dict(name=f"{site['name']} at JD {jd:.8f}", kind="native_ephemeris",
                                               jd_tt=float(jd), longitude=site["longitude"], latitude=site["latitude"],
                                               geometry=g, includeEarth=True, order=order, expected=actual))
        print(f"Audited {site['name']}: {len(epochs)} epochs", flush=True)
    output = {}
    for name, values in groups.items():
        actual, expected = np.asarray(values["actual"]), np.asarray(values["expected"])
        output[name] = dict(cases=len(actual), eclipse_cases_excluded=values["eclipses_excluded"],
                            components={k:metric(actual[:,j], expected[:,j], values["descriptions"])
                                        for j,k in enumerate(COMPONENTS)})
    return output, golden_candidates


def quadrature_audit(response):
    base = scalar_geometry(geometry(2462570.5, 90.0, 0.0))
    # Synthetic source settings isolate the integration, including the bright
    # crescent both above and below the horizon. They are not dated scenes.
    cases = []
    for elevation in (-0.6, -0.3, -0.267, -0.15, 0, 0.15, 0.267, 0.3, 0.6, 1, 2):
        g = {**base, "sun_elevation_deg":float(elevation), "sun_radius_deg":0.267, "source_separation_deg":90.0}
        cases.append(dict(name=f"Solar center {elevation:g} deg", source="sun", geometry=g, includeEarth=False))
    for phase in (0, 30, 60, 90, 120, 150, 170, 179):
        for elevation in (-2, -1.1, -1, -0.8, -0.3, 0, 0.3, 0.8, 1, 1.1, 2):
            for beta in (0, math.pi/2, math.pi):
                g = {**base, "sun_elevation_deg":-70.0, "earth_elevation_deg":float(elevation),
                     "earth_phase_deg":float(phase), "earth_radius_deg":1.0,
                     "earth_bright_limb_rad":float(beta), "source_separation_deg":90.0}
                cases.append(dict(name=f"Earth phase {phase:g}, center {elevation:g}, limb {beta:g}",
                                  source="earth", geometry=g, includeEarth=True))
    low, high, labels = [], [], []
    golden = []
    for i, case in enumerate(cases):
        a = response.light(case["geometry"], case["includeEarth"], order=12)
        b = response.light(case["geometry"], case["includeEarth"], order=24)
        low.append([a[k] for k in COMPONENTS]); high.append([b[k] for k in COMPONENTS])
        labels.append(dict(case=case["name"], source=case["source"], geometry=case["geometry"]))
        if i < 11 or (case["geometry"]["earth_elevation_deg"] == 0 and case["geometry"]["earth_phase_deg"] in (0, 90, 179)):
            golden.append(dict(name=case["name"], kind="synthetic_quadrature", geometry=case["geometry"],
                               includeEarth=case["includeEarth"], order=12, expected=a))
    low, high = np.asarray(low), np.asarray(high)
    return dict(cases=len(cases), low_order=12, reference_order=24,
                evidence="Quadrature-order sensitivity for the same exported response tables; order 24 is a numerical comparator, not exact truth.",
                components={k:metric(low[:,j], high[:,j], labels) for j,k in enumerate(COMPONENTS)}), golden


def compact_settings(settings):
    keys = ("quality", "albedo", "increment_tolerance", "order_cap", "order_converged",
            "solar_disk_radius_deg", "sun_strips", "convergence_min_elevation_deg",
            "angular_order_per_interval", "azimuth_samples", "continuation_runner", "continued_orders")
    output = {k:settings[k] for k in keys if k in settings}
    orders = settings.get("channel_orders", [])
    if orders:
        output["channel_order_range"] = [min(orders), max(orders)]
        output["channels"] = len(orders)
    return output


def historical_solar_regression(response, transfer, solved_path=HISTORICAL_SOLVED, curve_path=HISTORICAL_CURVE):
    """Sample only saved historical elevations, at the same fixed source size.

    The published atlas flux and the moment readout are separate references.
    The sea curve has a further direct-beam interpolation convention, recorded
    below.  Their differences are regression observations rather than a universal
    acceptance threshold or an independent observational validation.
    """
    solved = json.loads(Path(solved_path).read_text())
    curve_product = json.loads(Path(curve_path).read_text())
    if solved.get("schema") != "terluna.illumination.solved-spherical-sky/1":
        raise ValueError("Unsupported historical solved-sky product")
    if curve_product.get("schema") != "terluna.research.sea-lighting-calendar/1":
        raise ValueError("Unsupported historical sea-calendar product")
    moon = solved["worlds"]["moon_1.2atm"]
    settings = moon["settings"]
    radius = float(settings["solar_disk_radius_deg"])
    base = scalar_geometry(geometry(2462570.5, 0.0, 0.0))
    base.update(sun_radius_deg=radius, sun_distance_factor=1.0,
                sun_distance_m=float(solved["producer"]["constants"]["AU"]),
                source_separation_deg=90.0)
    def evaluate(elevations):
        elevations = np.asarray(elevations, dtype=float)
        if not np.all(np.isfinite(elevations)) or not np.all(np.diff(elevations) > 0):
            raise ValueError("Historical elevations must be finite and strictly increasing")
        return [response.light({**base,"sun_elevation_deg":float(e)}, include_earth=False, order=12)
                for e in elevations]
    native_rows = sorted(moon["samples"], key=lambda row: row["sun_deg"])
    elevations = np.array([r["sun_deg"] for r in native_rows])
    values = evaluate(elevations)
    labels = [dict(sun_elevation_deg=float(e)) for e in elevations]
    direct = np.array([r["direct_horizontal_lux"] for r in native_rows])
    moment = np.array([r["moment_diffuse_horizontal_lux"] for r in native_rows])
    atlas = np.array([r["total_horizontal_lux"] for r in native_rows])
    predicted_direct = np.array([r["solar_direct"] for r in values])
    predicted_diffuse = np.array([r["solar_diffuse"] for r in values])
    predicted_total = np.array([r["solar"] for r in values])
    published = dict(samples=len(elevations), supported_elevation_span_deg=[float(elevations.min()),float(elevations.max())],
                     components=dict(direct=metric(predicted_direct,direct,labels),
                                     diffuse_against_moments=metric(predicted_diffuse,moment,labels),
                                     total_against_moments=metric(predicted_total,direct+moment,labels),
                                     total_against_published_atlas=metric(predicted_total,atlas,labels)),
                     curve=dict(elevation_deg=elevations.tolist(), calendar_direct_lux=predicted_direct.tolist(),
                                calendar_diffuse_lux=predicted_diffuse.tolist(), calendar_solar_lux=predicted_total.tolist(),
                                historical_direct_lux=direct.tolist(), historical_moment_diffuse_lux=moment.tolist(),
                                historical_atlas_total_lux=atlas.tolist()))
    curve = curve_product["sky"]["ground_lux_by_sun_elevation"]
    sea_elevations = np.asarray(curve["elevation_deg"], dtype=float)
    sea_reference = np.asarray(curve["lux"], dtype=float)
    if sea_reference.shape != sea_elevations.shape or not np.all(np.isfinite(sea_reference)) or np.any(sea_reference < 0):
        raise ValueError("Historical solar curve requires matching nonnegative finite lux")
    sea_values = np.array([r["solar"] for r in evaluate(sea_elevations)])
    sea_labels = [dict(sun_elevation_deg=float(e)) for e in sea_elevations]
    sea_bands = []
    for name, low, high in (("deep_night",-90.000001,-60), ("deep_twilight",-60,-30),
                           ("twilight",-30,0), ("daylight",0,90.000001)):
        mask = (sea_elevations >= low) & (sea_elevations < high)
        if np.any(mask):
            sea_bands.append(dict(band=name, **metric(sea_values[mask], sea_reference[mask],
                                                     [label for label, keep in zip(sea_labels, mask) if keep])))
    sea = dict(samples=len(sea_elevations), supported_elevation_span_deg=[float(sea_elevations.min()),float(sea_elevations.max())],
               total=metric(sea_values,sea_reference,sea_labels), bands=sea_bands,
               curve=dict(elevation_deg=sea_elevations.tolist(), calendar_solar_lux=sea_values.tolist(),
                          historical_solar_lux=sea_reference.tolist()),
               historical_last_scattering_order_fraction_max=curve_product["sky"]["last_scattering_order_fraction_max"])
    return dict(evidence="Sampled regression between distinct numerical producers on their saved elevation nodes; no universal pass/fail tolerance.",
                normalization=dict(sun_distance="1 AU, inverse-square factor fixed to one", solar_disk_radius_deg=radius,
                                   earth_included=False, disk_quadrature_order=12),
                historical=dict(solved_product_schema=solved["schema"], solved_producer=solved["producer"],
                                sea_curve_schema=curve_product["schema"], sea_curve_producer=curve_product["producer"],
                                transport_settings=compact_settings(settings),
                                disk_method="12 strips across the visible solar disk for attenuation; the diffuse scattering phase uses the disk-center direction.",
                                convergence_rule="Surface increment test applies to elevations at least -24 degrees, together with a whole-field maximum test; retained tolerance 2e-5. One retained channel continued beyond the initial order cap to order 104.",
                                published_readout="Published diffuse flux integrates the directional atlas; the product separately retains moment-grid diffuse flux. Both are compared.",
                                sea_curve_readout="Moment diffuse flux linearly interpolated in elevation, plus a log-interpolated direct-normal beam times sin(elevation), set to zero at and below zero elevation. No stellar/airglow floor is included in this solar-only curve."),
                calendar=dict(producer=transfer["producer"], transport_settings=compact_settings(transfer["transport_settings"]),
                              disk_method="Point-source transport convolved over a finite uniform disk after export; direct and diffuse source-elevation responses are both integrated.",
                              convergence_rule="The point-source surface increment test spans the full -90 to +90 degree source-elevation grid, with the recorded tolerance and whole-field maximum test."),
                interpretation="A difference includes source-disk treatment, retained scattering orders, readout and interpolation changes. Agreement would be numerical consistency, not independent physical validation. No reference is extrapolated beyond its saved elevation support.",
                published_solved_sky=published, sea_calendar_solar_curve=sea)


def validate(transfer_path=TRANSFER, output_path=OUTPUT):
    transfer_path, output_path = Path(transfer_path), Path(output_path)
    if not transfer_path.exists():
        raise FileNotFoundError(f"Run the transport export before this audit: {transfer_path}")
    data = json.loads(transfer_path.read_text())
    response = Response(data)
    # A stale export can still be internally interpolated, but cannot substantiate
    # the current implementation. Fail explicitly rather than re-pinning it.
    for name, expected in data["producer"]["files"].items():
        if digest(ROOT/name) != expected:
            raise ValueError(f"Transfer producer changed since export: {name}")
    records, epochs, native, directions, ranged, range_mask = read_reference()
    comparison, golden = geometry_audit(response, epochs, native, directions, ranged, range_mask, order=12)
    quadrature, disk_golden = quadrature_audit(response)
    historical = historical_solar_regression(response, data)
    golden.extend(disk_golden)
    # Exact new-Earth boundary and toggled Earthlight are separate portable cases.
    for include_earth, phase in ((True,180.0), (False,60.0)):
        g = {**golden[0]["geometry"], "earth_phase_deg":phase}
        golden.append(dict(name=f"Earth enabled={include_earth}, phase={phase:g}", kind="synthetic_boundary",
                           geometry=g, includeEarth=include_earth, order=12,
                           expected=response.light(g, include_earth)))
    files = ["illumination/calendar/validate.py", "illumination/calendar/disk.py", "illumination/calendar/geometry.py",
             "illumination/calendar/evaluator.mjs", "geography/lunar_ephemeris.py"]
    result = dict(schema=SCHEMA, evidence="Independent-reader numerical audit and propagation of sampled JPL direction/range residuals through a conditional atmosphere/source model.",
                  producer=dict(files={name:digest(ROOT/name) for name in files}),
                  inputs={str(p.relative_to(ROOT)):digest(p) for p in (transfer_path, REFERENCE, REFERENCE_INPUTS,
                                                                     HISTORICAL_SOLVED, HISTORICAL_CURVE)},
                  time_scale="TT; no UTC or future leap-second assumption enters this comparison.",
                  frame="MOON_ME, east-positive longitude; local topocentric geometry subtracts the lunar observer position.",
                  ephemeris_samples=len(records), geometric_range_samples=int(range_mask.sum()),
                  span_jd_tt=[float(epochs.min()),float(epochs.max())], sites=SITES,
                  conventions=dict(direction_only="Both ephemerides use identical native Earth/Moon and Sun/Moon ranges; only the directions differ. This also covers rows without independent ranges.",
                                   directions_and_geometric_ranges="Only the 2,507 epochs with both independently sampled JPL geometric ranges; directions and ranges are changed together.",
                                   earth_source="Phase/bright-limb orientation are recomputed from each topocentric geometry. Normalized Lambert spatial weights multiply the exported measured spectral phase response exactly once.",
                                   dilution="Solar inverse square from Sun/observer; Earth solid-angle ratio to the reference distance, multiplied by Sun/Earth inverse square.",
                                   comparator="JPL observer coordinates are apparent and light-time corrected; native vectors are instantaneous. Residuals retain that convention difference.",
                                   interpolation="Linear in native source elevation and empirical phase nodes, matching evaluator.mjs. Disk quadrature orders 12 and 24.",
                                   eclipses="Potential Sun/Earth angular-disk overlaps are excluded from error summaries; eclipse attenuation is not solved.",
                                   error_reading="Observed maxima at saved samples only; RMS describes this deliberately nonuniform sample, not a climatological average. Small source-relative errors can shift grazing brightness thresholds greatly. This audit does not assign a uniform timing accuracy.",
                                   boundary="No observational validation of the atmosphere or future Earth cloud/surface pattern; spatial Earth spectrum is uniform. Terrain, refraction, polarization and finite-distance variation within the atmospheric shell remain outside this response."),
                  ephemeris_to_light=comparison, disk_quadrature=quadrature,
                  historical_solar_regression=historical,
                  golden_light=golden,
                  golden_reading_rule="Feed geometry, includeEarth and order to lightAt. native_ephemeris cases also permit end-to-end geometry+light comparison at their TT Julian date and coordinates.")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    summary = {name:{k:row["max_abs_lux"] for k,row in group["components"].items() if k in ("solar","earth","total")}
               for name,group in comparison.items()}
    print(json.dumps(dict(output=str(output_path), bytes=output_path.stat().st_size,
                         ephemeris_cases=summary, golden_cases=len(golden)), indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transfer",type=Path,default=TRANSFER)
    parser.add_argument("--output",type=Path,default=OUTPUT)
    args=parser.parse_args()
    validate(args.transfer,args.output)


if __name__ == "__main__":
    main()
