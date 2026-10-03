"""The numbers behind the review of the wave studies (review.md, 2026-10-03).

Each quantity the review quotes is recomputed here from committed products and
from the 60-day atmospheric archive that drives the waves:

- the seas' sizes in the 28% atlas;
- the GCM's open-sea roughness against Charnock's law at lunar and Earth gravity;
- the drag coefficient of that law against SWAN's Wu drag at the strip winds;
- near-surface winds over the equatorial seas in the GCM and in CM1;
- the height of the GCM's lowest model level;
- the inverse-barometer tilt of each sea under the GCM's surface pressure;
- a first estimate of the monthly equilibrium tide in each sea;
- what the computed waves look like: wavelength, speed and steepness.

Run from the repository root with `python -m climate.waves.review_checks`.
The archive is research/runs/waves/cycle/atmosphere.npz on the research
drive; the script stops if it or a committed product is missing.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from climate.waves.model import ROOT, sha256
from climate.waves.weather import drag_coefficient, stress_equivalent_wind
from shared.constants import (ANOMALISTIC_MONTH_DAYS, DATA, DRACONIC_MONTH_DAYS, EARTH_GM,
                              EARTH_MOON_DISTANCE, MOON_EQUATOR_TO_ECLIPTIC_DEG, MOON_LOVE_H2,
                              MOON_LOVE_K2, MOON_ORBIT_ECCENTRICITY, MOON_ORBIT_INCLINATION_DEG,
                              MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY, SUN_GM,
                              SYNODIC_MONTH_DAYS, AU)

HERE = Path(__file__).resolve().parent
ATLAS_JSON = ROOT / "geography/results/atlas.json"
ATLAS_GRID = ROOT / "geography/products/atlas_28pct_4ppd.npz"
CM1_FORCING = HERE / "results/forcing.json"
COASTAL_HISTORY = HERE / "results/coastal_history.json"
CYCLE = HERE / "results/cycle.json"
WATER = ROOT / "shared/scenarios/waves.json"
ARCHIVE = ROOT / "research/runs/waves/cycle/atmosphere.npz"
OUTPUT = HERE / "results/review_checks.json"

GCM = DATA["model_closures"]["exoplasim_3_4_2"]
CHARNOCK, Z0_FLOOR, KAPPA = GCM["charnock_parameter"], GCM["sea_roughness_floor_m"], GCM["von_karman_constant"]
NEARSIDE, SPA, SMYTHII = 433, 1417, 874
HEADLAND = (93.4, 2.34)  # the eastern Smythii headland of the shore studies (lon, lat)
DAY = 86400.0


def charnock_drag(u10, gravity, alpha=CHARNOCK, floor=Z0_FLOOR, kappa=KAPPA, height=10.0):
    """Neutral drag coefficient at `height` for ExoPlaSim's open-sea law, z0 = max(alpha u*^2/g, floor)."""
    u10 = np.asarray(u10, dtype=float)
    ustar = 0.035 * u10
    for _ in range(200):
        z0 = np.maximum(alpha * ustar**2 / gravity, floor)
        ustar = kappa * u10 / np.log(height / z0)
    return (ustar / np.where(u10 > 0, u10, 1.0))**2


def neutral_wind(ustar, gravity, alpha=CHARNOCK, floor=Z0_FLOOR, kappa=KAPPA, height=10.0):
    """The neutral wind at `height` that the same law gives for a friction velocity."""
    z0 = np.maximum(alpha * np.asarray(ustar)**2 / gravity, floor)
    return np.asarray(ustar) / kappa * np.log(height / z0)


def deep_water(period_s, gravity):
    """Wavelength (m), phase speed and group speed (m/s) of a deep-water wave."""
    length = gravity * period_s**2 / (2 * math.pi)
    speed = gravity * period_s / (2 * math.pi)
    return length, speed, speed / 2


def earth_equivalent(length_m):
    """Period and speed of the deep-water Earth wave with the same wavelength."""
    return math.sqrt(2 * math.pi * length_m / STANDARD_GRAVITY), math.sqrt(STANDARD_GRAVITY * length_m / (2 * math.pi))


def unit_vectors(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)), axis=-1)


def first_tide_estimate(points_lat, points_lon, weights, days=365.25, step_days=0.125):
    """Monthly equilibrium tide (m) at points of one enclosed sea, its area mean removed.

    Earth's distance follows the Keplerian eccentricity; the sub-Earth point
    moves by the leading optical librations, 2e sin M in longitude and
    -(I + i) sin F in latitude. The Sun's sub-solar point circles the equator
    once a synodic month. Degree-2 potential, Love factor 1 + k2 - h2.
    """
    t = np.arange(0.0, days, step_days)
    mean_anomaly = 2 * math.pi * t / ANOMALISTIC_MONTH_DAYS
    latitude_argument = 2 * math.pi * t / DRACONIC_MONTH_DAYS
    e = MOON_ORBIT_ECCENTRICITY
    earth_lon = np.degrees(2 * e * np.sin(mean_anomaly))
    earth_lat = -(MOON_EQUATOR_TO_ECLIPTIC_DEG + MOON_ORBIT_INCLINATION_DEG) * np.sin(latitude_argument)
    distance_factor = (1 - e * np.cos(mean_anomaly))**-3  # (a/r)^3 to first order in the orbit
    sun_lon = 180.0 - 360.0 * t / SYNODIC_MONTH_DAYS
    love = 1 + MOON_LOVE_K2 - MOON_LOVE_H2
    earth = EARTH_GM * MOON_RADIUS**2 / (MOON_SURFACE_GRAVITY * EARTH_MOON_DISTANCE**3) * love
    sun = SUN_GM * MOON_RADIUS**2 / (MOON_SURFACE_GRAVITY * AU**3) * love
    points = unit_vectors(points_lat, points_lon)
    p2 = lambda c: 1.5 * c * c - 0.5
    height = (earth * distance_factor[:, None] * p2(unit_vectors(earth_lat, earth_lon) @ points.T)
              + sun * p2(unit_vectors(np.zeros_like(sun_lon), sun_lon) @ points.T))
    height -= (height * weights).sum(axis=1, keepdims=True) / weights.sum()
    return height, dict(earth_coefficient_m=earth, sun_coefficient_m=sun, love_factor=love,
                        samples=len(t), span_days=days, step_days=step_days)


def require(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing input: {path}. Restore it before running the review checks.")
    return path


def seas():
    atlas = json.loads(require(ATLAS_JSON).read_text())
    water = atlas["water_share"]
    rows = []
    for body in atlas["bodies"]:
        rows.append(dict(label=body["label"], names=body["water_names"][:4], area_km2=body["area_km2"],
                         share_of_moon=body["share"], share_of_water=body["share"] / water,
                         mean_depth_m=body["mean_depth_m"], max_depth_m=body["max_depth_m"]))
    return dict(water_share_of_moon=water, small_lakes_share_of_moon=atlas["small_lakes_share"], bodies=rows)


def archive_fields(path):
    with np.load(require(path), allow_pickle=False) as d:
        fields = {k: d[k] for k in d.files if k != "metadata"}
        meta = json.loads(d["metadata"].item())
    return fields, meta


def gcm_roughness(a, meta):
    sea = a["land_fraction"][0] < 0.5
    ustar = np.sqrt(np.hypot(a["stress_eastward_Pa"], a["stress_northward_Pa"]) / a["stress_closure_density_kg_m3"])
    z0, us = a["roughness_m"][:, sea], ustar[:, sea]
    above = us > np.sqrt(Z0_FLOOR * meta["gravity_m_s2"] / CHARNOCK) * 2  # well above the roughness floor
    def ratio(g):
        r = z0[above] / (CHARNOCK * us[above]**2 / g)
        return dict(median=float(np.median(r)), p5=float(np.percentile(r, 5)), p95=float(np.percentile(r, 95)))
    return dict(samples=int(above.sum()), friction_velocity_threshold_m_s=float(np.sqrt(Z0_FLOOR * meta["gravity_m_s2"] / CHARNOCK) * 2),
                ratio_to_charnock_at_lunar_gravity=ratio(meta["gravity_m_s2"]),
                ratio_to_charnock_at_earth_gravity=ratio(STANDARD_GRAVITY),
                charnock_parameter=CHARNOCK, roughness_floor_m=Z0_FLOOR)


def drag_table(winds):
    rows = []
    for u in winds:
        wu = float(drag_coefficient(u))
        lunar, earth = float(charnock_drag(u, MOON_SURFACE_GRAVITY)), float(charnock_drag(u, STANDARD_GRAVITY))
        rows.append(dict(wind_10m_m_s=u, wu=wu, charnock_lunar=lunar, charnock_earth=earth,
                         lunar_charnock_to_wu=lunar / wu))
    return rows


def equatorial_winds(a, cm1, density):
    lat = a["latitude_deg"]
    rows = np.abs(lat) < 3.0
    sea = (a["land_fraction"][0] < 0.5) & rows[:, None]
    tx, ty = a["stress_eastward_Pa"][:, sea], a["stress_northward_Pa"][:, sea]
    ue, ve = stress_equivalent_wind(tx, ty, density)
    ustar = np.sqrt(np.hypot(tx, ty) / a["stress_closure_density_kg_m3"][:, sea])
    lowest = np.hypot(a["eastward_lowest_m_s"][:, sea], a["northward_lowest_m_s"][:, sea])
    q = (50, 90, 99)
    pct = lambda x: dict(zip(("p50", "p90", "p99"), np.percentile(x, q).round(3).tolist()))
    stats = cm1["water_wind_statistics"]
    cm1_row = dict(p50=round(stats["median"], 3), p90=round(stats["p90"], 3), p99=round(stats["p99"], 3))
    neutral = neutral_wind(ustar, MOON_SURFACE_GRAVITY)
    return dict(gcm_rows_deg=lat[rows].round(3).tolist(), gcm_sea_cells=int(sea.sum()), samples=int(tx.size),
                gcm_stress_equivalent_wu=pct(np.hypot(ue, ve)), gcm_neutral_10m_charnock=pct(neutral),
                gcm_lowest_level=pct(lowest), cm1_equatorial_ring_10m=cm1_row,
                gcm_neutral_to_cm1={k: round(pct(neutral)[k] / cm1_row[k], 3) for k in cm1_row},
                cm1_span_days=cm1["span_days"])


def lowest_level(a, meta):
    sea = a["land_fraction"][0] < 0.5
    temperature = float(a["surface_temperature_K"][:, sea].mean())
    scale_height = meta["gas_constant_J_kg_K"] * temperature / meta["gravity_m_s2"]
    return dict(sigma=meta["sigma_lowest"], mean_sea_surface_temperature_K=temperature,
                scale_height_m=scale_height, height_m=-scale_height * math.log(meta["sigma_lowest"]))


def sea_labels_at(lat_deg, lon_deg, grid):
    """Atlas water label at each (lat, lon), from the nearest 4-pixel-per-degree node."""
    rows = np.clip(np.round((89.875 - lat_deg) * 4).astype(int), 0, grid.shape[0] - 1)
    cols = np.round((np.mod(lon_deg, 360.0) - 0.125) * 4).astype(int) % grid.shape[1]
    return grid[rows, cols]


def inverse_barometer(a, labels, water_density):
    lat, lon = a["latitude_deg"], a["longitude_deg"]
    lon2, lat2 = np.meshgrid(lon, lat)
    sea = a["land_fraction"][0] < 0.5
    cell_label = sea_labels_at(lat2, lon2, labels)
    out = {}
    for body in (NEARSIDE, SPA, SMYTHII):
        cells = sea & (cell_label == body)
        if cells.sum() < 2:
            out[str(body)] = dict(cells=int(cells.sum()))
            continue
        w = np.cos(np.radians(lat2[cells]))
        p = a["surface_pressure_Pa"][:, cells]
        anomaly = p - (p * w).sum(axis=1, keepdims=True) / w.sum()
        level = -anomaly / (water_density * MOON_SURFACE_GRAVITY)
        span = level.max(axis=0) - level.min(axis=0)
        out[str(body)] = dict(cells=int(cells.sum()), range_median_m=float(np.median(span)),
                              range_max_m=float(span.max()))
    return out


def tide_first_estimates(grid):
    labels, lat, lon = grid["water_label"], grid["lat_deg"], grid["lon_deg"]
    out = {}
    for body in (NEARSIDE, SPA, SMYTHII):
        rows, cols = np.where(labels[::4, ::4] == body)  # every fourth node: a 1-degree sample
        plat, plon = lat[::4][rows], lon[::4][cols]
        height, settings = first_tide_estimate(plat, plon, np.cos(np.radians(plat)))
        span = height.max(axis=0) - height.min(axis=0)
        out[str(body)] = dict(nodes=int(len(rows)), range_median_m=float(np.median(span)),
                              range_p95_m=float(np.percentile(span, 95)), range_max_m=float(span.max()))
        if body == SMYTHII:
            # The headland is one more point of the same sea; its level is relative to the sea's mean.
            pts_lat, pts_lon = np.append(plat, HEADLAND[1]), np.append(plon, HEADLAND[0])
            weights = np.append(np.cos(np.radians(plat)), 0.0)
            series, _ = first_tide_estimate(pts_lat, pts_lon, weights)
            out["headland"] = dict(longitude_deg=HEADLAND[0], latitude_deg=HEADLAND[1],
                                   range_m=float(series[:, -1].max() - series[:, -1].min()))
    out["settings"] = settings
    return out


def smythii_development(atlas_seas, density, water_density):
    """Smythii's strongest basin-mean forcing against the fully developed limit, and its wind setup."""
    cycle = json.loads(require(CYCLE).read_text())
    hours = np.array(cycle["forcing"]["monthly_time_days"]) * 24
    stress = np.array(cycle["forcing"]["monthly_basin_mean_stress_Pa"])
    start, end = cycle["clock"]["report_interval_hours"]
    window = (hours >= start) & (hours <= end)
    peak = float(stress[window].max())
    wind = float(stress_equivalent_wind(np.array([peak]), np.array([0.0]), density)[0][0])
    body = next(b for b in atlas_seas["bodies"] if b["label"] == SMYTHII)
    length = 2 * math.sqrt(body["area_km2"] * 1e6 / math.pi)  # diameter of a disc of the sea's area
    return dict(peak_basin_mean_stress_Pa=peak, stress_equivalent_wind_m_s=wind,
                fully_developed_hs_m=0.24 * wind**2 / MOON_SURFACE_GRAVITY,
                fully_developed_rule="g Hs / U10^2 = 0.24, the Pierson-Moskowitz limit",
                computed_maximum_hs_m=cycle["statistics"]["max_hs_m"],
                wind_setup_m=peak * length / (water_density * MOON_SURFACE_GRAVITY * body["mean_depth_m"]),
                wind_setup_rule="tau L / (rho g h) across a disc-equivalent sea of the mean depth")


def appearance():
    coast = json.loads(require(COASTAL_HISTORY).read_text())["stations"]["west_face"]
    cycle = json.loads(require(CYCLE).read_text())["statistics"]["peak_record"]
    cases = [("western face, strongest arrival", coast["maximum_hs_m"], coast["mean_period_at_peak_s"]),
             ("basin peak, mean period", cycle[3], cycle[5]),
             ("basin peak, peak period", cycle[3], cycle[4])]
    rows = []
    for name, hs, period in cases:
        length, speed, group = deep_water(period, MOON_SURFACE_GRAVITY)
        earth_period, earth_speed = earth_equivalent(length)
        rows.append(dict(case=name, hs_m=hs, period_s=period, wavelength_m=length, phase_speed_m_s=speed,
                         group_speed_m_s=group, steepness_hs_over_wavelength=hs / length,
                         earth_period_same_wavelength_s=earth_period, earth_speed_same_wavelength_m_s=earth_speed))
    return dict(cases=rows, speed_ratio_same_wavelength=math.sqrt(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY),
                fall_time_ratio_same_height=math.sqrt(STANDARD_GRAVITY / MOON_SURFACE_GRAVITY))


def build(archive=ARCHIVE):
    cm1 = json.loads(require(CM1_FORCING).read_text())
    water_density = json.loads(require(WATER).read_text())["water_density_kg_m3"]
    density = cm1["smythii_marginis_equatorial_density_statistics"]["mean"]  # the coupled SWAN build's air
    a, meta = archive_fields(archive)
    with np.load(require(ATLAS_GRID), allow_pickle=False) as g:
        grid = {k: g[k] for k in ("water_label", "lat_deg", "lon_deg")}
    stats = cm1["water_wind_statistics"]
    atlas_seas = seas()
    producer = [Path(__file__), ROOT / "climate/waves/weather.py", ROOT / "shared/constants.json"]
    return dict(
        schema="terluna.climate.wave-review-checks/1",
        producer=dict(domain="climate", files={str(p.relative_to(ROOT)): sha256(p) for p in producer}),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in
                (ATLAS_JSON, ATLAS_GRID, CM1_FORCING, COASTAL_HISTORY, CYCLE, WATER, archive)},
        evidence=("Recomputations from committed wave, climate and geography products and the 60-day "
                  "GCM surface archive. The tide is a first estimate: Keplerian distance, the two leading "
                  "optical librations and a circular solar term, in equilibrium with each sea's mean removed; "
                  "the geography tide product supersedes it."),
        reading_rule=("Percentiles pool every open-sea cell and archive time named in each section. Ranges are "
                      "maximum minus minimum over the sampled span at each node, summarised over the sea's nodes."),
        seas=atlas_seas,
        smythii_development=smythii_development(atlas_seas, density, water_density),
        gcm_sea_roughness=gcm_roughness(a, meta),
        drag_against_wu=drag_table([stats["median"], stats["p90"], stats["p99"], 7.0, 10.0]),
        equatorial_sea_winds=equatorial_winds(a, cm1, density),
        gcm_lowest_level=lowest_level(a, meta),
        inverse_barometer=dict(water_density_kg_m3=water_density,
                               seas=inverse_barometer(a, grid["water_label"], water_density)),
        monthly_tide_first_estimate=tide_first_estimates(grid),
        appearance=appearance(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    record = build(args.archive)
    args.output.write_text(json.dumps(record, indent=1) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
