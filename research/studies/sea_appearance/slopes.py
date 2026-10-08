"""The sea's slopes at the six coasts of the lighting calendar, hour by hour through the shore month.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.sea_appearance.slopes

How water reflects the sky, the Sun and the Earth follows from the distribution of its surface slopes. Two parts
make them up here:

- the waves the wave runs resolve, down to about 1 m (0.497 Hz at lunar gravity): their slope covariance from
  the directional spectra, at the nearside shores from the nearside run's restart files at the end of each
  48-hour segment (exact spectra, illumination/water_surface/hotfile.py), at the Smythii headland's east face from
  the shore history's hourly spectra as the optical comfort study read them;
- the shorter waves: the short-wave regime of Elfouhaily et al. (1997) at lunar gravity, above the wave runs'
  highest frequency (illumination/water_surface/short_waves.py), driven hour by hour by the friction velocity of
  the GCM's surface stress at each coast, interpolated as the wave coupling interpolates it, and aligned with
  that stress.

The South Pole-Aitken sea has no wave run, so only its short waves are known. Earth's reference is Cox and
Munk's (1954) clean sea under the same 10 m wind, which comes from the GCM's own neutral surface layer
(Charnock's roughness at lunar gravity).
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from climate.waves.nearside import periodic_weights
from climate.waves.weather import apply_weights
from illumination.water_surface.hotfile import read_hotfile
from illumination.water_surface.model import slope_moments, wavenumber
from illumination.water_surface.short_waves import capillary_scales, neutral_wind, short_wave_slopes
from shared.constants import (EXOPLASIM_CHARNOCK, EXOPLASIM_SEA_ROUGHNESS_FLOOR, EXOPLASIM_VON_KARMAN,
                              MOON_SURFACE_GRAVITY, STANDARD_GRAVITY)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CALENDAR = HERE / "results/lighting_calendar.json"
RESULT = HERE / "results/sea_slopes.json"
SERIES = HERE / "results/sea_slopes.npz"
ATMOSPHERE = ROOT / "research/runs/waves/cycle/atmosphere.npz"
NEARSIDE = ROOT / "research/runs/waves/nearside/dt300"
CUBE = NEARSIDE / "sea.npz"
SMYTHII = ROOT / "research/runs/optical_comfort/water_slope_history.json"
SMYTHII_RECORD = ROOT / "research/studies/optical_comfort/results/water_slopes.json"
SMYTHII_TABLE = ROOT / "research/runs/waves/shore_history/coast_s4_dt900_0_1419/sites.tbl"
SMYTHII_STATION = "east_face"
RUN_START = datetime(2000, 1, 1)
GRAVITY = MOON_SURFACE_GRAVITY


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def quantiles(values):
    values = np.asarray(values, float)
    if not len(values):
        return None
    return dict(zip(("p10", "p50", "p90", "maximum"), np.round(np.quantile(values, [.1, .5, .9, 1]), 6).tolist()))


def refine(f, e):
    """Halve each frequency interval, interpolating the variance linearly (as the optical comfort study does)."""
    x = np.empty(2 * len(f) - 1); x[::2] = f; x[1::2] = (f[1:] + f[:-1]) / 2
    y = np.empty((len(x), e.shape[1])); y[::2] = e; y[1::2] = (e[1:] + e[:-1]) / 2
    return x, y


def resolved_covariance(frequency, direction, variance, depth):
    f, e = frequency, variance
    for _ in range(2):
        f, e = refine(f, e)
    return slope_moments(f, direction, e, depth, GRAVITY)["covariance"]


def rotated(along, across, toward_deg):
    """East-north covariance of slopes with variances along and across a direction (degrees from east)."""
    a = np.radians(toward_deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[along * c * c + across * s * s, (along - across) * c * s],
                     [(along - across) * c * s, along * s * s + across * c * c]])


def cox_munk_clean(u10):
    """Cox and Munk's (1954) clean-sea mean square slope, 0.003 + 5.12e-3 W, with W the wind at 12.5 m on Earth."""
    u10 = np.asarray(u10, float)
    u_star = 0.04 * u10
    for _ in range(30):
        roughness = np.maximum(EXOPLASIM_CHARNOCK * u_star ** 2 / STANDARD_GRAVITY, EXOPLASIM_SEA_ROUGHNESS_FLOOR)
        u_star = EXOPLASIM_VON_KARMAN * u10 / np.log(10 / roughness)
    return 0.003 + 5.12e-3 * u10 * np.log(12.5 / roughness) / np.log(10 / roughness)


def friction(hours, coasts):
    """The GCM's surface stress at each coast, as the wave coupling interpolates it: bilinear over the cells that
    stay open water through the record, linear in time between its three-hourly snapshots."""
    with np.load(ATMOSPHERE, allow_pickle=False) as d:
        times = d["time_s"] / 3600
        lat, lon = d["latitude_deg"][::-1], d["longitude_deg"]
        water = np.all((d["land_fraction"][:, ::-1] < .5) & (d["sea_ice_fraction"][:, ::-1] < .01), axis=0)
        fields = {k: d[k][:, ::-1] for k in ("stress_eastward_Pa", "stress_northward_Pa",
                                               "stress_closure_density_kg_m3")}
    if not np.allclose(np.diff(times), 3) or times[0] != 0:
        raise ValueError("A complete three-hourly record is required")
    out = {}
    for short, _, lon_c, lat_c in coasts:
        indices, weights, fallback = periodic_weights(lon, lat, water, np.array([[lon_c]]), np.array([[lat_c]]))
        tx, ty, rho = (np.interp(hours, times, apply_weights(fields[k], indices, weights)[:, 0, 0])
                       for k in ("stress_eastward_Pa", "stress_northward_Pa", "stress_closure_density_kg_m3"))
        u_star = np.sqrt(np.hypot(tx, ty) / rho)
        out[short] = dict(u_star=u_star, toward_deg=np.degrees(np.arctan2(ty, tx)) % 360,
                          u10=neutral_wind(u_star, GRAVITY, EXOPLASIM_CHARNOCK, EXOPLASIM_SEA_ROUGHNESS_FLOOR,
                                           EXOPLASIM_VON_KARMAN),
                          nearest_open_water_fallback=bool(fallback.any()))
    return out


def nearside_resolved(coasts, first, last):
    """Slope covariance of the resolved waves at the nearside shores, from every restart file in the month."""
    with np.load(CUBE, allow_pickle=False) as z:
        nodes = np.c_[z["lon"], z["lat"]]
        node = {short: int(np.argmin(np.hypot(*(nodes - (lon, lat)).T))) for short, _, lon, lat in coasts}
        depth = {short: float(z["depth"][i]) for short, i in node.items()}
        cube = dict(hours=z["hours"], tp=z["tp"][:, list(node.values())], node=node)
    rows, hashes = {short: [] for short, *_ in coasts}, {}
    for segment in sorted(NEARSIDE.glob("s*_*")):
        end = int(segment.name.split("_")[1])
        if not first <= end <= last:
            continue
        recorded = json.loads((segment / "run.json").read_text())["output_sha256"]["end.hot"]
        if sha256(segment / "end.hot") != recorded:
            raise ValueError(f"Restart file differs from its run record: {segment.name}")
        hashes[str((segment / "end.hot").relative_to(ROOT))] = recorded
        spectra = read_hotfile(segment / "end.hot", nodes=[(lon, lat) for _, _, lon, lat in coasts])
        if spectra["time"] != RUN_START + timedelta(hours=end):
            raise ValueError(f"Restart time differs from the segment's end: {segment.name}")
        for i, (short, _, lon, lat) in enumerate(coasts):
            if not np.allclose(spectra["lonlat"][i], (lon, lat), atol=1e-3):
                raise ValueError(f"No grid point at the shore {short}")
            rows[short].append((end, resolved_covariance(spectra["frequency"], spectra["direction"],
                                                         spectra["variance"][i], depth[short])))
    return rows, depth, cube, hashes, float(spectra["frequency"][-1])


def smythii_resolved(first, last):
    record = json.loads(SMYTHII_RECORD.read_text())["producer"]["history"]
    if sha256(SMYTHII) != record["sha256"]:
        raise ValueError("The Smythii slope history differs from the optical comfort study's record")
    history = json.loads(SMYTHII.read_text())
    rows = [(r["hour"], np.array(r["covariance_east_north"]), r["depth_m"]) for r in history["rows"]
            if r["site"] == SMYTHII_STATION and first <= r["hour"] <= last]
    table = np.loadtxt(SMYTHII_TABLE).reshape(-1, 4, 12)
    station = json.loads((SMYTHII_TABLE.parent / "product.json").read_text())["sites"]
    column = list(station).index(SMYTHII_STATION)
    peak = dict(zip(np.rint(table[:, column, 0] / 3600).astype(int), table[:, column, 4]))
    return rows, peak, station[SMYTHII_STATION], record["sha256"]


def peak_speed(period, depth):
    """Phase speed of the resolved sea's peak; infinite where no peak is known."""
    if not np.isfinite(period) or period <= 0:
        return np.inf
    k = wavenumber(np.array([1 / period]), depth, GRAVITY)[0]
    return float(2 * np.pi / period / k)


def build():
    calendar = json.loads(CALENDAR.read_text())
    hours = np.array(calendar["time_hours"], float)
    coasts = [(s["short"], s["name"], *s["lon_lat"]) for s in calendar["sites"]]
    first, last = hours[0], hours[-1]
    wind = friction(hours, coasts)
    nearside = coasts[:4]
    near_rows, near_depth, cube, hot_hashes, top_frequency = nearside_resolved(nearside, first, last)
    smythii_rows, smythii_peak, smythii_site, smythii_hash = smythii_resolved(first, last)
    k_m, c_m = capillary_scales(GRAVITY)
    deep_cutoff = (2 * np.pi * top_frequency) ** 2 / GRAVITY
    series, summary = dict(hours=hours), []
    for short, name, lon, lat in coasts:
        w = wind[short]
        if short in near_depth:
            depth = near_depth[short]
            column = list(cube["node"]).index(short)
            periods = np.interp(hours, cube["hours"], cube["tp"][:, column])
            resolved = near_rows[short]
            source = "nearside run restart files, every 48 hours"
        elif short == "Smythii headland":
            depth = smythii_rows[0][2]
            periods = np.array([smythii_peak.get(int(h), np.nan) for h in hours])
            resolved = [(h, c) for h, c, _ in smythii_rows]
            source = f"Smythii shore history, hourly, station {SMYTHII_STATION} at {smythii_site}"
        else:
            depth, periods, resolved, source = np.inf, np.full(len(hours), np.nan), [], None
        cutoff = deep_cutoff if not np.isfinite(depth) else float(wavenumber(np.array([top_frequency]), depth, GRAVITY)[0])
        short_cov, level = np.zeros((len(hours), 2, 2)), np.zeros(len(hours))
        for i, (u_star, toward, period) in enumerate(zip(w["u_star"], w["toward_deg"], periods)):
            s = short_wave_slopes(u_star, GRAVITY, cutoff,
                                  peak_speed=peak_speed(period, depth if np.isfinite(depth) else 1e4))
            short_cov[i], level[i] = rotated(s["along"], s["across"], toward), s["level"]
        earth = cox_munk_clean(w["u10"])
        short_mss = np.trace(short_cov, axis1=1, axis2=2)
        resolved_hours = np.array([h for h, _ in resolved], float)
        resolved_cov = np.array([c for _, c in resolved]).reshape(-1, 2, 2)
        at = np.searchsorted(hours, resolved_hours)
        total_cov = resolved_cov + short_cov[at] if len(at) else resolved_cov
        total_mss = np.trace(total_cov, axis1=1, axis2=2)
        key = short.replace(" ", "_").lower()
        series.update({f"{key}/u_star": w["u_star"], f"{key}/u10": w["u10"], f"{key}/toward_deg": w["toward_deg"],
                       f"{key}/short_covariance": short_cov.reshape(-1, 4), f"{key}/short_level": level,
                       f"{key}/earth_cox_munk_mss": earth, f"{key}/resolved_hours": resolved_hours,
                       f"{key}/resolved_covariance": resolved_cov.reshape(-1, 4)})
        summary.append(dict(
            name=name, short=short, lon_lat=[lon, lat], depth_m=None if not np.isfinite(depth) else round(depth, 2),
            resolved_waves=source, cutoff_wavenumber_rad_m=round(cutoff, 4),
            cutoff_wavelength_m=round(2 * np.pi / cutoff, 3),
            gcm_nearest_open_water_fallback=w["nearest_open_water_fallback"],
            friction_velocity_m_s=quantiles(w["u_star"]), wind_10m_m_s=quantiles(w["u10"]),
            hours_without_short_waves=int(np.sum(level == 0)),
            short_wave_mss=quantiles(short_mss),
            resolved_mss=quantiles(np.trace(resolved_cov, axis1=1, axis2=2)), resolved_samples=len(resolved),
            total_mss=quantiles(total_mss),
            earth_cox_munk_mss_same_wind=quantiles(earth),
            total_over_earth_same_wind=quantiles(total_mss / earth[at]) if len(at) else None))
    SERIES.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(SERIES, **{k: np.asarray(v, np.float32 if k != "hours" else float) for k, v in series.items()})
    own = [Path(__file__), ROOT / "illumination/water_surface/short_waves.py",
           ROOT / "illumination/water_surface/hotfile.py", ROOT / "illumination/water_surface/model.py"]
    # The wave coupling's interpolation helpers; the constants they read serve other paths of those modules.
    files = own + [ROOT / "climate/waves/nearside.py", ROOT / "climate/waves/weather.py"]
    return dict(
        schema="terluna.research.sea-slopes/1",
        evidence=("Surface slope statistics at six coasts through the shore month: the resolved waves' slope "
                  "covariance from the wave runs' directional spectra, and the short waves from the unified "
                  "spectrum of Elfouhaily et al. (1997) at lunar gravity, driven by the GCM's friction velocity. "
                  "The short-wave law is an Earth fit; it carries over if the lunar seas' short waves follow the "
                  "same balance of wind and capillarity, with u*/c_m inside the fitted range."),
        reading_rule=("Covariances are of the surface gradient, east and north, dimensionless. Mean square slope "
                      "is their trace. Short waves start at the wave runs' highest frequency, 0.497 Hz, and run "
                      "through the capillary roll-off; they align with the GCM's stress. Resolved and short-wave "
                      "slopes add at the hours with spectra. Hours count from the first atmospheric snapshot, as in "
                      "the lighting calendar. The GCM's three-hourly stress is a 170 km mean: gusts, sea breezes and "
                      "surface films are absent, so the light-wind hours carry most uncertainty."),
        producer=dict(lane="research", runner=str(Path(__file__).relative_to(ROOT)),
                      files={str(p.relative_to(ROOT)): sha256(p) for p in files}, constants=constants_used(own)),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (CALENDAR, ATMOSPHERE, CUBE, SMYTHII_TABLE)} | {
            str(SMYTHII.relative_to(ROOT)): smythii_hash} | hot_hashes,
        series=dict(file=str(SERIES.relative_to(ROOT)), sha256=sha256(SERIES)),
        gravity_m_s2=GRAVITY, capillary_wavenumber_rad_m=round(k_m, 3), capillary_wavelength_m=round(2 * np.pi / k_m, 5),
        minimum_phase_speed_m_s=round(c_m, 5), short_wave_threshold_u_star_m_s=round(c_m / np.e, 5),
        earth_capillary_wavelength_m=round(2 * np.pi / capillary_scales(STANDARD_GRAVITY)[0], 5),
        earth_minimum_phase_speed_m_s=round(capillary_scales(STANDARD_GRAVITY)[1], 5),
        time_hours=[float(first), float(last)], coasts=summary)


def main():
    product = build()
    RESULT.write_text(json.dumps(product, indent=1, allow_nan=False) + "\n")
    for c in product["coasts"]:
        print(c["short"], "| u*", c["friction_velocity_m_s"], "| no short waves", c["hours_without_short_waves"],
              "| short", c["short_wave_mss"], "| resolved", c["resolved_mss"], "| total", c["total_mss"],
              "| Earth", c["earth_cox_munk_mss_same_wind"])


if __name__ == "__main__":
    main()
