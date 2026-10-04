"""A lighting calendar of the seas: the Sun, the Earth and the light at six coasts through a month.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.sea_appearance.lighting

Hour by hour through the shore month of the wave studies (climate/waves/results/shore_month.json, whose tide month
matches the GCM's Sun), at the four shores of that month, the eastern Smythii headland of the shore studies and a
South Pole–Aitken coast by Mare Ingenii:

- the Sun's and the Earth's elevation and azimuth, from geography/lunar_ephemeris.py, checked against JPL Horizons
  within 0.05 degrees; the Earth's direction is topocentric (its parallax reaches 0.26 degrees);
- the Earth's phase angle and lit fraction, and its earthlight above the air by wavelength: Glenar et al.'s
  spectrum of the whole Earth scaled to Robinson et al.'s (2025) observed visual phase curve, geometric albedo
  0.242 (illumination/earthlight/model.py), lit by the unshielded sunlight of the Earth control;
- the clear-sky light on the ground from the Sun: the solved spherical sky's direct beam, interpolated between its
  stored solar elevations, and its diffuse light, read from the cached scattering solution, which spans the Sun from
  90 degrees above the horizon to 90 below;
- the earthlight on the ground through the same air, channel by channel: the solved sky's direct beam and
  diffuse light for a source at the Earth's elevation, weighted by the earthlight's spectrum;
- starlight and airglow as a placeholder of 0.001 lux, a clear moonless night on Earth;
- which source dominates, and the mode of vision for an 18% grey surface against CIE 191:2010's mesopic range,
  0.005 to 5 cd/m2.

It also maps, for every water cell at 1-degree spacing, how much of 2026-2045 the Earth spends above the horizon.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from climate.waves.cycle import CYCLE_HOURS
from geography import lunar_ephemeris as ephemeris
from illumination.earthlight import model as earthlight
from illumination.earthlight.fetch_inputs import INPUTS as EARTHLIGHT_INPUTS
from shared.constants import AU, EARTH_MOON_DISTANCE, MOON_RADIUS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE / "results/lighting_calendar.json"
MAP = HERE / "results/earth_over_seas.npz"
SHORE_MONTH = ROOT / "climate/waves/results/shore_month.json"
SKY = ROOT / "illumination/sky/results/solved_sky.json"
SPHERICAL = ROOT / "research/runs/optical_comfort/spherical"
ATLAS = ROOT / "geography/products/atlas_28pct_4ppd.npz"
ATLAS_JSON = ROOT / "geography/results/atlas.json"
FLOOR_LUX = 1e-3                 # placeholder: starlight and airglow of a clear moonless night on Earth
GREY = 0.18                      # reference reflectance for the mode of vision
MESOPIC_CD_M2 = (0.005, 5.0)     # CIE 191:2010
JD_START, YEARS = 2461041.5, 19.0
# The Smythii headland of the shore studies (geography/tides.py) and the South Pole–Aitken sea's coast nearest
# Mare Ingenii's centre (found in the atlas by nearest_coast); the four shores come from the shore month.
HEADLAND = ("Eastern Smythii headland", "Smythii headland", 93.4, 2.34)
INGENII = (163.2, -33.7)
SPA_BODY = 1417


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit(lon_deg, lat_deg):
    lon, lat = np.radians(lon_deg), np.radians(lat_deg)
    return np.stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)), axis=-1)


def local_frame(lon_deg, lat_deg):
    """East, north and up at a site, in the Moon's body frame (x toward 0 E 0 N, z north)."""
    lon, lat = math.radians(lon_deg), math.radians(lat_deg)
    east = np.array([-math.sin(lon), math.cos(lon), 0.0])
    north = np.array([-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat)])
    return east, north, unit(lon_deg, lat_deg)


def horizontal(vectors, frame):
    """Elevation and azimuth (degrees, azimuth clockwise from north) of direction vectors (..., 3)."""
    east, north, up = frame
    v = vectors / np.linalg.norm(vectors, axis=-1, keepdims=True)
    return np.degrees(np.arcsin(np.clip(v @ up, -1, 1))), np.degrees(np.arctan2(v @ east, v @ north)) % 360


def geometry(jd, lon_deg, lat_deg):
    """The Sun and the Earth over a site, the Earth's phase and the sunlight that falls on it."""
    g = ephemeris.earth_and_sun(jd)
    frame = local_frame(lon_deg, lat_deg)
    earth = g["earth"] * g["earth_distance_m"][..., None] - MOON_RADIUS * frame[2]
    sun_el, sun_az = horizontal(g["sun"], frame)
    earth_el, earth_az = horizontal(earth, frame)
    to_sun = g["sun"] * g["sun_distance_m"][..., None] - g["earth"] * g["earth_distance_m"][..., None]
    to_moon = -g["earth"]
    cos_alpha = np.sum(to_sun * to_moon, axis=-1) / np.linalg.norm(to_sun, axis=-1)
    alpha = np.arccos(np.clip(cos_alpha, -1, 1))
    return dict(sun_elevation=sun_el, sun_azimuth=sun_az, earth_elevation=earth_el, earth_azimuth=earth_az,
                earth_phase_angle=np.degrees(alpha), earth_lit_fraction=(1 + np.cos(alpha)) / 2,
                earth_distance_m=g["earth_distance_m"],
                earth_sunlight_factor=(AU / np.linalg.norm(to_sun, axis=-1)) ** 2)


class Sky:
    """Clear-sky light on the ground from a source at any elevation, from the solved spherical sky."""

    def __init__(self, sky=SKY, directory=SPHERICAL):
        product = json.loads(Path(sky).read_text())
        moon = product["worlds"]["moon_1.2atm"]
        path = Path(directory) / "moon_1.2atm_standard.npz"
        with np.load(path, allow_pickle=False) as z:
            self.angles = np.degrees(z["a"])
            down = z["moments"][0, :, :, 5] @ z["xyz"][:, 1]
            last = z["last_order"][0, :, :, 5] @ z["xyz"][:, 1]
            self.above_air = float(z["xyz"][:, 1].sum())
            # Per channel and unit source: diffuse light and the direct beam on the ground, by source elevation.
            self.channel_y, self.channel_band, self.channel_energy = z["xyz"][:, 1], z["band"], z["energy"]
            self.diffuse_channels = z["moments"][0, :, :, 5]
            self.beam_angles, self.direct_channels = np.degrees(z["sa"]), z["beam"][0, :, :, 1]
        self.diffuse_table = down
        self.last_order_fraction = float(np.max(last / np.maximum(down, 1e-300)))
        rows = sorted(moon["samples"], key=lambda r: r["sun_deg"])
        self.sample_deg = np.array([r["sun_deg"] for r in rows])
        self.direct_normal = np.array([r["direct_normal_lux"] for r in rows])
        with np.load(Path(directory) / "earth_control_standard.npz", allow_pickle=False) as z:
            self.unfiltered_above_air = float(z["xyz"][:, 1].sum())
            self.unshielded_band, self.unshielded_energy = z["band"], z["energy"]
        stored = np.array([r["moment_diffuse_horizontal_lux"] for r in rows])
        self.readout_check = float(np.max(np.abs(self.diffuse(self.sample_deg) / stored - 1)))
        self.inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (Path(sky), path,
                       Path(directory) / "earth_control_standard.npz")}

    def diffuse(self, elevation):
        return np.interp(elevation, self.angles, self.diffuse_table)

    def direct(self, elevation):
        """Direct beam on a horizontal surface, log-interpolated between the stored solar elevations."""
        e = np.asarray(elevation, float)
        lit = self.direct_normal > 0
        normal = np.exp(np.interp(e, self.sample_deg[lit], np.log(self.direct_normal[lit])))
        return np.where(e > 0, normal * np.sin(np.radians(np.maximum(e, 0))), 0.0)

    def ground(self, elevation):
        return self.direct(elevation) + self.diffuse(elevation)

    def ground_channels(self, elevation):
        """Light on the ground per channel and unit source, for sources at the given elevations (n, channels)."""
        e = np.atleast_1d(np.asarray(elevation, float))
        diffuse = np.stack([np.interp(e, self.angles, self.diffuse_channels[:, c])
                            for c in range(self.diffuse_channels.shape[1])], -1)
        direct = np.stack([np.interp(e, self.beam_angles, self.direct_channels[:, c])
                           for c in range(self.direct_channels.shape[1])], -1)
        return diffuse + direct


class Earthlight:
    """The Earth's light above the air in each of the sky's channels, per unit of that channel's shielded sunlight.

    Bands share a weight across their channels: within a 10 nm band the earthlight and the sunlight keep the same
    spectral shape, and the channels split a band by absorption strength.
    """

    def __init__(self, sky, step_deg=0.5):
        bands = np.unique(sky.channel_band, axis=0)
        if not np.array_equal(bands, np.unique(sky.unshielded_band, axis=0)):
            raise ValueError("The shielded and unshielded skies use different bands")
        shielded = np.array([sky.channel_energy[(sky.channel_band == b).all(1)].sum() for b in bands])
        unshielded = np.array([sky.unshielded_energy[(sky.unshielded_band == b).all(1)].sum() for b in bands])
        self.phase = np.arange(0.0, 180.0 + step_deg / 2, step_deg)
        table = np.array([earthlight.calibrated_irradiance(bands, unshielded, a, EARTH_MOON_DISTANCE)
                          for a in self.phase])
        self.per_unit = table / shielded
        self.channel_index = np.array([int(np.flatnonzero((bands == b).all(1))[0]) for b in sky.channel_band])
        self.reference_solid_angle = earthlight.solid_angle(EARTH_MOON_DISTANCE)

    def weights(self, phase_deg, distance_m, sunlight_factor):
        """Per-channel weights (n, channels) for the Earth at the given phase angles and distances."""
        per_band = np.stack([np.interp(phase_deg, self.phase, self.per_unit[:, b])
                             for b in range(self.per_unit.shape[1])], -1)
        scale = earthlight.solid_angle(distance_m) / self.reference_solid_angle * sunlight_factor
        return per_band[:, self.channel_index] * np.asarray(scale)[:, None]


def regime(lux):
    luminance = GREY * np.asarray(lux) / math.pi
    return np.where(luminance >= MESOPIC_CD_M2[1], "photopic",
                    np.where(luminance >= MESOPIC_CD_M2[0], "mesopic", "scotopic"))


def nearest_coast(labels, body, lon_deg, lat_deg, ppd=4):
    """The water cell of a body, beside land, nearest a point (cell centres from 90 N and 0 E)."""
    water = labels == body
    land = labels == 0
    shore = water & (np.roll(land, 1, 0) | np.roll(land, -1, 0) | np.roll(land, 1, 1) | np.roll(land, -1, 1))
    rows, cols = np.nonzero(shore)
    lat, lon = 90 - (rows + 0.5) / ppd, (cols + 0.5) / ppd
    k = int(np.argmax(unit((lon + 180) % 360 - 180, lat) @ unit(lon_deg, lat_deg)))
    return float((lon[k] + 180) % 360 - 180), float(lat[k])


def earth_over_seas(stride=4, step_days=0.25):
    """Share of 2026-2045 the Earth's centre stands above each water cell's horizon, at 1-degree spacing."""
    with np.load(ATLAS, allow_pickle=False) as z:
        labels = z["water_label"][stride // 2::stride, stride // 2::stride]
        lat, lon = z["lat_deg"][stride // 2::stride], z["lon_deg"][stride // 2::stride]
    jd = JD_START + np.arange(0.0, YEARS * 365.25, step_days)
    g = ephemeris.earth_and_sun(jd)
    earth, distance = g["earth"], g["earth_distance_m"][:, None]
    rows, cols = np.nonzero(labels > 0)
    up = unit(lon[cols], lat[rows])
    share, mean, low, high = (np.zeros(len(rows)) for _ in range(4))
    for k in range(0, len(rows), 400):
        # Topocentric elevation of the Earth's centre: sin(el) = (d cos(z) - R) / |d e - R u|.
        cos_z = earth @ up[k:k + 400].T
        el = np.degrees(np.arcsin(np.clip((distance * cos_z - MOON_RADIUS) /
                                          np.sqrt(distance**2 - 2 * MOON_RADIUS * distance * cos_z + MOON_RADIUS**2), -1, 1)))
        share[k:k + 400], mean[k:k + 400] = (el > 0).mean(axis=0), el.mean(axis=0)
        low[k:k + 400], high[k:k + 400] = el.min(axis=0), el.max(axis=0)
    return dict(lat=lat, lon=lon, rows=rows, cols=cols, labels=labels[rows, cols], share=share, mean=mean,
                low=low, high=high, samples=len(jd), step_days=step_days)


def summarise_map(grid):
    bodies = {b["label"]: b for b in json.loads(ATLAS_JSON.read_text())["bodies"]}
    weight = np.cos(np.radians(grid["lat"][grid["rows"]]))
    out = []
    listed = [b for b in set(grid["labels"].tolist()) if b in bodies]      # the atlas names bodies above its size floor
    for body in sorted(listed, key=lambda b: -bodies[b]["area_km2"])[:12]:
        m = grid["labels"] == body
        w = weight[m] / weight[m].sum()
        out.append(dict(body=int(body), names=bodies[body]["water_names"][:3], area_km2=bodies[body]["area_km2"],
                        share_of_area_earth_always_up=float(w @ (grid["share"][m] >= 1)),
                        share_of_area_earth_rises_and_sets=float(w @ ((grid["share"][m] > 0) & (grid["share"][m] < 1))),
                        share_of_area_earth_never_up=float(w @ (grid["share"][m] <= 0)),
                        mean_share_of_time_earth_up=float(w @ grid["share"][m])))
    return out


def build():
    month = json.loads(SHORE_MONTH.read_text())
    tide = month["tide_month"]
    hours = np.array(month["time_hours"], float)
    middle_hour = 1.5 * CYCLE_HOURS
    jd = tide["middle_jd_tt"] + (hours - middle_hour) / 24
    with np.load(ATLAS, allow_pickle=False) as z:
        spa = nearest_coast(z["water_label"], SPA_BODY, *INGENII)
    sites = [(s["name"], s["short"], *s["lon_lat"]) for s in month["shores"]]
    sites += [HEADLAND, ("South Pole–Aitken sea, by Mare Ingenii", "Ingenii coast", *spa)]
    sky = Sky()
    earth_light = Earthlight(sky)
    records = []
    for name, short, lon, lat in sites:
        g = geometry(jd, lon, lat)
        sun = sky.ground(g["sun_elevation"])
        weights = earth_light.weights(g["earth_phase_angle"], g["earth_distance_m"], g["earth_sunlight_factor"])
        earth_above_air = weights @ sky.channel_y
        earth = np.einsum("tc,tc,c->t", weights, sky.ground_channels(g["earth_elevation"]), sky.channel_y)
        total = sun + earth + FLOOR_LUX
        sources = np.array(["sky lit by the Sun", "the Earth", "stars and airglow (placeholder)"])
        dominant = sources[np.argmax(np.stack([sun, earth, np.full_like(sun, FLOOR_LUX)]), axis=0)]
        modes = regime(total)
        night = g["sun_elevation"] < -0.27
        records.append(dict(
            name=name, short=short, lon_lat=[lon, lat],
            hours_by_mode={m: int(np.sum(modes == m)) for m in ("photopic", "mesopic", "scotopic")},
            hours_by_dominant_source={str(s): int(np.sum(dominant == s)) for s in sources},
            night_hours=int(night.sum()),
            darkest_lux=float(total.min()), brightest_night_lux=float(total[night].max()) if night.any() else None,
            earth_elevation_range_deg=[float(g["earth_elevation"].min()), float(g["earth_elevation"].max())],
            hours_earth_up=int(np.sum(g["earth_elevation"] > 0)),
            series=dict(sun_elevation_deg=g["sun_elevation"].round(2).tolist(),
                        sun_azimuth_deg=g["sun_azimuth"].round(1).tolist(),
                        earth_elevation_deg=g["earth_elevation"].round(2).tolist(),
                        earth_azimuth_deg=g["earth_azimuth"].round(1).tolist(),
                        earth_lit_fraction=g["earth_lit_fraction"].round(3).tolist(),
                        earth_above_air_lux=earth_above_air.round(3).tolist(),
                        ground_lux_from_sun=np.round(sun, 5).tolist(),
                        ground_lux_from_earth=np.round(earth, 5).tolist())))
    grid = earth_over_seas()
    MAP.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(MAP, lat_deg=grid["lat"], lon_deg=grid["lon"], rows=grid["rows"].astype(np.int16),
                        cols=grid["cols"].astype(np.int16), labels=grid["labels"].astype(np.int32),
                        share_earth_up=grid["share"].astype(np.float32), mean_elevation_deg=grid["mean"].astype(np.float32),
                        lowest_elevation_deg=grid["low"].astype(np.float32), highest_elevation_deg=grid["high"].astype(np.float32))
    elevations = np.arange(-90, 90.25, 0.5)
    files = [Path(__file__), ROOT / "geography/lunar_ephemeris.py", ROOT / "climate/waves/cycle.py",
             ROOT / "illumination/earthlight/model.py"]
    return dict(
        schema="terluna.research.sea-lighting-calendar/1",
        evidence=("Geometry from a truncated lunar theory checked against JPL Horizons; clear-sky light from the solved "
                  "spherical sky of the design atmosphere (horizontally uniform, no refraction, no clouds); earthlight "
                  "from a model spectrum of the whole Earth (Glenar et al. 2019) scaled to its observed visual phase "
                  "curve (Robinson et al. 2025), carried through the same air by wavelength. Starlight and airglow are "
                  "a placeholder."),
        reading_rule=("time_hours count from the first atmospheric snapshot of the wave runs; the calendar dates are the "
                      "shore month's tide month. Elevations in degrees above the horizon, azimuths clockwise from north. "
                      "ground_lux_* is horizontal illuminance at the surface, clear sky. The mode of vision applies to an "
                      "18% grey surface lit by the total including the 0.001-lux placeholder."),
        producer=dict(domain="research", files={str(p.relative_to(ROOT)): sha256(p) for p in files},
                      constants=constants_used(files)),
        inputs=dict(shore_month_sha256=sha256(SHORE_MONTH), atlas_sha256=sha256(ATLAS), **sky.inputs,
                    **{f"illumination/earthlight/inputs/{f['name']}": f["sha256"]
                       for f in EARTHLIGHT_INPUTS.manifest()["files"]}),
        sky=dict(shielded_sunlight_above_air_lux=sky.above_air, unfiltered_sunlight_above_air_lux=sky.unfiltered_above_air,
                 diffuse_readout_relative_difference_max=sky.readout_check,
                 last_scattering_order_fraction_max=sky.last_order_fraction,
                 ground_lux_by_sun_elevation=dict(elevation_deg=elevations.tolist(),
                                                  lux=np.round(sky.ground(elevations), 6).tolist())),
        month=dict(start_tt=tide["start_tt"], end_tt=tide["end_tt"], middle_jd_tt=tide["middle_jd_tt"]),
        time_hours=hours.round(3).tolist(), sites=records,
        earth_over_seas=dict(path=str(MAP.relative_to(ROOT)), sha256=sha256(MAP), samples=grid["samples"],
                             step_days=grid["step_days"], seas=summarise_map(grid)))


def main():
    record = build()
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(record, separators=(",", ":"), allow_nan=False) + "\n")
    print(RESULT)


if __name__ == "__main__":
    main()
