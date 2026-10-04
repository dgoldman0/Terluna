"""What each rendering frame shows, gathered from the study's products for describing the scenes.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.sea_appearance.scenes

The four rendering coasts at the six moments of the regime results. Each frame has a camera 2 m above the water
(CAMERAS: where, and why there) and a view (VIEWS: the azimuth, clockwise from north, and the pitch of a
65-degree lens on a 16:9 frame). For a 1920 by 1080 frame the runner gathers:
- where the Sun, the Earth and the bright stars stand in it, the disks' colours through the air, the Earth's phase;
- the sky's and the sea's colour and brightness across the frame, read from the regime panoramas
  (research/runs/sea_appearance/regimes.npz), with the glitter of both disks from the reflection model;
- the coast's skyline from LOLA heights at the moment's tide (geography/coast_terrain.py), its distance, the side
  the light falls on and the haze of the air in front of it;
- the sea state from the hour's wave spectrum, and the water's own colour.
Colours are given as sRGB at one exposure per frame (1.25 times the 95th percentile of the sky's luminance in the
frame, highlights rolling off from 80% of white), without adapting to the light of the moment, as the study's
figures show them. The runner writes results/scenes.json; scenes.md describes each scene in words from it. The
renderer in visualization/reference-renderer/seas/ draws the same frames.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from geography import coast_terrain, tides
from illumination.stars import sky as star_sky
from illumination.water_column import model as water_optics
from illumination.water_surface import reflection
from illumination.water_surface.hotfile import read_hotfile
from illumination.water_surface.model import iter_swan, slope_moments, wavenumber
from illumination.water_surface.short_waves import capillary_scales, short_wave_level, short_wave_slopes
from research.studies.sea_appearance import lighting
from shared.constants import AU, EARTH_RADIUS, MOON_RADIUS, MOON_SURFACE_GRAVITY, SUN_RADIUS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE / "results/scenes.json"
REGIMES = HERE / "results/regimes.json"
CALENDAR = HERE / "results/lighting_calendar.json"
SLOPE_SERIES = HERE / "results/sea_slopes.npz"
PANORAMAS = ROOT / "research/runs/sea_appearance/regimes.npz"
SHORE_MONTH = ROOT / "climate/waves/results/shore_month.json"
TIDES = ROOT / "geography/products/tides_28pct_4ppd.npz"
TERRAIN = ROOT / "research/runs/sea_appearance/terrain"
NEARSIDE = ROOT / "research/runs/waves/nearside/dt300"
SMYTHII_HISTORY = ROOT / "research/runs/waves/shore_history/coast_s4_dt900_0_1419"
WATERS = ROOT / "biosphere/living_water/waters.json"
SOLUTION = lighting.SPHERICAL / "moon_1.2atm_standard.npz"
EYE_HEIGHT_M = 2.0
WATER = "productive_coast"
WIDTH, HEIGHT, FOV_DEG = 1920, 1080, 65.0
SHOULDER = 0.8
D65 = (0.3127, 0.3290)
XYZ_TO_SRGB = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])

# Where each camera stands, metres east and north of the study's coast point, and why.
CAMERAS = {
    "W Procellarum": dict(offset_m=(0.0, 0.0), spectrum="nearside",
                          note="at the wave runs' shore point, 1.7 km from the nearest islet"),
    "Smythii headland": dict(offset_m=(-5500.0, 650.0), spectrum="west_face",
                             note="1 km off the headland's west face, beside the shore history's west-face station"),
    "S Nubium": dict(offset_m=(-1500.0, 500.0), spectrum="nearside",
                     note="1.5 km off the west shore of the island on which the wave grid's shore point falls"),
    "Ingenii coast": dict(offset_m=(0.0, 0.0), spectrum="borrowed",
                          note="at the study's coast point, 5.4 km off the coast"),
}
# The view of each frame: azimuth clockwise from north and pitch, both in degrees.
VIEWS = {
    ("W Procellarum", "noon"): (205.0, -3.0), ("W Procellarum", "low Sun"): (262.0, 1.0),
    ("W Procellarum", "after sunset"): (268.0, 3.0), ("W Procellarum", "evening"): (92.0, 3.0),
    ("W Procellarum", "deep evening"): (285.0, 3.0), ("W Procellarum", "darkest"): (94.0, 4.0),
    ("Smythii headland", "noon"): (270.0, -3.0), ("Smythii headland", "low Sun"): (268.0, 1.0),
    ("Smythii headland", "after sunset"): (268.0, 3.0), ("Smythii headland", "evening"): (268.0, 3.0),
    ("Smythii headland", "deep evening"): (266.0, 3.0), ("Smythii headland", "darkest"): (264.5, 2.0),
    ("S Nubium", "noon"): (180.0, -3.0), ("S Nubium", "low Sun"): (270.0, 1.0),
    ("S Nubium", "after sunset"): (266.0, 3.0), ("S Nubium", "evening"): (260.0, 3.0),
    ("S Nubium", "deep evening"): (246.0, 3.0), ("S Nubium", "darkest"): (137.0, 3.0),
    ("Ingenii coast", "noon"): (180.0, -2.0), ("Ingenii coast", "low Sun"): (268.0, 1.0),
    ("Ingenii coast", "after sunset"): (262.0, 3.0), ("Ingenii coast", "evening"): (255.0, 3.0),
    ("Ingenii coast", "deep evening"): (236.0, 3.0), ("Ingenii coast", "darkest"): (180.0, 3.0),
}
TERRAIN_NAMES = {"W Procellarum": "w_procellarum", "Smythii headland": "smythii_headland", "S Nubium": "s_nubium",
                 "Ingenii coast": "ingenii_coast"}
LAND_SOIL = {"W Procellarum": "12001", "Smythii headland": "12001", "S Nubium": "12001", "Ingenii coast": "62231"}
# Names for the brightest stars, by J2000 place in degrees (matched within 0.3 degrees).
NAMED_STARS = {
    "Sirius": (101.287, -16.716), "Canopus": (95.988, -52.696), "Arcturus": (213.915, 19.182),
    "Alpha Centauri": (219.902, -60.834), "Vega": (279.234, 38.784), "Capella": (79.172, 45.998),
    "Rigel": (78.634, -8.202), "Procyon": (114.825, 5.225), "Achernar": (24.429, -57.237),
    "Betelgeuse": (88.793, 7.407), "Hadar": (210.956, -60.373), "Altair": (297.696, 8.868),
    "Acrux": (186.650, -63.099), "Aldebaran": (68.980, 16.509), "Antares": (247.352, -26.432),
    "Spica": (201.298, -11.161), "Pollux": (116.329, 28.026), "Fomalhaut": (344.413, -29.622),
    "Deneb": (310.358, 45.280), "Mimosa": (191.930, -59.689), "Regulus": (152.093, 11.967),
    "Adhara": (104.656, -28.972), "Castor": (113.650, 31.888), "Shaula": (263.402, -37.104),
    "Bellatrix": (81.283, 6.350), "Elnath": (81.573, 28.608), "Miaplacidus": (138.300, -69.717),
    "Alnilam": (84.053, -1.202), "Alnitak": (85.190, -1.943), "Mintaka": (83.002, -0.299),
    "Alnair": (332.058, -46.961), "Alioth": (193.507, 55.960), "Dubhe": (165.932, 61.751),
    "Mirfak": (51.081, 49.861), "Peacock": (306.412, -56.735), "Alkaid": (206.885, 49.313),
    "Polaris": (37.955, 89.264), "Avior": (125.628, -59.509), "Sargas": (264.330, -42.998),
    "Menkalinan": (89.882, 44.948), "Alhena": (99.428, 16.399), "Atria": (252.166, -69.028),
    "Wezen": (107.098, -26.393), "Kaus Australis": (276.043, -34.385), "Alphard": (141.897, -8.659),
    "Hamal": (31.793, 23.462), "Nunki": (283.816, -26.297), "Mirach": (17.433, 35.621),
    "Alpheratz": (2.097, 29.090), "Menkar": (45.570, 4.090), "Diphda": (10.897, -17.987),
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit(elevation_deg, azimuth_deg):
    e, a = np.radians(elevation_deg), np.radians(azimuth_deg)
    return np.stack(np.broadcast_arrays(np.cos(e) * np.sin(a), np.cos(e) * np.cos(a), np.sin(e)), -1)


def camera_site(calendar, coast):
    """East longitude and latitude of a coast's camera."""
    site = next(s for s in calendar["sites"] if s["short"] == coast)
    x, y = CAMERAS[coast]["offset_m"]
    lat, lon = coast_terrain.destination(*site["lon_lat"], np.array(x), np.array(y))
    return float((lon + 180) % 360 - 180), float(lat)


# --------------------------------------------------------------------------------------------- sea spectra
def nearside_spectrum(lon, lat, hour):
    """The nearest restart file's directional spectrum at a nearside shore (as the slope study reads them)."""
    best = None
    for segment in sorted(NEARSIDE.glob("s*_*")):
        end = int(segment.name.split("_")[1])
        if best is None or abs(end - hour) < abs(best[0] - hour):
            best = (end, segment)
    end, segment = best
    if abs(end - hour) > 24:
        raise ValueError("No restart file within a day of the moment")
    spectra = read_hotfile(segment / "end.hot", nodes=[(lon, lat)])
    with np.load(NEARSIDE / "sea.npz", allow_pickle=False) as z:
        node = int(np.argmin(np.hypot(z["lon"] - lon, z["lat"] - lat)))
        depth = float(z["depth"][node])
    return dict(frequency=spectra["frequency"], direction=spectra["direction"], variance=spectra["variance"][0],
                depth=depth, source=f"nearside wave run, restart file at hour {end} ({end - hour:+.0f} h)",
                file=str((segment / "end.hot").relative_to(ROOT)))


def smythii_spectrum(hour, station="west_face"):
    """The Smythii shore history's hourly spectrum at a station."""
    product = json.loads((SMYTHII_HISTORY / "product.json").read_text())
    names = list(product["sites"])
    column = names.index(station)
    bulk = np.loadtxt(SMYTHII_HISTORY / "sites.tbl").reshape(-1, len(names), len(product["columns"]))
    ti, di = product["columns"].index("time_s"), product["columns"].index("depth_m")
    origin = None
    for record in iter_swan(SMYTHII_HISTORY / "sites.spc"):
        if origin is None:
            origin = record["time"]
        h = (record["time"] - origin).total_seconds() / 3600 + product["start_hour"]
        if abs(h - hour) < 0.01:
            depth = float(bulk[np.flatnonzero(np.abs(bulk[:, column, ti] - h * 3600) < 1)[0], column, di])
            return dict(frequency=record["frequency"], direction=record["direction"],
                        variance=record["variance"][column], depth=depth,
                        source=f"Smythii shore history, {station.replace('_', '-')} station, hour {hour:.0f}",
                        file=str((SMYTHII_HISTORY / "sites.spc").relative_to(ROOT)))
    raise ValueError(f"No Smythii spectrum at hour {hour}")


def borrowed_spectrum(month, hours, toward_deg):
    """For a sea without a wave run: the nearside shores' restart spectrum whose slope is the month's median,
    turned to the local wind; deep water."""
    candidates = []
    with np.load(NEARSIDE / "sea.npz", allow_pickle=False) as z:
        lon_n, lat_n, depth_n = z["lon"], z["lat"], z["depth"]
    shores = month["shores"]
    for segment in sorted(NEARSIDE.glob("s*_*")):
        end = int(segment.name.split("_")[1])
        if not hours[0] <= end <= hours[-1]:
            continue
        spectra = read_hotfile(segment / "end.hot", nodes=[tuple(s["lon_lat"]) for s in shores])
        for i, s in enumerate(shores):
            node = int(np.argmin(np.hypot(lon_n - s["lon_lat"][0], lat_n - s["lon_lat"][1])))
            k = wavenumber(spectra["frequency"], float(depth_n[node]), MOON_SURFACE_GRAVITY)
            e = spectra["variance"][i]
            mss = float(np.trapezoid(e.sum(axis=1) * (360 / e.shape[1]) * k ** 2, spectra["frequency"]))
            candidates.append((mss, end, s["short"], spectra["frequency"], spectra["direction"], e,
                               str((segment / "end.hot").relative_to(ROOT))))
    median = float(np.median([c[0] for c in candidates]))
    mss, end, shore, f, direction, e, path = min(candidates, key=lambda c: abs(c[0] - median))
    rad = np.radians(direction)
    weight = e.sum(axis=0)
    mean_dir = math.degrees(math.atan2(np.sum(weight * np.sin(rad)), np.sum(weight * np.cos(rad))))
    turned = (direction + (toward_deg - mean_dir)) % 360
    order = np.argsort(turned)
    return dict(frequency=f, direction=turned[order], variance=e[:, order], depth=1000.0,
                source=(f"borrowed: the {shore} restart spectrum at hour {end}, the month's median resolved slope, "
                        f"turned {toward_deg - mean_dir:+.0f} degrees to the local wind, as deep water"),
                file=path)


def sea_spectrum(coast, hour, calendar, month, hours, toward_deg):
    camera = CAMERAS[coast]
    site = next(s for s in calendar["sites"] if s["short"] == coast)
    if camera["spectrum"] == "nearside":
        return nearside_spectrum(*site["lon_lat"], hour)
    if camera["spectrum"] == "west_face":
        return smythii_spectrum(hour)
    return borrowed_spectrum(month, hours, toward_deg)


def refine(f, e):
    x = np.empty(2 * len(f) - 1); x[::2] = f; x[1::2] = (f[1:] + f[:-1]) / 2
    y = np.empty((len(x), e.shape[1])); y[::2] = e; y[1::2] = (e[1:] + e[:-1]) / 2
    return x, y


def rotated(along, across, toward_deg):
    a = math.radians(toward_deg)
    c, s = math.cos(a), math.sin(a)
    return np.array([[along * c * c + across * s * s, (along - across) * c * s],
                     [(along - across) * c * s, along * s * s + across * c * c]])


def sea_state(spec, u_star, toward_deg, u10):
    """Bulk description of the hour's sea, and the slope covariance the glitter needs."""
    g = MOON_SURFACE_GRAVITY
    f, theta, e, depth = spec["frequency"], spec["direction"], spec["variance"], spec["depth"]
    m0 = float(np.trapezoid(e.sum(axis=1) * (360 / e.shape[1]), f))
    fp = float(f[np.argmax(e.sum(axis=1))])
    kp = float(wavenumber(np.array([fp]), depth, g)[0])
    weight = np.trapezoid(e, f, axis=0)
    rad = np.radians(theta)
    c, s = np.sum(weight * np.cos(rad)), np.sum(weight * np.sin(rad))
    mean_toward = math.degrees(math.atan2(s, c)) % 360
    spread = math.degrees(math.sqrt(2 * max(0.0, 1 - math.hypot(c, s) / weight.sum())))
    ff, ee = f, e
    for _ in range(2):
        ff, ee = refine(ff, ee)
    resolved = slope_moments(ff, theta, ee, depth, g)["covariance"]
    k_from = float(wavenumber(np.array([f[-1]]), depth, g)[0])
    short = short_wave_slopes(u_star, g, k_from, peak_speed=2 * math.pi * fp / kp)
    short_cov = rotated(short["along"], short["across"], toward_deg)
    k_m, _ = capillary_scales(g)
    compass = lambda math_deg: (90 - math_deg) % 360          # Cartesian 'toward' to compass bearing
    return dict(
        source=spec["source"], file=spec["file"], depth_m=round(spec["depth"], 1),
        significant_height_m=round(4 * math.sqrt(m0), 3), peak_period_s=round(1 / fp, 2),
        peak_wavelength_m=round(2 * math.pi / kp, 2), waves_travel_toward_compass_deg=round(compass(mean_toward), 1),
        directional_spread_deg=round(spread, 1), steepness_hs_over_wavelength=round(4 * math.sqrt(m0) * kp / (2 * math.pi), 4),
        wind_10m_m_s=round(u10, 2), wind_toward_compass_deg=round(compass(toward_deg), 1),
        friction_velocity_m_s=round(u_star, 4), short_waves_present=bool(short_wave_level(u_star, g) > 0),
        capillary_wavelength_cm=round(200 * math.pi / k_m, 2),
        mean_square_slope=dict(longer_than_1m=round(float(np.trace(resolved)), 5),
                               short=round(float(np.trace(short_cov)), 5),
                               total=round(float(np.trace(resolved + short_cov)), 5)),
        rms_tilt_deg=round(math.degrees(math.atan(math.sqrt(float(np.trace(resolved + short_cov))))), 2)), \
        resolved + short_cov


# ------------------------------------------------------------------------------------------------- the frame
class Frame:
    def __init__(self, azimuth, pitch):
        self.f = unit(pitch, azimuth)
        self.r = np.array([math.cos(math.radians(azimuth)), -math.sin(math.radians(azimuth)), 0.0])
        self.u = np.cross(self.r, self.f)
        self.th = math.tan(math.radians(FOV_DEG) / 2)
        self.aspect = WIDTH / HEIGHT

    def direction(self, x, y):
        d = (self.f[None] + ((2 * x / WIDTH - 1) * self.th)[:, None] * self.r[None]
             + ((1 - 2 * y / HEIGHT) * self.th / self.aspect)[:, None] * self.u[None])
        return d / np.linalg.norm(d, axis=1, keepdims=True)

    def project(self, d):
        d = np.atleast_2d(d)
        z = d @ self.f
        x = (d @ self.r / z / self.th + 1) / 2 * WIDTH
        y = (1 - d @ self.u / z / self.th * self.aspect) / 2 * HEIGHT
        return np.where(z > 0, x, np.nan), np.where(z > 0, y, np.nan)


def panorama_sampler(panorama, elevations, azimuths):
    el_asc = elevations[::-1]
    pan = panorama[::-1]

    def sample(el, az):
        el = np.clip(el, el_asc[0], el_asc[-1])
        i = np.clip(np.searchsorted(el_asc, el) - 1, 0, len(el_asc) - 2)
        u = ((el - el_asc[i]) / (el_asc[i + 1] - el_asc[i]))[..., None]
        q = (az - azimuths[0]) % 360
        j = np.floor(q).astype(int) % 360
        v = (q - np.floor(q))[..., None]
        j1 = (j + 1) % 360
        return ((1 - u) * ((1 - v) * pan[i, j] + v * pan[i, j1]) + u * ((1 - v) * pan[i + 1, j] + v * pan[i + 1, j1]))
    return sample


def to_hex(xyz, white):
    """sRGB for XYZ at an exposure: luminance above 80% of white rolls off toward white, keeping chromaticity."""
    xyz = np.asarray(xyz, float)
    y = xyz[1] / white
    rolled = y if y <= SHOULDER else SHOULDER + (1 - SHOULDER) * (1 - math.exp(-(y - SHOULDER) / (1 - SHOULDER)))
    lin = np.clip((xyz / white * (rolled / y if y > 0 else 0)) @ XYZ_TO_SRGB.T, 0, 1)
    enc = np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * lin ** (1 / 2.4) - 0.055)
    return "#" + "".join(f"{int(round(255 * c)):02x}" for c in enc)


def chromaticity(xyz):
    t = float(np.sum(xyz))
    return [round(float(xyz[0]) / t, 4), round(float(xyz[1]) / t, 4)] if t > 0 else None


def colour_name(srgb):
    """Plain colour words for a display colour (#rrggbb): hue, saturation and lightness as seen on screen."""
    r, g, b = (int(srgb[i:i + 2], 16) / 255 for i in (1, 3, 5))
    hi, lo = max(r, g, b), min(r, g, b)
    sat = 0.0 if hi == 0 else (hi - lo) / hi
    if hi < 0.06:
        return "black"
    if sat < 0.07:
        tone = "white" if hi > 0.92 else "light grey" if hi > 0.7 else "grey" if hi > 0.4 else "dark grey"
        return tone
    if hi == r:
        hue = (60 * (g - b) / (hi - lo)) % 360
    elif hi == g:
        hue = 60 * (b - r) / (hi - lo) + 120
    else:
        hue = 60 * (r - g) / (hi - lo) + 240
    names = [(12, "red"), (28, "red-orange"), (42, "orange"), (52, "amber"), (64, "golden yellow"), (75, "yellow"),
             (100, "yellow-green"), (150, "green"), (172, "aqua"), (195, "cyan"), (218, "sky blue"), (245, "blue"),
             (275, "violet"), (320, "magenta-pink"), (345, "rose"), (361, "red")]
    name = next(n for limit, n in names if hue < limit)
    if sat < 0.2:
        name = ("greyish " if hi < 0.75 else "very pale ") + name
    elif sat < 0.4:
        name = ("muted " if hi < 0.6 else "pale ") + name
    elif sat > 0.8 and hi > 0.6:
        name = "vivid " + name
    if hi < 0.3:
        name = "very dark " + name
    elif hi < 0.55:
        name = "dark " + name
    return name


# ------------------------------------------------------------------------------------------------- the coast
def skyline(height, axis, camera_xy, offset, azimuths):
    """Per azimuth: the highest elevation angle of land above the water (deg), its distance (m) and height (m),
    and the nearest land's distance along the line."""
    dist = np.arange(50.0, 1.05e5, 50.0)
    out = np.full((len(azimuths), 4), np.nan)
    dx = axis[1] - axis[0]
    for i, az in enumerate(np.radians(azimuths)):
        x = camera_xy[0] + dist * math.sin(az)
        y = camera_xy[1] + dist * math.cos(az)
        u, v = (x - axis[0]) / dx, (y - axis[0]) / dx
        inside = (u >= 0) & (v >= 0) & (u < len(axis) - 1) & (v < len(axis) - 1)
        if not inside.any():
            continue
        ui, vi = u[inside], v[inside]
        i0, j0 = ui.astype(int), vi.astype(int)
        fu, fv = ui - i0, vi - j0
        z = ((1 - fu) * (1 - fv) * height[j0, i0] + fu * (1 - fv) * height[j0, i0 + 1]
             + (1 - fu) * fv * height[j0 + 1, i0] + fu * fv * height[j0 + 1, i0 + 1]) + offset
        d = dist[inside]
        land = z > 0
        if not land.any():
            continue
        angle = np.degrees(np.arctan2(z - EYE_HEIGHT_M - d ** 2 / (2 * MOON_RADIUS), d))
        angle = np.where(land, angle, -90.0)
        k = int(np.argmax(angle))
        out[i] = [angle[k], d[k], z[k], d[land][0]]
    return out


# --------------------------------------------------------------------------------------------------- builder
class Inputs:
    def __init__(self):
        self.regimes = json.loads(REGIMES.read_text())
        self.calendar = json.loads(CALENDAR.read_text())
        self.month = json.loads(SHORE_MONTH.read_text())
        self.hours = np.array(self.calendar["time_hours"], float)
        self.jd0 = self.month["tide_month"]["middle_jd_tt"]
        with np.load(SLOPE_SERIES) as z:
            self.series = {k: z[k] for k in z.files}
        with np.load(PANORAMAS, allow_pickle=False) as z:
            self.panoramas = {k: z[k] for k in z.files}
        self.sky = lighting.Sky()
        with np.load(SOLUTION, allow_pickle=False) as z:
            self.wavelength, self.xyz = z["wavelength"], z["xyz"].astype(float)
            self.extinction = z["scattering"][0] + z["absorption"][0]          # per metre at the ground
        self.sky.channel_wavelength = self.wavelength
        self.earth_light = lighting.Earthlight(self.sky)
        self.waters = json.loads(WATERS.read_text())
        self.n = float(reflection.refractive_index(550.0))
        with np.load(TIDES, allow_pickle=False) as z:
            self.tides = {k: z[k] for k in ("lat_deg", "lon_deg", "rows", "cols", "basis")}
        days = np.arange(0.0, tides.YEARS * 365.25, tides.STEP_DAYS)
        self.tide_mean = tides.forcing_coefficients(tides.JD_START + days).mean(axis=0)
        self.stars = star_sky.catalogue()

    def jd(self, hour):
        return self.jd0 + (hour - 1.5 * lighting.CYCLE_HOURS) / 24

    def tide(self, lon, lat, jd):
        t = self.tides
        lat_c, lon_c = t["lat_deg"][t["rows"]], t["lon_deg"][t["cols"]]
        pixel = int(np.argmax(lighting.unit(lon_c, lat_c) @ lighting.unit(lon, lat)))
        return float(((tides.forcing_coefficients(np.array([jd])) - self.tide_mean) @ t["basis"][:, pixel])[0])

    def beam(self, elevation):
        """Direct normal irradiance per channel and unit source at a source elevation (zero below the horizon)."""
        if elevation <= -0.3:
            return np.zeros(self.sky.direct_channels.shape[1])
        horizontal = np.array([np.interp(elevation, self.sky.beam_angles, self.sky.direct_channels[:, c])
                               for c in range(self.sky.direct_channels.shape[1])])
        return horizontal / math.sin(math.radians(max(elevation, 0.5)))


def hours_from_sunset(series, index):
    sun = np.asarray(series["sun_elevation_deg"])
    n = len(sun)
    if sun[index] > 0:
        for step in range(1, n):                              # until the next sunset
            if sun[(index + step - 1) % n] > 0 >= sun[(index + step) % n]:
                return -float(step)
    for step in range(1, n):                                  # since the last sunset
        if sun[(index - step) % n] > 0 >= sun[(index - step + 1) % n]:
            return float(step - 1)
    return None


def scene(inputs, index, record):
    coast, moment, hour = record["coast"], record["moment"], record["hour"]
    lon, lat = camera_site(inputs.calendar, coast)
    jd = inputs.jd(hour)
    g = lighting.geometry(np.array([jd]), lon, lat)
    sun_el, sun_az = float(g["sun_elevation"][0]), float(g["sun_azimuth"][0])
    earth_el, earth_az = float(g["earth_elevation"][0]), float(g["earth_azimuth"][0])
    weights = inputs.earth_light.weights(g["earth_phase_angle"], g["earth_distance_m"],
                                         g["earth_sunlight_factor"])[0]
    azimuth, pitch = VIEWS[(coast, moment)]
    frame = Frame(azimuth, pitch)
    # The sea state and its slopes at that hour.
    key = coast.replace(" ", "_").lower()
    i = int(np.argmin(np.abs(inputs.hours - hour)))
    u_star, toward = float(inputs.series[f"{key}/u_star"][i]), float(inputs.series[f"{key}/toward_deg"][i])
    u10 = float(inputs.series[f"{key}/u10"][i])
    spec = sea_spectrum(coast, hour, inputs.calendar, inputs.month, inputs.hours, toward)
    waves, covariance = sea_state(spec, u_star, toward, u10)
    view_compass = azimuth
    relative = (waves["waves_travel_toward_compass_deg"] - view_compass + 540) % 360 - 180
    waves["travel_relative_to_view_deg"] = round(relative, 1)
    # The coast at the tide's level.
    tide = inputs.tide(lon, lat, jd)
    with np.load(TERRAIN / f"{TERRAIN_NAMES[coast]}.npz", allow_pickle=False) as z:
        height, axis = z["height_m"].astype(float), z["x_m"].astype(float)
    camera_xy = CAMERAS[coast]["offset_m"]
    half_h = math.degrees(math.atan(frame.th / frame.aspect))
    azimuths = (azimuth + np.arange(-40.0, 40.01, 0.25)) % 360
    line = skyline(height, axis, camera_xy, -tide, azimuths)
    dip = math.degrees(math.sqrt(2 * EYE_HEIGHT_M / MOON_RADIUS))
    # A coarse sampling of the frame: sky, land and sea, and the exposure.
    xs, ys = np.meshgrid(np.linspace(5, WIDTH - 5, 192), np.linspace(5, HEIGHT - 5, 108))
    d = frame.direction(xs.ravel(), ys.ravel())
    el = np.degrees(np.arcsin(np.clip(d[:, 2], -1, 1)))
    az = np.degrees(np.arctan2(d[:, 0], d[:, 1])) % 360
    sample = panorama_sampler(inputs.panoramas["xyz"][index].astype(float), inputs.panoramas["elevation_deg"],
                              inputs.panoramas["azimuth_deg"])
    xyz = sample(el, az)
    rel = (az - azimuth + 540) % 360 - 180
    top = np.interp(rel, np.arange(-40.0, 40.01, 0.25), np.nan_to_num(line[:, 0], nan=-90.0))
    is_land = (el <= top) & (top > -dip) & (el > -dip - 0.05)      # between the horizon and the skyline
    is_sea = (el <= -dip) & ~is_land
    is_sky = ~is_land & ~is_sea
    white = 1.25 * float(np.percentile(xyz[is_sky, 1], 95)) if is_sky.any() else 1.25 * float(np.percentile(xyz[:, 1], 95))

    def point(label, x, y):
        dd = frame.direction(np.array([float(x)]), np.array([float(y)]))[0]
        e = math.degrees(math.asin(dd[2]))
        a = math.degrees(math.atan2(dd[0], dd[1])) % 360
        value = sample(np.array([e]), np.array([a]))[0]
        return dict(label=label, x=int(x), y=int(y), elevation_deg=round(e, 2), azimuth_deg=round(a, 1),
                    luminance_cd_m2=round(float(value[1]), 5), chromaticity_xy=chromaticity(value),
                    colour=colour_name(to_hex(value, white)), srgb=to_hex(value, white))

    def row_of(elevation, x=WIDTH / 2):
        """Pixel row of an elevation at a column (by bisection on the frame's rays)."""
        lo, hi = -2000.0, 3000.0
        for _ in range(50):
            mid = (lo + hi) / 2
            e = math.degrees(math.asin(frame.direction(np.array([float(x)]), np.array([mid]))[0][2]))
            lo, hi = (mid, hi) if e > elevation else (lo, mid)
        return (lo + hi) / 2

    horizon_y = row_of(-dip)
    sky_points, sea_points = [], []
    for label, e in (("top of frame", None), ("15 degrees up", 15.0), ("8 degrees up", 8.0), ("3 degrees up", 3.0),
                     ("1 degree up", 1.0), ("just above the horizon", 0.25)):
        for side, x in (("left edge", 40), ("centre", WIDTH / 2), ("right edge", WIDTH - 40)):
            y = 30 if e is None else row_of(e, x)
            if 0 <= y < HEIGHT:
                sky_points.append(point(f"{label}, {side}", x, y))
    for label, e in (("just below the horizon", -0.4), ("1 degree down", -1.0), ("3 degrees down", -3.0),
                     ("8 degrees down", -8.0), ("bottom of frame", None)):
        for side, x in (("left edge", 40), ("centre", WIDTH / 2), ("right edge", WIDTH - 40)):
            y = HEIGHT - 30 if e is None else row_of(e, x)
            if 0 <= y < HEIGHT:
                sea_points.append(point(f"{label}, {side}", x, y))
    # The disks.
    out_disks = {}
    sun_dir, earth_dir = unit(sun_el, sun_az), unit(earth_el, earth_az)
    sun_radius = math.asin(SUN_RADIUS / (AU * float(g["earth_sunlight_factor"][0]) ** -0.5))
    earth_radius = math.asin(EARTH_RADIUS / float(g["earth_distance_m"][0]))
    pixel_deg = FOV_DEG / WIDTH
    for name, direction, radius, channels, lit in (
            ("sun", sun_dir, sun_radius, inputs.beam(sun_el), 1.0),
            ("earth", earth_dir, earth_radius, inputs.beam(earth_el) * weights, float(g["earth_lit_fraction"][0]))):
        x, y = frame.project(direction)
        normal = channels @ inputs.xyz
        omega = 2 * math.pi * (1 - math.cos(radius)) * max(lit, 1e-4)
        info = dict(elevation_deg=round(float(math.degrees(math.asin(direction[2]))), 2),
                    azimuth_deg=round(float(math.degrees(math.atan2(direction[0], direction[1])) % 360), 1),
                    diameter_deg=round(math.degrees(2 * radius), 3),
                    diameter_px=round(math.degrees(2 * radius) / pixel_deg, 1),
                    above_horizon=bool(direction[2] > -math.sin(radius)),
                    in_frame=bool(np.isfinite(x[0]) and 0 <= x[0] < WIDTH and 0 <= y[0] < HEIGHT),
                    x=None if not np.isfinite(x[0]) else round(float(x[0]), 1),
                    y=None if not np.isfinite(y[0]) else round(float(y[0]), 1))
        if normal[1] > 0:
            info.update(normal_illuminance_lux=round(float(normal[1]), 6),
                        disk_luminance_cd_m2=round(float(normal[1] / omega), 3),
                        times_white=round(float(normal[1] / omega / white), 1),
                        colour_chromaticity=chromaticity(normal),
                        colour=colour_name(to_hex(normal / normal[1] * 0.75 * white, white)),
                        disk_srgb=to_hex(normal / normal[1] * 0.75 * white, white))
            if name == "earth":
                above = weights @ inputs.xyz
                info["above_the_air_colour"] = colour_name(to_hex(above / above[1] * 0.75 * white, white))
        if name == "earth":
            info.update(lit_fraction=round(float(g["earth_lit_fraction"][0]), 3),
                        phase_angle_deg=round(float(g["earth_phase_angle"][0]), 1))
            # Which way the lit side faces in the frame: the Sun's direction across the Earth's disk.
            t = sun_dir - (sun_dir @ earth_dir) * earth_dir
            if np.linalg.norm(t) > 1e-9 and info["x"] is not None:
                ahead = earth_dir + 0.01 * t / np.linalg.norm(t)
                x2, y2 = frame.project(ahead)
                angle = math.degrees(math.atan2(-(y2[0] - y[0]), x2[0] - x[0]))
                info["lit_side_toward_deg"] = round(angle % 360, 0)
                info["lit_side_toward"] = ["right", "upper right", "top", "upper left", "left", "lower left",
                                           "bottom", "lower right"][int(((angle % 360) + 22.5) // 45) % 8]
        out_disks[name] = info
    # Glitter of both disks on the frame's sea, from the reflection model with this hour's slopes.
    sea_d = d[is_sea]
    views = -sea_d
    glitter = {}
    reflected_sea = xyz[is_sea]
    for name, direction, radius, channels in (("sun", sun_dir, sun_radius, inputs.beam(sun_el)),
                                              ("earth", earth_dir, earth_radius, inputs.beam(earth_el) * weights)):
        normal = channels @ inputs.xyz
        if normal[1] <= 0 or not len(views) or direction[2] < -0.02:
            continue
        k = reflection.disk_glint(direction, radius, views, covariance, inputs.n)
        lum = k * normal[1]
        strong = lum > 0.5 * reflected_sea[:, 1]
        if not strong.any():
            continue
        px, py = xs.ravel()[is_sea][strong], ys.ravel()[is_sea][strong]
        peak = int(np.argmax(lum))
        glitter[name] = dict(peak_luminance_cd_m2=round(float(lum[peak]), 3),
                             peak_x=int(xs.ravel()[is_sea][peak]), peak_y=int(ys.ravel()[is_sea][peak]),
                             x_range=[int(px.min()), int(px.max())], y_range=[int(py.min()), int(py.max())],
                             share_of_frame_sea=round(float(strong.mean()), 4),
                             colour=colour_name(to_hex(normal / normal[1] * 0.75 * white, white)),
                             times_white=round(float(lum[peak] / white), 1))
    # The water's own light seen straight down.
    downwelling = inputs.sky.ground_channels(sun_el)[0] + weights * inputs.sky.ground_channels(earth_el)[0]
    w = inputs.waters["waters"][WATER]
    iop = water_optics.water(inputs.wavelength, w["chlorophyll_mg_m3"]["best"], w["dissolved_organic_440_per_m"]["best"],
                             w["fines_g_m3"]["best"], w["soil"])
    leaving = (water_optics.remote_sensing_reflectance(iop["a"], iop["bb"]) * downwelling) @ inputs.xyz
    # The coast in the frame.
    lam, refl = water_optics.soil_reflectance(LAND_SOIL[coast])
    albedo = np.interp(inputs.wavelength, lam, refl)
    diffuse = np.array([np.interp(sun_el, inputs.sky.angles, inputs.sky.diffuse_channels[:, c])
                        for c in range(inputs.sky.diffuse_channels.shape[1])])
    diffuse += weights * np.array([np.interp(earth_el, inputs.sky.angles, inputs.sky.diffuse_channels[:, c])
                                   for c in range(inputs.sky.diffuse_channels.shape[1])])
    land = []
    in_frame = np.abs(np.arange(-40.0, 40.01, 0.25)) <= FOV_DEG / 2 + 2
    seen = in_frame & np.isfinite(line[:, 0]) & (line[:, 0] > -dip + 0.002)
    segments, start = [], None
    for j, flag in enumerate(seen):
        if flag and start is None:
            start = j
        if (not flag or j == len(seen) - 1) and start is not None:
            segments.append((start, j if flag else j - 1))
            start = None
    for a, b in segments:
        part = line[a:b + 1]
        k = int(np.nanargmax(part[:, 0]))
        az_from, az_to = azimuths[a], azimuths[b]
        x_from = float(frame.project(unit(0.0, az_from))[0][0])
        x_to = float(frame.project(unit(0.0, az_to))[0][0])
        peak_dir = unit(float(part[k, 0]), azimuths[a + k])
        px, py = frame.project(peak_dir)
        distance = float(part[k, 1])
        # Which face of the land the camera sees, and how lit it is.
        to_camera = (azimuths[a + k] + 180) % 360

        def facing(el_, az_):
            # Cosine of incidence on a slope tilted 20 degrees toward the camera.
            e_, d_ = math.radians(el_), math.radians(az_ - to_camera)
            return math.cos(math.radians(20)) * math.sin(e_) + math.sin(math.radians(20)) * math.cos(e_) * math.cos(d_)

        sun_cos = facing(sun_el, sun_az) if sun_el > 0 else 0.0
        earth_cos = facing(earth_el, earth_az) if earth_el > 0 else 0.0
        sun_side, earth_side = sun_cos > 0.1, earth_cos > 0.1
        direct = inputs.beam(sun_el) * max(sun_cos, 0.0) + inputs.beam(earth_el) * weights * max(earth_cos, 0.0)
        radiance = albedo * (diffuse * 0.7 + direct) / math.pi
        transmit = np.exp(-inputs.extinction * distance)
        horizon_sky = sample(np.array([0.5]), np.array([azimuths[a + k]]))[0]
        t_xyz = (transmit[:, None] * inputs.xyz).sum(0) / inputs.xyz.sum(0)
        apparent = t_xyz * (radiance @ inputs.xyz) + (1 - t_xyz) * horizon_sky
        nearest = int(np.nanargmin(part[:, 3]))
        mean_level = skyline(height, axis, camera_xy, 0.0, np.array([azimuths[a + nearest]]))[0]
        land.append(dict(
            shore_now_km=round(float(part[nearest, 3]) / 1e3, 3),
            shore_at_mean_level_km=None if not np.isfinite(mean_level[3]) else round(float(mean_level[3]) / 1e3, 3),
            shore_azimuth_deg=round(float(azimuths[a + nearest]), 1),
            x_from=round(x_from), x_to=round(x_to), azimuth_from_deg=round(float(az_from), 1),
            azimuth_to_deg=round(float(az_to), 1), highest_elevation_deg=round(float(part[k, 0]), 3),
            highest_y=round(float(py[0]), 1), highest_x=round(float(px[0]), 1),
            height_above_horizon_px=round(float(horizon_y - py[0]), 1),
            highest_point_distance_km=round(distance / 1e3, 2), highest_point_height_m=round(float(part[k, 2])),
            nearest_shore_km=round(float(np.nanmin(part[:, 3])) / 1e3, 2),
            farthest_skyline_km=round(float(np.nanmax(part[:, 1])) / 1e3, 2),
            lit_face_toward_camera=bool(sun_side), earthlit_face_toward_camera=bool(earth_side),
            haze_transmittance_550nm=round(float(np.interp(550, inputs.wavelength[np.argsort(inputs.wavelength)],
                                                           transmit[np.argsort(inputs.wavelength)])), 3),
            haze_transmittance_450nm=round(float(np.interp(450, inputs.wavelength[np.argsort(inputs.wavelength)],
                                                           transmit[np.argsort(inputs.wavelength)])), 3),
            apparent_luminance_cd_m2=round(float(apparent[1]), 5), apparent_colour=colour_name(to_hex(apparent, white)),
            apparent_srgb=to_hex(apparent, white)))
    # The stars: in frame, above the terrain and the sea horizon, by brightness, with the eye's limit there.
    stars = inputs.stars
    body = star_sky.body_directions(stars[:, 0], stars[:, 1], jd)
    east, north, up = lighting.local_frame(lon, lat)
    local = np.stack([body @ east, body @ north, body @ up], -1)
    s_el = np.degrees(np.arcsin(np.clip(local[:, 2], -1, 1)))
    s_az = np.degrees(np.arctan2(local[:, 0], local[:, 1])) % 360
    sx, sy = frame.project(local)
    inside = np.isfinite(sx) & (sx >= 0) & (sx < WIDTH) & (sy >= 0) & (sy < HEIGHT) & (s_el > 0.3)
    s_rel = (s_az - azimuth + 540) % 360 - 180
    terrain_top = np.interp(s_rel, np.arange(-40.0, 40.01, 0.25), np.nan_to_num(line[:, 0], nan=-90.0))
    inside &= s_el > terrain_top
    listed = []
    for j in np.flatnonzero(inside)[np.argsort(stars[inside, 2])]:
        background = float(sample(np.array([s_el[j]]), np.array([s_az[j]]))[0][1])
        limit = float(star_sky.naked_eye_limit(background))
        transmit = inputs.beam(s_el[j]) @ inputs.xyz[:, 1] / inputs.xyz[:, 1].sum()
        dimmed = stars[j, 2] - 2.5 * math.log10(max(transmit, 1e-12))
        name = next((n for n, (ra, dec) in NAMED_STARS.items()
                     if math.degrees(math.acos(min(1.0, math.cos(math.radians(dec)) * math.cos(stars[j, 1])
                                                    * math.cos(stars[j, 0] - math.radians(ra))
                                                    + math.sin(math.radians(dec)) * math.sin(stars[j, 1])))) < 0.3),
                    None)
        pixel_sr = math.radians(FOV_DEG / WIDTH) ** 2
        point = star_sky.V0_LUX * 10 ** (-0.4 * dimmed) / (4 * pixel_sr)
        listed.append(dict(name=name, v_above_air=round(float(stars[j, 2]), 2), v_seen=round(dimmed, 2),
                           b_minus_v=round(float(stars[j, 3]), 2), x=round(float(sx[j])), y=round(float(sy[j])),
                           elevation_deg=round(float(s_el[j]), 2), naked_eye_limit=round(limit, 2),
                           eye_sees_it=bool(dimmed < limit),
                           long_exposure_shows_it=bool(point > 0.2 * background and point > 0.05 * white)))
    eye = [s for s in listed if s["eye_sees_it"]]
    photo = [s for s in listed if s["long_exposure_shows_it"]]
    total = record["light_on_ground_lux"]
    series = next(s for s in inputs.calendar["sites"] if s["short"] == coast)["series"]
    return dict(
        coast=coast, moment=moment, hour=hour, day_of_month=round((hour - inputs.hours[0]) / 24, 2),
        hours_from_sunset=hours_from_sunset(series, int(np.argmin(np.abs(inputs.hours - hour)))),
        camera=dict(lon_lat_deg=[round(lon, 5), round(lat, 5)], note=CAMERAS[coast]["note"],
                    eye_height_m=EYE_HEIGHT_M, view_azimuth_deg=azimuth, pitch_deg=pitch,
                    horizontal_fov_deg=FOV_DEG, vertical_fov_deg=round(2 * half_h, 2), frame=[WIDTH, HEIGHT],
                    horizon_row_at_centre=round(horizon_y, 1), horizon_dip_deg=round(dip, 4)),
        light=dict(ground_lux=total, vision=str(lighting.regime(total)), white_cd_m2=round(white, 5),
                   share_of_frame=dict(sky=round(float(is_sky.mean()), 3), land=round(float(is_land.mean()), 3),
                                       sea=round(float(is_sea.mean()), 3))),
        sun=out_disks["sun"], earth=out_disks["earth"], sky=sky_points, sea=sea_points, glitter=glitter,
        water_leaving=dict(luminance_cd_m2=round(float(leaving[1]), 4), chromaticity_xy=chromaticity(leaving),
                           colour=colour_name(to_hex(leaving / max(leaving[1], 1e-30) * 0.5 * white, white)),
                           srgb=to_hex(leaving, white)),
        waves=waves, tide_m=round(tide, 3), land=land,
        stars=dict(in_frame=len(listed), eye_sees=len(eye), long_exposure_shows=len(photo),
                   brightest=listed[:14]))


def build():
    inputs = Inputs()
    scenes = [scene(inputs, i, r) for i, r in enumerate(inputs.regimes["scenes"])]
    files = [Path(__file__), ROOT / "illumination/stars/sky.py", ROOT / "geography/coast_terrain.py"]
    return dict(
        schema="terluna.research.sea-scenes/1",
        evidence=("What each of the 24 rendering frames holds, gathered from the study's products: the regime "
                  "panoramas (statistical sea, no coastline) for the colours of sky and sea, LOLA terrain at the "
                  "tide's level for the coast, the hour's wave spectrum, the reflection model for the glitter and "
                  "the bright-star catalogue. The coast's colour is a placeholder: bare Apollo soil, lit by the "
                  "clear sky and hazed by Beer-Lambert extinction toward the horizon sky."),
        reading_rule=("Pixel coordinates are for a 1920 by 1080 frame from the top left; azimuths clockwise from "
                      "north; elevations from the eye's horizontal. sRGB colours are at the frame's exposure "
                      "(white_cd_m2 shown white, highlights rolling off from 80%), on a display with D65 white, "
                      "unadapted. Wave directions are compass bearings the waves travel toward."),
        producer=dict(lane="research", runner=str(Path(__file__).relative_to(ROOT)),
                      files={str(p.relative_to(ROOT)): sha256(p) for p in files}, constants=constants_used(files[:1])),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (REGIMES, CALENDAR, SLOPE_SERIES, PANORAMAS, TIDES,
                                                             WATERS)},
        scenes=scenes)


def main():
    product = build()
    RESULT.write_text(json.dumps(product, indent=1, allow_nan=False) + "\n")
    for s in product["scenes"]:
        print(s["coast"], s["moment"], "white", s["light"]["white_cd_m2"], "land", len(s["land"]),
              "stars eye", s["stars"]["eye_sees"], flush=True)


if __name__ == "__main__":
    main()
