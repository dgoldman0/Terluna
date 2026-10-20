"""Living light on the Open Moon: what bioluminescence can give, what people need outdoors at night, and the skyglow.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.living_light.run    # ignored research/runs/living_light/

One screen in three parts on committed products and cited values:
- emitters: the light and the cost of luminous organisms and engineered plants, from measured quantum yields and
  photon outputs, turned into luminance, path illuminance, ATP and carbon per photon and per square metre through a
  354-hour night, set against the plant model's night stores (biosphere/canopy/results/plant.json), and taken to the
  physical limits;
- needs: orientation, walking, faces and reading against each setting's natural night light, from the calendar's
  evaluator (illumination/calendar) on its transfer product, checked against the provisioning branch's daily-life
  calendars and main's per-place Earthlight, and over the ways-of-living homes;
- skyglow: the land under a light-polluted sky a person, from ecology's people screen, carried to the Moon's air by a
  single-scattering skyglow model and to each setting's natural sky, with twilight-aware lamps and living light.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np

from illumination.calendar.geometry import geometry_from_vectors
from illumination.calendar.validate import Response
from shared.constants import (AU, AVOGADRO, EARTH_MOON_DISTANCE, EARTH_RADIUS, JULIAN_DAY, MOON_RADIUS, PLANCK,
                              SPEED_OF_LIGHT, STANDARD_GRAVITY, SUN_RADIUS, SYNODIC_MONTH_DAYS)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.living-light/1'
OUT = HERE / 'results' / 'living_light.json'
SOURCES = 'research/studies/living_light/sources.json'

INPUTS = {
    'illumination/calendar/results/transfer.json': 'terluna.illumination.calendar-transfer/1',
    'illumination/sky/results/solved_sky.json': 'terluna.illumination.solved-spherical-sky/1',
    'illumination/aerosol/results/haze.json': 'terluna.illumination.haze/1',
    'climate/results/gcm/global_winds_A28_dim5_moon.json': 'terluna.climate.gcm-global-winds/1',
    'biosphere/canopy/results/plant.json': 'terluna.biosphere.plant-carbon-cycle/1',
    'biosphere/ecology/results/people.json': 'terluna.biosphere.ecology-people/1',
    'biosphere/ecology/results/taxa.json': 'terluna.biosphere.ecology-taxa/1',
}
# The calendar's evaluator, called live as the provisioning branch's daily-life study calls it.
IMPORTED_CODE = ('illumination/calendar/geometry.py', 'illumination/calendar/validate.py', 'illumination/calendar/disk.py')

# Files from the former domain branches, bound by blob and SHA-256 and read in-tree since the joint integration.
MAIN_COMMIT = None
PROVISIONING_COMMIT = None
CROSS_BRANCH = {
    'research/studies/joint_synthesis/results/places.json': dict(
        branch='main', commit=MAIN_COMMIT, blob='3aa6d0939faf502433c03105b831a9994ab3315d',
        sha256='9707e4d1c7332ba352b90dbf77c66fde080ec5e0fe24e57b54067c61ecdcfc8d',
        schema='terluna.research.joint-synthesis-places/1'),
    'research/studies/daily_life/results/daily_life.json': dict(
        branch='main', commit=PROVISIONING_COMMIT, blob='bb3b666527ee6fc3cac21fe50ac97e74e0c1f28c',
        sha256='47dec4b1fb14bf6353628cc11183bddace1041fd0ee1cfcc22565628565042f0',
        schema='terluna.research.daily-life/1'),
    'research/studies/ways_of_living/results/ways_of_living.json': dict(
        branch='main', commit=PROVISIONING_COMMIT, blob='879dd4194f60547ba40ba5cecf2659dc978de2f8',
        sha256='9ef2b91e2d80edf9703036d011045f83655d10d4fccc39968fafa738c0d38148',
        schema='terluna.research.ways-of-living/1'),
}

# ---- Time and the Sun -----------------------------------------------------------------------------------------------
CYCLE_H = SYNODIC_MONTH_DAYS * JULIAN_DAY / 3600.0       # 708.73 h from noon to noon
NIGHT_H = CYCLE_H / 2                                     # the Sun's centre below the horizon, mean geometry
NIGHT_S = NIGHT_H * 3600.0
DEG_PER_H = 360.0 / CYCLE_H
SUN_RADIUS_DEG = math.degrees(math.asin(SUN_RADIUS / AU))
SAMPLES = 7200                                            # 0.098 h apart through one cycle
EARTH_SAMPLES = 720                                       # the Earth's light changes slowly; interpolated between


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_product(path, schema):
    data = json.loads(Path(path).read_text())
    if data.get('schema') != schema:
        raise ValueError(f'Unexpected schema for {path}: {data.get("schema")}')
    return data


def cross_branch_bytes(path):
    """A bound file: in-tree when its hash matches, as every pin does since the joint integration; else the Git
    object at a pinned commit."""
    spec = CROSS_BRANCH[path]
    local = ROOT / path
    if local.exists():
        raw = local.read_bytes()
        if hashlib.sha256(raw).hexdigest() == spec['sha256']:
            return raw, 'in-tree'
    try:
        raw = subprocess.run(['git', '-C', str(ROOT), 'show', f"{spec['commit']}:{path}"], check=True,
                             capture_output=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None, 'unavailable'
    if hashlib.sha256(raw).hexdigest() != spec['sha256']:
        return None, 'hash mismatch'
    return raw, 'git object at the pinned commit'


def cross_branch_json(path):
    raw, source = cross_branch_bytes(path)
    if raw is None:
        return None, source
    data = json.loads(raw)
    if data.get('schema') != CROSS_BRANCH[path]['schema']:
        raise ValueError(f'Unexpected schema for {path}')
    return data, source


def plain(x, digits=4):
    """JSON-ready copy: numpy scalars and arrays to Python, tuples to lists, floats to `digits` significant figures
    (None keeps them whole)."""
    if isinstance(x, dict):
        return {str(k): plain(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v, digits) for v in x]
    if isinstance(x, np.ndarray):
        return plain(x.tolist(), digits)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, (float, np.floating)):
        x = float(x)
        if not math.isfinite(x):
            raise ValueError('non-finite value in the product')
        return x if digits is None else float(f'{x:.{digits}g}')
    return x


# ==== Natural night light ============================================================================================
class Sky:
    """Clear-sky light on level ground in the mean geometry, through the calendar's evaluator (as daily_life uses it)."""

    def __init__(self, transfer):
        self.response = Response(transfer)
        self.elevation = np.linspace(-90.0, 90.0, 9001)
        base = dict(sun_radius_deg=SUN_RADIUS_DEG, sun_distance_factor=1.0, earth_phase_deg=180.0,
                    source_separation_deg=180.0, earth_radius_deg=0.0)
        parts = [self.response.light({**base, 'sun_elevation_deg': float(e)}, include_earth=False)
                 for e in self.elevation]
        self.solar = np.array([p['solar'] for p in parts])
        self.solar_diffuse = np.array([p['solar_diffuse'] for p in parts])
        if np.any(np.diff(self.solar) <= 0):
            raise ValueError("The Sun's light must rise with its elevation")

    def sun(self, elevation_deg):
        return np.interp(elevation_deg, self.elevation, self.solar)

    def sun_diffuse(self, elevation_deg):
        return np.interp(elevation_deg, self.elevation, self.solar_diffuse)

    def earth(self, lat, lon, hour_angle_deg):
        """The Earth's light (total and diffuse) for the Sun at the given hour angles from local noon."""
        n = hour_angle_deg.size
        sun_lon = np.radians(lon - hour_angle_deg)            # the subsolar point moves west through the cycle
        g = dict(earth=np.tile([1.0, 0.0, 0.0], (n, 1)), earth_distance_m=np.full(n, EARTH_MOON_DISTANCE),
                 sun=np.stack([np.cos(sun_lon), np.sin(sun_lon), np.zeros(n)], axis=-1), sun_distance_m=np.full(n, AU))
        geo = geometry_from_vectors(g, lon, lat)
        parts = [self.response.light({k: float(v[i]) for k, v in geo.items()}) for i in range(n)]
        total = np.array([p['earth'] for p in parts])
        diffuse = np.array([p['earth_diffuse'] for p in parts])
        return total, diffuse, float(geo['earth_elevation_deg'][0])


def hours_below(t, lux, level):
    """Hours per cycle below a level, by log-linear interpolation between samples (daily_life's crossing rule)."""
    y, line = np.log(np.maximum(lux, 1e-30)), math.log(level)
    y0, y1 = y, np.roll(y, -1)
    dt = t[1] - t[0]
    below = np.where((y0 < line) & (y1 < line), dt, 0.0)
    fall = (y0 >= line) & (y1 < line)
    rise = (y0 < line) & (y1 >= line)
    frac = np.where(fall | rise, (y0 - line) / np.where(y0 != y1, y0 - y1, 1.0), 0.0)
    below = below + np.where(fall, (1 - frac) * dt, 0.0) + np.where(rise, frac * dt, 0.0)
    return float(below.sum())


def night_light(sky, lat, lon, levels):
    """One place's natural light through the cycle: the night's darkest moment, its parts, and hours below levels."""
    hour_angle = np.arange(SAMPLES) * 360.0 / SAMPLES
    t = hour_angle / DEG_PER_H
    elevation = np.degrees(np.arcsin(np.cos(math.radians(lat)) * np.cos(np.radians(hour_angle))))
    sun, sun_diffuse = sky.sun(elevation), sky.sun_diffuse(elevation)
    coarse = np.arange(EARTH_SAMPLES) * 360.0 / EARTH_SAMPLES
    e_total, e_diffuse, earth_elevation = sky.earth(lat, lon, coarse)
    wrap = lambda y: np.interp(hour_angle, np.append(coarse, 360.0), np.append(y, y[0]))
    earth, earth_diffuse = wrap(e_total), wrap(e_diffuse)
    total = sun + earth
    night = (t >= 90.0 / DEG_PER_H) & (t < 270.0 / DEG_PER_H)
    i = int(np.argmin(np.where(night, total, np.inf)))
    mid = SAMPLES // 2
    return dict(lat=lat, lon=lon, earth_elevation_deg=earth_elevation,
                midnight_lux=float(total[mid]), darkest_lux=float(total[i]),
                darkest_hours_after_sunset=float(t[i] - 90.0 / DEG_PER_H),
                darkest_parts=dict(sun_diffuse_lux=float(sun_diffuse[i]), earth_lux=float(earth[i]),
                                   earth_diffuse_lux=float(earth_diffuse[i]),
                                   sun_elevation_deg=float(elevation[i])),
                hours_below={f'{x:g}': hours_below(t, total, x) for x in levels})


# ==== Skyglow in the Moon's air and Earth's =========================================================================
# A single-scattering model of the zenith brightness a ground source makes at a distance, after Garstang (1986): light
# leaves the source, scatters once toward the observer's zenith line, and is attenuated along both legs on a sphere.
# Only ratios between the two worlds and between thresholds are used, so the source's absolute flux cancels.
EARTH_SURFACE_DENSITY = 1.225             # kg/m3, U.S. Standard Atmosphere 1976
EARTH_SURFACE_PA = 101325.0               # Pa, U.S. Standard Atmosphere 1976
SKY_DISTANCES_M = np.geomspace(50.0, 800e3, 64)


def scattering_layers(world, column, haze=None):
    """Scattering layers (extinction per m at height z, phase, single-scattering albedo) for one world's air."""
    sigma = column['sigma_m2_per_kg']
    if world == 'moon':
        z_m = np.asarray(column['height_km'], float) * 1e3
        rho = np.asarray(column['density_kg_m3'], float)
        top_h = column['top_scale_height_m']

        def rho_at(z):
            inside = np.interp(z, z_m, rho)
            above = rho[-1] * np.exp(-(z - z_m[-1]) / top_h)
            return np.where(z <= z_m[-1], inside, above)
    else:
        h_e = column['scale_height_m']

        def rho_at(z):
            return EARTH_SURFACE_DENSITY * np.exp(-z / h_e)
    layers = [dict(beta=lambda z, f=rho_at: sigma * f(np.maximum(z, 0.0)), g=None, ssa=1.0)]
    if haze is not None and haze['aod'] > 0:
        b0, hh = haze['aod'] / haze['scale_height_m'], haze['scale_height_m']
        layers.append(dict(beta=lambda z, b0=b0, hh=hh: b0 * np.exp(-np.maximum(z, 0.0) / hh), g=haze['asymmetry'],
                           ssa=haze['single_scattering_albedo']))
    return layers


def zenith_skyglow(distances_m, radius_m, layers, top_m, pattern='lambert', n_h=500, n_path=160):
    """Zenith radiance at ground distances from a ground source of unit flux into the upper hemisphere (per sr).

    pattern 'lambert': light reflected from lit ground, intensity sin(elevation)/pi (fully shielded lamps);
    'uniform': light spread evenly over the upper hemisphere, 1/(2 pi) (lamps that also shine upward).
    """
    hs = np.geomspace(1.0, top_m, n_h)
    grid = np.concatenate([[0.0], hs])
    tau_out = 0.0                                             # the vertical leg, from each height down to the ground
    for L in layers:
        b = L['beta'](grid)
        tau_out = tau_out + np.cumsum(0.5 * (b[1:] + b[:-1]) * np.diff(grid))
    u = (np.arange(n_path) + 0.5) / n_path                    # midpoints along each slant path
    source = np.array([0.0, radius_m])
    out = []
    for d in distances_m:
        phi = d / radius_m
        zen = np.array([math.sin(phi), math.cos(phi)])
        v = (radius_m + hs)[:, None] * zen[None, :] - source[None, :]
        s = np.linalg.norm(v, axis=1)
        up = v[:, 1] / s                                      # sine of the ray's elevation at the source
        pts = source[None, None, :] + u[None, :, None] * v[:, None, :]
        alt = np.linalg.norm(pts, axis=2) - radius_m
        blocked = (up <= 0) | np.any(alt < -1.0, axis=1)
        tau_in = sum(L['beta'](alt).sum(axis=1) * s / n_path for L in layers)
        cospsi = -(v @ zen) / s
        scatter = 0.0
        for L in layers:
            if L['g'] is None:
                phase = 0.75 * (1 + cospsi ** 2) / (4 * math.pi)
            else:
                g = L['g']
                phase = (1 - g * g) / (4 * math.pi * (1 + g * g - 2 * g * cospsi) ** 1.5)
            scatter = scatter + L['ssa'] * L['beta'](hs) * phase
        intensity = up / math.pi if pattern == 'lambert' else np.full_like(up, 1 / (2 * math.pi))
        integrand = np.where(blocked, 0.0, intensity / s ** 2 * np.exp(-tau_in - tau_out) * scatter)
        out.append(float(np.trapezoid(integrand, hs)))
    return np.array(out)


# ==== The settings: places, latitude rows and the ways-of-living homes ===============================================
# Latitude bands of the ways-of-living homes (their bands_dry_land_km2), each read at the middles of five equal-area
# slices; the never-dark latitude closes the middle band.
BANDS = dict(equatorial=(0.0, 15.0), tropical=(15.0, 30.0), middle=(30.0, 48.56), high=(48.56, 90.0))
# Where the Earth stands, by the homes' earth_up_shares (ways of living: high is 30 degrees and more): read at these
# Earth elevations; the limb's libration is read with the Earth on the horizon, the far side on the 180 meridian.
EARTH_CLASSES = dict(high=(30.0, 50.0, 70.0), low=(12.0, 20.0), rising_and_setting=(0.0,), never=None)
HOMES = ('summit_metropolis', 'coastal_and_lake_towns', 'wetland_edge_towns', 'canopy_towns', 'farming_country',
         'highland_towns', 'elevated_districts', 'nearside_towns', 'limb_towns', 'far_side_places', 'high_latitude_towns')


BAND_SLICES = 5


def band_latitudes(lo, hi, n=BAND_SLICES):
    """The latitudes at the middle of n equal-area slices of a band (area uniform in sin latitude)."""
    a, b = math.sin(math.radians(lo)), math.sin(math.radians(hi))
    return [math.degrees(math.asin(a + (k + 0.5) / n * (b - a))) for k in range(n)]


def earth_longitude(lat, elevation):
    """East longitude at which the mean-geometry Earth (over 0, 0) stands at this elevation, at this latitude."""
    c = math.sin(math.radians(elevation)) / math.cos(math.radians(lat))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def zenith_ratio_table(solved):
    """Zenith luminance over the sky's mean luminance (diffuse horizontal / pi) by source elevation, Moon's solved sky."""
    samples = sorted(solved['worlds']['moon_1.2atm']['samples'], key=lambda s: s['sun_deg'])
    e = np.array([s['sun_deg'] for s in samples])
    k = np.array([s['zenith_luminance_cd_m2'] / (s['diffuse_horizontal_lux'] / math.pi) for s in samples])
    return e, k


def natural_zenith_cd_m2(parts, earth_elevation, ratio, twilight_k, floor_cd_m2):
    """The natural sky's zenith luminance at a night's darkest moment.

    Twilight is all diffuse light; its zenith share of the sky's mean luminance is twilight_k (the solved sky gives
    0.65 at the deepest sampled depression, 30 degrees). Earthlight's diffuse part takes the solved sky's ratio at the
    Earth's elevation. floor_cd_m2 stands for starlight, zodiacal light and airglow, which the transfer product omits.
    """
    e, k = ratio
    earth_k = float(np.interp(max(earth_elevation, 0.0), e, k))
    return (twilight_k * parts['sun_diffuse_lux'] / math.pi + earth_k * parts['earth_diffuse_lux'] / math.pi
            + floor_cd_m2)


def disk_profile(distances_m, brightness, radius_m, xs_m, n_rho=1200):
    """Zenith brightness at distances xs from the centre of a uniformly lit disk town of the given radius, per unit of
    its total flux, from a point source's curve (flat geometry across the town; the curve carries the curvature).
    Below the curve's first distance the brightness is carried as 1/distance, its near-field form."""
    d0, b0 = distances_m[0], brightness[0]

    def point(rho):
        inner = b0 * d0 / np.maximum(rho, 1e-3)
        outer = np.exp(np.interp(np.log(np.maximum(rho, d0)), np.log(distances_m), np.log(brightness)))
        return np.where(rho < d0, inner, outer)

    out = []
    for x in xs_m:
        rho = np.concatenate([np.geomspace(1.0, max(2.0, x + radius_m), n_rho)])
        if x <= 1e-9:
            arc = np.where(rho <= radius_m, 2 * math.pi * rho, 0.0)
        else:
            inside = rho <= radius_m - x
            cross = (rho > abs(radius_m - x)) & (rho < radius_m + x)
            c = np.clip((rho ** 2 + x ** 2 - radius_m ** 2) / (2 * rho * x), -1.0, 1.0)
            arc = np.where(inside, 2 * math.pi * rho, np.where(cross, 2 * rho * np.arccos(c), 0.0))
        out.append(float(np.trapezoid(point(rho) * arc, rho)) / (math.pi * radius_m ** 2))
    return np.array(out)


def polluted_radius(xs_m, profile, level):
    """The distance from a town's centre out to which its skyglow exceeds a level (0 when it never does)."""
    if level >= profile[0]:
        return 0.0
    tail = np.maximum(profile[1:], 1e-300)
    if level <= tail[-1]:
        return float(xs_m[-1])
    return float(np.exp(np.interp(-math.log(level), -np.log(tail), np.log(xs_m[1:]))))


# ==== Photometry =====================================================================================================
# CIE 1924 photopic V(lambda) (the CIE 1931 2-degree y-bar) and CIE 1951 scotopic V'(lambda), 380-780 nm at 5 nm, as
# tabulated by the Colour & Vision Research Laboratory (CVRL) and distributed with HDR-VDP 3.0.7 (sources.json).
LAMBDA_NM = np.arange(380.0, 781.0, 5.0)
V_PHOTOPIC = np.array([
    3.9e-05, 6.4e-05, 0.00012, 0.000217, 0.000396, 0.00064, 0.00121, 0.00218, 0.004, 0.0073, 0.0116, 0.01684, 0.023,
    0.0298, 0.038, 0.048, 0.06, 0.0739, 0.09098, 0.1126, 0.13902, 0.1693, 0.20802, 0.2586, 0.323, 0.4073, 0.503,
    0.6082, 0.71, 0.7932, 0.862, 0.91485, 0.954, 0.9803, 0.99495, 1, 0.995, 0.9786, 0.952, 0.9154, 0.87, 0.8163,
    0.757, 0.6949, 0.631, 0.5668, 0.503, 0.4412, 0.381, 0.321, 0.265, 0.217, 0.175, 0.1382, 0.107, 0.0816, 0.061,
    0.04458, 0.032, 0.0232, 0.017, 0.01192, 0.00821, 0.005723, 0.004102, 0.002929, 0.002091, 0.001484, 0.001047,
    0.00074, 0.00052, 0.0003611, 0.0002492, 0.0001719, 0.00012, 8.48e-05, 6e-05, 4.24e-05, 3e-05, 2.12e-05,
    1.499e-05])
V_SCOTOPIC = np.array([
    0.000589, 0.001108, 0.002209, 0.00453, 0.00929, 0.01852, 0.03484, 0.0604, 0.0966, 0.1436, 0.1998, 0.2625, 0.3281,
    0.3931, 0.455, 0.513, 0.567, 0.62, 0.676, 0.734, 0.793, 0.851, 0.904, 0.949, 0.982, 0.998, 0.997, 0.975, 0.935,
    0.88, 0.811, 0.733, 0.65, 0.564, 0.481, 0.402, 0.3288, 0.2639, 0.2076, 0.1602, 0.1212, 0.0899, 0.0655, 0.0469,
    0.03315, 0.02312, 0.01593, 0.01088, 0.00737, 0.00497, 0.003335, 0.002235, 0.001497, 0.001005, 0.000677, 0.000459,
    0.0003129, 0.0002146, 0.000148, 0.0001026, 7.15e-05, 5.01e-05, 3.533e-05, 2.501e-05, 1.78e-05, 1.273e-05,
    9.14e-06, 6.6e-06, 4.78e-06, 3.482e-06, 2.546e-06, 1.87e-06, 1.379e-06, 1.022e-06, 7.6e-07, 5.67e-07, 4.25e-07,
    3.196e-07, 2.413e-07, 1.829e-07, 1.39e-07])
KM_PHOTOPIC = 683.0                       # lm/W at 540 THz (555 nm), the SI definition of the candela
KM_SCOTOPIC = 1700.0                      # lm/W at 507 nm, CIE scotopic maximum
GREY = 0.18                               # the 18% grey surface the light calendar's modes of vision use
MESOPIC_CD_M2 = (0.005, 5.0)              # CIE 191:2010, the mesopic range of adaptation luminance


def photon_spectrum(peak_nm, fwhm_nm):
    """Photon spectrum on the CIE grid, normalised to one photon: a Gaussian in wavelength (an approximation)."""
    sd = fwhm_nm / (2 * math.sqrt(2 * math.log(2)))
    w = np.exp(-0.5 * ((LAMBDA_NM - peak_nm) / sd) ** 2)
    return w / w.sum()


def per_photon(peak_nm, fwhm_nm):
    """Mean photon energy (J) and photopic and scotopic lumens carried per photon per second."""
    n = photon_spectrum(peak_nm, fwhm_nm)
    energy = PLANCK * SPEED_OF_LIGHT / (LAMBDA_NM * 1e-9)
    return dict(photon_j=float(n @ energy), lm_s_photopic=float(KM_PHOTOPIC * (n * energy) @ V_PHOTOPIC),
                lm_s_scotopic=float(KM_SCOTOPIC * (n * energy) @ V_SCOTOPIC))


def photons_for_lux(lux, peak_nm, fwhm_nm, scotopic=False):
    """Photon irradiance (photons per second per m2) that gives an illuminance, photopic or scotopic lux."""
    p = per_photon(peak_nm, fwhm_nm)
    return lux / (p['lm_s_scotopic'] if scotopic else p['lm_s_photopic'])


def lambert_luminance(photons_m2_s, peak_nm, fwhm_nm):
    """Photopic and scotopic luminance (cd/m2) of a Lambertian surface emitting this photon exitance."""
    p = per_photon(peak_nm, fwhm_nm)
    return dict(photopic=photons_m2_s * p['lm_s_photopic'] / math.pi, scotopic=photons_m2_s * p['lm_s_scotopic'] / math.pi)


# ==== The cost of living light =======================================================================================
C_G_PER_MOL = 12.011                      # g per mol of carbon, its conventional atomic weight
GLUCOSE_KJ_PER_G_C = 15.7 * 180.156 / 72.064   # FAO (2003)'s 15.7 kJ per g of monosaccharide (as the lunar-cycle
#                                         ecology study records it), per g of its carbon: 39.2 kJ
# ATP equivalents of each cofactor, from Hinkle (2005)'s P/O ratios of about 2.5 for NADH and 1.5 for succinate
# (2.3 and 1.4 if ATP synthase moves 10/3 protons): NAD(P)H counts 2.3-2.5, and an acetyl unit burnt in the citric
# acid cycle 3 NADH + 1 FADH2 + 1 GTP; respiring glucose gives 10 NADH + 2 FADH2 + 4 ATP and GTP (derived).
P_O_NADH = (2.3, 2.5)
P_O_FADH2 = (1.4, 1.5)


def atp_values(i):
    """ATP equivalents of the cofactors and carbon units, at the low (0) or high (1) P/O ratios."""
    nadh, fadh = P_O_NADH[i], P_O_FADH2[i]
    acetyl = 3 * nadh + fadh + 1
    pyruvate = acetyl + nadh
    glucose = 10 * nadh + 2 * fadh + 4
    return dict(ATP=1.0, ATP_to_AMP=2.0, NADPH=nadh, NADH=nadh, acetyl_CoA=acetyl, pyruvate=pyruvate,
                glucose=glucose, per_mol_C=glucose / 6)


H3H_UNCOUPLING = 0.3                     # NADPH spent without hydroxylating hispidin (Tong et al. 2020)
# What one turnover of each luciferin cycle consumes and returns (sources.json gives each enzyme's reaction). Positive
# counts are consumed, negative returned. The firefly's luciferin is regenerated from oxyluciferin with one cysteine
# (Gomi & Kajiyama 2001); cysteine from serine and sulfate is costed below by its pathway.
PATHWAYS = {
    'fungal': dict(  # caffeic acid cycle (Kotlobay et al. 2018): HispS, H3H, Luz, CPH; ACC makes the malonyl-CoA
        consumed=dict(ATP_to_AMP=1, ATP=4, NADPH=1 / (1 - H3H_UNCOUPLING), acetyl_CoA=2), returned=dict(pyruvate=1),
        note='HispS adenylates caffeic acid and adds two malonyl-CoA (two ATP at acetyl-CoA carboxylase, two more at '
             'ATP-citrate lyase for the cytosol\'s acetyl-CoA); H3H uses '
             'NADPH and O2, wasting 30% of the NADPH; Luz releases CO2 and light; CPH returns caffeic acid and pyruvate'),
    'firefly': dict(  # luciferase and luciferin regeneration
        consumed=dict(ATP_to_AMP=1, cysteine=1), returned=dict(),
        note='luciferase spends ATP to AMP and pyrophosphate; oxyluciferin returns to luciferin through '
             '2-cyano-6-hydroxybenzothiazole and one cysteine'),
    'bacterial': dict(  # lux: luciferase, flavin reductase and the LuxCDE fatty-acid reductase
        consumed=dict(NADPH=2, ATP_to_AMP=1), returned=dict(),
        note='FMNH2 from NAD(P)H; the fatty acid back to aldehyde by LuxCDE with ATP to AMP and NADPH'),
}


def cysteine_atp(v):
    """Cysteine from 3-phosphoglycerate and sulfate (derived from the pathway's stoichiometry): sulfate activated by ATP
    to AMP, eight electrons to sulfide (four NAD(P)H), the acetyl of O-acetylserine lost as acetate and reactivated by
    ATP to AMP, serine's amino group assimilated through glutamine synthetase and glutamate synthase (one ATP and one
    NAD(P)H), and the serine skeleton valued as the 3-phosphoglycerate it diverts (one ATP and a pyruvate) less the
    NADH made on the way to serine. The two carbons that leave as thioglycolate are not credited back."""
    return (v['ATP_to_AMP'] + 4 * v['NADPH'] + v['ATP_to_AMP'] + (v['ATP'] + v['NADPH'])
            + (1 + v['pyruvate'] - v['NADH']))


def turnover_cost(pathway, i):
    """ATP equivalents consumed by one luciferin turnover, at the low (0) or high (1) P/O ratios."""
    v = atp_values(i)
    v = dict(v, cysteine=cysteine_atp(v))
    spec = PATHWAYS[pathway]
    return sum(n * v[k] for k, n in spec['consumed'].items()) - sum(n * v[k] for k, n in spec['returned'].items())


# ==== People's light and the polluted sky ==========================================================================
# Ecology's people screen puts 791-1,420 m2 of land a person under a light-polluted sky with the lunar design's
# shielded lamps (a tenth of Earth's 7,910-14,200 m2 a person in Europe and the US; Falchi et al. 2016 over the World
# Bank's land and people). Those areas are Earth's: Earth's air and Earth's natural night sky. Here each is carried to
# the Moon as a town of N people whose polluted zone on Earth covers N times the area a person, read through the
# single-scattering model in both worlds' air, against each setting's natural sky, with the town's light scaled by
# twilight-aware lamps and by living light.
TOWN_PEOPLE = (1e4, 1e5, 1e6)            # towns read; the metropolis is read on its own at its 100 million
METROPOLIS_PEOPLE = 1e8                  # the summit metropolis brief's about 100 million
SETTLED_M2 = dict(town=100.0, metropolis=25.0)   # ways of living: a town at 10,000 a km2, the metropolis at 40,000
XS_M = np.concatenate([[0.0], np.geomspace(20.0, 790e3, 160)])


def grid_points():
    """(band, latitude, class, Earth elevation, longitude, weight in the class) for every reading of the homes' grid.

    Within a class the readings are weighted by the area each stands for, which goes as the cosine of the Earth's
    elevation."""
    points = []
    for band, (lo, hi) in BANDS.items():
        for lat in band_latitudes(lo, hi):
            for cls, elevations in EARTH_CLASSES.items():
                if elevations is None:
                    points.append((band, lat, cls, None, 180.0, 1.0))
                    continue
                uses = []
                for e in elevations:
                    e_use = min(e, 90.0 - lat - 1.0)
                    if cls == 'high' and e_use < 30.0:
                        e_use = min(30.0, 90.0 - lat - 1.0)
                    uses.append(e_use)
                norm = sum(math.cos(math.radians(e)) for e in uses)
                for e_use in uses:
                    points.append((band, lat, cls, e_use, earth_longitude(lat, e_use),
                                   math.cos(math.radians(e_use)) / norm))
    return points


# ==== People across the settings =====================================================================================
LUNAR_PEOPLE_BILLION = (4.0, 6.0, 8.0)   # the author's working range on the Moon, 2026-10-11 (the provisioning
#                                          branch's decisions register at b82e0dc, People on the Moon and beyond)
ALOFT_SHARE = 0.5                        # 2-4 billion of them aloft (the same register): half
# The author's direction of the same day, land residents "mostly in the drier highlands, the summit metropolis round its
# tower above all", read the way provisioning builds its mixes: four times as dense on the highland homes (a guess at
# the mapping; the distribution is the author's to choose).
DERIVED_MIXES = dict(highlands=dict(base='spread', weight=4.0, homes=('highland_towns', 'elevated_districts')))


def mix_shares(ways, mix):
    """Each surface home's share of the people outside the metropolis and the undersea domes, by a provisioning mix
    (ways of living: the total25 case's all-surface split, which fixes the metropolis at its planned 100 million), or by
    a mix derived from one of them."""
    spec = DERIVED_MIXES.get(mix)
    people = ways['cases']['total25']['by_aerial_share']['0']['surface'][spec['base'] if spec else mix]['people_billion']
    if spec:
        people = {k: v * (spec['weight'] if k in spec['homes'] else 1.0) for k, v in people.items()}
    free = {k: v for k, v in people.items() if k not in ('summit_metropolis', 'undersea_communities')}
    total = sum(free.values())
    return {k: v / total for k, v in free.items()}, people['summit_metropolis']


# ==== Living light: the measured systems ==============================================================================
# Emission peaks: fungal luciferases 520-530 nm (Kotlobay et al. 2018; Perry et al. 2025), Photinus pyralis luciferase
# 560-566 nm (Niwa et al. 2010, as Niwa et al. 2023 tabulate), bacterial luciferase about 490 nm (Brodl et al. 2018),
# the sea's luminous plankton about 475 nm (as Brodl et al. 2018 give for marine bioluminescence). The widths are
# assumptions for the Gaussian stand-ins.
CHEMISTRIES = {
    'fungal pathway': dict(pathway='fungal', quantum_yield=(0.00096, 0.017), peak_nm=525.0, fwhm_nm=80.0,
                           basis='Mycena luciferases with 3-hydroxyhispidin: 0.017 (gill) and 0.00096 (stipe), '
                                 'Teranishi 2018; the fungal quantum yield is not yet established (Kusuma et al. 2025)'),
    'firefly chemistry': dict(pathway='firefly', quantum_yield=(0.41, 0.48), peak_nm=560.0, fwhm_nm=70.0,
                              basis='Photinus pyralis luciferase: 0.41 +/- 0.074 (Ando et al. 2008); 0.48 (Niwa et '
                                    'al. 2010, as Niwa et al. 2023 report)'),
    'bacterial chemistry': dict(pathway='bacterial', quantum_yield=(0.10, 0.17), peak_nm=490.0, fwhm_nm=80.0,
                                basis='0.17 photons per aldehyde (Shimomura, Johnson & Kohama 1972); about 0.1-0.16 '
                                      '(Kaku et al. 2021, citing Hastings & Nealson 1977)'),
    'thermodynamic floor': dict(pathway='floor', quantum_yield=(1.0,), peak_nm=525.0, fwhm_nm=80.0,
                                basis='one photon per luciferin and nothing lost'),
}
BEST_BEETLE_QY = 0.61                    # Pyrearinus termitilluminans luciferase at 539 nm, the highest beetle quantum
#                                          yield measured (Niwa et al. 2010, as Niwa et al. 2023 tabulate)
PER_CM2_MIN = 1e4 / 60.0                 # photons per minute per cm2 to photons per second per m2
EMITTERS = [
    dict(key='Nicotiana benthamiana flowers, 2020', kind='engineered plant',
         organism='N. benthamiana carrying the fungal pathway (Mitiouchkina et al. 2020, Supplementary Table 1, mean '
                  'of 5)', photons_m2_s=6.47e10 * PER_CM2_MIN, peak_nm=525.0, fwhm_nm=80.0, chemistry='fungal',
         quantum_yield=CHEMISTRIES['fungal pathway']['quantum_yield'], source='mitiouchkina_2020', measured=True,
         derivation='6.47e10 photons a minute per cm2'),
    dict(key='Nicotiana tabacum flowers, 2020', kind='engineered plant',
         organism='N. tabacum carrying the fungal pathway (Mitiouchkina et al. 2020, Supplementary Table 1, mean of 5)',
         photons_m2_s=3.17e10 * PER_CM2_MIN, peak_nm=525.0, fwhm_nm=80.0, chemistry='fungal',
         quantum_yield=CHEMISTRIES['fungal pathway']['quantum_yield'], source='mitiouchkina_2020', measured=True,
         derivation='3.17e10 photons a minute per cm2'),
    dict(key='Nicotiana leaves, 2020', kind='engineered plant',
         organism='N. benthamiana leaves (Mitiouchkina et al. 2020, Supplementary Table 1, mean of 19); tobacco '
                  'leaves give 4.31e9', photons_m2_s=5.58e9 * PER_CM2_MIN, peak_nm=525.0, fwhm_nm=80.0,
         chemistry='fungal', quantum_yield=CHEMISTRIES['fungal pathway']['quantum_yield'], source='mitiouchkina_2020', measured=True,
         derivation='5.58e9 photons a minute per cm2'),
    dict(key='enhanced tobacco flowers, 2023', kind='engineered plant',
         organism='tobacco with the enhanced fungal pathway, brightest flowers (Zheng et al. 2023)',
         photons_m2_s=3e11 * PER_CM2_MIN, peak_nm=525.0, fwhm_nm=80.0, chemistry='fungal',
         quantum_yield=CHEMISTRIES['fungal pathway']['quantum_yield'], source='zheng_2023', measured=True,
         derivation='a maximum of 3e11 photons a minute per cm2'),
    dict(key='improved pathway, 2024, tenfold', kind='engineered plant',
         organism='the improved pathway of Shakhova et al. (2024), at the low end of its one to two orders of magnitude '
                  'over the 2020 pathway, applied to the 2020 N. benthamiana flowers (derived)',
         photons_m2_s=10 * 6.47e10 * PER_CM2_MIN, peak_nm=525.0, fwhm_nm=80.0, chemistry='fungal',
         quantum_yield=CHEMISTRIES['fungal pathway']['quantum_yield'], source='shakhova_2024', measured=False,
         derivation='ten times 6.47e10 photons a minute per cm2'),
    dict(key='improved pathway, 2024, hundredfold', kind='engineered plant',
         organism='the same at the high end, a hundred times the 2020 flowers (derived)',
         photons_m2_s=100 * 6.47e10 * PER_CM2_MIN, peak_nm=525.0, fwhm_nm=80.0, chemistry='fungal',
         quantum_yield=CHEMISTRIES['fungal pathway']['quantum_yield'], source='shakhova_2024', measured=False,
         derivation='a hundred times 6.47e10 photons a minute per cm2'),
    dict(key='milky sea', kind='natural',
         organism='luminous bacteria at the sea surface, at the satellite\'s adjusted detection threshold (Miller et al. '
                  '2005); the brightest patches seen since approached crescent- to quarter-moon light (Miller et al. 2021)',
         photons_m2_s=1.4e11 * 1e4, peak_nm=490.0, fwhm_nm=80.0, chemistry='bacterial',
         quantum_yield=CHEMISTRIES['bacterial chemistry']['quantum_yield'], source='miller_2005', measured=True,
         derivation='1.4e11 photons per cm2 a second (1.8e-4 W/m2/sr)'),
]
# Per-organism outputs that are not surfaces: what it takes to glow at a given luminance, and a flash's store.
FLASHES = {
    'Pyrodinium bahamense': dict(peak_photons_s=2.93e8, photons_per_flash=3.35e8, peak_nm=475.0, fwhm_nm=40.0,
                                 sources=['letendre_2025', 'letendre_2024']),
    'Lingulodinium polyedra': dict(peak_photons_s=1.11e8, photons_per_flash=1.17e8, peak_nm=475.0, fwhm_nm=40.0,
                                   sources=['letendre_2025', 'letendre_2024']),
    'Pyrocystis fusiformis': dict(peak_photons_s=1.59e10, photons_per_flash=None, peak_nm=475.0, fwhm_nm=40.0,
                                  sources=['letendre_2025']),
    'Noctiluca scintillans': dict(peak_photons_s=9.0e10, photons_per_flash=2.5e11, peak_nm=475.0, fwhm_nm=40.0,
                                  sources=['letendre_2024']),
}
BACTERIA_PHOTONS_S = (1e3, 2.0e4)        # a cell in a squid's light organ (about 1e3; Perez et al. 2010) to
#                                          Photobacterium leiognathi at its brightest (2.0e4; Dunlap 1985)
FIREFLY_FLASH = dict(photons=(1e8, 1e11), lantern_mm2=(0.74, 3.52), sources=['silver_2026', 'lopez_palafox_2020'])
ENZYME_PHOTONS_S = dict(NanoLuc=8.6, firefly_luciferase=0.2)   # photons a second per molecule at saturating substrate
#                                                                (Klein et al. 2023)
TAXA_USED = ('G28', 'M09', 'G39', 'C07', 'C12')


# ==== What people need outdoors at night =============================================================================
# Levels are horizontal illuminance on the ground. A face takes semi-cylindrical illuminance; under a uniform night
# sky it gets half the ground's light (derived), so the face levels are doubled to the ground.
SEMICYLINDRICAL_OVER_HORIZONTAL = 0.5
NEEDS = {
    'orientation': dict(lux=0.05, task='moving about open ground by eye, as under a summer full moon on Earth',
                        basis='full moonlight on Earth is typically 0.05-0.1 lux in temperate summers, 0.05-0.2 lux '
                              'at most full moons and at most about 0.3 lux (Kyba, Mohar & Posch 2017)'),
    'colour': dict(lux=0.32, task='telling colours',
                   basis='paper samples were sorted into their colour groups above 0.32 lux; at 0.1-0.01 lux only reds '
                         'and oranges were named reliably (Pokorny et al. 2006)'),
    'walking': dict(lux=1.0, task='walking safely: seeing trip hazards',
                    basis='1.0 lux lets pedestrians of all ages detect a 10-mm trip hazard under any lamp, 0.22-0.93 lux '
                          'by lamp and age (Fotios & Uttley 2018); detection reaches its plateau at 2.0 lux (Uttley et '
                          'al. 2017); the lowest pedestrian class, P6, asks 2.0 lux on average and 0.4 lux at least '
                          '(CIE 115:2010, as Fotios & Gibbons 2018 tabulate it)'),
    'faces at 4 m': dict(lux=0.8 / SEMICYLINDRICAL_OVER_HORIZONTAL, semicylindrical_lux=0.8,
                         task='recognising a face at 4 m',
                         basis='semi-cylindrical illuminance of 0.8 lux at 4 m and 2.7 lux at 10 m (Caminada & van '
                               'Bommel 1980, as Fotios, Yang & Cheal 2015 report it)'),
    'faces at 10 m': dict(lux=2.7 / SEMICYLINDRICAL_OVER_HORIZONTAL, semicylindrical_lux=2.7,
                          task='recognising a face at 10 m', basis='as for faces at 4 m'),
    'reading': dict(lux=8.0, task='reading ordinary print',
                    basis='at a mesopic 2 cd/m2 people read at full speed with a larger critical print size (Penaloza et '
                          'al. 2025): about 8 lux on paper of reflectance 0.8 (derived)'),
    'work': dict(lux=100.0, task='close work and reading small print outdoors',
                 basis='100 lux for reading addresses outdoors and 500 lux for reading and writing indoors (EN 12464-2 '
                       'and EN 12464-1, as the Zumtobel lighting handbook reprints them)'),
}
OUTDOOR_LUX = NEEDS['walking']['lux']   # twilight-aware outdoor lamps add only what the natural light lacks of this


# ==== Geometry: what a glowing surface gives a path ==================================================================
def hedge_view_factor(half_width_m, height_m):
    """Configuration factor from a point on a path's centre line to two long glowing hedges, one on each side.

    For a horizontal element and a long strip seen between angles phi1 and phi2 from the zenith in the plane across
    the path, the factor is (sin phi2 - sin phi1) / 2; a hedge runs from the horizon (90 degrees) to its top."""
    top = math.atan2(half_width_m, height_m)
    return 2 * (1 - math.sin(top)) / 2


GEOMETRIES = {
    'closed canopy overhead': dict(view_factor=1.0, note='a glowing underside that fills the sky seen from the path'),
    'canopy with a third open': dict(view_factor=2 / 3, note='an arbour or forest edge with a third of the sky open'),
    'hedges 1 m tall on a 2 m path': dict(view_factor=hedge_view_factor(1.0, 1.0), note='both sides glowing'),
    'hedges 2 m tall on a 3 m path': dict(view_factor=hedge_view_factor(1.5, 2.0), note='both sides glowing'),
}


def glow_cost(photons_m2_s, quantum_yield, pathway, hours, i, delivered=1.0):
    """ATP equivalents and carbon (g C per m2) for a glow held for a number of hours (low or high P/O ratios, i), when
    a share `delivered` of the photons made leaves the tissue toward where it is wanted."""
    v = atp_values(i)
    turnovers = photons_m2_s / quantum_yield / delivered
    atp_per_photon = turnover_cost(pathway, i) / quantum_yield / delivered
    mol_c_per_mol_photons = atp_per_photon / v['per_mol_C']
    mol_photons = photons_m2_s * hours * 3600.0 / AVOGADRO
    return dict(atp_per_photon=atp_per_photon, mol_c_per_mol_photons=mol_c_per_mol_photons,
                turnovers_m2_s=turnovers, g_c_m2=mol_photons * mol_c_per_mol_photons * C_G_PER_MOL)


def floor_mol_c_per_mol_photons(photon_j):
    """The thermodynamic floor: one photon per luciferin and the photon's energy paid at the energy of respiring
    carbohydrate, with nothing lost."""
    return photon_j * AVOGADRO / (GLUCOSE_KJ_PER_G_C * 1e3 * C_G_PER_MOL)


# Light levels at which hours below are counted, lux on level ground (the needs below read from these).
LEVELS_LUX = tuple(sorted({0.001, 0.003, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0,
                           10.0, 20.0, 50.0, 100.0, 200.0, 500.0} | {v['lux'] for v in NEEDS.values()}))


# ==== Assembling the natural light ===================================================================================
def moon_class_areas(n=721):
    """Each latitude band's share of the Moon's area, and within it the share of each Earth class (mean geometry)."""
    out = {}
    lon = np.radians(np.linspace(-180.0, 180.0, n)[:-1])
    for band, (lo, hi) in BANDS.items():
        lat = np.radians(np.linspace(lo, hi, 181))
        w = np.cos(lat)
        e = np.degrees(np.arcsin(np.cos(lat)[:, None] * np.cos(lon)[None, :]))
        cls = dict(high=e >= 30.0, low=(e >= 8.0) & (e < 30.0), rising_and_setting=np.abs(e) < 8.0, never=e <= -8.0)
        total = (w[:, None] * np.ones_like(e)).sum()
        out[band] = dict(area_share=(math.sin(math.radians(hi)) - math.sin(math.radians(lo))),
                         classes={k: float((w[:, None] * v).sum() / total) for k, v in cls.items()})
    return out


def natural_light(sky, daily, places):
    """The places and latitude rows of the daily-life study recomputed, the homes' grid, and the checks."""
    rows = {}
    worst = dict(lux_beyond_print=0.0, hours_beyond_print=0.0, note='largest deviation from the daily-life calendars '
                 'beyond their printed precision (lux to 4 decimals, hours to 2)')
    lv = daily['levels']
    for name, c in daily['light']['calendars'].items():
        r = night_light(sky, c['lat'], c['lon'], LEVELS_LUX + (lv['dusk'], lv['dark'], lv['scotopic']))
        hb = r['hours_below']
        for ours, theirs in ((r['darkest_lux'], c['darkest_lux']), (r['midnight_lux'], c['midnight_lux'])):
            worst['lux_beyond_print'] = max(worst['lux_beyond_print'], abs(ours - theirs) - 5e-5)
        for lvl, field in ((lv['dusk'], 'below_dusk_h'), (lv['dark'], 'below_dark_h'), (lv['scotopic'], 'below_scotopic_h')):
            worst['hours_beyond_print'] = max(worst['hours_beyond_print'], abs(hb[f'{lvl:g}'] - c[field]) - 5e-3)
        rows[name] = dict(r, hours_below={k: v for k, v in hb.items() if float(k) in LEVELS_LUX})
    main = {}
    for name, p in places['places'].items():
        r = night_light(sky, p['lat'], p['lon'], (0.1,))
        main[name] = dict(main_darkest_lux=p['night']['natural_darkest_lux'], darkest_lux=r['darkest_lux'],
                          main_midnight_lux=p['night']['natural_midnight_lux'], midnight_lux=r['midnight_lux'],
                          main_hours_below_0_1_lux=p['night']['hours_below_0_1_lux_natural'],
                          hours_below_0_1_lux=r['hours_below']['0.1'])
    grid, cache = [], {}
    for band, lat, cls, e, lon, w_class in grid_points():
        key = (round(lat, 6), round(lon, 6))
        if key not in cache:
            cache[key] = night_light(sky, lat, lon, LEVELS_LUX)
        grid.append(dict(cache[key], band=band, earth_class=cls, earth_elevation_target_deg=e, class_weight=w_class))
    return dict(places=rows, main_places=main, grid=grid, checks_against_daily_life=worst,
                classes_by_band=moon_class_areas())


# ==== Assembling the skyglow =========================================================================================
def columns(transfer, winds):
    """The two worlds' air for the skyglow model: the design run's global-mean density profile for the Moon, an
    exponential atmosphere with the 1976 standard's surface density and pressure for Earth, both scattering at the
    mass cross-section of the solved 1.2-atm column at 550 nm."""
    atm = transfer['atmosphere']['column']
    sigma = atm['rayleigh_optical_depth_550nm'] / atm['air_column_kg_m2']
    z, rho = np.array(winds['height_km'], float), np.array(winds['density_mean_kg_m3'], float)
    top_h = (z[-1] - z[-2]) * 1e3 / math.log(rho[-2] / rho[-1])
    moon = dict(sigma_m2_per_kg=sigma, height_km=z, density_kg_m3=rho, top_scale_height_m=top_h)
    earth = dict(sigma_m2_per_kg=sigma, scale_height_m=EARTH_SURFACE_PA / STANDARD_GRAVITY / EARTH_SURFACE_DENSITY)
    return moon, earth


def column_tau(layers, top_m):
    z = np.geomspace(1.0, top_m, 6000)
    return float(sum(np.trapezoid(L['beta'](z), z) for L in layers))


def local_slope(distances, curve, lo, hi):
    """d log B / d log D between two distances."""
    a, b = (np.exp(np.interp(math.log(x), np.log(distances), np.log(curve))) for x in (lo, hi))
    return math.log(b / a) / math.log(hi / lo)


# Design guesses for people's light (ranges, labelled once here and in the README):
DESIGN_GUESSES = dict(
    residual_share=(0.1, 0.5),          # upward light from windows, vehicles and task lights that no outdoor rule
    #                                     governs; street lights give 13% of a city's light after midnight on Earth
    #                                     (Kyba et al. 2021), and the lunar design governs every outdoor lamp
    living_light_share=(0.3, 0.7),      # share of the outdoor lighting service that living light takes over
    living_light_uplight=(0.01, 0.1),   # living light's upward flux per unit of that service, against shielded lamps
    delivered_share=(0.1, 0.9),         # share of the photons a glowing plant makes that reach the path: up to 0.9 from
    #                                     glowing cells on leaves' undersides beneath a reflecting layer, 0.1 from glow
    #                                     made inside green leaf tissue, which absorbs it and shines it upward as well
)
# Assumptions about the natural sky (ranges):
TWILIGHT_K = (0.4, 0.8)                 # the zenith's share of the night sky's mean luminance below 30 degrees of
#                                         depression (solved sky: 0.65 at 30 degrees; deeper is unsampled)
NATURAL_SKY_CD_M2 = 174e-6              # Falchi et al. 2016's natural zenith sky at solar minimum
AIRGLOW_SHARE = 0.48                    # airglow's share of Earth's natural sky (Duriscoe & Hung 2016)
FALCHI_EXTINCTION_MAG = 0.33             # Earth's vertical extinction at sea level in Falchi et al.'s atlas (K = 1)


def natural_floor(moon_tau):
    """Starlight, zodiacal and galactic light and airglow, which the light products omit: Earth's natural zenith
    (Falchi et al. 2016) at most; at least its extraterrestrial part, 52% (Duriscoe & Hung 2016), dimmed by the Moon's
    extra optical depth over Earth's (no airglow)."""
    earth_tau = FALCHI_EXTINCTION_MAG / (2.5 * math.log10(math.e))
    return ((1 - AIRGLOW_SHARE) * math.exp(-(moon_tau - earth_tau)) * NATURAL_SKY_CD_M2, NATURAL_SKY_CD_M2)
POLLUTED_SHARE = 14.0 / 174.0           # Falchi et al. 2016: light-polluted above 14 ucd/m2 of artificial brightness
# Earth's comparison air carries the haze of Falchi et al. 2016's atlas: aerosol clarity K = 1, a vertical extinction
# of 0.33 magnitudes at sea level, less the molecular part (derived); its scale height and phase are assumptions.
EARTH_AEROSOL_SHAPE = dict(scale_height_m=1.5e3, asymmetry=0.7, single_scattering_albedo=0.95)


def moon_haze(haze, case):
    """The design Moon's night haze from the haze product, averaged over its regions by area (one case)."""
    rows = [(r['area_share'], r['cases'][case]['night']) for r in haze['regions'].values()]
    w = sum(a for a, _ in rows)
    mean = lambda k: sum(a * n[k] for a, n in rows) / w
    return dict(aod=mean('aod'), scale_height_m=mean('mixed_layer_km') * 1e3, asymmetry=mean('asymmetry'),
                single_scattering_albedo=mean('single_scattering_albedo'))


def earth_aerosol(earth_col):
    """Falchi et al.'s atlas haze: its vertical extinction in magnitudes as optical depth, less the molecules'."""
    tau_total = FALCHI_EXTINCTION_MAG / (2.5 * math.log10(math.e))
    return dict(EARTH_AEROSOL_SHAPE, aod=tau_total - column_tau(scattering_layers('earth', earth_col), 120e3))


def skyglow_curves(transfer, winds, haze):
    moon_col, earth_col = columns(transfer, winds)
    aerosol = earth_aerosol(earth_col)
    cases = dict(
        rayleigh=(scattering_layers('moon', moon_col), scattering_layers('earth', earth_col)),
        hazy=(scattering_layers('moon', moon_col, moon_haze(haze, 'central')),
              scattering_layers('earth', earth_col, aerosol)))
    curves = {}
    for name, (lm, le) in cases.items():
        curves[name] = dict(moon=zenith_skyglow(SKY_DISTANCES_M, MOON_RADIUS, lm, top_m=500e3),
                            earth=zenith_skyglow(SKY_DISTANCES_M, EARTH_RADIUS, le, top_m=120e3))
    model = dict(
        moon_tau_550=column_tau(cases['rayleigh'][0], 500e3), earth_tau_550=column_tau(cases['rayleigh'][1], 120e3),
        moon_haze_central=moon_haze(haze, 'central'), earth_aerosol=aerosol,
        earth_rayleigh_slope_10_50_km=local_slope(SKY_DISTANCES_M, curves['rayleigh']['earth'], 10e3, 50e3),
        moon_rayleigh_slope_10_50_km=local_slope(SKY_DISTANCES_M, curves['rayleigh']['moon'], 10e3, 50e3),
        earth_hazy_slope_10_35_km=local_slope(SKY_DISTANCES_M, curves['hazy']['earth'], 10e3, 35e3),
        moon_hazy_slope_10_35_km=local_slope(SKY_DISTANCES_M, curves['hazy']['moon'], 10e3, 35e3),
        earth_hazy_slope_10_50_km=local_slope(SKY_DISTANCES_M, curves['hazy']['earth'], 10e3, 50e3),
        moon_hazy_slope_10_50_km=local_slope(SKY_DISTANCES_M, curves['hazy']['moon'], 10e3, 50e3),
        brightness_ratio_moon_over_earth={f'{d / 1e3:g} km': {k: float(np.exp(np.interp(math.log(d), np.log(SKY_DISTANCES_M),
                                          np.log(c['moon'] / c['earth'])))) for k, c in curves.items()}
                                          for d in (1e3, 5e3, 20e3, 100e3)})
    return curves, model


class Towns:
    """Disk profiles of towns in both worlds' air, calibrated on Earth to a polluted area of lit_m2 a person."""

    def __init__(self, curves, lit_m2):
        self.rows = {}
        for case, c in curves.items():
            for people in TOWN_PEOPLE + (METROPOLIS_PEOPLE,):
                settled = SETTLED_M2['metropolis' if people == METROPOLIS_PEOPLE else 'town']
                a = math.sqrt(settled * people / math.pi)
                pe = disk_profile(SKY_DISTANCES_M, c['earth'], a, XS_M)
                pm = disk_profile(SKY_DISTANCES_M, c['moon'], a, XS_M)
                r_target = math.sqrt(lit_m2 * people / math.pi)
                threshold = float(np.exp(np.interp(math.log(r_target), np.log(XS_M[1:]),
                                                   np.log(np.maximum(pe[1:], 1e-300)))))
                self.rows[(case, people)] = dict(earth=pe, moon=pm, threshold=threshold,
                                                 r_earth=polluted_radius(XS_M, pe, threshold), town_radius=a,
                                                 centre_over_threshold_moon=float(pm[0] / threshold))

    def factor(self, case, people, nu, q):
        """The Moon's polluted area over the Earth's for the same town, natural-sky ratio nu and light share q."""
        row = self.rows[(case, people)]
        if q <= 0:
            return 0.0
        return (polluted_radius(XS_M, row['moon'], nu * row['threshold'] / q) / row['r_earth']) ** 2


def flux_share(darkest_lux, outdoor_lux, design, residual, living_share=0.0, uplight=1.0):
    """The share of the shielded design's upward light still emitted at the night's darkest moment.

    shielded: every lamp all night, as people.py's lunar design assumes; twilight-aware: outdoor lamps add only what
    the natural light lacks of the outdoor level, while the residual (windows, vehicles, task lights) stays; living
    light then takes over a share of the outdoor service with a fraction of the lamps' upward flux."""
    if design == 'shielded':
        return 1.0
    lamps = max(0.0, 1.0 - darkest_lux / outdoor_lux)
    outdoor = lamps * (1.0 - living_share * (1.0 - uplight)) if design == 'living' else lamps
    return residual + (1.0 - residual) * outdoor


def point_factors(point, towns, ratio, outdoor_lux, floors):
    """Polluted-area factors at one reading of the grid: low and high over the zenith share, the town sizes, the air
    cases and the design guesses, for each design."""
    parts, e = point['darkest_parts'], point['earth_elevation_deg']
    nus = [natural_zenith_cd_m2(parts, e, ratio, k, f) / NATURAL_SKY_CD_M2 for k in TWILIGHT_K for f in floors]
    r_lo, r_hi = DESIGN_GUESSES['residual_share']
    p_lo, p_hi = DESIGN_GUESSES['living_light_share']
    u_lo, u_hi = DESIGN_GUESSES['living_light_uplight']
    qs = dict(shielded=[1.0],
              twilight=[flux_share(point['darkest_lux'], outdoor_lux, 'twilight', r) for r in (r_lo, r_hi)],
              living=[flux_share(point['darkest_lux'], outdoor_lux, 'living', r, p, u)
                      for r in (r_lo, r_hi) for p in (p_lo, p_hi) for u in (u_lo, u_hi)])
    out = dict(natural_zenith_over_earth=[min(nus), max(nus)])
    for design, q_list in qs.items():
        by_town = {}
        for people in TOWN_PEOPLE:
            vals = [towns.factor(case, people, nu, q) for case in ('rayleigh', 'hazy') for nu in nus for q in q_list]
            by_town[people] = [min(vals), max(vals)]
        out[design] = dict(flux_share=[min(q_list), max(q_list)],
                           area_factor=[min(v[0] for v in by_town.values()), max(v[1] for v in by_town.values())],
                           area_factor_by_town=by_town)
    return out


def home_weights(ways, home, grid):
    """Weights of the grid's readings for one home: its land by band, its Earth classes and the band's slices."""
    phys = ways['settings'][home]['physical']
    bands = phys['bands_dry_land_km2']
    total = sum(bands.values())
    shares = phys['earth_up_shares']
    weights = []
    for i, g in enumerate(grid):
        wb = bands[g['band']] / total
        wc = shares[g['earth_class']]
        if wb > 0 and wc > 0:
            weights.append((i, wb * wc / BAND_SLICES * g['class_weight']))
    s = sum(w for _, w in weights)
    return [(i, w / s) for i, w in weights]


def lit_shares(per_person, metropolis_m2, ways, dry_land_km2, mix, billions, aloft_share, aloft_m2):
    """Share of the dry land under a light-polluted sky for a number of lunar residents (billions), by a mix.

    Surface residents fill the homes by the mix; each home's lit land is capped at its own land. Residents aloft light
    aloft_m2 each (their towns spread over the whole Moon, a share of it over dry land). Returns [low, high]."""
    shares, metro_billion = mix_shares(ways, mix)
    land = {h: ways['settings'][h]['home_dry_land_km2'] for h in shares}
    out = []
    for end in (0, 1):
        n_surface = billions * (1 - aloft_share) - metro_billion
        lit = metro_billion * 1e9 * metropolis_m2[end] / 1e6
        for h, s in shares.items():
            lit += min(land[h], n_surface * s * 1e9 * per_person[h][end] / 1e6)
        lit += billions * aloft_share * 1e9 * aloft_m2[end] / 1e6
        out.append(min(1.0, lit / dry_land_km2))
    return out


def polluted_through_the_night(sky, towns, lat, people, k, floor, case='rayleigh', step_h=2.0):
    """A far-side town's polluted area through one night as the twilight fades and returns, against its darkest-hour
    area: hours with any polluted sky, and the night's mean area over the darkest hour's (twilight only, no Earth)."""
    t = np.arange(0.0, NIGHT_H + 1e-9, step_h)
    hour_angle = 90.0 + t * DEG_PER_H
    e = np.degrees(np.arcsin(np.cos(math.radians(lat)) * np.cos(np.radians(hour_angle))))
    nu = (k * sky.sun_diffuse(e) / math.pi + floor) / NATURAL_SKY_CD_M2
    f = np.array([towns.factor(case, people, x, 1.0) for x in nu])
    peak = f.max()
    return dict(hours_polluted=float(step_h * np.count_nonzero(f > 0)),
                mean_over_darkest=float(f.mean() / peak) if peak > 0 else 0.0,
                hours_above_half_of_darkest=float(step_h * np.count_nonzero(f > 0.5 * peak)) if peak > 0 else 0.0)


def forest_floor_darkest(point, canopy):
    """The night's darkest light on a forest floor (leaf area index 5): the daily-life study's floor shares, by the
    source's elevation for the Earth's light and its twilight share for the Sun's."""
    shares = canopy['floor_share_by_sun_deg']
    e = np.array(sorted(float(k) for k in shares))
    s = np.array([shares[f'{int(x)}' if float(x).is_integer() else str(x)] for x in e])
    parts = point['darkest_parts']
    earth_share = float(np.interp(max(point['earth_elevation_deg'], 1.0), e, s))
    return parts['earth_lux'] * earth_share + parts['sun_diffuse_lux'] * canopy['twilight_floor_share']


# ==== The product ====================================================================================================
def skyglow_section(nat, curves, model, ways, people, solved, sky, outdoor_lux):
    """Per-person polluted land by home and design, the metropolis, residents aloft and the shares by mix."""
    ratio = zenith_ratio_table(solved)
    lit = dict(lunar_design=people['needs']['practices']['lunar_design']['lit_sky_m2'],
               earth_practice=people['needs']['practices']['earth_practice']['lit_sky_m2'])
    dry_km2 = people['land']['dry_land_km2']
    grid = nat['grid']
    floors = natural_floor(model['moon_tau_550'])
    out = dict(model=model, earth_lit_sky_m2=lit, designs={}, towns={})
    for practice, (lo_m2, hi_m2) in lit.items():
        towns = {0: Towns(curves, lo_m2), 1: Towns(curves, hi_m2)}
        # Every town size is calibrated to the same Earth area a person, which sets its upward light a person: the
        # town's light is the threshold brightness over its profile's value per unit light at the calibrated radius.
        out['towns'][practice] = {f'{case} {people_n:.0e}': dict(
            town_radius_km=towns[0].rows[(case, people_n)]['town_radius'] / 1e3,
            earth_polluted_radius_km=[towns[i].rows[(case, people_n)]['r_earth'] / 1e3 for i in (0, 1)],
            centre_over_earth_threshold_moon=[towns[i].rows[(case, people_n)]['centre_over_threshold_moon'] for i in (0, 1)],
            upward_light_a_person_over_a_town_of_1e4=[
                towns[i].rows[(case, TOWN_PEOPLE[0])]['threshold'] * TOWN_PEOPLE[0]
                / (towns[i].rows[(case, people_n)]['threshold'] * people_n) for i in (0, 1)])
            for case in ('rayleigh', 'hazy') for people_n in TOWN_PEOPLE + (METROPOLIS_PEOPLE,)}
        factors = {i: [point_factors(g, towns[i], ratio, outdoor_lux, floors) for g in grid] for i in (0, 1)}
        designs = {}
        for design in ('shielded', 'twilight', 'living'):
            homes = {}
            for home in HOMES:
                if home == 'summit_metropolis':
                    continue
                w = home_weights(ways, home, grid)
                homes[home] = [sum(wi * factors[0][i][design]['area_factor'][0] for i, wi in w) * lo_m2,
                               sum(wi * factors[1][i][design]['area_factor'][1] for i, wi in w) * hi_m2]
            # the metropolis on its own: its 100 million at its own density, read at the summit
            summit = nat['places']['summit']
            nus = [natural_zenith_cd_m2(summit['darkest_parts'], summit['earth_elevation_deg'], ratio, k, f)
                   / NATURAL_SKY_CD_M2 for k in TWILIGHT_K for f in floors]
            qs = point_factors(summit, towns[0], ratio, outdoor_lux, floors)[design]['flux_share']
            metro, metro_town_light = [], []
            for i, m2 in ((0, lo_m2), (1, hi_m2)):
                vals = [towns[i].factor(case, METROPOLIS_PEOPLE, nu, q) for case in ('rayleigh', 'hazy') for nu in nus
                        for q in qs]
                metro.append((min(vals) if i == 0 else max(vals)) * m2)
                # the same 100 million, each lighting the sky as a resident of a calibrated town of 10,000 does
                rows = towns[i].rows
                vals = [towns[i].factor(case, METROPOLIS_PEOPLE, nu, q * rows[(case, METROPOLIS_PEOPLE)]['threshold']
                                        * METROPOLIS_PEOPLE / (rows[(case, TOWN_PEOPLE[0])]['threshold']
                                                               * TOWN_PEOPLE[0]))
                        for case in ('rayleigh', 'hazy') for nu in nus for q in qs]
                metro_town_light.append((min(vals) if i == 0 else max(vals)) * m2)
            # residents aloft: their towns over the whole Moon by area, the dry land's share of what they light
            moon_w = []
            for i, g in enumerate(grid):
                band = nat['classes_by_band'][g['band']]
                moon_w.append((i, band['area_share'] * band['classes'][g['earth_class']] / BAND_SLICES
                               * g['class_weight']))
            s = sum(w for _, w in moon_w)
            land_share = people['land']['dry_land']
            aloft = [land_share * sum(w / s * factors[0][i][design]['area_factor'][0] for i, w in moon_w) * lo_m2,
                     land_share * sum(w / s * factors[1][i][design]['area_factor'][1] for i, w in moon_w) * hi_m2]
            shares = {}
            for mix in list(ways['assumptions']['mixes']) + list(DERIVED_MIXES):
                shares[mix] = {f'{b:g}': dict(
                    all_on_the_land=lit_shares(homes, metro, ways, dry_km2, mix, b, 0.0, aloft),
                    half_aloft=lit_shares(homes, metro, ways, dry_km2, mix, b, ALOFT_SHARE, aloft),
                    all_on_the_land_metropolis_lit_like_a_town=lit_shares(homes, metro_town_light, ways, dry_km2, mix,
                                                                          b, 0.0, aloft))
                    for b in LUNAR_PEOPLE_BILLION}
                for b_key, row in shares[mix].items():
                    row['m2_per_person_all_on_the_land'] = [x * dry_km2 * 1e6 / (float(b_key) * 1e9)
                                                            for x in row['all_on_the_land']]
            by_town = {}
            for people_n in TOWN_PEOPLE:
                rows = {}
                for home in HOMES:
                    if home == 'summit_metropolis':
                        continue
                    w = home_weights(ways, home, grid)
                    rows[home] = [sum(wi * factors[0][i][design]['area_factor_by_town'][people_n][0] for i, wi in w)
                                  * lo_m2,
                                  sum(wi * factors[1][i][design]['area_factor_by_town'][people_n][1] for i, wi in w)
                                  * hi_m2]
                spread_share = {f'{b:g}': lit_shares(rows, metro, ways, dry_km2, 'spread', b, 0.0, aloft)
                                for b in LUNAR_PEOPLE_BILLION}
                by_town[f'{people_n:.0e}'] = dict(per_person_m2_by_home=rows, spread_shares_all_on_the_land=spread_share)
            land_total, land_dark = 0.0, 0.0
            for home in HOMES:
                if home == 'summit_metropolis':
                    continue
                land = ways['settings'][home]['home_dry_land_km2']
                land_total += land
                land_dark += land * sum(wi for i, wi in home_weights(ways, home, grid)
                                        if factors[1][i][design]['area_factor'][1] > 0)
            designs[design] = dict(per_person_m2_by_home=homes, metropolis_m2_per_person=metro,
                                   metropolis_dry_land_share=[x * METROPOLIS_PEOPLE / 1e6 / dry_km2 for x in metro],
                                   metropolis_m2_per_person_lit_like_a_town_of_1e4=metro_town_light,
                                   metropolis_dry_land_share_lit_like_a_town_of_1e4=[
                                       x * METROPOLIS_PEOPLE / 1e6 / dry_km2 for x in metro_town_light],
                                   aloft_m2_per_person=aloft, shares_of_dry_land=shares, by_town_size=by_town,
                                   dry_land_share_where_towns_light_their_sky=land_dark / land_total)
        out['designs'][practice] = designs
        out[f'grid_factors_{practice}'] = [dict(band=g['band'], lat=g['lat'], earth_class=g['earth_class'],
                                               darkest_lux=g['darkest_lux'],
                                               natural_zenith_over_earth=factors[0][i]['natural_zenith_over_earth'],
                                               area_factor={d: [factors[0][i][d]['area_factor'][0],
                                                                factors[1][i][d]['area_factor'][1]]
                                                            for d in ('shielded', 'twilight', 'living')})
                                          for i, g in enumerate(grid)]
    calibrated = [Towns(curves, m2) for m2 in lit['lunar_design']]
    night = {}
    for lat in (0.0, 5.4125, 15.0, 30.0):
        night[f'{lat:g}'] = {}
        for n in (1e5, METROPOLIS_PEOPLE):
            runs = [polluted_through_the_night(sky, t, lat, n, k, f, case) for t in calibrated
                    for case in ('rayleigh', 'hazy') for k in TWILIGHT_K for f in floors]
            night[f'{lat:g}'][f'{n:.0e}'] = {key: [min(r[key] for r in runs), max(r[key] for r in runs)]
                                             for key in runs[0]}
    out['through_the_night'] = night
    # The ground light skyglow adds at the polluted boundary, 8% of the natural zenith over pi/k, where towns pollute;
    # apart, the dark belt's readings (any hour below 0.01 lux).
    edge, edge_belt = [], []
    for i, g in enumerate(grid):
        f = point_factors(g, calibrated[1], ratio, outdoor_lux, floors)
        if f['shielded']['area_factor'][1] > 0:
            for nu in f['natural_zenith_over_earth']:
                lux = [POLLUTED_SHARE * nu * NATURAL_SKY_CD_M2 * math.pi / k for k in TWILIGHT_K]
                edge += lux
                if g['darkest_lux'] < MELATONIN_LUX[0]:
                    edge_belt += lux
    out['ground_lux_at_the_polluted_boundary'] = [min(edge), max(edge)]
    out['ground_lux_at_the_polluted_boundary_in_the_dark_belt'] = [min(edge_belt), max(edge_belt)]
    # The natural night sky at the darkest moment against Earth's, by band and Earth class, and the readings' latitudes
    # where towns of up to a million still light their sky above Falchi et al.'s level (shielded lamps).
    by_class = {}
    for g in out['grid_factors_lunar_design']:
        row = by_class.setdefault(g['band'], {}).setdefault(g['earth_class'], dict(
            natural_zenith_over_earth=[math.inf, 0.0], latitudes_where_towns_light_their_sky=[]))
        lo, hi = g['natural_zenith_over_earth']
        row['natural_zenith_over_earth'] = [min(row['natural_zenith_over_earth'][0], lo),
                                            max(row['natural_zenith_over_earth'][1], hi)]
        if g['area_factor']['shielded'][1] > 0 and g['lat'] not in row['latitudes_where_towns_light_their_sky']:
            row['latitudes_where_towns_light_their_sky'].append(g['lat'])
    out['natural_sky_by_band_and_class'] = by_class
    return out


def night_stores(plant):
    """The canopy model's night store and growth per m2 of ground, equatorial and 30-degree bands, both sides."""
    out = {}
    for side in ('near', 'far'):
        for band in ('0', '30'):
            for strategy in ('earth_like', 'idle_quarter'):
                r = plant['results'][side][band][strategy]
                out[f'{side} {band} {strategy}'] = dict(store_g_c_m2=r['store'], growth_g_c_m2_cycle=r['growth'],
                                                        respiration_in_the_dark_g_c_m2=r['respiration_in_the_dark'],
                                                        dark_hours=r['dark_hours'])
    return out


def emitter_rows(emitters, stores, needs_lux):
    """Each luminous organism's or engineered plant's light: luminance, what it gives a path, its cost."""
    store_lo = min(v['store_g_c_m2'] for v in stores.values())
    store_hi = max(v['store_g_c_m2'] for v in stores.values())
    growth_lo = min(v['growth_g_c_m2_cycle'] for v in stores.values())
    rows = []
    for e in emitters:
        q = e['photons_m2_s']
        lum = lambert_luminance(q, e['peak_nm'], e['fwhm_nm'])
        p = per_photon(e['peak_nm'], e['fwhm_nm'])
        row = dict(key=e['key'], organism=e['organism'], kind=e['kind'], measured=e['measured'], photons_m2_s=q,
                   peak_nm=e['peak_nm'],
                   fwhm_nm=e['fwhm_nm'], luminance_cd_m2=lum, scotopic_over_photopic=p['lm_s_scotopic'] / p['lm_s_photopic'],
                   optical_w_m2=q * p['photon_j'], source=e['source'], derivation=e.get('derivation', ''),
                   path_lux={g: math.pi * lum['photopic'] * v['view_factor'] for g, v in GEOMETRIES.items()},
                   brightness_needed_for={k: lux / (math.pi * lum['photopic']) for k, lux in needs_lux.items()})
        if e.get('chemistry') in PATHWAYS and e.get('quantum_yield'):
            costs = [glow_cost(q, qy, e['chemistry'], NIGHT_H, i, d) for qy in e['quantum_yield'] for i in (0, 1)
                     for d in DESIGN_GUESSES['delivered_share']]
            g = [c['g_c_m2'] for c in costs]
            by_share = {f'{d:g}': [glow_cost(q, qy, e['chemistry'], NIGHT_H, i, d)['g_c_m2']
                                   for qy in e['quantum_yield'] for i in (0, 1)]
                        for d in DESIGN_GUESSES['delivered_share']}
            row['night_cost'] = dict(g_c_m2_354_h=[min(g), max(g)],
                                     g_c_m2_354_h_by_delivered_share={k: [min(v), max(v)] for k, v in by_share.items()},
                                     share_of_night_store=[min(g) / store_hi, max(g) / store_lo],
                                     share_of_cycle_growth=[min(g) / max(v['growth_g_c_m2_cycle'] for v in stores.values()),
                                                            max(g) / growth_lo],
                                     atp_per_photon=[min(c['atp_per_photon'] for c in costs),
                                                     max(c['atp_per_photon'] for c in costs)])
        rows.append(row)
    return rows


def path_targets(needs_lux, stores, chemistries, hours_by_need):
    """The glow a closed glowing canopy needs to give each need level on the path below it, and its cost through the
    whole 354-hour night and through the hours the summit's night stays below that level."""
    store_lo = min(v['store_g_c_m2'] for v in stores.values())
    store_hi = max(v['store_g_c_m2'] for v in stores.values())
    out = {}
    for need, lux in needs_lux.items():
        row = {}
        for name, c in chemistries.items():
            q = photons_for_lux(lux, c['peak_nm'], c['fwhm_nm'])
            entry = dict(photons_m2_s=q, luminance_cd_m2=lux / math.pi)
            for label, hours in (('night', NIGHT_H), ('summit_hours', hours_by_need[need])):
                by_share = {}
                if c['pathway'] == 'floor':
                    mol_c = floor_mol_c_per_mol_photons(per_photon(c['peak_nm'], c['fwhm_nm'])['photon_j'])
                    g = [q * hours * 3600.0 / AVOGADRO * mol_c * C_G_PER_MOL] * 2
                else:
                    by_share = {f'{d:g}': [glow_cost(q, qy, c['pathway'], hours, i, d)['g_c_m2']
                                           for qy in c['quantum_yield'] for i in (0, 1)]
                                for d in DESIGN_GUESSES['delivered_share']}
                    g = [x for v in by_share.values() for x in v]
                entry[label] = dict(hours=hours, g_c_m2=[min(g), max(g)],
                                    share_of_night_store=[min(g) / store_hi, max(g) / store_lo])
                if by_share:
                    entry[label]['by_delivered_share'] = {
                        k: dict(g_c_m2=[min(v), max(v)], share_of_night_store=[min(v) / store_hi, max(v) / store_lo])
                        for k, v in by_share.items()}
            entry['night_mean_w_m2'] = [x * GLUCOSE_KJ_PER_G_C * 1e3 / NIGHT_S for x in entry['night']['g_c_m2']]
            row[name] = entry
        out[need] = row
    return out


def per_organism(needs_lux):
    """Organisms that are not glowing surfaces: how many dinoflagellate cells or bacteria under each square metre it
    takes to glow at a luminance, a firefly's flash against an engineered flower, and enzymes at full speed."""
    flashes = {}
    for name, f in FLASHES.items():
        lm = per_photon(f['peak_nm'], f['fwhm_nm'])['lm_s_photopic']
        flashes[name] = dict(peak_photons_s=f['peak_photons_s'], photons_per_flash=f['photons_per_flash'],
                             cells_m2_flashing_together_for={f'{lum:g} cd/m2': lum * math.pi / lm / f['peak_photons_s']
                                                            for lum in (0.01, 0.1, 1.0)}, sources=f['sources'])
    sea = next(e for e in EMITTERS if e['key'] == 'milky sea')
    flower = EMITTERS[0]['photons_m2_s'] / 1e6                # photons a second per mm2
    walk = photons_for_lux(needs_lux['walking'], 525.0, 80.0)
    per_umol = {k: v * AVOGADRO * 1e-6 for k, v in ENZYME_PHOTONS_S.items()}
    return dict(
        dinoflagellates=flashes,
        bacteria_m2_for_the_milky_sea=[sea['photons_m2_s'] / x for x in BACTERIA_PHOTONS_S[::-1]],
        firefly_flash_photons_per_mm2=[FIREFLY_FLASH['photons'][0] / FIREFLY_FLASH['lantern_mm2'][1],
                                       FIREFLY_FLASH['photons'][1] / FIREFLY_FLASH['lantern_mm2'][0]],
        firefly_flash_in_seconds_of_the_2020_flower=[FIREFLY_FLASH['photons'][0] / FIREFLY_FLASH['lantern_mm2'][1] / flower,
                                                     FIREFLY_FLASH['photons'][1] / FIREFLY_FLASH['lantern_mm2'][0] / flower],
        enzyme_photons_m2_s_per_umol_m2=per_umol,
        enzyme_over_walking_light_per_umol_m2={k: v / walk for k, v in per_umol.items()},
        firefly_luciferase_lux_per_umol_m2_closed_canopy=per_umol['firefly_luciferase']
        * per_photon(560.0, 70.0)['lm_s_photopic'])


def energy_limited_lux(stores, chemistries, store_shares=(0.01, 0.1)):
    """The path illuminance under a closed glowing canopy that a share of the night store buys through a 354-hour
    night, by chemistry: the energy limit of living light (the floor is the thermodynamic one)."""
    out = {}
    stores_lo_hi = (min(v['store_g_c_m2'] for v in stores.values()), max(v['store_g_c_m2'] for v in stores.values()))

    def path_lux(per, share, lm_s):
        lux = [share * store / C_G_PER_MOL / x * AVOGADRO / NIGHT_S * lm_s for store in stores_lo_hi for x in per]
        return [min(lux), max(lux)]

    for name, c in chemistries.items():
        p = per_photon(c['peak_nm'], c['fwhm_nm'])
        if c['pathway'] == 'floor':
            out[name] = {f'{share:g} of the store': path_lux([floor_mol_c_per_mol_photons(p['photon_j'])], share,
                                                             p['lm_s_photopic']) for share in store_shares}
            continue
        per = {f'{d:g}': [turnover_cost(c['pathway'], i) / qy / d / atp_values(i)['per_mol_C']
                          for qy in c['quantum_yield'] for i in (0, 1)] for d in DESIGN_GUESSES['delivered_share']}
        every = [x for v in per.values() for x in v]
        row = {f'{share:g} of the store': path_lux(every, share, p['lm_s_photopic']) for share in store_shares}
        row['by_delivered_share'] = {k: {f'{share:g} of the store': path_lux(v, share, p['lm_s_photopic'])
                                         for share in store_shares} for k, v in per.items()}
        out[name] = row
    return out


def marker_contrast(marker_cd_m2, background_lux, albedo=GREY):
    """Weber contrast of a glowing waymark against ground of this albedo under the natural light."""
    ground = albedo * background_lux / math.pi
    return (marker_cd_m2 - ground) / ground if ground > 0 else float('inf')


def marker_needed_cd_m2(background_lux, contrast, albedo=GREY):
    """A waymark's luminance for a given contrast against lit ground."""
    return albedo * background_lux / math.pi * (1 + contrast)


def hours_at(point, lux):
    """Hours per cycle the natural light at a reading stays below a level (interpolated in log level)."""
    levels = np.array(LEVELS_LUX)
    hours = np.array([point['hours_below'][f'{x:g}'] for x in LEVELS_LUX])
    if lux <= levels[0]:
        return float(hours[0])
    if lux >= levels[-1]:
        return float(hours[-1])
    return float(np.interp(math.log(lux), np.log(levels), hours))


def needs_by_setting(nat, ways, needs_lux, markers):
    """Each need against each setting's natural night: hours below its level, the lamp light it lacks at the darkest
    moment, and whether waymarks glowing at each of the markers' luminances stand out from the darkest ground
    (a Weber contrast of at least 1 against 18% grey)."""
    def verdicts(points):
        out = {}
        for need, lux in needs_lux.items():
            hours = [(w, hours_at(p, lux)) for w, p in points]
            lack = [(w, max(0.0, lux - p['darkest_lux'])) for w, p in points]
            out[need] = dict(level_lux=lux, hours_below_mean=sum(w * h for w, h in hours),
                             hours_below_range=[min(h for _, h in hours), max(h for _, h in hours)],
                             share_never_below=sum(w for w, h in hours if h == 0.0),
                             lamp_lux_at_darkest_range=[min(x for _, x in lack), max(x for _, x in lack)])
        return out

    places = {}
    for name, p in nat['places'].items():
        places[name] = dict(lat=p['lat'], lon=p['lon'], earth_elevation_deg=p['earth_elevation_deg'],
                            darkest_lux=p['darkest_lux'], midnight_lux=p['midnight_lux'],
                            waymark_contrast_at_darkest={k: marker_contrast(v, p['darkest_lux'])
                                                         for k, v in markers.items()},
                            needs=verdicts([(1.0, p)]))
    homes, land_never = {}, {need: 0.0 for need in needs_lux}
    land_total = 0.0
    for home in HOMES:
        if home == 'summit_metropolis':
            pts = [(1.0, nat['places']['summit'])]
        else:
            pts = [(w, nat['grid'][i]) for i, w in home_weights(ways, home, nat['grid'])]
        v = verdicts(pts)
        land = ways['settings'][home]['home_dry_land_km2']
        land_total += land
        for need in needs_lux:
            land_never[need] += land * v[need]['share_never_below']
        homes[home] = dict(dry_land_km2=land, darkest_lux_range=[min(p['darkest_lux'] for _, p in pts),
                                                                 max(p['darkest_lux'] for _, p in pts)],
                           waymark_visible_share={k: sum(w for w, p in pts if marker_contrast(v, p['darkest_lux'])
                                                         >= 1.0) for k, v in markers.items()},
                           needs=v)
    return dict(places=places, homes=homes,
                dry_land_share_where_the_night_never_falls_below={k: v / land_total for k, v in land_never.items()})


# Night-life thresholds already in the ecology register (biosphere/ecology/literature.md): melatonin in fishes and
# rodents suppressed at 0.01-0.03 lux (Grubisic et al. 2019); full moonlight on Earth 0.05-0.3 lux (Kyba et al. 2017).
MELATONIN_LUX = (0.01, 0.03)
FULL_MOON_LUX = (0.05, 0.3)
FLEET_ROLLED_LUX = 0.02                  # the ring fleet's night half rolled near 80 degrees with a Moon-facing scatter near
#                                          1e-5: about 0.02 lux at a far-side midnight (this branch's decisions register,
#                                          the fleet's night-zone entry of 2026-10-08, from the joint synthesis)


def night_life(nat, ways, canopy, calibrated, ratio):
    """Dark hours for night life by home, the forest floor's darkest light, and where a town's skyglow alone puts
    melatonin-level light on the ground."""
    homes = {}
    for home in HOMES:
        if home == 'summit_metropolis':
            pts = [(1.0, nat['places']['summit'])]
        else:
            pts = [(w, nat['grid'][i]) for i, w in home_weights(ways, home, nat['grid'])]
        homes[home] = dict(
            hours_below_melatonin_levels=[sum(w * hours_at(p, x) for w, p in pts) for x in MELATONIN_LUX],
            share_with_any_hour_below_0_01_lux=sum(w for w, p in pts if hours_at(p, MELATONIN_LUX[0]) > 0),
            forest_floor_darkest_lux=[min(forest_floor_darkest(p, canopy) for _, p in pts),
                                      max(forest_floor_darkest(p, canopy) for _, p in pts)])
    # Skyglow's own ground light, about pi L / k from a sky of zenith luminance L, against the melatonin levels: the
    # zenith luminance it takes, in Earth thresholds (14 ucd/m2), and how far from a town's centre the sky gets there.
    out_towns = {}
    need = [x * k / math.pi / (POLLUTED_SHARE * NATURAL_SKY_CD_M2) for x in MELATONIN_LUX for k in TWILIGHT_K]
    for people in TOWN_PEOPLE + (METROPOLIS_PEOPLE,):
        rows = [t.rows[(case, people)] for t in calibrated for case in ('rayleigh', 'hazy')]
        radii = [polluted_radius(XS_M, row['moon'], n * row['threshold']) / 1e3 for row in rows for n in need]
        out_towns[f'{people:.0e}'] = dict(
            centre_in_earth_thresholds=[min(r['centre_over_threshold_moon'] for r in rows),
                                        max(r['centre_over_threshold_moon'] for r in rows)],
            melatonin_level_in_earth_thresholds=[min(need), max(need)],
            radius_km_where_skyglow_alone_reaches_it=[min(radii), max(radii)],
            town_radius_km=rows[0]['town_radius'] / 1e3)
    belt = [g['darkest_lux'] for g in nat['grid'] if hours_at(g, MELATONIN_LUX[0]) > 0]
    belt.append(nat['places']['summit']['darkest_lux'])
    land = {home: ways['settings'][home]['home_dry_land_km2'] for home in HOMES}
    belt_km2 = sum(land[h] * homes[h]['share_with_any_hour_below_0_01_lux'] for h in HOMES)
    requirement = dict(dark_belt_km2=belt_km2, dark_belt_dry_land_share=belt_km2 / sum(land.values()),
                       dark_belt_darkest_lux=[min(belt), max(belt)],
                       falchi_share_of_the_dark_belt_darkest_lux=[POLLUTED_SHARE * min(belt), POLLUTED_SHARE * max(belt)],
                       fleet_rolled_lux=FLEET_ROLLED_LUX,
                       fleet_over_dark_belt_darkest=[FLEET_ROLLED_LUX / max(belt), FLEET_ROLLED_LUX / min(belt)])
    return dict(homes=homes, skyglow_ground_light=out_towns, melatonin_lux=list(MELATONIN_LUX),
                full_moon_lux=list(FULL_MOON_LUX), dark_night_requirement=requirement)


EVIDENCE = (
    'A screen on committed products and cited values. The natural night is the calendar\'s evaluator on its transfer '
    'product in the mean geometry (clear sky, level ground, no libration, refraction, terrain, cloud or haze), '
    'reproducing the provisioning branch\'s daily-life calendars and main\'s per-place Earthlight. Living light is '
    'measured photon outputs and quantum yields from the literature, turned into luminance and path light by CIE '
    'photometry with Gaussian stand-ins for the emission spectra, and costed from each luciferin cycle\'s stoichiometry '
    'at Hinkle\'s P/O ratios, without the cost of making the enzymes. Needs are Earth lighting standards and '
    'experiments. Skyglow is a single-scattering model of a uniform disk town in the Moon\'s and Earth\'s air, '
    'calibrated on ecology\'s Earth-statistics areas a person; multiple scattering, spectra, towns\' shapes and their '
    'overlaps are not modelled. No organism, lamp or town here exists or has been tested.')
READING_RULE = (
    'Light is lux on level ground and cd/m2 for surfaces (photopic unless named scotopic). Hours are per 708.73-hour '
    'cycle. [low, high] pairs bound each quantity over the ranges named in design_guesses and in the cited values. '
    'emitters.organisms gives each luminous organism\'s or engineered plant\'s light and its cost through a 354-hour '
    'night; emitters.targets the glow and the cost a closed glowing canopy needs for each need on the path below it; '
    'emitters.limits the path light a share of the night store buys. needs gives, per place and per ways-of-living '
    'home, the hours below each need\'s level and the lamp light lacking at the darkest moment. skyglow gives the area '
    'a person under a light-polluted sky at the night\'s darkest moment by home and design, in m2, and the share of '
    'the dry land by provisioning\'s mixes at 4, 6 and 8 billion lunar residents, all on the land or half aloft. '
    'Values computed here are screen values; arithmetic on cited values is derived; design guesses are ranges.')


def results():
    """The whole product as a dict (deterministic; the stored product equals it)."""
    data = {p: load_product(ROOT / p, s) for p, s in INPUTS.items()}
    bound, sources = {}, {}
    for path in CROSS_BRANCH:
        bound[path], sources[path] = cross_branch_json(path)
        if bound[path] is None:
            raise RuntimeError(f'A bound cross-branch input is unreadable: {path} ({sources[path]})')
    places = bound['research/studies/joint_synthesis/results/places.json']
    daily = bound['research/studies/daily_life/results/daily_life.json']
    ways = bound['research/studies/ways_of_living/results/ways_of_living.json']
    transfer = data['illumination/calendar/results/transfer.json']
    solved = data['illumination/sky/results/solved_sky.json']
    people = data['biosphere/ecology/results/people.json']
    sky = Sky(transfer)
    nat = natural_light(sky, daily, places)
    curves, model = skyglow_curves(transfer, data['climate/results/gcm/global_winds_A28_dim5_moon.json'],
                                   data['illumination/aerosol/results/haze.json'])
    stores = night_stores(data['biosphere/canopy/results/plant.json'])
    needs_lux = {k: v['lux'] for k, v in NEEDS.items()}
    sky_section = skyglow_section(nat, curves, model, ways, people, solved, sky, OUTDOOR_LUX)
    calibrated = [Towns(curves, m2) for m2 in people['needs']['practices']['lunar_design']['lit_sky_m2']]
    emitters = dict(
        photometry=dict(km_photopic_lm_w=KM_PHOTOPIC, km_scotopic_lm_w=KM_SCOTOPIC,
                        photons_per_s_m2_for_one_lux={name: photons_for_lux(1.0, c['peak_nm'], c['fwhm_nm'])
                                                      for name, c in CHEMISTRIES.items()}),
        pathways={k: dict(note=v['note'], atp_equivalents_per_turnover=[turnover_cost(k, i) for i in (0, 1)],
                          quantum_yield=next(c['quantum_yield'] for c in CHEMISTRIES.values() if c['pathway'] == k),
                          atp_equivalents_per_photon=[turnover_cost(k, i) / qy for i in (0, 1) for qy in
                                                      next(c['quantum_yield'] for c in CHEMISTRIES.values()
                                                           if c['pathway'] == k)])
                  for k, v in PATHWAYS.items()},
        chemistries={k: dict(v) for k, v in CHEMISTRIES.items()},
        cysteine_atp_equivalents=[cysteine_atp(dict(atp_values(i))) for i in (0, 1)],
        atp_per_mol_c=[atp_values(i)['per_mol_C'] for i in (0, 1)],
        floor_mol_c_per_mol_photons={name: floor_mol_c_per_mol_photons(per_photon(c['peak_nm'], c['fwhm_nm'])['photon_j'])
                                     for name, c in CHEMISTRIES.items()},
        night_stores=stores,
        organisms=emitter_rows(EMITTERS, stores, needs_lux),
        geometries={k: v['view_factor'] for k, v in GEOMETRIES.items()},
        targets=path_targets({k: v for k, v in needs_lux.items() if v <= 10.0}, stores, CHEMISTRIES,
                             {k: hours_at(nat['places']['summit'], v) for k, v in needs_lux.items()}),
        per_organism=per_organism(needs_lux),
        limits=energy_limited_lux(stores, dict(CHEMISTRIES, **{'best beetle luciferase': dict(
            CHEMISTRIES['firefly chemistry'], quantum_yield=(BEST_BEETLE_QY,), peak_nm=540.0)})),
        mesopic_at_walking_level={name: mesopic_luminance(GREY * needs_lux['walking'] / math.pi,
                                                          GREY * needs_lux['walking'] / math.pi
                                                          * per_photon(c['peak_nm'], c['fwhm_nm'])['lm_s_scotopic']
                                                          / per_photon(c['peak_nm'], c['fwhm_nm'])['lm_s_photopic'])
                                  / (GREY * needs_lux['walking'] / math.pi)
                                  for name, c in CHEMISTRIES.items() if c['pathway'] != 'floor'})
    plants = [e for e in EMITTERS if e['kind'] == 'engineered plant']
    lum = lambda e: lambert_luminance(e['photons_m2_s'], e['peak_nm'], e['fwhm_nm'])['photopic']
    waymarks = {'brightest measured plant': max(lum(e) for e in plants if e['measured']),
                'improved pathway at a hundredfold': max(lum(e) for e in plants if not e['measured'])}
    needs = needs_by_setting(nat, ways, needs_lux, waymarks)
    needs['waymarks_cd_m2'] = waymarks
    life = night_life(nat, ways, daily['volume']['under_canopy'], calibrated, zenith_ratio_table(solved))
    taxa = data['biosphere/ecology/results/taxa.json']
    registered = {e['id'] for key in ('groups', 'modules', 'communities') for e in taxa.get(key, [])}
    own = [HERE / 'run.py']
    product = dict(
        schema=SCHEMA,
        producer=dict(
            study='living_light',
            files={str(p.relative_to(ROOT)): digest(p) for p in own},
            imported_code={p: digest(ROOT / p) for p in IMPORTED_CODE},
            inputs={p: digest(ROOT / p) for p in list(INPUTS) + [SOURCES]},
            cross_branch_inputs={p: dict(branch=s['branch'], commit=s['commit'], blob=s['blob'], sha256=s['sha256'],
                                         read_from=sources[p],
                                         bind='Read in-tree once the branches are integrated; until then from Git '
                                              'at the commit.')
                                 for p, s in CROSS_BRANCH.items()},
            constants=constants_used([str(p) for p in own] + list(IMPORTED_CODE))),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        design_guesses=dict(DESIGN_GUESSES, outdoor_lux=OUTDOOR_LUX),
        natural_sky_assumptions=dict(twilight_zenith_share=TWILIGHT_K, floor_cd_m2=natural_floor(model['moon_tau_550']),
                                     earth_natural_cd_m2=NATURAL_SKY_CD_M2, polluted_share=POLLUTED_SHARE),
        needs_basis=NEEDS,
        natural=dict(places={k: {f: v[f] for f in ('lat', 'lon', 'earth_elevation_deg', 'midnight_lux', 'darkest_lux',
                                                    'darkest_hours_after_sunset', 'darkest_parts', 'hours_below')}
                             for k, v in nat['places'].items()},
                     main_places=nat['main_places'], checks_against_daily_life=nat['checks_against_daily_life'],
                     classes_by_band=nat['classes_by_band'],
                     grid=[{f: g[f] for f in ('band', 'lat', 'lon', 'earth_class', 'earth_elevation_deg', 'darkest_lux',
                                              'darkest_parts', 'hours_below')} for g in nat['grid']]),
        emitters=emitters, needs=needs, skyglow=sky_section, night_life=life,
        taxonomy=dict(identifiers_used=list(TAXA_USED), all_registered=all(t in registered for t in TAXA_USED)),
        no_claims=['No organism, glowing plant, lamp or town named here exists or has been tested.',
                   'Need levels are Earth standards and experiments with people adapted to Earth\'s nights.',
                   'Skyglow areas are Falchi et al.\'s astronomical measure carried to the Moon; ecological effects '
                   'begin at other levels.'])
    product['headlines'] = headlines(product)
    producer = plain(product.pop('producer'), None)
    return dict(schema=SCHEMA, producer=producer, **{k: v for k, v in plain(product).items() if k != 'schema'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'research/runs/living_light/living_light.json')
    args = parser.parse_args()
    product = results()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    print(args.output)



def headlines(product):
    """The numbers the README quotes, gathered in one place."""
    em, nd, sk = product['emitters'], product['needs'], product['skyglow']
    lunar = sk['designs']['lunar_design']
    shares = {design: {mix: {b: v['all_on_the_land'] for b, v in rows.items()}
                       for mix, rows in lunar[design]['shares_of_dry_land'].items()}
              for design in lunar}
    aloft = {design: {mix: {b: v['half_aloft'] for b, v in rows.items()}
                      for mix, rows in lunar[design]['shares_of_dry_land'].items()}
             for design in lunar}
    per_person = {design: dict(homes=lunar[design]['per_person_m2_by_home'],
                               metropolis=lunar[design]['metropolis_m2_per_person'],
                               aloft=lunar[design]['aloft_m2_per_person']) for design in lunar}
    brightest = max((r for r in em['organisms'] if r['kind'] == 'engineered plant' and r['measured']),
                    key=lambda r: r['photons_m2_s'])
    return dict(
        brightest_engineered_plant=dict(key=brightest['key'], photons_m2_s=brightest['photons_m2_s'],
                                        luminance_cd_m2=brightest['luminance_cd_m2'],
                                        closed_canopy_lux=brightest['path_lux']['closed canopy overhead']),
        targets=em['targets'], limits=em['limits'],
        dry_land_share_where_the_night_never_falls_below=nd['dry_land_share_where_the_night_never_falls_below'],
        skyglow_model=dict(moon_tau=sk['model']['moon_tau_550'], earth_tau=sk['model']['earth_tau_550'],
                           brightness_ratio=sk['model']['brightness_ratio_moon_over_earth']),
        per_person_m2=per_person, shares_all_on_the_land=shares, shares_half_aloft=aloft,
        through_the_night=sk['through_the_night'])


def mesopic_luminance(photopic_cd_m2, scotopic_cd_m2, iterations=40):
    """CIE 191:2010's MES2 mesopic luminance (as Goodman 2013 sets out its iteration): m from 0.5, then
    L_mes = [m L_p + (1 - m) L_s V'(555)] / [m + (1 - m) V'(555)] with V'(555) = 683/1700, m = 0.767 + 0.3334 log L_mes,
    m held to 0-1 (1 at 5 cd/m2 and above, 0 at 0.005 and below)."""
    v = KM_PHOTOPIC / KM_SCOTOPIC
    m = 0.5
    lmes = photopic_cd_m2
    for _ in range(iterations):
        lmes = (m * photopic_cd_m2 + (1 - m) * scotopic_cd_m2 * v) / (m + (1 - m) * v)
        m = min(1.0, max(0.0, 0.7670 + 0.3334 * math.log10(max(lmes, 1e-12))))
    return lmes


if __name__ == '__main__':
    main()
