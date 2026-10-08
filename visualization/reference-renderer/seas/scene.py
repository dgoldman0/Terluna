"""Scene inputs for the sea-scene path tracer at the sea-appearance study's four coasts.

Everything physical comes from the domains and the study: the moments and their sea states from the study's
products, the waves drawn from the hour's spectra (illumination/water_surface/realization.py), the sky and the
air in front of each surface from the solved spherical sky (illumination/sky/aerial.py), the Earth's light
(illumination/earthlight through the study's lighting), the water's own light (illumination/water_column with
the biosphere's productive-coast guess), the coast from LOLA heights above the geoid (geography/coast_terrain.py)
at the tide's level (the geography tide product), and the stars (illumination/stars). The renderer adds
only geometry: the camera, the explicit waves and the coast's shading.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from numba import njit, prange

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
import kernels as k  # noqa: E402
from geography import coast_terrain, tides  # noqa: E402
from illumination.sky import aerial  # noqa: E402
from illumination.sky.solved_transport import evaluate  # noqa: E402
from illumination.water_column import model as water_optics  # noqa: E402
from illumination.water_surface import realization  # noqa: E402
from illumination.water_surface.model import wavenumber  # noqa: E402
from illumination.water_surface.reflection import fresnel, refractive_index  # noqa: E402
from research.studies.sea_appearance import lighting, regimes, scenes  # noqa: E402
from research.studies.sea_appearance.scenes import CAMERAS, LAND_SOIL, TERRAIN_NAMES, VIEWS  # noqa: E402
from shared.constants import AU, EARTH_RADIUS, MOON_RADIUS, MOON_SURFACE_GRAVITY, SUN_RADIUS  # noqa: E402

STUDY = ROOT / "research/studies/sea_appearance"
REGIMES = STUDY / "results/regimes.json"
SLOPE_SERIES = STUDY / "results/sea_slopes.npz"
SHORE_MONTH = ROOT / "climate/waves/results/shore_month.json"
TIDES = ROOT / "geography/products/tides_28pct_4ppd.npz"
TERRAIN = ROOT / "research/runs/sea_appearance/terrain"
WATERS = ROOT / "biosphere/living_water/waters.json"
STARS = ROOT / "illumination/stars/bright_stars.json"
EYE_HEIGHT_M = 2.0
WATER = "productive_coast"
# Spectral bands for the coast's colour carried through the air (nm edges); narrower where the air's
# transmittance changes fastest.
BAND_EDGES = np.array([360.0, 430.0, 470.0, 510.0, 550.0, 590.0, 630.0, 680.0, 831.0])
# Sky map: elevation rows dense at the horizon, absolute azimuth columns.
SKY_EL = np.unique(np.r_[np.arange(-0.12, 0.1, 0.005), np.arange(0.1, 1.0, 0.02), np.arange(1.0, 5.0, 0.05),
                         np.arange(5.0, 20.0, 0.2), np.arange(20.0, 90.01, 0.5)])
SKY_REL_AZ = np.unique(np.r_[np.arange(0.0, 5.0, 0.1), np.arange(5.0, 20.0, 0.25), np.arange(20.0, 180.01, 1.0)])
SKY_AZ_STEP = 0.25
# The air in front of surfaces.
LUT_EL = np.array([-90.0, -40.0, -20.0, -10.0, -6.0, -4.0, -3.0, -2.0, -1.5, -1.0, -0.7, -0.5, -0.35, -0.25, -0.18,
                   -0.12, -0.08, -0.04, 0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0,
                   11.0, 15.0, 25.0, 90.0])
LUT_REL_AZ = np.array([0.0, 5.0, 10.0, 20.0, 30.0, 45.0, 60.0, 80.0, 100.0, 120.0, 150.0, 180.0])
LUT_D = np.array([0.0, 25.0, 50.0, 100.0, 200.0, 400.0, 700.0, 1000.0, 1500.0, 2000.0, 2600.0, 3500.0, 5000.0, 7000.0,
                  1e4, 1.4e4, 2e4, 2.8e4, 4e4, 5.6e4, 8e4, 1.1e5, 1.6e5, 3e5])
LUT_AZ_STEP = 5.0
SEA_BLOCK = 16
TERRAIN_BLOCK = 8
TEXTURE_STRIDE = 2           # the coast's light at 200 m, its shape at 100 m
HORIZON_AZIMUTHS = 72

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit(elevation_deg, azimuth_deg):
    e, a = math.radians(elevation_deg), math.radians(azimuth_deg)
    return np.array([math.cos(e) * math.sin(a), math.cos(e) * math.cos(a), math.sin(e)])


class Common:
    """What every scene shares: the solved sky, the study's light helpers and products."""

    def __init__(self):
        with np.load(regimes.SOLUTION, allow_pickle=False) as z:
            self.solution = {key: z[key] for key in z.files if key != "meta"}
            self.meta = json.loads(str(z["meta"]))
        self.xyz = np.ascontiguousarray(self.solution["xyz"], dtype=np.float64)
        self.wavelength = self.solution["wavelength"]
        self.band = np.clip(np.searchsorted(BAND_EDGES, self.wavelength) - 1, 0, len(BAND_EDGES) - 2)
        self.sky = lighting.Sky()
        self.sky.channel_wavelength = self.wavelength
        self.earth_light = lighting.Earthlight(self.sky)
        self.regimes = json.loads(REGIMES.read_text())
        self.calendar = json.loads(regimes.CALENDAR.read_text())
        self.month = json.loads(SHORE_MONTH.read_text())
        self.hours = np.array(self.calendar["time_hours"], float)
        self.jd0 = self.month["tide_month"]["middle_jd_tt"]
        with np.load(SLOPE_SERIES) as z:
            self.series = {key: z[key] for key in z.files}
        self.waters = json.loads(WATERS.read_text())
        self.n_water = float(refractive_index(550.0))
        self.inputs = {str(p.relative_to(ROOT)): sha256(p) for p in
                       (REGIMES, SLOPE_SERIES, SHORE_MONTH, TIDES, WATERS, STARS, regimes.CALENDAR)}
        self.inputs[str(regimes.SOLUTION.relative_to(ROOT))] = sha256(regimes.SOLUTION)
        with np.load(TIDES, allow_pickle=False) as z:
            self.tides = {key: z[key] for key in ("lat_deg", "lon_deg", "rows", "cols", "basis")}
        days = np.arange(0.0, tides.YEARS * 365.25, tides.STEP_DAYS)
        self.tide_mean = tides.forcing_coefficients(tides.JD_START + days).mean(axis=0)

    def jd(self, hour):
        return self.jd0 + (hour - 1.5 * lighting.CYCLE_HOURS) / 24

    def band_xyz(self, channels):
        """XYZ per spectral band (bands, 3) of per-channel radiance or irradiance."""
        out = np.zeros((len(BAND_EDGES) - 1, 3))
        for b in range(out.shape[0]):
            sel = self.band == b
            out[b] = channels[sel] @ self.xyz[sel]
        return out

    def tide(self, lon, lat, jd):
        """Water level (m) relative to its 2026-2045 mean at the tide grid's nearest sea cell."""
        t = self.tides
        lat_c, lon_c = t["lat_deg"][t["rows"]], t["lon_deg"][t["cols"]]
        u = lambda lo, la: np.stack((np.cos(np.radians(la)) * np.cos(np.radians(lo)),
                                     np.cos(np.radians(la)) * np.sin(np.radians(lo)), np.sin(np.radians(la))), -1)
        pixel = int(np.argmax(u(lon_c, lat_c) @ u(lon, lat)))
        return float((tides.forcing_coefficients(np.array([jd])) - self.tide_mean) @ t["basis"][:, pixel])


# --------------------------------------------------------------------------------------------- the sea state
def sea_state(common, coast, hour, seed):
    """The hour's explicit waves: tiles, filtered-slope pyramids and the residual covariance."""
    key = coast.replace(" ", "_").lower()
    i = int(np.argmin(np.abs(common.hours - hour)))
    u_star = float(common.series[f"{key}/u_star"][i])
    toward = float(common.series[f"{key}/toward_deg"][i])
    spec = scenes.sea_spectrum(coast, hour, common.calendar, common.month, common.hours, toward)
    g = MOON_SURFACE_GRAVITY
    f, e, depth = spec["frequency"], spec["variance"], spec["depth"]
    fp = float(f[np.argmax(e.sum(axis=1))])
    kp = float(wavenumber(np.array([fp]), depth, g)[0])
    peak_speed = 2 * math.pi * fp / kp
    k_from = float(wavenumber(np.array([f[-1]]), depth, g)[0])
    densities = [realization.swan_density(f, spec["direction"], e, depth, g),
                 realization.short_density(u_star, toward, g, k_from, peak_speed)]
    tiles = realization.realize(densities, depth, g, seed)
    residual = realization.residual_covariance(u_star, toward, g, tiles[-1].k_high, peak_speed)
    # Filtered-slope pyramids, flattened.
    pyramids = [realization.lean_pyramid(t.slope) for t in tiles]
    levels = max(len(p) for p in pyramids)
    offsets = np.zeros((len(tiles), levels), np.int64)
    sizes = np.ones((len(tiles), levels), np.int64)
    chunks, position = [], 0
    for ti, pyramid in enumerate(pyramids):
        for li, level in enumerate(pyramid):
            offsets[ti, li] = position
            sizes[ti, li] = level.shape[1]
            chunks.append(level.ravel())
            position += level.size
    lean = np.concatenate(chunks).astype(np.float32)
    heights = tiles[0].height
    n = heights.shape[0]
    blocks = heights.reshape(n // SEA_BLOCK, SEA_BLOCK, n // SEA_BLOCK, SEA_BLOCK).max(axis=(1, 3))
    # A block's cells reach the next blocks' first nodes.
    first = heights[::SEA_BLOCK, ::SEA_BLOCK]
    blocks = np.maximum.reduce([blocks, np.roll(first, -1, 0), np.roll(first, -1, 1),
                                np.roll(np.roll(first, -1, 0), -1, 1),
                                np.roll(heights[::SEA_BLOCK, :].reshape(n // SEA_BLOCK, n // SEA_BLOCK, SEA_BLOCK)
                                        .max(axis=2), -1, 0),
                                np.roll(heights[:, ::SEA_BLOCK].reshape(n // SEA_BLOCK, SEA_BLOCK, n // SEA_BLOCK)
                                        .max(axis=1), -1, 1)]).astype(np.float32)
    record = dict(
        spectrum=spec["source"], spectrum_file=spec["file"], depth_m=round(depth, 2),
        friction_velocity_m_s=round(u_star, 4), wind_toward_deg=round(toward, 1),
        wind_10m_m_s=round(float(common.series[f"{key}/u10"][i]), 2), peak_period_s=round(1 / fp, 2),
        significant_height_m=round(float(4 * np.std(heights)), 3), seed=seed,
        tiles=[dict(size_m=t.size_m, cells=t.cells, wavelengths_m=[round(2 * math.pi / t.k_high, 4),
               None if t.k_low == 0 else round(2 * math.pi / t.k_low, 3)], choppy=t.choppy,
               slope_covariance_expected=np.round(t.expected_covariance, 6).tolist(),
               slope_covariance_drawn=np.round(realization.sample_covariance(t), 6).tolist()) for t in tiles],
        residual_covariance=np.round(residual, 7).tolist(),
        mean_square_slope=round(float(sum(np.trace(realization.sample_covariance(t)) for t in tiles)
                                      + np.trace(residual)), 5))
    return dict(heights=heights, blocks=blocks, hmin=float(heights.min()), hmax=float(heights.max()),
                dx=tiles[0].spacing_m, lean=lean, offsets=offsets, sizes=sizes,
                spacing=np.array([t.spacing_m for t in tiles]),
                levels=np.array([len(p) for p in pyramids], np.int64),
                residual=np.array([residual[0, 0], residual[1, 1], residual[0, 1]]), record=record)


# ----------------------------------------------------------------------------------------------- the light
def sky_field(common, source_el, weights, height=EYE_HEIGHT_M):
    inscatter, _, _ = aerial.view(common.solution, source_el, SKY_EL, SKY_REL_AZ, [1e9], height_m=height,
                                  weights=weights)
    return inscatter[:, :, 0, :]


def absolute(field, rel_grid, source_az, azimuths):
    """A field over azimuth from its source, laid onto absolute azimuths."""
    rel = np.abs((azimuths - source_az + 180) % 360 - 180)
    j = np.clip(np.searchsorted(rel_grid, rel) - 1, 0, len(rel_grid) - 2)
    v = ((rel - rel_grid[j]) / (rel_grid[j + 1] - rel_grid[j]))
    v = v.reshape((1, -1) + (1,) * (field.ndim - 2))
    return (1 - v) * field[:, j] + v * field[:, j + 1]


def light(common, geometry, weights):
    """Sky map, the air's LUTs, the disks and the water's own light for one moment."""
    sun_el, sun_az = geometry["sun_el"], geometry["sun_az"]
    earth_el, earth_az = geometry["earth_el"], geometry["earth_az"]
    earth_on = earth_el > -30
    azimuths = np.arange(0.0, 360.0, SKY_AZ_STEP)
    sky = absolute(sky_field(common, sun_el, None), SKY_REL_AZ, sun_az, azimuths)
    if earth_on:
        sky = sky + absolute(sky_field(common, earth_el, weights), SKY_REL_AZ, earth_az, azimuths)
    floor = lighting.FLOOR_LUX / np.pi * np.array([regimes.D65[0] / regimes.D65[1], 1.0,
                                                    (1 - sum(regimes.D65)) / regimes.D65[1]])
    sky = sky + floor
    # The air in front of surfaces, from the eye.
    lut_az = np.arange(0.0, 360.0, LUT_AZ_STEP)
    s_sun, t_ch, _ = aerial.view(common.solution, sun_el, LUT_EL, LUT_REL_AZ, LUT_D, height_m=EYE_HEIGHT_M)
    lut_s = absolute(s_sun, LUT_REL_AZ, sun_az, lut_az)
    if earth_on:
        s_earth, _, _ = aerial.view(common.solution, earth_el, LUT_EL, LUT_REL_AZ, LUT_D, height_m=EYE_HEIGHT_M,
                                    weights=weights)
        lut_s = lut_s + absolute(s_earth, LUT_REL_AZ, earth_az, lut_az)
    # Transmittance per band, weighted by daylight's spectrum within each band, and per XYZ for the sea.
    daylight = common.sky.ground_channels(30.0)[0]
    lut_tb = np.zeros((len(LUT_EL), len(LUT_D), len(BAND_EDGES) - 1))
    lut_t = np.zeros((len(LUT_EL), len(LUT_D), 3))
    for b in range(lut_tb.shape[2]):
        sel = common.band == b
        w = daylight[sel] * common.xyz[sel].sum(axis=1)
        lut_tb[:, :, b] = t_ch[:, :, sel] @ w / w.sum()
    for c in range(3):
        w = daylight * common.xyz[:, c]
        lut_t[:, :, c] = t_ch @ w / w.sum()
    # The disks through the air.
    sun_radius = math.asin(SUN_RADIUS / (AU * geometry["sunlight_factor"] ** -0.5))
    sun_normal = regimes.normal_beam(common.sky, sun_el) @ common.xyz
    sun_omega = 2 * math.pi * (1 - math.cos(sun_radius))
    sun = np.zeros(10)
    sun[0:3], sun[3] = unit(sun_el, sun_az), sun_radius
    sun[4:7], sun[7:10] = sun_normal, sun_normal / sun_omega
    earth_radius = math.asin(EARTH_RADIUS / geometry["earth_distance_m"])
    earth_omega = 2 * math.pi * (1 - math.cos(earth_radius))
    earth_normal = (regimes.normal_beam(common.sky, earth_el) * weights) @ common.xyz
    lit = max(geometry["earth_lit_fraction"], 1e-4)
    earth = np.zeros(12)
    earth[0:3], earth[3] = unit(earth_el, earth_az), earth_radius
    earth[4:7] = unit(sun_el, sun_az)
    earth[7:10] = earth_normal / (lit * earth_omega)
    earth[10], earth[11] = earth_omega, 1.0 if earth_el > -earth_radius * 180 / math.pi else 0.0
    # The water's own light (as the regime study computes it).
    downwelling = common.sky.ground_channels(sun_el)[0] + weights * common.sky.ground_channels(earth_el)[0]
    w = common.waters["waters"][WATER]
    iop = water_optics.water(common.wavelength, w["chlorophyll_mg_m3"]["best"],
                             w["dissolved_organic_440_per_m"]["best"], w["fines_g_m3"]["best"], w["soil"])
    rrs = water_optics.remote_sensing_reflectance(iop["a"], iop["bb"])
    water = (rrs * downwelling) @ common.xyz
    # Coarse sky per band for lighting the coast: radiance per channel by elevation and absolute azimuth.
    coarse_el = np.r_[np.arange(0.5, 10.0, 1.0), np.arange(11.0, 90.0, 3.0)]
    coarse_rel = np.arange(0.0, 180.1, 5.0)
    sky_ch = evaluate(common.solution, common.meta, np.array([sun_el]), coarse_el, coarse_rel)[0]
    coarse_az = np.arange(2.5, 360.0, 5.0)
    sky_b = absolute(sky_ch, coarse_rel, sun_az, coarse_az)                     # (el, az, channels)
    if earth_on:
        sky_e = evaluate(common.solution, common.meta, np.array([earth_el]), coarse_el, coarse_rel)[0] * weights
        sky_b = sky_b + absolute(sky_e, coarse_rel, earth_az, coarse_az)
    sky_bands = np.stack([[common.band_xyz(sky_b[i, j]) for j in range(len(coarse_az))]
                          for i in range(len(coarse_el))])                      # (el, az, bands, 3)
    sky_bands += (floor / len(BAND_EDGES[:-1]))[None, None, None, :]
    sun_bands = common.band_xyz(regimes.normal_beam(common.sky, sun_el))
    earth_bands = common.band_xyz(regimes.normal_beam(common.sky, earth_el) * weights)
    return dict(sky=sky.astype(np.float64), azimuths=azimuths, lut_s=lut_s, lut_t=lut_t, lut_tb=lut_tb, sun=sun,
                earth=earth, water=water, coarse_el=coarse_el, coarse_az=coarse_az, sky_bands=sky_bands,
                sun_bands=sun_bands, earth_bands=earth_bands, floor=floor)


# ------------------------------------------------------------------------------------------------- the coast
@njit(parallel=True, cache=True)
def coast_light(h, dx, stride, inv2r, offset, sky_bands, coarse_el, coarse_az, sun, sun_bands, sun_radius, earth,
                earth_bands, earth_radius, albedo, n_az):
    """Lambertian radiance (bands, XYZ) of the coast on a grid every stride nodes: direct light where the
    terrain's horizon lets the disks through, the sky above that horizon, and a single bounce off the slopes
    below it lit as this cell is."""
    ny, nx = h.shape
    oy, ox = (ny - 1) // stride + 1, (nx - 1) // stride + 1
    nb = sky_bands.shape[2]
    out = np.zeros((oy, ox, nb, 3), np.float32)
    steps = np.empty(64)
    d = dx
    for s in range(64):
        steps[s] = d
        d *= 1.09
    for jj in prange(oy):
        j = jj * stride
        for ii in range(ox):
            i = ii * stride
            z0 = h[j, i] + offset
            if z0 < -2.0:
                continue
            # Normal from central differences.
            gx = (h[j, min(i + 1, nx - 1)] - h[j, max(i - 1, 0)]) / (dx * (min(i + 1, nx - 1) - max(i - 1, 0)))
            gy = (h[min(j + 1, ny - 1), i] - h[max(j - 1, 0), i]) / (dx * (min(j + 1, ny - 1) - max(j - 1, 0)))
            nn = math.sqrt(gx * gx + gy * gy + 1)
            n0, n1, n2 = -gx / nn, -gy / nn, 1 / nn
            # The horizon's elevation in each azimuth sector.
            horizon = np.full(n_az, -0.5)
            for a in range(n_az):
                az = math.radians((a + 0.5) * 360.0 / n_az)
                sx, sy = math.sin(az), math.cos(az)
                best = -0.5
                for s in range(64):
                    dist = steps[s]
                    x = i + sx * dist / dx
                    y = j + sy * dist / dx
                    if x < 0 or y < 0 or x > nx - 1 or y > ny - 1:
                        break
                    xi, yi = int(x), int(y)
                    fx, fy = x - xi, y - yi
                    x1, y1 = min(xi + 1, nx - 1), min(yi + 1, ny - 1)
                    zz = ((1 - fx) * (1 - fy) * h[yi, xi] + fx * (1 - fy) * h[yi, x1] + (1 - fx) * fy * h[y1, xi]
                          + fx * fy * h[y1, x1]) + offset - dist * dist * inv2r
                    el = math.degrees(math.atan2(zz - z0, dist))
                    if el > best:
                        best = el
                horizon[a] = best
            irr = np.zeros((nb, 3))
            hidden = 0.0
            for ei in range(len(coarse_el)):
                el = coarse_el[ei]
                lo = el - (coarse_el[ei] - coarse_el[ei - 1]) / 2 if ei > 0 else 0.0
                hi = el + (coarse_el[ei + 1] - coarse_el[ei]) / 2 if ei + 1 < len(coarse_el) else 90.0
                band = (math.sin(math.radians(hi)) - math.sin(math.radians(lo))) * math.radians(360.0 / len(coarse_az))
                ce, se = math.cos(math.radians(el)), math.sin(math.radians(el))
                for ai in range(len(coarse_az)):
                    az = math.radians(coarse_az[ai])
                    cosine = n0 * ce * math.sin(az) + n1 * ce * math.cos(az) + n2 * se
                    if cosine <= 0:
                        continue
                    sector = int(coarse_az[ai] / (360.0 / n_az)) % n_az
                    if el > horizon[sector]:
                        for b in range(nb):
                            for c in range(3):
                                irr[b, c] += sky_bands[ei, ai, b, c] * cosine * band
                    else:
                        hidden += cosine * band
            for src, bands, radius in ((sun, sun_bands, sun_radius), (earth, earth_bands, earth_radius)):
                el = math.degrees(math.asin(src[2]))
                az = math.degrees(math.atan2(src[0], src[1])) % 360.0
                sector = int(az / (360.0 / n_az)) % n_az
                seen = min(max((el - horizon[sector]) / (2 * math.degrees(radius)) + 0.5, 0.0), 1.0)
                cosine = n0 * src[0] + n1 * src[1] + n2 * src[2]
                if seen > 0 and cosine > 0:
                    for b in range(nb):
                        for c in range(3):
                            irr[b, c] += bands[b, c] * cosine * seen
            for b in range(nb):
                bounce = 1.0 + albedo[b] * hidden / math.pi
                for c in range(3):
                    out[jj, ii, b, c] = albedo[b] * irr[b, c] * bounce / math.pi
    return out


def coast(common, coast_name, camera_xy, tide_m, light_, sun_radius, earth_radius):
    """The coast's heights round the camera at the tide's level, its trace blocks and its light."""
    with np.load(TERRAIN / f"{TERRAIN_NAMES[coast_name]}.npz", allow_pickle=False) as z:
        h = z["height_m"].astype(np.float64)
        x = z["x_m"].astype(float)
        meta = json.loads(str(z["metadata"]))
    dx = float(x[1] - x[0])
    x0, y0 = x[0] - camera_xy[0], x[0] - camera_xy[1]
    offset = -tide_m                         # heights above the water at this moment
    n = h.shape[0]
    cells = n - 1
    nb = cells // TERRAIN_BLOCK
    blocks = np.full((nb, nb), -1e9)
    for bj in range(nb):
        for bi in range(nb):
            blocks[bj, bi] = h[bj * TERRAIN_BLOCK:(bj + 1) * TERRAIN_BLOCK + 1,
                               bi * TERRAIN_BLOCK:(bi + 1) * TERRAIN_BLOCK + 1].max()
    lam, refl = water_optics.soil_reflectance(LAND_SOIL[coast_name])
    albedo_ch = np.interp(common.wavelength, lam, refl)
    albedo = np.array([np.average(albedo_ch[common.band == b], weights=common.xyz[common.band == b, 1] + 1e-12)
                       for b in range(len(BAND_EDGES) - 1)])
    # The coast's own shading uses positions relative to the camera for the Moon's curvature.
    hc = h.copy()
    yy, xx = np.meshgrid(x - camera_xy[1], x - camera_xy[0], indexing="ij")
    texture = coast_light(hc, dx, TEXTURE_STRIDE, 0.5 / MOON_RADIUS, offset, light_["sky_bands"], light_["coarse_el"],
                          light_["coarse_az"], light_["sun"][0:3], light_["sun_bands"], sun_radius,
                          light_["earth"][0:3], light_["earth_bands"], earth_radius, albedo, HORIZON_AZIMUTHS)
    record = dict(source=meta["evidence"], terrain_sha256=sha256(TERRAIN / f"{TERRAIN_NAMES[coast_name]}.npz"),
                  tide_m=round(tide_m, 3), land_soil=water_optics.SOILS[LAND_SOIL[coast_name]]["name"],
                  land_albedo_by_band=np.round(albedo, 4).tolist())
    return dict(h=h.astype(np.float32), blocks=blocks.astype(np.float32), x0=x0, y0=y0, dx=dx, offset=offset,
                hmax=float(h.max() + offset), texture=texture.astype(np.float64), tex_dx=dx * TEXTURE_STRIDE,
                record=record)


# --------------------------------------------------------------------------------------------------- a scene
def build(common, coast_name, moment, width, height, fov_deg=65.0, seed=20380207):
    s = next(x for x in common.regimes["scenes"] if x["coast"] == coast_name and x["moment"] == moment)
    hour = s["hour"]
    site = next(x for x in common.calendar["sites"] if x["short"] == coast_name)
    camera = CAMERAS[coast_name]
    lat, lon = coast_terrain.destination(*site["lon_lat"], np.array(camera["offset_m"][0]),
                                         np.array(camera["offset_m"][1]))
    lon = float((lon + 180) % 360 - 180)
    lat = float(lat)
    jd = common.jd(hour)
    g = lighting.geometry(np.array([jd]), lon, lat)
    geometry = dict(sun_el=float(g["sun_elevation"][0]), sun_az=float(g["sun_azimuth"][0]),
                    earth_el=float(g["earth_elevation"][0]), earth_az=float(g["earth_azimuth"][0]),
                    earth_phase_deg=float(g["earth_phase_angle"][0]),
                    earth_lit_fraction=float(g["earth_lit_fraction"][0]),
                    earth_distance_m=float(g["earth_distance_m"][0]),
                    sunlight_factor=float(g["earth_sunlight_factor"][0]))
    weights = common.earth_light.weights(g["earth_phase_angle"], g["earth_distance_m"],
                                         g["earth_sunlight_factor"])[0]
    light_ = light(common, geometry, weights)
    sea = sea_state(common, coast_name, hour, seed)
    tide = common.tide(lon, lat, jd)
    coast_ = coast(common, coast_name, camera["offset_m"], tide, light_, light_["sun"][3], light_["earth"][3])
    azimuth, pitch = VIEWS[(coast_name, moment)]
    f = unit(pitch, azimuth)
    right = np.array([math.cos(math.radians(azimuth)), -math.sin(math.radians(azimuth)), 0.0])
    up = np.cross(right, f)
    eye_z = float(k.surface_height(sea["heights"], 0.0, 0.0, 0.0, 0.0, sea["dx"], True)) + EYE_HEIGHT_M
    camera_vec = np.r_[[0.0, 0.0, eye_z], f, right, up]
    params = np.zeros(24)
    params[k.P_R] = MOON_RADIUS
    params[k.P_SEA_DX], params[k.P_SEA_HMIN], params[k.P_SEA_HMAX] = sea["dx"], sea["hmin"], sea["hmax"]
    params[k.P_N] = common.n_water
    params[k.P_F0] = float(fresnel(1.0, common.n_water))
    params[k.P_TER_X0], params[k.P_TER_Y0], params[k.P_TER_DX] = coast_["x0"], coast_["y0"], coast_["dx"]
    params[k.P_TER_OFFSET], params[k.P_TER_HMAX], params[k.P_HAS_TERRAIN] = coast_["offset"], coast_["hmax"], 1.0
    params[k.P_TEX_X0], params[k.P_TEX_Y0], params[k.P_TEX_DX] = coast_["x0"], coast_["y0"], coast_["tex_dx"]
    params[k.P_SKY_AZ_STEP], params[k.P_LUT_AZ_STEP] = SKY_AZ_STEP, LUT_AZ_STEP
    params[k.P_PIXEL] = 2 * math.tan(math.radians(fov_deg) / 2) / width     # a pixel's angle at the centre
    params[k.P_INTERREFLECT] = 4 * sea["dx"]
    params[k.P_SEA_BLOCK], params[k.P_TER_BLOCK] = SEA_BLOCK, TERRAIN_BLOCK
    params[k.P_TAN_HALF], params[k.P_ASPECT] = math.tan(math.radians(fov_deg) / 2), width / height
    params[k.P_SUN_TEST] = math.sin(math.radians(15.0))
    arrays = dict(params=params, sea_h=sea["heights"], sea_blocks=sea["blocks"], lean=sea["lean"],
                  lean_offsets=sea["offsets"], lean_sizes=sea["sizes"], lean_spacing=sea["spacing"],
                  lean_levels=sea["levels"], residual=sea["residual"], ter_h=coast_["h"],
                  ter_blocks=coast_["blocks"], tex=coast_["texture"], sky=light_["sky"],
                  sky_el=SKY_EL.astype(np.float64), lut_s=light_["lut_s"], lut_t=light_["lut_t"],
                  lut_tb=light_["lut_tb"], lut_el=LUT_EL, lut_d=LUT_D, sun=light_["sun"], earth=light_["earth"],
                  water=light_["water"])
    record = dict(coast=coast_name, moment=moment, hour=hour, date_tt=jd, camera_lon_lat=[round(lon, 5), round(lat, 5)],
                  camera_note=camera["note"], eye_height_m=EYE_HEIGHT_M, view_azimuth_deg=azimuth,
                  view_pitch_deg=pitch, horizontal_fov_deg=fov_deg, width=width, height=height,
                  sun=dict(elevation_deg=round(geometry["sun_el"], 3), azimuth_deg=round(geometry["sun_az"], 2)),
                  earth=dict(elevation_deg=round(geometry["earth_el"], 3), azimuth_deg=round(geometry["earth_az"], 2),
                             phase_angle_deg=round(geometry["earth_phase_deg"], 2),
                             lit_fraction=round(geometry["earth_lit_fraction"], 3)),
                  sea=sea["record"], land=coast_["record"], water=WATER,
                  water_leaving_nadir_cd_m2=round(float(light_["water"][1]), 5))
    return camera_vec, arrays, record, geometry, light_
