"""Data for the globe's appearance mode: surface colour, relief, clouds, ground light and haze.

Called by render.py. Every input is a domain product, and each piece says what kind it is:

- computed: the atmosphere (Rayleigh coefficients and scale height from the illumination domain's
  engine-atmosphere product), the ground's direct and diffuse light against Sun elevation (the
  Open Moon clear-sky atlas, 55 scattering orders), and land, sea and relief at 16 pixels per degree
  (LOLA above the GRAIL geoid through the geography domain, at the atlas product's sea level), with
  the geography domain's rivers and rain-fed lakes above sea level;
- informed: cloud cover from the climate domain's ExoPlaSim climatology (run A), a low and a high deck
  each moved with the Sun through its hour-angle composite, placed at the deck's cover-weighted
  height and carried east by the deck's cover-weighted mean wind; light scattered more than once, an
  isotropic source scaled from the sky atlas's diffuse light so that the air seen straight down
  matches Eddington's conservative-scattering reflectance;
- guesstimate: surface cover from rainfall, soil moisture, nearness to water and rivers, height, hollows
  and slope (forest, woodland, grassland, bare highland or basaltic soil, rock), river width from discharge,
  water colour from depth, and the
  cloud decks' shapes (noise the page evaluates, thresholded to the run's cover) and opacity. The
  biosphere and climate work replace it.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / 'illumination' / 'sky' / 'products' / 'engine_atmosphere.json'
SKY_ATLAS = ROOT / 'illumination' / 'sky' / 'data' / 'moon_atlas.npz'
CLIMATOLOGY = ROOT / 'climate' / 'gcm' / 'products' / 'climatology_A.npz'
DRAINAGE = ROOT / 'geography' / 'products' / 'drainage_28pct_16ppd.npz'
DRAINAGE_SCHEMA = 'terluna.geography.drainage/1'
# Channel width from discharge: Earth's relation, about 7.2 m per (m3/s)^0.5, widened by 1.6 because flow in
# lunar gravity is slower (speed with the square root of g) and needs a larger cross-section.
RIVER_WIDTH = 7.2 * 1.6
LUMINANCE = np.array([0.2126, 0.7152, 0.0722])
LUT_ELEVATIONS = np.linspace(-90.0, 90.0, 121)   # twilight light reaches far past the terminator
MOON_KM = 1737.4
NORMAL_STRENGTH = 10.0                           # relief is baked into the normal map at this exaggeration

# Guesstimated reflectances (linear sRGB), after typical Earth surfaces.
FOREST = np.array([0.030, 0.058, 0.026])
SAVANNA = np.array([0.070, 0.088, 0.042])
GRASS = np.array([0.135, 0.132, 0.075])
HIGHLAND_SOIL = np.array([0.34, 0.31, 0.25])
BASALT_SOIL = np.array([0.16, 0.10, 0.07])
ROCK = np.array([0.20, 0.19, 0.18])
DEEP_WATER = np.array([0.008, 0.018, 0.036])
MID_WATER = np.array([0.012, 0.027, 0.046])
TURBID_WATER = np.array([0.034, 0.048, 0.042])
RIVER_WATER = np.array([0.085, 0.078, 0.058])   # sediment-laden, as rivers crossing fresh ground run
CLOUD_ALBEDO = 0.72                 # of the part of a deck that intercepts light
LOW_OPACITY, HIGH_OPACITY = 0.9, 0.4  # guesstimates: the low deck holds the model's cloud water, the high deck little
CLOUD_ROWS = 90                     # 2-degree rows for the cloud maps
CALIBRATION_FLOOR_DEG = 10.0        # below this Sun height the scattering calibration holds its 10-degree value

# Cloud noise, evaluated by the page and replicated here to measure its distribution. The low deck is
# warped fractal noise with features from about 400 km down to a few km; the high deck is drawn out
# four times longer east-west than north-south.
CLOUD_NOISE = dict(lowFrequency=4.0, lowOctaves=7, lowWarp=0.35, highFrequency=3.0, highOctaves=6,
                   highWarp=0.5, highStretch=4.0, gain=0.55, lacunarity=2.03)
QUANTILES = 33


def smoothstep(a, b, x):
    t = np.clip((x - a) / (b - a), 0.0, 1.0).astype(np.float32)
    return t * t * (3 - 2 * t)


def fbm(shape, seed, octaves=7, base=(6, 3), gain=0.55, stretch=1.0):
    """Fractal value noise in float32, seamless in longitude: each octave is a random lattice enlarged with
    cubic interpolation. stretch > 1 elongates features east-west."""
    import cv2
    rng = np.random.default_rng(seed)
    h, w = shape
    out = np.zeros(shape, np.float32)
    amp, total = 1.0, 0.0
    for o in range(octaves):
        nx, ny = int(base[0] * 2 ** o), max(2, int(round(base[1] * 2 ** o * stretch)))
        if nx > w // 2 or ny > h // 2:
            break
        lattice = rng.random((ny, nx), dtype=np.float32)
        pad = np.concatenate([lattice[:, -2:], lattice, lattice[:, :2]], axis=1)
        big = cv2.resize(pad, (int(round(pad.shape[1] * w / nx)), h), interpolation=cv2.INTER_CUBIC)
        x0 = int(round(2 * w / nx))
        out += amp * big[:, x0:x0 + w]
        total += amp
        amp *= gain
    return out / total


def regrid(field, lat_src, lon_src, lat_dst, lon_dst):
    """Bilinear interpolation from the model grid (north-first latitudes, east 0-360) to the texture grid."""
    order = np.argsort(lat_src)
    ls, f = lat_src[order], field[order]
    lon_ext = np.concatenate([lon_src - 360.0, lon_src, lon_src + 360.0])
    f_ext = np.concatenate([f, f, f], axis=1)
    la = np.clip(lat_dst, ls[0], ls[-1])
    i = np.clip(np.searchsorted(ls, la) - 1, 0, ls.size - 2)
    ty = ((la - ls[i]) / (ls[i + 1] - ls[i]))[:, None]
    j = np.clip(np.searchsorted(lon_ext, lon_dst % 360.0) - 1, 0, lon_ext.size - 2)
    tx = ((lon_dst % 360.0 - lon_ext[j]) / (lon_ext[j + 1] - lon_ext[j]))[None, :]
    a = f_ext[i][:, j] * (1 - tx) + f_ext[i][:, j + 1] * tx
    b = f_ext[i + 1][:, j] * (1 - tx) + f_ext[i + 1][:, j + 1] * tx
    return a * (1 - ty) + b * ty


def wrapped(filter_, a, pad=32):
    """Apply an image filter as if longitude wrapped around: pad east and west from the far side, then crop."""
    out = filter_(np.concatenate([a[:, -pad:], a, a[:, :pad]], axis=1))
    return out[:, pad:pad + a.shape[1]]


def gradients(a, dy):
    """North-south and east-west differences per metre of arc (east-west before the cos(latitude) factor),
    central and wrapping in longitude."""
    gy = np.gradient(a, axis=0) / dy
    gx = (np.roll(a, -1, axis=1) - np.roll(a, 1, axis=1)) / (2.0 * dy)
    return gy, gx


def climate_field(clim, name, shape):
    """A climatology field on the texture grid: regridded at 4 pixels per degree, then enlarged."""
    import cv2
    lat4 = 90.0 - (np.arange(720) + 0.5) / 4.0
    lon4 = (np.arange(1440) + 0.5) / 4.0 - 180.0
    f = regrid(clim[name], clim['lat'], clim['lon'], lat4, lon4).astype(np.float32)
    s = shape[1] // f.shape[1]
    pad = np.concatenate([f[:, -1:], f, f[:, :1]], axis=1)
    return cv2.resize(pad, (pad.shape[1] * s, shape[0]), interpolation=cv2.INTER_LINEAR)[:, s:s + shape[1]]


def drainage_grids(shape, path=DRAINAGE):
    """Discharge (m3/s) along channels and the level of rain-fed lakes above sea level on the texture grid
    (180 W at the left edge), from the geography domain's drainage product; None when it is absent."""
    if not Path(path).is_file():
        return None
    z = np.load(path)
    meta = json.loads(str(z['metadata']))
    if meta.get('schema') != DRAINAGE_SCHEMA or meta['grid']['rows'] != shape[0]:
        raise SystemExit(f'Expected {DRAINAGE_SCHEMA} on a {shape[0]}-row grid in {path}')
    q = np.zeros(shape[0] * shape[1], np.float32)
    q[z['channel_cell']] = z['channel_discharge_m3s']
    lake = np.full(shape[0] * shape[1], np.nan, np.float32)
    lake[z['lake_cell']] = z['lake_level_m']
    roll = shape[1] // 2
    return dict(discharge=np.roll(q.reshape(shape), roll, axis=1), lake_level=np.roll(lake.reshape(shape), roll, axis=1),
                metadata=meta)


def surface_colour(height, level, clim, ppd, drainage=None):
    """Guesstimated surface reflectance (linear sRGB, float32), the water mask (sea and lakes) and the share
    of each cell that river channels cover, 180 W at the left edge."""
    import cv2
    height = height.astype(np.float32)
    rows, cols = height.shape
    lat = 90.0 - (np.arange(rows) + 0.5) / ppd
    km = MOON_KM * np.radians(1.0 / ppd)
    sea = height < level
    if drainage is not None:
        lake = np.isfinite(drainage['lake_level'])
        discharge = np.where(sea | lake, 0.0, drainage['discharge'])
        river = np.clip(RIVER_WIDTH * np.sqrt(discharge) / (km * 1000.0), 0.0, 1.0).astype(np.float32)
        surface_level = np.where(lake, drainage['lake_level'], level).astype(np.float32)
    else:
        lake = np.zeros_like(sea)
        discharge = river = None
        surface_level = np.float32(level)
    water = sea | lake

    rain = climate_field(clim, 'pr_mm_day', height.shape)
    wet = rain
    wet += 1.5 * smoothstep(0.02, 0.2, climate_field(clim, 'mrso_m', height.shape))
    pad = cols // 8
    tiled = np.concatenate([water[:, -pad:], water, water[:, :pad]], axis=1)
    shore_km = cv2.distanceTransform((~tiled).astype(np.uint8), cv2.DIST_L2, 5)[:, pad:pad + cols] * km
    wet += 1.4 * np.exp(-shore_km / 120.0)
    del tiled, shore_km
    land_height = np.maximum(height, level)
    wet += 0.8 * np.exp(-(land_height - level) / 400.0)
    curvature = wrapped(lambda a: cv2.Laplacian(cv2.GaussianBlur(a, (0, 0), 3.0), cv2.CV_32F, ksize=3), land_height)
    wet += 0.6 * np.clip(curvature / (curvature[~water].std() + 1e-6), -2.5, 2.5)   # wetter in hollows
    del curvature
    wet += 2.5 * (fbm(height.shape, 11, octaves=9, base=(24, 12)) - 0.5)
    if discharge is not None:                        # floodplains: wetter ground along larger rivers
        wet += 0.9 * np.minimum(wrapped(lambda a: cv2.GaussianBlur(a, (0, 0), 2.5), np.log10(1.0 + discharge).astype(np.float32)), 2.0)
    veg = smoothstep(0.4, 3.2, wet)
    dense = smoothstep(3.2, 5.5, wet)
    del wet
    basaltic = smoothstep(level + 900.0, level + 150.0, height) * smoothstep(
        0.35, 0.65, fbm(height.shape, 17, octaves=8, base=(18, 9)))
    gy, gx = gradients(land_height, km * 1000.0)
    rocky = 0.5 * smoothstep(0.35, 0.7, np.hypot(gx / np.maximum(np.cos(np.radians(lat)), 0.05)[:, None], gy))
    del gy, gx, land_height
    grass = smoothstep(0.0, 0.35, veg)
    wooded = smoothstep(0.35, 1.0, veg)
    del veg
    depth = np.where(water, surface_level - height, 0.0)
    shallow = smoothstep(250.0, 30.0, depth)
    deep = smoothstep(400.0, 1500.0, depth)
    del depth
    colour = np.empty(height.shape + (3,), np.float32)
    for c in range(3):
        bare = HIGHLAND_SOIL[c] * (1 - basaltic) + BASALT_SOIL[c] * basaltic
        land = (bare * (1 - grass) + GRASS[c] * grass) * (1 - wooded) + (SAVANNA[c] * (1 - dense) + FOREST[c] * dense) * wooded
        land = land * (1 - rocky) + ROCK[c] * rocky
        sea = (MID_WATER[c] * (1 - deep) + DEEP_WATER[c] * deep) * (1 - shallow) + TURBID_WATER[c] * shallow
        if river is not None:
            land = land * (1 - river) + RIVER_WATER[c] * river
        colour[..., c] = np.where(water, sea, land)
    return colour, water, river, surface_level


def relief_normals(height, level, ppd):
    """Tangent-space normals of the surface with water flat at its level, relief baked at NORMAL_STRENGTH."""
    surf = np.maximum(height, level).astype(np.float32)
    lat = 90.0 - (np.arange(surf.shape[0]) + 0.5) / ppd
    dy = MOON_KM * 1000.0 * np.radians(1.0 / ppd)
    gy, gx = gradients(surf, dy)
    n = np.stack([-(gx / np.maximum(np.cos(np.radians(lat)), 1e-3)[:, None]) * NORMAL_STRENGTH,
                  gy * NORMAL_STRENGTH, np.ones_like(surf)], axis=-1)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    return np.round((n + 1) / 2 * 255).astype(np.uint8)


def block_mean(a, k):
    """Mean over k x k blocks."""
    r, c = a.shape[0] // k, a.shape[1] // k
    return a[:r * k, :c * k].reshape(r, k, c, k).mean(axis=(1, 3))


def srgb_encode(x):
    x = np.clip(x, 0.0, 1.0)
    return np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1 / 2.4) - 0.055)


def engine_atmosphere():
    d = json.loads(ENGINE.read_text())
    return d['atmospheres']['moon'], d


def ground_light_lut(sun_rgb_lux):
    """Direct (on a horizontal surface) and diffuse irradiance from the clear-sky atlas, relative to the
    top-of-atmosphere sunlight in each channel, at LUT_ELEVATIONS."""
    a = np.load(SKY_ATLAS)
    suns = a['suns']
    direct = np.array([np.interp(LUT_ELEVATIONS, suns, a['direct_horizontal'][:, c]) for c in range(3)]).T
    diffuse = np.array([np.interp(LUT_ELEVATIONS, suns, a['diffuse'][:, c]) for c in range(3)]).T
    return np.clip(direct, 0, None) / sun_rgb_lux, np.clip(diffuse, 0, None) / sun_rgb_lux


def eddington_reflectance(tau, mu0):
    """Reflectance of a conservatively scattering plane layer over a black surface (Eddington)."""
    return (tau + (2.0 / 3.0 - mu0) * (1.0 - np.exp(-tau / mu0))) / (4.0 / 3.0 + tau)


def single_scatter_nadir(tau, mu0):
    """Single-scattering Rayleigh reflectance of the same layer, seen straight down."""
    phase = 0.75 * (1.0 + mu0 * mu0)
    return phase / (4.0 * (mu0 + 1.0)) * (1.0 - np.exp(-tau * (1.0 / mu0 + 1.0)))


def diffuse_source_factor(moon, diffuse, elevations=LUT_ELEVATIONS):
    """Per-channel factor k for light scattered more than once, treated as an isotropic source of k times the
    sky atlas's diffuse light at the scattering point's Sun height. k is set so that the air seen straight
    down matches Eddington's conservative-scattering reflectance of the same column, and holds its
    10-degree value for lower Suns and beyond the terminator."""
    tau = np.array(moon['rayleigh_scattering_fitted_per_km']) * moon['rayleigh_scale_height_km']
    held = np.maximum(elevations, CALIBRATION_FLOOR_DEG)
    mu0 = np.sin(np.radians(held))[:, None]
    sky = np.array([np.interp(held, elevations, diffuse[:, c]) for c in range(3)]).T
    higher_orders = (eddington_reflectance(tau[None, :], mu0) - single_scatter_nadir(tau[None, :], mu0)) * mu0
    return higher_orders / (sky * (1.0 - np.exp(-tau[None, :])))


def rows_to_uniform(field, lat_src, rows=CLOUD_ROWS):
    """Interpolate a north-first model field (lat x any) to uniform rows, north first."""
    centres = 90.0 - (np.arange(rows) + 0.5) * 180.0 / rows
    order = np.argsort(lat_src)
    return np.stack([np.interp(centres, lat_src[order], field[order, k]) for k in range(field.shape[1])], axis=1)


def deck_heights_km(clim, h_km):
    """Cover-weighted height of the low and high decks, heights from sigma with the atmosphere's scale height."""
    return tuple(deck_mean(clim, -h_km * np.log(clim['sigma'])))


def deck_mean(clim, values):
    """Cover-weighted mean of a per-layer quantity over the low and the high deck."""
    cover = clim['cl_layer_mean']
    low = clim['sigma'] > clim['deck_sigma']
    return [float((values * cover)[m].sum() / cover[m].sum()) for m in (low, ~low)]


# The page's cloud noise, replicated with the same integer hash so that its distribution can be measured.
def _hash(x):
    x = x.astype(np.uint32)
    x ^= x >> np.uint32(16)
    x *= np.uint32(0x7FEB352D)
    x ^= x >> np.uint32(15)
    x *= np.uint32(0x846CA68B)
    x ^= x >> np.uint32(16)
    return x


def _lattice(c, seed):
    u = lambda v: v.astype(np.int64).astype(np.uint32)
    h = _hash(u(c[..., 2]) + np.uint32(seed))
    h = _hash(u(c[..., 1]) + h)
    return _hash(u(c[..., 0]) + h).astype(np.float64) / 4294967295.0


def value_noise(x, seed):
    i = np.floor(x)
    f = x - i
    u = f * f * f * (f * (f * 6.0 - 15.0) + 10.0)
    c = i.astype(np.int64)
    corner = lambda dx, dy, dz: _lattice(c + np.array([dx, dy, dz]), seed)
    x00 = corner(0, 0, 0) + (corner(1, 0, 0) - corner(0, 0, 0)) * u[..., 0]
    x10 = corner(0, 1, 0) + (corner(1, 1, 0) - corner(0, 1, 0)) * u[..., 0]
    x01 = corner(0, 0, 1) + (corner(1, 0, 1) - corner(0, 0, 1)) * u[..., 0]
    x11 = corner(0, 1, 1) + (corner(1, 1, 1) - corner(0, 1, 1)) * u[..., 0]
    y0 = x00 + (x10 - x00) * u[..., 1]
    y1 = x01 + (x11 - x01) * u[..., 1]
    return y0 + (y1 - y0) * u[..., 2]


def cloud_fbm(p, octaves, seed, gain=CLOUD_NOISE['gain'], lacunarity=CLOUD_NOISE['lacunarity']):
    s, a, n = 0.0, 1.0, 0.0
    for i in range(octaves):
        s = s + a * value_noise(p, seed + i * 7919)
        n += a
        a *= gain
        p = p * lacunarity
    return s / n


def cloud_field(n, layer, noise=CLOUD_NOISE):
    """The page's cloud field for unit vectors n (..., 3) in the Moon's frame, before any drift."""
    if layer == 0:
        w = np.stack([value_noise(n * 2.0, 11), value_noise(n * 2.0 + 5.2, 12), value_noise(n * 2.0 + 9.7, 13)], -1) - 0.5
        return cloud_fbm((n + noise['lowWarp'] * w) * noise['lowFrequency'], noise['lowOctaves'], 1)
    q = n * np.array([1.0, noise['highStretch'], 1.0])
    w = np.stack([value_noise(q * 1.5, 21), value_noise(q * 1.5 + 3.1, 22), value_noise(q * 1.5 + 7.3, 23)], -1) - 0.5
    return cloud_fbm((q + noise['highWarp'] * w) * noise['highFrequency'], noise['highOctaves'], 2)


def cloud_quantiles(samples=200_000, seed=3):
    """Quantiles of each deck's field over the sphere, so the page can turn it into a uniform variable and
    threshold it at the run's cover."""
    rng = np.random.default_rng(seed)
    n = rng.normal(size=(samples, 3))
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    q = np.linspace(0.0, 1.0, QUANTILES)
    return [np.quantile(cloud_field(n, layer), q).round(5).tolist() for layer in (0, 1)]


def build(height, level, ppd):
    """Appearance data for the globe page: textures (arrays, 180 W at the left edge) and parameters.

    height: metres above the geoid at ppd pixels per degree, 0 E at the left edge as LOLA has it."""
    moon, engine = engine_atmosphere()
    colour_rel = np.array(moon['sun_colour_linear_srgb'])
    sun_rgb_lux = moon['sun_illuminance_lux'] * colour_rel / float(LUMINANCE @ colour_rel)
    clim = dict(np.load(CLIMATOLOGY))
    clim_meta = json.loads(str(clim.pop('metadata')))
    clim['deck_sigma'] = clim_meta['deck_sigma']
    h = np.roll(height, height.shape[1] // 2, axis=1)
    drainage = drainage_grids(h.shape)
    colour, water, river, surface_level = surface_colour(h, level, clim, ppd, drainage)
    albedo = np.empty(colour.shape, np.uint8)
    for c in range(3):
        albedo[..., c] = np.round(srgb_encode(colour[..., c]) * 255)
    del colour
    half = block_mean(np.maximum(h, surface_level), 2)
    normal = relief_normals(half, level, ppd // 2)
    wet_share = water.astype(np.float32) if river is None else np.maximum(water.astype(np.float32), river)
    water_mask = np.round(wet_share * 255).astype(np.uint8)
    del half, wet_share

    # Cloud: the Sun-relative composites of each deck and its geographic pattern relative to its zonal mean.
    sun = np.dstack([rows_to_uniform(np.clip(clim[k], 0, 1), clim['lat']) for k in ('low_sun', 'high_sun', 'clt_sun')])
    rows = 90.0 - (np.arange(CLOUD_ROWS) + 0.5) * 180.0 / CLOUD_ROWS
    cols = -180.0 + (np.arange(2 * CLOUD_ROWS) + 0.5) * 180.0 / CLOUD_ROWS
    geo = []
    for k in ('low_mean', 'high_mean'):
        zonal = np.maximum(clim[k].mean(axis=1, keepdims=True), 1e-3)
        geo.append(np.clip(regrid(clim[k] / zonal, clim['lat'], clim['lon'], rows, cols), 0.0, 2.55))
    geo = np.dstack([geo[0], geo[1], np.zeros_like(geo[0])])
    low_km, high_km = deck_heights_km(clim, moon['rayleigh_scale_height_km'])
    low_wind, high_wind = deck_mean(clim, clim['ua_layer_mean'])

    direct, diffuse = ground_light_lut(sun_rgb_lux)
    source = diffuse_source_factor(moon, diffuse)
    light = np.stack([direct, diffuse, source]).astype(np.float32)          # 3 rows x elevations x RGB
    textures = dict(albedo=albedo, normal=normal, water=water_mask,
                    cloud_sun=np.round(sun * 255).astype(np.uint8), cloud_geo=np.round(geo * 100).astype(np.uint8))
    area = lambda f: float((f * np.cos(np.radians(clim['lat']))[:, None]).sum()
                           / (np.cos(np.radians(clim['lat'])).sum() * f.shape[1]))
    params = dict(
        tau0=(np.array(moon['rayleigh_scattering_fitted_per_km']) * moon['rayleigh_scale_height_km']).round(4).tolist(),
        scaleHeightKm=moon['rayleigh_scale_height_km'], radiusKm=moon['bottom_radius_km'],
        topKm=moon['bottom_radius_km'] + moon['atmosphere_height_km'], sunColour=moon['sun_colour_linear_srgb'],
        sunLux=moon['sun_illuminance_lux'],
        lutElevations=[float(LUT_ELEVATIONS[0]), float(LUT_ELEVATIONS[-1])],
        light=[[[float(f'{v:.4g}') for v in rgb] for rgb in row] for row in light],
        surfaceKmPerTexel=round(MOON_KM * np.radians(1.0 / ppd), 3), normalStrength=NORMAL_STRENGTH,
        clouds=dict(lowKm=round(low_km, 1), highKm=round(high_km, 1), lowOpacity=LOW_OPACITY, highOpacity=HIGH_OPACITY,
                    albedo=CLOUD_ALBEDO, geoScale=100.0, lowCover=round(area(clim['low_mean']), 3),
                    highCover=round(area(clim['high_mean']), 3), totalCover=round(area(clim['clt_mean']), 3),
                    seaIce=round(float(clim['sic'].max()), 4), snowM=round(float(clim['snd_m'].max()), 4),
                    lowWindMs=round(low_wind, 2), highWindMs=round(high_wind, 2), noise=CLOUD_NOISE,
                    quantiles=cloud_quantiles()),
        climatology=dict(run=clim_meta['run'], years=clim_meta['years'], schema=clim_meta['schema']),
        legend=[['Forest', FOREST], ['Woodland, savanna', SAVANNA], ['Grassland, scrub', GRASS],
                ['Bare highland soil', HIGHLAND_SOIL], ['Bare basaltic soil', BASALT_SOIL], ['Rock', ROCK],
                ['River water, sediment-laden', RIVER_WATER], ['Shallow, silty water', TURBID_WATER],
                ['Deep water', DEEP_WATER]])
    params['legend'] = [[name, '#' + ''.join(f'{int(round(v * 255)):02x}' for v in srgb_encode(np.asarray(rgb)))]
                        for name, rgb in params['legend']]
    sources = dict(engine=ENGINE, sky_atlas=SKY_ATLAS, climatology=CLIMATOLOGY)
    if drainage is not None:
        sources['drainage'] = DRAINAGE
        params['drainage'] = dict(schema=drainage['metadata']['schema'], riverWidthPerRootDischarge=RIVER_WIDTH)
    return textures, params, sources, dict(sea=h < level, water=water)
