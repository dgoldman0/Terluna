"""Data for the globe's appearance mode: surface colour, clouds, ground light and haze.

Called by render.py. Every input is a domain product, and each piece says what kind it is:

- computed: the atmosphere (Rayleigh coefficients and scale height from the illumination domain's
  engine-atmosphere product) and the ground's direct and diffuse light against Sun elevation (the
  Open Moon clear-sky atlas, 55 scattering orders);
- informed: cloud cover from the climate domain's ExoPlaSim climatology (run A), a low and a high deck
  each moved with the Sun through its hour-angle composite and placed at the deck's cover-weighted
  height; light scattered more than once, an isotropic source scaled from the sky atlas's diffuse light
  so that the air seen straight down matches Eddington's conservative-scattering reflectance;
- guesstimate: surface cover from rainfall, soil moisture, nearness to water, height and slope (forest,
  woodland, grassland, bare highland or basaltic soil, rock), water colour from depth, and the cloud
  decks' opacity. The biosphere and climate work replace it.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / 'illumination' / 'sky' / 'products' / 'engine_atmosphere.json'
SKY_ATLAS = ROOT / 'illumination' / 'sky' / 'data' / 'moon_atlas.npz'
CLIMATOLOGY = ROOT / 'climate' / 'gcm' / 'products' / 'climatology_A.npz'
LUMINANCE = np.array([0.2126, 0.7152, 0.0722])
LUT_ELEVATIONS = np.linspace(-90.0, 90.0, 121)   # twilight light reaches far past the terminator

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
CLOUD_ALBEDO = 0.72                 # of the part of a deck that intercepts light
LOW_OPACITY, HIGH_OPACITY = 0.9, 0.4  # guesstimates: the low deck holds the model's cloud water, the high deck little
CLOUD_ROWS = 90                     # 2-degree rows for the cloud maps
CALIBRATION_FLOOR_DEG = 10.0        # below this Sun height the scattering calibration holds its 10-degree value


def smoothstep(a, b, x):
    t = np.clip((x - a) / (b - a), 0.0, 1.0)
    return t * t * (3 - 2 * t)


def fbm(shape, seed, octaves=7, base=(6, 3), gain=0.55, stretch=1.0):
    """Fractal value noise, seamless in longitude; stretch > 1 elongates features east-west."""
    rng = np.random.default_rng(seed)
    h, w = shape
    out = np.zeros(shape)
    amp, total = 1.0, 0.0
    for o in range(octaves):
        nx, ny = int(base[0] * 2 ** o), max(2, int(base[1] * 2 ** o * stretch))
        grid = rng.random((ny + 1, nx))
        y = np.linspace(0, ny, h, endpoint=False)[:, None]
        x = np.linspace(0, nx, w, endpoint=False)[None, :]
        y0, x0 = np.floor(y).astype(int), np.floor(x).astype(int)
        fy, fx = y - y0, x - x0
        fy, fx = fy * fy * (3 - 2 * fy), fx * fx * (3 - 2 * fx)
        y1, x1 = np.minimum(y0 + 1, ny), (x0 + 1) % nx
        x0 = x0 % nx
        v = (grid[y0, x0] * (1 - fx) + grid[y0, x1] * fx) * (1 - fy) + (grid[y1, x0] * (1 - fx) + grid[y1, x1] * fx) * fy
        out += amp * v
        total += amp
        amp *= gain
    return out / total


def uniformise(x):
    """Rank-transform to a uniform 0..1 distribution, so a threshold at 1 - c covers a share c."""
    r = np.empty(x.size)
    r[np.argsort(x, axis=None)] = np.linspace(0.0, 1.0, x.size)
    return r.reshape(x.shape)


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


def surface_colour(height, level, lat, lon, clim):
    """Guesstimated surface reflectance (linear sRGB) and the water mask, 180 W at the left edge."""
    depth = level - height
    water = depth > 0
    from scipy.ndimage import distance_transform_edt
    rain = regrid(clim['pr_mm_day'], clim['lat'], clim['lon'], lat, lon)
    soil = regrid(clim['mrso_m'], clim['lat'], clim['lon'], lat, lon)
    patch = fbm(height.shape, 11, octaves=6, base=(24, 12)) - 0.5
    km_per_px = 1737.4 * np.radians(180.0 / height.shape[0])
    tiled = np.concatenate([water, water, water], axis=1)
    shore_km = distance_transform_edt(~tiled)[:, height.shape[1]:2 * height.shape[1]] * km_per_px
    coast = np.exp(-shore_km / 120.0)
    lowland = np.exp(-np.clip(height - level, 0, None) / 400.0)
    wet = rain + 1.5 * smoothstep(0.02, 0.2, soil) + 1.4 * coast + 0.8 * lowland + 2.5 * patch
    veg = smoothstep(0.4, 3.2, wet)
    dense = smoothstep(3.2, 5.5, wet)
    basaltic = smoothstep(level + 900.0, level + 150.0, height) * smoothstep(0.35, 0.65, fbm(height.shape, 17, octaves=5, base=(18, 9)))
    bare = HIGHLAND_SOIL * (1 - basaltic[..., None]) + BASALT_SOIL * basaltic[..., None]
    gy, gx = np.gradient(height)
    dy = 1737.4e3 * np.radians(180.0 / height.shape[0])
    slope = np.hypot(gx / (dy * np.maximum(np.cos(np.radians(lat))[:, None], 0.05)), gy / dy)
    rocky = smoothstep(0.16, 0.34, slope)[..., None]
    grassy = bare * (1 - smoothstep(0.0, 0.35, veg))[..., None] + GRASS * smoothstep(0.0, 0.35, veg)[..., None]
    wooded = SAVANNA * (1 - dense[..., None]) + FOREST * dense[..., None]
    land = grassy * (1 - smoothstep(0.35, 1.0, veg))[..., None] + wooded * smoothstep(0.35, 1.0, veg)[..., None]
    land = land * (1 - 0.5 * rocky) + ROCK * 0.5 * rocky
    shallow = smoothstep(250.0, 30.0, depth)[..., None]
    deep = smoothstep(400.0, 1500.0, depth)[..., None]
    sea = MID_WATER * (1 - deep) + DEEP_WATER * deep
    sea = sea * (1 - shallow) + TURBID_WATER * shallow
    colour = np.where(water[..., None], sea, land)
    return colour, water


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
    sigma, cover = clim['sigma'], clim['cl_layer_mean']
    z = -h_km * np.log(sigma)
    low = sigma > clim['deck_sigma']
    return float((z * cover)[low].sum() / cover[low].sum()), float((z * cover)[~low].sum() / cover[~low].sum())


def build(height, level, lat, ppd):
    """Appearance data for the globe page: textures (arrays) and parameters (JSON-ready)."""
    moon, engine = engine_atmosphere()
    colour_rel = np.array(moon['sun_colour_linear_srgb'])
    sun_rgb_lux = moon['sun_illuminance_lux'] * colour_rel / float(LUMINANCE @ colour_rel)
    clim = dict(np.load(CLIMATOLOGY))
    clim_meta = json.loads(str(clim.pop('metadata')))
    clim['deck_sigma'] = clim_meta['deck_sigma']
    lon = (np.arange(height.shape[1]) + 0.5) / ppd - 180.0
    h = np.roll(height, height.shape[1] // 2, axis=1)
    colour, water = surface_colour(h, level, lat, lon, clim)
    albedo = np.round(srgb_encode(colour) * 255).astype(np.uint8)

    # Cloud: per-deck noise (low decks lumpy, high decks drawn out east-west), the Sun-relative
    # composites and the geographic pattern of each deck relative to its zonal mean.
    low_noise = uniformise(fbm(h.shape, 5, octaves=8, base=(24, 12), gain=0.6))
    high_noise = uniformise(fbm(h.shape, 7, octaves=7, base=(10, 10), gain=0.62, stretch=2.2))
    noise = np.dstack([low_noise, high_noise, water.astype(float)])      # blue: the water mask, for sunglint
    sun = np.dstack([rows_to_uniform(np.clip(clim[k], 0, 1), clim['lat']) for k in ('low_sun', 'high_sun', 'clt_sun')])
    rows = 90.0 - (np.arange(CLOUD_ROWS) + 0.5) * 180.0 / CLOUD_ROWS
    cols = -180.0 + (np.arange(2 * CLOUD_ROWS) + 0.5) * 180.0 / CLOUD_ROWS
    geo = []
    for k in ('low_mean', 'high_mean'):
        zonal = np.maximum(clim[k].mean(axis=1, keepdims=True), 1e-3)
        geo.append(np.clip(regrid(clim[k] / zonal, clim['lat'], clim['lon'], rows, cols), 0.0, 2.55))
    geo = np.dstack([geo[0], geo[1], np.zeros_like(geo[0])])
    low_km, high_km = deck_heights_km(clim, moon['rayleigh_scale_height_km'])

    direct, diffuse = ground_light_lut(sun_rgb_lux)
    source = diffuse_source_factor(moon, diffuse)
    light = np.stack([direct, diffuse, source]).astype(np.float32)          # 3 rows x elevations x RGB
    textures = dict(albedo=albedo, cloud_noise=np.round(noise * 255).astype(np.uint8),
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
        clouds=dict(lowKm=round(low_km, 1), highKm=round(high_km, 1), lowOpacity=LOW_OPACITY, highOpacity=HIGH_OPACITY,
                    albedo=CLOUD_ALBEDO, geoScale=100.0, lowCover=round(area(clim['low_mean']), 3),
                    highCover=round(area(clim['high_mean']), 3), totalCover=round(area(clim['clt_mean']), 3),
                    seaIce=round(float(clim['sic'].max()), 4), snowM=round(float(clim['snd_m'].max()), 4)),
        climatology=dict(run=clim_meta['run'], years=clim_meta['years'], schema=clim_meta['schema']),
        legend=[['Forest', FOREST], ['Woodland, savanna', SAVANNA], ['Grassland, scrub', GRASS],
                ['Bare highland soil', HIGHLAND_SOIL], ['Bare basaltic soil', BASALT_SOIL], ['Rock', ROCK],
                ['Shallow, silty water', TURBID_WATER], ['Deep water', DEEP_WATER]])
    params['legend'] = [[name, '#' + ''.join(f'{int(round(v * 255)):02x}' for v in srgb_encode(np.asarray(rgb)))]
                        for name, rgb in params['legend']]
    return textures, params, dict(engine=ENGINE, sky_atlas=SKY_ATLAS, climatology=CLIMATOLOGY)
