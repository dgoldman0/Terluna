"""Data for the globe's appearance mode: surface colour, relief, clouds, ground light and haze.

Called by render.py. Every input is a domain product, and each piece says what kind it is:

- computed: the atmosphere (Rayleigh coefficients and scale height from the illumination domain's
  engine-atmosphere product), the ground's direct and diffuse light against Sun elevation (the
  Open Moon clear-sky atlas, 55 scattering orders), and land, sea and relief at 16 pixels per degree
  (LOLA above the GRAIL geoid through the geography domain, at the atlas product's sea level), with
  the geography domain's rivers and rain-fed lakes above sea level;
- informed: cloud cover from the climate domain's ExoPlaSim climatology (the design case, A28_corrected_970
  years 10-19, since 16 October; run A before), a low and a high deck
  each moved with the Sun through its hour-angle composite, placed at the deck's cover-weighted
  height and carried east by the deck's cover-weighted mean wind; light scattered more than once, an
  isotropic source scaled from the sky atlas's diffuse light so that the air seen straight down
  matches Eddington's conservative-scattering reflectance;
- guesstimate: surface cover from rainfall, soil moisture, nearness to water and rivers, height, hollows
  and slope (forest, woodland, grassland, bare highland or basaltic soil, rock), river width from discharge,
  water colour from depth, and the
  cloud decks' shapes (noise the page evaluates, thresholded to the run's cover) and opacity. The
  biosphere and climate work replace it.

On the turned Moon (7 November 2026) render.py hands build() the heights already carried into the turned frame
(climate.gcm.turned's rotation), a climatology of a run on the turned Moon, which is in that frame already, and the
rotation, with which drainage_grids carries the rivers and lakes, routed on today's ground, into the same frame.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / 'illumination' / 'sky' / 'products' / 'engine_atmosphere.json'
SKY_ATLAS = ROOT / 'illumination' / 'sky' / 'data' / 'moon_atlas.npz'
CLIMATOLOGY = ROOT / 'climate' / 'gcm' / 'products' / 'climatology_A28_corrected_970.npz'   # the design case
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
# True colour (the design surface's globe, --moisture). The cover follows the design climate's dryness: the moisture
# map's aridity, the run's rain over its demand, in each of the run's land cells. The colours are Earth's surfaces
# under the same dryness, as reflectance in linear sRGB; a guesstimate, labelled so on the page.
MOISTURE = ROOT / 'research' / 'studies' / 'atlas_round' / 'results' / 'moisture_map_final.json'
MOISTURE_SCHEMA = 'terluna.biosphere.ecology-moisture-map/1'
WATER_COLOURS = ROOT / 'research' / 'studies' / 'sea_appearance' / 'results' / 'water_colours.json'
GREEN_HERB = np.array([0.050, 0.095, 0.035])    # green grass and shrubs
DRY_HERB = np.array([0.20, 0.17, 0.10])         # grass gone to straw
BASALT_ROCK = np.array([0.08, 0.075, 0.07])
# Weathered soil. Under the air space weathering ends, and oxygen and water rust the rock's iron, so the bare ground
# takes the colours of Earth's soils on like rock under like dryness. The terrain's upkeep (register, 13 October)
# carries rock dust back uphill at 1-42 mm a century, renewing the top metre every 2,400 to 100,000 years, so the
# topsoil stays moderately weathered. Each end of the rock's iron has a typical dry soil colour (Munsell notation)
# at three aridities, with its xyY under Illuminant C from the Munsell renotation (Newhall, Nickerson and Judd 1943;
# RIT's real.dat). The feldspathic highland's iron, 4-6% FeO, matches Earth's upper continental crust (about 5%), so
# its soils redden as Earth's do.
SOIL_FEO = (4.5, 20.0)               # FeO (weight per cent) of the highland and the mare basalt ends
SOIL_MUNSELL = dict(
    highland=((0.05, '7.5YR 6/4', 'light brown, a desert soil with lime in it', (0.3860, 0.3652, 30.05)),
              (0.3, '7.5YR 5/4', 'brown, a grassland soil with some organic matter', (0.3991, 0.3714, 19.77)),
              (1.0, '5YR 5/6', 'yellowish red, a warm humid soil', (0.4420, 0.3808, 19.77))),
    basalt=((0.05, '5YR 4/4', 'reddish brown, a desert soil on basalt', (0.4187, 0.3679, 12.00)),
            (0.3, '5YR 3/4', 'dark reddish brown, a dry grassland soil on basalt', (0.4376, 0.3715, 6.555)),
            (1.0, '2.5YR 3/6', 'dark red, terra roxa on basalt', (0.4954, 0.3692, 6.555))))
TITANIUM_DARKENING = 0.97            # per weight per cent of TiO2 above 2% in the basalt: ilmenite, opaque, outlasts weathering
# Plant cover against aridity: least + (most - least) / (1 + (half / aridity)^steep). Earth's drylands lose their
# cover fast below an aridity of 0.2 (Berdugo et al. 2020, Science 367, 787) and its humid forests close their
# canopy; the curve gives 2.5% at 0.02, 17% at 0.1, 46% at 0.2, 72% at 0.33, 86% at 0.5, 92% at 0.65 and 97% at 1.
COVER = dict(least=0.02, most=1.0, half=0.22, steep=2.2)
TREES = dict(half=0.6, steep=3.0)    # the trees' share of the cover: 11% at 0.3, 56% at 0.65 (woodland), 82% at 1
GREEN = dict(half=0.35, steep=2.0)   # the green share of the grass: a quarter at 0.2, half at 0.35, 89% at 1
# Within the climate's cells (display rules, guesstimates): ground along larger rivers and in hollows wetter, steep
# slopes drier, all as a factor on the aridity, and patchiness drawn as noise.
FLOODPLAIN_LN = 0.6                  # times the blurred log10 discharge, at most 2: up to 3.3 times as wet
HOLLOW_LN = 0.15                     # times the standardized curvature, clipped at 2.5
SLOPE_LN = 0.4                       # taken off where the ground steepens from 0.2 to 0.6
PATCH_LN = 0.25                      # the noise's standard deviation in log aridity
SHELF_M = (150.0, 20.0)              # water shallower than these takes the productive coast's colour, fully at 20 m
PLUME_KM_PER_ROOT = 0.3              # a plume reaches about 0.3 km per (m3/s)^0.5 of sediment-bearing flow: 9.5 km at 1,000 m3/s
# Falls (the falls study, research/studies/atlas_round/falls.py), display guesstimates: a river's water turns to foam as
# its step steepens from FOAM_SLOPE[0] to FOAM_SLOPE[1], and spray whitens the ground and water round the falls of
# MIST_GW or more, out to about MIST_KM per cube root of their GW (5 km round 42 GW), at most MIST_MAX of each pixel.
WHITE_WATER = np.array([0.55, 0.58, 0.60])
MIST_WHITE = np.array([0.62, 0.64, 0.66])
FOAM_SLOPE = (0.03, 0.2)
MIST_GW = 2.0
MIST_KM = 1.5
MIST_MAX = 0.5
CLOUD_ALBEDO = 0.72                 # of the part of a deck that intercepts light
LOW_OPACITY, HIGH_OPACITY = 0.9, 0.4  # guesstimates: the low deck holds the model's cloud water, the high deck little
CLOUD_ROWS = 90                     # 2-degree rows for the cloud maps
# The climate model's fields are smoothed over this share of its grid spacing before they set the ground cover
# (a display choice): bilinear interpolation of its 5.5-degree cells leaves kinks along the grid lines, and its
# grid-scale ripples, finer than the spectral model resolves, crossed the cover's thresholds as straight east-west
# and north-south edges. At 0.7 of a cell, ripples two cells long keep a tenth of their contrast and waves 30 degrees
# long about 70%.
CLIMATE_SMOOTHING_CELLS = 0.7
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


def table_noise(table, x, y, z):
    """Value noise in 3D from a 256-cubed table of random values, wrapping, with the quintic fade the page uses."""
    flat = table.reshape(-1)
    xi, yi, zi = np.floor(x), np.floor(y), np.floor(z)
    fade = lambda t: (t * t * t * (t * (t * 6.0 - 15.0) + 10.0)).astype(np.float32)
    u, v, w = fade(x - xi), fade(y - yi), fade(z - zi)
    i0, j0, k0 = (xi.astype(np.int64) & 255), (yi.astype(np.int64) & 255), (zi.astype(np.int64) & 255)
    i0, i1 = i0 << 16, ((i0 + 1) & 255) << 16
    j0, j1 = j0 << 8, ((j0 + 1) & 255) << 8
    k1 = (k0 + 1) & 255
    lerp = lambda a, b, t: a + t * (b - a)
    near = lerp(lerp(flat[i0 + j0 + k0], flat[i1 + j0 + k0], u), lerp(flat[i0 + j1 + k0], flat[i1 + j1 + k0], u), v)
    far = lerp(lerp(flat[i0 + j0 + k1], flat[i1 + j0 + k1], u), lerp(flat[i0 + j1 + k1], flat[i1 + j1 + k1], u), v)
    return lerp(near, far, w)


def sphere_fbm(shape, seed, octaves, base, spread, gain=0.55):
    """Fractal value noise evaluated in 3D at the sphere's points, so that it is round everywhere, the poles
    included, and seamless in longitude: base cells round the equator at the first octave and twice as many at
    each next, down to cells two pixels across. Float32 with mean 0.5 and standard deviation `spread`."""
    rng = np.random.default_rng(seed)
    table = rng.random((256, 256, 256), dtype=np.float32)
    rows, cols = shape
    ppd = cols / 360.0
    freqs = [base / (2.0 * np.pi) * 2 ** o for o in range(octaves)]
    freqs = [f for f in freqs if np.degrees(1.0 / f) * ppd >= 2.0]
    shifts = rng.random((len(freqs), 3)) * 256.0
    amps = gain ** np.arange(len(freqs))
    lon = np.radians(-180.0 + (np.arange(cols) + 0.5) / ppd)
    out = np.empty(shape, np.float32)
    for r0 in range(0, rows, 96):
        lat = np.radians(90.0 - (np.arange(r0, min(rows, r0 + 96)) + 0.5) / ppd)[:, None]
        x, y, z = np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat) + 0.0 * lon
        acc = np.zeros(x.shape)
        for f, (a, b, c), amp in zip(freqs, shifts, amps):
            acc += amp * table_noise(table, x * f + a, y * f + b, z * f + c)
        out[r0:r0 + 96] = acc
    return (0.5 + (out - out.mean()) * (spread / out.std())).astype(np.float32)


def sphere_blur(a, sigma, lat):
    """A Gaussian blur that is round on the ground: sigma rows north-south and, in each row, sigma / cos(latitude)
    columns east-west, wrapping round, so that it keeps its shape up to the poles. In float64."""
    from scipy import fft
    from scipy.ndimage import gaussian_filter1d
    out = gaussian_filter1d(a.astype(np.float64), sigma, axis=0, mode='nearest')
    across = sigma / np.maximum(np.cos(np.radians(lat)), 1e-3)
    response = np.exp(-2.0 * np.pi ** 2 * (across[:, None] * np.fft.rfftfreq(a.shape[1])[None, :]) ** 2)
    return fft.irfft(fft.rfft(out, axis=1) * response, n=a.shape[1], axis=1)


def sphere_laplacian(a, lat):
    """The Laplacian on the ground in pixel units, from a field blurred with sphere_blur: the north-south second
    difference plus the east-west one scaled by 1 / cos^2(latitude) for the columns' closer spacing."""
    c2 = np.maximum(np.cos(np.radians(lat)), 1e-3)[:, None] ** 2
    north, south = np.vstack([a[:1], a[:-1]]), np.vstack([a[1:], a[-1:]])
    return (north + south - 2.0 * a) + (np.roll(a, 1, axis=1) + np.roll(a, -1, axis=1) - 2.0 * a) / c2


def shore_distance_km(water, ppd):
    """Distance on the ground to the nearest water, in km: the nearest 0.25-degree cell holding any water, found
    among the cells' points on the sphere (so that it is round up to the poles) and spread back to the grid."""
    import cv2
    from scipy.spatial import cKDTree
    k = max(1, ppd // 4)
    rows, cols = water.shape[0] // k, water.shape[1] // k
    wet = water[:rows * k, :cols * k].reshape(rows, k, cols, k).any(axis=(1, 3))
    lat = np.radians(90.0 - (np.arange(rows) + 0.5) * k / ppd)[:, None]
    lon = np.radians(-180.0 + (np.arange(cols) + 0.5) * k / ppd)[None, :]
    points = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat) + 0.0 * lon], axis=-1)
    chord, _ = cKDTree(points[wet]).query(points[~wet])
    dist = np.zeros((rows, cols), np.float32)
    dist[~wet] = 2.0 * np.arcsin(np.minimum(chord / 2.0, 1.0)) * MOON_KM
    pad = np.concatenate([dist[:, -1:], dist, dist[:, :1]], axis=1)
    big = cv2.resize(pad, (pad.shape[1] * k, water.shape[0]), interpolation=cv2.INTER_LINEAR)
    return big[:, k:k + water.shape[1]]


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


def gradients(a, dy, lat):
    """North-south and east-west slopes per metre of ground, central and wrapping in longitude. Each east-west
    difference spans about the same ground as a north-south one, two rows: towards the poles, where the grid's
    columns crowd together, it reaches across as many columns as make that up (at most an eighth of the way
    round), so that noise between near neighbours does not swell into spokes."""
    gy = np.gradient(a, axis=0) / dy
    cos = np.maximum(np.cos(np.radians(lat)), 1e-6)
    gx = np.empty(a.shape, np.float32)
    for i, (row, c) in enumerate(zip(a, cos)):
        m = int(np.clip(np.round(1.0 / c), 1, a.shape[1] // 8))
        gx[i] = (np.roll(row, -m) - np.roll(row, m)) / (2.0 * m * dy * c)
    return gy, gx


def climate_field(clim, name, shape):
    """A climatology field on the texture grid (grid_field)."""
    return grid_field(clim[name], clim['lat'], clim['lon'], shape)


def turn_texture(a, rotation, left_lon=-180.0, nearest=False, chunk=256):
    """An equirectangular texture (rows from 90 N, columns eastward from left_lon) carried into the turned frame: each
    pixel of the result takes the texture at its place today, bilinear or nearest."""
    import cv2
    from climate.gcm.turned import to_today
    rows, cols = a.shape[:2]
    ppd, pad = cols / 360.0, 2
    extra = [(0, 0)] * (a.ndim - 2)
    src = np.pad(np.pad(a, [(0, 0), (pad, pad)] + extra, mode='wrap'), [(pad, pad), (0, 0)] + extra, mode='edge')
    lon_t = left_lon + (np.arange(cols) + 0.5) / ppd
    out = np.empty_like(a)
    for r in range(0, rows, chunk):
        lat_t = 90.0 - (np.arange(r, min(r + chunk, rows)) + 0.5) / ppd
        lat0, lon0 = to_today(*np.meshgrid(lat_t, lon_t, indexing='ij'), rotation)
        map_y = ((90.0 - lat0) * ppd - 0.5 + pad).astype(np.float32)
        map_x = (((lon0 - left_lon) % 360.0) * ppd - 0.5 + pad).astype(np.float32)
        out[r:r + chunk] = cv2.remap(src, map_x, map_y, cv2.INTER_NEAREST if nearest else cv2.INTER_LINEAR,
                                     borderMode=cv2.BORDER_REPLICATE)
    return out


def turn_channels(q, rotation, left_lon=-180.0):
    """Channel discharge carried into the turned frame: every channel cell lands in the turned pixel that holds its
    centre, and every turned pixel whose centre lies in a channel cell today takes that cell, keeping the larger
    discharge, so channels stay continuous where the turned pixels are smaller or larger than today's."""
    from climate.gcm.turned import to_turned
    rows, cols = q.shape
    ppd = cols / 360.0
    r, c = np.nonzero(q > 0)
    lat_t, lon_t = to_turned(90.0 - (r + 0.5) / ppd, left_lon + (c + 0.5) / ppd, rotation)
    rt = np.clip(np.floor((90.0 - lat_t) * ppd).astype(int), 0, rows - 1)
    ct = np.floor(((lon_t - left_lon) % 360.0) * ppd).astype(int) % cols
    out = turn_texture(q, rotation, left_lon, nearest=True)
    np.maximum.at(out, (rt, ct), q[r, c])
    return out


def drainage_grids(shape, path=DRAINAGE, rotation=None):
    """Discharge (m3/s) along channels and the level of rain-fed lakes above sea level on the texture grid
    (180 W at the left edge), from the geography domain's drainage product, carried into the turned frame when a
    rotation is given; None when the product is absent."""
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
    q, lake = np.roll(q.reshape(shape), roll, axis=1), np.roll(lake.reshape(shape), roll, axis=1)
    if rotation is not None:
        q, lake = turn_channels(q, rotation), turn_texture(lake, rotation, nearest=True)
    return dict(discharge=q, lake_level=lake, metadata=meta)


def land_discharge(routing_path):
    """The discharge each cell carries off land since the last lake or sea (m3/s, flat on the routing's grid): lakes
    and seas trap their rivers' sediment, so their outflows run clear and only the runoff gathered below them carries
    sediment on. From the routing product's receivers and discharge: each dry cell's own runoff is its discharge less
    what flows into it, gathered downstream as the routing gathers discharge, with every water cell passing none on.
    The near-side sea's spill, which the routing puts in at the cell past the neck's crest, is sea water and carries
    none."""
    with np.load(routing_path, allow_pickle=False) as z:
        receiver = z['receiver'].astype(np.int64)
        q = z['discharge_m3s']
        water = (z['sea'] >= 0) | z['routing_lake'] | z['outlet']
        water[z['lake_cell']] = True
        neck = json.loads(str(z['metadata'])).get('neck')
    n = receiver.size
    donor = receiver != np.arange(n)
    inflow = np.bincount(receiver[donor], weights=q[donor], minlength=n)
    if neck:
        inflow[neck['spill_cell']] += neck['flow_m3s']
    carried = np.where(water, 0.0, np.maximum(q - inflow, 0.0))
    waiting = np.bincount(receiver[donor], minlength=n)
    front = np.flatnonzero((waiting == 0) & donor)
    while front.size:
        to = receiver[front]
        np.add.at(carried, to, np.where(water[front], 0.0, carried[front]))
        np.subtract.at(waiting, to, 1)
        to = np.unique(to)
        front = to[(waiting[to] == 0) & donor[to]]
    carried[water] = 0.0
    return carried


def surface_grids(shape, path, rotation=None, routing=None):
    """Every water body's surface and the rivers' discharge on the texture grid (180 W at the left edge) from a
    surface that holds each water cell's level (the atlas round's globe surface), carried into the turned frame when a
    rotation is given: the drainage product's form, with every sea, lake and parted part of a sea as a water body at
    its own level. With the routing product the surface came from, also each channel's discharge off land since the
    last lake or sea (land_discharge), as `sediment`."""
    z = np.load(path)
    meta = json.loads(str(z['metadata']))
    if meta['grid']['rows'] != shape[0] or meta['grid']['cols'] != shape[1]:
        raise SystemExit(f'Expected a {shape[0]}-row grid in {path}')
    q = np.zeros(shape[0] * shape[1], np.float32)
    q[z['channel_cell']] = z['channel_discharge_m3s']
    level = np.full(shape[0] * shape[1], np.nan, np.float32)
    level[z['water_cell']] = z['water_level_m']
    roll = shape[1] // 2
    q, level = np.roll(q.reshape(shape), roll, axis=1), np.roll(level.reshape(shape), roll, axis=1)
    if rotation is not None:
        q, level = turn_channels(q, rotation), turn_texture(level, rotation, nearest=True)
    out = dict(discharge=q, lake_level=level, metadata=meta)
    if routing is not None:
        s = np.zeros(shape[0] * shape[1], np.float32)
        s[z['channel_cell']] = land_discharge(routing)[z['channel_cell']]
        s = np.roll(s.reshape(shape), roll, axis=1)
        out['sediment'] = turn_channels(s, rotation) if rotation is not None else s
    return out


def surface_colour(height, level, clim, ppd, drainage=None, sea=None):
    """Guesstimated surface reflectance (linear sRGB, float32), the water mask (sea and lakes) and the share
    of each cell that river channels cover, 180 W at the left edge. sea, when given, replaces the ground below
    `level` as the sea (an empty mask where every water body comes as a lake at its own level)."""
    height = height.astype(np.float32)
    rows, cols = height.shape
    lat = 90.0 - (np.arange(rows) + 0.5) / ppd
    km = MOON_KM * np.radians(1.0 / ppd)
    sea = height < level if sea is None else sea
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
    wet += 1.4 * np.exp(-shore_distance_km(water, ppd) / 120.0)
    land_height = np.maximum(height, level)
    wet += 0.8 * np.exp(-(land_height - level) / 400.0)
    curvature = sphere_laplacian(sphere_blur(land_height, 3.0, lat), lat).astype(np.float32)
    wet += 0.6 * np.clip(curvature / (curvature[~water].std() + 1e-6), -2.5, 2.5)   # wetter in hollows
    del curvature
    # The noise keeps the spread (0.146) the earlier grid noise had, so that the cover keeps its balance.
    wet += 2.5 * (sphere_fbm(height.shape, 11, octaves=9, base=24, spread=0.146) - 0.5)
    if discharge is not None:                        # floodplains: wetter ground along larger rivers
        wet += 0.9 * np.minimum(sphere_blur(np.log10(1.0 + discharge), 2.5, lat).astype(np.float32), 2.0)
    veg = smoothstep(0.4, 3.2, wet)
    dense = smoothstep(3.2, 5.5, wet)
    del wet
    basaltic = smoothstep(level + 900.0, level + 150.0, height) * smoothstep(
        0.35, 0.65, sphere_fbm(height.shape, 17, octaves=8, base=18, spread=0.127))
    gy, gx = gradients(land_height, km * 1000.0, lat)
    rocky = 0.5 * smoothstep(0.35, 0.7, np.hypot(gx, gy))
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


def plant_cover(aridity):
    """The share of the ground plants cover at an aridity (rain over demand), by COVER."""
    c = COVER
    return c['least'] + (c['most'] - c['least']) / (1.0 + (c['half'] / aridity) ** c['steep'])


def share_of(aridity, curve):
    """A share rising with aridity through half at curve['half']."""
    return 1.0 / (1.0 + (curve['half'] / aridity) ** curve['steep'])


XYZ_TO_SRGB = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])
BRADFORD = np.array([[0.8951, 0.2664, -0.1614], [-0.7502, 1.7135, 0.0367], [0.0389, -0.0685, 1.0296]])


def munsell_srgb(xyY):
    """A Munsell renotation colour's reflectance in linear sRGB: its xyY under Illuminant C, Y scaled by 0.975 to the
    perfect diffuser (the renotation's Y is relative to magnesium oxide), carried to D65 by the Bradford transform."""
    x, y, Y = xyY
    Y = Y * 0.975 / 100.0
    white = lambda wx, wy: np.array([wx / wy, 1.0, (1.0 - wx - wy) / wy])
    adapt = np.linalg.inv(BRADFORD) @ np.diag((BRADFORD @ white(0.31271, 0.32902)) / (BRADFORD @ white(0.31006, 0.31616))) @ BRADFORD
    return np.clip(XYZ_TO_SRGB @ (adapt @ np.array([x / y * Y, Y, (1.0 - x - y) / y * Y])), 0.0, None)


def soil_colour(aridity, feo, tio2=0.0):
    """Bare weathered ground (linear sRGB, float32, a trailing axis of 3): each end's soil colour interpolated in log
    aridity between its three (held beyond them), the two ends mixed by the iron (FeO, weight per cent) between theirs,
    and the basalt end darkened by its titanium."""
    ln = np.log(np.asarray(aridity, np.float32))
    ends = []
    for end in ('highland', 'basalt'):
        knots = SOIL_MUNSELL[end]
        at = np.log([k[0] for k in knots])
        rgb = np.array([munsell_srgb(k[3]) for k in knots])
        ends.append(np.stack([np.interp(ln, at, rgb[:, c]) for c in range(3)], axis=-1).astype(np.float32))
    f = np.clip((np.asarray(feo, np.float32) - SOIL_FEO[0]) / (SOIL_FEO[1] - SOIL_FEO[0]), 0.0, 1.0)[..., None]
    dark = (TITANIUM_DARKENING ** np.maximum(np.asarray(tio2, np.float32) - 2.0, 0.0))[..., None]
    return (ends[0] * (1 - f) + ends[1] * dark * f).astype(np.float32)


def ground_colour(aridity, soil, mafic=0.0, rocky=0.0):
    """Reflectance of the ground (linear sRGB, float32, a trailing axis of 3) at an aridity: plants over bare ground in
    the share Earth's land holds under the same dryness, the trees' share and the grass's green rising with the rain;
    the bare ground the weathered soil (soil_colour), and rock in its place on steep ground (rocky the share of the
    bare ground; mafic the basalt's share of the rock)."""
    a = np.asarray(aridity, np.float32)[..., None]
    m = np.asarray(mafic, np.float32)[..., None]
    r = np.asarray(rocky, np.float32)[..., None]
    bare = np.asarray(soil, np.float32) * (1 - r) + (ROCK * (1 - m) + BASALT_ROCK * m) * r
    green = share_of(a, GREEN)
    trees = share_of(a, TREES)
    canopy = FOREST * trees + (DRY_HERB * (1 - green) + GREEN_HERB * green) * (1 - trees)
    cover = plant_cover(a)
    return (bare * (1 - cover) + canopy * cover).astype(np.float32)


def grid_field(values, lat_src, lon_src, shape):
    """A field on the climate model's grid carried to the texture grid: regridded at 4 pixels per degree, smoothed
    over CLIMATE_SMOOTHING_CELLS of the model's grid spacing, then enlarged."""
    import cv2
    lat4 = 90.0 - (np.arange(720) + 0.5) / 4.0
    lon4 = (np.arange(1440) + 0.5) / 4.0 - 180.0
    f = regrid(values, lat_src, lon_src, lat4, lon4)
    spacing = float(np.abs(np.diff(lat_src)).mean())
    f = sphere_blur(f, CLIMATE_SMOOTHING_CELLS * spacing * 4.0, lat4).astype(np.float32)
    s = shape[1] // f.shape[1]
    pad = np.concatenate([f[:, -1:], f, f[:, :1]], axis=1)
    return cv2.resize(pad, (pad.shape[1] * s, shape[0]), interpolation=cv2.INTER_LINEAR)[:, s:s + shape[1]]


def moisture_record(path=MOISTURE):
    record = json.loads(Path(path).read_text())
    if record.get('schema') != MOISTURE_SCHEMA:
        raise SystemExit(f'Expected {MOISTURE_SCHEMA} in {path}')
    return record


def moisture_field(record, shape):
    """The moisture map's log aridity on the texture grid: each of the run's land cells at its centre, every other cell
    given its nearest land cell's (as the map assigns the design surface's land there), then regridded and smoothed as
    the climate's fields are. Both are on the run's grid, in the globe's frame."""
    from scipy.spatial import cKDTree
    lat = np.asarray(record['grid']['lat'], float)
    lon = np.asarray(record['grid']['lon'], float)
    cells = record['cells']
    ln = np.full((lat.size, lon.size), np.nan)
    ln[np.asarray(cells['row']), np.asarray(cells['col'])] = np.log(np.asarray(cells['aridity'], float))
    la, lo = np.meshgrid(np.radians(lat), np.radians(lon), indexing='ij')
    points = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)
    land = np.isfinite(ln)
    _, nearest = cKDTree(points[land]).query(points[~land])
    ln[~land] = ln[land][nearest]
    return grid_field(ln, lat, lon, shape)


def water_albedos(path=WATER_COLOURS):
    """Each designed water's reflectance (linear sRGB) as the page lights it: pi times its remote-sensing reflectance,
    weighted in each channel by sunlight (the illumination domain's solar spectrum, the spectrum that sets the page's
    sunlight), for the open sea, the productive coast and the river mouth at their mean chlorophyll."""
    from illumination.sky.atmospheres import SW, SF, MOON
    from illumination.sky.colour_matching import colour_weights
    d = json.loads(Path(path).read_text())
    lam = np.asarray(d['display_wavelength_nm'], float)
    weights = (np.interp(lam, SW, SF) * MOON.visible_filter)[:, None] * colour_weights(lam)
    case = lambda water: next(c for c in d['cases']
                              if c['water'] == water and c['phase'] == 'mean' and c['soil'] == 'developed')
    return {k: (np.pi * np.asarray(case(k)['rrs_display_per_sr']) @ weights / weights.sum(axis=0)).astype(np.float32)
            for k in ('open_sea', 'productive_coast', 'river_mouth')}


def plume_share(discharge, water, lat, km):
    """The share of each water pixel a river's sediment plume colours: round each mouth (a channel pixel beside the
    water), out to about PLUME_KM_PER_ROOT km per (m3/s)^0.5 of its discharge, Gaussian, larger mouths in classes a
    decade wide."""
    beside = np.roll(water, 1, axis=1) | np.roll(water, -1, axis=1)
    beside[1:] |= water[:-1]
    beside[:-1] |= water[1:]
    mouth = np.where((discharge > 0) & ~water & beside, discharge, 0.0)
    plume = np.zeros(water.shape, np.float32)
    for low in (10.0, 100.0, 1000.0, 10000.0):
        high = low * 10.0 if low < 10000.0 else np.inf
        these = ((mouth >= low) & (mouth < high)).astype(np.float32)
        if not these.any():
            continue
        sigma = PLUME_KM_PER_ROOT * np.sqrt(low * np.sqrt(10.0)) / km
        spread = sphere_blur(these, sigma, lat) * (2.0 * np.pi * sigma * sigma)
        plume = np.maximum(plume, np.clip(spread, 0.0, 1.0).astype(np.float32))
    return plume * water


def falls_grids(shape, record, product, rotation=None):
    """Foam and mist on the texture grid (180 W at the left edge), in the globe's frame, from the falls study: each
    steep step's foam share by its steepness, carried as channels are, and the mist round each listed reach of MIST_GW
    or more, Gaussian, its width growing with the cube root of the reach's power, peaking at MIST_MAX."""
    from scipy.ndimage import maximum_filter
    rows, cols = shape
    with np.load(product, allow_pickle=False) as z:
        cell, slope, reach = z['steep_cell'], z['steep_slope'], z['steep_reach']
    foam = np.zeros(rows * cols, np.float32)
    foam[cell] = smoothstep(FOAM_SLOPE[0], FOAM_SLOPE[1], slope)
    roll = cols // 2
    foam = np.roll(foam.reshape(shape), roll, axis=1)
    foam = turn_channels(foam, rotation) if rotation is not None else foam
    ppd = cols / 360.0
    lat = 90.0 - (np.arange(rows) + 0.5) / ppd
    km = MOON_KM * np.radians(1.0 / ppd)
    mist = np.zeros(shape, np.float32)
    powers = [x['power_gw'] for x in record['reaches']]
    for lo, hi in ((MIST_GW, 6.0), (6.0, 18.0), (18.0, np.inf)):
        these = [i for i, pw in enumerate(powers) if lo <= pw < hi]
        if not these:
            continue
        src = np.zeros(rows * cols, np.float32)
        src[cell[np.isin(reach, these)]] = 1.0
        src = np.roll(src.reshape(shape), roll, axis=1)
        if rotation is not None:
            src = turn_channels(src, rotation)
        sigma = MIST_KM * np.cbrt(np.sqrt(lo * (hi if np.isfinite(hi) else 2.0 * lo))) / km
        spread = sphere_blur(src, sigma, lat) * (2.0 * np.pi * sigma * sigma)
        mist = np.maximum(mist, MIST_MAX * np.clip(spread, 0.0, 1.0).astype(np.float32))
    return dict(foam=foam, mist=mist)


def true_colour(height, level, ppd, drainage, moisture, rock, waters, falls=None):
    """The ground in true colour (linear sRGB, float32), the water mask and the share of each pixel river channels
    cover, 180 W at the left edge, on a surface holding every water body at its own level (surface_grids). The cover
    follows the moisture map's aridity, wetter along larger rivers and in hollows and drier on steep slopes, with noise
    for patchiness (ground_colour), over the soil weathered from the rock (rock: FeO and TiO2 in weight per cent and
    the basalt's share, on the texture grid; soil_colour); the waters take the designed waters' colours: open sea, the
    productive coast's over shelves and in shallow lakes, and the river mouth's in rivers and where rivers carry
    sediment into a lake or sea (drainage['sediment'], the runoff gathered since the last lake or sea)."""
    height = height.astype(np.float32)
    rows, _ = height.shape
    lat = 90.0 - (np.arange(rows) + 0.5) / ppd
    km = MOON_KM * np.radians(1.0 / ppd)
    water = np.isfinite(drainage['lake_level'])
    surface_level = np.where(water, drainage['lake_level'], level).astype(np.float32)
    discharge = np.where(water, 0.0, drainage['discharge']).astype(np.float32)
    river = np.clip(RIVER_WIDTH * np.sqrt(discharge) / (km * 1000.0), 0.0, 1.0).astype(np.float32)
    ln_aridity = moisture_field(moisture, height.shape)
    ground = np.where(water, surface_level, height)
    curvature = sphere_laplacian(sphere_blur(ground, 3.0, lat), lat).astype(np.float32)
    ln_aridity += HOLLOW_LN * np.clip(curvature / (curvature[~water].std() + 1e-6), -2.5, 2.5)
    del curvature
    ln_aridity += FLOODPLAIN_LN * np.minimum(sphere_blur(np.log10(1.0 + discharge), 2.5, lat).astype(np.float32), 2.0)
    gy, gx = gradients(ground, km * 1000.0, lat)
    slope = np.hypot(gx, gy)
    del gy, gx, ground
    ln_aridity -= SLOPE_LN * smoothstep(0.2, 0.6, slope)
    ln_aridity += PATCH_LN * (sphere_fbm(height.shape, 11, octaves=9, base=24, spread=0.146) - 0.5) / 0.146
    aridity = np.exp(ln_aridity)
    del ln_aridity
    colour = ground_colour(aridity, soil_colour(aridity, rock['feo'], rock['tio2']), rock['mafic'],
                           0.6 * smoothstep(0.35, 0.7, slope))
    del aridity, slope
    shelf = smoothstep(SHELF_M[0], SHELF_M[1], np.where(water, surface_level - height, 0.0))[..., None]
    sediment = np.where(water, 0.0, drainage['sediment'] if 'sediment' in drainage else discharge)
    plume = plume_share(sediment, water, lat, km)[..., None]
    del sediment
    sea = (waters['open_sea'] * (1 - shelf) + waters['productive_coast'] * shelf) * (1 - plume) + waters['river_mouth'] * plume
    del shelf, plume
    river_colour = waters['river_mouth']
    if falls is not None:                            # white water where the rivers drop steeply
        river_colour = river_colour * (1 - falls['foam'][..., None]) + WHITE_WATER * falls['foam'][..., None]
    colour = colour * (1 - river[..., None]) + river_colour * river[..., None]
    colour = np.where(water[..., None], sea, colour).astype(np.float32)
    if falls is not None:                            # and spray round the largest falls
        colour = colour * (1 - falls['mist'][..., None]) + MIST_WHITE * falls['mist'][..., None]
    return colour.astype(np.float32), water, river, surface_level


LEGEND_ARIDITY = ((2.0, 'Wet forest'), (0.8, 'Humid forest'), (0.57, 'Dry woodland'), (0.33, 'Semi-arid grassland'),
                  (0.1, 'Arid shrubland'), (0.02, 'Desert, highland soil'))
HIGHLAND_ROCK = (SOIL_FEO[0], 0.3)   # FeO and TiO2 (weight per cent) of the legend's highland and basalt
MARE_ROCK = (SOIL_FEO[1], 1.5)


def true_legend(waters):
    """The true-colour legend: the ground on highland rock at the moisture map's class middles, basalt desert, and the
    three waters."""
    ground = lambda a, rock, mafic: ground_colour(a, soil_colour(a, *rock), mafic)
    return ([[name, ground(a, HIGHLAND_ROCK, 0.0)] for a, name in LEGEND_ARIDITY]
            + [['Desert, basalt soil', ground(0.02, MARE_ROCK, 1.0)]]
            + [['Open sea', waters['open_sea']], ['Shelf seas, shallow lakes', waters['productive_coast']],
               ['Rivers, sediment plumes', waters['river_mouth']]])


def enlarge(f, shape):
    """A coarse field (rows from 90 N, a whole number of texture pixels per cell) enlarged to the texture grid,
    bilinear, wrapping in longitude."""
    import cv2
    s = shape[1] // f.shape[1]
    pad = np.concatenate([f[:, -1:], f, f[:, :1]], axis=1).astype(np.float32)
    return cv2.resize(pad, (pad.shape[1] * s, shape[0]), interpolation=cv2.INTER_LINEAR)[:, s:s + shape[1]]


def rock_fields(shape, ppd, rotation=None):
    """FeO and TiO2 (weight per cent) and the basalt's share on the texture grid (180 W at the left edge), in the
    globe's frame: the maria where LROC maps them, with Lunar Prospector's oxides there but at least 15% FeO (mare
    basalt soils hold 14-22%; the instrument's 5-degree footprint dilutes the small maria), and elsewhere the oxides
    the highlands hold once the maria are left out of the smoothing."""
    from geography import mare_map, surface_composition as composition
    mare = mare_map.basalt(ppd)
    coarse = block_mean(mare.astype(np.float32), ppd) > 0.5
    whole, highland = composition.grid(1), composition.grid(1, exclude=coarse)
    feo = np.where(mare, np.maximum(enlarge(whole['feo'], shape), 15.0), enlarge(highland['feo'], shape))
    tio2 = np.where(mare, enlarge(whole['tio2'], shape), enlarge(highland['tio2'], shape))
    fields = dict(feo=feo.astype(np.float32), tio2=np.maximum(tio2, 0.0).astype(np.float32), mafic=mare.astype(np.float32))
    if rotation is not None:
        fields = {k: turn_texture(v, rotation) for k, v in fields.items()}
    return fields


def relief_normals(height, ppd):
    """Tangent-space normals of the ground, sea floors and lake beds included, relief baked at NORMAL_STRENGTH.
    The globe's appearance mode flattens them under water, and its map mode under the sea."""
    surf = height.astype(np.float32)
    lat = 90.0 - (np.arange(surf.shape[0]) + 0.5) / ppd
    dy = MOON_KM * 1000.0 * np.radians(1.0 / ppd)
    gy, gx = gradients(surf, dy, lat)
    n = np.stack([-gx * NORMAL_STRENGTH, gy * NORMAL_STRENGTH, np.ones_like(surf)], axis=-1)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    return np.round((n + 1) / 2 * 255).astype(np.uint8)


SKY_AZIMUTHS = 16                    # directions searched for the horizon
SKY_REACH_KM = 50.0                  # how far along each
SKY_STEPS = 20                       # distances along each, spaced evenly in their logarithm from half a pixel


def sky_view(height, ppd, water=None):
    """The share of an even sky's diffuse light that reaches level ground (float32, 1 in the open): the mean over
    SKY_AZIMUTHS directions of cos^2 of the horizon's elevation, the horizon found out to SKY_REACH_KM along each with
    the Moon's curvature taken off (the ground falls away by d^2 / 2R, 0.7 km at 50 km). Rows from 90 N, columns
    wrapping in longitude; water, where given, sees the whole sky."""
    import cv2
    h = height.astype(np.float32)
    rows, cols = h.shape
    km = MOON_KM * np.radians(1.0 / ppd)
    lat = 90.0 - (np.arange(rows) + 0.5) / ppd
    per_col = (1.0 / (km * np.maximum(np.cos(np.radians(lat)), 0.02))).astype(np.float32)[:, None]
    yy, xx = np.mgrid[0:rows, 0:cols].astype(np.float32)
    seen = np.zeros_like(h)
    for k in range(SKY_AZIMUTHS):
        north, east = np.cos(2.0 * np.pi * k / SKY_AZIMUTHS), np.sin(2.0 * np.pi * k / SKY_AZIMUTHS)
        steepest = np.zeros_like(h)
        for d in np.geomspace(0.5 * km, SKY_REACH_KM, SKY_STEPS):
            map_x = np.mod(xx + np.float32(d * east) * per_col, np.float32(cols))
            map_y = np.clip(yy - np.float32(d * north / km), 0.0, rows - 1.0)
            ahead = cv2.remap(h, map_x, map_y, cv2.INTER_LINEAR, borderMode=cv2.BORDER_WRAP)
            np.maximum(steepest, ((ahead - h) / 1000.0 - d * d / (2.0 * MOON_KM)) / d, out=steepest)
        seen += 1.0 / (1.0 + steepest * steepest)
    seen /= SKY_AZIMUTHS
    if water is not None:
        seen[water] = 1.0
    return seen


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


def build(height, level, ppd, climatology=CLIMATOLOGY, drainage_path=DRAINAGE, rotation=None, surface_path=None,
          moisture_path=None, routing_path=None, falls_paths=None):
    """Appearance data for the globe page: textures (arrays, 180 W at the left edge) and parameters.

    height: metres above the geoid at ppd pixels per degree, 0 E at the left edge as LOLA has it, in the globe's frame;
    climatology: in the same frame; rotation: the turned frame's, which carries the drainage product into it.
    surface_path: a surface holding every water body at its own level (surface_grids), in place of the sea below
    `level` and the drainage product; `level` then stays the heights' reference for the ground's look.
    moisture_path: with a surface, the moisture map of the climatology's run and years, which sets the ground in true
    colour (true_colour) in place of the guessed cover types; routing_path, the routing product the surface was drawn
    from, whose receivers give the rivers' sediment-bearing flow."""
    moon, engine = engine_atmosphere()
    colour_rel = np.array(moon['sun_colour_linear_srgb'])
    sun_rgb_lux = moon['sun_illuminance_lux'] * colour_rel / float(LUMINANCE @ colour_rel)
    clim = dict(np.load(climatology))
    clim_meta = json.loads(str(clim.pop('metadata')))
    clim['deck_sigma'] = clim_meta['deck_sigma']
    h = np.roll(height, height.shape[1] // 2, axis=1)
    truly = surface_path is not None and moisture_path is not None
    if truly:
        from geography import fetch_surface_inputs as geography_inputs, mare_map, surface_composition
        moisture = moisture_record(moisture_path)
        made = moisture['producer']
        if made['run'] != clim_meta['run'] or [made['years'][0], made['years'][-1]] != [clim_meta['years'][0], clim_meta['years'][-1]]:
            raise SystemExit(f"The moisture map is of {made['run']} years {made['years']}, the clouds of {clim_meta['run']} "
                             f"years {clim_meta['years']}: build it on the same run and years")
        drainage = surface_grids(h.shape, surface_path, rotation, routing=routing_path)
        rock = rock_fields(h.shape, ppd, rotation)
        waters = water_albedos()
        falls = None
        if falls_paths is not None:
            falls_record = json.loads(Path(falls_paths[0]).read_text())
            falls = falls_grids(h.shape, falls_record, falls_paths[1], rotation)
        colour, water, river, _ = true_colour(h, level, ppd, drainage, moisture, rock, waters, falls)
        del falls
        weight = np.cos(np.radians(90.0 - (np.arange(h.shape[0]) + 0.5) / ppd))[:, None]
        land = (~water) * weight
        basalt_land = float(((rock['mafic'] > 0.5) * land).sum() / land.sum())
        feo_land = float((rock['feo'] * land).sum() / land.sum())
        del rock, weight, land
    elif surface_path is not None:
        drainage = surface_grids(h.shape, surface_path, rotation)
        colour, water, river, _ = surface_colour(h, level, clim, ppd, drainage, sea=np.zeros(h.shape, bool))
    else:
        drainage = drainage_grids(h.shape, drainage_path, rotation)
        colour, water, river, _ = surface_colour(h, level, clim, ppd, drainage)
    albedo = np.empty(colour.shape, np.uint8)
    for c in range(3):
        albedo[..., c] = np.round(srgb_encode(colour[..., c]) * 255)
    del colour
    if truly:                        # the full relief, and the sky each place sees
        normal = relief_normals(h, ppd)
        sky = np.round(sky_view(np.where(water, drainage['lake_level'], h), ppd, water) * 255).astype(np.uint8)
    else:
        normal = relief_normals(block_mean(h, 2), ppd // 2)
    wet_share = water.astype(np.float32) if river is None else np.maximum(water.astype(np.float32), river)
    water_mask = np.round(wet_share * 255).astype(np.uint8)
    del wet_share

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
    if truly:
        textures['sky'] = sky
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
    sources = dict(engine=ENGINE, sky_atlas=SKY_ATLAS, climatology=Path(climatology))
    if truly:
        params['legend'] = true_legend(waters)
        params['legendTitle'] = 'Ground in true colour'
        params['legendNote'] = ('Each in daylight, from Earth\'s land and soils under the same dryness: aridity (rain '
                                'over the water the air could draw off) ' + ', '.join(f'{a:g}' for a, _ in LEGEND_ARIDITY)
                                + ' on highland rock, and basalt desert. From space the air veils them in blue.')
        params['cover'] = dict(moisture=dict(schema=moisture['schema'], run=made['run'], years=made['years']),
                               cover=COVER, trees=TREES, green=GREEN, basaltLandShare=round(basalt_land, 4),
                               landFeO=round(feo_land, 2),
                               soils={end: [[a, m, name] for a, m, name, _ in knots] for end, knots in SOIL_MUNSELL.items()},
                               waters={k: [round(float(v), 5) for v in a] for k, a in waters.items()})
        sources.update(moisture=Path(moisture_path), mare=geography_inputs.path(mare_map.ARCHIVE),
                       composition=geography_inputs.path(surface_composition.TABLE), water_colours=WATER_COLOURS)
        if routing_path is not None:
            sources['routing'] = Path(routing_path)
        if falls_paths is not None:
            sources.update(falls=Path(falls_paths[0]), falls_grid=Path(falls_paths[1]))
            params['cover']['falls'] = dict(foamSlope=list(FOAM_SLOPE), mistGW=MIST_GW, mistKm=MIST_KM, mistMax=MIST_MAX)
    params['legend'] = [[name, '#' + ''.join(f'{int(round(v * 255)):02x}' for v in srgb_encode(np.asarray(rgb)))]
                        for name, rgb in params['legend']]
    if surface_path is not None:
        sources['surface'] = Path(surface_path)
        params['drainage'] = dict(schema=drainage['metadata']['schema'], riverWidthPerRootDischarge=RIVER_WIDTH)
        return textures, params, sources, dict(sea=np.isfinite(drainage['lake_level']), water=water)
    if drainage is not None:
        sources['drainage'] = Path(drainage_path)
        params['drainage'] = dict(schema=drainage['metadata']['schema'], riverWidthPerRootDischarge=RIVER_WIDTH)
    return textures, params, sources, dict(sea=h < level, water=water)
