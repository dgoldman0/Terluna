"""The food systems' zones on the 28% atlas: land by zone, its height and distance to water, and the seas by depth.

A grid screen on the ignored grid products (the atlas and drainage grids and the design climatology). Zones partition
the dry land, the land outside the seas and the rain-fed lakes:
- highlands: ground 4 km and more above the sea, at every latitude;
- rain belt: within 30 degrees of the equator, below 4 km;
- middle latitudes: 30-60 degrees below 4 km, split into wet land (the design run's rain at least 0.5 mm a day, the
  aerosol study's wet-land threshold) and the fog desert (less);
- polar land: poleward of 60 degrees below 4 km.
The rain-fed fill reads ecology's latitude bands, so the grid also gives each band's dry land, and the highland part
of it, at every height.
"""
from __future__ import annotations

import json
import math

import numpy as np
from scipy.spatial import cKDTree

from shared.constants import MOON_RADIUS

ATLAS = 'geography/products/atlas_28pct_4ppd.npz'
DRAINAGE = 'geography/products/drainage_28pct_16ppd.npz'
CLIMATOLOGY = 'climate/gcm/products/climatology_A28_dim5_moon.npz'
GRID_SCHEMAS = {ATLAS: 'terluna.geography.atlas-grid/1', DRAINAGE: 'terluna.geography.drainage/1',
                CLIMATOLOGY: 'terluna.climate.gcm-climatology/1'}
HIGHLAND_KM = 4.0                 # ways of living's highlands: ground 4 km and more above the sea
RAIN_BELT_DEG = 30.0              # ecology: 83% of the rain-fed land lies within 30 degrees
POLAR_DEG = 60.0                  # the aerosol study's polar dry land starts at 60 degrees
WET_MM_DAY = 0.5                  # the aerosol study's and the GCM README's wet land: at least 0.5 mm of rain a day
BIG_WATER_KM2 = 1000.0            # seas and great lakes (ways of living); smaller waters are lakes and wet margins
ECO_BANDS = ((0.0, 15.0), (15.0, 30.0), (30.0, 45.0), (45.0, 60.0), (60.0, 75.0), (75.0, 90.01))
DEPTHS_M = (10.0, 20.0, 50.0, 100.0, 200.0)
ZONES = ('highlands', 'rain_belt', 'middle_wet', 'fog_desert', 'polar')


class MissingInput(RuntimeError):
    """An ignored grid product or a bound input on another branch is not available here."""


def load_npz(root, rel):
    path = root / rel
    if not path.exists():
        raise MissingInput(f'{rel} is not present; regenerate it in its domain')
    with np.load(path, allow_pickle=False) as z:
        data = {k: z[k] for k in z.files}
    meta = json.loads(str(data.pop('metadata')))
    if meta.get('schema') != GRID_SCHEMAS[rel]:
        raise ValueError(f'Unexpected schema in {rel}: {meta.get("schema")}')
    return data, meta


def pixel_areas(lat_deg, ppd, ncols):
    """km2 of each cell of a regular latitude-longitude grid on the Moon's sphere."""
    top = np.radians(lat_deg + 0.5 / ppd)
    bottom = np.radians(lat_deg - 0.5 / ppd)
    row = (MOON_RADIUS / 1e3) ** 2 * math.radians(1.0 / ppd) * (np.sin(top) - np.sin(bottom))
    return np.repeat(row[:, None], ncols, axis=1)


def regrid(clim, name, lat, lon):
    """Bilinear interpolation of a T21 field onto the atlas grid (as ways of living and geography's drainage do)."""
    f = clim[name].astype(float)
    order = np.argsort(clim['lat'])
    ls, fs = clim['lat'][order].astype(float), f[order]
    lon_src = np.concatenate([clim['lon'] - 360.0, clim['lon'], clim['lon'] + 360.0]).astype(float)
    fs = np.concatenate([fs, fs, fs], axis=1)
    la = np.clip(lat, ls[0], ls[-1])
    i = np.clip(np.searchsorted(ls, la) - 1, 0, ls.size - 2)
    ty = ((la - ls[i]) / (ls[i + 1] - ls[i]))[:, None]
    j = np.clip(np.searchsorted(lon_src, lon) - 1, 0, lon_src.size - 2)
    tx = ((lon - lon_src[j]) / (lon_src[j + 1] - lon_src[j]))[None, :]
    top = fs[i][:, j] * (1 - tx) + fs[i][:, j + 1] * tx
    bottom = fs[i + 1][:, j] * (1 - tx) + fs[i + 1][:, j + 1] * tx
    return top * (1 - ty) + bottom * ty


def weighted_quantiles(values, weights, qs):
    order = np.argsort(values, kind='stable')
    v, w = values[order], weights[order]
    cw = np.cumsum(w) / w.sum()
    return [float(v[min(np.searchsorted(cw, q), v.size - 1)]) for q in qs]


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def zone_screen(root):
    """Dry land by zone and by ecology band, heights above the sea, distance to the seas and great lakes, rain, and
    the seas' area by depth. Raises MissingInput without the ignored grids."""
    atlas, meta = load_npz(root, ATLAS)
    drainage, _ = load_npz(root, DRAINAGE)
    clim, _ = load_npz(root, CLIMATOLOGY)
    h = atlas['height_m'].astype(float)
    label = atlas['water_label']
    lat, lon = atlas['lat_deg'].astype(float), atlas['lon_deg'].astype(float)
    ppd = int(round(1.0 / abs(lat[0] - lat[1])))
    area = pixel_areas(lat, ppd, lon.size)
    sea_level = float(meta['sea_level_m'])
    water = label > 0
    rows16, cols16 = 180 * 16, 360 * 16
    lake = np.zeros(rows16 * cols16, bool)
    lake[drainage['lake_cell']] = True
    factor = rows16 // h.shape[0]
    f_lake = lake.reshape(rows16, cols16).reshape(h.shape[0], factor, h.shape[1], factor).mean(axis=(1, 3))
    f_lake = np.where(water, 0.0, f_lake)
    dry = np.where(water, 0.0, area) * (1.0 - f_lake)
    LA, LO = np.meshgrid(lat, lon, indexing='ij')
    alat = np.abs(LA)
    hsl_km = (h - sea_level) / 1e3
    rain = regrid(clim, 'pr_mm_day', lat, lon)
    high = hsl_km >= HIGHLAND_KM
    masks = {
        'highlands': high & ~water,
        'rain_belt': ~high & ~water & (alat < RAIN_BELT_DEG),
        'middle_wet': ~high & ~water & (alat >= RAIN_BELT_DEG) & (alat < POLAR_DEG) & (rain >= WET_MM_DAY),
        'fog_desert': ~high & ~water & (alat >= RAIN_BELT_DEG) & (alat < POLAR_DEG) & (rain < WET_MM_DAY),
        'polar': ~high & ~water & (alat >= POLAR_DEG),
    }
    body_area = np.bincount(label.ravel(), weights=area.ravel())
    big = body_area >= BIG_WATER_KM2
    big[0] = False
    seas = big[label] & water
    xyz = unit(LA, LO)
    tree = cKDTree(xyz[seas])
    zones = {}
    for name in ZONES:
        m = masks[name]
        w = dry[m]
        chord, _ = tree.query(xyz[m])
        km = 2.0 * MOON_RADIUS / 1e3 * np.arcsin(np.minimum(chord / 2.0, 1.0))
        zones[name] = dict(
            dry_land_km2=float(w.sum()),
            height_above_sea_km_mean=float((hsl_km[m] * w).sum() / w.sum()),
            height_above_sea_km_p10_p50_p90=weighted_quantiles(hsl_km[m], w, (0.1, 0.5, 0.9)),
            km_to_sea_or_great_lake_p10_p50_p90=weighted_quantiles(km, w, (0.1, 0.5, 0.9)),
            km_to_sea_or_great_lake_mean=float((km * w).sum() / w.sum()),
            rain_mm_day=float((rain[m] * w).sum() / w.sum()),
            abs_lat_deg_mean=float((alat[m] * w).sum() / w.sum()))
    bands = []
    for lo, hi in ECO_BANDS:
        m = (alat >= lo) & (alat < hi) & ~water
        bands.append(dict(abs_lat_deg=[lo, min(hi, 90.0)], dry_land_km2=float(dry[m].sum()),
                          highland_dry_land_km2=float(dry[m & high].sum())))
    depth = sea_level - h
    sea_by_depth = {f'{d:g}': float(area[seas & (depth > 0.0) & (depth <= d)].sum()) for d in DEPTHS_M}
    return dict(moon_km2=float(area.sum()), dry_land_km2=float(dry.sum()), zones=zones, bands=bands,
                seas_and_great_lakes_km2=float(area[seas].sum()), sea_km2_by_depth_m=sea_by_depth,
                water_share=float(area[water].sum() / area.sum()), sea_level_m=sea_level,
                fog_desert_all_heights_km2=float(dry[~water & (alat >= RAIN_BELT_DEG) & (alat < POLAR_DEG) &
                                                     (rain < WET_MM_DAY)].sum()))
