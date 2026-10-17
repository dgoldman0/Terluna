"""The flyover's terrain, regenerated from geography's inputs and bound to the island study's product and its atlas.

Heights are LOLA's LDEM_16 radii at their native pixel centres less the degree-200 GRAIL GL0420A geoid, evaluated on
those centres as the island study evaluates it, less the common water level the island study's product gives. A
global surround at 4 px/deg comes the same way from LDEM_4. The island is the connected land holding the product's
summit. The arrays build in about a second, so nothing is cached and no array is kept in Git.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import label

from geography import fetch_inputs
from geography.topography import Grid, geoid, read_ldem
from water_depth import SOUNDINGS

ROOT = Path(__file__).resolve().parents[2]
ISLAND = ROOT / 'research/studies/jura_island/results/shores.json'
ATLAS = ROOT / 'geography/results/atlas.json'
ISLAND_SCHEMA = 'terluna.research.jura-shores/1'
ATLAS_SCHEMA = 'terluna.geography.atlas/1'
INPUTS = ('ldem_16.img', 'ldem_4.img', 'jggrx_0420a_sha.tab')
STUDY_INPUTS = ('ldem_16.img', 'jggrx_0420a_sha.tab')        # the inputs the island study measured
GEOID_DEGREE = 200                                             # as the island study and the atlas take it
REGION = dict(south=27.0, north=64.0, west=293.0, east=352.0)  # the regional crop, degrees north and east


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    path = Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def products(island_path=None, atlas_path=None):
    """The island study's product and the atlas it was measured against, checked: their schemas, the atlas's bytes
    against the sha256 the island study recorded, and their two sea levels against each other."""
    island_path, atlas_path = Path(island_path or ISLAND), Path(atlas_path or ATLAS)
    island, atlas = json.loads(island_path.read_text()), json.loads(atlas_path.read_text())
    if island.get('schema') != ISLAND_SCHEMA:
        raise ValueError(f'Unsupported island product: {island.get("schema")}')
    if atlas.get('schema') != ATLAS_SCHEMA:
        raise ValueError(f'Unsupported atlas product: {atlas.get("schema")}')
    atlas_hash = digest(atlas_path)
    if atlas_hash != island['producer']['inputs'].get('geography/results/atlas.json'):
        raise ValueError('The atlas differs from the one the island study measured against')
    if island['sea_level_m'] != atlas['sea_level_m']:
        raise ValueError('The island study and the atlas give different sea levels')
    return dict(sea_level_m=float(island['sea_level_m']), water_share=float(atlas['water_share']),
                pixels_per_degree=int(island['pixels_per_degree']), island=island['island'],
                study_inputs=island['producer']['inputs'],
                hashes={relative(island_path): digest(island_path), relative(atlas_path): atlas_hash})


def check_inputs(study_inputs):
    """The sha256 of the LOLA and GRAIL files, checked against geography's manifest and, for the two the island study
    read, against the hashes its product recorded. Any difference stops the build."""
    manifest = {f['name']: f for f in fetch_inputs.manifest()['files']}
    hashes = {}
    for name in INPUTS:
        path = fetch_inputs.path(name)
        value = digest(path)
        if value != manifest[name]['sha256'] or path.stat().st_size != manifest[name]['bytes']:
            raise ValueError(f'Terrain source bytes differ from geography/inputs.json: {name}')
        key = f'geography/inputs/{name}'
        if name in STUDY_INPUTS and study_inputs.get(key) != value:
            raise ValueError(f'Terrain source bytes differ from those the island study measured: {name}')
        hashes[key] = value
    return hashes


def regional_grid(ppd):
    """Native pixel centres over the regional crop, north first, in degrees east."""
    lat = REGION['north'] - (np.arange(round((REGION['north'] - REGION['south']) * ppd)) + 0.5) / ppd
    lon = REGION['west'] + (np.arange(round((REGION['east'] - REGION['west']) * ppd)) + 0.5) / ppd
    return Grid(lat, lon, ppd)


def island_mask(lat, lon, height, summit_lat, summit_lon):
    """The 8-connected land holding the pixel nearest the summit (longitudes in one convention throughout)."""
    land = height > 0
    parts, _ = label(land, structure=np.ones((3, 3)))
    i, j = int(np.argmin(np.abs(lat - summit_lat))), int(np.argmin(np.abs(lon - summit_lon)))
    if not land[i, j]:
        raise ValueError('The island product\'s summit lies under water in the regenerated terrain')
    return parts == parts[i, j]


@dataclass
class Terrain:
    """Heights (m) above the common water level. Latitudes run north first; longitudes are degrees east, negative to
    the west, as the film's keys give them."""
    lat: np.ndarray
    lon: np.ndarray
    z: np.ndarray
    island: np.ndarray
    lat4: np.ndarray
    lon4: np.ndarray
    z4: np.ndarray
    sea_level_m: float
    water_share: float
    pixels_per_degree: int
    summit: tuple            # longitude, latitude and height (m) from the island study's product
    hashes: dict             # the products and inputs read, by sha256

    def heights(self, lon, lat):
        """Bilinear heights (m) at points inside the regional crop."""
        f = RegularGridInterpolator((self.lat[::-1], self.lon), self.z[::-1], bounds_error=True)
        return f(np.column_stack([np.ravel(lat), np.ravel(lon)])).reshape(np.shape(lon))

    def cell_area_km2(self):
        return Grid(self.lat, self.lon, self.pixels_per_degree).cell_area_m2 / 1e6

    def texture_bounds(self):
        """The regional crop's outer edges (west, east, south, north; degrees), half a pixel beyond its outermost
        pixel centres, so that each pixel centre falls on the centre of its texel in the depth texture."""
        half = 0.5 / self.pixels_per_degree
        return self.lon.min() - half, self.lon.max() + half, self.lat.min() - half, self.lat.max() + half

    def surround(self, step=2):
        """The global surround at every second 4 px/deg sample, with its first column repeated a turn to the east so
        that its mesh closes at 0 deg E: longitudes, latitudes and heights (m)."""
        z = self.z4[::step, ::step]
        return np.append(self.lon4[::step], self.lon4[0] + 360.0), self.lat4[::step], np.concatenate([z, z[:, :1]], 1)


def texture_uv(lon, lat, bounds):
    """The depth texture's coordinates of points: u from the west edge, v from the north edge."""
    west, east, south, north = bounds
    return np.stack([(lon - west) / (east - west), (north - lat) / (north - south)], axis=-1)


def build(island_path=None, atlas_path=None):
    p = products(island_path, atlas_path)
    ppd = p['pixels_per_degree']
    if ppd != 16:
        raise ValueError('The island product must stand at LOLA\'s native 16 px/deg')
    hashes = dict(p['hashes'], **check_inputs(p['study_inputs']))
    grid = regional_grid(ppd)
    rows = np.round((90 - grid.lat_deg) * ppd - 0.5).astype(int)
    cols = np.round(grid.lon_deg * ppd - 0.5).astype(int)
    raw = np.memmap(fetch_inputs.path('ldem_16.img'), dtype='<i2', mode='r', shape=(180 * ppd, 360 * ppd))
    z = np.asarray(raw[np.ix_(rows, cols)], dtype=float) * 0.5 - geoid(grid, max_degree=GEOID_DEGREE) - p['sea_level_m']
    radius4, grid4 = read_ldem(4)
    z4 = radius4 - geoid(grid4, max_degree=GEOID_DEGREE) - p['sea_level_m']
    s = p['island']
    summit = (s['summit_lon_deg'] - 360.0, s['summit_lat_deg'], float(s['summit_m']))
    lon = grid.lon_deg - 360.0
    return Terrain(lat=grid.lat_deg, lon=lon, z=z, island=island_mask(grid.lat_deg, lon, z, summit[1], summit[0]),
                   lat4=grid4.lat_deg, lon4=grid4.lon_deg - 360.0, z4=z4, sea_level_m=p['sea_level_m'],
                   water_share=p['water_share'], pixels_per_degree=ppd, summit=summit, hashes=hashes)


def soundings(terrain):
    """The named points' depths below the water (m), sampled bilinearly, with their values rounded to 10 m."""
    out = []
    for name, lon, lat in SOUNDINGS:
        depth = -float(terrain.heights(lon, lat))
        out.append(dict(name=name, longitude_deg=lon, latitude_deg=lat, depth_m=depth,
                        display_depth_m=int(np.round(depth / 10) * 10)))
    return out


def references(terrain):
    """Every reference point the panels name: longitude, latitude, the value shown (m) and its sense."""
    refs = {s['name']: (s['longitude_deg'], s['latitude_deg'], s['display_depth_m'], 'deep')
            for s in soundings(terrain)}
    lon, lat, height = terrain.summit
    refs['Jura summit'] = (lon, lat, int(round(height)), 'high')
    return refs
