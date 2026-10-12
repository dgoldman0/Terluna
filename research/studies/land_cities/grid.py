"""The 28% atlas on its 4 px/deg grid, as the lava-tube and highland screens read it.

Dry land is atlas land less the drainage estimate's rain-fed lakes, counted as the share of each atlas cell that
the 16 px/deg lake cells cover, as the ways-of-living study counts it. Sides follow that study's rule: nearside where
the Earth stands above the horizon beyond the optical libration's reach, far side where it stays below it, limb
between. Fields on the GCM's T21 grid are carried onto the atlas bilinearly, periodic in longitude.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from shared.constants import LIBRATION_LATITUDE_DEG, LIBRATION_LONGITUDE_DEG, MOON_RADIUS

RADIUS_KM = MOON_RADIUS / 1e3
LAKE_PPD = 16                 # the drainage grid's pixels per degree
WET_CELL_SHARE = 0.5          # an atlas cell at least half under rain-fed lakes counts as standing water


def pixel_areas(lat_deg, ppd, ncols):
    """Exact spherical cell areas (km2) of a cylindrical grid with pixel-centred latitudes."""
    top = np.radians(np.asarray(lat_deg) + 0.5 / ppd)
    bottom = np.radians(np.asarray(lat_deg) - 0.5 / ppd)
    row = RADIUS_KM ** 2 * math.radians(1.0 / ppd) * (np.sin(top) - np.sin(bottom))
    return np.repeat(row[:, None], ncols, axis=1)


def great_circle_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    c = np.clip(np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(np.radians(np.asarray(lon2) - lon1)),
                -1.0, 1.0)
    return RADIUS_KM * np.arccos(c)


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def chord_km(chord):
    return 2.0 * RADIUS_KM * np.arcsin(np.minimum(np.asarray(chord) / 2.0, 1.0))


def regrid(field, src_lat, src_lon, lat, lon):
    """Bilinear interpolation of a [lat][lon] field onto a grid of latitudes and longitudes (east, degrees),
    periodic in longitude and held at the outermost rows beyond them."""
    f = np.asarray(field, float)
    src_lat = np.asarray(src_lat, float)
    order = np.argsort(src_lat)
    ls, fs = src_lat[order], f[order]
    src_lon = np.asarray(src_lon, float)
    lon_src = np.concatenate([src_lon - 360.0, src_lon, src_lon + 360.0])
    fs = np.concatenate([fs, fs, fs], axis=1)
    la = np.clip(np.asarray(lat, float), ls[0], ls[-1])
    i = np.clip(np.searchsorted(ls, la) - 1, 0, ls.size - 2)
    ty = ((la - ls[i]) / (ls[i + 1] - ls[i]))[:, None]
    lo = np.asarray(lon, float) % 360.0
    j = np.clip(np.searchsorted(lon_src, lo) - 1, 0, lon_src.size - 2)
    tx = ((lo - lon_src[j]) / (lon_src[j + 1] - lon_src[j]))[None, :]
    top = fs[i][:, j] * (1 - tx) + fs[i][:, j + 1] * tx
    bottom = fs[i + 1][:, j] * (1 - tx) + fs[i + 1][:, j + 1] * tx
    return top * (1 - ty) + bottom * ty


def sides(LA, LO):
    """0 nearside, 1 limb, 2 far side, by where the Earth stands against the libration's reach."""
    cos_t = np.cos(np.radians(LA)) * np.cos(np.radians(LO))
    earth_elev = 90.0 - np.degrees(np.arccos(np.clip(cos_t, -1.0, 1.0)))
    lon_c = np.radians(((LO + 180.0) % 360.0) - 180.0)
    az = np.arctan2(np.sin(lon_c) * np.cos(np.radians(LA)), np.sin(np.radians(LA)))
    reach = np.hypot(LIBRATION_LONGITUDE_DEG * np.sin(az), LIBRATION_LATITUDE_DEG * np.cos(az))
    return np.where(earth_elev > reach, 0, np.where(earth_elev < -reach, 2, 1))


def lake_mask16(drainage):
    """The drainage estimate's rain-fed lake cells as a 16 px/deg mask, and each cell's lake level (m above the
    geoid; NaN off the lakes)."""
    rows, cols = 180 * LAKE_PPD, 360 * LAKE_PPD
    lake = np.zeros(rows * cols, bool)
    lake[drainage['lake_cell']] = True
    level = np.full(rows * cols, np.nan, np.float32)
    level[drainage['lake_cell']] = drainage['lake_level_m']
    return lake.reshape(rows, cols), level.reshape(rows, cols)


@dataclass
class Atlas:
    lat: np.ndarray
    lon: np.ndarray
    LA: np.ndarray
    LO: np.ndarray
    area: np.ndarray          # km2 per cell
    height_m: np.ndarray      # ground above the geoid
    sea_level_m: float
    hsl_km: np.ndarray        # ground above the sea
    water: np.ndarray         # seas and atlas lakes (every closed depression below the common level)
    water_label: np.ndarray
    f_lake: np.ndarray        # share of a land cell under rain-fed lakes
    lake_level_m: np.ndarray  # mean surface of the rain-fed lakes in a cell, above the geoid (NaN without lakes)
    land: np.ndarray          # km2 of land per cell
    dry: np.ndarray           # km2 of dry land per cell
    side: np.ndarray
    ppd: int

    def wet(self):
        """Standing water: the seas and atlas lakes, and cells at least half under rain-fed lakes."""
        return self.water | (self.f_lake >= WET_CELL_SHARE)

    def water_level_m(self):
        """The water's surface above the geoid on wet cells: the common level on the seas and atlas lakes, the lakes'
        own level on rain-fed lake cells."""
        return np.where(self.water, self.sea_level_m, self.lake_level_m)

    def km_to_water(self):
        """Great-circle km from every land cell to the nearest standing water (0 on water), and the index of that
        wet cell in the flattened grid."""
        wet = self.wet()
        xyz = unit(self.LA, self.LO)
        chord, idx = cKDTree(xyz[wet]).query(xyz[~wet])
        out = np.zeros(self.LA.shape)
        out[~wet] = chord_km(chord)
        nearest = np.arange(wet.size).reshape(wet.shape)
        flat_wet = np.flatnonzero(wet.ravel())
        nearest[~wet] = flat_wet[idx]
        return out, nearest

    def cell(self, lat, lon):
        """The cell holding a point; a point on a cell edge falls in the cell to its south and east, as the 16 px/deg
        lookups take it."""
        north, west = self.lat[0] + 0.5 / self.ppd, self.lon[0] - 0.5 / self.ppd
        i = int(np.clip(np.floor((north - lat) * self.ppd), 0, self.lat.size - 1))
        j = int(np.floor((lon % 360.0 - west) * self.ppd)) % self.lon.size
        return i, j


def build_atlas(atlas_grid, atlas_meta, drainage):
    h = atlas_grid['height_m'].astype(float)
    lat, lon = atlas_grid['lat_deg'].astype(float), atlas_grid['lon_deg'].astype(float)
    ppd = int(round(1.0 / abs(lat[0] - lat[1])))
    area = pixel_areas(lat, ppd, lon.size)
    sea = float(atlas_meta['sea_level_m'])
    water_label = atlas_grid['water_label']
    water = water_label > 0
    lake16, level16 = lake_mask16(drainage)
    factor = lake16.shape[0] // h.shape[0]
    blocks = (h.shape[0], factor, h.shape[1], factor)
    f_lake = lake16.reshape(blocks).mean(axis=(1, 3))
    f_lake = np.where(water, 0.0, f_lake)
    with np.errstate(invalid='ignore'):
        level_sum = np.nansum(level16.reshape(blocks).astype(float), axis=(1, 3))
        level = np.where(f_lake > 0.0, level_sum / np.maximum(lake16.reshape(blocks).sum(axis=(1, 3)), 1), np.nan)
    land = np.where(water, 0.0, area)
    LA, LO = np.meshgrid(lat, lon, indexing='ij')
    return Atlas(lat=lat, lon=lon, LA=LA, LO=LO, area=area, height_m=h, sea_level_m=sea, hsl_km=(h - sea) / 1e3,
                 water=water, water_label=water_label, f_lake=f_lake, lake_level_m=level, land=land,
                 dry=land * (1.0 - f_lake), side=sides(LA, LO), ppd=ppd)
