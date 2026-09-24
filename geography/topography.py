"""Lunar topography above an equipotential surface, from LOLA radii and the GRAIL geoid.

A water surface at rest follows an equipotential, not a height above a sphere.
The geoid here is the GRAIL GL0420A gravity field (degrees 2 to n_max, fully
normalized coefficients), plus Earth's static tidal potential scaled by 1 + k2
(the model's degree-2 terms exclude the permanent tide) and the Moon's
centrifugal potential. Heights are LOLA radius minus geoid radius, relative to
the 1737.4 km reference sphere; the geoid has zero mean over the sphere, so the
height datum is a constant offset that does not affect areas or volumes.

LOLA heights are in the mean-Earth/polar-axis frame and GL0420A in the DE421
principal-axis frame; the ~0.02 degree offset is ignored, below the ~50 km
resolution of a degree-200 geoid.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

from shared.constants import EARTH_GM, EARTH_MOON_DISTANCE, MOON_GM, MOON_RADIUS, SIDEREAL_MONTH_DAYS
from geography import fetch_inputs

LOVE_K2 = 0.0248          # GL0420A label: held fixed in the gravity solution


@dataclass(frozen=True)
class Grid:
    """Simple cylindrical grid with pixel-centred latitudes (north first) and east longitudes."""
    lat_deg: np.ndarray
    lon_deg: np.ndarray
    pixels_per_degree: int

    @property
    def cell_area_m2(self):
        """Exact area of each cell on the reference sphere (lat x lon)."""
        half = 0.5 / self.pixels_per_degree
        top = np.radians(self.lat_deg + half)
        bottom = np.radians(self.lat_deg - half)
        band = MOON_RADIUS**2 * np.radians(1.0 / self.pixels_per_degree) * (np.sin(top) - np.sin(bottom))
        return np.repeat(band[:, None], self.lon_deg.size, axis=1)


def read_ldem(pixels_per_degree: int):
    """LOLA LDEM radius grid as height (m) above the 1737.4 km sphere, and its grid."""
    name = {4: 'ldem_4.img', 16: 'ldem_16.img'}[pixels_per_degree]
    lines, samples = 180 * pixels_per_degree, 360 * pixels_per_degree
    raw = np.fromfile(fetch_inputs.path(name), dtype='<i2')
    if raw.size != lines * samples:
        raise ValueError(f'{name}: expected {lines}x{samples} samples')
    height = raw.reshape(lines, samples).astype(float) * 0.5
    lat = 90.0 - (np.arange(lines) + 0.5) / pixels_per_degree
    lon = (np.arange(samples) + 0.5) / pixels_per_degree
    return height, Grid(lat, lon, pixels_per_degree)


def read_gravity(max_degree: int):
    """GL0420A coefficients C[n, m], S[n, m] (fully normalized), GM (m^3/s^2), reference radius (m)."""
    rows = fetch_inputs.path('jggrx_0420a_sha.tab').read_text().splitlines()
    head = [v.strip() for v in rows[0].split(',')]
    ref_radius, gm = float(head[0]) * 1e3, float(head[1]) * 1e9
    if int(head[5]) != 1:
        raise ValueError('Expected fully normalized coefficients')
    c = np.zeros((max_degree + 1, max_degree + 1))
    s = np.zeros_like(c)
    for row in rows[1:]:
        v = row.split(',')
        n, m = int(v[0]), int(v[1])
        if n > max_degree:
            break
        c[n, m], s[n, m] = float(v[2]), float(v[3])
    return c, s, gm, ref_radius


def legendre_column(m: int, max_degree: int, t, u):
    """Fully normalized (4-pi, no Condon-Shortley phase) P_nm(t) for n = m..max_degree.

    t = sin(latitude), u = cos(latitude). Standard stable column recursion.
    """
    p = np.zeros((max_degree - m + 1, t.size))
    pmm = np.ones_like(t)
    for k in range(1, m + 1):
        pmm = pmm * u * math.sqrt((2 * k + 1) / (2 * k)) if k > 1 else pmm * u * math.sqrt(3.0)
    p[0] = pmm
    if max_degree > m:
        p[1] = math.sqrt(2 * m + 3) * t * pmm
    for n in range(m + 2, max_degree + 1):
        a = math.sqrt((2 * n - 1) * (2 * n + 1) / ((n - m) * (n + m)))
        b = math.sqrt((2 * n + 1) * (n + m - 1) * (n - m - 1) / ((n - m) * (n + m) * (2 * n - 3)))
        p[n - m] = a * t * p[n - m - 1] - b * p[n - m - 2]
    return p


def geoid(grid: Grid, max_degree: int = 200, radius: float = MOON_RADIUS, tides: bool = True):
    """Geoid height (m) above the reference sphere on a grid.

    Degree 0 and 1 are excluded (mean radius and centre of mass), so the gravity
    part has zero mean; the tidal and centrifugal terms are written as zero-mean
    degree-2 functions.
    """
    c, s, gm, ref = read_gravity(max_degree)
    lat, lon = np.radians(grid.lat_deg), np.radians(grid.lon_deg)
    t, u = np.sin(lat), np.cos(lat)
    scale = (ref / radius) ** np.arange(max_degree + 1)
    a = np.zeros((lat.size, max_degree + 1))
    b = np.zeros_like(a)
    for m in range(max_degree + 1):
        p = legendre_column(m, max_degree, t, u)
        n = np.arange(m, max_degree + 1)
        keep = n >= 2
        a[:, m] = (scale[n][keep] * c[n[keep], m]) @ p[keep]
        b[:, m] = (scale[n][keep] * s[n[keep], m]) @ p[keep]
    mm = np.arange(max_degree + 1)
    height = radius * (a @ np.cos(np.outer(mm, lon)) + b @ np.sin(np.outer(mm, lon)))
    if tides:
        gamma = gm / radius**2
        # Earth's static tide about the mean sub-Earth point (0 N, 0 E).
        cos_psi = np.outer(np.cos(lat), np.cos(lon))
        tide = (1 + LOVE_K2) * EARTH_GM * radius**2 / EARTH_MOON_DISTANCE**3 * 0.5 * (3 * cos_psi**2 - 1)
        omega = 2 * math.pi / (SIDEREAL_MONTH_DAYS * 86400.0)
        spin = 0.5 * omega**2 * radius**2 * (np.cos(lat)**2 - 2.0 / 3.0)
        height = height + (tide + spin[:, None]) / gamma
    return height


def height_above_geoid(pixels_per_degree: int = 16, max_degree: int = 200):
    """LOLA topography relative to the geoid (m), the geoid itself, and the grid."""
    radius_height, grid = read_ldem(pixels_per_degree)
    coarse = min(pixels_per_degree, 4)
    _, cgrid = read_ldem(coarse) if coarse != pixels_per_degree else (None, grid)
    n_coarse = geoid(cgrid, max_degree)
    if coarse != pixels_per_degree:
        # The degree-200 geoid is smooth on the 0.25-degree grid; interpolate it.
        from scipy.interpolate import RegularGridInterpolator
        lat_c = np.concatenate([[90.0], cgrid.lat_deg, [-90.0]])
        pad = np.vstack([n_coarse[:1].mean(1, keepdims=True).repeat(n_coarse.shape[1], 1), n_coarse,
                         n_coarse[-1:].mean(1, keepdims=True).repeat(n_coarse.shape[1], 1)])
        lon_c = np.concatenate([cgrid.lon_deg - 360.0, cgrid.lon_deg, cgrid.lon_deg + 360.0])
        pad = np.hstack([pad, pad, pad])
        f = RegularGridInterpolator((lat_c[::-1], lon_c), pad[::-1])
        la, lo = np.meshgrid(grid.lat_deg, grid.lon_deg, indexing='ij')
        n_full = f(np.stack([la.ravel(), lo.ravel()], 1)).reshape(la.shape)
    else:
        n_full = n_coarse
    return radius_height - n_full, n_full, grid
