"""Recover photopic angular shapes from the committed, packed clear-sky atlas.

The packer preserves luminance before half-float quantization. Its panorama
contains diffuse radiance; the solar beam is in separate tables. Only the upper
hemisphere is used here. Surface light comes from the current column product.
"""
from __future__ import annotations

import base64
from dataclasses import dataclass
import gzip
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .run import ROOT

ATLAS = ROOT / 'immersion/assets/sky/atmosphere.json'
Y_WEIGHTS = np.array([0.2126, 0.7152, 0.0722])
WORLD_ATLAS = {'moon': 'moon_no_ozone', 'earth': 'earth'}


def direction(elevation_deg, azimuth_deg):
    e, a = np.radians([elevation_deg, azimuth_deg])
    return np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])


@lru_cache(maxsize=4)
def sky_quadrature(order=128):
    """Gauss integration in elevation, periodic midpoint integration in azimuth."""
    x, w = np.polynomial.legendre.leggauss(order)
    e = (x + 1) * math.pi / 4
    az = (np.arange(4 * order) + 0.5) * 2 * math.pi / (4 * order)
    directions = np.stack(np.broadcast_arrays(
        np.cos(e)[:, None] * np.cos(az), np.cos(e)[:, None] * np.sin(az),
        np.sin(e)[:, None]), axis=-1)
    weights = np.broadcast_to((w * math.pi / 4 * np.cos(e))[:, None]
                              * 2 * math.pi / len(az), directions.shape[:-1])
    return directions.reshape(-1, 3), weights.ravel()


class AngularSky:
    """Bilinear luminance in sin(elevation) and symmetric solar-relative azimuth.

    The normalization is the exact horizontal integral of this interpolant.
    Returned radiance therefore has unit horizontal diffuse illuminance.
    """
    def __init__(self, mu, az, luminance):
        mu, az, luminance = map(np.asarray, (mu, az, luminance))
        if (luminance.shape != (len(mu), len(az)) or
                not np.all(np.isfinite(luminance)) or np.any(luminance < 0) or
                not np.all(np.diff(mu) > 0) or not np.all(np.diff(az) > 0) or
                not np.allclose([mu[0], mu[-1], az[0], az[-1]], [0, 1, 0, math.pi])):
            raise ValueError('Invalid sky grid')
        # Reflection symmetry doubles the integral over 0 <= azimuth <= pi.
        ring = 2 * np.sum((luminance[:, :-1] + luminance[:, 1:]) / 2 * np.diff(az), axis=1)
        a, b = mu[:-1], mu[1:]
        self.horizontal = float(np.sum((b - a) / 6 * ((2 * a + b) * ring[:-1]
                                                      + (a + 2 * b) * ring[1:])))
        if self.horizontal <= 0:
            raise ValueError('Positive sky flux required')
        self.maximum = float(luminance.max() / self.horizontal)
        self.interpolator = RegularGridInterpolator((mu, az), luminance / self.horizontal)

    def __call__(self, directions):
        d = np.asarray(directions)
        if np.any(d[..., 2] < -1e-10):
            raise ValueError('Sky lookup is only defined above the horizon')
        mu = np.clip(d[..., 2], 0, 1)
        az = np.abs(np.arctan2(d[..., 1], d[..., 0]))
        return self.interpolator(np.stack((mu, az), axis=-1))


def uniform_sky():
    return AngularSky(np.array([0., 1.]), np.array([0., math.pi]), np.ones((2, 2)))


class PackedAtlas:
    def __init__(self, path=ATLAS):
        self.data = json.loads(Path(path).read_text())
        if (self.data.get('version') != 1 or
                self.data.get('units') != 'photopic-weighted linear sRGB; sky is cd/m2-equivalent, irradiance lux-equivalent'):
            raise ValueError('Unsupported packed-atlas format or units')
        self.maps = {}
        for name in WORLD_ATLAS.values():
            w = self.data['worlds'][name]
            if (w['height'], w['width']) != (41, 33):
                raise ValueError('Unrecognized packed angular grid; check the producing bake')
            if not np.all(np.diff(w['suns']) > 0):
                raise ValueError('Solar elevations must increase')
            raw = gzip.decompress(base64.b64decode(w['data'], validate=True))
            rgb = np.frombuffer(raw, dtype='<f2').astype(float).reshape(len(w['suns']), 41, 33, 3)
            rgb *= np.asarray(w['scale'])[:, None, None, None]
            if np.any(rgb < 0) or not np.all(np.isfinite(rgb)):
                raise ValueError('Invalid packed radiance')
            self.maps[name] = rgb @ Y_WEIGHTS
        self.mu = np.sin(np.linspace(0, 1, 21) ** 2 * math.pi / 2)
        self.az = np.linspace(0, math.pi, 33)

    def get(self, world, sun_deg):
        name = WORLD_ATLAS[world]
        w = self.data['worlds'][name]
        if sun_deg not in w['suns']:
            raise ValueError('This study uses exact stored solar elevations')
        index = w['suns'].index(sun_deg)
        sky = AngularSky(self.mu, self.az, self.maps[name][index, 20:])
        return sky, dict(
            atlas_world=name, sun_deg=sun_deg,
            stored_diffuse_lux=float(np.asarray(w['diffuse'][index]) @ Y_WEIGHTS),
            reconstructed_diffuse_lux=sky.horizontal,
            reconstructed_over_stored=sky.horizontal / float(np.asarray(w['diffuse'][index]) @ Y_WEIGHTS),
            source_archive_sha256=next(p['sha256'] for p in self.data['provenance']
                                      if p['file'] == name + '_atlas.npz'),
        )


@dataclass
class LightField:
    sky: AngularSky
    sun_deg: float
    direct_h: float
    diffuse_h: float

    def __post_init__(self):
        if not 0 < self.sun_deg <= 90 or min(self.direct_h, self.diffuse_h) < 0:
            raise ValueError('Invalid daylight field')
        if not np.all(np.isfinite([self.sun_deg, self.direct_h, self.diffuse_h])):
            raise ValueError('Nonfinite daylight field')

    @property
    def sun(self):
        return direction(self.sun_deg, 0)

    @property
    def dni(self):
        return self.direct_h / math.sin(math.radians(self.sun_deg))

    def luminance(self, rays):
        return self.sky(rays) * self.diffuse_h

    def open_plane(self, normal, albedo=0.1, order=128):
        rays, weights = sky_quadrature(order)
        sky = float(np.sum(self.luminance(rays) * np.maximum(rays @ normal, 0) * weights))
        ground = albedo * (self.direct_h + self.diffuse_h) * (1 - normal[2]) / 2
        direct = self.dni * max(float(normal @ self.sun), 0)
        return dict(sky_lux=sky, ground_lux=ground, sun_direct_lux=direct,
                    ambient_lux=sky + ground, total_lux=sky + ground + direct)
