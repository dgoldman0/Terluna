"""Clear molecular paths to elevated features in spherical, layered air.

This is a direct-beam calculation, not cloud radiance or a scattering solver.
Layer intersections are analytic; existing domain products supply the optics.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np

from atmosphere.radiative_convective.optics import rayleigh_optical_depth
from illumination.sky.atmospheres import MOON_ZERO, O3, SF, SW, rayleigh
from illumination.sky.colour_matching import cmf
from illumination.surface_light.model import COLUMNS, build, column_facts, sunlight
from shared.constants import AU, MOON_RADIUS, SUN_RADIUS

SUN_RADIUS_DEG = math.degrees(math.asin(SUN_RADIUS / AU))


@dataclass
class MolecularColumn:
    radius_m: float
    height_m: np.ndarray
    extinction_m1: np.ndarray
    wavelength_nm: np.ndarray

    def __post_init__(self):
        self.height_m = np.asarray(self.height_m, dtype=float)
        self.extinction_m1 = np.asarray(self.extinction_m1, dtype=float)
        self.wavelength_nm = np.asarray(self.wavelength_nm, dtype=float)
        if (not np.isfinite(self.radius_m) or self.radius_m <= 0 or
                self.height_m.ndim != 1 or len(self.height_m) < 2 or self.height_m[0] != 0 or
                not np.isfinite(self.height_m).all() or np.any(np.diff(self.height_m) <= 0) or
                self.wavelength_nm.ndim != 1 or len(self.wavelength_nm) == 0 or
                not np.isfinite(self.wavelength_nm).all() or np.any(self.wavelength_nm <= 0) or
                np.any(np.diff(self.wavelength_nm) <= 0) or
                self.extinction_m1.shape != (len(self.height_m) - 1, len(self.wavelength_nm)) or
                not np.isfinite(self.extinction_m1).all() or np.any(self.extinction_m1 < 0)):
            raise ValueError('Invalid spherical molecular column')

    def path(self, height_m, mu, length_m=math.inf):
        """Lengths in shells and solid-body blockage, for a ray or finite segment.

        mu is the outward radial direction cosine at the starting point. A
        finite ray may point below horizontal without intersecting the ground.
        Tangency alone has zero solid-body path and is not blockage.
        """
        if (not np.isfinite([height_m, mu]).all() or not 0 <= height_m <= self.height_m[-1] or
                not -1 <= mu <= 1 or math.isnan(length_m) or length_m < 0):
            raise ValueError('Invalid optical path')
        r = self.radius_m + height_m
        b2, centre = r * r * (1 - mu * mu), -r * mu
        edges = self.radius_m + self.height_m
        # Length inside each nested sphere, clipped to the requested segment.
        discriminant = edges * edges - b2
        roots = np.sqrt(np.maximum(discriminant, 0))
        inside = np.maximum(0, np.minimum(centre + roots, length_m) - np.maximum(centre - roots, 0))
        inside[discriminant <= 0] = 0
        if mu >= 0:
            # An outward ray starting outside the solid body cannot enter it;
            # near the horizon, subtracting large squares can leave roundoff.
            inside[0] = 0
        return np.maximum(np.diff(inside), 0), bool(inside[0] > 1e-7)

    def transmission(self, height_m, mu, length_m=math.inf):
        lengths, blocked = self.path(height_m, mu, length_m)
        if blocked:
            return np.zeros(len(self.wavelength_nm))
        return np.exp(-lengths @ self.extinction_m1)

    def sun_transmission(self, height_m, sun_elevation_deg, order=32):
        """Uniform finite solar disk, small-angle vertical-strip quadrature.

        The neglected transverse disk curvature and centre-normal projection
        are O(solar angular radius squared); refraction/limb darkening omitted.
        """
        if not np.isfinite(sun_elevation_deg) or not -90 <= sun_elevation_deg <= 90:
            raise ValueError('Invalid solar elevation')
        offsets, weights = solar_strips(order)
        values = [self.transmission(height_m, math.sin(math.radians(sun_elevation_deg + d)))
                  for d in offsets]
        return weights @ np.asarray(values)


@lru_cache(maxsize=8)
def solar_strips(order):
    if not isinstance(order, int) or order < 2:
        raise ValueError('Solar disk order must be an integer >= 2')
    theta = np.arange(1, order + 1) * math.pi / (order + 1)
    return SUN_RADIUS_DEG * np.cos(theta), 2 * np.sin(theta) ** 2 / (order + 1)


def solved_moon():
    """Current moist/solved column with its layer-mean Rayleigh and ozone.

    Rayleigh layer columns are preserved exactly. Ozone uses the existing
    coarse visible cross sections. Other absorption bands are omitted.
    """
    col, air, trace, hashes = build(COLUMNS['moon_1.2atm'])
    beta = rayleigh_optical_depth(col, 1e7 / SW) / np.diff(col['z_m'])[:, None]
    beta += trace['O3'][:, None] * 1e6 * O3[None, :]
    return MolecularColumn(MOON_RADIUS, col['z_m'], beta, SW), dict(
        column_facts=column_facts(col, air, trace), stored_profiles=hashes,
        density='Piecewise constant layer extinction; vertical molecular columns preserved exactly.',
        omitted='H2O/O2 and other absorption, aerosols, cloud extinction, refraction, all in-scattering.')


def atlas_proxy(step_m=1000.):
    """The zero-ozone atlas exponential atmosphere, in conservative thin shells."""
    if not np.isfinite(step_m) or step_m <= 0:
        raise ValueError('Positive layer spacing required')
    atm = MOON_ZERO
    z = np.linspace(0, atm.top - atm.radius_m, math.ceil((atm.top - atm.radius_m) / step_m) + 1)
    average = atm.scale_height_m * (np.exp(-z[:-1] / atm.scale_height_m)
                                    - np.exp(-z[1:] / atm.scale_height_m)) / np.diff(z)
    return MolecularColumn(MOON_RADIUS, z, average[:, None] * rayleigh(SW, atm), SW)


def photopic_weights(design=False):
    """Existing 10-nm solar spectrum and CMF fit, trapezoidal lux weights."""
    width = np.gradient(SW)
    width[[0, -1]] *= .5
    transmission, _ = sunlight('design' if design else 'unfiltered')
    source = SF * (transmission(SW) if transmission is not None else 1)
    return 683 * width * cmf(SW)[:, 1] * source


def cloud_view(radius_m, height_m, surface_range_m, observer_height_m=1.6):
    """Equatorial great-circle section; positive range is west, toward sunset."""
    if (not np.isfinite([radius_m, height_m, surface_range_m, observer_height_m]).all() or
            radius_m <= 0 or height_m < 0 or observer_height_m < 0 or
            abs(surface_range_m) >= math.pi * radius_m):
        raise ValueError('Invalid cloud-view geometry')
    alpha = surface_range_m / radius_m
    x = (radius_m + height_m) * math.sin(alpha)
    z = (radius_m + height_m) * math.cos(alpha) - (radius_m + observer_height_m)
    distance = math.hypot(x, z)
    if distance == 0:
        raise ValueError('Feature and observer coincide')
    return dict(arc_angle_deg=math.degrees(alpha), ray_length_m=distance, mu=z / distance,
                apparent_elevation_deg=math.degrees(math.atan2(z, abs(x))),
                ray_x=x / distance, ray_z=z / distance)


def sunset_depression_deg(radius_m, height_m, surface_range_m=0.):
    """Observer's Sun depression when the cloud's solar centre meets its horizon."""
    if (not np.isfinite([radius_m, height_m, surface_range_m]).all() or
            radius_m <= 0 or height_m < 0 or abs(surface_range_m) >= math.pi * radius_m):
        raise ValueError('Invalid sunset geometry')
    return math.degrees(surface_range_m / radius_m + math.acos(radius_m / (radius_m + height_m)))
