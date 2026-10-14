"""Optional external gravity and eclipse geometry for reduced fleet replay.

All supplied positions share a Moon-centred, nonrotating coordinate frame.
The third-body force subtracts the same body's acceleration at the origin;
it is not an extra central attraction toward the external body's absolute
position. No ephemeris is opened, fetched or substituted by this module.
"""
from __future__ import annotations

from functools import lru_cache
from math import asin, isfinite

import numpy as np

from shared import constants as K


def _position(value, name):
    result = np.asarray(value, dtype=float)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise ValueError(f'{name} must be a finite three-vector')
    return result


def _positive_scalar(value, name):
    if np.ndim(value) != 0:
        raise ValueError(f'{name} must be a finite positive scalar')
    result = float(value)
    if not isfinite(result) or result <= 0:
        raise ValueError(f'{name} must be a finite positive scalar')
    return result


def differential_point_acceleration(points, perturbers):
    """Sum GM*((R-x)/|R-x|³-R/|R|³), retaining the shape of points.

    `points` has final dimension three. Each perturber is a dictionary with
    `position_m` and positive `gm_m3_s2`. The exact point-mass expression is
    evaluated with log1p/expm1 for far sources to avoid catastrophic loss from
    subtracting nearly equal origin and spacecraft accelerations. Coincidence
    with a source or a source at the coordinate origin is rejected.

    This is differential gravity in a frame whose origin shares each listed
    body's point-mass acceleration. Other perturbers, lunar harmonics and
    relativistic effects are not implicitly included.
    """
    original = np.asarray(points, dtype=float)
    if original.ndim == 0 or original.shape[-1] != 3 or not np.isfinite(original).all():
        raise ValueError('points must be finite with final dimension three')
    coordinates = original.reshape(-1, 3)
    result = np.zeros_like(coordinates)
    for source in perturbers:
        if not isinstance(source, dict) or 'position_m' not in source or 'gm_m3_s2' not in source:
            raise ValueError('Each perturber needs position_m and gm_m3_s2')
        position = _position(source['position_m'], 'Perturber position_m')
        gm = _positive_scalar(source['gm_m3_s2'], 'gm_m3_s2')
        radius = float(np.linalg.norm(position))
        if radius == 0:
            raise ValueError('External perturber cannot coincide with the coordinate origin')
        displacement = position-coordinates
        distance = np.linalg.norm(displacement, axis=1)
        if (distance == 0).any():
            raise ValueError('Point-mass gravity is singular at a perturber')
        normalized = coordinates/radius
        unit = position/radius
        delta = np.sum(normalized*normalized, axis=1)-2*(normalized@unit)
        far = np.abs(delta) < .5
        factor_minus_one = np.expm1(-1.5*np.log1p(delta[far]))
        result[far] += gm/radius**3*(position*factor_minus_one[:, None]-
                                    coordinates[far]*(1+factor_minus_one[:, None]))
        # Near a body the direct formula does not subtract near-equal fields.
        result[~far] += gm*(displacement[~far]/distance[~far, None]**3-position/radius**3)
    if not np.isfinite(result).all():
        raise ValueError('Gravity parameters overflowed the differential field')
    return result.reshape(original.shape)


def sphere_occlusion(centres, source_directions, spheres):
    """Opaque sphere transmission for forward rays, returned as (F,S) zeros/ones.

    Directions point from each face toward the source and must be unit length.
    Spheres contain `position_m` and positive `radius_m`. A ray crossing or
    tangentially touching a sphere at positive distance is blocked. Rays that
    leave a sphere surface immediately are unblocked; points strictly inside
    are blocked. The source is treated as distant, so callers must not supply
    occulting objects beyond a finite nearby light source.

    Disk integration, penumbra sampling and face-area sampling belong to the
    optical caller. There is no atmosphere, diffraction or light-time model.
    """
    faces = np.asarray(centres, dtype=float)
    directions = np.asarray(source_directions, dtype=float)
    if faces.ndim == 1:
        faces = faces[None, :]
    if directions.ndim == 1:
        directions = directions[None, :]
    if faces.ndim != 2 or faces.shape[1] != 3 or not np.isfinite(faces).all():
        raise ValueError('centres must be a finite (F,3) array')
    if (directions.ndim != 2 or directions.shape[1] != 3 or not len(directions) or
            not np.isfinite(directions).all() or
            not np.allclose(np.linalg.norm(directions, axis=1), 1., rtol=1e-10, atol=1e-12)):
        raise ValueError('source_directions must contain finite unit vectors')
    transmission = np.ones((len(faces), len(directions)))
    for sphere in spheres:
        if not isinstance(sphere, dict) or 'position_m' not in sphere or 'radius_m' not in sphere:
            raise ValueError('Each sphere needs position_m and radius_m')
        position = _position(sphere['position_m'], 'Sphere position_m')
        radius = _positive_scalar(sphere['radius_m'], 'radius_m')
        offset = position-faces
        longitudinal = offset@directions.T
        perpendicular = offset[:, None, :]-longitudinal[:, :, None]*directions[None, :, :]
        perpendicular_squared = np.sum(perpendicular*perpendicular, axis=2)
        hit = perpendicular_squared <= radius**2
        half_chord = np.sqrt(np.maximum(radius**2-perpendicular_squared, 0.))
        positive_intersection = longitudinal+half_chord > 0.
        transmission[hit & positive_intersection] = 0.
    return transmission


def lunar_environment(ephemeris, *, solar_disk_rings: int = 1):
    """Return cached at(t) callbacks using an already supplied ephemeris object.

    Accepts the existing Ephemeris.state API or InterpolatedEphemeris.at API.
    The payload directly supports Fleet.source and its perturbing_bodies and
    occluders callbacks. Source direction and angular radius are evaluated at
    the Moon origin, so solar parallax across an extended fleet is omitted.
    Earth and Sun supply differential point gravity; Earth and Moon are opaque
    spherical occultors. Other planets and nonspherical gravity remain absent.

    Missing files, out-of-range epochs and failed upstream interpolation are
    propagated to the caller. There is no circular-orbit or synthetic fallback.
    Cache entries must be treated as read-only by consumers.
    """
    if not callable(getattr(ephemeris, 'at', None)) and not callable(getattr(ephemeris, 'state', None)):
        raise ValueError('An already supplied state or interpolated ephemeris is required')
    if isinstance(solar_disk_rings, bool) or not isinstance(solar_disk_rings, int) or solar_disk_rings < 0:
        raise ValueError('solar_disk_rings must be a nonnegative integer')

    @lru_cache(maxsize=16)
    def at(t):
        if np.ndim(t) != 0 or not isfinite(float(t)):
            raise ValueError('Environment time must be a finite scalar')
        time = float(t)
        if callable(getattr(ephemeris, 'at', None)):
            positions = ephemeris.at(time)['positions']
            sun = _position(positions['sun'], 'Ephemeris Sun position')
            earth = _position(positions['earth'], 'Ephemeris Earth position')
        else:
            moon = _position(ephemeris.state('moon', time)[0], 'Ephemeris Moon position')
            sun = _position(ephemeris.state('sun', time)[0], 'Ephemeris Sun position')-moon
            earth = _position(ephemeris.state('earth', time)[0], 'Ephemeris Earth position')-moon
        distance = float(np.linalg.norm(sun))
        if distance <= K.SUN_RADIUS:
            raise ValueError('Moon-origin observer must be outside the Sun')
        if np.linalg.norm(earth) <= K.EARTH_RADIUS:
            raise ValueError('Earth centre cannot contain the Moon-origin observer')
        # Local import avoids a cycle when optics imports sphere_occlusion.
        from .flexible_optics import solar_disk_quadrature
        centre_direction = sun/distance
        directions, weights = solar_disk_quadrature(centre_direction, asin(K.SUN_RADIUS/distance), solar_disk_rings)
        # Normalize directional intensity so a plane normal to the central
        # ray receives the requested inverse-square irradiance.
        pressure = (K.SOLAR_CONSTANT/K.SPEED_OF_LIGHT*(K.AU/distance)**2/
                    float(weights@(directions@centre_direction)))
        return dict(directions=directions, weights=weights, pressure_Pa=pressure,
                    source_directions=directions, source_weights=weights,
                    perturbing_bodies=[dict(name='earth', position_m=earth, gm_m3_s2=K.EARTH_GM),
                                       dict(name='sun', position_m=sun, gm_m3_s2=K.SUN_GM)],
                    occluders=[dict(name='moon', position_m=np.zeros(3), radius_m=K.MOON_RADIUS),
                               dict(name='earth', position_m=earth, radius_m=K.EARTH_RADIUS)],
                    assumptions='Moon-origin parallel solar rays; uniform solar disk; Earth/Sun differential point gravity; spherical opaque Earth/Moon; no parallax, light-time, lunar harmonics or other planets')
    return at
