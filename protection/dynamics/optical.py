"""Finite-Sun occultation and momentum bounds for a Sun-facing aperture."""
from __future__ import annotations
import numpy as np
from shared import constants as K


def length(x):
    return np.linalg.norm(x, axis=-1)


def unit(x):
    return x / length(x)[..., None]


def disk_overlap(r1, r2, separation):
    """Area of two disks; angular-plane radii are sufficient at these angles."""
    r1, r2, d = np.broadcast_arrays(r1, r2, separation)
    safe = np.maximum(d, np.finfo(float).tiny)
    c1 = np.clip((d*d + r1*r1 - r2*r2) / (2 * safe * r1), -1, 1)
    c2 = np.clip((d*d + r2*r2 - r1*r1) / (2 * safe * r2), -1, 1)
    radicand = np.maximum((-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2), 0)
    area = r1*r1*np.arccos(c1) + r2*r2*np.arccos(c2) - 0.5*np.sqrt(radicand)
    area = np.where(d <= np.abs(r1-r2), np.pi*np.minimum(r1, r2)**2, area)
    return np.where(d >= r1+r2, 0, area)


def solar_visibility(sun_vector, earth_vector):
    """Uniform-brightness solar disk blocked by Earth's spherical solid body."""
    ds, de = length(sun_vector), length(earth_vector)
    if np.any(de <= K.EARTH_RADIUS):
        raise ValueError("Point is inside Earth")
    rs, re = np.arcsin(K.SUN_RADIUS/ds), np.arcsin(K.EARTH_RADIUS/de)
    separation = np.arctan2(length(np.cross(sun_vector, earth_vector)),
                           np.sum(sun_vector*earth_vector, axis=-1))
    blocked = disk_overlap(rs, re, separation) / (np.pi*rs*rs)
    return np.where(de < ds, np.clip(1-blocked, 0, 1), 1)


def arriving_ray_vectors(sun_vector, earth_vector, sun_velocity, earth_velocity, sun_acceleration, earth_acceleration):
    """Retarded source and blocker vectors for photons reaching the screen now.

    Taylor-expand barycentric body motion to second order. Iterate the solar
    travel time; the blocker time is its along-ray distance/c. Exact kernel
    evaluations validate this approximation separately. Observer aberration
    is common to both rays; forces are expressed in the barycentric axes.
    """
    tau=length(sun_vector)/K.SPEED_OF_LIGHT
    for _ in range(2):
        sun=sun_vector-sun_velocity*tau[...,None]+0.5*sun_acceleration*tau[...,None]**2
        tau=length(sun)/K.SPEED_OF_LIGHT
    tau_e=np.maximum(np.sum(earth_vector*unit(sun),axis=-1),0)/K.SPEED_OF_LIGHT
    earth=earth_vector-earth_velocity*tau_e[...,None]+0.5*earth_acceleration*tau_e[...,None]**2
    return sun,earth


def optical_projection(required, sun_vector, sigma, diverted_fraction=1.0, visibility=1.0):
    """Closest acceleration in the relaxed photon-momentum ball.

    A diverted fraction q carries incoming momentum b*e, with b=q*S/(c*sigma).
    Any passive outgoing direction, or mixtures of directions, gives a ball
    centred at b*e with radius b. Its sunward support is zero. Unchanged forward
    transmission contributes no momentum. This relaxes real angular/spectral
    constraints, so the residual is a LOWER BOUND, not an implemented sail law.
    """
    sigma = np.asarray(sigma)
    diverted_fraction = np.asarray(diverted_fraction)
    if np.any(sigma <= 0) or np.any((diverted_fraction < 0) | (diverted_fraction > 1)):
        raise ValueError("Positive areal mass and a fractional diverted spectrum required")
    flux = K.SOLAR_CONSTANT * (K.AU/length(sun_vector))**2 * np.asarray(visibility)
    b = diverted_fraction * flux / (K.SPEED_OF_LIGHT*sigma)
    e = -unit(sun_vector)
    center = b[..., None]*e
    delta = required-center
    mag = length(delta)
    fraction = np.divide(b, mag, out=np.ones_like(mag), where=mag > 0)
    optical = center + delta*np.minimum(fraction, 1)[..., None]
    return optical, required-optical


def axial_optical(sun_vector, sigma, coefficient, visibility=1.0):
    flux = K.SOLAR_CONSTANT*(K.AU/length(sun_vector))**2*np.asarray(visibility)
    return -unit(sun_vector)*(coefficient*flux/(K.SPEED_OF_LIGHT*sigma))[..., None]


def covering_radius(distance, solar_distance, protected_radius, transverse_offset=0):
    """Conservative radius in a plane normal to the Moon–Sun direction.

    Straight rays through the full solar disk intersect the screen plane.
    A displaced screen centre adds its transverse offset to the aperture.
    """
    distance = np.asarray(distance)
    if np.any(distance <= 0) or np.any(distance >= solar_distance):
        raise ValueError("The screen must lie strictly between Moon and Sun")
    return protected_radius*(1-distance/solar_distance) + K.SUN_RADIUS*distance/solar_distance + transverse_offset


def light_time_allowance(distance, moon_barycentric_velocity, sun_barycentric_velocity):
    """Conservative first-order shift from target and source motion.

    The moving Moon reaches the shadow after d/c. The retarded Sun shifts the
    incident ray direction by order v_sun/c. Full speed bounds their transverse
    components. The separate 50 m formation margin exceeds the omitted
    acceleration terms for the present 20,000–180,000 km screen family.
    """
    return np.asarray(distance)/K.SPEED_OF_LIGHT*(length(moon_barycentric_velocity)+length(sun_barycentric_velocity))
