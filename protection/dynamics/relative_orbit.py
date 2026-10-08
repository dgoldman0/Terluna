"""Relative orbits of shield tiles about the Moon.

A formation repeats only when its members share an orbital energy. Tiles whose
osculating semi-major axes differ by da drift apart along-track by about
3*pi*da per revolution. The linear Hill-Clohessy-Wiltshire solution about a
circular reference, written in quasi-nonsingular relative orbital elements
(da, dlambda, de_x, de_y, di_x, di_y, all divided by the reference semi-major
axis), gives the repeating relative orbits (da = 0) and their closed-form
minimum separation in the radial-normal plane (D'Amico and Montenbruck 2006,
J. Guid. Control Dyn. 29, 554-563). With de parallel to di the radial
separation peaks where the cross-track separation vanishes.

Frames: R is radial, T completes the right-handed set in the orbit plane (the
along-track direction for a circular orbit), N is along the orbit normal.
Relative velocities are rates seen in the rotating RTN frame. The argument of
latitude u is measured in the reference orbit plane from a chosen origin; the
elements depend on that origin only through the phases of de and di.
"""
from __future__ import annotations

import numpy as np

from shared import constants as K


def semi_major_axis(states, mu=K.MOON_GM):
    """Osculating two-body semi-major axis of Moon-relative states (..., 6)."""
    states = np.asarray(states, float)
    r = np.linalg.norm(states[..., :3], axis=-1)
    v2 = np.sum(states[..., 3:6]**2, axis=-1)
    return 1/(2/r-v2/mu)


def period(a, mu=K.MOON_GM):
    return 2*np.pi*np.sqrt(np.asarray(a, float)**3/mu)


def drift_per_revolution(a, reference, mu=K.MOON_GM):
    """Along-track arc gained on the reference during one reference period.

    Negative values trail: a larger orbit has a longer period. For small da
    the value approaches -3*pi*da.
    """
    a = np.asarray(a, float)
    return 2*np.pi*reference*((reference/a)**1.5-1)


def energy_matching_burns(states, target, mu=K.MOON_GM, iterations=4):
    """Velocity changes along each velocity that set the semi-major axis to target.

    target is a scalar or one value per state. Vis-viva gives
    da = 2 a^2 v dv / mu; a few Newton steps converge to rounding.
    """
    s = np.array(states, float)
    total = np.zeros_like(s[..., 3:6])
    for _ in range(iterations):
        a = semi_major_axis(s, mu)
        v = np.linalg.norm(s[..., 3:6], axis=-1)
        dv = mu*(target-a)/(2*a**2*v)
        step = (dv/v)[..., None]*s[..., 3:6]
        s[..., 3:6] += step
        total += step
    return total


def rtn_frame(reference):
    """Rows R, T, N of the reference state's orbital frame."""
    r, v = np.asarray(reference[:3], float), np.asarray(reference[3:6], float)
    radial = r/np.linalg.norm(r)
    normal = np.cross(r, v)
    normal /= np.linalg.norm(normal)
    return np.vstack([radial, np.cross(normal, radial), normal])


def rtn_relative(reference, states):
    """Relative states (..., 6) in the rotating RTN frame of the reference state."""
    reference = np.asarray(reference, float)
    states = np.asarray(states, float)
    frame = rtn_frame(reference)
    r, v = reference[:3], reference[3:6]
    omega = np.cross(r, v)/np.dot(r, r)
    dr = states[..., :3]-r
    dv = states[..., 3:6]-v-np.cross(omega, dr)
    return np.concatenate([dr@frame.T, dv@frame.T], axis=-1)


def hcw_elements(relative, a, n, u):
    """Relative orbital elements from RTN relative states at argument of latitude u."""
    x, y, z, xd, yd, zd = np.moveaxis(np.asarray(relative, float), -1, 0)
    da = 4*x/a+2*yd/(a*n)
    radial = 3*x/a+2*yd/(a*n)
    rate = xd/(a*n)
    c, s = np.cos(u), np.sin(u)
    return np.stack([da, y/a-2*rate, radial*c+rate*s, radial*s-rate*c,
                     z/a*s+zd/(a*n)*c, -z/a*c+zd/(a*n)*s], axis=-1)


def hcw_state(elements, a, n, u0, dt):
    """RTN relative state dt after the epoch of the elements, where u = u0."""
    da, dl, dex, dey, dix, diy = np.moveaxis(np.asarray(elements, float), -1, 0)
    u = u0+n*np.asarray(dt, float)
    c, s = np.cos(u), np.sin(u)
    return np.stack([a*(da-dex*c-dey*s),
                     a*(dl-1.5*n*da*dt+2*(dex*s-dey*c)),
                     a*(dix*s-diy*c),
                     a*n*(dex*s-dey*c),
                     a*(-1.5*n*da+2*n*(dex*c+dey*s)),
                     a*n*(dix*c+diy*s)], axis=-1)


def minimum_rn_separation(de, di, a):
    """Smallest radial-normal distance over one orbit of a relative orbit with da = 0.

    de, di (..., 2) are relative eccentricity and inclination vectors. In the
    linear solution x = -a|de|cos(u-phi) and z = a|di|sin(u-theta); the
    minimum is exact for that motion.
    """
    de, di = np.asarray(de, float), np.asarray(di, float)
    A = a*np.linalg.norm(de, axis=-1)
    B = a*np.linalg.norm(di, axis=-1)
    gap = np.arctan2(de[..., 1], de[..., 0])-np.arctan2(di[..., 1], di[..., 0])
    C = np.sqrt(np.maximum(A**4+B**4-2*A**2*B**2*np.cos(2*gap), 0.))
    return np.sqrt(np.maximum((A**2+B**2-C)/2, 0.))


def normal_thickness(states):
    """Spread of positions along the mean orbit normal of a set of states (N, 6)."""
    states = np.asarray(states, float)
    centre = states.mean(axis=0)
    normal = rtn_frame(centre)[2]
    offsets = (states[:, :3]-centre[:3])@normal
    return float(offsets.max()-offsets.min())
