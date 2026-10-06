"""Impulsive keeping of nested rings relative to a reference ring.

Each ring's representative is held to designed relative orbital elements
about the reference ring: its semi-major-axis step, the reference's
eccentricity vector and its own designed plane. Rings slide past one another,
so the along-track element is left free. Once per orbit two tangential burns
half an orbit apart restore the semi-major-axis step and the eccentricity
vector, and one normal burn restores the plane (the near-circular impulsive
solution of D'Amico and Montenbruck's relative-orbit control). The linear
model's errors at an eccentricity of about 0.02 are corrected on later orbits.
"""
from __future__ import annotations

import numpy as np

from shared import constants as K
from .relative_orbit import hcw_elements, rtn_frame, rtn_relative, semi_major_axis


def argument_of_latitude(reference, origin):
    """Angle of the reference state from origin projected into its orbit plane, toward the motion."""
    R, T, N = rtn_frame(reference)
    o = origin-np.dot(origin, N)*N
    o /= np.linalg.norm(o)
    return float(np.arctan2(np.dot(np.cross(o, R), N), np.dot(o, R)))


def elements(reference, states, origin, mu=K.MOON_GM):
    """Relative orbital elements of states about the reference, with u measured from origin."""
    a = semi_major_axis(reference[None], mu)[0]
    n = np.sqrt(mu/a**3)
    u = argument_of_latitude(reference, origin)
    return hcw_elements(rtn_relative(reference, states), a, n, u), a, n, u


def geometric_elements(reference, states, origin, mu=K.MOON_GM):
    """Exact relative elements about the reference: (da, 0, dex, dey, dix, diy), divided by a.

    Semi-major axes come from vis-viva; the eccentricity-vector difference is
    projected on the reference plane's axes (origin direction and the axis
    90 degrees ahead); the plane difference follows z = a(dix sin u - diy cos u)
    from the unit orbit normals. These carry no linearisation about a circular
    reference, so offsets of tens of kilometres stay free of eccentricity
    artefacts.
    """
    states = np.asarray(states, float)
    a = semi_major_axis(states, mu)
    a0 = semi_major_axis(reference[None], mu)[0]
    h = np.cross(states[:, :3], states[:, 3:6])
    e = np.cross(states[:, 3:6], h)/mu-states[:, :3]/np.linalg.norm(states[:, :3], axis=1)[:, None]
    hr = np.cross(reference[:3], reference[3:6])
    er = np.cross(reference[3:6], hr)/mu-reference[:3]/np.linalg.norm(reference[:3])
    N = hr/np.linalg.norm(hr)
    ox = origin-np.dot(origin, N)*N
    ox /= np.linalg.norm(ox)
    oy = np.cross(N, ox)
    dh = h/np.linalg.norm(h, axis=1)[:, None]-N
    de = e-er
    n0 = np.sqrt(mu/a0**3)
    return (np.stack([(a-a0)/a0, np.zeros(len(a)), de@ox, de@oy, -(dh@oy), dh@ox], axis=1), a0, n0,
            argument_of_latitude(reference, origin))


def plan(error, a, n):
    """Burns that remove an element error (target minus current, divided by a).

    Returns a list of (phase, (dv_radial, dv_along, dv_normal)) with phases as
    arguments of latitude. Two tangential burns at phase p and p+pi change the
    semi-major axis by their sum and the eccentricity vector along p by their
    difference; a normal burn at the phase of the inclination correction turns
    the plane.
    """
    da, _, dex, dey, dix, diy = error
    burns = []
    de = np.hypot(dex, dey)
    p = np.arctan2(dey, dex) if de > 0 else 0.
    t1 = n*a*(da+de)/4
    t2 = n*a*(da-de)/4
    burns += [(p, (0., t1, 0.)), (p+np.pi, (0., t2, 0.))]
    di = np.hypot(dix, diy)
    if di > 0:
        burns.append((np.arctan2(diy, dix), (0., 0., n*a*di)))
    return burns


def apply(state, impulse):
    """Add an RTN impulse (radial, along-track, normal) to one state."""
    frame = rtn_frame(state)
    out = np.array(state, float)
    out[3:] += np.asarray(impulse, float)@frame
    return out


def gravity_shift(env, t, state, dt, step=15.):
    """A state carried dt along its orbit under the environment's gravity (fixed-step RK4).

    Over the minutes between neighbouring rings this follows the tide's
    short-period changes of the osculating orbit, which a two-body shift holds
    fixed; the sail force over that time moves the elements by metres.
    """
    count = max(1, int(np.ceil(abs(dt)/step)))
    h = dt/count
    y = np.array(state, float)

    def f(time, x):
        return np.r_[x[3:], env.gravity(time, x[None, :3])[0]]
    for i in range(count):
        s = t+i*h
        k1 = f(s, y)
        k2 = f(s+h/2, y+h/2*k1)
        k3 = f(s+h/2, y+h/2*k2)
        k4 = f(s+h, y+h*k3)
        y = y+h/6*(k1+2*k2+2*k3+k4)
    return y


class ContinuousKeeping:
    """Low-thrust feedback holding each ring to its neighbour toward a free reference ring.

    Ring k is compared with its neighbour on the reference side, carried along
    that neighbour's orbit to ring k's argument of latitude under the
    environment's gravity (gravity_shift), so the tide's short-period swing of
    the osculating orbit is matched as well.
    Neighbouring rings slide by only kilometres per orbit, so the shift stays
    short and the comparison avoids the tide's short-period differences between
    distant phases, which reach about a kilometre after two orbits for rings
    eight apart. Errors are geometric_elements against the same quantities at
    the start. The errors (semi-major-axis step, eccentricity vector, plane)
    decay with time constant tau through Gauss's near-circular equations: the
    acceleration is the transposed influence matrix times the scaled error, which
    reduces the squared error monotonically in the linear model. The along-track
    element is free, the reference ring flies free, and accelerations are capped
    at limit (m/s^2).
    """
    def __init__(self, model, start, reference, origin, tau, limit=1e-3, mu=K.MOON_GM):
        self.model, self.reference, self.origin = model, reference, origin
        self.tau, self.limit, self.mu = tau, limit, mu
        self.targets = self.relative(0., np.asarray(start, float))

    def neighbour(self, k):
        return k-1 if k > self.reference else k+1

    def relative(self, t, states):
        """Elements of each ring about its phase-matched neighbour (the reference row stays zero)."""
        out = np.zeros((len(states), 6))
        for k in range(len(states)):
            if k == self.reference:
                continue
            m = self.neighbour(k)
            h = np.cross(states[m, :3], states[m, 3:6])
            rate = np.linalg.norm(h)/np.dot(states[m, :3], states[m, :3])
            ahead = np.arctan2(np.dot(h/np.linalg.norm(h), np.cross(states[m, :3], states[k, :3])),
                               np.dot(states[m, :3], states[k, :3]))
            matched = gravity_shift(self.model.env, t, states[m], ahead/rate)
            out[k] = geometric_elements(matched, states[k][None], self.origin, self.mu)[0][0]
        return out

    def control(self, t, states):
        a = semi_major_axis(states[self.reference][None], self.mu)[0]
        n = np.sqrt(self.mu/a**3)
        error = (self.targets-self.relative(t, states))*a
        error[:, 1] = 0.
        out = np.zeros((len(states), 3))
        for k in range(len(states)):
            if k == self.reference:
                continue
            u = argument_of_latitude(states[k], self.origin)
            c, s = np.cos(u), np.sin(u)
            # Rows: d(a da), d(a dex), d(a dey), d(a dix), d(a diy) per unit (aR, aT, aN), each divided by n.
            B = np.array([[0., 2., 0.], [s, 2*c, 0.], [-c, 2*s, 0.], [0., 0., c], [0., 0., s]])/n
            rtn = error[k, [0, 2, 3, 4, 5]]@B/(self.tau*np.sum(B*B, axis=0))
            norm = np.linalg.norm(rtn)
            if norm > self.limit:
                rtn *= self.limit/norm
            out[k] = rtn@rtn_frame(states[k])
        return out

    def acceleration(self, t, states):
        return self.model.acceleration(t, states)+self.control(t, states)
