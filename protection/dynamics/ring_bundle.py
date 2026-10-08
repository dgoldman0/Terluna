"""Bundles of nested ring segments with per-tile attitudes.

A ring is one natural orbit carrying a string of tiles; rings step outward in
radius and slide past one another. Each tile keeps its own radial-facing
attitude, referenced to its velocity and tilted about the orbit normal so that
neighbours on a ring overlap like shingles. Tiles on one ring start as
time-shifted copies of the ring's centre trajectory, so the tide and the
sail force act on them alike.

The sail force reuses the compiled finite-Sun shadow geometry of
eclipse_parallel and fast_parallel, evaluated in the bundle centroid's frame;
attitudes across a bundle differ by well under a degree. Incidence and force
direction use each tile's own normal.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree

from shared import constants as K
from .collection import require_uneclipsed_tiles
from .cycling import ray_geometry
from .eclipse_parallel import moments as masked_moments
from .fast_parallel import moments as parallel_moments
from .fleet import sun_points
from .natural_pattern import finite_gravity, retarded_sun
from .optical import length, unit

SIDE, CLEAR = 10000., 9890.


def tile_frames(states, tilt):
    """Per-tile axes (..., 3, 3), columns u, v, n: radial-facing from each velocity, tilted about the normal."""
    states = np.asarray(states, float)
    v = unit(states[..., 3:6])
    n = unit(np.cross(states[..., :3], states[..., 3:6]))
    r = np.cross(v, n)
    s, c = np.sin(tilt), np.cos(tilt)
    return np.stack([s*r+c*v, n, c*r-s*v], axis=-1)


def centroid_frame(states, tilt):
    return tile_frames(np.asarray(states, float).mean(axis=0), tilt)


def eclipse_bodies(q, sample):
    """Moon and Earth blocker vectors when any tile's envelope meets a body shadow, else None."""
    rays = ray_geometry(q, sample)
    try:
        require_uneclipsed_tiles(rays)
        return None
    except ValueError as exc:
        if 'Body eclipse' not in str(exc):
            raise
    return np.stack([rays['moon'], rays['earth']], axis=1)


def visible_moments(env, t, q, frame, source, sample, bodies=False):
    """Visible clear area and first moments per tile for one solar source point."""
    if bodies is False:
        bodies = eclipse_bodies(q, sample)
    if bodies is None:
        return parallel_moments(q, frame, source)
    return masked_moments(q, frame, source, bodies, [K.MOON_RADIUS, K.EARTH_RADIUS])


class RingBundleForces:
    """Gravity with finite-square quadrature plus the sail force on per-tile normals."""
    def __init__(self, env, tilt, mass_ratio, suns=8, rotation=.317):
        self.env, self.tilt, self.suns, self.rotation = env, tilt, suns, rotation
        self.mass = env.sigma*SIDE**2*np.atleast_1d(np.asarray(mass_ratio, float))

    def light(self, t, states):
        q = states[:, :3]
        frames = tile_frames(states, self.tilt)
        common = centroid_frame(states, self.tilt)
        sample = self.env.at(t)
        normals = frames[:, :, 2]
        force = np.zeros_like(q)
        first = np.zeros(len(q))
        bare = np.zeros(len(q))
        bodies = eclipse_bodies(q, sample)
        for source in sun_points(retarded_sun(sample), self.suns, self.rotation):
            area = visible_moments(self.env, t, q, common, source, sample, bodies)[:, 0]
            ray = source-q
            cosine = np.sum(unit(ray)*normals, axis=1)
            density = K.SOLAR_CONSTANT*(K.AU/length(ray))**2/self.suns*np.abs(cosine)
            first += density*area
            bare += density*CLEAR**2
            pressure = -2*self.env.central_fraction/K.SPEED_OF_LIGHT*density*cosine
            force += (pressure*area)[:, None]*normals
        return force, first, bare

    def acceleration(self, t, states):
        frames = tile_frames(states, self.tilt)
        force, _, _ = self.light(t, states)
        gravity = finite_gravity(self.env, t, states[:, :3], frames[:, :, 0], frames[:, :, 1])
        return gravity+force/self.mass[:, None]


def propagate(acceleration, states, t0, t1, rtol=1e-10, atol=1e-3, max_step=600.):
    n = len(states)

    def rhs(t, y):
        s = y.reshape(n, 6)
        return np.hstack([s[:, 3:], acceleration(t, s)]).ravel()
    sol = solve_ivp(rhs, (t0, t1), np.asarray(states, float).ravel(), method='DOP853', rtol=rtol,
                    atol=atol, max_step=max_step, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return (lambda t: sol.sol(t).reshape(n, 6)), sol


def time_track(acceleration, centre, lags, rtol=1e-10, atol=1e-3, max_step=600.):
    """States one trajectory reaches after each lag (s); the centre is its state at epoch zero.

    Positive lags integrate forward from zero and negative lags backward, so a zero lag returns the
    centre itself and every lag keeps its time relative to the epoch.
    """
    centre = np.asarray(centre, float)
    lags = np.asarray(lags, float)
    states = np.repeat(centre[None], len(lags), axis=0)
    for side in (lags > 0, lags < 0):
        if side.any():
            end = lags[side][np.argmax(np.abs(lags[side]))]
            at, _ = propagate(acceleration, centre[None], 0., end+np.sign(end), rtol, atol, max_step)
            states[side] = [at(lag)[0] for lag in lags[side]]
    return states


def build(r0, centre_track, rings, per_ring, pitch, step, height_pitch, sun_plane, normal, start_u):
    """Initial states of a ring bundle about a reference orbit plane.

    r0: the middle ring's radius. sun_plane: unit Sun direction projected into the reference plane; normal:
    the plane's unit normal; start_u: the bundle centre's angle from sun_plane
    toward the motion. Ring k (0..rings-1) has radius r0+(k-mid)*step and peaks
    (k-mid)*height_pitch out of plane at the sub-solar point. centre_track(state,
    lags) returns the states reached from a ring centre after the given lags;
    tiles are those states, so each ring is one trajectory shifted in time.
    """
    w = np.cross(normal, sun_plane)
    mid = (rings-1)/2
    states, ring, slot = [], [], []
    for k in range(rings):
        radius = r0+(k-mid)*step
        tilt = np.arcsin((k-mid)*height_pitch/radius)
        s_k = np.cos(tilt)*sun_plane+np.sin(tilt)*normal
        speed = np.sqrt(K.MOON_GM/radius)
        centre = np.r_[radius*(np.cos(start_u)*s_k+np.sin(start_u)*w),
                       speed*(-np.sin(start_u)*s_k+np.cos(start_u)*w)]
        lags = (np.arange(per_ring)-(per_ring-1)/2)*pitch/speed
        states.extend(centre_track(centre, lags))
        ring.extend([k]*per_ring)
        slot.extend(range(per_ring))
    return np.array(states), np.array(ring), np.array(slot)


def stumpff(z):
    """Stumpff functions C(z) and S(z) of the universal-variable two-body solution."""
    z = np.asarray(z, float)
    c, s = np.empty_like(z), np.empty_like(z)
    pos, neg = z > 1e-8, z < -1e-8
    small = ~(pos | neg)
    r = np.sqrt(z[pos])
    c[pos], s[pos] = (1-np.cos(r))/z[pos], (r-np.sin(r))/r**3
    r = np.sqrt(-z[neg])
    c[neg], s[neg] = (np.cosh(r)-1)/-z[neg], (np.sinh(r)-r)/r**3
    w = z[small]
    c[small], s[small] = 1/2-w/24+w*w/720, 1/6-w/120+w*w/5040
    return c, s


def kepler_shift(states, dt, mu=K.MOON_GM, iterations=10):
    """Two-body states (..., 6) advanced by dt (...), from the universal Kepler equation."""
    states = np.asarray(states, float)
    r0, v0 = states[..., :3], states[..., 3:6]
    dt = np.broadcast_to(np.asarray(dt, float), r0.shape[:-1])
    rn = np.linalg.norm(r0, axis=-1)
    vr = np.sum(r0*v0, axis=-1)/rn
    alpha = 2/rn-np.sum(v0*v0, axis=-1)/mu
    root = np.sqrt(mu)
    chi = root*np.abs(alpha)*dt
    for _ in range(iterations):
        z = alpha*chi**2
        c, s = stumpff(z)
        f = rn*vr/root*chi**2*c+(1-alpha*rn)*chi**3*s+rn*chi-root*dt
        df = rn*vr/root*chi*(1-z*s)+(1-alpha*rn)*chi**2*c+rn
        chi = chi-f/df
    z = alpha*chi**2
    c, s = stumpff(z)
    f, g = 1-chi**2/rn*c, dt-chi**3/root*s
    r = f[..., None]*r0+g[..., None]*v0
    r_n = np.linalg.norm(r, axis=-1)
    fd, gd = root/(r_n*rn)*(alpha*chi**3*s-chi), 1-chi**2/r_n*c
    return np.concatenate([r, fd[..., None]*r0+gd[..., None]*v0], axis=-1)


class ContinuousRings:
    """One representative tile per ring; its ring and the neighbouring rings appear as rotated copies.

    A complete ring is one trajectory shifted in time, so a tile's neighbours on
    the ring and on the adjacent rings are the representatives' states advanced
    or delayed by whole tile lags. Over lags of minutes the two-body motion
    carries that shift along the ring's eccentric orbit; the tide and sail force
    move a tile by metres in that time. The copies cast shadows and enter
    clearance and overlap checks; only the representatives are propagated. Rings
    0 and the last are the bundle's edges.
    """
    def __init__(self, env, tilt, mass_ratio, lags, reach=3, suns=8, rotation=.317):
        self.forces = RingBundleForces(env, tilt, mass_ratio, suns, rotation)
        self.env, self.tilt, self.lags, self.reach = env, tilt, np.asarray(lags, float), reach

    def copies(self, states):
        """Rotated copies around each representative: (states, ring, index) with index 0 the representative."""
        normals = unit(np.cross(states[:, :3], states[:, 3:6]))
        rate = np.linalg.norm(np.cross(states[:, :3], states[:, 3:6]), axis=1)/np.sum(states[:, :3]**2, axis=1)
        wanted = set()
        for k in range(len(states)):
            wanted.update((k, j) for j in range(-self.reach, self.reach+1) if j)
            for m in (k-1, k+1):
                if 0 <= m < len(states):
                    # Phase of representative k seen from ring m, in ring m's tile lags.
                    rel = np.arctan2(np.dot(normals[m], np.cross(states[m, :3], states[k, :3])),
                                     np.dot(states[m, :3], states[k, :3]))
                    j0 = int(np.round(rel/(rate[m]*self.lags[m])))
                    wanted.update((m, j) for j in range(j0-self.reach, j0+self.reach+1) if j)
        pairs = np.array(sorted(wanted))
        m, j = pairs[:, 0], pairs[:, 1]
        copies = kepler_shift(states[m], j*self.lags[m])
        return (np.concatenate([states, copies]), np.r_[np.arange(len(states)), m],
                np.r_[np.zeros(len(states), int), j])

    def acceleration(self, t, states):
        everything, ring, index = self.copies(states)
        force, _, _ = self.forces.light(t, everything)
        frames = tile_frames(states, self.tilt)
        gravity = finite_gravity(self.env, t, states[:, :3], frames[:, :, 0], frames[:, :, 1])
        return gravity+force[:len(states)]/self.forces.mass[:len(states), None]


def own_attitude_clearance(states, tilt, side=SIDE):
    """Smallest surface distance between nearly parallel squares, each pair in its mean frame."""
    q = states[:, :3]
    frames = tile_frames(states, tilt)
    pairs = cKDTree(q).query_pairs(np.sqrt(2)*side+500., output_type='ndarray')
    if len(pairs) == 0:
        return np.inf, (-1, -1)
    f = frames[pairs[:, 0]]+frames[pairs[:, 1]]
    f /= np.linalg.norm(f, axis=1)[:, None, :]
    d = np.einsum('ki,kij->kj', q[pairs[:, 0]]-q[pairs[:, 1]], f)
    gap = np.linalg.norm(np.maximum(np.abs(d)-[side, side, 0.], 0.), axis=1)
    k = int(np.argmin(gap))
    return float(gap[k]), (int(pairs[k, 0]), int(pairs[k, 1]))
