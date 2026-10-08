"""Returning solar-sail trajectories and conservative shadow service geometry.

This module is separate from the pinned held-aperture experiment. Its optical
actuator is an ideal, steerable specular membrane, not the relaxed momentum
ball used by that experiment. No electric thrust is added to propagation.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from scipy.optimize import brentq

from shared import constants as K
from .ephemeris import BODY_GM
from .model import normalized_derivatives
from .optical import length, unit, disk_overlap, arriving_ray_vectors


def lunar_frame(samples):
    """Instantaneous Earth–Moon orbital triad and its inertial derivative."""
    p, v, a = (-samples["positions"]["earth"],
               -samples["velocities"]["earth"], -samples["earth_relative_a"])
    x, xd, _ = normalized_derivatives(p, v, a)
    h = np.cross(p, v)
    hd = np.cross(p, a)
    z = unit(h)
    zd = (hd-z*np.sum(z*hd, axis=-1)[..., None])/length(h)[..., None]
    y, yd = np.cross(z, x), np.cross(zd, x)+np.cross(z, xd)
    return np.stack([x, y, z], axis=-1), np.stack([xd, yd, zd], axis=-1)


def to_inertial(state, frame, frame_d):
    q = np.einsum("...ij,...j->...i", frame, state[..., :3])
    v = (np.einsum("...ij,...j->...i", frame, state[..., 3:])+
         np.einsum("...ij,...j->...i", frame_d, state[..., :3]))
    return np.concatenate([q, v], axis=-1)


def from_inertial(state, frame, frame_d):
    q = np.einsum("...ji,...j->...i", frame, state[..., :3])
    v = np.einsum("...ji,...j->...i", frame,
                  state[..., 3:]-np.einsum("...ij,...j->...i", frame_d, q))
    return np.concatenate([q, v], axis=-1)


def visible_sun(sun, earth, moon, quadrature=96):
    """Uniform finite solar disk, occulted by the union of Earth and Moon.

    Exact two-circle overlap is used unless both silhouettes overlap the Sun
    and one another. That uncommon three-disk case uses Gauss–Legendre strips
    across the solar disk; an occulted ray is counted once. Angular-plane
    circles neglect spherical projection corrections of order solar angle².
    """
    sun, earth, moon = np.broadcast_arrays(sun, earth, moon)
    shape = sun.shape[:-1]
    s, e, m = (x.reshape(-1, 3) for x in (sun, earth, moon))
    ds, de, dm = length(s), length(e), length(m)
    if np.any(de <= K.EARTH_RADIUS) or np.any(dm <= K.MOON_RADIUS):
        raise ValueError("Sail intersects a solid primary")
    u = unit(s)
    rs = np.arcsin(K.SUN_RADIUS/ds)
    re, rm = np.arcsin(K.EARTH_RADIUS/de), np.arcsin(K.MOON_RADIUS/dm)
    se = np.arctan2(length(np.cross(s, e)), np.sum(s*e, axis=-1))
    sm = np.arctan2(length(np.cross(s, m)), np.sum(s*m, axis=-1))
    ae = disk_overlap(rs, re, se)*(de < ds)
    am = disk_overlap(rs, rm, sm)*(dm < ds)
    sep = np.arctan2(length(np.cross(e, m)), np.sum(e*m, axis=-1))
    both = (ae > 0) & (am > 0) & (sep < re+rm)
    blocked = ae+am
    if np.any(both):
        nodes, weights = np.polynomial.legendre.leggauss(quadrature)
        for i in np.flatnonzero(both):
            if ae[i] >= np.pi*rs[i]**2*(1-1e-12) or am[i] >= np.pi*rs[i]**2*(1-1e-12):
                blocked[i] = np.pi*rs[i]**2
                continue
            tangent = e[i]-u[i]*np.dot(e[i], u[i])
            if length(tangent) < 1e-10:
                tangent = np.cross(u[i], [0., 0., 1.])
            b1 = unit(tangent)
            b2 = np.cross(u[i], b1)
            centers = []
            for body, angle in ((e[i], se[i]), (m[i], sm[i])):
                v = body-u[i]*np.dot(body, u[i])
                direction = v/max(length(v), 1e-100)
                centers.append(angle*np.array([np.dot(direction, b1), np.dot(direction, b2)]))
            # Split at all circle x-extrema to avoid integrating over kinks.
            edges = [-rs[i], rs[i]]
            for (cx, _), radius in zip(centers, (re[i], rm[i])):
                edges += list(np.clip([cx-radius, cx+radius], -rs[i], rs[i]))
            area = 0.
            for left, right in zip(sorted(set(edges))[:-1], sorted(set(edges))[1:]):
                x = (left+right)/2+(right-left)/2*nodes
                sy = np.sqrt(np.maximum(rs[i]**2-x*x, 0))
                intervals = []
                for (cx, cy), radius in zip(centers, (re[i], rm[i])):
                    dy = np.sqrt(np.maximum(radius**2-(x-cx)**2, 0))
                    lo, hi = np.maximum(-sy, cy-dy), np.minimum(sy, cy+dy)
                    hi = np.maximum(lo, hi)
                    intervals.append((lo, hi))
                (l1, h1), (l2, h2) = intervals
                union = h1-l1+h2-l2-np.maximum(np.minimum(h1, h2)-np.maximum(l1, l2), 0)
                area += (right-left)/2*np.dot(weights, union)
            blocked[i] = area
    return np.clip(1-blocked/(np.pi*rs*rs), 0, 1).reshape(shape)


def ray_geometry(q, sample):
    """Retarded incident rays and future lunar intercept of the central ray.

    Source/blocker motion is expanded to second order in barycentric axes.
    Target motion during the tile-to-Moon light time is included. Bounds on
    the finite-Sun footprint classify a tile as central whenever any ray can
    reach the solid lunar disk; mixed central/annulus tiles get the stricter
    central spectrum. This does not resolve seams between finite-size tiles.
    """
    sun, earth = arriving_ray_vectors(
        sample["positions"]["sun"]-q, sample["positions"]["earth"]-q,
        sample["sun_v"], sample["earth_v"], sample["sun_a"], sample["earth_a"])
    u = unit(sun)
    moon = -q
    tau = np.maximum(np.sum(moon*u, axis=-1), 0)/K.SPEED_OF_LIGHT
    moon_blocker = moon-sample["moon_v"]*tau[..., None]+.5*sample["moon_a"]*tau[..., None]**2
    d = np.sum(q*u, axis=-1)
    tau_forward = np.maximum(d, 0)/K.SPEED_OF_LIGHT
    future_q = q-sample["moon_v"]*tau_forward[..., None]-.5*sample["moon_a"]*tau_forward[..., None]**2
    d = np.sum(future_q*u, axis=-1)
    transverse = future_q-u*d[..., None]
    spread = np.maximum(d, 0)*K.SUN_RADIUS/length(sun)
    return {"sun": sun, "earth": earth, "moon": moon_blocker,
            "sunward": d, "transverse": transverse, "spread": spread,
            "central": (d > 0) & (length(transverse)-spread <= K.MOON_RADIUS)}


def sail_basis(sun):
    e = -unit(sun)
    eps = np.deg2rad(K.EARTH_OBLIQUITY_DEG)
    north = np.array([0., -np.sin(eps), np.cos(eps)])
    b1 = unit(np.cross(north, e))
    return e, b1, np.cross(e, b1)


def specular_sail(sun, angles, sigma, fraction, visibility):
    """Ideal spectral mirror: a = 2 f S/(c sigma) cos²(alpha) n.

    Alpha is the sail-normal cone from the outgoing solar direction; clock
    angle turns the normal about that direction. Unreflected light transmits.
    Angles must stay within the passive, illuminated hemisphere. Force is per
    physical membrane mass, and projected useful area is cos(alpha).
    """
    alpha, clock = np.asarray(angles)[..., 0], np.asarray(angles)[..., 1]
    if np.any(np.abs(alpha) > np.pi/2+1e-12) or sigma <= 0:
        raise ValueError("Invalid sail cone or areal mass")
    e, b1, b2 = sail_basis(sun)
    cosine = np.cos(alpha)
    normal = (e*cosine[..., None]+np.sin(alpha)[..., None]*
              (b1*np.cos(clock)[..., None]+b2*np.sin(clock)[..., None]))
    magnitude = (2*np.asarray(fraction)*K.SOLAR_CONSTANT*(K.AU/length(sun))**2*
                 np.asarray(visibility)/(K.SPEED_OF_LIGHT*sigma)*cosine**2)
    return magnitude[..., None]*normal, cosine


def outgoing_clearance(q, rays, angles, protected_radius):
    """Conservative reflected central-ray clearance of a protected sphere.

    A specular reflection preserves the solar angular radius. Reject a state
    if its reflected beam cone can intersect the protected sphere. Positive
    clearance is angular (radians); a positive value certifies all reflected
    solar rays miss the sphere in the frozen geometry. Add a light-time motion
    margin separately when interpreting the result.
    """
    e, b1, b2 = sail_basis(rays["sun"])
    alpha, clock = angles[..., 0], angles[..., 1]
    n = e*np.cos(alpha)[..., None]+np.sin(alpha)[..., None]*(
        b1*np.cos(clock)[..., None]+b2*np.sin(clock)[..., None])
    outgoing = e-2*np.cos(alpha)[..., None]*n
    target = -unit(q)
    separation = np.arctan2(length(np.cross(outgoing, target)),
                           np.sum(outgoing*target, axis=-1))
    result = (separation-np.arcsin(np.minimum(protected_radius/length(q), 1.))-
              np.arcsin(K.SUN_RADIUS/length(rays["sun"])))
    return np.where(length(q) <= protected_radius, -np.pi, result)


def service_fraction(trace, protected_radius):
    """Projected membrane area whose incident rays cross the protected disk.

    This is a ray-area/time inventory statistic for an infinitesimal patch,
    not a fraction of the disk covered by one patch. Eclipsed incident rays
    and reflected rays that could return to the protected sphere are excluded.
    """
    spread = np.maximum(trace["spread"], 1.)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        overlap = disk_overlap(protected_radius, spread, length(trace["transverse"]))/(np.pi*spread**2)
    safe = outgoing_clearance(trace["state"][..., :3], trace, trace["angles"], protected_radius) > 0
    return (overlap*(trace["sunward"] > 0)*trace["projected"]*trace["visibility"]*safe)


class CyclingEnvironment:
    """Vectorized bounded ephemeris interpolator for many shooting evaluations."""
    def __init__(self, samples, sigma, central_fraction, spectrum="full_off_core"):
        self.sigma, self.central_fraction = sigma, central_fraction
        if spectrum not in ("filter_only", "full_off_core"):
            raise ValueError("Unknown optical spectrum")
        self.spectrum = spectrum
        self.names = list(BODY_GM)
        self.gm = np.array([BODY_GM[n] for n in self.names])
        t = samples["t"]
        self.limits = t[0], t[-1]
        p = np.stack([samples["positions"][n] for n in self.names], axis=-2)
        v = np.stack([samples["velocities"][n] for n in self.names], axis=-2)
        self.positions = CubicHermiteSpline(t, p, v, extrapolate=False)
        self.other = {k: CubicSpline(t, samples[k], extrapolate=False) for k in
                      ("moon_a", "moon_v", "sun_v", "earth_v", "sun_a", "earth_a")}
        self.frame, self.frame_d = lunar_frame(samples)
        self.frames = CubicHermiteSpline(t, self.frame, self.frame_d, extrapolate=False)

    def at(self, t):
        if np.any(np.asarray(t) < self.limits[0]) or np.any(np.asarray(t) > self.limits[1]):
            raise ValueError("Time outside ephemeris support")
        p = self.positions(t)
        return {"positions": {n: p[..., i, :] for i, n in enumerate(self.names)},
                **{k: s(t) for k, s in self.other.items()}}

    def gravity(self, t, q):
        p = self.positions(t)
        delta = p-q[..., None, :]
        return (-K.MOON_GM*q/length(q)[..., None]**3-self.other["moon_a"](t)+
                np.sum(self.gm[:, None]*delta/length(delta)[..., None]**3, axis=-2))

    def reflected_fraction(self, central):
        if self.spectrum == "filter_only":
            return np.full_like(central, self.central_fraction, dtype=float)
        return np.where(central, self.central_fraction, 1.)

    def acceleration(self, t, q, angles, enabled=True, diagnostics=False):
        gravity = self.gravity(t, q)
        if not enabled and not diagnostics:
            return gravity
        rays = ray_geometry(q, self.at(t))
        visibility = visible_sun(rays["sun"], rays["earth"], rays["moon"])
        fraction = self.reflected_fraction(rays["central"])
        sail, projected = specular_sail(rays["sun"], angles, self.sigma, fraction, visibility)
        if not enabled:
            sail = np.zeros_like(gravity)
        if diagnostics:
            return {**rays, "gravity": gravity, "sail": sail, "visibility": visibility,
                    "fraction": fraction, "projected": projected}
        return gravity+sail


def jacobi_value(local, distance):
    """Instantaneous circular-problem diagnostic, not an ephemeris invariant."""
    mu = K.MOON_GM/(K.EARTH_GM+K.MOON_GM)
    distance = np.asarray(distance)
    speed = np.sqrt((K.EARTH_GM+K.MOON_GM)/distance)
    p, v = local[..., :3]/distance[..., None], local[..., 3:]/speed[..., None]
    earth = p+np.array([1., 0., 0.])
    return ((p[..., 0]+1-mu)**2+p[..., 1]**2+
            2*((1-mu)/length(earth)+mu/length(p))-np.sum(v*v, axis=-1))


def service_coast_angles(state, rays, sample, protected_radius, feather_width=1e6):
    """Face the Sun in service; feather with reflection directed off the Moon."""
    margin = length(state[..., :3])/K.SPEED_OF_LIGHT*length(sample["moon_v"])+50.
    outside = (length(rays["transverse"])-protected_radius-rays["spread"]-margin)/feather_width
    f = np.where(rays["sunward"] > 0, np.clip(outside, 0., 1.), 1.)
    _, b1, b2 = sail_basis(rays["sun"])
    tangent = -rays["transverse"]
    alpha = np.pi/2*f*f*(3-2*f)
    clock = np.arctan2(np.sum(tangent*b2, axis=-1), np.sum(tangent*b1, axis=-1))
    return np.stack([alpha, clock], axis=-1)


class JacobiTacker:
    """Bounded optical feedback using the rotating-frame energy gradient.

    The exact CR3BP identity is dC/dt = -2 v_rot dot a. The ephemeris is not
    circular, so C is only a control diagnostic. When C is too small, choose
    the specular normal minimizing sail work; when too large, maximize work.
    The cone is blended from face-on to that optimum as the error grows.
    This heuristic must be judged by propagated trajectories, not the identity.
    """
    def __init__(self, env, target, width=.01, mode="tacking", protected_radius=None):
        self.env, self.target, self.width, self.mode = env, target, width, mode
        self.protected_radius = protected_radius

    def control(self, t, state, rays=None, sample=None):
        if sample is None:
            sample = self.env.at(t)
        if rays is None:
            rays = ray_geometry(state[:3], sample)
        if self.mode == "coast":
            return np.array([np.pi/2, 0.])
        if self.mode == "service_coast":
            # Face-on through the required ray corridor, smoothly feather over
            # 1000 km outside it. The patch supplies no off-duty power here.
            alpha, clock = service_coast_angles(state, rays, sample, self.protected_radius)
        elif self.mode == "face_on":
            alpha, clock = 0., 0.
        else:
            frame, fd = self.env.frames(t), self.env.frames(t, 1)
            local = from_inertial(state, frame, fd)
            distance = length(sample["positions"]["earth"])
            error = self.target-jacobi_value(local, distance)
            velocity = frame@local[3:]
            e, b1, b2 = sail_basis(rays["sun"])
            ve = np.dot(velocity, e)
            v1, v2 = np.dot(velocity, b1), np.dot(velocity, b2)
            vt = np.hypot(v1, v2)
            radical = np.sqrt(9*ve*ve+8*vt*vt)
            minimize = error >= 0
            optimum = np.arctan2(max(radical+(3*ve if minimize else -3*ve), 0.), 4*vt)
            alpha = np.tanh(abs(error)/self.width)*optimum
            clock = np.arctan2(v2, v1)+(np.pi if minimize else 0.)
        angles = np.array([alpha, clock])
        if self.protected_radius is not None:
            margin = length(state[:3])/K.SPEED_OF_LIGHT*length(sample["moon_v"])+50.
            if outgoing_clearance(state[:3], rays, angles, self.protected_radius+margin) <= 0:
                # Reflect radially away from the Moon. This is a passive normal
                # and guarantees positive clearance for these distant orbits.
                e, b1, b2 = sail_basis(rays["sun"])
                delta = e-unit(state[:3])
                if length(delta) < 1e-12:
                    angles = np.array([np.pi/2, 0.])
                else:
                    normal = unit(delta)
                    angles = np.array([np.arccos(np.clip(np.dot(normal, e), 0, 1)),
                                       np.arctan2(np.dot(normal, b2), np.dot(normal, b1))])
        return angles

    def derivative(self, t, state):
        sample = self.env.at(t)
        rays = ray_geometry(state[:3], sample)
        angles = self.control(t, state, rays, sample)
        if abs(np.cos(angles[0])) < 1e-12:
            return np.r_[state[3:], self.env.gravity(t, state[:3])]
        visibility = visible_sun(rays["sun"], rays["earth"], rays["moon"])
        fraction = self.env.reflected_fraction(rays["central"])
        sail, _ = specular_sail(rays["sun"], angles, self.env.sigma, fraction, visibility)
        return np.r_[state[3:], self.env.gravity(t, state[:3])+sail]


def propagate_tacker(env, initial, duration, start=0., width=.01, mode="tacking",
                     step=1800., rtol=2e-9, protected_radius=None, target=None, max_step=None):
    local = from_inertial(initial, env.frames(start), env.frames(start, 1))
    if target is None:
        target = jacobi_value(local, length(env.at(start)["positions"]["earth"]))
    controller = JacobiTacker(env, target, width, mode, protected_radius)
    times = np.linspace(start, start+duration, int(np.ceil(duration/step))+1)
    def collision(t, state):
        # This study requires a screen outside the protected atmosphere.
        return length(state[:3])-(protected_radius or K.MOON_RADIUS)*1.01
    collision.terminal, collision.direction = True, -1
    def escape(t, state):
        return 300e6-length(state[:3])
    escape.terminal, escape.direction = True, -1
    sol = solve_ivp(controller.derivative, [start, start+duration], initial,
                    method="DOP853", t_eval=times, rtol=rtol,
                    atol=np.array([.01]*3+[1e-6]*3), max_step=max_step or step,
                    events=[collision, escape])
    if not sol.success:
        raise RuntimeError(sol.message)
    t, state = sol.t, sol.y.T
    angles = np.array([controller.control(ti, y) for ti, y in zip(t, state)])
    diagnostics = env.acceleration(t, state[:, :3], angles, diagnostics=True)
    local = from_inertial(state, env.frames(t), env.frames(t, 1))
    jacobi = jacobi_value(local, length(env.at(t)["positions"]["earth"]))
    return {"t": t, "state": state, "angles": angles, "local": local,
            "jacobi": jacobi, "jacobi_target": target,
            "completed": sol.status == 0, "mode": mode, **diagnostics}


def cr3bp_derivative(t, state, mu):
    x, y, z, vx, vy, vz = state
    p1, p2 = np.array([x+mu, y, z]), np.array([x-1+mu, y, z])
    gravity = -(1-mu)*p1/length(p1)**3-mu*p2/length(p2)**3
    return [vx, vy, vz, x+2*vy+gravity[0], y-2*vx+gravity[1], gravity[2]]


def dro_seed(radius_m):
    """Symmetry-shoot a planar, lunar distant-retrograde CR3BP seed."""
    mu = K.MOON_GM/(K.EARTH_GM+K.MOON_GM)
    distance = K.EARTH_MOON_DISTANCE
    time_unit = np.sqrt(distance**3/(K.EARTH_GM+K.MOON_GM))
    rho = radius_m/distance
    def crossing(t, y, mu):
        return y[1]
    crossing.direction, crossing.terminal = 1, True
    def half(vy):
        s = np.array([1-mu+rho, 0, 0, 0, vy, 0])
        sol = solve_ivp(cr3bp_derivative, [0, 15], s, args=(mu,),
                        events=crossing, rtol=2e-10, atol=2e-12,
                        max_step=.05, dense_output=False)
        events = sol.t_events[0]
        candidates = [i for i, t in enumerate(events) if t > 1e-5 and sol.y_events[0][i, 4] > 0]
        if not candidates:
            raise ValueError("No retrograde half-orbit crossing")
        i = candidates[0]
        return sol.y_events[0][i, 3], events[i]
    nominal = -np.sqrt(mu/rho)-rho
    grid = np.linspace(1.7*nominal, .45*nominal, 25)
    bracket = None
    last = None
    for vy in grid:
        f, _ = half(vy)
        if last is not None and f*last[1] < 0:
            bracket = (last[0], vy)
            break
        last = vy, f
    if bracket is None:
        raise ValueError("No DRO symmetry bracket")
    vy = brentq(lambda v: half(v)[0], *bracket, xtol=1e-12)
    _, duration = half(vy)
    initial = np.array([radius_m, 0, 0, 0, vy*distance/time_unit, 0])
    return initial, 2*duration*time_unit
