"""Moon-centred point-mass dynamics with ephemeris reference acceleration."""
from __future__ import annotations
import numpy as np
from shared import constants as K
from .ephemeris import BODY_GM
from .optical import length, solar_visibility, covering_radius, optical_projection, light_time_allowance, arriving_ray_vectors


def normalized_derivatives(p, v, a):
    r = length(p)
    n = p/r[..., None]
    rd = np.sum(n*v, axis=-1)
    rdd = (np.sum(v*v+p*a, axis=-1)-rd*rd)/r
    nd = (v-n*rd[..., None])/r[..., None]
    ndd = (a-n*rdd[..., None]-2*nd*rd[..., None])/r[..., None]
    return n, nd, ndd


def angular_derivatives(p, v, a):
    """Longitude and two time derivatives after projection to the ecliptic."""
    eps = np.deg2rad(K.EARTH_OBLIQUITY_DEG)
    def xy(z):
        return z[..., 0], z[..., 1]*np.cos(eps)+z[..., 2]*np.sin(eps)
    x, y = xy(p); vx, vy = xy(v); ax, ay = xy(a)
    r2 = x*x+y*y
    cross = x*vy-y*vx
    return np.arctan2(y, x), cross/r2, (x*ay-y*ax)/r2 - cross*2*(x*vx+y*vy)/r2**2


def geometry(samples: dict, harmonics: int = 3):
    u, ud, udd = normalized_derivatives(samples["positions"]["sun"], samples["velocities"]["sun"], samples["sun_relative_a"])
    eps = np.deg2rad(K.EARTH_OBLIQUITY_DEG)
    north = np.array([0, -np.sin(eps), np.cos(eps)])
    v, vd, vdd = normalized_derivatives(np.cross(north, u), np.cross(north, ud), np.cross(north, udd))
    w = np.cross(u, v); wd = np.cross(ud, v)+np.cross(u, vd)
    wdd = np.cross(udd, v)+2*np.cross(ud, vd)+np.cross(u, vdd)
    p, pd, pdd = angular_derivatives(-samples["positions"]["earth"], -samples["velocities"]["earth"], -samples["earth_relative_a"])
    s, sd, sdd = angular_derivatives(samples["positions"]["sun"], samples["velocities"]["sun"], samples["sun_relative_a"])
    phi, phid, phidd = p-s, pd-sd, pdd-sdd
    basis = [np.ones_like(phi)]; bd = [np.zeros_like(phi)]; bdd = [np.zeros_like(phi)]
    for n in range(1, harmonics+1):
        c, sn = np.cos(n*phi), np.sin(n*phi)
        basis.extend([c, sn]); bd.extend([-n*phid*sn, n*phid*c])
        bdd.extend([-n*phidd*sn-(n*phid)**2*c, n*phidd*c-(n*phid)**2*sn])
    return {"frame": np.stack([u, v, w], axis=1),
            "frame_d": np.stack([ud, vd, wd], axis=1),
            "frame_dd": np.stack([udd, vdd, wdd], axis=1),
            "basis": np.stack(basis, axis=1), "basis_d": np.stack(bd, axis=1),
            "basis_dd": np.stack(bdd, axis=1), "phase": phi}


def target(geo: dict, coefficients: np.ndarray):
    """Trajectory coefficients in km in the moving (sunward, transverse) frame."""
    coords = geo["basis"]@coefficients.T*1000
    cd = geo["basis_d"]@coefficients.T*1000
    cdd = geo["basis_dd"]@coefficients.T*1000
    q = np.einsum("ni,nij->nj", coords, geo["frame"])
    qd = np.einsum("ni,nij->nj", cd, geo["frame"])+np.einsum("ni,nij->nj", coords, geo["frame_d"])
    qdd = (np.einsum("ni,nij->nj", cdd, geo["frame"])+
           2*np.einsum("ni,nij->nj", cd, geo["frame_d"])+
           np.einsum("ni,nij->nj", coords, geo["frame_dd"]))
    return {"q": q, "v": qd, "a": qdd, "coords": coords}


def relative_gravity(q, positions, moon_a):
    r = length(q)
    if np.any(r <= K.MOON_RADIUS):
        raise ValueError("Point is inside Moon")
    g = -K.MOON_GM*q/r[..., None]**3 - moon_a
    for name, gm in BODY_GM.items():
        delta = positions[name]-q
        g = g+gm*delta/length(delta)[..., None]**3
    return g


def moon_acceleration_residual(samples):
    g = np.zeros_like(samples["moon_a"])
    for name, gm in BODY_GM.items():
        p = samples["positions"][name]
        g += gm*p/length(p)[..., None]**3
    return g-samples["moon_a"]


def evaluate(samples, geo, coefficients, sigma=0.05, protected_radii=4, fraction=1.0):
    trajectory = target(geo, coefficients)
    q = trajectory["q"]
    required = trajectory["a"]-relative_gravity(q, samples["positions"], samples["moon_a"])
    sun,earth = arriving_ray_vectors(samples["positions"]["sun"]-q,samples["positions"]["earth"]-q,
                                    samples["sun_v"],samples["earth_v"],samples["sun_a"],samples["earth_a"])
    visibility = solar_visibility(sun, earth)
    optical, residual = optical_projection(required, sun, sigma, fraction, visibility)
    d = trajectory["coords"][:, 0]
    offset = length(trajectory["coords"][:, 1:])
    radius = covering_radius(d, length(samples["positions"]["sun"]), protected_radii*K.MOON_RADIUS, offset)
    ray_margin = light_time_allowance(d, samples["moon_v"], samples["sun_v"])
    radius = radius+ray_margin
    return {**trajectory, "required": required, "optical": optical, "residual": residual,
            "visibility": visibility, "covering_radius": radius, "light_time_allowance": ray_margin, "sun": sun,
            "earth_clearance": length(earth)-K.EARTH_RADIUS,
            "sunward_required": np.sum(required*geo["frame"][:, 0], axis=-1)}
