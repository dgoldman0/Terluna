"""The Moon's orientation under the Open Moon's water, at once and through time.

A synchronously rotating Moon keeps its axis of least moment of inertia (A) toward Earth and its axis of greatest
moment (C) along the orbit's normal. Water placed off those axes adds products of inertia that tilt them, and the
Moon follows its new principal axes. The present moments come from lunar laser ranging (beta = (C - A)/B and
gamma = (B - A)/C, as tabulated by Murphy 2013) and C/MR^2 = 0.393112 (Williams et al. 2014). A surface load's
inertia is its degree-2 quadratic form; a geoid change N carries the whole redistribution, the solid Moon's response
included, through Delta I = -(2 M R^2 / 3) Q[N / R]. The spin axis's tilt follows Cassini state 2,
tan eps = sin i / (cos i - x) with x = (3/2)(n / |node rate|)(M_E / (M_E + M_M))((C - A)/C) cos eps, the orbit
average of the potential (3/4) n^2 (C - A) sin^2 eps (Wisdom 2006). The hydrostatic part of today's figure,
J2 = (5/6) k_f q and C22 = (1/4) k_f q, takes k_f from the Darwin-Radau relation
C/MR^2 = (2/3)(1 - (2/5) sqrt(5 / (1 + k_f) - 1)).

Through time the water's degree-2 inertia acts through 1 + k'_2(t) of a viscoelastic Moon (viscoelastic.py), the
seas and lakes as a load applied at once and the pores' water as it arrives, and the rotational and tidal bulge
re-forms about the moving axes through the tidal k_2(t). Both are convolutions of a load history with a relaxation
spectrum, integrated exactly for loads linear within each time step; the axes at the end of a step enter the
bulge's forcing, so each step is iterated to agreement.

The frame of a turned Moon is the rotation matrix whose rows are its axes a (toward Earth), b and c (north) in
today's body frame (x toward 0 N 0 E, z north). Latitude and longitude convert between the frames through it.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from shared.constants import (EARTH_GM, GRAVITATIONAL_CONSTANT, JULIAN_DAY, MOON_GM, MOON_ORBIT_INCLINATION_DEG,
                              MOON_RADIUS, SIDEREAL_MONTH_DAYS)

HERE = Path(__file__).resolve().parent
_SOURCES = {s["key"]: s["values_used"] for s in json.loads((HERE / "interior_sources.json").read_text())["sources"]}
BETA = _SOURCES["moments_llr"]["beta"]
GAMMA = _SOURCES["moments_llr"]["gamma"]
C_MR2 = _SOURCES["moments_llr"]["c_mr2"]
R_C = _SOURCES["moments_llr"]["r_km"] * 1e3       # radius of C/MR^2 and of GRAIL's coefficients
NODE_DAYS = _SOURCES["nodal_period"]["node_days"]
MASS = MOON_GM / GRAVITATIONAL_CONSTANT
EARTH_SHARE = EARTH_GM / (EARTH_GM + MOON_GM)      # M_E / (M_E + M_M)
MR2 = MASS * R_C**2
C_MOMENT = C_MR2 * MR2
B_MOMENT = C_MOMENT * (1 + GAMMA) / (1 + BETA)
A_MOMENT = B_MOMENT - GAMMA * C_MOMENT
BASE = np.diag([A_MOMENT, B_MOMENT, C_MOMENT])
SPIN_Q = (2 * math.pi / (SIDEREAL_MONTH_DAYS * JULIAN_DAY)) ** 2 * R_C**3 / MOON_GM


# ------------------------------------------------------------- loads and moments

def grid_geometry(lat_deg, lon_deg):
    """Exact cell solid angles and unit vectors (3, nlat, nlon) of a regular latitude-longitude grid of centres."""
    lat = np.radians(np.asarray(lat_deg, float))[:, None]
    lon = np.radians(np.asarray(lon_deg, float))[None, :]
    dlat = abs(float(np.radians(lat_deg[1] - lat_deg[0])))
    dlon = abs(float(np.radians(lon_deg[1] - lon_deg[0])))
    weight = (np.sin(np.minimum(lat + dlat / 2, np.pi / 2)) - np.sin(np.maximum(lat - dlat / 2, -np.pi / 2))) * dlon
    weight = weight * np.ones_like(lon)
    units = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat) * np.ones_like(lon)])
    return weight, units


def degree2(field, lat_deg, lon_deg):
    """Degree-1 vector and degree-2 quadratic form of a field on the sphere: f = a . n + n^T Q n + ..."""
    weight, units = grid_geometry(lat_deg, lon_deg)
    a = 3 / (4 * np.pi) * np.einsum("ijk,jk->i", units, field * weight)
    second = np.einsum("ijk,ljk->il", units, units * (field * weight)[None])
    Q = 15 / (8 * np.pi) * (second - np.eye(3) * np.trace(second) / 3)
    return a, Q


def inertia_of_load(sigma, lat_deg, lon_deg, radius=MOON_RADIUS):
    """Traceless inertia change (kg m^2) and mass (kg) of a surface load sigma (kg/m^2) at `radius`."""
    weight, _ = grid_geometry(lat_deg, lon_deg)
    _, Q = degree2(sigma, lat_deg, lon_deg)
    return -radius**4 * (8 * np.pi / 15) * Q, float(np.sum(sigma * weight) * radius**2)


def inertia_of_geoid(N, lat_deg, lon_deg, radius=MOON_RADIUS):
    """Inertia change (kg m^2) carried by a geoid change N (m) on the sphere of `radius`."""
    _, Q = degree2(N / radius, lat_deg, lon_deg)
    return -(2 * MASS * radius**2 / 3) * Q


def principal(dI, base=BASE):
    """Moments A <= B <= C and the frame (rows a, b, c): a toward Earth's side (x > 0), c north."""
    values, vectors = np.linalg.eigh(base + dI)
    a = vectors[:, 0] * np.sign(vectors[0, 0])
    c = vectors[:, 2] * np.sign(vectors[2, 2])
    return values, np.array([a, np.cross(c, a), c])


def fluid_love(c_mr2):
    """Darwin-Radau: C/MR^2 = (2/3)(1 - (2/5) sqrt(5 / (1 + k_f) - 1)), solved for k_f."""
    return 5 / ((2.5 * (1 - 1.5 * c_mr2)) ** 2 + 1) - 1


def bulge_shape(frame):
    """Traceless inertia (kg m^2) of the rotational and tidal bulge per unit Love number, aligned with `frame`:
    J2 = (5/6) q and C22 = q / 4 per unit k."""
    J2, C22 = 5 / 6 * SPIN_Q, SPIN_Q / 4
    local = np.diag([-(J2 / 3 + 2 * C22), -(J2 / 3 - 2 * C22), 2 * J2 / 3]) * MR2
    return frame.T @ local @ frame


def cassini_tilt(c_moment, a_moment):
    """Spin axis's angle to the ecliptic's normal (deg) in Cassini state 2 of a synchronous rotator."""
    i = math.radians(MOON_ORBIT_INCLINATION_DEG)
    k = 1.5 * (NODE_DAYS / SIDEREAL_MONTH_DAYS) * EARTH_SHARE * (c_moment - a_moment) / c_moment
    eps = i
    for _ in range(60):
        eps = math.atan2(math.sin(i), math.cos(i) - k * math.cos(eps))
    return math.degrees(eps) - MOON_ORBIT_INCLINATION_DEG


def describe(dI, base=BASE):
    """The pole's move and its direction, the point beneath Earth, the lock and C - A, the tilt and the free
    libration in longitude of the Moon with inertia base + dI."""
    values, frame = principal(dI, base)
    a, b, c = frame
    A, B, C = values
    return dict(pole_deg=math.degrees(math.acos(min(1.0, c[2]))), pole_toward_lon_deg=math.degrees(math.atan2(c[1], c[0])),
                sub_earth_lat_deg=math.degrees(math.asin(a[2])), sub_earth_lon_deg=math.degrees(math.atan2(a[1], a[0])),
                face_move_deg=math.degrees(math.acos(min(1.0, a[0]))),
                b_minus_a_change=(B - A) / (B_MOMENT - A_MOMENT) - 1, c_minus_a_change=(C - A) / (C_MOMENT - A_MOMENT) - 1,
                tilt_to_ecliptic_deg=cassini_tilt(C, A),
                libration_days=SIDEREAL_MONTH_DAYS / math.sqrt(3 * (B - A) / C), frame=frame)


# ------------------------------------------------------------- frames

def unit(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat) * np.ones_like(lon)], axis=-1)


def to_turned(lat_deg, lon_deg, rotation):
    """Today's latitude and longitude (deg) to the turned Moon's, `rotation` the frame's rows in today's frame."""
    p = unit(np.asarray(lat_deg, float), np.asarray(lon_deg, float)) @ rotation.T
    return np.degrees(np.arcsin(np.clip(p[..., 2], -1, 1))), np.degrees(np.arctan2(p[..., 1], p[..., 0]))


def to_today(lat_deg, lon_deg, rotation):
    """The turned Moon's latitude and longitude (deg) back to today's."""
    p = unit(np.asarray(lat_deg, float), np.asarray(lon_deg, float)) @ rotation
    return np.degrees(np.arcsin(np.clip(p[..., 2], -1, 1))), np.degrees(np.arctan2(p[..., 1], p[..., 0]))


def frame_record(frame):
    """The turned frame as a data product: rotation matrix, poles and the point beneath Earth in today's frame."""
    a, _, c = frame
    lat = lambda v: math.degrees(math.asin(v[2]))
    lon = lambda v: math.degrees(math.atan2(v[1], v[0]))
    return dict(rotation_today_to_turned=frame.tolist(),
                north_pole=dict(lat_deg=lat(c), lon_deg=lon(c)), south_pole=dict(lat_deg=lat(-c), lon_deg=lon(-c)),
                sub_earth_point=dict(lat_deg=lat(a), lon_deg=lon(a)))


# ------------------------------------------------------------- through time

def _relax(state, target_start, target_end, rates, dt):
    """Advance mode states x_j, dx_j/dt = r_j (target(t) - x_j), over dt with the target linear in time."""
    rd = rates * dt
    decay = np.exp(-rd)
    phi = np.where(rd > 1e-12, -np.expm1(-rd) / np.maximum(rd, 1e-300), 1.0 - rd / 2)
    ramp = np.where(rd > 1e-12, 1.0 - phi, rd / 2)         # (dt - (1 - e^-rd)/r) / dt
    slope = target_end - target_start
    return (state * decay[:, None, None] + target_start[None] * (1 - decay)[:, None, None]
            + slope[None] * ramp[:, None, None])


def evolve(times, loads, tide_love, rates, base=BASE, bulge=True, iterations=4):
    """The orientation through time.

    times: years, ascending from 0. loads: pairs (load, love), load(t) a rigid degree-2 inertia (3 x 3) at time t,
    taken linear between the times given, and love its (elastic, amplitudes) of k'_2 on `rates`. tide_love: the same
    for the tidal k_2, through which the rotational and tidal bulge re-forms about the moving axes.
    Returns describe() at each time."""
    k_tide, a_tide = tide_love
    _, frame0 = principal(np.zeros((3, 3)), base)
    shape0 = bulge_shape(frame0)
    histories = [np.array([np.asarray(load(t), float) for t in times]) for load, _ in loads]
    states = [np.zeros((len(rates), 3, 3)) for _ in loads]
    x_tide = np.zeros((len(rates), 3, 3))
    rows = []
    prev_shape = shape0
    for i, t in enumerate(times):
        dt = t - times[i - 1] if i else 0.0
        dI_load = np.zeros((3, 3))
        for j, ((_, (k_e, amps)), history) in enumerate(zip(loads, histories)):
            if i:
                states[j] = _relax(states[j], history[i - 1], history[i], rates, dt)
            dI_load = dI_load + (1 + k_e) * history[i] + np.einsum("j,jkl->kl", amps, states[j])
        shape = prev_shape
        x_t = x_tide
        for _ in range(iterations):
            dI = dI_load
            if bulge:
                x_t = _relax(x_tide, prev_shape - shape0, shape - shape0, rates, dt) if i else x_tide
                dI = dI_load + k_tide * (shape - shape0) + np.einsum("j,jkl->kl", a_tide, x_t)
            _, frame = principal(dI, base)
            shape = bulge_shape(frame)
        x_tide, prev_shape = x_t, shape
        row = describe(dI, base)
        row["years"] = float(t)
        rows.append(row)
    return rows
