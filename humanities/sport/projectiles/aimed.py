"""Aimed projectiles on the Moon beside Earth: how an arrow, an air-rifle pellet, a small-bore bullet, a paintball and a
clay target fly in the design air, what holding still and taking recoil ask of a shooter, and whether a paintball can
be dodged once it is seen.

    python -m humanities.sport.projectiles.aimed     # writes humanities/sport/projectiles/results/projectiles.json

The sport study's stage for aimed projectiles (8 November 2026). The projectiles fly as points with quadratic drag,
a = -k |v_air| v_air - g z with k = rho C_D A / 2m, in the design run's air at the sea-level datum (state.RHO_MOON)
beside Earth's, integrated by classical Runge-Kutta; the clay target also carries lift as a spinning disc. Shooters
are the sport study's two bodies (state.BODIES); a target dodging a paintball moves under the reach engine's floor law
(reach.profile), starting a reaction time after it sees the shot.

Values and their evidence (sources.json beside this file):
- Paintball: 17.3 mm (.68 calibre) and 3.2 g (the mass assumed, typical), fired at the 91.44 m/s (300 ft/s) limit the
  markers' ASTM standard sets; a smooth sphere's drag below the drag crisis, C_D 0.47 (Reynolds number about 1.2e5
  at the muzzle in either air).
- Arrow: an outdoor target arrow with a 5.5 mm shaft, 19.4 g, at 60 m/s (assumed, typical of a recurve); C_D 3.1 on
  the shaft's section, Park's figure for an arrow yawed 0.75 degrees, with 1.5, his figure for an aligned arrow whose
  boundary layer stays laminar, as the low case.
- Air-rifle pellet: 4.5 mm, 0.53 g, 175 m/s, and a wadcutter's C_D of 0.65 (all assumed, typical of 10 m match
  pellets).
- Small-bore bullet: the .22 LR standard-velocity round, 40 grain (2.59 g), 5.7 mm, at 326 m/s; C_D 0.30, fitted to
  the maker's table (CCI: 1070 ft/s at the muzzle, 977 at 50 yards, 908 at 100); subsonic in both airs.
- Clay target: 110 mm, 105 g, planform 0.0095 m^2, thrown at 30 m/s along its own plane, its attitude held by its
  spin (73 rad/s); C_L = 0.1 + 1.9 alpha (held at its 30-degree value beyond) and C_D = 0.2 + 0.8 C_L^2, the starting
  values of Andert, Freudenthal and Levedag (2016), checked here against their measured skeet throws.
- Heights: a muzzle or release 1.4 m above the floor, a target's chest 1.2 m (the dance body, assumed); the skeet
  houses' launch heights, 3.15 m and 1.42 m, from the same paper.
- Recoil: each projectile's momentum, plus the powder gases' for the firearms (about 0.16 N s for the .22 and 1.6 N s
  for a 24 g trap load fired at 400 m/s; assumed from powder charges near 1200 m/s); the gun's mass moves with the
  shooter; the impulse acts 0.45 m above the centre of mass (a shoulder, assumed).
- Sway: a standing body's intrinsic ankle stiffness is 91% (+-23%) of the gravitational toppling stiffness m g h on
  Earth (Loram and Lakie 2002), too little to stand without active control; the same stiffness is carried to the Moon.
- Dodging: a paintball aimed at the chest misses once the body's centre has moved the torso's half-width (0.20 m for
  the dance body, scaled by stature) plus the ball's radius off its line.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from humanities.dance import envelope as dance
from humanities.sport import envelope as env
from humanities.sport import reach, state
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'projectiles.json'
SCHEMA = 'terluna.humanities.sport-projectiles/1'
FILES = ('humanities/sport/projectiles/aimed.py', 'humanities/sport/state.py', 'humanities/sport/reach.py',
         'humanities/sport/envelope.py', 'humanities/dance/envelope.py', 'humanities/dance/body/model.py',
         'humanities/dance/body/pose.py', 'habitation/low_gravity/movement.py')
INPUTS = (dance.WINDS, dance.POSES)

AIRS = (('moon', state.G_M, state.RHO_MOON), ('earth', state.G_E, state.RHO_EARTH))
MUZZLE_M = 1.4
CHEST_M = 1.2
PROJECTILES = {
    'paintball': dict(kg=0.0032, d_m=0.0173, cd=0.47, m_s=91.44, distances_m=(10, 20, 30, 40, 50, 60)),
    'arrow': dict(kg=0.0194, d_m=0.0055, cd=3.1, cd_low=1.5, m_s=60.0, distances_m=(18, 30, 50, 70, 90)),
    'air-rifle pellet': dict(kg=0.00053, d_m=0.0045, cd=0.65, m_s=175.0, distances_m=(10,)),
    'small-bore bullet': dict(kg=0.002592, d_m=0.0057, cd=0.30, m_s=326.1, distances_m=(25, 50, 100)),
}
CLAY = dict(kg=0.105, d_m=0.11, area_m2=0.0095, cl0=0.1, cla=1.9, cd0=0.2, k=0.8, alpha_max_rad=math.radians(30),
            m_s=30.0)
SKEET = {'high house, long': dict(height_m=3.15, pitch_rad=0.187, measured_m=(61.08, 63.18)),
         'low house, long': dict(height_m=1.42, pitch_rad=0.285, measured_m=(68.22, 69.54))}
SHOOTERS = {                               # projectile kg, m/s, gases N s, gun or bow kg (assumed)
    'paintball marker': (0.0032, 91.44, 0.0, 1.2),
    'recurve bow': (0.0194, 60.0, 0.0, 1.5),
    'air rifle': (0.00053, 175.0, 0.0, 5.0),
    'small-bore rifle': (0.002592, 326.1, 0.16, 5.0),
    'shotgun, 24 g trap load': (0.024, 400.0, 1.6, 3.5),
}
SHOULDER_ABOVE_COM_M = 0.45
ANKLE_STIFFNESS_SHARE = (0.91, 0.68, 1.14)  # of m g h on Earth: Loram and Lakie's mean and one standard deviation
TORSO_HALF_M = 0.20
DODGE_POSTURES = ('ready', 'upright', 'moving, foot down', 'moving, mid-bound', 'aloft')
MID_BOUND_S = 0.25                        # a bounding runner waits this long for its next stance (half a 5 cm bound)
DUCK_M = 0.5                              # sinking the head behind cover that reaches the chest


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16]


def r(x, n=3):
    return None if x is None or not math.isfinite(float(x)) else round(float(x), n)


def drag_k(p, rho, cd=None):
    """k = rho C_D A / 2m per metre: the speed's e-folding distance in a straight flight is 1/k."""
    return rho * (p['cd'] if cd is None else cd) * math.pi * p['d_m'] ** 2 / 8 / p['kg']


# --- flights -----------------------------------------------------------------------------------------------------------

def fly(p, g, rho, elev_deg, h0=MUZZLE_M, wind=0.0, cd=None, dt=None, t_max=120.0, x_max=None):
    """Flights at the projectile's muzzle speed from h0, one per elevation (degrees), in a crosswind `wind` blowing
    along +y: the times t (n,) and x, y, z and speed (m, n), NaN once a flight has reached the floor or passed x_max."""
    el = np.radians(np.atleast_1d(np.asarray(elev_deg, float)))
    k = drag_k(p, rho, cd)
    v0 = p['m_s']
    dt = dt or min(0.01, 0.5 / (k * v0 + 1e-9), 0.2 / v0 * 5)
    s = np.zeros((el.size, 6))
    s[:, 2], s[:, 3], s[:, 5] = h0, v0 * np.cos(el), v0 * np.sin(el)

    def f(q):
        vx, vy, vz = q[:, 3], q[:, 4] - wind, q[:, 5]
        sp = np.sqrt(vx * vx + vy * vy + vz * vz)
        return np.stack([q[:, 3], q[:, 4], q[:, 5], -k * sp * vx, -k * sp * vy, -g - k * sp * vz], 1)

    rows, alive, t = [s.copy()], np.ones(el.size, bool), 0.0
    while alive.any() and t < t_max:
        k1 = f(s)
        k2 = f(s + dt / 2 * k1)
        k3 = f(s + dt / 2 * k2)
        k4 = f(s + dt * k3)
        s = s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        rec = s.copy()
        rec[~alive] = np.nan
        rows.append(rec)
        alive &= (s[:, 2] > 0) & ((s[:, 0] < x_max) if x_max else True)
    S = np.stack(rows, 1)
    speed = np.linalg.norm(S[..., 3:], axis=2)
    return np.arange(S.shape[1]) * dt, S[..., 0], S[..., 1], S[..., 2], speed


def at_distance(t, x, y, z, v, D):
    """Each flight's first crossing of x = D, interpolated: time, y, z and speed (NaN where it never gets there)."""
    m = x.shape[0]
    out = np.full((m, 4), np.nan)
    for i in range(m):
        xi = x[i]
        j = np.flatnonzero(np.nan_to_num(xi, nan=-1.0) >= D)
        if j.size == 0 or j[0] == 0:
            continue
        j = j[0]
        w = (D - xi[j - 1]) / (xi[j] - xi[j - 1])
        out[i] = [t[j - 1] + w * (t[j] - t[j - 1])] + [a[i, j - 1] + w * (a[i, j] - a[i, j - 1]) for a in (y, z, v)]
    return out


def table(p, g, rho, cd=None):
    """At each of the projectile's distances: a horizontal shot's flight time, speed and drop below the bore, and the
    shot aimed to reach chest height there (its elevation, flight time, arrival speed and energy, and the drift a
    crosswind of 1 m/s gives it)."""
    rows = []
    t, x, y, z, v = fly(p, g, rho, [0.0], cd=cd, x_max=max(p['distances_m']) + 1)
    flat = {D: at_distance(t, x, y, z, v, D)[0] for D in p['distances_m']}
    grid = np.linspace(-8.0, 25.0, 1321)
    t, x, y, z, v = fly(p, g, rho, grid, cd=cd, x_max=max(p['distances_m']) + 1)
    for D in p['distances_m']:
        a = at_distance(t, x, y, z, v, D)
        ok = np.isfinite(a[:, 2])
        row = dict(distance_m=D, flat=dict(flight_s=r(flat[D][0]), speed_m_s=r(flat[D][3], 1),
                                           drop_m=r(MUZZLE_M - flat[D][2])))
        zs, els = a[ok, 2], grid[ok]
        hit = None
        for i in range(len(zs) - 1):
            if (zs[i] - CHEST_M) * (zs[i + 1] - CHEST_M) <= 0:
                hit = els[i] + (CHEST_M - zs[i]) / (zs[i + 1] - zs[i]) * (els[i + 1] - els[i])
                break
        if hit is not None:
            tt, xx, yy, zz, vv = fly(p, g, rho, [hit], cd=cd, x_max=D + 1)
            b = at_distance(tt, xx, yy, zz, vv, D)[0]
            tw, xw, yw, zw, vw = fly(p, g, rho, [hit], cd=cd, wind=1.0, x_max=D + 1)
            c = at_distance(tw, xw, yw, zw, vw, D)[0]
            row['aimed'] = dict(elevation_deg=r(hit, 3), flight_s=r(b[0]), speed_m_s=r(b[3], 1),
                                energy_J=r(0.5 * p['kg'] * b[3] ** 2, 2), drift_per_m_s_wind_m=r(c[1]))
        rows.append(row)
    return rows


def farthest(p, g, rho, cd=None):
    """The longest flight from the muzzle's height, over elevations, to the floor: the range a stray shot can carry."""
    els = np.arange(1.0, 60.01, 1.0)
    t, x, y, z, v = fly(p, g, rho, els, cd=cd, dt=0.01, t_max=300.0)
    land = np.array([np.nanmax(xi) for xi in x])
    i = int(np.argmax(land))
    fine = np.arange(max(els[i] - 1, 0.1), els[i] + 1.001, 0.05)
    t, x, y, z, v = fly(p, g, rho, fine, cd=cd, dt=0.005, t_max=300.0)
    land = np.array([np.nanmax(xi) for xi in x])
    j = int(np.argmax(land))
    alive = np.isfinite(x[j])
    return dict(elevation_deg=r(fine[j], 2), range_m=r(land[j], 1), flight_s=r(t[alive][-1], 2),
                apex_m=r(np.nanmax(z[j]), 1))


def projectile(name, p):
    out = dict(kg=p['kg'], diameter_m=p['d_m'], drag_coefficient=p['cd'], muzzle_m_s=p['m_s'])
    for world, g, rho in AIRS:
        k = drag_k(p, rho)
        out[world] = dict(k_per_m=r(k, 6), speed_e_folding_m=r(1 / k, 1), terminal_m_s=r(math.sqrt(g / k), 2),
                          table=table(p, g, rho), farthest=farthest(p, g, rho))
        if 'cd_low' in p:
            out[world]['low_drag'] = dict(drag_coefficient=p['cd_low'], table=table(p, g, rho, cd=p['cd_low']),
                                          farthest=farthest(p, g, rho, cd=p['cd_low']))
    return out


# --- the clay target ----------------------------------------------------------------------------------------------------

def clay_flight(g, rho, height_m, pitch_rad, v0=CLAY['m_s'], dt=0.002, t_max=60.0):
    """A spinning disc thrown along its own plane at `pitch_rad`, its attitude held by its spin, flown in the vertical
    plane with lift across the velocity and drag along it: range, flight time and apex."""
    m, A = CLAY['kg'], CLAY['area_m2']
    x, z, vx, vz, t, apex = 0.0, height_m, v0 * math.cos(pitch_rad), v0 * math.sin(pitch_rad), 0.0, height_m

    def acc(vx, vz):
        gam = math.atan2(vz, vx)
        alpha = min(pitch_rad - gam, CLAY['alpha_max_rad'])
        cl = CLAY['cl0'] + CLAY['cla'] * alpha
        cd = CLAY['cd0'] + CLAY['k'] * cl * cl
        q = 0.5 * rho * (vx * vx + vz * vz) * A / m
        return (-q * cd * math.cos(gam) - q * cl * math.sin(gam), -q * cd * math.sin(gam) + q * cl * math.cos(gam) - g)

    while z > 0 and t < t_max:
        a1 = acc(vx, vz)
        a2 = acc(vx + dt / 2 * a1[0], vz + dt / 2 * a1[1])
        a3 = acc(vx + dt / 2 * a2[0], vz + dt / 2 * a2[1])
        a4 = acc(vx + dt * a3[0], vz + dt * a3[1])
        nvx = vx + dt / 6 * (a1[0] + 2 * a2[0] + 2 * a3[0] + a4[0])
        nvz = vz + dt / 6 * (a1[1] + 2 * a2[1] + 2 * a3[1] + a4[1])
        x, z = x + dt * (vx + nvx) / 2, z + dt * (vz + nvz) / 2
        vx, vz, t = nvx, nvz, t + dt
        apex = max(apex, z)
    return dict(range_m=x, flight_s=t, apex_m=apex)


def clays():
    out = {}
    for name, h in SKEET.items():
        row = dict(measured_on_earth_m=list(h['measured_m']))
        for world, g, rho in AIRS:
            f = clay_flight(g, rho, h['height_m'], h['pitch_rad'])
            row[world] = dict(range_m=r(f['range_m'], 1), flight_s=r(f['flight_s'], 2), apex_m=r(f['apex_m'], 1))
        # the throwing speed that carries a lunar clay as far as the measured Earth throw
        want = sum(h['measured_m']) / 2
        lo, hi = 2.0, CLAY['m_s']
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if clay_flight(state.G_M, state.RHO_MOON, h['height_m'], h['pitch_rad'], v0=mid)[
                'range_m'] < want else (lo, mid)
        f = clay_flight(state.G_M, state.RHO_MOON, h['height_m'], h['pitch_rad'], v0=hi)
        row['moon_for_the_same_range'] = dict(throw_m_s=r(hi, 2), flight_s=r(f['flight_s'], 2), apex_m=r(f['apex_m'], 1))
        out[name] = row
    return out


# --- the shooter --------------------------------------------------------------------------------------------------------

def standing_pitch_inertia(body):
    """The standing body's inertia about its ankles in pitch: the model's whole-body inertia, arms down, about the
    centre of mass, scaled to the body, plus m h^2."""
    I_com = dance.whole_body_inertia()[2][1, 1] * (body.mass_kg / env.MASS_KG) * (body.stature_m / env.STATURE_M) ** 2
    return I_com + body.mass_kg * body.com_height_m ** 2


def recoil(body, g, grip=state.MU):
    """Each shot's push back on the shooter: the body and gun take the impulse's speed, the feet's grip takes it
    away again, and its moment about the centre of mass rocks the body back until the floor's push at the support's
    edge stops it (the throw's balance in envelope.throws, for an impulse)."""
    M = body.mass_kg
    Ia = standing_pitch_inertia(body)
    out = {}
    for name, (kg, v, gas, gun) in SHOOTERS.items():
        J = kg * v + gas
        dv = J / (M + gun)
        L = J * SHOULDER_ABOVE_COM_M
        omega = L / Ia
        alpha = M * g * body.fore_aft_reach_m / Ia
        out[name] = dict(impulse_N_s=r(J, 3), speed_change_m_s=r(dv, 4), grip_stops_it_s=r(dv / (grip * g), 3),
                         rock_back_deg=r(math.degrees(omega * omega / (2 * alpha)), 3), rock_stops_s=r(omega / alpha, 3))
    return out


def sway(body, g):
    """The standing shooter as an inverted pendulum on its ankles: with Earth's intrinsic ankle stiffness K, the net
    stiffness K - m g h decides whether the body stands still without active control, and sets the period it sways
    at when it does."""
    m, h = body.mass_kg, body.com_height_m
    Ia = standing_pitch_inertia(body)
    rows = {}
    for share in ANKLE_STIFFNESS_SHARE:
        K = share * m * state.G_E * h
        net = K - m * g * h
        rows[f'stiffness_{share}_of_earth_mgh'] = dict(
            share_of_toppling=r(K / (m * g * h), 2), passively_stable=net > 0,
            sway_period_s=r(2 * math.pi * math.sqrt(Ia / net), 2) if net > 0 else None)
    return rows


# --- dodging a paintball ------------------------------------------------------------------------------------------------

def aimed_times(p, g, rho, distances):
    """The flight time of the shot aimed to reach chest height at each distance, from one fan of elevations: the
    elevation that brings the fan to chest height is interpolated, and the time with it (None beyond reach)."""
    grid = np.linspace(-8.0, 25.0, 1321)
    t, x, y, z, v = fly(p, g, rho, grid, x_max=max(distances) + 1)
    out = []
    for D in distances:
        a = at_distance(t, x, y, z, v, D)
        ok = np.isfinite(a[:, 2])
        zs, ts = a[ok, 2], a[ok, 0]
        hit = None
        for i in range(len(zs) - 1):
            if (zs[i] - CHEST_M) * (zs[i + 1] - CHEST_M) <= 0:
                w = (CHEST_M - zs[i]) / (zs[i + 1] - zs[i])
                hit = ts[i] + w * (ts[i + 1] - ts[i])
                break
        out.append(hit)
    return out


def paintball_curve(g, rho, far_m=120):
    """The aimed paintball's flight time against distance, every metre from 5 m (None beyond reach)."""
    Ds = list(range(5, far_m + 1))
    return Ds, aimed_times(PROJECTILES['paintball'], g, rho, Ds)


def flight_at(curve, D):
    Ds, ts = curve
    pts = [(d, t) for d, t in zip(Ds, ts) if t is not None]
    if not pts or D > pts[-1][0]:
        return None
    return float(np.interp(D, [d for d, _ in pts], [t for _, t in pts]))


_PROFILES = {}


def escape_m(body, posture, t, grip, g):
    """How far the body's centre can be off the line of a shot it saw, t seconds after the shot left the muzzle: the
    reach engine's profile across the line, the target facing the shooter."""
    move = t - state.REACTION_S
    if posture == 'aloft' or move <= 0:
        return 0.0
    if posture.startswith('moving'):
        st = 'moving'
        nxt = MID_BOUND_S if posture.endswith('mid-bound') else 0.0
    else:
        st, nxt = posture, 0.0
    key = (body.key, st, nxt, grip, g)
    if key not in _PROFILES:
        _PROFILES[key] = reach.profile(body, st, (0.0, 1.0), grip, g, next_stance_s=nxt, facing=0.0)[:2]
    tt, x = _PROFILES[key]
    return float(np.interp(move, tt, x))


def dodge(grip=state.MU):
    """For each world and posture: the escape the body makes before the ball arrives at each distance, and the least
    distance from which a seen shot can be dodged; and how long an aloft player is a sure hit."""
    out = {}
    for world, g, rho in AIRS:
        curve = paintball_curve(g, rho)
        rows = {}
        for key, body in state.BODIES.items():
            need = TORSO_HALF_M * body.stature_m / env.STATURE_M + PROJECTILES['paintball']['d_m'] / 2
            flights = {D: r(flight_at(curve, D)) for D in (10, 20, 30, 40, 50, 60, 80)}
            by = {}
            for posture in DODGE_POSTURES:
                esc = {D: (escape_m(body, posture, tf, grip, g) if tf else None) for D, tf in flights.items()}
                least = None
                if posture != 'aloft':
                    lo, hi = 5.0, 80.0
                    if escape_m(body, posture, flight_at(curve, hi) or 0, grip, g) >= need:
                        for _ in range(25):
                            mid = 0.5 * (lo + hi)
                            tf = flight_at(curve, mid)
                            lo, hi = (mid, hi) if (tf is None or escape_m(body, posture, tf, grip, g) < need) else (lo, mid)
                        least = hi
                by[posture] = dict(escape_m={str(D): r(e) for D, e in esc.items()}, least_dodge_m=r(least, 1))
            rows[key] = dict(need_m=r(need), flight_s={str(D): tf for D, tf in flights.items()}, postures=by)
        jumps = []
        best = state.BODIES['dancer'].take_off_m_s(g) ** 2 / (2 * g)
        for rise in [x for x in (0.3, 1.0, 2.0) if x < best] + [best]:
            air = 2 * math.sqrt(2 * rise / g)
            jumps.append(dict(rise_m=r(rise, 2), airtime_s=r(air, 2),
                              sure_hit_window_s={str(D): (r(max(0.0, air - tf), 2) if air and tf else None)
                                                 for D, tf in rows['dancer']['flight_s'].items()}))
        duck = {f'{DUCK_M}_m': r(env.dip_s(DUCK_M, g), 2)}
        out[world] = dict(by_body=rows, aloft_sure_hit=jumps, duck_by_gravity_s=duck)
    return out


# --- what a game needs --------------------------------------------------------------------------------------------------

def for_games():
    """The constants a game engine needs to fly a paintball and to move a player within these limits."""
    p = PROJECTILES['paintball']
    dancer = state.BODIES['dancer']
    worlds = {}
    for world, g, rho in AIRS:
        worlds[world] = dict(gravity_m_s2=r(g, 4), air_density_kg_m3=r(rho, 4), paintball_k_per_m=r(drag_k(p, rho), 6),
                             floor_push_limit_m_s2={str(mu): r(mu * g, 3) for mu in state.GRIPS},
                             full_take_off_m_s=r(dancer.take_off_m_s(g), 3),
                             crouch_take_off_m_s=r(dancer.take_off_m_s(g, dancer.crouch_push_m), 3))
    return dict(
        paintball=dict(kg=p['kg'], diameter_m=p['d_m'], drag_coefficient=p['cd'], muzzle_m_s=p['m_s'],
                       law='a = -k |v - w| (v - w) - g z, k = rho C_D (pi d^2 / 4) / (2 m), w the wind'),
        player=dict(reaction_s=state.REACTION_S, sprint_profile=dict(F0_N_kg=state.SPRINT_F0_N_KG,
                                                                     V0_m_s=state.SPRINT_V0_M_S),
                    top_speed_m_s=state.V_MAX_M_S, hands_up_m=dancer.reach_up_m, hands_out_m=dancer.reach_out_m,
                    floor_law='horizontal push per kg <= min(mu g, (g/h)(lead + x), F0 (1 - v/V0)) (reach.py)',
                    jump_timing='state.Body.jump_timing: a dip only as deep as the take-off needs'),
        worlds=worlds)


def results():
    return dict(
        schema=SCHEMA,
        producer=dict(domain='humanities', files={f: digest(f) for f in FILES},
                      constants=constants_used([ROOT / f for f in FILES]),
                      inputs={str(p.relative_to(ROOT)): digest(p.relative_to(ROOT)) for p in INPUTS}),
        evidence='Point-mass flights with stated drag (and a disc\'s lift for the clay) in the design run\'s air at '
                 'the shared gravity beside Earth\'s, with specifications from a marker manual, a maker\'s ballistic '
                 'table, arrow and clay-target aerodynamics papers and stated assumptions (sources.json); shooters '
                 'and dodging targets from the sport study\'s bodies and reach engine. Screens, not measurements: '
                 'nothing was fired, thrown or dodged.',
        heights_m=dict(muzzle=MUZZLE_M, chest=CHEST_M),
        projectiles={name: projectile(name, p) for name, p in PROJECTILES.items()},
        clay=clays(),
        shooter={key: dict(recoil={w: recoil(b, g) for w, g, _ in AIRS}, sway={w: sway(b, g) for w, g, _ in AIRS})
                 for key, b in state.BODIES.items()},
        paintball_dodge=dodge(),
        for_games=for_games())


def main() -> int:
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
