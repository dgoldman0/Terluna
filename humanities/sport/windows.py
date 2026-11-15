"""The ball windows: what a ball allows. How fast it can be thrown from a stance, where it goes, how long it takes, and
what can catch it, on the Moon in the design run's air beside Earth.

    from humanities.sport import windows
    from humanities.sport.state import BODIES, G_M, RHO_MOON, State
    fan = windows.passes(BODIES['dancer'], State(posture='split'), 'split', 'handball (size 3)', (15.0, 0.0), G_M,
                         RHO_MOON, grip=0.8)

The contests weigh a ball's flight against the reach engine (reach.py): a pass is open when a receiver's hands can be
on its path before a defender's. This module gives the flights and the limits a body puts on them.

The model and its choices:
- Flights. A ball is a point with quadratic drag, a = -k |v| v - g z^, with k = rho C_D A / 2m (envelope.drag_area:
  a smooth sphere below the drag crisis, the shuttlecock from its terminal speed), in still air, integrated by
  classical Runge-Kutta. A flight ends when the ball's centre reaches the floor (z = 0; its radius is ignored), found
  by interpolating within the last step.
- Releases. A stance holds a throw while the floor can carry the throw's push without turning or sliding the thrower
  (envelope.throws): the hand pushes the ball to its release speed over envelope.THROW_S, so the hand force is
  F = m v / THROW_S; it acts at the lever, the release's height above the centre of mass, and the floor's push N can
  sit as far behind the centre of mass as the stance's reach, so F lever <= N reach; the feet hold F <= mu N; and the
  arm's own fastest release, RELEASE_MAX_M_S, caps both. The lever defaults to Body.release_above_com_m, the
  envelope's overhand release (1.0 m for the dance body); envelope.AIR_LEVERS_M gives 0.3 m for a chest pass, whose
  short lever lets the grip and the arm bind instead. The stances: 'feet_under' (the reach is fore_aft_reach, N the
  weight), 'split' (split_reach, the weight), 'landing' (fore_aft_reach, N the landing's mean push for a landing at
  landing_v_down over envelope.STROKE_M) and 'aloft' (no floor: the arm alone, and the throw's angular momentum
  m v lever, its line taken horizontal where it is largest, turns the body until the landing,
  envelope.attitude_error_deg).
- Catches. The same balance for a ball caught at speed v over the hand's give (state.CATCH_S, as a throw gives it;
  a longer give takes a smaller push). The arm's limit is taken to hold for a catch as for a throw (assumed).
- Passes. The thrower's own velocity, and its vertical speed when aloft, add to the release; the stance limits the
  release speed relative to the hand. The ball is aimed along the floor plan at the target; a moving thrower's ball
  drifts off that line, and the fan reports the drift.
"""
from __future__ import annotations

import math

import numpy as np

from humanities.sport import envelope, reach
from humanities.sport.state import CATCH_S, MU, RELEASE_MAX_M_S

STANCES = ('feet_under', 'split', 'landing', 'aloft')


def drag_k(ball, rho):
    """The drag constant per unit mass, k = rho C_D A / 2m, for a ball named in envelope.BALLS."""
    return 0.5 * rho * envelope.drag_area(ball) / envelope.BALLS[ball][0]


def _acc(v, k, g):
    a = -k * np.sqrt((v * v).sum(-1, keepdims=True)) * v
    a[..., 2] -= g
    return a


def _step(p, v, k, g, dt):
    a1 = _acc(v, k, g)
    a2 = _acc(v + 0.5 * dt * a1, k, g)
    a3 = _acc(v + 0.5 * dt * a2, k, g)
    a4 = _acc(v + dt * a3, k, g)
    return (p + dt * v + dt * dt / 6 * (a1 + a2 + a3), v + dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4))


def fly(ball, p0, v0, g, rho, dt=0.004, t_max=30.0):
    """Flights from p0 ((3,) or (m, 3)) at the initial velocities v0 ((m, 3)) until each reaches the floor or t_max:
    the times t (n,), positions p and velocities v (m, n, 3; NaN after landing), and the landing times land_s, places
    land_p and speeds land_m_s (NaN for a flight still up at t_max). Memory grows as m times the longest flight."""
    v = np.array(v0, float).reshape(-1, 3)
    m = v.shape[0]
    p = np.broadcast_to(np.asarray(p0, float), (m, 3)).astype(float)
    k = drag_k(ball, rho)
    land_s, land_v = np.full(m, np.nan), np.full(m, np.nan)
    land_p = np.full((m, 3), np.nan)
    up = (p[:, 2] > 0) | ((p[:, 2] == 0) & (v[:, 2] > 0))
    P, V = [p.copy()], [v.copy()]
    t = 0.0
    with np.errstate(invalid='ignore'):
        while up.any() and t < t_max - 1e-9:
            pn, vn = _step(p, v, k, g, dt)
            hit = up & (pn[:, 2] <= 0)
            if hit.any():
                f = p[hit, 2] / (p[hit, 2] - pn[hit, 2])
                land_s[hit] = t + f * dt
                land_p[hit] = p[hit] + f[:, None] * (pn[hit] - p[hit])
                land_v[hit] = np.linalg.norm(v[hit] + f[:, None] * (vn[hit] - v[hit]), axis=1)
                up &= ~hit
                pn[~up], vn[~up] = np.nan, np.nan
            p, v, t = pn, vn, t + dt
            P.append(p.copy())
            V.append(v.copy())
    return dict(t=np.arange(len(P)) * dt, p=np.stack(P, 1), v=np.stack(V, 1), land_s=land_s, land_p=land_p,
                land_m_s=land_v)


def _support(body, stance, g, landing_v_down):
    """The floor's push N and how far behind the centre of mass it can sit, in a stance."""
    if stance == 'feet_under':
        return body.mass_kg * g, body.fore_aft_reach_m
    if stance == 'split':
        return body.mass_kg * g, body.split_reach_m
    if stance == 'landing':
        return body.mass_kg * (landing_v_down ** 2 / (2 * envelope.STROKE_M) + g), body.fore_aft_reach_m
    raise ValueError(f'unknown stance {stance!r}; one of {STANCES}')


def _limit(body, stance, ball, g, grip, push_s, landing_v_down, lever, descent_s):
    m = envelope.BALLS[ball][0]
    lever = body.release_above_com_m if lever is None else float(lever)
    if stance == 'aloft':
        v = RELEASE_MAX_M_S
        dL = m * v * lever
        turn = {k: (envelope.attitude_error_deg(dL, body.inertia_kg_m2[k], g) if descent_s is None
                    else math.degrees(dL / body.inertia_kg_m2[k] * descent_s)) for k in ('layout', 'tuck')}
        return dict(stance=stance, m_s=v, limited_by='arm', lever_m=lever, angular_momentum_N_m_s=dL,
                    turn_by_landing_deg=turn, thrower_speed_change_m_s=m * v / body.mass_kg)
    N, reach_m = _support(body, stance, g, landing_v_down)
    lims = dict(turning=N * reach_m * push_s / (m * lever), sliding=grip * N * push_s / m, arm=RELEASE_MAX_M_S)
    by = min(lims, key=lims.get)
    return dict(stance=stance, m_s=lims[by], limited_by=by, lever_m=lever, turning_m_s=lims['turning'],
                sliding_m_s=lims['sliding'], arm_m_s=RELEASE_MAX_M_S, floor_push_N=N, reach_m=reach_m,
                hand_force_N=m * lims[by] / push_s)


def release_limit(body, stance, ball, g, grip, landing_v_down=2.5, lever=None, descent_s=None):
    """The fastest release a stance holds without turning or sliding the thrower (the module's balance), as a dict:
    m_s, what limits it ('turning', 'sliding' or 'arm') and the terms. Aloft, the arm's limit with the angular
    momentum the throw leaves and the turn it makes by the landing (layout and tuck; the descent of envelope's apex
    jump unless descent_s is given)."""
    return _limit(body, stance, ball, g, grip, envelope.THROW_S, landing_v_down, lever, descent_s)


def catch_limit(body, stance, ball, g, grip, give_s=CATCH_S, landing_v_down=2.5, lever=None, descent_s=None):
    """The fastest ball a stance can catch without turning or sliding the catcher: the same balance with the ball's
    momentum taken over the give, F = m v / give_s, at the lever (the catch's height above the centre of mass)."""
    return _limit(body, stance, ball, g, grip, give_s, landing_v_down, lever, descent_s)


def release_height(body, thrower, lever=None):
    """The release's height: the centre of mass's (the posture's on the floor, the arc's when aloft) plus the lever."""
    lever = body.release_above_com_m if lever is None else float(lever)
    com = thrower.z if thrower.posture == 'aloft' and thrower.z is not None else body.com_m(thrower.posture)
    return com + lever


def passes(body, thrower, stance, ball, target, g, rho, grip, speeds=None, elevations=None, lever=None,
           sample_dt=0.05, dt=0.004):
    """The fan of flights from the thrower's release toward a target on the floor plan. speeds are release speeds
    relative to the hand (default nine up to the stance's limit), elevations in degrees (default -10 to 75 by 5).
    For each (speed, elevation): when the ball passes the target's distance along the line (pass_s, NaN if it lands
    short), its height and speed there and its sideways drift off the line; where it lands; and the path sampled every
    sample_dt (path (S, E, n, 3), NaN after landing, at times t) for interception checks with reach.reach_along."""
    lim = release_limit(body, stance, ball, g, grip, lever=lever)['m_s']
    speeds = lim * np.linspace(0.2, 1.0, 9) if speeds is None else np.atleast_1d(np.asarray(speeds, float))
    elev = np.arange(-10.0, 76.0, 5.0) if elevations is None else np.atleast_1d(np.asarray(elevations, float))
    px, py = float(thrower.p[0]), float(thrower.p[1])
    z0 = release_height(body, thrower, lever)
    L = math.hypot(target[0] - px, target[1] - py)
    ux, uy = (target[0] - px) / L, (target[1] - py) / L
    S, E = np.meshgrid(speeds, np.radians(elev), indexing='ij')
    vz0 = thrower.vz if thrower.posture == 'aloft' else 0.0
    v0 = np.stack([S * np.cos(E) * ux + thrower.v[0], S * np.cos(E) * uy + thrower.v[1], S * np.sin(E) + vz0], -1)
    f = fly(ball, (px, py, z0), v0.reshape(-1, 3), g, rho, dt=dt)
    rel = f['p'] - np.array([px, py, 0.0])
    along = rel[..., 0] * ux + rel[..., 1] * uy
    with np.errstate(invalid='ignore'):
        past = along >= L
    i = np.argmax(past, 1)
    has = past.any(1) & (i > 0)
    r = np.arange(len(i))
    j = np.maximum(i - 1, 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        w = np.where(has, (L - along[r, j]) / (along[r, i] - along[r, j]), np.nan)
    at = f['p'][r, j] + w[:, None] * (f['p'][r, i] - f['p'][r, j])
    vel = f['v'][r, j] + w[:, None] * (f['v'][r, i] - f['v'][r, j])
    step = max(1, int(round(sample_dt / dt)))
    shape = S.shape
    land = f['land_p'] - np.array([px, py, 0.0])
    return dict(speeds_m_s=speeds, elevations_deg=elev, limit_m_s=lim, held=speeds <= lim + 1e-9,
                release_m=(px, py, z0), distance_m=L,
                pass_s=np.where(has, (j + w) * dt, np.nan).reshape(shape), height_m=at[:, 2].reshape(shape),
                speed_m_s=np.linalg.norm(vel, axis=1).reshape(shape),
                drift_m=(-(at[:, 0] - px) * uy + (at[:, 1] - py) * ux).reshape(shape),
                land_s=f['land_s'].reshape(shape), land_m=(land[:, 0] * ux + land[:, 1] * uy).reshape(shape),
                t=f['t'][::step], path=f['p'][:, ::step].reshape(shape + (-1, 3)))


def _shoot(ball, p0, V, T, g, rho, n=400):
    """Positions and velocities at exactly T for initial velocities V (m, 3), with no floor; and the lowest height."""
    k = drag_k(ball, rho)
    p = np.broadcast_to(np.asarray(p0, float), V.shape).astype(float)
    v = np.array(V, float)
    low = p[:, 2].copy()
    for _ in range(n):
        p, v = _step(p, v, k, g, T / n)
        low = np.minimum(low, p[:, 2])
    return p, v, low


def lob(ball, p0, target, flight_s, g, rho, tol=1e-4, iterations=20):
    """The release velocity that brings the ball from p0 to the point `target` (x, y, z) in flight_s, with drag:
    Newton's method on the shot, from the vacuum's answer, with the Jacobian by finite differences. Returns the
    velocity, its speed and elevation, the arrival velocity and speed, the miss and the lowest height on the way (a
    negative one means the shot passes through the floor)."""
    p0, tgt = np.asarray(p0, float), np.asarray(target, float)
    T = float(flight_s)
    v = (tgt - p0) / T + np.array([0.0, 0.0, 0.5 * g * T])
    for _ in range(iterations):
        eps = 1e-5 * max(1.0, float(np.linalg.norm(v)))
        P, Vel, low = _shoot(ball, p0, np.vstack([v, v + eps * np.eye(3)]), T, g, rho)
        miss = P[0] - tgt
        if np.linalg.norm(miss) < tol:
            break
        v = v - np.linalg.solve(((P[1:] - P[0]) / eps).T, miss)
    P, Vel, low = _shoot(ball, p0, v[None, :], T, g, rho)
    return dict(v0=v, speed_m_s=float(np.linalg.norm(v)),
                elevation_deg=math.degrees(math.atan2(v[2], math.hypot(*v[:2]))), arrives=Vel[0],
                arrives_m_s=float(np.linalg.norm(Vel[0])), miss_m=float(np.linalg.norm(P[0] - tgt)),
                lowest_m=float(low[0]))


def catch_band(body, posture, g):
    """The heights the hands can catch at: standing (Body.hand_band at the standing height), in the posture itself, and
    with the posture's best jump (its full take-off), with that jump's delay to take-off and to the apex."""
    best = reach.jumps(body, posture, g)
    return dict(standing_m=body.hand_band(body.com_height_m), in_posture_m=body.hand_band(body.com_m(posture)),
                jump_top_m=best['hands_top_m'], take_off_s=best['take_off_s'], apex_s=best['apex_s'], jump=best)


def stance_figures(body, ball, g, grip=MU, lever=None):
    """The fastest release and catch by stance, as one table."""
    return {s: dict(release_m_s=release_limit(body, s, ball, g, grip, lever=lever)['m_s'],
                    catch_m_s=catch_limit(body, s, ball, g, grip, lever=lever)['m_s']) for s in STANCES}
