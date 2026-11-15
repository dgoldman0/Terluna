"""The reach engine: where a player's hands can be, and when, on a floor at lunar gravity beside Earth.

    from humanities.sport import reach
    from humanities.sport.state import BODIES, G_M, State
    r = reach.reach_time(BODIES['dancer'], State(posture='ready'), (0.0, 3.0, 1.2), grip=0.8, g=G_M)

The second round of the sport study asks who can be where, and when, against the time a ball takes to get there. The
engine answers for one player from a State (state.py): the centre of mass moves over the floor under the floor's push,
and the hands reach the band Body.hand_band around it, out to Body.reach_out_m horizontally (a cylinder: the arms'
sweep inside it is not modelled). Times run from the moment the player starts to act; the contests add REACTION_S to a
cue themselves. Grip and gravity are arguments everywhere, so the contests sweep both.

The model and its choices:
1. On the floor, per direction. Along a unit direction u in the drift frame (the place p + v t the player would drift
   to), x_u(t) is the most distance the centre of mass can cover by t, from rest in that frame.
   - Lean, then grip: x'' = min(mu g, (g/h)(lead0(u) + x), muscle_push(s)). The first term is the floor's grip, the
     second the inverted pendulum about the support's far edge (envelope.lean_start), the third the legs' sprint
     profile at the speed s = v.u + x' along u (state.muscle_push). h is the posture's centre-of-mass height
     (Body.com_m) and lead0(u) the support's farthest extent behind the centre of mass along u (lead0). On the Moon
     the grip binds below about 6-7 m/s; on Earth the muscles bind from the start at grips of 0.8 and above.
   - A drop step first, where lead0(u) is short of split_reach: no push for STEP_S while a foot goes split_reach
     behind, then the same law from there at the ready height. Whichever covers more at t counts. From upright this
     takes the body down to the ready height within the step, which on the Moon is quicker than gravity allows (a
     0.16 m sink and stop takes about 0.5 s), so the upright drop step is flattered; at its own height it would reach
     the first 0.1 m sideways about 0.04 s later. Even so it beats the lean alone there on the Moon.
   - Moving: the feet are placed as the push needs, x'' = min(mu g, muscle_push(s)) from next_stance_s on (a stride's
     average push).
   - The speed along u, v.u + x', stays at or below V_MAX_M_S.
   The law has no closed form once the muscles join the minimum, so each profile is integrated (classical Runge-Kutta
   at DT, exact for the grip's constant push) and read between steps by cubic Hermite interpolation. The reach set at
   t is the drift point plus x_u(t) along every u. This assumes the player heads for the point from the start, which
   is the double integrator's optimum when the push is the same in every direction: exact for 'moving' while the grip
   binds and the cap does not, and an inner estimate where the support, the cap or the muscles make the push depend
   on direction. A point is in reach when it lies within x_u(t) + reach_out of the drift point along its own ray.
2. Jumps. The take-off follows the contract's timing (Body.jump_timing): a jump sinks only as deep as its take-off
   speed vz needs (at least MIN_DIP_M from upright, none from the ready crouch or split stance up to the crouch's best
   take-off), then pushes; a moving player rebounds at the next stance in one stride's stance. The legs leave the
   floor at the take-off height, push_off - dip above standing, and land when the centre of mass is back at the
   standing height (State.landing_s), so an arc returns to its take-off height after 2 vz/g and lands a little later.
   The push adds horizontal speed up to k vz along u, with k = min(mu, lead0(u)/h) on the floor postures
   (envelope.jumps' lateral hop) and mu when moving. The player may travel on the floor first: the take-off comes at
   t_f + dip + push, with the floor's speed at t_f held through the dip and push (assumed; the push's own horizontal
   build-up over its stroke is left out, a slight understatement). Aloft, the arc is fixed and the hands ride
   Body.hand_band(z). The search runs over the time in the air and vz on small grids (VZ_SHARES of the full take-off,
   the crouch's best take-off without a dip, and per point the slowest take-off that reaches its height), vectorised.
3. Aloft and landing. An aloft player drifts until State.landing_s. The landing's stance, 2 STROKE_M / v_down
   (envelope.braking_landing), can change the travel by up to mu (v_down + g T) in any direction, at a steady push
   over the stance (the muscles' cap is not applied to it: the legs take it as a strut); from then on the player is
   moving, and may rebound out of the landing.
4. Down. A fallen player gets up first (envelope.floor_and_falls' stand-up time from 0.45 m) and is then ready. Until
   up, the hands reach from the floor to the standing band's low end within reach_out (assumed).
5. Heights. On the floor the hands reach the standing band up to reach_up (rising out of the crouch is taken as
   free); below the posture's band they reach a point by sinking, no faster than envelope.dip_s allows from the
   start (the scoop; the sink is assumed to run alongside the travel).

Margins are signed distances to spare: the smaller of the horizontal spare along the point's ray and the vertical spare
inside the band. Where both the floor and a jump reach a point, the engine answers with the floor (no commitment). The
jump search's grids make a jump's margin a slight underestimate (a few centimetres at most where the window is
narrow); reach_along reads the floor profiles from a direction-time table (N_DIR directions, DT steps) built once per
state, and reach_at and reach_time use the same table, so all three agree.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from humanities.dance import envelope as dance
from humanities.sport import envelope
from humanities.sport.state import (MIN_DIP_M, MU, SPRINT_F0_N_KG, SPRINT_V0_M_S, STEP_S, V_MAX_M_S, Reach, State)

MODES = ('ground', 'jump', 'aloft', 'landing')    # reach_along's mode codes index this
DT = 0.01                                         # the profiles' integration and table step, s
SPAN_S = 12.0                                     # the table's span; beyond it the last speed holds
N_DIR = 120                                       # the table's direction bins (3 degrees)
VZ_SHARES = (0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 1.0)  # take-off speeds searched, as shares of the full take-off
N_TAU = 14                                        # times in the air searched, plus the apex and the band's centring
GET_UP_FROM_M = 0.45                              # the stand-up the envelope times (floor_and_falls)
V_DOWN_MIN = 0.05                                 # a landing's least downward speed, so its stance stays finite
BIG = 1e3


# --- the support, the push and the profiles ------------------------------------------------------------------------

def lead0(body, posture, ux, uy, facing=0.0):
    """The support's farthest extent behind the centre of mass along (ux, uy): where the floor's push can sit to lean
    the body that way. The support is a rectangle in the body's frame (forward = facing): upright ±fore_aft_reach by
    ±side_reach, ready ±fore_aft_reach by ±ready_half_stance, split from -split_reach to +fore_aft_reach forward by
    ±side_reach. 'down' and 'moving' read as ready (the stance a player gets up into; a moving jump's k is mu)."""
    c, s = math.cos(facing), math.sin(facing)
    uf, us = ux * c + uy * s, -ux * s + uy * c
    back, front, side = body.fore_aft_reach_m, body.fore_aft_reach_m, body.ready_half_stance_m
    if posture == 'upright':
        side = body.side_reach_m
    elif posture == 'split':
        back, side = body.split_reach_m, body.side_reach_m
    return np.where(uf >= 0, back * uf, -front * uf) + side * np.abs(us)


def muscle(v_along):
    """state.muscle_push vectorised: the legs' horizontal push per kilogram at a speed along the push."""
    return SPRINT_F0_N_KG * np.maximum(0.0, 1.0 - np.asarray(v_along, float) / SPRINT_V0_M_S)


def _run(lead, h, mu, g, v_along, n, x0=0.0, v0=0.0, dt=DT):
    """x and x' on s = 0, dt, ..., n dt under x'' = min(mu g, (g/h)(lead + x), muscle(v_along + x')), the speed held at
    or below V_MAX_M_S - v_along: classical Runge-Kutta on (x, x'), exact while one constant push binds."""
    lead, va = np.broadcast_arrays(np.asarray(lead, float), np.asarray(v_along, float))
    cap = np.clip(V_MAX_M_S - va, 0.0, None)
    x = np.zeros(lead.shape) + x0
    v = np.minimum(np.zeros(lead.shape) + v0, cap)
    X, V = np.empty(lead.shape + (n + 1,)), np.empty(lead.shape + (n + 1,))
    X[..., 0], V[..., 0] = x, v
    k = g / h

    def acc(x, v):
        a = np.minimum(np.minimum(mu * g, k * (lead + x)), muscle(va + v))
        return np.where(v >= cap, 0.0, a)

    for i in range(n):
        a1 = acc(x, v)
        a2 = acc(x + 0.5 * dt * v, v + 0.5 * dt * a1)
        a3 = acc(x + 0.5 * dt * v + 0.25 * dt * dt * a1, v + 0.5 * dt * a2)
        a4 = acc(x + dt * v + 0.5 * dt * dt * a2, v + dt * a3)
        x = x + dt * v + dt * dt / 6 * (a1 + a2 + a3)
        v = np.minimum(v + dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4), cap)
        X[..., i + 1], V[..., i + 1] = x, v
    return X, V


def _hermite(X, V, s, dt=DT):
    """X and V (rows on the step dt from 0) read at the times s (a 1-D array): cubic Hermite for x, linear for x';
    before 0 the row's start, beyond the end its last speed held."""
    n = X.shape[-1]
    f = np.maximum(s, 0.0) / dt
    i = np.minimum(np.floor(f), n - 2).astype(np.intp)
    u = f - i
    w = np.minimum(u, 1.0)
    x0, x1, v0, v1 = X[..., i], X[..., i + 1], V[..., i], V[..., i + 1]
    x = ((1 + 2 * w) * (1 - w) ** 2 * x0 + w * (1 - w) ** 2 * dt * v0 + w * w * (3 - 2 * w) * x1
         + w * w * (w - 1) * dt * v1 + (u - w) * dt * v1)
    return x, v0 + (v1 - v0) * w


def get_up_s(body, g):
    """A fallen player's time to stand: envelope.floor_and_falls' rise of 0.45 m that ends at rest."""
    v = body.take_off_m_s(g)
    return envelope.rise_and_stop_s(GET_UP_FROM_M, v * v / (2 * body.push_off_m), g)


def sink_m(t, g):
    """How far the centre of mass can sink and stop within t from rest: envelope.dip_s inverted (it grows as the square
    root of the depth)."""
    return (np.maximum(np.asarray(t, float), 0.0) / envelope.dip_s(1.0, g)) ** 2


def landing(body, state, grip, g):
    """An aloft player's landing: when the feet touch (State.landing_s), the downward speed then, the stance's length
    2 STROKE_M / v_down and the change of travel it allows, mu (v_down + g T) (envelope.braking_landing)."""
    T = state.landing_s(body, g)
    v_down = max(g * T - state.vz, V_DOWN_MIN)
    stance = 2 * envelope.STROKE_M / v_down
    return dict(landing_s=T, v_down_m_s=v_down, stance_s=stance, change_m_s=grip * (v_down + g * stance))


def _profiles(body, st, ux, uy, t, grip, g, steps=True):
    """x_u(t) and x_u'(t) for the state along each direction (ux, uy) (1-D arrays) at the times t (1-D): arrays of shape
    (directions, times)."""
    ux, uy = np.atleast_1d(np.asarray(ux, float)), np.atleast_1d(np.asarray(uy, float))
    t = np.atleast_1d(np.asarray(t, float))
    va = float(st.v[0]) * ux + float(st.v[1]) * uy
    X, V = np.zeros((ux.size, t.size)), np.zeros((ux.size, t.size))
    tmax = float(t.max()) if t.size else 0.0
    p = st.posture
    if p == 'aloft':
        L = landing(body, st, grip, g)
        cap = np.clip(V_MAX_M_S - va, 0.0, None)[:, None]
        a_s = L['change_m_s'] / L['stance_s']
        s = np.clip(t - L['landing_s'], 0.0, L['stance_s'])[None, :]
        sc = cap / a_s                                       # the cap within the stance
        X = np.where(s > sc, 0.5 * a_s * sc * sc + cap * (s - sc), 0.5 * a_s * s * s)
        V = np.minimum(a_s * s, cap)
        t1 = L['landing_s'] + L['stance_s']
        if tmax > t1:
            n = int(math.ceil((tmax - t1) / DT)) + 1
            xs, vs = _run(BIG, body.com_height_m, grip, g, va, n, _end(a_s, L['stance_s'], cap),
                          np.minimum(L['change_m_s'], cap[:, 0]))
            xa, vb = _hermite(xs, vs, t - t1)
            after = (t > t1)[None, :]
            X, V = np.where(after, xa, X), np.where(after, vb, V)
        return X, V
    if p == 'moving':
        plans = [(st.next_stance_s, BIG, body.com_height_m, None)]
    else:
        t0 = get_up_s(body, g) if p == 'down' else 0.0
        q = 'ready' if p == 'down' else p
        lead = lead0(body, q, ux, uy, st.facing)
        plans = [(t0, lead, body.com_m(q), None)]
        if steps:
            plans.append((t0 + STEP_S, body.split_reach_m, body.ready_com_m, lead < body.split_reach_m - 1e-12))
    for t0, lead, h, mask in plans:
        if tmax <= t0:
            continue
        n = int(math.ceil((tmax - t0) / DT)) + 1
        xs, vs = _run(np.broadcast_to(lead, ux.shape), h, grip, g, va, n)
        xi, vi = _hermite(xs, vs, t - t0)
        live = (t > t0)[None, :] & (xi > X)
        if mask is not None:
            live &= mask[:, None]
        X, V = np.where(live, xi, X), np.where(live, vi, V)
    return X, V


def _end(a_s, stance, cap):
    """The distance a landing's stance covers at its steady push, its speed capped."""
    sc = cap[:, 0] / a_s
    return np.where(stance > sc, 0.5 * a_s * sc * sc + cap[:, 0] * (stance - sc), 0.5 * a_s * stance * stance)


def profile(body, posture, u, grip, g, v_along=0.0, next_stance_s=0.0, facing=0.0, t=None, steps=True, z=None,
            vz=0.0):
    """Arrays t, x and x' for one direction u (an angle in radians, or a 2-vector): the most distance the centre of mass
    can cover along u in the drift frame by each t, and its speed there relative to the drift. v_along is the drift
    velocity's component along u (it sets the speed cap and the muscles' speed); steps=False leaves the drop step out
    (the lean alone); z and vz give the arc for 'aloft'."""
    if np.ndim(u) == 0:
        ux, uy = math.cos(u), math.sin(u)
    else:
        ux, uy = np.asarray(u, float) / np.hypot(*u)
    t = np.arange(0.0, 10.0 + DT / 2, DT) if t is None else np.asarray(t, float)
    st = State(v=(v_along * ux, v_along * uy), posture=posture, facing=facing, z=z, vz=vz, next_stance_s=next_stance_s)
    x, v = _profiles(body, st, ux, uy, t, grip, g, steps=steps)
    return t, x[0], v[0]


def reach_set(body, state, t, grip, g, n=72):
    """An (n, 2) boundary of where the centre of mass can be at time t: the drift point plus x_u(t) along n
    directions."""
    th = np.arange(n) * (2 * math.pi / n)
    x, _ = _profiles(body, state, np.cos(th), np.sin(th), [float(t)], grip, g)
    px, py = state.p[0] + state.v[0] * t, state.p[1] + state.v[1] * t
    return np.stack([px + x[:, 0] * np.cos(th), py + x[:, 0] * np.sin(th)], -1)


# --- jumps -----------------------------------------------------------------------------------------------------------

def take_off_height(body):
    """The centre of mass's height as the feet leave the floor, on the toes: push_off - dip above standing."""
    return body.com_height_m + body.push_off_m - body.dip_m


def airtime(body, vz, g):
    """From take-off at vz until the centre of mass is back at the standing height (State.landing_s's landing)."""
    dz = take_off_height(body) - body.com_height_m
    return (vz + np.sqrt(np.asarray(vz, float) ** 2 + 2 * g * dz)) / g


def prep_s(body, posture, vz, g):
    """Body.jump_timing vectorised over vz: the dip's seconds plus the push's seconds (inf where vz is out of reach).
    'down' jumps as ready once up, and 'aloft' as moving out of its landing; the caller adds those delays."""
    vz = np.asarray(vz, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        net = dance.F0_N_KG * (1 - vz / (2 * dance.V0_M_S)) - g
        s = np.where(net > 0, vz * vz / (2 * np.where(net > 0, net, 1.0)), np.inf)
        ok = (vz > 0) & (s <= body.push_off_m + 1e-9)
        c = envelope.dip_s(1.0, g)
        q = {'down': 'ready', 'aloft': 'moving'}.get(posture, posture)
        if q == 'moving':
            out = np.full(vz.shape, envelope.STRIDE_S)
        elif q in ('ready', 'split'):
            dip = np.maximum(0.0, s - body.crouch_push_m)
            out = c * np.sqrt(dip) + (body.crouch_push_m + dip + s) / vz
        else:
            lift = body.push_off_m - body.dip_m
            dip = np.maximum(MIN_DIP_M, s - lift)
            out = c * np.sqrt(dip) + (dip + lift + s) / vz
        return np.where(ok, out, np.inf)


def push_ratio(body, posture, ux, uy, grip, facing=0.0):
    """The horizontal speed a take-off's push can add per unit of vz along (ux, uy): min(mu, lead0/h) from a floor
    posture (envelope.jumps' lateral hop), mu from a gait or a landing."""
    if posture in ('moving', 'aloft'):
        return np.full(np.shape(ux), float(grip))
    q = 'ready' if posture == 'down' else posture
    return np.minimum(grip, lead0(body, q, ux, uy, facing) / body.com_m(q))


def _ready_s(body, st, g):
    """When the jump's preparation can start: after getting up, after the landing, or at the next stance."""
    if st.posture == 'down':
        return get_up_s(body, g)
    if st.posture == 'aloft':
        return st.landing_s(body, g)
    return st.next_stance_s if st.posture == 'moving' else 0.0


def _jump_unit(st, dip, sideways):
    q = st.posture
    if q in ('moving', 'aloft'):
        return 'rebound_jump'
    if q == 'upright':
        return 'standing_jump'
    return 'side_hop' if sideways and dip <= 1e-9 else 'crouch_jump'


def jumps(body, posture, g, vz=None, next_stance_s=0.0):
    """A jump from a posture (or a State) at take-off speed vz (default the full take-off): its dip, its delay to
    take-off from the moment the player acts, its airtime, its apex and the hands' highest reach there."""
    st = posture if isinstance(posture, State) else State(posture=posture, next_stance_s=next_stance_s)
    vz = body.take_off_m_s(g) if vz is None else float(vz)
    q = {'down': 'ready', 'aloft': 'moving'}.get(st.posture, st.posture)
    timing = body.jump_timing(q, vz, g)
    if timing is None:
        return None
    dip, dip_s, push_s = timing
    lift = _ready_s(body, st, g) + dip_s + push_s
    z_to = take_off_height(body)
    return dict(unit=_jump_unit(st, dip, False), take_off_m_s=vz, dip_m=dip, dip_s=dip_s, push_s=push_s,
                take_off_s=lift, airtime_s=float(airtime(body, vz, g)), apex_s=lift + vz / g,
                rise_m=vz * vz / (2 * g),
                hands_top_m=z_to + vz * vz / (2 * g) + body.reach_up_m - body.com_height_m)


# --- the evaluator ---------------------------------------------------------------------------------------------------

def _key(state):
    """The state the floor profiles depend on (its place does not matter), as a hashable cache key."""
    return State(p=(0.0, 0.0), v=(float(state.v[0]), float(state.v[1])), posture=state.posture,
                 facing=float(state.facing), z=None if state.z is None else float(state.z), vz=float(state.vz),
                 next_stance_s=float(state.next_stance_s))


@lru_cache(maxsize=16)
def _table(body, key, grip, g):
    th = np.arange(N_DIR) * (2 * math.pi / N_DIR)
    t = np.arange(int(round(SPAN_S / DT)) + 1) * DT
    X, V = _profiles(body, key, np.cos(th), np.sin(th), t, grip, g)
    return np.ascontiguousarray(X).ravel(), np.ascontiguousarray(V).ravel(), t.size


def _lookup(tab, ang, t):
    """The table's x and x' at headings ang and times t (broadcast): linear in time and in direction; beyond the span
    the last speed holds."""
    X, V, nt = tab
    fa = np.mod(ang, 2 * math.pi) * (N_DIR / (2 * math.pi))
    ia = np.floor(fa)
    wa = fa - ia
    ia = ia.astype(np.intp) % N_DIR
    ib = (ia + 1) % N_DIR
    ft = np.maximum(t, 0.0) / DT
    it = np.minimum(np.floor(ft), nt - 2)
    wt = ft - it
    w = np.minimum(wt, 1.0)
    it = it.astype(np.intp)
    i0, i1 = ia * nt + it, ib * nt + it
    v0 = V[i0] + (V[i0 + 1] - V[i0]) * w
    v1 = V[i1] + (V[i1 + 1] - V[i1]) * w
    x0 = X[i0] + (X[i0 + 1] - X[i0]) * w + (wt - w) * DT * V[i0 + 1]
    x1 = X[i1] + (X[i1 + 1] - X[i1]) * w + (wt - w) * DT * V[i1 + 1]
    return x0 + (x1 - x0) * wa, v0 + (v1 - v0) * wa


def _evaluate(body, state, grip, g, t, q, details=False):
    """Margins and mode codes for the hands at points q (n, 3) at times t (n,); with details, the best jump's terms."""
    t = np.asarray(t, float).ravel()
    q = np.asarray(q, float).reshape(-1, 3)
    tab = _table(body, _key(state), float(grip), float(g))
    dx = q[:, 0] - (state.p[0] + state.v[0] * t)
    dy = q[:, 1] - (state.p[1] + state.v[1] * t)
    qz = q[:, 2]
    d = np.hypot(dx, dy)
    ang = np.where(d > 1e-9, np.arctan2(dy, dx), state.facing)
    R, top = body.reach_out_m, body.reach_up_m
    a_low, a_up = body.com_height_m - body.reach_low_m, body.reach_up_m - body.com_height_m
    x, _ = _lookup(tab, ang, t)
    p = state.posture
    if p == 'aloft':
        T = state.landing_s(body, g)
        zc = state.z + state.vz * t - 0.5 * g * t * t
        air = np.minimum(R - d, np.minimum(qz - (zc - a_low), zc + a_up - qz))
        lo = np.maximum(body.reach_low_m - sink_m(t - T, g), 0.0)
        land = np.minimum(x + R - d, np.minimum(qz - lo, top - qz))
        m0, c0 = np.where(t <= T, air, land), np.where(t <= T, 2, 3)
    else:
        t0 = get_up_s(body, g) if p == 'down' else 0.0
        low0 = body.hand_band(body.com_m('ready' if p == 'down' else p))[0]
        lo = np.maximum(low0 - sink_m(t - t0, g), 0.0)
        hi = np.full(t.shape, top)
        if p == 'down':
            up = t >= t0
            lo, hi = np.where(up, lo, 0.0), np.where(up, hi, body.reach_low_m)
        m0, c0 = np.minimum(x + R - d, np.minimum(qz - lo, hi - qz)), np.zeros(t.shape, int)
    m1, info = _jumps(body, state, grip, g, t, ang, d, qz, tab, details)
    jump = (m0 < 0) & (m1 > m0)
    return np.where(jump, m1, m0), np.where(jump, 1, c0), info


def _jumps(body, state, grip, g, t, ang, d, qz, tab, details):
    """The best jump for each point: over the take-off speed and the time in the air, the take-off at t - tau after
    floor travel to t_f = t - tau - prep."""
    n = t.size
    full = body.take_off_m_s(g)
    z_to = take_off_height(body)
    a_low, a_up = body.com_height_m - body.reach_low_m, body.reach_up_m - body.com_height_m
    ready = _ready_s(body, state, g)
    vz = np.broadcast_to(np.array(VZ_SHARES) * full, (n, len(VZ_SHARES)))
    apex_vz = np.sqrt(2 * g * np.maximum(qz - a_up - z_to, 0.0))
    extra = [apex_vz, 1.15 * apex_vz]
    if state.posture in ('ready', 'split', 'down'):
        extra.append(np.full(n, body.take_off_m_s(g, body.crouch_push_m)))
    vz = np.clip(np.concatenate([vz, np.stack(extra, 1)], 1), 0.02 * full, full)       # (n, J)
    prep = prep_s(body, state.posture, vz, g)
    tau_max = np.minimum(airtime(body, vz, g), t[:, None] - ready - prep)              # (n, J)
    centre = qz[:, None] - 0.5 * (a_up - a_low) - z_to
    sq = np.sqrt(np.maximum(vz * vz - 2 * g * centre, 0.0))
    special = np.stack([(vz - sq) / g, (vz + sq) / g, vz / g], -1)
    tau = np.concatenate([tau_max[..., None] * np.linspace(0.0, 1.0, N_TAU), special], -1)
    tau = np.clip(tau, 0.0, np.maximum(tau_max, 0.0)[..., None])                      # (n, J, K)
    tf = t[:, None, None] - tau - prep[..., None]
    ux, uy = np.cos(ang), np.sin(ang)
    k = push_ratio(body, state.posture, ux, uy, grip, state.facing)
    xf, vf = _lookup(tab, ang[:, None, None], tf)
    D = xf + vf * (prep[..., None] + tau) + (k[:, None] * vz)[..., None] * tau
    zc = z_to + vz[..., None] * tau - 0.5 * g * tau * tau
    m = np.minimum(D + body.reach_out_m - d[:, None, None],
                   np.minimum(qz[:, None, None] - zc + a_low, zc + a_up - qz[:, None, None]))
    m = np.where((tau_max >= 0)[..., None] & np.isfinite(prep)[..., None], m, -np.inf).reshape(n, -1)
    i = m.argmax(1)
    best = m[np.arange(n), i]
    info = None
    if details:
        j, kk = np.unravel_index(i, tau.shape[1:])
        r = np.arange(n)
        info = dict(vz=vz[r, j], tau=tau[r, j, kk], prep=prep[r, j], tf=tf[r, j, kk], k=k)
    return best, info


# --- the answers -----------------------------------------------------------------------------------------------------

def _travel(body, st, ang, t, grip, g):
    """The vocabulary units for the floor travel up to t along heading ang."""
    ux, uy = math.cos(ang), math.sin(ang)
    p = st.posture
    x = float(_profiles(body, st, ux, uy, [t], grip, g)[0][0, 0])
    if p == 'aloft':
        L = landing(body, st, grip, g)
        if t <= L['landing_s']:
            return ()
        return ('braking_landing',) + (('bound',) if t > L['landing_s'] + L['stance_s'] + 1e-9 else ())
    if p == 'moving':
        return ('bound',) if x > 1e-3 else ()
    t0 = get_up_s(body, g) if p == 'down' else 0.0
    q = 'ready' if p == 'down' else p
    pre = ('fall_and_rise',) if p == 'down' and t >= t0 else ()
    if x <= 1e-3:
        return pre
    c, s = math.cos(st.facing), math.sin(st.facing)
    stride = 'shuffle' if abs(-ux * s + uy * c) > abs(ux * c + uy * s) else 'bound'
    lean = float(_profiles(body, st, ux, uy, [t], grip, g, steps=False)[0][0, 0])
    if x > lean + 1e-9:
        first, start, L, h = 'drop_step', t0 + STEP_S, body.split_reach_m, body.ready_com_m
    else:
        first, start, L, h = 'lean_start', t0, float(lead0(body, q, ux, uy, st.facing)), body.com_m(q)
    end = start + math.sqrt(h / g) * math.acosh(max(grip * h / L, 1.0))
    return pre + (first,) + ((stride,) if t > end + 1e-9 else ())


def reach_at(body, state, point, t, grip, g):
    """Whether the hands can be at point (x, y, z) at exactly time t, on the floor in the standing band (or a scoop
    below it) or by a jump timed to it, and how: a Reach with the mode, the units, the time committed to an arc and,
    for a jump, the take-off."""
    t = float(t)
    m, code, info = _evaluate(body, state, grip, g, [t], [point], details=True)
    m, mode = float(m[0]), MODES[int(code[0])]
    px, py = state.p[0] + state.v[0] * t, state.p[1] + state.v[1] * t
    dx, dy = point[0] - px, point[1] - py
    ang = math.atan2(dy, dx) if math.hypot(dx, dy) > 1e-9 else state.facing
    T = state.landing_s(body, g)
    if mode == 'aloft':
        return Reach(ok=m >= 0, seconds=t, margin_m=m, mode=mode, moves=(), committed_s=t, take_off_s=0.0)
    if mode in ('ground', 'landing'):
        moves = _travel(body, state, ang, t, grip, g)
        q = 'ready' if state.posture == 'down' else ('moving' if state.posture == 'aloft' else state.posture)
        if point[2] < body.hand_band(body.com_m(q))[0]:
            moves += ('scoop',)
        return Reach(ok=m >= 0, seconds=t, margin_m=m, mode=mode, moves=moves,
                     committed_s=T if mode == 'landing' else 0.0, take_off_s=0.0 if mode == 'landing' else None)
    tau, tf, vz = float(info['tau'][0]), float(info['tf'][0]), float(info['vz'][0])
    ready = _ready_s(body, state, g)
    moves = _travel(body, state, ang, tf, grip, g) if tf > ready + 1e-6 or state.posture == 'aloft' else ()
    if state.posture == 'down' and 'fall_and_rise' not in moves:
        moves = ('fall_and_rise',) + moves
    q = {'down': 'ready', 'aloft': 'moving'}.get(state.posture, state.posture)
    timing = body.jump_timing(q, vz, g)
    c, s = math.cos(state.facing), math.sin(state.facing)
    ux, uy = math.cos(ang), math.sin(ang)
    sideways = abs(-ux * s + uy * c) > abs(ux * c + uy * s) and tf <= ready + 1e-6
    moves += (_jump_unit(state, timing[0] if timing else 0.0, sideways),)
    return Reach(ok=m >= 0, seconds=t, margin_m=m, mode='jump', moves=moves,
                 committed_s=tau + (T if state.posture == 'aloft' else 0.0), take_off_s=t - tau)


def reach_time(body, state, point, grip, g, t_max=10.0, dt=0.01):
    """The earliest time the hands can be at a fixed point: a scan on dt, then bisection to a millisecond's
    fraction. Unreachable within t_max: ok False, seconds inf, and the closest the scan came as the margin."""
    ts = np.arange(int(math.floor(t_max / dt + 1e-9)) + 1) * dt
    q = np.broadcast_to(np.asarray(point, float), (ts.size, 3))
    m, code, _ = _evaluate(body, state, grip, g, ts, q)
    ok = m >= 0
    if not ok.any():
        i = int(np.argmax(m))
        return Reach(ok=False, seconds=math.inf, margin_m=float(m[i]), mode=MODES[int(code[i])])
    i = int(np.argmax(ok))
    lo, hi = (float(ts[i - 1]), float(ts[i])) if i > 0 else (0.0, 0.0)
    while hi - lo > 2e-4:
        mid = 0.5 * (lo + hi)
        if _evaluate(body, state, grip, g, [mid], [point])[0][0] >= 0:
            hi = mid
        else:
            lo = mid
    return reach_at(body, state, point, hi, grip, g)


def reach_along(body, state, times, points, grip, g, chunk=500):
    """reach_at over many (time, point) pairs at once, as arrays: ok, margin_m and the mode code (an index into MODES:
    the floor where it reaches, else the best of the rest). 10,000 points take a fraction of a second."""
    points = np.asarray(points, float).reshape(-1, 3)
    times = np.broadcast_to(np.asarray(times, float).ravel(), (points.shape[0],))
    margin, code = np.empty(points.shape[0]), np.empty(points.shape[0], np.int8)
    for s in range(0, points.shape[0], chunk):
        m, c, _ = _evaluate(body, state, grip, g, times[s:s + chunk], points[s:s + chunk])
        margin[s:s + chunk], code[s:s + chunk] = m, c
    return margin >= 0, margin, code
