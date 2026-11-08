"""A chase: a pursuer against an evader on a floor, on the Moon beside Earth, and what decides it.

    from humanities.sport import contest_chase
    contest_chase.screen()      # the contest's results; contests.py writes them into results/contests.json

Two players start still, GAP_M apart along x, centred on an open floor or in a square court of 10, 20 or 40 m, and
both start on the same cue. Each is a centre of mass under the reach engine's floor law for a player in a gait
(reach.py, 'moving'): it pushes in any direction at up to min(mu g, the legs' sprint profile at its speed along the
push), and no push takes it past V_MAX_M_S; a wall stops a centre of mass WALL_M short of it and takes away the
velocity into it. The pursuer tags the evader when its hands reach the evader's centre of mass: within reach_out on
the floor plan and inside the pursuer's standing band in height (state.Body.hand_band). Time steps of DT.

Each player sees the other as it was REACTION_S ago, carried forward at the velocity it had then; an evader in the air
is therefore read exactly until it lands, which is how the pursuer sees each take-off and reacts REACTION_S later.
The strategies, kept simple:
- the 'extrapolating' pursuer heads for the evader's position extrapolated by the time to go, the earliest time its
  reach covers that place, pushing toward the place from its own drift point (the double integrator's heading, as in
  reach.py); the 'anticipating' pursuer does the same on the floor, but whenever it sees the evader in the air it sends
  itself to the evader's landing point at the landing time, so as to cover a cut either way;
- the evader runs ('bound') to one of two plans. 'evade': it pushes along whichever of HEADINGS directions keeps the
  most distance beyond the pursuer's reach a horizon ahead (the court's crossing time at its push, 0.8-3 s), with
  places beyond the walls clipped and penalised. 'circuit' (in a court): it runs a loop CIRCUIT_INSET_M inside the
  walls' reach at the speed its turning allows there, v^2 / R = CIRCUIT_SHARE of its push, steering toward the loop's
  tangent; after a tactic it carries on round the loop in whichever sense it is then turning. On an open floor both
  plans flee straight away;
- its tactics, used whenever the pursuer comes within TRIGGER_M (centre to centre), at most once per REFRACT_S of
  running: a cut ('bounce_cut' through CUT_DEG at its speed, with the bounce envelope.cut_bounce needs on the way in
  and out: the bounce in is the tell, the bounce out commits); a hop ('rebound_jump' after one stride's stance, rising
  HOP_RISE_M or as high as the take-off allows, its push adding mu vz along the evader's best heading, then the
  landing's 'braking_landing' stance); a feint (the cut's bounce in, then no cut).

What the contest reports for each court, world and grip: the time to the tag (None if none within T_MAX_S); for each
tactic against the plain run, the change in that time, and the evader's time in the air; how often a wall stopped
either player; the closest approach on an open floor, where equal players should never close the gap; and the bodies'
match-ups, the larger body chasing the dance body and fleeing it.

Assumptions, beside the shared ones in state.py: both players start already in a gait (posture 'moving', next stance
now), so the lean and the drop step are left out of the start; a run's strides are averaged into the floor law, and
only a cut's or a hop's bounces are modelled as flights; a cut's stance turns the velocity at a steady rate; the
pursuer stays on the floor and its hands reach only its standing band; a player plans its own reach as the grip's
parabola capped at V_MAX_M_S (the muscles' fall-off is left to the motion itself); the court's walls are hard.
"""
from __future__ import annotations

import math

from humanities.sport import envelope as env
from humanities.sport import reach, state

COURTS_M = (None, 10.0, 20.0, 40.0)       # None: an open floor
GAP_M = 4.0
DT = 0.01
T_MAX_S = 60.0
WALL_M = 0.25
TRIGGER_M = 2.0
REFRACT_S = 0.5
CUT_DEG = 60.0
CUT_MIN_M_S = 1.0                         # slower than this, the evader does not cut
HOP_RISE_M = 1.0
CIRCUIT_INSET_M = 1.0                     # the loop's radius is the walls' limit for a centre of mass less this
CIRCUIT_SHARE = 0.8                       # the share of the push the loop's turn takes at its steady speed
CIRCUIT_GAIN = 1.0                        # 1/s: how hard the loop's steering pulls back onto its radius
CIRCUIT_TAU = 0.2                         # s: the time the steering takes to close a velocity error
HEADINGS = 16
CLIP_WEIGHT = 2.0                         # the evader's penalty per metre a planned place lies beyond the walls
T_GO = tuple(0.05 * i for i in range(1, 61)) + tuple(3.0 + 0.25 * i for i in range(1, 13))
TACTICS = ('run', 'cut', 'hop', 'feint')
MATCHES = (                               # reader, plan, tactic, pursuer, evader
    ('extrapolating', 'evade', 'run', 'dancer', 'dancer'),
    ('extrapolating', 'circuit', 'run', 'dancer', 'dancer'),
    ('extrapolating', 'circuit', 'cut', 'dancer', 'dancer'),
    ('extrapolating', 'circuit', 'hop', 'dancer', 'dancer'),
    ('extrapolating', 'circuit', 'feint', 'dancer', 'dancer'),
    ('anticipating', 'circuit', 'cut', 'dancer', 'dancer'),
    ('anticipating', 'circuit', 'feint', 'dancer', 'dancer'),
    ('extrapolating', 'circuit', 'run', 'larger', 'dancer'),
    ('extrapolating', 'circuit', 'run', 'dancer', 'larger'),
    ('extrapolating', 'evade', 'run', 'larger', 'dancer'),
)
WORKERS = 3
_COS = tuple(math.cos(2 * math.pi * k / HEADINGS) for k in range(HEADINGS))
_SIN = tuple(math.sin(2 * math.pi * k / HEADINGS) for k in range(HEADINGS))


def planned_reach_m(T, a):
    """How far a player plans to get in T from its drift point: the grip's parabola, then V_MAX_M_S."""
    tc = state.V_MAX_M_S / a
    return 0.5 * a * T * T if T <= tc else 0.5 * a * tc * tc + state.V_MAX_M_S * (T - tc)


def cut_bounce_m_s(speed, angle_deg, grip, g):
    """The vertical speed each way of the bounce a cut through one stance needs (envelope.cut_bounce, unrounded)."""
    dv = 2 * speed * math.sin(math.radians(angle_deg) / 2)
    return max(0.0, (dv / grip - g * env.STRIDE_S) / 2)


def _push(vx, vy, dx, dy, a_max, dt):
    """One step of a push along the unit (dx, dy) at a_max, no push taking the speed past V_MAX_M_S."""
    nvx, nvy = vx + a_max * dx * dt, vy + a_max * dy * dt
    s, lim = math.hypot(nvx, nvy), max(state.V_MAX_M_S, math.hypot(vx, vy))
    return (nvx * lim / s, nvy * lim / s) if s > lim else (nvx, nvy)


def _floor_a(vx, vy, dx, dy, grip, g):
    """The floor law's push for a player in a gait: the grip, or the legs' sprint profile at the speed along it."""
    return min(grip * g, state.muscle_push(vx * dx + vy * dy))


def _inside(x, y, lim):
    """A place clipped into the court (unchanged on an open floor)."""
    return (x, y) if lim is None else (min(max(x, -lim), lim), min(max(y, -lim), lim))


def _wall(x, y, vx, vy, lim):
    hit = False
    if lim is not None:
        if x > lim or x < -lim:
            x, vx, hit = math.copysign(lim, x), (min(vx, 0.0) if x > 0 else max(vx, 0.0)), True
        if y > lim or y < -lim:
            y, vy, hit = math.copysign(lim, y), (min(vy, 0.0) if y > 0 else max(vy, 0.0)), True
    return x, y, vx, vy, hit


def _best_heading(ex, ey, evx, evy, qx, qy, qvx, qvy, horizon, own_m, their_m, lim):
    """The evader's heading: the most distance beyond the pursuer's planned reach at the horizon, places beyond the
    walls clipped and penalised. Returns (k, score)."""
    bx, by = ex + evx * horizon, ey + evy * horizon
    px, py = qx + qvx * horizon, qy + qvy * horizon
    best, kb = -1e9, 0
    for k in range(HEADINGS):
        cx, cy = bx + own_m * _COS[k], by + own_m * _SIN[k]
        clip = 0.0
        if lim is not None:
            lx, ly = min(max(cx, -lim), lim), min(max(cy, -lim), lim)
            clip, cx, cy = math.hypot(cx - lx, cy - ly), lx, ly
        s = math.hypot(cx - px, cy - py) - their_m - CLIP_WEIGHT * clip
        if s > best:
            best, kb = s, k
    return kb, best


def _circuit(ex, ey, evx, evy, radius, speed, sense, a_max):
    """The loop's steering: toward the tangent at `speed`, pulled back onto `radius`; the push's direction and size."""
    r = math.hypot(ex, ey)
    ux, uy = (ex / r, ey / r) if r > 1e-6 else (1.0, 0.0)
    k = CIRCUIT_GAIN * (radius - r)
    dvx, dvy = -uy * sense * speed + k * ux - evx, ux * sense * speed + k * uy - evy
    n = math.hypot(dvx, dvy)
    return (1.0, 0.0, 0.0) if n < 1e-9 else (dvx / n, dvy / n, min(a_max, n / CIRCUIT_TAU))


def contest(g, grip, court=None, tactic='run', reader='extrapolating', pursuer='dancer', evader='dancer', dt=DT,
            t_max=T_MAX_S, gap=GAP_M, trace=False, plan='circuit'):
    """One chase. Returns its summary (and, with trace, the states step by step)."""
    P, E = state.BODIES[pursuer], state.BODIES[evader]
    lim = None if court is None else court / 2 - WALL_M
    a_p, a_e = min(grip * g, state.SPRINT_F0_N_KG), min(grip * g, state.SPRINT_F0_N_KG)
    horizon = 2.0 if court is None else min(3.0, max(0.8, math.sqrt(court / a_e)))
    own_m, their_m = planned_reach_m(horizon, a_e), planned_reach_m(horizon, a_p) + P.reach_out_m
    loop = plan == 'circuit' and lim is not None
    radius = (lim - CIRCUIT_INSET_M) if loop else 0.0
    loop_speed = min(state.V_MAX_M_S, math.sqrt(CIRCUIT_SHARE * a_e * radius)) if loop else 0.0
    sense = 1.0
    reach_p = [planned_reach_m(T, a_p) + P.reach_out_m for T in T_GO]
    n_delay = int(round(state.REACTION_S / dt))
    R = n_delay * dt
    lo, hi = P.reach_low_m, P.reach_up_m
    z_to = reach.take_off_height(E)
    px, py, pvx, pvy = -gap / 2, 0.0, 0.0, 0.0
    ex, ey, evx, evy = gap / 2, 0.0, 0.0, 0.0
    ez, evz = None, 0.0
    hist_p, hist_e = [(px, py, pvx, pvy)], [(ex, ey, evx, evy, None, 0.0)]
    mode, left, seg = 'run', 0.0, []          # the evader's segment, its time left, and the segments queued
    a_seg = (0.0, 0.0)
    next_tactic, used, air_s, walls_p, walls_e = 0.0, 0, 0.0, 0, 0
    min_gap, tag, log = gap, None, []
    steps = int(round(t_max / dt))
    for i in range(1, steps + 1):
        t = i * dt
        # --- what each sees: the other REACTION_S ago, carried forward ------------------------------------------
        ox, oy, ovx, ovy, oz, ovz = hist_e[max(i - 1 - n_delay, 0)]
        qx, qy, qpx, qpy = hist_p[max(i - 1 - n_delay, 0)]
        qx, qy = qx + qpx * R, qy + qpy * R
        # --- the pursuer ----------------------------------------------------------------------------------------
        dx = dy = 0.0
        if reader == 'anticipating' and oz is not None:
            t_land = (ovz + math.sqrt(max(ovz * ovz + 2 * g * (oz - E.com_height_m), 0.0))) / g
            T = t_land - R
            if T > 0:
                tx, ty = _inside(ox + ovx * t_land, oy + ovy * t_land, lim)
                dx, dy = tx - px - pvx * T, ty - py - pvy * T
        if dx == 0.0 and dy == 0.0:
            cx, cy = ox + ovx * R, oy + ovy * R
            for T, rp in zip(T_GO, reach_p):
                tx, ty = _inside(cx + ovx * T, cy + ovy * T, lim)    # the evader cannot be beyond the walls
                dx, dy = tx - px - pvx * T, ty - py - pvy * T
                if dx * dx + dy * dy <= rp * rp:
                    break
        n = math.hypot(dx, dy)
        if n < 1e-9:
            dx, dy, n = ovx - pvx, ovy - pvy, math.hypot(ovx - pvx, ovy - pvy) or 1.0
        dx, dy = dx / n, dy / n
        pa = _floor_a(pvx, pvy, dx, dy, grip, g)
        pvx, pvy = _push(pvx, pvy, dx, dy, pa, dt)
        px, py = px + pvx * dt, py + pvy * dt
        px, py, pvx, pvy, hit = _wall(px, py, pvx, pvy, lim)
        walls_p += hit
        # --- the evader -------------------------------------------------------------------------------------------
        if mode == 'run' and t >= next_tactic and tactic != 'run' and math.hypot(ex - qx, ey - qy) < TRIGGER_M:
            speed = math.hypot(evx, evy)
            if tactic == 'hop':
                seg, used = [('prep', env.STRIDE_S), ('hop', 0.0)], used + 1
            elif speed >= CUT_MIN_M_S:
                th = math.atan2(evy, evx)
                outs = [(math.cos(th + s * math.radians(CUT_DEG)), math.sin(th + s * math.radians(CUT_DEG)))
                        for s in (1, -1)]
                scores = [_best_heading(ex, ey, speed * ux, speed * uy, qx, qy, qpx, qpy, horizon, 0.0, their_m,
                                        lim)[1] for ux, uy in outs]
                ux, uy = outs[0] if scores[0] >= scores[1] else outs[1]
                w = cut_bounce_m_s(speed, CUT_DEG, grip, g)
                turn = ('stance', env.STRIDE_S, (speed * ux - evx) / env.STRIDE_S, (speed * uy - evy) / env.STRIDE_S)
                if tactic == 'cut':
                    seg = ([('air', w)] if w > 0 else []) + [turn] + ([('air', w)] if w > 0 else [])
                elif w > 0:
                    seg = [('air', w)]
                used += bool(seg)
            if seg:
                mode = None
        if mode is None:                          # start the next queued segment
            if not seg:
                mode, next_tactic = 'run', t + REFRACT_S
                turning = ex * evy - ey * evx
                sense = math.copysign(1.0, turning) if abs(turning) > 1e-6 else sense
            else:
                kind, *rest = seg.pop(0)
                if kind == 'air':
                    mode, ez, evz = 'air', E.com_height_m, rest[0]
                elif kind == 'stance':
                    mode, left, a_seg = 'stance', rest[0], (rest[1], rest[2])
                elif kind == 'prep':
                    mode, left = 'prep', rest[0]
                else:                             # the hop leaves the floor along the best heading
                    k, _ = _best_heading(ex, ey, evx, evy, qx, qy, qpx, qpy, horizon, own_m, their_m, lim)
                    vz = min(math.sqrt(2 * g * HOP_RISE_M), E.take_off_m_s(g))
                    evx, evy = evx + grip * vz * _COS[k], evy + grip * vz * _SIN[k]
                    mode, ez, evz = 'air', z_to, vz
                    seg = [('land', 0.0)] + seg
        if mode == 'air':
            evz_new = evz - g * dt
            ez, evz = ez + 0.5 * (evz + evz_new) * dt, evz_new
            ex, ey = ex + evx * dt, ey + evy * dt
            air_s += dt
            if ez <= E.com_height_m:
                v_down, ez, evz, mode = -evz, None, 0.0, None
                if seg and seg[0][0] == 'land':
                    seg.pop(0)
                    mode, left = 'land', 2 * env.STROKE_M / max(v_down, 0.05)
                    a_land = grip * (max(v_down, 0.05) / left + g)
        elif mode == 'stance':
            evx, evy = evx + a_seg[0] * dt, evy + a_seg[1] * dt
            ex, ey = ex + evx * dt, ey + evy * dt
        else:                                     # running, preparing a hop, or in a landing's stance
            if loop:
                ux, uy, _ = _circuit(ex, ey, evx, evy, radius, loop_speed, sense, 1.0)
                ea = a_land if mode == 'land' else _floor_a(evx, evy, ux, uy, grip, g)
                ux, uy, ea = _circuit(ex, ey, evx, evy, radius, loop_speed, sense, ea)
            else:
                k, _ = _best_heading(ex, ey, evx, evy, qx, qy, qpx, qpy, horizon, own_m, their_m, lim)
                ux, uy = _COS[k], _SIN[k]
                ea = a_land if mode == 'land' else _floor_a(evx, evy, ux, uy, grip, g)
            evx, evy = _push(evx, evy, ux, uy, ea, dt)
            ex, ey = ex + evx * dt, ey + evy * dt
        if mode in ('stance', 'prep', 'land'):
            left -= dt
            if left <= 1e-9:
                mode = None if seg else 'run'
                if mode == 'run':
                    next_tactic = t + REFRACT_S
                    turning = ex * evy - ey * evx
                    sense = math.copysign(1.0, turning) if abs(turning) > 1e-6 else sense
        ex, ey, evx, evy, hit = _wall(ex, ey, evx, evy, lim)
        walls_e += hit
        hist_p.append((px, py, pvx, pvy))
        hist_e.append((ex, ey, evx, evy, ez, evz))
        d = math.hypot(ex - px, ey - py)
        min_gap = min(min_gap, d)
        if trace:
            log.append((t, px, py, pvx, pvy, pa, ex, ey, evx, evy, ez, mode))
        z = E.com_height_m if ez is None else ez
        if d <= P.reach_out_m and lo <= z <= hi:
            tag = t
            break
    out = dict(tag_s=None if tag is None else round(tag, 2), tactics_used=used, evader_air_s=round(air_s, 2),
               loop=dict(radius_m=round(radius, 2), speed_m_s=round(loop_speed, 2)) if loop else None,
               wall_stops=dict(pursuer=walls_p, evader=walls_e), closest_m=round(min_gap, 2))
    if trace:
        out['trace'] = log
    return out


def _case(args):
    world, g, grip, court, reader, plan, tactic, pursuer, evader = args
    out = contest(g, grip, court, tactic, reader, pursuer, evader, plan=plan)
    return dict(world=world, grip=grip, court_m=court, reader=reader, plan=plan, tactic=tactic, pursuer=pursuer,
                evader=evader, **out)


def _effects(cases):
    """Against the plain run in the same court, world and grip: how each tactic and reader moves the tag."""
    base = {(c['world'], c['grip'], c['court_m'], c['plan']): c['tag_s'] for c in cases
            if c['tactic'] == 'run' and c['pursuer'] == c['evader'] == 'dancer'}
    out = []
    for c in cases:
        if c['tactic'] == 'run' and c['pursuer'] == c['evader'] == 'dancer':
            continue
        b = base[(c['world'], c['grip'], c['court_m'], c['plan'])]
        if b is None or c['tag_s'] is None:
            change = None if b is None and c['tag_s'] is None else ('escapes' if c['tag_s'] is None else 'tagged')
        else:
            change = round(c['tag_s'] - b, 2)
        out.append(dict(world=c['world'], grip=c['grip'], court_m=c['court_m'], reader=c['reader'],
                        plan=c['plan'], tactic=c['tactic'], pursuer=c['pursuer'], evader=c['evader'], run_tag_s=b,
                        tag_s=c['tag_s'], change_s=change))
    return out


def screen(worlds=None, grips=None, courts=COURTS_M, workers=WORKERS):
    worlds = worlds or list(state.WORLDS)
    grips = grips or state.GRIPS
    jobs = [(w, g, grip, court, *m) for w, g in worlds for grip in grips for court in courts for m in MATCHES]
    if workers > 1:
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context('fork')) as pool:
            cases = list(pool.map(_case, jobs))
    else:
        cases = [_case(j) for j in jobs]
    return dict(gap_m=GAP_M, reaction_s=state.REACTION_S, t_max_s=T_MAX_S, dt_s=DT, trigger_m=TRIGGER_M,
                cut_deg=CUT_DEG, hop_rise_m=HOP_RISE_M, units=dict(run='bound', cut='bounce_cut',
                                                                   hop=['rebound_jump', 'braking_landing'],
                                                                   feint='bounce_cut'),
                cases=cases, effects=_effects(cases))
