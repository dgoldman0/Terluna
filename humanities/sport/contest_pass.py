"""A pass against an interceptor: whether a ball gets from one player to another past a third, on the Moon beside
Earth, and how that contest is decided.

    from humanities.sport import contest_pass
    contest_pass.screen()       # the contest's results; contests.py writes them into results/contests.json

Three players start still in the ready crouch on an open floor: a passer with the ball at the origin, a receiver D
metres along +x facing the passer, and an interceptor at one place of a grid beside the pass line, facing back along
the line toward the passer's end.
The passer throws from a stance (a chest pass with the feet under the body, or overhand from a long split stance) at
any speed the stance holds (windows.release_limit) and any elevation, aimed at the receiver or into the space a few
degrees to either side, and the ball flies in the design air (windows.fly). Receiver and interceptor both start a
reaction time after the release (state.REACTION_S). A pass is caught at the first point of its path the receiver's
hands can reach in time (reach.reach_along) with a ball the receiver's chest can hold, and it is safe when the
interceptor's hands can reach no point of its path before that catch.

What the contest reports for each world, grip, distance and stance:
- the shadow: the share of the grid from which the interceptor makes every pass unsafe, which is the floor one player
  controls around a pass line, and its area;
- outside the shadow, the best safe pass's margin (how long after the catch the interceptor could first have touched
  the ball) and whether it goes flat or as a lob;
- how the interceptor has to decide: the shadow when it reacts to the release, and when it guesses the pass and starts
  moving HEAD_STARTS_S before the release; the head start a place needs is how early the interceptor must commit to
  control it, and floor controlled only by guessing is where a fake pays;
- a fake's worth: the shadow when a fake has drawn the interceptor into a jump that left the floor BAIT_S before the
  real release, set against the floor the interceptor controlled by guessing.

Assumptions, beside the shared ones in state.py: the players are the dance body unless named (the larger body
intercepts in one variant); a pass ends where the ball reaches the floor; a touch is an interception; the receiver
gives with the ball over 0.30 s, which holds twice what the 0.15 s give holds.
"""
from __future__ import annotations

import math

import numpy as np

from humanities.sport import envelope as env
from humanities.sport import reach, state, windows

BALL = env.THROW_BALL
DISTANCES_M = (6.0, 10.0, 15.0)
STANCES = {'chest pass, feet under': ('feet_under', env.AIR_LEVERS_M['chest pass, near the centre of mass']),
           'overhand, split stance': ('split', None)}
ELEVATIONS_DEG = tuple(range(-5, 66, 10))
AZIMUTHS_DEG = (0.0, 8.0, -8.0)
SPEEDS = 4                                # release speeds from MIN_SPEED to the stance's limit
MIN_SPEED_M_S = 3.0
GRID_DX_M = 1.5                           # the interceptor's grid: along the line, and out from it on one side (the
GRID_DY_M = 0.5                           # passes aim to both sides, so the shadow is twice the one side's)
GRID_Y_M = tuple(0.25 + GRID_DY_M * i for i in range(10))
SAMPLE_S = 0.03
HEAD_STARTS_S = (0.3, 1.0)                # an interceptor guessing the pass starts this long before the release
WORKERS = 3                               # processes for the sweep (each well under 1 GB)
BAIT_S = 0.20                             # a fake draws a jump that leaves the floor this long before the release
BAIT_RISE_M = 0.7                         # or the crouch's best, where that is lower (Earth: 0.22 m)
GIVE_S = 0.30


# --- adapters to the libraries (the only calls this module makes on them) --------------------------------------------

def _fly(v0, p0, g, rho):
    """Flights of BALL from p0 for each initial velocity, sampled every SAMPLE_S: times (n,), points (m, n, 3) with NaN
    after the floor."""
    out = windows.fly(BALL, p0, v0, g, rho, dt=SAMPLE_S, t_max=12.0)
    t, xyz = np.asarray(out['t']), np.asarray(out['p'])
    step = max(1, int(round(SAMPLE_S / (t[1] - t[0]))))
    return t[::step], xyz[:, ::step]


def _along(body, st, times, points, grip, g):
    """reach.reach_along: whether the hands can be at each point at its time, and the mode used."""
    out = reach.reach_along(body, st, np.asarray(times, float), np.asarray(points, float), grip, g)
    if isinstance(out, dict):
        return np.asarray(out['ok'], bool), np.asarray(out.get('mode'))
    ok, _, mode = out
    return np.asarray(ok, bool), np.asarray(mode)


# --- one contest -------------------------------------------------------------------------------------------------------

def _passes(passer, stance, lever, g, rho, grip, distance):
    """The passer's fan: velocities and release point, within the stance's hold."""
    kind = stance
    limit = windows.release_limit(passer, kind, BALL, g, grip, lever=lever if lever else passer.release_above_com_m)
    limit = float(limit['m_s'] if isinstance(limit, dict) else limit)
    posture = 'split' if kind == 'split' else 'ready'
    h0 = passer.com_m(posture) + (lever if lever else passer.release_above_com_m)
    speeds = np.linspace(MIN_SPEED_M_S, max(limit, MIN_SPEED_M_S), SPEEDS)
    v0 = []
    for s in speeds:
        for el in ELEVATIONS_DEG:
            for az in AZIMUTHS_DEG:
                e, a = math.radians(el), math.radians(az)
                v0.append((s * math.cos(e) * math.cos(a), s * math.cos(e) * math.sin(a), s * math.sin(e)))
    meta = np.array([(s, el, az) for s in speeds for el in ELEVATIONS_DEG for az in AZIMUTHS_DEG])
    return np.array(v0), np.array([0.0, 0.0, h0]), limit, meta


def _first(ok, valid):
    """Index of the first True along each row (n if none)."""
    hit = ok & valid
    n = hit.shape[1]
    return np.where(hit.any(axis=1), hit.argmax(axis=1), n)


def contest(g, rho, grip, distance, stance, passer=None, receiver=None, interceptor=None, grid=None):
    passer = passer or state.BODIES['dancer']
    receiver = receiver or state.BODIES['dancer']
    interceptor = interceptor or state.BODIES['dancer']
    kind, lever = STANCES[stance]
    v0, p0, limit, meta = _passes(passer, kind, lever, g, rho, grip, distance)
    t, path = _fly(v0, p0, g, rho)
    m, n = path.shape[:2]
    valid = np.isfinite(path[..., 2])
    speed_at = np.zeros((m, n))
    speed_at[:, 1:] = np.linalg.norm(np.diff(np.nan_to_num(path), axis=1), axis=2) / (t[1] - t[0])
    hold = windows.catch_limit(receiver, 'feet_under', BALL, g, grip, give_s=GIVE_S,
                               lever=env.AIR_LEVERS_M['chest pass, near the centre of mass'])
    hold = float(hold['m_s'] if isinstance(hold, dict) else hold)

    # the receiver, reacting at the release
    rstate = state.State(p=(distance, 0.0), posture='ready', facing=math.pi)
    tr = np.repeat((t - state.REACTION_S)[None, :], m, axis=0)
    ok_r, _ = _along(receiver, rstate, tr.ravel(), np.nan_to_num(path).reshape(-1, 3), grip, g)
    ok_r = ok_r.reshape(m, n) & (tr >= 0) & (speed_at <= hold)
    k_r = _first(ok_r, valid)
    caught = k_r < n
    t_catch = np.where(caught, t[np.minimum(k_r, n - 1)], np.inf)

    grid = np.array(grid if grid is not None else [(x, y) for x in np.arange(-1.0, distance + 2.0 + 1e-9, GRID_DX_M)
                                                   for y in GRID_Y_M])
    keep = np.flatnonzero(caught)
    if keep.size == 0:
        return dict(limit_m_s=round(limit, 2), passes=int(m), caught_by_receiver=0, cells=[],
                    summary=_summary([], distance))
    # only caught passes matter, and only their paths up to the latest catch
    n_cut = int(min(n, k_r[keep].max() + 1))
    path, valid, t_catch, meta = path[keep, :n_cut], valid[keep, :n_cut], t_catch[keep], meta[keep]
    t, (m, n) = t[:n_cut], (keep.size, n_cut)
    # every interceptor faces back along the line, toward the passer's end, so one reach table serves the grid
    rel = np.nan_to_num(path)[None] - np.concatenate([grid, np.zeros((len(grid), 1))], axis=1)[:, None, None, :]
    istate = state.State(p=(0.0, 0.0), posture='ready', facing=math.pi)
    vz0 = min(math.sqrt(2 * g * BAIT_RISE_M), interceptor.take_off_m_s(g, interceptor.crouch_push_m))
    z_to = interceptor.com_height_m + interceptor.push_off_m - interceptor.dip_m
    aloft = state.State(p=(0.0, 0.0), posture='aloft', facing=math.pi,
                        z=z_to + vz0 * BAIT_S - 0.5 * g * BAIT_S ** 2, vz=vz0 - g * BAIT_S)
    first = {}
    leads = {'react': (istate, -state.REACTION_S), **{f'head_start_{h:g}_s': (istate, h) for h in HEAD_STARTS_S},
             'baited': (aloft, 0.0)}
    for name, (st, lead) in leads.items():
        ti = np.broadcast_to(t + lead, rel.shape[:3])
        ok, _ = _along(interceptor, st, ti.ravel(), rel.reshape(-1, 3), grip, g)
        ok = ok.reshape(rel.shape[:3]) & (ti >= 0) & valid[None]
        k = np.where(ok.any(axis=2), ok.argmax(axis=2), n)
        first[name] = np.where(k < n, t[np.minimum(k, n - 1)], np.inf)        # (cells, passes)
    cells = []
    for c, q in enumerate(grid):
        row = dict(x=float(q[0]), y=float(q[1]))
        for name in first:
            row[name] = _cell(first[name][c] - t_catch, t_catch, meta)
        cells.append(row)
    timing = interceptor.jump_timing('ready', vz0, g)
    return dict(limit_m_s=round(limit, 2), passes=int(len(v0)), caught_by_receiver=int(keep.size),
                bait_tell_s=None if timing is None else round(timing[1] + timing[2], 3), cells=cells,
                summary=_summary(cells, distance))


def _cell(lead, t_catch, meta):
    """One interceptor place: lead is, per pass, how long after the receiver's catch the interceptor could first touch
    the ball (inf if never). A pass is safe with at least CONTEST_S to spare and contested within it."""
    safe = lead >= state.CONTEST_S
    contested = (lead > -state.CONTEST_S) & ~safe
    if not safe.any():
        return dict(safe=False, contested_share=round(float(contested.mean()), 3))
    i = int(np.argmax(np.where(safe, np.minimum(lead, 99.0), -np.inf)))
    return dict(safe=True, safe_share=round(float(safe.mean()), 3), contested_share=round(float(contested.mean()), 3),
                never_touched=bool(np.isinf(lead[i])), speed_m_s=round(float(meta[i, 0]), 1),
                elevation_deg=int(meta[i, 1]), azimuth_deg=int(meta[i, 2]), catch_s=round(float(t_catch[i]), 2))


def _summary(cells, distance):
    out = {}
    if not cells:
        return dict(no_pass_caught=True)
    names = [k for k in cells[0] if k not in ('x', 'y')]
    for name in names:
        shadow = [c for c in cells if not c[name]['safe']]
        open_ = [c for c in cells if c[name]['safe']]
        out[name] = dict(shadow_share=round(len(shadow) / len(cells), 3),
                         shadow_m2=round(2 * len(shadow) * GRID_DX_M * GRID_DY_M, 1),
                         median_safe_share=(round(float(np.median([c[name]['safe_share'] for c in open_])), 3)
                                            if open_ else None),
                         best_never_touched_share=(round(sum(c[name]['never_touched'] for c in open_) / len(open_), 2)
                                                   if open_ else None),
                         lob_share=(round(sum(c[name]['elevation_deg'] >= 30 for c in open_) / len(open_), 2)
                                    if open_ else None),
                         median_catch_s=(round(float(np.median([c[name]['catch_s'] for c in open_])), 2)
                                         if open_ else None))
    first = HEAD_STARTS_S[0]
    out['guess_only_share'] = round(out[f'head_start_{first:g}_s']['shadow_share'] - out['react']['shadow_share'], 3)
    out['bait_opens_share'] = round(out['react']['shadow_share'] - out['baited']['shadow_share'], 3)
    # a fake pays against a guessing interceptor: the floor its guess controlled that a drawn jump gives up
    out['fake_undoes_share'] = round(out[f'head_start_{first:g}_s']['shadow_share'] - out['baited']['shadow_share'], 3)
    return out


def _case(args):
    world, g, rho, grip, d, stance, interceptor = args
    c = contest(g, rho, grip, d, stance, interceptor=state.BODIES[interceptor])
    row = dict(world=world, grip=grip, distance_m=d, stance=stance, limit_m_s=c['limit_m_s'], passes=c['passes'],
               caught_by_receiver=c['caught_by_receiver'], bait_tell_s=c.get('bait_tell_s'), summary=c['summary'])
    if interceptor != 'dancer':
        row['interceptor'] = interceptor
    return row


def screen(worlds=None, grips=None, distances=DISTANCES_M, workers=WORKERS):
    worlds = worlds or [(w, g, state.RHO_MOON if w == 'moon' else state.RHO_EARTH) for w, g in state.WORLDS]
    grips = grips or state.GRIPS
    jobs = [(w, g, rho, grip, d, stance, 'dancer') for w, g, rho in worlds for grip in grips for d in distances
            for stance in STANCES]
    # the larger body intercepting, at the central grip, 10 m, the chest pass
    jobs += [(w, g, rho, state.MU, 10.0, 'chest pass, feet under', 'larger') for w, g, rho in worlds]
    if workers > 1:
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context('fork')) as pool:
            cases = list(pool.map(_case, jobs))
    else:
        cases = [_case(j) for j in jobs]
    return dict(ball=BALL, reaction_s=state.REACTION_S, contest_s=state.CONTEST_S, head_starts_s=list(HEAD_STARTS_S),
                bait=dict(rise_m=BAIT_RISE_M, before_s=BAIT_S), receiver_give_s=GIVE_S, cases=cases)
