"""A lobbed ball contested by two players: who gets hands on it first, on the Moon beside Earth, and how that is
decided.

    from humanities.sport import contest_lob
    contest_lob.screen()        # the contest's results; contests.py writes them into results/contests.json

A launcher at the origin lobs a ball from 2 m up to land at L, LOB_M along +x, at one of three elevations, with the
speed that lands it there in the design air (windows.fly). Two players start still in the ready crouch, facing back
toward the launcher (along -x): B at a fixed place (3 m short of L on the line, or 3 m beside L), and A at each place of a grid around L.
Both start a reaction time after the release (state.REACTION_S), and each takes the ball at the first point of its
path the hands can reach in time (reach.reach_along), on the floor or by a jump. The earlier hands win; within
CONTEST_S of each other the ball is contested, and goes to hands and bodies.

What the contest reports for each world, grip, elevation and B's place:
- the shares of A's places where A wins, B wins or the ball is contested;
- how players first touch the ball, on the floor or by a jump;
- the slack: how much later than its reaction a player could still have started and touched the ball as early, which
  is the time it waits under the ball and the time in which players can watch each other and change their plans;
- in a contested ball, who holds position: the player with the more slack arrives first under the ball;
- the crowd: how far from the spot under the earliest touch, in each direction, a player can start and still contest
  the ball, and the area that gives, which is how many players a lob gathers;
- the launcher's control: the share of A's places where the elevation alone changes who gets the ball.

Assumptions, beside the shared ones in state.py: both players are the dance body unless named; the ball is the
handball of the throws; a catch is a touch with the hands, with no check of the hold (lobs arrive slowly).
"""
from __future__ import annotations

import math

import numpy as np

from humanities.sport import envelope as env
from humanities.sport import reach, state, windows

BALL = env.THROW_BALL
LOB_M = 12.0
RELEASE_M = env.RELEASE_HEIGHT_M
ELEVATIONS_DEG = (45.0, 60.0, 70.0)
B_PLACES = {'3 m short, on the line': (-3.0, 0.0), '3 m beside': (0.0, 3.0)}
GRID_M = tuple(float(v) for v in range(-8, 9))
CONTEST_S = state.CONTEST_S
SAMPLE_S = 0.04
DELAYS_S = (0.25, 0.5, 1.0, 2.0, 3.0)
WORKERS = 3
CROWD_DIRECTIONS = (0, 45, 90, 135, 180, 225, 270, 315)   # degrees from +x, the launcher behind at 180


def _fly(v0, g, rho):
    out = windows.fly(BALL, np.array([0.0, 0.0, RELEASE_M]), np.atleast_2d(v0), g, rho, dt=SAMPLE_S, t_max=20.0)
    t, xyz = np.asarray(out['t']), np.asarray(out['p'])
    step = max(1, int(round(SAMPLE_S / (t[1] - t[0]))))
    return t[::step], xyz[:, ::step]


def _along(body, st, times, points, grip, g):
    out = reach.reach_along(body, st, np.asarray(times, float), np.asarray(points, float), grip, g)
    if isinstance(out, dict):
        return np.asarray(out['ok'], bool), np.asarray(out.get('mode'))
    ok, _, mode = out
    return np.asarray(ok, bool), np.asarray(mode)


def lob(g, rho, elevation_deg):
    """The lob at this elevation that lands LOB_M away: its speed, times and path (NaN after the floor)."""
    e = math.radians(elevation_deg)
    speeds = np.linspace(3.0, 40.0, 371)
    v0 = np.stack([speeds * math.cos(e), np.zeros_like(speeds), speeds * math.sin(e)], axis=1)
    t, xyz = _fly(v0, g, rho)
    last = np.where(np.isfinite(xyz[..., 0]), xyz[..., 0], -np.inf).max(axis=1)
    i = int(np.argmin(np.abs(last - LOB_M)))
    lo, hi = speeds[max(i - 1, 0)], speeds[min(i + 1, len(speeds) - 1)]
    fine = np.linspace(lo, hi, 101)
    v0 = np.stack([fine * math.cos(e), np.zeros_like(fine), fine * math.sin(e)], axis=1)
    t, xyz = _fly(v0, g, rho)
    last = np.where(np.isfinite(xyz[..., 0]), xyz[..., 0], -np.inf).max(axis=1)
    j = int(np.argmin(np.abs(last - LOB_M)))
    path = xyz[j]
    n = int(np.isfinite(path[:, 2]).sum())
    return dict(speed_m_s=float(fine[j]), t=t[:n], path=path[:n])


def earliest(body, place, t, path, grip, g, start=state.REACTION_S):
    """The first sample of the path the player at `place`, facing back toward the launcher, can reach, starting `start` after the
    release: (index or None, its mode)."""
    st = state.State(p=(0.0, 0.0), posture='ready', facing=math.pi)
    rel = path - np.array([place[0], place[1], 0.0])
    tt = t - start
    ok, mode = _along(body, st, tt, rel, grip, g)
    ok &= tt >= 0
    if not ok.any():
        return None, None
    k = int(np.argmax(ok))
    return k, (mode[k] if mode is not None and mode.shape else None)


def read_window(body, place, t, path, k, grip, g):
    """How much later than the reaction the player could start and still take the ball at sample k or before."""
    best = 0.0
    for d in DELAYS_S:
        j, _ = earliest(body, place, t[:k + 1], path[:k + 1], grip, g, start=state.REACTION_S + d)
        if j is None:
            break
        best = d
    return best


def crowd_radius(body, t, path, t_best, grip, g, directions=CROWD_DIRECTIONS, far_m=40.0):
    """How far from the first touch's spot (the floor under the ball at the earliest touch), in each direction, a player
    can start and still touch the ball within CONTEST_S of that touch (bisection on the distance, to 0.1 m)."""
    L = path[int(np.argmin(np.abs(t - t_best))), :2]
    out = []
    for a in directions:
        u = np.array([math.cos(math.radians(a)), math.sin(math.radians(a))])
        lo, hi = 0.0, far_m
        k, _ = earliest(body, L, t, path, grip, g)
        if k is None or t[k] > t_best + CONTEST_S:
            out.append(0.0)
            continue
        while hi - lo > 0.1:
            mid = 0.5 * (lo + hi)
            k, _ = earliest(body, L + mid * u, t, path, grip, g)
            lo, hi = (mid, hi) if (k is not None and t[k] <= t_best + CONTEST_S) else (lo, mid)
        out.append(round(lo, 1))
    return out


def contest(g, rho, grip, elevation_deg, b_place, a_body=None, b_body=None):
    a_body = a_body or state.BODIES['dancer']
    b_body = b_body or state.BODIES['dancer']
    shot = lob(g, rho, elevation_deg)
    t, path = shot['t'], shot['path']
    L = np.array([LOB_M, 0.0])
    bq = L + np.array(B_PLACES[b_place])
    kb, mb = earliest(b_body, bq, t, path, grip, g)
    tb = t[kb] if kb is not None else math.inf
    slack_b = read_window(b_body, bq, t, path, kb, grip, g) if kb is not None else None
    cells = []
    for dx in GRID_M:
        for dy in GRID_M:
            q = L + np.array([dx, dy])
            if np.hypot(*(q - bq)) < 0.5:
                continue
            ka, ma = earliest(a_body, q, t, path, grip, g)
            ta = t[ka] if ka is not None else math.inf
            if not (math.isfinite(ta) or math.isfinite(tb)):
                who = 'neither'
            elif abs(ta - tb) < CONTEST_S:
                who = 'contested'
            else:
                who = 'A' if ta < tb else 'B'
            row = dict(dx=dx, dy=dy, who=who, a_s=None if ka is None else round(float(ta), 2))
            if ka is not None:
                row['mode'] = reach.MODES[int(ma)]
                row['slack_s'] = read_window(a_body, q, t, path, ka, grip, g)
            if who == 'contested' and slack_b is not None:
                diff = row['slack_s'] - slack_b
                row['position'] = 'even' if abs(diff) < CONTEST_S else ('A' if diff > 0 else 'B')
            cells.append(row)
    n = len(cells)
    reached = [c for c in cells if 'slack_s' in c]
    t_best = min([tb] + [c['a_s'] for c in cells if c['a_s'] is not None])
    radius = crowd_radius(a_body, t, path, t_best, grip, g)
    contested = [c for c in cells if c['who'] == 'contested']
    return dict(speed_m_s=round(shot['speed_m_s'], 2), flight_s=round(float(t[-1]), 2),
                apex_m=round(float(np.nanmax(path[:, 2])), 2),
                b=dict(touch_s=None if kb is None else round(float(tb), 2),
                       mode=None if mb is None else reach.MODES[int(mb)], slack_s=slack_b),
                shares={w: round(sum(c['who'] == w for c in cells) / n, 3) for w in ('A', 'B', 'contested', 'neither')},
                earliest_touch_s=round(float(t_best), 2),
                first_touch_spot_m=[round(float(v), 1) for v in path[int(np.argmin(np.abs(t - t_best))), :2]],
                crowd_radius_m=dict(by_direction_deg=dict(zip(CROWD_DIRECTIONS, radius)),
                                    median=round(float(np.median(radius)), 1)),
                crowd_area_m2=round(float(np.mean(np.square(radius)) * math.pi), 1),
                touch_by_jump=(round(sum(c['mode'] == 'jump' for c in reached) / len(reached), 2) if reached else None),
                median_slack_s=(round(float(np.median([c['slack_s'] for c in reached])), 2) if reached else None),
                position_in_contests=({k: round(sum(c.get('position') == k for c in contested) / len(contested), 2)
                                       for k in ('A', 'B', 'even')} if contested else None),
                cells=cells)


def _case(args):
    world, g, rho, grip, b_place = args
    by_el = {el: contest(g, rho, grip, el, b_place) for el in ELEVATIONS_DEG}
    cases = [dict(world=world, grip=grip, elevation_deg=el, b_place=b_place, **{k: v for k, v in c.items() if k != 'cells'})
             for el, c in by_el.items()]
    # the launcher's control: where the elevation alone changes who gets the ball
    whos = np.array([[c['who'] for c in by_el[el]['cells']] for el in ELEVATIONS_DEG])
    control = dict(world=world, grip=grip, b_place=b_place, share=round(float((whos != whos[0]).any(axis=0).mean()), 3))
    return cases, control


def screen(worlds=None, grips=None, workers=WORKERS):
    worlds = worlds or [(w, g, state.RHO_MOON if w == 'moon' else state.RHO_EARTH) for w, g in state.WORLDS]
    grips = grips or state.GRIPS
    jobs = [(w, g, rho, grip, b) for w, g, rho in worlds for grip in grips for b in B_PLACES]
    if workers > 1:
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context('fork')) as pool:
            done = list(pool.map(_case, jobs))
    else:
        done = [_case(j) for j in jobs]
    return dict(ball=BALL, lob_m=LOB_M, release_m=RELEASE_M, reaction_s=state.REACTION_S, contest_s=CONTEST_S,
                cases=[c for cases, _ in done for c in cases], launcher_control=[c for _, c in done])
