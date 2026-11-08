"""Parkour moves for a lunar obstacle game: the surface moves of surfaces.py written in the athletic vocabulary's form,
joined to its units, and set against the floor on the Moon beside Earth.

    python -m humanities.sport.parkour     # writes humanities/sport/results/parkour.json (about a minute)

A paintball game's arena is built from obstacles, and on the Moon obstacles are the fast way to move. The floor's
grip caps a stance's change of velocity at mu g on average, so a cut at speed is a visible jump; a push off a wall or
a box presses the foot or the hand into it, so its grip grows with the push. surfaces.py times the moves that use
surfaces (vaults, a foot's push off a wall, running up a wall, catching a ledge and climbing it, a precision landing,
a drop and roll, a bar's swing and release, a slide into cover and pulling down behind it) and routes a player across
an arena with them. This module writes those moves as units (vocabulary/parkour_units.json) in the form of the
athletic vocabulary (vocabulary.py), with the tell, commitment and recovery a contest reads, adds one of its own (the
rise from low that a slide or a pull-down ends in), and evaluates every unit for both bodies on both worlds at each of
the floor's grips. Two postures join the athletic ones: hang (from a ledge or a bar) and down (low behind cover).

The product also holds what the moves change: the wall's push against the floor's cut, a sustained run round a curved
wall, and the quickest ways from cover to cover across two arenas (surfaces.quickest_route) beside floor travel alone
and an open floor. The surface moves' assumptions (the arms' pull, the legs' yield, the walls' grip, the clearances)
are stated in surfaces.py; the vocabulary's in vocabulary.py and state.py.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from humanities.dance import envelope as dance
from humanities.sport import envelope as env
from humanities.sport import state, surfaces
from humanities.sport import vocabulary as V
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
UNITS_JSON = HERE / 'vocabulary' / 'parkour_units.json'
OUT = HERE / 'results' / 'parkour.json'
SCHEMA = 'terluna.humanities.sport-parkour/1'
FILES = ('humanities/sport/parkour.py', 'humanities/sport/vocabulary/parkour_units.json', 'humanities/sport/surfaces.py',
         'humanities/sport/vocabulary.py', 'humanities/sport/vocabulary/units.json', 'humanities/sport/reach.py',
         'humanities/sport/state.py', 'humanities/sport/envelope.py', 'humanities/dance/envelope.py',
         'humanities/dance/body/model.py', 'humanities/dance/body/pose.py', 'habitation/low_gravity/movement.py')
INPUTS = (dance.WINDS, dance.POSES)
UNITS = surfaces.PARKOUR_UNITS + ('rise_from_low',)
ARENAS = {'a 1 m barrier across the way': (surfaces.BARRIER, surfaces.BARRIER_START, surfaces.BARRIER_GOAL),
          'bunkers, a low box, a 3 m wall and a bar': (surfaces.EXAMPLE, surfaces.EXAMPLE_START, surfaces.EXAMPLE_GOAL)}
CURVED_WALL_RADII_M = (3.0, 4.7, 10.0)
WORKERS = 3


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16]


def r(x, n=3):
    return None if x is None or not math.isfinite(float(x)) else round(float(x), n)


def load():
    return json.loads(UNITS_JSON.read_text())


def physics(body, g, mu):
    """The vocabulary's names for one body, world and grip, with the surface moves added as functions of the move's
    own parameters."""
    ns = V.physics(body, g, mu)
    ns.update(
        vault=lambda h, length, v, kind: surfaces.vault(body, h, length, v, g, kind=kind),
        tic_tac=lambda v, angle: surfaces.tic_tac(body, v, angle, g, grip=mu),
        floor_cut=lambda v, turn: surfaces.floor_cut(v, turn, mu, g),
        wall_run=lambda v, steps, push: surfaces.wall_run(body, v, g, steps=steps, push=push),
        cat_leap=lambda v_down: surfaces.cat_leap(body, g, v_down=v_down),
        climb_up=lambda: surfaces.climb_up(body, g),
        precision_jump=lambda gap: surfaces.precision_jump(body, gap, g),
        drop_roll=lambda level: surfaces.drop_roll(body, level, g, grip=mu),
        lache=lambda v: surfaces.lache(body, g, catch_m_s=v),
        slide=lambda v: surfaces.slide(body, v, g, grip=mu),
        pull_down=lambda depth: surfaces.pull_down(body, depth, g),
        slide_com=surfaces.SLIDE_COM_M * body.stature_m / env.STATURE_M, int=int, bool=bool)
    return ns


def times(unit, body, g, mu, given=None):
    """A unit's phases, three times and figures for one body, world and grip (vocabulary.times with these names)."""
    names = physics(body, g, mu)
    names.update({k: v['default'] for k, v in unit.get('params', {}).items()})
    names.update(given or {})
    for k, expr in unit.get('let', {}).items():
        names[k] = eval(expr, {'__builtins__': {}}, names)
    phases = [dict(kind=p['kind'], seconds=V.evaluate(p['duration'], names)) for p in unit['phases']]
    out = dict(duration_s=sum(p['seconds'] for p in phases), phases=phases, tell_s=V.evaluate(unit['tell_s'], names),
               commits_s=V.evaluate(unit['commits_s'], names), recover_s=V.evaluate(unit['recover_s'], names))
    out['figures'] = {k: V.evaluate(v, names) for k, v in unit.get('figures', {}).items()}
    return out


def figures(unit):
    out = {}
    for key, body in state.BODIES.items():
        out[key] = {}
        for world, g in state.WORLDS:
            out[key][world] = {}
            for mu in state.GRIPS:
                t = times(unit, body, g, mu)
                row = dict(duration_s=r(t['duration_s']), tell_s=r(t['tell_s']), commits_s=r(t['commits_s']),
                           recover_s=r(t['recover_s']), phases=[[p['kind'], r(p['seconds'])] for p in t['phases']])
                row.update({k: r(v) for k, v in t['figures'].items()})
                out[key][world][str(mu)] = row
    return out


def joined_graph():
    """Transitions over the athletic units and the parkour units together (vocabulary.transitions)."""
    units = V.load()['units'] + load()['units']
    return units, V.transitions(units)


def wall_against_floor():
    """A foot's push off a wall against the floor's cut through the same turn, at 3 m/s meeting the wall 45 degrees
    from its normal, for the dance body."""
    body = state.BODIES['dancer']
    out = {}
    for world, g in state.WORLDS:
        out[world] = {}
        for mu in state.GRIPS:
            w = surfaces.tic_tac(body, 3.0, 45.0, g, grip=mu)
            f = w['floor_cut']
            out[world][str(mu)] = dict(
                wall_change_m_s=r(w['kept']['velocity_change_m_s'], 2), wall_largest_change_m_s=r(w['most']['velocity_change_m_s'], 2),
                turn_deg=r(w['kept']['turn_deg'], 1), largest_turn_deg=r(w['most']['turn_deg'], 1), contact_s=r(w['contact_s']),
                floor_one_stance_m_s=r(f['one_stance_m_s'], 2), floor_cut_bounce_rise_m=r(f['bounce_rise_m'], 2),
                floor_cut_total_s=r(f['total_s'], 2),
                wall_over_floor=r(w['kept']['velocity_change_m_s'] / f['one_stance_m_s'], 1))
    return out


def curved_wall():
    """The least speed at which friction holds a runner's weight on a curved wall of radius R: mu v^2 / R >= g."""
    return {world: {str(R): {str(mu): r(surfaces.wall_of_death_m_s(g, R, mu), 2) for mu in state.GRIPS}
                    for R in CURVED_WALL_RADII_M} for world, g in state.WORLDS}


def _route(args):
    name, world, mu, body_key = args
    arena, start, goal = ARENAS[name]
    g = dict(state.WORLDS)[world]
    out = surfaces.quickest_route(state.BODIES[body_key], state.State(p=start, posture='ready'), goal, arena, mu, g)
    best = out['best']
    row = dict(arena=name, world=world, grip=mu, body=body_key, open_floor_s=r(out['open_floor']['seconds'], 2),
               floor_round_s=r(out['floor_only']['seconds'], 2) if out['floor_only']['ok'] else None,
               with_surface_moves_s=r(best['seconds'], 2),
               moves=[m['unit'] for m in best['moves'] if m['unit'] != 'floor'])
    return row


def routes(workers=WORKERS):
    jobs = [(name, world, mu, 'dancer') for name in ARENAS for world, _ in state.WORLDS for mu in state.GRIPS]
    jobs += [(name, world, state.MU, 'larger') for name in ARENAS for world, _ in state.WORLDS]
    if workers > 1:
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context('fork')) as pool:
            return list(pool.map(_route, jobs))
    return [_route(j) for j in jobs]


def results():
    units, graph = joined_graph()
    parkour = load()['units']
    starts = V.reachable(graph, 'ready_crouch')
    back = {k: [a for a in graph if k in graph[a]] for k in graph}
    return dict(
        schema=SCHEMA,
        producer=dict(domain='humanities', files={f: digest(f) for f in FILES},
                      constants=constants_used([ROOT / f for f in FILES]),
                      inputs={str(p.relative_to(ROOT)): digest(p.relative_to(ROOT)) for p in INPUTS}),
        evidence='Screens of the surface moves (surfaces.py: pushes on the jump\'s force-velocity profile, friction '
                 'on surfaces, an assumed arm pull, the reach engine\'s floor law) for two representative bodies at the '
                 'shared gravity and Earth\'s, written as vocabulary units and joined to the athletic units; the routes '
                 'come from a beam search over up to three surface moves. Nothing was solved on the whole body, '
                 'measured in low gravity or run.',
        grips=list(state.GRIPS), bodies={k: b.name for k, b in state.BODIES.items()},
        units=[dict(id=u['id'], name=u['name'], entry=u['entry'], exit=u['exit'], successors=graph[u['id']],
                    by_body=figures(u)) for u in parkour],
        transitions=dict(units=len(units), parkour_units=len(parkour),
                         reach_all_from_ready_crouch=len(starts) == len(units),
                         all_lead_back_to_ready_crouch=len(V.reachable(back, 'ready_crouch')) == len(units)),
        wall_against_floor=wall_against_floor(), curved_wall_least_m_s=curved_wall(), routes=routes())


def main() -> int:
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
