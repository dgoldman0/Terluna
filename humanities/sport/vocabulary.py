"""The athletic vocabulary: the moves a lunar sport can be made of, written as units like the dance's.

    python -m humanities.sport.vocabulary     # writes humanities/sport/results/vocabulary.json at once

The second round of the sport study (8 November 2026) writes the moves before any game. Each unit is a repeatable
action of play in the dance grammar's form (humanities/dance/grammar.py): it enters from a state and leaves in another,
through phases timed by the envelope's mechanics. A state here is a posture (state.POSTURES), a momentum (still, slow
or travelling) and the ball (held or free). To the dance's form each unit adds the three times a contest turns on:
- the tell: how long an opponent can see the move coming before it can no longer be changed;
- the commitment: how long it then holds the body to its course (an arc, a push, a contact);
- the recovery: how long after it ends before the next quick action can start.
Every time is evaluated for both bodies, on the Moon and on Earth, at each of the floor's grips (state.GRIPS). A unit
can follow another wherever the one's exit state is an entry state of the other, which gives the vocabulary's
transitions.

Expressions in vocabulary/units.json read the names physics() binds for one body, world and grip: the body's
dimensions, the take-off speeds, and the envelope's mechanics as functions (a dip no faster than g, a rise that stops
no faster than g, the lean's inverted pendulum, the grip on a stance's vertical impulse, a throw's moment balance).
A unit's `let` names are evaluated first, in order. Two assumptions are this module's own: a throwing arm or a
kicking leg winds up in 0.20 s before its push or swing (the muscles' clock, either world), and a kicked ball leaves
at 1.2 times the foot's speed (instep kicks, from general knowledge).
"""
from __future__ import annotations

import hashlib
import json
import math
from itertools import product
from pathlib import Path

from humanities.dance import envelope as dance
from humanities.sport import envelope as env
from humanities.sport import state
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
UNITS_JSON = HERE / 'vocabulary' / 'units.json'
OUT = HERE / 'results' / 'vocabulary.json'
SCHEMA = 'terluna.humanities.sport-vocabulary/1'
FILES = ('humanities/sport/vocabulary.py', 'humanities/sport/vocabulary/units.json', 'humanities/sport/state.py',
         'humanities/sport/envelope.py', 'humanities/dance/envelope.py', 'humanities/dance/body/model.py',
         'humanities/dance/body/pose.py', 'habitation/low_gravity/movement.py')
INPUTS = (dance.WINDS, dance.POSES)

WINDUP_S = 0.20                           # a throwing arm or a kicking leg winds up (assumed)
KICK_BALL_RATIO = 1.2                     # a kicked ball's speed over the foot's (instep kicks; general knowledge)
LEG_M = 0.9                               # the kicking leg, hip to foot (envelope: 15 rad/s gives about 13 m/s)
KICK_CONTACT_M = 0.15                     # the ball's centre above the floor at a kick
CHEST_LEVER_M = env.AIR_LEVERS_M['chest pass, near the centre of mass']
HANDS_ABOVE_COM_M = 0.40                  # a shove's hands above the crouched centre of mass (assumed)
STILL_M_S = 0.15                          # below this a body counts as still
MOMENTA = ('still', 'slow', 'travelling')
BALLS = ('held', 'free')


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16]


def r(x, n=3):
    return round(float(x), n)


def physics(body, g, mu):
    """The names a unit's expressions can read, for one body on one world at one grip."""
    v_full = body.take_off_m_s(g)
    v_crouch = body.take_off_m_s(g, body.crouch_push_m)
    a_up = v_full ** 2 / (2 * body.push_off_m)
    m = body.mass_kg
    h, hr = body.com_height_m, body.ready_com_m
    ball_kg = env.BALLS[env.THROW_BALL][0]
    kick_kg = env.BALLS['football (size 5)'][0]
    _, leg_I = env.leg_inertia_about_hip()
    sl, sm = body.stature_m / env.STATURE_M, m / env.MASS_KG
    leg_I *= sm * sl * sl

    def takeoff(rise):
        return math.sqrt(2 * g * max(rise, 0.0))

    def rise_of(vz):
        return vz * vz / (2 * g)

    def airtime(vz):
        return 2 * vz / g

    def push_s(vz, stroke):
        """A push over `stroke` at the mean speed vz/2."""
        return 2 * stroke / max(vz, 1e-6)

    def dip_s(depth):
        return env.dip_s(depth, g)

    def rise_s(height):
        """Raising the centre of mass and arriving at rest: quick at the push, gravity's at the stop."""
        return env.rise_and_stop_s(height, a_up, g)

    def lean_s(distance, reach, height):
        """The inverted pendulum about the support's edge, from rest (envelope.lean_start, unrounded)."""
        return math.sqrt(height / g) * math.acosh(1 + distance / reach)

    def lean_v(distance, reach, height):
        tau = math.sqrt(height / g)
        return reach / tau * math.sinh(math.acosh(1 + distance / reach))

    def start_s(distance, lead, height, cap=state.V_MAX_M_S, dt=5e-4):
        """Covering `distance` from rest: x'' = min(mu g, (g/height)(lead + x), the legs' push at that speed), the lean
        about the support's edge `lead` behind the centre of mass until it gives the grip's limit, with the speed
        capped at `cap` (stepped at dt, semi-implicit)."""
        x = v = t = 0.0
        while x < distance:
            a = min(mu * g, g / height * (lead + x), state.muscle_push(v))
            v = min(v + a * dt, cap)
            x += v * dt
            t += dt
        return t - (x - distance) / max(v, 1e-9)

    def jump(posture, vz, part):
        """A take-off at vz from a posture (Body.jump_timing): part 0 the dip in m, 1 the dip's seconds, 2 the push's."""
        out = body.jump_timing(posture, vz, g)
        return math.inf if out is None else out[part]

    def grip_s(dv):
        return dv / (mu * g)

    def grip_m(dv):
        return dv * dv / (2 * mu * g)

    def stride_dv(rise):
        """The most one stride's stance changes the travel: grip on its vertical impulse, landing plus take-off."""
        return mu * (2 * takeoff(rise) + g * env.STRIDE_S)

    def cut_w(speed, angle_deg):
        """The vertical speed each way a cut through one stance needs (envelope.cut_bounce)."""
        dv = 2 * speed * math.sin(math.radians(angle_deg) / 2)
        return max(0.0, (dv / mu - g * env.STRIDE_S) / 2)

    def landing_s(v_down):
        return 2 * env.STROKE_M / max(v_down, 1e-6)

    def landing_shed(speed, v_down):
        return min(speed, mu * (v_down + g * landing_s(v_down)))

    def landing_push(v_down):
        """The floor's mean push in a landing absorbed over the stroke, per kilogram."""
        return v_down * v_down / (2 * env.STROKE_M) + g

    def pivot_s(angle_deg, half_stance):
        """Turning the body about the vertical on both feet, half the time speeding the turn and half stopping it, with
        the feet's friction couple mu m g half_stance against the body's inertia about its long axis."""
        tau = mu * m * g * half_stance
        return 2 * math.sqrt(math.radians(angle_deg) * body.inertia_kg_m2['twist'] / tau)

    def hold_m_s(reach, lever, push, kg=ball_kg, t=env.THROW_S):
        """The fastest release (or catch) whose hand force the support holds: the moment balance of envelope.throws
        (floor push times reach against hand force times lever), the grip (hand force at most mu times the push) and
        the arm's own limit."""
        by_moment = push * reach * t / (kg * lever) if lever > 0 else math.inf
        return min(by_moment, mu * push * t / kg, state.RELEASE_MAX_M_S)

    def turn_deg(dL, inertia, t_left):
        """The attitude an angular momentum dL leaves by the landing, t_left later."""
        return math.degrees(dL / inertia * t_left)

    def topple_s(a0, a1):
        return dance.topple_seconds(a0, a1, g, h)

    toss_a = env.TOSS_FORCE_N / m - g
    toss_vz = math.sqrt(2 * toss_a * env.TOSS_STROKE_M) if toss_a > 0 else 0.0
    foot_m_s = env.KICK_RAD_S * LEG_M
    ns = dict(
        g=g, g_earth=env.G_E, mu=mu, m=m, h=h, h_ready=hr, half_stance=body.ready_half_stance_m, side=body.side_reach_m,
        fore_aft=body.fore_aft_reach_m, split=body.split_reach_m, push_off=body.push_off_m,
        crouch_push=body.crouch_push_m, dip=body.dip_m, reach_up=body.reach_up_m, reach_out=body.reach_out_m,
        reach_low=body.reach_low_m, lever=body.release_above_com_m, chest_lever=CHEST_LEVER_M,
        hands_above=HANDS_ABOVE_COM_M, v_full=v_full, v_crouch=v_crouch, v_max=state.V_MAX_M_S,
        walk_max=math.sqrt(g * h), step_s=state.STEP_S, stride_s=env.STRIDE_S, stroke=env.STROKE_M,
        throw_s=env.THROW_S, catch_s=state.CATCH_S, windup_s=WINDUP_S, ball_kg=ball_kg, kick_kg=kick_kg,
        I_layout=body.inertia_kg_m2['layout'], I_tuck=body.inertia_kg_m2['tuck'], leg_I=leg_I,
        kick_rad_s=env.KICK_RAD_S, kick_build_s=env.KICK_BUILD_S, foot_m_s=foot_m_s,
        kick_ball_m_s=KICK_BALL_RATIO * foot_m_s, kick_lever=h - KICK_CONTACT_M,
        apex_descent_s=env.APEX_TAKE_OFF_M_S / g, toss_vz=toss_vz, toss_stroke=env.TOSS_STROKE_M,
        takeoff=takeoff, rise_of=rise_of, airtime=airtime, push_s=push_s, dip_s=dip_s, rise_s=rise_s, lean_s=lean_s,
        lean_v=lean_v, start_s=start_s, jump=jump, grip_s=grip_s, grip_m=grip_m, stride_dv=stride_dv, cut_w=cut_w, landing_s=landing_s,
        landing_shed=landing_shed, landing_push=landing_push, pivot_s=pivot_s, hold_m_s=hold_m_s, turn_deg=turn_deg,
        topple_s=topple_s, min=min, max=max, sqrt=math.sqrt, sin=math.sin, cos=math.cos, radians=math.radians,
        pi=math.pi)
    return ns


def load():
    return json.loads(UNITS_JSON.read_text())


def evaluate(expr, names):
    if isinstance(expr, (int, float)):
        return float(expr)
    return float(eval(expr, {'__builtins__': {}}, names))


def times(unit, body, g, mu, given=None):
    """A unit's phases, its three times and its own figures for one body, world and grip."""
    names = physics(body, g, mu)
    names.update({k: v['default'] for k, v in unit.get('params', {}).items()})
    names.update(given or {})
    for k, expr in unit.get('let', {}).items():
        names[k] = evaluate(expr, names)
    phases = [dict(kind=p['kind'], seconds=evaluate(p['duration'], names), overlaps=bool(p.get('overlaps')))
              for p in unit['phases']]
    out = dict(duration_s=sum(p['seconds'] for p in phases if not p['overlaps']), phases=phases,
               tell_s=evaluate(unit['tell_s'], names), commits_s=evaluate(unit['commits_s'], names),
               recover_s=evaluate(unit['recover_s'], names))
    out['figures'] = {k: evaluate(v, names) for k, v in unit.get('figures', {}).items()}
    return out


# --- states and transitions ------------------------------------------------------------------------------------------

def _listed(x, every):
    if x in ('any', None):
        return list(every)
    return list(x) if isinstance(x, list) else [x]


def entry_states(unit):
    e = unit['entry']
    return list(product(_listed(e['posture'], state.POSTURES), _listed(e.get('momentum'), MOMENTA),
                        _listed(e.get('ball'), BALLS)))


def exit_state(unit, entry):
    x = unit['exit']
    posture, momentum, ball = entry
    return (posture if x['posture'] == 'same' else x['posture'],
            momentum if x.get('momentum', 'same') == 'same' else x['momentum'],
            ball if x.get('ball', 'same') == 'same' else x['ball'])


def transitions(units):
    """Each unit's successors: the units with an entry state among its exit states."""
    entries = {u['id']: set(entry_states(u)) for u in units}
    out = {}
    for u in units:
        exits = {exit_state(u, s) for s in entries[u['id']]}
        out[u['id']] = [v['id'] for v in units if exits & entries[v['id']]]
    return out


def reachable(graph, start):
    seen, todo = {start}, [start]
    while todo:
        for nxt in graph[todo.pop()]:
            if nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    return seen


# --- the product -----------------------------------------------------------------------------------------------------

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


def readings(units):
    """Moon beside Earth for the dance body at the central grip: the times a contest reads, and their ratio."""
    body = state.BODIES['dancer']
    rows = []
    for u in units:
        m, e = (times(u, body, g, state.MU) for _, g in state.WORLDS)
        rows.append(dict(unit=u['id'], **{f'{k}_moon': r(m[k], 2) for k in ('duration_s', 'tell_s', 'commits_s')},
                         **{f'{k}_earth': r(e[k], 2) for k in ('duration_s', 'tell_s', 'commits_s')}))
    return dict(body='dancer', grip=state.MU, rows=rows,
                no_tell=[u['id'] for u in units if times(u, body, state.G_M, state.MU)['tell_s'] == 0],
                committed_over_1_s_on_the_moon=[u['id'] for u in units
                                                if times(u, body, state.G_M, state.MU)['commits_s'] > 1.0])


def results():
    vocab = load()
    units = vocab['units']
    graph = transitions(units)
    return dict(
        schema=SCHEMA,
        producer=dict(domain='humanities', files={f: digest(f) for f in FILES},
                      constants=constants_used([ROOT / f for f in FILES]),
                      inputs={str(p.relative_to(ROOT)): digest(p.relative_to(ROOT)) for p in INPUTS}),
        evidence='Physics screens of each move\'s timing from the envelope\'s closed forms, for two representative '
                 'bodies at the shared gravity, with every assumption stated beside its value (units.json, '
                 'vocabulary.py, state.py, envelope.py). The units are drawn from the practice of court and field '
                 'sports, athletics, gymnastics and the martial arts, and designed for the Moon; none is measured in '
                 'low gravity or played.',
        grips=list(state.GRIPS), bodies={k: b.name for k, b in state.BODIES.items()},
        units=[dict(id=u['id'], name=u['name'], family=u['family'], entry=u['entry'], exit=u['exit'],
                    successors=graph[u['id']], by_body=figures(u)) for u in units],
        transitions=dict(units=len(units), edges=sum(len(v) for v in graph.values()),
                         reach_all_from_ready_crouch=len(reachable(graph, 'ready_crouch')) == len(units)),
        readings=readings(units))


def main() -> int:
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
