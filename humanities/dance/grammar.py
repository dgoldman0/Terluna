"""The grammar: how the vocabulary's units join into phrases a dancer could perform on the Moon.

    python -m humanities.dance.grammar      # composes the example phrases into humanities/dance/results/phrases.json

A phrase is a sequence of placed units. Each unit enters from a state (which foot or feet bear the weight, how fast
and which way the body travels, how fast it turns, where the sleeves are) and leaves in another. The grammar places a
unit only where its entry fits the state the previous unit left; where it does not, it inserts the linking step a
dancer would use (a weight change, strides that build or spend speed, a braking turn, a gather of the sleeves), each
sized by the envelope's physics. Durations come from that physics too: flights from their rise, pushes and landings
from their speeds, descents from gravity, falls from the inverted pendulum, turn preparations from the floor's grip,
waits for fabric from its settling speed.

Two clocks run through every phrase. Gravity's clock (falls, flights, fabric, the pendulum of a swinging leg) runs
2.45 times slower on the Moon than on Earth; the muscles' clock (an arm's whip, the feet's patter, the push of a
take-off) keeps Earth's tempo. The grammar keeps both: quick articulation can sit inside a slow arc, and the body can
rise as fast as on Earth but never sink faster than the Moon pulls it.

Devices from choreographic practice (repetition, mirroring, augmentation, diminution, accumulation, retrograde of
order) act on motifs; each result is placed and checked again. The evaluation reports what a choreographer would look
at (length, variety, contrast of effort and level, time aloft, the lunar share, travel and curvature, stops) as
measures for comparison; the author and dancers judge the art.
"""
from __future__ import annotations

import copy
import json
import math
import random
from dataclasses import dataclass, field, replace
from pathlib import Path

from humanities.dance import envelope as env

HERE = Path(__file__).resolve().parent
VOCAB = HERE / 'vocabulary'
OUT = HERE / 'results' / 'phrases.json'
SCHEMA = 'terluna.humanities.dance-phrases/1'
PRODUCER_FILES = ('humanities/dance/grammar.py', 'humanities/dance/envelope.py', 'humanities/dance/vocabulary/units.json',
                  'humanities/dance/vocabulary/poses.json', 'humanities/dance/body/model.py')

G = env.G_M
MU = env.FRICTION
COM_H = None            # set lazily from the body model


def _com_height():
    global COM_H
    if COM_H is None:
        COM_H = float(env.whole_body_inertia()[1][2])
    return COM_H


# --- physics available to the vocabulary's duration expressions ---------------------------------------------------

def takeoff(rise_m):
    return math.sqrt(2 * G * rise_m)


def flight_s(rise_m):
    return 2 * math.sqrt(2 * rise_m / G)


def push_s(v):
    """A push from a crouch whose depth grows with the speed it gives (ankles alone for a small hop, the whole leg
    for a full jump), covered at the mean speed v/2."""
    depth = min(env.H_PO_M, 0.03 + 0.11 * v)
    return 2 * depth / max(v, 1e-3)


def absorb_s(v, depth):
    return 2 * depth / max(v, 1e-3)


def descend_s(depth, share=0.7):
    """Lowering the centre of mass by depth and stopping, with the downward acceleration at most share * g: the body
    can never sink faster than gravity pulls it."""
    return 2 * math.sqrt(depth / (share * G))


def topple_s(a0, a1):
    return env.topple_seconds(a0, a1, G, _com_height())


def fall_speed(angle_deg, a0=2.0):
    """Speed of the centre of mass after toppling from a0 to angle_deg about the ankle."""
    L = _com_height()
    w2 = 2 * G / L * (math.cos(math.radians(a0)) - math.cos(math.radians(angle_deg)))
    return L * math.sqrt(max(w2, 0.0))


def fall_grip(angle_deg, a0=2.0):
    """The most grip a fall as one piece asks of the floor on its way from a0 to angle_deg: the ratio of the floor's
    sideways push to its upward push, for the whole body pivoting about its foot from rest. The body's own moment of
    inertia about its centre of mass slows the fall's sideways acceleration, so a dancer needs less grip than a point
    mass on a rod (tan(angle)); the ratio is the same at any gravity."""
    L = _com_height()
    k2 = env.whole_body_inertia()[2][1, 1] / env.MASS_KG
    r = L * L / (L * L + k2)
    c0 = math.cos(math.radians(a0))
    worst = 0.0
    for i in range(1, 201):
        th = math.radians(a0 + (angle_deg - a0) * i / 200)
        s, c = math.sin(th), math.cos(th)
        worst = max(worst, r * s * (3 * c - 2 * c0) / (1 - r * (s * s + 2 * c * (c0 - c))))
    return worst


def shift_s(across_m, half_stance_m=0.10, half_foot_m=0.04):
    """Carrying the centre of mass from between the feet onto one foot without a step, as fast as the floor allows. The
    body is an inverted pendulum pushed sideways from the far foot's outer edge and braked from the near foot's outer
    edge (x'' = (g/h)(x - p), p at the support's edges), coming to rest over the standing foot. Gravity sets the pace:
    a 10 cm shift takes 1.7 s on the Moon and 0.7 s on Earth."""
    w = math.sqrt(G / _com_height())
    edge = half_stance_m + half_foot_m

    def meet(a):          # pushed for a/w, braked for b/w: positions where the two arcs meet
        b = math.asinh(edge * math.sinh(a) / (edge - across_m))
        return edge * (math.cosh(a) - 1) - (edge - (edge - across_m) * math.cosh(b)), b

    lo, hi = 0.0, 6.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if meet(mid)[0] < 0 else (lo, mid)
    return (lo + meet(lo)[1]) / w


def prep_s(rev_per_s, arms_open=0):
    """Preparation that builds a turn's angular momentum with the floor's grip, for a turn with the arms gathered (the
    default) or open in second (arms_open=1), which nearly triples the body's inertia about the vertical."""
    I = env.whole_body_inertia(arms_side=bool(arms_open))[2][2, 2]
    tau = MU * env.MASS_KG * G * env.FEET_APART_M / 2
    return I * 2 * math.pi * rev_per_s / tau


def stride_dv(stride_rise_m=0.08):
    """Speed one bounding stride can add or take away. Its stance must reverse the vertical speed v_z as well as
    carry the weight, so the floor presses with m (g + 2 v_z / T) on average and friction allows mu (g T + 2 v_z):
    about 0.8 m/s a stride for the orbit run's 8 cm strides, where a walking step's mu g T gives a third of that."""
    vz = takeoff(stride_rise_m)
    return MU * (G * push_s(vz) + 2 * vz)


PHYSICS = dict(takeoff=takeoff, flight_s=flight_s, push_s=push_s, absorb_s=absorb_s, descend_s=descend_s,
               topple_s=topple_s, fall_speed=fall_speed, prep_s=prep_s, shift_s=shift_s, max=max, min=min,
               sqrt=math.sqrt)


def evaluate(expr, params):
    if isinstance(expr, (int, float)):
        return float(expr)
    return float(eval(expr, {'__builtins__': {}}, {**PHYSICS, **params}))


# --- vocabulary ------------------------------------------------------------------------------------------------------

def load_vocabulary():
    units = json.loads((VOCAB / 'units.json').read_text())
    poses = json.loads((VOCAB / 'poses.json').read_text())
    return {u['id']: u for u in units['units']}, poses


UNITS, POSES = load_vocabulary()

LUNAR_NATIVE = {uid for uid, u in UNITS.items() if u['lunar'].startswith('native')}


# --- state ------------------------------------------------------------------------------------------------------------

@dataclass
class State:
    support: str = 'both'          # both, L, R, none, floor
    speed: float = 0.0             # horizontal speed, m/s
    heading: float = 0.0           # direction of travel, degrees anticlockwise from the stage's +x
    spin: float = 0.0              # revolutions a second about the vertical
    sleeves: str = 'gathered'      # gathered, aloft, trailing
    sleeves_down_at: float = 0.0   # when cast sleeves are back within reach
    t: float = 0.0
    x: float = 0.0
    y: float = 0.0
    aloft: bool = False            # the unit ended in the air: the next stance is a landing
    last_kind: str = ''            # the kind of the phase the unit ended in (a shift ends on both feet)


@dataclass
class Placed:
    unit: str
    side: str                       # the gesture side G
    params: dict
    t0: float
    duration: float
    phases: list
    entry: State
    exit: State
    connector: bool = False
    note: str = ''
    dropped_first: bool = False     # its first stance continues the stance it entered on
    dropped_landing: bool = False   # the next unit's landing replaces its own

    def describe(self):
        u = UNITS[self.unit]
        return f'{u["name"]} ({"G=" + self.side})'


def other(side):
    return {'L': 'R', 'R': 'L'}[side]


def resolve(support, side):
    """A phase's support in absolute terms for gesture side `side`."""
    return {'S': other(side), 'G': side}.get(support, support)


def params_for(unit_id, given):
    u = UNITS[unit_id]
    p = {k: v['default'] for k, v in u['params'].items() if 'default' in v}
    p.update(given or {})
    for k, spec in u['params'].items():
        if isinstance(p.get(k), (int, float)) and 'min' in spec:
            p[k] = min(max(p[k], spec['min']), spec['max'])
    return p


def expand_phases(u, p, side):
    """The unit's phases in order with numeric durations, repeats unrolled, supports made absolute; repeated strides
    alternate the feet."""
    out = []
    s = side
    for ph in u['phases']:
        if ph['kind'] == 'repeat':
            n = int(evaluate(ph['count'], p))
            for _ in range(n):
                for sub in ph['of']:
                    q = dict(sub)
                    q['support'] = resolve(sub.get('support', 'none'), s)
                    q['duration'] = evaluate(sub['duration'], p)
                    out.append(q)
                s = other(s)
        else:
            q = dict(ph)
            q['support'] = resolve(ph.get('support', 'none'), side)
            q['duration'] = evaluate(ph['duration'], p) if not isinstance(ph['duration'], (int, float)) else float(ph['duration'])
            if q.get('shapes') == 'along path':
                # The arms pass through the path's positions at equal intervals of the phase.
                path = p['path']
                q['shapes'] = [{'at': k / max(len(path) - 1, 1), 'arms': a} for k, a in enumerate(path)]
            out.append(q)
    return out


def exit_state(u, p, side, entry, t_end, phases):
    ex = u['exit']
    st = replace(entry, t=t_end)
    sup = ex.get('support', 'same')
    if sup == 'same':
        pass
    elif sup == 'alternating':
        last = [ph for ph in phases if ph['kind'] == 'stance']
        st.support = last[-1]['support'] if last else entry.support
        if phases and phases[-1]['kind'] == 'flight' and st.support in ('L', 'R'):
            st.support = other(st.support)       # a run ends in the air and lands on the other foot
    else:
        st.support = resolve(sup, side)
    mom = ex.get('momentum', 'same')
    if mom == 'still':
        st.speed, st.spin = 0.0, 0.0
    elif mom == 'travelling':
        st.speed = p.get('speed_m_s', max(entry.speed, 0.6))
        st.spin = 0.0
    elif mom == 'turning':
        st.spin = p.get('rev_per_s', entry.spin)
        # A turn on a pivot cannot drift: the dancer travels only round the circle the unit names, once in its time.
        st.speed = math.pi * p.get('circle_m', 0.0) / max(t_end - entry.t, 1e-9)
    if 'turn_deg' in p:
        st.heading = (entry.heading + p['turn_deg'] * (1 if side == 'L' else -1)) % 360
    if 'sleeves' in ex:
        st.sleeves = ex['sleeves'] if ex['sleeves'] != 'same' else entry.sleeves
        if st.sleeves == 'aloft':
            vm = math.sqrt(2 * env.FABRIC_KG_M2['silk chiffon'] * G / (env.air_density(0.0) * env.CN_FLAT))
            st.sleeves_down_at = t_end + 1.5 / vm
    # Travel: the speed integrated over the unit along its heading (a screen of floor use, not a path plan).
    dist = 0.5 * (entry.speed + st.speed) * (t_end - entry.t)
    if 'radius_m' in p and p.get('radius_m'):
        dtheta = math.degrees(dist / p['radius_m']) * (1 if side == 'L' else -1)
        st.heading = (entry.heading + dtheta) % 360
    h = math.radians(0.5 * (entry.heading + st.heading))
    st.x, st.y = entry.x + dist * math.cos(h), entry.y + dist * math.sin(h)
    return st


def continues_stance(unit_id, side, st: State):
    """Whether the unit's first phase is a stance on the one foot already bearing the weight from a stance (not a
    landing): then that stance is its push, and a dancer never adds a second stance on the same foot. A weight shift
    ends with both feet still down, so a jump after it keeps its own take-off on the foot the weight moved onto."""
    if st.aloft or st.support not in ('L', 'R') or st.last_kind not in ('stance', 'turn'):
        return False
    first = UNITS[unit_id]['phases'][0]
    if first['kind'] == 'repeat':
        first = first['of'][0]
    return first['kind'] == 'stance' and resolve(first.get('support', 'none'), side) == st.support


def place(unit_id, side, given, entry, connector=False, note='', drop_first=False, drop_landing=False):
    """Place a unit from an entry state. drop_first: its first stance is the one the dancer is already standing in.
    drop_landing: the unit ends in the air because the next unit (a carve) is its landing."""
    u = UNITS[unit_id]
    p = params_for(unit_id, given)
    phases = expand_phases(u, p, side)
    if drop_first and phases and phases[0]['kind'] == 'stance':
        phases = phases[1:]
    if drop_landing and len(phases) >= 2 and phases[-1]['kind'] == 'stance' and phases[-2]['kind'] == 'flight':
        phases = phases[:-1]
    t = entry.t
    for ph in phases:
        ph['t0'] = t
        t += ph['duration']
    ex = exit_state(u, p, side, entry, t, phases)
    if drop_landing:
        ex.support = 'none'
    ex.aloft = bool(phases) and phases[-1]['kind'] == 'flight'
    ex.last_kind = phases[-1]['kind'] if phases else entry.last_kind
    return Placed(unit_id, side, p, entry.t, t - entry.t, phases, entry, ex, connector, note, drop_first, drop_landing)


# --- entry checks and the linking steps that satisfy them ---------------------------------------------------------

def entry_problems(unit_id, side, st: State, given=None):
    u = UNITS[unit_id]
    need = u['entry']
    probs = []
    sup = need.get('support', 'any')
    allowed = sup if isinstance(sup, list) else [sup]
    allowed = [resolve(a, side) for a in allowed]
    if 'any' not in allowed and 'same' not in allowed and st.support not in allowed:
        probs.append(('support', allowed))
    lo, hi = need.get('speed_m_s', [0.0, 99.0])
    if need.get('momentum') == 'still':
        hi = min(hi, 0.15)
    p = params_for(unit_id, given)
    if unit_id in ('grand_jete', 'skim'):          # these carry their approach speed; orbit runs build their own
        lo = max(lo, 0.6 * p.get('speed_m_s', lo))
    if st.speed < lo - 1e-9:
        probs.append(('slow', lo))
    if st.speed > hi + 1e-9:
        probs.append(('fast', hi))
    if st.spin > 0.05 and u['family'] not in ('turn',):
        probs.append(('spinning', 0.0))
    sl = need.get('sleeves')
    if sl == 'gathered' and st.sleeves != 'gathered':
        probs.append(('sleeves', 'gathered'))
    if isinstance(sl, list) and st.sleeves not in sl:
        probs.append(('nothing to gather', sl))        # no linking step makes this unit meaningful here
    return probs


def motion_links(probs, unit_id, side, st: State, given=None):
    """Linking steps for the body's motion and the garment, in a dancer's order: stop the spin, gather the sleeves,
    then build or spend speed."""
    out = []
    s = st
    kinds = {k for k, _ in probs}
    if 'spinning' in kinds:
        # Spend the spin through the feet: a braking step takes the turn's preparation time in reverse.
        pl = place('pas_de_bourree', side if s.support != side else other(side), {'step_s': max(0.3, prep_s(s.spin) / 3)}, s,
                   connector=True, note='brake the turn')
        pl.exit.spin = 0.0
        out.append(pl)
        s = pl.exit
    if 'sleeves' in kinds:
        wait = max(1.0, s.sleeves_down_at - s.t)
        pl = place('sleeve_recover', side, {'seconds': min(3.0, wait)}, s, connector=True,
                   note=f'wait {wait:.1f} s for the sleeves to settle, then gather')
        out.append(pl)
        s = pl.exit
    target = None
    for k, v in probs:
        if k in ('slow', 'fast'):
            target = v
    if target is not None:
        want = params_for(unit_id, given).get('speed_m_s', target)
        if 'fast' in kinds:
            want = target
        n = max(2, math.ceil(abs(want - s.speed) / stride_dv()))
        lead = other(s.support) if s.support in ('L', 'R') else side
        # Run on from a stance: that stance pushes into the first flight, and one more stride does the speed's work.
        cont = continues_stance('orbit_run', lead, s)
        n += 1 if cont else 0
        pl = place('orbit_run', lead, {'speed_m_s': want, 'strides': n, 'radius_m': max(2.0, want ** 2 / (0.5 * MU * G))}, s,
                   connector=True, note=f'{n} strides to {"build" if want > s.speed else "spend"} speed to {want:.2g} m/s',
                   drop_first=cont)
        pl.exit.speed = want
        out.append(pl)
        s = pl.exit
    return out, s


def support_links(unit_id, side, st: State):
    """Linking steps that put the weight where the unit starts: a rise from the floor, a weight change, a close."""
    out = []
    s = st
    need = UNITS[unit_id]['entry'].get('support', 'any')
    allowed = need if isinstance(need, list) else [need]
    allowed = [resolve(a, side) for a in allowed]
    if s.support == 'floor' and 'floor' not in allowed and 'any' not in allowed:
        pl = place('rise_from_floor', side, {}, s, connector=True, note='rise to standing')
        out.append(pl)
        s = pl.exit
    if 'any' not in allowed and s.support not in allowed:
        want = allowed[0]
        pl = None
        if want in ('L', 'R') and s.support == 'both' and s.speed < 0.15:
            # Standing still on both feet: the weight is carried across onto the standing foot in place.
            pl = place('weight_shift', other(want), {}, s, connector=True, note=f'weight onto {want} in place')
        elif want in ('L', 'R') and s.support in ('both', other(want)):
            pl = place('glide_step', other(want), {'speed_m_s': max(0.3, s.speed), 'step_m': 0.3}, s, connector=True,
                       note=f'weight onto {want}')
            pl.exit.support = want
        elif want == 'both':
            # A closing step that also stops the travel: friction brakes it at mu g at most, so the step is long enough
            # to spend the speed (v^2 / (2 mu g)), covered at the mean speed.
            v = s.speed
            braking = v > 0.15
            step = max(0.25, v * v / (2 * MU * G)) if braking else 0.25
            # Out of the air the step lands on the foot the travel comes down on; standing on one foot, the dancer
            # steps onto the other.
            foot = s.support if s.support in ('L', 'R') else other(side)
            pl = place('glide_step', other(foot) if s.aloft else foot,
                       {'speed_m_s': max(0.2, v / 2) if braking else max(0.3, v), 'step_m': step},
                       s, connector=True, note='close to both feet' + (f', stopping {v:.2g} m/s' if braking else ''))
            pl.exit.support = 'both'
            pl.exit.speed = 0.0
        if pl:
            out.append(pl)
            s = pl.exit
    return out, s


@dataclass
class Phrase:
    name: str
    placed: list = field(default_factory=list)

    @property
    def state(self):
        return self.placed[-1].exit if self.placed else State()

    def add(self, unit_id, side='R', flexible_side=True, **given):
        """Place a unit: first the links its motion needs (spin, sleeves, speed), then the side, taken on the other
        foot when only the supporting foot is wrong (as a dancer would), then any weight change still needed."""
        # A landing unit (the carve) replaces the landing of the jump before it, which then ends in the air.
        need = UNITS[unit_id]['entry'].get('support')
        if need in ('none', ['none']) and self.placed and self.state.support != 'none' and not self.placed[-1].connector:
            prev = self.placed[-1]
            if len(prev.phases) >= 2 and prev.phases[-1]['kind'] == 'stance' and prev.phases[-2]['kind'] == 'flight':
                self.placed[-1] = place(prev.unit, prev.side, prev.params, prev.entry, prev.connector, prev.note,
                                        drop_first=prev.dropped_first, drop_landing=True)
        st = self.state
        fits = lambda sd, state: not any(k == 'support' for k, _ in entry_problems(unit_id, sd, state, given))
        probs = [p for p in entry_problems(unit_id, side, st, given) if p[0] != 'support']
        links, st = motion_links(probs, unit_id, side, st, given)
        if flexible_side and not fits(side, st) and fits(other(side), st):
            side = other(side)
        more, st = support_links(unit_id, side, st)
        self.placed.extend(links + more)
        self.placed.append(place(unit_id, side, given, st, drop_first=continues_stance(unit_id, side, st)))
        return self

    @property
    def duration(self):
        return self.placed[-1].exit.t if self.placed else 0.0


# --- devices -------------------------------------------------------------------------------------------------------

def motif(*items):
    """A motif: a list of (unit_id, side, params) to be placed in order."""
    return [(u, s, dict(p)) for u, s, p in items]


def mirror(m):
    return [(u, other(s), p) for u, s, p in m]


def augment(m, k=1.4):
    out = []
    for u, s, p in m:
        q = dict(p)
        for key in ('rise_m', 'revolutions', 'hold_s', 'under_s'):
            if key in q or key in UNITS[u]['params']:
                base = q.get(key, UNITS[u]['params'].get(key, {}).get('default'))
                if isinstance(base, (int, float)):
                    q[key] = base * k
        out.append((u, s, q))
    return out


def diminish(m, k=0.7):
    return augment(m, k)


def accumulate(m):
    out = []
    for i in range(1, len(m) + 1):
        out.extend(m[:i])
    return out


def retrograde(m):
    return list(reversed(m))


def build(name, sequence):
    ph = Phrase(name)
    for u, s, p in sequence:
        ph.add(u, s, **p)
    return ph


# --- evaluation and the score ----------------------------------------------------------------------------------------

def evaluate_phrase(ph: Phrase):
    T = ph.duration
    aloft = sum(p['duration'] for pl in ph.placed for p in pl.phases if p['kind'] == 'flight')
    long_flights = sum(1 for pl in ph.placed for p in pl.phases if p['kind'] == 'flight' and p['duration'] > 1.0)
    native = sum(pl.duration for pl in ph.placed if pl.unit in LUNAR_NATIVE)
    efforts = [UNITS[pl.unit]['effort'] for pl in ph.placed]
    contrast = sum(sum(a[k] != b[k] for k in a) for a, b in zip(efforts, efforts[1:]))
    families = {UNITS[pl.unit]['family'] for pl in ph.placed}
    travel = sum(0.5 * (pl.entry.speed + pl.exit.speed) * pl.duration for pl in ph.placed)
    turning = sum(abs(((pl.exit.heading - pl.entry.heading + 180) % 360) - 180) for pl in ph.placed)
    stops = sum(1 for pl in ph.placed[1:-1] if pl.exit.speed < 0.15 and pl.exit.spin < 0.05
                and UNITS[pl.unit]['family'] not in ('level', 'upper', 'fabric'))
    levels = set()
    for pl in ph.placed:
        for p in pl.phases:
            levels.add({'flight': 'air', 'fall': 'falling'}.get(p['kind'], 'floor' if p['support'] in ('floor', 'kneel') else 'standing'))
    flights = [(pl.t0 + p['t0'] - pl.t0, p['duration']) for pl in ph.placed for p in pl.phases if p['kind'] == 'flight']
    climax = max(flights, key=lambda f: f[1])[0] / T if flights else None
    linking = sum(pl.duration for pl in ph.placed if pl.connector)
    # Revolutions turned in the air, by axis.
    aerial = {}
    for pl in ph.placed:
        for p in pl.phases:
            if p['kind'] != 'flight':
                continue
            for key, axis in (('somersault_rev', 'somersault'), ('cartwheel_rev', 'cartwheel'), ('twist_rev', 'twist'),
                              ('turn_rev', 'turn')):
                if p.get(key) is not None:
                    v = evaluate(p[key], pl.params) if isinstance(p[key], str) else float(p[key])
                    aerial[axis] = round(aerial.get(axis, 0.0) + v, 2)
    return dict(duration_s=round(T, 2), units=len([p for p in ph.placed if not p.connector]),
                connectors=len([p for p in ph.placed if p.connector]), connector_share=round(linking / T, 3),
                families=sorted(families),
                time_aloft_s=round(aloft, 2), share_aloft=round(aloft / T, 3), flights_over_1_s=long_flights,
                lunar_native_share=round(native / T, 3), effort_changes=contrast, levels=sorted(levels),
                travel_m=round(travel, 1), heading_change_deg=round(turning), stops=stops,
                climax_at=round(climax, 2) if climax is not None else None, aerial_revolutions=aerial)


def score(ph: Phrase, bpm=60.0):
    """A count score: the start of each unit in counts at the given tempo, with its support, physics and garment."""
    beat = 60.0 / bpm
    rows = []
    for pl in ph.placed:
        u = UNITS[pl.unit]
        flights = [p['duration'] for p in pl.phases if p['kind'] == 'flight']
        rows.append(dict(
            count=round(pl.t0 / beat + 1, 2), time_s=round(pl.t0, 2), seconds=round(pl.duration, 2),
            unit=u['name'], unit_id=pl.unit, gesture_side=pl.side, connector=pl.connector, note=pl.note,
            params={k: (round(v, 3) if isinstance(v, float) else v) for k, v in pl.params.items()},
            support_in=pl.entry.support, support_out=pl.exit.support,
            speed_in_m_s=round(pl.entry.speed, 2), speed_out_m_s=round(pl.exit.speed, 2),
            flight_s=[round(f, 2) for f in flights], sleeves=pl.exit.sleeves, effort=u['effort'],
            phases=[dict(kind=p['kind'], support=p.get('support'), t0=round(p['t0'], 3),
                         seconds=round(p['duration'], 3), shapes=p.get('shapes', []),
                         turn_rev=p.get('turn_rev'), somersault_rev=p.get('somersault_rev'),
                         cartwheel_rev=p.get('cartwheel_rev'), twist_rev=p.get('twist_rev'), body=p.get('body'))
                    for p in pl.phases]))
    return rows


def score_text(ph: Phrase, bpm=60.0):
    lines = [f'{ph.name}: {ph.duration:.1f} s at {bpm:.0f} beats a minute']
    for r in score(ph, bpm):
        mark = '  (link)' if r['connector'] else ''
        fl = f'  aloft {"+".join(str(f) for f in r["flight_s"])} s' if r['flight_s'] else ''
        lines.append(f'  {r["count"]:6.2f}  {r["time_s"]:6.2f} s  {r["unit"]:<26s} G={r["gesture_side"]} '
                     f'{r["support_in"]:>5s}->{r["support_out"]:<5s} {r["speed_in_m_s"]:.1f}->{r["speed_out_m_s"]:.1f} m/s'
                     f'{fl}{mark}{"  " + r["note"] if r["note"] else ""}')
    return '\n'.join(lines)


# --- the generator -------------------------------------------------------------------------------------------------------

FOLLOWS = {   # which families a choreographer most often lets follow each, as weights (a starting palette)
    'travel': {'flight': 4, 'turn': 2, 'level': 2, 'fabric': 1, 'travel': 1},
    'flight': {'flight': 1, 'travel': 2, 'turn': 2, 'level': 2, 'fabric': 1, 'floor': 1},
    'turn': {'travel': 2, 'level': 2, 'fabric': 2, 'floor': 1, 'flight': 1},
    'level': {'travel': 2, 'flight': 2, 'turn': 2, 'fabric': 1, 'upper': 1},
    'floor': {'level': 2, 'upper': 2, 'travel': 1},
    'upper': {'level': 2, 'travel': 2, 'fabric': 2, 'turn': 1},
    'fabric': {'level': 2, 'travel': 2, 'turn': 1, 'upper': 1},
}


def generate(seed, units=8, palette=None, name=None):
    """A phrase drawn from the family transition weights, placed through the grammar; reproducible from its seed."""
    rng = random.Random(seed)
    palette = palette or [u for u in UNITS if UNITS[u].get('realisable', True)]
    fam_of = {u: UNITS[u]['family'] for u in palette}
    seq = []
    fam = 'upper'
    side = 'R'
    for _ in range(units):
        w = FOLLOWS.get(fam, {})
        fams = [f for f in w if any(fam_of[u] == f for u in palette)]
        fam = rng.choices(fams, weights=[w[f] for f in fams])[0]
        probe = build('probe', seq) if seq else Phrase('probe')
        choices = [u for u in palette if fam_of[u] == fam and (not seq or u != seq[-1][0])
                   and not any(k == 'nothing to gather' for k, _ in entry_problems(u, side, probe.state))]
        if not choices:
            continue
        uid = rng.choice(choices)
        side = side if rng.random() < 0.6 else other(side)
        seq.append((uid, side, {}))
    return build(name or f'generated {seed}', seq), seq


def best_generated(seeds=range(40), units=8):
    """Among generated phrases, the one with the most lunar-native time and effort contrast and the fewest stops."""
    best = None
    for s in seeds:
        ph, seq = generate(s, units)
        e = evaluate_phrase(ph)
        climax = abs((e['climax_at'] if e['climax_at'] is not None else 0.0) - 0.65)
        key = (e['lunar_native_share'] + 0.02 * e['effort_changes'] + 0.05 * e['flights_over_1_s'] - 0.1 * e['stops']
               - 0.5 * e['connector_share'] - 0.3 * climax)
        if best is None or key > best[0]:
            best = (key, s, ph, seq)
    return best


# --- the example phrases ---------------------------------------------------------------------------------------------------

def first_phrase():
    """Under the Canopy: the first lunar phrase, composed by hand. The sleeves go up, the dancer falls under them,
    gathers speed on a curve, hangs a split in a long jeté, steers out of the landing, whirls until the skirt stands
    out, and spirals down as the sleeves come home."""
    return build('Under the Canopy', motif(
        ('port_de_bras', 'R', {'seconds': 3.0, 'path': ['bras_bas', 'first', 'fifth']}),
        ('canopy', 'R', {'under_s': 3.0}),
        ('long_fall', 'R', {'to_deg': 35}),
        ('orbit_run', 'L', {'speed_m_s': 1.6, 'strides': 3, 'radius_m': 4.0}),
        ('grand_jete', 'R', {'rise_m': 1.2, 'speed_m_s': 1.6}),
        ('carve_landing', 'R', {'turn_deg': 45}),
        ('whirl', 'R', {'revolutions': 3, 'rev_per_s': 0.45}),
        ('spiral_descent', 'R', {'revolutions': 1.0}),
    ))


def suspension_study():
    """A motif and its development: a skim into a flight phrase, then the same motif mirrored, augmented, and
    resolved into a long penché."""
    m = motif(('skim', 'R', {'hops': 4}), ('flight_phrase', 'R', {'rise_m': 1.2, 'speed_m_s': 0.8}))
    seq = m + mirror(m) + augment(m, 1.4) + motif(('long_penche', 'R', {'hold_s': 2.5}))
    return build('Suspension Study (motif, mirror, augmentation)', seq)


def air_study():
    """Air Study: the aerial figures the Moon's long flights open, composed by hand. A run into a star wheel, a
    hanging salto and a corkscrew on the spot, a sleeve flight that leaves the sleeves standing, the highest jump laid
    out flat and spinning, and a spiral down to a kneel as the sleeves come home."""
    return build('Air Study', motif(
        ('orbit_run', 'R', {'speed_m_s': 0.8, 'strides': 2, 'radius_m': 3.0}),
        ('star_wheel', 'R', {'rise_m': 1.4, 'speed_m_s': 0.8}),
        ('hanging_salto', 'R', {'rise_m': 1.6}),
        ('corkscrew', 'R', {'rise_m': 1.6, 'twists': 1}),
        ('sleeve_flight', 'R', {'rise_m': 1.4, 'speed_m_s': 0.8}),
        ('flat_spin', 'R', {'rise_m': 2.3}),
        ('spiral_descent', 'R', {'revolutions': 1.0}),
    ))


def canopy_and_salto():
    """Under the Canopy with the hanging salto: the first phrase up to its whirl, then the aerial figure that the
    motion engine realised first, as its ending. The sleeves go up and the dancer falls under them, gathers speed,
    hangs a split in the long jete and steers out of it, runs the speed down into the whirl, closes to both feet and
    turns over backward in a salto that hangs upside down, landing still."""
    return build('Under the Canopy, with the hanging salto', motif(
        ('port_de_bras', 'R', {'seconds': 3.0, 'path': ['bras_bas', 'first', 'fifth']}),
        ('canopy', 'R', {'under_s': 3.0}),
        ('long_fall', 'R', {'to_deg': 35}),
        ('orbit_run', 'L', {'speed_m_s': 1.6, 'strides': 3, 'radius_m': 4.0}),
        ('grand_jete', 'R', {'rise_m': 1.2, 'speed_m_s': 1.6}),
        ('carve_landing', 'R', {'turn_deg': 45}),
        ('whirl', 'R', {'revolutions': 3, 'rev_per_s': 0.45}),
        ('hanging_salto', 'R', {'rise_m': 1.6}),
    ))


def results():
    phrases = [first_phrase(), suspension_study(), air_study(), canopy_and_salto()]
    key, seed, gen, _ = best_generated()
    gen.name = f'Generated, seed {seed} (best of 40 by lunar share, contrast, linking and climax)'
    phrases.append(gen)
    return dict(
        schema=SCHEMA,
        producer=dict(domain='humanities', files={f: env.digest(f) for f in PRODUCER_FILES},
                      inputs={'humanities/dance/results/envelope.json': env.digest('humanities/dance/results/envelope.json')}),
        evidence='Phrases placed by the grammar from the vocabulary, with durations and entry checks from the '
                 'envelope\'s physics: screens of feasibility and structure. The motion engine realises them with '
                 'whole-body dynamics; artistic worth is for the author and dancers to judge.',
        vocabulary=dict(units=len(UNITS), lunar_native=sorted(LUNAR_NATIVE),
                        families={f: sorted(u for u in UNITS if UNITS[u]['family'] == f)
                                  for f in sorted({u['family'] for u in UNITS.values()})}),
        phrases=[dict(name=p.name, evaluation=evaluate_phrase(p), score=score(p), text=score_text(p)) for p in phrases])


def main() -> int:
    OUT.parent.mkdir(exist_ok=True)
    r = results()
    OUT.write_text(json.dumps(r, indent=1) + '\n')
    for p in r['phrases']:
        print(p['text'])
        print('  ', p['evaluation'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
