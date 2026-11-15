"""What the dance's solved whole-body motions did when they landed, took off and turned through a stance.

    python -m humanities.sport.whole_body     # writes humanities/sport/results/whole_body.json (needs the drive)

The envelope (envelope.py) sets bounds by closed forms; this record measures solved motions of the motion engine, the
55 kg body model on the dance's vinyl floor (friction 0.5) at lunar gravity, and on Earth for a check, each from its
motion file and per-node report on the research drive:
- the whole phrase *Under the Canopy, with the hanging salto* as solved on 6 November 2026, before landings had a
  stroke;
- the hanging salto as a short phrase of its own and the run, grand jete and carve landing of the longer phrase, solved
  on 7 November with the engine's landing impulses costed on what they do to the centre of mass and held inside the
  friction cone: first with the engine's Euler flights, then the salto with Runge-Kutta flights in two stages, landing
  flat-footed and then on the balls of the feet, rolling down into a deep bend and rising out of it;
- a check landing solved on 8 November on Earth and on the Moon, a sauté that touches down on the balls of the feet,
  sinks on them into a deep bend and rises, without and with the engine's cost on landing forces (weight 100), and the
  hanging salto landing the same way with the cost.

For every stance that a flight ends and another flight starts (strides, the grand jete's take-off and its carve
landing), every landing from a flight, the hanging salto's standing take-off and each motion's last state, it reads the
centre of mass's velocity before and after, the floor's forces through the stance, the change of velocity the forces do
not account for, how far the centre of mass sinks after touchdown and where it ends over the feet. The engine stops each
landing foot with an instant impulse (realise.py), which the motion file's forces leave out, so that unaccounted change
is the landing impulse; its sideways share against its vertical share shows whether it stayed inside the floor's grip.
Velocities come from fits to the free-flight frames on either side of a stance, or central differences of the 60 frames
a second the motion file carries inside one; the flights' own free fall checks them.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.humanities.sport-whole-body/1'
OUT = HERE / 'results' / 'whole_body.json'
FILES = ('humanities/sport/whole_body.py',)
DRIVE = Path('/media/projectspace/terluna-research/humanities/dance-runs')
LANDING = DRIVE / 'landing_2026-11-07' / 'motion'
LOADS = DRIVE / 'landing_force_2026-11-08' / 'motion'
MOTIONS = {
    # the whole phrase before landings were costed: rigid landings
    'phrase_2026_11_06': DRIVE / 'motion_apart' / 'canopy_and_salto',
    # the landing impulse costed (weight 1e3) with Euler flights: the optimiser floats down on the integration error
    'salto_landing_cost_euler_2026_11_07': LANDING / 'hanging_salto_stroke',
    'leap_landing_cost_euler_2026_11_07': LANDING / 'leap_stroke',
    # Runge-Kutta flights, solved in two stages (the landing cost, weight 1e2, added to a solution without it): the
    # salto landing flat-footed with straight legs, then touching down on the balls of the feet and rolling down
    'salto_flat_landing_rk4_2026_11_07': LANDING / 'salto_flatfoot_staged_rk4_flight_w100_stage2',
    'salto_ball_landing_rk4_2026_11_07': LANDING / 'salto_staged_rk4_flight_w100_stage2',
    # The check landing on the balls of the feet into a deep bend, on Earth and on the Moon, without the cost on landing
    # forces and with it (weight 100), and the salto landing the same way with the cost (8 November)
    'check_earth_2026_11_08': LOADS / 'check2_earth_stage2',
    'check_earth_landing_cost_2026_11_08': LOADS / 'check2_earth_stage3_w100',
    'check_moon_2026_11_08': LOADS / 'check2_moon_stage2',
    'check_moon_landing_cost_2026_11_08': LOADS / 'check2_moon_stage3_w100',
    'salto_on_the_balls_landing_cost_2026_11_08': LOADS / 'salto_balls_stage3_w100',
}
MOTION = MOTIONS['phrase_2026_11_06'].with_suffix('.npz')        # the reference motion the tests need
MASS_KG = 55.0
FLOOR_FRICTION = 0.5


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def r(x, n=3):
    return round(float(x), n)


def paths(stem):
    return stem.with_suffix('.npz'), stem.parent / (stem.name + '_report.json')


def load(stem):
    motion, report = paths(stem)
    with np.load(motion, allow_pickle=False) as z:
        m = {k: z[k] for k in ('time', 'com', 'contact_on', 'contact_force', 'seg_pos', 'segments')}
        m['meta'] = json.loads(str(z['meta']))
    return m, json.loads(report.read_text())


def spans(on):
    """(start, end) frame ranges, end exclusive, of consecutive frames where `on` holds."""
    out, i, n = [], 0, len(on)
    while i < n:
        if on[i]:
            j = i
            while j < n and on[j]:
                j += 1
            out.append((i, j))
            i = j
        else:
            i += 1
    return out


def velocity(com, t):
    return np.gradient(com, t, axis=0)


def flight_velocity(m, frames, at):
    """The centre of mass's velocity at time `at` from a fit to free-flight frames alone (straight lines along, a
    parabola up), so that a landing's impulse in the next frame does not enter it."""
    tt, c = m['time'][frames], m['com'][frames]
    return np.array([np.polyval(np.polyder(np.polyfit(tt, c[:, k], 2 if k == 2 else 1)), at) for k in range(3)])


def velocity_before(m, a, k=6):
    """Velocity at the last airborne frame before a stance starting at frame a."""
    on = m['contact_on'].any(axis=1)
    frames = [i for i in range(max(a - k, 0), a) if not on[i]]
    return flight_velocity(m, frames, m['time'][a - 1])


def velocity_after(m, b, k=6):
    """Velocity at the first airborne frame after a stance ending before frame b."""
    on = m['contact_on'].any(axis=1)
    frames = [i for i in range(b, min(b + k, len(on))) if not on[i]]
    return flight_velocity(m, frames, m['time'][b])


def flight_checks(m, v):
    """Each flight's measured vertical acceleration (a straight-line fit of its vertical velocity) against -g."""
    g = m['meta']['gravity']
    out = []
    for a, b in spans(~m['contact_on'].any(axis=1)):
        if b - a < 12:
            continue
        sel = slice(a + 2, b - 2)
        slope = np.polyfit(m['time'][sel], v[sel, 2], 1)[0]
        out.append(dict(t=[r(m['time'][a], 3), r(m['time'][b - 1], 3)], vertical_acceleration_m_s2=r(slope, 3),
                        share_of_gravity=r(-slope / g, 3)))
    return out


def stance_record(m, v, a, b, unit, g, ends_in_flight=True):
    """One stance, frames a..b-1 bearing load: velocities either side, the floor's impulse, the impulse the forces
    leave out, the friction bound on the whole change, the peak push and the turn of the travel. The velocity before
    comes from the flight that ends at the stance; the velocity after, from the next flight, or inside the stance
    (central differences) for a window that ends on the floor."""
    t, F = m['time'], m['contact_force']
    dt = float(t[1] - t[0])
    i0, i1 = max(a - 1, 0), min(b, len(t) - 1)          # the last frame before the stance, the first after it
    v_in = velocity_before(m, a)
    v_out = velocity_after(m, b) if ends_in_flight else v[i1]
    N = F[a:b, :, 2].sum(axis=1)
    H = F[a:b, :, :2].sum(axis=1)
    forced = np.concatenate([H.sum(axis=0), [N.sum()]]) * dt / MASS_KG - np.array([0.0, 0.0, g * (t[i1] - t[i0])])
    total = v_out - v_in
    jump = total - forced                                # the change the stored forces leave out: the landing impulses
    vert_impulse = N.sum() * dt / MASS_KG + max(jump[2], 0.0)
    dvh = float(np.linalg.norm(total[:2]))
    hin, hout = math.atan2(v_in[1], v_in[0]), math.atan2(v_out[1], v_out[0])
    turn = math.degrees((hout - hin + math.pi) % (2 * math.pi) - math.pi) if min(np.linalg.norm(v_in[:2]), np.linalg.norm(v_out[:2])) > 0.2 else None
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(N > 5.0, np.linalg.norm(H, axis=1) / N, 0.0)
    return dict(
        unit=unit, t=[r(t[a]), r(t[b - 1])], stance_s=r(t[i1] - t[i0]),
        before=dict(along_m_s=r(np.linalg.norm(v_in[:2]), 2), vertical_m_s=r(v_in[2], 2)),
        after=dict(along_m_s=r(np.linalg.norm(v_out[:2]), 2), vertical_m_s=r(v_out[2], 2)),
        travel_change_m_s=r(dvh, 2), travel_turn_deg=None if turn is None else r(turn, 0),
        floor_vertical_impulse_per_kg_m_s=r(N.sum() * dt / MASS_KG, 2),
        landing_impulse_m_s=dict(along=r(np.linalg.norm(jump[:2]), 2), vertical=r(jump[2], 2),
                                 along_over_vertical=r(np.linalg.norm(jump[:2]) / jump[2], 2) if jump[2] > 0.3 else None),
        friction_bound_on_travel_change_m_s=r(FLOOR_FRICTION * vert_impulse, 2),
        share_of_bound_used=r(dvh / (FLOOR_FRICTION * vert_impulse), 2) if vert_impulse > 0 else None,
        peak_push=dict(N=r(N.max(), 0), lunar_bw=r(N.max() / (MASS_KG * MOON_SURFACE_GRAVITY), 1),
                       earth_bw=r(N.max() / (MASS_KG * STANDARD_GRAVITY), 2)),
        most_friction_used=r(ratio.max(), 2))


def salto_take_off_and_landing(m, v, rep, stances, on, g):
    """The hanging salto's standing take-off, from its own stance in the report (in the longer phrase the whirl, the
    pas de bourree and the glide step stand on the floor before it without a flight) from the start of the dip to
    lift-off, and the first 0.4 s of its landing."""
    t, com = m['time'], m['com']
    t0 = rep['nodes'][0]['t']                           # the report's time of the motion file's zero
    salto = next(a for a in rep['meta']['aerial'] if a['unit'].startswith('hanging_salto'))
    t_off, t_land = (x - t0 for x in salto['t'])
    st = [n['t'] - t0 for n in rep['nodes'] if n['unit'].startswith('hanging_salto') and n['phase'] == 'stance' and n['t'] - t0 < t_off]
    a = int(np.searchsorted(t, st[0]))
    b = int(np.searchsorted(t, t_off - 1e-9))
    while b < len(t) and on[b]:
        b += 1                                           # the first airborne frame
    dip_i = a + int(np.argmin(com[a:b, 2]))
    onset = dip_i
    while onset > a and v[onset - 1, 2] < -0.02:
        onset -= 1                                       # back to where the descent began
    standing = dict(unit='hanging_salto, standing take-off', stance_from_s=r(st[0]), descent_starts_s=r(t[onset]),
                    dip_m=r(com[onset, 2] - com[dip_i, 2]), dip_lowest_at_s=r(t[dip_i]),
                    descent_s=r(t[dip_i] - t[onset]), dip_to_lift_off_s=r(t[b] - t[dip_i]),
                    onset_to_lift_off_s=r(t[b] - t[onset]), lift_off_vertical_m_s=r(velocity_after(m, b)[2], 2),
                    peak_push=dict(N=r(m['contact_force'][a:b, :, 2].sum(axis=1).max(), 0),
                                   lunar_bw=r(m['contact_force'][a:b, :, 2].sum(axis=1).max() / (MASS_KG * MOON_SURFACE_GRAVITY), 1)))
    kl = min(k for k, (a, b) in enumerate(stances) if t[a] > t_land - 0.05)
    a, b = stances[kl]
    landing = stance_record(m, v, a, min(b, a + int(round(0.4 / float(t[1] - t[0])))),
                            'hanging_salto, landing (first 0.4 s)', g, ends_in_flight=False)
    return standing, landing


def measure(stem):
    m, rep = load(stem)
    g = m['meta']['gravity']
    assert min(abs(g - MOON_SURFACE_GRAVITY), abs(g - STANDARD_GRAVITY)) < 1e-9     # the Moon, or Earth for a check
    t, com = m['time'], m['com']
    v = velocity(com, t)
    on = m['contact_on'].any(axis=1)
    stances = spans(on)
    # The motion file's time starts at zero at the first node solved; the report keeps the phrase's own times.
    t0 = rep['nodes'][0]['t']
    unit_at = {round(n['t'] - t0, 3): n['unit'] for n in rep['nodes']}

    def unit_of(time):
        keys = np.array(list(unit_at))
        return unit_at[float(keys[np.argmin(np.abs(keys - time))])]

    between = []
    for k, (a, b) in enumerate(stances):
        flight_before = a > 0 and not on[a - 1]
        flight_after = b < len(t) and not on[b]
        if flight_before and flight_after and t[b] - t[a - 1] <= 1.0:      # strides, take-offs and landings
            between.append(stance_record(m, v, a, b, f'{unit_of(t[a - 1])} into {unit_of(t[b])}', g))
    standing = landing = None
    if any(n['unit'].startswith('hanging_salto') and n['phase'] == 'flight' for n in rep['nodes']):
        standing, landing = salto_take_off_and_landing(m, v, rep, stances, on, g)
    landings = []
    for a, b in spans(~on):
        if a == 0 or b >= len(t) or b - a < 12:
            continue                                     # a flight with a landing after it
        e = b
        while e < len(t) and on[e] and t[e] - t[b] < 2.0:
            e += 1
        rec = stance_record(m, v, b, e, f'{unit_of(t[a])} landing', g, ends_in_flight=e < len(t) and not on[e])
        low = b + int(np.argmin(com[b:e, 2]))
        rec['stroke'] = dict(com_at_touchdown_m=r(com[b, 2]), lowest_m=r(com[low, 2]), sinks_m=r(com[b, 2] - com[low, 2]),
                             lowest_after_s=r(t[low] - t[b]))
        landings.append(rec)
    return dict(motion=dict(file=str(paths(stem)[0]), phrase=rep['meta'].get('phrase'), span_s=rep['meta']['span_s'],
                            iterations=rep['meta'].get('iterations'), dynamic_gap=rep['meta'].get('dynamic_gap')),
                flights=flight_checks(m, v), stances_between_flights=between, landings=landings,
                salto_standing_take_off=standing, salto_landing=landing, finish=finish(m, v))


def finish(m, v):
    """Where the body ends at the last frame: how far its centre of mass stands, in the horizontal, from the midpoint
    of the two ankles (the feet's segment frames; the motion file carries no contact positions), how hard the floor
    still pushes and how fast the body still moves."""
    F = m['contact_force'][-1]
    names = [str(n) for n in m['segments']]
    w = F[:, 2]
    out = dict(com_speed_m_s=r(np.linalg.norm(v[-1]), 3))
    if w.sum() > 1.0:
        feet = [names.index(f'{s}_foot') for s in 'LR']
        mid = m['seg_pos'][-1][feet].mean(axis=0)
        out.update(com_over_feet_m=r(np.linalg.norm(m['com'][-1, :2] - mid[:2]), 3), floor_push_N=r(w.sum(), 1))
    return out


def results():
    present = {k: stem for k, stem in MOTIONS.items() if paths(stem)[0].is_file()}
    inputs = {}
    for stem in present.values():
        for path in paths(stem):
            inputs[str(path)] = sha(path)
    return dict(
        schema=SCHEMA,
        producer=dict(domain='humanities', files={f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest()[:16] for f in FILES},
                      constants=constants_used([ROOT / f for f in FILES]), inputs=inputs),
        evidence='Measurements of solved whole-body motions (multiple-shooting optimal control with rigid contacts, '
                 'friction cones on the stance forces, torque bounds and landing impulses, steps of 25 ms by Euler or, '
                 'in flight, Runge-Kutta; since 8 November a cost on landing forces) for one 55 kg body on a 0.5 '
                 'friction floor at lunar gravity, and on Earth for a check; feasible motions of that model, not '
                 'measured human movement.',
        body=dict(mass_kg=MASS_KG, floor_friction=FLOOR_FRICTION),
        motions={k: measure(stem) for k, stem in present.items()})


def main() -> int:
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
