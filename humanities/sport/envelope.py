"""What an athlete's body can do at lunar gravity in the Open Moon's air: the moves a sport can be made of.

    python -m humanities.sport.envelope     # writes humanities/sport/results/envelope.json at once

The sport study starts where the dance did, from the body before any rules. Each figure sets the Moon beside Earth for
the same body. The families follow what play asks of a body: getting going, stopping and turning on the ground;
jumping and the commitment a jump makes; the flight; landing; throwing, catching, kicking and striking; contact
between players; getting down and getting up. Balls fly in the design run's air at the sea-level datum.

Bodies. The dance's body model (humanities/dance/body/model.py: 1.64 m, 55 kg, de Leva's female proportions) and a
larger athlete of 1.80 m and 80 kg, which is that body scaled: lengths by stature, masses by mass, inertias by mass
times length squared, the push-off by stature and every reach by stature, with the same force-velocity profile per
kilogram. The larger body stands in for a male body model, which the dance work leaves for later (its question 2).

Assumptions, each where its figures are made:
- Grip. An athletic shoe on a sports floor takes a friction coefficient of 0.8, an assumed value (the earlier sports
  screens used 0.8; no measurement was retrieved this round); the dance's vinyl floor (0.5) and a high-grip 1.2 sit
  beside it where grip decides the figure.
- Balance. Standing, the floor's push can move 0.13 m to the side of the centre of mass (half the hip width and half a
  foot) and 0.10 m ahead of or behind it; a long split stance puts a rear foot 0.50 m behind it; a wide ready stance
  sets the feet 0.70 m apart with the centre of mass at 0.75 m. Placements for the dance body, scaled for the larger.
- Strides and landings. A running stride's stance lasts 0.33 s, as in the solved phrase's strides; a landing's legs
  absorb it over a 0.4 m stroke; the legs brake a drop into a crouch at up to 15 m/s^2 on either world (assumed,
  muscle-limited).
- Throws. The hand pushes the ball for 0.15 s and releases an overhand throw 1.0 m above the standing centre of mass
  (assumed). A kick swings the leg at 15 rad/s at impact, built over 0.15 s (assumed).
- Balls. Masses and diameters are the governing bodies' sizes, from general knowledge, not retrieved this round; drag is
  a smooth sphere's below the drag crisis (C_D 0.47); the shuttlecock's drag area comes from its Earth terminal speed of
  about 6.7 m/s. A soft ball of 0.10 kg and 0.20 m is an illustrative choice.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from habitation.low_gravity.movement import bank_deg
from humanities.dance import envelope as dance
from humanities.dance.body import pose
from humanities.dance.body.model import MASS_KG, STATURE_M, spec as body_spec
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.humanities.sport-envelope/1'
OUT = HERE / 'results' / 'envelope.json'
FILES = ('humanities/sport/envelope.py', 'humanities/dance/envelope.py', 'humanities/dance/body/model.py',
         'humanities/dance/body/pose.py', 'habitation/low_gravity/movement.py')
INPUTS = (dance.WINDS, dance.POSES)

G_M = MOON_SURFACE_GRAVITY
G_E = STANDARD_GRAVITY
WORLDS = (('moon', G_M), ('earth', G_E))
RHO_EARTH = dance.RHO_EARTH

MU = 0.8                                  # athletic shoe on a sports floor (assumed)
MU_RANGE = (0.5, 0.8, 1.2)                # the dance floor, the central value, a high-grip floor
LARGE = dict(stature_m=1.80, mass_kg=80.0)
SIDE_REACH_M = 0.13                       # the floor's push beside the standing centre of mass
FORE_AFT_REACH_M = 0.10                   # ahead of or behind it, feet under the body
SPLIT_REACH_M = 0.50                      # to a rear foot in a long split stance (assumed)
READY_HALF_STANCE_M = 0.35                # a wide ready stance, feet 0.70 m apart (assumed)
READY_COM_M = 0.75                        # the centre of mass in that crouch (assumed)
STRIDE_S = 0.33                           # a running stride's stance
STROKE_M = 0.4                            # a landing's stroke
BRAKE_M_S2 = 15.0                         # the legs' braking of a drop into a crouch (assumed)
SHAPE_S = 0.75                            # one held shape in a flight (the dance envelope's 0.7-0.8 s)
APEX_TAKE_OFF_M_S = 2.3                   # the jump a contact or a release in the air happens at the top of
THROW_S = 0.15                            # the hand's push on a thrown or caught ball (assumed)
RELEASE_ABOVE_COM_M = 1.0                 # an overhand release above the standing centre of mass (assumed)
AIR_LEVERS_M = {'chest pass, near the centre of mass': 0.3, "overhead, at arm's length": 0.8}
KICK_RAD_S = 15.0                         # the kicking leg's swing at impact (assumed)
KICK_BUILD_S = 0.15
CD_SPHERE = 0.47
SHUTTLE_EARTH_TERMINAL_M_S = 6.7
BALLS = {                                 # mass kg, diameter m
    'basketball (size 7)': (0.62, 0.239),
    'handball (size 3)': (0.45, 0.19),
    'football (size 5)': (0.43, 0.22),
    'volleyball': (0.27, 0.21),
    'soft ball (illustrative)': (0.10, 0.20),
    'shuttlecock': (0.0052, None),
}
THROW_BALL = 'handball (size 3)'
RELEASE_M_S = (5.0, 10.0, 15.0)
RELEASE_HEIGHT_M = 2.0                    # a ball's release above the floor for its flight
TOSS_FORCE_N = dance.THROW_FORCE_N        # a two-arm push under a teammate, over the arms' stroke (the dance's partner)
TOSS_STROKE_M = dance.THROW_STROKE_M


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16]


def r(x, n=3):
    return round(float(x), n)


# --- the bodies ---------------------------------------------------------------------------------------------------

def leg_inertia_about_hip(side='R'):
    """One leg's inertia about its hip joint for a swing in the sagittal plane (axis y), from the body model's thigh,
    shank, foot and toes in the rest pose."""
    sp = body_spec()
    hip = np.array(next(s for s in sp['segments'] if s['name'] == f'{side}_thigh')['origin'], float)
    total = mass = 0.0
    for s in sp['segments']:
        if s['name'] in (f'{side}_thigh', f'{side}_shank', f'{side}_foot', f'{side}_toes'):
            c = np.array(s['origin'], float) + np.array(s['com'], float)
            d = c - hip
            total += s['inertia_kg_m2'][1] + s['mass_kg'] * (d[0] ** 2 + d[2] ** 2)
            mass += s['mass_kg']
    return mass, total


def bodies():
    """The two bodies' masses, centres of mass, inertias in the flight shapes and reaches."""
    M, com, _ = dance.whole_body_inertia()
    I = {k: pose.inertia(s)[2] for k, s in dance.AERIAL_SHAPES.items()}
    leg_m, leg_I = leg_inertia_about_hip()
    base = dict(stature_m=STATURE_M, mass_kg=r(M, 2), com_height_m=r(com[2]),
                inertia_kg_m2=dict(layout=r(I['layout, arms overhead'][1, 1]), tuck=r(I['tuck'][1, 1]),
                                   twist=r(I['layout, arms gathered'][2, 2])),
                push_off_m=dance.H_PO_M, side_reach_m=SIDE_REACH_M, fore_aft_reach_m=FORE_AFT_REACH_M,
                split_reach_m=SPLIT_REACH_M, ready_half_stance_m=READY_HALF_STANCE_M, ready_com_m=READY_COM_M,
                release_above_com_m=RELEASE_ABOVE_COM_M, leg=dict(mass_kg=r(leg_m, 2), inertia_about_hip_kg_m2=r(leg_I)))
    sl, sm = LARGE['stature_m'] / STATURE_M, LARGE['mass_kg'] / MASS_KG
    large = dict(stature_m=LARGE['stature_m'], mass_kg=LARGE['mass_kg'], com_height_m=r(com[2] * sl),
                 inertia_kg_m2={k: r(v * sm * sl ** 2) for k, v in base['inertia_kg_m2'].items()},
                 push_off_m=r(dance.H_PO_M * sl), side_reach_m=r(SIDE_REACH_M * sl),
                 fore_aft_reach_m=r(FORE_AFT_REACH_M * sl), split_reach_m=r(SPLIT_REACH_M * sl),
                 ready_half_stance_m=r(READY_HALF_STANCE_M * sl), ready_com_m=r(READY_COM_M * sl),
                 release_above_com_m=r(RELEASE_ABOVE_COM_M * sl),
                 leg=dict(mass_kg=r(leg_m * sm, 2), inertia_about_hip_kg_m2=r(leg_I * sm * sl ** 2)))
    return {'dancer, 1.64 m and 55 kg': base, 'larger athlete, 1.80 m and 80 kg (scaled)': large}


def take_off(body, g):
    return dance.samozino_take_off(g, h=body['push_off_m'])


# --- on the ground: getting going, stopping, turning --------------------------------------------------------------

def lean_start(distance, g, h, reach):
    """A standing body leaning off its support before any step: an inverted pendulum x'' = (g/h)(x - p) with the
    floor's push held at the support's edge p = -reach, so x = reach (cosh(t/tau) - 1) and v = (reach/tau) sinh(t/tau).
    Returns the time to carry the centre of mass `distance` sideways and the speed it then has."""
    tau = math.sqrt(h / g)
    t = tau * math.acosh(1 + distance / reach)
    return dict(seconds=r(t), speed_m_s=r(reach / tau * math.sinh(t / tau), 2))


def cut_bounce(speed, angle_deg, mu, g, stance_s=STRIDE_S):
    """The vertical speed a cut through one stance needs: turning the travel by `angle_deg` at constant speed changes
    the velocity by 2 v sin(angle/2), and the floor's grip gives at most mu times the vertical impulse, which carries
    the landing's vertical speed w, the take-off's w and the weight over the stance: mu (2 w + g T) >= 2 v sin(a/2)."""
    dv = 2 * speed * math.sin(math.radians(angle_deg) / 2)
    w = max(0.0, (dv / mu - g * stance_s) / 2)
    return dict(velocity_change_m_s=r(dv), vertical_speed_each_way_m_s=r(w), bounce_rise_m=r(w * w / (2 * g)),
                bounce_airtime_s=r(2 * w / g))


def braking_landing(speed, v_down, mu, g, stroke=STROKE_M):
    """A landing that brakes travel: the legs absorb v_down over the stroke in T = 2 d / v_down, the floor's vertical
    impulse per kilogram is v_down + g T, and its grip takes at most mu times that off the travel; the rest slides
    off at mu g."""
    T = 2 * stroke / v_down
    shed = min(speed, mu * (v_down + g * T))
    left = speed - shed
    return dict(stance_s=r(T), shed_in_landing_m_s=r(shed), left_m_s=r(left),
                distance_in_landing_m=r((speed + left) / 2 * T), slide_after_m=r(left * left / (2 * mu * g)),
                total_m=r((speed + left) / 2 * T + left * left / (2 * mu * g)),
                mean_push_body_weights=dict(lunar=r((v_down * v_down / (2 * stroke) + g) / G_M, 2),
                                            earth=r((v_down * v_down / (2 * stroke) + g) / G_E, 2)))


def ground(body):
    h = body['com_height_m']
    out = dict(friction=MU)
    out['balance_clock_s'] = {w: r(math.sqrt(h / g)) for w, g in WORLDS}
    out['lean_start'] = [dict(distance_m=d, **{w: lean_start(d, g, h, body['side_reach_m']) for w, g in WORLDS})
                         for d in (0.1, 0.3)]
    out['friction_limited_start'] = [dict(speed_m_s=v, **{w: dict(seconds=r(v / (MU * g), 2), metres=r(v * v / (2 * MU * g), 2))
                                                             for w, g in WORLDS}) for v in (2.0, 3.0, 4.0)]
    out['stride_change_m_s'] = {f'vertical_{vz}_m_s_each_way': {w: r(MU * (2 * vz + g * STRIDE_S), 2) for w, g in WORLDS}
                                for vz in (0.4, 0.6)}
    out['stops'] = []
    for v in (2.0, 3.0, 4.0):
        row = dict(speed_m_s=v, sliding_m={w: r(v * v / (2 * MU * g), 2) for w, g in WORLDS})
        for vd in (1.5, 2.5):
            row[f'braking_landing_from_{vd}_m_s_down'] = {w: braking_landing(v, vd, MU, g) for w, g in WORLDS}
        out['stops'].append(row)
    out['steady_turn_radius_m'] = [dict(speed_m_s=v, **{f'mu_{mu}': {w: r(v * v / (mu * g), 2) for w, g in WORLDS}
                                                         for mu in MU_RANGE}) for v in (2.0, 3.0, 4.0)]
    out['bank_for_a_5_m_circle_deg'] = {f'{v}_m_s': r(bank_deg(v, 5.0, G_M, MU), 1) for v in (2.0, 3.0, 4.0)}
    out['cuts_through_one_stance'] = [dict(speed_m_s=v, angle_deg=a, moon=cut_bounce(v, a, MU, G_M),
                                           earth=cut_bounce(v, a, MU, G_E))
                                      for v in (2.0, 3.0, 4.0) for a in (30, 60, 90, 180)]
    return out


# --- jumps and commitment ------------------------------------------------------------------------------------------

def dip_s(depth, g, brake=BRAKE_M_S2):
    """Shortest time to lower the centre of mass by `depth` from rest and stop: it can fall no faster than g, then the
    legs brake at up to `brake`, so it falls d1 = depth brake / (g + brake) and brakes over the rest."""
    d1 = depth * brake / (g + brake)
    return math.sqrt(2 * d1 / g) + math.sqrt(2 * (depth - d1) / brake)


def rise_and_stop_s(height, a_up, g):
    """Shortest time to raise the centre of mass by `height` and arrive at rest without leaving the floor: the legs
    push up at a_up, but slowing the rise is gravity's alone (the floor cannot pull), so the body coasts up the last
    d2 = height a_up / (a_up + g)."""
    d2 = height * a_up / (a_up + g)
    return math.sqrt(2 * (height - d2) / a_up) + math.sqrt(2 * d2 / g)


def jumps(body):
    out = {}
    for w, g in WORLDS:
        v = take_off(body, g)
        out[w] = dict(take_off_m_s=r(v), rise_m=r(v * v / (2 * g)), airtime_s=r(2 * v / g),
                      push_s=r(2 * body['push_off_m'] / v))
    depth = 0.29                                       # the solved salto's dip
    out['countermovement'] = dict(
        dip_m=depth, dip_s={w: r(dip_s(depth, g)) for w, g in WORLDS},
        dip_and_push_s={w: r(dip_s(depth, g) + out[w]['push_s']) for w, g in WORLDS},
        note='a standing jump announces itself by its dip; a jump out of a landing (a rebound) needs none')
    ratio = min(MU, body['ready_half_stance_m'] / body['ready_com_m'])
    out['lateral_hop'] = dict(push_ratio=r(ratio), limited_by='stance' if ratio < MU else 'grip', hops=[
        dict(vertical_m_s=vz, sideways_m_s=r(ratio * vz),
             **{w: dict(rise_m=r(vz * vz / (2 * g)), airtime_s=r(2 * vz / g), covered_m=r(ratio * vz * 2 * vz / g))
                for w, g in WORLDS}) for vz in (1.0, 1.5, 2.3)])
    out['commitment'] = [dict(rise_m=h, take_off_m_s=r(math.sqrt(2 * G_M * h)), airtime_s=r(2 * math.sqrt(2 * h / G_M)),
                              shapes_in_flight=r(2 * math.sqrt(2 * h / G_M) / SHAPE_S, 1),
                              earth_airtime_same_take_off_s=r(2 * math.sqrt(2 * G_M * h) / G_E))
                         for h in (0.3, 0.6, 1.0, 1.6, round(out['moon']['rise_m'], 2))]
    return out


def take_off_turning(body, v, g):
    """The dance envelope's estimate of the somersault a take-off can give, for this body: the mean push over the
    push-off times the lean past the foot (0.25 m for the dance body, scaled by stature). It bounds what the floor
    can give; the solved hanging salto took 14.7 N m s of its push (salto variants, 2 November)."""
    s = body['stature_m'] / STATURE_M
    t = 2 * body['push_off_m'] / v
    push = body['mass_kg'] * (v / t + g)
    return dict(push_s=r(t), mean_push_N=r(push, 1), somersault_N_m_s=r(push * dance.AERIAL_LEAN_M * s * t, 1))


def flight(body):
    out = {}
    for w, g in WORLDS:
        v = take_off(body, g)
        T = 2 * v / g
        out[w] = dict(take_off_m_s=r(v), airtime_s=r(T), shapes=r(T / SHAPE_S, 1),
                      take_off_can_give=take_off_turning(body, v, g),
                      one_somersault_needs_N_m_s={k: r(I * 2 * math.pi / T, 1) for k, I in body['inertia_kg_m2'].items()
                                                  if k in ('layout', 'tuck')},
                      somersault_rev_s=r(1 / T, 2))
    out['note'] = 'the centre of mass keeps its arc: shape changes orientation and turn rate; contact alone changes the arc'
    return out


# --- landings --------------------------------------------------------------------------------------------------------

def attitude_error_deg(dL, inertia, g, take_off=APEX_TAKE_OFF_M_S):
    """Rotation by the landing that an angular momentum dL picked up at the top of a jump leaves: dL / I over the
    descent's take_off / g."""
    return math.degrees(dL / inertia * take_off / g)


def landings(body):
    I = body['inertia_kg_m2']
    rows = []
    for v in (1.5, 2.5, 3.0):
        row = dict(landing_m_s=v)
        for d in (0.2, 0.4, 0.6):
            row[f'stroke_{d}_m'] = dict(stance_s=r(2 * d / v), **{w: dict(
                mean_push_own_bw=r((v * v / (2 * d) + g) / g, 1), mean_push_earth_bw=r((v * v / (2 * d) + g) / G_E, 2),
                brakes_or_turns_m_s=r(MU * (v + g * 2 * d / v), 2)) for w, g in WORLDS})
        rows.append(row)
    drops = [dict(drop_m=h, lands_m_s=r(math.sqrt(2 * G_M * h), 2), like_an_earth_drop_of_m=r(h * G_M / G_E, 2))
             for h in (1.0, 2.0, 3.0, 5.0)]
    errors = [dict(angular_momentum_N_m_s=dL, **{w: dict(layout_deg=r(attitude_error_deg(dL, I['layout'], g), 0),
                                                         tuck_deg=r(attitude_error_deg(dL, I['tuck'], g), 0))
                                                 for w, g in WORLDS}) for dL in (1.0, 3.0, 5.0, 10.0)]
    return dict(loads=rows, drops=drops, attitude_after_a_bump_at_the_top=dict(
        take_off_m_s=APEX_TAKE_OFF_M_S, descent_s={w: r(APEX_TAKE_OFF_M_S / g) for w, g in WORLDS}, rows=errors))


# --- balls and throws -----------------------------------------------------------------------------------------------

def drag_area(name):
    m, d = BALLS[name]
    if d is None:
        return 2 * m * G_E / (RHO_EARTH * SHUTTLE_EARTH_TERMINAL_M_S ** 2)
    return CD_SPHERE * math.pi * d * d / 4


def ball_flights(name, v0, g, rho, h0=RELEASE_HEIGHT_M, dt=2e-3):
    """Every release angle from 5 to 85 degrees at once, with quadratic drag, from h0 to the floor (semi-implicit
    Euler). Returns the best range with its angle, flight time and arrival speed."""
    m, _ = BALLS[name]
    k = 0.5 * rho * drag_area(name) / m
    ang = np.radians(np.arange(5.0, 85.5, 1.0))
    vx, vz = v0 * np.cos(ang), v0 * np.sin(ang)
    x, z = np.zeros_like(ang), np.full_like(ang, h0)
    t = 0.0
    done = np.zeros(ang.shape, bool)
    land = np.zeros_like(ang)
    tland = np.zeros_like(ang)
    vland = np.zeros_like(ang)
    while not done.all():
        s = np.hypot(vx, vz)
        vx = vx - k * s * vx * dt
        vz = vz - (g + k * s * vz) * dt
        x, z, t = x + vx * dt, z + vz * dt, t + dt
        hit = (~done) & (z <= 0)
        land[hit], tland[hit], vland[hit] = x[hit], t, np.hypot(vx[hit], vz[hit])
        done |= hit
        if t > 120:
            break
    i = int(np.argmax(land))
    return dict(best_angle_deg=r(np.degrees(ang[i]), 0), range_m=r(land[i], 1), flight_s=r(tland[i], 2),
                arrives_m_s=r(vland[i], 1))


def balls(rho):
    out = {}
    for name, (m, d) in BALLS.items():
        A = drag_area(name)
        term = {w: r(math.sqrt(2 * m * g / (rr * A)), 2) for (w, g), rr in zip(WORLDS, (rho, RHO_EARTH))}
        out[name] = dict(mass_kg=m, diameter_m=d, drag_area_m2=r(A, 5), terminal_m_s=term, throws=[
            dict(release_m_s=v, moon=ball_flights(name, v, G_M, rho), earth=ball_flights(name, v, G_E, RHO_EARTH),
                 moon_without_air=ball_flights(name, v, G_M, 0.0)) for v in RELEASE_M_S])
    return out


def throws(body, rho):
    """The throw's push back on the thrower against what the support can hold, and its turn in the air."""
    m_ball = BALLS[THROW_BALL][0]
    M = body['mass_kg']
    lever = body['release_above_com_m']
    h_release = body['com_height_m'] + lever
    rows = []
    for v in RELEASE_M_S:
        F = m_ball * v / THROW_S
        row = dict(release_m_s=v, hand_force_N=r(F, 1))
        for w, g in WORLDS:
            N = M * g
            row[w] = dict(friction_share_braced=r(F / (MU * N), 2),
                          reach_needed_recoiling_m=r(F * lever / N, 2), reach_needed_braced_m=r(F * h_release / N, 2),
                          held_feet_under=F * lever <= N * body['fore_aft_reach_m'])
        rows.append(row)
    # The fastest release whose push the support holds without turning the thrower: the floor's push N at the reach
    # carries the moment F lever, so v = N reach t / (m_ball lever). In a landing (2.5 m/s over the stroke) the floor
    # pushes with the landing's mean push instead of the weight.
    held = {}
    for stance, reach, landing in (('feet under the body', body['fore_aft_reach_m'], False),
                                   ('long split stance', body['split_reach_m'], False),
                                   ('in a landing, feet under', body['fore_aft_reach_m'], True)):
        held[stance] = {}
        for w, g in WORLDS:
            N = M * (2.5 ** 2 / (2 * STROKE_M) + g) if landing else M * g
            held[stance][w] = r(N * reach * THROW_S / (m_ball * lever), 1)
    air = []
    for label, arm in AIR_LEVERS_M.items():
        for v in (5.0, 10.0):
            dL = m_ball * v * arm
            air.append(dict(throw=label, release_m_s=v, angular_momentum_N_m_s=r(dL, 2),
                            **{w: dict(layout_deg=r(attitude_error_deg(dL, body['inertia_kg_m2']['layout'], g), 0),
                                       tuck_deg=r(attitude_error_deg(dL, body['inertia_kg_m2']['tuck'], g), 0))
                               for w, g in WORLDS}))
    leg_L = body['leg']['inertia_about_hip_kg_m2'] * KICK_RAD_S
    kick = dict(leg_angular_momentum_N_m_s=r(leg_L, 1), mean_moment_N_m=r(leg_L / KICK_BUILD_S, 0),
                one_foot_can_hold_N_m={w: r(M * g * body['fore_aft_reach_m'], 1) for w, g in WORLDS},
                note='what one foot cannot hold, the arms and trunk must take by turning against the leg')
    return dict(ball=THROW_BALL, ball_kg=m_ball, push_s=THROW_S, release_above_com_m=lever, standing=rows,
                fastest_release_held_m_s=held, in_the_air=air,
                catches='a catch takes the same impulse in reverse: the tables hold for a ball caught at that speed',
                kick=kick)


# --- contact between players --------------------------------------------------------------------------------------

def contact(bods):
    names = list(bods)
    masses = {n: bods[n]['mass_kg'] for n in names}
    pairs = [(names[0], names[0]), (names[0], names[1]), (names[1], names[1])]
    bumps = []
    for a, b in pairs:
        ma, mb = masses[a], masses[b]
        for v in (1.0, 2.0, 3.0):
            dva, dvb = mb / (ma + mb) * v, ma / (ma + mb) * v
            row = dict(bodies=[a, b], closing_m_s=v, speed_change_m_s=[r(dva, 2), r(dvb, 2)],
                       energy_J=r(0.5 * ma * mb / (ma + mb) * v * v, 0))
            for w, g in WORLDS:
                row[w] = dict(slide_m=[r(dva ** 2 / (2 * MU * g), 2), r(dvb ** 2 / (2 * MU * g), 2)],
                              to_absorb_by_grip_s=[r(dva / (MU * g), 2), r(dvb / (MU * g), 2)])
            bumps.append(row)
    pushes = {n: {w: r(MU * masses[n] * g, 0) for w, g in WORLDS} for n in names}
    air_exchange = []
    for a, b in pairs:
        ma, mb = masses[a], masses[b]
        for v in (0.5, 1.0):
            dva = mb / (ma + mb) * v
            air_exchange.append(dict(bodies=[a, b], closing_m_s=v, speed_change_m_s=r(dva, 2),
                                     landing_moves_m={w: r(dva * APEX_TAKE_OFF_M_S / g, 2) for w, g in WORLDS}))
    tosses = {}
    for n in names:
        for w, g in WORLDS:
            a = TOSS_FORCE_N / masses[n] - g                 # the flyer's net upward acceleration over the stroke
            v = math.sqrt(2 * a * TOSS_STROKE_M) if a > 0 else None
            tosses.setdefault(n, {})[w] = (dict(take_off_m_s=r(v), rise_m=r(v * v / (2 * g)), airtime_s=r(2 * v / g))
                                           if v else None)
    return dict(friction=MU, ground_collisions=bumps, standing_push_N=pushes,
                bumps_in_the_air_change_the_arcs=air_exchange,
                toss=dict(push_N=TOSS_FORCE_N, stroke_m=TOSS_STROKE_M, flyers=tosses),
                note='a bump in the air also turns the bodies; landings holds the attitude it leaves at the landing')


# --- getting down, getting up, falling ------------------------------------------------------------------------------

def floor_and_falls(body):
    h = body['com_height_m']
    out = dict(topple_s={f'{a}_to_{b}_deg': {w: r(dance.topple_seconds(a, b, g, h), 2) for w, g in WORLDS}
                         for a, b in ((2, 15), (2, 35), (5, 90))})
    out['get_down_s'] = {f'{d}_m': {w: r(dip_s(d, g)) for w, g in WORLDS} for d in (0.3, 0.45)}
    up = {}
    for w, g in WORLDS:
        v = take_off(body, g)
        a_up = v * v / (2 * body['push_off_m'])
        up[w] = {f'{d}_m': r(rise_and_stop_s(d, a_up, g)) for d in (0.3, 0.45)}
        up[w]['push_acceleration_m_s2'] = r(a_up, 1)
    out['stand_up_and_stop_s'] = up
    out['note'] = ('a change of height that must end at rest takes about twice as long as on Earth at either end: a body '
                   'sinks no faster than g, and a rise slows no faster than g; one that ends in motion (a jump, a drop '
                   'into a landing) is quick at its muscle end')
    return out


def results():
    rho = dance.air_density(0.0)
    bods = bodies()
    per_body = {n: dict(ground=ground(b), jumps=jumps(b), flight=flight(b), landings=landings(b), throws=throws(b, rho),
                        floor_and_falls=floor_and_falls(b)) for n, b in bods.items()}
    return dict(
        schema=SCHEMA,
        producer=dict(domain='humanities', files={f: digest(f) for f in FILES},
                      constants=constants_used([ROOT / f for f in FILES]),
                      inputs={str(p.relative_to(ROOT)): digest(p.relative_to(ROOT)) for p in INPUTS}),
        evidence='Closed-form and one-dimensional mechanics for two representative bodies (the dance body model and the '
                 'same body scaled to 1.80 m and 80 kg) at the shared gravity and the design run\'s air, with stated '
                 'parameter choices beside each figure. A screen of what a body can do for sport; the whole-body motion '
                 'engine tests it, and the whole-body record (whole_body.json) measures what its solved motions did.',
        gravity=dict(moon=G_M, earth=G_E, ratio=r(G_E / G_M, 4)),
        air_density_kg_m3=dict(moon_sea_level=r(rho, 4), earth=RHO_EARTH),
        bodies=bods, by_body=per_body, balls=balls(rho), contact=contact(bods),
        microgravity=dict(
            ground='no weight presses the feet: no grip, no stop, no turn without a surface to push against',
            jumps='every push sets an endless drift; nothing returns a player to play',
            balls='a ball flies straight until the air or a wall stops it; there is no arc to read',
            falls='there is no falling, landing or getting up, and no down to play against'))


def main() -> int:
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
