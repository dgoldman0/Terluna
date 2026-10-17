"""Building for people at a sixth of Earth's gravity: guards, stairs, handholds, ramps and movement halls.

    python -m habitation.low_gravity.movement    # writes habitation/low_gravity/results/low_gravity.json at once

Closed-form mechanics at the shared lunar gravity (GM/R^2 from shared/constants.json), set against Earth's building
codes and measurements of how people move. The cited values are held here, each beside its source in sources.json,
and a test checks them against that file. The design note (README.md) reads the product and marks which figures are
derived and which are design choices.

Assumptions, each used where its output is made:
- A hop or jump leaves the ground at its Earth take-off speed, so it rises g_E/g_M times as high after take-off and
  stays aloft that many times as long. The legs lift a sixth of the weight, so the same effort can leave faster: the
  same leg work over the push-off bounds the take-off from above (jump_bounds), and Earth's take-off speed sets the
  jump figures low.
- Steady running is the exception: its vertical take-off speed scales with the square root of gravity (Polet et al.
  2018), so its bounce stays about Earth's height while its flight lasts sqrt(g_E/g_M) times as long.
- Falls follow the daily-life study's quadratic drag, which ends a long lunar fall at 17.9 m/s; on Earth the short
  drops compared here are taken without drag.
- Shoes and tyres grip a dry floor with a friction coefficient of 0.6 (as in the daily-life study).
- A body's standing centre of mass sits at 0.55 of its stature.
- A hand holds as much as its grip strength, and arms and body take up a caught fall over 0.3 m.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.habitation.low-gravity/1'
OUT = HERE / 'results' / 'low_gravity.json'
FILES = ('habitation/low_gravity/movement.py',)

G_MOON, G_EARTH = MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
RATIO = G_EARTH / G_MOON                  # 6.04: how much higher the same take-off speed lifts a body

# Earth's codes (sources.json). Heights and drops in metres.
EARTH_GUARD_M = 1.1                      # Approved Document K, diagram 3.1: landings, balconies, edges (UK)
EARTH_GUARD_IBC_M = 1.067                # IBC 1015.3: 42 inches
UNGUARDED_DROP_M = {                     # the largest drop each code leaves without a guard
    'UK, buildings other than dwellings (two risers, or 380 mm)': 0.38,
    'UK, dwellings': 0.60,
    'IBC 1015.2 (30 inches)': 0.762,
}
EARTH_STAIR = dict(general_access_rise_m=(0.150, 0.170), general_access_going_m=(0.250, 0.400),
                   two_r_plus_g_m=(0.550, 0.700), ibc_rise_m=(0.102, 0.178), ibc_tread_min_m=0.279)
EARTH_FLIGHT_RISE_M = {'UK general access, 12 risers of 170 mm': 12 * 0.170, 'IBC 1011.8 (12 feet)': 3.658}
EARTH_HEADROOM_M = {'UK escape routes': 2.0, 'IBC (80 inches)': 2.032}
EARTH_HANDRAIL_M = {'UK (Approved Document K)': (0.900, 1.000), 'ADA 505.4 (34-38 inches)': (0.865, 0.965)}
EARTH_RAMP = dict(ada_max_slope=1 / 12, ada_max_rise_m=0.760, ada_landing_m=1.525,
                  adk_max_going_m={'1:12': 2.0, '1:15': 5.0, '1:20': 10.0})
OPENING_SPHERE_M = 0.100                 # UK 100 mm; IBC 1015.4, 4 inches

# Measurements of people (sources.json).
MARKOVIC_M = dict(                       # Markovic et al. 2004, 93 male physical-education students
    countermovement_after_take_off=0.353,  # flight time on a contact mat: the rise after take-off
    abalakow_without_arm_swing=0.462,    # waist belt and tape, from standing
    abalakow_with_arm_swing=0.522,       # the same, with an arm swing
)
# The belt reads the push-off's extension as well as the flight: for the same countermovement jump it reads 0.109 m
# above the mat. The extension is the same on both worlds; only the rise after take-off scales.
PUSH_OFF_EXTENSION_M = MARKOVIC_M['abalakow_without_arm_swing'] - MARKOVIC_M['countermovement_after_take_off']
JUMPS_EARTH_M = {                        # the rise after take-off
    'countermovement jump': MARKOVIC_M['countermovement_after_take_off'],
    'jump with an arm swing': MARKOVIC_M['abalakow_with_arm_swing'] - PUSH_OFF_EXTENSION_M,
}
WALK_RUN_MOON_M_S = 1.42                 # preferred walk-run transition at lunar gravity, De Witt et al. 2014
GRIP_KGF = {                             # Dodds et al. 2014, twelve British studies
    "women's peak median": 31.0, "men's peak median": 51.0,
    "women's weak grip (T-score -2.5)": 16.0, "men's weak grip (T-score -2.5)": 27.0,
}
TERMINAL_M_S = 17.9                      # the daily-life study's end of a long lunar fall, with drag

# Assumptions (module docstring) and the example bodies.
FRICTION = 0.6
ARM_GIVE_M = 0.3
BODY_KG = 75.0                           # a round adult mass for the worked examples
COM_SHARE = 0.55                         # standing centre of mass over stature
STATURES_M = (1.91, 1.75, 1.62, 1.40, 1.20, 1.00)
EXAMPLE_ADULT_M = 1.75
PUSH_OFF_DEPTH_M = (0.35, 0.45)          # the centre of mass's travel through a jump's push-off
TALL_ADULT_M = 1.9                       # stature for headroom in the movement halls
FOOT_M = 0.3                             # a shoe's length, added to stopping distances on landings

# The design choices this note makes (README.md gives the reasons).
GRIP_SHARE = 0.5                         # design turns, landings and ramps use half the grip
LUNAR_GUARD_M = 2.2                      # the summit drawings' railing, kept
HEADROOM_HOP_EARTH_M = 0.10
EVERYDAY_STAIR = dict(rise_m=0.150, going_m=0.400, risers_a_bound=(2, 4, 6))
BOUNDING_STAIR = dict(rise_m=0.300, going_m=0.600, risers_a_bound=(1, 2, 4))
FLIGHT_RISE_CHOICE_M = 12.0
LANDING_M = dict(everyday=dict(depth_m=3.0, arrival_m_s=WALK_RUN_MOON_M_S),
                 bounding=dict(depth_m=4.5, arrival_m_s=2.0))
RAMP_LANDING_M = 3.0
FOOTHOLDS_M = {'a kerb': 0.15, 'a bench seat or a planter\'s rim': 0.45, 'a free-standing handrail': 1.0}
HALL_JUMP = 'jump with an arm swing'
HALL_ARMS_M = 1.0                        # raised arms and clearance above the head at a jump's top
STEPPING_HOP_M_S = dict(along=1.4, up=1.0)


def rise(take_off_m_s, g=G_MOON):
    """Height a body's centre of mass rises after leaving the ground at this vertical speed."""
    return take_off_m_s ** 2 / (2 * g)


def take_off(rise_m, g=G_EARTH):
    """Vertical take-off speed that lifts the centre of mass by rise_m."""
    return math.sqrt(2 * g * rise_m)


def airtime(take_off_m_s, g=G_MOON):
    """Time from take-off back to the same level."""
    return 2 * take_off_m_s / g


def same_work_take_off(earth_take_off_m_s, depth_m):
    """The lunar take-off speed the same leg work gives over a push-off of this depth: the legs' work pays the speed
    and the lift through the push-off, W/m = v^2/2 + g d, on both worlds. Muscle force falls with speed, so this is an
    upper bound."""
    return math.sqrt(earth_take_off_m_s ** 2 + 2 * (G_EARTH - G_MOON) * depth_m)


def impact_speed(height_m, g=G_MOON, terminal=TERMINAL_M_S):
    """Impact speed from a fall of height_m with quadratic drag."""
    return terminal * math.sqrt(1 - math.exp(-2 * g * height_m / terminal ** 2))


def drop_for_impact(speed_m_s, g=G_MOON, terminal=TERMINAL_M_S):
    """The fall with drag that lands at this speed; the inverse of impact_speed."""
    return -terminal ** 2 / (2 * g) * math.log(1 - speed_m_s ** 2 / terminal ** 2)


def fall_seconds(height_m, g=G_MOON, terminal=TERMINAL_M_S):
    """Time a fall of height_m lasts with quadratic drag: h = (vt^2/g) ln cosh(g t / vt)."""
    return terminal / g * math.acosh(math.exp(g * height_m / terminal ** 2))


def drop(height_m):
    """A lunar fall with drag: impact speed, the Earth drop with the same impact, and the time it lasts."""
    v = impact_speed(height_m)
    return dict(drop_m=height_m, impact_m_s=round(v, 2), like_earth_m=round(v ** 2 / (2 * G_EARTH), 3),
                seconds=round(fall_seconds(height_m), 2))


def guard_height(hop_earth_m, earth_guard_m=EARTH_GUARD_M):
    """The hop rule: the guard that keeps the margin Earth's guard has at the top of the same hop, one that rises
    hop_earth_m after take-off on Earth: the same take-off lifts the body (RATIO - 1) * hop_earth_m higher here."""
    return earth_guard_m + (RATIO - 1) * hop_earth_m


def energy_guard(com_m, earth_guard_m=EARTH_GUARD_M):
    """The energy rule: the guard whose top costs the same lift per kilogram as Earth's, g_M (H - c) = g_E (H_E - c)."""
    return com_m + RATIO * (earth_guard_m - com_m)


def lift_speed(top_m, com_m, g=G_MOON):
    """The impact speed whose whole kinetic energy pays the lift of the centre of mass from com_m over top_m."""
    return math.sqrt(2 * g * max(top_m - com_m, 0.0))


def foothold_reach(height_m, take_off_m_s, g=G_MOON):
    """How far a foothold of this height has to stand from a guard to give a hop at this take-off speed no higher
    crossing of the guard's plane than a hop from the floor at the guard. The paths of a hop at speed v reach no higher
    than y0 + v^2/2g - g x^2 / 2v^2 at a distance x, so the foothold's height s is spent at x = v sqrt(2 s / g)."""
    return take_off_m_s * math.sqrt(2 * height_m / g)


def stop_m(speed_m_s, g=G_MOON, friction=FRICTION, slope=0.0):
    """Braking distance from speed_m_s on a floor of this slope (rise over run), going down it."""
    theta = math.atan(slope)
    return speed_m_s ** 2 / (2 * g * (friction * math.cos(theta) - math.sin(theta)))


def turn_radius(speed_m_s, g=G_MOON, friction=FRICTION):
    """The tightest level turn this grip allows at this speed."""
    return speed_m_s ** 2 / (friction * g)


def bank_deg(speed_m_s, radius_m, g=G_MOON, friction=FRICTION):
    """Banking a floor needs so that this grip holds a turn of this radius at this speed:
    tan(beta + atan mu) = v^2/(g r)."""
    return max(0.0, math.degrees(math.atan(speed_m_s ** 2 / (g * radius_m)) - math.atan(friction)))


def runaway_run_out(slope, friction=FRICTION):
    """Level run-out per metre of ramp that a wheel needs to stop after rolling the whole run unbraked, braking on
    the level landing at the foot with this grip: (2 g sin(theta) L) / (2 friction g) = L sin(theta) / friction.
    Gravity cancels, since the speed gained and the braking both scale with it."""
    return math.sin(math.atan(slope)) / friction


def stair(rise_m, going_m, risers_a_bound=(2, 4, 6)):
    """A stair at lunar gravity against Earth's general-access stair at its steepest (170 mm rise). A bound clears
    k risers, its top at the landing tread's level."""
    earth_rise = EARTH_STAIR['general_access_rise_m'][1]
    bounds = []
    for k in risers_a_bound:
        up = take_off(k * rise_m, G_MOON)
        t_top = up / G_MOON
        bounds.append(dict(risers=k, rise_m=round(k * rise_m, 3), take_off_up_m_s=round(up, 2),
                           seconds_to_top=round(t_top, 2), speed_along_m_s=round(k * going_m / t_top, 2),
                           like_earth_hop_m=round(k * rise_m / RATIO, 3)))
    return dict(
        rise_m=rise_m, going_m=going_m, pitch_deg=round(math.degrees(math.atan(rise_m / going_m)), 1),
        two_r_plus_g_m=round(2 * rise_m + going_m, 3),
        work_per_riser_j_kg=round(G_MOON * rise_m, 3),
        earth_rise_with_the_same_work_m=round(rise_m / RATIO, 4),
        step_down=dict(seconds=round(math.sqrt(2 * rise_m / G_MOON), 2),
                       landing_m_s=round(math.sqrt(2 * G_MOON * rise_m), 2),
                       like_earth_step_m=round(rise_m / RATIO, 4)),
        earth_general_access=dict(rise_m=earth_rise, work_per_riser_j_kg=round(G_EARTH * earth_rise, 3),
                                  step_down_seconds=round(math.sqrt(2 * earth_rise / G_EARTH), 2),
                                  step_down_landing_m_s=round(math.sqrt(2 * G_EARTH * earth_rise), 2)),
        bounds=bounds)


def arrest_m(grip_n, hands, body_kg=BODY_KG, g=G_MOON, give_m=ARM_GIVE_M):
    """The fall a grip can stop over the arms' give: hands * grip = m g (1 + h / give). Negative: it cannot hold the
    body's weight at all."""
    return give_m * (hands * grip_n / (body_kg * g) - 1)


def results():
    design_grip = GRIP_SHARE * FRICTION
    sized_hop = (LUNAR_GUARD_M - EARTH_GUARD_M) / (RATIO - 1)       # the hop the 2.2 m guard keeps Earth's margin at
    sized_take_off = take_off(sized_hop)
    adult_com = COM_SHARE * EXAMPLE_ADULT_M

    hops = []
    for name, h in [('a hop that rises 5 cm on Earth', 0.05), ('a hop that rises 10 cm on Earth', 0.10),
                    ('the hop the 2.2 m guard is sized for', sized_hop)] + list(JUMPS_EARTH_M.items()):
        v = take_off(h)
        hops.append(dict(name=name, earth_rise_m=round(h, 4), take_off_m_s=round(v, 2), moon_rise_m=round(rise(v), 2),
                         earth_airtime_s=round(airtime(v, G_EARTH), 2), moon_airtime_s=round(airtime(v), 2),
                         guard_keeping_earths_margin_at_its_top_m=round(guard_height(h), 2)))
    jump_bounds = {}
    for name, h in JUMPS_EARTH_M.items():
        v = take_off(h)
        extension = PUSH_OFF_EXTENSION_M if name == 'jump with an arm swing' else 0.0
        jump_bounds[name] = dict(
            earth_take_off_m_s=round(v, 2),
            same_take_off_speed=dict(take_off_m_s=round(v, 2), moon_rise_m=round(extension + rise(v), 2)),
            same_leg_work=[dict(push_off_m=d, take_off_m_s=round(same_work_take_off(v, d), 2),
                                moon_rise_m=round(extension + rise(same_work_take_off(v, d)), 2))
                           for d in PUSH_OFF_DEPTH_M],
            rise_measured_from='standing' if extension else 'take-off')

    falls = [drop(h) for h in (1.0, 2.0, 3.0, 10.0, 30.0)]
    thresholds = {}
    for name, h in UNGUARDED_DROP_M.items():
        v = math.sqrt(2 * G_EARTH * h)
        moon_h = drop_for_impact(v)
        thresholds[name] = dict(earth_drop_m=h, impact_m_s=round(v, 2), moon_drop_same_impact_m=round(moon_h, 2),
                                moon_seconds=round(fall_seconds(moon_h), 2),
                                earth_seconds=round(math.sqrt(2 * h / G_EARTH), 2))
    strictest = drop_for_impact(math.sqrt(2 * G_EARTH * min(UNGUARDED_DROP_M.values())))
    loosest = drop_for_impact(math.sqrt(2 * G_EARTH * max(UNGUARDED_DROP_M.values())))

    energy = []
    for s in STATURES_M:
        c = COM_SHARE * s
        energy.append(dict(stature_m=s, com_m=round(c, 3),
                           earth_lift_j_kg=round(G_EARTH * (EARTH_GUARD_M - c), 2),
                           earth_lift_speed_m_s=round(lift_speed(EARTH_GUARD_M, c, G_EARTH), 2),
                           moon_guard_same_lift_m=round(energy_guard(c), 2),
                           moon_lift_over_2_2_m_j_kg=round(G_MOON * (LUNAR_GUARD_M - c), 2),
                           moon_lift_speed_over_2_2_m_m_s=round(lift_speed(LUNAR_GUARD_M, c), 2)))
    # The lowest centre of mass for which 2.2 m costs at least Earth's lift.
    com_kept = (RATIO * EARTH_GUARD_M - LUNAR_GUARD_M) / (RATIO - 1)
    child_rows = [r for r in energy if r['moon_guard_same_lift_m'] > LUNAR_GUARD_M]
    guards = dict(
        rules=dict(
            standing_and_leaning='Tipping over a guard turns on where the centre of mass stands against its top, which '
                                 'gravity leaves alone: Earth\'s height keeps Earth\'s margin.',
            hops='H = H_Earth + (g_E/g_M - 1) j keeps the margin Earth\'s guard has at the top of the same hop, '
                 'one that rises j after take-off on Earth.',
            stumbling_and_running='H = c + (g_E/g_M) (H_Earth - c) costs the same lift per kilogram over the top as '
                                  'Earth\'s guard, so that the same impact speed carries the centre of mass over it on '
                                  'both worlds if the same share of the impact\'s energy goes into the lift.'),
        by_hop=[dict(earth_hop_m=j, moon_rise_m=round(RATIO * j, 2), uk_based_m=round(guard_height(j), 2),
                     ibc_based_m=round(guard_height(j, EARTH_GUARD_IBC_M), 2))
                for j in (0.0, 0.05, 0.10, 0.20, round(sized_hop, 4)) + tuple(JUMPS_EARTH_M.values())],
        sized_hop=dict(guard_m=LUNAR_GUARD_M, earth_hop_m=round(sized_hop, 4), take_off_m_s=round(sized_take_off, 2),
                       share_of_countermovement_jump=round(sized_hop / JUMPS_EARTH_M['countermovement jump'], 2),
                       example_adult=dict(stature_m=EXAMPLE_ADULT_M, com_m=round(adult_com, 3),
                                          com_over_top_at_the_hops_top_m=dict(
                                              earth=round(adult_com + sized_hop - EARTH_GUARD_M, 3),
                                              moon=round(adult_com + RATIO * sized_hop - LUNAR_GUARD_M, 3))),
                       hands_rise_m=round(rise(sized_take_off), 2)),
        by_stature=energy,
        energy_kept=dict(guard_m=LUNAR_GUARD_M, lowest_com_m=round(com_kept, 3),
                         shortest_stature_m=round(com_kept / COM_SHARE, 2),
                         children_need_m=[min(r['moon_guard_same_lift_m'] for r in child_rows),
                                          max(r['moon_guard_same_lift_m'] for r in child_rows)],
                         for_statures_m=[min(r['stature_m'] for r in child_rows),
                                         max(r['stature_m'] for r in child_rows)]),
        push_and_impact=dict(
            body_kg=BODY_KG,
            standing_push_n=dict(moon=round(FRICTION * BODY_KG * G_MOON), earth=round(FRICTION * BODY_KG * G_EARTH)),
            impact_energy_j={f'{v:g} m/s': round(0.5 * BODY_KG * v ** 2) for v in (WALK_RUN_MOON_M_S, 2.0, 3.0)}),
        bands=dict(kind='design choice, for the author\'s review', rows=[
            dict(drops_m=[0.0, round(strictest, 2)], guard=None,
                 why='lands no harder than the strictest Earth code\'s unguarded 380 mm'),
            dict(drops_m=[round(strictest, 2), round(loosest, 2)], options_m=[EARTH_GUARD_M, LUNAR_GUARD_M],
                 why='lands like the 0.38-0.76 m drops Earth\'s codes guard in some places and not others; Earth\'s '
                     '1.1 m keeps Earth\'s margin for standing and leaning and costs stumbling and running a sixth of '
                     'Earth\'s lift, 2.2 m keeps the longer drops\' guard'),
            dict(drops_m=[round(loosest, 2), None], guard_m=LUNAR_GUARD_M,
                 why='keeps Earth\'s lift for every adult whose centre of mass stands at the lowest_com_m of '
                     'energy_kept or higher, and Earth\'s margin at the top of the sized hop'),
            dict(drops_m=[round(loosest, 2), None], where='where young children play',
                 options=['guards of children_need_m in energy_kept', 'play kept back from the edge'])]),
        footholds=dict(take_off_m_s=round(sized_take_off, 2), rule='the guard\'s height counts from any surface a foot '
                       'can use that stands closer to it than reach_m', reach_m={
                           name: dict(height_m=s, reach_m=round(foothold_reach(s, sized_take_off), 2))
                           for name, s in FOOTHOLDS_M.items()}),
        opening_sphere_m=OPENING_SPHERE_M)

    everyday, bounding = stair(**EVERYDAY_STAIR), stair(**BOUNDING_STAIR)
    flight = {name: dict(earth_rise_m=round(h, 3), moon_rise_same_tumble_m=round(h * RATIO, 1))
              for name, h in EARTH_FLIGHT_RISE_M.items()}
    stops = {f'{v:g} m/s': dict(moon_m=round(stop_m(v), 2), earth_m=round(stop_m(v, G_EARTH), 2),
                                moon_at_half_grip_m=round(stop_m(v, friction=design_grip), 2),
                                moon_turn_radius_m=round(turn_radius(v), 2),
                                moon_turn_radius_at_half_grip_m=round(turn_radius(v, friction=design_grip), 2),
                                earth_turn_radius_m=round(turn_radius(v, G_EARTH), 2))
             for v in (1.0, WALK_RUN_MOON_M_S, 2.0, 3.0)}
    headroom = {name: round(h + (RATIO - 1) * HEADROOM_HOP_EARTH_M, 2) for name, h in EARTH_HEADROOM_M.items()}
    landings = {kind: dict(spec, needed_m=round(stop_m(spec['arrival_m_s'], friction=design_grip) + FOOT_M, 2))
                for kind, spec in LANDING_M.items()}
    stairs = dict(everyday=dict(everyday, kind='design choice'), bounding=dict(bounding, kind='design choice'),
                  riser_landing_like_earth_170mm_m=round(EARTH_STAIR['general_access_rise_m'][1] * RATIO, 2),
                  flight_limits=flight, flight_rise_choice_m=FLIGHT_RISE_CHOICE_M,
                  landings=dict(landings, grip_share=GRIP_SHARE, kind='design choice'),
                  headroom_m=dict(headroom, hop_earth_m=HEADROOM_HOP_EARTH_M, choice_m=2.5, kind='design choice'))

    grips = {name: dict(kgf=f, newtons=round(f * G_EARTH),
                        lunar_body_kg_one_hand=round(f * RATIO),
                        arrest_m=dict(moon_one_hand=round(arrest_m(f * G_EARTH, 1), 2),
                                      moon_two_hands=round(arrest_m(f * G_EARTH, 2), 2),
                                      earth_two_hands=round(arrest_m(f * G_EARTH, 2, g=G_EARTH), 2)))
             for name, f in GRIP_KGF.items()}
    handholds = dict(body_kg=BODY_KG, give_m=ARM_GIVE_M, weight_n=dict(moon=round(BODY_KG * G_MOON),
                                                                       earth=round(BODY_KG * G_EARTH)),
                     grips=grips, handrail_heights_m=EARTH_HANDRAIL_M)

    run_rise = design_grip * RAMP_LANDING_M
    by_slope = {}
    for name, s in (('1:12', 1 / 12), ('1:15', 1 / 15), ('1:20', 1 / 20)):
        theta = math.atan(s)
        by_slope[name] = dict(
            share_of_grip=round(s / FRICTION, 3), run_out_per_m_of_run=round(runaway_run_out(s), 4),
            run_for_the_landing_m=round(run_rise / math.sin(theta), 1),
            adk_max_going_m=EARTH_RAMP['adk_max_going_m'][name],
            stop_on_the_ramp_from_walking_m=dict(moon=round(stop_m(WALK_RUN_MOON_M_S, slope=s), 2),
                                                 earth=round(stop_m(WALK_RUN_MOON_M_S, G_EARTH, slope=s), 2)),
            earth_slope_with_the_same_climb=f'1:{1 / math.tan(math.asin(math.sin(theta) / RATIO)):.0f}')
    steep = {f'{deg:g} degrees': dict(share_of_grip_up_or_down=round(math.tan(math.radians(deg)) / FRICTION, 2),
                                      stop_going_down_from_walking_m=dict(
                                          moon=round(stop_m(WALK_RUN_MOON_M_S, slope=math.tan(math.radians(deg))), 1),
                                          earth=round(stop_m(WALK_RUN_MOON_M_S, G_EARTH,
                                                             slope=math.tan(math.radians(deg))), 1)))
             for deg in (10.0, 18.4, 25.0, 30.0)}
    ramps = dict(by_slope=by_slope, steep_paths=steep,
                 steepest_within_half_grip_deg=round(math.degrees(math.atan(design_grip)), 1),
                 landing_m=RAMP_LANDING_M, grip_share=GRIP_SHARE, run_rise_for_the_landing_m=round(run_rise, 2),
                 earth=dict(EARTH_RAMP,
                            ada_landing_stops_runs_rising_m=round(EARTH_RAMP['ada_landing_m'] * FRICTION, 2)),
                 lunar_slope_with_the_climb_of_earth_1_12_deg=round(math.degrees(math.asin(
                     RATIO * math.sin(math.atan(1 / 12)))), 1),
                 kind='design choice for the landing and the runs; the gradients are Earth\'s')

    hall_v = take_off(JUMPS_EARTH_M[HALL_JUMP])
    hall_rise = PUSH_OFF_EXTENSION_M + rise(hall_v)
    need = TALL_ADULT_M + hall_rise + HALL_ARMS_M
    bound_rises = [PUSH_OFF_EXTENSION_M + rise(same_work_take_off(hall_v, d)) for d in PUSH_OFF_DEPTH_M]
    cmj_v = take_off(JUMPS_EARTH_M['countermovement jump'])
    cmj_bound = [same_work_take_off(cmj_v, d) for d in PUSH_OFF_DEPTH_M]
    halls = dict(
        jump=HALL_JUMP, jump_rise_from_standing_m=round(hall_rise, 2), head_top_m=round(TALL_ADULT_M + hall_rise, 2),
        needed_m=round(need, 2), clear_height_m=math.ceil(need * 2) / 2, kind='design choice, provisional',
        needed_at_the_same_leg_work_m=[round(TALL_ADULT_M + r + HALL_ARMS_M, 2) for r in bound_rises],
        running_hop=dict(along_m_s=3.0, up_m_s=round(cmj_v, 2), airtime_s=round(airtime(cmj_v), 2),
                         length_m=round(3.0 * airtime(cmj_v), 1),
                         earth_length_m=round(3.0 * airtime(cmj_v, G_EARTH), 1),
                         at_the_same_leg_work=dict(up_m_s=[round(v, 2) for v in cmj_bound],
                                                   length_m=[round(3.0 * airtime(v), 1) for v in cmj_bound])),
        grip_share=GRIP_SHARE,
        banking_deg={f'{v:g} m/s at {r:g} m': round(bank_deg(v, r, friction=design_grip), 1)
                     for v, r in ((WALK_RUN_MOON_M_S, 1.0), (2.0, 2.0), (3.0, 3.0), (3.0, 5.0), (5.0, 10.0))},
        level_turn_radius_m={f'{v:g} m/s': round(turn_radius(v, friction=design_grip), 2)
                             for v in (WALK_RUN_MOON_M_S, 2.0, 3.0)},
        running_airtime_ratio=round(math.sqrt(RATIO), 2))

    u, w = STEPPING_HOP_M_S['along'], STEPPING_HOP_M_S['up']
    stepping = dict(along_m_s=u, up_m_s=w, moon_length_m=round(u * airtime(w), 2),
                    earth_length_m=round(u * airtime(w, G_EARTH), 2), moon_rise_m=round(rise(w), 2),
                    like_earth_hop_m=round(rise(w, G_EARTH), 3), kind='design choice')

    return dict(
        schema=SCHEMA,
        producer=dict(domain='habitation', files={f: digest(ROOT / f) for f in FILES}, constants=constants_used(FILES)),
        evidence=('Closed-form mechanics at the shared lunar gravity, set against Earth\'s building codes and '
                  'measurements of people moving on Earth and in simulated lunar gravity (sources.json). Hops keep '
                  'Earth\'s take-off speeds, which sets jumps low, with the same leg work as an upper bound; steady '
                  'running follows Polet et al. (2018); falls take the daily-life study\'s drag; a centre of mass sits '
                  'at 0.55 of stature; grip on a dry floor is 0.6 and a hand holds its grip strength: assumptions, '
                  'stated in movement.py. Entries marked "design choice" are this note\'s choices, for the author\'s '
                  'review. No human-factors test at lunar gravity stands behind any of them; the longest stay on the '
                  'Moon is 75 hours.'),
        reading_rule=('Lengths in metres, speeds in m/s, times in seconds, energies per kilogram of body in J/kg. '
                      '"like_earth" is the Earth drop or step with the same impact speed. A hop\'s rise is its rise '
                      'after take-off unless marked from standing. Guard heights are measured from the walking '
                      'surface, along stairs from the nosings\' line, and from any foothold closer than its reach. '
                      'Physical limits use the full grip of 0.6; design turns, landings and ramps use half of it.'),
        gravity=dict(moon_m_s2=round(G_MOON, 4), earth_m_s2=G_EARTH, ratio=round(RATIO, 3)),
        hops=hops, jump_bounds=jump_bounds, falls=falls, unguarded_drops=thresholds, guards=guards, stairs=stairs,
        stops=stops, handholds=handholds, ramps=ramps, movement_halls=halls, stepping_stones=stepping)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def main() -> int:
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, ensure_ascii=False) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
