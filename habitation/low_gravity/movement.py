"""Building for people at a sixth of Earth's gravity: guards, stairs, handholds, ramps and movement halls.

    python -m habitation.low_gravity.movement    # writes habitation/low_gravity/results/low_gravity.json at once

Closed-form mechanics at the shared lunar gravity (GM/R^2 from shared/constants.json), set against Earth's building
codes and measurements of how people move, both recorded in sources.json. The design note (README.md) reads the
product and marks which figures are derived and which are design choices.

Assumptions, each used where its output is made:
- A hop or jump leaves the ground with the same speed as on Earth (the daily-life study's "same push"), so it rises
  g_E/g_M times as high and lasts that many times as long, and a hop of a given height on Earth stands for its
  take-off speed.
- Steady running is the exception: its vertical take-off speed scales with the square root of gravity (Polet et al.
  2018), so its bounce stays about Earth's height while its flight lasts sqrt(g_E/g_M) times as long.
- Shoes and tyres grip a dry floor with a friction coefficient of 0.6 (as in the daily-life study).
- A hand holds as much as its grip strength, and arms and body take up a caught fall over 0.3 m.
- Falls are taken without air drag. Below 30 m drag lowers the lunar impact speed by under 7%; long falls end at
  the daily-life study's 17.9 m/s.
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
RATIO = G_EARTH / G_MOON                  # 6.04: how much higher the same push lifts a body

# Earth's codes (sources.json). Heights and drops in metres.
EARTH_GUARD_M = 1.1                      # Approved Document K, diagram 3.1: landings, balconies, edges (UK)
EARTH_GUARD_IBC_M = 1.067                # IBC 1015.3: 42 inches
EARTH_FLIGHT_GUARD_M = 0.9               # Approved Document K, diagram 3.1: along flights of stairs
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
                  uk_gradient_by_going_m={'1:12': 2.0, '1:15': 5.0, '1:20': 10.0})
OPENING_SPHERE_M = 0.100                 # UK 100 mm; IBC 1015.4, 4 inches

# Measurements of people (sources.json).
JUMPS_EARTH_M = {                        # rise of the centre of mass in a standing jump, Markovic et al. 2004
    'countermovement jump': 0.353,       # 93 physical-education students, men, flight time on a contact mat
    'jump with an arm swing': 0.522,     # the same men, Abalakow jump with arm swing
}
WALK_RUN_MOON_M_S = 1.42                 # preferred walk-run transition at lunar gravity, De Witt et al. 2014
GRIP_KGF = {                             # Dodds et al. 2014, twelve British studies
    "women's peak median": 31.0, "men's peak median": 51.0,
    "women's weak grip (T-score -2.5)": 16.0, "men's weak grip (T-score -2.5)": 27.0,
}

# Assumptions (module docstring) and the example body.
FRICTION = 0.6
ARM_GIVE_M = 0.3
BODY_KG = 75.0                           # a round adult mass for the worked examples
TALL_ADULT_M = 1.9                       # stature for headroom in the movement halls
FOOT_M = 0.3                             # a shoe's length, added to stopping distances on landings

# The design choices this note makes (README.md gives the reasons).
GUARD_HOP_EARTH_M = (2.2 - EARTH_GUARD_M) / (RATIO - 1)   # the hop the summit drawings' 2.2 m guard covers
HEADROOM_HOP_EARTH_M = 0.10
EVERYDAY_STAIR = dict(rise_m=0.150, going_m=0.400, risers_a_bound=(2, 4, 6))
BOUNDING_STAIR = dict(rise_m=0.300, going_m=0.600, risers_a_bound=(1, 2, 4))
FLIGHT_RISE_CHOICE_M = 12.0
LANDING_M = 3.0
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


def drop(height_m, g=G_MOON):
    """A fall without drag: impact speed, the Earth drop with the same impact, and the time it lasts."""
    v = math.sqrt(2 * g * height_m)
    return dict(drop_m=height_m, impact_m_s=round(v, 2), like_earth_m=round(v ** 2 / (2 * G_EARTH), 3),
                seconds=round(math.sqrt(2 * height_m / g), 2))


def guard_height(hop_earth_m, earth_guard_m=EARTH_GUARD_M):
    """The guard that keeps Earth's margin over the centre of mass at the top of a hop that rises hop_earth_m on
    Earth: the same take-off lifts the body (RATIO - 1) * hop_earth_m higher at lunar gravity."""
    return earth_guard_m + (RATIO - 1) * hop_earth_m


def stop_m(speed_m_s, g=G_MOON, friction=FRICTION, slope=0.0):
    """Braking distance from speed_m_s on a floor of this slope (rise over run), going down it."""
    theta = math.atan(slope)
    return speed_m_s ** 2 / (2 * g * (friction * math.cos(theta) - math.sin(theta)))


def turn_radius(speed_m_s, g=G_MOON, friction=FRICTION):
    """The tightest level turn grip allows at this speed."""
    return speed_m_s ** 2 / (friction * g)


def bank_deg(speed_m_s, radius_m, g=G_MOON, friction=FRICTION):
    """Banking a floor needs so that grip holds a turn of this radius at this speed: tan(beta + atan mu) = v^2/(g r)."""
    return max(0.0, math.degrees(math.atan(speed_m_s ** 2 / (g * radius_m)) - math.atan(friction)))


def runaway_run_out(slope, friction=FRICTION):
    """Level run-out per metre of ramp that a wheel needs to stop after rolling the whole run unbraked, braking on
    the level landing at the foot with its full grip: (2 g sin(theta) L) / (2 friction g) = L sin(theta) / friction.
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
    hops = []
    for name, h in [('a hop that rises 5 cm on Earth', 0.05), ('a hop that rises 10 cm on Earth', 0.10),
                    ('the hop the 2.2 m guard covers', GUARD_HOP_EARTH_M)] + list(JUMPS_EARTH_M.items()):
        v = take_off(h)
        hops.append(dict(name=name, earth_rise_m=round(h, 4), take_off_m_s=round(v, 2), moon_rise_m=round(rise(v), 2),
                         earth_airtime_s=round(airtime(v, G_EARTH), 2), moon_airtime_s=round(airtime(v), 2),
                         guard_above_its_top_m=round(guard_height(h), 2)))

    falls = [drop(h) for h in (1.0, 2.0, 3.0, 4.6, 10.0, 30.0)]
    thresholds = {name: dict(earth_drop_m=h, moon_drop_same_impact_m=round(h * RATIO, 2),
                             moon_seconds=round(math.sqrt(2 * h * RATIO / G_MOON), 2),
                             earth_seconds=round(math.sqrt(2 * h / G_EARTH), 2))
                  for name, h in UNGUARDED_DROP_M.items()}
    strictest = min(UNGUARDED_DROP_M.values()) * RATIO
    loosest = max(UNGUARDED_DROP_M.values()) * RATIO

    guards = dict(
        rule='H = H_earth + (g_E/g_M - 1) j: Earth\'s guard, raised by the extra height a hop that rises j on Earth '
             'gains at lunar gravity, keeps Earth\'s margin over the centre of mass at the hop\'s top',
        by_hop=[dict(earth_hop_m=j, uk_based_m=round(guard_height(j), 2),
                     ibc_based_m=round(guard_height(j, EARTH_GUARD_IBC_M), 2)) for j in (0.0, 0.05, 0.10, 0.15, 0.20,
                                                                                       0.25, 0.30, 0.353, 0.522)],
        summit_drawings=dict(height_m=2.2, covers_earth_hop_m=round(GUARD_HOP_EARTH_M, 3),
                             take_off_m_s=round(take_off(GUARD_HOP_EARTH_M), 2),
                             share_of_countermovement_jump=round(
                                 GUARD_HOP_EARTH_M / JUMPS_EARTH_M['countermovement jump'], 2)),
        hands_rise_in_that_hop_m=round(rise(take_off(GUARD_HOP_EARTH_M)), 2),
        push_and_impact=dict(
            body_kg=BODY_KG,
            standing_push_n=dict(moon=round(FRICTION * BODY_KG * G_MOON), earth=round(FRICTION * BODY_KG * G_EARTH)),
            impact_energy_j={f'{v:g} m/s': round(0.5 * BODY_KG * v ** 2) for v in (WALK_RUN_MOON_M_S, 2.0, 3.0)}),
        classes=[
            dict(drops_m=[0.0, round(strictest, 2)], guard_m=None, kind='design choice',
                 why='lands no harder than the strictest Earth code\'s unguarded 380 mm'),
            dict(drops_m=[round(strictest, 2), round(loosest, 2)], guard_m=EARTH_GUARD_M, kind='design choice',
                 why='lands like the 0.38-0.76 m drops Earth\'s codes guard in some places and not in others'),
            dict(drops_m=[round(loosest, 2), None], guard_m=2.2, kind='design choice',
                 why='keeps Earth\'s margin over the centre of mass at the top of a hop that rises 0.22 m on Earth, '
                     'with a top that gives no grip'),
        ],
        opening_sphere_m=OPENING_SPHERE_M)

    everyday, bounding = stair(**EVERYDAY_STAIR), stair(**BOUNDING_STAIR)
    flight = {name: dict(earth_rise_m=round(h, 3), moon_rise_same_tumble_m=round(h * RATIO, 1))
              for name, h in EARTH_FLIGHT_RISE_M.items()}
    stops = {f'{v:g} m/s': dict(moon_m=round(stop_m(v), 2), earth_m=round(stop_m(v, G_EARTH), 2),
                                moon_turn_radius_m=round(turn_radius(v), 2),
                                earth_turn_radius_m=round(turn_radius(v, G_EARTH), 2))
             for v in (1.0, WALK_RUN_MOON_M_S, 2.0, 3.0)}
    headroom = {name: round(h + (RATIO - 1) * HEADROOM_HOP_EARTH_M, 2) for name, h in EARTH_HEADROOM_M.items()}
    stairs = dict(everyday=dict(everyday, kind='design choice'), bounding=dict(bounding, kind='design choice'),
                  riser_landing_like_earth_170mm_m=round(EARTH_STAIR['general_access_rise_m'][1] * RATIO, 2),
                  flight_limits=flight, flight_rise_choice_m=FLIGHT_RISE_CHOICE_M,
                  landing_m=LANDING_M, landing_need_from_2_m_s_m=round(stop_m(2.0) + FOOT_M, 2),
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

    by_slope = {}
    for name, s in (('1:12', 1 / 12), ('1:15', 1 / 15), ('1:20', 1 / 20)):
        theta = math.atan(s)
        out = runaway_run_out(s)
        by_slope[name] = dict(
            share_of_grip=round(s / FRICTION, 3), run_out_per_m_of_run=round(out, 4),
            longest_run_for_the_landing_m=round(LANDING_M / out, 1),
            rise_of_that_run_m=round(LANDING_M / out * math.sin(theta), 2),
            stop_on_the_ramp_from_walking_m=dict(moon=round(stop_m(WALK_RUN_MOON_M_S, slope=s), 2),
                                                 earth=round(stop_m(WALK_RUN_MOON_M_S, G_EARTH, slope=s), 2)),
            earth_slope_with_the_same_climb=f'1:{1 / math.tan(math.asin(math.sin(theta) / RATIO)):.0f}')
    steep = {f'{deg:g} degrees': dict(share_of_grip=round(math.tan(math.radians(deg)) / FRICTION, 2),
                                      stop_from_walking_m=dict(
                                          moon=round(stop_m(WALK_RUN_MOON_M_S, slope=math.tan(math.radians(deg))), 1),
                                          earth=round(stop_m(WALK_RUN_MOON_M_S, G_EARTH,
                                                             slope=math.tan(math.radians(deg))), 1)))
             for deg in (10.0, 18.4, 25.0, 30.0)}
    ramps = dict(by_slope=by_slope, steep_paths=steep, landing_m=LANDING_M,
                 run_rise_for_the_landing_m=round(LANDING_M * FRICTION, 2),
                 earth=dict(EARTH_RAMP,
                            ada_landing_stops_runs_rising_m=round(EARTH_RAMP['ada_landing_m'] * FRICTION, 2)),
                 lunar_slope_with_the_climb_of_earth_1_12_deg=round(math.degrees(math.asin(
                     RATIO * math.sin(math.atan(1 / 12)))), 1),
                 kind='design choice for the landing and the runs; the gradients are Earth\'s')

    hall_jump = JUMPS_EARTH_M[HALL_JUMP]
    hall_rise = rise(take_off(hall_jump))
    cmj_v = take_off(JUMPS_EARTH_M['countermovement jump'])
    need = TALL_ADULT_M + hall_rise + HALL_ARMS_M
    halls = dict(
        jump=HALL_JUMP, jump_rise_m=round(hall_rise, 2), head_top_m=round(TALL_ADULT_M + hall_rise, 2),
        needed_m=round(need, 2), clear_height_m=math.ceil(need * 2) / 2, kind='design choice',
        running_hop=dict(along_m_s=3.0, up_m_s=round(cmj_v, 2), airtime_s=round(airtime(cmj_v), 2),
                         length_m=round(3.0 * airtime(cmj_v), 1),
                         earth_length_m=round(3.0 * airtime(cmj_v, G_EARTH), 1)),
        banking_deg={f'{v:g} m/s at {r:g} m': round(bank_deg(v, r), 1)
                     for v, r in ((WALK_RUN_MOON_M_S, 1.0), (2.0, 2.0), (3.0, 3.0), (3.0, 5.0), (5.0, 10.0))},
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
                  'Earth\'s take-off speeds, steady running follows Polet et al. (2018), grip on a dry floor is 0.6 '
                  'and a hand holds its grip strength: assumptions, stated in movement.py. Entries marked "design '
                  'choice" are this note\'s choices. No human-factors test at lunar gravity stands behind any of '
                  'them; the longest stay on the Moon is 75 hours.'),
        reading_rule=('Lengths in metres, speeds in m/s, times in seconds, energies per kilogram of body in J/kg. '
                      '"like_earth" is the Earth drop or step with the same impact speed. Guard heights are measured '
                      'from the walking surface or, along stairs, from the nosings\' line. Falls are without drag.'),
        gravity=dict(moon_m_s2=round(G_MOON, 4), earth_m_s2=G_EARTH, ratio=round(RATIO, 3)),
        hops=hops, falls=falls, unguarded_drops=thresholds, guards=guards, stairs=stairs, stops=stops,
        handholds=handholds, ramps=ramps, movement_halls=halls, stepping_stones=stepping)


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
