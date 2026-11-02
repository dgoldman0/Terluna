"""What a dancer's body, a partner and a garment can do at lunar gravity in the Open Moon's air.

    python -m humanities.dance.envelope     # writes humanities/dance/results/envelope.json at once

The first screen of the dance domain. Each figure sets the Moon beside Earth, and where it helps beside microgravity,
for the same body: the 1.64 m, 55 kg trained female dancer of the body model (humanities/dance/body/model.py). The
Moon's gravity is the shared GM/R^2; its air is the design run's at the sea-level datum (1.40 kg/m^3) and at the summit
metropolis (1.17 kg/m^3), read from the GCM's global winds product; Earth's air is the standard 1.225 kg/m^3.

Assumptions, each where its figures are made:
- Jumps. The legs' push follows Samozino's linear force-velocity model of the lower limb, F = F0 (1 - v/V0) per
  kilogram, over a push-off of h_po; F0 and V0 are the profile measured in female handball players (Palinkas et al.
  2025), whose countermovement jumps rose 27-28 cm, and h_po = 0.36 m is assumed. The same muscles pushing against a
  sixth of the weight give the lunar take-off. Earth's take-off speed and the same leg work bracket it below and above.
- Rotation. The body model's inertia in two poses (standing, arms down; arms overhead) gives the layout and twist
  moments; the tuck and pike take Yeadon's ratios to the layout. A turn on the floor gets its torque from two feet
  pushing tangentially at the friction limit, 0.40 m apart.
- Figures in the air. The body model posed in the vocabulary's shapes (body/pose.py) gives each figure's inertia, and
  the computed tuck and pike shares are reported beside Yeadon's. A take-off gives a somersault or cartwheel the mean
  push over the push-off times an assumed lean of 0.25 m past the foot, and a spin about the vertical the floor's grip
  on that push times an assumed lever of 0.30 m. Sleeves are the garment solver's 1.5 m of chiffon.
- Friction between dance shoe and floor: 0.5 for steps and turns, an assumed value for leather and suede on a vinyl
  dance floor, for which no measurement was found (sources.json, open); the habitation note's 0.6 stays its own.
- Balance. A standing dancer topples like an inverted pendulum about the ankle with the centre of mass at its standing
  height.
- Air. Drag areas: the body alone across a vertical flow 0.20 m^2, with a flared full skirt 0.80 m^2, with large
  sleeves or a cloak spread 2.0 m^2 (assumed); flat fabric falls face-on with a normal-force coefficient of 1.17.
- Fabrics: silk chiffon 26.1 g/m^2 and georgette 34.3 g/m^2 as measured by Fuss (2024); crepe de chine, satin and
  linen are stated brackets. The garment solver (garment/) carries the full cloth mechanics.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from habitation.low_gravity.movement import stop_m, turn_radius, bank_deg
from humanities.dance.body.model import MASS_KG, STATURE_M, spec as body_spec
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.humanities.dance-envelope/1'
OUT = HERE / 'results' / 'envelope.json'
FILES = ('humanities/dance/envelope.py', 'humanities/dance/body/model.py', 'humanities/dance/body/pose.py',
         'habitation/low_gravity/movement.py')
POSES = ROOT / 'humanities' / 'dance' / 'vocabulary' / 'poses.json'
WINDS = ROOT / 'climate' / 'results' / 'gcm' / 'global_winds_A28_dim5_moon.json'

G_M = MOON_SURFACE_GRAVITY
G_E = STANDARD_GRAVITY
RHO_EARTH = 1.225
SUMMIT_KM = 11.6                          # the summit metropolis's ground (research/studies/sky_fleet)

# Jumping: Samozino's force-velocity model, set to an Earth countermovement jump of a trained woman.
F0_N_KG = 29.0                            # maximal force per kilogram at zero speed (Palinkas et al. 2025: 28-29)
V0_M_S = 3.05                             # unloaded extension speed (Palinkas et al. 2025: 3.0-3.1 m/s)
H_PO_M = 0.36                             # push-off distance, crouch to take-off (assumed; no number found)
# Rotation: Yeadon's tuck and pike as shares of the layout's somersault inertia.
TUCK_SHARE = 0.30
PIKE_SHARE = 0.55
FEET_APART_M = 0.40                       # fourth position's preparation for a turn
PIVOT_RADIUS_M = 0.02                     # effective radius of the ball of the foot's contact on demi-pointe
FRICTION = 0.5                            # dance shoe on a vinyl dance floor
SPIN_FRICTION = 0.4                       # turning friction on the ball of the foot
# Air and cloth
DRAG_AREA_M2 = {'body alone, vertical flow': 0.20, 'a flared full skirt': 0.80, 'large sleeves or a cloak spread': 2.0}
CN_FLAT = 1.17
FABRIC_KG_M2 = {'silk chiffon': 0.026, 'silk georgette': 0.034, 'crepe de chine': 0.065, 'silk satin': 0.090,
                'linen': 0.150}          # chiffon and georgette measured (Fuss 2024); the rest stated brackets
SKIRT = dict(waist_radius_m=0.12, length_m=0.70)
# Partnering
THROW_FORCE_N = 350.0                     # two-arm push beyond the thrower's own support, over the throw's stroke
THROW_STROKE_M = 0.60
LIFT_ARM_FORCE_N = 250.0                  # a strong one-arm overhead press


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16]


def air_density(height_km: float) -> float:
    w = json.loads(WINDS.read_text())
    assert w['schema'] == 'terluna.climate.gcm-global-winds/1'
    return float(np.interp(height_km, w['height_km'], w['density_mean_kg_m3']))


# --- jumps ---------------------------------------------------------------------------------------------------------

def samozino_take_off(g, f0=F0_N_KG, v0=V0_M_S, h=H_PO_M):
    """Take-off speed in Samozino's model: v^2 = 2 h (F0 (1 - v/(2 V0)) - g), the mean force at the mean speed v/2."""
    b = f0 * h / v0
    c = -2 * h * (f0 - g)
    return (-b + math.sqrt(b * b - 4 * c)) / 2


def ode_take_off(g, f0=F0_N_KG, v0=V0_M_S, h=H_PO_M, n=20000):
    """The same push integrated exactly: v dv/dx = F0 (1 - v/V0) - g from rest over h (a check on the mean-force
    form; it counts the force at each speed rather than at the mean)."""
    v, dx = 1e-6, h / n
    for _ in range(n):
        a = f0 * (1 - v / v0) - g
        v = math.sqrt(max(v * v + 2 * a * dx, 1e-12))
    return v


def flight(v_up, g):
    return dict(take_off_m_s=round(v_up, 3), rise_m=round(v_up ** 2 / (2 * g), 3), airtime_s=round(2 * v_up / g, 3))


def jumps():
    earth = samozino_take_off(G_E)
    moon = samozino_take_off(G_M)
    same_work = math.sqrt(earth ** 2 + 2 * (G_E - G_M) * H_PO_M)
    return dict(
        model=dict(F0_N_kg=F0_N_KG, V0_m_s=V0_M_S, push_off_m=H_PO_M),
        earth=flight(earth, G_E),
        moon_force_velocity=flight(moon, G_M),
        moon_earth_take_off=flight(earth, G_M),
        moon_same_work=flight(same_work, G_M),
        ode_check=dict(earth_m_s=round(ode_take_off(G_E), 3), moon_m_s=round(ode_take_off(G_M), 3)),
        choreographic=[dict(share_of_max_rise=s, **flight(math.sqrt(s) * moon, G_M)) for s in (0.25, 0.5, 0.75)])


# --- inertia ---------------------------------------------------------------------------------------------------------

def whole_body_inertia(arms_up=False, arms_side=False):
    """The whole body's inertia about its centre of mass, standing, from the body model's segments (parallel-axis
    sums of each segment's principal moments). arms_up mirrors the arm segments about the shoulder joints' height;
    arms_side lays them out along +-y at shoulder height (second position)."""
    sp = body_spec()
    by = {s['name']: s for s in sp['segments']}

    def world_origin(name):
        s = by[name]
        return np.array(s['origin'], float)

    items = []
    for s in sp['segments']:
        m = s['mass_kg']
        c = world_origin(s['name']) + np.array(s['com'])
        I = np.diag(s['inertia_kg_m2'])
        arm = any(k in s['name'] for k in ('upper_arm', 'forearm', 'hand'))
        if arm and (arms_up or arms_side):
            side = 1 if s['name'].startswith('L') else -1
            sh = world_origin(f'{"L" if side > 0 else "R"}_upper_arm')
            drop = sh[2] - c[2]
            if arms_up:
                c = np.array([c[0], c[1], sh[2] + drop])
            else:
                c = np.array([c[0], sh[1] + side * drop, sh[2]])
                I = np.diag([I[2, 2], I[0, 0], I[1, 1]])      # the long axis turns from z to y
        items.append((m, c, I))
    M = sum(m for m, _, _ in items)
    com = sum(m * c for m, c, _ in items) / M
    total = np.zeros((3, 3))
    for m, c, I in items:
        d = c - com
        total += I + m * (np.dot(d, d) * np.eye(3) - np.outer(d, d))
    return M, com, total


def rotations(jump):
    _, com, I_down = whole_body_inertia()
    _, _, I_up = whole_body_inertia(arms_up=True)
    _, _, I_side = whole_body_inertia(arms_side=True)
    layout = I_up[1, 1]
    I = dict(somersault_layout_arms_up=layout, somersault_layout_arms_down=I_down[1, 1],
             somersault_pike=PIKE_SHARE * layout, somersault_tuck=TUCK_SHARE * layout,
             twist_arms_in=I_down[2, 2], twist_arms_side=I_side[2, 2])
    flights = dict(moon_half_height=jump['choreographic'][1]['airtime_s'],
                   moon_max=jump['moon_force_velocity']['airtime_s'], earth_max=jump['earth']['airtime_s'])
    out = dict(inertia_kg_m2={k: round(v, 3) for k, v in I.items()}, com_height_m=round(float(com[2]), 3),
               flights_s=flights, figures=[])
    for name, axis, turns in (('one layout somersault', 'somersault_layout_arms_up', 1),
                              ('one tucked somersault', 'somersault_tuck', 1),
                              ('a full twist, arms in', 'twist_arms_in', 1),
                              ('a full turn, arms in second', 'twist_arms_side', 1),
                              ('a triple twist, arms in', 'twist_arms_in', 3)):
        row = dict(figure=name, inertia_kg_m2=round(I[axis], 3))
        for where, t in flights.items():
            w = 2 * math.pi * turns / t
            row[where] = dict(rad_s=round(w, 2), angular_momentum_N_m_s=round(I[axis] * w, 1))
        out['figures'].append(row)
    return out


# --- the floor: traction, turns, balance -----------------------------------------------------------------------------

def traction():
    rows = []
    for v in (1.0, 1.5, 2.0, 3.0):
        rows.append(dict(speed_m_s=v,
                         build_or_stop_m=dict(moon=round(stop_m(v, G_M, FRICTION), 2),
                                              earth=round(stop_m(v, G_E, FRICTION), 2)),
                         tightest_turn_radius_m=dict(moon=round(turn_radius(v, G_M, FRICTION), 2),
                                                     earth=round(turn_radius(v, G_E, FRICTION), 2)),
                         bank_for_a_2_m_circle_deg=round(bank_deg(v, 2.0, G_M, FRICTION), 1)))
    return dict(friction=FRICTION, peak_acceleration_m_s2=dict(moon=round(FRICTION * G_M, 3),
                                                               earth=round(FRICTION * G_E, 3)), speeds=rows)


def landing_redirect(v_down, absorb_m, g):
    """Horizontal speed a landing can add or take away: the stance lasts T = 2 d / v as the legs absorb v over d,
    the floor pushes N = m (v / T + g) on average, and friction allows mu N T / m."""
    T = 2 * absorb_m / v_down
    n_per_kg = v_down / T + g
    return dict(stance_s=round(T, 3), mean_load_body_weights=round(n_per_kg / g, 2),
                redirect_m_s=round(FRICTION * n_per_kg * T, 3))


def turns_on_the_floor(rot):
    weight = MASS_KG * G_M
    tau = dict(moon=FRICTION * MASS_KG * G_M * FEET_APART_M / 2, earth=FRICTION * MASS_KG * G_E * FEET_APART_M / 2)
    I = rot['inertia_kg_m2']['twist_arms_in']
    rows = []
    for rps in (0.5, 1.0, 2.0):
        L = I * 2 * math.pi * rps
        row = dict(revolutions_per_s=rps, angular_momentum_N_m_s=round(L, 2))
        for where, g in (('moon', G_M), ('earth', G_E)):
            decay = SPIN_FRICTION * MASS_KG * g * PIVOT_RADIUS_M
            row[where] = dict(preparation_s=round(L / tau[where], 2), spin_decay_s=round(L / decay, 1),
                              revolutions_before_stopping=round(L / decay * rps / 2, 1))
        rows.append(row)
    return dict(turning_torque_N_m={k: round(v, 1) for k, v in tau.items()}, lunar_weight_N=round(weight, 1),
                pirouettes=rows)


def topple_seconds(theta0_deg, theta1_deg, g, length):
    """Inverted pendulum about the ankle: theta'' = (g/L) sin theta, from rest at theta0 to theta1."""
    th, w, t, dt = math.radians(theta0_deg), 0.0, 0.0, 1e-4
    while th < math.radians(theta1_deg):
        w += (g / length) * math.sin(th) * dt
        th += w * dt
        t += dt
    return t


def balance(rot):
    L = rot['com_height_m']
    out = dict(com_height_m=L, time_constant_s=dict(moon=round(math.sqrt(L / G_M), 3),
                                                     earth=round(math.sqrt(L / G_E), 3)))
    for a, b in ((2, 15), (2, 45), (5, 90)):
        out[f'from_{a}_to_{b}_deg_s'] = dict(moon=round(topple_seconds(a, b, G_M, L), 2),
                                             earth=round(topple_seconds(a, b, G_E, L), 2))
    return out


# --- air -----------------------------------------------------------------------------------------------------------

def drag_flight(v_up, g, rho, drag_area, mass=MASS_KG, dt=1e-4):
    """A vertical jump with quadratic drag, launched at v_up and landing at the launch height."""
    z, v, t, apex = 0.0, v_up, 0.0, 0.0
    k = 0.5 * rho * drag_area / mass
    while True:
        a = -g - k * v * abs(v)
        v += a * dt
        z += v * dt
        t += dt
        apex = max(apex, z)
        if z < 0 and v < 0:
            return dict(apex_m=round(apex, 3), airtime_s=round(t, 3), landing_m_s=round(-v, 3))


def air(jump):
    rho = dict(moon_sea_level=air_density(0.0), moon_summit=air_density(SUMMIT_KM), earth=RHO_EARTH)
    v_up = jump['choreographic'][1]['take_off_m_s']
    flights = []
    for name, area in DRAG_AREA_M2.items():
        flights.append(dict(garment=name, drag_area_m2=area,
                            moon=drag_flight(v_up, G_M, rho['moon_sea_level'], area),
                            moon_vacuum=flight(v_up, G_M),
                            earth=drag_flight(v_up, G_E, RHO_EARTH, area)))
    terminal = {}
    for name, area in (('body spread face-down', 0.70),) + tuple(DRAG_AREA_M2.items()):
        terminal[name] = dict(moon=round(math.sqrt(2 * MASS_KG * G_M / (rho['moon_sea_level'] * area)), 1),
                              earth=round(math.sqrt(2 * MASS_KG * G_E / (RHO_EARTH * area)), 1))
    return dict(density_kg_m3={k: round(v, 4) for k, v in rho.items()}, jump_take_off_m_s=v_up,
                jumps_with_drag=flights, terminal_speed_m_s=terminal)


def fabric(rho_moon):
    rows = []
    for name, ma in FABRIC_KG_M2.items():
        vm = math.sqrt(2 * ma * G_M / (rho_moon * CN_FLAT))
        ve = math.sqrt(2 * ma * G_E / (RHO_EARTH * CN_FLAT))
        rows.append(dict(fabric=name, areal_density_kg_m2=ma,
                         face_on_fall_m_s=dict(moon=round(vm, 3), earth=round(ve, 3)),
                         seconds_to_settle_1_5_m=dict(moon=round(1.5 / vm, 1), earth=round(1.5 / ve, 1)),
                         air_over_weight_at_1_m_s=dict(moon=round(0.5 * rho_moon * CN_FLAT / (ma * G_M), 2),
                                                       earth=round(0.5 * RHO_EARTH * CN_FLAT / (ma * G_E), 2))))
    r, l = SKIRT['waist_radius_m'], SKIRT['length_m']

    def rpm(theta_deg, g):
        th = math.radians(theta_deg)
        return 60 / (2 * math.pi) * math.sqrt(g * math.tan(th) / (r + 2 * l / 3 * math.sin(th)))

    flare = [dict(flare_deg=a, rpm=dict(moon=round(rpm(a, G_M), 1), earth=round(rpm(a, G_E), 1)))
             for a in (30, 45, 60, 75)]
    return dict(normal_force_coefficient=CN_FLAT, fabrics=rows,
                bending_length_ratio_moon_to_earth=round((G_E / G_M) ** (1 / 3), 3),
                skirt=dict(**SKIRT, rigid_strip_flare=flare))


# --- partnering --------------------------------------------------------------------------------------------------------

def partnering():
    m = MASS_KG
    out = dict(partner_weight_N=dict(moon=round(m * G_M, 1), earth=round(m * G_E, 1)),
               one_arm_overhead_hold=dict(moon=LIFT_ARM_FORCE_N >= m * G_M, earth=LIFT_ARM_FORCE_N >= m * G_E))
    for where, g in (('moon', G_M), ('earth', G_E)):
        a = THROW_FORCE_N / m - g
        out[f'throw_{where}'] = flight(math.sqrt(2 * a * THROW_STROKE_M), g) if a > 0 else None
    # A catch that keeps Earth's landing load: the catcher's arms and legs take the partner's speed over 0.6 m.
    return out


# --- figures in the air -------------------------------------------------------------------------------------------------

AERIAL_LEAN_M = 0.25      # the centre of mass's lean past the pushing foot at an acrobatic take-off (assumed)
SPIN_LEVER_M = 0.30       # the lever of the floor's sideways push about the centre of mass in a spinning take-off (assumed)
SLEEVE_LENGTH_M = 1.5     # the garment solver's water sleeves
AERIAL_SHAPES = {         # vocabulary shapes (poses.json) for the figures, over the motion engine's carriage
    'layout, arms overhead': dict(legs='layout_air', arms='overhead_straight'),
    'layout, arms gathered': dict(legs='layout_air', arms='tucked_in'),
    'tuck': dict(legs='tuck_air', arms='tucked_in', torso='contraction'),
    'pike': dict(legs='pike_air', arms='second'),
    'star': dict(legs='star_air', arms='high_v'),
    'flat spread': dict(legs='straddle_air', arms='second'),
}


def take_off_turning(v, g):
    """How much angular momentum a take-off at speed v can give the body, from the mean push over the push-off: the
    push's lean past the foot turns a somersault or cartwheel; the floor's grip, pushing sideways on a lever, turns the
    body about the vertical."""
    t = 2 * H_PO_M / v
    push = MASS_KG * (v / t + g)
    return dict(push_s=round(t, 3), mean_push_N=round(push, 1),
                somersault_or_cartwheel_N_m_s=round(push * AERIAL_LEAN_M * t, 1),
                about_the_vertical_N_m_s=round(FRICTION * push * SPIN_LEVER_M * t, 1))


def aerials(jump, rho_moon):
    """Rotations in the air that the Moon's long flights open, with the body model's own inertia in each shape: a
    twisting somersault, a cartwheel in a star, a flat spin, a somersault that hangs upside down, and sleeves left
    standing above a falling body. Each needs the angular momentum its rotations take over the flight, against what the
    take-off can give; Earth's figures are for its best standing jump."""
    from humanities.dance.body import pose
    I = {k: pose.inertia(s)[2] for k, s in AERIAL_SHAPES.items()}
    moon, high = jump['choreographic'][1], jump['choreographic'][2]
    earth = jump['earth']
    flights = dict(moon=(moon['take_off_m_s'], moon['airtime_s'], G_M), moon_high=(high['take_off_m_s'], high['airtime_s'], G_M),
                   earth=(earth['take_off_m_s'], earth['airtime_s'], G_E))
    cap = {k: take_off_turning(v, g) for k, (v, _, g) in flights.items()}

    def need(inertia, turns, where):
        return round(inertia * 2 * math.pi * turns / flights[where][1], 1)

    out = dict(lean_m=AERIAL_LEAN_M, spin_lever_m=SPIN_LEVER_M, flights_s={k: round(t, 3) for k, (_, t, _) in flights.items()},
               take_off=cap, inertia_kg_m2={k: dict(cartwheel=round(v[0, 0], 3), somersault=round(v[1, 1], 3),
                                                    twist=round(v[2, 2], 3)) for k, v in I.items()},
               tuck_and_pike_shares=dict(tuck=round(I['tuck'][1, 1] / I['layout, arms overhead'][1, 1], 3),
                                         pike=round(I['pike'][1, 1] / I['layout, arms overhead'][1, 1], 3),
                                         assumed=dict(tuck=TUCK_SHARE, pike=PIKE_SHARE)))
    # A layout somersault with a full twist: the somersault's angular momentum drives the twist once the body tilts
    # out of its plane (twist rate L sin(tilt) / I_twist), so a 1:1 corkscrew needs sin(tilt) = I_twist / I_somersault.
    Is, It = I['layout, arms gathered'][1, 1], I['layout, arms gathered'][2, 2]
    out['corkscrew'] = dict(
        tilt_deg=round(math.degrees(math.asin(It / Is)), 1),
        twist_rev_s={k: round(1 / flights[k][1], 2) for k in ('moon', 'earth')},
        angular_momentum_N_m_s={k: need(Is, 1, k) for k in ('moon', 'earth')})
    out['star_wheel'] = dict(rev_s={k: round(1 / flights[k][1], 2) for k in ('moon', 'earth')},
                             angular_momentum_N_m_s={k: need(I['star'][0, 0], 1, k) for k in ('moon', 'earth')})
    out['flat_spin'] = dict(angular_momentum_N_m_s={k: need(I['flat spread'][0, 0], 1, k) for k in flights},
                            spread_is_the_largest_axis=bool(np.argmax(np.linalg.eigvalsh(I['flat spread'])) == 2
                                                            and I['flat spread'][0, 0] >= I['flat spread'][1:, 1:].diagonal().max()))
    # A somersault that tucks to turn over, opens when upside down and tucks again: open 0-0.1 of the flight, tuck to
    # 0.3, open to 0.7 (the hang, centred on the inversion), tuck to 0.9, open to land.
    Io, It_ = I['layout, arms overhead'][1, 1], I['tuck'][1, 1]
    program = ((0.1, Io), (0.2, It_), (0.4, Io), (0.2, It_), (0.1, Io))
    hang = {}
    for k in ('moon', 'earth'):
        T = flights[k][1]
        L = 2 * math.pi / (T * sum(f / i for f, i in program))
        hang[k] = dict(angular_momentum_N_m_s=round(L, 1), hang_s=round(0.4 * T, 2),
                       hang_sweep_deg=round(math.degrees(L * 0.4 * T / Io), 1))
    out['hanging_salto'] = dict(open_over_tuck=round(Io / It_, 2), **hang)
    # Sleeves cast up at the take-off: from the apex the body falls freely while chiffon sinks at most at its face-on
    # speed, so the sleeves stand above the hands until the landing, up to their length, then settle.
    ma = FABRIC_KG_M2['silk chiffon']
    sleeve = {}
    for k, rho in (('moon', rho_moon), ('earth', RHO_EARTH)):
        v, T, g = flights[k]
        vt = math.sqrt(2 * ma * g / (rho * CN_FLAT))
        rise = v * v / (2 * g)
        gap = min(max(rise - vt * T / 2, 0.0), SLEEVE_LENGTH_M)
        sleeve[k] = dict(body_descent_s=round(T / 2, 2), landing_m_s=round(v, 2), sleeve_sink_m_s=round(vt, 3),
                         sleeves_above_hands_at_landing_m=round(gap, 2), settle_after_landing_s=round(gap / vt, 1))
    out['sleeve_flight'] = sleeve
    return out


def stage_for_a_leap(jump):
    v_up = jump['choreographic'][1]['take_off_m_s']
    rows = []
    for v_along in (1.0, 1.5, 2.0):
        t = 2 * v_up / G_M
        rows.append(dict(speed_along_m_s=v_along, approach_m=round(stop_m(v_along, G_M, FRICTION), 2),
                         flight_m=round(v_along * t, 2), run_out_m=round(stop_m(v_along, G_M, FRICTION), 2),
                         total_m=round(2 * stop_m(v_along, G_M, FRICTION) + v_along * t, 2),
                         earth_same_speeds_total_m=round(2 * stop_m(v_along, G_E, FRICTION)
                                                         + v_along * 2 * v_up / G_E, 2)))
    return dict(take_off_up_m_s=v_up, rows=rows)


def results():
    j = jumps()
    rot = rotations(j)
    a = air(j)
    return dict(
        schema=SCHEMA,
        producer=dict(domain='humanities', files={f: digest(f) for f in FILES},
                      constants=constants_used([ROOT / f for f in FILES]),
                      inputs={str(p.relative_to(ROOT)): digest(p.relative_to(ROOT)) for p in (WINDS, POSES)}),
        evidence='Closed-form and one-dimensional mechanics for one representative dancer (the body model) at the '
                 'shared gravity and the design run\'s air, with stated parameter choices beside each figure. A '
                 'screen of capabilities and limits; the whole-body motion engine and the garment solver test them.',
        body=dict(stature_m=STATURE_M, mass_kg=MASS_KG),
        gravity=dict(moon=G_M, earth=G_E, ratio=round(G_E / G_M, 4)),
        jumps=j, rotations=rot, traction=traction(),
        landings={f'{v}_m_s_over_0_4_m': dict(moon=landing_redirect(v, 0.4, G_M), earth=landing_redirect(v, 0.4, G_E))
                  for v in (1.5, 2.5)},
        turns=turns_on_the_floor(rot), balance=balance(rot), air=a, fabric=fabric(a['density_kg_m3']['moon_sea_level']),
        partnering=partnering(), stage_for_a_leap=stage_for_a_leap(j),
        aerials=aerials(j, a['density_kg_m3']['moon_sea_level']),
        microgravity=dict(
            traction='none: no weight presses the feet, so every push needs a surface to push against',
            flight='unending: a push sets a straight drift that only a surface or the air stops',
            falls_and_balance='absent: no toppling, no landing, no down to dance against',
            fabric='unsettled: cloth drifts until the air stops it, with no direction to hang or fall'))


def main() -> int:
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
