"""How fabrics and garment cuts act on the Moon: a closed-form screen of every fabric in the library against Earth.

    python -m humanities.dance.garment.screen          # writes humanities/dance/results/garments.json

Cloth under gravity in air answers to three groups of numbers, and the Moon moves each one differently.

1. Hanging, swinging, falling and settling are driven by gravity and resisted by the cloth's inertia, the air's
   added mass and its drag. Written in the time unit sqrt(L/g), a hanging panel's swing has no g left in it: its
   equation keeps only the fabric's mass against the air's (m / rho L) and the swing's shape. A lunar garment left to
   itself therefore moves as the same garment does on Earth, slowed by sqrt(g_E / g_M) = 2.46, the gravity clock the
   dancer's own flights keep, and only the 14% denser lunar air changes the shape of the motion.
2. Stiffness against weight sets the drape. Peirce's bending length c = (B / m g)^(1/3) is 1.82 times longer on the
   Moon, so a fabric drapes there as one with six times its bending rigidity drapes on Earth: folds are fewer and
   broader, and crisp fabrics hold their shapes.
3. Motions the dancer drives at the speed of muscle (a turn, a run, an arm's sweep) meet the air and the centrifugal
   pull at Earth's speeds while weight is a sixth: against the weight, the centrifugal pull counts 6.0 times more
   than on Earth and the air's push 6.9 times more (the 6.0 times the denser air's 1.14). Between the two forced
   loads the weight drops out: the ratio of a fabric's centrifugal pull to the air's drag in a turn is its drag length
   2m / (rho C_N) over its radius, whatever the gravity, so once the pull outweighs the weight a lunar skirt's shape
   no longer depends on how fast it turns.

The screen gives each fabric's figures under all three. For circle skirts in steady turns (full and half circles) it
gives a meridian string without air and with skin friction alone (a skirt without folds), and an index of the folds'
drag against the centrifugal pull, which holds neither the rate of turn nor gravity, with the weighted hem that brings
it down. A string with fold drag (validate.chain_estimate) finds no steady shape on the Moon: its iteration wanders
between a wound and an open skirt, so the trial runs (humanities/dance/garment/trials.py) decide that case.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from humanities.dance.envelope import G_E, G_M, RHO_EARTH, air_density
from humanities.dance.garment.cloth import Settings
from humanities.dance.garment.garments import BAND_CLEARANCE
from humanities.dance.garment.materials import BENDING_STATUS, FABRICS, bending_length
from humanities.dance.motion_io import body_spec
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.humanities.dance-garments/1'
OUT = ROOT / 'humanities' / 'dance' / 'results' / 'garments.json'
FILES = ('humanities/dance/garment/screen.py', 'humanities/dance/garment/materials.py',
         'humanities/dance/garment/cloth.py', 'humanities/dance/garment/garments.py', 'humanities/dance/envelope.py')
INPUTS = ('humanities/dance/body/body.json', 'humanities/dance/garment/sources.json',
          'climate/results/gcm/global_winds_A28_dim5_moon.json')

C_N = Settings.c_normal                   # the solver's normal-pressure and skin-friction coefficients
C_T = Settings.c_tangent
SKIRT_LENGTH_M = 0.85                     # the solver's skirts
PANEL = dict(length_m=0.85, width_m=0.40, release_deg=30.0)   # a hanging panel let go from 30 degrees
FOLD_SCALE_M = 0.20                       # a fold's width, for the air's added mass
GARMENT_SCALE_M = 1.0                     # a whole garment's
TRAVEL_M_S = (1.0, 2.0, 3.0)              # a walk, the orbit run, a bounding run
TURN_RPM = (15.0, 27.0, 60.0)             # a slow turn, the lunar whirl (0.45 rev/s), an Earth whirl
HEM_KG_M = (0.0, 0.025, 0.05, 0.1, 0.2)   # a weighted hem's mass per metre of hem


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16]


def worlds():
    return dict(moon=dict(g=G_M, rho=air_density(0.0)), earth=dict(g=G_E, rho=RHO_EARTH))


def waist_radius():
    ab = {s['name']: s for s in body_spec()['segments']}['abdomen']
    return ab['capsules'][0]['radius'] + BAND_CLEARANCE


# --- 1. the gravity clock: falling, settling and swinging -----------------------------------------------------------

def face_on_fall(m, g, rho, factor=1.0):
    """Terminal speed of a flat sheet falling face-on."""
    return math.sqrt(2 * m * g / (rho * C_N * factor))


def fall_time(d, vt, g):
    """Time to fall d from rest with quadratic drag: d = (vt^2 / g) ln cosh(g t / vt)."""
    return vt / g * math.acosh(math.exp(g * d / vt ** 2))


def added_mass(rho, a, factor=1.0):
    """Added mass per unit area of a strip of width a moving normal to itself: pi rho a^2 / 4 per unit span. An open
    weave or net lets air through and carries less; it is taken to carry the same share as the pressure it takes
    (factor), an assumption."""
    return math.pi * rho * a / 4 * factor


def swing(m, g, rho, factor=1.0, L=PANEL['length_m'], W=PANEL['width_m'], theta0_deg=PANEL['release_deg'], n=20000):
    """A rigid panel hanging from its top edge, let go at theta0, swinging face-on through still air (a Morison-type
    estimate with the strip's potential-flow added mass and the steady normal-pressure coefficient):
    (m + m_a) L^3/3 theta'' = -m g L^2/2 sin(theta) - rho C_N L^4/8 |theta'| theta'. In the time unit sqrt(L/g) the
    equation holds no g, so the Moon's swing is Earth's slowed by the gravity clock, with only the air density changed;
    that holds for any drag and added-mass coefficients that depend on the swing's shape (its Keulegan-Carpenter
    number) rather than its speed. Returns the air's added mass over the fabric's, the time to the first passage
    through the vertical, and the largest swing back past it, which is an upper bound: oscillating plates at these
    amplitudes take drag coefficients several times the steady one (Keulegan and Carpenter), so they swing back less."""
    ma = added_mass(rho, W, factor)
    I = (m + ma) * L ** 3 / 3
    k = rho * C_N * factor * L ** 4 / 8
    G = m * g * L ** 2 / 2

    def rhs(a, w):
        return w, (-G * math.sin(a) - k * abs(w) * w) / I

    th, om = math.radians(theta0_deg), 0.0
    dt = 8.0 * math.sqrt(L / g) * math.sqrt(1 + ma / m) / n
    first, back = None, 0.0
    for i in range(1, n + 1):
        a1, b1 = rhs(th, om)
        a2, b2 = rhs(th + 0.5 * dt * a1, om + 0.5 * dt * b1)
        a3, b3 = rhs(th + 0.5 * dt * a2, om + 0.5 * dt * b2)
        a4, b4 = rhs(th + dt * a3, om + dt * b3)
        th += dt * (a1 + 2 * a2 + 2 * a3 + a4) / 6
        om += dt * (b1 + 2 * b2 + 2 * b3 + b4) / 6
        if first is None and th < 0:
            first = i * dt
        if first is not None:
            if om >= 0 and th < 0:                      # the first swing back has reached its far point
                back = -th
                break
    return dict(added_mass_ratio=round(ma / m, 2),
                first_through_vertical_s=None if first is None else round(first, 2),
                largest_swing_back_deg=round(math.degrees(back), 1))


# --- 3. forced motion: travel and turns -----------------------------------------------------------------------------

def trail_angle(m, g, rho, U, factor=1.0):
    """A straight panel trailing from its leading edge in travel at U: the normal pressure 1/2 rho C_N U^2 sin^2(a)
    holds its weight's normal share m g cos(a) at an angle a below the horizontal (the solver's pressure law; a thin
    plate's attached flow at shallow incidence pushes harder and holds it flatter, so the angle is an upper bound)."""
    q = 0.5 * rho * C_N * factor * U ** 2
    w = m * g
    # q sin^2 a = w cos a  ->  q (1 - c^2) = w c  ->  c = (-w + sqrt(w^2 + 4 q^2)) / (2 q)
    c = (-w + math.sqrt(w * w + 4 * q * q)) / (2 * q)
    return math.degrees(math.acos(min(1.0, c)))


def meridian(m, g, rho, omega, cut=2 * math.pi, hem_kg_m=0.0, L=SKIRT_LENGTH_M, r_w=None, factor=1.0, folds=True,
             n=300, iters=6000, relax=0.15):
    """A circle skirt turning steadily with the waist in still air, as one meridian string per radian of the
    dancer's turn (extending humanities/dance/garment/validate.chain_estimate to other cuts and a weighted hem).
    Along the string, at slant length s from the waist, the cloth per radian of turn spans rc(s) = r_w + cut s / 2pi
    of pattern width and carries m rc of mass; the hem adds hem_kg_m rc(L). Loads: weight, centrifugal pull and the
    air's circumferential force. With folds, the cloth's excess width over the circle it sits on, rc against R, lies
    in a zig-zag whose faces meet the flow at beta, cos(beta) = R / rc (an upper bound on fold drag); without folds,
    skin friction only (a lower bound). The string lies along the sum of the loads below each point. Where a cut
    narrower than a full circle would open past the cone it makes when laid smooth (dR/ds > cut / 2pi), its hoop
    goes taut: the string is held at that slope, keeping the loads' circumferential share."""
    r_w = waist_radius() if r_w is None else r_w
    s = (np.arange(n) + 0.5) * L / n
    ds = L / n
    rc = r_w + cut * s / (2 * math.pi)
    lam = m * rc
    hem_m = hem_kg_m * (r_w + cut * L / (2 * math.pi))
    slope_max = cut / (2 * math.pi)
    X = np.c_[r_w + 0.3 * slope_max * s, np.zeros(n), -s * 0.95]
    res = np.inf
    t = None
    for it in range(iters):
        Rh = np.hypot(X[:, 0], X[:, 1])
        f = np.zeros((n, 3))
        f[:, 2] = -lam * g
        f[:, 0] = lam * omega ** 2 * X[:, 0]
        f[:, 1] = lam * omega ** 2 * X[:, 1]
        cosb = np.clip(Rh / rc, 0, 1)
        sinb = np.sqrt(1 - cosb ** 2)
        ceff = (C_N * factor * sinb ** 3 + C_T * cosb ** 3) if folds else np.full(n, C_T)
        vel = omega * np.c_[-X[:, 1], X[:, 0], np.zeros(n)]
        sp = np.linalg.norm(vel, axis=1)
        f += -0.5 * rho * ceff[:, None] * rc[:, None] * sp[:, None] * vel
        Fb = np.cumsum((f * ds)[::-1], axis=0)[::-1]
        tip = X[-1]
        Fb += np.array([hem_m * omega ** 2 * tip[0], hem_m * omega ** 2 * tip[1], -hem_m * g])
        t = Fb / np.linalg.norm(Fb, axis=1, keepdims=True)
        if slope_max < 1.0:
            rhat = X[:, :2] / np.maximum(Rh, 1e-12)[:, None]
            that = np.c_[-rhat[:, 1], rhat[:, 0]]
            tr = (t[:, :2] * rhat).sum(1)
            tt = (t[:, :2] * that).sum(1)
            over = tr > slope_max
            tr2 = np.where(over, slope_max, tr)
            rest = np.sqrt(np.maximum(0.0, 1 - tr2 ** 2))
            tt2 = np.where(over, np.clip(tt, -rest, rest), tt)
            tz2 = np.where(over, -np.sqrt(np.maximum(0.0, 1 - tr2 ** 2 - tt2 ** 2)), t[:, 2])
            t = np.where(over[:, None], np.c_[tr2[:, None] * rhat + tt2[:, None] * that, tz2], t)
        Xn = np.array([r_w, 0.0, 0.0]) + np.cumsum(t * ds, axis=0) - 0.5 * t * ds
        res = np.abs(Xn - X).max()
        X = (1 - relax) * X + relax * Xn
        if res < 1e-9:
            break
    hem = X[-1] + 0.5 * ds * t[-1]
    R = math.hypot(hem[0], hem[1])
    return dict(hem_radius_m=round(R, 3), drop_m=round(float(-hem[2]), 3),
                flare_deg=round(math.degrees(math.atan2(R - r_w, -hem[2])), 1),
                lag_deg=round(-math.degrees(math.atan2(hem[1], hem[0])), 1), converged=bool(res < 1e-6))


def wind_up_index(m, rho, cut=2 * math.pi, hem_kg_m=0.0, factor=1.0, L=SKIRT_LENGTH_M, r_w=None):
    """How hard the air's drag on a turning skirt's folds pulls against its centrifugal pull: the folds' circumferential
    drag over the radial pull, summed down the skirt with the skirt open to the width of its own cut (rc(s), flat for a
    full circle, the smooth cone for a part circle) and its folds' faces square to the flow, an upper bound:
    W = (rho C_N / 2) int rc^3 ds / (m int rc^2 ds + hem_kg_m rc(L)^2). Both loads grow as the square of the rate of
    turn and neither holds g, so W is the same at any rate and on either world. Above 1 the folds' drag outweighs the
    pull and the skirt winds round the body; a fold-free skirt meets skin friction only."""
    r_w = waist_radius() if r_w is None else r_w
    s = np.linspace(0.0, L, 2001)
    rc = r_w + cut * s / (2 * math.pi)
    trap = lambda f: float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(s)))
    drag = 0.5 * rho * C_N * factor * trap(rc ** 3)
    return drag / (m * trap(rc ** 2) + hem_kg_m * rc[-1] ** 2)


def hem_for_index(m, rho, target=1.0, cut=2 * math.pi, factor=1.0, L=SKIRT_LENGTH_M, r_w=None):
    """The weighted hem (kg per metre of hem) that brings the wind-up index down to target."""
    r_w = waist_radius() if r_w is None else r_w
    s = np.linspace(0.0, L, 2001)
    rc = r_w + cut * s / (2 * math.pi)
    trap = lambda f: float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(s)))
    need = (0.5 * rho * C_N * factor * trap(rc ** 3) / target - m * trap(rc ** 2)) / rc[-1] ** 2
    return max(0.0, need)


# --- the screen -----------------------------------------------------------------------------------------------------

def fabric_rows(W):
    rows = []
    for name, mat in FABRICS.items():
        m, f = mat.areal_density, mat.normal_factor
        row = dict(fabric=name, areal_density_kg_m2=m, bending_rigidity_uN_m=round(mat.bending_rigidity * 1e6, 3),
                   normal_factor=f)
        for wn, w in W.items():
            g, rho = w['g'], w['rho']
            vt = face_on_fall(m, g, rho, f)
            row[wn] = dict(
                bending_length_cm=round(100 * bending_length(mat, g), 2),
                face_on_fall_m_s=round(vt, 3),
                fall_0_5_m_s=round(fall_time(0.5, vt, g), 2),
                fall_1_5_m_s=round(fall_time(1.5, vt, g), 2),
                drag_length_cm=round(100 * 2 * m / (rho * C_N * f), 2),
                added_mass_over_fabric=dict(fold=round(added_mass(rho, FOLD_SCALE_M, f) / m, 1),
                                            garment=round(added_mass(rho, GARMENT_SCALE_M, f) / m, 1)),
                swing=swing(m, g, rho, f),
                trail_deg={f'{U:g}_m_s': round(trail_angle(m, g, rho, U, f), 1) for U in TRAVEL_M_S})
        row['earth_bending_rigidity_for_the_moon_drape_uN_m'] = round(
            mat.bending_rigidity * W['earth']['g'] / W['moon']['g'] * 1e6, 2)
        # falling and trailing turn on weight against the air, m g / rho: a lunar fabric falls and trails as one
        # lighter by (g_E / g_M)(rho_M / rho_E) = 6.9 does on Earth
        row['falls_and_trails_as_on_earth_kg_m2'] = round(
            m * W['moon']['g'] / W['moon']['rho'] * W['earth']['rho'] / W['earth']['g'], 4)
        row['bending_status'] = BENDING_STATUS.get(name, 'unstated')
        rows.append(row)
    return rows


def turn_rows(W):
    """Circle skirts in steady turns: the shape without air (weight against centrifugal pull), the shape with skin
    friction alone (a fold-free skirt, the lower bound of the air's drag), and the wind-up index of the folds' drag."""
    rows = []
    for name, mat in FABRICS.items():
        m, f = mat.areal_density, mat.normal_factor
        for cut_name, cut in (('full circle', 2 * math.pi), ('half circle', math.pi)):
            r = dict(fabric=name, cut=cut_name,
                     wind_up_index=dict(moon=round(wind_up_index(m, W['moon']['rho'], cut, factor=f), 2),
                                        earth=round(wind_up_index(m, W['earth']['rho'], cut, factor=f), 2)))
            for wn, w in W.items():
                g, rho = w['g'], w['rho']
                r[wn] = {f'{rpm:g}_rpm': dict(no_air=meridian(m, g, 0.0, rpm * math.pi / 30, cut=cut),
                                              skin_friction=meridian(m, g, rho, rpm * math.pi / 30, cut=cut, factor=f,
                                                                     folds=False))
                         for rpm in TURN_RPM}
            rows.append(r)
    return rows


def hem_rows(W):
    rows = []
    for name, mat in FABRICS.items():
        m, f = mat.areal_density, mat.normal_factor
        rho = W['moon']['rho']
        rows.append(dict(fabric=name, cut='full circle',
                         wind_up_index_moon={f'{h * 1000:g}_g_m': round(wind_up_index(m, rho, hem_kg_m=h, factor=f), 2)
                                             for h in HEM_KG_M},
                         hem_kg_m_for_index={f'{x:g}': round(hem_for_index(m, rho, x, factor=f), 3)
                                             for x in (1.0, 0.5)}))
    return rows


def results():
    W = worlds()
    return dict(
        schema=SCHEMA,
        producer=dict(domain='humanities', files={f: digest(f) for f in FILES},
                      constants=constants_used([ROOT / f for f in FILES]),
                      inputs={p: digest(p) for p in INPUTS}),
        evidence='Closed-form and one-dimensional estimates for the fabric library in materials.py (values and their '
                 'standing in sources.json) at the shared gravity and the design run\'s sea-level air, against Earth. '
                 'Falls use the solver\'s normal-pressure law; swings add a strip\'s potential-flow added mass '
                 '(a Morison-type estimate); turns use a meridian string in still air without air and with skin '
                 'friction alone (a fold-free skirt), and an upper-bound index of the folds\' drag against the '
                 'centrifugal pull. A screen for choosing fabrics and cuts; the cloth solver\'s trial runs test it.',
        worlds={k: dict(gravity_m_s2=v['g'], air_kg_m3=round(v['rho'], 4)) for k, v in W.items()},
        scalings=dict(gravity_clock=round(math.sqrt(G_E / G_M), 3),
                      bending_length_moon_over_earth=round((G_E / G_M) ** (1 / 3), 3),
                      forced_load_over_weight_moon_over_earth=round(G_E / G_M * W['moon']['rho'] / RHO_EARTH, 2)),
        coefficients=dict(normal=C_N, skin_friction=C_T),
        settings=dict(skirt_length_m=SKIRT_LENGTH_M, waist_radius_m=round(waist_radius(), 4), panel=PANEL,
                      fold_scale_m=FOLD_SCALE_M, garment_scale_m=GARMENT_SCALE_M, travel_m_s=TRAVEL_M_S,
                      turn_rpm=TURN_RPM, hem_kg_m=HEM_KG_M),
        fabrics=fabric_rows(W),
        turns=turn_rows(W),
        weighted_hems=hem_rows(W))


def main() -> int:
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
