"""Raindrop breakup at lunar gravity: the size scale that carries the rain limits of CM1's Morrison scheme to the Moon.

    python -m atmosphere.rain.breakup        # about a minute on one thread
    # -> atmosphere/rain/results/rain_breakup.json

CM1's Morrison scheme gives rain an exponential spectrum and limits its mean diameter D_m = 1/lambda twice. Rain's
self-collection is offset by collisional breakup above a cutoff D_b of 300 um, by the factor
2 - exp(2300 m^-1 (D_m - D_b)) (Verlinde and Cotton 1993), so that D_m settles at D_b + ln 2 / 2300 m^-1 = 0.60 mm;
and D_m is held below 2.8 mm. The rain's mean fall speeds are held below 9.1 m/s in the scheme's reference air. All
three are Earth's (shared/constants.json, model_closures).

Two drops that collide coalesce or come apart as the collision's kinetic energy against surface tension decides. The
outcome follows the collision's Weber number, We = CKE / S_c, its kinetic energy over the surface energy of the
coalesced drop: the pair coalesces with efficiency exp(-1.15 We) (Straub et al. 2010), the liquid's viscosity entering
only through an Ohnesorge number of a few thousandths. At lunar gravity the same drops fall more slowly and meet more
gently. A pair of drops on another world collides as an Earth pair of the same size ratio does when it is s times that
pair's size, s set by the Weber number through each world's terminal speeds, which the single-drop solver gives
(drops.py). Weighted by the breakup each Earth pair carries in Morrison's spectra at the cutoff and at the
equilibrium, the pairs' s gives the size scale of the breakup law on that world.

At equal air density s is the ratio of the worlds' capillary lengths sqrt(sigma / (rho_w g)), to within 0.2% at the
equilibrium and 4% at the cutoff: a drop's shape follows that length, its inertial fall the speed
(g sigma / rho_w)^(1/4) (rho_w / rho_air)^(1/2) it implies, and the Weber numbers of its collisions both. Denser air
slows the drops and softens their collisions, so s grows with the air's density as (rho_air / rho_air,Earth)^alpha:
at the equilibrium alpha is 0.69-0.71 for gravities from the Moon's to Earth's at 1.2 times Earth's density. The cap
on D_m bounds thin rain sorted by its fall, where breakup is too slow to act; it follows the largest drop that holds
together, the stability screen, which scales with the capillary length alone. The cap on the mean fall speed follows
that drop's terminal speed, which under the same similarity scales as (g / 9.81)^(1/4). Low and List's (1982) fit of
the coalescence efficiency, read with its dimensional constants at lunar sizes and speeds, makes the Moon's drops
break at smaller sizes: its exponent carries an absolute length of about 5 mm, which on Earth's terminal speeds stands
in for the Weber number, and lunar pairs lie beyond the 0.4-4.6 mm it was measured over. Its weights check the Weber
law's, and its literal reading is given beside them.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from shared.constants import (CM1_GRAVITY, CM1_RAIN_BREAKUP_CUTOFF, CM1_RAIN_BREAKUP_RATE, CM1_RAIN_FALL_CAP,
                              CM1_RAIN_FALL_COEFFICIENT, CM1_RAIN_FALL_EXPONENT, CM1_RAIN_MEAN_DIAMETER_CAP,
                              MOON_SURFACE_GRAVITY, STANDARD_GRAVITY)
from shared.provenance import constants_used
from atmosphere.radiative_convective.thermodynamics import DryAir, earthlike_air
from atmosphere.rain import drops, product as rain_product, targets as climate_targets

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.rain-breakup/1'
OUT = HERE / 'results' / 'rain_breakup.json'
CODE = ('atmosphere/rain/breakup.py', 'atmosphere/rain/drops.py', 'atmosphere/rain/product.py',
        'atmosphere/rain/targets.py', 'atmosphere/radiative_convective/thermodynamics.py')
# Straub, Beheng, Seifert, Schlottke and Weigand (2010, part II): colliding raindrops coalesce with efficiency
# exp(-1.15 We), We = CKE / S_c.
STRAUB_DECAY = 1.15
# Low and List (1982a): E = a (1 + p)^-2 exp(-b sigma E_T^2 / S_c) for E_T = CKE + S_T - S_c below 5 uJ, 0 above,
# fitted to pairs of 0.395-4.6 mm and 0.395-1.8 mm.
LOW_LIST_A, LOW_LIST_B, LOW_LIST_TOP_J = 0.778, 2.61e6, 5.0e-6
LOW_LIST_LARGE_M, LOW_LIST_SMALL_M = (0.395e-3, 4.6e-3), (0.395e-3, 1.8e-3)
EQUILIBRIUM_M = CM1_RAIN_BREAKUP_CUTOFF + math.log(2.0) / CM1_RAIN_BREAKUP_RATE   # where breakup balances, 0.60 mm
SPECTRA = (('cutoff', CM1_RAIN_BREAKUP_CUTOFF), ('equilibrium', EQUILIBRIUM_M))
TABLE_POINTS, TABLE_SMALLEST_M = 320, 0.02e-3     # terminal speeds on a log grid up to just below the screen
SCREEN_MARGIN = 0.998
PAIR_POINTS, PAIR_SMALLEST_M = 200, 0.1e-3        # Earth's drops in the spectra, from drizzle to the screen
LINE_SAMPLES, BISECTIONS = 160, 60
ATMOSPHERE_PA = rain_product.ATMOSPHERE_PA
CO2_PPM = rain_product.CO2_PPM
# The general check: gravities (m/s2) and surface pressures (atm) of the Moon's composition of air at the resort's 2 m
# temperature and humidity; 3.7 m/s2 is a gravity between the Moon's and Earth's, near Mars's.
GENERAL = ((STANDARD_GRAVITY, 1.2), (STANDARD_GRAVITY, 1.5), (MOON_SURFACE_GRAVITY, 1.0), (MOON_SURFACE_GRAVITY, 1.2),
           (MOON_SURFACE_GRAVITY, 1.5), (3.7, 1.0), (3.7, 1.2))
ADOPTED = ('resort', 'surface')                    # the target and height the lunar values are given at

EVIDENCE = ('A runnable calculation: the single-drop solver\'s terminal speeds, which reproduce Earth\'s measured fall '
            'speeds to 0.20 m/s, in each world\'s air, with the collision outcome taken to follow the Weber number '
            'alone, as Straub et al. (2010) found for raindrop pairs, the liquid\'s viscosity and the gas film '
            'between the drops left out. The weights are Morrison\'s exponential spectra on Earth. The lunar values '
            'are Earth\'s breakup carried by that similarity; no laboratory has measured lunar collisions, and no '
            'CM1 run has tested the values.')
READING_RULE = ('Sizes in mm. scale gives, for each target and height (surface: 2 m; cloud_base: the parcel cloud base of '
                'each world\'s column), the size factor s that carries Morrison\'s spectrum at the cutoff (0.30 mm) '
                'and at the equilibrium (0.60 mm) from Earth to the Moon: the breakup-weighted geometric mean over '
                'drop pairs of the ratio of the sizes at which each world\'s pairs of one size ratio reach one Weber '
                'number, with its 10, 50 and 90% points and the share of the breakup in Earth pairs the Moon cannot '
                'match below its stability screen (taken at the Moon\'s largest Weber number for that ratio, a lower '
                'bound on s). low_list_weights is the same mean with Low and List\'s energetic breakup as the weights; '
                'low_list_literal reads their fit with its dimensional constants on both worlds for the pairs within '
                'its measured sizes, with the share of the breakup it covers. capillary is the ratio of capillary '
                'lengths. alpha = ln(s / capillary) / ln(density ratio). general repeats the scale at other gravities '
                'and pressures. lunar_limits gives the values adopted at the resort\'s 2 m air as the rain-limits patch '
                'of climate/crm realizes them at the design pressure, and its namelist settings.')
ASSUMPTIONS = [
    'Drops collide at their terminal speeds in still air, with the geometric kernel pi/4 (D_L + D_s)^2 |V_L - V_s| '
    '(collision efficiency 1), between 0.1 mm and the stability screen.',
    'Each pair\'s outcome is set by its Weber number and size ratio; the impact parameter averages alike on both '
    'worlds. The Ohnesorge number of a 1 mm drop, 0.004 on Earth and 0.0025 at the lunar size, and the 20% denser '
    'gas film are left out.',
    'Morrison\'s exponential spectra stand for the rain whose breakup the scheme\'s law was fitted to; the law\'s two '
    'lengths, the cutoff and its e-folding length, scale together by the equilibrium\'s s.',
    'The Moon\'s air is the design air at 1.2 atm and Earth\'s is at 1 atm, at each target\'s 2 m temperature and '
    'humidity (rain_drops.json\'s airs).',
]


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Speeds:
    """Terminal speeds in still air at one height of a column on a log grid of sizes up to just below its stability
    screen, with the air's density and water's surface tension there."""

    def __init__(self, column, altitude_m=0.0):
        state = column.state(altitude_m)
        self.gravity = column.gravity_m_s2
        self.density = state['density_kg_m3']
        self.sigma = drops.surface_tension(state['temperature_k'])
        self.screen = drops.stable_diameter(column, altitude_m)
        self.diameter = np.geomspace(TABLE_SMALLEST_M, SCREEN_MARGIN * self.screen, TABLE_POINTS)
        self.speed = np.array([drops.terminal_speed(d, column, altitude_m) for d in self.diameter])
        self._log = (np.log(self.diameter), np.log(self.speed))

    def __call__(self, diameter_m):
        return np.exp(np.interp(np.log(diameter_m), *self._log))

    @property
    def capillary_length(self):
        return math.sqrt(self.sigma / (self.gravity * (drops.WATER_DENSITY - self.density)))


def collision(speeds, large, small):
    """Pairs of drops falling at their terminal speeds: the collision's kinetic energy CKE, the coalesced drop's surface
    energy S_c, the two drops' surface energy S_T (J) and the Weber number CKE / S_c."""
    reduced = math.pi / 12 * drops.WATER_DENSITY * large**3 * small**3 / (large**3 + small**3)
    cke = reduced * (speeds(large) - speeds(small))**2
    s_c = math.pi * speeds.sigma * (large**3 + small**3)**(2 / 3)
    s_t = math.pi * speeds.sigma * (large**2 + small**2)
    return cke, s_c, s_t, cke / s_c


def straub(weber):
    """Coalescence efficiency by Straub et al. (2010): exp(-1.15 We)."""
    return np.exp(-STRAUB_DECAY * weber)


def low_list(speeds, large, small):
    """Low and List's (1982a) coalescence efficiency and its limit at no collision energy, a (1 + p)^-2."""
    cke, s_c, s_t, _ = collision(speeds, large, small)
    energy = cke + s_t - s_c
    ceiling = LOW_LIST_A * (1 + small / large)**-2
    return np.where(energy < LOW_LIST_TOP_J, ceiling * np.exp(-LOW_LIST_B * speeds.sigma * energy**2 / s_c), 0.0), ceiling


def rising_size(speeds, ratio, weber):
    """For pairs of size ratio p = D_s / D_L, the smallest large drop on a world whose pair reaches a Weber number:
    along D_s = p D_L the Weber number rises with size to a largest value below the stability screen and may fall
    beyond it, so the size is sought on the rising branch. Returns the size (nan where the line never reaches the Weber
    number), the size of the line's largest Weber number and that number."""
    lo = np.log(np.maximum(TABLE_SMALLEST_M / ratio, TABLE_SMALLEST_M) * 1.0001)
    hi = np.full_like(lo, math.log(SCREEN_MARGIN * speeds.screen * 0.9999))
    line = np.exp(lo[:, None] + (hi - lo)[:, None] * np.linspace(0.0, 1.0, LINE_SAMPLES))
    w = collision(speeds, line, ratio[:, None] * line)[3]
    rows, top = np.arange(len(ratio)), np.argmax(w, axis=1)

    def at(log_size):
        return collision(speeds, np.exp(log_size), ratio * np.exp(log_size))[3]

    # the line's largest Weber number, sampled finely between the samples beside the largest sample (the speeds'
    # interpolation leaves the top of the line flat and slightly uneven, so a search for one maximum can stop short)
    a = np.log(line[rows, np.maximum(top - 1, 0)])
    b = np.log(line[rows, np.minimum(top + 1, LINE_SAMPLES - 1)])
    fine = a[:, None] + (b - a)[:, None] * np.linspace(0.0, 1.0, LINE_SAMPLES)
    w_fine = collision(speeds, np.exp(fine), ratio[:, None] * np.exp(fine))[3]
    top_log = fine[rows, np.argmax(w_fine, axis=1)]
    peak = np.maximum(w_fine.max(axis=1), w[rows, top])
    a, b = lo.copy(), top_log
    for _ in range(BISECTIONS):
        middle = 0.5 * (a + b)
        lower = at(middle) < weber
        a, b = np.where(lower, middle, a), np.where(lower, b, middle)
    size = np.where(weber <= w[:, 0], np.exp(lo), np.exp(0.5 * (a + b)))
    return np.where(weber <= peak * (1.0 + 1e-4), size, np.nan), np.exp(top_log), peak


def pair_scale(world, earth, large, small):
    """Size factor s for Earth pairs: the smallest pair of the same size ratio on the world reaching the Earth pair's
    Weber number over the smallest such pair on Earth. Pairs the world cannot match take the size of its largest Weber
    number at that ratio, and are flagged."""
    ratio = small / large
    weber = collision(earth, large, small)[3]
    on_earth, earth_top, _ = rising_size(earth, ratio, weber)
    on_earth = np.where(np.isfinite(on_earth), on_earth, earth_top)    # at the sampled peak of Earth's own line
    on_world, top_size, _ = rising_size(world, ratio, weber)
    reached = np.isfinite(on_world)
    return np.where(reached, on_world, top_size) / on_earth, reached


def breakup_weights(earth, mean_diameter_m):
    """Earth's drop pairs (large, small) in Morrison's exponential spectrum of a mean diameter, from 0.1 mm to the
    stability screen, each weighted by the breakup it carries: the geometric kernel times both drops' numbers times the
    share that does not coalesce, by Straub et al.'s law and by Low and List's energetic breakup (their limit at no
    energy less their efficiency)."""
    d = np.geomspace(PAIR_SMALLEST_M, SCREEN_MARGIN * 0.999 * earth.screen, PAIR_POINTS)
    number = np.gradient(np.log(d)) * d * np.exp(-d / mean_diameter_m)
    i, j = np.triu_indices(PAIR_POINTS, k=1)
    large, small = d[j], d[i]
    kernel = math.pi / 4 * (large + small)**2 * np.abs(earth(large) - earth(small)) * number[j] * number[i]
    weber = collision(earth, large, small)[3]
    efficiency, ceiling = low_list(earth, large, small)
    return large, small, dict(straub=kernel * (1.0 - straub(weber)), low_list=kernel * (ceiling - efficiency))


def scale(world, earth, mean_diameter_m, law='straub'):
    """The breakup-weighted size factor carrying an Earth spectrum of a mean diameter to the world: geometric mean, 10,
    50 and 90% points, the share of the weight the world cannot match, and the breaking pairs' mean sizes."""
    large, small, weights = breakup_weights(earth, mean_diameter_m)
    w = weights[law]
    keep = w > 1e-9 * w.max()
    large, small, w = large[keep], small[keep], w[keep]
    s, reached = pair_scale(world, earth, large, small)
    log_s = np.log(s)
    order = np.argsort(log_s)
    cumulative = np.cumsum(w[order]) / w.sum()

    def point(f):
        return float(np.exp(log_s[order][min(int(np.searchsorted(cumulative, f)), len(order) - 1)]))

    return dict(scale=float(np.exp(np.sum(w * log_s) / w.sum())), p10=point(0.1), p50=point(0.5), p90=point(0.9),
                unreached_share=float(w[~reached].sum() / w.sum()),
                breaking_large_mm=float(np.sum(w * large) / w.sum() * 1e3),
                breaking_small_mm=float(np.sum(w * small) / w.sum() * 1e3))


def low_list_literal(world, earth, mean_diameter_m):
    """Low and List's fit read with its dimensional constants on both worlds: for the Earth pairs within its measured
    sizes that sometimes coalesce, the factor s at which the world's pair of the same size ratio has the Earth pair's
    efficiency, weighted by their energetic breakup in the spectrum; and the share of that breakup the check covers."""
    large, small, weights = breakup_weights(earth, mean_diameter_m)
    w = weights['low_list']
    target = low_list(earth, large, small)[0]
    inside = ((large >= LOW_LIST_LARGE_M[0]) & (large <= LOW_LIST_LARGE_M[1]) & (small >= LOW_LIST_SMALL_M[0])
              & (small <= LOW_LIST_SMALL_M[1]))
    use = inside & (target > 0.0) & (w > 0.0)
    large, small, target, weight = large[use], small[use], target[use], w[use]
    a = np.full(large.shape, math.log(0.5))
    b = np.log(SCREEN_MARGIN * 0.9999 * world.screen / large)
    above = low_list(world, np.exp(a) * large, np.exp(a) * small)[0] >= target
    below = low_list(world, np.exp(b) * large, np.exp(b) * small)[0] <= target
    bracketed = above & below
    for _ in range(BISECTIONS):                         # the efficiency falls as the pair grows at a fixed ratio
        middle = 0.5 * (a + b)
        high = low_list(world, np.exp(middle) * large, np.exp(middle) * small)[0] > target
        a, b = np.where(high, middle, a), np.where(high, b, middle)
    log_s = 0.5 * (a + b)
    covered = weight[bracketed].sum() / w.sum()
    return dict(scale=float(np.exp(np.sum(weight[bracketed] * log_s[bracketed]) / weight[bracketed].sum())),
                breakup_share_covered=float(covered))


def comparison(world, earth):
    """The scale at the cutoff and at the equilibrium by both laws' weights, the capillary ratio and alpha."""
    density_ratio = world.density / earth.density
    capillary = world.capillary_length / earth.capillary_length
    out = dict(density_kg_m3=dict(world=world.density, earth=earth.density), density_ratio=density_ratio,
               capillary=capillary, stability_screen_mm=dict(world=world.screen * 1e3, earth=earth.screen * 1e3))
    for name, dm in SPECTRA:
        row = scale(world, earth, dm)
        row['low_list_weights'] = scale(world, earth, dm, law='low_list')['scale']
        row['low_list_literal'] = low_list_literal(world, earth, dm)
        row['over_capillary'] = row['scale'] / capillary
        row['alpha'] = (math.log(row['over_capillary']) / math.log(density_ratio)
                        if abs(density_ratio - 1.0) > 1e-3 else None)
        out[name] = row
    return out


def general_column(gravity, atmospheres, air):
    pressure = atmospheres * ATMOSPHERE_PA
    return drops.Column(gravity_m_s2=gravity, surface_pressure_pa=pressure, surface_temperature_k=air['air_2m_k'],
                        surface_relative_humidity=air['relative_humidity'], air=earthlike_air(pressure, CO2_PPM))


def mass_weighted_speed(speeds, mean_diameter_m):
    """The mass-weighted mean terminal speed of an exponential spectrum, up to the stability screen."""
    d = np.geomspace(1e-5, speeds.diameter[-1], 2000)
    weight = np.gradient(d) * d**3 * np.exp(-d / mean_diameter_m)
    return float(np.sum(weight * speeds(d)) / np.sum(weight))


def fall_speed_cap(cols, lunar_equilibrium_m):
    """The largest drops' terminal speeds against the similarity's (g / 9.81)^(1/4), and Morrison's mean fall speed of
    rain at the lunar equilibrium against the caps."""
    moon, earth_g = Speeds(cols['moon']), Speeds(cols['moon_air_earth_gravity'])
    largest = dict(moon=float(moon.speed[-1]), moon_air_earth_gravity=float(earth_g.speed[-1]),
                   earth=float(Speeds(cols['earth']).speed[-1]))
    density = moon.density
    factor = rain_product.cm1_density_factor(density)
    gravity_factor = rain_product.cm1_gravity_factor(MOON_SURFACE_GRAVITY)
    terluna_cap = CM1_RAIN_FALL_CAP * gravity_factor * factor
    similar_cap = CM1_RAIN_FALL_CAP * (MOON_SURFACE_GRAVITY / CM1_GRAVITY)**0.25 * factor
    morrison = (CM1_RAIN_FALL_COEFFICIENT * gravity_factor * factor * math.gamma(4 + CM1_RAIN_FALL_EXPONENT) / 6
                * lunar_equilibrium_m**CM1_RAIN_FALL_EXPONENT)
    return dict(largest_drop_speed_m_s=largest,
                gravity_ratio_in_the_same_air=largest['moon'] / largest['moon_air_earth_gravity'],
                similarity_quarter_power=(MOON_SURFACE_GRAVITY / CM1_GRAVITY)**0.25,
                at_2m_air=dict(density_kg_m3=density, terluna_cap_m_s=terluna_cap, similarity_cap_m_s=similar_cap,
                               lunar_equilibrium_mm=lunar_equilibrium_m * 1e3,
                               morrison_mass_weighted_speed_m_s=morrison,
                               solver_mass_weighted_speed_m_s=mass_weighted_speed(moon, lunar_equilibrium_m)))


def lunar_limits(row, cols):
    """The values adopted at the resort's 2 m air: the scale of the collisional breakup at the equilibrium, its
    exponent on the air's density, the cap on the mean diameter by the capillary length and the fall-speed cap, as the
    rain-limits patch realizes them from its namelist at CM1's gravity of reference and the design pressure."""
    alpha = row['equilibrium']['alpha']
    rounded = round(alpha, 3)
    capillary = math.sqrt(CM1_GRAVITY / MOON_SURFACE_GRAVITY)          # the patch's (9.81 / g)^(1/2)
    pressure_ratio = rain_product.MOON_PRESSURE_PA / ATMOSPHERE_PA      # CM1's density ratio at equal temperature
    collisional = capillary * pressure_ratio**rounded
    cutoff = CM1_RAIN_BREAKUP_CUTOFF * collisional
    rate = CM1_RAIN_BREAKUP_RATE / collisional
    equilibrium = cutoff + math.log(2.0) / rate
    return dict(
        scale_capillary=capillary, scale_collisional=collisional, alpha=alpha,
        derived_scale_collisional=row['equilibrium']['scale'],
        earth=dict(breakup_cutoff_mm=CM1_RAIN_BREAKUP_CUTOFF * 1e3, breakup_rate_per_m=CM1_RAIN_BREAKUP_RATE,
                   equilibrium_mm=EQUILIBRIUM_M * 1e3, mean_diameter_cap_mm=CM1_RAIN_MEAN_DIAMETER_CAP * 1e3,
                   fall_speed_cap_m_s=CM1_RAIN_FALL_CAP,
                   terluna_fall_speed_cap_m_s=CM1_RAIN_FALL_CAP * rain_product.cm1_gravity_factor(MOON_SURFACE_GRAVITY)),
        moon=dict(breakup_cutoff_mm=cutoff * 1e3, breakup_rate_per_m=rate, equilibrium_mm=equilibrium * 1e3,
                  mean_diameter_cap_mm=CM1_RAIN_MEAN_DIAMETER_CAP * capillary * 1e3,
                  fall_speed_cap_m_s=CM1_RAIN_FALL_CAP * (MOON_SURFACE_GRAVITY / CM1_GRAVITY)**0.25),
        namelist=dict(brkup_gscale=1, brkup_rhoexp=rounded,
                      brkup_rhor='the case\'s surface pressure over 1013.25 hPa, CM1\'s ratio of its air\'s density to '
                                 'Earth\'s at the same temperature (1.2 at the design pressure)',
                      brkup_db=CM1_RAIN_BREAKUP_CUTOFF),
        fall_speed=fall_speed_cap(cols, equilibrium))


def build():
    data, binding = climate_targets.load()
    air = rain_product.target_air(climate_targets.read(data))
    product = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, assumptions=ASSUMPTIONS,
                   earth_limits=dict(breakup_cutoff_mm=CM1_RAIN_BREAKUP_CUTOFF * 1e3,
                                     breakup_rate_per_m=CM1_RAIN_BREAKUP_RATE, equilibrium_mm=EQUILIBRIUM_M * 1e3,
                                     mean_diameter_cap_mm=CM1_RAIN_MEAN_DIAMETER_CAP * 1e3,
                                     fall_speed_cap_m_s=CM1_RAIN_FALL_CAP),
                   scale={})
    for key in rain_product.TARGETS:
        cols = rain_product.columns(air[key])
        moon, earth = cols['moon'], cols['earth']
        product['scale'][key] = {
            'surface': comparison(Speeds(moon), Speeds(earth)),
            'cloud_base': dict(comparison(Speeds(moon, moon.cloud_base_m), Speeds(earth, earth.cloud_base_m)),
                               heights_m=dict(moon=moon.cloud_base_m, earth=earth.cloud_base_m))}
    resort = air['resort']
    earth = Speeds(rain_product.columns(resort)['earth'])
    product['general'] = [dict(gravity_m_s2=g, pressure_atm=p, **comparison(Speeds(general_column(g, p, resort)), earth))
                          for g, p in GENERAL]
    target, height = ADOPTED
    product['lunar_limits'] = lunar_limits(product['scale'][target][height], rain_product.columns(resort))
    product['producer'] = dict(domain='atmosphere', files={p: sha256(ROOT / p) for p in CODE}, **binding,
                               constants=constants_used([ROOT / p for p in CODE]))
    return product


def main(argv=None) -> int:
    product = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    for key, rows in product['scale'].items():
        for height, row in rows.items():
            print(f"{key:9s} {height:10s} capillary {row['capillary']:.4f}  cutoff {row['cutoff']['scale']:.3f}  "
                  f"equilibrium {row['equilibrium']['scale']:.3f} (alpha {row['equilibrium']['alpha']:.3f}, "
                  f"unreached {row['equilibrium']['unreached_share']:.2f})")
    for row in product['general']:
        print(f"g {row['gravity_m_s2']:6.3f} {row['pressure_atm']:.1f} atm: over capillary {row['cutoff']['over_capillary']:.3f}"
              f" / {row['equilibrium']['over_capillary']:.3f}")
    moon = product['lunar_limits']['moon']
    print(f"lunar limits: cutoff {moon['breakup_cutoff_mm']:.3f} mm, rate {moon['breakup_rate_per_m']:.0f}/m, "
          f"equilibrium {moon['equilibrium_mm']:.3f} mm, cap {moon['mean_diameter_cap_mm']:.2f} mm, "
          f"fall-speed cap {moon['fall_speed_cap_m_s']:.2f} m/s")
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
