#!/usr/bin/env python3
"""Snow in the Open Moon's sky winter: how snowflakes fall, grow and melt at a sixth of Earth's gravity, beside the
raindrop solver (atmosphere/rain).

    python -m atmosphere.snow.flakes            # seconds: results/snow.json

A snowflake here is an aggregate of ice crystals. Its mass follows Brown and Francis (1995), m = 0.0185 D^1.9 kg for a
maximum dimension D in metres; its projected area is 0.35 of the circle across D, the area ratio large aggregates of
dendrites and plates approach (0.33-0.41; Leinonen and co-authors 2026); and it falls by Heymsfield and Westbrook's
(2010) law: the modified Best number X* = 8 m g rho / (pi Ar^0.5 eta^2), the Reynolds number
Re = (d0^2/4) ((1 + 4 X*^0.5 / (d0^2 C0^0.5))^0.5 - 1)^2 with d0 = 8.0 and C0 = 0.35, and the speed eta Re / (rho D).
Gravity enters through the Best number alone. The flake grows from the vapour by diffusion to a capacitance of 0.25 D
(Westbrook, Hogan and Illingworth 2008), ventilated as Hall and Pruppacher (1976) give it for ice, against the
saturation pressures of Murphy and Koop (2005); it sweeps up the cloud's small crystals with Lin, Farley and Orville's
(1983) sticking efficiency, exp(0.09 (T - 273.15)); and below the freezing level it melts by the heat it takes from the
air less what its wet surface's evaporation carries away (Pruppacher and Klett 2010, chapter 16), at the same
capacitance and ventilation. The air's viscosity, conductivity and vapour diffusivity are the raindrop solver's
(drops.transport) for each world's air.

A winter column is a straight temperature profile through the freezing level, its pressure hydrostatic. Its cloud runs
from a top temperature down to the freezing level, at liquid saturation (a mixed-phase cloud, where snow grows fastest)
and holding small crystals of 0.2 mm at a set mixing ratio; below the freezing level the air holds a set relative
humidity. A flake starts at the cloud top as a 0.2 mm crystal, falls at its terminal speed, grows to the freezing
level, then melts, keeping its size and mass while its ice turns to water. The Moon's column is the resort's sky ladder
at 301 K (research/studies/climate_targets, read as atmosphere/rain/targets.py reads it): the freezing level, the
pressure there and the lapse from the ladder's levels about it. Earth's is a winter column with its freezing level at
850 hPa and the standard 6.5 K a km. Breakup of large flakes, riming and the melting flake's collapse are left out.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from shared.constants import GAS_CONSTANT, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from atmosphere.radiative_convective.thermodynamics import DryAir, MOLAR_MASS, TRIPLE_POINT_K, earthlike_air
from atmosphere.rain import drops, targets

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.snow/1'
OUT = HERE / 'results' / 'snow.json'
CODE = ('atmosphere/snow/flakes.py', 'atmosphere/rain/drops.py', 'atmosphere/rain/targets.py')
T0 = 273.15                                  # K, the melting point
R_VAPOUR = GAS_CONSTANT / MOLAR_MASS['H2O']
LATENT_FUSION = 3.337e5                      # J/kg at 0 °C
MASS_COEFFICIENT, MASS_EXPONENT = 0.0185, 1.9          # Brown and Francis (1995), SI
AREA_RATIO = 0.35                            # large aggregates (Leinonen and co-authors 2026: 0.33-0.41)
DELTA0, C0 = 8.0, 0.35                       # Heymsfield and Westbrook (2010)
CAPACITANCE = 0.25                           # times D, aggregates (Westbrook, Hogan and Illingworth 2008)
CRYSTAL_M = 2e-4                             # the cloud's small crystals and the flake's start, 0.2 mm
EARTH_PRESSURE_PA, EARTH_LAPSE_K_M = 101325.0, 6.5e-3
EARTH_FREEZING_PA = 85000.0                  # a winter freezing level near 1.5 km
LADDER_SPAN_KM = (20.0, 40.0)                # the resort's ladder levels the lapse is fitted over


def ice_saturation_pressure(t):
    """Pa over ice (Murphy and Koop 2005, their eq. 7), for 110 K and up."""
    return math.exp(9.550426 - 5723.265 / t + 3.53068 * math.log(t) - 0.00728332 * t)


def liquid_saturation_pressure(t):
    """Pa over liquid water, supercooled too (Murphy and Koop 2005, their eq. 10), 123-332 K."""
    return math.exp(54.842763 - 6763.22 / t - 4.210 * math.log(t) + 0.000367 * t
                    + math.tanh(0.0415 * (t - 218.8)) * (53.878 - 1331.22 / t - 9.44523 * math.log(t) + 0.014025 * t))


def flake_mass(d):
    return MASS_COEFFICIENT * d ** MASS_EXPONENT


def flake_size(m):
    return (m / MASS_COEFFICIENT) ** (1.0 / MASS_EXPONENT)


class Column:
    """A winter column: a straight temperature profile through the freezing level at height z_freeze (m), pressure
    hydrostatic from the pressure there, in a world's gravity and dry air."""

    def __init__(self, name, gravity, air, freezing_height_m, freezing_pressure_pa, lapse_k_m):
        self.name, self.gravity, self.air = name, gravity, air
        self.z_freeze, self.p_freeze, self.lapse = freezing_height_m, freezing_pressure_pa, lapse_k_m
        self.exponent = gravity * air.molar_mass / (GAS_CONSTANT * lapse_k_m)

    def temperature(self, z):
        return T0 - self.lapse * (z - self.z_freeze)

    def height_of(self, t):
        return self.z_freeze + (T0 - t) / self.lapse

    def state(self, z, relative_humidity=None):
        """The air at height z: temperature, pressure, density, viscosity, vapour diffusivity, conductivity and the
        vapour pressure, at liquid saturation when no relative humidity is given."""
        t = self.temperature(z)
        p = self.p_freeze * (t / T0) ** self.exponent
        e = liquid_saturation_pressure(t) * (1.0 if relative_humidity is None else relative_humidity)
        viscosity, diffusivity, conductivity, _ = drops.transport(t, p, e / p, self.air.fractions)
        density = (p - e) * self.air.molar_mass / (GAS_CONSTANT * t) + e / (R_VAPOUR * t)
        return dict(t=t, p=p, rho=density, eta=viscosity, dv=diffusivity, ka=conductivity, e=e)


def fall_speed(d, mass, gravity, s):
    """Terminal speed (m/s) and Reynolds number of a particle of maximum dimension d and mass mass, by Heymsfield and
    Westbrook (2010)."""
    x = 8.0 * mass * gravity * s['rho'] / (math.pi * math.sqrt(AREA_RATIO) * s['eta'] ** 2)
    re = DELTA0 ** 2 / 4.0 * (math.sqrt(1.0 + 4.0 * math.sqrt(x) / (DELTA0 ** 2 * math.sqrt(C0))) - 1.0) ** 2
    return s['eta'] * re / (s['rho'] * d), re


def ventilation(re, s):
    """Hall and Pruppacher (1976) for ice: X = Sc^(1/3) Re^(1/2)."""
    x = (s['eta'] / (s['rho'] * s['dv'])) ** (1.0 / 3.0) * math.sqrt(re)
    return 1.0 + 0.14 * x * x if x < 1.0 else 0.86 + 0.28 * x


def deposition_rate(d, re, s):
    """kg/s gained from the vapour by a flake of size d (negative when it sublimates)."""
    t = s['t']
    ls = drops.latent_heat(T0) + LATENT_FUSION                   # sublimation, held at its 0 °C value
    esi = ice_saturation_pressure(t)
    resistance = (ls / (R_VAPOUR * t) - 1.0) * ls / (s['ka'] * t) + R_VAPOUR * t / (s['dv'] * esi)
    return 4.0 * math.pi * CAPACITANCE * d * ventilation(re, s) * (s['e'] / esi - 1.0) / resistance


def melting_rate(d, re, s):
    """kg/s of ice melted on a flake of size d whose wet surface stands at 0 °C (zero when the air takes no heat)."""
    rho_surface = liquid_saturation_pressure(T0) / (R_VAPOUR * T0)
    rho_air = s['e'] / (R_VAPOUR * s['t'])
    heat = 4.0 * math.pi * CAPACITANCE * d * ventilation(re, s) * (
        s['ka'] * (s['t'] - T0) + drops.latent_heat(T0) * s['dv'] * (rho_air - rho_surface))
    return max(heat, 0.0) / LATENT_FUSION


def grow(column, top_c):
    """A 0.2 mm crystal from the cloud top (top_c, °C) to the freezing level, growing from the vapour alone: its size,
    mass and speed there, the hours taken, the cloud's depth and the air's mass in it (kg/m2)."""
    z_top = column.height_of(T0 + top_c)

    def rates(z, y):
        m = y[0]
        d = flake_size(m)
        s = column.state(z)
        v, re = fall_speed(d, m, column.gravity, s)
        return [-deposition_rate(d, re, s) / v, -1.0 / v]

    m0 = flake_mass(CRYSTAL_M)
    out = solve_ivp(rates, (z_top, column.z_freeze), [m0, 0.0], rtol=1e-8, atol=[m0 * 1e-9, 1e-3], max_step=200.0)
    m, hours = out.y[0, -1], out.y[1, -1] / 3600.0
    d = flake_size(m)
    air_mass = (column.state(column.z_freeze)['p'] - column.state(z_top)['p']) / column.gravity
    return dict(cloud_top_c=top_c, cloud_depth_m=float(z_top - column.z_freeze), cloud_air_kg_m2=float(air_mass),
                size_mm=1e3 * d, mass_mg=1e6 * m,
                speed_m_s=fall_speed(d, m, column.gravity, column.state(column.z_freeze))[0],
                hours=float(hours))


def crossing_hours(column, top_c, size_mm):
    """Hours a flake of a set size takes to fall from the cloud top (top_c, °C) to the freezing level."""
    d = size_mm * 1e-3
    zs = np.linspace(column.height_of(T0 + top_c), column.z_freeze, 400)
    v = np.array([fall_speed(d, flake_mass(d), column.gravity, column.state(z))[0] for z in zs])
    return float(np.trapezoid(1.0 / v, -zs) / 3600.0)


def melt(column, size_mm, relative_humidity):
    """A flake of a set size falling from the freezing level through air at a set relative humidity: the depths below
    it (m) where half and all of its ice has melted, and the hours to melt, or None if it reaches the ground."""
    d = size_mm * 1e-3
    m = flake_mass(d)

    def rates(z, y):
        s = column.state(z, relative_humidity)
        v, re = fall_speed(d, m, column.gravity, s)
        return [melting_rate(d, re, s) / v, -1.0 / v]

    half = lambda z, y: y[0] - 0.5 * m
    gone = lambda z, y: y[0] - 1e-3 * m
    gone.terminal, gone.direction = True, -1
    out = solve_ivp(rates, (column.z_freeze, 0.0), [m, 0.0], events=[half, gone], rtol=1e-8, atol=[m * 1e-9, 1e-3],
                    max_step=20.0)
    depth = lambda k: float(column.z_freeze - out.t_events[k][0]) if len(out.t_events[k]) else None
    hours = float(out.y_events[1][0][1]) / 3600.0 if len(out.t_events[1]) else None
    return dict(size_mm=size_mm, relative_humidity=relative_humidity, half_melted_below_m=depth(0),
                melted_below_m=depth(1), hours_to_melt=hours)


def moon_column():
    """The resort's sky ladder about its freezing level, with the climate targets' binding."""
    data, binding = targets.load()
    ladder = targets.read(data)['resort']['climate']['ladder']
    z = np.array([x['height_km'] for x in ladder]) * 1e3
    t = np.array([x['air_c'] for x in ladder]) + T0
    p = np.array([x['pressure_atm'] for x in ladder]) * EARTH_PRESSURE_PA
    near = (z >= LADDER_SPAN_KM[0] * 1e3) & (z <= LADDER_SPAN_KM[1] * 1e3)
    lapse = -np.polyfit(z[near], t[near], 1)[0]
    above = int(np.argmax(t < T0))
    f = (t[above - 1] - T0) / (t[above - 1] - t[above])
    z_f = z[above - 1] + f * (z[above] - z[above - 1])
    p_f = math.exp(math.log(p[above - 1]) + f * (math.log(p[above]) - math.log(p[above - 1])))
    air = earthlike_air(1.2 * EARTH_PRESSURE_PA, 400.0)
    return Column('moon', MOON_SURFACE_GRAVITY, air, z_f, p_f, lapse), binding


def earth_column():
    air = DryAir(EARTH_PRESSURE_PA)
    return Column('earth', STANDARD_GRAVITY, air, 1500.0, EARTH_FREEZING_PA, EARTH_LAPSE_K_M)


def speeds(column, sizes_mm=(0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0)):
    """Terminal speeds (m/s) of flakes of each size in the air at the column's freezing level."""
    s = column.state(column.z_freeze)
    return {f'{d:g} mm': fall_speed(d * 1e-3, flake_mass(d * 1e-3), column.gravity, s)[0] for d in sizes_mm}


def analyse():
    moon, binding = moon_column()
    earth = earth_column()
    tops, sizes, humidities = (-10.0, -15.0, -20.0), (2.0, 5.0, 10.0, 20.0), (0.9, 0.7)
    growth = {c.name: [grow(c, top) for top in tops] for c in (moon, earth)}
    crossing = {c.name: {f'{top:g} C': {f'{s:g} mm': crossing_hours(c, top, s) for s in sizes} for top in tops}
                for c in (moon, earth)}
    melting = {c.name: [melt(c, s, rh) for s in sizes for rh in humidities] for c in (moon, earth)}
    describe = lambda c: dict(gravity_m_s2=c.gravity, freezing_height_m=c.z_freeze, freezing_pressure_pa=c.p_freeze,
                              lapse_k_per_km=c.lapse * 1e3, air_density_at_freezing_kg_m3=c.state(c.z_freeze)['rho'])
    files = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16] for p in CODE}
    return dict(
        schema=SCHEMA, columns=dict(moon=describe(moon), earth=describe(earth)),
        fall_speeds_m_s=dict(moon=speeds(moon), earth=speeds(earth)), vapour_growth=growth, crossing_hours=crossing,
        melting=melting,
        producer=dict(domain='atmosphere', files=files, inputs=binding['inputs']),
        evidence=('Single-particle physics with published laws fitted on Earth (mass, area, fall law, capacitance, '
                  'ventilation, sticking), carried to lunar gravity and air through the quantities they depend on; '
                  'idealised columns and clouds. A first estimate of the sky winter\'s snow, not a cloud model.'),
        reading_rule=('Sizes are maximum dimensions. vapour_growth starts a 0.2 mm crystal at each cloud top and '
                      'grows it from the vapour alone to the freezing level, with the cloud\'s depth and the air\'s '
                      'mass in it, which at one mixing ratio of cloud ice sets the ice a falling flake can sweep up; '
                      'crossing_hours are the hours a flake of a set size takes from the cloud top to the freezing '
                      'level; melting gives, for flakes of set sizes falling from the freezing level through air at a '
                      'set relative humidity, the depths below it where half and all of the ice has melted. The '
                      'Moon\'s column is the resort\'s sky ladder at 301 K; Earth\'s a winter column with its '
                      'freezing level at 850 hPa.'))


def main() -> int:
    out = analyse()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1, default=float) + '\n')
    for name, c in out['columns'].items():
        print(name, {k: round(float(v), 4) for k, v in c.items()})
    for name, v in out['fall_speeds_m_s'].items():
        print(name, 'fall speeds m/s', {k: round(float(x), 3) for k, x in v.items()})
    for name, rows in out['vapour_growth'].items():
        for r in rows:
            print(f"{name} vapour growth from {r['cloud_top_c']:+.0f} C: cloud {r['cloud_depth_m'] / 1e3:.1f} km, "
                  f"{r['cloud_air_kg_m2']:.0f} kg/m2 of air; crystal {r['size_mm']:.2f} mm {r['mass_mg']:.3f} mg at "
                  f"{r['speed_m_s']:.2f} m/s after {r['hours']:.2f} h")
    for name, tops in out['crossing_hours'].items():
        print(name, 'crossing hours', {k: {s: round(h, 2) for s, h in v.items()} for k, v in tops.items()})
    for name, rows in out['melting'].items():
        for r in rows:
            print(f"{name} melting {r['size_mm']:g} mm at rh {r['relative_humidity']}: "
                  f"half {r['half_melted_below_m']}, "
                  f"all {r['melted_below_m']} m below the freezing level, {r['hours_to_melt']} h")
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
