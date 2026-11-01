"""Single raindrops in the lower air: fall speed with drag and shape, ventilated evaporation, survival below cloud.

A drop is an oblate spheroid in its equilibrium shape (Green 1975), falling under gravity and buoyancy against the
drag of a sphere (Clift and Gauvin 1970) raised for the spheroid's surface (Loth 2008), and trading vapour and heat
with the air by ventilated diffusion (Pruppacher and Klett 2010). The equations are those of Loftus and Wordsworth
(2021) and their code, rainprops (MIT licence, kept in vendor/rainprops), adapted on 11 October 2026 for the Open
Moon's rain and brought into the repository on 17 October.

    from atmosphere.rain import drops
    from atmosphere.radiative_convective.thermodynamics import earthlike_air
    column = drops.Column(gravity_m_s2=1.6242, surface_pressure_pa=121590.0, surface_temperature_k=300.0,
                          surface_relative_humidity=0.6, air=earthlike_air(121590.0, 400.0))
    drops.terminal_speed(0.002, column)      # m/s in still air at the ground
    drops.fall(0.002, column)                # from the cloud base to the ground

Sizes are volume-equivalent diameters in metres, and vertical air speed is positive upward. The air below the cloud
is hydrostatic under constant gravity, cools linearly with height (by default at the dry adiabat g/cp, with the
repository's heat capacity of the moist air) and holds a uniform vapour mole fraction, so the cloud base is where that
air saturates: a parcel estimate for a well-mixed layer. The rain leaves the air's temperature and vapour unchanged,
and drops fall independently. The water is pure liquid, with the ground between 273.15 and 313.15 K and the cloud base
at or above the triple point.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, replace

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from shared.constants import GAS_CONSTANT
from atmosphere.radiative_convective.thermodynamics import (CP_VAPOUR, DryAir, MOLAR_MASS, TRIPLE_POINT_K,
                                                            saturation_pressure)

WATER_MOLAR_MASS = MOLAR_MASS['H2O']
# Liquid water near ambient temperature: its density at 20 C, held fixed so that a drop's mass rho pi D^3 / 6 is the
# conserved quantity, and its heat capacity.
WATER_DENSITY = 998.2                     # kg/m3
WATER_HEAT_CAPACITY = 4184.0              # J/kg/K
GROUND_RANGE_K = (273.15, 313.15)         # the warm-water domain at the ground
CUTOFF_DIAMETER = 1e-6                    # m: a drop this small has evaporated
CUTOFF_MASS = WATER_DENSITY * math.pi * CUTOFF_DIAMETER**3 / 6
LARGEST_DIAMETER = 0.03                   # m: the solver's upper bound on size

# Each gas: molar mass (kg/mol), Lennard-Jones collision diameter (angstrom) and well depth over Boltzmann's constant
# (K), and the coefficients of cp = c0 + c1 t + c2 t^2 + c3 t^3 (kJ/kg/K, t = T / 1000 K), as rainprops tabulates them
# (its data/gas_properties.csv, after Reid, Prausnitz and Sherwood 1977 and Borgnakke and Sonntag 2009). Argon, about
# 1% of the dry air, takes nitrogen's transport, heat capacity and molar mass here, which puts this heat capacity of
# the moist air 0.8% above thermodynamics.py's. It enters the Prandtl number of a drop's heat exchange; the column's
# dry adiabat takes the repository's (moist_heat_capacity), and the dry air's own molar mass sets the density.
GASES = {
    'N2': (0.0280134, 3.798, 71.4, (1.11, -0.48, 0.96, -0.42)),
    'O2': (0.031998, 3.467, 106.7, (0.88, -0.0001, 0.54, -0.33)),
    'CO2': (0.04401, 3.941, 195.2, (0.45, 1.67, -1.27, 0.39)),
    'H2O': (WATER_MOLAR_MASS, 2.641, 809.1, (1.79, 0.107, 0.586, -0.2)),
}


def surface_tension(temperature_k):
    """Surface tension of water against its vapour (N/m): the IAPWS (1994) form of Vargaftik et al. (1983)."""
    tau = 1.0 - temperature_k / 647.096
    return 0.2358 * tau**1.256 * (1.0 - 0.625 * tau)


def latent_heat(temperature_k):
    """Latent heat of vaporization (J/kg), a linear fit within 0.01% of IAPWS-95 from 0 to 40 C."""
    return 2.501e6 - 2369.0 * (temperature_k - 273.15)


def liquid_saturation_pressure(temperature_k):
    """Saturation vapour pressure over liquid water (Pa), Wagner and Pruss (2002), from the triple point upward."""
    if not temperature_k >= TRIPLE_POINT_K:
        raise ValueError(f'{temperature_k:.2f} K is below the triple point; the drops are liquid water only')
    return float(saturation_pressure(temperature_k))


def transport(temperature_k, pressure_pa, vapour_fraction, dry_fractions):
    """Viscosity (Pa s), vapour diffusivity (m2/s), thermal conductivity (W/m/K) and heat capacity (J/kg/K).

    Each gas's Chapman-Enskog viscosity with the Neufeld collision integral, Eucken's conductivity, Wilke's mixing
    rule for both, and the vapour's binary diffusion through each gas combined by Blanc's law with the vapour's own
    self-diffusion term, as rainprops computes them (Reid, Prausnitz and Sherwood 1977).
    """
    fractions = dict(dry_fractions)
    fractions['N2'] = fractions.get('N2', 0.0) + fractions.pop('Ar', 0.0)
    names = [n for n in fractions if fractions[n] > 0] + ['H2O']
    x = np.array([fractions[n] * (1 - vapour_fraction) for n in names[:-1]] + [vapour_fraction])
    molar = np.array([GASES[n][0] for n in names])
    sigma = np.array([GASES[n][1] for n in names])
    well = np.array([GASES[n][2] for n in names])
    coefficients = np.array([GASES[n][3] for n in names])
    t_star = temperature_k / well
    omega = 1.16145 / t_star**0.14874 + 0.52487 / np.exp(0.77320 * t_star) + 2.16178 / np.exp(2.43787 * t_star)
    viscosity_i = 1e-7 * 26.69 * np.sqrt(molar * 1e3 * temperature_k) / sigma**2 / omega
    theta = temperature_k / 1000.0
    cp_i = coefficients @ np.array([1.0, theta, theta**2, theta**3]) * 1000.0
    cv_i = cp_i - GAS_CONSTANT / molar
    conductivity_i = viscosity_i * cv_i / 4.0 * (9.0 * cp_i / cv_i - 5.0)
    phi = (1 + np.sqrt(viscosity_i[:, None] / viscosity_i[None, :]) * (molar[None, :] / molar[:, None])**0.25)**2
    phi /= 4.0 / np.sqrt(2.0) * np.sqrt(1 + molar[:, None] / molar[None, :])
    weight = phi @ x
    viscosity = np.sum(x * viscosity_i / weight)
    conductivity = np.sum(x * conductivity_i / weight)
    heat_capacity = np.sum(x * molar * cp_i) / np.sum(x * molar)
    pair_sigma = 0.5 * (sigma + GASES['H2O'][1])
    t_star = temperature_k / np.sqrt(well * GASES['H2O'][2])
    omega_d = (1.06036 / t_star**0.15610 + 0.193 / np.exp(0.47635 * t_star)
               + 1.03587 / np.exp(1.52996 * t_star) + 1.76474 / np.exp(3.89411 * t_star))
    binary = (1.858e-3 * np.sqrt(1 / (molar * 1e3) + 1 / (WATER_MOLAR_MASS * 1e3)) / pair_sigma**2
              * temperature_k**1.5 / (pressure_pa / 101325.0) / omega_d * 1e-4)
    diffusivity = 1.0 / np.sum(x / binary)
    return float(viscosity), float(diffusivity), float(conductivity), float(heat_capacity)


def moist_heat_capacity(vapour_mole_fraction, air):
    """The heat capacity at constant pressure of moist air (J/kg/K) as the repository takes it
    (atmosphere/radiative_convective/thermodynamics.py): the dry air's, from each gas's molar heat capacity, and the
    vapour's, weighted by mass."""
    x = float(vapour_mole_fraction)
    q = x * WATER_MOLAR_MASS / (x * WATER_MOLAR_MASS + (1 - x) * air.molar_mass)
    return (1 - q) * air.cp + q * CP_VAPOUR


@dataclass(frozen=True)
class Column:
    """The air a drop falls through.

    surface_pressure_pa is the total pressure at the ground; air supplies the dry composition and molar mass (its own
    pressure is not used). lapse_k_m defaults to the dry adiabat g/cp with the repository's heat capacity of the moist
    air at the ground (moist_heat_capacity), and cloud_base_m to the height where that air saturates. updraft_m_s is a
    uniform vertical air speed, positive upward.
    """
    gravity_m_s2: float
    surface_pressure_pa: float
    surface_temperature_k: float
    surface_relative_humidity: float
    air: DryAir
    lapse_k_m: float | None = None
    cloud_base_m: float | None = None
    updraft_m_s: float = 0.0

    def __post_init__(self):
        for name in ('gravity_m_s2', 'surface_pressure_pa', 'surface_temperature_k'):
            value = getattr(self, name)
            if not (math.isfinite(value) and value > 0):
                raise ValueError(f'{name} must be finite and positive')
        if not GROUND_RANGE_K[0] <= self.surface_temperature_k <= GROUND_RANGE_K[1]:
            raise ValueError('The ground lies outside the warm-water domain, 273.15-313.15 K')
        if not 0 < self.surface_relative_humidity <= 1:
            raise ValueError('Relative humidity at the ground must lie in (0, 1]')
        if not math.isfinite(self.updraft_m_s):
            raise ValueError('The updraft must be finite')
        self.air.validate()
        vapour = (self.surface_relative_humidity * liquid_saturation_pressure(self.surface_temperature_k)
                  / self.surface_pressure_pa)
        if vapour >= 0.1:
            raise ValueError('Vapour above a tenth of the air lies outside the dilute-vapour domain')
        object.__setattr__(self, '_vapour', vapour)
        object.__setattr__(self, '_molar_mass', (1 - vapour) * self.air.molar_mass + vapour * WATER_MOLAR_MASS)
        cp = moist_heat_capacity(vapour, self.air)
        lapse = self.gravity_m_s2 / cp if self.lapse_k_m is None else float(self.lapse_k_m)
        if not (math.isfinite(lapse) and lapse >= 0):
            raise ValueError('The lapse rate must be finite and non-negative')
        object.__setattr__(self, 'lapse_k_m', lapse)
        if self.cloud_base_m is None:
            object.__setattr__(self, 'cloud_base_m', self._saturation_height())
        if not (math.isfinite(self.cloud_base_m) and self.cloud_base_m >= 0):
            raise ValueError('The cloud base must be finite and non-negative')
        if self.surface_temperature_k - lapse * self.cloud_base_m < TRIPLE_POINT_K - 1e-8:
            raise ValueError('A cloud base below freezing needs ice and melting, which this model leaves out')

    def _saturation_height(self):
        if self.surface_relative_humidity == 1:
            return 0.0
        if self.lapse_k_m == 0:
            raise ValueError('An isothermal column needs its cloud base given')
        top = (self.surface_temperature_k - TRIPLE_POINT_K) / self.lapse_k_m

        def excess(z):
            t, p = self._temperature_pressure(z)
            return self._vapour * p / liquid_saturation_pressure(t) - 1

        if top <= 0 or excess(top) < 0:
            raise ValueError('The air does not saturate above freezing; give a cloud base or another column')
        return float(brentq(excess, 0.0, top))

    def _temperature_pressure(self, z):
        t = self.surface_temperature_k - self.lapse_k_m * z
        exponent = self.gravity_m_s2 * self._molar_mass / GAS_CONSTANT
        if self.lapse_k_m == 0:
            return t, self.surface_pressure_pa * math.exp(-exponent * z / t)
        return t, self.surface_pressure_pa * (t / self.surface_temperature_k)**(exponent / self.lapse_k_m)

    def state(self, altitude_m=0.0):
        """The air at a height: temperature, pressure, vapour, density and the transport properties."""
        z = float(altitude_m)
        if not math.isfinite(z):
            raise ValueError('The altitude must be finite')
        t, p = self._temperature_pressure(z)
        viscosity, diffusivity, conductivity, cp = transport(t, p, self._vapour, self.air.fractions)
        vapour_pressure = self._vapour * p
        return dict(altitude_m=z, temperature_k=t, pressure_pa=p, vapour_pressure_pa=vapour_pressure,
                    relative_humidity=vapour_pressure / liquid_saturation_pressure(t),
                    density_kg_m3=p * self._molar_mass / (GAS_CONSTANT * t), viscosity_pa_s=viscosity,
                    vapour_diffusivity_m2_s=diffusivity, thermal_conductivity_w_m_k=conductivity,
                    heat_capacity_j_kg_k=cp)

    def describe(self):
        """The column's inputs and its derived lapse and cloud base, for a product."""
        out = asdict(self)
        out['air'] = dict(fractions=dict(self.air.fractions), molar_mass_kg_mol=self.air.molar_mass)
        out['vapour_mole_fraction'] = self._vapour
        return out


def transport_adiabat(column):
    """The column with its dry adiabat at the transport model's heat capacity at the ground, which carries argon as
    nitrogen and stands 0.8% above the repository's, and its parcel cloud base: the recovered solver's column of
    11 October."""
    cp = column.state(0.0)['heat_capacity_j_kg_k']
    return replace(column, lapse_k_m=column.gravity_m_s2 / cp, cloud_base_m=None)


def axis_ratio(diameter_m, gravity_m_s2, air_density_kg_m3, temperature_k):
    """Height over width of a drop in its equilibrium oblate shape (Green 1975; Loftus and Wordsworth 2021, eq. 2)."""
    radius = 0.5 * diameter_m
    capillary = math.sqrt(surface_tension(temperature_k) / (gravity_m_s2 * (WATER_DENSITY - air_density_kg_m3)))
    if radius / capillary < 1e-4:
        return 1.0

    def mismatch(q):
        return q**(-1 / 6) * math.sqrt(max(0.0, q**-2 - 2 * q**(-1 / 3) + 1)) - radius / capillary

    return brentq(mismatch, 0.01, 1.0, xtol=2e-12)


def drag(diameter_m, speed_m_s, state, gravity_m_s2, drop_temperature_k, deform=True):
    """Drag force (N) on a drop at an air-relative speed, with its axis ratio and Reynolds number.

    Clift and Gauvin's sphere drag, times Loth's (2008) factor for an oblate spheroid's surface, on the spheroid's
    frontal area (Loftus and Wordsworth 2021, eqs. 5-8). deform=False holds the drop spherical. C_D v^2 is written
    without its 1/v so that a drop at rest starts to fall.
    """
    density, viscosity = state['density_kg_m3'], state['viscosity_pa_s']
    ratio = axis_ratio(diameter_m, gravity_m_s2, density, drop_temperature_k) if deform else 1.0
    area = math.pi * (diameter_m / 2)**2 * ratio**(-2 / 3)
    if ratio < 0.999:
        eccentricity = math.sqrt(1 - ratio**2)
        surface = (ratio**(-2 / 3) / 2
                   + ratio**(4 / 3) / (4 * eccentricity) * math.log((1 + eccentricity) / (1 - eccentricity)))
        shape = 1 + 1.5 * math.sqrt(max(surface - 1, 0.0)) + 6.7 * (surface - 1)
    else:
        shape = 1.0
    v = abs(speed_m_s)
    reynolds = density * v * diameter_m / viscosity
    cd_v2 = (24 * viscosity / (density * diameter_m) * v * (1 + 0.15 * reynolds**(229 / 333))
             + 0.42 * v * v / (1 + 4.25e4 * max(reynolds, 1e-30)**-1.16)) * shape
    return 0.5 * density * area * cd_v2, ratio, reynolds


def terminal_speed(diameter_m, column, altitude_m=0.0, drop_temperature_k=None, deform=True):
    """Fall speed (m/s, downward, relative to the air) at which drag and buoyancy balance the drop's weight."""
    d = float(diameter_m)
    if not (math.isfinite(d) and 0 < d <= LARGEST_DIAMETER):
        raise ValueError('The diameter must lie in (0, 0.03] m')
    state = column.state(altitude_m)
    t = state['temperature_k'] if drop_temperature_k is None else drop_temperature_k
    weight = (WATER_DENSITY - state['density_kg_m3']) * math.pi * d**3 / 6 * column.gravity_m_s2
    return brentq(lambda v: drag(d, v, state, column.gravity_m_s2, t, deform)[0] - weight, 0.0, 100.0, xtol=1e-10)


def stable_diameter(column, altitude_m=0.0):
    """The largest drop the force balance lets hold together, 2 sqrt(3 sigma / (2 g drho)) (Loftus and Wordsworth 2021,
    eq. 22): a screening scale; collisions and turbulence break drops below it."""
    s = column.state(altitude_m)
    return 2 * math.sqrt(1.5 * surface_tension(s['temperature_k'])
                         / (column.gravity_m_s2 * (WATER_DENSITY - s['density_kg_m3'])))


def ventilation(reynolds, number):
    """Ventilation factor for a falling drop (Pruppacher and Klett 2010, eqs. 13-57 and 13-58); number is the
    Schmidt number for vapour and the Prandtl number for heat."""
    x = math.sqrt(max(reynolds, 0.0)) * number**(1 / 3)
    return 1 + 0.108 * x * x if x < 1.4 else 0.78 + 0.308 * x


def drop_rates(y, column, evaporation=True):
    """Rates of change of a falling drop's state y = (height m, downward speed over the ground m/s, mass kg,
    temperature K): its fall, its acceleration under gravity, buoyancy and drag, its ventilated evaporation and its
    heat balance, with the air-relative speed setting drag and ventilation."""
    z, u, m, td = y
    m = max(m, CUTOFF_MASS * 0.05)
    d = (6 * m / (math.pi * WATER_DENSITY))**(1 / 3)
    air = column.state(max(z, -1.0))
    relative = u + column.updraft_m_s
    force, _, reynolds = drag(d, relative, air, column.gravity_m_s2, td)
    du = column.gravity_m_s2 * (1 - air['density_kg_m3'] / WATER_DENSITY) - math.copysign(force, relative) / m
    if not evaporation:
        return [-u, du, 0.0, 0.0]
    schmidt = air['viscosity_pa_s'] / (air['vapour_diffusivity_m2_s'] * air['density_kg_m3'])
    prandtl = air['viscosity_pa_s'] * air['heat_capacity_j_kg_k'] / air['thermal_conductivity_w_m_k']
    dm = (2 * math.pi * d * air['vapour_diffusivity_m2_s'] * ventilation(reynolds, schmidt) * WATER_MOLAR_MASS
          / GAS_CONSTANT * (air['vapour_pressure_pa'] / air['temperature_k']
                            - liquid_saturation_pressure(td) / td))
    heat = 2 * math.pi * d * air['thermal_conductivity_w_m_k'] * ventilation(reynolds, prandtl) * (
        air['temperature_k'] - td)
    return [-u, du, dm, (heat + latent_heat(td) * dm) / (m * WATER_HEAT_CAPACITY)]


def fall(initial_diameter_m, column, sample_heights_m=None, initial_speed_m_s=None, max_time_s=20000.0,
         rtol=2e-6, evaporation=True):
    """Follow one drop from the cloud base to the ground: height, downward speed over the ground, mass and temperature.

    By default the drop leaves the cloud base at its terminal speed less the updraft; initial_speed_m_s releases it
    with a given downward speed, as a drip from a leaf at rest. The fall ends at the ground, when the drop shrinks to
    1 um ('evaporated'), when it stops falling ('lofted') or when it outgrows the stable diameter ('unstable'). A
    drop that evaporates or stops leaves nothing on the ground; a 1 um remnant counts as liquid still aloft.
    """
    d0 = float(initial_diameter_m)
    if not (math.isfinite(d0) and 1e-5 <= d0 <= LARGEST_DIAMETER):
        raise ValueError('The initial diameter must lie between 10 um and 30 mm')
    if not (math.isfinite(max_time_s) and max_time_s > 0):
        raise ValueError('max_time_s must be positive')
    z0 = column.cloud_base_m
    start = column.state(z0)
    v_terminal = terminal_speed(d0, column, z0)
    u0 = v_terminal - column.updraft_m_s if initial_speed_m_s is None else float(initial_speed_m_s)
    if not math.isfinite(u0):
        raise ValueError('The initial speed must be finite')
    if initial_speed_m_s is not None and u0 < 0:
        raise ValueError('A drop must start downward or at rest')
    m0 = WATER_DENSITY * math.pi * d0**3 / 6
    heights = np.linspace(z0, 0, 33) if sample_heights_m is None else np.asarray(sample_heights_m, dtype=float)
    if (heights.ndim != 1 or not np.all(np.isfinite(heights)) or np.any(np.diff(heights) > 0) or np.any(heights < 0)
            or np.any(heights > z0 + 1e-7)):
        raise ValueError('Sample heights must be finite, descending and within [0, cloud base]')
    result = dict(initial_diameter_m=d0, initial_mass_kg=m0, initial_terminal_speed_m_s=v_terminal,
                  initial_speed_m_s=u0, surface_diameter_m=0.0, surface_mass_kg=0.0, surface_speed_m_s=0.0,
                  surface_terminal_speed_m_s=0.0, survived=False, fall_time_s=0.0, evaporated_mass_fraction=0.0,
                  remaining_airborne_mass_fraction=1.0, status='pending')
    if d0 > stable_diameter(column, z0):
        result.update(status='unstable_initial', trajectory=None)
        return result
    if u0 <= 0 and initial_speed_m_s is None:
        result.update(status='lofted_initial', trajectory=None)
        return result
    cutoff_mass = CUTOFF_MASS

    def rates(t, y):
        return drop_rates(y, column, evaporation)

    def ground(t, y):
        return y[0]

    def evaporated(t, y):
        return y[2] - cutoff_mass

    def lofted(t, y):
        return y[1]

    def unstable(t, y):
        return stable_diameter(column, max(y[0], 0.0)) - (max(y[2], 0.0) * 6 / (WATER_DENSITY * math.pi))**(1 / 3)

    for event in (ground, evaporated, lofted, unstable):
        event.terminal, event.direction = True, -1
    # A drop released at rest accelerates downward: its speed crosses zero upward, which the lofted event ignores.
    y0 = [z0, u0, m0, start['temperature_k']]
    if z0 == 0:
        solution, times, final, status = None, np.array([0.0]), np.array(y0), 'surface'
    else:
        solution = solve_ivp(rates, (0, max_time_s), y0, method='Radau', rtol=rtol,
                             atol=[1e-4, 2e-7, max(m0 * 1e-10, 1e-20), 1e-6], max_step=30.0, dense_output=True,
                             events=(ground, evaporated, lofted, unstable))
        if not solution.success:
            raise RuntimeError(solution.message)
        times, final = solution.t, solution.y[:, -1]
        names = ('surface', 'evaporated', 'lofted', 'unstable')
        status = next((n for n, hits in zip(names, solution.t_events) if len(hits)), 'time_limit')
    z, u, m, td = final
    survived = bool(status == 'surface' and u > 0)
    kept = max(m, 0.0) / m0
    final_diameter = float((max(m, 0.0) * 6 / (WATER_DENSITY * math.pi))**(1 / 3))
    result.update(status=status, survived=survived, fall_time_s=float(times[-1]), evaporated_mass_fraction=1 - kept,
                  remaining_airborne_mass_fraction=0.0 if survived else kept, final_altitude_m=float(max(z, 0.0)),
                  final_drop_temperature_k=float(td), final_diameter_m=final_diameter,
                  function_evaluations=0 if solution is None else int(solution.nfev))
    if survived:
        result.update(surface_diameter_m=final_diameter, surface_mass_kg=float(m), surface_speed_m_s=float(u),
                      surface_terminal_speed_m_s=terminal_speed(final_diameter, column, 0.0, td),
                      surface_air_relative_speed_m_s=float(u + column.updraft_m_s))
    track = {k: [] for k in ('altitude_m', 'diameter_m', 'speed_m_s', 'terminal_speed_m_s', 'time_s',
                             'drop_temperature_k')}
    for h in heights:
        track['altitude_m'].append(float(h))
        if h < z - 1e-5:
            for k in track:
                if k != 'altitude_m':
                    track[k].append(None)
            continue
        if abs(h - z0) < 1e-8:
            t, y = 0.0, np.array(y0)
        elif solution is None:
            t, y = 0.0, final
        else:
            t = (brentq(lambda s: solution.sol(s)[0] - h, 0, float(times[-1]), xtol=1e-7) if h > z + 1e-7
                 else float(times[-1]))
            y = solution.sol(t)
        d = float((max(y[2], 0.0) * 6 / (WATER_DENSITY * math.pi))**(1 / 3))
        track['diameter_m'].append(d)
        track['speed_m_s'].append(float(y[1]))
        track['terminal_speed_m_s'].append(terminal_speed(d, column, h, float(y[3])))
        track['time_s'].append(float(t))
        track['drop_temperature_k'].append(float(y[3]))
    result['trajectory'] = track
    return result
