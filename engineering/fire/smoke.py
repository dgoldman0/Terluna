"""Fire and smoke in buildings under any gravity and air: the standard fire-engineering correlations, with gravity
and the ambient air written in.

Fire engineering's correlations were fitted in Earth laboratories. Their physics is buoyancy, so they follow Froude
scaling (Quintiere 1989): two fires in spaces of the same shape and size behave alike when their dimensionless heat
release Q* = Q / (rho cp T sqrt(g) L^(5/2)) is the same, with times scaling as sqrt(L/g) and velocities as
sqrt(g L). In the same space, a fire in weaker gravity or thinner air therefore acts as a larger fire would on
Earth, played more slowly: at lunar gravity in air of Earth's density, as a fire of about 2.5 times the heat release
over about 2.5 times the time. `froude` gives the factors. The Earth correlations below are carried to any gravity
and air through them, or written with gravity explicit where the source already is:

- flame height (Heskestad), with the air's oxygen, which sets how much air a flame must take in to burn its fuel;
- the axisymmetric plume's mass flow (NFPA 92);
- the ceiling jet of a fire growing as Q = alpha t^2, and the response of smoke detectors, heat detectors and
  sprinklers to it (Heskestad and Delichatsios, as in NFPA 72 Annex B, and the RTI model);
- a two-zone model of smoke filling a closed space, and the exhaust or natural vent area that holds the smoke above
  a clear height;
- compartment fires: the hot layer's temperature and the flashover threshold (McCaffrey, Quintiere and Harkleroad),
  the ventilation-controlled heat release (EN 1991-1-2 Annex E) with the openings' air supply scaled to the gravity
  and air, and the larger supply that wind blowing through failed openings can give;
- the stack pressure along a shaft, and the pressure at the foot of a water column.

Scaling carries the buoyant flow. It does not carry the fire itself: how fast a fire grows and how much heat it
releases per unit area at another gravity and oxygen fraction are inputs (design_fires.py). Radiation, heat lost to
walls, soot, and the change from turbulent to laminar flow at low Grashof numbers do not follow Froude scaling, so
the results are first-order, as the correlations themselves are on Earth.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

from atmosphere.radiative_convective import thermodynamics as th
from shared.constants import GAS_CONSTANT, STANDARD_GRAVITY

HEAT_PER_OXYGEN_KJ_KG = 13_100.0   # heat released per kg of oxygen consumed, nearly the same for most fuels (Huggett)
WATER_DENSITY = 1000.0


@dataclass(frozen=True)
class Air:
    """Still air around a fire: gravity, pressure, temperature and dry composition."""
    gravity_m_s2: float
    pressure_pa: float
    temperature_k: float
    dry: th.DryAir = field(default_factory=lambda: th.DryAir(101325.0))

    @property
    def density(self) -> float:
        return self.pressure_pa * self.dry.molar_mass / (GAS_CONSTANT * self.temperature_k)

    @property
    def cp_kj(self) -> float:
        return self.dry.cp / 1000.0

    @property
    def oxygen_mole_fraction(self) -> float:
        return self.dry.fractions['O2']

    @property
    def oxygen_mass_fraction(self) -> float:
        return self.oxygen_mole_fraction * th.MOLAR_MASS['O2'] / self.dry.molar_mass

    @property
    def heat_per_air_kj_kg(self) -> float:
        """Heat released per kg of air burned: the oxygen it carries times 13.1 MJ per kg of oxygen."""
        return HEAT_PER_OXYGEN_KJ_KG * self.oxygen_mass_fraction


EARTH_REFERENCE = Air(STANDARD_GRAVITY, 101325.0, 293.15)   # the laboratory air of the correlations


def froude(air: Air, reference: Air = EARTH_REFERENCE) -> dict:
    """Factors that carry a correlation fitted in `reference` air to `air`, in a space of the same size.

    heat: a fire of heat release Q in `air` has the same Q* as a fire of heat x Q in the reference. The reference
    correlation's times, velocities, mass flows and temperature rises for that fire convert back to `air` when
    multiplied by time, velocity, mass_flow and temperature."""
    def b(a):
        return a.density * a.cp_kj * a.temperature_k * np.sqrt(a.gravity_m_s2)
    return dict(heat=b(reference) / b(air),
                time=float(np.sqrt(reference.gravity_m_s2 / air.gravity_m_s2)),
                velocity=float(np.sqrt(air.gravity_m_s2 / reference.gravity_m_s2)),
                mass_flow=air.density * np.sqrt(air.gravity_m_s2) / (reference.density * np.sqrt(reference.gravity_m_s2)),
                temperature=air.temperature_k / reference.temperature_k)


def flame_height_m(q_kw, diameter_m, air: Air):
    """Mean flame height of a fire of heat release q (kW) on a base of the given diameter (Heskestad):
    L = D (15.6 N^(1/5) - 1.02), N = cp T Q^2 / (g rho^2 (dHc/r)^3 D^5), with dHc/r the heat released per kg of air
    burned. In Earth's air this is L = 0.235 Q^(2/5) - 1.02 D. Weaker gravity and thinner air lengthen the flame, and
    so does air with less oxygen, of which the flame must take in more to burn the same fuel."""
    q, d = np.asarray(q_kw, float), np.asarray(diameter_m, float)
    n = air.cp_kj * air.temperature_k * q ** 2 / (air.gravity_m_s2 * air.density ** 2 * air.heat_per_air_kj_kg ** 3
                                                   * d ** 5)
    return d * (15.6 * n ** 0.2 - 1.02)


def _plume_earth(q, z):
    """NFPA 92's axisymmetric plume in Earth air: kg/s at height z (m) above the fire's base for a convective heat
    release q (kW), above and within the flame (limiting height 0.166 q^(2/5))."""
    if z > 0.166 * q ** 0.4:
        return 0.071 * q ** (1 / 3) * z ** (5 / 3) + 0.0018 * q
    return 0.032 * q ** 0.6 * z


def plume_mass_flow_kg_s(qc_kw, height_m, air: Air):
    """Mass flow of an axisymmetric plume at a height above the fire's base, for a convective heat release qc (kW):
    NFPA 92's m = 0.071 Qc^(1/3) z^(5/3) + 0.0018 Qc above the flame and 0.032 Qc^(3/5) z within it, carried to the
    air by Froude scaling."""
    f = froude(air)
    q = np.asarray(qc_kw, float) * f['heat']
    z = np.asarray(height_m, float)
    m = np.where(z > 0.166 * q ** 0.4, 0.071 * q ** (1 / 3) * z ** (5 / 3) + 0.0018 * q, 0.032 * q ** 0.6 * z)
    return m * f['mass_flow']


# Ceiling jet of a t-squared fire: Evans and Stroup (1985, equation 2), from Heskestad and Delichatsios's tests with
# wood cribs. NFPA 72 Annex B lowered the heat release behind those tests by a factor of 1.67, which turns the arrival
# and temperature constants 0.954, 0.188 and 0.313 into 0.861, 0.146 and 0.242 for a fire's actual heat release.
HD_ARRIVAL, HD_OFFSET, HD_SLOPE = 0.954, 0.188, 0.313
HD_VELOCITY, HD_VELOCITY_EXPONENT = 0.59, -0.63
HD_HEAT_CORRECTION = 1.67
SMOKE_DETECTOR_RISE_K = 13.0   # the ceiling-jet temperature rise used as a surrogate for smoke detection (NFPA 72)


def hd_constants(correction: float = HD_HEAT_CORRECTION) -> tuple[float, float, float]:
    """The arrival, offset and slope constants for a fire's actual heat release: the dimensionless time scales as
    alpha^(1/5) and the temperature as alpha^(-2/5), so the arrival scales as correction^(-1/5) and the others as
    correction^(-1/2)."""
    return HD_ARRIVAL * correction ** -0.2, HD_OFFSET * correction ** -0.5, HD_SLOPE * correction ** -0.5


def _growth_scale(alpha_kw_s2, air: Air):
    """A alpha, with A = g / (cp T rho): the combination through which gravity and the air enter (m^4/s^5)."""
    return air.gravity_m_s2 / (air.cp_kj * air.temperature_k * air.density) * np.asarray(alpha_kw_s2, float)


def ceiling_jet(time_s, alpha_kw_s2, ceiling_m, radius_m, air: Air):
    """Temperature rise (K) and velocity (m/s) of the ceiling jet at a radius from a fire growing as Q = alpha t^2
    (total heat release, kW), with the ceiling at the given height above the fuel."""
    arrival, offset, slope = hd_constants()
    t, h, r = (np.asarray(v, float) for v in (time_s, ceiling_m, radius_m))
    aa = _growth_scale(alpha_kw_s2, air)
    rh = r / h
    t_star = t * aa ** 0.2 * h ** -0.8
    rise_star = np.clip((t_star - arrival * (1 + rh)) / (offset + slope * rh), 0.0, None) ** (4 / 3)
    rise = rise_star * aa ** 0.4 * (air.temperature_k / air.gravity_m_s2) * h ** -0.6
    velocity = HD_VELOCITY * rh ** HD_VELOCITY_EXPONENT * np.sqrt(rise_star) * aa ** 0.2 * h ** 0.2
    return rise, velocity


def gas_rise_time_s(rise_k, alpha_kw_s2, ceiling_m, radius_m, air: Air):
    """When the ceiling jet at a radius first warms by rise_k: the response of a detector with no thermal lag, and
    with SMOKE_DETECTOR_RISE_K the usual surrogate for a smoke detector."""
    arrival, offset, slope = hd_constants()
    h, r = np.asarray(ceiling_m, float), np.asarray(radius_m, float)
    aa = _growth_scale(alpha_kw_s2, air)
    rh = r / h
    rise_star = np.asarray(rise_k, float) / (aa ** 0.4 * (air.temperature_k / air.gravity_m_s2) * h ** -0.6)
    t_star = arrival * (1 + rh) + (offset + slope * rh) * rise_star ** 0.75
    return t_star / (aa ** 0.2 * h ** -0.8)


def sensor_rise(times_s, gas_rise_k, velocity_m_s, rti):
    """Temperature rise of a sensing element with response time index rti ((m s)^1/2) in a gas stream:
    dTd/dt = sqrt(u) (Tg - Td) / RTI, stepped exactly over each interval with the gas held at its end value."""
    t, g, u = (np.asarray(v, float) for v in (times_s, gas_rise_k, velocity_m_s))
    out = np.zeros_like(g)
    for k in range(1, t.size):
        decay = np.exp(-np.sqrt(u[k]) * (t[k] - t[k - 1]) / rti)
        out[k] = g[k] + (out[k - 1] - g[k]) * decay
    return out


def activation_time_s(rise_k, rti, alpha_kw_s2, ceiling_m, radius_m, air: Air, t_max_s: float = 3600.0,
                      step_s: float = 0.25) -> float:
    """Time for a heat detector or sprinkler bulb (rti, (m s)^1/2) to warm by rise_k above the ambient in a growing
    fire's ceiling jet; inf if not within t_max_s."""
    t = np.arange(0.0, t_max_s + step_s, step_s)
    gas, u = ceiling_jet(t, alpha_kw_s2, ceiling_m, radius_m, air)
    hit = np.nonzero(sensor_rise(t, gas, u, rti) >= rise_k)[0]
    return float(t[hit[0]]) if hit.size else float('inf')


def smoke_filling(hrr_kw, area_m2: float, ceiling_m: float, air: Air, clear_m: float = 2.0,
                  convective: float = 0.7, loss: float = 0.3, hot_limit_k: float = 180.0,
                  t_max_s: float = 3600.0, step_s: float = 0.5) -> dict:
    """Smoke filling a closed space from a fire at floor level, as a two-zone model: the plume carries air into a hot
    layer under the ceiling, which keeps the convective heat less the share `loss` given up to the ceiling and walls.

    hrr_kw is the fire's total heat release as a function of time (s); `convective` of it drives the plume. Returns
    the times the layer's base falls to clear_m (smoke at head height) and its temperature rise reaches hot_limit_k
    (a layer at about 200 C radiates 2.5 kW/m2 to the people below), inf when not within t_max_s; `untenable_s` is
    the earlier of the two."""
    f = froude(air)
    rho, cp, t_amb = air.density, air.cp_kj, air.temperature_k
    z, mass, heat = ceiling_m, 0.0, 0.0
    t_clear = t_hot = float('inf')
    rise_at_clear = None
    t = 0.0
    while t < t_max_s and (t_clear == float('inf') or t_hot == float('inf')):
        q = convective * float(hrr_kw(t + 0.5 * step_s))
        t += step_s
        if q > 0.0:
            mass += _plume_earth(q * f['heat'], max(z, 0.05)) * f['mass_flow'] * step_s
            heat += (1.0 - loss) * q * step_s
        if mass <= 0.0:
            continue
        rise = heat / (cp * mass)
        z = ceiling_m - mass * (t_amb + rise) / (rho * t_amb) / area_m2
        if z <= clear_m and t_clear == float('inf'):
            t_clear, rise_at_clear = t, rise
        if rise >= hot_limit_k and t_hot == float('inf'):
            t_hot = t
        if z <= 0.0:
            break
    governed = None if min(t_clear, t_hot) == float('inf') else ('smoke at head height' if t_clear <= t_hot
                                                                 else 'hot layer overhead')
    return dict(clear_s=t_clear, hot_s=t_hot, untenable_s=min(t_clear, t_hot), governed_by=governed,
                layer_rise_at_clear_k=None if rise_at_clear is None else round(float(rise_at_clear), 1))


def exhaust_to_hold(qc_kw, clear_m, air: Air, loss: float = 0.3) -> dict:
    """The steady exhaust that holds a smoke layer's base at clear_m above a fire of convective heat release qc: the
    plume's mass flow at that height (NFPA 92's steady method), the layer's temperature rise, and its volume flow."""
    m = float(plume_mass_flow_kg_s(qc_kw, clear_m, air))
    rise = (1.0 - loss) * float(qc_kw) / (air.cp_kj * m)
    volume = m * (air.temperature_k + rise) / (air.density * air.temperature_k)
    return dict(mass_kg_s=m, layer_rise_k=rise, volume_m3_s=volume)


def natural_vent_area_m2(mass_kg_s, layer_depth_m, layer_rise_k, air: Air, discharge: float = 0.6):
    """Area of roof or wall vents through which a hot layer of this depth and temperature rise drives mass_kg_s out
    by its own buoyancy, with ample inlets below: m = Cd A rho sqrt(2 g d theta T) / (T + theta)."""
    t = air.temperature_k
    return mass_kg_s * (t + layer_rise_k) / (discharge * air.density
                                             * np.sqrt(2.0 * air.gravity_m_s2 * layer_depth_m * layer_rise_k * t))


# McCaffrey, Quintiere and Harkleroad (1981): the hot layer of a room with one opening, 6.85 (Q^2/(Ao sqrt(Ho) hk AT))^(1/3)
# in Earth air (Q kW, Ao m2, Ho m, hk kW/m2/K, AT m2); flashover at a rise of about 500 K.
MQH_EARTH = 6.85
FLASHOVER_RISE_K = 500.0


def _mqh_groups(opening_area_m2, opening_height_m, air: Air):
    return np.sqrt(air.gravity_m_s2) * air.cp_kj * air.density * opening_area_m2 * np.sqrt(opening_height_m)


def _mqh_coefficient() -> float:
    """The dimensionless coefficient C in dT/T = C (Q/(b T))^(2/3) (hk AT/b)^(-1/3), b = sqrt(g) cp rho Ao sqrt(Ho),
    set so that the general form gives 6.85 in Earth air."""
    a = EARTH_REFERENCE
    b = np.sqrt(a.gravity_m_s2) * a.cp_kj * a.density
    return MQH_EARTH / (a.temperature_k * (b * a.temperature_k) ** (-2 / 3) * b ** (1 / 3))


def layer_rise_k(q_kw, opening_area_m2, opening_height_m, lining_kw_m2_k, surface_area_m2, air: Air):
    """Hot-layer temperature rise in a room with one opening (MQH), written with gravity and the air explicit."""
    b = _mqh_groups(opening_area_m2, opening_height_m, air)
    t = air.temperature_k
    return _mqh_coefficient() * t * (np.asarray(q_kw, float) / (b * t)) ** (2 / 3) \
        * (lining_kw_m2_k * surface_area_m2 / b) ** (-1 / 3)


def flashover_kw(opening_area_m2, opening_height_m, lining_kw_m2_k, surface_area_m2, air: Air,
                 rise_k: float = FLASHOVER_RISE_K):
    """The heat release that takes a room's hot layer to flashover (MQH inverted); 610 (hk AT Ao sqrt(Ho))^(1/2) kW
    in Earth air, and less at lower gravity, where the opening lets the hot gas out more slowly."""
    b = _mqh_groups(opening_area_m2, opening_height_m, air)
    t = air.temperature_k
    return b * t * (rise_k / (_mqh_coefficient() * t)) ** 1.5 * (lining_kw_m2_k * surface_area_m2 / b) ** 0.5


def ventilation_limit_kw(opening_area_m2, opening_height_m, air: Air):
    """The largest heat release a room's openings can feed with air: EN 1991-1-2's 0.10 m Hu Av sqrt(heq) MW for
    wood (m = 0.8, Hu = 17.5 MJ/kg: 1.4 MW per m^(5/2) in Earth air), with the inflow of air scaled by rho sqrt(g)
    and the heat per kg of air by the air's oxygen."""
    ref = EARTH_REFERENCE
    earth = 0.10 * 0.8 * 17.5e3 * np.asarray(opening_area_m2, float) * np.sqrt(opening_height_m)
    inflow = air.density * np.sqrt(air.gravity_m_s2) / (ref.density * np.sqrt(ref.gravity_m_s2))
    return earth * inflow * air.oxygen_mass_fraction / ref.oxygen_mass_fraction


def wind_ventilation_limit_kw(opening_area_m2, wind_m_s, air: Air, discharge: float = 0.6,
                              pressure_coefficient_difference: float = 1.0):
    """The heat release that wind blowing through a burning space can feed: equal failed openings on the windward and
    leeward faces in series (effective area A / sqrt(2)), driven by a pressure difference of dCp x 1/2 rho v^2, with
    all the oxygen that passes burned. Wind does not weaken with gravity, so where a space's openings face the wind,
    this can exceed the buoyant limit many times over, as in wind-driven fires in Earth's tall buildings."""
    flow = discharge * np.asarray(opening_area_m2, float) / np.sqrt(2.0) * np.sqrt(pressure_coefficient_difference) \
        * air.density * np.asarray(wind_m_s, float)
    return flow * air.heat_per_air_kj_kg


def stack_pressure_pa(height_m, outside: Air, inside_temperature_k: float):
    """Stack effect: the pressure difference that builds along a shaft of this height, g h (rho_out - rho_in), with
    the inside air at the same pressure; positive when the inside is warmer and air is drawn in at the bottom."""
    inside = outside.density * outside.temperature_k / inside_temperature_k
    return outside.gravity_m_s2 * np.asarray(height_m, float) * (outside.density - inside)


def water_pressure_pa(height_m, gravity_m_s2: float):
    """Pressure at the foot of a water column of this height."""
    return WATER_DENSITY * gravity_m_s2 * np.asarray(height_m, float)
