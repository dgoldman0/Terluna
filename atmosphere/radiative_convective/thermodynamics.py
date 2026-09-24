"""Moist thermodynamics and the prescribed column of the inverse climate calculation.

A column is fixed by its surface temperature: a saturated pseudoadiabat
(non-dilute, so it stays valid when water vapour is a large fraction of the
air) up to a prescribed stratospheric temperature, then an isothermal
stratosphere holding the tropopause water mixing ratio. Humidity follows a
chosen relative-humidity profile on that temperature structure, as in the
standard one-dimensional habitable-zone calculations (Kasting 1988; Kopparapu et
al. 2013; Koll & Cronin 2018). Gravity falls as GM/r^2, and the dry air is
specified by its surface partial pressure and composition, so the same code
serves the Moon and an Earth control.

This is a prescribed structure for radiative calculations, not a solved
radiative-convective equilibrium: the stratospheric temperature and the
humidity profile are inputs.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import math
import numpy as np

from shared.constants import (AVOGADRO, BOLTZMANN, GAS_CONSTANT, MOON_GM, MOON_RADIUS,
                              EARTH_GM, EARTH_RADIUS)

# Molar masses (kg/mol) and ideal-gas heat capacities near 300 K (J/mol/K),
# NIST Chemistry WebBook. Held constant with temperature.
MOLAR_MASS = {'N2': 0.0280134, 'O2': 0.0319988, 'Ar': 0.039948, 'CO2': 0.0440095, 'H2O': 0.01801528}
MOLAR_CP = {'N2': 29.12, 'O2': 29.38, 'Ar': 20.786, 'CO2': 37.12, 'H2O': 33.58}
R_VAPOUR = GAS_CONSTANT / MOLAR_MASS['H2O']
CP_VAPOUR = MOLAR_CP['H2O'] / MOLAR_MASS['H2O']
TRIPLE_POINT_K = 273.16

# Wagner & Pruss (2002), J. Phys. Chem. Ref. Data 31, 387, eq. 2.5: saturation
# pressure over liquid water, 273.16 K to the critical point.
_TC, _PC = 647.096, 22.064e6
_WP = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
_WP_EXP = (1.0, 1.5, 3.0, 3.5, 4.0, 7.5)
# Wagner, Riethmann, Feistel & Harvey (2011), J. Phys. Chem. Ref. Data 40,
# 043103, eq. 6: sublimation pressure over hexagonal ice below the triple point.
_PT = 611.657
_ICE_A = (-21.2144006, 27.3203819, -6.10598130)
_ICE_B = (0.00333333333, 1.20666667, 1.70333333)


def _ln_esat_liquid(t):
    tau = 1 - t / _TC
    return np.log(_PC) + _TC / t * sum(a * tau**b for a, b in zip(_WP, _WP_EXP))


def _dln_esat_liquid(t):
    tau = 1 - t / _TC
    series = sum(a * tau**b for a, b in zip(_WP, _WP_EXP))
    dseries = sum(-a * b * tau**(b - 1) / _TC for a, b in zip(_WP, _WP_EXP))
    return -_TC / t**2 * series + _TC / t * dseries


def _ln_esat_ice(t):
    theta = t / TRIPLE_POINT_K
    return np.log(_PT) + sum(a * theta**b for a, b in zip(_ICE_A, _ICE_B)) / theta


def _dln_esat_ice(t):
    theta = t / TRIPLE_POINT_K
    s = sum(a * theta**b for a, b in zip(_ICE_A, _ICE_B))
    ds = sum(a * b * theta**(b - 1) for a, b in zip(_ICE_A, _ICE_B))
    return (ds / theta - s / theta**2) / TRIPLE_POINT_K


def _check(t):
    t = np.asarray(t, dtype=float)
    if np.any(t < 100) or np.any(t >= _TC):
        raise ValueError('Temperature outside 100 K to the critical point')
    return t


def saturation_pressure(t):
    """Saturation vapour pressure (Pa) over liquid at or above 273.16 K, over ice below."""
    t = _check(t)
    return np.exp(np.where(t >= TRIPLE_POINT_K, _ln_esat_liquid(np.maximum(t, TRIPLE_POINT_K)),
                           _ln_esat_ice(np.minimum(t, TRIPLE_POINT_K))))


def latent_heat(t):
    """Clausius-Clapeyron latent heat (J/kg) consistent with saturation_pressure.

    L = R_v T^2 dln(e_s)/dT for an ideal vapour and negligible condensate volume:
    vaporization at or above the triple point, sublimation below it.
    """
    t = _check(t)
    d = np.where(t >= TRIPLE_POINT_K, _dln_esat_liquid(np.maximum(t, TRIPLE_POINT_K)),
                 _dln_esat_ice(np.minimum(t, TRIPLE_POINT_K)))
    return R_VAPOUR * t**2 * d


def _esat_and_latent(t):
    """Scalar fast path used inside the adiabat integration."""
    if t >= TRIPLE_POINT_K:
        return math.exp(_ln_esat_liquid(t)), R_VAPOUR * t * t * _dln_esat_liquid(t)
    return math.exp(_ln_esat_ice(t)), R_VAPOUR * t * t * _dln_esat_ice(t)


@dataclass(frozen=True)
class Planet:
    name: str
    radius_m: float
    gm_m3_s2: float

    @property
    def surface_gravity(self):
        return self.gm_m3_s2 / self.radius_m**2

    def gravity(self, z):
        return self.gm_m3_s2 / (self.radius_m + np.asarray(z, dtype=float))**2


MOON = Planet('moon', MOON_RADIUS, MOON_GM)
EARTH = Planet('earth', EARTH_RADIUS, EARTH_GM)


@dataclass(frozen=True)
class DryAir:
    """Non-condensing air: surface partial pressure and dry mole fractions."""
    surface_pressure_pa: float
    fractions: dict = field(default_factory=lambda: {'N2': 0.7808, 'O2': 0.2095, 'Ar': 0.0093, 'CO2': 0.0004})

    def validate(self):
        if not self.surface_pressure_pa > 0:
            raise ValueError('Dry surface pressure must be positive')
        if set(self.fractions) - set(MOLAR_MASS) or 'H2O' in self.fractions:
            raise ValueError('Dry air holds N2, O2, Ar and CO2 only')
        if any(v < 0 for v in self.fractions.values()) or abs(sum(self.fractions.values()) - 1) > 1e-9:
            raise ValueError('Dry mole fractions must be non-negative and sum to one')

    @property
    def molar_mass(self):
        return sum(x * MOLAR_MASS[s] for s, x in self.fractions.items())

    @property
    def gas_constant(self):
        return GAS_CONSTANT / self.molar_mass

    @property
    def cp(self):
        """Specific heat at constant pressure, J/kg/K."""
        return sum(x * MOLAR_CP[s] for s, x in self.fractions.items()) / self.molar_mass


def earthlike_air(surface_pressure_pa, co2_ppm, oxygen_partial_pa=21_227.0, argon_ratio=0.0093 / 0.7808):
    """Dry air holding Earth's oxygen partial pressure, Earth's Ar/N2 ratio and a CO2 amount.

    The default oxygen partial pressure is 0.2095 of one standard atmosphere, so a
    1.2-atm column carries 17.46% O2, as in the repository's 17.5% design case.
    """
    o2 = oxygen_partial_pa / surface_pressure_pa
    co2 = co2_ppm * 1e-6
    rest = 1 - o2 - co2
    if rest <= 0:
        raise ValueError('Oxygen and CO2 exceed the dry pressure')
    n2 = rest / (1 + argon_ratio)
    return DryAir(surface_pressure_pa, {'N2': n2, 'O2': o2, 'Ar': rest - n2, 'CO2': co2})


def manabe_wetherald(rh_surface=0.77, floor=0.02):
    """Relative humidity RH_s (p/p_s - floor)/(1 - floor) (Manabe & Wetherald 1967)."""
    def rh(p, ps):
        return np.clip(rh_surface * (p / ps - floor) / (1 - floor), 0.0, 1.0)
    rh.label = f'manabe_wetherald_{rh_surface:g}'
    return rh


def uniform_rh(value):
    def rh(p, ps):
        return np.full_like(np.asarray(p, dtype=float), value)
    rh.label = f'uniform_{value:g}'
    return rh


def adiabat_slope(t, p, air: DryAir):
    """dlnT/dlnp on the saturated pseudoadiabat, total pressure p.

    Ding & Pierrehumbert (2016), ApJ 822, 24, eq. 12, with the latent heat taken
    from the Clausius-Clapeyron slope of the saturation curve. It reduces to the
    familiar dilute moist adiabat when vapour is scarce and to the saturation
    curve itself for a pure-vapour atmosphere.
    """
    ec, lv = _esat_and_latent(t)
    pa = p - ec
    if pa <= 0:
        return R_VAPOUR * t / lv
    ra, cpa = air.gas_constant, air.cp
    r_sat = ra / R_VAPOUR * ec / pa
    alpha, beta, gamma = lv / (ra * t), lv / (R_VAPOUR * t), lv / (cpa * t)
    moist = (1 + (CP_VAPOUR / cpa + (beta - 1) * gamma) * r_sat) / (1 + alpha * r_sat)
    return 1.0 / (ec / p * beta + pa / p * cpa / ra * moist)


def moist_adiabat(ts, ps, air: DryAir, t_stop=100.0, x_max=40.0, step=0.002):
    """Saturated pseudoadiabat from (ts, ps) upward, by fixed-step RK4 in x = ln(ps/p).

    Integration stops once T falls below t_stop. Returns (x, T) samples; use
    adiabat_temperature to evaluate them at chosen pressures.
    """
    h = step
    xs, lnt = [0.0], [math.log(ts)]

    def f(x, y):
        return -adiabat_slope(math.exp(y), ps * math.exp(-x), air)

    x, y = 0.0, math.log(ts)
    while x < x_max and math.exp(y) > t_stop:
        k1 = f(x, y); k2 = f(x + h / 2, y + h / 2 * k1)
        k3 = f(x + h / 2, y + h / 2 * k2); k4 = f(x + h, y + h * k3)
        y += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        x += h
        xs.append(x); lnt.append(y)
    return np.asarray(xs), np.exp(np.asarray(lnt))


def adiabat_temperature(adiabat, ps, p):
    """Interpolate an adiabat from moist_adiabat at pressures p (log-linear)."""
    xs, ts = adiabat
    x = np.log(ps / np.asarray(p, dtype=float))
    return np.exp(np.interp(x, xs, np.log(ts)))


@dataclass(frozen=True)
class ColumnSpec:
    planet: Planet
    air: DryAir
    surface_temperature_k: float
    stratosphere_k: float = 200.0
    humidity: object = None          # callable rh(p, ps); default Manabe-Wetherald 0.77
    top_pressure_pa: float = 1.0
    troposphere_levels: int = 70     # between the surface and 1% of it
    upper_levels: int = 20           # from there to the top

    def validate(self):
        self.air.validate()
        if not self.stratosphere_k <= self.surface_temperature_k < _TC:
            raise ValueError('Need T_strat <= T_s below the critical point')
        if self.top_pressure_pa <= 0 or self.troposphere_levels < 8 or self.upper_levels < 2:
            raise ValueError('Invalid vertical grid')


def build_column(spec: ColumnSpec) -> dict:
    """Levels (surface first) and layers of the prescribed column.

    Returns pressure p_pa, temperature t_k, water mole fraction x_h2o, height z_m
    and gravity g at levels; for layers, the mid-pressure, mid-temperature,
    layer water fraction, layer total molecular column (molecules/cm^2) and mass
    (kg/m^2). The dry-air composition is uniform; CO2, N2, O2 and Ar fractions
    in moist air are their dry fractions times (1 - x_h2o).
    """
    spec.validate()
    rh = spec.humidity or manabe_wetherald()
    air, ts = spec.air, spec.surface_temperature_k
    es_s = float(saturation_pressure(ts))
    ps = air.surface_pressure_pa + float(rh(np.array([1.0]), 1.0)[0]) * es_s
    p_mid = 0.01 * ps
    if spec.top_pressure_pa >= p_mid:
        raise ValueError('Top pressure must lie above 1% of the surface pressure')
    lower = np.exp(np.linspace(math.log(ps), math.log(p_mid), spec.troposphere_levels + 1))
    upper = np.exp(np.linspace(math.log(p_mid), math.log(spec.top_pressure_pa), spec.upper_levels + 1))[1:]
    p = np.concatenate([lower, upper])
    adiabat = moist_adiabat(ts, ps, air, t_stop=0.999 * spec.stratosphere_k,
                            x_max=math.log(ps / spec.top_pressure_pa))
    xs, ta = adiabat
    if ta[-1] < spec.stratosphere_k:
        # Tropopause: where the adiabat reaches the stratospheric temperature; made a level.
        k = int(np.argmax(ta < spec.stratosphere_k))
        w = (ta[k - 1] - spec.stratosphere_k) / (ta[k - 1] - ta[k])
        p_trop = ps * math.exp(-((1 - w) * xs[k - 1] + w * xs[k]))
        if np.min(np.abs(np.log(p / p_trop))) > 1e-6:
            p = np.sort(np.append(p, p_trop))[::-1]
    else:
        p_trop = spec.top_pressure_pa
    t_adiabat = adiabat_temperature(adiabat, ps, p)
    strat = p <= p_trop * (1 + 1e-12)
    t = np.where(strat, spec.stratosphere_k, t_adiabat)
    t[0] = ts
    t_trop = float(adiabat_temperature(adiabat, ps, p_trop))
    x_trop = float(rh(np.array([p_trop]), ps)[0]) * float(saturation_pressure(t_trop)) / p_trop
    x_h2o = np.where(strat, x_trop, rh(p, ps) * saturation_pressure(np.maximum(t, 100.0)) / p)
    x_h2o = np.clip(x_h2o, 0.0, 1.0)

    col = column_from_levels(spec.planet, air, p, t, x_h2o)
    col.update(surface_pressure_pa=ps, tropopause_pa=p_trop, tropopause_h2o=x_trop,
               humidity=getattr(rh, 'label', 'custom'), spec=spec)
    return col


def column_from_levels(planet: Planet, air: DryAir, p, t, x_h2o, constant_gravity=False) -> dict:
    """Hydrostatic heights and layer columns for given levels (surface first).

    Layers are bounded by consecutive levels; their pressure is the log-pressure
    midpoint and their mass dp/g at mid-height (g = GM/r^2, or the surface value
    when constant_gravity is set, as in plane-parallel codes).
    """
    p, t, x_h2o = (np.asarray(v, dtype=float) for v in (p, t, x_h2o))
    if not (p.ndim == t.ndim == x_h2o.ndim == 1 and p.size == t.size == x_h2o.size >= 2):
        raise ValueError('Levels must be matching one-dimensional arrays')
    if np.any(np.diff(p) >= 0):
        raise ValueError('Pressure must decrease upward from the surface')
    m_dry = air.molar_mass
    molar = m_dry * (1 - x_h2o) + MOLAR_MASS['H2O'] * x_h2o
    r_specific = GAS_CONSTANT / molar
    dphi = -0.5 * (r_specific[1:] * t[1:] + r_specific[:-1] * t[:-1]) * np.diff(np.log(p))
    phi = np.concatenate([[0.0], np.cumsum(dphi)])
    g0, rad = planet.surface_gravity, planet.radius_m
    if constant_gravity:
        z = phi / g0
        gravity = lambda height: np.full_like(np.asarray(height, dtype=float), g0)
    else:
        if np.any(phi >= g0 * rad):
            raise ValueError('Column extends beyond the potential of the planet')
        z = rad * phi / (g0 * rad - phi)
        gravity = planet.gravity
    pl = np.sqrt(p[1:] * p[:-1])
    tl = 0.5 * (t[1:] + t[:-1])
    xl = 0.5 * (x_h2o[1:] + x_h2o[:-1])
    zl = 0.5 * (z[1:] + z[:-1])
    mass = (p[:-1] - p[1:]) / gravity(zl)
    molar_l = m_dry * (1 - xl) + MOLAR_MASS['H2O'] * xl
    column = mass / molar_l * AVOGADRO * 1e-4           # molecules per cm^2
    return dict(p_pa=p, t_k=t, x_h2o=x_h2o, z_m=z, g=gravity(z), surface_pressure_pa=float(p[0]),
                layer_p_pa=pl, layer_t_k=tl, layer_x_h2o=xl, layer_z_m=zl,
                layer_mass_kg_m2=mass, layer_column_cm2=column, layer_dz_cm=(z[1:] - z[:-1]) * 100.0,
                dry_fractions=dict(air.fractions), planet=planet)


def layer_number_density(col) -> np.ndarray:
    """Total molecules per cm^3 at layer mid-points."""
    return col['layer_p_pa'] / (BOLTZMANN * col['layer_t_k']) * 1e-6


def species_column(col, species) -> np.ndarray:
    """Molecules per cm^2 of a species in each layer."""
    if species == 'H2O':
        return col['layer_column_cm2'] * col['layer_x_h2o']
    return col['layer_column_cm2'] * (1 - col['layer_x_h2o']) * col['dry_fractions'].get(species, 0.0)
