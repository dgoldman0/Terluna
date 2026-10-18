"""A one-dimensional lake through the lunar month: heat, mixing, ice and dissolved oxygen.

The lake is a column of fixed layers over a basin's hypsography (the area at each depth). Each step, of an hour or
two, it takes the sunlight and the infrared from the sky, the air's temperature and humidity and the wind, and:

- heats the column: the sunlight's near-infrared share in the top few decimetres and its visible share by the
  water's own attenuation (Paulson & Simpson 1977), the light reaching the sloping floor heating the water beside it,
  and the infrared, sensible and latent heat at the surface, with the bulk transfer the GCM's own surface layer uses
  (PlaSim's neutral coefficient over Charnock's roughness at the configured gravity, Louis's stability factors and,
  over warm water, Miller et al.'s free convection carried to the configured gravity), paired as PlaSim pairs it when
  the forcing gives the wind's height: the coefficient at the height of the GCM's lowest level, with that level's wind,
  potential temperature and humidity; and the Moon's heat flow through the floor;
- mixes it: statically unstable water overturns, and the energy of that convection and of the wind's stirring
  entrains the water below the mixed layer against its stratification, with the General Lake Model's energy balance
  and its coefficients (Hipsey et al. 2019, eq. 45; C_K 0.2, C_W 0.23, C_T 0.51); below the mixed layer heat and
  dissolved matter diffuse with a diffusivity set by the stratification, a mixing efficiency times a dissipation
  rate over the squared buoyancy frequency (Osborn 1980), up to a configured cap, and across statically unstable
  interfaces at a configured value no smaller than the cap, every buoyancy term with the configured gravity; the
  dissipation comes from the wind (a share of the wind's energy flux at 10 m, P10 = rho_air C10 U10^3, over the mass
  of water below the mixed layer, as Wuest, Piepke & Van Senden 2000 measured it), or is fixed, or the diffusivity
  follows the Lake number (Yeates & Imberger 2003);
- freezes, when the surface reaches freezing, into ice that grows and melts by its own surface balance and
  conduction (a zero-layer ice model in Semtner's manner), shading the water and closing it to the air;
- carries dissolved oxygen: exchange with the air at the saturation the Benson and Krause fit gives (Garcia &
  Gordon 1992, with Benson & Krause's pressure and moisture factor) for the air's own oxygen fraction, at a transfer
  velocity from the surface-renewal model of Tedford et al. (2014) or from Cole & Caraco (1998); gross photosynthesis
  where light reaches; respiration in the euphotic water, the water below and the sediment, limited by the oxygen
  present, at rates the supply of organic matter fixes or with an Arrhenius dependence on temperature from a stated
  reference; the euphotic water's share either at a fixed rate or from a labile pool that the day's production feeds
  and that decays with a stated e-folding time, so respiration runs down through a long night; and the oxygen debt
  of the demand an anoxic layer cannot meet, reduced substances that take oxygen back where it returns;
- tracks the water's age since it last touched the surface, and the dissolved products of respiration in one chosen
  column, which make the water denser and so can hold deep water against an overturn: aerobic respiration's CO2
  (0.284e-3 per g/kg, as Lake Kivu's density is reckoned), and the anaerobic products the debt stands for, its CO2
  and methane's lightening (-1.25e-3 per g/kg) at a configured density per unit of debt, with an optional capped
  term for the dissolved iron and bicarbonate of iron reduction. The solute is carried and mixed in every column
  whether or not it reaches the density.

Concentrations are g/m3 (equal to mg/L); heat is per unit volume at a fixed heat capacity, so the column's heat
changes only by the fluxes through its surface and floor. Buoyancy follows Kell's (1975) density of pure water with a
linear salinity term.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg.lapack import dgtsv

from shared.constants import (BOLTZMANN, ELEMENTARY_CHARGE, EXOPLASIM_CHARNOCK, EXOPLASIM_SEA_ROUGHNESS_FLOOR,
                              EXOPLASIM_VON_KARMAN, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY, STEFAN_BOLTZMANN)

# Water and ice
RHO0 = 998.0                       # reference density for heat capacity and buoyancy, kg/m3
CP_WATER = 4186.0                  # J/kg/K
RHO_CP = RHO0 * CP_WATER
KAPPA_HEAT = 1.4e-7                # molecular thermal diffusivity, m2/s
DIFFUSIVITY_O2 = 2.1e-9            # molecular diffusivity of O2 in water near 20 C, m2/s
SALINITY_DENSITY = 0.8             # kg/m3 per g/kg, the linear salinity term
FREEZING_C = 0.0
RHO_ICE, LATENT_FUSION, K_ICE = 917.0, 3.34e5, 2.2
ICE_ALBEDO, ICE_ATTENUATION = 0.5, 1.5          # bare ice: reflection, and attenuation of what passes (1/m)
WATER_EMISSIVITY = 0.98            # PlaSim's emissivity of open water (radmod.f90)
SUN_AT_PAR = 2.06                  # micromol photons per J of sunlight: 0.45 of it photosynthetically active at 4.57
# Oxygen solubility (USGS Office of Water Quality Technical Memorandum 2011.03, equations 9-14)
MG_PER_ML_O2 = 1.42905             # mg per mL of O2 as a real gas at STP
O2_REFERENCE_FRACTION = 0.20946    # dry air's O2 mole fraction in Benson & Krause's and Garcia & Gordon's fits
GG_A = (2.00907, 3.22014, 4.05010, 4.94457, -0.256847, 3.88767)
GG_B = (-6.24523e-3, -7.37614e-3, -1.03410e-2, -8.17083e-3)
GG_C0 = -4.88682e-7
SCHMIDT_O2_20C = 1800.6 - 120.10 * 20.0 + 3.7818 * 400.0 - 0.047608 * 8000.0   # Wanninkhof (1992), fresh water
RENEWAL_DEPTH_M = 0.15             # depth of the dissipation rate in the surface-renewal model (Tedford et al. 2014)
BOLTZMANN_EV = BOLTZMANN / ELEMENTARY_CHARGE         # eV/K, for respiration's Arrhenius dependence
# Dissolved gases' effect on density, as Lake Kivu's is reckoned (Schmid et al. 2004): kg/m3 of density for each g/m3
# dissolved (the coefficients per g/kg times the water's 1,000 kg/m3)
KIVU_BETA_CO2 = 0.284e-3
KIVU_BETA_CH4 = -1.25e-3
M_CO2, M_CH4 = 44.01, 16.04
# PlaSim's Louis factors (fluxmod.f90: vdiff_b, vdiff_c, vdiff_d)
LOUIS_B = LOUIS_C = LOUIS_D = 5.0
# The free-convection transfer velocity over water in unstable air, 0.0016 m/s per K^(1/3) of the virtual temperature
# difference (Miller et al. 1992, as PlaSim's fluxmod.f90 carries it at Earth's gravity); turbulent free convection
# scales as (g dT)^(1/3), so it is carried to the configured gravity by (g / g_Earth)^(1/3).
FREE_CONVECTION = 0.0016


# --- water ----------------------------------------------------------------------------------------------------------

def density(t_c, salinity=0.0):
    """Density of fresh water, kg/m3: Kell (1975) for pure water at one atmosphere and a linear salinity term."""
    t = np.asarray(t_c, dtype=float)
    pure = ((((((-280.54253e-12 * t + 105.56302e-9) * t - 46.170461e-6) * t - 7.9870401e-3) * t + 16.945176) * t
             + 999.83952) / (1.0 + 16.879850e-3 * t))
    return pure + SALINITY_DENSITY * salinity


def viscosity(t_c):
    """Dynamic viscosity of water, Pa s (Vogel's equation)."""
    return 2.414e-5 * 10.0 ** (247.8 / (np.asarray(t_c, dtype=float) + 273.15 - 140.0))


def schmidt_o2(t_c):
    """Schmidt number of O2 in fresh water: Wanninkhof's (1992) value at 20 C carried to other temperatures by the
    Stokes-Einstein scaling of the diffusivity (D proportional to T over the viscosity), which holds past his fit's
    0-30 C."""
    t = np.asarray(t_c, dtype=float)
    return (SCHMIDT_O2_20C * (viscosity(t) / viscosity(20.0)) ** 2 * (density(20.0) / density(t))
            * (293.15 / (t + 273.15)))


def _rho(t, salinity=0.0):
    """density() for one temperature, in plain floats for the time loop."""
    return ((((((-280.54253e-12 * t + 105.56302e-9) * t - 46.170461e-6) * t - 7.9870401e-3) * t + 16.945176) * t
             + 999.83952) / (1.0 + 16.879850e-3 * t)) + SALINITY_DENSITY * salinity


_TABLE_T = np.linspace(-5.0, 60.0, 65001)              # 0.001 K steps; linear interpolation errs by ~1e-9 kg/m3


def _density_table(salinity):
    """density() on a fine temperature table, for the time loop's arrays."""
    return density(_TABLE_T, salinity)


def _mu(t):
    return 2.414e-5 * 10.0 ** (247.8 / (t + 273.15 - 140.0))


_MU20 = 2.414e-5 * 10.0 ** (247.8 / (20.0 + 273.15 - 140.0))


def _schmidt(t):
    mu = _mu(t)
    return SCHMIDT_O2_20C * (mu / _MU20) ** 2 * (_rho(20.0) / _rho(t)) * (293.15 / (t + 273.15))


def _expansion(t_c):
    """Thermal expansion coefficient of fresh water (1/K) from Kell's density."""
    return -(_rho(t_c + 0.05) - _rho(t_c - 0.05)) / (0.1 * _rho(t_c))


# --- oxygen at saturation -------------------------------------------------------------------------------------------

def vapour_pressure_atm(t_c):
    """Water's vapour pressure (atm), Benson & Krause (1980) as in the USGS memorandum's equation 11."""
    tk = np.asarray(t_c, dtype=float) + 273.15
    return np.exp(11.8571 - 3840.70 / tk - 216961.0 / tk ** 2)


def oxygen_saturation(t_c, pressure_atm=1.0, salinity=0.0, o2_fraction=O2_REFERENCE_FRACTION):
    """Dissolved oxygen in equilibrium with moist air, mg/L.

    Garcia & Gordon's (1992) fit to Benson & Krause's data (mL/L at one atmosphere of moist air, converted at
    1.42905 mg/mL), Benson & Krause's (1984) pressure factor with the water's vapour pressure and the virial term,
    and the dry air's O2 mole fraction relative to Earth's 0.20946. pressure_atm is the total pressure of the moist
    air at the surface. The fits are stated for 0-40 C, 0-40 g/kg and 0.5-1.1 atm; the pressure factor is
    thermodynamic, and at 1.2 atm the next correction it omits, the Poynting term for dissolved O2, is 3e-4.
    """
    t = np.asarray(t_c, dtype=float)
    ts = np.log((298.15 - t) / (273.15 + t))
    a0, a1, a2, a3, a4, a5 = GG_A
    b0, b1, b2, b3 = GG_B
    ln_c = (a0 + ts * (a1 + ts * (a2 + ts * (a3 + ts * (a4 + ts * a5))))
            + salinity * (b0 + ts * (b1 + ts * (b2 + ts * b3))) + GG_C0 * salinity ** 2)
    c0 = MG_PER_ML_O2 * np.exp(ln_c)
    u = vapour_pressure_atm(t)
    theta = 0.000975 - 1.426e-5 * t + 6.436e-8 * t ** 2
    p = np.asarray(pressure_atm, dtype=float)
    fp = (p - u) * (1.0 - theta * p) / ((1.0 - u) * (1.0 - theta))
    return c0 * fp * np.asarray(o2_fraction, dtype=float) / O2_REFERENCE_FRACTION


def _saturation(t, p, salinity, o2_fraction):
    """oxygen_saturation() for one temperature and pressure, in plain floats for the time loop."""
    ts = math.log((298.15 - t) / (273.15 + t))
    a0, a1, a2, a3, a4, a5 = GG_A
    b0, b1, b2, b3 = GG_B
    ln_c = (a0 + ts * (a1 + ts * (a2 + ts * (a3 + ts * (a4 + ts * a5))))
            + salinity * (b0 + ts * (b1 + ts * (b2 + ts * b3))) + GG_C0 * salinity * salinity)
    tk = t + 273.15
    u = math.exp(11.8571 - 3840.70 / tk - 216961.0 / (tk * tk))
    theta = 0.000975 - 1.426e-5 * t + 6.436e-8 * t * t
    return (MG_PER_ML_O2 * math.exp(ln_c) * (p - u) * (1.0 - theta * p) / ((1.0 - u) * (1.0 - theta))
            * o2_fraction / O2_REFERENCE_FRACTION)


# --- air-water exchange ---------------------------------------------------------------------------------------------

def saturation_humidity(t_k, p_pa):
    """Specific humidity at saturation over water (Bolton 1980), as the climate studies take it."""
    e = 611.2 * math.exp(17.67 * (t_k - 273.15) / (t_k - 29.65))
    return 0.622 * e / (p_pa - 0.378 * e)


def neutral_coefficient(u, z, gravity):
    """PlaSim's neutral transfer coefficient over water, (kappa / ln(z/z0 + 1))^2, with z0 by Charnock's relation at
    the given gravity and PlaSim's floor (seamod.f90 and fluxmod.f90)."""
    z0 = EXOPLASIM_SEA_ROUGHNESS_FLOOR
    for _ in range(4):
        cn = (EXOPLASIM_VON_KARMAN / math.log(z / z0 + 1.0)) ** 2
        z0 = max(EXOPLASIM_CHARNOCK * cn * u * u / gravity, EXOPLASIM_SEA_ROUGHNESS_FLOOR)
    return (EXOPLASIM_VON_KARMAN / math.log(z / z0 + 1.0)) ** 2, z0


def reduce_wind(u_level, z_level, z=10.0, gravity=MOON_SURFACE_GRAVITY):
    """The wind at height z under a wind u_level at z_level, on the neutral log profile over water that PlaSim's
    surface layer assumes between the ground and its lowest level."""
    u_level = float(u_level)
    if u_level <= 0.0:
        return 0.0
    cn, z0 = neutral_coefficient(u_level, z_level, gravity)
    u_star = math.sqrt(cn) * u_level
    return u_star / EXOPLASIM_VON_KARMAN * math.log(z / z0 + 1.0)


@dataclass(frozen=True)
class Air:
    """The air over the lake: its dry gas constant and heat capacity, and the reference heights."""
    gas_constant: float                 # J/kg/K of the dry air
    cp: float                           # J/kg/K
    gravity: float = MOON_SURFACE_GRAVITY
    wind_height_m: float = 10.0
    min_wind: float = 0.5               # m/s floor on the wind in the transfer, as PlaSim's zumin floors its own
    free_convection_gravity: float = None   # gravity for the free-convection term; None: the configured gravity


def bulk_fluxes(ts_c, air_k, q_air, wind, p_pa, air: Air, z=None):
    """Sensible and latent heat leaving the water (W/m2), evaporation (kg/m2/s) and the wind stress (N/m2).

    The neutral coefficient at the wind's height z (air.wind_height_m if not given) with Charnock's roughness, and
    PlaSim's factors on the bulk Richardson number between the surface and that height, with the configured gravity:
    Louis's for momentum, and for heat Louis's in stable air and Miller et al.'s free convection over water in unstable
    air. The wind, air temperature and humidity are taken at the same height: PlaSim (fluxmod.f90) pairs its lowest
    level's wind and humidity and that level's potential temperature at the surface pressure (its tas, dtsa) with the
    coefficient at the level's height, znl, about 0.86 km over the Moon's design air.
    """
    tsk = ts_c + 273.15
    u = max(wind, air.min_wind)
    qs = saturation_humidity(tsk, p_pa)
    q_air = min(q_air, saturation_humidity(air_k, p_pa))
    tv_air, tv_s = air_k * (1.0 + 0.608 * q_air), tsk * (1.0 + 0.608 * qs)
    rho_air = p_pa / (air.gas_constant * tv_air)
    z = air.wind_height_m if z is None else float(z)
    cn, z0 = neutral_coefficient(u, z, air.gravity)
    ri = air.gravity * z * (tv_air - tv_s) / (0.5 * (tv_air + tv_s) * u * u)
    if ri <= 0.0:
        denom = 1.0 + 3.0 * LOUIS_C * LOUIS_B * math.sqrt(-ri * (z / z0 + 1.0)) * cn
        fm = 1.0 - 2.0 * LOUIS_B * ri / denom
        g_fc = air.gravity if air.free_convection_gravity is None else air.free_convection_gravity
        w_fc = FREE_CONVECTION * (g_fc / STANDARD_GRAVITY) ** (1.0 / 3.0) * max(tv_s - tv_air, 0.0) ** (1.0 / 3.0)
        fh = (1.0 + (w_fc / (u * cn)) ** 1.25) ** 0.8          # PlaSim's form over water (Miller et al. 1992)
    else:
        root = math.sqrt(1.0 + LOUIS_D * ri)
        fh = 1.0 / (1.0 + 3.0 * LOUIS_B * ri * root)
        fm = 1.0 / (1.0 + 2.0 * LOUIS_B * ri / root)
    ch = cn * fh * u * rho_air
    latent_heat = 2.501e6 - 2370.0 * ts_c
    evaporation = ch * (qs - q_air)
    return ch * air.cp * (tsk - air_k), latent_heat * evaporation, evaporation, rho_air * cn * fm * u * u


def transfer_velocity(model, wind10, u_star_water, cooling_buoyancy, t_c):
    """Gas transfer velocity of O2 (m/s).

    'renewal': the surface-renewal model k = c (nu eps)^(1/4) Sc^(-1/2) with c = 0.5 and the dissipation at 0.15 m of
    Tedford et al. (2014), 0.56 u*^3/(kappa z) + 0.77 B when the surface cools (B the buoyancy flux, which carries the
    gravity) and 0.6 u*^3/(kappa z) when it warms, as Erkkila et al. (2018) and Kiuru et al. (2019) apply it.
    'cole_caraco': k600 = 2.07 + 0.215 U10^1.7 cm/h (Cole & Caraco 1998), carried to O2 with Sc^-2/3 below 3.7 m/s
    and Sc^-1/2 above (Jahne et al. 1987).
    """
    t = float(t_c)
    sc = _schmidt(t)
    if model == 'renewal':
        eps_s = u_star_water ** 3 / (EXOPLASIM_VON_KARMAN * RENEWAL_DEPTH_M)
        eps = 0.56 * eps_s + 0.77 * cooling_buoyancy if cooling_buoyancy > 0.0 else 0.6 * eps_s
        nu = _mu(t) / _rho(t)
        return 0.5 * (nu * eps) ** 0.25 * sc ** -0.5
    if model == 'cole_caraco':
        k600 = (2.07 + 0.215 * wind10 ** 1.7) / 3.6e5
        return k600 * (sc / 600.0) ** (-2.0 / 3.0 if wind10 < 3.7 else -0.5)
    raise ValueError(f'unknown transfer model {model!r}')


# --- the basin ------------------------------------------------------------------------------------------------------

@dataclass
class Basin:
    """Layers between interfaces at depths z (N+1, from 0 at the surface), with the lake's area A at each interface."""
    z: np.ndarray
    area: np.ndarray

    def __post_init__(self):
        self.z = np.asarray(self.z, dtype=float)
        self.area = np.asarray(self.area, dtype=float)
        if self.z[0] != 0.0 or np.any(np.diff(self.z) <= 0) or np.any(np.diff(self.area) > 1e-9 * self.area[0]):
            raise ValueError('interfaces must deepen from 0 and the area must not grow with depth')
        self.dz = np.diff(self.z)
        self.volume = 0.5 * (self.area[:-1] + self.area[1:]) * self.dz
        self.sediment = self.area[:-1] - self.area[1:]
        self.sediment[-1] += self.area[-1]                  # the floor's flat part lies in the bottom layer
        self.centre = 0.5 * (self.z[:-1] + self.z[1:])
        self.depth = float(self.z[-1])

    @property
    def n(self):
        return self.dz.size

    @classmethod
    def pan(cls, depth, layer_m, area_m2):
        n = max(1, int(round(depth / layer_m)))
        return cls(np.linspace(0.0, depth, n + 1), np.full(n + 1, float(area_m2)))

    @classmethod
    def from_hypsography(cls, depth, area_m2, depth_fraction, area_fraction, top_m=0.5, top_depth=20.0,
                         growth=1.08, max_m=8.0):
        """Fine layers to top_depth, then growing by growth to max_m; the area at each interface interpolated in the
        hypsography (the share of the surface area deeper than each fraction of the deepest point)."""
        z = [0.0]
        dz = top_m
        while z[-1] < depth - 1e-9:
            if z[-1] >= top_depth:
                dz = min(max_m, dz * growth)
            z.append(min(depth, z[-1] + dz))
        if len(z) > 2 and z[-1] - z[-2] < 0.3 * dz:
            z.pop(-2)
        z = np.array(z)
        frac = np.interp(z / depth, depth_fraction, area_fraction)
        frac = np.minimum.accumulate(np.maximum(frac, 0.0))
        return cls(z, frac * area_m2)


# --- parameters -----------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Physics:
    gravity: float = MOON_SURFACE_GRAVITY
    salinity: float = 0.0               # g/kg
    albedo: float = 0.069               # PlaSim's open water (seamod.f90 albsea)
    near_ir_share: float = 0.58         # Paulson & Simpson (1977) water type I
    near_ir_scale_m: float = 0.35
    attenuation: float = 0.57           # visible light's attenuation, 1/m
    c_k: float = 0.2                    # GLM's convective overturn efficiency (Hipsey et al. 2019, table 1)
    c_w: float = 0.23                   # wind stirring
    c_t: float = 0.51                   # unsteady turbulence (acceleration of entrained water)
    mixing_efficiency: float = 0.2      # Osborn's (1980) Gamma
    dissipation: float = 1e-10          # turbulent dissipation below the mixed layer, W/kg
    max_diffusivity: float = 1e-4       # m2/s: the cap on stratified water's diffusivity
    unstable_diffusivity: float = 1e-4  # m2/s across statically unstable interfaces below the mixed layer, raised to
                                        # the cap where the cap is higher, so stable water never mixes faster
    deep_mixing: str = 'dissipation'    # the treatment below the mixed layer: 'dissipation', Osborn's rule on the fixed
                                        # dissipation above; 'wind', Osborn's rule on the share wind_share of the
                                        # wind's energy flux P10 spread over the water below the mixed layer;
                                        # 'lake_number', Yeates & Imberger's (2003) diffusivity from the Lake number
    wind_share: float = 0.003           # 'wind': the share of P10 reaching the stratified water as turbulent kinetic
                                        # energy, 0.3 +/- 0.1% (Wuest, Piepke & Van Senden 2000)
    lake_number_constant: float = 200.0     # 'lake_number': K = c K_M (N^2/N_max^2)/L_N (Yeates & Imberger 2003)
    lake_number_z_t_m: float = 10.0     # the Lake number's thermocline: the strongest stratification below this
    heat_flow: float = 0.015            # W/m2 through the floor
    full_mixing: bool = False           # a slab: the whole column mixed every step (a check against a slab ocean)
    solute_density: float = 0.0         # kg/m3 of density per g/m3 of O2 respired aerobically (its CO2)
    solute_column: int = 0              # the oxygen column whose respiration sets the solutes
    debt_density: float = 0.0           # kg/m3 per g/m3 of that column's debt: the anaerobic products it stands for
    iron_density_per_debt: float = 0.0  # kg/m3 per g/m3 of debt from iron reduction's dissolved iron and bicarbonate,
    iron_density_cap: float = 0.0       # up to this, kg/m3 (the dissolved iron a ferruginous lake holds); the cap acts
                                        # on the debt, so a mixture of capped and uncapped water is up to twice as
                                        # dense in iron as the mean of its parts: a bracket carried with the debt


@dataclass(frozen=True)
class Scenario:
    """One oxygen column: what grows and respires, how fast gas crosses the surface, and how warmth acts on demand."""
    name: str
    gpp: float                          # gross production, mean over the cycle, g O2 per m2 of lake surface per day
    attenuation: float                  # its water's visible attenuation for the light photosynthesis uses, 1/m
    sediment: float                     # sediment demand, g O2 per m2 of floor per day
    deep_water: float                   # water-column demand below the euphotic zone, g O2 per m2 per day
    euphotic_share: float = 0.8         # share of production respired in the euphotic water itself
    deep_scale_m: float = 25.0          # thickness below the euphotic zone over which the deep demand acts
    transfer: str = 'renewal'
    activation_ev: float = 0.0          # respiration's Arrhenius energy (0: demand fixed by supply)
    half_saturation: float = 0.5        # g/m3: respiration slows as oxygen falls toward it
    debt_rate: float = 1.0              # per day: reduced substances taking oxygen back where it returns
    light_saturation: float = 150.0     # micromol photons/m2/s: tanh light response
    labile_days: float = None           # e-folding of the euphotic water's labile pool (None: a fixed rate)
    reference_c: float = 20.0           # the temperature at which the Arrhenius dependence leaves the rates as given


@dataclass
class Forcing:
    """Step-by-step forcing, each array one value per step of dt seconds."""
    dt: float
    sw_down: np.ndarray                 # W/m2
    lw_down: np.ndarray                 # W/m2
    air_k: np.ndarray                   # air temperature near the surface, K
    q_air: np.ndarray                   # specific humidity
    wind10: np.ndarray                  # m/s at 10 m (Cole & Caraco's transfer; the bulk exchange's if no level)
    pressure: np.ndarray                # surface pressure, Pa
    o2_fraction: float                  # the dry air's O2 mole fraction
    wind_level: np.ndarray = None       # m/s at level_height, where the bulk exchange pairs wind, air and humidity
    level_height: np.ndarray = None     # m: that height, the GCM's lowest level (None: the air's wind_height_m)

    def __len__(self):
        return int(np.asarray(self.sw_down).size)


# --- the model ------------------------------------------------------------------------------------------------------

def light_shape(basin: Basin, physics: Physics):
    """Share of the absorbed sunlight that heats each layer; the light reaching the floor heats the water over it."""
    z, a = basin.z, basin.area
    level = (physics.near_ir_share * np.exp(-z / physics.near_ir_scale_m)
             + (1.0 - physics.near_ir_share) * np.exp(-physics.attenuation * z))
    through = level * a
    share = through[:-1] - through[1:]
    share[-1] = through[-2]                              # whatever reaches the last layer stays in it
    return share / a[0]


def layer_light(basin: Basin, attenuation):
    """Each layer's mean visible light as a share of the light just below the surface."""
    z = basin.z
    return (np.exp(-attenuation * z[:-1]) - np.exp(-attenuation * z[1:])) / (attenuation * basin.dz)


def demand_profile(basin: Basin, s: Scenario, euphotic=None):
    """Respiration's volumetric rate in each layer (g/m3/day) at the reference temperature.

    The euphotic water (to the 1% light depth) respires the euphotic share of production; the deep demand acts over
    deep_scale_m below it per m2 of the euphotic zone's base, its part that finds no water there falling on the
    sediment; the sediment demand acts on every square metre of floor. With a labile pool (s.labile_days) the euphotic
    share comes from the pool as the run goes, and the profile leaves it out unless euphotic is True."""
    z_eu = min(4.6 / s.attenuation, basin.depth)
    top, bottom, centre = basin.z[:-1], basin.z[1:], basin.centre
    a0 = basin.area[0]
    rate = np.zeros(basin.n)
    include = (s.labile_days is None) if euphotic is None else bool(euphotic)
    euphotic = np.clip((np.minimum(bottom, z_eu) - top) / basin.dz, 0.0, 1.0)
    v_eu = float((euphotic * basin.volume).sum())
    if include:
        rate += euphotic * s.euphotic_share * s.gpp * a0 / v_eu
    lower = min(z_eu + s.deep_scale_m, basin.depth)
    deep = np.clip((np.minimum(bottom, lower) - np.maximum(top, z_eu)) / basin.dz, 0.0, 1.0)
    v_deep = float((deep * basin.volume).sum())
    a_eu = float(np.interp(z_eu, basin.z, basin.area))
    volumetric = s.deep_water / s.deep_scale_m
    rate += deep * volumetric
    missing = max(0.0, s.deep_water * a_eu - volumetric * v_deep)          # g/day with no water to act in
    floor_in_range = basin.sediment * ((centre >= z_eu) & (centre < lower))
    if missing > 0.0:
        weight = floor_in_range if floor_in_range.sum() > 0 else basin.sediment * (np.arange(basin.n) == basin.n - 1)
        rate += missing * weight / weight.sum() / basin.volume
    rate += s.sediment * basin.sediment / basin.volume
    return rate


def anaerobic_density(methane_share):
    """kg/m3 of density for each g/m3 of debt (oxygen equivalents of reduced products): its methane, a share of its
    oxygen equivalents at 64 g a mole, and its CO2, a mole with each mole of methane and a mole for each 32 g of the
    rest, at Lake Kivu's coefficients. When the debt is oxidised its products become aerobic respiration's CO2."""
    co2 = M_CO2 * (methane_share / 64.0 + (1.0 - methane_share) / 32.0)
    ch4 = M_CH4 * methane_share / 64.0
    return KIVU_BETA_CO2 * co2 + KIVU_BETA_CH4 * ch4


DEEP_MIXING = ('dissipation', 'wind', 'lake_number')   # the deep-mixing treatments implemented


def drag_coefficient_10m(u10):
    """The wind's drag coefficient at 10 m over a lake, C10 = tau/(rho_air U10^2), as Wuest & Lorke (2003) give it
    for light and moderate winds and Simstrat takes it (Gaudard et al. 2017, eq. 5)."""
    u = np.maximum(np.asarray(u10, dtype=float), 0.0)
    return np.where(u <= 0.1, 6.215e-2,
                    np.where(u <= 3.85, 4.4e-3 * np.maximum(u, 0.1) ** -1.15,
                             -7.12e-7 * u * u + 7.387e-5 * u + 6.605e-4))


def wind_energy_flux(u10, rho_air):
    """The wind's energy flux at 10 m, P10 = rho_air C10 U10^3 (W/m2), the reference flux of Wuest, Piepke & Van
    Senden (2000)."""
    u = np.asarray(u10, dtype=float)
    return rho_air * drag_coefficient_10m(u) * u ** 3


def wind_dissipation(p10, hypolimnion_kg_m2, physics: Physics):
    """The mean dissipation below the mixed layer (W/kg) from the wind: the share of P10 that reaches the stratified
    water, spread over its mass per square metre of surface."""
    return physics.wind_share * p10 / hypolimnion_kg_m2 if hypolimnion_kg_m2 > 0.0 else 0.0


def diffusivity(n2, physics: Physics, dissipation=None):
    """The turbulent diffusivity below the mixed layer (m2/s): Osborn's mixing efficiency times the dissipation over
    the squared buoyancy frequency, capped, and across statically unstable interfaces the unstable value or the cap,
    whichever is larger, so that no stable interface mixes faster than an unstable one. The dissipation is the fixed
    one, or the wind's, given."""
    if physics.deep_mixing not in DEEP_MIXING:
        raise NotImplementedError(f'deep mixing {physics.deep_mixing!r}: only {DEEP_MIXING} is implemented')
    eps = physics.dissipation if dissipation is None else float(dissipation)
    n2 = np.asarray(n2, dtype=float)
    k = physics.mixing_efficiency * eps / np.maximum(n2, 1e-30)
    unstable = max(physics.unstable_diffusivity, physics.max_diffusivity)
    return np.where(n2 > 0.0, np.minimum(k, physics.max_diffusivity), unstable)


def lake_number_diffusivity(n2, physics: Physics, lake_number):
    """Yeates & Imberger's (2003) diffusivity below the mixed layer (m2/s), calibrated on Lake Kinneret and two
    reservoirs: c K_M (N^2/N_max^2)/L_N with c = 200 and K_M the molecular diffusivity of heat, so it grows as the
    Lake number falls and is largest at the thermocline; capped, with the unstable rule of diffusivity()."""
    n2 = np.asarray(n2, dtype=float)
    n2_max = float(n2.max()) if n2.size else 0.0
    if n2_max > 0.0 and lake_number > 0.0 and math.isfinite(lake_number):
        k = physics.lake_number_constant * KAPPA_HEAT * np.maximum(n2, 0.0) / n2_max / lake_number
    else:
        k = np.zeros_like(n2)
    unstable = max(physics.unstable_diffusivity, physics.max_diffusivity)
    return np.where(n2 > 0.0, np.minimum(k, physics.max_diffusivity), unstable)


@dataclass
class State:
    t: np.ndarray                       # temperature, C
    oxygen: np.ndarray                  # (N, M) g/m3
    debt: np.ndarray                    # (N, M) g O2 per m3 owed by reduced substances
    age: np.ndarray                     # days since surface contact
    ice: float = 0.0                    # thickness, m
    ice_surface_c: float = 0.0
    mixed_depth: float = 0.0
    solute: np.ndarray = None           # g/m3 of O2 respired, net of production, since the water left the air
    pool: np.ndarray = None             # (N, M) g O2 per m3 of labile organic matter awaiting respiration

    def __post_init__(self):
        if self.solute is None:
            self.solute = np.zeros_like(self.t)
        if self.pool is None:
            self.pool = np.zeros_like(self.oxygen)

    def copy(self):
        return State(self.t.copy(), self.oxygen.copy(), self.debt.copy(), self.age.copy(), self.ice,
                     self.ice_surface_c, self.mixed_depth, self.solute.copy(), self.pool.copy())


def initial_state(basin: Basin, scenarios, t_c, oxygen=None):
    n, m = basin.n, len(scenarios)
    t = np.full(n, float(t_c)) if np.isscalar(t_c) else np.asarray(t_c, dtype=float).copy()
    ox = np.full((n, m), float(oxygen if oxygen is not None else 8.0)) if oxygen is None or np.isscalar(oxygen) \
        else np.asarray(oxygen, dtype=float).copy()
    return State(t, ox, np.zeros((n, m)), np.zeros(n))


class Lake:
    def __init__(self, basin: Basin, physics: Physics, scenarios, air: Air):
        self.basin, self.physics, self.scenarios, self.air = basin, physics, list(scenarios), air
        self.heat_share = light_shape(basin, physics)
        self.light = np.stack([layer_light(basin, s.attenuation) for s in self.scenarios], axis=1)
        self.demand = np.stack([demand_profile(basin, s) for s in self.scenarios], axis=1) / 86400.0
        self.activation = np.array([s.activation_ev for s in self.scenarios]) / BOLTZMANN_EV
        self.kinetic = bool(np.any(self.activation != 0.0))
        self.reference_k = np.array([s.reference_c for s in self.scenarios]) + 273.15
        self.labile_k = np.array([0.0 if s.labile_days is None else 1.0 / (s.labile_days * 86400.0)
                                  for s in self.scenarios])
        self.feed_share = np.array([0.0 if s.labile_days is None else s.euphotic_share for s in self.scenarios])
        self.pooled = bool(np.any(self.labile_k > 0.0))
        self.debt_dense = bool(physics.debt_density or physics.iron_density_per_debt)
        self.half_sat = np.array([s.half_saturation for s in self.scenarios])
        self.debt_rate = np.array([s.debt_rate for s in self.scenarios]) / 86400.0
        self.ik = np.array([s.light_saturation for s in self.scenarios])
        self.models = sorted({s.transfer for s in self.scenarios})
        self.model_index = np.array([self.models.index(s.transfer) for s in self.scenarios])
        self.production_scale = None
        self._cum_volume = np.cumsum(basin.volume)
        self._centre_of_volume = float((basin.centre * basin.volume).sum() / basin.volume.sum())
        self._rho_table = _density_table(physics.salinity)

    # production is scaled so that its mean over the forcing is each scenario's areal rate
    def calibrate_production(self, forcing: Forcing):
        b = self.basin
        par0 = np.asarray(forcing.sw_down, dtype=float) * (1.0 - self.physics.albedo) * SUN_AT_PAR
        levels = np.linspace(0.0, max(par0.max(), 1.0), 400)
        scale = np.zeros(len(self.scenarios))
        for m, s in enumerate(self.scenarios):
            g = np.tanh(levels[:, None] * self.light[None, :, m] / s.light_saturation) @ b.volume
            mean = float(np.interp(par0, levels, g).mean())
            scale[m] = s.gpp / 86400.0 * b.area[0] / mean if mean > 0 else 0.0
        self.production_scale = scale
        return scale

    def run(self, forcing: Forcing, state: State, record_every=1, profiles_at=()):
        b, ph, air = self.basin, self.physics, self.air
        n_layers, m_cols = b.n, len(self.scenarios)
        if self.production_scale is None:
            self.calibrate_production(forcing)
        dt = float(forcing.dt)
        g = ph.gravity
        a0 = b.area[0]
        v, dz, zc, ztop, zbot = b.volume, b.dz, b.centre, b.z[:-1], b.z[1:]
        v0 = v[0]
        cum_v = np.cumsum(v)
        heat_share = self.heat_share
        geothermal = ph.heat_flow * b.sediment / (RHO_CP * v) * dt
        heat_vec = heat_share * a0 / (RHO_CP * v) * dt          # K per (W/m2 of absorbed sunlight) per step
        surface_k = a0 / (RHO_CP * v0) * dt
        rho_table = self._rho_table
        exchange_ad = b.area[1:-1] / np.diff(b.centre) * dt
        # interface coefficients: area / distance between layer centres
        a_int = b.area[1:-1]
        d_int = np.diff(zc)
        k_mol_t, k_mol_x = KAPPA_HEAT, DIFFUSIVITY_O2
        sw, lw, ta, qa, u10, pa = (np.asarray(x, dtype=float) for x in (forcing.sw_down, forcing.lw_down,
                                                                         forcing.air_k, forcing.q_air,
                                                                         forcing.wind10, forcing.pressure))
        p_atm = pa / 101325.0
        uz = u10 if forcing.wind_level is None else np.asarray(forcing.wind_level, dtype=float)
        zl = None if forcing.level_height is None else np.asarray(forcing.level_height, dtype=float)
        o2_fraction = forcing.o2_fraction
        steps = len(forcing)
        rec_n = (steps + record_every - 1) // record_every
        out = dict(surface_t=np.zeros(rec_n), bottom_t=np.zeros(rec_n), mixed_depth=np.zeros(rec_n),
                   ice=np.zeros(rec_n), surface_o2=np.zeros((rec_n, m_cols)), bottom_o2=np.zeros((rec_n, m_cols)),
                   deep_o2=np.zeros((rec_n, m_cols)), saturation=np.zeros(rec_n), bottom_age=np.zeros(rec_n),
                   bottom_debt=np.zeros((rec_n, m_cols)), bottom_solute=np.zeros(rec_n), u_star=np.zeros(rec_n),
                   dissipation=np.zeros(rec_n), interfaces_below=0, interfaces_capped=0, interfaces_held_by_cap=0,
                   interfaces_unstable=0,
                   heat_in=0.0,
                   oxygen_in=np.zeros(m_cols),
                   produced=np.zeros(m_cols), respired=np.zeros(m_cols), mean_o2=np.zeros((rec_n, m_cols)),
                   profiles={})
        deep_cut = float(cum_v[-1] * 0.9)
        deep_mask = (cum_v - 0.5 * v) >= deep_cut           # the deepest tenth of the volume, by layer centres
        if not deep_mask.any():
            deep_mask[-1] = True
        v_deep = v[deep_mask]
        t, ox, debt, age, sol, pool = state.t, state.oxygen, state.debt, state.age, state.solute, state.pool
        beta, col = ph.solute_density, ph.solute_column
        ice, ice_ts, h_prev = state.ice, state.ice_surface_c, state.mixed_depth
        prod_scale = self.production_scale[None, :]
        light = self.light
        sigma = STEFAN_BOLTZMANN
        profile_set = set(int(i) for i in profiles_at)
        zero_columns = np.zeros((n_layers, m_cols))
        pooled = self.pooled
        stack = np.zeros((n_layers, (3 if pooled else 2) * m_cols + 2))
        n2_factor = g / (RHO0 * d_int)
        mixing = ph.deep_mixing
        if mixing not in DEEP_MIXING:
            raise NotImplementedError(f'deep mixing {mixing!r}: only {DEEP_MIXING} is implemented')
        cap_reached = ph.max_diffusivity * (1.0 - 1e-12)
        cap_slower = ph.max_diffusivity * dt / (d_int * d_int) < 1.0      # the cap mixes less than a step can
        labile_k, feed_share = self.labile_k[None, :], self.feed_share[None, :]
        reference_k = self.reference_k[None, :]
        w_deep = np.where(deep_mask, v, 0.0) / v[deep_mask].sum()
        w_mean = v / cum_v[-1]
        half_sat, debt_rate, ik = self.half_sat[None, :], self.debt_rate[None, :], self.ik[None, :]
        light_ik = light / ik
        for n in range(steps):
            # --- surface heat ---------------------------------------------------------------------------------
            t0 = float(t[0])
            z_n = None if zl is None else zl[n]
            sens, lat_heat, evap, tau = bulk_fluxes(t0 if ice <= 0.0 else ice_ts, ta[n], qa[n], uz[n], pa[n], air,
                                                    z_n)
            sw_in = sw[n] * (1.0 - ph.albedo)
            if ice > 0.0:
                sw_water, open_water = self._ice_step_light(sw[n], ice), False
            else:
                sw_water, open_water = sw_in, True
            if open_water:
                surface = WATER_EMISSIVITY * (lw[n] - sigma * (t0 + 273.15) ** 4) - sens - lat_heat
                t += sw_water * heat_vec + geothermal
                t[0] += surface * surface_k
                out['heat_in'] += (sw_water + surface) * a0 * dt
                if t[0] < FREEZING_C:                      # freezing: the deficit becomes ice
                    deficit = (FREEZING_C - t[0]) * RHO_CP * v0
                    ice = deficit / (RHO_ICE * LATENT_FUSION * a0)
                    t[0] = FREEZING_C
                    ice_ts = FREEZING_C
            else:
                ice, ice_ts, q_atm, to_water = self._ice_step(ice, ice_ts, sw[n], lw[n], ta[n], qa[n], uz[n],
                                                              pa[n], t, dt, z_n)
                t[0] += to_water * surface_k
                t += sw_water * heat_vec + geothermal
                out['heat_in'] += (sw_water + q_atm) * a0 * dt
            out['heat_in'] += ph.heat_flow * a0 * dt
            # --- mixing ---------------------------------------------------------------------------------------
            u_star_w = math.sqrt(tau / RHO0) if open_water else 0.0
            if ph.full_mixing:
                k_full, frac, h_prev = n_layers, 0.0, b.depth
                t[:] = float((t * v).sum() / cum_v[-1])
                ox[:] = (ox * v[:, None]).sum(0) / cum_v[-1]
                debt[:] = (debt * v[:, None]).sum(0) / cum_v[-1]
                pool[:] = (pool * v[:, None]).sum(0) / cum_v[-1]
                age[:] = 0.0
                sol[:] = 0.0
            else:
                k_full, frac, h_prev = self._mix(t, ox, debt, age, u_star_w, dt, h_prev, g, sol, pool)
            # --- diffusion below the mixed layer --------------------------------------------------------------
            rho = np.interp(t, _TABLE_T, rho_table)
            if beta:
                rho = rho + beta * sol
            if self.debt_dense:
                rho = rho + self._debt_density(debt[:, col])
            n2 = np.diff(rho) * n2_factor
            if mixing == 'wind':                           # a share of the wind's energy over the water below
                below_kg_m2 = RHO0 * float(cum_v[-1] - cum_v[k_full - 1]) / a0 if k_full < n_layers else 0.0
                p10 = float(wind_energy_flux(u10[n], pa[n] / (air.gas_constant * ta[n]))) if open_water else 0.0
                eps_h = wind_dissipation(p10, below_kg_m2, ph)
                k_turb = diffusivity(n2, ph, eps_h)
            elif mixing == 'lake_number':
                eps_h = math.nan
                k_turb = lake_number_diffusivity(n2, ph, self.lake_number(rho, u_star_w * u_star_w, g))
            else:
                eps_h = ph.dissipation
                k_turb = diffusivity(n2, ph)
            first = max(k_full - 1, 0)                     # the interfaces from the mixed layer's base down
            if first < n_layers - 1:
                below = n2[first:]
                out['interfaces_below'] += below.size
                out['interfaces_unstable'] += int((below <= 0.0).sum())
                capped = (below > 0.0) & (k_turb[first:] >= cap_reached)
                out['interfaces_capped'] += int(capped.sum())
                # held: the cap slower than the step's own mixing across the interface, so it limits the exchange
                out['interfaces_held_by_cap'] += int((capped & cap_slower[first:]).sum())
            _diffuse(t, (k_turb + k_mol_t) * exchange_ad, v)
            stack[:, :m_cols] = ox
            stack[:, m_cols:2 * m_cols] = debt
            if pooled:
                stack[:, 2 * m_cols:3 * m_cols] = pool
            stack[:, -2] = age
            stack[:, -1] = sol
            _diffuse(stack, (k_turb + k_mol_x) * exchange_ad, v)
            ox[:], debt[:], age[:], sol[:] = stack[:, :m_cols], stack[:, m_cols:2 * m_cols], stack[:, -2], stack[:, -1]
            if pooled:
                pool[:] = stack[:, 2 * m_cols:3 * m_cols]
            # --- biology --------------------------------------------------------------------------------------
            par0 = sw_water * SUN_AT_PAR                   # photosynthetically active light below the surface
            if par0 > 0.0:
                production = prod_scale * np.tanh(par0 * light_ik)     # g/m3/s
            else:
                production = zero_columns
            demand = self.demand * np.exp(self.activation * (1.0 / reference_k - 1.0 / (t[:, None] + 273.15))) \
                if self.kinetic else self.demand
            if pooled:                                     # the day's production feeds the labile pool, which decays
                pool[:] = (pool + dt * feed_share * production) / (1.0 + dt * labile_k)
                demand = demand + labile_k * pool
            gained = ox + dt * production
            ox_new = gained / (1.0 + dt * demand / (half_sat + np.maximum(ox, 0.0)))
            used = gained - ox_new
            debt += np.maximum(dt * demand - used, 0.0)
            oxidised = np.minimum(debt * (1.0 - 1.0 / (1.0 + dt * debt_rate * ox_new / (half_sat + ox_new))), ox_new)
            debt -= oxidised
            ox_new -= oxidised
            out['produced'] += dt * (v @ production)
            out['respired'] += v @ (used + oxidised)
            sol += used[:, col] + oxidised[:, col] - dt * production[:, col]   # carried whether or not it is dense
            ox[:] = ox_new
            # --- exchange with the air ------------------------------------------------------------------------
            saturation = _saturation(float(t[0]), p_atm[n], ph.salinity, o2_fraction)
            if open_water:
                cooling = sens + lat_heat - WATER_EMISSIVITY * (lw[n] - sigma * (t[0] + 273.15) ** 4) \
                    - sw_water * float(heat_share[0])
                buoyancy = g * _expansion(float(t[0])) * max(cooling, 0.0) / RHO_CP
                kk = [transfer_velocity(mdl, u10[n], u_star_w, buoyancy, t[0]) for mdl in self.models]
                kcol = np.array(kk)[self.model_index]
                kk_full = max(1, k_full)
                v_ml = cum_v[kk_full - 1]
                relax = 1.0 - np.exp(-kcol * a0 * dt / v_ml)
                change = (saturation - ox[:kk_full]) * relax[None, :]
                ox[:kk_full] += change
                out['oxygen_in'] += v[:kk_full] @ change
                sol[:kk_full] *= 1.0 - relax[col]              # the products leave as CO2 at the gas's own pace
            # --- age ------------------------------------------------------------------------------------------
            age += dt / 86400.0
            if open_water:
                age[:max(1, k_full)] = 0.0
            # --- record ---------------------------------------------------------------------------------------
            if n % record_every == 0:
                r = n // record_every
                out['surface_t'][r] = t[0]
                out['bottom_t'][r] = t[-1]
                out['mixed_depth'][r] = h_prev
                out['ice'][r] = ice
                out['surface_o2'][r] = ox[0]
                out['bottom_o2'][r] = ox[-1]
                out['deep_o2'][r] = w_deep @ ox
                out['mean_o2'][r] = w_mean @ ox
                out['saturation'][r] = saturation
                out['bottom_age'][r] = age[-1]
                out['bottom_debt'][r] = debt[-1]
                out['bottom_solute'][r] = sol[-1]
                out['u_star'][r] = u_star_w
                out['dissipation'][r] = eps_h
            if n in profile_set:
                out['profiles'][n] = dict(t=t.copy(), oxygen=ox.copy(), debt=debt.copy(), age=age.copy(),
                                          solute=sol.copy(), pool=pool.copy())
        state.ice, state.ice_surface_c, state.mixed_depth = ice, ice_ts, h_prev
        return out

    # --- pieces ---------------------------------------------------------------------------------------------------
    def lake_number(self, rho, u_star2, g):
        """The column's Lake number (Imberger & Patterson 1990, in the form lake_oxygen reads it): L_N = g S z_T /
        (rho0 u*^2 A0^(1/2) z_cv), S the Schmidt work per unit g and area, z_T the strongest stratification below
        physics.lake_number_z_t_m, z_cv the centre of volume; infinite in calm water."""
        if u_star2 <= 0.0:
            return math.inf
        b = self.basin
        v = b.volume
        zcv = self._centre_of_volume
        rbar = float((rho * v).sum() / v.sum())
        work = float(((b.centre - zcv) * (rho - rbar) * v).sum()) / b.area[0]
        n2 = np.diff(rho) / np.diff(b.centre)
        below = b.centre[1:] > self.physics.lake_number_z_t_m
        zt = float(b.centre[1:][below][np.argmax(n2[below])]) if below.any() else float(b.centre[-1])
        return g * work * zt / (RHO0 * u_star2 * math.sqrt(b.area[0]) * zcv)

    def _debt_density(self, d):
        """The density the debt of the coupled column adds (kg/m3): its anaerobic products, and iron reduction's
        dissolved iron and bicarbonate up to the cap. The cap acts on the debt given, so where the mixed layer joins
        capped and uncapped water the mixture's iron follows its mean debt and can reach twice the mean of its parts'
        iron: the iron term is a bracket on the density, carried with the debt."""
        ph = self.physics
        out = ph.debt_density * d
        if ph.iron_density_per_debt:
            out = out + np.minimum(ph.iron_density_per_debt * d, ph.iron_density_cap)
        return out

    def _mix(self, t, ox, debt, age, u_star, dt, h_prev, g, sol=None, pool=None):
        """Convective overturn, then wind and convective energy entraining the water below (GLM's eq. 45), the last
        layer partly; returns the fully mixed layers, the share of the next and the mixed layer's depth.

        Entraining layer by layer from the surface, the mixed layer just above layer j is the mix of every layer above
        it, so each layer's cost follows from cumulative means and the energy reaches as far as their running sum."""
        b, ph = self.basin, self.physics
        v, dz, zc, zbot, ztop, cv = b.volume, b.dz, b.centre, b.z[1:], b.z[:-1], self._cum_volume
        n = b.n
        table = self._rho_table
        beta = ph.solute_density if sol is not None else 0.0
        rho = np.interp(t, _TABLE_T, table)
        heat = np.cumsum(t * v)
        tbar = heat / cv
        rbar = np.interp(tbar, _TABLE_T, table)
        if beta:
            rho = rho + beta * sol
            rbar = rbar + beta * np.cumsum(sol * v) / cv
        if self.debt_dense:                           # the mixture's debt sets the mixture's density
            dcol = debt[:, ph.solute_column]
            rho = rho + self._debt_density(dcol)
            rbar = rbar + self._debt_density(np.cumsum(dcol * v) / cv)
        stable = rbar[:-1] < rho[1:]
        k = int(np.argmax(stable)) + 1 if stable.any() else n
        released = g * float(((rbar[k - 1] - rho[:k]) * v[:k] * zc[:k]).sum())
        w3 = max(released, 0.0) / (RHO0 * b.area[0] * dt)
        u3 = u_star ** 3
        energy = (ph.c_k * w3 + ph.c_w * u3) * dt
        frac = 0.0
        j = k
        if k < n and energy > 0.0:
            # the velocity scale as the paper's eq. 45 writes it; GLM's code (glm_mixer.c) takes w*^3 + (C_W/C_K) u*^3,
            # which makes q^2 2.9 times larger under wind alone and changes the deepening by under 1% (a check run
            # against this code's tests' stratifications)
            q2 = (w3 + ph.c_w * u3) ** (2.0 / 3.0)
            gp = np.maximum(g * (rho[k:] - rbar[k - 1:-1]) / RHO0, 0.0)
            new = np.clip((zbot[k:] - np.maximum(ztop[k:], h_prev)) / dz[k:], 0.0, 1.0)
            cost = (gp * zbot[k - 1:-1] + ph.c_t * q2 * new) * dz[k:]
            spent = np.cumsum(cost)
            j = k + int(np.searchsorted(spent, energy, side='right'))
            if j < n:
                left = energy - (float(spent[j - k - 1]) if j > k else 0.0)
                frac = left / float(cost[j - k]) if cost[j - k] > 0.0 else 1.0
                if frac >= 1.0:
                    j, frac = j + 1, 0.0
        k_full = j
        v_ml = float(cv[k_full - 1])
        if frac > 0.0 and k_full < n:
            vf = frac * v[k_full]
            w = v[:k_full]
            t_new = (float(heat[k_full - 1]) + vf * float(t[k_full])) / (v_ml + vf)
            t[k_full] = frac * t_new + (1.0 - frac) * t[k_full]
            t[:k_full] = t_new
            for x in self._tracers(ox, debt, sol, pool):
                mean = (w @ x[:k_full] + vf * x[k_full]) / (v_ml + vf)
                x[k_full] = frac * mean + (1.0 - frac) * x[k_full]
                x[:k_full] = mean
            mean_age = (float(w @ age[:k_full]) + vf * float(age[k_full])) / (v_ml + vf)
            age[k_full] = frac * mean_age + (1.0 - frac) * age[k_full]
            age[:k_full] = mean_age
            depth = float(ztop[k_full] + frac * dz[k_full])
        else:
            w = v[:k_full]
            t[:k_full] = float(heat[k_full - 1]) / v_ml
            for x in self._tracers(ox, debt, sol, pool):
                x[:k_full] = (w @ x[:k_full]) / v_ml
            age[:k_full] = float(w @ age[:k_full]) / v_ml
            depth = float(zbot[k_full - 1])
        return k_full, frac, depth

    def _tracers(self, ox, debt, sol, pool):
        """The dissolved and suspended matter the mixing carries with the water."""
        out = [ox, debt]
        if sol is not None:
            out.append(sol)
        if pool is not None and self.pooled:
            out.append(pool)
        return out

    def _ice_step_light(self, sw_down, ice):
        return sw_down * (1.0 - ICE_ALBEDO) * math.exp(-ICE_ATTENUATION * ice)

    def _ice_step(self, ice, ts_ice, sw, lw, ta, qa, u10, pa, t, dt, z=None):
        """A zero-layer ice step: the ice surface's temperature from its balance with conduction, growth or melt at
        the base from the conduction and the top water's heat above freezing, melt at the top once the surface
        reaches freezing. Returns the thickness, the surface temperature, the heat the ice takes from the air
        (W/m2, the column's gain) and the heat the top layer gains from the ice (W/m2)."""
        absorbed = sw * (1.0 - ICE_ALBEDO) * (1.0 - math.exp(-ICE_ATTENUATION * ice))
        h = max(ice, 1e-3)

        def from_air(ts):
            s, lh, _, _ = bulk_fluxes(ts, ta, qa, u10, pa, self.air, z)
            return absorbed + WATER_EMISSIVITY * (lw - STEFAN_BOLTZMANN * (ts + 273.15) ** 4) - s - lh

        if from_air(FREEZING_C) > 0.0:                 # the surface reaches freezing and melts at the top
            ts = FREEZING_C
            q_atm = from_air(FREEZING_C)
            top_melt = q_atm * dt / (RHO_ICE * LATENT_FUSION)
        else:
            lo, hi = -80.0, FREEZING_C
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                if from_air(mid) + K_ICE * (FREEZING_C - mid) / h > 0.0:
                    lo = mid
                else:
                    hi = mid
            ts = 0.5 * (lo + hi)
            q_atm = -K_ICE * (FREEZING_C - ts) / h      # what the air takes is what conduction brings up
            top_melt = 0.0
        conduction = K_ICE * (FREEZING_C - ts) / h
        v0, a0 = self.basin.volume[0], self.basin.area[0]
        excess = max(float(t[0]) - FREEZING_C, 0.0) * RHO_CP * v0 / (a0 * dt)
        new = ice + (conduction - excess) * dt / (RHO_ICE * LATENT_FUSION) - top_melt
        to_water = -excess
        if new <= 0.0:                                 # it melts away, and the heat left over warms the water
            to_water += -new * RHO_ICE * LATENT_FUSION / dt
            new = 0.0
        return new, ts, q_atm, to_water


def _diffuse(x, c, v):
    """Implicit, conservative diffusion of x (a column or columns) through the interfaces, where c holds the volume
    exchanged per unit difference over the step at each interface (LAPACK's dgtsv)."""
    if v.size < 2:
        return
    diag = v.copy()
    diag[:-1] += c
    diag[1:] += c
    rhs = x * (v[:, None] if x.ndim == 2 else v)
    _, _, _, sol, info = dgtsv(-c, diag, -c, rhs, overwrite_b=1)
    if info != 0:
        raise RuntimeError(f'dgtsv failed ({info})')
    x[...] = sol


def heat_content(basin: Basin, state: State):
    """Heat of the column relative to ice at freezing (J): the water's sensible heat with the ice's latent heat
    taken off."""
    return float(RHO_CP * (state.t * basin.volume).sum() - RHO_ICE * LATENT_FUSION * state.ice * basin.area[0])


def oxygen_content(basin: Basin, state: State):
    return (state.oxygen * basin.volume[:, None]).sum(0)
