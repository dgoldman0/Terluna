"""Screening estimates behind the conservation gates.

Each function is an order-of-magnitude estimate with its source and assumptions in
the docstring. They date each loss in the transformation and size the heritage
enclosures; domain models supersede them as they are built (surface processes in
geography and climate, frost and heating in atmosphere, structures in engineering,
the appearance from Earth in illumination).
"""
from __future__ import annotations
import math

from shared.constants import GAS_CONSTANT, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

ATM = 101325.0
WATER_DENSITY = 1000.0
AIR_MOLAR_MASS = 0.02897          # kg/mol, Earth-like air
REGOLITH_GRAIN_DENSITY = 3100.0   # kg/m3, lunar soil grains (Carrier et al., Lunar Sourcebook ch. 9)
BASALT_DENSITY = 2700.0


def pressure_at_depth(depth_m, surface_pa, g=MOON_SURFACE_GRAVITY, density=WATER_DENSITY):
    """Absolute pressure under water: surface air pressure plus the water column."""
    return surface_pa + density * g * depth_m


def earth_dive_equivalent(depth_m, surface_pa):
    """Depth of an Earth sea-level dive reaching the same absolute pressure."""
    return (pressure_at_depth(depth_m, surface_pa) - ATM) / (WATER_DENSITY * STANDARD_GRAVITY)


def saltation_threshold(pressure_pa, temperature_k=270.0, g=MOON_SURFACE_GRAVITY,
                        grain_density=REGOLITH_GRAIN_DENSITY, roughness_m=1e-4, height_m=10.0):
    """Lowest threshold friction speed for saltation and the matching wind at 10 m.

    Shao and Lu (2000), J. Geophys. Res. 105(D17), 22437: u*t^2 = A_N ((rho_p/rho_a) g d + gamma/(rho_a d)),
    A_N = 0.0123, gamma = 3e-4 kg/s2, minimised over grain size d. A logarithmic wind profile over a
    smooth bare surface (roughness 0.1 mm) gives the 10 m wind.
    """
    rho_a = pressure_pa * AIR_MOLAR_MASS / (GAS_CONSTANT * temperature_k)
    gamma, a_n = 3e-4, 0.0123
    d = math.sqrt(gamma / (grain_density * g))
    u_star = math.sqrt(a_n * ((grain_density / rho_a) * g * d + gamma / (rho_a * d)))
    return dict(pressure_atm=pressure_pa / ATM, grain_um=d * 1e6, u_star_m_s=u_star,
                wind_10m_m_s=u_star / 0.4 * math.log(height_m / roughness_m))


FROST = {  # triple point (K, Pa) and a sublimation enthalpy (J/mol) for Clausius-Clapeyron
    'N2': (63.15, 12520.0, 6900.0),
    'O2': (54.36, 146.3, 8200.0),
}


def frost_pressure(species, temperature_k):
    """Vapour pressure over the solid (Pa), extrapolated from the triple point with a constant
    sublimation enthalpy. For N2 it gives about 1 Pa at 37 K, matching Pluto's surface pressure over
    its nitrogen ice (New Horizons); for O2 the enthalpy is uncertain by about 15%."""
    t3, p3, h = FROST[species]
    return p3 * math.exp(-h / GAS_CONSTANT * (1.0 / temperature_k - 1.0 / t3))


def warming_depth_m(years, diffusivity_m2_s):
    """Depth a surface temperature change reaches by diffusion, sqrt(kappa t)."""
    return math.sqrt(diffusivity_m2_s * years * 3.15576e7)


def dissolution_years(thickness_m, rate_mol_m2_s, molar_volume_m3=0.120 / 2800.0):
    """Years for dissolution to remove a layer at a surface-normalised rate. The default molar volume is
    basaltic glass per mole of silicon (about 50% SiO2 by mass, density 2,800 kg/m3). Gislason and Oelkers
    (2003) measured about 3e-11 mol Si m-2 s-1 near 25 C and pH 7; field rates run 10^2 to 10^5 times slower
    (White and Brantley 2003)."""
    return thickness_m / (rate_mol_m2_s * molar_volume_m3 * 3.15576e7)


def runoff_scaling(g=MOON_SURFACE_GRAVITY, g_ref=STANDARD_GRAVITY):
    """Ratios of lunar to terrestrial sediment behaviour for the same discharge and slope
    (Darcy-Weisbach flow, Shields entrainment, Meyer-Peter Mueller bedload, Stokes settling)."""
    r = g / g_ref
    return dict(flow_depth=r ** (-1 / 3), bed_shear_stress=r ** (2 / 3), shields_number=r ** (-1 / 3),
                stream_power=r, bedload_far_above_threshold=r ** 0.5 * (r ** (-1 / 3)) ** 1.5,
                settling_speed=r, suspension_ease=r ** (-2 / 3), cohesive_slope_height=1 / r)


def raindrop_energy_ratio(diameter_m, air_density=None, g=MOON_SURFACE_GRAVITY, surface_pa=None,
                          temperature_k=295.0, earth_air_density=1.2):
    """Impact energy per kilogram of rain on the Moon relative to Earth for one drop size (terminal speed
    from a sphere drag law: Schiller-Naumann below Re 1000, Cd 0.44 above)."""
    if air_density is None:
        air_density = (surface_pa or 1.2 * ATM) * AIR_MOLAR_MASS / (GAS_CONSTANT * temperature_k)

    def terminal(gravity, rho_a, mu=1.8e-5):
        v = 1.0
        for _ in range(200):
            re = rho_a * v * diameter_m / mu
            cd = 24 / re * (1 + 0.15 * re ** 0.687) if re < 1000 else 0.44
            v = math.sqrt(4 * gravity * diameter_m * (WATER_DENSITY - rho_a) / (3 * cd * rho_a))
        return v
    return (terminal(g, air_density) / terminal(STANDARD_GRAVITY, earth_air_density)) ** 2


def dome_shell_m(net_pressure_pa, radius_m, modulus_pa, allowable_pa, knockdown=0.2, poisson=0.3):
    """Hemispherical shell thickness for a net external pressure: the larger of membrane stress
    (p R / 2t) and classical buckling 2E/sqrt(3(1-nu^2)) (t/R)^2 reduced by a knockdown factor."""
    classical = 2.0 / math.sqrt(3 * (1 - poisson ** 2))
    buckling = radius_m * math.sqrt(net_pressure_pa / (knockdown * classical * modulus_pa))
    membrane = net_pressure_pa * radius_m / (2 * allowable_pa)
    return max(buckling, membrane)


def rock_equivalent_m(pressure_pa, g=MOON_SURFACE_GRAVITY, density=BASALT_DENSITY):
    """Thickness of rock whose weight equals a pressure: the ballast a floor would need against uplift."""
    return pressure_pa / (density * g)


def earth_tower_equivalent_m(height_m, g=MOON_SURFACE_GRAVITY):
    """Height of an Earth tower with the same self-weight stress as a lunar tower."""
    return height_m * g / STANDARD_GRAVITY


def sediment_return_watts(lowering_m_per_myr, land_share, lift_m, area_m2, rock_density=2600.0, g=MOON_SURFACE_GRAVITY):
    """Ideal power to carry eroded rock back uphill: land lowering times land area, lifted lift_m."""
    mass_per_s = lowering_m_per_myr / 1e6 / 3.15576e7 * land_share * area_m2 * rock_density
    return mass_per_s * g * lift_m, mass_per_s * 3.15576e7


def kardashev_watts(k):
    """Sagan's interpolation K = (log10 P - 6) / 10, P in watts."""
    return 10 ** (10 * k + 6)


def moonlight_factors(bond_albedo_new, phase_integral_new=1.22, geometric_albedo_now=0.12,
                      quarter_share_now=(1 / 15, 1 / 11), quarter_share_new=0.15):
    """Brightness of a terraformed Moon relative to today's, at full and at quarter phase.
    New geometric albedo = Bond albedo / phase integral (Earth-like 1.22, Robinson 2026); today's quarter
    Moon is 1/11-1/15 of full (opposition surge; Krisciunas and Schaefer 1991), an Earth-like planet's
    about 15% (Mallama 2017)."""
    full = bond_albedo_new / phase_integral_new / geometric_albedo_now
    return dict(full=full, quarter=(full * quarter_share_new / quarter_share_now[1],
                                    full * quarter_share_new / quarter_share_now[0]))
