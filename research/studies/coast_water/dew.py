"""The surface's water balance through the Moon's day and night: the plants' demand by Penman-Monteith, dew on a canopy,
on an exposed leaf and on bare soil from the night's energy balance, and fog where the air at the ground saturates.

Every function takes arrays of the air at about 2 m (K), its specific humidity (kg/kg), the surface pressure (Pa), the
sky's downward longwave (W/m2) and the wind at 2 m (m/s), so the same code runs on the GCM's means, on the CM1 rings'
nights and on an Earth night for the check.

- **Demand.** FAO-56's Penman-Monteith (Allen et al. 1998, equation 3) for the short reference surface, its
  aerodynamic resistance 208/u2 (s/m) with u2 held at 0.5 m/s or more, its albedo of 0.23 and its emission at the air's
  temperature, with the hourly convention for periods of day and of night: surface resistance 50 s/m by day and 200 by
  night (ASCE-EWRI 2005) and the soil taking 0.1 of the net radiation by day and giving back 0.5 of the loss by night
  (FAO-56 equations 45 and 46). The psychrometric constant and the air's density follow the local pressure. Night
  condensation is left to the dew balance, so the demand is held at zero or more.
- **Dew on a canopy.** The canopy as one surface, grass-like: it loses ε(L↓ - σTc^4) to the sky, draws half that loss
  back from the soil at night as the demand's convention does, and takes heat and vapour from the air at 2 m through
  the aerodynamic resistance 208/u2 and its leaves' boundary layers in series. Dew forms where the surface falls below
  the air's dew point.
- **Dew on an exposed leaf.** A single leaf held in the open air at the top of a plant, its upper face to the sky and
  its lower face to surroundings at the air's temperature, with no soil heat, exchanging through its own boundary layer
  alone. Its thin boundary layer holds it close to the air's temperature, so it gathers dew only in nearly saturated
  air. Its conductances are a flat plate's laminar ones (Nusselt and Sherwood numbers 0.664 Re^1/2 Pr^1/3 and Sc^1/3)
  times 1.4 for outdoor air, which at sea level on Earth is Campbell and Norman's 0.189 (u/d)^1/2 mol m-2 s-1 for
  heat, with the Moon's air at its own pressure, and free convection at lunar gravity where the wind fails.
- **Dew on bare soil.** The ground at the GCM's own temperature takes vapour through the bare ground's aerodynamic
  resistance, ln(2/0.001) ln(2/0.0001) / (k^2 u2): an upper bound, since the latent heat released would warm it.
- **Fog.** The air at the ground saturates where its humidity reaches saturation at the ground's temperature; fog
  water reaches the plants at its liquid water content times a deposition velocity.

Saturation is Bolton's over water, as climate_targets takes it.
"""
from __future__ import annotations

import numpy as np

from research.studies.climate_targets import extract as CT
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY, STEFAN_BOLTZMANN

LV = 2.45e6                 # J/kg, FAO-56's latent heat near 20 C
CP = 1013.0                 # J/kg/K, FAO-56 (equation 8's 1.013e-3 MJ/kg/K)
RD = CT.RD                  # the GCM's dry-air gas constant for the design air
EPSILON = 0.622
KARMAN = 0.41
REFERENCE_ALBEDO = 0.23
RA_REFERENCE = 208.0        # s/m times u2: FAO-56's short reference surface
U2_FLOOR = 0.5              # m/s, FAO-56's lower limit in the equation
RS_DAY, RS_NIGHT = 50.0, 200.0
G_DAY, G_NIGHT = 0.1, 0.5   # the soil's share of the net radiation by day and of its loss by night
SOIL_Z0M, SOIL_Z0H = 0.001, 0.0001
LEAF_EMISSIVITY = 0.97
OUTDOOR_FACTOR = 1.4        # Campbell and Norman: outdoor leaf conductances against a laminar plate's
LEAF_WIDTHS_M = (0.02, 0.05, 0.10)
CHARACTERISTIC = 0.72       # Campbell and Norman: a leaf's characteristic dimension over its width
CANOPY_LAI = 3.0


def saturation_pressure(t):
    return CT.es_water(t)


def saturation_humidity(t, p):
    return CT.qsat(t, p)


def slope(t):
    """d e_s / dT (Pa/K) for Bolton's form."""
    return saturation_pressure(t) * 17.67 * 243.5 / (t - 29.65) ** 2


def vapour_pressure(q, p):
    return q * p / (EPSILON + (1.0 - EPSILON) * q)


def dew_point(q, p):
    """The dew point (K) of air at humidity q and pressure p, inverting Bolton's form."""
    e = np.maximum(vapour_pressure(q, p), 1e-3)
    x = np.log(e / 611.2)
    return (243.5 * x / (17.67 - x)) + 273.15


def density(t, q, p):
    return p / (RD * t * (1.0 + 0.608 * q))


def air_properties(t, p):
    """Kinematic viscosity, thermal diffusivity and vapour diffusivity (m2/s) of the air at t (K) and p (Pa):
    Sutherland's viscosity, the US Standard Atmosphere's conductivity and a vapour diffusivity of 21.2 mm2/s at 0 C
    and 101.325 kPa going as T^1.75 / p."""
    mu = 1.458e-6 * t ** 1.5 / (t + 110.4)
    k = 2.64638e-3 * t ** 1.5 / (t + 245.4 * 10.0 ** (-12.0 / t))
    rho = p / (RD * t)
    return mu / rho, k / (rho * CP), 21.2e-6 * (t / 273.15) ** 1.75 * (101325.0 / p)


def leaf_conductances(t, p, wind, width, dt, gravity=MOON_SURFACE_GRAVITY):
    """Heat and vapour conductances (m/s) of one face of a leaf of the given width in wind (m/s), with free convection
    at the leaf-air difference dt (K) combined as the cube root of the sum of the cubes."""
    nu, kappa, dv = air_properties(t, p)
    d = CHARACTERISTIC * width
    re = np.maximum(wind, 1e-6) * d / nu
    forced_nu = OUTDOOR_FACTOR * 0.664 * np.sqrt(re) * (nu / kappa) ** (1.0 / 3.0)
    forced_sh = OUTDOOR_FACTOR * 0.664 * np.sqrt(re) * (nu / dv) ** (1.0 / 3.0)
    gr = gravity * np.abs(dt) / t * d ** 3 / nu ** 2
    free_nu = 0.54 * (gr * nu / kappa) ** 0.25
    free_sh = 0.54 * (gr * nu / dv) ** 0.25
    heat = np.cbrt(forced_nu ** 3 + free_nu ** 3) * kappa / d
    vapour = np.cbrt(forced_sh ** 3 + free_sh ** 3) * dv / d
    return heat, vapour


def wind_2m(u):
    return np.maximum(u, U2_FLOOR)


def demand(net_sun, sky, t, q, p, u2, day, resistance=(RS_DAY, RS_NIGHT), soil=(G_DAY, G_NIGHT)):
    """Penman-Monteith demand (kg/m2/s, mm/s) of the short reference surface: net_sun the sunlight it absorbs at its
    own albedo (W/m2), sky the downward longwave, day a mask of daylit periods; the surface resistance (s/m) and the
    soil's share of the net radiation by day and by night."""
    q = np.minimum(q, saturation_humidity(t, p))
    rn = net_sun + sky - STEFAN_BOLTZMANN * t ** 4
    g = np.where(day, soil[0] * rn, soil[1] * rn)
    u = wind_2m(u2)
    ra = RA_REFERENCE / u
    rs = np.where(day, resistance[0], resistance[1])
    gamma = CP * p / (EPSILON * LV)
    vpd = saturation_pressure(t) - vapour_pressure(q, p)
    rho = density(t, q, p)
    flux = (slope(t) * (rn - g) + rho * CP * vpd / ra) / (slope(t) + gamma * (1.0 + rs / ra))
    return np.maximum(flux, 0.0) / LV


def _bisect(balance, lo, hi, steps=48):
    """The root of a balance that is positive below it and negative above, by bisection between lo and hi (K)."""
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        f = balance(mid)
        lo, hi = np.where(f > 0, mid, lo), np.where(f > 0, hi, mid)
    return 0.5 * (lo + hi)


def canopy_dew(t, q, p, sky, u2, width=LEAF_WIDTHS_M[1], soil_share=G_NIGHT, lai=CANOPY_LAI,
               gravity=MOON_SURFACE_GRAVITY):
    """Dew on a grass-like canopy at night (kg/m2/s of ground) and its surface temperature (K)."""
    q = np.minimum(q, saturation_humidity(t, p))
    rho = density(t, q, p)
    ra = RA_REFERENCE / wind_2m(u2)

    def fluxes(tc):
        heat, vapour = leaf_conductances(t, p, u2, width, t - tc, gravity)
        rh, rv = ra + 1.0 / (2.0 * lai * heat), ra + 1.0 / (2.0 * lai * vapour)
        dew = np.maximum(rho * (q - saturation_humidity(tc, p)), 0.0) / rv
        return rh, dew

    def balance(tc):
        rh, dew = fluxes(tc)
        loss = LEAF_EMISSIVITY * (sky - STEFAN_BOLTZMANN * tc ** 4)
        return (1.0 - soil_share) * loss + rho * CP * (t - tc) / rh + LV * dew

    tc = _bisect(balance, t - 25.0, t + 5.0)
    return fluxes(tc)[1], tc


def leaf_dew(t, q, p, sky, wind, width=LEAF_WIDTHS_M[1], gravity=MOON_SURFACE_GRAVITY):
    """Dew on an exposed leaf (kg/m2/s of one face's area, gathered on both faces) and its temperature (K)."""
    q = np.minimum(q, saturation_humidity(t, p))
    rho = density(t, q, p)

    def fluxes(tl):
        heat, vapour = leaf_conductances(t, p, wind, width, t - tl, gravity)
        dew = 2.0 * vapour * np.maximum(rho * (q - saturation_humidity(tl, p)), 0.0)
        return heat, dew

    def balance(tl):
        heat, dew = fluxes(tl)
        upper = LEAF_EMISSIVITY * (sky - STEFAN_BOLTZMANN * tl ** 4)
        lower = LEAF_EMISSIVITY * STEFAN_BOLTZMANN * (t ** 4 - tl ** 4)
        return upper + lower + 2.0 * heat * rho * CP * (t - tl) + LV * dew

    tl = _bisect(balance, t - 25.0, t + 5.0)
    return fluxes(tl)[1], tl


def soil_dew(t_ground, q, p, u2, t_air):
    """Dew on bare ground at the GCM's own ground temperature (kg/m2/s), with no warming by the latent heat it
    releases: an upper bound."""
    q = np.minimum(q, saturation_humidity(t_air, p))
    ra = np.log(2.0 / SOIL_Z0M) * np.log(2.0 / SOIL_Z0H) / (KARMAN ** 2 * wind_2m(u2))
    return np.maximum(density(t_air, q, p) * (q - saturation_humidity(t_ground, p)), 0.0) / ra


def saturated_at_ground(t_ground, q, p):
    """Whether the air in contact with the ground would saturate at the ground's temperature."""
    return q >= saturation_humidity(t_ground, p)


def earth_night(rh, deficit_w_m2, hours, t=288.15, u2=1.0, p=101325.0):
    """Dew (mm) on the canopy and on the exposed leaf through an Earth night at relative humidity rh, the sky
    deficit (W/m2) below the air's own emission, Earth's gravity: the check against measured dew."""
    t, rh = np.asarray(t, float), np.asarray(rh, float)
    q = EPSILON * rh * saturation_pressure(t) / (p - (1.0 - EPSILON) * rh * saturation_pressure(t))
    sky = STEFAN_BOLTZMANN * t ** 4 - deficit_w_m2
    c, _ = canopy_dew(t, q, p, sky, u2, gravity=STANDARD_GRAVITY)
    leaf, _ = leaf_dew(t, q, p, sky, u2, gravity=STANDARD_GRAVITY)
    return c * hours * 3600.0, leaf * hours * 3600.0
