/* Physical and astronomical constants shared by every domain (constants.json).
 * Python reads the same file through constants.py. Derived values are computed
 * here, never stored.
 */
import data from './constants.json' with { type: 'json' };

export const DATA = data;
export const BOLTZMANN = data.physics.boltzmann_J_K;
export const AVOGADRO = data.physics.avogadro_per_mol;
export const STEFAN_BOLTZMANN = data.physics.stefan_boltzmann_W_m2_K4;
export const PLANCK = data.physics.planck_J_s;
export const SPEED_OF_LIGHT = data.physics.speed_of_light_m_s;
export const STANDARD_GRAVITY = data.physics.standard_gravity_m_s2;

export const MOON_RADIUS = data.moon.radius_m;
export const MOON_GM = data.moon.gm_m3_s2;
export const MOON_SURFACE_GRAVITY = MOON_GM / MOON_RADIUS ** 2;
export const SYNODIC_MONTH_DAYS = data.moon.synodic_month_days;
export const SIDEREAL_MONTH_DAYS = data.moon.sidereal_month_days;

export const EARTH_RADIUS = data.earth.radius_m;
export const EARTH_GM = data.earth.gm_m3_s2;
export const EARTH_MOON_DISTANCE = data.earth.moon_mean_distance_m;

export const ANOMALISTIC_MONTH_DAYS = data.moon.anomalistic_month_days;
export const DRACONIC_MONTH_DAYS = data.moon.draconic_month_days;
export const LIBRATION_LONGITUDE_DEG = data.moon.optical_libration_longitude_deg;
export const LIBRATION_LATITUDE_DEG = data.moon.optical_libration_latitude_deg;
export const EARTH_SIDEREAL_DAY_S = data.earth.sidereal_day_s;
export const EARTH_OBLIQUITY_DEG = data.earth.obliquity_deg;
export const EARTH_GEOMETRIC_ALBEDO = data.earth.geometric_albedo_visible;

export const SOLAR_CONSTANT = data.sun.solar_constant_W_m2;
export const SUN_RADIUS = data.sun.radius_m;
export const AU = data.sun.au_m;

// Rounded lunar gravity still used by the sky solver, column model and immersion.
export const LEGACY_MOON_GRAVITY = data.legacy.moon_surface_gravity_rounded_m_s2.value;
