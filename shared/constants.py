"""Physical and astronomical constants shared by every domain.

Values come from constants.json, which constants.js also reads, so Python and the
browser use identical numbers. Derived values are computed here, never stored.
"""
from __future__ import annotations

import json
from pathlib import Path

DATA = json.loads(Path(__file__).with_name("constants.json").read_text())

BOLTZMANN = DATA["physics"]["boltzmann_J_K"]
AVOGADRO = DATA["physics"]["avogadro_per_mol"]
GAS_CONSTANT = BOLTZMANN * AVOGADRO
STEFAN_BOLTZMANN = DATA["physics"]["stefan_boltzmann_W_m2_K4"]
STANDARD_GRAVITY = DATA["physics"]["standard_gravity_m_s2"]

MOON_RADIUS = DATA["moon"]["radius_m"]
MOON_GM = DATA["moon"]["gm_m3_s2"]
MOON_SURFACE_GRAVITY = MOON_GM / MOON_RADIUS**2
SYNODIC_MONTH_DAYS = DATA["moon"]["synodic_month_days"]
SIDEREAL_MONTH_DAYS = DATA["moon"]["sidereal_month_days"]
MOON_EQUATOR_TO_ECLIPTIC_DEG = DATA["moon"]["equator_to_ecliptic_deg"]

EARTH_RADIUS = DATA["earth"]["radius_m"]
EARTH_GM = DATA["earth"]["gm_m3_s2"]
EARTH_MOON_DISTANCE = DATA["earth"]["moon_mean_distance_m"]

SOLAR_CONSTANT = DATA["sun"]["solar_constant_W_m2"]
SUN_RADIUS = DATA["sun"]["radius_m"]
AU = DATA["sun"]["au_m"]

# Rounded lunar gravity still used by the sky solver, column model and immersion.
LEGACY_MOON_GRAVITY = DATA["legacy"]["moon_surface_gravity_rounded_m_s2"]["value"]
