"""The Moon's settings for the ocean models.

Gravity, radius and rotation come from shared/constants.json; the seas'
water is the held fresh water of the sea-chemistry study (0.14 g/kg) at the
seas' median surface warmth at 301 K (29.6 C, the sea-oxygen study), whose
density the lake model's Kell (1975) formula gives. The tide runs use that
density only for the water's self-attraction and for stresses and energy;
the depth-averaged flow itself does not depend on it.

Bed friction is quadratic, tau = rho Cd |u| u, with Cd = 0.0025, the value
geography's linear tide check also assumes, and a bracket of 0.001 to 0.005
for beds smoother or rougher than the usual open-shelf value. Cd from a
logarithmic boundary layer holds no gravity, so it carries to the Moon as it
stands; Manning's n, which hides g, would not.
"""
from __future__ import annotations

import math

from biosphere.ecology import lake
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, SIDEREAL_MONTH_DAYS

GRAVITY = MOON_SURFACE_GRAVITY                       # m/s^2
RADIUS = MOON_RADIUS                                 # m
ROTATION_PERIOD_S = SIDEREAL_MONTH_DAYS * 86400.0    # one turn against the stars
OMEGA = 2.0 * math.pi / ROTATION_PERIOD_S            # rad/s

SEA_SALINITY_G_KG = 0.14     # the held fresh seas (research/studies/sea_chemistry)
SEA_SURFACE_C = 29.6         # the seas' median surface at 301 K (research/studies/sea_oxygen)
WATER_DENSITY = float(lake.density(SEA_SURFACE_C, SEA_SALINITY_G_KG))   # kg/m^3

DRAG = 0.0025                # quadratic bed drag coefficient, central
DRAG_BRACKET = (0.001, 0.005)
LINEAR_FRICTION_M_S = 2.5e-5  # geography/tide_dynamics.py: Cd 0.0025 times a 1 cm/s tidal flow

DAY_S = 86400.0


def wave_speed(depth_m):
    """Long-wave speed sqrt(g h), m/s."""
    return math.sqrt(GRAVITY * depth_m)


def coriolis(lat_deg):
    """f = 2 Omega sin(latitude), 1/s."""
    return 2.0 * OMEGA * math.sin(math.radians(lat_deg))
