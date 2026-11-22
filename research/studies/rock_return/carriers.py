"""The living return of rock, part by part: what a carrier spends, what it can reach, whether the loop holds itself,
and what keeps the returned rock in place. Pure functions; run.py puts the design surface's figures through them.

A carrier here is an animal that takes settled sediment from where water left it and carries it back to the ground it
came from, as part of how it lives: a flyer that feeds or gathers nest earth at a delta or a hollow and nests on high
ground, or a diver that brings a floor's sediment up to the shallows. Flight follows the sky ecology study's rule,
weight over glide ratio over muscle efficiency for each metre flown, with the climb paid in full.
"""
from __future__ import annotations

import math

import numpy as np

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

CALCITE_G_PER_MOL = 100.087
NITROGEN_G_PER_MOL = 14.007
GRAIN_KG_M3 = 2600.0
WATER_KG_M3 = 1000.0


def payload_cost_j_per_kg(haul_m, lift_m, glide, efficiency, gravity=MOON_SURFACE_GRAVITY):
    """Metabolic energy to fly one kilogram a distance over the ground and up a height."""
    if glide <= 0 or not 0 < efficiency <= 1:
        raise ValueError('a positive glide ratio and an efficiency in (0, 1] are required')
    return gravity * haul_m / (glide * efficiency) + gravity * lift_m / efficiency


def delivered_cost_j_per_kg(haul_m, lift_m, glide, efficiency, payload_over_body, glides_home=1.0,
                            gravity=MOON_SURFACE_GRAVITY):
    """Metabolic energy for each kilogram of rock delivered: the load and the carrier's own body flown out and up,
    and the body flown back level over the share of the way home that a glide from the height gained does not cover."""
    if payload_over_body <= 0 or not 0 <= glides_home <= 1:
        raise ValueError('a positive load and a share in [0, 1] are required')
    out = payload_cost_j_per_kg(haul_m, lift_m, glide, efficiency, gravity)
    back = (1.0 - glides_home) * gravity * haul_m / (glide * efficiency)
    return out + (out + back) / payload_over_body


def lift_in_water_j_per_kg(depth_m, efficiency, gravity=MOON_SURFACE_GRAVITY, grain=GRAIN_KG_M3, water=WATER_KG_M3):
    """Metabolic energy to raise a kilogram of mineral grains through water: its weight less its buoyancy."""
    return gravity * (1.0 - water / grain) * depth_m / efficiency


def pressure_atm(depth_m, gravity=MOON_SURFACE_GRAVITY, surface_atm=1.0, water=WATER_KG_M3):
    return surface_atm + water * gravity * depth_m / 101325.0


def depth_at_earths_pressure(earth_depth_m, surface_atm=1.0):
    """The depth in the Moon's water whose pressure is that of a depth in Earth's: what a diver's body bears."""
    return earth_depth_m * STANDARD_GRAVITY / MOON_SURFACE_GRAVITY + (1.0 - surface_atm) * 101325.0 / (WATER_KG_M3 * MOON_SURFACE_GRAVITY)


def reach(on_dry_ground, in_water, plain, plain_and_front, floor_within_a_dive):
    """The shares of the rock that each kind of carrier can reach, given what settles on dry ground and in standing
    water, the share of the water's part that the delta plain keeps, the share the plain and the delta front keep
    together, and the share of the floors beyond that lies within a dive."""
    if not 0 <= plain <= plain_and_front <= 1 or not 0 <= floor_within_a_dive <= 1:
        raise ValueError('shares out of order')
    flyers = on_dry_ground + in_water * plain
    fronts = in_water * (plain_and_front - plain)
    floors = in_water * (1.0 - plain_and_front)
    return dict(flyers_and_waders=flyers, with_divers_at_the_delta_fronts=flyers + fronts,
                with_divers_on_the_floors_within_a_dive=flyers + fronts + floors * floor_within_a_dive,
                left_on_the_deeper_floors=floors * (1.0 - floor_within_a_dive))


def carriers_for(rock_kg_yr, body_kg, payload_over_body, trips_a_year):
    """How many carriers of one body deliver a yearly mass, each carrying its load on so many trips."""
    return rock_kg_yr / (body_kg * payload_over_body * trips_a_year)


def hours_a_year(trips, haul_m, speed_m_s, handling_s):
    return trips * (2.0 * haul_m / speed_m_s + handling_s) / 3600.0


# ---- staying put -----------------------------------------------------------------------------------------------
def crest_slope(lowering_m_yr, reach_m, kappa_m2_yr):
    """The slope soil creep needs just below a line where carriers set the rock down, if that line feeds the ground for
    reach_m on each side and the ground is losing lowering_m_yr: the flux to pass is the loss times the reach, and
    creep passes kappa times the slope."""
    return lowering_m_yr * reach_m / kappa_m2_yr


def reach_for_slope(slope, lowering_m_yr, kappa_m2_yr):
    """How far to each side a line of returned rock can feed the ground without creep needing more than a slope."""
    return slope * kappa_m2_yr / lowering_m_yr


def calcite_by_urea_kg(nitrogen_kg):
    """Calcite that bacteria splitting urea can lay down: one carbonate for each urea, two nitrogen atoms."""
    return nitrogen_kg * CALCITE_G_PER_MOL / (2.0 * NITROGEN_G_PER_MOL)


# ---- the loop --------------------------------------------------------------------------------------------------
def loop_equilibrium(supply, duty, growth, death, half):
    """The standing deposit and the carriers at which the loop rests. Sediment arrives at supply; carriers take
    duty * B * S / (S + half) of it a year; they multiply at growth * S / (S + half) and die at death."""
    if not growth > death > 0:
        raise ValueError('the carriers must be able to outgrow their deaths')
    deposit = half * death / (growth - death)
    carriers = supply * growth / (duty * death)
    return deposit, carriers


def loop_eigenvalues(supply, duty, growth, death, half):
    """The two rates at which a small disturbance of the resting loop changes (1/yr); both real parts are negative
    for every positive setting, so the loop returns to rest."""
    s, b = loop_equilibrium(supply, duty, growth, death, half)
    d = half / (s + half) ** 2
    jac = np.array([[-duty * b * d, -duty * s / (s + half)], [b * growth * d, growth * s / (s + half) - death]])
    return np.linalg.eigvals(jac)


def loop_run(supply, duty, growth, death, half, deposit, carriers, years, dt=0.05):
    """The deposit and the carriers through time from a starting state, with supply a number or a function of the
    year. Returns the years, the deposit, the carriers and what the carriers return each year."""
    n = int(round(years / dt))
    t = np.arange(n + 1) * dt
    s_out, b_out, r_out = np.empty(n + 1), np.empty(n + 1), np.empty(n + 1)
    s, b = float(deposit), float(carriers)
    feed = supply if callable(supply) else (lambda year: supply)

    def rates(s, b, year):
        take = s / (s + half) if s > 0 else 0.0
        return feed(year) - duty * b * take, b * (growth * take - death)
    for i in range(n + 1):
        s_out[i], b_out[i] = s, b
        r_out[i] = duty * b * (s / (s + half) if s > 0 else 0.0)
        if i == n:
            break
        k1 = rates(s, b, t[i])
        k2 = rates(s + 0.5 * dt * k1[0], b + 0.5 * dt * k1[1], t[i] + 0.5 * dt)
        k3 = rates(s + 0.5 * dt * k2[0], b + 0.5 * dt * k2[1], t[i] + 0.5 * dt)
        k4 = rates(s + dt * k3[0], b + dt * k3[1], t[i] + dt)
        s = max(s + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0, 0.0)
        b = max(b + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0, 0.0)
    return t, s_out, b_out, r_out


def years_until_within(t, series, target, tolerance=0.05):
    """The first year after which a series stays within a tolerance of its target."""
    off = np.abs(series - target) > tolerance * abs(target)
    if not off.any():
        return 0.0
    last = int(np.flatnonzero(off)[-1])
    return float(t[last + 1]) if last + 1 < t.size else math.inf


def open_after(span_yr, fill_years_by_share, returned_share):
    """The share of the lakes' area still open after a span, when carriers return a share of what the lakes keep, from
    the years by which each share of the area fills with nothing returned (read on a logarithmic scale of years)."""
    if not 0 <= returned_share < 1:
        return 1.0
    pairs = sorted((float(q), y) for q, y in fill_years_by_share.items() if y)
    shares = np.array([q for q, _ in pairs])
    years = np.log10([y for _, y in pairs])
    t = math.log10(span_yr * (1.0 - returned_share))
    return 1.0 - float(np.interp(t, years, shares, left=0.0, right=shares[-1]))
