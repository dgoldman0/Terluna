#!/usr/bin/env python3
"""The high cloud's ice at a world's gravity: a slab of cold air rising through a layer of the upper air, freezing
its haze droplets into crystals that grow and fall out of it, for the GCM's factors on its cold cloud.

A slab of pressure thickness dp_slab rises at a steady pressure velocity omega (Pa/s) from ice saturation at the
bottom of a layer to its top, cooling along its dry adiabat and warmed by the heat its ice releases. Its sulphate
solution droplets freeze homogeneously at Koop and co-authors' (2000) rate J(da_w), da_w = (e - e_si) / e_sl, the
droplets' water activity less ice's at the same temperature, each droplet freezing at J V per second. The crystals
grow from the vapour by diffusion to their capacitance, with the vapour's gas-kinetic limit at a deposition
coefficient alpha, D* = D / (1 + 4 D / (alpha v r)), and the conduction of their latent heat (Pruppacher and Klett
2010, chapter 13), against Murphy and Koop's (2005) saturation over ice and supercooled water. They fall at Heymsfield
and Westbrook's (2010) speed and leave the slab through its base, a fraction rho g u dt / dp_slab each step. A
crystal is a solid ice sphere up to 97 µm, where Brown and Francis's (1995) mass law for larger ice meets the sphere's
and takes over; its projected area ratio is 1 up to that size and falls as a power of the size to the aggregates'
0.35 at 1 mm (flakes.AREA_RATIO), its capacitance from 0.5 D to the aggregates' 0.25 D alike. The air's viscosity,
vapour diffusivity and conductivity are the raindrop solver's (drops.transport).

Gravity reaches the crystals through their fall speed and the slab's mass: equal pressure velocities and slabs on two
worlds cool their air alike and condense the same vapour per kilogram, while on the Moon a crystal falls more slowly
and a pascal of slab holds six times the air. Each run returns the slab's ice over the cloudy part of its rise: the
mean ice mixing ratio, the effective radius 3 q_i / (4 rho_ice A), A the crystals' projected area per kilogram of air,
the number of crystals and the share of the condensed ice that fell out. Left out: aggregation of the crystals,
heterogeneous ice nuclei, turbulence and waves in the cooling, and ice falling in from above.
"""
from __future__ import annotations

import math

import numpy as np

from shared.constants import GAS_CONSTANT
from atmosphere.radiative_convective.thermodynamics import MOLAR_MASS
from atmosphere.rain import drops
from atmosphere.snow import flakes as sf

RHO_ICE = 917.0                              # kg/m3
R_VAPOUR = sf.R_VAPOUR
LATENT_SUBLIMATION = drops.latent_heat(sf.T0) + sf.LATENT_FUSION   # J/kg, held at its 0 °C value as flakes holds it
SPHERE_LIMIT_M = (6.0 * sf.MASS_COEFFICIENT / (math.pi * RHO_ICE)) ** (1.0 / (3.0 - sf.MASS_EXPONENT))   # 97 µm
AGGREGATE_M = 1.0e-3                         # where the area ratio reaches the aggregates' 0.35
AREA_EXPONENT = math.log(sf.AREA_RATIO) / math.log(AGGREGATE_M / SPHERE_LIMIT_M)
KOOP = (-906.7, 8502.0, -26924.0, 29180.0)   # log10 J (cm-3 s-1) as a cubic in da_w, Koop and co-authors (2000)
KOOP_RANGE = (0.26, 0.34)                    # the fit's range of da_w
HAZE_RADIUS_M = 0.15e-6                      # the solution droplets' radius as they freeze
HAZE_PER_KG = 5.0e8                          # sulphate droplets per kilogram of air, about 150 a cubic centimetre aloft


def koop_rate(da_w):
    """Homogeneous freezing rate (per m3 of solution per second) at a water-activity difference da_w."""
    x = np.minimum(np.asarray(da_w, float), KOOP_RANGE[1])
    log_j = KOOP[0] + x * (KOOP[1] + x * (KOOP[2] + x * KOOP[3]))
    return np.where(x < KOOP_RANGE[0], 0.0, 1e6 * 10.0 ** log_j)


def crystal_mass(d):
    d = np.asarray(d, float)
    return np.where(d <= SPHERE_LIMIT_M, RHO_ICE * math.pi / 6.0 * d ** 3, sf.MASS_COEFFICIENT * d ** sf.MASS_EXPONENT)


def crystal_size(m):
    m = np.asarray(m, float)
    sphere = (6.0 * m / (math.pi * RHO_ICE)) ** (1.0 / 3.0)
    return np.where(sphere <= SPHERE_LIMIT_M, sphere, (m / sf.MASS_COEFFICIENT) ** (1.0 / sf.MASS_EXPONENT))


def area_ratio(d):
    d = np.asarray(d, float)
    return np.clip((d / SPHERE_LIMIT_M) ** AREA_EXPONENT, sf.AREA_RATIO, 1.0)


def capacitance(d):
    """Capacitance (m): 0.5 D for the spheres, falling with the area ratio to the aggregates' 0.25 D."""
    d = np.asarray(d, float)
    share = np.log(area_ratio(d)) / math.log(sf.AREA_RATIO)            # 0 for spheres, 1 for aggregates
    return d * (0.5 - 0.25 * share)


def fall_speed(d, m, gravity, s):
    """Terminal speed (m/s) and Reynolds number by Heymsfield and Westbrook (2010), with each crystal's area ratio."""
    x = 8.0 * m * gravity * s['rho'] / (math.pi * np.sqrt(area_ratio(d)) * s['eta'] ** 2)
    re = sf.DELTA0 ** 2 / 4.0 * (np.sqrt(1.0 + 4.0 * np.sqrt(x) / (sf.DELTA0 ** 2 * math.sqrt(sf.C0))) - 1.0) ** 2
    return s['eta'] * re / (s['rho'] * d), re


def growth_coefficient(d, re, s, alpha):
    """kg/s a crystal gains per unit of ice supersaturation (dm/dt = k (S_i - 1)), gas-kinetic limit included."""
    t = s['t']
    c = capacitance(d)
    ls = LATENT_SUBLIMATION
    esi = sf.ice_saturation_pressure(t)
    speed_v = math.sqrt(8.0 * R_VAPOUR * t / math.pi)
    speed_a = math.sqrt(8.0 * s['r_air'] * t / math.pi)
    dv = s['dv'] / (1.0 + 4.0 * s['dv'] / (alpha * speed_v * c))
    ka = s['ka'] / (1.0 + 4.0 * s['ka'] / (s['rho'] * s['cp'] * speed_a * c))   # thermal accommodation 1
    x = (s['eta'] / (s['rho'] * s['dv'])) ** (1.0 / 3.0) * np.sqrt(re)
    vent = np.where(x < 1.0, 1.0 + 0.14 * x * x, 0.86 + 0.28 * x)
    resistance = (ls / (R_VAPOUR * t) - 1.0) * ls / (ka * t) + R_VAPOUR * t / (dv * esi)
    return 4.0 * math.pi * c * vent / resistance


def _lag(z):
    """1 - (1 - exp(-z)) / z, accurate for small z: the share of a step's supply deposited within the step."""
    return z / 2.0 - z * z / 6.0 + z ** 3 / 24.0 if z < 1e-3 else (z + math.expm1(-z)) / z


class Air:
    """A world's dry air and gravity, giving the slab's state at a temperature, pressure and vapour pressure."""

    def __init__(self, gravity, air):
        self.gravity, self.air = gravity, air
        self.r_air = GAS_CONSTANT / air.molar_mass

    def state(self, t, p, e):
        eta, dv, ka, cp = drops.transport(t, p, e / p, self.air.fractions)
        rho = (p - e) / (self.r_air * t) + e / (R_VAPOUR * t)
        return dict(t=t, p=p, e=e, rho=rho, eta=eta, dv=dv, ka=ka, cp=cp, r_air=self.r_air)


def slab(world, t_bottom, p_bottom, p_top, omega, dp_slab, alpha=0.5, haze_per_kg=HAZE_PER_KG,
         haze_radius_m=HAZE_RADIUS_M, max_cohorts=600):
    """A slab of air rising at omega (Pa/s, the magnitude) from ice saturation at (t_bottom, p_bottom) to p_top:
    its ice over the cloudy part of the rise, where the ice holds more than a thousandth of the slab's first vapour,
    and what it did along the way."""
    eps = MOLAR_MASS['H2O'] / world.air.molar_mass
    t, p = t_bottom, p_bottom
    q = eps * sf.ice_saturation_pressure(t) / p                        # vapour, kg/kg (q << 1 aloft)
    q_start = q
    haze = haze_per_kg
    v_haze = 4.0 / 3.0 * math.pi * haze_radius_m ** 3
    m0 = RHO_ICE * v_haze                                              # a frozen droplet's mass
    n = np.empty(0)                                                    # crystals per kg of air in each cohort
    m = np.empty(0)                                                    # each cohort's crystal mass
    born = np.empty(0)
    fallen, frozen_total, elapsed, nucleating = 0.0, 0.0, 0.0, 0.0
    sums = dict(time=0.0, q=0.0, area=0.0, number=0.0)
    s_max = 1.0
    while p > p_top * (1.0 + 1e-12):
        state = world.state(t, p, q * p / eps)
        esi = sf.ice_saturation_pressure(t)
        qsi = eps * esi / p
        s_max = max(s_max, q / qsi)
        da_w = (q * p / eps - esi) / sf.liquid_saturation_pressure(t)
        rate = float(koop_rate(da_w)) * v_haze                         # each droplet's freezing rate, 1/s
        cp, ls = state['cp'], LATENT_SUBLIMATION
        kappa = world.r_air / cp
        supply = qsi * omega / p * (ls * kappa / (R_VAPOUR * t) - 1.0)   # excess vapour the cooling makes, kg/kg/s
        dt = 2.0 if da_w > 0.25 else 60.0                              # short while the haze can freeze
        if n.size:
            d = crystal_size(m)
            k = growth_coefficient(d, fall_speed(d, m, world.gravity, state)[1], state, alpha)
            relative = k * abs(q / qsi - 1.0) / m                      # each cohort's growth per second
            dt = min(dt, max(0.5, 0.3 / float(relative.max()))) if relative.max() > 0 else dt
        dt = min(dt, (p - p_top) / omega)
        if rate > 0.0 and haze > 0.0:
            new = haze * -math.expm1(-rate * dt)
            haze -= new
            frozen_total += new
            nucleating += dt
            if n.size and elapsed - born[-1] < 120.0 and m[-1] < 10.0 * m0:
                m[-1] = (m[-1] * n[-1] + m0 * new) / (n[-1] + new)
                n[-1] += new
            else:
                n, m, born = np.append(n, new), np.append(m, m0), np.append(born, elapsed)
        deposited = 0.0
        if n.size:
            d = crystal_size(m)
            k = growth_coefficient(d, fall_speed(d, m, world.gravity, state)[1], state, alpha)
            uptake = float(np.sum(n * k)) / qsi                        # 1/s: the excess's relaxation rate
            z = uptake * dt
            # the excess x relaxes toward supply / uptake: what deposits over the step, exact as uptake -> 0
            deposited = max((q - qsi) * -math.expm1(-z) + supply * dt * _lag(z), -float(np.sum(n * m)))
            m = np.maximum(m + deposited * k / (uptake * qsi), 0.0)    # shared by each crystal's demand
        q = max(q - deposited, 0.0)
        p_new = p - omega * dt
        t = t * (p_new / p) ** kappa + ls * deposited / cp           # its dry adiabat and the heat of its ice
        p = p_new
        if n.size:                                                     # fall-out through the slab's base
            d = crystal_size(np.maximum(m, m0))
            u, _ = fall_speed(d, np.maximum(m, m0), world.gravity, state)
            keep = np.exp(-state['rho'] * world.gravity * u * dt / dp_slab)
            fallen += float(np.sum(n * (1.0 - keep) * m))
            n = n * keep
            alive = (m > 0.5 * m0) & (n > 1e-9 * float(n.sum()))
            q += float(np.sum(n[~alive] * m[~alive]))                  # what is left of a vanishing cohort
            n, m, born = n[alive], m[alive], born[alive]
            if n.size > max_cohorts:                                   # merge neighbours in mass
                order = np.argsort(m)
                n, m, born = n[order], m[order], born[order]
                pairs = n.size // 2
                nn = n[:2 * pairs].reshape(-1, 2).sum(1)
                mm = (n[:2 * pairs] * m[:2 * pairs]).reshape(-1, 2).sum(1) / nn
                bb = born[:2 * pairs].reshape(-1, 2).max(1)
                n, m, born = (np.concatenate([nn, n[2 * pairs:]]), np.concatenate([mm, m[2 * pairs:]]),
                              np.concatenate([bb, born[2 * pairs:]]))
        elapsed += dt
        qi = float(np.sum(n * m))
        if qi > 1e-3 * q_start:
            d = crystal_size(m)
            sums['time'] += dt
            sums['q'] += qi * dt
            sums['area'] += float(np.sum(n * area_ratio(d) * math.pi / 4.0 * d * d)) * dt
            sums['number'] += float(n.sum()) * dt
    tc = sums['time']
    qi_mean = sums['q'] / tc if tc else 0.0
    area_mean = sums['area'] / tc if tc else 0.0
    ice_end = float(np.sum(n * m))
    return dict(gravity_m_s2=world.gravity, t_bottom_k=t_bottom, t_top_k=t, p_bottom_pa=p_bottom, p_top_pa=p_top,
                omega_pa_s=omega, slab_pa=dp_slab, alpha=alpha, hours=elapsed / 3600.0,
                cloudy_share=tc / elapsed if elapsed else 0.0, ice_kg_kg=qi_mean,
                effective_radius_um=1e6 * 3.0 * qi_mean / (4.0 * RHO_ICE * area_mean) if area_mean else 0.0,
                crystals_per_kg=sums['number'] / tc if tc else 0.0, frozen_per_kg=frozen_total,
                fallen_share=fallen / max(fallen + ice_end, 1e-30), vapour_start_kg_kg=q_start,
                vapour_end_kg_kg=q, ice_end_kg_kg=ice_end, fallen_kg_kg=fallen,
                largest_supersaturation=s_max - 1.0, nucleating_hours=nucleating / 3600.0)
