"""The river deltas screen's physics, with gravity, rotation and the waters written in.

One set of functions serves the Moon and the Earth control. Each relation carried from Earth is used in its
dimensionless form with g explicit; where an Earth fit hides g inside a coefficient (Chezy's C, Manning's n, the
CERC and Ashton-Murray constants), it is rewritten here so that g shows. Sources and how far each was read are in
sources.json; design guesses are marked where they are set, in run.py.
"""
from __future__ import annotations

import math

import numpy as np

KAPPA = 0.41                     # von Karman's constant (the standard 0.40-0.41)
RHO_W0 = 999.974950              # Tanaka et al. (2001): pure water's density at its maximum, kg/m3
CERC_K = 0.39                    # CERC's coefficient with the significant height (Shore Protection Manual 1984)
BREAKER_INDEX = 0.78             # breaking height over depth (McCowan's solitary-wave limit, as CERC is used)
# Ferguson & Church (2004): C1 = 18 (Stokes) and C2 = 1.0 for natural sieved grains, 0.4 for smooth spheres and
# 1.2 for very angular grains.
FC_C1, FC_C2 = 18.0, 1.0
# Bankfull Shields number of single-thread sand-bed rivers, 17 Rep^-1/2 (Trampush et al. 2014's compilation, after
# Wilkerson & Parker 2011's finding that it falls as Rep^-1/2), as Hayden et al. (2019) carry it to Mars.
SAND_SHIELDS_A, SAND_SHIELDS_N = 17.0, -0.5
# Parker et al. (2007), gravel beds: H = 0.382 Q^0.4 / g^0.2, B = 4.63 Q^0.4 / g^0.2 Qhat^0.0667,
# S = 0.101 Qhat^-0.344, with Qhat = Q / (sqrt(g D) D^2).
PARKER_GRAVEL = dict(depth=0.382, width=4.63, width_n=0.0667, slope=0.101, slope_n=-0.344)
THRESHOLD_RATIO = 1.2            # bankfull over critical stress of the channel's most resistant material (Dunne &
                                 # Jerolmack 2020: their 18 sections near 1.2, none above 1.5)
PLUNGE_FROUDE = 0.5              # densimetric Froude number at the plunge line (Akiyama & Stefan 1984; Parker &
                                 # Toniolo 2007; 0.4-0.8 reported)
CM_COEFF = 0.17                  # Crosato & Mosselman (2009)


def water_density(t_c):
    """Pure water's density (kg/m3) at t_c degrees C (Tanaka et al. 2001)."""
    t = np.asarray(t_c, float)
    return RHO_W0 * (1.0 - (t - 3.983035) ** 2 * (t + 301.797) / (522528.9 * (t + 69.34881)))


def kinematic_viscosity(t_c):
    """Water's kinematic viscosity (m2/s): Vogel's fit mu = 2.414e-5 * 10**(247.8 / (T - 140)) Pa s over density."""
    t = np.asarray(t_c, float)
    return 2.414e-5 * 10.0 ** (247.8 / (t + 273.15 - 140.0)) / water_density(t)


def density(t_c, salt_g_kg=0.0, sediment_kg_m3=0.0, haline=0.8, rho_s=2650.0):
    """Water with dissolved salts (haline kg/m3 per g/kg) and suspended sediment of grain density rho_s."""
    rho = water_density(t_c)
    return rho + haline * np.asarray(salt_g_kg, float) + np.asarray(sediment_kg_m3, float) * (1.0 - rho / rho_s)


def settling_velocity(d, g, rho_s, nu, rho=1000.0, c1=FC_C1, c2=FC_C2):
    """Ferguson & Church (2004): w = R g D^2 / (C1 nu + sqrt(0.75 C2 R g D^3)), Stokes for fine grains."""
    r = rho_s / rho - 1.0
    d = np.asarray(d, float)
    return r * g * d * d / (c1 * nu + np.sqrt(0.75 * c2 * r * g * d ** 3))


def particle_reynolds(d, g, r, nu):
    return np.sqrt(r * g * d) * d / nu


def shields_critical(rep):
    """Brownlie's (1981) fit of the Shields curve, as written by Parker (2004): 0.5 [0.22 Rep^-0.6 +
    0.06 x 10^(-7.7 Rep^-0.6)]."""
    x = np.asarray(rep, float) ** -0.6
    return 0.5 * (0.22 * x + 0.06 * 10.0 ** (-7.7 * x))


def sand_bankfull_shields(rep):
    """Loose-bank scaling: the Earth compilation's bankfull Shields number at equal Rep, at least 1.2 critical."""
    rep = np.asarray(rep, float)
    return np.maximum(SAND_SHIELDS_A * rep ** SAND_SHIELDS_N, THRESHOLD_RATIO * shields_critical(rep))


def engelund_hansen(tau_star, cf):
    """Engelund & Hansen's (1967) total bed-material load, q* = (0.05 / Cf) tau*^(5/2)."""
    return 0.05 / cf * np.asarray(tau_star, float) ** 2.5


def unit_load(tau_star, cf, d, g, r):
    """Bed-material load per unit width, m3/s per m of grains: q* sqrt(R g D) D."""
    return engelund_hansen(tau_star, cf) * np.sqrt(r * g * d) * d


def parker_gravel(q, d, g):
    """Parker et al. (2007): bankfull depth, width and slope of a single-thread gravel river."""
    qhat = q / (np.sqrt(g * d) * d * d)
    base = np.asarray(q, float) ** 0.4 / g ** 0.2
    p = PARKER_GRAVEL
    return p['depth'] * base, p['width'] * base * qhat ** p['width_n'], p['slope'] * qhat ** p['slope_n']


def loose_channel(q, d, g, r, nu, cf, gravel_d=0.002):
    """Loose-bank (alluvial) scaling at bankfull. Width from Parker et al.'s (2007) dimensionless width,
    B = 4.63 Q^0.4 g^-0.2 Qhat^0.0667 (fitted on gravel; used here for sand beds as well); for sand beds the
    bankfull Shields number at equal Rep sets depth times slope and the friction coefficient the velocity, so that
    continuity gives the depth; gravel beds take Parker et al.'s depth and slope. Returns depth, width, velocity,
    slope and bed stress (Pa)."""
    q, d = np.asarray(q, float), np.asarray(d, float)
    hg, b, sg = parker_gravel(q, d, g)
    ts = sand_bankfull_shields(particle_reynolds(d, g, r, nu))
    hs = ts * r * d
    u = np.sqrt(g * hs / cf)
    h = q / (b * u)
    s = hs / h
    gravel = d >= gravel_d
    h, s = np.where(gravel, hg, h), np.where(gravel, sg, s)
    u = np.where(gravel, q / (b * hg), u)
    return dict(depth=h, width=b, velocity=u, slope=s, stress=1000.0 * g * h * s)


def threshold_channel(loose, tau_bed_c, tau_bank_c, bank_height_max):
    """Threshold channel (Dunne & Jerolmack 2020): bankfull stress 1.2 times the critical stress of the most
    resistant material lining it, bed or bank, in pascals that do not change with gravity. At the loose channel's
    discharge and slope, depth scales with the stress, velocity with its square root (the same friction
    coefficient) and width by continuity: H_t / H_l = tau_t / tau_l, B_t / B_l = (tau_l / tau_t)^(3/2). Channels
    deeper than a steep cohesive bank stands (bank_height_max) are flagged; Earth's big rivers keep taller banks
    by their slopes and layering, so the flag is not a cap."""
    tau = THRESHOLD_RATIO * np.maximum(tau_bed_c, tau_bank_c)
    k = tau / loose['stress']
    h = loose['depth'] * k
    return dict(depth=h, width=loose['width'] * k ** -1.5, velocity=loose['velocity'] * np.sqrt(k),
                slope=loose['slope'], stress=tau, bank_limited=h > bank_height_max)


def bank_height_limit(cohesion_pa, g, rho_bank=1800.0, stability_number=6.0):
    """Height a cohesive bank stands before it fails, N c / (rho g): Taylor's stability number of about 6 for a
    steep bank of weak friction (3.85 for a vertical cut). Lunar banks stand six times taller."""
    return stability_number * np.asarray(cohesion_pa, float) / (rho_bank * g)


def rouse(ws, u_star):
    return ws / (KAPPA * u_star)


def plunge_depth(q_unit, g_reduced, froude=PLUNGE_FROUDE):
    """Depth at which an inflow plunges: (q^2 / (g' Fp^2))^(1/3)."""
    return (np.asarray(q_unit, float) ** 2 / (np.maximum(g_reduced, 1e-30) * froude ** 2)) ** (1.0 / 3.0)


def rossby_radius(g_reduced, thickness, lat_deg, omega):
    """Internal Rossby radius of a plume, sqrt(g' h) / |f| with f = 2 Omega sin(latitude)."""
    f = 2.0 * omega * np.abs(np.sin(np.radians(lat_deg)))
    return np.sqrt(np.maximum(g_reduced, 0.0) * thickness) / np.maximum(f, 1e-12)


ANGLE_FACTOR = max(math.cos(x) ** 1.2 * math.sin(x) for x in np.linspace(0.0, math.pi / 2.0, 20001))


def longshore_max(h0, t, g, rho_s, rho=1000.0, k=CERC_K, gamma=BREAKER_INDEX):
    """Maximum potential alongshore transport (kg/s) of deep-water waves h0 (significant, m) of period t (s).

    CERC, I = K P_l, with the immersed weight transport I = (rho_s - rho) g Q, recast into deep-water terms with
    shoaling and refraction over shore-parallel contours (Ashton & Murray 2006): Q = K2 H0^(12/5) T^(1/5)
    cos^(6/5)(phi) sin(phi), where here K2 = K rho rho_s / (8 (rho_s - rho)) 2 pi (4 pi)^(-6/5) g^(3/5) gamma^(-2/5)
    holds gravity explicitly; the angle factor is its maximum, 0.4685 at 42.4 degrees."""
    k2 = k * rho * rho_s / (8.0 * (rho_s - rho)) * 2.0 * math.pi * (4.0 * math.pi) ** -1.2 * g ** 0.6 * gamma ** -0.4
    return k2 * np.asarray(h0, float) ** 2.4 * np.asarray(t, float) ** 0.2 * ANGLE_FACTOR


def wave_height_from_power(power_w_m, t, g, rho=1000.0):
    """Deep-water significant height carrying a mean wave power P (W/m): P = rho g^2 H^2 T / (64 pi)."""
    return np.sqrt(64.0 * math.pi * np.asarray(power_w_m, float) / (rho * g * g * np.asarray(t, float)))


def fetch_waves(u10, fetch_m, g):
    """JONSWAP fetch-limited growth with the fully developed limit (Shore Protection Manual 1984 forms):
    g H / U^2 = 0.0016 (g F / U^2)^(1/2) to 0.2433; g Tp / U = 0.2857 (g F / U^2)^(1/3) to 8.134."""
    x = g * np.asarray(fetch_m, float) / u10 ** 2
    h = np.minimum(1.6e-3 * np.sqrt(x), 0.2433) * u10 ** 2 / g
    tp = np.minimum(0.2857 * np.cbrt(x), 8.134) * u10 / g
    return h, tp


def tidal_discharge(amplitude, width, slope, period_s):
    """Tidal discharge amplitude at a mouth (m3/s): the channel's surface over the reach where its bed rises by the
    tidal amplitude, a / S, filled and drained sinusoidally, omega a^2 B / S (after Nienhuis et al. 2018, 2020)."""
    return 2.0 * math.pi / period_s * np.asarray(amplitude, float) ** 2 * width / np.maximum(slope, 1e-7)


def bar_mode(width, slope, q, d, r, cf, g, b_exp):
    """Crosato & Mosselman (2009): m^2 = 0.17 g (b - 3) B^3 S / (sqrt(Delta D50) C Q), Chezy's C = sqrt(g / Cf)."""
    chezy = np.sqrt(g / cf)
    m2 = CM_COEFF * g * (b_exp - 3.0) * width ** 3 * slope / (np.sqrt(r * d) * chezy * q)
    return np.sqrt(np.maximum(m2, 0.0))


def mouth_bar_distance(width, velocity, d, r, g, k=3.0, pi_ref=300.0):
    """Distance to the mouth bar: k mouth widths at Earth's typical jet, scaled by Edmonds & Slingerland's (2007)
    one-fifth power of the jet's momentum over the grains' weight, Pi = U^2 / (R g D)."""
    return k * width * (velocity ** 2 / (r * g * d) / pi_ref) ** 0.2


def avulsion_time(depth, belt_width, length, kept_kg_yr, dry_density, superelevation=1.0):
    """Years to build one superelevation over the channel belt along the avulsion length (Jerolmack & Mohrig 2007:
    T_A = h / v_A, with v_A the near-channel aggradation the kept load makes)."""
    return superelevation * depth * belt_width * length * dry_density / np.maximum(kept_kg_yr, 1e-30)


def front_failure_height(cohesion_pa, slope_deg, g, rho_sat=1600.0, rho=1000.0):
    """Height of an infinite submerged mud slope at which cohesion alone stops holding it: c / (gamma' sin b cos b)."""
    b = np.radians(slope_deg)
    return cohesion_pa / ((rho_sat - rho) * g * np.sin(b) * np.cos(b))


def lognormal_fractions(edges_m, median_m, sigma_ln):
    """Mass fractions of a lognormal grain-size distribution between the bin edges."""
    from math import erf
    z = (np.log(edges_m) - math.log(median_m)) / (sigma_ln * math.sqrt(2.0))
    cdf = 0.5 * (1.0 + np.array([erf(v) for v in z]))
    return np.diff(cdf)
