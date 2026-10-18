"""Dissolved gases in the lakes: solubility, the carbonate system, the hydrostatic bubble point, sediment ebullition
and a rising bubble, each with the gravity explicit.

- Solubility as Bunsen coefficients (mL of gas at STP per mL of water per atm of the gas): Weiss (1970) for N2, O2
  and Ar, Yamamoto et al. (1976) for CH4, both in fresh water; Plummer & Busenberg (1982) for CO2 (mol/kg/atm), the
  constants the carbon cycle study uses, with carbonic acid's K1 and K2 and Davies activities.
- The bubble point: water whose dissolved gases' partial pressures, with the vapour, sum past the local pressure (the
  surface pressure plus rho g z) can form bubbles; the gas saturation ratio is that sum over the local pressure, the
  criterion behind the hazard assessments of Lake Kivu (Schmid et al. 2005) and Lake Nyos.
- Sediment ebullition: the steady porewater balance of Bazhin (2003, 2010) as Langenegger et al. (2019) apply it, with
  production falling off exponentially in the sediment and N2 held at its value at the sediment surface (no
  stripping), so that bubbles form below the depth where methane reaches the critical concentration and every mole
  made below it leaves as gas: the ebullition share is exp(-u), with (1 + u) exp(-u) = 1 - F*/P.
- A rising bubble of CH4, N2 and O2 exchanging with the water around it: Tomiyama's (1998) drag for contaminated
  bubbles at the local gravity, Froessling's mass transfer, the ideal gas and the Laplace pressure.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from biosphere.ecology import lake as L

STP_MOLAR_VOLUME_L = 22.414           # L/mol, ideal gas at 0 C and 1 atm
ATM = 101325.0
GAS_R = 8.314462618
SURFACE_TENSION = 0.072               # N/m
BUNSEN = {                            # A1, A2, A3 in ln(beta) = A1 + A2 (100/T) + A3 ln(T/100), fresh water
    'N2': (-59.6274, 85.7661, 24.3696),       # Weiss (1970)
    'O2': (-58.3877, 85.8079, 23.8439),       # Weiss (1970)
    'Ar': (-55.6578, 82.0262, 22.5929),       # Weiss (1970)
    'CH4': (-67.1962, 99.1624, 27.9015),      # Yamamoto et al. (1976)
}
SCHMIDT_20C = {                       # Wanninkhof (1992), fresh water at 20 C
    'CH4': 1897.8 - 114.28 * 20 + 3.2902 * 400 - 0.039061 * 8000,
    'CO2': 1911.1 - 118.11 * 20 + 3.4527 * 400 - 0.041320 * 8000,
    'N2': 1970.7 - 131.45 * 20 + 4.1390 * 400 - 0.052106 * 8000,
    'O2': 1800.6 - 120.10 * 20 + 3.7818 * 400 - 0.047608 * 8000,
    'Ar': 1759.7 - 117.37 * 20 + 3.6959 * 400 - 0.046527 * 8000,
}
# Partial molar volumes behind the density effect of dissolved gas (Lake Kivu's coefficients, Schmid et al. 2004):
# CO2 adds 0.284e-3 and CH4 takes 1.25e-3 of the density per g/kg.
BETA_CO2_PER_G_KG, BETA_CH4_PER_G_KG = 0.284e-3, -1.25e-3


def solubility(gas, t_c):
    """mol per m3 of water per atm of the gas's partial pressure."""
    if gas == 'CO2':
        return co2_kh(t_c) * L.density(t_c)
    a1, a2, a3 = BUNSEN[gas]
    tk = np.asarray(t_c, dtype=float) + 273.15
    beta = np.exp(a1 + a2 * 100.0 / tk + a3 * np.log(tk / 100.0))
    return beta * 1000.0 / STP_MOLAR_VOLUME_L


def schmidt(gas, t_c):
    """Schmidt number in fresh water: Wanninkhof's value at 20 C carried by the Stokes-Einstein scaling, as the lake
    model carries O2's."""
    t = np.asarray(t_c, dtype=float)
    return (SCHMIDT_20C[gas] * (L.viscosity(t) / L.viscosity(20.0)) ** 2 * (L.density(20.0) / L.density(t))
            * (293.15 / (t + 273.15)))


def diffusivity(gas, t_c):
    t = np.asarray(t_c, dtype=float)
    return L.viscosity(t) / L.density(t) / schmidt(gas, t)


# --- the carbonate system ---------------------------------------------------------------------------------------

def co2_kh(t_c):
    """CO2's solubility, mol/kg/atm (Plummer & Busenberg 1982)."""
    tk = np.asarray(t_c, dtype=float) + 273.15
    return 10 ** (108.3865 + 0.01985076 * tk - 6919.53 / tk - 40.45154 * np.log10(tk) + 669365.0 / tk ** 2)


def carbonate_constants(t_c):
    tk = np.asarray(t_c, dtype=float) + 273.15
    lt = np.log10(tk)
    k1 = 10 ** (-356.3094 - 0.06091964 * tk + 21834.37 / tk + 126.8339 * lt - 1684915.0 / tk ** 2)
    k2 = 10 ** (-107.8871 - 0.03252849 * tk + 5151.79 / tk + 38.92561 * lt - 563713.9 / tk ** 2)
    kw = 10 ** -(4470.99 / tk - 6.0875 + 0.01706 * tk)
    davies = 0.4883 + 8.074e-4 * (tk - 273.15)
    return k1, k2, kw, davies


def gammas(t_c, ionic):
    _, _, _, a = carbonate_constants(t_c)
    r = math.sqrt(ionic)
    f = r / (1.0 + r) - 0.3 * ionic
    return 10 ** (-a * f), 10 ** (-4.0 * a * f)


def speciate(dic, ta, t_c, ionic=1e-3):
    """[H+] activity's pH and dissolved CO2 (mol/kg) for DIC and alkalinity (mol/kg), by bisection in pH; arrays
    broadcast."""
    dic, ta, t = np.broadcast_arrays(np.asarray(dic, float), np.asarray(ta, float), np.asarray(t_c, float))
    k1, k2, kw, _ = carbonate_constants(t)
    g1, g2 = gammas(t, ionic) if np.ndim(ionic) == 0 else (1.0, 1.0)
    k1, k2, kw = k1 / (g1 * g1), k2 / g2, kw / (g1 * g1)
    lo, hi = np.full(dic.shape, -14.0), np.full(dic.shape, -2.0)       # log10 of [H+]
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        h = 10.0 ** mid
        alk = dic * (k1 * h + 2.0 * k1 * k2) / (h * h + k1 * h + k1 * k2) + kw / h - h
        high = alk > ta                                  # too little H+: move up
        lo = np.where(high, mid, lo)
        hi = np.where(high, hi, mid)
    h = 10.0 ** (0.5 * (lo + hi))
    co2 = dic * h * h / (h * h + k1 * h + k1 * k2)
    return -np.log10(h * g1), co2


def alkalinity_at(ph, pco2_atm, t_c, ionic=1e-3):
    """Alkalinity (mol/kg) of water at a pH (activity) in equilibrium with pco2_atm."""
    k1, k2, kw, _ = carbonate_constants(t_c)
    g1, g2 = gammas(t_c, ionic)
    k1, k2, kw = k1 / (g1 * g1), k2 / g2, kw / (g1 * g1)
    h = 10.0 ** (-ph) / g1
    co2 = co2_kh(t_c) * pco2_atm
    return float(co2 * (k1 / h + 2.0 * k1 * k2 / (h * h)) + kw / h - h), float(co2 * (1 + k1 / h + k1 * k2 / h ** 2))


# --- pressure and the bubble point ------------------------------------------------------------------------------

def local_pressure_atm(z_m, surface_atm, gravity, rho=L.RHO0):
    return surface_atm + rho * gravity * np.asarray(z_m, dtype=float) / ATM


def gas_pressure_atm(t_c, ch4=0.0, co2_aq=0.0, n2=None, ar=None, o2=0.0):
    """The sum of the dissolved gases' partial pressures and the vapour (atm): CH4 and N2/Ar/O2 in mol/m3, dissolved
    CO2 in mol/kg."""
    t = np.asarray(t_c, dtype=float)
    total = L.vapour_pressure_atm(t) + np.asarray(ch4) / solubility('CH4', t) + np.asarray(co2_aq) / co2_kh(t)
    total = total + np.asarray(o2) / solubility('O2', t)
    if n2 is not None:
        total = total + np.asarray(n2) / solubility('N2', t)
    if ar is not None:
        total = total + np.asarray(ar) / solubility('Ar', t)
    return total


# --- sediment ebullition ----------------------------------------------------------------------------------------

def ebullition_share(production, depth_m, t_c, surface_atm, gravity, n2_atm, bottom_ch4=0.0, porosity=0.9,
                     decay_per_m=30.0):
    """Share of the sediment's methane production (mol/m2/s) that leaves as bubbles, in the steady porewater balance
    with production a exp(-b z) and N2 fixed at its partial pressure at the sediment surface (Bazhin 2010;
    Langenegger et al. 2019). Bubbles form where the methane's partial pressure reaches the local pressure less the
    N2 and the vapour."""
    if production <= 0.0:
        return 0.0
    p = float(local_pressure_atm(depth_m, surface_atm, gravity)) - float(L.vapour_pressure_atm(t_c))
    c_crit = float(solubility('CH4', t_c)) * max(p - n2_atm, 0.0)
    tortuosity = 1.0 - math.log(porosity ** 2)
    d_sed = float(diffusivity('CH4', t_c)) / tortuosity
    f_star = porosity * d_sed * decay_per_m * max(c_crit - bottom_ch4, 0.0)
    ratio = f_star / production
    if ratio >= 1.0:
        return 0.0
    target = 1.0 - ratio
    u = brentq(lambda x: (1.0 + x) * math.exp(-x) - target, 0.0, 200.0, xtol=1e-12)
    return math.exp(-u)


# --- a rising bubble --------------------------------------------------------------------------------------------

def rise_velocity(d, t_c, gravity):
    """Terminal speed (m/s) of a contaminated bubble of diameter d: Tomiyama's (1998) drag for contaminated systems,
    C_D = max(24/Re (1 + 0.15 Re^0.687), 8/3 Eo/(Eo + 4)), solved with v = (4 g d / 3 C_D)^0.5."""
    rho, mu = float(L.density(t_c)), float(L.viscosity(t_c))
    eo = gravity * rho * d * d / SURFACE_TENSION
    v = 0.2
    for _ in range(60):
        re = max(rho * v * d / mu, 1e-9)
        cd = max(24.0 / re * (1.0 + 0.15 * re ** 0.687), 8.0 / 3.0 * eo / (eo + 4.0))
        v_new = math.sqrt(4.0 * gravity * d / (3.0 * cd))
        if abs(v_new - v) < 1e-10:
            break
        v = 0.5 * (v + v_new)
    return v


def bubble_survival(depth_m, diameter_m, t_c, surface_atm, gravity, water, x_ch4=0.8, interface='rigid'):
    """Share of a bubble's methane that reaches the surface when it leaves the sediment at depth_m with diameter_m and
    methane mole fraction x_ch4 (the rest N2). water gives the dissolved CH4, N2 and O2 (mol/m3) it meets, uniform.
    interface 'rigid' takes Froessling's transfer for a contaminated surface, 'clean' Higbie's penetration,
    k = 2 (D v / pi d)^0.5, the faster bound. Returns the share and the rise time (s)."""
    tk = t_c + 273.15
    gases = ('CH4', 'N2', 'O2')
    beta = np.array([float(solubility(g, t_c)) for g in gases])
    diff = np.array([float(diffusivity(g, t_c)) for g in gases])
    sc = np.array([float(schmidt(g, t_c)) for g in gases])
    cw = np.array([water.get(g, 0.0) for g in gases])
    pv = float(L.vapour_pressure_atm(t_c)) * ATM
    rho, mu = float(L.density(t_c)), float(L.viscosity(t_c))

    def pressure(z):
        return surface_atm * ATM + rho * gravity * max(z, 0.0)

    r0 = 0.5 * diameter_m
    p0 = pressure(depth_m) + 2.0 * SURFACE_TENSION / r0
    n0 = (p0 - pv) * (4.0 / 3.0 * math.pi * r0 ** 3) / (GAS_R * tk)
    y0 = np.array([x_ch4 * n0, (1.0 - x_ch4) * n0, 0.0, depth_m])

    def radius(n, z):
        r = r0
        for _ in range(30):                              # the Laplace pressure ties the radius to itself
            vol = n * GAS_R * tk / (pressure(z) + 2.0 * SURFACE_TENSION / r - pv)
            r_new = (3.0 * vol / (4.0 * math.pi)) ** (1.0 / 3.0)
            if abs(r_new - r) < 1e-12:
                break
            r = r_new
        return r

    def rhs(_, y):
        n = np.maximum(y[:3], 0.0)
        tot = n.sum()
        if tot <= 1e-30:
            return np.zeros(4)
        r = radius(tot, y[3])
        d = 2.0 * r
        v = rise_velocity(d, t_c, gravity)
        re = rho * v * d / mu
        if interface == 'clean':
            k = 2.0 * np.sqrt(diff * v / (math.pi * d))
        else:
            k = diff * (2.0 + 0.6 * math.sqrt(re) * sc ** (1.0 / 3.0)) / d
        p_gas = (pressure(y[3]) + 2.0 * SURFACE_TENSION / r - pv) / ATM
        flux = 4.0 * math.pi * r * r * k * (beta * (n / tot) * p_gas - cw)
        return np.concatenate([-flux, [-v]])

    surface = lambda _, y: y[3]
    surface.terminal, surface.direction = True, -1
    gone = lambda _, y: y[:3].sum() - 1e-9 * n0
    gone.terminal = True
    sol = solve_ivp(rhs, (0.0, 1e6), y0, events=(surface, gone), rtol=1e-7, atol=1e-16 * max(n0, 1e-12),
                    max_step=max(depth_m / 0.05 / 400.0, 0.5))
    share = max(float(sol.y[0, -1]), 0.0) / y0[0] if sol.status == 1 and sol.t_events[0].size else 0.0
    return share, float(sol.t[-1])
