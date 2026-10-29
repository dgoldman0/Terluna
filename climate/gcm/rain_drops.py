"""PlaSim's falling rain evaporated from its drops' sizes and fall speeds at lunar gravity, in place of the Earth-tuned
constant gamma (register, 29 October 2026; the record is the section "PlaSim's rain from its drops" in README.md
beside this file).

PlaSim's mkrain (rainmod.f90) carries the precipitation each column forms down through its layers, turns it to snow
where the air is below freezing and back to rain where it is above, and with nevapprec 1 evaporates in each unsaturated
layer a share gamma of the layer's saturation deficit each step, at most the precipitation there. gamma, 0.01, was set
for Earth and knows nothing of the drops. With nevapprec 2 (this patch) the rain carries, beside its mass flux F, the
flux of its drops' number F_N, and its spectrum is Morrison's exponential one, N0 exp(-lambda D), whose slope follows
from the two: F_N/F = lambda^3 V_n(lambda) / (pi rho_w V_m(lambda)). Through each layer, in midpoint steps of at most
rdrop_step metres:

- evaporation, Morrison's equation as CM1's lunar build applies it: the rain loses k F per metre of fall,
      k = 2 rho D_v lambda^4 W (q_s - q) / (rho_w V_m AB),  W = F1R/lambda^2 + F2R (dens rho/mu)^(1/2) Sc^(1/3) T3,
  with AB PlaSim's own latent-heat correction, its drops' number falling with their mass (NSUBR), so that
  dF_N/ds = -k F_N V_m/V_n;
- collisions, Morrison's self-collection with the lunar breakup: dF_N/ds = -K dum (F/V_m)(F_N/V_n), K the lunar rate
  from the drops' collision kernels (2.919 m3/kg/s), dum = 1 below the breakup cutoff D_b and 2 - exp(KB (1/lambda -
  D_b)) above it, D_b and KB CM1's Earth values scaled as the rain limits scale them, by sqrt(9.81/g) times the column's
  surface pressure over 1013.25 hPa to the power 0.687; only on rain of at least 1e-8 kg/kg;
- V_m, V_n and the ventilation integral T3 from atmosphere/rain's lunar fall law (rain_fall_law.json, the table CM1's
  rain law reads), times Morrison's density factor dens = (RHOSU/rho)^0.54 and its square root, the mean speeds capped at
  9.1 (g/9.81)^(1/4) m/s times dens, lambda held between the lunar cap on the mean diameter, 2.8 mm sqrt(9.81/g), and
  20 um.

New rain enters each layer at its middle. Rain melted from snow there comes as new drops of mean diameter rdrop_melt;
rain formed as rain there comes as new drops of rdrop_warm, or with nrdrop_join 1 joins the rain already falling through
the layer, as accretion adds to drops without adding drops. What the layer evaporates is at most what would saturate it
over the step, as gamma's form is. Snow sublimates with gamma as before; mkrain's phase rule, which makes each layer's
precipitation all rain or all snow, is unchanged.

The scheme was tested on CM1's design boxes before it was built (rules and outcomes on the research drive,
session-tools/rain_physics): carrying the size down a column from the boxes' own mean air with CM1's own sizes for the
new rain gives CM1's column evaporation within 1.01-1.12 on the sea box, the land box and the inflow day. With every new
drop at one size, 0.465 mm, the mean of CM1's own size where the sea box's rain forms, it gives 1.17-1.31, outside the
test's 1.25 on the land box: the size new rain takes is open, and it waits on the snow work for the rain melted from
lunar snow. CM1 makes that rain from ice it sizes by Earth's limits.

`mkrain` is the patched mkrain's Python copy with nevapprec 2, columns as rows and levels top first, for the tests.
RAIN_DROP_PATCHES put the scheme into rainmod.f90 behind nevapprec; with nevapprec 1 the routine is as before.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from shared.constants import (CM1_GRAVITY, CM1_RAIN_BREAKUP_CUTOFF, CM1_RAIN_BREAKUP_RATE, CM1_RAIN_DENSITY_EXPONENT,
                              CM1_RAIN_DENSITY_REFERENCE_PRESSURE, CM1_RAIN_DENSITY_REFERENCE_TEMPERATURE,
                              CM1_RAIN_FALL_CAP, CM1_RAIN_MEAN_DIAMETER_CAP, CM1_RAIN_MEAN_DIAMETER_FLOOR,
                              CM1_RAIN_SELF_COLLECTION_FLOOR, CM1_RAIN_VENTILATION_F1, CM1_RAIN_VENTILATION_F2,
                              CM1_DRY_AIR_GAS_CONSTANT, CM1_VAPOUR_DIFFUSIVITY_COEFFICIENT,
                              CM1_VAPOUR_DIFFUSIVITY_EXPONENT, CM1_VISCOSITY_COEFFICIENT, CM1_VISCOSITY_SUTHERLAND,
                              CM1_WATER_DENSITY, STANDARD_ATMOSPHERE)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FALL_LAW = ROOT / 'atmosphere' / 'rain' / 'results' / 'rain_fall_law.json'
BREAKUP = ROOT / 'atmosphere' / 'rain' / 'results' / 'rain_breakup.json'
MARKER = 'Terluna: rain from its drops'
RHOSU = CM1_RAIN_DENSITY_REFERENCE_PRESSURE / (CM1_DRY_AIR_GAS_CONSTANT * CM1_RAIN_DENSITY_REFERENCE_TEMPERATURE)
INVERSE_POINTS = 2001             # the slope's inverse table, uniform in ln lambda across the scheme's range
EXPONENT_LIMIT = 50.0             # the most a step's exponent may change the rain's water or drops (e^50)
TINY_FLUX = 1.0e-30               # kg/m2/s: rain below it is gone
# The size new rain takes where CM1's sea box forms it: the mean of its rain-mass-weighted 1/lambda over the levels and
# outputs of days 29.5-59.0 where its rain's growth outran the rain arriving from above, weighted by that growth.
CM1_SEA_BIRTH_SIZE_M = 0.465e-3


def products(fall_law: Path = FALL_LAW, breakup: Path = BREAKUP) -> dict:
    """The lunar fall law's table and self-collection rate, and the breakup's density exponent."""
    law = json.loads(fall_law.read_text())
    lim = json.loads(breakup.read_text())['lunar_limits']['namelist']
    t = law['table']
    table = np.array([t['ln_mass_weighted'], t['ln_number_weighted'], t['ln_ventilation']], float)
    assert table.shape == (3, t['points'])
    return dict(x0=float(t['ln_lambda_first']), dx=float(t['ln_lambda_step']), table=table,
                ksc=float(law['namelist']['rain_selfcoll']), rhoexp=float(lim['brkup_rhoexp']))


PRODUCTS = products()
BIRTH_RULES = {
    'one_size': dict(nrdrop_join=0, rdrop_melt=CM1_SEA_BIRTH_SIZE_M, rdrop_warm=CM1_SEA_BIRTH_SIZE_M),
    'join': dict(nrdrop_join=1, rdrop_melt=CM1_SEA_BIRTH_SIZE_M, rdrop_warm=CM1_SEA_BIRTH_SIZE_M),
}
PARAMETERS = dict(BIRTH_RULES['one_size'], rdrop_ksc=PRODUCTS['ksc'], rdrop_rhoexp=PRODUCTS['rhoexp'],
                  rdrop_step=1000.0)


def lookup(lnlam, p=PRODUCTS):
    """The table's three moments at ln lambda (mass-weighted and number-weighted speeds in the reference air, the
    ventilation integral), linear in ln lambda and held at the ends, as CM1's TERLUNA_RAIN reads them."""
    x = (lnlam - p['x0']) / p['dx']
    n = p['table'].shape[1]
    i = min(max(int(x), 0), n - 2)
    w = min(max(x - i, 0.0), 1.0)
    t = p['table']
    return tuple(math.exp((1.0 - w) * t[m, i] + w * t[m, i + 1]) for m in range(3))


def limits(g: float) -> dict:
    """The lunar limits that do not depend on the column: the slope's range and the cap on the mean speeds."""
    scap = math.sqrt(CM1_GRAVITY / g)
    return dict(lnlam_min=math.log(1.0 / (CM1_RAIN_MEAN_DIAMETER_CAP * scap)),
                lnlam_max=math.log(1.0 / CM1_RAIN_MEAN_DIAMETER_FLOOR),
                vcap=CM1_RAIN_FALL_CAP * (g / CM1_GRAVITY) ** 0.25)


def inverse(g: float, p=PRODUCTS) -> dict:
    """ln(F_N/F) on INVERSE_POINTS slopes uniform in ln lambda over the scheme's range (the speeds capped), from which
    the slope of a rain of given fluxes is read back; the density factor cancels."""
    lim = limits(g)
    lnl = np.array([lim['lnlam_min'] + (lim['lnlam_max'] - lim['lnlam_min']) * i / (INVERSE_POINTS - 1)
                    for i in range(INVERSE_POINTS)])
    lnr = np.empty(INVERSE_POINTS)
    for i, x in enumerate(lnl):
        vm, vn, _ = lookup(x, p)
        lnr[i] = 3.0 * x + math.log(min(vn, lim['vcap'])) - math.log(min(vm, lim['vcap'])) \
            - math.log(math.pi * CM1_WATER_DENSITY)
    assert np.all(np.diff(lnr) > 0)
    return dict(lim, lnl=lnl, lnr=lnr, step=(lim['lnlam_max'] - lim['lnlam_min']) / (INVERSE_POINTS - 1))


def slope_of(r: float, inv: dict) -> float:
    """ln lambda for ln(F_N/F) = r: a search of the inverse table, linear between its points, held at its ends."""
    lnr, lnl = inv['lnr'], inv['lnl']
    if r <= lnr[0]:
        return float(lnl[0])
    if r >= lnr[-1]:
        return float(lnl[-1])
    lo, hi = 0, len(lnr) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if lnr[mid] <= r:
            lo = mid
        else:
            hi = mid
    return float(lnl[lo] + (r - lnr[lo]) / (lnr[hi] - lnr[lo]) * (lnl[hi] - lnl[lo]))


def number_per_flux(d: float, inv: dict) -> float:
    """F_N/F of rain whose mean diameter is d, its slope held within the scheme's range."""
    x = min(max(-math.log(d), inv['lnlam_min']), inv['lnlam_max'])
    i = min(int((x - inv['lnl'][0]) / inv['step']), INVERSE_POINTS - 2)
    w = (x - inv['lnl'][i]) / inv['step']
    return math.exp((1.0 - w) * inv['lnr'][i] + w * inv['lnr'][i + 1])


def rates(f, n, air, inv, prm, p=PRODUCTS):
    """The rain's loss of mass per metre, a, and of drops per metre, b, at fluxes f (kg/m2/s) and n (per m2 per s)."""
    if f <= TINY_FLUX or n <= 0.0:
        return 0.0, 0.0
    lnlam = slope_of(math.log(max(n, 1e-300) / f), inv)
    lam = math.exp(lnlam)
    vm0, vn0, t3 = lookup(lnlam, p)
    dens = (RHOSU / air['rho']) ** CM1_RAIN_DENSITY_EXPONENT
    vm = min(dens * vm0, inv['vcap'] * dens)
    vn = min(dens * vn0, inv['vcap'] * dens)
    w = CM1_RAIN_VENTILATION_F1 / lam**2 + CM1_RAIN_VENTILATION_F2 * math.sqrt(dens * air['rho'] / air['mu']) \
        * air['sc'] ** (1.0 / 3.0) * t3
    a = 0.0
    if air['q'] < air['qs']:
        a = 2.0 * air['rho'] * air['dv'] * lam**4 * w * (air['qs'] - air['q']) / (CM1_WATER_DENSITY * vm * air['ab'])
    d = 1.0 / lam
    dum = 1.0 if d < air['db'] else 2.0 - math.exp(air['kb'] * (d - air['db']))
    c = 0.0
    if f / vm / air['rho'] >= CM1_RAIN_SELF_COLLECTION_FLOOR:
        c = prm['rdrop_ksc'] * dum * (f / vm) / vn
    return a, a * vm / vn + c


def rdropevap(jlev, zq, zt, zqsat, zcor, zpra, zmelt, zprc, zprl, zprn, zdqdtl, zdtdt, s, c, prm, inv):
    """One layer of every column: the rain's evaporation and drops through it (rdropevap in Fortran, line for line).
    zprc, zprl (sigma-weighted mixing ratio per second), zprn (per m2 per s) and the layer's tendency columns change in
    place."""
    for j in range(zprc.size):
        mass = s['dp'][j] / c['ga']
        ftot = (zprc[j] + zprl[j]) * mass
        if ftot <= 0.0:
            zprn[j] = 0.0
            continue
        farr = min(max(zpra[j], 0.0), zprc[j] + zprl[j]) * mass
        bm = min(max(zmelt[j], 0.0) * mass, ftot - farr)
        bw = max(ftot - farr - bm, 0.0)
        if farr <= 0.0:
            zprn[j] = 0.0
        tv = zt[j] * (1.0 + (1.0 / c['rdbrv'] - 1.0) * zq[j])
        p = s['dp'][j] * s['sigma'][jlev]
        rho = p / (c['gascon'] * tv)
        dz = c['gascon'] * tv / c['ga'] * s['dsigma'][jlev] / s['sigma'][jlev]
        lcpe = c['ALV'] / (c['acpd'] * (1.0 + c['ADV'] * s['dq'][j, jlev]))
        ab = 1.0 + lcpe * c['ra2'] * (c['TMELT'] - c['ra4']) * zqsat[j] * zcor[j] / (zt[j] - c['ra4']) ** 2
        mu = CM1_VISCOSITY_COEFFICIENT * zt[j] ** 1.5 / (zt[j] + CM1_VISCOSITY_SUTHERLAND)
        dv = CM1_VAPOUR_DIFFUSIVITY_COEFFICIENT * zt[j] ** CM1_VAPOUR_DIFFUSIVITY_EXPONENT / p
        scol = math.sqrt(CM1_GRAVITY / c['ga']) * (s['dp'][j] / STANDARD_ATMOSPHERE) ** prm['rdrop_rhoexp']
        air = dict(q=zq[j], qs=zqsat[j], rho=rho, ab=ab, mu=mu, dv=dv, sc=mu / (rho * dv),
                   db=CM1_RAIN_BREAKUP_CUTOFF * scol, kb=CM1_RAIN_BREAKUP_RATE / scol)
        f, n = farr, zprn[j]
        if f > 0.0 and n <= 0.0:
            n = f * number_per_flux(prm['rdrop_warm'], inv)
        ratio = n / f if f > 0.0 else 0.0            # the rain's last drops per unit of its flux
        nstep = max(1, math.ceil(0.5 * dz / prm['rdrop_step']))
        ds = 0.5 * dz / nstep
        for half in (1, 2):
            if half == 2:                                  # the layer's new rain, at its middle
                if bm > 0.0:
                    n = n + bm * number_per_flux(prm['rdrop_melt'], inv)
                if bw > 0.0 and not (prm['nrdrop_join'] == 1 and f > 0.0):
                    n = n + bw * number_per_flux(prm['rdrop_warm'], inv)
                f = f + bm + bw
                if f > 0.0:
                    ratio = n / f
            for _ in range(nstep):
                if f <= TINY_FLUX:
                    f, n = 0.0, 0.0
                    break
                a1, b1 = rates(f, n, air, inv, prm)
                fm = f * math.exp(-min(0.5 * a1 * ds, EXPONENT_LIMIT))
                if fm <= TINY_FLUX:                       # gone within the half step
                    f, n = 0.0, 0.0
                    break
                nm = n * math.exp(-max(min(0.5 * b1 * ds, EXPONENT_LIMIT), -EXPONENT_LIMIT))
                a2, b2 = rates(fm, nm, air, inv, prm)
                f = f * math.exp(-min(a2 * ds, EXPONENT_LIMIT))
                n = n * math.exp(-max(min(b2 * ds, EXPONENT_LIMIT), -EXPONENT_LIMIT))
                if f > TINY_FLUX:
                    ratio = n / f
        evap = max(ftot - f, 0.0) / mass
        cap = max(zqsat[j] - zq[j], 0.0) * s['dsigma'][jlev] / c['deltsec2'] / ab
        if evap > cap:                                    # the layer saturates: the rest keeps its last size
            n = (ftot - cap * mass) * ratio
            evap = cap
        fc = zprc[j] / (zprc[j] + zprl[j])
        zprc[j] = max(zprc[j] - evap * fc, 0.0)
        zprl[j] = max(zprl[j] - evap * (1.0 - fc), 0.0)
        zprn[j] = n
        dq = evap / s['dsigma'][jlev]
        zdtdt[j] = zdtdt[j] - dq * lcpe
        zdqdtl[j] = zdqdtl[j] + dq


def mkrain(s: dict, c: dict, prm: dict = PARAMETERS) -> dict:
    """The patched mkrain with nevapprec 2 on columns (rows) and levels (top first): s holds dt and dq (with the ground
    after the levels), dtdt, dqdt, dp, sigma, dsigma and the precipitation formed by level (dprcl, dprscl, dprll,
    dprsll), and its tendencies are added in place; returns the precipitation at the ground (m/s of water) and the
    drops' number flux there."""
    n, nlev = s['dprcl'].shape
    inv = inverse(c['ga'])
    zprc, zprl, zprsc, zprsl, zprn = (np.zeros(n) for _ in range(5))
    zdtdt, zdqdts, zdqdtl = (np.zeros((n, nlev)) for _ in range(3))
    for jlev in range(nlev):
        zpra = zprc + zprl
        zprc = zprc + s['dprcl'][:, jlev] - s['dprscl'][:, jlev]
        zprl = zprl + s['dprll'][:, jlev] - s['dprsll'][:, jlev]
        zprsc = zprsc + s['dprscl'][:, jlev]
        zprsl = zprsl + s['dprsll'][:, jlev]
        zpsa = zprsc + zprsl
        # rain and snow by the layer's temperature, and in the lowest layer by the ground's as well (mkrain's where
        # and elsewhere: the snow melts where the air is above freezing, the rain freezes where it is not)
        conversions = [(s['dt'][:, jlev], False)] + ([(s['dt'][:, nlev], True)] if jlev == nlev - 1 else [])
        for tcol, add in conversions:
            warm = (zprsc + zprsl > 0.0) & (tcol > c['TMELT'])
            cold = ~warm & (zprc + zprl > 0.0) & (tcol <= c['TMELT'])
            zlcp = (c['ALS'] - c['ALV']) / (c['acpd'] * (1.0 + c['ADV'] * s['dq'][:, jlev]))
            base = zdtdt[:, jlev] if add else 0.0
            zdtdt[:, jlev] = np.where(warm, -(zprsc + zprsl) / s['dsigma'][jlev] * zlcp + base,
                                      np.where(cold, (zprc + zprl) / s['dsigma'][jlev] * zlcp + base, zdtdt[:, jlev]))
            zprc, zprl, zprsc, zprsl = (np.where(warm, zprsc + zprc, np.where(cold, 0.0, zprc)),
                                        np.where(warm, zprsl + zprl, np.where(cold, 0.0, zprl)),
                                        np.where(warm, 0.0, np.where(cold, zprc + zprsc, zprsc)),
                                        np.where(warm, 0.0, np.where(cold, zprl + zprsl, zprsl)))
        zq = s['dq'][:, jlev] + s['dqdt'][:, jlev] * c['deltsec2']
        zt = s['dt'][:, jlev] + s['dtdt'][:, jlev] * c['deltsec2']
        zqsat = c['rdbrv'] * c['ra1'] * np.exp(c['ra2'] * (zt - c['TMELT']) / (zt - c['ra4'])) \
            / (s['dp'] * s['sigma'][jlev])
        zqsat = np.minimum(zqsat, c['rdbrv'])
        zcor = 1.0 / (1.0 - (1.0 / c['rdbrv'] - 1.0) * zqsat)
        zqsat = zqsat * zcor
        zqsat = np.where(zqsat < 0.0, 0.0, zqsat)
        zmelt = np.maximum(zpsa - zprsc - zprsl, 0.0)
        rdropevap(jlev, zq, zt, zqsat, zcor, zpra, zmelt, zprc, zprl, zprn, zdqdtl[:, jlev], zdtdt[:, jlev], s, c,
                  prm, inv)
        for key in ('l', 'c'):
            zprs = zprsl if key == 'l' else zprsc
            sub = (zq < zqsat) & (zprs > 0.0)
            zlcpe = c['ALS'] / (c['acpd'] * (1.0 + c['ADV'] * s['dq'][:, jlev]))
            zdqdt = np.minimum(c['gamma'] * (zqsat - zq) * s['dsigma'][jlev] / c['deltsec2']
                               / (1.0 + zlcpe * c['ra2'] * (c['TMELT'] - c['ra4']) * zqsat * zcor
                                  / (zt - c['ra4']) ** 2), zprs)
            zdqdt = np.where(sub, zdqdt, 0.0)
            zprs[:] = zprs - zdqdt
            zdqdt = zdqdt / s['dsigma'][jlev]
            zdtdt[:, jlev] = np.where(sub, zdtdt[:, jlev] - zdqdt * zlcpe, zdtdt[:, jlev])
            zdqdts[:, jlev] = np.where(sub, zdqdts[:, jlev] + zdqdt, zdqdts[:, jlev])
    s['dqdt'][:, :nlev] += zdqdtl + zdqdts
    s['dtdt'][:, :nlev] += zdtdt
    to_ms = s['dp'] / c['ga'] / 1000.0
    return dict(dprc=(zprc + zprsc) * to_ms, dprl=(zprl + zprsl) * to_ms, dprs=(zprsc + zprsl) * to_ms,
                number=zprn, rain_evaporation=zdqdtl, snow_sublimation=zdqdts)


# ----------------------------------------------------------------------------------------------- the Fortran patch

def _fortran_list(values, per_line=5) -> str:
    rows = [', '.join(values[i:i + per_line]) for i in range(0, len(values), per_line)]
    return ' &\n     &  ' + ', &\n     &  '.join(rows) + ' &\n     & '


def declarations(prm: dict = PARAMETERS, p=PRODUCTS) -> str:
    """rainmod's new variables: the namelist's settings, the constants as parameters, the fall law's table and the
    slope's inverse table."""
    n = p['table'].shape[1]
    values = [repr(float(p['table'][m, i])) for i in range(n) for m in range(3)]
    return (
        f"      integer :: nrdrop_join = {prm['nrdrop_join']} ! {MARKER} (nevapprec 2): 1 joins rain formed as rain to "
        "the rain falling\n"
        f"      real :: rdrop_melt = {prm['rdrop_melt']!r} ! mean diameter (m) of the drops of rain melted from snow\n"
        f"      real :: rdrop_warm = {prm['rdrop_warm']!r} ! mean diameter (m) of the drops of rain formed as rain\n"
        f"      real :: rdrop_ksc = {prm['rdrop_ksc']!r} ! rain's self-collection at lunar gravity (m3/kg/s)\n"
        f"      real :: rdrop_rhoexp = {prm['rdrop_rhoexp']!r} ! the breakup's scale: exponent on the air's density "
        "ratio\n"
        f"      real :: rdrop_step = {prm['rdrop_step']!r} ! the longest midpoint step through a layer (m)\n"
        f"      integer, parameter :: nrdtab = {n} ! the lunar fall law (atmosphere/rain, rain_fall_law.json)\n"
        f"      real, parameter :: rdx0 = {p['x0']!r}, rddx = {p['dx']!r} ! its ln lambda (lambda per m)\n"
        "      real, parameter :: rdtab(3,nrdtab) = reshape((/" + _fortran_list(values) + "/), (/3,nrdtab/))\n"
        f"      real, parameter :: rdf1 = {CM1_RAIN_VENTILATION_F1!r}, rdf2 = {CM1_RAIN_VENTILATION_F2!r} ! Morrison's "
        "ventilation\n"
        f"      real, parameter :: rdrhow = {CM1_WATER_DENSITY!r}, rdrhosu = {RHOSU!r}, rdexpd = "
        f"{CM1_RAIN_DENSITY_EXPONENT!r}\n"
        f"      real, parameter :: rddv0 = {CM1_VAPOUR_DIFFUSIVITY_COEFFICIENT!r}, rddvx = "
        f"{CM1_VAPOUR_DIFFUSIVITY_EXPONENT!r}, rdmu0 = {CM1_VISCOSITY_COEFFICIENT!r}, rdmus = "
        f"{CM1_VISCOSITY_SUTHERLAND!r}\n"
        f"      real, parameter :: rdbreak = {CM1_RAIN_BREAKUP_CUTOFF!r}, rdkbreak = {CM1_RAIN_BREAKUP_RATE!r} ! CM1's "
        "breakup, Earth's\n"
        f"      real, parameter :: rdcap0 = {CM1_RAIN_MEAN_DIAMETER_CAP!r}, rdfloor0 = "
        f"{CM1_RAIN_MEAN_DIAMETER_FLOOR!r}, rdvcap0 = {CM1_RAIN_FALL_CAP!r}\n"
        f"      real, parameter :: rdgref = {CM1_GRAVITY!r}, rdpref = {STANDARD_ATMOSPHERE!r}, rdcfloor = "
        f"{CM1_RAIN_SELF_COLLECTION_FLOOR!r}\n"
        f"      integer, parameter :: nrdinv = {INVERSE_POINTS} ! the slope's inverse table\n"
        f"      real, parameter :: rdxlim = {EXPONENT_LIMIT!r}, rdtiny = {TINY_FLUX!r} ! a step's largest exponent; "
        "rain gone\n"
        "      real :: rdinvl(nrdinv), rdinvr(nrdinv), rdlmin, rdlmax, rdlstep, rdvcap\n"
        "      integer :: nrdinit = 0\n")


NAMELIST_NAMES = ('nrdrop_join', 'rdrop_melt', 'rdrop_warm', 'rdrop_ksc', 'rdrop_rhoexp', 'rdrop_step')

ROUTINES = """
!     ====================
!     SUBROUTINE RDROPINIT
!     ====================

      subroutine rdropinit
      use rainmod
!
!     Terluna: rain from its drops: the slope's inverse table, ln(F_N/F) on nrdinv slopes uniform in ln lambda across
!     the lunar range (the cap on the mean diameter rdcap0 sqrt(rdgref/g), the floor rdfloor0), the mean speeds capped
!     at rdvcap0 (g/rdgref)**(1/4); the density factor cancels in the ratio
!
      integer :: ji
      real :: zvm,zvn,zv3
      rdlmin=LOG(1./(rdcap0*SQRT(rdgref/ga)))
      rdlmax=LOG(1./rdfloor0)
      rdlstep=(rdlmax-rdlmin)/REAL(nrdinv-1)
      rdvcap=rdvcap0*(ga/rdgref)**0.25
      do ji=1,nrdinv
       rdinvl(ji)=rdlmin+(rdlmax-rdlmin)*REAL(ji-1)/REAL(nrdinv-1)
       call rdroptab(rdinvl(ji),zvm,zvn,zv3)
       rdinvr(ji)=3.*rdinvl(ji)+LOG(MIN(zvn,rdvcap))-LOG(MIN(zvm,rdvcap))-LOG(4.*ATAN(1.)*rdrhow)
      enddo
      nrdinit=1
      return
      end subroutine rdropinit

!     ===================
!     SUBROUTINE RDROPTAB
!     ===================

      subroutine rdroptab(plnl,pvm,pvn,pv3)
      use rainmod
!
!     Terluna: the lunar fall law's mass-weighted and number-weighted mean speeds in the reference air and its
!     ventilation integral at ln lambda plnl, linear in ln lambda and held at the table's ends (CM1's TERLUNA_RAIN)
!
      real :: plnl,pvm,pvn,pv3
      real :: zx,zw
      integer :: ji
      zx=(plnl-rdx0)/rddx
      ji=MIN(MAX(INT(zx),0),nrdtab-2)
      zw=MIN(MAX(zx-REAL(ji),0.),1.)
      pvm=EXP((1.-zw)*rdtab(1,ji+1)+zw*rdtab(1,ji+2))
      pvn=EXP((1.-zw)*rdtab(2,ji+1)+zw*rdtab(2,ji+2))
      pv3=EXP((1.-zw)*rdtab(3,ji+1)+zw*rdtab(3,ji+2))
      return
      end subroutine rdroptab

!     =================
!     FUNCTION RDROPSLOPE
!     =================

      real function rdropslope(pr)
      use rainmod
!
!     Terluna: ln lambda of rain whose ln(F_N/F) is pr: a search of the inverse table, linear between its points and
!     held at its ends
!
      real :: pr
      integer :: jlo,jhi,jmid
      if(pr <= rdinvr(1)) then
       rdropslope=rdinvl(1)
      elseif(pr >= rdinvr(nrdinv)) then
       rdropslope=rdinvl(nrdinv)
      else
       jlo=1
       jhi=nrdinv
       do while(jhi-jlo > 1)
        jmid=(jlo+jhi)/2
        if(rdinvr(jmid) <= pr) then
         jlo=jmid
        else
         jhi=jmid
        endif
       enddo
       rdropslope=rdinvl(jlo)+(pr-rdinvr(jlo))/(rdinvr(jhi)-rdinvr(jlo))*(rdinvl(jhi)-rdinvl(jlo))
      endif
      return
      end function rdropslope

!     ================
!     FUNCTION RDROPNPER
!     ================

      real function rdropnper(pd)
      use rainmod
!
!     Terluna: F_N/F of rain whose mean diameter is pd (m), its slope held within the scheme's range
!
      real :: pd,zx,zw
      integer :: ji
      zx=MIN(MAX(-LOG(pd),rdlmin),rdlmax)
      ji=MIN(INT((zx-rdinvl(1))/rdlstep),nrdinv-2)+1
      zw=(zx-rdinvl(ji))/rdlstep
      rdropnper=EXP((1.-zw)*rdinvr(ji)+zw*rdinvr(ji+1))
      return
      end function rdropnper

!     =====================
!     SUBROUTINE RDROPRATES
!     =====================

      subroutine rdroprates(pf,pn,pq,pqs,prho,pab,pmu,pdv,psc,pdb,pkb,pa,pb)
      use rainmod
!
!     Terluna: the rain's loss of water per metre of fall pa and of drops per metre pb, at mass flux pf (kg/m2/s) and
!     number flux pn (per m2 per s): Morrison's evaporation, with the drops' number falling with their mass, and his
!     self-collection with breakup (rain_drops.rates)
!
      real :: pf,pn,pq,pqs,prho,pab,pmu,pdv,psc,pdb,pkb,pa,pb
      real :: zlnl,zlam,zvm0,zvn0,zv3,zdens,zvm,zvn,zw,zd,zdum,zc
      real :: rdropslope
      pa=0.
      pb=0.
      if(pf <= rdtiny .or. pn <= 0.) return
      zlnl=rdropslope(LOG(MAX(pn,1.e-300)/pf))
      zlam=EXP(zlnl)
      call rdroptab(zlnl,zvm0,zvn0,zv3)
      zdens=(rdrhosu/prho)**rdexpd
      zvm=MIN(zdens*zvm0,rdvcap*zdens)
      zvn=MIN(zdens*zvn0,rdvcap*zdens)
      zw=rdf1/zlam**2+rdf2*SQRT(zdens*prho/pmu)*psc**(1./3.)*zv3
      if(pq < pqs) pa=2.*prho*pdv*zlam**4*zw*(pqs-pq)/(rdrhow*zvm*pab)
      zd=1./zlam
      if(zd < pdb) then
       zdum=1.
      else
       zdum=2.-EXP(pkb*(zd-pdb))
      endif
      zc=0.
      if(pf/zvm/prho >= rdcfloor) zc=rdrop_ksc*zdum*(pf/zvm)/zvn
      pb=pa*zvm/zvn+zc
      return
      end subroutine rdroprates

!     ====================
!     SUBROUTINE RDROPEVAP
!     ====================

      subroutine rdropevap(jlev,pq,pt,pqsat,pcor,ppra,pmelt,pprc,pprl,pprn,pdqdt,pdtdt)
      use rainmod
!
!     Terluna: rain from its drops (nevapprec 2; climate/gcm/rain_drops.py, its rdropevap line for line). Through
!     layer jlev of every column the rain arriving from above (ppra, its drops pprn) falls half the layer; the layer's
!     new rain joins it at the middle, pmelt of it melted from snow as drops of rdrop_melt and the rest formed as rain
!     as drops of rdrop_warm or, with nrdrop_join 1, joining the rain falling; and the whole falls the lower half, in
!     midpoint steps of at most rdrop_step, losing water and drops (rdroprates). At most what would saturate the layer
!     over the step evaporates, shared between the convective and the large-scale rain by their fluxes.
!
      integer :: jlev
      real :: pq(NHOR),pt(NHOR),pqsat(NHOR),pcor(NHOR),ppra(NHOR),pmelt(NHOR)
      real :: pprc(NHOR),pprl(NHOR),pprn(NHOR),pdqdt(NHOR),pdtdt(NHOR)
      integer :: jhor,jhalf,jsub,nsub
      real :: zmass,zftot,zfarr,zbm,zbw,ztv,zp,zrho,zdz,zlcpe,zab,zmu,zdv,zsc,zscol,zdb,zkb
      real :: zf,zn,zds,za1,zb1,za2,zb2,zfm,znm,zevap,zcap,zfc,zdq,zrat
      real :: rdropnper
      if(nrdinit == 0) call rdropinit
      do jhor=1,NHOR
       zmass=dp(jhor)/ga
       zftot=(pprc(jhor)+pprl(jhor))*zmass
       if(zftot <= 0.) then
        pprn(jhor)=0.
        cycle
       endif
       zfarr=MIN(MAX(ppra(jhor),0.),pprc(jhor)+pprl(jhor))*zmass
       zbm=MIN(MAX(pmelt(jhor),0.)*zmass,zftot-zfarr)
       zbw=MAX(zftot-zfarr-zbm,0.)
       if(zfarr <= 0.) pprn(jhor)=0.
       ztv=pt(jhor)*(1.+(1./rdbrv-1.)*pq(jhor))
       zp=dp(jhor)*sigma(jlev)
       zrho=zp/(gascon*ztv)
       zdz=gascon*ztv/ga*dsigma(jlev)/sigma(jlev)
       zlcpe=ALV/(acpd*(1.+ADV*dq(jhor,jlev)))
       zab=1.+zlcpe*ra2*(TMELT-ra4)*pqsat(jhor)*pcor(jhor)/(pt(jhor)-ra4)**2
       zmu=rdmu0*pt(jhor)**1.5/(pt(jhor)+rdmus)
       zdv=rddv0*pt(jhor)**rddvx/zp
       zsc=zmu/(zrho*zdv)
       zscol=SQRT(rdgref/ga)*(dp(jhor)/rdpref)**rdrop_rhoexp
       zdb=rdbreak*zscol
       zkb=rdkbreak/zscol
       zf=zfarr
       zn=pprn(jhor)
       if(zf > 0. .and. zn <= 0.) zn=zf*rdropnper(rdrop_warm)
       zrat=0.
       if(zf > 0.) zrat=zn/zf
       nsub=MAX(1,CEILING(0.5*zdz/rdrop_step))
       zds=0.5*zdz/REAL(nsub)
       do jhalf=1,2
        if(jhalf == 2) then
         if(zbm > 0.) zn=zn+zbm*rdropnper(rdrop_melt)
         if(zbw > 0. .and. .not.(nrdrop_join == 1 .and. zf > 0.)) zn=zn+zbw*rdropnper(rdrop_warm)
         zf=zf+zbm+zbw
         if(zf > 0.) zrat=zn/zf
        endif
        do jsub=1,nsub
         if(zf <= rdtiny) then
          zf=0.
          zn=0.
          exit
         endif
         call rdroprates(zf,zn,pq(jhor),pqsat(jhor),zrho,zab,zmu,zdv,zsc,zdb,zkb,za1,zb1)
         zfm=zf*EXP(-MIN(0.5*za1*zds,rdxlim))
         if(zfm <= rdtiny) then
          zf=0.
          zn=0.
          exit
         endif
         znm=zn*EXP(-MAX(MIN(0.5*zb1*zds,rdxlim),-rdxlim))
         call rdroprates(zfm,znm,pq(jhor),pqsat(jhor),zrho,zab,zmu,zdv,zsc,zdb,zkb,za2,zb2)
         zf=zf*EXP(-MIN(za2*zds,rdxlim))
         zn=zn*EXP(-MAX(MIN(zb2*zds,rdxlim),-rdxlim))
         if(zf > rdtiny) zrat=zn/zf
        enddo
       enddo
       zevap=MAX(zftot-zf,0.)/zmass
       zcap=MAX(pqsat(jhor)-pq(jhor),0.)*dsigma(jlev)/deltsec2/zab
       if(zevap > zcap) then   ! the layer saturates: the rest keeps its last size
        zn=(zftot-zcap*zmass)*zrat
        zevap=zcap
       endif
       zfc=pprc(jhor)/(pprc(jhor)+pprl(jhor))
       pprc(jhor)=MAX(pprc(jhor)-zevap*zfc,0.)
       pprl(jhor)=MAX(pprl(jhor)-zevap*(1.-zfc),0.)
       pprn(jhor)=zn
       zdq=zevap/dsigma(jlev)
       pdtdt(jhor)=pdtdt(jhor)-zdq*zlcpe
       pdqdt(jhor)=pdqdt(jhor)+zdq
      enddo
      return
      end subroutine rdropevap
"""

EVAPORATION_BLOCK = f"""!*     {MARKER} (nevapprec 2): the rain's evaporation from its drops' sizes and fall speeds
!      (rdropevap); the snow sublimates with gamma as below
!
       if(nevapprec == 2) then
         zq(:)=dq(:,jlev)+dqdt(:,jlev)*deltsec2
         zt(:)=dt(:,jlev)+dtdt(:,jlev)*deltsec2
         zqsat(:)=rdbrv*ra1*exp(ra2*(zt(:)-TMELT)/(zt(:)-ra4))          &
     &         /(dp(:)*sigma(jlev))
         zqsat(:)=AMIN1(zqsat(:),rdbrv)
         zcor(:)=1./(1.-(1./rdbrv-1.)*zqsat(:))
         zqsat(:)=zqsat(:)*zcor(:)
         where (zqsat(:) < 0.0)
           zqsat(:) = 0.0
         endwhere
         zmelt(:)=MAX(zpsa(:)-zprsc(:)-zprsl(:),0.)
         call rdropevap(jlev,zq,zt,zqsat,zcor,zpra,zmelt,zprc,zprl,zprn,zdqdtl(:,jlev),zdtdt(:,jlev))
         where((zq(:).lt.zqsat(:)).AND.zprsl(:) > zeps)
           zlcpe(:)=ALS/(acpd*(1.+ADV*dq(:,jlev)))
           zdqdt(:)=AMIN1(gamma*(zqsat(:)-zq(:))*dsigma(jlev)/deltsec2  &
     &                   /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                 &
     &                   *zqsat(:)*zcor(:)/(zt(:)-ra4)**2)              &
     &                   ,zprsl(:))
           zprsl(:)=zprsl(:)-zdqdt(:)
           zdqdt(:)=zdqdt(:)/dsigma(jlev)
           zdtdt(:,jlev)=zdtdt(:,jlev)-zdqdt(:)*zlcpe(:)
           zdqdts(:,jlev)=zdqdts(:,jlev)+zdqdt(:)
         endwhere
         where((zq(:).lt.zqsat(:)).AND.zprsc(:) > zeps)
           zlcpe(:)=ALS/(acpd*(1.+ADV*dq(:,jlev)))
           zdqdt(:)=AMIN1(gamma*(zqsat(:)-zq(:))*dsigma(jlev)/deltsec2  &
     &                   /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                 &
     &                   *zqsat(:)*zcor(:)/(zt(:)-ra4)**2)              &
     &                   ,zprsc(:))
           zprsc(:)=zprsc(:)-zdqdt(:)
           zdqdt(:)=zdqdt(:)/dsigma(jlev)
           zdtdt(:,jlev)=zdtdt(:,jlev)-zdqdt(:)*zlcpe(:)
           zdqdts(:,jlev)=zdqdts(:,jlev)+zdqdt(:)
         endwhere
       endif
!
"""

NEVAPPREC_LINE = '      integer :: nevapprec= 1 ! switch for evaporation of precip.(1/0=yes/no)\n'
RESET_LINE = '!     reset defaults (according to general setup... tuning)\n'
BROADCAST_LINE = '      call mpbci(nevapprec)\n'
DECLARE_LINE = '      real :: zdqdtl(NHOR,NLEV)         ! evaporation rate\n'
LOOP_LINES = '      do jlev=1,NLEV\n       zprc(:)=zprc(:)+dprcl(:,jlev)-dprscl(:,jlev)\n'
SNOW_LINE = '       zprsl(:)=zprsl(:)+dprsll(:,jlev)\n'
EVAPORATION_LINE = '!*     evaporation/sublimation of rain/snow\n'
END_LINE = '      end subroutine mkrain\n'


def edits(prm: dict = PARAMETERS) -> list:
    """The edits of rainmod.f90, each (old, new) on a line no other Terluna patch touches."""
    return [
        (NEVAPPREC_LINE, NEVAPPREC_LINE + declarations(prm)),
        (RESET_LINE, f"      namelist/rainmod_nl/{','.join(NAMELIST_NAMES)}   ! {MARKER}\n" + RESET_LINE),
        (BROADCAST_LINE, BROADCAST_LINE + '      call mpbci(nrdrop_join)\n'
         + ''.join(f'      call mpbcr({name})\n' for name in NAMELIST_NAMES[1:])),
        (DECLARE_LINE, DECLARE_LINE + f'      real :: zpra(NHOR),zpsa(NHOR),zmelt(NHOR),zprn(NHOR) ! {MARKER} '
                                      '(nevapprec 2)\n'),
        (LOOP_LINES, f'      zprn(:)=0.   ! {MARKER}: the drops\' number flux (per m2 per s), none above the top\n'
                     '      do jlev=1,NLEV\n'
                     '       if(nevapprec == 2) zpra(:)=zprc(:)+zprl(:)   ! the rain arriving from above\n'
                     '       zprc(:)=zprc(:)+dprcl(:,jlev)-dprscl(:,jlev)\n'),
        (SNOW_LINE, SNOW_LINE + '       if(nevapprec == 2) zpsa(:)=zprsc(:)+zprsl(:)   ! the snow before the '
                                "layer's melting\n"),
        (EVAPORATION_LINE, EVAPORATION_BLOCK + EVAPORATION_LINE),
        (END_LINE, END_LINE + ROUTINES),
    ]


RAIN_DROP_EDITS = edits()
RAIN_DROP_PATCHES = [('rainmod.f90', MARKER, RAIN_DROP_EDITS)]


def setting(rule: str) -> dict:
    """The scheme as a run's configuration records it, for a birth rule of BIRTH_RULES."""
    if rule not in BIRTH_RULES:
        raise ValueError(f'the rain from its drops takes a birth rule of {sorted(BIRTH_RULES)}, not {rule!r}')
    prm = dict(PARAMETERS, **BIRTH_RULES[rule])
    return dict(rule=rule, nevapprec=2, namelist={k: prm[k] for k in NAMELIST_NAMES},
                fall_law=str(FALL_LAW.relative_to(ROOT)), breakup=str(BREAKUP.relative_to(ROOT)),
                snow='sublimates with gamma, as PlaSim takes it',
                code='climate/gcm/rain_drops.py (mkrain, rdropevap, rates; rdropevap and rdroprates in Fortran)')
