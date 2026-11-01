"""PlaSim's vertical diffusion and dry adjustment, one column at a time in Python, and a nonlocal boundary layer for it
(27 October 2026; the record is the section "A nonlocal boundary layer for the Moon" in README.md beside this file).

PlaSim mixes its air in two routines. Its vertical diffusion (fluxmod.f90, vdiff) gives heat, vapour and momentum the
local diffusivity K = l^2 |dV/dz| F(Ri), with Blackadar's mixing length l = lambda kappa z / (lambda + kappa z) (lambda
vdiff_lamm for momentum and sqrt(7.5) vdiff_lamm for heat and vapour, heights counted at 250 K), the Louis-ECHAM
stability functions and a least wind-speed difference between levels, and solves the diffusion implicitly in sigma,
which keeps the column's sum of sigma-thickness times temperature, vapour and wind. It heats the air by the kinetic
energy the diffusion of momentum dissipates. Its dry adjustment (rainmod.f90, mkdca), later in the step, mixes potential
temperature and vapour evenly through any layers where the dry potential temperature falls with height. `vdiff` and
`dca` are those routines on one column, step for step.

The nonlocal boundary layer (`nonlocal_column`) follows Holtslag and Boville (1993) as CAM3 implements it (Collins and
others 2004, section 4.11.2; CAM's hb_diff.F90), with the explicit entrainment flux of Hong, Noh and Dudhia (2006, the YSU
scheme). Its depth h comes from the bulk Richardson number between the lowest level and the levels above (Vogelezang
and Holtslag's form, with B u*^2 in the shear), with the thermals' excess b (w'theta_v')_0 / w_m in unstable air; with an
explicit entrainment flux h is the inversion's base, where that number reaches zero (YSU), and 0.3 in stable air
(CAM3). Inside h the diffusivity is kappa w z (1 - z/h)^2, w the velocity scale from u* and w*: w_m = (u*^3 + 0.6
w*^3)^(1/3) for momentum and w_m / Pr for heat and vapour, the Prandtl number from surface-layer similarity at 0.1 h;
in stable air u* / phi_h(z/L). A countergradient flux K gamma, gamma = a w* (w'c')_0 / (w_m^2 h), carries heat and
vapour, and an entrainment flux at h, -0.15 times the surface buoyancy flux in free convection, carries each across
the inversion in proportion to its jump between the levels bracketing h, falling as (z/h)^3 below. Inside h the
diffusivity is the larger of the profile's and the local one; above h the local scheme runs alone. The nonlocal fluxes
enter explicitly before vdiff's implicit solve, as in Holtslag and Boville's appendix B, so the column's heat, water and
momentum are kept exactly as vdiff keeps them; a vapour profile the explicit step would leave negative keeps its own.

Arrays run top first, as PlaSim's levels do: layer k lies between half levels sigmah[k-1] (above; 0 at the top) and
sigmah[k] (below; 1 at the ground), and the interface between layers k and k+1 is sigmah[k], k < n - 1.
"""
from __future__ import annotations

import math

import numpy as np

from shared.constants import (EXOPLASIM_DIFFUSION_HEIGHT_TEMPERATURE, EXOPLASIM_DIFFUSION_STABILITY,
                              EXOPLASIM_MIXING_LENGTH, EXOPLASIM_VON_KARMAN, EXOPLASIM_WIND_FLOOR, STANDARD_GRAVITY)

KARMAN = EXOPLASIM_VON_KARMAN                    # fluxmod's vonkarman, 0.4
LAMBDA_EARTH_M = EXOPLASIM_MIXING_LENGTH         # vdiff_lamm, 160 m: Blackadar's length for momentum, set for Earth
WIND_FLOOR = EXOPLASIM_WIND_FLOOR                # zumin, 1 m/s
STABILITY = EXOPLASIM_DIFFUSION_STABILITY        # vdiff_b = vdiff_c = vdiff_d = 5
HEIGHT_TEMPERATURE = EXOPLASIM_DIFFUSION_HEIGHT_TEMPERATURE   # ztscal, 250 K


def lunar_length(gravity: float) -> float:
    """vdiff_lamm stretched by the ratio of the gravities, the rule climate/crm/cm1_run.py applies to CM1's l_inf: a
    column stretched by g_E / g keeps its velocities and pressures and lengthens its lengths by that ratio (a scaling
    assumption the CM1 boxes share, tested here as a sensitivity)."""
    return LAMBDA_EARTH_M * STANDARD_GRAVITY / gravity


# The nonlocal scheme's constants, each dimensionless, with its source (Holtslag and Boville 1993, HB93; CAM3's
# technical description, Collins and others 2004, and CAM's hb_diff.F90; Hong, Noh and Dudhia 2006, YSU).
PARAMETERS = dict(
    b_excess=8.5,           # thermals' excess b (w'theta_v')_0 / w_m: HB93 (3.12), (A14), b = 0.85 d; CAM fak
    a_countergradient=7.2,  # gamma = a w* (w'c')_0 / (w_m^2 h): HB93 (3.10), (A14), a = 0.85^2 d, d = 10; CAM fakn
    eps_surface=0.1,        # the surface layer's share of h, where Pr is matched: HB93 (A13); CAM sffrac
    beta_m=15.0,            # phi_m = (1 - 15 z/L)^(-1/3), unstable (Dyer 1974): HB93 (A8); CAM betam
    beta_h=15.0,            # phi_h = (1 - 15 z/L)^(-1/2), unstable: HB93 (A6); CAM betah
    beta_s=5.0,             # phi_h = 1 + 5 z/L (z/L <= 1), 5 + z/L above, stable: HB93 (A3), (A5); CAM betas
    ricr_stable=0.3,        # critical bulk Richardson number, stable or neutral surface flux: CAM3 (4.471); CAM ricr
    ricr_unstable=0.0,      # ..., unstable: YSU (A12), h the inversion's base where the entrainment flux is explicit
    b_ustar=100.0,          # B u*^2 in the bulk number's shear: Vogelezang and Holtslag (1996) as CAM3 (4.471); CAM fac
    a_entrain=0.15,         # (w'theta_v')_h = -0.15 (theta_v / g) w_e^3 / h: YSU (A9), Noh and others (2003)
    c_entrain_ustar=5.0,    # w_e^3 = w*^3 + 5 u*^3: YSU (A8), (A9)
)
TINY = 1.0e-36            # CAM's least squared wind difference in the bulk number (hb_diff.F90, tiny): a numerical guard


def plasim_constants(c: dict) -> dict:
    """The constants the routines read: gravity ga, gascon, akap (gascon / acpd unless given), rdbrv, acpd and ADV."""
    out = dict(ga=float(c['ga']), gascon=float(c['gascon']), rdbrv=float(c['rdbrv']), acpd=float(c['acpd']),
               ADV=float(c.get('ADV', 0.0)))
    out['akap'] = float(c.get('akap', out['gascon'] / out['acpd']))
    return out


def full_sigma(sigmah) -> np.ndarray:
    """PlaSim's full levels, midway in sigma between the half levels that bound each layer (plasim.f90, initpm)."""
    sigmah = np.asarray(sigmah, float)
    return 0.5 * (np.r_[0.0, sigmah[:-1]] + sigmah)


def heights(t, q, sigma, sigmah, g: float, gascon: float, rdbrv: float):
    """Heights (m) above the ground of the full levels and of each layer's lower half level, from the layers' virtual
    temperatures (hydrostatic in sigma; the top layer's upper boundary, sigma 0, is left out)."""
    tv = np.asarray(t, float) * (1.0 + (1.0 / rdbrv - 1.0) * np.asarray(q, float))
    sigma, sigmah = np.asarray(sigma, float), np.asarray(sigmah, float)
    n = tv.size
    zh = np.zeros(n)
    for k in range(n - 1, 0, -1):
        zh[k - 1] = zh[k] + gascon * tv[k] / g * math.log(sigmah[k] / sigmah[k - 1])
    zf = zh + gascon * tv / g * np.log(sigmah / sigma)
    return zf, zh


def local_coefficients(t, q, u, v, sigma, sigmah, dtime: float, c: dict, lamm: float = LAMBDA_EARTH_M,
                       umin: float = WIND_FLOOR, b: float = STABILITY, cc: float = STABILITY, d: float = STABILITY):
    """vdiff's coefficients on one column: for each interface k (between layers k and k+1) the sigma-form diffusivities
    for momentum and for heat and vapour (dimensionless, dtime K (d sigma/dz)^2 / (sigma_k+1 - sigma_k)), and the
    interface's d sigma / dz (1/m, g sigmah / (R T_h) with vdiff's mean temperature T_h), as fluxmod.f90 computes them."""
    c = plasim_constants(c)
    ga, gascon, akap, rdbrv = c['ga'], c['gascon'], c['akap'], c['rdbrv']
    t, q, u, v = (np.asarray(a, float) for a in (t, q, u, v))
    sigma, sigmah = np.asarray(sigma, float), np.asarray(sigmah, float)
    dsigma = np.diff(np.r_[0.0, sigmah])
    n = t.size
    zlamh = lamm * math.sqrt(3.0 * d * 0.5)
    zkonst1 = ga * dtime / gascon
    zkonst2 = zkonst1 * ga / gascon
    zkonst3 = zkonst2 * ga / gascon
    zskap = sigma ** akap
    zskaph = sigmah ** akap
    zabsu = np.sqrt(u ** 2 + v ** 2)
    km, kh, to_z = np.zeros(n - 1), np.zeros(n - 1), np.zeros(n - 1)
    for k in range(n - 1):
        jp = k + 1
        zzlev = -gascon * HEIGHT_TEMPERATURE * math.log(sigmah[k]) / ga
        zdz = gascon * HEIGHT_TEMPERATURE * math.log(sigma[jp] / sigma[k]) / ga
        zmixh = zlamh * KARMAN * zzlev / (zlamh + KARMAN * zzlev)
        zmixm = lamm * KARMAN * zzlev / (lamm + KARMAN * zzlev)
        zmixh2, zmixm2 = zmixh * zmixh, zmixm * zmixm
        zrdsig = 1.0 / (sigma[jp] - sigma[k])
        zfac = zkonst3 * sigmah[k] * sigmah[k] * sigmah[k] * zrdsig
        zzkfac = math.sqrt(((((zzlev + zdz) / zzlev) ** (1.0 / 3.0) - 1.0) / zdz) ** 3 / zzlev)
        zrthl = (dsigma[jp] + dsigma[k]) / (t[k] * dsigma[jp] + t[jp] * dsigma[k])
        zdvds = max(umin, abs(zabsu[k] - zabsu[jp])) * zrdsig
        zqf1 = 1.0 + (1.0 / rdbrv - 1.0) * q[k]
        zqf2 = 1.0 + (1.0 / rdbrv - 1.0) * q[jp]
        zdthds = (t[k] * zqf1 / zskap[k] - t[jp] * zqf2 / zskap[jp]) * zrdsig
        zri = zskaph[k] * gascon * zdthds / (sigmah[k] * zdvds * zdvds)
        if zri <= 0.0:
            zdenomh = 1.0 + 3.0 * cc * b * math.sqrt(-1.0 * zri) * zzkfac * zmixh2
            zdenomm = 1.0 + 3.0 * cc * b * math.sqrt(-1.0 * zri) * zzkfac * zmixm2
            zrifh = 1.0 - 3.0 * b * zri / zdenomh
            zrifm = 1.0 - 2.0 * b * zri / zdenomm
        else:
            zdenom = math.sqrt(1 + d * zri)
            zrifh = 1.0 / (1.0 + 3.0 * b * zri * zdenom)
            zrifm = 1.0 / (1.0 + 2.0 * b * zri / zdenom)
        kh[k] = zfac * zmixh2 * zrifh * zdvds * zrthl * zrthl * zrthl
        km[k] = zfac * zmixm2 * zrifm * zdvds * zrthl * zrthl * zrthl
        to_z[k] = ga * sigmah[k] * zrthl / gascon
    return km, kh, to_z


def nonlocal_column(t, q, u, v, sigma, sigmah, wth0: float, wq0: float, ustar: float, km_local, kh_local, to_z,
                    dtime: float, c: dict, prm=PARAMETERS) -> dict:
    """The nonlocal boundary layer on one column, given the state (top first), the surface's kinematic fluxes of
    potential temperature (K m/s) and vapour (kg/kg m/s), upward, the friction velocity (m/s) and vdiff's sigma-form
    local diffusivities and d sigma / dz at the interfaces. Returns the depth h (m) and the velocity scales; the
    interfaces' diffusivities in sigma form with the profile's inside h (the larger of the two); the explicit
    countergradient and entrainment increments of temperature (K) and vapour (kg/kg) over dtime; and the kinematic
    fluxes they carry (theta in K m/s, vapour in kg/kg m/s, upward, at the interfaces)."""
    c = plasim_constants(c)
    g, gascon, akap, rdbrv = c['ga'], c['gascon'], c['akap'], c['rdbrv']
    t, q, u, v = (np.asarray(a, float) for a in (t, q, u, v))
    sigma, sigmah = np.asarray(sigma, float), np.asarray(sigmah, float)
    km_local, kh_local, to_z = (np.asarray(a, float) for a in (km_local, kh_local, to_z))
    n = t.size
    nb = n - 1                                                  # the lowest level
    dsigma = np.diff(np.r_[0.0, sigmah])
    eps_v = 1.0 / rdbrv - 1.0
    theta = t / sigma ** akap
    thv = theta * (1.0 + eps_v * q)
    zf, zh = heights(t, q, sigma, sigmah, g, gascon, rdbrv)
    thv_n = thv[nb]
    wthv0 = wth0 * (1.0 + eps_v * q[nb]) + eps_v * theta[nb] * wq0
    buoy = g / thv_n * wthv0                                    # the surface buoyancy flux, m2/s3
    unstable = wthv0 > 0.0
    ust3 = ustar ** 3
    kappa, eps = KARMAN, prm['eps_surface']

    def depth(thref, ricr):
        """CAM's pblintd: the bulk number from the lowest level up, h linear in the number between levels."""
        ri_below = 0.0
        for k in range(nb - 1, -1, -1):
            vvk = max((u[k] - u[nb]) ** 2 + (v[k] - v[nb]) ** 2 + prm['b_ustar'] * ustar * ustar, TINY)
            ri = g * (thv[k] - thref) * (zf[k] - zf[nb]) / (thv_n * vvk)
            if ri >= ricr:
                if ri - ri_below <= 0.0:
                    return zf[k + 1]
                return zf[k + 1] + (ricr - ri_below) / (ri - ri_below) * (zf[k] - zf[k + 1])
            ri_below = ri
        return zf[0]

    out = dict(unstable=bool(unstable), wthv0=wthv0, buoyancy=buoy, z_full=zf, z_half=zh[:-1], ustar=ustar)
    if unstable:
        h0 = depth(thv_n, prm['ricr_unstable'])
        wm0 = (ust3 + prm['beta_m'] * eps * kappa * buoy * h0) ** (1.0 / 3.0)
        excess = prm['b_excess'] * wthv0 / wm0
        h = depth(thv_n + excess, prm['ricr_unstable'])
        wstar3 = buoy * h
        wstar = wstar3 ** (1.0 / 3.0)
        phim_inv = (ust3 + prm['beta_m'] * eps * kappa * buoy * h) ** (1.0 / 3.0)   # u* / phi_m(eps h / L), = w_m
        phih_inv = math.sqrt(ust3 + prm['beta_h'] * eps * kappa * buoy * h)         # u*^(3/2) / phi_h(eps h / L)
        wm = phim_inv
        pr = math.sqrt(ustar) * phim_inv / phih_inv + eps * kappa * prm['a_countergradient'] * wstar / wm
        gamma_t = prm['a_countergradient'] * wstar * wth0 / (wm * wm * h)
        gamma_q = prm['a_countergradient'] * wstar * wq0 / (wm * wm * h)
        out.update(h_first=h0, excess=excess, wstar=wstar, wm=wm, pr=pr)
    else:
        h = depth(thv_n, prm['ricr_stable'])
        wstar = wm = 0.0
        pr = 1.0
        gamma_t = gamma_q = 0.0
        out.update(h_first=h, excess=0.0, wstar=0.0, wm=0.0, pr=1.0)
    out['h'] = h
    zi = zh[:-1]                                                # the interfaces' heights
    inside = zi < h
    km_pbl, kh_pbl = np.zeros(n - 1), np.zeros(n - 1)
    for k in np.flatnonzero(inside):
        z = zi[k]
        shape = z * (1.0 - z / h) ** 2
        if unstable:
            km_pbl[k] = kappa * wm * shape
            kh_pbl[k] = km_pbl[k] / pr
        else:
            bneg = -buoy                                        # >= 0 in stable or neutral air
            if ustar <= 0.0:
                wt = 0.0
            elif z * kappa * bneg <= ust3:                      # z/L <= 1: u* / (1 + 5 z/L)
                wt = ustar * ust3 / (ust3 + prm['beta_s'] * z * kappa * bneg)
            else:                                               # z/L > 1: u* / (5 + z/L)
                wt = ustar * ust3 / (prm['beta_s'] * ust3 + z * kappa * bneg)
            km_pbl[k] = kappa * wt * shape
            kh_pbl[k] = km_pbl[k]
    fac = dtime * to_z * to_z / (sigma[1:] - sigma[:-1])        # kinematic diffusivity to sigma form
    km, kh = km_local.copy(), kh_local.copy()
    km[inside] = np.maximum(km_local[inside], km_pbl[inside] * fac[inside])
    kh[inside] = np.maximum(kh_local[inside], kh_pbl[inside] * fac[inside])
    # the countergradient fluxes, K gamma with the diffusivity used, at the interfaces inside h (kinematic)
    flux_cg_t, flux_cg_q = np.zeros(n - 1), np.zeros(n - 1)
    if unstable:
        k_used = kh[inside] / fac[inside]
        flux_cg_t[inside] = k_used * gamma_t
        flux_cg_q[inside] = k_used * gamma_q
    # the entrainment flux: at the interface between the levels bracketing h, falling as (z/h)^3 below
    flux_en_t, flux_en_q = np.zeros(n - 1), np.zeros(n - 1)
    out.update(entrainment_velocity=0.0, inversion_interface=None)
    below = np.flatnonzero(zf <= h)
    jb = int(below[0]) if below.size else nb                    # the highest level not above h
    if unstable and jb > 0:
        ja = jb - 1                                             # the level just above h
        went3 = wstar3 + prm['c_entrain_ustar'] * ust3
        went = went3 ** (1.0 / 3.0)
        flux_thv_h = -prm['a_entrain'] * thv_n / g * went3 / h
        dthv = thv[ja] - thv[jb]
        we = went if dthv <= 0.0 else min(-flux_thv_h / dthv, went)
        ft = -we * (theta[ja] - theta[jb])
        fq = -we * (q[ja] - q[jb])
        flux_en_t[ja], flux_en_q[ja] = ft, fq
        for k in range(ja + 1, n - 1):
            r3 = (zi[k] / h) ** 3
            flux_en_t[k], flux_en_q[k] = ft * r3, fq * r3
        out.update(entrainment_velocity=we, inversion_interface=ja, entrainment_flux_thv=flux_thv_h,
                   inversion_jump_thv=dthv)
    # explicit increments: each interface's flux into the layer below it out of the layer above, in sigma form
    shk = sigmah[:-1] ** akap
    g_t = shk * dtime * to_z * (flux_cg_t + flux_en_t)          # temperature: theta's flux times the interface's Exner
    g_q = dtime * to_z * (flux_cg_q + flux_en_q)
    dt_inc = (np.r_[g_t, 0.0] - np.r_[0.0, g_t]) / dsigma
    dq_inc = (np.r_[g_q, 0.0] - np.r_[0.0, g_q]) / dsigma
    kept = bool(np.any(q + dq_inc < 0.0))
    if kept:                                                    # HB93 appendix B: such a profile keeps its own
        dq_inc = np.zeros(n)
    out.update(inside=inside, km=km, kh=kh, km_profile=km_pbl, kh_profile=kh_pbl, gamma_t=gamma_t, gamma_q=gamma_q,
               flux_countergradient_t=flux_cg_t, flux_countergradient_q=flux_cg_q, flux_entrainment_t=flux_en_t,
               flux_entrainment_q=flux_en_q, dt_explicit=dt_inc, dq_explicit=dq_inc, vapour_kept=kept)
    return out


def _solve(x, k, dsigma):
    """vdiff's elimination for momentum and vapour: dsigma_j (x'_j - x_j) = k_j (x'_j+1 - x'_j) - k_j-1 (x'_j - x'_j-1)."""
    n = x.size
    zebs, xn = np.zeros(n - 1), np.zeros(n)
    zebs[0] = k[0] / (dsigma[0] + k[0])
    xn[0] = dsigma[0] * x[0] / (dsigma[0] + k[0])
    for j in range(1, n - 1):
        denom = dsigma[j] + k[j] + k[j - 1] * (1.0 - zebs[j - 1])
        zebs[j] = k[j] / denom
        xn[j] = (x[j] * dsigma[j] + k[j - 1] * xn[j - 1]) / denom
    xn[n - 1] = (x[n - 1] * dsigma[n - 1] + k[n - 2] * xn[n - 2]) / (dsigma[n - 1] + k[n - 2] * (1.0 - zebs[n - 2]))
    for j in range(n - 2, -1, -1):
        xn[j] = xn[j] + zebs[j] * xn[j + 1]
    return xn


def _solve_t(t, kp, dsigma, zskap):
    """vdiff's elimination for temperature, diffusing potential temperature: kp is the heat diffusivity times
    sigmah^kappa at each interface."""
    n = t.size
    zebs, tn = np.zeros(n - 1), np.zeros(n)
    zebs[0] = kp[0] / (dsigma[0] + kp[0] / zskap[0])
    tn[0] = dsigma[0] * t[0] / (dsigma[0] + kp[0] / zskap[0])
    for j in range(1, n - 1):
        denom = dsigma[j] + (kp[j] + kp[j - 1] * (1.0 - zebs[j - 1] / zskap[j - 1])) / zskap[j]
        zebs[j] = kp[j] / denom
        tn[j] = (t[j] * dsigma[j] + kp[j - 1] / zskap[j - 1] * tn[j - 1]) / denom
    tn[n - 1] = (t[n - 1] * dsigma[n - 1] + kp[n - 2] * tn[n - 2] / zskap[n - 2]) \
        / (dsigma[n - 1] + kp[n - 2] / zskap[n - 1] * (1.0 - zebs[n - 2] / zskap[n - 2]))
    for j in range(n - 2, -1, -1):
        tn[j] = tn[j] + zebs[j] * tn[j + 1] / zskap[j + 1]
    return tn


def vdiff(t, q, u, v, sigma, sigmah, dtime: float, c: dict, lamm: float = LAMBDA_EARTH_M, umin: float = WIND_FLOOR,
          q_heat=None, ndheat: bool = True, surface: dict | None = None, prm=PARAMETERS) -> dict:
    """PlaSim's vertical diffusion on one column (fluxmod.f90, vdiff), from the state after the tendencies so far
    (t, q, u, v, top first), over the implicit step dtime (PlaSim's deltsec2): the tendencies of temperature (the
    diffusion's and the kinetic energy's dissipation, which vdiff adds separately), vapour and wind, and the sigma-form
    diffusivities used. q_heat is the vapour of the moist heat capacity in the dissipation's heating (PlaSim's dq, the
    state before this step's tendencies; q when not given). With surface, dict(wth0=, wq0=, ustar=) (kinematic surface
    fluxes, upward, and the friction velocity), the nonlocal boundary layer acts first (nonlocal_column)."""
    c = plasim_constants(c)
    akap, acpd, adv = c['akap'], c['acpd'], c['ADV']
    t, q, u, v = (np.asarray(a, float).copy() for a in (t, q, u, v))
    sigma, sigmah = np.asarray(sigma, float), np.asarray(sigmah, float)
    dsigma = np.diff(np.r_[0.0, sigmah])
    q_heat = q if q_heat is None else np.asarray(q_heat, float)
    km, kh, to_z = local_coefficients(t, q, u, v, sigma, sigmah, dtime, c, lamm, umin)
    t0, q0 = t.copy(), q.copy()
    pbl = None
    if surface is not None:
        pbl = nonlocal_column(t, q, u, v, sigma, sigmah, surface['wth0'], surface['wq0'], surface['ustar'], km, kh,
                              to_z, dtime, c, prm)
        km, kh = pbl['km'], pbl['kh']
        t = t + pbl['dt_explicit']
        q = q + pbl['dq_explicit']
    un, vn = _solve(u, km, dsigma), _solve(v, km, dsigma)
    dudt, dvdt = (un - u) / dtime, (vn - v) / dtime
    dtdt_ke = np.zeros_like(t)
    if ndheat:
        absu = np.sqrt(u ** 2 + v ** 2)
        ke = 0.5 * absu * absu
        ken = _solve(ke, km, dsigma)
        dtdt_ke = -((un * un - u * u + vn * vn - v * v) / dtime - (ken - ke) / dtime) * 0.5 / acpd / (1.0 + adv * q_heat)
    qn = _solve(q, kh, dsigma)
    dqdt = (qn - q0) / dtime
    zskap = sigma ** akap
    tn = _solve_t(t, kh * sigmah[:-1] ** akap, dsigma, zskap)
    dtdt = (tn - t0) / dtime
    return dict(dtdt=dtdt + dtdt_ke, dtdt_diffusion=dtdt, dtdt_dissipation=dtdt_ke, dqdt=dqdt, dudt=dudt, dvdt=dvdt,
                km=km, kh=kh, to_z=to_z, pbl=pbl)


def dca(t, q, sigma, sigmah, akap: float, dtime: float):
    """PlaSim's dry convective adjustment (rainmod.f90, mkdca) on one column: the lowest pair of layers whose dry
    potential temperature falls with height, extended upward through the layers colder than their mixture, takes the
    mean potential temperature by sigma-thickness times Exner and the mean vapour by sigma-thickness, until no such pair
    remains. As in PlaSim, a column still unstable after 2n + 1 passes is left as it was. Returns the tendencies (K/s,
    kg/kg/s) and the number of passes that adjusted the column (0: it was stable)."""
    t, q = np.asarray(t, float), np.asarray(q, float)
    sigma = np.asarray(sigma, float)
    dsigma = np.diff(np.r_[0.0, np.asarray(sigmah, float)])
    n = t.size
    zskap = sigma ** akap
    zth = t / zskap
    zqn = q.copy()
    jiter = 0
    while True:
        icon = False
        ibase, itop = -1, n - 1
        zsum1 = zsum2 = zsumq1 = zsumq2 = ztht = zqt = 0.0
        for k in range(n - 1):                                  # the last (lowest) unstable pair wins
            kp = k + 1
            if zth[k] < zth[kp]:
                ibase, itop = kp, k - 1
                zsum1 = zskap[kp] * dsigma[kp] + zskap[k] * dsigma[k]
                zsum2 = zskap[kp] * dsigma[kp] * zth[kp] + zskap[k] * dsigma[k] * zth[k]
                zsumq1 = dsigma[kp] + dsigma[k]
                zsumq2 = dsigma[kp] * zqn[kp] + dsigma[k] * zqn[k]
                ztht = zsum2 / zsum1
                zqt = zsumq2 / zsumq1
                icon = True
        if not icon:
            break
        for k in range(n - 3, -1, -1):                          # Fortran's NLEM-1 down to 1
            if itop == k and zth[k] < ztht:
                itop = k - 1
                zsum1 += zskap[k] * dsigma[k]
                zsum2 += zskap[k] * dsigma[k] * zth[k]
                zsumq1 += dsigma[k]
                zsumq2 += dsigma[k] * zqn[k]
                ztht = zsum2 / zsum1
                zqt = zsumq2 / zsumq1
        for k in range(n):
            if itop < k <= ibase:
                zth[k] = ztht
                zqn[k] = zqt
        jiter += 1
        if jiter > 2 * n:
            return np.zeros(n), np.zeros(n), -jiter
    return (zth * zskap - t) / dtime, (zqn - q) / dtime, jiter
