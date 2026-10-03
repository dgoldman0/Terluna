"""Non-inductive charge separation when graupel meets ice crystals or snow and they bounce apart.

A graupel particle of diameter D_g falling at V_g(D_g) sweeps a smaller particle of diameter D_x falling at V_x(D_x);
the pair collides at rate (pi/4)(D_g + D_x)^2 |V_g - V_x| per unit of both number densities, bounces apart with the
separation efficiency, and each bounce leaves charge dq on the graupel and -dq on the other particle. For two
exponential size distributions N0 exp(-lam D) the charging rate per unit volume is

    rho_dot = N_g N_x sum_ij w_i w_j K(u_i / lam_g, u_j / lam_x),   K = (pi/4)(D_g + D_x)^2 |V_g - V_x| E_sep dq,

with N = N0 / lam the number per unit volume and (u, w) Gauss-Laguerre nodes and weights. The charge per bounce comes
from a laboratory law, each a function of the crystal diameter, the impact speed, the temperature and the cloud water
or the rime accretion rate. The laws follow their implementation in WRF-ELEC (MicroTed/wrf4-elec at commit e43041b,
elec/module_mp_nssl_2mom_elec.F, in the public domain under the WRF notice; Mansell et al. 2005, J. Geophys. Res.
110, D12101; Fierro et al. 2013, Mon. Wea. Rev. 141, 2390):

- 'saunders_peck': Saunders and Peck (1998) as WRF-ELEC's default (nssl_isaund = 12, routine saund6 with idelq = 1
  and the Brooks et al. 1997 critical rime accretion rate above -15 C). dq = B D^a |dV|^b q(RAR, T) fC, with B, a, b
  by crystal size and sign (Mansell et al. 2005 Table 1, from Brooks et al. 1997) and RAR = E_cw LWC V_graupel in
  g m-2 s-1. No charge below RAR = 0.1 or below -32.47 C. The laboratory data behind it were taken at 3-14 m/s.
- 'takahashi': Takahashi (1978) through WRF-ELEC's lookup table (Wojcik 1994; elec/takahashi.txt, 0 to -30 C, 0.01 to
  30 g/m3 of cloud water, linear in both) scaled by alpha = min(5 (D/100 um)^2 |dV|/(8 m/s), 10) (routine TAKA). The
  laboratory data were taken at 9 m/s.
- 'low_speed_kumar' and 'low_speed_avila': the charge per bounce measured near the lunar impact speeds, held fixed:
  Pradeep Kumar et al. (2024, J. Earth Syst. Sci. 133, 226) Eq. 3 for 2-mm graupel at 1.2 m/s, +(1.392 C^2 - 3.149 C
  + 8.108) fC for cloud water C of 0.1-1.72 g/m3 from -7 to -18.5 C; and the top of Avila et al.'s (2013, J. Geophys.
  Res. Atmos. 118, 6680) +0.01-0.2 fC at 1-3 m/s from -7 to -13 C. They bracket what slow impacts may carry.

The per-bounce charge is capped at 500 fC and the charging rate at 3 nC m-3 s-1, as in WRF-ELEC.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np

NODES, WEIGHTS = np.polynomial.laguerre.laggauss(24)
FEMTO = 1.0e-15
DQ_CAP_C = 500.0 * FEMTO                    # WRF-ELEC's cap on the charge per bounce (delqxxa, delqnxa)
RATE_CAP_C_M3_S = 3.0e-9                    # and on the charging rate (scxacymax)
RIMING_EFFICIENCY = 0.7                     # graupel's collection of cloud droplets, Morrison's ECI in the CM1 runs


def pair_rate(n_g, lam_g, n_x, lam_x, speed_g, speed_x, transfer, e_sep=1.0):
    """Charging rate (C m-3 s-1, positive when graupel gains positive charge) of graupel against another ice type,
    for cells given by arrays of the same shape: number densities n (m-3), slopes lam (1/m), fall-speed functions
    speed(d) -> m/s taking the diameter array (broadcast against the cells), the charge per bounce
    transfer(d_g, d_x, dv) -> C, and the separation efficiency (a number or a function of d_x)."""
    n_g, lam_g, n_x, lam_x = (np.asarray(a, float)[..., None, None] for a in (n_g, lam_g, n_x, lam_x))
    d_g = NODES[:, None] / lam_g
    d_x = NODES[None, :] / lam_x
    dv = np.abs(speed_g(d_g) - speed_x(d_x))
    eff = e_sep(d_x) if callable(e_sep) else e_sep
    kernel = 0.25 * np.pi * (d_g + d_x) ** 2 * dv * eff * transfer(d_g, d_x, dv)
    weight = WEIGHTS[:, None] * WEIGHTS[None, :]
    return (n_g * n_x * (kernel * weight).sum(axis=(-2, -1), keepdims=True))[..., 0, 0]


def collision_rate(n_g, lam_g, n_x, lam_x, speed_g, speed_x):
    """Graupel collisions with the other type per unit volume per second (m-3 s-1)."""
    return pair_rate(n_g, lam_g, n_x, lam_x, speed_g, speed_x, lambda d_g, d_x, dv: np.ones_like(dv))


def power_law(dq_ref_c: float, d_ref_m: float, v_ref_m_s: float, a: float, b: float):
    """A charge-per-bounce law dq = dq_ref (D_x / d_ref)^a (dv / v_ref)^b, for sensitivity to the speed exponent b."""
    def transfer(d_g, d_x, dv):
        return dq_ref_c * (d_x / d_ref_m) ** a * (dv / v_ref_m_s) ** b
    return transfer


def rar_critical(t_c):
    """The critical rime accretion rate (g m-2 s-1) above which graupel charges positively, as WRF-ELEC's default:
    Brooks et al. (1997) above -15 C, Saunders and Peck's (1998) polynomial from -15 to -23.7 C, its cubic
    continuation below (to zero at -37 C) or the polynomial where that is larger, 0.1 below -33 C, and never under 0.1."""
    t = np.asarray(t_c, float)
    poly = 1.0 + t * (7.9262e-2 + t * (4.4847e-2 + t * (7.4754e-3 + t * (5.4686e-4 + t * (1.6737e-5 + t * 1.7613e-7)))))
    cubic = 3.39608 * (1.0 - np.abs(((t + 23.7) / (-23.7 + 37.0)) ** 3))
    rarc = np.where(t > -15.0, np.clip(-1.47 - 0.2 * t, 0.0, 3.29), np.where(t > -23.7, poly, cubic))
    rarc = np.maximum(rarc, 0.1)
    rarc = np.where(t <= -23.7, np.maximum(0.0, poly), rarc)            # idelq = 1: the polynomial below -23.7 C
    return np.where(t < -33.0, 0.1, rarc)


def saunders_peck_charge_fc(rar, t_c):
    """The charge function q(RAR, T) (fC) of WRF-ELEC's default Saunders-Peck scheme, before the size and speed
    factor: positive above the critical rate, a parabola through zero at 0.1 and at the critical rate below it, none
    at or below 0.1 or at -32.47 C and colder, falling linearly to that from -25 C."""
    rar, t = np.asarray(rar, float), np.asarray(t_c, float)
    rarc = rar_critical(t)
    fac = np.where(t > -25.0, 1.0, np.clip(1.0 + (t + 25.0) / (32.47 - 25.0), 0.0, 1.0))
    positive = 6.74 * (rar - rarc)
    span = np.where(np.abs(rarc - 0.1) > 1e-9, rarc - 0.1, 1e-9)
    negative = np.minimum(0.0, np.minimum(1.0, 0.6 * np.abs(rarc - 0.1)) * 6.5
                          * (-1.0 + 4.0 / span ** 2 * (rar - (rarc + 0.1) / 2.0) ** 2))
    q = np.where(rar > rarc, positive, negative)
    return np.where((rar <= 0.1) | (t <= -32.47), 0.0, fac * q)


def size_speed_factor(d_x, dv, positive):
    """B D^a |dV|^b of Mansell et al. (2005) Table 1 (D in m, dV in m/s): by crystal size, separately for charge of
    each sign on the graupel."""
    um = np.asarray(d_x, float) * 1.0e6
    dv = np.abs(np.asarray(dv, float))
    pos = np.where(um < 155.0, 4.9e13 * d_x ** 3.76, np.where(um <= 452.0, 4.0e6 * d_x ** 1.9, 52.8 * d_x ** 0.44))
    neg = np.where(um < 253.0, 5.24e8 * d_x ** 2.54, 24.0 * d_x ** 0.5)
    return np.where(positive, pos * dv ** 2.5, neg * dv ** 2.8)


def saunders_peck(t_c, lwc_g_m3, speed_g, other=None):
    """The Saunders-Peck charge per bounce, transfer(d_g, d_x, dv) in C, for cells at temperature t_c with cloud
    water lwc_g_m3; each graupel size rimes at RAR = E_cw LWC V_g(d_g)."""
    def transfer(d_g, d_x, dv):
        rar = RIMING_EFFICIENCY * lwc_g_m3 * speed_g(d_g)
        q = saunders_peck_charge_fc(rar, t_c)
        return np.clip(size_speed_factor(d_x, dv, q > 0.0) * q * FEMTO, -DQ_CAP_C, DQ_CAP_C)
    return transfer


def takahashi_table(path):
    """Takahashi's charge per collision (fC) from WRF-ELEC's table: the cloud water values (g/m3, rising) and the
    charges at 0, -1, ..., -30 C for each."""
    values = np.array(Path(path).read_text().split(), float).reshape(-1, 32)
    values = values[values[:, 0] > 0.0]                     # the first record is the header of temperatures, 0 to -30
    order = np.argsort(values[:, 0])
    return values[order, 0], values[order, 1:]


def takahashi_charge_fc(lwc_g_m3, t_c, table):
    """Takahashi's charge (fC) at cloud water lwc_g_m3 and temperature t_c, interpolated linearly in both as the TAKA
    routine does, the -30 C values used colder and rolled off to none at -40 C; none below 0.01 g/m3."""
    cwc, charge = table
    lwc = np.minimum(np.asarray(lwc_g_m3, float), 30.0)
    t = np.asarray(t_c, float)
    tt = np.clip(-t, 0.0, 30.0)
    lo = np.clip(np.floor(tt).astype(int), 0, 29)
    frac = tt - lo
    by_t = (1.0 - frac[..., None]) * charge.T[lo] + frac[..., None] * charge.T[lo + 1]   # (..., len(cwc))
    k = np.clip(np.searchsorted(cwc, lwc, side='right') - 1, 0, len(cwc) - 2)
    w = np.clip((lwc - cwc[k]) / (cwc[k + 1] - cwc[k]), 0.0, 1.0)
    q = (1.0 - w) * np.take_along_axis(by_t, k[..., None], -1)[..., 0] + w * np.take_along_axis(by_t, k[..., None] + 1, -1)[..., 0]
    fac = np.where(t > -30.0, 1.0, np.clip(1.0 - ((t + 30.0) / 10.0) ** 2, 0.0, 1.0))
    return np.where(lwc >= 0.01, fac * q, 0.0)


def takahashi(table):
    """The Takahashi charge per bounce with WRF-ELEC's size and speed factor, as a law."""
    def law(t_c, lwc_g_m3, speed_g=None, other=None):
        q = takahashi_charge_fc(lwc_g_m3, t_c, table)

        def transfer(d_g, d_x, dv):
            alpha = np.minimum(5.0 * (d_x / 100.0e-6) ** 2 * np.abs(dv) / 8.0, 10.0)
            return np.clip(alpha * q * FEMTO, -DQ_CAP_C, DQ_CAP_C)
        return transfer
    return law


def low_speed_kumar(t_c, lwc_g_m3, speed_g=None, other=None):
    """Pradeep Kumar et al. (2024) Eq. 3 at 1.2 m/s: +(1.392 C^2 - 3.149 C + 8.108) fC per bounce for cloud water C of
    0.1-1.72 g/m3 (held at the ends) from -7 to -18.5 C, none outside that range."""
    c = np.clip(np.asarray(lwc_g_m3, float), 0.1, 1.72)
    inside = (np.asarray(t_c, float) <= -7.0) & (np.asarray(t_c, float) >= -18.5)
    q = np.where(inside, 1.392 * c ** 2 - 3.149 * c + 8.108, 0.0) * FEMTO

    def transfer(d_g, d_x, dv):
        return q * np.ones_like(dv)
    return transfer


def low_speed_avila(t_c, lwc_g_m3, speed_g=None, other=None):
    """The top of Avila et al.'s (2013) charges at 1-3 m/s: +0.2 fC per bounce from -7 to -13 C, none outside."""
    t = np.asarray(t_c, float)
    q = np.where((t <= -7.0) & (t >= -13.0), 0.2, 0.0) * FEMTO

    def transfer(d_g, d_x, dv):
        return q * np.ones_like(dv)
    return transfer


def laws(input_dir=None) -> dict:
    """The laws by name; Takahashi's needs its table among the recorded inputs."""
    from atmosphere.electricity import fetch_inputs
    out = dict(saunders_peck=saunders_peck)
    try:
        out['takahashi'] = takahashi(takahashi_table(fetch_inputs.path('takahashi.txt')))
    except FileNotFoundError:
        pass
    out.update(low_speed_kumar=low_speed_kumar, low_speed_avila=low_speed_avila)
    return out
