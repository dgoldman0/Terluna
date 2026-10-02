"""Non-inductive charge separation when graupel meets ice crystals or snow and they bounce apart.

A graupel particle of diameter D_g falling at V_g(D_g) sweeps a smaller particle of diameter D_x falling at V_x(D_x);
the pair collides at rate (pi/4)(D_g + D_x)^2 |V_g - V_x| per unit of both number densities, bounces apart with the
separation efficiency, and each bounce leaves charge dq on the graupel and -dq on the other particle. For two
exponential size distributions N0 exp(-lam D) the charging rate per unit volume is

    rho_dot = N_g N_x sum_ij w_i w_j K(u_i / lam_g, u_j / lam_x),   K = (pi/4)(D_g + D_x)^2 |V_g - V_x| E_sep dq,

with N = N0 / lam the number per unit volume and (u, w) Gauss-Laguerre nodes and weights. The charge per bounce comes
from a laboratory law (`laws`), each a function of the crystal diameter, the impact speed, the temperature and the
rime accretion rate.
"""
from __future__ import annotations
import numpy as np

NODES, WEIGHTS = np.polynomial.laguerre.laggauss(24)


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
