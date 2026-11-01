"""Böhm hydrodynamics for prescribed porous spherical aggregate envelopes.

Vectorized translation of McSnow revision 1e428d7e901cbc5b91f2040e575f736ec5963a74.
Only alpha=1, rough ice, and 0<q<=1 are supported. The caller supplies mass,
diameter, gravity and local air properties; this module does not evolve shape.
See boehm.md for source discrepancies, comparison scope and verification.

Derived from McSnow: Copyright (c) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH,
MeteoSwiss. Authors: Slavko Brdar, Axel Seifert, Christoph Siewert,
Markus Karrer and Jan-Niklas Welss. SPDX-License-Identifier: BSD-3-Clause.
The retained license is references/mcsnow_bsd_3_clause.txt.
"""
from __future__ import annotations

import numpy as np

SOURCE_REVISION = '1e428d7e901cbc5b91f2040e575f736ec5963a74'
ROUGH_ICE_X0 = 2.8e6
# The retained Fortran source leaves these literals unkinded (default REAL).
# Preserve their float32 rounding for source parity, rather than silently changing it.
_SOURCE_SHAPE_QUADRATIC = float(np.float32(9.18))
_SOURCE_SHAPE_CUBIC = float(np.float32(3.53))


def _positive(name, value):
    value = np.asarray(value, dtype=float)
    if not np.isfinite(value).all() or (value <= 0).any():
        raise ValueError(f'{name} must be finite and positive')
    return value


def _shape(phi, q):
    """McSnow hydrodynamicValuesBoehm, also used for the reduced two-body state."""
    k = np.minimum(np.minimum(np.maximum(.82 + .18*phi, .85), .37 + .63/phi),
                   1.33/(np.maximum(np.log(phi), 0.) + 1.19))
    gamma = np.maximum(1., np.minimum(1.98, 3.76-8.41*phi+_SOURCE_SHAPE_QUADRATIC*phi**2-_SOURCE_SHAPE_CUBIC*phi**3))
    cdp = np.maximum(.292*k*gamma, .492-.2/np.sqrt(phi))
    cdp *= np.maximum(1., q*(1.46*q-.46))
    return k, gamma, cdp


def particle_properties(D, m, g, state, q, alpha=1, sigma=0):
    """Return broadcast numeric arrays for rough, alpha=1 porous aggregates.

    D is outer diameter [m], m mass [kg], g gravity [m/s2], q the effective
    projected-area/circumscribing-circle ratio. state['rho'] is kg/m3 and
    state['eta'] is dynamic viscosity [Pa s]. sigma is the source's internal
    lognormal speed-width parameter, not ambient turbulence. All returned
    fields are arrays (scalars are zero-dimensional arrays).
    """
    D, m, g, rho, eta, q, alpha, sigma = np.broadcast_arrays(
        _positive('D', D), _positive('m', m), _positive('g', g),
        _positive('rho', state['rho']), _positive('eta', state['eta']),
        _positive('q', q), np.asarray(alpha, dtype=float), np.asarray(sigma, dtype=float))
    if (q > 1).any():
        raise ValueError('this adapter requires 0 < q <= 1')
    if not (alpha == 1).all():
        raise ValueError('this adapter supports fall-axis aspect alpha=1 only')
    if not np.isfinite(sigma).all() or (sigma < 0).any():
        raise ValueError('sigma must be finite and nonnegative')
    k, gamma, cdp = _shape(alpha, q)
    X = 8*m*g*rho/(np.pi*eta**2*q**.25)
    # This algebra avoids overflow in (X/X0)^2 while retaining source limits.
    ratio = X/ROUGH_ICE_X0
    with np.errstate(over='ignore'):
        transition = 1. - 1./(1. + ratio**2)
    cdp_prime = cdp*(1. + .6*transition)
    increment = cdp_prime/(6*k)*np.sqrt(X/cdp_prime)
    beta = increment/(np.sqrt(1. + increment) + 1.)
    re0 = 6*k/cdp_prime*beta**2
    gamma_small = (4.5*k**2-cdp)/(4*cdp)
    Re = re0*(1 + 2*beta*np.exp(-beta*gamma_small)/((2+beta)*(1+beta)))
    v = eta*Re/(rho*D)
    r = D/2
    r_e = r*np.sqrt(q)
    delta = 3.60*r/np.sqrt(Re*gamma)
    CDI = X/Re**2-24*k/Re
    result = dict(D=D, m=m, g=g, q=q, alpha=alpha, sigma=sigma, rho=rho, eta=eta,
                  r=r, r_e=r_e, alpha_e=r/r_e, q_e=q, j=np.full(D.shape, 2.),
                  k=k, delta=delta, CDI=CDI, v=v, Re=Re, X=X)
    if not all(np.isfinite(x).all() for x in result.values()):
        raise FloatingPointError('nonfinite Böhm particle state; inputs exceed numerical range')
    return {key: np.asarray(value) for key, value in result.items()}


def _velocity_sign(v1, v2, sigma):
    log_ratio = np.log(v1/v2)
    deterministic = sigma < 1e-10
    safe_sigma = np.where(deterministic, 1., sigma)
    arg = (log_ratio + safe_sigma**2)/(2*safe_sigma)
    # Fortran SIGN(1, LOG(1)) is +1, unlike numpy.sign(0).
    sign = np.where(log_ratio >= 0, 1., -1.)
    return np.where(deterministic, sign,
                    np.tanh(2/np.sqrt(np.pi)*arg + .1012*arg**3))


def collision_efficiency(left, right, state, *, diagnostics=False):
    """Broadcast source hydrodynamic efficiency; no sticking or outer area.

    For a pair matrix, pass each left property with [:, None] and each right
    property with [None, :]. Equal-speed pairs return the source efficiency;
    the caller's differential-settling kernel still vanishes through |dv|.
    diagnostics=True returns Ec and intermediate arrays, including anomaly
    flags. Raw finite efficiencies outside [0,1] are retained, never clipped.
    """
    rho, eta = _positive('rho', state['rho']), _positive('eta', state['eta'])
    for properties in (left, right):
        if not np.array_equal(np.broadcast_arrays(properties['rho'], rho)[0],
                              np.broadcast_arrays(properties['rho'], rho)[1]):
            raise ValueError('particle density differs from collision air state')
        if not np.array_equal(np.broadcast_arrays(properties['eta'], eta)[0],
                              np.broadcast_arrays(properties['eta'], eta)[1]):
            raise ValueError('particle viscosity differs from collision air state')
    if not np.equal(left['g'], right['g']).all():
        raise ValueError('particles must use the same gravity')
    faster_left = np.asarray(left['v']) > np.asarray(right['v'])
    keys = ('v', 'm', 'r', 'r_e', 'alpha', 'alpha_e', 'q_e', 'j', 'k', 'delta', 'CDI', 'sigma')
    one = {k: np.where(faster_left, left[k], right[k]) for k in keys}
    two = {k: np.where(faster_left, right[k], left[k]) for k in keys}
    v1, v2, j = one['v'], two['v'], one['j']
    vi1 = v1*(2+j)/4*np.minimum(one['alpha'], 1.)
    mstar = one['m']*two['m']/(one['m']+two['m'])
    rstar = one['r']*two['r']/(one['r_e']+two['r_e'])
    astar = np.pi*rstar**2
    phi1 = _velocity_sign(v1, v2, one['sigma'])
    phi2 = _velocity_sign(v2, v1, two['sigma'])
    vstar = np.maximum(np.abs(v1*phi1+v2*phi2), 1e-10)
    delta_star = one['delta']*np.sqrt(v1/vstar)+two['delta']*np.sqrt(v2/vstar)
    alpha_star = rstar*(one['alpha_e']/one['r_e']+two['alpha_e']/two['r_e'])
    qstar = rstar*(one['q_e']/one['r_e']+two['q_e']/two['r_e'])
    kstar, _, cdpstar = _shape(alpha_star, qstar)
    Restar = 2*rstar*vstar*rho/eta
    cd0star = 4.5*kstar**2*np.maximum(alpha_star, 1.)
    beta = np.sqrt(Restar*cdpstar/(6*kstar))
    gamma = (cd0star-cdpstar)/(4*cdpstar)
    cdistar = cdpstar*(1-4/beta*np.expm1(-gamma*beta))
    G = 6*kstar*eta/(rho*rstar*vstar*cdistar)
    F = np.sqrt(G**2+one['CDI']*vi1**2/(cdistar*vstar**2))
    H = 2*mstar/(rho*astar*cdistar*delta_star)
    fh = F/H
    small = fh < 10
    # Evaluate each branch on its own mask: the unused hyperbolic branch can
    # overflow for large particles. The source's >=10 approximation is kept.
    safe_fh = np.where(small, fh, 0.)
    low = H*np.log(np.cosh(safe_fh)+(1+G)/F*np.sinh(safe_fh))-G
    high = H*np.log((F+G)/(2*F))+F-G
    base = np.where(small, low, high)
    boundary = base**(2/j)
    b = 3*v1/(one['r_e']+one['delta'])
    k2x = np.where(two['alpha_e'] < 1, .57+.43*two['alpha_e'], two['k']**1.15)
    cx = two['m']/(3*np.pi*two['r_e']*eta*k2x)
    cy = two['m']/(3*np.pi*two['r_e']*eta*two['k'])
    bcy = b*cy
    active = bcy > .5
    delx = np.sqrt(2*b*cx/j+1)
    dely = np.sqrt(np.maximum(2*bcy-1, 0.))
    # atan2 reproduces both source angular branches and the bcy=1 limit.
    tdelta = cy/np.where(active, dely, 1.)*np.arctan2(dely, bcy-1)
    ratio_t = tdelta/cx
    x = delx*ratio_t
    # Scaled exponential form of cosh(x)+sinh(x)/delx avoids overflow.
    log_term = ((delx-1)*ratio_t
                + np.log((1+1/delx)/2 + (1-1/delx)/2*np.exp(-2*x)))
    potential_high = np.exp(-j*log_term)
    potential_low = (2*delx/(1+delx)*np.exp(-(delx-1)*ratio_t))**j
    potential = np.where(active, np.where(bcy >= 1, potential_high, potential_low), 0.)
    Ec = boundary*potential
    if not np.isfinite(Ec).all():
        raise FloatingPointError('nonfinite Böhm collision efficiency; inspect source applicability')
    if diagnostics:
        return dict(Ec=np.asarray(Ec), boundary_efficiency=boundary,
                    potential_efficiency=potential, vstar=vstar, F=F, G=G, H=H,
                    F_over_H=fh, b_cy=bcy, Re_star=Restar, alpha_star=alpha_star,
                    q_star=qstar, outside_unit_interval=(Ec < 0) | (Ec > 1))
    return np.asarray(Ec)
