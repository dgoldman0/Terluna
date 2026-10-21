"""The simplified Betts-Miller convection scheme as the Terluna patch puts it into PlaSim (exoplasim_run.SBM_CODE), one
column at a time in Python, for its tests.

The scheme follows Frierson (2007, J. Atmos. Sci. 64, 1959-1976). A parcel lifted from the lowest level, dry and then
moist adiabatically with PlaSim's saturation (the ice's latent heat below freezing, as Kuo's parcel), and mixing in a
share 1 - exp(-eps ln(p_below/p)) of the surrounding air in each layer it rises through, gives the reference
temperature from the lowest level to the top of its buoyancy; the reference vapour is rh of saturation there. Where the
parcel is buoyant somewhere above its start, with its inhibition below cin when that is set, the temperature and vapour
relax toward the references over tau. The rain is the vapour removed; the heating is shifted by one amount through the
convecting layers so that it releases the rain's latent heat, keeping the column's moist enthalpy. A buoyant column
whose vapour lies below its reference, or whose inhibition holds it back, is left to PlaSim's shallow scheme. Arrays
run top first, as PlaSim's levels do.
"""
from __future__ import annotations

import math

import numpy as np

# PlaSim's constants as the tests take them (pumamod and the run's planet); the Fortran harness is given the same.
PLASIM = dict(ALV=2.5008e6, ALS=2.8345e6, acpd=1004.6, ADV=0.8608, TMELT=273.16, ra1=610.78, ra2=17.2693882,
              ra4=35.86, rdbrv=0.622, gascon=287.05, ga=1.62)


def saturation(t, p, c=PLASIM):
    """PlaSim's saturation specific humidity (kg/kg) at temperature t (K) and pressure p (Pa)."""
    q = c['rdbrv'] * c['ra1'] * math.exp(c['ra2'] * (t - c['TMELT']) / (t - c['ra4'])) / p
    q = min(q, c['rdbrv'])
    return q / (1.0 - (1.0 / c['rdbrv'] - 1.0) * q)


def condense(t, q, p, c=PLASIM):
    """A supersaturated parcel condensed to saturation in two steps, as Kuo's parcel is."""
    for _ in range(2):
        qs = c['rdbrv'] * c['ra1'] * math.exp(c['ra2'] * (t - c['TMELT']) / (t - c['ra4'])) / p
        qs = min(qs, c['rdbrv'])
        cor = 1.0 / (1.0 - (1.0 / c['rdbrv'] - 1.0) * qs)
        qs *= cor
        if qs >= q:
            break
        lcp = (c['ALS'] if t < c['TMELT'] else c['ALV']) / (c['acpd'] * (1.0 + c['ADV'] * q))
        dqsdt = c['ra2'] * (c['TMELT'] - c['ra4']) * qs * cor / ((t - c['ra4']) ** 2)
        rain = (q - qs) / (1.0 + lcp * dqsdt)
        t, q = t + lcp * rain, q - rain
    return t, q


def column(t, q, sigma, dsigma, ps, tau=7200.0, rh=0.7, eps=0.0, cin_limit=0.0, step=1800.0, akap=0.2857,
           c=PLASIM):
    """One column's tendencies by the scheme: temperature (K/s) and vapour (kg/kg/s) by level, the rain formed at each
    level (kg/m2/s), the convection's top level (index, or None), the parcel's inhibition (J/kg), whether the
    column is left to the shallow scheme and the parcel's temperature, the reference. step is PlaSim's deltsec2,
    the shortest relaxation time."""
    t, q = np.asarray(t, float), np.asarray(q, float)
    sigma, dsigma = np.asarray(sigma, float), np.asarray(dsigma, float)
    n = t.size
    p, m = ps * sigma, ps * dsigma / c['ga']
    cp = c['acpd'] * (1.0 + c['ADV'] * q)
    out = dict(dtdt=np.zeros(n), dqdt=np.zeros(n), rain=np.zeros(n), top=None, cin=0.0, shallow=False,
               reference_t=None)
    tp = np.zeros(n)
    tn, qn = condense(t[-1], q[-1], p[-1], c)
    tp[-1] = tn
    for k in range(n - 2, -1, -1):
        tn *= (sigma[k] / sigma[k + 1]) ** akap
        if eps > 0.0:
            mix = 1.0 - math.exp(-eps * math.log(sigma[k + 1] / sigma[k]))
            tn, qn = (1.0 - mix) * tn + mix * t[k], (1.0 - mix) * qn + mix * q[k]
        tn, qn = condense(tn, qn, p[k], c)
        tp[k] = tn
    out['reference_t'] = tp.copy()
    buoyant = [k for k in range(n - 2, -1, -1) if tp[k] > t[k]]
    if not buoyant:
        return out
    lfc = buoyant[0]
    top = lfc
    for k in range(lfc - 1, -1, -1):
        if tp[k] <= t[k]:
            break
        top = k
    cin = sum(c['gascon'] * max(0.0, t[k] - tp[k]) * dsigma[k] / sigma[k] for k in range(lfc + 1, n))
    out.update(top=top, cin=cin)
    tau = max(tau, step)
    layers = range(top, n)
    dtr = {k: (tp[k] - t[k]) / tau for k in layers}
    dqr = {k: (rh * saturation(tp[k], p[k], c) - q[k]) / tau for k in layers}
    pq = -sum(dqr[k] * m[k] for k in layers)
    pt = sum(cp[k] * dtr[k] * m[k] for k in layers)
    cpm = sum(cp[k] * m[k] for k in layers)
    if pq <= 0.0 or (cin_limit > 0.0 and cin > cin_limit):
        out['shallow'] = True
        return out
    corr = (c['ALV'] * pq - pt) / cpm
    cond = sum(max(0.0, -dqr[k]) * m[k] for k in layers)
    for k in layers:
        out['dtdt'][k] = dtr[k] + corr
        out['dqdt'][k] = dqr[k]
        if dqr[k] < 0.0:
            out['rain'][k] = -dqr[k] * pq / cond * m[k]
    return out
