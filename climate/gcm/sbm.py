"""The simplified Betts-Miller convection scheme as the Terluna patch puts it into PlaSim (exoplasim_run.SBM_CODE), one
column at a time in Python, for its tests.

The scheme follows Frierson (2007, J. Atmos. Sci. 64, 1959-1976). A parcel lifted from the lowest level, dry and then
moist adiabatically with PlaSim's saturation (the ice's latent heat below freezing, as Kuo's parcel), and mixing in a
share 1 - exp(-eps ln(p_below/p)) of the surrounding air in each layer it rises through, gives the reference
temperature from the lowest level to the top of its buoyancy; the reference vapour is rh of saturation there. Where the
parcel is buoyant somewhere above its start, with its inhibition below cin when that is set and, with the CAPE trigger
set, the CAPE of a second parcel mixing eps_trigger per e-folding of pressure at least cape_limit, the temperature and
vapour relax toward the references over tau. The rain is the vapour removed; the heating is shifted by one amount
through the convecting layers so that it releases the rain's latent heat, keeping the column's moist enthalpy. A
buoyant column whose vapour lies below its reference, or whose inhibition or trigger holds it back, is left to PlaSim's
shallow scheme. Arrays run top first, as PlaSim's levels do.

Two options of the cloud-depth patch (exoplasim_run.SBM_DEPTH_PATCHES), off by default:
- the depth test (deep_pa > 0): the column convects deep only where the parcel's cloud, from the condensation level of
  the lowest level's air to the scheme's top level, is more than deep_pa thick; a thinner one goes to the shallow
  scheme, as the IFS separates deep from shallow convection at 200 hPa (documentation Cy36r1, Part IV, section 6.4);
- the mixing within layers (steps > 1, with eps > 0): the parcel crosses each layer in that many equal steps of log
  pressure, each lifting it dry adiabatically, mixing in its step's share of the surrounding air interpolated there in
  log pressure, and condensing the mixture, the last step ending at the level itself. One step is the scheme as before.
"""
from __future__ import annotations

import math

import numpy as np

# PlaSim's constants as the tests take them (pumamod and the run's planet); the Fortran harness is given the same.
PLASIM = dict(ALV=2.5008e6, ALS=2.8345e6, acpd=1004.6, ADV=0.8608, TMELT=273.16, ra1=610.78, ra2=17.2693882,
              ra4=35.86, rdbrv=0.622, gascon=287.05, ga=1.62)
LCL_ITERATIONS = 24                # halvings of the condensation level's bracket in log pressure, as the Fortran's


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


def parcel(t, q, sigma, p, eps, akap, c=PLASIM, steps=1, ps=None):
    """The scheme's parcel lifted from the lowest level, dry and then moist adiabatically, mixing in a share
    1 - exp(-eps ln(p_below/p)) of the surrounding air in each layer it rises through: its temperature by level. With
    steps > 1 and eps > 0 the mixing is integrated within each layer in that many equal steps of log pressure, the
    surrounding air interpolated linearly in log pressure; ps, the surface pressure, gives the steps' pressures."""
    tp = np.zeros(t.size)
    tn, qn = condense(t[-1], q[-1], p[-1], c)
    tp[-1] = tn
    for k in range(t.size - 2, -1, -1):
        if eps > 0.0 and steps > 1:
            sub = math.log(sigma[k + 1] / sigma[k]) / float(steps)
            mix = 1.0 - math.exp(-eps * sub)
            shrink, lift = math.exp(-sub), math.exp(-akap * sub)       # each step's sigma ratio and dry lift
            sig0 = sigma[k + 1]
            for s in range(1, steps + 1):
                if s < steps:
                    sig1 = sig0 * shrink
                    tn *= lift
                else:                                                  # the last step ends at the level itself
                    sig1 = sigma[k]
                    tn *= (sig1 / sig0) ** akap
                w = float(s) / float(steps)
                tn = (1.0 - mix) * tn + mix * ((1.0 - w) * t[k + 1] + w * t[k])
                qn = (1.0 - mix) * qn + mix * ((1.0 - w) * q[k + 1] + w * q[k])
                tn, qn = condense(tn, qn, ps * sig1, c)
                sig0 = sig1
        else:
            tn *= (sigma[k] / sigma[k + 1]) ** akap
            if eps > 0.0:
                mix = 1.0 - math.exp(-eps * math.log(sigma[k + 1] / sigma[k]))
                tn, qn = (1.0 - mix) * tn + mix * t[k], (1.0 - mix) * qn + mix * q[k]
            tn, qn = condense(tn, qn, p[k], c)
        tp[k] = tn
    return tp


def condensation_level(t0, q0, p0, p_top, akap, c=PLASIM, iterations=LCL_ITERATIONS):
    """The pressure (Pa) at which air of temperature t0 and vapour q0 at p0, lifted dry adiabatically with its vapour
    kept, reaches PlaSim's saturation: p0 if it is saturated there, p_top if it is still unsaturated at p_top, and
    otherwise found by bisection in log pressure between the two."""
    if q0 >= saturation(t0, p0, c):
        return p0
    if q0 < saturation(t0 * (p_top / p0) ** akap, p_top, c):
        return p_top
    lo, hi = math.log(p_top), math.log(p0)
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        pm = math.exp(mid)
        if q0 >= saturation(t0 * (pm / p0) ** akap, pm, c):
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (lo + hi))


def trigger_cape(t, tp, sigma, dsigma, c=PLASIM):
    """A parcel's CAPE (J/kg) as the scheme counts its inhibition: the gas constant times (tp - t) dsigma/sigma, summed
    from its level of free convection, the lowest above its start where it is buoyant, up through the layers it stays
    buoyant; nil for a parcel buoyant nowhere above its start."""
    lfc = next((k for k in range(t.size - 2, -1, -1) if tp[k] > t[k]), None)
    cape = 0.0
    if lfc is not None:
        for k in range(lfc, -1, -1):
            if tp[k] <= t[k]:
                break
            cape += c['gascon'] * (tp[k] - t[k]) * dsigma[k] / sigma[k]
    return cape


def column(t, q, sigma, dsigma, ps, tau=7200.0, rh=0.7, eps=0.0, cin_limit=0.0, step=1800.0, akap=0.2857,
           c=PLASIM, cape_limit=0.0, eps_trigger=0.0, deep_pa=0.0, steps=1):
    """One column's tendencies by the scheme: temperature (K/s) and vapour (kg/kg/s) by level, the rain formed at each
    level (kg/m2/s), the convection's top level (index, or None), the parcel's inhibition (J/kg), whether the
    column is left to the shallow scheme and the parcel's temperature, the reference; with the CAPE trigger set
    (cape_limit > 0), the trigger parcel's CAPE too, and with the depth test set (deep_pa > 0, Pa), the cloud's base
    and depth (Pa). steps > 1 integrates the parcel's mixing within each layer. step is PlaSim's deltsec2, the shortest
    relaxation time."""
    t, q = np.asarray(t, float), np.asarray(q, float)
    sigma, dsigma = np.asarray(sigma, float), np.asarray(dsigma, float)
    n = t.size
    p, m = ps * sigma, ps * dsigma / c['ga']
    cp = c['acpd'] * (1.0 + c['ADV'] * q)
    out = dict(dtdt=np.zeros(n), dqdt=np.zeros(n), rain=np.zeros(n), top=None, cin=0.0, shallow=False,
               reference_t=None)
    tp = parcel(t, q, sigma, p, eps, akap, c, steps=steps, ps=ps)
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
    if deep_pa > 0.0:                                       # the depth test: a thin cloud goes to the shallow scheme
        base = condensation_level(t[-1], q[-1], p[-1], p[top], akap, c)
        out.update(base=base, depth=base - p[top])
        if base - p[top] <= deep_pa:
            out['shallow'] = True
            return out
    if cape_limit > 0.0:
        out['trigger_cape'] = trigger_cape(t, parcel(t, q, sigma, p, eps_trigger, akap, c), sigma, dsigma, c)
    tau = max(tau, step)
    layers = range(top, n)
    dtr = {k: (tp[k] - t[k]) / tau for k in layers}
    dqr = {k: (rh * saturation(tp[k], p[k], c) - q[k]) / tau for k in layers}
    pq = -sum(dqr[k] * m[k] for k in layers)
    pt = sum(cp[k] * dtr[k] * m[k] for k in layers)
    cpm = sum(cp[k] * m[k] for k in layers)
    if pq <= 0.0 or (cin_limit > 0.0 and cin > cin_limit) or (cape_limit > 0.0 and out['trigger_cape'] < cape_limit):
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
