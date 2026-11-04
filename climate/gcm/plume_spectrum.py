"""A spectral plume model of a convecting column's free troposphere, for the Open Moon at lunar gravity (4 November
2026): a spread of entraining plumes from the lowest level's air, each rising to its own neutral level, and the
environment's equilibrium temperature that the spread sets, with PlaSim's own thermodynamics (sbm.condense: the ice's
latent heat below 273.16 K, the condensate dropped) so that its result can serve the simplified Betts-Miller scheme
(sbm, sbm_redesign2) as a reference. Arrays run top first, as PlaSim's levels do.

The formulation follows the buoyancy sorting of Zhou and Xie (2019, J. Atmos. Sci. 76, 2801-2814): plumes with
distinct entrainment rates, each positively buoyant up to its own level of neutral buoyancy, with a one-to-one relation
between a plume's entrainment and that level. At each level the environment stands at the temperature of the plume
whose neutral level it is; the less entraining plumes, warmer, rise on through it and the more entraining ones have
stopped below. The environment's lapse rate is therefore smaller than any single plume's, and its deficit below the
undiluted parcel can stop growing, or shrink, with height (the C shape). With one plume the environment is that plume
everywhere, which is Singh and O'Gorman's (2013) zero-buoyancy plume: sbm_redesign2.zbp_reference.

The closure (the relation between entrainment and neutral level) is given here as the distribution of the plumes' tops:
G(s), the share of the base mass flux whose plumes have stopped by s, s the environment's height above the
condensation level of the lowest level's air (metres, or e-foldings of pressure), from counted cloud tops
(tops_closure) or from a convective mass-flux profile (tops_from_mass_flux). The plumes are ordered from the most
to the least entraining, so the plume neutral at s is the one at quantile G(s) of the base mass flux counted from the
most entraining end (quantile matching). Each plume keeps its own mass flux from its base to its top, its lateral
detrainment equal to its entrainment, and leaves the column at its top; the spread's mass flux is then 1 - G(s) of the
base's, its fractional detrainment the tops' outflow (dG/dz) / (1 - G) plus the lateral exchange, and Romps's (2014)
relative humidity delta / (delta + gamma) follows, gamma = -d ln q*/dz of the environment. The closure is this module's construction from the published
summaries of Zhou and Xie's model; their own way of fixing the relation was not read (the paper's text was not
reachable). Arakawa and Schubert's (1974) quasi-equilibrium closure for the mass flux is not used.

Each plume rises dry adiabatically with its vapour kept and undiluted to the condensation level of the lowest level's
air, crossing a layer wholly below it in one step (sbm.parcel's way); every other layer is crossed in nsub equal steps
of log pressure, each lifting all plumes dry adiabatically to the step's end and condensing them, then finding the
environment there (the closure, between the two plumes that bracket its quantile, linearly in their index) and mixing
into each plume a share 1 - exp(-e dz) of environment air at the environment's temperature and (rh + x (1 - rh)) of
its saturation: dz the step's metres above the base from the plume's own virtual temperature, and e the plume's rate
per metre ('constant') or a / z with z its height above the ground at the step's middle ('inverse_z', zbp_reference's
form). x is Romps's (2014) share of freshly detrained, saturated cloud in the entrained air, which turns the dilution
eps into (1 - x) eps. The temperature recorded at a level is the plumes' after condensing there and before mixing, and
the environment's is the closure's choice among them.

Each plume's moist static energy h = acpd T + g z_h + ALV q (z_h its height from its own dry temperature, so that the
dry lift keeps acpd T + g z_h to third order in the step) changes only by entrainment, mix (h_env - h), and by the
scheme's own condensation convention, acpd lcp rain - ALV rain (the vapour's heat capacity in cp and the ice's latent
heat below freezing); both are summed, with the dry lift's third-order departure, and the remainder, rounding, is
returned per plume as the budget's residual.
"""
from __future__ import annotations

import math

import numpy as np

from climate.gcm import sbm

FAMILIES = ('constant', 'inverse_z')
NSUB = 16                          # steps a layer above the condensation level; zbp_reference's (within 0.2 K of 64)
N_PLUMES = 32                      # nodes of a spectrum built by log_spectrum


def saturation(t, p, c=sbm.PLASIM):
    """sbm.saturation for an array of temperatures (K) at one pressure or an array of pressures (Pa)."""
    t = np.asarray(t, float)
    q = c['rdbrv'] * c['ra1'] * np.exp(c['ra2'] * (t - c['TMELT']) / (t - c['ra4'])) / p
    q = np.minimum(q, c['rdbrv'])
    return q / (1.0 - (1.0 / c['rdbrv'] - 1.0) * q)


def condense(t, q, p, c=sbm.PLASIM):
    """sbm.condense for an array of plumes at pressure p (Pa): their temperature and vapour after the two condensing
    steps, and the change of acpd T + ALV q (J/kg) that the scheme's own convention makes in each, the latent heat
    taken over cp with the vapour's heat capacity and the ice's below freezing."""
    t, q = np.array(t, float), np.array(q, float)
    defect = np.zeros_like(t)
    for _ in range(2):
        qs = c['rdbrv'] * c['ra1'] * np.exp(c['ra2'] * (t - c['TMELT']) / (t - c['ra4'])) / p
        qs = np.minimum(qs, c['rdbrv'])
        cor = 1.0 / (1.0 - (1.0 / c['rdbrv'] - 1.0) * qs)
        qs = qs * cor
        wet = qs < q
        if not np.any(wet):
            break
        lcp = np.where(t < c['TMELT'], c['ALS'], c['ALV']) / (c['acpd'] * (1.0 + c['ADV'] * q))
        dqsdt = c['ra2'] * (c['TMELT'] - c['ra4']) * qs * cor / ((t - c['ra4']) ** 2)
        rain = np.where(wet, (q - qs) / (1.0 + lcp * dqsdt), 0.0)
        defect += (c['acpd'] * lcp - c['ALV']) * rain
        t, q = t + lcp * rain, q - rain
    return t, q, defect


def log_spectrum(lo, hi, n=N_PLUMES):
    """n plumes' entrainment parameters spaced evenly in log between lo and hi, most entraining first, with equal base
    mass flux (a spread uniform in the logarithm of the rate)."""
    params = np.exp(np.linspace(math.log(hi), math.log(lo), n)) if n > 1 else np.array([float(hi)])
    return params, np.full(n, 1.0 / n)


def tops_closure(heights, shares):
    """The tops' distribution G(s) from points: s above the condensation level (metres or e-foldings) and the share of
    the base mass flux stopped by s, nondecreasing from 0 to 1, linear between the points, 0 below the first and 1 above
    the last."""
    s, g = np.asarray(heights, float), np.asarray(shares, float)
    if np.any(np.diff(s) <= 0.0) or np.any(np.diff(g) < 0.0) or g[0] < 0.0 or g[-1] > 1.0:
        raise ValueError('the tops need rising heights and nondecreasing shares within 0..1')

    def closure(x):
        return float(np.interp(x, s, g, left=0.0, right=1.0))
    closure.points = (s.tolist(), g.tolist())
    return closure


def tops_from_mass_flux(heights, mass_flux, base):
    """The tops' distribution from a convective mass-flux profile (heights above the ground, m, rising; kg/m2/s), as a
    cloud-resolving model's updraughts give it: no plume stops below the height where the mass flux peaks, and above
    it the share stopped is 1 - M(z) / M_peak, kept nondecreasing and reaching 1 where the mass flux ends; heights
    above base (m, the condensation level's height)."""
    z, m = np.asarray(heights, float), np.asarray(mass_flux, float)
    k = int(np.argmax(m))
    g = np.clip(np.maximum.accumulate(1.0 - m[k:] / m[k]), 0.0, 1.0)
    end = next((i for i in range(g.size) if g[i] >= 0.995), g.size - 1)
    s = np.r_[0.0, z[k:k + end + 1] - base]
    shares = np.r_[0.0, g[:end + 1]]
    shares[-1] = 1.0
    keep = np.r_[True, np.diff(s) > 0.0]
    return tops_closure(s[keep], shares[keep])


def equilibrium(t_src, q_src, sigma, ps, c, akap, params, family='constant', weights=None, tops=None,
                tops_units='m', rh=0.70, x=0.0, nsub=NSUB, base=None) -> dict:
    """The environment's equilibrium temperature by level set by a spread of plumes from the lowest level's air
    (temperature t_src, vapour q_src at sigma[-1]), as the module's docstring states it. params: the plumes'
    entrainment, per metre ('constant') or the coefficient a of a / z ('inverse_z'); weights: their base mass-flux
    shares (equal by default); tops: the closure G(s) (tops_closure), s in tops_units 'm' (the environment's metres
    above the condensation level) or 'efold' (its e-foldings of pressure above it), needed for more than one plume; rh:
    the entrained air's relative humidity, a number or one value a level (top first, interpolated in log pressure);
    x: Romps's share of detrained cloud in the entrained air; base: the condensation level (Pa) of the lowest level's
    air, found as sbm.condensation_level finds it when None.

    Returns, by level: t_env, the environment's temperature (K); z_env, its height above the ground (m); quantile, the
    share of the base mass flux stopped below the level (G); param_env, the entrainment parameter of the plume neutral
    there; mass_flux, the spread's mass flux over the base's (1 - G); eps_bulk, the mass-flux-weighted entrainment of
    the plumes still rising (per metre); detrainment, the spread's fractional detrainment (per metre); rh_romps,
    Romps's (1 - x) delta / ((1 - x) delta + gamma). By plume (most entraining first): params, weights, t_plumes
    (plume by level), top_s (where its quantile is reached, in tops_units), top_z (m above the ground), top_p (Pa) and
    the moist static energy budget to the top level (J/kg): budget_entrainment, budget_convention, budget_lift (the dry
    lift's third-order departure, largest in the one-step layers below the base) and budget_residual, the rest. Also
    base (Pa), z_base (m), sorting (the largest error, K, of the plumes' order against the environment: a plume still
    rising colder than it, or a stopped one warmer) and steps. param_env and eps_bulk for 'inverse_z' take a and a / z
    with z the environment's height."""
    sigma = np.asarray(sigma, float)
    n = sigma.size
    p = ps * sigma
    params = np.atleast_1d(np.asarray(params, float))
    if family not in FAMILIES:
        raise ValueError(f'family is one of {FAMILIES}, not {family!r}')
    order = np.argsort(-params, kind='stable')
    params = params[order]
    weights = np.full(params.size, 1.0 / params.size) if weights is None else np.asarray(weights, float)[order]
    weights = weights / weights.sum()
    m = params.size
    if m > 1 and tops is None:
        raise ValueError('a spread of plumes needs the tops distribution (tops_closure)')
    if tops_units not in ('m', 'efold'):
        raise ValueError(f"tops_units is 'm' or 'efold', not {tops_units!r}")
    mid = np.cumsum(weights) - 0.5 * weights                   # each plume's quantile, most entraining first
    index = np.arange(m, dtype=float)
    if base is None:
        base = sbm.condensation_level(t_src, q_src, p[-1], p[0], akap, c)
    lb = math.log(base / ps)
    scale = c['gascon'] / c['ga']
    epsv = 1.0 / c['rdbrv'] - 1.0
    rh_levels = np.full(n, float(rh)) if np.ndim(rh) == 0 else np.asarray(rh, float)
    lnsig = np.log(sigma)
    inverse_z = family == 'inverse_z'

    def pick(temps, quantile):
        """The environment: the temperature of the plume at that quantile, between its two neighbours."""
        if m == 1:
            return temps[0], 0.0
        pos = float(np.interp(quantile, mid, index))
        i0 = min(int(pos), m - 2)
        f = pos - i0
        return temps[i0] * (1.0 - f) + temps[i0 + 1] * f, pos

    t_pl = np.zeros((m, n))
    t_env, z_env = np.zeros(n), np.zeros(n)
    quant, pos_lv = np.zeros(n), np.zeros(n)
    tn, qn, _ = condense(np.full(m, float(t_src)), np.full(m, float(q_src)), p[-1], c)
    t_pl[:, -1] = tn
    t_env[-1] = tn[0]
    z = scale * tn * (1.0 + epsv * qn) * math.log(1.0 / sigma[-1])     # each plume's height above the ground
    zh = scale * tn * math.log(1.0 / sigma[-1])                         # its height from its dry temperature (h)
    h_src = c['acpd'] * tn + c['ga'] * zh + c['ALV'] * qn
    ent, conv, lift = np.zeros(m), np.zeros(m), np.zeros(m)
    te, qe = float(tn[0]), float(qn[0])
    ze = float(z[0])
    z_env[-1] = ze
    z_base = None                                       # the environment's height at the condensation level
    steps = 0
    for k in range(n - 2, -1, -1):
        x1, x0 = lnsig[k + 1], lnsig[k]
        if min(x1, lb) - x0 <= 0.0:                     # wholly below the base: one dry step, as sbm.parcel
            t_old = tn
            tn = tn * (sigma[k] / sigma[k + 1]) ** akap
            t_lift = tn
            tn, qn, d = condense(tn, qn, p[k], c)
            conv += d
            z = z + scale * 0.5 * (t_old + tn) * (1.0 + epsv * qn) * (x1 - x0)
            zh = zh + scale * 0.5 * (t_old + t_lift) * (x1 - x0)
            lift += c['acpd'] * (t_lift - t_old) + c['ga'] * scale * 0.5 * (t_old + t_lift) * (x1 - x0)
            te_new, qe_new = float(tn[0]), float(qn[0])
            ze += scale * 0.5 * (te + te_new) * (1.0 + epsv * 0.5 * (qe + qe_new)) * (x1 - x0)
            te, qe = te_new, qe_new
            t_pl[:, k], t_env[k], z_env[k] = tn, te, ze
            steps += 1
            continue
        width = (x1 - x0) / nsub
        xa = x1
        for s in range(1, nsub + 1):
            xb = x0 if s == nsub else x1 - s * width
            t_old = tn
            tn = tn * math.exp(akap * (xb - xa))
            t_lift = tn
            pb = ps * math.exp(xb)
            tn, qn, d = condense(tn, qn, pb, c)
            conv += d
            dz_step = scale * 0.5 * (t_old + tn) * (1.0 + epsv * qn) * (xa - xb)
            zh = zh + scale * 0.5 * (t_old + t_lift) * (xa - xb)
            lift += c['acpd'] * (t_lift - t_old) + c['ga'] * scale * 0.5 * (t_old + t_lift) * (xa - xb)
            ze_new = ze + scale * te * (1.0 + epsv * qe) * (xa - xb)          # the environment's height, predicted
            if z_base is None and xb < lb:
                z_base = ze + (ze_new - ze) * (xa - lb) / (xa - xb)
            above = 0.0 if z_base is None else (ze_new - z_base if tops_units == 'm' else lb - xb)
            g = 0.0 if m == 1 else tops(max(above, 0.0))
            te_new, pos = pick(tn, g)
            rh_b = float(np.interp(xb, lnsig[::-1], rh_levels[::-1]))
            qe_new = (rh_b + x * (1.0 - rh_b)) * float(saturation(te_new, pb, c))
            mixed = min(xa, lb) - xb
            if mixed > 0.0:
                dz = dz_step * mixed / (xa - xb)
                e = params / (z + 0.5 * dz_step) if inverse_z else params
                mix = 1.0 - np.exp(-e * dz)
                hp = c['acpd'] * tn + c['ga'] * zh + c['ALV'] * qn
                he = c['acpd'] * te_new + c['ga'] * zh + c['ALV'] * qe_new
                ent += mix * (he - hp)
                t_rec = tn
                tn = (1.0 - mix) * tn + mix * te_new
                qn = (1.0 - mix) * qn + mix * qe_new
            else:
                t_rec = tn
            z = z + dz_step
            ze += scale * 0.5 * (te + te_new) * (1.0 + epsv * 0.5 * (qe + qe_new)) * (xa - xb)
            te, qe = te_new, qe_new
            xa = xb
            steps += 1
        t_pl[:, k], t_env[k], z_env[k] = t_rec, te, ze
        quant[k], pos_lv[k] = g, pos
    h_end = c['acpd'] * tn + c['ga'] * zh + c['ALV'] * qn
    residual = h_end - h_src - ent - conv - lift
    # Diagnostics by level: the plume neutral there, the spread's mass flux, its entrainment and detrainment.
    z_b = z_env[-1] if z_base is None else z_base
    param_env = np.interp(pos_lv, index, params) if m > 1 else np.full(n, params[0])
    rising = mid[None, :] > quant[:, None]                      # (level, plume): still rising above the level
    mass = 1.0 - quant
    w = rising * weights[None, :]
    rate = params[None, :] / np.maximum(z_env[:, None], 1.0) if inverse_z else np.broadcast_to(params, (n, m))
    eps_bulk = np.where(w.sum(1) > 0.0, (w * rate).sum(1) / np.maximum(w.sum(1), 1e-300), 0.0)
    above_base = z_env > z_b
    dg = np.gradient(quant, z_env)
    detr = np.where(above_base & (mass > 1e-9), eps_bulk + dg / np.maximum(mass, 1e-9), np.nan)
    qs_env = saturation(t_env, p, c)
    gamma = -np.gradient(np.log(qs_env), z_env)
    dd = (1.0 - x) * detr
    rh_romps = np.where(above_base, dd / (dd + gamma), np.nan)
    sorting = 0.0
    for k in range(n - 1):
        if not above_base[k] or m == 1:
            continue
        up = rising[k]
        if np.any(up):
            sorting = max(sorting, float(np.max(t_env[k] - t_pl[up, k])))
        if np.any(~up):
            sorting = max(sorting, float(np.max(t_pl[~up, k] - t_env[k])))
    if m > 1:
        s_lv = (z_env - z_b) if tops_units == 'm' else (lb - lnsig)
        top_s = np.array([_inverse(tops, u) for u in mid])
        top_z = np.interp(top_s, s_lv[::-1], z_env[::-1])
        top_p = np.exp(np.interp(top_z, z_env[::-1], np.log(p)[::-1]))
    else:
        top_s = top_z = top_p = np.array([np.nan])
    return dict(t_env=t_env, z_env=z_env, quantile=quant, param_env=param_env, mass_flux=mass, eps_bulk=eps_bulk,
                detrainment=detr, rh_romps=rh_romps, params=params, weights=weights, t_plumes=t_pl, top_s=top_s,
                top_z=top_z, top_p=top_p, budget_residual=residual, budget_entrainment=ent, budget_convention=conv,
                budget_lift=lift, base=base, z_base=z_b, sorting=sorting, steps=steps)


def _inverse(tops, u, lo=0.0, hi=None):
    """The s at which the tops' distribution reaches u: from its points when it has them (tops_closure), otherwise by
    bisection between lo and hi (200 km by default) for any nondecreasing callable."""
    pts = getattr(tops, 'points', None)
    if pts is not None:
        s, g = np.asarray(pts[0]), np.asarray(pts[1])
        if u <= g[0]:
            return float(s[0])
        if u >= g[-1]:
            return float(s[-1])
        return float(np.interp(u, g, s))
    hi = 2.0e5 if hi is None else hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if tops(mid) < u:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def deficit(t_ref, t, q, sigma, ps, c, akap, nsub=NSUB):
    """The finely integrated undiluted parcel from the lowest level (sbm.diluted_parcel at dilution 0 in nsub steps a
    layer, the dilution study's yardstick) minus t_ref, by level (K)."""
    fine = sbm.diluted_parcel(np.asarray(t, float), np.asarray(q, float), np.asarray(sigma, float), ps, 0.0, akap, c,
                              max_mix=sbm.MAX_MIX, min_steps=nsub)
    return fine - np.asarray(t_ref, float)
