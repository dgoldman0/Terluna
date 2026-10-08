"""The state the upper air settles to under the light the solar shield lets through, from the limb tracing.

atmosphere/middle_atmosphere/limb_heat.py traces sunlight below 175 nm through the shield's window stack, its annulus
films, gaps over its aperture and open sky beyond it, along slant paths through the Open Moon's air. It tabulates the
heat each part leaves in the thermosphere against the heat that swells the air, for several protected radii, and
writes them to results/limb_heat.json. A state is the first heat, counting up from zero, at which the heat the swollen
air takes, with the sky's Lyman-alpha glow, equals the heat that swells it; where none lies in the tables the air runs
away past their top. The author adopted this count of the UV transmission through gaps for O1 on 2026-10-07
(research/decisions.md), with FISM2's measured rise of the X-rays at solar maximum.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

import numpy as np

from atmosphere.middle_atmosphere import escape

ROOT = Path(__file__).resolve().parents[2]
LIMB = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results' / 'limb_heat.json'
LIMB_SCHEMA = 'terluna.atmosphere.limb-heat/3'
FILM = 'annulus_4_um'           # the ring fleet's annulus film (decisions.md, 2026-10-07)
HEATING_SHARE = 0.1             # a protected radius covers the heating when light outside it gives at most this share
SHAPE_X = np.round(np.arange(0.0, 40.0001, 0.05), 4)   # log pressures above the base for a state's heating shape
_PRODUCT = {}


def limb_product():
    """The limb tracing's product, checked against the schema this module reads."""
    if 'product' not in _PRODUCT:
        product = json.loads(LIMB.read_text())
        if product['schema'] != LIMB_SCHEMA:
            raise ValueError(f"unexpected limb-heat schema {product['schema']}")
        _PRODUCT['product'] = product
    return _PRODUCT['product']


def entries(case):
    """The tabulated heats of one case (shield and treatment joined, as in limb_heat.CASES)."""
    return limb_product()['cases'][case]['by_profile_heat']


def protected_radii():
    return [float(x) for x in limb_product()['design']['summary_radii_R']]


def activity_of(activity):
    """An activity as a dict: the loss response's quiet Sun or solar maximum by name, or a dict as given."""
    if isinstance(activity, str):
        from atmosphere.loss_response.model import ACTIVITY    # read late: the loss response imports this module
        return ACTIVITY[activity]
    return activity


def scales(activity, bands=None):
    """Each band's factor over the quiet Sun for an activity: its own 'bands' (a measured spectrum) if it has them;
    otherwise FISM2's rise of the X-rays by band below 10 nm (escape.xray_scale) and the activity's single factor on
    the ultraviolet above."""
    activity = activity_of(activity)
    bands = limb_product()['design']['bands_nm'] if bands is None else bands
    if 'bands' in activity:
        out = np.asarray(activity['bands'], float)
        if len(out) != len(bands) - 1:
            raise ValueError('a measured spectrum needs one factor per band')
        return out
    return np.array([float(escape.xray_scale(.5 * (lo + hi), activity['xray'])) if hi <= 10.0 else activity['uv']
                     for lo, hi in zip(bands[:-1], bands[1:])])


def band_factors(wavelength_nm, activity, bands=None):
    """The factor over the quiet Sun at each wavelength: its band's (see scales)."""
    bands = limb_product()['design']['bands_nm'] if bands is None else bands
    index = np.clip(np.searchsorted(bands, np.asarray(wavelength_nm, float), side='right') - 1, 0, len(bands) - 2)
    return scales(activity, bands)[index]


def parts(rows, film, activity, glow, protected_R, bands=None):
    """Heat (W/m^2 of lunar surface) each part of the light leaves in the thermosphere against the heat that swells
    the air, at one protected radius and solar activity (a name or a dict, see scales): the window stack over the
    window, the film over the annulus, sunlight beyond the aperture, the sky's glow (W/m^2, already scaled for
    activity), and the whole band through gaps over the aperture per unit of transmission. The bands default to the
    limb product's."""
    rows = [e for e in rows if 'radii' in e and not e['outflow_limit']]
    q = np.array([e['heat_W_m2'] for e in rows])
    key = f'{float(protected_R):g}'
    s = scales(activity, bands)

    def get(source, zone):
        return np.array([float(np.dot(e['radii'][key][source][zone], s)) for e in rows])
    return q, dict(window=get('window', 'window'), annulus=get(film, 'annulus'),
                   beyond_aperture=get('open', 'outside'), glow=np.full(len(q), glow), gap_unit=get('open', 'aperture'))


def _share_below(quantiles, levels, x):
    """A part's cumulative share of heat below each log pressure x, from its quantiles."""
    q = np.asarray(quantiles, float)
    xp = np.r_[0.0, q, q[-1] + max(q[-1] - q[-2], 1e-3)]
    xp = np.maximum.accumulate(xp) + np.arange(len(xp)) * 1e-9          # strictly rising
    return np.interp(x, xp, np.r_[0.0, levels, 1.0])


def shape(rows, film, activity, glow, glow_shape, protected_R, gaps, heat, design=None):
    """The cumulative share of a state's heat below each log pressure above the base (on SHAPE_X): each part of the
    traced light where the limb tables put it, at the state's heat between the tabulated heats and mixed over its groups
    of bands by their heat at the activity, and the sky's glow (glow, W/m^2) where O2 absorbs it (glow_shape, on
    SHAPE_X)."""
    design = limb_product()['design'] if design is None else design
    bands, groups = design['bands_nm'], design['shape_groups_nm']
    levels = np.asarray(design['shape_quantiles'], float)
    key = f'{float(protected_R):g}'
    if float(protected_R) not in [float(x) for x in design['shape_radii_R']]:
        raise ValueError(f'the limb tables follow the heat in height for {design["shape_radii_R"]} lunar radii only')
    rows = [e for e in rows if 'radii' in e and not e['outflow_limit']]
    q = np.array([e['heat_W_m2'] for e in rows])
    s = scales(activity, bands)
    mids = .5 * (np.asarray(bands[:-1]) + np.asarray(bands[1:]))
    group = np.searchsorted(groups, mids, side='right') - 1
    n = len(groups) - 1
    j = int(np.clip(np.searchsorted(q, heat), 1, len(q) - 1))
    w = float(np.clip((heat - q[j - 1]) / (q[j] - q[j - 1]), 0.0, 1.0))
    pieces = (('window', 'window', 'window', 1.0), (film, 'annulus', f'annulus_{key}', 1.0),
              ('open', 'outside', f'outside_{key}', 1.0), ('open', 'aperture', f'aperture_{key}', gaps))
    below, total = np.zeros_like(SHAPE_X), 0.0
    for i, weight in ((j - 1, 1.0 - w), (j, w)):
        for source, zone, place, factor in pieces:
            heats = np.bincount(group, weights=np.asarray(rows[i]['radii'][key][source][zone]) * s, minlength=n) * factor
            for g, quantiles in enumerate(rows[i]['shapes'][place][source]):
                if heats[g] > 0 and quantiles is not None:
                    below += weight * heats[g] * _share_below(quantiles, levels, SHAPE_X)
                    total += weight * heats[g]
    below += glow * np.asarray(glow_shape)
    total += glow
    if total <= 0:
        raise ValueError('a state with no heat has no heating shape')
    out = np.clip(below / total, 0.0, 1.0)
    out[0], out[-1] = 0.0, 1.0
    return np.maximum.accumulate(out)


def first_state(q, p, gaps):
    """The first heat, counting up from zero, at which the heat the swollen air takes equals the heat that swells it,
    with a grey transmission `gaps` through gaps over the aperture; None if the heat runs past the tables' top."""
    fine = np.linspace(0., q[-1], 4001)
    taken = sum(np.interp(fine, q, p[k]) for k in ('window', 'annulus', 'beyond_aperture', 'glow'))
    excess = taken + gaps * np.interp(fine, q, p['gap_unit']) - fine
    settled = np.nonzero(excess <= 0)[0]
    if not len(settled):
        return None
    i = settled[0]
    return float(fine[i - 1] + excess[i - 1] * (fine[i] - fine[i - 1]) / (excess[i - 1] - excess[i]))


def shares(q, p, gaps, heat):
    """Each part's heat at a state (W/m^2), and the share of the heat that comes from beyond the aperture."""
    out = {k: float(np.interp(heat, q, v)) for k, v in p.items() if k != 'gap_unit'}
    out['gaps'] = float(gaps * np.interp(heat, q, p['gap_unit']))
    out['beyond_share'] = out['beyond_aperture'] / heat if heat > 0 else 0.
    return out


def allowed(q, p, heat):
    """Largest transmission through gaps whose state stays at or below the heat; None if the films, beyond-aperture
    light and glow alone pass it."""
    state = first_state(q, p, 0.)
    if state is None or state > heat:
        return None
    if (first_state(q, p, .3) or math.inf) <= heat:
        return .3
    lo, hi = math.log(1e-9), math.log(.3)
    for _ in range(60):
        mid = .5 * (lo + hi)
        s = first_state(q, p, math.exp(mid))
        lo, hi = (mid, hi) if s is not None and s <= heat else (lo, mid)
    return math.exp(.5 * (lo + hi))


def heating_radius(rows, film, activity, glow, gaps, share=HEATING_SHARE):
    """The smallest protected radius (lunar radii) at which light from beyond the aperture gives at most `share` of
    the heat of the state it settles to, interpolated between the tabulated radii; None if no tabulated radius does
    or the air runs away at every one."""
    radii = protected_radii()
    values = []
    for x in radii:
        q, p = parts(rows, film, activity, glow, x)
        state = first_state(q, p, gaps)
        values.append(math.inf if state is None else shares(q, p, gaps, state)['beyond_share'])
    if values[0] <= share:
        return radii[0]
    for i in range(1, len(radii)):
        if values[i] <= share:
            v0, v1 = values[i - 1], values[i]
            if not math.isfinite(v0):
                return radii[i]
            return radii[i - 1] + (v0 - share) * (radii[i] - radii[i - 1]) / (v0 - v1)
    return None

