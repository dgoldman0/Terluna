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

ROOT = Path(__file__).resolve().parents[2]
LIMB = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results' / 'limb_heat.json'
LIMB_SCHEMA = 'terluna.atmosphere.limb-heat/2'
FILM = 'annulus_4_um'           # the ring fleet's annulus film (decisions.md, 2026-10-07)
HEATING_SHARE = 0.1             # a protected radius covers the heating when light outside it gives at most this share
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


def parts(rows, film, activity, glow, protected_R):
    """Heat (W/m^2 of lunar surface) each part of the light leaves in the thermosphere against the heat that swells
    the air, at one protected radius: the window stack over the window, the film over the annulus, sunlight beyond the
    aperture, the sky's glow (W/m^2, already scaled for activity), and the whole band through gaps over the aperture
    per unit of transmission."""
    rows = [e for e in rows if 'radii' in e and not e['outflow_limit']]
    q = np.array([e['heat_W_m2'] for e in rows])
    key = f'{float(protected_R):g}'

    def get(source, zone):
        return np.array([e['radii'][key][source][zone][activity] for e in rows])
    return q, dict(window=get('window', 'window'), annulus=get(film, 'annulus'),
                   beyond_aperture=get('open', 'outside'), glow=np.full(len(q), glow), gap_unit=get('open', 'aperture'))


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

