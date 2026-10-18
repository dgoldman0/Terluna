"""What one saved lake column shows: how far the nights mix it, how much of its floor they reach, how often it turns
over, how much oxygen its floor and deep water keep, and what an overturn does to its surface.

Each column is a case of ladder.py, saved by columns.py: the step-by-step mixed depth, the Sun's height and the
surface, bottom and deepest-tenth oxygen of every oxygen column, and the whole column about every day. The floor of
each layer is the area the basin loses across it (lake.py's Basin.sediment), so a layer's share of the floor is its
share of the lake's bed in plan. A night runs from sunset to the next sunset, so its mixing takes in the cooling that
goes on after dawn.
"""
from __future__ import annotations

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_oxygen import run as R

HYPOXIC, ANOXIC = R.HYPOXIC, R.ANOXIC
AFTER_DAYS = R.AFTER_OVERTURN_DAYS
QUANTILES = (0.1, 0.5, 0.9)


def q(x, qs=QUANTILES, nd=3):
    x = np.asarray(x, dtype=float)
    return [round(float(np.quantile(x, v)), nd) for v in qs] if x.size else None


def basin_of(rec):
    return L.Basin(np.asarray(rec['z'], float), np.asarray(rec['area'], float))


def night_reach(rec):
    """Each night's deepest mixing (m), with the morning after it, as the water goes on cooling past sunrise: from
    each sunset to the next at the lake, the step ranges with it. A lunar cycle holds one night, so this is
    lake_oxygen's deepest mixing of each cycle grouped by its night."""
    mu = np.asarray(rec['mu'], float)
    md = np.asarray(rec['mixed_depth'], float)
    up = mu > 0.0
    sunset = np.flatnonzero(up[:-1] & ~up[1:]) + 1
    ns = [(int(a), int(b)) for a, b in zip(sunset[:-1], sunset[1:])]
    return np.array([md[a:b].max() for a, b in ns]), ns


def overturn_cycles(rec, bsn):
    """The lunar cycles at the lake, sunrise to sunrise, and whether the mixed layer reached the bottom layer in each
    (lake_oxygen's complete overturn)."""
    mu = np.asarray(rec['mu'], float)
    md = np.asarray(rec['mixed_depth'], float)
    cyc = R.cycles(mu)
    full = np.array([md[a:b].max() >= bsn.z[-2] - 1e-6 for a, b in cyc], dtype=bool)
    return cyc, full


def runs_of(flags, value):
    """Lengths of the runs of consecutive entries equal to value."""
    out, n = [], 0
    for f in flags:
        if f == value:
            n += 1
        elif n:
            out.append(n)
            n = 0
    if n:
        out.append(n)
    return out


def episodes(rec, bsn):
    """Overturn episodes: runs of consecutive lunar cycles with a complete overturn, as (first step, last step)."""
    cyc, full = overturn_cycles(rec, bsn)
    out, start = [], None
    for (a, b), f in zip(cyc, full):
        if f and start is None:
            start = a
        if not f and start is not None:
            out.append((start, a))
            start = None
    if start is not None:
        out.append((start, cyc[-1][1]))
    return out


def floor_by_depth(bsn, depths):
    """The share of the floor lying within each depth."""
    sed = bsn.sediment / bsn.sediment.sum()
    return [round(float(sed[bsn.centre <= d].sum()), 4) for d in depths]


def column(rec, summary, m, names):
    """The statistics of oxygen column m (index into the run's 24) and of the column's mixing. An overturn is a
    lunar cycle in which the mixed layer reaches the bottom layer; consecutive ones make an episode."""
    bsn = basin_of(rec)
    dt = float(rec['dt'])
    years = len(rec['mu']) * dt / (365.25 * 86400.0)
    sed = bsn.sediment / bsn.sediment.sum()
    centre = bsn.centre
    reach, ns = night_reach(rec)
    to_floor = reach >= bsn.z[-2] - 1e-6                             # the nights that mix the lake to its floor
    reached = reach[:, None] >= centre[None, :]                      # nights x layers
    layer_share = reached.mean(0)
    floor_each_night = reached.astype(float) @ sed
    ox = np.asarray(rec['oxygen'][:, :, m], float)                   # days x layers
    oxic, anoxic = (ox >= HYPOXIC).mean(0), (ox < ANOXIC).mean(0)
    eps = episodes(rec, bsn)
    cyc, full = overturn_cycles(rec, bsn)
    cycle_years = float(np.mean([b - a for a, b in cyc])) * dt / (365.25 * 86400.0)
    quiet = runs_of(full, False)
    so = np.asarray(rec['surface_o2'][:, m], float)
    md = np.asarray(rec['mixed_depth'], float)
    events = []                                  # each lunar cycle with a complete overturn, from its overturn on
    for (a, b), f in zip(cyc, full):
        if not f:
            continue
        a = a + int(np.argmax(md[a:b] >= bsn.z[-2] - 1e-6))
        seg = so[a:b]
        below = seg < HYPOXIC
        back = np.flatnonzero(~below[int(np.argmin(seg)):]) if below.any() else np.array([0])
        events.append(dict(day=round(a * dt / 86400.0, 1), surface_min=round(float(seg.min()), 3),
                           hours_below_2=round(float(below.sum() * dt / 3600.0), 1),
                           hours_below_0_5=round(float((seg < ANOXIC).sum() * dt / 3600.0), 1),
                           hours_from_minimum_back_above_2=round(float(back[0] * dt / 3600.0), 1)
                           if back.size else None))
    o = summary['oxygen'][names[m]]
    out = dict(
        years=round(years, 2), layers=bsn.n, depth_m=float(bsn.depth),
        mean_depth_m=round(float(bsn.volume.sum() / bsn.area[0]), 2),
        nights=len(ns), night_mixing_m=q(reach, (0.1, 0.25, 0.5, 0.75, 0.9)),
        # the ordinary nights, those that stop short of the floor: how deep a night mixes where no overturn takes it
        night_mixing_short_of_the_floor_m=q(reach[~to_floor], (0.1, 0.25, 0.5, 0.75, 0.9)),
        nights_mixing_to_the_floor=round(float(to_floor.mean()), 4),
        floor_ventilated_each_night=q(floor_each_night),
        floor_ventilated_in_half_the_nights=round(float(sed[layer_share >= 0.5].sum()), 4),
        floor_ventilated_in_nine_nights_of_ten=round(float(sed[layer_share >= 0.9].sum()), 4),
        floor_oxic_share_of_time=round(float(oxic @ sed), 4), floor_anoxic_share_of_time=round(float(anoxic @ sed), 4),
        floor_oxic_nine_tenths_of_the_time=round(float(sed[oxic >= 0.9].sum()), 4),
        complete_overturn_cycles=summary['mixing']['complete_overturns'], cycles=summary['cycles'],
        complete_overturn_share=summary['mixing']['complete_overturn_share'],
        overturn_episodes=len(eps), episodes_per_decade=round(len(eps) / years * 10.0, 2),
        years_between_episodes=round(float(np.diff([a for a, _ in eps]).mean() * dt / (365.25 * 86400.0)), 2)
        if len(eps) > 1 else None,
        longest_without_overturn_years=round(max(quiet) * cycle_years, 2) if quiet else 0.0,
        median_stretch_without_overturn_years=round(float(np.median(quiet)) * cycle_years, 2) if quiet else 0.0,
        deepest_tenth=dict(mean=o['deepest_tenth']['mean'], hypoxic_share=o['deepest_tenth']['hypoxic_share'],
                           anoxic_share=o['deepest_tenth']['anoxic_share']),
        bottom=dict(mean=o['bottom']['mean'], hypoxic_share=o['bottom']['hypoxic_share'],
                    anoxic_share=o['bottom']['anoxic_share']),
        bottom_t_c=summary['bottom_t_c'], surface_t_c=summary['surface_t_c'], ice_max_m=summary['ice_max_m'],
        whole_lake_mean=o['whole_lake_mean'], after_overturns=events)
    if 'nights' in o:
        n = o['nights']
        out['shallows_nights'] = dict(surface_minimum=n['surface_minimum'], hypoxic_share=n['hypoxic_share'],
                                      anoxic_share=n['anoxic_share'],
                                      surface_minimum_share_of_saturation=n['surface_minimum_share_of_saturation'],
                                      bottom_minimum=n['bottom_minimum'],
                                      hours_bottom_hypoxic=n['hours_bottom_hypoxic'],
                                      hours_bottom_anoxic=n['hours_bottom_anoxic'])
    return out, dict(centre=centre, sediment=sed, layer_night_share=layer_share, oxic=oxic, anoxic=anoxic)


def brackets(summary, trophic):
    """The deepest tenth's mean oxygen and the shallows' median night across the transfer and kinetics variants."""
    out = {}
    for key in summary['oxygen']:
        if not key.startswith(trophic + '/'):
            continue
        o = summary['oxygen'][key]
        out[key.split('/', 1)[1]] = dict(deepest_tenth_mean=o['deepest_tenth']['mean'],
                                         night_surface_minimum_median=o['nights']['surface_minimum'][1]
                                         if 'nights' in o else None)
    return out


def giant(rec, summary, m, names, sat):
    """A giant's fifty years: whether the nights or a cold spell ever mix it to the floor, and how fast its deepest
    tenth loses the oxygen it started with."""
    bsn = basin_of(rec)
    dt = float(rec['dt'])
    t_yr = (np.arange(len(rec['mu'])) + 0.5) * dt / (365.25 * 86400.0)
    deep = np.asarray(rec['deep_o2'][:, m], float)
    reach, ns = night_reach(rec)
    eps = episodes(rec, bsn)
    first = lambda level: (round(float(t_yr[np.argmax(deep < level)]), 2) if (deep < level).any() else None)
    # the loss from the oxygenated start until the deepest tenth turns hypoxic, or over the whole run if it does not
    end = first(HYPOXIC) or float(t_yr[-1])
    window = (t_yr > 0.5) & (t_yr <= end)
    slope = float(np.polyfit(t_yr[window][::24], deep[window][::24], 1)[0]) if window.sum() > 48 else float('nan')
    decade = [round(float(deep[(t_yr >= y) & (t_yr < y + 10)].mean()), 3) for y in range(0, int(t_yr[-1]) + 1, 10)
              if ((t_yr >= y) & (t_yr < y + 10)).sum() > 0]
    o = summary['oxygen'][names[m]]
    bottom_t = np.asarray(rec['bottom_t'], float)
    cyc, full = overturn_cycles(rec, bsn)
    cycle_years = float(np.mean([b - a for a, b in cyc])) * dt / (365.25 * 86400.0)
    quiet = runs_of(full, False)
    return dict(years=round(float(t_yr[-1]), 1), depth_m=float(bsn.depth), layers=bsn.n, ice_max_m=summary['ice_max_m'],
                mean_depth_m=round(float(bsn.volume.sum() / bsn.area[0]), 1),
                complete_overturn_cycles=int(full.sum()), cycles=len(cyc),
                longest_without_overturn_years=round(max(quiet) * cycle_years, 2) if quiet else 0.0,
                complete_overturn_episodes=len(eps), night_mixing_m=q(reach, (0.5, 0.9, 0.99)),
                deepest_night_mixing_m=round(float(reach.max()), 1),
                deepest_tenth_start=round(float(deep[0]), 3), deepest_tenth_end=round(float(deep[-1]), 3),
                deepest_tenth_by_decade=decade, loss_g_m3_per_year=round(-slope, 4),
                loss_fitted_over_years=[0.5, round(end, 2)],
                years_to_hypoxia=first(HYPOXIC), years_to_anoxia=first(ANOXIC),
                years_from_saturation_to_2_mg_l_at_the_loss=round((sat - HYPOXIC) / -slope, 1) if slope < 0 else None,
                bottom_t_c_first_and_last_year=[round(float(bottom_t[t_yr < 1].mean()), 3),
                                                round(float(bottom_t[t_yr > t_yr[-1] - 1].mean()), 3)],
                whole_lake_mean=o['whole_lake_mean'])
