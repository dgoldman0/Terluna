"""What the saved sea columns show: each deep column's mixing, deep water and oxygen; the shelves' and bays' nights;
the two bounds a basin's columns set on its horizontal circulation; and whether connected basins share deep water.

Each column is a case of run.py saved by columns.py: lake_oxygen's summarise() of its step-by-step record and the
whole column about every day. The deepest tenth is the deepest tenth of the basin's volume by layer centres, as the
lake model reckons it; the floor of each layer is the area the basin loses across it (lake.py's Basin.sediment).
"""
from __future__ import annotations

import math

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_oxygen import run as R
from research.studies.sea_oxygen.columns import EPISODE_DAYS
from shared.constants import MOON_RADIUS, SIDEREAL_MONTH_DAYS

HYPOXIC, ANOXIC = R.HYPOXIC, R.ANOXIC
TROPHIC = tuple(R.TROPHIC)
CENTRAL = '/renewal/supply'
YEAR_DAYS = 365.25


def rnd(x, nd=3):
    return None if x is None else round(float(x), nd)


def basin_of(z):
    return L.Basin(np.asarray(z['z'], float), np.asarray(z['area'], float))


def deepest_tenth(bsn):
    v = bsn.volume
    cum = np.cumsum(v)
    mask = (cum - 0.5 * v) >= 0.9 * cum[-1]
    if not mask.any():
        mask[-1] = True
    return np.where(mask, v, 0.0) / v[mask].sum()


def bracket_names(trophic):
    """The columns of one trophic state that bound its central one, as lake_oxygen brackets its answers: the other
    transfer velocity and the labile pool's fast and slow e-folding; the warm kinetics' transient stands beside them."""
    labile = [f'{trophic}/renewal/supply/labile_{R.LABILE_DAYS[k]:g}d' for k in ('fast', 'slow')]
    return [f'{trophic}/{tr}/supply' for tr in R.TRANSFER] + labile, [f'{trophic}/{tr}/warm' for tr in R.TRANSFER]


def bracket(res, trophic, field=('deepest_tenth', 'mean')):
    """The smallest and largest of a figure across the bounding columns of one trophic state, and across its warm
    transient."""
    def get(name):
        v = res['oxygen'][name]
        for f in field:
            v = v[f]
        return v
    bounds, warm = bracket_names(trophic)
    b, w = [get(n) for n in bounds], [get(n) for n in warm]
    return dict(bounds=[rnd(min(b)), rnd(max(b))], warm_transient=[rnd(min(w)), rnd(max(w))])


def profile_figures(path, names, sills=None):
    """Figures only the daily columns and the step-by-step series hold: the deepest tenth's age and temperature, the
    shares of the volume and of the floor under hypoxic and anoxic water, the time-mean oxygen profile's lowest point,
    the temperature at given depths (the sills to other basins), when in the decade each complete overturn begins, and
    each year's mean oxygen in the deepest tenth of each trophic state's central column."""
    with np.load(path) as z:
        bsn = basin_of(z)
        t = np.asarray(z['t'], float)
        age = np.asarray(z['age'], float)
        wd = deepest_tenth(bsn)
        out = dict(deepest_tenth_age_years=dict(mean=rnd(float((age @ wd).mean()) / YEAR_DAYS, 2),
                                                max=rnd(float((age @ wd).max()) / YEAR_DAYS, 2)),
                   column_age_years=dict(mean=rnd(float((age @ (bsn.volume / bsn.volume.sum())).mean()) / YEAR_DAYS,
                                                  2)),
                   deepest_tenth_t_c=dict(mean=rnd(float((t @ wd).mean())), min=rnd(float((t @ wd).min())),
                                          max=rnd(float((t @ wd).max()))),
                   surface_t_min_c=rnd(float(np.min(z['surface_t']))),
                   days=int(t.shape[0]), layers=int(bsn.n), oxygen={})
        v, sed = bsn.volume, bsn.sediment
        # each night's deepest mixing and the floor it reaches, as the lake kinds study reckons them
        md, mu = np.asarray(z['mixed_depth'], float), np.asarray(z['mu'], float)
        reach = np.array([md[a:b].max() for a, b in R.nights(mu)])
        share = (reach[:, None] >= bsn.centre[None, :]).mean(0)
        floor = sed / sed.sum()
        out['night_mixing_m'] = [rnd(np.quantile(reach, q), 1) for q in (0.1, 0.5, 0.9, 0.99)]
        out['floor_reached_in_half_the_nights'] = rnd(float(floor[share >= 0.5].sum()), 4)
        out['floor_reached_in_nine_nights_of_ten'] = rnd(float(floor[share >= 0.9].sum()), 4)
        # when in the decade the mixed layer reaches the bottom layer, steps within EPISODE_DAYS one overturn, as the
        # coarse passes count them (columns.recorded_passes())
        dt = float(z['dt'])
        full = np.flatnonzero(md >= float(bsn.z[-2]) - 1e-6)
        starts = full[np.r_[True, np.diff(full) > EPISODE_DAYS * 86400.0 / dt]] if full.size else full
        out['overturns_begin_years'] = [rnd(n * dt / 86400.0 / YEAR_DAYS, 2) for n in starts]
        year = int(round(YEAR_DAYS * 86400.0 / dt))
        deep_o2 = np.asarray(z['deep_o2'], float)
        out['deepest_tenth_yearly'] = {t: [rnd(deep_o2[i:i + year, names.index(t + CENTRAL)].mean())
                                           for i in range(0, deep_o2.shape[0] - year + 1, year)] for t in TROPHIC}
        oxygen = z['oxygen']
        for trophic in TROPHIC:
            m = names.index(trophic + CENTRAL)
            ox = np.asarray(oxygen[:, :, m], float)
            mean_profile = ox.mean(0)
            k = int(np.argmin(mean_profile))
            out['oxygen'][trophic] = dict(
                volume_hypoxic_share=rnd(float(((ox < HYPOXIC) @ v).mean() / v.sum()), 4),
                volume_anoxic_share=rnd(float(((ox < ANOXIC) @ v).mean() / v.sum()), 4),
                floor_hypoxic_share=rnd(float(((ox < HYPOXIC) @ sed).mean() / sed.sum()), 4),
                floor_anoxic_share=rnd(float(((ox < ANOXIC) @ sed).mean() / sed.sum()), 4),
                lowest_mean_g_m3=rnd(float(mean_profile[k])), lowest_at_m=rnd(float(bsn.centre[k]), 1),
                column_mean_g_m3=rnd(float((ox @ (v / v.sum())).mean())))
        if sills:
            out['at_sills'] = {}
            for other, depth in sills.items():
                if depth is None:
                    continue
                k = int(np.clip(np.searchsorted(bsn.centre, depth), 0, bsn.n - 1))
                out['at_sills'][other] = dict(depth_m=depth, layer_centre_m=rnd(float(bsn.centre[k]), 1),
                                              t_min_c=rnd(float(t[:, k].min())), t_mean_c=rnd(float(t[:, k].mean())))
    return out


def deep_case(res, path, sills=None):
    """A deep column's figures."""
    m = res['mixing']
    names = res['case']['scenarios']
    out = dict(name=res['case']['name'], basin=res['case']['basin'], column=res['case']['column'],
               cell=dict(lat=rnd(res['case']['lat'], 2), lon=rnd(res['case']['lon'], 2), row=res['case']['row'],
                         column=res['case']['t21_column']),
               depth_m=res['case']['depth_m'], layers=res['case']['layers'], mean_depth_m=res['case']['mean_depth_m'],
               forcing=res['case']['forcing'], cycles=res['cycles'], years=res['case']['years'],
               complete_overturns=m['complete_overturns'], years_between=m['years_between'],
               mixing_m=m['deepest_m'], bottom_age_days=m['bottom_age_days'], bottom_t_c=res['bottom_t_c'],
               surface_t_c=res['surface_t_c'], ice_max_m=res['ice_max_m'],
               spin_up=long_term(res['settling']['spin_up'], names),
               drift=dict(bottom_t_first_last_two_years=[res['settling']['bottom_t_first_two_years'],
                                                         res['settling']['bottom_t_last_two_years']]),
               oxygen={})
    for trophic in TROPHIC:
        o = res['oxygen'][trophic + CENTRAL]
        entry = dict(deepest_tenth=dict(mean=o['deepest_tenth']['mean'], p5=o['deepest_tenth']['p5'],
                                        min=o['deepest_tenth']['min'],
                                        hypoxic_share=o['deepest_tenth']['hypoxic_share'],
                                        anoxic_share=o['deepest_tenth']['anoxic_share']),
                     deepest_tenth_bracket=bracket(res, trophic), bottom_mean=o['bottom']['mean'],
                     whole_column_mean=o['whole_lake_mean'], surface=o['surface'])
        if 'overturns' in o:
            entry['after_overturns'] = dict(column_mean=o['overturns']['lake_mean_after'],
                                            surface_minimum=o['overturns']['surface_minimum_after'])
        i = names.index(trophic + CENTRAL)
        first, last = res['settling']['deep_o2_first_two_years'][i], res['settling']['deep_o2_last_two_years'][i]
        entry['deepest_tenth_first_last_two_years'] = [first, last]
        out['oxygen'][trophic] = entry
    out['profiles'] = profile_figures(path, names, sills)
    out['deep_mixing'] = res.get('deep_mixing')
    out['lake_number'] = res.get('lake_number')
    out['physics'] = res['case'].get('physics')
    return out


LONG_HALF = 0.5          # the long state is read from the last half of the spin-up's passes
TREND_PASSES = 20        # and its trends from the last twenty


def long_term(spun, names):
    """The deep water's long state from the spin-up's passes (each ten years of the record in 4-step means): over
    the last half of them, the share of passes with a complete overturn and the overturns a century, each trophic
    state's deepest-tenth oxygen and its time under 2 and 0.5 g/m3, the deepest tenth's age at the passes' ends, and
    the trends over the last twenty passes."""
    by = spun['by_pass']
    n = spun['passes']
    fine = spun.get('at_the_analysis_step') or {}
    half = slice(int(n * (1.0 - LONG_HALF)), n)
    tail = slice(max(0, n - TREND_PASSES), n)
    slope = lambda ys: float(np.polyfit(np.arange(len(ys)), ys, 1)[0]) if len(ys) > 2 else None
    over = np.asarray(by['overturns'][half], float)
    age = np.asarray(by['deep_age_days_end'][half], float) / YEAR_DAYS
    out = dict(passes=n, years=n * 10, read_from_pass=half.start, started_from=spun.get('started_from'),
               at_the_analysis_step=dict(passes=fine.get('passes'), settled=fine.get('settled'),
                                         limit=fine.get('limit'), bottom_age_days_by_pass=fine.get(
                                             'bottom_age_days_by_pass'), deep_o2_by_pass=fine.get('deep_o2_by_pass'),
                                         overturned_by_pass=fine.get('overturned_by_pass')),
               passes_with_a_complete_overturn=rnd(float((over > 0).mean()), 3),
               overturns_a_century=rnd(float(over.sum()) / (over.size / 10.0), 2),
               last_pass_with_a_complete_overturn=int(max([i for i, x in enumerate(by['overturns']) if x > 0],
                                                          default=-1)),
               deepest_tenth_age_years=dict(mean=rnd(age.mean(), 1), max=rnd(age.max(), 1),
                                            trend_per_pass=rnd(slope(list(np.asarray(by['deep_age_days_end'][tail],
                                                                                     float) / YEAR_DAYS)), 2)),
               deepest_tenth_t_c=dict(mean=rnd(np.mean(by['deep_t_end'][half])),
                                      trend_k_per_pass=rnd(slope(by['deep_t_end'][tail]), 5)),
               deepest_mixing_m=dict(median=rnd(np.median(by['deepest_mixing_m'][half]), 1),
                                     min=rnd(np.min(by['deepest_mixing_m'][half]), 1)),
               oxygen={})
    for trophic in TROPHIC:
        m = names.index(trophic + CENTRAL)
        o2 = np.asarray([x[m] for x in by['deep_o2'][half]])
        out['oxygen'][trophic] = dict(
            deepest_tenth_mean=rnd(o2.mean()), deepest_tenth_least=rnd(min(x[m] for x in by['deep_o2_min'][half])),
            hypoxic_share=rnd(np.mean([x[m] for x in by['deep_hypoxic'][half]]), 3),
            anoxic_share=rnd(np.mean([x[m] for x in by['deep_anoxic'][half]]), 3),
            trend_per_pass=rnd(slope([x[m] for x in by['deep_o2'][tail]]), 4),
            first_anoxic_pass=next((i for i, x in enumerate(by['deep_anoxic']) if x[m] > 0.5), None))
    return out


def pan_case(res):
    """A pan's nights: each night's lowest surface oxygen and, in the 10 m pans, the floor's."""
    out = dict(name=res['case']['name'], basin=res['case']['basin'], column=res['case']['column'],
               depth_m=res['case']['depth_m'], cell=dict(lat=rnd(res['case']['lat'], 2), lon=rnd(res['case']['lon'], 2)),
               forcing=res['case']['forcing'], surface_t_c=res['surface_t_c'], ice_max_m=res['ice_max_m'],
               production_g_m2_day=res['budget']['production_g_m2_day'], nights={})
    for trophic in TROPHIC:
        n = res['oxygen'][trophic + CENTRAL]['nights']
        bounds, warm = bracket_names(trophic)
        lows = [res['oxygen'][name]['nights']['surface_minimum'][1] for name in bounds]
        lows_warm = [res['oxygen'][name]['nights']['surface_minimum'][1] for name in warm]
        bottom = res['oxygen'][trophic + CENTRAL]['bottom']
        out['nights'][trophic] = dict(surface_minimum=n['surface_minimum'],
                                      share_of_saturation=n['surface_minimum_share_of_saturation'],
                                      before_dawn=n['surface_before_dawn'], hypoxic_share=n['hypoxic_share'],
                                      anoxic_share=n['anoxic_share'], bottom_minimum=n['bottom_minimum'],
                                      hours_bottom_hypoxic=n['hours_bottom_hypoxic'],
                                      hours_bottom_anoxic=n['hours_bottom_anoxic'],
                                      bottom_time_hypoxic_share=bottom['hypoxic_share'],
                                      bottom_time_anoxic_share=bottom['anoxic_share'],
                                      bracket_of_medians=[rnd(min(lows)), rnd(max(lows))],
                                      warm_transient_medians=[rnd(min(lows_warm)), rnd(max(lows_warm))])
    out['physics'] = res['case'].get('physics')
    return out


def pan_bracket(pans, suffix='_fixed_dissipation'):
    """Each bracketed pan beside its central run: under the wind's share of P10 spread over the thin water below its
    mixed layer, and under the first round's fixed dissipation and cap, the nights' surface minimum (median), the share
    of nights under 2 g/m3 at the surface, the floor's hours under 2 and 0.5 g/m3 a night (median and 90th
    percentile) and the floor's time under them."""
    out = {}
    for name, v in pans.items():
        if not name.endswith(suffix) or name[:-len(suffix)] not in pans:
            continue
        base = pans[name[:-len(suffix)]]
        entry = {}
        for trophic in TROPHIC:
            row = {}
            for label, p in (('wind', base), ('fixed_dissipation', v)):
                n = p['nights'][trophic]
                row[label] = dict(surface_minimum_median=n['surface_minimum'][1],
                                  surface_hypoxic_nights=n['hypoxic_share'],
                                  bottom_minimum_median=n['bottom_minimum'][1],
                                  bottom_hypoxic_hours=[n['hours_bottom_hypoxic'][1], n['hours_bottom_hypoxic'][2]],
                                  bottom_anoxic_hours=[n['hours_bottom_anoxic'][1], n['hours_bottom_anoxic'][2]],
                                  bottom_time_hypoxic_share=n['bottom_time_hypoxic_share'],
                                  bottom_time_anoxic_share=n['bottom_time_anoxic_share'])
            entry[trophic] = row
        out[name[:-len(suffix)]] = dict(depth_m=v['depth_m'], basin=v['basin'], nights=entry)
    return out


def lake_pans(nights_of):
    """The lakes' 3 m pans in the same terms, from lake_oxygen's summaries (its run_case())."""
    out = {}
    for band, res in nights_of.items():
        entry = {}
        for trophic in TROPHIC:
            n = res['oxygen'][trophic + CENTRAL]['nights']
            entry[trophic] = dict(surface_minimum=n['surface_minimum'], hypoxic_share=n['hypoxic_share'],
                                  anoxic_share=n['anoxic_share'],
                                  share_of_saturation=n['surface_minimum_share_of_saturation'])
        out[band] = entry
    return out


def demand_deepest_tenth(case_res, path):
    """The deepest tenth's reference oxygen demand (g/m3 a year) for each trophic state's central column, from the lake
    model's own demand profile on the column's layers: the demand the central columns keep, fixed by supply, at any
    warmth (warm_factor() gives how far the warm kinetics would raise it)."""
    with np.load(path) as z:
        bsn = basin_of(z)
    wd = deepest_tenth(bsn)
    out = {}
    for s in R.scenarios():
        if s.name.endswith(CENTRAL):
            out[s.name.split('/')[0]] = rnd(float(L.demand_profile(bsn, s) @ wd) * YEAR_DAYS, 4)
    return out


COLD_END_MARGINS_K = (0.5, 1.0, 2.0)


def cold_end_shares(cells_of_basin, coldest, margins=COLD_END_MARGINS_K):
    """The share of a basin's water in cells whose coldest 3-day mean lies within each margin of the coldest column's
    cell's (coldest, C), the cell setting.RULE takes among those holding at least 1,000 km2 of the basin."""
    total = sum(c['area_km2'] for c in cells_of_basin)
    return {f'{m:g}': rnd(sum(c['area_km2'] for c in cells_of_basin if c['min'] <= coldest + m) / total, 4)
            for m in margins}


def warm_factor(t_c):
    """How many times the warm kinetics' demand exceeds the reference demand at a deep water's temperature: ecosystem
    respiration's activation energy from the temperature the rates were measured at (lake_oxygen's KINETICS['warm']
    and SOURCE_C); the central columns keep the reference demand, fixed by supply, whatever the warmth."""
    ev = R.KINETICS['warm']
    return float(np.exp(ev / L.BOLTZMANN_EV * (1.0 / (R.SOURCE_C + 273.15) - 1.0 / (t_c + 273.15))))


def screen(interval, shares, coldest, demand):
    """The screen of exchange on one interval between the coldest column's renewals (years): at each margin's share
    of the basin, the interval stretched over that share, and the deep water's mean oxygen, the coldest column's less
    the deepest tenth's reference demand over half the added interval (a linear decline, floored at zero)."""
    out = dict(renewal_interval_coldest_years=rnd(interval, 2), margins={})
    for margin, share in shares.items():
        interval_x = interval / share if interval and share else None
        oxygen = {}
        for trophic in TROPHIC:
            o_cold = coldest['oxygen'][trophic]['deepest_tenth']['mean']
            r = demand.get(trophic)
            oxygen[trophic] = None if interval_x is None or r is None else \
                rnd(max(0.0, o_cold - r * 0.5 * (interval_x - interval)))
        out['margins'][margin] = dict(cold_end_share=share, renewal_interval_with_exchange_years=rnd(interval_x, 2),
                                      deepest_tenth_mean_with_exchange=oxygen)
    return out


def overturn_counts(case):
    """A column's complete overturns, the mixed layer reaching the bottom layer: in its saved decade, as overturns
    (steps within EPISODE_DAYS one, as the coarse passes count them) and their rate a century, and as the lunar cycles
    an overturn touches (lake_oxygen's complete_overturns; an overturn spanning two cycles counts in both); and the
    coarse passes' last half's overturns a century."""
    episodes = len(case['profiles']['overturns_begin_years'])
    return dict(saved_decade=episodes, saved_decade_a_century=rnd(episodes * 100.0 / case['years'], 1),
                saved_decade_cycles=case['complete_overturns'], cycles=case['cycles'],
                coarse_a_century=case['spin_up']['overturns_a_century'])


def bounds(basin, typical, coldest, shares):
    """What a basin's typical and coldest columns bound, and what horizontal exchange would do to them: their deep
    water's temperature, age and oxygen in the saved decade, and their overturns in the saved decade and in the coarse
    passes' last half. The coldest column makes the basin's coldest deep water and sets the floor to its temperature;
    the typical column sets it from above; each gives the renewal its own weather gives a column of the whole basin.
    With exchange the coldest end's water sinks under the rest: the deep water takes its temperature, which the typical
    part's surface reaches or never reaches (typical_surface_reaches_cold_deep_water), and only the cold end renews it,
    the share of the basin's area whose coldest 3-day means lie within a margin of the coldest cell's. The screen
    spreads the coldest column's renewal over that share and takes the deepest tenth's reference demand over the added
    interval; it reads the interval from the saved decade's overturns, as the oxygen is read, with the coarse passes'
    beside it (coarse). It scales the two columns and holds no circulation. The reference demand is the central
    columns' (fixed by supply); the warm kinetics would raise it warm_kinetics_factor times at the deep water's warmth."""
    t_typ, t_cold = typical['profiles']['deepest_tenth_t_c']['mean'], coldest['profiles']['deepest_tenth_t_c']['mean']
    n = len(coldest['profiles']['overturns_begin_years'])
    interval = coldest['years'] / n if n else None
    per_century = coldest['spin_up']['overturns_a_century']
    interval_coarse = 100.0 / per_century if per_century else None
    demand = coldest.get('deep_demand_g_m3_yr') or {}
    out = dict(basin=basin, deep_t_c=dict(typical=t_typ, coldest=t_cold, difference_k=rnd(t_typ - t_cold)),
               overturns=dict(typical=overturn_counts(typical), coldest=overturn_counts(coldest)),
               deep_age_years=dict(typical=typical['profiles']['deepest_tenth_age_years']['mean'],
                                   coldest=coldest['profiles']['deepest_tenth_age_years']['mean']),
               cold_end_area_share=shares,
               typical_surface_min_c=typical['profiles']['surface_t_min_c'],
               typical_surface_reaches_cold_deep_water=bool(typical['profiles']['surface_t_min_c'] < t_cold),
               coldest_surface_c=coldest['surface_t_c'],
               warm_kinetics_factor=dict(typical=rnd(warm_factor(t_typ), 2), coldest=rnd(warm_factor(t_cold), 2)),
               oxygen={trophic: dict(typical=dict(mean=typical['oxygen'][trophic]['deepest_tenth']['mean'],
                                                  anoxic_share=typical['oxygen'][trophic]['deepest_tenth'][
                                                      'anoxic_share']),
                                     coldest=dict(mean=coldest['oxygen'][trophic]['deepest_tenth']['mean'],
                                                  anoxic_share=coldest['oxygen'][trophic]['deepest_tenth'][
                                                      'anoxic_share']),
                                     deepest_tenth_reference_demand_g_m3_yr=demand.get(trophic))
                       for trophic in TROPHIC},
               screen=dict(read_from='the saved decade\'s overturns',
                           **screen(interval, shares, coldest, demand),
                           coarse=dict(read_from='the coarse passes\' last half',
                                       **screen(interval_coarse, shares, coldest, demand))))
    return out


OMEGA = 2.0 * math.pi / (SIDEREAL_MONTH_DAYS * 86400.0)   # rad/s: the Moon turns once a sidereal month


def internal_wave_speed(path, gravity):
    """The first internal mode's speed through a column's saved decade, taken as two layers: each day
    c = (g' h)^(1/2), h the day's mean mixed depth and g' = g (rho_deep - rho_mixed) / rho0 between the mixed layer
    and the deepest tenth (the lake model's density, as the Lake number takes it); days with no density step, the
    column overturned, are left out. The 10th, 50th and 90th percentiles over the days, with the density step's
    temperature difference and the mixed depth at the median."""
    with np.load(path) as z:
        bsn = basin_of(z)
        t = np.asarray(z['t'], float)
        steps = np.asarray(z['step'], int)
        md = np.asarray(z['mixed_depth'], float)
        every = int(z['every'])
    v, wd = bsn.volume, deepest_tenth(bsn)
    rho = L.density(t)
    c, dt_k, hs = [], [], []
    for k, n in enumerate(steps):
        h = float(md[max(0, n - every + 1):n + 1].mean())
        above = bsn.centre < h
        if not above.any():
            continue
        rho_mixed = float((rho[k][above] * v[above]).sum() / v[above].sum())
        g_reduced = gravity * (float(rho[k] @ wd) - rho_mixed) / L.RHO0
        if g_reduced <= 0.0:
            continue
        c.append(math.sqrt(g_reduced * h))
        dt_k.append(float((t[k][above] * v[above]).sum() / v[above].sum()) - float(t[k] @ wd))
        hs.append(h)
    c, dt_k, hs = np.asarray(c), np.asarray(dt_k), np.asarray(hs)
    return dict(c_m_s=[rnd(np.quantile(c, q), 3) for q in (0.1, 0.5, 0.9)],
                temperature_step_k_median=rnd(np.median(dt_k), 2), mixed_depth_m_median=rnd(np.median(hs), 1),
                days=int(c.size), of_days=int(steps.size))


def rossby_radius_km(c, lat_deg):
    """The internal Rossby radius (km) of a wave speed c (m/s) at a latitude: c/|f| with f = 2 Omega sin(lat), or near
    the equator the equatorial radius (c/beta)^(1/2), beta = 2 Omega cos(lat)/R, where that is the smaller."""
    lat = math.radians(lat_deg)
    f = 2.0 * OMEGA * abs(math.sin(lat))
    beta = 2.0 * OMEGA * math.cos(lat) / MOON_RADIUS
    mid = c / f if f > 0.0 else math.inf
    return min(mid, math.sqrt(c / beta)) / 1000.0


def rotation(waves, basin):
    """Whether the Moon's rotation shapes a basin's response to the wind: the internal Rossby radius at the median wave
    speed of its typical column, at the latitudes holding the basin's area's 5th, 50th and 95th percentiles (and at the
    median latitude for the 90th percentile speed, the deeper mixed layers), against the basin's own scale, the square
    root of its area. A basin several radii across answers the wind by Ekman transport,
    upwelling along its shores within a radius of them and Kelvin and Poincare waves; one about a radius across tilts
    as a lake does."""
    c, c90 = waves['c_m_s'][1], waves['c_m_s'][2]
    lats = basin['lat_quantiles_by_area']
    radii = {q: rnd(rossby_radius_km(c, lat), 0) for q, lat in lats.items()}
    scale = math.sqrt(basin['area_km2'])
    return dict(wave=waves, latitudes=lats, rossby_radius_km=radii,
                rossby_radius_km_at_the_90th_percentile_speed=rnd(rossby_radius_km(c90, lats['0.5']), 0),
                basin_scale_km=rnd(scale, 0), scale_in_radii=rnd(scale / radii['0.5'], 1) if radii['0.5'] else None)


def connected(entries, pairs):
    """Whether a basin's water at its sill to another basin of the same sea is ever colder than that basin's deep
    water, so that it could flow over the sill and sink into it. pairs: (source case, receiving case, the receiving
    basin's key, the sill's depth)."""
    out = []
    for a, b, b_key, depth in pairs:
        if depth is None:
            continue
        at = entries[a]['profiles'].get('at_sills', {}).get(b_key)
        if not at:
            continue
        deep_b = entries[b]['profiles']['deepest_tenth_t_c']['mean']
        out.append(dict(source=a, into=b, sill_m=depth, source_t_min_at_sill_c=at['t_min_c'],
                        source_t_mean_at_sill_c=at['t_mean_c'], deep_t_c_of_receiver=deep_b,
                        colder_than_receivers_deep_water=bool(at['t_min_c'] < deep_b),
                        margin_k=rnd(deep_b - at['t_min_c'])))
    return out


def saturation_at(st, share):
    """The seas' oxygen at saturation in a target between the two runs, as lake_oxygen's saturation() reckons each of
    its targets, at a share of the way: the seas' median and band temperatures at 1.2 atm, and the distribution of
    every water cell's 3-day means at its own pressure."""
    p_design = R.DESIGN_PRESSURE_PA / 101325.0
    lunar = lambda t, s=None: float(L.oxygen_saturation(t, p_design, st.salinity if s is None else s, st.o2_fraction))
    runs = st.targets['runs']
    mix = lambda a, b: a + share * (b - a)
    design, test = runs['design']['seas'], runs['test']['seas']
    median_c = mix(runs['design']['sea_median_c'], runs['test']['sea_median_c'])
    out = dict(share_of_the_way=share, sea_median_c=rnd(median_c, 2), mg_l=rnd(lunar(median_c)),
               seawater_mg_l=rnd(lunar(median_c, 35.0)),
               warm_99th_c=rnd(mix(design['all']['p99'], test['all']['p99']), 2),
               warm_99th_mg_l=rnd(lunar(mix(design['all']['p99'], test['all']['p99']))), bands={})
    for band in design['bands']:
        t_band = mix(design['bands'][band]['median'], test['bands'][band]['median'])
        out['bands'][band] = dict(median_c=rnd(t_band, 2), mg_l=rnd(lunar(t_band)))
    stats = {}
    for name in ('lock', 'hot'):
        ts = st.f[f'{name}__waters_ts'] - 273.15
        ps = st.f[f'{name}__waters_ps'] * 100.0 / 101325.0
        c = L.oxygen_saturation(ts, ps, st.salinity, st.o2_fraction)
        w = np.broadcast_to(st.f[f'{name}__waters_weight'], c.shape)
        sea = np.broadcast_to(st.f[f'{name}__waters_surface_m'] < 50.0, c.shape)
        stats[name] = dict(median=R.weighted_quantile(c.ravel(), w.ravel(), 0.5),
                           p1=R.weighted_quantile(c.ravel(), w.ravel(), 0.01),
                           p99=R.weighted_quantile(c.ravel(), w.ravel(), 0.99),
                           seas=R.weighted_quantile(c[sea], w[sea], 0.5),
                           seas_pressure_atm=float(np.average(ps[sea], weights=w[sea])))
    lw, hw = stats['lock'], stats['hot']
    out['three_day_means'] = dict(median=rnd(mix(lw['median'], hw['median'])), p1=rnd(mix(lw['p1'], hw['p1'])),
                                  p99=rnd(mix(lw['p99'], hw['p99'])), seas_median=rnd(mix(lw['seas'], hw['seas'])),
                                  seas_pressure_atm=rnd(lw['seas_pressure_atm'], 4))
    return out


RUN_PASSES = 3           # a last pass that overturns after this many passes without one renews a column cut off


def analysis_step(case):
    """What a column's passes at the analysis step say of its saved decade, which follows them. A column that settled
    gives its long state. One that did not leans one way: still losing oxygen in its coupled column, its oxygen is an
    upper value; still gaining, a lower value; its last pass overturning after a run of passes without one, its decade
    is the column just renewed, its oxygen an upper value and its age a lower value; ageing with no overturn at all, its
    age is a lower value. Each pass is the ten-year record; the coupled column is the trophic state whose physics the
    column runs (its deepest tenth's mean oxygen by pass), and the bottom water's age is its mean over the pass."""
    fine = case['spin_up']['at_the_analysis_step']
    over = [bool(x) for x in fine['overturned_by_pass']]
    o2 = [float(x) for x in fine['deep_o2_by_pass']]
    age = [float(x) / YEAR_DAYS for x in fine['bottom_age_days_by_pass']]
    runs, run = [], 0
    for x in over:
        run = 0 if x else run + 1
        runs.append(run)
    per_pass = case['years']
    out = dict(passes=fine['passes'], settled=bool(fine['settled']), passes_overturning=int(sum(over)),
               longest_run_without_a_complete_overturn_years=rnd(max(runs) * per_pass, 1),
               last_pass_overturned=over[-1],
               bottom_age_years=dict(first=rnd(age[0], 1), oldest=rnd(max(age), 1), last=rnd(age[-1], 1)),
               coupled=case['physics']['coupled'],
               coupled_deepest_tenth_o2=dict(first=rnd(o2[0]), least=rnd(min(o2)), most=rnd(max(o2)), last=rnd(o2[-1]),
                                             last_change=rnd(o2[-1] - o2[-2]) if len(o2) > 1 else None))
    change = o2[-1] - o2[-2] if len(o2) > 1 else 0.0
    if fine['settled']:
        out.update(reading='settled', oxygen_is=None, age_is=None)
    elif over[-1] and len(runs) > 1 and runs[-2] >= RUN_PASSES:
        out.update(reading=(f'its decade follows its one overturn at the analysis step, after {runs[-2]} passes '
                            f'without one, and holds the column just renewed'), oxygen_is='an upper value',
                   age_is='a lower value')
    elif not any(over):
        out.update(reading='no overturn at the analysis step: cut off from the surface and ageing',
                   oxygen_is='its deep water anoxic' if max(o2) < ANOXIC else 'an upper value', age_is='a lower value')
    elif change < -R.SETTLED_O2:
        out.update(reading='still losing oxygen', oxygen_is='an upper value', age_is=None)
    elif change > R.SETTLED_O2:
        out.update(reading='still gaining oxygen', oxygen_is='a lower value', age_is=None)
    else:
        out.update(reading='unsettled in its age or solute', oxygen_is=None, age_is=None)
    return out


def _deep_answer(case):
    """One column's answer. Its deep water's oxygen is the saved decade's, after the passes at the analysis step: the
    deepest tenth's mean with its bracket (the other transfer and the labile pool's e-folding) and its time hypoxic and
    anoxic, the shares of the volume and floor under hypoxic water; analysis_step says whether those passes settled
    and, where they did not, which way the decade's figures lean. The overturns are the saved decade's, counted as the
    coarse passes count them (steps within EPISODE_DAYS one) and as the lunar cycles they touch, and a century the
    coarse passes' last half's; the long view is the coarse passes' last half: the deepest tenth's age and its oxygen
    there, which the coarse steps put higher than the analysis step does."""
    lt = case['spin_up']
    prof = case['profiles']
    out = dict(settled_at_the_analysis_step=lt['at_the_analysis_step']['settled'],
               passes_at_the_analysis_step=lt['at_the_analysis_step']['passes'],
               analysis_step=analysis_step(case),
               decade_overturns=len(prof['overturns_begin_years']),
               decade_overturns_a_century=rnd(len(prof['overturns_begin_years']) * 100.0 / case['years'], 1),
               decade_cycles_overturning=case['complete_overturns'],
               overturns_a_century=lt['overturns_a_century'],
               deep_age_years=prof['deepest_tenth_age_years']['mean'],
               long_deep_age_years=lt['deepest_tenth_age_years']['mean'],
               deep_t_c=prof['deepest_tenth_t_c']['mean'], night_mixing_m=prof.get('night_mixing_m'),
               floor_reached_in_half_the_nights=prof.get('floor_reached_in_half_the_nights'), oxygen={})
    for trophic in TROPHIC:
        d = case['oxygen'][trophic]
        p = prof['oxygen'][trophic]
        o = lt['oxygen'][trophic]
        out['oxygen'][trophic] = dict(deepest_tenth_mean=d['deepest_tenth']['mean'],
                                      bracket=d['deepest_tenth_bracket']['bounds'],
                                      hypoxic_share=d['deepest_tenth']['hypoxic_share'],
                                      anoxic_share=d['deepest_tenth']['anoxic_share'],
                                      volume_hypoxic_share=p['volume_hypoxic_share'],
                                      floor_hypoxic_share=p['floor_hypoxic_share'],
                                      decade_first_last_two_years=d['deepest_tenth_first_last_two_years'],
                                      long_deepest_tenth_mean=o['deepest_tenth_mean'],
                                      long_anoxic_share=o['anoxic_share'])
    return out


def answers(deep, variants, cells, stands_for):
    """Each basin's deep water: its typical column at the central wind share, the variants run at its typical cell
    (the wind share's bracket, and at the two largest basins the sensitivity rows and the other form of the deep
    mixing) keyed typical_<variant>, and its coldest column at the central share."""
    out = {}
    for key, ch in cells.items():
        entry = dict(stands_for=stands_for[key]['share_of_the_seas'], one_column=ch['one_column'],
                     typical=_deep_answer(deep[f'{key}_typical']))
        for name in sorted(variants):
            if name.startswith(f'{key}_typical_'):
                entry[name[len(key) + 1:]] = _deep_answer(variants[name])
        if f'{key}_coldest' in deep:
            entry['coldest'] = _deep_answer(deep[f'{key}_coldest'])
        out[key] = entry
    return out


SEAS_VIEWS = ('typical', 'typical_wind_share_low', 'typical_wind_share_high', 'typical_shih_bouffard', 'coldest')


def seas_wide(answers, views=SEAS_VIEWS):
    """The seas' deep water by area, each run basin standing for its share of the seas (seas.json's stands_for): for
    each trophic state the area-weighted mean of the deepest tenth's decade mean, and the shares of the seas' area
    whose deepest tenth averages under 2 and under 0.5 g/m3. Each view takes one column of every basin, the typical
    one where the basin lacks it (a basin with one column has no coldest; the other form ran at two basins), and
    says how much of the seas' area its own column covers and which of the columns it takes did not settle at the
    analysis step, with the way each leans (analysis_step()), so that a view holding an upper value is one."""
    out = {}
    for view in views:
        w, own, rows, unsettled = [], 0.0, {t: [] for t in TROPHIC}, {}
        for key, e in answers.items():
            pick = e.get(view)
            own += e['stands_for'] if pick is not None else 0.0
            pick = pick or e['typical']
            w.append(e['stands_for'])
            for t in TROPHIC:
                rows[t].append(pick['oxygen'][t]['deepest_tenth_mean'])
            if not pick['analysis_step']['settled']:
                unsettled[key] = dict(reading=pick['analysis_step']['reading'],
                                      oxygen_is=pick['analysis_step']['oxygen_is'], stands_for=e['stands_for'])
        w = np.asarray(w, float)
        entry = dict(own_columns_share=rnd(own / w.sum(), 4), unsettled_columns=unsettled)
        for t in TROPHIC:
            o = np.asarray(rows[t], float)
            entry[t] = dict(mean=rnd(float(o @ w / w.sum())), area_under_2=rnd(float(w[o < HYPOXIC].sum() / w.sum()), 4),
                            area_under_0_5=rnd(float(w[o < ANOXIC].sum() / w.sum()), 4))
        out[view] = entry
    return out


def forms(answers, form='typical_shih_bouffard'):
    """The two forms of the deep mixing at the basins where both ran, at their typical cells: Osborn's rule at the
    central wind share with the share's bracket beside it, and the other form at the central share; the decade's
    overturns, the deepest tenth's age and temperature and each trophic state's deepest-tenth oxygen and its time
    anoxic, and each figure's range between the forms."""
    out = {}
    for key, e in answers.items():
        if form not in e:
            continue
        rows = dict(osborn=e['typical'], osborn_wind_share_low=e.get('typical_wind_share_low'),
                    osborn_wind_share_high=e.get('typical_wind_share_high'), shih_bouffard=e[form])
        entry = {label: None if r is None else dict(
            decade_overturns=r['decade_overturns'], decade_overturns_a_century=r['decade_overturns_a_century'],
            decade_cycles_overturning=r['decade_cycles_overturning'], overturns_a_century=r['overturns_a_century'],
            deep_age_years=r['deep_age_years'], long_deep_age_years=r['long_deep_age_years'], deep_t_c=r['deep_t_c'],
            settled_at_the_analysis_step=r['settled_at_the_analysis_step'], analysis_step=r['analysis_step'],
            oxygen={t: dict(deepest_tenth_mean=r['oxygen'][t]['deepest_tenth_mean'],
                            anoxic_share=r['oxygen'][t]['anoxic_share'],
                            volume_hypoxic_share=r['oxygen'][t]['volume_hypoxic_share']) for t in TROPHIC})
            for label, r in rows.items()}
        a, b = e['typical'], e[form]
        entry['range'] = dict(
            deep_age_years=sorted([a['deep_age_years'], b['deep_age_years']]),
            oxygen={t: sorted([a['oxygen'][t]['deepest_tenth_mean'], b['oxygen'][t]['deepest_tenth_mean']])
                    for t in TROPHIC})
        out[key] = entry
    return out
