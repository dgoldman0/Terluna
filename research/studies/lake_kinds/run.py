"""The Moon's lakes sorted by kind at the resort's 301 K, with what each kind needs.

    python -m research.studies.lake_kinds.ladder      # the lake runs at 301 K (about 2.5 hours on two cores)
    python -m research.studies.lake_kinds.run         # the census and the product (a few minutes)

It reads the lake runs ladder.py saved in research/runs/lake_kinds (their index ladder.json holds the sha256 of every
saved file, which this checks), the waters' oxygen study's product, the anoxic depths study's gas physics and remedies
(run on this study's 301 K columns), the census's lake set (census.py: geography's drainage on LOLA, the climatology
it routed and groundwater stage 2), the sea chemistry's left-alone path, the springs' water and climate_targets, and
writes results/lake_kinds.json:
- the 301 K resort and the proof of its construction at 300 K;
- the ladder: lakes with floors from 5 to 352 m, each trophic state with its own physics, and what their nights, their
  overturns and their floors do; the depth thresholds of the kinds; the giants' fifty years; the checks;
- the census: every body of the drainage product by kind, band, size, outlet and groundwater, with the floor each
  kind's nights ventilate and its area by floor depth for the methane study to multiply;
- the closed lakes' path to soda and the screen for meromixis, and the springs' density against the lakes';
- what each kind needs, sized, and the overturn events of the nutrient-poor lakes with the refuges their fish need.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
import time
from pathlib import Path
from unittest import mock

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_anoxia import gases as G
from research.studies.lake_anoxia import run as A
from research.studies.lake_kinds import analysis as AN
from research.studies.lake_kinds import census as CE
from research.studies.lake_kinds import columns as C
from research.studies.lake_kinds import ladder as LD
from research.studies.lake_kinds import setting as S
from research.studies.lake_oxygen import run as R
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'lake_kinds.json'
SCHEMA = 'terluna.research.lake-kinds/1'
SEA_CHEMISTRY, SPRINGS = CE.SEA_CHEMISTRY, CE.SPRINGS
TROPHIC = LD.TROPHIC
BANDS = LD.BAND_NAMES
CENTRAL_TROPHIC = 'oligotrophic'       # the state the land's dissolved phosphorus holds the lakes in (lake_anoxia),
                                       # which the census's central sorting takes; the others are sorted beside it
# The kinds by the deepest point: floors the nights reach, then the edges the brief and the ladder set.
KIND_EDGES_M = dict(middle=50.0, deep=1000.0)
KINDS = ('self_ventilated', 'middle', 'deep', 'giant')
KIND_NAMES = dict(self_ventilated='self-ventilated: deepest point within the ventilated depth, the floor the nights '
                                  'reach on half of them',
                  middle='the middle lakes: beyond the ventilated depth to 50 m',
                  deep='deep: beyond the ventilated depth and 50 m, to 1,000 m',
                  giant='giants: deeper than 1,000 m')
FLOOR_EDGES_M = A.DEPTH_EDGES          # lake_anoxia's floor classes, for the methane study's fluxes by floor depth
SALINITY_DENSITY, MARGINS_KG_M3, FRESHENING = CE.SALINITY_DENSITY, CE.MARGINS_KG_M3, CE.FRESHENING
EARTH_MOORINGS = dict(lake='Lake Erie, its moorings across the central basin', moorings=10, sensors_per_mooring=(6, 8),
                      lake_area_km2=25700.0, source='NOAA OxyMoorON; Michigan Sea Grant for the area')
MOORING_KM2 = 2000.0                   # proposed design: a mooring for each further 2,000 km2 of a large lake
SMALL_KM2 = 1.0                        # lakes under this take profiles from boats or a sensor string
OXYGENATION_SHARES = (0.01, 0.1)       # shares of the deep lakes' area chosen for deep-water life, to size the power


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def dumps(obj) -> str:
    text = json.dumps(obj, indent=1, allow_nan=False)
    return re.sub(r'\[\s+([-0-9.eE,\s]+?|null|true|false)\s+\]', lambda m: '[' + ' '.join(m.group(1).split()) + ']',
                  text)


def rnd(x, sig=4):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, np.ndarray):
        return [rnd(v, sig) for v in x.tolist()]
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        x = x.item()
    if isinstance(x, bool) or not isinstance(x, float):
        return x
    if not math.isfinite(x):
        return None
    if x == 0.0:
        return 0.0
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


# --- the runs ---------------------------------------------------------------------------------------------------

class Runs:
    """ladder.py's index, every saved summary and column checked against it."""

    def __init__(self):
        self.index = json.loads(LD.INDEX.read_text())
        if self.index['schema'] != LD.SCHEMA:
            raise SystemExit(f"{LD.INDEX} is {self.index['schema']}")
        self.summaries = {}
        for name, entry in self.index['cases'].items():
            path = C.CASES / entry['summary']
            if sha256(path) != entry['summary_sha256']:
                raise SystemExit(f'{path.name} does not match ladder.json; run ladder.py again')
            self.summaries[name] = json.loads(path.read_text())
        self.pans = {}
        for band, entry in self.index['pans'].items():
            path = C.CASES / entry['summary']
            if sha256(path) != entry['summary_sha256']:
                raise SystemExit(f'{path.name} does not match ladder.json')
            self.pans[band] = json.loads(path.read_text())
        self.checked = set()

    def record(self, name, keys=None):
        """The case's saved passes, joined, with each file checked against its sha256 the first time it is read;
        keys limits the arrays read (the basin, the step and the time step always come)."""
        entry = self.index['cases'][name]
        parts = []
        for f, h in entry['files'].items():
            path = C.CASES / f
            if f not in self.checked:
                if sha256(path) != h:
                    raise SystemExit(f'{f} does not match ladder.json; run ladder.py again')
                self.checked.add(f)
            with np.load(path) as z:
                want = z.files if keys is None else [k for k in z.files if k in set(keys) | {'z', 'area', 'dt',
                                                                                             'every', 'step', 'mu'}]
                parts.append({k: z[k] for k in want})
        if len(parts) == 1:
            return parts[0]
        out = dict(parts[0])
        joined = ('mu', 'pressure') + C.SERIES + C.COLUMN_SERIES + ('t', 'oxygen', 'debt', 'age', 'solute')
        for k in joined:
            if k in out:
                out[k] = np.concatenate([p[k] for p in parts])
        n = [len(p['mu']) for p in parts]
        out['step'] = np.concatenate([p['step'] + sum(n[:i]) for i, p in enumerate(parts)])
        return out


class Gases(A.Setting):
    """lake_anoxia's setting, its gas physics and remedies read on this study's columns."""

    def load(self, name):
        return dict(np.load(C.CASES / name))


def own_column(summary, trophic):
    names = summary['case']['scenarios']
    return names.index(f'{trophic}/renewal/supply'), names


# --- the ladder -------------------------------------------------------------------------------------------------

def bubble_years(gas, depth, t_c, trophic):
    """lake_anoxia's isolated hypolimnion: its trophic state's floor and deep-water demand over a hypolimnion reaching
    from the floor to 15 m below the surface, or 50 m at most, 5 m at least, at the lake's own deep temperature."""
    thick = min(max(depth - 15.0, 5.0), 50.0)
    p = R.TROPHIC[trophic]
    rate = (p['sediment'] + p['deep_water']) / thick * 365.25
    return A.isolated_years(gas, depth, t_c, rate)


def event_summary(events):
    """The overturning lunar months of a run: how many, how many took the surface below 2 mg/L, and how low and for
    how long."""
    if not events:
        return dict(events=0, events_below_2=0)
    smin = np.array([e['surface_min'] for e in events])
    h2 = np.array([e['hours_below_2'] for e in events])
    h05 = np.array([e['hours_below_0_5'] for e in events])
    below = smin < R.HYPOXIC
    return dict(events=len(events), events_below_2=int(below.sum()),
                surface_min=dict(lowest=float(smin.min()), median=float(np.median(smin))),
                hours_below_2=dict(median_of_those_below=float(np.median(h2[below])) if below.any() else 0.0,
                                   most=float(h2.max())),
                hours_below_0_5=dict(most=float(h05.max())))


def ladder(runs, gas, eta='central'):
    """Each band, floor depth and trophic state, with its own physics, at one share of the wind's energy below the
    mixed layer (the depths and trophic states the ladder ran at it)."""
    out = {}
    layers = {}
    for b, band in enumerate(BANDS):
        out[band] = {}
        for depth in LD.LADDER_M:
            out[band][f'{depth:g}'] = {}
            for trophic in TROPHIC:
                name = LD.ladder_name(b, depth, trophic, eta)
                if name not in runs.summaries:
                    continue
                s = runs.summaries[name]
                m, names = own_column(s, trophic)
                rec = runs.record(name)
                stats, lay = AN.column(rec, s, m, names)
                stats['brackets'] = AN.brackets(s, trophic)
                if depth >= 16.0:
                    stats['years_to_bubble_point_if_isolated'] = bubble_years(gas, depth, s['bottom_t_c']['mean'],
                                                                              trophic)
                stats['budget'] = s['budget']
                stats['deep_mixing'] = s['deep_mixing']
                spin = s['settling']['spin_up']
                fine = spin.get('at_the_analysis_step') or {}
                stats['spin_up_passes'] = spin['passes']
                stats['passes_at_the_analysis_step'] = dict(passes=fine.get('passes'), settled=fine.get('settled'))
                events = stats.pop('after_overturns')
                stats['after_overturns'] = event_summary(events)
                out[band][f'{depth:g}'][trophic] = stats
                layers[(band, depth, trophic)] = dict(lay, events=events)
                del rec
    return out, layers


def lake_oxygen_convention(runs):
    """The 352 m lakes at 301 K with lake_oxygen's convention, every trophic state's oxygen under mesotrophic physics,
    beside the waters' oxygen study's resort at 300 K."""
    product = json.loads(R.OUT.read_text())
    out = {}
    for b, band in enumerate(BANDS):
        s = runs.summaries[LD.ladder_name(b, R.DEEP_M, 'mesotrophic')]
        old = product['runs']['deep'][f'resort/{band}']
        row = dict(complete_overturns=s['mixing']['complete_overturns'], cycles=s['cycles'],
                   mixing_m=s['mixing']['deepest_m'], bottom_age_mean_days=s['mixing']['bottom_age_days']['mean'],
                   bottom_t_c=s['bottom_t_c']['mean'],
                   at_300k=dict(complete_overturns=old['mixing']['complete_overturns'],
                                mixing_m=old['mixing']['deepest_m'],
                                bottom_age_mean_days=old['mixing']['bottom_age_days']['mean'],
                                bottom_t_c=old['bottom_t_c']['mean']), traits={})
        for trophic in TROPHIC:
            o, p = s['oxygen'][f'{trophic}/renewal/supply'], old['oxygen'][f'{trophic}/renewal/supply']
            row['traits'][trophic] = dict(
                deepest_tenth_mean=o['deepest_tenth']['mean'], anoxic_share=o['deepest_tenth']['anoxic_share'],
                lake_mean_after_overturns=(o.get('overturns') or {}).get('lake_mean_after'),
                at_300k=dict(deepest_tenth_mean=p['deepest_tenth']['mean'],
                             anoxic_share=p['deepest_tenth']['anoxic_share'],
                             lake_mean_after_overturns=(p.get('overturns') or {}).get('lake_mean_after')))
        out[band] = row
    pans = {}
    for band, s in runs.pans.items():
        old = product['runs']['shallow'][f'resort/{band}']
        pans[band] = {t: dict(night_minimum=s['oxygen'][f'{t}/renewal/supply']['nights']['surface_minimum'],
                              hypoxic_nights=s['oxygen'][f'{t}/renewal/supply']['nights']['hypoxic_share'],
                              anoxic_nights=s['oxygen'][f'{t}/renewal/supply']['nights']['anoxic_share'],
                              at_300k=dict(night_minimum=old['oxygen'][f'{t}/renewal/supply']['nights'][
                                  'surface_minimum'],
                                  hypoxic_nights=old['oxygen'][f'{t}/renewal/supply']['nights']['hypoxic_share']))
                      for t in TROPHIC}
    return dict(deep_352m=out, pans_3m=pans)


SHARE_OF_NIGHTS = dict(half=0.5, quarter=0.25, three_quarters=0.75)   # the ventilated depth and its bracket


def floor_crossing(depths, shares, level):
    """The floor depth at which the share of nights mixing a lake to its floor falls through level, interpolated in the
    log of the depth between the ladder's rungs, with its bound where the rungs do not bracket it: 'at least' the
    deepest rung if that still passes level, 'at most' the shallowest if none does."""
    d = np.log(np.asarray(depths, dtype=float))
    s = np.asarray(shares, dtype=float)
    below = np.flatnonzero(s < level)
    if not below.size:
        return float(np.exp(d[-1])), 'at least'
    j = int(below[0])
    if j == 0:
        return float(np.exp(d[0])), 'at most'
    f = (s[j - 1] - level) / (s[j - 1] - s[j])
    return float(np.exp(d[j - 1] + f * (d[j] - d[j - 1]))), None


def thresholds(lad, base=None):
    """The depths that sort the lakes at 301 K. The ventilated depth is the floor depth to which the nights mix a lake
    to its floor on half the nights, with a quarter and three quarters of them beside it, interpolated in the log of
    the depth between the ladder's rungs. The wind's mixing below the night's mixed layer keeps a shallower lake's
    deep water near the mixed layer's density, so the nights reach the floors of lakes far deeper than they mix a
    352 m lake on its ordinary nights; those, the nights that stop short of its floor, are given beside it with every
    night's mixing. base fills the rungs a ladder lacks (the share's bracket runs only the 352 m lakes) from the
    central ladder."""
    out = dict(ventilated_depth_m={}, ordinary_night_mixing_352m_m={}, night_mixing_all_nights_352m_m={},
               polymictic_depth_m={}, monthly_overturn_depth_m={}, by_depth={},
               rungs_from_the_central_ladder=[] if base is None else
               [f'{d:g}' for d in LD.LADDER_M if not lad['0-15'][f'{d:g}']])
    present = [t for t in TROPHIC if t in lad['0-15']['352']]

    def row(band, d, t):
        r = lad[band][f'{d:g}'].get(t)
        return r if r is not None else base[band][f'{d:g}'][t]

    depths = list(LD.LADDER_M)
    names = ('p10', 'p25', 'median', 'p75', 'p90')
    for band in BANDS:
        out['ventilated_depth_m'][band] = {}
        for t in present:
            shares = [row(band, d, t)['nights_mixing_to_the_floor'] for d in depths]
            z = {}
            for key, level in SHARE_OF_NIGHTS.items():
                z[key], bound = floor_crossing(depths, shares, level)
                if bound:
                    z[f'{key}_bound'] = bound
            out['ventilated_depth_m'][band][t] = z
        r352 = {t: row(band, 352.0, t) for t in present}
        out['ordinary_night_mixing_352m_m'][band] = {t: dict(zip(names, r['night_mixing_short_of_the_floor_m'] or
                                                                 r['night_mixing_m'])) for t, r in r352.items()}
        out['night_mixing_all_nights_352m_m'][band] = {t: dict(zip(names, r['night_mixing_m']))
                                                       for t, r in r352.items()}
        poly, monthly = {}, {}
        for t in present:
            ok = [d for d in depths if row(band, d, t)['nights_mixing_to_the_floor'] >= 0.5]
            poly[t] = max(ok) if ok else None
            ok = [d for d in depths if (row(band, d, t)['complete_overturn_share'] or 0) >= 0.5]
            monthly[t] = max(ok) if ok else None
        out['polymictic_depth_m'][band] = poly
        out['monthly_overturn_depth_m'][band] = monthly
        out['by_depth'][band] = {f'{d:g}': {t: dict(
            nights_mixing_to_the_floor=row(band, d, t)['nights_mixing_to_the_floor'],
            complete_overturn_share=row(band, d, t)['complete_overturn_share'],
            episodes_per_decade=row(band, d, t)['episodes_per_decade'],
            floor_ventilated_in_half_the_nights=row(band, d, t)['floor_ventilated_in_half_the_nights'],
            floor_oxic_share_of_time=row(band, d, t)['floor_oxic_share_of_time'],
            deepest_tenth_mean=row(band, d, t)['deepest_tenth']['mean'],
            deepest_tenth_anoxic_share=row(band, d, t)['deepest_tenth']['anoxic_share'])
            for t in present} for d in depths}
    return out


def giants(runs, st):
    out = {}
    for name, entry in runs.index['cases'].items():
        if not name.startswith('giant_'):
            continue
        s = runs.summaries[name]
        trophic = entry['case']['coupled'].split('/')[0]
        m, names = own_column(s, trophic)
        rec = runs.record(name, keys=('mixed_depth', 'deep_o2', 'bottom_t', 't'))
        sat = C.saturated(st, float(np.asarray(rec['t'][0], float).mean()))
        g = AN.giant(rec, s, m, names, sat)
        g['budget'] = s['budget']
        spin = s['settling']['spin_up']
        fine = spin.get('at_the_analysis_step') or {}
        g['spin_up'] = dict(passes=spin['passes'], settled=spin['settled'], limit=spin.get('limit'),
                            at_the_analysis_step=dict(passes=fine.get('passes'), settled=fine.get('settled'),
                                                      limit=fine.get('limit'),
                                                      overturned_by_pass=fine.get('overturned_by_pass')),
                            passes_with_an_overturn=int(sum(spin.get('overturned_by_pass', []))),
                            last_pass_with_an_overturn=max([i for i, o in enumerate(spin.get('overturned_by_pass', []))
                                                            if o], default=None),
                            bottom_age_days_last_two=spin.get('bottom_age_days_by_pass', [])[-2:])
        out[name] = dict(g, band=BANDS[entry['case']['band']], trophic=trophic, depth_m=entry['case']['depth_m'],
                         eta=entry['case'].get('eta', 'central'), deep_mixing=s['deep_mixing'])
        del rec
    return out


def meromictic(runs, st, gas):
    out = {}
    for name, entry in runs.index['cases'].items():
        if not name.startswith('meromictic_'):
            continue
        s = runs.summaries[name]
        case = entry['case']
        trophic = case['coupled'].split('/')[0]
        m, names = own_column(s, trophic)
        rec = runs.record(name, keys=('mixed_depth', 'deep_o2', 'bottom_t', 't', 'oxygen', 'debt', 'solute'))
        bsn = AN.basin_of(rec)
        sat = C.saturated(st, float(np.asarray(rec['t'][0], float).mean()))
        g = AN.giant(rec, s, m, names, sat)
        below = bsn.centre > case['step_depth_m']
        sol = np.asarray(rec['solute'], float)
        step = s['case']['step_solute_g_m3']
        kept = (sol[:, below].mean(1) - sol[:, ~below].mean(1)) / step      # the step's contrast left, by day
        ox = np.asarray(rec['oxygen'][:, :, m], float)
        debt = np.asarray(rec['debt'][:, :, m], float)
        days = np.asarray(rec['step'], float) * float(rec['dt']) / 86400.0
        deep_ox = ox[:, below] @ bsn.volume[below] / bsn.volume[below].sum()
        deep_debt = debt[:, below] @ bsn.volume[below] / bsn.volume[below].sum()
        first = lambda x, level: round(float(days[np.argmax(x < level)] / 365.25), 2) if (x < level).any() else None
        t_deep = float(np.asarray(rec['t'][:, below], float).mean())
        out[name] = dict(band=BANDS[case['band']], trophic=trophic, depth_m=case['depth_m'],
                         step_kg_m3=case['step_kg_m3'], step_depth_m=case['step_depth_m'],
                         step_contrast_left=[round(float(kept[days < 365.25].mean()), 3),
                                             round(float(kept[days > days[-1] - 365.25].mean()), 3)],
                         complete_overturn_episodes=g['complete_overturn_episodes'],
                         night_mixing_m=g['night_mixing_m'], deepest_night_mixing_m=g['deepest_night_mixing_m'],
                         below_step_years_to_hypoxia=first(deep_ox, R.HYPOXIC),
                         below_step_years_to_anoxia=first(deep_ox, R.ANOXIC),
                         below_step_debt_g_m3_by_decade=[round(float(deep_debt[(days >= y * 365.25) &
                                                                               (days < (y + 10) * 365.25)].mean()), 2)
                                                         for y in range(0, 50, 10)],
                         below_step_t_c_first_and_last_year=[
                             round(float(np.asarray(rec['t'][:, below], float)[days < 365.25].mean()), 3),
                             round(float(np.asarray(rec['t'][:, below], float)[days > days[-1] - 365.25].mean()), 3)],
                         years_to_bubble_point_if_isolated=bubble_years(gas, case['depth_m'], t_deep, trophic),
                         ice_max_m=g['ice_max_m'], budget=s['budget'])
        del rec
    return out


def checks(runs, lad):
    """The ladder's sensitivity to the step and the layers at 0-15, and the 352 m lakes' to the wind's share."""
    out = {}

    def row(name, depth):
        s = runs.summaries[name]
        m, names = own_column(s, 'mesotrophic')
        stats, _ = AN.column(runs.record(name), s, m, names)
        return dict(night_mixing_median_m=stats['night_mixing_m'][2],
                    nights_mixing_to_the_floor=stats['nights_mixing_to_the_floor'],
                    episodes_per_decade=stats['episodes_per_decade'],
                    floor_oxic_share_of_time=stats['floor_oxic_share_of_time'],
                    deepest_tenth_mean={t: s['oxygen'][f'{t}/renewal/supply']['deepest_tenth']['mean']
                                        for t in TROPHIC}, bottom_t_c=s['bottom_t_c']['mean'])

    def central(depth):
        s = runs.summaries[LD.ladder_name(0, depth, 'mesotrophic')]
        st = lad['0-15'][f'{depth:g}']['mesotrophic']
        return dict(night_mixing_median_m=st['night_mixing_m'][2],
                    nights_mixing_to_the_floor=st['nights_mixing_to_the_floor'],
                    episodes_per_decade=st['episodes_per_decade'],
                    floor_oxic_share_of_time=st['floor_oxic_share_of_time'],
                    deepest_tenth_mean={t: s['oxygen'][f'{t}/renewal/supply']['deepest_tenth']['mean']
                                        for t in TROPHIC}, bottom_t_c=s['bottom_t_c']['mean'])
    out['hour_steps_352m'] = dict(check=row('check_352m_0-15_mesotrophic_hour_steps', 352.0), central=central(352.0))
    out['two_hour_steps_30m'] = dict(check=row('check_30m_0-15_mesotrophic_two_hour_steps', 30.0),
                                     central=central(30.0))
    out['quarter_metre_layers_12m'] = dict(check=row('check_12m_0-15_mesotrophic_quarter_metre_layers', 12.0),
                                           central=central(12.0))
    out['wind_share_352m'] = {band: {eta: row(LD.ladder_name(b, 352.0, 'mesotrophic', eta), 352.0)
                                     for eta in LD.ETAS} for b, band in enumerate(BANDS)}
    return out


# --- the census -------------------------------------------------------------------------------------------------

def kinds_of(deepest, z_vent):
    """0 self-ventilated, 1 middle, 2 deep, 3 giant."""
    kind = np.full(deepest.shape, 3)
    kind[deepest <= KIND_EDGES_M['deep']] = 2
    kind[deepest <= KIND_EDGES_M['middle']] = 1
    kind[deepest <= z_vent] = 0
    return kind


def tally(mask, area, extra=None):
    out = dict(bodies=int(mask.sum()), area_km2=float(area[mask].sum()))
    if extra:
        for k, v in extra.items():
            out[k] = v
    return out


def census_tables(B, d, b, kind, vent_km2, below_km3, z_by_band):
    """Counts, areas and floors of every kind, by band and size."""
    n = B['n']
    area = B['area_km2']
    band = B['band']
    out = dict(ventilated_depth_m=z_by_band, kinds={}, by_band={}, by_size={}, floor_by_depth_km2={})
    for k, key in enumerate(KINDS):
        m = kind == k
        out['kinds'][key] = dict(name=KIND_NAMES[key], bodies=int(m.sum()), area_km2=float(area[m].sum()),
                                 area_share=float(area[m].sum() / area.sum()),
                                 count_share=float(m.sum() / n), volume_km3=float(B['volume_km3'][m].sum()),
                                 ventilated_floor_km2=float(vent_km2[m].sum()),
                                 ventilated_floor_share=float(vent_km2[m].sum() / max(area[m].sum(), 1e-9)),
                                 ventilated_floor_share_of_each_body=[float(x) for x in np.quantile(
                                     vent_km2[m] / np.maximum(area[m], 1e-12), (0.1, 0.5, 0.9))] if m.any() else None,
                                 volume_below_ventilated_depth_km3=float(below_km3[m].sum()),
                                 mean_depth_m=float(B['volume_km3'][m].sum() * 1e3 / max(area[m].sum(), 1e-9)),
                                 deepest_median_m=float(np.median(B['deepest_m'][m])) if m.any() else None,
                                 closed_bodies=int((m & B['closed']).sum()),
                                 closed_part_area_km2=float(B['closed_area_km2'][m].sum()))
        out['by_band'][key] = {CE.band_name(i): tally(m & (band == i), area) for i in range(4)}
        out['by_size'][key] = {CE.size_name(i): tally(m & (B['size'] == i), area)
                               for i in range(len(CE.SIZE_EDGES_KM2) - 1)}
        cells = m[b['body'] - 1]
        hist = np.histogram(d['depth'][cells], FLOOR_EDGES_M, weights=d['area'][cells])[0] / 1e6
        out['floor_by_depth_km2'][key] = [float(x) for x in hist]
    return out


def census(st, thr, runs, thr_eta=None):
    """Every body sorted by kind under each trophic state's own ventilated depth, at the wind's central share and, for
    the states the ladder ran at it, at its bracket; the central state's sorting at the central share carries the
    needs, the springs and the closed lakes."""
    d, b, B, E, GW, cp, springs = CE.load()
    area = B['area_km2']
    band = B['band']
    zv = {t: np.array([thr['ventilated_depth_m'][BANDS[min(i, 2)]][t]['half'] for i in range(4)]) for t in TROPHIC}
    total = dict(bodies=B['n'], area_km2=float(area.sum()), closed_bodies=int(B['closed'].sum()),
                 closed_area_km2=float(area[B['closed']].sum()), volume_km3=float(B['volume_km3'].sum()))
    out = dict(total=total, central_trophic=CENTRAL_TROPHIC, floor_depth_edges_m=[float(x) for x in FLOOR_EDGES_M],
               by_trophic={})
    keep = {}
    for t in TROPHIC:
        z_body = zv[t][band]
        kind = kinds_of(B['deepest_m'], z_body)
        vent_km2, below_km3 = CE.floor_within(z_body[b['body'] - 1], d, b)
        z_names = {(BANDS[i] if i < 3 else '45-90 (as 30-45)'): float(zv[t][i]) for i in range(4)}
        out['by_trophic'][t] = census_tables(B, d, b, kind, vent_km2, below_km3, z_names)
        keep[t] = dict(kind=kind, z_body=z_body, vent_km2=vent_km2, below_km3=below_km3)
    c = keep[CENTRAL_TROPHIC]
    for key in ('kinds', 'ventilated_depth_m', 'by_band', 'by_size', 'floor_by_depth_km2'):
        out[key] = out['by_trophic'][CENTRAL_TROPHIC][key]
    if thr_eta:
        out['by_trophic_eta'] = {}
        for eta, th in thr_eta.items():
            out['by_trophic_eta'][eta] = {}
            for t in th['ventilated_depth_m']['0-15']:
                z = np.array([th['ventilated_depth_m'][BANDS[min(i, 2)]][t]['half'] for i in range(4)])
                z_body = z[band]
                kind = kinds_of(B['deepest_m'], z_body)
                vent_km2, below_km3 = CE.floor_within(z_body[b['body'] - 1], d, b)
                names = {(BANDS[i] if i < 3 else '45-90 (as 30-45)'): float(z[i]) for i in range(4)}
                out['by_trophic_eta'][eta][t] = census_tables(B, d, b, kind, vent_km2, below_km3, names)
    out['by_share_of_nights'] = {}
    for key in ('quarter', 'three_quarters'):
        zq = np.array([thr['ventilated_depth_m'][BANDS[min(i, 2)]][CENTRAL_TROPHIC][key] for i in range(4)])[band]
        kq = kinds_of(B['deepest_m'], zq)
        out['by_share_of_nights'][key] = dict(ventilated_depth_m=sorted({float(x) for x in zq}),
                                              self_ventilated=tally(kq == 0, area))
    # the census at the ladder's own rungs: the bodies whose deepest point lies within each rung's floor depth, by
    # band, which brackets the self-ventilated kind where the ventilated depth is interpolated between two rungs
    out['within_each_rung_by_band'] = {
        CE.band_name(i): {f'{z:g}': tally((band == i) & (B['deepest_m'] <= z), area) for z in LD.LADDER_M}
        for i in range(4)}
    out['within_each_rung_by_band']['all'] = {f'{z:g}': tally(B['deepest_m'] <= z, area) for z in LD.LADDER_M}
    return out, dict(d=d, b=b, B=B, E=E, GW=GW, cp=cp, springs=springs, **c)


# --- closed lakes -----------------------------------------------------------------------------------------------

def closed_lakes(cz):
    """The closed parts by kind: how long each kind takes along the sea chemistry's left-alone path, and which could
    hold a chemocline once soda."""
    cp = cz['cp']
    has, area = cp['has'], cp['area_km2']
    kind = kinds_of(cp['deepest_m'], cz['z_body'])
    s_star = {f'{f:g}/{m:g}': m / (SALINITY_DENSITY * f) for f in FRESHENING for m in MARGINS_KG_M3}
    out = dict(parts=int(has.sum()), area_km2=float(area[has].sum()), by_kind={})
    for k, key in enumerate(KINDS):
        m = has & (kind == k)
        if not m.any():
            continue
        w = area[m]
        row = dict(parts=int(m.sum()), area_km2=float(w.sum()),
                   mean_depth_m=float((cp['depth_m'][m] * w).sum() / w.sum()))
        for label in ('by_scaling', 'by_balance'):
            f = cp[label][m]
            row[label] = dict(
                years_to_1_g_kg={c: CE.weighted_median(f * cp['milestones'][c]['one_g_kg'], w)
                                 for c in ('high', 'central', 'low')},
                years_to_earth_salinity={c: CE.weighted_median(f * cp['milestones'][c]['earth_salinity'], w)
                                         for c in ('high', 'central', 'low')},
                years_to_hold_a_chemocline={k2: CE.weighted_median(f * CE.years_to(cp['path'], v), w)
                                            for k2, v in s_star.items()})
        out['by_kind'][key] = row
    deep_closed = has & (kind >= 1)
    out['meromixis_candidates'] = dict(
        rule='closed parts whose deepest point lies beyond the ventilated depth: once their salinity passes the '
             'salinity at which a wet spell\'s freshening of the mixed layer outweighs the margin, they can hold a '
             'chemocline through the cold spells',
        bodies=int(deep_closed.sum()), area_km2=float(area[deep_closed].sum()),
        by_kind={key: dict(parts=int((has & (kind == k)).sum()), area_km2=float(area[has & (kind == k)].sum()))
                 for k, key in enumerate(KINDS) if k >= 1})
    return out


# --- springs ----------------------------------------------------------------------------------------------------

def spring_fed(cz, thr, lad, water, gi):
    """Which lakes the groundwater could hold meromictic, at most: their inflow all reaching the deep layer across the
    floor, spread over it, against the time between the overturns that would flush it. The springs' density excess
    is water's one temperature row, which main() builds at the 352 m oligotrophic lakes' mean bottom temperature
    (spring_water(..., temps=(t_deep,)))."""
    B, GW, kind = cz['B'], cz['GW'], cz['kind']
    area = B['area_km2']
    floor = B['area_km2']
    deep_floor = floor - cz['vent_km2']
    v_deep = cz['below_km3']
    rows = {k for v in water.values() for k in v['at']}
    if len(rows) != 1:
        raise ValueError(f'spring_fed takes the springs\' water at one temperature, the deep water\'s; got {rows}')
    t_key = rows.pop()
    excess = min(v['at'][t_key]['density_excess_kg_m3'][0] for v in water.values())
    excess_hi = max(v['at'][t_key]['density_excess_kg_m3'][1] for v in water.values())
    gain_m_yr = np.where(GW['resolved'], np.nan_to_num(GW['gain_m_yr']), np.nan_to_num(GW['by_shore_m_yr']))
    inflow_deep_km3 = gain_m_yr * deep_floor / 1000.0                 # km3 a year into the deep layer (m * km2)
    gain = np.where(v_deep > 0, inflow_deep_km3 / np.maximum(v_deep, 1e-12), 0.0)    # share of the deep layer a year
    # the time the springs have before an overturn flushes the deep layer: the longest stretch without one in the run
    # of the body's band at the next depth at or below its deepest point, the ladder's lakes and then the giants run
    # in that band, for the central trophic state at the wind's central share; beyond the deepest run, that run's
    interval = np.zeros(B['n'])
    band_of = np.minimum(B['band'], 2)
    for b, band in enumerate(BANDS):
        rows = [(d, lad[band][f'{d:g}'][CENTRAL_TROPHIC]['longest_without_overturn_years']) for d in LD.LADDER_M]
        rows += [(g['depth_m'], g['longest_without_overturn_years']) for g in gi.values()
                 if g['band'] == band and g['trophic'] == CENTRAL_TROPHIC and g['eta'] == 'central']
        rows.sort()
        depths = np.array([r[0] for r in rows])
        quiet = np.array([r[1] for r in rows])
        nearest = np.clip(np.searchsorted(depths, B['deepest_m']), 0, depths.size - 1)
        interval = np.where(band_of == b, quiet[nearest], interval)
    out = dict(rule='inflow per area of the 15 km lakes (or, for a body too small for that grid, its band\'s inflow '
                    'per metre of shore along its own shore), spread over the floor, the part below the ventilated '
                    'depth entering the deep layer at the deep water\'s temperature with the springs\' density '
                    'excess; a candidate builds the margin within the longest stretch its kind of lake goes '
                    'without an overturn',
               density_excess_kg_m3=[excess, excess_hi], margins_kg_m3=list(MARGINS_KG_M3), by_margin={})
    for margin in MARGINS_KG_M3:
        years = np.where(gain > 0, margin / np.maximum(gain * excess, 1e-30), np.inf)
        years_hi = np.where(gain > 0, margin / np.maximum(gain * excess_hi, 1e-30), np.inf)
        cand = (kind >= 1) & (years <= interval)
        cand_hi = (kind >= 1) & (years_hi <= interval)
        out['by_margin'][f'{margin:g}'] = dict(
            candidates=tally(cand, area), resolved=tally(cand & GW['resolved'], area),
            by_shore=tally(cand & ~GW['resolved'], area),
            with_silica=tally(cand_hi, area),
            by_kind={key: tally(cand & (kind == k), area) for k, key in enumerate(KINDS)},
            median_years_to_margin_resolved_deep_lakes=float(np.median(years[(kind >= 2) & GW['resolved']])))
        cz[f'spring_candidates_{margin:g}'] = cand
    out['groundwater'] = dict(GW['check'], resolved_bodies=int(GW['resolved'].sum()),
                              resolved_area_km2=float(area[GW['resolved']].sum()),
                              gain_m_yr_resolved=[float(x) for x in np.nanquantile(GW['gain_m_yr'], (0.1, 0.5, 0.9))],
                              gain_m_yr_by_shore=[float(x) for x in np.nanquantile(GW['by_shore_m_yr'],
                                                                                    (0.1, 0.5, 0.9))],
                              unresolved_bodies=int((~GW['resolved']).sum()),
                              unresolved_in_seeping_cells=int(((~GW['resolved']) & (GW['seeping_around_share'] > 0.5))
                                                              .sum()))
    return out


# --- needs ------------------------------------------------------------------------------------------------------

def remedies_301(gas, runs):
    """lake_anoxia's remedies and degassing pipes on the three 352 m lakes at 301 K (mesotrophic physics, its
    convention), labelled as its resort."""
    deep = {}
    for b, band in enumerate(BANDS):
        name = LD.ladder_name(b, R.DEEP_M, 'mesotrophic')
        entry = runs.index['cases'][name]
        f = list(entry['files'])[0]
        case = dict(target='resort', band=band, every=int(np.load(C.CASES / f)['every']))
        deep[name] = A.deep_lake(gas, f, case)
    temps = dict(resort=float(np.mean([v['bottom_t_c'] for v in deep.values()])))
    with mock.patch.object(A, 'TARGET_ORDER', ('resort',)):
        rem = A.remedies(gas, deep, temps)['targets']['resort']
        deg = A.degassing(gas, deep)['resort']
    return deep, temps, rem, deg


def circulation_kw_per_km2(depth):
    """lake_anoxia's circulation: Lorenzen & Fast's air per km2 compressed isothermally to the depth at lunar gravity,
    at its compressor efficiency (kW per km2)."""
    ratio = float(G.local_pressure_atm(depth, A.SURFACE_ATM, MOON_SURFACE_GRAVITY)) / A.SURFACE_ATM
    q = A.CIRCULATION_AIR_M3_MIN_KM2 / 60.0
    return A.SURFACE_ATM * G.ATM * q * math.log(ratio) / A.COMPRESSOR_EFFICIENCY / 1000.0


def permissible_load(mean_depth_m, q_s, settling=1.0):
    """Vollenweider's load (g P per m2 a year) holding a lake at the oligotrophic 10 mg/m3, as lake_anoxia takes it,
    with its settling at Earth's gravity (1) or carried to the Moon's (lake_anoxia's SETTLING); 0 cancels the
    retention, the most anoxic floors can do at steady state, by returning all that settles."""
    return A.OECD_TP['oligotrophic'] * q_s * (1.0 + settling * math.sqrt(mean_depth_m / q_s)) / 1000.0


RETENTION = dict(A.SETTLING, anoxic_floors=0.0)   # lake_anoxia's settling brackets, and the retention at most cancelled


def needs(cen, cz, lad, thr, rem, deg, closed, springs_fed, mero, anoxia):
    kinds = cen['kinds']
    B, kind = cz['B'], cz['kind']
    area = B['area_km2']
    q_s = anoxia['phosphorus']['hydraulic_load_m_yr']
    dissolved = anoxia['phosphorus']['loads_g_m2_yr']['dissolved']
    out = {}
    for key in KINDS:
        z = kinds[key]['mean_depth_m']
        out[key] = dict(bodies=kinds[key]['bodies'], area_km2=kinds[key]['area_km2'],
                        oligotrophic_load_g_p_m2_yr={k: permissible_load(z, q_s, f) for k, f in RETENTION.items()},
                        dissolved_load_over_it={k: [x / permissible_load(z, q_s, f) for x in dissolved]
                                                for k, f in RETENTION.items()})
    # moorings: one at the deepest point of every lake deeper than the ventilated depth and at least SMALL_KM2, and one
    # more for each further MOORING_KM2 of a large lake; the smaller ones take profiles from boats or a sensor string
    deep_enough = kind >= 1
    big = deep_enough & (area >= SMALL_KM2)
    moorings = np.where(big, 1 + np.floor(np.maximum(area - MOORING_KM2, 0.0) / MOORING_KM2), 0)
    out['monitoring'] = dict(
        rule=f'one profiling mooring at the deepest point of every lake deeper than the ventilated depth and of at '
             f'least {SMALL_KM2:g} km2, and one more for each further {MOORING_KM2:g} km2; smaller ones profiled from '
             f'boats or a sensor string; the self-ventilated lakes watched from orbit and the shore',
        earth=EARTH_MOORINGS, moorings=int(moorings.sum()),
        by_kind={key: int(moorings[kind == k].sum()) for k, key in enumerate(KINDS)},
        by_size={CE.size_name(i): int(moorings[B['size'] == i].sum()) for i in range(len(CE.SIZE_EDGES_KM2) - 1)},
        small_deep_lakes_profiled=tally(deep_enough & ~big, area),
        watched_only=tally(kind == 0, area))
    # the deep kind and the giants together, each trophic state's rate per km2 (from the 352 m lakes at 301 K) taken
    # for both kinds
    deep_area = float(area[kind >= 2].sum())
    oxy = rem['traits']
    out['oxygenation_for_deep_water_life'] = dict(
        kw_per_km2={t: oxy[t]['oxygenation']['kw_per_km2'] for t in TROPHIC},
        o2_t_km2_day={t: oxy[t]['oxygenation']['o2_t_km2_day'] for t in TROPHIC},
        deep_and_giant_lakes_area_km2=deep_area,
        rate_taken_for='the deep kind and the giants alike',
        gw_if_chosen={f'{s:g}': {t: [k * deep_area * s / 1e6 for k in oxy[t]['oxygenation']['kw_per_km2']]
                                 for t in TROPHIC} for s in OXYGENATION_SHARES},
        deep_t_c=rem['deep_t_c'], warm_transient_factor=rem['warm_factor'])
    middle_area = float(area[kind == 1].sum())
    out['middle'].update(
        circulation_kw_per_km2={f'{d:g}': circulation_kw_per_km2(d) for d in (20.0, 30.0, 50.0)},
        circulation_gw_if_all=circulation_kw_per_km2(50.0) * middle_area / 1e6,
        held_meromictic={k: dict(years_to_anoxia=v['below_step_years_to_anoxia'],
                                 years_to_bubble_point=v['years_to_bubble_point_if_isolated'],
                                 step_contrast_left=v['step_contrast_left'],
                                 complete_overturn_episodes=v['complete_overturn_episodes']) for k, v in mero.items()})
    cand_area = closed['meromixis_candidates']['area_km2']
    spring = springs_fed['by_margin'][f'{MARGINS_KG_M3[1]:g}']['candidates']
    out['meromictic_candidates'] = dict(
        closed_deep=dict(bodies=closed['meromixis_candidates']['bodies'], area_km2=cand_area),
        spring_fed=spring,
        degassing_per_km2={t: dict(water_m3_s=deg[t]['water_m3_s_per_km2'], pumping_kw=deg[t]['pumping_kw_per_km2'],
                                   methane_kw=deg[t]['methane_energy_kw_per_km2'],
                                   self_lifting=deg[t]['self_lifting']) for t in TROPHIC},
        degassing_gw_if_all={t: deg[t]['pumping_kw_per_km2'] * (cand_area + spring['area_km2']) / 1e6
                             for t in TROPHIC})
    out['built_lakes'] = built_lakes(lad, thr, rem)
    return out


def built_lakes(lad, thr, rem):
    """The design rule for reservoirs and park lakes, which the geography set does not hold."""
    zv = {band: thr['ventilated_depth_m'][band] for band in BANDS}
    rows = {}
    for depth in (5.0, 8.0, 12.0, 16.0, 20.0, 30.0, 50.0):
        rows[f'{depth:g}'] = {band: {t: dict(
            nights_mixing_to_the_floor=lad[band][f'{depth:g}'][t]['nights_mixing_to_the_floor'],
            floor_oxic_share_of_time=lad[band][f'{depth:g}'][t]['floor_oxic_share_of_time'],
            night_surface_minimum_median=(lad[band][f'{depth:g}'][t].get('shallows_nights') or {}).get(
                'surface_minimum', [None, None])[1],
            hypoxic_nights=(lad[band][f'{depth:g}'][t].get('shallows_nights') or {}).get('hypoxic_share'),
            years_to_bubble_point_if_isolated=lad[band][f'{depth:g}'][t].get('years_to_bubble_point_if_isolated'))
            for t in ('oligotrophic', 'new_reservoir')} for band in BANDS}
    new = R.TROPHIC['new_reservoir']
    oxygenation = rem['traits']['new_reservoir']['oxygenation']
    return dict(ventilated_depth_m=zv, ladder=rows,
                first_decade_oxygen_consumed_t_km2_day=new['sediment'] + new['deep_water'],
                first_decade_oxygen_delivered_t_km2_day=oxygenation['o2_t_km2_day'],
                first_decade_kw_per_km2=oxygenation['kw_per_km2'],
                circulation_kw_per_km2={f'{d:g}': circulation_kw_per_km2(d) for d in (20.0, 30.0, 50.0)})


# --- overturn events --------------------------------------------------------------------------------------------

def overturn_events(runs, lad, layers):
    """The nutrient-poor lakes' overturns at 301 K: how often the surface falls below 2 mg/L and for how long, and
    the oxygen the separate shallows keep at the same times."""
    out = {}
    for b, band in enumerate(BANDS):
        rows = {}
        pans = {}
        for depth in (5.0, 8.0):
            name = LD.ladder_name(b, depth, 'oligotrophic')
            rec = runs.record(name, keys=('surface_o2', 'bottom_o2'))
            m, _ = own_column(runs.summaries[name], 'oligotrophic')
            pans[depth] = (np.asarray(rec['surface_o2'][:, m], float), np.asarray(rec['bottom_o2'][:, m], float),
                           float(rec['dt']))
        for depth in LD.LADDER_M:
            if depth < 16.0:
                continue
            st = lad[band][f'{depth:g}']['oligotrophic']
            ev = layers[(band, depth, 'oligotrophic')]['events']
            name = LD.ladder_name(b, depth, 'oligotrophic')
            dt = float(runs.summaries[name]['case']['step_s'])
            refuge = []
            for e in ev:
                t0 = e['day'] * 86400.0
                for pd, (so, bo, pdt) in pans.items():
                    i0 = int(t0 / pdt)
                    i1 = i0 + int(R.AFTER_OVERTURN_DAYS * 86400.0 / pdt)
                    refuge.append(dict(pan_m=pd, surface_min=round(float(so[i0:i1].min()), 3),
                                       bottom_min=round(float(bo[i0:i1].min()), 3)))
            below = [e for e in ev if e['surface_min'] < R.HYPOXIC]
            low = [r for r, e in zip(refuge[0::2], ev) if e['surface_min'] < R.HYPOXIC] + \
                  [r for r, e in zip(refuge[1::2], ev) if e['surface_min'] < R.HYPOXIC]
            rows[f'{depth:g}'] = dict(
                episodes_per_decade=st['episodes_per_decade'], events=len(ev), events_below_2=len(below),
                surface_min=[e['surface_min'] for e in below],
                hours_below_2=[e['hours_below_2'] for e in below],
                hours_below_0_5=[e['hours_below_0_5'] for e in below],
                hours_from_minimum_back_above_2=[e['hours_from_minimum_back_above_2'] for e in below],
                shallows_lowest_while_below_2={f'{pd:g} m': dict(
                    surface=min([r['surface_min'] for r in low if r['pan_m'] == pd], default=None),
                    bottom=min([r['bottom_min'] for r in low if r['pan_m'] == pd], default=None))
                    for pd in (5.0, 8.0)},
                step_s=dt)
        out[band] = rows
    return out


# --- the product ------------------------------------------------------------------------------------------------

def main(argv=None) -> int:
    tic = time.time()
    runs = Runs()
    st = S.Setting301()
    gas = Gases()
    print('runs and settings', f'{time.time() - tic:.0f} s', flush=True)
    lad, layers = ladder(runs, gas)
    thr = thresholds(lad)
    lad_eta = {eta: ladder(runs, gas, eta)[0] for eta in LD.ETAS}
    thr_eta = {eta: thresholds(x, base=lad) for eta, x in lad_eta.items()}
    print('ladder', f'{time.time() - tic:.0f} s', flush=True)
    conv = lake_oxygen_convention(runs)
    gi = giants(runs, st)
    me = meromictic(runs, st, gas)
    ch = checks(runs, lad)
    print('giants, meromictic and checks', f'{time.time() - tic:.0f} s', flush=True)
    cen, cz = census(st, thr, runs, thr_eta)
    print('census', f'{time.time() - tic:.0f} s', flush=True)
    cl = closed_lakes(cz)
    chem = json.loads(SEA_CHEMISTRY.read_text())
    held = next(r for r in chem['runs'] if all(r[k] == v for k, v in R.HELD_FRESH.items()))
    t_deep = float(np.mean([lad[band]['352'][CENTRAL_TROPHIC]['bottom_t_c']['mean'] for band in BANDS]))
    water = CE.spring_water(cz['springs'], float(held['rows'][-1]['lake_salinity_g_kg']),
                            chem['setting']['si_per_eq'], temps=(t_deep,))
    sf = spring_fed(cz, thr, lad, water, gi)
    deep301, temps, rem, deg = remedies_301(gas, runs)
    anoxia = json.loads(A.OUT.read_text())
    nd = needs(cen, cz, lad, thr, rem, deg, cl, sf, me, anoxia)
    ev = overturn_events(runs, lad, layers)
    print('closed lakes, springs, needs, events', f'{time.time() - tic:.0f} s', flush=True)
    inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (R.OUT, R.LAKES, A.OUT, R.TARGETS, SEA_CHEMISTRY,
                                                             SPRINGS, CE.OUT)}
    inputs.update(CE.inputs())
    inputs[str(LD.INDEX.relative_to(ROOT))] = sha256(LD.INDEX)
    code = [Path(__file__), HERE / 'analysis.py', HERE / 'census.py', HERE / 'columns.py', HERE / 'ladder.py',
            HERE / 'meromictic.py', HERE / 'setting.py', Path(L.__file__), Path(R.__file__), Path(A.__file__),
            Path(G.__file__), ROOT / 'research' / 'studies' / 'lake_oxygen' / 'lakes.py']
    out = dict(
        schema=SCHEMA,
        producer=dict(files={str(p.resolve().relative_to(ROOT)): sha256(p) for p in code}, inputs=inputs,
                      constants=dict(MOON_SURFACE_GRAVITY=MOON_SURFACE_GRAVITY, STANDARD_GRAVITY=STANDARD_GRAVITY)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        setting=dict(resort=S.describe(st), reproduction_at_300k=runs.index['reproduction'],
                     kinds=KIND_NAMES, kind_edges_m=KIND_EDGES_M, central_trophic=CENTRAL_TROPHIC,
                     margins_kg_m3=MARGINS_KG_M3, freshening=FRESHENING),
        ladder=lad, thresholds=thr, ladder_eta=lad_eta, thresholds_eta=thr_eta,
        lake_oxygen_convention=conv, giants=gi, meromictic=me, checks=ch,
        census=cen, closed_lakes=cl, spring_water=water, spring_fed=sf,
        anoxia_at_301k=dict(deep_t_c=temps['resort'], remedies=rem, degassing=deg,
                            deep_lakes={k: dict(overturn_events=v['overturn_events'],
                                                anoxic_floor_days_per_year={t: v['traits'][t]['anoxic_floor'][
                                                    'days_per_year'] for t in TROPHIC},
                                                saturation_ratio_max={t: v['traits'][t]['saturation_ratio_max']
                                                                      for t in TROPHIC}) for k, v in deep301.items()}),
        needs=nd, overturn_events=ev, seconds=round(time.time() - tic, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(rnd(out)) + '\n')
    print('wrote', OUT.relative_to(ROOT), f'in {time.time() - tic:.0f} s')
    return 0


EVIDENCE = (
    'Runnable calculations on one lake model (biosphere/ecology/lake.py, as corrected after the independent review '
    'of 18 October) forced from one GCM\'s climate at T21: the resort at 301 K is the lock\'s weather moved by '
    'climate_targets\' share of the way on the ground\'s temperature, the construction that reproduces the waters\' '
    'oxygen study\'s 300 K runs exactly; each lake is a band lake\'s hypsography scaled to its floor, with oxygen '
    'demand from bounded scenarios of Earth\'s lakes, and the stirring below the night\'s mixing set by the share of '
    'the wind\'s energy Earth\'s lakes pass to their stratified water (Wuest, Piepke & Van Senden 2000), 0.3% with '
    'the bracket 0.2-0.4% on the 352 m lakes and a giant, the flux on Wuest & Lorke\'s lake drag (the drag behind the '
    'measured share is open). The census is arithmetic on '
    'geography\'s drainage product (a first estimate from the climate run A28_dim5_moon, which regenerates after the '
    'design case), LOLA at 16 pixels a degree and groundwater stage 2 re-solved from its own heads; the closed lakes\' '
    'soda follows the sea chemistry\'s box model scaled by each lake\'s turnover; the remedies are lake_anoxia\'s '
    'Earth-based sizes at 301 K; the moorings and the reservoirs\' rule are proposed designs. No lunar lake has been '
    'measured.')
READING_RULE = (
    'Oxygen in g/m3 (mg/L), areas in km2 of lake surface, depths in m. ladder is keyed band, floor depth (m) and '
    'trophic state, each with its own physics and its own renewal/supply column, with the wind\'s deep mixing at its '
    'central share, and ladder_eta the 352 m lakes\' oligotrophic and mesotrophic water at its bracket (low 0.2%, '
    'high 0.4%); each case\'s deep_mixing gives the dissipation the wind gave and the shares of the interfaces at '
    'and held back by the cap; night_mixing_m gives the 10th, '
    '25th, 50th, 75th and 90th percentiles of each night\'s deepest mixing; floor shares are of the lake\'s bed in '
    'plan; episodes group consecutive lunar cycles with a complete overturn. thresholds.ventilated_depth_m is the '
    'floor depth to which the nights mix a lake to its floor on half the nights (half), and on a quarter and three '
    'quarters of them, interpolated in the log of the depth between the ladder\'s rungs, with a bound where the rungs '
    'do not bracket it; ordinary_night_mixing_352m_m gives the 352 m lakes\' nights that stop short of the floor and '
    'night_mixing_all_nights_352m_m every night\'s. census.kinds sort every body by its deepest point against its '
    'band\'s ventilated depth for oligotrophic water at the wind\'s central share, and the edges 50 and 1,000 m '
    '(by_share_of_nights at the quarter and three quarters); census.by_trophic and census.by_trophic_eta sort them '
    'under each state\'s '
    'own ventilated depth and at the share\'s bracket; census.within_each_rung_by_band counts the bodies whose deepest '
    'point lies within each of the ladder\'s floor depths, by band, which brackets the self-ventilated kind between '
    'the rungs the ventilated depth is interpolated across; floor_by_depth_km2 gives each kind\'s floor in '
    'lake_anoxia\'s classes (floor_depth_edges_m), to multiply by the methane study\'s fluxes per area of floor. '
    'Closed bodies are those with at least half their cells in depressions that do not spill. Years to salinities '
    'and chemoclines are per closed part, area-weighted medians over a kind; triples run [high, central, low] '
    'weathering where named. spring_water is at the 352 m oligotrophic lakes\' mean bottom temperature, the one '
    'spring_fed reads. giants and meromictic hold the fifty saved years after each one\'s spin-up. Powers in kW per '
    'km2 of lake and GW in total; needs.oxygenation_for_deep_water_life takes each trophic state\'s rate per km2 for '
    'the deep kind and the giants together; built_lakes gives a new reservoir\'s first-decade oxygen as consumed and '
    'as delivered, up to three times it as the devices\' stirring raises the demand. Phosphorus loads in g P per m2 a '
    'year, with lake_anoxia\'s settling brackets and anoxic_floors, the retention cancelled, the most anoxic floors '
    'can do at steady state.')


if __name__ == '__main__':
    sys.exit(main())
