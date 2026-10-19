"""The lakes' and seas' methane at the resort's 301 K, with each lake's own flux and the methane-eating biology sized.

    python -m research.studies.lake_methane.run --checks     # lake_anoxia's saved columns: the reproduction
    python -m research.studies.lake_methane.run --seas       # the seas' columns
    python -m research.studies.lake_methane.run --targets    # lake_anoxia's regenerated lock, resort and hot world
    python -m research.studies.lake_methane.run --work       # one of several processes filling the cache
    python -m research.studies.lake_methane.run              # all of it, then results/lake_methane.json
    python -m research.studies.lake_methane.run --cache-only # the product from the cache, tracing no case

The lake model is not run here. The tracer (tracer.py, through cases.py) carries methane through columns the model
saved: lake_anoxia's 352 m lakes and its 1,925 m lake in the lock, the resort at 300 K and the hot world, the lake
kinds' columns at 301 K (the ladder of floors from 5 to 352 m, the giants to 3,500 m and the lakes held meromictic)
and the seas' columns. It writes results/lake_methane.json (schema terluna.research.lake-methane/1).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_anoxia import gases as G
from research.studies.lake_anoxia import run as A
from research.studies.lake_methane import biology as B
from research.studies.lake_methane import cases as CS
from research.studies.lake_methane import columns as C
from research.studies.lake_methane import tracer as T
from research.studies.lake_oxygen import run as R

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'research' / 'runs' / 'lake_methane'
OUT = HERE / 'results' / 'lake_methane.json'
SCHEMA = 'terluna.research.lake-methane/1'
ANOXIA = ROOT / 'research' / 'studies' / 'lake_anoxia' / 'results' / 'lake_anoxia.json'
LAKES = ROOT / 'research' / 'studies' / 'lake_oxygen' / 'results' / 'lakes.json'
SEAS_JSON = ROOT / 'research' / 'studies' / 'sea_oxygen' / 'results' / 'seas.json'
CHEMISTRY = ROOT / 'research' / 'studies' / 'sea_chemistry' / 'results' / 'sea_chemistry.json'
SEAS_PRODUCT = ROOT / 'research' / 'studies' / 'sea_oxygen' / 'results' / 'sea_oxygen.json'
CENSUS = ROOT / 'research' / 'studies' / 'lake_kinds' / 'results' / 'lake_kinds.json'
LEDGERS, GAS_FORCING, TARGETS = A.LEDGERS, A.GAS_FORCING, A.TARGETS
TRAITS = A.TRAITS
TG = 1e12
M_CH4 = 16.04
CASES = CS.CASES


def sha256(path) -> str:
    return C.sha256(path)


def repo_path(path) -> str:
    """A path as the repository names it (relative to its root), for the product."""
    p = Path(path)
    return str(p.relative_to(ROOT)) if p.is_absolute() and str(p).startswith(str(ROOT)) else str(p)


def dumps(obj) -> str:
    return A.dumps(obj)


def rnd(x, sig=4):
    return A.rnd(x, sig)


trace_case = CS.trace_case
summarise = CS.summarise

# --- the lane sets ------------------------------------------------------------------------------------------------

COMMUNITIES = (('earth', {}), ('designed', B.FAST), ('poor', B.POOR), ('none', dict(aerobic=False)))


def community_lanes(columns, shares=True, aom=None, extra=()):
    """For each oxygen column: Earth's lake community, the designed and poor ends and no methanotrophs; the methane
    share's low and high ends; the anaerobic pathways' combined sets at their acceptors' supply; and extra sets."""
    lanes = []
    for j in columns:
        for label, kw in COMMUNITIES:
            lanes.append(T.Lane(j, label=label, **dict(B.EARTH, **kw)))
        if shares:
            for s_ in ('low', 'high'):
                lanes.append(T.Lane(j, share=B.SHARES[s_], label=f'share/{s_}', **B.EARTH))
        if aom is not None:
            lanes.extend(l for l in B.aom_lanes(j, aom) if l.label in ('aom/all/central', 'aom/all/high'))
        for make in extra:
            lanes.extend(make(j))
    return lanes


def kinetics_lanes(j):
    """(a) each kinetic parameter at the ends of Earth's ranges, the others at Earth's central values."""
    return [T.Lane(j, label=f'kinetics/{k}', **dict(B.EARTH, **kw)) for k, kw in B.KINETICS.items()
            if k not in ('central', 'designed', 'poor', 'none')]


def overturn_lanes(j):
    """(b) the stock the water keeps against the community's response time."""
    out = []
    for cells in B.STOCKS_CELLS_ML:
        for mu in B.GROWTH_PER_DAY:
            kw = dict(B.EARTH, seed_c_um=cells * B.UM_C_PER_CELL_ML, mu20_per_h=mu / 24.0)
            out.append(T.Lane(j, label=f'overturn/{cells:g}/{mu:g}', **kw))
    return out


def aom_detail(supply):
    """(c) each anaerobic pathway alone at its acceptor's supply, and together."""
    return lambda j: [l for l in B.aom_lanes(j, supply) if l.label not in ('aom/all/central', 'aom/all/high')]


# --- lake_anoxia's regenerated columns: the three targets beside its two rules ------------------------------------

def anoxia_targets(st, product, only=None):
    """lake_anoxia's 352 m lakes and its 1,925 m lake, regenerated on the final lake model, traced: the central wind
    share with the communities and the share's ends, the share's brackets with Earth's community and none, the 1,925
    m lake through its fifty years."""
    meta = C.anoxia_meta(C.ANOXIA_RUNS)
    if C.model_version(meta['producer']['files'].get('biosphere/ecology/lake.py')) != 'current':
        raise SystemExit("lake_anoxia's saved columns are not on the current lake model")
    profiles = C.ANOXIA_RUNS / 'profiles.json'
    out = dict(profiles=repo_path(profiles), profiles_sha256=sha256(profiles),
               deep={}, deepest={})
    four = range(len(TRAITS))
    for name, case in meta['cases'].items():
        if only is not None and not only(name):
            continue
        col = C.anoxia_column(st, name, meta, folder=C.ANOXIA_RUNS)
        if case['kind'] == 'deep':
            central = case.get('eta', 'central') == 'central'
            lanes = community_lanes(four, shares=central) if central else \
                [T.Lane(j, label=k, **dict(B.EARTH, **kw)) for j in four for k, kw in COMMUNITIES
                 if k in ('earth', 'none')]
            res = trace_case(col, lanes, f'target_{name}', cyclic=True)
            out['deep'][name] = dict(target=case['target'], band=case['band'], eta=case.get('eta', 'central'),
                                     case=res['case'], lanes=res['lanes'], passes=res['passes'],
                                     settled=res['settled'])
        else:
            lanes = [T.Lane(j, label=k, **dict(B.EARTH, **kw)) for j in four for k, kw in COMMUNITIES]
            res = trace_case(col, lanes, f'target_{name}', cyclic=False, start_empty=True,
                             chunks=case.get('passes', 5))
            out['deepest'][name] = dict(target=case['target'], eta=case.get('eta', 'central'),
                                        decades=[dict(lanes=p['lanes'], overturns=len(p['case']['overturns']))
                                                 for p in res['parts']], case=res['case'])
        del col
    if only is None:
        out['beside_lake_anoxia'] = two_rule_table(out['deep'], product)
    return out


def two_rule_table(deep, product):
    """The tracer's emission beside lake_anoxia's two release rules, by target and trophic state, the mean of the
    three bands at the central wind share: its surface flux (between overturns and at them) against lake_anoxia's
    pelagic flux under the net-rise and the stock rules, and the bubbles from each lake's own floors against its
    bubbles from the Moon's floors."""
    out = {}
    for target in A.TARGET_ORDER:
        runs = [v for v in deep.values() if v['target'] == target and v['eta'] == 'central']
        rows = {}
        for t in TRAITS:
            def pick(label):
                return [next(l for l in r['lanes'] if l['label'] == label and l['column'].startswith(t + '/'))
                        for r in runs]
            entry = {}
            for label in ('earth', 'designed', 'poor', 'none', 'share/low', 'share/high'):
                ls = pick(label)
                entry[label] = dict(
                    surface_mol_m2_yr=float(np.mean([l['diffusive_mol_m2_yr'] for l in ls])),
                    at_overturns_mol_m2_yr=float(np.mean([l['emitted_at_overturns_mol_m2_yr'] for l in ls])),
                    bubbles_own_floors_mol_m2_yr=float(np.mean([l['bubbles_to_air_mol_m2_yr'] for l in ls])),
                    made_mol_m2_yr=float(np.mean([l['made_anoxic_mol_m2_yr'] for l in ls])),
                    oxycline_oxidised_share=[l['oxycline_oxidised_share'] for l in ls])
            m = product['methane'][target]['traits'][t]
            rule = m['by_release_rule']
            entry['lake_anoxia'] = dict(
                made_mol_m2_yr=m['anoxic_methane_made_mol_m2_yr']['central'],
                bubbles_moon_floors_mol_m2_yr=m['ebullition_mol_m2_yr']['central'],
                net_rise_pelagic_mol_m2_yr=rule['net_rise']['pelagic_mol_m2_yr']['central'],
                stock_pelagic_mol_m2_yr=rule['stock']['pelagic_mol_m2_yr']['central'],
                net_rise_released_mol_m2_yr=rule['net_rise']['released_at_overturns_mol_m2_yr'],
                stock_released_mol_m2_yr=rule['stock']['released_at_overturns_mol_m2_yr'])
            rows[t] = entry
        out[target] = rows
    return out


# --- the lake kinds' columns at 301 K -----------------------------------------------------------------------------

LADDER_DETAIL_M = (30.0, 100.0, 352.0)        # rungs with the overturn grid and each anaerobic pathway
KINETICS_M = (352.0,)                         # rungs with each kinetic parameter at its ends
DETAIL_BAND = 0                               # the biology's sweeps on the 0-15 band, 62% of the lake area
GIANT_PASSES = 10                             # a giant's deep store turns over slowly: up to a century of its decade


def kinds_cases(st301, index, supply, only=None):
    """Every case of the kinds' 301 K index traced: the ladder's rungs; the giants, spun up until settled and saved
    for five passes of the record, from their last saved decade passed over until the methane settles; the lakes held
    meromictic, which start oxygenated with their chemocline and run fifty years, by decades from that start; and the
    wind share's brackets and the checks with Earth's community and none."""
    out = dict(index=repo_path(C.KINDS_RUNS / 'ladder.json'), index_sha256=sha256(C.KINDS_RUNS / 'ladder.json'),
               cases={})
    gs = gas_setting(st301)
    for name, entry in index['cases'].items():
        if only is not None and not only(name):
            continue
        case = entry['case']
        depth = float(case['depth_m'])
        central = case.get('eta', 'central') == 'central' and not name.startswith('check')
        giant = name.startswith('giant')
        col = C.kinds_column(st301, name, last_pass(index, name) if giant else index)
        if name.startswith('ladder') and central:
            extra = []
            if case['band'] == DETAIL_BAND and depth in KINETICS_M:
                extra.append(kinetics_lanes)
            if case['band'] == DETAIL_BAND and depth in LADDER_DETAIL_M:
                extra.extend([overturn_lanes, aom_detail(supply)])
            lanes = community_lanes([0], aom=supply, extra=extra)
            res = trace_case(col, lanes, f'kinds_{name}', cyclic=True)
        elif giant and central:
            detail = [aom_detail(supply)] if depth == 1925.0 and case['band'] == DETAIL_BAND else []
            lanes = community_lanes([0], aom=supply, extra=detail)
            res = trace_case(col, lanes, f'kinds_{name}', cyclic=True, max_passes=GIANT_PASSES)
        elif name.startswith('meromictic'):
            lanes = community_lanes([0], aom=supply, extra=[aom_detail(supply)])
            res = trace_case(col, lanes, f'kinds_{name}', cyclic=False, start_empty=True, chunks=case['passes'])
        else:
            lanes = [T.Lane(0, label=k, **dict(B.EARTH, **kw)) for k, kw in COMMUNITIES if k in ('earth', 'none')]
            res = trace_case(col, lanes, f'kinds_{name}', cyclic=True, max_passes=GIANT_PASSES if giant else 4)
        if name.startswith(('giant', 'meromictic')):
            res = dict(res, floors=floor_table(gs, col, res, lanes))
        out['cases'][name] = dict(case=case, summary=entry['summary'], files=entry['files'], result=res)
        del col
    return out


def floor_table(gs, col, res, lanes):
    """The floor water's gas saturation ratio (1 is the bubble point) from each lane's own methane at its traced
    record's close (a giant's settled cycle) or each decade's (a meromictic lake from its oxygenated start)."""
    parts = res.get('parts') or [res]
    out = {}
    for j, l in enumerate(lanes):
        if l.kind != 'methane' or str(j) not in parts[-1].get('closing_mol_m3', {}):
            continue
        rows = []
        for part in parts:
            ch4 = part['closing_mol_m3'][str(j)][-1]
            ratio, p_gas, p_local = floor_saturation(gs, col, part['closing_snapshot'], l.column, ch4, l.share)
            rows.append(dict(ratio=ratio, gas_atm=p_gas, floor_atm=p_local, ch4_mol_m3=ch4))
        out[l.label] = rows
    return out


def last_pass(index, name):
    """The index with one case trimmed to its last saved pass (a giant's settled last decade)."""
    entry = index['cases'][name]
    last = sorted(entry['files'])[-1]
    trimmed = dict(index)
    trimmed['cases'] = dict(index['cases'])
    trimmed['cases'][name] = dict(entry, files={last: entry['files'][last]})
    return trimmed


def lane(res_lanes, label):
    return next((l for l in res_lanes if l['label'] == label), None)


def long_table(kinds, prefix):
    """The giants' or the meromictic lakes' decades from their oxygenated start, per lane."""
    out = {}
    for name, c in kinds['cases'].items():
        if not name.startswith(prefix):
            continue
        r = c['result']
        decades = []
        for p in r['parts']:
            decades.append(dict(overturns=len(p['case']['overturns']),
                                lanes={l['label']: dict(emitted=l['emitted_mol_m2_yr'], made=l['made_anoxic_mol_m2_yr'],
                                                        at_overturns=l['emitted_at_overturns_mol_m2_yr'],
                                                        bubbles=l['bubbles_to_air_mol_m2_yr'],
                                                        oxidised=l['oxidised_mol_m2_yr'],
                                                        aom=l['aom_floor_mol_m2_yr'] + l['aom_water_mol_m2_yr'],
                                                        inventory_end=l['inventory_mol_m2']['end'])
                                       for l in p['lanes']}))
        out[name] = dict(case={k: v for k, v in c['case'].items() if k != 'start'}, decades=decades,
                         area_km2=r['case']['area_km2'], mean_depth_m=r['case']['mean_depth_m'])
    return out


# --- the Moon's lakes by rung ---------------------------------------------------------------------------------------

LADDER_M = (5.0, 8.0, 12.0, 16.0, 20.0, 30.0, 50.0, 100.0, 352.0)
GIANTS_M = (1000.0, 1925.0, 3500.0)
BAND_EDGES = ((0.0, 15.0), (15.0, 30.0), (30.0, 90.0))


def rung_edges(rungs):
    r = np.asarray(rungs, float)
    return np.concatenate([[0.0], np.sqrt(r[:-1] * r[1:]), [np.inf]])


def lake_bodies():
    """Geography's rain-fed lakes on LOLA as lake_oxygen's lakes reads them: each body's area, deepest point, centre."""
    from research.studies.lake_oxygen import lakes as K
    d = K.lake_depths()
    b = K.bodies(d)
    return dict(area=b['area'][1:], deepest=b['deepest'][1:], lat=b['lat'][1:])


def rung_table(bodies, rungs):
    """Each rung's lakes: those whose deepest point lies nearest it (between the geometric midpoints), their area and
    its shares by latitude band (the 30-45 band's lake standing for 30-90)."""
    edges = rung_edges(rungs)
    total = bodies['area'].sum()
    out = {}
    for i, r in enumerate(rungs):
        m = (bodies['deepest'] >= edges[i]) & (bodies['deepest'] < edges[i + 1])
        a = bodies['area'][m]
        lat = np.abs(bodies['lat'][m])
        bands = [float(a[(lat >= lo) & (lat < hi)].sum() / max(a.sum(), 1.0)) for lo, hi in BAND_EDGES]
        out[f'{r:g}'] = dict(deepest_m=[float(edges[i]), float(edges[i + 1])], bodies=int(m.sum()),
                             area_km2=float(a.sum() / 1e6), area_share=float(a.sum() / total), band_shares=bands)
    return out


def flux_by_rung(ladder, rungs_table, labels):
    """Each rung's flux per m2 of lake surface, its bands weighted by its lakes' area in each band (a band not run
    takes the 0-15 band's), for each trophic state and lane; the giants are rungs too."""
    out = {}
    for rung, info in rungs_table.items():
        w = info['band_shares']
        entry = {}
        for t in TRAITS:
            runs = {C.BANDS.index(b): v['lanes'] for b, v in ladder.get(rung, {}).get(t, {}).items()}
            if not runs:
                continue
            per = {}
            for label in labels:
                vals = {i: v[label]['emitted'] for i, v in runs.items() if label in v}
                if not vals:
                    continue
                full = [vals.get(i, vals.get(0, next(iter(vals.values())))) for i in range(3)]
                per[label] = float(sum(wi * fi for wi, fi in zip(w, full)))
            entry[t] = dict(bands_run=sorted(C.BANDS[i] for i in runs), emitted_mol_m2_yr=per)
        out[rung] = entry
    return out


def moon_lakes(by_rung, rungs_table, area_km2, trait_fallback=True):
    """The Moon's lakes, Tg a year: each rung's flux over its lakes' area; a rung lacking a trophic state (the
    giants' new reservoirs) takes the deepest rung run in that state, and says so."""
    out = {}
    for t in TRAITS:
        per_label, borrowed = {}, []
        for rung, info in rungs_table.items():
            entry = by_rung.get(rung, {}).get(t)
            src = rung
            if entry is None and trait_fallback:
                deeper = [r for r in by_rung if t in by_rung[r]]
                src = max(deeper, key=float) if deeper else None
                entry = by_rung[src][t] if src else None
                borrowed.append(dict(rung=rung, takes=src))
            if entry is None:
                continue
            for label, f in entry['emitted_mol_m2_yr'].items():
                per_label[label] = per_label.get(label, 0.0) + f * M_CH4 * info['area_km2'] * 1e6 / TG
        out[t] = dict(tg_yr=per_label, borrowed=borrowed)
    return out


def uniform_scale_up(by_rung, area_km2):
    """lake_anoxia's way: every lake at the 352 m lakes' flux."""
    return {t: {label: f * M_CH4 * area_km2 * 1e6 / TG for label, f in by_rung['352'][t]['emitted_mol_m2_yr'].items()}
            for t in TRAITS if t in by_rung.get('352', {})}


CENSUS_BANDS = ('0-15', '15-30', '30-45')      # the census reads its 45-90 band as its 30-45, as BAND_EDGES' last does
SPLITS = ('at_100m', 'census', 'at_352m')


def ventilation_split(bodies, by_rung, rungs_table, area_km2, census=None):
    """The Moon's lakes, Tg a year, with the lakes whose deepest point lies between the 100 m and the 352 m rungs split
    at the ventilated depth (the floor the nights reach on half of them): those within it take the 100 m rung's
    band-weighted flux, those beyond it the 352 m rung's, and every other lake its nearest rung's as in moon_lakes.
    The split is set at the lake kinds census's ventilated depth for each trophic state and band (census), and at the
    two rungs for its bracket (the census's self-ventilated area at its central trophic state, 226,500-629,700 km2
    against 403,700). Each split's self-ventilated area, every lake shallower than its ventilated depth, is given beside
    the census's."""
    census = census or json.loads(CENSUS.read_text())['census']
    rungs = [float(r) for r in rungs_table]
    edges = rung_edges(rungs)
    i, j = rungs.index(100.0), rungs.index(352.0)
    deep, a = bodies['deepest'], bodies['area']
    band = np.digitize(np.abs(bodies['lat']), [15.0, 30.0])              # 0-15, 15-30, 30-90
    inside = (deep >= edges[i]) & (deep < edges[j + 1])
    nearest = moon_lakes(by_rung, rungs_table, area_km2)
    out = dict(deepest_m=[float(edges[i]), float(edges[j + 1])], ventilated_depth_m={}, self_ventilated_km2={},
               census_self_ventilated_km2=dict(at_100m=census['within_each_rung_by_band']['all']['100']['area_km2'],
                                               census=census['kinds']['self_ventilated']['area_km2'],
                                               at_352m=census['within_each_rung_by_band']['all']['352']['area_km2'],
                                               trophic=census['central_trophic']),
               moon_tg_yr={s_: {} for s_ in SPLITS})
    for t in TRAITS:
        vent = census['by_trophic'][t]['ventilated_depth_m']
        depths = dict(at_100m=[100.0] * 3, census=[float(vent[b]) for b in CENSUS_BANDS], at_352m=[352.0] * 3)
        out['ventilated_depth_m'][t] = depths
        f100, f352 = (by_rung[r][t]['emitted_mol_m2_yr'] for r in ('100', '352'))
        a100, a352 = (rungs_table[r]['area_km2'] * 1e6 for r in ('100', '352'))
        for split, d in depths.items():
            dv = np.asarray(d)[band]
            within, beyond = float(a[inside & (deep < dv)].sum()), float(a[inside & (deep >= dv)].sum())
            out['self_ventilated_km2'].setdefault(split, {})[t] = float(a[deep < dv].sum() / 1e6)
            tg = {}
            for label, total in nearest[t]['tg_yr'].items():
                if label not in f100 or label not in f352:
                    continue
                swap = (f100[label] * (within - a100) + f352[label] * (beyond - a352)) * M_CH4 / TG
                tg[label] = total + swap
            out['moon_tg_yr'][split][t] = dict(tg_yr=tg, borrowed=nearest[t]['borrowed'],
                                               within_km2=within / 1e6, beyond_km2=beyond / 1e6)
    out['nearest_rung_tg_yr'] = nearest
    return out


def meromictic_flux(mero):
    """A meromictic lake's flux per km2 (t a year), each decade from its oxygenated start and the last, per lane."""
    out = {}
    for name, m in mero.items():
        last = m['decades'][-1]['lanes']
        out[name] = dict(depth_m=m['case']['depth_m'], step_kg_m3=m['case'].get('step_kg_m3'),
                         trophic=m['case']['coupled'].split('/')[0],
                         t_km2_yr={label: v['emitted'] * M_CH4 for label, v in last.items()},
                         by_decade_t_km2_yr={label: [d['lanes'][label]['emitted'] * M_CH4 for d in m['decades']]
                                             for label in last},
                         made_t_km2_yr={label: v['made'] * M_CH4 for label, v in last.items()})
    return out


AOM_SETS = {
    'aom/all/central': 'every pathway at the high-erosion end: the floors\' iron and manganese at the earliest erosion '
                       '(8.3 Gt a year, the soil loss the phosphorus particles imply) with highland soils, nitrate at '
                       'Lake Constance\'s geometric-mean actual rate, sulphate at its renewal',
    'aom/all/high': 'every pathway at the high-erosion end with mare soils and nitrate at Lake Constance\'s highest '
                    'actual rate',
    'aom/iron/latest_highland': 'the floors\' iron and manganese at the latest erosion, 0.71 Gt a year (the '
                                'conservation study\'s 10 m a million years), highland soils',
    'aom/iron/earliest_highland': 'the floors\' iron and manganese at the earliest erosion, highland soils',
    'aom/iron/earliest_mare': 'the floors\' iron and manganese at the earliest erosion, mare soils',
    'aom/nitrate/low': 'nitrate and nitrite at Lake Constance\'s lowest actual rate',
    'aom/nitrate/high': 'nitrate and nitrite at its highest actual rate',
    'aom/sulphate': 'sulphate at 1.4 uM, renewed by the throughflow'}

AOM_LABELS = ('aom/all/central', 'aom/all/high', 'aom/iron/latest_highland', 'aom/iron/earliest_highland',
              'aom/iron/earliest_mare', 'aom/nitrate/low', 'aom/nitrate/high', 'aom/sulphate')


def biology_tables(ladder, kinds=None):
    """(a) the 352 m lakes' emission under each kinetic variant; (b) under each stock and response time, and at each
    overturn the shares of the stored methane oxidised and emitted and the oxygen the lake holds for it; (c) each
    anaerobic pathway's cut on every rung that carries it. Each lane is set beside Earth's community on the same bands:
    the sweeps run on the 0-15 band, the combined anaerobic sets on every band."""
    band0 = C.BANDS[DETAIL_BAND]

    def pair(rung, t, label, key='emitted', only=None):
        """The lane's flux and Earth's community's, each the mean over the bands that carry the lane (or the one)."""
        runs = ladder.get(rung, {}).get(t, {})
        have = [v['lanes'] for b, v in runs.items() if label in v['lanes'] and (only is None or b == only)]
        if not have:
            return None, None
        return float(np.mean([v[label][key] for v in have])), float(np.mean([v['earth'][key] for v in have]))

    out = dict(band=band0, kinetics={}, overturns={}, anaerobic={})
    for t in TRAITS:
        row = {k: pair('352', t, f'kinetics/{k}', only=band0)[0] for k in B.KINETICS
               if k not in ('central', 'designed', 'poor', 'none')}
        row.update({k: pair('352', t, k, only=band0)[0] for k in ('earth', 'designed', 'poor', 'none')})
        out['kinetics'][t] = row
        grid = {}
        for rung in ('30', '100', '352'):
            grid[rung] = dict(earth=pair(rung, t, 'earth', only=band0)[0],
                              grid={f'{cells:g}': {f'{mu:g}': pair(rung, t, f'overturn/{cells:g}/{mu:g}',
                                                                   only=band0)[0]
                                                   for mu in B.GROWTH_PER_DAY} for cells in B.STOCKS_CELLS_ML})
        out['overturns'][t] = grid
        cuts = {}
        for rung in ladder:
            row = {}
            for label in AOM_LABELS:
                v, e = pair(rung, t, label)
                if v is not None:
                    row[label] = dict(emitted=v, earth=e, cut_share=(1.0 - v / e) if e > 0 else None)
            if row:
                cuts[rung] = row
        out['anaerobic'][t] = cuts
    if kinds is not None:
        out['overturn_events'] = overturn_table(kinds)
    return out


def giants_reading(ladder, forms, rungs_table, census=None):
    """The giants' deep oxygen as the lake kinds study saved it (the deepest tenth's mean, mg/L, by decade of its
    record), what that makes of the traced methane, and for the 3,500 m nutrient-poor giant, which the kinds study
    finds anoxic within its record's next 13 years, the Moon-wide methane it would add if it settled anoxic as the
    1,925 m nutrient-poor giant does under Shih and Bouffard's weaker deep mixing (their change per m2, Earth's
    community, over the 3,500 m rung's lakes)."""
    census = census or json.loads(CENSUS.read_text())
    out = {}
    for name, g in census['giants'].items():
        o2 = g['deepest_tenth_by_decade']
        anoxic = max(o2) < 0.05
        change = o2[-1] - o2[0]
        if anoxic:
            reading, methane = 'anoxic deep water throughout', 'as traced, its store still shrinking slightly'
        elif g.get('years_to_anoxia'):
            reading, methane = f"losing oxygen, anoxic within {g['years_to_anoxia']:g} years", 'a lower value'
        elif change > 0.01:
            reading, methane = 'gaining oxygen slowly', 'as traced or slightly an upper value'
        else:
            reading, methane = 'settled', 'as traced'
        out[name] = dict(deepest_tenth_o2_by_decade_mg_l=o2, years_to_anoxia=g.get('years_to_anoxia'),
                         reading=reading, methane=methane)
    shih = forms['cases'].get('shih_301k_giant_1925m_0-15_oligotrophic')
    osb = ladder.get('1925', {}).get('oligotrophic', {}).get('0-15')
    if shih is not None and osb is not None:
        d = next(l for l in shih['result']['lanes'] if l['label'] == 'earth')['emitted_mol_m2_yr'] - \
            osb['lanes']['earth']['emitted']
        out['if_the_3500m_nutrient_poor_giant_settled_anoxic'] = dict(
            change_mol_m2_yr=d, moon_tg_yr=d * M_CH4 * rungs_table['3500']['area_km2'] * 1e6 / TG)
    return out


def overturn_table(kinds, rungs=(30.0, 100.0, 352.0)):
    """(b) the overturns of the rungs that carry the overturn grid (the 0-15 band, the central wind share): how often
    they come, the methane each lane has stored when one begins and the oxygen the overturned lake then holds against
    the 64 g a mole that store needs, and what leaves within the overturns' windows (their first full mixing to thirty
    days past their last) a year: emitted, across the surface and by bubbles, and oxidised in the water."""
    out = {}
    for name, c in kinds['cases'].items():
        case = c['case']
        if not name.startswith('ladder') or case.get('eta', 'central') != 'central' or \
                case['band'] != DETAIL_BAND or float(case['depth_m']) not in rungs:
            continue
        r = c['result']
        ov = r['case']['overturns']
        held = [e['o2_held_g_m2'][0] for e in ov]               # the case's one traced column
        rows = {}
        for l in r['lanes']:
            stored = [e['stored_mol_m2'] for e in l['events']]
            rows[l['label']] = dict(
                stored_per_overturn_mol_m2=float(np.mean(stored)) if stored else 0.0,
                o2_needed_over_held=float(np.sum(stored) * 64.0 / np.sum(held)) if held and np.sum(held) > 0 else None,
                emitted_at_overturns_mol_m2_yr=l['emitted_at_overturns_mol_m2_yr'],
                oxidised_at_overturns_mol_m2_yr=l['oxidised_at_overturns_mol_m2_yr'],
                emitted_mol_m2_yr=l['emitted_mol_m2_yr'])
        out.setdefault(f"{float(case['depth_m']):g}", {})[case['coupled'].split('/')[0]] = dict(
            overturns=len(ov), years=r['case']['years'], overturns_per_year=len(ov) / r['case']['years'],
            o2_held_per_overturn_g_m2=float(np.mean(held)) if held else None, lanes=rows)
    return out


# --- the reproduction and the seas' part --------------------------------------------------------------------------

def anoxia_reproduction(st, product, folder=None):
    """Where the tracer's physics coincides with lake_anoxia's: its methane made, read from the daily columns as
    lake_anoxia reads them (the saved oxygen, its rule for the water counted); the debt's clones against the model's
    own debt, on the mixed layer's oxygen step by step; the water's age; the overturns and the oxygen they leave."""
    folder = folder or C.ANOXIA_RUNS
    meta = C.anoxia_meta(folder)
    rule = T.lake_anoxia_made_below(product)
    four = range(len(TRAITS))
    profiles = Path(folder) / 'profiles.json'
    out = dict(profiles=repo_path(profiles), profiles_sha256=sha256(profiles),
               product_schema=product['schema'], made_below_g_m3=rule, lakes={})
    for name, case in meta['cases'].items():
        if case['kind'] != 'deep' or case.get('eta', 'central') != 'central':
            continue
        col = C.anoxia_column(st, name, meta, folder=folder)
        lanes = [T.Lane(0, kind='age')] + [T.Lane(j, label='earth', **B.EARTH) for j in four]
        rep = trace_case(col, lanes, f'repro_{name}', cyclic=False, ml_oxygen=False, clones=True, made_below=rule,
                         oxygen_budget=False)
        clones = trace_case(col, [T.Lane(j, kind='debt') for j in four], f'clones_{name}', cyclic=False,
                            clones=True, made_below=rule, oxygen_budget=False)
        theirs = product['deep_lakes'][name]
        rows = {}
        for t, row in zip(TRAITS, rep['lanes']):
            them = theirs['traits'][t]
            debt = clones['clones'][f'debt/{t}/{C.CENTRAL}']
            rows[t] = dict(made_mol_m2_yr=row['made_anoxic_mol_m2_yr'],
                           lake_anoxia_made_mol_m2_yr=them['methane_made_mol_m2_yr']['central'],
                           debt_inventory_mean_mol_m2=debt['clone_mean_mol_m2'],
                           model_inventory_mean_mol_m2=debt['model_mean_mol_m2'],
                           lake_anoxia_inventory_mean_mol_m2=them['inventory_mean_mol_m2'],
                           o2_held_after_overturns_g_m2=[e['o2_held_g_m2'][TRAITS.index(t)]
                                                         for e in rep['case']['overturns']],
                           lake_anoxia_o2_after_overturns_g_m2=them['oxygen_after_overturns_g_m2'])
        out['lakes'][name] = dict(target=case['target'], band=case['band'], overturns=len(rep['case']['overturns']),
                                  lake_anoxia_overturns=theirs['overturn_events'], age=rep['clones']['age'],
                                  transfer_u_star_rel_err=u_star_error(col), traits=rows)
        del col
    return out


def u_star_error(col):
    """The rebuilt forcing's friction velocity against the one the run saved, the largest relative error."""
    if col.u_star is None or col.u_star_rebuilt is None:
        return None
    u = np.asarray(col.u_star, float)
    return float(np.max(np.abs(np.asarray(col.u_star_rebuilt) - u) / np.maximum(u, 1e-9)))


def seas_part(ladder=None):
    """The seas' central columns and, once their index is complete, the variants, with their totals; given the lakes'
    ladder, the seas' shelves to 200 m take its rungs."""
    part = seas_cases()
    part['totals'] = seas_totals(part, ladder=ladder)
    index, _ = C.seas_index()
    if index.get('complete'):
        part['variants'] = seas_cases(variants=True)
        part['variants_table'] = seas_variants_table(part)
        part['forms'] = seas_forms_range(part, ladder=ladder)
    part['directions'] = seas_directions(part)
    part['seas_product_sha256'] = sha256(SEAS_PRODUCT)
    if part.get('variants'):
        part['variant_totals'] = seas_variant_totals(part, ladder=ladder)
    return part


def ladder_table(kinds, prefixes=('ladder',)):
    """Each rung's (and giant's) emission per m2 of lake surface by band and trophic state, and its parts, for each
    lane, at the central wind share; the giants read from their settled last decade."""
    out = {}
    for name, c in kinds['cases'].items():
        if not name.startswith(prefixes) or c['case'].get('eta', 'central') != 'central':
            continue
        depth = f"{c['case']['depth_m']:g}"
        band = C.BANDS[c['case']['band']]
        trophic = c['case']['coupled'].split('/')[0]
        r = c['result']
        out.setdefault(depth, {}).setdefault(trophic, {})[band] = dict(
            overturns=len(r['case']['overturns']), years=r['case']['years'], settled=r['settled'],
            mixed_depth_m=r['case']['mixed_depth_m'],
            lanes={l['label']: dict(emitted=l['emitted_mol_m2_yr'], at_overturns=l['emitted_at_overturns_mol_m2_yr'],
                                    bubbles=l['bubbles_to_air_mol_m2_yr'], surface=l['diffusive_mol_m2_yr'],
                                    made=l['made_anoxic_mol_m2_yr'], floors=l['floors_mol_m2_yr'],
                                    oxidised=l['oxidised_mol_m2_yr'], oxycline_share=l['oxycline_oxidised_share'],
                                    aom=l['aom_floor_mol_m2_yr'] + l['aom_water_mol_m2_yr'],
                                    inventory=l['inventory_mol_m2']['mean'],
                                    first_order_per_h=l['first_order_max_per_h'],
                                    store_change_share=d['change_share'] if d else None)
                   for l, d in zip(r['lanes'], store_drift(r) or [None] * len(r['lanes']))})
    return out
def gas_setting(st):
    """lake_anoxia's carbonate system and air gases for the floor's bubble point: the held fresh seas' lake pH (sea
    chemistry), the design air's CO2 and the lake's alkalinity, N2 and Ar at saturation."""
    chem = json.loads(A.CHEMISTRY.read_text())
    held = next(r for r in chem['runs'] if all(r[k] == v for k, v in R.HELD_FRESH.items()))
    lake_ph = float(held['rows'][-1]['lake_ph'])
    pco2 = A.DESIGN_CO2_PPM * 1e-6 * A.SURFACE_ATM
    ta = G.alkalinity_at(lake_ph, pco2, 21.15, A.LAKE_IONIC)[0]
    fr = T.earthlike_air(R.DESIGN_PRESSURE_PA, R.DESIGN_CO2_PPM).fractions
    return SimpleNamespace(pco2=pco2, ta=ta, lake_ph=lake_ph, n2_frac=float(fr['N2']), ar_frac=float(fr['Ar']),
                           salinity=st.salinity, o2_fraction=st.o2_fraction, gravity=st.gravity)


def floor_saturation(gs, col, snapshot, j, ch4, share):
    """The floor water's gas saturation ratio (1 is the bubble point) with the tracer's methane: CO2 from the oxygen the
    water has used and the debt's carbon (a mole with each mole of methane, one for each 32 g of the rest), in the
    lake's carbonate system as lake_anoxia reckons it, N2 and Ar at the surface's saturation, and the vapour."""
    t = float(col.t[snapshot, -1])
    ox = max(float(col.oxygen[snapshot, -1, j]), 0.0)
    debt = max(float(col.debt[snapshot, -1, j]), 0.0)
    sat = float(L.oxygen_saturation(t, A.SURFACE_ATM, gs.salinity, gs.o2_fraction))
    added = max(sat - ox, 0.0) / 32.0 + debt * (share / 64.0 + (1.0 - share) / 32.0)
    base = float(A.equilibrium_dic(gs, np.array(t)))
    _, co2 = G.speciate(base + added / 1000.0, gs.ta, t, A.LAKE_IONIC)
    n2, ar = A.air_gases(gs, t)
    p_gas = float(G.gas_pressure_atm(t, ch4=ch4, co2_aq=float(co2), n2=float(n2), ar=float(ar), o2=ox / 32.0))
    p_local = float(G.local_pressure_atm(col.basin.z[-1], A.SURFACE_ATM, gs.gravity))
    return p_gas / p_local, p_gas, p_local


def giant_summary(res, col, gs):
    """A giant's decades from its oxygenated start, per lane: what its anoxic water made, what it emitted (between
    overturns and at them), what it stored, and its floor's approach to the bubble point from the tracer's own
    methane at the last two decades' close."""
    parts = res['parts']
    out = {}
    for j, row in enumerate(res['lanes']):
        label = row['label']
        decades, inv_prev = [], 0.0
        for part in parts:
            r = part['lanes'][j]
            years = part['case']['years']
            inv = r['inventory_mol_m2']['end']
            decades.append(dict(made_mol_m2_yr=r['made_anoxic_mol_m2_yr'], emitted_mol_m2_yr=r['emitted_mol_m2_yr'],
                                at_overturns_mol_m2_yr=r['emitted_at_overturns_mol_m2_yr'],
                                bubbles_to_air_mol_m2_yr=r['bubbles_to_air_mol_m2_yr'],
                                invaded_mol_m2_yr=r.get('invaded_mol_m2_yr', 0.0),
                                oxidised_mol_m2_yr=r['oxidised_mol_m2_yr'], stored_mol_m2_yr=(inv - inv_prev) / years,
                                inventory_end_mol_m2=inv, overturns=len(part['case']['overturns'])))
            inv_prev = inv
        sat = []
        for part in parts[-2:]:
            ch4 = part['closing_mol_m3'][str(j)][-1]
            sat.append(floor_saturation(gs, col, part['closing_snapshot'], col.names.index(row['column']), ch4,
                                        row['share']))
        years_last = parts[-1]['case']['years']
        rate = (sat[-1][0] - sat[0][0]) / years_last if len(sat) == 2 else 0.0
        out[label] = dict(column=row['column'], decades=decades,
                          fifty_years=dict(made_mol_m2_yr=float(np.mean([d['made_mol_m2_yr'] for d in decades])),
                                           emitted_mol_m2_yr=float(np.mean([d['emitted_mol_m2_yr']
                                                                             for d in decades]))),
                          floor_saturation_ratio=sat[-1][0], floor_gas_pressure_atm=sat[-1][1],
                          floor_pressure_atm=sat[-1][2], floor_saturation_rate_per_yr=rate,
                          years_to_bubble_point=(1.0 - sat[-1][0]) / rate if rate > 0 else None)
    return out


# --- the seas ------------------------------------------------------------------------------------------------------

SEA_VARIANT_SUFFIXES = ('_wind_share_low', '_wind_share_high', '_lake_number', '_fixed_dissipation', '_solute_none',
                        '_oligotrophic_physics', '_shih_bouffard')
SEA_COMMUNITIES = (('earth', {}), ('designed', B.FAST), ('poor', B.POOR), ('none', dict(aerobic=False)))


def sea_lanes(supply, deep):
    """Each trophic state's methane under Earth's lake community, the designed and the poor ends and none; in the deep
    columns also the methane share's low and high ends. The held seas carry no sulphate and their floors take a
    thirtieth of the lakes' river iron a square metre (biology.acceptor_supply), so their anaerobic pathways are given
    as that supply's bound rather than traced."""
    lanes = []
    for j in range(len(TRAITS)):
        for label, kw in SEA_COMMUNITIES:
            lanes.append(T.Lane(j, label=label, **dict(B.EARTH, **kw)))
        if deep:
            for s_ in ('low', 'high'):
                lanes.append(T.Lane(j, share=B.SHARES[s_], label=f'share/{s_}', **B.EARTH))
    return lanes


def is_sea_variant(name):
    return any(name.endswith(x) for x in SEA_VARIANT_SUFFIXES)


def seas_cases(variants=False):
    """The seas' columns traced, central (or, with variants, the variants), each case cached."""
    from research.studies.sea_oxygen import setting as SS
    index, digest = C.seas_index()
    st_sea = SS.SeaSetting()
    supply = B.acceptor_supply(waters='seas')
    out = dict(index=repo_path(C.SEAS_RUNS / 'columns.json'), index_sha256=digest, lake_py=index['lake_py']['sha256'],
               complete=bool(index.get('complete')), acceptors=supply, columns={})
    for name, entry in index['cases'].items():
        if is_sea_variant(name) != variants:
            continue
        col = C.seas_column(st_sea, name, index)
        deep = entry['kind'] == 'deep'
        if not variants:
            # the central columns were traced first under the earlier settle rule, every lane to 0.5%: they keep it
            res = trace_case(col, sea_lanes(supply, deep), f'sea_{name}', cyclic=True, max_passes=4, settle_rule='all')
        else:
            form = col.form != 'osborn'
            lanes = [T.Lane(j, label=k, **dict(B.EARTH, **kw)) for j in range(len(TRAITS))
                     for k, kw in SEA_COMMUNITIES if form or k in ('earth', 'none')]
            res = trace_case(col, lanes, f'sea_{name}', cyclic=True, max_passes=4)
        out['columns'][name] = dict(basin=entry['basin']['key'], kind=entry['kind'], column=entry['case']['column'],
                                    pan_m=entry['case'].get('pan_m'), stands_for=entry['basin']['stands_for'],
                                    lat=entry.get('lat'),
                                    diffusivity=col.form,
                                    physics=entry['physics'], files=entry['files'], summary=entry['summary'],
                                    summary_sha256=entry['summary_sha256'], case=res['case'], lanes=res['lanes'],
                                    passes=res['passes'], settled=res['settled'], bubbles=res['bubbles'],
                                    store=store_drift(res))
        del col
    return out


def store_drift(res):
    """Each lane's methane store over the last pass of a cyclic case (mol m-2 at its start and end, and the change a
    year): a store still shrinking makes the lane's emission an upper value, one still growing a lower one."""
    if not res.get('passes'):
        return None
    last, years = res['passes'][-1], res['case']['years']
    return [dict(label=l['label'], column=l['column'], start_mol_m2=a, end_mol_m2=b, change_mol_m2_yr=(b - a) / years,
                 change_share=(b - a) / a if a > 0 else None)
            for l, a, b in zip(res['lanes'], last['start_mol_m2'], last['end_mol_m2'])]


def lane_of(column, trophic, label):
    return next(l for l in column['lanes'] if l['label'] == label and l['column'].startswith(trophic + '/'))


SHELF_RUNGS = ((30.0, 'pan10'), (math.sqrt(30.0 * 50.0), '30'), (math.sqrt(50.0 * 100.0), '50'),
               (math.sqrt(100.0 * 352.0), '100'), (200.0, '352'))


def seas_shelves(seas):
    """Each basin's cumulative shelf shares of the water it stands for: its own water at its own shares (seas.json's
    run basins), and the water it stands for beyond its own (the other regions and the crater seas) at the shares the
    seas' whole leaves once every basin's own shelves are taken out, so that the seas' total keeps seas.json's."""
    total = seas['seas']
    area = total['area_km2']
    whole = {float(k): v for k, v in total['shelf_area_share'].items()}
    basins = {b['key']: b for b in seas['run_basins']}
    own = sum(b['stands_for']['own_km2'] for b in basins.values())
    rest = {d: (whole[d] * area - sum(b['shelf_area_share'][f'{d:g}'] * b['stands_for']['own_km2']
                                      for b in basins.values())) / (area - own) for d in whole}
    out = {}
    for key, b in basins.items():
        a_own = b['stands_for']['own_km2']
        a_rest = b['stands_for']['other_regions_km2'] + b['stands_for']['crater_seas_km2']
        out[key] = {d: (b['shelf_area_share'][f'{d:g}'] * a_own + rest[d] * a_rest) / (a_own + a_rest)
                    for d in whole}
    return out, rest


def shelf_at(shares, depth):
    """A cumulative shelf share at a depth, interpolated in the logarithm of depth between seas.json's."""
    xs = sorted(shares)
    if depth <= xs[0]:
        return shares[xs[0]] * depth / xs[0]
    if depth >= xs[-1]:
        return shares[xs[-1]]
    return float(np.interp(math.log(depth), np.log(xs), [shares[x] for x in xs]))


def sea_parts(shares):
    """A basin's water by local depth: shallower than 5 m at the 3 m pan, 5-30 m at the 10 m pan, 30-200 m at the
    lakes' ladder rungs nearest its depth (the geometric midpoints between rungs), deeper at the deep column."""
    parts, below = dict(pan3=shelf_at(shares, 5.0)), shelf_at(shares, 5.0)
    for edge, name in SHELF_RUNGS:
        upto = shelf_at(shares, edge)
        parts[name] = parts.get(name, 0.0) + upto - below
        below = upto
    parts['deep'] = 1.0 - below
    return parts


def sea_band(lat):
    a = abs(float(lat))
    return C.BANDS[0] if a < 15.0 else (C.BANDS[1] if a < 30.0 else C.BANDS[2])


def rung_flux(ladder, rung, t, band, label):
    """A ladder rung's emission in a band (the 0-15 band's where that band was not run), or None."""
    runs = ladder.get(rung, {}).get(t, {})
    v = runs.get(band) or runs.get(C.BANDS[0])
    return None if v is None or label not in v['lanes'] else v['lanes'][label]['emitted']


def pan_flux(pan, ladder, t, band, label, rung):
    """A pan's emission; for a lane the pans were not run with (the methane share's ends), the pan's under Earth's
    community times the ratio the lakes' rung of the same depth shows between that lane and Earth's community."""
    if pan is None:
        return None
    lane_ = next((l for l in pan['lanes'] if l['label'] == label and l['column'].startswith(t + '/')), None)
    if lane_ is not None:
        return lane_['emitted_mol_m2_yr']
    earth = lane_of(pan, t, 'earth')['emitted_mol_m2_yr']
    a, b = rung_flux(ladder, rung, t, band, label), rung_flux(ladder, rung, t, band, 'earth')
    return earth * (a / b if a is not None and b else 1.0)


def basin_per_m2(parts, deep, pan3, pan10, ladder, t, band, label):
    """A basin's emission per m2 of the water it stands for, by its parts; without the ladder (the seas alone) the
    shelves deeper than 30 m take the deep column's flux."""
    f = deep * parts['deep']
    f += pan3 * parts['pan3'] + pan10 * parts['pan10']
    for _, rung in SHELF_RUNGS[1:]:
        r = rung_flux(ladder, rung, t, band, label) if ladder else None
        f += (deep if r is None else r) * parts[rung]
    return f


def seas_totals(part, seas_json=None, ladder=None):
    """The seas' methane, Tg a year, weighted by what each basin's columns stand for: each basin's water by local depth
    (sea_parts, each basin's own shelf shares), its typical deep column beyond 200 m; beside it the deep columns over
    all the water and each basin's coldest columns in place of its typical ones."""
    seas_doc = json.loads((seas_json or SEAS_JSON).read_text())
    seas = seas_doc['seas']
    area = seas['area_km2'] * 1e6
    shelves, rest = seas_shelves(seas_doc)
    cols = part['columns']
    basins = sorted({c['basin'] for c in cols.values()}, key=lambda b: -next(
        c['stands_for']['share_of_the_seas'] for c in cols.values() if c['basin'] == b))
    labels = [l['label'] for l in next(iter(cols.values()))['lanes'] if l['column'].startswith(TRAITS[0] + '/')]
    out = dict(area_km2=seas['area_km2'], shelf_share=seas['shelf_area_share'], shelf_share_beyond_own=rest,
               ladder_for_shelves=bool(ladder), basins={}, totals={})
    tg = lambda mol_m2_yr, a: mol_m2_yr * M_CH4 * a / TG
    for b in basins:
        mine = {k: c for k, c in cols.items() if c['basin'] == b}
        w = next(iter(mine.values()))['stands_for']['share_of_the_seas']
        typ, cold = mine.get(f'{b}_typical'), mine.get(f'{b}_coldest')
        p3, p10 = mine.get(f'{b}_typical_pan3m'), mine.get(f'{b}_typical_pan10m')
        p3c = mine.get(f'{b}_coldest_pan3m') or p3
        band = sea_band((typ or next(iter(mine.values())))['lat'])
        parts = sea_parts(shelves[b])
        rows = {}
        for t in TRAITS:
            entry = {}
            for label in labels:
                deep_lane = next((l for l in typ['lanes'] if l['label'] == label and l['column'].startswith(t + '/')),
                                 None) if typ else None
                if deep_lane is None:
                    continue
                f_typ = deep_lane['emitted_mol_m2_yr']
                f3 = pan_flux(p3, ladder or {}, t, band, label, '5')
                f10 = pan_flux(p10, ladder or {}, t, band, label, '12')
                row = dict(deep_mol_m2_yr=f_typ, pan3_mol_m2_yr=f3, pan10_mol_m2_yr=f10,
                           columns_alone_tg_yr=tg(f_typ, w * area))
                if f3 is not None and f10 is not None:
                    per = basin_per_m2(parts, f_typ, f3, f10, ladder, t, band, label)
                    lower = basin_per_m2(parts, f_typ, f3, f10, None, t, band, label)
                    pans = f3 * parts['pan3'] + f10 * parts['pan10']
                    row.update(with_shelves_mol_m2_yr=per, with_shelves_tg_yr=tg(per, w * area),
                               shelves_tg_yr=tg(per - f_typ * parts['deep'], w * area),
                               to_30m_tg_yr=tg(pans, w * area),
                               from_30_to_200m_tg_yr=tg(per - pans - f_typ * parts['deep'], w * area),
                               with_shelves_deep_columns_beyond_30m_tg_yr=tg(lower, w * area))
                if cold is not None:
                    c_lane = next((l for l in cold['lanes'] if l['label'] == label and l['column'].startswith(t + '/')),
                                  None)
                    if c_lane is not None:
                        f3c = pan_flux(p3c, ladder or {}, t, band, label, '5')
                        row['coldest_deep_mol_m2_yr'] = c_lane['emitted_mol_m2_yr']
                        if f3c is not None and f10 is not None:
                            row['coldest_with_shelves_tg_yr'] = tg(basin_per_m2(
                                parts, c_lane['emitted_mol_m2_yr'], f3c, f10, ladder, t, band, label), w * area)
                entry[label] = row
            rows[t] = entry
        out['basins'][b] = dict(stands_for=w, band=band, parts=parts, columns=sorted(mine), traits=rows)
    for t in TRAITS:
        entry = {}
        for label in labels:
            rows = [out['basins'][b]['traits'][t][label] for b in basins if label in out['basins'][b]['traits'][t]]
            entry[label] = dict(
                with_shelves_tg_yr=float(sum(r.get('with_shelves_tg_yr', r['columns_alone_tg_yr']) for r in rows)),
                with_shelves_deep_columns_beyond_30m_tg_yr=float(sum(r.get(
                    'with_shelves_deep_columns_beyond_30m_tg_yr', r['columns_alone_tg_yr']) for r in rows)),
                columns_alone_tg_yr=float(sum(r['columns_alone_tg_yr'] for r in rows)),
                shelves_tg_yr=float(sum(r.get('shelves_tg_yr', 0.0) for r in rows)),
                to_30m_tg_yr=float(sum(r.get('to_30m_tg_yr', 0.0) for r in rows)),
                from_30_to_200m_tg_yr=float(sum(r.get('from_30_to_200m_tg_yr', 0.0) for r in rows)),
                coldest_with_shelves_tg_yr=float(sum(r.get('coldest_with_shelves_tg_yr',
                                                           r.get('with_shelves_tg_yr', r['columns_alone_tg_yr']))
                                                     for r in rows)))
        out['totals'][t] = entry
    return out


def central_of(name):
    """The central column a variant varies: its name without the variant's suffix."""
    for x in SEA_VARIANT_SUFFIXES:
        if name.endswith(x):
            return name[:-len(x)]
    return name


def seas_variants_table(part):
    """Each variant column's emission beside its central column's, per trophic state and lane: the surface flux, the
    bubbles, the whole and the whole's ratio to the central's."""
    out = {}
    for name, v in part.get('variants', {}).get('columns', {}).items():
        base = part['columns'].get(central_of(name))
        if base is None:
            continue
        rows = {}
        for l in v['lanes']:
            t = l['column'].split('/')[0]
            b = lane_of(base, t, l['label'])
            rows.setdefault(t, {})[l['label']] = dict(
                emitted=l['emitted_mol_m2_yr'], surface=l['diffusive_mol_m2_yr'], bubbles=l['bubbles_to_air_mol_m2_yr'],
                made=l['made_anoxic_mol_m2_yr'], central_emitted=b['emitted_mol_m2_yr'],
                central_surface=b['diffusive_mol_m2_yr'], central_bubbles=b['bubbles_to_air_mol_m2_yr'],
                central_made=b['made_anoxic_mol_m2_yr'],
                ratio=l['emitted_mol_m2_yr'] / b['emitted_mol_m2_yr'] if b['emitted_mol_m2_yr'] > 0 else None)
        out[name] = dict(central=central_of(name), basin=v['basin'], kind=v['kind'], diffusivity=v['diffusivity'],
                         physics=v['physics'], settled=v['settled'], traits=rows)
    return out


def seas_variant_totals(part, seas_json=None, ladder=None):
    """The seas' Tg a year with each basin's typical deep column and pans replaced by a variant wherever the variant
    was run (the central columns elsewhere), on the same partition by local depth, for each community the variants
    carry: the wind share's brackets over every basin, and the other variants over the basins and pans they were run
    for."""
    seas_doc = json.loads((seas_json or SEAS_JSON).read_text())
    area = seas_doc['seas']['area_km2'] * 1e6
    cols, var = part['columns'], part.get('variants', {}).get('columns', {})
    basins = part['totals']['basins']
    out = {}
    for suffix in SEA_VARIANT_SUFFIXES:
        if not any(n.endswith(suffix) for n in var):
            continue
        labels = sorted({l['label'] for n, v in var.items() if n.endswith(suffix) for l in v['lanes']})
        replaced = sorted(n for n in var if n.endswith(suffix))
        rows = {}
        for t in TRAITS:
            entry = {}
            for label in labels:
                total = 0.0
                for b, row in basins.items():
                    pick = lambda base: var.get(base + suffix) or cols.get(base)
                    deep, p3, p10 = (pick(f'{b}_typical{x}') for x in ('', '_pan3m', '_pan10m'))
                    if deep is None:
                        continue
                    f = lane_of(deep, t, label)['emitted_mol_m2_yr']
                    if p3 is not None and p10 is not None:
                        f = basin_per_m2(row['parts'], f, lane_of(p3, t, label)['emitted_mol_m2_yr'],
                                         lane_of(p10, t, label)['emitted_mol_m2_yr'], ladder, t, row['band'], label)
                    total += f * M_CH4 * row['stands_for'] * area / TG
                entry[label] = total
            rows[t] = entry
        out[suffix.lstrip('_')] = dict(replaced=replaced, tg_yr=rows)
    return out


METHANE_IS = {'an upper value': 'a lower value', 'a lower value': 'an upper value',
              'its deep water anoxic': 'a lower value'}


def seas_directions(part, product=None):
    """The seas' product's reading of each column at its analysis step (answers.<basin>.<column>.analysis_step) where
    it is not settled, and what that makes of the traced methane: oxygen an upper value makes the methane a lower one
    and the reverse; a column cut off from the surface and still ageing, its methane still gathering, a lower one."""
    product = product or json.loads(SEAS_PRODUCT.read_text())
    out = {}
    both = list(part['columns'].items()) + list(part.get('variants', {}).get('columns', {}).items())
    for name, c in both:
        b = c['basin']
        a = product['answers'].get(b, {}).get(name[len(b) + 1:], {}).get('analysis_step')
        if a is None or a.get('reading') == 'settled':
            continue
        out[name] = dict(reading=a['reading'], oxygen_is=a['oxygen_is'], age_is=a['age_is'],
                         passes=a['passes'], passes_overturning=a['passes_overturning'],
                         methane_is=METHANE_IS.get(a['oxygen_is']))
    return out


SEA_FORM_LABELS = ('earth', 'designed', 'poor', 'none')


def seas_forms_range(part, seas_json=None, ladder=None):
    """The seas' methane under Shih and Bouffard's form beside Osborn's, Tg a year. The basins whose typical deep
    column was run under both forms take their own traced rows for their deep water; every other basin's deep water
    takes its own emission plus the change the form makes per m2 in the run basins' columns, in the same trophic state
    and community (their mean weighted by what they stand for, and each alone for the bracket), never below nothing;
    the shelves to 200 m stay as under Osborn's, as the lakes' rungs within the ventilated depth do."""
    shih = {v['basin']: v for v in part.get('variants', {}).get('columns', {}).values()
            if v['diffusivity'] == 'shih_bouffard' and v['kind'] == 'deep' and v['column'] == 'typical'}
    if not shih:
        return None
    seas_doc = json.loads((seas_json or SEAS_JSON).read_text())
    area = seas_doc['seas']['area_km2'] * 1e6
    cols, totals = part['columns'], part['totals']
    w = {b: shih[b]['stands_for']['share_of_the_seas'] for b in shih}
    out = dict(run_basins={b: dict(stands_for=w[b], column=f'{b}_typical_shih_bouffard') for b in sorted(shih)},
               traits={})
    for t in TRAITS:
        rows = {}
        for label in SEA_FORM_LABELS:
            change = {b: lane_of(v, t, label)['emitted_mol_m2_yr'] - lane_of(cols[f'{b}_typical'], t, label)[
                'emitted_mol_m2_yr'] for b, v in shih.items()}
            weighted = sum(w[b] * change[b] for b in shih) / sum(w.values())
            schemes = dict(weighted=weighted, **{f'by_{b}': change[b] for b in sorted(shih)})
            est = {}
            for scheme, d in schemes.items():
                total = 0.0
                for b, row in totals['basins'].items():
                    typ = cols.get(f'{b}_typical')
                    if typ is None:
                        continue
                    if b in shih:
                        deep = lane_of(shih[b], t, label)['emitted_mol_m2_yr']
                    else:
                        deep = max(lane_of(typ, t, label)['emitted_mol_m2_yr'] + d, 0.0)
                    cell = row['traits'][t][label]
                    p3, p10 = cell.get('pan3_mol_m2_yr'), cell.get('pan10_mol_m2_yr')
                    per = basin_per_m2(row['parts'], deep, p3, p10, ladder, t, row['band'], label) \
                        if p3 is not None else deep
                    total += per * M_CH4 * row['stands_for'] * area / TG
                est[scheme] = total
            osb = totals['totals'][t][label]['with_shelves_tg_yr']
            vals = [osb] + list(est.values())
            rows[label] = dict(osborn_tg_yr=osb, shih_tg_yr=est['weighted'],
                               shih_bracket_tg_yr={k: v for k, v in est.items() if k != 'weighted'},
                               change_mol_m2_yr=dict(weighted=weighted, **change),
                               range_tg_yr=[min(vals), max(vals)])
        out['traits'][t] = rows
    return out


def air_shim(absorbed):
    """lake_anoxia's air() with its ledger, the reduced-gas forcing and one target's absorbed sunlight, so that the
    steady methane and the dimming are its arithmetic exactly."""
    st = SimpleNamespace(ledgers=json.loads(LEDGERS.read_text()),
                         forcing=json.loads(GAS_FORCING.read_text())['moon_at_design'],
                         absorbed={t: absorbed for t in A.TARGET_ORDER})
    return st


def air_for(absorbed, tg_by_trait):
    """The air's steady methane and the added dimming for the lakes' (and seas') Tg a year by trait, each with its
    low, central and high case, through lake_anoxia's air()."""
    st = air_shim(absorbed)
    rows = {t: dict(moon_tg_yr=tg_by_trait[t]) for t in tg_by_trait}
    budget = {target: dict(traits=rows) for target in A.TARGET_ORDER}
    out = A.air(st, budget)
    hold = out['soil_sink_tg_yr_per_ppm'] * out['held_ppm']                # Tg a year the soils take at 1.87 ppm
    traits = out['targets']['resort']
    for row in traits.values():
        for case in row.values():
            # the waters' room under the held methane, Tg a year, at the ledger's low- and high-cost ends: what the
            # soils take at the held methane less every source left (negative, the removal that holds it)
            case['lakes_uncut']['room_under_held_tg_yr'] = tuple(hold - x for x in case['lakes_uncut']['sources_tg_yr'])
    return dict(ledger_lakes_tg_yr=out['ledger_lakes_tg_yr'], held_ppm=out['held_ppm'],
                held_forcing_w_m2=out['held_forcing_w_m2'], soil_sink_tg_yr_per_ppm=out['soil_sink_tg_yr_per_ppm'],
                absorbed_sunlight_w_m2=absorbed, traits=traits)


def absorbed_sunlight(share):
    """lake_anoxia's reading: between the design case's 211.2 W/m2 and the warm Moon's 225.4 by the share of the
    way."""
    return 211.2 + share * (225.4 - 211.2)




# --- the deep mixing's form: Shih and Bouffard's beside Osborn's ----------------------------------------------------

def forms_cases(st, only=None):
    """lake_oxygen's mixing-forms columns traced: the 352 m 0-15 lake at 301 K in each trophic state (its own
    physics) and the 1,925 m giant from it (its settled last decade), under Shih and Bouffard's form; the same lake in
    the lock, the resort at 300 K and the hot world under it (the four central columns their mesotrophic physics
    carries); Osborn's rule with the wind's dissipation bounded, in the lock."""
    index, digest = C.mixing_forms_index()
    out = dict(index=repo_path(C.MIXING_FORMS_RUNS / 'mixing_forms_index.json'), index_sha256=digest, cases={})
    for name, entry in index['cases'].items():
        if only is not None and not only(name):
            continue
        use = last_pass(index, name) if 'giant' in name else index
        col = C.mixing_forms_column(name, use)
        lanes = community_lanes(range(len(col.names)), shares=False)
        res = trace_case(col, lanes, f'forms_{name}', cyclic=True, max_passes=GIANT_PASSES if 'giant' in name else 4)
        out['cases'][name] = dict(case=entry['case'], form=col.form, columns=col.names, result=res)
        del col
    return out


def forms_ratios(forms, ladder, targets):
    """Shih and Bouffard's emission against Osborn's for the same lake: at 301 K the 352 m 0-15 lake and the 1,925 m
    giant against the kinds' ladder, each trophic state with its own physics; in the three targets the 352 m 0-15 lake
    against lake_anoxia's regenerated column. Per lane, the surface flux, the bubbles and the whole."""
    out = dict(at_301k={}, targets={})
    for t in TRAITS:
        for rung, name in (('352', f'shih_301k_352m_0-15_{t}'), ('1925', f'shih_301k_giant_1925m_0-15_{t}')):
            shih = forms['cases'].get(name)
            osb = ladder.get(rung, {}).get(t, {}).get('0-15')
            if shih is None or osb is None:
                continue
            rows = {}
            for l in shih['result']['lanes']:
                o = osb['lanes'].get(l['label'])
                if o is None:
                    continue
                rows[l['label']] = dict(
                    shih=dict(emitted=l['emitted_mol_m2_yr'], surface=l['diffusive_mol_m2_yr'],
                              bubbles=l['bubbles_to_air_mol_m2_yr'], made=l['made_anoxic_mol_m2_yr']),
                    osborn=dict(emitted=o['emitted'], surface=o['surface'], bubbles=o['bubbles'], made=o['made']),
                    surface_ratio=(l['diffusive_mol_m2_yr'] / o['surface']) if o['surface'] > 0 else None)
            out['at_301k'].setdefault(t, {})[rung] = rows
    for target in A.TARGET_ORDER:
        name = f'shih_{target}_352m_0-15'
        shih = forms['cases'].get(name)
        osb = targets['deep'].get(f'deep_{target}_0-15')
        if shih is None or osb is None:
            continue
        rows = {}
        for t in TRAITS:
            for label in ('earth', 'designed', 'poor', 'none'):
                s_ = next((l for l in shih['result']['lanes']
                           if l['label'] == label and l['column'].startswith(t + '/')), None)
                o = next((l for l in osb['lanes'] if l['label'] == label and l['column'].startswith(t + '/')), None)
                if s_ is None or o is None:
                    continue
                rows.setdefault(t, {})[label] = dict(shih=s_['emitted_mol_m2_yr'], osborn=o['emitted_mol_m2_yr'],
                                                     shih_surface=s_['diffusive_mol_m2_yr'],
                                                     osborn_surface=o['diffusive_mol_m2_yr'])
        out['targets'][target] = rows
    return out


def shih_scaled_by_rung(by_rung, ratios):
    """Shih and Bouffard's flux by rung, carried from its rows: each rung's emission plus the change the form makes per
    m2 where it was run, in the same trophic state and lane, the 352 m lake's for the 352 m rung and the 1,925 m
    giant's for the giants, never below nothing. The rungs to 100 m, within the ventilated depth the nights reach, are
    as under Osborn's. A change carried by addition stays bounded where Osborn's flux is near nothing, as it is in a
    ventilated column, which a ratio would not."""
    out = {}
    for rung, entry in by_rung.items():
        src = '1925' if float(rung) >= 1000.0 else ('352' if float(rung) >= 352.0 else None)
        out[rung] = {}
        for t, e in entry.items():
            rows = ratios['at_301k'].get(t, {}).get(src, {}) if src else {}
            per, changes = {}, {}
            for label, f in e['emitted_mol_m2_yr'].items():
                row = rows.get(label)
                d = row['shih']['emitted'] - row['osborn']['emitted'] if row else 0.0
                per[label] = max(f + d, 0.0)
                changes[label] = d
            out[rung][t] = dict(emitted_mol_m2_yr=per, change_mol_m2_yr=changes, change_from=src if rows else None)
    return out


def air_by_community(absorbed, tg):
    """The air for each community: lake_anoxia's air() with the central (Earth's community) as its central case, the
    designed community as its low and the poor one as its high, per trophic state; tg: {trait: {label: Tg a year}}."""
    rows = {t: dict(low=v['designed'], central=v['earth'], high=v['poor']) for t, v in tg.items()}
    return air_for(absorbed, rows)


def resort_absorbed(st301):
    return absorbed_sunlight(st301.share['resort'])


def producer_files():
    from research.studies.lake_kinds import setting as KS
    from research.studies.sea_oxygen import columns as SC
    from research.studies.sea_oxygen import setting as SS
    from research.studies.lake_oxygen import mixing_forms as MF
    paths = [HERE / f for f in ('run.py', 'cases.py', 'tracer.py', 'columns.py', 'biology.py')] + \
        [Path(L.__file__), Path(G.__file__), Path(A.__file__), Path(R.__file__), Path(MF.__file__), Path(KS.__file__),
         Path(SS.__file__), Path(SC.__file__)]
    return {str(p.resolve().relative_to(ROOT)): sha256(p) for p in paths}


def inputs_of():
    paths = [ANOXIA, C.ANOXIA_RUNS / 'profiles.json', C.KINDS_RUNS / 'ladder.json', C.SEAS_RUNS / 'columns.json',
             LAKES, SEAS_JSON, SEAS_PRODUCT, CENSUS, CHEMISTRY, LEDGERS, GAS_FORCING, TARGETS, HERE / 'sources.json']
    return {str(p.relative_to(ROOT)): sha256(p) for p in paths}


# --- tracing in several processes ---------------------------------------------------------------------------------

CLAIMS = RUNS / 'claims'


def claim(name):
    """A process's claim on a case: the first to create its claim file traces it, and the others pass it by."""
    CLAIMS.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(CLAIMS / f'{name}.claim', os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, 'w') as f:
        f.write(f'{os.getpid()} {time.strftime("%Y-%m-%d %H:%M:%S")}\n')
    return True


def pending():
    """Every case the full run traces, the slowest first: lake_anoxia's 1,925 m lakes, the giants and the meromictic
    lakes, the mixing forms' giants and the rest of them, the ladder with its detailed rungs and its deepest first,
    lake_anoxia's 352 m lakes; the reproduction is traced by the full run itself."""
    meta = C.anoxia_meta(C.ANOXIA_RUNS)
    kinds = C.kinds_index()['cases']
    forms = C.mixing_forms_index()[0]['cases']
    order = [('targets', n) for n, c in meta['cases'].items() if c['kind'] != 'deep']
    order += [('kinds', n) for n in kinds if n.startswith(('giant', 'meromictic'))]
    order += [('forms', n) for n in forms if 'giant' in n] + [('forms', n) for n in forms if 'giant' not in n]
    rest = [n for n in kinds if not n.startswith(('giant', 'meromictic'))]
    detail = lambda n: kinds[n]['case']['band'] == DETAIL_BAND and \
        float(kinds[n]['case']['depth_m']) in LADDER_DETAIL_M + KINETICS_M
    order += [('kinds', n) for n in sorted(rest, key=lambda n: (not detail(n), -float(kinds[n]['case']['depth_m'])))]
    order += [('targets', n) for n, c in meta['cases'].items() if c['kind'] == 'deep']
    return order


def work():
    """One of several processes filling the cache: each pending case is traced by the process that claims it, with
    the same calls the full run makes, so the full run then finds every case cached."""
    st = R.Setting()
    product = json.loads(ANOXIA.read_text())
    from research.studies.lake_kinds import setting as KS
    st301 = KS.Setting301()
    index = C.kinds_index()
    supply = B.acceptor_supply(waters='lakes')
    mine = 0
    for kind, name in pending():
        if not claim(f'{kind}_{name}'):
            continue
        only = {name}.__contains__
        if kind == 'targets':
            anoxia_targets(st, product, only=only)
        elif kind == 'kinds':
            kinds_cases(st301, index, supply, only=only)
        else:
            forms_cases(st, only=only)
        mine += 1
    print('claimed', mine, 'cases', flush=True)


EVIDENCE = (
    'Arithmetic on one lake model\'s saved columns: biosphere/ecology/lake.py (the final model, wind-driven deep '
    'mixing) run by lake_anoxia, the lake kinds study and the sea oxygen study, whose daily columns and step-by-step '
    'series the tracer reads without running the model. Methane is a tracer of its own, made by lake_anoxia\'s rule in'
    ' anoxic water and on every floor, moved by the model\'s own mixing and diffusivity, bubbled by lake_anoxia\'s '
    'porewater balance and rising bubble, oxidised by Michaelis-Menten methanotrophs whose rates are Earth\'s measured'
    ' ranges carried to engineered communities as design variables, its oxygen accounted against the model\'s '
    'reoxidation of its debt, and exchanged with the air at the model\'s own transfer velocity. The transport is '
    'checked against the model\'s own age and debt, the kinetics against an analytic profile, the bubbles against '
    'lake_anoxia\'s physics and Lake Soppen. The anaerobic pathways are held to their acceptors\' supply from the '
    'land\'s erosion and Earth\'s rates. No lunar lake or sea has been measured; the air\'s steady methane is '
    'lake_anoxia\'s arithmetic on the domain synthesis\'s ledger.')
READING_RULE = (
    'Fluxes in mol per m2 of lake (or sea) surface a year unless named per m2 of floor; Tg of CH4 a year Moon-wide; t'
    ' a year per km2 for the meromictic lakes. Lanes: earth is Earth\'s lake community at Earth\'s central rates '
    '(biology.EARTH), designed and poor the ends of Earth\'s ranges a design reaches for or must allow, none no '
    'methanotrophs (its traced transient; all it makes would leave once its store filled), share/low and share/high '
    'the methane share\'s 0.3 and 0.8, kinetics/* one parameter at an end of its range, overturn/<cells per mL>/<growth'
    ' per day> the stock the water keeps and the response time, aom/* the anaerobic pathways at their acceptors\' '
    'supply. emitted is diffusive and bubbles together; at_overturns within an overturn\'s window (its first full '
    'mixing to thirty days past its last); overturns_summed sums a lane\'s overturns. passes give each pass\'s stores '
    'and, where recorded, its emission and methane made. Rungs are the kinds\' ladder and giants by each lake\'s '
    'deepest point; kinds.moon_tg_yr splits the lakes between the 100 m and 352 m rungs at the lake kinds census\'s '
    'ventilated depth, and kinds.ventilation gives the split at either rung and by the nearest rung; kinds.giants '
    'reads each giant\'s deep oxygen. Seas by basin weight each basin by what its columns stand for, its water by '
    'local depth from its own shelf shares: the 3 m pan to 5 m, the 10 m pan to 30 m, the lakes\' ladder rungs at the '
    'basin\'s band to 200 m (upper-leaning; with_shelves_deep_columns_beyond_30m gives the deep columns there '
    'instead) and the deep column beyond; aom/all/central is every anaerobic pathway at the high-erosion end '
    '(setting.anaerobic_sets); '
    'seas.directions carries the seas\' study\'s reading of its unsettled columns and what it makes of their methane. '
    'Shih and Bouffard\'s figures add the change the form makes per m2 where it was run. air gives lake_anoxia\'s air()'
    ' for lakes and seas together: lakes_uncut is the waters at these fluxes with every other source cut 90%, and '
    'room_under_held_tg_yr the Tg a year the waters could add before the held methane is passed (negative, the '
    'removal that holds it); lakes_and_seas_share carries the methane share\'s 0.3, 0.5 and 0.8 as low, central and '
    'high.')


def lean(obj):
    """The product's copy of the traced summaries: each lane's overturn-by-overturn list summed (with its count) and the
    closing columns left in the cache; everything else as traced."""
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            summed = isinstance(v, list) and (not v or (isinstance(v[0], dict) and 'stored_mol_m2' in v[0]))
            if k == 'events' and summed:
                out['overturns_summed'] = dict({x: float(sum(e[x] for e in v)) for x in (v[0] if v else {})},
                                               count=len(v))
            elif k != 'closing_mol_m3':
                out[k] = lean(v)
        return out
    if isinstance(obj, list):
        return [lean(x) for x in obj]
    return obj


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--checks', action='store_true', help="lake_anoxia's saved columns: the reproduction")
    ap.add_argument('--seas', action='store_true', help="the seas' columns")
    ap.add_argument('--targets', action='store_true', help="lake_anoxia's regenerated columns, the three targets")
    ap.add_argument('--work', action='store_true', help='one of several processes filling the cache')
    ap.add_argument('--cache-only', action='store_true', help='the product from the cache; stop at a case not cached')
    args = ap.parse_args(argv)
    if args.cache_only:
        os.environ[R.CACHE_ONLY] = '1'
    tic = time.time()
    if args.work:
        work()
        print('worked', f'{time.time() - tic:.0f} s', flush=True)
        return 0
    st = R.Setting()
    if args.checks:
        part = anoxia_reproduction(st, json.loads(ANOXIA.read_text()))
        (RUNS / 'checks.json').write_text(dumps(rnd(part)) + '\n')
        print('wrote', RUNS / 'checks.json', f'in {time.time() - tic:.0f} s')
        return 0
    if args.seas:
        part = seas_part()
        (RUNS / 'seas.json').write_text(dumps(rnd(part)) + '\n')
        print('wrote', RUNS / 'seas.json', f'in {time.time() - tic:.0f} s')
        return 0
    product = json.loads(ANOXIA.read_text())
    reproduction = anoxia_reproduction(st, product)
    targets = anoxia_targets(st, product)
    if args.targets:
        (RUNS / 'targets.json').write_text(dumps(rnd(targets)) + '\n')
        print('wrote', RUNS / 'targets.json', f'in {time.time() - tic:.0f} s')
        return 0
    from research.studies.lake_kinds import setting as KS
    st301 = KS.Setting301()
    index = C.kinds_index()
    supply = B.acceptor_supply(waters='lakes')
    kinds = kinds_cases(st301, index, supply)
    ladder = ladder_table(kinds, prefixes=('ladder', 'giant'))
    mero = long_table(kinds, 'meromictic')
    bodies = lake_bodies()
    rungs = rung_table(bodies, LADDER_M + GIANTS_M)
    labels = ('earth', 'designed', 'poor', 'none', 'share/low', 'share/high', 'aom/all/central', 'aom/all/high')
    by_rung = flux_by_rung(ladder, rungs, labels)
    area = float(bodies['area'].sum() / 1e6)
    vent = ventilation_split(bodies, by_rung, rungs, area)
    lakes_tg = vent['moon_tg_yr']['census']
    forms = forms_cases(st)
    ratios = forms_ratios(forms, ladder, targets)
    shih_by_rung = shih_scaled_by_rung(by_rung, ratios)
    vent_shih = ventilation_split(bodies, shih_by_rung, rungs, area)
    lakes_tg_shih = vent_shih['moon_tg_yr']['census']
    seas = seas_part(ladder)
    seas_tg = {t: {label: v['with_shelves_tg_yr'] for label, v in seas['totals']['totals'][t].items()}
               for t in TRAITS}
    both = {t: {label: lakes_tg[t]['tg_yr'].get(label, 0.0) + seas_tg[t].get(label, 0.0)
                for label in ('earth', 'designed', 'poor')} for t in TRAITS}
    seas_tg_shih = {t: {label: v['shih_tg_yr'] for label, v in seas['forms']['traits'][t].items()}
                    for t in TRAITS} if seas.get('forms') else seas_tg
    both_shih = {t: {label: lakes_tg_shih[t]['tg_yr'].get(label, 0.0) + seas_tg_shih[t].get(label, 0.0)
                     for label in ('earth', 'designed', 'poor')} for t in TRAITS}
    # Earth's community at the methane share's ends (0.3 and 0.8) and its centre (0.5), lakes and seas together
    both_share = {t: dict(low=lakes_tg[t]['tg_yr']['share/low'] + seas_tg[t]['share/low'],
                          central=lakes_tg[t]['tg_yr']['earth'] + seas_tg[t]['earth'],
                          high=lakes_tg[t]['tg_yr']['share/high'] + seas_tg[t]['share/high']) for t in TRAITS}
    giants = giants_reading(ladder, forms, rungs)
    # Earth's community with every anaerobic pathway at the high-erosion end (biology.aom_lanes' all/central, the
    # earliest erosion with highland soils) and with mare soils (all/high), the seas as traced
    both_aom = {t: dict(low=lakes_tg[t]['tg_yr'].get('aom/all/high', 0.0) + seas_tg[t].get('earth', 0.0),
                        central=lakes_tg[t]['tg_yr'].get('aom/all/central', 0.0) + seas_tg[t].get('earth', 0.0),
                        high=lakes_tg[t]['tg_yr'].get('earth', 0.0) + seas_tg[t].get('earth', 0.0)) for t in TRAITS}
    # the ventilated depth's bracket with Earth's community: the split at 100 m, at the census's depth and at 352 m
    both_vent = {t: {case: vent['moon_tg_yr'][split][t]['tg_yr']['earth'] + seas_tg[t]['earth']
                     for case, split in zip(('low', 'central', 'high'), SPLITS)} for t in TRAITS}
    absorbed = resort_absorbed(st301)
    out = dict(
        schema=SCHEMA,
        producer=dict(files=producer_files(), inputs=inputs_of(),
                      constants=dict(MOON_SURFACE_GRAVITY=L.MOON_SURFACE_GRAVITY)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        setting=dict(ground_k=st301.ground_k, share_of_the_way=st301.share['resort'], absorbed_sunlight_w_m2=absorbed,
                     methane_share=B.SHARES, communities=dict(earth=B.EARTH, designed=dict(B.EARTH, **B.FAST),
                                                              poor=dict(B.EARTH, **B.POOR)),
                     made_below_g_m3=T.ANOXIC, settled=dict(share=CS.SETTLED, mol_m2=CS.SETTLED_MOL_M2),
                     acceptors_lakes=supply, anaerobic_sets=AOM_SETS, overturn_window_days=T.OVERTURN_AFTER_DAYS),
        reproduction=reproduction,
        targets=targets,
        kinds=dict(index=kinds['index'], index_sha256=kinds['index_sha256'], ladder=ladder, meromictic=mero,
                   rungs=rungs, flux_by_rung=by_rung, lake_area_km2=area, moon_tg_yr=lakes_tg, ventilation=vent,
                   uniform_scale_up_tg_yr=uniform_scale_up(by_rung, area),
                   meromictic_t_km2_yr=meromictic_flux(mero), giants=giants,
                   files={name: c['files'] for name, c in kinds['cases'].items()}),
        seas=seas,
        forms=dict(index=forms['index'], index_sha256=forms['index_sha256'], cases=forms['cases'], ratios=ratios,
                   shih_flux_by_rung=shih_by_rung, shih_moon_tg_yr=lakes_tg_shih, shih_ventilation=vent_shih),
        air=dict(lakes_and_seas=air_by_community(absorbed, both),
                 lakes_and_seas_share=air_for(absorbed, both_share), lakes_and_seas_share_tg_yr=both_share,
                 lakes_and_seas_shih=air_by_community(absorbed, both_shih),
                 lakes=air_by_community(absorbed, {t: lakes_tg[t]['tg_yr'] for t in TRAITS}),
                 lakes_and_seas_aom=air_for(absorbed, both_aom),
                 lakes_and_seas_ventilation=air_for(absorbed, both_vent), lakes_and_seas_ventilation_tg_yr=both_vent,
                 lakes_and_seas_tg_yr=both, lakes_and_seas_shih_tg_yr=both_shih, lakes_and_seas_aom_tg_yr=both_aom,
                 seas_tg_yr=seas_tg, seas_shih_tg_yr=seas_tg_shih),
        biology=biology_tables(ladder, kinds),
        seconds=round(time.time() - tic, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(rnd(lean(out))) + '\n')
    print('wrote', OUT.relative_to(ROOT), f'in {time.time() - tic:.0f} s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
