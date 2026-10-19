"""Tracing a saved case, cached: the tracer's passes over a case, their summary per lane, and the cache key.

A cyclic record (a settled cycle, as lake_anoxia's and the kinds' lakes and the seas' columns are) is passed over until
the methane's inventory at its end settles; a record that starts from an oxygenated column (the giants, the meromictic
lakes) runs once, in parts, each summarised and its closing column kept. Each case is cached in
research/runs/lake_methane/cases under a key of the code that traces it, the saved files' contents, the column's
diffusivity form, the lanes and the options. The files enter by their sha256 alone, so a producer that renames its
files (the seas' study re-keys its columns by hard links) leaves the key standing.

The key holds this module, the tracer, the readers, the biology, the lake model and lake_anoxia's gas physics whole,
and of lake_oxygen's and lake_anoxia's runners only what the tracing uses, by name: the source of each function and
the value of each constant, as this process loaded them (read once, at import). An edit elsewhere in those files leaves
the cached cases standing; any change to what the tracing runs re-keys them. run.py's analysis is not in the key. A
cache-only run (lake_oxygen's cache_only()) traces no case and stops at the first one not cached.
"""
from __future__ import annotations

import dataclasses
import hashlib
import inspect
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_anoxia import gases as G
from research.studies.lake_anoxia import run as A
from research.studies.lake_methane import biology as B
from research.studies.lake_methane import columns as C
from research.studies.lake_methane import tracer as T
from research.studies.lake_oxygen import mixing_forms as MF
from research.studies.lake_oxygen import run as R

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASES = ROOT / 'research' / 'runs' / 'lake_methane' / 'cases'
SETTLED = 0.005                           # a pass settles when the inventories move less than this share,
SETTLED_MOL_M2 = 1e-3                     # or less than this, mol m-2, from the pass before
WHOLE = (Path(__file__), HERE / 'tracer.py', HERE / 'columns.py', HERE / 'biology.py', Path(L.__file__),
         Path(G.__file__))
USED = {
    R: dict(code=('kd', 'scenarios', 'window_humidity', 'windows', 'resort_windows', 'hourly', 'physics',
                  'forcing_for', 'site_info'),
            values=('DESIGN_PRESSURE_PA', 'DESIGN_CO2_PPM', 'HELD_FRESH', 'TARGET_NAMES', 'STEPS_PER_MEAN',
                    'STEPS_PER_MEAN_DEEP', 'LOWEST_FLOOR_MS', 'SOLAR_HOUR_BINS', 'ANOXIC', 'TROPHIC', 'TRANSFER',
                    'KINETICS', 'SOURCE_C', 'LABILE_DAYS', 'METHANE_SHARE', 'IRON', 'WIND_SHARE', 'DEEP_MIXING',
                    'CAP_WIND', 'DISSIPATION', 'CAP', 'UNSTABLE_DIFFUSIVITY', 'SOLUTE', 'COUPLED')),
    A: dict(code=(),
            values=('HALF_SATURATION', 'SURFACE_ATM', 'BUBBLE_MM', 'BUBBLE_CH4_FRACTION', 'SEDIMENT_DECAY_PER_M',
                    'ANOXIC')),
    MF: dict(code=('shih_bouffard',), values=('NU', 'RE_MOLECULAR', 'RE_ENERGETIC')),
}


def _fingerprint():
    h = hashlib.sha256()
    for p in WHOLE:
        h.update(p.name.encode())
        h.update(C.sha256(p).encode())
    for module, used in USED.items():
        h.update(module.__name__.encode())
        for name in used['code']:
            if hasattr(module, name):
                h.update(name.encode())
                h.update(inspect.getsource(getattr(module, name)).encode())
        for name in used['values']:
            if hasattr(module, name):
                h.update(name.encode())
                h.update(json.dumps(getattr(module, name), sort_keys=True, default=repr).encode())
    return h.hexdigest()


FINGERPRINT = _fingerprint()              # read once, as the modules load


def case_key(col, lanes, extra):
    h = hashlib.sha256()
    h.update(FINGERPRINT.encode())
    files = col.meta.get('files') or {col.meta.get('file'): col.meta.get('sha256')}
    h.update(json.dumps([files[k] for k in sorted(files)]).encode())      # the contents in their passes' order
    h.update(json.dumps(getattr(col, 'form', None) or 'osborn').encode())
    h.update(json.dumps([dataclasses.asdict(l) for l in lanes], sort_keys=True).encode())
    h.update(json.dumps(extra, sort_keys=True, default=str).encode())
    return h.hexdigest()[:20]


def settled_between(prev, now, lanes, rule):
    """Whether two passes' closing inventories (mol m-2, the methane lanes in order) agree: under 'all' every lane to
    SETTLED of its inventory; under 'oxidising' the lanes whose methanotrophs oxidise, to SETTLED or SETTLED_MOL_M2.
    A lane without oxidation keeps gathering methane in anoxic water until its store is full, when all it makes
    leaves: its steady emission is what it makes, and it is read so."""
    prev, now = np.asarray(prev, float), np.asarray(now, float)
    if rule == 'all':
        return bool(np.all(np.abs(now - prev) <= SETTLED * np.maximum(np.abs(prev), 1e-6)))
    keep = np.array([l.aerobic for l in lanes if l.kind == 'methane'], bool)
    if not keep.any():
        keep[:] = True
    tol = np.maximum(SETTLED * np.abs(prev), SETTLED_MOL_M2)
    return bool(np.all((np.abs(now - prev) <= tol)[keep]))


def pass_fluxes(tr, out, col):
    """One pass's emission (across the surface and by bubbles) and methane made, mol per m2 of surface a year, for
    each methane lane in order: how the reported fluxes moved from pass to pass."""
    start, stop = out['steps']
    years = (stop - start) * col.dt / T.YEAR
    bud = out['budget']
    per = lambda k: bud[k][:, tr.is_m].sum(0) / tr.a0 / years
    return dict(emitted_mol_m2_yr=(per('evaded') + per('bubbles_air')).tolist(),
                made_mol_m2_yr=(per('made_water') + per('made_floor')).tolist())


def trace_case(col, lanes, name, cyclic=True, max_passes=4, ml_oxygen=True, clones=False, start_empty=False,
               made_below=T.ANOXIC, oxygen_budget=True, oxygen_rule='proportional', chunks=1, settle_rule='oxidising'):
    """The tracer through a case (module docstring), cached; returns the summary."""
    extra = dict(cyclic=cyclic, max_passes=max_passes, ml_oxygen=ml_oxygen, clones=clones, settled=SETTLED,
                 start_empty=start_empty, made_below=made_below, oxygen_budget=oxygen_budget,
                 oxygen_rule=oxygen_rule, chunks=chunks)
    if settle_rule != 'all':
        extra.update(settle_rule=settle_rule, settled_mol_m2=SETTLED_MOL_M2)
    key = case_key(col, lanes, extra)
    path = CASES / f'{name}_{key}.json'
    if path.exists():
        return json.loads(path.read_text())
    if R.cache_only():
        R.uncached(path)
    tic = time.time()
    probe = T.Tracer(col, lanes)
    bub = None
    if any(l.kind == 'methane' and l.ebullition for l in lanes):
        bub = T.Bubbles(col.basin, col.t.astype(float).mean(0), col.physics.gravity, probe.floor_production_max())
    tr = T.Tracer(col, lanes, bubbles=bub, mixed_layer_oxygen=ml_oxygen, made_below=made_below,
                  oxygen_budget=oxygen_budget, oxygen_rule=oxygen_rule)
    state = tr.initial(0)
    if start_empty:
        for k in ('c', 'model'):
            state[k][:, tr.is_m] = 0.0
    passes, parts = [], []
    m = tr.is_m
    first = int(col.snap_step[0]) + 1
    if cyclic:
        for p in range(max_passes):
            inv0 = tr.inventory(state['c'])
            out = tr.run(state, start=first if p == 0 else 0, profiles=clones)
            inv1 = tr.inventory(state['c'])
            passes.append(dict(start_mol_m2=(inv0[m] / tr.a0).tolist(), end_mol_m2=(inv1[m] / tr.a0).tolist(),
                               **pass_fluxes(tr, out, col)))
            if p >= 1 and settled_between(passes[-2]['end_mol_m2'], passes[-1]['end_mol_m2'], lanes, settle_rule):
                break
        res = summarise(tr, out, col, lanes, clones)
        res['settled'] = bool(len(passes) > 1 and settled_between(passes[-2]['end_mol_m2'],
                                                                   passes[-1]['end_mol_m2'], lanes, settle_rule))
        if settle_rule != 'all':                 # the closing column, for the floor's approach to its bubble point
            res['closing_mol_m3'] = {str(j): state['c'][:, j].tolist() for j in np.flatnonzero(m)}
            res['closing_snapshot'] = int(col.snap_step.size - 1)
    else:
        edges = np.linspace(first, col.steps, chunks + 1).round().astype(int)
        for a_, b_ in zip(edges[:-1], edges[1:]):
            out = tr.run(state, start=int(a_), stop=int(b_), profiles=clones)
            part = summarise(tr, out, col, lanes, clones)
            part['closing_mol_m3'] = {str(j): state['c'][:, j].tolist() for j in np.flatnonzero(m)}
            part['closing_snapshot'] = int(np.searchsorted(col.snap_step, b_ - 1, side='right')) - 1
            parts.append(part)
        res = parts[-1] if chunks == 1 else dict(case=parts[-1]['case'], lanes=parts[-1]['lanes'], parts=parts)
        res['settled'] = None
    res['passes'] = passes
    res['seconds'] = round(time.time() - tic, 1)
    res['bubbles'] = None if bub is None else dict(deepest_bubbling_floor_m=bub.deepest_m,
                                                   survival_at_10m=float(np.interp(10.0, col.basin.centre,
                                                                                   bub.survival)))
    CASES.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res) + '\n')
    print(f'  {name}: {len(lanes)} lanes, {len(passes) or chunks} {"passes" if cyclic else "parts"}, '
          f'{res["seconds"]:.0f} s', flush=True)
    return res


def summarise(tr, out, col, lanes, clones):
    """Per lane, per m2 of lake surface and per year of the last pass: what was made, what left by each road, the
    inventory, each overturn's store and release; per case the overturns, the oxygen they leave, the setting."""
    start, stop = out['steps']
    years = (stop - start) * col.dt / T.YEAR
    a0, v = tr.a0, tr.v
    bud, ev = out['budget'], out['events']
    rows = []
    for j, lane in enumerate(lanes):
        if lane.kind != 'methane':
            continue
        f = lambda k, ph=None: float((bud[k][:, j].sum() if ph is None else bud[k][ph, j]) / a0 / years)
        made = f('made_water') + f('made_floor')
        emitted = [f('evaded', ph) + f('bubbles_air', ph) for ph in (0, 1)]
        ox_b, ev_b = f('oxidised', 0), f('evaded', 0)
        inv = out['snap_inventory'][:, j] / a0
        stock = out['mean_b'][:, j]
        k1 = tr.q20[j] * stock / tr.kc[j] * 3600.0                                   # per hour at 20 C
        events = []
        for e, (a_, b_) in enumerate(tr.events):
            g = lambda k: float(ev['budget'][k][e, j] / a0)
            events.append(dict(stored_mol_m2=float(ev['inventory_before'][e, j] / a0),
                               emitted_mol_m2=g('evaded') + g('bubbles_air'), oxidised_mol_m2=g('oxidised'),
                               aom_mol_m2=g('aom_water'), made_mol_m2=g('made_water') + g('made_floor')))
        rows.append(dict(
            label=lane.label, column=col.names[lane.column], share=lane.share,
            made_anoxic_mol_m2_yr=made, floors_mol_m2_yr=f('floor_total'), bubbles_mol_m2_yr=f('bubbles'),
            emitted_mol_m2_yr=emitted[0] + emitted[1], emitted_between_mol_m2_yr=emitted[0],
            emitted_at_overturns_mol_m2_yr=emitted[1], diffusive_mol_m2_yr=f('evaded'),
            bubbles_to_air_mol_m2_yr=f('bubbles_air'), invaded_mol_m2_yr=f('invaded'),
            oxidised_mol_m2_yr=f('oxidised'),
            oxidised_between_mol_m2_yr=ox_b, oxidised_at_overturns_mol_m2_yr=f('oxidised', 1),
            floor_surface_oxidised_mol_m2_yr=f('floor_surface_oxidised'), aom_floor_mol_m2_yr=f('aom_floor'),
            aom_water_mol_m2_yr=f('aom_water'),
            oxycline_oxidised_share=ox_b / (ox_b + ev_b) if ox_b + ev_b > 0 else None,
            o2_used_g_m2_yr=f('o2_used_g'), o2_model_methane_g_m2_yr=f('o2_model_g'),
            inventory_mol_m2=dict(mean=float(inv.mean()), max=float(inv.max()), end=float(inv[-1])),
            stock_max_um_c=float(stock.max() * 1e3), first_order_max_per_h=float(k1.max()),
            stock_cells_ml_max=float(stock.max() * 1e3 / B.UM_C_PER_CELL_ML),
            events=events))
    # the case: its overturns and the oxygen the overturned lake holds (lake_anoxia's reading: the column just after
    # each event begins, whole-lake mean)
    events = []
    for e, (a_, b_) in enumerate(tr.events):
        k = min(int(np.searchsorted(col.snap_step, a_)) + 1, col.snap_step.size - 1)
        held = [float((col.oxygen[k, :, c].astype(float) * v).sum() / a0) for c in range(len(col.names))]
        events.append(dict(first_step=a_, last_step=b_, year=a_ * col.dt / T.YEAR, o2_held_g_m2=held))
    case = dict(name=col.name, years=years, steps=stop - start, depth_m=col.basin.depth, area_km2=a0 / 1e6,
                mean_depth_m=float(v.sum() / a0), layers=col.basin.n, columns=col.names,
                overturns=events, bottom_t_c=float(col.t[:, -1].astype(float).mean()),
                surface_t_c=float(col.surface_t.mean()),
                mixed_depth_m=[float(np.quantile(col.mixed_depth, q)) for q in (0.5, 0.9, 0.99)],
                transfer_velocity_m_day=float(col.k_o2.mean() * T.DAY * tr.k_ratio))
    res = dict(case=case, lanes=rows)
    if clones:
        res['clones'] = clone_checks(tr, out, col, lanes)
    return res


def clone_checks(tr, out, col, lanes):
    """The age and debt clones against the model's own saved age and debt."""
    snap, steps = out['snap_c'], out['snap_steps']
    idx = np.searchsorted(col.snap_step, steps)
    v, a0 = tr.v, tr.a0
    res = {}
    for j, lane in enumerate(lanes):
        if lane.kind == 'age':
            model = col.age[idx].astype(float)
            err = np.abs(snap[:, :, j] - model)
            res['age'] = dict(mean_abs_days=float(err.mean()), p999_days=float(np.quantile(err, 0.999)),
                              max_days=float(err.max()),
                              lake_mean_max_days=float(np.abs((snap[:, :, j] - model) @ v / v.sum()).max()),
                              lake_mean_days=float((model @ v).mean() / v.sum()))
        if lane.kind == 'debt':
            model = (col.debt[idx][:, :, lane.column].astype(float) * lane.share / T.CH4_TO_O2) @ v / a0
            clone = snap[:, :, j] @ v / a0
            res[f'debt/{col.names[lane.column]}'] = dict(model_mean_mol_m2=float(model.mean()),
                                                         clone_mean_mol_m2=float(clone.mean()),
                                                         model_max_mol_m2=float(model.max()),
                                                         clone_max_mol_m2=float(clone.max()))
    return res


