"""A run's settled climate estimated at a nearby global temperature by pattern scaling, written as data products beside
the run's own.

    climate/gcm/.venv/bin/python -m climate.gcm.pattern_scaling BASE:45-59 --target RUN:35-49 \\
        --pair WARM:45-59/COOL:4-8 --pair WARM:3-8/COOL:35-49

Pattern scaling (Santer et al. 1990; Mitchell 2003; Tebaldi and Arblaster 2014) carries a model's climate to another
global mean temperature by a pattern of change per kelvin of global warming, measured on runs that differ in their
forcing alone. Here the pattern is uniform over land and over open water: each field's mean change per kelvin of the
global surface temperature over the land that both runs of a pair hold and over the water both hold, from pairs of runs
of the same physics at two sunlights. Where a field's maps of change are mostly the noise of short windows, as on the
Moon's T21 grid, a uniform pattern is what the pairs support. A field changes by the mean of the pairs' values where
they agree in sign and is carried unchanged where they disagree, its response unresolved. Amounts (rain, evaporation,
soil water, the column's vapour, snow, the wind's friction velocity and speed) change by a fraction per kelvin,
temperatures and covers by an amount per kelvin, covers kept within 0 and 1. The surface temperature's two values are
scaled together so that the global mean moves by exactly the step asked; the rest keep their measured values.

The base window's climatology (climate.gcm.climatology) and surface wind (climate.gcm.surface_wind) are scaled by the
step from its global mean surface temperature to the target window's, and written as climatology_<run>_y<a>-<b>_<label>
.npz and surface_wind_<run>_y<a>-<b>_<label>.npz in products/, in their own schemas with the run's metadata and an
estimate block, so that every reader of a climatology reads them unchanged. The record, results/gcm/pattern_scaling_
<run>_y<a>-<b>_<label>.json, gives each pair's means and responses and the values applied. A window without a product of
its own is read from the run's year files and written as climatology_<run>_y<a>-<b>.npz or surface_wind_<run>_y<a>-<b>
.npz first.

Santer, B. D., Wigley, T. M. L., Schlesinger, M. E. & Mitchell, J. F. B. (1990). Developing climate scenarios from
equilibrium GCM results. Max-Planck-Institut fur Meteorologie Report 47. Mitchell, T. D. (2003). Pattern scaling: an
examination of the accuracy of the technique for describing future climates. Climatic Change 60, 217-242. Tebaldi, C.
& Arblaster, J. M. (2014). Pattern scaling: its strengths and limitations, and an update on the latest model
simulations. Climatic Change 122, 459-471.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PRODUCTS = HERE / 'products'
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-pattern-scaling/1'

# How each field of a climatology and of a surface wind product moves with the global temperature: 'ratio' by a
# fraction per kelvin, 'add' by an amount per kelvin, 'cover' as 'add' kept within 0 and 1. Fields not named (the
# grid, the land mask, the hour angles, the layers, the daylit counts) are carried unchanged.
CLIMATOLOGY_KINDS = {
    'ts_k': 'add', 'tas_land_sun': 'add', 'tas_land_sun_near': 'add', 'tas_land_sun_far': 'add',
    'pr_mm_day': 'ratio', 'evap_mm_day': 'ratio', 'mrso_m': 'ratio', 'prw_kg_m2': 'ratio', 'snd_m': 'ratio',
    'clt_mean': 'cover', 'low_mean': 'cover', 'high_mean': 'cover', 'sic': 'cover', 'clt_sun': 'cover',
    'low_sun': 'cover', 'high_sun': 'cover', 'cl_layer_mean': 'cover', 'ua_layer_mean': 'add'}
SURFACE_WIND_KINDS = {'ustar_day_ms': 'ratio', 'wind10_day_ms': 'ratio', 'ustar_all_ms': 'ratio',
                      'wind10_all_ms': 'ratio', 'sic': 'cover'}
TEMPERATURE = 'ts_k'
TINY = 1e-12
EVIDENCE = ('An estimate, not a run: the base window\'s fields moved by a uniform pattern of change per kelvin of global '
            'surface temperature, over land and over open water, measured on pairs of runs of the same physics at two '
            'sunlights (climate.gcm.pattern_scaling). The pattern carries the land\'s and the water\'s mean responses; '
            'where they fall within each is the base run\'s own. ')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def weights(lat, lon):
    """Each cell's area weight on a regular-in-longitude grid: the cosine of its latitude."""
    return np.cos(np.radians(np.asarray(lat, float)))[:, None] * np.ones((1, np.asarray(lon).size))


def global_mean(field, lat, lon) -> float:
    w = weights(lat, lon)
    return float((np.asarray(field, float) * w).sum() / w.sum())


def on_grid(field, lat, lon) -> bool:
    field = np.asarray(field)
    return field.ndim >= 2 and field.shape[-2:] == (np.asarray(lat).size, np.asarray(lon).size)


def mean_over(field, weight, mask=None) -> float:
    """The weighted mean of a field over a mask (both broadcast to its shape), its NaN left out."""
    f = np.asarray(field, float)
    w = np.broadcast_to(np.asarray(weight, float) * (1.0 if mask is None else np.asarray(mask, float)), f.shape)
    w = np.where(np.isfinite(f), w, 0.0)
    total = w.sum()
    return float((np.where(w > 0, f, 0.0) * w).sum() / total) if total > 0 else math.nan


def domain_means(field, lat, lon, land, water):
    """A field's means over the land and the water when it lies on the grid, or over all of it otherwise (a composite
    by latitude and hour, cos-latitude weighted, or a profile by layer)."""
    f = np.asarray(field, float)
    if on_grid(f, lat, lon):
        w = weights(lat, lon)
        return {'land': mean_over(f, w, land), 'water': mean_over(f, w, water)}
    if f.ndim == 2 and f.shape[0] == np.asarray(lat).size:
        return {'all': mean_over(f, np.cos(np.radians(np.asarray(lat, float)))[:, None])}
    return {'all': mean_over(f, np.ones(f.shape))}


def pair_response(warm, cool, kinds, delta=None):
    """Each field's change per kelvin of global surface temperature between two windows of runs on the same grid, over
    the land both hold and the water both hold: a fraction per kelvin for 'ratio' fields, an amount for the others.
    delta is the warm window's global mean surface temperature less the cool one's, from ts_k unless given."""
    lat, lon = warm['lat'], warm['lon']
    if delta is None:
        delta = global_mean(warm[TEMPERATURE], lat, lon) - global_mean(cool[TEMPERATURE], lat, lon)
    if abs(delta) < 0.1:
        raise ValueError(f'the pair differs by {delta:.3f} K in global surface temperature, too little to scale by')
    land = (np.asarray(warm['lsm']) > 0.5) & (np.asarray(cool['lsm']) > 0.5)
    water = (np.asarray(warm['lsm']) < 0.5) & (np.asarray(cool['lsm']) < 0.5)
    fields = {}
    for name, kind in kinds.items():
        if name not in warm or name not in cool:
            continue
        mw, mc = domain_means(warm[name], lat, lon, land, water), domain_means(cool[name], lat, lon, land, water)
        fields[name] = {}
        for d in mw:
            if kind == 'ratio':
                value = (mw[d] / mc[d] - 1.0) / delta if abs(mc[d]) > TINY and math.isfinite(mw[d]) else math.nan
            else:
                value = (mw[d] - mc[d]) / delta
            fields[name][d] = dict(warm=mw[d], cool=mc[d], per_K=value)
    return dict(delta_K=delta, land_share=float((weights(lat, lon) * land).sum() / weights(lat, lon).sum()),
                water_share=float((weights(lat, lon) * water).sum() / weights(lat, lon).sum()), fields=fields)


def combine(pairs):
    """The response applied to each field and domain: the mean of the pairs' values where all are finite and agree in
    sign, zero (the field carried unchanged) where they do not."""
    out = {}
    names = set.intersection(*(set(p['fields']) for p in pairs))
    for name in sorted(names):
        out[name] = {}
        for d in pairs[0]['fields'][name]:
            values = [p['fields'][name][d]['per_K'] for p in pairs]
            finite = all(math.isfinite(v) for v in values)
            agree = finite and (all(v > 0 for v in values) or all(v < 0 for v in values))
            out[name][d] = dict(per_K=float(np.mean(values)) if agree else 0.0, pairs=values, resolved=bool(agree))
    return out


def normalise_temperature(applied, fields, name=TEMPERATURE):
    """The surface temperature's land and water values scaled together so that, on the base window's land and water,
    the global mean moves by one kelvin per kelvin; returns the factor."""
    lat, lon = fields['lat'], fields['lon']
    w = weights(lat, lon)
    land = np.asarray(fields['lsm']) > 0.5
    share = float((w * land).sum() / w.sum())
    r = applied[name]
    if not (r['land']['resolved'] and r['water']['resolved']):
        raise ValueError('the pairs do not resolve the surface temperature\'s response')
    moved = share * r['land']['per_K'] + (1.0 - share) * r['water']['per_K']
    factor = 1.0 / moved
    for d in ('land', 'water'):
        r[d]['measured_per_K'] = r[d]['per_K']
        r[d]['per_K'] = r[d]['per_K'] * factor
    return factor


def apply(fields, kinds, applied, delta):
    """The fields moved by delta kelvin of global surface temperature: on the grid each land cell (the base window's mask
    above one half) by the land's response and each water cell by the water's, elsewhere by the field's one response."""
    lat, lon = fields['lat'], fields['lon']
    land = np.asarray(fields['lsm']) > 0.5
    out = dict(fields)
    for name, kind in kinds.items():
        if name not in fields or name not in applied:
            continue
        f = np.asarray(fields[name], float)
        r = applied[name]
        k = r['all']['per_K'] if 'all' in r else np.where(land, r['land']['per_K'], r['water']['per_K'])
        moved = f * np.maximum(0.0, 1.0 + k * delta) if kind == 'ratio' else f + k * delta
        if kind == 'cover':
            moved = np.where(np.isfinite(moved), np.clip(moved, 0.0, 1.0), moved)
        out[name] = moved
    return out


# ---------------------------------------------------------------- windows and products

def parse_window(text):
    """'RUN:a-b' as (run, a, b)."""
    run, span = text.rsplit(':', 1)
    first, last = (int(v) for v in span.split('-'))
    return run, first, last


def window_path(prefix, run, first, last):
    return PRODUCTS / f'{prefix}_{run}_y{first}-{last}.npz'


def read_product(path, run, first, last):
    with np.load(path) as z:
        meta = json.loads(str(z['metadata']))
        if meta.get('run') != run or list(meta.get('years')) != [first, last] or 'estimate' in meta:
            raise SystemExit(f'{path} is not {run} years {first}-{last} as the run gave them')
        return {k: np.asarray(z[k], float) for k in z.files if k != 'metadata'}, meta


def window(prefix, run, first, last, log=print):
    """A window's climatology ('climatology') or surface wind ('surface_wind') product, written from the run's year
    files first when there is none; returns its fields, metadata and path."""
    path = window_path(prefix, run, first, last)
    if not path.exists():
        from climate.gcm import climatology as cl, surface_wind as sw
        from climate.gcm.exoplasim_run import RUNS
        log(f'reading {run} years {first}-{last} for its {prefix.replace("_", " ")}')
        progress = json.loads((RUNS / run / 'progress.json').read_text())
        if prefix == 'climatology':
            data = cl.climatology(run, first, last)
            meta = dict(schema=cl.SCHEMA, run=run, years=[first, last], configuration=progress['configuration'],
                        deck_sigma=cl.DECK_SIGMA, producer=dict(domain='climate', files={'gcm/climatology.py': sha256(cl.__file__)[:16]}),
                        evidence=cl.EVIDENCE, reading_rule=cl.READING_RULE)
        else:
            data, gascon = sw.surface_wind(run, first, last)
            meta = dict(schema=sw.SCHEMA, run=run, years=[first, last], configuration=progress['configuration'],
                        gas_constant_J_kg_K=gascon, day_cosine=sw.DAY_COSINE,
                        producer=dict(domain='climate', files={'gcm/surface_wind.py': sha256(sw.__file__)[:16]}),
                        evidence=sw.EVIDENCE, reading_rule=sw.READING_RULE)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, metadata=json.dumps(meta), **{k: np.asarray(v).astype(np.float32) for k, v in data.items()})
    fields, meta = read_product(path, run, first, last)
    return fields, meta, path


def label_of(target_k):
    return f'scaled{target_k:.2f}K'


def write_estimate(path, fields, meta, estimate):
    meta = dict(meta)
    meta['estimate'] = estimate
    meta['evidence'] = EVIDENCE + meta.get('evidence', '')
    np.savez_compressed(path, metadata=json.dumps(meta), **{k: np.asarray(v).astype(np.float32) for k, v in fields.items()})
    return path


def rounded(x, digits=6):
    if isinstance(x, dict):
        return {k: rounded(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rounded(v, digits) for v in x]
    if isinstance(x, float):
        return None if not math.isfinite(x) else round(x, digits)
    return x


def build(base, target, pairs, log=print):
    """The estimate: base, target and each pair's windows as (run, first, last); returns the record."""
    run, first, last = base
    clim, clim_meta, clim_path = window('climatology', *base, log=log)
    wind, wind_meta, wind_path = window('surface_wind', *base, log=log)
    goal, _, goal_path = window('climatology', *target, log=log)
    base_k = global_mean(clim[TEMPERATURE], clim['lat'], clim['lon'])
    target_k = global_mean(goal[TEMPERATURE], goal['lat'], goal['lon'])
    delta = target_k - base_k
    measured, measured_wind, used = [], [], []
    sunlight = lambda m: m['configuration']['model'].get('sunlight_scale')                 # noqa: E731
    for warm_w, cool_w in pairs:
        warm, warm_meta, warm_path = window('climatology', *warm_w, log=log)
        cool, cool_meta, cool_path = window('climatology', *cool_w, log=log)
        response = pair_response(warm, cool, CLIMATOLOGY_KINDS)
        warm_s, _, warm_s_path = window('surface_wind', *warm_w, log=log)
        cool_s, _, cool_s_path = window('surface_wind', *cool_w, log=log)
        response_wind = pair_response(warm_s, cool_s, SURFACE_WIND_KINDS, delta=response['delta_K'])
        measured.append(response)
        measured_wind.append(response_wind)
        used.append(dict(warm=dict(run=warm_w[0], years=list(warm_w[1:]), sunlight_scale=sunlight(warm_meta),
                                   global_surface_temperature_K=global_mean(warm[TEMPERATURE], warm['lat'], warm['lon'])),
                         cool=dict(run=cool_w[0], years=list(cool_w[1:]), sunlight_scale=sunlight(cool_meta),
                                   global_surface_temperature_K=global_mean(cool[TEMPERATURE], cool['lat'], cool['lon'])),
                         products={str(p.name): sha256(p) for p in (warm_path, cool_path, warm_s_path, cool_s_path)}))
        log(f'pair {warm_w[0]}:{warm_w[1]}-{warm_w[2]} / {cool_w[0]}:{cool_w[1]}-{cool_w[2]}: '
            f'{response["delta_K"]:+.3f} K global')
    applied = combine(measured)
    factor = normalise_temperature(applied, clim)
    applied_wind = combine(measured_wind)
    label = label_of(target_k)
    estimate = dict(method='uniform pattern scaling over land and open water (climate.gcm.pattern_scaling)', label=label,
                    base=dict(run=run, years=[first, last], global_surface_temperature_K=base_k),
                    target=dict(run=target[0], years=list(target[1:]), global_surface_temperature_K=target_k),
                    delta_K=delta, pairs=[dict(warm=u['warm'], cool=u['cool']) for u in used],
                    record=str(record_path(run, first, last, label).relative_to(HERE.parent.parent)))
    out_clim = PRODUCTS / f'climatology_{run}_y{first}-{last}_{label}.npz'
    out_wind = PRODUCTS / f'surface_wind_{run}_y{first}-{last}_{label}.npz'
    scaled = apply(clim, CLIMATOLOGY_KINDS, applied, delta)
    write_estimate(out_clim, scaled, clim_meta, dict(estimate, responses=rounded(
        {k: {d: v['per_K'] for d, v in r.items()} for k, r in applied.items()})))
    scaled_wind = apply(wind, SURFACE_WIND_KINDS, applied_wind, delta)
    write_estimate(out_wind, scaled_wind, wind_meta, dict(estimate, responses=rounded(
        {k: {d: v['per_K'] for d, v in r.items()} for k, r in applied_wind.items()})))
    lat, lon = clim['lat'], clim['lon']
    w = weights(lat, lon)
    land = clim['lsm'] > 0.5
    over = lambda f, m: float((np.asarray(f, float) * w * m).sum() / (w * m).sum())   # noqa: E731
    summary = {}
    for name in ('ts_k', 'pr_mm_day', 'evap_mm_day', 'mrso_m', 'prw_kg_m2', 'clt_mean'):
        summary[name] = dict(land=[over(clim[name], land), over(scaled[name], land)],
                             water=[over(clim[name], ~land), over(scaled[name], ~land)])
    pe = lambda f: f['pr_mm_day'] - f['evap_mm_day']                                     # noqa: E731
    summary['rain_less_evaporation_mm_day'] = dict(land=[over(pe(clim), land), over(pe(scaled), land)])
    summary['global_surface_temperature_K'] = [base_k, global_mean(scaled[TEMPERATURE], lat, lon)]
    return dict(
        schema=SCHEMA, producer=dict(domain='climate', files={'gcm/pattern_scaling.py': sha256(__file__)[:16]}),
        evidence=('An estimate by pattern scaling, not a run: the base window\'s fields moved by a uniform pattern of '
                  'change per kelvin of global surface temperature over land and over open water, measured on pairs of '
                  'runs of the same physics at two sunlights. A pair\'s window still warming or cooling gives a '
                  'transient response, which the land follows within a year or two of the global temperature. Not '
                  'validated against a run at the target.'),
        reading_rule=('Temperatures in K, rain and evaporation in mm/day, soil water in m, vapour in kg/m2, covers as '
                      'fractions, winds in m/s. Each pair gives each field\'s mean over the land both runs hold and the '
                      'water both hold (or over all of it for the composites and profiles), warm and cool, and per_K: '
                      'a fraction per kelvin for the ratio fields, an amount per kelvin for the others. applied is '
                      'what the estimate moved each field by per kelvin, zero where the pairs disagree in sign '
                      '(resolved false); the surface temperature\'s values are scaled by temperature_factor so the '
                      'global mean moves by delta_K exactly. summary gives [base, estimate] means over the base '
                      'window\'s land and water.'),
        estimate={k: v for k, v in estimate.items() if k != 'record'},
        temperature_factor=factor,
        pairs=[dict(u, delta_K=m['delta_K'], land_share=m['land_share'], water_share=m['water_share'],
                    climatology=rounded(m['fields']), surface_wind=rounded(mw['fields']))
               for u, m, mw in zip(used, measured, measured_wind)],
        applied=dict(climatology=rounded(applied), surface_wind=rounded(applied_wind)),
        summary=rounded(summary),
        products={str(p.relative_to(HERE.parent.parent)): sha256(p) for p in (out_clim, out_wind)},
        inputs={str(p.relative_to(HERE.parent.parent)): sha256(p) for p in (clim_path, wind_path, goal_path)})


def record_path(run, first, last, label):
    return RESULTS / f'pattern_scaling_{run}_y{first}-{last}_{label}.json'


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('base', help='the window to scale, RUN:first-last')
    p.add_argument('--target', required=True, help='the window whose global mean surface temperature the estimate holds')
    p.add_argument('--pair', action='append', required=True, help='WARM:first-last/COOL:first-last, two or more')
    args = p.parse_args(argv)
    pairs = [tuple(parse_window(w) for w in text.split('/')) for text in args.pair]
    if len(pairs) < 2:
        raise SystemExit('give at least two pairs, so that a response the pairs disagree on is left unapplied')
    base = parse_window(args.base)
    record = build(base, parse_window(args.target), pairs, log=lambda *a: print(*a, flush=True))
    path = record_path(*base, record['estimate']['label'])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=1) + '\n')
    s = record['summary']
    print(f"global surface temperature {s['global_surface_temperature_K'][0]:.3f} -> {s['global_surface_temperature_K'][1]:.3f} K "
          f"(step {record['estimate']['delta_K']:+.3f} K, temperature factor {record['temperature_factor']:.3f})")
    for name in ('ts_k', 'pr_mm_day', 'evap_mm_day', 'mrso_m', 'rain_less_evaporation_mm_day'):
        print(f"  {name}: land {s[name]['land'][0]:.4f} -> {s[name]['land'][1]:.4f}")
    unresolved = [f'{k}.{d}' for k, r in {**record['applied']['climatology'], **record['applied']['surface_wind']}.items()
                  for d, v in r.items() if not v['resolved']]
    print('  unresolved, carried unchanged:', ', '.join(unresolved) or 'none')
    print(f'wrote {path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
