"""A run's settled climate estimated on a changed surface from runs already made, written as data products beside the
run's own.

    climate/gcm/.venv/bin/python -m climate.gcm.surface_change BASE:45-59:scaled300.31K --mask <mask record> \\
        --pair LAND_HERE:35-49/WATER_HERE:45-59:scaled300.31K --tag resurfaced

When the water on the Moon changes, so does the surface the climate model runs on. The model holds each T21 cell as
all land or all water (the cells with the least land become sea until the sea covers the surface's water share), so
two settled windows of the same physics at one global temperature whose surfaces differ show what a change does where
it happens, in the cells that are land in one and water in the other. On the Moon's month-long day the land heats and
the rain falls over it: such a cell rains about a third as much as water as it does as land, and as water it rains
much the same whatever it rained as land, while it evaporates at about the open water's rate for its latitude.

Lakes form where it rains, so carrying a climate onto a new surface cell by whole cell sets the lakes and the rain
chasing each other: a cell that turns to water loses its rain, its lake shrinks, the cell turns back to land, and the
rain returns. The estimate therefore moves each cell by the change of its water fraction, from the surface the run
stood on to the new one, which varies smoothly with the lakes:

- a cell's rain changes by the change of its water fraction times the difference between its rain as water and as
  land. For a cell the run held as land, its rain as land is the run's and its rain as water the line a + b x (its
  rain as land), fitted over the pairs' switched cells in both directions; for a cell the run held as water, its
  rain as water is the run's and its rain as land the mean of the run's land cells beside it (of its latitude row,
  or of all the land, where it has none beside it);
- the evaporation of a cell the run held as land changes by the change of its water fraction times the difference
  between the open water's mean evaporation at its latitude and its own; a cell the run held as water keeps its
  evaporation, which also keeps the open water's rate by latitude as the run gave it;
- the rain that the new water no longer takes falls on the land. In every cell the share that is land has its rain
  as land raised by one factor, and in the cells the run held as land the evaporation follows that rain along a
  curve of Turc's and Pike's form, E = P / (1 + (P / Em)^n)^(1/n), whose exponent n is fitted on the base window's
  own land cells and whose ceiling Em is the cell's own, so that dry land evaporates nearly all of the added rain
  and wet land sheds it. The factor is the one at which the whole surface's rain less evaporation is what it was in
  the base window: the model's own closure of its water;
- every other field (temperatures, clouds, soil water, winds) and the land mask are the base window's. On the
  surface the run stood on the estimate is the run's own climate.

Where the displaced rain falls is the estimate's assumption, not a measurement: the pairs' swaps nearly cancel, so
they do not show it. The rain over water in the runs is nearly the same whatever surrounds it, which leaves the land.
A uniform factor over all the land is the plainest placement; rain displaced from a large new lake would in truth
fall mostly on the land round it. Moving a cell in proportion to its water fraction is an interpolation between the
two states the model can hold it in, which no run tests. The pairs' masks are compared cell by cell on the model's
grid in each run's own frame, which may differ by a degree where a cell spans five and a half.

The products are climatology_<run>_y<a>-<b>_<label>_<tag>.npz and surface_wind_<run>_y<a>-<b>_<label>_<tag>.npz in
products/ (the surface wind is the base window's, carried under the new label), in their own schemas with the run's
metadata and an estimate block, so that every reader of a climatology reads them unchanged. The record is
results/gcm/surface_change_<run>_y<a>-<b>_<label>_<tag>.json.

Turc, L. (1954). Le bilan d'eau des sols: relation entre les precipitations, l'evaporation et l'ecoulement. Annales
Agronomiques 5, 491-569. Pike, J. G. (1964). The estimation of annual run-off from meteorological data in a tropical
climate. Journal of Hydrology 2, 116-123.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from climate.gcm import boundary

HERE = Path(__file__).resolve().parent
PRODUCTS = HERE / 'products'
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-surface-change/1'
MASK_SCHEMA = 'terluna.research.atlas-round-design-mask/1'
RAIN, EVAPORATION, MASK = 'pr_mm_day', 'evap_mm_day', 'lsm'
EXPONENTS = np.arange(1.0, 8.01, 0.05)      # the curve's exponents tried
CEILINGS = np.arange(1.0, 8.01, 0.05)       # and its ceilings, mm/day, for the fit of the exponent
FIT_LATITUDE = 45.0                         # the exponent is fitted on land within this latitude, where one ceiling serves
EVIDENCE = ('An estimate, not a run: the base window\'s rain and evaporation carried onto a changed surface, each cell '
            'moved by the change of its water fraction times the contrast between land and water measured where '
            'settled runs of the same physics hold different masks, the displaced rain put on the land and the '
            'surface\'s water budget closed as the model closes it (climate.gcm.surface_change). Temperatures, clouds, '
            'soil water, winds and the land mask are the base window\'s. ')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def weights(lat, lon):
    """Each cell's area weight on a regular-in-longitude grid: the cosine of its latitude."""
    return np.cos(np.radians(np.asarray(lat, float)))[:, None] * np.ones((1, np.asarray(lon).size))


def parse_window(text):
    """'RUN:a-b' or 'RUN:a-b:label' as (run, a, b, label or None)."""
    parts = text.split(':')
    if len(parts) not in (2, 3):
        raise SystemExit(f'{text}: a window is RUN:first-last or RUN:first-last:label')
    first, last = (int(v) for v in parts[1].split('-'))
    return parts[0], first, last, (parts[2] if len(parts) == 3 else None)


def product_path(prefix, run, first, last, label=None):
    return PRODUCTS / (f'{prefix}_{run}_y{first}-{last}' + (f'_{label}' if label else '') + '.npz')


def read_product(prefix, run, first, last, label=None):
    """A window's product, the run's own or an estimate under a label; returns its fields, metadata and path."""
    path = product_path(prefix, run, first, last, label)
    with np.load(path) as z:
        meta = json.loads(str(z['metadata']))
        if meta.get('run') != run or list(meta.get('years')) != [first, last] or \
                meta.get('estimate', {}).get('label') != label:
            raise SystemExit(f'{path} is not {run} years {first}-{last}' + (f' estimated as {label}' if label else ''))
        return {k: np.asarray(z[k], float) for k in z.files if k != 'metadata'}, meta, path


# ---------------------------------------------------------------- the contrast between land and water

def switched(a, b):
    """The cells that are land in one window and water in the other, from both directions: each cell's rain and
    evaporation as land and as water, its row and its area weight."""
    la, lb = np.asarray(a[MASK]) > 0.5, np.asarray(b[MASK]) > 0.5
    w = weights(a['lat'], a['lon'])
    rows = np.arange(la.shape[0])[:, None] * np.ones(la.shape[1], int)
    out = {k: [] for k in ('rain_land', 'rain_water', 'evaporation_land', 'evaporation_water', 'row', 'weight')}
    for land, water, m in ((a, b, la & ~lb), (b, a, lb & ~la)):
        out['rain_land'].append(land[RAIN][m])
        out['rain_water'].append(water[RAIN][m])
        out['evaporation_land'].append(land[EVAPORATION][m])
        out['evaporation_water'].append(water[EVAPORATION][m])
        out['row'].append(rows[m])
        out['weight'].append(w[m])
    return {k: np.concatenate(v) for k, v in out.items()}


def water_rain_fit(cells):
    """The line rain as water = a + b x rain as land over switched cells, by area-weighted least squares, with the
    spread about it and the ratio of the means."""
    x, y, w = cells['rain_land'], cells['rain_water'], cells['weight']
    xm, ym = np.average(x, weights=w), np.average(y, weights=w)
    b = float(np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2))
    a = float(ym - b * xm)
    spread = float(np.sqrt(np.average((y - a - b * x) ** 2, weights=w)))
    return dict(a_mm_day=a, b=b, spread_mm_day=spread, cells=int(x.size), rain_as_land_mm_day=float(xm),
                rain_as_water_mm_day=float(ym), ratio_of_means=float(ym / xm),
                evaporation_as_land_mm_day=float(np.average(cells['evaporation_land'], weights=w)),
                evaporation_as_water_mm_day=float(np.average(cells['evaporation_water'], weights=w)))


def zonal_water_evaporation(fields):
    """The open water's mean evaporation by latitude row, rows without water taking their neighbours' by latitude."""
    water = np.asarray(fields[MASK]) < 0.5
    ev = np.where(water, fields[EVAPORATION], np.nan)
    has = water.any(axis=1)
    zonal = np.full(water.shape[0], np.nan)
    zonal[has] = np.nanmean(ev[has], axis=1)
    lat = np.asarray(fields['lat'], float)
    order = np.argsort(lat[has])
    return np.interp(lat, lat[has][order], zonal[has][order])


# ---------------------------------------------------------------- the land's evaporation on its rain

def curve(rain, ceiling, n):
    """Turc's and Pike's form: the evaporation of land that rains `rain` under a ceiling `ceiling`."""
    rain = np.asarray(rain, float)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        return np.where(np.isfinite(ceiling), rain / (1.0 + (rain / ceiling) ** n) ** (1.0 / n), rain)


def fit_exponent(fields, within=FIT_LATITUDE):
    """The curve's exponent, with one ceiling, that best gives the evaporation of the window's land cells within
    `within` degrees of the equator from their rain (area-weighted least squares over a grid of both)."""
    land = np.asarray(fields[MASK]) > 0.5
    lat = np.asarray(fields['lat'], float)[:, None] * np.ones((1, land.shape[1]))
    m = land & (np.abs(lat) <= within) & (fields[RAIN] > 0)
    p, e, w = fields[RAIN][m], fields[EVAPORATION][m], weights(fields['lat'], fields['lon'])[m]
    best = None
    for n in EXPONENTS:
        for ceiling in CEILINGS:
            miss = float(np.average((curve(p, ceiling, n) - e) ** 2, weights=w))
            if best is None or miss < best[0]:
                best = (miss, float(n), float(ceiling))
    return dict(n=best[1], ceiling_mm_day=best[2], spread_mm_day=math.sqrt(best[0]), cells=int(m.sum()))


def ceilings(rain, evaporation, n):
    """Each land cell's own ceiling, the one that puts its rain and evaporation on the curve; infinite where the cell
    evaporates all its rain or more, which then evaporates the same share of whatever it rains."""
    rain, evaporation = np.asarray(rain, float), np.asarray(evaporation, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = rain / evaporation
        return np.where((evaporation > 0) & (ratio > 1.0 + 1e-9), rain / np.maximum(ratio ** n - 1.0, 1e-300) ** (1.0 / n), np.inf)


def land_evaporation(rain_new, rain, evaporation, n):
    """The land's evaporation on a new rain: along each cell's own curve, or the same share of the rain where the
    cell has no ceiling."""
    em = ceilings(rain, evaporation, n)
    with np.errstate(divide='ignore', invalid='ignore'):
        share = np.where(rain > 0, evaporation / rain, 1.0)
    return np.where(np.isfinite(em), curve(rain_new, em, n), np.asarray(rain_new, float) * share)


# ---------------------------------------------------------------- the estimate

def beside(mask):
    """The cells with a marked cell among their eight neighbours (longitude wraps), and how many."""
    count = np.zeros(mask.shape, int)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                moved = np.roll(mask, dj, axis=1)
                if di == 1:
                    moved = np.vstack([np.zeros((1, mask.shape[1]), bool), moved[:-1]])
                elif di == -1:
                    moved = np.vstack([moved[1:], np.zeros((1, mask.shape[1]), bool)])
                count += moved
    return count


def neighbours_mean(field, mask):
    """Each cell's mean of a field over its marked neighbours, NaN where it has none."""
    total = np.zeros(mask.shape)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                m = np.roll(mask, dj, axis=1)
                f = np.roll(np.where(mask, field, 0.0), dj, axis=1)
                if di == 1:
                    m, f = np.vstack([np.zeros((1, mask.shape[1]), bool), m[:-1]]), np.vstack([np.zeros((1, mask.shape[1])), f[:-1]])
                elif di == -1:
                    m, f = np.vstack([m[1:], np.zeros((1, mask.shape[1]), bool)]), np.vstack([f[1:], np.zeros((1, mask.shape[1]))])
                total += f
    n = beside(mask)
    return np.where(n > 0, total / np.maximum(n, 1), np.nan)


def carry(fields, wet_run, wet_new, fit, n):
    """The base window's fields on a surface whose cells' water fractions have changed from wet_run, the surface the
    run stood on, to wet_new: each cell's rain and evaporation moved by the change of its water fraction times the
    contrast between that cell as water and as land, the land's rain raised by one factor with its evaporation
    following, and the surface's rain less evaporation held. Returns the fields and a summary."""
    lat, lon = fields['lat'], fields['lon']
    w = weights(lat, lon)
    land = np.asarray(fields[MASK]) > 0.5                              # the cells the run held as land
    wet_run, wet_new = np.asarray(wet_run, float), np.asarray(wet_new, float)
    if wet_new.shape != land.shape or wet_run.shape != land.shape:
        raise ValueError(f'the water fractions are {wet_run.shape} and {wet_new.shape}, the window\'s cells {land.shape}')
    rain, evaporation = np.asarray(fields[RAIN], float), np.asarray(fields[EVAPORATION], float)
    held = float(((rain - evaporation) * w).sum())                    # the base window's rain less evaporation
    rows = np.arange(land.shape[0])[:, None] * np.ones(land.shape[1], int)
    # each cell's rain as land and as water: the run's own where it held the cell so, the line's or the land beside
    # it where it did not
    near = neighbours_mean(rain, land)
    row = np.array([np.average(rain[i][land[i]], weights=w[i][land[i]]) if land[i].any() else np.nan
                    for i in range(land.shape[0])])
    beside_land = np.where(np.isfinite(near), near, np.where(np.isfinite(row[rows]), row[rows],
                                                             float(np.average(rain[land], weights=w[land]))))
    rain_land = np.where(land, rain, beside_land)
    rain_water = np.where(land, np.maximum(fit['a_mm_day'] + fit['b'] * rain, 0.0), rain)
    change = wet_new - wet_run
    moved = np.maximum(rain + change * (rain_water - rain_land), 0.0)
    # evaporation moves in the cells the run held as land, towards the open water's at their latitude; the cells it
    # held as water keep theirs, which is also what gives the open water's rate by latitude
    zonal = zonal_water_evaporation(fields)
    steamed = np.where(land, np.maximum(evaporation + change * (zonal[rows] - evaporation), 0.0), evaporation)
    # the land's rain raised by one factor until the surface's rain less evaporation is the base window's: in every
    # cell the share that is land takes the factor on its rain as land, and in the cells the run held as land its
    # evaporation follows along the cell's curve
    dry = 1.0 - wet_new

    def with_factor(factor):
        p = moved + dry * rain_land * (factor - 1.0)
        e = steamed + np.where(land, dry * (land_evaporation(rain * factor, rain, evaporation, n) - evaporation), 0.0)
        return p, e

    def left(factor):
        p, e = with_factor(factor)
        return float(((p - e) * w).sum()) - held

    lo, hi = 0.2, 5.0
    if left(lo) > 0.0 or left(hi) < 0.0:
        raise ValueError('no factor on the land\'s rain between 0.2 and 5 closes the water budget')
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if left(mid) < 0.0 else (lo, mid)
    factor = 0.5 * (lo + hi)
    new_rain, new_evaporation = with_factor(factor)
    out = dict(fields)
    out[RAIN], out[EVAPORATION] = new_rain, new_evaporation
    total = w.sum()
    mean = lambda f, m: float((f * w)[m].sum() / w[m].sum()) if m.any() else None       # noqa: E731
    more, less = change > 0.1, change < -0.1
    summary = dict(
        water_share=[float((wet_run * w).sum() / total), float((wet_new * w).sum() / total)],
        mean_change_of_water_fraction=float((np.abs(change) * w).sum() / total),
        cells=dict(more_water_by_a_tenth=int(more.sum()), less_water_by_a_tenth=int(less.sum())),
        land_rain_factor=factor,
        more_water=dict(rain_mm_day=[mean(rain, more), mean(new_rain, more)],
                        evaporation_mm_day=[mean(evaporation, more), mean(new_evaporation, more)]),
        less_water=dict(rain_mm_day=[mean(rain, less), mean(new_rain, less)],
                        evaporation_mm_day=[mean(evaporation, less), mean(new_evaporation, less)]),
        held_as_land=dict(rain_mm_day=[mean(rain, land), mean(new_rain, land)],
                          evaporation_mm_day=[mean(evaporation, land), mean(new_evaporation, land)]),
        held_as_water=dict(rain_mm_day=[mean(rain, ~land), mean(new_rain, ~land)],
                           evaporation_mm_day=[mean(evaporation, ~land), mean(new_evaporation, ~land)]),
        surface=dict(rain_mm_day=[float((rain * w).sum() / total), float((new_rain * w).sum() / total)],
                     evaporation_mm_day=[float((evaporation * w).sum() / total), float((new_evaporation * w).sum() / total)]),
        rain_less_evaporation_left_mm_day=float(((new_rain - new_evaporation) * w).sum() / total - held / total))
    return out, summary


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


def record_path(run, first, last, label):
    return RESULTS / f'surface_change_{run}_y{first}-{last}_{label}.json'


def build(base, mask_path, pairs, tag, log=print):
    """The estimate: base as (run, first, last, label), the new mask's record, the pairs' windows as the same and
    the tag added to the base's label; returns the record."""
    run, first, last, base_label = base
    clim, clim_meta, clim_path = read_product('climatology', *base)
    wind, wind_meta, wind_path = read_product('surface_wind', *base)
    mask_record = json.loads(Path(mask_path).read_text())
    if mask_record.get('schema') != MASK_SCHEMA:
        raise SystemExit(f'{mask_path} is not a design mask record')
    if mask_record['run'] != run or not np.allclose(mask_record['lat'], clim['lat'], atol=1e-3) or \
            not np.allclose(mask_record['lon'], clim['lon'], atol=1e-3):
        raise SystemExit(f'{mask_path} is not a mask on {run}\'s grid and frame')
    wet_new = 1.0 - np.asarray(mask_record['land_fraction'], float)
    surface = PRODUCTS / clim_meta['configuration']['frame']['surface']      # the surface the run stood on
    run_surface = boundary.read_t21(surface)
    if not np.allclose(run_surface['lat'], clim['lat'], atol=1e-3):
        raise SystemExit(f'{surface} is not on the window\'s grid')
    wet_run = 1.0 - run_surface['land']
    measured, pooled = [], []
    for one, other in pairs:
        a, _, a_path = read_product('climatology', *one)
        b, _, b_path = read_product('climatology', *other)
        w = weights(a['lat'], a['lon'])
        mean_k = lambda f: float((f['ts_k'] * w).sum() / w.sum())                         # noqa: E731
        cells = switched(a, b)
        pooled.append(cells)
        measured.append(dict(windows=[dict(run=x[0], years=[x[1], x[2]], label=x[3], global_surface_temperature_K=mean_k(f))
                                      for x, f in ((one, a), (other, b))],
                             products={p.name: sha256(p) for p in (a_path, b_path)}, **water_rain_fit(cells)))
        log(f'pair {one[0]}:{one[1]}-{one[2]} / {other[0]}:{other[1]}-{other[2]}: {measured[-1]["cells"]} switched cells, '
            f'rain as water {measured[-1]["a_mm_day"]:.2f} + {measured[-1]["b"]:.3f} x rain as land')
    fit = water_rain_fit({k: np.concatenate([c[k] for c in pooled]) for k in pooled[0]})
    exponent = fit_exponent(clim)
    scaled, summary = carry(clim, wet_run, wet_new, fit, exponent['n'])
    label = f'{base_label}_{tag}' if base_label else tag
    estimate = dict(method='the land-water contrast of settled runs on different masks, applied to a new mask '
                           '(climate.gcm.surface_change)', label=label,
                    base=dict(run=run, years=[first, last], label=base_label, estimate=clim_meta.get('estimate')),
                    mask=dict(record=str(mask_path), sha256=sha256(mask_path)),
                    water_rain=dict(a_mm_day=fit['a_mm_day'], b=fit['b']), land_curve_exponent=exponent['n'],
                    land_rain_factor=summary['land_rain_factor'],
                    record=str(record_path(run, first, last, label).relative_to(HERE.parent.parent)))
    out_clim = product_path('climatology', run, first, last, label)
    out_wind = product_path('surface_wind', run, first, last, label)
    write_estimate(out_clim, scaled, clim_meta, estimate)
    write_estimate(out_wind, wind, wind_meta, dict(estimate, note='the base window\'s surface wind, carried unchanged'))
    return dict(
        schema=SCHEMA, producer=dict(domain='climate', files={'gcm/surface_change.py': sha256(__file__)[:16]}),
        evidence=('An estimate of the base window\'s rain and evaporation on a changed surface, not a run. The '
                  'contrast between land and water is measured where settled windows of the same physics hold '
                  'different masks, and applied in proportion to the change of each cell\'s water fraction; where '
                  'the displaced rain falls is assumed (on the land, by one factor); the land\'s evaporation follows '
                  'its rain along a curve fitted on the base window\'s own land. Not validated against a run on the '
                  'new surface.'),
        reading_rule=('Rain and evaporation in mm/day. pairs gives each pair of windows, the cells that are land in '
                      'one and water in the other and the line rain as water = a + b x rain as land fitted over '
                      'them, with their mean rain and evaporation as land and as water; water_rain is the line over '
                      'all the pairs\' cells, which the estimate applies. land_curve is the exponent of Turc\'s and '
                      'Pike\'s form fitted, with one ceiling, on the base window\'s land within fit_latitude degrees '
                      'of the equator; each land cell then keeps its own ceiling. summary gives the surface\'s water '
                      'share on the run\'s surface and the new one, the mean change of a cell\'s water fraction, the '
                      'cells with a tenth more or less water, the factor on the land\'s rain, [base, estimate] means '
                      'over the cells named (held_as_land and held_as_water by the run\'s mask), and what the '
                      'surface\'s rain less evaporation is left differing from the base window\'s.'),
        estimate={k: v for k, v in estimate.items() if k != 'record'},
        pairs=rounded(measured), water_rain=rounded(fit), land_curve=rounded(dict(exponent, fit_latitude=FIT_LATITUDE)),
        summary=rounded(summary),
        products={str(p.relative_to(HERE.parent.parent)): sha256(p) for p in (out_clim, out_wind)},
        inputs={str(p.relative_to(HERE.parent.parent)): sha256(p) for p in (clim_path, wind_path)},
        ignored_inputs={str(surface.relative_to(HERE.parent.parent)): sha256(surface)})


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('base', help='the window to carry onto the new surface, RUN:first-last[:label]')
    p.add_argument('--mask', required=True, type=Path,
                   help='the new surface\'s record on the run\'s grid and frame (the atlas round\'s design mask), '
                        'whose land_fraction gives each cell\'s dry share')
    p.add_argument('--pair', action='append', required=True,
                   help='ONE:first-last[:label]/OTHER:first-last[:label]: windows at one temperature on different masks')
    p.add_argument('--tag', default='resurfaced', help='added to the base window\'s label')
    args = p.parse_args(argv)
    pairs = [tuple(parse_window(w) for w in text.split('/')) for text in args.pair]
    base = parse_window(args.base)
    record = build(base, args.mask, pairs, args.tag, log=lambda *a: print(*a, flush=True))
    path = record_path(*base[:3], record['estimate']['label'])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=1) + '\n')
    s = record['summary']
    print(f"water {s['water_share'][0]:.4f} -> {s['water_share'][1]:.4f} of the surface; {s['cells']['more_water_by_a_tenth']} "
          f"cells with a tenth more water, {s['cells']['less_water_by_a_tenth']} with a tenth less")
    print(f"land rain x {s['land_rain_factor']:.4f}; surface rain {s['surface']['rain_mm_day'][0]:.4f} -> "
          f"{s['surface']['rain_mm_day'][1]:.4f} mm/day")
    print(f'wrote {path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
