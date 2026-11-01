"""How deep the snow zone of the design-climate storm boxes stands above the freezing level.

The two 192 km CM1 boxes on the design climate (the cold-cloud run at 298.5 K), over land at 0 N, 0 E and over the
Near-side Sea at 5 S, 45 W, save their fields every three hours for 59 days (crm-runs on the project-space drive).
For every saved column holding snow above its freezing level this measures the freezing level (where the air first
falls through 0 C going up), the top of the snow (the highest level holding at least 1e-5 kg of snow per kg of air,
with 1e-6 as a sensitivity), the heights above the freezing level below which half and nine-tenths of the column's
snow mass lie above that level, the temperature at the snow's top, and the same top for cloud ice.  Columns are
counted alone and weighted by the snow mass they hold above the freezing level.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from atmosphere.snow import turbulence as t
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = ROOT / 'climate/crm/runs'
PRODUCT = HERE / 'results/snow_zone_depth.json'
CASES = ('box_0e_design', 'box_near_side_sea')
THRESHOLDS = (1e-5, 1e-6)        # kg snow per kg air
FIELDS = ('th', 'prs', 'qs', 'qi')
FREEZING_K = 273.15


def case_layout(name):
    path = RUNS / name
    shape, fields, length = t.layout(path / 'cm1out_s.ctl')
    _, mf, ml = t.layout(path / 'cm1out_metadata.ctl')
    metadata = np.fromfile(path / 'cm1out_metadata.dat', dtype='<f4').reshape(-1, ml)
    lookup = {int(row[mf['nwrite'][0]]): float(row[mf['mtime'][0]]) for row in metadata}
    numbers = sorted(int(p.name[8:14]) for p in path.glob('cm1out_t*_s.dat'))
    interfaces = np.loadtxt(path / 'input_grid_z')
    return dict(path=path, shape=shape, fields=fields, length=length, numbers=numbers,
                times=np.array([lookup[n] for n in numbers]), interfaces=interfaces,
                z=(interfaces[:-1] + interfaces[1:]) / 2, dz=np.diff(interfaces),
                case=json.loads((path / 'case.json').read_text()))


def read_frame(item, number, digest):
    path = item['path'] / f'cm1out_t{number:06d}_s.dat'
    if path.stat().st_size != item['length'] * 4:
        raise ValueError(f'GrADS size mismatch: {path}')
    raw = np.memmap(path, dtype='<f4', mode='r', shape=(item['length'],))
    out = {}
    s = item['shape']
    for name in FIELDS:
        offset, nz = item['fields'][name]
        field = np.asarray(raw[offset:offset + nz * s['y'] * s['x']], dtype='<f4').reshape(nz, s['y'], s['x']).copy()
        if not np.all(np.isfinite(field)):
            raise ValueError(f'Invalid {name}: {path}')
        digest.update(name.encode() + b'\0' + field.tobytes(order='C'))
        out[name] = field.astype(float).reshape(nz, -1)
    return out


def column_measures(item, frame, threshold):
    """Per column: freezing level, snow top and mass quantiles above it, the top's temperature, the ice top."""
    z, dz = item['z'], item['dz']
    T = frame['th'] * (frame['prs'] / t.CM1_P00) ** (CM1_DRY_AIR_GAS_CONSTANT / t.CM1_CP)
    rho = frame['prs'] / (CM1_DRY_AIR_GAS_CONSTANT * T)
    rows = []
    for c in range(T.shape[1]):
        Tc = T[:, c]
        below = np.flatnonzero((Tc[:-1] >= FREEZING_K) & (Tc[1:] < FREEZING_K))
        if Tc[0] < FREEZING_K or not below.size:
            continue
        k = below[0]
        z0 = z[k] + (Tc[k] - FREEZING_K) / (Tc[k] - Tc[k + 1]) * (z[k + 1] - z[k])
        above = z > z0
        snow = above & (frame['qs'][:, c] >= threshold)
        if not snow.any():
            continue
        top = int(np.flatnonzero(snow).max())
        mass = np.where(above, rho[:, c] * frame['qs'][:, c] * dz, 0.)
        cumulative = np.cumsum(mass) / mass.sum()
        h = z - z0
        ice = above & (frame['qi'][:, c] >= threshold)
        rows.append((z0, z[top] - z0, Tc[top] - FREEZING_K, float(np.interp(.5, cumulative, h)),
                     float(np.interp(.9, cumulative, h)), mass.sum(),
                     z[int(np.flatnonzero(ice).max())] - z0 if ice.any() else np.nan))
    return np.array(rows).reshape(-1, 7)


def quantiles(values, weights=None):
    v = np.asarray(values, float)
    keep = np.isfinite(v)
    v = v[keep]
    if not v.size:
        return None
    if weights is None:
        q = np.percentile(v, (10, 50, 90))
    else:
        w = np.asarray(weights, float)[keep]
        order = np.argsort(v)
        c = np.cumsum(w[order]) / w.sum()
        q = np.interp((.1, .5, .9), c, v[order])
    return dict(p10=float(q[0]), p50=float(q[1]), p90=float(q[2]))


def analyse(stride=1):
    cases = {}
    for name in CASES:
        item = case_layout(name)
        digest = hashlib.sha256()
        per = {thr: [] for thr in THRESHOLDS}
        columns = 0
        for number in item['numbers'][::stride]:
            frame = read_frame(item, number, digest)
            columns += frame['th'].shape[1]
            for thr in THRESHOLDS:
                per[thr].append(column_measures(item, frame, thr))
        out = {}
        for thr in THRESHOLDS:
            a = np.concatenate(per[thr])
            labels = ('freezing_level_km', 'snow_top_above_freezing_km', 'snow_top_temperature_C',
                      'snow_half_mass_above_freezing_km', 'snow_nine_tenths_mass_above_freezing_km',
                      'snow_mass_above_freezing_kg_m2', 'ice_top_above_freezing_km')
            scale = (1e-3, 1e-3, 1., 1e-3, 1e-3, 1., 1e-3)
            heavy = a[:, 5] >= np.percentile(a[:, 5], 90) if len(a) else np.zeros(0, bool)
            out[f'{thr:g}_kg_kg'] = dict(
                snowy_columns=int(len(a)), snowy_fraction=float(len(a) / columns),
                by_column={lab: quantiles(a[:, i] * s) for i, (lab, s) in enumerate(zip(labels, scale))},
                weighted_by_snow_mass={lab: quantiles(a[:, i] * s, a[:, 5]) for i, (lab, s) in enumerate(zip(labels, scale))
                                       if i != 5},
                heaviest_tenth={lab: quantiles(a[heavy, i] * s) for i, (lab, s) in enumerate(zip(labels, scale))})
        cases[name] = dict(purpose=item['case']['configuration'].get('purpose', '')[:240],
                           frames_read=len(item['numbers'][::stride]), stride=stride, cadence_s=float(np.diff(item['times'])[0]),
                           days=[float(item['times'][0] / 86400), float(item['times'][-1] / 86400)],
                           columns_examined=columns, fields_sha256=digest.hexdigest(), thresholds=out)
    return dict(schema='terluna.atmosphere.snow-zone-depth/1',
                evidence=('Diagnosis of the archived design-climate CM1 boxes: heights of resolved snow above each '
                          "column's freezing level; no new simulation."),
                reading_rule=('Heights are km above each column\'s freezing level; by_column counts every snowy column '
                              'alike, weighted_by_snow_mass weights each by its snow above the freezing level; the 1e-05 '
                              'threshold is the reference and 1e-06 a sensitivity.'),
                definitions=dict(freezing_level='lowest height where the air falls through 0 C going up',
                                 snow='CM1 snow mixing ratio qs at or above the threshold',
                                 weighting='columns weighted by their snow mass above the freezing level'),
                cases=cases, producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stride', type=int, default=1)
    parser.add_argument('--dry', action='store_true', help='print without writing the product')
    args = parser.parse_args()
    product = analyse(args.stride)
    if not args.dry:
        PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
    for name, case in product['cases'].items():
        r = case['thresholds']['1e-05_kg_kg']
        print(name, 'snowy columns', r['snowy_columns'], 'fraction %.3f' % r['snowy_fraction'])
        for key in ('freezing_level_km', 'snow_top_above_freezing_km', 'snow_top_temperature_C',
                    'snow_half_mass_above_freezing_km', 'snow_nine_tenths_mass_above_freezing_km', 'ice_top_above_freezing_km'):
            print(f'  {key:42s} column {r["by_column"][key]}  mass-weighted {r["weighted_by_snow_mass"].get(key)}')


if __name__ == '__main__':
    main()
