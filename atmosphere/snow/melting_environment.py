"""The air that snow falls through below the freezing level in the design-climate storm boxes.

The two CM1 boxes on the design climate save their fields every three hours: the land box at 0 N, 0 E
(box_0e_design, its second lunar day) and the Near-side Sea box at 5 S, 45 W (box_near_side_sea_day3, its third
lunar day, run on from the second's restart).  For every saved column carrying snow through its freezing level, this
keeps the column's temperature, pressure, vapour pressure, cloud condensate, rain and CM1's own melting and
rain-evaporation tendencies on the model levels from the ground to about 32 km, with the freezing level, the snow
there and the local hour angle.  Columns are read as the boxes save them; nothing is simulated.

    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting_environment      # a few minutes, one thread

The arrays go to the research drive (ENV_DIR); melting.py reads them and records their hashes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from atmosphere.snow import turbulence as t
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT, CM1_VAPOUR_GAS_CONSTANT, SYNODIC_MONTH_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = ROOT / 'climate/crm/runs'
ENV_DIR = ROOT / 'research/runs/snow_agent_3/environment'
FREEZING_K = 273.15
SNOW_AT_FREEZING = 1e-5          # kg snow per kg air at the first level above the freezing level
TOP_M = 32e3                     # levels kept: ground to the first level above this height
FIELDS = ('th', 'prs', 'qv', 'qc', 'qi', 'qs', 'qr', 'tt_melt', 'tt_evar')
DAY_S = SYNODIC_MONTH_DAYS * 86400.
# Each box over the lunar day it is read for (days from the start of its run).
CASES = {
    'land': dict(run='box_0e_design', days=(SYNODIC_MONTH_DAYS, 2 * SYNODIC_MONTH_DAYS)),
    'sea': dict(run='box_near_side_sea_day3', days=(2 * SYNODIC_MONTH_DAYS, 3 * SYNODIC_MONTH_DAYS)),
}


def hour_angle(start_deg, lon_deg, t_s):
    """Local hour angle in degrees (0 noon, -90 sunrise, 90 sunset), as the climate branch's box analyses take it:
    the box starts with local time start_deg at 0 E and the Sun moves 360 degrees in one solar day."""
    return float((start_deg + lon_deg + 360. * t_s / DAY_S + 180.) % 360. - 180.)


def case_layout(run):
    path = RUNS / run
    shape, fields, length = t.layout(path / 'cm1out_s.ctl')
    _, mf, ml = t.layout(path / 'cm1out_metadata.ctl')
    metadata = np.fromfile(path / 'cm1out_metadata.dat', dtype='<f4').reshape(-1, ml)
    lookup = {int(row[mf['nwrite'][0]]): float(row[mf['mtime'][0]]) for row in metadata}
    numbers = sorted(int(p.name[8:14]) for p in path.glob('cm1out_t*_s.dat'))
    lines = (path / 'cm1out_s.ctl').read_text().splitlines()
    i = next(k for k, line in enumerate(lines) if line.startswith('zdef'))
    nz = int(lines[i].split()[1])
    zh = np.array([float(x) for line in lines[i + 1:i + 1 + nz] for x in line.split()][:nz]) * 1e3
    return dict(path=path, shape=shape, fields=fields, length=length, numbers=numbers,
                times=np.array([lookup[n] for n in numbers]), zh=zh,
                case=json.loads((path / 'case.json').read_text()))


def read_frame(item, number, digest, levels):
    path = item['path'] / f'cm1out_t{number:06d}_s.dat'
    if path.stat().st_size != item['length'] * 4:
        raise ValueError(f'GrADS size mismatch: {path}')
    raw = np.memmap(path, dtype='<f4', mode='r', shape=(item['length'],))
    s, out = item['shape'], {}
    for name in FIELDS:
        offset, nz = item['fields'][name]
        field = np.asarray(raw[offset:offset + nz * s['y'] * s['x']], dtype='<f4').reshape(nz, -1)[:levels].copy()
        if not np.all(np.isfinite(field)) or np.any(field < -9e5):
            raise ValueError(f'Invalid {name}: {path}')
        digest.update(name.encode() + b'\0' + field.tobytes(order='C'))
        out[name] = field.astype(float)
    return out


def snowy_columns(item, frame):
    """Columns with snow at the first level above their freezing level: that level's index, the freezing height
    (where the air first falls through 0 C going up, interpolated) and the snow's mass there, rho qs."""
    zh = item['zh'][:frame['th'].shape[0]]
    T = frame['th'] * (frame['prs'] / t.CM1_P00) ** (CM1_DRY_AIR_GAS_CONSTANT / t.CM1_CP)
    rho = frame['prs'] / (CM1_DRY_AIR_GAS_CONSTANT * T)
    warm = T >= FREEZING_K
    cross = warm[:-1] & ~warm[1:]
    has = cross.any(axis=0) & warm[0]
    k = np.argmax(cross, axis=0)
    cols = np.flatnonzero(has)
    k = k[cols]
    above = k + 1
    qs = frame['qs'][above, cols]
    keep = qs >= SNOW_AT_FREEZING
    cols, k = cols[keep], k[keep]
    Tk, Tk1 = T[k, cols], T[k + 1, cols]
    z0 = zh[k] + (Tk - FREEZING_K) / (Tk - Tk1) * (zh[k + 1] - zh[k])
    return cols, k, z0, rho[k + 1, cols] * frame['qs'][k + 1, cols], T


def extract(name):
    spec = CASES[name]
    item = case_layout(spec['run'])
    levels = int(np.argmax(item['zh'] > TOP_M)) + 1
    record = item['case']['configuration']
    start, lon = record['start_hour_angle_deg'], record['site']['lon_deg']
    eps = CM1_DRY_AIR_GAS_CONSTANT / CM1_VAPOUR_GAS_CONSTANT
    digest = hashlib.sha256()
    keep = {k: [] for k in ('T', 'p', 'e', 'cond', 'qr', 'qs', 'melt', 'evar')}
    meta = {k: [] for k in ('time_s', 'hour_deg', 'z0_m', 'snow_kg_m3', 'frame', 'column')}
    frames = columns = 0
    lo, hi = spec['days'][0] * 86400. - 1., spec['days'][1] * 86400. + 1.
    for number, time_s in zip(item['numbers'], item['times']):
        if not lo <= time_s <= hi:
            continue
        frame = read_frame(item, number, digest, levels)
        frames += 1
        columns += frame['th'].shape[1]
        cols, k, z0, snow, T = snowy_columns(item, frame)
        if not len(cols):
            continue
        qv = frame['qv'][:, cols]
        p = frame['prs'][:, cols]
        keep['T'].append(T[:, cols].T.astype('<f4'))
        keep['p'].append(p.T.astype('<f4'))
        keep['e'].append((qv / (eps + qv) * p).T.astype('<f4'))
        keep['cond'].append((frame['qc'][:, cols] + frame['qi'][:, cols]).T.astype('<f4'))
        keep['qr'].append(frame['qr'][:, cols].T.astype('<f4'))
        keep['qs'].append(frame['qs'][:, cols].T.astype('<f4'))
        keep['melt'].append(frame['tt_melt'][:, cols].T.astype('<f4'))
        keep['evar'].append(frame['tt_evar'][:, cols].T.astype('<f4'))
        meta['time_s'].append(np.full(len(cols), time_s))
        meta['hour_deg'].append(np.full(len(cols), hour_angle(start, lon, time_s)))
        meta['z0_m'].append(z0)
        meta['snow_kg_m3'].append(snow)
        meta['frame'].append(np.full(len(cols), number))
        meta['column'].append(cols)
    arrays = {k: np.concatenate(v) for k, v in keep.items()}
    arrays.update({k: np.concatenate(v) for k, v in meta.items()})
    arrays['zh_m'] = item['zh'][:levels]
    ENV_DIR.mkdir(parents=True, exist_ok=True)
    out = ENV_DIR / f'{name}.npz'
    np.savez(out, **arrays)
    summary = dict(run=spec['run'], days=list(spec['days']), frames=frames, columns_examined=columns,
                   snowy_at_freezing=int(len(arrays['z0_m'])), levels=levels, fields=list(FIELDS),
                   fields_sha256=digest.hexdigest(), file=str(out.relative_to(ROOT)),
                   file_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
                   start_hour_angle_deg=start, lon_deg=lon, surface=record['site']['surface'])
    (ENV_DIR / f'{name}.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


def load(name):
    """The kept columns of one box and the extraction's record."""
    record = json.loads((ENV_DIR / f'{name}.json').read_text())
    with np.load(ENV_DIR / f'{name}.npz') as z:
        data = {k: z[k] for k in z.files}
    return data, record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=tuple(CASES), action='append')
    args = parser.parse_args()
    for name in args.case or CASES:
        s = extract(name)
        print(name, s['frames'], 'frames;', s['snowy_at_freezing'], 'of', s['columns_examined'],
              'columns carry snow through the freezing level', flush=True)


if __name__ == '__main__':
    main()
