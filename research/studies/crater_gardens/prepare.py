"""Snapshot a consistent regional water and climate context for the crater-garden count."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from research.studies.crater_gardens.count import DEFAULT_WORK, sha256

CLIMATE = 'climate/gcm/products/climatology_A28R_301k_design_l15_snow_s938_turned_refit_precession_y10-49_today.npz'
ATLAS = 'geography/products/atlas_28pct_4ppd.npz'


def prepare(source_root, work_dir, reuse_hydrology=False):
    source_root, work_dir = Path(source_root), Path(work_dir)
    hydro = work_dir / 'hydrology'
    climate_file, atlas_file = source_root / CLIMATE, source_root / ATLAS
    with np.load(climate_file, allow_pickle=False) as z:
        cm = json.loads(str(z['metadata']))
        climate = {k: z[k] for k in ('lat', 'lon', 'pr_mm_day', 'ts_k')}
    if cm['schema'] != 'terluna.climate.gcm-climatology/1' or cm.get('frame') != 'today':
        raise ValueError('Use the corrected design climate carried to today\'s coordinates')
    if not reuse_hydrology:
        subprocess.run([sys.executable, '-m', 'geography.drainage', '--climatology', str(climate_file),
                        '--out', str(hydro)], cwd=source_root, check=True)
    drainage = hydro / 'products/drainage_28pct_16ppd.npz'
    summary = json.loads((hydro / 'results/drainage.json').read_text())
    if summary['sources'].get(str(climate_file.relative_to(source_root))) != sha256(climate_file):
        raise ValueError('The drainage must have been computed on the same climate product')
    with np.load(atlas_file, allow_pickle=False) as z:
        am = json.loads(str(z['metadata']))
        common = z['water_label'] != 0
    with np.load(drainage, allow_pickle=False) as z:
        dm = json.loads(str(z['metadata']))
        grid = dm['grid']
        if dm['schema'] != 'terluna.geography.drainage/1':
            raise ValueError('Unexpected drainage schema')
        if grid['rows'] % common.shape[0] or grid['cols'] % common.shape[1]:
            raise ValueError('Water grids must align in the present-day lunar frame')
        water = common.repeat(grid['rows'] // common.shape[0], axis=0).repeat(grid['cols'] // common.shape[1], axis=1).astype(np.uint8)
        water.ravel()[z['lake_cell']] |= 2
    meta = dict(schema='terluna.research.crater-gardens-context/1', frame='today',
                evidence='The atlas\'s common-level water at 4 pixels/degree, repeated categorically onto the '
                         '16 pixels/degree rain-fed lake grid. The latter is recomputed with the same corrected '
                         'design climate used for the regional rain and mean ground-temperature exclusions.',
                reading_rule='water has bit 1 for common-level water and bit 2 for rain-fed lakes; zero means '
                             'neither mask covers the cell. Water rows run north to south, columns east from '
                             '0 degrees. Climate fields retain their own Gaussian grid and are sampled bilinearly. '
                             'The common-water mask resolves about 7.6 km and the lakes about 1.9 km at the equator.',
                climate=dict(run=cm['run'], years=cm['years'], frame=cm['frame'], carried=cm['carried']),
                lake_summary=summary['lakes'],
                sources={str(atlas_file): sha256(atlas_file), str(climate_file): sha256(climate_file),
                         str(drainage): sha256(drainage)},
                hydrology_producer=summary['producer'],
                producer=dict(file='research/studies/crater_gardens/prepare.py', sha256=sha256(__file__),
                              source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=source_root, text=True).strip()))
    target = work_dir / 'context.npz'
    np.savez_compressed(target, metadata=json.dumps(meta), water=water, **climate)
    print(json.dumps(dict(file=str(target), bytes=target.stat().st_size, sha256=sha256(target)), indent=2))
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path('/home/kir/Documents/Projects/terluna'))
    parser.add_argument('--work-dir', type=Path, default=DEFAULT_WORK)
    parser.add_argument('--reuse-hydrology', action='store_true')
    args = parser.parse_args()
    prepare(args.source_root, args.work_dir, args.reuse_hydrology)
