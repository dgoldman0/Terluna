"""What the coast-water study reads from the GCM's own output: the design case and the hot world's 20-year test, the
two runs climate_targets builds the resort from, over their years 10-19, and the design case's two years written
every 6.5 hours.

    climate/gcm/.venv/bin/python -m research.studies.coast_water.extract     # about a minute; needs netCDF4

For each run it keeps each cell's mean rain and evaporation over the years, and for every output mean the fields the
night's surface balance and the plants' demand need: the Sun's cosine, the air at 2 m, the ground, the surface
pressure, the surface's net longwave (all-sky and clear), its net sunlight, its sensible and latent heat, the cloud
cover, and the lowest level's humidity, wind and air. The humidity is kept as written, the value at each mean's
first step (ExoPlaSim 3.4.2, plasim.f90 line 3218), and run.py takes each mean's own from its start and its end as
climate_targets does. It writes research/runs/coast_water/gcm_extract.npz (ignored, on the research drive) and a
sidecar gcm_extract.json with the sha256 and size of every file it read.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from research.studies.climate_targets import extract as CT

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / 'research' / 'runs' / 'coast_water' / 'gcm_extract.npz'
SIDECAR = OUT.with_suffix('.json')
SCHEMA = 'terluna.research.coast-water-gcm-extract/1'
RUNS = dict(design=('A28_corrected_970', (10, 19)), test=('A28_corrected_970_lw_ice', (10, 19)),
            eddies=('A28_corrected_970_eddies', (0, 1)))
MEANS = ('pr', 'evap')                                    # each cell's mean over the years, m/s of water
SERIES = ('czen', 'tas', 'ts', 'ps', 'rls', 'rlscs', 'rss', 'hfss', 'hfls', 'clt')
LOWEST = ('hus', 'ua', 'va', 'ta')                        # the lowest model level, about 0.9 km up


def read(name, years):
    """One run's fields over its years, with the hashes of the files read and its sun window."""
    import netCDF4
    files, series, sums, lsm, n = {}, {k: [] for k in SERIES + LOWEST}, {k: 0.0 for k in MEANS}, None, 0
    lat = lon = None
    for year in range(years[0], years[1] + 1):
        path = CT.output(name, year)
        files[str(path.relative_to(ROOT))] = dict(sha256=CT.sha256(path), bytes=path.stat().st_size)
        with netCDF4.Dataset(path) as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            mask = np.asarray(d['lsm'][:], float) > 0.5
            if not (mask == mask[0]).all() or (lsm is not None and not (mask[0] == lsm).all()):
                raise ValueError(f'{name}: the land mask changes')
            lsm = mask[0]
            for k in MEANS:
                v = np.asarray(d[k][:], float)
                sums[k] = sums[k] + v.sum(axis=0)
            n += d['pr'].shape[0]
            for k in SERIES:
                series[k].append(np.asarray(d[k][:], np.float32))
            for k in LOWEST:
                series[k].append(np.asarray(d[k][:, -1], np.float32))
    progress, progress_file = CT.record(name)
    files.update({p: dict(sha256=s, bytes=(ROOT / p).stat().st_size) for p, s in progress_file.items()})
    out = {f'{k}_mean': (sums[k] / n).astype(np.float64) for k in MEANS}
    out.update({k: np.concatenate(v) for k, v in series.items()})
    out['lsm'] = lsm
    meta = dict(run=name, years=list(years), means=n, sun_window_rad=float(CT.sun_window(progress['configuration'])),
                writes_per_lunar_day=progress['configuration']['model']['writes_per_lunar_day'])
    return out, files, meta, lat, lon


def extract():
    arrays, files, runs, grid = {}, {}, {}, None
    for key, (name, years) in RUNS.items():
        out, f, meta, lat, lon = read(name, years)
        if grid is None:
            grid = (lat, lon)
        elif not (np.array_equal(grid[0], lat) and np.array_equal(grid[1], lon)):
            raise ValueError('the runs differ in their grid')
        arrays.update({f'{key}_{k}': v for k, v in out.items()})
        files.update(f)
        runs[key] = meta
    b = CT.boundary()
    order = [int(np.argmin(np.abs(b['lat'] - x))) for x in grid[0]]          # the boundary rows in the output's order
    files[str(CT.BOUNDARY.relative_to(ROOT))] = dict(sha256=CT.sha256(CT.BOUNDARY), bytes=CT.BOUNDARY.stat().st_size)
    arrays.update(lat=grid[0], lon=grid[1], land_fraction=b['land'][order], land_elevation_m=b['elevation'][order],
                  cell_area_fraction=b['area'][order])
    return arrays, dict(schema=SCHEMA, runs=runs, files=files,
                        reading_rule=('Arrays are keyed run_field. Series run along the output means, years in order, '
                                      'as (means, 32, 64) from the north; hus, ua, va and ta are the lowest level, '
                                      'hus as written (the value at each mean\'s first step). pr_mean and evap_mean '
                                      'are m/s of water over every mean of the years; ExoPlaSim writes evaporation '
                                      'as negative. lsm is the model\'s own land mask, land_fraction, '
                                      'land_elevation_m and cell_area_fraction the T21 boundary file\'s, in the '
                                      'output\'s row order.'))


def main(argv=None) -> int:
    arrays, sidecar = extract()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, **arrays)
    sidecar['npz'] = dict(path=str(OUT.relative_to(ROOT)), sha256=CT.sha256(OUT), bytes=OUT.stat().st_size)
    SIDECAR.write_text(json.dumps(sidecar, indent=1) + '\n')
    for key, meta in sidecar['runs'].items():
        print(key, meta['run'], meta['means'], 'means')
    return 0


if __name__ == '__main__':
    sys.exit(main())
