"""The temperatures of the water each large river brings and of the sea it enters, from the project's lake model.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.river_oxygen.waters   # about 35 minutes

A river of the design surface reaches its sea a few kilometres to a few hundred below a lake large enough to set its
temperature (mouths.py), so the water it brings is that lake's surface water. For each river kept one by one this
module runs two water columns of biosphere/ecology/lake.py on the design years' weather, as the cold waters study
runs its ponds and lake (research/studies/cold_waters/run.py, whose forcing and readers it uses):

- the source lake, at its own level and to its mean depth (at most COLUMN_M), under the weather of the land cell
  nearest it;
- the sea beside the mouth, a column of COLUMN_M at the sea's level, under the weather of the water cell nearest the
  mouth.

Each column's temperature profile is kept every PROFILE_EVERY steps (about six hours), with its surface temperature,
its temperature at its floor and its mixed depth step by step, on the research drive (river-oxygen/columns/), one
file for each distinct column; results/waters.json records them with their hashes. The design years are the run A28R_301k_design_l15_snow_s934, years 10-19, on today's Moon with the
earlier atlas's waters: the design surface's own climate has no record written every three days, so the contrast
between a lake and the sea beside it is read on the earlier design climate, 0.7 K warmer in the mean.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.cold_waters import run as CW
from shared.constants import MOON_SURFACE_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MOUTHS = HERE / 'results' / 'mouths.json'
OUT = HERE / 'results' / 'waters.json'
COLUMNS = Path('/media/projectspace/terluna-research/river-oxygen/columns')
SCHEMA = 'terluna.research.river-oxygen-waters/1'
COLUMN_M = 40.0                          # a sea column, and the most a lake column is given
LAYER_M = 1.0
SHALLOW_LAYERS = 20                      # a lake shallower than COLUMN_M is run in this many layers, at least 0.25 m
SPIN_MEANS = 240                         # 3-day means left out before reading: two years, as the cold waters' lake
PROFILE_EVERY = 6                        # steps between kept temperature profiles (a step is 3,575 s)
EARLIER_SEA_LEVEL_M = -1654.2            # the design years' sea level above the dry Moon's geoid (geography's atlas)
SEARCH_DEG = 8.0                         # how far a column looks for a cell of its own kind (land for a lake, water for a sea)


def sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def cell_for(lat, lon, grid_lat, grid_lon, water, want_water):
    """The model cell nearest a place, of the wanted kind where one lies within SEARCH_DEG of arc, else the nearest."""
    la, lo = np.radians(grid_lat)[:, None], np.radians(grid_lon)[None, :]
    a, b = np.radians(lat), np.radians(lon)
    arc = np.degrees(np.arccos(np.clip(np.sin(a) * np.sin(la) + np.cos(a) * np.cos(la) * np.cos(lo - b), -1.0, 1.0)))
    kind = np.where(water == want_water, arc, np.inf)
    j, i = np.unravel_index(np.argmin(kind), arc.shape)
    if not np.isfinite(kind[j, i]) or kind[j, i] > SEARCH_DEG:
        j, i = np.unravel_index(np.argmin(arc), arc.shape)
    return int(j), int(i), bool(water[j, i] == want_water)


def column_key(j, i, level_m, depth_m):
    return f'j{j:02d}_i{i:02d}_level{int(round(level_m / 10.0)) * 10:+06d}_depth{depth_m:05.1f}'


def plan(mouths, data, heights):
    """The distinct columns the rivers need, and each river's pair of them."""
    water = data['lsm'][0] < 0.5
    columns, rivers = {}, {}
    for sea, group in mouths['rivers'].items():
        sea_level = mouths['seas'][sea]['level_m']
        rivers[sea] = []
        for m in group['kept']:
            src = m['course']['source']
            entry = dict(cell=m['cell'])
            for role, lat, lon, level, depth, want_water in (
                    ('lake', src['lat'], src['lon_east'], src['level_m'],
                     float(np.clip(src.get('mean_depth_m', COLUMN_M), 1.0, COLUMN_M)), False),
                    ('sea', m['lat'], m['lon_east'], sea_level, COLUMN_M, True)):
                j, i, matched = cell_for(lat, lon, data['lat'], data['lon'], water, want_water)
                depth = round(depth, 1) if depth < COLUMN_M else COLUMN_M
                key = column_key(j, i, level, depth)
                above_sea = level - EARLIER_SEA_LEVEL_M
                columns.setdefault(key, dict(j=j, i=i, lat=float(data['lat'][j]), lon_east=float(data['lon'][i]),
                                             cell_is_water=bool(water[j, i]), level_above_the_earlier_sea_m=above_sea,
                                             cell_ground_above_the_earlier_sea_m=float(heights[j, i]),
                                             lift_m=above_sea - float(heights[j, i]), depth_m=depth))
                entry[role] = key
                entry[f'{role}_cell_of_its_kind'] = matched
            rivers[sea].append(entry)
    return columns, rivers


def run_column(data, spec, lapse):
    forcing, air_k, window = CW.site_forcing(data, spec['j'], spec['i'], spec['lift_m'], lapse)
    layer = LAYER_M if spec['depth_m'] >= SHALLOW_LAYERS * LAYER_M else max(0.25, spec['depth_m'] / SHALLOW_LAYERS)
    rd, cp = CW.air_constants()
    basin = L.Basin.pan(spec['depth_m'], layer, 1e4)                       # the water and its setting as CW.run_water's
    scenario = L.Scenario('temperature', gpp=1.6, attenuation=1.7 / 3.0, sediment=0.3, deep_water=0.3)
    lake = L.Lake(basin, L.Physics(deep_mixing='wind'), [scenario], L.Air(rd, cp, MOON_SURFACE_GRAVITY))
    state = L.initial_state(basin, [scenario], float(np.mean(air_k)) - 273.15)
    start = SPIN_MEANS * CW.STEPS_PER_MEAN
    kept = range(start, len(forcing), PROFILE_EVERY)
    out = lake.run(forcing, state, profiles_at=kept)
    series = {k: np.asarray(out[k], dtype=np.float32)[start:] for k in ('surface_t', 'bottom_t', 'mixed_depth')}
    series['profile_t'] = np.array([out['profiles'][n]['t'] for n in kept], dtype=np.float32)
    series['layer_top_m'] = np.asarray(basin.z[:-1], dtype=np.float32)
    series['layer_bottom_m'] = np.asarray(basin.z[1:], dtype=np.float32)
    series['pressure_pa'] = np.asarray(forcing.pressure, dtype=np.float32)[start:]
    series['air_c'] = (np.asarray(air_k, dtype=np.float32) - 273.15)[start:]
    return series, float(forcing.dt), layer


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    p.add_argument('--limit', type=int, default=None, help='run at most this many new columns (a trial)')
    args = p.parse_args(argv)
    started = time.time()
    log = lambda *a: print(*a, flush=True)                                     # noqa: E731
    mouths = json.loads(MOUTHS.read_text())
    data = CW.read_design_years()
    heights, topo = CW.topography_m(data['lat'], data['lon'])
    land = data['lsm'][0] > 0.5
    lapse = float(np.mean(data['lapse'][:, land]))
    columns, rivers = plan(mouths, data, heights)
    log(f'{len(columns)} distinct columns for {sum(len(v) for v in rivers.values())} rivers; lapse {lapse * 1e3:.3f} K/km')
    COLUMNS.mkdir(parents=True, exist_ok=True)
    done = 0
    for key, spec in sorted(columns.items()):
        path = COLUMNS / f'{key}.npz'
        if not path.exists():
            if args.limit is not None and done >= args.limit:
                continue
            t0 = time.time()
            series, dt, layer = run_column(data, spec, lapse)
            np.savez_compressed(path, dt_s=dt, layer_m=layer, **series)
            done += 1
            log(f'{key}: {time.time() - t0:.0f} s; surface {series["surface_t"].mean():.2f} C '
                f'({series["surface_t"].min():.2f} to {series["surface_t"].max():.2f}), floor {series["bottom_t"].mean():.2f} C')
        with np.load(path) as z:
            spec.update(file=path.name, sha256=sha256(path), dt_s=float(z['dt_s']), layer_m=float(z['layer_m']),
                        surface_mean_c=float(z['surface_t'].mean()), surface_lowest_c=float(z['surface_t'].min()),
                        surface_highest_c=float(z['surface_t'].max()), floor_mean_c=float(z['bottom_t'].mean()),
                        floor_lowest_c=float(z['bottom_t'].min()), mixed_depth_median_m=float(np.median(z['mixed_depth'])),
                        mixed_depth_p95_m=float(np.percentile(z['mixed_depth'], 95)))
    missing = [k for k, v in columns.items() if 'file' not in v]
    record = dict(
        schema=SCHEMA,
        producer=dict(files={str(Path(f).resolve().relative_to(ROOT)): sha256(f)
                             for f in (__file__, CW.__file__, CW.L.__file__)}),
        inputs={str(MOUTHS.relative_to(ROOT)): sha256(MOUTHS), str(CW.LAKE_OXYGEN.relative_to(ROOT)): sha256(CW.LAKE_OXYGEN)},
        ignored_inputs={str(topo): sha256(topo), **{str(f): sha256(f) for f in CW.year_files()}},
        evidence=('Runs of the project\'s one-dimensional lake model on one climate model\'s 3-day means at T21, made '
                  'linear in time, as the cold waters study runs its lake: no lake or sea like these has been measured. '
                  'A lake far above or below its cell\'s mean ground takes the cell\'s air along the land\'s mean lapse '
                  'rate. A column stands for the surface water a river draws from its lake and for the sea beside the '
                  'mouth; neither feels the other, and the sea column is 40 m of a sea hundreds of metres deep, so its '
                  'floor gives the coldest water the nights\' mixing makes there.'),
        reading_rule=('Temperatures in C, depths in m. columns[key] gives each distinct column: its model cell, whether '
                      'that cell is water in the design years\' run, its level and the lift above the cell\'s mean ground, '
                      'its depth, and its surface\'s and floor\'s temperatures after two years\' spin-up. rivers[sea] '
                      'lists, in mouths.json\'s order, each river\'s lake column and sea column, and whether each stands '
                      'on a cell of its own kind (land for a lake, water for a sea). The series are on the research '
                      'drive under river-oxygen/columns/, one .npz for each column: profile_t holds the temperature of '
                      'each layer (layer_top_m to layer_bottom_m) every sixth step, and surface_t, bottom_t and '
                      'mixed_depth every step.'),
        design_years=dict(run=CW.RUN, years=list(CW.YEARS)), lapse_k_per_km=lapse * 1e3,
        settings=dict(column_m=COLUMN_M, layer_m=LAYER_M, spin_means=SPIN_MEANS, profile_every=PROFILE_EVERY, search_deg=SEARCH_DEG,
                      earlier_sea_level_m=EARLIER_SEA_LEVEL_M),
        columns=columns, rivers=rivers, columns_not_yet_run=missing, seconds=round(time.time() - started))
    args.out.write_text(json.dumps(record, indent=1) + '\n')
    log(f'wrote {args.out}; {done} columns run now, {len(missing)} not yet run, {time.time() - started:.0f} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
