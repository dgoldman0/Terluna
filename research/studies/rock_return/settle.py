"""Does the ground settle under the return as designed, or do its hollows go on filling? Five hundred million years.

    /media/projectspace/terluna-research/venvs/optics/bin/python -m research.studies.rock_return.settle   # about an hour on three cores; wants numba

process.py stops at a hundred million years, where the hollows of lake country have lost a fifth of their volume
with no clear trend and the hollows of dry crater country a quarter, at a steady rate. This module runs the reading
the register takes, rock shared by wear with all of it reached and set down within 20 km, for five hundred million
years, and at every fifty million counts the hollows (hollows.py's census), the relief, the slopes, the forms kept
and how much rock is being cut, which says whether the land's wear is winding down or holding.

The dry region is run both ways with water, since the model treats rock that slides into a hollow differently in
the two: onto a lake it is held as rock come to rest and offered to the carriers again; onto a dry floor it is
ground. Running the dry southern highlands with lakes in their hollows shows what clearing the floors would do.
Writes results/settle.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from research.studies.rock_return import evolve as EV
from research.studies.rock_return import hollows as HO

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'settle.json'
SCHEMA = 'terluna.research.rock-return-settle/1'
YEARS = tuple(5.0e7 * k for k in range(1, 11))
RUNS = dict(korolev_with_lakes=('korolev', 'wet'), southern_highlands_dry=('southern_highlands', 'dry'),
            southern_highlands_with_lakes=('southern_highlands', 'wet'))


def one(job):
    label, region, water = job
    started = time.time()
    g, h0 = EV.grid_of(region), EV.start_of(region)
    g['rng'] = np.random.default_rng(20261122)
    shape = g['shape']
    inside = EV.interior(region)
    land = inside & region['eroding'].ravel()
    rate, _ = EV.calibrate(region, water)
    start, old = HO.census(h0, g, inside.reshape(shape))
    relief0 = EV.local_relief(h0, shape)[land].mean()
    f0 = EV.forms(h0, shape)[inside]
    h, pending, t, cut_before = h0.copy(), None, 0.0, 0.0
    total = dict(cut=0.0, set_down=0.0, gone=0.0, slid=0.0, laid=0.0)
    rows = {}
    for mark in YEARS:
        while t < mark - 1e-6:
            h, moved, pending = EV.step(h, g, water, 'worn', rate, EV.STEP_YR, pending=pending)
            for k in total:
                total[k] += moved[k]
            t += EV.STEP_YR
        row, _ = HO.census(h, g, inside.reshape(shape), old)
        z = h.reshape(shape)
        fall = np.maximum(np.abs(np.diff(z, axis=0, prepend=z[:1])) / g['north_south'], np.abs(np.diff(z, axis=1, prepend=z[:, :1])) / g['east_west'])
        row.update(volume_over_start=row['volume_km3'] / start['volume_km3'], relief_over_start=float(EV.local_relief(h, shape)[land].mean() / relief0),
                   forms_kept=float(np.corrcoef(f0, EV.forms(h, shape)[inside])[0, 1]),
                   mean_fall_to_a_neighbour_deg=float(np.degrees(np.arctan(fall.ravel()[land])).mean()),
                   highest_m=float(h[inside].max()), mean_change_m=float((h - h0)[inside].mean()),
                   rock_cut_km3_per_myr=(total['cut'] - cut_before) / 1e9 / 50.0, rock_gone_km3=total['gone'] / 1e9, rock_slid_km3=total['slid'] / 1e9)
        cut_before = total['cut']
        rows[f'{mark:g}'] = row
    return label, dict(water=water, rate=rate, start=dict(start, highest_m=float(h0[inside].max())), years=rows, seconds=round(time.time() - started))


def build(workers=3, log=print):
    from multiprocessing import Pool
    started = time.time()
    regions = EV.load_regions(('korolev', 'southern_highlands'))
    small = dict(ground=np.zeros((8, 8), np.float32), water=np.full((8, 8), np.nan, np.float32), sea=np.zeros((8, 8), bool),
                 eroding=np.ones((8, 8), bool), lowering=np.zeros((8, 8), np.float32), east_west_m=1.0, north_south_m=1.0)
    for water in EV.WATERS:
        EV.step(EV.start_of(small) + np.arange(64.0), EV.grid_of(small), water, 'even', 1.0e-7, 1.0)
    jobs = [(label, regions[name], water) for label, (name, water) in RUNS.items()]
    results = {}
    with Pool(workers) as pool:
        for label, row in pool.imap_unordered(one, jobs):
            results[label] = row
            log(f'{label}: {row["seconds"]} s; volume over start {[round(r["volume_over_start"], 2) for r in row["years"].values()]}, {time.time() - started:.0f} s')
    placement = json.loads(EV.PLACEMENT.read_text())
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.rock_return.settle',
                      code_sha256={p.name: EV.sha256(p) for p in (Path(__file__), HERE / 'evolve.py', HERE / 'hollows.py')}, compiled=EV.COMPILED),
        inputs={'research/studies/rock_return/results/placement.json': EV.sha256(EV.PLACEMENT)},
        ignored_inputs={str(EV.REGIONS): placement['product']['sha256']},
        evidence=('evolve.py\'s landscape model run for five hundred million years under the reading the register '
                  'takes, rock shared by wear, all of it reached and set down within 20 km. Over such a span the '
                  'model\'s fixed rain, its one strength of rock and its open edge weigh more than over a hundred '
                  'million years; it shows a tendency and predicts no place. Not validated.'),
        reading_rule=('runs: for each region and water, at every fifty million years: hollows.py\'s census (number by '
                      'class of size, share of the ground in a hollow, volume below the lips, hollows on old and new '
                      'ground); volume_over_start; relief_over_start; forms_kept; the mean fall to a neighbour; the '
                      'highest ground; rock_cut_km3_per_myr over the fifty million years before, which falls if the '
                      'land\'s wear is winding down.'),
        setting=dict(years=list(YEARS), runs={k: list(v) for k, v in RUNS.items()}, rule='worn', range_m=EV.RANGE_M, step_yr=EV.STEP_YR,
                     critical_deg=EV.CRITICAL_DEG, creep_m2_yr=EV.CREEP_M2_YR),
        runs=results, seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--workers', type=int, default=3)
    args = p.parse_args(argv)
    out = build(workers=args.workers, log=lambda *a: print(*a, flush=True))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
