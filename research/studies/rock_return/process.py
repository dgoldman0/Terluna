"""What the return as designed does to the ground: where it puts the rock, and how high it builds.

    /media/projectspace/terluna-research/venvs/optics/bin/python -m research.studies.rock_return.process   # about eight minutes on four cores; wants numba

evolve.py sets several rules of placement side by side. This module follows the one the design holds: carriers
take the rock where it comes to rest and set it down on freshly worn ground round that place, most of it within
20 km. The design does not say whether ground can get more than it lost, so the placement is run in two readings:
shared in proportion to the wear, and offered so but taken only up to what each piece of ground lost. It runs each region under the water the design gives it (Korolev's country stands
mostly under lakes, the southern highlands are dry) for a hundred million years, with the share of the rock the
carriers can reach as the case: all of it, the share that divers at the delta fronts add where rivers bring sand,
and the share flyers and waders reach alone.

It answers in plan and in height. In plan: the map of where the ground has risen and fallen. In height: the
spread of heights at the start and after, the change of the ground by where it stood on the starting map (crests,
slopes, floors), the change of the highest summits and of the rims of the largest hollows, and the most any ground
has risen. Pictures go to the research drive. Writes results/process.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

from research.studies.rock_return import evolve as EV

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'process.json'
SCHEMA = 'terluna.research.rock-return-process/1'
SNAPSHOTS = EV.DRIVE / 'process_snapshots.npz'
YEARS = (1.0e6, 3.0e6, 1.0e7, 3.0e7, 1.0e8)
DESIGN = dict(korolev=dict(water='wet', reaches=dict(all_of_it=1.0, with_divers_at_the_delta_fronts=0.95, flyers_and_waders_alone=0.42)),
              southern_highlands=dict(water='dry', reaches=dict(all_of_it=1.0)))
READINGS = dict(shared_by_wear='worn', no_more_than_it_lost='healed')   # two readings of setting rock down on freshly worn ground
PERCENTS = (0.1, 1.0, 10.0, 50.0, 90.0, 99.0, 99.9)
SUMMITS = 30
SUMMIT_CELLS = 5                     # a summit is the highest ground within this many cells


def heights(h, mask):
    v = h[mask]
    return dict(percentiles={f'{p:g}%': float(np.percentile(v, p)) for p in PERCENTS}, highest=float(v.max()), lowest=float(v.min()),
                mean=float(v.mean()), spread=float(v.std()))


def positions(h0, region):
    """Where each cell stood on the starting map: a crest, a floor or between by its height against the ground
    round it, and whether it is among the steepest quarter of the land."""
    shape = region['ground'].shape
    form = EV.forms(h0, shape)
    inside = EV.interior(region)
    s = form[inside].std()
    relief = EV.local_relief(h0, shape)
    steep = relief >= np.percentile(relief[inside], 75)
    return dict(crests=inside & (form > s), floors=inside & (form < -s), between=inside & (np.abs(form) <= s), steepest_quarter=inside & steep)


def by_position(h, h0, where):
    out = {}
    for name, mask in where.items():
        d = (h - h0)[mask]
        out[name] = dict(share_of_ground=float(mask.sum() / where['crests'].size), mean_change_m=float(d.mean()),
                         tenth_m=float(np.percentile(d, 10)), ninetieth_m=float(np.percentile(d, 90)))
    return out


def summits_of(h0, region):
    """The highest distinct summits of the starting map inside the measured part."""
    shape = region['ground'].shape
    g = h0.reshape(shape)
    top = (g == ndimage.maximum_filter(g, size=2 * SUMMIT_CELLS + 1, mode='nearest')) & EV.interior(region).reshape(shape)
    cells = np.flatnonzero(top.ravel())
    return cells[np.argsort(h0[cells])[::-1][:SUMMITS]]


def summit_change(h, h0, cells, shape):
    """Each summit's change: the highest ground within two cells of it now, less its height at the start."""
    now = ndimage.maximum_filter(h.reshape(shape), size=5, mode='nearest').ravel()[cells]
    d = now - h0[cells]
    return dict(median_m=float(np.median(d)), least_m=float(d.min()), most_m=float(d.max()), start_heights_m=[float(h0[cells].min()), float(h0[cells].max())])


def one(job):
    name, reading, label, reach, region, water = job
    started = time.time()
    shots, books, rate = EV.run(region, water, READINGS[reading], years=YEARS, reach=reach, seed=20261122)
    g, h0 = EV.grid_of(region), EV.start_of(region)
    shape = g['shape']
    inside = EV.interior(region)
    followed = EV.followed_hollows(region)
    start_marks = [EV.hollow_now(h0, hol, g) for hol in followed]
    where = positions(h0, region)
    tops = summits_of(h0, region)
    rims0 = np.array([np.median(h0.reshape(shape)[hol['ring']]) for hol in followed])
    rows = {}
    for mark, h in shots.items():
        change = (h - h0)[inside]
        rims = np.array([np.median(h.reshape(shape)[hol['ring']]) for hol in followed]) - rims0
        rows[f'{mark:g}'] = dict(EV.measure(h, h0, region, followed, start_marks), heights=heights(h, inside), by_where_it_stood=by_position(h, h0, where),
                                 summits=summit_change(h, h0, tops, shape),
                                 rims_of_the_followed_hollows_m=dict(median=float(np.median(rims)), least=float(rims.min()), most=float(rims.max())),
                                 change_m={f'{p:g}%': float(np.percentile(change, p)) for p in PERCENTS} | dict(most=float(change.max()), least=float(change.min())),
                                 share_risen_over_100_m=float((change > 100.0).mean()), share_fallen_over_100_m=float((change < -100.0).mean()),
                                 rock={k: v / 1e9 for k, v in books[mark].items()})
    return name, reading, label, dict(reach=reach, water=water, rule=READINGS[reading], rate=rate, start=dict(heights=heights(h0, inside)), years=rows,
                             seconds=round(time.time() - started)), {f'{mark:g}': h.astype(np.float32) for mark, h in shots.items()}


def picture(name, region, water, shots, row, path, reading):
    """The designed return on one region: the map at the start and after, where the ground rose and fell, and the
    spread of heights."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    shape = region['ground'].shape
    h0 = EV.start_of(region).reshape(shape)
    inside = EV.interior(region).reshape(shape)
    light = LightSource(azdeg=315, altdeg=35)
    lo, hi = np.percentile(h0, [1, 99])
    extent = [0, shape[1] * region['east_west_m'] / 1e3, shape[0] * region['north_south_m'] / 1e3, 0]
    fig, axes = plt.subplots(2, 3, figsize=(19, 12.5))

    def shaded(ax, h, title):
        rgb = light.shade(np.clip(h, lo - (hi - lo), hi + (hi - lo)), cmap=plt.cm.gist_earth, vmin=lo, vmax=hi, blend_mode='overlay',
                          vert_exag=6.0, dx=region['east_west_m'], dy=region['north_south_m'])
        ax.imshow(rgb, extent=extent, aspect='equal')
        ax.set_title(title, fontsize=12)
        ax.set_xlabel('km')
    shaded(axes[0, 0], h0, 'the map at the start')
    shaded(axes[0, 1], shots['3e+07'].reshape(shape), 'after 30 million years')
    shaded(axes[0, 2], shots['1e+08'].reshape(shape), 'after 100 million years')
    diff = shots['1e+08'].reshape(shape) - h0
    im = axes[1, 0].imshow(diff, extent=extent, aspect='equal', cmap='RdBu', vmin=-600, vmax=600)
    axes[1, 0].set_title('where the ground rose (blue) and fell (red) in 100 million years, m', fontsize=12)
    axes[1, 0].set_xlabel('km')
    fig.colorbar(im, ax=axes[1, 0], shrink=0.8)
    ax = axes[1, 1]
    bins = np.linspace(np.percentile(h0[inside], 0.05), np.percentile(h0[inside], 99.95) + 500.0, 120)
    for key, label, colour in (('start', 'at the start', 'black'), ('3e+07', '30 million years', 'tab:blue'), ('1e+08', '100 million years', 'tab:red')):
        h = h0 if key == 'start' else shots[key].reshape(shape)
        ax.hist(h[inside], bins=bins, histtype='step', color=colour, label=label, linewidth=1.6)
    ax.set_yscale('log')
    ax.set_xlabel('height of the ground, m')
    ax.set_ylabel('cells')
    ax.set_title('how the ground is spread over height', fontsize=12)
    ax.legend()
    ax = axes[1, 2]
    order = np.argsort(h0[inside])
    start = h0[inside][order]
    edges = np.linspace(0, start.size, 41).astype(int)
    mid = [start[a:b].mean() for a, b in zip(edges[:-1], edges[1:])]
    for key, label, colour in (('3e+07', '30 million years', 'tab:blue'), ('1e+08', '100 million years', 'tab:red')):
        d = (shots[key].reshape(shape) - h0)[inside][order]
        ax.plot(mid, [d[a:b].mean() for a, b in zip(edges[:-1], edges[1:])], color=colour, label=f'mean, {label}')
        ax.fill_between(mid, [np.percentile(d[a:b], 10) for a, b in zip(edges[:-1], edges[1:])],
                        [np.percentile(d[a:b], 90) for a, b in zip(edges[:-1], edges[1:])], color=colour, alpha=0.15)
    ax.axhline(0.0, color='grey', linewidth=0.8)
    ax.set_xlabel('height at the start, m')
    ax.set_ylabel('change of the ground, m (band: a tenth to nine tenths)')
    ax.set_title('change by the height the ground started at', fontsize=12)
    ax.legend()
    s = row['years']['1e+08']['summits']
    fig.suptitle(f'{name.replace("_", " ")}, {region["centre"][0]:g}° N {region["centre"][1]:g}° E, '
                 f'{"lakes in every hollow" if water == "wet" else "dry hollows"}: rock returned to freshly worn ground, '
                 f'{"shared by wear" if reading == "shared_by_wear" else "no ground getting more than it lost"}, all of it reached. '
                 f'Highest summits changed {s["least_m"]:+.0f} to {s["most_m"]:+.0f} m in 100 million years', fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def build(workers=4, log=print):
    from multiprocessing import Pool
    started = time.time()
    regions = EV.load_regions(tuple(DESIGN))
    small = dict(ground=np.zeros((8, 8), np.float32), water=np.full((8, 8), np.nan, np.float32), sea=np.zeros((8, 8), bool),
                 eroding=np.ones((8, 8), bool), lowering=np.zeros((8, 8), np.float32), east_west_m=1.0, north_south_m=1.0)
    for water in EV.WATERS:
        EV.step(EV.start_of(small) + np.arange(64.0), EV.grid_of(small), water, 'even', 1.0e-7, 1.0)
    jobs = [(name, reading, label, reach, regions[name], d['water']) for name, d in DESIGN.items() for reading in READINGS
            for label, reach in d['reaches'].items()]
    results, shots = {}, {}
    with Pool(workers) as pool:
        for name, reading, label, row, snaps in pool.imap_unordered(one, jobs):
            results.setdefault(name, {}).setdefault(reading, {})[label] = row
            shots[(name, reading, label)] = snaps
            log(f'{name} {reading} {label}: {row["seconds"]} s, forms kept {row["years"]["1e+08"]["forms_kept"]:.3f} at 100 million years, {time.time() - started:.0f} s')
    pictures, arrays = {}, {}
    for name, d in DESIGN.items():
        for reading in READINGS:
            path = EV.DRIVE / f'{name}_as_designed_{reading}.png'
            picture(name, regions[name], d['water'], shots[(name, reading, 'all_of_it')], results[name][reading]['all_of_it'], path, reading)
            pictures[f'{name}, {reading}'] = str(path)
            for label in d['reaches']:
                for mark, h in shots[(name, reading, label)].items():
                    arrays[f'{name}__{reading}__{label}__{mark}'] = h.reshape(regions[name]['ground'].shape)
    np.savez_compressed(SNAPSHOTS, **arrays)
    placement = json.loads(EV.PLACEMENT.read_text())
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.rock_return.process',
                      code_sha256={p.name: EV.sha256(p) for p in (Path(__file__), HERE / 'evolve.py')}, compiled=EV.COMPILED),
        inputs={'research/studies/rock_return/results/placement.json': EV.sha256(EV.PLACEMENT)},
        ignored_inputs={str(EV.REGIONS): placement['product']['sha256']},
        snapshots=dict(path=str(SNAPSHOTS), sha256=EV.sha256(SNAPSHOTS)), pictures=pictures,
        evidence=('evolve.py\'s landscape model on two squares of the design surface, with the placement the design '
                  'holds: rock set down on freshly worn ground round where it came to rest, about 86% of it within 20 km, in '
                  'two readings: shared in proportion to the wear, or offered so and taken only up to what each cell lost. Each cell is 1.9 km, so where on a slope the carriers set the rock '
                  'is below what the model sees, and the returned rock wears as the ground does. The shares the '
                  'carriers reach are the study\'s for lakes whose rivers bring sand. Not validated.'),
        reading_rule=('regions: for each region under the water the design gives it, each reading of the placement '
                      '(shared_by_wear, no_more_than_it_lost) and each share of the rock the carriers reach, at each year: evolve.py\'s measures; heights, the spread of the ground over '
                      'height (percentiles, the highest and the lowest); by_where_it_stood, the change of the ground '
                      'that was crest, floor or between on the starting map, and of the steepest quarter; summits, '
                      'the change of the thirty highest summits; the change of the rims of the largest hollows; '
                      'change_m, the change of the ground by share of the land, with the most any cell rose or fell; '
                      'rock, km3 cut, laid in hollows, set down and gone.'),
        setting=dict(years=list(YEARS), design=DESIGN, readings=READINGS, range_m=EV.RANGE_M, summits=SUMMITS, summit_cells=SUMMIT_CELLS, percents=list(PERCENTS)),
        regions=results, seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args(argv)
    out = build(workers=args.workers, log=lambda *a: print(*a, flush=True))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
