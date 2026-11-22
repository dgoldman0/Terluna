"""What the return as designed does to the ground: where it puts the rock, and how high it builds.

    /media/projectspace/terluna-research/venvs/optics/bin/python -m research.studies.rock_return.process   # about half an hour on six cores; wants numba

evolve.py sets several rules of placement side by side. This module follows the one the design holds: carriers
take the rock where it comes to rest and set it down on freshly worn ground round that place, most of it within
20 km. The design does not say whether ground can get more than it lost, so the placement is run in two readings:
shared in proportion to the wear, and offered so but taken only up to what each piece of ground lost.

Each region runs under the water the design gives it (Korolev's country stands mostly under lakes, the southern
highlands are dry) for a hundred million years, on evolve.py's model with soil creep and with ground that fails
where it is steeper than 35 degrees. The cases are the share of the rock the carriers reach, a longer range, and
loose ground that fails at 30 degrees. Two checks ask whether the model's own grain and step shape the answer: one
crater's surroundings run on the design surface's cells and on cells half their size, and a region run with half
the time step.

It answers in plan and in height. In plan: where the ground has risen and fallen. In height: the spread of heights
at the start and after, the change of the ground by where it stood on the starting map, the change of the highest
summits and of the rims of the largest hollows, the most any ground has risen, and how much of the ground stands
at the failing slope. Pictures go to the research drive. Writes results/process.json.
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
SCHEMA = 'terluna.research.rock-return-process/2'
SNAPSHOTS = EV.DRIVE / 'process_snapshots.npz'
YEARS = (1.0e6, 3.0e6, 1.0e7, 3.0e7, 1.0e8)
READINGS = dict(shared_by_wear='worn', no_more_than_it_lost='healed')   # two readings of setting rock down on freshly worn ground
DESIGN = dict(
    southern_highlands=dict(water='dry', cases=dict(all_of_it={}, all_of_it_on_loose_ground=dict(critical_deg=EV.LOOSE_DEG))),
    korolev=dict(water='wet', cases=dict(all_of_it={}, with_divers_at_the_delta_fronts=dict(reach=0.95), flyers_and_waders_alone=dict(reach=0.42),
                                         all_of_it_within_50_km=dict(range_m=EV.FAR_RANGE_M), all_of_it_on_loose_ground=dict(critical_deg=EV.LOOSE_DEG))))
WINDOWS = dict(southern_highlands=(slice(170, 330), slice(150, 370), 248), korolev=(slice(60, 220), slice(190, 350), 150))   # a closer look: rows, columns, the row of the cut
CHECK_YEARS = (1.0e7, 3.0e7)
CHECK_MARGIN = 24                    # cells of the design surface left out along the edges of the crater's window
PERCENTS = (0.1, 1.0, 10.0, 50.0, 90.0, 99.0, 99.9)
SUMMITS = 30
SUMMIT_CELLS = 5                     # a summit is the highest ground within this many cells
CAP_M = 4000.0                       # pictures draw the ground no higher than this above the window's highest starting ground


def heights(h, mask):
    v = h[mask]
    return dict(percentiles={f'{p:g}%': float(np.percentile(v, p)) for p in PERCENTS}, highest=float(v.max()), lowest=float(v.min()),
                mean=float(v.mean()), spread=float(v.std()))


def steepness(h, region):
    """Each cell's steepest fall to a neighbour, in degrees, and the share of the measured land within two degrees of
    the failing slope or over it."""
    g = EV.grid_of(region)
    z = h.reshape(g['shape'])
    best = np.zeros_like(z)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy or dx:
                d = np.hypot(dy * g['north_south'], dx * g['east_west'])
                fall = (z - np.roll(np.roll(z, dy, 0), dx, 1)) / d
                best = np.maximum(best, fall)
    return np.degrees(np.arctan(best)).ravel()


def positions(h0, region):
    """Where each cell stood on the starting map: a crest, a floor or between by its height against the ground
    round it, and whether it is among the steepest quarter of the land."""
    shape = region['ground'].shape
    scale = region.get('finer', 1)
    form = EV.forms(h0, shape, scale)
    inside = EV.interior(region)
    s = form[inside].std()
    relief = EV.local_relief(h0, shape, scale=scale)
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
    size = 2 * SUMMIT_CELLS * region.get('finer', 1) + 1
    top = (g == ndimage.maximum_filter(g, size=size, mode='nearest')) & EV.interior(region).reshape(shape)
    cells = np.flatnonzero(top.ravel())
    return cells[np.argsort(h0[cells])[::-1][:SUMMITS]]


def summit_change(h, h0, cells, shape, scale=1):
    """Each summit's change: the highest ground within two cells of it now, less its height at the start."""
    now = ndimage.maximum_filter(h.reshape(shape), size=4 * scale + 1, mode='nearest').ravel()[cells]
    d = now - h0[cells]
    return dict(median_m=float(np.median(d)), least_m=float(d.min()), most_m=float(d.max()), start_heights_m=[float(h0[cells].min()), float(h0[cells].max())])


def account(region, water, rule, years, **settings):
    """Run a region and measure it at each year."""
    started = time.time()
    shots, books, rate = EV.run(region, water, rule, years=years, seed=20261122, **settings)
    g, h0 = EV.grid_of(region), EV.start_of(region)
    shape = g['shape']
    scale = region.get('finer', 1)
    inside = EV.interior(region)
    land = inside & region['eroding'].ravel()
    followed = EV.followed_hollows(region)
    start_marks = [EV.hollow_now(h0, hol, g) for hol in followed]
    where = positions(h0, region)
    tops = summits_of(h0, region)
    rims0 = np.array([np.median(h0.reshape(shape)[hol['ring']]) for hol in followed])
    failing = settings.get('critical_deg', EV.CRITICAL_DEG)
    steep0 = steepness(h0, region)
    rows = {}
    for mark, h in shots.items():
        change = (h - h0)[inside]
        steep = steepness(h, region)
        row = dict(EV.measure(h, h0, region, followed, start_marks), heights=heights(h, inside), by_where_it_stood=by_position(h, h0, where),
                   summits=summit_change(h, h0, tops, shape, scale),
                   change_m={f'{p:g}%': float(np.percentile(change, p)) for p in PERCENTS} | dict(most=float(change.max()), least=float(change.min())),
                   share_risen_over_100_m=float((change > 100.0).mean()), share_fallen_over_100_m=float((change < -100.0).mean()),
                   share_risen_over_300_m=float((change > 300.0).mean()), share_risen_over_1000_m=float((change > 1000.0).mean()),
                   share_of_land_at_the_failing_slope=float((steep[land] > failing - 2.0).mean()),
                   mean_slope_deg=float(steep[land].mean()),
                   rock={k: v / 1e9 for k, v in books[mark].items()})
        if followed:
            rims = np.array([np.median(h.reshape(shape)[hol['ring']]) for hol in followed]) - rims0
            row['rims_of_the_followed_hollows_m'] = dict(median=float(np.median(rims)), least=float(rims.min()), most=float(rims.max()))
        rows[f'{mark:g}'] = row
    return dict(settings=settings, rule=rule, water=water, rate=rate,
                start=dict(heights=heights(h0, inside), share_of_land_at_the_failing_slope=float((steep0[land] > failing - 2.0).mean()),
                           mean_slope_deg=float(steep0[land].mean())),
                years=rows, seconds=round(time.time() - started)), {f'{mark:g}': h.astype(np.float32) for mark, h in shots.items()}


def window_of(region, rows, cols, finer=1):
    """A crater's surroundings cut from a region, on the region's cells or on cells finer times smaller (the same
    ground drawn finer: no detail is added)."""
    def cut(a, order):
        w = a[rows, cols]
        return w if finer == 1 else ndimage.zoom(w.astype(np.float64), finer, order=order)
    ground = cut(region['ground'], 1).astype(np.float32)
    flag = lambda a: cut(a.astype(np.float64), 0) > 0.5                                 # noqa: E731
    return dict(ground=ground, water=np.full(ground.shape, np.nan, np.float32), sea=flag(region['sea']), eroding=flag(region['eroding']),
                lowering=cut(region['lowering'], 1).astype(np.float32), east_west_m=region['east_west_m'] / finer,
                north_south_m=region['north_south_m'] / finer, centre=region['centre'], margin=CHECK_MARGIN * finer, finer=finer)


def one(job):
    kind, name, reading, label, region, water, years, settings = job
    row, snaps = account(region, water, READINGS[reading], years, **settings)
    return kind, name, reading, label, row, snaps


# ---- pictures ----------------------------------------------------------------------------------------------------
def shade(ax, w, lo, hi, extent, region, title):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    rgb = LightSource(azdeg=315, altdeg=35).shade(w, cmap=plt.cm.gist_earth, vmin=lo, vmax=hi + 1500.0, blend_mode='overlay', vert_exag=6.0,
                                                   dx=region['east_west_m'], dy=region['north_south_m'])
    ax.imshow(rgb, extent=extent, aspect='equal')
    ax.set_title(title, fontsize=12)
    ax.set_xlabel('km')


def picture(name, region, water, shots, row, path, reading):
    """The designed return on one region: the map at the start and after, where the ground rose and fell, and the
    spread of heights."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    shape = region['ground'].shape
    h0 = EV.start_of(region).reshape(shape)
    inside = EV.interior(region).reshape(shape)
    lo, hi = np.percentile(h0, [1, 99])
    top = h0.max() + CAP_M
    extent = [0, shape[1] * region['east_west_m'] / 1e3, shape[0] * region['north_south_m'] / 1e3, 0]
    fig, axes = plt.subplots(2, 3, figsize=(19, 12.5))
    shade(axes[0, 0], h0, lo, hi, extent, region, 'the map at the start')
    shade(axes[0, 1], np.minimum(shots['3e+07'].reshape(shape), top), lo, hi, extent, region, 'after 30 million years')
    shade(axes[0, 2], np.minimum(shots['1e+08'].reshape(shape), top), lo, hi, extent, region, 'after 100 million years')
    diff = shots['1e+08'].reshape(shape) - h0
    im = axes[1, 0].imshow(diff, extent=extent, aspect='equal', cmap='RdBu', vmin=-600, vmax=600)
    axes[1, 0].set_title('where the ground rose (blue) and fell (red) in 100 million years, m', fontsize=12)
    axes[1, 0].set_xlabel('km')
    fig.colorbar(im, ax=axes[1, 0], shrink=0.8)
    ax = axes[1, 1]
    bins = np.linspace(np.percentile(h0[inside], 0.05), h0[inside].max() + CAP_M, 140)
    for key, label, colour in (('start', 'at the start', 'black'), ('3e+07', '30 million years', 'tab:blue'), ('1e+08', '100 million years', 'tab:red')):
        h = h0 if key == 'start' else shots[key].reshape(shape)
        ax.hist(np.minimum(h[inside], bins[-1]), bins=bins, histtype='step', color=colour, label=label, linewidth=1.6)
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
    at = row['years']['1e+08']
    fig.suptitle(f'{name.replace("_", " ")}, {region["centre"][0]:g}° N {region["centre"][1]:g}° E, '
                 f'{"lakes in every hollow" if water == "wet" else "dry hollows"}: rock returned to freshly worn ground, '
                 f'{"shared by wear" if reading == "shared_by_wear" else "no ground getting more than it lost"}, all of it reached. '
                 f'Highest ground {row["start"]["heights"]["highest"]:.0f} m at the start, {at["heights"]["highest"]:.0f} m after 100 million years', fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def closer(name, region, shots, path, reading):
    """A crater or a shore from near: from above, slanted with heights drawn taller, and a cut through the ground."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    rows, cols, line = WINDOWS[name]
    ew, ns = region['east_west_m'] / 1e3, region['north_south_m'] / 1e3
    h0 = EV.start_of(region).reshape(region['ground'].shape)
    views = [('at the start', h0)] + [(f'after {t} million years', shots[k].reshape(h0.shape)) for t, k in ((30, '3e+07'), (100, '1e+08'))]
    w0 = h0[rows, cols]
    lo, hi = np.percentile(w0, [1, 99])
    top = w0.max() + CAP_M
    x, y = np.arange(w0.shape[1]) * ew, np.arange(w0.shape[0]) * ns
    X, Y = np.meshgrid(x, y)
    tall = 8.0
    fig = plt.figure(figsize=(19, 16))
    for k, (title, h) in enumerate(views):
        w = np.minimum(h[rows, cols], top)
        ax = fig.add_subplot(3, 3, k + 1)
        shade(ax, w, lo, hi, [0, x[-1], y[-1], 0], region, f'from above, {title}')
        ax.axhline((line - rows.start) * ns, color='white', linewidth=0.8, linestyle='--')
        ax = fig.add_subplot(3, 3, k + 4, projection='3d')
        rgb = LightSource(azdeg=315, altdeg=35).shade(w, cmap=plt.cm.gist_earth, vmin=lo, vmax=hi + 1500.0, blend_mode='soft', vert_exag=6.0,
                                                       dx=ew * 1e3, dy=ns * 1e3)
        ax.plot_surface(X, Y, w / 1e3, facecolors=rgb, rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
        ax.set_box_aspect((x[-1], y[-1], (top - w0.min()) / 1e3 * tall))
        ax.set_zlim(w0.min() / 1e3, top / 1e3)
        ax.view_init(elev=32, azim=-62)
        ax.set_title(f'slanted view, heights drawn {tall:g} times too tall, {title}')
        ax.set_xlabel('km')
        ax.set_ylabel('km')
        ax.set_zlabel('km')
    ax = fig.add_subplot(3, 1, 3)
    for (title, h), colour in zip(views, ('black', 'tab:blue', 'tab:red')):
        ax.plot(x, np.minimum(h[line, cols], top), color=colour, label=title, linewidth=1.6 if colour == 'black' else 1.2)
    ax.set_xlabel('km along the dashed line')
    ax.set_ylabel('height of the ground, m')
    ax.set_title('a cut through the ground along the dashed line, heights drawn about 25 times too tall')
    ax.legend()
    ax.grid(alpha=0.3)
    fig.suptitle(f'{name.replace("_", " ")}: a closer look, '
                 f'{"rock shared by wear" if reading == "shared_by_wear" else "ground taking no more than it lost"}', fontsize=15)
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def build(workers=6, log=print):
    from multiprocessing import Pool
    started = time.time()
    regions = EV.load_regions(tuple(DESIGN))
    small = dict(ground=np.zeros((8, 8), np.float32), water=np.full((8, 8), np.nan, np.float32), sea=np.zeros((8, 8), bool),
                 eroding=np.ones((8, 8), bool), lowering=np.zeros((8, 8), np.float32), east_west_m=1.0, north_south_m=1.0)
    for water in EV.WATERS:                                      # compile the kernels once, before the workers start
        EV.step(EV.start_of(small) + np.arange(64.0), EV.grid_of(small), water, 'even', 1.0e-7, 1.0)
    jobs = [('run', name, reading, label, regions[name], d['water'], YEARS, settings)
            for name, d in DESIGN.items() for reading in READINGS for label, settings in d['cases'].items()]
    south = regions['southern_highlands']
    rows, cols, _ = WINDOWS['southern_highlands']
    for reading in READINGS:
        for finer in (1, 2):
            jobs.append(('grain', f'cells_{finer}', reading, f'{finer}', window_of(south, rows, cols, finer), 'dry', CHECK_YEARS, {}))
        jobs.append(('step', 'southern_highlands', reading, 'half_the_step', south, 'dry', CHECK_YEARS, dict(dt=0.5 * EV.STEP_YR)))
    results, checks, shots = {}, dict(grain={}, step={}), {}
    with Pool(workers) as pool:
        for kind, name, reading, label, row, snaps in pool.imap_unordered(one, jobs):
            if kind == 'run':
                results.setdefault(name, {}).setdefault(reading, {})[label] = row
                shots[(name, reading, label)] = snaps
            else:
                checks[kind].setdefault(reading, {})[label] = row
            log(f'{kind} {name} {reading} {label}: {row["seconds"]} s, {time.time() - started:.0f} s')
    pictures, arrays = {}, {}
    for name, d in DESIGN.items():
        for reading in READINGS:
            snaps = shots[(name, reading, 'all_of_it')]
            path = EV.DRIVE / f'{name}_as_designed_{reading}.png'
            picture(name, regions[name], d['water'], snaps, results[name][reading]['all_of_it'], path, reading)
            near = EV.DRIVE / f'{name}_closer_{reading}.png'
            closer(name, regions[name], snaps, near, reading)
            pictures[f'{name}, {reading}'] = dict(whole=str(path), closer=str(near))
            for label in d['cases']:
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
        evidence=('evolve.py\'s landscape model on two squares of the design surface, with soil creep and with ground '
                  'that fails where it is steeper than 35 degrees, and the placement the design holds: rock set down '
                  'on freshly worn ground round where it came to rest, about 86% of it within 20 km, in two readings: '
                  'shared in proportion to the wear, or offered so and taken only up to what each cell lost. Each cell '
                  'is 1.9 km, so where on a slope the carriers set the rock is below what the model sees, and the '
                  'returned rock wears and fails as the ground does. The shares the carriers reach are the study\'s '
                  'for lakes whose rivers bring sand. Not validated.'),
        reading_rule=('regions: for each region under the water the design gives it, each reading of the placement '
                      '(shared_by_wear, no_more_than_it_lost) and each case (the share of the rock reached, the '
                      'longer range, loose ground failing at 30 degrees), at each year: evolve.py\'s measures; heights, '
                      'the spread of the ground over height; by_where_it_stood, the change of the ground that was '
                      'crest, floor or between on the starting map, and of the steepest quarter; summits, the change '
                      'of the thirty highest summits; the change of the rims of the largest hollows; change_m, the '
                      'change of the ground by share of the land; the shares risen over 100, 300 and 1,000 m; the '
                      'share of the land within two degrees of the failing slope; rock, km3. checks: grain, one '
                      'crater\'s surroundings in the southern highlands on the design surface\'s cells and on cells '
                      'half their size; step, the southern highlands with half the time step.'),
        setting=dict(years=list(YEARS), design=DESIGN, readings=READINGS, range_m=EV.RANGE_M, critical_deg=EV.CRITICAL_DEG, loose_deg=EV.LOOSE_DEG,
                     creep_m2_yr=EV.CREEP_M2_YR, step_yr=EV.STEP_YR, margin_cells=EV.MARGIN, check_years=list(CHECK_YEARS),
                     check_margin_cells=CHECK_MARGIN, summits=SUMMITS, summit_cells=SUMMIT_CELLS, percents=list(PERCENTS),
                     windows={k: [v[0].start, v[0].stop, v[1].start, v[1].stop, v[2]] for k, v in WINDOWS.items()}),
        regions=results, checks=checks, seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--workers', type=int, default=6)
    args = p.parse_args(argv)
    out = build(workers=args.workers, log=lambda *a: print(*a, flush=True))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
