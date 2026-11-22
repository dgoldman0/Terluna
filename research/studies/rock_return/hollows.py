"""Do hollows fill, and do new ones form, under the return as designed? A count on the runs process.py keeps.

    python -m research.studies.rock_return.hollows      # about a minute; reads process.py's snapshots on the research drive

process.py follows the ground under two readings of the placement and measures the largest hollows of the starting
map. Whether the land keeps a population of hollows is a different question: old hollows can fill while new ones
form. This module counts every hollow of each region at the start and after 30 and 100 million years, three ways
that need not agree: by number in each class of size, by the share of the ground that lies in a hollow, and by the
volume the hollows hold below their lips. It then sorts the later hollows by the ground they stand on: on ground
that was hollow at the start (an old hollow, or a piece of one), on ground that was not (a new hollow), or on both.
Writes results/hollows.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

from research.studies.rock_return import evolve as EV
from research.studies.rock_return import process as PR

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'hollows.json'
SCHEMA = 'terluna.research.rock-return-hollows/1'
PROCESS = HERE / 'results' / 'process.json'
CLASSES_KM2 = ((0.0, 10.0, 'under 10 km2'), (10.0, 100.0, '10-100 km2'), (100.0, 1000.0, '100-1,000 km2'), (1000.0, np.inf, 'over 1,000 km2'))
COUNTED_KM2 = 10.0                   # hollows smaller than this are near the model's grain and are left out of the sorting
IN_A_HOLLOW_M = 1.0                  # ground this far under its hollow's lip counts as in the hollow
YEARS = ('3e+07', '1e+08')


def census(h, g, inside, old=None):
    """The hollows of a ground whose cells reach the measured part: their number by class of size, the share of the
    measured ground in a hollow, the volume below the lips, and, against the ground that was hollow at the start,
    how many stand on old ground, on new ground or on both."""
    depth, labels, count = EV.hollows_of(h, g)
    depth = depth.reshape(g['shape'])
    index = np.arange(1, count + 1)
    km2 = g['cell_area'] / 1e6
    size = ndimage.sum(np.ones(g['shape']), labels, index=index) * km2
    reaches = ndimage.maximum(inside.astype(float), labels, index=index) > 0
    now = depth > IN_A_HOLLOW_M
    out = dict(number={label: int((reaches & (size >= lo) & (size < hi)).sum()) for lo, hi, label in CLASSES_KM2},
               share_of_ground=float(now[inside].mean()), volume_km3=float(depth[inside].sum() * km2 / 1e3))
    if old is not None:
        on_old = ndimage.mean(old.astype(float), labels, index=index)
        deepest = ndimage.maximum(depth, labels, index=index)
        counted = reaches & (size >= COUNTED_KM2)
        kinds = dict(on_old_ground=counted & (on_old >= 0.9), on_new_ground=counted & (on_old < 0.1),
                     on_both=counted & (on_old >= 0.1) & (on_old < 0.9))
        out['hollows_of_10_km2_or_more'] = {
            name: dict(number=int(pick.sum()), area_km2=float(size[pick].sum()),
                       median_depth_m=float(np.median(deepest[pick])) if pick.any() else None,
                       deepest_m=float(deepest[pick].max()) if pick.any() else None) for name, pick in kinds.items()}
        out['ground_in_a_hollow_now_and_not_at_the_start'] = float((now & ~old)[inside].mean())
        out['ground_in_a_hollow_at_the_start_and_not_now'] = float((old & ~now)[inside].mean())
    return out, now


def build(log=print):
    process = json.loads(PROCESS.read_text())
    snapshots = Path(process['snapshots']['path'])
    if EV.sha256(snapshots) != process['snapshots']['sha256']:
        raise SystemExit(f'{snapshots} does not match the snapshots {PROCESS} records')
    regions = EV.load_regions(tuple(PR.DESIGN))
    out = {}
    with np.load(snapshots, allow_pickle=False) as snap:
        for name, region in regions.items():
            g, h0 = EV.grid_of(region), EV.start_of(region)
            inside = EV.interior(region).reshape(g['shape'])
            start, old = census(h0, g, inside)
            start['hollows_of_10_km2_or_more'] = sum(v for k, v in start['number'].items() if k != CLASSES_KM2[0][2])
            out[name] = dict(water=PR.DESIGN[name]['water'], start=start)
            for reading in PR.READINGS:
                out[name][reading] = {}
                for year in YEARS:
                    row, _ = census(snap[f'{name}__{reading}__all_of_it__{year}'].astype(np.float64).ravel(), g, inside, old)
                    out[name][reading][year] = row
                    log(f'{name} {reading} {year}: {row["number"]}, {row["share_of_ground"]:.3f} of the ground, {row["volume_km3"]:.0f} km3')
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.rock_return.hollows',
                      code_sha256={p.name: EV.sha256(p) for p in (Path(__file__), HERE / 'evolve.py', HERE / 'process.py')}),
        inputs={'research/studies/rock_return/results/process.json': EV.sha256(PROCESS)},
        ignored_inputs={str(snapshots): process['snapshots']['sha256']},
        evidence=('A count on the ground of process.py\'s runs with all the rock reached: every hollow that the model\'s '
                  'flood finds, at the start and after 30 and 100 million years. A hollow is ground below the lip its '
                  'water would spill over, at 1.9 km a cell; hollows under 10 km2 are a few cells and near the model\'s '
                  'grain. Whether a hollow holds a lake depends on the region\'s water. Not validated.'),
        reading_rule=('For each region and reading: number, the hollows by class of size; share_of_ground, the share '
                      'of the measured ground in a hollow; volume_km3, the volume below the lips; and for hollows of '
                      '10 km2 or more, how many stand on ground that was hollow at the start (nine tenths or more of '
                      'their ground), on ground that was not (under a tenth) or on both, with their area and depths. '
                      'A starting hollow that is cut in two by rock set down across it counts afterwards as two '
                      'hollows on old ground.'),
        setting=dict(classes_km2=[c[2] for c in CLASSES_KM2], counted_km2=COUNTED_KM2, in_a_hollow_m=IN_A_HOLLOW_M, years=list(YEARS)),
        regions=out)


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
