"""A check of the design surface's water against its own heights: does the ground beside each lake hold it?

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python3 -m research.studies.atlas_round.lake_support   # about a minute

A lake's surface is level, so ground beside a lake that lies below its surface is under the lake too, or the lake
stands lower than drawn. The design surface (globe_surface.py) takes its lakes from the routing (geography/drainage.py's
route, through the design routing), which decides each hollow's floor by itself: every cell that drains to one lowest
point is one catchment, and the routing fills a catchment to the lowest pass on its way to a sea, the level at which
the whole hollow round it would run over. A hollow of several catchments is then filled piece by piece: a piece with
water to spare stands full to the outer pass, while the piece beside it, whose evaporation outruns its inflow, keeps a
small lake or none, with nothing between the two but the hollow's floor.

This check measures where that has happened, on the 1.9 km cells:

- the lake cells with dry ground beside them (one of the eight neighbours) lower than the lake's surface by more than
  DROP_M, and the bodies of touching lake cells that have such a shore;
- the dry ground below an adjoining lake's surface and joined to it across ground below that surface: what the lake
  would cover if it kept its level;
- for each pair of lake cell and lower dry neighbour, whether the two lie in one catchment or in two, on the surface
  and seas the routing used.

The seas are measured the same way and reported apart: the near-side sea's main part stands over its neck by design.
The check changes nothing; it is the measure a corrected lake rule has to bring to zero.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

from geography import drainage as dr
from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import globe_surface as GS
from research.studies.atlas_round import redrawn_drainage as RD
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'lake_support.json'
SCHEMA = 'terluna.research.atlas-round-lake-support/1'
CODE = ('research/studies/atlas_round/lake_support.py',)
DROP_M = 1.0                         # ground this far below a lake's surface counts as below it
NEIGHBOURS = [(dr_, dc) for dr_ in (-1, 0, 1) for dc in (-1, 0, 1) if dr_ or dc]


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def shifted(a, d_row, d_col, fill):
    """The grid moved so that each cell holds its neighbour's value: rows end at the poles (filled), columns wrap."""
    out = np.roll(a, -d_col, axis=1)
    if d_row:
        out = np.roll(out, -d_row, axis=0)
        if d_row > 0:
            out[-d_row:] = fill
        else:
            out[:-d_row] = fill
    return out


def unheld(h, level, source, drop=DROP_M):
    """The water cells of `source` with a dry neighbour lower than their surface by more than `drop`. level is NaN on
    dry ground."""
    dry = np.isnan(level)
    out = np.zeros(h.shape, bool)
    for d_row, d_col in NEIGHBOURS:
        out |= source & shifted(dry, d_row, d_col, False) & (shifted(h, d_row, d_col, np.inf) < level - drop)
    return out


def ground_below(h, level, source, drop=DROP_M):
    """The level that reaches each dry cell from the water cells of `source`, spreading from cell to neighbour across dry
    ground lower than that level by more than `drop`; the highest level where several reach. -inf where none does."""
    rows, cols = h.shape
    hf = h.ravel()
    dry = np.isnan(level).ravel()
    reach = np.full(hf.size, -np.inf, np.float32)
    front = np.flatnonzero(source.ravel())
    at = level.ravel()[front].astype(np.float32)
    while front.size:
        r, c = front // cols, front % cols
        found, found_level = [], []
        for d_row, d_col in NEIGHBOURS:
            rr = r + d_row
            ok = (rr >= 0) & (rr < rows)
            j = rr[ok] * cols + (c[ok] + d_col) % cols
            lv = at[ok]
            take = dry[j] & (hf[j] < lv - drop) & (reach[j] < lv)
            found.append(j[take])
            found_level.append(lv[take])
        j, lv = np.concatenate(found), np.concatenate(found_level)
        if not j.size:
            break
        order = np.lexsort((lv, j))
        j, lv = j[order], lv[order]
        last = np.r_[j[1:] != j[:-1], True]                 # the highest level reaching each cell
        front, at = j[last], lv[last]
        reach[front] = at
    return reach.reshape(rows, cols)


def measure(h, level, source, area, drop=DROP_M):
    """The counts and areas for one kind of water (source): its shore, the part of it not held, the bodies that part
    lies on, and the dry ground below and joined. Returns them with the level reaching each dry cell."""
    dry = np.isnan(level)
    shore = np.zeros(h.shape, bool)
    for d_row, d_col in NEIGHBOURS:
        shore |= source & shifted(dry, d_row, d_col, False)
    loose = unheld(h, level, source, drop)
    bodies, n = ndimage.label(source, structure=np.ones((3, 3), bool))
    marked = np.zeros(n + 1, bool)
    marked[np.unique(bodies[loose])] = True
    marked[0] = False
    reach = ground_below(h, level, source, drop)
    below = np.isfinite(reach)
    depth = (reach - h)[below]
    total = float(area.sum())
    water = float(area[source].sum())
    return dict(
        cells=int(source.sum()), area_km2=round(water), share_of_moon=round(water / total, 5),
        shore_cells=int(shore.sum()), shore_cells_not_held=int(loose.sum()),
        shore_share_not_held=round(float(loose.sum()) / max(int(shore.sum()), 1), 4),
        bodies=int(n), bodies_with_such_a_shore=int(marked.sum()),
        area_in_those_bodies_km2=round(float(area[marked[bodies]].sum())),
        area_share_in_those_bodies=round(float(area[marked[bodies]].sum()) / max(water, 1.0), 4),
        ground_below_and_joined_km2=round(float(area[below].sum())),
        ground_below_share_of_moon=round(float(area[below].sum()) / total, 5),
        ground_below_over_water_area=round(float(area[below].sum()) / max(water, 1.0), 4),
        ground_below_depth_m=(dict(median=round(float(np.median(depth))), p90=round(float(np.percentile(depth, 90))),
                                   largest=round(float(depth.max()))) if depth.size else None)), reach


def catchment_pairs(h, level, source, outlet, lat, cell_km, spills, drop=DROP_M):
    """The pairs of a water cell of `source` and a lower dry neighbour: how many lie in one catchment (the cells that
    drain to one lowest point) and how many in two, and in how many the water's own catchment overflows."""
    rows, cols = h.shape
    receiver = dr.d8_receivers(h.astype(np.float32), outlet, lat, cell_km)
    floor = dr.terminal(receiver)
    hf, lv, dry = h.ravel(), level.ravel(), np.isnan(level).ravel() & ~outlet.ravel()
    cells = np.flatnonzero(source.ravel())
    r, c = cells // cols, cells % cols
    same = other = overflowing = 0
    for d_row, d_col in NEIGHBOURS:
        rr = r + d_row
        ok = (rr >= 0) & (rr < rows)
        i = cells[ok]
        j = rr[ok] * cols + (c[ok] + d_col) % cols
        pair = dry[j] & (hf[j] < lv[i] - drop)
        i, j = i[pair], j[pair]
        same += int((floor[i] == floor[j]).sum())
        other += int((floor[i] != floor[j]).sum())
        overflowing += int(spills.ravel()[i].sum())
    return dict(pairs=same + other, in_one_catchment=same, in_two_catchments=other, the_lake_overflows=overflowing)


def build(log=print):
    started = time.time()
    with np.load(GS.PRODUCT, allow_pickle=False) as z:
        h = z['height_m'].astype(np.float32)
        cells, levels, bodies, spill = z['water_cell'], z['water_level_m'], z['water_body'], z['water_spills']
    rows, cols = h.shape
    level = np.full(rows * cols, np.nan, np.float32)
    level[cells] = levels
    body = np.full(rows * cols, -1, np.int8)
    body[cells] = bodies
    spills = np.zeros(rows * cols, bool)
    spills[cells] = spill
    level, body, spills = (a.reshape(rows, cols) for a in (level, body, spills))
    lat = 90.0 - (np.arange(rows) + 0.5) * 180.0 / rows
    cell_km = MOON_RADIUS / 1000.0 * np.radians(180.0 / rows)
    area = np.repeat((cell_km ** 2 * np.cos(np.radians(lat)))[:, None], cols, axis=1)
    lakes, reach = measure(h, level, body == 0, area)
    log(f"lakes: {lakes['shore_cells_not_held']:,} of {lakes['shore_cells']:,} shore cells not held; "
        f"{lakes['ground_below_and_joined_km2']:,} km2 of dry ground below and joined, {time.time() - started:.0f} s")
    seas, _ = measure(h, level, body > 0, area)
    log(f"seas: {seas['shore_cells_not_held']:,} of {seas['shore_cells']:,} shore cells not held, {time.time() - started:.0f} s")
    bands = {}
    for top in range(90, -90, -30):
        sl = slice(int((90 - top) * rows / 180), int((120 - top) * rows / 180))
        bands[f'{top} to {top - 30}'] = round(float(area[sl][np.isfinite(reach[sl])].sum()))
    # the catchments, on the surface and seas the routing used
    routing = DR.load(log)
    surf = RD.load_surface(log)
    outlet = (routing['sea'] >= 0).reshape(rows, cols)
    lake_level = np.full(rows * cols, np.nan, np.float32)
    lake_level[routing['lake_cell']] = routing['lake_level_m']
    lake_spills = np.zeros(rows * cols, bool)
    lake_spills[routing['lake_cell']] = routing['lake_spills']
    lake_level, lake_spills = lake_level.reshape(rows, cols), lake_spills.reshape(rows, cols)
    pairs = catchment_pairs(surf['h'], lake_level, np.isfinite(lake_level), outlet, surf['lat'], surf['cell_km'], lake_spills)
    log(f"pairs: {pairs['in_two_catchments']:,} of {pairs['pairs']:,} lie in two catchments, {time.time() - started:.0f} s")
    record = json.loads(GS.OUT.read_text())
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={str(GS.OUT.relative_to(ROOT)): sha256(GS.OUT), str(DR.OUT.relative_to(ROOT)): sha256(DR.OUT)},
        ignored_inputs={str(GS.PRODUCT): record['product']['sha256']},
        evidence=('A calculation on products: the design surface\'s water cells and levels compared with its own heights, '
                  'cell by cell, and the routing\'s lakes with the catchments of the surface it ran on. It measures a '
                  'defect of the lake rule and proposes no remedy; what a corrected rule does to the lakes\' area and '
                  'the rivers\' flows is not computed here.'),
        reading_rule=('lakes and seas: a shore cell is a water cell with dry ground among its eight neighbours; it is not '
                      'held where one of them lies more than drop_m below the water\'s surface. bodies are touching '
                      'water cells. ground_below_and_joined is the dry ground lower than an adjoining water surface and '
                      'joined to it across such ground. ground_below_by_latitude gives the lakes\' in km2 by 30-degree '
                      'band. catchment_pairs counts each lake cell and lower dry neighbour on the routing\'s surface: in '
                      'one catchment (cells draining to one lowest point) or two, and whether the lake overflows.'),
        parameters=dict(drop_m=DROP_M),
        lakes=lakes, seas=seas, ground_below_by_latitude_km2=bands, catchment_pairs=pairs,
        seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    record = build()
    OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    print(f"wrote {OUT}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
