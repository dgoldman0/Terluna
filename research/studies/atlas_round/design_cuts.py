"""The works' cuts on the design surface as trenches: the cells whose ground the cuts from Moscoviense take down to the
water they carry, the water's level in each and the flow through it.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.design_cuts   # about 6 minutes

The register's rows of 16 November divide Moscoviense's water at two weirs and send the rest north through a cut in its
northern ridge, to Humboldtianum and on through further cuts to Smythii-Marginis. governed_drainage.py finds those
routes (redrawn_drainage's CutSearch from Moscoviense: the route whose deepest cut is shallowest, and among those the one
whose cuts add up to the least), and the design routing (design_routing.py) puts each branch's water in past each cut;
neither keeps the cut's ground. This module finds the same routes on the design routing's climate and works and keeps
each trench. The channel runs at its water's level from the lake it leaves, over each pass it must cut, to the hollow it
falls into; before the first pass and after the last its way runs down the steepest slope to the water, and between
passes through the ground that stands least above the water. The globe surface (globe_surface.py) takes the trenches'
ground down to the water and draws the branches' flow through them. The trenches are kept apart from the routing so
that the routing's product, which the climate run on the design surface and the maps name by its hash, stays as it is.

The product, atlas-round/design_cuts/design_cuts_16ppd.npz on the research drive, indexes cells as the design routing
does (rows from 90 N, columns from 0 E, 16 per degree). results/design_cuts.json records its inputs, the design
routing's record among them, and its hash; load() gives the product back once it matches that record and the record's
inputs and code are current.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from geography import drainage as dr
from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'design_cuts.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/design_cuts')
PRODUCT = BULK / 'design_cuts_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-design-cuts/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-design-cuts-grid/1'
CODE = ('research/studies/atlas_round/design_cuts.py', 'research/studies/atlas_round/design_routing.py',
        'research/studies/atlas_round/final_drainage.py',
        'research/studies/atlas_round/governed_drainage.py', 'research/studies/atlas_round/redrawn_drainage.py',
        'research/studies/drainage_network/network.py', 'geography/drainage.py')
TRENCH_PAD = 40                      # cells round two passes in which the trench between them is sought
TRENCH_HIGH_M = 50.0                 # the way between passes keeps low: a step costs its length times 1 + height above water / this


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def downhill(cell, h_flat, receiver, level, limit=4000):
    """The cells from `cell` down the routing's steepest way while the ground stands above `level`."""
    out = []
    while h_flat[cell] > level and len(out) < limit:
        out.append(int(cell))
        nxt = int(receiver[cell])
        if nxt == cell:
            break
        cell = nxt
    return out


def low_way(a, b, h, level, cell_km, lat):
    """The 8-connected way from cell a to cell b through the ground that stands least above `level`, within TRENCH_PAD
    cells of both and wrapping in longitude: each step costs its length on the ground times 1 + the next cell's height
    above the water over TRENCH_HIGH_M."""
    import heapq
    rows, cols = h.shape
    ra, ca, rb, cb = a // cols, a % cols, b // cols, b % cols
    db = (cb - ca + cols // 2) % cols - cols // 2          # b's column from a's, the short way round
    r0, r1 = max(min(ra, rb) - TRENCH_PAD, 0), min(max(ra, rb) + TRENCH_PAD, rows - 1)
    c0, c1 = min(0, db) - TRENCH_PAD, max(0, db) + TRENCH_PAD
    goal, best, back = (rb, db), {(ra, 0): 0.0}, {}
    heap = [(0.0, ra, 0)]
    while heap:
        cost, r, dc = heapq.heappop(heap)
        if (r, dc) == goal:
            break
        if cost > best.get((r, dc), np.inf):
            continue
        for sr in (-1, 0, 1):
            for sc in (-1, 0, 1):
                nr, nc = r + sr, dc + sc
                if (sr == 0 and sc == 0) or not (r0 <= nr <= r1 and c0 <= nc <= c1):
                    continue
                step = cell_km * np.hypot(sr, sc * np.cos(np.radians(lat[nr])))
                rise = max(0.0, float(h[nr, (ca + nc) % cols]) - level)
                nxt = cost + step * (1.0 + rise / TRENCH_HIGH_M)
                if nxt < best.get((nr, nc), np.inf):
                    best[(nr, nc)], back[(nr, nc)] = nxt, (r, dc)
                    heapq.heappush(heap, (nxt, nr, nc))
    if goal not in best:
        raise RuntimeError(f'no way between the passes at cells {a} and {b} within {TRENCH_PAD} cells')
    way, node = [], goal
    while node != (ra, 0):
        way.append(node[0] * cols + (ca + node[1]) % cols)
        node = back[node]
    return way + [a]


def trenches(geo, state, gov, log=print):
    """The works' cuts as trenches: for each branch from Moscoviense (governed_drainage.branches_for's chains, found
    again by the same search), the channel runs at its water's level from the lake it leaves, over each pass it must
    cut, to the hollow it falls into, and the ground on that way above the water is taken down to it. Before the first
    pass and after the last the way runs down the steepest slope to the water; between two passes it is the way through
    the ground that stands least above the water (low_way). Returns the trench's cells, the water's level in each and the
    flow through it, the two branches' shares added where they share it."""
    started = time.time()
    levels, in_sea = RD.node_levels(geo, state)
    search = RD.CutSearch(geo, levels, in_sea, state['names'], GD.SOURCE)
    allowed = search.shallowest(only=(GD.HELD_BRANCH, GD.REST_BRANCH))
    h = geo['h']
    h_flat = h.ravel()
    receiver = state['r']['receiver']
    level_at, flow_at = {}, {}
    for target in (GD.HELD_BRANCH, GD.REST_BRANCH):
        chain, _ = search.route(target, allowed[target])
        segments, current = [], None
        for e, lev in chain:
            if current is None and search.wt_l[e] - lev >= RD.CUT_COUNTED_M:
                current = []                              # the channel starts at the first cut
            if current is not None:
                current.append((e, lev))
                y = int(search.dst[e])
                if (lev if lev < search.lev_l[y] else search.lev_l[y]) == search.lev_l[y]:
                    segments.append(current)              # it falls into a hollow and runs on by itself
                    current = None
        if current:
            segments.append(current)
        mine = {}
        for seg in segments:
            lev = seg[0][1]
            way = downhill(int(search.c_from[seg[0][0]]), h_flat, receiver, lev)
            for (e1, _), (e2, _) in zip(seg[:-1], seg[1:]):
                way += low_way(int(search.c_to[e1]), int(search.c_from[e2]), h, lev, geo['cell_km'], geo['lat'])
            for e, _ in seg:
                way += [int(search.c_from[e]), int(search.c_to[e])]
            way += downhill(int(search.c_to[seg[-1][0]]), h_flat, receiver, lev)
            for c in way:
                if h_flat[c] > lev:
                    mine[c] = min(mine.get(c, np.inf), lev)
        for c, lev in mine.items():
            level_at[c] = min(level_at.get(c, np.inf), lev)
            flow_at[c] = flow_at.get(c, 0.0) + float(gov['share'][target])
    cells = np.array(sorted(level_at), np.int32)
    level = np.array([level_at[c] for c in cells], np.float32)
    flow = np.array([flow_at[c] for c in cells], np.float32)
    deepest = float((h_flat[cells] - level).max()) if cells.size else 0.0
    log(f'the cuts\' trenches: {cells.size} cells, deepest {deepest:.1f} m, {time.time() - started:.0f} s')
    return dict(cell=cells, level=level, flow=flow, deepest_m=deepest)


def inputs():
    """The design routing's inputs and its record, with their hashes, and the climatology (ignored by Git)."""
    tracked, ignored = DR.inputs()
    tracked[DR.rel(DR.OUT)] = sha256(DR.OUT)
    return tracked, ignored


def build(log=print):
    """The trenches on the design routing's climate and works: the same steps as design_routing.route up to the governed
    shares, then trenches()."""
    started = time.time()
    fd = json.loads(FD.OUT.read_text())
    case = fd['cases'][DR.BUDGET]
    frame_path, record = FD.frame_record()
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    surf = RD.load_surface(log)
    runoff, lake_net, _, _ = FD.load(case['evaporation'])
    geo = RD.geometry(surf, runoff, lake_net, log)
    state = RD.design_state(geo, runoff, lake_net, log=log)
    branches, info = GD.branches_for(geo, state, rotation, log)
    gov = GD.govern(geo, state, runoff, lake_net, branches, log=log)
    trench = trenches(geo, state, gov, log)
    cols = geo['h'].shape[1]
    h_flat = geo['h'].ravel()
    depth = h_flat[trench['cell']] - trench['level']
    deepest = int(trench['cell'][np.argmax(depth)]) if depth.size else None
    area = np.repeat((geo['cell_km'] * 1e3) ** 2 * np.cos(np.radians(geo['lat'])), cols)
    shares = {k: round(float(gov['share'][k]), 1) for k in (GD.SALT_SEA, GD.HELD_BRANCH, GD.REST_BRANCH)}
    tracked, ignored = inputs()
    meta = dict(schema=PRODUCT_SCHEMA, grid=dict(pixels_per_degree=dr.PPD, rows=geo['h'].shape[0], cols=cols),
                common_level_m=float(geo['level0']), shares_m3s=shares, inputs=tracked, ignored_inputs=ignored,
                reading_rule=('Cells index rows from 90 N and columns from 0 E, 16 per degree, on today\'s Moon: cell = '
                              'row * 5760 + column. cut_cell lists the trenches\' cells, cut_level_m the water\'s level in '
                              'each above the dry Moon\'s geoid (the ground there stands above it until cut) and '
                              'cut_flow_m3s the branches\' flow through it, both branches\' shares added where they share '
                              'the trench.'),
                evidence=('The routes governed_drainage.py finds for the works of 16 November on the design climate\'s '
                          'closed budget, with their ground between 1.9 km cells: a design, not validated.'))
    BULK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PRODUCT, metadata=json.dumps(meta), cut_cell=trench['cell'], cut_level_m=trench['level'],
                        cut_flow_m3s=trench['flow'])
    summary = dict(trench_cells=int(trench['cell'].size), deepest_m=round(trench['deepest_m'], 1),
                   deepest_at=NW.place(deepest, cols, rotation) if deepest is not None else None,
                   recorded_cuts={k: v['cuts'][:3] for k, v in info.items()},
                   ground_taken_km3=round(float((np.maximum(depth, 0.0) * area[trench['cell']]).sum()) / 1e9, 1),
                   trench_km=round(float(np.sqrt(area[trench['cell']]).sum()) / 1e3),
                   flow_m3s=round(float(trench['flow'].max())) if trench['flow'].size else 0, shares_m3s=shares)
    log(f'cuts: {json.dumps(summary)[:800]}, {time.time() - started:.0f} s')
    return dict(schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}), inputs=tracked,
                ignored_inputs=ignored, evidence=meta['evidence'],
                reading_rule=('Depths and flows in m and m3/s; ground taken in km3; the trench\'s length summed over its '
                              'cells. product names the trenches on the research drive, which design_cuts.load() reads.'),
                summary=summary, product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)),
                seconds=round(time.time() - started))


def stale():
    """What in the record no longer matches: the product, a tracked input or the code. Empty when current."""
    if not OUT.exists() or not PRODUCT.exists():
        return ['no cuts yet']
    record = json.loads(OUT.read_text())
    tracked, ignored = inputs()
    out = [] if sha256(PRODUCT) == record['product']['sha256'] else [str(PRODUCT)]
    out += [k for k, v in {**tracked, **ignored}.items() if record['inputs'].get(k, record['ignored_inputs'].get(k)) != v]
    out += [f for f in CODE if record['producer']['files'].get(f) != sha256(ROOT / f)]
    return out


def load(log=print):
    """The trenches, once they match their record and the record's inputs and code are current."""
    behind = stale()
    if behind:
        raise SystemExit(f'the design cuts are stale ({", ".join(behind)}); '
                         'run python -m research.studies.atlas_round.design_cuts')
    with np.load(PRODUCT, allow_pickle=False) as z:
        out = {k: z[k] for k in z.files if k != 'metadata'}
        out['meta'] = json.loads(str(z['metadata']))
    log(f'the design cuts from {PRODUCT}')
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    record = build(lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1) + '\n')
    print(f'wrote {args.out}: {json.dumps(record["summary"])[:1200]}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
