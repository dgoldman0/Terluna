"""The design surface's routing on the redrawn floor: each cell's downstream neighbour, discharge and drained area, every
sea at its own level and the lakes, on the design climate's closed budget with the near-side sea parted at its neck.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.design_routing   # about 6 minutes

The design surface is routed once for each climate and kept, so that the maps built on it read one routing: the globe
surface (globe_surface.py), the satoyama map (satoyama_final.py) and the climate model's surface for the run on the
design surface (design_turned_surface.py) read it, and later maps can. The routing:

- the governed routing on the closed budget of the design climate (final_drainage.py's carried climatology, its works
  and shares; governed_drainage.govern): the salt sea and Humboldtianum held at the common level by the fixed weirs,
  and Moscoviense, Orientale and Smythii-Marginis full at their rims;
- the near-side sea parted (near_side_neck.py's steady state), its main part and its southern seas each at its own
  level;
- routed again with both parts at their levels, so that rivers run on over the floor the parted sea leaves dry, and
  with the main part's spill put in past the neck's crest, so that the river over the neck runs down to the southern
  seas.

The product, atlas-round/design_routing/design_routing_16ppd.npz on the research drive, indexes cells from 90 N and
0 E, 16 per degree, and holds for every cell its receiver, discharge and drained area, the sea over it and the
routing's outlet and lake flags, and for the lakes off the seas their cells, levels and whether each spills; its
metadata gives the seas' levels, the neck's spill and the climate. results/design_routing.json records its inputs and
its hash, and load() gives the product back once it matches that record and the record's inputs are current.
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
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import near_side_neck as N
from research.studies.atlas_round import redrawn_drainage as RD

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'design_routing.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/design_routing')
PRODUCT = BULK / 'design_routing_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-design-routing/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-design-routing-grid/1'
CODE = ('research/studies/atlas_round/design_routing.py', 'research/studies/atlas_round/final_drainage.py',
        'research/studies/atlas_round/near_side_neck.py', 'research/studies/atlas_round/governed_drainage.py',
        'research/studies/atlas_round/redrawn_drainage.py', 'research/studies/drainage_network/network.py',
        'geography/drainage.py')
BUDGET = 'closed_budget'
SEAS = ('near_side_main', 'near_side_south', 'south_pole_aitken', 'humboldtianum', 'moscoviense', 'orientale',
        'smythii_marginis')
UPPER_SEAS = ('moscoviense', 'orientale', 'smythii_marginis')       # full at their rims, so they spill over them
SPILL_SEARCH_DEG = 0.6               # the spill is put in at the lowest ground beside the main part this near the throat


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path) -> str:
    return str(Path(path).relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path)


def spill_cell(h_flat, main_water, cols, centre):
    """The cell past the neck's crest where the main part's spill leaves it: the lowest ground beside the main part's
    water within SPILL_SEARCH_DEG of the throat's centre."""
    rows = h_flat.size // cols
    i0 = int((90.0 - centre[0]) * dr.PPD)
    j0 = int((centre[1] % 360.0) * dr.PPD)
    k = int(SPILL_SEARCH_DEG * dr.PPD)
    r_ = np.arange(max(i0 - k, 0), min(i0 + k + 1, rows))
    c_ = np.arange(j0 - k, j0 + k + 1) % cols
    win = main_water.reshape(rows, cols)[np.ix_(r_, c_)]
    beside = ndimage.binary_dilation(win, structure=np.ones((3, 3), bool)) & ~win
    cells = (r_[:, None] * cols + c_[None, :])[beside]
    if cells.size == 0:
        raise SystemExit('no ground beside the main part near the throat')
    return int(cells[np.argmin(h_flat[cells])])


def governed_ops(state, gov, branches, salt_sea=GD.SALT_SEA, source=GD.SOURCE):
    """The governed routing's put-ins and take-outs (governed_drainage.govern's last round): the weir's share to the salt
    sea at Moscoviense's pass, the full seas' surplus at theirs, and each branch's share along its cuts."""
    entry, share, over = state['entry'], gov['share'], gov['over']
    ops = [(entry[source], share[salt_sea])] + [(entry[k], over[k]) for k in GD.FULL_SEAS]
    for target in (GD.HELD_BRANCH, GD.REST_BRANCH):
        unit_ops, _ = branches[target]
        ops += [(c, sign * share[target]) for c, sign in unit_ops]
    return ops


def inputs():
    """The records the routing follows from, with their hashes, and the climatology it reads (ignored by Git)."""
    frame_path, _ = FD.frame_record()
    return ({rel(p): sha256(p) for p in (FD.OUT, N.OUT, ROOT / RD.SURFACE, frame_path)},
            {rel(ROOT / FD.CLIMATOLOGY): sha256(ROOT / FD.CLIMATOLOGY)})


def route(log=print):
    """The design surface's routing: the governed routing on the closed budget, routed again with the near-side sea at
    its two levels and the main part's spill put in past the neck's crest. Returns the routing, the sea over each cell
    (an index into SEAS, -1 elsewhere) with the seas' levels, the neck's spill and the frame and climate it used."""
    started = time.time()
    fd = json.loads(FD.OUT.read_text())
    neck = json.loads(N.OUT.read_text())
    case = fd['cases'][BUDGET]
    frame_path, record = FD.frame_record()
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    surf = RD.load_surface(log)
    runoff, lake_net, _, _ = FD.load(case['evaporation'])
    geo = RD.geometry(surf, runoff, lake_net, log)
    level0, h = geo['level0'], geo['h']
    rows, cols = h.shape
    h_flat = h.ravel()
    area = np.repeat((geo['cell_km'] * 1e3) ** 2 * np.cos(np.radians(geo['lat'])), cols)
    state = RD.design_state(geo, runoff, lake_net, log=log)
    branches, _ = GD.branches_for(geo, state, rotation, log)
    gov = GD.govern(geo, state, runoff, lake_net, branches, log=log)
    log(f'governed shares {dict((k, round(v)) for k, v in gov["share"].items())}, {time.time() - started:.0f} s')

    # the near-side sea parted: each part's water at its steady level (near_side_neck.py)
    steady = neck['budgets'][BUDGET]['natural_neck']['steady']
    z_main, z_south = level0 + steady['main_level_m'], level0 + steady['south_level_m']
    st = N.setting(log)
    main_water = st['part']['main'] & (st['joined']['main'] < z_main)
    south_water = st['part']['south'] & (st['joined']['south'] < z_south)
    del st
    names, masks, levels = state['names'], state['masks'], state['levels']
    waters, sea_masks = [], []
    for k in names:
        if k == 'near_side':
            continue
        waters.append((masks[k], levels[k]))
        sea_masks.append((SEAS.index(k), masks[k], levels[k]))
    waters += [(main_water, z_main), (south_water, z_south)]
    sea_masks += [(0, main_water, z_main), (1, south_water, z_south)]
    names_r = [k for k in names if k != 'near_side'] + ['near_side_main', 'near_side_south']

    # the governed routing again over the parted sea, with the spill over the neck put in
    r0 = state['r']
    lake_full = ~r0['outlet'] & (h_flat < r0['level'][r0['label']])
    spill = spill_cell(h_flat, main_water, cols, neck['throat']['centre_today'])
    ops = governed_ops(state, gov, branches) + [(spill, steady['neck_flow_m3s'])]
    r, supply = GD.route_with(geo, runoff, lake_net, waters, names_r, ops, area, lake_full)
    log(f'routed over the parted sea: supplies {dict((k, round(v)) for k, v in supply.items())}, '
        f'{time.time() - started:.0f} s')

    # each sea over the cells below its level, in the order the globe surface drew them
    sea = np.full(h_flat.size, -1, np.int8)
    sea_levels = {}
    for i, mask, lev in sea_masks:
        sea[mask & (h_flat < lev)] = i
        sea_levels[SEAS[i]] = float(lev)
    return dict(r=r, sea=sea, sea_levels=sea_levels, level0=level0, rows=rows, cols=cols, area=area, h=h_flat,
                spill=spill, neck_flow=float(steady['neck_flow_m3s']), supply=supply, evaporation=case['evaporation'],
                frame=dict(record=rel(frame_path), sha256=sha256(frame_path), rotation_today_to_turned=rotation.tolist(),
                           state=record['frame']['state']))


def write(d, log=print):
    """The product and its record from route()'s routing."""
    started = time.time()
    r, sea, area, h = d['r'], d['sea'], d['area'], d['h']
    lake = r['lake'] & ~r['outlet'] & (sea < 0)
    cells = np.flatnonzero(lake).astype(np.int32)
    spills = ~r['closed'][r['label'][cells]]
    tracked, ignored = inputs()
    meta = dict(schema=PRODUCT_SCHEMA, grid=dict(pixels_per_degree=dr.PPD, rows=d['rows'], cols=d['cols']),
                common_level_m=d['level0'], seas=list(SEAS), levels_m=d['sea_levels'], upper_seas=list(UPPER_SEAS),
                neck=dict(spill_cell=d['spill'], flow_m3s=d['neck_flow']), frame=d['frame'],
                climate=dict(run=FD.RUN, years=FD.YEARS, climatology=FD.CLIMATOLOGY, budget=BUDGET,
                             open_water_evaporation_factor=d['evaporation']),
                inputs=tracked, ignored_inputs=ignored,
                reading_rule=('Cells index rows from 90 N and columns from 0 E, 16 per degree, on today\'s Moon: cell = '
                              'row * 5760 + column. receiver is the cell each cell drains to, discharge_m3s the mean '
                              'flow through it and drained_km2 the area upstream of it; sea is the index into seas of the '
                              'sea over the cell, -1 off the seas, each sea at levels_m above the dry Moon\'s geoid (the '
                              'common level common_level_m); outlet marks the routing\'s ends (the seas\' cells) and '
                              'routing_lake the cells its hollows\' water covers; lake_cell lists the cells under lakes '
                              'off the seas, with lake_level_m and lake_spills (whether the lake overflows). The upper '
                              'seas spill over their rims; the near-side sea is parted, neck.spill_cell the ground past '
                              'its crest where the main part\'s spill of neck.flow_m3s is put in.'),
                evidence=('The redrawn floor at the settled state, the seas at the levels the design climate\'s closed '
                          'budget gives them with the works of 16 November and the fixed weirs, and the routing of that '
                          'budget\'s rain over it with the near-side sea parted. Not validated against any observation.'))
    BULK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PRODUCT, metadata=json.dumps(meta), receiver=r['receiver'].astype(np.int32),
                        discharge_m3s=np.asarray(r['discharge'], np.float64),
                        drained_km2=(r['drained_m2'] / 1e6).astype(np.float32), sea=sea,
                        outlet=np.asarray(r['outlet'], bool), routing_lake=np.asarray(r['lake'], bool),
                        lake_cell=cells, lake_level_m=r['lake_level'][cells].astype(np.float32), lake_spills=spills)
    total = float(area.sum())
    q = r['discharge']
    into = np.flatnonzero((sea[r['receiver']] >= 0) & (sea < 0) & ~lake)
    summary = dict(
        seas={k: dict(level_above_common_m=round(v - d['level0'], 2),
                      area_km2=round(float(area[sea == SEAS.index(k)].sum()) / 1e6))
              for k, v in d['sea_levels'].items()},
        sea_share=round(float(area[sea >= 0].sum()) / total, 5),
        lake_cells=int(cells.size), lake_share=round(float(area[cells].sum()) / total, 5),
        closed_lake_share=round(float(area[cells[~spills]].sum()) / total, 5),
        rivers_into_the_seas_m3s=round(float(q[into].sum())),
        supplies_m3s={k: round(float(v)) for k, v in d['supply'].items()},
        neck_spill_m3s=round(d['neck_flow']))
    log(f'wrote {PRODUCT}, {time.time() - started:.0f} s')
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs=tracked, ignored_inputs=ignored, evidence=meta['evidence'],
        reading_rule=('Levels in m above the redrawn common level; areas in km2, flows in m3/s; shares of the Moon\'s '
                      'surface. product names the routing on the research drive, which design_routing.load() reads.'),
        summary=summary, product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)))


def stale():
    """What in the record no longer matches: the product, a tracked input or the code. Empty when current."""
    if not OUT.exists() or not PRODUCT.exists():
        return ['no routing yet']
    record = json.loads(OUT.read_text())
    tracked, ignored = inputs()
    out = [] if sha256(PRODUCT) == record['product']['sha256'] else [str(PRODUCT)]
    out += [k for k, v in {**tracked, **ignored}.items() if record['inputs'].get(k, record['ignored_inputs'].get(k)) != v]
    out += [f for f in CODE if record['producer']['files'].get(f) != sha256(ROOT / f)]
    return out


def load(log=print):
    """The routing product, once it matches its record and the record's inputs and code are current; arrays by name and
    its metadata under 'meta'."""
    behind = stale()
    if behind:
        raise SystemExit(f'the design routing is stale ({", ".join(behind)}); '
                         'run python -m research.studies.atlas_round.design_routing')
    with np.load(PRODUCT, allow_pickle=False) as z:
        out = {k: z[k] for k in z.files if k != 'metadata'}
        out['meta'] = json.loads(str(z['metadata']))
    log(f'the design routing from {PRODUCT}')
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    started = time.time()
    log = lambda *a: print(*a, flush=True)                                     # noqa: E731
    record = write(route(log), log)
    record['seconds'] = round(time.time() - started)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1) + '\n')
    log(f'wrote {args.out} in {record["seconds"]} s: {json.dumps(record["summary"])}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
