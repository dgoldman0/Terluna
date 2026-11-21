"""The design surface's routing on the redrawn floor: each cell's downstream neighbour, discharge and drained area, every
sea at its own level and the lakes, on the design climate's closed budget, the lakes filled by the corrected rule.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.design_routing   # about 6 minutes

The design surface is routed once for each climate and kept, so that the maps built on it read one routing: the globe
surface (globe_surface.py), the satoyama map (satoyama_final.py) and the climate model's surface for the run on the
design surface (design_turned_surface.py) read it, and later maps can. The routing:

- the weirs settled on the closed budget of the design climate (final_drainage.py's carried climatology and factor;
  lake_rule_balance.govern): the salt sea and Humboldtianum held at the common level, and Moscoviense, Orientale and
  Smythii-Marginis full at their rims;
- the near-side sea's parts, the routing and the levels set against each other (lake_rule_balance.settle_parts) until
  each part stands where its loss to the air meets the rivers that reach it: the main part over the neck's crest
  sending the southern part what the neck passes, on the ground past the crest while the southern part stands under
  it, or through the water once both stand over it;
- the lakes filled by the corrected rule, each hollow as one water (lake_rule.py), so that no lake stands above dry
  ground beside it, and the receivers joined through the lakes (hollow_courses.py), so that a river's course leads
  through its lakes to its sea.

Until 18 November the routing balanced each catchment by itself (research/studies/drainage_network/routing.py) and was
routed once with the parts at the levels near_side_neck.py gave.

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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'design_routing.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/design_routing')
PRODUCT = BULK / 'design_routing_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-design-routing/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-design-routing-grid/1'
CODE = ('research/studies/atlas_round/design_routing.py', 'research/studies/atlas_round/final_drainage.py',
        'research/studies/atlas_round/near_side_neck.py', 'research/studies/atlas_round/lake_rule_balance.py',
        'research/studies/atlas_round/lake_rule.py', 'research/studies/atlas_round/hollow_courses.py',
        'research/studies/atlas_round/design_works.py', 'research/studies/drainage_network/network.py',
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


def inputs():
    """The records the routing follows from, with their hashes, and the climatology it reads (ignored by Git)."""
    from research.studies.atlas_round import design_works as DW
    frame_path, _ = FD.frame_record()
    return ({rel(p): sha256(p) for p in (FD.OUT, N.OUT, DW.OUT, frame_path)},
            {rel(ROOT / FD.CLIMATOLOGY): sha256(ROOT / FD.CLIMATOLOGY), str(DW.PRODUCT): sha256(DW.PRODUCT)})


def route(log=print):
    """The design surface's routing: the weirs settled on the closed budget, then the near-side sea's parts, the
    routing and the levels set against each other, with the lakes filled by the corrected rule and the receivers
    joined through them. Returns the routing, the sea over each cell (an index into SEAS, -1 elsewhere) with the
    seas' levels, the neck's flow and where it is put in, the weirs' flows and the frame and climate it used."""
    from research.studies.atlas_round import design_works as DW
    from research.studies.atlas_round import hollow_courses as HC
    from research.studies.atlas_round import lake_rule as LR
    from research.studies.atlas_round import lake_rule_balance as LB
    started = time.time()
    fd = json.loads(FD.OUT.read_text())
    neck_record = json.loads(N.OUT.read_text())
    case = fd['cases'][BUDGET]
    frame_path, record = FD.frame_record()
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    d, st = DW.load(log)
    level0, h = d['level0'], d['h']
    rows, cols = h.shape
    h_flat = h.ravel()
    names = list(LB.WATERS)
    runoff, lake_net, _, _ = FD.load(case['evaporation'])
    # the weirs, begun from the shares final_drainage.py settled on the same climate
    begun = dict(share=case['weir_history'][-1]['shares_m3s'], over=case['full_seas_pass_on_m3s'])
    g = LB.govern(d, st, runoff, lake_net, start=begun, log=log)
    gov = g['gov']
    log(f'governed shares {dict((k, round(v)) for k, v in gov["share"].items())}, {time.time() - started:.0f} s')
    area = (((d['cell_km'] * 1000.0) ** 2 * np.maximum(np.cos(np.radians(d['lat'])), 1e-6))[:, None] * np.ones((1, cols))).ravel()
    net = lake_net.ravel() * dr.MM_DAY * area
    # the near-side sea's parts, the routing and the levels set against each other
    neck = LB.neck_from_record(st, neck_record)
    tab = LB.part_tables(d, st, -net, area)
    rivers = dict(main=float(g['supply'][LB.MAIN]), south=float(g['supply'][LB.SOUTH]))
    parts = LB.settle_parts(d, st, tab, neck, neck_record['throat']['centre_today'], runoff, lake_net, net, area,
                            g['fixed'], g['branches'], rivers, log=log)
    r = LR.hollow_routing(parts['s'], parts['terms'], parts['state'])
    c = HC.connect(parts['s'], parts['terms'], parts['state'], r, h.shape)
    water, levels, z = parts['water'], parts['levels'], parts['z']
    gathered = np.array([float(c['discharge'][water == i].sum()) for i in range(len(names))])
    if np.abs(gathered - parts['got']).max() > 1e-3:
        raise SystemExit('the rivers gathered along the joined courses do not bring the seas what the hollows give them')
    b = r['budget']
    gap = b['runoff_of_dry_ground'] + b['put_in'] + b['lakes_rain_less_evaporation'] - b['reaches_the_seas'] - b['kept_by_lakes_not_full']
    if abs(gap) > 1.0 or float(r['left'].min()) < -1e-6:
        raise SystemExit(f'the routing\'s budget does not close ({gap:.3g} m3/s) or a cut takes more than its lake holds')
    log(f'routed with the near-side sea at {z["main"] - level0:.2f} and {z["south"] - level0:.2f} m: supplies '
        f'{dict((k, round(float(v))) for k, v in zip(names, parts["got"]))}, {time.time() - started:.0f} s')
    joined = bool(z['south'] > neck.crest)               # both parts over the crest: one sea across the neck
    spill = None
    if z['neck_flow'] > 0.0 and not joined:
        spill = spill_cell(h_flat, water == LB.MAIN, cols, neck_record['throat']['centre_today'])
    sea = np.full(h_flat.size, -1, np.int8)
    sea_levels = {}
    for i, k in enumerate(names):
        sea[water == i] = SEAS.index(k)
        sea_levels[k] = float(levels[i])
    loss = {k: -float(net[water == i].sum()) for i, k in enumerate(names)}
    own = {k: round(float(v)) for k, v in zip(names, gov['own'])}
    weirs = dict(moscoviense_surplus_m3s=round(gov['surplus']), released_m3s={k: round(v) for k, v in gov['share'].items()},
                 arrived_m3s={k: round(float(g['supply'][names.index(k)] - gov['own'][names.index(k)]))
                              for k in (GD.SALT_SEA, GD.HELD_BRANCH, GD.REST_BRANCH)},
                 full_seas_pass_on_m3s={k: round(v) for k, v in gov['over'].items()},
                 held_seas_left_short_m3s={k: round(v) for k, v in gov['unmet'].items()},
                 branch_flows_m3s={k: [round(q, 1) for q in v] for k, v in zip(('humboldtianum_branch', 'smythii_marginis_branch'), parts['flows'])},
                 own_rivers_at_the_common_level_m3s=own)
    return dict(r=dict(receiver=c['receiver'], discharge=c['discharge'], drained_m2=c['drained_m2'], lake=r['lake'],
                       lake_level=r['lake_level'], lake_spills=r['lake_spills'], outlet=r['outlet'], crossing=c['crossing'],
                       keeps=~c['passes']),
                sea=sea, sea_levels=sea_levels, level0=level0, rows=rows, cols=cols, area=area, h=h_flat, spill=spill,
                neck_flow=float(z['neck_flow']), joined=joined, supply=dict(zip(names, (float(v) for v in parts['got']))),
                loss=loss, weirs=weirs, budget={k: round(v, 1) for k, v in b.items()},
                parts=dict(history=parts['history'], sills=parts['sills'],
                           left_out_of_balance_m3s={k: round(v, 1) for k, v in parts['short'].items()},
                           spills_over_a_sill_m3s={k: round(v, 1) for k, v in parts['spills'].items()},
                           crest_above_common_m=round(neck.crest - level0, 2)),
                evaporation=case['evaporation'],
                frame=dict(record=rel(frame_path), sha256=sha256(frame_path), rotation_today_to_turned=rotation.tolist(),
                           state=record['frame']['state']))


def write(d, log=print):
    """The product and its record from route()'s routing."""
    started = time.time()
    r, sea, area, h = d['r'], d['sea'], d['area'], d['h']
    lake = r['lake'] & ~r['outlet'] & (sea < 0)
    cells = np.flatnonzero(lake).astype(np.int32)
    spills = r['lake_spills'][cells]
    tracked, ignored = inputs()
    meta = dict(schema=PRODUCT_SCHEMA, grid=dict(pixels_per_degree=dr.PPD, rows=d['rows'], cols=d['cols']),
                common_level_m=d['level0'], seas=list(SEAS), levels_m=d['sea_levels'], upper_seas=list(UPPER_SEAS),
                neck=dict(spill_cell=d['spill'], flow_m3s=d['neck_flow'], joined=d['joined']), frame=d['frame'],
                lake_rule='each hollow filled as one water (lake_rule.py); receivers joined through the lakes (hollow_courses.py)',
                climate=dict(run=FD.RUN, years=FD.YEARS, climatology=FD.CLIMATOLOGY, budget=BUDGET,
                             open_water_evaporation_factor=d['evaporation']),
                inputs=tracked, ignored_inputs=ignored,
                reading_rule=('Cells index rows from 90 N and columns from 0 E, 16 per degree, on today\'s Moon: cell = '
                              'row * 5760 + column. receiver is the cell each cell drains to, discharge_m3s the mean '
                              'flow through it and drained_km2 the area upstream of it; sea is the index into seas of the '
                              'sea over the cell, -1 off the seas, each sea at levels_m above the dry Moon\'s geoid (the '
                              'common level common_level_m); outlet marks the routing\'s ends (the seas\' cells) and '
                              'routing_lake the cells its hollows\' water covers; lake_cell lists the cells under lakes '
                              'off the seas, with lake_level_m and lake_spills (whether the lake overflows). A lowest '
                              'point\'s receiver is the cell its lake\'s water goes on to, so the receivers lead from any '
                              'cell to a sea; under a lake the discharge is what has come in on the way to its exit, '
                              'and crossing marks the dry cells of a sill that a lake\'s water is carried across; keeps '
                              'marks the lowest points whose lakes pass nothing on, where a course\'s water ends. The '
                              'upper seas spill over their rims. neck.flow_m3s is what the near-side sea\'s main part '
                              'sends its southern part: where neck.joined is false the southern part stands under the '
                              'neck\'s crest and the flow is put in on the ground past it at neck.spill_cell; where it '
                              'is true both parts stand over the crest, one sea across the neck, and the cell is null.'),
                evidence=('The redrawn floor at the settled state, the seas at the levels the design climate\'s closed '
                          'budget gives them with the works of 16 November and the fixed weirs, and the routing of that '
                          'budget\'s rain over it, the lakes filled by the corrected rule. Not validated against any '
                          'observation.'))
    BULK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PRODUCT, metadata=json.dumps(meta), receiver=r['receiver'].astype(np.int32),
                        discharge_m3s=np.asarray(r['discharge'], np.float64),
                        drained_km2=(r['drained_m2'] / 1e6).astype(np.float32), sea=sea,
                        outlet=np.asarray(r['outlet'], bool), routing_lake=np.asarray(r['lake'], bool),
                        crossing=np.asarray(r['crossing'], bool), keeps=np.asarray(r['keeps'], bool),
                        lake_cell=cells, lake_level_m=r['lake_level'][cells].astype(np.float32), lake_spills=spills)
    total = float(area.sum())
    q = r['discharge']
    into = np.flatnonzero((sea[r['receiver']] >= 0) & (sea < 0) & ~lake)
    summary = dict(
        seas={k: dict(level_above_common_m=round(v - d['level0'], 2),
                      area_km2=round(float(area[sea == SEAS.index(k)].sum()) / 1e6),
                      loses_m3s=round(d['loss'][k]), receives_m3s=round(d['supply'][k]))
              for k, v in d['sea_levels'].items()},
        sea_share=round(float(area[sea >= 0].sum()) / total, 5),
        lake_cells=int(cells.size), lake_share=round(float(area[cells].sum()) / total, 5),
        closed_lake_share=round(float(area[cells[~spills]].sum()) / total, 5),
        rivers_into_the_seas_m3s=round(float(q[into].sum())),
        supplies_m3s={k: round(float(v)) for k, v in d['supply'].items()},
        neck_spill_m3s=round(d['neck_flow']), near_side_sea_joined_across_its_neck=d['joined'],
        open_water_evaporation_factor=d['evaporation'], weirs=d['weirs'], near_side=d['parts'], budget_m3s=d['budget'])
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
