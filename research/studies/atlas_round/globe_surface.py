"""The design surface the globe draws on the redrawn floor: the settled heights, every water body at its own level and
the rivers, on the design climate's closed budget with the near-side sea parted at its neck.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.globe_surface   # about 10 minutes

The design climate on the final surface is the run at sunlight 0.930 (register, 17 November), its closed budget the
one on which the near-side sea parts. The surface holds:

- the heights: LOLA with the settled case's U - N, as redrawn_drainage.load_surface gives them;
- the seas at their own levels: Moscoviense, Orientale and Smythii-Marginis full at their rims; the salt sea and
  Humboldtianum at the common level, where the fixed weirs, sized on this budget, hold them; and the near-side sea
  parted (near_side_neck.py's steady state), its main part 2.1 m over the neck's crest and Cognitum, Humorum, Nubium and
  southern Procellarum about 100 m lower;
- the lakes and rivers of the governed routing (final_drainage.py's works and shares) routed again with the near-side
  sea at its two levels, so that rivers run on over the floor the parted sea leaves dry, and with the main part's spill
  put in past the neck's crest, so that the river over the neck runs down to the southern seas.

It is written to the research drive in geography's drainage product's form (cells from 90 N and 0 E, 16 per degree),
with the heights, each water cell's level and body, and the channels, and summarised in results/globe_surface.json.
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
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'globe_surface.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/globe_surface')
PRODUCT = BULK / 'globe_surface_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-globe-surface/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-globe-surface-grid/1'
CODE = ('research/studies/atlas_round/globe_surface.py', 'research/studies/atlas_round/final_drainage.py',
        'research/studies/atlas_round/near_side_neck.py', 'research/studies/atlas_round/governed_drainage.py',
        'research/studies/atlas_round/redrawn_drainage.py', 'research/studies/drainage_network/network.py')
BUDGET = 'closed_budget'
CHANNEL_M3S = 1.0                    # discharge kept, as geography's drainage product keeps it
BODIES = ('lake', 'near-side sea', 'southern near-side seas', 'South Pole–Aitken Sea', 'Mare Humboldtianum',
          'Mare Moscoviense', 'Mare Orientale', 'Smythii–Marginis Sea')
BODY_OF = {'south_pole_aitken': 3, 'humboldtianum': 4, 'moscoviense': 5, 'orientale': 6, 'smythii_marginis': 7}
SPILL_SEARCH_DEG = 0.6               # the spill is put in at the lowest ground beside the main part this near the throat


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


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


def build(log=print):
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
    waters, body_masks = [], []
    for k in names:
        if k == 'near_side':
            continue
        waters.append((masks[k], levels[k]))
        body_masks.append((BODY_OF[k], masks[k], levels[k]))
    waters += [(main_water, z_main), (south_water, z_south)]
    body_masks += [(1, main_water, z_main), (2, south_water, z_south)]
    names_r = [k for k in names if k != 'near_side'] + ['near_side_main', 'near_side_south']

    # the governed routing again over the parted sea, with the spill over the neck put in
    r0 = state['r']
    lake_full = ~r0['outlet'] & (h_flat < r0['level'][r0['label']])
    spill = spill_cell(h_flat, main_water, cols, neck['throat']['centre_today'])
    ops = governed_ops(state, gov, branches) + [(spill, steady['neck_flow_m3s'])]
    r, supply = GD.route_with(geo, runoff, lake_net, waters, names_r, ops, area, lake_full)
    log(f'routed over the parted sea: supplies {dict((k, round(v)) for k, v in supply.items())}, {time.time() - started:.0f} s')

    # the water: each sea's cells at its level, then the lakes off the seas
    level = np.full(h_flat.size, np.nan, np.float32)
    body = np.zeros(h_flat.size, np.int8)
    for b, mask, lev in body_masks:
        cells = np.flatnonzero(mask & (h_flat < lev))
        level[cells] = lev
        body[cells] = b
    lake = r['lake'] & ~r['outlet'] & np.isnan(level)
    level[lake] = r['lake_level'][lake]
    wet = ~np.isnan(level)
    water_cell = np.flatnonzero(wet).astype(np.int32)
    q = r['discharge']
    keep = np.flatnonzero((q >= CHANNEL_M3S) & ~r['outlet'] & ~r['lake'] & ~wet).astype(np.int32)

    # the neck's river, from the spill down the routing
    path = RD.walk(r, spill)
    surface = np.where(wet, level, h_flat)
    neck_river = NW.long_profile(path, surface, q, geo['lat'], cols, geo['cell_km'], level0, r['lake'])
    neck_river.update(starts=NW.place(spill, cols, rotation), reaches=('southern near-side seas' if body[path[-1]] == 2
                                                                      else f'body {int(body[path[-1]])}'))

    # the summary: shares, each sea's area, depth and volume, the lakes, the rivers into the seas
    total = float(area.sum())
    depth = np.where(wet, np.nan_to_num(level, nan=0.0) - h_flat, 0.0)
    seas = {}
    for b in range(1, len(BODIES)):
        m = body == b
        a = float(area[m].sum())
        seas[BODIES[b]] = dict(level_m=round(float(level[m][0] - level0), 2) if m.any() else None,
                               area_km2=round(a / 1e6), share=round(a / total, 5),
                               mean_depth_m=round(float((depth[m] * area[m]).sum() / max(a, 1.0))),
                               greatest_depth_m=round(float(depth[m].max())) if m.any() else 0,
                               volume_km3=round(float((depth[m] * area[m]).sum()) / 1e9))
    lake_labels, n_lakes = ndimage.label((body == 0).reshape(rows, cols) & wet.reshape(rows, cols),
                                         structure=np.ones((3, 3), bool))
    lake_area = float(area[(body == 0) & wet].sum())
    into = np.flatnonzero((body[r['receiver']] > 0) & (body == 0) & ~wet)
    summary = dict(
        common_level_above_the_dry_geoid_m=round(level0, 2),
        water_share=round(float(area[wet].sum()) / total, 5), sea_share=round(float(area[body > 0].sum()) / total, 5),
        lake_share=round(lake_area / total, 5), lakes=int(n_lakes),
        lake_volume_km3=round(float((depth[(body == 0) & wet] * area[(body == 0) & wet]).sum()) / 1e9),
        rivers_into_the_seas_m3s=round(float(q[into].sum())), seas=seas,
        neck_spill_m3s=round(steady['neck_flow_m3s']), neck_river=neck_river,
        deepest_water_m=round(float(depth.max())), highest_ground_m=round(float(h_flat.max() - level0)),
        lowest_ground_m=round(float(h_flat.min() - level0)))
    log(f'summary {json.dumps(summary)[:1500]}, {time.time() - started:.0f} s')

    BULK.mkdir(parents=True, exist_ok=True)
    meta = dict(schema=PRODUCT_SCHEMA, grid=dict(pixels_per_degree=dr.PPD, rows=rows, cols=cols),
                common_level_m=level0, bodies=list(BODIES),
                frame=dict(record=str(frame_path.relative_to(ROOT)), sha256=sha256(frame_path),
                           rotation_today_to_turned=rotation.tolist(), state=record['frame']['state'],
                           pole_move_deg=record['frame'].get('pole_move_deg'),
                           sub_earth_point=record['frame'].get('sub_earth_point')),
                climate=dict(run=FD.RUN, years=FD.YEARS, budget=BUDGET, open_water_evaporation_factor=case['evaporation']),
                reading_rule=('Cells index rows from 90 N and columns from 0 E, 16 per degree, on today\'s Moon: cell = '
                              'row * 5760 + column. height_m is the ground above the dry Moon\'s geoid as redrawn at the '
                              'settled state; water_cell lists every cell under water, water_level_m its surface above '
                              'the same geoid and water_body its body (an index into bodies, 0 a lake off the seas); '
                              'channel_discharge_m3s is the mean discharge through channel_cell, on dry ground.'),
                evidence=('The redrawn floor at the settled state, the seas at the levels the design climate\'s closed '
                          'budget gives them with the works of 16 November and the fixed weirs, and the routing of that '
                          'budget\'s rain over it with the near-side sea parted. Not validated against any observation.'))
    np.savez_compressed(PRODUCT, metadata=json.dumps(meta), height_m=h.astype(np.float32), water_cell=water_cell,
                        water_level_m=level[water_cell].astype(np.float32), water_body=body[water_cell],
                        channel_cell=keep, channel_discharge_m3s=q[keep].astype(np.float32))
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (FD.OUT, N.OUT, ROOT / RD.SURFACE, frame_path)},
        ignored_inputs={rel(ROOT / FD.CLIMATOLOGY): sha256(ROOT / FD.CLIMATOLOGY)},
        evidence=meta['evidence'],
        reading_rule=('Levels and heights in m above the redrawn common level; areas in km2, volumes in km3, flows in '
                      'm3/s; shares of the Moon\'s surface. seas are the seven bodies at their own levels; lakes are '
                      'the routing\'s lakes off the seas; neck_river is the main part\'s spill run down the routing from '
                      'the ground past the crest. product names the grid the globe reads, on the research drive.'),
        summary=summary, product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)),
        seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    record = build(log=lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1) + '\n')
    print(f'wrote {args.out} in {record["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
