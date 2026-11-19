"""The design surface the globe draws on the redrawn floor: the settled heights, every water body at its own level and
the rivers, on the design climate's closed budget with the near-side sea parted at its neck.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.globe_surface   # about a minute

The design climate is final_drainage.py's (register, 18 November: the run on the design surface carried by pattern
scaling to 300.31 K), its closed budget the one on which the near-side sea parts. The surface holds:

- the heights: LOLA with the settled case's U - N, as redrawn_drainage.load_surface gives them;
- the seas at their own levels: Moscoviense, Orientale and Smythii-Marginis full at their rims; the salt sea and
  Humboldtianum at the common level, where the fixed weirs, sized on this budget, hold them; and the near-side sea
  parted (near_side_neck.py's steady state), its main part over the neck's crest and Cognitum, Humorum, Nubium and
  southern Procellarum lower;
- the lakes and rivers of the design routing, the governed routing routed again with the near-side sea at its two
  levels, so that rivers run on over the floor the parted sea leaves dry, and with the main part's spill put in past the
  neck's crest, so that the river over the neck runs down to the southern seas;
- the works' cuts (register, 16 November; design_cuts.py): the trenches from Moscoviense towards Humboldtianum and
  Smythii-Marginis,
  whose ground is taken down to the water they carry, with the branches' flow drawn through them as a river. The
  routing puts that flow in past each trench; the trench's own channel, a few metres to tens of metres of water, is not
  cut below the water's level.

The routing is design_routing.py's product, which this module reads after design_routing.load() has checked it against
its record; the routing is computed there once for each climate and nothing is routed here. The surface is written to
the research drive in geography's drainage product's form (cells from 90 N and 0 E, 16 per degree), with the heights,
each water cell's level and body, and the channels, and summarised in results/globe_surface.json.
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
from research.studies.atlas_round import design_cuts as DC
from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'globe_surface.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/globe_surface')
PRODUCT = BULK / 'globe_surface_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-globe-surface/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-globe-surface-grid/1'
CODE = ('research/studies/atlas_round/globe_surface.py', 'research/studies/atlas_round/design_routing.py',
        'research/studies/atlas_round/design_cuts.py',
        'research/studies/atlas_round/redrawn_drainage.py', 'research/studies/drainage_network/network.py')
BUDGET = DR.BUDGET
CHANNEL_M3S = 1.0                    # discharge kept, as geography's drainage product keeps it
BODIES = ('lake', 'near-side sea', 'southern near-side seas', 'South Pole–Aitken Sea', 'Mare Humboldtianum',
          'Mare Moscoviense', 'Mare Orientale', 'Smythii–Marginis Sea')
BODY_OF = {'south_pole_aitken': 3, 'humboldtianum': 4, 'moscoviense': 5, 'orientale': 6, 'smythii_marginis': 7}


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(log=print):
    started = time.time()
    routing = DR.load(log)
    meta = routing['meta']
    frame_path, record = FD.frame_record()
    if meta['frame']['sha256'] != sha256(frame_path):
        raise SystemExit(f'the design routing was made in another frame than {frame_path.name}')
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    surf = RD.load_surface(log)
    h, level0, lat, cell_km = surf['h'], surf['level0'], surf['lat'], surf['cell_km']
    if abs(level0 - meta['common_level_m']) > 1e-6:
        raise SystemExit('the design routing\'s common level is not the surface\'s')
    rows, cols = h.shape
    h_flat = h.ravel()
    area = np.repeat((cell_km * 1e3) ** 2 * np.cos(np.radians(lat)), cols)
    r = dict(receiver=routing['receiver'], outlet=routing['outlet'], lake=routing['routing_lake'])
    q = routing['discharge_m3s']

    # the water: each sea's cells at its level, then the lakes off the seas
    sea = routing['sea']
    on_sea = sea >= 0
    level = np.full(h_flat.size, np.nan, np.float32)
    body = np.zeros(h_flat.size, np.int8)
    level[on_sea] = np.array([meta['levels_m'][k] for k in meta['seas']], np.float32)[sea[on_sea]]
    body[on_sea] = sea[on_sea] + 1
    lake = routing['lake_cell']
    level[lake] = routing['lake_level_m']
    spills = np.zeros(h_flat.size, bool)
    spills[lake] = routing['lake_spills']
    for k in meta['upper_seas']:                      # the upper seas spill over their rims; the lower seas are terminal
        spills[body == meta['seas'].index(k) + 1] = True
    wet = ~np.isnan(level)
    water_cell = np.flatnonzero(wet).astype(np.int32)

    # the works' cuts: each trench's ground taken down to the water it carries; a lake standing higher that the trench
    # crosses drains to that water, its floor below staying wet and the rest drying; the branches' flow runs through
    cuts = DC.load(log)
    cut, cut_level = cuts['cut_cell'], cuts['cut_level_m']
    cut_depth = np.maximum(h_flat[cut] - cut_level, 0.0)
    h_flat[cut] = np.minimum(h_flat[cut], cut_level)
    crossed = wet[cut] & (body[cut] == 0) & (np.nan_to_num(level[cut], nan=-np.inf) > cut_level)
    drained_lakes, drained_km2 = 0, 0.0
    if crossed.any():
        lakes_lab, _ = ndimage.label((wet & (body == 0)).reshape(rows, cols), structure=np.ones((3, 3), bool))
        lakes_lab = lakes_lab.ravel()
        for lab_ in np.unique(lakes_lab[cut[crossed]]):
            cells = np.flatnonzero(lakes_lab == lab_)
            new = float(cut_level[crossed & (lakes_lab[cut] == lab_)].min())
            drained_lakes += 1
            drained_km2 += float(area[cells[h_flat[cells] >= new]].sum()) / 1e6
            level[cells] = np.where(h_flat[cells] < new, np.float32(new), np.float32(np.nan))
        wet = ~np.isnan(level)
        water_cell = np.flatnonzero(wet).astype(np.int32)
    in_cut = np.zeros(h_flat.size, bool)
    in_cut[cut] = True
    cut_pieces = []                                  # each piece of trench at one water level: its deepest point
    for lev in np.unique(cut_level):
        piece = cut_level == lev
        k = int(np.argmax(np.where(piece, cut_depth, -1.0)))
        cut_pieces.append(dict(level_m=round(float(lev - level0), 1), deepest_m=round(float(cut_depth[k]), 1),
                               cells=int(piece.sum()), km=round(float(np.sqrt(area[cut[piece]]).sum()) / 1e3),
                               flow_m3s=round(float(cuts['cut_flow_m3s'][piece].max())),
                               at=NW.place(int(cut[k]), cols, rotation)))
    cut_pieces.sort(key=lambda x: -x['deepest_m'])
    q = q.copy()
    q[cut] = np.maximum(q[cut], cuts['cut_flow_m3s'])
    keep = np.flatnonzero((q >= CHANNEL_M3S) & ~r['outlet'] & (~r['lake'] | in_cut) & ~wet).astype(np.int32)

    # the neck's river, from the spill down the routing; none where both parts of the near-side sea stand over the
    # neck's crest and the flow between them passes through the sea's water
    spill = meta['neck']['spill_cell']
    neck_river = None
    if spill is not None:
        path = RD.walk(r, spill)
        surface = np.where(wet, level, h_flat)
        neck_river = NW.long_profile(path, surface, q, lat, cols, cell_km, level0, r['lake'])
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
        closed_lake_share=round(float(area[(body == 0) & wet & ~spills].sum()) / total, 5),
        lake_volume_km3=round(float((depth[(body == 0) & wet] * area[(body == 0) & wet]).sum()) / 1e9),
        rivers_into_the_seas_m3s=round(float(routing['discharge_m3s'][into].sum())), seas=seas,   # the cuts' flow leaves a sea
        neck_spill_m3s=round(meta['neck']['flow_m3s']), neck_river=neck_river,
        cuts=dict(trench_cells=int(cut.size), deepest_m=round(float(cut_depth.max()), 1) if cut.size else 0.0,
                  ground_taken_km3=round(float((cut_depth * area[cut]).sum()) / 1e9, 1),
                  lakes_drained=drained_lakes, lake_drained_km2=round(drained_km2),
                  pieces=cut_pieces,
                  flow_m3s=round(float(cuts['cut_flow_m3s'].max())) if cut.size else 0),
        deepest_water_m=round(float(depth.max())), highest_ground_m=round(float(h_flat.max() - level0)),
        lowest_ground_m=round(float(h_flat.min() - level0)))
    log(f'summary {json.dumps(summary)[:1500]}, {time.time() - started:.0f} s')

    BULK.mkdir(parents=True, exist_ok=True)
    climate = meta['climate']
    product_meta = dict(
        schema=PRODUCT_SCHEMA, grid=dict(pixels_per_degree=dr.PPD, rows=rows, cols=cols),
        common_level_m=level0, bodies=list(BODIES),
        frame=dict(record=str(frame_path.relative_to(ROOT)), sha256=sha256(frame_path),
                   rotation_today_to_turned=rotation.tolist(), state=record['frame']['state'],
                   pole_move_deg=record['frame'].get('pole_move_deg'),
                   sub_earth_point=record['frame'].get('sub_earth_point')),
        climate=dict(run=climate['run'], years=climate['years'], budget=climate['budget'],
                     open_water_evaporation_factor=climate['open_water_evaporation_factor']),
        reading_rule=('Cells index rows from 90 N and columns from 0 E, 16 per degree, on today\'s Moon: cell = '
                      'row * 5760 + column. height_m is the ground above the dry Moon\'s geoid as redrawn at the '
                      'settled state; water_cell lists every cell under water, water_level_m its surface above '
                      'the same geoid, water_body its body (an index into bodies, 0 a lake off the seas) and '
                      'water_spills whether its water overflows (a lake that spills, an upper sea at its rim); '
                      'channel_discharge_m3s is the mean discharge through channel_cell, on dry ground. The heights '
                      'carry the works\' cuts, taken down to the water their trenches carry.'),
        evidence=('The redrawn floor at the settled state, the seas at the levels the design climate\'s closed '
                  'budget gives them with the works of 16 November and the fixed weirs, and the routing of that '
                  'budget\'s rain over it with the near-side sea parted. Not validated against any observation.'))
    np.savez_compressed(PRODUCT, metadata=json.dumps(product_meta), height_m=h.astype(np.float32),
                        water_cell=water_cell, water_level_m=level[water_cell].astype(np.float32),
                        water_body=body[water_cell], water_spills=spills[water_cell],
                        channel_cell=keep, channel_discharge_m3s=q[keep].astype(np.float32))
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)   # noqa: E731
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (DR.OUT, DC.OUT, ROOT / RD.SURFACE, frame_path)},
        ignored_inputs={str(DR.PRODUCT): sha256(DR.PRODUCT), str(DC.PRODUCT): sha256(DC.PRODUCT)},
        evidence=product_meta['evidence'],
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
