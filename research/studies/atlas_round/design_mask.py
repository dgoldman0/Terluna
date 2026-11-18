"""The land mask the climate model would run on for a routing's seas and lakes, in the frame of the design climate's run.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python3 -m research.studies.atlas_round.design_mask   # under a minute

The design climate is a run on the design surface: the design routing's seas and lakes turned into the Moon's settled
frame and averaged onto the model's T21 grid, where the cells with the least land become sea until the sea covers the
surface's own water share (design_turned_surface.py, turned_surface.py, exoplasim_run.land_mask). When the routing's
waters change, the mask the model would run on changes with them. This module makes that mask for a routing's waters
by the same steps, in the frame the design climate's run stands in, without writing the model's surface files, and
sets it against the mask the run held. climate.gcm.surface_change carries the run's climate onto it.

The frame is the run's own, so that the mask lies on the grid the run's fields lie on. The settled turn under the new
waters differs from it by the weight the waters move; design_settled_turn.py gives that turn when a run is made.

By default the waters are the design routing's (design_routing.load()); --product takes another routing in the form
lake_rule_balance.py keeps. The record, results/design_mask.json, holds the mask and each cell's land fraction.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from climate.gcm import turned
from research.studies.atlas_round import design_routing as DRt
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import turned_surface as TS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'design_mask.json'
SCHEMA = 'terluna.research.atlas-round-design-mask/1'
CODE = ('research/studies/atlas_round/design_mask.py', 'research/studies/atlas_round/turned_surface.py',
        'climate/gcm/turned.py', 'climate/gcm/boundary.py')
GRID = 'gaussian_T21'


def rel(path) -> str:
    return str(Path(path).relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path)


def waters(product=None, log=print):
    """A routing's waters on the 16 px/deg grid: the sea over each cell (an index into levels, -1 off the seas), the
    seas' levels, the lakes' cells and levels, and what the routing is. The design routing's by default; a path names
    a routing in lake_rule_balance.py's form."""
    if product is None:
        r = DRt.load(log)
        meta = r['meta']
        return dict(sea=r['sea'], levels=[meta['levels_m'][k] for k in meta['seas']], lake_cell=r['lake_cell'],
                    lake_level=r['lake_level_m'].astype(float), common_level=float(meta['common_level_m']),
                    source=dict(path=str(DRt.PRODUCT), sha256=DRt.sha256(DRt.PRODUCT)))
    with np.load(product, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        return dict(sea=z['water'], levels=list(meta['levels_m']), lake_cell=z['lake_cell'],
                    lake_level=z['lake_level_m'].astype(float), common_level=float(meta['common_level_m']),
                    source=dict(path=str(product), sha256=RD.sha256(product)))


def surface(h, level0, w):
    """The waters on the redrawn floor in boundary.lakes_product's form on the 4 px/deg grid, as
    design_turned_surface.lakes_surface makes it, with the seas' fraction."""
    c = TS.coarse(h, w['sea'].reshape(h.shape), w['lake_cell'], w['lake_level'], w['levels'])
    rows, cols = c['water_fraction'].shape
    lat = 90.0 - (np.arange(rows) + 0.5) / (RD.dr.PPD / TS.FACTOR)
    lon = (np.arange(cols) + 0.5) / (RD.dr.PPD / TS.FACTOR)
    weights = np.cos(np.radians(lat))[:, None] * np.ones(cols)
    share = lambda a: float((a * weights).sum() / weights.sum())                # noqa: E731
    wet = c['water_fraction']
    prod = dict(lat_deg=lat, lon_deg=lon, water_fraction=wet, mean_height_m=c['mean_height_m'],
                mean_water_depth_m=c['mean_water_depth_m'],
                water_surface_m=np.where(wet > 0, c['water_surface_m'], level0), meta=dict(level_m=level0),
                sha256=w['source']['sha256'][:16], drainage_sha256=w['source']['sha256'][:16], routing_run=FD.RUN,
                sea_share=share(c['sea_fraction']), lake_share=share(wet - c['sea_fraction']))
    return prod, c['sea_fraction']


def mask_in_frame(prod, sea, rotation):
    """The T21 land fraction of a surface carried into a frame, and PlaSim's binary mask from it (1 land, 0 sea)."""
    moved = turned.turned_product(prod, sea, np.asarray(rotation, float))
    fields = TS.regrid_signed(moved, GRID)
    summary = TS.t21_summary(fields)
    return fields, summary, moved


def build(product=None, log=print):
    started = time.time()
    surf = RD.load_surface(log)
    h, level0 = surf['h'], float(surf['level0'])
    w = waters(product, log)
    if abs(level0 - w['common_level']) > RD.RIM_CHECK_M:
        raise SystemExit('the routing\'s common level is not the redrawn floor\'s')
    frame_path, frame = FD.frame_record()
    rotation = frame['frame']['rotation_today_to_turned']
    prod, sea = surface(h, level0, w)
    fields, summary, moved = mask_in_frame(prod, sea, rotation)
    with np.load(ROOT / FD.TURNED_CLIMATOLOGY) as z:
        lat, lon, run_mask = np.asarray(z['lat'], float), np.asarray(z['lon'], float), np.asarray(z['lsm'], float) > 0.5
    if not (np.allclose(fields['lat'], lat, atol=1e-3) and np.allclose(fields['lon'] % 360.0, lon % 360.0, atol=1e-3)):
        raise SystemExit('the model\'s grid is not the climatology\'s')
    mask = summary['mask'] > 0.5
    area = np.asarray(fields['cell_area_fraction'], float)
    share = lambda m: round(float(area[m].sum() / area.sum()), 5)                # noqa: E731
    log(f'water {prod["sea_share"] + prod["lake_share"]:.5f} of the surface today (seas {prod["sea_share"]:.5f}), T21 '
        f'water {summary["water_share"]:.4f}; against the run\'s mask {int((run_mask & ~mask).sum())} cells to water '
        f'and {int((~run_mask & mask).sum())} to land, {time.time() - started:.0f} s')
    return dict(
        schema=SCHEMA, producer=dict(files={f: RD.sha256(ROOT / f) for f in CODE}),
        run=FD.RUN, years=FD.YEARS, frame=dict(record=rel(frame_path), sha256=RD.sha256(frame_path)),
        routing=w['source'], inputs={rel(ROOT / RD.SURFACE): RD.sha256(ROOT / RD.SURFACE)},
        ignored_inputs={FD.TURNED_CLIMATOLOGY: RD.sha256(ROOT / FD.TURNED_CLIMATOLOGY)},
        evidence=('A calculation: the routing\'s seas and lakes on the redrawn floor averaged onto the atlas\'s 4 '
                  'pixel/degree grid, carried into the frame the design climate\'s run stands in and averaged onto the '
                  'model\'s T21 grid, with PlaSim\'s binary mask made from the land fractions as the runner makes it. '
                  'No surface files are written and no run is made on it.'),
        reading_rule=('lat (north first) and lon (east) are the T21 grid\'s in the run\'s frame; mask is 1 for land and '
                      '0 for sea, land_fraction each cell\'s dry share. shares are of the Moon\'s surface: today on the '
                      '4 pixel/degree grid, t21 on the model\'s grid, mask_land the binary mask\'s land. '
                      'against_the_run counts the cells whose mask differs from the run\'s own and their share of '
                      'the surface.'),
        shares=dict(today=dict(seas=round(prod['sea_share'], 5), lakes=round(prod['lake_share'], 5)),
                    turned=dict(seas=round(moved['sea_share'], 5), lakes=round(moved['lake_share'], 5)),
                    t21_water=round(summary['water_share'], 5), mask_land=round(summary['mask_land_share'], 5)),
        against_the_run=dict(to_water=int((run_mask & ~mask).sum()), to_land=int((~run_mask & mask).sum()),
                             share_to_water=share(run_mask & ~mask), share_to_land=share(~run_mask & mask),
                             run_mask_land=share(run_mask)),
        lat=[round(float(v), 6) for v in lat], lon=[round(float(v), 6) for v in lon],
        mask=mask.astype(int).tolist(), land_fraction=np.round(fields['land_fraction'], 5).tolist(),
        seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--product', type=Path, help='a routing in lake_rule_balance.py\'s form, in place of the design routing')
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    record = build(args.product, log=lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record) + '\n')
    print(f'wrote {args.out} in {record["seconds"]} s: {record["shares"]}, against the run {record["against_the_run"]}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
