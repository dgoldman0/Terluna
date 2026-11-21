"""The satoyama map on the design surface: satoyama_candidates.py's screen on the design surface's drainage and the
design climate, the run on the design surface at sunlight 0.936 over years 45-59 carried by pattern scaling to 300.31 K
(final_drainage.py, an estimate).

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.satoyama_final   # about a minute

The first map (satoyama_candidates.py) stood on the redrawn floor's design routing, with the near-side ocean, the salt
sea and Humboldtianum held at the common level, and on the present design climate's rain. This one runs the same
screen, through satoyama_candidates.py's own functions and criteria, on:

- the design routing (design_routing.py), which the globe surface draws: every sea at its own level, the near-side sea
  parted at its neck, and the governed lakes and rivers on the design climate's closed budget, so that the floor the
  parted sea leaves dry carries side valleys of its own;
- the design climate's rain and surface temperature carried to today's ground (final_drainage.py);
- the run's settled frame for the turned latitudes (the design surface's).

The screen reads each cell's receiver and drained area from the design routing's product, its discharge in the form
satoyama_candidates.py reads (float32, as the governed routing stores it). The record is results/satoyama_final.json,
in the first map's schema, with the first map's totals beside its own.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from geography import drainage as dr
from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import satoyama_candidates as SC

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'satoyama_final.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/satoyama_final')
CODE = ('research/studies/atlas_round/satoyama_final.py', 'research/studies/atlas_round/satoyama_candidates.py',
        'research/studies/atlas_round/design_routing.py', 'research/studies/atlas_round/final_drainage.py',
        'research/studies/atlas_round/redrawn_drainage.py', 'geography/drainage.py')
SEAS = DR.SEAS                          # globe_surface's bodies 1-7, in order


def sha256(path) -> str:
    return SC.sha256(path)


def load_routing(log=print):
    """The design routing in the form satoyama_candidates.py reads: its arrays with the discharge as float32."""
    routing = DR.load(log)
    routing['discharge_m3s'] = routing['discharge_m3s'].astype(np.float32)
    return routing


def first_map():
    """The first map's totals, from its record."""
    old = json.loads(SC.OUT.read_text())
    c, m = old['candidates'], old['methane']
    return dict(record=str(SC.OUT.relative_to(ROOT)), sha256=sha256(SC.OUT), valleys=old['valleys'],
                land_km2=old['land_km2'], candidates=c['count'], candidates_km2=c['km2'],
                share_of_a_fifth_of_the_land=c['share_of_a_fifth_of_the_land'],
                flooded_possible_km2=c['flooded_possible_km2'],
                methane_at_the_possible_area_tg_yr=m['at_the_possible_area_tg_yr'],
                satoumi=old['satoumi'],
                bands={b['turned_latitude']: dict(candidates_km2=b['candidates_km2'],
                                                  share_of_the_hybrid_fifth=b['share_of_the_hybrid_fifth'])
                       for b in old['bands']})


def build(log=print):
    started = time.time()
    clim_path = FD.carry()                         # the run's climatology on today's ground, current with its sources
    routing = load_routing(log)
    surf = RD.load_surface(log)
    h, level0, lat, cell_km = surf['h'], surf['level0'], surf['lat'], surf['cell_km']
    if abs(level0 - routing['meta']['common_level_m']) > 1e-6:
        raise SystemExit('the routing\'s common level is not the surface\'s')
    rows, cols = h.shape
    clim = dict(np.load(clim_path))
    clim.pop('metadata')
    climate = dict(pr=dr.climate_on_grid(clim, 'pr_mm_day', rows, cols), ts=dr.climate_on_grid(clim, 'ts_k', rows, cols))
    del clim
    frame_path, record = FD.frame_record()
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    c = SC.candidates(h, lat, cell_km, routing, climate, rotation, log)
    del climate
    area = c['cell_area']
    sat = SC.satoumi(h, level0, routing, area, lat, cols, rotation)
    out = SC.record(c, routing, level0, area, sat, rotation)
    BULK.mkdir(parents=True, exist_ok=True)
    SC.picture(c, routing, h, BULK / 'candidates_map.png')
    pos = np.full(rows * cols, -1, np.int32)
    members = c['unit_cell'] >= 0
    pos[members] = np.searchsorted(c['ids'], c['unit_cell'][members]).astype(np.int32)
    np.savez_compressed(BULK / 'satoyama_units.npz', unit=pos, candidate=c['ok'], outlet_cell=c['ids'].astype(np.int32),
                        flooded_possible_km2=c['possible'].astype(np.float32), woodland_km2=c['wood'].astype(np.float32),
                        valley_floor_km2=c['valley'].astype(np.float32), area_km2=c['area'].astype(np.float32))
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)       # noqa: E731
    return dict(
        schema=SC.SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (DR.OUT, FD.OUT, frame_path)},
        ignored_inputs={str(DR.PRODUCT): sha256(DR.PRODUCT), rel(clim_path): sha256(clim_path)},
        routing=dict(record=rel(DR.OUT), path=str(DR.PRODUCT), seas=routing['meta']['seas'],
                     levels_above_the_common_level_m={k: round(v - level0, 2)
                                                      for k, v in routing['meta']['levels_m'].items()},
                     climate=routing['meta']['climate']),
        evidence=('A screen on products: side valleys of the design surface\'s drainage (every sea at its own level, the '
                  'near-side sea parted, the governed lakes and rivers on the design climate\'s closed budget) at 1.9 '
                  'km cells, tested against criteria set from Earth\'s satoyama, tank cascades and terraces '
                  f'(satoyama.md), with the rain and surface temperature at T21 of {FD.CLIMATE}, carried to '
                  'today\'s ground, and the IPCC 2019 Refinement\'s methane factors for rice. Slopes at 1.9 km are '
                  'gentler than the ground within a cell, so the valley floors and the woodland are first estimates. '
                  'Not validated against any observation.'),
        reading_rule=('As satoyama_candidates.json: areas in km2, rain in mm a year, temperatures the climate model\'s '
                      'surface temperature in C, methane in Tg a year; places on today\'s Moon and the turned Moon in '
                      'the run\'s settled frame. first_map gives the first map\'s totals beside these.'),
        criteria=dict(valley_km2=list(SC.AREA_KM2), gentle_slope=SC.GENTLE, steep_deg=SC.STEEP_DEG,
                      valley_floor_share_min=SC.VALLEY_SHARE_MIN, woodland_share_min=SC.WOODLAND_SHARE_MIN,
                      rain_mm=list(SC.RAIN_MM), season_days=SC.SEASON_DAYS, need_mm=list(SC.NEED_MM),
                      irrigated_share_max=SC.IRRIGATED_SHARE_MAX),
        **out, first_map=first_map(), seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    product = build(log=lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(product, indent=1, ensure_ascii=False) + '\n')
    c, f = product['candidates'], product['first_map']
    print(f"{c['count']:,} candidates over {c['km2']:,} km2 (first map {f['candidates']:,}, {f['candidates_km2']:,} km2); "
          f"flooded possible {c['flooded_possible_km2']:,} km2 (first map {f['flooded_possible_km2']:,})")
    print('bands', json.dumps([(b['turned_latitude'], b['candidates_km2'], b['share_of_the_hybrid_fifth']) for b in product['bands']]))
    print('methane', product['methane']['at_the_possible_area_tg_yr'], 'satoumi', json.dumps(product['satoumi']))
    print(f'wrote {args.out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
