"""The drylands on the atlas round's final surface, at the two sunlights the night of 16 November ran, against the
corrected design climate on the atlas's surface, and how much of them the rivers' green land would take.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.drylands   # several minutes

The land is read as the moisture map reads it (biosphere/ecology/moisture_map.py): the run's rain over Priestley and
Taylor's demand from the means over the ten tenths of the lunar day, UNEP's classes of that aridity index, each land
cell of the run's own T21 grid weighted by its area. The windows:

- the final surface at sunlight 0.930, years 35-49, settled near 300.3 K;
- the final surface at sunlight 0.938, years 1-8, still warming from 301.1 to 301.7 K toward about 302 K, so its land
  is a first look, not a settled climate;
- the corrected design run on the atlas's surface (sunlight 0.938, 301.0 K), years 10-19, the design years' window,
  in its own elastic frame, as the reference.

The rivers' green land (river_losses.py's river_green.npz on the research drive): each river cell's water surface and a
belt of plants drinking its bed's seepage, at the narrow and the wide ends of their ranges, and the floodplains within
5 m of their river that it can keep wet, carried into the settled frame and onto the final surface's runs' T21 grid.
The green land counted is what falls on the run's land cells of each class; the cells keep their class, since the
corridors are narrower than a cell.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import moisture_map as bm
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import river_losses as RL
from research.studies.atlas_round import turned_surface as TS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'drylands.json'
SCHEMA = 'terluna.research.atlas-round-drylands/1'
CODE = ('research/studies/atlas_round/drylands.py', 'research/studies/atlas_round/river_losses.py',
        'biosphere/ecology/moisture_map.py')
RUNS = Path('/media/projectspace/terluna-research/gcm-runs')
CACHE = Path('/media/projectspace/terluna-research/atlas-round/drylands/cache')
REFERENCE_CACHE = Path('/media/projectspace/terluna-research/session-tools/precession/peek/cache')   # read only
WINDOWS = (('final_surface_0.930', 'A28S_301k_design_l15_snow_s930_settled', (35, 49), True),
           ('final_surface_0.938_first_years', 'A28S_301k_design_l15_snow_s938_settled', (1, 8), True),
           ('corrected_design_atlas_surface', 'A28R_301k_design_l15_snow_s938_turned_refit_precession', (10, 19), False))
CLASSES = [name for name, _, _ in bm.UNEP]


def shares(c):
    """The land's area and its share in each UNEP class, with the run's land rain and demand (mm/day)."""
    area = bm.cell_areas_km2(c['lat'], c['lon'])
    land = c['land']
    cls = bm.unep_class(c['ai'])
    total = float(area[land].sum())
    out = {name: round(float(area[land & (cls == k)].sum() / total), 4) for k, name in enumerate(CLASSES)}
    return dict(land_km2=round(total, 0), classes=out,
                drylands=round(sum(out[n] for n in CLASSES[:4]), 4),
                arid_and_semi_arid=round(out['arid'] + out['semi-arid'], 4),
                land_rain_mm_day=round(float((c['rain'] * area)[land].sum() / total), 3),
                land_demand_mm_day=round(float((c['pet'] * area)[land].sum() / total), 3)), cls, area


def ground(run, years):
    prog = json.loads((RUNS / run / 'progress.json').read_text())
    t = [y['surface_k'] for y in prog['years'] if years[0] <= y['year'] <= years[1]]
    return round(float(np.mean(t)), 2)


def green_cells(rotation):
    """The rivers' green land on the T21 grid of the settled frame (km2 per cell): the narrow and the wide readings of
    the channels and the belt, and the floodplains kept wet."""
    with np.load(RL.GREEN, allow_pickle=False) as d:
        meta = json.loads(str(d['metadata']))
        cols, ppd = meta['grid']['cols'], meta['grid']['pixels_per_degree']
        rc, rl, rq = d['river_cell'], d['river_length_km'].astype(float), d['river_flow_m3s'].astype(float)
        pc, pa, pk = d['plain_cell'], d['plain_area_km2'].astype(float), d['plain_kept_share'].astype(float)
    lat_of = lambda cells: 90.0 - (cells // cols + 0.5) / ppd
    lon_of = lambda cells: (cells % cols + 0.5) / ppd
    out = {}
    rr, rcol = RL.t21_cells(lat_of(rc), lon_of(rc), rotation)
    prow, pcol = RL.t21_cells(lat_of(pc), lon_of(pc), rotation)
    plains = np.zeros((32, 64))
    np.add.at(plains, (prow, pcol), pa * pk)
    for label, coef, belt in (('narrow', RL.WIDTH[0], RL.RIPARIAN_KM[0]), ('wide', RL.WIDTH[1], RL.RIPARIAN_KM[1])):
        corridor = rl * (RL.channel_width_m(rq, coef) / 1000.0 + belt)
        grid = np.zeros((32, 64))
        np.add.at(grid, (rr, rcol), corridor)
        out[label] = grid + plains
    out['floodplains'] = plains
    return out


def build(log=print):
    started = time.time()
    frame = json.loads(TS.record_path(RD.sha256(TS.ST.OUT)).read_text())
    rotation = np.asarray(frame['frame']['rotation_today_to_turned'])
    green = green_cells(rotation)
    windows = {}
    bm.RUNS = RUNS
    for label, run, years, final in WINDOWS:
        bm.CACHE = REFERENCE_CACHE if not final else CACHE
        c = bm.read_run(run, years, humidity=False)
        summary, cls, area = shares(c)
        summary.update(run=run, years=list(years), ground_k=ground(run, years))
        if final:
            land = c['land']
            total = float(area[land].sum())
            summary['rivers_green_land'] = {}
            for reading in ('narrow', 'wide', 'floodplains'):
                g = green[reading]
                by = {name: round(float(g[land & (cls == k)].sum()), 0) for k, name in enumerate(CLASSES)}
                dry = sum(by[n] for n in CLASSES[:3])
                summary['rivers_green_land'][reading] = dict(
                    km2_by_class=by, on_hyper_arid_arid_and_semi_arid_km2=round(dry, 0),
                    share_of_the_land=round(dry / total, 5),
                    arid_and_semi_arid_left=round(summary['arid_and_semi_arid']
                                                  - (by['arid'] + by['semi-arid']) / total, 4))
        windows[label] = summary
        log(f'{label}: {summary["ground_k"]} K, drylands {summary["drylands"]}, arid and semi-arid '
            f'{summary["arid_and_semi_arid"]}, classes {summary["classes"]}')
    return dict(
        schema=SCHEMA, producer=dict(files={f: RD.sha256(ROOT / f) for f in CODE}),
        inputs={str(p.relative_to(ROOT)): RD.sha256(p) for p in (RL.OUT, TS.record_path(RD.sha256(TS.ST.OUT)))},
        ignored_inputs={str(RL.GREEN): RD.sha256(RL.GREEN)},
        windows=windows,
        evidence=('GCM output read as the moisture map reads it (rain over Priestley-Taylor demand, UNEP classes, land '
                  'by area), and the rivers\' green land from river_losses.py carried onto the runs\' grid. The 0.938 '
                  'window is still warming; the reference stands on the atlas\'s surface in the elastic frame. Not '
                  'validated against observation.'),
        reading_rule=('classes are shares of each window\'s own land; drylands are the four dry classes and arid_and_'
                      'semi_arid their middle two; rivers_green_land gives the channels, belts and floodplains (narrow '
                      'and wide readings, and the floodplains alone) falling on land of each class, and the arid and '
                      'semi-arid share left once that land is taken out.'),
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
