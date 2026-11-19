"""The design climate on the design's own surface: the lakes and the rain set against each other until they agree.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python3 -m research.studies.atlas_round.design_surface_climate   # about 4 minutes a round

The design climate was run on a surface, and the routing of its rain makes another: lakes form where it rains, and
the corrected lake rule holds more of the Moon under them. On the Moon's month-long day the rain falls over land, so a
cell that becomes lake rains less and the land rains more (climate.gcm.surface_change, which measures both on runs
already made). The climate the design routes and the surface the routing makes therefore have to be found together.
No run is made for it (register, 18 November): the climate is estimated from the runs there are.

The surface is held as each T21 cell's water fraction W, in the frame the design climate's run stands in. From W the
climate is estimated (surface_change, the run's settled window pattern-scaled to the design temperature as its base);
on that climate the design's settings are made and its rain routed (lake_rule_balance.py's steps: the factor that
closes the budget, the weirs, the near-side sea's parts, the lakes by the corrected rule), which gives the routing's
own water fractions R(W) (design_mask.py's steps). W then moves part of the way to R(W), and the round is made again
until the two agree: a climate whose lakes are the lakes its own rain fills. Moving the whole way each round, or
carrying the climate by the model's all-or-nothing mask, sets the lakes and the rain chasing each other; the lakes
themselves would not, since they fill over centuries while the rain follows them within a season.

It begins from the surface the run stood on, where the estimate is the run's own climate. The climate's products are
written under the design climate's label each round (climatology_<run>_y<a>-<b>_<label>.npz, its surface wind, and
the climatology on today's ground), so that when the rounds have settled the design chain reads them as they stand:
final_drainage.py, near_side_neck.py and design_routing.py are run after this module. Each round's state is kept on
the research drive (atlas-round/design_surface_climate/), and a run takes up where the last stopped.

The record, results/design_surface_climate.json, gives each round and the settled surface;
results/design_surface.json holds that surface in design_mask.py's form, which the climate's estimate cites.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from climate.gcm import boundary
from climate.gcm import surface_change as SC
from geography import drainage as dr
from research.studies.atlas_round import design_mask as DM
from research.studies.atlas_round import design_works as DW
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import lake_rule as LR
from research.studies.atlas_round import lake_rule_balance as LB
from research.studies.atlas_round import near_side_neck as N
from research.studies.atlas_round import redrawn_drainage as RD

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'design_surface_climate.json'
SURFACE = HERE / 'results' / 'design_surface.json'
WORK = Path('/media/projectspace/terluna-research/atlas-round/design_surface_climate')
STATE = WORK / 'state.json'
SCHEMA = 'terluna.research.atlas-round-design-surface-climate/1'
CODE = ('research/studies/atlas_round/design_surface_climate.py', 'research/studies/atlas_round/lake_rule_balance.py',
        'research/studies/atlas_round/lake_rule.py', 'research/studies/atlas_round/design_mask.py',
        'research/studies/atlas_round/design_works.py', 'research/studies/atlas_round/final_drainage.py',
        'research/studies/atlas_round/near_side_neck.py', 'climate/gcm/surface_change.py')
BASE_LABEL = FD.BASE_LABEL           # the run's settled window at the design temperature (climate.gcm.pattern_scaling)
LABEL = FD.LABEL                     # the design climate's label
TAG = LABEL[len(BASE_LABEL) + 1:]
BASE = (FD.RUN, FD.YEARS[0], FD.YEARS[1], BASE_LABEL)
PAIRS = ((('A28S_301k_design_l15_snow_s930_settled', 35, 49, None), BASE),)   # settled at one temperature on two surfaces
STEP = 0.5                           # the share of the way the surface moves towards the routing's each round
SETTLED_MEAN = 0.004                 # settled when the cells' water fractions differ by less than this on average
SETTLED_CELLS = 12                   # and no more cells than this differ by APART or more
APART = 0.1
ROUNDS = 12


def rel(path) -> str:
    return str(Path(path).relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path)


def run_surface():
    """The surface the design climate's run stood on, as each T21 cell's land fraction, with the grid, the frame's
    record and rotation, read from the base window's product."""
    path = SC.product_path('climatology', *BASE)
    with np.load(path) as z:
        meta = json.loads(str(z['metadata']))
        lat, lon = np.asarray(z['lat'], float), np.asarray(z['lon'], float)
    frame = meta['configuration']['frame']
    surface = boundary.read_t21(SC.PRODUCTS / frame['surface'])
    record = json.loads((ROOT / frame['record']).read_text())
    return dict(land=surface['land'], lat=lat, lon=lon, frame_record=frame['record'],
                rotation=np.asarray(record['frame']['rotation_today_to_turned'], float),
                surface=str((SC.PRODUCTS / frame['surface']).relative_to(ROOT)))


def surface_record(land_fraction, run, note):
    """A surface in design_mask.py's form, as climate.gcm.surface_change reads one."""
    return dict(schema=DM.SCHEMA, run=FD.RUN, years=FD.YEARS, frame=dict(record=run['frame_record']),
                reading_rule=('lat (north first) and lon (east) are the T21 grid\'s in the run\'s frame; land_fraction is '
                              'each cell\'s dry share.'),
                note=note, lat=[round(float(v), 6) for v in run['lat']], lon=[round(float(v), 6) for v in run['lon']],
                land_fraction=np.round(land_fraction, 5).tolist())


def routed(d, st, neck, centre, run, log=print):
    """The design's settings made and its rain routed on the climate final_drainage.py now reads, and the water
    fractions of the surface that routing makes on the T21 grid in the run's frame. Returns the fractions and the
    round's figures."""
    factor, closing = LB.corrected_closing(d, st, log)
    runoff, lake_net, _, _ = FD.load(factor)
    cols = d['h'].shape[1]
    area = (((d['cell_km'] * 1000.0) ** 2 * np.maximum(np.cos(np.radians(d['lat'])), 1e-6))[:, None] * np.ones((1, cols))).ravel()
    net = lake_net.ravel() * dr.MM_DAY * area
    g = LB.govern(d, st, runoff, lake_net, log=log)
    rivers = dict(main=float(g['supply'][LB.MAIN]), south=float(g['supply'][LB.SOUTH]))
    tab = LB.part_tables(d, st, -net, area)
    parts = LB.settle_parts(d, st, tab, neck, centre, runoff, lake_net, net, area, g['fixed'], g['branches'], rivers, log=log)
    r = LR.hollow_routing(parts['s'], parts['terms'], parts['state'])
    water, levels = parts['water'], parts['levels']
    lake = r['lake'] & (water < 0)
    cells = np.flatnonzero(lake)
    prod, sea = DM.surface(d['h'], d['level0'], dict(sea=water, levels=levels, lake_cell=cells,
                                                     lake_level=r['lake_level'][cells], source=dict(sha256='0' * 64)))
    fields, _, _ = DM.mask_in_frame(prod, sea, run['rotation'])
    total = float(area.sum())
    seas = LB.seas_record(list(LB.WATERS), water, levels, d['level0'], area, net, parts['got'], parts['spills'])
    figures = dict(
        open_water_evaporation_factor=factor, surplus_at_factor_1_m3s=closing['at_factor_1']['surplus_m3s'],
        left_over_at_the_factor_m3s=closing['at_its_factor']['surplus_m3s'],
        released_m3s={k: round(v) for k, v in g['gov']['share'].items()},
        neck_flow_m3s=round(parts['z']['neck_flow']), near_side_sea_joined=bool(parts['z']['south'] > neck.crest),
        seas={k: dict(level_above_common_m=v['level_above_common_m'], area_km2=v['area_km2'], receives_m3s=v['receives_m3s'])
              for k, v in seas.items()},
        seas_share=round(float(area[water >= 0].sum()) / total, 5), lakes_share=round(float(area[lake].sum()) / total, 5),
        lakes_not_overflowing_km2=round(float(area[lake & ~r['lake_spills']].sum()) / 1e6))
    return 1.0 - np.asarray(fields['land_fraction'], float), figures


def build(log=print, restart=False):
    started = time.time()
    run = run_surface()
    weight = np.cos(np.radians(run['lat']))[:, None] * np.ones((1, run['lon'].size))
    share = lambda f: float((f * weight).sum() / weight.sum())                   # noqa: E731
    if not (ROOT / FD.TURNED_CLIMATOLOGY).exists():      # the design's modules read the run's frame from the climate's
        SURFACE.write_text(json.dumps(surface_record(run['land'], run, 'the surface the run stood on')) + '\n')   # product
        SC.build(BASE, SURFACE, PAIRS, TAG, log=log)
    d, st = DW.load(log)
    neck_record = json.loads(N.OUT.read_text())
    neck = LB.neck_from_record(st, neck_record)
    centre = neck_record['throat']['centre_today']
    wet = 1.0 - run['land']                              # the surface the climate is estimated on: the run's own first
    rounds = []
    WORK.mkdir(parents=True, exist_ok=True)
    key = dict(code={f: RD.sha256(ROOT / f) for f in CODE}, works=RD.sha256(DW.PRODUCT),
               base=RD.sha256(SC.product_path('climatology', *BASE)))
    if STATE.exists() and not restart:
        kept = json.loads(STATE.read_text())
        if kept['key'] == key:
            wet, rounds = np.asarray(kept['wet']), kept['rounds']
            log(f'taken up after {len(rounds)} rounds')
    settled = bool(rounds) and rounds[-1]['settled']
    estimate = None
    for k in range(len(rounds), ROUNDS):
        if settled:
            break
        SURFACE.write_text(json.dumps(surface_record(1.0 - wet, run, f'the surface of round {k}')) + '\n')
        estimate = SC.build(BASE, SURFACE, PAIRS, TAG, log=log)
        FD.carry(force=True)
        got, figures = routed(d, st, neck, centre, run, log)
        gap = got - wet
        this = dict(round=k, water_share_of_the_surface=dict(estimated_on=round(share(wet), 5), routed=round(share(got), 5)),
                    mean_difference_of_water_fraction=round(share(np.abs(gap)), 5),
                    cells_apart=int((np.abs(gap) >= APART).sum()), largest_difference=round(float(np.abs(gap).max()), 3),
                    land_rain_factor=round(estimate['summary']['land_rain_factor'], 5),
                    surface_rain_mm_day=round(estimate['summary']['surface']['rain_mm_day'][1], 5), **figures)
        this['settled'] = bool(this['mean_difference_of_water_fraction'] < SETTLED_MEAN and this['cells_apart'] <= SETTLED_CELLS)
        rounds.append(this)
        settled = this['settled']
        log(f'round {k}: {json.dumps({kk: vv for kk, vv in this.items() if kk != "seas"})}, {time.time() - started:.0f} s')
        if not settled:
            wet = wet + STEP * gap
        STATE.write_text(json.dumps(dict(key=key, wet=wet.tolist(), rounds=rounds)) + '\n')
    if not settled:
        raise SystemExit(f'the surface and the routing still differ after {len(rounds)} rounds: {rounds[-1]}')
    # the settled surface and its climate, written once more so that the record names the files as they stand
    SURFACE.write_text(json.dumps(surface_record(1.0 - wet, run, 'the settled surface: the routing of the climate estimated '
                                                 'on it makes it again within the record\'s differences')) + '\n')
    estimate = SC.build(BASE, SURFACE, PAIRS, TAG, log=log)
    path = SC.record_path(*BASE[:3], LABEL)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(estimate, indent=1) + '\n')
    carried = FD.carry(force=True)
    last = rounds[-1]
    return dict(
        schema=SCHEMA, producer=dict(files=key['code']),
        inputs={rel(p): RD.sha256(p) for p in (DW.OUT, N.OUT, ROOT / run['frame_record'])},
        ignored_inputs={rel(SC.product_path('climatology', *BASE)): key['base'], str(DW.PRODUCT): key['works'],
                        run['surface']: RD.sha256(ROOT / run['surface'])},
        evidence=('An estimate, not a run: the design climate\'s settled window carried onto the surface its own routing '
                  'makes, by climate.gcm.surface_change, with that surface and the routing set against each other until '
                  'they agree. The routing fills the lakes by the corrected rule and makes the design\'s settings as '
                  'final_drainage.py and design_routing.py make them. The climate\'s response to the surface is '
                  'measured on two settled runs and its placement of the displaced rain is assumed; nothing is '
                  'validated against a run on the settled surface.'),
        reading_rule=('rounds: each round\'s surface as its water share, estimated on and routed; the mean difference of '
                      'the cells\' water fractions between the two and the cells differing by apart or more; the factor '
                      'on the land\'s rain and the surface\'s mean rain in the estimate; and the routing\'s figures (the '
                      'factor that closes its budget, the weirs\' releases, the near-side sea, each sea\'s level above '
                      'the common level, area and supply, the seas\' and the lakes\' shares of the Moon). The surface '
                      'moves step of the way to the routing\'s each round and is settled within settled_mean and '
                      'settled_cells. climate names the estimate\'s products and record.'),
        parameters=dict(step=STEP, settled_mean=SETTLED_MEAN, settled_cells=SETTLED_CELLS, apart=APART),
        label=LABEL, base=dict(run=BASE[0], years=[BASE[1], BASE[2]], label=BASE_LABEL),
        pairs=[[dict(run=w[0], years=[w[1], w[2]], label=w[3]) for w in pair] for pair in PAIRS],
        rounds=rounds, settled=dict(rounds=len(rounds), **{k: last[k] for k in (
            'water_share_of_the_surface', 'mean_difference_of_water_fraction', 'cells_apart', 'land_rain_factor',
            'surface_rain_mm_day', 'open_water_evaporation_factor', 'seas_share', 'lakes_share')}),
        surface=dict(record=rel(SURFACE), sha256=RD.sha256(SURFACE)),
        climate=dict(record=rel(path), products=estimate['products'], on_todays_ground=rel(carried),
                     summary=estimate['summary']),
        seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--restart', action='store_true', help='begin again from the surface the run stood on')
    args = p.parse_args(argv)
    record = build(log=lambda *a: print(time.strftime('%H:%M:%S'), *a, flush=True), restart=args.restart)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(f'wrote {OUT} in {record["seconds"]} s: settled in {record["settled"]["rounds"]} rounds, {record["settled"]}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
