"""The atlas round's drainage on the design climate: the design climate's settled years carried onto today's ground
and routed over the redrawn floor with the works chosen on 16 November.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.final_drainage   # about 20 minutes

The design climate (register, 18 November) is the climate run on the design surface (A28D at sunlight 0.936, years 45-59,
301.00 K at the ground), carried by pattern scaling to the 300.31 K at which the run at sunlight 0.930 on the final
surface settled (A28S, years 35-49): climate.gcm.pattern_scaling moves its climatology and surface wind by the land's and
the water's mean responses per kelvin, measured on two pairs of runs, into climatology_<run>_y45-59_scaled300.31K.npz,
an estimate. A second estimate carries that onto the surface the design's own routing makes, where the corrected lake
rule holds more of the Moon under lakes than the run's surface did: design_surface_climate.py sets the surface and
the routing against each other, and climate.gcm.surface_change writes the climate under the label
scaled300.31K_resurfaced, which this module reads. Until 18 November this module read the run at 0.930 on the final
surface itself. The run keeps its fields
in the settled turned frame of the surface it ran on, the record its climatology names (design_turned_surface.py's for
the design surface), so the climatology is carried onto today's latitude and longitude as climate.gcm.turned_climatology
carries one, with that record's rotation, and written beside it as ..._today.npz.

The routing fills its lakes by the corrected rule (lake_rule.py: each hollow as one water) and sets the design's
settings as lake_rule_balance.py worked them out (register, 18 November): the factor on the open water's evaporation
that closes the routing's budget with every sea at the common level, one step from 1; then, on each budget, the design
state with the three upper seas full and spilling and the weirs set until the salt sea and Humboldtianum receive their
need, a cut taking no more than its lake can spare. The floor's fixed waters and the works' cells are design_works.py's.
The near-side sea is routed at the common level and the level its balance would bring it to is reported;
near_side_neck.py and design_routing.py part it. The record sets the numbers beside those of the earlier design
climate on the turned Moon at 8 degrees (governed_drainage.json), which stays as it is. Each budget's governed routing
goes to the research drive in the governed routing's form, its receivers joined through the lakes (hollow_courses.py).
Until 18 November the module routed with governed_drainage.py's and redrawn_drainage.py's functions on the earlier
rule, which balanced each catchment by itself.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from climate.gcm import turned_climatology as TC
from geography import drainage as dr
from research.studies.atlas_round import governed_drainage as GD
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'final_drainage.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/final_drainage')
SCHEMA = 'terluna.research.atlas-round-final-drainage/1'
PRODUCT_SCHEMA = GD.PRODUCT_SCHEMA
CODE = ('research/studies/atlas_round/final_drainage.py', 'research/studies/atlas_round/lake_rule_balance.py',
        'research/studies/atlas_round/lake_rule.py', 'research/studies/atlas_round/hollow_courses.py',
        'research/studies/atlas_round/design_works.py', 'research/studies/drainage_network/network.py',
        'geography/drainage.py', 'climate/gcm/turned_climatology.py', 'climate/gcm/turned.py')
SEAS = ('moscoviense', 'orientale', 'smythii_marginis', 'near_side', 'south_pole_aitken', 'humboldtianum')
SEAS_FIXED = ('moscoviense', 'orientale', 'smythii_marginis', 'south_pole_aitken', 'humboldtianum')   # design_works' fixed waters
RUN = 'A28D_301k_design_l15_snow_s936_settled'
YEARS = [45, 59]
BASE_LABEL = 'scaled300.31K'                # climate.gcm.pattern_scaling's estimate at the run at 0.930's 300.31 K
LABEL = f'{BASE_LABEL}_resurfaced'          # and climate.gcm.surface_change's on the surface the design's routing makes
BASE_WINDOW = f'{RUN}_y{YEARS[0]}-{YEARS[1]}_{BASE_LABEL}'
WINDOW = f'{RUN}_y{YEARS[0]}-{YEARS[1]}_{LABEL}'
BASE_TURNED_CLIMATOLOGY = f'climate/gcm/products/climatology_{BASE_WINDOW}.npz'
TURNED_CLIMATOLOGY = f'climate/gcm/products/climatology_{WINDOW}.npz'
CLIMATOLOGY = f'climate/gcm/products/climatology_{WINDOW}_today.npz'
SURFACE_WIND = f'climate/gcm/products/surface_wind_{WINDOW}.npz'
ESTIMATE = f'climate/results/gcm/pattern_scaling_{BASE_WINDOW}.json'      # the step to the design temperature
SURFACE_ESTIMATE = f'climate/results/gcm/surface_change_{WINDOW}.json'    # the step onto the design's own surface
CLIMATE = ('the design climate: the climate run on the design surface (sunlight 0.936, years 45-59, 301.00 K at the '
           'ground) carried by pattern scaling to the 300.31 K of the run at 0.930 on the final surface and by the '
           'contrast between land and water onto the surface its own routing makes, an estimate')
BUDGETS = ('routing_budget', 'closed_budget')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frame_record():
    """The settled turned frame the run stands in: the record its climatology names, unchanged since the run read it."""
    with np.load(ROOT / TURNED_CLIMATOLOGY) as z:
        frame = json.loads(str(z['metadata']))['configuration']['frame']
    path = ROOT / frame['record']
    if sha256(path)[:len(frame['sha256'])] != frame['sha256']:
        raise SystemExit(f'{frame["record"]} has changed since {RUN} ran in it')
    return path, json.loads(path.read_text())


def routing_path(budget):
    return BULK / f'final_routing_{budget}.npz'


def surface_ratios():
    """What the step onto the design's own surface does to each T21 cell's rain and evaporation: the design climate's
    over the run's at the design temperature, in the run's frame (1 where the base has none). The maps that read the
    run's years one by one move each year by these after the step to the design temperature."""
    with np.load(ROOT / BASE_TURNED_CLIMATOLOGY) as b, np.load(ROOT / TURNED_CLIMATOLOGY) as e:
        out = {}
        for key, name in (('pr', 'pr_mm_day'), ('evap', 'evap_mm_day')):
            base, new = np.asarray(b[name], float), np.asarray(e[name], float)
            out[key] = np.where(np.abs(base) > 1e-9, new / np.where(np.abs(base) > 1e-9, base, 1.0), 1.0)
        return out


# ---------------------------------------------------------------- the climate on today's ground

def carry(force=False):
    """Write the run's climatology on today's latitude and longitude (turned_climatology's method and form) unless it is
    there already, carried from the climatology it names; returns its path."""
    out = ROOT / CLIMATOLOGY
    src = ROOT / TURNED_CLIMATOLOGY
    path, record = frame_record()
    if out.exists() and not force:
        meta = json.loads(str(np.load(out)['metadata']))
        carried = meta.get('carried', {})
        if carried.get('climatology', {}).get('sha256') == sha256(src) and \
                carried.get('frame_record', {}).get('sha256') == sha256(path):
            return out
    z = np.load(src)
    meta = json.loads(str(z['metadata']))
    if meta.get('schema') != TC.SCHEMA or meta.get('run') != RUN or list(meta.get('years')) != YEARS or \
            meta.get('estimate', {}).get('label') != LABEL:
        raise SystemExit(f'{src} is not the climatology of {RUN} years {YEARS} estimated as {LABEL}')
    clim = {k: z[k] for k in z.files if k != 'metadata'}
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    fields = TC.to_today(clim, rotation)
    meta.update(frame='today', carried=dict(
        method=('each cell of today\'s grid takes the turned field at its turned place, bilinear on the turned grid; '
                'the Sun composites are left out (climate.gcm.turned_climatology.to_today)'),
        climatology=dict(file=src.name, sha256=sha256(src)),
        frame_record=dict(path=str(path.relative_to(ROOT)), sha256=sha256(path), state=record['frame']['state']),
        producer={'research/studies/atlas_round/final_drainage.py': sha256(Path(__file__))[:16],
                  'climate/gcm/turned_climatology.py': sha256(Path(TC.__file__))[:16]}))
    np.savez_compressed(out, metadata=json.dumps(meta), lat=clim['lat'], lon=clim['lon'], **fields)
    return out


def load(evaporation=1.0, rain=1.0):
    """The run's rain, evaporation, runoff and open-water balance (mm/day) on the 16 pixel/degree grid, as the network
    study's loader takes them from a climatology (network.load), with the open water's evaporation scaled."""
    z = np.load(ROOT / CLIMATOLOGY)
    meta = json.loads(str(z['metadata']))
    if meta.get('run') != RUN or list(meta.get('years')) != YEARS or meta.get('frame') != 'today' or \
            meta.get('estimate', {}).get('label') != LABEL:
        raise SystemExit(f'{CLIMATOLOGY} is not {RUN} years {YEARS} estimated as {LABEL} on today\'s ground')
    clim = {k: z[k] for k in z.files if k != 'metadata'}
    rows, cols = 180 * dr.PPD, 360 * dr.PPD
    lat = 90.0 - (np.arange(rows) + 0.5) / dr.PPD
    precip = (dr.climate_on_grid(clim, 'pr_mm_day', rows, cols) * np.float32(rain)).astype(np.float32)
    evap = dr.climate_on_grid(clim, 'evap_mm_day', rows, cols).astype(np.float32)
    runoff = np.maximum(precip - evap, 0.0).astype(np.float32)
    lake_net = (precip - np.float32(evaporation) * dr.lake_evaporation(clim, lat)[:, None].astype(np.float32)).astype(np.float32)
    return runoff, lake_net, precip, evap


# ---------------------------------------------------------------- the routing

def near_side_would_settle(d, st, supply, loss_field, area):
    """The level the near-side sea's balance would bring it to on a supply, as governed_drainage.near_side_level gives
    it: where its loss to the air meets the supply, up to its rim, the lowest pass between it and another sea."""
    from research.studies.atlas_round import lake_rule_balance as LB
    from research.studies.atlas_round import near_side_neck as N
    level0, hf = d['level0'], d['h'].ravel()
    joined = st['joined_main']                           # the level each cell joins the sea's water at, from Imbrium
    rim = min(float(joined[int(np.flatnonzero(d['water'] == i)[np.argmin(hf[d['water'] == i])])])
              for i in range(len(SEAS_FIXED)))
    table = LB.by_level(NW.Table(joined.astype(np.float64), area, loss_field, hf, rim + 1.0))
    held = N.level_for(table, supply, rim - NW.BELOW_RIM_M)
    at = table.at(rim - NW.BELOW_RIM_M if held is None else held)
    return dict(reaches_its_rim=held is None, level_m=round((rim if held is None else held) - level0, 1),
                area_km2=round(at[0] / 1e6), area_at_the_common_level_km2=round(table.at(level0)[0] / 1e6),
                rim_above_common_m=round(rim - level0, 1))


def governed(d, st, evaporation, start=None, log=print):
    """The weirs set on a budget's climate with the near-side sea at the common level (lake_rule_balance.govern), the
    routing with its courses joined through the lakes (hollow_courses.connect), and the record's figures for it."""
    from research.studies.atlas_round import hollow_courses as HC
    from research.studies.atlas_round import lake_rule_balance as LB
    runoff, lake_net, _, _ = load(evaporation)
    g = LB.govern(d, st, runoff, lake_net, start=start, log=log, routing=True)
    gov, supply, water, levels, r = g['gov'], g['supply'], g['water'], g['levels'], g['r']
    c = HC.connect(g['s'], g['terms'], g['state'], r, d['h'].shape)
    area = r['area']
    loss_field = -r['net']
    names = list(LB.WATERS)
    own, loss = gov['own'], gov['loss']
    seas = {}
    for k in SEAS:                                       # the near-side sea whole, its two parts at one level
        idx = [LB.MAIN, LB.SOUTH] if k == 'near_side' else [names.index(k)]
        cells = np.isin(water, idx)
        got, lost, mine = (float(sum(v[i] for i in idx)) for v in (supply, loss, own))
        seas[k] = dict(level_m=round(levels[idx[0]] - d['level0'], 2), area_km2=round(float(area[cells].sum()) / 1e6),
                       loss_m3s=round(lost), own_rivers_m3s=round(mine), supply_m3s=round(got), balance_m3s=round(got - lost))
    near = float(supply[LB.MAIN] + supply[LB.SOUTH])
    seas['near_side']['would_settle'] = near_side_would_settle(d, st, near, loss_field, area)
    seas['near_side']['parts_rivers_m3s'] = dict(main=round(float(supply[LB.MAIN])), south=round(float(supply[LB.SOUTH])))
    lake = r['lake'] & (water < 0)
    total = float(area.sum())
    record = dict(
        evaporation=evaporation, moscoviense_surplus_m3s=round(gov['surplus']),
        released_m3s={k: round(v) for k, v in gov['share'].items()},
        arrived_m3s={k: round(float(supply[names.index(k)] - own[names.index(k)])) for k in (GD.SALT_SEA, GD.HELD_BRANCH, GD.REST_BRANCH)},
        full_seas_pass_on_m3s={k: round(v) for k, v in gov['over'].items()},
        held_seas_left_short_m3s={k: round(v) for k, v in gov['unmet'].items()},
        branch_flows_m3s={k: [round(q, 1) for q in v] for k, v in gov['flows'].items()},
        seas=seas, lakes=dict(km2=round(float(area[lake].sum()) / 1e6), share_of_moon=round(float(area[lake].sum()) / total, 5),
                              passing_nothing_on_km2=round(float(area[lake & ~r['lake_spills']].sum()) / 1e6)),
        weirs_set_in_routings=len(gov['history']), weir_history=gov['history'])
    sea = np.full(water.size, -1, np.int8)
    for i, k in enumerate(SEAS):
        sea[np.isin(water, [LB.MAIN, LB.SOUTH] if k == 'near_side' else [names.index(k)])] = i
    routing = dict(receiver=c['receiver'], discharge=c['discharge'], drained_m2=c['drained_m2'], sea=sea, lake=lake,
                   lake_level=r['lake_level'], lake_spills=r['lake_spills'],
                   levels={k: float(levels[LB.MAIN if k == 'near_side' else names.index(k)]) for k in SEAS},
                   shares=gov['share'])
    return record, routing


def product(routing, level0, shape, budget_label, evaporation, path):
    """A budget's governed routing in the governed routing's form, with the climate it was routed on."""
    cells = np.flatnonzero(routing['lake']).astype(np.int32)
    meta = dict(schema=PRODUCT_SCHEMA, common_level_m=level0, seas=list(SEAS), levels_m=routing['levels'],
                grid=dict(pixels_per_degree=dr.PPD, rows=shape[0], cols=shape[1]),
                shares_m3s={k: round(v, 1) for k, v in routing['shares'].items()},
                climate=dict(run=RUN, years=YEARS, climatology=CLIMATOLOGY, budget=budget_label,
                             open_water_evaporation_factor=evaporation),
                reading_rule=('As design_routing.npz: cells index rows from 90 N and columns from 0 E, 16 per degree, on '
                              'today\'s Moon; receiver, discharge_m3s and drained_km2 per cell, sea the index into seas, '
                              'lake_cell with lake_level_m and lake_spills. The lakes are filled by the corrected rule '
                              '(lake_rule.py) and the receivers joined through them (hollow_courses.py). The weirs and '
                              'cuts of 16 November divide Moscoviense\'s water on the climate named, and the branches\' '
                              'water runs in the rivers it joins.'))
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, metadata=json.dumps(meta), receiver=routing['receiver'].astype(np.int32),
                        discharge_m3s=routing['discharge'].astype(np.float32),
                        drained_km2=(routing['drained_m2'] / 1e6).astype(np.float32), sea=routing['sea'],
                        lake_cell=cells, lake_level_m=routing['lake_level'][cells].astype(np.float32),
                        lake_spills=routing['lake_spills'][cells])
    return dict(path=str(path), sha256=sha256(path))


def against_the_present_climate(cases):
    """The governed numbers that move from the present design climate's routing to this one."""
    old = json.loads(GD.OUT.read_text())['cases']
    out = {}
    for label in BUDGETS:
        new, was = cases[label], old[label]
        out[label] = dict(
            moscoviense_surplus_m3s=[was['moscoviense_surplus_m3s'], new['moscoviense_surplus_m3s']],
            released_m3s={k: [was['released_m3s'][k], new['released_m3s'][k]] for k in new['released_m3s']},
            full_seas_pass_on_m3s={k: [was['full_seas_pass_on_m3s'][k], new['full_seas_pass_on_m3s'][k]]
                                   for k in new['full_seas_pass_on_m3s']},
            seas={k: dict(loss_m3s=[was['seas'][k]['loss_m3s'], new['seas'][k]['loss_m3s']],
                          own_rivers_m3s=[was['seas'][k]['own_rivers_m3s'], new['seas'][k]['own_rivers_m3s']],
                          balance_m3s=[was['seas'][k]['balance_m3s'], new['seas'][k]['balance_m3s']])
                  for k in new['seas']},
            near_side_would_settle_m=[was['seas']['near_side']['would_settle']['level_m'],
                                      new['seas']['near_side']['would_settle']['level_m']])
    return out


def build(log=print):
    from research.studies.atlas_round import design_works as DW
    from research.studies.atlas_round import lake_rule_balance as LB
    started = time.time()
    clim_path = carry()
    log(f'climatology on today\'s ground: {clim_path.name}, {time.time() - started:.0f} s')
    d, st = DW.load(log)
    frame_path, _ = frame_record()
    factor, closing = LB.corrected_closing(d, st, log)
    log(f'the routing\'s budget closes at {factor}, {time.time() - started:.0f} s')
    cases, products = {}, {}
    for label, evaporation in zip(BUDGETS, (1.0, factor)):
        record, routing = governed(d, st, evaporation, log=log)
        products[label] = product(routing, d['level0'], d['h'].shape, label, evaporation, routing_path(label))
        cases[label] = record
        log(f'{label}: near side {record["seas"]["near_side"]}, {time.time() - started:.0f} s')
        del routing
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (DW.OUT, frame_path, GD.OUT, ROOT / ESTIMATE)},
        ignored_inputs={rel(clim_path): sha256(clim_path), TURNED_CLIMATOLOGY: sha256(ROOT / TURNED_CLIMATOLOGY),
                        str(DW.PRODUCT): sha256(DW.PRODUCT)},
        climatology=dict(run=RUN, years=YEARS, estimate=LABEL, frame=rel(frame_path)),
        evidence=(f'A calculation on products: the routing of {CLIMATE}, carried from its settled turned frame onto '
                  'today\'s ground, over the redrawn floor with the works chosen on 16 November, the lakes filled by '
                  'the corrected rule (each hollow as one water) and each weir set by settling again until the sea it '
                  'serves receives its need. Not validated against any observation.'),
        reading_rule=('Levels in m above the redrawn common level, flows in m3/s, areas in km2; cases on the routing\'s '
                      'own budget and on the one its open-water evaporation factor closes. budget gives the routing '
                      'with every sea at the common level at factor 1 and at the factor that one step gives, where '
                      'surplus_m3s is what the budget still leaves over. In each case the near-side sea stands at the '
                      'common level, released is what each weir or branch lets go and arrived what reaches its sea '
                      'beyond its own rivers; would_settle gives the level the near-side sea\'s balance would bring it '
                      'to, up to its rim. against_the_present_climate pairs each number as [the earlier design climate '
                      'at 8 degrees on the earlier lake rule, this routing]. routing names each budget\'s governed '
                      'routing on the research drive.'),
        budget=closing['at_factor_1'], budget_at_its_factor=closing['at_its_factor'],
        cases=cases, against_the_present_climate=against_the_present_climate(cases), routing=products,
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
