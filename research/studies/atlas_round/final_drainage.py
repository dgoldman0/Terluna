"""The atlas round's drainage on the climate run on the final surface: the 0.930 run's settled years carried onto
today's ground and routed over the redrawn floor with the works chosen on 16 November.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.final_drainage   # about 20 minutes

The run on the final surface (A28S at sunlight 0.930, the third step's restart) settled at 300.3 K at the ground; its
years 35-49 are averaged by climate.gcm.climatology into climatology_<run>_y35-49.npz. The run keeps its fields in the
settled turned frame of the surface it ran on (turned_surface.py's record), so the climatology is carried onto today's
latitude and longitude as climate.gcm.turned_climatology carries one, with that record's rotation, and written beside
it as ..._y35-49_today.npz.

The routing is governed_drainage.py's, unchanged, on this climate: the network study's loader and budget closure as
network.py has them (the open water's evaporation raised until the routing's budget closes, with every sea at the
common level), the design state with the three upper seas full and spilling, the cuts' routes found again on the new
lakes, and the weirs set again on each budget until the salt sea and Humboldtianum receive their need. The near-side
ocean is routed at the common level and the level its balance would bring it to is reported, as before. The record sets
the numbers beside those of the present design climate on the turned Moon at 8 degrees (governed_drainage.json), which
stays as it is. Each budget's routing goes to the research drive in the governed routing's form.
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
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import turned_surface as TS
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'final_drainage.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/final_drainage')
SCHEMA = 'terluna.research.atlas-round-final-drainage/1'
PRODUCT_SCHEMA = GD.PRODUCT_SCHEMA
CODE = ('research/studies/atlas_round/final_drainage.py', 'research/studies/atlas_round/governed_drainage.py',
        'research/studies/atlas_round/redrawn_drainage.py', 'research/studies/drainage_network/network.py',
        'research/studies/drainage_network/routing.py', 'geography/drainage.py', 'climate/gcm/turned_climatology.py',
        'climate/gcm/turned.py')
RUN = 'A28S_301k_design_l15_snow_s930_settled'
YEARS = [35, 49]
TURNED_CLIMATOLOGY = f'climate/gcm/products/climatology_{RUN}_y{YEARS[0]}-{YEARS[1]}.npz'
CLIMATOLOGY = f'climate/gcm/products/climatology_{RUN}_y{YEARS[0]}-{YEARS[1]}_today.npz'
BUDGETS = ('routing_budget', 'closed_budget')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frame_record():
    """The settled turned frame the run stands in: the final surface's record, which the run names."""
    path = TS.record_path(RD.sha256(TS.ST.OUT))
    return path, json.loads(path.read_text())


def routing_path(budget):
    return BULK / f'final_routing_{budget}.npz'


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
    if meta.get('schema') != TC.SCHEMA or meta.get('run') != RUN or list(meta.get('years')) != YEARS:
        raise SystemExit(f'{src} is not the climatology of {RUN} years {YEARS}')
    if meta['configuration']['frame']['record'] != str(path.relative_to(ROOT)):
        raise SystemExit(f'{RUN} ran in {meta["configuration"]["frame"]["record"]}, not {path.relative_to(ROOT)}')
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
    if meta.get('run') != RUN or list(meta.get('years')) != YEARS or meta.get('frame') != 'today':
        raise SystemExit(f'{CLIMATOLOGY} is not {RUN} years {YEARS} on today\'s ground')
    clim = {k: z[k] for k in z.files if k != 'metadata'}
    rows, cols = 180 * dr.PPD, 360 * dr.PPD
    lat = 90.0 - (np.arange(rows) + 0.5) / dr.PPD
    precip = (dr.climate_on_grid(clim, 'pr_mm_day', rows, cols) * np.float32(rain)).astype(np.float32)
    evap = dr.climate_on_grid(clim, 'evap_mm_day', rows, cols).astype(np.float32)
    runoff = np.maximum(precip - evap, 0.0).astype(np.float32)
    lake_net = (precip - np.float32(evaporation) * dr.lake_evaporation(clim, lat)[:, None].astype(np.float32)).astype(np.float32)
    return runoff, lake_net, precip, evap


def budget(geo, runoff, lake_net, precip, evap, area):
    """The routing's budget with every sea at the common level and the open water's evaporation as loaded, and the
    factor on that evaporation that closes it (network.network's budget)."""
    level0, h, lat, cell_km = geo['level0'], geo['h'], geo['lat'], geo['cell_km']
    at_common = [(geo['joined'][k] < level0, level0) for k in geo['seas']]
    r, rivers = NW.route_at(h, lat, cell_km, runoff, lake_net, at_common)
    sea_cells = np.logical_or.reduce([m for m, _ in at_common])
    lake = r['lake'] & ~sea_cells
    land = ~sea_cells & ~r['lake']
    mm = dr.MM_DAY
    loss = -lake_net.ravel().astype(np.float64) * mm * area
    out = dict(
        land_rain_less_evaporation_m3s=round(float(((precip - evap).ravel() * mm * area)[land].sum())),
        routing_runoff_m3s=round(float((runoff.ravel() * mm * area)[land].sum())),
        lakes_rain_less_evaporation_m3s=round(float((lake_net.ravel() * mm * area)[lake].sum())),
        seas_loss_m3s=round(float(loss[sea_cells].sum())),
        delivered_to_the_seas_m3s=round(sum(rivers)),
        lakes_km2=round(float(area[lake].sum()) / 1e6),
        closed_lakes_km2=round(float(area[lake & r['closed'][r['label']]].sum()) / 1e6))
    out['surplus_m3s'] = out['delivered_to_the_seas_m3s'] - out['seas_loss_m3s']
    open_water = float((NW.evap_lake_of(lake_net, precip).ravel() * mm * area)[sea_cells | lake].sum())
    out['open_water_evaporation_m3s'] = round(open_water)
    out['open_water_evaporation_factor_closing_it'] = round(1.0 + out['surplus_m3s'] / open_water, 3)
    return out


# ---------------------------------------------------------------- the routing

def product(gov, state, geo, budget_label, evaporation, path):
    """A budget's governed routing in the governed routing's form, with the climate it was routed on."""
    r = gov['r']
    lake = r['lake'] & ~r['outlet']
    cells = np.flatnonzero(lake).astype(np.int32)
    meta = dict(schema=PRODUCT_SCHEMA, common_level_m=geo['level0'], seas=state['names'],
                levels_m={k: state['levels'][k] for k in state['names']},
                grid=dict(pixels_per_degree=dr.PPD, rows=geo['h'].shape[0], cols=geo['h'].shape[1]),
                shares_m3s={k: round(v, 1) for k, v in gov['share'].items()},
                climate=dict(run=RUN, years=YEARS, climatology=CLIMATOLOGY, budget=budget_label,
                             open_water_evaporation_factor=evaporation),
                reading_rule=('As design_routing.npz: cells index rows from 90 N and columns from 0 E, 16 per degree, on '
                              'today\'s Moon; receiver, discharge_m3s and drained_km2 per cell, sea the index into seas, '
                              'lake_cell with lake_level_m and lake_spills. The weirs and cuts of 16 November divide '
                              'Moscoviense\'s water on the climate named, and the branches\' water runs in the rivers '
                              'it joins.'))
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, metadata=json.dumps(meta), receiver=r['receiver'].astype(np.int32),
                        discharge_m3s=r['discharge'].astype(np.float32),
                        drained_km2=(r['drained_m2'] / 1e6).astype(np.float32), sea=state['sea_of'],
                        lake_cell=cells, lake_level_m=r['lake_level'][cells].astype(np.float32),
                        lake_spills=~r['closed'][r['label'][cells]])
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
    started = time.time()
    clim_path = carry()
    log(f'climatology on today\'s ground: {clim_path.name}, {time.time() - started:.0f} s')
    surf = RD.load_surface(log)
    runoff, lake_net, precip, evap = load(1.0)
    geo = RD.geometry(surf, runoff, lake_net, log)
    rims = RD.rims_against_the_first_step(geo, surf['case'])
    if not all(v['agree'] for v in rims.values()):
        raise SystemExit('the rims do not reproduce the first step\'s design case')
    frame_path, record = frame_record()
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    cols = geo['h'].shape[1]
    area = np.repeat((geo['cell_km'] * 1e3) ** 2 * np.cos(np.radians(geo['lat'])), cols)
    closing = budget(geo, runoff, lake_net, precip, evap, area)
    factor = closing['open_water_evaporation_factor_closing_it']
    log(f'the routing\'s budget: {closing}, {time.time() - started:.0f} s')
    del runoff, lake_net, precip, evap
    cases, products = {}, {}
    branches = info = None
    for label, evaporation in zip(BUDGETS, (1.0, factor)):
        runoff, lake_net, _, _ = load(evaporation)
        state = RD.design_state(geo, runoff, lake_net, log=log)
        design = RD.design_record(state, geo, rotation)
        if branches is None:                          # the cuts' routes stand on the floor and the routing's lakes
            branches, info = GD.branches_for(geo, state, rotation, log)
        gov = GD.govern(geo, state, runoff, lake_net, branches, log=log)
        loss_field = -lake_net.ravel().astype(np.float64) * dr.MM_DAY * area
        rec = GD.governed_record(geo, state, gov, info, rotation, loss_field, area)
        rec['evaporation'] = evaporation
        rec['design_state_without_the_works'] = {k: {kk: v[kk] for kk in ('level_m', 'area_km2', 'loss_m3s',
                                                                          'own_rivers_m3s', 'supply_m3s', 'balance_m3s')}
                                                 for k, v in design.items()}
        rec['lakes_the_works_fill'] = GD.new_lakes(state['r'], gov['r'], area, cols, rotation)
        products[label] = product(gov, state, geo, label, evaporation, routing_path(label))
        cases[label] = rec
        log(f'{label}: near side {rec["seas"]["near_side"]}, {time.time() - started:.0f} s')
        del state, gov, runoff, lake_net
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (ROOT / RD.SURFACE, frame_path, GD.OUT)},
        ignored_inputs={rel(clim_path): sha256(clim_path), TURNED_CLIMATOLOGY: sha256(ROOT / TURNED_CLIMATOLOGY),
                        str(surf['grids']): surf['case']['grids']['sha256']},
        climatology=dict(run=RUN, years=YEARS, frame=rel(frame_path)),
        evidence=('A calculation on products: the drainage network study\'s routing of the climate run on the final '
                  'surface (sunlight 0.930, years 35-49, 300.3 K at the ground), carried from its settled turned frame '
                  'onto today\'s ground, over the redrawn floor with the works chosen on 16 November, each weir set by '
                  'routing again until the sea it serves receives its need. The cuts\' routes are found between 1.9 km '
                  'cells. Not validated against any observation.'),
        reading_rule=('As governed_drainage.json: levels in m above the redrawn common level, flows in m3/s, areas in '
                      'km2; cases on the routing\'s own budget and on the one its open-water evaporation factor closes. '
                      'design_state_without_the_works is the same climate with the upper seas spilling over their '
                      'natural passes and no weirs or cuts. against_the_present_climate pairs each number as [present '
                      'design climate at 8 degrees, this run]. routing names each budget\'s governed routing on the '
                      'research drive.'),
        budget=closing, rims_against_the_first_step=rims,
        branches={k: dict(cuts=v['cuts'], opens_onto=v['onto'], cuts_least_in_sum=v['least_in_sum']) for k, v in info.items()},
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
