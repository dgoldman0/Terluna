"""Where the design surface's rivers drop steeply: the falls and cascades of the design routing, ranked by the power of
the falling water, for the globe's marks and white water.

    ulimit -v 10000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.falls   # about a minute

Each step of the design routing from a cell to its receiver drops from the cell's water surface, the ground on dry
cells and the level on a lake or sea, to the receiver's, on the design surface the globe draws (globe_surface.py: the
redrawn floor with the works' cuts, every water body at its level). Steps on rivers of at least MIN_Q steeper than
STEEP join along the flow into reaches, each starting where no steep step flows into it and running down while the
steps stay steep. Each reach gets its drop, its length, its mean flow, its steepest step and the power of its falling
water (density times lunar gravity times flow times drop), the place in today's frame and the globe's, the water it
ends in, and the named feature beside it (the IAU gazetteer: the smallest crater of 20 km or more within three quarters
of its diameter, else the nearest). The river that carries Moscoviense's overflow over its natural pass to the salt sea
is traced from its source and marked on each reach it holds.

The 1.9 km cells average the ground: a reach's steepest step is the steepest 1.9 km, and whether the water falls sheer
or cascades is finer than the heights. results/falls.json holds the reaches, largest first; the product on the research
drive holds every steep step's cell, steepness and flow for the globe's white water.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from geography import nomenclature as nm
from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import globe_surface as GS
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'falls.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/falls')
PRODUCT = BULK / 'falls_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-falls/1'
CODE = ('research/studies/atlas_round/falls.py', 'research/studies/atlas_round/globe_surface.py',
        'research/studies/atlas_round/design_routing.py')
MIN_Q = 100.0                        # m3/s: the rivers searched
STEEP = 0.03                         # a step dropping 30 m a km or more joins a reach
LISTED = 40                          # reaches kept in the record, largest power first
G_MOON, RHO = 1.62, 1000.0
MOON_KM = 1737.4
NAMED_KM = 20.0                      # named features at least this wide name a reach


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reaches(receiver, q, surface, water, cols, ppd=16, min_q=MIN_Q, steep=STEEP):
    """The steep steps (cell, slope) and the reaches they join into, as dicts of cells, drop, length, mean flow,
    steepest slope and power, on a grid of `cols` columns at `ppd` cells a degree (rows from 90 N)."""
    n = receiver.size
    idx = np.arange(n)
    lat = 90.0 - (idx // cols + 0.5) / ppd
    lon = (idx % cols + 0.5) / ppd
    rlat, rlon = lat[receiver], lon[receiver]
    d_km = MOON_KM * np.radians(np.hypot(rlat - lat, ((rlon - lon + 180.0) % 360.0 - 180.0) * np.cos(np.radians(lat))))
    drop = surface - surface[receiver]
    moving = (receiver != idx) & (q >= min_q) & (d_km > 0) & ~(water & water[receiver])
    slope = np.where(moving, drop / np.maximum(d_km * 1000.0, 1.0), 0.0)
    is_steep = moving & (slope >= steep)
    fed = np.zeros(n, bool)
    fed[receiver[is_steep]] = True
    heads = np.flatnonzero(is_steep & ~fed)
    seen = np.zeros(n, bool)
    out = []
    for cell in heads[np.argsort(-q[heads])]:
        path = [int(cell)]
        while is_steep[path[-1]] and not seen[path[-1]] and len(path) < 400:
            seen[path[-1]] = True
            path.append(int(receiver[path[-1]]))
        if len(path) < 2:
            continue
        steps = np.array(path[:-1])
        total = float(surface[path[0]] - surface[path[-1]])
        flow = float(q[steps].mean())
        out.append(dict(cells=path, drop_m=total, length_km=float(d_km[steps].sum()), flow_m3s=flow,
                        steepest=float(slope[steps].max()), power_gw=RHO * G_MOON * flow * total / 1e9))
    out.sort(key=lambda x: -x['power_gw'])
    steep_cells = np.flatnonzero(is_steep)
    return steep_cells, slope[steep_cells], out


def namer(features=None):
    """The named feature beside a place: the smallest crater or other feature of NAMED_KM or more whose three quarters
    of a diameter reach it, else the nearest such feature's name with 'near'."""
    feats = [f for f in (features if features is not None else nm.features()) if f.get('diameter_km', 0) >= NAMED_KM]
    la_f = np.array([f['lat'] for f in feats])
    lo_f = np.array([f['lon'] for f in feats])
    dia = np.array([f['diameter_km'] for f in feats])

    def name(lat, lon):
        d = MOON_KM * np.radians(np.hypot(la_f - lat, ((lo_f - lon + 180.0) % 360.0 - 180.0) * np.cos(np.radians(lat))))
        inside = d < 0.75 * dia
        if inside.any():
            k = np.flatnonzero(inside)[np.argmin(dia[inside])]
            return dict(name=feats[k]['name'], beside=True, km=round(float(d[k])))
        k = int(np.argmin(d))
        return dict(name=feats[k]['name'], beside=False, km=round(float(d[k])))
    return name


def build(log=print):
    started = time.time()
    routing = DR.load(log)
    meta = routing['meta']
    with np.load(GS.PRODUCT, allow_pickle=False) as z:
        smeta = json.loads(str(z['metadata']))
        h = z['height_m'].astype(np.float64).ravel()
        wcell, wlevel, wbody = z['water_cell'], z['water_level_m'], z['water_body']
    cols = smeta['grid']['cols']
    rotation = np.asarray(smeta['frame']['rotation_today_to_turned'], float)
    surface = h.copy()
    surface[wcell] = np.maximum(wlevel, h[wcell])
    water = np.zeros(h.size, bool)
    water[wcell] = True
    body = np.zeros(h.size, np.int8)
    body[wcell] = wbody
    receiver = routing['receiver'].astype(np.int64)
    q = routing['discharge_m3s']
    steep_cells, steep_slope, found = reaches(receiver, q, surface, water, cols)
    log(f'{len(found)} steep reaches on rivers of {MIN_Q:.0f} m3/s or more, {time.time() - started:.0f} s')

    # the river over Moscoviense's natural pass: from the cell where the weir's share is put in, down to the sea
    source = None
    share = None
    try:
        gd = json.loads((HERE / 'results' / 'final_drainage.json').read_text())
        share = gd['cases'][DR.BUDGET]['weir_history'][-1]['shares_m3s']['south_pole_aitken']
    except (OSError, KeyError, IndexError):
        pass
    big = [x for x in found if share and abs(x['flow_m3s'] - share) / share < 0.2]
    if big:
        cell = big[0]['cells'][0]
        donors = np.flatnonzero(receiver == cell)
        while donors.size and q[donors].max() >= 0.8 * q[cell]:
            cell = int(donors[np.argmax(q[donors])])
            donors = np.flatnonzero(receiver == cell)
        source = cell
    great = []
    if source is not None:
        k = source
        while not water[k] or body[k] == 0:
            great.append(k)
            if receiver[k] == k or len(great) > 50000:
                break
            k = int(receiver[k])
        great.append(k)
    on_great = set(great)
    name = namer()
    bodies = smeta['bodies']
    listed = []
    for x in found[:LISTED]:
        start, end = x['cells'][0], x['cells'][-1]
        at = NW.place(start, cols, rotation)
        listed.append(dict(
            drop_m=round(x['drop_m']), length_km=round(x['length_km'], 1), flow_m3s=round(x['flow_m3s']),
            steepest=round(x['steepest'], 3), steepest_deg=round(float(np.degrees(np.arctan(x['steepest']))), 1),
            power_gw=round(x['power_gw'], 2), start=start, end=end, at=at,
            into=(bodies[body[end]] if water[end] and body[end] > 0 else 'a lake' if water[end] else 'land'),
            beside=name(*at['today']), on_moscoviense_river=start in on_great))
    great_record = None
    if great:
        def step_km(a, b):
            la, lb = 90.0 - (a // cols + 0.5) / 16.0, 90.0 - (b // cols + 0.5) / 16.0
            dlon = ((b % cols - a % cols) / 16.0 + 180.0) % 360.0 - 180.0
            return MOON_KM * np.radians(np.hypot(lb - la, dlon * np.cos(np.radians(la))))
        dist = sum(step_km(a, b) for a, b in zip(great[:-1], great[1:]))
        great_record = dict(source=NW.place(great[0], cols, rotation), mouth=NW.place(great[-1], cols, rotation),
                            flow_m3s=round(float(q[great[0]])), length_km=round(float(dist)),
                            falls_m=round(float(surface[great[0]] - surface[great[-1]])),
                            into=bodies[body[great[-1]]] if water[great[-1]] else 'land')
    BULK.mkdir(parents=True, exist_ok=True)
    reach_of = np.full(h.size, -1, np.int16)
    for i, x in enumerate(found[:LISTED]):
        reach_of[x['cells'][:-1]] = i
    np.savez_compressed(PRODUCT, metadata=json.dumps(dict(schema=SCHEMA + '-grid', grid=smeta['grid'], min_q_m3s=MIN_Q,
                                                         steep=STEEP)),
                        steep_cell=steep_cells.astype(np.int32), steep_slope=steep_slope.astype(np.float32),
                        steep_flow_m3s=q[steep_cells].astype(np.float32), steep_reach=reach_of[steep_cells])
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)   # noqa: E731
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (DR.OUT, GS.OUT)},
        ignored_inputs={str(DR.PRODUCT): sha256(DR.PRODUCT), str(GS.PRODUCT): sha256(GS.PRODUCT)},
        evidence=('Steps of the design routing over the design surface between 1.9 km cells: where the large rivers '
                  'drop steeply, not whether they fall sheer, which is finer than the heights.'),
        reading_rule=('Drops in m, lengths in km, flows in m3/s, power in GW of falling water at lunar gravity; '
                      'steepest is the steepest 1.9 km step as a grade. at.today and at.turned are latitude and '
                      'longitude in today\'s frame and the globe\'s. product holds every steep step for the white water, '
                      'steep_reach the index into reaches of the listed reach each belongs to (-1 for none).'),
        settings=dict(min_q_m3s=MIN_Q, steep=STEEP, named_km=NAMED_KM),
        summary=dict(steep_reaches=len(found), steep_steps=int(steep_cells.size),
                     over_1_gw=sum(x['power_gw'] >= 1.0 for x in found),
                     over_500_m=sum(x['drop_m'] >= 500.0 for x in found), moscoviense_river=great_record),
        reaches=listed, product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)), seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    record = build(lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1) + '\n')
    print(json.dumps(record['summary']))
    for x in record['reaches'][:12]:
        print(f"{x['power_gw']:6.1f} GW {x['flow_m3s']:7,d} m3/s {x['drop_m']:5d} m / {x['length_km']:4.1f} km "
              f"into {x['into']:24s} {'beside' if x['beside']['beside'] else 'near'} {x['beside']['name']}"
              f"{'  (Moscoviense river)' if x['on_moscoviense_river'] else ''}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
