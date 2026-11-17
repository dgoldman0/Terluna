"""Dryland rivers on the final surface: where gravity could carry water through the arid and semi-arid land of the
climate run on the final surface to an outlet, from the full seas' surplus and from the lakes of the large rivers, with
each route's cut, its length through dryland, the green corridor it could keep, that corridor's loss, and what the
near-side sea gives up.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.dryland_rivers   # a few minutes

The land's aridity is the moisture map's reading of the run at sunlight 0.930, years 35-49 (rain over Priestley and
Taylor's demand on the run's T21 grid, drylands.py's window), carried onto each 1.9 km cell of today's ground at its
place in the run's settled frame. A source is a full sea at its rim (Moscoviense, Orientale, Smythii-Marginis) or a
lake through which at least Q_SOURCE flows in the closed budget's governed routing (final_drainage.py). A route leaves
a source over one pass of the catchment graph into a neighbouring catchment: a channel at the source's level, cut down
through the pass and the ground beyond it until that ground falls below the source's level, after which the water runs
down the routing by itself. Its new river is the stretch from the pass to the first river already carrying Q_JOINS,
or to the sea or closed lake it ends in, whichever comes first; its dryland is the arid and semi-arid land on that
stretch.

Each route carries Q_ROUTE (less where its source has less to give). Its green corridor is river_losses.py's: the
channel's water surface, widened for lunar gravity, and a belt of plants drinking the bed's seepage, at the narrow and
the wide ends of their ranges, over the new river's dryland; its loss is that corridor at the open water's net loss on
the closed budget (the climate's rain less the open water's evaporation, where that is negative). The ranking is green
land per metre of cut, the cut being the deepest the channel needs. Moscoviense's weirs are fixed crests (register,
17 November), each taking a fixed share of its overflow, and Orientale and Smythii-Marginis pass on all they receive, so
a route's loss falls on the sea its water was heading for: the near-side sea for routes ending in it, in Orientale or in
Smythii-Marginis; the salt sea or Humboldtianum for routes ending in them; and the three by Moscoviense's shares for routes
ending in Moscoviense. A route into a closed lake takes its whole flow from the sea its source's water reached. Each
sea's fall is its level's response to that much less supply: the parted near-side sea's south part's from
near_side_neck.py, the salt sea's and Humboldtianum's from their own levels' tables at their need.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import moisture_map as bm
from climate.gcm import turned_climatology as TC
from climate.gcm.turned import to_turned
from geography import drainage as dr
from research.studies.atlas_round import drylands as DL
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import river_losses as RL
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'dryland_rivers.json'
SCHEMA = 'terluna.research.atlas-round-dryland-rivers/1'
CODE = ('research/studies/atlas_round/dryland_rivers.py', 'research/studies/atlas_round/final_drainage.py',
        'research/studies/atlas_round/river_losses.py', 'research/studies/atlas_round/redrawn_drainage.py',
        'biosphere/ecology/moisture_map.py')
BUDGET = 'closed_budget'
FULL_SEAS = ('moscoviense', 'orientale', 'smythii_marginis')
Q_SOURCE = 1000.0                    # m3/s through a lake for it to be a source
Q_ROUTE = 500.0                      # m3/s a route carries, less where its source has less to give
Q_JOINS = 1000.0                     # m3/s: a river this large already keeps its own corridor
CUT_MAX_M = 300.0                    # the deepest cut considered
DRY_MIN_KM = 30.0                    # the least arid and semi-arid land on a new river's stretch
LISTED = 30
SAME_GAP_KM = 25.0                   # passes nearer than this to a better route's are the same gap
NEAR_SIDE = 'near_side'
PASS_ON = ('near_side', 'orientale', 'smythii_marginis')   # seas that pass on all they receive to the near side
HELD = ('south_pole_aitken', 'humboldtianum')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------- the land's aridity on today's ground

def aridity_on_ground(rows, cols, rotation, log=print):
    """The UNEP class of each 1.9 km cell (0 hyper-arid to 4 humid): the run's aridity index on its T21 grid, bilinear
    at the cell's place in the run's frame. Returns the classes (uint8, rows x cols) and the T21 reading."""
    bm.RUNS = DL.RUNS
    bm.CACHE = DL.CACHE
    label, run, years, _ = DL.WINDOWS[0]
    c = bm.read_run(run, years, humidity=False)
    lat_src, lon_src, ai = np.asarray(c['lat'], float), np.asarray(c['lon'], float), np.asarray(c['ai'], float)
    if lat_src[0] < lat_src[-1]:
        lat_src, ai = lat_src[::-1], ai[::-1]
    out = np.empty((rows, cols), np.uint8)
    lon = (np.arange(cols) + 0.5) / dr.PPD
    for i in range(0, rows, 180):
        la = 90.0 - (np.arange(i, min(i + 180, rows)) + 0.5) / dr.PPD
        g_lat, g_lon = np.meshgrid(la, lon, indexing='ij')
        t_lat, t_lon = to_turned(g_lat, g_lon, rotation)
        out[i:i + la.size] = bm.unep_class(TC.sample(ai, lat_src, lon_src, t_lat, t_lon))
    log(f'aridity on today\'s ground from {run} years {years}')
    return out, dict(run=run, years=list(years), window=label)


# ---------------------------------------------------------------- the routing's stretches

def step_km(receiver, rows, cols, lat, cell_km):
    return RL.flow_length_km(receiver, rows, cols, lat, cell_km)


def downstream_sum(receiver, value, stop, rounds=40):
    """For each cell, the sum of `value` from the cell down the routing to the first stop cell (the stop's own value
    left out), by doubling: each round adds what lies beyond the cells already summed."""
    total = np.where(stop, 0.0, value).astype(np.float64)
    nxt = np.where(stop, np.arange(receiver.size), receiver).astype(np.int64)
    for _ in range(rounds):
        nn = nxt[nxt]
        if np.array_equal(nn, nxt):
            break
        total = total + np.where(stop[nxt], 0.0, total[nxt])
        nxt = np.where(stop[nxt], nxt, nn)
    return total, nxt


def stretches(r, h_flat, dry, wet, step, outlet):
    """For each cell, the new river's stretch below it: its dryland (km of arid and semi-arid land) and its length (km)
    down to the first cell of a river carrying Q_JOINS or an outlet, and that end cell."""
    stop = (r['discharge'] >= Q_JOINS) | outlet
    dry_km, end = downstream_sum(r['receiver'], np.where(dry & ~wet, step, 0.0), stop)
    length_km, _ = downstream_sum(r['receiver'], step, stop)
    return dry_km, length_km, end


# ---------------------------------------------------------------- the sources and their routes

def water_surface(meta, sea, lake_cell, lake_level, h_flat):
    """Each cell's water level: its sea's, its lake's, NaN on dry ground."""
    surface = np.full(h_flat.size, np.nan)
    for i, k in enumerate(meta['seas']):
        surface[sea == i] = meta['levels_m'][k]
    surface[lake_cell] = lake_level
    return surface


def routes(geo_graph, h_flat, surface, r, sea, names, term, dry_km, length_km, end, rate_m2s, log=print):
    """Every route over one pass from a source's water into the next catchment, its far side cut until the ground falls
    below the source's level, kept where its new river crosses at least DRY_MIN_KM of arid and semi-arid land; the
    shallowest to each end from each source."""
    label, sinks, (u, v, w, cu, cv) = geo_graph
    n_nodes = sinks.size
    pit = sinks.astype(np.int64)
    level = np.where(np.isnan(surface[pit]), h_flat[pit], surface[pit])
    wet = ~np.isnan(surface[pit])
    outflow = np.zeros(n_nodes)
    np.maximum.at(outflow, label.ravel(), r['discharge'])
    node_sea = sea[pit]
    full = {names.index(k) for k in FULL_SEAS}
    is_full = np.isin(node_sea, list(full))
    is_lake = wet & (node_sea < 0) & (outflow >= Q_SOURCE)
    source = is_full | is_lake
    a = np.concatenate([u, v]).astype(np.int64)
    b = np.concatenate([v, u]).astype(np.int64)
    ww = np.concatenate([w, w]).astype(np.float64)
    ca = np.concatenate([cu, cv]).astype(np.int64)
    cb = np.concatenate([cv, cu]).astype(np.int64)
    cut0 = ww - level[a]
    keep = source[a] & (cut0 > RD.CUT_COUNTED_M) & (cut0 <= CUT_MAX_M)
    same = wet[b] & (np.abs(level[b] - level[a]) < 0.01)                     # the same water over two catchments
    keep &= ~same & ~(wet[b] & (level[b] >= level[a]))                       # water standing higher beyond: no way on
    keep &= ~((node_sea[b] >= 0) & (node_sea[b] == node_sea[a]))
    idx = np.flatnonzero(keep)
    log(f'{source.sum():,} sources ({is_full.sum():,} catchments under full seas, {is_lake.sum():,} lakes), '
        f'{idx.size:,} passes within {CUT_MAX_M:g} m')
    found = {}
    rec = r['receiver']
    for e in idx:
        lev = level[a[e]]
        c = int(cb[e])
        top = ww[e]
        for _ in range(5000):                         # the far side cut until the ground falls below the source's level
            if h_flat[c] <= lev or not np.isnan(surface[c]):
                break
            top = max(top, float(h_flat[c]))
            nxt = int(rec[c])
            if nxt == c:
                break
            c = nxt
        else:
            continue
        if top - lev > CUT_MAX_M or dry_km[c] < DRY_MIN_KM:
            continue
        group = ('sea', int(node_sea[a[e]])) if is_full[a[e]] else ('lake', round(float(lev), 1), int(term[pit[a[e]]]))
        key = (group, int(end[c]))
        cut = top - lev
        if key not in found or cut < found[key]['cut_m']:
            found[key] = dict(cut_m=cut, source_node=int(a[e]), source_level=float(lev), source_pit=int(pit[a[e]]),
                              pass_cell=int(cb[e]),
                              shore_cell=int(ca[e]), free_cell=c, end_cell=int(end[c]),
                              source_outflow_m3s=float(outflow[a[e]]), full_sea=bool(is_full[a[e]]),
                              source_sea=int(node_sea[a[e]]), source_terminal=int(term[pit[a[e]]]))
    log(f'{len(found):,} routes with {DRY_MIN_KM:g} km of dryland or more')
    return list(found.values())


def walk_to(rec, start, stop_cell, limit=200000):
    """The cells from `start` down the routing to `stop_cell` (or to an outlet)."""
    path, c = [int(start)], int(start)
    for _ in range(limit):
        if c == stop_cell:
            break
        n = int(rec[c])
        if n == c:
            break
        c = n
        path.append(c)
    return np.array(path, np.int64)


def own_rivers(fd_case, rec, sea, discharge, cols, names, source_sea, source_pit, full_sea):
    """The cells of the river a source's own water runs down: a lake's from its outlet; a full sea's overflow from the
    first cell past its pass (and Moscoviense's branches from where each runs on by itself), each as the governed
    routing leads it."""
    starts = []
    if not full_sea:
        starts.append(int(rec[source_pit]))
    else:
        k = names[source_sea]
        rivers = fd_case['rivers']
        places = [rivers[f'{k}_over_its_natural_pass']['leaves']['today']] if f'{k}_over_its_natural_pass' in rivers else []
        for key, v in rivers.items():
            if key.startswith(f'{k}_to_'):
                places += [seg['starts']['today'] for seg in v['segments']]
        for la, lo in places:
            i, j = int((90.0 - la) * dr.PPD), int((lo % 360.0) * dr.PPD)
            box = [(a % (180 * dr.PPD)) * cols + (b % cols) for a in range(i - 3, i + 4) for b in range(j - 3, j + 4)]
            box = [x for x in box if sea[x] < 0] or box
            starts.append(max(box, key=lambda x: discharge[x]))
    cells = set()
    for c in starts:
        cells.update(walk_to(rec, c, -1).tolist())
    return cells


def classify(route, rec, surface, sea, names, discharge, cols, fd_case):
    """Whether a route's new river runs back into its own source (a loop with no outlet of its own), rejoins the river
    its source's water already runs down (a side channel), or reaches another river, a sea or a closed lake."""
    path = walk_to(rec, route['free_cell'], route['end_cell'])
    if route['full_sea']:
        loop = bool(np.any(sea[path] == route['source_sea']))
    else:
        lev = route['source_level']
        loop = bool(np.any((sea[path] < 0) & (np.abs(np.nan_to_num(surface[path], nan=-1e9) - lev) < 0.01)))
    if loop:
        return 'loops back into its source'
    own = own_rivers(fd_case, rec, sea, discharge, cols, names, route['source_sea'], route['source_pit'], route['full_sea'])
    if route['end_cell'] in own:
        return 'rejoins its source\'s own river'
    return 'reaches another river, a sea or a closed lake'


def held_seas_falls(graph, h_flat, sea, names, area, loss_field, level0, case, step=1000.0):
    """How far the salt sea and Humboldtianum fall for each m3/s less supply at their need (network.Table on the level
    at which each cell joins the water over the sea's deepest cell, as redrawn_drainage.lake_at reads it)."""
    label, sinks, edges = graph
    out = {}
    for k in HELD:
        cells = np.flatnonzero(sea == names.index(k))
        seed = int(cells[np.argmin(h_flat[cells])])
        joined = NW.joining_level(h_flat, label, sinks.size, edges, seed)
        table = NW.Table(joined, area, loss_field, h_flat, level0 + 1.0)
        need = case['seas'][k]['loss_m3s']
        lev = table.level_for(need - step, level0 + 1e-6)
        out[k] = round((level0 - (level0 if lev is None else lev)) / step, 4)
    return out


def gap_km(a, b, cols):
    la1, lo1 = 90.0 - (a // cols + 0.5) / dr.PPD, (a % cols + 0.5) / dr.PPD
    la2, lo2 = 90.0 - (b // cols + 0.5) / dr.PPD, (b % cols + 0.5) / dr.PPD
    dlon = ((lo2 - lo1 + 180.0) % 360.0) - 180.0
    return math.radians(math.hypot(la2 - la1, dlon * math.cos(math.radians(0.5 * (la1 + la2))))) * 1737.4


# ---------------------------------------------------------------- the record

def place(cell, cols, rotation):
    return NW.place(int(cell), cols, rotation)


def where_it_ends(cell, r, sea, names, term):
    """What a new river reaches: a sea, a closed lake, or a large river (and the sea or closed lake that river ends in)."""
    t = int(term[cell])
    goes = names[sea[t]] if sea[t] >= 0 else 'a closed lake'
    if sea[cell] >= 0:
        return dict(reaches=names[sea[cell]], ends_in=names[sea[cell]])
    if r['receiver'][cell] == cell:
        return dict(reaches='a closed lake', ends_in='a closed lake')
    return dict(reaches=f'a river of {r["discharge"][cell]:,.0f} m3/s', ends_in=goes)


def bearers(sea_name, flow, shares):
    """Which seas give up a flow lost on its way to `sea_name` (m3/s each): the near-side sea for itself and the seas
    that pass all they receive on to it, a held sea for itself, the three by Moscoviense's fixed shares for Moscoviense."""
    if sea_name in PASS_ON:
        return {NEAR_SIDE: flow}
    if sea_name in HELD:
        return {sea_name: flow}
    if sea_name == 'moscoviense':
        return {NEAR_SIDE: flow * shares['smythii_marginis'], 'south_pole_aitken': flow * shares['south_pole_aitken'],
                'humboldtianum': flow * shares['humboldtianum']}
    return {}


def evaluate(route, r, sea, names, term, dry_km, length_km, rate_m2s, cols, rotation, falls, shares, level0=0.0):
    """A route's new river, green corridor and loss, and which seas give up what and fall how far. `falls` holds each
    bearing sea's fall per m3/s of supply (the near-side sea's for its parted south part), `shares` Moscoviense's."""
    c = route['free_cell']
    q = min(Q_ROUTE, 0.5 * route['source_outflow_m3s']) if not route['full_sea'] else Q_ROUTE
    out = dict(cut_m=round(route['cut_m'], 1), source=('full sea ' + names[route['source_sea']]) if route['full_sea']
               else f'lake {route["source_level"] - level0:+.0f} m against the common level, passing {route["source_outflow_m3s"]:,.0f} m3/s',
               pass_at=place(route['pass_cell'], cols, rotation), shore_at=place(route['shore_cell'], cols, rotation),
               flow_m3s=round(q), new_river_km=round(float(length_km[c])), dryland_km=round(float(dry_km[c])),
               **where_it_ends(route['end_cell'], r, sea, names, term), ends_at=place(route['end_cell'], cols, rotation))
    if route['full_sea']:
        own = names[route['source_sea']]                  # an upper sea's water: Moscoviense's shares or on to the near side
    else:
        t = int(route['source_terminal'])
        own = names[sea[t]] if sea[t] >= 0 else 'a closed lake'
    out['source_water_ends_in'] = own
    green, loss = {}, {}
    for reading, coef, belt in (('narrow', RL.WIDTH[0], RL.RIPARIAN_KM[0]), ('wide', RL.WIDTH[1], RL.RIPARIAN_KM[1])):
        width_m = float(RL.channel_width_m(np.array([q]), coef)[0]) + belt * 1e3
        lost = width_m * float(rate_m2s[c])
        share = min(1.0, q / lost) if lost > 0 else 1.0     # a river that loses all it carries greens only so far
        green[reading] = round(out['dryland_km'] * width_m / 1e3 * share)
        loss[reading] = round(min(lost, q))
    out['green_km2'] = green
    out['loss_m3s'] = loss
    out['green_km2_per_m_of_cut'] = round(green['narrow'] / max(route['cut_m'], 1.0), 1)
    ends = out['ends_in']
    cost, lowers = {}, {}
    for reading in ('narrow', 'wide'):
        given = bearers(own, q, shares) if ends == 'a closed lake' else bearers(ends, loss[reading], shares)
        cost[reading] = {k: round(v, 1) for k, v in given.items()}
        lowers[reading] = {k: round(v * falls[k], 2) for k, v in given.items()}
    out['seas_give_up_m3s'] = cost
    out['seas_fall_m'] = lowers
    return out


def drylands_by_band(cls, area, land, rotation, rows, cols):
    """The arid and semi-arid land (km2) by 15 degrees of turned latitude, on today's ground at 1.9 km cells."""
    dry = land & ((cls == 1) | (cls == 2))
    out = {}
    lon = (np.arange(cols) + 0.5) / dr.PPD
    for i in range(0, rows, 180):
        la = 90.0 - (np.arange(i, min(i + 180, rows)) + 0.5) / dr.PPD
        g_lat, g_lon = np.meshgrid(la, lon, indexing='ij')
        t_lat, _ = to_turned(g_lat, g_lon, rotation)
        band = (np.floor(np.abs(t_lat) / 15.0) * 15).astype(int)
        d = dry[i * cols:(i + la.size) * cols].reshape(la.size, cols)
        a = area[i * cols:(i + la.size) * cols].reshape(la.size, cols)
        for b in np.unique(band[d]):
            out[f'{b}-{b + 15}'] = out.get(f'{b}-{b + 15}', 0.0) + float(a[d & (band == b)].sum()) / 1e6
    return {k: round(v) for k, v in sorted(out.items(), key=lambda x: int(x[0].split('-')[0]))}


def build(log=print):
    started = time.time()
    fd = json.loads(FD.OUT.read_text())
    neck = json.loads((HERE / 'results' / 'near_side_neck.json').read_text())
    frame_path, record = FD.frame_record()
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    surf = RD.load_surface(log)
    h, level0, lat, cell_km = surf['h'], surf['level0'], surf['lat'], surf['cell_km']
    rows, cols = h.shape
    h_flat = h.ravel()
    with np.load(FD.routing_path(BUDGET), allow_pickle=False) as d:
        meta = json.loads(str(d['metadata']))
        r = dict(receiver=d['receiver'].astype(np.int64), discharge=d['discharge_m3s'].astype(np.float64))
        sea, lake_cell, lake_level = d['sea'], d['lake_cell'].astype(np.int64), d['lake_level_m'].astype(np.float64)
    names = meta['seas']
    surface = water_surface(meta, sea, lake_cell, lake_level, h_flat)
    wet = ~np.isnan(surface)
    idx = np.arange(h_flat.size)
    outlet = (r['receiver'] == idx) | (sea >= 0)
    term = dr.terminal(r['receiver'])
    cls, reading = aridity_on_ground(rows, cols, rotation, log)
    cls = cls.ravel()
    dry = (cls == 1) | (cls == 2)
    step = step_km(r['receiver'], rows, cols, lat, cell_km)
    _, lake_net, _, _ = FD.load(fd['cases'][BUDGET]['evaporation'])
    rate = np.maximum(-lake_net.ravel().astype(np.float64), 0.0) * dr.MM_DAY      # m/s the open water loses beyond rain
    area = np.repeat((cell_km * 1e3) ** 2 * np.cos(np.radians(lat)), cols)
    loss_field = -lake_net.ravel().astype(np.float64) * dr.MM_DAY * area
    del lake_net
    dry_km, length_km, end = stretches(r, h_flat, dry, wet, step, outlet)
    stop = (r['discharge'] >= Q_JOINS) | outlet
    rate_m2s, _ = downstream_sum(r['receiver'], np.where(dry & ~wet, step * 1e3 * rate, 0.0), stop)
    log(f'stretches, {time.time() - started:.0f} s')
    graph = RD.catchment_graph(h, lat, cell_km)
    found = routes(graph, h_flat, surface, r, sea, names, term, dry_km, length_km, end, rate_m2s, log)
    case = fd['cases'][BUDGET]
    held_falls = held_seas_falls(graph, h_flat, sea, names, area, loss_field, level0, case)
    log(f'the held seas fall {held_falls} m for each m3/s less, {time.time() - started:.0f} s')
    del graph
    rows_ = neck['budgets'][BUDGET]['natural_neck']['with_extra_loss']
    fall = (rows_[-1]['south_level_m'] - rows_[0]['south_level_m']) / (rows_[-1]['extra_loss_m3s'] - rows_[0]['extra_loss_m3s'])
    fall = -fall
    falls = {NEAR_SIDE: fall, **{k: held_falls[k] for k in HELD}}
    shares = {k: v / case['moscoviense_surplus_m3s'] for k, v in case['released_m3s'].items()}
    kinds = [classify(x, r['receiver'], surface, sea, names, r['discharge'], cols, case) for x in found]
    loops = sum(k.startswith('loops') for k in kinds)
    pairs = []
    for x, kind in zip(found, kinds):
        if kind.startswith('loops'):
            continue
        e = evaluate(x, r, sea, names, term, dry_km, length_km, rate_m2s, cols, rotation, falls, shares, level0)
        e['kind'] = kind
        pairs.append((e, x['pass_cell']))
    pairs.sort(key=lambda p: -p[0]['green_km2_per_m_of_cut'])
    evaluated, kept_passes, merged = [], [], 0
    for e, cell in pairs:
        if any(gap_km(cell, k, cols) < SAME_GAP_KM for k in kept_passes):
            merged += 1
            continue
        evaluated.append(e)
        kept_passes.append(cell)
    log(f'{loops} loops dropped, {merged} routes in a gap already listed, {len(evaluated)} routes')
    land = ~wet
    by_source = {}
    for x in evaluated:
        by_source.setdefault(x['source'] if x['source'].startswith('full sea') else 'lakes of the large rivers', []).append(x)
    top = evaluated[:LISTED]
    best_each_source = {k: v[:5] for k, v in by_source.items()}
    summary = dict(
        routes=len(evaluated), loops_dropped=loops, merged_into_a_gap_listed=merged,
        routes_by_kind={k: sum(1 for x in evaluated if x['kind'] == k) for k in sorted({x['kind'] for x in evaluated})},
        all_routes_green_km2={k: sum(x['green_km2'][k] for x in evaluated) for k in ('narrow', 'wide')},
        all_routes_loss_m3s={k: sum(x['loss_m3s'][k] for x in evaluated) for k in ('narrow', 'wide')},
        all_routes_seas_give_up_m3s={k: {s_: round(sum(x['seas_give_up_m3s'][k].get(s_, 0.0) for x in evaluated), 1)
                                         for s_ in (NEAR_SIDE,) + HELD} for k in ('narrow', 'wide')},
        all_routes_seas_fall_m={k: {s_: round(sum(x['seas_fall_m'][k].get(s_, 0.0) for x in evaluated), 2)
                                    for s_ in (NEAR_SIDE,) + HELD} for k in ('narrow', 'wide')},
        routes_by_source={k: len(v) for k, v in by_source.items()},
        listed_green_km2={k: sum(x['green_km2'][k] for x in top) for k in ('narrow', 'wide')},
        listed_loss_m3s={k: sum(x['loss_m3s'][k] for x in top) for k in ('narrow', 'wide')},

        arid_and_semi_arid_land_km2=round(float(area[land & dry].sum()) / 1e6),
        land_km2=round(float(area[land].sum()) / 1e6),
        falls_m_per_1000_m3s={k: round(v * 1000.0, 1) for k, v in falls.items()},
        moscoviense_fixed_shares={k: round(v, 4) for k, v in shares.items()})
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (FD.OUT, HERE / 'results' / 'near_side_neck.json', ROOT / RD.SURFACE, frame_path)},
        ignored_inputs={rel(FD.routing_path(BUDGET)): sha256(FD.routing_path(BUDGET))},
        aridity=reading, budget=BUDGET,
        settings=dict(q_source_m3s=Q_SOURCE, q_route_m3s=Q_ROUTE, q_joins_m3s=Q_JOINS, cut_max_m=CUT_MAX_M,
                      dryland_least_km=DRY_MIN_KM, width_m_per_sqrt_m3s=list(RL.WIDTH), riparian_km=list(RL.RIPARIAN_KM)),
        evidence=('A calculation on products: the closed budget\'s governed routing of the run at sunlight 0.930 on the '
                  'final surface (final_drainage.py), the moisture map\'s aridity of that run carried onto today\'s '
                  'ground, the redrawn floor\'s catchment graph between 1.9 km cells, and river_losses.py\'s corridor. '
                  'A route\'s cut is the deepest its channel needs; tunnels, costs and the wider floodplains are left '
                  'out. Not validated against any observation.'),
        reading_rule=('Each route leaves a source over one pass; its new river runs from where the cut ends to the first '
                      'river of at least q_joins or its outlet, and dryland_km is the arid and semi-arid land on that '
                      'stretch. green_km2 is the channel and the belt over that land at the narrow and wide readings, '
                      'cut short where the river would lose all it carries; loss_m3s what that corridor loses to the air '
                      'beyond its rain. seas_give_up is that loss on the sea the route\'s water was heading for, with '
                      'Moscoviense\'s weirs fixed (the near-side sea for routes into it, Orientale or Smythii-Marginis; '
                      'a held sea for routes into it; Moscoviense\'s shares for routes into Moscoviense), or the whole '
                      'flow where it ends in a closed lake; seas_fall is how far that lowers each (the near-side sea\'s '
                      'parted south part). kind says whether the new river rejoins the river '
                      'its source\'s water already runs down (a side channel through dry land) or reaches another river, '
                      'a sea or a closed lake; routes running back into their own source are dropped, and of passes '
                      'within 25 km of each other only the best is listed. Places are on today\'s Moon and the settled '
                      'turned one.'),
        summary=summary, drylands_by_turned_latitude_km2=drylands_by_band(cls, area, land, rotation, rows, cols),
        routes=top, best_each_source=best_each_source, seconds=round(time.time() - started))


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
