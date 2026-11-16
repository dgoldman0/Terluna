"""The atlas round's second step: the drainage and the river network on the redrawn floor.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.redrawn_drainage     # about 40 minutes

The first step (redrawn_surface.py) redrew the surface under its water at the crust's settled state. This step routes
the design years' rain over the redrawn floor:

- the heights are LOLA's at 16 pixels per degree above the dry Moon's geoid with the settled case's change of the
  ground against the sea surface, U - N, added, and the common level is the redrawn one, below which 28% of the
  surface lies;
- the drainage network study's analysis runs unchanged on them: its geometry, rims, fills, chain, history and long
  profiles, on the design years' rain, with the open water's evaporation raised until the routing's budget closes,
  and at the ends of the rain range the climate lane carries;
- the design state the round's frame sets is then routed: Moscoviense, Orientale and Smythii-Marginis full at their
  rims, the near-side ocean, the south sea and Humboldtianum at the common level, and every other hollow and lake
  filling or drying by its own balance as geography's drainage has them. Each full sea's surplus enters the land at
  the first cell past the pass its spill path crosses and is routed down from there with the rest, filling the hollows
  on its way, so each sea's supply counts the overflows of the seas above it;
- the gravity routes by which a full sea's water could go elsewhere: each pass from its water into a neighbouring
  catchment whose water runs to another sea, with the cut its level would need there, and the route to each sea that
  needs the shallowest deepest cut, the water filling every hollow it passes and leaving each over its lowest pass or
  through a cut. These are the ground for dividing Moscoviense's water and for the channel to Humboldtianum.

The climate is the design years' on the turned Moon at 8 degrees, as the network study takes it; the round's run on
the final surface replaces it.
"""
from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
import time
from pathlib import Path

import numpy as np

from geography import crustal_loading as CL
from geography import drainage as dr
from geography import fetch_inputs
from geography import topography as tp
from research.studies.atlas_round import redrawn_surface as RS
from research.studies.drainage_network import network as NW
from research.studies.south_sea_ocean import inflows as I
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'redrawn_drainage.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/redrawn_drainage')
DESIGN_ROUTING = BULK / 'design_routing.npz'
CHECKPOINT = BULK / 'redrawn_drainage_so_far.json'
SCHEMA = 'terluna.research.atlas-round-redrawn-drainage/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-design-routing/1'
SURFACE = 'research/studies/atlas_round/results/redrawn_surface.json'
CODE = ('research/studies/atlas_round/redrawn_drainage.py', 'research/studies/atlas_round/redrawn_surface.py',
        'research/studies/drainage_network/network.py', 'research/studies/drainage_network/routing.py',
        'geography/drainage.py', 'geography/topography.py', 'geography/crustal_loading.py',
        'research/studies/south_sea_ocean/inflows.py')
UPPER = RS.UPPER
LOWER = ('near_side', 'south_pole_aitken', 'humboldtianum')
RIM_CHECK_M = 0.05                      # the rims found here agree with the first step's within this
WAYS_A_SEA = 8                          # the lowest single cuts to distinct river mouths listed for each full sea
WAYS_DEEPEST_M = 1500.0                 # and no deeper than this
SAME_GAP_KM = 25.0                      # passes closer than this to a lower one listed are taken as the same gap
WAYS_POOL = 40                          # the passes measured for each full sea before the shallowest cuts are listed
CUT_GRID_M = (0.0, 25.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 650.0, 800.0, 1000.0, 1250.0, 1500.0, 2000.0, 3000.0)
CUT_PRECISION_M = 2.0                   # the bisection on a route's deepest cut stops within this
LABELS_CAP = 3_000_000                  # the least-sum route search gives up beyond this many labels
CUTS_LISTED = 12                        # the deepest cuts of a route listed
CUT_COUNTED_M = 1.0                     # a cut shallower than this counts as none on a route


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def distance_km(a, b, cols, lat):
    """Great-circle distance between two cells of the 16 px/deg grid."""
    (ra, ca), (rb, cb) = divmod(int(a), cols), divmod(int(b), cols)
    p1, p2 = math.radians(lat[ra]), math.radians(lat[rb])
    dl = math.radians((cb - ca) / dr.PPD)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return MOON_RADIUS / 1000.0 * math.acos(max(-1.0, min(1.0, c)))


# ---------------------------------------------------------------- the redrawn floor

def load_surface(log=print):
    """LOLA's heights above the dry Moon's geoid with the settled case's U - N added (float32, 16 px/deg), the redrawn
    common level that holds 28% of the surface below it, the grid, and the first step's record."""
    started = time.time()
    if fetch_inputs.restore()['status'] != 'PASS':
        raise SystemExit('BLOCKED: topography inputs missing; run python -m geography.fetch_inputs --download')
    surface = json.loads((ROOT / SURFACE).read_text())
    case = surface['cases'][RS.DESIGN]
    grids = Path(case['grids']['path'])
    if sha256(grids) != case['grids']['sha256']:
        raise SystemExit(f'{grids} does not match the grids {SURFACE} records')
    atlas = json.loads((ROOT / NW.ATLAS).read_text())
    height, _, g = tp.height_above_geoid(dr.PPD, atlas['grid']['geoid_degree'])
    with np.load(grids, allow_pickle=False) as d:
        meta = json.loads(str(d['metadata']))
        dev = d['U'].astype(float) - d['N'].astype(float)
    eff = height + CL.upsample(dev, CL.TRANSFORM_FACTOR)
    del height, dev
    row_area = g.cell_area_m2[:, 0].copy()
    level = CL.level_for_area(eff, row_area, CL.SHARE, meta['level_m'])
    if abs(level - meta['level_m']) > RIM_CHECK_M:
        raise SystemExit(f'the common level {level:.3f} m differs from the first step\'s {meta["level_m"]:.3f} m')
    h = eff.astype(np.float32)
    del eff
    rows, cols = h.shape
    lat = 90.0 - (np.arange(rows) + 0.5) / dr.PPD
    cell_km = MOON_RADIUS / 1000.0 * np.radians(1.0 / dr.PPD)
    log(f'redrawn floor: common level {level - atlas["sea_level_m"]:+.2f} m against the atlas\'s, {time.time() - started:.0f} s')
    return dict(h=h, level0=level, lat=lat, cell_km=cell_km, row_area=row_area, atlas=atlas, surface=surface,
                case=case, grids=grids)


def catchment_graph(h, lat, cell_km):
    """NW.catchment_graph with each catchment's pit kept: the label of every cell's catchment, the pits (cell of each
    node, in node order) and the lowest pass between each pair of catchments that touch."""
    rows, cols = h.shape
    idx2 = np.arange(rows * cols, dtype=np.int32).reshape(rows, cols)
    receiver = dr.d8_receivers(h, np.zeros(h.shape, bool), lat, cell_km)
    sinks = np.flatnonzero(receiver == np.arange(receiver.size)).astype(np.int64)
    node_of = np.full(receiver.size, -1, np.int32)
    node_of[sinks] = np.arange(sinks.size, dtype=np.int32)
    label = node_of[dr.terminal(receiver)]
    del receiver, node_of
    edges = dr.passes(label.reshape(rows, cols), h.astype(np.float32), idx2)
    return label, sinks, edges


def geometry(surf, runoff, lake_net, log=print):
    """The drainage network study's geometry (network.geometry) on the redrawn heights and common level, with the
    catchment graph kept for the routes out of the full seas."""
    started = time.time()
    h, level0, lat, cell_km = surf['h'], surf['level0'], surf['lat'], surf['cell_km']
    rows, cols = h.shape
    flat = h.ravel()
    lon = (np.arange(cols) + 0.5) / dr.PPD
    body = I.label_waters(h < level0)
    named, _ = I.name_seas(body, surf['row_area'], lat, ((lon + 180.0) % 360.0) - 180.0, surf['atlas']['bodies'])
    seas = list(I.SEAS.values())
    label, sinks, edges = catchment_graph(h, lat, cell_km)
    n_nodes = sinks.size
    seeds, joined = {}, {}
    for k in seas:
        cells = np.flatnonzero(body.ravel() == named[k]['label'])
        seeds[k] = int(cells[np.argmin(flat[cells])])
        joined[k] = NW.joining_level(flat, label, n_nodes, edges, seeds[k])
    del body
    log(f'{n_nodes:,} catchments and the seas\' joining levels, {time.time() - started:.0f} s')
    targets, spill_paths = NW.spill_targets(h, lat, cell_km, runoff, lake_net, seas, joined, seeds, level0)
    log(f'spill targets {targets}, {time.time() - started:.0f} s')
    passes = {(a, b): float(joined[a][seeds[b]]) for a in seas for b in seas if a != b}
    return dict(level0=level0, h=h, lat=lat, cell_km=cell_km, seas=seas, joined=joined, seeds=seeds, targets=targets,
                spill_paths=spill_paths, passes=passes, catchments=n_nodes, graph=(label, sinks, edges))


def rims_against_the_first_step(geo, case, upper=UPPER):
    """Each full sea's rim here and in the first step's design case, both above the redrawn common level."""
    seas, passes, level0 = geo['seas'], geo['passes'], geo['level0']
    out = {}
    for k in upper:
        here = min(passes[(k, j)] for j in seas if j != k) - level0
        there = case['seas'][k]['surface_above_the_common_level_m'] + NW.BELOW_RIM_M
        out[k] = dict(here_m=round(here, 3), first_step_m=round(there, 3), agree=abs(here - there) <= RIM_CHECK_M)
    return out


# ---------------------------------------------------------------- the design state

def design_state(geo, runoff, lake_net, upper=UPPER, lower=LOWER, rounds=4, log=print):
    """The full seas at their rims and the others at the common level, each full sea's surplus put in at the first cell
    past the pass its spill path crosses and routed down with the rivers. Returns the routing with the surpluses in,
    each sea's level, area, loss, own rivers and supply, and each surplus's entry cell and the path it runs."""
    started = time.time()
    h, level0, lat, cell_km = geo['h'], geo['level0'], geo['lat'], geo['cell_km']
    seas, joined, passes = geo['seas'], geo['joined'], geo['passes']
    rows, cols = h.shape
    area = np.repeat((cell_km * 1e3) ** 2 * np.maximum(np.cos(np.radians(lat)), 1e-6), cols)
    rims = {k: min(passes[(k, j)] for j in seas if j != k) for k in upper}
    levels = {**{k: rims[k] - NW.BELOW_RIM_M for k in upper}, **{k: level0 for k in lower}}
    names = list(upper) + list(lower)
    masks = {k: joined[k] < levels[k] for k in names}
    loss_field = -np.asarray(lake_net, np.float64).ravel() * dr.MM_DAY * area
    loss = {k: float(loss_field[masks[k]].sum()) for k in names}
    waters = [(masks[k], levels[k]) for k in names]
    sea_of = np.full(h.size, -1, np.int8)
    for i, k in enumerate(names):
        sea_of[masks[k]] = i
    r, got = NW.route_at(h, lat, cell_km, runoff, lake_net, waters)
    own = dict(zip(names, got))
    lake_full = ~r['outlet'] & (h.ravel() < r['level'][r['label']])
    entry = {}
    for k in upper:
        path, _ = geo['spill_paths'][k]
        leave = np.flatnonzero(joined[k][path] >= rims[k] - NW.BELOW_RIM_M)
        entry[k] = int(path[leave[0]])
    supply, put_in = dict(own), {}
    for _ in range(rounds):
        add_runoff = np.asarray(runoff, np.float64).copy()
        add_lake = np.asarray(lake_net, np.float64).copy()
        put_in, direct = {}, {}
        for k in upper:
            surplus = supply[k] - loss[k]
            if surplus <= 0.0:
                continue
            put_in[k] = surplus
            c = entry[k]
            if sea_of[c] >= 0:                           # the pass opens straight onto another sea
                direct[names[sea_of[c]]] = direct.get(names[sea_of[c]], 0.0) + surplus
                continue
            i, j = divmod(c, cols)
            target = add_lake if lake_full[c] else add_runoff
            target[i, j] += surplus / (dr.MM_DAY * area[c])
        del r
        r, got = NW.route_at(h, lat, cell_km, add_runoff, add_lake, waters)
        del add_runoff, add_lake
        new = dict(zip(names, got))
        for k, q in direct.items():
            new[k] += q
        moved = max(abs(new[k] - supply[k]) for k in upper)
        supply = new
        if moved < 1.0:
            break
    else:
        raise RuntimeError('the full seas\' overflows did not settle')
    paths = {}
    for k in upper:
        path = walk(r, entry[k])
        paths[k] = (names[sea_of[path[-1]]] if sea_of[path[-1]] >= 0 else None, path)
    log(f'design state: supplies {dict((k, round(v)) for k, v in supply.items())}, {time.time() - started:.0f} s')
    return dict(r=r, names=names, levels=levels, masks=masks, area={k: float(area[masks[k]].sum()) for k in names},
                loss=loss, own=own, supply=supply, put_in=put_in, entry=entry, paths=paths, sea_of=sea_of)


def design_record(state, geo, rotation):
    """The design state's numbers: each sea's level against the common level, area, loss to the air beyond the rain on
    it, its own rivers, what the full seas above send it, and its balance; each full sea's overflow with its long
    profile; and what each sea at the common level needs from the full seas to stay there."""
    r, level0, h, lat, cell_km = state['r'], geo['level0'], geo['h'], geo['lat'], geo['cell_km']
    cols = h.shape[1]
    surface = np.where(r['lake'], np.nan_to_num(r['lake_level'], nan=0.0), h.ravel())
    out = {}
    for k in state['names']:
        e = dict(level_m=round(state['levels'][k] - level0, 2), area_km2=round(state['area'][k] / 1e6),
                 loss_m3s=round(state['loss'][k]), own_rivers_m3s=round(state['own'][k]),
                 from_the_full_seas_m3s=round(state['supply'][k] - state['own'][k]),
                 supply_m3s=round(state['supply'][k]), balance_m3s=round(state['supply'][k] - state['loss'][k]))
        if k in state['entry']:
            target, path = state['paths'][k]
            profile = NW.long_profile(path, surface, r['discharge'], lat, cols, cell_km, level0, r['lake'])
            e.update(passes_on_m3s=round(state['put_in'].get(k, 0.0)), leaves=NW.place(state['entry'][k], cols, rotation),
                     reaches=target, overflow=profile)
        else:
            e['needs_from_the_full_seas_to_stay_m3s'] = round(max(state['loss'][k] - state['own'][k], 0.0))
        out[k] = e
    return out


# ---------------------------------------------------------------- the ways a full sea's water could leave

def node_levels(geo, state):
    """The water level of each catchment of the graph in the design state: a sea's level where the catchment's pit lies
    under a sea, otherwise the level its hollow fills to before it spills (the routing's flood)."""
    label, sinks, _ = geo['graph']
    r = state['r']
    level = np.empty(sinks.size)
    in_sea = state['sea_of'][sinks]
    for i, k in enumerate(state['names']):
        level[in_sea == i] = state['levels'][k]
    land = in_sea < 0
    level[land] = r['level'][r['label'][sinks[land]]]
    return level, in_sea


def walk(r, cell):
    """The cells from one cell down the routing to the sea it reaches."""
    path, seen = [int(cell)], {int(cell)}
    while not r['outlet'][cell]:
        cell = int(r['receiver'][cell])
        if cell in seen:
            break
        seen.add(cell)
        path.append(cell)
    return np.array(path, dtype=np.int64)


def single_cuts(geo, state, rotation, upper=UPPER, count=WAYS_A_SEA, deepest=WAYS_DEEPEST_M, same_gap_km=SAME_GAP_KM,
                pool=WAYS_POOL):
    """For each full sea, the passes on the divide round the land that drains to it into land whose water runs to
    another sea, one for each river mouth beyond: the deepest cut a channel at the sea's level needs from its shore over
    the pass and down the far side until the water can run on by itself (where the far side's hollows spill below the
    sea's level), the cut's length, and the river the water then joins. The lowest way to each sea is listed, then the
    lowest of the rest."""
    r, names, h, lat, cell_km, level0 = state['r'], state['names'], geo['h'], geo['lat'], geo['cell_km'], geo['level0']
    rows, cols = h.shape
    hf = h.ravel()
    mouth = dr.terminal(r['receiver'])
    goes = state['sea_of'][mouth]                         # the sea each cell's water reaches
    free_at = r['level'][r['label']]                      # the level water must reach in a cell's hollow to run on
    idx2 = np.arange(rows * cols, dtype=np.int64).reshape(rows, cols)
    surface = np.where(r['lake'], np.nan_to_num(r['lake_level'], nan=0.0), hf)
    profile = lambda path: NW.long_profile(path, surface, r['discharge'], lat, cols, cell_km, level0, r['lake'])
    out = {}
    for k in upper:
        i_k = names.index(k)
        level_k = state['levels'][k]
        shed = (goes == i_k).reshape(rows, cols)
        c_in, c_out = [], []
        for dr_, dc in ((0, 1), (1, 0), (1, 1), (1, -1)):
            nb = dr.shifted(idx2, dr_, dc, -1)
            cross = (nb >= 0) & (shed != dr.shifted(shed, dr_, dc, False))
            a_, b_ = idx2[cross], nb[cross]
            a_in = shed[cross]
            c_in.append(np.where(a_in, a_, b_))
            c_out.append(np.where(a_in, b_, a_))
        c_in, c_out = np.concatenate(c_in), np.concatenate(c_out)
        height = np.maximum(hf[c_in], hf[c_out]).astype(float)
        m = mouth[c_out]
        reach = goes[c_out]
        keep = (reach >= 0) & (reach != i_k) & (height - level_k <= deepest)
        c_in, c_out, height, m, reach = c_in[keep], c_out[keep], height[keep], m[keep], reach[keep]
        order = np.lexsort((height, m))
        first = np.ones(order.size, bool)
        first[1:] = m[order][1:] != m[order][:-1]
        best = order[first]
        best = best[np.argsort(height[best], kind='stable')].tolist()
        chosen, seen_seas = [], set()
        for j in best:                                   # the lowest pass to each sea, then the lowest of the rest
            if reach[j] not in seen_seas:
                chosen.append(j)
                seen_seas.add(int(reach[j]))
        for j in best:                                   # a pass near one already taken is the same gap
            if len(chosen) >= max(pool, len(seen_seas)):
                break
            if j not in chosen and all(distance_km(c_out[j], c_out[i], cols, lat) > same_gap_km for i in chosen):
                chosen.append(j)
        ways = []
        for j in chosen:
            inside, beyond = walk(r, c_in[j]), walk(r, c_out[j])
            free = np.flatnonzero(free_at[beyond] <= level_k)
            cut_beyond = beyond[:free[0] + 1] if free.size else beyond
            top = max(height[j], float(hf[inside].max()), float(hf[cut_beyond].max()))
            p_out = profile(beyond)
            ways.append(dict(cut_m=round(top - level_k, 1),
                             pass_above_the_common_level_m=round(float(height[j] - level0), 1),
                             at=NW.place(int(c_out[j]), cols, rotation), reaches=names[reach[j]],
                             cut_length_km=round(profile(inside)['length_km'] + profile(cut_beyond)['length_km']),
                             from_the_shore_km=profile(inside)['length_km'],
                             mouth=NW.place(int(m[j]), cols, rotation), river_length_km=p_out['length_km'],
                             lakes_passed=p_out['lakes_passed'], steep_km=p_out['steep_km'],
                             discharge_now_at_the_pass_and_mouth_m3s=p_out['discharge_m3s']))
        ways.sort(key=lambda x: x['cut_m'])
        listed, seas_listed = [], set()
        for x in ways:                                   # the shallowest cut to each sea, then the shallowest others
            if x['reaches'] not in seas_listed:
                listed.append(x)
                seas_listed.add(x['reaches'])
        listed += [x for x in ways if x not in listed][:max(0, count - len(listed))]
        out[k] = sorted(listed, key=lambda x: x['cut_m'])
    return out


def least_cut_routes(geo, state, levels, in_sea, source, rotation, grid=CUT_GRID_M, precision=CUT_PRECISION_M):
    """From one full sea's water, the gravity route to each other sea whose deepest cut is the shallowest. A channel
    carries the water at its level; where it enters a catchment whose hollow spills lower, the water falls into that
    hollow and runs on at the hollow's level, and where it crosses a pass above its level the pass is cut down to it.
    The shallowest deepest cut to each sea is found by bisection on the cut allowed, each trial a search that keeps the
    water as high as it can (which only ever helps later passes). Returns, for each sea, the deepest cut, the cuts
    deeper than CUT_COUNTED_M with their places, the seas passed on the way and the route's length through its
    passes; None where no cut within the grid reaches it."""
    label, sinks, (u, v, w, cu, cv) = geo['graph']
    names, h, lat, level0 = state['names'], geo['h'], geo['lat'], geo['level0']
    cols = h.shape[1]
    n = sinks.size
    src = np.concatenate([u, v]).astype(np.int64)
    dst = np.concatenate([v, u]).astype(np.int64)
    wt = np.concatenate([w, w]).astype(float)
    c_from = np.concatenate([cu, cv]).astype(np.int64)
    c_to = np.concatenate([cv, cu]).astype(np.int64)
    order = np.argsort(src, kind='stable')
    src, dst, wt, c_from, c_to = (x[order] for x in (src, dst, wt, c_from, c_to))
    start = np.searchsorted(src, np.arange(n + 1)).tolist()
    dst_l, wt_l, lev_l = dst.tolist(), wt.tolist(), levels.tolist()
    i_s = names.index(source)
    starts = np.flatnonzero(in_sea == i_s).tolist()
    level_s = float(levels[starts[0]])
    targets = {k: np.flatnonzero(in_sea == i).tolist() for i, k in enumerate(names) if k != source}

    def spread(allowed):
        best, prev = [-math.inf] * n, [-1] * n
        for x in starts:
            best[x] = level_s
        heap = [(-level_s, x) for x in starts]
        heapq.heapify(heap)
        while heap:
            neg, x = heapq.heappop(heap)
            lev = -neg
            if lev < best[x]:
                continue
            limit = lev + allowed
            for e in range(start[x], start[x + 1]):
                if wt_l[e] > limit:
                    continue
                y = dst_l[e]
                ly = lev if lev < lev_l[y] else lev_l[y]
                if ly > best[y]:
                    best[y], prev[y] = ly, e
                    heapq.heappush(heap, (-ly, y))
        return best, prev

    reached = lambda best, k: any(best[x] > -math.inf for x in targets[k])
    bracket, below = {}, -1e-3
    for allowed in grid:
        best, _ = spread(allowed)
        for k in targets:
            if k not in bracket and reached(best, k):
                bracket[k] = (below, allowed)
        below = allowed
        if len(bracket) == len(targets):
            break
    def cleanest(allowed, k, cap=LABELS_CAP):
        """Among routes whose cuts all lie within `allowed`, one whose cuts add up to the least: a search on the cuts'
        sum whose labels carry the water's level, a catchment keeping each label no other there beats on both."""
        store = [(x, 0.0, level_s, -1, -1) for x in starts]          # node, cuts so far, water level, label before, edge
        labels = {x: [(0.0, level_s, i)] for i, x in enumerate(starts)}
        heap = [(0.0, -level_s, i) for i in range(len(starts))]
        heapq.heapify(heap)
        dead, goal = set(), set(targets[k])
        while heap:
            total, _, i = heapq.heappop(heap)
            if i in dead:
                continue
            x, _, lev, _, _ = store[i]
            if x in goal:
                chain = []
                while store[i][3] >= 0:
                    chain.append((store[i][4], store[store[i][3]][2]))
                    i = store[i][3]
                return chain[::-1]
            limit = lev + allowed
            for e in range(start[x], start[x + 1]):
                if wt_l[e] > limit:
                    continue
                y = dst_l[e]
                c = wt_l[e] - lev
                nt = total + (c if c > 0.0 else 0.0)
                ly = lev if lev < lev_l[y] else lev_l[y]
                here = labels.get(y, [])
                if any(t2 <= nt and l2 >= ly for t2, l2, _ in here):
                    continue
                kept = []
                for t2, l2, j in here:
                    if nt <= t2 and ly >= l2:
                        dead.add(j)
                    else:
                        kept.append((t2, l2, j))
                kept.append((nt, ly, len(store)))
                labels[y] = kept
                store.append((y, nt, ly, i, e))
                heapq.heappush(heap, (nt, -ly, len(store) - 1))
                if len(store) > cap:
                    return None
        return None

    out = {k: None for k in targets}
    for k, (lo, hi) in bracket.items():
        while hi - lo > precision:
            mid = 0.5 * (lo + hi)
            if reached(spread(mid)[0], k):
                hi = mid
            else:
                lo = mid
        chain = cleanest(hi, k)
        tidy = chain is not None
        if not tidy:                                         # the search grew too large: the highest-water route
            best, prev = spread(hi)
            x, edges_back = max(targets[k], key=lambda x: best[x]), []
            while prev[x] >= 0 and in_sea[x] != i_s:
                edges_back.append(prev[x])
                x = int(src[prev[x]])
            chain = [(e, best[int(src[e])]) for e in edges_back[::-1]]
        depth = [max(0.0, wt_l[e] - lev) for e, lev in chain]
        cuts = [dict(depth_m=round(d, 1), pass_above_the_common_level_m=round(wt_l[e] - level0, 1),
                     water_above_the_common_level_m=round(lev - level0, 1),
                     at=NW.place(int(c_to[e]), cols, rotation)) for (e, lev), d in zip(chain, depth) if d >= CUT_COUNTED_M]
        passed = sorted({names[in_sea[int(dst[e])]] for e, _ in chain[:-1] if in_sea[int(dst[e])] >= 0} - {source, k})
        edges_only = [e for e, _ in chain]
        length = sum(distance_km(c_to[a], c_from[b], cols, lat) for a, b in zip(edges_only[:-1], edges_only[1:])) + \
            sum(distance_km(c_from[e], c_to[e], cols, lat) for e in edges_only)
        out[k] = dict(deepest_cut_m=round(max(depth, default=0.0), 1), cuts_total_m=round(sum(depth), 1),
                      cuts_over_1_m=len(cuts), cuts=sorted(cuts, key=lambda c: -c['depth_m'])[:CUTS_LISTED],
                      cuts_least_in_sum=tidy, seas_passed=passed, passes=len(chain),
                      length_through_the_passes_km=round(length))
    return out


def lake_at(geo, sea, supply_m3s, loss_field, area, level0, cap_m=1.0):
    """The level and area at which a sea at the common level or below settles on a given supply: its loss to the air
    meets the supply (network.Table on its joining levels)."""
    table = NW.Table(geo['joined'][sea], area, loss_field, geo['h'].ravel(), level0 + cap_m)
    held = table.level_for(supply_m3s, level0 + 1e-6)
    lev = level0 if held is None else held
    a, _, v = table.at(lev)
    return dict(supply_m3s=round(supply_m3s), level_m=round(lev - level0, 1), area_km2=round(a / 1e6),
                volume_km3=round(v / 1e9), held_at_the_common_level=held is None)


# ---------------------------------------------------------------- the run

def design_product(state, geo, path=DESIGN_ROUTING):
    """The design state's routing on the routing's own budget, for the satoyama map and the globe: every cell's
    receiver, mean discharge and drained area, the lakes with their levels and each cell's sea."""
    r = state['r']
    lake = r['lake'] & ~r['outlet']
    cells = np.flatnonzero(lake).astype(np.int32)
    meta = dict(schema=PRODUCT_SCHEMA, common_level_m=geo['level0'], seas=state['names'],
                levels_m={k: state['levels'][k] for k in state['names']}, grid=dict(pixels_per_degree=dr.PPD,
                rows=geo['h'].shape[0], cols=geo['h'].shape[1]),
                reading_rule=('Cells index rows from 90 N and columns from 0 E, 16 per degree: cell = row * 5760 + column, '
                              'on today\'s Moon. receiver is the cell each cell drains to (a full hollow\'s pit to the '
                              'cell past its pass, a sea\'s cells to themselves); discharge_m3s the mean flow through '
                              'the cell with the full seas\' overflows; drained_km2 the land that drains through it; '
                              'sea the index into seas of the sea over the cell, -1 elsewhere; lake_cell the cells under '
                              'lakes off the seas, lake_level_m their surface above the dry Moon\'s geoid as redrawn.'))
    BULK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, metadata=json.dumps(meta), receiver=r['receiver'].astype(np.int32),
                        discharge_m3s=r['discharge'].astype(np.float32),
                        drained_km2=(r['drained_m2'] / 1e6).astype(np.float32), sea=state['sea_of'],
                        lake_cell=cells, lake_level_m=r['lake_level'][cells].astype(np.float32),
                        lake_spills=~r['closed'][r['label'][cells]])
    return dict(path=str(path), sha256=sha256(path))


def design_piece(geo, rain, evaporation, rotation, area, full, log=print, write_product=True):
    """The design state on one rain and evaporation, with the routes out of the full seas, the rivers' long profiles and
    the product when `full` and `write_product`."""
    _, _, _, _, _, runoff, lake_net, _, _ = NW.load(rain, evaporation)
    state = design_state(geo, runoff, lake_net, log=log)
    record = dict(evaporation=evaporation, seas=design_record(state, geo, rotation))
    loss_field = -lake_net.ravel().astype(np.float64) * dr.MM_DAY * area
    h_need = record['seas']['humboldtianum']['needs_from_the_full_seas_to_stay_m3s']
    sm_spare = record['seas']['smythii_marginis']['passes_on_m3s']
    own_h = state['own']['humboldtianum']
    record['humboldtianum_on_a_channel'] = [
        lake_at(geo, 'humboldtianum', own_h + q, loss_field, area, geo['level0'])
        for q in sorted({0.0, 1000.0, 2500.0, 5000.0, float(min(h_need, sm_spare))})]
    product = None
    if full:
        started = time.time()
        levels, in_sea = node_levels(geo, state)
        record['single_cuts'] = single_cuts(geo, state, rotation)
        log(f'single cuts, {time.time() - started:.0f} s')
        record['least_cut_routes'] = {k: least_cut_routes(geo, state, levels, in_sea, k, rotation) for k in UPPER}
        log(f'least-cut routes, {time.time() - started:.0f} s')
        record['rivers'] = NW.river_profiles(state['r'], state['masks'], geo['h'], geo['lat'], geo['cell_km'],
                                             geo['level0'], rotation)
        if write_product:
            product = design_product(state, geo)
    return record, product


def build(rains=NW.RAINS, log=print, resume=False):
    """The network study's cases and the design state on the redrawn floor. Each piece is kept on the research drive
    as it is finished (CHECKPOINT), and `resume` takes the finished pieces from there, the geometry being found again."""
    started = time.time()
    so_far = json.loads(CHECKPOINT.read_text()) if resume and CHECKPOINT.exists() else {}

    def keep(key, value):
        so_far[key] = value
        BULK.mkdir(parents=True, exist_ok=True)
        CHECKPOINT.write_text(json.dumps(so_far, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
        log(f'kept {key}, {time.time() - started:.0f} s')

    surf = load_surface(log)
    _, _, _, _, _, runoff, lake_net, _, _ = NW.load(1.0)
    geo = geometry(surf, runoff, lake_net, log)
    del runoff, lake_net
    rims = rims_against_the_first_step(geo, surf['case'])
    log(f'rims against the first step: {rims}')
    if not all(v['agree'] for v in rims.values()):
        raise SystemExit('the rims do not reproduce the first step\'s design case')
    rotation = np.array(json.loads((ROOT / NW.FRAME).read_text())['frame']['rotation_today_to_turned'])
    rows, cols = geo['h'].shape
    area = np.repeat((geo['cell_km'] * 1e3) ** 2 * np.cos(np.radians(geo['lat'])), cols)
    if 'design_routing_budget' not in so_far:
        # the routing is the same each time; a product already written and unchanged is kept, so what reads it stays bound
        kept_product = 'design_routing' in so_far and DESIGN_ROUTING.exists() and \
            sha256(DESIGN_ROUTING) == so_far['design_routing']['sha256']
        record, product = design_piece(geo, rains[0], 1.0, rotation, area, True, log, write_product=not kept_product)
        if product is not None:
            keep('design_routing', product)
        keep('design_routing_budget', record)
    if 'case_0' not in so_far:
        keep('case_0', NW.network(rains[0], geo, log))
    factor = so_far['case_0']['budget']['open_water_evaporation_factor_closing_it']
    if 'design_closed_budget' not in so_far:
        keep('design_closed_budget', design_piece(geo, rains[0], factor, rotation, area, False, log)[0])
    for i, (rain, evaporation) in enumerate([(rains[0], factor)] + [(rain, 1.0) for rain in rains[1:]], start=1):
        if f'case_{i}' not in so_far:
            keep(f'case_{i}', NW.network(rain, geo, log, evaporation=evaporation))
    cases = [so_far[f'case_{i}'] for i in range(len(rains) + 1)]
    design = dict(routing_budget=so_far['design_routing_budget'], closed_budget=so_far['design_closed_budget'])
    product = so_far['design_routing']
    level_atlas = surf['atlas']['sea_level_m']
    inputs = [ROOT / SURFACE, ROOT / NW.ATLAS, ROOT / NW.DRAINAGE, ROOT / NW.FRAME]
    rel = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    return dict(
        schema=SCHEMA,
        producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in inputs},
        ignored_inputs={rel(ROOT / NW.CLIMATOLOGY): sha256(ROOT / NW.CLIMATOLOGY), str(surf['grids']): surf['case']['grids']['sha256'],
                        **{rel(fetch_inputs.path(n)): sha256(fetch_inputs.path(n)) for n in ('ldem_16.img', 'jggrx_0420a_sha.tab')}},
        climatology=dict(run=NW.RUN, years=NW.YEARS),
        evidence=('A calculation on products: the drainage network study\'s routing and analysis of the design years\' '
                  'rain over LOLA\'s 16 pixel/degree heights redrawn at the crust\'s settled state (the first step\'s '
                  'design case), with the design state routed with the three upper seas full at their rims and their '
                  'overflows run down to the seas below. The climate is the present design climate on the turned Moon '
                  'at 8 degrees, as the atlas\'s waters made it, and the routing\'s budget differs from the climate '
                  'model\'s as the network study found. Routes and cuts are found between 1.9 km cells. Not validated '
                  'against any observation.'),
        reading_rule=('Levels in m above the redrawn common level unless named otherwise; flows in m3/s; areas in km2 and '
                      'volumes in km3. cases are the drainage network study\'s records (network.json\'s form) on the '
                      'redrawn floor: the design years\' rain, the same with the open water\'s evaporation raised to '
                      'close the routing\'s budget, and the rain range\'s ends. design is the round\'s state on the '
                      'design years\' rain, on the routing\'s own budget and on the closed one: the three upper seas '
                      'full and spilling, the other three at the common level; balance is a sea\'s supply less its '
                      'loss to the air beyond the rain on it, and a sea at the common level with a negative balance '
                      'would fall unless the full seas send it needs_from_the_full_seas_to_stay. single_cuts and '
                      'least_cut_routes are gravity routes for a full sea\'s water, a cut being how far a pass stands '
                      'above the water that would run through it.'),
        common_level=dict(above_the_dry_geoid_m=round(geo['level0'], 2), against_the_atlas_m=round(geo['level0'] - level_atlas, 2)),
        rims_against_the_first_step=rims,
        catchments=geo['catchments'], spill_targets=geo['targets'],
        seeds_today={k: NW.place(c, cols)['today'] for k, c in geo['seeds'].items()},
        cases=cases, design=design, design_routing=product, seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--rain', type=float, nargs='*', default=list(NW.RAINS), help='factors on the design years\' rain')
    p.add_argument('--out', type=Path, default=OUT)
    p.add_argument('--resume', action='store_true', help='take the pieces already finished from the research drive')
    args = p.parse_args(argv)
    record = build(tuple(args.rain), log=lambda *a: print(*a, flush=True), resume=args.resume)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    print(f'wrote {args.out} in {record["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
