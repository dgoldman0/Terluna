"""The design surface's sustained rivers: the design climate's fifteen settled years routed one at a time on the design
routing's network, without storage and with the lakes and the upper seas carrying water from one year to the next; the
rivers' losses in the drylands on the design climate; and the major rivers classed as sustained or intermittent.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.sustained_rivers   # about 10 minutes

The first run also builds the design routing's structure, about 12 minutes, and keeps it on the research drive; --resume
takes the years already routed from the drive's checkpoint.

The design routing (design_routing.py) draws every river at the long-term mean flow of the design climate's closed
budget. This module asks which rivers keep flowing through the climate's dry years, in four parts:

- the years without storage. Each settled year of the design climate (the run on the design surface at sunlight 0.936,
  years 45-59) is read from the run's year sums (moisture_map.cached_year in drylands.py's cache), moved by the
  estimate's step as moisture_final.adjusted moves them, carried onto today's ground with the settled frame's rotation
  as final_drainage.carry carries the mean, and loaded as final_drainage.load loads it, with the closed budget's factor
  on the open water's evaporation. It is routed on the design routing's structure held fixed: the floor, the seas at
  their design levels, the lakes as the design fills them and the works where the design puts them. Each water passes
  on that year's own surplus: a lake standing full or an upper sea passes on what the year brings it beyond its loss,
  a hollow the design leaves below its pass passes nothing, Moscoviense's fixed crests divide its overflow in the
  design's fixed shares, and the near-side sea's main part sends over its neck the design's spill moved by the year's
  departure in its balance. Each year is also routed as a steady state of its own, every hollow filling and spilling
  as the routing has it, to show what the hollows held below their passes would add;
- the years with storage. Every lake and the three upper seas carry their volume from one year to the next, starting
  from the design state, in STEPS steps a year. A lake is a body of the routing's catchments standing at one level,
  whose outlet holds its volume. Above its spill level a water passes on what critical flow passes over an outlet as
  wide as Earth's alluvial channels at the water's design flow (river_losses.WIDTH, widened for lunar gravity, at both
  ends of the range); below it, it passes nothing and evaporates from the area its level covers. The main part's spill
  is critical flow over the neck's throat as near_side_neck.json measures it;
- the losses. river_losses.py's estimate (the channels, the plants along them and the floodplains, in the drylands) is
  run again on the design routing and the design climate, and each loss is taken from the river where it falls and
  carried down to the seas;
- the classes. Every reach of the rivers carrying at least THRESHOLD_M3S on the design routing is sustained when it
  carries that much in every year with the wide outlets' storage, and intermittent otherwise, with its mean, lowest-year
  and second-lowest-year flows in each way of routing the years, the years it carries that much and its flow after the
  losses.

The product, atlas-round/sustained_rivers/sustained_rivers_16ppd.npz on the research drive, holds each river cell's
flows by year and its sustained discharge and intermittent flag; results/sustained_rivers.json records it.
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
from research.studies.atlas_round import design_routing as DRt
from research.studies.atlas_round import drylands as DL
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import moisture_final as MF
from research.studies.atlas_round import near_side_neck as N
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import river_losses as RL
from research.studies.drainage_network import network as NW
from research.studies.drainage_network.routing import flood_from
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, MOON_SURFACE_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'sustained_rivers.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/sustained_rivers')
PRODUCT = BULK / 'sustained_rivers_16ppd.npz'
REACHES = BULK / 'sustained_rivers_reaches.json'
WORK = BULK / 'work'                            # pieces kept as they finish, so a rerun can resume
STRUCTURE = WORK / 'design_structure.npz'
TS_PRODUCTS = ROOT / 'climate' / 'gcm' / 'products'      # turned_surface.PRODUCTS, where the GCM surfaces are
SCHEMA = 'terluna.research.atlas-round-sustained-rivers/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-sustained-rivers-grid/1'
CODE = ('research/studies/atlas_round/sustained_rivers.py', 'research/studies/atlas_round/design_routing.py',
        'research/studies/atlas_round/final_drainage.py', 'research/studies/atlas_round/governed_drainage.py',
        'research/studies/atlas_round/redrawn_drainage.py', 'research/studies/atlas_round/near_side_neck.py',
        'research/studies/atlas_round/moisture_final.py', 'research/studies/atlas_round/river_losses.py',
        'research/studies/drainage_network/network.py', 'research/studies/drainage_network/routing.py',
        'geography/drainage.py', 'climate/gcm/turned_climatology.py', 'climate/gcm/pattern_scaling.py',
        'biosphere/ecology/moisture_map.py')
STRUCTURE_CODE = CODE[1:7] + CODE[8:11]          # what the design structure is built by
RUN, YEARS = FD.RUN, tuple(FD.YEARS)
BUDGET = DRt.BUDGET
WATERS = ('moscoviense', 'orientale', 'smythii_marginis', 'south_pole_aitken', 'humboldtianum', 'near_side_main',
          'near_side_south')                     # design_routing.route's waters, in its order
OPS = ('salt_sea_weir', 'orientale', 'smythii_marginis', 'humboldtianum_branch', 'smythii_marginis_branch', 'neck')
THRESHOLD_M3S = 100.0                # a major river: the falls' and the losses' threshold
LARGE_M3S = 1000.0                   # the large rivers, counted apart
TRACK_M3S = 1.0                      # the river cells kept year by year: the globe's channels
DRY_M3S = 1.0                        # a river carrying less than this in a year has stopped flowing
LAKE_BREAK_KM2 = 100.0               # a lake this large ends a reach; smaller ones lie on its course
STEPS = 24                           # steps a year of the storage
Q_FLOOR = 1.0                        # m3/s: the flow an outlet is sized at where the design passes on less
DEPTHS_M = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0,
            8192.0, 16384.0)          # the levels below a lake's spill level its area and loss are tabled at
G = MOON_SURFACE_GRAVITY
YEAR_S = JULIAN_YEAR_DAYS * JULIAN_DAY
MM_DAY = dr.MM_DAY


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def rel(path) -> str:
    return str(Path(path).relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path)


# ---------------------------------------------------------------- the design structure

def design_structure(log=print, refresh=False):
    """The design routing's fixed structure, built as design_routing.route builds it: the redrawn floor's heights, the
    water over each cell (an index into WATERS, -1 elsewhere) with each water's level, and the governed routing's
    put-ins and take-outs with the cell, the kind (OPS) and the sign of each, whether design_routing places it in a
    lake, and its flow in the design. Kept on the research drive under the design routing record's hash and its
    code's, since building it takes the governed routing's five minutes."""
    key = dict(record=sha256(DRt.OUT), code={f: sha256(ROOT / f) for f in STRUCTURE_CODE})
    if STRUCTURE.exists() and not refresh:
        with np.load(STRUCTURE, allow_pickle=False) as z:
            meta = json.loads(str(z['metadata']))
            if meta['key'] == key:
                log(f'the design structure from {STRUCTURE}')
                return dict(meta=meta, **{k: z[k] for k in z.files if k != 'metadata'})
    behind = DRt.stale()
    if behind:
        raise SystemExit(f'the design routing is stale ({", ".join(behind)})')
    started = time.time()
    fd = json.loads(FD.OUT.read_text())
    neck = json.loads(N.OUT.read_text())
    case = fd['cases'][BUDGET]
    frame_path, record = FD.frame_record()
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    surf = RD.load_surface(log)
    runoff, lake_net, _, _ = FD.load(case['evaporation'])
    geo = RD.geometry(surf, runoff, lake_net, log)
    level0, h = geo['level0'], geo['h']
    cols = h.shape[1]
    h_flat = h.ravel()
    state = RD.design_state(geo, runoff, lake_net, log=log)
    branches, _ = GD.branches_for(geo, state, rotation, log)
    gov = GD.govern(geo, state, runoff, lake_net, branches, log=log)
    if any(onto is not None for _, onto in branches.values()):
        raise SystemExit('a branch opens straight onto a sea, which this module does not carry')
    steady = neck['budgets'][BUDGET]['natural_neck']['steady']
    z_main, z_south = level0 + steady['main_level_m'], level0 + steady['south_level_m']
    st = N.setting(log)
    main_water = st['part']['main'] & (st['joined']['main'] < z_main)
    south_water = st['part']['south'] & (st['joined']['south'] < z_south)
    del st
    names, masks, levels = state['names'], state['masks'], state['levels']
    water = np.full(h_flat.size, -1, np.int8)
    level_of = []
    for k in names:
        if k == 'near_side':
            continue
        water[masks[k]] = WATERS.index(k)
        level_of.append(levels[k])
    water[main_water] = WATERS.index('near_side_main')
    water[south_water] = WATERS.index('near_side_south')
    level_of += [z_main, z_south]
    if [k for k in names if k != 'near_side'] + ['near_side_main', 'near_side_south'] != list(WATERS):
        raise SystemExit('design_routing.route\'s waters are not in WATERS\' order')
    r0 = state['r']
    lake_full0 = ~r0['outlet'] & (h_flat < r0['level'][r0['label']])
    spill = DRt.spill_cell(h_flat, main_water, cols, neck['throat']['centre_today'])
    entry, share, over = state['entry'], gov['share'], gov['over']
    ops = [(entry['moscoviense'], 'salt_sea_weir', 1.0, share[GD.SALT_SEA])]
    ops += [(entry[k], k, 1.0, over[k]) for k in GD.FULL_SEAS]
    for target, kind in ((GD.HELD_BRANCH, 'humboldtianum_branch'), (GD.REST_BRANCH, 'smythii_marginis_branch')):
        ops += [(c, kind, float(sign), sign * share[target]) for c, sign in branches[target][0]]
    ops.append((spill, 'neck', 1.0, float(steady['neck_flow_m3s'])))
    if DRt.governed_ops(state, gov, branches) + [(spill, steady['neck_flow_m3s'])] != [(c, q) for c, _, _, q in ops]:
        raise SystemExit('the put-ins and take-outs are not design_routing.route\'s')
    surplus = float(gov['surplus'])
    meta = dict(key=key, level0=float(level0), cell_km=float(geo['cell_km']), levels=[float(v) for v in level_of],
                rotation=rotation.tolist(), frame=dict(record=rel(frame_path), sha256=sha256(frame_path)),
                evaporation=float(case['evaporation']), surplus=surplus,
                fractions={k: float(share[k]) / surplus for k in (GD.SALT_SEA, GD.HELD_BRANCH, GD.REST_BRANCH)},
                over={k: float(v) for k, v in over.items()}, neck_flow=float(steady['neck_flow_m3s']),
                loss={k: float(v) for k, v in state['loss'].items()}, own={k: float(v) for k, v in state['own'].items()},
                seconds=round(time.time() - started))
    arrays = dict(h=h, lat=geo['lat'], water=water,
                  op_cell=np.array([c for c, _, _, _ in ops], np.int64),
                  op_kind=np.array([OPS.index(k) for _, k, _, _ in ops], np.int8),
                  op_sign=np.array([s for _, _, s, _ in ops], float),
                  op_lake=np.array([bool(lake_full0[c]) for c, _, _, _ in ops]),
                  op_design=np.array([q for _, _, _, q in ops], float))
    WORK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(STRUCTURE, metadata=json.dumps(meta), **arrays)
    log(f'the design structure built in {meta["seconds"]} s and kept in {STRUCTURE}')
    return dict(meta=meta, **arrays)


# ---------------------------------------------------------------- the routing on the fixed structure

class Network:
    """The design routing's network held fixed, and its routing on any year's water. The graph is routing.route_levels'
    on the design's floor and waters: each cell's steepest-descent neighbour, the catchments of the pits with the seas'
    outlets as roots of the flood, each catchment's spill level and the pass it spills through. replay() accumulates a
    year's runoff over it exactly as route_levels' loop does, the generations of that loop kept; cascade() does the same
    for the catchments alone, which is all the seas' balances need."""

    def __init__(self, h, water, levels, lat, cell_km, log=print):
        started = time.time()
        rows, cols = h.shape
        n = rows * cols
        self.rows, self.cols, self.n = rows, cols, n
        self.water = water
        self.n_waters = len(levels)
        outlet = water >= 0
        outlet_level = np.zeros(n)
        outlet_level[outlet] = np.asarray(levels, float)[water[outlet]]
        idx2 = np.arange(n, dtype=np.int32).reshape(rows, cols)
        receiver = dr.d8_receivers(h.astype(np.float32), outlet.reshape(rows, cols), lat, cell_km)
        sink = (receiver == idx2.ravel()) & ~outlet
        sinks = np.flatnonzero(sink).astype(np.int32)
        node_of = np.full(n, -1, np.int32)
        lv, groups = np.unique(outlet_level[outlet], return_inverse=True)
        n_out = max(1, lv.size)
        node_of[outlet] = groups
        node_of[sinks] = np.arange(n_out, sinks.size + n_out, dtype=np.int32)
        term = dr.terminal(receiver)
        label = node_of[term]
        n_nodes = sinks.size + n_out
        h_pass = np.where(outlet.reshape(rows, cols), -np.inf, h).astype(np.float32)
        edges = dr.passes(label.reshape(rows, cols), h_pass, idx2)
        level, spill_from, spill_to = flood_from(n_nodes, edges, list(enumerate(lv.tolist())))
        if np.isnan(level).any():
            raise RuntimeError(f'{int(np.isnan(level).sum())} catchments have no path to an outlet')
        del edges
        hp = h_pass.ravel()
        pass_height = np.full(n_nodes, -np.inf)
        pass_height[n_out:] = np.maximum(hp[spill_from[n_out:]], hp[spill_to[n_out:]])
        del h_pass, hp
        self.area = (((cell_km * 1000.0) ** 2 * np.maximum(np.cos(np.radians(lat)), 1e-6))[:, None]
                     * np.ones((1, cols))).ravel()
        hf = h.ravel()
        self.h = hf
        self.outlet, self.sink, self.sinks, self.node_of, self.label = outlet, sink, sinks, node_of, label
        self.level, self.spill_to, self.n_out, self.n_nodes = level, spill_to, n_out, n_nodes
        self.lake_full = ~outlet & (hf < level[label])
        final = receiver.copy()
        final[sinks] = spill_to[n_out:]
        self.receiver = final
        self.d8 = receiver
        # each catchment's parent, root (the outlet its water reaches) and depth below the root
        parent = np.arange(n_nodes, dtype=np.int64)
        parent[n_out:] = label[spill_to[n_out:]]
        self.parent = parent
        root, depth = parent.copy(), (parent != np.arange(n_nodes)).astype(np.int64)
        while True:                                  # pointer jumping: depth counts the hops to the current pointer
            nxt = root[root]
            if np.array_equal(nxt, root):
                break
            depth = depth + depth[root]
            root = nxt
        self.root, self.depth = root, depth
        # the lakes as bodies: a catchment whose pass lies under its parent's spill level is part of its parent's water
        # when that stands full, and the body's outlet is the catchment that spills over a pass of its own
        inner = np.zeros(n_nodes, bool)
        inner[n_out:] = (parent[n_out:] >= n_out) & (pass_height[n_out:] < level[parent[n_out:]])
        body = np.where(inner, parent, np.arange(n_nodes))
        while True:
            nxt = body[body]
            if np.array_equal(nxt, body):
                break
            body = nxt
        self.inner, self.body = inner, body
        # the water each cell's own path ends in (seas), and the water a catchment's spill reaches when its parent is a sea
        self.term_sea = np.where(outlet[term], water[term], -1).astype(np.int8)
        self.direct = np.flatnonzero((self.term_sea >= 0) & ~outlet)
        self.spill_sea = np.full(n_nodes, -1, np.int8)
        to_sea = np.arange(n_out, n_nodes)[parent[n_out:] < n_out]
        self.spill_sea[to_sea] = self.term_sea[spill_to[to_sea]]
        self.outlet_levels = lv
        # route_levels' generations over the cells
        donor = ~outlet
        indeg = np.bincount(final[donor], minlength=n).astype(np.int32)
        front = np.flatnonzero((indeg == 0) & donor)
        pieces, bounds = [], [0]
        while front.size:
            pieces.append(front.astype(np.int32))
            bounds.append(bounds[-1] + front.size)
            r = final[front]
            np.subtract.at(indeg, r, 1)
            cand = np.unique(r)
            front = cand[(indeg[cand] == 0) & donor[cand]]
        self.order = np.concatenate(pieces)
        self.bounds = np.asarray(bounds, np.int64)
        self.order_to = final[self.order]
        self.order_sink = sink[self.order]
        self.has_sink = np.add.reduceat(self.order_sink.astype(np.int64), self.bounds[:-1]) > 0
        del pieces, indeg
        log(f'network: {n_nodes - n_out:,} catchments, {self.bounds.size - 1} generations of cells, '
            f'{int(depth.max())} of catchments, {time.time() - started:.0f} s')

    def generations(self, nodes):
        """The catchments among `nodes` grouped by depth, deepest first: each group's parents lie in later groups."""
        nodes = np.asarray(nodes, np.int64)
        nodes = nodes[nodes >= self.n_out]
        d = self.depth[nodes]
        order = np.argsort(-d, kind='stable')
        nodes, d = nodes[order], d[order]
        cut = np.flatnonzero(np.diff(d)) + 1
        return np.split(nodes, cut)

    def sources(self, runoff_mm, lake_net_mm, ops=()):
        """The runoff on each cell and each catchment's lake balance (m3/s) from runoff and open-water balance fields
        (mm/day), with flows put in or taken out at cells as governed_drainage.route_with puts them: into a hollow's
        balance where `in_lake`, on the ground otherwise. ops: (cell, flow, in_lake)."""
        r64 = np.asarray(runoff_mm, np.float64).ravel().copy()
        l64 = np.asarray(lake_net_mm, np.float64).ravel().copy()
        for cell, q, in_lake in ops:
            if q == 0.0:
                continue
            if in_lake:
                l64[cell] += q / (MM_DAY * self.area[cell])
            elif q > 0.0:
                r64[cell] += q / (MM_DAY * self.area[cell])
            else:
                raise ValueError(f'a take-out at cell {cell} lies under no hollow\'s water')
        runoff = np.where(self.outlet | self.lake_full, 0.0, np.maximum(r64, 0.0) * MM_DAY * self.area)
        lake = np.bincount(self.label[self.lake_full], weights=(l64 * MM_DAY * self.area)[self.lake_full],
                           minlength=self.n_nodes)
        return runoff, lake

    def with_ops(self, base, runoff_mm, ops):
        """sources() with put-ins and take-outs, from the same fields' sources without them (`base`): each put-in cell's
        runoff taken again as route_with and route_levels take it, and each take-out from its hollow's balance."""
        runoff, lake = base[0].copy(), base[1].copy()
        cells = {}
        for cell, q, in_lake in ops:
            if q == 0.0:
                continue
            if in_lake:
                if self.lake_full[cell]:
                    lake[self.label[cell]] += q
            elif q > 0.0:
                cells.setdefault(cell, []).append(q)
            else:
                raise ValueError(f'a take-out at cell {cell} lies under no hollow\'s water')
        flat = np.asarray(runoff_mm).ravel()
        for cell, qs in cells.items():
            if self.outlet[cell] or self.lake_full[cell]:
                continue
            v = float(flat[cell])
            for q in qs:
                v += q / (MM_DAY * self.area[cell])
            runoff[cell] = max(v, 0.0) * MM_DAY * self.area[cell]
        return runoff, lake

    def cascade(self, runoff, lake, passed_on=None, groups=None, shut=None):
        """The catchments' inflows and what each passes on (m3/s), and what each water receives: route_levels' balance
        on the catchments alone, each passing on its inflow and lake balance where positive, or `passed_on` where
        given; the catchments marked `shut` pass nothing."""
        inflow = np.bincount(self.label, weights=runoff, minlength=self.n_nodes)
        out = np.zeros(self.n_nodes)
        for nodes in (groups if groups is not None else self.all_groups()):
            q = np.maximum(inflow[nodes] + lake[nodes], 0.0) if passed_on is None else passed_on[nodes]
            if shut is not None:
                q = np.where(shut[nodes], 0.0, q)
            out[nodes] = q
            np.add.at(inflow, self.parent[nodes], q)
        return inflow, out, self.supply(runoff, out)

    def all_groups(self):
        if not hasattr(self, '_groups'):
            self._groups = self.generations(np.arange(self.n_out, self.n_nodes))
        return self._groups

    def supply(self, runoff, out):
        """What each water receives: the runoff of the land draining straight into it and what the catchments spilling
        into it pass on."""
        direct = np.bincount(self.term_sea[self.direct], weights=runoff[self.direct], minlength=self.n_waters)
        into = np.flatnonzero(self.spill_sea >= 0)
        return direct + np.bincount(self.spill_sea[into], weights=out[into], minlength=self.n_waters)

    def supply_of(self, q):
        """What each water receives in a replay's discharge: the flow into its cells (route_at's count)."""
        if not hasattr(self, '_water_cells'):
            self._water_cells = [np.flatnonzero(self.water == i) for i in range(self.n_waters)]
        return np.array([float(q[c].sum()) for c in self._water_cells])

    def replay(self, runoff, lake, passed_on=None, loss=None, shut=None):
        """Each cell's discharge (m3/s): route_levels' accumulation over the kept generations, each pit passing on its
        inflow and lake balance where positive, or `passed_on` where given, the pits of catchments marked `shut` passing
        nothing; `loss` (m3/s per cell) is taken from what a cell passes on, never more than it carries. Returns the
        discharge and what each catchment passed on."""
        acc = runoff.copy()
        out_node = np.zeros(self.n_nodes)
        order, to, at_sink, b = self.order, self.order_to, self.order_sink, self.bounds
        for g in range(b.size - 1):
            i0, i1 = b[g], b[g + 1]
            cells = order[i0:i1]
            out = acc[cells]
            if loss is not None:
                out = np.maximum(out - loss[cells], 0.0)
            if self.has_sink[g]:
                s = at_sink[i0:i1]
                k = self.node_of[cells[s]]
                inflow = acc[cells[s]]
                out[s] = np.maximum(inflow + lake[k], 0.0) if passed_on is None else passed_on[k]
                if shut is not None:
                    out[s] = np.where(shut[k], 0.0, out[s])
                out_node[k] = out[s]
            np.add.at(acc, to[i0:i1], out)
        return acc, out_node


# ---------------------------------------------------------------- the years' climate

def year_sums():
    """The settled years' sums (moisture_map.read_year's, from drylands.py's cache) and the estimate's step with the
    moisture map's responses per kelvin (moisture_final.estimate_responses)."""
    bm.RUNS = DL.RUNS
    bm.CACHE = str(DL.CACHE)
    record, _, applied = MF.estimate_responses()
    years = {y: bm.cached_year(RUN, y, True) for y in range(YEARS[0], YEARS[1] + 1)}
    return years, applied, float(record['estimate']['delta_K'])


def year_climate(year, applied, delta, rotation):
    """One year's rain, evaporation and land on today's latitude and longitude: the year's sums moved by the estimate's
    step (moisture_final.adjusted), each field the mean over the year's records as climate.gcm.climatology takes it
    (evaporation turned positive upward), carried from the settled frame as climate.gcm.turned_climatology.to_today
    carries a climatology."""
    y = MF.adjusted(year, applied, delta)
    count = y['cnt'].sum(axis=0)
    pr = y['acc']['pr'].sum(axis=0) / count * bm.MM
    ev = y['acc']['evap'].sum(axis=0) / count * bm.MM
    ev = -ev if ev.mean() < 0 else ev
    fields = dict(pr_mm_day=pr.astype(np.float32), evap_mm_day=ev.astype(np.float32),
                  lsm=np.asarray(y['land'], float).astype(np.float32))
    return carried(fields, y['lat'], y['lon'], rotation)


def carried(fields, lat, lon, rotation):
    """Fields of the settled frame's grid on the same grid in today's frame (turned_climatology.to_today's sampling),
    the grid's coordinates in single precision as the climatology products hold them."""
    lat32, lon32 = np.asarray(lat, np.float32), np.asarray(lon, np.float32)
    lat, lon = np.asarray(lat32, float), np.asarray(lon32, float)
    grid_lat, grid_lon = np.meshgrid(lat, lon, indexing='ij')
    lat_t, lon_t = to_turned(grid_lat, grid_lon, rotation)
    out = {k: TC.sample(v, lat, lon, lat_t, lon_t).astype(np.float32) for k, v in fields.items()}
    out.update(lat=lat32, lon=lon32)
    return out


def loaded(clim, evaporation, rows, cols):
    """A climate's runoff, open-water balance, rain and evaporation (mm/day) on the 16 pixel/degree grid, as
    final_drainage.load takes them from the design climatology on today's ground."""
    lat = 90.0 - (np.arange(rows) + 0.5) / dr.PPD
    precip = (dr.climate_on_grid(clim, 'pr_mm_day', rows, cols) * np.float32(1.0)).astype(np.float32)
    evap = dr.climate_on_grid(clim, 'evap_mm_day', rows, cols).astype(np.float32)
    runoff = np.maximum(precip - evap, 0.0).astype(np.float32)
    lake_net = (precip - np.float32(evaporation) * dr.lake_evaporation(clim, lat)[:, None].astype(np.float32)).astype(np.float32)
    return runoff, lake_net, precip, evap


# ---------------------------------------------------------------- a year without storage

def op_flows(s, mosc, orientale, smythii, neck):
    """Each put-in's and take-out's flow from the outflows that feed it: Moscoviense's overflow divided by its fixed
    crests in the design's shares (register, 17 November), Orientale's and Smythii-Marginis's whole, and the neck's
    spill."""
    f = s['meta']['fractions']
    unit = dict(salt_sea_weir=mosc * f[GD.SALT_SEA], humboldtianum_branch=mosc * f[GD.HELD_BRANCH],
                smythii_marginis_branch=mosc * f[GD.REST_BRANCH], orientale=orientale, smythii_marginis=smythii,
                neck=neck)
    return s['op_sign'] * np.array([unit[k] for k in OPS])[s['op_kind']]


def ops_list(s, flows):
    return list(zip(s['op_cell'].tolist(), np.asarray(flows, float).tolist(), s['op_lake'].tolist()))


def water_losses(net, lake_net_mm):
    """Each water's loss to the air beyond the rain on it (m3/s) at its design level, as redrawn_drainage.design_state
    counts a sea's loss."""
    field = -np.asarray(lake_net_mm, np.float64).ravel() * MM_DAY * net.area
    return net.supply_of(field)


def without_storage(net, s, runoff_mm, lake_net_mm, reference=None, river_loss=None, shut=None):
    """One year routed on the fixed structure with no water carried in or out: each upper sea passes on what its rivers
    bring it beyond its loss that year, Moscoviense's by the fixed shares, and the main part of the near-side sea sends
    over its neck the design's spill moved by its balance's departure from the reference year (`reference`, the design
    climate's own routing; the design's spill itself when None). The catchments marked `shut`, the hollows the design
    leaves below their passes, pass nothing: held as the design fills them, a wet year cannot raise them to their passes.
    `river_loss` (m3/s a cell) is taken from the rivers where it falls, so a lake or sea below receives less. Returns
    the discharge, what each catchment passes on, each water's supply and loss, the outflows and the main part's
    supply."""
    iM, iO, iS, iN = (WATERS.index(k) for k in ('moscoviense', 'orientale', 'smythii_marginis', 'near_side_main'))
    loss = water_losses(net, lake_net_mm)
    base = net.sources(runoff_mm, lake_net_mm)

    def routed(mosc, orientale, smythii, neck, whole=False):
        runoff, lake = net.with_ops(base, runoff_mm, ops_list(s, op_flows(s, mosc, orientale, smythii, neck)))
        if river_loss is None and not whole:
            return None, None, net.cascade(runoff, lake, shut=shut)[2]
        q, out = net.replay(runoff, lake, loss=river_loss, shut=shut)
        return q, out, net.supply_of(q)

    _, _, own = routed(0.0, 0.0, 0.0, 0.0)
    mosc, orientale = max(own[iM] - loss[iM], 0.0), max(own[iO] - loss[iO], 0.0)
    _, _, sup = routed(mosc, orientale, 0.0, 0.0)
    smythii = max(sup[iS] - loss[iS], 0.0)
    _, _, sup = routed(mosc, orientale, smythii, 0.0)
    neck = s['meta']['neck_flow']
    if reference is not None:
        neck = max(neck + (sup[iN] - reference['main_supply']) - (loss[iN] - reference['loss'][iN]), 0.0)
    q, out, supply = routed(mosc, orientale, smythii, neck, whole=True)
    return dict(q=q, out=out, supply=supply, loss=loss, main_supply=float(sup[iN]),
                outflow=dict(moscoviense=mosc, orientale=orientale, smythii_marginis=smythii, neck=neck))


# ---------------------------------------------------------------- the years with storage

def neck_rating(level0):
    """The main part's spill over the neck against its level: critical flow over the throat (near_side_neck.Neck from
    near_side_neck.json's throat), as the steady state passes it, choked; returns the crest and a table of heads above
    it (m) and flows (m3/s)."""
    rec = json.loads(N.OUT.read_text())
    crest = level0 + rec['crest']['four_connected_m']
    by = rec['throat']['by_level']
    neck = N.Neck(crest, [crest] + [level0 + b['level_m'] for b in by], [0.0] + [b['section_m2'] for b in by],
                  [0.0] + [b['width_km'] * 1e3 for b in by], rec['sill']['length_km'] * 1e3)
    heads = np.concatenate([np.linspace(0.0, 2.0, 81)[1:], np.linspace(2.0, 30.0, 113)[1:]])
    return crest, np.concatenate([[0.0], heads]), np.array([0.0] + [neck.critical(crest + e) for e in heads])


def solve_head(area, c, s):
    """The head h >= 0 over an outlet at the end of a step: area * h + c * h^1.5 = s, s > 0 (an implicit step of a
    water whose outflow over the step is c * h^1.5 / dt). Newton's method on u = h^0.5 from an upper bound, which a
    convex rising function approaches from above."""
    u = np.minimum(np.sqrt(s / np.maximum(area, 1e-30)), np.cbrt(s / np.maximum(c, 1e-30)))
    for _ in range(60):
        f = area * u * u + c * u ** 3 - s
        step = f / np.maximum(2.0 * area * u + 3.0 * c * u * u, 1e-300)
        u = np.maximum(u - step, 0.0)
        if np.all(np.abs(step) <= 1e-12 * np.maximum(u, 1e-12)):
            break
    return u * u


def row_interp(x, xs, ys):
    """Linear interpolation of each row's table (xs falling along the row) at that row's x, held at the ends."""
    j = np.clip((xs >= x[:, None]).sum(axis=1) - 1, 0, xs.shape[1] - 2)
    r = np.arange(x.size)
    x0, x1, y0, y1 = xs[r, j], xs[r, j + 1], ys[r, j], ys[r, j + 1]
    t = np.clip(np.where(x0 > x1, (x0 - x) / np.where(x0 > x1, x0 - x1, 1.0), 0.0), 0.0, 1.0)
    return y0 + t * (y1 - y0)


class Storage:
    """The lakes' and the upper seas' storage from step to step on the fixed network. Each lake is a body of catchments
    standing at one level (Network.body); its outlet catchment holds the body's volume against the spill level, and the
    catchments inside it pass their water on as the routing has them. Above the spill level the body passes on critical
    flow over an outlet of width `width` (m per (m3/s)^0.5, river_losses.channel_width_m at the design's outflow, at
    least Q_FLOOR) at its head over the water the design holds at that level (the catchments that spill, joined to the
    outlet through catchments that spill); below it, it passes nothing and gains or loses the rain and evaporation of the
    outlet catchment's cells its level covers, tabled at DEPTHS_M. Moscoviense, Orientale and Smythii-Marginis hold their volume over their
    rims the same way, Moscoviense's outflow divided by its fixed shares; the main part of the near-side sea spills over
    the neck by neck_rating, its balance less the fixed offset (main_part_balance) by which the design routing brings it
    more than the neck study's steady state, so that the design climate keeps the design's spill. The salt sea,
    Humboldtianum and the southern seas pass nothing on. Each step routes the
    catchments upstream first, in stages that follow the seas' order (Moscoviense and Orientale; the salt sea,
    Humboldtianum and Smythii-Marginis; the main part; the southern seas), so each put-in is known before the
    catchments it feeds."""

    def __init__(self, net, s, design, width, log=print):
        started = time.time()
        self.net, self.s, self.width = net, s, width
        m = s['meta']
        n_out = net.n_out
        h = net.h.astype(np.float64)
        # the lakes' tables: area, volume against the spill level, and the cells' order for the balance below it
        cells = np.flatnonzero(net.lake_full)
        lab = net.label[cells].astype(np.int64)
        hc = h[cells]
        order = np.lexsort((hc, lab))
        cells, lab, hc = cells[order], lab[order], hc[order]
        a = net.area[cells]
        lo = float(hc.min()) - 1.0
        span = float(max(np.nanmax(net.level), hc.max()) - lo) + 1.0
        key = lab * span + (hc - lo)
        nodes = np.arange(net.n_nodes)
        start = np.searchsorted(lab, nodes, side='left')
        end = np.searchsorted(lab, nodes, side='right')
        lake_nodes = nodes[(nodes >= n_out) & (end > start) & ~net.inner]       # each body's outlet catchment
        depth = np.asarray(DEPTHS_M)
        z = net.level[lake_nodes][:, None] - depth[None, :]
        pos = np.searchsorted(key, lake_nodes[:, None] * span + (z - lo))
        pos = np.clip(pos, start[lake_nodes][:, None], end[lake_nodes][:, None])
        cum_a = np.concatenate([[0.0], np.cumsum(a)])
        cum_ah = np.concatenate([[0.0], np.cumsum(a * hc)])
        s0 = start[lake_nodes][:, None]
        area_j = cum_a[pos] - cum_a[s0]
        vol_j = z * area_j - (cum_ah[pos] - cum_ah[s0])
        self.cells, self.pos, self.s0 = cells, pos, s0
        self.lake_nodes = lake_nodes
        self.slot = np.full(net.n_nodes, -1, np.int64)
        self.slot[lake_nodes] = np.arange(lake_nodes.size)
        # the head stands over the body's water the design holds at the spill level: each catchment that spills, joined
        # to the outlet through catchments that spill too (a closed one inside holds a smaller lake of its own, lower)
        node_area = np.bincount(net.label[net.lake_full], weights=net.area[net.lake_full], minlength=net.n_nodes)
        spills = design['out'] > 0.0
        joined = np.zeros(net.n_nodes, bool)
        for group in reversed(net.all_groups()):                 # the outlets first, then the catchments inside
            joined[group] = spills[group] & np.where(net.inner[group], joined[net.parent[group]], True)
        self.area_head = np.maximum(np.bincount(net.body[joined], weights=node_area[joined], minlength=net.n_nodes),
                                    node_area)[lake_nodes]
        self.vol_full = vol_j[:, 0]
        self.vtab = vol_j - vol_j[:, :1]                       # volume against the spill level, falling to -vol_full
        # the outlets: critical flow c_out * h^1.5 over a channel as wide as the design's outflow makes it
        q_design = design['out'][lake_nodes]
        self.k_out = RL.channel_width_m(np.maximum(q_design, Q_FLOOR), width) * math.sqrt(G) * (2.0 / 3.0) ** 1.5
        # the upper seas over their rims, and the main part over the neck
        iw = {k: WATERS.index(k) for k in WATERS}
        self.sea_area = np.array([net.area[net.water == i].sum() for i in range(net.n_waters)])
        self.sea_q = {'moscoviense': design['outflow']['moscoviense'], 'orientale': design['outflow']['orientale'],
                      'smythii_marginis': design['outflow']['smythii_marginis']}
        self.sea_k = {k: float(RL.channel_width_m(np.array([max(q, Q_FLOOR)]), width)[0]) * math.sqrt(G)
                      * (2.0 / 3.0) ** 1.5 for k, q in self.sea_q.items()}
        self.crest, self.neck_head, self.neck_q = neck_rating(m['level0'])
        # the design routing brings the main part more than the neck study's steady state took (its rivers run on over
        # the floor the parted sea leaves dry); held as a fixed offset, so the design climate keeps the design's spill
        iN = WATERS.index('near_side_main')
        self.main_offset = float(design['main_supply'] - design['loss'][iN] - design['outflow']['neck'])
        # the stages: each catchment by the water its own water reaches
        stage_of_water = dict(moscoviense=0, orientale=0, smythii_marginis=1, south_pole_aitken=1, humboldtianum=1,
                              near_side_main=2, near_side_south=3)
        stage_of_outlet = {}
        for k, i in iw.items():
            node = int(net.node_of[np.flatnonzero(net.water == i)[0]])
            if stage_of_outlet.setdefault(node, stage_of_water[k]) != stage_of_water[k]:
                raise SystemExit(f'{k} shares an outlet with a water of another stage')
        stage = np.full(net.n_nodes, -1, np.int64)
        for node, st_ in stage_of_outlet.items():
            stage[(net.root == node) & (nodes >= n_out)] = st_
        self.stage = stage
        self.groups = [net.generations(np.flatnonzero(stage == i)) for i in range(4)]
        self.to_sea = [np.flatnonzero((stage == i) & (net.spill_sea >= 0)) for i in range(4)]
        # where each put-in lands: a catchment's balance, a water straight, or nowhere (route_levels drops it)
        oc, ol = s['op_cell'], s['op_lake']
        lands_lake = ol & net.lake_full[oc]
        lands_ground = ~ol & ~net.outlet[oc] & ~net.lake_full[oc]
        target = np.where(lands_lake | lands_ground, net.label[oc], -1)
        self.op_node = np.where(target >= n_out, target, -1)
        self.op_sea = np.where((target >= 0) & (target < n_out), net.term_sea[oc], -1)
        kinds_after = {0: ('salt_sea_weir', 'humboldtianum_branch', 'smythii_marginis_branch', 'orientale'),
                       1: ('smythii_marginis',), 2: ('neck',)}
        for st_, kinds in kinds_after.items():
            for kind in kinds:
                sel = s['op_kind'] == OPS.index(kind)
                bad = sel & (self.op_node >= 0) & (stage[np.maximum(self.op_node, 0)] <= st_)
                if bad.any():
                    raise SystemExit(f'a put-in or take-out of {kind} lands in a catchment routed before it')
        self.iw = iw
        log(f'storage at width {width:g}: {lake_nodes.size:,} lakes, {sum(len(g) for g in self.groups)} groups, '
            f'{time.time() - started:.0f} s')

    def initial(self, design, lake_design, inflow_design):
        """The design state: each spilling lake at the head its design outflow needs, each closed lake where its loss
        meets its inflow, the upper seas at the heads their design overflows need and the main part where the neck
        passes the design's spill."""
        q = design['out'][self.lake_nodes]
        head = np.where(q > 0.0, (q / self.k_out) ** (2.0 / 3.0), 0.0)
        vol = self.area_head * head
        closed = q <= 0.0
        if closed.any():
            ntab = self.table(lake_design)[closed]
            need = -inflow_design[self.lake_nodes][closed]
            vol[closed] = self._level_for(ntab, self.vtab[closed], need)
        self.vol = vol
        self.sea_vol = {k: self.sea_area[self.iw[k]] * (self.sea_q[k] / self.sea_k[k]) ** (2.0 / 3.0)
                        for k in self.sea_q}
        e = float(np.interp(design['outflow']['neck'], self.neck_q, self.neck_head))
        self.main_vol = self.sea_area[self.iw['near_side_main']] * e

    @staticmethod
    def _level_for(ntab, vtab, need):
        """The volume at which each lake's balance below its spill level meets its inflow: the first tabled level, going
        down, at which the lake's loss no longer exceeds the inflow, interpolated in volume; at the spill level where
        the full lake's loss does not exceed it, empty where no level's does."""
        ok = ntab >= need[:, None]
        j = np.argmax(ok, axis=1)
        r = np.arange(need.size)
        jm = np.maximum(j - 1, 0)
        n0, n1, v0, v1 = ntab[r, jm], ntab[r, j], vtab[r, jm], vtab[r, j]
        t = np.where(n1 != n0, (need - n0) / np.where(n1 != n0, n1 - n0, 1.0), 1.0)
        out = v0 + np.clip(t, 0.0, 1.0) * (v1 - v0)
        out = np.where(j == 0, 0.0, out)
        return np.where(ok.any(axis=1), out, vtab[:, -1])

    def table(self, lake_net_mm):
        """Each lake's balance (m3/s, rain less evaporation) over the cells below each tabled level, for a year."""
        w = np.asarray(lake_net_mm, np.float64).ravel()[self.cells] * MM_DAY * self.net.area[self.cells]
        cum = np.concatenate([[0.0], np.cumsum(w)])
        return cum[self.pos] - cum[self.s0]

    def sea_step(self, name, vol, inflow, loss, dt):
        """An upper sea's volume over its rim and its outflow over a step."""
        a, k = self.sea_area[self.iw[name]], self.sea_k[name]
        s = vol + dt * (inflow - loss)
        if s <= 0.0:
            return s, 0.0
        h = float(solve_head(np.array([a]), np.array([dt * k]), np.array([s]))[0])
        return a * h, (s - a * h) / dt

    def main_step(self, vol, inflow, loss, dt):
        """The main part's volume over the neck's crest and its spill over a step (bisection on the tabled rating)."""
        a = self.sea_area[self.iw['near_side_main']]
        s = vol + dt * (inflow - loss)
        if s <= 0.0:
            return s, 0.0
        lo, hi = 0.0, s / a
        for _ in range(80):
            e = 0.5 * (lo + hi)
            if a * e + dt * float(np.interp(e, self.neck_head, self.neck_q)) > s:
                hi = e
            else:
                lo = e
        e = 0.5 * (lo + hi)
        return a * e, (s - a * e) / dt

    def year(self, runoff_mm, lake_net_mm, steps=STEPS):
        """One year in `steps` steps; returns the mean of what each catchment passed on, the waters' mean outflows and
        supplies, and their losses."""
        net, s = self.net, self.s
        dt = YEAR_S / steps
        runoff, lake = net.sources(runoff_mm, lake_net_mm)
        local = np.bincount(net.label, weights=runoff, minlength=net.n_nodes)
        direct = np.bincount(net.term_sea[net.direct], weights=runoff[net.direct], minlength=net.n_waters)
        loss = water_losses(net, lake_net_mm)
        ntab = self.table(lake_net_mm)
        slot = self.slot
        total = np.zeros(net.n_nodes)
        totals = dict(moscoviense=0.0, orientale=0.0, smythii_marginis=0.0, neck=0.0)
        supply_total = np.zeros(net.n_waters)
        for _ in range(steps):
            inflow = local.copy()
            sea_in = direct.copy()
            out = np.zeros(net.n_nodes)
            flows = dict(moscoviense=0.0, orientale=0.0, smythii_marginis=0.0, neck=0.0)

            def put(kinds):
                q = op_flows(s, flows['moscoviense'], flows['orientale'], flows['smythii_marginis'], flows['neck'])
                sel = np.isin(s['op_kind'], [OPS.index(k) for k in kinds])
                node, sea = self.op_node[sel], self.op_sea[sel]
                np.add.at(inflow, node[node >= 0], q[sel][node >= 0])
                np.add.at(sea_in, sea[sea >= 0], q[sel][sea >= 0])

            for st_ in range(4):
                if st_ == 1:
                    put(('salt_sea_weir', 'humboldtianum_branch', 'smythii_marginis_branch'))
                elif st_ == 2:
                    put(('orientale', 'smythii_marginis'))
                elif st_ == 3:
                    put(('neck',))
                for nodes in self.groups[st_]:
                    i = slot[nodes]
                    q = np.maximum(inflow[nodes] + lake[nodes], 0.0)
                    has = i >= 0
                    if has.any():
                        j = i[has]
                        v = self.vol[j]
                        bal = np.where(v >= 0.0, lake[nodes[has]], 0.0)
                        below = v < 0.0
                        if below.any():
                            bal[below] = row_interp(v[below], self.vtab[j[below]], ntab[j[below]])
                        sv = v + dt * (inflow[nodes[has]] + bal)
                        up = sv > 0.0
                        new = np.maximum(sv, -self.vol_full[j])
                        qh = np.zeros(j.size)
                        if up.any():
                            ah = self.area_head[j[up]]
                            hh = solve_head(ah, dt * self.k_out[j[up]], sv[up])
                            new[up] = ah * hh
                            qh[up] = (sv[up] - new[up]) / dt
                        self.vol[j] = new
                        q[has] = qh
                    out[nodes] = q
                    np.add.at(inflow, net.parent[nodes], q)
                into = self.to_sea[st_]
                sea_in += np.bincount(net.spill_sea[into], weights=out[into], minlength=net.n_waters)
                if st_ == 0:
                    for k in ('moscoviense', 'orientale'):
                        self.sea_vol[k], flows[k] = self.sea_step(k, self.sea_vol[k], sea_in[self.iw[k]],
                                                                  loss[self.iw[k]], dt)
                elif st_ == 1:
                    k = 'smythii_marginis'
                    self.sea_vol[k], flows[k] = self.sea_step(k, self.sea_vol[k], sea_in[self.iw[k]],
                                                              loss[self.iw[k]], dt)
                elif st_ == 2:
                    iN = self.iw['near_side_main']
                    self.main_vol, flows['neck'] = self.main_step(self.main_vol, sea_in[iN] - self.main_offset,
                                                                  loss[iN], dt)
            total += out
            for k in totals:
                totals[k] += flows[k]
            supply_total += sea_in
        return dict(out=total / steps, outflow={k: v / steps for k, v in totals.items()},
                    supply=supply_total / steps, loss=loss, sources=(runoff, lake))

    def field(self, runoff_mm, lake_net_mm, result):
        """A year's discharge in each cell from its mean outflows (a replay with each pit passing on its mean)."""
        o = result['outflow']
        ops = ops_list(self.s, op_flows(self.s, o['moscoviense'], o['orientale'], o['smythii_marginis'], o['neck']))
        runoff, lake = self.net.with_ops(result['sources'], runoff_mm, ops)
        q, _ = self.net.replay(runoff, lake, passed_on=result['out'])
        return q


# ---------------------------------------------------------------- the drylands' losses on the design climate

def dryland_losses(net, s, q, lake_cell, budgets, frame, log=print):
    """river_losses.py's estimate run again on the design routing's discharge `q` and the design climate's open-water
    balance on each budget (`budgets`: label to lake_net, mm/day): the dryland rivers of at least river_losses.Q_RIVER,
    their channels at both width coefficients, the belt of plants drinking their beds' seepage at both widths, the
    floodplains by height above the nearest river, the connected floodplains at river_losses.HAND_MAIN with the area
    each river can keep wet, and those carried onto the design surface's T21 grid against PlaSim's mask (`frame`: the
    design surface's record). Returns the record, in river_losses.json's form, and each cell's loss (m3/s) on the closed
    budget at the low end (the narrow channel and belt) and the high end (the wide channel and belt and the floodplains
    kept wet), each taken from the river cell it lies on or beside."""
    from climate.gcm import boundary as bd
    from research.studies.atlas_round import turned_surface as TS
    rows, cols = net.rows, net.cols
    lat, cell_km = np.asarray(s['lat'], float), s['meta']['cell_km']
    h = s['h']
    receiver = net.receiver.astype(np.int64)
    area_km2 = np.repeat(cell_km ** 2 * np.cos(np.radians(lat)), cols)
    land = net.water < 0
    land[lake_cell] = False
    length = RL.flow_length_km(receiver, rows, cols, lat, cell_km)
    river = land & (q >= RL.Q_RIVER)
    near = RL.nearest_river(river.reshape(rows, cols))
    record = dict(budgets={})
    main = None
    for label, lake_net in budgets.items():
        loss = np.maximum(-np.asarray(lake_net).ravel().astype(float), 0.0) * MM_DAY
        dry = loss > 0.0
        dry_river = river & dry
        out = dict(rivers={}, floodplains={})
        for qmin in RL.Q_CLASSES:
            m = dry_river & (q >= qmin)
            out['rivers'][f'{qmin:.0f}'] = dict(length_km=round(float(length[m].sum()), 0),
                                                largest_flow_m3s=round(float(q[m].max()), 0) if m.any() else None)
        lengths, losses = length[dry_river], loss[dry_river]
        out['channel'] = {}
        for coef in RL.WIDTH:
            width_km = RL.channel_width_m(q[dry_river], coef) / 1000.0
            out['channel'][f'{coef:g}'] = dict(area_km2=round(float((width_km * lengths).sum()), 0),
                                               evaporation_m3s=round(float((width_km * lengths * 1e6 * losses).sum()), 1))
        out['riparian_belt'] = {f'{b:g}': dict(area_km2=round(float((b * lengths).sum()), 0),
                                               evaporation_m3s=round(float((b * lengths * 1e6 * losses).sum()), 1))
                                for b in RL.RIPARIAN_KM}
        out['dry_land_share_of_land'] = round(float(area_km2[dry & land].sum() / area_km2[land].sum()), 4)
        for hand in RL.HAND_M:
            plain = RL.floodplain(h, near, dry, land, lat, cell_km, hand)
            out['floodplains'][f'{hand:g}'] = dict(
                area_km2=round(float(area_km2[plain].sum()), 0),
                evaporation_if_all_wet_m3s=round(float((area_km2[plain] * 1e6 * loss[plain]).sum()), 1))
            if label == BUDGET and hand == RL.HAND_MAIN:
                main = (plain, loss, dry)
        record['budgets'][label] = out
        log(f'{label}: dryland rivers {out["rivers"]["100"]["length_km"]:.0f} km, channels '
            f'{out["channel"]["3"]["evaporation_m3s"]}-{out["channel"]["15"]["evaporation_m3s"]} m3/s')
    from scipy import ndimage
    plain, loss, dry = main
    labels, count = ndimage.label(plain.reshape(rows, cols), structure=np.ones((3, 3), bool))
    labels = labels.ravel()
    cells = np.flatnonzero(labels)
    lab = labels[cells]
    area = np.bincount(lab, weights=area_km2[cells], minlength=count + 1)
    need = np.bincount(lab, weights=area_km2[cells] * 1e6 * loss[cells], minlength=count + 1)
    inflow = np.zeros(count + 1)
    np.maximum.at(inflow, lab, q[near[cells]])
    mean_loss = need / np.maximum(area * 1e6, 1e-30)
    capped = RL.capped_area_km2(area, mean_loss, inflow)
    keep = capped[lab] / np.maximum(area[lab], 1e-30)
    rotation = np.asarray(frame['frame']['rotation_today_to_turned'], float)
    t21 = bd.read_t21(TS.PRODUCTS / frame['surface'])
    cell_area = t21['area'] / t21['area'].sum() * 4.0 * np.pi * (MOON_RADIUS / 1e3) ** 2
    cut = float(t21['land'][t21['mask'] == 0].max())
    lon = (np.arange(cols) + 0.5) / dr.PPD
    tr, tc = RL.t21_cells(lat[cells // cols], lon[cells % cols], rotation)
    added = np.zeros(t21['land'].shape)
    np.add.at(added, (tr, tc), area_km2[cells] * keep)
    fraction = added / cell_area
    flips = (t21['mask'] > 0) & (t21['land'] - fraction < cut)
    biggest = np.argsort(area[1:])[::-1][:RL.TOP] + 1
    listed = []
    for k in biggest:
        members = cells[lab == k]
        r_, c_ = members // cols, members % cols
        la, lo = float(lat[r_].mean()), float(np.degrees(np.angle(np.exp(1j * np.radians(lon[c_])).mean())))
        lt, ln = to_turned(la, lo, rotation)
        listed.append(dict(area_km2=round(float(area[k]), 0), today=[round(la, 2), round(lo, 2)],
                           turned=[round(float(lt), 2), round(float(ln), 2)], river_flow_m3s=round(float(inflow[k]), 0),
                           net_evaporation_mm_year=round(float(mean_loss[k] * 1000.0 * RL.SECONDS_PER_YEAR), 0),
                           kept_wet_km2=round(float(capped[k]), 0),
                           evaporation_m3s=round(float(min(need[k], inflow[k])), 1)))
    record.update(
        floodplains_at=dict(hand_m=RL.HAND_MAIN, budget=BUDGET, count=int(count), area_km2=round(float(area[1:].sum()), 0),
                            kept_wet_km2=round(float(capped[1:].sum()), 0),
                            evaporation_m3s=round(float(np.minimum(need[1:], inflow[1:]).sum()), 1), largest=listed),
        gcm=dict(surface=frame['surface'], water_cells=int((t21['mask'] == 0).sum()),
                 mask_land_fraction_at_most=round(cut, 4), largest_added_water_fraction=round(float(fraction.max()), 4),
                 cells_with_added_water_over_1_percent=int((fraction > 0.01).sum()),
                 land_cells_that_would_turn_to_water=int(flips.sum())))
    # each loss taken from the river where it falls: the channel and belt on their own cells, a floodplain's from the
    # river cell nearest each of its cells
    dry_river = np.flatnonzero(river & dry)
    base = length[dry_river] * 1e6 * loss[dry_river]
    width = {c: RL.channel_width_m(q[dry_river], c) / 1000.0 for c in RL.WIDTH}
    low, high = np.zeros(net.n), np.zeros(net.n)
    low[dry_river] = (width[RL.WIDTH[0]] + RL.RIPARIAN_KM[0]) * base
    high[dry_river] = (width[RL.WIDTH[-1]] + RL.RIPARIAN_KM[-1]) * base
    np.add.at(high, near[cells], area_km2[cells] * 1e6 * loss[cells] * keep)
    return record, dict(low=low, high=high)


# ---------------------------------------------------------------- the major rivers' reaches

def lake_bodies(wet, rows, cols, area):
    """The lakes off the seas as bodies of eight-connected water cells, joined across the 0/360 seam, with each body's
    area (km2); body 0 is no lake."""
    from scipy import ndimage
    lab, n = ndimage.label(wet.reshape(rows, cols), structure=np.ones((3, 3), bool))
    parent = np.arange(n + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for dr_ in (-1, 0, 1):                              # the last column's cells touch the first column's
        a = lab[:, -1]
        b = np.roll(lab[:, 0], -dr_) if dr_ else lab[:, 0]
        if dr_ == 1:
            b[-1] = 0
        elif dr_ == -1:
            b[0] = 0
        for x, y in zip(a[(a > 0) & (b > 0)], b[(a > 0) & (b > 0)]):
            rx, ry = find(int(x)), find(int(y))
            if rx != ry:
                parent[max(rx, ry)] = min(rx, ry)
    roots = np.array([find(i) for i in range(n + 1)])
    lab = roots[lab.ravel()]
    return lab, np.bincount(lab, weights=area / 1e6, minlength=n + 1)


def next_down(cells, q, wet, outlet, water, receiver, sink, node_of, passes_on, body, body_km2, threshold,
              break_km2=LAKE_BREAK_KM2):
    """Where each major river cell's water goes next: the next major river cell, through the lakes smaller than
    `break_km2` that lie on its course and pass water on; or the end of the river, into a sea (its index), a lake of at
    least `break_km2` or a lake that passes nothing on (its body), or onto land where it carries less than the
    threshold. Returns (kind, value) by cell, kind one of 'river', 'sea', 'lake', 'closed', 'below'."""
    out = {}
    for c in cells.tolist():
        d = int(receiver[c])
        for _ in range(200000):
            if outlet[d]:
                out[c] = ('sea', int(water[d]))
                break
            if not wet[d]:
                out[c] = ('river', d) if q[d] >= threshold else ('below', d)
                break
            b = int(body[d])
            if body_km2[b] >= break_km2:
                out[c] = ('lake', b)
                break
            if sink[d] and passes_on[node_of[d]] <= 0.0:
                out[c] = ('closed', b)
                break
            d = int(receiver[d])
        else:
            raise RuntimeError(f'the water from cell {c} finds no end')
    return out


def stems(cells, q, down):
    """The major rivers as stems: from each head (a major river cell no other feeds) down to where the stem ends in a
    sea, a lake or the land, or joins a river that carries more. At a confluence the stem carrying the most flow goes on.
    Returns each stem's cells, head first, and how it ends."""
    donors = {}
    for c in cells.tolist():
        kind, d = down[c]
        if kind == 'river':
            donors.setdefault(d, []).append(c)
    main = {d: max(cs, key=lambda c: q[c]) for d, cs in donors.items()}
    out = []
    for head in (c for c in cells.tolist() if c not in donors):
        path = [head]
        while True:
            kind, d = down[path[-1]]
            if kind == 'river' and main[d] == path[-1]:
                path.append(d)
                continue
            out.append(dict(cells=path, ends=(('joins', d) if kind == 'river' else (kind, d))))
            break
    return out


def classify(flows, threshold):
    """Each cell's yearly flows (years by cells) read against a threshold: the mean, the lowest and second-lowest
    years, the years at or above the threshold, whether it is sustained (at or above it in every year) and whether it
    stops flowing (under DRY_M3S) in any year."""
    ranked = np.sort(flows, axis=0)
    return dict(mean=flows.mean(axis=0), lowest=ranked[0], second=ranked[1],
                years=(flows >= threshold).sum(axis=0), sustained=ranked[0] >= threshold, stops=ranked[0] < DRY_M3S)


# ---------------------------------------------------------------- the run

def year_files():
    """The year sums the run reads: the design climate's settled years and the windows of the estimate's pairs, in
    drylands.py's cache."""
    est = json.loads((ROOT / FD.ESTIMATE).read_text())['estimate']
    windows = [(RUN, list(YEARS))] + [(p[w]['run'], p[w]['years']) for p in est['pairs'] for w in ('warm', 'cool')]
    return sorted({DL.CACHE / f'{run}.{y:05d}.1.npz' for run, (a, b) in windows for y in range(a, b + 1)})


def ops_of(s, outflow):
    return ops_list(s, op_flows(s, outflow['moscoviense'], outflow['orientale'], outflow['smythii_marginis'],
                                outflow['neck']))


def year_summary(supply, loss, outflow, q, into):
    """A year's waters: what each receives and loses (m3/s), the upper seas' and the neck's outflows and the rivers'
    flow into the seas."""
    return dict(supply_m3s={k: round(float(v), 1) for k, v in zip(WATERS, supply)},
                loss_m3s={k: round(float(v), 1) for k, v in zip(WATERS, loss)},
                outflow_m3s={k: round(float(v), 1) for k, v in outflow.items()},
                rivers_into_the_seas_m3s=round(float(q[into].sum()), 1))


def checkpoint(path, done, flows, summaries, stores):
    """The years finished so far, with the storage's state, kept on the research drive."""
    WORK.mkdir(parents=True, exist_ok=True)
    state = {f'vol_{k}': st.vol for k, st in stores.items()}
    meta = dict(done=done, summaries=summaries,
                seas={k: dict(sea_vol=st.sea_vol, main_vol=st.main_vol) for k, st in stores.items()})
    np.savez(path, metadata=json.dumps(meta), **{f'flow_{k}': v for k, v in flows.items()}, **state)


def run_years(net, s, reference, design_in, closed, track, into, land, log=print, resume=False):
    """The settled years in order: each routed without storage ('none', the hollows the design leaves closed held
    shut), as a steady state of its own ('steady', those hollows filling and spilling in a wet year as a steady routing
    has them), and with the storage at both outlet widths ('narrow', 'wide'), the storage starting from the design state.
    Returns each scenario's flows at the tracked cells (years by cells), each year's summaries, the storage's check on
    the design climate, the years' mean climate on the climate model's grid and each year's runoff from the land
    (m3/s, the design's dry ground)."""
    m = s['meta']
    rotation = np.asarray(m['rotation'], float)
    rows, cols = net.rows, net.cols
    years = list(range(YEARS[0], YEARS[1] + 1))
    sums, applied, delta = year_sums()
    runoff_d, lake_d = design_in
    stores, steady = {}, {}
    for label, width in zip(('narrow', 'wide'), RL.WIDTH):
        st = Storage(net, s, reference, width, log)
        st.initial(reference, lake_d, reference['inflow'])
        t0 = time.time()
        check = st.year(runoff_d, lake_d)
        lakes = st.lake_nodes
        moved = np.abs(check['out'][lakes] - reference['out'][lakes])
        steady[label] = dict(
            buffering=buffering(st, reference), seconds_a_year=round(time.time() - t0, 1),
            largest_change_of_a_lake_outflow_m3s=round(float(moved.max()), 3),
            outflows_m3s={k: [round(reference['outflow'][k], 1), round(v, 1)] for k, v in check['outflow'].items()})
        st.initial(reference, lake_d, reference['inflow'])
        stores[label] = st
        log(f'storage at width {width:g} on the design climate: {steady[label]}')
    path = WORK / 'years.npz'
    flows = {k: np.zeros((len(years), track.size), np.float32) for k in ('none', 'steady', 'narrow', 'wide')}
    summaries = {k: [] for k in flows}
    done = 0
    if resume and path.exists():
        with np.load(path, allow_pickle=False) as z:
            meta = json.loads(str(z['metadata']))
            done = meta['done']
            summaries = meta['summaries']
            for k in flows:
                flows[k] = z[f'flow_{k}']
            for k, st in stores.items():
                st.vol = z[f'vol_{k}']
                st.sea_vol = meta['seas'][k]['sea_vol']
                st.main_vol = meta['seas'][k]['main_vol']
        log(f'resuming after {done} years')
    clim_sum, land_runoff = None, []
    for i, y in enumerate(years):
        t0 = time.time()
        clim = year_climate(sums[y], applied, delta, rotation)
        fields = ('pr_mm_day', 'evap_mm_day', 'lsm')
        clim_sum = ({k: clim[k].astype(np.float64) for k in fields} if clim_sum is None
                    else {k: clim_sum[k] + clim[k] for k in fields})
        runoff_y, lake_y, _, _ = loaded(clim, m['evaporation'], rows, cols)
        land_runoff.append(float((np.maximum(runoff_y.ravel().astype(np.float64), 0.0) * MM_DAY * net.area)[land].sum()))
        if i < done:
            continue
        for k, shut in (('none', closed), ('steady', None)):
            r = without_storage(net, s, runoff_y, lake_y, reference, shut=shut)
            flows[k][i] = r['q'][track]
            summaries[k].append(dict(year=y, **year_summary(r['supply'], r['loss'], r['outflow'], r['q'], into)))
            del r
        for k, st in stores.items():
            res = st.year(runoff_y, lake_y)
            q = st.field(runoff_y, lake_y, res)
            flows[k][i] = q[track]
            summaries[k].append(dict(year=y, **year_summary(res['supply'], res['loss'], res['outflow'], q, into)))
            del q
        checkpoint(path, i + 1, flows, summaries, stores)
        log(f'year {y}: {time.time() - t0:.0f} s; without storage Moscoviense {summaries["none"][-1]["outflow_m3s"]}')
    mean = {k: (v / len(years)).astype(np.float32) for k, v in clim_sum.items()}
    mean.update(lat=clim['lat'], lon=clim['lon'])
    return flows, summaries, steady, mean, (applied, delta), land_runoff


def main_part_balance(reference):
    """The main part of the near-side sea on the design climate: what the design routing brings it and loses at its
    level, against the neck study's steady state, whose rivers were counted with the near-side sea whole."""
    iN = WATERS.index('near_side_main')
    neck = json.loads(N.OUT.read_text())['budgets'][BUDGET]
    balance = reference['main_supply'] - reference['loss'][iN]
    return dict(routing_brings_m3s=round(reference['main_supply']), neck_study_brings_m3s=neck['supply_m3s']['main'],
                loss_at_its_level_m3s=round(float(reference['loss'][iN])), routing_balance_m3s=round(balance),
                design_spill_m3s=round(reference['outflow']['neck']),
                offset_held_m3s=round(balance - reference['outflow']['neck'], 1),
                note=('the design routing routes the rivers again over the floor the parted sea leaves dry, and brings the '
                      'main part more than the neck study counted with the sea whole; both years\' procedures keep the '
                      'design\'s spill on the design climate and move it by each year\'s departure'))


def threshold_scale(net, land, runoff_mm, frame):
    """What the threshold means on this surface: the land's mean runoff on the design climate, the catchment that
    yields THRESHOLD_M3S at it, and the climate model's land cells' area, over which each year's rain is resolved."""
    from climate.gcm import boundary as bd
    land_km2 = float(net.area[land].sum()) / 1e6
    runoff = float((np.maximum(np.asarray(runoff_mm, np.float64).ravel(), 0.0) * MM_DAY * net.area)[land].sum())
    t21 = bd.read_t21(TS_PRODUCTS / frame['surface'])
    cell = t21['area'] / t21['area'].sum() * 4.0 * np.pi * (MOON_RADIUS / 1e3) ** 2
    per_km2 = runoff / land_km2
    return dict(land_km2=round(land_km2), land_runoff_m3s=round(runoff),
                runoff_mm_a_year=round(per_km2 / 1e6 * YEAR_S * 1000.0, 1),
                catchment_yielding_the_threshold_km2=round(THRESHOLD_M3S / per_km2),
                climate_land_cell_km2_median=round(float(np.median(cell[t21['mask'] > 0]))))


def ratio_spread(a, b):
    """The spread of a over b where b is positive: median and 5th and 95th percentiles."""
    ok = b > 0
    r = a[ok] / b[ok]
    return dict(median=round(float(np.median(r)), 4), p5=round(float(np.percentile(r, 5)), 4),
                p95=round(float(np.percentile(r, 95)), 4), cells=int(ok.sum()))


def buffering(st, reference):
    """How long each water's storage holds its outflow: the e-folding time A / (dQ/dh) at the design state, in years,
    for the upper seas, the main part and the lakes passing on THRESHOLD_M3S or more."""
    out = {}
    for k, q in st.sea_q.items():
        head = (q / st.sea_k[k]) ** (2.0 / 3.0)
        out[k] = dict(head_m=round(head, 2), years=round(st.sea_area[st.iw[k]] * head / (1.5 * q) / YEAR_S, 2))
    e = float(np.interp(reference['outflow']['neck'], st.neck_q, st.neck_head))
    slope = float(np.interp(e + 0.05, st.neck_head, st.neck_q) - np.interp(e - 0.05, st.neck_head, st.neck_q)) / 0.1
    out['near_side_main'] = dict(head_m=round(e, 2), years=round(st.sea_area[st.iw['near_side_main']] / slope / YEAR_S, 2))
    q = reference['out'][st.lake_nodes]
    big = q >= THRESHOLD_M3S
    head = (q[big] / st.k_out[big]) ** (2.0 / 3.0)
    tau = st.area_head[big] * head / (1.5 * q[big]) / YEAR_S
    nodes = st.lake_nodes[big]
    rotation = np.asarray(st.s['meta']['rotation'], float)
    longest = [dict(at=NW.place(int(st.net.sinks[nodes[i] - st.net.n_out]), st.net.cols, rotation),
                    area_km2=round(float(st.area_head[big][i]) / 1e6), outflow_m3s=round(float(q[big][i])),
                    head_m=round(float(head[i]), 2), years=round(float(tau[i]), 2)) for i in np.argsort(-tau)[:10]]
    out['lakes_passing_on_the_threshold'] = dict(
        count=int(big.sum()), over_a_year=int((tau >= 1.0).sum()), over_a_tenth_of_a_year=int((tau >= 0.1).sum()),
        median_years=round(float(np.median(tau)), 4), longest=longest)
    return out


def build(log=print, resume=False):
    started = time.time()
    s = design_structure(log)
    m = s['meta']
    rotation = np.asarray(m['rotation'], float)
    net = Network(s['h'], s['water'], m['levels'], s['lat'], m['cell_km'], log)
    rows, cols = net.rows, net.cols
    routing = DRt.load(log)
    q_design = np.asarray(routing['discharge_m3s'], float)
    if not (np.array_equal(routing['receiver'], net.receiver) and np.array_equal(routing['outlet'], net.outlet)):
        raise SystemExit('the network is not the design routing\'s')
    wet = routing['routing_lake'] & ~net.outlet
    into = np.flatnonzero((net.water[net.receiver] >= 0) & (net.water < 0) & ~wet)
    runoff_d, lake_d, _, _ = FD.load(m['evaporation'])
    runoff, lake = net.sources(runoff_d, lake_d, ops_list(s, s['op_design']))
    q_exact, _ = net.replay(runoff, lake)
    reference = without_storage(net, s, runoff_d, lake_d)
    runoff, lake = net.sources(runoff_d, lake_d, ops_of(s, reference['outflow']))
    reference['inflow'] = net.cascade(runoff, lake)[0]
    del runoff, lake
    closed = (reference['out'] <= 0.0) & (np.arange(net.n_nodes) >= net.n_out)
    land = ~wet & ~net.outlet
    major_d = land & (q_design >= THRESHOLD_M3S)
    design_check = dict(
        replay_largest_difference_m3s=float(np.abs(q_exact - q_design).max()),
        years_procedure_largest_difference_m3s=round(float(np.abs(reference['q'] - q_design).max()), 3),
        years_procedure_outflows_m3s={k: [round(v, 1) for v in (dict(moscoviense=m['surplus'], **m['over'],
                                                                     neck=m['neck_flow'])[k], reference['outflow'][k])]
                                      for k in reference['outflow']},
        rivers_into_the_seas_m3s=round(float(q_design[into].sum())),
        main_part=main_part_balance(reference))
    del q_exact
    log(f'design check: {design_check}, {time.time() - started:.0f} s')
    track = np.flatnonzero(land & (q_design >= TRACK_M3S))
    flows, summaries, steady, clim_mean, (applied, delta), land_runoff = run_years(
        net, s, reference, (runoff_d, lake_d), closed, track, into, land, log, resume)
    log(f'the years routed, {time.time() - started:.0f} s')

    # the years' mean against the design: the climate's part (the years' mean climate routed) and the routing's
    runoff_m, lake_m, _, _ = loaded(clim_mean, m['evaporation'], rows, cols)
    mean_routed = without_storage(net, s, runoff_m, lake_m, reference, shut=closed)
    del lake_m
    at = np.flatnonzero(major_d[track])
    mean_of = lambda k, f: float(np.mean([f(y) for y in summaries[k]]))                     # noqa: E731
    q_at = q_design[track][at]
    against_design = dict(
        major_cells=int(at.size),
        years_mean_over_design=ratio_spread(flows['none'].mean(axis=0)[at], q_at),
        mean_climate_routed_over_design=ratio_spread(mean_routed['q'][track][at], q_at),
        years_mean_over_mean_climate_routed=ratio_spread(flows['none'].mean(axis=0)[at], mean_routed['q'][track][at]),
        steady_years_mean_over_design=ratio_spread(flows['steady'].mean(axis=0)[at], q_at),
        storage_years_mean_over_design={k: ratio_spread(flows[k].mean(axis=0)[at], q_at) for k in ('narrow', 'wide')},
        rivers_into_the_seas_m3s=dict(design=round(float(q_design[into].sum())),
                                      mean_climate_routed=round(float(mean_routed['q'][into].sum())),
                                      **{f'years_mean_{k}': round(mean_of(k, lambda y: y['rivers_into_the_seas_m3s']))
                                         for k in summaries}),
        outflows_m3s={o: dict(design=round(reference['outflow'][o]), mean_climate_routed=round(mean_routed['outflow'][o]),
                              **{f'years_mean_{k}': round(mean_of(k, lambda y: y['outflow_m3s'][o])) for k in summaries})
                      for o in reference['outflow']},
        land_runoff_m3s=dict(
            note=('the routing\'s runoff from the design\'s dry ground: each 1.9 km cell\'s rain less evaporation where '
                  'positive, so a year wetter where that straddles nothing adds runoff and a drier one cannot take any '
                  'away'),
            design=round(float((np.maximum(runoff_d.ravel().astype(np.float64), 0) * MM_DAY * net.area)[land].sum())),
            mean_climate=round(float((np.maximum(runoff_m.ravel().astype(np.float64), 0) * MM_DAY * net.area)[land].sum())),
            years_mean=round(float(np.mean(land_runoff))), lowest_year=round(min(land_runoff)),
            highest_year=round(max(land_runoff)), by_year=[round(v) for v in land_runoff]))
    del mean_routed, runoff_m

    # the losses on the design climate
    frame_path, frame = FD.frame_record()
    against_design['threshold_scale'] = threshold_scale(net, land, runoff_d, frame)
    lake_routing_budget = FD.load(1.0)[1]
    loss_record, loss_fields = dryland_losses(net, s, q_design, routing['lake_cell'],
                                              {'routing_budget': lake_routing_budget, BUDGET: lake_d}, frame, log)
    del lake_routing_budget
    after = {}
    for end, field in loss_fields.items():
        r = without_storage(net, s, runoff_d, lake_d, reference, river_loss=field, shut=closed)
        after[end] = dict(q=r['q'], supply=r['supply'], outflow=r['outflow'],
                          taken_m3s=round(float(field.sum()), 1), into=round(float(r['q'][into].sum()), 1))
        del r
    loss_record['routed'] = dict(
        rule=('each loss taken from the river cell it falls on (a floodplain\'s from the nearest river cell), never more '
              'than the cell carries, and the year routed again: a lake or sea below receives less, Moscoviense passes '
              'on less by its fixed shares and the neck less by the main part\'s balance'),
        ends=dict(low='the narrow channel and belt (width coefficient 3, belt 0.5 km)',
                  high='the wide channel and belt (15, 2 km) and the floodplains kept wet at 5 m'),
        taken_m3s={k: v['taken_m3s'] for k, v in after.items()},
        rivers_into_the_seas_m3s=dict(none=round(float(reference['q'][into].sum()), 1),
                                      **{k: v['into'] for k, v in after.items()}),
        waters_receive_m3s={w: dict(none=round(float(reference['supply'][i]), 1),
                                    **{k: round(float(v['supply'][i]), 1) for k, v in after.items()})
                            for i, w in enumerate(WATERS)},
        outflows_m3s={o: dict(none=round(reference['outflow'][o], 1),
                              **{k: round(v['outflow'][o], 1) for k, v in after.items()}) for o in reference['outflow']})
    earlier = json.loads(RL.OUT.read_text())
    loss_record['against_the_earlier_climate'] = dict(
        note=('[river_losses.json on the present design climate at 8 degrees with the governed routing, this estimate on '
              'the design climate with the design routing]; the connected floodplains are the earlier record\'s on the '
              'routing\'s own budget and these on the closed budget, whose routing the design routing is'),
        **{label: dict(dryland_rivers_100_km=[earlier['budgets'][label]['rivers']['100']['length_km'],
                                              loss_record['budgets'][label]['rivers']['100']['length_km']],
                       channel_m3s=[[earlier['budgets'][label]['channel'][f'{c:g}']['evaporation_m3s'],
                                     loss_record['budgets'][label]['channel'][f'{c:g}']['evaporation_m3s']] for c in RL.WIDTH],
                       riparian_belt_m3s=[[earlier['budgets'][label]['riparian_belt'][f'{b:g}']['evaporation_m3s'],
                                           loss_record['budgets'][label]['riparian_belt'][f'{b:g}']['evaporation_m3s']]
                                          for b in RL.RIPARIAN_KM],
                       dry_land_share_of_land=[earlier['budgets'][label]['dry_land_share_of_land'],
                                               loss_record['budgets'][label]['dry_land_share_of_land']])
           for label in ('routing_budget', BUDGET)},
        floodplains_kept_wet_km2=[earlier['floodplains_at']['kept_wet_km2'], loss_record['floodplains_at']['kept_wet_km2']],
        floodplains_m3s=[earlier['floodplains_at']['evaporation_m3s'], loss_record['floodplains_at']['evaporation_m3s']])
    q_after = {k: v['q'][track] for k, v in after.items()}
    taken_high = np.maximum(reference['q'][track] - q_after['high'], 0.0)
    del after, loss_fields
    log(f'losses: {loss_record["routed"]["taken_m3s"]}, {time.time() - started:.0f} s')

    # the classes
    stats = {k: classify(v, THRESHOLD_M3S) for k, v in flows.items()}
    stats_large = {k: classify(v, LARGE_M3S) for k, v in flows.items()}
    lake_body, body_km2 = lake_bodies(wet, rows, cols, net.area)
    major = track[q_design[track] >= THRESHOLD_M3S]
    down = next_down(major, q_design, wet, net.outlet, net.water, net.receiver, net.sink, net.node_of, reference['out'],
                     lake_body, body_km2, THRESHOLD_M3S)
    stem_list = stems(major, q_design, down)
    slot = np.full(net.n, -1, np.int64)
    slot[track] = np.arange(track.size)
    length = RL.flow_length_km(net.receiver.astype(np.int64), rows, cols, np.asarray(s['lat'], float), m['cell_km'])
    dry = (np.asarray(lake_d).ravel() < 0.0)
    from_lake = np.zeros(net.n, bool)
    from_lake[net.receiver[wet]] = True
    from research.studies.atlas_round import falls as FA
    name = FA.namer()
    reach_rows, cell_reach = [], np.full(track.size, -1, np.int64)
    for i, stem in enumerate(stem_list):
        cells = np.asarray(stem['cells'])
        j = slot[cells]
        primary = stats['wide']['sustained'][j]
        cut = np.flatnonzero(np.diff(primary.astype(int))) + 1
        for piece in np.split(np.arange(cells.size), cut):
            c = cells[piece]
            e = int(c[-1])
            je = slot[e]
            last = piece[-1] == cells.size - 1
            ends = stem['ends'] if last else ('changes', int(cells[piece[-1] + 1]))
            row = dict(stem=i, cells=int(c.size), length_km=round(float(length[c].sum()), 1),
                       dryland_km=round(float(length[c][dry[c]].sum()), 1),
                       head=NW.place(int(c[0]), cols, rotation), end=NW.place(e, cols, rotation),
                       leaves_a_lake=bool(from_lake[c[0]] and piece[0] == 0),
                       ends=dict(kind=ends[0], **({'water': WATERS[ends[1]]} if ends[0] == 'sea' else
                                                  {'lake_km2': round(float(body_km2[ends[1]]))} if ends[0] in ('lake', 'closed')
                                                  else {})),
                       design_m3s=round(float(q_design[e]), 1),
                       sustained=bool(stats['wide']['sustained'][je]),
                       stops=bool(stats['wide']['stops'][je]))
            for k in ('none', 'steady', 'narrow', 'wide'):
                row[k] = dict(mean_m3s=round(float(stats[k]['mean'][je]), 1), lowest_m3s=round(float(stats[k]['lowest'][je]), 1),
                              second_lowest_m3s=round(float(stats[k]['second'][je]), 1),
                              years_at_the_threshold=int(stats[k]['years'][je]), sustained=bool(stats[k]['sustained'][je]))
            row['after_losses_m3s'] = dict(low=round(float(q_after['low'][je]), 1), high=round(float(q_after['high'][je]), 1),
                                           lowest_year_high=round(float(max(stats['wide']['lowest'][je] - taken_high[je], 0.0)), 1))
            cell_reach[j[piece]] = len(reach_rows)
            reach_rows.append(row)
    log(f'{len(stem_list)} rivers in {len(reach_rows)} reaches, {time.time() - started:.0f} s')
    return dict(s=s, net=net, routing=routing, track=track, flows=flows, summaries=summaries, steady=steady,
                stats=stats, stats_large=stats_large, stems=stem_list, reaches=reach_rows, cell_reach=cell_reach,
                q_design=q_design, q_after=q_after, taken_high=taken_high, design_check=design_check,
                against_design=against_design, loss_record=loss_record, frame_path=frame_path,
                t21=TS_PRODUCTS / frame['surface'], applied=applied,
                delta=delta, reference=reference, name=name, rotation=rotation, started=started)


def river_row(reach, name):
    """A river's line in the record: where it runs and ends, beside which named feature, and its flows."""
    end = reach['end']['today']
    return dict(end=reach['end'], beside=name(*end), ends=reach['ends'], length_km=reach['length_km'],
                design_m3s=reach['design_m3s'], **{k: reach[k] for k in ('none', 'narrow', 'wide', 'after_losses_m3s')},
                leaves_a_lake=reach['leaves_a_lake'])


def write(b, log=print):
    """The product, the reaches' table and the record from build()'s results."""
    s, net, track, stats, flows = b['s'], b['net'], b['track'], b['stats'], b['flows']
    m = s['meta']
    years = list(range(YEARS[0], YEARS[1] + 1))
    q_track = b['q_design'][track]
    major = q_track >= THRESHOLD_M3S
    large = q_track >= LARGE_M3S
    intermittent = major & ~stats['wide']['sustained']
    length = RL.flow_length_km(net.receiver.astype(np.int64), net.rows, net.cols, np.asarray(s['lat'], float),
                               m['cell_km'])[track]
    evidence = ('A calculation on products: the design climate\'s fifteen settled years (an estimate carried by pattern '
                'scaling to 300.31 K), each routed on the design routing\'s structure held fixed, without storage and '
                'with the lakes and the upper seas carrying their volume from year to year over outlets sized by Earth\'s '
                'channel widths; the rivers\' losses in the drylands by river_losses.py\'s estimate. Fifteen years is a '
                'short record, and no figure is validated against an observation.')
    reading_rule = ('Cells index rows from 90 N and columns from 0 E, 16 per degree, on today\'s Moon: cell = row * 5760 + '
                    'column. cell lists the river cells on dry ground carrying at least track_m3s on the design routing; '
                    'flow_m3s is each settled year\'s mean discharge (years by cells) with the lakes\' and upper seas\' '
                    'storage over wide outlets, flow_narrow_m3s over narrow ones and flow_without_storage_m3s with none. '
                    'sustained_m3s is the lowest year\'s flow with storage over wide outlets, the sustained discharge; '
                    'major marks the cells carrying threshold_m3s on the design routing and intermittent those of them '
                    'under it in at least one year with that storage. after_losses_*_m3s are the design climate\'s flows '
                    'after the drylands\' losses at the low and high ends; reach indexes the reaches\' table.')
    meta = dict(schema=PRODUCT_SCHEMA, grid=dict(pixels_per_degree=dr.PPD, rows=net.rows, cols=net.cols), years=years,
                track_m3s=TRACK_M3S, threshold_m3s=THRESHOLD_M3S, outlet_width_coefficients=dict(narrow=RL.WIDTH[0],
                                                                                                wide=RL.WIDTH[-1]),
                design_routing=dict(path=str(DRt.PRODUCT), sha256=sha256(DRt.PRODUCT)), reaches=str(REACHES),
                evidence=evidence, reading_rule=reading_rule)
    BULK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PRODUCT, metadata=json.dumps(meta), cell=track.astype(np.int32),
                        design_m3s=q_track.astype(np.float32), flow_m3s=flows['wide'], flow_narrow_m3s=flows['narrow'],
                        flow_without_storage_m3s=flows['none'], sustained_m3s=stats['wide']['lowest'].astype(np.float32),
                        sustained_without_storage_m3s=stats['none']['lowest'].astype(np.float32), major=major,
                        intermittent=intermittent, after_losses_low_m3s=b['q_after']['low'].astype(np.float32),
                        after_losses_high_m3s=b['q_after']['high'].astype(np.float32),
                        reach=b['cell_reach'].astype(np.int32))
    REACHES.write_text(json.dumps(dict(schema=SCHEMA + '#reaches', threshold_m3s=THRESHOLD_M3S, reaches=b['reaches']),
                                  indent=0) + '\n')
    log(f'wrote {PRODUCT} and {REACHES}')

    # the classes, by river (each stem as it reaches its end) and by length
    reaches = b['reaches']
    ends = [r for r in reaches if r['ends']['kind'] != 'changes']
    by_class = {}
    for k in ('none', 'steady', 'narrow', 'wide'):
        sus = stats[k]['sustained']
        by_class[k] = dict(
            rivers_sustained_to_their_end=sum(r[k]['sustained'] for r in ends),
            rivers_intermittent_at_their_end=sum(not r[k]['sustained'] for r in ends),
            rivers_stopping_at_their_end_in_a_year=sum(r[k]['lowest_m3s'] < DRY_M3S for r in ends),
            km_sustained=round(float(length[major & sus].sum())), km_intermittent=round(float(length[major & ~sus].sum())),
            km_stopping=round(float(length[major & stats[k]['stops']].sum())),
            large_km_sustained_at_their_threshold=round(float(length[large & b['stats_large'][k]['sustained']].sum())),
            large_km_under_it_in_a_year=round(float(length[large & ~b['stats_large'][k]['sustained']].sum())))
    for k in ('narrow', 'wide'):                     # where the storage's refilling leaves less than none
        short = major & (stats[k]['lowest'] < stats['none']['lowest'] - 0.5)
        by_class[k]['km_whose_lowest_year_falls_below_none'] = round(float(length[short].sum()))
        by_class[k]['largest_shortfall_against_none_m3s'] = round(float(
            (stats['none']['lowest'] - stats[k]['lowest'])[major].max()), 1)
    leaving = [r for r in ends if any(x['leaves_a_lake'] for x in reaches if x['stem'] == r['stem'])]
    name = b['name']
    largest = sorted((r for r in ends if r['wide']['sustained']), key=lambda r: -r['design_m3s'])[:15]
    weakest = sorted((r for r in ends if not r['wide']['sustained']), key=lambda r: -r['design_m3s'])[:25]
    held = sorted((r for r in ends if r['wide']['sustained'] and not r['none']['sustained']), key=lambda r: -r['design_m3s'])
    into_seas = sorted((r for r in ends if r['ends']['kind'] == 'sea'), key=lambda r: -r['design_m3s'])[:15]
    bands = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))
    by_band = {f'{a}-{b}': {k: sum(1 for r in ends if a <= abs(r['end']['turned'][0]) < b and not r[k]['sustained'])
                            for k in ('none', 'wide', 'narrow')} | dict(rivers=sum(1 for r in ends if a <= abs(r['end']['turned'][0]) < b))
               for a, b in bands}
    summaries = b['summaries']
    scale = b['against_design'].pop('threshold_scale')
    record = dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (DRt.OUT, FD.OUT, N.OUT, RL.OUT, ROOT / FD.ESTIMATE, b['frame_path'])},
        ignored_inputs={**{str(DRt.PRODUCT): sha256(DRt.PRODUCT), FD.CLIMATOLOGY: sha256(ROOT / FD.CLIMATOLOGY),
                           str(STRUCTURE): sha256(STRUCTURE), rel(b['t21']): sha256(b['t21'])},
                        **{str(f): sha256(f) for f in year_files()}},
        evidence=evidence,
        reading_rule=('Flows in m3/s, lengths in km, areas in km2; a year is a settled year of the design climate (model '
                      'years 45-59). Each year is routed four ways: none, with no water carried between years and the '
                      'hollows the design leaves below their passes held there, the lower bound on a dry year\'s flow '
                      'for the waters the design holds full; steady, each year a steady state of its own, those hollows '
                      'filling and spilling in a wet year, shown to measure what holding them takes away; narrow and wide, '
                      'with the lakes\' and the upper seas\' storage carried from year to year over outlets of width 3 '
                      'and 15 m per (m3/s)^0.5 widened for lunar gravity, wide holding the least. A reach is part of a '
                      'major river (threshold_m3s or more on the design routing) between confluences, lakes of '
                      'lake_break_km2 or more and changes of class, its flows read at its downstream end; it is '
                      'sustained when it carries the threshold in every year with the wide outlets\' storage and '
                      'intermittent otherwise. A river is a stem from its head to where it ends or joins a larger river, '
                      'classed by its end. after_losses gives the design climate\'s flow after the drylands\' losses at '
                      'the low and high ends, and lowest_year_high the wide storage\'s lowest year less the high end\'s '
                      'loss of the mean year. product holds the cells, reaches the full table.'),
        limits=[
            'The climate is an estimate: the run at sunlight 0.936 carried by pattern scaling from 301.00 K to 300.31 K, '
            'each year keeping its own departures from the run\'s mean; a run at the design temperature would differ by '
            'about the noise of fifteen years\' means.',
            'Fifteen years is a short record: the lowest year is about the driest year in fifteen and the second lowest '
            'the driest in seven or eight; longer records hold drier years, so a river sustained here can fall under the '
            'threshold in a rarer drought. The years follow one another in the model\'s order.',
            'Without storage each year passes on only its own surplus: the lower bound on a dry year\'s flow for the '
            'waters the design holds full. With storage a lake that fell below its pass in a dry year refills before it '
            'spills again, so in the years after a drought its river can carry less than without storage '
            '(by_storage.*.km_whose_lowest_year_falls_below_none).',
            'The storage\'s outlets are not designed or measured: each passes critical flow over a section as wide as '
            'Earth\'s alluvial channels at its design flow, at both ends of their range, so the head a lake holds over '
            'its pass, and with it the water it carries into a dry year, is an assumption; a pass the water has not cut '
            'would be wider and hold less. A lake\'s area over its pass is held at its area at the pass.',
            f'Each year is a year\'s mean rain and evaporation in {STEPS} steps: the seasons and the storms within a '
            f'year are left out, and a year\'s rain is resolved over the climate model\'s land cells of about '
            f'{scale["climate_land_cell_km2_median"]:,} km2, the routing\'s 1.9 km cells taking it bilinearly.',
            'The losses are river_losses.py\'s estimate at the ends of its ranges, taken from the mean year\'s flow; the '
            'lowest year after losses deducts the mean year\'s loss.',
            'The routing is the design routing\'s steady fill and spill on 1.9 km cells, with its own approximations '
            '(each catchment balanced alone, as route_levels balances it). Nothing is validated against an '
            'observation.'],
        settings=dict(threshold_m3s=THRESHOLD_M3S, large_m3s=LARGE_M3S, track_m3s=TRACK_M3S, dry_m3s=DRY_M3S,
                      lake_break_km2=LAKE_BREAK_KM2, steps_a_year=STEPS, outlet_width_coefficients=list(RL.WIDTH),
                      lunar_widening=round(RL.LUNAR_WIDENING, 3), outlet_sized_at_least_m3s=Q_FLOOR,
                      lake_table_depths_m=list(DEPTHS_M)),
        climate=dict(run=RUN, years=list(YEARS), estimate=FD.LABEL, delta_K=round(b['delta'], 4),
                     responses_per_K={k: {d: round(b['applied'][k][d]['per_K'], 6) for d in ('land', 'water')}
                                      for k in ('pr', 'evap')},
                     open_water_evaporation_factor=m['evaporation'], frame=m['frame']),
        design_check=b['design_check'], against_the_design=b['against_design'],
        storage=b['steady'],
        years={k: v for k, v in summaries.items()},
        losses=b['loss_record'],
        classes=dict(threshold_m3s=THRESHOLD_M3S, threshold_scale=scale,
                     rivers=len(ends), reaches=len(reaches),
                     major_km=round(float(length[major].sum())), large_km=round(float(length[large].sum())),
                     by_storage=by_class, rivers_leaving_a_lake=len(leaving),
                     rivers_leaving_a_lake_sustained=dict(none=sum(r['none']['sustained'] for r in leaving),
                                                          wide=sum(r['wide']['sustained'] for r in leaving)),
                     intermittent_by_turned_latitude=by_band,
                     largest_into_the_seas=[river_row(r, name) for r in into_seas],
                     largest_sustained=[river_row(r, name) for r in largest],
                     intermittent_largest=[river_row(r, name) for r in weakest],
                     sustained_only_with_storage=[river_row(r, name) for r in held[:25]]),
        product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)),
        reaches=dict(path=str(REACHES), sha256=sha256(REACHES)),
        seconds=round(time.time() - b['started']))
    return record


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    p.add_argument('--resume', action='store_true', help='take the years already routed from the research drive')
    args = p.parse_args(argv)
    log = lambda *a: print(*a, flush=True)                                     # noqa: E731
    record = write(build(log, resume=args.resume), log)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1) + '\n')
    log(f'wrote {args.out} in {record["seconds"]} s: {json.dumps(record["classes"]["by_storage"])}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
