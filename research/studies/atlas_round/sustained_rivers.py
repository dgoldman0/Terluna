"""The design surface's sustained rivers: the design climate's fifteen settled years routed one at a time on the design
routing's network, without storage and with the lakes and the upper seas carrying water from one year to the next; the
rivers' losses in the drylands on the design climate; and the major rivers classed as sustained or intermittent.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.sustained_rivers   # about 10 minutes

--resume takes the years already routed from the research drive's checkpoint, where it was kept on this design routing.

The design routing (design_routing.py) draws every river at the long-term mean flow of the design climate's closed
budget. This module asks which rivers keep flowing through the climate's dry years, in four parts:

- the years without storage. Each settled year of the design climate (the run on the design surface at sunlight 0.936,
  years 45-59) is read from the run's year sums (moisture_map.cached_year in drylands.py's cache), moved by the
  estimate's two steps as moisture_final.adjusted moves them (to the design temperature, then onto the design's own
  surface), carried onto today's ground with the settled frame's rotation as final_drainage.carry carries the mean,
  and loaded as final_drainage.load loads it, with the closed budget's factor on the open water's evaporation. It is
  routed on the design routing's structure held fixed: the floor, the seas at their design levels, the lakes as the
  corrected rule fills them in the design (lake_rule.py: each hollow as one water, the courses joined through the
  lakes by hollow_courses.py) and the works where the design puts them. Each water passes on that year's own
  surplus: a lake standing full or an upper sea passes on what the year brings it beyond its loss, a lake the design
  leaves below its pass passes nothing, and Moscoviense's fixed crests divide its overflow in the design's fixed
  shares. Where the near-side sea stands parted, its main part sends over its neck the design's spill moved by the
  year's departure in its balance; where the design holds it as one sea across the neck, as it does since
  18 November, the flow between its parts passes through its water and both are held at their levels;
- the years with storage. Every lake and the three upper seas carry their volume from one year to the next, starting
  from the design state, in STEPS steps a year. A lake is a body of hollows standing as one water, whose volume is
  held at the lowest point it leaves by. Above its pass a water passes on what critical flow passes over an outlet as
  wide as Earth's alluvial channels at the water's design flow (river_losses.WIDTH, widened for lunar gravity, at both
  ends of the range); below it, it passes nothing and evaporates from the area its level covers, so a lake the design
  leaves below its pass can rise in wet years and pass water on;
- the losses. river_losses.py's estimate (the channels, the plants along them and the floodplains, in the drylands) is
  run again on the design routing and the design climate, and each loss is taken from the river where it falls and
  carried down to the seas;
- the classes. Every reach of the rivers carrying at least THRESHOLD_M3S on the design routing is sustained when it
  carries that much in every year with the wide outlets' storage, and intermittent otherwise, with its mean, lowest-year
  and second-lowest-year flows in each way of routing the years, the years it carries that much and its flow after the
  losses.

The design routing's product keeps the joined receivers and the flows; the years need the hollows as well, so the
routing is made again in memory from the same records (design_state) and must give the product back exactly. Until
18 November the network was the earlier lake rule's, each catchment balanced alone, and each year was also routed as a
steady state of its own with every hollow filling and spilling; on the corrected rule such a state moves the lakes'
shores, so its rivers would not lie on the design's cells, and the storage's years carry what it showed.

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
from research.studies.atlas_round import river_losses as RL
from research.studies.drainage_network import network as NW
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
        'research/studies/atlas_round/design_works.py', 'research/studies/atlas_round/near_side_neck.py',
        'research/studies/atlas_round/lake_rule.py', 'research/studies/atlas_round/lake_rule_balance.py',
        'research/studies/atlas_round/hollow_courses.py', 'research/studies/atlas_round/moisture_final.py',
        'research/studies/atlas_round/river_losses.py', 'research/studies/drainage_network/network.py',
        'geography/drainage.py', 'climate/gcm/turned_climatology.py', 'climate/gcm/pattern_scaling.py',
        'biosphere/ecology/moisture_map.py')
# The structure of the design on the earlier lake rule (work/design_structure.npz) was kept under these files'
# hashes; lake_rule.py, lake_rule_balance.py and design_works.py read it and its key. This module no longer builds it.
STRUCTURE_CODE = ('research/studies/atlas_round/design_routing.py', 'research/studies/atlas_round/final_drainage.py',
                  'research/studies/atlas_round/governed_drainage.py', 'research/studies/atlas_round/redrawn_drainage.py',
                  'research/studies/atlas_round/near_side_neck.py', 'research/studies/atlas_round/moisture_final.py',
                  'research/studies/drainage_network/network.py', 'research/studies/drainage_network/routing.py',
                  'geography/drainage.py')
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


# ---------------------------------------------------------------- the design routing made again

def design_state(log=print):
    """The design routing made again in memory, as design_routing.route makes it, so that the years can be routed on
    its hollows: the product keeps the joined receivers and the flows, and the years need the hollows' tree and which
    hollows stand full as well. The weirs are settled again from final_drainage.py's shares (lake_rule_balance.govern,
    the routing's own reckoning), the near-side sea's parts are taken at the levels the product records, and the
    hollows are settled once on them; the result must be the product's receivers and flows. Returns `s`, the structure
    the years' procedures read (the floor, the water over each cell as an index into WATERS with each water's level,
    and the works' cells with the kind, the sign and the design flow of each and the part of its branch's share each
    cut carries), the routing's pieces (the hollows' structure, their state, the routing and the joined courses), the
    product and the design climate's runoff and open-water balance."""
    from research.studies.atlas_round import design_works as DW
    from research.studies.atlas_round import hollow_courses as HC
    from research.studies.atlas_round import lake_rule as LR
    from research.studies.atlas_round import lake_rule_balance as LB
    started = time.time()
    if tuple(LB.WATERS) != WATERS:
        raise SystemExit('lake_rule_balance\'s waters are not in WATERS\' order')
    routing = DRt.load(log)
    made = json.loads(DRt.OUT.read_text())['summary']
    if made['near_side']['spills_over_a_sill_m3s']:
        raise SystemExit('a part of the near-side sea spills over a sill into a basin, which this module does not carry')
    case = json.loads(FD.OUT.read_text())['cases'][BUDGET]
    frame_path, record = FD.frame_record()
    d, st = DW.load(log)
    runoff, lake_net, _, _ = FD.load(case['evaporation'])
    begun = dict(share=case['weir_history'][-1]['shares_m3s'], over=case['full_seas_pass_on_m3s'])
    g = LB.govern(d, st, runoff, lake_net, start=begun, log=log)
    gov = g['gov']
    meta_r = routing['meta']
    neck = meta_r['neck']
    joined = bool(neck['joined'])
    water, levels = LB.parted_waters(d, st, meta_r['levels_m']['near_side_main'], meta_r['levels_m']['near_side_south'])
    hollows = LB.structure_of(d, water, levels)
    terms = LR.hollow_terms(hollows, runoff, lake_net)
    over_the_neck = [] if neck['spill_cell'] is None else [(int(neck['spill_cell']), float(neck['flow_m3s']))]
    state, flows, sources = LB.branch_state(hollows, terms, g['fixed'] + over_the_neck, g['branches'],
                                            LB.roots_of(hollows['tree']))
    r = LR.hollow_routing(hollows, terms, state)
    courses = HC.connect(hollows, terms, state, r, d['h'].shape)
    if not np.array_equal(courses['receiver'], routing['receiver']):
        raise SystemExit('the routing made again does not join its courses as the design routing\'s product does')
    apart = float(np.abs(courses['discharge'] - routing['discharge_m3s']).max())
    if apart > 1e-6:
        raise SystemExit(f'the routing made again differs from the design routing\'s product by {apart:.3g} m3/s')
    # the works in the product's order, each with its design flow and the part of its branch's share it carries
    share, over, surplus = gov['share'], gov['over'], float(gov['surplus'])
    unit = dict(salt_sea_weir=share[GD.SALT_SEA], humboldtianum_branch=share[GD.HELD_BRANCH],
                smythii_marginis_branch=share[GD.REST_BRANCH], **{k: over[k] for k in GD.FULL_SEAS})
    along = dict(humboldtianum_branch=iter(flows[0]), smythii_marginis_branch=iter(flows[1]))
    ops = []
    for (cell, sign, _), kind in zip(d['ops'], d['op_kind']):
        q = float(next(along[kind])) if kind in along else float(unit[kind])
        ops.append((int(cell), kind, float(sign), q, q / unit[kind] if unit[kind] > 0.0 else 0.0))
    ops += [(cell, 'neck', 1.0, q, 1.0) for cell, q in over_the_neck]
    if sorted((c, sign * q) for c, _, sign, q, _ in ops) != sorted((int(c), float(q)) for c, q in sources):
        raise SystemExit('the works\' flows are not the ones the hollows were settled on')
    meta = dict(level0=float(d['level0']), cell_km=float(d['cell_km']), levels=[float(v) for v in levels],
                rotation=np.asarray(d['rotation'], float).tolist(),
                frame=dict(record=rel(frame_path), sha256=sha256(frame_path)), evaporation=float(case['evaporation']),
                surplus=surplus,
                fractions={k: float(share[k]) / surplus for k in (GD.SALT_SEA, GD.HELD_BRANCH, GD.REST_BRANCH)},
                over={k: float(v) for k, v in over.items()}, neck_flow=float(neck['flow_m3s']), joined=joined,
                made_again_largest_difference_m3s=apart)
    s = dict(meta=meta, h=d['h'], lat=d['lat'], water=water,
             op_cell=np.array([c for c, _, _, _, _ in ops], np.int64),
             op_kind=np.array([OPS.index(k) for _, k, _, _, _ in ops], np.int8),
             op_sign=np.array([x for _, _, x, _, _ in ops], float),
             op_design=np.array([x * q for _, _, x, q, _ in ops], float),
             op_part=np.array([p for _, _, _, _, p in ops], float))
    log(f'the design routing made again ({apart:.3g} m3/s from its product at most), {time.time() - started:.0f} s')
    return dict(s=s, hollows=hollows, state=state, r=r, courses=courses, routing=routing, runoff=runoff,
                lake_net=lake_net)


# ---------------------------------------------------------------- the routing on the fixed structure

def bodies_of(tree, full):
    """The hollow whose water each hollow is part of (lake_rule.hollow_routing's and hollow_courses.connect's): two
    hollows that both stand full to the sill between them are one water with the larger hollow they make, and so on up
    while that one stands full beside a full one."""
    n_all = tree['parent'].size
    top = np.where(full & (tree['parent'] >= 0) & full[np.maximum(tree['sibling'], 0)], tree['parent'], np.arange(n_all))
    while True:
        nxt = top[top]
        if np.array_equal(nxt, top):
            return top
        top = nxt


class Network:
    """The design routing's network held fixed, and its routing on any year's water. The graph is the corrected rule's
    (lake_rule.py) with the courses joined through the lakes (hollow_courses.py): each cell's steepest-descent
    neighbour, each lowest point sending its water on to the cell the design routing sends it to, and the lakes as the
    design fills them. A lake is a body, the hollows standing full to the sills between them as one water; it leaves
    by one lowest point, where its whole balance is taken, and its other lowest points pass on what reaches them
    (`through`). The bodies the design leaves below their passes pass nothing (`closed`). replay() gathers a year's
    water over the cells as hollow_courses.connect does; cascade() does the same over the catchments alone, which is
    all the seas' balances need."""

    def __init__(self, hollows, state, routing, courses, water, n_waters, shape, log=print):
        started = time.time()
        rows, cols = shape
        n = rows * cols
        tree, n_out, n_nodes = hollows['tree'], hollows['n_out'], hollows['n_nodes']
        outlet, sink, sinks, label = hollows['outlet'], hollows['sink'], hollows['sinks'], hollows['label']
        self.rows, self.cols, self.n = rows, cols, n
        self.water, self.n_waters = water, n_waters
        self.outlet, self.sink, self.sinks, self.label = outlet, sink, sinks, label
        self.n_out, self.n_nodes = n_out, n_nodes
        self.node_of = np.where(outlet | sink, label, -1).astype(np.int32)
        self.area, self.h, self.d8 = hollows['area'], hollows['hf'], hollows['receiver']
        final = np.asarray(courses['receiver'], np.int64)
        self.receiver = final
        # the bodies: the hollow each lowest point's lake is part of, and the lowest point each body leaves by, the
        # one whose water goes on to a cell outside the body
        full = np.asarray(state['full'], bool)
        top = bodies_of(tree, full)
        leaves = np.arange(n_out, n_nodes)
        hollow = top[leaves]
        to = label[final[sinks]].astype(np.int64)
        leaves_body = (to < n_out) | (top[to] != hollow)
        exit_of = np.full(tree['parent'].size, -1, np.int64)
        exit_of[hollow[leaves_body]] = leaves[leaves_body]
        if np.unique(hollow[leaves_body]).size != int(leaves_body.sum()) or (exit_of[hollow] < 0).any():
            raise RuntimeError('a lake does not leave by one lowest point')
        self.through = np.zeros(n_nodes, bool)
        self.through[n_out:] = ~leaves_body
        self.inner = self.through
        self.closed = np.zeros(n_nodes, bool)
        self.closed[n_out:] = leaves_body & ~full[hollow]
        self.leaving = np.zeros(n_nodes, bool)
        self.leaving[n_out:] = leaves_body & full[hollow]
        self.body = np.arange(n_nodes, dtype=np.int64)
        self.body[n_out:] = exit_of[hollow]
        self.level = np.full(n_nodes, np.nan)                # the pass a body overflows at
        self.level[n_out:] = tree['spill'][hollow]
        self.design_out = np.zeros(n_nodes)                  # what the design's state sends on from each lake
        self.design_out[n_out:] = np.where(self.leaving[n_out:], np.where(tree['parent'][hollow] < 0, state['out'][hollow],
                                                                          state['over'][hollow]), 0.0)
        # the lakes as the design fills them, and every cell a body's water covers when it stands full to its pass
        lake = np.asarray(routing['lake'], bool) & ~outlet
        self.lake_full = lake
        self.lake_cells = np.flatnonzero(lake)
        self.lake_node = exit_of[routing['body'][self.lake_cells]]
        if (self.lake_node < 0).any():
            raise RuntimeError('a lake cell lies in no lake')
        node = exit_of[top[hollows['owner']]]
        self.cover_cells, self.cover_node = hollows['under'][node >= 0], node[node >= 0]
        # each catchment's parent (where its water goes next), root (the outlet it reaches) and depth below the root
        parent = np.arange(n_nodes, dtype=np.int64)
        parent[n_out:] = to
        self.parent = parent
        root, depth = parent.copy(), (parent != np.arange(n_nodes)).astype(np.int64)
        while True:                                  # pointer jumping: depth counts the hops to the current pointer
            nxt = root[root]
            if np.array_equal(nxt, root):
                break
            depth = depth + depth[root]
            root = nxt
        self.root, self.depth = root, depth
        # the water each cell's own path ends in (seas), and the water a catchment's lake reaches when it goes on to one
        terminal = hollows['terminal']
        self.term_sea = np.where(outlet[terminal], water[terminal], -1).astype(np.int8)
        self.direct = np.flatnonzero((self.term_sea >= 0) & ~outlet)
        self.spill_sea = np.full(n_nodes, -1, np.int8)
        self.into_seas = leaves[to < n_out]
        self.spill_sea[self.into_seas] = self.term_sea[final[sinks[self.into_seas - n_out]]]
        # hollow_courses.connect's generations over the cells
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
        if self.order.size != int(donor.sum()):
            raise RuntimeError('the joined receivers do not all lead to a sea')
        self.bounds = np.asarray(bounds, np.int64)
        self.order_to = final[self.order]
        self.order_sink = sink[self.order]
        self.has_sink = np.add.reduceat(self.order_sink.astype(np.int64), self.bounds[:-1]) > 0
        del pieces, indeg
        log(f'network: {n_nodes - n_out:,} catchments in {int(leaves_body.sum()):,} lakes, {int(self.closed.sum()):,} of '
            f'them below their passes, {self.bounds.size - 1} generations of cells, {int(depth.max())} of catchments, '
            f'{time.time() - started:.0f} s')

    def generations(self, nodes):
        """The catchments among `nodes` grouped by depth, deepest first: each group's parents lie in later groups."""
        nodes = np.asarray(nodes, np.int64)
        nodes = nodes[nodes >= self.n_out]
        d = self.depth[nodes]
        order = np.argsort(-d, kind='stable')
        nodes, d = nodes[order], d[order]
        cut = np.flatnonzero(np.diff(d)) + 1
        return np.split(nodes, cut)

    def sources(self, runoff_mm, lake_net_mm):
        """A year's water on the fixed lakes, from runoff and open-water balance fields (mm/day), taken as
        lake_rule.hollow_terms takes them: the runoff that starts on each dry cell (run, m3/s), each lake's rain less
        evaporation at the lowest point it leaves by (lake), what each catchment's dry ground gathers (local) and what
        the dry ground draining straight to each water yields (direct)."""
        run = np.maximum(np.asarray(runoff_mm).ravel(), 0.0) * MM_DAY * self.area
        run[self.outlet | self.lake_full] = 0.0
        net = np.asarray(lake_net_mm).ravel()[self.lake_cells] * MM_DAY * self.area[self.lake_cells]
        lake = np.bincount(self.lake_node, weights=net, minlength=self.n_nodes)
        local = np.bincount(self.label, weights=run, minlength=self.n_nodes)
        direct = np.bincount(self.term_sea[self.direct], weights=run[self.direct], minlength=self.n_waters)
        return dict(run=run, lake=lake, local=local, direct=direct)

    def lands(self, cells):
        """Where a flow put in at each cell goes: the catchment whose lowest point it runs down to, or the water it
        reaches straight; -1 for the other."""
        cells = np.asarray(cells, np.int64)
        k = self.label[cells].astype(np.int64)
        return np.where(k >= self.n_out, k, -1), np.where(k < self.n_out, self.term_sea[cells], -1)

    def passed(self, nodes, inflow, lake, passed_on):
        """What the catchments `nodes` pass on with `inflow` reaching their lowest points: a lake's other lowest
        points all that reaches them, the one it leaves by its inflow and balance where positive (or `passed_on`
        where given), and nothing from a lake the design leaves below its pass."""
        if passed_on is None:
            q = np.maximum(inflow + lake[nodes], 0.0)
            q[self.closed[nodes]] = 0.0
        else:
            q = passed_on[nodes]
        return np.where(self.through[nodes], inflow, q)

    def cascade(self, src, ops=(), passed_on=None, groups=None):
        """The catchments' inflows and what each passes on (m3/s), and what each water receives, with the flows `ops`
        ((cell, m3/s) pairs, negative for a take-out) put in where they land: the routing's balance on the catchments
        alone."""
        inflow, sea = src['local'].copy(), src['direct'].copy()
        if len(ops):
            node, water = self.lands([c for c, _ in ops])
            q = np.array([q for _, q in ops], float)
            np.add.at(inflow, node[node >= 0], q[node >= 0])
            np.add.at(sea, water[water >= 0], q[water >= 0])
        out = np.zeros(self.n_nodes)
        for nodes in (groups if groups is not None else self.all_groups()):
            q = self.passed(nodes, inflow[nodes], src['lake'], passed_on)
            out[nodes] = q
            np.add.at(inflow, self.parent[nodes], q)
        into = self.into_seas
        return inflow, out, sea + np.bincount(self.spill_sea[into], weights=out[into], minlength=self.n_waters)

    def all_groups(self):
        if not hasattr(self, '_groups'):
            self._groups = self.generations(np.arange(self.n_out, self.n_nodes))
        return self._groups

    def supply_of(self, q):
        """What each water receives in a replay's discharge: the flow into its cells."""
        if not hasattr(self, '_water_cells'):
            self._water_cells = [np.flatnonzero(self.water == i) for i in range(self.n_waters)]
        return np.array([float(q[c].sum()) for c in self._water_cells])

    def replay(self, src, ops=(), passed_on=None, loss=None):
        """Each cell's discharge (m3/s): hollow_courses.connect's gathering over the kept generations with the flows
        `ops` put in at their cells, each lowest point passing on what passed() gives it; `loss` (m3/s per cell) is
        taken from what a cell passes on, never more than it carries. Under a lake the discharge is what has come in
        on the way to its exit, and at the lowest point a spilling lake leaves by, what the lake passes on. Returns
        the discharge and what each catchment passed on."""
        acc = src['run'].copy()
        put = {}
        for cell, q in ops:
            put[cell] = put.get(cell, 0.0) + q
        for cell, q in put.items():
            acc[cell] += q
        lake = src['lake']
        out_node = np.zeros(self.n_nodes)
        order, to, at_sink, b = self.order, self.order_to, self.order_sink, self.bounds
        for g in range(b.size - 1):
            i0, i1 = b[g], b[g + 1]
            cells = order[i0:i1]
            out = acc[cells]
            if loss is not None:
                out = out - np.minimum(loss[cells], np.maximum(out, 0.0))
            if self.has_sink[g]:
                s = at_sink[i0:i1]
                pits = cells[s]
                k = self.node_of[pits]
                q = self.passed(k, acc[pits], lake, passed_on)
                out[s] = q
                out_node[k] = q
                leaving = self.leaving[k]
                acc[pits[leaving]] = q[leaving]
            np.add.at(acc, to[i0:i1], out)
        return acc, out_node


# ---------------------------------------------------------------- the years' climate

def year_sums():
    """The settled years' sums (moisture_map.read_year's, from drylands.py's cache), the estimate's step to the design
    temperature with the moisture map's responses per kelvin (moisture_final.estimate_responses), and its step onto
    the design's own surface as each cell's ratio of rain and of evaporation (final_drainage.surface_ratios)."""
    bm.RUNS = DL.RUNS
    bm.CACHE = str(DL.CACHE)
    record, _, applied = MF.estimate_responses()
    years = {y: bm.cached_year(RUN, y, True) for y in range(YEARS[0], YEARS[1] + 1)}
    return years, applied, float(record['estimate']['delta_K']), FD.surface_ratios()


def year_climate(year, applied, delta, rotation, ratios=None):
    """One year's rain, evaporation and land on today's latitude and longitude: the year's sums moved by the estimate's
    two steps (moisture_final.adjusted), each field the mean over the year's records as climate.gcm.climatology takes
    it (evaporation turned positive upward), carried from the settled frame as
    climate.gcm.turned_climatology.to_today carries a climatology."""
    y = MF.adjusted(year, applied, delta, ratios)
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

def op_flows(s, mosc, orientale, smythii, neck=0.0):
    """Each put-in's and take-out's flow from the outflows that feed it: Moscoviense's overflow divided by its fixed
    crests in the design's shares (register, 17 November), each cut along a branch carrying the part of the branch's
    share it carries in the design (a cut whose lake could not spare the whole share carries less), Orientale's and
    Smythii-Marginis's whole, and the neck's spill where the near-side sea stands parted."""
    f = s['meta']['fractions']
    unit = dict(salt_sea_weir=mosc * f[GD.SALT_SEA], humboldtianum_branch=mosc * f[GD.HELD_BRANCH],
                smythii_marginis_branch=mosc * f[GD.REST_BRANCH], orientale=orientale, smythii_marginis=smythii,
                neck=neck)
    return s['op_sign'] * s['op_part'] * np.array([unit[k] for k in OPS])[s['op_kind']]


def ops_list(s, flows):
    return list(zip(s['op_cell'].tolist(), np.asarray(flows, float).tolist()))


def water_losses(net, lake_net_mm):
    """Each water's loss to the air beyond the rain on it (m3/s) at its design level, as redrawn_drainage.design_state
    counts a sea's loss."""
    field = -np.asarray(lake_net_mm, np.float64).ravel() * MM_DAY * net.area
    return net.supply_of(field)


def without_storage(net, s, runoff_mm, lake_net_mm, reference=None, river_loss=None):
    """One year routed on the fixed structure with no water carried in or out: each lake the design holds full and
    each upper sea passes on what its rivers bring it beyond its loss that year, Moscoviense's by the fixed shares, and
    the lakes the design leaves below their passes pass nothing: a wet year cannot raise them to their passes within
    the year. Where the near-side sea stands parted, its main part sends over the neck the design's spill moved by its
    balance's departure from the reference year (`reference`, the design climate's own routing; the design's spill
    itself when None); where it is one sea across the neck, the flow between its parts passes through its water and
    is no river. `river_loss` (m3/s a cell) is taken from the rivers where it falls, so a lake or sea below receives
    less. Returns the discharge, what each catchment passes on, each water's supply and loss, the outflows, the main
    part's supply and the year's sources."""
    m = s['meta']
    iM, iO, iS, iN = (WATERS.index(k) for k in ('moscoviense', 'orientale', 'smythii_marginis', 'near_side_main'))
    loss = water_losses(net, lake_net_mm)
    src = net.sources(runoff_mm, lake_net_mm)

    def routed(mosc, orientale, smythii, neck, whole=False):
        ops = ops_list(s, op_flows(s, mosc, orientale, smythii, neck))
        if river_loss is None and not whole:
            return None, None, net.cascade(src, ops)[2]
        q, out = net.replay(src, ops, loss=river_loss)
        return q, out, net.supply_of(q)

    _, _, own = routed(0.0, 0.0, 0.0, 0.0)
    mosc, orientale = max(own[iM] - loss[iM], 0.0), max(own[iO] - loss[iO], 0.0)
    _, _, sup = routed(mosc, orientale, 0.0, 0.0)
    smythii = max(sup[iS] - loss[iS], 0.0)
    _, _, sup = routed(mosc, orientale, smythii, 0.0)
    outflow = dict(moscoviense=mosc, orientale=orientale, smythii_marginis=smythii)
    neck = 0.0
    if not m['joined']:
        neck = m['neck_flow']
        if reference is not None:
            neck = max(neck + (sup[iN] - reference['main_supply']) - (loss[iN] - reference['loss'][iN]), 0.0)
        outflow['neck'] = neck
    q, out, supply = routed(mosc, orientale, smythii, neck, whole=True)
    return dict(q=q, out=out, supply=supply, loss=loss, main_supply=float(sup[iN]), outflow=outflow, sources=src)


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
    """The lakes' and the upper seas' storage from step to step on the fixed network. Each lake is a body, the hollows
    standing as one water (Network.body); the lowest point it leaves by holds the body's volume against the pass it
    overflows at, and its other lowest points pass on what reaches them. Above the pass the body passes on critical
    flow over an outlet of width `width` (m per (m3/s)^0.5, river_losses.channel_width_m at the design's outflow, at
    least Q_FLOOR) at its head over its water's area at the pass; below it, it passes nothing and gains or loses the
    rain and evaporation of the cells its level covers, tabled at DEPTHS_M over all the cells its water covers when
    full, the ground it covers or uncovers ceasing or starting to yield runoff as the corrected rule has it. A lake
    the design leaves below its pass starts where its loss meets what reaches it, rises in wet years and passes water
    on once it reaches its pass. Moscoviense, Orientale and Smythii-Marginis hold their volume over
    their rims the same way, Moscoviense's outflow divided by its fixed shares. Where the near-side sea stands parted,
    its main part spills over the neck by neck_rating, its balance less the fixed offset left between the design
    routing's balance and its spill, so that the design climate keeps the design's spill; where it is one sea across
    the neck, both parts are held at their levels. The salt sea, Humboldtianum and the southern seas pass nothing on.
    Each step routes the catchments upstream first, in stages that follow the seas' order (Moscoviense and Orientale;
    the salt sea, Humboldtianum and Smythii-Marginis; the main part; the southern part), so each put-in is known
    before the catchments it feeds."""

    def __init__(self, net, s, design, width, log=print):
        started = time.time()
        self.net, self.s, self.width = net, s, width
        m = s['meta']
        n_out = net.n_out
        h = net.h.astype(np.float64)
        # the lakes' tables: area, volume against the pass, and the cells' order for the balance below it
        cells, lab = net.cover_cells, net.cover_node
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
        lake_nodes = nodes[(nodes >= n_out) & (end > start)]       # the lowest point each body leaves by
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
        self.filled = np.flatnonzero(net.lake_full[cells])       # the cells the design's lakes cover, among the tabled
        self.filled_slot = self.slot[lab[self.filled]]
        self.area_head = area_j[:, 0]                          # the head stands over the body's water at its pass
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
        self.parted = not m['joined']
        if self.parted:
            # what the design routing leaves between the main part's balance and its spill is held as a fixed offset,
            # so the design climate keeps the design's spill
            self.crest, self.neck_head, self.neck_q = neck_rating(m['level0'])
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
        # where each put-in lands: a catchment's lowest point, or a water straight
        self.op_node, self.op_sea = net.lands(s['op_cell'])
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

    def initial(self, design, runoff_design, lake_design, inflow_design):
        """The design state: each spilling lake at the head its design outflow needs, each closed lake where its loss
        meets its inflow, the upper seas at the heads their design overflows need and, where the near-side sea stands
        parted, the main part where the neck passes the design's spill."""
        q = design['out'][self.lake_nodes]
        head = np.where(q > 0.0, (q / self.k_out) ** (2.0 / 3.0), 0.0)
        vol = self.area_head * head
        closed = q <= 0.0
        if closed.any():
            ntab = self.table(runoff_design, lake_design)[closed]
            need = -inflow_design[self.lake_nodes][closed]
            vol[closed] = self._level_for(ntab, self.vtab[closed], need)
        self.vol = vol
        self.sea_vol = {k: self.sea_area[self.iw[k]] * (self.sea_q[k] / self.sea_k[k]) ** (2.0 / 3.0)
                        for k in self.sea_q}
        self.main_vol = 0.0
        if self.parted:
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

    def table(self, runoff_mm, lake_net_mm):
        """Each lake's balance (m3/s) with its water at each tabled level, for a year, beside the inflow the network
        gathers from the ground the design leaves dry: the rain less evaporation of the cells below the level, less
        the runoff those cells yield while dry (lake_rule.hollow_terms' `want`), with the runoff of all the cells the
        design's lake covers. At the design's level that is the lake's rain less evaporation; a lake that rises takes
        the runoff of the ground it covers, and one that falls gets the runoff of the bed it uncovers."""
        cells = self.cells
        area = self.net.area[cells]
        dry = np.maximum(np.asarray(runoff_mm).ravel()[cells], 0.0) * MM_DAY * area
        net = np.asarray(lake_net_mm).ravel()[cells] * MM_DAY * area
        cum = np.concatenate([[0.0], np.cumsum(net - dry)])
        back = np.bincount(self.filled_slot, weights=dry[self.filled], minlength=self.lake_nodes.size)
        return cum[self.pos] - cum[self.s0] + back[:, None]

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
        src = net.sources(runoff_mm, lake_net_mm)
        local, direct, lake = src['local'], src['direct'], src['lake']
        loss = water_losses(net, lake_net_mm)
        ntab = self.table(runoff_mm, lake_net_mm)
        slot = self.slot
        total = np.zeros(net.n_nodes)
        kinds_out = ('moscoviense', 'orientale', 'smythii_marginis') + (('neck',) if self.parted else ())
        totals = {k: 0.0 for k in kinds_out}
        supply_total = np.zeros(net.n_waters)
        for _ in range(steps):
            inflow = local.copy()
            sea_in = direct.copy()
            out = np.zeros(net.n_nodes)
            flows = {k: 0.0 for k in kinds_out}

            def put(kinds):
                q = op_flows(s, flows['moscoviense'], flows['orientale'], flows['smythii_marginis'], flows.get('neck', 0.0))
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
                    q = np.where(net.through[nodes], inflow[nodes], np.maximum(inflow[nodes] + lake[nodes], 0.0))
                    has = i >= 0
                    if has.any():
                        j = i[has]
                        v = self.vol[j]
                        bal = ntab[j, 0].copy()                  # at its pass or over it: the whole water's balance
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
                elif st_ == 2 and self.parted:
                    iN = self.iw['near_side_main']
                    self.main_vol, flows['neck'] = self.main_step(self.main_vol, sea_in[iN] - self.main_offset,
                                                                  loss[iN], dt)
            total += out
            for k in totals:
                totals[k] += flows[k]
            supply_total += sea_in
        return dict(out=total / steps, outflow={k: v / steps for k, v in totals.items()},
                    supply=supply_total / steps, loss=loss, sources=src)

    def field(self, runoff_mm, lake_net_mm, result):
        """A year's discharge in each cell from its mean outflows (a replay with each lake passing on its mean)."""
        q, _ = self.net.replay(result['sources'], ops_of(self.s, result['outflow']), passed_on=result['out'])
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
                                outflow.get('neck', 0.0)))


def year_summary(supply, loss, outflow, q, into):
    """A year's waters: what each receives and loses (m3/s), the upper seas' and the neck's outflows and the rivers'
    flow into the seas."""
    return dict(supply_m3s={k: round(float(v), 1) for k, v in zip(WATERS, supply)},
                loss_m3s={k: round(float(v), 1) for k, v in zip(WATERS, loss)},
                outflow_m3s={k: round(float(v), 1) for k, v in outflow.items()},
                rivers_into_the_seas_m3s=round(float(q[into].sum()), 1))


def checkpoint(path, done, flows, summaries, stores, key):
    """The years finished so far, with the storage's state, kept on the research drive under the design routing
    record's hash."""
    WORK.mkdir(parents=True, exist_ok=True)
    state = {f'vol_{k}': st.vol for k, st in stores.items()}
    meta = dict(done=done, summaries=summaries, design_routing=key,
                seas={k: dict(sea_vol=st.sea_vol, main_vol=st.main_vol) for k, st in stores.items()})
    np.savez(path, metadata=json.dumps(meta), **{f'flow_{k}': v for k, v in flows.items()}, **state)


def run_years(net, s, reference, design_in, track, into, land, log=print, resume=False):
    """The settled years in order: each routed without storage ('none', the lakes the design leaves below their passes
    passing nothing) and with the storage at both outlet widths ('narrow', 'wide'), the storage starting from the
    design state. Returns each scenario's flows at the tracked cells (years by cells), each year's summaries, the
    storage's check on the design climate, the years' mean climate on the climate model's grid, the estimate's steps
    and each year's runoff from the land (m3/s, the design's dry ground)."""
    m = s['meta']
    rotation = np.asarray(m['rotation'], float)
    rows, cols = net.rows, net.cols
    years = list(range(YEARS[0], YEARS[1] + 1))
    sums, applied, delta, ratios = year_sums()
    runoff_d, lake_d = design_in
    stores, steady = {}, {}
    for label, width in zip(('narrow', 'wide'), RL.WIDTH):
        st = Storage(net, s, reference, width, log)
        st.initial(reference, runoff_d, lake_d, reference['inflow'])
        t0 = time.time()
        check = st.year(runoff_d, lake_d)
        lakes = st.lake_nodes
        moved = np.abs(check['out'][lakes] - reference['out'][lakes])
        steady[label] = dict(
            buffering=buffering(st, reference), seconds_a_year=round(time.time() - t0, 1),
            largest_change_of_a_lake_outflow_m3s=round(float(moved.max()), 3),
            outflows_m3s={k: [round(reference['outflow'][k], 1), round(v, 1)] for k, v in check['outflow'].items()})
        st.initial(reference, runoff_d, lake_d, reference['inflow'])
        stores[label] = st
        log(f'storage at width {width:g} on the design climate: {steady[label]}')
    path = WORK / 'years.npz'
    key = sha256(DRt.OUT)
    flows = {k: np.zeros((len(years), track.size), np.float32) for k in ('none', 'narrow', 'wide')}
    summaries = {k: [] for k in flows}
    done = 0
    if resume and path.exists():
        with np.load(path, allow_pickle=False) as z:
            meta = json.loads(str(z['metadata']))
            if meta.get('design_routing') != key or z['flow_none'].shape != flows['none'].shape:
                raise SystemExit(f'{path} holds the years of another design routing; run without --resume')
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
        clim = year_climate(sums[y], applied, delta, rotation, ratios)
        fields = ('pr_mm_day', 'evap_mm_day', 'lsm')
        clim_sum = ({k: clim[k].astype(np.float64) for k in fields} if clim_sum is None
                    else {k: clim_sum[k] + clim[k] for k in fields})
        runoff_y, lake_y, _, _ = loaded(clim, m['evaporation'], rows, cols)
        land_runoff.append(float((np.maximum(runoff_y.ravel().astype(np.float64), 0.0) * MM_DAY * net.area)[land].sum()))
        if i < done:
            continue
        r = without_storage(net, s, runoff_y, lake_y, reference)
        flows['none'][i] = r['q'][track]
        summaries['none'].append(dict(year=y, **year_summary(r['supply'], r['loss'], r['outflow'], r['q'], into)))
        del r
        for k, st in stores.items():
            res = st.year(runoff_y, lake_y)
            q = st.field(runoff_y, lake_y, res)
            flows[k][i] = q[track]
            summaries[k].append(dict(year=y, **year_summary(res['supply'], res['loss'], res['outflow'], q, into)))
            del q, res
        checkpoint(path, i + 1, flows, summaries, stores, key)
        log(f'year {y}: {time.time() - t0:.0f} s; without storage {summaries["none"][-1]["outflow_m3s"]}')
    mean = {k: (v / len(years)).astype(np.float32) for k, v in clim_sum.items()}
    mean.update(lat=clim['lat'], lon=clim['lon'])
    return flows, summaries, steady, mean, (applied, delta), land_runoff


def near_side_balance(s, reference):
    """The near-side sea's parts on the design climate as the years' procedure routes it: what the rivers bring each
    part and what it loses at its level, against the flow the design sends from the main part to the southern one
    (through the sea's water where it stands as one sea across the neck, over the crest onto the ground past it
    where it stands parted, the spill then counted among the southern part's rivers)."""
    m = s['meta']
    iN, iS = WATERS.index('near_side_main'), WATERS.index('near_side_south')
    supply, loss, neck = reference['supply'], reference['loss'], m['neck_flow']
    through_water = neck if m['joined'] else 0.0
    return dict(joined_across_the_neck=m['joined'], neck_flow_m3s=round(neck, 1),
                main=dict(rivers_m3s=round(float(supply[iN]), 1), loses_m3s=round(float(loss[iN]), 1),
                          left_m3s=round(float(supply[iN] - loss[iN] - neck), 1)),
                south=dict(rivers_m3s=round(float(supply[iS]), 1), loses_m3s=round(float(loss[iS]), 1),
                           left_m3s=round(float(supply[iS] + through_water - loss[iS]), 1)),
                note=('left: what a part receives, with the neck\'s flow taken from the main part and given to the '
                      'southern one, beyond what it loses; the design routing settles each part within 0.2% of what '
                      'it receives'))


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
    for the upper seas, the main part where the near-side sea stands parted, and the lakes passing on THRESHOLD_M3S
    or more."""
    out = {}
    for k, q in st.sea_q.items():
        head = (q / st.sea_k[k]) ** (2.0 / 3.0)
        out[k] = dict(head_m=round(head, 2), years=round(st.sea_area[st.iw[k]] * head / (1.5 * q) / YEAR_S, 2))
    if st.parted:
        e = float(np.interp(reference['outflow']['neck'], st.neck_q, st.neck_head))
        slope = float(np.interp(e + 0.05, st.neck_head, st.neck_q) - np.interp(e - 0.05, st.neck_head, st.neck_q)) / 0.1
        out['near_side_main'] = dict(head_m=round(e, 2),
                                     years=round(st.sea_area[st.iw['near_side_main']] / slope / YEAR_S, 2))
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
    ds = design_state(log)
    s, routing = ds['s'], ds['routing']
    m = s['meta']
    rotation = np.asarray(m['rotation'], float)
    net = Network(ds['hollows'], ds['state'], ds['r'], ds['courses'], s['water'], len(WATERS), s['h'].shape, log)
    runoff_d, lake_d = ds['runoff'], ds['lake_net']
    del ds
    rows, cols = net.rows, net.cols
    q_design = np.asarray(routing['discharge_m3s'], float)
    wet = routing['routing_lake'] & ~net.outlet
    into = np.flatnonzero((net.water[net.receiver] >= 0) & (net.water < 0) & ~wet)
    # the network with the design's own lake outflows must give the design routing's flows back, and the years'
    # procedure on the design climate the same within the rounding of each lake's balance
    q_exact, _ = net.replay(net.sources(runoff_d, lake_d), ops_list(s, s['op_design']), passed_on=net.design_out)
    reference = without_storage(net, s, runoff_d, lake_d)
    reference['inflow'] = net.cascade(reference['sources'], ops_of(s, reference['outflow']))[0]
    land = ~wet & ~net.outlet
    major_d = land & (q_design >= THRESHOLD_M3S)
    designed = dict(moscoviense=m['surplus'], **m['over'], neck=m['neck_flow'])
    design_check = dict(
        routing_made_again_largest_difference_m3s=m['made_again_largest_difference_m3s'],
        replay_largest_difference_m3s=float(np.abs(q_exact - q_design).max()),
        years_procedure_largest_difference_m3s=round(float(np.abs(reference['q'] - q_design).max()), 3),
        years_procedure_outflows_m3s={k: [round(designed[k], 1), round(v, 1)] for k, v in reference['outflow'].items()},
        rivers_into_the_seas_m3s=round(float(q_design[into].sum())),
        lakes=dict(count=int((~net.through[net.n_out:]).sum()), below_their_passes=int(net.closed.sum()),
                   lowest_points_inside_larger_lakes=int(net.through.sum())),
        near_side=near_side_balance(s, reference))
    del q_exact
    log(f'design check: {design_check}, {time.time() - started:.0f} s')
    track = np.flatnonzero(land & (q_design >= TRACK_M3S))
    flows, summaries, steady, clim_mean, (applied, delta), land_runoff = run_years(
        net, s, reference, (runoff_d, lake_d), track, into, land, log, resume)
    log(f'the years routed, {time.time() - started:.0f} s')

    # the years' mean against the design: the climate's part (the years' mean climate routed) and the routing's
    runoff_m, lake_m, _, _ = loaded(clim_mean, m['evaporation'], rows, cols)
    mean_routed = without_storage(net, s, runoff_m, lake_m, reference)
    del lake_m
    at = np.flatnonzero(major_d[track])
    mean_of = lambda k, f: float(np.mean([f(y) for y in summaries[k]]))                     # noqa: E731
    q_at = q_design[track][at]
    against_design = dict(
        major_cells=int(at.size),
        years_mean_over_design=ratio_spread(flows['none'].mean(axis=0)[at], q_at),
        mean_climate_routed_over_design=ratio_spread(mean_routed['q'][track][at], q_at),
        years_mean_over_mean_climate_routed=ratio_spread(flows['none'].mean(axis=0)[at], mean_routed['q'][track][at]),
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
        r = without_storage(net, s, runoff_d, lake_d, reference, river_loss=field)
        after[end] = dict(q=r['q'], supply=r['supply'], outflow=r['outflow'],
                          taken_m3s=round(float(field.sum()), 1), into=round(float(r['q'][into].sum()), 1))
        del r
    loss_record['routed'] = dict(
        rule=('each loss taken from the river cell it falls on (a floodplain\'s from the nearest river cell), never more '
              'than the cell carries, and the year routed again: a lake or sea below receives less and Moscoviense '
              'passes on less by its fixed shares'),
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
            for k in ('none', 'narrow', 'wide'):
                row[k] = dict(mean_m3s=round(float(stats[k]['mean'][je]), 1), lowest_m3s=round(float(stats[k]['lowest'][je]), 1),
                              second_lowest_m3s=round(float(stats[k]['second'][je]), 1),
                              years_at_the_threshold=int(stats[k]['years'][je]), sustained=bool(stats[k]['sustained'][je]))
            row['after_losses_m3s'] = dict(low=round(float(q_after['low'][je]), 1), high=round(float(q_after['high'][je]), 1),
                                           lowest_year_high=round(float(max(stats['wide']['lowest'][je] - taken_high[je], 0.0)), 1))
            cell_reach[j[piece]] = len(reach_rows)
            reach_rows.append(row)
    log(f'{len(stem_list)} rivers in {len(reach_rows)} reaches, {time.time() - started:.0f} s')
    from research.studies.atlas_round import design_works as DW
    return dict(s=s, net=net, routing=routing, track=track, flows=flows, summaries=summaries, steady=steady,
                works=DW.OUT, works_product=DW.PRODUCT,
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
                'scaling to 300.31 K and then onto the design\'s own surface), each routed on the design routing\'s '
                'structure held fixed, the lakes as the corrected rule fills them, without storage and with the lakes '
                'and the upper seas carrying their volume from year to year over outlets sized by Earth\'s channel '
                'widths; the rivers\' losses in the drylands by river_losses.py\'s estimate. Fifteen years is a short '
                'record, and no figure is validated against an observation.')
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
    for k in ('none', 'narrow', 'wide'):
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
        inputs={rel(p): sha256(p) for p in (DRt.OUT, FD.OUT, N.OUT, b['works'], RL.OUT, ROOT / FD.ESTIMATE,
                                            ROOT / FD.SURFACE_ESTIMATE, b['frame_path'])},
        ignored_inputs={**{str(DRt.PRODUCT): sha256(DRt.PRODUCT), str(b['works_product']): sha256(b['works_product']),
                           **{f: sha256(ROOT / f) for f in (FD.CLIMATOLOGY, FD.BASE_TURNED_CLIMATOLOGY,
                                                            FD.TURNED_CLIMATOLOGY)},
                           rel(b['t21']): sha256(b['t21'])},
                        **{str(f): sha256(f) for f in year_files()}},
        evidence=evidence,
        reading_rule=('Flows in m3/s, lengths in km, areas in km2; a year is a settled year of the design climate (model '
                      'years 45-59). Each year is routed three ways: none, with no water carried between years and the '
                      'lakes the design leaves below their passes passing nothing, the lower bound on a dry year\'s flow '
                      'for the waters the design holds full; narrow and wide, with the lakes\' and the upper seas\' '
                      'storage carried from year to year over outlets of width 3 and 15 m per (m3/s)^0.5 widened for '
                      'lunar gravity, wide holding the least, a lake below its pass rising in wet years and passing '
                      'water on once it reaches it. A reach is part of a '
                      'major river (threshold_m3s or more on the design routing) between confluences, lakes of '
                      'lake_break_km2 or more and changes of class, its flows read at its downstream end; it is '
                      'sustained when it carries the threshold in every year with the wide outlets\' storage and '
                      'intermittent otherwise. A river is a stem from its head to where it ends or joins a larger river, '
                      'classed by its end. after_losses gives the design climate\'s flow after the drylands\' losses at '
                      'the low and high ends, and lowest_year_high the wide storage\'s lowest year less the high end\'s '
                      'loss of the mean year. product holds the cells, reaches the full table.'),
        limits=[
            'The climate is an estimate in two steps: the run at sunlight 0.936 carried by pattern scaling from 301.00 K '
            'to 300.31 K, each year keeping its own departures from the run\'s mean, and then each cell\'s rain and '
            'evaporation moved by the ratio the step onto the design\'s own surface gives the mean '
            '(climate.gcm.surface_change, from existing runs); a run on that surface at the design temperature would '
            'differ by more than the noise of fifteen years\' means where the lakes grew.',
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
            'The structure is held fixed: the lakes keep the design\'s shores in the years without storage, a lake that '
            'reaches its pass with storage sends its water on by the course the design gives it and not into the larger '
            'hollow it would fill with its neighbour, and each cut carries the design\'s part of its branch\'s share '
            'whether or not its lake can spare it that year. Where the near-side sea is one sea across its neck, both '
            'parts are held at their design levels.',
            'The routing is the design routing\'s on 1.9 km cells, each hollow filled as one water (lake_rule.py). '
            'Nothing is validated against an observation.'],
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
