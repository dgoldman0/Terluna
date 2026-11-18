"""A corrected lake rule tried on the design surface: every hollow filled as one water, and what that changes.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python3 -m research.studies.atlas_round.lake_rule   # about 4 minutes

The routing the design surface is drawn from (research/studies/drainage_network/routing.py's route_levels, with
geography/drainage.py's lake rule) balances each catchment, the cells draining to one lowest point, by itself, and fills
a catchment that has water to spare to the lowest pass on its way to a sea. A hollow with several lowest points is then
filled piece by piece: one piece stands full to the hollow's outer pass while the piece beside it keeps a small lake or
none, with nothing between the two but the hollow's floor (lake_support.py measures it: dry ground below an adjoining
lake over 9% of the Moon). Water does not stand so. The full piece pours over the sill between the two, and either the
second piece takes all of it and the river ends there, or both fill to the sill, join, and rise together as one lake
until its loss to the air meets what comes in or it reaches the next pass.

route_hollows is the routing with that rule. The catchments are joined into the tree of the hollows they make: taking
the passes between catchments from the lowest up, two hollows that a pass joins make a larger hollow, whose water
covers both once each has filled to the pass, and a hollow that a pass joins to ground already draining to a sea
overflows there. Each hollow needs a supply to stand full: its open water's loss beyond the rain on it, and the runoff
its drowned ground no longer yields. A hollow whose supply is less holds the lake that its supply keeps, lowest ground
first, and passes nothing on; one whose supply is more stands full and sends the rest over its pass, into the hollow
beside it while that one is not full, and into the larger hollow the two make once it is. Flows put in at cells (the
overflows of the full seas, the works' branches, the spill over the near-side sea's neck) join the supply of the
hollow they land in, and flows taken out leave it. Ground that stays dry yields its runoff whether or not it lies
below a pass (route_levels yields none below a hollow's outer pass); floor_runoff=False keeps route_levels' way, which
parts the rule's effect from this one's.

The module routes the design climate once on the design routing's own floor, seas and works (the structure
sustained_rivers.py keeps), with the seas held at their design levels and the works' flows at their design values,
and compares the result with the design routing: the lakes' area, what each sea receives against what it loses, and
the flow of each named waterway where it ends and at its head. It first routes with route_levels on the same inputs,
which must give the design routing's discharge back. Nothing in the design is changed: the seas' levels, the weirs'
shares and the factor on the open water's evaporation that closes the design's budget were all set on the present
rule, and setting them again on a corrected rule is the step this measure is for.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from climate.gcm.turned import to_turned
from geography import drainage as dr
from research.studies.atlas_round import design_routing as DRt
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import lake_support as LS
from research.studies.atlas_round import sustained_rivers as SR
from research.studies.drainage_network.routing import route_levels

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'lake_rule.json'
WATERWAYS = HERE / 'results' / 'waterways.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/lake_rule')
PRODUCT = BULK / 'lake_rule_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-lake-rule/1'
CODE = ('research/studies/atlas_round/lake_rule.py', 'research/studies/drainage_network/routing.py',
        'geography/drainage.py')
MM_DAY = dr.MM_DAY
CLOSED_LISTED = 12                   # the largest lakes that pass nothing on, listed
ON_A_COURSE_KM2 = 100.0              # a lake that passes nothing on is listed on a waterway's course from this size
LOST_SHARE = 0.05                    # a waterway's course is followed where its mouth loses more than this of its flow
TRACE_LISTED = 8                     # the first lakes listed where a lost waterway's water is followed


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------- the hollows' tree

def hollow_tree(n_nodes, n_out, edges, out_level):
    """The tree of the hollows the catchments make. Nodes below n_out are the seas' outlets, the rest catchments; edges
    are the lowest passes between touching nodes (geography.drainage.passes). Taking the passes from the lowest up,
    a pass between two hollows makes a larger hollow of them, and a pass onto ground that already drains to a sea
    is where a hollow overflows for good. Returns arrays over the hollows (the catchments keep their node numbers and
    the larger hollows follow): parent (-1 for a hollow that overflows towards a sea), spill (the level it overflows
    at), land (the cell its overflow lands on, across the pass), sibling (-1 without one) and the two hollows each
    larger one is made of; and the hollows that overflow towards a sea, in the order they were reached."""
    u, v, w, cu, cv = edges
    lev_u = np.where(u < n_out, out_level[np.minimum(u, n_out - 1)], -np.inf)
    lev_v = np.where(v < n_out, out_level[np.minimum(v, n_out - 1)], -np.inf)
    cross = np.maximum(w.astype(np.float64), np.maximum(lev_u, lev_v))      # a pass onto a sea is no lower than the sea
    order = np.argsort(cross, kind='stable')
    order = order[~((u < n_out) & (v < n_out))[order]]
    U, V, CU, CV, X = (a[order].tolist() for a in (u, v, cu, cv, cross))
    group = list(range(n_nodes))
    top = [-1] * n_out + list(range(n_out, n_nodes))    # the hollow a group has grown into; -1 once it drains to a sea
    parent = [-1] * n_nodes
    spill = [np.nan] * n_nodes
    land = [-1] * n_nodes
    first, second, reached = [], [], []
    for k in range(len(U)):
        a, b = U[k], V[k]
        while group[a] != a:
            group[a] = group[group[a]]
            a = group[a]
        while group[b] != b:
            group[b] = group[group[b]]
            b = group[b]
        if a == b:
            continue
        ta, tb = top[a], top[b]
        group[a] = b
        if ta < 0 and tb < 0:
            continue
        if ta < 0 or tb < 0:                             # one side drains to a sea: the other overflows onto it
            x, cell = (tb, CU[k]) if ta < 0 else (ta, CV[k])
            spill[x], land[x] = X[k], cell
            reached.append(x)
            top[b] = -1
            continue
        p = len(parent)                                  # two hollows joined by this pass make a larger one
        parent[ta] = parent[tb] = p
        spill[ta] = spill[tb] = X[k]
        land[ta], land[tb] = CV[k], CU[k]                # each one's overflow lands on the other's side
        parent.append(-1)
        spill.append(np.nan)
        land.append(-1)
        first.append(ta)
        second.append(tb)
        top[b] = p
    n_all = len(parent)
    parent, spill, land = np.asarray(parent, np.int64), np.asarray(spill, float), np.asarray(land, np.int64)
    first, second = np.asarray(first, np.int64), np.asarray(second, np.int64)
    sibling = np.full(n_all, -1, np.int64)
    sibling[first], sibling[second] = second, first
    kids = np.full((n_all, 2), -1, np.int64)
    kids[n_nodes:, 0], kids[n_nodes:, 1] = first, second
    if np.isnan(spill[n_out:]).any():
        raise RuntimeError(f'{int(np.isnan(spill[n_out:]).sum())} hollows have no way to a sea')
    return dict(parent=parent, spill=spill, land=land, sibling=sibling, kids=kids,
                reached=np.asarray(reached, np.int64), n_nodes=n_nodes, n_out=n_out)


def first_covered_by(tree, leaf, height):
    """For cells of catchments `leaf` at `height`: the smallest hollow whose water, standing full, covers each cell (the
    catchment itself where the cell lies below its pass), or -1 where no hollow's water reaches it."""
    parent, spill = tree['parent'], tree['spill']
    n_all = parent.size
    up = [np.where(parent >= 0, parent, np.arange(n_all))]
    while True:
        nxt = up[-1][up[-1]]
        if np.array_equal(nxt, up[-1]):
            break
        up.append(nxt)
    root = up[-1]
    node = leaf.astype(np.int64).copy()
    node[height >= spill[root[leaf]]] = -1               # above the outermost pass: never under water
    climb = np.flatnonzero((node >= 0) & (height >= spill[leaf]))
    y, hh = node[climb], height[climb]
    for table in reversed(up):                           # the highest hollow round the cell that it stands clear of
        c = table[y]
        y = np.where(spill[c] <= hh, c, y)
    node[climb] = parent[y]
    return node


# ---------------------------------------------------------------- the hollows' balance

def settle(tree, supply, need):
    """Each hollow's state on its supply. supply: the water (m3/s) each catchment gathers from its own ground and from
    what is put in; need: what each hollow needs beyond the hollows it is made of to stand full to its pass. A hollow
    is full when it is a catchment, or made of two full hollows, and has been given its need; what it is given beyond
    that goes over its pass: into the hollow beside it while that one is not full (to the catchment its overflow lands
    in), into the larger hollow the two make once it is, and from a hollow that overflows towards a sea to the
    catchment or sea its overflow lands on. Returns full (bool), held (the water each hollow keeps: its need when full,
    less when not), out (what leaves a hollow that overflows towards a sea) and over (what a full hollow has sent over
    its pass into the hollow beside it)."""
    parent, sibling, kids = tree['parent'].tolist(), tree['sibling'].tolist(), tree['kids'].tolist()
    n_nodes, n_out = tree['n_nodes'], tree['n_out']
    n_all = len(parent)
    land_leaf = tree['land_leaf'].tolist()
    held = np.zeros(n_all)
    held[:supply.size] = supply
    held, need = held.tolist(), need.tolist()
    full = bytearray(n_all)
    out, over = [0.0] * n_all, [0.0] * n_all
    # hollows that overflow towards a sea, the last reached first: each lands on ground reached before it, so every
    # hollow has all that comes from upstream before it is settled; within one, the smaller hollows first
    root = np.where(tree['parent'] >= 0, tree['parent'], np.arange(n_all))
    while True:
        nxt = root[root]
        if np.array_equal(nxt, root):
            break
        root = nxt
    rank = np.full(n_all, -1, np.int64)
    rank[tree['reached'][::-1]] = np.arange(tree['reached'].size)
    nodes = np.arange(n_out, n_all)
    nodes = nodes[np.lexsort((nodes, rank[root[nodes]]))]
    turn = rank[root].tolist()
    for x in nodes.tolist():
        if full[x]:
            continue
        if x >= n_nodes and not (full[kids[x][0]] and full[kids[x][1]]):
            continue
        stack = [(x, 0.0)]
        while stack:
            y, q = stack.pop()
            if not full[y]:
                held[y] += q
                if held[y] < need[y]:
                    continue
                q = held[y] - need[y]
                held[y] = need[y]
                full[y] = 1
            p = parent[y]
            if p < 0:                                    # over the pass towards a sea
                out[y] += q
                t = land_leaf[y]
                if t >= n_out:
                    if turn[t] <= turn[y]:
                        raise RuntimeError('an overflow lands in a hollow settled before it')
                    held[t] += q
            elif full[sibling[y]]:                       # both full: the water is the larger hollow's
                stack.append((p, q))
            elif q > 0.0:                                # over the sill into the hollow beside it
                over[y] += q
                stack.append((land_leaf[y], q))
    return np.frombuffer(full, np.uint8).astype(bool), np.asarray(held), np.asarray(out), np.asarray(over)


# ---------------------------------------------------------------- the routing

def accumulate(receiver, donor, start):
    """Flow gathered down the receivers: each donor cell passes all it holds to its receiver. Returns each cell's
    flow, its own included."""
    acc = start.astype(np.float64).copy()
    indeg = np.bincount(receiver[donor], minlength=receiver.size).astype(np.int32)
    front = np.flatnonzero((indeg == 0) & donor)
    while front.size:
        r = receiver[front]
        np.add.at(acc, r, acc[front])
        np.subtract.at(indeg, r, 1)
        cand = np.unique(r)
        front = cand[(indeg[cand] == 0) & donor[cand]]
    return acc


def route_hollows(h, outlet, lat_deg, cell_km, runoff_mm_day, lake_net_mm_day, outlet_level, sources=(), floor_runoff=True):
    """route_levels with every hollow filled as one water (see the module's text). sources: (cell, m3/s) pairs put in
    (or, negative, taken out) at cells. Returns a dict of flat arrays: receiver (each cell's steepest-descent
    neighbour, itself at a lowest point and on a sea), discharge (m3/s through each cell: dry ground's runoff and the
    overflows, gathered downhill; on a sea's cells, what arrives there), own_flow (what starts at each cell), lake,
    lake_level (NaN off the lakes), lake_spills (whether the cell's lake overflows), body (the hollow whose water a lake cell is part of), label, area, net (the open
    water's balance, m3/s), outlet; the hollows' tables under 'tree', 'full', 'held', 'need', 'out', 'over', 'left'
    (what a lake not full keeps beyond its last whole cell) and 'cover' (the hollow whose water first covers each cell); and
    the budget's pieces."""
    rows, cols = h.shape
    n = rows * cols
    idx2 = np.arange(n, dtype=np.int32).reshape(rows, cols)
    receiver = dr.d8_receivers(h.astype(np.float32), outlet, lat_deg, cell_km)
    flat_out = outlet.ravel()
    sink = (receiver == idx2.ravel()) & ~flat_out
    sinks = np.flatnonzero(sink).astype(np.int32)
    node_of = np.full(n, -1, np.int32)
    levels, groups = np.unique(np.asarray(outlet_level, dtype=float).ravel()[flat_out], return_inverse=True)
    n_out = max(1, levels.size)
    node_of[flat_out] = groups
    node_of[sinks] = np.arange(n_out, sinks.size + n_out, dtype=np.int32)
    label = node_of[dr.terminal(receiver)]
    n_nodes = sinks.size + n_out
    h_pass = np.where(outlet, -np.inf, h).astype(np.float32)
    edges = dr.passes(label.reshape(rows, cols), h_pass, idx2)
    out_level = levels if levels.size else np.array([float(h.min())])
    tree = hollow_tree(n_nodes, n_out, edges, out_level)
    del edges, h_pass
    land = tree['land']
    tree['land_leaf'] = np.where(land >= 0, label[np.maximum(land, 0)], -1)

    hf = h.ravel().astype(np.float64)
    area = (((cell_km * 1000.0) ** 2 * np.maximum(np.cos(np.radians(lat_deg)), 1e-6))[:, None] * np.ones((1, cols))).ravel()
    ground = np.flatnonzero(~flat_out)
    hollow = ground[label[ground] >= n_out]              # ground that drains straight to a sea lies in no hollow
    cover = np.full(n, -1, np.int64)
    cover[hollow] = first_covered_by(tree, label[hollow], hf[hollow])
    run = np.where(flat_out, 0.0, np.maximum(runoff_mm_day.ravel(), 0.0) * MM_DAY * area)
    if not floor_runoff:
        run[cover >= 0] = 0.0
    net = lake_net_mm_day.ravel() * MM_DAY * area
    n_all = tree['parent'].size
    under = np.flatnonzero(cover >= 0)
    order = np.lexsort((hf[under], cover[under]))
    under = under[order]                                 # the cells a hollow's water can cover, by hollow and height
    owner = cover[under]
    want = run[under] - net[under]                       # what keeps a cell under water: its loss, and its lost runoff
    need = np.bincount(owner, weights=want, minlength=n_all)
    supply = np.bincount(label[ground], weights=run[ground], minlength=n_nodes)
    put = np.zeros(n)
    for cell, q in sources:
        put[int(cell)] += q
        if not flat_out[int(cell)]:
            supply[label[int(cell)]] += q
    supply[:n_out] = 0.0
    full, held, out, over = settle(tree, supply, need)

    # the lakes: all a full hollow's water can cover, and in a hollow not full the lowest cells its supply keeps
    kids = tree['kids']
    leaf = np.arange(n_all) < n_nodes
    made = leaf.copy()
    made[n_nodes:] = full[kids[n_nodes:, 0]] & full[kids[n_nodes:, 1]]
    made[:n_out] = False
    start = np.searchsorted(owner, np.arange(n_all + 1))
    wet_sorted = full[owner]
    level = np.full(n_all, np.nan)                       # the water's level in a hollow that is not full
    level[n_nodes:] = np.where(made[n_nodes:] & ~full[n_nodes:], tree['spill'][kids[n_nodes:, 0]], np.nan)
    left = np.zeros(n_all)                               # supply a lake not full keeps beyond its last whole cell
    for x in np.flatnonzero(made & ~full).tolist():
        a, b = start[x], start[x + 1]
        left[x] = held[x]
        if a == b:
            continue
        reach = np.maximum.accumulate(np.cumsum(want[a:b]))
        k = int(np.searchsorted(reach, held[x], side='right'))
        if k:
            wet_sorted[a:a + k] = True
            level[x] = hf[under[a + k - 1]]
            left[x] = held[x] - float(np.cumsum(want[a:a + k])[-1])
    lake = np.zeros(n, bool)
    lake[under] = wet_sorted
    # a full hollow's water stands at its pass, or with the hollow beside it at the level of the larger one they make
    top = np.where(full & (tree['parent'] >= 0) & full[np.maximum(tree['sibling'], 0)], tree['parent'], np.arange(n_all))
    while True:
        nxt = top[top]
        if np.array_equal(nxt, top):
            break
        top = nxt
    stands = np.where(full, np.where(full[top], tree['spill'][top], level[top]), level)
    lake_level = np.full(n, np.nan)
    lake_level[under] = np.where(wet_sorted, stands[owner], np.nan)
    spills = np.zeros(n, bool)
    spills[under] = wet_sorted & full[top[owner]]
    body = np.full(n, -1, np.int64)                      # the hollow whose water each lake cell is part of
    body[under] = np.where(wet_sorted, top[owner], -1)

    # the rivers: dry ground's runoff, what is put in, and each overflow where it lands, gathered downhill
    flow = np.where(lake | flat_out, 0.0, run) + put
    roots = np.flatnonzero(full & (tree['parent'] < 0) & (np.arange(n_all) >= n_out))
    np.add.at(flow, land[roots], out[roots])
    sills = np.flatnonzero(full & (tree['parent'] >= 0) & ~full[np.maximum(tree['sibling'], 0)])
    np.add.at(flow, land[sills], over[sills])
    discharge = accumulate(receiver, ~flat_out & ~sink, flow)
    budget = dict(
        runoff_of_dry_ground=float(run[~lake & ~flat_out].sum()), put_in=float(put.sum()),
        lakes_rain_less_evaporation=float(net[lake].sum()), reaches_the_seas=float(discharge[flat_out].sum()),
        kept_by_lakes_not_full=float(left.sum()))
    return dict(receiver=receiver, label=label, discharge=discharge, lake=lake, lake_level=lake_level,
                lake_spills=spills, area=area, outlet=flat_out, sinks=sinks, outlet_nodes=n_out, tree=tree, full=full,
                held=held, need=need, out=out, over=over, cover=cover, body=body, made=made, left=left, own_flow=flow, net=net,
                budget=budget)


# ---------------------------------------------------------------- the design surface

def design_inputs(log=print):
    """The design routing's last routing as design_routing.route makes it, from the structure sustained_rivers.py
    keeps: the floor, the waters and their levels, the design climate's runoff and open-water balance on the closed
    budget, and the put-ins and take-outs at their design flows."""
    key = dict(record=sha256(DRt.OUT), code={f: sha256(ROOT / f) for f in SR.STRUCTURE_CODE})
    if not SR.STRUCTURE.exists():
        raise SystemExit(f'{SR.STRUCTURE} is missing; run python -m research.studies.atlas_round.sustained_rivers')
    with np.load(SR.STRUCTURE, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        if meta['key'] != key:
            raise SystemExit(f'{SR.STRUCTURE} is not the present design routing\'s; run '
                             'python -m research.studies.atlas_round.sustained_rivers, which builds it again')
        s = {k: z[k] for k in z.files if k != 'metadata'}
    log(f'the design structure from {SR.STRUCTURE}')
    h = s['h'].astype(np.float32)
    rows, cols = h.shape
    water = s['water']
    levels = [float(v) for v in meta['levels']]
    runoff, lake_net, _, _ = FD.load(meta['evaporation'])
    ops = [(int(c), float(q), bool(lk)) for c, q, lk in zip(s['op_cell'], s['op_design'], s['op_lake'])]
    return dict(h=h, lat=s['lat'], cell_km=float(meta['cell_km']), water=water, levels=levels, runoff=runoff,
                lake_net=lake_net, ops=ops, op_kind=[SR.OPS[k] for k in s['op_kind']], level0=float(meta['level0']),
                loss=meta['loss'], rotation=np.asarray(meta['rotation'], float))


def present_routing(d):
    """route_levels on the design's inputs as governed_drainage.route_with gives them to it: a flow put in on the
    ground as runoff, one put into or taken out of a full hollow in that catchment's balance."""
    h = d['h']
    rows, cols = h.shape
    area = (((d['cell_km'] * 1000.0) ** 2 * np.maximum(np.cos(np.radians(d['lat'])), 1e-6))[:, None] * np.ones((1, cols))).ravel()
    add_runoff = d['runoff'].astype(np.float64).copy()
    add_lake = d['lake_net'].astype(np.float64).copy()
    for cell, q, in_lake in d['ops']:
        i, j = divmod(cell, cols)
        (add_lake if in_lake else add_runoff)[i, j] += q / (MM_DAY * area[cell])
    outlet = (d['water'] >= 0).reshape(rows, cols)
    level = np.where(d['water'] >= 0, np.asarray(d['levels'])[np.maximum(d['water'], 0)], 0.0).reshape(rows, cols)
    return route_levels(h, outlet, d['lat'], d['cell_km'], add_runoff, add_lake, level)


def corrected_routing(d, floor_runoff=True):
    h = d['h']
    rows, cols = h.shape
    outlet = (d['water'] >= 0).reshape(rows, cols)
    level = np.where(d['water'] >= 0, np.asarray(d['levels'])[np.maximum(d['water'], 0)], 0.0).reshape(rows, cols)
    return route_hollows(h, outlet, d['lat'], d['cell_km'], d['runoff'], d['lake_net'], level,
                         sources=[(c, q) for c, q, _ in d['ops']], floor_runoff=floor_runoff)


def lakes_of(r, water, area_km2, total):
    """A routing's lakes off the seas: cells, area, share of the Moon, and the part in lakes that overflow."""
    lake = r['lake'] & (water < 0)
    spills = r['lake_spills'] & lake if 'lake_spills' in r else lake & ~r['closed'][r['label']]
    a = float(area_km2[lake].sum())
    return dict(cells=int(lake.sum()), area_km2=round(a), share_of_moon=round(a / total, 5),
                area_overflowing_km2=round(float(area_km2[spills].sum())),
                area_not_overflowing_km2=round(float(area_km2[lake & ~spills].sum())))


def unheld_shores(r, h, water, sea_levels, area_km2):
    """lake_support's measure on a routing: the lakes' shore cells with dry ground beside them more than
    lake_support.DROP_M below the lake's surface, and the dry ground below an adjoining lake and joined to it."""
    rows, cols = h.shape
    lake = (r['lake'] & (water < 0)).reshape(rows, cols)
    level = np.where(lake, r['lake_level'].reshape(rows, cols), np.nan).astype(np.float32)
    sea = (water >= 0).reshape(rows, cols)
    level[sea] = np.asarray(sea_levels, np.float32)[water[water >= 0]]     # a sea's cells are water, not dry ground
    loose = LS.unheld(h, level, lake)
    below = np.isfinite(LS.ground_below(h, level, lake))
    return dict(shore_cells_not_held=int(loose.sum()),
                ground_below_and_joined_km2=round(float(area_km2.reshape(rows, cols)[below].sum())))


def place(cell, cols, rotation):
    """A cell's centre on today's Moon and in the settled frame (the globe page's), as [latitude, longitude]."""
    lat, lon = 90.0 - (cell // cols + 0.5) / dr.PPD, (cell % cols + 0.5) / dr.PPD
    la, lo = to_turned(np.array([lat]), np.array([lon]), rotation)
    return dict(today=[round(float(lat), 2), round(float(lon), 2)],
                settled=[round(float(la[0]), 2), round(float(((lo[0] + 180.0) % 360.0) - 180.0), 2)])


def lake_bodies(r, h_flat, water, area_km2, was_lake, level0, cols, rotation):
    """The corrected routing's lakes off the seas as bodies (the hollow whose water each is): a table by body of the
    area, whether it overflows and its level, each cell's row in it, and the largest that pass nothing on."""
    lake = r['lake'] & (water < 0)
    cells = np.flatnonzero(lake)
    uniq, inv = np.unique(r['body'][cells], return_inverse=True)
    area = np.bincount(inv, weights=area_km2[cells])
    spills = np.zeros(uniq.size, bool)
    spills[inv] = r['lake_spills'][cells]
    level = np.zeros(uniq.size)
    level[inv] = r['lake_level'][cells]
    row = np.full(h_flat.size, -1, np.int64)
    row[cells] = inv
    by_body = np.argsort(inv, kind='stable')
    start = np.searchsorted(inv[by_body], np.arange(uniq.size + 1))

    def describe(k):
        own = cells[by_body[start[k]:start[k + 1]]]
        deepest = int(own[np.argmin(h_flat[own])])
        return dict(area_km2=round(float(area[k])), newly_under_water_km2=round(float(area_km2[own][~was_lake[own]].sum())),
                    level_above_common_m=round(float(level[k] - level0)), depth_m=round(float(level[k] - h_flat[deepest])),
                    deepest_at=place(deepest, cols, rotation))

    closed = np.argsort(-np.where(spills, 0.0, area))[:CLOSED_LISTED]
    largest = [dict(passes_water_on=bool(spills[k]), **describe(int(k))) for k in np.argsort(-area)[:CLOSED_LISTED]]
    return dict(row=row, area=area, spills=spills, describe=describe, count=int(uniq.size), ids=uniq, largest=largest,
                passing_nothing_on=int((~spills).sum()), largest_passing_nothing_on=[describe(int(k)) for k in closed])


def trace(r, cell, limit=200_000):
    """The water at a cell followed on a route_hollows routing: down the ground to a lake, over each overflowing lake's
    pass to where its overflow lands, and on. Returns the lakes it enters in order, each as (the hollow whose water it
    is, what that lake passes on in m3/s), and how it ends: at a sea ('sea', its cell), in a lake that passes nothing
    on ('lake', the last listed), or at a lowest point that holds no lake ('dry', its cell)."""
    receiver, lake, spills, body, outlet = r['receiver'], r['lake'], r['lake_spills'], r['body'], r['outlet']
    tree, out, over = r['tree'], r['out'], r['over']
    c = int(cell)
    lakes = []
    for _ in range(limit):
        if outlet[c]:
            return lakes, ('sea', c)
        if lake[c]:
            top = int(body[c])
            if not spills[c]:
                lakes.append((top, 0.0))
                return lakes, ('lake', top)
            lakes.append((top, float(out[top] if tree['parent'][top] < 0 else over[top])))
            c = int(tree['land'][top])
            continue
        nxt = int(receiver[c])
        if nxt == c:
            return lakes, ('dry', c)
        c = nxt
    return lakes, ('unsettled', c)


def course(receiver, head, mouth, limit=50_000):
    """The cells from a waterway's head down the design routing to its mouth."""
    path = [int(head)]
    while path[-1] != mouth and len(path) < limit:
        nxt = int(receiver[path[-1]])
        if nxt == path[-1]:
            break
        path.append(nxt)
    return np.asarray(path)


def waterway_cells(q, cols, record):
    """Each named waterway's mouth and head cells on the design routing: the cell within a cell of the recorded place
    whose discharge is nearest the recorded flow."""
    out = []
    for x in record['waterways']:
        cells = {}
        for end, flow in (('mouth', x['flow_mouth_m3s']), ('head', x['flow_head_m3s'])):
            lat, lon = x[end]['today']
            i0, j0 = int((90.0 - lat) * dr.PPD), int((lon % 360.0) * dr.PPD)
            best = None
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    c = (i0 + di) * cols + (j0 + dj) % cols
                    miss = abs(float(q[c]) - flow)
                    if best is None or miss < best[0]:
                        best = (miss, c)
            cells[end] = best
        out.append(cells)
    return out


def build(log=print):
    started = time.time()
    d = design_inputs(log)
    h, water = d['h'], d['water']
    rows, cols = h.shape
    design = DRt.load(log)
    present = present_routing(d)
    present = {k: present[k] for k in ('discharge', 'lake', 'lake_level', 'closed', 'label', 'area')}
    differs = float(np.abs(present['discharge'] - design['discharge_m3s']).max())
    design_receiver = design['receiver']
    del design
    log(f'route_levels on the kept structure: largest difference from the design routing {differs:.3g} m3/s, '
        f'{time.time() - started:.0f} s')
    if differs > 1e-3:
        raise SystemExit('the kept structure does not give the design routing back')
    area_km2 = present['area'] / 1e6
    total = float(area_km2.sum())
    names = list(SR.WATERS)
    masks = [water == i for i in range(len(names))]
    record = json.loads(WATERWAYS.read_text())
    ends = waterway_cells(present['discharge'], cols, record)
    cases = {}
    kept = {}
    for key, floor_runoff in (('corrected', True), ('corrected_without_floor_runoff', False)):
        r = corrected_routing(d, floor_runoff)
        b = r['budget']
        gap = b['runoff_of_dry_ground'] + b['put_in'] + b['lakes_rain_less_evaporation'] - b['reaches_the_seas'] - b['kept_by_lakes_not_full']
        log(f'{key}: lakes {lakes_of(r, water, area_km2, total)["share_of_moon"]:.4f} of the Moon, '
            f'{b["reaches_the_seas"]:,.0f} m3/s reach the seas, budget closes within {gap:.3g} m3/s, {time.time() - started:.0f} s')
        if floor_runoff:                                 # the corrected routing looked at more closely
            hf = h.ravel()
            bodies = lake_bodies(r, hf, water, area_km2, present['lake'], d['level0'], cols, d['rotation'])
            tree = r['tree']
            n_all = tree['parent'].size
            root = np.where(tree['parent'] >= 0, tree['parent'], np.arange(n_all))
            while True:
                nxt = root[root]
                if np.array_equal(nxt, root):
                    break
                root = nxt
            short = np.bincount(root, weights=r['need'], minlength=n_all) - np.bincount(root, weights=r['held'], minlength=n_all)
            hollow_of = {}
            for (cell, q, _), kind in zip(d['ops'], d['op_kind']):
                if q < 0.0:
                    top = int(root[r['label'][cell]])
                    hollow_of[cell] = dict(
                        outer_pass_above_common_m=round(float(tree['spill'][top] - d['level0']), 1), stands_full=bool(r['full'][top]),
                        overflows_m3s=round(float(r['out'][top]), 1),
                        short_of_standing_full_m3s=(0.0 if r['full'][top] else round(float(short[top]), 1)))
            bands = {}
            lake_now = r['lake'] & (water < 0)
            lake_was = present['lake'] & (water < 0)
            for top_lat in range(90, -90, -30):
                part = slice(int((90 - top_lat) * rows / 180) * cols, int((120 - top_lat) * rows / 180) * cols)
                bands[f'{top_lat} to {top_lat - 30}'] = dict(present_km2=round(float(area_km2[part][lake_was[part]].sum())),
                                                           corrected_km2=round(float(area_km2[part][lake_now[part]].sum())))
            changed = dict(newly_lake_km2=round(float(area_km2[lake_now & ~lake_was].sum())),
                           no_longer_lake_km2=round(float(area_km2[lake_was & ~lake_now].sum())),
                           lakes_volume_km3=dict(
                               present=round(float(((present['lake_level'] - hf) * area_km2)[lake_was].sum() / 1000.0)),
                               corrected=round(float(((r['lake_level'] - hf) * area_km2)[lake_now].sum() / 1000.0))))
            on_course = []
            for x, e in zip(record['waterways'], ends):
                mouth, head = e['mouth'][1], e['head'][1]
                lost = 1.0 - float(r['discharge'][mouth]) / max(float(present['discharge'][mouth]), 1e-9)
                found, drop = [], None
                if lost > LOST_SHARE:
                    path = course(design_receiver, head, mouth)
                    rows_on = bodies['row'][path]
                    seen = []
                    for k in rows_on[rows_on >= 0].tolist():
                        if k not in seen:
                            seen.append(k)
                    found = [bodies['describe'](k) for k in seen
                             if not bodies['spills'][k] and bodies['area'][k] >= ON_A_COURSE_KM2]
                    # where the corrected routing loses most beyond what the design loses, between dry cells of
                    # the design course, and the water above that place followed lake by lake
                    dry = path[~lake_now[path] & (water[path] < 0)]
                    if dry.size > 1:
                        flow, was = r['discharge'][dry], present['discharge'][dry]
                        k = int(np.argmax((flow[:-1] - flow[1:]) - (was[:-1] - was[1:])))
                        lakes, (kind, where) = trace(r, dry[k])
                        rows_in = np.searchsorted(bodies['ids'], [t for t, _ in lakes])
                        drop = dict(above_m3s=round(float(flow[k])), below_m3s=round(float(flow[k + 1])),
                                    design_above_m3s=round(float(was[k])), design_below_m3s=round(float(was[k + 1])),
                                    above_at=place(int(dry[k]), cols, d['rotation']), water_above_ends=kind,
                                    lakes_it_enters=len(lakes),
                                    first_lakes=[dict(passes_on_m3s=round(q, 1), **bodies['describe'](int(row)))
                                                 for (_, q), row in list(zip(lakes, rows_in))[:TRACE_LISTED]])
                        if kind in ('dry', 'sea'):
                            drop['ends_at'] = place(int(where), cols, d['rotation'])
                        elif kind == 'lake':
                            drop['ends_in'] = bodies['describe'](int(rows_in[-1]))
                on_course.append(dict(lost_share=round(lost, 4), lakes_passing_nothing_on=found, largest_fall=drop))
            closer = dict(lakes_as_bodies=bodies['count'], bodies_passing_nothing_on=bodies['passing_nothing_on'],
                          largest=bodies['largest'],
                          largest_passing_nothing_on=bodies['largest_passing_nothing_on'], lakes_by_latitude=bands, **changed)
            del bodies
        cases[key] = {k: r[k] for k in ('discharge', 'lake', 'lake_level', 'lake_spills')}
        kept[key] = dict(budget_m3s={k: round(v, 1) for k, v in b.items()}, budget_gap_m3s=round(gap, 3),
                         hollows=dict(catchments=int(r['tree']['n_nodes'] - r['tree']['n_out']),
                                      made_of_two=int(r['tree']['parent'].size - r['tree']['n_nodes']),
                                      full=int(r['full'].sum()), holding_a_lake_not_full=int((r['made'] & ~r['full'] & (r['held'] > 0)).sum())))
        del r
    out_cases = {}
    for key, r in (('present', present), *cases.items()):
        q = r['discharge']
        supplies = {k: float(q[m].sum()) for k, m in zip(names, masks)}
        lakes = lakes_of(r, water, area_km2, total)
        lakes.update(unheld_shores(r, h, water, d['levels'], area_km2))
        out_cases[key] = dict(lakes=lakes, seas_receive_m3s={k: round(v) for k, v in supplies.items()},
                              all_seas_receive_m3s=round(sum(supplies.values())), **kept.get(key, {}))
        log(f'{key}: {lakes}, {time.time() - started:.0f} s')
    # the works' take-outs: is the lake each leaves from still full?
    take_outs = []
    for (cell, q, _), kind in zip(d['ops'], d['op_kind']):
        if q < 0.0:
            take_outs.append(dict(kind=kind, cell=cell, at=place(cell, cols, d['rotation']), flow_m3s=round(-q, 1),
                                  its_hollow_corrected=hollow_of[cell],
                                  **{key: dict(lake=bool(r['lake'][cell]), overflows=bool(r['lake_spills'][cell]),
                                               level_m=(None if not r['lake'][cell] else round(float(r['lake_level'][cell]) - d['level0'], 1)))
                                     for key, r in cases.items()},
                                  present_level_m=round(float(present['lake_level'][cell]) - d['level0'], 1)))
    waterways = []
    for x, e, lost in zip(record['waterways'], ends, on_course):
        row = dict(name=x.get('name'), anchor=x['anchor'], ends_in=x['ends_in'], fed_by=x.get('fed_by'),
                   mouth_loses_share=lost['lost_share'], lakes_passing_nothing_on_its_course=lost['lakes_passing_nothing_on'],
                   largest_fall_on_its_course=lost['largest_fall'])
        for end in ('mouth', 'head'):
            miss, c = e[end]
            row[end] = dict(present_m3s=round(float(present['discharge'][c])),
                            **{key + '_m3s': round(float(r['discharge'][c])) for key, r in cases.items()},
                            **{key + '_under_a_lake': bool(r['lake'][c]) for key, r in cases.items()},
                            recorded_m3s=x[f'flow_{end}_m3s'])
        waterways.append(row)
    BULK.mkdir(parents=True, exist_ok=True)
    r = cases['corrected']
    cells = np.flatnonzero(r['lake'] & (water < 0)).astype(np.int32)
    np.savez_compressed(PRODUCT, metadata=json.dumps(dict(schema=SCHEMA + '-grid', grid=dict(pixels_per_degree=dr.PPD, rows=rows, cols=cols),
                                                          common_level_m=d['level0'])),
                        discharge_m3s=r['discharge'], lake_cell=cells, lake_level_m=r['lake_level'][cells].astype(np.float32),
                        lake_spills=r['lake_spills'][cells])
    structure = SR.STRUCTURE
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={str(DRt.OUT.relative_to(ROOT)): sha256(DRt.OUT), str(WATERWAYS.relative_to(ROOT)): sha256(WATERWAYS)},
        ignored_inputs={str(structure): sha256(structure), str(DRt.PRODUCT): sha256(DRt.PRODUCT),
                        FD.CLIMATOLOGY: sha256(ROOT / FD.CLIMATOLOGY)},
        evidence=('A calculation on the design routing\'s own inputs: its floor, its seas at their design levels, the '
                  'design climate\'s closed budget and the works\' flows at their design values, routed once with a lake '
                  'rule that fills each hollow as one water, beside the present rule on the same inputs. The seas\' '
                  'levels, the weirs\' shares and the factor on the open water\'s evaporation are the design\'s, set on '
                  'the present rule; they are not set again here, so the figures show what the rule changes before '
                  'anything is balanced anew. Nothing is validated against an observation.'),
        reading_rule=('cases: present is route_levels (the design routing, reproduced within '
                      'control_largest_difference_m3s); corrected fills each hollow as one water, with dry ground '
                      'yielding its runoff wherever it lies; corrected_without_floor_runoff keeps the present rule\'s '
                      'way of yielding no runoff below a hollow\'s outer pass, so the difference between the two '
                      'corrected cases is that convention\'s. lakes are off the seas; shore_cells_not_held and '
                      'ground_below_and_joined_km2 are lake_support.py\'s measures. seas_receive_m3s is what the rivers '
                      'and overflows bring each sea, to set against seas_lose_m3s (the open water\'s loss beyond the rain '
                      'on it; the near-side sea\'s is the whole sea\'s on the design state, before it was parted). '
                      'budget_m3s: dry ground\'s runoff, what is put in and the lakes\' rain less evaporation, against '
                      'what reaches the seas and what lakes not full keep beyond their last whole cell. waterways: each '
                      'named waterway\'s flow at its mouth and head cells on the design routing, by case, and whether '
                      'the cell lies under a lake in the corrected cases; for a waterway whose mouth loses more than '
                      'lost_share of its flow, the lakes of on_a_course_km2 or more that pass nothing on along its '
                      'design course from head to mouth, and the place between dry cells of that course where the '
                      'corrected flow falls most beyond the design\'s own fall there, with the water above it followed '
                      'lake by lake (lake_rule.trace: the first trace_listed lakes it enters, each with what it passes '
                      'on, and how it ends). corrected_closer: the corrected lakes as bodies (a hollow\'s '
                      'one water), the largest of all and the largest that pass nothing on with their place on today\'s Moon and in the settled '
                      'frame, the lakes\' area by 30-degree band of today\'s latitude, the ground newly and no longer '
                      'under a lake, and the lakes\' volume (water level less ground, over the lakes\' cells). take_outs: the lake each of the works\' branches leaves from, by case (level above '
                      'the common level), and in the corrected case the hollow it lies in: its outer pass, whether it '
                      'stands full, what it overflows, and the flow it is short of standing full.'),
        control_largest_difference_m3s=differs,
        seas_lose_m3s={k: round(v) for k, v in d['loss'].items()},
        parameters=dict(closed_listed=CLOSED_LISTED, on_a_course_km2=ON_A_COURSE_KM2, lost_share=LOST_SHARE,
                        trace_listed=TRACE_LISTED),
        cases=out_cases, corrected_closer=closer, waterways=waterways, take_outs=take_outs,
        product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)), seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    log = lambda *a: print(*a, flush=True)                                     # noqa: E731
    record = build(log)
    OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    log(f'wrote {OUT} in {record["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
