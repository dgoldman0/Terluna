"""Geography's routing with seas standing at different levels.

geography.drainage.route takes every outlet as one water at the highest outlet cell's height, which suits the atlas,
where every sea stands at the common level. The drainage network study needs seas at their own levels at once: a
hollow on the slope of a sea that has fallen must spill towards that sea's level, not another's. This is
geography.drainage.route with that one change: each outlet cell belongs to the water of its own level, each such
water is a root of the flood at its level, and the roots keep their levels.

The routing's other pieces (d8_receivers, terminal, passes, the fill-and-spill balance) are geography's own, called
from geography.drainage. The change stays in this study for now because six products elsewhere bind the hash of
geography/drainage.py (biosphere's terrain axis, the lake census, the south sea study's first stage, the storm-water
study, the terrain-through-time routing and resources' sea phosphorus), and the atlas round, which rebuilds
geography's drainage and those products, is the place to fold it in. test_network.py checks that with one level it
routes exactly as geography's routing does.
"""
from __future__ import annotations

import heapq

import numpy as np

from geography.drainage import MM_DAY, d8_receivers, passes, terminal


def flood_from(n_nodes, edges, roots):
    """Each catchment's spill level, the minimax pass height on its way to the roots, with the pass it spills
    through (the cell on its own side and the cell on its parent's side). roots: (node, level) pairs; each root
    keeps its level, and a catchment takes the lowest minimax over all of them."""
    u, v, w, cu, cv = edges
    src = np.concatenate([u, v]); dst = np.concatenate([v, u]); wt = np.concatenate([w, w])
    own = np.concatenate([cv, cu]); other = np.concatenate([cu, cv])     # cells on dst's and src's side
    order = np.argsort(src, kind='stable')
    src, dst, wt, own, other = src[order], dst[order], wt[order], own[order], other[order]
    start = np.searchsorted(src, np.arange(n_nodes + 1)).tolist()
    dst, wt, own, other = dst.tolist(), wt.tolist(), own.tolist(), other.tolist()
    level = np.full(n_nodes, np.nan)
    spill_from = np.full(n_nodes, -1, np.int64)
    spill_to = np.full(n_nodes, -1, np.int64)
    done = bytearray(n_nodes)
    heap = []
    for k, lev in roots:
        k, lev = int(k), float(lev)
        done[k] = 1
        level[k] = lev
        for e in range(start[k], start[k + 1]):
            if not done[dst[e]]:
                heap.append((max(wt[e], lev), dst[e], own[e], other[e]))
    heapq.heapify(heap)
    while heap:
        lev, x, c_own, c_other = heapq.heappop(heap)
        if done[x]:
            continue
        done[x] = 1
        level[x], spill_from[x], spill_to[x] = lev, c_own, c_other
        for k in range(start[x], start[x + 1]):
            y = dst[k]
            if not done[y]:
                heapq.heappush(heap, (max(wt[k], lev), y, own[k], other[k]))
    return level, spill_from, spill_to


def route_levels(h, outlet, lat_deg, cell_km, runoff_mm_day, lake_net_mm_day, outlet_level):
    """geography.drainage.route with each outlet cell standing at the level of its own water (outlet_level, an
    array of h's shape): outlets of one level form one outlet node, a root of the flood at that level. Returns the
    same dict as geography's routing, with outlet_nodes, the number of outlet nodes that lead the node tables."""
    rows, cols = h.shape
    n = rows * cols
    idx2 = np.arange(n, dtype=np.int32).reshape(rows, cols)
    receiver = d8_receivers(h.astype(np.float32), outlet, lat_deg, cell_km)
    flat_out = outlet.ravel()
    idx = idx2.ravel()
    sink = (receiver == idx) & ~flat_out
    sinks = np.flatnonzero(sink).astype(np.int32)
    node_of = np.full(n, -1, np.int32)
    levels, groups = np.unique(np.asarray(outlet_level, dtype=float).ravel()[flat_out], return_inverse=True)
    n_out = max(1, levels.size)
    node_of[flat_out] = groups
    node_of[sinks] = np.arange(n_out, sinks.size + n_out, dtype=np.int32)
    label = node_of[terminal(receiver)]
    n_nodes = sinks.size + n_out
    h_pass = np.where(outlet, -np.inf, h).astype(np.float32)
    edges = passes(label.reshape(rows, cols), h_pass, idx2)
    level, spill_from, spill_to = flood_from(n_nodes, edges, list(enumerate(levels.tolist())))
    if np.isnan(level).any():
        raise RuntimeError(f'{int(np.isnan(level).sum())} catchments have no path to an outlet')

    hf = h.ravel()
    area = ((cell_km * 1000.0) ** 2 * np.maximum(np.cos(np.radians(lat_deg)), 1e-6))[:, None] * np.ones((1, cols))
    area = area.ravel()
    lake_full = ~flat_out & (hf < level[label])
    runoff = np.where(flat_out | lake_full, 0.0, np.maximum(runoff_mm_day.ravel(), 0.0) * MM_DAY * area)
    lake_net = np.bincount(label[lake_full], weights=(lake_net_mm_day.ravel() * MM_DAY * area)[lake_full],
                           minlength=n_nodes)

    final = receiver.copy()
    final[sinks] = spill_to[n_out:]
    donor = ~flat_out
    indeg = np.bincount(final[donor], minlength=n).astype(np.int32)
    acc = runoff.copy()
    drained = area * np.where(flat_out, 0.0, 1.0)
    inflow = np.zeros(n_nodes)
    passed_on = np.zeros(n_nodes)
    front = np.flatnonzero((indeg == 0) & donor)
    while front.size:
        out = acc[front].copy()
        basin = drained[front].copy()
        at_sink = sink[front]
        if at_sink.any():
            k = node_of[front[at_sink]]
            inflow[k] = out[at_sink]
            out[at_sink] = np.maximum(out[at_sink] + lake_net[k], 0.0)
            passed_on[k] = out[at_sink]
            basin[at_sink] *= out[at_sink] > 0.0
        r = final[front]
        np.add.at(acc, r, out)
        np.add.at(drained, r, basin)
        np.subtract.at(indeg, r, 1)
        cand = np.unique(r)
        front = cand[(indeg[cand] == 0) & donor[cand]]

    # Depressions that cannot overflow keep the lake whose evaporation balances their inflow.
    closed = np.zeros(n_nodes, bool)
    closed[n_out:] = (inflow[n_out:] + lake_net[n_out:] < 0.0) & (lake_net[n_out:] < 0.0)
    loss = np.bincount(label[lake_full], weights=(-lake_net_mm_day.ravel() * MM_DAY * area)[lake_full], minlength=n_nodes)
    lake = lake_full.copy()
    cells = np.flatnonzero(lake_full & closed[label])
    if cells.size:
        order = np.lexsort((hf[cells], label[cells]))
        cells = cells[order]
        lab = label[cells]
        rate = np.where(loss[lab] > 0, (-lake_net_mm_day.ravel() * MM_DAY * area)[cells], 0.0)
        cum = np.cumsum(rate)
        group_start = np.r_[0, np.flatnonzero(lab[1:] != lab[:-1]) + 1]
        base = np.repeat(cum[group_start] - rate[group_start], np.diff(np.r_[group_start, cells.size]))
        lake[cells] = (cum - base) <= inflow[lab]
    lake_level = np.where(lake, level[label], np.nan)
    if cells.size:
        top = np.full(n_nodes, -np.inf)
        wet = cells[lake[cells]]
        np.maximum.at(top, label[wet], hf[wet])
        lake_level[cells] = np.where(lake[cells], top[label[cells]], np.nan)
    return dict(receiver=final, label=label, level=level, discharge=acc, drained_m2=drained, lake=lake,
                lake_level=lake_level, closed=closed, inflow=inflow, passed_on=passed_on, lake_net=lake_net,
                sinks=sinks, spill_from=spill_from, spill_to=spill_to, area=area, outlet=flat_out, outlet_nodes=n_out)
