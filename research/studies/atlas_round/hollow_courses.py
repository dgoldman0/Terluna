"""Courses through the lakes of a routing made by the corrected lake rule: each lowest point given the cell its water
goes on to, so that following the receivers from any cell reaches a sea, with the flow and the drained area gathered
along them.

lake_rule.py's routing keeps each cell's steepest-descent neighbour, and a lowest point is its own receiver: a lake's
overflow is a flow added where it lands. The maps that read the design routing follow a river down through its lakes
to its sea, find the waters a region drains to and trace a river up its largest donor, which needs the receivers
joined across each pass as research/studies/drainage_network/routing.py's are. connect() joins them from the
hollows' tree and their state:

- a hollow that overflows towards a sea sends its water to the cell its overflow lands on;
- two hollows that stand full to the sill between them are one lake: the one whose side the larger hollow's pass
  lies on carries the water on, and the other's goes over the sill into it;
- a full hollow beside one that is not full sends its water over the sill into that one;
- a hollow that is not full keeps its water. Its lowest point is still given the cell its water would go on to if it
  filled, so that every course reaches a sea; nothing flows along that step.

A lake of several lowest points is crossed from one to the next and left at the lowest point of the catchment its
pass lies in. The flow gathered along the joined receivers is the routing's own on dry ground. Under a lake it is
what has come in so far on the way to the lake's exit, and at the exit's lowest point it is what the lake passes on.
Two full hollows joined at their sill are one lake, yet the sill's own cell can stand dry between them: the water
that has come in is carried across it (`crossing` marks such cells), so a river's course is unbroken through the
lake. The drained area passes a lowest point where water does.
"""
from __future__ import annotations

import numpy as np

NEIGHBOURS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))


def roots_of(parent):
    """The hollow that overflows towards a sea round each hollow."""
    root = np.where(parent >= 0, parent, np.arange(parent.size))
    while True:
        nxt = root[root]
        if np.array_equal(nxt, root):
            return root
        root = nxt


def exit_cells(s, shape):
    """The cell on its own side of the pass each hollow overflows at: for a hollow joined to another at a pass, the
    cell the other's overflow lands on; for one that overflows towards a sea, the cell of its ground beside the
    cell its overflow lands on across the lowest crossing."""
    tree, label, hf = s['tree'], s['label'], s['hf']
    parent, sibling, land = tree['parent'], tree['sibling'], tree['land']
    n_out = tree['n_out']
    rows, cols = shape
    out = np.full(parent.size, -1, np.int64)
    joined = np.flatnonzero(parent >= 0)
    out[joined] = land[sibling[joined]]
    root = roots_of(parent)
    tops = np.flatnonzero((parent < 0) & (np.arange(parent.size) >= n_out))
    onto = land[tops]
    best = np.full(tops.size, np.inf)
    i, j = onto // cols, onto % cols
    for di, dj in NEIGHBOURS:
        ii, jj = i + di, (j + dj) % cols
        ok = (ii >= 0) & (ii < rows)
        c = np.where(ok, ii, 0) * cols + jj
        mine = ok & (label[c] >= n_out) & (root[np.maximum(label[c], 0)] == tops)
        crossing = np.where(mine, np.maximum(hf[c], hf[onto]) + 1e-9 * hf[c], np.inf)
        better = crossing < best
        best[better] = crossing[better]
        out[tops[better]] = c[better]
    if (out[tops] < 0).any():
        raise RuntimeError('a hollow\'s overflow lands on a cell that none of its ground lies beside')
    return out


def connect(s, terms, state, r, shape):
    """The joined receivers of a routing (lake_rule.hollow_structure, hollow_terms, hollow_state and hollow_routing's
    returns) on a grid of `shape`, with the flow and the drained area gathered along them. Returns receiver,
    discharge (m3/s), drained_m2, passes (whether water goes on from each cell) and crossing (the dry cells of a sill
    that a lake's water is carried across)."""
    tree = s['tree']
    parent, sibling, kids, land = tree['parent'], tree['sibling'], tree['kids'], tree['land']
    n_out, n_nodes = tree['n_out'], tree['n_nodes']
    n_all = parent.size
    label, sinks, outlet = s['label'], s['sinks'], s['outlet']
    full, out, over = r['full'], r['out'], r['over']
    # each hollow's catchments as a run of places, so that a catchment's hollow among two is told at once
    size = np.ones(n_all, np.int64)
    for x in range(n_nodes, n_all):
        size[x] = size[kids[x, 0]] + size[kids[x, 1]]
    lo = np.zeros(n_all, np.int64)
    tops = np.flatnonzero((parent < 0) & (np.arange(n_all) >= n_out))
    lo[tops] = np.concatenate([[0], np.cumsum(size[tops])[:-1]])
    for x in range(n_all - 1, n_nodes - 1, -1):
        a, b = kids[x]
        lo[a], lo[b] = lo[x], lo[x] + size[a]
    exits = exit_cells(s, shape)
    # where each hollow's water goes on to, and the cell on its own side it leaves by: its own pass's where it goes
    # over its own pass, the larger hollow's where it carries the larger hollow's water on
    on = np.full(n_all, -1, np.int64)
    leaves_by = np.full(n_all, -1, np.int64)
    on[tops], leaves_by[tops] = land[tops], exits[tops]
    lo_, size_, label_ = lo.tolist(), size.tolist(), label
    for x in range(n_all - 1, n_nodes - 1, -1):
        a, b = int(kids[x, 0]), int(kids[x, 1])
        at = lo_[int(label_[leaves_by[x]])]
        if not lo_[x] <= at < lo_[x] + size_[x]:
            raise RuntimeError('a hollow\'s water leaves by a cell that is not its own')
        c, o = (a, b) if at < lo_[a] + size_[a] else (b, a)      # c holds the cell the hollow's water leaves by
        if full[a] != full[b]:                           # the full one spills over the sill into the other, which
            y, z = (a, b) if full[a] else (b, a)         # keeps its water: its course goes on as the larger one's
            on[y], leaves_by[y] = land[y], exits[y]
            on[z], leaves_by[z] = on[x], (leaves_by[x] if z == c else exits[z])
        else:                                            # one lake to the sill, or neither full
            on[c], leaves_by[c] = on[x], leaves_by[x]
            on[o], leaves_by[o] = land[o], exits[o]
    receiver = s['receiver'].astype(np.int64).copy()
    receiver[sinks] = on[n_out:n_nodes]
    # what each lowest point passes on: its lake's outflow at the lake's exit, what has come in so far elsewhere in
    # a lake of several, nothing from a lake that is not full
    top = np.where(full & (parent >= 0) & full[np.maximum(sibling, 0)], parent, np.arange(n_all))
    while True:
        nxt = top[top]
        if np.array_equal(nxt, top):
            break
        top = nxt
    leaves = np.arange(n_out, n_nodes)
    body = top[leaves]
    at_exit = label[leaves_by[body]] == leaves
    spills = full[body]
    outflow = np.where(spills, np.where(parent[body] < 0, out[body], over[body]), 0.0)
    through = ~at_exit                                   # carried on unchanged inside a lake of several lowest points
    n = outlet.size
    pit = np.full(n, -1, np.int64)
    pit[sinks] = np.arange(sinks.size)
    flow = np.where(r['lake'] | outlet, 0.0, terms['run']) + state['put']
    acc = flow.astype(np.float64).copy()
    drained = s['area'] * np.where(outlet, 0.0, 1.0)
    passes = np.ones(n, bool)
    donor = ~outlet
    indeg = np.bincount(receiver[donor], minlength=n).astype(np.int32)
    front = np.flatnonzero((indeg == 0) & donor)
    done = 0
    while front.size:
        done += front.size
        send, basin = acc[front].copy(), drained[front].copy()
        k = pit[front]
        here = k >= 0
        if here.any():
            kk = k[here]
            send[here] = np.where(through[kk], send[here], outflow[kk])
            leaving = ~through[kk] & spills[kk]          # a lake's exit carries what the lake passes on
            acc[front[here][leaving]] = send[here][leaving]
            goes = send[here] > 0.0
            basin[here] *= goes
            passes[front[here]] = goes
        to = receiver[front]
        np.add.at(acc, to, send)
        np.add.at(drained, to, basin)
        np.subtract.at(indeg, to, 1)
        cand = np.unique(to)
        front = cand[(indeg[cand] == 0) & donor[cand]]
    if done != int(donor.sum()):
        raise RuntimeError('the joined receivers do not all lead to a sea')
    # the dry cells of a sill that a lake's water crosses between two of its lowest points
    crossing = np.zeros(n, bool)
    wet = r['lake'] | outlet
    for cell in on[n_out:n_nodes][through & (acc[sinks] > 0.0)].tolist():
        while not wet[cell]:
            crossing[cell] = True
            cell = int(receiver[cell])
    return dict(receiver=receiver, discharge=acc, drained_m2=drained, passes=passes, crossing=crossing)
