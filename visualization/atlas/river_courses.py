"""The rivers' courses for the globe: smooth lines that taper with the flow for the far view, and water laid in the
valleys of the 474 m ground for the close views.

    python3 visualization/atlas/river_courses.py --waterways research/studies/atlas_round/results/waterways.json \\
        --surface <globe_surface_16ppd.npz> --out <folder>          # river_lines.json, seconds per river

The page's first lines (render.river_marks) joined the centres of the routing's 1.9 km cells and crossed every lake
straight, at one width. Here each named waterway (waterways.py) is followed again on the design routing from its mouth
up its main donors, with the flow at every cell, and parted into channel reaches and lake crossings:

- a channel reach runs from shore to shore (the midpoints between its end cells and the water beside them), smoothed
  on the sphere over SMOOTH_KM, and carries the flow at every point, so that its width can follow the square root of
  the flow as the page's width rule does (appearance.RIVER_WIDTH);
- a lake crossing follows the lake's deep water from where the river enters to where it leaves (the least-cost way
  over the lake's cells, shallow water costing more), to be drawn faintly: the lake itself is the water there.

How wide the page draws a reach is a display rule (display_half_width_px): the computed width where it spans more
than the floor, and otherwise a floor in pixels that grows with the flow, which is an exaggeration at every zoom
where the water is narrower than a pixel. The patch for the page draws these lines.

For the close views, fine_rivers lays every river of FINE_FLOW or more on a window of the 474 m ground: the routed
1.9 km cells fix which cells a river crosses and its flow; within them the course takes the least-cost way over the
fine heights, low ground costing least, joined from cell to receiver only, so that no course leaves its routed cells;
tributaries end on the river they join; each course is smoothed and drawn at its computed width as the share of each
pixel its water covers: its area in a pixel is its length there times its width, so a width under a pixel is a share
of the pixel.
"""
from __future__ import annotations

import argparse
import heapq
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import appearance  # noqa: E402

MOON_KM = appearance.MOON_KM
HEAD_M3S = 100.0                      # a waterway is followed upstream while it carries this, as waterways.py does
SMOOTH_KM, SAMPLE_KM, PASSES = 9.0, 1.0, 3   # the far lines: a moving average over this span, this many times
LAKE_SHALLOW_COST, LAKE_DEPTH_M = 4.0, 20.0  # a lake cell costs 1 + this / (1 + depth / LAKE_DEPTH_M) to cross
# Display rules for the far view (exaggerations where the water is narrower): the least half width drawn, in pixels,
# at HEAD_M3S and at FLOOR_FLOW, growing between them with the square root of the flow.
FLOOR_HALF_PX, FLOOR_FLOW = (0.3, 1.25), 30000.0
FINE_FLOW = 30.0                      # the close views draw rivers from this flow up (m3/s)
K = 4                                 # fine cells along a routed cell's side
VALLEY_COST, RELIEF_M = 4.0, 5.0      # a fine cell costs 1 + VALLEY_COST * (its height above its block's lowest) /
#                                       (the block's relief + RELIEF_M): flat blocks cost alike, so courses run straight
SMOOTH_PX = 2.0                       # the fine courses' smoothing along their length (fine cells, about 1 km)
STEP, ACROSS = 0.05, 0.25             # the fine rivers' water is laid down at points this far apart along a course
#                                       and across it (fine cells)


def unit(lat, lon):
    la, lo = np.radians(np.asarray(lat, float)), np.radians(np.asarray(lon, float))
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def lat_lon(u):
    u = np.asarray(u, float)
    return np.degrees(np.arcsin(np.clip(u[..., 2], -1.0, 1.0))), np.degrees(np.arctan2(u[..., 1], u[..., 0]))


def routing_arrays(path):
    """The design routing as the courses read it: receivers, flows, the cells under water and their levels."""
    with np.load(path, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        rows, cols = meta['grid']['rows'], meta['grid']['cols']
        receiver = z['receiver'].astype(np.int64)
        q = z['discharge_m3s'].astype(np.float32)
        drained = z['drained_km2']
        sea = z['sea']
        water = z['routing_lake'] | (sea >= 0)
        water[z['lake_cell']] = True
        level = np.full(rows * cols, np.nan, np.float32)
        level[z['lake_cell']] = z['lake_level_m']
        for k, name in enumerate(meta['seas']):
            level[sea == k] = meta['levels_m'][name]
    return dict(rows=rows, cols=cols, receiver=receiver, q=q, drained=drained, water=water, level=level, meta=meta)


def main_donors(receiver, q):
    """For every cell, the cell draining into it that carries the most water, or -1 (as waterways.main_donors)."""
    order = np.argsort(q, kind='stable')
    order = order[(receiver[order] >= 0) & (receiver[order] != order)]
    donor = np.full(receiver.size, -1, np.int64)
    donor[receiver[order]] = order
    return donor


def follow(mouth, head, donor, q, least=HEAD_M3S):
    """A waterway's cells from its head down to its mouth: from the mouth up the main donors to the recorded head,
    or to where the flow falls under `least`."""
    cells = [int(mouth)]
    while cells[-1] != head:
        up = int(donor[cells[-1]])
        if up < 0 or q[up] < least:
            break
        cells.append(up)
    return np.asarray(cells[::-1], np.int64)


def pieces(wet):
    """Runs of channel and of lake along a course, as (kind, first, last) with both ends counted."""
    out, i = [], 0
    while i < len(wet):
        j = i
        while j + 1 < len(wet) and wet[j + 1] == wet[i]:
            j += 1
        out.append(('lake' if wet[i] else 'channel', i, j))
        i = j + 1
    return out


def smooth_line(u, values=None, window_km=SMOOTH_KM, sample_km=SAMPLE_KM, passes=PASSES):
    """A chain of unit vectors made a smooth line on the sphere: resampled every sample_km, then a moving average
    over window_km, `passes` times, its two ends kept. Returns the points, their distance along the line (km) and the
    values carried to them."""
    u = np.asarray(u, float)
    step = np.arccos(np.clip((u[1:] * u[:-1]).sum(axis=1), -1.0, 1.0)) * MOON_KM
    along = np.r_[0.0, np.cumsum(step)]
    n = max(int(np.ceil(along[-1] / sample_km)), 1) + 1
    t = np.linspace(0.0, along[-1], n)
    p = np.stack([np.interp(t, along, u[:, k]) for k in range(3)], axis=1)
    half = min(max(int(round(window_km / sample_km / 2.0)), 1), max((n - 1) // 2, 0))
    if half:
        kernel = np.ones(2 * half + 1) / (2 * half + 1)
        for _ in range(passes):
            padded = np.vstack([np.repeat(p[:1], half, axis=0), p, np.repeat(p[-1:], half, axis=0)])
            s = np.stack([np.convolve(padded[:, k], kernel, mode='valid') for k in range(3)], axis=1)
            s[0], s[-1] = p[0], p[-1]
            p = s
    p /= np.linalg.norm(p, axis=1, keepdims=True)
    return p, t, (np.interp(t, along, np.asarray(values, float)) if values is not None else None)


def lake_path(a, b, water, depth, rows, cols, margin=40):
    """The cells from a to b over water, favouring deep water: the least-cost way through the water cells of the box
    round the two (with a margin), a cell costing 1 + LAKE_SHALLOW_COST / (1 + depth / LAKE_DEPTH_M). Falls back to
    the two cells alone where the box's water does not join them."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    ra, ca = divmod(int(a), cols)
    rb, cb = divmod(int(b), cols)
    db = (cb - ca + cols // 2) % cols - cols // 2
    r0, r1 = max(min(ra, rb) - margin, 0), min(max(ra, rb) + margin, rows - 1)
    c0, c1 = min(0, db) - margin, max(0, db) + margin
    idx = np.arange(r0, r1 + 1)[:, None] * cols + ((ca + np.arange(c0, c1 + 1)) % cols)[None, :]
    wet = water[idx]
    cost = 1.0 + LAKE_SHALLOW_COST / (1.0 + np.clip(np.nan_to_num(depth[idx], nan=0.0), 0.0, None) / LAKE_DEPTH_M)
    nr, nc = idx.shape
    pid = np.arange(nr * nc).reshape(nr, nc)
    I, J, W = [], [], []
    for dr, dc, length in ((0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.414), (1, -1, 1.414)):
        p = (slice(0, nr - dr), slice(max(0, -dc), nc - max(0, dc)))
        s = (slice(dr, nr), slice(max(0, dc), nc + min(0, dc)))
        ok = wet[p] & wet[s]
        I.append(pid[p][ok])
        J.append(pid[s][ok])
        W.append((0.5 * length * (cost[p] + cost[s]))[ok])
    graph = coo_matrix((np.concatenate(W), (np.concatenate(I), np.concatenate(J))), shape=(nr * nc, nr * nc)).tocsr()
    start, end = (ra - r0) * nc + (0 - c0), (rb - r0) * nc + (db - c0)
    dist, pred = dijkstra(graph, directed=False, indices=start, return_predecessors=True)
    if not np.isfinite(dist[end]):
        return np.asarray([a, b], np.int64)
    out = [end]
    while out[-1] != start:
        out.append(int(pred[out[-1]]))
    return idx.ravel()[np.asarray(out[::-1])]


def cell_units(cells, cols, ppd=16):
    cells = np.asarray(cells)
    return unit(90.0 - (cells // cols + 0.5) / ppd, (cells % cols + 0.5) / ppd)


def river_lines(waterways, routing, height, rotation, log=None):
    """Each named waterway as the page draws it: its pieces from the head down to the mouth, a channel reach as
    points with the flow at each (m3/s), a lake crossing as the line through its deep water, all in the settled
    frame. `height` is the design surface's (flat, on the routing's grid), for the lakes' depths."""
    rows, cols, q, water, receiver = (routing[k] for k in ('rows', 'cols', 'q', 'water', 'receiver'))
    donor = main_donors(receiver, q)
    depth = routing['level'] - np.asarray(height, np.float32).ravel()
    at = lambda place: int((90.0 - place[0]) * 16) * cols + int((place[1] % 360.0) * 16)   # noqa: E731
    mid = lambda a, b: (lambda m: m / np.linalg.norm(m))(cell_units([a], cols)[0] + cell_units([b], cols)[0])   # noqa: E731
    turn = lambda u: np.asarray(u) @ np.asarray(rotation, float).T                          # noqa: E731
    out = []
    for x in waterways:
        if not x.get('name'):
            continue
        cells = follow(at(x['mouth']['today']), at(x['head']['today']), donor, q)
        runs = pieces(water[cells])
        parts = []
        for kind, i, j in runs:
            run = cells[i:j + 1]
            if kind == 'lake' and j > i:
                run = lake_path(run[0], run[-1], water, depth, rows, cols)
            u = cell_units(run, cols)
            flow = q[run].astype(float)
            # a reach ends at the shores: the midpoints to the cells before and after it along the course
            before = cells[i - 1] if i > 0 else None
            after = cells[j + 1] if j + 1 < len(cells) else (int(receiver[cells[j]]) if kind == 'channel' else None)
            if before is not None:
                u, flow = np.vstack([mid(before, run[0])[None], u]), np.r_[flow[0], flow]
            if after is not None and after >= 0 and after != cells[j]:
                u, flow = np.vstack([u, mid(run[-1], after)[None]]), np.r_[flow, flow[-1]]
            if len(u) < 2:
                continue
            p, km, f = smooth_line(u, flow)
            keep = np.unique(np.r_[np.arange(0, len(p), 2), len(p) - 1])     # a point every 2 km
            la, lo = lat_lon(turn(p[keep]))
            parts.append(dict(kind=kind, km=round(float(km[-1]), 1),
                              points=[[round(float(a), 4), round(float(b), 4), round(float(c), 1)]
                                      for a, b, c in zip(la, lo, f[keep])]))
        out.append(dict(name=x['name'], anchor=x['anchor'], flow=x['flow_mouth_m3s'], pieces=parts))
        if log:
            log(f"{x['name']}: {len(cells)} cells, {sum(p['kind'] == 'channel' for p in parts)} reaches, "
                f"{sum(p['kind'] == 'lake' for p in parts)} lake crossings")
    return out


def display_half_width_px(flow, km_per_px):
    """Half the width a channel is drawn at, in pixels: the computed width (appearance.RIVER_WIDTH times the square
    root of the flow) at the view's scale, and at least a floor that grows with the square root of the flow from
    FLOOR_HALF_PX[0] at HEAD_M3S to FLOOR_HALF_PX[1] at FLOOR_FLOW. Wherever the floor is the larger, the line is
    wider than the water: a display exaggeration. Returns the half width and the exaggeration factor."""
    flow = np.asarray(flow, float)
    true_half = 0.5 * appearance.RIVER_WIDTH * np.sqrt(flow) / 1000.0 / km_per_px
    s = np.clip((np.sqrt(flow) - np.sqrt(HEAD_M3S)) / (np.sqrt(FLOOR_FLOW) - np.sqrt(HEAD_M3S)), 0.0, 1.0)
    floor = FLOOR_HALF_PX[0] + (FLOOR_HALF_PX[1] - FLOOR_HALF_PX[0]) * s
    half = np.maximum(true_half, floor)
    return half, half / np.maximum(true_half, 1e-9)


# ----------------------------------------------------------------------------------------- the close views

def fine_rivers(h, wet, lat, r_first, c_first, routing, sediment=None, least=FINE_FLOW):
    """The rivers of `least` or more over a window of the fine ground. h (m) and wet are the window's fine heights and
    water, lat its rows' latitudes, r_first and c_first its first fine row and column on the globe (multiples of K, as
    its sides are); routing is routing_arrays' and sediment the flow each routed cell carries off land since the last
    lake (appearance.land_discharge), where given. Returns the share of each fine pixel river water covers, the
    largest river's flow through each pixel, the sediment-bearing share of that flow, and the courses."""
    from scipy.ndimage import gaussian_filter1d
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    R, C = h.shape
    if r_first % K or c_first % K or R % K or C % K:
        raise ValueError('the window must be made of whole routed cells')
    cols16 = routing['cols']
    nr, nc = R // K, C // K
    i0, j0 = r_first // K, c_first // K
    cell = (i0 + np.arange(nr))[:, None] * cols16 + ((j0 + np.arange(nc)) % cols16)[None, :]
    q = routing['q'][cell].ravel().astype(np.float64)
    water16 = routing['water'][cell].ravel()
    silt = (np.clip(np.asarray(sediment)[cell].ravel() / np.maximum(q, 1e-9), 0.0, 1.0) if sediment is not None
            else np.zeros(nr * nc))
    rcv_g = routing['receiver'][cell]
    ri, rj = rcv_g // cols16 - i0, (rcv_g % cols16 - j0) % cols16
    inside = (rcv_g >= 0) & (rcv_g != cell) & (ri >= 0) & (ri < nr) & (rj < nc)
    rcv = np.where(inside, ri * nc + rj, -1).ravel()
    has = rcv >= 0
    safe = np.where(has, rcv, 0)
    river = (q >= least) & ~water16
    terminal = np.zeros(nr * nc, bool)
    terminal[rcv[river & has]] = True
    terminal &= ~river
    upper = water16 & has & river[safe]                    # a lake's cell that passes its water to a river
    empty = dict(share=np.zeros(h.shape, np.float32), flow=np.zeros(h.shape, np.float32),
                 sediment=np.zeros(h.shape, np.float32), courses=[], unplaced=0)
    if not river.any():
        return empty

    # the fine graph: pixels of the rivers' cells and of the cells they end in, joined both ways within a cell and one
    # way from a cell to its receiver, so that a course can only run down its routed cells (a lake that takes one
    # river in and sends another out is the end of the first and never a way back up the second)
    cid = (np.arange(R)[:, None] // K) * nc + (np.arange(C)[None, :] // K)
    node = (river | terminal)[cid]
    blocks = h.reshape(nr, K, nc, K)
    low = np.repeat(np.repeat(blocks.min(axis=(1, 3)), K, axis=0), K, axis=1)
    high = np.repeat(np.repeat(blocks.max(axis=(1, 3)), K, axis=0), K, axis=1)
    cost = 1.0 + VALLEY_COST * (h - low) / (high - low + RELIEF_M)
    cost[wet] = 1.0
    pid = np.arange(R * C).reshape(R, C)
    I, J, W = [], [], []
    for dr, dc, length in ((0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.414), (1, -1, 1.414)):
        p = (slice(0, R - dr), slice(max(0, -dc), C - max(0, dc)))
        s = (slice(dr, R), slice(max(0, dc), C + min(0, dc)))
        a, b = cid[p], cid[s]
        both = node[p] & node[s]
        w = 0.5 * length * (cost[p] + cost[s])
        down = both & ((a == b) | (rcv[a] == b))           # water may pass from p's pixel to s's
        up = both & ((a == b) | (rcv[b] == a))             # and from s's to p's
        I += [pid[p][down], pid[s][up]]
        J += [pid[s][down], pid[p][up]]
        W += [w[down], w[up]]
    # searched from the ends upstream: the graph's edges turned round
    graph = coo_matrix((np.concatenate(W), (np.concatenate(J), np.concatenate(I))), shape=(R * C, R * C)).tocsr()
    wet_block = wet.reshape(nr, K, nc, K).any(axis=(1, 3)).ravel()
    leaves = river & ~has                                  # a river that leaves the window or ends on dry ground
    source = leaves[cid] | (terminal[cid] & (wet | ~wet_block[cid]))
    dist, pred, _ = dijkstra(graph, directed=True, indices=pid[source], return_predecessors=True, min_only=True)

    # the courses: the largest first, each from its head down to its end or to the river it joins
    donors = {}
    for d in np.flatnonzero((river | upper) & has).tolist():
        donors.setdefault(int(rcv[d]), []).append(d)
    heap = [(-q[c], int(c)) for c in np.flatnonzero(river & ~(has & river[safe]))]
    heapq.heapify(heap)
    placed, courses, unplaced = {}, [], 0
    hf, cidf = h.ravel(), cid.ravel()
    while heap:
        _, c = heapq.heappop(heap)
        chain, lake = [c], -1
        while True:
            ds = donors.get(chain[-1], [])
            if not ds:
                break
            best = max(ds, key=lambda d: q[d])
            for d in ds:
                if d != best and river[d]:
                    heapq.heappush(heap, (-q[d], d))
            if not river[best]:
                lake = best                                # the lake cell whose water the river carries off
                break
            chain.append(best)
        head = chain[-1]
        hr, hc = divmod(head, nc)
        block = pid[hr * K:(hr + 1) * K, hc * K:(hc + 1) * K].ravel()
        block = block[np.isfinite(dist[block])]
        if not block.size:
            unplaced += 1
            continue
        first = None
        start = int(block[np.argmin(hf[block])])           # a head on dry ground: its cell's lowest place
        if lake >= 0:                                      # a head at a lake: where its cell meets the lake's, the cheapest way down
            lr, lc = divmod(lake, nc)
            by, bx = block // C, block % C
            cy, cx = np.clip(by, lr * K, lr * K + K - 1), np.clip(bx, lc * K, lc * K + K - 1)
            beside = np.maximum(np.abs(by - cy), np.abs(bx - cx)) <= 1
            if beside.any():
                k = int(np.flatnonzero(beside)[np.argmin(dist[block[beside]])])
                start, first = int(block[k]), int(cy[k] * C + cx[k])
        own, path = set(chain), [start]
        while True:
            nxt = int(pred[path[-1]])
            if nxt < 0:
                break
            cn = int(cidf[nxt])
            if cn in placed and cn not in own:             # the river it joins: end on that river's course
                pts = np.asarray(placed[cn])
                y, x = divmod(path[-1], C)
                path.append(int(pts[np.argmin((pts // C - y) ** 2 + (pts % C - x) ** 2)]))
                break
            path.append(nxt)
        for cn, pix in zip(cidf[path[:-1]].tolist(), path[:-1]):
            placed.setdefault(cn, []).append(pix)
        if int(cidf[path[-1]]) in own:
            placed.setdefault(int(cidf[path[-1]]), []).append(int(path[-1]))
        if first is not None:
            path = [first] + path
        if len(path) < 2:
            continue
        path = np.asarray(path)
        cells = cidf[path]
        cells[-1] = cells[-2] if cells[-1] not in own else cells[-1]
        y, x = (path // C).astype(float), (path % C).astype(float)
        if len(path) > 4:
            ys, xs = gaussian_filter1d(y, SMOOTH_PX, mode='nearest'), gaussian_filter1d(x, SMOOTH_PX, mode='nearest')
            ys[[0, -1]], xs[[0, -1]] = y[[0, -1]], x[[0, -1]]
            y, x = ys, xs
        # the flow at each point: its routed cell's, and beyond the river's own cells (the lake it leaves, the water
        # it ends in) the nearest river cell's along the course
        ok = river[cells]
        if not ok.any():
            continue
        last = np.maximum.accumulate(np.where(ok, np.arange(len(cells)), -1))
        cells = cells[np.where(last >= 0, last, int(np.argmax(ok)))]
        courses.append(dict(y=y, x=x, flow=q[cells], silt=silt[cells]))

    # laid on the grid: points every STEP of a fine cell along each course, and across its width where it is wider
    # than ACROSS, each carrying its piece of the water's area to the pixel it falls in, so that the area is the
    # course's length times its computed width whatever the width
    dy = MOON_KM * 1000.0 * np.radians(1.0 / (16 * K))
    area, silted = np.zeros(R * C), np.zeros(R * C)
    flow_map = np.zeros(R * C, np.float32)
    for k in courses:
        y, x = k['y'], k['x']
        seg = np.hypot(np.diff(x), np.diff(y))
        along = np.r_[0.0, np.cumsum(seg)]
        if along[-1] <= 0.0:
            continue
        n = int(np.ceil(along[-1] / STEP))
        t = (np.arange(n) + 0.5) * along[-1] / n
        xs, ys = np.interp(t, along, x), np.interp(t, along, y)
        fs, ss = np.interp(t, along, k['flow']), np.interp(t, along, k['silt'])
        i = np.clip(np.searchsorted(along, t) - 1, 0, len(seg) - 1)
        tx, ty = np.diff(x)[i] / np.maximum(seg[i], 1e-9), np.diff(y)[i] / np.maximum(seg[i], 1e-9)
        dxm = dy * np.cos(np.radians(np.interp(ys, np.arange(R), lat)))
        thick = appearance.RIVER_WIDTH * np.sqrt(fs) * np.hypot(tx * dxm, ty * dy) / (dxm * dy)   # fine cells across
        piece = thick * along[-1] / n                                                          # fine cells of water
        parts = np.maximum(np.ceil(thick / ACROSS).astype(int), 1)
        for m in np.unique(parts).tolist():
            sel = parts == m
            for j in range(m):
                off = ((j + 0.5) / m - 0.5) * thick[sel]
                px = np.rint(xs[sel] - ty[sel] * off).astype(np.int64)
                py = np.rint(ys[sel] + tx[sel] * off).astype(np.int64)
                ok = (px >= 0) & (px < C) & (py >= 0) & (py < R)
                at = py[ok] * C + px[ok]
                area += np.bincount(at, weights=piece[sel][ok] / m, minlength=R * C)
                silted += np.bincount(at, weights=(piece[sel] * ss[sel])[ok] / m, minlength=R * C)
                np.maximum.at(flow_map, at, fs[sel][ok].astype(np.float32))
    share = np.minimum(area, 1.0).astype(np.float32).reshape(R, C)
    silt_map = (silted / np.maximum(area, 1e-9)).astype(np.float32).reshape(R, C)
    flow_map = flow_map.reshape(R, C)
    return dict(share=share, flow=flow_map, sediment=silt_map, courses=courses, unplaced=unplaced)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--waterways', type=Path, required=True)
    ap.add_argument('--surface', type=Path, required=True, help='the globe surface product (its heights give the lakes\' depths)')
    ap.add_argument('--routing', type=Path, help='the design routing product; by default the one its record names')
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args(argv)
    started = time.time()
    if a.routing is None:
        from research.studies.atlas_round import design_routing as DRt
        a.routing = Path(json.loads(DRt.OUT.read_text())['product']['path'])
    routing = routing_arrays(a.routing)
    with np.load(a.surface, allow_pickle=False) as z:
        height = z['height_m']
    rotation = routing['meta']['frame']['rotation_today_to_turned']
    record = json.loads(a.waterways.read_text())
    lines = river_lines(record['waterways'], routing, height, rotation, log=lambda s: print(s, flush=True))
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / 'river_lines.json').write_text(json.dumps(dict(
        schema='terluna.visualization.globe-river-lines/1',
        note=('Each named waterway from its head down to its mouth in the settled frame: channel reaches as points '
              '[latitude, longitude, flow in m3/s], smoothed over %.0f km, from shore to shore; lake crossings along '
              'the lake\'s deep water, to be drawn faintly. The page widens a reach to a floor in pixels '
              '(display_half_width_px): an exaggeration wherever the water is narrower.' % SMOOTH_KM),
        settings=dict(smooth_km=SMOOTH_KM, floor_half_px=FLOOR_HALF_PX, floor_flow_m3s=FLOOR_FLOW, head_m3s=HEAD_M3S),
        rivers=lines)) + '\n')
    print(f'{len(lines)} rivers written, {time.time() - started:.0f} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
