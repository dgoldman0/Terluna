"""The near-side sea's neck across western Procellarum on the climate run on the final surface: what the neck carries,
when the sea would part there, and three ways to treat it (let it part, a ship lane, a strait), each with the monthly
tide through the neck.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.near_side_neck   # a few minutes

The near-side sea's narrowest join (redrawn_surface.py: the level at which the near-side maria's seeds join, four-
connected) is a short sill in western Procellarum. Below it the sea parts into the south part (Cognitum, Humorum,
Nubium and southern Procellarum, the maria last to join) and the main part (Imbrium's). On final_drainage.py's routing
of the 0.930 run (each budget's governed routing on the research drive) this module finds:

- the parts: each cell of the near-side sea is given to the part whose water it joins below the crest (the network
  study's joining levels on the catchment graph, as the drainage's own tables), or, where it joins only above the
  crest, to the nearer part; each part's loss to the air beyond the rain at any level (network.Table on its cells) and
  the rivers it receives (the governed routing's mouths into the near-side sea, given to the part their sea cell is in);
- the neck: its crest, its throat at each level (the narrowest wetted section of the water joining the parts, a
  minimum cut through the shallow water between them, four-connected, each cell weighed by its depth times its mean
  side), and its sill (the length of the least-excavation route between the parts' water 3 m under the crest);
- the flow through it: the south part's supply less its loss, which the main part sends; the two parts' levels when
  that flow must pass the neck, with quadratic bed drag over the sill and the jet's velocity head lost beyond it,
  capped by critical flow over the throat; and how much more the near-side sea could lose before the south part falls
  half a metre and a metre behind and before the neck chokes;
- what happens past that with no work: the main part standing over the crest and spilling to the south part, which
  settles where its loss meets its rivers and the spill;
- a ship lane and a strait: the least excavation for a channel at each bed level between the parts' water deep enough
  already (eight-connected cells, in a window about the throat), the lane's volume at 300 m wide and the strait's width
  for the whole southern flow at a set drop, both with side slopes of 1 in 3, and the margin each gives;
- the monthly tide through the neck: a two-basin model in which each part's water follows geography's equilibrium tide
  at its cells' turned places plus a level of its own, the sea's volume held, and the neck passes the flow the head
  between the parts drives through its sections, with the southern flow on it.

All flows are first estimates. The bed drag is the ocean models' central coefficient (0.0025); the sill's sections are
taken as the throat's over the sill's length, which overstates its friction; the throat is found on 1.9 km cells. The
tide model leaves out the water's self-attraction (9% on the near side's seas), the parts' own seiches and the neck's
inertia, which at monthly periods is small.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.csgraph import breadth_first_order, dijkstra, maximum_flow

from climate.gcm.turned import to_turned
from climate.ocean import forcing as F
from climate.ocean import lunar
from geography import crustal_loading as CL
from geography import drainage as dr
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.drainage_network import network as NW

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'near_side_neck.json'
SCHEMA = 'terluna.research.atlas-round-near-side-neck/1'
CODE = ('research/studies/atlas_round/near_side_neck.py', 'research/studies/atlas_round/final_drainage.py',
        'research/studies/atlas_round/redrawn_drainage.py', 'research/studies/drainage_network/network.py',
        'geography/crustal_loading.py', 'climate/ocean/forcing.py')
SOUTH = ('Cognitum', 'Humorum', 'Nubium')
G = lunar.GRAVITY
DRAG = lunar.DRAG
EXIT_HEADS = 1.0                     # velocity heads lost where the neck's flow spreads into the water beyond
SILL_BELOW_CREST_M = 3.0             # the sill's length: the route between the parts' water this far under the crest
THROAT_ABOVE_CREST_M = (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 45.0, 69.0)
LANE_WIDTH_M = 300.0
LANE_DEPTHS_M = (10.0, 15.0, 20.0)
SIDE_SLOPE = 3.0                     # horizontal over vertical, each side
STRAIT_DROPS_M = (0.5, 1.0)
STRAIT_BELOW_CREST_M = (3.0, 11.0)   # the shared level a strait holds under the crest
STRAIT_WATER_M = 18.0                # the strait's water depth at that level
WINDOW_DEG = (8.0, 10.0)             # half-height and half-width of the channel search's window
EXTRA_LOSS_M3S = (0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0)
TIDE_DAYS = 2 * 365.25
TIDE_STEP_S = 3600.0
TIDE_SPINUP_DAYS = 60.0


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------- the near side and its parts

def setting(log=print):
    """The redrawn floor, the crest, the parts' cells and the rolled near-side box the throat is found on."""
    surf = RD.load_surface(log)
    h, level0, lat, cell_km = surf['h'], surf['level0'], surf['lat'], surf['cell_km']
    rows, cols = h.shape
    flat = h.ravel()
    mares = {k: CL.deepest_near(h, *v) for k, v in CL.NEARSIDE_MARIA.items()}
    join = CL.nearside_join(h, level0, mares)
    label, sinks, edges = RD.catchment_graph(h, lat, cell_km)
    n_nodes = sinks.size
    seed_m = mares['Imbrium'][0] * cols + mares['Imbrium'][1]
    seed_s = min((mares[k][0] * cols + mares[k][1] for k in SOUTH), key=lambda c: flat[c])
    joined_m = NW.joining_level(flat, label, n_nodes, edges, seed_m)
    joined_s = NW.joining_level(flat, label, n_nodes, edges, seed_s)
    crest = float(joined_m[seed_s])                  # the minimax pass between the parts on the catchment graph
    del label, sinks, edges
    area = np.repeat((cell_km * 1e3) ** 2 * np.cos(np.radians(lat)), cols)
    south_core = joined_s < crest - 1e-3
    main_core = joined_m < crest - 1e-3
    near = joined_m < level0 + 200.0                 # the near-side sea's water up to 200 m over the common level
    shift = cols // 2                                # the near side contiguous, away from the 0/360 seam
    cores = np.roll((south_core | main_core).reshape(rows, cols), shift, axis=1)
    _, (ii, jj) = ndimage.distance_transform_edt(~cores, return_indices=True)
    nearer_south = np.roll(np.roll(south_core.reshape(rows, cols), shift, axis=1)[ii, jj], -shift, axis=1).ravel()
    del ii, jj, cores
    south = near & (south_core | (~main_core & nearer_south))
    main = near & ~south
    log(f'crest {crest - level0:+.2f} m on the catchment graph, {join["join_level_m"] - level0:+.2f} m four-connected')
    return dict(surf=surf, h=h, flat=flat, level0=level0, lat=lat, cell_km=cell_km, rows=rows, cols=cols,
                mares=mares, join4=float(join['join_level_m']), crest=crest, seeds=dict(main=seed_m, south=seed_s),
                joined=dict(main=joined_m, south=joined_s), part=dict(main=main, south=south), area=area)


def tables(st, loss):
    """Each part's network.Table (area, loss, volume against level) on a loss field, and both parts held at one level
    ('joined', as a strait would hold them)."""
    out = {}
    both = np.full(st['flat'].size, np.inf)
    for k in ('main', 'south'):
        joined = np.where(st['part'][k], st['joined'][k], np.inf).astype(np.float64)
        out[k] = NW.Table(joined, st['area'], loss, st['flat'], st['level0'] + 1000.0)
        both = np.where(st['part'][k], joined, both)
    out['joined'] = NW.Table(both, st['area'], loss, st['flat'], st['level0'] + 1000.0)
    return out


def supplies(st, routing_path):
    """The rivers each part receives in a governed routing: its mouths into the near-side sea, by the part of the sea
    cell each one enters."""
    with np.load(routing_path, allow_pickle=False) as d:
        meta = json.loads(str(d['metadata']))
        receiver, discharge, sea = d['receiver'], d['discharge_m3s'].astype(np.float64), d['sea']
    ns = meta['seas'].index('near_side')
    mouths = np.flatnonzero((sea < 0) & (sea[receiver] == ns))
    into = receiver[mouths]
    q = discharge[mouths]
    south = st['part']['south'][into]
    return dict(main=float(q[~south].sum()), south=float(q[south].sum()), mouths=int(mouths.size),
                total=float(q.sum()))


# ---------------------------------------------------------------- the throat and the channels, on the rolled box

def box(st):
    """The near side rolled so it is contiguous (crustal_loading.nearside_join's box), with the parts under the crest."""
    cols, ppd = st['cols'], dr.PPD
    shift = cols // 2
    r0, r1 = int((90 - 75) * ppd), int((90 + 50) * ppd)
    c0, c1 = int(80 * ppd), int(280 * ppd)
    roll = lambda a: np.roll(a.reshape(st['rows'], cols), shift, axis=1)[r0:r1, c0:c1]
    sub = roll(st['h']).astype(np.float64)
    lab, _ = ndimage.label(sub < st['join4'] - 0.05)
    pts = {k: (i - r0, (j + shift) % cols - c0) for k, (i, j) in st['mares'].items()}
    a_id = {lab[pts[k]] for k in SOUTH} - {0}
    b_id = lab[pts['Imbrium']]
    if len(a_id) != 1 or b_id == 0 or b_id in a_id:
        raise SystemExit(f'the parts under the crest are not as expected: south {a_id}, main {b_id}')
    lat = st['lat'][r0:r1]
    lon = (np.arange(c0, c1) + 0.5) / ppd - 180.0
    return dict(sub=sub, A=lab == a_id.pop(), B=lab == b_id, lat=lat, lon=lon, pts=pts,
                w_row=st['cell_km'] * 1e3 * np.sqrt(np.cos(np.radians(lat))))


def throat(bx, level):
    """The smallest wetted section of the water joining the parts at a level: a minimum cut through the shallow water
    between them, each cell weighed by its depth times its mean side. Returns the section, its top width, its deepest
    water, its centre and the cut's cells, or None where the parts do not join at that level."""
    sub, A, B = bx['sub'], bx['A'], bx['B']
    lab, _ = ndimage.label(sub < level)
    pa, pb = bx['pts'][SOUTH[0]], bx['pts']['Imbrium']
    if lab[pa] == 0 or lab[pa] != lab[pb]:
        return None
    wet = lab == lab[pb]
    S = wet & ~A & ~B
    rows, cols = sub.shape
    idx = -np.ones(sub.shape, np.int64)
    cells = np.flatnonzero(S.ravel())
    n = cells.size
    idx.ravel()[cells] = np.arange(n)
    depth = (level - sub).ravel()[cells]
    width = np.repeat(bx['w_row'], cols)[cells]
    cap = depth * width
    unit = max(1.0, float(cap.sum()) / 4e8)
    c = np.maximum(np.round(cap / unit), 1).astype(np.int64)
    inf = int(c.sum() + 1)
    src, snk = 2 * n, 2 * n + 1
    r_, c_, v_ = [np.arange(n) * 2], [np.arange(n) * 2 + 1], [c]
    ci, cj = np.divmod(cells, cols)
    for di, dj in ((0, 1), (1, 0), (0, -1), (-1, 0)):
        ni, nj = ci + di, cj + dj
        ok = (ni >= 0) & (ni < rows) & (nj >= 0) & (nj < cols)
        u = np.flatnonzero(ok)
        nb_i, nb_j = ni[u], nj[u]
        nb = idx[nb_i, nb_j]
        s_nb = nb >= 0
        r_.append(u[s_nb] * 2 + 1); c_.append(nb[s_nb] * 2); v_.append(np.full(s_nb.sum(), inf))
        a_nb = A[nb_i, nb_j]
        r_.append(np.full(a_nb.sum(), src)); c_.append(u[a_nb] * 2); v_.append(np.full(a_nb.sum(), inf))
        b_nb = B[nb_i, nb_j]
        r_.append(u[b_nb] * 2 + 1); c_.append(np.full(b_nb.sum(), snk)); v_.append(np.full(b_nb.sum(), inf))
    g = sparse.coo_matrix((np.concatenate(v_), (np.concatenate(r_), np.concatenate(c_))), shape=(2 * n + 2, 2 * n + 2)).tocsr()
    g.sum_duplicates()
    g.data = np.minimum(g.data, inf)
    g = g.astype(np.int32)
    res = maximum_flow(g, src, snk)
    fp = sparse.csr_matrix(res.flow).astype(np.int64).maximum(0)
    resid = (g.astype(np.int64) - fp + fp.T).tocsr()
    resid.data[resid.data < 0] = 0
    resid.eliminate_zeros()
    seen = np.zeros(2 * n + 2, bool)
    seen[breadth_first_order(resid, src, directed=True, return_predecessors=False)] = True
    cut = np.flatnonzero(seen[0:2 * n:2] & ~seen[1:2 * n:2])
    ti, tj = np.divmod(cells[cut], cols)
    return dict(section_m2=float((depth[cut] * width[cut]).sum()), width_m=float(width[cut].sum()),
                deepest_m=float(depth[cut].max()), cells=(ti, tj),
                centre=[round(float(bx['lat'][int(round(ti.mean()))]), 2), round(float(bx['lon'][int(round(tj.mean()))]), 2)])


def section_at(bx, cut_cells, level):
    """A throat's wetted section and top width (m2, m) at a water level, on its own cells."""
    ti, tj = cut_cells
    d = np.maximum(level - bx['sub'][ti, tj], 0.0)
    w = bx['w_row'][ti]
    return float((d * w).sum()), float(w[d > 0].sum())


def channel(bx, bed, centre, level):
    """The least excavation (cut below `bed`, summed along the way) for a channel between the parts' water already that
    deep, eight-connected in a window about `centre`: its length, the stretch dug, the deepest and mean cut, the
    volume per metre of bed width and the side slopes' volume."""
    sub, A, B, lat, lon = bx['sub'], bx['A'], bx['B'], bx['lat'], bx['lon']
    ic = int(np.argmin(np.abs(lat - centre[0])))
    jc = int(np.argmin(np.abs(lon - centre[1])))
    wi, wj = int(WINDOW_DEG[0] * dr.PPD), int(WINDOW_DEG[1] * dr.PPD)
    i0, i1, j0, j1 = max(ic - wi, 0), min(ic + wi, sub.shape[0]), max(jc - wj, 0), min(jc + wj, sub.shape[1])
    win = sub[i0:i1, j0:j1]
    wr, wc = win.shape
    cell_m = (lon[1] - lon[0]) * math.pi / 180.0 * lunar.RADIUS
    dx = cell_m * np.cos(np.radians(lat[i0:i1]))
    dy = cell_m
    cutd = np.maximum(win - bed, 0.0).ravel()
    I, J = np.mgrid[0:wr, 0:wc]
    I, J = I.ravel(), J.ravel()
    r_, c_, v_ = [], [], []
    for di, dj in ((0, 1), (1, 0), (1, 1), (1, -1)):
        ni, nj = I + di, J + dj
        ok = (ni < wr) & (nj >= 0) & (nj < wc)
        u, v = I[ok] * wc + J[ok], ni[ok] * wc + nj[ok]
        step = np.sqrt((0.5 * (dx[I[ok]] + dx[ni[ok]]) * abs(dj)) ** 2 + (dy * di) ** 2)
        wgt = 0.5 * (cutd[u] + cutd[v]) * step + 1e-6 * step
        r_ += [u, v]; c_ += [v, u]; v_ += [wgt, wgt]
    n = wr * wc
    src = np.flatnonzero((A[i0:i1, j0:j1] & (win <= bed)).ravel())
    snk = np.flatnonzero((B[i0:i1, j0:j1] & (win <= bed)).ravel())
    if src.size == 0 or snk.size == 0:
        return None
    r_.append(np.full(src.size, n)); c_.append(src); v_.append(np.full(src.size, 1e-9))
    g = sparse.coo_matrix((np.concatenate(v_), (np.concatenate(r_), np.concatenate(c_))), shape=(n + 1, n + 1)).tocsr()
    dist, pred = dijkstra(g, directed=True, indices=n, return_predecessors=True)
    end = int(snk[np.argmin(dist[snk])])
    path = [end]
    while pred[path[-1]] != n and pred[path[-1]] >= 0:
        path.append(int(pred[path[-1]]))
    path = np.array(path[::-1])
    pi, pj = np.divmod(path, wc)
    steps = np.sqrt((np.diff(pj) * dx[pi[:-1]]) ** 2 + (np.diff(pi) * dy) ** 2)
    mid = 0.5 * (cutd[path][1:] + cutd[path][:-1])
    dug = mid > 0
    ground = 0.5 * (win.ravel()[path][1:] + win.ravel()[path][:-1])
    edge = pi.min() == 0 or pi.max() == wr - 1 or pj.min() == 0 or pj.max() == wc - 1
    return dict(route_km=round(float(steps.sum()) / 1e3, 1), dug_km=round(float(steps[dug].sum()) / 1e3, 1),
                deepest_cut_m=round(float(cutd[path].max()), 1),
                mean_cut_where_dug_m=round(float((mid[dug] * steps[dug]).sum() / max(steps[dug].sum(), 1.0)), 1),
                m3_per_m_of_bed_width=float((mid * steps).sum()),
                side_slopes_m3=float((SIDE_SLOPE * mid ** 2 * steps).sum()),
                through_land_km=round(float(steps[ground > level].sum()) / 1e3, 1),
                starts=[round(float(lat[i0 + pi[0]]), 2), round(float(lon[j0 + pj[0]]), 2)],
                ends=[round(float(lat[i0 + pi[-1]]), 2), round(float(lon[j0 + pj[-1]]), 2)],
                touches_the_window_edge=bool(edge))


# ---------------------------------------------------------------- the flow through the neck

class Neck:
    """The neck's hydraulics: the throat's section and top width at any level (interpolated between the levels its
    cut was found at), the sill's length, and channels beside it (bed level, bed width, length) that add their own
    conveyance."""

    def __init__(self, crest, levels, sections, widths, sill_m, channels=()):
        self.crest, self.sill_m = crest, sill_m
        self.lv, self.sec, self.wid = np.asarray(levels), np.asarray(sections), np.asarray(widths)
        self.channels = list(channels)

    def natural(self, level):
        if level <= self.crest:
            return 0.0, 0.0
        return float(np.interp(level, self.lv, self.sec)), float(np.interp(level, self.lv, self.wid))

    def lowest(self):
        """The lowest level the main part can stand at and still send water: the crest or the deepest channel's bed."""
        return min([self.crest] + [bed for bed, _, _ in self.channels])

    def section(self, upstream, downstream):
        """The wetted section (m2) the flow passes through at two levels: the throat's and each channel's, each taken
        at the higher of the two levels' mean and its critical water level."""
        out = 0.0
        if upstream > self.crest:
            out += self.natural(max(0.5 * (upstream + downstream), self.crest + 2.0 / 3.0 * (upstream - self.crest)))[0]
        for bed, width, _ in self.channels:
            if upstream > bed:
                y = max(0.5 * (upstream + downstream), bed + 2.0 / 3.0 * (upstream - bed)) - bed
                out += width * y + SIDE_SLOPE * y * y
        return out

    def conveyance(self, upstream, downstream):
        """The flow (m3/s) a metre of head drives: Q = K * sqrt(head), with the bed's drag over the sill (and each
        channel's length) and the jet's velocity head lost beyond it. Each section is taken at the higher of the two
        levels' mean and the critical water level over its bed (two thirds of the way up from its crest), where the
        water at the throat stands when the far side falls away."""
        out = 0.0
        if upstream > self.crest:
            y = max(0.5 * (upstream + downstream), self.crest + 2.0 / 3.0 * (upstream - self.crest))
            a, b = self.natural(y)
            if a > 0.0 and b > 0.0:
                out += a / math.sqrt(DRAG * self.sill_m / (G * a / b) + EXIT_HEADS / (2 * G))
        for bed, width, length in self.channels:
            if upstream <= bed:
                continue
            y = max(0.5 * (upstream + downstream), bed + 2.0 / 3.0 * (upstream - bed)) - bed
            a = width * y + SIDE_SLOPE * y * y
            r = a / (width + 2.0 * y * math.sqrt(1.0 + SIDE_SLOPE ** 2))
            out += a / math.sqrt(DRAG * length / (G * r) + EXIT_HEADS / (2 * G))
        return out

    def critical(self, upstream):
        """The most the throat and channels pass from water standing at `upstream` (critical flow at each)."""
        out = 0.0
        if upstream > self.crest:
            ys = np.linspace(self.crest, upstream, 200)[1:-1]
            out += float(np.max(np.interp(ys, self.lv, self.sec) * np.sqrt(2 * G * (upstream - ys))))
        for bed, width, length in self.channels:
            e = upstream - bed
            if e > 0.0:
                y = np.linspace(0.0, e, 200)[1:-1]
                a = width * y + SIDE_SLOPE * y * y
                out += float(np.max(a * np.sqrt(2 * G * (e - y))))
        return out

    def flow(self, main_level, south_level):
        """The flow from the main part to the south part at two levels (negative the other way): the smaller of what
        friction lets the head drive and what critical flow passes."""
        sign = 1.0 if main_level >= south_level else -1.0
        hi, lo = max(main_level, south_level), min(main_level, south_level)
        q = self.conveyance(hi, lo) * math.sqrt(hi - lo)
        return sign * min(q, self.critical(hi))


def level_for(table, supply, below):
    """network.Table.level_for by a search on the running greatest loss (the first level whose loss reaches the
    supply), kept on the table: the same level, found without scanning the cells each time."""
    if not hasattr(table, 'runmax'):
        table.runmax = np.maximum.accumulate(table.loss)
    n = int(np.searchsorted(table.c, below))
    i = int(np.searchsorted(table.runmax[:n], supply, side='left'))
    return float(table.c[i]) if i < n else None


def balance(tab, supply, neck, extra=0.0, cap=None):
    """The two parts' levels and the neck's flow at steady state, the main part's supply less `extra`: the flow Q at
    which the main part's level (its loss meeting its supply less Q) and the south part's (its loss meeting its rivers
    and Q) drive Q through the neck. None where the main part cannot stand over the crest."""
    sm, ss = supply['main'] - extra, supply['south']
    top = cap if cap is not None else neck.crest + 1000.0
    lev = lambda t, q: level_for(t, q, top)                      # None: the water reaches the cap
    q_hi = sm - tab['main'].at(neck.lowest())[1]               # the most the main part can spare and still send it
    if q_hi <= 0.0:
        return None

    def levels(q):
        zm = lev(tab['main'], sm - q)
        zs = lev(tab['south'], ss + q)
        return (top if zm is None else zm), (top if zs is None else zs)

    lo, hi = 0.0, q_hi
    for _ in range(80):
        q = 0.5 * (lo + hi)
        zm, zs = levels(q)
        if neck.flow(zm, zs) > q:
            lo = q
        else:
            hi = q
    q = 0.5 * (lo + hi)
    zm, zs = levels(q)
    return dict(flow_m3s=q, main_level=zm, south_level=zs, head_m=zm - zs,
                choked=bool(q >= 0.999 * neck.critical(zm)))


# ---------------------------------------------------------------- the monthly tide through the neck

def part_monomials(st, levels, frame_rotation):
    """Each part's area-mean of geography's forcing monomials at its water's turned places, and its area, at each
    part's level."""
    out = {}
    for k in ('main', 'south'):
        cells = np.flatnonzero(st['part'][k] & (st['joined'][k] < levels[k]))
        lat = 90.0 - (cells // st['cols'] + 0.5) / dr.PPD
        lon = ((cells % st['cols'] + 0.5) / dr.PPD + 180.0) % 360.0 - 180.0
        t_lat, t_lon = to_turned(lat, lon, frame_rotation)
        w = st['area'][cells]
        m = np.zeros(20)
        for i in range(0, cells.size, 200000):
            m += (F.monomials(t_lat[i:i + 200000], t_lon[i:i + 200000]) * w[i:i + 200000, None]).sum(axis=0)
        out[k] = dict(mean=m / w.sum(), area=float(w.sum()))
    return out


def tide(neck, steady, means, neck_monomials, days=TIDE_DAYS, step=TIDE_STEP_S, spinup=TIDE_SPINUP_DAYS):
    """The two-basin tide. Each part's water stands at the equilibrium tide's pattern plus a level of its own, the
    sea's volume held; the south part's own level moves by what the neck passes beyond the steady southern flow, and
    the neck passes what the head between the parts' own levels drives (Neck.flow at the levels on the neck's two
    sides). Returns, after a spin-up, the tidal exchange against the open neck's, the neck's flow and the throat's
    mean speed, and the south part's monthly range against the open neck's."""
    t = np.arange(0.0, (days + spinup) * 86400.0, step)
    coef = F.ephemeris_series(t)
    zm_bar, zs_bar = coef @ means['main']['mean'], coef @ means['south']['mean']
    z_neck = coef @ neck_monomials
    am, as_ = means['main']['area'], means['south']['area']
    sea = (am * zm_bar + as_ * zs_bar) / (am + as_)
    open_s = zs_bar - sea                            # the south part's mean tidal level with the neck wide open
    q0, zm0, zs0 = steady['flow_m3s'], steady['main_level'], steady['south_level']
    s_level = open_s[0]                              # the south part's mean tidal level, started open
    south = np.empty_like(t)
    q = np.empty_like(t)
    speed = np.empty_like(t)
    for n in range(t.size):
        c_s = s_level - zs_bar[n]
        c_m = -(as_ * s_level) / am - zm_bar[n]
        up, down = zm0 + z_neck[n] + c_m, zs0 + z_neck[n] + c_s
        q[n] = neck.flow(up, down)
        speed[n] = abs(q[n]) / max(neck.section(max(up, down), min(up, down)), 1.0)
        south[n] = s_level
        s_level += (q[n] - q0) * step / as_
    keep = t >= spinup * 86400.0
    month = 27.554550 * 86400.0
    ranges, open_ranges = [], []
    for a in np.arange(spinup * 86400.0, t[-1] - month, month):
        sel = (t >= a) & (t < a + month)
        ranges.append(float(np.ptp(south[sel])))
        open_ranges.append(float(np.ptp(open_s[sel])))
    open_q = np.gradient(open_s, t) * as_
    tidal = q[keep] - q0
    return dict(
        exchange_open_neck_largest_m3s=round(float(np.abs(open_q[keep]).max())),
        exchange_largest_m3s=round(float(np.abs(tidal).max())), exchange_rms_m3s=round(float(np.sqrt(np.mean(tidal ** 2)))),
        neck_flow_m3s=dict(least=round(float(q[keep].min())), greatest=round(float(q[keep].max())), steady=round(q0)),
        throat_mean_speed_largest_m_s=round(float(speed[keep].max()), 2),
        throat_mean_speed_median_m_s=round(float(np.median(speed[keep])), 2),
        throat_section_steady_m2=round(neck.section(zm0, zs0)),
        south_monthly_range_m=round(float(np.median(ranges)), 2),
        south_open_neck_monthly_range_m=round(float(np.median(open_ranges)), 2),
        parts_equilibrium_difference_monthly_range_m=round(float(np.median(
            [np.ptp((zs_bar - zm_bar)[(t >= a) & (t < a + month)]) for a in np.arange(spinup * 86400.0, t[-1] - month, month)])), 2))


# ---------------------------------------------------------------- the record

def rel(st, z):
    return round(float(z) - st['level0'], 2)


def threshold(tab, sup, neck, test, steps=40):
    """The least extra loss from the main part's rivers at which `test(state)` holds (0 where it holds already; a state
    in which the main part cannot stand over the crest counts as holding), by doubling and bisection."""
    s0 = balance(tab, sup, neck, 0.0)
    if s0 is None or test(s0):
        return 0
    lo, hi = 0.0, 100.0
    while True:
        s1 = balance(tab, sup, neck, hi)
        if s1 is None or test(s1):
            break
        lo, hi = hi, hi * 2.0
        if hi > 2e5:
            return None
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        s1 = balance(tab, sup, neck, mid)
        if s1 is None or test(s1):
            hi = mid
        else:
            lo = mid
    return round(hi)


def added_water(st, tab, sup, neck, drop):
    """The water to add to the main part's rivers (by gravity, from the seas above) at which the south part stands
    within `drop` of the main part and the neck flows unchoked, with the parts' state then; None beyond 200,000 m3/s."""
    joined = lambda s: s is not None and s['head_m'] <= drop and not s['choked']
    s0 = balance(tab, sup, neck, 0.0)
    if joined(s0):
        return dict(added_m3s=0, **state_record(st, tab, s0))
    lo, hi = 0.0, 1000.0
    while not joined(balance(tab, sup, neck, -hi)):
        lo, hi = hi, hi * 2.0
        if hi > 2e5:
            return None
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if joined(balance(tab, sup, neck, -mid)):
            hi = mid
        else:
            lo = mid
    return dict(added_m3s=round(hi), **state_record(st, tab, balance(tab, sup, neck, -hi)))


def margins(st, tab, sup, neck):
    """The extra loss the near-side sea can take from the main part's rivers before the south part stands half a metre
    and a metre under the main part, and before the neck chokes; and the parts at each extra loss in the table."""
    base = balance(tab, sup, neck)
    out = dict(steady=None if base is None else state_record(st, tab, base))
    if base is None:
        return out
    falls = lambda drop: (lambda s: s['head_m'] - base['head_m'] >= drop)
    out['extra_loss_before_m3s'] = dict(
        south_half_a_metre_behind=threshold(tab, sup, neck, falls(0.5)),
        south_a_metre_behind=threshold(tab, sup, neck, falls(1.0)),
        neck_chokes=threshold(tab, sup, neck, lambda s: s['choked']))
    rows = []
    for extra in EXTRA_LOSS_M3S:
        s = balance(tab, sup, neck, extra)
        rows.append(dict(extra_loss_m3s=extra, **({} if s is None else state_record(st, tab, s))))
    out['with_extra_loss'] = rows
    out['added_water_to_stay_one_sea'] = {f'within_{d:g}_m': added_water(st, tab, sup, neck, d) for d in (0.5, 1.0)}
    return out


def state_record(st, tab, s):
    return dict(main_level_m=rel(st, s['main_level']), south_level_m=rel(st, s['south_level']),
                head_m=round(s['head_m'], 2), neck_flow_m3s=round(s['flow_m3s']), choked=s['choked'],
                south_area_km2=round(tab['south'].at(s['south_level'])[0] / 1e6),
                main_area_km2=round(tab['main'].at(s['main_level'])[0] / 1e6))


def strait_width(bed, length, level, need, drop):
    """The bed width of a strait (bed level, length) that carries `need` at a drop `drop` with its water at `level`."""
    lo, hi = 1.0, 20000.0
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if Neck(level + 1000.0, [0.0, 1.0], [0.0, 0.0], [0.0, 0.0], 1.0, [(bed, mid, length)]).flow(level + 0.5 * drop, level - 0.5 * drop) < need:
            lo = mid
        else:
            hi = mid
    return hi


def build(log=print):
    started = time.time()
    fd = json.loads(FD.OUT.read_text())
    st = setting(log)
    level0 = st['level0']
    frame_path, record = FD.frame_record()
    rotation = np.asarray(record['frame']['rotation_today_to_turned'], float)
    bx = box(st)
    crest = st['join4']
    found = []
    for above in THROAT_ABOVE_CREST_M:
        t_ = throat(bx, crest + above)
        if t_ is not None:
            found.append((crest + above, t_))
    log(f'throats at {len(found)} levels, {time.time() - started:.0f} s')
    centre = min(found, key=lambda x: abs(x[0] - (crest + 3.0)))[1]['centre']
    sill = channel(bx, crest - SILL_BELOW_CREST_M, centre, crest)
    neck = Neck(crest, [crest] + [lv for lv, _ in found], [0.0] + [t_['section_m2'] for _, t_ in found],
                [0.0] + [t_['width_m'] for _, t_ in found], sill['dug_km'] * 1e3)
    c_lat, c_lon = centre
    t_lat, t_lon = to_turned(np.array([c_lat]), np.array([c_lon]), rotation)
    neck_mono = F.monomials(t_lat, t_lon)[0]
    out = dict(
        crest=dict(four_connected_m=rel(st, crest), catchment_graph_m=rel(st, st['crest'])),
        throat=dict(centre_today=centre, centre_turned=[round(float(t_lat[0]), 2), round(float(t_lon[0]), 2)],
                    by_level=[dict(level_m=rel(st, lv), section_m2=round(t_['section_m2']), width_km=round(t_['width_m'] / 1e3, 1),
                                   deepest_m=round(t_['deepest_m'], 1), centre=t_['centre']) for lv, t_ in found]),
        sill=dict(length_km=sill['dug_km'], route=sill),
        budgets={})
    for label in FD.BUDGETS:
        case = fd['cases'][label]
        evaporation = case['evaporation']
        _, lake_net, _, _ = FD.load(evaporation)
        loss = -lake_net.ravel().astype(np.float64) * dr.MM_DAY * st['area']
        del lake_net
        tab = tables(st, loss)
        sup = supplies(st, FD.routing_path(label))
        settle = case['seas']['near_side']['would_settle']
        b = dict(supply_m3s=dict(main=round(sup['main']), south=round(sup['south']), mouths=sup['mouths'],
                                 near_side_record=case['seas']['near_side']['supply_m3s']),
                 near_side_would_settle=settle)
        at_crest = {k: round(tab[k].at(crest)[1]) for k in ('main', 'south')}
        b['loss_at_the_crest_m3s'] = at_crest
        b['south_needs_at_the_crest_m3s'] = round(at_crest['south'] - sup['south'])
        if settle['reaches_its_rim']:
            b['note'] = 'the near-side sea rises to its rim on this budget; the neck lies far under its water'
            out['budgets'][label] = b
            log(f'{label}: reaches its rim')
            continue
        m = margins(st, tab, sup, neck)
        b['natural_neck'] = m
        base = balance(tab, sup, neck)
        if base is None:
            b['note'] = 'the main part cannot stand over the crest on this budget'
            out['budgets'][label] = b
            continue
        shared = base['main_level']
        lanes = []
        for depth in LANE_DEPTHS_M:
            route = channel(bx, shared - depth, centre, shared)
            if route is None:
                continue
            volume = route['m3_per_m_of_bed_width'] * LANE_WIDTH_M + route['side_slopes_m3']
            with_lane = Neck(neck.crest, neck.lv, neck.sec, neck.wid, neck.sill_m,
                             [(shared - depth, LANE_WIDTH_M, max(route['dug_km'], 0.1) * 1e3)])
            lm = margins(st, tab, sup, with_lane)
            lanes.append(dict(depth_under_the_level_m=depth, bed_m=rel(st, shared - depth), width_m=LANE_WIDTH_M,
                              volume_million_m3=round(volume / 1e6, 1), route=route,
                              extra_loss_before_m3s=lm.get('extra_loss_before_m3s')))
        b['ship_lanes'] = lanes
        z_one = level_for(tab['joined'], sup['main'] + sup['south'], level0 + 200.0)
        b['one_level_m'] = None if z_one is None else rel(st, z_one)
        held = (crest - 3.0, crest - 11.0) if z_one is None or z_one >= crest else (z_one, z_one - 10.0, z_one - 20.0)
        straits, strait_necks = [], []
        for z in held:
            bed = z - STRAIT_WATER_M
            route = channel(bx, bed, centre, z)
            if route is None:
                continue
            need = tab['south'].at(z)[1] - sup['south']
            whole_extra = sup['main'] + sup['south'] - tab['main'].at(z)[1] - tab['south'].at(z)[1]
            length = max(route['dug_km'], 0.1) * 1e3
            for drop in STRAIT_DROPS_M:
                w = strait_width(bed, length, z, need, drop)
                volume = route['m3_per_m_of_bed_width'] * w + route['side_slopes_m3']
                straits.append(dict(shared_level_m=rel(st, z), under_the_crest_m=round(crest - z, 1), bed_m=rel(st, bed),
                                    water_m=STRAIT_WATER_M, drop_m=drop, south_need_m3s=round(need),
                                    bed_width_m=round(w), volume_km3=round(volume / 1e9, 3),
                                    mean_speed_m_s=round(need / (w * STRAIT_WATER_M + SIDE_SLOPE * STRAIT_WATER_M ** 2), 2),
                                    extra_loss_it_holds_m3s=max(round(whole_extra), 0), route=route))
                strait_necks.append((straits[-1], Neck(neck.crest, neck.lv, neck.sec, neck.wid, neck.sill_m,
                                                       [(bed, w, length)])))
        b['straits'] = straits
        levels = {'main': base['main_level'], 'south': base['south_level']}
        b['tide'] = dict(natural_neck=tide(neck, base, part_monomials(st, levels, rotation), neck_mono))
        if lanes:
            deepest = lanes[-1]
            with_lane = Neck(neck.crest, neck.lv, neck.sec, neck.wid, neck.sill_m,
                             [(deepest['bed_m'] + level0, LANE_WIDTH_M, max(deepest['route']['dug_km'], 0.1) * 1e3)])
            base_l = balance(tab, sup, with_lane)
            if base_l is not None:
                lv = {'main': base_l['main_level'], 'south': base_l['south_level']}
                b['tide'][f'with_the_{deepest["depth_under_the_level_m"]:g}_m_lane'] = tide(
                    with_lane, base_l, part_monomials(st, lv, rotation), neck_mono)
        for spec, with_strait in strait_necks:
            if spec['drop_m'] != STRAIT_DROPS_M[0] or spec['extra_loss_it_holds_m3s'] > 0:
                continue                                  # the strait that holds the one level with no extra loss
            base_s = balance(tab, sup, with_strait)
            if base_s is not None:
                lv = {'main': base_s['main_level'], 'south': base_s['south_level']}
                b['tide']['with_the_strait_at_the_one_level'] = dict(
                    steady=state_record(st, tab, base_s),
                    **tide(with_strait, base_s, part_monomials(st, lv, rotation), neck_mono))
        out['budgets'][label] = b
        log(f'{label}: {json.dumps(b)[:1500]}, {time.time() - started:.0f} s')
    rel_ = lambda p: str(Path(p).relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel_(p): sha256(p) for p in (FD.OUT, ROOT / RD.SURFACE, frame_path)},
        ignored_inputs={rel_(FD.routing_path(k)): sha256(FD.routing_path(k)) for k in FD.BUDGETS},
        evidence=('A calculation on products: the redrawn floor\'s heights between 1.9 km cells, the drainage of the '
                  'climate run on the final surface (final_drainage.py) and geography\'s equilibrium tide forcing; '
                  'the neck\'s hydraulics are first estimates with the ocean models\' central bed drag and the throat\'s '
                  'section taken over the whole sill. Not validated against any observation or ocean model run.'),
        reading_rule=('Levels in m above the redrawn common level; flows in m3/s. main and south are the near-side sea\'s '
                      'parts beyond the neck. natural_neck.steady is the two parts at steady state with the neck as it '
                      'is; extra_loss_before is how much more the near-side sea could lose from the main part\'s rivers '
                      'before the south part stands half a metre or a metre further under the main part, or the neck '
                      'chokes (critical flow: the parts are separate, the main part spilling to the south one). Lanes '
                      'are 300 m wide at their bed, straits as wide as the south part\'s need at the shared level and '
                      'drop requires; volumes include side slopes of 1 in 3. tide is the two-basin model over two '
                      'years after a 60-day spin-up.'),
        **out, seconds=round(time.time() - started))


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
