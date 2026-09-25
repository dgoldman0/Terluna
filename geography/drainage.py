"""Rivers and lakes on the Open Moon at the scenario water share.

    python -m geography.drainage      # results/drainage.json and products/drainage_<pct>pct_16ppd.npz

Water is routed over LOLA at 16 pixels per degree (1.9 km) above the GRAIL geoid, with the atlas
product's sea level. Every cell below that level is sea or a lake at the common level, as the atlas
fills them; these are the outlets. Above it:

1. Each cell drains to the steepest of its eight neighbours.
2. Cells that drain nowhere are the floors of closed depressions. A priority flood over the graph of
   catchments and their lowest shared passes gives each depression its spill level, the lowest pass on
   its way to an outlet, and the pass it spills through.
3. Runoff is the climate run's precipitation less evaporation over land, where positive. It accumulates
   downstream. A depression fills to its spill level and overflows when its inflow and the rain on its
   lake exceed the lake's open-water evaporation (the run's evaporation over sea at that latitude);
   otherwise it holds the smaller lake whose evaporation balances its inflow and passes nothing on.

The balance is a steady annual mean at the climate model's resolution (T21, about 170 km). Groundwater,
infiltration losses and travel time come later.
"""
from __future__ import annotations
import hashlib
import heapq
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA = 'terluna.geography.drainage/1'
PPD = 16
CLIMATOLOGY = ROOT / 'climate' / 'gcm' / 'products' / 'climatology_A.npz'
CLIMATOLOGY_SCHEMA = 'terluna.climate.gcm-climatology/1'
CHANNEL_M3S = 1.0                         # discharge kept in the product
MM_DAY = 1.0 / 1000.0 / 86400.0           # mm/day to m/s
DIRECTIONS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
EVIDENCE = ('Steady annual routing of the climate run\'s runoff (precipitation less evaporation over land, T21) '
            'over LOLA 16 pixel/degree heights above the GRAIL geoid, with fill-and-spill lakes decided by each '
            'depression\'s water balance against the run\'s open-water evaporation. A first estimate of where '
            'rivers run and which closed depressions hold lakes; channels are grid-resolved (1.9 km), and '
            'groundwater, infiltration and travel time come later.')
READING_RULE = ('Cells index rows from 90 N and columns from 0 E, 16 per degree, as LOLA: cell = row * 5760 + '
                'column. channel_discharge_m3s[k] is the mean discharge through cell channel_cell[k]. lake_cell '
                'lists land cells under rain-fed lakes above sea level, with lake_level_m their surface above the '
                'geoid and lake_spills whether that lake overflows.')


def shifted(a, dr, dc, fill):
    """out[r, c] = a[r + dr, (c + dc) % cols]; rows beyond the poles get fill."""
    out = np.roll(a, -dc, axis=1) if dc else a
    if dr == -1:
        out = np.vstack([np.full((1, a.shape[1]), fill, a.dtype), out[:-1]])
    elif dr == 1:
        out = np.vstack([out[1:], np.full((1, a.shape[1]), fill, a.dtype)])
    return out


def d8_receivers(h, outlet, lat_deg, cell_km):
    """Index of each cell's steepest-descent neighbour (itself when none is lower, and for outlets)."""
    rows, cols = h.shape
    idx = np.arange(rows * cols, dtype=np.int32).reshape(rows, cols)
    dx = cell_km * np.maximum(np.cos(np.radians(lat_deg)), 1e-3)[:, None]
    best = np.zeros(h.shape, np.float32)
    receiver = idx.copy()
    for dr, dc in DIRECTIONS:
        slope = (h - shifted(h, dr, dc, np.inf)) / np.sqrt((dr * cell_km) ** 2 + (dc * dx) ** 2)
        better = slope > best
        best = np.where(better, slope, best)
        receiver = np.where(better, shifted(idx, dr, dc, -1), receiver)
    receiver[outlet] = idx[outlet]
    return receiver.ravel()


def terminal(receiver):
    """The cell each cell's path ends at, by pointer jumping."""
    root = receiver.copy()
    while True:
        nxt = root[root]
        if np.array_equal(nxt, root):
            return root
        root = nxt


def passes(label, h_pass, idx):
    """Lowest pass between every pair of touching catchments: (u, v, height, cell in u, cell in v), u < v."""
    parts = []
    for dr, dc in ((0, 1), (1, 0), (1, 1), (1, -1)):
        lb = shifted(label, dr, dc, -1)
        m = (lb >= 0) & (lb != label)
        a, b = label[m], lb[m]
        w = np.maximum(h_pass[m], shifted(h_pass, dr, dc, np.inf)[m])
        ca, cb = idx[m], shifted(idx, dr, dc, -1)[m]
        swap = a > b
        u, v = np.where(swap, b, a), np.where(swap, a, b)
        cu, cv = np.where(swap, cb, ca), np.where(swap, ca, cb)
        parts.append(_lowest(u, v, w, cu, cv))
    return _lowest(*(np.concatenate(x) for x in zip(*parts)))


def _lowest(u, v, w, cu, cv):
    key = u.astype(np.int64) * (int(v.max(initial=0)) + 1) + v
    order = np.lexsort((w, key))
    key = key[order]
    first = np.ones(key.size, bool)
    first[1:] = key[1:] != key[:-1]
    keep = order[first]
    return u[keep], v[keep], w[keep], cu[keep], cv[keep]


def flood(n_nodes, edges, sea_level):
    """Spill level of every catchment (the minimax pass height on its way to node 0, the outlets) and the pass
    it spills through: the cell on its own side and the cell on its parent's side."""
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
    heap = [(float(sea_level), 0, -1, -1)]
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


def route(h, outlet, lat_deg, cell_km, runoff_mm_day, lake_net_mm_day):
    """Route runoff over a latitude-longitude grid (longitude periodic) with fill-and-spill lakes.

    h: heights (m); outlet: cells that take water in (sea); runoff_mm_day: generated on land;
    lake_net_mm_day: precipitation less open-water evaporation on a lake surface. Returns a dict of flat
    arrays and per-depression tables."""
    rows, cols = h.shape
    n = rows * cols
    idx2 = np.arange(n, dtype=np.int32).reshape(rows, cols)
    receiver = d8_receivers(h.astype(np.float32), outlet, lat_deg, cell_km)
    flat_out = outlet.ravel()
    idx = idx2.ravel()
    sink = (receiver == idx) & ~flat_out
    sinks = np.flatnonzero(sink).astype(np.int32)
    node_of = np.full(n, -1, np.int32)
    node_of[flat_out] = 0
    node_of[sinks] = np.arange(1, sinks.size + 1, dtype=np.int32)
    label = node_of[terminal(receiver)]
    n_nodes = sinks.size + 1
    h_pass = np.where(outlet, -np.inf, h).astype(np.float32)
    edges = passes(label.reshape(rows, cols), h_pass, idx2)
    sea = float(np.max(h[outlet])) if outlet.any() else float(h.min())
    level, spill_from, spill_to = flood(n_nodes, edges, sea)
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
    final[sinks] = spill_to[1:]
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
    closed[1:] = (inflow[1:] + lake_net[1:] < 0.0) & (lake_net[1:] < 0.0)
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
                sinks=sinks, spill_from=spill_from, spill_to=spill_to, area=area, outlet=flat_out)


def lake_evaporation(clim, lat_deg):
    """Open-water evaporation (mm/day) at each latitude: the run's zonal mean over sea cells."""
    sea = clim['lsm'] < 0.5
    ev = np.where(sea, clim['evap_mm_day'], np.nan)
    zonal = np.nanmean(ev, axis=1)
    lat = clim['lat']
    ok = np.isfinite(zonal)
    order = np.argsort(lat[ok])
    return np.interp(lat_deg, lat[ok][order], zonal[ok][order])


def climate_on_grid(clim, name, rows, cols):
    """A climatology field on the 16 pixel/degree grid (rows from 90 N, columns from 0 E), bilinear."""
    lat = 90.0 - (np.arange(rows) + 0.5) / PPD
    lon = (np.arange(cols) + 0.5) / PPD
    f = clim[name]
    order = np.argsort(clim['lat'])
    ls, fs = clim['lat'][order], f[order]
    lon_src = np.concatenate([clim['lon'] - 360.0, clim['lon'], clim['lon'] + 360.0])
    fs = np.concatenate([fs, fs, fs], axis=1)
    i = np.clip(np.searchsorted(ls, np.clip(lat, ls[0], ls[-1])) - 1, 0, ls.size - 2)
    ty = ((np.clip(lat, ls[0], ls[-1]) - ls[i]) / (ls[i + 1] - ls[i]))[:, None].astype(np.float32)
    j = np.clip(np.searchsorted(lon_src, lon) - 1, 0, lon_src.size - 2)
    tx = ((lon - lon_src[j]) / (lon_src[j + 1] - lon_src[j]))[None, :].astype(np.float32)
    top = fs[i][:, j].astype(np.float32) * (1 - tx) + fs[i][:, j + 1].astype(np.float32) * tx
    bottom = fs[i + 1][:, j].astype(np.float32) * (1 - tx) + fs[i + 1][:, j + 1].astype(np.float32) * tx
    return top * (1 - ty) + bottom * ty


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compute():
    from geography import atlas as ga
    from geography import topography as tp
    from shared.constants import MOON_RADIUS
    atlas_path = HERE / 'results' / 'atlas.json'
    atlas = json.loads(atlas_path.read_text())
    if atlas.get('schema') != ga.SCHEMA:
        raise SystemExit(f'Expected {ga.SCHEMA} in {atlas_path}; run python -m geography.atlas')
    clim = dict(np.load(CLIMATOLOGY))
    clim_meta = json.loads(str(clim.pop('metadata')))
    if clim_meta.get('schema') != CLIMATOLOGY_SCHEMA:
        raise SystemExit(f'Expected {CLIMATOLOGY_SCHEMA} in {CLIMATOLOGY}')
    level = float(atlas['sea_level_m'])
    h, _, grid = tp.height_above_geoid(PPD, atlas['grid']['geoid_degree'])
    h = h.astype(np.float32)
    rows, cols = h.shape
    lat = 90.0 - (np.arange(rows) + 0.5) / PPD
    cell_km = MOON_RADIUS / 1000.0 * np.radians(1.0 / PPD)
    outlet = h < level
    precip = climate_on_grid(clim, 'pr_mm_day', rows, cols)
    runoff = np.maximum(precip - climate_on_grid(clim, 'evap_mm_day', rows, cols), 0.0)
    lake_net = precip - lake_evaporation(clim, lat)[:, None].astype(np.float32)
    del precip
    r = route(h, outlet, lat, cell_km, runoff, lake_net)
    return atlas, atlas_path, clim_meta, level, h, lat, cell_km, runoff, r


def summarise(atlas, level, h, lat, r, top=12):
    """Totals, the largest rivers by discharge at the mouth and the largest lakes as connected water bodies."""
    from scipy import ndimage
    rows, cols = h.shape
    area = r['area']
    total = float(area.sum())
    q = r['discharge']
    lake = r['lake']
    place = lambda row, col: (round(float(90.0 - (row + 0.5) / PPD), 2),
                              round(float(((col + 0.5) / PPD + 180.0) % 360.0 - 180.0), 2))
    mouths = np.flatnonzero(r['outlet'] & (q > 0))
    order = mouths[np.argsort(q[mouths])[::-1]]
    rivers = [dict(mouth=place(c // cols, c % cols), discharge_m3s=round(float(q[c]), 1),
                   basin_km2=round(float(r['drained_m2'][c]) / 1e6)) for c in order[:top]]
    body, n_bodies = ndimage.label(lake.reshape(rows, cols), structure=np.ones((3, 3), bool))
    body = body.ravel()
    wet = np.flatnonzero(lake)
    depth = r['lake_level'][wet] - h.ravel()[wet]
    b_area = np.bincount(body[wet], weights=area[wet], minlength=n_bodies + 1)
    b_volume = np.bincount(body[wet], weights=area[wet] * depth, minlength=n_bodies + 1)
    b_depth = np.zeros(n_bodies + 1)
    np.maximum.at(b_depth, body[wet], depth)
    b_level = np.full(n_bodies + 1, -np.inf)
    np.maximum.at(b_level, body[wet], r['lake_level'][wet])
    b_row = np.bincount(body[wet], weights=area[wet] * (wet // cols), minlength=n_bodies + 1)
    b_closed = np.bincount(body[wet], weights=r['closed'][r['label'][wet]].astype(float), minlength=n_bodies + 1)
    lakes = []
    for k in np.argsort(b_area)[::-1][:top]:
        cells = wet[body[wet] == k]
        col = np.angle(np.exp(1j * 2 * np.pi * (cells % cols) / cols).mean()) % (2 * np.pi) / (2 * np.pi) * cols
        lakes.append(dict(centre=place(b_row[k] / b_area[k], col), area_km2=round(float(b_area[k]) / 1e6),
                          volume_km3=round(float(b_volume[k]) / 1e9), max_depth_m=round(float(b_depth[k])),
                          level_above_sea_m=round(float(b_level[k] - level)), spills=bool(b_closed[k] == 0)))
    closed_bodies = int(np.count_nonzero(b_closed[1:] > 0))
    return dict(
        schema=SCHEMA, share=atlas['share'], sea_level_m=level, grid=dict(pixels_per_degree=PPD, rows=rows, cols=cols),
        depressions=int(r['level'].size - 1),
        lakes=dict(bodies=int(n_bodies), closed_bodies=closed_bodies, area_share=round(float(area[lake].sum()) / total, 4),
                   volume_km3=round(float(b_volume.sum()) / 1e9), median_depth_m=round(float(np.median(depth))),
                   largest=lakes),
        rivers=dict(discharge_to_sea_m3s=round(float(q[r['outlet']].sum())), channels_over_100_m3s_km=round(
            float((q[~r['outlet'] & ~lake] >= 100.0).sum()) * float(np.sqrt(area.mean())) / 1000.0), largest=rivers),
        evidence=EVIDENCE, reading_rule=READING_RULE)


def main(argv=None) -> int:
    atlas, atlas_path, clim_meta, level, h, lat, cell_km, runoff, r = compute()
    summary = summarise(atlas, level, h, lat, r)
    producer = dict(domain='geography', files={f'geography/{p.name}': digest(p)[:16] for p in
                                              (Path(__file__), HERE / 'topography.py', HERE / 'atlas.py')})
    sources = {str(atlas_path.relative_to(ROOT)): digest(atlas_path), str(CLIMATOLOGY.relative_to(ROOT)): digest(CLIMATOLOGY)}
    summary.update(producer=producer, sources=sources, climatology=dict(run=clim_meta['run'], years=clim_meta['years']))
    q = r['discharge']
    keep = np.flatnonzero((q >= CHANNEL_M3S) & ~r['outlet'] & ~r['lake']).astype(np.int32)
    lake_cells = np.flatnonzero(r['lake']).astype(np.int32)
    out = HERE / 'products' / f"drainage_{round(100 * atlas['share'])}pct_{PPD}ppd.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    meta = {k: summary[k] for k in ('schema', 'share', 'sea_level_m', 'grid', 'evidence', 'reading_rule', 'producer', 'sources')}
    np.savez_compressed(out, metadata=json.dumps(meta), channel_cell=keep, channel_discharge_m3s=q[keep].astype(np.float32),
                        lake_cell=lake_cells, lake_level_m=r['lake_level'][lake_cells].astype(np.float32),
                        lake_spills=~r['closed'][r['label'][lake_cells]])
    (HERE / 'results' / 'drainage.json').write_text(json.dumps(summary, indent=2) + '\n')
    lk = summary['lakes']
    print(f"{summary['depressions']:,} closed depressions; {lk['bodies']:,} rain-fed lakes ({lk['closed_bodies']:,} without "
          f"an outlet) over {lk['area_share']:.2%} of the surface, {lk['volume_km3']:,} km3; "
          f"{summary['rivers']['discharge_to_sea_m3s']:,} m3/s reaches the sea -> {out.name}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
