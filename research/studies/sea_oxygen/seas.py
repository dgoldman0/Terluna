"""The held fresh seas' basins: the common-level waters on LOLA, their bodies, depths and latitudes, how their deep
basins connect, and the basins the sea oxygen study runs.

    python -m research.studies.sea_oxygen.seas        # about three minutes; reads LOLA, the GRAIL geoid and the atlas

Geography's drainage takes every cell of LOLA's 16 pixel/degree grid whose height above the GRAIL geoid lies below
the atlas's common level (-1,654.2 m at the 28% share) as sea, the outlet of its rivers; drainage.py computes those
heights with topography.height_above_geoid, and this does the same. The seas' bodies are those cells' connected
regions, eight neighbours, joined across the 0/360 seam and the poles by geography's own labeller (basins.seas, as
the atlas labels its bodies). For each body this gives its area, volume, mean depth, deepest point and latitude span.

How the seas connect: a large sea is several deep basins joined over sills. The sill between two basins is the depth
of the deepest water path between their deepest points (the largest depth d at which both lie in one connected body
of water deeper than d); below it each basin's deep water is its own. The four seas over 100,000 km2 are divided
into their basins at PARTITION_M, where the Near-side Sea divides into its four main parts, each basin taking the
sea's shallower water nearest it through the sea (its region grows cell by cell through the shallows). The other
bodies are their own basins.

The basins the study runs (RULE): every basin of at least RUN_MIN_KM2, which takes in the largest seas and the
deepest point; the basin nearest each pole among those of at least POLAR_MIN_KM2; one crater sea, the body of
10,000-100,000 km2 whose deepest point lies nearest the area-weighted median of all bodies under 100,000 km2, among
those within 10 degrees of their median latitude; and one shallow sea, the largest region of a large sea whose
deepest point lies above SHALLOW_SEA_DEEPEST_M. Each run basin is listed with its hypsography, its shelves, the T21
cells that hold its water, which extract.py reads, and the seas' area it stands for: its own, each other region of
its sea nearest it in mean depth and, for the crater sea, every body under 100,000 km2 not run. Writes
results/seas.json.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
DRAINAGE = ROOT / 'geography' / 'results' / 'drainage.json'
OUT = HERE / 'results' / 'seas.json'
SCHEMA = 'terluna.research.sea-oxygen-seas/1'
PPD = 16
T21_ROWS, T21_COLS = 32, 64
LARGE_KM2 = 100000.0          # the seas divided into basins
PARTITION_M = 250.0           # where the Near-side Sea divides into its four main parts
RUN_MIN_KM2 = 100000.0        # every basin this large is run
POLAR_MIN_KM2 = 10000.0       # and the basin nearest each pole among those this large
CRATER_SEA_KM2 = (10000.0, 100000.0)
CRATER_SEA_LAT_WINDOW = 10.0
SHALLOW_SEA_DEEPEST_M = 500.0  # the shallow seas: regions of the large seas whose deepest point lies above this
SPLIT_DEPTHS_M = (100.0, 250.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0)
SPLIT_MIN_KM2 = 5000.0        # parts listed at each depth
SILL_TOLERANCE_M = 2.0        # bisection's precision
DEEPER_SPLITS_M = (500.0, 1000.0, 1500.0, 2000.0)
SUB_BASIN_MIN_KM2 = 50000.0   # deeper parts whose sills are found within the run basins
SHELF_DEPTHS_M = (2.0, 5.0, 10.0, 30.0, 100.0, 200.0)
LAT_BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90))
SIZE_CLASSES_KM2 = ((0, 10), (10, 100), (100, 1000), (1000, 1e4), (1e4, 1e5), (1e5, 1e8))
WATER_CODES = ('ME', 'OC', 'SI', 'LC', 'PA')


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def weighted_quantile(x, w, q):
    i = np.argsort(x, kind='stable')
    c = np.cumsum(w[i])
    return float(np.interp(q * c[-1], c - 0.5 * w[i], x[i]))


# --- the grid -------------------------------------------------------------------------------------------------------

def heights():
    """Heights above the geoid on LOLA's 16 pixel/degree grid (m), as geography's drainage computes them, the common
    level, the cells' areas (m2), the grid and the atlas."""
    from geography import topography as tp
    atlas = json.loads(ATLAS.read_text())
    h, _, grid = tp.height_above_geoid(PPD, atlas['grid']['geoid_degree'])
    return h.astype(np.float32), float(atlas['sea_level_m']), grid.cell_area_m2, grid, atlas


def label8(mask):
    """Eight-neighbour regions of a mask spanning the whole width, joined across the 0/360 seam as geography's
    basins.seas joins them (a cell in column 0 with the cell in the last column of the same row)."""
    from scipy import ndimage
    lab, count = ndimage.label(mask, structure=np.ones((3, 3), bool))
    uf = np.arange(count + 1)

    def root(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x
    for x, y in zip(lab[:, 0], lab[:, -1]):
        if x and y:
            rx, ry = root(x), root(y)
            if rx != ry:
                uf[ry] = rx
    roots = np.array([root(x) for x in range(count + 1)])
    return roots[lab]


def t21_lat():
    """T21's Gaussian latitudes, north first, as PlaSim's grid and the boundary file hold them."""
    x, _ = np.polynomial.legendre.leggauss(T21_ROWS)
    return np.degrees(np.arcsin(x))[::-1]


def t21_cells(rows, cols, grid):
    """The T21 cell (row, column) nearest each 16 pixel/degree cell's centre."""
    lat = t21_lat()
    mid = 0.5 * (lat[:-1] + lat[1:])                     # descending
    r = np.searchsorted(-mid, -grid.lat_deg[rows])
    c = np.rint(grid.lon_deg[cols] / (360.0 / T21_COLS)).astype(int) % T21_COLS
    return r, c


# --- bodies ---------------------------------------------------------------------------------------------------------

def body_table(lab, depth, area, grid):
    """Every body's area, volume, deepest point and its cell, latitude span and centroid."""
    n = int(lab.max()) + 1
    flat = lab.ravel()
    wet = np.flatnonzero(flat)
    b = flat[wet]
    a = area.ravel()[wet]
    d = depth.ravel()[wet].astype(float)
    rows, cols = wet // lab.shape[1], wet % lab.shape[1]
    lat, lon = grid.lat_deg[rows], grid.lon_deg[cols]
    out = dict(area=np.bincount(b, a, n), volume=np.bincount(b, a * d, n))
    deepest = np.zeros(n)
    np.maximum.at(deepest, b, d)
    out['deepest'] = deepest
    lo, hi = np.full(n, 90.0), np.full(n, -90.0)
    np.minimum.at(lo, b, lat)
    np.maximum.at(hi, b, lat)
    out['lat_min'], out['lat_max'] = lo, hi
    la, lg = np.radians(lat), np.radians(lon)
    x = np.bincount(b, a * np.cos(la) * np.cos(lg), n)
    y = np.bincount(b, a * np.cos(la) * np.sin(lg), n)
    z = np.bincount(b, a * np.sin(la), n)
    out['lat'] = np.degrees(np.arctan2(z, np.hypot(x, y)))
    out['lon'] = np.degrees(np.arctan2(y, x)) % 360.0
    out['abs_lat'] = np.bincount(b, a * np.abs(lat), n) / np.maximum(out['area'], 1.0)
    order = np.lexsort((-d, b))
    first = np.r_[0, np.flatnonzero(np.diff(b[order])) + 1]
    seed = np.zeros(n, np.int64)
    seed[b[order][first]] = wet[order][first]
    out['seed'] = seed
    return out


def numbered_by_area(lab, area):
    """The labels renumbered 1, 2, ... from the largest body."""
    n = int(lab.max()) + 1
    sizes = np.bincount(lab.ravel(), area.ravel(), n)
    sizes[0] = -1.0
    ids = np.flatnonzero(sizes > 0)
    remap = np.zeros(n, np.int32)
    remap[ids[np.argsort(-sizes[ids])]] = np.arange(1, ids.size + 1, dtype=np.int32)
    return remap[lab]


def features():
    from geography import nomenclature
    every = nomenclature.features()
    return [f for f in every if f['code'] in WATER_CODES], [f for f in every if f['code'] == 'AA']


def feature_cell(f, grid):
    r = int(np.clip(np.round((90.0 - f['lat']) * PPD - 0.5), 0, grid.lat_deg.size - 1))
    c = int(np.round((f['lon'] % 360.0) * PPD - 0.5)) % grid.lon_deg.size
    return r, c


def names_in(region_lab, k, water, grid):
    """The IAU water names whose centres lie in region k, largest first."""
    out = []
    for f in sorted(water, key=lambda f: -f['diameter_km']):
        r, c = feature_cell(f, grid)
        if region_lab[r, c] == k:
            out.append(f['name'])
    return out


# --- connections ----------------------------------------------------------------------------------------------------

def window(lab, k):
    """Rows holding body k (the whole width, so that the seam join applies)."""
    rows = np.flatnonzero((lab == k).any(1))
    return slice(int(rows[0]), int(rows[-1]) + 1)


def splits(lab, k, depth, area, grid, water):
    """Body k's parts deeper than each of SPLIT_DEPTHS_M: their number and the largest ones."""
    from scipy import ndimage
    w = window(lab, k)
    body = lab[w] == k
    lat_w = np.broadcast_to(grid.lat_deg[w][:, None], body.shape)
    out = []
    for d in SPLIT_DEPTHS_M:
        sub = label8(body & (depth[w] > d))
        ids = np.unique(sub[sub > 0])
        if not ids.size:
            break
        a = ndimage.sum(area[w], sub, ids) / 1e6
        parts = []
        full = np.zeros(lab.shape, sub.dtype)
        full[w] = sub
        for i in np.argsort(a)[::-1]:
            if a[i] < SPLIT_MIN_KM2:
                break
            m = sub == ids[i]
            parts.append(dict(area_km2=round(float(a[i])),
                              mean_lat=round(float(np.average(lat_w[m], weights=area[w][m])), 1),
                              deepest_m=round(float(depth[w][m].max())), names=names_in(full, ids[i], water, grid)[:4]))
        out.append(dict(depth_m=d, parts=int(ids.size), area_km2=round(float(a.sum())), largest=parts))
    return out


class Sills:
    """Sills between cells of one body, by bisection on the depth, each pair found once."""

    def __init__(self, lab, depth):
        self.lab, self.depth, self.known = lab, depth, {}

    def __call__(self, k, a, b):
        pair = (k,) + tuple(sorted((int(a), int(b))))
        if pair in self.known:
            return self.known[pair]
        w = window(self.lab, k)
        cols = self.lab.shape[1]
        ra, ca = divmod(int(a), cols)
        rb, cb = divmod(int(b), cols)
        ra, rb = ra - w.start, rb - w.start
        body = self.lab[w] == k
        dw = self.depth[w]

        def joined(d):
            s = label8(body & (dw > d))
            return s[ra, ca] > 0 and s[ra, ca] == s[rb, cb]
        lo, hi = 0.0, float(min(dw[ra, ca], dw[rb, cb]))
        out = None
        if joined(lo):
            while hi - lo > SILL_TOLERANCE_M:
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if joined(mid) else (lo, mid)
            out = round(lo, 1)
        self.known[pair] = out
        return out


def partition(lab, large, depth):
    """Basins of the large bodies: their regions deeper than PARTITION_M, each grown cell by cell (eight neighbours,
    across the seam) through the body's shallower water until the body is covered; shallows the growth cannot reach
    (joined to their body only across a pole) become regions of their own. Returns the regions' labels, 0 outside
    the large bodies, and the growth's steps."""
    from scipy import ndimage
    rows = np.flatnonzero(np.isin(lab, large).any(1))
    w = slice(int(rows[0]), int(rows[-1]) + 1)
    in_large = np.isin(lab[w], large)
    regions = label8(in_large & (depth[w] > PARTITION_M)).astype(np.int32)
    todo = in_large & (regions == 0)
    steps = 0
    while todo.any():
        grown = ndimage.maximum_filter(regions, size=3, mode=('nearest', 'wrap'))
        take = todo & (grown > 0)
        if not take.any():
            break
        regions[take] = grown[take]
        todo &= ~take
        steps += 1
    if todo.any():
        extra = label8(todo)
        regions[todo] = extra[todo] + regions.max()
    full = np.zeros(lab.shape, np.int32)
    full[w] = regions
    return full, steps


# --- a basin --------------------------------------------------------------------------------------------------------

HYPSOGRAPHY_POINTS = 41


def hypsography(depth, area, deepest):
    """The basin's area deeper than each of 41 fractions of its deepest point, as shares of its surface area: lake_oxygen's
    lakes.py's hypsography() unrounded. That one rounds the shares to five decimals, which in a basin of millions of
    km2 sets the deepest points' shares to zero and leaves the lake model layers without volume."""
    f = np.linspace(0.0, 1.0, HYPSOGRAPHY_POINTS)
    order = np.argsort(depth)
    d, a = depth[order], area[order]
    below = a.sum() - np.concatenate([[0.0], np.cumsum(a)])      # area deeper than each sorted depth
    frac = np.interp(f * deepest, np.concatenate([[0.0], d]), below) / a.sum()
    frac[0], frac[-1] = 1.0, 0.0
    return [float(x) for x in f], [float(x) for x in np.minimum.accumulate(frac)]


def describe(cells, depth, area, grid, water, craters, region_lab, k, seed_cell):
    """One basin's figures from its cells (flat indices)."""
    from research.studies.lake_oxygen import lakes as K
    cols = grid.lon_deg.size
    r, c = cells // cols, cells % cols
    d = depth.ravel()[cells].astype(float)
    a = area.ravel()[cells]
    lat, lon = grid.lat_deg[r], grid.lon_deg[c]
    tot = a.sum()
    sr, sc = divmod(int(seed_cell), cols)
    f, frac = hypsography(d, a, float(d.max()))
    tr, tc = t21_cells(r, c, grid)
    key = tr * T21_COLS + tc
    ids, inv = np.unique(key, return_inverse=True)
    cell_area = np.bincount(inv, a)
    cell_volume = np.bincount(inv, a * d)
    t21 = [dict(row=int(i // T21_COLS), column=int(i % T21_COLS), area_km2=round(float(x) / 1e6, 1),
                mean_depth_m=round(float(v / x), 1)) for i, x, v in zip(ids, cell_area, cell_volume)]
    t21.sort(key=lambda x: -x['area_km2'])
    la, lg = np.radians(lat), np.radians(lon)
    x, y, z = (a * np.cos(la) * np.cos(lg)).sum(), (a * np.cos(la) * np.sin(lg)).sum(), (a * np.sin(la)).sum()
    names = names_in(region_lab, k, water, grid)
    return dict(area_km2=round(float(tot) / 1e6, 1), volume_km3=round(float((a * d).sum()) / 1e9, 1),
                mean_depth_m=round(float((a * d).sum() / tot), 1), deepest_m=round(float(d.max()), 1),
                deepest_at=dict(lat=round(float(grid.lat_deg[sr]), 3), lon_east=round(float(grid.lon_deg[sc]), 3)),
                centroid=dict(lat=round(math.degrees(math.atan2(z, math.hypot(x, y))), 3),
                              lon_east=round(math.degrees(math.atan2(y, x)) % 360.0, 3)),
                lat_span=[round(float(lat.min()), 2), round(float(lat.max()), 2)],
                mean_abs_lat=round(float(np.average(np.abs(lat), weights=a)), 2),
                lat_quantiles_by_area={str(q): round(weighted_quantile(lat, a, q), 2) for q in (0.05, 0.5, 0.95)},
                depth_quantiles_by_area_m={str(q): round(weighted_quantile(d, a, q), 1)
                                           for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
                shelf_area_share={f'{s:g}': round(float(a[d < s].sum() / tot), 5) for s in SHELF_DEPTHS_M},
                names=names[:8],
                crater=None if names else K.crater_of(float(grid.lat_deg[sr]), float(grid.lon_deg[sc]), craters),
                hypsography=dict(depth_fraction=f, area_fraction=frac), t21_cells=t21, seed_cell=int(seed_cell))


def key_of(entry, used):
    """A short key from the basin's largest water name, its crater or its place."""
    if entry['names']:
        base = entry['names'][0].split()[-1].lower()
    elif entry['crater']:
        base = entry['crater']['name'].split()[0].lower().replace("'", '')
    else:
        lat, lon = entry['centroid']['lat'], entry['centroid']['lon_east']
        base = f"{abs(lat):.0f}{'n' if lat >= 0 else 's'}_{lon:.0f}e"
    k, i = base, 2
    while k in used:
        k, i = f'{base}_{i}', i + 1
    used.add(k)
    return k


# --- the whole ------------------------------------------------------------------------------------------------------

def summary(lab, depth, area, grid, t, sea_level):
    wet = lab > 0
    a_all = area[wet]
    d_all = depth[wet].astype(float)
    total = float(a_all.sum())
    lat_all = np.broadcast_to(grid.lat_deg[:, None], lab.shape)[wet]
    bands = {}
    for lo, hi in LAT_BANDS:
        m = (np.abs(lat_all) >= lo) & (np.abs(lat_all) < hi)
        bands[f'{lo}-{hi}'] = dict(area_share=round(float(a_all[m].sum() / total), 4),
                                   volume_share=round(float((a_all[m] * d_all[m]).sum() / (a_all * d_all).sum()), 4))
    sizes = t['area'][1:] / 1e6
    by_size = {}
    for lo, hi in SIZE_CLASSES_KM2:
        m = (sizes >= lo) & (sizes < hi)
        by_size[f'{lo:g}-{hi:g}'] = dict(
            bodies=int(m.sum()), area_share=round(float(sizes[m].sum() / sizes.sum()), 4),
            volume_share=round(float(t['volume'][1:][m].sum() / t['volume'][1:].sum()), 4),
            median_deepest_m=round(float(np.median(t['deepest'][1:][m])), 1) if m.any() else None)
    small = (sizes < LARGE_KM2)
    w = t['area'][1:][small]
    crater_seas = dict(bodies=int(small.sum()), area_km2=round(float(w.sum()) / 1e6),
                       area_share=round(float(w.sum() / t['area'][1:].sum()), 4),
                       deepest_quantiles_by_area_m={str(q): round(weighted_quantile(t['deepest'][1:][small], w, q), 1)
                                                    for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
                       mean_depth_quantiles_by_area_m={str(q): round(weighted_quantile(
                           t['volume'][1:][small] / t['area'][1:][small], w, q), 1) for q in (0.1, 0.5, 0.9)},
                       abs_lat_quantiles_by_area={str(q): round(weighted_quantile(
                           np.abs(t['lat'][1:][small]), w, q), 1) for q in (0.1, 0.5, 0.9)})
    return dict(sea_level_m=sea_level, cells=int(a_all.size), bodies=int(lab.max()),
                area_km2=round(total / 1e6), area_share_of_moon=round(total / float(area.sum()), 5),
                volume_km3=round(float((a_all * d_all).sum()) / 1e9),
                mean_depth_m=round(float((a_all * d_all).sum() / total), 1),
                local_depth_quantiles_by_area_m={str(q): round(weighted_quantile(d_all, a_all, q), 1)
                                                 for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
                shelf_area_share={f'{s:g}': round(float(a_all[d_all < s].sum() / total), 5) for s in SHELF_DEPTHS_M},
                bands=bands, by_size_km2=by_size, crater_seas=crater_seas)


def basins_of(lab, regions, large, t, area, depth, grid):
    """Every basin of at least POLAR_MIN_KM2: the large bodies' regions and the other bodies whole."""
    cols = lab.shape[1]
    flat_reg = regions.ravel()
    out = []
    order = np.argsort(flat_reg, kind='stable')
    bounds = np.searchsorted(flat_reg[order], np.arange(int(regions.max()) + 2))
    for r in range(1, int(regions.max()) + 1):
        cells = order[bounds[r]:bounds[r + 1]]
        if not cells.size or area.ravel()[cells].sum() / 1e6 < POLAR_MIN_KM2:
            continue
        seed = cells[np.argmax(depth.ravel()[cells])]
        out.append(dict(kind='region', id=r, body=int(lab.ravel()[cells[0]]), cells=cells, seed=int(seed)))
    flat_lab = lab.ravel()
    for k in range(1, t['area'].size):
        if k in large or t['area'][k] / 1e6 < POLAR_MIN_KM2:
            continue
        out.append(dict(kind='body', id=k, body=k, cells=np.flatnonzero(flat_lab == k), seed=int(t['seed'][k])))
    for b in out:
        a = area.ravel()[b['cells']]
        lat = grid.lat_deg[b['cells'] // cols]
        b['area'] = float(a.sum())
        b['lat'] = float(np.average(lat, weights=a))
        b['abs_lat'] = float(np.average(np.abs(lat), weights=a))
        b['deepest'] = float(depth.ravel()[b['seed']])
    return out


def choose(basins, stats):
    chosen = {}
    for b in basins:
        if b['area'] / 1e6 >= RUN_MIN_KM2:
            chosen[(b['kind'], b['id'])] = 'a basin of at least 100,000 km2'
    for sign, pole in ((1, 'north'), (-1, 'south')):
        best = max((b for b in basins if sign * b['lat'] > 0), key=lambda b: b['abs_lat'])
        chosen.setdefault((best['kind'], best['id']), f'the basin nearest the {pole} pole')
    cs = stats['crater_seas']
    median_deep, median_lat = cs['deepest_quantiles_by_area_m']['0.5'], cs['abs_lat_quantiles_by_area']['0.5']
    pool = [b for b in basins if b['kind'] == 'body' and CRATER_SEA_KM2[0] <= b['area'] / 1e6 < CRATER_SEA_KM2[1]
            and abs(b['abs_lat'] - median_lat) <= CRATER_SEA_LAT_WINDOW]
    best = min(pool, key=lambda b: abs(b['deepest'] - median_deep))
    chosen.setdefault((best['kind'], best['id']), 'the crater seas\' representative')
    shallow = [b for b in basins if b['kind'] == 'region' and b['deepest'] < SHALLOW_SEA_DEEPEST_M]
    best = max(shallow, key=lambda b: b['area'])
    chosen.setdefault((best['kind'], best['id']), 'the shallow seas\' representative')
    return chosen


def main(argv=None) -> int:
    from scipy import ndimage
    from geography import basins as bs
    from geography import fetch_inputs
    from geography import nomenclature
    tic = time.time()
    h, level, area, grid, atlas = heights()
    depth = (level - h).astype(np.float32)
    _, lab = bs.seas(h, area, level, min_area_fraction=0.0)
    del h
    lab = numbered_by_area(lab, area)
    t = body_table(lab, depth, area, grid)
    print(f'{t["area"].size - 1} bodies in {time.time() - tic:.0f} s', flush=True)
    water, craters = features()
    stats = summary(lab, depth, area, grid, t, level)
    large = [k for k in range(1, t['area'].size) if t['area'][k] / 1e6 >= LARGE_KM2]
    bodies = []
    for k in range(1, min(t['area'].size, 41)):
        entry = dict(body=k, area_km2=round(float(t['area'][k]) / 1e6, 1),
                     volume_km3=round(float(t['volume'][k]) / 1e9, 1),
                     mean_depth_m=round(float(t['volume'][k] / t['area'][k]), 1),
                     deepest_m=round(float(t['deepest'][k]), 1),
                     lat_span=[round(float(t['lat_min'][k]), 2), round(float(t['lat_max'][k]), 2)],
                     centroid=dict(lat=round(float(t['lat'][k]), 3), lon_east=round(float(t['lon'][k]), 3)),
                     names=names_in(lab, k, water, grid)[:6])
        near = min(atlas['bodies'], key=lambda x: math.hypot(
            x['centroid'][0] - entry['centroid']['lat'],
            ((x['centroid'][1] - entry['centroid']['lon_east'] + 180.0) % 360.0) - 180.0))
        if abs(near['area_km2'] - entry['area_km2']) / max(near['area_km2'], 1.0) < 0.25:
            entry['atlas_area_km2'], entry['atlas_deepest_m'] = near['area_km2'], near['max_depth_m']
        bodies.append(entry)
    connections = {str(k): dict(body=k, names=names_in(lab, k, water, grid)[:6],
                                splits=splits(lab, k, depth, area, grid, water)) for k in large}
    print(f'splits at {time.time() - tic:.0f} s', flush=True)
    regions, steps = partition(lab, large, depth)
    print(f'partition in {steps} steps at {time.time() - tic:.0f} s', flush=True)
    basins = basins_of(lab, regions, large, t, area, depth, grid)
    chosen = choose(basins, stats)
    sills = Sills(lab, depth)
    used = set()
    run = []
    for b in sorted(basins, key=lambda b: -b['area']):
        if (b['kind'], b['id']) not in chosen:
            continue
        entry = describe(b['cells'], depth, area, grid, water, craters, regions if b['kind'] == 'region' else lab,
                         b['id'], b['seed'])
        entry['key'] = key_of(entry, used)
        entry.update(body=b['body'], basin_of=('a region of a sea over 100,000 km2' if b['kind'] == 'region'
                                               else 'a whole body'), chosen_because=chosen[(b['kind'], b['id'])])
        run.append((entry, b))
        print(f"  {entry['key']}: {entry['area_km2']:.0f} km2, {entry['deepest_m']:.0f} m, lat {entry['lat_span']}, "
              f"{entry['chosen_because']}", flush=True)
    for entry, b in run:
        entry['sills_m'] = {e['key']: sills(b['body'], b['seed'], x['seed'])
                            for e, x in run if x['body'] == b['body'] and x is not b}
        entry['deeper_parts'] = []
        if b['kind'] != 'region':
            continue
        inside = np.zeros(lab.shape, bool)
        inside.ravel()[b['cells']] = True
        for d in DEEPER_SPLITS_M:
            comp = label8(inside & (depth > d))
            cid = np.unique(comp[comp > 0])
            if not cid.size:
                break
            a = ndimage.sum(area, comp, cid) / 1e6
            big = cid[a >= SUB_BASIN_MIN_KM2]
            if big.size < 2:
                continue
            seeds = [int(np.ravel_multi_index(ndimage.maximum_position(depth, comp, c), lab.shape)) for c in big]
            parts = [dict(names=names_in(comp, c, water, grid)[:3], area_km2=round(float(a[cid == c][0])),
                          deepest_m=round(float(depth.ravel()[s]), 1)) for c, s in zip(big, seeds)]
            pairs = [dict(a=i, b=j, sill_m=sills(b['body'], seeds[i], seeds[j]))
                     for i in range(len(seeds)) for j in range(i + 1, len(seeds))]
            entry['deeper_parts'].append(dict(split_at_m=d, parts=parts, sills=pairs))
        del inside
    print(f'sills at {time.time() - tic:.0f} s', flush=True)
    # the seas' area each run basin stands for: its own; each other region of a large sea goes to the run basin of
    # the same sea nearest it in mean depth, and the bodies under 100,000 km2 to the crater seas' representative
    stands = {e['key']: dict(own_km2=e['area_km2'], other_regions_km2=0.0, crater_seas_km2=0.0) for e, _ in run}
    flat_reg, flat_lab = regions.ravel(), lab.ravel()
    reg_area = np.bincount(flat_reg, area.ravel())
    reg_mean = np.bincount(flat_reg, (area * depth).ravel()) / np.maximum(reg_area, 1.0)
    reg_body = np.zeros(reg_area.size, np.int64)
    reg_body[flat_reg] = flat_lab
    run_regions = {b['id'] for e, b in run if b['kind'] == 'region'}
    other = dict(area_km2=0.0, by_deepest_m={})
    for r in range(1, reg_area.size):
        if r in run_regions or reg_area[r] <= 0:
            continue
        same = [(e, b) for e, b in run if b['body'] == reg_body[r] and b['kind'] == 'region']
        e, _ = min(same, key=lambda eb: abs(eb[0]['mean_depth_m'] - reg_mean[r]))
        stands[e['key']]['other_regions_km2'] += reg_area[r] / 1e6
        other['area_km2'] += reg_area[r] / 1e6
    reg_deep = np.zeros(reg_area.size)
    np.maximum.at(reg_deep, flat_reg, depth.ravel())
    unrun = np.array([r for r in range(1, reg_area.size) if r not in run_regions and reg_area[r] > 0])
    for lo, hi in ((0, 500), (500, 1000), (1000, 2000), (2000, 1e4)):
        m = (reg_deep[unrun] >= lo) & (reg_deep[unrun] < hi)
        other['by_deepest_m'][f'{lo:g}-{hi:g}'] = dict(
            area_km2=round(float(reg_area[unrun][m].sum()) / 1e6),
            mean_depth_m=round(float((reg_mean[unrun][m] * reg_area[unrun][m]).sum() / max(reg_area[unrun][m].sum(),
                                                                                          1.0)), 1))
    other['area_km2'] = round(other['area_km2'])
    other['area_share_of_the_seas'] = round(other['area_km2'] / stats['area_km2'], 4)
    crater_key = next(e['key'] for e, b in run if chosen[(b['kind'], b['id'])] == 'the crater seas\' representative')
    run_bodies = {b['id'] for e, b in run if b['kind'] == 'body'}
    for k in range(1, t['area'].size):
        if k in large or k in run_bodies:
            continue
        stands[crater_key]['crater_seas_km2'] += t['area'][k] / 1e6
    total = sum(sum(v.values()) for v in stands.values())
    for e, _ in run:
        s = stands[e['key']]
        s.update({k: round(v, 1) for k, v in s.items()})
        s['share_of_the_seas'] = round(sum(v for k, v in s.items() if k.endswith('km2')) / stats['area_km2'], 5)
        e['stands_for'] = s
    regions_summary = {}
    for k in large:
        inside = [b for b in basins if b['body'] == k and b['kind'] == 'region']
        regions_summary[str(k)] = dict(regions=int(np.unique(regions[lab == k]).size),
                                       basins_over_10000_km2=len(inside),
                                       their_area_share=round(sum(b['area'] for b in inside) / float(t['area'][k]), 4))
    inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (ATLAS, DRAINAGE)}
    for name in ('ldem_16.img', 'ldem_4.img', 'jggrx_0420a_sha.tab', nomenclature.ARCHIVE):
        inputs[f'geography/inputs/{name}'] = sha256(fetch_inputs.path(name))
    drainage = json.loads(DRAINAGE.read_text())
    code = [Path(__file__), ROOT / 'geography' / 'topography.py', ROOT / 'geography' / 'basins.py',
            ROOT / 'geography' / 'nomenclature.py', ROOT / 'research' / 'studies' / 'lake_oxygen' / 'lakes.py']
    out = dict(
        schema=SCHEMA,
        producer=dict(files={str(p.resolve().relative_to(ROOT)): sha256(p) for p in code}, inputs=inputs),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        checks=dict(drainage_sea_level_m=drainage['sea_level_m'], atlas_sea_level_m=atlas['sea_level_m'],
                    atlas_share=atlas['share'], share_on_this_grid=stats['area_share_of_moon'],
                    stands_for_total_km2=round(total, 1), partition_growth_steps=steps),
        seas=stats, bodies=bodies, connections=connections, regions=regions_summary, other_regions=other,
        selection=dict(large_km2=LARGE_KM2, partition_m=PARTITION_M, run_min_km2=RUN_MIN_KM2,
                       polar_min_km2=POLAR_MIN_KM2, crater_sea_km2=list(CRATER_SEA_KM2),
                       crater_sea_lat_window=CRATER_SEA_LAT_WINDOW, rule=RULE),
        t21=dict(lat=[round(float(x), 6) for x in t21_lat()], lon=[i * 360.0 / T21_COLS for i in range(T21_COLS)],
                 assignment='each 16 pixel/degree cell goes to the T21 cell with the nearest Gaussian latitude and '
                            'the nearest longitude'),
        run_basins=[e for e, _ in run], seconds=round(time.time() - tic, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT), f'in {time.time() - tic:.0f} s')
    return 0


RULE = ('The four seas over 100,000 km2 are divided into basins at 250 m, each taking the shallower water nearest it '
        'through the sea; other bodies are their own basins. Run: every basin of at least 100,000 km2 (the largest '
        'seas and the deepest point); the basin nearest each pole among those of at least 10,000 km2; and the '
        'crater seas\' representative, the body of 10,000-100,000 km2 within 10 degrees of the crater seas\' median '
        'latitude whose deepest point lies nearest their median deepest point, both by area over the bodies under '
        '100,000 km2; and the shallow seas\' representative, the largest region of a large sea whose deepest point '
        'lies above 500 m.')
EVIDENCE = ('Hydrostatic geometry: LOLA at 16 pixels a degree above the GRAIL geoid (degree 200), as geography\'s '
            'drainage takes it, with every cell below the atlas\'s common level of the 28% share counted as sea and '
            'the bodies labelled by geography\'s basins.seas, joined across the seam and the poles. Depths are the '
            'level less the ground at 1.9 km cells; a sill is the deepest water path between two deepest points, '
            'to 2 m. The basins of the large seas are a partition at 250 m; their boundaries in the shallows follow '
            'the nearest water and are a construction. Crustal loading, sediment and shoreline change are not '
            'counted.')
READING_RULE = ('Depths in m below the common level, areas in km2, latitudes north, longitudes east. bodies lists '
                'the 40 largest bodies with the atlas\'s areas where they match. connections gives each sea over '
                '100,000 km2 split at a ladder of depths: the parts deeper than each depth and the largest of them. '
                'run_basins gives each basin the study runs: its hypsography (the share of its surface area deeper '
                'than each fraction of its deepest point), the share of its area shallower than each shelf depth, '
                'the T21 cells holding its water with their areas, the sills to the other run basins of its sea, '
                'the sills between its own deeper parts where it is a region of a large sea, and the seas\' area it '
                'stands for (its own, the other regions of its sea nearest it in depth, and for the crater seas\' '
                'representative every body under 100,000 km2 not run).')


if __name__ == '__main__':
    sys.exit(main())
