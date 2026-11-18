"""The atlas's regions and districts, a first draft on the design surface (17 November 2026).

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 <python with scipy> -m research.studies.atlas_round.districts

The register defines a district by significance and wholeness: a region that reads as a whole and meaningful place and
stands apart from the others (16 November). districts.md beside this file turns that into rules; this module draws a first
draft by them on the design surface, today's grid at 8 pixels per degree (3.8 km):

- Centres (rule 1): the named craters of FEATURE_KM or more on land, the larger lakes, and fill-in points about
  SEED_KM apart where neither stands; every island gets one.
- Edges (rules 2 and 9): districts grow from their centres over the ground by a watershed, so that they meet on its
  ridges and rim crests; crossing the edge of a moisture zone (the moisture map's classes, smoothed as the climate's
  fields are) or of the mare basalt costs as much as climbing MOISTURE_EDGE_M or ROCK_EDGE_M. The seas bound them.
- Character and standing apart (rules 3 and 8): neighbours that share their moisture zone, their rock and the waters
  they drain to, with a boundary standing less than WEAK_EDGE_M above the higher of their mean ground, merge while the
  whole stays within MAX_KM2.
- Size (rule 11): districts under MIN_KM2 join the neighbour they share the most edge with; districts over MAX_KM2 split
  round new centres along the ground's own ridges. Islands under MIN_KM2 stay whole down to ISLET_KM2.
- Regions (rule 10): the same watershed from the basins of REGION_FEATURE_KM or more and fill-in points about
  REGION_SEED_KM apart, with each district given to the region holding most of it, so that districts never cross a
  region's edge; each sea is a region of its own.
- Coastal seas (the author's choice of 17 November): the sea within COAST_KM of land is cut into coastal sea districts,
  each the water nearest to a stretch of shore, joined up to COAST_MAX_KM2, and linked to the land districts that face it
  and to those whose rivers reach it (the hybrid the author chose for satoumi); the open water beyond is its sea's
  region without districts.

Codes stand in for names until the names come from the Terlunan language's roots (the author's choice): regions L01...
on land and S1... at sea, districts by region and number. Descriptions are short and computed from the facts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from climate.gcm.turned import to_turned
from geography import mare_map
from geography import nomenclature as nm
from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import design_cuts as DC
from research.studies.atlas_round import globe_surface as GS
from research.studies.atlas_round import moisture_final as MF

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'districts.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/districts')
PRODUCT = BULK / 'districts_8ppd.npz'
SCHEMA = 'terluna.research.atlas-round-districts/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-districts-grid/1'
FALLS = HERE / 'results' / 'falls.json'
HERITAGE = ROOT / 'research' / 'studies' / 'conservation' / 'results' / 'heritage_register.json'
CODE = ('research/studies/atlas_round/districts.py', 'research/studies/atlas_round/district_names.py')

PPD = 8
MOON_KM = 1737.4
MIN_KM2, MAX_KM2, ISLET_KM2 = 5_000.0, 150_000.0, 3_000.0
FEATURE_KM = 60.0                    # named craters whose floors centre districts
REGION_FEATURE_KM = 300.0            # basins whose floors centre regions
SEED_KM, REGION_SEED_KM = 200.0, 750.0
LAKE_SEED_KM2 = 1_000.0
FLOOR_SHARE = 0.3                    # a centre's disc reaches this share of its feature's diameter from the middle
MOISTURE_EDGE_M, ROCK_EDGE_M = 600.0, 400.0
WEAK_EDGE_M = 200.0
COAST_KM, COAST_MAX_KM2 = 60.0, 40_000.0
STRAIT_MASS_KM2 = 20_000.0           # land masses this large or larger make a coastal sea a bay or a strait
REGION_MIN_KM2, REGION_MAX_KM2 = 250_000.0, 1_200_000.0
SMOOTH_DEG = 0.7 * 5.625             # the moisture map is smoothed over 0.7 of the climate run's T21 spacing
SALT_SEAS = ('South Pole–Aitken Sea',)
# The gazetteer's named maria, oceans, bays, lakes and marshes each lend their name's meaning to one unit (rule 13,
# districts.md); a feature whose centre lies on open water or in a unit already taken goes to the nearest free unit
# within FEATURE_REACH of its diameter, FEATURE_REACH_KM at most. Features the names table pins by name join them.
FEATURE_TYPES = ('Mare, maria', 'Oceanus, oceani', 'Sinus, sinūs', 'Lacus, lacūs', 'Palus, paludes')
FEATURE_REACH, FEATURE_REACH_KM = 0.25, 150.0


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path) -> str:
    p = Path(path)
    return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)


# ----------------------------------------------------------------------------------------------- the grid

def grid(ppd=PPD):
    rows, cols = 180 * ppd, 360 * ppd
    lat = 90.0 - (np.arange(rows) + 0.5) / ppd
    lon = (np.arange(cols) + 0.5) / ppd
    area = (MOON_KM * np.radians(1.0 / ppd)) ** 2 * np.cos(np.radians(lat))
    return lat, lon, area


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def to_latlon(v):
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.degrees(np.arcsin(np.clip(v[..., 2], -1, 1))), np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360.0


def neighbours_differ(a):
    """Cells whose value differs from a 4-neighbour's (longitude wrapping)."""
    d = np.zeros(a.shape, bool)
    d[:-1] |= a[:-1] != a[1:]
    d[1:] |= a[1:] != a[:-1]
    d |= a != np.roll(a, 1, axis=1)
    d |= a != np.roll(a, -1, axis=1)
    return d


# ------------------------------------------------------------------------------------------------- inputs

def surface_classes(log=print):
    """The design surface at PPD: ground height, and each cell land (0), lake (1) or sea, with the sea's body index; the
    surface's metadata. Each 2 x 2 block of the 16 px/deg surface is sea where at least two of its cells are, else lake
    where at least two are, else land."""
    with np.load(GS.PRODUCT, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        height = z['height_m'].astype(np.float32)
        cells, body = z['water_cell'], z['water_body']
    rows16, cols16 = height.shape
    b16 = np.full(rows16 * cols16, -1, np.int8)
    b16[cells] = body
    b16 = b16.reshape(rows16, cols16)
    k = rows16 // (180 * PPD)
    blocks = lambda a: a.reshape(rows16 // k, k, cols16 // k, k)
    sea_n = blocks(b16 >= 1).sum(axis=(1, 3))
    lake_n = blocks(b16 == 0).sum(axis=(1, 3))
    kind = np.where(sea_n >= 2, 2, np.where(lake_n >= 2, 1, 0)).astype(np.int8)
    sea_body = np.zeros(kind.shape, np.int8)
    for b in range(1, len(meta['bodies'])):
        n = blocks(b16 == b).sum(axis=(1, 3))
        sea_body = np.where((kind == 2) & (n >= 2) & (sea_body == 0), b, sea_body)
    sea_body = np.where((kind == 2) & (sea_body == 0), 1, sea_body)
    h = blocks(height).mean(axis=(1, 3))
    log(f'surface at {PPD} px/deg: land {np.mean(kind == 0):.3f}, lakes {np.mean(kind == 1):.3f}, seas {np.mean(kind == 2):.3f} of the cells')
    return h, kind, sea_body, meta, b16


def routing_fields(b16, log=print):
    """From the design routing at 16 px/deg: the water each cell finally drains to (the surface's body index of its
    terminal cell, 0 for a lake, -1 for a dry hollow) and the flow through it, both taken at PPD as the block's
    largest flow and that cell's end."""
    with np.load(DR.PRODUCT, allow_pickle=False) as z:
        receiver = z['receiver'].astype(np.int64)
        q = z['discharge_m3s'].astype(np.float32)
    end = receiver.copy()
    for _ in range(40):
        nxt = end[end]
        if np.array_equal(nxt, end):
            break
        end = nxt
    flat = b16.ravel()
    to = np.where(flat[end] >= 0, flat[end], -1).astype(np.int8)
    rows16, cols16 = b16.shape
    k = rows16 // (180 * PPD)
    qb = q.reshape(rows16 // k, k, cols16 // k, k)
    q8 = qb.max(axis=(1, 3))
    arg = qb.transpose(0, 2, 1, 3).reshape(rows16 // k, cols16 // k, k * k).argmax(axis=2)
    r = np.arange(rows16 // k)[:, None] * k + arg // k
    c = np.arange(cols16 // k)[None, :] * k + arg % k
    to8 = to.reshape(rows16, cols16)[r, c]
    log('routing: ends found for every cell')
    return to8, q8, end, to


def moisture_classes(rotation, log=print):
    """The moisture map's classes on today's grid: each cell's place in the run's turned frame, and there the log
    aridity of the run's land cells (a water cell taking its nearest land cell's, as the map assigns them) averaged with
    Gaussian weights over SMOOTH_DEG, then classed by the map's bounds."""
    from scipy.spatial import cKDTree
    record = json.loads(MF.OUT.read_text())
    lat_g = np.asarray(record['grid']['lat'], float)
    lon_g = np.asarray(record['grid']['lon'], float)
    cells = record['cells']
    ln = np.full((lat_g.size, lon_g.size), np.nan)
    ln[np.asarray(cells['row']), np.asarray(cells['col'])] = np.log(np.asarray(cells['aridity'], float))
    pts = unit(*np.meshgrid(lat_g, lon_g, indexing='ij'))
    land = np.isfinite(ln)
    _, nearest = cKDTree(pts[land]).query(pts[~land])
    ln[~land] = ln[land][nearest]
    tree = cKDTree(pts.reshape(-1, 3))
    lat, lon, _ = grid()
    out = np.empty((lat.size, lon.size), np.int8)
    bounds = [c['aridity'][1] for c in sorted(record['classes'], key=lambda c: c['index'])][:-1]
    sigma = np.radians(SMOOTH_DEG)
    values = ln.ravel()
    for r0 in range(0, lat.size, 120):
        la, lo = np.meshgrid(lat[r0:r0 + 120], lon, indexing='ij')
        lat_t, lon_t = to_turned(la, lo, rotation)
        p = unit(lat_t, lon_t).reshape(-1, 3)
        d, i = tree.query(p, k=12)
        w = np.exp(-0.5 * (d / sigma) ** 2)
        a = np.exp((w * values[i]).sum(axis=1) / w.sum(axis=1))
        out[r0:r0 + 120] = np.searchsorted(bounds, a, side='right').reshape(la.shape)
    names = [c['name'] for c in sorted(record['classes'], key=lambda c: c['index'])]
    log(f'moisture zones on today\'s grid: {", ".join(f"{n} {np.mean(out == i):.2f}" for i, n in enumerate(names))}')
    return out, names


# ------------------------------------------------------------------------------------------------ centres

def disc_cells(lat0, lon0, radius_km, lat, lon):
    """Rows and columns of the grid's cells within radius_km of a place (longitude wrapping)."""
    ang = radius_km / MOON_KM
    r0 = max(0, int((90.0 - lat0 - np.degrees(ang)) * PPD) - 1)
    r1 = min(lat.size, int((90.0 - lat0 + np.degrees(ang)) * PPD) + 2)
    width = np.degrees(ang) / max(np.cos(np.radians(min(abs(lat0) + np.degrees(ang), 89.9))), 1e-3)
    if width >= 180.0:
        cols = np.arange(lon.size)
    else:
        c0 = int(((lon0 - width) % 360.0) * PPD) - 1
        cols = (c0 + np.arange(int(2 * width * PPD) + 4)) % lon.size
    rr, cc = np.meshgrid(np.arange(r0, r1), cols, indexing='ij')
    v = unit(lat[rr], lon[cc])
    inside = v @ unit(np.array(lat0), np.array(lon0)) >= np.cos(ang)
    return rr[inside], cc[inside]


def fibonacci(spacing_km):
    n = int(round(4 * np.pi * MOON_KM ** 2 / (0.866 * spacing_km ** 2)))
    i = np.arange(n) + 0.5
    lat = np.degrees(np.arcsin(1 - 2 * i / n))
    lon = (np.degrees(np.pi * (1 + 5 ** 0.5) * i)) % 360.0
    return lat, lon


def place_centres(on_land, kind, lat, lon, min_km, seed_km, lake_km2=None, log=print, islands=True):
    """Markers for the watershed: discs on the floors of named craters of min_km or more (the smaller painted over the
    larger, so that a crater inside a basin keeps its own), the lakes of lake_km2 or more not mostly under a disc,
    fill-in points about seed_km apart farther than 0.6 seed_km from any marker, and one point on each island left
    without one. Returns the marker grid and each marker's (kind, diameter or area, latitude, longitude)."""
    from scipy import ndimage
    from scipy.spatial import cKDTree
    markers = np.zeros(on_land.shape, np.int32)
    info = [None]
    feats = [f for f in nm.features() if f.get('diameter_km', 0) >= min_km
             and f.get('type') in ('Crater, craters', 'Satellite Feature')]
    feats.sort(key=lambda f: -f['diameter_km'])
    for f in feats:
        lo = f['lon'] % 360.0
        r, c = int((90.0 - f['lat']) * PPD), int(lo * PPD) % lon.size
        if not on_land[min(r, lat.size - 1), c]:
            continue
        rr, cc = disc_cells(f['lat'], lo, FLOOR_SHARE * f['diameter_km'], lat, lon)
        keep = on_land[rr, cc]
        if not keep.any():
            continue
        info.append(('crater', float(f['diameter_km']), float(f['lat']), lo, f['name']))
        markers[rr[keep], cc[keep]] = len(info) - 1
    if lake_km2:
        lakes, n = ndimage.label(kind == 1, structure=np.ones((3, 3)))
        _, _, area = grid()
        sizes = ndimage.sum(np.broadcast_to(area[:, None], kind.shape), lakes, np.arange(1, n + 1))
        for i in np.flatnonzero(sizes >= lake_km2) + 1:
            m = lakes == i
            if (markers[m] > 0).mean() > 0.5:
                continue
            rr, cc = np.nonzero(m)
            v = unit(lat[rr], lon[cc]).mean(axis=0)
            la, lo = to_latlon(v)
            info.append(('lake', float(sizes[i - 1]), float(la), float(lo), None))
            markers[m & (markers == 0)] = len(info) - 1
    taken = np.argwhere(markers > 0)
    tree = cKDTree(unit(lat[taken[:, 0]], lon[taken[:, 1]])) if taken.size else None
    sl, so = fibonacci(seed_km)
    r = np.clip(((90.0 - sl) * PPD).astype(int), 0, lat.size - 1)
    c = (so * PPD).astype(int) % lon.size
    ok = on_land[r, c] & (markers[r, c] == 0)
    if tree is not None:
        d, _ = tree.query(unit(lat[r], lon[c]))
        ok &= d * MOON_KM > 0.6 * seed_km
    for rr, cc in zip(r[ok], c[ok]):
        info.append(('ground', 0.0, float(lat[rr]), float(lon[cc]), None))
        markers[rr, cc] = len(info) - 1
    comps, n = ndimage.label(on_land, structure=np.ones((3, 3)))
    comps = join_wrapped(comps)
    have = np.unique(comps[markers > 0])
    for i in (np.setdiff1d(np.unique(comps[comps > 0]), have) if islands else []):
        rr, cc = np.nonzero(comps == i)
        v = unit(lat[rr], lon[cc]).mean(axis=0)
        k = int(np.argmax(unit(lat[rr], lon[cc]) @ v))
        info.append(('island', 0.0, float(lat[rr[k]]), float(lon[cc[k]]), None))
        markers[rr[k], cc[k]] = len(info) - 1
    kinds = [x[0] for x in info[1:]]
    log(f'centres: {kinds.count("crater")} craters of {min_km:.0f} km or more, {kinds.count("lake")} lakes, '
        f'{kinds.count("ground")} fill-in points, {kinds.count("island")} islands')
    return markers, info


def join_wrapped(labels):
    """Labels of components that touch across the grid's seam made one."""
    left, right = labels[:, 0], labels[:, -1]
    pairs = {(a, b) for a, b in zip(left, right) if a and b and a != b}
    if not pairs:
        return labels
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            x = parent[x]
        return x
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    lut = np.arange(labels.max() + 1)
    for x in parent:
        lut[x] = find(x)
    return lut[labels]


# ---------------------------------------------------------------------------------------------- watershed

def cost_surface(h, kind, moisture, basalt):
    """The ground the districts grow over: height in metres, more where a moisture zone or the mare basalt changes; the
    seas a barrier."""
    c = h - h[kind < 2].min() + 1.0
    c = c + MOISTURE_EDGE_M * neighbours_differ(moisture) + ROCK_EDGE_M * neighbours_differ(basalt)
    c = np.clip(c, 1, 65000)
    c[kind == 2] = 65535
    return c.astype(np.uint16)


def watershed(cost, markers, pad=None):
    """scipy's watershed by image foresting, with the grid wrapped in longitude."""
    from scipy import ndimage
    pad = pad or cost.shape[1] // 10
    cp = np.concatenate([cost[:, -pad:], cost, cost[:, :pad]], axis=1)
    mp = np.concatenate([markers[:, -pad:], markers, markers[:, :pad]], axis=1).astype(np.int32)
    out = ndimage.watershed_ift(cp, mp, structure=np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]))
    return out[:, pad:pad + cost.shape[1]]


# --------------------------------------------------------------------------------------------- the rules

def stats(labels, h, moisture, basalt, to8, area):
    """Per label: area, mean ground, the moisture zones' areas, the basalt share and the area draining to each water."""
    n = int(labels.max()) + 1
    w = np.broadcast_to(area[:, None], labels.shape)
    lab = labels.ravel()
    s = dict(area=np.bincount(lab, w.ravel(), n))
    s['h'] = np.bincount(lab, (h * w).ravel(), n) / np.maximum(s['area'], 1e-9)
    s['moist'] = np.bincount(lab * 6 + moisture.ravel().astype(np.int64), w.ravel(), n * 6).reshape(n, 6)
    s['basalt'] = np.bincount(lab, (basalt * w).ravel(), n) / np.maximum(s['area'], 1e-9)
    s['to'] = np.bincount(lab * 9 + (to8.ravel().astype(np.int64) + 1), w.ravel(), n * 9).reshape(n, 9)
    return s


def adjacency(labels, h):
    """Pairs of neighbouring labels (both above zero) with their shared edge's length in cell sides and the mean height
    of the higher cell along it."""
    pairs = []
    for a, b, ha, hb in ((labels[:, :], np.roll(labels, -1, axis=1), h, np.roll(h, -1, axis=1)),
                         (labels[:-1], labels[1:], h[:-1], h[1:])):
        m = (a != b) & (a > 0) & (b > 0)
        lo, hi = np.minimum(a[m], b[m]), np.maximum(a[m], b[m])
        pairs.append(np.stack([lo, hi, np.maximum(ha[m], hb[m])], axis=1))
    p = np.concatenate(pairs)
    key = p[:, 0].astype(np.int64) * (int(labels.max()) + 1) + p[:, 1].astype(np.int64)
    uk, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    hsum = np.bincount(inv, p[:, 2])
    n1 = int(labels.max()) + 1
    return np.stack([uk // n1, uk % n1], axis=1), cnt, hsum / cnt


def relabel(labels, lut):
    out = lut[labels]
    keep = np.unique(out[out > 0])
    remap = np.zeros(int(out.max()) + 1, np.int32)
    remap[keep] = np.arange(1, keep.size + 1)
    return remap[out]


def merge_pass(labels, h, moisture, basalt, to8, area, island_of, log=print):
    """One round of the size and standing-apart rules: each label merges at most once."""
    s = stats(labels, h, moisture, basalt, to8, area)
    pairs, length, edge_h = adjacency(labels, h)
    n = int(labels.max()) + 1
    lut = np.arange(n)
    used = np.zeros(n, bool)
    a_ = s['area']
    # too small: into the neighbour sharing the longest edge among those the join keeps within MAX_KM2, else the longest
    small = [i for i in range(1, n) if a_[i] < MIN_KM2]
    best = {}
    for (a, b), L in zip(pairs, length):
        for x, y in ((a, b), (b, a)):
            if a_[x] < MIN_KM2:
                fits = a_[x] + a_[y] <= MAX_KM2
                score = (fits, L)
                if score > best.get(x, ((False, 0), 0))[0]:
                    best[x] = (score, y)
    best = {x: (s[1], y) for x, (s, y) in best.items()}
    merged = 0
    for x in sorted(small, key=lambda i: a_[i]):
        if x in best and not used[x] and not used[best[x][1]]:
            lut[x] = best[x][1]
            used[x] = used[best[x][1]] = True
            merged += 1
    # neighbours alike, with a weak edge between them
    zone = s['moist'].argmax(axis=1)
    rock = s['basalt'] >= 0.5
    water = s['to'].argmax(axis=1)
    weak = []
    for (a, b), L, eh in zip(pairs, length, edge_h):
        if used[a] or used[b]:
            continue
        rise = eh - max(s['h'][a], s['h'][b])
        if (zone[a] == zone[b] and rock[a] == rock[b] and water[a] == water[b] and rise < WEAK_EDGE_M
                and a_[a] + a_[b] <= MAX_KM2 and island_of[a] == island_of[b]):
            weak.append((rise, a, b))
    alike = 0
    for rise, a, b in sorted(weak):
        if not used[a] and not used[b]:
            lut[b] = a
            used[a] = used[b] = True
            alike += 1
    log(f'  merge pass: {n - 1} districts, {merged} too small joined, {alike} alike pairs joined')
    return relabel(labels, lut), merged + alike


CLIMB_M = 100.0                      # in a split, each 100 m climbed between cells counts as one more cell's length


def geodesic_parts(cells_r, cells_c, h, lat, seeds):
    """Each of a district's cells given to the seed it is nearest by a path over its cells (4-neighbours, longitude
    wrapping), a step costing its length on the ground times 1 + the height climbed or descended over CLIMB_M, so that
    parts meet on the steep lines between their centres."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    cols = h.shape[1]
    index = {(r, c): k for k, (r, c) in enumerate(zip(cells_r.tolist(), cells_c.tolist()))}
    a, b, w = [], [], []
    step_ns = MOON_KM * np.radians(1.0 / PPD)
    for k, (r, c) in enumerate(zip(cells_r.tolist(), cells_c.tolist())):
        for rr, cc, step in ((r, (c + 1) % cols, step_ns * np.cos(np.radians(lat[r]))), (r + 1, c, step_ns)):
            j = index.get((rr, cc))
            if j is not None:
                a.append(k)
                b.append(j)
                w.append(step * (1.0 + abs(float(h[r, c]) - float(h[rr, cc])) / CLIMB_M))
    n = len(cells_r)
    graph = coo_matrix((w, (a, b)), shape=(n, n)).tocsr()
    _, _, source = dijkstra(graph, directed=False, indices=seeds, return_predecessors=True, min_only=True)
    position = {s: j for j, s in enumerate(seeds)}
    return np.array([position.get(int(s), 0) for s in source])


def split_large(labels, cost, area, lat, lon, log=print, h=None):
    """Districts over MAX_KM2 split into parts of about 0.6 MAX_KM2: centres spread by k-means on the sphere, and each
    cell given to its nearest centre over the ground (geodesic_parts), so that the parts come out balanced and meet on
    the district's walls and scarps."""
    s_area = np.bincount(labels.ravel(), np.broadcast_to(area[:, None], labels.shape).ravel())
    big = [i for i in range(1, int(labels.max()) + 1) if s_area[i] > MAX_KM2]
    if not big:
        return labels, 0
    out = labels.copy()
    nxt = int(labels.max()) + 1
    for i in big:
        rr, cc = np.nonzero(labels == i)
        v = unit(lat[rr], lon[cc])
        k = max(2, int(np.ceil(s_area[i] / (0.6 * MAX_KM2))))
        rng = np.random.default_rng(i)
        centres = v[rng.choice(len(v), k, replace=False)]
        for _ in range(20):
            g = np.argmax(v @ centres.T, axis=1)
            centres = np.stack([v[g == j].mean(axis=0) if (g == j).any() else centres[j] for j in range(k)])
            centres /= np.linalg.norm(centres, axis=1, keepdims=True)
        seeds = sorted({int(np.argmax(v @ c)) for c in centres})
        part = geodesic_parts(rr, cc, h, lat, seeds)
        for j in range(1, len(seeds)):
            sel = part == j
            out[rr[sel], cc[sel]] = nxt
            nxt += 1
    log(f'  split {len(big)} districts over {MAX_KM2:,.0f} km2')
    return relabel(out, np.arange(int(out.max()) + 1)), len(big)


def apply_rules(labels, h, moisture, basalt, to8, area, cost, lat, lon, islands, log=print):
    for round_ in range(12):
        changed = 0
        while True:
            island_of = np.zeros(int(labels.max()) + 1, np.int64)
            island_of[labels.ravel()] = islands.ravel()
            labels, n = merge_pass(labels, h, moisture, basalt, to8, area, island_of, log)
            changed += n
            if n == 0:
                break
        labels, n = split_large(labels, cost, area, lat, lon, log, h=h)
        if n == 0:
            break
    return labels


def island_territories(on_land, area, lat, lon, log=print):
    """Each landmass's territory (rule 15, districts.md). A landmass of REGION_MIN_KM2 or more is a mainland, its own
    territory. An island within COAST_KM of a mainland's shore, or of an island already in a territory, lies in that
    territory's waters and joins it, nearest first, so that a chain of islands joins with its first link. The islands
    left lie in international waters; those whose waters touch (within twice COAST_KM of each other) form a group, each
    group of ISLET_KM2 or more a territory of its own. A smaller group is rocks: its islets join the nearest district
    and claim no waters (territory -1). Returns the territory of every land cell (0 at sea), the landmasses (the
    component grid), the order in which islands joined and the shore cell each joined through, the groups'
    territories and the rocks' cells."""
    from scipy import ndimage
    from scipy.spatial import cKDTree
    comps, _ = ndimage.label(on_land, structure=np.ones((3, 3)))
    comps = join_wrapped(comps)
    n = int(comps.max()) + 1
    size = np.bincount(comps.ravel(), np.broadcast_to(area[:, None], comps.shape).ravel(), n)
    ids = np.unique(comps[comps > 0])
    terr = np.zeros(n, np.int64)
    for c in ids:
        if size[c] >= REGION_MIN_KM2:
            terr[c] = c
    shore = on_land & neighbours_differ(on_land)
    sr, sc = np.nonzero(shore)
    scomp = comps[sr, sc]
    v = unit(lat[sr], lon[sc])
    chord = lambda km: 2.0 * np.sin(km / MOON_KM / 2.0)                          # noqa: E731
    joined, through = [], {}
    pending = [c for c in ids if not terr[c]]
    while pending:
        have = np.flatnonzero(terr[scomp] > 0)
        if not have.size:
            break
        tree = cKDTree(v[have])
        moved = []
        for c in pending:
            mine = np.flatnonzero(scomp == c)
            d, k = tree.query(v[mine], distance_upper_bound=chord(COAST_KM))
            ok = np.isfinite(d)
            if ok.any():
                j = have[k[ok][np.argmin(d[ok])]]
                moved.append((float(d[ok].min()), c, j))
        if not moved:
            break
        for _, c, j in sorted(moved):
            terr[c] = terr[scomp[j]]
            joined.append(int(c))
            through[int(c)] = (int(sr[j]), int(sc[j]))
        pending = [c for c in pending if not terr[c]]
    groups = []
    if pending:
        parent = {c: c for c in pending}

        def root(c):
            while parent[c] != c:
                parent[c] = parent[parent[c]]
                c = parent[c]
            return c

        pts, owner = [], []
        for c in pending:                                   # at most a few hundred shore points per island
            mine = np.flatnonzero(scomp == c)
            mine = mine[::max(1, mine.size // 300)]
            pts.append(v[mine])
            owner.append(np.full(mine.size, c))
        pts, owner = np.concatenate(pts), np.concatenate(owner)
        for a, b in cKDTree(pts).query_pairs(chord(2.0 * COAST_KM), output_type='ndarray'):
            ra, rb = root(int(owner[a])), root(int(owner[b]))
            if ra != rb:
                parent[rb] = ra
        roots = sorted({root(c) for c in pending})
        for g, r in enumerate(roots):
            members = [c for c in pending if root(c) == r]
            if sum(size[c] for c in members) < ISLET_KM2:
                for c in members:
                    terr[c] = -1                        # rocks
                continue
            tid = n + g
            groups.append(tid)
            for c in members:
                terr[c] = tid
    territory = terr[comps]
    territory[~on_land] = 0
    rocks = territory == -1
    log(f'territories: {sum(1 for c in ids if size[c] >= REGION_MIN_KM2)} mainlands, {len(joined)} islands in their '
        f'waters, {len(groups)} island groups in international waters, {int(np.unique(comps[rocks]).size)} rocks')
    return territory, comps, joined, through, groups, rocks


def absorb_islets(labels, on_land, area, lat, lon, territory=None):
    """Islands under ISLET_KM2 join the nearest district across the water within their own territory (rule 15); in an
    island group of islets alone, the islets together make one district."""
    from scipy import ndimage
    from scipy.spatial import cKDTree
    comps, n = ndimage.label(on_land, structure=np.ones((3, 3)))
    comps = join_wrapped(comps)
    w = np.broadcast_to(area[:, None], labels.shape)
    sizes = np.bincount(comps.ravel(), w.ravel())
    tiny = np.flatnonzero(sizes < ISLET_KM2)
    tiny = tiny[tiny > 0]
    if not tiny.size:
        return labels
    is_tiny = np.isin(comps, tiny)
    out = labels.copy()
    if territory is None:
        territory = np.where(on_land, 1, 0)
    nxt = int(labels.max()) + 1
    for t in np.unique(territory[is_tiny]):
        mine = is_tiny & (territory == t)
        big = on_land & ~is_tiny & ((territory == t) if t != -1 else True)     # rocks join the nearest district
        tr, tc = np.nonzero(mine)
        if big.any():
            br, bc = np.nonzero(big)
            _, i = cKDTree(unit(lat[br], lon[bc])).query(unit(lat[tr], lon[tc]))
            out[tr, tc] = labels[br[i], bc[i]]
        else:
            out[tr, tc] = nxt
            nxt += 1
    return out


# ------------------------------------------------------------------------------------------- coastal seas

def coastal_seas(labels, kind, sea_body, h, area, lat, lon, log=print, region_of=None, rocks=None):
    """Coastal sea districts, the territorial waters (rule 12): the sea within COAST_KM of land, each cell with the
    land district its nearest shore belongs to, the pieces joined (within one sea and one region) into districts of
    MIN_KM2 to COAST_MAX_KM2. Rocks claim no waters. Returns the pieces, the facing districts, each piece's sea and its
    region."""
    from scipy.spatial import cKDTree
    land = kind < 2
    shore = land & neighbours_differ(land)
    if rocks is not None:
        shore &= ~rocks
    sr, sc = np.nonzero(shore)
    tree = cKDTree(unit(lat[sr], lon[sc]))
    wr, wc = np.nonzero(kind == 2)
    d, i = tree.query(unit(lat[wr], lon[wc]))
    km = 2 * np.arcsin(np.clip(d / 2, 0, 1)) * MOON_KM
    near = km <= COAST_KM
    faces = np.zeros(labels.shape, np.int32)
    faces[wr[near], wc[near]] = labels[sr[i[near]], sc[i[near]]]
    piece = np.zeros(labels.shape, np.int64)
    key = sea_body.astype(np.int64) * 100_000 + faces
    piece[faces > 0] = key[faces > 0]
    uk, inv = np.unique(piece[faces > 0], return_inverse=True)
    pieces = np.zeros(labels.shape, np.int32)
    pieces[faces > 0] = inv + 1
    body_of = np.concatenate([[0], uk // 100_000])
    reg_of = np.concatenate([[0], region_of[uk % 100_000]]) if region_of is not None else np.zeros(body_of.shape, np.int64)
    for _ in range(30):
        w = np.broadcast_to(area[:, None], labels.shape)
        a_ = np.bincount(pieces.ravel(), w.ravel())
        pairs, length, _ = adjacency(pieces, h)
        n = int(pieces.max()) + 1
        lut = np.arange(n)
        used = np.zeros(n, bool)
        order = sorted(zip(pairs, length), key=lambda x: min(a_[x[0][0]], a_[x[0][1]]))
        joined = 0
        for (a, b), L in order:
            if used[a] or used[b] or body_of[a] != body_of[b] or reg_of[a] != reg_of[b]:
                continue
            if min(a_[a], a_[b]) < MIN_KM2 or a_[a] + a_[b] <= COAST_MAX_KM2 * 0.5:
                lut[b] = a
                used[a] = used[b] = True
                joined += 1
        if not joined:
            break
        keep_body, keep_reg = body_of.copy(), reg_of.copy()
        pieces = lut[pieces]
        keep = np.unique(pieces[pieces > 0])
        remap = np.zeros(n, np.int32)
        remap[keep] = np.arange(1, keep.size + 1)
        body_of = np.concatenate([[0], keep_body[keep]])
        reg_of = np.concatenate([[0], keep_reg[keep]])
        pieces = remap[pieces]
    log(f'coastal seas: {int(pieces.max())} districts within {COAST_KM:.0f} km of the shore')
    return pieces, faces, body_of, reg_of


# ------------------------------------------------------------------------------------------------ regions

def regions(labels, on_land, kind, cost, area, lat, lon, log=print, territory=None, comps=None, joined=(), through=None,
            groups=()):
    """Land regions by the same watershed from the great basins and wide fill-in points, each district given to the
    region holding most of it; regions under REGION_MIN_KM2 join the neighbour they share the most edge with. With the
    territories (island_territories, rule 15), an island in a territory's waters takes the region of the land it joined
    through, and each island group in international waters makes a region of its own, exempt from the minimum.
    Returns each district's region, the centres' details and the island groups' regions."""
    # without island markers, an island falls to a region across the water as the flood crosses the sea last
    markers, info = place_centres(on_land, kind, lat, lon, REGION_FEATURE_KM, REGION_SEED_KM, log=log, islands=False)
    reg = watershed(cost, markers)
    reg[~on_land] = 0
    n = int(labels.max()) + 1
    w = np.broadcast_to(area[:, None], labels.shape).ravel()
    m = int(reg.max()) + 1
    table = np.bincount(labels.ravel().astype(np.int64) * m + reg.ravel(), w, n * m).reshape(n, m)
    table[:, 0] = 0
    of = table.argmax(axis=1)
    of[0] = 0
    exempt = set()
    if territory is not None:
        def majority(values):
            key = labels[on_land].astype(np.int64) * (int(values.max()) + 1) + values[on_land]
            uk, cnt = np.unique(key, return_counts=True)
            best = {}
            for k_, c_ in zip(uk, cnt):
                d_, v_ = divmod(int(k_), int(values.max()) + 1)
                if c_ > best.get(d_, (0, 0))[0]:
                    best[d_] = (c_, v_)
            out = np.zeros(n, np.int64)
            for d_, (_, v_) in best.items():
                out[d_] = v_
            return out
        t_of, c_of = majority(territory), majority(comps)
        nxt = int(of.max()) + 1
        for g in groups:                                    # an island group: a region of its own
            ds = np.flatnonzero(t_of == g)
            ds = ds[ds > 0]
            if ds.size:
                of[ds] = nxt
                exempt.add(nxt)
                nxt += 1
        for c in joined:                                    # an island in a territory's waters, nearest first
            ds = np.flatnonzero(c_of == c)
            ds = ds[ds > 0]
            if ds.size:
                of[ds] = of[labels[through[c]]]
    darea = np.bincount(labels.ravel(), w, n)
    dvec = np.zeros((n, 3))
    rr, cc = np.nonzero(labels > 0)
    np.add.at(dvec, labels[rr, cc], unit(lat[rr], lon[cc]) * np.broadcast_to(area[:, None], labels.shape)[rr, cc][:, None])
    dvec /= np.maximum(np.linalg.norm(dvec, axis=1, keepdims=True), 1e-12)
    from scipy.spatial import cKDTree

    def merge_small(of):
        """Regions under REGION_MIN_KM2 join the neighbour they share the most edge with among those the join keeps
        within REGION_MAX_KM2 (else the longest edge); an island region joins the nearest region across the water."""
        for _ in range(400):
            rmap = of[labels]
            ra = np.bincount(rmap.ravel(), w, int(of.max()) + 1)
            small = sorted((r for r in np.unique(of[1:]) if ra[r] < REGION_MIN_KM2 and r not in exempt),
                           key=lambda r: ra[r])
            if not small:
                return of
            r = small[0]
            pairs, length, _ = adjacency(rmap, np.zeros(rmap.shape, np.float32))
            cand = [((ra[r] + ra[b if a == r else a] <= REGION_MAX_KM2), L, b if a == r else a)
                    for (a, b), L in zip(pairs, length) if r in (a, b) and (b if a == r else a) not in exempt]
            if cand:
                of[of == r] = max(cand)[2]
                continue
            mine = np.argwhere(rmap == r)
            other = np.argwhere((rmap > 0) & (rmap != r) & ~np.isin(rmap, list(exempt)))
            tree = cKDTree(unit(lat[other[:, 0]], lon[other[:, 1]]))
            _, k = tree.query(unit(lat[mine[:, 0]], lon[mine[:, 1]]))
            of[of == r] = np.bincount(rmap[other[k, 0], other[k, 1]]).argmax()
        return of

    def split_big(of):
        """Regions over REGION_MAX_KM2 split by grouping their districts, kept whole, round centres seeded far apart and
        settled by area-weighted k-means."""
        nxt, count = int(of.max()) + 1, 0
        for r in np.unique(of[1:]):
            mine = np.flatnonzero(of == r)
            mine = mine[mine > 0]
            total = darea[mine].sum()
            if total <= REGION_MAX_KM2 or r in exempt:
                continue
            k = int(np.ceil(total / (0.7 * REGION_MAX_KM2)))
            v = dvec[mine]
            seeds = [int(np.argmin(v @ (v * darea[mine, None]).sum(axis=0)))]
            while len(seeds) < k:
                seeds.append(int(np.argmin(np.max(v @ v[seeds].T, axis=1))))
            centres = v[seeds]
            for _ in range(30):
                g = np.argmax(v @ centres.T, axis=1)
                centres = np.stack([(v[g == j] * darea[mine][g == j, None]).sum(axis=0) if (g == j).any()
                                    else centres[j] for j in range(k)])
                centres /= np.linalg.norm(centres, axis=1, keepdims=True)
            for j in range(1, k):
                of[mine[g == j]] = nxt
                nxt += 1
            count += 1
        return of, count

    for _ in range(8):
        of = merge_small(of)
        of, count = split_big(of)
        if not count:
            break
    uniq = np.unique(of[1:])
    remap = np.zeros(int(of.max()) + 1, np.int32)
    remap[uniq] = np.arange(1, uniq.size + 1)
    of = remap[of]
    island_regions = {int(remap[r]) for r in exempt if r < remap.size and remap[r]}
    log(f'regions on land: {uniq.size}, {len(island_regions)} of them island groups')
    return of, info, island_regions


# ---------------------------------------------------------------------------------------------- the draft

def shape_of(mask, kind, islands, la, lo, lat, lon, reach_km=120.0, rays=16):
    """Whether a coastal sea is a bay (land within reach_km along ten in sixteen rays cast from its label point) and
    whether it is a strait (it lies within COAST_KM of two separate land masses)."""
    hits = 0
    for k in range(rays):
        az = 2 * np.pi * k / rays
        for dkm in np.linspace(20.0, reach_km, 12):
            ang = dkm / MOON_KM
            la2 = np.degrees(np.arcsin(np.sin(np.radians(la)) * np.cos(ang)
                                       + np.cos(np.radians(la)) * np.sin(ang) * np.cos(az)))
            lo2 = lo + np.degrees(np.arctan2(np.sin(az) * np.sin(ang) * np.cos(np.radians(la)),
                                             np.cos(ang) - np.sin(np.radians(la)) * np.sin(np.radians(la2))))
            r = min(int((90.0 - la2) * PPD), lat.size - 1)
            c = int((lo2 % 360.0) * PPD) % lon.size
            if kind[r, c] < 2:
                hits += 1
                break
    from scipy import ndimage
    steps = int(np.ceil(COAST_KM / (MOON_KM * np.radians(1.0 / PPD)))) + 1
    rr, cc = np.nonzero(mask)
    r0, r1 = max(rr.min() - steps, 0), min(rr.max() + steps + 1, mask.shape[0])
    wide = cc.max() - cc.min() > mask.shape[1] // 2
    c0, c1 = (0, mask.shape[1]) if wide else (max(cc.min() - steps, 0), min(cc.max() + steps + 1, mask.shape[1]))
    near = ndimage.binary_dilation(mask[r0:r1, c0:c1], iterations=steps)
    masses = np.unique(islands[r0:r1, c0:c1][near & (kind[r0:r1, c0:c1] < 2)])
    return bool(hits >= rays * 10 / 16), bool((masses > 0).sum() >= 2)


def feature_units(labels, coast, lat, lon, pinned=()):
    """The one unit that carries each named mare, ocean, bay, lake and marsh of the IAU gazetteer, and each feature the
    names table pins by name: larger features first, each the unit holding its centre or the nearest free unit within
    FEATURE_REACH of its diameter (FEATURE_REACH_KM at most). Returns {('land' | 'coast', id): feature} with the
    feature's name, type, diameter and the distance from its centre to the unit's nearest cell."""
    feats = [f for f in nm.features() if f.get('type') in FEATURE_TYPES or f.get('name') in pinned]
    feats.sort(key=lambda f: (-f.get('diameter_km', 0.0), f['name']))
    out = {}
    for f in feats:
        lo0 = f['lon'] % 360.0
        reach = min(FEATURE_REACH * f.get('diameter_km', 0.0), FEATURE_REACH_KM)
        rr, cc = disc_cells(f['lat'], lo0, max(reach, 1.0), lat, lon)
        if not rr.size:
            continue
        d = np.degrees(np.arccos(np.clip(unit(lat[rr], lon[cc]) @ unit(np.array(f['lat']), np.array(lo0)), -1.0, 1.0)))
        km = np.radians(d) * MOON_KM
        best = None
        for k in np.argsort(km):
            r, c = int(rr[k]), int(cc[k])
            key = ('land', int(labels[r, c])) if labels[r, c] else ('coast', int(coast[r, c])) if coast[r, c] else None
            if key and key not in out:
                best = (key, float(km[k]))
                break
        if best:
            out[best[0]] = dict(name=f['name'], type=f['type'].split(',')[0], km=round(float(f.get('diameter_km', 0.0))),
                                distance_km=round(best[1]))
    return out


def anchors(labels, lat, lon):
    """Each label's label point: its cell farthest from its own edge (a pole of inaccessibility on the grid)."""
    from scipy import ndimage
    edge = neighbours_differ(labels) | (labels == 0)
    dist = ndimage.distance_transform_edt(~edge)
    n = int(labels.max()) + 1
    best = ndimage.maximum_position(dist, labels, np.arange(1, n))
    return {i + 1: (float(lat[r]), float(lon[c])) for i, (r, c) in enumerate(best)}


def code_order(items):
    """Codes in reading order: bands of 15 degrees from north to south, west to east within each."""
    return sorted(items, key=lambda x: (int((90.0 - x[1]) // 15), ((x[2] + 180.0) % 360.0), x[0]))


def build(log=print):
    started = time.time()
    lat, lon, area = grid()
    h, kind, sea_body, meta, b16 = surface_classes(log)
    rotation = np.asarray(meta['frame']['rotation_today_to_turned'], float)
    to8, q8, end16, to16 = routing_fields(b16, log)
    moisture, zone_names = moisture_classes(rotation, log)
    basalt = mare_map.basalt(PPD, left_lon=0.0).astype(bool)
    on_land = kind < 2
    cost = cost_surface(h, kind, moisture, basalt)
    markers, info = place_centres(on_land, kind, lat, lon, FEATURE_KM, SEED_KM, LAKE_SEED_KM2, log)
    labels = watershed(cost, markers)
    labels[~on_land] = 0
    base = labels.copy()
    from scipy import ndimage
    islands, _ = ndimage.label(on_land, structure=np.ones((3, 3)))
    islands = join_wrapped(islands)
    labels = relabel(labels, np.arange(int(labels.max()) + 1))
    log(f'watershed: {int(labels.max())} first districts, {time.time() - started:.0f} s')
    territory, comps, joined, through, groups, rocks = island_territories(on_land, area, lat, lon, log)
    labels = apply_rules(labels, h, moisture, basalt, to8, area, cost, lat, lon, islands, log)
    labels = absorb_islets(labels, on_land, area, lat, lon, territory)
    labels = relabel(labels, np.arange(int(labels.max()) + 1))
    log(f'land districts after the rules: {int(labels.max())}')
    region_of, region_info, island_regions = regions(labels, on_land, kind, cost, area, lat, lon, log,
                                                     territory=territory, comps=comps, joined=joined, through=through,
                                                     groups=groups)
    coast, faces, coast_body, coast_region = coastal_seas(labels, kind, sea_body, h, area, lat, lon, log,
                                                          region_of=region_of, rocks=rocks)

    # facts per land district
    s = stats(labels, h, moisture, basalt, to8, area)
    centre_kind = {}
    for i in range(1, int(labels.max()) + 1):
        ids, cnt = np.unique(base[labels == i], return_counts=True)
        top = ids[np.argmax(cnt)]
        centre_kind[i] = info[top] if 0 < top < len(info) else ('ground', 0.0, 0.0, 0.0, None)
    lakes, nl = ndimage.label(kind == 1, structure=np.ones((3, 3)))
    lake_counts = {}
    if nl:
        pair = np.unique(lakes.ravel().astype(np.int64) * (int(labels.max()) + 1) + labels.ravel())
        for p in pair:
            li, d = divmod(int(p), int(labels.max()) + 1)
            if li and d:
                lake_counts[d] = lake_counts.get(d, 0) + 1
    lake_area = np.bincount(labels.ravel(), (np.broadcast_to(area[:, None], labels.shape) * (kind == 1)).ravel(),
                            int(labels.max()) + 1)
    qmax = ndimage.maximum(q8, labels, np.arange(int(labels.max()) + 1)) if labels.max() else []
    p10 = ndimage.labeled_comprehension(h, labels, np.arange(int(labels.max()) + 1), lambda v: np.percentile(v, 10),
                                        float, 0.0)
    p90 = ndimage.labeled_comprehension(h, labels, np.arange(int(labels.max()) + 1), lambda v: np.percentile(v, 90),
                                        float, 0.0)
    falls = json.loads(FALLS.read_text()) if FALLS.exists() else {'reaches': []}
    cols16 = b16.shape[1]
    fall_at = {}
    for x in falls.get('reaches', []):
        r16, c16 = divmod(int(x['start']), cols16)
        d = int(labels[r16 // 2, c16 // 2])
        if d:
            fall_at.setdefault(d, []).append(float(x['power_gw']))
    cuts = DC.load(log=lambda *_: None)
    cut_at = set()
    for c in cuts['cut_cell']:
        r16, c16 = divmod(int(c), cols16)
        cut_at.add(int(labels[r16 // 2, c16 // 2]))
    heritage = json.loads(HERITAGE.read_text())['sites']
    sites_at = {}
    for site in heritage:
        r = min(int((90.0 - site['lat']) * PPD), lat.size - 1)
        c = int((site['lon'] % 360.0) * PPD) % lon.size
        d = int(labels[r, c])
        if d:
            sites_at.setdefault(d, []).append(site['name'])

    # links of the coast (the hybrid): the land districts each coastal sea faces, and those whose rivers end in it
    nc = int(coast.max())
    face_area = {}
    w = np.broadcast_to(area[:, None], labels.shape)
    m = coast > 0
    for cd, ld, a in zip(coast[m], faces[m], w[m]):
        face_area.setdefault(int(cd), {}).setdefault(int(ld), 0.0)
        face_area[int(cd)][int(ld)] += float(a)
    rows16, _ = b16.shape
    ends = end16.reshape(rows16, cols16)
    end8 = (ends // cols16) // 2, (ends % cols16) // 2
    end_coast = coast[end8]
    lab16 = np.repeat(np.repeat(labels, 2, axis=0), 2, axis=1)
    area16 = np.repeat(np.repeat(np.broadcast_to(area[:, None], labels.shape) / 4.0, 2, axis=0), 2, axis=1)
    ok = (lab16 > 0) & (end_coast > 0)
    key = end_coast[ok].astype(np.int64) * (int(labels.max()) + 1) + lab16[ok]
    uk, inv = np.unique(key, return_inverse=True)
    fed = np.bincount(inv, area16[ok])
    fed_by = {}
    for k_, a in zip(uk, fed):
        cd, ld = divmod(int(k_), int(labels.max()) + 1)
        fed_by.setdefault(cd, {})[ld] = float(a)
    mouth_q = ndimage.maximum(np.where(kind == 2, 0, q8), labels, np.arange(int(labels.max()) + 1))

    # codes: land regions in reading order, districts in reading order within each
    land_anchor = anchors(labels, lat, lon)
    region_cells = region_of[labels]
    region_anchor = anchors(region_cells, lat, lon)
    region_codes = {}
    for k, (r, la, lo) in enumerate(code_order([(r, *region_anchor[r]) for r in region_anchor])):
        region_codes[r] = f'L{k + 1:02d}'
    district_codes = {}
    for r in region_anchor:
        mine = [(d, *land_anchor[d]) for d in land_anchor if region_of[d] == r]
        for k, (d, la, lo) in enumerate(code_order(mine)):
            district_codes[d] = f'{region_codes[r]}.{k + 1}'
    bodies = meta['bodies']
    sea_codes = {b: f'S{b}' for b in range(1, len(bodies))}
    coast_anchor = anchors(coast, lat, lon)
    coast_codes = {}
    for b in sea_codes:
        mine = [(cd, *coast_anchor[cd]) for cd in coast_anchor if coast_body[cd] == b]
        for k, (cd, la, lo) in enumerate(code_order(mine)):
            coast_codes[cd] = f'{sea_codes[b]}.{k + 1}'

    from research.studies.atlas_round import district_names as DN
    carried = feature_units(labels, coast, lat, lon, pinned=DN.pinned_features())

    def place(la, lo):
        t = to_turned(la, lo, rotation)
        return dict(today=[round(la, 3), round(((lo + 180.0) % 360.0) - 180.0, 3)],
                    turned=[round(float(t[0]), 3), round(float(t[1]), 3)])

    def drains(i):
        row = s['to'][i]
        b = int(np.argmax(row)) - 1
        share = row.max() / max(row.sum(), 1e-9)
        if b >= 1:
            return sea_codes[b], round(float(share), 2)
        return ('its own lakes' if b == 0 else 'dry hollows'), round(float(share), 2)

    districts = []
    for i in sorted(district_codes, key=lambda d: district_codes[d]):
        ck = centre_kind[i]
        zone = zone_names[int(s['moist'][i].argmax())]
        zone_share = float(s['moist'][i].max() / max(s['moist'][i].sum(), 1e-9))
        relief = float(p90[i] - p10[i])
        landform = (f'ground in a {ck[1]:.0f} km crater' if ck[0] == 'crater' and ck[1] < 150 else
                    f'ground in a {ck[1]:.0f} km basin' if ck[0] == 'crater' else
                    'lake country' if ck[0] == 'lake' and lake_area[i] / s['area'][i] >= 0.05 else
                    'island' if ck[0] == 'island' else
                    'rugged highland' if relief > 2500 else 'hill country' if relief > 1200 else 'plain')
        to_code, to_share = drains(i)
        falls_here = sorted(fall_at.get(i, []), reverse=True)
        bits = [f"{zone[:1].upper() + zone[1:]} {landform}{' on mare basalt' if s['basalt'][i] >= 0.5 else ''}"]
        if lake_counts.get(i):
            share = lake_area[i] / s['area'][i]
            bits.append(f"{lake_counts[i]} lake{'s' if lake_counts[i] > 1 else ''} cover{'' if lake_counts[i] > 1 else 's'} "
                        f"{share:.0%} of it" if share >= 0.01 else
                        f"{lake_counts[i]} small lake{'s' if lake_counts[i] > 1 else ''}")
        if float(qmax[i]) >= 1000:
            bits.append(f"a river of {float(qmax[i]):,.0f} m³/s")
        if falls_here:
            bits.append(f"{len(falls_here)} steep reach{'es' if len(falls_here) > 1 else ''}, up to {falls_here[0]:.0f} GW")
        if i in cut_at:
            bits.append("the works' cut")
        if sites_at.get(i):
            bits.append('heritage: ' + '; '.join(sites_at[i]))
        drain_text = f"drains to {to_code}" if to_code.startswith('S') else f"drains to {to_code}"
        description = ', '.join(bits) + f"; {drain_text}."
        la, lo = land_anchor[i]
        coasts = sorted({coast_codes[cd] for cd, f in face_area.items() if i in f and cd in coast_codes}
                        | {coast_codes[cd] for cd, f in fed_by.items() if i in f and cd in coast_codes})
        districts.append(dict(code=district_codes[i], region=region_codes[region_of[i]], kind='land',
                              area_km2=round(float(s['area'][i])), at=place(la, lo), zone=zone,
                              zone_share=round(zone_share, 2), basalt_share=round(float(s['basalt'][i]), 2),
                              centre=ck[0], centre_km=round(ck[1]) if ck[0] == 'crater' else None,
                              centre_name=ck[4] if ck[0] == 'crater' else None,
                              centre_at=[round(ck[2], 2), round(((ck[3] + 180.0) % 360.0) - 180.0, 2)],
                              mean_height_m=round(float(s['h'][i]) - float(meta['common_level_m'])),
                              relief_m=round(relief), lakes=int(lake_counts.get(i, 0)),
                              lake_share=round(float(lake_area[i] / s['area'][i]), 3),
                              largest_flow_m3s=round(float(qmax[i])), drains_to=to_code, drains_share=to_share,
                              falls_gw=[round(g, 1) for g in falls_here[:3]], cut=i in cut_at,
                              heritage=sites_at.get(i, []), coasts=coasts, description=description,
                              feature=carried.get(('land', i))))
        districts[-1]['anchor'] = DN.anchor_of(districts[-1])
    sea_names = {b: bodies[b] for b in sea_codes}
    # bays and straits are judged against land masses of STRAIT_MASS_KM2 or more, so that islets make neither
    mass_area = np.bincount(islands.ravel(), np.broadcast_to(area[:, None], labels.shape).ravel())
    big_masses = np.flatnonzero(mass_area >= STRAIT_MASS_KM2)
    big_masses = big_masses[big_masses > 0]
    coastal = []
    depth = np.zeros(labels.shape, np.float32)
    lev = {}
    with np.load(GS.PRODUCT, allow_pickle=False) as z:
        wc, wl = z['water_cell'], z['water_level_m']
    lev16 = np.full(b16.size, np.nan, np.float32)
    lev16[wc] = wl
    lev8 = np.nanmean(lev16.reshape(b16.shape[0] // 2, 2, b16.shape[1] // 2, 2), axis=(1, 3))
    depth = np.where(kind == 2, lev8 - h, 0.0)
    for cd in sorted(coast_codes, key=lambda c: coast_codes[c]):
        mcd = coast == cd
        a = float(w[mcd].sum())
        dmean = float((depth[mcd] * w[mcd]).sum() / a)
        b = int(coast_body[cd])
        f_ = face_area.get(cd, {})
        facing = [district_codes[d] for d, _ in sorted(f_.items(), key=lambda x: -x[1]) if d in district_codes]
        feed = fed_by.get(cd, {})
        feeding = [district_codes[d] for d, x in sorted(feed.items(), key=lambda x: -x[1])
                   if d in district_codes and x >= 0.05 * float(s['area'][d])]
        inflow = max([float(mouth_q[d]) for d in feed if d < len(mouth_q)] or [0.0])
        la, lo = coast_anchor[cd]
        salt = sea_names[b] in SALT_SEAS
        shallow = 'shallow' if dmean < 200 else 'deep' if dmean > 1500 else 'shelving'
        enclosed, strait = shape_of(coast == cd, np.where(np.isin(islands, big_masses), kind, 2), islands, la, lo, lat,
                                    lon)
        form = 'strait' if strait else 'bay' if enclosed else 'coastal sea'
        description = (f"{'Salt' if salt else 'Fresh'} {form}, {shallow}, {dmean:,.0f} m deep on average; faces "
                       f"{', '.join(facing[:4]) or 'no district'}"
                       + (f"; rivers reach it from {', '.join(feeding[:4])}" if feeding else '') + '.')
        coastal.append(dict(code=coast_codes[cd], region=region_codes[int(coast_region[cd])], sea=sea_codes[b],
                            kind='coastal sea', area_km2=round(a),
                            at=place(la, lo), mean_depth_m=round(dmean), salt=salt, enclosed=enclosed, strait=strait,
                            faces=facing, fed_by=feeding, description=description, feature=carried.get(('coast', cd))))
        coastal[-1]['anchor'] = DN.anchor_of(coastal[-1])
    regions_out = []
    for r, code in sorted(region_codes.items(), key=lambda x: x[1]):
        mine = [d for d in districts if d['region'] == code]
        a = sum(d['area_km2'] for d in mine)
        zones = {}
        for d in mine:
            zones[d['zone']] = zones.get(d['zone'], 0) + d['area_km2']
        top = max(zones, key=zones.get)
        basaltic = sum(d['area_km2'] * d['basalt_share'] for d in mine) / max(a, 1)
        relief_mean = sum(d['area_km2'] * d['relief_m'] for d in mine) / max(a, 1)
        coastwise = sum(d['area_km2'] for d in mine if d['coasts']) / max(a, 1)
        shores = {x for c in coastal for x in c['faces']}
        shorewise = sum(d['area_km2'] for d in mine if d['code'] in shores) / max(a, 1)
        lakewise = sum(d['area_km2'] * d['lake_share'] for d in mine) / max(a, 1)
        la, lo = region_anchor[r]
        waters = sum(c['area_km2'] for c in coastal if c['region'] == code)
        regions_out.append(dict(code=code, kind='land', area_km2=a, districts=len(mine), at=place(la, lo),
                                island_group=r in island_regions, waters_km2=waters,
                                zones={k: round(v / max(a, 1), 3) for k, v in zones.items()},
                                relief_m=round(relief_mean), basalt_share=round(basaltic, 3),
                                coast_share=round(coastwise, 3), shore_share=round(shorewise, 3),
                                lake_share=round(lakewise, 3),
                                description=(f"An island group in international waters: " if r in island_regions else '')
                                            + f"{len(mine)} districts over {a:,} km², mostly {top}"
                                            + (', much of it mare basalt' if basaltic > 0.3 else '')
                                            + (f"; its waters {waters:,} km²" if waters else '') + '.'))
    sea_anchor = anchors(np.where(kind == 2, sea_body, 0).astype(np.int32), lat, lon)
    for b, code in sea_codes.items():
        mine = [c for c in coastal if c['sea'] == code]
        a = float(w[(kind == 2) & (sea_body == b)].sum())
        la, lo = sea_anchor.get(b, (0.0, 0.0))
        salt = sea_names[b] in SALT_SEAS
        regions_out.append(dict(code=code, kind='sea', area_km2=round(a), districts=len(mine), at=place(la, lo),
                                description=f"{'A salt' if salt else 'A fresh'} sea of {a:,.0f} km²: {len(mine)} coastal "
                                            f"seas, the waters of the regions round it, and beyond them international "
                                            f"waters, which belong to no region."))

    # names in the Moon's tongues (district_names), each region in its own; descriptions then read names for codes
    body_of_code = {code: bodies[b] for b, code in sea_codes.items()}
    for r in regions_out:
        if r['kind'] == 'sea':
            r['body'] = body_of_code[r['code']]
    units = dict(regions=regions_out, districts=districts, coastal_seas=coastal)
    DN.anchor_all(units)                                     # every anchor unique
    for x in regions_out + districts + coastal:
        x['description_codes'] = x['description']
    previous = DN.previous_names(json.loads(OUT.read_text())) if OUT.exists() else None
    DN.apply(units, previous=previous)                       # a redraw keeps the names whose anchors survive

    # the product: district and region at each cell, on today's grid
    district_grid = np.zeros(labels.shape, np.int16)
    region_grid = np.zeros(labels.shape, np.int16)
    index = {}
    for k, d in enumerate(districts + coastal):
        index[d['code']] = k + 1
    for i, code in district_codes.items():
        district_grid[labels == i] = index[code]
    for cd, code in coast_codes.items():
        district_grid[coast == cd] = index[code]
    for r, code in region_codes.items():
        region_grid[(labels > 0) & (region_of[labels] == r)] = int(code[1:])
    for b, code in sea_codes.items():
        region_grid[(kind == 2) & (sea_body == b) & (coast == 0)] = 100 + b      # international waters
    for cd, code in coast_codes.items():
        region_grid[coast == cd] = int(region_codes[int(coast_region[cd])][1:])   # a region's waters
    tracked = {rel(p): sha256(p) for p in (GS.OUT, DR.OUT, DC.OUT, MF.OUT, FALLS, HERITAGE)}
    ignored = {str(p): sha256(p) for p in (GS.PRODUCT, DR.PRODUCT, DC.PRODUCT)}
    pmeta = dict(schema=PRODUCT_SCHEMA, grid=dict(pixels_per_degree=PPD, rows=labels.shape[0], cols=labels.shape[1]),
                 reading_rule=('Rows from 90 N and columns from 0 E, 8 per degree, on today\'s Moon. district holds 1 + '
                               'the index of each cell\'s district in the record\'s list (land districts, then coastal seas), '
                               '0 on open water; region holds the region\'s number (L01 is 1) on its land and in its '
                               'waters, the coastal seas, and 100 + the sea\'s body in international waters.'),
                 codes=[d['code'] for d in districts + coastal])
    BULK.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PRODUCT, metadata=json.dumps(pmeta), district=district_grid, region=region_grid)
    sizes = np.array([d['area_km2'] for d in districts])
    summary = dict(land_regions=len(region_codes), island_group_regions=len(island_regions), sea_regions=len(sea_codes),
                   land_districts=len(districts),
                   coastal_sea_districts=len(coastal),
                   land_district_km2=dict(median=int(np.median(sizes)), p10=int(np.percentile(sizes, 10)),
                                          p90=int(np.percentile(sizes, 90)), smallest=int(sizes.min()),
                                          largest=int(sizes.max())),
                   under_min=int((sizes < MIN_KM2).sum()), over_max=int((sizes > MAX_KM2).sum()),
                   coastal_km2_median=int(np.median([c['area_km2'] for c in coastal])) if coastal else 0)
    log(f'districts: {json.dumps(summary)}, {time.time() - started:.0f} s')
    return dict(schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}), inputs=tracked,
                ignored_inputs=ignored,
                evidence=('A first draft by the rules of districts.md on the design surface, design climate and design '
                          'routing: computed units and facts with short computed descriptions; their significance '
                          'awaits written judgment and the author\'s review. Codes stand in for names.'),
                reading_rule=('Areas in km2, heights in metres above the common level; at.today and at.turned are each '
                              'unit\'s label point (latitude, longitude east). coasts lists the coastal sea districts a '
                              'land district faces or feeds; faces and fed_by list the land districts a coastal sea '
                              'faces and those whose rivers reach it with at least 5% of their area.'),
                parameters=dict(ppd=PPD, min_km2=MIN_KM2, max_km2=MAX_KM2, islet_km2=ISLET_KM2, feature_km=FEATURE_KM,
                                region_feature_km=REGION_FEATURE_KM, seed_km=SEED_KM, region_seed_km=REGION_SEED_KM,
                                lake_seed_km2=LAKE_SEED_KM2, moisture_edge_m=MOISTURE_EDGE_M, rock_edge_m=ROCK_EDGE_M,
                                weak_edge_m=WEAK_EDGE_M, coast_km=COAST_KM, coast_max_km2=COAST_MAX_KM2,
                                region_min_km2=REGION_MIN_KM2, region_max_km2=REGION_MAX_KM2, smooth_deg=SMOOTH_DEG),
                tongues={k: dict(label=t['label'], roots=t['roots'], order=t['order'], where=DN.WHERE[k])
                         for k, t in DN.TONGUES.items()},
                terlunan_coined=DN.TONGUES['terlunan']['coined'],
                summary=summary, regions=regions_out, districts=districts, coastal_seas=coastal,
                product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)), seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    record = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    print(f"wrote {OUT}: {json.dumps(record['summary'])}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
