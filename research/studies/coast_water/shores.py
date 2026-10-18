"""The shores at LOLA's native 16 px/deg: every water body and every land pixel near one, for the coast-water study.

    python -m research.studies.coast_water.shores      # about twenty seconds, 1.5 GB

Heights are LOLA at 16 px/deg above the degree-200 GRAIL geoid, as geography.drainage takes them
(geography.topography.height_above_geoid), against the 28% atlas's sea level. The water is the drainage product's:
every pixel below the sea level is common-level water, as the atlas fills it, and the rain-fed lakes above it are the
product's lake cells with their levels and whether each spills. Common-level bodies of at least the atlas's listing
size (geography.atlas.BODY_MIN_SHARE of the Moon) are seas and the smaller ones common-level lakes. Land bodies are
the 8-connected land with the rain-fed lakes in it, the largest the main land and the rest islands; labels join across
the longitude seam and across each pole.

For each land pixel within NEAR_KM of water it keeps the distance to the nearest water pixel's centre, the step to
that pixel in rows and columns, its height above that water's level, the water's class, the rain-fed lake that water
belongs to and the island it lies on, if any. For each land pixel with water beside it, it keeps the length of its
edges facing each class of water, its least height above the facing water of that class and its depth across those
edges. Distances come from Euclidean distance transforms in bands of one degree of latitude, each with the east-west
spacing at its centre, tiled across the seam. It also keeps each T21 cell's land and water at this resolution and the
rain-fed lakes' areas. The product goes to research/runs/coast_water/shores_16ppd.npz (ignored), with its inputs'
sha256s in its metadata.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from geography import atlas as GA
from geography import topography as tp
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / 'research' / 'runs' / 'coast_water' / 'shores_16ppd.npz'
SCHEMA = 'terluna.research.coast-water-shores/1'
CODE = ('research/studies/coast_water/shores.py', 'geography/topography.py', 'geography/atlas.py')
PPD = 16
ROWS, COLS = 180 * PPD, 360 * PPD
ATLAS = 'geography/results/atlas.json'
DRAINAGE = 'geography/products/drainage_28pct_16ppd.npz'
LOLA = ('geography/inputs/ldem_16.img', 'geography/inputs/ldem_4.img', 'geography/inputs/jggrx_0420a_sha.tab')
NEAR_KM = 12.0                   # land pixels kept: within this distance of a water pixel's centre
BAND_ROWS, MARGIN_ROWS = 16, 24  # distance transforms over a degree of latitude with a degree and a half either side
SEA, LAKE_COMMON, LAKE_OPEN, LAKE_CLOSED = 1, 2, 3, 4
WATER_CLASSES = {SEA: 'sea', LAKE_COMMON: 'common_level_lake', LAKE_OPEN: 'open_lake', LAKE_CLOSED: 'closed_lake'}


def digest(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def grid():
    """Each row's latitude, the north-south spacing (km), each row's east-west spacing at its centre and at its
    northern edge (and the southern pole's), and each row's pixel area (km2)."""
    lat = 90.0 - (np.arange(ROWS) + 0.5) / PPD
    dy = 2.0 * np.pi * MOON_RADIUS / 1000.0 / 360.0 / PPD
    dx = dy * np.cos(np.radians(lat))
    edges = 90.0 - np.arange(ROWS + 1) / PPD
    dx_edge = dy * np.abs(np.cos(np.radians(edges)))
    r = MOON_RADIUS / 1000.0
    area = r * r * np.radians(1.0 / PPD) * (np.sin(np.radians(edges[:-1])) - np.sin(np.radians(edges[1:])))
    return lat, dy, dx, dx_edge, area


def merge_labels(labels):
    """Labels joined across the longitude seam (8-connected) and across each pole, where a polar row's pixel touches
    the one half a turn round; renumbered 1..n by first appearance in sorted order, 0 kept for the background."""
    n, cols = int(labels.max()), labels.shape[1]
    pairs = []
    west, east = labels[:, 0], labels[:, -1]
    for shift in (-1, 0, 1):
        e = np.roll(east, shift)
        if shift == 1:
            e[0] = 0
        elif shift == -1:
            e[-1] = 0
        pairs.append(np.stack([west, e], 1))
    for row in (labels[0], labels[-1]):
        pairs.append(np.stack([row, np.roll(row, cols // 2)], 1))
    p = np.concatenate(pairs).astype(np.int64)
    p = p[(p[:, 0] > 0) & (p[:, 1] > 0)]
    graph = coo_matrix((np.ones(len(p)), (p[:, 0], p[:, 1])), shape=(n + 1, n + 1))
    _, comp = connected_components(graph, directed=False)
    names = np.unique(comp[1:])
    renumber = np.zeros(comp.max() + 1, np.int32)
    renumber[names] = np.arange(1, names.size + 1, dtype=np.int32)
    table = renumber[comp]
    table[0] = 0
    return table[labels]


def body_areas(labels, row_area):
    """Area (km2) of each label, row by row."""
    out = np.zeros(int(labels.max()) + 1)
    for r in range(labels.shape[0]):
        out += np.bincount(labels[r], minlength=out.size) * row_area[r]
    return out


def distances(water, lat, dy, band_rows=BAND_ROWS, margin_rows=MARGIN_ROWS):
    """Distance (km) from each pixel's centre to the nearest water pixel's centre and that pixel's flat index, by
    bands of latitude tiled three times round the seam; each band takes the east-west spacing at its centre."""
    rows, cols = water.shape
    dist = np.full((rows, cols), np.inf, np.float32)
    near = np.full((rows, cols), -1, np.int32)
    for r0 in range(0, rows, band_rows):
        r1 = min(rows, r0 + band_rows)
        a, b = max(0, r0 - margin_rows), min(rows, r1 + margin_rows)
        sub = water[a:b]
        if not sub.any():
            continue
        tiled = np.concatenate([sub, sub, sub], axis=1)
        sampling = (dy, dy * max(float(np.cos(np.radians(lat[r0:r1].mean()))), 1e-3))
        d, (ii, jj) = ndimage.distance_transform_edt(~tiled, sampling=sampling, return_indices=True)
        keep, mid = slice(r0 - a, r1 - a), slice(cols, 2 * cols)
        dist[r0:r1] = d[keep, mid]
        near[r0:r1] = ((ii[keep, mid] + a) * cols + (jj[keep, mid] % cols)).astype(np.int32)
    return dist, near


def shore_edges(keep, wc, lv, lake_body, h, dy, dx, dx_edge):
    """The land pixels in keep with water beside them, and by class of that water (sea, common-level lake, open and
    closed rain-fed lake, in columns) the length of their edges facing it (km), their depth across those edges (km,
    length-weighted), their least height above that water's level (m), and the rain-fed lake they face. East-west
    neighbours share a north-south edge of length dy and the pixel spans its row's dx across it; north-south
    neighbours share the row edge's east-west length dx_edge and the pixel spans dy across it."""
    rows, cols = keep.shape
    coastal = np.zeros((rows, cols), bool)
    for dr, dc in ((0, 1), (0, -1), (-1, 0), (1, 0)):
        coastal |= neighbour(wc, dr, dc, 0) > 0
    coastal &= keep
    cidx = np.flatnonzero(coastal)
    position = np.full(rows * cols, -1, np.int64)
    position[cidx] = np.arange(cidx.size)
    edge = np.zeros((cidx.size, 4), np.float32)
    across = np.zeros((cidx.size, 4), np.float32)
    rise = np.full((cidx.size, 4), np.nan, np.float32)
    facing_lake = np.zeros(cidx.size, np.int32)
    for dr, dc in ((0, 1), (0, -1), (-1, 0), (1, 0)):
        nb_class, nb_level = neighbour(wc, dr, dc, 0), neighbour(lv, dr, dc, np.nan)
        nb_lake = neighbour(lake_body, dr, dc, 0)
        rr, cc = np.nonzero(coastal & (nb_class > 0))
        p = position[rr * cols + cc]
        k = nb_class[rr, cc].astype(np.int64) - 1
        np.maximum.at(facing_lake, p, nb_lake[rr, cc])
        if dr == 0:
            length, depth = np.full(rr.size, dy), dx[rr]
        else:
            length, depth = (dx_edge[rr] if dr == -1 else dx_edge[rr + 1]), np.full(rr.size, dy)
        np.add.at(edge, (p, k), length.astype(np.float32))
        np.add.at(across, (p, k), (length * depth).astype(np.float32))
        z = (h[rr, cc] - nb_level[rr, cc]).astype(np.float32)
        np.fmin.at(rise, (p, k), z)
    across = np.where(edge > 0, across / np.maximum(edge, 1e-12), 0.0).astype(np.float32)
    return cidx, edge, across, rise, facing_lake


def t21_index(t21_lat):
    """Each 16 px/deg row's and column's T21 cell: the Gaussian rows bounded where the cumulative Gaussian weights put
    them (rows in the output's order, from the north), the columns centred on the output's longitudes."""
    _, weights = np.polynomial.legendre.leggauss(t21_lat.size)
    edges = np.concatenate([[-1.0], -1.0 + np.cumsum(weights)])               # sin(latitude), from the south
    lat = 90.0 - (np.arange(ROWS) + 0.5) / PPD
    k = np.clip(np.searchsorted(edges, np.sin(np.radians(lat))) - 1, 0, t21_lat.size - 1)
    south_first = np.sort(t21_lat)
    row = np.array([int(np.argmin(np.abs(t21_lat - south_first[i]))) for i in k], np.int16)
    lon = (np.arange(COLS) + 0.5) / PPD
    col = (np.floor((lon + 180.0 / 64.0) / (360.0 / 64.0)).astype(int) % 64).astype(np.int16)
    return row, col


def build():
    atlas = json.loads((ROOT / ATLAS).read_text())
    level = float(atlas['sea_level_m'])
    lat, dy, dx, dx_edge, row_area = grid()
    h, _, _ = tp.height_above_geoid(PPD, atlas['grid']['geoid_degree'])
    h = h.astype(np.float32)
    common = h < level
    with np.load(ROOT / DRAINAGE) as z:
        cells, lake_level, spills = z['lake_cell'].astype(np.int64), z['lake_level_m'], z['lake_spills']
    if common.ravel()[cells].any():
        raise ValueError('a rain-fed lake cell lies below the sea level')
    lake = np.zeros(ROWS * COLS, bool)
    lake[cells] = True
    lake = lake.reshape(ROWS, COLS)
    water = common | lake
    level_of = np.full(ROWS * COLS, np.nan, np.float32)
    level_of[common.ravel()] = level
    level_of[cells] = lake_level
    eight = np.ones((3, 3), bool)
    # common-level bodies: the seas and the common-level lakes
    body = merge_labels(ndimage.label(common, structure=eight)[0])
    area = body_areas(body, row_area)
    moon = 4.0 * np.pi * (MOON_RADIUS / 1000.0) ** 2
    sea_body = area >= GA.BODY_MIN_SHARE * moon
    sea_body[0] = False
    wclass = np.zeros(ROWS * COLS, np.int8)
    flat_common = common.ravel()
    wclass[flat_common] = np.where(sea_body[body.ravel()[flat_common]], SEA, LAKE_COMMON)
    wclass[cells] = np.where(spills, LAKE_OPEN, LAKE_CLOSED)
    seas = dict(count=int(sea_body.sum()), area_km2=float(area[sea_body].sum()),
                common_level_lakes=int((~sea_body[1:]).sum()),
                common_level_lake_km2=float(area[1:][~sea_body[1:]].sum()))
    del body
    # the rain-fed lakes, body by body
    lake_body = merge_labels(ndimage.label(lake, structure=eight)[0])
    lake_area = body_areas(lake_body, row_area)
    flat_lake_body = lake_body.ravel()
    closed_cells = np.bincount(flat_lake_body[cells], weights=(~spills).astype(float), minlength=lake_area.size)
    lake_closed = closed_cells > 0                          # as geography.drainage counts a body without an outlet
    rows_of = cells // COLS
    lake_row = np.bincount(flat_lake_body[cells], weights=rows_of.astype(float), minlength=lake_area.size)
    lake_n = np.bincount(flat_lake_body[cells], minlength=lake_area.size)
    lake_lat = 90.0 - (lake_row / np.maximum(lake_n, 1) + 0.5) / PPD
    # land bodies: the land with the rain-fed lakes in it
    land_body = merge_labels(ndimage.label(~common, structure=eight)[0])
    land_area = body_areas(land_body, row_area)
    land_area[0] = 0.0
    main = int(np.argmax(land_area))
    island = (land_body != main) & (land_body > 0)
    island_body = np.where(island, land_body, 0).astype(np.int32)
    del land_body
    island_land = body_areas(np.where(lake, 0, island_body), row_area)       # each island's land, its lakes left out
    island_land[0] = 0.0                                                      # everything off the islands
    island_km2 = float(island_land.sum())
    rows_index = np.broadcast_to(np.arange(ROWS)[:, None], (ROWS, COLS))
    weights = np.broadcast_to(row_area[:, None], (ROWS, COLS))
    angle = 2.0 * np.pi * (np.arange(COLS) + 0.5) / COLS
    flat_body = island_body.ravel()
    sel = flat_body > 0
    w_sel = weights.ravel()[sel]
    n_bodies = island_land.size
    mean_row = np.bincount(flat_body[sel], weights=w_sel * rows_index.ravel()[sel], minlength=n_bodies)
    total_w = np.bincount(flat_body[sel], weights=w_sel, minlength=n_bodies)
    cos_sum = np.bincount(flat_body[sel], weights=w_sel * np.broadcast_to(np.cos(angle), (ROWS, COLS)).ravel()[sel],
                          minlength=n_bodies)
    sin_sum = np.bincount(flat_body[sel], weights=w_sel * np.broadcast_to(np.sin(angle), (ROWS, COLS)).ravel()[sel],
                          minlength=n_bodies)
    island_row = mean_row / np.maximum(total_w, 1e-30)
    island_lat = 90.0 - (island_row + 0.5) / PPD
    island_lon = (np.degrees(np.arctan2(sin_sum, cos_sum)) % 360.0)
    # each island's ground above the sea: its highest pixel and its land's mean, the lakes on it left out
    on_land = sel & ~lake.ravel()
    above = h.ravel()[on_land] - level
    island_top = np.zeros(n_bodies, np.float32)
    np.maximum.at(island_top, flat_body[on_land], above.astype(np.float32))
    w_land = weights.ravel()[on_land]
    island_mean = (np.bincount(flat_body[on_land], weights=w_land * above, minlength=n_bodies)
                   / np.maximum(np.bincount(flat_body[on_land], weights=w_land, minlength=n_bodies), 1e-30))
    # distances from each land pixel to the nearest water
    dist, near = distances(water, lat, dy)
    keep = ~water & (dist <= NEAR_KM)
    idx = np.flatnonzero(keep)
    nearest = near.ravel()[idx]
    del near
    r_land, c_land = np.divmod(idx, COLS)
    r_w, c_w = np.divmod(nearest.astype(np.int64), COLS)
    step_c = (c_w - c_land) % COLS
    step_c = np.where(step_c > COLS // 2, step_c - COLS, step_c)
    pixels = dict(index=idx.astype(np.int32), distance_km=dist.ravel()[idx], step_rows=(r_w - r_land).astype(np.int16),
                  step_cols=step_c.astype(np.int16), height_m=(h.ravel()[idx] - level_of[nearest]).astype(np.float32),
                  water_class=wclass[nearest], lake=flat_lake_body[nearest].astype(np.int32),
                  island=island.ravel()[idx], island_body=flat_body[idx])
    del dist
    # edges facing water
    cidx, edge, across, rise, facing_lake = shore_edges(keep, wclass.reshape(ROWS, COLS),
                                                        level_of.reshape(ROWS, COLS), lake_body, h, dy, dx, dx_edge)
    coast = dict(index=cidx.astype(np.int32), edge_km=edge, across_km=across, rise_m=rise,
                 island=island.ravel()[cidx], lake=facing_lake, island_body=flat_body[cidx])
    return dict(level=level, lat=lat, dy=dy, dx=dx, row_area=row_area, wclass=wclass, pixels=pixels, coast=coast,
                seas=seas, island_km2=island_km2, lake_area=lake_area, lake_closed=lake_closed, lake_lat=lake_lat,
                island_land=(island & ~lake).ravel(), island_area=island_land, island_lat=island_lat,
                island_lon=island_lon, island_top=island_top, island_mean=island_mean.astype(np.float32))


def neighbour(a, dr, dc, fill):
    """Each pixel's neighbour dr rows down and dc columns east, across the seam; beyond a pole, fill."""
    out = np.roll(a, -dc, axis=1) if dc else a
    if dr == -1:
        out = np.vstack([np.full((1, a.shape[1]), fill, a.dtype), out[:-1]])
    elif dr == 1:
        out = np.vstack([out[1:], np.full((1, a.shape[1]), fill, a.dtype)])
    return out


def t21_areas(row, col, wclass, row_area, island_land):
    """Land, the islands' land and each class of water (km2) in every T21 cell at 16 px/deg."""
    out = {name: np.zeros(2048) for name in ['land', 'island_land'] + list(WATER_CLASSES.values())}
    cell_cols = col.astype(np.int64)
    w = wclass.reshape(ROWS, COLS)
    isl = island_land.reshape(ROWS, COLS)
    for r in range(ROWS):
        cell = int(row[r]) * 64 + cell_cols
        out['land'] += np.bincount(cell, weights=(w[r] == 0) * row_area[r], minlength=2048)
        out['island_land'] += np.bincount(cell, weights=isl[r] * row_area[r], minlength=2048)
        for k, name in WATER_CLASSES.items():
            out[name] += np.bincount(cell, weights=(w[r] == k) * row_area[r], minlength=2048)
    return {k: v.reshape(32, 64) for k, v in out.items()}


def main(argv=None) -> int:
    from research.studies.coast_water import extract as EX
    with np.load(EX.OUT) as z:
        t21_lat = z['lat'].copy()
    b = build()
    row, col = t21_index(t21_lat)
    cells = t21_areas(row, col, b['wclass'], b['row_area'], b['island_land'])
    meta = dict(schema=SCHEMA, sea_level_m=b['level'], near_km=NEAR_KM, pixels_per_degree=PPD,
                seas=b['seas'], island_km2=b['island_km2'], water_classes={str(k): v for k, v in WATER_CLASSES.items()},
                producer=dict(files={p: digest(ROOT / p) for p in CODE}),
                inputs={p: digest(ROOT / p) for p in (ATLAS, DRAINAGE) + LOLA},
                reading_rule=('pixel_* arrays run over the land pixels within near_km of water, index = row * 5760 + '
                              'column from 90 N and 0 E; distance_km is to the nearest water pixel\'s centre and '
                              'step_rows, step_cols the step to it; height_m is above that water\'s level; '
                              'water_class 1 sea, 2 common-level lake, 3 open rain-fed lake, 4 closed; lake is the '
                              'rain-fed lake body (0 for common-level water). coast_* arrays run over the land pixels '
                              'with water beside them, by class in columns (sea, common-level lake, open, closed): '
                              'edge_km the length of the edges facing it, across_km the pixel\'s depth across them, '
                              'rise_m the least height above the facing water of that class. cell_* are km2 in each '
                              'T21 cell (rows from the north), t21_row and t21_col each 16 px/deg row\'s and '
                              'column\'s cell. lake_* run over the rain-fed lake bodies.'))
    arrays = {f'pixel_{k}': v for k, v in b['pixels'].items()}
    arrays.update({f'coast_{k}': v for k, v in b['coast'].items()})
    arrays.update({f'cell_{k}': v for k, v in cells.items()})
    arrays.update(t21_row=row, t21_col=col, row_lat=b['lat'], row_dx_km=b['dx'], row_area_km2=b['row_area'],
                  lake_area_km2=b['lake_area'], lake_closed=b['lake_closed'], lake_lat=b['lake_lat'],
                  island_area_km2=b['island_area'], island_lat=b['island_lat'], island_lon=b['island_lon'],
                  island_top_m=b['island_top'], island_mean_m=b['island_mean'],
                  dy_km=np.array(b['dy']))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, metadata=json.dumps(meta), **arrays)
    print(f"{b['pixels']['index'].size:,} land pixels within {NEAR_KM:g} km of water; {b['coast']['index'].size:,} "
          f"beside it; seas {b['seas']}; islands {b['island_km2']:,.0f} km2")
    return 0


if __name__ == '__main__':
    sys.exit(main())
