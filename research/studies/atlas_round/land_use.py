"""A first sorting of the design surface's land into wild land, hybrid land and the land open to habitation.

    ulimit -v 8000000; python3 -m research.studies.atlas_round.land_use     # about 80 seconds on two processors

The register fixes what the three kinds of land are and what share of the land is kept wild and hybrid, and it lays out
one place, Jura Island; it holds no map of the kinds (land_use.md, part A). This module applies a short set of stated
rules to the design surface at 16 px/deg and counts what they give, for the Moon, for each moisture class, night regime
and band of turned latitude, and for each region and district of the redrawn atlas. Every rule is either a register
row, a study's application of one, or a proposal, and RULES says which. The variants differ only in what the register
leaves open: the wild share (half, or 30% beside the hybrid fifth) and where the satoyama-like mosaics belong.

The rules, in the order they are applied (each land cell takes the first kind that claims it):

1. Fixed places. The land within 32 km of the summit port is the metropolis, habitation. Jura Island takes the
   register's layout: its high arc, its outer coasts and its islets wild, the rest hybrid.
2. The satoyama-like candidate valleys are hybrid land (or, in the variant that places the managed biomes with
   habitation, the first of the rain belt's farm country).
3. Gardens fill the hybrid share from the dry classes, the land nearest a sea or great lake first.
4. In every moisture class, and within it in every night regime, the wild takes its share from the land ranked last
   for dense living; the rest is open to habitation and its farms.
5. Of the open land, the distribution C' at 8 billion people with the food blend fills the best-ranked land of each
   of its bands up to that band's need; the rest of the open land is not needed at that population.

Nothing here is a design: the rules are a first screen for the author's decisions, on products that are themselves
estimates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DRIVE = Path('/media/projectspace/terluna-research')
SURFACE = DRIVE / 'atlas-round' / 'globe_surface' / 'globe_surface_16ppd.npz'
SURFACE_RECORD = HERE / 'results' / 'globe_surface.json'
ROUTING = DRIVE / 'atlas-round' / 'design_routing' / 'design_routing_16ppd.npz'
MOISTURE = HERE / 'results' / 'moisture_map_final.json'
DISTRICTS = DRIVE / 'atlas-round' / 'corrected_rule' / 'districts_redrawn' / 'districts.json'
DISTRICTS_GRID = DRIVE / 'atlas-round' / 'corrected_rule' / 'districts_redrawn' / 'districts_8ppd.npz'
SATOYAMA = HERE / 'results' / 'satoyama_final.json'
SATOYAMA_UNITS = DRIVE / 'atlas-round' / 'satoyama_final' / 'satoyama_units.npz'
HERITAGE = ROOT / 'research' / 'studies' / 'conservation' / 'results' / 'heritage_register.json'
OUT = HERE / 'results' / 'land_use.json'
BULK = DRIVE / 'session-tools' / 'agents' / 'atlas_land_use'
PRODUCT = BULK / 'land_use_16ppd.npz'
SCHEMA = 'terluna.research.atlas-round-land-use/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-land-use-grid/1'
CODE = ('research/studies/atlas_round/land_use.py',)

WILD, HYBRID, OPEN = 1, 2, 3
KINDS = {WILD: 'wild', HYBRID: 'hybrid', OPEN: 'open'}
CLASS_KEYS = ('hyper_arid', 'arid', 'semi_arid', 'dry_subhumid', 'humid', 'wet')
DRY_MAX = 3                          # the drylands: hyper-arid, arid, semi-arid and dry sub-humid (the moisture map's classes)

# ---- the register's figures
HYBRID_SHARE = 0.20                  # register, 11 October: another 20% hybrid
WILD_SHARES = (0.30, 0.50)           # register, 11 October: at least half preferred, or 30% beside the hybrid fifth
SUMMIT = (5.4, -158.6)               # register, 27 September: the far side's summit, latitude and east longitude
METROPOLIS_KM = 32.0                 # register, 27 September: the metropolis lies within about 32 km of the tower

# ---- study figures the rules lean on
HIGH_M = 4000.0                      # land cities: the high ground stands 4 km or more above the sea
COMFORT_LAT = 30.0                   # land cities and the register's research practice: most comfort poleward of 30 deg
POLEWARD_FIRST_LAT = 45.0            # distributions: the dry highland towns go poleward of 45 deg first
HIGH_LATITUDE = 48.6                 # distributions, C': the high-latitude lowland homes lie poleward of 48.6 deg
GREAT_LAKE_KM2 = 1000.0              # ways of living: the great lakes are those of 1,000 km2 and more
SHORE_BELT_KM = 5.0                  # the belt beside a sea or great lake whose kinds each variant reports
SKY_BOAT_KM, SKY_FERRY_KM = 66.7, 180.0   # distributions: a four-seat sky boat's trip and the sky ferry's route

# C' at 8 billion people, 4 billion of them on the ground, fed by the blend (distributions and food systems, on the
# land of October, 23.4 million km2). Low and high ends in km2.
NEED = dict(
    metropolis=(2464.0, 2464.0),
    tower_districts=(8950.0, 16900.0),
    highland_towns=(115000.0, 193000.0),
    scattered_communities=(100000.0, 571000.0),
    fields={'30': (680000.0, 1860000.0), '50': (680000.0, 2170000.0)},   # by the wild share: kept 30% beside the parks, or half
    heritage_pasture=(298000.0, 595000.0),                             # 1 and 2 kg of heritage meat a person a year
    settled=(264000.0, 736000.0),
    receivers_share_of_settled=(0.36, 0.59))
SCATTERED_RAIN_BELT, SCATTERED_MIDDLE = 0.475, 0.252                    # shares of the scattered communities' land by band
OCTOBER_LAND_KM2 = 23.4e6
OCTOBER_FIFTH_KM2 = 4.68e6
LIT_SKY_SHARE_OF_OCTOBER_LAND = (0.00493, 0.172)   # distributions: C' as placed, the metropolis lit as a town / at Earth's area

# ---- Jura Island (register, 15 October), with the measures the layout needs
JURA_REGION = 'L09'
IRIDUM = (45.0101, -31.6653, 249.28998)   # the IAU gazetteer's Sinus Iridum: latitude, east longitude, diameter in km
JURA_HIGH_M = 1000.0                 # proposal: the high arc is the island's land more than 1 km above its sea
JURA_COAST_KM = 5.0                  # proposal: the outer coasts are the land within 5 km of the sea outside the bay

BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))
RANKS = {1: 'dry highlands poleward of 30 deg', 2: 'high-latitude lowland poleward of 48.6 deg',
         3: 'homestead belt: 30 to 48.6 deg, dry sub-humid to wet, below 4 km',
         4: 'rain belt: within 30 deg, humid and wet', 5: 'dry highlands within 30 deg', 6: 'other land'}

RULES = [
    dict(rule='The land within 32 km of the summit port (5.4 N, 158.6 W) is the metropolis: habitation.',
         state='register', source='decisions.md, The summit port and metropolis, 27 September 2026'),
    dict(rule='Jura Island: the high arc, the outer coasts and the islets are wild and the rest is hybrid, with one harbour '
              'town on the inner bay, which is smaller than a cell and is not drawn.',
         state='register, with proposed measures',
         source='decisions.md, Jura Island, 15 October 2026; the high arc as land more than 1 km above the sea and the '
                'outer coasts as land within 5 km of the sea outside the bay are proposals'),
    dict(rule='The satoyama-like candidate valleys are hybrid land.',
         state='register', source='decisions.md, 16 November 2026: much of the satoyama-like activity lies in the hybrid land'),
    dict(rule='Gardens fill the hybrid share from the drylands, the land nearest a sea or a great lake first.',
         state='register for the dry land, proposal for the order',
         source='decisions.md, 11 and 12 October 2026: megaparks and great gardens built round the lower water; the order '
                'by nearness to water follows Jura\'s garden shores and the resort study\'s cost of lifting water'),
    dict(rule='In every moisture class, and within it in every night regime, the wild keeps its share, taken from the land '
              'ranked last for dense living.',
         state='register for the share of the land and of the dark sky, study for holding it zone by zone, proposal for the '
               'night regimes and the order',
         source='decisions.md, 11 October 2026; the food and distributions studies hold the share in every zone; the night '
                'regimes are the moisture map\'s, so that the land with the darkest nights keeps its wild share too'),
    dict(rule='The land left is open to habitation and its farms, ranked: dry highlands poleward of 30 deg, the '
              'high-latitude lowland, the homestead belt, the rain belt, dry highlands within 30 deg, other land.',
         state='register for the direction, study for the bands',
         source='decisions.md, 11 October 2026 (the drier highlands; homesteads in the middle latitudes; the rest where '
                'it is better sustained); distributions (C\') and land cities for the bands and their limits'),
]


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def rel(path) -> str:
    try:
        return str(Path(path).relative_to(ROOT))
    except ValueError:
        return str(path)


# ---------------------------------------------------------------------------------------------- geometry

def grid(rows, cols):
    """Cell-centre latitudes (from 90 N) and east longitudes (from 0), the cell's height in km and each row's cell area."""
    lat = 90.0 - (np.arange(rows) + 0.5) * 180.0 / rows
    lon = (np.arange(cols) + 0.5) * 360.0 / cols
    cell_km = MOON_RADIUS / 1000.0 * np.radians(180.0 / rows)
    return lat, lon, cell_km, cell_km ** 2 * np.cos(np.radians(lat))


def turned(lat, lon, rotation):
    """Today's latitude and east longitude (deg) to the turned Moon's, longitude in [0, 360)."""
    la, lo = np.radians(np.asarray(lat, float)), np.radians(np.asarray(lon, float))
    p = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1) @ np.asarray(rotation, float).T
    return np.degrees(np.arcsin(np.clip(p[..., 2], -1, 1))), np.degrees(np.arctan2(p[..., 1], p[..., 0])) % 360.0


def t21_cell(tlat, tlon, lat_edges, n_lon):
    """The row and column of the climate model's cell that holds each turned latitude and longitude. Rows run north to
    south between lat_edges; column 0 is centred on longitude 0."""
    rows = np.searchsorted(-np.asarray(lat_edges, float)[1:-1], -np.asarray(tlat, float))
    step = 360.0 / n_lon
    cols = np.floor((np.asarray(tlon, float) + step / 2) / step).astype(np.int64) % n_lon
    return rows, cols


def great_circle_km(lat, lon, lat0, lon0):
    a, b = np.radians(np.asarray(lat, float)), np.radians(lat0)
    c = np.sin(a) * np.sin(b) + np.cos(a) * np.cos(b) * np.cos(np.radians(np.asarray(lon, float) - lon0))
    return MOON_RADIUS / 1000.0 * np.arccos(np.clip(c, -1, 1))


def distance_km(source, lat, cell_km, band_rows=80, reach_km=300.0):
    """The distance in km from each cell to the nearest cell of `source`, capped at reach_km. The grid is worked in bands
    of rows, each with the east-west width of a cell at the band's mean latitude, the columns wrapping round."""
    rows, cols = source.shape
    margin = int(np.ceil(reach_km / cell_km))
    out = np.full(source.shape, reach_km, np.float32)
    for r0 in range(0, rows, band_rows):
        r1 = min(r0 + band_rows, rows)
        a0, a1 = max(0, r0 - margin), min(rows, r1 + margin)
        width = max(float(np.cos(np.radians(lat[r0:r1])).mean()), 1e-3)
        pad = min(cols, int(np.ceil(reach_km / (cell_km * width))))
        block = np.concatenate([source[a0:a1, cols - pad:], source[a0:a1], source[a0:a1, :pad]], axis=1)
        if not block.any():
            continue
        d = ndimage.distance_transform_edt(~block, sampling=(cell_km, cell_km * width))
        out[r0:r1] = np.minimum(d[r0 - a0:r1 - a0, pad:pad + cols], reach_km)
    return out


def great_lakes(lake, row_area, min_km2=GREAT_LAKE_KM2):
    """The lake cells that belong to a lake of min_km2 or more: touching cells are one lake, the columns wrapping round."""
    labels, n = ndimage.label(lake, structure=np.ones((3, 3), bool))
    km2 = np.bincount(labels.ravel(), weights=np.repeat(row_area, lake.shape[1]), minlength=n + 1)
    parent = np.arange(n + 1)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    first, last = labels[:, 0], labels[:, -1]
    for a, b in ((first, last), (first[1:], last[:-1]), (first[:-1], last[1:])):
        both = (a > 0) & (b > 0)
        for i, j in np.unique(np.stack([a[both], b[both]], axis=1), axis=0):
            parent[find(int(i))] = find(int(j))
    root = np.array([find(i) for i in range(n + 1)])
    whole = np.bincount(root, weights=km2, minlength=n + 1)
    big = whole[root] >= min_km2
    big[0] = False
    return big[labels]


def outlets(cells, receiver, sea, steps=4000):
    """The cell each of `cells` drains to in the end: the first sea cell on its way, or the cell where its way stops."""
    p = np.asarray(cells, np.int64).copy()
    for _ in range(steps):
        nxt = receiver[p].astype(np.int64)
        move = (sea[p] < 0) & (nxt >= 0) & (nxt != p)
        if not move.any():
            break
        p[move] = nxt[move]
    return p


# ------------------------------------------------------------------------------------------- the sorting

def habitation_rank(moisture, turned_lat, up_m):
    """Each land cell's rank for dense living and its farms (RANKS), 1 the first taken."""
    lat = np.abs(turned_lat)
    dry, high = moisture <= DRY_MAX, up_m >= HIGH_M
    rank = np.full(moisture.shape, 6, np.int8)
    rank[dry & high & (lat < COMFORT_LAT)] = 5
    rank[(lat < COMFORT_LAT) & (moisture >= 4)] = 4
    rank[~high & (lat >= COMFORT_LAT) & (lat < HIGH_LATITUDE) & (moisture >= 3)] = 3
    rank[~high & (lat >= HIGH_LATITUDE)] = 2
    rank[dry & high & (lat >= COMFORT_LAT)] = 1
    return rank


def best_first(rank, turned_lat, rain, shore_km, first=None):
    """Each cell's place in the order in which habitation takes the land: by rank; within the dry highlands poleward of
    45 deg first; within the rain belt the wettest first; then the land nearest a sea or great lake. `first` marks
    cells put at the head of the rain belt (the candidate valleys, where the mosaics count as farm country)."""
    rank = rank.copy()
    second = np.zeros(rank.shape, np.float32)
    second[(rank == 1) & (np.abs(turned_lat) < POLEWARD_FIRST_LAT)] = 1.0
    wet = rank == 4
    second[wet] = -rain[wet]
    if first is not None:
        rank[first] = 4
        second[first] = -1e9
    place = np.empty(rank.size, np.int64)
    place[np.lexsort((shore_km, second, rank))] = np.arange(rank.size)
    return rank, place


def sort_land(area, moisture, stratum, fixed, candidate, shore_km, place, wild_share, hybrid_share, mosaics='hybrid'):
    """Each land cell's kind (WILD, HYBRID, OPEN) by rules 1 to 4, and how far inland the gardens reach.

    stratum numbers the groups of land that each keep the wild share (here a moisture class in a night regime). fixed
    holds the kinds already set (0 where none). mosaics says where the candidate valleys go: 'hybrid' (inside the hybrid
    share), 'beside' (hybrid, with the full share of gardens beside them) or 'habitation' (open land)."""
    kind = fixed.astype(np.int8).copy()
    total = float(area.sum())
    groups = int(stratum.max()) + 1
    group_km2 = np.bincount(stratum, weights=area, minlength=groups)
    dry_land = moisture <= DRY_MAX
    if mosaics in ('hybrid', 'beside'):
        kind[(kind == 0) & candidate] = HYBRID
    elif mosaics == 'habitation':
        kind[(kind == 0) & candidate] = OPEN
    else:
        raise ValueError(mosaics)
    target = hybrid_share * total
    if mosaics == 'beside':
        target += float(area[candidate & (kind == HYBRID)].sum())
    # the gardens: dry land nearest water first, leaving each group its wild share
    held = np.bincount(stratum[kind > WILD], weights=area[kind > WILD], minlength=groups)
    room = np.maximum((1.0 - wild_share) * group_km2 - held, 0.0)
    free = np.flatnonzero((kind == 0) & dry_land)
    free = free[np.argsort(shore_km[free], kind='stable')]
    a, m = area[free], stratum[free]
    within = np.zeros(free.size, bool)
    for c in np.unique(m):
        sel = m == c
        within[sel] = np.cumsum(a[sel]) <= room[c]
    take = within & (np.cumsum(np.where(within, a, 0.0)) <= target - float(area[kind == HYBRID].sum()))
    kind[free[take]] = HYBRID
    reach = float(shore_km[free[take]].max()) if take.any() else 0.0
    # the wild: each group's share, from the land ranked last for habitation
    for c in range(groups):
        in_c = stratum == c
        if not in_c.any():
            continue
        need = wild_share * group_km2[c] - float(area[in_c & (kind == WILD)].sum())
        free = np.flatnonzero(in_c & (kind == 0))
        free = free[np.argsort(-place[free], kind='stable')]
        kind[free[np.cumsum(area[free]) <= need]] = WILD
    kind[kind == 0] = OPEN
    return kind, reach


def place_people(area, kind, band, place, needs):
    """The open cells the distribution fills: in each band the best-placed open land up to the band's need. Returns the
    cells taken and, per band, the open land, the need, the land placed and the shortfall, in km2."""
    used = np.zeros(area.size, bool)
    out = {}
    for b, need in needs.items():
        free = np.flatnonzero((kind == OPEN) & (band == b))
        free = free[np.argsort(place[free], kind='stable')]
        take = np.cumsum(area[free]) <= need
        used[free[take]] = True
        have = float(area[free].sum())
        out[b] = dict(open_km2=round(have), need_km2=round(need), placed_km2=round(float(area[free[take]].sum())),
                      short_km2=round(max(need - have, 0.0)))
    return used, out


def band_needs(wild_share, end):
    """C' needs by band in km2 at the low (0) or high (1) end, on the figures of NEED: the dry highlands' tower districts
    and towns; the high-latitude lowland's scattered communities with those who leave the middle latitudes; the homestead
    belt's heritage pasture; the rain belt's fields with its scattered communities."""
    fields = NEED['fields']['50' if wild_share >= 0.5 else '30'][end]
    scattered = NEED['scattered_communities'][end]
    return {1: NEED['tower_districts'][end] + NEED['highland_towns'][end],
            2: (1.0 - SCATTERED_RAIN_BELT) * scattered,
            3: NEED['heritage_pasture'][end],
            4: fields + SCATTERED_RAIN_BELT * scattered}


def shares(area, kind, used, group, n):
    """Per group: land in km2 and the shares of wild, hybrid, open land and of the open land the distribution fills."""
    land = np.bincount(group, weights=area, minlength=n)
    out = [land]
    for k in (WILD, HYBRID, OPEN):
        out.append(np.bincount(group[kind == k], weights=area[kind == k], minlength=n))
    out.append(np.bincount(group[used], weights=area[used], minlength=n))
    return np.stack(out)


def share_rows(table, keys):
    rows = {}
    for i, key in enumerate(keys):
        land = float(table[0, i])
        if key is None or land <= 0:
            continue
        rows[key] = dict(land_km2=round(land), wild=round(float(table[1, i]) / land, 4),
                         hybrid=round(float(table[2, i]) / land, 4), open=round(float(table[3, i]) / land, 4),
                         filled=round(float(table[4, i]) / land, 4))
    return rows


# ------------------------------------------------------------------------------------------------ inputs

def load(log=print):
    """The design surface's land cells with everything the rules read, as flat arrays over the land cells."""
    started = time.time()
    with np.load(SURFACE, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        height = z['height_m'].astype(np.float32)
        water_cell, water_body = z['water_cell'], z['water_body']
    rows, cols = height.shape
    body = np.full(rows * cols, -1, np.int8)
    body[water_cell] = water_body
    body = body.reshape(rows, cols)
    land = body < 0
    lat, lon, cell_km, row_area = grid(rows, cols)
    li = np.flatnonzero(land)
    r, c = li // cols, li % cols
    area = row_area[r]
    log(f'land {area.sum() / 1e6:.2f} million km2 in {li.size:,} cells; lakes {row_area[(body == 0).nonzero()[0]].sum() / 1e6:.2f}, '
        f'seas {row_area[(body > 0).nonzero()[0]].sum() / 1e6:.2f} million km2, {time.time() - started:.0f} s')

    # the climate's cells in the turned frame
    rotation = np.asarray(meta['frame']['rotation_today_to_turned'], float)
    tlat = np.empty(li.size, np.float32)
    t_row = np.empty(li.size, np.int16)
    t_col = np.empty(li.size, np.int16)
    record = json.loads(MOISTURE.read_text())
    g_lat, g_lon = np.asarray(record['grid']['lat'], float), np.asarray(record['grid']['lon'], float)
    for s in range(0, li.size, 2_000_000):
        e = min(s + 2_000_000, li.size)
        la, lo = turned(lat[r[s:e]], lon[c[s:e]], rotation)
        tlat[s:e] = la
        t_row[s:e], t_col[s:e] = t21_cell(la, lo, record['grid']['lat_edges'], g_lon.size)
    cells = record['cells']
    c_row, c_col = np.asarray(cells['row']), np.asarray(cells['col'])
    own = np.zeros((g_lat.size, g_lon.size), bool)
    own[c_row, c_col] = True
    gla, glo = np.radians(np.repeat(g_lat, g_lon.size)), np.radians(np.tile(g_lon, g_lat.size))
    p = np.stack([np.cos(gla) * np.cos(glo), np.cos(gla) * np.sin(glo), np.sin(gla)], axis=1)
    flat_land = np.flatnonzero(own.ravel())
    nearest = flat_land[np.argmax(p @ p[flat_land].T, axis=1)]          # each cell's nearest land cell (itself on land)
    index_of = np.full(own.size, -1)
    index_of[c_row * g_lon.size + c_col] = np.arange(c_row.size)
    source = index_of[nearest]                                           # the entry of cells each T21 cell reads
    cell = t_row.astype(np.int64) * g_lon.size + t_col
    entry = source[cell]
    moisture = np.asarray(cells['moisture'], np.int8)[entry]
    night = np.asarray(cells['night_regime'], np.int8)[entry]
    rain = np.asarray(cells['rain_mm_yr'], np.float32)[entry]
    air = np.asarray(cells['air_c'], np.float32)[entry]
    warmest = np.asarray(cells['warmest_3day_air_c'], np.float32)[entry]
    assigned = ~own.ravel()[cell]
    log(f'moisture classes carried to the land cells, {assigned.mean():.3f} of the cells assigned from the nearest land '
        f'cell, {time.time() - started:.0f} s')

    # the atlas's regions and districts
    with np.load(DISTRICTS_GRID, allow_pickle=False) as z:
        k = rows // z['region'].shape[0]
        region = np.repeat(np.repeat(z['region'], k, axis=0), k, axis=1).ravel()[li].astype(np.int32)
        district = np.repeat(np.repeat(z['district'], k, axis=0), k, axis=1).ravel()[li].astype(np.int32)
    atlas = json.loads(DISTRICTS.read_text())

    # the satoyama-like candidate valleys
    with np.load(SATOYAMA_UNITS, allow_pickle=False) as z:
        unit = z['unit'][li]
        candidate = np.zeros(li.size, bool)
        in_unit = unit >= 0
        candidate[in_unit] = z['candidate'][unit[in_unit]]
    satoyama = json.loads(SATOYAMA.read_text())
    cand_km2 = float(area[candidate].sum())
    if abs(cand_km2 - satoyama['candidates']['km2']) > 0.01 * satoyama['candidates']['km2']:
        raise ValueError(f'the candidate valleys cover {cand_km2:,.0f} km2 here and '
                         f'{satoyama["candidates"]["km2"]:,} km2 in {SATOYAMA.name}: the units are not that record\'s')

    # nearness to a sea or a great lake
    lakes_big = great_lakes(body == 0, row_area)
    shore = distance_km((body > 0) | lakes_big, lat, cell_km).ravel()[li]
    log(f'great lakes {row_area[lakes_big.nonzero()[0]].sum() / 1e6:.2f} million km2; shores measured, {time.time() - started:.0f} s')

    with np.load(ROUTING, allow_pickle=False) as z:
        routing_meta = json.loads(str(z['metadata']))
        receiver, sea = z['receiver'], z['sea']
    up = height.ravel()[li] - np.float32(meta['common_level_m'])
    return dict(meta=meta, rows=rows, cols=cols, lat=lat, lon=lon, cell_km=cell_km, row_area=row_area, body=body, li=li,
                r=r, c=c, area=area, tlat=tlat, moisture=moisture, night=night, rain=rain, air=air, warmest=warmest,
                assigned=assigned, region=region, district=district, atlas=atlas, candidate=candidate, satoyama=satoyama,
                shore=shore, up=up, height=height.ravel()[li], receiver=receiver, sea=sea, routing_meta=routing_meta,
                moisture_record=record)


def summit_cells(d):
    """The land cells within METROPOLIS_KM of the summit, and what the surface holds there."""
    lon0 = SUMMIT[1] % 360.0
    km = great_circle_km(d['lat'][d['r']], d['lon'][d['c']], SUMMIT[0], lon0)
    inside = km <= METROPOLIS_KM
    a = d['area']
    circle = np.pi * METROPOLIS_KM ** 2
    return inside, dict(
        at=dict(lat=SUMMIT[0], lon_east=SUMMIT[1]), radius_km=METROPOLIS_KM, land_km2=round(float(a[inside].sum())),
        circle_km2=round(circle), water_share_of_the_circle=round(1.0 - float(a[inside].sum()) / circle, 3),
        mean_height_above_the_common_level_m=round(float(np.average(d['up'][inside], weights=a[inside]))),
        moisture_class=CLASS_KEYS[int(np.bincount(d['moisture'][inside], weights=a[inside]).argmax())],
        rain_mm_yr=round(float(np.average(d['rain'][inside], weights=a[inside]))),
        turned_latitude_deg=round(float(np.average(d['tlat'][inside], weights=a[inside])), 2),
        rank_without_the_rule=RANKS[int(np.bincount(habitation_rank(d['moisture'][inside], d['tlat'][inside], d['up'][inside]),
                                                    weights=a[inside]).argmax())])


def jura_layout(d, log=print):
    """Jura Island's kinds by the register's layout. Returns the fixed kinds over the land cells (0 elsewhere) and the
    island's measures."""
    number = int(JURA_REGION[1:])
    cols = d['cols']
    here = np.flatnonzero(d['region'] == number)
    r, c = d['r'][here], d['c'][here]
    r0, c0 = r.min(), c.min()
    box = np.zeros((r.max() - r0 + 1, c.max() - c0 + 1), bool)
    box[r - r0, c - c0] = True
    labels, n = ndimage.label(box, structure=np.ones((3, 3), bool))
    lab = labels[r - r0, c - c0]
    a = d['area'][here]
    km2 = np.bincount(lab, weights=a, minlength=n + 1)
    main = lab == int(km2.argmax())
    sea_level = float(d['routing_meta']['levels_m']['near_side_main'])
    above = d['height'][here] - sea_level
    end = outlets(d['li'][here], d['receiver'], d['sea'])
    end_lat, end_lon = d['lat'][end // cols], d['lon'][end % cols]
    to_bay = (great_circle_km(end_lat, end_lon, IRIDUM[0], IRIDUM[1] % 360.0) <= IRIDUM[2] / 2) & (d['sea'][end] >= 0)
    high = above > JURA_HIGH_M
    outer_coast = ~to_bay & (d['shore'][here] <= JURA_COAST_KM)
    wild = ~main | high | outer_coast
    fixed = np.zeros(d['area'].size, np.int8)
    fixed[here] = np.where(wild, WILD, HYBRID)
    total = float(a.sum())
    island = float(a[main].sum())
    measures = dict(
        region=JURA_REGION, land_km2=round(total), island_km2=round(island), islets=int(n - 1),
        islets_km2=round(total - island), summit_above_its_sea_m=round(float(above.max())),
        sea_level_above_the_geoid_m=round(sea_level, 1),
        draining_to_the_bay_share_of_the_island=round(float(a[main & to_bay].sum()) / island, 4),
        draining_to_a_circle_1_15_times_the_bays_share=round(float(a[main & (d['sea'][end] >= 0) & (
            great_circle_km(end_lat, end_lon, IRIDUM[0], IRIDUM[1] % 360.0) <= 1.15 * IRIDUM[2] / 2)].sum()) / island, 4),
        above_1_km_share_of_the_island=round(float(a[main & high].sum()) / island, 4),
        outer_coast_belt_share_of_the_island=round(float(a[main & outer_coast & ~high].sum()) / island, 4),
        named_garden_share_of_the_island=round(float(a[main & to_bay & ~high].sum()) / island, 4),
        unnamed_share_of_the_island=round(float(a[main & ~to_bay & ~high & ~outer_coast].sum()) / island, 4),
        hybrid_share_of_the_region=round(float(a[~wild].sum()) / total, 4),
        wild_share_of_the_region=round(float(a[wild].sum()) / total, 4),
        moisture=share_of(d['moisture'][here], a),
        note=('The register names the bay-facing slopes and shores (garden), the high arc, the outer coasts and the '
              'islets (wild) and a harbour town; the island\'s land that drains away from the bay, below the high arc '
              'and inland of the outer coast, is not named and is counted as hybrid here on the row\'s "mostly hybrid".'))
    log(f"Jura: {measures['island_km2']:,} km2 island, {measures['draining_to_the_bay_share_of_the_island']:.0%} draining to "
        f"the bay, {measures['hybrid_share_of_the_region']:.0%} hybrid by the layout")
    return fixed, measures


def share_of(index, area):
    km2 = np.bincount(index, weights=area, minlength=len(CLASS_KEYS))
    return {k: round(float(v) / float(area.sum()), 4) for k, v in zip(CLASS_KEYS, km2)}


def strata(d):
    """The groups of land that each keep the wild share: every moisture class within every night regime."""
    return d['moisture'].astype(np.int64) * len(d['moisture_record']['night_regimes']) + d['night']


# ------------------------------------------------------------------------------------------------- build

VARIANTS = (
    dict(key='A', title='30% wild beside the hybrid fifth, the mosaics in the hybrid land', wild=0.30, hybrid=HYBRID_SHARE,
         mosaics='hybrid', open_items='the register\'s compromise (11 October) and its row on the mosaics (16 November)'),
    dict(key='B', title='half wild beside the hybrid fifth, the mosaics in the hybrid land', wild=0.50, hybrid=HYBRID_SHARE,
         mosaics='hybrid', open_items='the register\'s preference for half wild, with the hybrid fifth of the comparison run'),
    dict(key='C', title='30% wild beside the hybrid fifth, the mosaics with habitation', wild=0.30, hybrid=HYBRID_SHARE,
         mosaics='habitation', open_items='the managed biomes placed in the land of full habitation'),
    dict(key='A+', title='30% wild, the mosaics beside a full fifth of gardens', wild=0.30, hybrid=HYBRID_SHARE,
         mosaics='beside', open_items='the hybrid fifth grown to hold the mosaics whole (16 November, open)'),
)


def run_variant(d, v, fixed, rank, log=print):
    a = d['area']
    first = d['candidate'] if v['mosaics'] == 'habitation' else None
    rank_v, place = best_first(rank, d['tlat'], d['rain'], d['shore'], first=first)
    kind, reach = sort_land(a, d['moisture'], strata(d), fixed, d['candidate'], d['shore'], place, v['wild'], v['hybrid'], v['mosaics'])
    band = np.where(rank_v <= 4, rank_v, 0).astype(np.int8)
    band[fixed == OPEN] = 0
    rooms = {}
    for end, name in ((0, 'low'), (1, 'high')):
        used, table = place_people(a, kind, band, place, band_needs(v['wild'], end))
        rooms[name] = {RANKS[b]: row for b, row in table.items()}
    used |= fixed == OPEN                                                 # the metropolis
    total = float(a.sum())
    open_km2 = float(a[kind == OPEN].sum())
    receivers = [NEED['settled'][i] * NEED['receivers_share_of_settled'][i] for i in (0, 1)]
    need = [NEED['metropolis'][i] + sum(band_needs(v['wild'], i).values()) for i in (0, 1)]
    groups = {}
    groups['moisture_class'] = share_rows(shares(a, kind, used, d['moisture'], len(CLASS_KEYS)), CLASS_KEYS)
    regimes = d['moisture_record']['night_regimes']
    groups['night_regime'] = share_rows(shares(a, kind, used, d['night'], len(regimes)), regimes)
    lat = np.abs(d['tlat'])
    band_index = np.searchsorted([b[1] for b in BANDS[:-1]], lat, side='right')
    groups['turned_latitude'] = share_rows(shares(a, kind, used, band_index, len(BANDS)), [f'{lo}-{hi}' for lo, hi in BANDS])
    groups['rank'] = share_rows(shares(a, kind, used, rank.astype(np.int64), 7), [None] + [RANKS[i] for i in range(1, 7)])
    codes = [None] * 200
    for reg in d['atlas']['regions']:
        if reg['kind'] == 'land':
            codes[int(reg['code'][1:])] = reg['code']
    region = np.where((d['region'] > 0) & (d['region'] < 100), d['region'], 0)
    groups['region'] = share_rows(shares(a, kind, used, region, 200), codes)
    outside = float(a[region == 0].sum())
    n_land = len(d['atlas']['districts'])
    dist = np.where((d['district'] > 0) & (d['district'] <= n_land), d['district'], 0)
    table = shares(a, kind, used, dist, n_land + 1)
    districts = {}
    for i, dd in enumerate(d['atlas']['districts']):
        land = float(table[0, i + 1])
        if land > 0:
            districts[dd['code']] = [round(float(table[k, i + 1]) / land, 3) for k in (1, 2, 3, 4)]
    # the same land sorted by whole districts: each district takes the kind that holds most of it
    d_land, d_kinds = table[0, 1:], table[1:4, 1:]
    has = d_land > 0
    lead = d_kinds.argmax(axis=0)
    mixed = has & (d_kinds.max(axis=0) <= 0.5 * d_land)
    whole = dict(
        shares={KINDS[k + 1]: round(float(d_land[has & (lead == k)].sum()) / float(d_land[has].sum()), 4) for k in range(3)},
        districts={KINDS[k + 1]: int((has & (lead == k)).sum()) for k in range(3)},
        no_kind_over_half=dict(districts=int(mixed.sum()), land_share=round(float(d_land[mixed].sum()) / float(d_land[has].sum()), 4)),
        land_in_its_districts_kind=round(float(d_kinds.max(axis=0)[has].sum()) / float(d_land[has].sum()), 4))
    hybrid = kind == HYBRID
    garden = hybrid & ~d['candidate'] & (fixed == 0)
    shore_belt = {}
    for name, part in (('drylands', d['moisture'] <= DRY_MAX), ('humid_and_wet', d['moisture'] > DRY_MAX)):
        belt = part & (d['shore'] <= SHORE_BELT_KM)
        km2 = float(a[belt].sum())
        shore_belt[name] = dict(km2=round(km2), **{KINDS[k]: round(float(a[belt & (kind == k)].sum()) / km2, 4) for k in KINDS})
    wild_by_class = {k: row['wild'] for k, row in groups['moisture_class'].items()}
    dark = groups['night_regime'].get('dark', {})
    short = {name: {b: row['short_km2'] for b, row in rows.items() if row['short_km2'] > 0} for name, rows in rooms.items()}
    fails, notes = [], []
    if min(wild_by_class.values()) < v['wild'] - 0.002:
        fails.append('the wild share in ' + ', '.join(k for k, s in wild_by_class.items() if s < v['wild'] - 0.002))
    if float(a[kind == WILD].sum()) / total < 0.5 - 0.002:
        notes.append('keeps the compromise\'s 30% of the land wild, where the register prefers at least half')
    if dark and dark['wild'] < v['wild'] - 0.002:
        fails.append(f"the wild share of the dark-night land ({dark['wild']:.0%})")
    lit = [x * OCTOBER_LAND_KM2 / total for x in LIT_SKY_SHARE_OF_OCTOBER_LAND]
    unlit = [round(1.0 - x - float(a[hybrid].sum()) / total, 4) for x in lit]
    if min(unlit) < v['wild'] - 0.002:
        fails.append(f'the share of the land under an unlit sky with every park lit ({min(unlit):.0%})')
    if abs(float(a[hybrid].sum()) / total - HYBRID_SHARE) > 0.01:
        fails.append(f"a hybrid share of one fifth ({float(a[hybrid].sum()) / total:.1%})")
    for name in ('low', 'high'):
        if short[name]:
            fails.append(f"room at the {name} end in " + '; '.join(f'{b} ({km2:,} km2 short)' for b, km2 in short[name].items()))
    out = dict(
        title=v['title'], wild_share_set=v['wild'], hybrid_share_set=v['hybrid'], mosaics=v['mosaics'], differs_in=v['open_items'],
        moon=dict(land_km2=round(total), wild=round(float(a[kind == WILD].sum()) / total, 4),
                  hybrid=round(float(a[hybrid].sum()) / total, 4), open=round(open_km2 / total, 4),
                  filled_at_the_high_end=round(float(a[used].sum()) / total, 4),
                  wild_km2=round(float(a[kind == WILD].sum())), hybrid_km2=round(float(a[hybrid].sum())), open_km2=round(open_km2)),
        hybrid=dict(candidate_valleys_km2=round(float(a[hybrid & d['candidate']].sum())), gardens_km2=round(float(a[garden].sum())),
                    gardens_reach_inland_km=round(reach, 1),
                    gardens_by_class={k: round(float(x)) for k, x in zip(CLASS_KEYS, np.bincount(d['moisture'][garden], weights=a[garden], minlength=len(CLASS_KEYS)))},
                    m2_a_resident_at_8_billion=round(float(a[hybrid].sum()) * 1e6 / 8e9)),
        land_within_5_km_of_a_sea_or_great_lake=shore_belt,
        room=dict(need_km2=[round(x) for x in need], need_share_of_the_land=[round(x / total, 4) for x in need],
                  receivers_km2=[round(x) for x in receivers],
                  need_with_receivers_share_of_the_land=[round((x + y) / total, 4) for x, y in zip(need, receivers)],
                  open_land_not_needed_share_of_the_land=[round((open_km2 - x - y) / total, 4) for x, y in zip(need, receivers)],
                  land_without_use_share=[round(1.0 - float(a[hybrid].sum()) / total - (x + y) / total, 4) for x, y in zip(need, receivers)],
                  by_band=rooms),
        dark_sky=dict(lit_by_settlement_share_of_the_land=[round(x, 4) for x in lit],
                      unlit_share_with_every_park_lit=unlit, wild_share_of_the_dark_night_land=dark.get('wild'),
                      filled_share_of_the_dark_night_land=dark.get('filled')),
        by=groups, districts=districts, by_whole_districts=whole, land_outside_the_regions_km2=round(outside),
        fails=fails, notes=notes)
    log(f"variant {v['key']}: wild {out['moon']['wild']:.3f}, hybrid {out['moon']['hybrid']:.3f}, open {out['moon']['open']:.3f}, "
        f"filled {out['moon']['filled_at_the_high_end']:.3f}; gardens reach {reach:.0f} km inland; fails: {fails or 'nothing'}")
    return out, kind, used


def access(d, kind, used):
    """How far the hybrid land lies from the open land the distribution fills, as shares of the hybrid land."""
    full = np.zeros(d['rows'] * d['cols'], bool)
    full[d['li'][used]] = True
    km = distance_km(full.reshape(d['rows'], d['cols']), d['lat'], d['cell_km']).ravel()[d['li']]
    h = kind == HYBRID
    a = d['area'][h]
    return dict(within_a_sky_boats_trip=round(float(a[km[h] <= SKY_BOAT_KM].sum() / a.sum()), 4),
                within_a_sky_ferrys_route=round(float(a[km[h] <= SKY_FERRY_KM].sum() / a.sum()), 4),
                beyond_300_km=round(float(a[km[h] >= 300.0].sum() / a.sum()), 4))


def heritage(d, kinds):
    """Each heritage site's place on the design surface: under water, or on land with its kind in each variant."""
    sites = json.loads(HERITAGE.read_text())['sites']
    cols = d['cols']
    where = np.full(d['rows'] * cols, -1, np.int64)
    where[d['li']] = np.arange(d['li'].size)
    bodies = d['meta']['bodies']
    out = []
    for s in sites:
        r = min(int((90.0 - s['lat']) * d['rows'] / 180.0), d['rows'] - 1)
        c = int((s['lon'] % 360.0) * cols / 360.0) % cols
        i = where[r * cols + c]
        row = dict(id=s['id'], name=s['name'], lat=s['lat'], lon=s['lon'], status_on_the_28_percent_atlas=s['status'])
        if i < 0:
            row['on_the_design_surface'] = 'under ' + bodies[int(d['body'][r, c])]
        else:
            row['on_the_design_surface'] = 'land'
            row['moisture_class'] = CLASS_KEYS[int(d['moisture'][i])]
            row['kind'] = {key: ('filled' if used[i] else KINDS[int(kind[i])]) for key, (kind, used) in kinds.items()}
        out.append(row)
    return out


def build(log=print):
    started = time.time()
    d = load(log)
    a = d['area']
    total = float(a.sum())
    rank = habitation_rank(d['moisture'], d['tlat'], d['up'])
    summit, summit_measures = summit_cells(d)
    jura_fixed, jura = jura_layout(d, log)
    fixed = jura_fixed.copy()
    fixed[summit] = OPEN
    variants, kinds = {}, {}
    grids = {}
    for v in VARIANTS:
        out, kind, used = run_variant(d, v, fixed, rank, log)
        out['access'] = access(d, kind, used)
        variants[v['key']] = out
        kinds[v['key']] = (kind, used)
        full = np.zeros(d['rows'] * d['cols'], np.int8)
        full[d['li']] = np.where(kind == OPEN, np.where(used, 3, 4), kind)
        grids[v['key']] = full.reshape(d['rows'], d['cols'])
        log(f"variant {v['key']} done, {time.time() - started:.0f} s")
    # what the general rules give Jura and the summit without their fixed rows, on the first variant
    free = np.zeros(a.size, np.int8)
    v = VARIANTS[0]
    _, place = best_first(rank, d['tlat'], d['rain'], d['shore'])
    kind0, _ = sort_land(a, d['moisture'], strata(d), free, d['candidate'], d['shore'], place, v['wild'], v['hybrid'], v['mosaics'])
    here = d['region'] == int(JURA_REGION[1:])
    jura['by_the_general_rules_alone'] = {KINDS[k]: round(float(a[here & (kind0 == k)].sum()) / float(a[here].sum()), 4) for k in KINDS}
    summit_measures['by_the_general_rules_alone'] = {KINDS[k]: round(float(a[summit & (kind0 == k)].sum()) / float(a[summit].sum()), 4) for k in KINDS}

    lat = np.abs(d['tlat'])
    dry = d['moisture'] <= DRY_MAX
    class_km2 = np.bincount(d['moisture'], weights=a, minlength=len(CLASS_KEYS))
    rank_km2 = np.bincount(rank, weights=a, minlength=7)
    cand = d['candidate']
    band_index = np.searchsorted([b[1] for b in BANDS[:-1]], lat, side='right')
    inventory = dict(
        land_km2=round(total), land_share_of_the_moon=round(total / float(d['row_area'].sum() * d['cols']), 4),
        by_moisture_class={k: dict(km2=round(float(x)), share=round(float(x) / total, 4)) for k, x in zip(CLASS_KEYS, class_km2)},
        moisture_map_shares_over_the_design_land=d['moisture_record']['shares']['over_the_design_land']['by_class'],
        assigned_share=round(float(a[d['assigned']].sum()) / total, 4),
        by_night_regime={k: round(float(x) / total, 4) for k, x in zip(d['moisture_record']['night_regimes'],
                                                                         np.bincount(d['night'], weights=a, minlength=4))},
        by_rank={RANKS[i]: dict(km2=round(float(rank_km2[i])), share=round(float(rank_km2[i]) / total, 4)) for i in range(1, 7)},
        high_ground_km2=round(float(a[d['up'] >= HIGH_M].sum())),
        high_ground_poleward_of_30_km2=round(float(a[(d['up'] >= HIGH_M) & (lat >= COMFORT_LAT)].sum())),
        high_ground_poleward_of_45_km2=round(float(a[(d['up'] >= HIGH_M) & (lat >= POLEWARD_FIRST_LAT)].sum())),
        night_regime_by_class_km2={n: {k: round(float(a[(d['night'] == j) & (d['moisture'] == i)].sum())) for i, k in enumerate(CLASS_KEYS)}
                                   for j, n in enumerate(d['moisture_record']['night_regimes'])},
        dry_high_ground_km2=round(float(a[dry & (d['up'] >= HIGH_M)].sum())),
        dry_high_ground_poleward_of_30_km2=round(float(a[dry & (d['up'] >= HIGH_M) & (lat >= COMFORT_LAT)].sum())),
        dry_high_ground_poleward_of_45_km2=round(float(a[dry & (d['up'] >= HIGH_M) & (lat >= POLEWARD_FIRST_LAT)].sum())),
        drylands_km2=round(float(a[dry].sum())),
        drylands_within_km_of_a_sea_or_great_lake={str(k): round(float(a[dry & (d['shore'] <= k)].sum())) for k in (5, 10, 25, 50, 100)},
        land_within_km_of_a_sea_or_great_lake={str(k): round(float(a[d['shore'] <= k].sum())) for k in (5, 10, 25, 50, 100)},
        candidate_valleys=dict(
            km2=round(float(a[cand].sum())), share_of_the_land=round(float(a[cand].sum()) / total, 4),
            share_of_a_fifth=round(float(a[cand].sum()) / (HYBRID_SHARE * total), 4),
            by_moisture_class={k: round(float(x)) for k, x in zip(CLASS_KEYS, np.bincount(d['moisture'][cand], weights=a[cand], minlength=len(CLASS_KEYS)))},
            by_turned_latitude={f'{lo}-{hi}': dict(km2=round(float(a[cand & (band_index == i)].sum())),
                                                  share_of_the_bands_fifth=round(float(a[cand & (band_index == i)].sum()) / (HYBRID_SHARE * float(a[band_index == i].sum())), 4))
                                for i, (lo, hi) in enumerate(BANDS)},
            woodland_km2=d['satoyama']['candidates']['woodland_km2'], valley_floor_km2=d['satoyama']['candidates']['valley_floor_km2'],
            flooded_possible_km2=d['satoyama']['candidates']['flooded_possible_km2']),
        mean_air_c_by_class={k: round(float(np.average(d['air'][d['moisture'] == i], weights=a[d['moisture'] == i])), 1) for i, k in enumerate(CLASS_KEYS)},
        warmest_3day_air_c_by_class={k: round(float(d['warmest'][d['moisture'] == i].max()), 1) for i, k in enumerate(CLASS_KEYS)})
    fifth = HYBRID_SHARE * total
    carried = dict(
        october_land_km2=OCTOBER_LAND_KM2, october_fifth_km2=OCTOBER_FIFTH_KM2, land_km2=round(total), fifth_km2=round(fifth),
        fifth_change_km2=round(fifth - OCTOBER_FIFTH_KM2), fifth_m2_a_resident_at_8_billion=round(fifth * 1e6 / 8e9),
        need_share_of_october_land=[round((NEED['settled'][i] + NEED['fields']['30'][i] + NEED['heritage_pasture'][i]) / OCTOBER_LAND_KM2, 4) for i in (0, 1)],
        need_share_of_the_land=[round((NEED['settled'][i] + NEED['fields']['30'][i] + NEED['heritage_pasture'][i]) / total, 4) for i in (0, 1)],
        note='Areas the studies give in km2 are carried unchanged; shares of the land are taken again on the present land.')

    BULK.mkdir(parents=True, exist_ok=True)
    surface_hash = sha256(SURFACE)
    surface_record = json.loads(SURFACE_RECORD.read_text())
    recorded = surface_record.get('product', {}).get('sha256')
    if recorded and recorded != surface_hash:
        raise ValueError(f'{SURFACE} is not the product {rel(SURFACE_RECORD)} records')
    meta = dict(schema=PRODUCT_SCHEMA, grid=dict(pixels_per_degree=16, rows=d['rows'], cols=d['cols']),
                reading_rule=('Rows from 90 N and columns from 0 E, 16 per degree, on today\'s Moon. Each variant\'s grid '
                              'holds 0 on water, 1 wild land, 2 hybrid land, 3 open land the distribution fills at its '
                              'high end, 4 open land beyond that need.'),
                variants={v['key']: v['title'] for v in VARIANTS}, surface_sha256=surface_hash)
    np.savez_compressed(PRODUCT, metadata=np.array(json.dumps(meta)), **{f'kind_{k.replace("+", "_plus")}': g for k, g in grids.items()})
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}),
        inputs={rel(p): sha256(p) for p in (SURFACE_RECORD, MOISTURE, SATOYAMA, HERITAGE)},
        ignored_inputs={str(p): (surface_hash if p == SURFACE else sha256(p)) for p in (SURFACE, ROUTING, DISTRICTS, DISTRICTS_GRID, SATOYAMA_UNITS)},
        evidence=('A screen on products: stated rules applied to the design surface\'s land at 16 px/deg, with the moisture '
                  'map\'s classes and night regimes at the climate model\'s cells of about 170 km, the satoyama map\'s '
                  'candidate valleys, the redrawn regions and districts and the design routing. The shares of wild and '
                  'hybrid land are the register\'s; where each kind lies follows rules of which several are proposals '
                  '(rules, parameters). The people\'s needs are the distributions and food studies\' figures for C\' at 8 '
                  'billion with the blend, worked on the land of October and carried in km2. The climate is an estimate '
                  'and nothing here is validated against a design of any place.'),
        reading_rule=('Shares are of land area. wild, hybrid and open add to 1 in every group; filled is the part of the '
                      'open land that C\' at 8 billion takes at its high end and is counted inside open. Regions and '
                      'districts carry the redrawn atlas\'s codes. room.by_band gives, for each of C\'s bands, the open '
                      'land in it, the need, the land placed and the shortfall, at the low and the high end. A variant\'s '
                      'fails lists what it does not meet of the register\'s shares and of the room for the people. '
                      'districts gives [wild, hybrid, open, filled] by district code. by_whole_districts gives what the '
                      'same sorting comes to when each district takes the one kind that holds most of its land: the '
                      'shares of the land, the districts by kind, those where no kind holds more than half, and the share '
                      'of the land that lies in its district\'s kind.'),
        parameters=dict(
            register=dict(hybrid_share=HYBRID_SHARE, wild_shares=WILD_SHARES, summit=SUMMIT, metropolis_km=METROPOLIS_KM),
            study=dict(high_ground_m=HIGH_M, comfort_latitude_deg=COMFORT_LAT, poleward_first_deg=POLEWARD_FIRST_LAT,
                       high_latitude_deg=HIGH_LATITUDE, great_lake_km2=GREAT_LAKE_KM2, sky_boat_km=SKY_BOAT_KM,
                       sky_ferry_km=SKY_FERRY_KM, need_km2=NEED, scattered_share_rain_belt=SCATTERED_RAIN_BELT,
                       scattered_share_middle=SCATTERED_MIDDLE, sinus_iridum=IRIDUM),
            proposal=dict(dry_classes=CLASS_KEYS[:DRY_MAX + 1], gardens_order='nearest a sea or great lake first',
                          wild_order='the land ranked last for habitation first, in each moisture class and night regime',
                          homestead_belt_classes=CLASS_KEYS[3:], rain_belt_classes=CLASS_KEYS[4:],
                          jura_high_arc_m=JURA_HIGH_M, jura_outer_coast_km=JURA_COAST_KM)),
        rules=RULES, ranks=RANKS, inventory=inventory, carried=carried, summit=summit_measures, jura=jura,
        variants=variants, heritage=heritage(d, kinds),
        limits=['The moisture classes, night regimes, rain and warmth are the climate model\'s cells, about 170 km across; '
                'a cell the model holds as water takes its nearest land cell\'s values, which errs dry in the lake districts.',
                'Distances to water are measured in bands of 5 degrees at each band\'s mean cell width and capped at 300 km; '
                'a great lake is one of 1,000 km2 or more at 1.9 km cells.',
                'The needs by band are bounds built from the studies\' printed ranges; the studies\' own corners pair the '
                'ranges differently, and their totals are in carried.',
                'A town, a crater garden town or a harbour is smaller than a cell: the open land is where such places may '
                'stand, and no site is chosen.'],
        product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)), seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    record = build()
    OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    print(f'wrote {OUT}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
