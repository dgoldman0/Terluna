"""Made sand seas on the design surface: where they can lie, the sand and haulage they take, and what the winds make
of them in ten thousand years (20 November 2026).

    # with the system Python, after sand.py has written its product (about four minutes on four capped CPUs, 5 GB):
    ulimit -v 12000000; OPENBLAS_NUM_THREADS=1 python3 -m research.studies.atlas_round.sand_seas

The paper is sand_seas.md beside this file. The sand study (sand.py, sand.md) found that the design's winds move sand
easily and that the land offers them too little for a sand sea within ten thousand years; the register's row of
20 November has one to three seas made in the hyper-arid land from the sand in the sediment that the terrain's upkeep
returns uphill, with the wild land's natural sand left open. This module sites them:

1. The ground (part A): the hyper-arid land as patches, with each patch's winds and relief.
2. The sand at hand (part B): where the upkeep's loop gathers durable sand, at the ends of the land's drainage in
   lakes and in dry hollows, and how far it lies from each site.
3. The seas (part C): each sea's outline by a stated rule on the ground and the winds, its sand, the years and the
   haulage to deliver it, its dunes, what its margins lose downwind and what it covers.
4. The map (part D): the sand study's open-sand case with the seas set down on it, as a product the globe's look
   reads in place of sand_16ppd.npz.

Every figure is a reckoning from the sand study's fields and inherits its evidence state: nothing here is observed.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from shared.constants import MOON_RADIUS

from . import sand as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'sand_seas.json'
BULK = S.DRIVE / 'atlas-round' / 'sand_seas'
PRODUCT = BULK / 'sand_seas_16ppd.npz'
HERITAGE = ROOT / 'research' / 'studies' / 'conservation' / 'results' / 'heritage_register.json'
SCHEMA = 'terluna.research.atlas-round-sand-seas/1'
CODE = ('research/studies/atlas_round/sand_seas.py', 'research/studies/atlas_round/sand.py')
ROWS, COLS, PPD, CELL_M = S.ROWS, S.COLS, S.PPD, S.CELL_M

HYPER_ARID = S.CLASS_EDGES[0]      # the moisture map's divide between the hyper-arid and the arid land (rain over demand, 0.05)
YEARS = 10000.0                    # the land's age in the sand study's central case; the seas are reckoned as standing that long
HAUL_KJ_PER_T_KM = 107.0           # rail freight's energy, as the terrain study prices its loop's haulage (33 billion t-km at 0.11 GW)
LOOP_RANGE_GW = (0.2, 11.0)        # the whole loop's power across the terrain study's erosion range
RADII_KM = (250.0, 500.0, 1000.0, 2000.0)
LEVEL_M = 60.0                     # ground counts as level where the heights within about 10 km spread by less than this
LEVEL_WINDOW = 11                  # cells (21 km)
TAPER_CELLS = 4.0                  # a sea thins to its edge over this many cells (7.6 km)
TAPER_EDGE = 0.15                  # share of the full thickness at the outermost cell
TAPER_SIDES = (0.35, 2.2)          # the taper's width on a sea's upwind and downwind sides, as shares of TAPER_CELLS
WATER_BUFFER_CELLS = 4             # a sea's outline keeps this far from standing water (7.6 km)
CHANNEL_BUFFER_CELLS = 3           # and this far from a channel carrying CHANNEL_M3S or more
CHANNEL_M3S = 20.0
CONVERGENCE_SIGMA = 5.0            # cells (9.5 km): the width the net drift's convergence is smoothed over for an outline
OUTLINE_SIGMA = 2.0                # cells (3.8 km): an outline is rounded off over this width, so that no single cell's height or score notches it
SURFACE_SIGMA = 8.0                # cells (15 km): a sea's own surface follows the ground smoothed over this width, so its sand fills hollows and thins over swells
THICKNESS_MOST = 2.5               # of a sea's mean thickness: the most sand a hollow under it holds
SECTORS = 8                        # a sea's sand is told apart by the side it comes from, in this many directions
EDGE_FROM = 0.8                    # of the hyper-arid divide's dryness: from here to the divide an outline by convergence is drawn back from the edge

# The seas. `seed` is a place on today's Moon inside the sea; `rule` says what draws its outline round the seed:
# 'convergence' takes the ground where the winds' net drift converges fastest, less `swell_m_per_kyr` for every mean
# thickness the ground stands above the ground round it (sand gathers in the lows and leaves the swells) and less up to
# `edge_m_per_kyr` towards the hyper-arid land's edge, from EDGE_FROM of the divide's dryness on, and 'floor'
# the lowest ground of a walled plain. `area_km2` is the outline's size, `thickness_m` the sand's mean equivalent
# thickness inside the taper and `deliver_years` the span its haulage plan is reckoned for.
SEAS = {
    'plain': dict(label='the basalt plain', seed=(23.2, 318.2), rule='convergence', area_km2=20000.0, thickness_m=25.0, reach_deg=9.0,
                  relief_m=150.0, deliver_years=5000.0, swell_m_per_kyr=0.10, edge_m_per_kyr=0.25),
    'walled': dict(label='the walled plain', seed=(70.6, 125.6), rule='floor', area_km2=5000.0, thickness_m=30.0, reach_deg=4.0,
                   relief_m=None, deliver_years=5000.0),
}
SIZES_KM2 = (5000.0, 10000.0, 20000.0, 40000.0)   # other sizes each sea's haulage is reckoned for, at its own thickness
DELIVER_YEARS = (1000.0, 2000.0, 5000.0, 10000.0)


# ---- the ground ---------------------------------------------------------------------------------------------------
def spread(height):
    """How far the heights within LEVEL_WINDOW cells spread (their standard deviation, m)."""
    from scipy import ndimage
    h = height.astype(np.float32)
    mean = ndimage.uniform_filter(h, LEVEL_WINDOW, mode='nearest')
    return np.sqrt(np.maximum(ndimage.uniform_filter(h * h, LEVEL_WINDOW, mode='nearest') - mean * mean, 0.0))


def convergence(rdp, az_deg, land):
    """The net drift's convergence: metres of sand in a thousand years that a flux at the wind's capacity would lay down
    (negative where it would take sand away), from the resultant drift (m2 a year) and its direction, smoothed over
    CONVERGENCE_SIGMA cells. Water counts as no drift."""
    from scipy import ndimage
    a = np.radians(np.nan_to_num(az_deg))
    east, north = np.where(land, rdp * np.sin(a), 0.0), np.where(land, rdp * np.cos(a), 0.0)
    dx = (CELL_M * np.cos(np.radians(S.fine_lat())))[:, None]
    div = (np.roll(east, -1, 1) - np.roll(east, 1, 1)) / (2.0 * np.maximum(dx, 50.0)) + (np.roll(north, 1, 0) - np.roll(north, -1, 0)) / (2.0 * CELL_M)
    div[0], div[-1] = 0.0, 0.0
    return (ndimage.gaussian_filter(np.where(land, -div, 0.0).astype(np.float32), CONVERGENCE_SIGMA, mode=('nearest', 'wrap')) * 1000.0).astype(np.float32)


def cell_of(lat, lon):
    """The row and column of a place on today's Moon."""
    return int(np.clip((90.0 - lat) * PPD, 0, ROWS - 1)), int((lon % 360.0) * PPD) % COLS


def distance_km(lat, lon, lat0, lon0):
    """Great-circle distance (km) from places to one place."""
    p, p0 = np.radians(lat), math.radians(lat0)
    d = np.radians(np.asarray(lon, float) - lon0)
    c = np.sin(p) * math.sin(p0) + np.cos(p) * math.cos(p0) * np.cos(d)
    return MOON_RADIUS / 1e3 * np.arccos(np.clip(c, -1.0, 1.0))


def patches(hyper, F, conv, rough, least_km2=3000.0):
    """The hyper-arid land as joined patches of least_km2 or more: where each lies, its winds and its relief."""
    from scipy import ndimage
    labels, count = ndimage.label(hyper, structure=np.ones((3, 3), bool))
    for r in range(ROWS):                                         # a patch across the seam at 0 E is one patch
        a, b = labels[r, 0], labels[r, -1]
        if a and b and a != b:
            labels[labels == b] = a
    ids = np.unique(labels[labels > 0])
    f = F.fine['central']
    out = []
    for i in ids:
        m = labels == i
        w = F.a2[m].astype(np.float64)
        if w.sum() < least_km2:
            continue
        la, lo = np.broadcast_to(F.lat_f[:, None], m.shape)[m], np.broadcast_to(F.lon_f[None, :], m.shape)[m]
        x = (np.cos(np.radians(la)) * np.cos(np.radians(lo)) * w).sum()
        y = (np.cos(np.radians(la)) * np.sin(np.radians(lo)) * w).sum()
        z = (np.sin(np.radians(la)) * w).sum()
        a = np.radians(f['az'][m])
        east, north = (f['rdp'][m] * np.sin(a) * w).sum() / w.sum(), (f['rdp'][m] * np.cos(a) * w).sum() / w.sum()
        mean = lambda field: float((field[m] * w).sum() / w.sum())   # noqa: E731
        out.append(dict(km2=round(float(w.sum())), at_today=[round(float(np.degrees(np.arctan2(z, np.hypot(x, y)))), 1), round(float(np.degrees(np.arctan2(y, x)) % 360.0), 1)],
                        latitudes=[round(float(la.min()), 1), round(float(la.max()), 1)], aridity=round(mean(F.aridity), 3),
                        drift_potential_m2_yr=round(mean(f['dp']), 1), net_over_total=round(float((f['rdp'][m] * w).sum() / max((f['dp'][m] * w).sum(), 1e-9)), 2),
                        net_towards_deg=round(float(np.degrees(np.arctan2(east, north)) % 360.0)), net_m2_yr=round(float(np.hypot(east, north)), 1),
                        convergence_m_per_kyr=round(mean(conv), 3), basalt_share=round(mean(F.mafic), 2),
                        level_share=round(float(w[rough[m] < LEVEL_M].sum() / w.sum()), 2), spread_m_median=round(float(np.median(rough[m]))),
                        height_m=[round(float(np.percentile(F.height[m], q))) for q in (10, 50, 90)]))
    return sorted(out, key=lambda p: -p['km2'])


# ---- the sand at hand ---------------------------------------------------------------------------------------------
def follow(receiver, stays):
    """Where each cell's way down ends: receivers followed from every cell to the first cell that stays (one marked in
    `stays`, or one that is its own receiver). The step is doubled at each turn until nothing moves, so a way of n
    cells takes about log2(n) turns."""
    nxt = np.where(stays, np.arange(receiver.size, dtype=receiver.dtype), receiver)
    for _ in range(64):
        new = nxt[nxt]
        if np.array_equal(new, nxt):
            break
        nxt = new
    return nxt


def ends_of_drainage(F, routing=None):
    """Each land cell's end of drainage: the first cell on its way down that stands beside standing water or is a dry
    hollow's lowest point (the sand study's mouths). Returns the end's cell for every cell (its own for water and for
    the ends themselves) and the mask of the ends."""
    routing = np.load(S.ROUTING, allow_pickle=False) if routing is None else routing
    n = ROWS * COLS
    flat_water = F.water.reshape(-1)
    receiver = np.clip(routing['receiver'], 0, n - 1).astype(np.int32)
    end = (~flat_water & flat_water[receiver]) | (routing['keeps'] & ~flat_water)
    return follow(receiver, flat_water | end), end


def sources(F, supply='central'):
    """Where the terrain's upkeep gathers sand: every end of the land's drainage with the durable sand its own ground
    sends there in a year. The yield is the sand study's (the terrain study's lowering by latitude and rain, four
    tenths sand-sized, half of that lasting); each land cell's sand goes to the first standing water or dry hollow
    floor on its way down, so no ground is counted twice and the ends together hold what the land sheds. An end is a
    lake's delta where the lake spills, a closed hollow's floor (a lake that keeps its water, or dry ground) or a
    sea's shore; the loop dredges the first, digs the second and leaves the third. Returns the ends' cells and places,
    the bulk cubic metres a year, the kind (0 lake delta, 1 closed hollow, 2 sea shore), the basalt's share of the
    sand by the rock of the ground it came from, and what the sand study's own rule (each end credited with up to its
    cap of the ground upstream) gives the same ends."""
    sup = S.SUPPLY[supply]
    routing = np.load(S.ROUTING, allow_pickle=False)
    surface = np.load(S.SURFACE, allow_pickle=False)
    n = ROWS * COLS
    flat_water = F.water.reshape(-1)
    nxt, end = ends_of_drainage(F, routing)
    shed = (F.yield_row[:, None] * F.rain_rule * sup['erosion'] * sup['durable'] * F.a2 * 1e6).reshape(-1)   # bulk m3 of durable sand a year off each cell
    land = ~flat_water
    reaches = land & end[nxt]
    at_end = np.bincount(nxt[reaches], weights=shed[reaches], minlength=n)
    dark = np.bincount(nxt[reaches], weights=(shed * F.mafic.reshape(-1))[reaches], minlength=n)
    body = np.full(n, -1, np.int16)
    spills = np.zeros(n, bool)
    body[surface['water_cell']] = surface['water_body']
    spills[surface['water_cell']] = surface['water_spills']
    receiver = np.clip(routing['receiver'], 0, n - 1)
    cells = np.flatnonzero(end)
    into = flat_water[receiver[cells]]
    kind = np.where(~into, 1, np.where(body[receiver[cells]] > 0, 2, np.where(spills[receiver[cells]], 0, 1))).astype(np.int8)
    rows, cols = cells // COLS, cells % COLS
    capped = (F.yield_row[rows] * F.rain_rule.reshape(-1)[cells] * sup['erosion'] * sup['durable']
              * np.minimum(routing['drained_km2'][cells], sup['catchment_cap_km2']) * 1e6)
    return dict(cell=cells, lat=F.lat_f[rows], lon=F.lon_f[cols], m3_yr=at_end[cells], kind=kind,
                basalt=dark[cells] / np.maximum(at_end[cells], 1e-9), land_yield_m3_yr=float(shed[land].sum()),
                unreached_m3_yr=float(shed[land & ~reaches].sum()), capped_m3_yr=capped)


def wind_offer(F, src, variant='central'):
    """The rivers' sand the wind is offered in a case of the sand study (km3 over the case's years, after the plants,
    the crust and the calm), by that study's rule of a capped catchment and by the ends' own ground: the measure of
    what the capped rule's double counting does to the sand study's supply."""
    var = S.VARIANTS[variant]
    kept = 0.0 if var['open_sand'] else F.kept_crust.reshape(-1)[src['cell']]
    take = 1.0 - kept * (1.0 - S.FRESH_TAKE)
    windy = np.clip(F.fine[var['subgrid']]['above'].reshape(-1)[src['cell']] / 0.01, 0.0, 1.0)
    factor = var['years'] * F.through_plants.reshape(-1)[src['cell']] * take * windy / 1e9
    return dict(capped_km3=round(float((src['capped_m3_yr'] * factor).sum()), 1), own_ground_km3=round(float((src['m3_yr'] * factor).sum()), 1))


def supply_totals(F, src):
    """The ends' sand by kind (km3 a year), what the land sheds, and the sand study's capped rule beside them."""
    names = ('lake deltas', 'closed hollows', 'sea shores')
    by = {names[k]: round(float(src['m3_yr'][src['kind'] == k].sum() / 1e9), 4) for k in range(3)}
    total = sum(by.values())
    return dict(km3_yr=by, share={k: round(v / total, 3) for k, v in by.items()}, ends={names[k]: int(((src['kind'] == k) & (src['m3_yr'] > 0)).sum()) for k in range(3)},
                the_loop_gathers_km3_yr=round(by['lake deltas'] + by['closed hollows'], 4), land_sheds_km3_yr=round(src['land_yield_m3_yr'] / 1e9, 4),
                unreached_km3_yr=round(src['unreached_m3_yr'] / 1e9, 5),
                basalt_share={names[k]: round(float((src['m3_yr'] * src['basalt'])[src['kind'] == k].sum() / max(src['m3_yr'][src['kind'] == k].sum(), 1e-9)), 3) for k in range(3)},
                sand_study_capped_rule=dict(ends_km3_yr=round(float(src['capped_m3_yr'].sum() / 1e9), 4),
                                            wind_offer_central=wind_offer(F, src, 'central'), wind_offer_opensand=wind_offer(F, src, 'opensand'),
                                            note='the sand study credits each end with up to 5,000 km2 of the ground upstream, so chains of lakes count the same ground more than once; '
                                                 'the wind\'s offer compares the rivers\' sand it leaves the wind under each rule'),
                note='durable sand, bulk volume, each land cell\'s sand counted once at the first standing water or dry hollow floor below it')


def haul(src, lat0, lon0, need_km3, years=DELIVER_YEARS, radii_km=RADII_KM):
    """Sand for one sea from the loop's ends nearest it (lake deltas and closed hollows). Within each radius: the sand
    a year, the years to gather `need_km3`, the mean distance, the haulage power while it lasts and the basalt's share
    of the sand. And for each span of years, the same for the nearest ends that gather the need within it (None where
    all the Moon's ends together cannot)."""
    use = src['kind'] < 2
    d = distance_km(src['lat'][use], src['lon'][use], lat0, lon0)
    order = np.argsort(d)
    d, rate, basalt = d[order], src['m3_yr'][use][order], src['basalt'][use][order]
    cum, cum_d, cum_b = np.cumsum(rate), np.cumsum(rate * d), np.cumsum(rate * basalt)
    tonnes = S.BULK_DENSITY / 1000.0                              # tonnes in a bulk cubic metre

    def within(k):
        if k < 0 or cum[k] <= 0:
            return None
        mean_km = cum_d[k] / cum[k]
        return dict(radius_km=round(float(d[k])), sand_km3_yr=round(float(cum[k] / 1e9), 5), years=round(float(need_km3 * 1e9 / cum[k])), mean_haul_km=round(float(mean_km)),
                    power_gw=round(float(cum[k] * tonnes * mean_km * HAUL_KJ_PER_T_KM * 1e3 / S.SECONDS_PER_YEAR / 1e9), 3),
                    energy_j=float(f'{need_km3 * 1e9 * tonnes * mean_km * HAUL_KJ_PER_T_KM * 1e3:.3g}'), basalt_share=round(float(cum_b[k] / cum[k]), 2),
                    share_of_the_loops_sand=round(float(cum[k] / cum[-1]), 3))

    by_years = {}
    for span in years:
        k = int(np.searchsorted(cum, need_km3 * 1e9 / span))
        by_years[str(int(span))] = within(k) if k < cum.size else None
    return dict(need_km3=round(need_km3, 1), by_radius=[within(int(np.searchsorted(d, r, side='right')) - 1) for r in radii_km], within_years=by_years)


# ---- the seas -----------------------------------------------------------------------------------------------------
def window(seed, reach_deg):
    """Rows and columns of a box of reach_deg each way round a seed (columns wrap)."""
    r, c = cell_of(*seed)
    half_r = int(reach_deg * PPD)
    half_c = int(min(reach_deg * PPD / max(math.cos(math.radians(seed[0])), 0.1), COLS // 2 - 1))
    rows = np.arange(max(r - half_r, 0), min(r + half_r + 1, ROWS))
    cols = np.arange(c - half_c, c + half_c + 1) % COLS
    return rows, cols, (r - rows[0], half_c)


def outline(seed, allowed, score, area_row, target_km2):
    """The joined ground round the seed whose score is highest, taken down the score until it holds target_km2: the
    cells with a score of t or more that join the seed through such cells, for the t that brings their area nearest
    the target. `allowed`, `score` and `area_row` are a window's; returns the mask and t."""
    from scipy import ndimage
    assert allowed[seed], 'the seed lies on ground the sea may not cover'
    area = np.broadcast_to(area_row[:, None], allowed.shape)
    lo, hi = float(score[allowed].min()), float(score[seed])
    best = None
    for _ in range(40):
        t = 0.5 * (lo + hi)
        labels, _n = ndimage.label(allowed & (score >= t))
        mask = labels == labels[seed]
        km2 = float(area[mask].sum())
        if best is None or abs(km2 - target_km2) < abs(best[2] - target_km2):
            best = (mask, t, km2)
        if km2 > target_km2:
            lo = t
        else:
            hi = t
        if hi - lo < 1e-6 * max(abs(hi), 1.0):
            break
    mask = ndimage.binary_fill_holes(best[0]) & allowed            # small islands of lower score inside the sea are sand too
    rounded = (ndimage.gaussian_filter(mask.astype(np.float32), OUTLINE_SIGMA) > 0.5) & allowed   # the outline rounded off over a few cells
    labels, _n = ndimage.label(rounded)
    if labels[seed]:
        mask = ndimage.binary_fill_holes(labels == labels[seed]) & allowed
    return mask, best[1]


def taper(mask, towards_deg=None, aspect=1.0):
    """A sea's thickness as a share of its full thickness: TAPER_EDGE at its outermost cells, rising to 1 over
    TAPER_CELLS cells inward. With the direction its net drift blows towards, the margin is narrow on the upwind side,
    where the wind strips a clean edge, and wide downwind, where sand laps out over the ground: TAPER_SIDES of
    TAPER_CELLS at the two. `aspect` is a column's width against a row's height (the cosine of the latitude), so that
    the widths are the same on the ground whichever way a margin faces."""
    from scipy import ndimage
    inward, nearest = ndimage.distance_transform_edt(mask, sampling=(1.0, aspect), return_indices=True)
    width = TAPER_CELLS
    if towards_deg is not None:
        rows, cols = np.indices(mask.shape)
        north, east = -(nearest[0] - rows).astype(np.float32), (nearest[1] - cols).astype(np.float32) * aspect   # from each cell to the nearest ground off the sea
        t = math.radians(towards_deg)
        downwind = (east * math.sin(t) + north * math.cos(t)) / np.maximum(np.hypot(east, north), 1e-6)      # 1 where that ground lies down the drift
        width = TAPER_CELLS * (TAPER_SIDES[0] + (TAPER_SIDES[1] - TAPER_SIDES[0]) * 0.5 * (1.0 + downwind))
    return np.where(mask, np.clip(TAPER_EDGE + (1.0 - TAPER_EDGE) * (inward - 1.0) / width, TAPER_EDGE, 1.0), 0.0).astype(np.float32)


def margin_loss(mask_full, thick, q_eff, az_deg, water, lat_f, area_row, years):
    """What a sea's downwind margin sheds in `years`: every cell of the sea whose net drift leads onto ground outside
    it passes sand on at the drift the cover and the shelter there allow, across its own width, for as long as it has
    sand to give. Returns the thin sand set down on the cells that receive it (m), the body's thickness less what
    left (m), the volume shed and the part of it that goes straight into water (km3)."""
    idx = np.flatnonzero(mask_full.reshape(-1))
    rows, cols = idx // COLS, idx % COLS
    dr_, dc, _stride = S.grid_step(az_deg.reshape(-1)[idx], lat_f[rows], 1)
    r2 = np.clip(rows + dr_, 0, ROWS - 1)
    target = r2 * COLS + (cols + dc) % COLS
    leaves = ~mask_full.reshape(-1)[target] & np.isfinite(az_deg.reshape(-1)[idx])
    area_m2 = area_row * 1e6
    width = np.sqrt(area_m2[rows])                                 # a cell's width across the drift, as the square root of its area
    volume = np.minimum(q_eff.reshape(-1)[idx] * width * years, thick.reshape(-1)[idx] * area_m2[rows]) * leaves
    left = thick.copy().reshape(-1)
    left[idx] -= (volume / area_m2[rows]).astype(np.float32)
    wet = water.reshape(-1)[target]
    thin = np.zeros(ROWS * COLS, np.float32)
    np.add.at(thin, target[leaves & ~wet], (volume / area_m2[r2])[leaves & ~wet].astype(np.float32))
    return thin.reshape(ROWS, COLS), left.reshape(ROWS, COLS), float(volume.sum() / 1e9), float(volume[wet].sum() / 1e9)


def by_side(src, lat0, lon0, radius_km):
    """The basalt's share of the loop's sand within radius_km of a place, by the side it lies on: SECTORS directions
    clockwise from north, each with its sand (bulk m3 a year). A side with no source takes the whole circle's share."""
    use = (src['kind'] < 2) & (distance_km(src['lat'], src['lon'], lat0, lon0) <= radius_km)
    east = (np.asarray(src['lon'][use]) - lon0 + 180.0) % 360.0 - 180.0
    bearing = np.degrees(np.arctan2(east * math.cos(math.radians(lat0)), np.asarray(src['lat'][use]) - lat0)) % 360.0
    sector = np.floor(bearing / (360.0 / SECTORS) + 0.5).astype(int) % SECTORS
    sand = np.bincount(sector, weights=src['m3_yr'][use], minlength=SECTORS)
    dark = np.bincount(sector, weights=(src['m3_yr'] * src['basalt'])[use], minlength=SECTORS)
    whole = float(dark.sum() / max(sand.sum(), 1e-9))
    return np.where(sand > 0, dark / np.maximum(sand, 1e-9), whole), sand, whole


def make_up(mask, lat_f, lon_f, lat0, lon0, shares, whole):
    """The basalt's share of a sea's sand cell by cell: at its middle the share of all its sand, and towards each edge
    the share of the sand that comes from that side, since each source's sand is set down on the side nearest it."""
    rows, cols = np.nonzero(mask)
    north = lat_f[rows] - lat0
    east = ((lon_f[cols] - lon0 + 180.0) % 360.0 - 180.0) * math.cos(math.radians(lat0))
    at = (np.degrees(np.arctan2(east, north)) % 360.0) / (360.0 / SECTORS)
    low = np.floor(at).astype(int)
    side = shares[low % SECTORS] * (1.0 - (at - low)) + shares[(low + 1) % SECTORS] * (at - low)
    out = np.hypot(north, east)
    out = out / max(float(out.max()), 1e-9)
    field = np.zeros(mask.shape, np.float32)
    field[rows, cols] = whole + (side - whole) * out
    return field


def dunes_of(net_over_total, drift_m2_yr, moving, thickness_m, years=YEARS):
    """The dunes the sand study's rules give a body of sand: the form by the net drift's share of the total, the crest
    spacing as the smaller of what the drift can gather in the years and three dune heights of sand, and the height."""
    form = 'crescent dunes and transverse ridges' if net_over_total >= 0.8 else 'linear dunes' if net_over_total >= 0.3 else 'star dunes'
    grown = math.sqrt(2.0 * S.CREST_ASPECT * max(drift_m2_yr * moving, 1e-6) * years)
    spacing = min(grown, 3.0 * S.CREST_ASPECT * thickness_m, S.MIXED_LAYER_M[0])
    return dict(form=form, crest_spacing_m=round(spacing), height_m=round(spacing / S.CREST_ASPECT),
                limited_by='the sand at hand' if spacing < grown else 'the drift and the years')


def lay_out(F, conv, rough, log=print):
    """Each sea's outline and sand on the grid, by its rule. Returns, per sea, the full-grid mask, the thickness (m)
    and the figures of its ground."""
    from scipy import ndimage
    surface = np.load(S.SURFACE, allow_pickle=False)
    channel = np.zeros(ROWS * COLS, bool)
    channel[surface['channel_cell'][surface['channel_discharge_m3s'] >= CHANNEL_M3S]] = True
    channel = channel.reshape(ROWS, COLS)
    out = {}
    for name, sea in SEAS.items():
        rows, cols, seed = window(sea['seed'], sea['reach_deg'])
        w = np.ix_(rows, cols)
        water, aridity, height = F.water[w], F.aridity[w], F.height[w].astype(np.float32)
        near_water = ndimage.binary_dilation(water, iterations=WATER_BUFFER_CELLS)
        near_channel = ndimage.binary_dilation(channel[w], iterations=CHANNEL_BUFFER_CELLS)
        allowed = ~water & ~near_water & ~near_channel & (aridity < HYPER_ARID)
        if sea['relief_m'] is not None:
            allowed &= rough[w] < sea['relief_m']
        proud = height - ndimage.gaussian_filter(height, SURFACE_SIGMA)   # how far the ground stands above the ground round it (m)
        if sea['rule'] == 'convergence':
            near_edge = np.clip((aridity / HYPER_ARID - EDGE_FROM) / (1.0 - EDGE_FROM), 0.0, 1.0)
            score = conv[w] - sea['swell_m_per_kyr'] * proud / sea['thickness_m'] - sea['edge_m_per_kyr'] * near_edge * near_edge * (3.0 - 2.0 * near_edge)
        else:
            score = -height
        mask, t = outline(seed, allowed, score, F.area_row[rows], sea['area_km2'])
        full = np.zeros((ROWS, COLS), bool)
        full[w] = mask
        drift = F.fine['central']
        a = np.radians(drift['az'][w][mask])
        towards = math.degrees(math.atan2(float((drift['rdp'][w][mask] * np.sin(a)).sum()), float((drift['rdp'][w][mask] * np.cos(a)).sum()))) % 360.0
        edge = taper(mask, towards, math.cos(math.radians(sea['seed'][0])))
        even = edge * sea['thickness_m']                          # the sand as an even sheet, and as it lies: deep in the hollows, thin on the swells
        lying = edge * np.clip(sea['thickness_m'] - proud, 0.0, THICKNESS_MOST * sea['thickness_m'])
        area = np.broadcast_to(F.area_row[rows][:, None], mask.shape)
        thick = np.zeros((ROWS, COLS), np.float32)
        thick[w] = lying * float((even * area).sum() / max(float((lying * area).sum()), 1e-9))
        out[name] = dict(mask=full, thick=thick, threshold=float(t))
        log(f'  {name}: outline of {float(F.a2[full].sum()):,.0f} km2 at a score of {t:.3f}, {float((thick * F.a2).sum() / 1e3):.0f} km3')
    return out


def ground_figures(F, sea, body, conv, rough, sites, regions, districts):
    """What one sea's ground is: its place, size, winds, relief, what it covers and what lies near it."""
    from scipy import ndimage
    m = body['mask']
    w = F.a2[m].astype(np.float64)
    f = F.fine['central']
    la, lo = np.broadcast_to(F.lat_f[:, None], m.shape)[m], np.broadcast_to(F.lon_f[None, :], m.shape)[m]
    x, y, z = (np.cos(np.radians(la)) * np.cos(np.radians(lo)) * w).sum(), (np.cos(np.radians(la)) * np.sin(np.radians(lo)) * w).sum(), (np.sin(np.radians(la)) * w).sum()
    lat0, lon0 = float(np.degrees(np.arctan2(z, np.hypot(x, y)))), float(np.degrees(np.arctan2(y, x)) % 360.0)
    mean = lambda field: float((field[m] * w).sum() / w.sum())    # noqa: E731
    a = np.radians(f['az'][m])
    east, north = (f['rdp'][m] * np.sin(a) * w).sum() / w.sum(), (f['rdp'][m] * np.cos(a) * w).sum() / w.sum()
    local = ndimage.maximum_filter(F.height, 3, mode='nearest')[m] - ndimage.minimum_filter(F.height, 3, mode='nearest')[m]   # relief within a 5.7 km square
    h = sea['thickness_m']
    extent = distance_km(la, lo, lat0, lon0)
    near = sorted(((round(float(distance_km(np.array([s['lat']]), np.array([s['lon'] % 360.0]), lat0, lon0)[0])), s['name']) for s in sites))[:3]
    reg = np.bincount(regions[m], weights=w)
    dis = np.bincount(districts[m], weights=w)
    return dict(label=sea['label'], rule=sea['rule'], seed_today=list(sea['seed']), centre_today=[round(lat0, 2), round(lon0, 2)],
                area_km2=round(float(w.sum())), sand_km3=round(float((body['thick'] * F.a2).sum() / 1e3), 1), thickness_m=h,
                reaches_km=round(float(extent.max())), latitudes=[round(float(la.min()), 2), round(float(la.max()), 2)],
                longitudes_east=[round(float(lo.min()), 2), round(float(lo.max()), 2)],
                height_m=[round(float(np.percentile(F.height[m], q))) for q in (10, 50, 90)], aridity=round(mean(F.aridity), 3),
                plant_cover=round(mean(F.cover), 3), drift_potential_m2_yr=round(mean(f['dp']), 1), net_m2_yr=round(float(np.hypot(east, north)), 1),
                net_towards_deg=round(float(np.degrees(np.arctan2(east, north)) % 360.0)),
                net_over_total=round(float((f['rdp'][m] * w).sum() / max((f['dp'][m] * w).sum(), 1e-9)), 2),
                convergence_m_per_kyr=round(mean(conv), 3), moving_share=round(mean(F.moving), 2), sheltered=round(mean(F.trap_index), 2),
                basalt_ground_share=round(mean(F.mafic), 2), level_share=round(float(w[rough[m] < LEVEL_M].sum() / w.sum()), 2),
                relief_in_5km=dict(under_the_sand=round(float(w[local <= h].sum() / w.sum()), 2), ridden_over=round(float(w[(local > h) & (local <= 3 * h)].sum() / w.sum()), 2),
                                   stands_through=round(float(w[local > 3 * h].sum() / w.sum()), 2)),
                nearest_heritage_sites_km=[dict(km=k, site=n) for k, n in near],
                regions={str(int(k)): round(float(v)) for k, v in enumerate(reg) if v > 0}, districts={str(int(k)): round(float(v)) for k, v in enumerate(dis) if v > 0})


def build(log=print):
    """The whole reckoning: the fields, the patches, the sources, the seas and the case of the map with the seas on it.
    Returns the record and the arrays for the product."""
    t0 = time.time()
    F = S.fields(log)
    f = F.fine['central']
    land = ~F.water
    conv = convergence(f['rdp'], f['az'], land)
    rough = spread(F.height)
    hyper = land & (F.aridity < HYPER_ARID)
    found = patches(hyper, F, conv, rough)
    log(f'  {len(found)} hyper-arid patches, {float(F.a2[hyper].sum()):,.0f} km2 ({time.time() - t0:.0f} s)')
    src = sources(F)
    totals = supply_totals(F, src)
    log(f'  sources: {totals["km3_yr"]} km3 a year; the land sheds {totals["land_sheds_km3_yr"]}; the capped rule {totals["sand_study_capped_rule"]}')
    bodies = lay_out(F, conv, rough, log)
    sites = json.loads(HERITAGE.read_text())['sites']
    units = np.load(S.REGIONS, allow_pickle=False)
    regions = np.repeat(np.repeat(units['region'], 2, axis=0), 2, axis=1).astype(np.int64)
    districts = np.repeat(np.repeat(units['district'], 2, axis=0), 2, axis=1).astype(np.int64)
    record = json.loads(S.REGION_RECORD.read_text())
    region_code = {k + 1: r['code'] for k, r in enumerate(record['regions'])}
    district_code = {k + 1: d['code'] for k, d in enumerate(record['districts'] + record['coastal_seas'])}

    # the margins' losses: each sea sheds sand downwind where its net drift leads off it, with the sand left open
    q_open = (f['rdp'] * F.moving * F.exposure).astype(np.float32)
    thin = np.zeros((ROWS, COLS), np.float32)
    thick = np.zeros((ROWS, COLS), np.float32)
    made = np.zeros((ROWS, COLS), np.uint8)
    seas = {}
    for k, (name, sea) in enumerate(SEAS.items(), start=1):
        body = bodies[name]
        g = ground_figures(F, sea, body, conv, rough, sites, regions, districts)
        g['regions'] = {region_code.get(int(i), '0'): v for i, v in g['regions'].items() if int(i) in region_code}
        g['districts'] = {district_code.get(int(i), '0'): v for i, v in g['districts'].items() if int(i) in district_code}
        shed, left, shed_km3, wet_km3 = margin_loss(body['mask'], body['thick'], q_open, f['az'], F.water, F.lat_f, F.area_row, YEARS)
        thin += shed
        thick += left
        made[body['mask']] = k
        g['margin'] = dict(shed_km3=round(shed_km3, 2), straight_into_water_km3=round(wet_km3, 3), share_of_the_sea=round(shed_km3 / max(g['sand_km3'], 1e-9), 4),
                           note='what the downwind margin passes on in ten thousand years at the net drift there, the sand left open; where it goes is the map\'s to say')
        g['dunes'] = dunes_of(g['net_over_total'], g['net_m2_yr'], g['moving_share'], sea['thickness_m'])
        g['haul'] = haul(src, g['centre_today'][0], g['centre_today'][1], g['sand_km3'])
        g['plan'] = dict(deliver_years=sea['deliver_years'], **(g['haul']['within_years'][str(int(sea['deliver_years']))] or {}))
        per_km2 = g['sand_km3'] / g['area_km2']
        g['other_sizes'] = [dict(area_km2=a, sand_km3=round(a * per_km2, 1), **{f'in_{int(y)}_years': v for y, v in
                                 haul(src, g['centre_today'][0], g['centre_today'][1], a * per_km2)['within_years'].items()}) for a in SIZES_KM2]
        g['id'] = k
        seas[name] = g
        log(f'  {name}: {g["area_km2"]:,} km2, {g["sand_km3"]} km3, sheds {shed_km3:.2f} km3 ({time.time() - t0:.0f} s)')

    # the map: the sand study's open-sand case on the ends' own ground, without the seas and with them set down on it
    n = ROWS * COLS
    rivers = np.zeros(n, np.float32)
    rivers[src['cell']] = src['m3_yr'] / (F.area_row[src['cell'] // COLS] * 1e6)
    rivers = rivers.reshape(ROWS, COLS)
    rule = 'each end of drainage given the sand of its own ground, down to the first standing water'
    held, bare = {}, {}
    var = dict(S.VARIANTS['opensand'], rests_on=f'the sand study\'s open-sand case with {rule}')
    open_sums = S.variant('open', var, F, lambda key, array: bare.__setitem__(key, array) if key in ('sand_share', 'sand_thickness_m') else None,
                          suffix='', rivers=rivers, log=log)
    var = dict(var, rests_on=var['rests_on'] + ', and the made seas set down on it: each sea\'s body lies where it was put, and what its downwind margin sheds in ten '
               'thousand years is carried as the land\'s own sand is')
    summary = S.variant('seas', var, F, lambda key, array: held.__setitem__(key, array), placed=dict(thin=thin, thick=thick), suffix='', rivers=rivers, log=log)
    open_share, open_thick = bare['sand_share'], bare['sand_thickness_m']
    for name, g in seas.items():
        m = bodies[name]['mask']
        lat0, lon0 = g['centre_today']
        rows, cols, _ = window((lat0, lon0), 20.0)
        w = np.ix_(rows, cols)
        area = F.a2[w]
        inside = m[w]
        share, form, spacing = held['sand_share'][w], held['form'][w], held['crest_spacing_m'][w]
        gained = np.maximum(share - open_share[w], 0.0) * ~inside
        dune = inside & (form >= 2)
        far = gained > 0.02
        by_form = {S.FORM_NAMES[i]: round(float(area[inside & (form == i)].sum())) for i in range(1, 6) if (inside & (form == i)).any()}
        g['as_mapped'] = dict(sand_cover_inside=round(float((share * area)[inside].sum() / area[inside].sum()), 3), forms_km2=by_form,
                              crest_spacing_m={k: round(float(np.percentile(spacing[dune], q))) if dune.any() else None for k, q in (('p10', 10), ('median', 50), ('p90', 90))},
                              plume_km2=round(float((gained * area).sum())),
                              plume_reaches_km=round(float(distance_km(np.broadcast_to(F.lat_f[rows][:, None], far.shape)[far], np.broadcast_to(F.lon_f[cols][None, :], far.shape)[far],
                                                                       lat0, lon0).max())) if far.any() else 0,
                              plume_sand_km3=round(float((np.maximum(held['sand_thickness_m'][w] - open_thick[w], 0.0) * ~inside * area).sum() / 1e3), 2))
    # which of each cell's sand is a made sea's, and that sand's basalt share by its plan's sources
    from scipy import ndimage
    part = np.where(held['sand_share'] > 0, np.clip((held['sand_share'] - open_share) / np.maximum(held['sand_share'], 1e-6), 0.0, 1.0), 0.0).astype(np.float32)
    part[made > 0] = 1.0
    nearest = ndimage.distance_transform_edt(made == 0, return_distances=False, return_indices=True)
    owner = made[nearest[0], nearest[1]]
    inside = np.zeros((ROWS, COLS), np.float32)
    for name, g in seas.items():
        shares, sand, whole = by_side(src, g['centre_today'][0], g['centre_today'][1], g['plan']['radius_km'])
        inside += make_up(bodies[name]['mask'], F.lat_f, F.lon_f, g['centre_today'][0], g['centre_today'][1], shares, whole)
        g['plan']['basalt_share_by_side'] = dict(clockwise_from_north=[round(float(v), 2) for v in shares],
                                                 sand_km3_yr=[round(float(v) / 1e9, 5) for v in sand])
    dark = np.where(part > 0, inside[nearest[0], nearest[1]], 0.0).astype(np.float32)   # a sea's fringe and plume take its nearest cell's sand
    del nearest, owner, inside
    height = np.where(held['form'] >= 2, held['crest_spacing_m'] / S.CREST_ASPECT, 0.0).astype(np.float32)
    arrays = dict(held, made_sea=made, dune_height_m=height, made_sand_part=part, made_sand_basalt=dark, sand_share__open=open_share, sand_thickness_m__open=open_thick)
    record = dict(years=YEARS, hyper_arid=dict(km2=round(float(F.a2[hyper].sum())), share_of_land=round(float(F.a2[hyper].sum() / F.land_km2), 4), patches=found),
                  supply=totals, seas=seas, map=dict(open=open_sums, seas=summary), land_km2=round(F.land_km2), seconds=round(time.time() - t0))
    return record, arrays, F


EVIDENCE = ('A siting reckoning on the sand study\'s fields, nothing in it observed. The winds, the cover, the supply by catchment and the carrying scheme are that '
            'study\'s and carry its limits (sand.json). The outlines follow stated rules on the 1.9 km grid; the sand each sea takes is its outline times a chosen '
            'thickness; the haulage is the terrain study\'s rail-freight energy over great-circle distances from the nearest ends of drainage; the dunes follow '
            'the sand study\'s rules of form and spacing, fitted on Earth.')
READING_RULE = ('Grids are 2,880 rows from 90 N by 5,760 columns from 0 E, 16 per degree, on today\'s Moon. The fields sand_share, active_share, form, '
                'crest_azimuth_deg, crest_spacing_m, drift_azimuth_deg and sand_thickness_m are the sand product\'s for one case, the open-sand case with the made seas '
                'set down on it and each end of drainage given the sand of its own ground, and are read as that product\'s central fields are. made_sea is 0 off the seas '
                'and each sea\'s number inside its outline (metadata.seas); dune_height_m is the crest spacing over the sand study\'s aspect of 20 where dunes stand; '
                'made_sand_part is the share of a cell\'s sand that a made sea brought (1 inside a sea, the gain over the same case without the seas outside) and '
                'made_sand_basalt the basalt\'s share of that sand by the sea\'s haulage plan; sand_share__open and sand_thickness_m__open are the same case without the seas.')


def write(log=print):
    BULK.mkdir(parents=True, exist_ok=True)
    record, arrays, _F = build(log)
    partial = PRODUCT.with_suffix('.partial.npz')
    writer = S.NpzWriter(partial)
    for key, array in arrays.items():
        writer.add(key, array)
    inputs = {str(p.relative_to(ROOT)): S.sha256(p) for p in (S.OUT, S.MOISTURE, HERITAGE)}
    ignored = {str(p): S.sha256(p) for p in (S.PRODUCT, S.SERIES, S.SURFACE, S.ROUTING, S.ROCK, S.REGIONS, S.REGION_RECORD)}
    producer = dict(files={c: S.sha256(ROOT / c) for c in CODE}, constants=dict(MOON_RADIUS=MOON_RADIUS))
    meta = dict(schema=SCHEMA, grid=dict(pixels_per_degree=PPD, rows=ROWS, cols=COLS), producer=producer, inputs=inputs, ignored_inputs=ignored,
                evidence=EVIDENCE, reading_rule=READING_RULE, forms=list(S.FORM_NAMES),
                seas=[dict(id=g['id'], name=name, label=g['label'], centre_today=g['centre_today'], area_km2=g['area_km2'], thickness_m=g['thickness_m'],
                           form=g['dunes']['form'], crest_spacing_m=g['dunes']['crest_spacing_m'], height_m=g['dunes']['height_m'],
                           net_towards_deg=g['net_towards_deg'], basalt_share=g['plan'].get('basalt_share', 0.0)) for name, g in record['seas'].items()])
    writer.add('metadata', np.array(json.dumps(meta)))
    writer.close()
    partial.replace(PRODUCT)
    out = dict(schema=SCHEMA, producer=producer, inputs=inputs, ignored_inputs=ignored, evidence=EVIDENCE, reading_rule=READING_RULE,
               parameters=dict(seas=SEAS, hyper_arid=HYPER_ARID, haul_kj_per_t_km=HAUL_KJ_PER_T_KM, loop_range_gw=LOOP_RANGE_GW, level_m=LEVEL_M,
                               level_window_cells=LEVEL_WINDOW, taper_cells=TAPER_CELLS, taper_edge=TAPER_EDGE, taper_sides=TAPER_SIDES, water_buffer_cells=WATER_BUFFER_CELLS,
                               channel_buffer_cells=CHANNEL_BUFFER_CELLS, channel_m3s=CHANNEL_M3S, convergence_sigma_cells=CONVERGENCE_SIGMA, outline_sigma_cells=OUTLINE_SIGMA,
                               bulk_density_kg_m3=S.BULK_DENSITY, supply=S.SUPPLY['central'], sizes_km2=SIZES_KM2, deliver_years=DELIVER_YEARS, radii_km=RADII_KM),
               **record, product=dict(path=str(PRODUCT), sha256=S.sha256(PRODUCT)),
               limits=['Everything the sand study lists as a limit holds here: one run\'s winds at 170 km and 3 days, the display cover curve, a capped catchment behind each end of drainage, Earth\'s dune rules.',
                       'A sea is set down whole and stands ten thousand years; the years of its building are not followed.',
                       'A sea\'s body stays where it is put: its own slow march downwind (the net drift over its thickness, a few kilometres to some tens in ten thousand years) is not drawn.',
                       'The margin sheds along the net drift only; the to-and-fro of winds from many sides, which fringes a sea with a sheet a few kilometres wide, enters as the map\'s last smoothing over 4.7 km.',
                       'Haulage is priced over great-circle distances at one energy per tonne-kilometre; no route, grade or carrier is designed.',
                       'The sand hauled is the durable half of the sand-sized grains; sorting it from the dredged sediment is not priced.'])
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    log(f'  wrote {OUT} and {PRODUCT} ({PRODUCT.stat().st_size / 1e6:.0f} MB)')
    return out


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.split('\n')[0]).parse_args(argv)
    record = write()
    for name, g in record['seas'].items():
        print(f'{name}: {g["area_km2"]:,} km2, {g["sand_km3"]} km3, {g["dunes"]["form"]} {g["dunes"]["crest_spacing_m"]} m apart; sheds {g["margin"]["shed_km3"]} km3')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
