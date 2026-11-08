"""The monthly tide in the shores' groundwater and where the sea's load settles, on the coast-water study's belts.

    python -m research.studies.coast_water.tide_terrain --download    # once, about 95 MB of LOLA rows
    python -m research.studies.coast_water.tide_run                   # about ten minutes, 1.5 GB

It recomputes the coast-water study's 301 K basis with its own functions (the deserts, the plants' uptake, the 16
px/deg shores and belts), gives every shore the monthly tide of its water (geography's tide grid for the seas; each
lake's own equilibrium tide as a closed basin of fixed volume, from the same forcing), and measures the belts with the
tide against the same belts without it. The tide's rise of the mean water table comes from the numerical cross-shore
model (tide.py), run for each of the coast-water study's cases over a grid of shore slope, tide amplitude and draw,
and checked against Ferris and Nielsen and against direct runs off the grid. Profiles at 100 m on LOLA's 256 px/deg rows across desert sea coasts at 45-60 degrees run the
cross-shore model on the real ground; load.py gives the load's minerals and rates. Writes results/coast_tide.json.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from research.studies.coast_water import belts as B
from research.studies.coast_water import load as L
from research.studies.coast_water import run as R
from research.studies.coast_water import shores as S
from research.studies.coast_water import tide as TD
from research.studies.coast_water import tide_terrain as TT
from shared.constants import JULIAN_DAY, MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'coast_tide.json'
SCHEMA = 'terluna.research.coast-tide/1'
CODE = ('research/studies/coast_water/tide_run.py', 'research/studies/coast_water/tide.py',
        'research/studies/coast_water/load.py', 'research/studies/coast_water/tide_terrain.py',
        'research/studies/coast_water/run.py', 'research/studies/coast_water/belts.py',
        'research/studies/coast_water/shores.py', 'research/studies/coast_water/extract.py',
        'research/studies/climate_targets/extract.py', 'research/studies/carbon_cycle/run.py',
        'research/studies/sea_chemistry/run.py', 'geography/tides.py', 'geography/coast_terrain.py',
        'geography/topography.py')
TIDES = 'geography/results/tides.json'
TIDE_GRID = 'geography/products/tides_28pct_4ppd.npz'
CHEMISTRY = 'research/studies/sea_chemistry/results/sea_chemistry.json'
COAST_WATER = 'research/studies/coast_water/results/coast_water.json'
TERRAIN_RECORD = 'research/studies/coast_water/tide_inputs.json'
INPUTS = {TIDES: 'terluna.geography.tides/1', CHEMISTRY: 'terluna.research.sea-chemistry/1',
          COAST_WATER: 'terluna.research.coast-water/1', TERRAIN_RECORD: 'terluna.research.coast-tide-inputs/1'}
LAKE_INPUTS = (S.DRAINAGE,) + S.LOLA
CACHE = 'research/runs/coast_water/tide/numerics.json'
YEAR_S = JULIAN_DAY * 365.25
GROUNDS = ('slow', 'central', 'sandy')
SEA_CLASSES = (0, 1)                      # sea and island coasts in SHORE_CLASSES
DESERT = 'rain_mask.every_cell'
EPS_TABLE = (0.01, 0.03, 0.1, 0.3, 0.6, 1.0, 2.0, 4.0, 10.0, 20.0)
RISE_CASES = ('slow_15m', 'slow_25m', 'central_15m', 'central_25m', 'sandy_15m', 'sandy_25m')
RISE_SLOPES = (0.006, 0.02, 0.07, 0.25)          # the sea coasts' 16 px/deg slopes run 0.005-0.40 (5th-95th)
RISE_AMPLITUDES = (0.25, 0.8, 1.6, 2.8)          # half their typical monthly ranges, 0.21-2.84 m
RISE_DRAWS_MM_DAY = (1.2, 2.0, 3.0)              # the uptake, 1.36-2.90 mm a day on the desert sea coasts
DISTANCES_M = (0, 25, 50, 100, 200, 400)
EVIDENCE = (
    'A runnable calculation on the coast-water study\'s own basis (its 301 K deserts, uptake, 16 px/deg shores and '
    'belts, recomputed with its functions), geography\'s tide grid and its forcing, the sea chemistry\'s held seas and '
    'LOLA\'s 256 px/deg rows. The tide\'s reach is linear theory (Ferris 1951; Nielsen 1990 for the sloping face) '
    'checked against a numerical cross-shore model (unconfined Boussinesq, implicit, the sea held at its level on '
    'the slope, a seepage face, the plants\' draw within the reach), and the belts with the tide take the model\'s '
    'rise of the mean water table, run for each case over a grid of slope, amplitude and draw and checked against '
    'direct runs off the grid. The lakes\' tides are '
    'each lake\'s equilibrium tide as a closed basin, without its own water\'s attraction. Specific yields are '
    'published ranges for basalt and sand; transmissivity is the coast-water study\'s. Beaches, sediment, waves, '
    'the capillary fringe\'s effect on the yield near the ground and the density of the fringe\'s water are left '
    'out. The load is the evaporation path of the sea chemistry\'s waters with its constants and a mass balance of '
    'the draw; no field measurement checks any of it.')
READING_RULE = (
    'Lengths in m (belts and distances) or km (shore), areas in km2, heights in m above the mean water level, tide '
    'amplitudes A are half the typical monthly range (the median over anomalistic months). eps = kA/slope with '
    'k = sqrt(w Sy/2T) for the 27.55-day line and c = E/2T. '
    'rise is the height the coast-water rule s W + c W^2 = R needs added to R to give the model\'s belt; a belt '
    '"always" within reach is where the water table at its lowest in the cycle stays within the reach. Cases are '
    'named by ground (slow, central, sandy crust as the coast-water study names them) and reach; specific yields '
    'are each ground\'s central value unless named. Deserts are the coast-water study\'s rain_mask.every_cell '
    'unless named. Deposits are mm per century of solid mineral, and of soil fully cemented at 30% porosity.')


def digest(path) -> str:
    h = hashlib.sha256()
    with open(ROOT / path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def central_yield(ground):
    lo, hi = TD.SPECIFIC_YIELD[ground]
    return math.sqrt(lo * hi)


# ----------------------------------------------------------------------------------------- the coast-water basis
def basis():
    data, sidecar = R.load()
    base = R.basis(data[R.TARGETS])
    f = base['share_of_the_way']
    g = np.load(ROOT / R.GCM)
    sh = np.load(ROOT / R.SHORES)
    grid = R.t21(g)
    shore_cells = {k: sh[f'cell_{k}'] for k in ('land', 'island_land') + tuple(S.WATER_CLASSES.values())}
    rain, masks, _ = R.deserts(g, grid, f, base, shore_cells, data[R.TARGETS])
    rings = {k: data[p] for k, p in R.RINGS.items()}
    ratios = R.ring_ratios(rings)
    vector = R.vector_speed_factor(g, sidecar['runs']['eddies']['sun_window_rad'], grid)
    demand = {}
    for key in ('design', 'test'):
        s = R.surface(g, key, sidecar, grid, ratios, vector)
        demand[key] = R.demand_by_cell(s)
        del s
    pet = {k: R.interpolate(demand['design'][k], demand['test'][k], f) for k in ('all', 'day', 'night')}
    uptake = np.maximum(pet['all'] - rain, 0.0) / 1000.0 / R.DAY_S
    table = data[R.TABLE]
    return dict(data=data, sidecar=sidecar, base=base, sh=sh, grid=grid, shore_cells=shore_cells, rain=rain,
                masks=masks, pet=pet, uptake=uptake, cases=R.cases(table))


# ------------------------------------------------------------------------------------------- the tides on shores
def unit(lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)), axis=-1)


def pixel_latlon(index):
    r, c = np.divmod(np.asarray(index, np.int64), S.COLS)
    return 90.0 - (r + 0.5) / S.PPD, (c + 0.5) / S.PPD


def tide_forcing():
    """The tide lines' forcing in geography's monomials, recovered from its grid (each node's line amplitudes are
    the forcing times the node's response basis), and the ratio of each node's typical monthly range to the rms of
    its lines, kappa, which turns a lake's lines into its typical range."""
    from geography import tides as GT
    z = np.load(ROOT / TIDE_GRID)
    basis_ = z['basis'].astype(np.float64)
    lines = z['line_amplitude'].astype(np.complex128)
    amp, _, rank, _ = np.linalg.lstsq(basis_.T, lines.T, rcond=1e-6)
    amp = amp.T
    residual = float(np.abs(amp @ basis_ - lines).max() / np.abs(lines).max())
    sigma = np.sqrt((np.abs(lines) ** 2).sum(axis=0) / 2.0)
    ratio = z['monthly_range_median_m'] / np.maximum(sigma, 1e-12)
    lat, lon = z['lat_deg'][z['rows']], z['lon_deg'][z['cols']]
    return dict(amp=amp, rank=int(rank), residual=residual, kappa=float(np.median(ratio)),
                kappa_p10_p90=[float(np.percentile(ratio, 10)), float(np.percentile(ratio, 90))],
                periods_days=z['line_period_days'], tree=cKDTree(unit(lat, lon)), node_range=z['monthly_range_median_m'],
                node_lines=lines, node_label=z['labels'], monomials=GT.monomial_values)


def sea_amplitude(forcing, lat, lon):
    """Half the typical monthly range (m) of the nearest sea node, and that node."""
    _, node = forcing['tree'].query(unit(lat, lon))
    return 0.5 * forcing['node_range'][node].astype(np.float64), node


def lake_bodies(sh):
    """The rain-fed lakes' and the common-level lakes' bodies at 16 px/deg, labelled as shores.build labels them,
    with each body's area-weighted mean of the tide's monomials; the common-level bodies of the seas' size are
    flagged."""
    from geography import atlas as GA
    from geography import topography as tp
    atlas = json.loads((ROOT / S.ATLAS).read_text())
    level = float(atlas['sea_level_m'])
    lat, dy, dx, dx_edge, row_area = S.grid()
    h, _, _ = tp.height_above_geoid(S.PPD, atlas['grid']['geoid_degree'])
    common = h < level
    del h
    with np.load(ROOT / S.DRAINAGE) as z:
        cells = z['lake_cell'].astype(np.int64)
    lake = np.zeros(S.ROWS * S.COLS, bool)
    lake[cells] = True
    lake = lake.reshape(S.ROWS, S.COLS)
    eight = np.ones((3, 3), bool)
    rain_fed = S.merge_labels(ndimage.label(lake, structure=eight)[0])
    common_body = S.merge_labels(ndimage.label(common, structure=eight)[0])
    area = S.body_areas(common_body, row_area)
    sea = area >= GA.BODY_MIN_SHARE * 4.0 * math.pi * (MOON_RADIUS / 1000.0) ** 2
    sea[0] = False
    lon = (np.arange(S.COLS) + 0.5) / S.PPD
    out = {}
    for name, labels in (('rain_fed', rain_fed), ('common', common_body)):
        n = int(labels.max()) + 1
        sums = np.zeros((n, 20))
        weight = np.zeros(n)
        for r in range(S.ROWS):
            row = labels[r]
            sel = row > 0
            if not sel.any():
                continue
            v = monomials(unit(np.full(sel.sum(), lat[r]), lon[sel]))
            w = row_area[r]
            np.add.at(sums, row[sel], v * w)
            weight += np.bincount(row[sel], minlength=n) * w
        out[name] = dict(labels=labels, mean=sums / np.maximum(weight, 1e-30)[:, None], area_km2=weight)
    out['common']['sea'] = sea
    return out


def monomials(points):
    """geography.tides' monomials of the unit position."""
    from geography import tides as GT
    return GT.monomial_values(points)


def lake_amplitude(forcing, means, body, lat, lon, chunk=200000):
    """Half the typical monthly range (m) of a closed lake's equilibrium tide at shore points: the forcing's lines
    on the monomials less their mean over the lake (means[body]), which a lake of fixed volume holds, turned into a
    typical range by kappa."""
    lat, lon, body = np.atleast_1d(lat), np.atleast_1d(lon), np.atleast_1d(body).astype(np.int64)
    out = np.zeros(lat.size)
    for s in range(0, lat.size, chunk):
        dv = monomials(unit(lat[s:s + chunk], lon[s:s + chunk])) - means[body[s:s + chunk]]
        lines = dv @ forcing['amp'].T
        out[s:s + chunk] = 0.5 * forcing['kappa'] * np.sqrt((np.abs(lines) ** 2).sum(axis=-1) / 2.0)
    return out


def shore_tides(px, sh, forcing, bodies):
    """Each shore pixel's tide amplitude A (m): for the pixels beside the water, by the class of the water each edge
    faces (columns sea, common-level lake, open and closed lake); for the pixels inland, that of their nearest water."""
    co, il = px['coast'], px['inland']
    cidx = sh['coast_index'].astype(np.int64)
    lat_c, lon_c = pixel_latlon(cidx)
    a_coast = np.zeros((cidx.size, 4))
    a_coast[:, 0] = sea_amplitude(forcing, lat_c, lon_c)[0]
    rf = bodies['rain_fed']
    lake = co['lake'].astype(np.int64)
    a_lake = lake_amplitude(forcing, rf['mean'], lake, lat_c, lon_c)
    a_coast[:, 2] = np.where(lake > 0, a_lake, 0.0)
    a_coast[:, 3] = a_coast[:, 2]
    cm = bodies['common']
    flat = cm['labels'].ravel()
    r, c = np.divmod(cidx, S.COLS)
    facing = np.zeros(cidx.size, np.int64)
    for dr, dc in ((0, 1), (0, -1), (-1, 0), (1, 0)):
        rr = np.clip(r + dr, 0, S.ROWS - 1)
        lab = flat[rr * S.COLS + (c + dc) % S.COLS].astype(np.int64)
        lab = np.where(cm['sea'][lab], 0, lab)
        facing = np.maximum(facing, lab)
    a_coast[:, 1] = np.where(facing > 0, lake_amplitude(forcing, cm['mean'], facing, lat_c, lon_c), 0.0)
    # inland pixels: their nearest water pixel
    pidx = sh['pixel_index'].astype(np.int64)
    inland = np.ones(pidx.size, bool)
    inland[np.isin(pidx, cidx)] = False
    pr, pc = np.divmod(pidx[inland], S.COLS)
    wr = np.clip(pr + sh['pixel_step_rows'][inland], 0, S.ROWS - 1)
    wc = (pc + sh['pixel_step_cols'][inland]) % S.COLS
    widx = wr * S.COLS + wc
    wlat, wlon = pixel_latlon(widx)
    water = il['water']
    a_in = np.zeros(widx.size)
    sea_sel = water == S.SEA
    a_in[sea_sel] = sea_amplitude(forcing, wlat[sea_sel], wlon[sea_sel])[0]
    rf_sel = np.isin(water, (S.LAKE_OPEN, S.LAKE_CLOSED)) & (il['lake'] > 0)
    a_in[rf_sel] = lake_amplitude(forcing, rf['mean'], il['lake'][rf_sel], wlat[rf_sel], wlon[rf_sel])
    cl_sel = water == S.LAKE_COMMON
    body = flat[widx[cl_sel]].astype(np.int64)
    a_in[cl_sel] = np.where(body > 0, lake_amplitude(forcing, cm['mean'], body, wlat[cl_sel], wlon[cl_sel]), 0.0)
    return a_coast, a_in


def lake_summary(forcing, bodies, sh):
    """The lakes' largest typical monthly ranges: each lake's greatest equilibrium range on its own boundary pixels
    (its pixels with land beside them), by size, and the lakes that reach a decimetre."""
    out = {}
    for name in ('rain_fed', 'common'):
        b = bodies[name]
        labels = b['labels']
        water = labels > 0
        edge = water & ~ndimage.binary_erosion(water, structure=np.ones((3, 3), bool), border_value=1)
        idx = np.flatnonzero(edge)
        lab = labels.ravel()[idx].astype(np.int64)
        if name == 'common':
            keep = ~b['sea'][lab]
            idx, lab = idx[keep], lab[keep]
        lat, lon = pixel_latlon(idx)
        rng = 2.0 * lake_amplitude(forcing, b['mean'], lab, lat, lon)
        top = np.zeros(labels.max() + 1)
        np.maximum.at(top, lab, rng)
        area = b['area_km2']
        real = np.flatnonzero(area > 0)
        if name == 'common':
            real = real[~b['sea'][real]]
        rows = {}
        for lo, hi in ((0, 10), (10, 100), (100, 1000), (1000, 10000), (10000, 1e7)):
            sel = real[(area[real] >= lo) & (area[real] < hi)]
            if sel.size:
                rows[f'{lo:g}-{hi:g}_km2'] = dict(count=int(sel.size), median_greatest_range_m=float(np.median(top[sel])),
                                                  greatest_range_m=float(top[sel].max()))
        deci = real[top[real] >= 0.1]
        big = real[np.argsort(-area[real])[:8]]
        centres = ndimage.center_of_mass(np.ones(labels.shape, np.uint8), labels, big) if big.size else []
        out[name] = dict(lakes=int(real.size), by_area=rows, reaching_a_decimetre=dict(
            count=int(deci.size), area_km2=float(area[deci].sum()),
            smallest_km2=float(area[deci].min()) if deci.size else None),
            largest=[dict(area_km2=float(area[k]), greatest_range_m=float(top[k]),
                          lat=float(90.0 - (y + 0.5) / S.PPD), lon=float((x + 0.5) / S.PPD))
                     for k, (y, x) in zip(big, centres)])
    return out


# ------------------------------------------------------------------------------------------------- the numerics
def ferris_table():
    """Ferris's damping length and lag for the leading lines in each ground, at each end of its specific yield."""
    out = {}
    lines = dict(m_prime=TD.MONTH_S, fortnight_2d=14.765294 * JULIAN_DAY)
    for ground, t in (('slow', 2.51846e-4), ('central', 5.17172e-4), ('sandy', 5.17172e-3)):
        for sy in TD.SPECIFIC_YIELD[ground]:
            row = {}
            for name, period in lines.items():
                k = TD.damping(period, sy, t)
                amp, lag = TD.ferris(np.array(DISTANCES_M[1:], float), 1.0, period, sy, t)
                row[name] = dict(damping_length_m=1.0 / k, swing_share={f'{d}': float(a) for d, a in
                                                                        zip(DISTANCES_M[1:], amp)},
                                 lag_days={f'{d}': float(l / JULIAN_DAY) for d, l in zip(DISTANCES_M[1:], lag)})
            out[f'{ground}_sy_{sy:g}'] = dict(transmissivity_m2_s=t, specific_yield=sy, **row)
    return out


def numerics(transmissivities):
    """The no-draw overheight against eps, and for each of the coast-water study's cases the model's rise on planar
    shores over a grid of slope, tide amplitude and draw (central specific yield). Cached by tide.py's hash and the
    grid."""
    key = digest('research/studies/coast_water/tide.py') + json.dumps(
        [EPS_TABLE, RISE_CASES, RISE_SLOPES, RISE_AMPLITUDES, RISE_DRAWS_MM_DAY, R.rounded(transmissivities)])
    cache = ROOT / CACHE
    if cache.exists():
        old = json.loads(cache.read_text())
        if old.get('key') == key:
            return old
    t0 = time.time()
    no_draw = [TD.overheight_run(e) for e in EPS_TABLE]
    rise = {}
    for name in RISE_CASES:
        ground, reach = name.split('_')[0], float(name.split('_')[1][:-1])
        t = transmissivities[ground]
        sy = central_yield(ground)
        k = TD.damping(TD.MONTH_S, sy, t)
        rows = []
        for s in RISE_SLOPES:
            for a in RISE_AMPLITUDES:
                for e in RISE_DRAWS_MM_DAY:
                    draw = e / 1000.0 / 86400.0
                    eps = k * a / s
                    run = TD.belt_run(s, t, sy, draw, reach, a, steps=160, tolerance=3e-3, max_windows=300)
                    rows.append(dict(slope=s, amplitude_m=a, draw_mm_day=e, eps=eps,
                                     rise_m=run['with']['effective_rise_m'],
                                     rise_without_m=run['without']['effective_rise_m'],
                                     belt_without_m=run['without']['edge_mean_m'],
                                     analytic_m=run['analytic_without_m'], belt_with_m=run['with']['edge_mean_m'],
                                     belt_always_m=run['with']['edge_always_m'],
                                     draw_without_m3_yr_per_km=run['without']['draw_m3_yr_per_km'],
                                     draw_with_m3_yr_per_km=run['with']['draw_m3_yr_per_km'],
                                     swing_half_range_m=run['with']['swing_half_range_m'],
                                     draw_share_reached_at_high_water=run['with']['draw_share_reached_at_high_water'],
                                     windows=run['with']['windows'], last_change_m=run['with']['last_change_m']))
        rise[name] = dict(transmissivity_m2_s=t, specific_yield=sy, reach_m=reach, damping_length_m=1.0 / k,
                          rows=rows)
    out = dict(key=key, no_draw=no_draw, rise=rise, seconds=time.time() - t0)
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(R.rounded(out), indent=1))
    return json.loads(cache.read_text())


def rise_tables(num):
    """The rise (m) of each case at any slope, amplitude and draw (m/s): its log linear in the logs between the grid's
    runs, so that a power law between them is exact, held at the grid's edges; below its least amplitude as the amplitude squared (as the runs scale there), and below its
    least draw towards the no-draw overheight at zero draw."""
    from scipy.interpolate import RegularGridInterpolator
    ls, la, le = (np.log(np.array(v)) for v in (RISE_SLOPES, RISE_AMPLITUDES, RISE_DRAWS_MM_DAY))
    eps_t = np.array([r['eps'] for r in num['no_draw']])
    f_t = np.array([r['overheight_over_amplitude'] for r in num['no_draw']])
    tables = {}
    for name, case in num['rise'].items():
        v = np.log(np.maximum([r['rise_m'] for r in case['rows']], 1e-4)).reshape(ls.size, la.size, le.size)
        grid = RegularGridInterpolator((ls, la, le), v)
        k = 1.0 / case['damping_length_m']

        def f(slope, amp, draw, grid=grid, k=k):
            slope = np.clip(np.asarray(slope, float), RISE_SLOPES[0], RISE_SLOPES[-1])
            amp = np.asarray(amp, float)
            e = np.asarray(draw, float) * 1000.0 * 86400.0
            a_in = np.clip(amp, RISE_AMPLITUDES[0], RISE_AMPLITUDES[-1])
            e_in = np.clip(e, RISE_DRAWS_MM_DAY[0], RISE_DRAWS_MM_DAY[-1])
            q = np.stack(np.broadcast_arrays(np.log(slope), np.log(a_in), np.log(e_in)), axis=-1)
            out = np.exp(grid(q))
            low = amp < RISE_AMPLITUDES[0]
            out = np.where(low, out * (amp / RISE_AMPLITUDES[0]) ** 2, out)
            weak = e < RISE_DRAWS_MM_DAY[0]
            if np.any(weak):
                eps = k * amp / slope
                none = amp * np.interp(np.log(np.maximum(eps, 1e-6)), np.log(eps_t), f_t, left=0.0)
                none = np.where(eps < eps_t[0], 0.5 * eps * amp, none)
                share = np.clip(e / RISE_DRAWS_MM_DAY[0], 0.0, 1.0)
                out = np.where(weak, none + share * (out - none), out)
            return np.where(amp > 1e-3, np.maximum(out, 0.0), 0.0)
        tables[name] = f
    return tables


# ----------------------------------------------------------------------------------------- the belts with the tide
def belts_reach(px, uptake_cell, transmissivity, reach_coast, reach_inland):
    """run.belts_for with a reach for each pixel (each edge of the pixels beside the water)."""
    c_cell = uptake_cell.ravel() / (2.0 * transmissivity)
    co = px['coast']
    c = c_cell[co['cell']][:, None]
    areas, widths, lengths = B.first_ring(co['edge'], co['across'], co['rise'], c, reach_coast)
    total = areas.sum(axis=1)
    scale = np.where(total > co['area'], co['area'] / np.maximum(total, 1e-30), 1.0)
    areas = areas * scale[:, None]
    cls = np.stack([np.where(co['island'], 1, 0), np.full(total.size, 2), np.full(total.size, 3),
                    np.full(total.size, 4)], axis=1)
    il = px['inland']
    share = B.inland(il['distance'], il['step_rows'], il['step_cols'], il['height'], c_cell[il['cell']],
                     reach_inland)
    il_area = share * il['area']
    il_cls = R.classify(il['water'], il['island'])
    belt = np.zeros((5, 2048))
    length = np.zeros((5, 2048))
    for k in range(4):
        np.add.at(belt, (cls[:, k], co['cell']), areas[:, k])
        np.add.at(length, (cls[:, k], co['cell']), lengths[:, k])
    np.add.at(belt, (il_cls, il['cell']), il_area)
    strips = dict(width=widths, length=lengths, cls=cls, cell=co['cell'])
    return dict(belt=belt, length=length, strips=strips, inland_area=il_area, inland_cls=il_cls, coast_areas=areas)


def pixel_rises(px, uptake_cell, rise, a_coast, a_in):
    """Each pixel's rise of the mean water table under its tide (m), from its case's table at its own slope (the
    slope belts.py gives it), amplitude and draw."""
    co, il = px['coast'], px['inland']
    up = uptake_cell.ravel()
    with np.errstate(divide='ignore', invalid='ignore'):
        slope_c = np.where(co['across'] > 0, np.maximum(np.nan_to_num(co['rise'], nan=0.0), 0.0)
                           / (0.5 * co['across'] * 1000.0), 0.0)
        n = np.maximum(np.maximum(np.abs(il['step_rows']), np.abs(il['step_cols'])), 1).astype(float)
        centre = il['distance'] - 0.5 * il['distance'] / n
        slope_i = np.maximum(il['height'], 0.0) / (centre * 1000.0)
    rc = rise(slope_c, a_coast, (up[co['cell']])[:, None])
    ri = rise(slope_i, a_in, up[il['cell']])
    return rc, ri


def compare(summary_without, summary_with, regions=(DESERT, 'log_rain_mask.every_cell', 'share_matched.every_cell')):
    """The belts with the tide against the belts without it, by region and shore class, with the bands."""
    out = {}
    for region in regions:
        a = summary_without['deserts'][region]
        b = summary_with['deserts'][region]
        name = next(iter(a))
        a, b = a[name], b[name]
        rows = {}
        for cls in R.SHORE_CLASSES:
            x, y = a['by_class'][cls], b['by_class'][cls]
            rows[cls] = dict(shore_km=x['shore_km'], belt_km2=[x['belt_km2'], y['belt_km2']],
                             belt_change=(y['belt_km2'] / x['belt_km2'] - 1.0) if x['belt_km2'] else None,
                             width_median_m=[x['width_m']['median'], y['width_m']['median']],
                             width_p90_m=[x['width_m']['p90'], y['width_m']['p90']],
                             supply_m3_yr_per_km=[x['supply_m3_yr_per_km'], y['supply_m3_yr_per_km']],
                             uptake_km3_yr=[x['uptake_km3_yr'], y['uptake_km3_yr']])
        bands = {k: dict(share_within_reach=[a['bands'][k]['share_within_reach'], b['bands'][k]['share_within_reach']],
                         belt_km2=[a['bands'][k]['belt_km2'], b['bands'][k]['belt_km2']]) for k in a['bands']}
        out[region] = dict(by_class=rows, bands=bands, belt_km2=[a['belt_km2'], b['belt_km2']],
                           share_within_reach=[a['share_within_reach'], b['share_within_reach']],
                           uptake_km3_yr=[a['uptake_km3_yr'], b['uptake_km3_yr']])
    alls = summary_without['by_class']['all'][next(iter(summary_without['by_class']['all']))]
    allt = summary_with['by_class']['all'][next(iter(summary_with['by_class']['all']))]
    out['all_shores'] = {cls: dict(belt_km2=[alls['by_class'][cls]['belt_km2'], allt['by_class'][cls]['belt_km2']],
                                   uptake_km3_yr=[alls['by_class'][cls]['uptake_km3_yr'],
                                                  allt['by_class'][cls]['uptake_km3_yr']])
                         for cls in R.SHORE_CLASSES}
    return out


def tide_distribution(px, a_coast, a_in, masks, grid):
    """The tides the shores carry: by class, the typical monthly range (2A) over shore length, Moon-wide and in the
    deserts, and the rises' spread."""
    co = px['coast']
    desert = masks[('rain_mask', 'every_cell')].ravel()
    cls = np.stack([np.where(co['island'], 1, 0), np.full(co['cell'].size, 2), np.full(co['cell'].size, 3),
                    np.full(co['cell'].size, 4)], axis=1)
    out = {}
    for region, cells in (('all', np.ones(2048, bool)), ('deserts', desert)):
        rows = {}
        for k, name in enumerate(R.SHORE_CLASSES):
            sel = (cls == k) & (co['edge'] > 0) & cells[co['cell']][:, None]
            w = co['edge'][sel]
            if w.sum() <= 0:
                continue
            rng = 2.0 * a_coast[sel]
            rows[name] = dict(shore_km=float(w.sum() * B.STAIRCASE),
                              range_m=dict(zip(('p10', 'median', 'p90'), R.weighted_percentiles(rng, w))),
                              share_over_decimetre=float(w[rng >= 0.1].sum() / w.sum()),
                              share_over_metre=float(w[rng >= 1.0].sum() / w.sum()))
        out[region] = rows
    return out


def moon_wide(b, forcing, bodies, tables):
    px = R.shore_pixels(b['sh'])
    a_coast, a_in = shore_tides(px, b['sh'], forcing, bodies)
    out = dict(tides=tide_distribution(px, a_coast, a_in, b['masks'], b['grid']), cases={})
    desert = b['masks'][('rain_mask', 'every_cell')].ravel()
    for name, case in b['cases'].items():
        t, reach = case['transmissivity_m2_s'], case['reach_m']
        base = R.belts_for(px, b['uptake'], t, reach)
        s0 = R.belt_summary({name: base}, b['uptake'], b['masks'], b['grid'], b['shore_cells'], {name: case})
        del base
        rc, ri = pixel_rises(px, b['uptake'], tables[name], a_coast, a_in)
        run = belts_reach(px, b['uptake'], t, reach + rc, reach + ri)
        s1 = R.belt_summary({name: run}, b['uptake'], b['masks'], b['grid'], b['shore_cells'], {name: case})
        del run
        entry = compare(s0, s1)
        co = px['coast']
        for k, label in ((0, 'sea'), (1, 'common_level_lake'), (2, 'open_lake'), (3, 'closed_lake')):
            sel = (co['edge'][:, k] > 0) & desert[co['cell']]
            entry[f'desert_{label}_rise_m'] = dict(zip(('p10', 'median', 'p90'), R.weighted_percentiles(
                rc[sel, k], co['edge'][sel, k])))
        out['cases'][name] = dict(specific_yield=central_yield(name.split('_')[0]), **entry)
    return out


# ------------------------------------------------------------------------------------- desert coasts at 100 m
SITE_HALF_M = 6000.0
SPACING_M = 100.0
TRANSECTS = 4
PROFILE_STEP_M = 10.0
PROFILE_LAND_M = 2500.0
PROFILE_SEA_M = 400.0


def sites(entries, b):
    """Up to two desert sea-coast sites in each slice: the 16 px/deg desert sea-coast pixels in the slice's rows,
    binned by degree of longitude, the two densest bins at least 10 degrees apart, each at its pixels' middle."""
    sh = b['sh']
    co = sh['coast_index'].astype(np.int64)
    lat, lon = pixel_latlon(co)
    cell = (sh['t21_row'][co // S.COLS].astype(np.int64) * 64 + sh['t21_col'][co % S.COLS].astype(np.int64))
    sea = sh['coast_edge_km'][:, 0] > 0
    out = []
    for e in entries:
        south, north = e['latitude_deg']
        half = 'east' if e['tile'].endswith('000_180') else 'west'
        for mask in ('rain_mask', 'log_rain_mask', 'share_matched'):
            m = b['masks'][(mask, 'every_cell')].ravel()[cell] & sea & (lat > south + 0.3) & (lat < north - 0.3)
            m &= (lon < 180.0) if half == 'east' else (lon >= 180.0)
            if m.sum() >= 20:
                break
        if m.sum() < 20:
            continue
        bins = np.floor(lon[m]).astype(int)
        counts = np.bincount(bins, minlength=361)
        chosen = []
        for j in np.argsort(-counts):
            if counts[j] < 10 or len(chosen) == 2:
                break
            if all(abs(j - k) >= 10 for k in chosen):
                chosen.append(int(j))
        for j in chosen:
            s = m & (np.floor(lon) == j)
            k = np.flatnonzero(s)[np.argmin(np.abs(lat[s] - np.median(lat[s])) + np.abs(lon[s] - np.median(lon[s])))]
            out.append(dict(slice=e['name'], mask=mask, lat=float(lat[k]), lon=float(lon[k]),
                            cell=int(cell[k]), desert_sea_pixels_in_degree=int(counts[j])))
    return out


def fine_belts_reach(height, spacing_m, c, reach):
    """run.fine_belts with a reach for each pixel of the grid."""
    water = height < 0.0
    dist, (ii, jj) = ndimage.distance_transform_edt(~water, sampling=spacing_m, return_indices=True)
    rows, cols = np.indices(height.shape)
    steps_r, steps_c = ii - rows, jj - cols
    land = ~water
    coastal = np.zeros_like(land)
    edge = np.zeros(height.shape)
    for dr, dc in ((0, 1), (0, -1), (1, 0), (-1, 0)):
        nb = np.zeros_like(water)
        src = water[max(dr, 0):height.shape[0] + min(dr, 0), max(dc, 0):height.shape[1] + min(dc, 0)]
        nb[max(-dr, 0):height.shape[0] + min(-dr, 0), max(-dc, 0):height.shape[1] + min(-dc, 0)] = src
        facing = land & nb
        coastal |= facing
        edge += facing * spacing_m / 1000.0
    pixel_area = (spacing_m / 1000.0) ** 2
    across = np.full(height.shape, spacing_m / 1000.0)
    areas, widths, lengths = B.first_ring(edge, across, height, c, reach)
    first = np.minimum(areas, pixel_area)
    inner = land & ~coastal
    inland_area = np.zeros(height.shape)
    inland_area[inner] = B.inland(dist[inner] / 1000.0, steps_r[inner], steps_c[inner], height[inner], c,
                                  reach[inner]) * pixel_area
    return np.where(coastal, first, 0.0) + inland_area, dist, coastal


def transects(height, axis, dist, coastal, count):
    """Cross-shore profiles from count shoreline pixels spread along the coast: inland along the distance field's
    gradient, seaward the opposite way, heights bilinear every PROFILE_STEP_M."""
    gy, gx = np.gradient(dist)
    idx = np.flatnonzero(coastal & (np.abs(axis)[None, :] < SITE_HALF_M - PROFILE_LAND_M / 2)
                         & (np.abs(axis)[:, None] < SITE_HALF_M - PROFILE_LAND_M / 2))
    if idx.size == 0:
        return []
    pick = idx[np.linspace(0, idx.size - 1, min(count, idx.size)).astype(int)]
    s = np.arange(-PROFILE_SEA_M, PROFILE_LAND_M + PROFILE_STEP_M / 2, PROFILE_STEP_M)
    out = []
    for p in pick:
        r, c = divmod(int(p), height.shape[1])
        d = np.array([gx[r, c], gy[r, c]])
        if not np.isfinite(d).all() or np.hypot(*d) == 0:
            continue
        d = d / np.hypot(*d)
        cols = c + s / SPACING_M * d[0]
        rows = r + s / SPACING_M * d[1]
        z = ndimage.map_coordinates(height.astype(float), [rows, cols], order=1, mode='nearest')
        # the shoreline: the last crossing below the mean water level before the land
        wet = np.flatnonzero(z[: int(PROFILE_SEA_M / PROFILE_STEP_M) + 20] < 0)
        if wet.size == 0:
            continue
        x = s - s[wet[-1]] - 0.5 * PROFILE_STEP_M
        out.append((x, z))
    return out


PROFILE_CAP = 120                 # windows a profile run may take; the drift over its last ten is recorded


def profile_run(x, z, transmissivity, sy, draw, reach, amplitude, rise):
    """The cross-shore model on a real profile without and with the tide, each to three of its belt's settling
    times and the geometric tail, or to the cap. The run with the tide starts from the run without it, raised by
    the planar table's rise at the profile's mean slope over its first kilometre. Each run is repeated to ten
    windows short of its end, and its drift is the change of the belt's edge and of the mean water table over the
    belt across those ten windows."""
    c = draw / (2.0 * transmissivity)
    span = min(1000.0, float(x[-1]))
    slope = max(float(np.interp(span, x, z)) / span, 1e-3)
    w0 = float(B.width(slope, c, reach))
    w_guess = float(B.width(slope, c, reach + 4.0 * amplitude))
    settle = max(4, int(math.ceil(1.2 * w_guess ** 2 * sy / transmissivity / TD.MONTH_S)))
    lift = float(np.ravel(rise(slope, amplitude, draw))[0])
    start = np.minimum(-c * (2.0 * w0 * np.clip(x, 0, w0) - np.clip(x, 0, w0) ** 2), np.maximum(z, 0.0))
    out = dict(slope=slope, planar_rise_m=lift)
    for label, a in (('without', 0.0), ('with', amplitude)):
        def go(cap, a=a, start=start):
            return TD.cross_shore(x, z, sy, transmissivity, TD.sinusoid(a), draw=draw, reach=reach,
                                  dt=TD.MONTH_S / 120, window=TD.MONTH_S, start=start, tolerance=5e-3,
                                  min_windows=min(settle, cap), max_windows=cap)
        r = go(PROFILE_CAP)
        r10 = go(max(r['windows'] - 10, 1))
        edge = TD.edge(x, z - r['mean'], reach)
        belt = (x >= 0) & (x <= edge)
        roots = (z - r['high']) <= TD.ROOTS_M
        out[label] = dict(edge_m=edge, edge_always_m=TD.edge(x, z - r['low'], reach), windows=r['windows'],
                          edge_drift_m=edge - TD.edge(x, z - r10['mean'], reach),
                          water_table_drift_m=float(np.mean((r['mean'] - r10['mean'])[belt])) if belt.any() else 0.0,
                          swing_half_range_100m_m=float(np.interp(100.0, x, 0.5 * (r['high'] - r['low']))),
                          draw_share_reached_at_high_water=float(r['draw'][roots].sum() / max(r['draw'].sum(), 1e-30)))
        start = np.minimum(r['mean'] + lift, np.maximum(z, 0.0))      # the tide's run: this belt raised
    return out


def profiles(b, forcing, tables):
    entries = TT.restore()
    out = dict(slices=[dict(name=e['name'], latitude_deg=e['latitude_deg'], tile=e['tile']) for e in entries],
               sites=[])
    by_name = {e['name']: e for e in entries}
    px = R.shore_pixels(b['sh'])
    for site in sites(entries, b):
        e = by_name[site['slice']]
        axis, height, _ = TT.local_grid(e, site['lon'], site['lat'], SITE_HALF_M, SPACING_M)
        height = height.astype(np.float64)
        e_cell = float(b['uptake'].ravel()[site['cell']])
        a = float(sea_amplitude(forcing, np.array([site['lat']]), np.array([site['lon']]))[0][0])
        rec = dict(site, uptake_mm_day=e_cell * 1000 * 86400, tide_range_m=2 * a,
                   land_share=float((height >= 0).mean()), cases={})
        for name in ('central_15m', 'sandy_25m'):
            case = b['cases'][name]
            ground = name.split('_')[0]
            t, reach = case['transmissivity_m2_s'], case['reach_m']
            sy = central_yield(ground)
            c = e_cell / (2.0 * t)
            zero = np.full(height.shape, reach)
            area0, dist, coastal = fine_belts_reach(height, SPACING_M, c, zero)
            # each pixel's rise at its slope from the shore, the square's water below the sea level taken as the sea
            slope = np.maximum(height, 0.0) / np.maximum(dist - 0.5 * SPACING_M, 0.5 * SPACING_M)
            rise = tables[name](slope, np.full(height.shape, a), np.full(height.shape, e_cell))
            area1, _, _ = fine_belts_reach(height, SPACING_M, c, reach + rise)
            # the 16 px/deg strips in the same square, with and without the tide
            inner = SITE_HALF_M
            co_lat, co_lon = pixel_latlon(b['sh']['coast_index'].astype(np.int64))
            x_, y_ = R.aeqd(site['lat'], site['lon'], co_lat, co_lon)
            in_c = (np.abs(x_) <= inner) & (np.abs(y_) <= inner)
            il_idx = px['inland']
            il_lat, il_lon = 90.0 - (il_idx['row'] + 0.5) / S.PPD, (il_idx['col'] + 0.5) / S.PPD
            xi, yi = R.aeqd(site['lat'], site['lon'], il_lat, il_lon)
            in_i = (np.abs(xi) <= inner) & (np.abs(yi) <= inner) & np.isin(il_idx['water'], (S.SEA, S.LAKE_COMMON))
            tile = R.subset(px, in_c, in_i)
            const = np.full((32, 64), e_cell)
            coarse0 = R.belts_for(tile, const, t, reach)
            a_c = np.zeros((int(in_c.sum()), 4))
            a_c[:, :2] = a                                   # the square's water below the sea level is the sea's
            rc, ri = pixel_rises(tile, const, tables[name], a_c, np.full(int(in_i.sum()), a))
            coarse1 = belts_reach(tile, const, t, reach + rc, reach + ri)
            runs = [profile_run(x, z, t, sy, e_cell, reach, a, tables[name])
                    for x, z in transects(height, axis, dist, coastal, TRANSECTS)]
            med = lambda label, key: float(np.median([q[label][key] for q in runs])) if runs else None
            unsettled = [dict(transect=i, **{f'{label}_{k}': q[label][k] for label in ('without', 'with')
                                             for k in ('windows', 'edge_drift_m', 'water_table_drift_m')})
                         for i, q in enumerate(runs) if any(q[label]['windows'] >= PROFILE_CAP and (
                             abs(q[label]['edge_drift_m']) > 1.0 or abs(q[label]['water_table_drift_m']) > 0.01)
                             for label in ('without', 'with'))]
            rec['cases'][name] = dict(
                belt_km2_100m=[float(area0.sum()), float(area1.sum())],
                belt_km2_16ppd=[float(coarse0['coast_areas'][:, :2].sum() + coarse0['inland_area'].sum()),
                                float(coarse1['coast_areas'][:, :2].sum() + coarse1['inland_area'].sum())],
                rise_m_100m=dict(zip(('p10', 'median', 'p90'), [float(np.percentile(rise[coastal], q))
                                                               for q in (10, 50, 90)])) if coastal.any() else None,
                transects=dict(count=len(runs), edge_without_m=med('without', 'edge_m'), edge_with_m=med('with', 'edge_m'),
                               edge_always_m=med('with', 'edge_always_m'),
                               swing_at_100m_m=med('with', 'swing_half_range_100m_m'),
                               draw_share_reached_at_high_water=med('with', 'draw_share_reached_at_high_water'),
                               largest_edge_drift_m=max((abs(q[l]['edge_drift_m']) for q in runs for l in ('without', 'with')), default=None),
                               largest_water_table_drift_m=max((abs(q[l]['water_table_drift_m']) for q in runs for l in ('without', 'with')),
                                                               default=None),
                               unsettled=unsettled, runs=runs))
        out['sites'].append(rec)
    return out


# --------------------------------------------------------------------------------------------- checks and load
def checks(transmissivities, uptake_median, day_share, forcing, tables):
    """Direct runs against the rise table, the two leading lines against one, the lunar day's draw, a finite
    aquifer, and Ferris's lag."""
    out = dict(direct=[], specific_yield=[])
    t_c = transmissivities['central']
    for ground, reach, slope, a, e in (('central', 15, 0.03, 1.85, uptake_median),
                                       ('central', 25, 0.1, 1.05, 1.3 * uptake_median),
                                       ('slow', 25, 0.01, 0.5, 0.7 * uptake_median),
                                       ('sandy', 25, 0.012, 1.85, uptake_median),
                                       ('sandy', 15, 0.04, 2.2, 1.2 * uptake_median),
                                       ('central', 15, 0.15, 2.4, uptake_median)):
        t = transmissivities[ground]
        sy = central_yield(ground)
        c = e / (2.0 * t)
        r = TD.belt_run(slope, t, sy, e, float(reach), a, steps=160, tolerance=3e-3)
        table = float(np.ravel(tables[f'{ground}_{reach}m'](slope, a, e))[0])
        out['direct'].append(dict(case=f'{ground}_{reach}m', slope=slope, amplitude_m=a, draw_mm_day=e * 1000 * 86400,
                                  rise_m=r['with']['effective_rise_m'], rise_table_m=table,
                                  belt_without_m=r['without']['edge_mean_m'], belt_with_m=r['with']['edge_mean_m'],
                                  belt_table_m=float(B.width(slope, c, reach + table)),
                                  belt_always_m=r['with']['edge_always_m']))
    # the specific yield's range on the desert sea coasts' median and gentle shores, central crust, 15 m
    for slope in (0.096, 0.018):
        for sy in TD.SPECIFIC_YIELD['central']:
            r = TD.belt_run(slope, t_c, sy, uptake_median, 15.0, 1.53, steps=160, tolerance=3e-3)
            out['specific_yield'].append(dict(slope=slope, specific_yield=sy, amplitude_m=1.53,
                                              belt_without_m=r['without']['edge_mean_m'],
                                              belt_with_m=r['with']['edge_mean_m'],
                                              swing_half_range_m=r['with']['swing_half_range_m']))
    # the two leading lines (M' and F) of a Frigoris coast node against a sinusoid of its typical range
    node = int(forcing['tree'].query(unit(np.array([56.0]), np.array([0.0])))[1][0])
    lines = forcing['node_lines'][:, node]
    per = forcing['periods_days'] * JULIAN_DAY
    a_typ = 0.5 * float(forcing['node_range'][node])
    sy = central_yield('central')
    beat = 1.0 / abs(1.0 / per[1] - 1.0 / per[0])
    two = TD.lines(np.abs(lines[:2]), per[:2], np.angle(lines[:2]))
    r2 = TD.belt_run(0.02, t_c, sy, uptake_median, 15.0, a_typ, sea=two, window=beat, steps=int(beat / 6.6 / 3600),
                     tolerance=3e-3, max_windows=12, start_rise=0.8 * a_typ)
    r1 = TD.belt_run(0.02, t_c, sy, uptake_median, 15.0, a_typ, steps=160, tolerance=3e-3)
    out['two_lines'] = dict(node_lat_lon=[56.0, 0.0], typical_range_m=2 * a_typ,
                            line_amplitudes_m=[float(v) for v in np.abs(lines[:2])], beat_days=beat / JULIAN_DAY,
                            rise_two_lines_m=r2['with']['effective_rise_m'], rise_one_line_m=r1['with']['effective_rise_m'],
                            belt_two_lines_m=r2['with']['edge_mean_m'], belt_one_line_m=r1['with']['edge_mean_m'],
                            always_two_lines_m=r2['with']['edge_always_m'], always_one_line_m=r1['with']['edge_always_m'],
                            windows=r2['with']['windows'])
    # the plants' draw through the lunar day, with and without the tide
    beat = 1.0 / abs(1.0 / TD.MONTH_S - 1.0 / TD.LUNAR_DAY_S)
    cycle = TD.lunar_day_draw(uptake_median, day_share)
    rows = {}
    for label, amp in (('no_tide', 0.0), ('tide', 1.85)):
        r = TD.belt_run(0.02, t_c, sy, uptake_median, 15.0, max(amp, 1e-6), draw_t=cycle, window=beat,
                        steps=int(beat / 3.0 / 3600), tolerance=3e-3, max_windows=30)
        rows[label] = dict(belt_mean_m=r['with']['edge_mean_m'], belt_always_m=r['with']['edge_always_m'],
                           swing_half_range_m=r['with']['swing_half_range_m'])
    out['lunar_day_draw'] = dict(day_share=day_share, beat_days=beat / JULIAN_DAY, slope=0.02, **rows,
                                 steady_belt_m=float(B.width(0.02, uptake_median / (2 * t_c), 15.0)))
    # a finite aquifer 100 m deep against the constant transmissivity
    rd = TD.belt_run(0.02, t_c, sy, uptake_median, 15.0, 1.85, steps=160, tolerance=3e-3, depth=100.0)
    out['finite_aquifer_100m'] = dict(belt_without_m=rd['without']['edge_mean_m'], belt_with_m=rd['with']['edge_mean_m'],
                                      rise_m=rd['with']['effective_rise_m'],
                                      finite_amplitude_overheight_m=float(TD.finite_amplitude_overheight(1.85, 100.0)))
    # Ferris's amplitude and lag on a near-vertical face, from the model's own record at 1/k
    k = TD.damping(TD.MONTH_S, sy, t_c)
    x, ground = TD.planar_grid(1.0, 0.1, k, 6.0 / k)
    at = int(np.searchsorted(x, 1.0 / k))
    level = TD.sinusoid(0.1)
    res = TD.cross_shore(x, ground, sy, t_c, level, dt=TD.MONTH_S / 480, window=TD.MONTH_S, tolerance=1e-6,
                         max_windows=8, record=np.array([at]))
    out['ferris_check'] = dict(x_m=float(x[at]), swing_share=float(0.5 * (res['high'][at] - res['low'][at]) / 0.1),
                               expected_share=math.exp(-k * x[at]),
                               lag_days=TD.lag(res['series'][:, 0], res['times'], TD.MONTH_S, level) / JULIAN_DAY,
                               expected_lag_days=k * x[at] / (2.0 * math.pi / TD.MONTH_S) / JULIAN_DAY)
    return out


def load_balance(b, cw, chemistry, num):
    seas_ = L.seas(chemistry)
    paths = {k: [L.evaporation_path(v, p) for p in L.SOIL_CO2_ATM] for k, v in seas_.items()}
    for k in paths:
        for p in paths[k]:
            p.pop('path', None)
    desert = cw['belts']['deserts'][DESERT]
    by_class = {}
    for case in ('central_15m', 'central_25m', 'sandy_25m'):
        rows = {}
        for cls, v in desert[case]['by_class'].items():
            if not v['belt_km2']:
                continue
            draw = v['uptake_km3_yr'] * 1e9 / (v['belt_km2'] * 1e6)          # m of water a year over the belt
            rows[cls] = dict(draw_m_yr=draw, **{k: L.deposits(s, draw) for k, s in seas_.items()},
                             leaching_m_yr_with_chloride=L.leaching_requirement(seas_['held_with_chloride'], draw))
        by_class[case] = rows
    rain = cw['demand']['desert_rain_mm_day'] * 365.25 / 1000.0
    flushed = {f"{r['slope']:g}": r['draw_share_reached_at_high_water'] for r in num['rise']['central_15m']['rows']
               if r['amplitude_m'] == 1.6 and r['draw_mm_day'] == 2.0}
    return dict(seas={k: {a: (v * 1000.0 if a in ('ca', 'mg', 'na', 'cl', 'so4', 'si', 'ta') else v)
                          for a, v in s.items()} for k, s in seas_.items()},
                evaporation_paths=paths, by_class=by_class, desert_rain_m_yr=rain,
                tide_reach_of_the_draw=dict(case='central_15m', amplitude_m=1.6, draw_mm_day=2.0, by_slope=flushed),
                constants=dict(soil_co2_atm=list(L.SOIL_CO2_ATM), cement_porosity=L.CEMENT_POROSITY,
                               root_zone_m=L.ROOT_ZONE_M, cec_mol_kg=list(L.CEC_MOL_KG), sodic_esp=L.SODIC_ESP,
                               saline_mol_m3=L.SALINE_MOL_M3, chloride_g_kg=L.CHLORIDE_G_KG))


def results():
    t0 = time.time()
    b = basis()
    data = {p: json.loads((ROOT / p).read_text()) for p in INPUTS}
    for p, schema in INPUTS.items():
        if data[p].get('schema') != schema:
            raise ValueError(f'{p}: schema {data[p].get("schema")} is not {schema}')
    cw = data[COAST_WATER]
    transmissivities = R.transmissivities(b['data'][R.TABLE])
    wl, desert = b['grid']['wl'], b['masks'][('rain_mask', 'every_cell')]
    uptake_median = R.weighted_percentiles(b['uptake'], wl * desert)[1]
    dm = wl * desert
    day_share = float((dm * b['pet']['day']).sum() / (dm * b['pet']['all']).sum())
    num = numerics(transmissivities)
    tables = rise_tables(num)
    print(f'{time.time() - t0:.0f} s: grids ready', file=sys.stderr, flush=True)
    forcing = tide_forcing()
    bodies = lake_bodies(b['sh'])
    lakes = lake_summary(forcing, bodies, b['sh'])
    mw = moon_wide(b, forcing, bodies, tables)
    print(f'{time.time() - t0:.0f} s: Moon-wide belts', file=sys.stderr, flush=True)
    del bodies
    prof = profiles(b, forcing, tables)
    print(f'{time.time() - t0:.0f} s: profiles', file=sys.stderr, flush=True)
    chk = checks(transmissivities, uptake_median, day_share, forcing, tables)
    load = load_balance(b, cw, data[CHEMISTRY], num)
    terrain = json.loads((ROOT / TERRAIN_RECORD).read_text())
    out = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
               basis=dict(target_k=b['base']['target_k'], share_of_the_way=b['base']['share_of_the_way'],
                          desert_uptake_median_mm_day=uptake_median * 1000 * 86400,
                          desert_day_share_of_demand=day_share, transmissivity_m2_s=transmissivities,
                          specific_yield=TD.SPECIFIC_YIELD, month_days=TD.MONTH_S / JULIAN_DAY),
               ferris=ferris_table(),
               numerics=dict(no_draw=num['no_draw'], rise=num['rise'], seconds=num['seconds']),
               tide_forcing=dict(rank=forcing['rank'], residual=forcing['residual'], kappa=forcing['kappa'],
                                 kappa_p10_p90=forcing['kappa_p10_p90']),
               lakes=lakes, belts=mw, profiles=prof, checks=chk, load=load,
               seconds=time.time() - t0)
    ignored = (R.GCM, R.GCM_SIDECAR, R.SHORES, TIDE_GRID) + LAKE_INPUTS
    out['producer'] = dict(study='research/studies/coast_water', files={p: digest(p) for p in CODE},
                           inputs={p: digest(p) for p in tuple(INPUTS) + tuple(R.INPUTS)},
                           ignored_inputs={p: digest(p) for p in ignored},
                           terrain_slices={e['name']: e.get('selected_sha256') for e in terrain['slices']},
                           gcm_files=b['sidecar']['files'])
    return out


def main(argv=None) -> int:
    product = R.rounded(results())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    c = product['belts']['cases']['central_15m'][DESERT]
    print('desert belts without and with the tide, central 15 m:', c['belt_km2'], c['share_within_reach'],
          'seconds', product['seconds'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
