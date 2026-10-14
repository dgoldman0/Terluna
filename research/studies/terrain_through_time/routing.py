"""The terrain study's geography: lakes, hollows, relief and channels from the drainage routing, and how they fill.

    python -m research.studies.terrain_through_time.routing   # results/terrain_routing.json, 3-6 minutes, 6 GB

The author asked on 13 October 2026 how much of the young Moon's terrain, its lakes, ponds and landforms, lasts
without upkeep and what could rebuild it. This runner repeats geography's drainage routing (geography.drainage.compute
with the climatology the committed drainage product names), checks that it reproduces geography/results/drainage.json
as resources' waters product does, and on it builds:

- where the land's rock comes off: in proportion to its local relief, the spread of heights within about 5 km of each
  cell, as Ahnert (1970) found denudation rising with relief, and to rain raised to the BQART model's exponent on
  water discharge (Syvitski & Milliman 2007), with the conservation study's total;
- the units that catch sediment: the spilling lakes as resources' retention bodies, and each closed depression with
  its lake or dry floor; for each its capacity, lake area, water inflow and the unit that takes what it passes, through
  the spilling lake's outlet or the closed depression's pass once it is full;
- each unit's fill through time from the rock eroded on the land that drains to it first, with Siyam's form of Brune's
  curves at lunar settling keeping sediment in the spilling lakes as the sea phosphorus study keeps particles, closed
  depressions keeping all of it, accretion on the lake floors, and full units passing everything on;
- the lift and distance from each eroding cell to where its sediment first settles, for carrying it back;
- the channels by slope and by discharge times slope, for where dam builders' dams would hold;
- the rivers' hardness relative to the weathering law's mean, from the land's temperature through Dessert's law, with
  rain and evaporation on the lakes diluting or concentrating it, for travertine.

The land loses the conservation study's 0.71, 3.55 and 7.1 Gt of rock a year, its 10, 50 and 100 m a million years
over the land. The per-place results
(each unit's fill times, the channels, local relief) go to research/runs/resources/terrain, ignored by Git, for the
atlas round; research/studies/terrain_through_time/run.py reads the summary.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'terrain_routing.json'
GRIDS = ROOT / 'research' / 'runs' / 'resources' / 'terrain' / 'terrain_16ppd.npz'
SCHEMA = 'terluna.research.terrain-routing/1'
DRAINAGE = 'geography/results/drainage.json'
ATLAS = 'geography/results/atlas.json'
GATES = 'research/studies/conservation/results/gates.json'
INPUTS = (DRAINAGE, ATLAS, GATES)
CODE = ('research/studies/terrain_through_time/routing.py', 'resources/sea_phosphorus_waters.py',
        'research/studies/carbon_cycle/run.py', 'geography/drainage.py', 'geography/topography.py',
        'geography/atlas.py')
IGNORED = ('climate/gcm/products/climatology_A28_dim5_moon.npz', 'geography/inputs/ldem_16.img',
           'geography/inputs/ldem_16.lbl', 'geography/inputs/ldem_4.img', 'geography/inputs/ldem_4.lbl',
           'geography/inputs/jggrx_0420a_sha.tab', 'geography/inputs/jggrx_0420a_sha.lbl')
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
GRAVITY_RATIO = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
ROCK_DENSITY = 2600.0                    # kg/m3, the conservation study's eroded rock (screens.sediment_return_watts)
RELIEF_HALF_WIDTH_KM = 5.0               # local relief over about 10 km (Montgomery & Brandon 2002's window)
WETLAND_DEPTH_M = 2.0                    # design guess: open water gives way to marsh where under about 2 m deep
HEADS_SHARE = 1.0 / 3.0                  # the heads variant: sediment back on each catchment's highest third by area
WATER_EXPONENT = 0.31                    # BQART's exponent on water discharge (Syvitski & Milliman 2007)
# Accretion on lake floors from what the lakes grow and settle themselves (organic matter, shells, diatoms and
# marl), mm a year: about a deep oligotrophic lake's (Baikal, 0.14), northern peat's long-term growth (0.41) and
# Earth's natural lakes' median before 1900 with their catchments' sediment (0.92); sources in sources.json.
ACCRETION_MM_YR = dict(low=0.1, central=0.4, high=1.0)
# Each case pairs the conservation study's erosion with the sea chemistry's dry density of what settles and the
# accretion: fast erosion with loose sediment and fast accretion fills earliest, slow with dense and slow latest.
CASES = dict(
    earliest=dict(erosion_m_myr=100.0, dry_density_kg_m3=500.0, accretion_mm_yr=ACCRETION_MM_YR['high']),
    central=dict(erosion_m_myr=50.0, dry_density_kg_m3=1000.0, accretion_mm_yr=ACCRETION_MM_YR['central']),
    latest=dict(erosion_m_myr=10.0, dry_density_kg_m3=1500.0, accretion_mm_yr=ACCRETION_MM_YR['low']),
    central_rock_only=dict(erosion_m_myr=50.0, dry_density_kg_m3=1000.0, accretion_mm_yr=0.0),
    central_by_runoff=dict(erosion_m_myr=50.0, dry_density_kg_m3=1000.0, accretion_mm_yr=ACCRETION_MM_YR['central'],
                           pattern='relief_and_runoff'))
STEPS_PER_DECADE = 20
END_YEARS = 1e9
SURVIVAL_YEARS = (100, 300, 1000, 3000, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9)
DECADES = (10, 100, 1000, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9)
FIFTHS = tuple(round(10 ** (k / 5), 6) for k in range(5, 46))   # five steps a decade, 10 to 1e9 years
AREA_CLASSES_KM2 = (0.0, 10.0, 100.0, 1000.0, 10000.0, math.inf)
LAT_BANDS_DEG = (0.0, 15.0, 30.0, 45.0, 60.0, 90.0)
SLOPE_BINS_PERCENT = (0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 10.0, 15.0, 23.0, math.inf)
LOG_QS_STEP = 0.05                       # decades of discharge times slope (m3/s) in the channel histogram
LOG_HARDNESS_STEP = 0.02                 # decades of hardness relative to the law's mean
SLOPE_FLOOR = 1e-5
EVIDENCE = ('The geography drainage routing rerun with its committed climatology and checked against the committed '
            'drainage product, on LOLA at 16 pixels a degree (1.9 km). Erosion is the conservation study\'s rock '
            'flux spread by local relief (Ahnert 1970) and rain to BQART\'s exponent (Syvitski & Milliman 2007), with '
            'relief alone and runoff kept as sensitivities; sediment settles on Siyam\'s (2000) form of '
            'Brune\'s (1953) median curve at lunar settling in spilling lakes and wholly in closed depressions, with '
            'design-guess accretion on lake floors; units fill in log-spaced steps of a twentieth of a decade with '
            'constant erosion and no uplift. Hardness is Dessert et al.\'s (2003) temperature term over the climate '
            'run\'s surface temperature, relative to the runoff-weighted mean. Features below the 1.9 km grid are '
            'outside it.')
READING_RULE = ('Units: 0 is the common-level waters; spilling lakes are resources\' retention bodies; each closed '
                'depression is one unit, a closed lake where it holds water and a dry hollow where not. A lake is '
                'lost when its remaining capacity falls below its area times 2 m or half its water, whichever is '
                'less, and full when none remains. Survival counts lakes still open after each year; lakes are '
                'classed by their area at the start. Rock rates are kg a year for each Gt a year the land loses, and '
                'lifts in m and distances in km are averages weighted by eroded mass. Channel lengths '
                'are km of 16 pixel/degree channel cells (mean discharge at least 1 m3/s, not in lakes), binned by '
                'slope in percent and by log10 of discharge times slope in m3/s; hardness is the channel water\'s '
                'dissolved load relative to the weathering law at the runoff-weighted mean temperature, before the '
                'young soils\' factor. Latitude bands are by absolute latitude.')


_T0 = time.time()


def log(stage):
    print(f'{time.time() - _T0:7.1f} s  {stage}', file=sys.stderr, flush=True)


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rnd(x, sig=6):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        x = x.item()
    if isinstance(x, float) and not math.isfinite(x):
        return None                      # JSON keeps no infinities: an open bin edge or a time never reached
    if isinstance(x, bool) or not isinstance(x, float) or x == 0.0:
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


def finite(x):
    """JSON keeps no infinities: a time never reached within the run is None."""
    return None if x is None or not math.isfinite(x) else x


# --- Grid helpers --------------------------------------------------------------------------------------------------

def latlon(cells, cols, ppd):
    row, col = np.divmod(np.asarray(cells, np.int64), cols)
    return 90.0 - (row + 0.5) / ppd, (col + 0.5) / ppd


def arc_m(a, b, cols, ppd):
    """Great-circle distance (m) between cells a and b."""
    la1, lo1 = (np.radians(x) for x in latlon(a, cols, ppd))
    la2, lo2 = (np.radians(x) for x in latlon(b, cols, ppd))
    s = (np.sin(0.5 * (la2 - la1)) ** 2 + np.cos(la1) * np.cos(la2) * np.sin(0.5 * (lo2 - lo1)) ** 2)
    return 2.0 * MOON_RADIUS * np.arcsin(np.sqrt(np.clip(s, 0.0, 1.0)))


def local_relief(surface, lat_deg, cell_km):
    """Highest less lowest surface within RELIEF_HALF_WIDTH_KM of each cell, north-south and east-west (m)."""
    from scipy.ndimage import maximum_filter1d, minimum_filter1d
    rows, cols = surface.shape
    hi = np.empty_like(surface)
    lo = np.empty_like(surface)
    for i in range(rows):
        w = RELIEF_HALF_WIDTH_KM / (cell_km * max(math.cos(math.radians(lat_deg[i])), 1e-3))
        w = int(min(cols // 2 - 1, max(1, round(w))))
        hi[i] = maximum_filter1d(surface[i], 2 * w + 1, mode='wrap')
        lo[i] = minimum_filter1d(surface[i], 2 * w + 1, mode='wrap')
    half = max(1, round(RELIEF_HALF_WIDTH_KM / cell_km))
    hi = maximum_filter1d(hi, 2 * half + 1, axis=0, mode='nearest')
    lo = minimum_filter1d(lo, 2 * half + 1, axis=0, mode='nearest')
    return hi - lo


def carry(r, source, passes):
    """Carry a flux down the routing's receivers in its own topological order. A depression's floor passes on the
    share passes[node] of what reaches it: all of it from a spilling lake, none from a closed depression."""
    final, outlet = r['receiver'], r['outlet']
    n = final.size
    node_of = np.full(n, -1, np.int64)
    node_of[r['sinks']] = np.arange(1, r['sinks'].size + 1)
    acc = source.astype(np.float64).copy()
    donor = ~outlet
    indeg = np.bincount(final[donor], minlength=n).astype(np.int64)
    front = np.flatnonzero((indeg == 0) & donor)
    while front.size:
        out = acc[front].copy()
        k = node_of[front]
        at = k > 0
        out[at] *= passes[k[at]]
        dest = final[front]
        np.add.at(acc, dest, out)
        np.subtract.at(indeg, dest, 1)
        cand = np.unique(dest)
        front = cand[(indeg[cand] == 0) & donor[cand]]
    return acc


# --- Units -----------------------------------------------------------------------------------------------------------

def units(r, shape, sea_level):
    """The units that catch sediment, each cell's first unit downstream, and each unit's parent."""
    from geography import drainage as gd
    from resources import sea_phosphorus_waters as spw
    rows, cols = shape
    n = rows * cols
    lake, label, area, outlet, hf = r['lake'], r['label'], r['area'], r['outlet'], r['height']
    final, level = r['receiver'], r['level']
    body, n_body, outflow, volume, tau, spills, census = spw.lake_bodies(r, shape)
    spill_node = r['passed_on'] > 0
    spill_node[0] = False
    n_nodes = level.size
    closed_nodes = np.flatnonzero(~spill_node)
    closed_nodes = closed_nodes[closed_nodes > 0]
    unit_of_node = np.full(n_nodes, -1, np.int64)
    unit_of_node[closed_nodes] = n_body + 1 + np.arange(closed_nodes.size)
    n_units = n_body + 1 + closed_nodes.size
    spilling_cell = lake & spill_node[label]
    closed_lake_cell = lake & ~spill_node[label]
    stop = outlet | lake
    stop[r['sinks'][~spill_node[1:]]] = True       # closed floors; a spilling floor at its pass height passes on
    unit_stop = np.full(n, -1, np.int64)
    unit_stop[spilling_cell] = body[spilling_cell]
    closed_stop = stop & ~outlet & ~spilling_cell
    unit_stop[closed_stop] = unit_of_node[label[closed_stop]]
    unit_stop[outlet] = 0
    if np.any(unit_stop[stop] < 0):
        raise ValueError('A stop cell has no unit')
    idx = np.arange(n, dtype=final.dtype)
    nxt = np.where(stop, idx, final)
    first = gd.terminal(nxt)
    unit_first = unit_stop[first]
    if np.any(unit_first < 0):
        raise ValueError('A cell drains to no unit')

    parent = np.full(n_units, -1, np.int64)
    parent[0] = 0
    out_all = r['discharge'].copy()
    out_all[r['sinks']] = np.maximum(out_all[r['sinks']] + r['lake_net'][1:], 0.0)
    exits = np.flatnonzero((body > 0) & ~outlet & (body[final] != body))
    if exits.size:
        down = unit_first[final[exits]]
        order = np.lexsort((-out_all[exits], body[exits]))
        b = body[exits][order]
        first_of = np.r_[True, b[1:] != b[:-1]]
        parent[b[first_of]] = down[order][first_of]
    no_exit = np.flatnonzero(parent[1:n_body + 1] < 0) + 1
    if no_exit.size:
        raise ValueError(f'{no_exit.size} spilling lakes have no outlet')
    parent[unit_of_node[closed_nodes]] = unit_first[r['spill_to'][closed_nodes]]

    depth = np.full(n_units, -1, np.int64)
    depth[0] = 0
    children = [[] for _ in range(n_units)]
    for u in range(1, n_units):
        children[parent[u]].append(u)
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for y in children[x]:
            depth[y] = depth[x] + 1
            queue.append(y)
    if np.any(depth < 0):
        raise ValueError(f'{int((depth < 0).sum())} units do not reach the common level')

    below = ~outlet & (hf < level[label])
    hollow = np.bincount(label[below], weights=((level[label] - hf) * area)[below], minlength=n_nodes)
    lake_m2 = np.zeros(n_units)
    lake_m3 = np.zeros(n_units)
    capacity = np.zeros(n_units)
    inflow = np.zeros(n_units)
    lake_m2[:n_body + 1] = np.bincount(body[spilling_cell], weights=area[spilling_cell], minlength=n_body + 1)
    lake_m3[:n_body + 1] = volume
    capacity[:n_body + 1] = volume
    inflow[:n_body + 1] = outflow * YEAR_S
    depth_cl = np.where(closed_lake_cell, r['lake_level'] - hf, 0.0)
    cl_m2 = np.bincount(label[closed_lake_cell], weights=area[closed_lake_cell], minlength=n_nodes)
    cl_m3 = np.bincount(label[closed_lake_cell], weights=(depth_cl * area)[closed_lake_cell], minlength=n_nodes)
    uc = unit_of_node[closed_nodes]
    lake_m2[uc] = cl_m2[closed_nodes]
    lake_m3[uc] = cl_m3[closed_nodes]
    capacity[uc] = hollow[closed_nodes]
    lake_m2[0] = lake_m3[0] = capacity[0] = 0.0
    closed = np.zeros(n_units, bool)
    closed[uc] = True

    # Where each unit lies: its lake, or its hollow when dry.
    where = np.full(n, -1, np.int64)
    where[spilling_cell] = body[spilling_cell]
    dry_hollow = below & ~lake & (unit_of_node[label] >= 0)
    where[closed_lake_cell | dry_hollow] = unit_of_node[label[closed_lake_cell | dry_hollow]]
    w = where >= 0
    lat, lon = latlon(np.flatnonzero(w), cols, 16)
    wu, wa = where[w], area[w]
    tot = np.bincount(wu, weights=wa, minlength=n_units)
    lat_c = np.bincount(wu, weights=wa * lat, minlength=n_units) / np.maximum(tot, 1e-30)
    z = np.bincount(wu, weights=wa * np.cos(np.radians(lon)), minlength=n_units) + 1j * np.bincount(
        wu, weights=wa * np.sin(np.radians(lon)), minlength=n_units)
    lon_c = np.degrees(np.angle(z)) % 360.0
    stop_level = np.where(outlet, sea_level, np.where(lake, r['lake_level'], hf))
    kind = np.where(closed, np.where(lake_m2 > 0, 2, 3), 1)
    kind[0] = 0
    return dict(n_body=n_body, n_units=n_units, parent=parent, depth=depth, closed=closed, kind=kind,
                capacity=capacity, lake_m2=lake_m2, lake_m3=lake_m3, inflow_m3_yr=inflow, lat=lat_c, lon=lon_c,
                unit_first=unit_first, first=first, stop_level=stop_level, census=census, tau=tau)


# --- Fill through time ---------------------------------------------------------------------------------------------

def evolve(u, local, case, beta):
    """Each unit's fill through time: rock in from its own land and from upstream units, kept on Brune's curve in
    spilling lakes (all of it in closed units), with accretion on the lake floor, until full; then passed on whole.
    Returns the years each lake is lost and each unit full, the sea's receipt and the lakes' capacity left."""
    n = u['n_units']
    rate = local * case['rock_gt_yr']
    rho = case['dry_density_kg_m3']
    accretion = case['accretion_mm_yr'] * 1e-3 * u['lake_m2']
    left = u['capacity'].copy()
    lake = u['lake_m2'] > 0
    thr = np.minimum(WETLAND_DEPTH_M * u['lake_m2'], 0.5 * u['lake_m3'])
    t_lost = np.where(lake, np.inf, np.nan)
    t_full = np.full(n, np.inf)
    t_full[0] = np.nan
    closed, inflow, parent = u['closed'], u['inflow_m3_yr'], u['parent']
    order = np.argsort(-u['depth'], kind='stable')
    d = u['depth'][order]
    cuts = np.flatnonzero(np.diff(d)) + 1
    groups = [g for g in np.split(order, cuts) if g.size and u['depth'][g[0]] > 0]
    edges = np.concatenate([[0.0], 10 ** (np.arange(0, round(math.log10(END_YEARS) * STEPS_PER_DECADE) + 1)
                                          / STEPS_PER_DECADE)])
    sea, lakes_left = [], []
    lake_cap = u['capacity'][lake].sum()
    for t0, t1 in zip(edges[:-1], edges[1:]):
        dt = t1 - t0
        m_in = rate * dt
        for g in groups:
            m = m_in[g]
            L = left[g]
            full = L <= 0.0
            te = np.where(closed[g], 1.0,
                          np.exp(-beta * inflow[g] / np.maximum(L, 1e-9) / GRAVITY_RATIO))
            te[full] = 0.0
            vol = te * m / rho + accretion[g] * dt * ~full
            frac = np.where(vol > L, L / np.maximum(vol, 1e-300), 1.0)
            frac[full] = 0.0
            kept = te * m * frac
            new = np.where(full, 0.0, np.maximum(L - vol, 0.0))
            th = thr[g]
            cross = ~full & (L > th) & (new <= th) & lake[g]
            t_lost[g[cross]] = t0 + dt * (L - th)[cross] / vol[cross]
            filled = ~full & (new <= 0.0)
            t_full[g[filled]] = t0 + dt * frac[filled]
            left[g] = new
            np.add.at(m_in, parent[g], m - kept)
        sea.append(m_in[0] / dt)
        lakes_left.append(left[lake].sum() / lake_cap)
    sea = np.array(sea)
    return t_lost, t_full, edges, sea, np.array(lakes_left), np.cumsum(sea * np.diff(edges))


def at_edges(edges, series, years):
    """The value at the step that ends at each year (the steps end on every decade)."""
    out = {}
    for y in years:
        i = int(np.argmin(np.abs(np.log10(edges[1:]) - math.log10(y))))
        out[f'{y:g}'] = float(series[i])
    return out


def classes(u):
    a = u['lake_m2'] / 1e6
    names = []
    for lo, hi in zip(AREA_CLASSES_KM2[:-1], AREA_CLASSES_KM2[1:]):
        names.append((f'{lo:g}_to_{hi:g}_km2' if math.isfinite(hi) else f'over_{lo:g}_km2', (a > lo) & (a <= hi)))
    return names


def bands(lat):
    out = []
    for lo, hi in zip(LAT_BANDS_DEG[:-1], LAT_BANDS_DEG[1:]):
        out.append((f'{lo:g}_{hi:g}', (np.abs(lat) >= lo) & (np.abs(lat) < hi if hi < 90 else np.abs(lat) <= hi)))
    return out


def weighted_median(t, w):
    ok = np.isfinite(t) | np.isinf(t)
    t, w = t[ok], w[ok]
    if not t.size or w.sum() <= 0:
        return None
    o = np.argsort(t)
    c = np.cumsum(w[o]) / w.sum()
    return finite(float(t[o][np.searchsorted(c, 0.5)]))


def survival(u, t_lost):
    lake = u['lake_m2'] > 0
    out = {}
    for name, m in [('all', lake)] + [(n, lake & m) for n, m in classes(u)]:
        a = u['lake_m2'][m] / 1e6
        t = t_lost[m]
        out[name] = dict(
            lakes=int(m.sum()), area_km2=float(a.sum()),
            open_count={f'{y:g}': int((t > y).sum()) for y in SURVIVAL_YEARS},
            open_area_km2={f'{y:g}': float(a[t > y].sum()) for y in SURVIVAL_YEARS},
            median_life_yr=weighted_median(t, np.ones(t.size)),
            area_weighted_median_life_yr=weighted_median(t, a))
    by_band = {}
    for name, m in bands(u['lat']):
        mm = lake & m
        a = u['lake_m2'][mm] / 1e6
        t = t_lost[mm]
        by_band[name] = dict(lakes=int(mm.sum()), area_km2=float(a.sum()),
                             open_area_km2={f'{y:g}': float(a[t > y].sum()) for y in SURVIVAL_YEARS},
                             area_weighted_median_life_yr=weighted_median(t, a))
    return out, by_band


# --- The run ---------------------------------------------------------------------------------------------------------

def compute():
    from geography import drainage as gd
    from resources import sea_phosphorus_waters as spw
    from research.studies.carbon_cycle import run as cc
    committed = json.loads((ROOT / DRAINAGE).read_text())
    gates = json.loads((ROOT / GATES).read_text())
    climatology = ROOT / IGNORED[0]
    if sha256(climatology) != committed['sources'][IGNORED[0]]:
        raise ValueError('Climatology does not match the drainage product')
    rates = sorted(float(k.split()[0]) for k in gates['sediment_return'])
    if rates != sorted({c['erosion_m_myr'] for c in CASES.values()}):
        raise ValueError(f'Erosion cases {rates} differ from the conservation study\'s')
    atlas_, _, _, level, h, lat_row, cell_km, runoff_mm, r = gd.compute(climatology)
    summary = gd.summarise(atlas_, level, h, lat_row, r)
    log('routed')
    reproduced = dict(
        discharge_to_sea_m3s=[summary['rivers']['discharge_to_sea_m3s'], committed['rivers']['discharge_to_sea_m3s']],
        lake_area_share=[summary['lakes']['area_share'], committed['lakes']['area_share']],
        lake_volume_km3=[summary['lakes']['volume_km3'], committed['lakes']['volume_km3']],
        lake_bodies=[summary['lakes']['bodies'], committed['lakes']['bodies']])
    for key, (mine, theirs) in reproduced.items():
        if mine != theirs:
            raise ValueError(f'Routing does not reproduce the drainage product: {key} {mine} != {theirs}')
    rows, cols = h.shape
    shape = h.shape
    hf = h.ravel().astype(np.float64)
    r = dict(r, height=hf, shape=shape)
    outlet, lake, area, final = r['outlet'], r['lake'], r['area'], r['receiver']
    n = hf.size

    # Relief, water and erosion. The land's rock comes off in proportion to its local relief (Ahnert 1970) and to
    # the water that moves it, rain to BQART's discharge exponent (Syvitski & Milliman 2007), so that storms on dry
    # land still move a little; the total is the conservation study's. Relief alone and runoff in place of rain are
    # kept as sensitivities of where the rock first settles.
    from geography.drainage import climate_on_grid
    clim = dict(np.load(climatology))
    clim.pop('metadata')
    rain = np.maximum(climate_on_grid(clim, 'pr_mm_day', rows, cols).ravel().astype(np.float64), 0.0) * 365.25
    wet = np.maximum(runoff_mm.ravel().astype(np.float64), 0.0) * 365.25
    surface = np.where(outlet, level, np.where(lake, r['lake_level'], hf)).reshape(shape).astype(np.float32)
    relief = local_relief(surface, lat_row, cell_km).ravel().astype(np.float64)
    del surface
    eroding = ~outlet & ~lake
    mean_relief = float((relief * area)[eroding].sum() / area[eroding].sum())

    def per_gt(weight):
        """kg of rock a year from each cell for each Gt a year from the land."""
        w = np.where(eroding, weight * area, 0.0)
        return w / w.sum() * 1e12
    patterns = dict(relief_and_rain=per_gt(relief * rain ** WATER_EXPONENT), relief_alone=per_gt(relief),
                    relief_and_runoff=per_gt(relief * wet ** WATER_EXPONENT))
    rock = patterns['relief_and_rain']

    u = units(r, shape, level)
    log('units')
    local = np.bincount(u['unit_first'], weights=rock, minlength=u['n_units'])

    # Carrying it back: lift to where it first settles, to the catchment's highest third, and distance.
    cells = np.flatnonzero(eroding)
    mass = rock[cells]
    uf = u['unit_first'][cells]
    stop_h = u['stop_level'][u['first'][cells]]
    lift = np.maximum(hf[cells] - stop_h, 0.0)
    # Each catchment's highest third: a unit's own land, or for land draining straight to the common level, the land
    # reaching it at one river mouth.
    key = np.where(uf == 0, u['n_units'] + u['first'][cells], uf)
    _, catchment = np.unique(key, return_inverse=True)
    o = np.lexsort((-hf[cells], catchment))
    a_sorted = area[cells][o]
    c_sorted = catchment[o]
    start = np.r_[0, np.flatnonzero(c_sorted[1:] != c_sorted[:-1]) + 1]
    size = np.diff(np.r_[start, o.size])
    cs = np.cumsum(a_sorted)
    cum = cs - np.repeat(cs[start] - a_sorted[start], size)
    top = cum <= HEADS_SHARE * np.repeat(np.add.reduceat(a_sorted, start), size)
    top[start] = True
    n_c = int(catchment.max()) + 1
    heads_a = np.bincount(c_sorted[top], weights=a_sorted[top], minlength=n_c)
    heads_h = np.bincount(c_sorted[top], weights=(hf[cells][o] * a_sorted)[top], minlength=n_c) / \
        np.maximum(heads_a, 1e-30)
    lift_heads = np.maximum(heads_h[catchment] - stop_h, 0.0)
    haul = arc_m(cells, u['first'][cells], cols, 16)
    total = float(mass.sum())
    first_kind = u['kind'][uf]
    kinds = ('common_level_waters', 'spilling_lakes', 'closed_lakes', 'dry_hollows')
    carrying = dict(
        mean_lift_m=float((mass * lift).sum() / total),
        mean_lift_to_heads_m=float((mass * lift_heads).sum() / total),
        heads_share_of_land=float(heads_a.sum() / area[eroding].sum()),
        mean_distance_km=float((mass * haul).sum() / total / 1e3),
        tonne_km_yr_per_gt_yr=float((mass / 1e3 * haul / 1e3).sum()),
        first_settles_in={name: {k: float(pat[cells][first_kind == i].sum() / 1e12) for i, k in enumerate(kinds)}
                          for name, pat in patterns.items()})

    # Units' census and fill through time.
    beta = spw.BRUNE_BETA['median']
    log('carrying back')
    lake_u = u['lake_m2'] > 0
    census = dict(
        spilling_lakes=u['n_body'], closed_lakes=int((u['kind'] == 2).sum()), dry_hollows=int((u['kind'] == 3).sum()),
        lake_area_km2=float(u['lake_m2'].sum() / 1e6), lake_water_km3=float(u['lake_m3'].sum() / 1e9),
        capacity_km3=dict(spilling_lakes=float(u['capacity'][u['kind'] == 1].sum() / 1e9),
                          closed_lakes=float(u['capacity'][u['kind'] == 2].sum() / 1e9),
                          dry_hollows=float(u['capacity'][u['kind'] == 3].sum() / 1e9)),
        tree_depth=int(u['depth'].max()),
        drainage_census=u['census'],
        by_class={name: dict(lakes=int((lake_u & m).sum()), area_km2=float(u['lake_m2'][lake_u & m].sum() / 1e6),
                             water_km3=float(u['lake_m3'][lake_u & m].sum() / 1e9),
                             own_rock_share=float(local[lake_u & m].sum() / 1e12))
                  for name, m in classes(u)})
    fills, grids = {}, {}
    for name, case in CASES.items():
        case = dict(case, rock_gt_yr=gates['sediment_return'][f"{case['erosion_m_myr']:g} m/Myr"]['rock_gt_per_year'],
                    pattern=case.get('pattern', 'relief_and_rain'))
        own = np.bincount(u['unit_first'], weights=patterns[case['pattern']], minlength=u['n_units'])
        t_lost, t_full, edges, sea, left, sea_cum = evolve(u, own, case, beta)
        log(f'fill {name}')
        surv, by_band = survival(u, t_lost)
        dry = u['kind'] == 3
        fills[name] = dict(
            case=case, survival=surv, by_band=by_band,
            lake_capacity_left=at_edges(edges, left, DECADES),
            lakes_capacity_one_percent_left_yr=finite(float(edges[1:][np.argmax(left <= 0.01)])
                                                      if np.any(left <= 0.01) else math.inf),
            to_common_level_gt_yr=at_edges(edges, sea / 1e12, DECADES),
            to_common_level_cumulative_gt=at_edges(edges, sea_cum / 1e12, FIFTHS),
            dry_hollows_full={f'{y:g}': int((t_full[dry] <= y).sum()) for y in DECADES})
        grids[f't_lost_{name}'] = t_lost.astype(np.float32)
        grids[f't_full_{name}'] = t_full.astype(np.float32)

    # Upkeep at the start: where the land's rock settles while every lake catches.
    t0_kept = local.copy()
    kept_by_lake = np.zeros(u['n_units'])
    m_in = local.copy()
    order = np.argsort(-u['depth'], kind='stable')
    for x in order:
        if x == 0:
            continue
        te = 1.0 if u['closed'][x] else math.exp(-beta * u['inflow_m3_yr'][x] / max(u['capacity'][x], 1e-9)
                                                    / GRAVITY_RATIO)
        kept_by_lake[x] = te * m_in[x]
        m_in[u['parent'][x]] += (1.0 - te) * m_in[x]
    kept_by_lake[0] = m_in[0]
    k = kept_by_lake[lake_u]
    ks = np.sort(k)[::-1]
    kc = np.cumsum(ks)
    upkeep = dict(
        kept_share=dict(lakes=float(k.sum() / total), dry_hollows=float(kept_by_lake[u['kind'] == 3].sum() / total),
                        common_level_waters=float(kept_by_lake[0] / total)),
        lakes_holding_half_of_it=int(np.searchsorted(kc, 0.5 * kc[-1]) + 1),
        lakes_holding_nine_tenths=int(np.searchsorted(kc, 0.9 * kc[-1]) + 1),
        by_class={name: float(kept_by_lake[lake_u & m].sum() / total) for name, m in classes(u)})
    del t0_kept

    log('fills and upkeep')
    # Channels: slope, discharge times slope, length and hardness.
    ts = climate_on_grid(clim, 'ts_k', rows, cols).ravel().astype(np.float64)
    runoff = np.where(outlet | (~outlet & (hf < r['level'][r['label']])), 0.0,
                      np.maximum(runoff_mm.ravel(), 0.0) * gd.MM_DAY * area)
    b = cc.DESSERT_B
    t_ref = float(np.log((runoff * np.exp(b * (ts - 300.0))).sum() / runoff.sum()) / b + 300.0)
    weight = np.exp(b * (ts - t_ref))
    passes = (r['passed_on'] > 0).astype(np.float64)
    load = carry(r, runoff * weight, passes)
    log('hardness carried')
    q = r['discharge']
    chan = np.flatnonzero((q >= gd.CHANNEL_M3S) & ~outlet & ~lake)
    rcv = final[chan]
    length = arc_m(chan, rcv, cols, 16)
    rcv_level = np.where(outlet[rcv], level, np.where(lake[rcv], r['lake_level'][rcv], hf[rcv]))
    slope = np.maximum((hf[chan] - rcv_level) / np.maximum(length, 1.0), SLOPE_FLOOR)
    qs = q[chan] * slope
    hardness = load[chan] / q[chan]
    lat_ch, _ = latlon(chan, cols, 16)
    s_bin = np.searchsorted(np.array(SLOPE_BINS_PERCENT[1:-1]), 100.0 * slope, side='right')
    q_bin = np.floor(np.log10(qs) / LOG_QS_STEP).astype(np.int64)
    h_bin = np.floor(np.log10(hardness) / LOG_HARDNESS_STEP).astype(np.int64)
    channel_bands = {}
    for name, m in bands(lat_ch):
        key = s_bin[m] * 100000 + (q_bin[m] + 50000)
        uk, inv = np.unique(key, return_inverse=True)
        km = np.bincount(inv, weights=length[m] / 1e3)
        hk, hinv = np.unique(h_bin[m], return_inverse=True)
        hkm = np.bincount(hinv, weights=length[m] / 1e3)
        channel_bands[name] = dict(
            length_km=float(length[m].sum() / 1e3),
            by_slope_and_log_qs=[[int(k // 100000), int(k % 100000 - 50000), float(v)] for k, v in zip(uk, km)],
            by_log_hardness=[[int(k), float(v)] for k, v in zip(hk, hkm)])
    land = eroding & (runoff > 0)
    w_bin = np.floor(np.log10(weight[land]) / LOG_HARDNESS_STEP).astype(np.int64)
    lat_land, _ = latlon(np.flatnonzero(land), cols, 16)
    spring_bands = {}
    for name, m in bands(lat_land):
        wk, winv = np.unique(w_bin[m], return_inverse=True)
        spring_bands[name] = [[int(k), float(v)] for k, v in zip(wk, np.bincount(winv, weights=area[land][m] / 1e6))]
    land_lat, _ = latlon(np.flatnonzero(eroding), cols, 16)
    relief_bands = {}
    for name, m in bands(land_lat):
        a = area[eroding][m]
        rr = relief[eroding][m]
        o2 = np.argsort(rr)
        c2 = np.cumsum(a[o2]) / a.sum()
        relief_bands[name] = dict(land_km2=float(a.sum() / 1e6), mean_relief_m=float((rr * a).sum() / a.sum()),
                                  mean_rain_mm_yr=float((rain[eroding][m] * a).sum() / a.sum()),
                                  rock_share=float(rock[eroding][m].sum() / 1e12),
                                  lowering_m_myr_per_gt_yr=float(rock[eroding][m].sum() / ROCK_DENSITY / a.sum() * 1e6),
                                  relief_percentiles_m={f'p{p}': float(rr[o2][np.searchsorted(c2, p / 100.0)])
                                                        for p in (10, 50, 90)})

    grids.update(
        unit_parent=u['parent'].astype(np.int32), unit_depth=u['depth'].astype(np.int32),
        unit_kind=u['kind'].astype(np.int8), unit_lat=u['lat'].astype(np.float32),
        unit_lon=u['lon'].astype(np.float32), unit_lake_m2=u['lake_m2'].astype(np.float32),
        unit_lake_m3=u['lake_m3'].astype(np.float32), unit_capacity_m3=u['capacity'].astype(np.float32),
        unit_inflow_m3_yr=u['inflow_m3_yr'].astype(np.float32),
        unit_own_rock_kg_yr_per_gt_yr=local.astype(np.float32),
        unit_kept_at_start_kg_yr_per_gt_yr=kept_by_lake.astype(np.float32),
        cell_unit_first=u['unit_first'].astype(np.int32).reshape(shape),
        cell_local_relief_m=np.clip(np.round(relief), 0, 65535).astype(np.uint16).reshape(shape),
        channel_cell=chan.astype(np.int32), channel_discharge_m3s=q[chan].astype(np.float32),
        channel_slope=slope.astype(np.float32), channel_hardness=hardness.astype(np.float32))
    product = dict(
        schema=SCHEMA,
        evidence=EVIDENCE, reading_rule=READING_RULE,
        reproduced_drainage=reproduced,
        grid=dict(pixels_per_degree=16, cell_km=cell_km, rows=rows, cols=cols),
        land=dict(eroding_km2=float(area[eroding].sum() / 1e6), mean_local_relief_m=mean_relief,
                  lowering_m_myr_per_gt_yr=float(1e12 / ROCK_DENSITY / area[eroding].sum() * 1e6),
                  by_band=relief_bands),
        erosion_cases_m_myr=rates,
        carrying_back=carrying,
        units=census,
        fill=fills,
        upkeep_at_start=upkeep,
        channels=dict(total_km=float(length.sum() / 1e3), slope_bins_percent=list(SLOPE_BINS_PERCENT),
                      log_qs_step=LOG_QS_STEP, log_hardness_step=LOG_HARDNESS_STEP, by_band=channel_bands),
        springs=dict(reference_temperature_k=t_ref, log_hardness_step=LOG_HARDNESS_STEP,
                     land_km2_by_log_hardness=spring_bands),
        parameters=dict(rock_density_kg_m3=ROCK_DENSITY, relief_half_width_km=RELIEF_HALF_WIDTH_KM,
                        wetland_depth_m=WETLAND_DEPTH_M, heads_share=HEADS_SHARE, brune_beta=beta,
                        gravity_ratio=GRAVITY_RATIO, accretion_mm_yr=ACCRETION_MM_YR,
                        steps_per_decade=STEPS_PER_DECADE, slope_floor=SLOPE_FLOOR))
    return product, grids


def main(argv=None) -> int:
    product, grids = compute()
    product = rnd(product)
    product['producer'] = dict(
        files={p: sha256(ROOT / p) for p in CODE},
        inputs={p: sha256(ROOT / p) for p in INPUTS},
        ignored_inputs={p: sha256(ROOT / p) for p in IGNORED},
        constants=constants_used([ROOT / p for p in CODE]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    GRIDS.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(GRIDS, metadata=json.dumps(dict(schema=SCHEMA + '-grids', product=str(OUT.relative_to(ROOT)),
                                                        reading_rule=READING_RULE)), **grids)
    print(OUT.relative_to(ROOT))
    print(GRIDS)
    return 0


if __name__ == '__main__':
    sys.exit(main())
