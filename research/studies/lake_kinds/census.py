"""The census's lake set: geography's rain-fed lakes as the waters' oxygen study labels them, with each body's floor
by depth, its outlet, its soda path, the groundwater it gains and the springs' water against the lakes'.

    python -m research.studies.lake_kinds.census      # about fifteen seconds and 0.6 GB: results/lake_census.json

Nothing here runs the lake model: the census is geometry on LOLA, hydrology from geography's products and chemistry
from the sea chemistry and springs studies, so it holds whatever the model finds. run.py sorts its bodies into kinds
with the depths the lake model sets.

- Bodies: lake_oxygen's lakes.py labels the drainage product's lake cells on LOLA's 16 pixel/degree grid as
  eight-connected bodies, the drainage summary's 27,865; each cell's depth is its lake's level less LOLA's height above
  the GRAIL geoid, and its area is the cell's.
- Outlets: the drainage product marks each lake cell whose depression does not spill. The drainage summary counts a
  body holding any such cell as closed, its 2,495; 2,287 of them lie mostly in depressions that do not spill, and the
  other 208 are open lakes in which nested pits keep a lower level of their own, the outer lake's water left off
  them by the routing's one level for each depression. The census takes the mostly closed bodies as the closed
  lakes, their closed cells a lake of their own, and the others as open.
- Evaporation and rain: geography's drainage takes each lake's open-water evaporation as the climate run's zonal mean
  over its seas at the lake's latitude, and the run's rain at the cell (A28_dim5_moon, years 20-29, the run it routed).
- Soda: each closed part runs the sea chemistry's left-alone path of the fresh seas at its own turnover, the sea
  chemistry's scaling (its depth over its evaporation, against the waters' 1,908 years) and a balance variant (the
  inflow the drainage balances, the evaporation less the rain on it, carrying the land water's load).
- Groundwater: stage 2's mature water table at the central crust (geography/groundwater_table.py) re-solved by its own
  solver from its saved heads, which must give its product's exchange with the lakes; the flow from land into each of
  its 1,061 lakes on the 15 km grid, per area of that lake, is what each body whose cells fall in it gains. A body too
  small to appear on the 15 km grid gains what its band's 15 km lakes gain for each metre of shore, the median over
  the gaining ones, along its own shore: the inflow through a lake's bed is the aquifer's flow towards it across its
  shoreline, so per area it grows as a lake gets smaller. Shores are taken as circles of the lake's area in both.
- Springs: the springs study's mature, central springs by band against the held lakes' water (sea chemistry), at the
  lake model's linear salinity term and Lake Kivu's coefficient for dissolved CO2.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'lake_census.json'
SCHEMA = 'terluna.research.lake-kinds-census/1'
CLIMATOLOGY = ROOT / 'climate' / 'gcm' / 'products' / 'climatology_A28_dim5_moon.npz'
SEA_CHEMISTRY = ROOT / 'research' / 'studies' / 'sea_chemistry' / 'results' / 'sea_chemistry.json'
SPRINGS = ROOT / 'research' / 'studies' / 'springs' / 'results' / 'springs.json'
HELD_FRESH = dict(corner='central', delivery='fresh', sinks='managed', sediment='kept_open')    # lake_oxygen's
GROUNDWATER_GRIDS = ROOT / 'research' / 'runs' / 'resources' / 'groundwater' / 'groundwater_table_2ppd.npz'
GROUNDWATER = ROOT / 'geography' / 'results' / 'groundwater_table.json'
DRAINAGE_SUMMARY = ROOT / 'geography' / 'results' / 'drainage.json'
PPD = 16
COARSE = 8                                   # 16 pixel/degree cells to a 2 pixel/degree cell
THRESHOLDS_M = (2.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0)
SIZE_EDGES_KM2 = (0.0, 1.0, 10.0, 100.0, 1000.0, 1e4, np.inf)     # lake_oxygen's lakes.py classes
BAND_EDGES_DEG = (0.0, 15.0, 30.0, 45.0, 90.01)
DEEPEST_EDGES_M = (0.0, 2.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0,
                   100.0, 150.0, 250.0, 352.0, 500.0, 700.0, 1000.0, 1500.0, 1925.0, 2500.0, 3500.0, 5000.0, 1e5)
FLOOR_EDGES_M = (0.0, 2.0, 5.0, 10.0, 20.0, 40.0, 70.0, 100.0, 150.0, 250.0, 400.0, 700.0, 1e5)   # lake_anoxia's
YEAR_S = 365.25 * 86400.0
SALINITY_DENSITY = 0.8                       # kg/m3 per g/kg, the lake model's linear salinity term (lake.py)
SPRING_TEMPERATURES_C = (22.0, 25.0, 28.0)   # about the deep water's range at the resort
MARGINS_KG_M3 = (0.003, 0.03)                # lake_oxygen's solute brackets: CO2's 0.003 already moved its deep lake's
                                             # mixing, and calcite's 0.03 shortened its deeper nights' from 55 to 20 m
FRESHENING = (0.01, 0.05)                    # a wet spell's freshening of a closed lake's mixed layer (see README)
RENEWAL_DEPTHS_M = (10.0, 12.0, 14.0, 16.0, 20.0, 25.0)    # floors below which the springs' share is followed
MIN_DEPTH_M = 0.1                            # lake_anoxia's floor on depth: a closed pan's level sits on its top cell
MIN_NET_EVAPORATION_M_YR = 0.01              # a closed lake balanced by rain alone still takes its catchment's trickle
MOSTLY = 0.5                                 # a body is a closed lake when this share of it lies in closed depressions


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def band_name(i):
    lo, hi = BAND_EDGES_DEG[i], BAND_EDGES_DEG[i + 1]
    return f'{lo:g}-{min(hi, 90.0):g}'


def size_name(i):
    lo, hi = SIZE_EDGES_KM2[i], SIZE_EDGES_KM2[i + 1]
    return f'{lo:g}-{hi:g} km2' if np.isfinite(hi) else f'over {lo:g} km2'


def lake_set():
    """Every lake cell and body: lake_oxygen's labelling, with the drainage product's outlets."""
    from research.studies.lake_oxygen import lakes as K
    d = K.lake_depths()
    b = K.bodies(d)
    z = np.load(K.DRAINAGE)
    if not np.array_equal(z['lake_cell'], d['cells']):
        raise SystemExit('the drainage product\'s lake cells are not the ones lake_oxygen labels')
    d['spills'] = z['lake_spills'].astype(bool)
    d['depth'] = np.maximum(d['depth'], 0.0)             # 1,420 cells sit a fraction of a millimetre above their level
    return d, b


def bodies(d, b):
    """Per-body table (index k for body k + 1): size, depth, place, outlet, and the floor within each threshold."""
    n = b['n']
    body, a, depth = b['body'], d['area'], d['depth']
    count = np.bincount(body, None, n + 1)[1:]
    closed_cells = np.bincount(body, (~d['spills']).astype(float), n + 1)[1:]
    closed_area = np.bincount(body, a * (~d['spills']), n + 1)[1:]
    area = b['area'][1:]
    lat = b['lat'][1:]
    out = dict(n=n, cells=count, area_km2=area / 1e6, volume_km3=b['volume'][1:] / 1e9,
               mean_depth_m=b['volume'][1:] / area, deepest_m=b['deepest'][1:], lat=lat, lon=b['lon'][1:],
               level_above_sea_m=b['level'][1:] - d['sea_level'], any_closed=closed_cells > 0,
               closed=closed_area / area >= MOSTLY, closed_area_share=closed_area / area,
               band=np.clip(np.searchsorted(BAND_EDGES_DEG, np.abs(lat), side='right') - 1, 0, len(BAND_EDGES_DEG) - 2),
               size=np.clip(np.searchsorted(SIZE_EDGES_KM2, area / 1e6, side='right') - 1, 0,
                            len(SIZE_EDGES_KM2) - 2),
               floor_within_km2={}, volume_below_km3={})
    for t in THRESHOLDS_M:
        out['floor_within_km2'][t] = np.bincount(body, a * (depth <= t), n + 1)[1:] / 1e6
        out['volume_below_km3'][t] = np.bincount(body, a * np.maximum(depth - t, 0.0), n + 1)[1:] / 1e9
    # the closed part of each body: its cells in depressions that do not spill, a lake of their own
    shut = ~d['spills']
    deepest_closed = np.zeros(n + 1)
    np.maximum.at(deepest_closed, body[shut], depth[shut])
    out['closed_area_km2'] = closed_area / 1e6
    out['closed_volume_km3'] = np.bincount(body, a * depth * shut, n + 1)[1:] / 1e9
    out['closed_deepest_m'] = deepest_closed[1:]
    return out


def floor_within(t_by_cell, d, b):
    """Floor area (km2) of each body within a depth that may differ from cell to cell, and the volume below it."""
    n = b['n']
    a, depth = d['area'], d['depth']
    return (np.bincount(b['body'], a * (depth <= t_by_cell), n + 1)[1:] / 1e6,
            np.bincount(b['body'], a * np.maximum(depth - t_by_cell, 0.0), n + 1)[1:] / 1e9)


def evaporation(d, b):
    """Each body's open-water evaporation and the rain on it (m a year, area means), as geography's drainage takes
    them from the climate run it routed."""
    from geography import drainage as G
    clim = dict(np.load(CLIMATOLOGY))
    clim.pop('metadata')
    n = b['n']
    a = d['area']
    evap = G.lake_evaporation(clim, d['lat']) * 365.25 / 1000.0
    rows, cols = 180 * PPD, 360 * PPD
    rain_grid = G.climate_on_grid(clim, 'pr_mm_day', rows, cols)
    rain = rain_grid.ravel()[d['cells']].astype(float) * 365.25 / 1000.0
    del rain_grid
    w = np.bincount(b['body'], a, n + 1)[1:]
    return dict(evaporation_m_yr=np.bincount(b['body'], a * evap, n + 1)[1:] / w,
                rain_m_yr=np.bincount(b['body'], a * rain, n + 1)[1:] / w)


def groundwater(d, b):
    """Stage 2's mature central water table, re-solved from its saved heads, and the groundwater each body gains."""
    from geography import groundwater_table as gt
    from geography import groundwater_uptake as gu
    from research.studies.springs import run as SP
    s = gt.coarsen(gt.setting())
    stage1 = json.loads((ROOT / gt.INPUTS[1]).read_text())
    p = stage1['cases']['central']['params']
    rech = gt.recharge(s, p['recharge_share'], gu.ARID_RECHARGE_SHARE)
    grids = np.load(GROUNDWATER_GRIDS)
    land = ~(s['sea'] | s['lake'])
    depth = grids['depth_m_mature'].astype(float)
    h0 = np.where(land & np.isfinite(depth), s['ground'] - np.nan_to_num(depth), s['ground'])
    r = gt.solve(s, rech, p['perm_1km'], gu.PERM_CAP_M2, gu.PERM_EXPONENT, h0)
    lab, n_lakes = gt.lakes_of(s)
    summary = gt.summarise(s, r, rech, lab, n_lakes)
    product = json.loads(GROUNDWATER.read_text())['cases']['central']['states']['mature']
    check = dict(labels_as_saved=bool(np.array_equal(lab, grids['lake_label'])),
                 iterations=r['iterations'], into_lakes_km3_yr=summary['into_lakes_km3_yr'],
                 product_into_lakes_km3_yr=product['into_lakes_km3_yr'],
                 gaining=summary['lakes']['gaining'], product_gaining=product['lakes']['gaining'],
                 losing=summary['lakes']['losing'], product_losing=product['lakes']['losing'],
                 gain_p50_m_yr=summary['lakes']['gain_m_yr_percentiles']['p50'],
                 product_gain_p50_m_yr=product['lakes']['gain_m_yr_percentiles']['p50'])
    area = gt.geometry(s)[3]
    a_face, b_face, flow = SP.faces(s, r['h'], p['perm_1km'])        # m3/s from a to b
    lake_cell = s['lake'].ravel()
    free = land.ravel()
    gross = np.zeros(lake_cell.size)
    m = free[a_face] & lake_cell[b_face] & (flow > 0)
    np.add.at(gross, b_face[m], flow[m])
    m = free[b_face] & lake_cell[a_face] & (flow < 0)
    np.add.at(gross, a_face[m], -flow[m])
    lake_area = np.bincount(lab.ravel(), (area * s['lake']).ravel(), n_lakes + 1)
    gross_rate = np.bincount(lab.ravel(), gross, n_lakes + 1) / np.maximum(lake_area, 1.0) * YEAR_S   # m/yr
    net_rate = np.bincount(lab.ravel(), r['into_fixed'].ravel(), n_lakes + 1) / np.maximum(lake_area, 1.0) * YEAR_S
    check['gross_into_lakes_km3_yr'] = float(gross.sum() * YEAR_S / 1e9)
    # each 16 pixel/degree lake cell's 15 km cell
    r16, c16 = np.divmod(d['cells'], 360 * PPD)
    idx = (r16 // COARSE) * s['cols'] + (c16 // COARSE)
    cell_lake = lab.ravel()[idx]
    a = d['area']
    nb = b['n']
    on_grid = np.bincount(b['body'], a * (cell_lake > 0), nb + 1)[1:]
    total = np.bincount(b['body'], a, nb + 1)[1:]
    gain = np.bincount(b['body'], a * gross_rate[cell_lake] * (cell_lake > 0), nb + 1)[1:] / np.maximum(on_grid, 1.0)
    net = np.bincount(b['body'], a * net_rate[cell_lake] * (cell_lake > 0), nb + 1)[1:] / np.maximum(on_grid, 1.0)
    # bodies the 15 km grid does not hold: their band's inflow per metre of shore, along their own shore
    lat2 = np.repeat(s['lat'][:, None], s['cols'], axis=1)
    lake_lat = np.bincount(lab.ravel(), (np.abs(lat2) * area * s['lake']).ravel(), n_lakes + 1) / np.maximum(lake_area,
                                                                                                         1.0)
    shore = 2.0 * np.sqrt(np.pi * lake_area)                                         # m, a circle of the lake's area
    per_shore = gross_rate * lake_area / np.maximum(shore, 1.0)                      # m3 a year per m of shore
    lake_band = np.clip(np.searchsorted(BAND_EDGES_DEG, lake_lat, side='right') - 1, 0, len(BAND_EDGES_DEG) - 2)
    gaining = (gross_rate > 0.001) & (np.arange(n_lakes + 1) > 0)
    band_shore = np.array([float(np.median(per_shore[gaining & (lake_band == i)]))
                           if (gaining & (lake_band == i)).any() else 0.0 for i in range(len(BAND_EDGES_DEG) - 1)])
    lat_b = np.abs(b['lat'][1:])
    body_band = np.clip(np.searchsorted(BAND_EDGES_DEG, lat_b, side='right') - 1, 0, len(BAND_EDGES_DEG) - 2)
    body_area = total
    by_shore = band_shore[body_band] * 2.0 * np.sqrt(np.pi * body_area) / np.maximum(body_area, 1.0)    # m/yr
    seeping = grids['seepage_mm_yr_mature'].astype(float).ravel()[idx] > 0.0
    near_seep = np.bincount(b['body'], a * seeping, nb + 1)[1:] / np.maximum(total, 1.0)
    resolved = on_grid / total >= 0.5
    check['band_inflow_per_shore_m3_yr_m'] = {band_name(i): float(band_shore[i]) for i in range(band_shore.size)}
    return dict(resolved=resolved, on_grid_share=on_grid / total, gain_m_yr=np.where(resolved, gain, np.nan),
                net_m_yr=np.where(resolved, net, np.nan), by_shore_m_yr=np.where(resolved, np.nan, by_shore),
                seeping_around_share=near_seep, check=check, lakes_on_grid=int(n_lakes),
                lake_gross_m_yr=gross_rate[1:], lake_net_m_yr=net_rate[1:], lake_area_km2=lake_area[1:] / 1e6)


def inputs():
    """The files the census reads, with their sha256."""
    from research.studies.lake_oxygen import lakes as K
    from geography import fetch_inputs
    paths = [K.DRAINAGE, K.DRAINAGE_SUMMARY, K.ATLAS, CLIMATOLOGY, GROUNDWATER, GROUNDWATER_GRIDS,
             ROOT / 'geography' / 'results' / 'groundwater_uptake.json', ROOT / 'geography' / 'products' /
             'atlas_28pct_4ppd.npz']
    out = {str(p.relative_to(ROOT)): sha256(p) for p in paths}
    for name in ('ldem_16.img', 'ldem_4.img', 'jggrx_0420a_sha.tab'):
        out[f'geography/inputs/{name}'] = sha256(fetch_inputs.path(name))
    return out


# --- soda -------------------------------------------------------------------------------------------------------

def weighted_median(x, w):
    """The w-weighted median of x."""
    i = np.argsort(x)
    c = np.cumsum(w[i])
    return float(x[i][np.searchsorted(c, 0.5 * c[-1])])


def soda_paths():
    """The sea chemistry's left-alone fresh seas (natural sinks, kept open) by weathering corner: the milestones and
    the central path, with the waters' turnover time its closed-lake scaling divides by."""
    chem = json.loads(SEA_CHEMISTRY.read_text())
    s = chem['setting']
    turnover = (s['sea_kg'] + s['lake_kg']) / s['discharge_kg_yr']
    runs = {r['corner']: r for r in chem['runs']
            if r['delivery'] == 'fresh' and r['sinks'] == 'natural' and r['sediment'] == 'kept_open'}
    milestones = {c: dict(one_g_kg=r['milestones_yr']['fresh_limit'], earth_salinity=r['milestones_yr'][
        'earth_salinity'], brine=r['milestones_yr']['brine_ceiling']) for c, r in runs.items()}
    path = [(row['year'], row['salinity_g_kg']) for row in runs['central']['rows'] if row['salinity_g_kg'] > 0]
    return turnover, milestones, path


def years_to(path, salinity):
    """The central left-alone path's year at a salinity, interpolated in logarithms."""
    years = np.array([p[0] for p in path if p[0] > 0])
    sal = np.maximum.accumulate(np.array([p[1] for p in path if p[0] > 0]))
    return float(10 ** np.interp(np.log10(salinity), np.log10(sal), np.log10(years)))


def closed_parts(d, b, B, springs):
    """Each body's closed cells as a lake of their own: area, depth, deepest point, evaporation and rain, and its
    turnover against the waters' by the sea chemistry's scaling and by the drainage balance."""
    from geography import drainage as Gd
    turnover, milestones, path = soda_paths()
    n = B['n']
    shut = ~d['spills']
    a = d['area'] * shut
    ca = np.bincount(b['body'], a, n + 1)[1:]
    clim = dict(np.load(CLIMATOLOGY))
    clim.pop('metadata')
    evap_cell = Gd.lake_evaporation(clim, d['lat']) * 365.25 / 1000.0
    rain_grid = Gd.climate_on_grid(clim, 'pr_mm_day', 180 * PPD, 360 * PPD)
    rain_cell = rain_grid.ravel()[d['cells']].astype(float) * 365.25 / 1000.0
    del rain_grid
    evap = np.bincount(b['body'], a * evap_cell, n + 1)[1:] / np.maximum(ca, 1.0)
    rain = np.bincount(b['body'], a * rain_cell, n + 1)[1:] / np.maximum(ca, 1.0)
    has = (ca > 0) & B['closed']
    area = B['closed_area_km2']
    depth = np.maximum(B['closed_volume_km3'] * 1e3 / np.maximum(area, 1e-12), MIN_DEPTH_M)
    w = springs['moon']['weathering']
    load_ratio = w['rivers_to_sea_alkalinity_meq_kg'] / w['land_water_alkalinity_meq_kg']
    net = np.maximum(evap - rain, MIN_NET_EVAPORATION_M_YR)
    return dict(has=has, area_km2=area, depth_m=depth, deepest_m=B['closed_deepest_m'], evaporation_m_yr=evap,
                rain_m_yr=rain, by_scaling=depth / np.maximum(evap, 1e-3) / turnover,
                by_balance=depth / net / turnover * load_ratio, turnover=turnover, milestones=milestones, path=path,
                load_ratio=load_ratio)


def soda_table(cp):
    """The closed parts by their deepest point: how long they take to 1 g/kg and to Earth's salinity, and to the
    salinity at which a wet spell's freshening of the mixed layer outweighs a margin."""
    has, area = cp['has'], cp['area_km2']
    s_star = {f'{f:g}/{m:g}': m / (SALINITY_DENSITY * f) for f in FRESHENING for m in MARGINS_KG_M3}
    out = []
    for lo, hi in zip(DEEPEST_EDGES_M[:-1], DEEPEST_EDGES_M[1:]):
        m = has & (cp['deepest_m'] > lo) & (cp['deepest_m'] <= hi) if lo > 0 else has & (cp['deepest_m'] <= hi)
        if not m.any():
            continue
        w = area[m]
        row = dict(deepest_m=[lo, hi], parts=int(m.sum()), area_km2=float(w.sum()),
                   mean_depth_m=float((cp['depth_m'][m] * w).sum() / w.sum()))
        for label in ('by_scaling', 'by_balance'):
            f = cp[label][m]
            row[label] = dict(
                years_to_1_g_kg={c: weighted_median(f * cp['milestones'][c]['one_g_kg'], w)
                                 for c in ('high', 'central', 'low')},
                years_to_earth_salinity={c: weighted_median(f * cp['milestones'][c]['earth_salinity'], w)
                                         for c in ('high', 'central', 'low')},
                area_share_past_1_g_kg_within={f'{y:g}': float(w[f * cp['milestones']['central']['one_g_kg'] <= y]
                                                              .sum() / w.sum()) for y in (1e3, 1e4, 1e5, 1e6)},
                years_to_hold_a_chemocline={k: weighted_median(f * years_to(cp['path'], v), w)
                                            for k, v in s_star.items()})
        out.append(row)
    return out, s_star


# --- springs ----------------------------------------------------------------------------------------------------

def spring_water(springs, lake_g_kg, si_per_eq, temps=SPRING_TEMPERATURES_C):
    """How much denser the springs' water is than the lakes' at the same temperature, by band, without and with the
    weathering's silica, and the warmth by which it would have to be cooler than the deep water to stay denser."""
    from biosphere.ecology import lake as L
    from research.studies.lake_anoxia import gases as G
    from research.studies.sea_chemistry import run as SC
    out = {}
    for band, x in springs['by_band'].items():
        if band not in ('0-15', '15-30', '30-45'):
            continue
        alk = x['alkalinity_meq_kg'] / 1000.0
        ca, mg = x['calcium_mmol_kg'] / 1000.0, x['magnesium_mmol_kg'] / 1000.0
        mono = max(alk - 2.0 * (ca + mg), 0.0)
        ions = ca * SC.cc.M_CA + mg * SC.cc.M_MG + mono * SC.cc.M_NA + alk * SC.M_HCO3            # g/kg
        silica = si_per_eq * alk * SC.M_SIO2
        row = dict(salinity_g_kg=[ions, ions + silica], lake_salinity_g_kg=lake_g_kg, at={})
        for t in temps:
            co2 = (x['pco2_uatm'] - springs['air_uatm']) * 1e-6 * float(G.co2_kh(t)) * 44.01        # g/kg
            low = SALINITY_DENSITY * (ions - lake_g_kg) + L.RHO0 * G.BETA_CO2_PER_G_KG * co2
            high = SALINITY_DENSITY * (ions + silica - lake_g_kg) + L.RHO0 * G.BETA_CO2_PER_G_KG * co2
            alpha = L._expansion(t)
            row['at'][f'{t:g}'] = dict(excess_co2_g_kg=co2, density_excess_kg_m3=[low, high],
                                       warmth_that_cancels_it_k=[low / (L.RHO0 * alpha), high / (L.RHO0 * alpha)])
        out[band] = row
    return out


def renewal(B, GW, d, b, excess):
    """What the groundwater does to the water below a floor depth: the share of that water it brings in a year if
    all its inflow across the floor below the depth stays there, and the years it takes to make it denser than the
    water above by each margin at the springs' least density excess."""
    gain = np.where(GW['resolved'], np.nan_to_num(GW['gain_m_yr']), np.nan_to_num(GW['by_shore_m_yr']))
    area = B['area_km2']
    out = {}
    for z in RENEWAL_DEPTHS_M:
        floor_below = area - floor_within(np.full(d['depth'].shape, z), d, b)[0]
        below = B['volume_below_km3'][z] if z in B['volume_below_km3'] else floor_within(
            np.full(d['depth'].shape, z), d, b)[1]
        share = np.where(below > 0, gain * floor_below / 1000.0 / np.maximum(below, 1e-12), 0.0)    # a year
        row = []
        for lo, hi in zip(DEEPEST_EDGES_M[:-1], DEEPEST_EDGES_M[1:]):
            m = (B['deepest_m'] > max(lo, z)) & (B['deepest_m'] <= hi)
            if not m.any():
                continue
            years = {f'{mg:g}': np.where(share[m] > 0, mg / np.maximum(share[m] * excess, 1e-30), np.inf)
                     for mg in MARGINS_KG_M3}
            row.append(dict(deepest_m=[lo, hi], bodies=int(m.sum()), area_km2=float(area[m].sum()),
                            renewed_share_a_year=[float(x) for x in np.quantile(share[m], (0.1, 0.5, 0.9))],
                            years_to_margin={k: dict(quantiles=[float(x) for x in np.quantile(v, (0.1, 0.5, 0.9))],
                                                     within_1_yr=int((v <= 1.0).sum()),
                                                     within_10_yr=int((v <= 10.0).sum()),
                                                     within_50_yr=int((v <= 50.0).sum()))
                                             for k, v in years.items()}))
        out[f'{z:g}'] = row
    return out


# --- the tables -------------------------------------------------------------------------------------------------

def tables(B, d, b, GW, cp):
    """Bodies, area and floors by deepest point, band and size; the cumulative count and area by deepest point."""
    area, band, size, deepest = B['area_km2'], B['band'], B['size'], B['deepest_m']
    n_bands, n_sizes = len(BAND_EDGES_DEG) - 1, len(SIZE_EDGES_KM2) - 1
    by_deepest = []
    for lo, hi in zip(DEEPEST_EDGES_M[:-1], DEEPEST_EDGES_M[1:]):
        m = (deepest > lo) & (deepest <= hi) if lo > 0 else deepest <= hi
        cells = m[b['body'] - 1]
        floor = np.histogram(d['depth'][cells], FLOOR_EDGES_M, weights=d['area'][cells])[0] / 1e6
        gain = np.where(GW['resolved'], GW['gain_m_yr'], GW['by_shore_m_yr'])
        by_deepest.append(dict(
            deepest_m=[lo, hi], bodies=int(m.sum()), area_km2=float(area[m].sum()),
            volume_km3=float(B['volume_km3'][m].sum()), closed_bodies=int((m & B['closed']).sum()),
            closed_part_area_km2=float(B['closed_area_km2'][m].sum()),
            floor_within_km2={f'{t:g}': float(B['floor_within_km2'][t][m].sum()) for t in THRESHOLDS_M},
            volume_below_km3={f'{t:g}': float(B['volume_below_km3'][t][m].sum()) for t in THRESHOLDS_M},
            floor_by_local_depth_km2=[float(x) for x in floor],
            by_band={band_name(i): dict(bodies=int((m & (band == i)).sum()),
                                        area_km2=float(area[m & (band == i)].sum())) for i in range(n_bands)},
            by_size={size_name(i): dict(bodies=int((m & (size == i)).sum()),
                                        area_km2=float(area[m & (size == i)].sum())) for i in range(n_sizes)},
            groundwater=dict(resolved_bodies=int((m & GW['resolved']).sum()),
                             gain_m_yr=[float(x) for x in np.quantile(gain[m], (0.1, 0.5, 0.9))] if m.any() else None,
                             seeping_around=int((m & ~GW['resolved'] & (GW['seeping_around_share'] > 0.5)).sum()))))
    cumulative = {}
    for i in range(n_bands + 1):
        mb = np.ones(area.shape, bool) if i == n_bands else band == i
        cumulative['all' if i == n_bands else band_name(i)] = {
            f'{e:g}': [int((mb & (deepest <= e)).sum()), float(area[mb & (deepest <= e)].sum())]
            for e in DEEPEST_EDGES_M[1:-1]}
    by_size = {size_name(i): dict(bodies=int((size == i).sum()), area_km2=float(area[size == i].sum()),
                                  deepest_m=[float(x) for x in np.quantile(deepest[size == i], (0.1, 0.5, 0.9))])
               for i in range(n_sizes)}
    by_band = {band_name(i): dict(bodies=int((band == i).sum()), area_km2=float(area[band == i].sum()),
                                  closed_bodies=int(((band == i) & B['closed']).sum()))
               for i in range(n_bands)}
    largest = np.argsort(area)[::-1][:12]
    return dict(by_deepest=by_deepest, cumulative_bodies_and_area_by_deepest=cumulative, by_size=by_size,
                by_band=by_band, largest=[dict(area_km2=float(area[k]), deepest_m=float(deepest[k]),
                                               mean_depth_m=float(B['mean_depth_m'][k]), lat=float(B['lat'][k]),
                                               lon=float(B['lon'][k]), closed=bool(B['closed'][k]))
                                          for k in largest])


def summarise(B, d, b, E, GW, cp, springs, chem):
    held = next(r for r in chem['runs'] if all(r[k] == v for k, v in HELD_FRESH.items()))
    water = spring_water(springs, float(held['rows'][-1]['lake_salinity_g_kg']), chem['setting']['si_per_eq'])
    excess = min(min(v['density_excess_kg_m3'][0] for v in row['at'].values()) for row in water.values())
    soda, s_star = soda_table(cp)
    drainage = json.loads(DRAINAGE_SUMMARY.read_text())['lakes']
    has = cp['has']
    nested = B['any_closed'] & ~B['closed']
    total = dict(bodies=B['n'], closed_bodies=int(B['closed'].sum()), area_km2=float(B['area_km2'].sum()),
                 volume_km3=float(B['volume_km3'].sum()),
                 closed_bodies_area_km2=float(B['area_km2'][B['closed']].sum()),
                 closed_part_area_km2=float(cp['area_km2'][has].sum()),
                 bodies_with_any_closed_cell=int(B['any_closed'].sum()),
                 their_area_km2=float(B['area_km2'][B['any_closed']].sum()),
                 wholly_closed_bodies=int((B['closed_area_share'] >= 0.999).sum()),
                 open_bodies_with_nested_closed_pits=dict(bodies=int(nested.sum()),
                                                          area_km2=float(B['area_km2'][nested].sum()),
                                                          pits_km2=float(B['closed_area_km2'][nested].sum())),
                 drainage_summary=dict(bodies=drainage['bodies'], closed_bodies=drainage['closed_bodies'],
                                       volume_km3=drainage['volume_km3'], median_depth_m=drainage['median_depth_m']),
                 median_local_depth_m=float(np.median(d['depth'])))
    gain = np.where(GW['resolved'], GW['gain_m_yr'], GW['by_shore_m_yr'])
    return dict(
        total=total, **tables(B, d, b, GW, cp),
        closed=dict(parts=int(has.sum()), area_km2=float(cp['area_km2'][has].sum()),
                    mean_depth_m=float(B['closed_volume_km3'][has].sum() * 1e3 / cp['area_km2'][has].sum()),
                    evaporation_m_yr=weighted_median(cp['evaporation_m_yr'][has], cp['area_km2'][has]),
                    rain_m_yr=weighted_median(cp['rain_m_yr'][has], cp['area_km2'][has]),
                    net_evaporation_m_yr=[float(x) for x in np.quantile((cp['evaporation_m_yr'] - cp['rain_m_yr'])[has],
                                                                        (0.1, 0.5, 0.9))],
                    turnover_of_the_waters_yr=cp['turnover'], land_water_over_rivers_load=1.0 / cp['load_ratio'],
                    milestones_of_the_left_alone_seas_yr=cp['milestones'],
                    salinity_holding_a_chemocline_g_kg=s_star, by_deepest=soda),
        groundwater=dict(GW['check'], resolved_bodies=int(GW['resolved'].sum()),
                         resolved_area_km2=float(B['area_km2'][GW['resolved']].sum()),
                         gain_m_yr_resolved=[float(x) for x in np.nanquantile(GW['gain_m_yr'], (0.1, 0.5, 0.9))],
                         gain_m_yr_by_shore=[float(x) for x in np.nanquantile(GW['by_shore_m_yr'], (0.1, 0.5, 0.9))],
                         unresolved_in_seeping_cells=int((~GW['resolved'] & (GW['seeping_around_share'] > 0.5)).sum()),
                         inflow_km3_yr=dict(resolved=float((np.nan_to_num(GW['gain_m_yr']) * B['area_km2']).sum()
                                                           / 1000.0),
                                            by_shore=float((np.nan_to_num(GW['by_shore_m_yr']) * B['area_km2']).sum()
                                                           / 1000.0)),
                         renewal_below=renewal(B, GW, d, b, excess), spring_density_excess_used_kg_m3=excess),
        spring_water=water,
        evaporation_m_yr=dict(open=weighted_median(E['evaporation_m_yr'][~B['closed']], B['area_km2'][~B['closed']]),
                              closed=weighted_median(E['evaporation_m_yr'][B['closed']], B['area_km2'][B['closed']])),
        rain_m_yr=dict(open=weighted_median(E['rain_m_yr'][~B['closed']], B['area_km2'][~B['closed']]),
                       closed=weighted_median(E['rain_m_yr'][B['closed']], B['area_km2'][B['closed']])),
        gain_m_yr_all=[float(x) for x in np.quantile(gain, (0.1, 0.5, 0.9))])


def load():
    """Everything the census computes, for this module's product and run.py's kinds."""
    d, b = lake_set()
    B = bodies(d, b)
    E = evaporation(d, b)
    GW = groundwater(d, b)
    springs_json = json.loads(SPRINGS.read_text())
    mature = springs_json['cases']['central']['states']['mature']['variants']['central']
    springs = dict(by_band=mature['by_band'], moon=mature['moon'],
                   air_uatm=springs_json['setting']['air_co2_atm'] * 1e6)
    cp = closed_parts(d, b, B, springs)
    return d, b, B, E, GW, cp, springs


def rnd(x, sig=4):
    import math
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        x = x.item()
    if isinstance(x, bool) or not isinstance(x, float):
        return x
    if not math.isfinite(x):
        return None
    return 0.0 if x == 0.0 else round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


EVIDENCE = (
    'Geometry and hydrology on committed products: geography\'s drainage of the climate run A28_dim5_moon over LOLA at '
    '16 pixels a degree above the GRAIL geoid (a first estimate; the lakes are regenerated once the design case '
    'settles), labelled as lake_oxygen labels it; groundwater stage 2\'s mature central water table re-solved from its '
    'own saved heads, reproducing its product; the sea chemistry\'s box model scaled to each closed lake\'s turnover; '
    'the springs\' modelled water. No lake model runs here, and no lunar lake has been measured.')
READING_RULE = (
    'Areas in km2 of lake surface, volumes in km3, depths in m below each lake\'s level. by_deepest bins the bodies by '
    'their deepest point, (lo, hi]; floor_within_km2 and volume_below_km3 are keyed by a floor depth, the area of bed '
    'within it and the water below it; floor_by_local_depth_km2 follows lake_anoxia\'s classes (floor_depth_edges_m). '
    'cumulative gives [bodies, area] with the deepest point at or within each depth, by band. closed.by_deepest bins '
    'the closed parts by their own deepest point, with years as area-weighted medians over the parts, keyed by '
    'weathering corner (high, central, low), and by the freshening share and margin for the chemocline. groundwater '
    'gains are m a '
    'year over each body\'s area; renewal_below gives, for each floor depth, the share of the water below it that the '
    'inflow across the floor below it brings in a year, and the years to make it denser by each margin (kg/m3).')


def main(argv=None) -> int:
    import time
    tic = time.time()
    d, b, B, E, GW, cp, springs = load()
    chem = json.loads(SEA_CHEMISTRY.read_text())
    out = dict(schema=SCHEMA,
               producer=dict(files={str(p.relative_to(ROOT)): sha256(p) for p in (
                   Path(__file__).resolve(), ROOT / 'research' / 'studies' / 'lake_oxygen' / 'lakes.py',
                   ROOT / 'geography' / 'groundwater_table.py', ROOT / 'geography' / 'drainage.py',
                   ROOT / 'research' / 'studies' / 'springs' / 'run.py')},
                             inputs=dict(inputs(), **{str(p.relative_to(ROOT)): sha256(p)
                                                      for p in (SEA_CHEMISTRY, SPRINGS)})),
               evidence=EVIDENCE, reading_rule=READING_RULE,
               setting=dict(deepest_edges_m=list(DEEPEST_EDGES_M), thresholds_m=list(THRESHOLDS_M),
                            floor_depth_edges_m=list(FLOOR_EDGES_M), size_edges_km2=[x if np.isfinite(x) else None
                                                                                    for x in SIZE_EDGES_KM2],
                            band_edges_deg=list(BAND_EDGES_DEG), salinity_density_kg_m3_per_g_kg=SALINITY_DENSITY,
                            margins_kg_m3=list(MARGINS_KG_M3), freshening=list(FRESHENING),
                            spring_temperatures_c=list(SPRING_TEMPERATURES_C)),
               census=summarise(B, d, b, E, GW, cp, springs, chem), seconds=None)
    out['seconds'] = round(time.time() - tic, 1)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rnd(out), indent=1, allow_nan=False) + '\n')
    print('wrote', OUT.relative_to(ROOT), f'in {time.time() - tic:.0f} s')
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
