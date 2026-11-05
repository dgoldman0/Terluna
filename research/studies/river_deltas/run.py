"""How the Open Moon's river deltas form and behave, on the young Moon and the mature one: a screen of every mouth.

    python3 -m research.studies.river_deltas.run        # results/river_deltas.json; the mouths' table on the drive

Every river mouth the terrain study found (research/runs/resources/terrain/terrain_16ppd.npz), into the seas, the
small waters at their level, the rain-fed lakes and the closed (soda) lakes, is screened for two states of the Moon
side by side, at lunar gravity:

- the river at its mouth: mean and bankfull discharge, and the bankfull channel by loose-bank scaling (Earth's
  dimensionless sand closure; Parker et al. 2007 for gravel) and as a threshold channel whose stress is held near
  the banks' own critical stress in pascals (Dunne & Jerolmack 2020);
- the load: the terrain study's rock at each mouth, split into bed material and wash load by the Rouse number at
  bankfull with lunar settling (Ferguson & Church 2004), and the fines' flocculation in each water;
- the plume through the lunar day: floating, level or plunging in each tenth of the day, the plunge depth from a
  densimetric criterion with lunar g', and the plume's internal Rossby radius at the Moon's rotation;
- what shapes the delta: river, wave and tidal sediment fluxes in Nienhuis et al.'s (2020) framework with gravity
  written in;
- the delta's scales: mouth bar, backwater, channel switching, bar mode, the front's failure height, the plain the
  monthly tide floods; and what share of the load the delta keeps.

The same functions run with Earth's gravity, rotation, seawater, semidiurnal tides and waves on four Earth deltas of
known type and on the lunar mouths. The flows are geography's drainage on the older climate run (A28_dim5_moon,
years 20-29), bracketed; the plume's temperatures come from the design years (A28R_301k_design_l15_snow_s934,
years 10-19). Latitude and the climate inputs are parameters, so the screen reruns when the turned Moon's design
run and the atlas round's shorelines arrive.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import resource
import time
from pathlib import Path

import numpy as np

from research.studies.river_deltas import physics as P
from shared.constants import (ANOMALISTIC_MONTH_DAYS, EARTH_SIDEREAL_DAY_S, JULIAN_DAY, JULIAN_YEAR_DAYS,
                              MOON_RADIUS, MOON_SURFACE_GRAVITY, SIDEREAL_MONTH_DAYS, STANDARD_GRAVITY)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'river_deltas.json'
DRIVE = Path('/media/projectspace/terluna-research/river-deltas')
TABLE = DRIVE / 'mouths.csv.gz'
SCHEMA = 'terluna.research.river-deltas/1'
GRIDS = ROOT / 'research' / 'runs' / 'resources' / 'terrain' / 'terrain_16ppd.npz'
DESIGN_YEARS = Path('/media/projectspace/terluna-research/session-tools/design_run/design_years.json')
INPUTS = dict(
    drainage=ROOT / 'geography' / 'results' / 'drainage.json',
    lakes=ROOT / 'geography' / 'products' / 'drainage_28pct_16ppd.npz',
    atlas=ROOT / 'geography' / 'products' / 'atlas_28pct_4ppd.npz',
    atlas_bodies=ROOT / 'geography' / 'results' / 'atlas.json',
    tides=ROOT / 'geography' / 'products' / 'tides_28pct_4ppd.npz',
    groundwater=ROOT / 'geography' / 'results' / 'groundwater_table.json',
    climatology=ROOT / 'climate' / 'gcm' / 'products' / 'climatology_A28_dim5_moon.npz',
    nearside_waves=ROOT / 'climate' / 'waves' / 'results' / 'nearside.json',
    smythii_waves=ROOT / 'climate' / 'waves' / 'results' / 'cycle.json',
    cold_waters=ROOT / 'research' / 'studies' / 'cold_waters' / 'results' / 'cold_waters.json',
    terrain=GRIDS, design_years=DESIGN_YEARS)
CODE = ('research/studies/river_deltas/run.py', 'research/studies/river_deltas/physics.py')

ROWS, COLS, PPD = 2880, 5760, 16
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
RHO = 1000.0
CHANNEL_M3S = 1.0                         # the drainage's channels: mean discharge of 1 m3/s and more
TERRAIN_GT_YR = dict(earliest=7.1, central=3.55, latest=0.71)   # the conservation study's land rock, Gt a year
# The young land yields more rock than the mature cover (design bracket): the storm-water study's bare dryland soils
# detach 10-23 times what crusted ones do (3.5-5.0 against 0.22-0.24 Gt a year, upper bounds) and 2.6-12.6 times
# the terrain study's central lowering; Earth's plowed fields erode 1-2 orders of magnitude faster than ground
# under native vegetation (Montgomery 2007); new planting cuts it. Paired with the terrain cases latest to earliest.
YOUNG_YIELD = dict(earliest=12.0, central=6.0, latest=3.0)
FLOW_CASES = ('low', 'central', 'high')
# Bankfull over mean discharge (design guesses, the terrain study's dam builders' floods): steady rivers 4, flashy 10;
# the young Moon's humid rivers lose baseflow to the filling crust, 6.
BANKFULL_RATIO = dict(young=dict(humid=6.0, dry=10.0), mature=dict(humid=4.0, dry=10.0))
INTERMITTENCY = 0.1                       # share of the year at bankfull that carries the bed-material load (design)
CF_SAND = 0.004                           # bankfull friction coefficient of sand beds with dunes (design, 0.002-0.006)
ROUSE_SPLIT = 0.11                        # wash load below it: Earth's 62.5 um at u* = 0.07 m/s (calibration)
GRAVEL_D = 0.002
# Grain sizes. The young Moon's regolith: median 72 um, poorly sorted (Carrier 2003: medians 40-130 um, Cu 9-30,
# 90% of the mass between 10 um and 1 mm), here lognormal with sigma 1.5 in ln D; mineral density 3,100 kg/m3.
# The mature candidate soils are standard textures from sand to clay built from clay, silt and sand modes.
MODES = dict(clay=(0.8e-6, 0.7), silt=(12e-6, 0.8), sand=(150e-6, 0.6))
SOILS = dict(regolith=dict(lognormal=(72e-6, 1.5), rho_s=3100.0),
             sand=dict(clay=0.05, silt=0.10, sand=0.85, rho_s=2700.0),
             sandy_loam=dict(clay=0.10, silt=0.25, sand=0.65, rho_s=2700.0),
             loam=dict(clay=0.20, silt=0.40, sand=0.40, rho_s=2700.0),
             clay_loam=dict(clay=0.33, silt=0.34, sand=0.33, rho_s=2700.0),
             clay=dict(clay=0.55, silt=0.30, sand=0.15, rho_s=2700.0))
EDGES = np.geomspace(1e-7, 1.6e-2, 61)
SIZES = np.sqrt(EDGES[1:] * EDGES[:-1])
# Candidate soils by climate zone, from the soils note's WRB groups read as textures (a proposal, not a decision).
ZONE_SOILS = {'hyper-arid': ('sand', 'sandy_loam'), 'arid': ('sand', 'sandy_loam', 'loam'),
              'semi-arid': ('sandy_loam', 'loam', 'clay'), 'dry sub-humid': ('loam', 'clay_loam', 'clay'),
              'humid': ('loam', 'clay_loam', 'clay')}
CLASSES = ('hyper-arid', 'arid', 'semi-arid', 'dry sub-humid', 'humid')
BANDS = ('0-15', '15-30', '30-45', '45-60', '60-91')
# Banks. Critical stress of bare soil, WEPP's 2.67 + 0.065 clay% - 0.058 VFS% Pa (Alberts et al. 1995, as the
# storm-water study uses it); roots raise it (design: young bare to newly planted 1-1.5, mature 1.5-6 with 3).
ROOT_FACTOR = dict(young=(1.0, 1.0, 1.5), mature=(1.5, 3.0, 6.0))
# Bank cohesion for the banks' standing height (design bracket): lunar soil in place 0.1-1 kPa (Carrier et al. in the
# Lunar Sourcebook) to 3 kPa newly planted; rooted developed soils 5-20 kPa.
BANK_COHESION_PA = dict(young=(300.0, 1000.0, 3000.0), mature=(5000.0, 10000.0, 20000.0))
# Waters. Seas (sea chemistry): young, at 1,000 years, 0.02-0.06 g/kg with calcium and magnesium 0.07-0.21 mmol/kg;
# mature, held fresh at 0.142 g/kg with 0.84 mmol/kg of calcium. Rain-fed lakes 0.02 and 0.04 g/kg. Soda lakes (the
# lake census's closed lakes): young 0.1-35 g/kg, mature 1-300. Rivers: alkalinity 1.26 meq/L at the catchment rate
# (terrain study), about 0.11 g/kg with 0.63 mmol/L of calcium and magnesium, the young soils 3.2 times that.
WATERS = dict(
    young=dict(sea_salt=0.04, sea_div=0.14, lake_salt=0.02, lake_div=0.10, soda_salt=(0.1, 1.0, 35.0),
               river_salt=0.34, river_div=2.0, organic=(0.0, 0.0, 0.2)),
    mature=dict(sea_salt=0.142, sea_div=0.84, lake_salt=0.04, lake_div=0.25, soda_salt=(1.0, 35.0, 300.0),
                river_salt=0.11, river_div=0.63, organic=(0.3, 0.6, 0.9)))
HALINE = 0.8                              # kg/m3 per g/kg, carbonate waters (the lake census)
HALINE_SEAWATER = 0.75                    # kg/m3 per g/kg, Earth seawater at 25 C (1023.3 against 997.0 at 35)
CCC_DIVALENT_MMOL = (0.3, 3.0)            # calcium coagulating smectite and illite clays (search summaries; open)
CCC_SODA_NA_MMOL = (30.0, 300.0)          # sodium at high pH, edges negative (search summary: about 100; open)
FLOC_BELOW_M = 16e-6                      # clay and fine silt form flocs
FLOC_SETTLING_EARTH = (1e-4, 3e-4, 1e-3)  # m/s, Earth's flocs; Stokes-like, so times g / g_Earth (design bracket)
# Temperatures (design estimates): land air by aridity class and tenth of the lunar day (design years), carried to
# the mouth's level along cold_waters' lapse rate from the land around it, less the catchment's height above the
# mouth (a share of the local relief); the seas' surface at the resort's median 29.6 C with a 1 K swing, their deep
# water at the coldest surface, 25-29 C; lakes the land air at their level plus an offset.
SEA_SURFACE_C, SEA_SWING_K, SEA_DEEP_C = (25.0, 29.6, 30.0), 0.5, (25.0, 26.0, 29.0)
LAKE_OFFSET_K = (0.0, 1.2, 2.4)
CATCHMENT_SHARE = (0.0, 0.5, 1.0)
LAND_BOX_DEG = 5.6                        # a T21 cell
WIND_10M = dict(moon=2.6, earth=7.0)      # m/s: the seas' median (wave study); Earth's oceans' mean (design)
DIRECTION_RESULTANT = 0.5                 # waves' directional steadiness where no wave product exists (design)
STRAIT_ZONES = {'western channels': dict(lon=(35.0, 44.5), lat=(-11.5, -0.5)),
                'Crisium passage': dict(lon=(53.0, 62.0), lat=(6.0, 11.5))}
# The channel bracket from the Earth control (screened large sand rivers 0.8-3.6 times too wide, here 0.5-2, and
# 3-10 times too steep), with the routing grid's slope, the raw shore terrain, as the steep side for young fans.
GEOMETRY_BRACKET = {'slope_third': dict(slope=1 / 3), 'slope_tenth': dict(slope=0.1),
                    'half_width_slope_third': dict(width=0.5, slope=1 / 3), 'half_width_slope_tenth': dict(width=0.5, slope=0.1),
                    'double_width_slope_third': dict(width=2.0, slope=1 / 3),
                    'double_width_slope_tenth': dict(width=2.0, slope=0.1), 'grid_slope': dict(slope='grid')}
STRAIT_SPEED = 2.25                       # m/s, the western channels' races (climate/ocean, 2.2-2.3)
MUD = dict(cohesion_pa=(500.0, 2000.0, 5000.0), rho_sat=1600.0, slopes_deg=(1.0, 5.0))   # design bracket
SAND_REPOSE_DEG = 32.0
BELT_WIDTHS = 10.0                        # channel belt in channel widths (design)
DRY_DENSITY = dict(sand=1600.0, mud=1000.0)
SEA_LEVEL_RATES_M_YR = dict(uptake_unbalanced=0.17, sag_central=0.015, sag_warm=0.41)
EARTH = dict(g=STANDARD_GRAVITY, omega=2.0 * math.pi / EARTH_SIDEREAL_DAY_S, tide_period=12.4206012 * 3600.0,
             sea_salt=35.0, haline_sea=HALINE_SEAWATER, wind='earth')
MOON = dict(g=MOON_SURFACE_GRAVITY, omega=2.0 * math.pi / (SIDEREAL_MONTH_DAYS * JULIAN_DAY),
            tide_period=ANOMALISTIC_MONTH_DAYS * JULIAN_DAY, sea_salt=None, haline_sea=HALINE, wind='moon')
EVIDENCE = ('A screen: closed-form relations evaluated per mouth on 1-D arrays, no morphodynamic model. Mouths, '
            'rock and delta areas are the terrain study\'s routing on geography\'s drainage (16 px/deg, older '
            'climate run A28_dim5_moon); flows are bracketed for the design climate; temperatures are design '
            'estimates from the design years\' land air by aridity class; soils are candidates, not choices; waves '
            'come from the wave study where it computed them and from fetch-limited growth elsewhere.')
READING_RULE = ('States: young (regolith, bare banks, own-catchment rock times the young yield factor, drier rivers, '
                'soft young seas) and mature (each candidate soil, rooted banks, the terrain study\'s cases, held '
                'seas). Shares of mouths are by count, of discharge by mean flow, of load by rock. Plume shares are of '
                'the ten tenths of the lunar day, at mean flow and at bankfull. Dominance follows Nienhuis et al. '
                '(2020): the largest of the river\'s, the waves\' and the tide\'s sediment flux. Lengths in m unless '
                'named in km, times in years, discharges m3/s, loads kg/s or Mt a year as named.')


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def rnd(x, sig=4):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, (bool, np.bool_)) or x is None or isinstance(x, str):
        return bool(x) if isinstance(x, np.bool_) else x
    if isinstance(x, (int, np.integer)):
        return int(x)
    x = float(x)
    if not math.isfinite(x):
        return None
    return float(f'{x:.{sig}g}') if x else 0.0


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    s = np.sin(0.5 * (p2 - p1)) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(0.5 * np.radians(lon2 - lon1)) ** 2
    return 2.0 * MOON_RADIUS / 1e3 * np.arcsin(np.sqrt(np.clip(s, 0.0, 1.0)))


def nearest(lat, lon, nlat, nlon):
    """Index of the nearest node and its distance in km for each point (points in modest numbers)."""
    idx = np.empty(lat.size, int)
    dist = np.empty(lat.size)
    for i0 in range(0, lat.size, 256):
        d = haversine_km(lat[i0:i0 + 256, None], lon[i0:i0 + 256, None], nlat[None, :], nlon[None, :])
        idx[i0:i0 + 256] = np.argmin(d, axis=1)
        dist[i0:i0 + 256] = d[np.arange(d.shape[0]), idx[i0:i0 + 256]]
    return idx, dist


# --- Grains -----------------------------------------------------------------------------------------------------------

def grain_fractions(name):
    spec = SOILS[name]
    if 'lognormal' in spec:
        f = P.lognormal_fractions(EDGES, *spec['lognormal'])
    else:
        f = sum(spec[k] * P.lognormal_fractions(EDGES, *MODES[k]) for k in ('clay', 'silt', 'sand'))
    return f / f.sum()


def bank_critical_pa(fractions, root):
    """WEPP's bare-soil critical shear from the texture's clay (< 2 um) and very fine sand (50-100 um), times roots."""
    clay = 100.0 * fractions[SIZES < 2e-6].sum()
    vfs = 100.0 * fractions[(SIZES >= 50e-6) & (SIZES < 100e-6)].sum()
    return max(2.67 + 0.065 * clay - 0.058 * vfs, 0.1) * root


def flocculated_share(divalent_mmol, soda_na_mmol, organic):
    """Share of the fines in flocs: salt coagulation between the critical concentrations' bracket (log-linear),
    calcium and magnesium in fresh water or sodium in the soda lakes, combined with organic flocculation."""
    lo, hi = CCC_DIVALENT_MMOL
    salt = np.clip(np.log(np.maximum(divalent_mmol, 1e-9) / lo) / math.log(hi / lo), 0.0, 1.0)
    lo, hi = CCC_SODA_NA_MMOL
    soda = np.clip(np.log(np.maximum(soda_na_mmol, 1e-9) / lo) / math.log(hi / lo), 0.0, 1.0)
    salt = np.maximum(salt, soda)
    return 1.0 - (1.0 - salt) * (1.0 - organic)


# --- Inputs -----------------------------------------------------------------------------------------------------------

def design_climate():
    d = json.loads(DESIGN_YEARS.read_text())['new']
    ar = d['aridity']
    air = np.array([[t['air_c'] for t in ar[c]['by_tenth']] for c in CLASSES])          # (class, tenth)
    hour = np.array([t['hour_deg'] for t in ar[CLASSES[0]]['by_tenth']])
    runoff = sum(ar[c]['land_share'] * np.mean([t['rain_mm_day'] - t['evaporation_mm_day'] for t in ar[c]['by_tenth']])
                 for c in CLASSES)
    weight = np.array([[ar[c]['land_share'] * ar[c]['latitude_shares'][b] for b in BANDS] for c in CLASSES])
    return dict(run=d['run'], air=air, hour=hour, runoff_mm_day=float(runoff), band_class=np.argmax(weight, axis=0),
                band_weights=weight / weight.sum(axis=0))


def old_runoff_mm_day():
    c = np.load(INPUTS['climatology'])
    land = c['lsm'] > 0.5
    lat = c['lat']
    w = np.cos(np.radians(lat))[:, None] * np.ones(land.shape)
    pr, ev = c['pr_mm_day'], c['evap_mm_day']
    net = pr - np.abs(ev)
    return float((net * w)[land].sum() / w[land].sum())


def band_index(lat):
    a = np.abs(lat)
    return np.minimum(np.searchsorted(np.array([15.0, 30.0, 45.0, 60.0]), a, side='right'), 4)


def load_mouths():
    """Every mouth with a resolved river: the terminal channel cells of the drainage (1 m3/s and more) next to the
    standing water cell where the land's rock first meets that water, with what that water is."""
    z = np.load(GRIDS)
    mouths = z['delta_mouth_cell'].astype(np.int64)
    chan = z['channel_cell'].astype(np.int64)
    q = z['channel_discharge_m3s'].astype(float)
    slope = z['channel_slope'].astype(float)
    rows, cols = chan // COLS, chan % COLS
    best = np.full(chan.size, -1)
    nbq = np.zeros(chan.size)
    for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)):
        rr = rows + dr
        ok = (rr >= 0) & (rr < ROWS)
        nb = np.where(ok, rr * COLS + (cols + dc) % COLS, 0)
        k = np.minimum(np.searchsorted(mouths, nb), mouths.size - 1)
        is_m = ok & (mouths[k] == nb)
        best = np.where((best < 0) & is_m, k, best)
        j = np.minimum(np.searchsorted(chan, nb), chan.size - 1)
        nbq = np.maximum(nbq, np.where(ok & (chan[j] == nb), q[j], 0.0))
    terminal = np.flatnonzero((best >= 0) & (nbq <= q * (1.0 + 1e-6)))
    k = best[terminal]
    q_m = np.bincount(k, weights=q[terminal], minlength=mouths.size)
    qmax = np.zeros(mouths.size)
    np.maximum.at(qmax, k, q[terminal])
    main = terminal[q[terminal] >= qmax[k]]                  # each mouth's largest river
    land_cell = np.full(mouths.size, -1)
    s_m = np.zeros(mouths.size)
    land_cell[best[main]] = chan[main]
    s_m[best[main]] = slope[main]
    river = np.flatnonzero(q_m >= CHANNEL_M3S)
    cell = mouths[river]
    r_, c_ = cell // COLS, cell % COLS
    lat = 90.0 - (r_ + 0.5) / PPD
    lon = (c_ + 0.5) / PPD
    lon = np.where(lon > 180.0, lon - 360.0, lon)
    body = z['delta_body_unit'][river].astype(np.int64)
    kind_u = z['unit_kind']
    n_units = kind_u.size
    lake = body < n_units
    kind = np.where(lake, np.where(kind_u[np.minimum(body, n_units - 1)] == 2, 'soda', 'lake'), 'common')
    lake_m2 = np.where(lake, z['unit_lake_m2'][np.minimum(body, n_units - 1)], np.nan)
    lake_m3 = np.where(lake, z['unit_lake_m3'][np.minimum(body, n_units - 1)], np.nan)
    relief = z['cell_local_relief_m'].ravel()[np.maximum(land_cell[river], 0)].astype(float)
    m = dict(mouth_cell=cell, land_cell=land_cell[river], lat=lat, lon=lon, q=q_m[river], slope_grid=s_m[river],
             rock_per_gt=z['delta_rock_kg_yr_per_gt_yr'][river].astype(float),
             depth0=z['delta_depth_at_mouth_m'][river].astype(float), kind=kind, body=body,
             lake_m2=lake_m2, lake_m3=lake_m3, relief=relief,
             unit_inflow=np.where(lake, z['unit_inflow_m3_yr'][np.minimum(body, n_units - 1)], np.nan))
    for case in TERRAIN_GT_YR:
        for y in (1000, 10000):
            m[f'delta_{case}_{y}'] = z[f'delta_{case}_{y}_km2'][river].astype(float)
    checks = dict(mouths_all=int(mouths.size), river_mouths=int(river.size),
                  terminal_channel_cells=int(terminal.size), channel_cells=int(chan.size),
                  rock_at_river_mouths_share=float(m['rock_per_gt'].sum() / z['delta_rock_kg_yr_per_gt_yr'].sum()))
    return m, checks


def attach_waters(m, dc):
    """Receiving water by kind and sea, water level, depth, tide, waves, the land around for the air's level."""
    atlas = np.load(INPUTS['atlas'])
    alat, alon = atlas['lat_deg'], atlas['lon_deg']
    h4 = atlas['height_m'].astype(float)
    wl = atlas['water_label']
    drainage = json.loads(INPUTS['drainage'].read_text())
    sea_level = float(drainage['sea_level_m'])
    tides = np.load(INPUTS['tides'])
    label_grid = np.full(h4.shape, -1, np.int64)
    label_grid[tides['rows'], tides['cols']] = tides['labels']
    node_grid = np.full(h4.shape, -1, np.int64)
    node_grid[tides['rows'], tides['cols']] = np.arange(tides['rows'].size)
    rng = tides['monthly_range_median_m'].astype(float)
    lon360 = np.where(m['lon'] < 0.0, m['lon'] + 360.0, m['lon'])
    alon360 = np.where(alon < 0.0, alon + 360.0, alon)
    order = np.argsort(alon360)
    r4 = np.argmin(np.abs(alat[None, :] - m['lat'][:, None]), axis=1)
    c4 = order[np.minimum(np.searchsorted(alon360[order], lon360), alon.size - 1)]
    common = m['kind'] == 'common'
    sea_label = np.full(m['lat'].size, -1)
    node = np.full(m['lat'].size, -1)
    for rad in (0, 1, 2, 3):
        todo = common & (sea_label < 0)
        if not todo.any():
            break
        for dr in range(-rad, rad + 1):
            for dcol in range(-rad, rad + 1):
                rr = np.clip(r4 + dr, 0, h4.shape[0] - 1)
                cc = (c4 + dcol) % h4.shape[1]
                hit = todo & (sea_label < 0) & (label_grid[rr, cc] >= 0)
                sea_label[hit] = label_grid[rr, cc][hit]
                node[hit] = node_grid[rr, cc][hit]
    # A sea is one connected water at 16 px/deg, so every mouth of a connected common-level water whose mouths
    # mostly meet one sea's tide nodes belongs to that sea; its tide is the nearest node of that sea.
    for comp in np.unique(m['body'][common]):
        members = common & (m['body'] == comp)
        found = sea_label[members & (sea_label >= 0)]
        if found.size:
            sea_label[members] = np.bincount(found).argmax()
    tlat, tlon = tides['lat_deg'][tides['rows']], tides['lon_deg'][tides['cols']]
    tlon = np.where(tlon > 180.0, tlon - 360.0, tlon)
    for lab in np.unique(sea_label[(sea_label >= 0) & (node < 0)]):
        sel = np.flatnonzero((sea_label == lab) & (node < 0))
        cand = np.flatnonzero(tides['labels'] == lab)
        i, _ = nearest(m['lat'][sel], m['lon'][sel], tlat[cand], tlon[cand])
        node[sel] = cand[i]
    kind = np.where(common & (sea_label >= 0), 'sea', np.where(common, 'small', m['kind']))
    m['kind'] = kind
    m['sea_label'] = sea_label
    m['tide_range'] = np.where(node >= 0, rng[np.maximum(node, 0)], 0.0)
    # Areas of the common-level waters from the atlas's labels (the nearest labelled water within 0.75 degrees, else
    # 10 km2, below the 4 px/deg atlas); lakes from the terrain units.
    cell_km2 = (MOON_RADIUS / 1e3 * math.radians(0.25)) ** 2 * np.cos(np.radians(alat))
    lab_area = np.bincount(wl.ravel().clip(0), weights=np.repeat(cell_km2, wl.shape[1]), minlength=wl.max() + 1)
    here = wl[r4, c4].astype(np.int64)
    for rad in (1, 2, 3):
        for dr in range(-rad, rad + 1):
            for dcol in range(-rad, rad + 1):
                cand = wl[np.clip(r4 + dr, 0, wl.shape[0] - 1), (c4 + dcol) % wl.shape[1]]
                here = np.where((here <= 0) & (cand > 0), cand, here)
    m['water_km2'] = np.where(np.isfinite(m['lake_m2']), m['lake_m2'] / 1e6,
                              np.where(here > 0, lab_area[here.clip(0)], 10.0))
    sea_depth = np.maximum(sea_level - h4[r4, c4], 0.0)
    lake_depth = np.where(np.isfinite(m['lake_m3']), m['lake_m3'] / np.maximum(m['lake_m2'], 1.0), 0.0)
    m['depth_avail'] = np.maximum(m['depth0'], np.where(np.isfinite(m['lake_m2']), lake_depth, sea_depth))
    # Water level at the mouth: lakes from the drainage product's lake cells, the common level elsewhere.
    lk = np.load(INPUTS['lakes'])
    lc, ll = lk['lake_cell'].astype(np.int64), lk['lake_level_m'].astype(float)
    o = np.argsort(lc)
    j = np.minimum(np.searchsorted(lc[o], m['mouth_cell']), lc.size - 1)
    in_lake = lc[o][j] == m['mouth_cell']
    m['level'] = np.where(in_lake & np.isfinite(m['lake_m2']), ll[o][j], sea_level)
    # Mean height of the land around each mouth over a T21 cell, where the air the design years give sits.
    from scipy.ndimage import uniform_filter
    land = (h4 > sea_level) & (wl <= 0)
    size = int(round(LAND_BOX_DEG * 4))
    num = uniform_filter(np.where(land, h4, 0.0), size=size, mode=('nearest', 'wrap'))
    den = uniform_filter(land.astype(float), size=size, mode=('nearest', 'wrap'))
    m['land_box_m'] = np.where(den[r4, c4] > 0.01, num[r4, c4] / np.maximum(den[r4, c4], 1e-9), m['level'])
    m['band'] = band_index(m['lat'])
    m['zone'] = dc['band_class'][m['band']]
    m['humid'] = np.array([CLASSES[z] == 'humid' for z in m['zone']])
    strait = np.zeros(m['lat'].size, bool)
    for zone in STRAIT_ZONES.values():
        strait |= ((kind == 'sea') & (m['lon'] >= zone['lon'][0]) & (m['lon'] <= zone['lon'][1])
                   & (m['lat'] >= zone['lat'][0]) & (m['lat'] <= zone['lat'][1]))
    m['strait'] = strait
    return dict(sea_level=sea_level, names=sea_names(), tide_nodes=int(rng.size))


def sea_names():
    bodies = json.loads(INPUTS['atlas_bodies'].read_text())['bodies']
    named = {433: 'Near-side Sea', 1417: 'South Pole-Aitken Sea', 874: 'Smythii-Marginis'}
    return {int(b['label']): named.get(int(b['label']), (b.get('water_names') or [str(b['label'])])[0])
            for b in bodies}


def attach_waves(m):
    """Waves at each mouth: the wave study's nearside coastal power and Smythii-Marginis' nodes; elsewhere
    fetch-limited growth at lunar gravity with the seas' median wind over the receiving water's size."""
    g = MOON['g']
    n = m['lat'].size
    h0, t0, rd = np.zeros(n), np.zeros(n), np.full(n, DIRECTION_RESULTANT)
    source = np.full(n, 'fetch', dtype=object)
    hf, tf = P.fetch_waves(WIND_10M['moon'], np.sqrt(m['water_km2']) * 1e3, g)
    h0[:], t0[:] = hf, tf / 1.2
    near = json.loads(INPUTS['nearside_waves'].read_text())['maps']
    cn = np.array(near['coastal_nodes'])
    nlat, nlon = np.array(near['lat'])[cn], np.array(near['lon'])[cn]
    sel = np.flatnonzero(m['sea_label'] == 433)
    if sel.size:
        i, d = nearest(m['lat'][sel], m['lon'][sel], nlat, nlon)
        power = np.array(near['coastal_mean_power_W_m'])[i]
        tm = np.array(near['energy_mean_tm01_s'])[cn][i]
        h0[sel] = P.wave_height_from_power(power, tm, g)
        t0[sel] = tm
        rd[sel] = np.array(near['direction_resultant'])[cn][i]
        source[sel] = 'nearside SWAN coast'
    cyc = json.loads(INPUTS['smythii_waves'].read_text())
    sel = np.flatnonzero(m['sea_label'] == 874)
    if sel.size:
        glat, glon = np.array(cyc['grid']['node_latitude_deg']), np.array(cyc['grid']['node_longitude_deg'])
        glon = np.where(glon > 180.0, glon - 360.0, glon)
        i, _ = nearest(m['lat'][sel], m['lon'][sel], glat, glon)
        cols = cyc['series']['offshore_columns']
        tcol = [k for k, c in enumerate(cols) if 'tm01' in c.lower()]
        tm = float(np.median(np.array(cyc['series']['offshore'])[:, tcol[0]])) if tcol else 10.0
        h0[sel] = np.array(cyc['statistics']['node_rms_hs_m'])[i]
        t0[sel] = tm
        rd[sel] = np.array(cyc['statistics']['node_direction_resultant'])[i]
        source[sel] = 'Smythii-Marginis SWAN nodes'
    m['wave_h'], m['wave_t'], m['wave_dir'], m['wave_source'] = h0, t0, rd, source


def temperatures(m, dc, lapse_k_km, catch_share, lake_offset, sea_surface, sea_deep):
    """River, surface and deep water temperatures (C) at each mouth through the ten tenths of the lunar day."""
    air = dc['air'][m['zone']] + lapse_k_km * (m['land_box_m'] - m['level'])[:, None] / 1e3
    river = air - lapse_k_km * catch_share * m['relief'][:, None] / 1e3
    swing = SEA_SWING_K * np.cos(np.radians(dc['hour']))[None, :]
    sea = m['kind'] == 'sea'
    surface = np.where(sea[:, None], sea_surface + swing, air + lake_offset + swing)
    deep = np.where(sea, sea_deep, (air + lake_offset).min(axis=1) - SEA_SWING_K)
    return river, surface, deep


# --- The screen -------------------------------------------------------------------------------------------------------

def screen(m, state, soil, env, rock_kg_yr, q_mean, root=None, organic=None, waters=None, temps=None,
           tau_bank=None, fractions=None, rho_s=None, geometry=None):
    """One state and soil over every mouth: channel, load, plume, shaping balance, scales and keep share."""
    g = env['g']
    w = dict(WATERS[state]) if waters is None else waters
    fractions = grain_fractions(soil) if fractions is None else fractions
    rho_s = SOILS[soil]['rho_s'] if rho_s is None else rho_s
    r = rho_s / RHO - 1.0
    root = ROOT_FACTOR[state][1] if root is None else root
    organic = w['organic'][1] if organic is None else organic
    n = q_mean.size
    ratio = np.where(m['humid'], BANKFULL_RATIO[state]['humid'], BANKFULL_RATIO[state]['dry']) \
        if 'bankfull_ratio' not in m else m['bankfull_ratio']
    q_bf = q_mean * ratio
    t_river, t_surf, t_deep = temps
    nu = P.kinematic_viscosity(np.mean(t_river, axis=1))
    ws = P.settling_velocity(SIZES[None, :], g, rho_s, nu[:, None])                  # (mouth, size)
    bank = (np.array([bank_critical_pa(fractions, root)]) * np.ones(n)) if tau_bank is None else tau_bank
    cohesion = BANK_COHESION_PA[state][1] if 'bank_cohesion' not in m else m['bank_cohesion']
    # Bed material and the channel together: the Rouse number at bankfull sets the split, the bed's median sets
    # the bed's threshold and the loose closure; a few passes settle them.
    bed = np.ones((n, SIZES.size), bool)
    for _ in range(6):
        fb = (fractions[None, :] * bed).sum(axis=1)
        cum = np.cumsum(fractions[None, :] * bed, axis=1)
        k50 = np.argmax(cum >= 0.5 * np.maximum(fb, 1e-12)[:, None], axis=1)
        d50 = np.where(fb > 1e-6, SIZES[k50], SIZES[np.flatnonzero(fractions > 1e-6)[-1]])
        loose = P.loose_channel(q_bf, d50, g, r, nu, CF_SAND, gravel_d=GRAVEL_D)
        tau_bed = P.shields_critical(P.particle_reynolds(d50, g, r, nu)) * RHO * r * g * d50
        ch = P.threshold_channel(loose, tau_bed, bank, P.bank_height_limit(cohesion, g))
        u_star = np.sqrt(ch['stress'] / RHO)
        bed = P.rouse(ws, u_star[:, None]) >= ROUSE_SPLIT
    f_bed = (fractions[None, :] * bed).sum(axis=1)
    # The bed's capacity at bankfull against the bed material it is given: above one the channel can cut down and
    # pass its sand on to the mouth, below one it builds bars and spreads.
    supply = np.maximum(rock_kg_yr * f_bed, 1.0) / (INTERMITTENCY * YEAR_S * rho_s)
    capacity = ch['width'] * P.unit_load(ch['stress'] / (RHO * r * g * d50), CF_SAND, d50, g, r) / supply
    ch['narrow'] = ch['bank_limited']
    # The channel the classifications read. The Earth control found the screened channels of large sand rivers
    # about twice too wide and 3-10 times too steep, so a geometry bracket multiplies the width and the slope (or
    # takes the routing grid's slope) at the screened velocity, depth following by continuity; the bed's stress,
    # its split and the threshold physics stay as screened.
    geo = dict(width=ch['width'], depth=ch['depth'], slope=ch['slope'], velocity=ch['velocity'])
    if geometry is not None:
        wf = geometry.get('width', 1.0)
        geo = dict(width=ch['width'] * wf, depth=ch['depth'] / wf, velocity=ch['velocity'],
                   slope=m['slope_grid'] if geometry.get('slope') == 'grid' else ch['slope'] * geometry.get('slope', 1.0))
    h_mean = geo['depth'] * (q_mean / q_bf) ** (2.0 / 3.0)
    u_mean = q_mean / (geo['width'] * h_mean)
    # The load in the water: flux-weighted concentration, at mean flow and bankfull through a rating of exponent 1.5.
    c_flux = rock_kg_yr / (q_mean * YEAR_S)
    susp = 1.0 - 0.5 * f_bed                                   # the wash load and half the bed material in suspension
    out = dict(q_mean=q_mean, q_bf=q_bf, width=geo['width'], depth=geo['depth'], velocity=geo['velocity'],
               slope=geo['slope'], stress=ch['stress'], bank_pa=bank, narrow=ch['narrow'],
               loose_width=loose['width'], loose_depth=loose['depth'], loose_velocity=loose['velocity'],
               loose_slope=loose['slope'], loose_stress=loose['stress'], d50_bed=d50, bed_share=f_bed,
               capacity_ratio=capacity, conc_kg_m3=c_flux, rock_mt_yr=rock_kg_yr / 1e9)
    # Fines in flocs, in the river and in the receiving water.
    if 'soda_na' in m:
        soda_na = m['soda_na']
    else:
        soda_na = np.where(m['kind'] == 'soda', w['soda_salt'][1] * 1000.0 / 84.0, 0.0)
    recv_div = np.where(m['kind'] == 'sea', w['sea_div'], np.where(m['kind'] == 'small', w['sea_div'], w['lake_div']))
    if 'recv_div' in m:
        recv_div = m['recv_div']
    floc = flocculated_share(np.maximum(recv_div, w['river_div']), soda_na, organic)
    w_floc = FLOC_SETTLING_EARTH[1] * g / STANDARD_GRAVITY
    fine = SIZES < FLOC_BELOW_M
    ws_eff = np.where(fine[None, :], (1.0 - floc[:, None]) * ws + floc[:, None] * np.maximum(w_floc, ws), ws)
    # Plume through the lunar day.
    if 'recv_salt' in m:
        s_recv = m['recv_salt']
    else:
        s_recv = np.where(m['kind'] == 'soda', w['soda_salt'][1], np.where(np.isin(m['kind'], ('sea', 'small')),
                          w['sea_salt'] if env['sea_salt'] is None else env['sea_salt'], w['lake_salt']))
    haline = np.where(np.isin(m['kind'], ('sea', 'small')), env['haline_sea'], HALINE)
    shares = {}
    for name, qq, cc in (('mean', q_mean, c_flux * ratio ** -0.25), ('bankfull', q_bf, c_flux * ratio ** 0.25)):
        rho_r = P.density(t_river, w['river_salt'], (cc * susp)[:, None], HALINE, rho_s)
        rho_s_ = P.density(t_surf, s_recv[:, None], 0.0, haline[:, None])
        drho = rho_r - rho_s_
        gr = g * np.abs(drho) / rho_s_
        hp = P.plunge_depth((qq / geo['width'])[:, None], gr)
        plunging = (drho > 1e-3) & (hp < m['depth_avail'][:, None])
        floating = drho < -1e-3
        level = ~plunging & ~floating
        rho_deep = P.density(t_deep, s_recv, 0.0, haline)
        floor = plunging & (rho_r > rho_deep[:, None])
        shares[name] = dict(floating=floating.mean(axis=1), level=level.mean(axis=1), plunging=plunging.mean(axis=1),
                            to_floor=floor.mean(axis=1),
                            plunge_depth=np.where(plunging.any(axis=1), np.nanmean(np.where(plunging, hp, np.nan),
                                                                                   axis=1), np.nan),
                            drho_mean=drho.mean(axis=1), g_reduced=gr.mean(axis=1))
    for name, s in shares.items():
        for k, v in s.items():
            out[f'plume_{name}_{k}'] = v
    # Rossby radius of the plume against its reach: the floating layer as deep as the river, or the underflow as
    # thick as its plunge depth; the reach is where half the wash load's median grain or floc has settled.
    sm = shares['mean']
    thick = np.where(sm['plunging'] > 0.5, np.nan_to_num(sm['plunge_depth'], nan=h_mean), h_mean)
    rossby = P.rossby_radius(sm['g_reduced'], thick, m['lat'], env['omega'])
    wash = ~bed
    wcum = np.cumsum(fractions[None, :] * wash, axis=1)
    fw = wcum[:, -1]
    kw = np.argmax(wcum >= 0.5 * np.maximum(fw, 1e-12)[:, None], axis=1)
    ws_wash = ws_eff[np.arange(n), kw]
    reach = math.log(2.0) * 0.5 * u_mean * h_mean / np.maximum(ws_wash, 1e-9)
    out.update(rossby_km=rossby / 1e3, plume_reach_km=reach / 1e3, floc_share=floc,
               coastal_current=(reach > rossby) & (sm['floating'] > 0.5))
    # What the delta keeps: its bed material and what of the wash load settles within its own reach, unless the
    # plume plunges and carries it down the front as an underflow.
    r_delta = np.maximum(np.sqrt(2.0 * m['delta_km2'] * 1e6 / math.pi), 500.0)
    t_reach = r_delta / np.maximum(0.5 * u_mean, 1e-3)
    settled = 1.0 - np.exp(-ws_eff * t_reach[:, None] / np.maximum(h_mean, 0.1)[:, None])
    wash_kept = (fractions[None, :] * wash * settled).sum(axis=1)
    keep = f_bed + wash_kept * (1.0 - sm['plunging'])
    kept_kg_yr = keep * rock_kg_yr
    out.update(keep_share=keep, keep_with_underflow=f_bed + wash_kept, delta_radius_km=r_delta / 1e3)
    # What shapes the delta.
    q_river = rock_kg_yr / YEAR_S
    h0, t0 = m['wave_h'], m['wave_t']
    if env['wind'] == 'earth' and 'earth_waves' not in m:
        hf, tf = P.fetch_waves(WIND_10M['earth'], np.sqrt(m['water_km2']) * 1e3, g)
        h0, t0 = hf, tf / 1.2
    q_wave = 2.0 * P.longshore_max(h0, t0, g, rho_s) * m['wave_dir']
    qw_tide = P.tidal_discharge(0.5 * m['tide_range'], geo['width'], geo['slope'], env['tide_period'])
    q_tide = qw_tide / q_mean * q_river
    gentle = P.tidal_discharge(0.5 * m['tide_range'], geo['width'], 0.1 * geo['slope'], env['tide_period'])
    flux = np.stack([q_river, q_wave, q_tide])
    dominant = np.array(['river', 'wave', 'tide'])[np.argmax(flux, axis=0)]
    swept = m['strait'] & (STRAIT_SPEED > geo['velocity'])
    dominant = np.where(swept, 'tide', dominant)
    out.update(q_river_kg_s=q_river, q_wave_kg_s=q_wave, q_tide_kg_s=q_tide, tidal_discharge=qw_tide,
               tide_ratio=qw_tide / q_mean, tide_ratio_gentle=gentle / q_mean,
               dominant=dominant, strait_swept=swept, wave_h=h0, wave_t=t0)
    # Scales.
    backwater = geo['depth'] / geo['slope']
    dens = f_bed * DRY_DENSITY['sand'] + (1.0 - f_bed) * DRY_DENSITY['mud']
    # Channels switch within the backwater reach on low-gradient deltas (Ganti et al. 2016) and over the delta's
    # own extent on steep fan deltas, whose backwater reach is metres.
    avulsion_length = np.maximum(backwater, r_delta)
    avulsion = P.avulsion_time(geo['depth'], BELT_WIDTHS * geo['width'], avulsion_length, kept_kg_yr, dens)
    bmode = P.bar_mode(geo['width'], geo['slope'], q_bf, d50, r, CF_SAND, g, np.where(d50 >= GRAVEL_D, 10.0, 4.0))
    pattern = np.where(bmode < 1.5, 'single-thread', np.where(bmode < 2.5, 'transitional', 'braided'))
    out.update(mouth_bar_km=P.mouth_bar_distance(geo['width'], geo['velocity'], d50, r, g) / 1e3,
               backwater_km=backwater / 1e3, avulsion_yr=avulsion, aggradation_m_yr=geo['depth'] / avulsion,
               bar_mode=bmode, pattern=pattern, tide_flood_km=m['tide_range'] / geo['slope'] / 1e3,
               underflow_channels=shares['bankfull']['plunging'] >= 0.5)
    return out


# --- Aggregates -------------------------------------------------------------------------------------------------------

def share_by(values, cats, weights):
    tot = weights.sum()
    return {str(c): float(weights[values == c].sum() / tot) if tot > 0 else 0.0 for c in cats}


def quantiles(x, w=None, qs=(0.1, 0.5, 0.9)):
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    if not ok.any():
        return None
    x = x[ok]
    w = np.ones_like(x) if w is None else np.asarray(w, float)[ok]
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    return [float(x[o][min(np.searchsorted(c, q), x.size - 1)]) for q in qs]


def summarise(m, out, sel=None):
    """Counts, flux-weighted shares and quantiles (10th, 50th, 90th) over the mouths selected: geometry and scales
    by mouth, the keep share and the shaping balance by load. Discharge weights would count a river once for every
    lake of a chain it runs through, so shares by load stand beside shares by mouth."""
    sel = np.ones(m['lat'].size, bool) if sel is None else sel
    if not sel.any():
        return dict(mouths=0)
    rock = out['rock_mt_yr'][sel]
    one = np.ones(rock.size)
    pick = {k: (v[sel] if isinstance(v, np.ndarray) and v.shape[:1] == sel.shape else v) for k, v in out.items()}
    plume = lambda flow, w: {k: float((pick[f'plume_{flow}_{k}'] * w).sum() / w.sum())
                             for k in ('floating', 'level', 'plunging', 'to_floor')}
    q = pick['q_mean']
    return dict(
        mouths=int(sel.sum()), discharge_m3s=float(q.sum()), rock_mt_yr=float(rock.sum()),
        dominance_by_count=share_by(pick['dominant'], ('river', 'wave', 'tide'), one),
        dominance_by_load=share_by(pick['dominant'], ('river', 'wave', 'tide'), rock),
        tide_ratio=quantiles(pick['tide_ratio']), tide_ratio_gentle_slope=quantiles(pick['tide_ratio_gentle']),
        tide_dominated_if_slope_tenfold_gentler_by_load=float((rock * (pick['tide_ratio_gentle'] > np.maximum(
            1.0, pick['q_wave_kg_s'] / np.maximum(pick['q_river_kg_s'], 1e-30)))).sum() / rock.sum()),
        wave_over_river=quantiles(pick['q_wave_kg_s'] / np.maximum(pick['q_river_kg_s'], 1e-30)),
        plume_mean_flow_by_load=plume('mean', rock), plume_bankfull_by_load=plume('bankfull', rock),
        plume_mean_flow_by_mouth=plume('mean', one), plume_bankfull_by_mouth=plume('bankfull', one),
        plunge_depth_m=quantiles(pick['plume_bankfull_plunge_depth']),
        density_excess_kg_m3=quantiles(pick['plume_mean_drho_mean']),
        pattern_by_mouth=share_by(pick['pattern'], ('single-thread', 'transitional', 'braided'), one),
        pattern_by_load=share_by(pick['pattern'], ('single-thread', 'transitional', 'braided'), rock),
        keep_share_by_load=float((pick['keep_share'] * rock).sum() / rock.sum()),
        keep_share_quantiles_by_load=quantiles(pick['keep_share'], rock),
        keep_with_underflow_by_load=float((pick['keep_with_underflow'] * rock).sum() / rock.sum()),
        bed_share_by_load=float((pick['bed_share'] * rock).sum() / rock.sum()),
        floc_share=quantiles(pick['floc_share'])[1],
        mean_discharge_m3s=quantiles(q), bankfull_discharge_m3s=quantiles(pick['q_bf']),
        width_m=quantiles(pick['width']), depth_m=quantiles(pick['depth']),
        velocity_m_s=quantiles(pick['velocity']), slope=quantiles(pick['slope']),
        bankfull_stress_pa=quantiles(pick['stress']), bank_critical_pa=quantiles(pick['bank_pa']),
        bed_d50_um=[v * 1e6 for v in quantiles(pick['d50_bed'])],
        loose_width_m=quantiles(pick['loose_width']), loose_depth_m=quantiles(pick['loose_depth']),
        loose_velocity_m_s=quantiles(pick['loose_velocity']), loose_slope=quantiles(pick['loose_slope']),
        loose_stress_pa=quantiles(pick['loose_stress']),
        threshold_over_loose_depth=quantiles(pick['depth'] / pick['loose_depth']),
        threshold_over_loose_width=quantiles(pick['width'] / pick['loose_width']),
        taller_than_steep_banks_by_mouth=float(pick['narrow'].mean()),
        capacity_ratio_by_load=quantiles(pick['capacity_ratio'], rock),
        concentration_kg_m3=quantiles(pick['conc_kg_m3']),
        mouth_bar_km=quantiles(pick['mouth_bar_km']), backwater_km=quantiles(pick['backwater_km']),
        avulsion_yr=quantiles(pick['avulsion_yr']), avulsion_yr_by_load=quantiles(pick['avulsion_yr'], rock),
        aggradation_m_yr_by_load=quantiles(pick['aggradation_m_yr'], rock),
        tide_flood_km=quantiles(np.where(pick['tide_flood_km'] > 0, pick['tide_flood_km'], np.nan)),
        rossby_km=quantiles(pick['rossby_km']), plume_reach_km=quantiles(pick['plume_reach_km']),
        reach_over_rossby=quantiles(pick['plume_reach_km'] / np.maximum(pick['rossby_km'], 1e-9)),
        coastal_current_by_load=float((pick['coastal_current'] * rock).sum() / rock.sum()),
        underflow_channels_by_mouth=float(pick['underflow_channels'].mean()),
        underflow_channels_by_load=float((pick['underflow_channels'] * rock).sum() / rock.sum()),
        strait_swept=int(pick['strait_swept'].sum()))


# --- Earth deltas -----------------------------------------------------------------------------------------------------

# Earth deltas of known type (recalled typical values, not re-read in this study; see sources.json): mean discharge,
# sediment load, bed-material texture, bank stress, waves, tide, receiving water and temperatures.
EARTH_DELTAS = dict(
    mississippi=dict(width_obs=1000.0, slope_obs=3e-05, pattern_obs='single-thread', type='river', plume='floating', q=18400.0, load_mt=200.0, texture=(0.30, 0.50, 0.20),
                     lat=29.0, wave_h=1.0, wave_t=6.0, wave_dir=0.5, tide_range=0.4, tide_period=24.84 * 3600.0,
                     salt=35.0, t_river=26.0, t_surface=29.0, t_deep=20.0, depth=60.0, delta_km2=25000.0,
                     bank_pa=10.0, ratio=2.0, haline=HALINE_SEAWATER),
    sao_francisco=dict(width_obs=800.0, slope_obs=0.0001, pattern_obs=None, type='wave', plume='floating', q=2800.0, load_mt=6.0, texture=(0.05, 0.15, 0.80),
                       lat=-10.5, wave_h=1.5, wave_t=7.0, wave_dir=0.7, tide_range=2.0, tide_period=12.42 * 3600.0,
                       salt=35.0, t_river=27.0, t_surface=27.0, t_deep=20.0, depth=50.0, delta_km2=800.0,
                       bank_pa=5.0, ratio=3.0, haline=HALINE_SEAWATER),
    fly=dict(width_obs=1500.0, slope_obs=2e-05, pattern_obs='single-thread', type='tide', plume='floating', q=6000.0, load_mt=100.0, texture=(0.30, 0.50, 0.20), lat=-8.5,
             wave_h=0.7, wave_t=6.0, wave_dir=0.5, tide_range=3.5, tide_period=12.42 * 3600.0, salt=35.0,
             t_river=27.0, t_surface=28.0, t_deep=20.0, depth=30.0, delta_km2=8000.0, bank_pa=10.0, ratio=2.0,
             haline=HALINE_SEAWATER),
    rhone_geneva=dict(width_obs=100.0, slope_obs=0.0015, pattern_obs=None, type='river', plume='plunging', q=300.0, load_mt=2.5, texture=(0.10, 0.50, 0.40), lat=46.4,
                      wave_h=0.3, wave_t=3.0, wave_dir=0.5, tide_range=0.0, tide_period=12.42 * 3600.0, salt=0.2,
                      t_river=9.0, t_surface=20.0, t_deep=6.0, depth=150.0, delta_km2=20.0, bank_pa=5.0, ratio=2.0,
                      haline=HALINE))


def earth_deltas():
    env = dict(EARTH)
    rows = {}
    for name, d in EARTH_DELTAS.items():
        clay, silt, sand = d['texture']
        f = sum(s * P.lognormal_fractions(EDGES, *MODES[k]) for s, k in ((clay, 'clay'), (silt, 'silt'), (sand, 'sand')))
        f = f / f.sum()
        one = lambda v: np.array([v], dtype=float)
        m = dict(lat=one(d['lat']), kind=np.array(['sea' if d['salt'] > 5 else 'lake']), humid=np.array([True]),
                 bankfull_ratio=one(d['ratio']), wave_h=one(d['wave_h']), wave_t=one(d['wave_t']),
                 wave_dir=one(d['wave_dir']), tide_range=one(d['tide_range']), strait=np.array([False]),
                 depth_avail=one(d['depth']), delta_km2=one(d['delta_km2']), water_km2=one(1e6), earth_waves=True,
                 recv_salt=one(d['salt']), recv_div=one(60.0 if d['salt'] > 5 else 1.5), soda_na=one(0.0),
                 bank_cohesion=one(BANK_COHESION_PA['mature'][1]))
        env['tide_period'] = d['tide_period']
        temps = (np.full((1, 10), d['t_river']), np.full((1, 10), d['t_surface']), one(d['t_deep']))
        waters = dict(river_salt=0.2, river_div=1.0, organic=(0.3, 0.3, 0.3), sea_div=60.0, lake_div=1.5,
                      soda_salt=(0, 0, 0), sea_salt=d['salt'], lake_salt=d['salt'])
        o = screen(m, 'mature', 'loam', dict(env, haline_sea=d['haline']), one(d['load_mt'] * 1e9), one(d['q']),
                   organic=0.3, waters=waters, temps=temps, tau_bank=one(d['bank_pa']), fractions=f, rho_s=2650.0)
        plume = max(('floating', 'level', 'plunging'), key=lambda k: float(o[f'plume_mean_{k}'][0]))
        # The same fluxes and bar criterion on the observed channel, to separate the framework from the geometry.
        qw_obs = P.tidal_discharge(0.5 * d['tide_range'], d['width_obs'], d['slope_obs'], d['tide_period'])
        flux_obs = (float(o['q_river_kg_s'][0]), float(o['q_wave_kg_s'][0]), qw_obs / d['q'] * float(o['q_river_kg_s'][0]))
        m_obs = float(P.bar_mode(d['width_obs'], d['slope_obs'], d['q'] * d['ratio'], float(o['d50_bed'][0]), 1.65,
                                 CF_SAND, env['g'], 4.0))
        rows[name] = dict(expected_type=d['type'], screened_type=str(o['dominant'][0]), expected_plume=d['plume'],
                          screened_plume=plume, q_river_kg_s=float(o['q_river_kg_s'][0]),
                          q_wave_kg_s=float(o['q_wave_kg_s'][0]), q_tide_kg_s=float(o['q_tide_kg_s'][0]),
                          width_m=float(o['width'][0]), depth_m=float(o['depth'][0]), slope=float(o['slope'][0]),
                          backwater_km=float(o['backwater_km'][0]), avulsion_yr=float(o['avulsion_yr'][0]),
                          rossby_km=float(o['rossby_km'][0]), plume_reach_km=float(o['plume_reach_km'][0]),
                          bar_mode=float(o['bar_mode'][0]), keep_share=float(o['keep_share'][0]),
                          pattern=str(o['pattern'][0]), observed_width_m=d['width_obs'], observed_slope=d['slope_obs'],
                          type_on_observed_channel=('river', 'wave', 'tide')[int(np.argmax(flux_obs))],
                          q_tide_on_observed_channel_kg_s=flux_obs[2], bar_mode_on_observed_channel=m_obs,
                          expected_pattern=d['pattern_obs'], plume_mean_plunging=float(o['plume_mean_plunging'][0]),
                          plunge_depth_m=float(o['plume_bankfull_plunge_depth'][0]))
    return rows


def seawater_plunge_threshold(t_c=25.0, salt=35.0, rho_s=2650.0):
    """Sediment (kg/m3) that makes fresh river water as dense as seawater at the same temperature."""
    drho = P.density(t_c, salt, 0.0, HALINE_SEAWATER) - P.density(t_c, 0.2, 0.0)
    return float(drho / (1.0 - P.water_density(t_c) / rho_s))


# --- The run ----------------------------------------------------------------------------------------------------------

def run():
    started = time.time()
    dc = design_climate()
    m, checks = load_mouths()
    waters = attach_waters(m, dc)
    attach_waves(m)
    lapse = json.loads(INPUTS['cold_waters'].read_text())['lapse_k_per_km']['central']
    drainage = json.loads(INPUTS['drainage'].read_text())
    n = m['lat'].size
    # The flows: the drainage's older climate, bracketed by the design years' land runoff against the old run's.
    ratio = dc['runoff_mm_day'] / old_runoff_mm_day()
    flow = dict(low=min(0.75, ratio), central=1.0, high=max(1.25, ratio))
    # The young Moon's drier rivers: the crust keeps the difference between the mature and the 1,000-year
    # groundwater seepage in each band, against the rivers' flow there.
    gw = json.loads(INPUTS['groundwater'].read_text())['cases']['central']['states']
    gw_bands = ('0-15', '15-30', '30-45', '45-60', '60-75', '75-90')
    gband = np.minimum((np.abs(m['lat']) // 15).astype(int), 5)
    to_common = np.isin(m['kind'], ('sea', 'small'))         # each river reaches the common level once
    river_km3 = np.bincount(gband[to_common], weights=m['q'][to_common] * YEAR_S / 1e9, minlength=6)
    young_flows = {}
    for stage in ('young_1000_yr', 'young_10000_yr'):
        loss = np.array([gw['mature']['by_band'][b]['seepage_km3_yr'] - gw[stage]['by_band'][b]['seepage_km3_yr']
                         for b in gw_bands])
        young_flows[stage] = np.clip(1.0 - loss / np.maximum(river_km3, 1e-9), 0.1, 1.0)
    young_flow = young_flows['young_1000_yr']
    checks.update(
        sea_discharge_m3s=float(m['q'][m['kind'] == 'sea'].sum()),
        common_level_discharge_m3s=float(m['q'][np.isin(m['kind'], ('sea', 'small'))].sum()),
        drainage_discharge_to_sea_m3s=float(drainage['rivers']['discharge_to_sea_m3s']),
        lake_mouths_over_inflow=quantiles(m['q'][np.isfinite(m['unit_inflow'])] * YEAR_S
                                          / np.maximum(m['unit_inflow'][np.isfinite(m['unit_inflow'])], 1.0)),
        largest_rivers=largest_check(m, drainage))
    temps = {}
    for key, k in (('warm_rivers', 0), ('central', 1), ('cool_rivers', 2)):
        temps[key] = temperatures(m, dc, lapse, CATCHMENT_SHARE[k], LAKE_OFFSET_K[k], SEA_SURFACE_C[k],
                                  SEA_DEEP_C[k])
    results, tables = {}, {}
    for state in ('young', 'mature'):
        soils = ('regolith',) if state == 'young' else ('sand', 'sandy_loam', 'loam', 'clay_loam', 'clay')
        q_state = m['q'] * (young_flow[gband] if state == 'young' else 1.0)
        for soil in soils:
            res = {}
            for case, gt in TERRAIN_GT_YR.items():
                factor = YOUNG_YIELD[case] if state == 'young' else 1.0
                year = 1000 if state == 'young' else 10000
                m['delta_km2'] = m[f'delta_{case}_{year}'] * factor
                rock = m['rock_per_gt'] * gt * factor
                out = screen(m, state, soil, MOON, rock, q_state, temps=temps['central'])
                res[case] = by_water(m, out, waters['names'], brief=case != 'central')
                if case == 'central':
                    tables[(state, soil)] = out
                    for fc in ('low', 'high'):
                        o2 = screen(m, state, soil, MOON, rock, q_state * flow[fc], temps=temps['central'])
                        res[f'flow_{fc}'] = by_water(m, o2, waters['names'], brief=True)
                    for tk in ('warm_rivers', 'cool_rivers'):
                        o3 = screen(m, state, soil, MOON, rock, q_state, temps=temps[tk])
                        res[tk] = by_water(m, o3, waters['names'], brief=True)
                    if soil in ('regolith', 'sand', 'loam', 'clay'):
                        bracket = {}
                        for label, geom in GEOMETRY_BRACKET.items():
                            o5 = screen(m, state, soil, MOON, rock, q_state, temps=temps['central'], geometry=geom)
                            bracket[label] = by_water(m, o5, waters['names'], brief=True)
                        res['channel_bracket'] = bracket
                    for k, label in ((0, 'banks_weak'), (2, 'banks_strong')):
                        o4 = screen(m, state, soil, MOON, rock, q_state, temps=temps['central'],
                                    root=ROOT_FACTOR[state][k])
                        res[label] = by_water(m, o4, waters['names'], brief=True)
            results[f'{state}/{soil}'] = res
    # The lunar mouths with Earth's constants, seawater in the seas and Earth's waves.
    earth_on_moon = {}
    for state, soil in (('young', 'regolith'), ('mature', 'loam')):
        m['delta_km2'] = m[f'delta_central_{1000 if state == "young" else 10000}'] * (
            YOUNG_YIELD['central'] if state == 'young' else 1.0)
        rock = m['rock_per_gt'] * TERRAIN_GT_YR['central'] * (YOUNG_YIELD['central'] if state == 'young' else 1.0)
        q_state = m['q'] * (young_flow[gband] if state == 'young' else 1.0)
        out = screen(m, state, soil, EARTH, rock, q_state, temps=temps['central'])
        earth_on_moon[f'{state}/{soil}'] = by_water(m, out, waters['names'], brief=True)
    agg = results['young/regolith']['central']['seas'].get('aggradation_m_yr_by_load')
    sea_level = dict(rates_m_yr=SEA_LEVEL_RATES_M_YR, young_sea_mouths_aggradation_m_yr_by_load=agg,
                     note=('The crust takes the seas\' water at a 17 cm a year fall if delivery stopped; the sag lowers '
                           'the deep floors 1.5 cm a year under the central mantle and up to 41 cm a year under the '
                           'warm profile in the first millennium (geography README). Where delivery lags, the young '
                           'deltas build into a moving sea.'))
    product = dict(
        schema=SCHEMA, producer=dict(files={c: sha256(ROOT / c) for c in CODE}), young_sea_level=sea_level,
        evidence=EVIDENCE, reading_rule=READING_RULE,
        inputs={k: dict(path=str(p), sha256=sha256(p)) for k, p in INPUTS.items()},
        climate=dict(flows_run='A28_dim5_moon years 20-29 (geography drainage)', temperatures_run=dc['run'],
                     design_land_runoff_mm_day=dc['runoff_mm_day'], flow_bracket=flow,
                     young_flow_factor_by_band={k: dict(zip(gw_bands, v.tolist())) for k, v in young_flows.items()},
                     rivers_to_common_level_km3_yr_by_band=dict(zip(gw_bands, river_km3.tolist())),
                     lapse_k_per_km=lapse, zone_by_band={b: CLASSES[c] for b, c in zip(BANDS, dc['band_class'])}),
        checks=checks, waters=dict(seas=waters['names'], sea_level_m=waters['sea_level'],
                                   mouths_by_kind={k: int((m['kind'] == k).sum())
                                                   for k in ('sea', 'small', 'lake', 'soda')},
                                   wave_sources={str(s): int((m['wave_source'] == s).sum())
                                                 for s in np.unique(m['wave_source'])}),
        states=results, earth_constants_on_lunar_mouths=earth_on_moon, earth_deltas=earth_deltas(),
        constants=constants_block(), largest_mouths=largest(m, tables),
        settings=settings_block(), seconds=time.time() - started)
    write_table(m, tables)
    product['peak_memory_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    product['seconds'] = time.time() - started
    return rnd(product)


def by_water(m, out, names, brief=False):
    kinds = dict(all=np.ones(m['lat'].size, bool), seas=m['kind'] == 'sea', small_waters=m['kind'] == 'small',
                 rain_fed_lakes=m['kind'] == 'lake', soda_lakes=m['kind'] == 'soda')
    res = {k: summarise(m, out, s) for k, s in kinds.items()}
    if not brief:
        res['by_sea'] = {names.get(int(lab), str(lab)): summarise(m, out, m['sea_label'] == lab)
                         for lab in np.unique(m['sea_label'][m['sea_label'] >= 0])}
        res['by_zone'] = {CLASSES[z]: summarise(m, out, m['zone'] == z) for z in np.unique(m['zone'])}
    else:
        keep = ('mouths', 'dominance_by_count', 'dominance_by_load', 'plume_mean_flow_by_load',
                'plume_bankfull_by_load', 'plume_mean_flow_by_mouth', 'keep_share_by_load', 'pattern_by_mouth',
                'pattern_by_load', 'width_m', 'depth_m', 'velocity_m_s', 'slope', 'avulsion_yr',
                'avulsion_yr_by_load', 'tide_ratio', 'tide_flood_km', 'backwater_km', 'mouth_bar_km', 'rossby_km',
                'plume_reach_km', 'coastal_current_by_load', 'density_excess_kg_m3')
        res = {k: {kk: v.get(kk) for kk in keep} for k, v in res.items()}
    return res


def largest_check(m, drainage):
    rows = []
    for r in drainage['rivers']['largest']:
        lat, lon = r['mouth']
        d = haversine_km(lat, lon, m['lat'], m['lon'])
        near = d < 30.0
        q = float(m['q'][near].max()) if near.any() else 0.0
        rows.append(dict(mouth=r['mouth'], drainage_m3s=r['discharge_m3s'], screened_m3s=q))
    return rows


def largest(m, tables, top=15):
    order = np.argsort(m['q'])[::-1][:top]
    try:
        from research.studies.terrain_through_time.routing import place_names
        names = place_names(m['lat'][order], m['lon'][order])
    except Exception:                                       # the gazetteer is an ignored input
        names = [None] * order.size
    rows = []
    for j, i in enumerate(order):
        row = dict(place=names[j], lat=float(m['lat'][i]), lon_east=float(m['lon'][i]), into=str(m['kind'][i]),
                   mean_discharge_m3s=float(m['q'][i]), tide_range_m=float(m['tide_range'][i]))
        for (state, soil), o in tables.items():
            if soil not in ('regolith', 'loam'):
                continue
            row[state] = dict(soil=soil, rock_mt_yr=float(o['rock_mt_yr'][i]), q_mean=float(o['q_mean'][i]),
                              width_m=float(o['width'][i]), depth_m=float(o['depth'][i]),
                              slope=float(o['slope'][i]), dominant=str(o['dominant'][i]),
                              plume_mean_plunging=float(o['plume_mean_plunging'][i]),
                              plume_mean_floating=float(o['plume_mean_floating'][i]),
                              keep_share=float(o['keep_share'][i]), pattern=str(o['pattern'][i]),
                              avulsion_yr=float(o['avulsion_yr'][i]), backwater_km=float(o['backwater_km'][i]),
                              rossby_km=float(o['rossby_km'][i]), plume_reach_km=float(o['plume_reach_km'][i]))
        rows.append(row)
    return rows


def constants_block():
    out = {}
    for name, env in (('moon', MOON), ('earth', EARTH)):
        out[name] = dict(
            g=env['g'], omega_s=env['omega'], tide_period_days=env['tide_period'] / JULIAN_DAY,
            mud_front_failure_m={f'{s:g}_deg': [float(P.front_failure_height(c, s, env['g'], MUD['rho_sat']))
                                                for c in MUD['cohesion_pa']] for s in MUD['slopes_deg']},
            sand_front_slope_deg=SAND_REPOSE_DEG,
            settling_mm_s={f'{d * 1e6:g}_um': float(P.settling_velocity(d, env['g'], 2700.0, 0.9e-6) * 1e3)
                           for d in (2e-6, 10e-6, 62.5e-6, 200e-6, 1e-3)},
            floc_settling_mm_s=[v * env['g'] / STANDARD_GRAVITY * 1e3 for v in FLOC_SETTLING_EARTH])
    out['seawater_plunge_threshold_kg_m3'] = seawater_plunge_threshold()
    out['sea_level_rates_m_yr'] = SEA_LEVEL_RATES_M_YR
    return out


def settings_block():
    return dict(young_yield=YOUNG_YIELD, terrain_gt_yr=TERRAIN_GT_YR, bankfull_ratio=BANKFULL_RATIO,
                intermittency=INTERMITTENCY, cf_sand=CF_SAND, rouse_split=ROUSE_SPLIT, root_factor=ROOT_FACTOR,
                waters=WATERS, ccc_divalent_mmol=CCC_DIVALENT_MMOL, ccc_soda_na_mmol=CCC_SODA_NA_MMOL,
                floc_settling_earth_m_s=FLOC_SETTLING_EARTH, sea_surface_c=SEA_SURFACE_C, sea_deep_c=SEA_DEEP_C,
                lake_offset_k=LAKE_OFFSET_K, catchment_share=CATCHMENT_SHARE, wind_10m=WIND_10M,
                direction_resultant=DIRECTION_RESULTANT, strait_speed=STRAIT_SPEED, mud=MUD,
                belt_widths=BELT_WIDTHS, dry_density=DRY_DENSITY, zone_soils=ZONE_SOILS,
                bank_cohesion_pa=BANK_COHESION_PA)


def write_table(m, tables):
    DRIVE.mkdir(parents=True, exist_ok=True)
    cols = ('q_mean', 'q_bf', 'rock_mt_yr', 'width', 'depth', 'velocity', 'slope', 'stress', 'loose_width',
            'loose_depth', 'loose_slope', 'd50_bed', 'bed_share', 'conc_kg_m3', 'floc_share', 'plume_mean_floating',
            'plume_mean_level', 'plume_mean_plunging', 'plume_bankfull_plunging', 'plume_bankfull_plunge_depth',
            'rossby_km', 'plume_reach_km', 'q_river_kg_s', 'q_wave_kg_s', 'q_tide_kg_s', 'dominant', 'keep_share',
            'mouth_bar_km', 'backwater_km', 'avulsion_yr', 'bar_mode', 'pattern', 'tide_flood_km')
    keys = [k for k in tables if k[1] in ('regolith', 'loam')]
    with gzip.open(TABLE, 'wt', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['mouth_cell', 'lat', 'lon_east', 'into', 'sea_label', 'zone', 'mean_discharge_m3s',
                     'grid_slope', 'tide_range_m', 'wave_h_m', 'wave_t_s', 'wave_source']
                    + [f'{s}_{c}' for s, _ in keys for c in cols])
        for i in range(m['lat'].size):
            row = [int(m['mouth_cell'][i]), f"{m['lat'][i]:.4f}", f"{m['lon'][i]:.4f}", m['kind'][i],
                   int(m['sea_label'][i]), CLASSES[m['zone'][i]], f"{m['q'][i]:.4g}", f"{m['slope_grid'][i]:.3g}",
                   f"{m['tide_range'][i]:.3g}", f"{m['wave_h'][i]:.3g}", f"{m['wave_t'][i]:.3g}",
                   m['wave_source'][i]]
            for key in keys:
                o = tables[key]
                row += [o[c][i] if isinstance(o[c][i], str) else f'{float(o[c][i]):.4g}' for c in cols]
            wr.writerow(row)


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    product = run()
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(OUT, TABLE)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
