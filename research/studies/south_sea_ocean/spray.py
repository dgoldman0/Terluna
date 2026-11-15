"""The salt a south sea at seawater's strength would send into the air, and where it would come down.

    /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.south_sea_ocean.spray
                                                              # results/spray.json; about five minutes, under 1 GB

A fresh sea's spray carries almost no salt; a salt sea's carries 35 g for each kilogram of drops. This module is a
first reading of how much, on the design years' own winds, and of where the rain and the ground take it out:

- the source: the sea-spray function the aerosol study's marine note assembled (Gong 2003 as printed by Jaegle et al.
  2011, with their weight for the water's warmth), driven by the 10 m wind of the corrected design run's year 10 as
  written every 6.5 hours (the record the littoral study measured its wind raise on). Lunar gravity enters as that
  note proposes: the wind is replaced by the Earth wind that puts the same energy into the waves (the same air
  density times friction velocity cubed), and the yield at that energy is carried as a low, central and high case,
  with the note's gentler wind law as a fourth;
- the fate: each 6.5 hours' salt from each of the sea's model cells is followed as a parcel on the record's winds of
  the lowest two levels (the lowest 1.9 km) round the sphere, losing salt to the ground at a deposition speed plus
  each size's own settling at lunar gravity, and to the rain at a scavenging ratio times the record's rain in the
  cell under it. Three cases of removal run from slow to fast. What reaches the ground is kept by model cell, wet and
  dry, with the parcels' time in the air.

It is a diagnostic on a T21 record: the model has no coast finer than 170 km, the parcels do not spread or rise out
of the lowest levels, and the spray of the surf and the coarsest drops, which fall within tens of kilometres, are not
in the source. Every lunar factor of the source is an extrapolation from Earth scalings (the marine note,
`atmospheric-electricity/aerosol_marine.md` on the research drive, sections 1.6 and 7.1).
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from research.studies.coast_water import load as L
from research.studies.south_sea_ocean import run as R
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'spray.json'
MAPS = ROOT / 'research' / 'runs' / 'south_sea_ocean' / 'spray_maps.npz'
SCHEMA = 'terluna.research.south-sea-spray/1'
RECORD_RUN = R.RUN + '_eddies'               # year 0 is the design run's year 10, written every 13 steps
RECORD = Path('/media/projectspace/terluna-research/gcm-runs') / RECORD_RUN / 'model' / 'MOST.00000.nc'
CODE = ('research/studies/south_sea_ocean/spray.py', 'research/studies/south_sea_ocean/run.py',
        'research/studies/coast_water/load.py', 'climate/gcm/turned.py')
WRITE_S = 13 * 1800.0                        # the record's interval, 6.5 hours

# --- the source ---------------------------------------------------------------------------------------------------
KARMAN, CHARNOCK = 0.41, 0.0144              # the marine note's log profile and Charnock roughness (Hayes et al. 2013)
AIR_VISCOSITY = 1.85e-5                      # Pa s, the marine note's for both airs
EARTH_AIR = 1.18                             # kg/m3, the marine note's Earth air
THETA = 30.0                                 # Gong's shape parameter
R80_EDGES_UM = (0.035, 0.5, 1.0, 2.0, 4.0, 7.0, 10.0, 14.0, 20.0)   # radius at 80% humidity; Gong's range
SALT_DENSITY = 2165.0                        # kg/m3, dry sea salt taken as halite's
WET_DENSITY = 1200.0                         # kg/m3 at 80% humidity (the marine note, section 6)
SST_FIT_TOP_C = 30.0                         # Jaegle et al.'s data end here; the weight is held at this value above
_IONS = R.ion_grams_per_kg()
CHLORIDE_SHARE = _IONS['cl'] / sum(_IONS.values())    # 0.55 of sea salt's mass
LAWS = dict(
    earth_law='Gong\'s function at the lunar wind itself, with no lunar change',
    low='the central case times 0.34: the film drops\' yield for each unit of wave energy at a sixth of the gravity',
    central='Gong\'s function at the Earth wind of equal energy flux into the waves',
    high='the central case times 2.46: air entrained by slower breakers at a sixth of the gravity',
    gentle_wind_law='Jaegle et al.\'s fitted wind law, 25.5 U^2.07 in place of 1.373 U^3.41, at the equal-energy wind')
MULTIPLIER = dict(low=0.34, central=1.0, high=2.46)

# --- the fate -----------------------------------------------------------------------------------------------------
LEVELS = 2                                   # the lowest levels whose mean wind carries the salt (to 1.9 km)
DEPTH_WATER_M, DEPTH_LAND_M = 1700.0, 5000.0  # the mixed layers CM1's rings give over water and over land by day
FATES = dict(                                # scavenging ratio by mass, and deposition speed besides settling (m/s)
    slow=dict(scavenging_ratio=100.0, deposition_m_s=0.0001),
    central=dict(scavenging_ratio=300.0, deposition_m_s=0.0005),
    fast=dict(scavenging_ratio=1000.0, deposition_m_s=0.003))
SUBSTEPS = 2
DROP = 1.0e-3                                # a parcel is ended when this share of its salt is left
MAX_YEARS = 6                                # the year's winds repeat until the parcels are spent
YEAR_S = JULIAN_YEAR_DAYS * JULIAN_DAY


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def gong_shape(r80_um):
    """Gong's (2003) spectrum without its wind factor: particles per m2 per s per um of r80 for each unit of
    1.373 U^3.41 (Jaegle et al. 2011, Eq. 1)."""
    r = np.asarray(r80_um, float)
    a = 4.7 * (1.0 + THETA * r) ** (-0.017 * r ** -1.44)
    b = (0.433 - np.log10(r)) / 0.433
    return r ** -a * (1.0 + 0.057 * r ** 3.45) * 10.0 ** (1.607 * np.exp(-b * b))


def bins(edges_um=R80_EDGES_UM, points=200):
    """For each size bin, the number and the dry salt mass (kg) it carries for each unit of the wind factor, and its
    mass-weighted mean r80 (um). The dry radius is half of r80."""
    number, mass, radius = [], [], []
    for lo, hi in zip(edges_um[:-1], edges_um[1:]):
        r = np.geomspace(lo, hi, points)
        n = gong_shape(r)
        m = n * SALT_DENSITY * 4.0 / 3.0 * math.pi * (0.5 * r * 1e-6) ** 3
        number.append(np.trapezoid(n, r))
        mass.append(np.trapezoid(m, r))
        radius.append(np.trapezoid(m * r, r) / np.trapezoid(m, r))
    return np.array(number), np.array(mass), np.array(radius)


def sst_weight(t_c):
    """Jaegle et al.'s (2011, Eq. 4) weight for the water's warmth, held at its value at the top of their data."""
    t = np.minimum(np.asarray(t_c, float), SST_FIT_TOP_C)
    return 0.3 + 0.1 * t - 0.0076 * t ** 2 + 0.00021 * t ** 3


def friction_velocity(u10, gravity, air_density):
    """Friction velocity (m/s) under a 10 m wind from the log profile over water whose roughness is Charnock's plus
    the smooth-flow term (Hayes et al. 2013, Eq. 2.4)."""
    u = np.maximum(np.asarray(u10, float), 1e-6)
    nu = AIR_VISCOSITY / air_density
    ust = 0.035 * u
    for _ in range(80):
        z0 = CHARNOCK * ust ** 2 / gravity + 0.11 * nu / ust
        ust = 0.5 * ust + 0.5 * KARMAN * u / np.log(10.0 / z0)
    return ust


def equivalent_wind(u10, gravity, air_density):
    """The Earth 10 m wind that puts the same energy flux into the waves, air density times friction velocity cubed."""
    grid = np.linspace(0.02, 80.0, 4000)
    earth = EARTH_AIR * friction_velocity(grid, STANDARD_GRAVITY, EARTH_AIR) ** 3
    here = air_density * friction_velocity(u10, gravity, air_density) ** 3
    return np.interp(here, earth, grid)


def wind_factor(u10, law, gravity=MOON_SURFACE_GRAVITY, air_density=1.404):
    """The source's wind factor (the multiple of Gong's shape) under each law."""
    u10 = np.asarray(u10, float)
    if law == 'earth_law':
        return 1.373 * u10 ** 3.41
    ueq = equivalent_wind(u10, gravity, air_density)
    if law == 'gentle_wind_law':
        return 25.5 * ueq ** 2.07
    return MULTIPLIER[law] * 1.373 * ueq ** 3.41


def settling(r80_um, gravity=MOON_SURFACE_GRAVITY, pressure_pa=121590.0, temperature_k=301.0):
    """Stokes settling speed (m/s) of a drop of radius r80 at 80% humidity, with Cunningham's slip correction."""
    d = 2.0 * np.asarray(r80_um, float) * 1e-6
    free_path = 0.0665e-6 * (101325.0 / pressure_pa) * (temperature_k / 293.0)
    kn = 2.0 * free_path / d
    slip = 1.0 + kn * (1.257 + 0.4 * np.exp(-1.1 / kn))
    return gravity * WET_DENSITY * d ** 2 * slip / (18.0 * AIR_VISCOSITY)


# --- the grid -----------------------------------------------------------------------------------------------------
def latitude_edges(glat):
    """The Gaussian rows' edges (degrees, north first), halfway between neighbours, with the poles."""
    g = np.asarray(glat, float)
    return np.concatenate([[90.0], 0.5 * (g[:-1] + g[1:]), [-90.0]])


def cell_area_m2(glat, nlon):
    e = np.radians(latitude_edges(glat))
    return MOON_RADIUS ** 2 * (np.sin(e[:-1]) - np.sin(e[1:])) * (2.0 * math.pi / nlon)


def cell_of(x, edges, lon0, nlon):
    """Row and column of unit vectors on the model grid."""
    lat = np.degrees(np.arcsin(np.clip(x[:, 2], -1.0, 1.0)))
    lon = np.degrees(np.arctan2(x[:, 1], x[:, 0]))
    i = np.clip(np.searchsorted(-edges, -lat, side='right') - 1, 0, edges.size - 2)
    j = np.rint(np.mod(lon - lon0, 360.0) / (360.0 / nlon)).astype(int) % nlon
    return i, j


def unit(lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)), axis=-1)


def step(x, u, v, dt):
    """Move unit vectors with an eastward and northward wind (m/s) for dt seconds along the great circle the wind
    points down."""
    lon = np.arctan2(x[:, 1], x[:, 0])
    sin_lat = np.clip(x[:, 2], -1.0, 1.0)
    cos_lat = np.sqrt(1.0 - sin_lat ** 2)
    east = np.stack((-np.sin(lon), np.cos(lon), np.zeros_like(lon)), axis=1)
    north = np.stack((-sin_lat * np.cos(lon), -sin_lat * np.sin(lon), cos_lat), axis=1)
    speed = np.hypot(u, v)
    angle = speed * dt / MOON_RADIUS
    heading = (u[:, None] * east + v[:, None] * north) / np.maximum(speed, 1e-30)[:, None]
    out = x * np.cos(angle)[:, None] + heading * np.sin(angle)[:, None]
    return out / np.linalg.norm(out, axis=1, keepdims=True)


def neighbours(mask):
    """The cells of a mask together with every cell touching one of them (the columns join across the seam)."""
    out = mask.copy()
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            shifted = np.roll(mask, dj, axis=1)
            if di == 1:
                shifted = np.vstack([np.zeros((1, mask.shape[1]), bool), shifted[:-1]])
            elif di == -1:
                shifted = np.vstack([shifted[1:], np.zeros((1, mask.shape[1]), bool)])
            out |= shifted
    return out


def follow(emission_kg, rows, cols, u, v, rain_kg_m2_s, water, glat, glon, settle_m_s, air_density, fate,
           seed=20261115, substeps=SUBSTEPS, drop=DROP, max_years=MAX_YEARS):
    """Follow each write's salt from each emitting cell to the ground. emission_kg[t, k, b] is the salt of size bin b
    leaving cell (rows[k], cols[k]) during write t. Returns the salt deposited in each model cell by the rain and
    dry, the time-integrated salt in the air over each cell (kg s), each bin's wet and dry totals, the mass-weighted
    mean distance between where salt rose and where it fell, and what was still in the air at the end."""
    nt, nk, nb = emission_kg.shape
    nlat, nlon = water.shape
    edges, lon0 = latitude_edges(glat), float(glon[0])
    depth = np.where(water, DEPTH_WATER_M, DEPTH_LAND_M)
    k_fix = (np.asarray(settle_m_s, float)[None, None, :] + fate['deposition_m_s']) / depth[:, :, None]
    wet_per_rain = fate['scavenging_ratio'] / (air_density * depth)
    rng = np.random.default_rng(seed)
    dt = WRITE_S / substeps
    wet = np.zeros(nlat * nlon)
    dry = np.zeros(nlat * nlon)
    airborne = np.zeros(nlat * nlon)
    by_bin = np.zeros((nb, 2))
    travel = 0.0
    x = np.zeros((0, 3))
    x0 = np.zeros((0, 3))
    m = np.zeros((0, nb))
    m0 = np.zeros(0)
    half_lon = 0.5 * 360.0 / nlon
    for n in range(nt * max_years):
        t = n % nt
        if n < nt:
            lat = rng.uniform(edges[rows + 1], edges[rows])
            lon = glon[cols] + rng.uniform(-half_lon, half_lon, nk)
            new = unit(lat, lon)
            x, x0 = np.concatenate([x, new]), np.concatenate([x0, new])
            m = np.concatenate([m, emission_kg[t]])
            m0 = np.concatenate([m0, emission_kg[t].sum(axis=1)])
        elif x.shape[0] == 0:
            break
        for _ in range(substeps):
            i, j = cell_of(x, edges, lon0, nlon)
            flat = i * nlon + j
            k_wet = (wet_per_rain[i, j] * rain_kg_m2_s[t, i, j])[:, None]
            k_dry = k_fix[i, j]
            total = k_wet + k_dry
            lost = m * -np.expm1(-total * dt)
            share_wet = k_wet / total
            # the salt in the air over the step: the integral of m exp(-k t)
            airborne += np.bincount(flat, weights=(lost / total).sum(axis=1), minlength=nlat * nlon)
            lw, ld = (lost * share_wet).sum(axis=1), (lost * (1.0 - share_wet)).sum(axis=1)
            wet += np.bincount(flat, weights=lw, minlength=nlat * nlon)
            dry += np.bincount(flat, weights=ld, minlength=nlat * nlon)
            by_bin[:, 0] += (lost * share_wet).sum(axis=0)
            by_bin[:, 1] += (lost * (1.0 - share_wet)).sum(axis=0)
            travel += float(((lw + ld) * np.arccos(np.clip((x * x0).sum(axis=1), -1.0, 1.0))).sum()) * MOON_RADIUS
            m = m - lost
            x = step(x, u[t, i, j], v[t, i, j], dt)
        left = m.sum(axis=1)
        spent = left < drop * m0
        if spent.any():                                 # the last thousandth falls where the parcel stands
            i, j = cell_of(x[spent], edges, lon0, nlon)
            dry += np.bincount(i * nlon + j, weights=left[spent], minlength=nlat * nlon)
            by_bin[:, 1] += m[spent].sum(axis=0)
            travel += float((left[spent] * np.arccos(np.clip((x[spent] * x0[spent]).sum(axis=1), -1.0, 1.0))).sum()
                            ) * MOON_RADIUS
            x, x0, m, m0 = x[~spent], x0[~spent], m[~spent], m0[~spent]
    emitted = float(emission_kg.sum())
    return dict(wet=wet.reshape(nlat, nlon), dry=dry.reshape(nlat, nlon), airborne=airborne.reshape(nlat, nlon),
                by_bin=by_bin, travel_m=travel / max(emitted, 1e-300), still_airborne_kg=float(m.sum()),
                emitted_kg=emitted)


# --- the record and the reading ----------------------------------------------------------------------------------
def read_record(path=RECORD, levels=LEVELS):
    import netCDF4
    with netCDF4.Dataset(path) as nc:
        get = lambda k: np.array(nc[k][:], dtype=np.float64)
        levp = get('levp')
        weight = np.diff(levp)[-levels:]
        ua = np.array(nc['ua'][:, -levels:], dtype=np.float64)
        va = np.array(nc['va'][:, -levels:], dtype=np.float64)
        w = (weight / weight.sum())[None, :, None, None]
        return dict(lat=get('lat'), lon=get('lon'), lsm=np.array(nc['lsm'][0], dtype=float), uas=get('uas'),
                    vas=get('vas'), u=(ua * w).sum(axis=1), v=(va * w).sum(axis=1), pr=get('pr'), ps=get('ps'),
                    tas=get('tas'), ts=get('ts'), top_sigma=float(levp[-levels - 1]))


def saline_years(chloride_kg_ha_yr):
    """Years for a chloride fall to bring a root zone's soil water to the coast-water study's saline threshold if
    nothing flushes it."""
    mol_m2_yr = chloride_kg_ha_yr * 0.1 / R.sc.M_CL
    return L.SALINE_MOL_M3 * L.CEMENT_POROSITY * L.ROOT_ZONE_M / mol_m2_yr


def quantiles(x, w, qs=(0.1, 0.5, 0.9)):
    x, w = np.asarray(x, float).ravel(), np.asarray(w, float).ravel()
    order = np.argsort(x)
    cum = np.cumsum(w[order])
    return [float(x[order][min(np.searchsorted(cum, q * cum[-1]), x.size - 1)]) for q in qs]


def build() -> dict:
    atlas = json.loads((ROOT / R.ATLAS).read_text())
    grid = np.load(ROOT / R.ATLAS_GRID)
    rotation = np.array(json.loads((ROOT / R.FRAME).read_text())['frame']['rotation_today_to_turned'])
    clim = dict(np.load(ROOT / R.CLIMATOLOGY))
    rec = read_record()
    glat, glon = rec['lat'], rec['lon']
    if not (np.allclose(glat, clim['lat']) and np.array_equal(rec['lsm'] > 0.5, clim['lsm'] > 0.5)):
        raise SystemExit('The 6.5-hour record is not on the design run\'s grid and coast')
    nt = rec['uas'].shape[0]
    water = rec['lsm'] < 0.5
    cell_m2 = np.repeat(cell_area_m2(glat, glon.size)[:, None], glon.size, axis=1)
    labels, alat, alon = grid['water_label'], grid['lat_deg'], grid['lon_deg']
    bodies = {b['label']: b for b in atlas['bodies']}
    seas = {name: R.on_model_grid(labels == label, alat, alon, rotation, glat, glon) * 1e6
            for label, name in R.SEAS.items()}
    own = seas[R.SEA_NAME] * water
    atlas_m2 = float(bodies[R.SEA_LABEL]['area_km2']) * 1e6
    own *= atlas_m2 / own.sum()                         # the sea's whole area, read at the model's water cells
    rows, cols = np.nonzero(own > 0)
    area = own[rows, cols]

    # the wind and the source
    speed = np.hypot(rec['uas'], rec['vas'])[:, rows, cols]
    rho = (rec['ps'] * 100.0 / (287.04 * rec['tas']))[:, rows, cols]
    air_density = float((rho * area).sum() / (area.sum() * nt))
    warmth = sst_weight(rec['ts'][:, rows, cols] - 273.15)
    number, mass, radius = bins()
    year_factor = YEAR_S / (nt * WRITE_S)               # the record's 357.5 days to a Julian year
    source = {}
    for law in LAWS:
        factor = wind_factor(speed, law, air_density=air_density) * warmth        # [t, k]
        kg = float((factor * area).sum() * WRITE_S * mass.sum()) * year_factor
        source[law] = dict(what=LAWS[law], salt_mt_yr=kg / 1e9, salt_g_m2_yr=kg * 1e3 / atlas_m2,
                           chloride_mt_yr=kg * CHLORIDE_SHARE / 1e9,
                           particles_m2_s=float((factor * area).sum() / (area.sum() * nt) * number.sum()))
    central_factor = wind_factor(speed, 'central', air_density=air_density) * warmth
    emission = central_factor[:, :, None] * area[None, :, None] * WRITE_S * mass[None, None, :]
    by_speed = {}
    for lo, hi in ((0, 2), (2, 3), (3, 4), (4, 5), (5, 99)):
        sel = (speed >= lo) & (speed < hi)
        by_speed[f'{lo}-{hi} m/s'] = dict(share_of_hours=float((sel * area).sum() / (area.sum() * nt)),
                                          share_of_salt=float((central_factor * sel * area).sum()
                                                              / (central_factor * area).sum()))
    wind = dict(mean_m_s=float((speed * area).sum() / (area.sum() * nt)),
                percentiles_m_s_50_90_99=quantiles(speed, np.broadcast_to(area, speed.shape), (0.5, 0.9, 0.99)),
                highest_m_s=float(speed.max()), air_density_kg_m3=air_density,
                equivalent_earth_wind_m_s={str(u): float(equivalent_wind(float(u), MOON_SURFACE_GRAVITY, air_density))
                                           for u in (2, 3, 5, 7, 10)},
                warmth_weight_mean=float((warmth * area).sum() / (area.sum() * nt)), by_speed=by_speed)

    # the fate, on the central source
    settle = settling(radius)
    rain = rec['pr'] * 1000.0                           # m/s of water to kg/m2/s
    land = ~water
    desert = land & (clim['pr_mm_day'] < 0.5)
    sea_share = {name: np.minimum(a / cell_m2, 1.0) for name, a in seas.items()}
    all_seas = np.minimum(sum(sea_share.values()), 1.0)
    beside = land & neighbours(seas[R.SEA_NAME] > 0)
    rain_year_kg_m2 = rain.mean(axis=0) * YEAR_S
    fates, maps = {}, {}
    for name, fate in FATES.items():
        t0 = time.time()
        f = follow(emission, rows, cols, rec['u'], rec['v'], rain, water, glat, glon, settle, air_density, fate)
        total = f['wet'] + f['dry']
        scale = year_factor                              # to a Julian year, like the source
        dep = total * scale / cell_m2                    # kg of salt per m2 a year
        wet_dep = f['wet'] * scale / cell_m2
        cl_ha = dep * CHLORIDE_SHARE * 1e4               # kg of chloride per hectare a year
        landish = 1.0 - all_seas                         # each cell's share outside the six seas
        to_seas = {s: float((total * sea_share[s]).sum() / total.sum()) for s in seas}
        in_rain = np.where(rain_year_kg_m2 > 0, wet_dep * CHLORIDE_SHARE / np.maximum(rain_year_kg_m2, 1e-30), 0.0)
        in_rain_mg_l = in_rain * 1e6                     # kg per kg of rain to mg per litre
        air_ug_m3 = f['airborne'] / (nt * WRITE_S) / cell_m2 / np.where(water, DEPTH_WATER_M, DEPTH_LAND_M) * 1e9
        w_land = cell_m2 * land
        fates[name] = dict(
            **fate, emitted_mt_yr=f['emitted_kg'] * scale / 1e9,
            deposited_over_emitted=float((total.sum() + f['still_airborne_kg']) / f['emitted_kg']),
            still_in_the_air_share=f['still_airborne_kg'] / f['emitted_kg'],
            share_by_rain=float(f['wet'].sum() / total.sum()),
            mean_days_in_the_air=float(f['airborne'].sum() / f['emitted_kg'] / JULIAN_DAY),
            mean_distance_km=f['travel_m'] / 1000.0,
            share_falling=dict(back_on_the_south_sea=to_seas[R.SEA_NAME],
                               on_the_other_seas={s: v for s, v in to_seas.items() if s != R.SEA_NAME},
                               on_land_and_lakes=float((total * landish).sum() / total.sum())),
            chloride_kg_ha_yr=dict(
                on_the_south_sea=float((cl_ha * seas[R.SEA_NAME]).sum() / seas[R.SEA_NAME].sum()),
                land_beside_the_sea=float((cl_ha * cell_m2 * beside).sum() / (cell_m2 * beside).sum()),
                land_mean=float((cl_ha * w_land).sum() / w_land.sum()),
                land_p10_p50_p90=quantiles(cl_ha[land], cell_m2[land]), land_highest=float(cl_ha[land].max()),
                desert_land_mean=float((cl_ha * cell_m2 * desert).sum() / (cell_m2 * desert).sum())),
            chloride_in_the_rain_mg_l=dict(
                land_mean=float((wet_dep * CHLORIDE_SHARE * w_land).sum() / (rain_year_kg_m2 * w_land).sum() * 1e6),
                land_p10_p50_p90=quantiles(in_rain_mg_l[land], cell_m2[land]),
                land_highest=float(in_rain_mg_l[land].max()),
                land_beside_the_sea=float((wet_dep * CHLORIDE_SHARE * cell_m2 * beside).sum()
                                          / (rain_year_kg_m2 * cell_m2 * beside).sum() * 1e6)),
            years_to_a_saline_root_zone_where_no_rain_flushes=dict(
                rule='coast_water.load: soil water at 4 dS/m, 40 mol of chloride per m3, in 0.3 m of soil water',
                land_beside_the_sea=saline_years(float((cl_ha * cell_m2 * beside).sum() / (cell_m2 * beside).sum())),
                desert_land_mean=saline_years(float((cl_ha * cell_m2 * desert).sum() / (cell_m2 * desert).sum()))),
            salt_in_the_air_ug_m3=dict(
                over_the_south_sea=float((air_ug_m3 * seas[R.SEA_NAME]).sum() / seas[R.SEA_NAME].sum()),
                land_beside_the_sea=float((air_ug_m3 * cell_m2 * beside).sum() / (cell_m2 * beside).sum()),
                land_mean=float((air_ug_m3 * w_land).sum() / w_land.sum()),
                land_highest=float(air_ug_m3[land].max())),
            by_size=[dict(r80_um=[float(a), float(b)], share_of_salt=float(mass[k] / mass.sum()),
                          settling_m_s=float(settle[k]),
                          share_by_rain=float(f['by_bin'][k, 0] / max(f['by_bin'][k].sum(), 1e-300)))
                     for k, (a, b) in enumerate(zip(R80_EDGES_UM[:-1], R80_EDGES_UM[1:]))],
            seconds=round(time.time() - t0, 1))
        maps[f'{name}_chloride_kg_ha_yr'] = cl_ha.astype(np.float32)
        maps[f'{name}_chloride_in_rain_mg_l'] = in_rain_mg_l.astype(np.float32)
        maps[f'{name}_salt_in_air_ug_m3'] = air_ug_m3.astype(np.float32)
    MAPS.parent.mkdir(parents=True, exist_ok=True)
    np.savez(MAPS, lat=glat, lon=glon, lsm=rec['lsm'], south_sea_m2=seas[R.SEA_NAME], **maps)

    ocean_salt_kg = R.salt(float(bodies[R.SEA_LABEL]['volume_km3']), 1.0)['salt_kg']
    central = fates['central']
    exported = {law: s['salt_mt_yr'] * (1.0 - central['share_falling']['back_on_the_south_sea'])
                for law, s in source.items()}
    return dict(
        schema=SCHEMA,
        evidence=('A diagnostic on one model year at T21: an Earth-fitted sea-spray source function carried to '
                  'lunar gravity by untested scalings (low, central and high cases), driven by 6.5-hour mean winds, '
                  'and parcels followed on the lowest levels\' winds with first-order removal. No measurement of '
                  'spray at reduced gravity exists. The surf\'s spray and drops above 20 um of r80, which fall near '
                  'the shore, are outside the source, and the model resolves no coast.'),
        reading_rule=('source gives the sea\'s yearly salt into the air under each law; the fate cases follow the '
                      'central source, and their shares apply to any multiple of it. Chloride is 55% of the salt. '
                      'chloride_kg_ha_yr and chloride_in_the_rain_mg_l are area-weighted over the model\'s land '
                      'cells (lakes included); land_beside_the_sea are the model\'s land cells that hold part of '
                      'the sea or touch a cell that does. salt_in_the_air_ug_m3 spreads the salt in the air over the '
                      'mixed layer '
                      '(1.7 km over water, 5 km over land). Multiply by the source cases\' ratios to the central '
                      'one for the low, high and gentle-law sources.'),
        record=dict(run=RECORD_RUN, what='the design run\'s year 10 written every 6.5 hours', writes=nt,
                    days=nt * WRITE_S / JULIAN_DAY, levels=LEVELS, top_sigma=rec['top_sigma']),
        sea=dict(name=R.SEA_NAME, area_km2=atlas_m2 / 1e6, model_cells=int(rows.size),
                 land_cells_beside_it=int(beside.sum()), land_beside_it_km2=float((cell_m2 * beside).sum() / 1e6)),
        wind=wind, source=source,
        sizes=dict(r80_edges_um=list(R80_EDGES_UM), share_of_salt=[float(x) for x in mass / mass.sum()],
                   share_of_number=[float(x) for x in number / number.sum()],
                   settling_m_s=[float(x) for x in settle]),
        fate=dict(depth_water_m=DEPTH_WATER_M, depth_land_m=DEPTH_LAND_M, cases=fates),
        against_the_sea=dict(
            salt_in_the_sea_kg=ocean_salt_kg,
            leaving_the_sea_for_good_mt_yr=exported,
            years_to_carry_off_a_thousandth={law: ocean_salt_kg * 1e-3 / (v * 1e9) for law, v in exported.items()}),
        chloride_share_of_salt=CHLORIDE_SHARE,
        inputs={**{p: sha256(ROOT / p) for p in (R.ATLAS, R.ATLAS_GRID, R.FRAME, R.CLIMATOLOGY)},
                str(RECORD): dict(bytes=RECORD.stat().st_size, sha256=sha256(RECORD))},
        producer=dict(files={p: sha256(ROOT / p) for p in CODE}))


def main(argv=None) -> int:
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    for law, s in out['source'].items():
        print(f"{law:16s} {s['salt_mt_yr']:8.2f} Mt of salt a year, {s['salt_g_m2_yr']:.3f} g/m2")
    for name, f in out['fate']['cases'].items():
        print(name, f"{f['mean_days_in_the_air']:.1f} days, {f['mean_distance_km']:.0f} km, by rain "
              f"{f['share_by_rain']:.2f}, back on the sea {f['share_falling']['back_on_the_south_sea']:.2f}, land "
              f"{f['chloride_kg_ha_yr']['land_mean']:.4f} kg Cl/ha/yr, rain {f['chloride_in_the_rain_mg_l']['land_mean']:.4f} "
              f"mg/L, {f['seconds']} s")
    return 0


if __name__ == '__main__':
    sys.exit(main())
