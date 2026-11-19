"""Wind-blown sand on the design surface: what moves it, where it could gather and a first map (19 November 2026).

    # once, with the climate environment (reads the run's NetCDF years; about a minute):
    ulimit -v 6000000; climate/gcm/.venv/bin/python -m research.studies.atlas_round.sand reduce
    # then, with the system Python (about four minutes on two capped CPUs, 4 GB):
    ulimit -v 6000000; OPENBLAS_NUM_THREADS=1 python3 -m research.studies.atlas_round.sand

The paper is sand.md beside this file. `reduce` keeps, from the design climate's run (sunlight 0.936, years 45-59), the
series a sand reckoning needs: each 3-day mean of the 10 m wind, the surface stress, the lowest level's wind, the air's
pressure and temperature at the ground, the Sun's height, the rain and the sensible heat. Everything after it reads
that file.

The reckoning has four steps, each with its evidence state in the record this module writes:

1. The physics (part A of the paper): Shao and Lu's (2000) fluid threshold and the flux law of Kok et al. (2012, their
   equation 2.34) at lunar gravity in the design's air, for grains of feldspar, glass and basalt.
2. The winds (part B): the run's 10 m wind is a wind over ground as rough as forest (PlaSim's land roughness, 2 m). It
   is carried to a sand surface at the lowest level's wind, and the part the run cannot hold (the turbulence of the
   sunlit hours, the circulations below 170 km, the storms' outflows) enters as a stated range taken from the
   cloud-resolving rings and the highland box. Each cell's yearly drift is summed by direction.
3. Supply and holding (parts C and D): the sand the land's own mantle, its rivers' ends and its dry shores offer, and
   the plants, crusts, lag and wetting that hold it.
4. The map (parts E and F): the offered sand is carried downwind on the 1.9 km grid as far as the drift takes it in
   the land's age, kept where the ground shelters it, lost where it meets water and held where it meets cover.

The climate run's grid is the turned Moon's; the product is on today's grid. Places turn with
`rotation_today_to_turned` (the design surface's frame), and a wind's direction turns with its place (`turn_azimuth`).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from climate.gcm.turned import to_today, to_turned
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'sand.json'
DRIVE = Path('/media/projectspace/terluna-research')
BULK = DRIVE / 'atlas-round' / 'sand'
SERIES = BULK / 'winds_A28D_y45-59.npz'
PRODUCT = BULK / 'sand_16ppd.npz'
RUN = DRIVE / 'gcm-runs' / 'A28D_301k_design_l15_snow_s936_settled'
YEARS = (45, 59)
SURFACE = DRIVE / 'atlas-round' / 'globe_surface' / 'globe_surface_16ppd.npz'
ROUTING = DRIVE / 'atlas-round' / 'design_routing' / 'design_routing_16ppd.npz'
ROCK = DRIVE / 'atlas-round' / 'fine_surface' / 'look_fields.npz'
REGIONS = DRIVE / 'atlas-round' / 'corrected_rule' / 'districts_redrawn' / 'districts_8ppd.npz'
REGION_RECORD = DRIVE / 'atlas-round' / 'corrected_rule' / 'districts_redrawn' / 'districts.json'
MOISTURE = HERE / 'results' / 'moisture_map_final.json'
SCHEMA = 'terluna.research.atlas-round-sand/1'
SERIES_SCHEMA = 'terluna.research.atlas-round-sand-winds/1'
CODE = ('research/studies/atlas_round/sand.py',)
PPD, ROWS, COLS = 16, 2880, 5760
SECONDS_PER_YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY

# ---- Part A: grains, air and the transport laws ---------------------------------------------------------------
G_MOON = MOON_SURFACE_GRAVITY            # m/s2
G_EARTH = STANDARD_GRAVITY
AIR_DENSITY = 1.35                       # kg/m3: the run's own surface air over land (ps / (R T), mean 1.354); the record uses 1.40-1.44 at sea level
AIR_DENSITY_EARTH = 1.225
VISCOSITY = 1.86e-5                      # Pa s at 29 C (the aerosol note's value)
VON_KARMAN = 0.4
SHAO_LU_A = 0.111                        # Shao and Lu (2000) as printed in Kok et al. (2012) eq. 2.8; fitted on Earth
SHAO_LU_GAMMA = (1.65e-4, 3.0e-4, 5.0e-4)   # N/m: Shao and Lu's range for dry loose grains, central 3e-4; fitted on Earth
IMPACT_RATIO = (0.5, 0.82, 1.0)          # impact over fluid threshold: Bagnold's 0.82 on Earth carried over; range from the aerosol note
FLUX_C = 5.0                             # C_DK of Kok et al. (2012) eq. 2.34, from Earth's 250 um sand (hop 0.1 m, speed gain 1 m/s)
GRAINS = {                               # particle densities, kg/m3 (handbook values, recalled; glass spans vesicular to dense)
    'feldspar': 2730.0, 'glass': 2500.0, 'basalt': 3000.0}
BULK_DENSITY = 1600.0                    # kg/m3 of loose dune sand (porosity about 0.4 for feldspar)
SAND_Z0 = (1e-4, 1e-3)                   # m: roughness of a sand surface at rest and in saltation, as the desert analysis brackets it
Z0_CENTRAL = 1e-3
MODEL_LAND_Z0 = 2.0                      # m: PlaSim's dz0land, the run's land roughness everywhere (landmod.f90; no surface file overrides it)
RING_LAND_Z0 = 0.1                       # m: the cloud-resolving cases' placeholder land (climate/crm/README.md)
DRAG_LENGTH_FACTOR = 2.2                 # saturation length over the drag length (rho_p / rho_a) d; Andreotti, Claudin and Pouliquen (2010), Earth and water-tank fit
ELEMENTARY_FACTOR = (12.0, 20.0)         # smallest dunes' wavelength over the saturation length (same source)
CREST_ASPECT = 20.0                      # crest spacing over dune height, Earth's dunes (Lancaster 1995, recalled)
MIXED_LAYER_M = (5100.0, 13200.0)        # the sunlit mixed layer: 5.1-5.2 km mean on the wet rings, 13.2 km on the design-ground ring

# ---- Part B: the run's winds and the part it cannot hold --------------------------------------------------------
LOWEST_LEVEL_M = 313.0                   # height of the run's lowest level (sigma 0.99415 in air near 300 K)
# The run's 10 m wind stands over a 2 m roughness. At the lowest level's wind the same air over a sand surface blows
# harder at 10 m. TO_SAND is that factor on the run's 10 m wind: low, central, high (see sand_wind_factor()).
SUBGRID = {                              # the unresolved part as a factor on the speed and a storm tail; see part B
    'none':    dict(scale=1.00, spread=None, storms=0.0, storm_gust=1.0),
    'low':     dict(scale=1.00, spread=4.0, storms=0.5, storm_gust=1.0),
    'central': dict(scale=1.30, spread=4.0, storms=1.0, storm_gust=1.25),
    'high':    dict(scale=2.00, spread=4.0, storms=2.0, storm_gust=1.5),
}
RING_TO_SAND = 1.2                       # a 10 m wind over the rings' 10 cm roughness carried to sand (1.15-1.27, see ring_to_sand())
STORM_GUST = dict(median=5.2, sigma_ln=0.174)   # per-storm highest 10 m wind on the wet 45 deg rings: median 5.2, 90th 6.4-6.7 m/s
STORM_TIME = 0.005                       # share of the time a land point lies in a storm's outflow at STORM_RAIN (rings: 0.6 storms a snapshot, 50 km of outflow each, over 5,000 km of land)
STORM_RAIN = 0.35                        # mm/day: the wet rings' land rain (0.31 and 0.40)
STORM_TIME_CAP = 0.05
DAY_COSINE = 0.1                         # a record counts as sunlit where its mean cosine of the Sun's zenith angle exceeds this
DIRECTIONS = 16


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def fluid_threshold(d, rho_p, rho_a=AIR_DENSITY, g=G_MOON, gamma=SHAO_LU_GAMMA[1], a_n=SHAO_LU_A):
    """Shao and Lu's (2000) friction velocity (m/s) at which wind alone starts grains of diameter d (m) moving."""
    d = np.asarray(d, float)
    return a_n * np.sqrt((rho_p - rho_a) / rho_a * g * d + gamma / (rho_a * d))


def easiest_grain(rho_p, rho_a=AIR_DENSITY, g=G_MOON, gamma=SHAO_LU_GAMMA[1], a_n=SHAO_LU_A):
    """The diameter (m) with the lowest fluid threshold and that threshold (m/s)."""
    d = math.sqrt(gamma / ((rho_p - rho_a) * g))
    return d, float(fluid_threshold(d, rho_p, rho_a, g, gamma, a_n))


def friction_velocity(u10, z0=Z0_CENTRAL):
    """u* (m/s) of a 10 m wind over roughness z0 on the neutral profile PlaSim uses, ln(1 + z / z0)."""
    return VON_KARMAN * np.asarray(u10, float) / math.log(1.0 + 10.0 / z0)


def wind_at_10m(ustar, z0=Z0_CENTRAL):
    return np.asarray(ustar, float) * math.log(1.0 + 10.0 / z0) / VON_KARMAN


def saltation_flux(ustar, ustar_it, rho_a=AIR_DENSITY, g=G_MOON, c=FLUX_C):
    """Mass of sand crossing a metre of width each second (kg/m/s): Q = C (rho_a / g) u*it (u*^2 - u*it^2) above the
    impact threshold, zero below it (Kok et al. 2012, eq. 2.34, after Duran et al. 2011)."""
    ustar = np.asarray(ustar, float)
    return np.where(ustar > ustar_it, c * rho_a / g * ustar_it * (ustar * ustar - ustar_it * ustar_it), 0.0)


def saturation_length(d, rho_p, rho_a=AIR_DENSITY):
    """The distance over which the flux adjusts to the wind, 2.2 drag lengths (m). Gravity does not enter it."""
    return DRAG_LENGTH_FACTOR * rho_p / rho_a * d


def hop(ustar, g=G_MOON, launch_factor=(0.5, 1.0)):
    """Height and length (m) of a ballistic hop without drag for grains leaving at launch_factor times u* upward and
    the same forward: an upper bound, since drag shortens it."""
    out = []
    for f in launch_factor:
        v = f * float(ustar) * 2.0   # lift-off speeds of order u* to 2 u* (Kok et al. 2012, their figure 2.12 scale)
        out.append((v * v / (2.0 * g), 2.0 * v * v / g))
    return out


def sand_wind_factor(z0_sand=Z0_CENTRAL, z0_model=MODEL_LAND_Z0, level=LOWEST_LEVEL_M, model_share=None):
    """The 10 m wind over sand as a multiple of the run's 10 m wind over its own roughness, the lowest level's wind
    held. With model_share None both profiles are neutral; otherwise model_share is the run's own 10 m wind over its
    lowest level's (0.52 at the median over land, unstable air bringing more wind down) and the sand surface takes the
    neutral share, so the factor is a lower one."""
    sand = math.log(1.0 + 10.0 / z0_sand) / math.log(1.0 + level / z0_sand)
    model = math.log(1.0 + 10.0 / z0_model) / math.log(1.0 + level / z0_model) if model_share is None else model_share
    return sand / model


def ring_to_sand(z0_sand=Z0_CENTRAL, z0_ring=RING_LAND_Z0, level=50.0):
    """The same for a 10 m wind over the cloud-resolving cases' 10 cm land, the wind at `level` held."""
    sand = math.log(1.0 + 10.0 / z0_sand) / math.log(1.0 + level / z0_sand)
    ring = math.log(1.0 + 10.0 / z0_ring) / math.log(1.0 + level / z0_ring)
    return sand / ring


TO_SAND = (1.4, 1.6, 2.06)               # low: the run's unstable share against a neutral sand profile; high: both neutral


# ---- the frame ----------------------------------------------------------------------------------------------------
def east_north(lat_deg, lon_deg):
    """Unit vectors to the east and to the north at a place, in the body frame its latitude and longitude belong to."""
    lat, lon = np.radians(np.asarray(lat_deg, float)), np.radians(np.asarray(lon_deg, float))
    east = np.stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)], axis=-1)
    north = np.stack([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat) * np.ones_like(lon)], axis=-1)
    return east, north


def vector_to_today(u_east, v_north, lat_t, lon_t, rotation):
    """A horizontal vector given by its eastward and northward parts at a turned-frame place, as a 3-vector in
    today's body frame (R^T applied to its turned-frame 3-vector)."""
    east, north = east_north(lat_t, lon_t)
    w = np.asarray(u_east, float)[..., None] * east + np.asarray(v_north, float)[..., None] * north
    return w @ np.asarray(rotation, float)


def azimuth_today(vec_today, lat_deg, lon_deg):
    """Eastward and northward parts, and the azimuth (degrees clockwise from north), of a body-frame 3-vector at a
    place on today's Moon."""
    east, north = east_north(lat_deg, lon_deg)
    e, n = (vec_today * east).sum(-1), (vec_today * north).sum(-1)
    return e, n, np.degrees(np.arctan2(e, n)) % 360.0


def turn_azimuth(azimuth_turned_deg, lat_t, lon_t, rotation):
    """A direction of travel given clockwise from the turned Moon's north at a turned-frame place, as the direction
    clockwise from today's north at the same point of the ground."""
    a = np.radians(np.asarray(azimuth_turned_deg, float))
    v = vector_to_today(np.sin(a), np.cos(a), lat_t, lon_t, rotation)
    lat, lon = to_today(lat_t, lon_t, rotation)
    return azimuth_today(v, lat, lon)[2]


# ---- the reduction (climate environment) ------------------------------------------------------------------------
SERIES_FIELDS = ('uas', 'vas', 'tauu', 'tauv', 'ps', 'tas', 'czen', 'pr', 'hfss', 'mrso')


def reduce_winds(log=print):
    """The run's design years as one file of 3-day means on its own grid (the turned Moon)."""
    import netCDF4
    BULK.mkdir(parents=True, exist_ok=True)
    out = {k: [] for k in SERIES_FIELDS + ('ual', 'val')}
    hashes, lat, lon, sigma = {}, None, None, None
    for year in range(YEARS[0], YEARS[1] + 1):
        path = RUN / 'model' / f'MOST.{year:05d}.nc'
        hashes[f'gcm-runs/{RUN.name}/model/{path.name}'] = sha256(path)
        with netCDF4.Dataset(path) as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            sigma = float(np.asarray(d['lev'][:], float).max())
            low = int(np.argmax(np.asarray(d['lev'][:], float)))
            for k in SERIES_FIELDS:
                out[k].append(np.asarray(d[k][:], np.float32))
            out['ual'].append(np.asarray(d['ua'][:, low], np.float32))
            out['val'].append(np.asarray(d['va'][:, low], np.float32))
            lsm = np.asarray(d['lsm'][:], np.float32).mean(axis=0)
        log(f'  year {year}')
    arrays = {k: np.concatenate(v) for k, v in out.items()}
    namelist = lambda name: (RUN / 'model' / name).read_text()
    meta = dict(schema=SERIES_SCHEMA, run=RUN.name, years=list(YEARS), records=int(arrays['uas'].shape[0]),
                frame='the run\'s own grid, the turned Moon at the settled state (the design surface\'s frame record)',
                lowest_sigma=sigma, inputs=hashes,
                fluxmod_namelist=namelist('fluxmod_namelist'), landmod_namelist=namelist('landmod_namelist'),
                planet_namelist=namelist('planet_namelist'),
                reading_rule='Arrays run over records (12 lunar days a model year, 10 means a lunar day, each 143 steps '
                             'of 30 minutes: 2.98 days), latitudes north first, longitudes east from 0. uas, vas the 10 m '
                             'wind (m/s) on Geleyn\'s profile over the land\'s 2 m roughness, tauu, tauv the surface stress '
                             '(Pa), ual, val the lowest level\'s wind (m/s), ps hPa, tas K, czen the mean cosine of the '
                             'Sun\'s zenith angle, pr m/s of water, hfss W/m2, mrso m; lsm the land mask.',
                evidence='ExoPlaSim 3.4.2 output means of the design climate\'s run at T21; nothing here is observed.')
    np.savez_compressed(SERIES, metadata=json.dumps(meta), lat=lat, lon=lon, lsm=lsm, **arrays)
    log(f'  wrote {SERIES} ({SERIES.stat().st_size / 1e6:.0f} MB)')
    return SERIES


def gas_constant(planet_namelist: str) -> float:
    for line in planet_namelist.splitlines():
        if line.strip().upper().startswith('GASCON'):
            return float(line.split('=')[1])
    raise ValueError('no GASCON in the planet namelist')


# ---- the drift ----------------------------------------------------------------------------------------------------
def weibull_nodes(shape, n=24):
    """n equally likely values of a Weibull variable of this shape, scaled to a mean of one."""
    p = (np.arange(n) + 0.5) / n
    x = (-np.log(1.0 - p)) ** (1.0 / shape)
    return x / x.mean()


def storm_speeds(gust=1.0, n=40, seed_nodes=12):
    """Equally likely 10 m winds over sand (m/s) within a storm's outflow: the storm's highest wind lognormal about
    the rings' median, the wind at a point uniform between half of it and all of it, carried to sand and to gusts."""
    from statistics import NormalDist
    z = np.array([NormalDist().inv_cdf((i + 0.5) / seed_nodes) for i in range(seed_nodes)])
    peak = STORM_GUST['median'] * np.exp(STORM_GUST['sigma_ln'] * z)
    part = 0.5 + 0.5 * (np.arange(n // seed_nodes + 1) + 0.5) / (n // seed_nodes + 1)
    return (peak[:, None] * part[None, :]).ravel() * RING_TO_SAND * gust


def drift(series, subgrid='central', to_sand=TO_SAND[1], z0=Z0_CENTRAL, ustar_ft=None, impact=IMPACT_RATIO[1],
          rho_p=GRAINS['feldspar']):
    """Each cell's yearly potential drift over bare loose sand, on the run's grid.

    Returns a dict of (lat, lon) arrays: dp_m2_yr the sand all winds together could carry (bulk m3 per metre of width a
    year), by_direction_m2_yr the same by the direction the sand travels (DIRECTIONS sectors clockwise from the turned
    Moon's north), resultant east and north parts, time_above the share of the time the wind exceeds the threshold,
    storm_m2_yr the storms' part (without direction), and the mean and percentiles of the 10 m wind over sand."""
    cfg = SUBGRID[subgrid]
    if ustar_ft is None:
        ustar_ft = easiest_grain(rho_p)[1]
    ustar_it = impact * ustar_ft
    u, v = series['uas'].astype(np.float64), series['vas'].astype(np.float64)
    speed = np.hypot(u, v) * to_sand * cfg['scale']
    sector = (np.floor((np.degrees(np.arctan2(u, v)) % 360.0) / (360.0 / DIRECTIONS) + 0.5).astype(int)) % DIRECTIONS
    nodes = np.ones(1) if cfg['spread'] is None else weibull_nodes(cfg['spread'])
    q = np.zeros(speed.shape)
    above = np.zeros(speed.shape)
    for x in nodes:
        us = friction_velocity(speed * x, z0)
        q += saltation_flux(us, ustar_it) / nodes.size
        above += (us > ustar_it) / nodes.size
    to_volume = SECONDS_PER_YEAR / BULK_DENSITY            # kg/m/s -> bulk m3 per metre a year
    records = speed.shape[0]
    by_dir = np.zeros((DIRECTIONS,) + speed.shape[1:])
    for k in range(DIRECTIONS):
        by_dir[k] = np.where(sector == k, q, 0.0).sum(axis=0) / records * to_volume
    ang = np.radians(np.arange(DIRECTIONS) * 360.0 / DIRECTIONS)
    east = (by_dir * np.sin(ang)[:, None, None]).sum(axis=0)
    north = (by_dir * np.cos(ang)[:, None, None]).sum(axis=0)
    rain = series['pr'].astype(np.float64).mean(axis=0) * 86400.0 * 1000.0          # mm/day
    storm_time = np.minimum(STORM_TIME * cfg['storms'] * rain / STORM_RAIN, STORM_TIME_CAP)
    storm_q = saltation_flux(friction_velocity(storm_speeds(cfg['storm_gust']), z0), ustar_it).mean() if cfg['storms'] else 0.0
    storm = storm_time * storm_q * to_volume
    return dict(dp_m2_yr=by_dir.sum(axis=0) + storm, by_direction_m2_yr=by_dir, east_m2_yr=east, north_m2_yr=north,
                storm_m2_yr=storm, time_above=above.mean(axis=0) + storm_time * (1.0 if cfg['storms'] else 0.0),
                speed_mean=speed.mean(axis=0), speed_p50=np.percentile(speed, 50, axis=0),
                speed_p90=np.percentile(speed, 90, axis=0), speed_p99=np.percentile(speed, 99, axis=0),
                ustar_it=ustar_it, ustar_ft=ustar_ft, rain_mm_day=rain)


def load_series():
    z = np.load(SERIES, allow_pickle=False)
    s = {k: z[k] for k in z.files if k != 'metadata'}
    s['metadata'] = json.loads(str(z['metadata']))
    return s


# ---- Parts C and D: supply and holding ----------------------------------------------------------------------------
COVER = dict(least=0.02, most=1.0, half=0.22, steep=2.2)   # the globe's cover curve (visualization/atlas/appearance.py COVER, a display rule), restated here
CLASS_EDGES = (0.05, 0.2, 0.5, 0.65, 1.0)                  # the moisture map's aridity classes (UNEP's, humid divided at 1)
CLASS_NAMES = ('hyper-arid', 'arid', 'semi-arid', 'dry sub-humid', 'humid', 'wet')
FORM_NAMES = ('none', 'sheet or streaks', 'crescent dunes and transverse ridges', 'linear dunes', 'star dunes',
              'dunes held under plants')
COVER_EFOLD = 0.05             # plant cover at which the sand flux falls by e: transport all but ends above 15% cover (Lancaster and Baas 1998, Owens Lake salt grass; search summary)
CRUST_FACTOR = 4.8             # a full biological crust's factor on the threshold (Rodriguez-Caballero et al. 2022, via the aerosol note)
CRUST_EDGE, CRUST_STEEP = 0.05, 4.0   # the kept crust's share of the open ground, 1 / (1 + (edge / aridity)^steep): half at the divide between the hyper-arid land, where the soils study keeps none (a crust grows only while wet), and the arid land, where it keeps one; a tenth at an aridity of 0.029, nine tenths at 0.087
FRESH_SHARE = 0.1              # share of crusted land's sand that lies fresh and loose at any time (floods, shores, disturbance); a bracket, 0.03-0.3
FRESH_TAKE = 0.3               # share of river-laid and shore sand the wind takes before a crust returns in 1-5 years; a bracket, 0.1-1
DRY_HOURS = 24.0               # hours a rained-on sand surface stays too damp to move (a bracket, 6-72 h)
STORM_RAIN_MM = 10.0           # rain a storm leaves at a point (Walnut Gulch's order, as the soils study reads the storm-water study)
LUNAR_DAY_H = 708.7
LUNAR_DAYS_PER_YEAR = 365.25 / 29.53059
EROSION_BANDS = ((0.0, 15.0), (15.0, 30.0), (30.0, 45.0), (45.0, 60.0), (60.0, 90.0))   # the terrain study's bands of latitude
EROSION_CENTRES_DEG = (7.5, 22.5, 37.5, 52.5, 75.0)
EROSION_M_PER_MYR = (79.0, 77.0, 57.0, 35.0, 25.0)   # its central lowering in each band, read here at the bands' middles and joined by straight lines between
EROSION_RANGE = (0.2, 2.0)     # its slowest and fastest cases against the central one (0.71 and 7.1 against 3.55 Gt a year)
EROSION_RAIN_POWER = 0.31      # the terrain study spreads its erosion as rain to this power (the BQART model's exponent on discharge); applied here to each mouth's rain against its band's mean
SAND_SIZE_SHARE = 0.4          # share of the eroded mantle between 63 um and 2 mm (Carrier 2003 through the soils study: median 40-130 um; river deltas: 41% bed material)
ROCK_DENSITY = 2600.0          # kg/m3, the terrain study's
BOWL_SHARE = 0.19              # share of the young land inside crater bowls of 10 m and more (storm-water study through the soils study)
BOWL_DEPTH_M = (0.3, 0.15)     # mean depth of those bowls on the young land and after ten thousand years of creep (terrain study: a 10 m bowl halves in 57-750 years, a 100 m bowl in 3,700-44,000); brackets
FULL_COVER_M = 1.5             # equivalent thickness at which free sand covers 63% of the ground (barchan fields to sand seas; Wasson and Hyde 1983, recalled)
DUNE_THICKNESS_M = 0.5         # a body of free sand thinner than this where the scheme sets it down lies as a sheet or streaks; a cell is drawn as dunes where most of its free sand, after the final sharing, came from deeper bodies
MOBILE_THICKNESS_M = (0.3, 30.0)   # thickness of the moving body: streaks to great dunes
TRAP_BASE, TRAP_GAIN = 0.05, 0.75  # share of a moving body a 1.9 km cell holds back at each step: open ground, and added in full shelter
SHELTER_REACH = 8              # cells looked upwind for higher ground (15 km)
SHELTER_SLOPE = (0.03, 0.10)   # rise over distance to the upwind skyline at which shelter begins and is complete (flow parts from the ground behind rims within about ten heights)
FLOOR_TPI_M, FLOOR_SLOPE = -150.0, 0.02   # a floor: this far below the mean ground within 10 km, and this gentle
MAX_STEPS = 200                # cell-to-cell moves a body of sand may make in the scheme (380 km along a meridian)
SIDE_SHARE = (0.02, 0.3)       # share of a moving body that steps to each side at a move, 0.5 (1 - net over total drift) within these limits: sand fans out where the winds vary (a proposal)
SOURCE_THICKNESS_M = 2.0       # loose sand a river end or a shore holds on a cell before it lies wider: the fan or beach it stands for is spread to this (a bracket, 1-5 m)
SETTLE_SIGMAS = (1.5, 3.0, 6.0, 12.0)   # cells: widths over which sand above that thickness is shared out, in turn
FINAL_SIGMA = 2.5              # cells (4.7 km): the width over which the carried sand and the bowls' sand are shared at the end, since the scheme places neither to a cell
CLIMATE_GRID = 4               # cells a degree of the turned-frame grid the climate's fields are carried on
CLIMATE_SMOOTHING_CELLS = 0.7  # of the run's grid spacing: the width the round smooths the climate's fields over (visualization/atlas/appearance.py, restated)
STEP_LIMIT, STEP_SHELTER = 0.03, 0.2   # the check on the map's smoothness: neighbouring cells differing by more than this in sand_share, away from water and from where shelter changes by more than this
CELL_M = 2.0 * math.pi * MOON_RADIUS / (360.0 * PPD)        # a cell's length along a meridian, 1,895 m
SUPPLY = {                     # what is offered; every entry a bracket's central value or its high end
    'central': dict(deflation_m=0.05,        # depth of raw mantle winnowed before a lag of fragments armours it (bracket 0.02-0.2 m)
                    durable=0.5,             # share of sand-sized grains that last as sand (mineral and rock fragments, beads) and do not break up as welded glass does (0.3-0.7)
                    erosion=1.0,             # factor on the terrain study's central lowering
                    catchment_cap_km2=5000.0,   # a mouth is credited with at most this much catchment: lakes upstream hold back the rest
                    shore_m2_yr=0.5,         # sand the waves and the monthly tide offer a downwind dry shore of a sea, bulk m3 per metre a year (0.1-3)
                    fetch_km=100.0),         # waters narrower than this offer in proportion to their width
    'ample': dict(deflation_m=0.2, durable=0.7, erosion=2.0, catchment_cap_km2=20000.0, shore_m2_yr=3.0, fetch_km=100.0),
}
VARIANTS = {
    'central':  dict(subgrid='central', open_sand=False, years=10000.0, supply='central',
                     rests_on='the run\'s winds carried to sand, the unresolved part at the wet 45 deg rings\' factor (1.3) with the storms\' tail; crusts kept on arid and semi-arid sand as on other open ground; central supply; the land ten thousand years after it dried'),
    'nogust':   dict(subgrid='none', open_sand=False, years=10000.0, supply='central',
                     rests_on='the central case on the run\'s 3-day means alone: no turbulence, no circulations below 170 km, no storms'),
    'highgust': dict(subgrid='high', open_sand=False, years=10000.0, supply='central',
                     rests_on='the central case with the unresolved part at the design-ground 45 deg rings\' factor (2.0), storms twice as frequent and their gusts half again as strong'),
    'opensand': dict(subgrid='central', open_sand=True, years=10000.0, supply='central',
                     rests_on='the central case with crust upkeep absent from sand (dune fields, sheets, bowl floors, river-laid and shore sand) in the arid and semi-arid land; plants alone hold it'),
    'young':    dict(subgrid='central', open_sand=False, years=1000.0, supply='central',
                     rests_on='the central case one thousand years after the land dried: a tenth of the rivers\' and shores\' sand, deeper bowls'),
    'ample':    dict(subgrid='central', open_sand=False, years=10000.0, supply='ample',
                     rests_on='the central case with every supply bracket at its high end: 0.2 m of mantle winnowed, seven tenths of the sand grains lasting, twice the erosion, mouths credited with 20,000 km2, shores offering 3 m3 a metre a year'),
}
RINGS = {   # the cloud-resolving cases whose 10 m winds over land part B rests on (climate/results/crm)
    'ring_ring_45n_lsw.json': '45 N, two-dimensional ring, land kept moist (evaporating 2.0 mm/day)',
    'ring_ring_45s_lsw.json': '45 S, two-dimensional ring, land kept moist',
    'ring_ring_45n_corrected.json': '45 N, two-dimensional ring on the design ground (wetness 0.014)',
    'ring_ring_45n_corrected_twilight.json': '45 N, the same with the curved air\'s light',
    'ring_ring_45s_corrected_twilight.json': '45 S, two-dimensional ring on the design ground (wetness 0.046)',
    'ring_ring_equator.json': 'equator, two-dimensional ring',
    'ring_ring_80n_lsw.json': '80 N, two-dimensional ring with the large-scale sinking',
    'ring_ring_80s_lsw.json': '80 S, two-dimensional ring with the large-scale sinking',
    'highland_box_highland.json': '44.7 S 246.1 E, three-dimensional box over the dry highlands, 100 km radius',
}
BANDS = ((0, 15), (15, 30), (30, 40), (40, 50), (50, 60), (60, 75), (75, 90))


def plant_cover(aridity):
    c = COVER
    return c['least'] + (c['most'] - c['least']) / (1.0 + (c['half'] / np.maximum(aridity, 1e-4)) ** c['steep'])


def class_of(aridity):
    return np.digitize(aridity, CLASS_EDGES).astype(np.int8)


def erosion_m_per_yr(lat_deg):
    """The terrain study's central lowering of the land (m of rock a year) at a latitude: its bands' rates at their
    middles, joined by straight lines."""
    return np.interp(np.abs(lat_deg), EROSION_CENTRES_DEG, EROSION_M_PER_MYR) * 1e-6


# ---- the grids ----------------------------------------------------------------------------------------------------
def fine_lat():
    return 90.0 - (np.arange(ROWS) + 0.5) / PPD


def fine_lon():
    return (np.arange(COLS) + 0.5) / PPD


def cell_area_km2():
    """Area of a cell in each row (km2)."""
    d = np.radians(1.0 / PPD)
    return (MOON_RADIUS / 1e3) ** 2 * d * d * np.cos(np.radians(fine_lat()))


def coarse_weights(lat_t, lon_t, lat_c, lon_c):
    """Bilinear weights of the run's grid at turned-frame places: rows and columns of the four cells and the shares
    along latitude and longitude. Beyond the last Gaussian latitude the edge row stands."""
    y = np.interp(-np.asarray(lat_t, float), -lat_c, np.arange(lat_c.size))
    y0 = np.clip(np.floor(y).astype(np.int64), 0, lat_c.size - 2)
    fy = np.clip(y - y0, 0.0, 1.0)
    x = (np.asarray(lon_t, float) % 360.0) / (360.0 / lon_c.size)
    x0 = np.floor(x).astype(np.int64) % lon_c.size
    return y0, x0, (x0 + 1) % lon_c.size, fy, x - np.floor(x)


def bilinear(field, w):
    y0, x0, x1, fy, fx = w
    return ((field[y0, x0] * (1 - fx) + field[y0, x1] * fx) * (1 - fy)
            + (field[y0 + 1, x0] * (1 - fx) + field[y0 + 1, x1] * fx) * fy)


def fill_from_land(field, known):
    """A field known on some of the run's cells spread onto the rest from the nearest known cells, cell by cell
    outward (longitude wraps)."""
    out, known = np.where(known, field, 0.0).astype(float), known.copy()
    while not known.all():
        total, count = np.zeros_like(out), np.zeros_like(out)
        for dr_, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            src = np.roll(out, (dr_, dc), axis=(0, 1))
            ok = np.roll(known, (dr_, dc), axis=(0, 1))
            if dr_ == -1:
                ok[-1] = False
            if dr_ == 1:
                ok[0] = False
            total += np.where(ok, src, 0.0)
            count += ok
        new = ~known & (count > 0)
        out[new] = total[new] / count[new]
        known |= new
    return out


def moisture_aridity(lat_c, lon_c):
    """The moisture map's aridity on the run's grid (its land cells; water cells from the nearest land), and the mask
    of the cells it classes."""
    c = json.loads(MOISTURE.read_text())['cells']
    a = np.zeros((lat_c.size, lon_c.size))
    have = np.zeros(a.shape, bool)
    a[np.asarray(c['row']), np.asarray(c['col'])] = np.asarray(c['aridity'], float)
    have[np.asarray(c['row']), np.asarray(c['col'])] = True
    return np.exp(fill_from_land(np.log(np.maximum(a, 1e-3)), have)), have


def read_surface():
    z = np.load(SURFACE, allow_pickle=False)
    meta = json.loads(str(z['metadata']))
    water = np.zeros(ROWS * COLS, bool)
    water[z['water_cell']] = True
    return z['height_m'], water.reshape(ROWS, COLS), meta


def step_offsets(az_deg, lat_deg, k):
    """Row and column offsets of the place k cell-lengths (1.9 km each) along an azimuth, the columns narrowing with
    latitude."""
    a = np.radians(az_deg)
    dr_ = -np.rint(k * np.cos(a)).astype(np.int32)
    dc = np.rint(k * np.sin(a) / np.maximum(np.cos(np.radians(lat_deg)), 0.02)).astype(np.int32)
    return dr_, dc


def grid_step(az_deg, lat_deg, step):
    """Row and column offsets of the step-th cell-to-cell move along an azimuth, and the metres a move covers. A move
    is one cell at most each way (columns narrow with latitude, so an eastward move is shorter near the poles), and
    successive moves alternate between neighbours so that their sum follows the direction."""
    a = np.radians(az_deg)
    vr, vc = -np.cos(a), np.sin(a) / np.maximum(np.cos(np.radians(lat_deg)), 0.02)
    norm = np.maximum(np.maximum(np.abs(vr), np.abs(vc)), 1e-9)
    ur, uc = vr / norm, vc / norm
    dr_ = (np.rint(step * ur) - np.rint((step - 1) * ur)).astype(np.int32)
    dc = (np.rint(step * uc) - np.rint((step - 1) * uc)).astype(np.int32)
    return dr_, dc, CELL_M / norm


def sphere_blur(a, sigma, lat_deg):
    """A Gaussian blur that is round on the ground, as the round's appearance.sphere_blur makes it: sigma rows north-south
    and, in each row, sigma / cos(latitude) columns east-west, wrapping round."""
    from scipy import fft
    from scipy.ndimage import gaussian_filter1d
    out = gaussian_filter1d(np.asarray(a, np.float32), sigma, axis=0, mode='nearest')
    across = sigma / np.maximum(np.cos(np.radians(lat_deg)), 1e-3)
    response = np.exp(-2.0 * np.pi ** 2 * (across[:, None] * np.fft.rfftfreq(out.shape[1])[None, :]) ** 2).astype(np.float32)
    return fft.irfft(fft.rfft(out, axis=1) * response, n=out.shape[1], axis=1).astype(np.float32)


def on_turned_grid(field, lat_c, lon_c):
    """A field of the run's grid on a grid of CLIMATE_GRID cells a degree on the turned Moon (rows from 90 N, columns
    from 0 E), smoothed over CLIMATE_SMOOTHING_CELLS of the run's spacing: the way the round's look carries the
    climate's fields (appearance.grid_field), so that no edge of a 170 km cell is left in them."""
    lat4 = 90.0 - (np.arange(180 * CLIMATE_GRID) + 0.5) / CLIMATE_GRID
    lon4 = (np.arange(360 * CLIMATE_GRID) + 0.5) / CLIMATE_GRID
    lon2, lat2 = np.meshgrid(lon4, lat4)
    f = bilinear(np.asarray(field, float), coarse_weights(lat2, lon2, lat_c, lon_c))
    return sphere_blur(f, CLIMATE_SMOOTHING_CELLS * float(np.abs(np.diff(lat_c)).mean()) * CLIMATE_GRID, lat4)


def turned_weights(lat_t, lon_t):
    """Bilinear weights of that turned-frame grid at turned-frame places."""
    rows, cols = 180 * CLIMATE_GRID, 360 * CLIMATE_GRID
    y = np.clip((90.0 - np.asarray(lat_t, float)) * CLIMATE_GRID - 0.5, 0.0, rows - 1.000001)
    y0 = np.floor(y).astype(np.int64)
    x = (np.asarray(lon_t, float) % 360.0) * CLIMATE_GRID - 0.5
    x0 = np.floor(x).astype(np.int64) % cols
    return y0, x0, (x0 + 1) % cols, y - y0, x - np.floor(x)


class Ground:
    """Sand shared among the land's cells by a blur that is round on the ground; none goes to water and none is lost."""

    def __init__(self, water):
        self.land = (~water).astype(np.float32)
        self.lat = fine_lat()
        self.area = cell_area_km2().astype(np.float32)[:, None]
        self.reach = {}

    def share(self, field, sigma):
        if sigma not in self.reach:
            self.reach[sigma] = np.maximum(sphere_blur(self.land, sigma, self.lat), 1e-3)
        out = np.maximum(self.land * sphere_blur(field / self.reach[sigma], sigma, self.lat), 0.0)
        had, has = float((field * self.area).sum()), float((out * self.area).sum())
        return (out * (had / has)).astype(np.float32) if has > 0 else out


def settle(free, ground, limit=SOURCE_THICKNESS_M):
    """Sand lying deeper than `limit` on a cell shared out over the land round it, wider at each turn, until none
    does: a river end lays its sand over its fan and a shore along its beach, and the cell stood for both."""
    for sigma in SETTLE_SIGMAS:
        excess = np.maximum(free - limit, 0.0)
        if float(excess.max()) < 0.05 * limit:
            break
        free = free - excess + ground.share(excess, sigma)
    return free.astype(np.float32)


def neighbour_steps(share, near_water, trap_index):
    """How smooth the map is: among pairs of neighbouring cells away from water and from the edges of shelter, the
    share whose sand_share differs by more than STEP_LIMIT, the same among pairs that hold sand, and the largest
    difference."""
    total = broken = sandy = sandy_broken = 0
    worst = 0.0
    for axis in (0, 1):
        pair = (lambda a: (a[:-1], a[1:])) if axis == 0 else (lambda a: (a, np.roll(a, -1, axis=1)))
        (a, b), (na, nb), (sa, sb) = pair(share), pair(near_water), pair(trap_index)
        ok = ~na & ~nb & (np.abs(sa - sb) <= STEP_SHELTER)
        d = np.abs(a - b)
        s = ok & ((a > 0) | (b > 0))
        total += int(ok.sum())
        broken += int((ok & (d > STEP_LIMIT)).sum())
        sandy += int(s.sum())
        sandy_broken += int((s & (d > STEP_LIMIT)).sum())
        worst = max(worst, float(d[ok].max()))
    return dict(pairs=total, share_of_pairs=round(broken / max(total, 1), 6), sandy_pairs=sandy,
                share_of_sandy_pairs=round(sandy_broken / max(sandy, 1), 5), largest_step=round(worst, 3))


def shelter(height, az_deg):
    """How far each cell lies in the lee of higher ground (0 open to 1 sheltered): the steepest rise to the ground
    within SHELTER_REACH cells upwind, against SHELTER_SLOPE."""
    out = np.zeros((ROWS, COLS), np.float32)
    lat = fine_lat()
    cols = np.arange(COLS, dtype=np.int32)[None, :]
    block = 360
    for r0 in range(0, ROWS, block):
        r1 = min(ROWS, r0 + block)
        az = az_deg[r0:r1]
        rows = np.arange(r0, r1, dtype=np.int32)[:, None]
        la = lat[r0:r1, None]
        best = np.full(az.shape, -1.0, np.float32)
        for k in range(1, SHELTER_REACH + 1):
            dr_, dc = step_offsets((az + 180.0) % 360.0, la, k)
            rise = (height[np.clip(rows + dr_, 0, ROWS - 1), (cols + dc) % COLS] - height[r0:r1]) / (k * CELL_M)
            best = np.maximum(best, rise.astype(np.float32))
        out[r0:r1] = np.clip((best - SHELTER_SLOPE[0]) / (SHELTER_SLOPE[1] - SHELTER_SLOPE[0]), 0.0, 1.0)
    return out


def floors(height):
    """1 on gentle ground well below its surroundings (crater and basin floors), else 0."""
    from scipy import ndimage
    mean = ndimage.uniform_filter(height, size=11, mode='nearest')
    gy, gx = np.gradient(height, CELL_M)
    gx = gx / np.maximum(np.cos(np.radians(fine_lat()))[:, None], 0.02).astype(np.float32)
    return ((height - mean < FLOOR_TPI_M) & (np.hypot(gx, gy) < FLOOR_SLOPE)).astype(np.float32)


def carry(free, q_eff, az_deg, trap, water, bowl, bowl_cap, years, side=None, log=print):
    """Carry the free sand (equivalent thickness, m) downwind cell by cell, in place. A body moves on while the drift
    could take it this far in `years` (a thin body travels faster than a thick one); each cell holds back TRAP
    of it at every step and fills its bowls first; sand that reaches water is lost; a cell keeps at most the thickest
    dunes' sand and passes the rest on. With `side`, that share of a moving body steps to each side at every move, so a
    body fans out where the winds vary. Returns the volume lost (km3)."""
    lat = fine_lat()
    area = cell_area_km2() * 1e6                                  # m2 by row
    flat_free, flat_bowl = free.reshape(-1), bowl.reshape(-1)
    flat_q, flat_az, flat_trap, flat_water = q_eff.reshape(-1), az_deg.reshape(-1), trap.reshape(-1), water.reshape(-1)
    flat_side = None if side is None else side.reshape(-1)
    lost = 0.0
    idx = np.flatnonzero(flat_free > 1e-4)
    for step in range(1, MAX_STEPS + 1):
        if idx.size == 0:
            break
        v = flat_free[idx].astype(np.float64)
        rows = idx // COLS
        dr_, dc, stride = grid_step(flat_az[idx], lat[rows], step)
        thick = np.clip(2.0 * v, *MOBILE_THICKNESS_M)
        r2 = rows + dr_
        go = (v > 1e-4) & (flat_q[idx] * years / (thick * stride) >= step) & np.isfinite(flat_az[idx]) & (r2 >= 0) & (r2 < ROWS)
        idx, v, rows, r2, dr_, dc = idx[go], v[go], rows[go], r2[go], dr_[go], dc[go]
        if idx.size == 0:
            break
        cols = idx % COLS
        moved = v * (1.0 - flat_trap[idx])                       # thickness leaving the cell
        volume = moved * area[rows]
        flat_free[idx] -= moved.astype(np.float32)
        t_rows, t_cols, t_vol = r2, (cols + dc) % COLS, volume
        if flat_side is not None:                                # a quarter turn to the left and to the right of the move
            s = flat_side[idx].astype(np.float64)
            t_rows = np.concatenate([r2, np.clip(rows - dc, 0, ROWS - 1), np.clip(rows + dc, 0, ROWS - 1)])
            t_cols = np.concatenate([t_cols, (cols + dr_) % COLS, (cols - dr_) % COLS])
            t_vol = np.concatenate([volume * (1.0 - 2.0 * s), volume * s, volume * s])
        target = t_rows * COLS + t_cols
        wet = flat_water[target]
        lost += float(t_vol[wet].sum())
        got, inverse = np.unique(target[~wet], return_inverse=True)
        arrive = np.bincount(inverse, weights=t_vol[~wet] / area[t_rows[~wet]], minlength=got.size)
        into = np.minimum(arrive, np.maximum(bowl_cap[got] - flat_bowl[got], 0.0))
        flat_bowl[got] += into.astype(np.float32)
        flat_free[got] += (arrive - into).astype(np.float32)
        idx = np.union1d(idx, got)                               # what stayed behind and what arrived may both move on
        if step % 40 == 0:
            log(f'      step {step}: {idx.size} cells hold moving sand')
    # no cell holds more than the thickest dunes' sand: what is over that spills on downwind, as a river's end or a
    # shore spreads its sand over more ground than one cell
    cap = MOBILE_THICKNESS_M[1]
    for turn in range(1, 401):
        over = np.flatnonzero(flat_free > cap)
        if over.size == 0:
            break
        rows, cols = over // COLS, over % COLS
        az = np.where(np.isfinite(flat_az[over]), flat_az[over], 0.0)
        dr_, dc, _ = grid_step(az, lat[rows], turn)
        r2 = np.clip(rows + dr_, 0, ROWS - 1)
        target = r2 * COLS + (cols + dc) % COLS
        volume = (flat_free[over].astype(np.float64) - cap) * area[rows]
        flat_free[over] = cap
        wet = flat_water[target]
        lost += float(volume[wet].sum())
        got, inverse = np.unique(target[~wet], return_inverse=True)
        flat_free[got] += np.bincount(inverse, weights=volume[~wet] / area[r2[~wet]], minlength=got.size).astype(np.float32)
    return lost / 1e9


class NpzWriter:
    """An .npz written one array at a time, so the product never stands whole in memory."""

    def __init__(self, path):
        import zipfile
        self.zf = zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED, allowZip64=True)

    def add(self, name, array):
        array = np.asarray(array)
        with self.zf.open(name + '.npy', 'w', force_zip64=True) as f:
            np.lib.format.write_array(f, np.ascontiguousarray(array) if array.ndim else array, allow_pickle=False)

    def close(self):
        self.zf.close()


def build(store, log=print):
    """Everything from the reduced winds to the map: arrays go to store(name, array) as they are finished; the sums
    come back."""
    t0 = time.time()
    s = load_series()
    lat_c, lon_c = s['lat'].astype(float), s['lon'].astype(float)
    land_c = s['lsm'] > 0.5
    height, water, surface_meta = read_surface()
    rotation = np.asarray(surface_meta['frame']['rotation_today_to_turned'], float)
    aridity_c, classed = moisture_aridity(lat_c, lon_c)
    to_sand = np.where(land_c, TO_SAND[1], 1.0)                 # over the run's water the 10 m wind already blows over a smooth surface
    lonc2, latc2 = np.meshgrid(lon_c, lat_c)
    drifts = {}
    for sub in SUBGRID:
        d = drift(s, sub, to_sand=to_sand[None])
        d['vector_today'] = vector_to_today(d['east_m2_yr'], d['north_m2_yr'], latc2, lonc2, rotation)
        drifts[sub] = d
    log(f'  coarse drift {time.time() - t0:.0f} s')

    used = sorted({v['subgrid'] for v in VARIANTS.values()})
    lat_f, lon_f = fine_lat(), fine_lon()
    area_row = cell_area_km2()
    aridity = np.empty((ROWS, COLS), np.float32)
    rain = np.empty((ROWS, COLS), np.float32)
    mafic = np.empty((ROWS, COLS), np.float32)
    rock = np.load(ROCK, allow_pickle=False)['mafic']
    fine = {sub: {k: np.empty((ROWS, COLS), np.float32) for k in ('dp', 'rdp', 'az', 'above')} for sub in used}
    smooth = lambda field: on_turned_grid(field, lat_c, lon_c)
    ln_aridity = smooth(np.log(aridity_c))
    ln_rain = smooth(np.log(np.maximum(drifts['central']['rain_mm_day'] * 365.0, 0.1)))
    carried = {sub: dict(vec=[smooth(drifts[sub]['vector_today'][..., i]) for i in range(3)], dp=smooth(drifts[sub]['dp_m2_yr']),
                         above=smooth(drifts[sub]['time_above'])) for sub in used}
    block = 180
    for r0 in range(0, ROWS, block):
        r1 = min(ROWS, r0 + block)
        lon2, lat2 = np.meshgrid(lon_f, lat_f[r0:r1])
        lat_t, lon_t = to_turned(lat2, lon2, rotation)
        w = turned_weights(lat_t, lon_t)
        aridity[r0:r1] = np.exp(bilinear(ln_aridity, w))
        rain[r0:r1] = np.exp(bilinear(ln_rain, w))
        rr = np.clip(np.floor((90.0 - lat_t) * PPD).astype(np.int64), 0, ROWS - 1)
        cc = np.floor((lon_t + 180.0) * PPD).astype(np.int64) % COLS
        mafic[r0:r1] = rock[rr, cc]                                # the rock's rows and columns are the turned Moon's, 180 W first
        for sub in used:
            d = carried[sub]
            vec = np.stack([bilinear(d['vec'][i], w) for i in range(3)], axis=-1)
            e, n, az = azimuth_today(vec, lat2, lon2)
            f = fine[sub]
            f['rdp'][r0:r1] = np.hypot(e, n)
            f['az'][r0:r1] = az
            f['dp'][r0:r1] = np.maximum(bilinear(d['dp'], w), 0.0)
            f['above'][r0:r1] = np.clip(bilinear(d['above'], w), 0.0, 1.0)
    del rock
    log(f'  fine fields {time.time() - t0:.0f} s')
    cover = plant_cover(aridity).astype(np.float32)
    classes = class_of(aridity)
    through_plants = np.exp(-cover / COVER_EFOLD).astype(np.float32)
    damp = np.clip(rain / LUNAR_DAYS_PER_YEAR / STORM_RAIN_MM * DRY_HOURS / LUNAR_DAY_H, 0.0, 0.9).astype(np.float32)
    moving = through_plants * (1.0 - damp)
    kept_crust = (1.0 / (1.0 + (CRUST_EDGE / np.maximum(aridity, 1e-4)) ** CRUST_STEEP)).astype(np.float32)   # the kept crust's share of the open ground
    del damp
    lee = shelter(height, fine['central']['az'])
    trap_index = np.maximum(lee, floors(height))
    trap = (TRAP_BASE + TRAP_GAIN * trap_index).astype(np.float32)
    exposure = ((1.0 - lee) ** 2).astype(np.float32)
    del lee
    log(f'  shelter and floors {time.time() - t0:.0f} s')

    routing = np.load(ROUTING, allow_pickle=False)
    flat_water = water.reshape(-1)
    mouth = (~flat_water) & flat_water[np.clip(routing['receiver'], 0, ROWS * COLS - 1)]
    mouth |= routing['keeps'] & ~flat_water
    drained = np.where(mouth, routing['drained_km2'], 0.0).astype(np.float32).reshape(ROWS, COLS)   # catchment behind each mouth, km2
    mouths_n = int(mouth.sum())
    del mouth
    from scipy import ndimage
    labels, count = ndimage.label(water)
    sizes = np.bincount(labels.reshape(-1), minlength=count + 1).astype(np.float64)
    rows2 = np.arange(ROWS, dtype=np.int32)[:, None]
    cols2 = np.arange(COLS, dtype=np.int32)[None, :]
    dr_, dc, _ = grid_step((fine['central']['az'] + 180.0) % 360.0, lat_f[:, None], 1)
    up_r, up_c = np.clip(rows2 + dr_, 0, ROWS - 1), (cols2 + dc) % COLS
    upwind_water_km = np.where(~water & water[up_r, up_c], np.sqrt(sizes * 3.59)[labels[up_r, up_c]], 0.0).astype(np.float32)   # width of the water a downwind shore faces
    shore_n = int((upwind_water_km > 0).sum())
    del labels, up_r, up_c, dr_, dc
    yield_row = erosion_m_per_yr(lat_f) * SAND_SIZE_SHARE * ROCK_DENSITY / BULK_DENSITY    # m of bulk sand-sized grains a year off the catchment
    means = []                                                # the land's mean rain in each of the terrain study's bands
    for lower, upper in EROSION_BANDS:
        band = (np.abs(lat_f) >= lower) & (np.abs(lat_f) < upper + 1e-9)
        weight = area_row[band][:, None] * ~water[band]
        means.append(float((rain[band] * weight).sum() / weight.sum()))
    band_rain = dict(centres_deg=list(EROSION_CENTRES_DEG), mean_mm_yr=[round(m, 1) for m in means])
    typical = np.interp(np.abs(lat_f), EROSION_CENTRES_DEG, means)
    rain_rule = ((np.maximum(rain, 0.1) / typical[:, None]) ** EROSION_RAIN_POWER).astype(np.float32)   # each mouth's rain against its latitude's, to the terrain study's power
    del rain
    ground = Ground(water)
    near_water = ndimage.binary_dilation(water, structure=np.ones((3, 3), bool))
    log(f'  sources {time.time() - t0:.0f} s: {mouths_n} mouths, {shore_n} downwind shore cells')

    region = np.repeat(np.repeat(np.load(REGIONS, allow_pickle=False)['region'], 2, axis=0), 2, axis=1).reshape(-1).astype(np.int64)
    region_list = json.loads(REGION_RECORD.read_text())['regions']
    a2 = np.broadcast_to(area_row[:, None].astype(np.float32), (ROWS, COLS))
    land_by_region = np.bincount(region, weights=(a2 * ~water).reshape(-1), minlength=200)
    land_km2 = float((a2 * ~water).sum())
    class_km2 = {CLASS_NAMES[k]: round(float(a2[(classes == k) & ~water].sum())) for k in range(6)}
    crest_min = ELEMENTARY_FACTOR[0] * saturation_length(easiest_grain(GRAINS['feldspar'])[0], GRAINS['feldspar'])
    summary = {}
    for name, var in VARIANTS.items():
        f, years, sup = fine[var['subgrid']], var['years'], SUPPLY[var['supply']]
        kept = np.zeros_like(kept_crust) if var['open_sand'] else kept_crust
        crust = (1.0 - kept * (1.0 - FRESH_SHARE)).astype(np.float32)
        take = (1.0 - kept * (1.0 - FRESH_TAKE)).astype(np.float32)
        windy = np.clip(f['above'] / 0.01, 0.0, 1.0)
        q_eff = (f['rdp'] * moving * crust * exposure).astype(np.float32)
        bowl_depth = BOWL_DEPTH_M[0] if years <= 1000.0 else BOWL_DEPTH_M[1]
        bowl_cap = np.where(water, 0.0, BOWL_SHARE * bowl_depth).astype(np.float32)
        loose = through_plants * crust                              # the share of in-place sand that stays loose sand and does not become soil under cover
        mantle = (sup['deflation_m'] * SAND_SIZE_SHARE * sup['durable'] * (1.0 - cover) * loose).astype(np.float32)
        mantle[water] = 0.0
        river = (yield_row[:, None] * rain_rule * sup['erosion'] * sup['durable'] * np.minimum(drained, sup['catchment_cap_km2']) / area_row[:, None]
                 * years * through_plants * take * windy).astype(np.float32)
        shore = (sup['shore_m2_yr'] * np.clip(upwind_water_km / sup['fetch_km'], 0.0, 1.0) * CELL_M / (area_row[:, None] * 1e6)
                 * years * through_plants * take * windy).astype(np.float32)
        supply_km3 = {k: round(float((v * a2).sum() / 1e3), 1) for k, v in (('mantle', mantle), ('rivers', river), ('shores', shore))}
        offered = mantle + settle(river + shore, ground)            # a river end's or a shore's sand lies over its fan or beach
        del mantle, river, shore, take, loose
        bowl = np.minimum(offered, bowl_cap)
        free = offered - bowl
        del offered
        log(f'  {name}: offered {supply_km3} km3')
        ratio = f['rdp'] / np.maximum(f['dp'], 1e-9)
        side = np.clip(0.5 * (1.0 - ratio), *SIDE_SHARE).astype(np.float32)
        lost = carry(free, q_eff, f['az'], trap, water, bowl, bowl_cap.reshape(-1), years, side, log)
        del side
        free = settle(free, ground)
        deep = free >= DUNE_THICKNESS_M                              # deep enough for dunes where it is set down
        piled = ground.share(np.where(deep, free, 0.0).astype(np.float32), FINAL_SIGMA)
        seeded = ground.share(deep.astype(np.float32), FINAL_SIGMA)
        free = ground.share(free, FINAL_SIGMA)
        bowl = np.minimum(ground.share(bowl, FINAL_SIGMA), bowl_cap)
        del deep
        share = np.clip(BOWL_SHARE * np.sqrt(bowl / np.maximum(bowl_cap, 1e-6)) + (1.0 - np.exp(-free / FULL_COVER_M)), 0.0, 1.0).astype(np.float32)
        share[water] = 0.0
        share[share < 0.005] = 0.0
        has = share > 0
        active = (moving * crust * windy).astype(np.float32)
        active[~has] = 0.0
        form = np.zeros((ROWS, COLS), np.uint8)
        dunes = has & (piled > 0.5 * free) & (free > 0.02)
        form[has] = 1
        form[dunes & (ratio >= 0.8)] = 2
        form[dunes & (ratio >= 0.3) & (ratio < 0.8)] = 3
        form[dunes & (ratio < 0.3)] = 4
        form[dunes & (active < 0.2)] = 5
        az = f['az']
        crest = np.full((ROWS, COLS), np.nan, np.float32)
        crest[form == 2] = (az[form == 2] + 90.0) % 180.0
        for k in (3, 5):
            crest[form == k] = az[form == k] % 180.0
        spacing = np.zeros((ROWS, COLS), np.float32)
        isdune = form >= 2
        body = piled[isdune] / np.maximum(seeded[isdune], 1e-3)      # thickness of the dune sand where it was set down
        grown = np.minimum(np.sqrt(2.0 * CREST_ASPECT * np.maximum(f['rdp'][isdune] * moving[isdune], 1e-6) * years),
                           3.0 * CREST_ASPECT * body)
        spacing[isdune] = np.clip(grown, crest_min, MIXED_LAYER_M[0])
        drift_az = np.where(has & (active > 0), az, np.nan).astype(np.float32)
        suffix = '' if name == 'central' else f'__{name}'
        thickness = (free + bowl).astype(np.float32)
        for key, array in (('sand_share', share), ('active_share', active), ('form', form), ('crest_azimuth_deg', crest),
                           ('crest_spacing_m', spacing), ('drift_azimuth_deg', drift_az), ('sand_thickness_m', thickness)):
            store(key + suffix, array)
        if name in ('central', 'nogust', 'highgust'):
            store('drift_potential_m2_yr' + suffix, np.where(water, 0.0, np.round(f['dp'], 1)).astype(np.float32))
        if name == 'central':
            store('drift_resultant_m2_yr', np.where(water, 0.0, np.round(f['rdp'], 1)).astype(np.float32))
            store('drift_direction_deg', np.where(water, np.nan, np.round(az)).astype(np.float32))
        sand_km2 = (share * a2).reshape(-1)
        moving_km2 = sand_km2 * active.reshape(-1)
        flat_form = form.reshape(-1).astype(np.int64)
        by_form = np.bincount(flat_form, weights=sand_km2, minlength=6)
        by_class = np.bincount(classes.reshape(-1).astype(np.int64), weights=sand_km2, minlength=6)
        r_sand = np.bincount(region, weights=sand_km2, minlength=200)
        r_move = np.bincount(region, weights=moving_km2, minlength=200)
        r_form = np.bincount(region * 6 + flat_form, weights=sand_km2, minlength=1200).reshape(200, 6)
        r_basalt = np.bincount(region, weights=sand_km2 * mafic.reshape(-1), minlength=200)
        table = []
        for i, r in enumerate(region_list):
            k = i + 1
            if r_sand[k] < 1.0 or land_by_region[k] <= 0:
                continue
            table.append(dict(code=r['code'], at_today=r['at']['today'], land_km2=round(float(land_by_region[k])),
                              sand_km2=round(float(r_sand[k])), sand_share=round(float(r_sand[k] / land_by_region[k]), 4),
                              moving_km2=round(float(r_move[k])), dunes_km2=round(float(r_form[k, 2:].sum())),
                              basalt_share=round(float(r_basalt[k] / r_sand[k]), 2),
                              main_form=FORM_NAMES[int(r_form[k, 1:].argmax()) + 1]))
        table.sort(key=lambda t: -t['sand_km2'])
        total = float(sand_km2.sum())
        summary[name] = dict(rests_on=var['rests_on'], years=years, subgrid=var['subgrid'], open_sand=var['open_sand'],
                             supply=var['supply'], sand_km2=round(total), share_of_land=round(total / land_km2, 5),
                             moving_km2=round(float(moving_km2.sum())), moving_share=round(float(moving_km2.sum() / max(total, 1e-9)), 4),
                             dunes_km2=round(float(by_form[2:].sum())),
                             basalt_share_of_sand=round(float((sand_km2 * mafic.reshape(-1)).sum() / max(total, 1e-9)), 3),
                             cells_over_half_sand_km2=round(float(a2[share > 0.5].sum())),
                             by_form_km2={FORM_NAMES[k]: round(float(by_form[k])) for k in range(1, 6)},
                             by_class_km2={CLASS_NAMES[k]: round(float(by_class[k])) for k in range(6)},
                             volume_km3=round(float((thickness * a2).sum() / 1e3), 1), in_bowls_km3=round(float((bowl * a2).sum() / 1e3), 1),
                             free_km3=round(float((free * a2).sum() / 1e3), 1), offered_km3=supply_km3, lost_to_water_km3=round(lost, 1),
                             crest_spacing_m=dict(median=round(float(np.median(spacing[isdune]))) if isdune.any() else None,
                                                  p90=round(float(np.percentile(spacing[isdune], 90))) if isdune.any() else None),
                             thickest_m=round(float(free.max()), 1), neighbour_steps=neighbour_steps(share, near_water, trap_index), regions=table)
        log(f'  {name}: sand {summary[name]["sand_km2"]} km2, moving {summary[name]["moving_km2"]}, dunes {summary[name]["dunes_km2"]}, lost {lost:.1f} km3, steps {summary[name]["neighbour_steps"]["share_of_sandy_pairs"]} ({time.time() - t0:.0f} s)')
        del share, active, form, crest, spacing, drift_az, thickness, free, bowl, q_eff, sand_km2, moving_km2, flat_form, ratio, piled, seeded
    store('basalt_share', np.where(water, 0.0, np.round(mafic, 2)).astype(np.float32))
    store('plant_cover', np.where(water, 0.0, np.round(cover, 2)).astype(np.float32))
    store('shelter', np.round(trap_index, 2).astype(np.float32))
    return dict(summary=summary, drifts=drifts, series=s, surface_meta=surface_meta, rotation=rotation, land_c=land_c,
                aridity_c=aridity_c, classed=classed, land_km2=land_km2, class_km2=class_km2, regions=region_list,
                counts=dict(mouths=mouths_n, downwind_shore_cells=shore_n), band_rain_mm_yr=band_rain, seconds=time.time() - t0)


# ---- the record ---------------------------------------------------------------------------------------------------
def physics():
    """Part A's figures."""
    out = dict(gravity_m_s2=G_MOON, air_density_kg_m3=AIR_DENSITY, earth=dict(gravity_m_s2=G_EARTH, air_density_kg_m3=AIR_DENSITY_EARTH))
    sizes_um = (20, 40, 63, 100, 130, 180, 264, 500, 1000, 2000)
    grains = {}
    for name, rho_p in GRAINS.items():
        d, u = easiest_grain(rho_p)
        lo, hi = easiest_grain(rho_p, gamma=SHAO_LU_GAMMA[0]), easiest_grain(rho_p, gamma=SHAO_LU_GAMMA[2])
        grains[name] = dict(
            density_kg_m3=rho_p, easiest_grain_um=round(d * 1e6), easiest_grain_um_range=[round(lo[0] * 1e6), round(hi[0] * 1e6)],
            fluid_threshold_m_s=round(u, 4), fluid_threshold_range_m_s=[round(lo[1], 4), round(hi[1], 4)],
            fluid_wind_10m_m_s={f'z0_{z:g}_m': round(float(wind_at_10m(u, z)), 2) for z in SAND_Z0},
            impact_wind_10m_m_s={f'ratio_{r:g}': [round(float(wind_at_10m(r * u, z)), 2) for z in SAND_Z0] for r in IMPACT_RATIO},
            threshold_by_size_m_s={str(sz): round(float(fluid_threshold(sz * 1e-6, rho_p)), 4) for sz in sizes_um},
            saturation_length_m={str(sz): round(float(saturation_length(sz * 1e-6, rho_p)), 2) for sz in (130, 264, 500)},
            smallest_dunes_m={str(sz): [round(float(k * saturation_length(sz * 1e-6, rho_p))) for k in ELEMENTARY_FACTOR] for sz in (130, 264, 500)})
    out['grains'] = grains
    dq, uq = easiest_grain(2650.0, AIR_DENSITY_EARTH, G_EARTH)
    out['earth']['quartz'] = dict(easiest_grain_um=round(dq * 1e6), fluid_threshold_m_s=round(uq, 4),
                                  fluid_wind_10m_m_s={f'z0_{z:g}_m': round(float(wind_at_10m(uq, z)), 2) for z in SAND_Z0},
                                  saturation_length_m_180um=round(float(saturation_length(180e-6, 2650.0, AIR_DENSITY_EARTH)), 2))
    uf = grains['feldspar']['fluid_threshold_m_s']
    out['threshold_ratio_moon_over_earth'] = round(uf / easiest_grain(GRAINS['feldspar'], AIR_DENSITY_EARTH, G_EARTH)[1], 3)
    flux = {}
    for u10 in (2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0):
        us = float(friction_velocity(u10))
        q = float(saltation_flux(us, IMPACT_RATIO[1] * uf))
        qe = float(saltation_flux(us, IMPACT_RATIO[1] * uq, AIR_DENSITY_EARTH, G_EARTH))
        flux[str(u10)] = dict(ustar_m_s=round(us, 3), moon_kg_m_s=round(q, 5), moon_m3_per_m_day=round(q * 86400.0 / BULK_DENSITY, 3),
                              earth_kg_m_s=round(qe, 5), earth_m3_per_m_day=round(qe * 86400.0 / BULK_DENSITY, 3))
    out['flux_by_10m_wind_over_sand'] = flux
    out['flux_ratio_at_equal_ustar_and_threshold'] = round(AIR_DENSITY / AIR_DENSITY_EARTH * G_EARTH / G_MOON, 2)
    out['hop_without_drag_m'] = {str(us): [dict(height=round(h, 3), length=round(l, 2)) for h, l in hop(us)] for us in (0.15, 0.25)}
    out['hop_without_drag_earth_m'] = {str(us): [dict(height=round(h, 3), length=round(l, 2)) for h, l in hop(us, G_EARTH)] for us in (0.25, 0.4)}
    out['wind_to_sand'] = dict(model_land_roughness_m=MODEL_LAND_Z0, lowest_level_m=LOWEST_LEVEL_M, factor_low_central_high=list(TO_SAND),
                               both_neutral=round(sand_wind_factor(), 3), sand_at_0p1_mm=round(sand_wind_factor(SAND_Z0[0]), 3),
                               model_share_0p52=round(sand_wind_factor(model_share=0.52), 3),
                               ring_to_sand=dict(used=RING_TO_SAND, level_50_m=round(ring_to_sand(), 3), level_313_m=round(ring_to_sand(level=313.0), 3)))
    out['largest_dunes_m'] = dict(mixed_layer_wet_rings=MIXED_LAYER_M[0], mixed_layer_design_ground_ring=MIXED_LAYER_M[1])
    return out


def weighted_quantiles(values, weights, qs):
    order = np.argsort(values)
    c = np.cumsum(weights[order])
    return [float(values[order][min(int(np.searchsorted(c, c[-1] * q / 100.0)), values.size - 1)]) for q in qs]


def wind_tables(b):
    """Part B's figures: the run's own winds, and for each gust case the wind over sand and the drift by band of
    latitude and by moisture class, weighted by area."""
    s, drifts, land = b['series'], b['drifts'], b['land_c']
    lat, lon = s['lat'].astype(float), s['lon'].astype(float)
    w = np.cos(np.radians(lat))[:, None] * np.ones((1, lon.size))
    sp = np.hypot(s['uas'], s['vas']).astype(np.float64)
    spl = np.hypot(s['ual'], s['val']).astype(np.float64)
    rho = s['ps'] * 100.0 / (gas_constant(s['metadata']['planet_namelist']) * s['tas'])
    ust = np.sqrt(np.hypot(s['tauu'], s['tauv']) / rho)
    ww = np.broadcast_to(w, sp.shape)
    groups = {f'{lo}-{hi}': land & ((np.abs(lat) >= lo) & (np.abs(lat) < hi))[:, None] for lo, hi in BANDS}
    cls = class_of(b['aridity_c'])
    groups.update({CLASS_NAMES[k]: land & b['classed'] & (cls == k) for k in range(6)})
    groups['all land'] = land
    run = {}
    for name, m in groups.items():
        if not m.any():
            continue
        q10 = weighted_quantiles(sp[:, m].ravel(), ww[:, m].ravel(), (50, 90, 99))
        ql = weighted_quantiles(spl[:, m].ravel(), ww[:, m].ravel(), (50, 90, 99))
        run[name] = dict(cells=int(m.sum()), wind_10m_p50_p90_p99_m_s=[round(x, 2) for x in q10], wind_10m_max_m_s=round(float(sp[:, m].max()), 2),
                         lowest_level_p50_p90_p99_m_s=[round(x, 2) for x in ql],
                         ustar_p50_p90_m_s=[round(x, 3) for x in weighted_quantiles(ust[:, m].ravel(), ww[:, m].ravel(), (50, 90))])
    share = sp[:, land] / np.maximum(spl[:, land], 1e-6)
    out = dict(run_over_its_own_roughness=run, run_10m_over_lowest_level_median=round(float(np.median(share)), 3),
               air_density_land_mean_kg_m3=round(float((rho.mean(axis=0) * w)[land].sum() / w[land].sum()), 3))
    eq = int(np.argmin(np.abs(lat)))
    lon_sun = np.degrees(np.angle((s['czen'][:, eq, :] * np.exp(1j * np.radians(lon))[None, :]).sum(axis=1))) % 360.0
    ha = np.broadcast_to((lon[None, None, :] - lon_sun[:, None, None] + 180.0) % 360.0 - 180.0, sp.shape)
    m = groups['40-50']
    edges = np.arange(-180, 181, 30)
    out['run_10m_wind_by_hour_angle_40_50'] = dict(
        hour_angle_deg=[float(x) for x in 0.5 * (edges[:-1] + edges[1:])],
        wind_10m_m_s=[round(float(sp[:, m][(ha[:, m] >= a) & (ha[:, m] < c)].mean()), 2) for a, c in zip(edges[:-1], edges[1:])],
        lowest_level_m_s=[round(float(spl[:, m][(ha[:, m] >= a) & (ha[:, m] < c)].mean()), 2) for a, c in zip(edges[:-1], edges[1:])],
        note='hour angle 0 at noon, positive in the afternoon; the Sun\'s longitude from the mean cosine of its zenith angle along the equator')
    cases = {}
    for sub, d in drifts.items():
        rdp = np.hypot(d['east_m2_yr'], d['north_m2_yr'])
        rows = {}
        for name, m in groups.items():
            if not m.any():
                continue
            x = w[m]
            mean = lambda a: float((a[m] * x).sum() / x.sum())
            rows[name] = dict(wind_over_sand_p50_p90_p99_m_s=[round(mean(d[k]), 2) for k in ('speed_p50', 'speed_p90', 'speed_p99')],
                              time_above_threshold=round(mean(d['time_above']), 4),
                              drift_potential_m2_yr=dict(mean=round(mean(d['dp_m2_yr']), 2), median=round(float(np.median(d['dp_m2_yr'][m])), 2),
                                                         highest=round(float(d['dp_m2_yr'][m].max()), 1)),
                              storms_part_m2_yr=round(mean(d['storm_m2_yr']), 2),
                              resultant_over_potential=round(float((rdp[m] * x).sum() / max(float((d['dp_m2_yr'][m] * x).sum()), 1e-12)), 2),
                              mobility_index=round(100.0 * mean(d['time_above']) / max(float((b['aridity_c'][m] * x).sum() / x.sum()), 1e-6), 0))
        cases[sub] = dict(settings=SUBGRID[sub], impact_threshold_m_s=round(float(d['ustar_it']), 4), by_group=rows)
    out['over_sand'] = cases
    return out


def pathways(b, top=12):
    """Part E at the run's scale: where the central case's resultant drift converges (what a saturated flux would lay
    down) and where it heads, by band."""
    d = b['drifts']['central']
    s = b['series']
    lat, lon = s['lat'].astype(float), s['lon'].astype(float)
    land = b['land_c']
    u, v = d['east_m2_yr'], d['north_m2_yr']
    phi = np.radians(lat)[:, None]
    dlon = np.radians(lon[1] - lon[0])
    du = (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2.0 * dlon)
    vc = v * np.cos(phi)
    dv = np.gradient(vc, np.radians(lat), axis=0)
    conv = -(du + dv) / (MOON_RADIUS * np.cos(phi))             # m of bulk sand a year, positive where the drift converges
    conv_land = np.where(land, conv, -np.inf)
    order = np.argsort(conv_land.ravel())[::-1][:top]
    lonc2, latc2 = np.meshgrid(lon, lat)
    lat_today, lon_today = to_today(latc2, lonc2, b['rotation'])
    az_today = azimuth_today(d['vector_today'], lat_today, lon_today)[2]
    regions8 = np.load(REGIONS, allow_pickle=False)['region']
    code = lambda la, lo: (lambda k: b['regions'][k - 1]['code'] if 1 <= k <= len(b['regions']) else None)(
        int(regions8[min(int((90.0 - la) * 8), 1439), int((lo % 360.0) * 8) % 2880]))
    places = [dict(lat_today=round(float(lat_today.ravel()[i]), 1), lon_today=round(float(lon_today.ravel()[i]) % 360.0, 1),
                   region=code(float(lat_today.ravel()[i]), float(lon_today.ravel()[i])),
                   convergence_m_per_1000_yr=round(float(conv.ravel()[i]) * 1000.0, 3),
                   drift_m2_yr=round(float(np.hypot(u, v).ravel()[i]), 1), aridity=round(float(b['aridity_c'].ravel()[i]), 3)) for i in order]
    bands = {}
    w = np.cos(phi) * np.ones_like(u)
    for hemi, sign in (('north', 1), ('south', -1)):
        for lo, hi in BANDS:
            m = land & ((sign * lat >= lo) & (sign * lat < hi))[:, None]
            if not m.any():
                continue
            e, n = float((u[m] * w[m]).sum()), float((v[m] * w[m]).sum())
            bands[f'{hemi} {lo}-{hi}'] = dict(towards_deg_turned_frame=round(math.degrees(math.atan2(e, n)) % 360.0),
                                              steadiness=round(math.hypot(e, n) / max(float((np.hypot(u, v)[m] * w[m]).sum()), 1e-12), 2),
                                              converging_share_of_land=round(float(w[m][conv[m] > 0].sum() / w[m].sum()), 2))
    return dict(strongest_convergence=places, by_band=bands,
                land_convergence_m_per_1000_yr=dict(p50=round(float(np.median(conv[land])) * 1000.0, 4), p90=round(float(np.percentile(conv[land], 90)) * 1000.0, 3),
                                                    highest=round(float(conv[land].max()) * 1000.0, 3)),
                note='the convergence is what a flux at the wind\'s capacity would lay down; the supply found here keeps the flux far below capacity')


def timescales(b):
    """Part E's times: how long the supply found here needs for a body of sand of a given size."""
    sup = SUPPLY['central']
    out = {}
    for lat in (50.0, 70.0):
        y = float(erosion_m_per_yr(lat)) * SAND_SIZE_SHARE * sup['durable'] * ROCK_DENSITY / BULK_DENSITY * 1e6   # m3 of sand per km2 a year
        out[f'river_sand_m3_per_km2_yr_at_{lat:g}'] = round(y, 2)
    dry_rain = 12.0                                            # mm a year, the moisture map's hyper-arid mean
    typical = float(np.interp(70.0, b['band_rain_mm_yr']['centres_deg'], b['band_rain_mm_yr']['mean_mm_yr']))
    rain_rule = (dry_rain / typical) ** EROSION_RAIN_POWER
    out['land_rain_at_70_mm_yr'] = round(typical, 1)
    out['rain_rule_at_12_mm_at_70'] = round(rain_rule, 3)
    river = float(erosion_m_per_yr(70.0)) * rain_rule * SAND_SIZE_SHARE * sup['durable'] * ROCK_DENSITY / BULK_DENSITY * 1e6 * sup['catchment_cap_km2']
    shore = sup['shore_m2_yr'] * 100e3
    dp = b['drifts']['central']['dp_m2_yr']
    typical = float(np.median(dp[b['land_c']]))
    for label, km3 in (('dune field, 100 km2 at 10 m', 1.0), ('sand sea, 10,000 km2 at 20 m', 200.0)):
        out[label] = dict(volume_km3=km3,
                          years_from_one_river_end=round(km3 * 1e9 / river), years_from_100_km_of_downwind_shore=round(km3 * 1e9 / shore),
                          years_at_the_winds_capacity_across_100_km=round(km3 * 1e9 / (typical * 100e3)))
    out['basis'] = dict(river_end_m3_yr=round(river), catchment_km2=sup['catchment_cap_km2'], shore_m3_yr=round(shore), capacity_m2_yr=round(typical, 1),
                        note='a river end in the hyper-arid land at 70 deg (12 mm of rain a year) with the capped catchment, all its durable sand taken; a downwind sea shore; the median central drift potential over land')
    return out


def ring_evidence():
    base = ROOT / 'climate' / 'results' / 'crm'
    out, hashes = {}, {}
    for name, what in RINGS.items():
        path = base / name
        d = json.loads(path.read_text())
        hashes[f'climate/results/crm/{name}'] = sha256(path)
        if 'wind_10m_m_s' in d:
            land = d['wind_10m_m_s']['land']
            entry = dict(what=what, land_10m_m_s={k: round(float(v), 2) for k, v in land.items()},
                         rain_land_mm_day=round(float(d['rain_mm_day']), 4) if not isinstance(d.get('rain_mm_day'), dict) else None)
            storms = d.get('storms') or {}
            if storms.get('gust_10m_max_m_s'):
                entry['storms'] = dict(count=storms['count'], per_snapshot=round(float(storms['per_snapshot']), 3),
                                       width_km_median=round(float(storms['width_km']['median']), 1),
                                       gust_10m_max_m_s={k: round(float(v), 2) for k, v in storms['gust_10m_max_m_s'].items()})
            hour = d.get('by_hour_angle', {})
            if hour.get('land', {}).get('wind_10m_m_s'):
                entry['by_hour_angle'] = dict(hour_angle_deg=hour['hour_angle_deg'], land_10m_m_s=hour['land']['wind_10m_m_s'],
                                              mixed_layer_km=hour['land'].get('mixed_layer_km'))
        else:
            entry = dict(what=what, land_10m_m_s={k: v for k, v in d['box']['wind_10m_m_s'].items()},
                         rain_mm_per_lunar_day=d['box']['rain_mm_per_lunar_day'], mean_height_m=d['box']['mean_height_m'])
        out[name] = entry
    return out, hashes


EVIDENCE = ('A first reckoning, nothing in it observed. The thresholds and the flux law are published relations fitted on Earth '
            'and carried to lunar gravity; no experiment at 1.62 m/s2 in dense air tests them. The winds are 3-day means of the '
            'design climate\'s run at T21 (about 170 km), used as the run gives them, carried from its 2 m land roughness to a sand '
            'surface by a stated factor, with the unresolved part as a stated range from the two-dimensional cloud-resolving rings '
            'at 45 deg and the highland box. The supply (mantle, river ends, shores) and the holding (plants, crusts, damp) are '
            'brackets from the project\'s studies and from Earth\'s literature. The map is a screening scheme on the 1.9 km grid, '
            'not a transport model: it places the offered sand, carries it downwind as far as the drift allows in the land\'s age, '
            'fanning out where the winds vary, and stops it at shelter, cover and water. Areas are what the scheme gives for each '
            'stated case.')
READING_RULE = ('Grids are 2,880 rows from 90 N by 5,760 columns from 0 E, 16 per degree, on today\'s Moon (the design surface\'s '
                'grid). sand_share is the share of the cell\'s ground under loose sand (0 over water); active_share the share of '
                'that sand that is bare and moves within a year; form 0 none, 1 sheet or streaks (and sand on the floors of small '
                'bowls), 2 crescent dunes and transverse ridges, 3 linear dunes, 4 star dunes, 5 dunes held under plants or crust; '
                'crest_azimuth_deg the direction the crests run, clockwise from today\'s north, between 0 and 180 (NaN for forms 0, '
                '1 and 4); crest_spacing_m the distance between crests (0 where no dunes); drift_azimuth_deg the direction the '
                'sand is carried towards, clockwise from today\'s north (NaN where nothing moves). These hold the central case; '
                'each variant has the same fields with the suffix __<variant>. Also kept: sand_thickness_m (the sand spread evenly '
                'over the cell, bulk metres); drift_potential_m2_yr (what all winds together could carry over bare loose sand, bulk '
                'm3 per metre of width a year, for the central case and for __nogust and __highgust), drift_resultant_m2_yr (the '
                'net of it, central) and drift_direction_deg (the direction that net travels, clockwise from today\'s north, on '
                'all land, central; NaN over water); '
                'basalt_share (the mare map\'s share of the rock under the cell), plant_cover (the globe\'s cover curve at the '
                'cell\'s aridity) and shelter (0 open, 1 in the lee of higher ground or on a floor).')


def write(log=print):
    BULK.mkdir(parents=True, exist_ok=True)
    partial = PRODUCT.with_suffix('.partial.npz')
    writer = NpzWriter(partial)
    b = build(writer.add, log)
    rings, ring_hashes = ring_evidence()
    inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (MOISTURE,)}
    inputs.update(ring_hashes)
    ignored = {str(p): sha256(p) for p in (SERIES, SURFACE, ROUTING, ROCK, REGIONS, REGION_RECORD)}
    producer = dict(files={c: sha256(ROOT / c) for c in CODE},
                    constants=dict(MOON_RADIUS=MOON_RADIUS, MOON_SURFACE_GRAVITY=MOON_SURFACE_GRAVITY, STANDARD_GRAVITY=STANDARD_GRAVITY,
                                   JULIAN_YEAR_DAYS=JULIAN_YEAR_DAYS))
    variants = {k: {kk: vv for kk, vv in v.items() if kk != 'regions'} for k, v in b['summary'].items()}
    meta = dict(schema=SCHEMA, grid=dict(pixels_per_degree=PPD, rows=ROWS, cols=COLS),
                frame=dict(grid='today\'s Moon, rows from 90 N, columns from 0 E', winds='the run\'s grid on the turned Moon, carried here with places and directions turned',
                           rotation_today_to_turned=b['surface_meta']['frame']['rotation_today_to_turned'], record=b['surface_meta']['frame']['record']),
                producer=producer, inputs=inputs, ignored_inputs=ignored, evidence=EVIDENCE, reading_rule=READING_RULE,
                forms=list(FORM_NAMES), variants={k: v['rests_on'] for k, v in VARIANTS.items()}, central='central',
                areas_km2={k: dict(sand=v['sand_km2'], moving=v['moving_km2'], dunes=v['dunes_km2']) for k, v in variants.items()})
    writer.add('metadata', np.array(json.dumps(meta)))
    writer.close()
    partial.replace(PRODUCT)
    record = dict(schema=SCHEMA, producer=producer, inputs=inputs, ignored_inputs=ignored, evidence=EVIDENCE, reading_rule=READING_RULE,
                  run=dict(name=RUN.name, years=list(YEARS), records=int(b['series']['uas'].shape[0]),
                           fluxmod_namelist=b['series']['metadata']['fluxmod_namelist'], landmod_namelist=b['series']['metadata']['landmod_namelist'],
                           note='the 10 m wind (uas, vas) is PlaSim\'s lowest level\'s wind (sigma 0.99415, about 313 m) brought down on '
                                'Geleyn\'s profile with the model\'s own stability factor over its land roughness of 2 m (NWIND10 1); the '
                                'stress (tauu, tauv) is its bulk transfer\'s with Beljaars\'s gust (NGUST 1); both are means over the 143 '
                                'steps of each output, 2.98 days'),
                  physics=physics(), winds=wind_tables(b), cloud_resolving=rings, pathways=pathways(b), timescales=timescales(b),
                  land=dict(km2=round(b['land_km2']), by_class_km2=b['class_km2'], band_rain_mm_yr=b['band_rain_mm_yr'], **b['counts']),
                  parameters=dict(subgrid=SUBGRID, supply=SUPPLY, to_sand=TO_SAND, ring_to_sand=RING_TO_SAND, storm_gust=STORM_GUST, storm_time=STORM_TIME,
                                  storm_rain_mm_day=STORM_RAIN, cover=COVER, cover_efold=COVER_EFOLD, crust_factor=CRUST_FACTOR,
                                  crust_edge=CRUST_EDGE, crust_steep=CRUST_STEEP, fresh_share=FRESH_SHARE, fresh_take=FRESH_TAKE,
                                  dry_hours=DRY_HOURS, erosion_m_per_myr=EROSION_M_PER_MYR, erosion_centres_deg=EROSION_CENTRES_DEG, erosion_range=EROSION_RANGE, erosion_rain_power=EROSION_RAIN_POWER, side_share=SIDE_SHARE, source_thickness_m=SOURCE_THICKNESS_M, settle_sigmas=SETTLE_SIGMAS, final_sigma=FINAL_SIGMA, climate_grid=CLIMATE_GRID, climate_smoothing_cells=CLIMATE_SMOOTHING_CELLS, step_limit=STEP_LIMIT, step_shelter=STEP_SHELTER, sand_size_share=SAND_SIZE_SHARE,
                                  bowl_share=BOWL_SHARE, bowl_depth_m=BOWL_DEPTH_M, full_cover_m=FULL_COVER_M, dune_thickness_m=DUNE_THICKNESS_M,
                                  mobile_thickness_m=MOBILE_THICKNESS_M, trap=[TRAP_BASE, TRAP_GAIN], shelter_reach_cells=SHELTER_REACH,
                                  shelter_slope=SHELTER_SLOPE, floor=[FLOOR_TPI_M, FLOOR_SLOPE], max_steps=MAX_STEPS, bulk_density_kg_m3=BULK_DENSITY,
                                  impact_ratio=IMPACT_RATIO, flux_c=FLUX_C, shao_lu=dict(a=SHAO_LU_A, gamma=SHAO_LU_GAMMA)),
                  variants=b['summary'], product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)),
                  limits=['The winds are one run\'s, at 170 km and 3 days; sea and lake breezes, slope winds and the storms\' gusts enter only through the stated factor and tail.',
                          'The cover is the globe\'s display curve at each cell\'s aridity; wetter ground along rivers and in hollows is left out.',
                          'Rivers\' sand is credited by catchment behind each mouth with a cap; lakes upstream that hold sand back are not followed.',
                          'Shore sand is a bracket per metre of downwind shore; no wave or tide model stands behind it.',
                          'The carrying scheme moves sand along one direction per cell, each case\'s resultant, for at most 200 moves from cell to cell (380 km along a meridian); it fans out by a stated share at each move and is shared over 4.7 km at the end.',
                          'The climate\'s fields reach the 1.9 km grid smoothed over 0.7 of the run\'s 170 km spacing, as the round\'s look carries them: nothing in the map is finer than that except what the ground gives (water, shelter, river ends and shores).',
                          'Dune forms follow the resultant over the potential and the sand\'s thickness by Earth\'s classes; crest spacings are an estimate from the sand at hand and the time.'],
                  seconds=round(b['seconds']))
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    log(f'  wrote {OUT} and {PRODUCT} ({PRODUCT.stat().st_size / 1e6:.0f} MB)')
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('command', nargs='?', default='build', choices=('build', 'reduce'))
    args = parser.parse_args(argv)
    if args.command == 'reduce':
        reduce_winds()
        return 0
    record = write()
    for name, v in record['variants'].items():
        print(f'{name}: sand {v["sand_km2"]:,} km2 ({100 * v["share_of_land"]:.2f}% of the land), moving {v["moving_km2"]:,}, dunes {v["dunes_km2"]:,}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
