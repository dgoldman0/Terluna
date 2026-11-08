"""The Open Moon's moisture map: the land classed by its water balance on the design climate of the turned Moon, with
the warmth, night light and settings beside each class, and the years, the GCM-CM1 gap and windows of years carried as
its uncertainty.

The water. Each T21 land cell of the turned run at the target (A28R_301k_design_l15_snow_s924_turned_refit, years 10-49)
takes its aridity, the rain over a Priestley-Taylor demand (alpha 1.26) on the cell's own mean net radiation at the
ground, its air at 2 m and its pressure, read as the design run's land readings read it: each field's means by tenth of
the lunar day, placed by each 3-day mean's hour angle (climate/crm/gcm_site.py), then averaged over the ten tenths. The
classes are UNEP's aridity classes (Middleton & Thomas 1997, as the Millennium Ecosystem Assessment's dryland chapter
tabulates them, Safriel & Adeel 2005, table 22.1), with the humid class divided where the rain meets the demand (aridity
1, below which a place has a yearly moisture deficit). The demand is radiation-based where the World Atlas's is
Thornthwaite's temperature-based one: on the Moon the warmth and the sunlight do not keep Earth's relation, so the
boundaries are UNEP's and the demand is the project's. The classes are named by their moisture, and beside each stands
the broad biome the dryland chapter gives the class on Earth, as a reference only. Which organisms each class holds on
the Moon is a design question, and none enters the map; a biome map comes when the moisture is joined by further axes.

The frame. The run lives in the turned Moon's frame, geology's orientation product at the elastic state: its cells,
their latitudes and the Earth above the sub-Earth point at 0 N 0 E are the turned Moon's. Geography's grids stand on
today's ground, so each pixel of the 28% atlas (4 px/deg) with its tides, and each rain-fed lake pixel (16 px/deg), is
counted in the cell that holds it on the turned Moon, carried by the orientation's rotation (climate/gcm/turned.py);
the pixels keep today's areas, which are the ground's. The lakes come from the drainage routed on the same years'
climatology carried back to today's ground: geography/turned_y10-49 for the forty years, geography/turned for 10-19.

The warmth. The 2 m air's mean, its coldest and warmest 3-day means and the interval minima PlaSim writes (mint)
recheck ecology finding 15: no frost, alpine, boreal or seasonal biome on natural ground. A cell whose mean air lies
below 24 C, the lower bound of Holdridge's tropical belt (recalled, Holdridge 1967), is marked upland.

The night. Each cell's darkest natural light comes from the living-light study's grid of the calendar's evaluator
(research/studies/living_light/results/living_light.json, natural.grid: clear sky, level ground, no libration),
interpolated in the cell's latitude and the Earth's mean elevation there; where the Earth stands below the horizon the
twilight alone counts. Night regimes: dark below 0.01 lux (Earth-like darkness, below the 0.01-0.03 lux at which light
at night suppresses melatonin in fishes and rodents, Grubisic et al. 2019, ecology finding 5), dim to 1 lux, and
Earthlit or twilit from 1 lux. Plants' dark hours (photosynthetic light under 1 umol/m2/s) come from the plant model's
three latitude bands (biosphere/canopy/results/plant.json), falling to none at 65 deg (finding 15).

The settings. The atlas, the rain-fed lakes and the tides give each cell its seas, lakes and closed lakes, its coast's
length and the monthly tidal range along it; the fresh-shore belts are the coast's length times the shore screen's
flat-shore widths, 0.6-2.5 km (ecology README, 15 October). Closed lakes in dry cells are the soda lakes' candidates
(register, 15 October). The soils enter as the soils study's zones by aridity class and the storm-water study's four
surfaces side by side; neither is chosen here (register, 3 November).

The uncertainty, each source kept apart and then combined:
- the years: the two-year spans' aridity gives the mean's standard error, and the band is two standard errors either
  side; 400 resamples of the run's years drawn with replacement give the shares' sampling spread and each cell's share
  of resamples that keep its class;
- the gap: CM1 rains 0.83 of the GCM's rain at the three gap boxes forced by the design years of 1 November, so the rain
  is scaled uniformly by 0.83, a sensitivity of the classes and never a bias-corrected climate; the boxes have not been
  forced by the turned years, so the scale is the earlier climate's; the aridity's own ratio at the two land boxes
  (0.90-0.92, rule R1b) is read beside it;
- combined: a cell is robust where its class holds from 0.83 times the band's low end to the band's high end;
- the timing: the GCM put 50% of the Near-side Sea box's rain into the night against CM1's 10% on the design years, so a
  cell where the GCM puts a quarter or more of its rain into the night carries a night share stated as 0.19-1 times the
  GCM's;
- the windows: the four decades 10-19, 20-29, 30-39 and 40-49 are classed alone, so the share of land whose class
  differs between independent ten-year windows of one climate measures how far a ten-year map moves by chance; years
  10-19 are mapped whole beside the forty years, with their own lakes;
- the hollows' store runs belong to the design run of 1 November and its land, so they do not enter.

    /home/kir/Documents/Projects/terluna/climate/gcm/.venv/bin/python -m biosphere.ecology.moisture_map
    ... -m biosphere.ecology.moisture_map --cache DIR     # the year files' sums kept, so the windows need no rereading
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import EARTH_MOON_DISTANCE, LEGACY_MOON_GRAVITY, MOON_RADIUS
from shared.provenance import constants_used
from climate.gcm.turned import read_frame, to_turned

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'results' / 'moisture_map.json'
RUNS = ROOT / 'climate' / 'gcm' / 'runs'
SCHEMA = 'terluna.biosphere.ecology-moisture-map/1'
DRIVE = Path('/media/projectspace/terluna-research/session-tools')
DESIGN_YEARS = DRIVE / 'design_run' / 'design_years.json'           # the established land readings, checked against
CLASS_SENSITIVITY = DRIVE / 'design_run' / 'class_sensitivity.json'

RUN, YEARS = 'A28R_301k_design_l15_snow_s924_turned_refit', (10, 49)        # the turned run at the target
COMPARE_YEARS = (10, 19)                # the turned design years, mapped whole beside the forty years
DECADES = ((10, 19), (20, 29), (30, 39), (40, 49))
DESIGN_RUN = 'A28R_301k_design_l15_snow_s934'          # the design run of 1 November, whose land readings check it
STORE_RUNS = ()                         # the hollows' store runs belong to the design run of 1 November
STORE_YEARS = (0, 9)            # ten years from the design run's year-9 restart, the design years' counterparts
SPAN = 2                        # years in each span of the years' spread

GEOGRAPHY = {'atlas': 'geography/results/atlas.json', 'atlas_grid': 'geography/products/atlas_28pct_4ppd.npz',
             'tides_grid': 'geography/products/tides_28pct_4ppd.npz'}
DRAINAGE = {(10, 49): 'geography/turned_y10-49', (10, 19): 'geography/turned'}     # each window's lakes
ORIENTATION = 'geology/results/moon_orientation.json'
TURNED_YEARS_LAND = DRIVE / 'turned_design' / 'turned_years_land_A28R_301k_design_l15_snow_s924_turned_refit.json'
LIVING_LIGHT = 'research/studies/living_light/results/living_light.json'
PLANT = 'biosphere/canopy/results/plant.json'
SCHEMAS = {LIVING_LIGHT: 'terluna.research.living-light/1', GEOGRAPHY['atlas']: 'terluna.geography.atlas/1'}

# The gap boxes' readings on the design years (session-tools/agents/gap_boxes/read_rules_*.out): lunar-day rain in
# mm/day, CM1 then the GCM; the land boxes' aridity (R1b); each box's share of rain from sunset to sunrise (R2a's
# third and fourth quarters).
GAP_RAIN = {'Near-side Sea, 5 S 45 W': (1.17442, 1.41976), 'land, 0 N 0 E': (4.97434, 6.00912),
            'dry box, 44.3 S 56.25 E': (0.51337, 0.60936)}
GAP_ARIDITY = {'land, 0 N 0 E': (1.73568, 1.89613), 'dry box, 44.3 S 56.25 E': (0.27348, 0.30562)}
GAP_NIGHT_SHARE = {'Near-side Sea, 5 S 45 W': (0.0948, 0.5039), 'land, 0 N 0 E': (0.0075, 0.0214),
                   'dry box, 44.3 S 56.25 E': (0.0578, 0.0246)}
GAP_SCALE = round(float(np.mean([c / g for c, g in GAP_RAIN.values()])), 2)          # 0.83
NIGHT_FLAG = 0.25                       # a quarter of the rain at night: the sea box's regime (0.50), not the land's
NIGHT_LOW_RATIO = round(GAP_NIGHT_SHARE['Near-side Sea, 5 S 45 W'][0] / GAP_NIGHT_SHARE['Near-side Sea, 5 S 45 W'][1], 2)

MM = 86400e3                            # m/s of water to mm/day
LATENT, CP_AIR = 2.5e6, 1004.0          # J/kg, J/kg/K, as the design run's land readings take them
KEYS = ('pr', 'evap', 'mrso', 'rss', 'rls', 'tas', 'ps', 'ts')
CENTRES = np.array([-180.0 + 36.0 * (k + 0.5) for k in range(10)])        # each tenth's hour angle, degrees past noon
NIGHT = np.where(np.abs(CENTRES) > 90.0, 1.0, np.where(np.abs(CENTRES) == 90.0, 0.5, 0.0))

UNEP = (('hyper-arid', 0.0, 0.05), ('arid', 0.05, 0.2), ('semi-arid', 0.2, 0.5), ('dry sub-humid', 0.5, 0.65),
        ('humid', 0.65, math.inf))
WATER_BOUNDS = (0.05, 0.2, 0.5, 0.65, 1.0)
# The moisture classes: UNEP's aridity classes with the humid class divided where the rain meets the demand. Beside each
# stands the biome the dryland chapter gives the class on Earth (Safriel & Adeel 2005, table 22.1), as a reference; which
# organisms each class holds on the Moon is a later design question, so none enters the map.
CLASSES = (
    dict(key='hyper_arid', name='hyper-arid', aridity=(0.0, 0.05), unep='hyper-arid', earth_biome='hyper-arid desert',
         earth_broad_biome='desert', evidence='published relation (UNEP class, MA 2005 table 22.1)'),
    dict(key='arid', name='arid', aridity=(0.05, 0.2), unep='arid', earth_biome='arid desert and shrubland',
         earth_broad_biome='desert', evidence='published relation (UNEP class, MA 2005 table 22.1)'),
    dict(key='semi_arid', name='semi-arid', aridity=(0.2, 0.5), unep='semi-arid', earth_biome='semi-arid grassland',
         earth_broad_biome='grassland', evidence='published relation (UNEP class, MA 2005 table 22.1)'),
    dict(key='dry_subhumid', name='dry sub-humid', aridity=(0.5, 0.65), unep='dry sub-humid',
         earth_biome='dry sub-humid woodland', earth_broad_biome='forest (mainly woodland)',
         evidence='published relation (UNEP class, MA 2005 table 22.1)'),
    dict(key='humid', name='humid', aridity=(0.65, 1.0), unep='humid', earth_biome='humid forest',
         earth_broad_biome='forest', evidence='published relation (UNEP humid class)'),
    dict(key='wet', name='wet', aridity=(1.0, math.inf), unep='humid', earth_biome='wet forest', earth_broad_biome='forest',
         evidence='published relation (aridity 1, where the rain meets the demand; MA 2005)'),
)
UPLAND_C = 24.0                         # Holdridge's tropical belt's lower bound, biotemperature 24 C (recalled)
LUX_DARK, LUX_BRIGHT = 0.01, 1.0
NIGHT_REGIMES = ('dark', 'dim', 'Earthlit', 'twilit')
SHORE_BELT_KM = (0.6, 2.5)              # flat shores' groundwater belts, the shore screen of 15 October
TIDE_CLASSES = (('microtidal', 0.0, 2.0), ('mesotidal', 2.0, 4.0), ('macrotidal', 4.0, math.inf))   # Davies (1964)
SOIL_ZONES = {   # the soils study of 3 November (session-tools/agents/soils/soils_maintenance.md), its table
    'hyper-arid': dict(soils='Leptosols, Regosols under a rock-fragment pavement; Calcisols where water gathers',
                       kept_by='the pavement left whole; traffic on set ways', ease='easy where left alone'),
    'arid': dict(soils='Calcisols, Arenosols of feldspar-and-glass sand, Durisols',
                 kept_by='a biological crust over the open ground; disturbed patches reseeded; grazers steered off',
                 ease='easy to moderate'),
    'semi-arid': dict(soils='Calcisols, Kastanozems, Vertisols on mare remnants and in hollows',
                      kept_by='plants and litter over about 35-60% of the ground, crusts between; hollows as planting '
                              'pits; rotated grazing', ease='moderate'),
    'dry sub-humid': dict(soils='Chernozems, Phaeozems, Luvisols, Lixisols, Vertisols',
                          kept_by='a nearly closed perennial cover of about 50-80% with litter; managed disturbance '
                                  'in place of fire', ease='moderate'),
    'humid': dict(soils='Andosol-like young soils, Cambisols, Nitisols on basalt, aging towards Lixisols and Acrisols; '
                        'Gleysols, Histosols, Fluvisols on wet ground',
                  kept_by='forest and wetland cover; nutrient return by the project\'s loops; conservation farming',
                  ease='easy to hold, sustained effort for fertility'),
}
STORM_WATER_SURFACES = {   # recharge, % of rain, Walnut Gulch's storms, young land, clay floors / sandy banks
    'source': 'the storm-water study on domain/hydrosphere, as the soils study tabulates it (section 3.3)',
    'sandy loam': {'arid': (1.3, 1.4), 'semi-arid': (2.1, 2.4), 'dry sub-humid': (3.1, 3.4)},
    'biological crust': {'arid': (1.5, 1.8), 'semi-arid': (2.4, 2.8), 'dry sub-humid': (3.5, 3.9)},
    'loam': {'arid': (5.9, 7.3), 'semi-arid': (9.0, 10.5), 'dry sub-humid': (12.1, 13.7)},
    'sealed loam': {'arid': (17.2, 21.5), 'semi-arid': (23.9, 28.2), 'dry sub-humid': (29.7, 34.3)},
}
EARTH_DRYLANDS = dict(source='Safriel & Adeel (2005), Millennium Ecosystem Assessment, Dryland Systems, table 22.1',
                      shares_of_land={'hyper-arid': 0.066, 'arid': 0.106, 'semi-arid': 0.152, 'dry sub-humid': 0.087},
                      drylands=0.413)


def sha(path, n=64):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:n]


# ==== Classes ========================================================================================================
def unep_class(ai):
    """UNEP's five aridity classes, 0 hyper-arid to 4 humid, each [lower, upper) as the design readings take them."""
    return np.digitize(np.asarray(ai, float), [b for _, b, _ in UNEP[1:]])


def moisture(ai):
    """The six moisture classes, 0 hyper-arid to 5 wet: UNEP's classes with the humid class split at aridity 1."""
    return np.digitize(np.asarray(ai, float), WATER_BOUNDS)


def night_regime(darkest_lux, earth_elevation_deg):
    """0 dark (< 0.01 lux at the darkest), 1 dim (< 1 lux), 2 Earthlit or 3 twilit (1 lux and more)."""
    lux, e = np.asarray(darkest_lux, float), np.asarray(earth_elevation_deg, float)
    return np.where(lux < LUX_DARK, 0, np.where(lux < LUX_BRIGHT, 1, np.where(e > 0.0, 2, 3)))


def tide_class(range_m):
    r = np.asarray(range_m, float)
    return np.where(np.isnan(r), -1, np.digitize(np.nan_to_num(r), [lo for _, lo, _ in TIDE_CLASSES[1:]]))


def bands(ai, ai_se, scale=GAP_SCALE):
    """The aridity's band ends: the low end the rain scaled by the gap's ratio less two standard errors of the
    mean over the years, the high end two standard errors above the run's own mean."""
    ai, se = np.asarray(ai, float), np.asarray(ai_se, float)
    return np.maximum(scale * (ai - 2.0 * se), 0.0), ai + 2.0 * se


def robustness(ai, ai_se, store_ai=(), scale=GAP_SCALE):
    """Each cell's class as run, at the gap's scale, across the years' band and combined, and on the store runs;
    robust where the class holds across the source's band."""
    ai = np.asarray(ai, float)
    run = moisture(ai)
    lo, hi = bands(ai, ai_se, scale)
    years_lo, years_hi = moisture(np.maximum(ai - 2.0 * np.asarray(ai_se, float), 0.0)), moisture(hi)
    out = dict(run=run, gap=moisture(scale * ai), years_low=years_lo, years_high=years_hi,
               combined_low=moisture(lo), combined_high=moisture(hi))
    out['robust_gap'] = out['gap'] == run
    out['robust_years'] = (years_lo == run) & (years_hi == run)
    out['robust'] = (out['combined_low'] == run) & (out['combined_high'] == run)
    stores = [moisture(s) for s in store_ai]
    out['store'] = stores
    out['robust_store'] = np.all([s == run for s in stores], axis=0) if stores else np.ones(run.shape, bool)
    out['robust_all'] = out['robust'] & out['robust_store']
    return out


def shares(classes, weight, n):
    w = np.asarray(weight, float) / np.sum(weight)
    return [float(w[np.asarray(classes) == k].sum()) for k in range(n)]


def margins(ai):
    """The aridity's distance to its class's lower and upper bounds, as fractions of the aridity (how far the rain
    must fall or rise, uniformly, to change the class)."""
    ai = np.asarray(ai, float)
    k = moisture(ai)
    lower = np.array((0.0,) + WATER_BOUNDS)[k]
    upper = np.array(WATER_BOUNDS + (math.inf,))[k]
    with np.errstate(divide='ignore', invalid='ignore'):
        down = np.where(lower > 0.0, 1.0 - lower / ai, np.inf)
        up = np.where(np.isfinite(upper), upper / ai - 1.0, np.inf)
    return down, up


# ==== The climate ====================================================================================================
def demand(m):
    """Priestley-Taylor's demand (mm/day) on the mean net radiation at the ground, the 2 m air and the pressure, as the
    design run's land readings compute it."""
    t = m['tas'] - 273.15
    es = 0.6108 * np.exp(17.27 * t / (t + 237.3))                                   # kPa
    s = 4098.0 * es / (t + 237.3) ** 2
    gamma = CP_AIR * m['ps'] * 100.0 / 1000.0 / (0.622 * LATENT)                     # ps in hPa; kPa/K
    return np.maximum(1.26 * s / (s + gamma) * (m['rss'] + m['rls']) / LATENT * 86400.0, 1e-3)


def read_year(run, year, humidity=True):
    """One year file's sums by tenth of the lunar day, its extremes and its night humidity at the lowest level."""
    import netCDF4
    from climate.crm import gcm_site as gs
    from climate.gcm import year_files as yf
    with netCDF4.Dataset(RUNS / run / 'model' / f'MOST.{year:05d}.nc') as d:
        lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
        land = np.asarray(d['lsm'][:], float)[0] > 0.5
        f = {k: np.asarray(d[k][:], float) for k in KEYS}
        czen = np.asarray(d['czen'][:], float)
        mint, maxt = np.asarray(d['mint'][:], float), np.asarray(d['maxt'][:], float)
        digest = hashlib.sha256()
        for a in (lat, lon, np.asarray(d['lsm'][:], float), czen, mint, maxt, *(f[k] for k in KEYS)):
            digest.update(np.ascontiguousarray(a).tobytes())
        hour = (lon[None, :] - gs.subsolar_longitude(czen, lat, lon)[:, None] + 180.0) % 360.0 - 180.0
        if humidity:
            keep = yf.records(d)
            q = (np.asarray(d[yf.HUS_MEAN][:, -1], float) if yf.HUS_MEAN in d.variables
                 else yf.interval_means(d['hus'][:, -1]))                          # the lowest level's 3-day means
            ta = np.asarray(d['ta'][:, -1], float)[keep]
            sigma = float(d['lev'][-1])
            p_hpa = sigma * f['ps'][keep]
            e = q * p_hpa / (0.622 + 0.378 * q)
            tc = ta - 273.15
            rh = e / (6.112 * np.exp(17.67 * tc / (tc + 243.5)))                    # Bolton (1980) over water
            night_rec = (np.abs(hour[keep]) > 90.0)[:, None, :]
            rh_sum = np.where(night_rec, rh, 0.0).sum(axis=0)
            rh_n = np.broadcast_to(night_rec, rh.shape).sum(axis=0)
            digest.update(np.ascontiguousarray(q).tobytes())
            digest.update(np.ascontiguousarray(ta).tobytes())
    tenth = np.floor((hour + 180.0) / 36.0).astype(int) % 10
    acc = {k: np.zeros((10,) + land.shape) for k in KEYS}
    cnt = np.zeros((10,) + land.shape)
    for k in range(10):
        sel = (tenth == k)[:, None, :]
        for v in KEYS:
            acc[v][k] += np.where(sel, f[v], 0.0).sum(axis=0)
        cnt[k] += np.broadcast_to(sel, f['pr'].shape).sum(axis=0)
    out = dict(lat=lat, lon=lon, land=land, acc=acc, cnt=cnt, tas_min=f['tas'].min(axis=0), tas_max=f['tas'].max(axis=0),
               mint=mint.min(axis=0), maxt=maxt.max(axis=0), digest=digest.hexdigest()[:32])
    if humidity:
        out.update(rh_sum=rh_sum, rh_n=rh_n)
    return out


def combine(years):
    """Means by tenth and their mean over the tenths, from read_year's sums over the given years."""
    acc = {k: sum(y['acc'][k] for y in years) for k in KEYS}
    cnt = sum(y['cnt'] for y in years)
    tenth = {k: v / np.maximum(cnt, 1) for k, v in acc.items()}
    return tenth, {k: v.mean(axis=0) for k, v in tenth.items()}


def aridity(mean):
    rain = mean['pr'] * MM
    pet = demand(mean)
    return rain, pet, rain / pet


CACHE = None                            # a folder for read_year's sums (--cache), so a rebuild need not reread the run


def cached_year(run, year, humidity):
    if CACHE is None:
        return read_year(run, year, humidity)
    path = Path(CACHE) / f'{run}.{year:05d}.{int(humidity)}.npz'
    if path.exists():
        z = np.load(path)
        out = {k: z[k] for k in z.files if not k.startswith('acc_') and k != 'digest'}
        out['acc'] = {k: z[f'acc_{k}'] for k in KEYS}
        out['digest'] = str(z['digest'])
        out['land'] = out['land'].astype(bool)
        return out
    out = read_year(run, year, humidity)
    path.parent.mkdir(parents=True, exist_ok=True)
    flat = {k: v for k, v in out.items() if k != 'acc'}
    flat.update({f'acc_{k}': v for k, v in out['acc'].items()})
    np.savez(path, **flat)
    return out


def read_run(run, years, humidity=True):
    data = [cached_year(run, y, humidity) for y in range(years[0], years[1] + 1)]
    tenth, mean = combine(data)
    rain, pet, ai = aridity(mean)
    spans = []
    for s in range(0, len(data) - SPAN + 1, SPAN):
        spans.append(aridity(combine(data[s:s + SPAN])[1])[2])
    height = np.loadtxt(RUNS / run / 'inputs' / 'moon_topography.sra', skiprows=1).ravel().reshape(data[0]['land'].shape)
    out = dict(run=run, years=list(years), lat=data[0]['lat'], lon=data[0]['lon'], land=data[0]['land'], tenth=tenth,
               mean=mean, rain=rain, pet=pet, ai=ai, spans=np.array(spans), height=height / LEGACY_MOON_GRAVITY,
               tas_min=np.min([y['tas_min'] for y in data], axis=0), tas_max=np.max([y['tas_max'] for y in data], axis=0),
               mint=np.min([y['mint'] for y in data], axis=0), maxt=np.max([y['maxt'] for y in data], axis=0),
               files={f'{run}/model/MOST.{y:05d}.nc': dict(fields_read_sha256=d['digest'], bytes=(RUNS / run / 'model' / f'MOST.{y:05d}.nc').stat().st_size)
                      for y, d in zip(range(years[0], years[1] + 1), data)})
    if humidity:
        out['night_rh'] = sum(y['rh_sum'] for y in data) / np.maximum(sum(y['rh_n'] for y in data), 1)
    out['data'] = data
    return out


BOOTSTRAP, SEED = 400, 20261104


def bootstrap(data, n=BOOTSTRAP, seed=SEED):
    """The aridity rebuilt from the run's years drawn with replacement, n times: the sampling spread of a mean over
    this many years, each year's sums composited as the full record is."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        pick = rng.integers(0, len(data), len(data))
        out.append(aridity(combine([data[i] for i in pick])[1])[2])
    return np.array(out)


# ==== The light ======================================================================================================
def earth_elevation(lat, lon):
    """The Earth's mean elevation (degrees) over the ground, the sub-Earth point at 0 N 0 E in the mean geometry, less
    the parallax of the Moon's radius at the Earth's mean distance (0.26 deg on the horizon), as the grid has it."""
    phi, lam = np.radians(lat), np.radians(lon)
    e = np.arcsin(np.cos(phi) * np.cos(lam))
    return np.degrees(e - np.arcsin(MOON_RADIUS / EARTH_MOON_DISTANCE * np.cos(e)))


def night_light(lat2, lon2):
    """The darkest natural light (lux) at each cell centre, from the living-light grid: linear in latitude and the
    Earth's elevation on log lux where the Earth is up, the twilight alone (the grid's never-rising rows) where not."""
    from scipy.interpolate import griddata
    grid = json.loads((ROOT / LIVING_LIGHT).read_text())['natural']['grid']
    lat_g = np.array([abs(p['lat']) for p in grid])
    e_g = np.array([p['earth_elevation_deg'] for p in grid])
    lux_g = np.log10(np.array([p['darkest_lux'] for p in grid]))
    never = np.array([p['earth_class'] == 'never' for p in grid])
    e = earth_elevation(lat2, lon2)
    alat = np.abs(lat2)
    order = np.argsort(lat_g[never])
    twilight = np.interp(alat, lat_g[never][order], lux_g[never][order])
    up = ~never
    pts = np.column_stack([lat_g[up], e_g[up]])
    lin = griddata(pts, lux_g[up], (alat, np.maximum(e, e_g[up].min())), method='linear')
    near = griddata(pts, lux_g[up], (alat, np.maximum(e, e_g[up].min())), method='nearest')
    earthlit = np.where(np.isnan(lin), near, lin)
    lux = np.where(e > e_g[up].min(), np.maximum(earthlit, twilight), twilight)
    return 10.0 ** lux, e


def dark_hours(lat2):
    """Hours a cycle with photosynthetic light under 1 umol/m2/s, from the plant model's bands (0, 30 and 60 deg) and
    none poleward of 65 deg (ecology finding 15)."""
    plant = json.loads((ROOT / PLANT).read_text())['results']['near']
    hours = []
    for band in ('0', '30', '60'):
        row = plant[band]
        hours.append(next(v['dark_hours'] for v in row.values() if isinstance(v, dict) and 'dark_hours' in v))
    return np.interp(np.abs(lat2), [0.0, 30.0, 60.0, 65.0], hours + [0.0]), dict(zip(('0', '30', '60'), hours))


# ==== The settings ===================================================================================================
def cell_edges(lat):
    mid = 0.5 * (lat[1:] + lat[:-1])
    return np.concatenate([[90.0], mid, [-90.0]])                                       # north to south


def cell_index(lat_px, lon_px, lat, lon):
    edges = cell_edges(lat)
    row = np.clip(np.searchsorted(-edges, -np.asarray(lat_px), side='right') - 1, 0, lat.size - 1)
    dlon = 360.0 / lon.size
    col = np.floor(((np.asarray(lon_px) - lon[0] + 0.5 * dlon) % 360.0) / dlon).astype(int) % lon.size
    return row, col


def cell_areas_km2(lat, lon):
    edges = np.radians(cell_edges(lat))
    band = (np.sin(edges[:-1]) - np.sin(edges[1:])) * (2.0 * math.pi / lon.size) * (MOON_RADIUS / 1e3) ** 2
    return band[:, None] * np.ones((1, lon.size))


def settings(lat, lon, rotation=None, drainage=DRAINAGE[YEARS]):
    """Seas, lakes, closed lakes, the coast's length and the tide along it, for each T21 cell. Geography's pixels
    stand on today's ground; with a rotation (today to the run's frame) each is counted in the cell that holds it in the
    run's frame, and keeps its area. drainage is the folder whose products hold the lakes routed on the run's
    climate."""
    shape = (lat.size, lon.size)
    r_km = MOON_RADIUS / 1e3
    z = np.load(ROOT / GEOGRAPHY['atlas_grid'])
    atlas = json.loads((ROOT / GEOGRAPHY['atlas']).read_text())
    sea_labels = np.array(sorted({int(b['label']) for b in atlas['bodies']}))
    lab = z['water_label']
    lat4, lon4 = z['lat_deg'], z['lon_deg']
    if not lat[0] > lat[-1]:
        raise ValueError('the GCM rows must run north to south')
    sea = np.isin(lab, sea_labels)
    low_lake = (lab > 0) & ~sea
    landpx = lab == 0
    step = math.radians(0.25)
    area4 = (r_km * step) ** 2 * np.cos(np.radians(lat4))[:, None] * np.ones((1, lon4.size))
    plat, plon = lat4[:, None] * np.ones((1, lon4.size)), np.ones((lat4.size, 1)) * lon4[None, :]
    if rotation is not None:
        plat, plon = to_turned(plat, plon, rotation)
    rr, cc = cell_index(plat, plon, lat, lon)
    flat = (rr * lon.size + cc).ravel()

    def total(weights):
        return np.bincount(flat, weights=np.asarray(weights, float).ravel(), minlength=lat.size * lon.size).reshape(shape)

    cell_area4 = total(area4)
    out = dict(sea_fraction=total(area4 * sea) / cell_area4, sea_level_lake_fraction=total(area4 * low_lake) / cell_area4,
               atlas_land_fraction=total(area4 * landpx) / cell_area4)
    # the coast: land-sea pixel edges, each credited to its land pixel's cell, with the sea pixel's monthly tide
    t = np.load(ROOT / GEOGRAPHY['tides_grid'])
    if not (np.array_equal(t['lat_deg'], lat4) and np.array_equal(t['lon_deg'], lon4)):
        raise ValueError('the tide grid must be the atlas grid')
    tide = np.full(lab.shape, np.nan)
    tide[t['rows'], t['cols']] = t['monthly_range_median_m']
    coast = np.zeros(lab.shape)
    tide_len = np.zeros(lab.shape)
    tide_w = np.zeros(lab.shape)
    ew = r_km * step                                                                 # a meridional pixel edge, km
    for shift in (1, -1):                                                            # east and west neighbours
        nb_sea = np.roll(sea, -shift, axis=1)
        nb_tide = np.roll(tide, -shift, axis=1)
        edge = landpx & nb_sea
        coast += edge * ew
        tide_len += np.where(edge & np.isfinite(nb_tide), nb_tide * ew, 0.0)
        tide_w += np.where(edge & np.isfinite(nb_tide), ew, 0.0)
    for shift in (1, -1):                                                            # the next row and the last
        nb_sea = np.zeros_like(sea)
        nb_tide = np.full(tide.shape, np.nan)
        lat_edge = lat4.copy()
        if shift == 1:
            nb_sea[:-1], nb_tide[:-1] = sea[1:], tide[1:]
            lat_edge[:-1] = 0.5 * (lat4[:-1] + lat4[1:])
        else:
            nb_sea[1:], nb_tide[1:] = sea[:-1], tide[:-1]
            lat_edge[1:] = 0.5 * (lat4[1:] + lat4[:-1])
        length = (r_km * step * np.cos(np.radians(lat_edge)))[:, None] * np.ones((1, lon4.size))
        edge = landpx & nb_sea
        coast += edge * length
        tide_len += np.where(edge & np.isfinite(nb_tide), nb_tide * length, 0.0)
        tide_w += np.where(edge & np.isfinite(nb_tide), length, 0.0)
    out['coast_km'] = total(coast)
    with np.errstate(invalid='ignore', divide='ignore'):
        with_tide = total(tide_w)
        out['tide_m'] = np.where(with_tide > 0, total(tide_len) / with_tide, np.nan)
    # rain-fed lakes above sea level (16 px/deg), and the closed ones among them
    dz = np.load(ROOT / drainage / 'products' / 'drainage_28pct_16ppd.npz')
    cells, spills = dz['lake_cell'].astype(np.int64), dz['lake_spills']
    lat16 = 90.0 - ((cells // 5760) + 0.5) / 16.0
    lon16 = ((cells % 5760) + 0.5) / 16.0
    a16 = (r_km * math.radians(1.0 / 16.0)) ** 2 * np.cos(np.radians(lat16))
    if rotation is not None:
        lat16, lon16 = to_turned(lat16, lon16, rotation)
    r16, c16 = cell_index(lat16, lon16, lat, lon)
    f16 = r16 * lon.size + c16
    box = cell_areas_km2(lat, lon)
    lakes = np.bincount(f16, weights=a16, minlength=lat.size * lon.size).reshape(shape)
    closed = np.bincount(f16, weights=a16 * ~spills, minlength=lat.size * lon.size).reshape(shape)
    out.update(rainfed_lake_fraction=lakes / box, closed_lake_fraction=closed / box, closed_lake_km2=closed,
               sea_level_lake_km2=out['sea_level_lake_fraction'] * box)
    return out


# ==== The map ========================================================================================================
def build(run=RUN, years=YEARS, store_runs=STORE_RUNS, store_years=STORE_YEARS, drainage=None, orientation=ORIENTATION):
    """The map of one window of years: the run's cells, carried to today's ground through the orientation product
    when one is given (None for a run on today's ground), with the lakes of the drainage folder."""
    drainage = drainage or DRAINAGE.get(tuple(years), DRAINAGE[YEARS])
    frame = read_frame(ROOT / orientation) if orientation else None
    c = read_run(run, years)
    lat, lon, land = c['lat'], c['lon'], c['land']
    lat2 = lat[:, None] * np.ones((1, lon.size))
    lon2 = np.ones((lat.size, 1)) * lon[None, :]
    w = np.cos(np.radians(lat2))
    wl = w * land / (w * land).sum()
    stores = [read_run(s, store_years, humidity=False) for s in store_runs]
    for s in stores:
        if not np.array_equal(s['land'], land):
            raise ValueError(f"{s['run']} has another land mask")
    ai_se = c['spans'].std(axis=0, ddof=1) / math.sqrt(c['spans'].shape[0])
    rob = robustness(c['ai'], ai_se, [s['ai'] for s in stores])
    lux, elev = night_light(lat2, lon2)
    regime = night_regime(lux, elev)
    hours, hour_bands = dark_hours(lat2)
    geo = settings(lat, lon, frame['rotation'] if frame else None, drainage)
    night_share = np.tensordot(NIGHT, c['tenth']['pr'], axes=(0, 0)) / np.maximum(c['tenth']['pr'].sum(axis=0), 1e-30)
    tas = c['mean']['tas'] - 273.15
    swing = c['tenth']['tas'].max(axis=0) - c['tenth']['tas'].min(axis=0)
    upland = tas < UPLAND_C
    down, up = margins(c['ai'])
    area = cell_areas_km2(lat, lon)
    ncom, names = len(CLASSES), [k['name'] for k in CLASSES]
    L = land

    def share_rows(classes):
        return dict(zip(names, (round(v, 4) for v in shares(classes[L], w[L], ncom))))

    by_class = {}
    for k, com in enumerate(CLASSES):
        cls = L & (rob['run'] == k)
        row = dict(share=round(float(wl[cls].sum()), 4), area_km2=round(float(area[cls].sum())),
                   robust_share=round(float(wl[cls & rob['robust']].sum()), 4),
                   robust_gap_share=round(float(wl[cls & rob['robust_gap']].sum()), 4),
                   robust_years_share=round(float(wl[cls & rob['robust_years']].sum()), 4),
                   robust_all_share=round(float(wl[cls & rob['robust_all']].sum()), 4),
                   share_at_gap_scale=round(float(wl[L & (rob['gap'] == k)].sum()), 4),
                   share_at_band_low=round(float(wl[L & (rob['combined_low'] == k)].sum()), 4),
                   share_at_band_high=round(float(wl[L & (rob['combined_high'] == k)].sum()), 4))
        if cls.any():
            ww = w * cls / (w * cls).sum()
            g = lambda a: float((np.nan_to_num(a) * ww).sum())
            row.update(rain_mm_yr=round(g(c['rain']) * 365.0), demand_mm_yr=round(g(c['pet']) * 365.0),
                       aridity=round(g(c['ai']), 3), air_c=round(g(tas), 2), night_rain_share=round(g(night_share), 3),
                       latitude_shares={f'{a}-{b}': round(float((ww * ((np.abs(lat2) >= a) & (np.abs(lat2) < b))).sum()), 3)
                                        for a, b in ((0, 15), (15, 30), (30, 45), (45, 60), (60, 91))},
                       night_regimes={n: round(float(ww[regime == j].sum()), 3) for j, n in enumerate(NIGHT_REGIMES)},
                       near_side_share=round(float(ww[elev > 0].sum()), 3))
        by_class[com['key']] = row
    boot = bootstrap(c['data'])
    boot_run = np.array([moisture(a) for a in boot])
    boot_gap = np.array([moisture(GAP_SCALE * a) for a in boot])
    boot_shares = np.array([shares(k[L], w[L], ncom) for k in boot_run])
    boot_gap_shares = np.array([shares(k[L], w[L], ncom) for k in boot_gap])
    keep_prob = (boot_run == rob['run'][None]).mean(axis=0)
    for k, com in enumerate(CLASSES):
        by_class[com['key']].update(
            share_sd_years=round(float(boot_shares[:, k].std(ddof=1)), 4),
            share_years_95=[round(float(v), 4) for v in np.quantile(boot_shares[:, k], [0.025, 0.975])],
            share_at_gap_scale_years_95=[round(float(v), 4) for v in np.quantile(boot_gap_shares[:, k], [0.025, 0.975])])
    one_se = (moisture(np.maximum(GAP_SCALE * (c['ai'] - np.asarray(ai_se)), 0.0)) == rob['run']) & \
             (moisture(c['ai'] + np.asarray(ai_se)) == rob['run'])
    sea_cell = ~land
    beside_sea = np.zeros_like(land)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            shifted_sea = np.roll(sea_cell, dc, axis=1)
            if dr == 1:
                shifted_sea = np.vstack([shifted_sea[1:], np.zeros((1, lon.size), bool)])
            elif dr == -1:
                shifted_sea = np.vstack([np.zeros((1, lon.size), bool), shifted_sea[:-1]])
            beside_sea |= shifted_sea
    beside_sea &= land
    transitions = {names[a]: {names[b]: round(float(wl[L & (rob['run'] == a) & (rob['combined_low'] == b)].sum()), 4)
                              for b in range(ncom) if (L & (rob['run'] == a) & (rob['combined_low'] == b)).any()}
                   for a in range(ncom)}
    unep_run = unep_class(c['ai'])
    unep_shares = dict(zip([u[0] for u in UNEP], (round(v, 4) for v in shares(unep_run[L], w[L], 5))))
    sens = {}
    for scale in (0.8, 0.83, 0.9, 1.0, 1.1, 1.2):
        k = unep_class(c['ai'] * scale)
        sens[f'{scale:.2f}'] = dict(land_changing_class=round(float(wl[L & (k != unep_run)].sum()), 4),
                                    class_shares=dict(zip([u[0] for u in UNEP], (round(v, 4) for v in shares(k[L], w[L], 5)))))
    aridity_ratio = float(np.mean([a / b for a, b in GAP_ARIDITY.values()]))
    at_aridity_ratio = moisture(c['ai'] * aridity_ratio)
    flagged_night = L & (night_share >= NIGHT_FLAG)
    p = w[L] / w[L].sum()
    simpson = 1.0 - sum(v * v for v in shares(rob['run'][L], w[L], ncom))
    combo = rob['run'] * len(NIGHT_REGIMES) + regime
    combo_shares = shares(combo[L], w[L], ncom * len(NIGHT_REGIMES))
    checks = dict(shares_sum=round(sum(by_class[k['key']]['share'] for k in CLASSES), 6),
                  gap_never_wetter=bool(np.all(rob['gap'][L] <= rob['run'][L])),
                  robust_within_gap_and_years=bool(np.all(~rob['robust'][L] | (rob['robust_gap'][L] & rob['robust_years'][L]))))
    if DESIGN_YEARS.exists() and run == DESIGN_RUN and tuple(years) == (10, 19):
        ref = json.loads(DESIGN_YEARS.read_text())['new']
        checks['design_years_aridity_shares_max_difference'] = max(
            abs(ref['aridity'][n]['land_share'] - float(wl[L & (unep_run == j)].sum())) for j, (n, _, _) in enumerate(UNEP))
        checks['design_years_land_rain_difference_mm_day'] = abs(ref['land_rain_mm_day'] - float((c['rain'] * wl).sum()))
        checks['design_years_json_sha256'] = sha(DESIGN_YEARS)
    if TURNED_YEARS_LAND.exists() and run == RUN and tuple(years) == COMPARE_YEARS:
        ref = json.loads(TURNED_YEARS_LAND.read_text())['turned_years']['unep_share']
        checks['turned_years_aridity_shares_max_difference'] = max(
            abs(ref[n] - float(wl[L & (unep_run == j)].sum())) for j, (n, _, _) in enumerate(UNEP))
        checks['turned_years_land_json_sha256'] = sha(TURNED_YEARS_LAND)
    if CLASS_SENSITIVITY.exists() and run == DESIGN_RUN and tuple(years) == (10, 19):
        ref = json.loads(CLASS_SENSITIVITY.read_text())['by_rain_scale']
        checks['class_sensitivity_max_difference'] = max(
            abs(ref[s]['class_shares'][n] - sens[f'{float(s):.2f}']['class_shares'][n]) for s in ref for n, _, _ in UNEP)
        checks['class_sensitivity_json_sha256'] = sha(CLASS_SENSITIVITY)
    frost = dict(coldest_3day_mean_air_c=round(float(c['tas_min'][L].min()) - 273.15, 2),
                 coldest_interval_minimum_air_c=round(float(c['mint'][L].min()) - 273.15, 2),
                 warmest_3day_mean_air_c=round(float(c['tas_max'][L].max()) - 273.15, 2),
                 warmest_interval_maximum_air_c=round(float(c['maxt'][L].max()) - 273.15, 2),
                 mean_air_c_range=[round(float(tas[L].min()), 2), round(float(tas[L].max()), 2)],
                 lunar_day_swing_k_median=round(float(np.median(swing[L])), 2),
                 upland_share=round(float(wl[L & upland].sum()), 4))
    hgt = c['height']
    slope = np.polyfit(hgt[L] / 1e3, tas[L], 1)[0]
    top = L & (hgt >= np.quantile(hgt[L], 0.9))
    bottom = L & (hgt <= np.quantile(hgt[L], 0.1))
    frost.update(air_lapse_k_per_km=round(float(slope), 3),
                 highest_tenth_minus_lowest_tenth_air_k=round(float(tas[top].mean() - tas[bottom].mean()), 2),
                 highest_tenth_height_km=round(float(hgt[top].mean()) / 1e3, 2),
                 lowest_tenth_height_km=round(float(hgt[bottom].mean()) / 1e3, 2))
    dry = L & (c['ai'] < 0.65)
    robust_dry = dry & (bands(c['ai'], ai_se)[1] < 0.65)
    soda = dict(closed_rainfed_lake_km2_in_dry_cells=round(float(geo['closed_lake_km2'][dry].sum())),
                closed_rainfed_lake_km2_in_robustly_dry_cells=round(float(geo['closed_lake_km2'][robust_dry].sum())),
                cells_with_closed_lakes_in_dry_cells=int((dry & (geo['closed_lake_km2'] > 0)).sum()),
                sea_level_lake_km2_in_dry_cells=round(float(geo['sea_level_lake_km2'][dry].sum())),
                closed_rainfed_lake_km2_on_land=round(float(geo['closed_lake_km2'][L].sum())),
                note=f'Closed rain-fed lakes are the non-spilling lakes of the drainage routed on these years\' climate '
                     f'({drainage}); the sea-level lakes are the atlas\'s small closed depressions below the seas\' '
                     'level. Where they stand in dry cells they are the soda lakes\' candidates.')
    coast = L & (geo['coast_km'] > 0)
    tclass = tide_class(geo['tide_m'])
    shore = dict(coast_km_on_land_cells=round(float(geo['coast_km'][L].sum())),
                 coast_km_on_dry_cells=round(float(geo['coast_km'][dry].sum())),
                 shore_belt_km2=[round(float(geo['coast_km'][L].sum()) * b) for b in SHORE_BELT_KM],
                 shore_belt_km2_dry_cells=[round(float(geo['coast_km'][dry].sum()) * b) for b in SHORE_BELT_KM],
                 coastal_cell_share_of_land=round(float(wl[coast].sum()), 4),
                 tide_classes_by_coast_km={n: round(float(geo['coast_km'][L & (tclass == j)].sum()))
                                           for j, (n, _, _) in enumerate(TIDE_CLASSES)},
                 note='Coast length counts the 4 px/deg land-sea pixel edges in each land cell, so it depends on that '
                      'resolution; the belts take the flat shores\' 0.6-2.5 km on every shore, an upper-end screen.')
    def col(a, digits=4, mask=L):
        a = np.asarray(a)[mask]
        if a.dtype == bool:
            return [bool(x) for x in a]
        if a.dtype.kind in 'iu':
            return [int(x) for x in a]
        return [None if not np.isfinite(x) else float(f'{x:.{digits}g}') for x in a.astype(float)]

    rows2, cols2 = np.meshgrid(np.arange(lat.size), np.arange(lon.size), indexing='ij')
    cells = dict(row=col(rows2), col=col(cols2), lat=col(lat2, 6), lon=col(lon2, 6), area_km2=col(area, 5), weight=col(wl, 6), height_m=col(c['height'], 5),
                 rain_mm_yr=col(c['rain'] * 365.0), demand_mm_yr=col(c['pet'] * 365.0), aridity=col(c['ai'], 5),
                 aridity_spans=[col(s, 5) for s in c['spans']], aridity_se=col(ai_se, 4),
                 aridity_store=[col(s['ai'], 5) for s in stores],
                 margin_down=col(np.minimum(down, 99.0), 4), margin_up=col(np.minimum(up, 99.0), 4),
                 moisture=col(rob['run']), moisture_gap=col(rob['gap']), moisture_band_low=col(rob['combined_low']),
                 moisture_band_high=col(rob['combined_high']), moisture_years_low=col(rob['years_low']),
                 moisture_years_high=col(rob['years_high']), moisture_store=[col(s) for s in rob['store']],
                 moisture_at_aridity_ratio=col(at_aridity_ratio), unep_class=col(unep_run),
                 robust_gap=col(rob['robust_gap']), robust_years=col(rob['robust_years']), robust=col(rob['robust']),
                 robust_store=col(rob['robust_store']), robust_all=col(rob['robust_all']),
                 air_c=col(tas), coldest_3day_air_c=col(c['tas_min'] - 273.15), warmest_3day_air_c=col(c['tas_max'] - 273.15),
                 lunar_day_swing_k=col(swing, 3), upland=col(upland), night_rain_share=col(night_share, 3),
                 night_rain_flag=col(night_share >= NIGHT_FLAG), beside_sea=col(beside_sea), class_kept_in_resamples=col(keep_prob, 3), night_rh_lowest_level=col(c['night_rh'], 3),
                 soil_water_mm=col(c['mean']['mrso'] * 1e3, 4), darkest_lux=col(lux, 4), earth_elevation_deg=col(elev, 4),
                 night_regime=col(regime), plant_dark_hours=col(hours, 4),
                 sea_fraction=col(geo['sea_fraction'], 3), sea_level_lake_fraction=col(geo['sea_level_lake_fraction'], 3),
                 rainfed_lake_fraction=col(geo['rainfed_lake_fraction'], 3),
                 closed_lake_fraction=col(geo['closed_lake_fraction'], 3), coast_km=col(geo['coast_km'], 4),
                 tide_m=col(geo['tide_m'], 3), tide_class=col(tclass))
    lakes = [f'{drainage}/results/drainage.json', f'{drainage}/products/drainage_28pct_16ppd.npz']
    inputs = {path: sha(ROOT / path) for path in [*GEOGRAPHY.values(), *lakes, LIVING_LIGHT, PLANT, *([orientation] if frame else [])]}
    for s in (c, *stores):
        for name, ident in s['files'].items():
            inputs[f'climate/gcm/runs/{name}'] = ident
        inputs[f"climate/gcm/runs/{s['run']}/inputs/moon_topography.sra"] = sha(RUNS / s['run'] / 'inputs' / 'moon_topography.sra')
    files = [HERE / 'moisture_map.py', ROOT / 'climate' / 'crm' / 'gcm_site.py', ROOT / 'climate' / 'gcm' / 'year_files.py',
             ROOT / 'climate' / 'gcm' / 'turned.py']
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='biosphere', model='biosphere/ecology/moisture_map.py',
                      files={str(f.relative_to(ROOT)): sha(f) for f in files},
                      constants=constants_used([str(f.relative_to(ROOT)) for f in files[:1]]),
                      run=run, years=list(years), store_runs=list(store_runs), store_years=list(store_years),
                      frame=(dict(orientation=orientation, sha256=frame['sha256'], tilt_deg=frame['tilt_deg']) if frame else None),
                      drainage=drainage, inputs=inputs, input_schemas={**SCHEMAS, lakes[0]: 'terluna.geography.drainage/1'}),
        evidence=('A zoning of the design climate by its water balance, not a model of vegetation: each T21 land cell '
                  '(about 170 km across at the equator) of the turned run at the target takes the moisture class its '
                  'aridity gives it, with its warmth, night light and settings beside it. The aridity is the GCM\'s rain '
                  'over a Priestley-Taylor demand without wind or humidity, over the years named; the classes are UNEP\'s, '
                  'the humid class divided at aridity 1, each with the broad biome the Millennium Ecosystem Assessment '
                  'gives the class on Earth beside it as a reference. Which organisms each class holds is a design '
                  'question. The cells are the turned Moon\'s; geography\'s pixels are counted in the cell that holds '
                  'them on the turned Moon through geology\'s orientation product. The gap\'s scale of 0.83 is a '
                  'sensitivity of the classes from three CM1 boxes forced by the design years of 1 November, never a '
                  'bias-corrected climate. The lakes come from the drainage routed on the same years\' climate; the '
                  'night light from the living-light grid in the mean geometry under a clear sky. Provisional while the '
                  'GCM-CM1 gap stands.'),
        reading_rule=('cells holds one entry per GCM land cell (land_binary_mask > 0.5) of the run, in columns: lat and '
                      'lon (degrees, east) are the T21 cell centres; weight is the cell\'s share of the land by cos(lat), '
                      'as the design readings weight it; area_km2 the cell\'s box area. aridity is rain over demand '
                      '(mm/yr, 365-day years) from the means over the years; aridity_spans the two-year spans; '
                      'aridity_se the mean\'s standard error from them; aridity_store the store runs\', in '
                      'store_runs order. moisture indexes classes (0 hyper-arid to 5 wet) as run; '
                      'moisture_gap at the rain times gap.scale; moisture_band_low and _high at the combined '
                      'band\'s ends, gap.scale times (aridity - 2 se) and aridity + 2 se. robust is true where the '
                      'class holds across that band, robust_gap and robust_years for each source alone, '
                      'robust_store on all three store runs, robust_all both. margin_down and margin_up are how far '
                      'the rain must fall or rise, as a fraction, to leave the class (99 where no bound). '
                      'class_kept_in_resamples is the share of the bootstrap\'s resamples of the years that keep the '
                      'cell\'s class; beside_sea marks land cells with a GCM water cell among their eight neighbours. night_regime '
                      'indexes night_regimes; tide_class indexes tide_classes (-1 no coast); shares are of the land.'),
        units=dict(rain='mm a year', demand='mm a year', air='deg C at 2 m', light='lux on level ground',
                   coast='km', areas='km2'),
        classes=[dict(c_, index=i, aridity=[c_['aridity'][0], c_['aridity'][1] if math.isfinite(c_['aridity'][1]) else None])
                     for i, c_ in enumerate(CLASSES)],
        night_regimes=list(NIGHT_REGIMES), night_regime_lux=dict(dark=LUX_DARK, bright=LUX_BRIGHT),
        tide_classes=[dict(name=n, range_m=[lo, hi if math.isfinite(hi) else None]) for n, lo, hi in TIDE_CLASSES],
        unep_classes=[dict(name=n, aridity=[lo, hi if math.isfinite(hi) else None]) for n, lo, hi in UNEP],
        gap=dict(scale=GAP_SCALE, box_rain_cm1_gcm_mm_day=GAP_RAIN,
                 box_ratios={k: round(a / b, 4) for k, (a, b) in GAP_RAIN.items()},
                 land_box_aridity_cm1_gcm=GAP_ARIDITY, aridity_ratio=round(aridity_ratio, 4),
                 night_share_cm1_gcm=GAP_NIGHT_SHARE, night_flag_share=NIGHT_FLAG, night_low_ratio=NIGHT_LOW_RATIO,
                 note='The rain scaled uniformly toward CM1\'s ratio at the three gap boxes is a sensitivity of the '
                      'classes, never a bias-corrected climate; three boxes do not give the gap\'s pattern.'),
        grid=dict(lat=[float(x) for x in lat], lon=[float(x) for x in lon],
                  lat_edges=[float(x) for x in cell_edges(lat)],
                  note='T21 Gaussian rows north to south; row edges are midpoints between rows; cells index row and col'),
        land=dict(cells=int(L.sum()), area_km2=round(float(area[L].sum())), share_of_moon=round(float(area[L].sum() / area.sum()), 4),
                  resolution='T21: 64 by 32 cells, 5.625 deg of longitude (170.6 km at the equator) by about 5.5 deg of '
                             'latitude (about 168 km)'),
        shares=dict(by_class=by_class, unep_as_run=unep_shares, unep_by_rain_scale=sens,
                    class_at_aridity_ratio=share_rows(at_aridity_ratio),
                    robust=round(float(wl[L & rob['robust']].sum()), 4),
                    robust_gap=round(float(wl[L & rob['robust_gap']].sum()), 4),
                    robust_years=round(float(wl[L & rob['robust_years']].sum()), 4),
                    robust_store=round(float(wl[L & rob['robust_store']].sum()), 4),
                    robust_all=round(float(wl[L & rob['robust_all']].sum()), 4),
                    land_changing_at_aridity_ratio=round(float(wl[L & (at_aridity_ratio != rob['run'])].sum()), 4),
                    store_runs={s: dict(land_changing=round(float(wl[L & (k != rob['run'])].sum()), 4),
                                        by_class=share_rows(k),
                                        difference_over_sampling_sd={n: round(float((v - by_class[cm['key']]['share']) / max(math.sqrt(2.0) * by_class[cm['key']]['share_sd_years'], 1e-9)), 2)
                                                                     for (n, v), cm in zip(share_rows(k).items(), CLASSES)})
                                for s, k in zip(store_runs, rob['store'])},
                    transitions_to_band_low=transitions,
                    years_spans_changing=round(float(wl[L & np.any(np.array([moisture(a) for a in c['spans']]) != rob['run'], axis=0)].sum()), 4),
                    night_regimes=dict(zip(NIGHT_REGIMES, (round(v, 4) for v in shares(regime[L], w[L], len(NIGHT_REGIMES))))),
                    night_rain_flagged=round(float(wl[flagged_night].sum()), 4),
                    night_rain_flagged_beside_sea=round(float(wl[flagged_night & beside_sea].sum()), 4),
                    night_rain_flagged_poleward_of_60=round(float(wl[flagged_night & (np.abs(lat2) >= 60.0)].sum()), 4),
                    beside_sea=round(float(wl[beside_sea].sum()), 4),
                    robust_at_one_se=round(float(wl[L & one_se].sum()), 4),
                    kept_in_95_percent_of_resamples=round(float(wl[L & (keep_prob >= 0.95)].sum()), 4),
                    kept_in_half_of_resamples=round(float(wl[L & (keep_prob >= 0.5)].sum()), 4),
                    bootstrap=dict(resamples=BOOTSTRAP, seed=SEED, note='years drawn with replacement from the run\'s '
                                   'years, each resample composited as the full record; shares\' spread is the '
                                   'sampling spread of a mean over this many years'),
                    class_by_night_regime={f'{names[i // len(NIGHT_REGIMES)]} | {NIGHT_REGIMES[i % len(NIGHT_REGIMES)]}': round(v, 4)
                                               for i, v in enumerate(combo_shares) if v > 0}),
        mix=dict(simpson_diversity=round(simpson, 4), effective_classes=round(1.0 / (1.0 - simpson), 3),
                 classes_over_5_percent=int(sum(by_class[k['key']]['share'] >= 0.05 for k in CLASSES)),
                 class_and_night_pairs_over_1_percent=int(sum(v >= 0.01 for v in combo_shares)),
                 drylands_share=round(float(wl[dry].sum()), 4), earth_drylands=EARTH_DRYLANDS),
        warmth=frost,
        light=dict(plant_dark_hours_by_band=hour_bands, source=LIVING_LIGHT),
        settings=dict(shores=shore, soda_lakes=soda, soil_zones=SOIL_ZONES, storm_water_surfaces=STORM_WATER_SURFACES,
                      night_humidity=dict(lowest_level_night_rh_median_dry=round(float(np.median(c['night_rh'][dry])), 3),
                                          lowest_level_night_rh_median_humid=round(float(np.median(c['night_rh'][L & ~dry])), 3),
                                          land_share_night_rh_over_0_9=round(float(wl[L & (c['night_rh'] >= 0.9)].sum()), 4),
                                          note='Relative humidity of the 3-day mean state at the lowest level, nights '
                                               '(|hour angle| > 90 deg) averaged; dew needs about 0.9 (coast-water study).')),
        checks=checks,
        cells=cells)
    return product


def compare_windows(product, compare, run=RUN, decades=DECADES):
    """How far the classes and the lakes move between windows of one run: the main window against a window mapped
    whole beside it, each with its own lakes, and the decades classed alone, against each other and the main window.
    The decades are independent samples of one climate, so the land that changes class between two of them measures how
    far a ten-year map moves by chance; the run's slow drift adds a little to it."""
    c, k = product['cells'], compare['cells']
    if (c['row'], c['col']) != (k['row'], k['col']):
        raise ValueError('the windows must share the run\'s land cells')
    rows, cols, w = np.array(c['row']), np.array(c['col']), np.array(c['weight'])
    main_cls, cmp_cls = np.array(c['moisture']), np.array(k['moisture'])
    dec = {f'{a}-{b}': moisture(read_run(run, (a, b))['ai'])[rows, cols] for a, b in decades}
    agree = lambda x, y: round(float(w[x == y].sum()), 4)                                           # noqa: E731
    class_shares = lambda x: {cl['key']: round(float(w[x == i].sum()), 4) for i, cl in enumerate(CLASSES)}  # noqa: E731
    names = list(dec)
    pairs = {f'{p} | {q}': agree(dec[p], dec[q]) for i, p in enumerate(names) for q in names[i + 1:]}
    area = np.array(c['area_km2'])
    lake_km2 = lambda cells, key: round(float((np.array(cells[key], float) * area).sum()))      # noqa: E731
    c['moisture_compare'] = k['moisture']
    c['moisture_decades'] = [[int(x) for x in dec[n]] for n in names]
    return dict(
        main_years=product['producer']['years'], compare_years=compare['producer']['years'], decades=names,
        main_against_compare=dict(
            land_keeping_class=agree(main_cls, cmp_cls), shares_main=class_shares(main_cls), shares_compare=class_shares(cmp_cls),
            drylands_main=product['mix']['drylands_share'], drylands_compare=compare['mix']['drylands_share'],
            rainfed_lake_km2_on_land=[lake_km2(c, 'rainfed_lake_fraction'), lake_km2(k, 'rainfed_lake_fraction')],
            closed_lake_km2_on_land=[lake_km2(c, 'closed_lake_fraction'), lake_km2(k, 'closed_lake_fraction')],
            rainfed_lake_fraction_mean_abs_difference=round(float(np.mean(np.abs(np.array(c['rainfed_lake_fraction'], float)
                                                                                 - np.array(k['rainfed_lake_fraction'], float)))), 5),
            compare_checks=compare['checks'], compare_drainage=compare['producer']['drainage']),
        decades_alone=dict(shares={n: class_shares(x) for n, x in dec.items()}, land_keeping_class_between=pairs,
                           between_mean=round(float(np.mean(list(pairs.values()))), 4),
                           between_range=[min(pairs.values()), max(pairs.values())],
                           against_main={n: agree(x, main_cls) for n, x in dec.items()}),
        note=('Both windows read the same run\'s land cells. The compare window overlaps the main one, so their '
              'agreement exceeds that of independent windows; the decades classed alone are independent samples, and '
              'the land that changes class between two of them is the chance movement of a ten-year map.'),
        reading_rule=('cells.moisture_compare holds each land cell\'s class in the compare window and cells.moisture_decades '
                      'its class in each decade, in decades order; shares and agreements are of the land by weight.'))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--run', default=RUN)
    ap.add_argument('--years', default=f'{YEARS[0]}-{YEARS[1]}')
    ap.add_argument('--compare-years', default=f'{COMPARE_YEARS[0]}-{COMPARE_YEARS[1]}',
                    help="a window mapped whole beside the main one, with the decades classed alone; 'none' for neither")
    ap.add_argument('--orientation', default=ORIENTATION, help="the run's frame; 'none' for a run on today's ground")
    ap.add_argument('--store-runs', nargs='*', default=list(STORE_RUNS))
    ap.add_argument('--store-years', default=f'{STORE_YEARS[0]}-{STORE_YEARS[1]}')
    ap.add_argument('--out', default=str(OUT))
    ap.add_argument('--cache', default=str(DRIVE / 'moisture_map' / 'cache'),
                    help="a folder for each year file's sums, so the windows and a rebuild need not reread the run")
    a = ap.parse_args(argv)
    global CACHE
    CACHE = a.cache
    if CACHE:
        Path(CACHE).mkdir(parents=True, exist_ok=True)
    years = tuple(int(v) for v in a.years.split('-'))
    store_years = tuple(int(v) for v in a.store_years.split('-'))
    orientation = None if a.orientation == 'none' else a.orientation
    product = build(a.run, years, tuple(a.store_runs), store_years, orientation=orientation)
    if a.compare_years != 'none':
        cy = tuple(int(v) for v in a.compare_years.split('-'))
        product['windows'] = compare_windows(product, build(a.run, cy, (), store_years, orientation=orientation), a.run)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(product, indent=1) + '\n')
    s = product['shares']
    print(f"{a.run} years {a.years}: {product['land']['cells']} land cells, robust {s['robust']:.3f} (gap {s['robust_gap']:.3f}, "
          f"years {s['robust_years']:.3f}, all {s['robust_all']:.3f})")
    for key, row in s['by_class'].items():
        print(f"  {key}: {row['share']:.3f} as run, {row['robust_share']:.3f} robust, {row['share_at_gap_scale']:.3f} at the gap's scale")
    print('  checks', json.dumps(product['checks']))
    print('  warmth', json.dumps(product['warmth']))
    print('  mix', json.dumps({k: v for k, v in product['mix'].items() if k != 'earth_drylands'}))
    if 'windows' in product:
        wd = product['windows']
        print('  windows', json.dumps({k: wd['main_against_compare'][k] for k in ('land_keeping_class', 'drylands_main', 'drylands_compare',
              'rainfed_lake_km2_on_land')}), json.dumps({k: wd['decades_alone'][k] for k in ('between_mean', 'between_range', 'against_main')}))
    print(f'wrote {out} sha256 {sha(out)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
