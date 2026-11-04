"""Where the Open Moon's fresh waters run coldest: highland lakes and shallow ponds through the lunar night.

A first estimate on the design years, A28R_301k_design_l15_snow_s934 years 10-19 (1,200 3-day means at T21). The
seas' deep water settles at the temperature of the coldest water their surface makes, 25-30 C at 301 K
(research/studies/sea_oxygen), so two places could still hold cooler water.

Highland lakes. Every rain-fed lake of geography's drainage product (16 px/deg) takes its T21 cell's air at 2 m
(the mean, the coldest 3-day mean and the coldest minimum within a mean) carried from the cell's mean ground (the
design run's topography) to the lake's level along the design years' own lapse rate over land, the land columns'
mean warming going down between 0.9 and 5 km above the ground (hypsometric heights on the sigma levels, virtual
temperature). A deep lake's floor water settles near its surface's coldest spell, which the 30 m lake below tests.

Ponds. The project's lake model (biosphere/ecology/lake.py) on pans of 0.5, 1 and 2 m, open to the 10 m wind and
sheltered (the wind cut to SHELTER of it, floored at the model's 0.5 m/s), through the ten years at three sites: the
coldest lake's, its air carried to the lake's level; the land cell with the coldest 3-day mean; and the land cell with
the median mean air. A 30 m pan at the coldest lake's site shows how cold a lake's floor gets. Forcing as the waters'
oxygen study builds it from the GCM's means (research/studies/lake_oxygen/run.py, windows): sunlight at the ground
rss / (1 - alb), sky infrared rls + sigma ts^4 over land, the air and the lowest level's humidity, the surface
pressure; here the 10 m wind with the model's 10 m coefficients, and every field linear in time between the means'
centres, which smears sunrise and sunset over a 3-day mean.

    /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.cold_waters.run   # minutes

writes results/cold_waters.json (schema terluna.research.cold-waters/1) with the sha256 of every input read.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = Path('/media/projectspace/terluna-research/gcm-runs')
RUN, YEARS = 'A28R_301k_design_l15_snow_s934', (10, 19)
DRAINAGE = ROOT / 'geography' / 'products' / 'drainage_28pct_16ppd.npz'
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
LAKE_OXYGEN = ROOT / 'research' / 'studies' / 'lake_oxygen' / 'results' / 'lake_oxygen.json'
OUT = HERE / 'results' / 'cold_waters.json'
SCHEMA = 'terluna.research.cold-waters/1'
PPD = 16
STEPS_PER_MEAN = 72                     # lake-model steps in each 3-day mean (3,575 s), as the waters' oxygen study
LAPSE_LEVELS_KM = (0.9, 5.0)            # the band whose mean warming going down sets the lapse rate
THRESHOLDS_C = (26.0, 24.0, 22.0, 20.0, 18.0, 15.0)
POND_DEPTHS_M = (0.5, 1.0, 2.0)
POND_LAYER_M = 0.1
LAKE_DEPTH_M, LAKE_LAYER_M = 30.0, 0.5
SHELTER = 0.3                           # a sheltered pond's share of the 10 m wind (banks, reeds, trees)
SPIN_MEANS = dict(pond=20, lake=240)    # 3-day means left out before reading: two lunar days; two years for the lake
STEFAN_BOLTZMANN = 5.670374419e-8


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def year_files():
    return [RUNS / RUN / 'model' / f'MOST.{y:05d}.nc' for y in range(YEARS[0], YEARS[1] + 1)]


def read_design_years():
    """The ten years' 3-day means of the fields both estimates need, and each land column's warming going down between
    the lapse band's heights (K/m), one value per mean."""
    import netCDF4
    keys = ('tas', 'mint', 'ts', 'rss', 'rls', 'alb', 'ps', 'uas', 'vas', 'lsm')
    out = {k: [] for k in keys + ('hus_low', 'lapse')}
    for path in year_files():
        with netCDF4.Dataset(path) as nc:
            lat, lon, sigma = (np.array(nc[k][:], dtype=float) for k in ('lat', 'lon', 'lev'))
            for k in keys:
                out[k].append(np.array(nc[k][:], dtype=np.float32))
            ta = np.array(nc['ta'][:], dtype=float)
            hus = np.array(nc['hus'][:], dtype=float)
            out['hus_low'].append(hus[:, -1].astype(np.float32))
            out['lapse'].append(lapse_rates(ta, hus, sigma).astype(np.float32))
    data = {k: np.concatenate(v) for k, v in out.items()}
    data.update(lat=lat, lon=lon)
    return data


def lapse_rates(ta, hus, sigma, gas_constant=None, gravity=MOON_SURFACE_GRAVITY, levels_km=LAPSE_LEVELS_KM):
    """Each column's warming going down (K/m) between two heights above its ground: the sigma levels' heights by the
    hypsometric equation on the virtual temperature (the lowest level from the ground, sigma = 1)."""
    rd = gas_constant or air_constants()[0]
    tv = ta * (1.0 + 0.608 * hus)                                       # (time, level, lat, lon), top level first
    log_s = np.log(sigma)
    z = np.empty_like(tv)
    z[:, -1] = rd / gravity * tv[:, -1] * (0.0 - log_s[-1])
    for k in range(len(sigma) - 2, -1, -1):
        z[:, k] = z[:, k + 1] + rd / gravity * 0.5 * (tv[:, k] + tv[:, k + 1]) * (log_s[k + 1] - log_s[k])
    lo, hi = (1000.0 * x for x in levels_km)
    # heights rise toward index 0; interpolate each column's temperature at the two heights
    zz = z[:, ::-1]
    tt = ta[:, ::-1]
    def at(height):
        above = np.argmax(zz >= height, axis=1)                         # first level at or above the height
        above = np.clip(above, 1, zz.shape[1] - 1)
        z1 = np.take_along_axis(zz, above[:, None] - 1, 1)[:, 0]
        z2 = np.take_along_axis(zz, above[:, None], 1)[:, 0]
        t1 = np.take_along_axis(tt, above[:, None] - 1, 1)[:, 0]
        t2 = np.take_along_axis(tt, above[:, None], 1)[:, 0]
        return t1 + (t2 - t1) * (height - z1) / (z2 - z1)
    return (at(lo) - at(hi)) / (hi - lo)


def air_constants():
    """The design air's gas constant and heat capacity, as the waters' oxygen study holds them."""
    s = json.loads(LAKE_OXYGEN.read_text())['setting']
    return float(s['air_gas_constant']), float(s['air_cp'])


def topography_m(lat, lon):
    """Each T21 cell's mean ground above sea level (m): the design run's surface geopotential over g."""
    path = RUNS / RUN / 'inputs' / 'moon_topography.sra'
    phi = np.loadtxt(path, skiprows=1).ravel().reshape(len(lat), len(lon))
    return phi / MOON_SURFACE_GRAVITY, path


# ---------------------------------------------------------------- highland lakes

def lake_bodies():
    """Geography's rain-fed lakes as connected bodies (eight neighbours): area (km2), level above sea level (m) and
    centre (degrees, longitude 0-360)."""
    from scipy import ndimage
    d = np.load(DRAINAGE, allow_pickle=False)
    cells = d['lake_cell'].astype(np.int64)
    level = d['lake_level_m'].astype(float)
    sea_level = float(json.loads(ATLAS.read_text())['sea_level_m'])
    rows, cols = 180 * PPD, 360 * PPD
    mask = np.zeros(rows * cols, bool)
    mask[cells] = True
    labels, n = ndimage.label(mask.reshape(rows, cols), structure=np.ones((3, 3), bool))
    body = labels.ravel()[cells]
    r, c = cells // cols, cells % cols
    lat = 90.0 - (r + 0.5) / PPD
    lon = (c + 0.5) / PPD
    half = 0.5 / PPD
    area = (MOON_RADIUS / 1e3) ** 2 * np.radians(1.0 / PPD) * (np.sin(np.radians(lat + half))
                                                                - np.sin(np.radians(lat - half)))
    a = np.bincount(body, area, n + 1)
    top = np.full(n + 1, -np.inf)
    np.maximum.at(top, body, level)
    b_lat = np.bincount(body, area * lat, n + 1) / np.maximum(a, 1e-12)
    x = np.bincount(body, area * np.cos(np.radians(lon)), n + 1)
    y = np.bincount(body, area * np.sin(np.radians(lon)), n + 1)
    keep = np.arange(1, n + 1)
    return dict(area_km2=a[keep], height_m=top[keep] - sea_level, lat=b_lat[keep],
                lon=np.degrees(np.arctan2(y[keep], x[keep])) % 360.0, sea_level_m=sea_level)


def nearest_cell(lat, lon, grid_lat, grid_lon):
    j = np.abs(np.asarray(lat)[:, None] - grid_lat[None, :]).argmin(1)
    dl = np.abs((np.asarray(lon)[:, None] - grid_lon[None, :] + 180.0) % 360.0 - 180.0)
    return j, dl.argmin(1)


def carried(air_k, cell_height_m, lake_height_m, lapse):
    """Air at a lake's level: the cell's air moved by the lapse rate (K/m) over the lake's height below the cell's
    mean ground (colder above it)."""
    return air_k + lapse * (cell_height_m - lake_height_m)


def highland_lakes(data, heights, lapse):
    lakes = lake_bodies()
    j, i = nearest_cell(lakes['lat'], lakes['lon'], data['lat'], data['lon'])
    tas = data['tas'].astype(float)
    mean_k, coldest_k = tas.mean(0)[j, i], tas.min(0)[j, i]
    minimum_k = data['mint'].astype(float).min(0)[j, i]
    cell_h = heights[j, i]
    out = {}
    for name, rate in lapse.items():
        out[name] = dict(mean_c=carried(mean_k, cell_h, lakes['height_m'], rate) - 273.15,
                         coldest_c=carried(coldest_k, cell_h, lakes['height_m'], rate) - 273.15,
                         minimum_c=carried(minimum_k, cell_h, lakes['height_m'], rate) - 273.15)
    return lakes, dict(j=j, i=i, cell_height_m=cell_h), out


def summarise_lakes(lakes, where, air):
    a, h = lakes['area_km2'], lakes['height_m']
    order = np.argsort(h)
    cum = np.cumsum(a[order]) / a.sum()
    def height_at(q):
        return float(h[order][np.searchsorted(cum, q)])
    central = air['central']
    below = {}
    for t in THRESHOLDS_C:
        sel = central['coldest_c'] < t
        sel_m = central['mean_c'] < t
        below[f'{t:g}'] = dict(lakes_coldest_spell_below=int(sel.sum()), area_km2_coldest_spell_below=float(a[sel].sum()),
                               lakes_mean_below=int(sel_m.sum()), area_km2_mean_below=float(a[sel_m].sum()),
                               lakes_coldest_spell_below_any_lapse=int(np.any([air[k]['coldest_c'] < t
                                                                               for k in air], axis=0).sum()))
    coldest = np.argsort(central['coldest_c'])[:20]
    return dict(
        lakes=int(a.size), area_km2=float(a.sum()),
        height_above_sea_m=dict(area_weighted_median=height_at(0.5), area_weighted_p90=height_at(0.9),
                                area_weighted_p99=height_at(0.99), highest=float(h.max()),
                                lakes_above_2km=int((h > 2000).sum()), lakes_above_4km=int((h > 4000).sum()),
                                area_km2_above_2km=float(a[h > 2000].sum())),
        air_at_the_lakes_c=dict(mean_area_weighted=float(np.sum(a * central['mean_c']) / a.sum()),
                                coldest_spell_area_weighted=float(np.sum(a * central['coldest_c']) / a.sum()),
                                coldest_spell_lowest=float(central['coldest_c'].min()),
                                minimum_lowest=float(central['minimum_c'].min())),
        below_threshold=below,
        coldest_lakes=[dict(lat=round(float(lakes['lat'][k]), 2), lon=round(float(lakes['lon'][k]), 2),
                            height_above_sea_m=round(float(h[k])), cell_ground_above_sea_m=round(float(where['cell_height_m'][k])),
                            area_km2=round(float(a[k]), 2),
                            mean_air_c=round(float(central['mean_c'][k]), 2),
                            coldest_spell_c=round(float(central['coldest_c'][k]), 2),
                            coldest_minimum_c=round(float(central['minimum_c'][k]), 2),
                            coldest_spell_c_lapse_bracket=[round(float(air[x]['coldest_c'][k]), 2) for x in ('low', 'high')])
                       for k in coldest])


# ---------------------------------------------------------------- ponds

def site_forcing(data, j, i, lift_m=0.0, lapse=0.0, shelter=1.0):
    """A cell's 3-day means as the lake model's step forcing; lift_m raises the water above the cell's mean ground,
    carrying the air and sky infrared up the lapse rate and the pressure up its scale height."""
    rd, _ = air_constants()
    g = MOON_SURFACE_GRAVITY
    tas = data['tas'][:, j, i].astype(float)
    alb = data['alb'][:, j, i].astype(float)
    sw = np.where(alb < 0.999, data['rss'][:, j, i] / np.maximum(1.0 - alb, 1e-3), 0.0)
    lw = data['rls'][:, j, i] + STEFAN_BOLTZMANN * data['ts'][:, j, i].astype(float) ** 4      # land: emissivity 1
    q_open = data['hus_low'][:, j, i].astype(float)
    q = q_open.copy()
    q[:-1] = 0.5 * (q_open[:-1] + q_open[1:])        # the output's humidity is the value opening each window
    p = data['ps'][:, j, i].astype(float) * 100.0
    wind = np.hypot(data['uas'][:, j, i], data['vas'][:, j, i]).astype(float) * shelter
    air = tas - lapse * lift_m
    if lift_m:
        lw = lw * (air / tas) ** 4
        p = p * np.exp(-g * lift_m / (rd * 0.5 * (air + tas)))
    qs = np.array([L.saturation_humidity(a, pp) for a, pp in zip(air, p)])
    q = np.minimum(q, qs)
    window = 1430 / 10 * 30.0 * 60.0                  # a 3-day mean: a tenth of the lunar day of 1,430 half-hour steps
    dt = window / STEPS_PER_MEAN
    centres = (np.arange(tas.size) + 0.5) * window
    t = (np.arange(tas.size * STEPS_PER_MEAN) + 0.5) * dt
    f = {k: np.interp(t, centres, v) for k, v in dict(sw=sw, lw=lw, air=air, q=q, u=wind, p=p).items()}
    forcing = L.Forcing(dt=dt, sw_down=np.maximum(f['sw'], 0.0), lw_down=f['lw'], air_k=f['air'], q_air=f['q'],
                        wind10=f['u'], pressure=f['p'], o2_fraction=0.2)
    return forcing, f['air'], window


def run_water(forcing, depth, layer, start_c):
    rd, cp = air_constants()
    basin = L.Basin.pan(depth, layer, 1e4)
    scenario = L.Scenario('temperature', gpp=1.6, attenuation=1.7 / 3.0, sediment=0.3, deep_water=0.3)
    lake = L.Lake(basin, L.Physics(deep_mixing='wind'), [scenario], L.Air(rd, cp, MOON_SURFACE_GRAVITY))
    state = L.initial_state(basin, [scenario], start_c)
    return lake.run(forcing, state)


def read_water(out, air_k, dt, window, spin_means, deep=False):
    """The water's surface (or floor) temperature after the spin-up against the air: mean, lowest, the lowest a
    3-day mean reaches, and the hours a lunar day spent below each threshold."""
    key = 'bottom_t' if deep else 'surface_t'
    x = np.asarray(out[key], dtype=float)
    start = spin_means * STEPS_PER_MEAN
    x, a = x[start:], np.asarray(air_k, dtype=float)[start:] - 273.15
    n = (x.size // STEPS_PER_MEAN) * STEPS_PER_MEAN
    means = x[:n].reshape(-1, STEPS_PER_MEAN).mean(1)
    lunar_days = x.size * dt / (10 * window)
    return dict(mean_c=float(x.mean()), lowest_c=float(x.min()), lowest_3day_mean_c=float(means.min()),
                air_mean_c=float(a.mean()), air_lowest_c=float(a.min()),
                hours_per_lunar_day_below={f'{t:g}': float((x < t).sum() * dt / 3600.0 / lunar_days) for t in THRESHOLDS_C})


def ponds(data, heights, lakes, where, air, lapse_central):
    land = data['lsm'][0] > 0.5
    tas = data['tas'].astype(float)
    coldest_cell = np.unravel_index(np.where(land, tas.min(0), np.inf).argmin(), land.shape)
    means = np.where(land, tas.mean(0), np.nan)
    median = np.nanmedian(means)
    typical_cell = np.unravel_index(np.nanargmin(np.abs(means - median)), land.shape)
    k = int(np.argmin(air['central']['coldest_c']))
    lake_cell = (int(where['j'][k]), int(where['i'][k]))
    lift = float(lakes['height_m'][k] - where['cell_height_m'][k])
    sites = {'the coldest lake': (lake_cell, lift), 'the coldest land cell': (coldest_cell, 0.0),
             'a typical land cell': (typical_cell, 0.0)}
    out = {}
    for name, ((j, i), lift_m) in sites.items():
        site = dict(lat=float(data['lat'][j]), lon=float(data['lon'][i]), cell_ground_above_sea_m=float(heights[j, i]),
                    water_above_cell_ground_m=lift_m, cases={})
        for shelter_name, shelter in (('open', 1.0), ('sheltered', SHELTER)):
            forcing, air_k, window = site_forcing(data, j, i, lift_m, lapse_central, shelter)
            start = float(np.mean(air_k)) - 273.15
            for depth in POND_DEPTHS_M:
                res = run_water(forcing, depth, POND_LAYER_M, start)
                site['cases'][f'{depth:g} m pond, {shelter_name}'] = read_water(res, air_k, forcing.dt, window,
                                                                               SPIN_MEANS['pond'])
            if name == 'the coldest lake' and shelter_name == 'open':
                res = run_water(forcing, LAKE_DEPTH_M, LAKE_LAYER_M, start)
                site['cases'][f'{LAKE_DEPTH_M:g} m lake, surface'] = read_water(res, air_k, forcing.dt, window,
                                                                                 SPIN_MEANS['lake'])
                site['cases'][f'{LAKE_DEPTH_M:g} m lake, floor'] = read_water(res, air_k, forcing.dt, window,
                                                                               SPIN_MEANS['lake'], deep=True)
        out[name] = site
    return out


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    started = time.time()
    data = read_design_years()
    heights, topo = topography_m(data['lat'], data['lon'])
    land = data['lsm'][0] > 0.5
    lapse_cells = data['lapse'].astype(float).mean(0)[land]
    lapse = dict(central=float(np.mean(data['lapse'][:, land])), low=float(np.percentile(lapse_cells, 10)),
                 high=float(np.percentile(lapse_cells, 90)))
    lakes, where, air = highland_lakes(data, heights, lapse)
    record = dict(
        schema=SCHEMA, producer='research/studies/cold_waters/run.py',
        evidence=('A first estimate from model results on the design years (A28R_301k_design_l15_snow_s934, years '
                  '10-19), read before the GCM and CM1 agree: the lapse rate and the night air are among what that gap '
                  'could move. The lakes take their cells\' air carried along the lapse rate; the ponds run the '
                  'project\'s lake model on the cells\' means made linear in time.'),
        reading_rule=('Air temperatures are 2 m air at each lake\'s own level; coldest spell is the coldest 3-day mean '
                      'over ten years, minimum the coldest minimum within a mean. Lapse rates are K per km of warming '
                      'going down between 0.9 and 5 km above the ground, over land: central the mean over every land '
                      'column and mean, low and high the 10th and 90th percentiles of the land cells\' means.'),
        lapse_k_per_km={k: round(v * 1000.0, 3) for k, v in lapse.items()},
        highland_lakes=summarise_lakes(lakes, where, air),
        ponds=ponds(data, heights, lakes, where, air, lapse['central']),
        settings=dict(pond_depths_m=POND_DEPTHS_M, pond_layer_m=POND_LAYER_M, lake_depth_m=LAKE_DEPTH_M,
                      shelter=SHELTER, spin_means=SPIN_MEANS, steps_per_mean=STEPS_PER_MEAN,
                      thresholds_c=THRESHOLDS_C),
        inputs={str(p.relative_to(RUNS)): sha256(p) for p in year_files() + [topo]}
        | {str(p.relative_to(ROOT)): sha256(p) for p in (DRAINAGE, ATLAS, LAKE_OXYGEN)},
        seconds=round(time.time() - started, 1))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(OUT)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
