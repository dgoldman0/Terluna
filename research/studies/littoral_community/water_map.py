"""The waters' warmth across the Open Moon at the design climate: every cell's open water and its shallows through
the lunar month.

The forcing of five shore sites (forcing.py) is carried to every T21 cell of the corrected design climate on the
turned Moon, A28R_301k_design_l15_snow_s938_turned_refit_precession, years 10-19, by the same steps: the fitted Sun,
the run's clear sky read from its end-of-window snapshots as a function of the Sun's elevation, each window's cloud
factor from its all-sky sunlight over that clear sky averaged through the window, the sunlight shaped within the
window, and the sky infrared, the pressure and the lowest level's air, humidity, wind and height, at which the bulk
exchange is paired as PlaSim pairs it, linear in time between the windows' centres. One
Sun serves every cell: the noon time and the period are fitted on the row nearest the equator and carried along the
longitude, the declination taken as zero.

The water is the forcing's slab, stepped for every cell and case at once: a mixed layer of a given depth takes the
sunlight it absorbs, its infrared exchange with the sky, and the sensible and latent heat of the lake model's bulk
exchange (biosphere/ecology/lake.py, written here for arrays and tested against it), and trades a stated share of its
water each day with the open water beside it. The cases:

- still water 1, 2 and 4 m deep, enclosed: ponds, lagoons, sheltered bays and lake margins that trade no water;
- water mixed to 10 m, enclosed: the body of a lake's upper water, the open water of the land's cells;
- shore water 2 m deep trading three tenths, or all, of its water a day with the open water's body: the run's own sea
  surface on its water cells, the 10 m water on its land cells;
- the sea surface the run itself carries, a mixed layer 50 m deep, on its water cells: the body of the sea's upper
  water.

A mixed layer has one temperature through its depth. Calm deep water does not: through the lunar day its top few
metres warm above the water below (the seas' oxygen study finds a day's mixed layer of 3-6 m in its columns). So the
body's figures are the surface's only where waves, currents and tides keep the surface mixed into it, and the still
layers 2 and 4 m deep show what a calm surface does. As a check the lake model's own column, the forcing's 20 m pan
under the cell's unsheltered wind, runs at the forcing's five sites, and its surface is recorded beside the layers
there (checks.sites, calm_column).

Every cell takes its own lowest level's wind unsheltered, the speed of each window's mean vector floored at PlaSim's
least wind and raised by the forcing's factor (forcing.level_wind, forcing.wind_raise), over water's own roughness: on
a land cell that is the wind a water wide enough for the air to settle over it would meet. A land cell's
waters stand at its mean ground level in the run. The statistics are over the hourly steps of the ten years after two
lunar days of spin-up: the mean, the 1st and 99th percentiles (from counts in quarter-kelvin bins), the extremes and
the hours of each lunar day above 34, 36 and 38 C, each rounded to a hundredth.

The shares weight each class by what the waters are: the seas by the area of the run's water cells, for their body
and for still layers there; the sea coast by
its length in each land cell (the moisture map's settings, from geography's atlas), each stretch taking the mean of
the neighbouring water cells' shore values, or the land cell's own where it has no water cell beside it; and the
rain-fed lakes by their area in each land cell (the moisture map's settings).

    /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.littoral_community.water_map  # minutes

writes results/water_map.json (schema terluna.research.littoral-water-map/1) with the sha256 of every input.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from biosphere.ecology.moisture_map import cell_areas_km2, cell_edges
from research.studies.littoral_community import forcing as F
from shared.constants import (EXOPLASIM_CHARNOCK, EXOPLASIM_SEA_ROUGHNESS_FLOOR, EXOPLASIM_VON_KARMAN,
                              MOON_SURFACE_GRAVITY, STANDARD_GRAVITY)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MOISTURE = ROOT / 'biosphere' / 'ecology' / 'results' / 'moisture_map.json'
FORCING = HERE / 'results' / 'forcing.json'
OUT = HERE / 'results' / 'water_map.json'
SCHEMA = 'terluna.research.littoral-water-map/1'
SB = 5.670374419e-8

# name, depth (m), share of its water traded with the open water each day
CASES = (('still_1m', 1.0, 0.0), ('still_2m', 2.0, 0.0), ('still_4m', 4.0, 0.0), ('mixed_10m', 10.0, 0.0),
         ('shore_2m_x0.3', 2.0, 0.3), ('shore_2m_x1', 2.0, 1.0))
OPEN_CASE = 3                               # the 10 m water, the land cells' open water
BIN0, BIN_W, NBINS = 5.0, 0.25, 240         # the counts' bins: 5 to 65 C by a quarter kelvin
THRESHOLDS_C = (34.0, 36.0, 38.0)
CLASS_EDGES_C = (28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0)
CLASS_NAMES = ('under 28 °C', '28–30 °C', '30–32 °C', '32–34 °C', '34–36 °C', '36–38 °C', '38–40 °C', '40 °C and more')
CENTRES = np.arange(-89.0, 89.01, 2.0)      # the clear-sky tables' elevations: forcing.run_clear_sky's bin centres
LATITUDE_BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))
# forcing.site_forcing shapes photosynthetic light from a clear-sky atlas; the water's heat reads none of it
NO_PAR = dict(elevation=np.array([-90.0, 90.0]), shortwave=np.array([1.0, 1.0]), par0=np.array([0.0, 0.0]))


# ---- the lake model's bulk exchange, for arrays ---------------------------------------------------------------------------
def saturation_humidity(t_k, p_pa):
    e = 611.2 * np.exp(17.67 * (t_k - 273.15) / (t_k - 29.65))
    return 0.622 * e / (p_pa - 0.378 * e)


def neutral(u, z, gravity):
    """lake.neutral_coefficient for arrays: PlaSim's neutral coefficient over water at height z and its roughness."""
    z0 = np.full(np.broadcast(u, z).shape, EXOPLASIM_SEA_ROUGHNESS_FLOOR)
    for _ in range(4):
        cn = (EXOPLASIM_VON_KARMAN / np.log(z / z0 + 1.0)) ** 2
        z0 = np.maximum(EXOPLASIM_CHARNOCK * cn * u * u / gravity, EXOPLASIM_SEA_ROUGHNESS_FLOOR)
    return (EXOPLASIM_VON_KARMAN / np.log(z / z0 + 1.0)) ** 2, z0


def reduce_wind(u_level, z_level, z=10.0, gravity=MOON_SURFACE_GRAVITY):
    """lake.reduce_wind for arrays of positive wind: the wind at z under u_level at z_level on the neutral profile
    over water."""
    cn, z0 = neutral(u_level, z_level, gravity)
    return np.sqrt(cn) * u_level / EXOPLASIM_VON_KARMAN * np.log(z / z0 + 1.0)


def bulk_fluxes(ts_c, air_k, q_air, wind, p_pa, air: L.Air, z=None):
    """lake.bulk_fluxes for arrays: sensible and latent heat leaving the water (W/m2), the wind, air and humidity
    paired at height z (the air's 10 m reference where none is given)."""
    tsk = ts_c + 273.15
    u = np.maximum(wind, air.min_wind)
    qs = saturation_humidity(tsk, p_pa)
    q_air = np.minimum(q_air, saturation_humidity(air_k, p_pa))
    tv_air, tv_s = air_k * (1.0 + 0.608 * q_air), tsk * (1.0 + 0.608 * qs)
    rho_air = p_pa / (air.gas_constant * tv_air)
    z = air.wind_height_m if z is None else z
    cn, _ = neutral(u, z, air.gravity)
    ri = air.gravity * z * (tv_air - tv_s) / (0.5 * (tv_air + tv_s) * u * u)
    g_fc = air.gravity if air.free_convection_gravity is None else air.free_convection_gravity
    w_fc = L.FREE_CONVECTION * (g_fc / STANDARD_GRAVITY) ** (1.0 / 3.0) * np.cbrt(np.maximum(tv_s - tv_air, 0.0))
    unstable = (1.0 + (w_fc / (u * cn)) ** 1.25) ** 0.8
    rs = np.maximum(ri, 0.0)
    stable = 1.0 / (1.0 + 3.0 * L.LOUIS_B * rs * np.sqrt(1.0 + L.LOUIS_D * rs))
    ch = cn * np.where(ri <= 0.0, unstable, stable) * u * rho_air
    latent = (2.501e6 - 2370.0 * ts_c) * ch * (qs - q_air)
    return ch * air.cp * (tsk - air_k), latent


# ---- the Sun, the tables and the windows ------------------------------------------------------------------------------------
def global_sun(data):
    """The noon time at longitude zero (hours from the first window's start) and the period, fitted on every cell of
    the row nearest the equator and carried along the longitude; with the spread of the cells' own noons about it."""
    lat, lon = data['lat'], data['lon']
    j = int(np.argmin(np.abs(lat)))
    fits = [F.fit_sun(data['czen'][:, j, i].astype(float), float(lat[j])) for i in range(lon.size)]
    period = float(np.median([f['period_h'] for f in fits]))
    at_zero = np.array([f['noon_h'] + lon[i] / 360.0 * period for i, f in enumerate(fits)])
    angle = 2.0 * np.pi * at_zero / period
    noon0 = float((np.arctan2(np.sin(angle).mean(), np.cos(angle).mean()) * period / (2.0 * np.pi)) % period)
    spread = (at_zero - noon0 + period / 2.0) % period - period / 2.0
    return dict(noon0_h=noon0, period_h=period, row_lat=float(lat[j]), noon_spread_h=float(np.abs(spread).max()),
                fit_rms_max=float(max(f['rms'] for f in fits)))


def elevation(lat2, noon2, period, t_h):
    """The Sun's elevation (degrees) over every cell at a time, the declination zero."""
    return np.degrees(np.arcsin(np.cos(np.radians(lat2)) * np.cos(2.0 * np.pi * (t_h - noon2) / period)))


def lookup(tables, elev):
    """Each cell's table read at its elevation, linear between the centres and constant beyond them (np.interp)."""
    x = np.clip((elev - CENTRES[0]) / (CENTRES[1] - CENTRES[0]), 0.0, CENTRES.size - 1.0)
    i0 = np.minimum(x.astype(int), CENTRES.size - 2)
    f = x - i0
    rows = np.arange(tables.shape[0])
    return tables[rows, i0] * (1.0 - f) + tables[rows, i0 + 1] * f


def between(k, n_w, w):
    """The two windows whose centres bracket step k of window w, and the later one's weight (np.interp's ends)."""
    x = (k + 0.5) / F.STEPS_PER_WINDOW - 0.5
    if x < 0.0:
        a, b, wt = w - 1, w, x + 1.0
    else:
        a, b, wt = w, w + 1, x
    if a < 0:
        return 0, 0, 0.0
    if b > n_w - 1:
        return n_w - 1, n_w - 1, 0.0
    return a, b, wt


def windows(data, sun):
    """Each window's fields for every cell (windows x cells), the clear-sky tables and the cloud factors."""
    n_w = data['tas'].shape[0]
    flat = lambda k: data[k].reshape(n_w, -1).astype(float)
    lat2 = np.repeat(data['lat'], data['lon'].size)
    lon2 = np.tile(data['lon'], data['lat'].size)
    noon2 = sun['noon0_h'] - lon2 / 360.0 * sun['period_h']
    tas, alb, ts = flat('tas'), flat('alb'), flat('ts')
    lit_ground = alb < 0.999
    sw_all = np.where(lit_ground, flat('rss') / np.maximum(1.0 - alb, 1e-3), 0.0)
    clear_end = np.where(lit_ground, flat('rsscs') / np.maximum(1.0 - alb, 1e-3), 0.0)
    q_open = flat('hus_low')
    q = q_open.copy()
    q[:-1] = 0.5 * (q_open[:-1] + q_open[1:])            # the output's humidity is the value opening each window
    p = flat('ps') * 100.0
    q = np.minimum(q, saturation_humidity(tas, p))
    # the lowest level's wind and height, as forcing.level_wind takes them at one cell
    z_low = F.air_constants()[0] * 0.5 * (flat('ta_low') + tas) / MOON_SURFACE_GRAVITY * math.log(1.0 / data['sigma_low'])
    level = np.maximum(np.hypot(flat('ua_low'), flat('va_low')), F.LEVEL_FLOOR_MS) * F.wind_raise()
    wind = reduce_wind(level, z_low)
    lw = flat('rls') + SB * ts ** 4
    # the run's clear sky by the Sun's elevation, each cell's own (forcing.run_clear_sky on its end-of-window values)
    elev_end = np.stack([elevation(lat2, noon2, sun['period_h'], (w + 1.0) * F.WINDOW_H) for w in range(n_w)])
    tables = np.zeros((lat2.size, CENTRES.size))
    for c in range(lat2.size):
        e, f = F.run_clear_sky(clear_end[:, c], elev_end[:, c])
        tables[c] = np.interp(CENTRES, e, f)
    # each window's mean clear sky and its cloud factor
    dt_h = F.WINDOW_H / F.STEPS_PER_WINDOW
    mean_clear = np.zeros((n_w, lat2.size))
    for w in range(n_w):
        for k in range(F.STEPS_PER_WINDOW):
            mean_clear[w] += lookup(tables, elevation(lat2, noon2, sun['period_h'], (w * F.STEPS_PER_WINDOW + k + 0.5) * dt_h))
    mean_clear /= F.STEPS_PER_WINDOW
    lit = mean_clear > F.CLEAR_FLOOR_W
    cloud_all = np.where(lit, sw_all, 0.0).sum(0) / np.maximum(np.where(lit, mean_clear, 0.0).sum(0), 1e-9)
    ratio = np.where(lit, sw_all / np.maximum(mean_clear, 1e-9), cloud_all[None, :])
    return dict(n_w=n_w, lat2=lat2, lon2=lon2, noon2=noon2, air=tas, q=q, p=p, wind=wind, level=level, z=z_low, lw=lw,
                sea_c=ts - 273.15,
                tables=tables, ratio=ratio, cloud_all=cloud_all)


# ---- the water --------------------------------------------------------------------------------------------------------------
def run(win, sun, sea, cases=CASES, spin_windows=F.SPIN_WINDOWS, progress=None):
    """Step every cell's cases through the windows and count their temperatures. Returns the counts
    (cases + 1, cells, bins), the last being the sea surface the run carries, with each one's sum, lowest and highest."""
    rd, cp, _, _ = F.air_constants()
    air = L.Air(rd, cp, MOON_SURFACE_GRAVITY)
    ph = L.Physics()
    n_w, n_c = win['n_w'], win['lat2'].size
    dt = F.WINDOW_H * 3600.0 / F.STEPS_PER_WINDOW
    dt_h = dt / 3600.0
    depth = np.array([c[1] for c in cases])[:, None]
    keep = np.exp(-np.array([c[2] for c in cases]) * dt / 86400.0)[:, None]
    trades = np.array([c[2] > 0.0 for c in cases])
    n_k = len(cases)
    temp = np.repeat((win['air'][0] - 273.15)[None, :], n_k, 0)
    counts = np.zeros((n_k + 1, n_c, NBINS), dtype=np.int32)
    total = np.zeros((n_k + 1, n_c))
    low, high = np.full((n_k + 1, n_c), np.inf), np.full((n_k + 1, n_c), -np.inf)
    ki, ci = np.meshgrid(np.arange(n_k + 1), np.arange(n_c), indexing='ij')
    steps = 0
    for w in range(n_w):
        for k in range(F.STEPS_PER_WINDOW):
            a, b, wt = between(k, n_w, w)
            mix = lambda x: x[a] * (1.0 - wt) + x[b] * wt
            air_k, q, p, level, z, lw, sea_c = (mix(win[key]) for key in ('air', 'q', 'p', 'level', 'z', 'lw', 'sea_c'))
            t_h = (w * F.STEPS_PER_WINDOW + k + 0.5) * dt_h
            sw = np.maximum(lookup(win['tables'], elevation(win['lat2'], win['noon2'], sun['period_h'], t_h)) * win['ratio'][w], 0.0)
            sens, latent = bulk_fluxes(temp, air_k[None, :], q[None, :], level[None, :], p[None, :], air, z[None, :])
            emitted = L.WATER_EMISSIVITY * (SB * (temp + 273.15) ** 4 - lw[None, :])
            net = (1.0 - ph.albedo) * sw[None, :] - emitted - sens - latent + ph.heat_flow
            temp = temp + net * dt / (L.RHO_CP * depth)
            open_c = np.where(sea, sea_c, temp[OPEN_CASE])
            temp[trades] = open_c[None, :] + (temp[trades] - open_c[None, :]) * keep[trades]
            if w >= spin_windows:
                both = np.concatenate([temp, sea_c[None, :]])
                bins = np.clip(((both - BIN0) / BIN_W).astype(int), 0, NBINS - 1)
                counts[ki, ci, bins] += 1
                total += both
                low, high = np.minimum(low, both), np.maximum(high, both)
                steps += 1
        if progress and (w + 1) % 200 == 0:
            progress(w + 1, n_w)
    return dict(counts=counts, total=total, low=low, high=high, steps=steps, dt_h=dt_h)


def percentile(counts, q):
    """The q-th percentile of each row of counts, linear within its bin."""
    cum = np.cumsum(counts, -1)
    target = q / 100.0 * cum[..., -1:]
    idx = np.argmax(cum >= target, -1)
    prev = np.take_along_axis(cum, np.maximum(idx - 1, 0)[..., None], -1)[..., 0]
    prev = np.where(idx > 0, prev, 0)
    here = np.take_along_axis(counts, idx[..., None], -1)[..., 0]
    frac = (target[..., 0] - prev) / np.maximum(here, 1)
    return BIN0 + (idx + frac) * BIN_W


def summarise(out):
    counts = out['counts']
    lunar_days = out['steps'] * out['dt_h'] / F.SYNODIC_H
    edges = BIN0 + np.arange(NBINS + 1) * BIN_W
    stats = dict(mean_c=out['total'] / out['steps'], p01_c=percentile(counts, 1.0), p99_c=percentile(counts, 99.0),
                 min_c=out['low'], max_c=out['high'])
    for x in THRESHOLDS_C:
        above = counts[..., edges[:-1] >= x].sum(-1)
        stats[f'hours_above_{x:g}'] = above * out['dt_h'] / lunar_days
    return stats


def classify(t_c):
    return np.searchsorted(np.array(CLASS_EDGES_C), np.asarray(t_c), side='right')


def shares(classes, weight):
    w = np.asarray(weight, float)
    total = w.sum()
    return [float(w[classes == k].sum() / total) if total > 0 else 0.0 for k in range(len(CLASS_NAMES))]


def coast_values(values2, sea2, land_rows, land_cols):
    """For each listed land cell, the mean of a field over its neighbouring water cells (eight neighbours, the
    longitude wrapping), or the cell's own value where no water cell lies beside it; and whether one did."""
    n_lat, n_lon = sea2.shape
    out, beside = np.zeros(len(land_rows)), np.zeros(len(land_rows), bool)
    for n, (j, i) in enumerate(zip(land_rows, land_cols)):
        got = [values2[jj, (i + di) % n_lon] for jj in (j - 1, j, j + 1) if 0 <= jj < n_lat for di in (-1, 0, 1)
               if sea2[jj, (i + di) % n_lon]]
        beside[n] = bool(got)
        out[n] = float(np.mean(got)) if got else float(values2[j, i])
    return out, beside


def calm_column(data, site, depth=F.OPEN_LAKE_M):
    """The lake model's own column at one of the forcing's sites: its still pan of the open lake's depth under the
    cell's unsheltered wind, by the site's own fitted Sun. Returns its surface's statistics after the spin-up, with
    the floor's mean."""
    here = dict(j=site['j'], i=site['i'], lat=site['lat'], lon=site['lon'], wind_share=1.0)
    forcing = F.site_forcing(data, here, site['sun'], NO_PAR)[0]
    out = F.run_water(forcing, depth, float(forcing.air_k[0] - 273.15))
    start = F.SPIN_WINDOWS * F.STEPS_PER_WINDOW
    surface = np.asarray(out['surface_t'], float)[start:]
    dt_h = forcing.dt / 3600.0
    lunar_days = surface.size * dt_h / F.SYNODIC_H
    row = dict(depth_m=depth, wind_mean_m_s=float(np.mean(forcing.wind10)), surface_mean_c=float(surface.mean()),
               surface_p99_c=float(np.percentile(surface, 99)), surface_max_c=float(surface.max()),
               floor_mean_c=float(np.mean(np.asarray(out['bottom_t'], float)[start:])))
    for x in THRESHOLDS_C:
        row[f'surface_hours_above_{x:g}'] = float((surface > x).sum() * dt_h / lunar_days)
    return {k: round(v, 2) for k, v in row.items()}


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    started = time.time()
    data = F.read_years()
    lat, lon = data['lat'], data['lon']
    shape = (lat.size, lon.size)
    sea2 = data['lsm'] <= 0.5
    sea = sea2.ravel()
    sun = global_sun(data)
    win = windows(data, sun)
    out = run(win, sun, sea, progress=lambda a, b: print(f'windows {a}/{b}', flush=True))
    # rounded here, so that the classes and shares below are those anyone counts from the product's cells
    stats = {k: np.round(v, 2) for k, v in summarise(out).items() if k != 'min_c'}
    names = [c[0] for c in CASES] + ['sea_surface']
    area2 = cell_areas_km2(lat, lon)
    coast2 = F.coastal(data['lsm'])
    air_mean = win['air'].mean(0) - 273.15

    # the settings that weight the shares: the moisture map's coast lengths and rain-fed lakes by land cell
    moist = json.loads(MOISTURE.read_text())
    mc = moist['cells']
    rows, cols = np.array(mc['row']), np.array(mc['col'])
    if not np.allclose(moist['grid']['lat'], lat) or not np.allclose(moist['grid']['lon'], lon):
        raise ValueError('the moisture map stands on another grid')
    if np.any(sea2[rows, cols]):
        raise ValueError("the moisture map's land cells include water cells of this run")
    coast_km = np.array(mc['coast_km'], float)
    lake_km2 = np.array(mc['rainfed_lake_fraction'], float) * np.array(mc['area_km2'], float)

    def field(name, key):
        return stats[key][names.index(name)].reshape(shape)

    share = {}
    for key in ('mean_c', 'p99_c'):
        s = {}
        s['seas_surface'] = shares(classify(field('sea_surface', key)[sea2]), area2[sea2])
        for name in ('still_1m', 'still_2m', 'still_4m', 'mixed_10m'):
            s[f'seas_{name}'] = shares(classify(field(name, key)[sea2]), area2[sea2])
        for name in ('still_1m', 'still_2m', 'still_4m', 'shore_2m_x0.3', 'shore_2m_x1'):
            v, _ = coast_values(field(name, key), sea2, rows, cols)
            s[f'sea_coast_{name}'] = shares(classify(v), coast_km)
        for name in ('still_1m', 'still_2m', 'still_4m', 'mixed_10m', 'shore_2m_x0.3'):
            s[f'lakes_{name}'] = shares(classify(field(name, key)[rows, cols]), lake_km2)
        share[key] = s
    _, beside = coast_values(field('still_2m', 'p99_c'), sea2, rows, cols)
    lat_cells = np.repeat(lat, lon.size).reshape(shape)
    wind_mean = win['wind'].mean(0).reshape(shape)
    by_latitude = []
    for lo, hi in LATITUDE_BANDS:
        band = (np.abs(lat_cells) >= lo) & (np.abs(lat_cells) < hi)
        row = dict(band_deg=[lo, hi])
        for label, mask in (('water_cells', band & sea2), ('land_cells', band & ~sea2)):
            mean = lambda a: round(float((a * area2)[mask].sum() / area2[mask].sum()), 2)
            row[label] = dict(still_2m_p99_c=mean(field('still_2m', 'p99_c')), still_2m_mean_c=mean(field('still_2m', 'mean_c')),
                              still_4m_p99_c=mean(field('still_4m', 'p99_c')),
                              still_2m_hours_above_36=mean(field('still_2m', 'hours_above_36')),
                              still_4m_hours_above_36=mean(field('still_4m', 'hours_above_36')),
                              shore_2m_x03_p99_c=mean(field('shore_2m_x0.3', 'p99_c')),
                              mixed_10m_p99_c=mean(field('mixed_10m', 'p99_c')),
                              air_mean_c=mean(air_mean.reshape(shape)), wind_mean_m_s=mean(wind_mean))
            if label == 'water_cells':
                row[label]['sea_surface_mean_c'] = mean(field('sea_surface', 'mean_c'))
                row[label]['sea_surface_p99_c'] = mean(field('sea_surface', 'p99_c'))
        by_latitude.append(row)

    # the five sites of the forcing product, by its own per-site Sun and series
    forcing = json.loads(FORCING.read_text())
    sites = {}
    for name, site in forcing['sites'].items():
        c = site['j'] * lon.size + site['i']
        own_noon = site['sun']['noon_h']
        here = (win['noon2'][c] - own_noon + sun['period_h'] / 2.0) % sun['period_h'] - sun['period_h'] / 2.0
        row = dict(noon_minus_site_fit_h=float(here))
        for mine, theirs in (('still_2m', None), ('shore_2m_x0.3', 'slab_2m_x0.3'), ('shore_2m_x1', 'slab_2m_x1')):
            if theirs and site['kind'] == 'sea shore':
                row[mine] = dict(p99_c=float(stats['p99_c'][names.index(mine), c]),
                                 forcing_p99_c=site['cases'][theirs]['water_p99_c'],
                                 mean_c=float(stats['mean_c'][names.index(mine), c]),
                                 forcing_mean_c=site['cases'][theirs]['water_mean_c'])
        row['still_2m_p99_c'] = float(stats['p99_c'][names.index('still_2m'), c])
        row['still_4m_p99_c'] = float(stats['p99_c'][names.index('still_4m'), c])
        row['body_p99_c'] = float(stats['p99_c'][names.index('sea_surface' if sea[c] else 'mixed_10m'), c])
        row['forcing_closed_2m_pan_p99_c'] = site['cases']['pan_2m']['water_p99_c']
        row['calm_column'] = calm_column(data, site)
        print('calm column', name, row['calm_column'], flush=True)
        sites[name] = row
    # the 10 m layer under the run's near-surface air against the run's own sea surface, on its water cells
    apart = (field('mixed_10m', 'mean_c') - field('sea_surface', 'mean_c'))[sea2]
    against_run = dict(mean_c=round(float((apart * area2[sea2]).sum() / area2[sea2].sum()), 2),
                       range_c=[float(apart.min()), float(apart.max())])

    r2 = lambda a: [round(float(x), 2) for x in np.asarray(a).ravel()]
    cells = dict(row=np.repeat(np.arange(lat.size), lon.size).tolist(), col=np.tile(np.arange(lon.size), lat.size).tolist(),
                 sea=sea.astype(int).tolist(), coastal_sea=coast2.ravel().astype(int).tolist(),
                 area_km2=r2(area2), air_mean_c=r2(air_mean), wind_mean_m_s=r2(win['wind'].mean(0)),
                 wind_level_mean_m_s=r2(win['level'].mean(0)),
                 run_wind10_vector_mean_m_s=r2(np.hypot(data['uas'], data['vas']).reshape(win['n_w'], -1).mean(0)),
                 cloud_factor=r2(win['cloud_all']))
    for n, name in enumerate(names):
        cells[name] = {k: r2(v[n]) for k, v in stats.items()}
    files_used = ['research/studies/littoral_community/water_map.py', 'research/studies/littoral_community/forcing.py',
                  'biosphere/ecology/lake.py', 'biosphere/ecology/moisture_map.py']
    sea_mean = field('sea_surface', 'mean_c')
    record = dict(
        schema=SCHEMA,
        producer=dict(lane='research', runner='research/studies/littoral_community/water_map.py',
                      files={f: F.sha256(ROOT / f) for f in files_used}, constants=constants_used(files_used)),
        evidence=('Model results on the corrected design climate (the run named in settings, years 10-19) at T21: '
                  'mixed slabs of water under each cell\'s sunlight, sky, air and 10 m wind, with the lake model\'s '
                  'bulk exchange paired at the run\'s lowest level as the waters\' oxygen study pairs it. The still '
                  'cases trade no water and take that level\'s wind over water, raised by the forcing\'s one factor '
                  'for the wind inside a 3-day mean, with no shelter, no sea breeze and no gust of the boundary '
                  'layer\'s eddies; the shore cases\' exchange is a '
                  'stated share, taken with the open water\'s body. A calm deep water\'s own surface late in the '
                  'lunar day lies above its body: at the five sites of checks the lake model\'s 20 m column puts it '
                  'between the still 4 m and 2 m layers. The climate model and CM1 differ on the air near the '
                  'ground, a gap the project carries as a range.'),
        reading_rule=('Temperatures in C of a mixed layer of the stated depth, over hourly steps of ten years, rounded '
                      'to a hundredth before the classes and shares are counted. The run\'s water cells are the seas '
                      'and the large lakes; there sea_surface is the run\'s own surface, a mixed layer 50 m deep: the '
                      'body of the upper water. On land cells the waters stand at the cell\'s mean ground level and '
                      'mixed_10m is the body of a lake\'s upper water. Read the body as the surface only where the '
                      'surface is kept mixed into it; a calm surface follows the still layers. hours_above_X '
                      'are hours of each lunar day (708.7 h). Classes are bins of 2 K from under 28 C to 40 C and '
                      'more. by_latitude gives area-weighted means over the cells of each band of |latitude|; '
                      'wind_level_mean_m_s is the lowest level\'s wind the waters take, floored and raised; '
                      'wind_mean_m_s is that wind at 10 m over water; run_wind10_vector_mean_m_s is the run\'s own '
                      '10 m wind, the speed of its mean vectors over each cell\'s own ground.'),
        settings=dict(run=F.RUN, years=list(F.YEARS), cases=[dict(name=c[0], depth_m=c[1], exchange_per_day=c[2]) for c in CASES],
                      open_water=dict(water_cells="the run's sea surface", land_cells='mixed_10m'),
                      spin_windows=F.SPIN_WINDOWS, steps_per_window=F.STEPS_PER_WINDOW, window_h=F.WINDOW_H,
                      pairing='the lowest level, as PlaSim pairs it', lowest_sigma=data['sigma_low'],
                      level_floor_m_s=F.LEVEL_FLOOR_MS, wind_raise=F.wind_raise(),
                      bins=dict(first_c=BIN0, width_k=BIN_W, count=NBINS), thresholds_c=list(THRESHOLDS_C)),
        classes=dict(edges_c=list(CLASS_EDGES_C), names=list(CLASS_NAMES)),
        grid=dict(lat=lat.tolist(), lon=lon.tolist(), lat_edges=cell_edges(lat).tolist(),
                  note="The run's own frame, the turned Moon's; the near side is centred on longitude 0."),
        sun=sun,
        summary=dict(
            seas=dict(area_km2=float(area2[sea2].sum()), surface_mean_c=float((sea_mean * area2)[sea2].sum() / area2[sea2].sum()),
                      surface_mean_range_c=[float(sea_mean[sea2].min()), float(sea_mean[sea2].max())]),
            coast=dict(length_km=float(coast_km.sum()), share_beside_a_water_cell=float(coast_km[beside].sum() / coast_km.sum())),
            lakes=dict(area_km2=float(lake_km2.sum())),
            by_latitude=by_latitude),
        shares=share,
        checks=dict(sites=sites, mixed_10m_minus_the_runs_sea_surface_means=against_run),
        cells=cells,
        inputs={str(p.relative_to(F.RUNS)): F.sha256(p) for p in F.year_files()}
        | {str(MOISTURE.relative_to(ROOT)): F.sha256(MOISTURE), str(FORCING.relative_to(ROOT)): F.sha256(FORCING)},
        design_air=dict(zip(('air_gas_constant', 'air_cp', 'o2_fraction', 'salinity_g_kg'), F.air_constants())),
        seconds=round(time.time() - started, 1))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(json.dumps(dict(sun=sun, summary=record['summary'], sites=sites), indent=1))
    print(OUT, record['seconds'], 's')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
