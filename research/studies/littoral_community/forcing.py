"""The littoral's physical setting through the design climate's lunar days: the water's temperature, its oxygen
exchange with the air and the light entering it, at shore sites of the corrected design climate.

The climate is the corrected design run on the turned geography, A28R_301k_design_l15_snow_s938_turned_refit_precession,
years 10-19, the design years (register, 8 November 2026): 1,200 3-day means at T21. Five sites stand for the
littoral's settings: four sea cells on a coast (a sea cell with land among its eight neighbours), one in each of the
latitude bands 0-10, 25-35, 40-50 and 55-65 degrees in the run's own frame, each the coastal cell with the band's median
ten-year surface temperature; and one land cell within 15 degrees of the equator with the median land temperature,
standing for the margins of the rain-fed lakes, 92% of whose lit floor lies within 30 degrees of the equator
(research/studies/underwater_light).

The Sun. The run's output windows (143 half-hour steps, 71.5 hours) do not divide its lunar day, so the Sun's phase
drifts from one lunar day's windows to the next. The noon time and the period are fitted, by least squares, to the
windows' mean cosine of the zenith angle (czen, night counted as zero), with the declination taken as zero (the
seasonal swing is 1.54 degrees). The fitted Sun then gives every step's elevation.

The light. The run writes its all-sky fluxes as window means and its clear-sky diagnostics (rsscs, rstcs) as the
values of each window's last step: against the fitted Sun, the clear-sky surface field follows the Sun at the
window's end (correlation 0.992-0.997 at the five sites, against 0.92 with the window's mean Sun) and falls to
twilight's 10-12 W/m2 when the Sun has set by then, where a mean over the window would not. So the run's clear sky is read as a function of the Sun's elevation, the medians of those end-of-window
values in 2-degree bins of the fitted Sun's elevation at each window's end, and each window's clear-sky mean is that
function averaged over the window's steps. A window's cloud factor is its all-sky downwelling shortwave over that
clear-sky mean; in windows too dark for a ratio (deep twilight and night) the ten years' ratio of the two is used.
Photosynthetic light just below the surface is the solved spherical sky's (research/runs/optical_comfort/spherical, as
research/studies/underwater_light reads it: the clear sky's photosynthetic photons below a flat surface at every solar
elevation, the long twilight included) times the window's cloud factor, clouds taken to dim photosynthetic light as
they dim the whole shortwave.

The water. The project's lake model (biosphere/ecology/lake.py) on pans of 1, 2 and 4 m, forced as the cold-waters and
waters' oxygen studies force it from the run's means (research/studies/cold_waters/run.py, site_forcing): the air and
the lowest level's humidity, the sky infrared rls + sigma ts^4, the surface pressure and the 10 m wind, each linear in
time between the windows' centres. The sunlight is sharpened: each window keeps the run's all-sky downwelling
shortwave, rss / (1 - alb), as its mean, shared over the window's steps in proportion to the run's clear-sky shortwave
at each step's elevation, so the water warms from dawn and cools from dusk instead of across a smeared 3-day mean. The
lake margin takes a third of its land cell's 10 m wind (SHELTER), as the cold-waters study's sheltered ponds do; the
open lake beside it takes the whole of it, a wind over land's roughness that a lake's fetch would raise, so the lake
site's daytime heat is an upper-side estimate. The oxygen transfer velocity at each step is the lake model's
surface-renewal form (Tedford et al. 2014) on the run's water friction velocity and the cooling buoyancy flux of the
surface's own heat loss (its infrared, sensible and latent heat at the recorded surface temperature), with Cole and
Caraco's (1998) wind form beside it.

A closed pan holds its heat, while a shore's water is changed by waves, wind-driven currents and the monthly tide.
So beside each pan runs a slab of the same depth that trades water with the open water at a stated share a day
(EXCHANGE_PER_DAY): its heat follows the same surface fluxes (the lake model's bulk exchange and emission, all the
absorbed sunlight heating the slab as the lake model lets the floor's light heat the water over it) and relaxes
toward the open water's temperature, the run's sea surface at the sea cells (its 50 m mixed layer, linear between
the windows' centres) and a 20 m pan's surface at the lake margin. The slab's oxygen transfer velocity takes the same
surface-renewal form on its own wind stress and cooling.

    /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.littoral_community.forcing  # minutes

writes the step series to /media/projectspace/terluna-research/littoral-community/forcing/<site>.npz and
results/forcing.json (schema terluna.research.littoral-forcing/1) with the sha256 of every input and series.
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
from research.studies.underwater_light import optics
from research.studies.underwater_light import run as UL
from shared.constants import MOON_SURFACE_GRAVITY, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = Path('/media/projectspace/terluna-research/gcm-runs')
RUN, YEARS = 'A28R_301k_design_l15_snow_s938_turned_refit_precession', (10, 19)
LAKE_OXYGEN = ROOT / 'research' / 'studies' / 'lake_oxygen' / 'results' / 'lake_oxygen.json'
BULK = Path('/media/projectspace/terluna-research/littoral-community/forcing')
OUT = HERE / 'results' / 'forcing.json'
SCHEMA = 'terluna.research.littoral-forcing/1'

WINDOW_STEPS = 143                      # the run's half-hour steps in each output window
WINDOW_H = WINDOW_STEPS * 0.5           # 71.5 hours
STEPS_PER_WINDOW = 72                   # lake-model steps in each window (3,575 s), as the waters' oxygen study
SPIN_WINDOWS = 20                       # windows left out before reading: two lunar days
BANDS = {'shore_0': (0.0, 10.0), 'shore_30': (25.0, 35.0), 'shore_45': (40.0, 50.0), 'shore_60': (55.0, 65.0)}
LAKE_BAND = (0.0, 15.0)
DEPTHS_M = (1.0, 2.0, 4.0)
LAYER_M = {1.0: 0.1, 2.0: 0.1, 4.0: 0.2, 20.0: 0.5}
SHELTER = 0.3                           # a sheltered margin's share of the 10 m wind (cold_waters' ponds)
CLEAR_FLOOR_W = 5.0                     # W/m2: below this window-mean clear shortwave the cloud factor is not read
CLOUD_CAP = 1.05                        # the largest cloud factor photosynthetic light takes from a window
PAR_DARK = 1.0                          # umol m-2 s-1, the underwater light study's darkness
THRESHOLDS_C = (32.0, 34.0, 36.0, 38.0)
EXCHANGE_PER_DAY = (0.1, 0.3, 1.0)      # the share of a slab's water the open water replaces each day
OPEN_LAKE_M = 20.0                      # the open lake beside the lake margin, a pan of this depth
SYNODIC_H = SYNODIC_MONTH_DAYS * 24.0
STEFAN_BOLTZMANN = 5.670374419e-8


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def year_files():
    return [RUNS / RUN / 'model' / f'MOST.{y:05d}.nc' for y in range(YEARS[0], YEARS[1] + 1)]


def read_years():
    """The ten years' windows of the fields the forcing needs (float32), with the grid."""
    import netCDF4
    keys = ('tas', 'ts', 'rss', 'rsscs', 'rls', 'alb', 'ps', 'uas', 'vas', 'czen')
    out = {k: [] for k in keys + ('hus_low',)}
    for path in year_files():
        with netCDF4.Dataset(path) as nc:
            lat, lon = np.array(nc['lat'][:], dtype=float), np.array(nc['lon'][:], dtype=float)
            lsm = np.array(nc['lsm'][0], dtype=float)
            for k in keys:
                out[k].append(np.array(nc[k][:], dtype=np.float32))
            out['hus_low'].append(np.array(nc['hus'][:, -1], dtype=np.float32))
    data = {k: np.concatenate(v) for k, v in out.items()}
    data.update(lat=lat, lon=lon, lsm=lsm)
    return data


def coastal(lsm):
    """Sea cells with land among their eight neighbours (longitude wraps; no wrap across the poles)."""
    land = lsm > 0.5
    near = np.zeros_like(land)
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            if dj or di:
                shifted = np.roll(land, di, 1)
                if dj > 0:
                    shifted = np.vstack([np.zeros_like(shifted[:1]), shifted[:-1]])
                elif dj < 0:
                    shifted = np.vstack([shifted[1:], np.zeros_like(shifted[:1])])
                near |= shifted
    return (~land) & near


def median_cell(data, mask):
    """The cell of a mask with the median ten-year mean surface temperature (the lower middle for an even count)."""
    js, iss = np.nonzero(mask)
    t = data['ts'][:, js, iss].astype(float).mean(0)
    order = np.argsort(t, kind='stable')
    k = order[(len(order) - 1) // 2]
    return int(js[k]), int(iss[k]), float(t[k] - 273.15), int(len(order))


def choose_sites(data):
    lat = data['lat']
    sites = {}
    coast = coastal(data['lsm'])
    for name, (lo, hi) in BANDS.items():
        band = (np.abs(lat) >= lo) & (np.abs(lat) < hi)
        j, i, t, n = median_cell(data, coast & band[:, None])
        sites[name] = dict(kind='sea shore', j=j, i=i, lat=float(lat[j]), lon=float(data['lon'][i]), wind_share=1.0,
                           cells_in_band=n, ten_year_surface_c=round(t, 3))
    land = data['lsm'] > 0.5
    band = (np.abs(lat) >= LAKE_BAND[0]) & (np.abs(lat) < LAKE_BAND[1])
    j, i, t, n = median_cell(data, land & band[:, None])
    sites['lake_0'] = dict(kind='lake margin', j=j, i=i, lat=float(lat[j]), lon=float(data['lon'][i]),
                           wind_share=SHELTER, cells_in_band=n, ten_year_surface_c=round(t, 3))
    return sites


# ---- the Sun ---------------------------------------------------------------------------------------------------------
def window_mean_cosine(lat_deg, t0_h, period_h, n_windows):
    """Each window's mean of max(cos zenith, 0) on half-hour steps, for a Sun at declination zero."""
    t = (np.arange(n_windows * WINDOW_STEPS) + 0.5) * 0.5
    c = math.cos(math.radians(lat_deg)) * np.cos(2.0 * math.pi * (t - t0_h) / period_h)
    return np.maximum(c, 0.0).reshape(n_windows, WINDOW_STEPS).mean(1)


def fit_sun(czen, lat_deg):
    """Noon time (hours from the first window's start) and period fitted to the windows' mean cosine."""
    from scipy.optimize import least_squares
    n = czen.size
    period0 = SYNODIC_H
    # a first noon from the brightest of the first lunar day's windows, then a scan of the phase at that period
    grid = np.arange(0.0, period0, 1.0)
    cost = [np.sum((window_mean_cosine(lat_deg, g, period0, min(n, 40)) - czen[:40]) ** 2) for g in grid]
    t0 = float(grid[int(np.argmin(cost))])
    fit = least_squares(lambda p: window_mean_cosine(lat_deg, p[0], p[1], n) - czen, x0=[t0, period0],
                        x_scale=[10.0, 0.1], xtol=1e-12, ftol=1e-12)
    resid = fit.fun
    return dict(noon_h=float(fit.x[0]), period_h=float(fit.x[1]), rms=float(np.sqrt(np.mean(resid ** 2))),
                max_abs=float(np.max(np.abs(resid))))


def elevation_at(lat_deg, t_h, sun):
    hour = (360.0 * (t_h - sun['noon_h']) / sun['period_h'] + 180.0) % 360.0 - 180.0
    return optics.sun_elevation(lat_deg, hour, 0.0), hour


# ---- the light ---------------------------------------------------------------------------------------------------------
def clear_sky():
    """The clear sky's shortwave on the water (W/m2) and photosynthetic photons just below a flat surface
    (umol m-2 s-1) on the light study's elevation grid."""
    light = UL.sky('moon')
    shortwave = (light['direct'] + light['diffuse']).sum(1)
    _, _, _, salinity = air_constants()
    open_sea = dict(chlorophyll_mg_m3=0.1, cdom_440=0.02, fines_g_m3=0.0, salinity_g_kg=salinity)
    column = optics.Column(light['wavelength'], np.array([0.0]), **open_sea)
    par0 = UL.table(light, column)[:, 0, 0]
    return dict(elevation=UL.ELEVATION.copy(), shortwave=shortwave, par0=par0)


def run_clear_sky(clear_down_end, elevation_end, bin_deg=2.0):
    """The run's clear-sky downwelling shortwave at the surface (W/m2) against the Sun's elevation: the medians of its
    end-of-window values in bins of the Sun's elevation at each window's end. Returns the bins' centres and medians."""
    edges = np.arange(-90.0, 90.0 + 0.5 * bin_deg, bin_deg)
    index = np.clip(np.digitize(elevation_end, edges) - 1, 0, edges.size - 2)
    centres = 0.5 * (edges[:-1] + edges[1:])
    filled = np.unique(index)
    return centres[filled], np.array([np.median(clear_down_end[index == b]) for b in filled])


# ---- the water ---------------------------------------------------------------------------------------------------------
def air_constants():
    """The design air's gas constant and heat capacity, its O2 mole fraction and the held seas' salinity, as the
    waters' oxygen study holds them; the product records the four values it read (design_air)."""
    s = json.loads(LAKE_OXYGEN.read_text())['setting']
    return float(s['air_gas_constant']), float(s['air_cp']), float(s['o2_fraction']), float(s['salinity_g_kg'])


def site_forcing(data, site, sun, sky, wind_share=None):
    """The lake model's step forcing at a site and the step series of light and the Sun; wind_share replaces the
    site's own share of the cell's 10 m wind."""
    rd, _, o2_fraction, _ = air_constants()
    j, i = site['j'], site['i']
    n_w = data['tas'].shape[0]
    tas = data['tas'][:, j, i].astype(float)
    alb = data['alb'][:, j, i].astype(float)
    sw_all = np.where(alb < 0.999, data['rss'][:, j, i] / np.maximum(1.0 - alb, 1e-3), 0.0)
    lw = data['rls'][:, j, i] + STEFAN_BOLTZMANN * data['ts'][:, j, i].astype(float) ** 4
    q_open = data['hus_low'][:, j, i].astype(float)
    q = q_open.copy()
    q[:-1] = 0.5 * (q_open[:-1] + q_open[1:])        # the output's humidity is the value opening each window
    p = data['ps'][:, j, i].astype(float) * 100.0
    wind = np.hypot(data['uas'][:, j, i], data['vas'][:, j, i]).astype(float) * (
        site['wind_share'] if wind_share is None else wind_share)
    qs = np.array([L.saturation_humidity(a, pp) for a, pp in zip(tas, p)])
    q = np.minimum(q, qs)
    dt = WINDOW_H * 3600.0 / STEPS_PER_WINDOW
    centres = (np.arange(n_w) + 0.5) * WINDOW_H * 3600.0
    t = (np.arange(n_w * STEPS_PER_WINDOW) + 0.5) * dt
    f = {k: np.interp(t, centres, v) for k, v in dict(lw=lw, air=tas, q=q, u=wind, p=p).items()}
    # the Sun at every step; the run's clear sky from its end-of-window values; each window's clouds and sunlight
    elevation, hour = elevation_at(site['lat'], t / 3600.0, sun)
    elevation_end, _ = elevation_at(site['lat'], (np.arange(n_w) + 1.0) * WINDOW_H, sun)
    rsscs = data['rsscs'][:, j, i].astype(float)
    clear_end = np.where(alb < 0.999, rsscs / np.maximum(1.0 - alb, 1e-3), 0.0)
    grid_e, grid_f = run_clear_sky(clear_end, elevation_end)
    clear_sw = np.interp(elevation, grid_e, grid_f)
    clear_par = np.interp(elevation, sky['elevation'], sky['par0'])
    w = np.repeat(np.arange(n_w), STEPS_PER_WINDOW)
    mean_clear = clear_sw.reshape(n_w, STEPS_PER_WINDOW).mean(1)
    lit = mean_clear > CLEAR_FLOOR_W
    cloud_all = float(sw_all[lit].sum() / mean_clear[lit].sum())
    ratio = np.where(lit, sw_all / np.maximum(mean_clear, 1e-9), cloud_all)
    cloud = np.where(lit, np.clip(ratio, 0.0, CLOUD_CAP), cloud_all)
    sw = np.maximum(clear_sw * ratio[w], 0.0)
    par = np.maximum(clear_par * cloud[w], 0.0)
    forcing = L.Forcing(dt=dt, sw_down=sw, lw_down=f['lw'], air_k=f['air'], q_air=f['q'], wind10=f['u'],
                        pressure=f['p'], o2_fraction=o2_fraction)
    sun_end = np.maximum(np.sin(np.radians(elevation_end)), 0.0)
    solved = np.interp(grid_e, sky['elevation'], sky['shortwave'])
    high = grid_e >= 10.0
    checks = dict(window_shortwave_kept=float(np.max(np.abs(sw.reshape(n_w, -1).mean(1)[lit] - sw_all[lit]))),
                  clear_sky_follows_window_end=float(np.corrcoef(clear_end, sun_end)[0, 1]),
                  clear_sky_follows_window_mean=float(np.corrcoef(clear_end, data['czen'][:, j, i].astype(float))[0, 1]),
                  clear_sky_after_sunset_w=float(np.median(clear_end[(elevation_end < -1.5) & (elevation_end > -12.0)])),
                  run_clear_over_solved_channels=[float((grid_f[high] / solved[high]).min()),
                                                  float((grid_f[high] / solved[high]).max())],
                  cloud_factor_ten_years=cloud_all,
                  lit_windows_above_cap=float(np.mean(ratio[lit] > CLOUD_CAP)))
    return forcing, dict(t_h=t / 3600.0, elevation=elevation, hour=hour, par0=par, cloud=cloud[w]), checks


def run_water(forcing, depth, start_c):
    rd, cp, _, salinity = air_constants()
    basin = L.Basin.pan(depth, LAYER_M[depth], 1e4)
    scenario = L.Scenario('temperature', gpp=1.6, attenuation=1.7 / 3.0, sediment=0.3, deep_water=0.3)
    lake = L.Lake(basin, L.Physics(deep_mixing='wind', salinity=salinity), [scenario],
                  L.Air(rd, cp, MOON_SURFACE_GRAVITY))
    state = L.initial_state(basin, [scenario], start_c)
    return lake.run(forcing, state)


def transfer(forcing, out):
    """The O2 transfer velocity (m/s) at each step: the surface-renewal form on the run's water friction velocity and
    its surface's cooling buoyancy flux, and Cole and Caraco's wind form."""
    rd, cp, _, _ = air_constants()
    air = L.Air(rd, cp, MOON_SURFACE_GRAVITY)
    ts = np.asarray(out['surface_t'], float)
    renewal = np.zeros(ts.size)
    cole = np.zeros(ts.size)
    for n in range(ts.size):
        sens, latent, _, _ = L.bulk_fluxes(float(ts[n]), float(forcing.air_k[n]), float(forcing.q_air[n]),
                                           float(forcing.wind10[n]), float(forcing.pressure[n]), air)
        loss = L.WATER_EMISSIVITY * (STEFAN_BOLTZMANN * (ts[n] + 273.15) ** 4 - forcing.lw_down[n]) + sens + latent
        b = MOON_SURFACE_GRAVITY * L._expansion(float(ts[n])) * loss / L.RHO_CP if loss > 0.0 else 0.0
        renewal[n] = L.transfer_velocity('renewal', float(forcing.wind10[n]), float(out['u_star'][n]), b, float(ts[n]))
        cole[n] = L.transfer_velocity('cole_caraco', float(forcing.wind10[n]), 0.0, 0.0, float(ts[n]))
    return renewal, cole


def open_sea_c(data, site, n_steps, dt):
    """The run's sea surface at a sea cell (C), linear in time between the windows' centres."""
    ts = data['ts'][:, site['j'], site['i']].astype(float) - 273.15
    centres = (np.arange(ts.size) + 0.5) * WINDOW_H * 3600.0
    return np.interp((np.arange(n_steps) + 0.5) * dt, centres, ts)


def slab(forcing, depth, exchange_per_day, open_c, start_c):
    """A well-mixed slab of water trading a share of its water with the open water each day: its temperature (C) and
    O2 transfer velocity (m/s) at every step."""
    rd, cp, _, _ = air_constants()
    air = L.Air(rd, cp, MOON_SURFACE_GRAVITY)
    ph = L.Physics()
    dt = float(forcing.dt)
    keep = math.exp(-exchange_per_day * dt / 86400.0)
    n = len(forcing)
    temp, k = np.zeros(n), np.zeros(n)
    t = float(start_c)
    for i in range(n):
        sens, latent, _, tau = L.bulk_fluxes(t, float(forcing.air_k[i]), float(forcing.q_air[i]),
                                             float(forcing.wind10[i]), float(forcing.pressure[i]), air)
        emitted = L.WATER_EMISSIVITY * (STEFAN_BOLTZMANN * (t + 273.15) ** 4 - forcing.lw_down[i])
        net = (1.0 - ph.albedo) * forcing.sw_down[i] - emitted - sens - latent + ph.heat_flow
        loss = emitted + sens + latent
        b = MOON_SURFACE_GRAVITY * L._expansion(t) * loss / L.RHO_CP if loss > 0.0 else 0.0
        k[i] = L.transfer_velocity('renewal', float(forcing.wind10[i]), math.sqrt(tau / L.RHO0), b, t)
        t = t + net * dt / (L.RHO_CP * depth)
        t = open_c[i] + (t - open_c[i]) * keep
        temp[i] = t
    return temp, k


def summarise_water(tw, k, dt_h):
    """Mean, spread, extremes and hours a lunar day above each threshold of a water's temperature, and its O2
    transfer velocity."""
    lunar_days = tw.size * dt_h / SYNODIC_H
    return dict(
        water_mean_c=float(tw.mean()), water_p01_c=float(np.percentile(tw, 1)), water_p99_c=float(np.percentile(tw, 99)),
        water_min_c=float(tw.min()), water_max_c=float(tw.max()),
        hours_per_lunar_day_above={f'{x:g}': float((tw > x).sum() * dt_h / lunar_days) for x in THRESHOLDS_C},
        k_m_day_mean=float(k.mean() * 86400.0), k_m_day_p10=float(np.percentile(k, 10) * 86400.0))


def summarise_site(common, dt_h):
    lunar_days = common['par0'].size * dt_h / SYNODIC_H
    par = common['par0']
    return dict(wind10_mean=float(common['wind10'].mean()), open_water_mean_c=float(common['open_c'].mean()),
                open_water_min_c=float(common['open_c'].min()), open_water_max_c=float(common['open_c'].max()),
                par0_day_mean_umol=float(par.mean()), par0_mol_m2_d=float(par.mean() * 0.0864),
                dark_hours_per_lunar_day=float((par < PAR_DARK).sum() * dt_h / lunar_days),
                cloud_ratio_mean=float(common['cloud'].mean()))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--sites', nargs='*', default=None, help='site names to run (default: all)')
    args = ap.parse_args(argv)
    started = time.time()
    data = read_years()
    sites = choose_sites(data)
    sky = clear_sky()
    BULK.mkdir(parents=True, exist_ok=True)
    record_sites, files = {}, {}
    for name, site in sites.items():
        if args.sites and name not in args.sites:
            continue
        czen = data['czen'][:, site['j'], site['i']].astype(float)
        sun = fit_sun(czen, site['lat'])
        forcing, light, checks = site_forcing(data, site, sun, sky)
        start = SPIN_WINDOWS * STEPS_PER_WINDOW
        dt_h = forcing.dt / 3600.0
        start_c = float(forcing.air_k[0] - 273.15)
        series, cases = {}, {}
        for depth in DEPTHS_M:
            out = run_water(forcing, depth, start_c)
            renewal, cole = transfer(forcing, out)
            tw = 0.5 * (np.asarray(out['surface_t']) + np.asarray(out['bottom_t']))
            series[f'pan_{depth:g}m_t'] = tw
            series[f'pan_{depth:g}m_surface_t'] = np.asarray(out['surface_t'])
            series[f'pan_{depth:g}m_k'] = renewal
            series[f'pan_{depth:g}m_k_cole_caraco'] = cole
            cases[f'pan_{depth:g}m'] = summarise_water(tw[start:], renewal[start:], dt_h) | dict(
                surface_max_c=float(np.max(out['surface_t'][start:])),
                k_cole_caraco_m_day_mean=float(np.mean(cole[start:]) * 86400.0))
        if site['kind'] == 'sea shore':
            open_c = open_sea_c(data, site, len(forcing), forcing.dt)
        else:
            open_forcing = site_forcing(data, site, sun, sky, wind_share=1.0)[0]
            open_c = np.asarray(run_water(open_forcing, OPEN_LAKE_M, start_c)['surface_t'], float)
        for depth in DEPTHS_M:
            for share in EXCHANGE_PER_DAY:
                tw, k = slab(forcing, depth, share, open_c, start_c)
                key = f'slab_{depth:g}m_x{share:g}'
                series[key + '_t'], series[key + '_k'] = tw, k
                cases[key] = summarise_water(tw[start:], k[start:], dt_h)
        common = dict(t_h=light['t_h'], elevation=light['elevation'], par0=light['par0'], cloud=light['cloud'],
                      wind10=np.asarray(forcing.wind10), air_c=np.asarray(forcing.air_k) - 273.15,
                      sw_down=np.asarray(forcing.sw_down), pressure_pa=np.asarray(forcing.pressure), open_c=open_c)
        arrays = {k: np.asarray(v[start:], dtype=np.float64 if k == 't_h' else np.float32)
                  for k, v in (common | series).items()}
        path = BULK / f'{name}.npz'
        np.savez_compressed(path, dt_s=np.float64(forcing.dt), **arrays)
        files[name] = dict(path=str(path), sha256=sha256(path))
        record_sites[name] = site | dict(sun=sun, light_checks=checks,
                                         setting=summarise_site({k: arrays[k] for k in common}, dt_h), cases=cases)
        print(name, json.dumps(record_sites[name]['setting']), flush=True)
        for key, c in cases.items():
            print(' ', key, round(c['water_mean_c'], 2), round(c['water_p01_c'], 2), round(c['water_p99_c'], 2),
                  round(c['water_max_c'], 2), round(c['k_m_day_mean'], 2), flush=True)
    files_used = ['research/studies/littoral_community/forcing.py', 'biosphere/ecology/lake.py',
                  'research/studies/underwater_light/run.py', 'research/studies/underwater_light/optics.py',
                  'illumination/water_column/model.py']
    record = dict(
        schema=SCHEMA,
        producer=dict(lane='research', runner='research/studies/littoral_community/forcing.py',
                      files={f: sha256(ROOT / f) for f in files_used}, constants=constants_used(files_used)),
        evidence=('Model results on the corrected design climate (the run named in settings, years 10-19), carried '
                  'into the project\'s lake model on shallow pans and into exchanging slabs at five sites, with the '
                  'solved clear sky\'s light scaled by the run\'s cloud factor by window. A first physical setting '
                  'for shore life at T21: the run\'s 3-day means hold no weather within a window, the pans stand '
                  'closed and the slabs\' exchange with the open water is a stated range, the lake site\'s wind is '
                  'its land cell\'s, and the climate model and CM1 differ on the air near the ground, a gap the '
                  'project carries as a range.'),
        reading_rule=('Times in hours from the first window kept (after the spin-up). A pan\'s water temperature is '
                      'the mean of its surface and floor layers; a slab is one mixed layer trading the stated share of '
                      'its water with the open water each day. par0 is photosynthetic photons just below the surface, '
                      'umol m-2 s-1. k is the O2 transfer velocity, m/s in the series and m/day in the summary.'),
        settings=dict(run=RUN, years=list(YEARS), window_h=WINDOW_H, steps_per_window=STEPS_PER_WINDOW,
                      spin_windows=SPIN_WINDOWS, depths_m=list(DEPTHS_M), layers_m={f'{k:g}': v for k, v in LAYER_M.items()},
                      shelter=SHELTER, clear_floor_w=CLEAR_FLOOR_W, cloud_cap=CLOUD_CAP, par_dark=PAR_DARK, bands=BANDS,
                      lake_band=LAKE_BAND, exchange_per_day=list(EXCHANGE_PER_DAY), open_lake_m=OPEN_LAKE_M),
        sites=record_sites,
        series=files,
        inputs={str(p.relative_to(RUNS)): sha256(p) for p in year_files()}
        | {str(UL.SKY['moon'].relative_to(ROOT)): sha256(UL.SKY['moon'])},
        design_air=dict(zip(('air_gas_constant', 'air_cp', 'o2_fraction', 'salinity_g_kg'), air_constants())),
        seconds=round(time.time() - started, 1))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(OUT, record['seconds'], 's')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
