"""Why the settled Near-side Sea box runs colder aloft than the GCM: each model's heat budget at the site by PlaSim's
layers over a lunar day, both longwave codes on both columns, and each model's air against the moist adiabat of its
own air near the ground.

    climate/gcm/.venv/bin/python -m climate.crm.box_upper_air    # ../results/crm/box_upper_air_near_side_sea.json

- **CM1** (box_near_side_sea_day3, its third lunar day, 236 diagnostic outputs 3 hours apart): box_budget's terms
  (radiation split into longwave and sunlight, all-sky and clear-sky as CM1 writes them; microphysics; turbulence;
  resolved transport; the mean vertical wind on the box's own air; the day-night table, reckoned from its table as the
  patch applies it; nudging; damping), each the temperature tendency at constant pressure, and the storage between
  the day's first and last outputs. Each term is summed over each of PlaSim's layers by pressure, the box's levels
  shared out by their overlap, in W/m2.
- **The GCM** (its six sea cells at the site): the rh 0.7 branch over years 10-14 and Kuo's A28R_301k_cold_moist over
  years 0-5, the run that forced the box. Its heating diagnostics (codes 22-29) are single steps at the end of each
  3-day mean, which fall evenly over the lunar day at the site, so their plain mean is the lunar day's. What the
  physics' sum (code 21) holds beyond the listed processes is mkrain's evaporation of falling rain. The circulation is
  the storage, the trend of the yearly means, less the physics' sum.
- **The longwave codes on one column**: CM1's RRTMG longwave, as compiled for rrtmg_edge.py from the lunar build's own
  objects (its module's source is the box build's byte for byte), and PlaSim's (plasim_longwave.lwr with the run's
  continuum, its fitted correction and its curved-air escape), each on the box's lunar-day mean column and on the
  GCM's, clear sky. The GCM's column is laid on the box's layers by log pressure for RRTMG, and the box's on PlaSim's
  levels for PlaSim. CM1's own curved-air escape, which its build adds after RRTMG, is reckoned from the box's
  terluna_twilight.txt at each column's temperatures.
- **Sunlight**: PlaSim's shortwave copy with water vapour alone (plasim_shortwave.net_flux, the run's absorption and
  Rayleigh settings) on both columns, a daily mean over six Gauss-Legendre hour angles at 5 S under the box's
  1,108 W/m2, to see how the column alone moves the sunlight absorbed aloft.
- **The moist adiabats**: each model's lunar-day mean air on the box's levels against the undiluted pseudo-adiabat of
  its own air near the ground (the box's lowest 0.5 km mixed, as box_onset's parcel; the GCM's lowest level, as its
  schemes lift it), and the box's air against zero-buoyancy plumes, the columns that plumes entraining eps per metre
  into air of a given relative humidity hold at their own temperature, saturated moist static energy falling as
  eps (1 - rh) L q* (the zero-buoyancy plume of Singh and O'Gorman 2013, not read here).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from climate.crm import box_budget as bb
from climate.crm import box_onset as bo
from climate.crm import cm1_output as co
from climate.crm import ring_analysis as ra

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
OUT = RESULTS / 'box_upper_air_near_side_sea.json'
SCHEMA = 'terluna.climate.crm-box-upper-air/1'
DAY = 86400.0
CASE, SPAN = 'box_near_side_sea_day3', (59.0, 88.5)
START_CASE = 'box_near_side_sea'                                  # its first output is the GCM's mean air
GCMS = (('A28R_301k_cold_moist', (0, 5), 'Kuo, the run that forced the box'),
        ('A28R_301k_cold_sbm_rh07', (10, 14), 'simplified Betts-Miller at rh 0.7'))
TERMS = ('dt_lw', 'dt_sw', 'dt_conv', 'dt_lsp', 'dt_vdif', 'dt_sfc', 'dt_dca', 'dt_nudge', 'dt_diab')
CM1_TERMS = ('longwave', 'sunlight', 'longwave_clear', 'sunlight_clear', 'microphysics', 'turbulence',
             'resolved_transport', 'mean_vertical_wind', 'day_night_table', 'nudging', 'damping', 'storage')
BANDS = (('0-14 km', 7, 10), ('14-35 km', 5, 7), ('35-48 km', 4, 5), ('48-63 km', 3, 4), ('63-110 km', 1, 3))
SEA_K = 302.0355                                                  # the box's sea (its case's ground_k)
CO2, O2 = 4.0e-4, 0.174579                                        # the lunar build's (its TERLUNA flags)
SUNLIGHT_W_M2 = 1108.02                                           # the box's var17
HEIGHTS_M = (18900.0, 29800.0, 41800.0, 55800.0)
YEAR_S = 17160 * 1800.0                                           # the runs' year, 17,160 half-hour steps
EPS_PLUMES = (0.0, 2.0e-5, 5.0e-5, 1.0e-4, 2.0e-4)


def gravity() -> float:
    return json.loads((ra.RUNS / CASE / 'case.json').read_text())['build']['gravity_m_s2']


# ------------------------------------------------------------------------------------------------- the two budgets

def cm1_budget(case: str = CASE, span=SPAN) -> dict:
    """The box's lunar-day mean of each heating term (K/s of temperature at constant pressure) by level, with its mean
    temperature, vapour, pressure and density, and the storage between its first and last outputs."""
    path0 = ra.RUNS / case
    geo = ra.case_geometry(path0)
    zh, zw = geo['zh'], geo['zw']
    text = (path0 / 'namelist.input').read_text()
    table = bb.read_day_night(path0 / 'terluna_lsadv.txt')
    var18, var19 = bb.namelist_value(text, 'var18'), bb.namelist_value(text, 'var19')
    eq = co.eqtset(path0)
    layout = ra.grads_layout(path0 / 'cm1out_diag_s.ctl')[1]
    acc, first, last = {}, None, None
    for path in sorted(path0.glob('cm1out_diag_*_s.dat')):
        d = bb.read_diag(path0, path, layout)
        if not span[0] * DAY + 1.0 < d['mtime'] <= span[1] * DAY + 1.0:
            continue
        ex = (d['prs'] / co.P00) ** (co.RD / co.CP)
        row = dict(longwave=co.applied_rate(d['lwten'], eq), sunlight=co.applied_rate(d['swten'], eq),
                   longwave_clear=co.applied_rate(d['lwtenc'], eq), sunlight_clear=co.applied_rate(d['swtenc'], eq),
                   microphysics=d['ptb_mp'] * ex, turbulence=(d['ptb_pbl'] + d['ptb_diss']) * ex,
                   resolved_transport=sum(d[k] for k in ('ptb_hadv', 'ptb_vadv', 'ptb_div', 'ptb_hidiff',
                                                         'ptb_vidiff')) * ex,
                   mean_vertical_wind=d['ptb_lsw'] * ex,
                   day_night_table=bb.day_night_heating(*table, zh, d['mtime'], var18, var19) * ex,
                   nudging=d['ptb_nudge'] * ex, damping=d['ptb_rdamp'] * ex, radiation=d['ptb_rad'] * ex,
                   t=d['t'], qv=d['qv'], p=d['prs'], rho=d['rho'])
        for k, v in row.items():
            acc.setdefault(k, []).append(np.asarray(v, float))
        first = first or (d['mtime'], np.asarray(d['t'], float))
        last = (d['mtime'], np.asarray(d['t'], float))
    out = {k: np.mean(v, axis=0) for k, v in acc.items()}
    out['storage'] = (last[1] - first[1]) / (last[0] - first[0])
    out.update(zh=zh, zw=zw, outputs=len(acc['t']), span_days=[first[0] / DAY, last[0] / DAY])
    return out


def gcm_budget(run: str, years, box: str = START_CASE) -> dict:
    """The GCM's lunar-day mean at the box's cells: each heating diagnostic in W/m2 by layer (top first), the rain's
    evaporation as the physics' sum less the listed processes, the storage from the yearly means' trend and the
    circulation as the storage less the physics; with its mean air."""
    import netCDF4
    from climate.crm import gcm_site as gs
    from climate.gcm import sbm_offline as so
    from climate.gcm import year_files as yf
    c = so.run_constants(run)
    record = json.loads((ra.RUNS / box / 'case.json').read_text())
    cells = [tuple(x) for x in record['gcm_site']['cells']]
    acc = {k: [] for k in TERMS + ('ta', 'hus', 'ps', 'ts', 'w')}
    yearly = []
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(HERE.parent / 'gcm' / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            idx = [(int(np.argmin(abs(lat - a))), int(np.argmin(abs(lon - b)))) for a, b in cells]
            f = yf.fields(d, TERMS + ('ta', 'hus', 'ps', 'ts'))       # the humidity's means: gcm/year_files.py
        for j, i in idx:
            for k in TERMS + ('ta', 'hus'):
                acc[k].append(f[k][:, :, j, i])
            acc['ps'].append(f['ps'][:, j, i] * 100.0)
            acc['ts'].append(f['ts'][:, j, i])
            acc['w'].append(np.full(f['ps'].shape[0], np.cos(np.radians(lat[j]))))
        yearly.append(np.mean([f['ta'][:, :, j, i].mean(axis=0) for j, i in idx], axis=0))
    a = {k: np.concatenate(v, axis=0) for k, v in acc.items()}
    w = a['w']
    mean = {k: np.average(a[k], axis=0, weights=w) for k in TERMS + ('ta', 'hus', 'ps', 'ts')}
    edges = gs.plasim_edges(sigma)
    mass = mean['ps'] * np.diff(edges) / c['ga']
    cp = c['acpd'] * (1.0 + c['ADV'] * mean['hus'])
    wm2 = {k: cp * mean[k] * mass for k in TERMS}
    wm2['rain_evaporation'] = wm2['dt_diab'] - sum(wm2[k] for k in TERMS[:-1])
    trend = np.polyfit(np.arange(years[0], years[1] + 1), np.array(yearly), 1)[0] / YEAR_S if len(yearly) > 1 else 0.0
    wm2['storage'] = cp * trend * mass
    wm2['circulation'] = wm2['storage'] - wm2['dt_diab']
    return dict(run=run, years=list(years), sigma=sigma, edges=edges, mean=mean, wm2=wm2, mass=mass, cp=cp,
                z=gs.level_heights(mean['ta'], mean['hus'], sigma, c['ga']), cell_means=int(a['ta'].shape[0]),
                constants=c)


def interfaces(p, zh, zw):
    """Pressures (Pa) at the box's w levels from its scalar levels', log-linear, the ground's carried below."""
    from climate.gcm import sbm_offline as so
    p = np.asarray(p, float)
    return so.w_level_pressures(p, float(np.exp(np.log(p[0]) + (np.log(p[0]) - np.log(p[1])) * (zh[0] - zw[0])
                                                / (zh[1] - zh[0]))), zh, zw)


def on_layers(values, pw, edges, per_level: bool = False, g: float | None = None) -> np.ndarray:
    """Values on the box's levels summed over each of PlaSim's layers (top first) by their overlap in pressure: a rate
    (K/s) as cp times it times the overlap's mass (W/m2), or with per_level a quantity already per level (W/m2)."""
    pw = np.asarray(pw, float)
    g = gravity() if g is None else g
    thick = pw[:-1] - pw[1:]
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        top, bot = a * pw[0], b * pw[0]
        overlap = np.clip(np.minimum(pw[:-1], bot) - np.maximum(pw[1:], top), 0.0, None)
        out.append(float(np.sum(values * overlap / thick)) if per_level
                   else float(np.sum(values * co.CP * overlap / g)))
    return np.array(out)


def layer_means(values, pw, edges) -> np.ndarray:
    """The pressure-weighted mean of a box profile over each of PlaSim's layers."""
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        top, bot = a * pw[0], b * pw[0]
        overlap = np.clip(np.minimum(pw[:-1], bot) - np.maximum(pw[1:], top), 0.0, None)
        out.append(float(np.sum(np.asarray(values) * overlap) / max(overlap.sum(), 1e-12)))
    return np.array(out)


def bands(per_layer) -> dict:
    """PlaSim's layers grouped (top-first indices [a, b))."""
    return {name: float(np.sum(np.asarray(per_layer)[a:b])) for name, a, b in BANDS}


# ------------------------------------------------------------------------------------------------- the radiation

def twilight_shells(case: str = CASE):
    """The box's curved-air table: shell edges (m), each shell's extra infrared escape (W/m2) and the temperatures
    it is given at."""
    rows = [r for r in (ra.RUNS / case / 'terluna_twilight.txt').read_text().split('\n') if r.strip()]
    ntw, nsh = (int(v) for v in rows[0].split())
    return (np.array(rows[1].split(), float), np.array(rows[2 + ntw].split(), float)[:nsh],
            np.array(rows[3 + ntw].split(), float)[:nsh])


def cm1_escape(t, zw, shells) -> np.ndarray:
    """CM1's extra infrared escape (W/m2, positive out) of each level at temperatures t, as its build adds it after
    RRTMG: each shell's share over the levels it overlaps, scaled with the levels' temperature to the fourth power."""
    tz, lw, tt = shells
    out = np.zeros(len(t))
    for k in range(len(t)):
        for s in range(lw.size):
            ov = min(zw[k + 1], tz[s + 1]) - max(zw[k], tz[s])
            if ov > 0.0:
                out[k] += lw[s] * ov / (tz[s + 1] - tz[s]) * (t[k] / tt[s]) ** 4
    return out


def rrtmg_longwave(columns, g: float) -> list:
    """CM1's RRTMG longwave, clear sky, from rrtmg_edge's driver on (interface p, interface t, layer p, layer t,
    layer vapour mixing ratio, ground temperature) columns: heating (K/day) by layer, ground first."""
    from climate.crm import ring_radiation as rr
    from climate.crm import rrtmg_edge as re_
    exe = re_.WORK / 'as_built' / 'lw_columns'
    if not exe.exists():
        exe = re_.build('as_built')
    nlay = columns[0][2].size
    fmt = lambda v: ' '.join(f'{float(x):.10e}' for x in v)
    rows = [f'{len(columns)} {nlay} {co.CP:.10e} {g:.10e}']
    for pi, ti, pl, tl, ql, ts in columns:
        x = rr.mole_fraction(np.asarray(ql), 0.0289644)
        rows += [fmt([ts]), fmt(pi / 100.0), fmt(ti), fmt(pl / 100.0), fmt(tl), fmt(x / (1.0 - x)),
                 fmt([CO2] * nlay), fmt([O2] * nlay)]
    link = re_.WORK / 'RRTMG_LW_DATA'
    if not link.exists():
        link.symlink_to(re_.rrtmg_data())
    out = subprocess.run([str(exe)], input='\n'.join(rows) + '\n', capture_output=True, text=True, check=True,
                         cwd=re_.WORK).stdout.split('\n')
    return [np.array([float(v) for v in out[3 * i].split()]) for i in range(len(columns))]


def rrtmg_column(p_layer, t_layer, q_layer, zh, zw, psfc, tsfc):
    """A column as CM1 hands it to RRTMG (ring_radiation.column_levels: the box's interfaces, the ground at its own
    temperature, six isothermal layers above the top to 1 Pa): interfaces, layer pressures, temperatures, vapour."""
    from climate.crm import ring_radiation as rr
    p, t, q = rr.column_levels(zw, zh, p_layer, t_layer, q_layer, psfc, t_ground=tsfc)
    pl = np.exp(0.5 * (np.log(p[1:]) + np.log(p[:-1])))
    pl[:len(p_layer)] = p_layer
    extra = p.size - 1 - len(t_layer)
    return (p, t, pl, np.concatenate([t_layer, np.full(extra, t[-1])]),
            np.concatenate([q_layer, np.full(extra, q[-1])]), tsfc)


def plasim_longwave(run: str, sigma, edges, ps, t, q, t_ground, escape: bool = True) -> np.ndarray:
    """PlaSim's clear-sky longwave (plasim_longwave.lwr with the run's settings and correction, and its curved-air
    escape unless escape is False) on a column on its levels: W/m2 by layer, top first."""
    from climate.gcm import plasim_longwave as pl
    s = pl.run_settings(run)
    twlw, twlt = pl.twilight_terms(run) if escape else (None, None)
    _, heat, _ = pl.lwr(sigma, ps, np.append(t, t_ground), q, 0.0, s['gravity'], s['acpd'], s['co2'], s['th2oc'],
                        twlw=twlw, twlt=twlt, correction=s['correction'])
    return s['acpd'] * (1.0 + pl.ADV * np.asarray(q)) * heat * ps * np.diff(edges) / s['gravity']


def plasim_water_sunlight(sigma, ps, t, q, g: float, lat_deg: float = -5.0) -> np.ndarray:
    """PlaSim's clear-sky shortwave copy with water vapour alone (the runs' absorption 1.292 in 'scale' form and
    Rayleigh scale 2.237, ocean albedo) on a column: the daily mean absorbed in each layer (W/m2), the Sun over the
    equator and SUNLIGHT_W_M2 at the top, six Gauss-Legendre hour angles over the sunlit half."""
    from climate.gcm import plasim_shortwave as psw
    col = dict(sigma=sigma, ps=ps, t=t, q=q, albedo=0.069, albedo_direct=psw.ocean_direct_albedo,
               rcoeff=psw.RCOEFF_PER_RAYSCALE * 2.237, solar1=0.502127834309164, gravity=g)
    x, w = np.polynomial.legendre.leggauss(6)
    hours, weights = 0.25 * np.pi * (x + 1.0), 0.25 * np.pi * w
    absorbed = np.zeros(len(sigma))
    for mu, wt in zip(np.cos(np.radians(lat_deg)) * np.cos(hours), weights):
        net = psw.net_flux(col, mu, 'scale', 1.292)
        absorbed += wt * mu * (net[:-1] - net[1:])
    return SUNLIGHT_W_M2 * absorbed / np.pi


# ------------------------------------------------------------------------------------------------- the adiabats

def saturated_temperature(hs, z, p, g, guess):
    """The temperature (K) of saturated air of moist static energy hs at height z and pressure p."""
    t = guess
    for _ in range(40):
        qs = bo.saturation(t, p)
        step = (co.CP * t + g * z + bo.LV * qs - hs) / (co.CP + bo.LV * (bo.saturation(t + 0.01, p) - qs) / 0.01)
        t -= step
        if abs(step) < 1e-7:
            break
    return t


def pseudo_adiabat(zh, p, t0, q0, z0, g) -> np.ndarray:
    """The undiluted pseudo-adiabat (K) of air t0, q0 at height z0, keeping its moist static energy, on levels zh."""
    h, qt = co.CP * t0 + g * z0 + bo.LV * q0, q0
    out = np.zeros(len(zh))
    for k, (z, pk) in enumerate(zip(zh, p)):
        out[k], q, _ = bo.saturate(h, qt, z, pk, g)
        qt = min(qt, q)
    return out


def zero_buoyancy(zh, p, t0, q0, z0, g, eps, rh) -> np.ndarray:
    """The column (K) a plume entraining eps per metre into air of relative humidity rh holds at its own temperature:
    undiluted until it saturates, then its saturated moist static energy falling as eps (1 - rh) L q* with height."""
    h = co.CP * t0 + g * z0 + bo.LV * q0
    out = np.zeros(len(zh))
    wet = False
    hs = h
    for k in range(len(zh)):
        if not wet:
            out[k], _, wet = bo.saturate(h, q0, zh[k], p[k], g)
            continue
        hs -= eps * (1.0 - rh[k - 1]) * bo.LV * bo.saturation(out[k - 1], p[k - 1]) * (zh[k] - zh[k - 1])
        out[k] = saturated_temperature(hs, zh[k], p[k], g, out[k - 1])
    return out


# ------------------------------------------------------------------------------------------------- the analysis

def analyse() -> dict:
    from climate.gcm import sbm_offline as so
    g = gravity()
    b = cm1_budget()
    zh, zw = b['zh'], b['zw']
    pw = interfaces(b['p'], zh, zw)
    shells = twilight_shells()
    gcms = [gcm_budget(run, years) for run, years, _ in GCMS]
    edges, sigma = gcms[0]['edges'], gcms[0]['sigma']
    box = {k: on_layers(b[k], pw, edges, g=g) for k in CM1_TERMS}
    box['escape'] = -on_layers(cm1_escape(b['t'], zw, shells), pw, edges, per_level=True)
    box['cloud_longwave'] = box['longwave'] - box['longwave_clear'] - box['escape']
    box['convection'] = box['microphysics'] + box['resolved_transport'] + box['turbulence']
    box['forcing'] = box['mean_vertical_wind'] + box['day_night_table']
    box['rest'] = box['storage'] - sum(box[k] for k in ('longwave', 'sunlight', 'microphysics', 'turbulence',
                                                        'resolved_transport', 'forcing', 'nudging', 'damping'))
    t_box, q_box = layer_means(b['t'], pw, edges), layer_means(b['qv'], pw, edges)
    # the forcing on the box's start air, the GCM's own: the table's lunar-day mean and the mean vertical wind on it
    case = ra.RUNS / CASE
    z_t, th_t = bb.read_day_night(case / 'terluna_lsadv.txt')
    table = np.interp(zh, z_t, th_t.mean(axis=0))
    wls = np.loadtxt(case / 'terluna_wls.txt', skiprows=1)
    w = np.interp(zh, wls[:, 0], wls[:, 1])
    start = ra.RUNS / START_CASE
    d0 = bb.read_diag(start, sorted(start.glob('cm1out_diag_*_s.dat'))[0], ra.grads_layout(start /
                                                                                         'cm1out_diag_s.ctl')[1])
    ex0 = (d0['prs'] / co.P00) ** (co.RD / co.CP)
    forcing_start = on_layers((table - w * np.gradient(d0['th'], zh)) * ex0, pw, edges, g=g)
    ex = (b['p'] / co.P00) ** (co.RD / co.CP)
    forcing_now = on_layers((table - w * np.gradient(b['t'] / ex, zh)) * ex, pw, edges, g=g)
    # the longwave codes on both columns, and the sunlight copy
    cm1_col = rrtmg_column(b['p'], b['t'], b['qv'], zh, zw, pw[0], SEA_K)
    cols, plasim, sun = [cm1_col], {}, {}
    cm1_on_plasim = so.on_levels(dict(t=b['t'], q=b['qv'] / (1.0 + b['qv']), p=b['p'], ps=pw[0]), 'level', sigma,
                                 edges)
    for gm in gcms:
        m = gm['mean']
        lnp = np.log(sigma * m['ps'])[::-1]
        pc = b['p'] * m['ps'] / pw[0]
        tg = np.interp(-np.log(pc), -lnp, m['ta'][::-1])
        qg = np.interp(-np.log(pc), -lnp, (m['hus'] / (1.0 - m['hus']))[::-1])
        cols.append(rrtmg_column(pc, tg, qg, zh, zw, m['ps'], float(m['ts'])))
        gm['on_box_levels'] = dict(t=tg, q=qg, p=pc)
        for escape in (True, False):
            plasim[(gm['run'], 'own', escape)] = plasim_longwave(gm['run'], sigma, edges, m['ps'], m['ta'], m['hus'],
                                                                 float(m['ts']), escape)
            plasim[(gm['run'], 'box', escape)] = plasim_longwave(gm['run'], sigma, edges, pw[0], cm1_on_plasim[0],
                                                                 cm1_on_plasim[1], SEA_K, escape)
        sun[gm['run']] = plasim_water_sunlight(sigma, m['ps'], m['ta'], m['hus'], g)
    sun['box'] = plasim_water_sunlight(sigma, pw[0], cm1_on_plasim[0], cm1_on_plasim[1], g)
    heats = rrtmg_longwave(cols, g)
    rrtmg = []
    for (pi, *_), hr in zip(cols, heats):
        energy = hr / DAY * co.CP * (pi[:-1] - pi[1:]) / g
        rrtmg.append(on_layers(energy, pi, edges, per_level=True))
    rrtmg_box_own_check = dict(
        z_m=[float(zh[k]) for k in [int(np.argmin(abs(zh - z))) for z in HEIGHTS_M]],
        offline_on_mean_column_k_day=[float(heats[0][int(np.argmin(abs(zh - z)))]) for z in HEIGHTS_M],
        box_clear_mean_k_day=[float(b['longwave_clear'][int(np.argmin(abs(zh - z)))] * DAY) for z in HEIGHTS_M])
    # the moist adiabats
    mass = b['rho'] * np.diff(zw)
    t0, q0 = bo.source(b['p'], b['t'], b['qv'], mass * bo.layer_share(zw, 500.0), 500.0)
    parcel_box = pseudo_adiabat(zh, b['p'], t0, q0, zh[0], g)
    rh_box = b['qv'] / bo.saturation(b['t'], b['p'])
    at = lambda a: [float(np.interp(z, zh, a)) for z in HEIGHTS_M]
    adiabats = dict(heights_m=list(HEIGHTS_M), box=dict(
        source=dict(t_k=float(t0), q_g_kg=float(q0 * 1e3), h_kj_kg=float((co.CP * t0 + g * zh[0] + bo.LV * q0) / 1e3)),
        air_k=at(b['t']), parcel_k=at(parcel_box), air_minus_parcel_k=at(b['t'] - parcel_box),
        relative_humidity=at(rh_box)), gcms={})
    for gm in gcms:
        m, c = gm['mean'], gm['constants']
        tl, ql, zl = float(m['ta'][-1]), float(m['hus'][-1] / (1.0 - m['hus'][-1])), float(gm['z'][-1])
        lev = gm['on_box_levels']
        parcel = pseudo_adiabat(zh, lev['p'], tl, ql, zl, g)
        adiabats['gcms'][gm['run']] = dict(
            source=dict(t_k=tl, q_g_kg=ql * 1e3, height_m=zl, h_kj_kg=(co.CP * tl + g * zl + bo.LV * ql) / 1e3),
            air_k=at(lev['t']), parcel_k=at(parcel), air_minus_parcel_k=at(lev['t'] - parcel),
            parcel_minus_box_parcel_k=at(parcel - parcel_box), air_minus_box_air_k=at(lev['t'] - b['t']))
    adiabats['zero_buoyancy'] = {
        f'{eps:g}/m': dict(eps_per_efolding=float(eps * co.RD * 270.0 / g),
                           box_air_minus_plume_box_rh_k=at(b['t'] - zero_buoyancy(zh, b['p'], t0, q0, zh[0], g, eps,
                                                                                    rh_box)),
                           box_air_minus_plume_rh07_k=at(b['t'] - zero_buoyancy(zh, b['p'], t0, q0, zh[0], g, eps,
                                                                                np.full(zh.size, 0.7))))
        for eps in EPS_PLUMES}
    # the record
    layer = lambda a: [float(x) for x in a]
    gcm_out = {}
    for gm, r_own in zip(gcms, rrtmg[1:]):
        W = gm['wm2']
        conv = W['dt_conv'] + W['dt_lsp'] + W['rain_evaporation'] + W['dt_vdif'] + W['dt_dca']
        gcm_out[gm['run']] = dict(
            years=gm['years'], cell_means=gm['cell_means'], height_m=layer(gm['z']), t_k=layer(gm['mean']['ta']),
            vapour_g_kg=layer(gm['mean']['hus'] * 1e3), surface_k=float(gm['mean']['ts']),
            w_m2={k: layer(W[k]) for k in ('dt_lw', 'dt_sw', 'dt_conv', 'dt_lsp', 'rain_evaporation', 'dt_vdif',
                                           'dt_sfc', 'dt_dca', 'dt_nudge', 'dt_diab', 'storage', 'circulation')},
            convection_w_m2=layer(conv), radiation_w_m2=layer(W['dt_lw'] + W['dt_sw']),
            bands_w_m2=dict(radiation=bands(W['dt_lw'] + W['dt_sw']), longwave=bands(W['dt_lw']),
                            sunlight=bands(W['dt_sw']), convection=bands(conv), circulation=bands(W['circulation'])),
            plasim_clear_longwave_w_m2=dict(own_column=layer(plasim[(gm['run'], 'own', True)]),
                                            own_column_without_escape=layer(plasim[(gm['run'], 'own', False)]),
                                            box_column=layer(plasim[(gm['run'], 'box', True)]),
                                            box_column_without_escape=layer(plasim[(gm['run'], 'box', False)])),
            rrtmg_clear_longwave_w_m2=dict(own_column=layer(r_own),
                                           cm1_escape_own_column=layer(-on_layers(
                                               cm1_escape(gm['on_box_levels']['t'], zw, shells), pw, edges,
                                               per_level=True))),
            plasim_water_sunlight_w_m2=dict(own_column=layer(sun[gm['run']])))
    return dict(
        schema=SCHEMA, case=CASE, span_days=b['span_days'], outputs=b['outputs'], sigma=layer(sigma),
        layers=[f'L{k}' for k in range(len(sigma))], bands={n: [f'L{k}' for k in range(a, c)] for n, a, c in BANDS},
        box=dict(t_k=layer(t_box), vapour_g_kg=layer(q_box * 1e3),
                 w_m2={k: layer(v) for k, v in box.items()},
                 bands_w_m2={k: bands(box[k]) for k in ('longwave', 'sunlight', 'convection', 'forcing', 'storage',
                                                        'longwave_clear', 'cloud_longwave', 'escape')},
                 forcing_on_start_air_w_m2=layer(forcing_start), forcing_on_day3_air_w_m2=layer(forcing_now),
                 rrtmg_clear_longwave_w_m2=layer(rrtmg[0]),
                 plasim_water_sunlight_w_m2=layer(sun['box']), rrtmg_offline_check=rrtmg_box_own_check),
        gcms=gcm_out, adiabats=adiabats, evidence=EVIDENCE, reading_rule=READING_RULE)


EVIDENCE = ('Diagnosis from existing output and offline radiation, no new runs. The box: one realization of its third '
            'lunar day (236 domain-mean outputs). The GCM: single-step heating snapshots spread evenly over the lunar '
            'day, five or six years at six sea cells, with the circulation as the remainder. The longwave codes run '
            'clear sky on lunar-day mean columns, so their heating differs from the mean of each model\'s own heating '
            '(nonlinearity, clouds); the sunlight copy holds water vapour alone and serves for the change between '
            'columns only. The zero-buoyancy plumes are a calculation on the box\'s mean sounding.')
READING_RULE = ('Every per-layer list follows PlaSim\'s layers top first (L0 at the top to L9 at the ground; sigma '
                'gives their full levels); bands group them (L9-L7 0-14 km, L6-L5 14-35, L4 35-48, L3 48-63, L2-L1 '
                '63-110 km at the site). W/m2 are the lunar-day mean heating of each layer\'s air. box.w_m2: CM1\'s '
                'terms laid on those layers by pressure; convection is microphysics plus resolved transport plus '
                'turbulence; forcing the mean vertical wind on the box\'s own air plus the day-night table; escape the '
                'curved air\'s extra infrared; cloud_longwave the all-sky less clear-sky longwave less the escape. '
                'gcms[run].w_m2: PlaSim\'s diagnostics by name, rain_evaporation the physics\' sum less the listed '
                'processes, circulation the storage less the physics; convection_w_m2 is convection, large-scale '
                'condensation, rain evaporation, vertical diffusion and dry adjustment together. The radiation '
                'entries are clear sky on lunar-day mean columns: rrtmg_* CM1\'s RRTMG without its escape (which '
                'cm1_escape_own_column gives), plasim_clear_longwave PlaSim\'s with its escape unless marked '
                'without; box_column is the box\'s air on PlaSim\'s levels. adiabats: temperatures at heights_m on the '
                'box\'s levels, each model\'s air against the undiluted pseudo-adiabat of its own air near the ground, '
                'and the box\'s air against zero-buoyancy plumes at each eps per metre (eps_per_efolding at 270 K).')


def compact(x, digits: int = 5):
    if isinstance(x, dict):
        return {str(k): compact(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [compact(v, digits) for v in x]
    if isinstance(x, (float, np.floating)):
        return None if not np.isfinite(x) else float(f'{float(x):.{digits}g}')
    if isinstance(x, np.integer):
        return int(x)
    return x


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args(argv)
    out = compact(analyse())
    out['producer'] = dict(domain='climate', files={
        f'{p.parent.name}/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
        for p in (Path(__file__).resolve(), HERE / 'box_budget.py', HERE / 'rrtmg_edge.py',
                  HERE / 'rrtmg_offline' / 'lw_columns.f90', HERE.parent / 'gcm' / 'plasim_longwave.py',
                  HERE.parent / 'gcm' / 'plasim_shortwave.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT.relative_to(HERE.parents[1])}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
