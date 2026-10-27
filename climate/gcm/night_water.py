"""The sea's water by process and height, quarter by quarter of the lunar day: what supplies the vapour the GCM
condenses over the Near-side Sea after dark, where it condenses and how much of it falls, against CM1's sea box.

    climate/gcm/.venv/bin/python -m climate.gcm.night_water        # under a minute: results/gcm/night_water_l15.json

The GCM rains about 2 mm/day after midnight at the Near-side Sea, where CM1's sea box rains 0.14. With the Betts-Miller
scheme's deep convection held off at night (the diagnostic run of 27 October) the rain stays, and large-scale
condensation forms it. This module follows the water through each layer at the site, in three GCM runs and in the box:
- kuo_l10, the Kuo run on ten levels whose budget at the site forced the box (A28R_301k_cold_moist, years 0-5);
- sbm_l15, the Betts-Miller run on fifteen levels with the depth test and the shallow cap
  (A28R_301k_cold_sbm_rh07_deep200_l15cap, years 0-1);
- sbm_l15_night_off, the same with deep convection held off wherever the Sun is down, a diagnostic that writes the
  vapour of every process (A28R_301k_cold_sbm_rh07_deep200_l15cap_diagnostic_night, years 0-1);
- sbm_l15_night_off_cloud, the same with condensate.py's cloud water and the humidity's true mean
  (A28R_301k_cold_sbm_rh07_deep200_l15cap_diagnostic_night_cloud, years 0-1), which writes its cloud's water and
  processes too;
- CM1's sea box over its settled third lunar day (box_near_side_sea_day3, days 59.125-88.5), forced by kuo_l10.

In the GCM, at the box's six sea cells. PlaSim writes each process's tendency as a single time step at the end of each
3-day mean (pbl_mixing.snapshot_hours), so the physics of a quarter is the mean of the steps that fall in it. Code 40,
which the runner names dq_dyn, is the physics' moistening as the spectral step applies it (spectrald adds the diabatic
tendencies and mkdqtgp writes them back to the grid), and code 41 the horizontal diffusion; code 40's global means by
layer match the physics' sum (code 31), and code 41 sums to nothing over the globe. The storage of each layer's vapour
and heat over a quarter comes from the 3-day means, fitted with the day's mean and four harmonics, each divided by the
damping its window gives it (crm/gcm_site.py), the humidity's means from year_files (each record of hus holds the
humidity at its interval's start, so each interval's mean is the average of the record and the next, and each year
file's first and last records are left out), and the circulation is the storage less the physics as applied: the
transport sideways and up or down with the horizontal diffusion, the heat's with the adiabatic warming of sinking air.
The column's physics over a quarter is exactly its evaporation less its rain, and its storage the change of its
precipitable water, from the means fitted the same way; the snapshots' column sum beside the physics measures their
noise, and the layers' summed storage beside the precipitable water's checks the humidity's means. Large-scale
condensation is also read from its heating (mklsp heats by ALV over the moist specific heat), which the runs without the
vapour by process need. The rain's evaporation (mkrain) is written apart only in the diagnostic run; elsewhere it stays
in the physics' sum.

In CM1, its domain diagnostics every three hours. Its vapour budget closes, layer by layer, as the prescribed
circulation (the day-night table, terluna_lsadv.txt, and the imposed mean vertical wind's advection, qvb_lsw; the
nudging is nil below its start), the resolved transport (the resolved and implicit-diffusion flux terms, qvb_vturbr
and qvb_vturbd), the boundary-layer scheme with the surface's evaporation (qvb_pbl) and the microphysics (qvb_mp). The
microphysics splits by CM1's own temperature tendencies into condensation, evaporation from cloud and from rain,
deposition and sublimation (tt_*), each taken to vapour with the latent heat of vaporisation or sublimation; their sum
matches qvb_mp within a few percent. The heat budget is grouped the same way, with the radiation (ptb_rad) and the
dissipation. The storage over a quarter is the change of the domain-mean vapour, and heat, between the quarter's ends.

The layers are the fifteen-level grid's (levels.py's l15) at the site: each GCM layer's heights by quarter from the
fifteen-level run's means. CM1's levels go into them by the share of each level's depth inside each layer, and the
ten-level run's layers by the share of each layer's height; both keep every total.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-night-water/1'
OUT = RESULTS / 'night_water_l15.json'
DAY = 86400.0
WATER_DENSITY = 1000.0
BOX = 'box_near_side_sea'                    # the case whose gcm_site cells are the site
CM1_CASE = 'box_near_side_sea_day3'
RUNS = (('kuo_l10', 'A28R_301k_cold_moist', (0, 5)),
        ('sbm_l15', 'A28R_301k_cold_sbm_rh07_deep200_l15cap', (0, 1)),
        ('sbm_l15_night_off', 'A28R_301k_cold_sbm_rh07_deep200_l15cap_diagnostic_night', (0, 1)),
        ('sbm_l15_night_off_cloud', 'A28R_301k_cold_sbm_rh07_deep200_l15cap_diagnostic_night_cloud', (0, 1)))
REFERENCE = 'sbm_l15'                        # the run whose layers are the bands
QUARTERS = (('sunrise-noon', -90.0, 0.0), ('noon-sunset', 0.0, 90.0), ('sunset-midnight', 90.0, 180.0),
            ('midnight-sunrise', -180.0, -90.0))
HARMONICS = 4
MEANS = ('ta', 'hus', 'ps', 'czen', 'pr', 'prc', 'evap', 'prw', 'cl')
SNAP_WATER = ('dq_diab', 'dq_fill', 'dq_dyn', 'dq_dyn_late', 'dq_sfc', 'dq_vdif', 'dq_lsp', 'dq_conv', 'dq_dca',
              'dcw_cond', 'dq_cwevap', 'dcw_rain')                    # the last three with condensate.py's cloud
SNAP_HEAT = ('dt_diab', 'dt_nudge', 'dt_sfc', 'dt_vdif', 'dt_lw', 'dt_sw', 'dt_lsp', 'dt_conv', 'dt_dca', 'dt_evap')
LV_CM1, LS_CM1 = 2.501e6, 2.834e6            # J/kg, CM1's constants.F
CM1_WATER = {'table': None, 'imposed_ascent': ('qvb_lsw',), 'nudging': ('qvb_nudge',),
             'resolved': ('qvb_vturbr', 'qvb_vturbd'), 'boundary_layer': ('qvb_pbl',), 'microphysics': ('qvb_mp',)}
CM1_HEAT = {'table': None, 'imposed_ascent': ('ptb_lsw',), 'nudging': ('ptb_nudge',),
            'resolved': ('ptb_vturbr', 'ptb_vturbd'), 'boundary_layer': ('ptb_pbl',), 'microphysics': ('ptb_mp',),
            'radiation': ('ptb_rad',), 'dissipation': ('ptb_diss',), 'damping': ('ptb_rdamp',)}
CM1_MICRO = (('condensation', 'tt_cond', LV_CM1), ('evaporation_from_cloud', 'tt_evac', LV_CM1),
             ('evaporation_from_rain', 'tt_evar', LV_CM1), ('deposition', 'tt_dep', LS_CM1),
             ('sublimation', 'tt_subl', LS_CM1))


def quarter_index(hour_deg) -> np.ndarray:
    """The index in QUARTERS of each hour angle (degrees, 0 = noon)."""
    h = (np.asarray(hour_deg, float) + 180.0) % 360.0 - 180.0
    out = np.full(h.shape, -1)
    for k, (_, a, b) in enumerate(QUARTERS):
        out[(h >= a) & (h < b)] = k
    return out


def overlap_matrix(src_edges, dst_edges) -> np.ndarray:
    """The share of each source layer (rows) lying in each target layer (columns), by height: layers between
    consecutive edges (m, the ground up). Each row sums to the share of its layer the targets cover."""
    s, d = np.asarray(src_edges, float), np.asarray(dst_edges, float)
    lo, hi = np.maximum(s[:-1, None], d[None, :-1]), np.minimum(s[1:, None], d[None, 1:])
    return np.clip(hi - lo, 0.0, None) / np.diff(s)[:, None]


def quarter_storage(centre_deg, values, window_deg: float, day_s: float, harmonics: int = HARMONICS,
                    undamp: bool = True) -> np.ndarray:
    """The rate of change (per day) of a quantity over each quarter, from means centred at the hour angles centre_deg
    that cover window_deg each: the series fitted with its mean and harmonics, each undamped for the window (or as
    they stand, for single steps), and its change from the quarter's start to its end over the quarter's length."""
    from climate.crm import gcm_site as gs
    coef = gs.harmonic_coefficients(centre_deg, values, harmonics)
    if undamp:
        coef = gs.undamped(coef, window_deg, harmonics)
    out = []
    for _, a, b in QUARTERS:
        change = gs.evaluate(coef, np.array([b]), harmonics)[0] - gs.evaluate(coef, np.array([a]), harmonics)[0]
        out.append(change / ((b - a) / 360.0 * day_s / DAY))
    return np.asarray(out)


def quarter_average(centre_deg, values, window_deg: float, harmonics: int = HARMONICS, points: int = 181) -> np.ndarray:
    """The mean over each quarter of a quantity from means centred at centre_deg, fitted and undamped as above."""
    from climate.crm import gcm_site as gs
    coef = gs.undamped(gs.harmonic_coefficients(centre_deg, values, harmonics), window_deg, harmonics)
    return np.asarray([gs.evaluate(coef, np.linspace(a, b, points), harmonics).mean(axis=0) for _, a, b in QUARTERS])


def saturation(t, p, c) -> np.ndarray:
    """PlaSim's saturation specific humidity (mklsp's formula) at temperature t (K) and pressure p (Pa)."""
    q = c['rdbrv'] * c['ra1'] * np.exp(c['ra2'] * (t - c['TMELT']) / (t - c['ra4'])) / p
    q = np.minimum(q, c['rdbrv'])
    return q / (1.0 - (1.0 / c['rdbrv'] - 1.0) * q)


def gcm_site_water(run: str, years) -> dict:
    """A run's vapour and heat by process at the site's cells by quarter (layers top first): mm/day of water and W/m2
    in each layer, the storage and circulation from the means, the column's exact physics, and the mean state."""
    import netCDF4
    from climate.crm import gcm_site as gs
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import solar_day_s
    from climate.gcm import pbl_mixing as pm
    from climate.gcm import sbm_offline as so
    from climate.gcm import year_files as yf
    record = json.loads((ra.RUNS / BOX / 'case.json').read_text())
    cells = [tuple(cc) for cc in record['gcm_site']['cells']]
    lon_site = record['configuration']['site']['lon_deg']
    c = so.run_constants(run)
    g = c['ga']
    snaps = {k: [] for k in SNAP_WATER + ('dq_evap',) + SNAP_HEAT + ('lsp_from_heating', 'conv_from_heating',
                                                                      'cloud_water')}
    means = {k: [] for k in ('water', 'enthalpy', 't', 'q', 'rh', 'cl', 'z_half', 'e', 'p', 'pc', 'prw')}
    hours, hours_snap, present = [], [], set()
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(HERE / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            idx = [(int(np.argmin(abs(lat - a))), int(np.argmin(abs(lon - b)))) for a, b in cells]
            jj, ii = np.array([j for j, _ in idx]), np.array([i for _, i in idx])
            w = np.cos(np.radians(lat[jj]))
            w = w / w.sum()
            names = [k for k in MEANS + SNAP_WATER + SNAP_HEAT + ('dq_evap', 'qcw') if k in d.variables]
            present |= set(names)
            got = yf.fields(d, names)                                   # the humidity's means: year_files
            f = {k: v[..., jj, ii] for k, v in got.items() if k != 'czen'}   # (time, [lev,] cells)
            hour = (lon_site - gs.subsolar_longitude(got['czen'], lat, lon) + 180.0) % 360.0 - 180.0
        edges = gs.plasim_edges(sigma)
        dsig = np.diff(edges)
        mean = lambda a: np.tensordot(a, w, axes=(-1, 0))                     # the cells' weighted mean
        ps = f['ps'] * 100.0                                                   # (time, cells), Pa
        mass = ps[:, None, :] * dsig[None, :, None] / g                        # kg/m2 by layer
        cpq = c['acpd'] * (1.0 + c['ADV'] * f['hus'])
        to_mm = 86.4e6                                                         # m/s of water to mm/day
        for k in SNAP_WATER:
            if k in f:
                snaps[k].append(mean(f[k]) * to_mm)
        if 'dq_evap' in f:
            snaps['dq_evap'].append(mean(f['dq_evap'] * mass / WATER_DENSITY) * to_mm)
        for k in SNAP_HEAT:
            if k in f:
                snaps[k].append(mean(f[k] * cpq * mass))
        snaps['lsp_from_heating'].append(mean(f['dt_lsp'] * cpq * mass / c['ALV']) * DAY)
        if 'qcw' in f:
            snaps['cloud_water'].append(mean(f['qcw'] * mass))                # kg/m2 by layer, at each step
        snaps['conv_from_heating'].append(mean(f['dt_conv'] * cpq * mass / c['ALV']) * DAY)
        means['water'].append(mean(f['hus'] * mass))
        means['enthalpy'].append(mean(cpq * f['ta'] * mass))
        means['t'].append(mean(f['ta']))
        means['q'].append(mean(f['hus']))
        p_full = sigma[None, :, None] * ps[:, None, :]
        means['rh'].append(mean(f['hus'] / saturation(f['ta'], p_full, c)))
        means['cl'].append(mean(f['cl']))
        zh = pm.half_level_heights(np.moveaxis(f['ta'], 1, -1), np.moveaxis(f['hus'], 1, -1), None, edges, g,
                                   c['gascon'], c['rdbrv'])                    # (time, cells, nlev-1), top first
        means['z_half'].append(mean(np.moveaxis(zh, -1, 1)))
        means['e'].append(-mean(f['evap']) * to_mm)
        means['p'].append(mean(f['pr']) * to_mm)
        means['pc'].append(mean(f['prc']) * to_mm)
        means['prw'].append(mean(f['prw']))
        hours.append(hour)
        hours_snap.append(pm.snapshot_hours(hour))
    hour, hour_snap = np.concatenate(hours), np.concatenate(hours_snap)
    s = {k: np.concatenate(v) for k, v in snaps.items() if v}
    m = {k: np.concatenate(v) for k, v in means.items()}
    step = float(np.median(np.diff(np.degrees(np.unwrap(np.radians(hour))))))
    day_s = solar_day_s()
    qs, qm = quarter_index(hour_snap), quarter_index(hour)
    snap_q = {k: np.stack([v[qs == n].mean(axis=0) for n in range(len(QUARTERS))]) for k, v in s.items()}
    mean_q = {k: np.stack([v[qm == n].mean(axis=0) for n in range(len(QUARTERS))]) for k, v in m.items()}
    storage_w = quarter_storage(hour, m['water'], step, day_s)               # kg/m2/day = mm/day, by layer
    storage_h = quarter_storage(hour, m['enthalpy'], step, day_s) / DAY      # W/m2 by layer
    exact = {k: quarter_average(hour, m[k], step) for k in ('e', 'p', 'pc')}
    storage_col = quarter_storage(hour, m['prw'], step, day_s)
    if 'cloud_water' in s:                                   # the cloud's storage from its single steps
        storage_cloud = quarter_storage(hour_snap, s['cloud_water'].sum(axis=1), step, day_s, undamp=False)
    else:
        storage_cloud = np.zeros(len(QUARTERS))
    out = dict(run=run, years=list(years), cells=[list(cc) for cc in cells], sigma=sigma.tolist(),
               snapshots=[int((qs == n).sum()) for n in range(len(QUARTERS))],
               means=[int((qm == n).sum()) for n in range(len(QUARTERS))], window_deg=step,
               writes_vapour_by_process='dq_lsp' in present, quarters={})
    for n, (name, _, _) in enumerate(QUARTERS):
        sq = {k: v[n] for k, v in snap_q.items()}
        vap = dict(storage=storage_w[n], physics_applied=sq['dq_dyn'], physics=sq['dq_diab'],
                   circulation=storage_w[n] - sq['dq_dyn'], horizontal_diffusion=sq['dq_dyn_late'],
                   large_scale_condensation_from_heating=-sq['lsp_from_heating'],
                   convective_heating_as_water=sq['conv_from_heating'], fill=sq['dq_fill'])
        for k, label in (('dq_sfc', 'surface'), ('dq_vdif', 'diffusion'), ('dq_lsp', 'large_scale_condensation'),
                         ('dq_conv', 'convection'), ('dq_dca', 'dry_adjustment'), ('dq_evap', 'rain_evaporation'),
                         ('dcw_cond', 'to_cloud'), ('dq_cwevap', 'cloud_evaporation'),
                         ('dcw_rain', 'cloud_to_precipitation')):
            if k in sq:
                vap[label] = sq[k]
        if 'dq_lsp' not in sq:                                                  # all the physics but condensation
            vap['physics_less_large_scale_condensation'] = sq['dq_diab'] + sq['lsp_from_heating']
        heat = dict(storage=storage_h[n], physics=sq['dt_diab'], circulation=storage_h[n] - sq['dt_diab'],
                    longwave=sq['dt_lw'], shortwave=sq['dt_sw'], large_scale_condensation=sq['dt_lsp'],
                    convection=sq['dt_conv'], mixing=sq['dt_vdif'] + sq['dt_dca'] + sq['dt_sfc'], nudge=sq['dt_nudge'])
        listed = sum(sq[k] for k in SNAP_HEAT[1:] if k in sq)
        heat['rain_evaporation' if 'dt_evap' in sq else 'rain_evaporation_by_residual'] = (
            sq['dt_evap'] if 'dt_evap' in sq else sq['dt_diab'] - listed)
        mq = {k: v[n] for k, v in mean_q.items()}
        column = dict(evaporation=float(exact['e'][n]), rain=float(exact['p'][n]),
                      convective_rain=float(exact['pc'][n]), large_scale_rain=float(exact['p'][n] - exact['pc'][n]),
                      physics_exact=float(exact['e'][n] - exact['p'][n]),
                      physics_applied_from_snapshots=float(sq['dq_dyn'].sum()),
                      storage=float(storage_col[n]), storage_from_layers=float(storage_w[n].sum()),
                      cloud_storage=float(storage_cloud[n]),
                      circulation_exact=float(storage_col[n] + storage_cloud[n] - (exact['e'][n] - exact['p'][n])))
        out['quarters'][name] = dict(vapour_mm_day={k: np.asarray(v, float) for k, v in vap.items()},
                                     heat_w_m2={k: np.asarray(v, float) for k, v in heat.items()},
                                     column_mm_day=column,
                                     state=dict(t_k=mq['t'], q_g_kg=mq['q'] * 1000.0, rh=mq['rh'], cloud=mq['cl'],
                                                **({'cloud_water_kg_m2': sq['cloud_water']}
                                                   if 'cloud_water' in sq else {})),
                                     z_half_m=mq['z_half'])
    return out


def bands_from(reference: dict) -> dict:
    """The bands by quarter: the reference run's layer edges (m, the ground up) at the site."""
    out = {}
    for name, q in reference['quarters'].items():
        zh = np.asarray(q['z_half_m'], float)[::-1]                               # the ground up, inner half levels
        out[name] = np.r_[0.0, zh, 2.0 * zh[-1] - zh[-2]]                        # the top layer given its depth below
    return out


def on_bands(values_top_first, z_half_top_first, edges) -> np.ndarray:
    """Layer totals (top first, as PlaSim holds them) laid on the bands (the ground up) by the share of each layer's
    height in each band; the top layer is taken as deep as the one below it."""
    v = np.asarray(values_top_first, float)[..., ::-1]
    zh = np.asarray(z_half_top_first, float)[::-1]
    src = np.r_[0.0, zh, 2.0 * zh[-1] - zh[-2]]
    return np.tensordot(v, overlap_matrix(src, edges), axes=(-1, 0))


def gcm_on_bands(site: dict, edges: dict, reference: bool) -> dict:
    """A run's layers on the bands, the ground up: its own order reversed for the reference run's grid, and laid by
    height for another grid."""
    out = {k: v for k, v in site.items() if k != 'quarters'}
    out['quarters'] = {}
    for name, q in site['quarters'].items():
        lay = (lambda v: np.asarray(v, float)[::-1]) if reference else (
            lambda v: on_bands(v, q['z_half_m'], edges[name]))
        state = {k: (np.asarray(v, float)[::-1] if reference else None) for k, v in q['state'].items()}
        out['quarters'][name] = dict(vapour_mm_day={k: lay(v) for k, v in q['vapour_mm_day'].items()},
                                     heat_w_m2={k: lay(v) for k, v in q['heat_w_m2'].items()},
                                     column_mm_day=q['column_mm_day'],
                                     state=state if reference else dict(
                                         own_layers_ground_up={k: np.asarray(v, float)[::-1] for k, v in
                                                               q['state'].items()},
                                         own_half_levels_m=np.asarray(q['z_half_m'], float)[::-1]))
    return out


def cm1_site_water(edges: dict, case_name: str = CM1_CASE) -> dict:
    """CM1's box by quarter on the bands (the ground up): its vapour (mm/day) and heat (W/m2) by process, the
    microphysics split, the storage and what closes the budget, the column's evaporation and rain, and the mean state
    (vapour, relative humidity, cloud fraction and condensate)."""
    from climate.crm import box_budget as bb
    from climate.crm import cm1_output as co
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    case = ra.RUNS / case_name
    record = json.loads((case / 'case.json').read_text())
    span = json.loads((HERE.parent / 'results' / 'crm' / f'box_onset_{case_name}.json').read_text())['span_days']
    start, lon = record['configuration']['start_hour_angle_deg'], record['configuration']['site']['lon_deg']
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    dz = np.diff(zw)
    text = (case / 'namelist.input').read_text()
    z_q, table_q = bb.read_day_night(case / 'terluna_lsadv.txt', 'qv')
    z_t, table_t = bb.read_day_night(case / 'terluna_lsadv.txt', 'theta')
    var18, var19 = bb.namelist_value(text, 'var18'), bb.namelist_value(text, 'var19')
    layout = ra.grads_layout(case / 'cm1out_diag_s.ctl')[1]
    rows = []
    for path in sorted(case.glob('cm1out_diag_*_s.dat')):
        d = bb.read_diag(case, path, layout)
        if not span[0] * DAY - 1.0 <= d['mtime'] <= span[1] * DAY + 1.0:
            continue
        mass = d['rho'] * dz
        exner = (d['prs'] / co.P00) ** (co.RD / co.CP)
        vap = {k: sum(d[n] for n in names) * mass * DAY for k, names in CM1_WATER.items() if names}
        vap['table'] = bb.day_night_heating(z_q, table_q, zh, d['mtime'], var18, var19) * mass * DAY
        heat = {k: sum(d[n] for n in names) * exner * co.CP * mass for k, names in CM1_HEAT.items() if names}
        heat['table'] = bb.day_night_heating(z_t, table_t, zh, d['mtime'], var18, var19) * exner * co.CP * mass
        micro = {k: -d[n] * exner * co.CP / lat * mass * DAY for k, n, lat in CM1_MICRO}
        rows.append(dict(t=d['mtime'] / DAY, hour=local_hour_angle(start, lon, d['mtime'], solar_day_s()), vap=vap,
                         heat=heat, micro=micro, water=d['qv'] * mass,
                         enthalpy=co.CP * d['th'] * exner * mass, mass=mass,
                         q=d['qv'] * mass, rh=d['rh'] * mass, cloud=d['cldfra'] * mass,
                         cloud_water=d['qc'] * mass, ice=(d['qi'] + d['qs'] + d['qg']) * mass,
                         rain_water=d['qr'] * mass,
                         e=d['qfx'] * DAY, p=d['prate'] * DAY))
    t = np.array([r['t'] for r in rows])
    hour = np.array([r['hour'] for r in rows])
    qi = quarter_index(hour)
    out = dict(case=case_name, span_days=list(span), outputs=len(rows), forced_by=record['gcm_site']['run'],
               outputs_by_quarter=[int((qi == n).sum()) for n in range(len(QUARTERS))], quarters={})
    unwrapped = np.degrees(np.unwrap(np.radians(hour)))
    for n, (name, a, b) in enumerate(QUARTERS):
        sel = qi == n
        share = overlap_matrix(zw, edges[name])
        lay = lambda v: np.asarray(v, float) @ share
        avg = lambda key, sub=None: (np.mean([r[key][sub] if sub else r[key] for r, s in zip(rows, sel) if s], axis=0))
        # the quarter's ends: the hour angles a and b on the unwrapped series, the stored quantities linear between
        span_q = unwrapped[sel]
        lo = span_q.min() - ((span_q.min() - a) % 360.0)
        hi = lo + (b - a)
        ends = lambda key: [np.array([np.interp(x, unwrapped, [r[key][k] for r in rows])
                                      for k in range(len(zh))]) for x in (lo, hi)]
        length = (t[-1] - t[0]) / (unwrapped[-1] - unwrapped[0]) * (b - a)       # days
        w0, w1 = ends('water')
        h0, h1 = ends('enthalpy')
        vap = {k: lay(avg('vap', k)) for k in rows[0]['vap']}
        vap['circulation'] = vap['table'] + vap['imposed_ascent'] + vap['nudging']
        vap['storage'] = lay((w1 - w0) / length)
        vap['rest'] = vap['storage'] - sum(vap[k] for k in ('circulation', 'resolved', 'boundary_layer',
                                                             'microphysics'))
        heat = {k: lay(avg('heat', k)) for k in rows[0]['heat']}
        heat['circulation'] = heat['table'] + heat['imposed_ascent'] + heat['nudging']
        heat['storage'] = lay((h1 - h0) / (length * DAY))
        heat['rest'] = heat['storage'] - sum(heat[k] for k in ('circulation', 'resolved', 'boundary_layer',
                                                               'microphysics', 'radiation', 'dissipation', 'damping'))
        micro = {k: lay(avg('micro', k)) for k in rows[0]['micro']}
        mass = lay(avg('mass'))
        out['quarters'][name] = dict(
            vapour_mm_day=vap, heat_w_m2=heat, microphysics_mm_day=micro,
            column_mm_day=dict(evaporation=float(avg('e')), rain=float(avg('p')),
                               physics_exact=float(avg('e') - avg('p')),
                               storage=float(vap['storage'].sum()), circulation=float(vap['circulation'].sum()),
                               rest=float(vap['rest'].sum())),
            state=dict(q_g_kg=lay(avg('q')) / mass * 1000.0, rh=lay(avg('rh')) / mass, cloud=lay(avg('cloud')) / mass,
                       cloud_water_kg_m2=lay(avg('cloud_water')), ice_kg_m2=lay(avg('ice')),
                       rain_water_kg_m2=lay(avg('rain_water'))),
            quarter_days=float(length), ends_hour_deg=[float(lo), float(hi)])
    return out


def rounded(x, digits=3):
    if isinstance(x, dict):
        return {k: rounded(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rounded(v, digits) for v in x]
    if isinstance(x, np.ndarray):
        return [rounded(float(v), digits) for v in np.ravel(x)] if x.ndim else rounded(float(x), digits)
    if isinstance(x, float):
        return round(x, digits) if math.isfinite(x) else None
    return x


EVIDENCE = ('Diagnostics of existing runs; no new run. GCM: one realization each of six years of the Kuo run that '
            'forced the box and two years each of the fifteen-level runs, the first two of a spin-up from the '
            'fifteen-level restart; the physics by process are single time steps at the end of each 3-day mean, about '
            "60 a quarter in two years and 180 in six, averaged over the six sea cells, so a layer's value carries "
            "their noise, which the column's snapshot sum against its exact physics measures; the storage and the "
            "column's exact physics come from the 3-day means fitted with four harmonics. CM1: one realization of the "
            "sea box's settled third lunar day, its domain means every three hours (about 59 a quarter), each quarter "
            'sampled once; the box loses 0.37 mm/day of vapour over that day, more than half above 9 km. The kuo_l10 '
            'run is the one the box was forced by; the fifteen-level runs have their own circulation, which no box has'
            ' met.')
READING_RULE = ("bands: the fifteen-level grid's layers at the site, the ground up, edges_m by quarter from the "
                "sbm_l15 run's means; every profile is a list over them, the ground up. Vapour in mm/day of water "
                'gained by each layer (positive moistening), heat in W/m2 gained (positive warming). GCM: storage from'
                ' the fitted means; physics_applied (code 40) as the spectral step applies it, physics (code 31) at '
                'the grid points; circulation = storage - physics_applied (the transport with the horizontal '
                'diffusion, and for heat the adiabatic warming); the processes as PlaSim writes them where the run '
                'does (surface, diffusion, large_scale_condensation, convection, dry_adjustment, rain_evaporation), '
                "with large-scale condensation also from its heating; column_mm_day: the quarter's exact evaporation, "
                "rain and their difference from the means against the snapshots' column sum, its storage from the "
                "precipitable water against the layers' sum, and the circulation as that storage less the exact "
                "physics. The ten-level run's layers are laid on the bands by height, its state on its own layers. "
                "CM1: circulation = table + imposed_ascent + nudging (the GCM's forcing), resolved (its own convection"
                " and eddies), boundary_layer (the scheme, the surface's evaporation entering at the ground), "
                'microphysics (net), rest = storage less them; microphysics_mm_day: the gross processes as vapour '
                '(condensation and deposition negative, the evaporations and sublimation positive); state: vapour '
                '(g/kg), relative humidity and cloud fraction as domain and mass means, and the condensate held '
                '(kg/m2). Quarters by hour angle from noon. ranges: the same summed over the bands 1-4, 5-7 and 8-9 '
                'from the ground (about 0-5.2, 5.2-11.8 and 11.8-23.4 km); spread: those sums and the column for each '
                'year of the fifteen-level runs and each two-year block of the Kuo run, their differences the measure '
                "of the snapshots' noise.")


GROUPS = ((0, 4), (4, 7), (7, 9))            # the bands from the ground: about 0-5.2, 5.2-11.8 and 11.8-23.4 km
BLOCKS = {'kuo_l10': ((0, 1), (2, 3), (4, 5)), 'sbm_l15': ((0, 0), (1, 1)), 'sbm_l15_night_off': ((0, 0), (1, 1)),
          'sbm_l15_night_off_cloud': ((0, 0), (1, 1))}


def group_sums(quarter: dict, edges) -> dict:
    """A quarter's vapour and heat by process summed over the height ranges GROUPS, keyed by their edges in km."""
    out = {}
    for a, b in GROUPS:
        key = f'{edges[a] / 1000.0:.1f}-{edges[b] / 1000.0:.1f} km'
        out[key] = {part: {k: float(np.sum(np.asarray(v, float)[a:b])) for k, v in quarter[part].items()}
                    for part in ('vapour_mm_day', 'heat_w_m2', 'microphysics_mm_day') if part in quarter}
    return out


def analyse() -> dict:
    sites = {label: gcm_site_water(run, years) for label, run, years in RUNS}
    edges = bands_from(sites[REFERENCE])
    same = lambda site: site['sigma'] == sites[REFERENCE]['sigma']
    gcm = {label: gcm_on_bands(site, edges, same(site)) for label, site in sites.items()}
    cm1 = cm1_site_water(edges)
    ranges = dict(cm1={q: group_sums(cm1['quarters'][q], edges[q]) for q in cm1['quarters']})
    for label, run, _ in RUNS:
        ranges[label] = {q: group_sums(gcm[label]['quarters'][q], edges[q]) for q in gcm[label]['quarters']}
    spread = {}
    for label, run, _ in RUNS:
        spread[label] = []
        for years in BLOCKS[label]:
            site = gcm_on_bands(gcm_site_water(run, years), edges, same(sites[label]))
            spread[label].append(dict(years=list(years), quarters={
                q: dict(ranges=group_sums(v, edges[q]), column_mm_day=v['column_mm_day'])
                for q, v in site['quarters'].items()}))
    return dict(schema=SCHEMA, quarters=[q[0] for q in QUARTERS],
                bands=dict(reference=REFERENCE, edges_m={k: v for k, v in edges.items()}),
                gcm=gcm, cm1=cm1, ranges=ranges, spread=spread, evidence=EVIDENCE, reading_rule=READING_RULE)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args(argv)
    out = analyse()
    files = (Path(__file__), HERE / 'pbl_mixing.py', HERE / 'sbm_offline.py', HERE.parent / 'crm' / 'gcm_site.py',
             HERE.parent / 'crm' / 'box_budget.py')
    out['producer'] = dict(domain='climate', files={f'{p.parent.name}/{p.name}': hashlib.sha256(p.read_bytes())
                                                    .hexdigest()[:16] for p in files})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rounded(out), indent=1) + '\n')
    print('wrote', args.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
