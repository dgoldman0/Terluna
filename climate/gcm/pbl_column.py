"""PlaSim's boundary layer in a single column on fifteen levels, against the CM1 boxes (27 October 2026; the record is
the section "A nonlocal boundary layer for the Moon" in README.md beside this file).

    climate/gcm/.venv/bin/python -m climate.gcm.pbl_column     # ../results/gcm/pbl_column_l15.json, minutes on one core

Each column takes a CM1 box's air on PlaSim's fifteen levels (levels.py's l15) at the start of one quarter of the lunar
day and runs two quarters with the box's own forcing, interpolated between its 3-hourly outputs: its radiation, the
day-night table and the mean ascent acting on the column's own air, its surface sensible and latent heat and friction
velocity, and its domain-mean wind, which the column holds as given. The second quarter is compared with CM1's means over
the same hours. The three columns are the sea's evening (from noon), the sea's night (from sunset) and the land's
afternoon (from sunrise).

A step follows PlaSim's physics in order: the surface fluxes into the lowest layer and the vertical diffusion
(pbl_scheme.vdiff, fluxmod.f90's vdiff as a copy, with or without the nonlocal boundary layer); the box's radiation; the
mass-flux convection scheme (massflux.column, as in PlaSim from the state after those tendencies, its boundary-layer
forcing those tendencies with the large-scale terms) with PlaSim's evaporation of its falling rain at the box's constant
(massflux.evaporate); the dry adjustment (pbl_scheme.dca, rainmod.f90's mkdca as a copy); then the large-scale terms and a
saturation adjustment whose condensate rains out. Four schemes run from the same air under the same forcing, the 2x2 of
Blackadar's length (Earth's 160 m, or stretched by the gravities to 966 m) and the boundary layer (PlaSim's local
diffusion alone, or the nonlocal scheme beneath it). The step is PlaSim's 1800 s, the implicit diffusion taken over it.

Each column reports, over the compared quarter and beside CM1: the mixed layer's depth by pbl_mixing's rule (the first
height above the lowest full level where the mean virtual potential temperature has risen 0.5 K), the vapour at the
full levels near 0.3 and 6.2 km and its drop between them, the virtual potential temperature's rise from 6.2 to 8.2 km
(the stable layer that caps CM1's night cloud), the heat and vapour the diffusion and the dry adjustment carry through
the half levels, the nonlocal scheme's depth and velocity scales, how often the dry adjustment acts, the convection's
type, tops and rain, and the column's water and heat budgets, the spin-up quarter's beside the compared one's. As a
sensitivity (one_quarter) the same columns also start from the box's air at the compared quarter's own start.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from climate.gcm import massflux as mf
from climate.gcm import pbl_scheme as ps

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
OUT = RESULTS / 'pbl_column_l15.json'
SCHEMA = 'terluna.climate.gcm-pbl-column/1'
LEVEL_SET = 'l15'
STEP_S = 1800.0                                     # PlaSim's time step at T21 (deltsec)
DAY = 86400.0
# (name, site, CM1 case, span in days, the quarter the column starts in, the quarter compared)
CASES = (('sea evening', 'sea', 'box_near_side_sea_day3', (59.0, 88.5), 'afternoon', 'evening'),
         ('sea night', 'sea', 'box_near_side_sea_day3', (59.0, 88.5), 'evening', 'night'),
         ('land afternoon', 'land', 'box_0e_design', (29.5, 59.0), 'morning', 'afternoon'))
SCHEMES = (('earth_local', 'earth', False), ('lunar_local', 'lunar', False), ('earth_nonlocal', 'earth', True),
           ('lunar_nonlocal', 'lunar', True))
PBL_QUARTERS = {'morning': 'sunrise-noon', 'afternoon': 'noon-sunset', 'evening': 'sunset-midnight',
                'night': 'midnight-sunrise'}           # massflux's quarter names and pbl_mixing's
LOW_KM, MID_KM, LID_KM = 0.3, 6.2, 8.2                 # the vapour drop's levels and the lid's top
REPORT_HALF_KM = (0.6, 1.8, 3.2, 5.2, 7.2, 9.3)      # half levels whose fluxes are reported


def grid() -> np.ndarray:
    from climate.gcm import levels as lv
    return np.asarray(lv.half_levels(LEVEL_SET), float)


def box_rows(case: str, span) -> dict:
    """The box's outputs as massflux.cm1_outputs gives them, each with its ground-relative domain-mean wind (u, v, top
    first) from the same diagnostic file."""
    from climate.crm import box_budget as bb
    from climate.crm import ring_analysis as ra
    box = mf.cm1_outputs(case, span)
    folder = ra.RUNS / case
    layout = ra.grads_layout(folder / 'cm1out_diag_s.ctl')[1]
    winds = {}
    for path in sorted(folder.glob('cm1out_diag_*_s.dat')):
        d = bb.read_diag(folder, path, layout)
        if not span[0] * DAY + 1.0 < d['mtime'] <= span[1] * DAY + 1.0:
            continue
        winds[round(float(d['mtime']) / DAY, 6)] = ((d['u'] + d['umove'])[::-1], (d['v'] + d['vmove'])[::-1])
    for r in box['rows']:
        r['u'], r['v'] = winds[round(r['day'], 6)]
    return box


def gridded_rows(box: dict, half) -> list:
    """Each output laid on the grid by layer means in pressure (massflux.to_grid), with its winds, table and ascent."""
    out = []
    for r in box['rows']:
        gr = mf.to_grid(r, half)
        for key in ('table_dtdt', 'table_dqdt', 'w_mean', 'u', 'v'):
            gr[key] = mf.to_grid(dict(r, t=r[key], q=r[key]), half)['t']
        gr.update(day=r['day'], hour=r['hour'], quarter=r['quarter'], ust=r['ust'], g=r['g'])
        out.append(gr)
    return out


def constants(g: float) -> dict:
    """PlaSim's constants as the massflux tests take them, at the box's gravity."""
    return dict(mf.PLASIM, ga=g)


def interface_flux(dtdt, dqdt, mass, c):
    """The heat (W/m2, PlaSim's dry heat capacity) and vapour (W/m2 of latent heat) a conservative process carries
    upward through each interface, from its tendencies summed from the top."""
    heat = np.cumsum(c['acpd'] * mass * dtdt)[:-1]
    vapour = np.cumsum(c['ALV'] * mass * dqdt)[:-1]
    return heat, vapour


def step(t, q, rows, times, now, sigma, half, pf, ph, mass, g, c, land, tag, lamm, nonlocal_):
    """One step of the column at time now: returns the new temperature and vapour and the step's record."""
    akap = c['gascon'] / c['acpd']
    j = int(np.clip(np.searchsorted(times, now, side='right') - 1, 0, len(times) - 2))
    f = (now - times[j]) / (times[j + 1] - times[j])
    lerp = lambda key: (1.0 - f) * np.asarray(rows[j][key], float) + f * np.asarray(rows[j + 1][key], float)
    hour = rows[j]['hour'] + f * ((rows[j + 1]['hour'] - rows[j]['hour'] + 180.0) % 360.0 - 180.0)
    shf, lhf, ust = float(lerp('shf')), float(lerp('lhf')), float(lerp('ust'))
    u, v = lerp('u'), lerp('v')
    z_full, z_half = mf.heights(t, q, pf, ph, g, c)
    # the surface fluxes into the lowest layer, as mkshfl and mkevap put them there before vdiff
    t1, q1 = t.copy(), q.copy()
    t1[-1] += STEP_S * shf / (c['acpd'] * mass[-1])
    q1[-1] += STEP_S * lhf / (c['ALV'] * mass[-1])
    surface = None
    if nonlocal_:
        rho = pf[-1] / (c['gascon'] * t1[-1] * (1.0 + (1.0 / c['rdbrv'] - 1.0) * q1[-1]))
        surface = dict(wth0=shf / (rho * c['acpd']), wq0=lhf / (rho * c['ALV']), ustar=ust)
    bl = ps.vdiff(t1, q1, u, v, sigma, half, STEP_S, c, lamm=lamm, q_heat=q, surface=surface)
    t2 = t1 + STEP_S * bl['dtdt']
    q2 = q1 + STEP_S * bl['dqdt']
    # the box's radiation, then the large-scale terms PlaSim's dynamics stand for, on the step's starting air
    t3 = t2 + STEP_S * lerp('rad_dtdt')
    w = lerp('w_mean')
    th = t * (1.0e5 / pf) ** akap
    dyn_t = lerp('table_dtdt') - w * np.gradient(th, z_full) * (pf / 1.0e5) ** akap
    dyn_q = lerp('table_dqdt') - w * np.gradient(q, z_full)
    omega = -(pf / (c['gascon'] * t)) * g * w
    # convection from the state after the tendencies so far, its boundary-layer forcing those and the large-scale terms
    bl_t = (t3 - t) / STEP_S + dyn_t
    bl_q = (q2 - q) / STEP_S + dyn_q
    zf3, zh3 = mf.heights(t3, q2, pf, ph, g, c)
    cv = mf.column(t3, q2, pf, ph, g, STEP_S, shf, lhf, bl_t, bl_q, omega, np.hypot(u, v), land, c=c, z_full=zf3,
                   z_half=zh3, melt='plasim')
    ev = mf.evaporate(cv['column'], cv['net_rain'], cv['net_snow'], cv['dtdt'], cv['dqdt'], mf.GAMMA[tag])
    t4 = t3 + STEP_S * (cv['dtdt'] + ev['dtdt'])
    q4 = q2 + STEP_S * (cv['dqdt'] + ev['dqdt'])
    # the dry adjustment
    dtd, dqd, passes = ps.dca(t4, q4, sigma, half, akap, STEP_S)
    t5 = t4 + STEP_S * dtd
    q5 = q4 + STEP_S * dqd
    # the large-scale terms and the condensation of what they saturate
    t6 = t5 + STEP_S * dyn_t
    q6 = np.maximum(q5 + STEP_S * dyn_q, 1.0e-9)
    t7, q7, ls_rain = mf.condense_large_scale(t6, q6, pf, mass, c)
    heat_bl, vap_bl = interface_flux(bl['dtdt_diffusion'], bl['dqdt'], mass, c)
    heat_dca, vap_dca = interface_flux(dtd, dqd, mass, c)
    rain = ev['rain'] + ev['snow'] + float(ls_rain.sum()) / STEP_S
    # where the water and the heat go over the step (kg/m2/s and W/m2): the surface's evaporation and the large-scale
    # terms bring water, the rain takes it, the rest is stored; the air's dry enthalpy takes the surface's sensible heat,
    # the radiation, the large-scale terms and the latent heat of the net condensation (the vapour lost to rain)
    budget = dict(evaporation=lhf / c['ALV'], large_scale_water=float(np.sum(mass * dyn_q)), rain=rain,
                  water_stored=float(np.sum(mass * (q7 - q))) / STEP_S, sensible=shf,
                  radiation=float(np.sum(mass * c['acpd'] * lerp('rad_dtdt'))),
                  large_scale_heat=float(np.sum(mass * c['acpd'] * dyn_t)),
                  heat_stored=float(np.sum(mass * c['acpd'] * (t7 - t))) / STEP_S)
    rec = dict(hour=hour, quarter=mf.quarter_name(hour), heat_bl=heat_bl, vap_bl=vap_bl, heat_dca=heat_dca,
               vap_dca=vap_dca, dca_passes=passes, conv_type=cv['type'], conv_top_m=cv['z_top'],
               conv_base_m=cv['z_base'], rain_conv=ev['rain'] + ev['snow'], rain_ls=float(ls_rain.sum()) / STEP_S,
               conv_dqdt=cv['dqdt'] + ev['dqdt'], conv_dtdt=cv['dtdt'] + ev['dtdt'], vapour_kept=False,
               budget=budget)
    p = bl['pbl']
    if p is not None:
        rec.update(h=p['h'], wstar=p['wstar'], wm=p['wm'], we=p['entrainment_velocity'], unstable=p['unstable'],
                   vapour_kept=p['vapour_kept'], excess=p['excess'])
    return t7, q7, rec


def rho_half(t, q, pf, ph, c):
    """Air density at the interfaces (kg/m3), from vdiff's mean temperature there; for reporting only."""
    tv = t * (1.0 + (1.0 / c['rdbrv'] - 1.0) * q)
    return ph[1:-1] / (c['gascon'] * 0.5 * (tv[:-1] + tv[1:]))


def run_column(rows, start_quarter: str, compare_quarter: str, site: str, scheme: tuple, half) -> dict:
    """The column from the box's first output in start_quarter through its last in compare_quarter, under one scheme
    (name, length 'earth' or 'lunar', nonlocal). Returns the compared quarter's means and each quarter's rain and
    convection."""
    name, length, nonlocal_ = scheme
    first = next(i for i, r in enumerate(rows) if r['quarter'] == start_quarter)
    last = max(i for i, r in enumerate(rows) if r['quarter'] == compare_quarter)
    sub = rows[first:last + 1]
    g = sub[0]['g']
    c = constants(g)
    lamm = ps.LAMBDA_EARTH_M if length == 'earth' else ps.lunar_length(g)
    t, q = sub[0]['t'].copy(), sub[0]['q'].copy()
    pf, ph = sub[0]['p_full'], sub[0]['p_half']
    sigma, sigmah = pf / ph[-1], ph[1:] / ph[-1]
    mass = np.diff(ph) / g
    times = np.array([r['day'] for r in sub]) * DAY
    land = site == 'land'
    acc = {}
    now = times[0]
    while now < times[-1] - 1.0:
        t, q, rec = step(t, q, sub, times, now, sigma, sigmah, pf, ph, mass, g, c, land, site, lamm, nonlocal_)
        a = acc.setdefault(rec['quarter'], dict(t=[], q=[], rec=[]))
        a['t'].append(t)
        a['q'].append(q)
        a['rec'].append(rec)
        now += STEP_S
    akap = c['gascon'] / c['acpd']
    out = dict(scheme=name, mixing_length_m=lamm, nonlocal_=nonlocal_, quarters={})
    for quarter, a in acc.items():
        recs = a['rec']
        tm, qm = np.mean(a['t'], axis=0), np.mean(a['q'], axis=0)
        zf, zh = mf.heights(tm, qm, pf, ph, g, c)
        theta = tm * (1.0e5 / pf) ** akap
        thv = theta * (1.0 + (1.0 / c['rdbrv'] - 1.0) * qm)
        convect = [r for r in recs if r['conv_type'] is not None]
        tops = [r['conv_top_m'] for r in convect if r['conv_top_m'] is not None]
        qd = dict(steps=len(recs), z_full_m=zf, z_half_m=zh[1:-1], theta_v_k=thv, q_g_kg=qm / (1.0 - qm) * 1000.0,
                  t_k=tm, heat_bl_w_m2=np.mean([r['heat_bl'] for r in recs], axis=0),
                  vapour_bl_w_m2=np.mean([r['vap_bl'] for r in recs], axis=0),
                  heat_dca_w_m2=np.mean([r['heat_dca'] for r in recs], axis=0),
                  vapour_dca_w_m2=np.mean([r['vap_dca'] for r in recs], axis=0),
                  dca_share=float(np.mean([r['dca_passes'] > 0 for r in recs])),
                  dca_failed=int(sum(r['dca_passes'] < 0 for r in recs)),
                  rain_conv_mm_day=float(np.mean([r['rain_conv'] for r in recs]) * DAY),
                  rain_ls_mm_day=float(np.mean([r['rain_ls'] for r in recs]) * DAY),
                  convective_share=len(convect) / len(recs),
                  types={k: sum(r['conv_type'] == k for r in recs) / len(recs) for k in ('deep', 'shallow', 'mid')},
                  top_median_km=float(np.median(tops)) / 1000.0 if tops else None,
                  top_p90_km=float(np.percentile(tops, 90)) / 1000.0 if tops else None,
                  conv_moistening_mm_day=np.mean([r['conv_dqdt'] for r in recs], axis=0) * mass * DAY,
                  conv_heating_w_m2=np.mean([r['conv_dtdt'] for r in recs], axis=0) * mass * c['acpd'],
                  water_mm_day={k: float(np.mean([r['budget'][k] for r in recs])) * DAY
                                for k in ('evaporation', 'large_scale_water', 'rain', 'water_stored')},
                  heat_w_m2={k: float(np.mean([r['budget'][k] for r in recs]))
                             for k in ('sensible', 'radiation', 'large_scale_heat', 'heat_stored')})
        if nonlocal_:
            qd.update(h_median_km=float(np.median([r['h'] for r in recs])) / 1000.0,
                      h_mean_km=float(np.mean([r['h'] for r in recs])) / 1000.0,
                      wstar_median=float(np.median([r['wstar'] for r in recs])),
                      we_median_mm_s=float(np.median([r['we'] for r in recs])) * 1000.0,
                      unstable_share=float(np.mean([r['unstable'] for r in recs])),
                      excess_median_k=float(np.median([r['excess'] for r in recs])),
                      vapour_kept_steps=int(sum(r['vapour_kept'] for r in recs)))
        out['quarters'][quarter] = qd
    return out


def nearest(z, target_m):
    return int(np.argmin(np.abs(np.asarray(z, float) - target_m)))


def measures(z_full, thv, q_g_kg) -> dict:
    """The comparison's measures on a profile at the column's full levels (top first): the mixed layer's depth by
    pbl_mixing's rule from the lowest full level, the vapour at the levels near 0.3 and 6.2 km and its drop, and the
    virtual potential temperature's rise from 6.2 to 8.2 km."""
    from climate.gcm import pbl_mixing as pm
    zf = np.asarray(z_full, float)[::-1]
    lo, mid, lid = nearest(zf, LOW_KM * 1000.0), nearest(zf, MID_KM * 1000.0), nearest(zf, LID_KM * 1000.0)
    th, qq = np.asarray(thv, float)[::-1], np.asarray(q_g_kg, float)[::-1]
    return dict(mixed_depth_km=_km(pm.rise_height(zf, th, float(zf[0]))), q_low_g_kg=float(qq[lo]),
                q_mid_g_kg=float(qq[mid]), vapour_drop_g_kg=float(qq[lo] - qq[mid]),
                lid_rise_k=float(th[lid] - th[mid]), levels_km=[round(float(zf[i]) / 1000.0, 2) for i in (lo, mid, lid)])


def _km(x):
    return None if x is None else x / 1000.0


def cm1_quarter(cm1: dict, quarter: str, z_full, z_half) -> dict:
    """CM1's means over the quarter at the column's levels: the measures, the profiles and the fluxes."""
    m = cm1['quarters'][PBL_QUARTERS[quarter]]
    zf, zhh = np.asarray(z_full, float), np.asarray(z_half, float)
    thv = np.interp(zf, cm1['zh_m'], m['theta_v_k'])
    qq = np.interp(zf, cm1['zh_m'], m['q_g_kg'])
    out = measures(zf, thv, qq)
    out['mixed_depth_km'] = _km(_cm1_depth(cm1, m, float(zf[-1])))
    out.update(theta_v_k=thv, q_g_kg=qq, heat_w_m2=np.interp(zhh, cm1['zw_m'], m['heat_w_m2']),
               vapour_w_m2=np.interp(zhh, cm1['zw_m'], m['vapour_w_m2']),
               vapour_resolved_w_m2=np.interp(zhh, cm1['zw_m'], m['vapour_resolved_w_m2']),
               heat_resolved_w_m2=np.interp(zhh, cm1['zw_m'], m['heat_resolved_w_m2']),
               hpbl_km=float(m['hpbl_m']) / 1000.0)
    return out


def _cm1_depth(cm1, m, z0):
    """CM1's mixed layer by the same rule on its own levels, from the column's lowest full level."""
    from climate.gcm import pbl_mixing as pm
    return pm.rise_height(cm1['zh_m'], m['theta_v_k'], z0)


def gap_closed(cm1_value, control, new):
    """The share of the control's gap to CM1 that the new scheme closes (1 closes it, 0 none, negative widens it)."""
    if cm1_value is None or control is None or new is None:
        return None
    gap = abs(control - cm1_value)
    return None if gap == 0.0 else 1.0 - abs(new - cm1_value) / gap


def analyse() -> dict:
    from climate.gcm import pbl_mixing as pm
    half = grid()
    boxes, cm1s, rains = {}, {}, {}
    results = {}
    for name, site, case, span, start_q, compare_q in CASES:
        if case not in boxes:
            boxes[case] = gridded_rows(box_rows(case, span), half)
            cm1s[case] = pm.cm1_box(case)
            rains[case] = mf.cm1_microphysics(case, span)
        rows = boxes[case]
        runs = {s[0]: run_column(rows, start_q, compare_q, site, s, half) for s in SCHEMES}
        ref = runs['earth_local']['quarters'][compare_q]
        cm1 = cm1_quarter(cm1s[case], compare_q, ref['z_full_m'], ref['z_half_m'])
        entry = dict(site=site, case=case, start_quarter=start_q, compared_quarter=compare_q,
                     cm1=dict(measures={k: cm1[k] for k in ('mixed_depth_km', 'q_low_g_kg', 'q_mid_g_kg',
                                                            'vapour_drop_g_kg', 'lid_rise_k', 'levels_km', 'hpbl_km')},
                              rain_ground_mm_day={q: rains[case][q]['ground_mm_day'] for q in (start_q, compare_q)},
                              theta_v_k=cm1['theta_v_k'], q_g_kg=cm1['q_g_kg'], heat_w_m2=cm1['heat_w_m2'],
                              vapour_w_m2=cm1['vapour_w_m2'], vapour_resolved_w_m2=cm1['vapour_resolved_w_m2'],
                              heat_resolved_w_m2=cm1['heat_resolved_w_m2']),
                     schemes={})
        for sname, run in runs.items():
            qd = run['quarters'][compare_q]
            sq = run['quarters'][start_q]
            entry['schemes'][sname] = dict(
                mixing_length_m=run['mixing_length_m'], nonlocal_=run['nonlocal_'],
                measures=measures(qd['z_full_m'], qd['theta_v_k'], qd['q_g_kg']),
                **{k: v for k, v in qd.items() if k not in ('z_full_m', 'z_half_m', 't_k')},
                z_full_m=qd['z_full_m'], z_half_m=qd['z_half_m'],
                spin_up_quarter=dict(rain_conv_mm_day=sq['rain_conv_mm_day'], rain_ls_mm_day=sq['rain_ls_mm_day'],
                                     convective_share=sq['convective_share'], types=sq['types'],
                                     top_median_km=sq['top_median_km'], dca_share=sq['dca_share'],
                                     water_mm_day=sq['water_mm_day'], heat_w_m2=sq['heat_w_m2'],
                                     **({'h_median_km': sq['h_median_km']} if 'h_median_km' in sq else {})))
        cm = entry['cm1']['measures']
        ctrl = entry['schemes']['earth_local']['measures']
        entry['gap_closed'] = {s: {k: gap_closed(cm[k], ctrl[k], entry['schemes'][s]['measures'][k])
                                   for k in ('mixed_depth_km', 'vapour_drop_g_kg', 'lid_rise_k')}
                               for s in ('lunar_local', 'earth_nonlocal', 'lunar_nonlocal')}
        # a sensitivity: the same columns started from the box's air at the compared quarter's own start
        short = {s[0]: run_column(rows, compare_q, compare_q, site, s, half)['quarters'][compare_q] for s in SCHEMES}
        short_m = {s: measures(q['z_full_m'], q['theta_v_k'], q['q_g_kg']) for s, q in short.items()}
        entry['one_quarter'] = dict(
            measures=short_m,
            gap_closed={s: {k: gap_closed(cm[k], short_m['earth_local'][k], short_m[s][k])
                            for k in ('mixed_depth_km', 'vapour_drop_g_kg', 'lid_rise_k')}
                        for s in ('lunar_local', 'earth_nonlocal', 'lunar_nonlocal')},
            rain_mm_day={s: q['rain_conv_mm_day'] + q['rain_ls_mm_day'] for s, q in short.items()},
            deep_share={s: q['types']['deep'] for s, q in short.items()})
        results[name] = entry
    night = results['sea night']
    stop = {s: dict(depth=night['gap_closed'][s]['mixed_depth_km'], drop=night['gap_closed'][s]['vapour_drop_g_kg'])
            for s in ('lunar_local', 'earth_nonlocal', 'lunar_nonlocal')}
    passes = {s: bool(v['depth'] is not None and v['drop'] is not None and v['depth'] >= 0.5 and v['drop'] >= 0.5)
              for s, v in stop.items()}
    lid = {s: night['schemes'][s]['measures']['lid_rise_k'] >= 0.5 * night['cm1']['measures']['lid_rise_k']
           for s in night['schemes']}
    return dict(schema=SCHEMA, level_set=LEVEL_SET, step_s=STEP_S, cases=results,
                stop_rule=dict(rule=('the new scheme closes at least half of the control column\'s gap to CM1 in the '
                                     'sea night\'s mixed depth and in its 0.3-6.2 km vapour drop'),
                               closed=stop, passes=passes, new_scheme='lunar_nonlocal',
                               passed=passes['lunar_nonlocal']),
                lid=dict(rule='the sea night\'s virtual potential temperature rise from 6.2 to 8.2 km stays at least half '
                              'of CM1\'s', kept=lid),
                parameters=ps.PARAMETERS, evidence=EVIDENCE, reading_rule=READING_RULE, producer=producer())


EVIDENCE = ('Single-column calculations on PlaSim\'s fifteen levels driven by the CM1 boxes\' own forcing (one '
            'realization each of the settled sea box and the 0 E land box on a 6 km grid); PlaSim\'s vertical diffusion '
            'and dry adjustment as Python copies of fluxmod.f90 and rainmod.f90, the nonlocal scheme as written in '
            'pbl_scheme.py, the mass-flux convection scheme of massflux.py. No GCM or CM1 run is made.')
READING_RULE = ('cases: each column\'s compared quarter (its second), by scheme: earth_local is PlaSim as it stands '
                '(160 m), lunar_local the length stretched by the gravities (966 m), earth_nonlocal and lunar_nonlocal '
                'the nonlocal scheme beneath each. measures: mixed_depth_km by pbl_mixing.rise_height from the lowest '
                'full level, q_low/q_mid at the full levels nearest 0.3 and 6.2 km (mixing ratio, g/kg), '
                'vapour_drop_g_kg their difference, lid_rise_k the virtual potential temperature at the level nearest '
                '8.2 km less that near 6.2 km, levels_km the levels used. Profiles top first at z_full_m; fluxes at the '
                'interfaces z_half_m (top first), upward, W/m2: heat_bl and vapour_bl by vdiff (diffusion with the '
                'nonlocal fluxes), heat_dca and vapour_dca by the dry adjustment; CM1\'s total (resolved and '
                'parameterized) at the same heights. dca_share: steps on which the dry adjustment acted. Convection: '
                'types and convective_share by step, tops (km) of convecting steps, conv_moistening_mm_day and '
                'conv_heating_w_m2 by layer. gap_closed: 1 - |new - CM1| / |earth_local - CM1|. one_quarter: a '
                'sensitivity, each column started from the box\'s air at its compared quarter\'s own start.')


def producer() -> dict:
    files = ('gcm/pbl_column.py', 'gcm/pbl_scheme.py', 'gcm/massflux.py', 'gcm/pbl_mixing.py', 'crm/box_budget.py')
    return dict(domain='climate', files={f: hashlib.sha256((HERE.parent / f).read_bytes()).hexdigest()[:16]
                                         for f in files})


def compact(x, digits: int = 4):
    if isinstance(x, dict):
        return {str(k): compact(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [compact(v, digits) for v in x]
    if isinstance(x, np.ndarray):
        return compact(x.tolist(), digits)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return None if not math.isfinite(float(x)) else float(f'{float(x):.{digits}g}')
    return x


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--output', type=Path, default=OUT)
    a = p.parse_args(argv)
    res = analyse()
    a.output.write_text(json.dumps(compact(res), indent=1) + '\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
