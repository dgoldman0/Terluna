"""An entrainment that follows the cloud's environment for the mass-flux scheme, tested against the CM1 boxes (27
October 2026; the record is massflux.md's section "Entrainment that follows the environment" and README.md's section of
the same date).

    climate/gcm/.venv/bin/python -m climate.gcm.massflux_entrainment   # ../results/gcm/massflux_entrainment_l15.json

The scheme's fixed profile mixes 0.8e-4 /m aloft and 2e-4 /m more decaying over 2 km above the cloud base, shallow
clouds twice that, whatever surrounds the cloud. The CM1 sea box's cloudy updraughts mix 1.2-1.6e-3 /m at night and
0.2-0.3e-3 /m in the afternoon. The form tested (massflux.PARAMETERS entrainment_form 'environment') replaces the
profile's part near the base with the IFS's entrainment (Cy40r1 (6.7); Bechtold and others 2008, 2014; de Rooy and
others 2013, (26)): eps = eps_bg + eps_env (1.3 - RH) (qs / qs_base)^3, RH the surrounding air's relative humidity and
qs, qs_base its saturation specific humidity at the level and at the cloud base, shallow clouds twice that, with the
IFS's eps_env = 1.75e-3 /m and the CM1 boxes' eps_bg = 0.8e-4 /m aloft. Three parts:

1. On CM1's own environments (cm1_check): each box's 3-D output over its settled lunar day, sampled as box_entrainment
   samples it, by quarter. The form is evaluated on the quarter's mean sounding at the levels where CM1's cloudy
   updraughts are sampled, with the cloud base at their lowest sampled level, and weighted by their mass flux as CM1's
   bulk rate is (below 16 km, where CM1's rate is diagnosed); beside it the caps of a plume mixing at the form's rates
   in the mean sounding (box_entrainment.plume_cap), against the columns' cloud tops.
2. The column test (columns): pbl_column's testbed on fifteen levels, PlaSim's local mixing at Earth's length, the
   mass-flux scheme with the fixed profile and with the form, each column started a quarter early: the sea's afternoon
   (from sunrise), evening (from noon) and night (from sunset) and the land's afternoon (from sunrise). The criteria
   (CRITERIA) were set before the run.
3. The running column from the GCM's air (gcm_air): massflux.scm on fifteen levels at the sea, started from the rh 0.7
   branch's site mean and run three lunar days, each form: whether the form still carries the GCM's warm upper air
   toward CM1's, as the fixed profile does (massflux.md, test 4).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from climate.gcm import massflux as mf
from climate.gcm import pbl_column as pc

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
OUT = RESULTS / 'massflux_entrainment_l15.json'
SCHEMA = 'terluna.climate.gcm-massflux-entrainment/1'
ENVIRONMENT = dict(mf.PARAMETERS, entrainment_form='environment')
FORMS = (('fixed', mf.PARAMETERS), ('environment', ENVIRONMENT))
BOXES = (('sea', 'box_near_side_sea_day3', (59.0, 88.5)), ('land', 'box_0e_design', (29.5, 59.0)))
QUARTERS = ('night', 'morning', 'afternoon', 'evening')
DEEP_PA = 20000.0                       # the scheme's deep_dp: a cloud thicker than this is deep
UPPER_SIGMA = (0.4368, 0.31685)          # the fifteen levels' full levels near 41.8 and 55.8 km at the sea site
# The columns (name, site, CM1 case, span, the quarter started in, the quarter compared)
CASES = (('sea afternoon', 'sea', 'box_near_side_sea_day3', (59.0, 88.5), 'morning', 'afternoon'),
         ('sea evening', 'sea', 'box_near_side_sea_day3', (59.0, 88.5), 'afternoon', 'evening'),
         ('sea night', 'sea', 'box_near_side_sea_day3', (59.0, 88.5), 'evening', 'night'),
         ('land afternoon', 'land', 'box_0e_design', (29.5, 59.0), 'morning', 'afternoon'))
SCHEME = ('earth_local', 'earth', False)  # PlaSim's own local mixing at Earth's length
# The pass criteria, set before the run (the brief of 27 October). The sea night's gaps are the testbed's own with the
# fixed profile (pbl_column_l15.json, earth_local: 4.87 km and 3.34 g/kg) against CM1's 6.48 km and 0.34 g/kg; CM1's
# rain at the ground over the boxes' settled days (2.27, 2.96 and 8.71 mm/day); the lid is CM1's 1.38 K.
CRITERIA = dict(
    night_tops_below_10km_share=0.75,        # of the night's steps, convecting or not, none topped above 10 km
    night_mixed_depth_km=0.5 * (4.87 + 6.48),
    night_vapour_drop_g_kg=0.5 * (3.34 + 0.34),
    night_lid_rise_k=0.5 * 1.38,
    night_rain_mm_day=0.5,
    sea_afternoon_rain_mm_day=(0.5 * 2.27, 1.5 * 2.27),
    sea_evening_rain_mm_day=2.0 / 3.0 * 2.96,
    land_afternoon_rain_mm_day=(0.5 * 8.71, 1.5 * 8.71),
    afternoon_tops_above_20km_share=0.5,     # of the convecting steps, in each afternoon column
    upper_air_tolerance_k=0.1)                # |form - CM1| at most |fixed - CM1| + this, each day column and height


# ------------------------------------------------------------------------------------------- 1. CM1's environments

def form_profile(t, q, p, k_base: int, kind: str, prm=ENVIRONMENT) -> np.ndarray:
    """The form's entrainment (per metre) at each level of a sounding (levels the ground up; temperature K, specific
    humidity, pressure Pa) for a cloud based at level k_base: massflux.entrainment with the surrounding air's relative
    humidity and its saturation at the level and at the base, PlaSim's saturation as the scheme takes it."""
    qs = np.array([mf.qsat(a, b) for a, b in zip(t, p)])
    rh = np.asarray(q, float) / qs
    return np.array([mf.entrainment(kind, 0.0, rh[k], 0.0, prm, qs[k], qs[k_base]) for k in range(len(t))])


def cm1_check() -> dict:
    """By box and quarter: CM1's bulk entrainment of its cloudy updraughts (qt and h), the form's bulk on the quarter's
    mean sounding at the same levels and weights, deep and shallow, which of the two the quarter's clouds are by the
    scheme's 200 hPa rule (base to the columns' median top), and the caps of plumes mixing at the form's rates and at
    CM1's bulk rate in the mean sounding."""
    from climate.crm import box_entrainment as be
    out = {}
    for site, case, span in BOXES:
        src = be.read_source(case, span)
        info, zh, zw = src['info'], np.asarray(src['info']['zh'], float), src['info']['zw']
        g = info['g']
        out[site] = dict(case=case, span_days=list(span), quarters={})
        for quarter in QUARTERS:
            group = src['groups'][quarter]
            if not group['rows']:
                continue
            pooled = be.pooled(group['rows'])
            e = be.entrainment(pooled, zh)['updraft']
            snd = be.sounding(pooled, zh, zw)
            q = snd['qv'] / (1.0 + snd['qv'])                       # CM1's mixing ratio as specific humidity
            measured = np.isfinite(e['eps_qt']) & (zh < be.LOWER_M) & (e['mass_flux'] > 0.0)
            tops = np.concatenate(group['tops']) if group['tops'] else np.zeros(0)
            entry = dict(cm1_bulk_qt_per_m=e['bulk_qt'], cm1_bulk_h_per_m=e['bulk_h'],
                         cloudy_share=float(tops.size / max(group['columns'], 1)),
                         top_median_m=float(np.median(tops)) if tops.size else None,
                         top_p90_m=float(np.percentile(tops, 90)) if tops.size else None)
            if measured.any():
                k_base = int(np.flatnonzero(measured)[0])
                p_top = float(np.interp(entry['top_median_m'], zh, snd['p'])) if tops.size else None
                depth = float(snd['p'][k_base] - p_top) if p_top is not None else None
                levels = []
                bulk = {}
                for kind in ('deep', 'shallow'):
                    prof = form_profile(snd['t'], q, snd['p'], k_base, kind)
                    bulk[kind] = be.bulk(zh, np.where(measured, prof, np.nan), e['mass_flux'])
                    if kind == 'deep':
                        qs = np.array([mf.qsat(a, b) for a, b in zip(snd['t'], snd['p'])])
                        for k in np.flatnonzero(measured):
                            levels.append(dict(z_m=float(zh[k]), mass_flux_kg_m2_s=float(e['mass_flux'][k]),
                                               cm1_eps_qt_per_m=float(e['eps_qt'][k]),
                                               cm1_eps_h_per_m=(None if not np.isfinite(e['eps_h'][k])
                                                                else float(e['eps_h'][k])),
                                               rh=float(q[k] / qs[k]), qs_ratio_cubed=float((qs[k] / qs[k_base]) ** 3),
                                               form_deep_per_m=float(prof[k])))
                kind = 'deep' if depth is not None and depth > DEEP_PA else 'shallow'
                entry.update(base_m=float(zh[k_base]), depth_to_median_top_pa=depth, kind=kind,
                             form_bulk_deep_per_m=bulk['deep'], form_bulk_shallow_per_m=bulk['shallow'],
                             form_bulk_per_m=bulk[kind], levels=levels)
                # plumes from the lowest 0.5 km mixed: at the form's rates (its base the plume's own condensation
                # level) and at CM1's bulk rate, condensate kept (reversible) and dropped
                undiluted = be.plume_cap(snd, zh, zw, 0.0, g)
                plumes = {}
                if undiluted['lcl_m'] is not None:
                    kb = int(np.argmin(np.abs(zh - undiluted['lcl_m'])))
                    for label, rate in (('form, deep', form_profile(snd['t'], q, snd['p'], kb, 'deep')),
                                        ('form, shallow', form_profile(snd['t'], q, snd['p'], kb, 'shallow')),
                                        ('CM1 bulk by qt', e['bulk_qt']), ('fixed profile, deep',
                                         fixed_profile(zh, kb, 'deep')), ('fixed profile, shallow',
                                                                         fixed_profile(zh, kb, 'shallow'))):
                        plumes[label] = {f: be.plume_cap(snd, zh, zw, rate, g, fallout)['cap_m']
                                         for f, fallout in (('reversible', False), ('pseudo_adiabatic', True))}
                entry['plume_caps_m'] = plumes
            out[site]['quarters'][quarter] = entry
    sea = out['sea']['quarters']
    out['contrast'] = dict(
        cm1_night_over_afternoon_qt=sea['night']['cm1_bulk_qt_per_m'] / sea['afternoon']['cm1_bulk_qt_per_m'],
        cm1_night_over_afternoon_h=sea['night']['cm1_bulk_h_per_m'] / sea['afternoon']['cm1_bulk_h_per_m'],
        form_night_over_afternoon=sea['night']['form_bulk_per_m'] / sea['afternoon']['form_bulk_per_m'],
        form_night_kind=sea['night']['kind'], form_afternoon_kind=sea['afternoon']['kind'],
        target=dict(night_per_m=(1.2e-3, 1.6e-3), afternoon_per_m=(0.2e-3, 0.3e-3)))
    return out


def fixed_profile(zh, k_base: int, kind: str, prm=mf.PARAMETERS) -> np.ndarray:
    """The fixed profile's entrainment (per metre) at heights zh for a cloud based at level k_base."""
    zb = float(zh[k_base])
    return np.array([mf.entrainment(kind, float(z), 0.8, zb, prm) for z in zh])


# ------------------------------------------------------------------------------------------------- 2. the columns

def column_entry(run: dict, compare_q: str, rows: list, cm1_m: dict) -> dict:
    """A column's compared quarter: the comparison's measures, rain at the ground, convection's types and tops by
    step, and the upper air at the levels near 41.8 and 55.8 km against CM1's mean over the same quarter."""
    qd = run['quarters'][compare_q]
    steps = qd['steps_record']
    tops = [s['top_m'] for s in steps if s['type'] is not None and s['top_m'] is not None]
    flux = [s for s in steps if s['type'] is not None and (s['mass_flux'] or 0.0) > 0.0]
    mb = np.array([s['mass_flux'] for s in flux]) if flux else np.zeros(0)
    n = len(steps)
    sigma = np.asarray(rows[0]['p_full'], float) / float(rows[0]['p_half'][-1])
    k_up = [int(np.argmin(np.abs(sigma - s))) for s in UPPER_SIGMA]
    cm1_t = np.mean([r['t'] for r in rows if r['quarter'] == compare_q], axis=0)
    zf = np.asarray(qd['z_full_m'], float)
    return dict(
        measures=pc.measures(qd['z_full_m'], qd['theta_v_k'], qd['q_g_kg']),
        rain_ground_mm_day=qd['rain_conv_mm_day'] + qd['rain_ls_mm_day'],
        rain_convective_mm_day=qd['rain_conv_mm_day'], rain_large_scale_mm_day=qd['rain_ls_mm_day'],
        steps=n, convecting_share=qd['convective_share'], types=qd['types'],
        tops_below_10km_share=sum(1 for s in steps if s['type'] is None or s['top_m'] is None
                                  or s['top_m'] < 10000.0) / n,
        convecting_tops_above_20km_share=(sum(1 for z in tops if z > 20000.0) / len(tops)) if tops else None,
        top_median_km=float(np.median(tops)) / 1000.0 if tops else None,
        top_p90_km=float(np.percentile(tops, 90)) / 1000.0 if tops else None,
        deep_with_mass_flux_share=sum(1 for s in flux if s['type'] == 'deep') / n,
        mass_flux_median_kg_m2_s=float(np.median(mb)) if mb.size else None,
        base_median_km=float(np.median([s['base_m'] for s in flux])) / 1000.0 if flux else None,
        flux_weighted_top_km=float(np.sum(mb * [s['top_m'] for s in flux]) / mb.sum()) / 1000.0 if flux else None,
        upper_air=[dict(sigma=float(sigma[k]), z_km=float(zf[k]) / 1000.0, t_k=float(qd['t_k'][k]),
                        cm1_t_k=float(cm1_t[k]), minus_cm1_k=float(qd['t_k'][k] - cm1_t[k])) for k in k_up],
        dca_share=qd['dca_share'], water_mm_day=qd['water_mm_day'], heat_w_m2=qd['heat_w_m2'],
        conv_heating_w_m2=qd['conv_heating_w_m2'], conv_moistening_mm_day=qd['conv_moistening_mm_day'],
        theta_v_k=qd['theta_v_k'], q_g_kg=qd['q_g_kg'], z_full_m=qd['z_full_m'], cm1=cm1_m)


def columns() -> dict:
    from climate.gcm import pbl_mixing as pm
    half = pc.grid()
    rows_by_case, cm1s, rains = {}, {}, {}
    out = {}
    for name, site, case, span, start_q, compare_q in CASES:
        if case not in rows_by_case:
            rows_by_case[case] = pc.gridded_rows(pc.box_rows(case, span), half)
            cm1s[case] = pm.cm1_box(case)
            rains[case] = mf.cm1_microphysics(case, span)
        rows = rows_by_case[case]
        forms = {}
        for form, prm in FORMS:
            run = pc.run_column(rows, start_q, compare_q, site, SCHEME, half, prm=prm, keep_steps=True)
            ref = run['quarters'][compare_q]
            cm1 = pc.cm1_quarter(cm1s[case], compare_q, ref['z_full_m'], ref['z_half_m'])
            cm1_m = {k: cm1[k] for k in ('mixed_depth_km', 'q_low_g_kg', 'q_mid_g_kg', 'vapour_drop_g_kg',
                                         'lid_rise_k', 'levels_km', 'hpbl_km')}
            entry = column_entry(run, compare_q, rows, cm1_m)
            sq = run['quarters'][start_q]
            entry['spin_up_quarter'] = dict(rain_ground_mm_day=sq['rain_conv_mm_day'] + sq['rain_ls_mm_day'],
                                            types=sq['types'], top_median_km=sq['top_median_km'])
            forms[form] = entry
        out[name] = dict(site=site, case=case, start_quarter=start_q, compared_quarter=compare_q,
                         cm1_rain_ground_mm_day={q: rains[case][q]['ground_mm_day'] for q in (start_q, compare_q)},
                         forms=forms)
    return out


def judge(cols: dict) -> dict:
    """The criteria (CRITERIA) on the form's columns, each with its value, its threshold and whether it is met."""
    c = CRITERIA
    night = cols['sea night']['forms']['environment']
    m = night['measures']
    res = dict(
        night_tops_below_10km=dict(value=night['tops_below_10km_share'], at_least=c['night_tops_below_10km_share']),
        night_mixed_depth=dict(value=m['mixed_depth_km'], at_least=c['night_mixed_depth_km']),
        night_vapour_drop=dict(value=m['vapour_drop_g_kg'], at_most=c['night_vapour_drop_g_kg']),
        night_lid_rise=dict(value=m['lid_rise_k'], at_least=c['night_lid_rise_k']),
        night_rain=dict(value=night['rain_ground_mm_day'], at_most=c['night_rain_mm_day']),
        sea_afternoon_rain=dict(value=cols['sea afternoon']['forms']['environment']['rain_ground_mm_day'],
                                within=c['sea_afternoon_rain_mm_day']),
        sea_evening_rain=dict(value=cols['sea evening']['forms']['environment']['rain_ground_mm_day'],
                              at_least=c['sea_evening_rain_mm_day']),
        land_afternoon_rain=dict(value=cols['land afternoon']['forms']['environment']['rain_ground_mm_day'],
                                 within=c['land_afternoon_rain_mm_day']))
    for name in ('sea afternoon', 'land afternoon'):
        v = cols[name]['forms']['environment']['convecting_tops_above_20km_share']
        res[f'{name.replace(" ", "_")}_tops_above_20km'] = dict(value=v, at_least=c['afternoon_tops_above_20km_share'])
    for name in ('sea afternoon', 'sea evening', 'land afternoon'):
        f = cols[name]['forms']
        for i, label in enumerate(('41.8 km', '55.8 km')):
            new = abs(f['environment']['upper_air'][i]['minus_cm1_k'])
            old = abs(f['fixed']['upper_air'][i]['minus_cm1_k'])
            res[f'{name.replace(" ", "_")}_upper_air_{label.replace(" ", "")}'] = dict(
                value=new, at_most=old + c['upper_air_tolerance_k'], fixed=old)
    for r in res.values():
        v = r['value']
        if v is None:
            r['met'] = False
        elif 'within' in r:
            r['met'] = bool(r['within'][0] <= v <= r['within'][1])
        elif 'at_least' in r:
            r['met'] = bool(v >= r['at_least'])
        else:
            r['met'] = bool(v <= r['at_most'])
    night_keys = [k for k in res if k.startswith('night_')]
    return dict(criteria=res, night_passes=all(res[k]['met'] for k in night_keys),
                day_passes=all(r['met'] for k, r in res.items() if k not in night_keys),
                passes=all(r['met'] for r in res.values()))


# -------------------------------------------------------------------------------- 3. the running column, GCM's air

def gcm_air(cycles: int = mf.SCM_CYCLES) -> dict:
    """massflux.scm on fifteen levels at the sea from the rh 0.7 branch's site mean, each form: the column's temperature
    less CM1's mean at 41.8 and 55.8 km (the nearest full levels) at the start and in each lunar day, and its rain by
    quarter in the last."""
    from climate.gcm import levels as lv
    grids = mf.grid_sets()
    box = dict(outputs=mf.cm1_outputs('box_near_side_sea_day3', (59.0, 88.5)))
    g = mf.gcm_site_columns(mf.GCM_RUN, mf.GCM_YEARS, 'box_near_side_sea')
    w = g['w'] / g['w'].sum()
    t10, q10 = np.tensordot(w, g['t'], 1), np.tensordot(w, g['q'], 1)
    start = dict(t=lv.remap(t10, grids['10'], grids['15']), q=lv.remap(q10, grids['10'], grids['15']))
    out = dict(run=mf.GCM_RUN, years=list(mf.GCM_YEARS), lunar_days=cycles, forms={})
    for form, prm in FORMS:
        s = mf.scm(box, 'sea', grids['15'], prm=prm, start=start, cycles=cycles)
        idx = [int(np.argmin(abs(s['z_m'] - z))) for z in (41800.0, 55800.0)]
        out['forms'][form] = dict(
            height_m=[float(s['z_m'][k]) for k in idx],
            start_minus_cm1_k=[float(s['t_start'][k] - s['t_cm1'][k]) for k in idx],
            mean_minus_cm1_k_by_lunar_day=[[float(h['t_mean'][k] - s['t_cm1'][k]) for k in idx] for h in s['cycles']],
            rain_mm_day_last_day=s['rain_mm_day'])
    return out


# ------------------------------------------------------------------------------------------------------ the record

LITERATURE = (
    dict(source='IFS Documentation Cy40r1 (2014), Part IV, sections 6.3-6.4', state='read',
         where='ecmwf.int elibrary 9204 (sha256 f9fec6f0246be4f3), text in session-tools/agents/climate-massflux/'
               'literature: (6.7) the entrainment, 1.75e-3 /m (1.3 - RH) (qsat(T)/qsat(T_base))^3, fe 1 deep, 2 '
               'shallow and mid-level, applied above the cloud base to positively buoyant convection; (6.8)-(6.9) the '
               'detrainment; the test parcels\' 0.4 times it (6.21)'),
    dict(source='Bechtold and others (2008), QJRMS 134, 1337-1351', state='read as ECMWF Technical Memorandum 556',
         where='elibrary 8011 (sha256 431b13baaf05e12e), section 2.1: eps = c0 Fe + c1 (qs - q)/qs Fe, '
               'Fe = (qs/qs_b)^2, c0 and c1 of order 1e-4 /m, the organised part for deep buoyant convection'),
    dict(source='Bechtold and others (2014), J. Atmos. Sci. 71, 734-753',
         state='read as ECMWF Technical Memorandum 705',
         where='elibrary 8018 (sha256 2e057dd36bfae7bd), appendix A: 1.8e-3 /m (1.3 - RH) fs, shallow twice deep, '
               '"as also supported by CRM data"'),
    dict(source='de Rooy and others (2013), QJRMS 139, 1-19', state='read: sections 2, 4.3-4.5 and 5',
         where='KNMI\'s copy (cdn.knmi.nl, sha256 d3b6d4fc8f609e90): (26), 1.8e-3 /m, applied where the updraught is '
               'buoyant, twice for a test cloud under 200 hPa; LES rates fall as clouds deepen (Fig. 8); Gregory\'s '
               '(2001) beta B/w^2 (beta 0.03 deep, 0.06 shallow) and beta/z give too little sensitivity in the ECMWF '
               'SCM; detrainment varies more than entrainment'),
    dict(source='Derbyshire and others (2004), QJRMS 130, 3055-3079', state='not read: Wiley refused the download',
         where='taken as de Rooy and others (2013, section 4.4) and Bechtold and others (2008) describe it'),
    dict(source='Gregory (2001), QJRMS 127; Neggers and others (2002), J. Atmos. Sci. 59; Simpson (1971), J. Atmos. '
                'Sci. 28 (0.2/R)', state='not read', where='taken as de Rooy and others (2013) give them'))
EVIDENCE = ('Offline single-column calculations, no new CM1 or GCM run: CM1\'s 3-D output of one realization each of '
            'the settled sea box and the 0 E land box on a 6 km grid, sampled as box_entrainment.py samples it; the '
            'columns in pbl_column.py\'s testbed driven by the boxes\' own forcing; the running column of massflux.py '
            'from the rh 0.7 branch\'s site mean. The form\'s check on CM1\'s environments evaluates it on mean '
            'soundings at the levels where CM1 samples its clouds, and CM1\'s bulk rates come from the bulk-plume '
            'method, which Romps (2010) found gives about half the directly measured mixing.')
READING_RULE = ('cm1_check: by box and quarter, CM1\'s bulk entrainment (per metre, mass-flux-weighted below 16 km, by '
                'qt and h) and the form\'s at the same levels and weights, deep and shallow, the quarter\'s kind by '
                'the scheme\'s 200 hPa rule from the sampled base to the columns\' median top; plume_caps_m: the cap '
                'of a plume from the lowest 0.5 km mixed through the quarter\'s mean sounding at each rate, reversible and '
                'pseudo-adiabatic. columns: each column\'s compared quarter with the fixed profile and the form: '
                'measures as pbl_column\'s, rain at the ground (convective after PlaSim\'s evaporation of falling rain '
                'at the box\'s constant, and large-scale), tops_below_10km_share of all steps (a step without '
                'convection counts as below), convecting_tops_above_20km_share of the convecting steps, '
                'deep_with_mass_flux_share of all steps, the cloud-base mass flux, base and mass-flux-weighted top of '
                'the steps that convect with a mass flux, upper_air at '
                'the full levels near 41.8 and 55.8 km against CM1\'s mean on the same levels over the quarter. '
                'judgement: CRITERIA, set before the run. gcm_air: the running column\'s temperature less CM1\'s mean '
                'at 41.8 and 55.8 km.')


def producer() -> dict:
    return dict(domain='climate', files={
        name: hashlib.sha256((HERE.parent / name).read_bytes()).hexdigest()[:16]
        for name in ('gcm/massflux_entrainment.py', 'gcm/massflux.py', 'gcm/pbl_column.py', 'gcm/pbl_scheme.py',
                     'crm/box_entrainment.py')})


def analyse(parts=('cm1', 'columns', 'gcm_air')) -> dict:
    out = dict(schema=SCHEMA, form=dict(entrainment_form='environment', eps_env_per_m=ENVIRONMENT['eps_env'],
                                        eps_bg_per_m=ENVIRONMENT['eps_bg'],
                                        shallow_factor=ENVIRONMENT['shallow_factor'],
                                        rule='eps_bg + eps_env (1.3 - RH) (qs / qs_base)^3, shallow clouds times '
                                             'shallow_factor; test parcels\' sub-cloud mixing unchanged'),
               criteria=CRITERIA, literature=LITERATURE, evidence=EVIDENCE, reading_rule=READING_RULE,
               producer=producer())
    if 'cm1' in parts:
        out['cm1_check'] = cm1_check()
    if 'columns' in parts:
        out['columns'] = columns()
        out['judgement'] = judge(out['columns'])
    if 'gcm_air' in parts:
        out['gcm_air'] = gcm_air()
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--parts', nargs='+', choices=('cm1', 'columns', 'gcm_air'), default=['cm1', 'columns', 'gcm_air'])
    p.add_argument('--output', type=Path, default=OUT)
    a = p.parse_args(argv)
    res = pc.compact(analyse(tuple(a.parts)))
    a.output.write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res.get('judgement', {}), indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
