"""The simplified Betts-Miller scheme offline: the rain it would make on the CM1 site boxes' own columns and on the
GCM's columns at the same sites, by quarter of the lunar day, at each setting of its reference humidity, its parcel's
mixing and its inhibition limit, beside the rain the boxes and the GCM made.

    climate/gcm/.venv/bin/python -m climate.gcm.sbm_offline      # 15 seconds: ../results/gcm/sbm_offline.json

The scheme is sbm.column, the tested Python copy of the Fortran routine, given each column as PlaSim would hold it and
the constants of the runs (the gas constant and gravity from the GCM run's planet namelist, ExoPlaSim's kappa of
0.286, PlaSim's saturation and latent heats). Its rain is the rate at once: the vapour its relaxation removes per
second, the rain formed before PlaSim's evaporation of falling rain (mkrain), in mm/day.

- The boxes' columns are CM1's domain means at each 3-hourly output over the second lunar day (cm1out_diag:
  temperature, specific humidity, pressure, surface pressure). The main mapping, 'level', interpolates temperature
  and vapour in log pressure to PlaSim's ten full levels at the box's surface pressure. Two others show what the
  mapping does: 'layer' takes each of PlaSim's layers' mean by mass, and 'cm1' runs the scheme on CM1's own 111
  levels. CM1's rain is its domain-mean rate at the same outputs, as its water budget reckons it.
- The GCM's columns are its 3-day means (ta, hus, ps) in each of its cells at the site (gcm_site.site_cells: six sea
  cells, four land), at the hour angle of the site from the means' subsolar longitude. Beside them stand the GCM's
  convective rain at the ground (prc, a mean) and the rain its scheme formed at the end of each mean (dt_conv, a
  snapshot of one step, 18 degrees of hour angle later), its column integral of cp times the heating over L.

A box's columns were made by CM1's convection, so the scheme meets them away from its own equilibrium: its rain at
once is the vapour above its reference over the relaxation time, while a model running the scheme keeps its column
near that reference and its rain follows what supplies the vapour. The comparison on CM1's air therefore shows when
the scheme would convect there and how far that air stands from the reference; on the GCM's own columns, held near
the reference, it shows what each setting does at once to the GCM's day and night.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np

from climate.gcm import sbm

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
OUT = RESULTS / 'sbm_offline.json'
SCHEMA = 'terluna.climate.gcm-sbm-offline/1'
DAY = 86400.0
KAPPA = 0.286                      # ExoPlaSim's p_exo.f90; the runs' planet namelist leaves it
RV, ACPV = 461.51, 1870.0          # plasimmod.f90
STEP_S = 3600.0                    # PlaSim's deltsec2 at the runs' 30-minute step
TAU_S = 7200.0
EPS = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
RH = (0.6, 0.7, 0.8)
CIN = ((0.7, 0.0, 10.0), (0.7, 0.0, 30.0), (0.7, 1.0, 10.0), (0.7, 1.0, 30.0))     # (rh, eps, cin limit J/kg)
MAPPINGS = ('level', 'layer', 'cm1')
QUARTERS = (('sunrise-noon', -90.0, 0.0), ('noon-sunset', 0.0, 90.0), ('sunset-midnight', 90.0, 180.0),
            ('midnight-sunrise', -180.0, -90.0))
BOXES = (('sea', 'box_near_side_sea'), ('land', 'box_0e_design'))
GCM_BRANCHES = (('A28R_301k_cold_sbm_rh07', 0.7), ('A28R_301k_cold_sbm_rh06', 0.6))   # each at its own rh
GCM_YEARS = (10, 14)
SPAN_DAYS = (29.5, 59.0)
# The pass test: the scheme's rain on CM1's columns near CM1's at the sea after midnight and in the morning, and in the
# afternoon and evening at both sites; near means within half a millimetre a day or half of CM1's rain.
PASS_QUARTERS = (('sea', 'midnight-sunrise'), ('sea', 'sunrise-noon'), ('sea', 'noon-sunset'),
                 ('sea', 'sunset-midnight'), ('land', 'noon-sunset'), ('land', 'sunset-midnight'))
PASS_MM_DAY, PASS_SHARE = 0.5, 0.5


def run_constants(run: str = GCM_BRANCHES[0][0]) -> dict:
    """PlaSim's constants as a run took them: its gas constant and gravity from its planet namelist, kappa, and the
    derived specific heat, gas-constant ratio and vapour's heat-capacity term, over sbm.PLASIM's saturation and
    latent heats."""
    text = (HERE / 'runs' / run / 'model' / 'planet_namelist').read_text()
    value = lambda key: float(re.search(rf'^\s*{key}\s*=\s*([-+0-9.eEdD]+)', text, re.M | re.I).group(1))
    gascon, ga = value('GASCON'), value('GA')
    acpd = gascon / KAPPA
    return dict(sbm.PLASIM, gascon=gascon, ga=ga, acpd=acpd, rdbrv=gascon / RV, ADV=ACPV / acpd - 1.0)


def quarter_of(hour_deg) -> np.ndarray:
    """The index in QUARTERS of each hour angle (degrees, 0 = noon, -180..180)."""
    h = (np.asarray(hour_deg, float) + 180.0) % 360.0 - 180.0
    out = np.full(h.shape, -1)
    for k, (_, a, b) in enumerate(QUARTERS):
        out[(h >= a) & (h < b)] = k
    return out


def by_quarter(values, quarter, weights=None) -> list:
    v, w = np.asarray(values, float), np.ones(len(values)) if weights is None else np.asarray(weights, float)
    return [float(np.sum(v[quarter == k] * w[quarter == k]) / np.sum(w[quarter == k])) for k in range(len(QUARTERS))]


def w_level_pressures(p, ps, zh, zw):
    """Pressure (Pa) at CM1's w levels zw from its scalar levels' pressures p (the ground up): log-linear between
    scalar levels, the surface pressure at the ground and the top level's lapse of log pressure carried to the top."""
    lnp = np.log(np.asarray(p, float))
    pw = np.exp(np.interp(zw, zh, lnp))
    pw[0] = ps
    pw[-1] = np.exp(lnp[-1] + (lnp[-1] - lnp[-2]) / (zh[-1] - zh[-2]) * (zw[-1] - zh[-1]))
    return pw


def layer_means(values, pw, edges_pa):
    """Means by mass of a profile held as layer values between the pressures pw (the ground up, falling) over the
    target layers between edges_pa (top first, rising), each over the part of it the profile covers."""
    v = np.asarray(values, float)
    top_first = np.asarray(pw, float)[::-1]
    held = np.concatenate([[0.0], np.cumsum((v * (np.asarray(pw)[:-1] - np.asarray(pw)[1:]))[::-1])])
    f = np.interp(edges_pa, top_first, held)
    covered = np.asarray(edges_pa[1:], float) - np.maximum(np.asarray(edges_pa[:-1], float), top_first[0])
    return (f[1:] - f[:-1]) / covered


def on_levels(col: dict, mapping: str, sigma, edges, zh=None, zw=None):
    """A CM1 column as the scheme takes it: temperature, vapour, sigma and layer thickness in sigma, top first."""
    ps = col['ps']
    if mapping == 'level':
        lnp = np.log(col['p'][::-1])
        x = np.log(np.asarray(sigma) * ps)
        return (np.interp(x, lnp, col['t'][::-1]), np.interp(x, lnp, col['q'][::-1]), np.asarray(sigma),
                np.diff(edges))
    pw = w_level_pressures(col['p'], ps, zh, zw)
    if mapping == 'layer':
        target = np.asarray(edges) * ps
        return layer_means(col['t'], pw, target), layer_means(col['q'], pw, target), np.asarray(sigma), np.diff(edges)
    if mapping == 'cm1':
        return col['t'][::-1], col['q'][::-1], col['p'][::-1] / ps, ((pw[:-1] - pw[1:]) / ps)[::-1]
    raise ValueError(f'mapping is one of {MAPPINGS}, not {mapping!r}')


def cm1_columns(case: str, span=SPAN_DAYS) -> dict:
    """A box's domain-mean columns at each diagnostic output over the span: temperature (K), specific humidity,
    pressure (Pa) by level (the ground up), surface pressure, and the box's rain (mm/day) and CAPE and inhibition
    (J/kg, CM1's), with each output's hour angle at the site."""
    from climate.crm import box_budget as bb
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    folder = ra.RUNS / case
    record = json.loads((folder / 'case.json').read_text())
    start, lon = record['configuration']['start_hour_angle_deg'], record['configuration']['site']['lon_deg']
    geo = ra.case_geometry(folder)
    layout = ra.grads_layout(folder / 'cm1out_diag_s.ctl')[1]
    cols = []
    for path in sorted(folder.glob('cm1out_diag_*_s.dat')):
        d = bb.read_diag(folder, path, layout)
        if not span[0] * DAY - 1.0 <= d['mtime'] <= span[1] * DAY + 1.0:
            continue
        cols.append(dict(day=d['mtime'] / DAY, hour=local_hour_angle(start, lon, d['mtime'], solar_day_s()),
                         t=d['t'], q=d['shv'], p=d['prs'], ps=d['psfc'], rain=d['prate'] * DAY, cape=d['cape'],
                         cin=d['cin']))
    return dict(case=case, columns=cols, zh=geo['zh'], zw=geo['zw'],
                inputs={'runner': record.get('runner'), 'gcm': record.get('gcm'), 'build': record['build']['label']})


def gcm_columns(run: str, years, box: str) -> dict:
    """The GCM's 3-day mean columns in each of its cells at a box's site (the cells the box's case records), with
    their weights (cos latitude), the hour angle of the site, the GCM's convective rain at the ground (prc, mm/day)
    and the rain its scheme formed at the end of the mean (mm/day, from the dt_conv snapshot). The humidity's means
    come from year_files, which leaves out each year file's first and last records."""
    import netCDF4
    from climate.crm import gcm_site as gs
    from climate.crm import ring_analysis as ra
    from climate.gcm import year_files as yf
    record = json.loads((ra.RUNS / box / 'case.json').read_text())
    cells = [tuple(c) for c in record['gcm_site']['cells']]
    lon_site = record['configuration']['site']['lon_deg']
    c = run_constants(run)
    out = dict(t=[], q=[], ps=[], weight=[], hour=[], year=[], prc=[], pr=[], formed=[])
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(HERE / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            edges = gs.plasim_edges(sigma)
            idx = [(int(np.argmin(abs(lat - a))), int(np.argmin(abs(lon - b)))) for a, b in cells]
            w = np.cos(np.radians([lat[j] for j, _ in idx]))
            f = yf.fields(d, ('ta', 'hus', 'ps', 'prc', 'pr', 'dt_conv', 'czen'))     # the humidity's means: year_files
            hour = (lon_site - gs.subsolar_longitude(f['czen'], lat, lon) + 180.0) % 360.0 - 180.0
        for n in range(f['ta'].shape[0]):
            for wk, (j, i) in zip(w, idx):
                ps = f['ps'][n, j, i] * 100.0
                q = f['hus'][n, :, j, i]
                m = ps * np.diff(edges) / c['ga']
                out['t'].append(f['ta'][n, :, j, i])
                out['q'].append(q)
                out['ps'].append(ps)
                out['weight'].append(wk)
                out['hour'].append(hour[n])
                out['year'].append(year)
                out['prc'].append(f['prc'][n, j, i] * 1000.0 * DAY)
                out['pr'].append(f['pr'][n, j, i] * 1000.0 * DAY)
                out['formed'].append(float((c['acpd'] * (1.0 + c['ADV'] * q) * f['dt_conv'][n, :, j, i] * m).sum())
                                     / c['ALV'] * DAY)
    return dict({k: np.asarray(v) for k, v in out.items()}, sigma=sigma, edges=edges, cells=cells)


def scheme(t, q, sigma, dsigma, ps, rh, eps, cin_limit, c, akap=KAPPA, cape_limit=0.0, eps_trigger=0.0) -> dict:
    """The scheme on one column: its rain at once (mm/day), whether it convects deep and rains, whether its parcel is
    buoyant anywhere, and the parcel's inhibition (J/kg); with the CAPE trigger set, held back below cape_limit."""
    o = sbm.column(t, q, sigma, dsigma, ps, tau=TAU_S, rh=rh, eps=eps, cin_limit=cin_limit, step=STEP_S, akap=akap, c=c,
                   cape_limit=cape_limit, eps_trigger=eps_trigger)
    return dict(rain=float(o['rain'].sum()) * DAY, deep=o['top'] is not None and not o['shallow'],
                buoyant=o['top'] is not None, cin=float(o['cin']))


def settings():
    """Every (rh, eps, cin limit) the test takes."""
    return [(rh, eps, 0.0) for rh in RH for eps in EPS] + list(CIN)


def summarise(results, quarter, weights=None) -> dict:
    rain = [r['rain'] for r in results]
    deep = [float(r['deep']) for r in results]
    buoyant = [float(r['buoyant']) for r in results]
    cin = np.array([r['cin'] if r['buoyant'] else np.nan for r in results])
    w = np.ones(len(results)) if weights is None else np.asarray(weights, float)
    cin_q = []
    for k in range(len(QUARTERS)):
        m = (quarter == k) & np.isfinite(cin)
        cin_q.append(float(np.sum(cin[m] * w[m]) / np.sum(w[m])) if m.any() else None)
    rain_q = by_quarter(rain, quarter, weights)
    return dict(rain_mm_day=rain_q, vapour_above_reference_kg_m2=[r * TAU_S / DAY for r in rain_q],
                convecting_share=by_quarter(deep, quarter, weights),
                buoyant_share=by_quarter(buoyant, quarter, weights),
                parcel_inhibition_j_kg=cin_q)


def passes(scheme_rain: dict, cm1_rain: dict) -> dict:
    """The pass test on one setting: each of PASS_QUARTERS within PASS_MM_DAY or PASS_SHARE of CM1's rain."""
    names = [q[0] for q in QUARTERS]
    rows = []
    for site, quarter in PASS_QUARTERS:
        s, b = scheme_rain[site][names.index(quarter)], cm1_rain[site][names.index(quarter)]
        rows.append(dict(site=site, quarter=quarter, scheme=s, cm1=b,
                         near=bool(abs(s - b) <= max(PASS_MM_DAY, PASS_SHARE * b))))
    return dict(passes=all(r['near'] for r in rows), quarters=rows)


def on_gcm_branch(run: str, rh: float, years, c: dict) -> dict:
    """The scheme at the branch's own reference humidity and each mixing on the branch's own cell columns at both
    sites, beside the branch's rain at the ground and its scheme's rain formed at the end of each mean."""
    out = dict(run=run, rh=rh, years=list(years), sites={}, settings=[])
    gcm = {site: gcm_columns(run, years, case) for site, case in BOXES}
    for site, g in gcm.items():
        quarter = quarter_of(g['hour'])
        snap = quarter_of(np.asarray(g['hour']) + 18.0)
        out['sites'][site] = dict(cells=[list(x) for x in g['cells']], cell_outputs=int(len(g['hour'])),
                                  prc_mm_day=by_quarter(g['prc'], quarter, g['weight']),
                                  formed_snapshot_mm_day=by_quarter(g['formed'], snap, g['weight']),
                                  formed_snapshot_share=by_quarter((g['formed'] > 0.05).astype(float), snap,
                                                                   g['weight']))
    for eps in EPS:
        row = dict(rh=rh, eps=eps, cin_limit_j_kg=0.0)
        for site, g in gcm.items():
            quarter = quarter_of(g['hour'])
            res = [scheme(t, q, g['sigma'], np.diff(g['edges']), ps, rh, eps, 0.0, c)
                   for t, q, ps in zip(g['t'], g['q'], g['ps'])]
            row[site] = summarise(res, quarter, g['weight'])
            rain = np.array([r['rain'] for r in res])
            per_year = [by_quarter(rain[g['year'] == y], quarter[g['year'] == y], g['weight'][g['year'] == y])
                        for y in range(years[0], years[1] + 1)]
            row[site]['rain_range_over_years_mm_day'] = [[float(np.min(v)), float(np.max(v))]
                                                        for v in np.array(per_year).T]
            row[site]['convecting_cell_outputs'] = [int(sum(r['deep'] for r, k in zip(res, quarter) if k == j))
                                                    for j in range(len(QUARTERS))]
        out['settings'].append(row)
    return out


# The CAPE trigger on the rh 0.7 branch's own columns (24 October): at most five settings (the CAPE in J/kg a parcel
# mixing eps_trigger per e-folding of pressure must have), beside the scheme without it and with the reference itself
# mixed at eps 1, as Stage A's screen took it. CM1's rain by quarter, the sea box's third lunar day and the land box's
# second, is the target; a setting passes if the sea's rain after midnight falls to at most half the branch's and
# within 0.5 mm/day of CM1's, while the afternoon and evening at both sites stay within a factor of two of CM1's. Each
# setting is judged on the offline rain itself and on the branch's own rain scaled quarter by quarter by the
# setting's offline rain over the scheme's without it.
TRIGGERS = ((50.0, 0.0), (100.0, 0.0), (200.0, 0.0), (5.0, 1.0), (5.0, 4.0))
TRIGGER_BRANCH = ('A28R_301k_cold_sbm_rh07', 0.7)
CM1_TARGET = dict(sea=[0.00, 2.23, 2.97, 0.14], land=[2.91, 8.59, 1.53, 0.08])


def trigger_passes(rain: dict, branch_night: float) -> dict:
    """The trigger test's rule on one set of rain by quarter at both sites."""
    names = [q[0] for q in QUARTERS]
    night = rain['sea'][names.index('midnight-sunrise')]
    rows = [dict(site='sea', quarter='midnight-sunrise', rain=night, cm1=CM1_TARGET['sea'][3],
                 near=bool(night <= 0.5 * branch_night and night <= CM1_TARGET['sea'][3] + 0.5))]
    for site in ('sea', 'land'):
        for quarter in ('noon-sunset', 'sunset-midnight'):
            r, b = rain[site][names.index(quarter)], CM1_TARGET[site][names.index(quarter)]
            rows.append(dict(site=site, quarter=quarter, rain=r, cm1=b, near=bool(0.5 * b <= r <= 2.0 * b)))
    return dict(passes=all(r['near'] for r in rows), quarters=rows)


def trigger_test(years=GCM_YEARS) -> dict:
    run, rh = TRIGGER_BRANCH
    c = run_constants(run)
    gcm = {site: gcm_columns(run, years, case) for site, case in BOXES}
    settings = [dict(label='without the trigger', eps=0.0, cape_limit=0.0, eps_trigger=0.0),
                dict(label='the reference mixed at eps 1', eps=1.0, cape_limit=0.0, eps_trigger=0.0)]
    settings += [dict(label=f'CAPE above {cape:g} J/kg, trigger mixing {mix:g}', eps=0.0, cape_limit=cape,
                      eps_trigger=mix) for cape, mix in TRIGGERS]
    branch = {site: by_quarter(g['pr'], quarter_of(g['hour']), g['weight']) for site, g in gcm.items()}
    held = {}                       # the CAPE the trigger parcel finds on the columns where the scheme rains without it
    for site, g in gcm.items():
        quarter = quarter_of(g['hour'])
        raining = np.array([scheme(t, q, g['sigma'], np.diff(g['edges']), ps, rh, 0.0, 0.0, c)['rain'] > 0.0
                            for t, q, ps in zip(g['t'], g['q'], g['ps'])])
        held[site] = {}
        for mix in (0.0, 1.0):
            cape = np.array([sbm.trigger_cape(t, sbm.parcel(t, q, g['sigma'], ps * g['sigma'], mix, KAPPA, c),
                                              g['sigma'], np.diff(g['edges']), c)
                             for t, q, ps in zip(g['t'], g['q'], g['ps'])])
            held[site][f'trigger mixing {mix:g}'] = [float(np.median(cape[raining & (quarter == k)]))
                                                     if (raining & (quarter == k)).any() else None
                                                     for k in range(len(QUARTERS))]
    rows = []
    for st in settings:
        row = dict(st)
        for site, g in gcm.items():
            quarter = quarter_of(g['hour'])
            res = [scheme(t, q, g['sigma'], np.diff(g['edges']), ps, rh, st['eps'], 0.0, c,
                          cape_limit=st['cape_limit'], eps_trigger=st['eps_trigger'])
                   for t, q, ps in zip(g['t'], g['q'], g['ps'])]
            row[site] = summarise(res, quarter, g['weight'])
        rows.append(row)
    base = rows[0]
    for row in rows:
        row['branch_scaled_mm_day'] = {
            site: [b * (r / s0) if s0 > 0.0 else 0.0
                   for b, r, s0 in zip(branch[site], row[site]['rain_mm_day'], base[site]['rain_mm_day'])]
            for site in ('sea', 'land')}
        row['test'] = dict(offline=trigger_passes({s: row[s]['rain_mm_day'] for s in ('sea', 'land')},
                                                  branch['sea'][3]),
                           branch_scaled=trigger_passes(row['branch_scaled_mm_day'], branch['sea'][3]))
    return dict(run=run, rh=rh, years=list(years), branch_rain_mm_day=branch, cm1_target_mm_day=CM1_TARGET,
                raining_columns_cape_median_j_kg=held, settings=rows,
                passing=[r['label'] for r in rows[2:] if r['test']['offline']['passes'] or
                         r['test']['branch_scaled']['passes']])


def analyse(span=SPAN_DAYS, branches=GCM_BRANCHES, years=GCM_YEARS) -> dict:
    c = run_constants(branches[0][0])
    boxes = {site: cm1_columns(case, span) for site, case in BOXES}
    import netCDF4
    from climate.crm import gcm_site as gs
    with netCDF4.Dataset(HERE / 'runs' / branches[0][0] / 'model' / f'MOST.{years[0]:05d}.nc') as d:
        sigma = np.asarray(d['lev'][:], float)
    edges = gs.plasim_edges(sigma)
    cm1 = {}
    for site, b in boxes.items():
        quarter = quarter_of([x['hour'] for x in b['columns']])
        cm1[site] = dict(case=b['case'], outputs=len(b['columns']),
                         span_days=[b['columns'][0]['day'], b['columns'][-1]['day']],
                         outputs_by_quarter=[int((quarter == k).sum()) for k in range(len(QUARTERS))],
                         rain_mm_day=by_quarter([x['rain'] for x in b['columns']], quarter),
                         cape_j_kg=by_quarter([x['cape'] for x in b['columns']], quarter),
                         inhibition_j_kg=by_quarter([x['cin'] for x in b['columns']], quarter), inputs=b['inputs'])
    on_boxes = {}
    for mapping in MAPPINGS:
        rows = []
        for rh, eps, cin_limit in settings():
            if mapping != 'level' and (rh != 0.7 or cin_limit):
                continue
            row = dict(rh=rh, eps=eps, cin_limit_j_kg=cin_limit)
            for site, b in boxes.items():
                quarter = quarter_of([x['hour'] for x in b['columns']])
                res = [scheme(*on_levels(x, mapping, sigma, edges, b['zh'], b['zw']), x['ps'], rh, eps, cin_limit, c)
                       for x in b['columns']]
                row[site] = summarise(res, quarter)
            if mapping == 'level':
                row['test'] = passes({s: row[s]['rain_mm_day'] for s in ('sea', 'land')},
                                     {s: cm1[s]['rain_mm_day'] for s in ('sea', 'land')})
            rows.append(row)
        on_boxes[mapping] = rows
    tests = dict(sbm.PLASIM, ga=c['ga'])                          # test_sbm.py's constants, for their effect
    on_boxes['level_test_constants'] = []
    for eps in (0.0, 1.0):
        row = dict(rh=0.7, eps=eps, cin_limit_j_kg=0.0)
        for site, b in boxes.items():
            quarter = quarter_of([x['hour'] for x in b['columns']])
            row[site] = summarise([scheme(*on_levels(x, 'level', sigma, edges), x['ps'], 0.7, eps, 0.0, tests,
                                          akap=tests['gascon'] / tests['acpd']) for x in b['columns']], quarter)
        on_boxes['level_test_constants'].append(row)
    on_gcm = [on_gcm_branch(run, rh, years, run_constants(run)) for run, rh in branches]
    return dict(schema=SCHEMA, quarters=[q[0] for q in QUARTERS], constants=dict(
                    gascon=c['gascon'], ga=c['ga'], kappa=KAPPA, acpd=c['acpd'], rdbrv=c['rdbrv'], ADV=c['ADV'],
                    tau_s=TAU_S, step_s=STEP_S),
                sigma=sigma.tolist(), cm1=cm1, on_boxes=on_boxes, on_gcm=on_gcm, trigger_test=trigger_test(years),
                pass_rule=dict(quarters=[list(q) for q in PASS_QUARTERS], within_mm_day=PASS_MM_DAY,
                               within_share_of_cm1=PASS_SHARE, mapping='level'),
                passing_settings=[dict(rh=r['rh'], eps=r['eps'], cin_limit_j_kg=r['cin_limit_j_kg'])
                                  for r in on_boxes['level'] if r['test']['passes']],
                evidence=('Offline diagnosis: the tested Python copy of the scheme applied at once to columns, with no '
                          'model integrating it. The boxes\' columns are CM1 domain means over one lunar day of one '
                          'run each, standing for a GCM cell, forced by Kuo\'s GCM; the GCM\'s are 3-day means of five '
                          'years of each of two branches, in which the scheme ran without mixing at rh 0.7 and 0.6. '
                          'The scheme\'s '
                          'instantaneous rate is set against CM1\'s realised rain and the GCM\'s rain.'),
                reading_rule=('Rain in mm/day by quarter of the lunar day at the site (sunrise-noon, noon-sunset, '
                              'sunset-midnight, midnight-sunrise); on the scheme\'s side the rain formed, before '
                              'PlaSim\'s evaporation of falling rain; vapour_above_reference_kg_m2 is that rain '
                              'times the '
                              'relaxation time, the vapour above the reference in the convecting layers, nil where the '
                              'scheme does not rain. convecting_share is the share of columns on which '
                              'the scheme convects deep and rains, buoyant_share the share on which its parcel is '
                              'buoyant somewhere, parcel_inhibition_j_kg the parcel\'s mean inhibition where buoyant. '
                              'on_boxes holds each setting on CM1\'s columns by mapping (level: interpolation to '
                              'PlaSim\'s full levels, the pass test\'s; layer: layer means by mass; cm1: CM1\'s own '
                              'levels; level_test_constants: the level mapping with test_sbm.py\'s constants in '
                              'place of the runs\'); on_gcm, for each branch, the scheme at the branch\'s own rh on '
                              'its own cell columns, '
                              'cell-weighted by cos latitude, '
                              'with the range of the five yearly values and the count of convecting cell-outputs. '
                              'formed_snapshot is the GCM\'s own scheme at the end of each mean, one step. '
                              'trigger_test: the CAPE trigger (sbm.column\'s cape_limit and eps_trigger) on the rh 0.7 '
                              'branch\'s cell columns beside the scheme without it and with its reference mixed at eps '
                              '1; branch_scaled multiplies the branch\'s own rain (pr) in each quarter by the '
                              'setting\'s '
                              'offline rain over the scheme\'s without the trigger; test applies the rule of '
                              'TRIGGERS\' comment to both; raining_columns_cape_median_j_kg is the median CAPE the '
                              'trigger parcel finds, undiluted or mixing at 1, on the columns where the scheme without '
                              'the trigger rains.'))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args(argv)
    out = analyse()
    out['producer'] = dict(domain='climate', files={f'gcm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'sbm.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    fmt = lambda v: ' '.join(f'{x:7.2f}' for x in v)
    for site in ('sea', 'land'):
        print(f"{site}: CM1 {fmt(out['cm1'][site]['rain_mm_day'])}   "
              f"inhibition {fmt(out['cm1'][site]['inhibition_j_kg'])}")
        for r in out['on_boxes']['level']:
            print(f"  rh {r['rh']:.1f} eps {r['eps']:4.1f} cin {r['cin_limit_j_kg']:4.0f}: "
                  f"{fmt(r[site]['rain_mm_day'])}"
                  f"   convecting {fmt(r[site]['convecting_share'])}")
        for b in out['on_gcm']:
            g = b['sites'][site]
            print(f"  {b['run']} prc {fmt(g['prc_mm_day'])}   formed (snapshots) {fmt(g['formed_snapshot_mm_day'])}")
            for r in b['settings']:
                print(f"    its columns rh {r['rh']:.1f} eps {r['eps']:4.1f}: {fmt(r[site]['rain_mm_day'])}"
                      f"   convecting {fmt(r[site]['convecting_share'])}   buoyant {fmt(r[site]['buoyant_share'])}")
    print('passing settings:', out['passing_settings'] or 'none')
    t = out['trigger_test']
    print('the CAPE trigger on', t['run'], 'years', t['years'], '; the branch\'s rain', t['branch_rain_mm_day'])
    for r in t['settings']:
        print(f"  {r['label']:42s} sea {fmt(r['sea']['rain_mm_day'])}  land {fmt(r['land']['rain_mm_day'])}  "
              f"branch scaled: sea {fmt(r['branch_scaled_mm_day']['sea'])}  "
              f"land {fmt(r['branch_scaled_mm_day']['land'])}"
              f"  passes {r['test']['offline']['passes']}/{r['test']['branch_scaled']['passes']}")
    print('passing triggers:', t['passing'] or 'none')
    print(f'wrote {OUT.relative_to(HERE.parents[1])}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
