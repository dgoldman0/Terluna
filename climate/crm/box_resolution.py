"""A sea box's held hour on three grids: whether the Near-side Sea box stays quiet because its 6 km columns cannot
resolve its boundary layer's eddies. Two tests use it, each read by rules set before its runs (TESTS):

- night (25 October 2026): box_near_side_sea_night_6km, _3km and _2km start from box_near_side_sea_day3's air at day
  84.5, 40 degrees past midnight, when that Kuo-forced box stayed quiet under a 6 km deep boundary layer;
- evening (28 October 2026): box_near_side_sea_l15_evening_6km, _3km and _2km start from box_near_side_sea_l15_day4's
  air at day 107.625, 142 degrees of hour angle, in the late evening's lull of the box forced by the fifteen-level GCM,
  when that GCM rains about 2.8 mm/day at the site.

    climate/gcm/.venv/bin/python -m climate.crm.box_resolution                  # the night test
    climate/gcm/.venv/bin/python -m climate.crm.box_resolution --test evening   # the evening test
                                       # ../results/crm/box_resolution_near_side_sea_night.json, _l15_evening.json

Each test's three boxes start from the source's air at one output averaged over its columns and hold that hour's
forcing for two Earth days; they differ in their grid alone. For each, and for the source box itself over the same two
days, this reads:

- every domain diagnostic (half-hourly in the boxes, 3-hourly in the source): the rain, from CM1's running domain total;
  CM1's own CAPE and inhibition, the means of its columns' 0.5 km mixed-layer parcels; the boundary layer's depth, the
  mean of its columns' bulk-Richardson heights; the strongest updraft, the domain maximum of w over the levels below
  16 km and below 60 km; the highest cloud, the highest level at which any column holds cloud water or ice of at least
  0.001 g/kg (CM1's own threshold there); and box_onset's parcel from the lowest 0.5 km of the domain-mean sounding,
  undiluted and mixing at eps 8;
- every 3-D output (6-hourly in the boxes, 3-hourly in the source): box_onset's snapshot measures, among them the
  99.9th percentiles of w below 16 km and at 16-60 km and the shares of columns by cloud top (0.01 g/kg), and where
  the cloud stands: the tops of the cloud below 40 km, the columns with cloud at 16-40 km, and the cloud above 40 km
  that stands apart from any below it, as a thin ice layer at 50-52 km does in all these runs at night;
- the domain-mean drift from the start, at fixed heights and over the GCM's lowest four layers, and the water budget
  over the second day and over both (box_budget.water).

The control's check and the reading of the three grids follow rules set before the runs (CONTROL and READING).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from climate.crm import box_budget as bb
from climate.crm import box_onset as bo
from climate.crm import cm1_output as co
from climate.crm import ring_analysis as ra

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-box-resolution/1'
DAY, HOUR = 86400.0, 3600.0
GRIDS = {'6 km': 'box_near_side_sea_night_6km', '3 km': 'box_near_side_sea_night_3km',
         '2 km': 'box_near_side_sea_night_2km'}
SOURCE = dict(case='box_near_side_sea_day3', first_day=84.5, last_day=86.5)
DAYS = (('first day', 0.0, 24.0), ('second day', 24.0, 48.0))         # hours from each run's start
SHALLOW_M, DEEP_M, UPPER_M = 16000.0, 40000.0, 60000.0
HEIGHTS_M = (1000.0, 5000.0, 16000.0, 40000.0, 60000.0, 100000.0)   # where the drift is tabulated
GCM_NIGHT_RAIN_MM_DAY = 1.45     # the GCM with its simplified Betts-Miller scheme at the site after midnight
# Set before the runs. The control reproduces the source box's quiet night if over its second day it rains under
# 0.3 mm/day (twice the source's mean after midnight, 0.14), its strongest updraft below 60 km averages under 3 m/s
# (the source's night: 1.5 m/s, its afternoon storms 10.1), the 99.9th percentile of its w below 16 km stays under 2 m/s
# in every snapshot (the source's night: 1.0-1.25 m/s) and under 5% of its columns, on average over its snapshots, hold
# cloud tops above 40 km (the source's night: 1%, its evening 63%).
CONTROL = dict(rain_mm_day=0.3, strongest_updraft_m_s=3.0, w_p999_below_16km_m_s=2.0, deep_top_share=0.05)
# A grid convects if over its second day it rains at least a third of the GCM's night rain (0.5 mm/day) with deep
# cells: its strongest updraft averaging at least 5 m/s or cloud tops above 40 km in at least 5% of its columns.
READING = dict(rain_mm_day=0.5, strongest_updraft_m_s=5.0, deep_top_share=0.05)
# The evening test (28 October 2026), set before its runs. Its source over the two days after the held hour rains
# 0.09-0.11 mm/day, its strongest updraft below 60 km averaging 1.0 m/s and the 99.9th percentile of its w below 16 km
# at most 0.9 m/s, while its inhibition falls from 10 to 3 J/kg; in the two days before, the sunset's storms rained
# 0.35-0.76 mm/day with updrafts averaging 1.8-2.3 m/s. The control reproduces the lull if over its second day it
# rains under 0.4 mm/day (twice the source's at 135-150 degrees, 0.2), its strongest updraft averages under 3 m/s and
# the 99.9th percentile of its w below 16 km stays under 2 m/s in every snapshot. The source carries the sunset storms'
# ice above 40 km over much of the domain, which a box started from mean air does not, so no cloud-top share enters.
EVENING_CONTROL = dict(rain_mm_day=0.4, strongest_updraft_m_s=3.0, w_p999_below_16km_m_s=2.0)
# A grid convects if over its second day it rains at least a third of the GCM's rain at that hour (2.8 mm/day, its
# 3-day means at 135-150 degrees) with deep cells: its strongest updraft averaging at least 5 m/s, or the 99.9th
# percentile of its w at 16-60 km averaging at least 1.5 m/s (the source's storms after midnight: 2.4-2.7 m/s).
EVENING_READING = dict(rain_mm_day=0.93, strongest_updraft_m_s=5.0, w_p999_16_60km_m_s=1.5)

NIGHT_EVIDENCE = ("Three CM1 boxes, one realization each, started from one domain-mean sounding of the sea box's third"
                  " lunar day with random perturbations and run two Earth days under that hour's forcing held fixed; "
                  'they differ in grid spacing (6, 3 and 2 km) and the 2 km box in its width (128 km against 192). The'
                  ' source box over the same two days, under its forcing by local time, stands beside them. A '
                  "numerical test of the grid's effect within CM1 at these settings, with its boundary-layer scheme on"
                  ' at every spacing.')
NIGHT_READING_RULE = ("Each grid's days are the 24 hours from its start and the next 24; the source's are days "
                      "84.5-85.5 and 85.5-86.5 of its run. Rain in mm/day from CM1's running domain total. Means over "
                      "the diagnostics after each day's start to its end (half-hourly in the boxes, 3-hourly in the "
                      "source): CM1's CAPE and inhibition (J/kg, its columns' 0.5 km mixed parcels), the boundary "
                      "layer's depth (hpbl, m), the strongest updraft (the domain maximum of w below 16 and below 60 "
                      'km, m/s; a domain of 4,096 columns samples more than one of 1,024), the highest cloud (m, any '
                      'column holding 0.001 g/kg of cloud water or ice) and the 0.5 km parcel of the domain-mean '
                      'sounding (box_onset, undiluted and at eps 8). snapshot: means over the 3-D outputs in the day '
                      "(6-hourly in the boxes, 3-hourly in the source) of box_onset's measures: the 99.9th percentiles"
                      ' of w, the highest cloud top and the shares of columns by cloud top (0.01 g/kg); clouds: the '
                      'mean and highest top of the cloud below 40 km, the share of columns with cloud at 16-40 km and '
                      "the share whose cloud above 40 km has none at 16-40 km beneath it. quiet: the control's "
                      'criteria, set before the runs (control_rule), on the second day; convects: rain of at least '
                      "reading_rule_limits' 0.5 mm/day with the strongest updraft averaging 5 m/s or cloud tops above "
                      '40 km in 5% of the columns. drift: the domain-mean temperature at fixed heights (K) and the '
                      "mass-weighted temperature (C) and vapour (g/kg) of the GCM's lowest four layers at the start, "
                      'after a day and at the end; water: box_budget.water over the second day and both, rest closes '
                      "the vapour and condensate budget (mm/day), and refill_of_empty_air the vapour CM1's "
                      'microphysics puts back in levels holding no water, where the prescribed table dries air above '
                      'about 55 km that has none to give, which the rest holds.')
EVENING_EVIDENCE = ('Three CM1 boxes, one realization each, started from one domain-mean sounding of the sea box '
                    "forced by the fifteen-level GCM, in its fourth lunar day's late evening (day 107.625, 142 degrees"
                    " of hour angle), with random perturbations, and run two Earth days under that hour's forcing held"
                    ' fixed; they differ in grid spacing (6, 3 and 2 km) and the 2 km box in its width (128 km against'
                    ' 192). The source box over the same two days, under its forcing by local time, stands beside '
                    "them. A numerical test of the grid's effect within CM1 at these settings, with its boundary-layer"
                    " scheme on at every spacing; the GCM's rain at that hour is its 3-day means at the site binned by"
                    ' hour angle.')
EVENING_READING_RULE = ("Each grid's days are the 24 hours from its start and the next 24; the source's are days "
                        "107.625-108.625 and 108.625-109.625 of its run. Rain in mm/day from CM1's running domain "
                        "total. Means over the diagnostics after each day's start to its end (half-hourly in the "
                        "boxes, 3-hourly in the source): CM1's CAPE and inhibition (J/kg, its columns' 0.5 km mixed "
                        "parcels), the boundary layer's depth (hpbl, m), the strongest updraft (the domain maximum of "
                        'w below 16 and below 60 km, m/s; a domain of 4,096 columns samples more than one of 1,024), '
                        'the highest cloud (m, any column holding 0.001 g/kg of cloud water or ice) and the 0.5 km '
                        'parcel of the domain-mean sounding (box_onset, undiluted and at eps 8). snapshot: means over '
                        "the 3-D outputs in the day (6-hourly in the boxes, 3-hourly in the source) of box_onset's "
                        'measures: the 99.9th percentiles of w, the highest cloud top and the shares of columns by '
                        'cloud top (0.01 g/kg); clouds: the mean and highest top of the cloud below 40 km, the share '
                        'of columns with cloud at 16-40 km and the share whose cloud above 40 km has none at 16-40 km '
                        "beneath it. quiet: the control's criteria, set before the runs (control_rule), on the second "
                        "day; convects: rain of at least reading_rule_limits' 0.93 mm/day, a third of the GCM's 2.8 at"
                        ' that hour, with the strongest updraft averaging 5 m/s or the 99.9th percentile of w at 16-60'
                        ' km averaging 1.5 m/s. drift: the domain-mean temperature at fixed heights (K) and the '
                        "mass-weighted temperature (C) and vapour (g/kg) of the GCM's lowest four layers at the start,"
                        ' after a day and at the end; water: box_budget.water over the second day and both, rest '
                        "closes the vapour and condensate budget (mm/day), and refill_of_empty_air the vapour CM1's "
                        'microphysics puts back in levels holding no water, where the prescribed table dries air above'
                        ' about 55 km that has none to give, which the rest holds.')
TESTS = {
    'night': dict(grids=GRIDS, source=SOURCE, control=CONTROL, reading=READING, quiet_spell='quiet sea nights',
                  gcm_key='gcm_night_rain_mm_day', gcm_rain_mm_day=GCM_NIGHT_RAIN_MM_DAY,
                  product='box_resolution_near_side_sea_night.json', evidence=NIGHT_EVIDENCE,
                  reading_rule=NIGHT_READING_RULE),
    'evening': dict(grids={'6 km': 'box_near_side_sea_l15_evening_6km', '3 km': 'box_near_side_sea_l15_evening_3km',
                           '2 km': 'box_near_side_sea_l15_evening_2km'},
                    source=dict(case='box_near_side_sea_l15_day4', first_day=107.625, last_day=109.625),
                    control=EVENING_CONTROL, reading=EVENING_READING, quiet_spell='quiet late evenings',
                    gcm_key='gcm_rain_mm_day', gcm_rain_mm_day=2.8,
                    product='box_resolution_near_side_sea_l15_evening.json', evidence=EVENING_EVIDENCE,
                    reading_rule=EVENING_READING_RULE),
}


def cloud_top_m(qcfrac, qifrac, zh):
    """The highest level at which any column holds cloud water or ice (CM1's fractions of columns by level), or None."""
    cloudy = (np.asarray(qcfrac) > 0.0) | (np.asarray(qifrac) > 0.0)
    return float(np.asarray(zh)[cloudy].max()) if cloudy.any() else None


def rain_over(rows, a_h, b_h):
    """Rain (mm/day) between hours a_h and b_h from CM1's running domain total (cm), linear between outputs."""
    hours = np.array([r['hour'] for r in rows])
    total = np.array([r['rain_total_mm'] for r in rows])
    if a_h < hours[0] - 1e-6 or b_h > hours[-1] + 1e-6:
        return None
    return float((np.interp(b_h, hours, total) - np.interp(a_h, hours, total)) / ((b_h - a_h) / 24.0))


def diag_rows(name: str, first_s: float = 0.0, last_s: float = np.inf) -> list:
    """Every domain diagnostic of a case between model times first_s and last_s, as hours from first_s."""
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    g = record['build']['gravity_m_s2']
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    dz = np.diff(zw)
    s_layout = ra.grads_layout(case / 'cm1out_diag_s.ctl')[1]
    w_layout = ra.grads_layout(case / 'cm1out_diag_w.ctl')[1]
    rows = []
    for path in sorted(case.glob('cm1out_diag_*_s.dat')):
        d = bb.read_diag(case, path, s_layout)
        if not first_s - 1.0 <= d['mtime'] <= last_s + 1.0:
            continue
        w = bb.read_diag(case, path.with_name(path.name.replace('_s.dat', '_w.dat')), w_layout)
        t = d['th'] * (d['prs'] / co.P00) ** (co.RD / co.CP)
        mass = d['rho'] * dz
        parcel = {f'eps {eps:g}': bo.parcel(zh, d['prs'], t, d['qv'], mass, dz, 500.0, eps, g) for eps in (0.0, 8.0)}
        rows.append(dict(hour=(d['mtime'] - first_s) / HOUR, rain_total_mm=10.0 * d['rain'], rain_mm_day=d['prate'] * DAY,
                         cm1_cape_j_kg=d['cape'], cm1_cin_j_kg=d['cin'], hpbl_m=d['hpbl'],
                         w_max_below_16km_m_s=float(w['wmax'][zw < SHALLOW_M].max()),
                         w_max_below_60km_m_s=float(w['wmax'][zw < UPPER_M].max()),
                         cloud_top_m=cloud_top_m(d['qcfrac'], d['qifrac'], zh),
                         parcel_cape_j_kg=parcel['eps 0']['cape'], parcel_cin_j_kg=parcel['eps 0']['cin'],
                         parcel_lfc_m=parcel['eps 0']['lfc_m'], parcel_cape_eps8_j_kg=parcel['eps 8']['cape'],
                         temperature_k={f'{z / 1000:g} km': float(np.interp(z, zh, d['t'])) for z in HEIGHTS_M},
                         layers={layer: dict(temperature_c=tc, vapour_g_kg=qv) for (layer, tc), qv in
                                 zip(bb.layer_means(d['t'] - 273.15, d['rho'], zw).items(),
                                     bb.layer_means(d['qv'] * 1000.0, d['rho'], zw).values())}))
    return rows


def cloud_layers(qc, qi, zh) -> dict:
    """Where a snapshot's cloud (liquid and ice of at least 0.01 g/kg, as box_onset counts it) stands, column by column:
    the highest and the mean top of the cloud below 40 km, the share of columns with any cloud at 16-40 km, and the
    share whose cloud above 40 km stands apart, with none at 16-40 km beneath it."""
    zh = np.asarray(zh, float)
    cloudy = (np.asarray(qc) + np.asarray(qi)) >= ra.CLOUD_KG_KG
    below = zh < DEEP_M
    levels = np.arange(zh.size)[:, None]
    lower = cloudy & below[:, None]
    has = lower.any(axis=0)
    tops = zh[np.where(lower, levels, 0).max(axis=0)][has]
    middle = cloudy[(zh >= SHALLOW_M) & below].any(axis=0)
    high = cloudy[~below].any(axis=0)
    return dict(top_below_40km_highest_m=float(tops.max()) if has.any() else None,
                top_below_40km_mean_m=float(tops.mean()) if has.any() else None,
                cloud_16_40km_share=float(middle.mean()), apart_above_40km_share=float((high & ~middle).mean()))


def snapshot_rows(name: str, first_s: float = 0.0, last_s: float = np.inf) -> list:
    """box_onset's measures of every 3-D output of a case between model times first_s and last_s, with where its cloud
    stands (cloud_layers)."""
    case = ra.RUNS / name
    tap = json.loads((case / 'case.json').read_text())['configuration']['output_s']
    zh = ra.case_geometry(case)['zh']
    rows = []
    for n in sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat')):
        t = (n - 1) * tap
        if first_s - 1.0 <= t <= last_s + 1.0:
            d = ra.read_snapshot(case, n)
            rows.append(dict(hour=(t - first_s) / HOUR, **bo.snapshot(case, n, zh), clouds=cloud_layers(d['qc'], d['qi'], zh)))
    return rows


def day_summary(rows, snaps, a_h: float, b_h: float) -> dict:
    """One day of a run: its rain from the running total, the means (and maxima) of the diagnostics after its start up
    to its end, and the means of its snapshots' measures over the same outputs."""
    inside = [r for r in rows if a_h + 1e-6 < r['hour'] <= b_h + 1e-6]
    shots = [s for s in snaps if a_h + 1e-6 < s['hour'] <= b_h + 1e-6]
    mean = lambda xs: float(np.mean(xs)) if xs else None
    values = lambda key: [r[key] for r in inside if r[key] is not None]
    out = dict(hours=[a_h, b_h], outputs=len(inside), snapshots=len(shots), rain_mm_day=rain_over(rows, a_h, b_h),
               **{k: mean(values(k)) for k in ('cm1_cape_j_kg', 'cm1_cin_j_kg', 'hpbl_m', 'parcel_cape_j_kg',
                                                'parcel_cin_j_kg', 'parcel_lfc_m', 'parcel_cape_eps8_j_kg',
                                                'w_max_below_16km_m_s', 'w_max_below_60km_m_s', 'cloud_top_m')},
               strongest_updraft_max_m_s=max(values('w_max_below_60km_m_s'), default=None),
               cloud_top_highest_m=max(values('cloud_top_m'), default=None))
    if shots:
        out['snapshot'] = dict(
            w_p999_below_16km_m_s=mean([s['w_p999_below_16km_m_s'] for s in shots]),
            w_p999_below_16km_max_m_s=max(s['w_p999_below_16km_m_s'] for s in shots),
            w_p999_16_60km_m_s=mean([s['w_p999_16_60km_m_s'] for s in shots]),
            w_max_m_s=mean([s['w_max_m_s'] for s in shots]),
            cloud_top_max_m=mean([s['cloud_top_max_m'] for s in shots if s['cloud_top_max_m'] is not None]),
            cloud_tops_share={k: mean([s['cloud_tops_share'][k] for s in shots]) for k in shots[0]['cloud_tops_share']},
            clouds={k: mean([s['clouds'][k] for s in shots if s['clouds'][k] is not None]) for k in shots[0]['clouds']},
            cloud_top_below_40km_highest_m=max((s['clouds']['top_below_40km_highest_m'] for s in shots
                                                if s['clouds']['top_below_40km_highest_m'] is not None), default=None),
            columns={k: mean([s['columns'][k] for s in shots]) for k in shots[0]['columns']})
    return out


def quiet(day: dict, control: dict = CONTROL) -> dict:
    """The control's criteria (a test's control rule) on one day of a run: each criterion's value, its limit and
    whether it holds."""
    s = day['snapshot']
    values = dict(rain_mm_day=day['rain_mm_day'], strongest_updraft_m_s=day['w_max_below_60km_m_s'],
                  w_p999_below_16km_m_s=s['w_p999_below_16km_max_m_s'],
                  deep_top_share=s['cloud_tops_share']['above_40km'])
    tests = {k: (values[k], limit) for k, limit in control.items()}
    out = {k: dict(value=v, below=limit, holds=bool(v < limit)) for k, (v, limit) in tests.items()}
    out['quiet'] = all(v['holds'] for v in out.values())
    return out


def convects(day: dict, rule: dict = READING) -> bool:
    """A test's deep convection on one day of a run: enough rain, with strong updrafts or deep cells (cloud tops above
    40 km in the night's rule, strong updrafts aloft in the evening's)."""
    s = day['snapshot']
    deep = (('deep_top_share' in rule and s['cloud_tops_share']['above_40km'] >= rule['deep_top_share'])
            or ('w_p999_16_60km_m_s' in rule and s['w_p999_16_60km_m_s'] >= rule['w_p999_16_60km_m_s']))
    strong = day['w_max_below_60km_m_s'] >= rule['strongest_updraft_m_s']
    return bool(day['rain_mm_day'] >= rule['rain_mm_day'] and (strong or deep))


def reading(control_quiet: bool, finer: dict, quiet_spell: str = 'quiet sea nights') -> str:
    """A test's outcome under the rule stated before the runs, from the control's check and the finer grids' second
    days (grid -> (quiet, convects))."""
    if not control_quiet:
        return (f'the control does not reproduce the source box\'s {quiet_spell}, so the method is not valid and the '
                'finer grids are not read')
    if all(c for _, c in finer.values()):
        return f'storms start on the finer grids but not at 6 km: the {quiet_spell} come from the grid'
    if all(q for q, _ in finer.values()):
        return f'none of the three convect: the {quiet_spell} hold up across the grids'
    return 'in between: ' + '; '.join(f"{g} {'convects' if c else 'stays quiet' if q else 'neither quiet nor convecting'}"
                                      for g, (q, c) in finer.items())


def drift(rows) -> dict:
    """The domain-mean air's drift from a run's start: the temperature at HEIGHTS_M and the mass-weighted temperature (C)
    and vapour (g/kg) of the GCM's lowest four layers (box_budget.layer_means) at the start, after the first day and at
    the end, with the change over each day."""
    first = rows[0]
    middle = min(rows, key=lambda r: abs(r['hour'] - DAYS[1][1]))
    last = rows[-1]
    entry = lambda a, b, c: dict(start=round(a, 4), day_1=round(b, 4), end=round(c, 4), first_day_change=round(b - a, 4),
                                 second_day_change=round(c - b, 4))
    return dict(hours=[first['hour'], middle['hour'], last['hour']],
                temperature_k={k: entry(first['temperature_k'][k], middle['temperature_k'][k], last['temperature_k'][k])
                               for k in first['temperature_k']},
                layers={layer: {q: entry(first['layers'][layer][q], middle['layers'][layer][q], last['layers'][layer][q])
                                for q in first['layers'][layer]} for layer in first['layers']})


def refilled(name: str, first_day: float, last_day: float) -> float:
    """The vapour (mm/day) CM1's microphysics puts back in the levels whose air holds neither vapour nor condensate
    (under 1e-8 kg/kg each), where a prescribed drying finds nothing to take: its budget counts the refill as
    microphysics, which the water budget's closure leaves out with the condensate it makes, so it stands in the rest."""
    case = ra.RUNS / name
    dz = np.diff(ra.case_geometry(case)['zw'])
    layout = ra.grads_layout(case / 'cm1out_diag_s.ctl')[1]
    sums = []
    for path in sorted(case.glob('cm1out_diag_*_s.dat')):
        d = bb.read_diag(case, path, layout)
        if not first_day * DAY - 1.0 <= d['mtime'] <= last_day * DAY + 1.0:
            continue
        empty = (d['qv'] < 1e-8) & (sum(d[k] for k in bb.CONDENSATE) < 1e-8)
        sums.append(float((d['qvb_mp'] * d['rho'] * dz)[empty].sum()) * DAY)
    return float(np.mean(sums))


def water(name: str, span_days=(0.0, 2.0)) -> dict:
    """box_budget.water over the second day of a span and over all of it, its main terms (mm/day), with the refill
    of empty air (refilled) that the rest holds."""
    keys = ('span_days', 'evaporation_mm_day', 'rain_mm_day', 'storage_mm_day', 'rest_mm_day',
            'boundary_layer_over_evaporation', 'column_mm_day')
    out = {}
    for label, start in (('second day', span_days[0] + 1.0), ('both days', span_days[0])):
        w = bb.water(name, start, span_days[1])
        out[label] = {k: w[k] for k in keys}
        out[label]['refill_of_empty_air_mm_day'] = refilled(name, start, span_days[1])
        out[label]['rest_less_refill_mm_day'] = w['rest_mm_day'] - out[label]['refill_of_empty_air_mm_day']
    return out


def series(rows) -> dict:
    """A run's diagnostics at whole hours as compact columns (the run folders hold every output)."""
    keys = ('hour', 'rain_mm_day', 'w_max_below_16km_m_s', 'w_max_below_60km_m_s', 'cloud_top_m', 'cm1_cape_j_kg',
            'cm1_cin_j_kg', 'parcel_cape_j_kg', 'parcel_cin_j_kg', 'hpbl_m')
    fmt = lambda v: None if v is None else round(float(v), 3)
    hourly = [r for r in rows if abs(r['hour'] - round(r['hour'])) < 1e-6]
    return {k: [fmt(r[k]) for r in hourly] for k in keys}


def run_record(name: str) -> dict:
    """What a box ran: its grid, time steps and wall time, and the size of its folder, its inputs (under 1 MB) and
    everything it wrote."""
    from climate.crm.cm1_run import stats
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    progress = json.loads((case / 'progress.json').read_text())
    s = stats(case)
    written = sum(p.stat().st_size for p in case.iterdir() if p.is_file() and not p.is_symlink())
    dt = s['dt'][s['dt'] > 0.0]
    return dict(case=name, grid=record['grid'], runner=record['runner'], executable=record['build']['executable_sha256'],
                held=record['held'], initial=record['initial'],
                wall_hours=round(sum(seg['wall_s'] for seg in progress['segments']) / HOUR, 3),
                threads=sorted({seg['threads'] for seg in progress['segments']}),
                time_step_s=dict(mean=float(dt.mean()), min=float(dt.min()), max=float(dt.max())),
                folder_mb=round(written / 1e6, 1))


def analyse(test: str = 'night') -> dict:
    t = TESTS[test]
    grids, finer = {}, {}
    for label, name in t['grids'].items():
        rows, snaps = diag_rows(name), snapshot_rows(name)
        days = {d: day_summary(rows, snaps, a, b) for d, a, b in DAYS}
        grids[label] = dict(run=run_record(name), days=days, quiet=quiet(days['second day'], t['control']),
                            convects=convects(days['second day'], t['reading']), drift=drift(rows), water=water(name),
                            series=series(rows))
    src = t['source']
    first_s, last_s = src['first_day'] * DAY, src['last_day'] * DAY
    rows, snaps = diag_rows(src['case'], first_s, last_s), snapshot_rows(src['case'], first_s, last_s)
    source = dict(src, days={d: day_summary(rows, snaps, a, b) for d, a, b in DAYS}, drift=drift(rows),
                  water=water(src['case'], (src['first_day'], src['last_day'])), series=series(rows))
    source['quiet'] = quiet(source['days']['second day'], t['control'])
    control = grids['6 km']['quiet']['quiet']
    for label in ('3 km', '2 km'):
        finer[label] = (grids[label]['quiet']['quiet'], grids[label]['convects'])
    return dict(schema=SCHEMA, grids=grids, source=source, control_quiet=control,
                reading=reading(control, finer, t['quiet_spell']), **{t['gcm_key']: t['gcm_rain_mm_day']},
                control_rule=t['control'], reading_rule_limits=t['reading'],
                evidence=t['evidence'], reading_rule=t['reading_rule'])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--test', choices=sorted(TESTS), default='night', help='which held hour (default: the night)')
    args = parser.parse_args(argv)
    out = analyse(args.test)
    out['producer'] = dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'box_onset.py', HERE / 'box_budget.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / TESTS[args.test]['product']
    path.write_text(json.dumps(out, indent=1) + '\n')
    fmt = lambda v, f='.2f': '-' if v is None else format(v, f)
    for label, e in [*out['grids'].items(), ('source', out['source'])]:
        for day, d in e['days'].items():
            s = d.get('snapshot', {})
            print(f"{label:6s} {day}: rain {fmt(d['rain_mm_day'])} mm/day; strongest updraft {fmt(d['w_max_below_60km_m_s'])} "
                  f"(below 16 km {fmt(d['w_max_below_16km_m_s'])}, max {fmt(d['strongest_updraft_max_m_s'])}) m/s; "
                  f"w p99.9 <16 km {fmt(s.get('w_p999_below_16km_m_s'))}; top {fmt((d['cloud_top_m'] or 0) / 1000, '.1f')} km, "
                  f"snapshot top {fmt((s.get('cloud_top_max_m') or 0) / 1000, '.1f')} km, deep share "
                  f"{fmt(s.get('cloud_tops_share', {}).get('above_40km'), '.3f')}; CAPE {fmt(d['cm1_cape_j_kg'], '.0f')} "
                  f"CIN {fmt(d['cm1_cin_j_kg'], '.1f')}; parcel {fmt(d['parcel_cape_j_kg'], '.0f')}/"
                  f"{fmt(d['parcel_cin_j_kg'], '.1f')}; hpbl {fmt(d['hpbl_m'], '.0f')} m")
        if 'quiet' in e:
            print(f"   quiet: {e['quiet']['quiet']}", {k: v['holds'] for k, v in e['quiet'].items() if k != 'quiet'})
    print('reading:', out['reading'])
    print(f'wrote {path.relative_to(HERE.parents[1])}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
