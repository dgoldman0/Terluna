"""A forced box's lunar days one after another, across its case and the cases run on from its restarts: the rain by
quarter of the lunar day, CM1's CAPE and inhibition by quarter, the domain-mean air's temperature at chosen heights
and its drift, the water budget's closure, and the share of the falling water evaporated with the constant gamma at
which PlaSim evaporates that share at the box's place.

    climate/gcm/.venv/bin/python -m climate.crm.box_days box_near_side_sea box_near_side_sea_day3
                                                                     # ../results/crm/box_days_box_near_side_sea.json

The cases are read in order, each from its own outputs after the last one's end, as a case set up by cm1_run's
setup_restart carries its source's numbering on. Each lunar day is a span between whole restart times (0-29.5,
29.5-59.0 and 59.0-88.5 days for boxes started at sunrise): its rain is the domain-mean rate at each 3-hourly output, as
the box's water budget reckons it, and its quarters follow the site's hour angle. The water budget (box_budget.water)
and the share evaporated (box_analysis.microphysics_budget) are taken over a span that runs to its case's end; gamma
comes from rain_evaporation.fit on the pooled GCM shares in ../results/gcm/rain_evaporation.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from climate.crm import box_budget as bb
from climate.crm import ring_analysis as ra

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-box-days/1'
DAY = 86400.0
SPANS = ((0.0, 29.5), (29.5, 59.0), (59.0, 88.5))
HEIGHTS_M = (40000.0, 60000.0)
QUARTERS = (('sunrise-noon', -90.0, 0.0), ('noon-sunset', 0.0, 90.0), ('sunset-midnight', 90.0, 180.0),
            ('midnight-sunrise', -180.0, -90.0))
REGIONS = {'sea': 'near_side_sea', 'land': 'cm1_box'}       # rain_evaporation's place for each surface


def outputs(cases) -> list:
    """Every diagnostic output of the cases in order, each case's after the last one's end: model day, hour angle at
    the site, rain (mm/day), CAPE and inhibition (J/kg), the temperature at HEIGHTS_M (K) and the case it came from."""
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    rows, last = [], -np.inf
    for case in cases:
        folder = ra.RUNS / case
        record = json.loads((folder / 'case.json').read_text())
        start, lon = record['configuration']['start_hour_angle_deg'], record['configuration']['site']['lon_deg']
        zh = ra.case_geometry(folder)['zh']
        layout = ra.grads_layout(folder / 'cm1out_diag_s.ctl')[1]
        mine = []
        for path in sorted(folder.glob('cm1out_diag_*_s.dat')):
            d = bb.read_diag(folder, path, layout)
            if d['mtime'] <= last + 1.0:
                continue
            mine.append(dict(day=d['mtime'] / DAY, hour=local_hour_angle(start, lon, d['mtime'], solar_day_s()),
                             rain=d['prate'] * DAY, cape=d['cape'], cin=d['cin'], case=case,
                             t_k=[float(np.interp(z, zh, d['t'])) for z in HEIGHTS_M]))
        rows += mine
        last = max(last, mine[-1]['day'] * DAY) if mine else last
    return rows


def quarters(rows, key) -> list:
    h = np.array([r['hour'] for r in rows])
    v = np.array([r[key] for r in rows], float)
    return [float(v[(h >= a) & (h < b)].mean()) if ((h >= a) & (h < b)).any() else None for _, a, b in QUARTERS]


def span_summary(rows, span) -> dict:
    """One lunar day: its rain and CM1's CAPE and inhibition by quarter, and the temperature at each height at its
    first and last output, both at sunrise, with their change per day and the day's mean."""
    inside = [r for r in rows if span[0] - 1e-6 <= r['day'] <= span[1] + 1e-6]
    if len(inside) < 2:
        return None
    day = np.array([r['day'] for r in inside])
    t = np.array([r['t_k'] for r in inside])
    return dict(span_days=[float(day[0]), float(day[-1])], outputs=len(inside), cases=sorted({r['case'] for r in inside}),
                rain_mm_day=float(np.mean([r['rain'] for r in inside])), rain_by_quarter_mm_day=quarters(inside, 'rain'),
                cape_by_quarter_j_kg=quarters(inside, 'cape'), inhibition_by_quarter_j_kg=quarters(inside, 'cin'),
                temperature_k={f'{z / 1000:g} km': dict(first=float(t[0, k]), last=float(t[-1, k]),
                                                        change_k_day=float((t[-1, k] - t[0, k]) / (day[-1] - day[0])),
                                                        mean=float(t[:, k].mean()))
                               for k, z in enumerate(HEIGHTS_M)})


def evaporated(case: str, from_day: float, surface: str) -> dict:
    """The share of the falling water evaporated over a case from from_day to its end, with the gamma at which
    PlaSim's pooled runs evaporate it at the box's place."""
    from climate.crm.box_analysis import microphysics_budget
    from climate.gcm import rain_evaporation as re_
    out = microphysics_budget(ra.RUNS / case, from_day)
    pooled = json.loads(re_.OUT.read_text())['pooled']
    out['gamma'] = re_.fit(pooled, REGIONS[surface], out['share_of_falling_water_evaporated'])
    out['gamma_region'] = REGIONS[surface]
    return out


def analyse(cases) -> dict:
    rows = outputs(cases)
    record = json.loads((ra.RUNS / cases[0] / 'case.json').read_text())
    surface = record['configuration']['site'].get('surface', 'land')
    days = [s for s in (span_summary(rows, span) for span in SPANS) if s is not None]
    for before, now in zip(days, days[1:]):                         # the lunar days' means one after another
        for k, v in now['temperature_k'].items():
            v['mean_change_from_day_before_k_day'] = (v['mean'] - before['temperature_k'][k]['mean']) / (
                np.mean(now['span_days']) - np.mean(before['span_days']))
    # the water budget and the share evaporated over each span that a single case runs to its end
    for summary in days:
        after = {r['case'] for r in rows if summary['span_days'][0] + 1e-6 < r['day'] <= summary['span_days'][1] + 1e-6}
        case = after.pop() if len(after) == 1 else None
        last = max(r['day'] for r in rows if r['case'] == case) if case else None
        if case and abs(last - summary['span_days'][1]) < 1e-6:
            w = bb.water(case, summary['span_days'][0])
            summary['water_mm_day'] = {k: w[k] for k in ('evaporation_mm_day', 'rain_mm_day', 'storage_mm_day',
                                                        'rest_mm_day', 'boundary_layer_over_evaporation')}
            summary['water_mm_day']['by_process'] = w['column_mm_day']
            summary['falling_water'] = evaporated(case, summary['span_days'][0], surface)
    start = rows[0]
    return dict(schema=SCHEMA, cases=list(cases), surface=surface, heights_m=list(HEIGHTS_M),
                start=dict(day=start['day'], temperature_k={f'{z / 1000:g} km': start['t_k'][k]
                                                            for k, z in enumerate(HEIGHTS_M)}),
                series=dict(day=[round(r['day'], 4) for r in rows],
                            temperature_k={f'{z / 1000:g} km': [round(r['t_k'][k], 3) for r in rows]
                                           for k, z in enumerate(HEIGHTS_M)}),
                lunar_days=days,
                evidence=('CM1 domain means at each 3-hourly output of one box run on from its own restarts, one '
                          'realization; the gamma from PlaSim\'s pooled runs at the box\'s place.'),
                reading_rule=('Each lunar day: rain in mm/day (domain-mean rate at the outputs), by quarter of the '
                              'lunar day at the site; CM1\'s CAPE and inhibition (J/kg, the mean of its columns\'); the '
                              'domain-mean temperature at each height at the first and last output (both at sunrise) '
                              'and its change per day between them, its mean over the lunar day and that mean\'s change '
                              'per day from the lunar day before; the first output of all is the GCM\'s mean air, '
                              'where the box starts; where one case runs over the whole day to its end, '
                              'its vapour budget (box_budget.water: rest closes it) and the share of the water leaving '
                              'its cloud that evaporated on the way down, with the gamma at which the GCM\'s pooled runs '
                              'evaporate that share at the box\'s place (null outside them).'))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('cases', nargs='+', help='the case, then the cases run on from its restarts, in order')
    args = parser.parse_args(argv)
    out = analyse(args.cases)
    out['producer'] = dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'box_budget.py',
                                                              HERE / 'box_analysis.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f'box_days_{args.cases[0]}.json'
    path.write_text(json.dumps(out, indent=1) + '\n')
    fmt = lambda v: ' '.join(f'{x:7.2f}' for x in v)
    for d in out['lunar_days']:
        temps = ', '.join(f"{k} {v['first']:.1f}->{v['last']:.1f} K ({v['change_k_day']:+.3f} K/day), mean "
                          f"{v['mean']:.1f}" for k, v in d['temperature_k'].items())
        print(f"days {d['span_days'][0]:.1f}-{d['span_days'][1]:.1f}: rain {d['rain_mm_day']:.2f} by quarter "
              f"{fmt(d['rain_by_quarter_mm_day'])}; inhibition {fmt(d['inhibition_by_quarter_j_kg'])}; {temps}")
        if 'water_mm_day' in d:
            w, f = d['water_mm_day'], d['falling_water']
            print(f"   water: evaporation {w['evaporation_mm_day']:.2f}, rain {w['rain_mm_day']:.2f}, storage "
                  f"{w['storage_mm_day']:.2f}, rest {w['rest_mm_day']:.3f}; falling water evaporated "
                  f"{f['share_of_falling_water_evaporated']:.3f}, gamma {f['gamma']}")
    print(f'wrote {path.relative_to(HERE.parents[1])}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
