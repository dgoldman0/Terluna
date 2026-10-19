"""Rain rates on the land in the three target climates, from the GCM: how much of the land's time is wet and how hard
it rains, in the 3-day means of the two runs behind the climate targets and in the design case's 6.5-hour years.

    climate/gcm/.venv/bin/python -m atmosphere.rain.rates      # about ten seconds; needs netCDF4 and the GCM runs
    # -> atmosphere/rain/results/rain_rates.json

The climate targets (research/studies/climate_targets, bound through targets.py) read the design case,
A28_corrected_970 near 295 K, and the 20-year test of the longwave refitted with ice optics, A28_corrected_970_lw_ice
near 307 K, over years 10-19 of each, and set each target between them, linear in the ground's temperature: the lock
is the design case, the hot world the test and the resort, at 301 K, 52% of the way between them. This product reads
the same output's total precipitation, `pr`, a true mean over each output window: ExoPlaSim 3.4.2 adds each step's
rain to its sum (outmod.f90, line 2517), divides the sum by the window's steps as it writes it (line 539) and resets
it for the next window (line 2395). The design case's 6.5-hour years, the two after its year 19, averaged eleven at a
time match its own 3-day means over years 10-19 in distribution, within 0.8% in the wet share and every quantile to
the 99th. Each cell is weighed by its land area, its land fraction times its area, as the targets weigh the land's
rain. A mean is wet at 0.1 mm/h or more. For each run it gives the wet share of the land's time and the wet means'
quantiles; for each target, those of the two runs interpolated as the targets interpolate, beside the climate
targets' land rain with the cut that target's own dimming makes.

A28_corrected_970_eddies carried the design case on from its year 19 for two years, writing every 6.5 hours, eleven
means to each 3-day mean. Its 6.5-hour means set beside the 3-day means their elevens make show what the 3-day means
smooth. Each target's 6.5-hour statistics take its 3-day ones times the design case's ratio of 6.5-hour to 3-day at
each quantile, and its wet share times the design case's ratio of wet shares: the design case's short weather, carried
to the warmer targets.

Every mean is a cell's mean, over a T21 cell about 170 km across; rain at a point under a shower falls harder and
for shorter times.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from atmosphere.rain import targets as climate_targets

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.rain-rates/1'
OUT = HERE / 'results' / 'rain_rates.json'
CODE = ('atmosphere/rain/rates.py', 'atmosphere/rain/targets.py')
RUNS_DIR = ROOT / 'climate' / 'gcm' / 'runs'
BOUNDARY = ROOT / 'climate' / 'gcm' / 'products' / 'moon_28pct_water_lakes_gaussian_T21.nc'
DESIGN, TEST, SHORT = 'A28_corrected_970', 'A28_corrected_970_lw_ice', 'A28_corrected_970_eddies'
YEARS, SHORT_YEARS = (10, 19), (0, 1)
TARGETS = climate_targets.KEYS          # the hot world, the resort and the lock (targets.py)
WET_MM_H = 0.1
QUANTILES = (0.1, 0.5, 0.9, 0.99)
MM_H_PER_M_S = 3.6e6

EVIDENCE = ('Statistics of the GCM\'s output: total precipitation in 3-day and 6.5-hour means over the land of a T21 '
            'grid, weighted by land area. The targets\' 3-day statistics interpolate two runs of different physics '
            'by the ground\'s temperature, as the climate targets do; their 6.5-hour statistics carry the design '
            'case\'s short weather to them by quantile ratios. Both runs take the design sunlight without the '
            'dimming a target needs; the climate targets\' cut of the land\'s mean rain for that dimming is carried '
            'beside each target\'s statistics.')
READING_RULE = ('Rates are mm/h of rain at the ground, means over a T21 cell and an output window. A mean is wet at '
                'wet_threshold_mm_h or more; wet.share is the share of the land\'s time in wet means and wet.p10-p99 '
                'and top are the wet means\' quantiles and largest, each weighted by land area. runs gives the two '
                'runs\' 3-day means over years 10-19; short gives the design case\'s 6.5-hour years, the 3-day means '
                'their elevens make and the ratio of the first to the second at each statistic; targets gives each '
                'target\'s 3-day statistics, interpolated with its share_of_the_way from the climate targets '
                '(held_in: their targets, or their resort band\'s rows), and its 6.5-hour statistics, those times the '
                'design case\'s ratios. land_rain_with_the_dimmings_own_cut_mm_day is the climate targets\' land rain '
                'with the dimming\'s own cut, [low, high], and dimmings_own_cut_percent the cut it takes from the '
                'mean, [least, most].')
ASSUMPTIONS = [
    'The land is every cell\'s land fraction, weighted by the cell\'s area, as the climate targets weigh the land\'s '
    'rain; the T21 grid cannot resolve the coasts.',
    'Each target between the two runs takes their statistics linear in the ground\'s temperature, as the climate '
    'targets take their climate; the hot world is the test run itself, the lock the design case and the resort the '
    'climate targets\' candidate at 301 K, where the author set it on 18 October.',
    'The short weather is the design case\'s, near 295 K, over two years: the warmer targets take its ratios of '
    '6.5-hour to 3-day rain.',
    'Both runs take the design sunlight, without the dimming the resort and the lock need. That dimming cuts the '
    'land\'s mean rain by dimmings_own_cut_percent in the climate targets; how the cut divides between how often and '
    'how hard it rains is open, and spread evenly over the rates it would lower each by that share.',
    'A cell\'s mean rain spreads a storm over 170 km and over the window: a point under a shower sees heavier, '
    'briefer rain than any of these means.',
]


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def output(run, year):
    return RUNS_DIR / run / 'model' / f'MOST.{year:05d}.nc'


def boundary():
    import netCDF4
    with netCDF4.Dataset(BOUNDARY) as d:
        return dict(lat=np.asarray(d['lat'][:], float), land=np.asarray(d['land_fraction'][:], float),
                    area=np.asarray(d['cell_area_fraction'][:], float))


def window_hours(run):
    """The length of one output mean (h) from the run's own settings."""
    model = json.loads((RUNS_DIR / run / 'progress.json').read_text())['configuration']['model']
    return model['timestep_min'] * model['steps_per_lunar_day'] / model['writes_per_lunar_day'] / 60.0


def read_rain(run, years, b, files):
    """The run's rain (mm/h) by mean and cell, and each cell's land area in the output's row order."""
    import netCDF4
    out, lat = [], None
    for year in range(years[0], years[1] + 1):
        path = output(run, year)
        files[str(path.relative_to(ROOT))] = sha256(path)
        with netCDF4.Dataset(path) as d:
            lat = np.asarray(d['lat'][:], float)
            out.append(np.asarray(d['pr'][:], float) * MM_H_PER_M_S)
    order = [int(np.argmin(np.abs(b['lat'] - x))) for x in lat]
    return np.concatenate(out), b['land'][order] * b['area'][order]


def wquantile(x, w, q):
    """The quantile q (0-1) of values x under weights w, linear between the weights' midpoints, as the climate
    targets take theirs."""
    i = np.argsort(x, kind='stable')
    x, w = x[i], w[i]
    c = np.cumsum(w)
    return float(np.interp(q * c[-1], c - 0.5 * w, x))


def statistics(rain, land_area):
    """The land's mean rain, the wet share of its time and the wet means' quantiles."""
    keep = land_area > 0
    x = rain[:, keep].ravel()
    w = np.broadcast_to(land_area[keep], rain[:, keep].shape).ravel()
    wet = x >= WET_MM_H
    out = dict(land_rain_mm_day=float((w * x).sum() / w.sum() * 24.0),
               wet=dict(share=float(w[wet].sum() / w.sum()),
                        **{f'p{round(100 * q)}': wquantile(x[wet], w[wet], q) for q in QUANTILES},
                        top=float(x[wet].max())))
    return out


def ratios(fine, coarse):
    return {k: fine['wet'][k] / coarse['wet'][k] for k in fine['wet']}


def build():
    data, binding = climate_targets.load()
    read = climate_targets.read(data)
    files = {str(BOUNDARY.relative_to(ROOT)): sha256(BOUNDARY)}
    b = boundary()
    runs = {}
    for run, role in ((DESIGN, 'the design case'), (TEST, 'the 20-year test, longwave refitted, ice optics')):
        rain, land_area = read_rain(run, YEARS, b, files)
        runs[run] = dict(role=role, years=list(YEARS), means=int(rain.shape[0]), window_hours=window_hours(run),
                         **statistics(rain, land_area))
    rain, land_area = read_rain(SHORT, SHORT_YEARS, b, files)
    hours = window_hours(SHORT)
    per = int(round(window_hours(DESIGN) / hours))
    whole = rain.shape[0] // per * per
    fine = statistics(rain[:whole], land_area)
    coarse = statistics(rain[:whole].reshape(-1, per, *rain.shape[1:]).mean(axis=1), land_area)
    short = dict(run=SHORT, years=list(SHORT_YEARS), window_hours=hours, means_to_a_3_day_mean=per,
                 means=int(whole), six_and_a_half_hours=fine, three_day_means_of_them=coarse,
                 ratio=ratios(fine, coarse), blocks_over_the_design_cases_own=ratios(coarse, runs[DESIGN]))
    out_targets = {}
    for key in TARGETS:
        target = read[key]
        f = target['share_of_the_way']
        a, c = runs[DESIGN], runs[TEST]
        three = dict(land_rain_mm_day=a['land_rain_mm_day'] + f * (c['land_rain_mm_day'] - a['land_rain_mm_day']),
                     wet={k: a['wet'][k] + f * (c['wet'][k] - a['wet'][k]) for k in a['wet']})
        six = dict(wet={k: three['wet'][k] * short['ratio'][k] for k in three['wet']})
        cut = target['land_rain_with_the_dimmings_own_cut_mm_day']
        out_targets[key] = dict(
            name_in_targets=target['name_in_targets'], held_in=target['held_in'], ground_k=target['ground_k'],
            share_of_the_way=f, three_day_means=three, six_and_a_half_hour_means=six,
            land_rain_with_the_dimmings_own_cut_mm_day=cut,
            dimmings_own_cut_percent=sorted(100.0 * (1.0 - x / target['land_rain_mm_day']) for x in cut))
    return dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, assumptions=ASSUMPTIONS,
                wet_threshold_mm_h=WET_MM_H, runs=runs, short=short, targets=out_targets,
                producer=dict(domain='atmosphere', files={p: sha256(ROOT / p) for p in CODE}, **binding, read=files))


def main(argv=None) -> int:
    result = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1, allow_nan=False) + '\n')
    s = result['short']
    print(f"design case, 6.5-hour means over 3-day: wet share {s['ratio']['share']:.2f}, median "
          f"{s['ratio']['p50']:.2f}, 99th percentile {s['ratio']['p99']:.2f}, top {s['ratio']['top']:.2f}")
    for key, t in result['targets'].items():
        three, six = t['three_day_means']['wet'], t['six_and_a_half_hour_means']['wet']
        cut = t['dimmings_own_cut_percent']
        print(f"{key:9s} land {t['three_day_means']['land_rain_mm_day']:.2f} mm/day; wet {three['share']:.1%} of "
              f"3-day means, median {three['p50']:.2f}, p99 {three['p99']:.2f} mm/h; 6.5-hour wet {six['share']:.1%}, "
              f"median {six['p50']:.2f}, p99 {six['p99']:.2f} mm/h; the dimming's cut {cut[0]:.1f}-{cut[1]:.1f}%")
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
