"""What the 3-day means leave out of the waters' warmth: the corrected design climate's year 10 written every 6.5
hours, against its own 3-day means.

The waters' forcing and map read the design years as 3-day means and rebuild what a mean hides: the speed of each
window's mean wind vector is floored at PlaSim's least wind and raised by one factor borrowed from an earlier run
(forcing.wind_raise), the air and the sky run linearly between the windows' centres, and the sunlight is shaped
within the window. The record `corrected` of wind_record.py is the design run's year 10 again, written every 13
half-hour steps, eleven records to each window. This runner puts the two side by side on every cell:

- the record against the design year's own file: each field the forcing reads, the eleven aggregated as the field
  itself is (their mean; the last for the clear-sky snapshots; the first for the humidity, which the run writes as
  each interval's opening value), against the 3-day record, and the two runs' year-end restarts, so that the two are
  shown to be one weather;
- the humidity at a year's ends, which no file holds: the fill forcing.window_humidity gives a year file's empty
  first record, against the record's values 6.5 hours inside the year;
- the map's waters (water_map.run: still water 1, 2 and 4 m deep, water mixed to 10 m, the 2 m shore) stepped through
  the year under each of these forcings, all by one Sun and one set of clear-sky tables (the record's own, 1,320
  values a cell):

      six_hourly         the 6.5-hour records as they are, the wind floored and unraised: the reference;
      means_raised       the 3-day means with the forcing's raise: the forcing's and the map's method;
      means_unraised     the 3-day means, the wind floored and unraised;
      means_true_speed   the 3-day means, each window taking the mean of its eleven floored speeds: what a perfect
                         account of the wind's speed alone would give;
      means_cell_factor  the 3-day means raised by each cell's own factor from this record (wind_record.floored);
      means_cell_excess  the 3-day means with each cell's own excess added in quadrature
                         (wind_record.quadrature_excess).

  The statistics are over the year after the map's two lunar days of spin-up, and each pass is reported as its
  difference from six_hourly.

The year is one of the ten design years, so its statistics are noisier than the map's; the differences between passes
are paired, the same weather under each. The two cell-by-cell raises are fitted on the year they are tried on.

    /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.littoral_community.inside_windows  # minutes

writes results/inside_windows.json (schema terluna.research.littoral-inside-windows/1) with the sha256 of every input.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology.moisture_map import cell_areas_km2
from research.studies.littoral_community import forcing as F
from research.studies.littoral_community import water_map as WM
from research.studies.littoral_community import wind_record as R
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'inside_windows.json'
SCHEMA = 'terluna.research.littoral-inside-windows/1'
PARENT_YEAR = 10                              # the design year the record repeats
PER = R.WRITES_PER_WINDOW                     # records to a 3-day window
FINE_STEPS = R.WRITE_STEPS                    # half-hour steps to a record
FINE_H = FINE_STEPS * 0.5
MEANS = ('tas', 'ts', 'rss', 'rls', 'alb', 'ps', 'uas', 'vas', 'czen', 'ta_low', 'ua_low', 'va_low')
LAST = ('rsscs',)                             # single steps at each interval's end
FIRST = ('hus_low',)                          # each interval's opening value
REFERENCE = 'six_hourly'
PASSES = ('six_hourly', 'means_raised', 'means_unraised', 'means_true_speed', 'means_cell_factor', 'means_cell_excess')
REPORT = ('still_2m', 'still_4m', 'mixed_10m', 'shore_2m_x0.3')
KEYS = ('mean_c', 'p99_c', 'hours_above_36')


def fine_file():
    return R.record_files('corrected')[0]


def parent_file(year=PARENT_YEAR):
    return F.RUNS / F.RUN / 'model' / f'MOST.{year:05d}.nc'


def restarts():
    """The record's year-end restart and the design run's at the same step."""
    return (F.RUNS / R.CORRECTED / 'model' / 'MOST_REST.00000', F.RUNS / F.RUN / 'model' / f'MOST_REST.{PARENT_YEAR:05d}')


def aggregate(fine, per=PER):
    """The record's fields brought to 3-day windows, each as the field itself is: the mean of a window's records, the
    last for a snapshot at the interval's end, the first for an interval's opening value."""
    out = {}
    for k in MEANS + LAST + FIRST:
        a = np.asarray(fine[k])
        n = a.shape[0] // per
        blocks = a[:n * per].reshape(n, per, *a.shape[1:]).astype(float)
        out[k] = blocks.mean(1) if k in MEANS else blocks[:, -1] if k in LAST else blocks[:, 0]
    return out


def agreement(fine, coarse, per=PER):
    """Each field's difference between the aggregated record and the 3-day record: the largest and the root mean
    square, each over the field's largest magnitude in the 3-day record (scale, in the field's own unit); the
    humidity without the year's first window, whose opening value neither file holds."""
    out = {}
    for k, a in aggregate(fine, per).items():
        b = np.asarray(coarse[k], dtype=float)
        if k in FIRST:
            a, b = a[1:], b[1:]
        scale = max(float(np.abs(b).max()), 1e-30)
        out[k] = dict(largest=float(np.abs(a - b).max() / scale), rms=float(np.sqrt(np.mean((a - b) ** 2)) / scale),
                      scale=scale)
    return out


def humidity_at_the_years_ends(fine, coarse, next_year_second):
    """The fill of a year file's empty first record against the record's values 6.5 hours inside the year, over the
    cells: at the year's start the fill is the next 3-day record (the series' first year), and at its end the mean of
    the last window's opening value and the next year's second record (every later year). Each as the root mean
    square and the mean of (fill - record) / record."""
    q_fine, q = np.asarray(fine['hus_low'], float), np.asarray(coarse['hus_low'], float)
    def apart(fill, near):
        rel = (fill - near) / near
        return dict(rms=round(float(np.sqrt(np.mean(rel ** 2))), 4), mean=round(float(rel.mean()), 4))
    return dict(as_written_over_the_next_record=round(float((q[0] / q[1]).mean()), 4),
                start_fill_against_6h5_inside=apart(q[1], q_fine[1]),
                end_fill_against_6h5_inside=apart(0.5 * (q[-1] + np.asarray(next_year_second, float)), q_fine[-1]),
                a_window_of_real_change=apart(q[2], q[1]))


def build_passes(fine, coarse, sun):
    """Every pass's windows, with the wind's cell-by-cell measures at the lowest level."""
    lat2, _, noon2 = WM.cell_suns(fine, sun)
    tables = WM.clear_tables(WM.clear_sky_end(fine), lat2, noon2, sun['period_h'], FINE_H)
    speed, vector = R.windows(fine['ua_low'], fine['va_low'], PER)
    inside, outside, factor = R.floored(speed, vector)
    excess = R.quadrature_excess(inside, outside)
    flat = lambda a: a.reshape(a.shape[0], -1)
    means = lambda **kw: WM.windows(coarse, sun, tables=tables, **kw)
    wins = dict(six_hourly=WM.windows(fine, sun, FINE_H, FINE_STEPS, wind_raise=1.0, tables=tables),
                means_raised=means(), means_unraised=means(wind_raise=1.0), means_true_speed=means(level=flat(inside)),
                means_cell_factor=means(level=flat(outside * factor)),
                means_cell_excess=means(level=flat(np.sqrt(outside ** 2 + excess ** 2))))
    return wins, dict(factor=factor, excess=excess)


def differences(stats, names, area, sea2, lat):
    """Each pass less the reference, for the reported cases and keys: area-weighted mean and root mean square over the
    water cells and the land cells, and the mean by latitude band."""
    shape = sea2.shape
    lat2 = np.abs(lat)[:, None] * np.ones(shape)
    out = {}
    for name in PASSES:
        if name == REFERENCE:
            continue
        rows = {}
        for case in REPORT:
            n = names.index(case)
            row = {}
            for key in KEYS:
                d = (stats[name][key][n] - stats[REFERENCE][key][n]).reshape(shape)
                entry = {}
                for label, mask in (('water_cells', sea2), ('land_cells', ~sea2)):
                    w = area[mask]
                    entry[label] = dict(mean=round(float((d[mask] * w).sum() / w.sum()), 3),
                                        rms=round(float(np.sqrt((d[mask] ** 2 * w).sum() / w.sum())), 3),
                                        by_band=[round(float((d * area)[mask & (lat2 >= lo) & (lat2 < hi)].sum()
                                                             / area[mask & (lat2 >= lo) & (lat2 < hi)].sum()), 3)
                                                 for lo, hi in WM.LATITUDE_BANDS])
                row[key] = entry
            rows[case] = row
        out[name] = rows
    return out


def levels(stats, names, area, sea2, lat):
    """The reference's and the map's method's own band means, for the reported cases and keys."""
    shape = sea2.shape
    lat2 = np.abs(lat)[:, None] * np.ones(shape)
    out = {}
    for name in (REFERENCE, 'means_raised'):
        out[name] = {}
        for case in REPORT:
            n = names.index(case)
            out[name][case] = {key: {label: [round(float((stats[name][key][n].reshape(shape) * area)[mask & (lat2 >= lo) & (lat2 < hi)].sum()
                                                         / area[mask & (lat2 >= lo) & (lat2 < hi)].sum()), 2)
                                             for lo, hi in WM.LATITUDE_BANDS]
                                     for label, mask in (('water_cells', sea2), ('land_cells', ~sea2))} for key in KEYS}
    return out


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    import netCDF4
    started = time.time()
    fine, coarse = F.read_years([fine_file()]), F.read_years([parent_file()])
    with netCDF4.Dataset(parent_file(PARENT_YEAR + 1)) as nc:
        next_second = np.array(nc['hus'][1, -1], dtype=float)
    lat, lon = coarse['lat'], coarse['lon']
    shape = (lat.size, lon.size)
    sea2 = coarse['lsm'] <= 0.5
    area = cell_areas_km2(lat, lon)
    agree = agreement(fine, coarse)
    mine, theirs = restarts()
    restart_identical = bool(mine.exists() and theirs.exists() and F.sha256(mine) == F.sha256(theirs))
    print('agreement', {k: f"{v['largest']:.1e}" for k, v in agree.items()}, 'restart identical', restart_identical, flush=True)
    sun = WM.global_sun(coarse)
    wins, wind = build_passes(fine, coarse, sun)
    names = [c[0] for c in WM.CASES] + ['sea_surface']
    stats = {}
    for name in PASSES:
        spin = F.SPIN_WINDOWS * (PER if name == REFERENCE else 1)
        stats[name] = WM.summarise(WM.run(wins[name], sun, sea2.ravel(), spin_windows=spin))
        print(name, 'done', round(time.time() - started, 1), 's', flush=True)
    mean = lambda a, m: round(float((a * area)[m].sum() / area[m].sum()), 3)
    level_mean = {name: dict(water=mean(wins[name]['level'].mean(0).reshape(shape), sea2),
                             land=mean(wins[name]['level'].mean(0).reshape(shape), ~sea2)) for name in PASSES}
    forcing = json.loads(WM.FORCING.read_text())
    sites = {}
    for site_name, site in forcing['sites'].items():
        c = site['j'] * lon.size + site['i']
        sites[site_name] = dict(kind=site['kind'], lat=site['lat'],
                                wind_factor=round(float(wind['factor'][site['j'], site['i']]), 3),
                                wind_excess_m_s=round(float(wind['excess'][site['j'], site['i']]), 3),
                                passes={name: {case: {key: round(float(stats[name][key][names.index(case), c]), 2) for key in KEYS}
                                               for case in ('still_2m', 'still_4m')} for name in PASSES})
    files_used = ['research/studies/littoral_community/inside_windows.py', 'research/studies/littoral_community/water_map.py',
                  'research/studies/littoral_community/forcing.py', 'research/studies/littoral_community/wind_record.py',
                  'biosphere/ecology/lake.py', 'biosphere/ecology/moisture_map.py']
    inputs = [fine_file(), parent_file(), parent_file(PARENT_YEAR + 1)]
    record = dict(
        schema=SCHEMA,
        producer=dict(lane='research', runner='research/studies/littoral_community/inside_windows.py',
                      files={f: F.sha256(ROOT / f) for f in files_used}, constants=constants_used(files_used)),
        evidence=("Model results on one year of the corrected design climate, its year 10, read twice: as the design "
                  "run's own 3-day means and as the same year run again and written every 6.5 hours, the two ending "
                  "on one restart. The waters are the map's mixed slabs under each forcing. The differences between "
                  'passes are paired, the same weather under each; the year is one of ten, and the 6.5-hour means '
                  'are themselves means. The cell-by-cell raises are fitted on the year they are tried on.'),
        reading_rule=("agreement gives each field's largest and root-mean-square difference between the record's "
                      "aggregated values and the 3-day record, over the field's largest magnitude (scale). passes gives each pass less six_hourly "
                      'for a case and a statistic: the area-weighted mean and root mean square over the water cells '
                      'and the land cells, and the mean by band of |latitude| (settings.bands). mean_c and p99_c are '
                      'in K, hours_above_36 in hours of each lunar day. levels gives the two named passes\' own band '
                      "means. wind_at_the_level is each pass's mean wind at the lowest level, m/s. A positive "
                      'difference means the pass runs warmer than the 6.5-hour forcing.'),
        settings=dict(record=R.CORRECTED, record_year=0, design_run=F.RUN, design_year=PARENT_YEAR, records_per_window=PER,
                      record_h=FINE_H, window_h=F.WINDOW_H, steps=dict(six_hourly=FINE_STEPS, means=F.STEPS_PER_WINDOW),
                      spin_windows=F.SPIN_WINDOWS, level_floor_m_s=F.LEVEL_FLOOR_MS, wind_raise=F.wind_raise(),
                      cases=[dict(name=c[0], depth_m=c[1], exchange_per_day=c[2]) for c in WM.CASES],
                      bands=[list(b) for b in WM.LATITUDE_BANDS], reference=REFERENCE, passes=list(PASSES)),
        sun=sun,
        checks=dict(agreement=agree, restart_identical=restart_identical,
                    humidity_at_the_years_ends=humidity_at_the_years_ends(fine, coarse, next_second)),
        wind_at_the_level=level_mean,
        cell_raises=dict(factor=dict(water=mean(wind['factor'], sea2), land=mean(wind['factor'], ~sea2)),
                         excess_m_s=dict(water=mean(wind['excess'], sea2), land=mean(wind['excess'], ~sea2))),
        passes=differences(stats, names, area, sea2, lat),
        levels=levels(stats, names, area, sea2, lat),
        sites=sites,
        inputs={str(p.relative_to(F.RUNS)): F.sha256(p) for p in inputs}
        | {str(WM.FORCING.relative_to(ROOT)): F.sha256(WM.FORCING)},
        seconds=round(time.time() - started, 1))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    for name, rows in record['passes'].items():
        print(name, {case: {key: (rows[case][key]['water_cells']['mean'], rows[case][key]['land_cells']['mean']) for key in KEYS}
                     for case in ('still_2m', 'mixed_10m')})
    print(OUT, record['seconds'], 's')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
