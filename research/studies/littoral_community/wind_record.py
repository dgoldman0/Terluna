"""The wind inside a 3-day mean, over the whole grid of the two records that show it.

The design years are kept as 3-day means, and the speed of a mean wind vector lies below the mean speed of the wind
that made it. Two records written every 13 half-hour steps (6.5 hours, eleven to each 3-day mean) show the difference:

- corrected: the corrected design run's year 10 again, run from that run's own year-9 restart with its writing
  interval alone changed (A28R_301k_design_l15_snow_s938_turned_refit_precession_eddies, year 0; domain/climate). It
  ends on the design run's year-10 restart byte for byte, so its weather is the design year's own
  (inside_windows.py checks its records against that year's 3-day means). It carries the wind at the level the
  waters' exchange is paired at, 313 m, and the run's 10 m wind over each cell's own ground.
- earlier: the two years of the earlier design case (A28_corrected_970_eddies, years 0-1, on today's ground and ten
  levels, the lowest about 0.9 km up), on which the waters' oxygen study measured the factor the littoral forcing and
  the water map raise their wind by (forcing.wind_raise: the mean of its six cells' cube ratios).

For each record and wind this runner writes what a raise rests on:

- the ratio of the 6.5-hour means' mean speed to the speed of their 3-day mean vectors, and the same in the cube, at
  every cell, summarised over the grid, by surface and by latitude band;
- the same ratio by the speed of the window's mean vector, with the root of the mean square speed the windows hold
  beyond their mean vectors: how a raise differs between calm windows and windy ones;
- for the lowest level, the forcing's own operation: the mean of the 6.5-hour speeds floored at PlaSim's least wind
  over the mean of the floored 3-day vector speeds (the factor that would stand in the borrowed one's place), and
  the excess that, added to the floored vector speed in quadrature, gives each cell's windows their mean speed.

The 6.5-hour means are themselves means: what the wind does inside 6.5 hours is not in either record.

    /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.littoral_community.wind_record  # a minute

writes results/wind_record.json (schema terluna.research.littoral-wind-record/2) with the sha256 of every input.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from biosphere.ecology.moisture_map import cell_areas_km2
from research.studies.littoral_community import forcing as F
from shared.constants import MOON_SURFACE_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'wind_record.json'
SCHEMA = 'terluna.research.littoral-wind-record/2'
CORRECTED = F.RUN + '_eddies'                # year 0 is the design run's year 10
RECORDS = dict(
    corrected=dict(run=CORRECTED, years=(0, 0),
                   what="the corrected design run's year 10, run again from its own year-9 restart and written "
                        'every 13 steps'),
    earlier=dict(run='A28_corrected_970_eddies', years=(0, 1),
                 what="two years of the earlier design case, on today's ground and ten levels, written every 13 "
                      'steps'))
WRITES_PER_WINDOW = 11                       # 6.5-hour writes to a 3-day mean: 110 and 10 to a model day
WRITE_STEPS = F.WINDOW_STEPS // WRITES_PER_WINDOW        # 13 half-hour steps
BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))
SPEED_EDGES = (0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0)


def record_files(key='earlier'):
    r = RECORDS[key]
    return [F.RUNS / r['run'] / 'model' / f'MOST.{y:05d}.nc' for y in range(r['years'][0], r['years'][1] + 1)]


def read_record(key='earlier'):
    """A record's wind at every write (float32): the lowest level's, and the run's 10 m wind where it writes one; with
    the grid, the run's land mask, the lowest level's sigma and its mean height over the grid's cells."""
    import netCDF4
    u, v, u10, v10, height = [], [], [], [], []
    for path in record_files(key):
        with netCDF4.Dataset(path) as nc:
            lat, lon = np.array(nc['lat'][:], dtype=float), np.array(nc['lon'][:], dtype=float)
            lsm = np.array(nc['lsm'][0], dtype=float)
            sigma = float(np.array(nc['lev'][:], dtype=float)[-1])
            u.append(np.array(nc['ua'][:, -1], dtype=np.float32))
            v.append(np.array(nc['va'][:, -1], dtype=np.float32))
            if 'uas' in nc.variables and 'vas' in nc.variables:
                u10.append(np.array(nc['uas'][:], dtype=np.float32))
                v10.append(np.array(nc['vas'][:], dtype=np.float32))
            mean_t = 0.5 * (np.array(nc['ta'][:, -1], dtype=float).mean(0) + np.array(nc['tas'][:], dtype=float).mean(0))
            height.append(F.air_constants()[0] * mean_t / MOON_SURFACE_GRAVITY * math.log(1.0 / sigma))
    out = dict(u=np.concatenate(u), v=np.concatenate(v), lat=lat, lon=lon, lsm=lsm, sigma=sigma,
               height_m=np.mean(height, 0))
    if u10:
        out.update(u10=np.concatenate(u10), v10=np.concatenate(v10))
    return out


def windows(u, v, per=WRITES_PER_WINDOW):
    """Each window's writes' speeds (windows x writes x cells...) and the speed of its mean vector."""
    n = (u.shape[0] // per) * per
    uu, vv = u[:n].reshape(-1, per, *u.shape[1:]).astype(float), v[:n].reshape(-1, per, *v.shape[1:]).astype(float)
    return np.hypot(uu, vv), np.hypot(uu.mean(1), vv.mean(1))


def ratios(speed, vector):
    """Per cell: the mean speed over the mean vector's speed, and the cube's ratio as the waters' oxygen study takes
    it, ((speed^3) / (vector^3))^(1/3) of the means."""
    linear = speed.mean((0, 1)) / vector.mean(0)
    cube = ((speed ** 3).mean((0, 1)) / (vector ** 3).mean(0)) ** (1.0 / 3.0)
    return linear, cube


def excess(speed, vector, mask=None):
    """The root of the mean square speed the windows hold beyond their mean vectors (m/s), over the windows and
    cells of a mask (windows x cells...)."""
    square = (speed ** 2).mean(1) - vector ** 2
    return float(np.sqrt(max(square.mean() if mask is None else square[mask].mean(), 0.0)))


def floored(speed, vector, floor=F.LEVEL_FLOOR_MS):
    """The forcing's own operation, per cell: each window's mean of the writes' speeds floored at PlaSim's least wind
    (windows x cells...), the floored speed of each window's mean vector, and the factor between their means, which a
    raise of the floored vector speed by one factor would need at that cell."""
    inside = np.maximum(speed, floor).mean(1)
    outside = np.maximum(vector, floor)
    return inside, outside, inside.mean(0) / outside.mean(0)


def quadrature_excess(inside, outside, highest=20.0):
    """Per cell, the excess e (m/s) that gives its windows their mean speed when added to the floored vector speed in
    quadrature: mean(sqrt(outside^2 + e^2)) = mean(inside), by bisection; zero where the windows hold none."""
    target = inside.mean(0)
    lo, hi = np.zeros(target.shape), np.full(target.shape, highest)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        low = np.sqrt(outside ** 2 + mid ** 2).mean(0) < target
        lo, hi = np.where(low, mid, lo), np.where(low, hi, mid)
    return np.where(target > outside.mean(0), 0.5 * (lo + hi), 0.0)


def describe(speed, vector, area, sea, lat):
    """A wind's statistics over the grid: the two ratios by cell, surface and latitude band, and by the speed of the
    window's mean vector, with the excess."""
    linear, cube = ratios(speed, vector)
    lat2 = np.abs(lat)[:, None] * np.ones(area.shape)
    mean = lambda a, m: round(float((a * area)[m].sum() / area[m].sum()), 3)
    spread = lambda a: dict(lowest=round(float(a.min()), 3), p05=round(float(np.percentile(a, 5)), 3),
                            median=round(float(np.median(a)), 3), p95=round(float(np.percentile(a, 95)), 3),
                            highest=round(float(a.max()), 3))
    over = lambda m: np.broadcast_to(m, vector.shape)
    everywhere = np.ones(area.shape, bool)
    grid = dict(mean_speed_over_vector=dict(area_mean=mean(linear, everywhere), water=mean(linear, sea), land=mean(linear, ~sea),
                                            cells=spread(linear)),
                cube_ratio=dict(area_mean=mean(cube, everywhere), water=mean(cube, sea), land=mean(cube, ~sea),
                                cells=spread(cube)),
                speed_m_s=dict(water=mean(speed.mean((0, 1)), sea), land=mean(speed.mean((0, 1)), ~sea)),
                vector_m_s=dict(water=mean(vector.mean(0), sea), land=mean(vector.mean(0), ~sea)),
                excess_rms_m_s=dict(all=round(excess(speed, vector), 3), water=round(excess(speed, vector, over(sea)), 3),
                                    land=round(excess(speed, vector, over(~sea)), 3)))
    bands = []
    for lo, hi in BANDS:
        band = (lat2 >= lo) & (lat2 < hi)
        row = dict(band_deg=[lo, hi])
        for label, m in (('water', band & sea), ('land', band & ~sea)):
            row[label] = dict(mean_speed_over_vector=mean(linear, m), cube_ratio=mean(cube, m),
                              excess_rms_m_s=round(excess(speed, vector, over(m)), 3))
        bands.append(row)
    by_speed = []
    window_speed = speed.mean(1)
    for lo, hi in zip(SPEED_EDGES, SPEED_EDGES[1:] + (np.inf,)):
        m = (vector >= lo) & (vector < hi)
        if not m.any():
            continue
        by_speed.append(dict(vector_m_s=[lo, None if np.isinf(hi) else hi], share_of_windows=round(float(m.mean()), 4),
                             mean_vector_m_s=round(float(vector[m].mean()), 3), mean_speed_m_s=round(float(window_speed[m].mean()), 3),
                             ratio=round(float(window_speed[m].mean() / vector[m].mean()), 3),
                             excess_rms_m_s=round(excess(speed, vector, m), 3)))
    return dict(grid=grid, by_latitude=bands, by_vector_speed=by_speed)


def describe_floored(speed, vector, area, sea, lat):
    """The forcing's operation at the lowest level over the grid: the factor and the quadrature excess by cell,
    surface and latitude band, and how often the wind stands under the floor."""
    inside, outside, factor = floored(speed, vector)
    e = quadrature_excess(inside, outside)
    lat2 = np.abs(lat)[:, None] * np.ones(area.shape)
    mean = lambda a, m: round(float((a * area)[m].sum() / area[m].sum()), 3)
    spread = lambda a, m: dict(p05=round(float(np.percentile(a[m], 5)), 3), median=round(float(np.median(a[m])), 3),
                               p95=round(float(np.percentile(a[m], 95)), 3))
    out = dict(floor_m_s=F.LEVEL_FLOOR_MS,
               factor=dict(water=mean(factor, sea), land=mean(factor, ~sea), water_cells=spread(factor, sea),
                           land_cells=spread(factor, ~sea)),
               quadrature_excess_m_s=dict(water=mean(e, sea), land=mean(e, ~sea), water_cells=spread(e, sea),
                                          land_cells=spread(e, ~sea)),
               writes_under_the_floor=dict(water=round(float((speed < F.LEVEL_FLOOR_MS)[:, :, sea].mean()), 4),
                                           land=round(float((speed < F.LEVEL_FLOOR_MS)[:, :, ~sea].mean()), 4)),
               windows_with_vector_under_the_floor=dict(water=round(float((vector < F.LEVEL_FLOOR_MS)[:, sea].mean()), 4),
                                                        land=round(float((vector < F.LEVEL_FLOOR_MS)[:, ~sea].mean()), 4)),
               by_latitude=[])
    for lo, hi in BANDS:
        band = (lat2 >= lo) & (lat2 < hi)
        out['by_latitude'].append(dict(band_deg=[lo, hi], water=dict(factor=mean(factor, band & sea), quadrature_excess_m_s=mean(e, band & sea)),
                                       land=dict(factor=mean(factor, band & ~sea), quadrature_excess_m_s=mean(e, band & ~sea))))
    return out


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    started = time.time()
    records, inputs = {}, {}
    for key, r in RECORDS.items():
        rec = read_record(key)
        lat, lon = rec['lat'], rec['lon']
        area = cell_areas_km2(lat, lon)
        sea = rec['lsm'] <= 0.5
        speed, vector = windows(rec['u'], rec['v'])
        row = dict(run=r['run'], years=list(r['years']), what=r['what'], lowest_sigma=rec['sigma'],
                   lowest_level_height_m=round(float((rec['height_m'] * area).sum() / area.sum()), 1),
                   windows=int(vector.shape[0]), water_share_of_area=round(float(area[sea].sum() / area.sum()), 4),
                   lowest_level=describe(speed, vector, area, sea, lat))
        row['lowest_level']['floored'] = describe_floored(speed, vector, area, sea, lat)
        if 'u10' in rec:
            row['at_10_m'] = describe(*windows(rec['u10'], rec['v10']), area, sea, lat)
        records[key] = row
        inputs.update({str(p.relative_to(F.RUNS)): F.sha256(p) for p in record_files(key)})
        print(key, json.dumps({k: row['lowest_level']['grid'][k] for k in ('mean_speed_over_vector', 'cube_ratio', 'excess_rms_m_s')}),
              flush=True)
    files_used = ['research/studies/littoral_community/wind_record.py', 'research/studies/littoral_community/forcing.py',
                  'biosphere/ecology/moisture_map.py']
    record = dict(
        schema=SCHEMA,
        producer=dict(lane='research', runner='research/studies/littoral_community/wind_record.py',
                      files={f: F.sha256(F.ROOT / f) for f in files_used}, constants=constants_used(files_used)),
        evidence=('Read from two climate-model records written every 13 half-hour steps (6.5 hours). corrected is the '
                  "corrected design run's year 10 again, from its own year-9 restart, ending on that run's year-10 "
                  "restart byte for byte: the weather inside one of the design years' own 3-day means, at the level "
                  "the waters' exchange is paired at and at 10 m over each cell's own ground. earlier is two years "
                  "of the earlier design case, on today's ground and ten levels, the record the forcing's raise was "
                  'measured on. One year and two years are short records; the 6.5-hour means are themselves means, '
                  'and what the wind does inside them is in neither.'),
        reading_rule=('mean_speed_over_vector is the mean speed of the 6.5-hour means over the mean speed of their '
                      "3-day mean vectors; cube_ratio is the waters' oxygen study's factor, the cube root of the mean "
                      'cubed speed over the mean cubed vector speed. excess_rms_m_s is the root of the mean square '
                      'speed the windows hold beyond their mean vectors. by_vector_speed groups windows of all cells '
                      "by their mean vector's speed. floored reads the forcing's operation at the lowest level: "
                      "factor is the mean of the 6.5-hour speeds floored at PlaSim's least wind over the mean of the "
                      'floored 3-day vector speeds, and quadrature_excess_m_s the excess that, added to the floored '
                      "vector speed in quadrature, gives a cell's windows their mean speed. Band and surface means "
                      'weight cells by area; the cells entries are percentiles over cells.'),
        settings=dict(writes_per_window=WRITES_PER_WINDOW, write_steps=WRITE_STEPS, write_h=WRITE_STEPS * 0.5,
                      window_h=F.WINDOW_H, bands=[list(b) for b in BANDS], speed_edges_m_s=list(SPEED_EDGES)),
        records=records,
        the_forcings_raise=dict(factor=F.wind_raise(), cells=F.wind_factors()),
        inputs=inputs | {str(F.LAKE_OXYGEN.relative_to(F.ROOT)): F.sha256(F.LAKE_OXYGEN)},
        seconds=round(time.time() - started, 1))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(OUT, record['seconds'], 's')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
