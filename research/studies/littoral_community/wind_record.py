"""The wind inside a 3-day mean, over the whole grid of the one record that shows it.

The corrected design climate is written as 3-day means, and the speed of a mean wind vector lies below the mean speed
of the wind that made it. The waters' oxygen study measured the difference at its six lake and water cells in the two
years of the earlier design case written every 6.4 hours (A28_corrected_970_eddies, years 0-1, on today's ground and
ten levels), and raises its wind by the cube's ratio there; the littoral forcing and the water map take the mean of
those six factors (forcing.wind_raise). This runner reads the same record at every cell, the lowest level's wind as
that study read it, and writes what the raise rests on:

- the ratio of the 6.4-hour means' mean speed to the speed of their 3-day mean vectors, and the same in the cube, at
  every cell, summarised over the grid, by surface and by latitude band;
- the same ratio by the speed of the window's mean vector, with the root of the mean square speed the windows hold
  beyond their mean vectors: how the raise would differ between calm windows and windy ones.

The record's wind is the lowest level's, about 0.9 km up, on a run whose ground and physics differ from the corrected
design climate's; the 10 m wind inside the corrected run's windows is not measured here.

    /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.littoral_community.wind_record  # a minute

writes results/wind_record.json (schema terluna.research.littoral-wind-record/1) with the sha256 of every input.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology.moisture_map import cell_areas_km2
from research.studies.littoral_community import forcing as F
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'wind_record.json'
SCHEMA = 'terluna.research.littoral-wind-record/1'
RECORD, RECORD_YEARS = 'A28_corrected_970_eddies', (0, 1)
WRITES_PER_WINDOW = 11                       # 6.4-hour writes to a 3-day mean: 110 and 10 to a lunar day
BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))
SPEED_EDGES = (0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0)


def record_files():
    return [F.RUNS / RECORD / 'model' / f'MOST.{y:05d}.nc' for y in range(RECORD_YEARS[0], RECORD_YEARS[1] + 1)]


def read_record():
    """The lowest level's wind at every write (float32), with the grid and the run's land mask."""
    import netCDF4
    u, v = [], []
    for path in record_files():
        with netCDF4.Dataset(path) as nc:
            lat, lon = np.array(nc['lat'][:], dtype=float), np.array(nc['lon'][:], dtype=float)
            lsm = np.array(nc['lsm'][0], dtype=float)
            sigma = float(np.array(nc['lev'][:], dtype=float)[-1])
            u.append(np.array(nc['ua'][:, -1], dtype=np.float32))
            v.append(np.array(nc['va'][:, -1], dtype=np.float32))
    return dict(u=np.concatenate(u), v=np.concatenate(v), lat=lat, lon=lon, lsm=lsm, sigma=sigma)


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


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    started = time.time()
    rec = read_record()
    lat, lon = rec['lat'], rec['lon']
    speed, vector = windows(rec['u'], rec['v'])
    linear, cube = ratios(speed, vector)
    area = cell_areas_km2(lat, lon)
    sea = rec['lsm'] <= 0.5
    lat2 = np.abs(lat)[:, None] * np.ones(area.shape)
    mean = lambda a, m: round(float((a * area)[m].sum() / area[m].sum()), 3)
    spread = lambda a: dict(lowest=round(float(a.min()), 3), p05=round(float(np.percentile(a, 5)), 3),
                            median=round(float(np.median(a)), 3), p95=round(float(np.percentile(a, 95)), 3),
                            highest=round(float(a.max()), 3))
    everywhere = np.ones(area.shape, bool)
    grid = dict(mean_speed_over_vector=dict(area_mean=mean(linear, everywhere), water=mean(linear, sea), land=mean(linear, ~sea),
                                            cells=spread(linear)),
                cube_ratio=dict(area_mean=mean(cube, everywhere), water=mean(cube, sea), land=mean(cube, ~sea),
                                cells=spread(cube)),
                speed_m_s=dict(water=mean(speed.mean((0, 1)), sea), land=mean(speed.mean((0, 1)), ~sea)),
                vector_m_s=dict(water=mean(vector.mean(0), sea), land=mean(vector.mean(0), ~sea)),
                excess_rms_m_s=dict(all=round(excess(speed, vector), 3),
                                    water=round(excess(speed, vector, np.broadcast_to(sea, vector.shape)), 3),
                                    land=round(excess(speed, vector, np.broadcast_to(~sea, vector.shape)), 3)))
    bands = []
    for lo, hi in BANDS:
        band = (lat2 >= lo) & (lat2 < hi)
        row = dict(band_deg=[lo, hi])
        for label, m in (('water', band & sea), ('land', band & ~sea)):
            row[label] = dict(mean_speed_over_vector=mean(linear, m), cube_ratio=mean(cube, m),
                              excess_rms_m_s=round(excess(speed, vector, np.broadcast_to(m, vector.shape)), 3))
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
    files_used = ['research/studies/littoral_community/wind_record.py', 'research/studies/littoral_community/forcing.py',
                  'biosphere/ecology/moisture_map.py']
    record = dict(
        schema=SCHEMA,
        producer=dict(lane='research', runner='research/studies/littoral_community/wind_record.py',
                      files={f: F.sha256(F.ROOT / f) for f in files_used}, constants=constants_used(files_used)),
        evidence=('Read from one climate-model record, the two years of the earlier design case written every 6.4 '
                  'hours: the lowest level\'s wind, about 0.9 km up, on today\'s ground and ten levels. It shows what '
                  '3-day means smooth away in that run; the 10 m wind inside the corrected design climate\'s windows '
                  'is not measured.'),
        reading_rule=('mean_speed_over_vector is the mean speed of the 6.4-hour means over the mean speed of their '
                      '3-day mean vectors; cube_ratio is the waters\' oxygen study\'s factor, the cube root of the mean '
                      'cubed speed over the mean cubed vector speed. excess_rms_m_s is the root of the mean square '
                      'speed the windows hold beyond their mean vectors. by_vector_speed groups windows of all cells '
                      'by their mean vector\'s speed. Band and surface means weight cells by area.'),
        settings=dict(run=RECORD, years=list(RECORD_YEARS), lowest_sigma=rec['sigma'], writes_per_window=WRITES_PER_WINDOW,
                      write_h=round(F.SYNODIC_H / 110.0, 3), windows=int(vector.shape[0]), bands=[list(b) for b in BANDS]),
        grid=grid, by_latitude=bands, by_vector_speed=by_speed,
        the_forcings_raise=dict(factor=F.wind_raise(), cells=F.wind_factors()),
        inputs={str(p.relative_to(F.RUNS)): F.sha256(p) for p in record_files()}
        | {str(F.LAKE_OXYGEN.relative_to(F.ROOT)): F.sha256(F.LAKE_OXYGEN)},
        seconds=round(time.time() - started, 1))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(json.dumps(dict(grid=grid, by_latitude=bands, by_vector_speed=by_speed), indent=1))
    print(OUT, record['seconds'], 's')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
