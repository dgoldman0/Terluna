"""The surface wind of an ExoPlaSim run, for the state of its seas, written as a data product for other lanes.

    climate/gcm/.venv/bin/python -m climate.gcm.surface_wind A:35-49     # products/surface_wind_A_y35-49.npz

A water surface reflects the Sun and the sky through the slopes of its waves, and its short waves grow with the
friction velocity of the wind's stress on it (illumination/water_surface/short_waves.py). For the chosen model years this
module reads the run's 3-day output means of the surface stress (tauu, tauv), the 10 m wind (uas, vas), the surface
pressure, the 2 m air temperature and the Sun's height, and gives for each cell of the model grid the friction velocity
u* = sqrt(|tau| / rho), with rho = ps / (R T) the air's density at the surface from the run's own gas constant (its
planet_namelist's GASCON), and the 10 m wind speed: their 10th, 50th and 90th percentiles over the output intervals in
which the Sun is up there (the interval's mean cosine of the Sun's zenith angle above DAY_COSINE), the median over all
intervals, the number of daylit intervals, and the land mask and sea ice. Each output is a 3-day mean, a tenth of a lunar
day: the size of a 3-day mean stress or wind falls short of the mean size where the wind turns within the interval, and
the spread of 3-day values is smaller than that of the instantaneous wind.
"""
from __future__ import annotations
import hashlib
import json
import re
import sys
import warnings
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.climate.gcm-surface-wind/1'
DAY_COSINE = 0.1          # an interval counts as daylit where its mean cosine of the Sun's zenith angle exceeds this
PERCENTILES = (10, 50, 90)
EVIDENCE = ('Percentiles of ExoPlaSim 3-day output means over the stated model years at T21 (32 x 64), each cell a 170 '
            'km mean. The stress comes from PlaSim\'s bulk surface transfer with the Terluna builds\' Beljaars gust in '
            'unstable air (README, "Beljaars\'s gust in PlaSim"); gusts below the model\'s scales, sea breezes and the '
            'turning of the wind within an interval are absent from it.')
READING_RULE = ('Latitudes north first as in the model; longitudes east, 0-360, in the run\'s own frame. ustar_day_ms[k] '
                'and wind10_day_ms[k] are the PERCENTILES[k] percentile of the friction velocity (m/s) and of the 10 m '
                'wind speed (m/s) over the daylit intervals; ustar_all_ms and wind10_all_ms the medians over all '
                'intervals; day_intervals the count of daylit intervals; lsm and sic the run\'s mean land mask and sea '
                'ice. Over land the stress is that of the land\'s own roughness.')


def gas_constant(model_dir: Path) -> float:
    """The run's specific gas constant (J/kg/K), GASCON in its planet_namelist."""
    found = re.search(r'^\s*GASCON\s*=\s*([0-9.eE+-]+)', (Path(model_dir) / 'planet_namelist').read_text(), re.M)
    if not found:
        raise ValueError(f'no GASCON in {model_dir}/planet_namelist')
    return float(found.group(1))


def friction_velocity(tauu, tauv, ps_hpa, tas_k, gascon):
    """u* = sqrt(|tau| / rho), rho = ps / (R T), from stress components (Pa), pressure (hPa) and temperature (K)."""
    rho = np.asarray(ps_hpa, float) * 100.0 / (gascon * np.asarray(tas_k, float))
    return np.sqrt(np.hypot(tauu, tauv) / rho)


def daylit_percentiles(values, czen, percentiles=PERCENTILES, threshold=DAY_COSINE):
    """Percentiles over time (axis 0) of values where the interval's mean Sun cosine exceeds the threshold, NaN where a
    cell has no daylit interval, and the count of daylit intervals."""
    day = np.asarray(czen) > threshold
    masked = np.where(day, values, np.nan)
    count = day.sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)          # cells with no daylit interval give NaN
        out = np.nanpercentile(masked, percentiles, axis=0)
    return out, count


def read(folder: str, first: int, last: int):
    import netCDF4
    from climate.gcm.exoplasim_run import RUNS
    model = RUNS / folder / 'model'
    gascon = gas_constant(model)
    keys = ('tauu', 'tauv', 'uas', 'vas', 'ps', 'tas', 'czen', 'lsm', 'sic')
    fields = {k: [] for k in keys}
    for year in range(first, last + 1):
        with netCDF4.Dataset(model / f'MOST.{year:05d}.nc') as d:
            for k in keys:
                fields[k].append(np.asarray(d[k][:], dtype=float))
            lat = np.asarray(d['lat'][:], dtype=float)
            lon = np.asarray(d['lon'][:], dtype=float)
    return {k: np.concatenate(v) for k, v in fields.items()}, lat, lon, gascon


def surface_wind(folder: str, first: int, last: int) -> dict:
    f, lat, lon, gascon = read(folder, first, last)
    ustar = friction_velocity(f['tauu'], f['tauv'], f['ps'], f['tas'], gascon)
    wind = np.hypot(f['uas'], f['vas'])
    ustar_day, count = daylit_percentiles(ustar, f['czen'])
    wind_day, _ = daylit_percentiles(wind, f['czen'])
    return dict(lat=lat, lon=lon, percentiles=np.array(PERCENTILES, float), ustar_day_ms=ustar_day,
                wind10_day_ms=wind_day, ustar_all_ms=np.median(ustar, axis=0), wind10_all_ms=np.median(wind, axis=0),
                day_intervals=count.astype(float), lsm=f['lsm'].mean(axis=0), sic=f['sic'].mean(axis=0)), gascon


def main(argv=None) -> int:
    args = list(argv or sys.argv[1:])
    if len(args) != 1 or ':' not in args[0]:
        raise SystemExit('Give one run and year range, for example A:35-49')
    folder, span = args[0].split(':')
    first, last = (int(v) for v in span.split('-'))
    data, gascon = surface_wind(folder, first, last)
    from climate.gcm.exoplasim_run import RUNS
    progress = json.loads((RUNS / folder / 'progress.json').read_text())
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    meta = dict(schema=SCHEMA, run=folder, years=[first, last], configuration=progress['configuration'],
                gas_constant_J_kg_K=gascon, day_cosine=DAY_COSINE,
                producer=dict(domain='climate', files={'gcm/surface_wind.py': digest(__file__)}),
                evidence=EVIDENCE, reading_rule=READING_RULE)
    out = HERE / 'products' / f'surface_wind_{folder}_y{first}-{last}.npz'
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, metadata=json.dumps(meta), **{k: v.astype(np.float32) for k, v in data.items()})
    sea = (data['lsm'] < 0.5) & (data['day_intervals'] > 0)
    w = np.repeat(np.cos(np.radians(data['lat']))[:, None], data['lon'].size, axis=1)[sea]
    mean = lambda a: float((a[sea] * w).sum() / w.sum())
    print(f"{out}: over the open-water cells where the Sun is up, median friction velocity {mean(data['ustar_day_ms'][1]):.3f} "
          f"m/s (10th percentile {mean(data['ustar_day_ms'][0]):.3f}, 90th {mean(data['ustar_day_ms'][2]):.3f}), median 10 m "
          f"wind {mean(data['wind10_day_ms'][1]):.2f} m/s; all hours {mean(data['ustar_all_ms']):.3f} m/s")
    return 0


if __name__ == '__main__':
    sys.exit(main())
