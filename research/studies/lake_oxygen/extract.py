"""What the study reads from the GCM's own output: the forcing over each chosen lake and the waters' surface.

    climate/gcm/.venv/bin/python -m research.studies.lake_oxygen.extract     # about a minute; needs netCDF4

For the design case (A28_corrected_970, the lock) and the 20-year test with the longwave refitted by height and ice
optics for the cold cloud (A28_corrected_970_lw_ice, the hot world), over years 10-19 of each, in their 3-day means
(10 to a lunar day):
- at the T21 cell under each lake of results/lakes.json, and at the nearest water cell in the same row: PlaSim's
  tas, which is not air at 2 m but the lowest level's temperature brought down its dry adiabat to the surface pressure,
  that level's potential temperature at the surface (fluxmod.f90's dtsa); the lowest level's air, humidity and wind
  (about 0.9 km up), the surface pressure, the surface's net shortwave and longwave, its albedo and temperature, the
  sensible and latent heat, the mean cosine of the Sun's zenith angle and PlaSim's land mask;
- at every water cell (land fraction under one half in the T21 boundary file, as climate_targets counts the seas):
  the surface temperature and pressure, with each cell's sea area and latitude;
- the subsolar longitude of each mean, from the Sun's height near the equator.
It reads the same fields at the lake cells from the two years of the design case written every 6.4 hours
(A28_corrected_970_eddies), which measure what the 3-day means smooth away. It writes
research/runs/lake_oxygen/forcing.npz and forcing.json (ignored, on the research drive), with the sha256 of every
file read; run.py binds them.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS_DIR = ROOT / 'climate' / 'gcm' / 'runs'
BOUNDARY = ROOT / 'climate' / 'gcm' / 'products' / 'moon_28pct_water_lakes_gaussian_T21.nc'
LAKES = HERE / 'results' / 'lakes.json'
OUT = ROOT / 'research' / 'runs' / 'lake_oxygen' / 'forcing.npz'
META = OUT.with_suffix('.json')
SCHEMA = 'terluna.research.lake-oxygen-forcing/1'
DESIGN, TEST, EDDIES = 'A28_corrected_970', 'A28_corrected_970_lw_ice', 'A28_corrected_970_eddies'
YEARS = (10, 19)
EDDY_YEARS = (0, 1)
SURFACE = ('tas', 'ps', 'rss', 'rls', 'alb', 'ts', 'hfss', 'hfls', 'czen', 'lsm', 'clt', 'pr')
LOWEST = ('ta', 'hus', 'ua', 'va')


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
                    elevation=np.asarray(d['land_elevation_m'][:], float),
                    water_surface=np.asarray(d['water_surface_m'][:], float),
                    area=np.asarray(d['cell_area_fraction'][:], float))


def subsolar(czen, lat, lon):
    """Longitude of greatest mean Sun height near the equator for each mean, as climate/gcm/climatology does."""
    rows = np.abs(lat) < 20.0
    return lon[np.argmax(czen[:, rows, :].mean(axis=1), axis=1)]


def cells(lat, lon, lsm):
    """Each chosen lake's T21 cell and the nearest water cell (PlaSim's mask) in the same row."""
    lakes = json.loads(LAKES.read_text())['chosen']
    out = []
    for k in lakes:
        j = int(np.argmin(np.abs(lat - k['lat'])))
        i = int(np.argmin(np.abs(((lon - k['lon_east'] + 180.0) % 360.0) - 180.0)))
        water = np.flatnonzero(lsm[j] < 0.5)
        if water.size == 0:
            raise SystemExit(f'no water cell in row {j}')
        gap = np.abs(((lon[water] - lon[i] + 180.0) % 360.0) - 180.0)
        out.append(dict(band=k['band'], row=j, column=int(i), water_column=int(water[np.argmin(gap)]),
                        lat=float(lat[j]), lon=float(lon[i]), water_lon=float(lon[water[np.argmin(gap)]])))
    return out


def read_run(run, years, sites, b, keep_waters=True):
    """The fields at the sites (lake and water cells) and the waters' surface over the years."""
    import netCDF4
    files, series = {}, {}
    lat = lon = sigma = None
    sub, czen_all, waters_ts, waters_ps = [], [], [], []
    for year in range(years[0], years[1] + 1):
        path = output(run, year)
        files[str(path.relative_to(ROOT))] = sha256(path)
        with netCDF4.Dataset(path) as d:
            lat = np.asarray(d['lat'][:], float)
            lon = np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            czen = np.asarray(d['czen'][:], float)
            sub.append(subsolar(czen, lat, lon))
            rows = [s['row'] for s in sites] + [s['row'] for s in sites]
            columns = [s['column'] for s in sites] + [s['water_column'] for s in sites]
            for k in SURFACE:
                v = np.asarray(d[k][:], float)
                series.setdefault(k, []).append(v[:, rows, columns])
            for k in LOWEST:
                v = np.asarray(d[k][:, -1, :, :], float)
                series.setdefault(k, []).append(v[:, rows, columns])
            if keep_waters:
                waters_ts.append(np.asarray(d['ts'][:], float))
                waters_ps.append(np.asarray(d['ps'][:], float))
    out = {k: np.concatenate(v) for k, v in series.items()}
    out['subsolar_lon'] = np.concatenate(sub)
    if keep_waters:
        order = [int(np.argmin(np.abs(b['lat'] - x))) for x in lat]
        land, area = b['land'][order], b['area'][order]
        sea = land < 0.5
        out['waters_ts'] = np.concatenate(waters_ts)[:, sea]
        out['waters_ps'] = np.concatenate(waters_ps)[:, sea]
        out['waters_weight'] = ((1.0 - land) * area)[sea]
        out['waters_lat'] = (lat[:, None] * np.ones_like(land))[sea]
        out['waters_surface_m'] = b['water_surface'][order][sea]
    return out, files, lat, lon, sigma


def main(argv=None) -> int:
    import netCDF4
    b = boundary()
    with netCDF4.Dataset(output(DESIGN, YEARS[0])) as d:
        lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
        lsm = np.asarray(d['lsm'][0], float)
    sites = cells(lat, lon, lsm)
    order = [int(np.argmin(np.abs(b['lat'] - x))) for x in lat]
    for s in sites:
        s['land_fraction'] = float(b['land'][order][s['row'], s['column']])
        s['elevation_m'] = float(b['elevation'][order][s['row'], s['column']])
        s['water_elevation_m'] = float(b['water_surface'][order][s['row'], s['water_column']])
    files = {str(BOUNDARY.relative_to(ROOT)): sha256(BOUNDARY), str(LAKES.relative_to(ROOT)): sha256(LAKES)}
    arrays, runs = {}, {}
    for name, run, years, waters in (('lock', DESIGN, YEARS, True), ('hot', TEST, YEARS, True),
                                     ('eddies', EDDIES, EDDY_YEARS, False)):
        data, f, _, _, sigma = read_run(run, years, sites, b, keep_waters=waters)
        files.update(f)
        for k, v in data.items():
            arrays[f'{name}__{k}'] = v.astype(np.float64)
        with netCDF4.Dataset(output(run, years[0])) as d:
            steps = int(np.diff(np.asarray(d['time'][:2], float))[0]) if d['time'].size > 1 else None
        progress = RUNS_DIR / run / 'progress.json'
        p = json.loads(progress.read_text())
        files[str(progress.relative_to(ROOT))] = sha256(progress)
        m = p['configuration']['model']
        runs[name] = dict(run=run, years=list(years), means=int(data['tas'].shape[0]), timestep_min=m['timestep_min'],
                          steps_per_lunar_day=m['steps_per_lunar_day'], writes_per_lunar_day=m['writes_per_lunar_day'],
                          steps_per_mean=m['steps_per_lunar_day'] // m['writes_per_lunar_day'], time_axis_steps=steps,
                          sunlight_scale=m['sunlight_scale'], gravity=p['configuration']['planet']['gravity'],
                          lowest_sigma=float(sigma[-1]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT, **arrays)
    meta = dict(schema=SCHEMA, sites=sites, runs=runs, lat=lat.tolist(), lon=lon.tolist(),
                columns=dict(site_order='the three lake cells, then the three water cells, in results/lakes.json order',
                             surface=list(SURFACE), lowest_level=list(LOWEST)),
                producer=dict(files={str(Path(__file__).resolve().relative_to(ROOT)): sha256(__file__)}, read=files),
                npz_sha256=sha256(OUT))
    META.write_text(json.dumps(meta, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT), 'and', META.relative_to(ROOT))
    for s in sites:
        print(s)
    return 0


if __name__ == '__main__':
    sys.exit(main())
