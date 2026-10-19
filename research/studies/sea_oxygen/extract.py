"""What the sea study reads from the GCM's own output: the weather over the seas' run basins.

    /home/kir/Documents/Projects/terluna/climate/gcm/.venv/bin/python -m research.studies.sea_oxygen.extract
    # about two minutes; needs netCDF4, so it runs in the GCM's environment

For the design case (A28_corrected_970, the lock) and the 20-year test with the longwave refitted by height and ice
optics (A28_corrected_970_lw_ice, the hot world), over years 10-19 of each, in their 3-day means, as the waters'
oxygen study's extract reads them (research/studies/lake_oxygen/extract.py, whose readers this imports):
- at every T21 cell holding water of a run basin in results/seas.json that PlaSim takes as water at the seas' level
  (its land mask is water and the boundary file's water surface lies within 50 m of the common level): the air at
  2 m, the lowest level's air, humidity and wind (about 0.9 km up), the surface pressure, the surface's net
  shortwave and longwave, its albedo and temperature, the sensible and latent heat, the mean cosine of the Sun's
  zenith angle, PlaSim's land mask, the cloud cover and the rain; a basin with no such cell takes the cell holding
  most of its water, flagged, so that run.py carries its air down to the seas' level as lake_oxygen carries its lakes';
- at every water cell (land fraction under one half in the boundary file, as climate_targets counts the seas): the
  surface temperature and pressure with each cell's sea area, latitude and water surface, as lake_oxygen's extract
  keeps them, for the reproduction of its figures;
- the subsolar longitude of each mean, from the Sun's height near the equator.
It reads the same fields at the same cells from the two years of the design case written every 6.4 hours
(A28_corrected_970_eddies), which measure what the 3-day means smooth away. The values are stored as the GCM wrote
them (float32). It writes research/runs/sea_oxygen/forcing.npz and forcing.json (ignored, on the research drive),
with the sha256 of every file read; run.py binds them.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from research.studies.lake_oxygen import extract as X

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SEAS = HERE / 'results' / 'seas.json'
OUT = ROOT / 'research' / 'runs' / 'sea_oxygen' / 'forcing.npz'
META = OUT.with_suffix('.json')
SCHEMA = 'terluna.research.sea-oxygen-forcing/1'
SEA_LEVEL_WITHIN_M = 50.0             # the boundary's water surface (m above the common level) of a sea cell


def sites_of(seas, lsm, b):
    """The cells the study reads: each run basin's water cells at the seas' level, or its largest cell."""
    cells = {}
    for basin in seas['run_basins']:
        water = [c for c in basin['t21_cells']
                 if lsm[c['row'], c['column']] < 0.5 and b['water_surface'][c['row'], c['column']] < SEA_LEVEL_WITHIN_M]
        take = water or [max(basin['t21_cells'], key=lambda c: c['area_km2'])]
        for c in take:
            s = cells.setdefault((c['row'], c['column']), dict(row=c['row'], column=c['column'], basins={}))
            s['basins'][basin['key']] = c['area_km2']
            s['land_cell'] = not water
    return [cells[k] for k in sorted(cells)]


def read_run(run, years, sites, b, keep_waters=True):
    """The fields at the sites and the waters' surface over the years, as the GCM wrote them."""
    import netCDF4
    files, series = {}, {}
    rows = np.array([s['row'] for s in sites])
    cols = np.array([s['column'] for s in sites])
    sub, waters_ts, waters_ps = [], [], []
    lat = lon = sigma = None
    for year in range(years[0], years[1] + 1):
        path = X.output(run, year)
        files[str(path.relative_to(ROOT))] = X.sha256(path)
        with netCDF4.Dataset(path) as d:
            lat = np.asarray(d['lat'][:], float)
            lon = np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            czen = np.asarray(d['czen'][:])
            sub.append(X.subsolar(czen.astype(float), lat, lon))
            for k in X.SURFACE:
                v = np.asarray(d[k][:])
                series.setdefault(k, []).append(v[:, rows, cols])
            for k in X.LOWEST:
                v = np.asarray(d[k][:, -1, :, :])
                series.setdefault(k, []).append(v[:, rows, cols])
            if keep_waters:
                waters_ts.append(np.asarray(d['ts'][:]))
                waters_ps.append(np.asarray(d['ps'][:]))
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
        out['waters_row'] = np.broadcast_to(np.arange(lat.size)[:, None], land.shape)[sea]
        out['waters_column'] = np.broadcast_to(np.arange(lon.size)[None, :], land.shape)[sea]
    return out, files, lat, lon, sigma


def main(argv=None) -> int:
    import netCDF4
    seas = json.loads(SEAS.read_text())
    b = X.boundary()
    with netCDF4.Dataset(X.output(X.DESIGN, X.YEARS[0])) as d:
        lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
        lsm = np.asarray(d['lsm'][0], float)
    if np.abs(lat - np.asarray(seas['t21']['lat'])).max() > 1e-3 or np.abs(lon - seas['t21']['lon']).max() > 1e-6:
        raise SystemExit('the GCM grid differs from the T21 grid seas.json assigns the basins to')
    order = [int(np.argmin(np.abs(b['lat'] - x))) for x in lat]
    bb = dict(b, water_surface=b['water_surface'][order], land=b['land'][order], elevation=b['elevation'][order])
    sites = sites_of(seas, lsm, bb)
    for i, s in enumerate(sites):
        s.update(index=i, lat=float(lat[s['row']]), lon=float(lon[s['column']]), lsm=float(lsm[s['row'], s['column']]),
                 land_fraction=float(bb['land'][s['row'], s['column']]),
                 elevation_m=float(bb['elevation'][s['row'], s['column']]),
                 water_surface_m=float(bb['water_surface'][s['row'], s['column']]))
    files = {str(X.BOUNDARY.relative_to(ROOT)): X.sha256(X.BOUNDARY), str(SEAS.relative_to(ROOT)): X.sha256(SEAS)}
    arrays, runs = {}, {}
    for name, run, years, waters in (('lock', X.DESIGN, X.YEARS, True), ('hot', X.TEST, X.YEARS, True),
                                     ('eddies', X.EDDIES, X.EDDY_YEARS, False)):
        data, f, _, _, sigma = read_run(run, years, sites, b, keep_waters=waters)
        files.update(f)
        for k, v in data.items():
            arrays[f'{name}__{k}'] = v
        with netCDF4.Dataset(X.output(run, years[0])) as d:
            steps = int(np.diff(np.asarray(d['time'][:2], float))[0]) if d['time'].size > 1 else None
        progress = X.RUNS_DIR / run / 'progress.json'
        p = json.loads(progress.read_text())
        files[str(progress.relative_to(ROOT))] = X.sha256(progress)
        m = p['configuration']['model']
        runs[name] = dict(run=run, years=list(years), means=int(data['tas'].shape[0]), timestep_min=m['timestep_min'],
                          steps_per_lunar_day=m['steps_per_lunar_day'], writes_per_lunar_day=m['writes_per_lunar_day'],
                          steps_per_mean=m['steps_per_lunar_day'] // m['writes_per_lunar_day'], time_axis_steps=steps,
                          sunlight_scale=m['sunlight_scale'], gravity=p['configuration']['planet']['gravity'],
                          lowest_sigma=float(sigma[-1]))
        print(f'  read {name}: {runs[name]["means"]} means at {len(sites)} cells', flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT, **arrays)
    meta = dict(schema=SCHEMA, sites=sites, runs=runs, lat=lat.tolist(), lon=lon.tolist(),
                columns=dict(site_order='the cells in sites, by row and column', surface=list(X.SURFACE),
                             lowest_level=list(X.LOWEST), stored='float32, as the GCM wrote them'),
                producer=dict(files={str(Path(p).resolve().relative_to(ROOT)): X.sha256(p)
                                     for p in (__file__, X.__file__)}, read=files),
                npz_sha256=X.sha256(OUT))
    META.write_text(json.dumps(meta, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT), 'and', META.relative_to(ROOT), f'({len(sites)} cells)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
