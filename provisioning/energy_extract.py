"""What the energy-options study reads from the design run: sunlight at the ground, cloud and the low winds.

    climate/gcm/.venv/bin/python -m provisioning.energy_extract
    # -> provisioning/results/energy_climate.json, in about a minute; needs netCDF4 and the run's output

The design climate is the resort at 301 K, and its design run is A28R_301k_design_l15_snow_s934 (ExoPlaSim at T21 with
15 levels, the sunlight scaled by 0.934), read over years 10-19 as the climate studies read it. The output is in
climate/gcm/runs (a link to the research drive, ignored by Git); every file read is recorded with its sha256 and size.

Each record is a mean over a tenth of the lunar day (143 steps of 30 minutes). From each the extract takes:
- sunlight arriving at the ground, the global horizontal irradiance a flat panel meets: the net surface shortwave
  (rss) over one minus the surface albedo (alb), all-sky and clear-sky (rsscs);
- sunlight arriving at the top of the air: the net (rst) less the outgoing (rsut), which the output signs negative;
- cloud cover (clt), the air at 2 m (tas), its water column (prw) and the winds of the lowest level (sigma 0.994,
  about 0.3 km up) with the surface stress.
Means are area-weighted by latitude row, separately for the model's land cells (lsm) and its water cells.

The lunar day: the Sun's place in each record is the longitude where the record's mean cosine of the Sun's zenith
angle (czen) peaks over the rows within 15 degrees of the equator, refined by a parabola through the peak and its two
neighbours. A cell's hour angle is its longitude less that place, and the land's sunlight is binned by hour angle in
10-degree bins over latitude bands of 15 degrees, north and south together. A record's mean spans 36 degrees of the
Sun's travel, so the profile is the day's shape smoothed over a tenth of the day.

The winds are 3-day means of the wind vector, so their speeds and cubes understate the instantaneous winds a turbine
meets; the gusts and storms within each mean are missing (reading rule).
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA = 'terluna.provisioning.energy-climate/1'
OUT = HERE / 'results' / 'energy_climate.json'
CODE = ('provisioning/energy_extract.py',)
RUN = 'A28R_301k_design_l15_snow_s934'
YEARS = tuple(range(10, 20))
RUN_DIR = ROOT / 'climate' / 'gcm' / 'runs' / RUN
HOUR_BIN_DEG = 10.0
BAND_DEG = 15.0
EQUATOR_ROWS_DEG = 15.0
CUT_IN_M_S = 3.0                   # a typical turbine's cut-in speed (published class value, for the share below it)

EVIDENCE = ("Means of the design run's own output over years 10-19, by latitude row and for the model's land and "
            "water cells: a GCM at T21 (about 170 km a cell) with its radiation validated against line-by-line columns "
            "and its clouds as the climate work leaves them. Sunlight at the ground is all-sky and clear-sky on a flat "
            "surface; the winds are 3-day means of the lowest level, which understate the instantaneous winds.")
READING_RULE = ("W/m2 are means over every hour of the lunar cycle, night included, unless labelled; by_row lists the "
                "model's 32 latitude rows from north to south; bands fold north and south and are named by their "
                "equatorward edge; the day profile's bins are hour angles in degrees from local noon, west of noon "
                "negative (morning); winds in m/s at sigma 0.994, about 0.3 km up; shares are of time.")


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def sun_longitude(czen, lat, lon):
    """Each record's subsolar longitude (degrees) from where czen peaks over the equatorial rows."""
    rows = np.abs(lat) < EQUATOR_ROWS_DEG
    m = czen[:, rows, :].mean(axis=1)                       # (time, lon)
    k = np.argmax(m, axis=1)
    n = lon.size
    y0, y1, y2 = (m[np.arange(m.shape[0]), (k + s) % n] for s in (-1, 0, 1))
    denom = y0 - 2 * y1 + y2
    shift = np.where(np.abs(denom) > 1e-12, 0.5 * (y0 - y2) / denom, 0.0)
    return (lon[k] + shift * (lon[1] - lon[0])) % 360.0


def read_year(path):
    import netCDF4
    with netCDF4.Dataset(path) as d:
        v = {k: np.asarray(d.variables[k][:], float) for k in
             ('lat', 'lon', 'rss', 'rsscs', 'alb', 'clt', 'czen', 'rst', 'rsut', 'tas', 'prw', 'tauu', 'tauv')}
        v['lsm'] = np.asarray(d.variables['lsm'][0], float)
        v['u'] = np.asarray(d.variables['ua'][:, -1], float)
        v['v'] = np.asarray(d.variables['va'][:, -1], float)
        v['sigma_lowest'] = float(np.asarray(d.variables['lev'][:])[-1])
    return v


def results():
    files = {}
    sums = None
    progress = RUN_DIR / 'progress.json'
    config = json.loads(progress.read_text())['configuration']
    files[str(progress.relative_to(ROOT))] = dict(sha256=sha256(progress), bytes=progress.stat().st_size)
    nbins = int(360 / HOUR_BIN_DEG)
    nbands = int(90 / BAND_DEG)
    for year in YEARS:
        path = RUN_DIR / 'model' / f'MOST.{year:05d}.nc'
        files[str(path.relative_to(ROOT))] = dict(sha256=sha256(path), bytes=path.stat().st_size)
        y = read_year(path)
        lat, lon = y['lat'], y['lon']
        land = y['lsm'] > 0.5
        ghi = y['rss'] / (1.0 - y['alb'])
        ghi_clear = y['rsscs'] / (1.0 - y['alb'])
        toa = y['rst'] - y['rsut']
        speed = np.hypot(y['u'], y['v'])
        dissipation = y['tauu'] * y['u'] + y['tauv'] * y['v']
        sun = sun_longitude(y['czen'], lat, lon)
        hour = (lon[None, :] - sun[:, None] + 180.0) % 360.0 - 180.0          # (time, lon), noon 0, morning negative
        if sums is None:
            nrow = lat.size
            sums = dict(n=0, lat=lat, lon=lon, land=land, sigma=y['sigma_lowest'],
                        rows={k: np.zeros(nrow) for k in ('ghi', 'ghi_clear', 'toa', 'clt', 'tas_sun', 'sun_w', 'prw',
                                                         'speed', 'cube', 'below', 'diss')},
                        land_rows={k: np.zeros(nrow) for k in ('ghi', 'ghi_clear', 'toa', 'clt', 'tas_sun', 'sun_w',
                                                              'speed', 'cube', 'below', 'diss')},
                        water_rows={k: np.zeros(nrow) for k in ('ghi', 'toa', 'clt', 'speed', 'cube', 'alb')},
                        profile=np.zeros((nbands, nbins)), profile_n=np.zeros((nbands, nbins)),
                        profile_clear=np.zeros((nbands, nbins)))
        s = sums
        nt = ghi.shape[0]
        s['n'] += nt
        nland = np.maximum(land.sum(axis=1), 1)
        nwater = np.maximum((~land).sum(axis=1), 1)

        def acc(target, field, mask=None, count=None):
            if mask is None:
                return target + field.sum(axis=0).mean(axis=1)
            return target + (field * mask[None]).sum(axis=0).sum(axis=1) / count
        R = s['rows']
        for key, field in (('ghi', ghi), ('ghi_clear', ghi_clear), ('toa', toa), ('clt', y['clt']), ('prw', y['prw']),
                           ('speed', speed), ('cube', speed ** 3), ('below', (speed < CUT_IN_M_S).astype(float)),
                           ('diss', dissipation), ('sun_w', ghi), ('tas_sun', y['tas'] * ghi)):
            R[key] = acc(R[key], field)
        L = s['land_rows']
        for key, field in (('ghi', ghi), ('ghi_clear', ghi_clear), ('toa', toa), ('clt', y['clt']), ('speed', speed),
                           ('cube', speed ** 3), ('below', (speed < CUT_IN_M_S).astype(float)), ('diss', dissipation),
                           ('sun_w', ghi), ('tas_sun', y['tas'] * ghi)):
            L[key] = acc(L[key], field, land, nland)
        W = s['water_rows']
        for key, field in (('ghi', ghi), ('toa', toa), ('clt', y['clt']), ('speed', speed), ('cube', speed ** 3),
                           ('alb', y['alb'])):
            W[key] = acc(W[key], field, ~land, nwater)
        # The land's day profile by band (north and south folded) and hour angle.
        weight = np.cos(np.radians(lat))
        band = np.minimum((np.abs(lat) // BAND_DEG).astype(int), nbands - 1)
        hb = ((hour + 180.0) // HOUR_BIN_DEG).astype(int) % nbins                 # (time, lon)
        for j in range(lat.size):
            cells = land[j]
            if not cells.any():
                continue
            idx = hb[:, cells].ravel()
            np.add.at(s['profile'][band[j]], idx, weight[j] * ghi[:, j, cells].ravel())
            np.add.at(s['profile_clear'][band[j]], idx, weight[j] * ghi_clear[:, j, cells].ravel())
            np.add.at(s['profile_n'][band[j]], idx, weight[j])
    s = sums
    n = s['n']
    lat = s['lat']
    weight = np.cos(np.radians(lat))
    rows = {k: v / n for k, v in s['rows'].items()}
    land_rows = {k: v / n for k, v in s['land_rows'].items()}
    water_rows = {k: v / n for k, v in s['water_rows'].items()}
    land_cells = s['land'].sum(axis=1)
    has_land = land_cells > 0

    def wmean(x, mask=None, w=None):
        w = weight if w is None else w
        m = np.ones_like(x, bool) if mask is None else mask
        return float((x * w)[m].sum() / w[m].sum())
    land_w = weight * land_cells                         # land area weights by row
    water_w = weight * (s['lon'].size - land_cells)
    rows['tas_sun'] = rows['tas_sun'] / np.maximum(rows['sun_w'], 1e-9)
    land_rows['tas_sun'] = land_rows['tas_sun'] / np.maximum(land_rows['sun_w'], 1e-9)
    by_row = dict(
        lat_deg=[round(float(x), 3) for x in lat],
        land_cells=[int(x) for x in land_cells],
        all=dict(ghi_w_m2=rows['ghi'], ghi_clear_w_m2=rows['ghi_clear'], toa_w_m2=rows['toa'], cloud=rows['clt'],
                 water_column_kg_m2=rows['prw'], wind_m_s=rows['speed'], wind_cube_m3_s3=rows['cube'],
                 share_below_cut_in=rows['below'], stress_dissipation_w_m2=rows['diss'],
                 air_k_weighted_by_sunlight=rows['tas_sun']),
        land=dict(ghi_w_m2=land_rows['ghi'], ghi_clear_w_m2=land_rows['ghi_clear'], toa_w_m2=land_rows['toa'],
                  cloud=land_rows['clt'], wind_m_s=land_rows['speed'], wind_cube_m3_s3=land_rows['cube'],
                  share_below_cut_in=land_rows['below'], stress_dissipation_w_m2=land_rows['diss'],
                  air_k_weighted_by_sunlight=land_rows['tas_sun']),
        water=dict(ghi_w_m2=water_rows['ghi'], toa_w_m2=water_rows['toa'], cloud=water_rows['clt'],
                   wind_m_s=water_rows['speed'], wind_cube_m3_s3=water_rows['cube']))
    for part in by_row.values():
        if isinstance(part, dict):
            for k, v in part.items():
                part[k] = [round(float(x), 5) if np.isfinite(x) else None for x in np.asarray(v)]
    profile = np.where(s['profile_n'] > 0, s['profile'] / np.maximum(s['profile_n'], 1e-12), np.nan)
    profile_clear = np.where(s['profile_n'] > 0, s['profile_clear'] / np.maximum(s['profile_n'], 1e-12), np.nan)
    nbands = profile.shape[0]
    bands = {}
    for b in range(nbands):
        lo = b * BAND_DEG
        rows_in = (np.abs(lat) >= lo) & (np.abs(lat) < lo + BAND_DEG) & has_land
        if not rows_in.any():
            continue
        bands[f'{lo:g}'] = dict(
            lat_deg=[lo, lo + BAND_DEG],
            land_rows=int(rows_in.sum()),
            land_ghi_w_m2=wmean(np.asarray(land_rows['ghi']), rows_in, land_w),
            land_ghi_clear_w_m2=wmean(np.asarray(land_rows['ghi_clear']), rows_in, land_w),
            land_toa_w_m2=wmean(np.asarray(land_rows['toa']), rows_in, land_w),
            land_cloud=wmean(np.asarray(land_rows['clt']), rows_in, land_w),
            land_wind_m_s=wmean(np.asarray(land_rows['speed']), rows_in, land_w),
            land_wind_cube_m3_s3=wmean(np.asarray(land_rows['cube']), rows_in, land_w),
            land_share_below_cut_in=wmean(np.asarray(land_rows['below']), rows_in, land_w),
            land_air_k_weighted_by_sunlight=wmean(np.asarray(land_rows['tas_sun']), rows_in, land_w),
            day_profile_ghi_w_m2=[round(float(x), 3) for x in profile[b]],
            day_profile_ghi_clear_w_m2=[round(float(x), 3) for x in profile_clear[b]])
    whole = dict(
        toa_w_m2=wmean(np.asarray(rows['toa'])),
        ghi_w_m2=wmean(np.asarray(rows['ghi'])),
        ghi_clear_w_m2=wmean(np.asarray(rows['ghi_clear'])),
        cloud=wmean(np.asarray(rows['clt'])),
        water_column_kg_m2=wmean(np.asarray(rows['prw'])),
        land=dict(ghi_w_m2=wmean(np.asarray(land_rows['ghi']), has_land, land_w),
                  ghi_clear_w_m2=wmean(np.asarray(land_rows['ghi_clear']), has_land, land_w),
                  toa_w_m2=wmean(np.asarray(land_rows['toa']), has_land, land_w),
                  cloud=wmean(np.asarray(land_rows['clt']), has_land, land_w),
                  wind_m_s=wmean(np.asarray(land_rows['speed']), has_land, land_w),
                  wind_cube_m3_s3=wmean(np.asarray(land_rows['cube']), has_land, land_w),
                  share_below_cut_in=wmean(np.asarray(land_rows['below']), has_land, land_w),
                  stress_dissipation_w_m2=wmean(np.asarray(land_rows['diss']), has_land, land_w),
                  air_k_weighted_by_sunlight=wmean(np.asarray(land_rows['tas_sun']), has_land,
                                                   land_w * np.asarray(land_rows['sun_w']))),
        water=dict(ghi_w_m2=wmean(np.asarray(water_rows['ghi']), None, water_w),
                   albedo=wmean(np.asarray(water_rows['alb']), None, water_w),
                   cloud=wmean(np.asarray(water_rows['clt']), None, water_w),
                   wind_m_s=wmean(np.asarray(water_rows['speed']), None, water_w),
                   wind_cube_m3_s3=wmean(np.asarray(water_rows['cube']), None, water_w)),
        stress_dissipation_w_m2=wmean(np.asarray(rows['diss'])),
        land_share_of_cells=float((land_w.sum()) / (weight.sum() * s['lon'].size)))
    return dict(
        schema=SCHEMA,
        producer=dict(domain='provisioning', files={p: sha256(ROOT / p) for p in CODE}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        run=dict(name=RUN, years=list(YEARS), records=n, sunlight_scale=config['model']['sunlight_scale'],
                 pressure_pa=config['pressure_pa'], land_albedo=config['model']['land_albedo'],
                 lowest_sigma=s['sigma'], files=files),
        method=dict(hour_bin_deg=HOUR_BIN_DEG, band_deg=BAND_DEG, equator_rows_deg=EQUATOR_ROWS_DEG,
                    cut_in_m_s=CUT_IN_M_S,
                    ghi='rss / (1 - alb)', toa='rst - rsut (rsut signed negative)',
                    hour_angle='longitude less the record subsolar longitude from czen over |lat| < 15 deg'),
        whole=whole, bands=bands, by_row=by_row)


def dumps(product) -> str:
    return json.dumps(product, indent=1) + '\n'


def main(argv=None) -> int:
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(product))
    print(json.dumps(product['whole'], indent=1))
    for k, b in product['bands'].items():
        print(k, round(b['land_ghi_w_m2'], 1), round(b['land_toa_w_m2'], 1), round(b['land_cloud'], 3),
              round(b['land_wind_m_s'], 2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
