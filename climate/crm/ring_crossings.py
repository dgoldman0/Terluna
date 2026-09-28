"""Where two rings along great circles cross: one place simulated in two vertical slices facing different ways,
under the same Sun, compared by local time.

    climate/gcm/.venv/bin/python -m climate.crm.ring_crossings ring_a ring_a_prime --from-day 29.5

writes the comparison to ../results/crm/crossings_<ring>_<ring>.json.

Two great circles cross at two opposite points. At each, the analysis takes each ring's columns within 100 km of
the crossing along the ring, land and water separately, in every 3-hour snapshot of the analysed span: the air at
2 m with its dewpoint and humidity, the dewpoint at the height of the GCM's lowest layer, rain, cloud, fog, cloud in
the flight band, the wind at 10 m and at 40 km, and whether the air is comfortable to be out in. Independent runs do not line up storm by storm, so the rings meet
through composites by local time and through their means over the day, the night and the whole span. Every ring
takes its local time from longitude, so each snapshot pairs the two patches at the same local time; the paired
differences are resampled a day at a time (a moving-block bootstrap) for the range the weather of the span allows.
The GCM design case's cells near the crossing give a third view.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from climate.crm import ring_analysis as ra

RUNS, RESULTS = ra.RUNS, ra.RESULTS
SCHEMA = 'terluna.climate.crm-crossings/1'
RADIUS_M = 100.0e3                    # a ring's columns this near a crossing, along the ring, make its patch
GCM_RADIUS_M = 150.0e3                # GCM cells (170 km wide at the equator) whose centres lie this near
BINS = 10                             # 36 degrees of local time, the bins of the GCM comparison
MIN_COLUMNS = 5                       # a surface is compared where each ring has at least this many columns of it
BLOCK_SNAPSHOTS = 8                   # one Earth day of 3-hour snapshots, the unit the bootstrap resamples
RESAMPLES = 2000
FLIGHT_WIND_M = 40000.0
COMFORT = {'strict': (18.0, 26.0, 15.0), 'loose': (17.0, 27.0, 16.0)}   # air at 2 m from, to (C); highest dewpoint (C)
VARIABLES = ('air_c', 'dewpoint_c', 'dewpoint_gcm_layer_c', 'humidity', 'rain_mm_h', 'cloud_cover', 'fog',
             'cloud_in_flight_band', 'wind_10m_m_s', 'wind_40_km_m_s', 'comfortable_strict', 'comfortable_loose')
GCM_RUN, GCM_YEARS = 'A28_dim5', (15, 24)
EVIDENCE = ('Two CM1 rings along great circles, compared where they cross (each ring\'s own summary says what a ring '
            'is). At a crossing both simulate the same ground under the same Sun, each in a vertical slice facing a '
            'different way, so their agreement shows how much a result depends on the direction a two-dimensional '
            'model has to choose; errors every ring shares, such as the missing inflow from the sides, agree and stay '
            'hidden. The analysed span is one lunar day after a lunar day of spin-up, so the ranges reflect the '
            'weather of that one day. The GCM design case (A28_dim5, years 15-24, 3-day means, cells about 170 km '
            'wide) is a third view.')
READING_RULE = ('A ring\'s patch is its columns within the radius of the crossing along the ring, land and water '
                'separately; a surface is compared where each ring has at least 5 columns of it. Values are patch means of '
                'snapshots every 3 model hours, composited by the local hour angle at the crossing (0 = noon, negative = '
                'morning) in 36-degree bins; day is when the Sun is up there. Differences are the first ring minus the '
                'second, snapshot by snapshot at the same local time; p5 and p95 bound 90% of the means of series rebuilt '
                'from one-day blocks drawn at random (a moving-block bootstrap, 2,000 draws). Comfortable means air at 2 m '
                'of 18-26 C with a dewpoint of at most 15 C (strict), or 17-27 C and 16 C (loose); hours per lunar day are '
                'the comfortable share of the patch through the span times the lunar day. The GCM\'s cells are those whose '
                'centres lie within 150 km, land and sea-level water separately. Its air is its lowest layer carried to '
                'the ground along a dry adiabat, but its vapour is that of the layer itself, centred about 0.86 km up '
                '(sigma 0.983); dewpoint_gcm_layer_c is the rings\' dewpoint at that pressure, the like-for-like '
                'comparison. Its comfortable hours use that dewpoint and 3-day means, which smooth out the warmest and '
                'coolest hours.')


def ring_path(record: dict) -> dict:
    """The great circle a ring runs along: its tilt and the longitude where it heads north across the equator,
    x = 0. The equatorial ring runs east from 0 E, a ring of no tilt."""
    if record.get('path'):
        return dict(tilt_deg=float(record['path']['tilt_deg']), node_deg=float(record['path']['node_deg']))
    if record.get('latitude_deg', 0.0) == 0.0:
        return dict(tilt_deg=0.0, node_deg=0.0)
    raise ValueError(f"{record.get('case')} runs along a circle of latitude, which is not a great circle")


def _frame(path: dict):
    """Unit vectors to a ring's x = 0 and to the point a quarter of the ring further east."""
    tilt, node = np.radians(path['tilt_deg']), np.radians(path['node_deg'])
    start = np.array([np.cos(node), np.sin(node), 0.0])
    ahead = np.array([-np.sin(node) * np.cos(tilt), np.cos(node) * np.cos(tilt), np.sin(tilt)])
    return start, ahead


def crossings(path_a: dict, path_b: dict) -> list:
    """The two opposite points where two great circles cross, in order of longitude: latitude, longitude and
    where each ring passes them, as a share of its length from x = 0."""
    line = np.cross(np.cross(*_frame(path_a)), np.cross(*_frame(path_b)))
    size = np.linalg.norm(line)
    if size < 1e-9:
        raise ValueError('the two rings run along one great circle')
    found = []
    for point in (line / size, -line / size):
        along = []
        for path in (path_a, path_b):
            start, ahead = _frame(path)
            along.append(float(np.round(np.arctan2(point @ ahead, point @ start) / (2.0 * np.pi), 12)) % 1.0)
        found.append(dict(lat_deg=float(np.round(np.degrees(np.arcsin(np.clip(point[2], -1.0, 1.0))), 9)) + 0.0,
                          lon_deg=float(np.round(np.degrees(np.arctan2(point[1], point[0])), 9)) % 360.0, along=along))
    return sorted(found, key=lambda c: c['lon_deg'])


def heading_deg(path: dict, along: float) -> float:
    """The direction a ring runs (eastward along x) at a point on it, clockwise from north."""
    start, ahead = _frame(path)
    s = 2.0 * np.pi * along
    point, velocity = np.cos(s) * start + np.sin(s) * ahead, -np.sin(s) * start + np.cos(s) * ahead
    lat, lon = np.arcsin(np.clip(point[2], -1.0, 1.0)), np.arctan2(point[1], point[0])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    return float(np.degrees(np.arctan2(velocity @ east, velocity @ north))) % 360.0


def along_ring(path: dict, lat_deg: float, lon_deg: float) -> float:
    """Where on a ring a point lies (or the point of the ring nearest it), as a share of its length from x = 0."""
    start, ahead = _frame(path)
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    point = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    return float(np.round(np.arctan2(point @ ahead, point @ start) / (2.0 * np.pi), 12)) % 1.0


def harmonic_fit(values, centres_deg, harmonics: int) -> np.ndarray:
    """A periodic series (bins along the first axis) refitted, level by level, with its mean and the first
    `harmonics` harmonics of the day, by least squares."""
    a = np.radians(np.asarray(centres_deg, dtype=float))
    basis = np.column_stack([np.ones_like(a)] + [f(n * a) for n in range(1, harmonics + 1) for f in (np.cos, np.sin)])
    coef, *_ = np.linalg.lstsq(basis, np.asarray(values, dtype=float), rcond=None)
    return basis @ coef


def day_night_advection(names, lat_deg: float, lon_deg: float, from_day: float, wavenumbers: int, bins: int,
                        harmonics: int, top_m: float, radius_m: float = RADIUS_M) -> dict:
    """The heating and moistening the planet-wide day-night circulation brings to a place on rings along great
    circles: on each ring, the advection along it by the flow of its longest waves (wavenumbers up to
    `wavenumbers` round the ring, which leaves out storms and storm systems), averaged over its columns within
    radius_m of the place, by local time in `bins` bins through the span from from_day; the rings averaged, each
    level refitted with the day's mean and first harmonics, and nothing above top_m."""
    from climate.crm.cm1_run import solar_day_s
    day_s = solar_day_s()
    edges = np.linspace(-180.0, 180.0, bins + 1)
    total_th, total_qv, used = 0.0, 0.0, []
    for name in names:
        case = RUNS / name
        geo = ra.case_geometry(case)
        record = geo['record']
        tap, start = record['configuration']['output_s'], record['configuration']['start_hour_angle_deg']
        nx, dx, zh = record['grid']['nx'], geo['dx'], geo['zh']
        cols = patch_columns(nx, along_ring(ring_path(record), lat_deg, lon_deg), radius_m, dx)
        k = np.fft.rfftfreq(nx, d=1.0 / nx)
        keep = (k <= wavenumbers)[None, :]
        slope = keep * (2j * np.pi * k / (nx * dx))[None, :]
        low = lambda f: np.fft.irfft(np.fft.rfft(f, axis=-1) * keep, n=nx, axis=-1)[:, cols]
        ddx = lambda f: np.fft.irfft(np.fft.rfft(f, axis=-1) * slope, n=nx, axis=-1)[:, cols]
        th, qv, count = np.zeros((bins, zh.size)), np.zeros((bins, zh.size)), np.zeros(bins)
        outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
        for n in [n for n in outputs if (n - 1) * tap >= from_day * 86400.0]:
            d = ra.read_snapshot(case, n)
            u = low(d['uinterp'])
            b = min(int(np.digitize(local_hour_angle(start, lon_deg, (n - 1) * tap, day_s), edges)) - 1, bins - 1)
            th[b] -= (u * ddx(d['th'])).mean(axis=1)
            qv[b] -= (u * ddx(d['qv'])).mean(axis=1)
            count[b] += 1
        if (count == 0).any():
            raise RuntimeError(f'{name}: no snapshot in some local-time bins after day {from_day}')
        total_th, total_qv = total_th + th / count[:, None], total_qv + qv / count[:, None]
        used.append(dict(case=name, columns=[int(c) for c in cols], snapshots=int(count.sum())))
    centres = 0.5 * (edges[1:] + edges[:-1])
    below = zh <= top_m
    fit = lambda v: np.where(below[None, :], harmonic_fit(v / len(names), centres, harmonics), 0.0)
    return dict(z_m=zh, hour_angle_deg=centres, theta_k_s=fit(total_th), qv_kg_kg_s=fit(total_qv), rings=used,
                wavenumbers=wavenumbers, shortest_wave_km=nx * dx / wavenumbers / 1000.0, harmonics=harmonics,
                top_m=top_m, from_day=from_day)


def patch_columns(nx: int, along: float, radius_m: float, dx: float) -> np.ndarray:
    """Columns whose centres lie within radius_m, along the ring, of the point at `along` (share of the length)."""
    length = nx * dx
    gap = ((np.arange(nx) + 0.5) * dx - along * length + 0.5 * length) % length - 0.5 * length
    return np.flatnonzero(np.abs(gap) <= radius_m)


def local_hour_angle(start_deg: float, lon_deg: float, t_s: float, day_s: float) -> float:
    """Hour angle (degrees, -180..180, 0 = noon) at a longitude at model time t: every ring starts the Sun with
    local time start_deg at 0 E."""
    return float((start_deg + lon_deg + 360.0 * t_s / day_s + 180.0) % 360.0 - 180.0)


def dewpoint_c(e):
    """Dewpoint (C) of air with vapour pressure e (Pa): the inverse of ring_analysis.saturation_pa."""
    g = np.log(np.maximum(e, 1.0e-3) / 611.2)
    return 243.5 * g / (17.67 - g)


def at_sigma(field, prs, psfc, sigma: float):
    """A section's values (nz, nx) where the pressure is sigma times the surface pressure in each column,
    interpolated linearly in log pressure between the levels either side."""
    target = sigma * psfc
    above = prs < target[None, :]
    k1 = np.maximum(np.where(above.any(axis=0), np.argmax(above, axis=0), prs.shape[0] - 1), 1)   # the first level above it
    k0, cols = k1 - 1, np.arange(prs.shape[1])
    w = np.log(prs[k0, cols] / target) / np.log(prs[k0, cols] / prs[k1, cols])
    return field[k0, cols] + w * (field[k1, cols] - field[k0, cols])


def gcm_lowest_sigma(run=GCM_RUN, year=GCM_YEARS[0]) -> float:
    """Sigma of the GCM's lowest full level, where its vapour lives."""
    import netCDF4
    with netCDF4.Dataset(ra.HERE.parent / 'gcm' / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
        return float(d['lev'][-1])


def comfortable(air_c, dew_c, band: str):
    low, high, dew_max = COMFORT[band]
    return (air_c >= low) & (air_c <= high) & (dew_c <= dew_max)


def block_bootstrap(values, block=BLOCK_SNAPSHOTS, resamples=RESAMPLES, seed=0) -> dict:
    """The mean of a series in time, with the 5th and 95th percentiles of the means of series rebuilt from runs
    of `block` consecutive values starting at random (wrapping round), which keeps the weather's persistence."""
    v = np.asarray(values, dtype=float)
    n = v.size
    if n == 0:
        return None
    starts = np.random.default_rng(seed).integers(0, n, size=(resamples, -(-n // block)))
    means = v[((starts[:, :, None] + np.arange(block)) % n).reshape(resamples, -1)[:, :n]].mean(axis=1)
    return dict(mean=float(v.mean()), p5=float(np.percentile(means, 5)), p95=float(np.percentile(means, 95)))


def patch_series(case: Path, patches: list, from_day: float, to_day: float | None, sigma: float) -> dict:
    """Model times from from_day (to to_day) and, per snapshot, the mean of every variable over each patch's
    columns. patches holds, for each crossing, the land and the water columns near it; sigma is the GCM's lowest
    level."""
    geo = ra.case_geometry(case)
    tap = geo['record']['configuration']['output_s']
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    use = [n for n in outputs if from_day * 86400.0 <= (n - 1) * tap <= (np.inf if to_day is None else to_day * 86400.0)]
    zh = geo['zh']
    near_ground = zh < ra.FOG_TOP_M
    flight = (zh >= ra.FLIGHT_BAND_M[0]) & (zh <= ra.FLIGHT_BAND_M[1])
    k40 = ra.level_index(zh, FLIGHT_WIND_M)
    out = dict(time_s=[], patches=[{surface: {key: [] for key in VARIABLES} for surface in cols} for cols in patches])
    for n in use:
        d = ra.read_snapshot(case, n)
        e2 = ra.vapour_pa(d['q2'], d['psfc'])
        air, dew = d['t2'] - 273.15, dewpoint_c(e2)
        cloudy = (d['qc'] + d['qi']) >= ra.CLOUD_KG_KG
        u40 = d['uinterp'][k40]
        layer_p = sigma * d['psfc']
        layer_dew = dewpoint_c(ra.vapour_pa(at_sigma(d['qv'], d['prs'], d['psfc'], sigma), layer_p))
        fields = dict(air_c=air, dewpoint_c=dew, dewpoint_gcm_layer_c=layer_dew, humidity=e2 / ra.saturation_pa(d['t2']), rain_mm_h=d['prate'] * 3600.0,
                      cloud_cover=cloudy.any(axis=0), fog=cloudy[near_ground].any(axis=0),
                      cloud_in_flight_band=cloudy[flight].any(axis=0), wind_10m_m_s=d['s10'],
                      wind_40_km_m_s=np.hypot(u40, d['vinterp'][k40]) if 'vinterp' in d else np.abs(u40),
                      comfortable_strict=comfortable(air, dew, 'strict'), comfortable_loose=comfortable(air, dew, 'loose'))
        out['time_s'].append((n - 1) * tap)
        for cols, store in zip(patches, out['patches']):
            for surface, idx in cols.items():
                for key in VARIABLES:
                    store[surface][key].append(float(np.mean(fields[key][idx])))
    return out


def compare(hour_angle, series: dict, lunar_day_h: float) -> dict:
    """One surface at one crossing, from one or two rings' patch means at the same local times: each ring's
    composites by local time, its means over the day, the night and the span, its comfortable hours and rain per
    lunar day, and with two rings the differences (first minus second) with their bootstrap ranges."""
    h = np.asarray(hour_angle, dtype=float)
    k = np.clip(np.digitize(h, np.linspace(-180.0, 180.0, BINS + 1)) - 1, 0, BINS - 1)
    day = np.cos(np.radians(h)) > 0.0
    names = list(series)
    binned = lambda v: [float(v[k == j].mean()) if (k == j).any() else None for j in range(BINS)]
    mean = lambda v: float(v.mean()) if v.size else None
    spans = lambda v: dict(all=mean(v), day=mean(v[day]), night=mean(v[~day]))
    out = dict(hour_angle_deg=[float(c) for c in (np.arange(BINS) + 0.5) * 360.0 / BINS - 180.0])
    for key in VARIABLES:
        values = {name: np.asarray(series[name][key], dtype=float) for name in names}
        row = dict(by_local_time={name: binned(v) for name, v in values.items()},
                   mean={name: spans(v) for name, v in values.items()})
        if len(names) == 2:
            diff = values[names[0]] - values[names[1]]
            by_bin = np.array([np.nan if b is None else b for b in binned(diff)])
            row['difference'] = dict(all=block_bootstrap(diff), day=block_bootstrap(diff[day]),
                                     night=block_bootstrap(diff[~day]),
                                     rms_by_local_time=float(np.sqrt(np.nanmean(by_bin ** 2))))
        out[key] = row
    out['comfortable_hours_per_lunar_day'] = {band: {name: float(np.mean(series[name][f'comfortable_{band}']) * lunar_day_h)
                                                     for name in names} for band in COMFORT}
    out['rain_mm_per_lunar_day'] = {name: float(np.mean(series[name]['rain_mm_h']) * lunar_day_h) for name in names}
    return out


def gcm_near(points: list, radius_m=GCM_RADIUS_M, run=GCM_RUN, years=GCM_YEARS) -> list:
    """The GCM design case near each crossing: its cells whose centres lie within radius_m, land and sea-level
    water separately, by local time in the same bins (3-day means): the air at 2 m, the dewpoint of its lowest
    layer, rain and cloud cover, and the share of 3-day means comfortable to be out in, as hours per lunar day."""
    import netCDF4
    from climate.gcm.climatology import subsolar_longitudes
    from climate.crm.cm1_run import planet, solar_day_s
    folder = ra.HERE.parent / 'gcm' / 'runs' / run / 'model'
    get = {k: [] for k in ('tas', 'pr', 'clt', 'czen', 'ps', 'hus')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            for k in get:
                get[k].append(np.asarray(d[k][:, -1] if k == 'hus' else d[k][:], dtype=float))   # the lowest layer's vapour
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            lsm = np.asarray(d['lsm'][:], float)[0] > 0.5
    f = {k: np.concatenate(v) for k, v in get.items()}
    radius = planet()['radius_m']
    height = np.loadtxt(folder.parent / 'inputs' / 'moon_topography.sra', skiprows=1).ravel().reshape(lsm.shape) / planet()['gravity_m_s2']
    q, p = f['hus'], f['ps'] * 100.0
    values = dict(air_c=f['tas'] - 273.15, dewpoint_c=dewpoint_c(q * p / (0.622 + 0.378 * q)), rain_mm_h=f['pr'] * 3.6e6,
                  cloud_cover=f['clt'])
    hour = (lon[None, :] - subsolar_longitudes(f['czen'], lat, lon)[:, None] + 180.0) % 360.0 - 180.0     # (t, lon)
    k = np.broadcast_to(np.clip(((hour + 180.0) / (360.0 / BINS)).astype(int), 0, BINS - 1)[:, None, :], f['tas'].shape)
    day = np.broadcast_to((np.cos(np.radians(hour)) > 0.0)[:, None, :], f['tas'].shape)
    lat2, lon2 = np.radians(np.meshgrid(lat, lon, indexing='ij'))
    lunar_day_h = solar_day_s() / 3600.0
    out = []
    for point in points:
        la, lo = np.radians(point['lat_deg']), np.radians(point['lon_deg'])
        angle = np.arccos(np.clip(np.sin(lat2) * np.sin(la) + np.cos(lat2) * np.cos(la) * np.cos(lon2 - lo), -1.0, 1.0))
        near = angle * radius <= radius_m
        entry = dict(run=run, years=list(years), radius_km=radius_m / 1000.0)
        for surface, cells in (('land', near & lsm), ('water', near & ~lsm & (height < 1.0))):
            if not cells.any():
                entry[surface] = None
                continue
            sel = np.broadcast_to(cells, f['tas'].shape)
            rows = dict(cells=[[float(lat[i]), float(lon[j])] for i, j in zip(*np.nonzero(cells))])
            for key, v in values.items():
                rows[key] = dict(by_local_time=[float(v[sel & (k == b)].mean()) for b in range(BINS)],
                                 mean=dict(all=float(v[sel].mean()), day=float(v[sel & day].mean()),
                                           night=float(v[sel & ~day].mean())))
            rows['comfortable_hours_per_lunar_day'] = {
                band: float(comfortable(values['air_c'][sel], values['dewpoint_c'][sel], band).mean() * lunar_day_h)
                for band in COMFORT}
            rows['rain_mm_per_lunar_day'] = float(values['rain_mm_h'][sel].mean() * lunar_day_h)
            entry[surface] = rows
        out.append(entry)
    return out


def analyse(names, from_day: float, radius_m: float = RADIUS_M, gcm: bool = True, to_day: float | None = None) -> dict:
    from climate.crm.cm1_run import solar_day_s
    geos = [ra.case_geometry(RUNS / name) for name in names]
    records = [g['record'] for g in geos]
    for key in ('start_hour_angle_deg', 'output_s'):
        if len({r['configuration'][key] for r in records}) != 1:
            raise ValueError(f'the rings differ in {key}')
    start = records[0]['configuration']['start_hour_angle_deg']
    paths = [ring_path(r) for r in records]
    points = crossings(*paths)
    patches = []                                                       # per ring, per crossing: surface -> columns
    for j, g in enumerate(geos):
        per = []
        for point in points:
            cols = patch_columns(g['record']['grid']['nx'], point['along'][j], radius_m, g['dx'])
            per.append({surface: cols[sel[cols]] for surface, sel in (('land', g['land']), ('water', ~g['land']))
                        if sel[cols].any()})
        patches.append(per)
    sigma = gcm_lowest_sigma()
    series = [patch_series(RUNS / name, per, from_day, to_day, sigma) for name, per in zip(names, patches)]
    times = sorted(set(series[0]['time_s']) & set(series[1]['time_s']))
    if not times:
        raise RuntimeError(f'no common output after day {from_day}')
    rows = [[s['time_s'].index(t) for t in times] for s in series]
    day_s = solar_day_s()
    near_gcm = gcm_near(points) if gcm else [None] * len(points)
    found = []
    for i, point in enumerate(points):
        entry = dict(lat_deg=point['lat_deg'], lon_deg=point['lon_deg'], rings={})
        for j, name in enumerate(names):
            nx = records[j]['grid']['nx']
            entry['rings'][name] = dict(column=int(point['along'][j] * nx) % nx,
                                        heading_deg=heading_deg(paths[j], point['along'][j]),
                                        **{f'{surface}_columns': int(len(patches[j][i].get(surface, ()))) for surface in ('land', 'water')})
        hour = [local_hour_angle(start, point['lon_deg'], t, day_s) for t in times]
        for surface in ('land', 'water'):
            have = [name for j, name in enumerate(names) if len(patches[j][i].get(surface, ())) >= MIN_COLUMNS]
            entry[surface] = compare(hour, {name: {key: np.asarray(series[j]['patches'][i][surface][key])[rows[j]]
                                                   for key in VARIABLES}
                                            for j, name in enumerate(names) if name in have},
                                     day_s / 3600.0) if have else None
        entry['gcm'] = near_gcm[i]
        found.append(entry)
    return dict(schema=SCHEMA, rings=list(names), evidence=EVIDENCE, reading_rule=READING_RULE,
                span_days=[times[0] / 86400.0, times[-1] / 86400.0], snapshots=len(times), radius_km=radius_m / 1000.0,
                gcm_lowest_sigma=sigma,
                paths={name: path for name, path in zip(names, paths)}, crossings=found)


def rounded(value, digits=3):
    if isinstance(value, float):
        return round(value, digits)
    if isinstance(value, dict):
        return {k: rounded(v, digits) for k, v in value.items()}
    if isinstance(value, list):
        return [rounded(v, digits) for v in value]
    return value


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('rings', nargs=2)
    parser.add_argument('--from-day', type=float, default=29.5)
    parser.add_argument('--to-day', type=float, help='the last model day to use (default: all output)')
    parser.add_argument('--radius-km', type=float, default=RADIUS_M / 1000.0)
    parser.add_argument('--output', type=Path, help='where to write the comparison (default: the results folder)')
    args = parser.parse_args(argv)
    result = analyse(args.rings, args.from_day, args.radius_km * 1000.0, to_day=args.to_day)
    result['producer'] = dict(domain='climate', files={'crm/ring_crossings.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    path = args.output or RESULTS / f'crossings_{args.rings[0]}_{args.rings[1]}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rounded(result), indent=1) + '\n')
    print(f'{path}: {result["snapshots"]} snapshots, days {result["span_days"][0]:.1f}-{result["span_days"][1]:.1f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
