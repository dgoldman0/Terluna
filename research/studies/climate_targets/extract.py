"""What the study reads from the GCM's own output: the design case and the 20-year test of the longwave refitted by
height with ice optics for the cold cloud, over years 10-19 of each, and the test's year-by-year progress.

    climate/gcm/.venv/bin/python -m research.studies.climate_targets.extract     # about a minute; needs netCDF4

It reads the ExoPlaSim output in climate/gcm/runs (ignored by Git) and the T21 boundary file in climate/gcm/products,
and writes research/runs/climate_targets/extract.json (ignored), with the sha256 of every file it read; run.py binds
it. The heat is in the 3-day means, each about a tenth of the lunar day, at two heights. One is the lowest model
level, about 0.9 km above the ground. The other is 2 m, where the output gives the air's temperature and the wet-bulb
temperature takes the lowest level's specific humidity, capped at saturation, as in air mixed through the lowest
kilometre. Every figure is weighted by area, its quantiles included: the land's by each land cell's land area (land
fraction above one half), the seas' by each sea cell's sea area. Each place's warmest and coolest tenth average every
cell's own warmest and coolest tenth of its 3-day means. The coasts are land cells with a sea cell among their eight
neighbours, joined across the longitude seam and kept apart across the poles. The land's wet-bulb figures counted cell
by cell, as the readout of 17 October counted them, stand beside the weighted ones under their own name.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / 'research' / 'runs' / 'climate_targets' / 'extract.json'
SCHEMA = 'terluna.research.climate-targets-extract/2'
RUNS_DIR = ROOT / 'climate' / 'gcm' / 'runs'
BOUNDARY = ROOT / 'climate' / 'gcm' / 'products' / 'moon_28pct_water_lakes_gaussian_T21.nc'
DESIGN, TEST = 'A28_corrected_970', 'A28_corrected_970_lw_ice'
YEARS = (10, 19)
BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))
HEIGHTS_KM = (0.9, 2, 5, 10, 15, 20, 25, 30, 35, 40)
DESERT_MM_DAY = 0.5                  # the desert measure of 14 October: land under 0.5 mm of rain a day
# Saturation over water (Bolton 1980) and the psychrometric balance; the dry gas constant and gravity as the GCM's.
RD, GRAVITY, CP = 288.3767, 1.6242188765608128, 1004.64
FREEZING_K = 273.15


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def output(run, year):
    return RUNS_DIR / run / 'model' / f'MOST.{year:05d}.nc'


def es_water(t):
    return 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))


def qsat(t, p):
    e = es_water(t)
    return 0.622 * e / (p - 0.378 * e)


def wet_bulb(t, q, p):
    """Psychrometric wet-bulb temperature (K) by bisection."""
    lo, hi = np.full_like(t, 200.0), t.copy()
    for _ in range(40):
        tw = 0.5 * (lo + hi)
        lv = 2.501e6 - 2370.0 * (tw - 273.15)
        f = qsat(tw, p) - q - (CP + 1850.0 * q) / lv * (t - tw)
        lo, hi = np.where(f < 0, tw, lo), np.where(f < 0, hi, tw)
    return 0.5 * (lo + hi)


def boundary():
    import netCDF4
    with netCDF4.Dataset(BOUNDARY) as d:
        return dict(lat=np.asarray(d['lat'][:], float), land=np.asarray(d['land_fraction'][:], float),
                    elevation=np.asarray(d['land_elevation_m'][:], float),
                    area=np.asarray(d['cell_area_fraction'][:], float))


def wquantile(x, w, q):
    """The quantile q (0-1) of values x under weights w, linear between the weights' midpoints."""
    i = np.argsort(x, kind='stable')
    x, w = x[i], w[i]
    c = np.cumsum(w)
    return float(np.interp(q * c[-1], c - 0.5 * w, x))


def wstats(x, w):
    return dict(median=wquantile(x, w, 0.5), p99=wquantile(x, w, 0.99), top=float(x.max()))


def share(x, w, threshold):
    """The weighted share of values above a threshold."""
    return float((w * (x > threshold)).sum() / w.sum())


def beside(mask):
    """Cells with a true cell among their eight neighbours or in themselves. The rows run from pole to pole, so
    neighbours join across the longitude seam and never across a pole."""
    out = np.zeros_like(mask)
    for di in (-1, 0, 1):
        rows = np.zeros_like(mask)
        if di == 0:
            rows[:] = mask
        elif di > 0:
            rows[di:] = mask[:-di]
        else:
            rows[:di] = mask[-di:]
        for dj in (-1, 0, 1):
            out |= np.roll(rows, dj, 1)
    return out


def freezing_height(ta, z):
    """Each column's first height above the ground (m) at which the air falls through freezing, linear between
    levels, or the lowest level's height where that level already freezes. The levels run from the top down."""
    lowest = ta.shape[1] - 1
    k = lowest - np.argmax((ta < FREEZING_K)[:, ::-1], axis=1)
    below = np.minimum(k + 1, lowest)
    t1, t2 = (np.take_along_axis(ta, i[:, None], 1)[:, 0] for i in (k, below))
    z1, z2 = (np.take_along_axis(z, i[:, None], 1)[:, 0] for i in (k, below))
    return z2 + (t2 - FREEZING_K) / np.maximum(t2 - t1, 1e-6) * (z1 - z2)


def read(run, b):
    """Means and distributions over the run's years 10-19."""
    import netCDF4
    files, toa, fields = {}, {}, {k: [] for k in ('ta', 'hus', 'ps', 'ts', 'tas', 'pr')}
    for year in range(YEARS[0], YEARS[1] + 1):
        path = output(run, year)
        files[str(path.relative_to(ROOT))] = sha256(path)
        with netCDF4.Dataset(path) as d:
            lat = np.asarray(d['lat'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            for k in ('rst', 'rsut', 'rlut', 'rltcs', 'rstcs', 'ntr'):
                toa.setdefault(k, []).append(np.asarray(d[k][:], float).mean(axis=0))
            for k in fields:
                fields[k].append(np.asarray(d[k][:], float))
    nodes, weights = np.polynomial.legendre.leggauss(lat.size)       # the T21 rows' Gaussian weights
    w = weights[np.argsort(np.argsort(np.sin(np.radians(lat))))][:, None] * np.ones(toa['rst'][0].shape)
    gmean = lambda x: float((x * w).sum() / w.sum())
    t = {k: gmean(np.mean(v, axis=0)) for k, v in toa.items()}
    order = [int(np.argmin(np.abs(b['lat'] - x))) for x in lat]      # the boundary rows in the output's order
    land, elevation, area = b['land'][order], b['elevation'][order], b['area'][order]
    ta, q = np.concatenate(fields['ta']), np.concatenate(fields['hus'])
    ps, ts, pr = np.concatenate(fields['ps']) * 100.0, np.concatenate(fields['ts']), np.concatenate(fields['pr'])
    tas = np.concatenate(fields['tas'])
    t0, q0, p0 = ta[:, -1], q[:, -1], ps * sigma[-1]
    # wet-bulb and air temperatures (C) at the two heights in every 3-day mean
    hot = dict(lowest_level=(wet_bulb(t0, q0, p0) - 273.15, t0 - 273.15),
               at_2_m=(wet_bulb(tas, np.minimum(q0, qsat(tas, ps)), ps) - 273.15, tas - 273.15))
    sst = ts - 273.15
    alat = np.abs(lat)[:, None] * np.ones(land.shape)
    island, sea = land > 0.5, land < 0.5
    coast = island & beside(sea)
    wl, ws = land * area, (1.0 - land) * area          # each cell's land and sea area
    on_land = np.broadcast_to(island, sst.shape)
    # rain and deserts by land area
    rain = pr.mean(axis=0) * 86400.0 * 1000.0
    dry = rain < DESERT_MM_DAY

    def land_share(mask):
        ww = wl * mask
        return dict(desert=float((ww * dry).sum() / ww.sum()), rain_mm_day=float((ww * rain).sum() / ww.sum()))

    # every cell's own warmest and coolest tenth of its 3-day means
    tenth = max(1, sst.shape[0] // 10)
    tenths = {}
    for height, pair in hot.items():
        tenths[height] = []
        for x in pair:
            s = np.sort(x, axis=0)
            tenths[height].append((s[-tenth:].mean(axis=0), s[:tenth].mean(axis=0)))

    def heat_at(height, mask):
        (tw, air), ((tw_warm, tw_cool), (air_warm, air_cool)) = hot[height], tenths[height]
        m = np.broadcast_to(mask, tw.shape)
        x, a, ww = tw[m], air[m], np.broadcast_to(wl, tw.shape)[m]
        cells = wl * mask
        mean = lambda f: float((cells * f).sum() / cells.sum())
        return dict(wet_bulb=wstats(x, ww), above_28=share(x, ww, 28.0), above_31=share(x, ww, 31.0),
                    air=wstats(a, ww), warmest_tenth=dict(wet_bulb=mean(tw_warm), air=mean(air_warm)),
                    coolest_tenth=dict(wet_bulb=mean(tw_cool), air=mean(air_cool)))

    def heat(mask):
        out = heat_at('lowest_level', mask)
        out['at_2_m'] = heat_at('at_2_m', mask)
        return out

    def seas(mask):
        m = np.broadcast_to(mask, sst.shape)
        return wstats(sst[m], np.broadcast_to(ws, sst.shape)[m])

    land_heat = heat(island)
    x = hot['lowest_level'][0][on_land]
    land_heat['counted_by_cell'] = dict(wet_bulb=dict(median=float(np.median(x)), p99=float(np.percentile(x, 99))),
                                        above_28=float((x > 28.0).mean()), above_31=float((x > 31.0).mean()))
    out = dict(top_of_air_w_m2=dict(incoming=t['rst'] - t['rsut'], reflected=-t['rsut'], absorbed=t['rst'],
                                     absorbed_clear_sky=t['rstcs'], outgoing_infrared=-t['rlut'],
                                     outgoing_infrared_clear_sky=-t['rltcs'], net=t['ntr']),
               ground_k=gmean(ts.mean(axis=0)), rain_mm_day=gmean(rain),
               land=dict(all=land_share(np.ones_like(island)),
                         bands={f'{lo}-{hi}': land_share((alat >= lo) & (alat < hi)) for lo, hi in BANDS}),
               heat=dict(land=land_heat, coasts=heat(coast),
                         bands={f'{lo}-{hi}': heat(island & (alat >= lo) & (alat < hi)) for lo, hi in BANDS},
                         below_1_km=heat(island & (elevation < 1000.0)),
                         from_4_km=heat(island & (elevation >= 4000.0))),
               seas=dict(all=seas(sea),
                         bands={f'{lo}-{hi}': seas(sea & (alat >= lo) & (alat < hi)) for lo, hi in BANDS}))
    # the sky over land: each column's air by height above its ground, by land area, hypsometric heights
    z = np.zeros_like(ta)
    z[:, -1] = RD * ta[:, -1] / GRAVITY * np.log(1.0 / sigma[-1])
    for k in range(sigma.size - 2, -1, -1):
        z[:, k] = z[:, k + 1] + RD * 0.5 * (ta[:, k] + ta[:, k + 1]) / GRAVITY * np.log(sigma[k + 1] / sigma[k])
    tl, zl, psl, al = ta[:, :, island], z[:, :, island], ps[:, island], wl[island]
    ladder = []
    for h in HEIGHTS_KM:
        hz = h * 1000.0
        i = np.clip((zl > hz).sum(axis=1) - 1, 0, sigma.size - 2)
        z1, z2 = (np.take_along_axis(zl, (i + j)[:, None, :], 1)[:, 0] for j in (0, 1))
        t1, t2 = (np.take_along_axis(tl, (i + j)[:, None, :], 1)[:, 0] for j in (0, 1))
        f = (hz - z2) / np.where(z1 == z2, 1.0, z1 - z2)
        th = t2 + f * (t1 - t2)
        lp = np.log(sigma[i + 1]) + f * (np.log(sigma[i]) - np.log(sigma[i + 1]))
        ph = np.exp(lp) * psl / 101325.0
        ladder.append(dict(height_km=h, air_c=float((al * th).sum() / al.sum() / th.shape[0]) - 273.15,
                           pressure_atm=float((al * ph).sum() / al.sum() / ph.shape[0])))
    out['sky_over_land'] = ladder
    out['freezing_level_over_land_km'] = wquantile(freezing_height(ta, z)[on_land],
                                                   np.broadcast_to(wl, sst.shape)[on_land], 0.5) / 1000.0
    return out, files


def progress(run):
    path = RUNS_DIR / run / 'progress.json'
    p = json.loads(path.read_text())
    return ([dict(year=y['year'], ground_k=y['surface_k'], top_of_air_net_w_m2=y['toa_net_w_m2'])
             for y in p['years']], {str(path.relative_to(ROOT)): sha256(path)})


def extract():
    b = boundary()
    files = {str(BOUNDARY.relative_to(ROOT)): sha256(BOUNDARY)}
    runs = {}
    for run in (DESIGN, TEST):
        runs[run], f = read(run, b)
        files.update(f)
    years, f = progress(TEST)
    files.update(f)
    return dict(schema=SCHEMA, years=list(YEARS), runs=runs, test_progress=years,
                producer=dict(files={str(Path(__file__).resolve().relative_to(ROOT)): sha256(__file__)}, read=files))


def main(argv=None) -> int:
    out = extract()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
