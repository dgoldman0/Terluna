"""What the study reads from the GCM's own output: the design case and the 20-year test of the longwave refitted by
height with ice optics for the cold cloud, over years 10-19 of each, and the test's year-by-year progress; the glow
pair on the corrected radiation, from its runs' yearly records; the design case's two years written every
6.5 hours; and the air over the summit, where the summit tower stands, in those runs, in the test's every year and in
the earlier design case whose air the metropolis brief reads.

    climate/gcm/.venv/bin/python -m research.studies.climate_targets.extract     # about a minute; needs netCDF4

It reads the ExoPlaSim output in climate/gcm/runs (ignored by Git) and the T21 boundary file in climate/gcm/products,
and writes research/runs/climate_targets/extract.json (ignored), with the sha256 and size of every file it read;
run.py binds it.

Land and water are the model's own mask (lsm) for the heat, the seas and the sky. The runner turns cells to water,
least land first, until the planned water share is met, so cells with a land fraction a little above one half run as
water. Land cells are weighted by their land area (land fraction times cell area) and water cells by their water
area, the quantiles included. The coasts are land cells with a water cell among their eight neighbours, joined across
the longitude seam and kept apart across the poles. The land's rain and deserts are counted two ways: every cell's
land with that cell's own rain, since the rain falls on a cell's land as on its water, and the model's own land cells
alone.

The heat is in the 3-day means, each about a tenth of the lunar day, at two heights. One is the lowest model level,
about 0.9 km above the ground. The other is 2 m, where the output's air is the lowest level's carried down the dry
adiabat to the ground's pressure (PlaSim's mkshfl, with its kappa of 0.286), and the wet-bulb temperature carries the
lowest level's specific humidity down with it, capped at saturation there: the reading of a well-mixed lowest
kilometre, as by day. The output's humidity is the humidity at each mean's first step (window_humidity), so a mean's
humidity is that of its start and its end, the next record's, and each run's last mean leaves the heat's statistics.
The wet-bulb temperature is the psychrometric one over water with Bolton's saturation, capped at the air's own
temperature where the air is supersaturated. Each place's warmest and coolest tenth are each cell's own tenth of its
3-day means ranked by one quantity, with the others averaged over the same means.

The day's 3-day means are those the Sun is up for more than half of. The output's czen is each mean's average cosine
of the Sun's zenith angle, and a mean is the day's when its czen exceeds half_sunlit() of the Sun's travel through the
mean times the cosine of the latitude; the cut is exact where the Sun's declination is zero, and the fixed orbit
swings it through the tilt, 1.54 degrees either way. Of the day's means, those with the Sun past its highest are the
afternoon's. The design case's 6.5-hour years (A28_corrected_970_eddies, two years from the design case's year-19
restart, eleven means to each 3-day mean) give the land's heat in those means and in the 3-day means their elevens
make, over the same two years, each 6.5-hour mean's humidity from its start and its end as in the 3-day output.

The summit's air is read at the T21 cells round it. The output carries no orography, so the model's own ground comes
from its sea-level pressure, which pyburn extrapolated from that ground at the standard atmosphere's 6.5 K/km
(ground_geopotential); the levels' geopotential rises from it hypsometrically with the virtual temperature, and
heights take inverse-square gravity. Each column is interpolated linearly in geopotential to a ladder of heights
above the model's sea level, then weighted over the cells: bilinearly between the four cell centres round the summit,
the cell holding it alone, and that cell with its eight neighbours by area.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from atmosphere.radiative_convective.thermodynamics import R_VAPOUR
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / 'research' / 'runs' / 'climate_targets' / 'extract.json'
SCHEMA = 'terluna.research.climate-targets-extract/5'
RUNS_DIR = ROOT / 'climate' / 'gcm' / 'runs'
BOUNDARY = ROOT / 'climate' / 'gcm' / 'products' / 'moon_28pct_water_lakes_gaussian_T21.nc'
DESIGN, TEST = 'A28_corrected_970', 'A28_corrected_970_lw_ice'
YEARS = (10, 19)
# The glow pair on the corrected radiation, from the design case's year 19 with 13 October's added dimming: the gases
# and the glow with 4.22% of added dimming, and with 3.52% (climate/gcm README).
GLOW_PAIR = ('A28_corrected_gases_glow', 'A28_corrected_gases_glow_held')
# The design case's years written every 6.5 hours, 110 means a lunar day, from its year 19.
EDDIES, EDDY_YEARS = 'A28_corrected_970_eddies', (0, 1)
START = f'climate/gcm/runs/{DESIGN}/model/MOST_REST.00019'           # the design case's year-19 restart
BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))             # the heat's and the seas' bands of latitude
RAIN_BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90))
HEIGHT_CLASSES = (('below_2_km', None, 2000.0), ('2_to_4_km', 2000.0, 4000.0), ('from_4_km', 4000.0, None))
HEIGHTS_KM = (0.9, 2, 5, 10, 15, 20, 25, 30, 35, 40)
DESERT_MM_DAY = 0.5                  # the desert measure of 14 October: land under 0.5 mm of rain a day
# Saturation over water (Bolton 1980) and the psychrometric balance; the dry gas constant and gravity as the GCM's.
RD, GRAVITY, CP = 288.3767, 1.6242188765608128, 1004.64
AKAP = 0.286                         # PlaSim's kappa (p_exo.f90), with which mkshfl carries the air down to 2 m
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
    """Psychrometric wet-bulb temperature (K) by bisection between 200 K and the air's temperature, so supersaturated
    air takes the air's own temperature."""
    lo, hi = np.full_like(t, 200.0), t.copy()
    for _ in range(40):
        tw = 0.5 * (lo + hi)
        lv = 2.501e6 - 2370.0 * (tw - 273.15)
        f = qsat(tw, p) - q - (CP + 1850.0 * q) / lv * (t - tw)
        lo, hi = np.where(f < 0, tw, lo), np.where(f < 0, hi, tw)
    return 0.5 * (lo + hi)


def window_humidity(written):
    """Each mean's specific humidity from the output's records along time. ExoPlaSim 3.4.2 computes its output
    humidity only at the steps where mod(nstep, nafter) == 0 (plasim.f90, gridpointd, line 3218) and adds that held
    value at every step of the mean (outmod.f90, outaccu, line 2556), so each record holds the humidity at its mean's
    first step while the temperatures, the pressure and the Sun's cosine are true means. A mean's humidity is the mean
    of its start and its end, the next record's; the last record has no next one and gives none."""
    return 0.5 * (written[:-1] + written[1:])


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


def wmean(x, w):
    return float((w * x).sum() / w.sum())


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
    levels: the lowest level's height where that level already freezes, and NaN where no level freezes. The levels
    run from the top down."""
    lowest = ta.shape[1] - 1
    cold = ta < FREEZING_K
    k = lowest - np.argmax(cold[:, ::-1], axis=1)          # the first freezing level from the ground up
    below = np.minimum(k + 1, lowest)
    t1, t2 = (np.take_along_axis(ta, i[:, None], 1)[:, 0] for i in (k, below))
    z1, z2 = (np.take_along_axis(z, i[:, None], 1)[:, 0] for i in (k, below))
    height = z2 + (t2 - FREEZING_K) / np.maximum(t2 - t1, 1e-6) * (z1 - z2)
    return np.where(cold.any(axis=1), height, np.nan)


def half_sunlit(window):
    """The mean of the Sun's cosine, over the cosine of the latitude, in a mean of the given width (radians of the
    Sun's hour angle) that the Sun is up for half of: with its declination zero, a mean centred on sunrise or sunset.
    A mean the Sun is up for more of has a larger one, and every mean wholly in the day one larger still."""
    return (1.0 - np.cos(window / 2.0)) / window


def sun_window(configuration):
    """The Sun's travel (radians of its hour angle) through one output mean: the mean's steps at the Sun's synodic
    rate, from the run's own rotation and year. The model's 1,430 steps of a lunar day make 29.79 days, longer than
    the Sun's 29.53, so a tenth of them takes the Sun through a little more than 36 degrees."""
    m, p = configuration['model'], configuration['planet']
    steps, rest = divmod(m['steps_per_lunar_day'], m['writes_per_lunar_day'])
    if rest:
        raise ValueError('the output means do not hold whole steps')
    days = steps * m['timestep_min'] / 1440.0
    return 2.0 * np.pi * days * (1.0 / p['rotationperiod'] - 1.0 / p['year'])


def ranked_tenths(x, others, k):
    """Each cell's warmest and coolest k means ranked by x, with every series in others averaged over the same
    means, cell by cell."""
    order = np.argsort(x, axis=0, kind='stable')
    return {end: {name: np.take_along_axis(v, idx, 0).mean(axis=0) for name, v in others.items()}
            for end, idx in (('warmest', order[-k:]), ('coolest', order[:k]))}


def heat_in(hot, wl, land, sel=None):
    """The land's heat at both heights in the means of its cells (land, times by cells) that a selection of them
    picks, weighted by land area, with the selection's share of the land's time."""
    ww = np.broadcast_to(wl, land.shape)
    m = land if sel is None else land & sel
    out = dict(share_of_land_time=float((ww * m).sum() / (ww * land).sum()))
    for height, (tw, air) in hot.items():
        x, a, w = tw[m], air[m], ww[m]
        out[height] = dict(wet_bulb=wstats(x, w), air=wstats(a, w), above_28=share(x, w, 28.0),
                           above_31=share(x, w, 31.0))
    return out


def land_and_water(masks):
    """The model's own land mask from each year's lsm, which must hold through the years."""
    lsm = np.concatenate(masks) > 0.5
    if not (lsm == lsm[0]).all():
        raise ValueError('the model\'s land mask changes from mean to mean')
    return lsm[0]


def read(run, b, configuration):
    """Means and distributions over the run's years 10-19."""
    import netCDF4
    files, toa = {}, {}
    fields = {k: [] for k in ('ta', 'hus', 'ps', 'ts', 'tas', 'pr', 'czen', 'lsm')}
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
    island = land_and_water(fields.pop('lsm'))
    sea = ~island
    ta, q = np.concatenate(fields['ta']), np.concatenate(fields['hus'])
    ps, ts, pr = np.concatenate(fields['ps']) * 100.0, np.concatenate(fields['ts']), np.concatenate(fields['pr'])
    tas, czen = np.concatenate(fields['tas']), np.concatenate(fields['czen'])
    del fields
    sst = ts - 273.15
    alat = np.abs(lat)[:, None] * np.ones(land.shape)
    coast = island & beside(sea)
    wl, ws = land * area, (1.0 - land) * area          # each cell's land and water area
    on_land_all = np.broadcast_to(island, sst.shape)
    lw_all = np.broadcast_to(wl, sst.shape)[on_land_all]
    # The heat. Each record's humidity is the humidity at its mean's first step (window_humidity), so a mean's own is
    # taken from its start and its end, the next record's; the last mean has no next record and leaves the heat's
    # statistics, which cover the other means.
    qm = window_humidity(q[:, -1])
    n = qm.shape[0]
    t0, p0, tas_h, ps_h = ta[:n, -1], ps[:n] * sigma[-1], tas[:n], ps[:n]
    # wet-bulb and air temperatures (C) at the two heights in every 3-day mean
    hot = dict(lowest_level=(wet_bulb(t0, qm, p0) - 273.15, t0 - 273.15),
               at_2_m=(wet_bulb(tas_h, np.minimum(qm, qsat(tas_h, ps_h)), ps_h) - 273.15, tas_h - 273.15))
    on_land = np.broadcast_to(island, t0.shape)
    lw = np.broadcast_to(wl, t0.shape)[on_land]         # the heat's land means' weights
    # rain and deserts on the land, counted two ways: every cell's land with that cell's own rain, and the model's own
    # land cells alone; the ground's temperature on the model's land
    rain = pr.mean(axis=0) * 86400.0 * 1000.0
    dry = rain < DESERT_MM_DAY
    ground = sst.mean(axis=0)

    def land_share(weights, mask, with_ground):
        ww = weights * mask
        out = dict(desert=float((ww * dry).sum() / ww.sum()), rain_mm_day=float((ww * rain).sum() / ww.sum()))
        if with_ground:
            out['ground_c'] = float((ww * ground).sum() / ww.sum())
        return out

    def elevations(lo, hi):
        return ((elevation >= lo) if lo is not None else True) & ((elevation < hi) if hi is not None else True)

    def rain_and_deserts(weights, with_ground=False):
        return dict(all=land_share(weights, np.ones_like(island), with_ground),
                    bands={f'{lo}-{hi}': land_share(weights, (alat >= lo) & (alat < hi), with_ground)
                           for lo, hi in RAIN_BANDS},
                    heights={name: land_share(weights, elevations(lo, hi), with_ground)
                             for name, lo, hi in HEIGHT_CLASSES})

    # the day's 3-day means, the Sun up for more than half of each, and among them the afternoon's, the Sun past its
    # highest: the next mean's czen below the last one's
    window = sun_window(configuration)
    day_all = czen > half_sunlit(window) * np.cos(np.radians(lat))[:, None]
    falling = np.empty_like(day_all)
    falling[1:-1] = czen[2:] < czen[:-2]
    falling[0], falling[-1] = czen[1] < czen[0], czen[-1] < czen[-2]
    day, afternoon = day_all[:n], (day_all & falling)[:n]
    # each cell's warmest and coolest tenth of its 3-day means, ranked by the wet-bulb temperature and by the air, with
    # the other quantities over the same means
    k = max(1, n // 10)
    sets = {}
    for height, (tw, air) in hot.items():
        others = dict(wet_bulb=tw, air=air, ground=sst[:n], day=day.astype(float), afternoon=afternoon.astype(float))
        sets[height] = dict(by_wet_bulb=ranked_tenths(tw, others, k), by_air=ranked_tenths(air, others, k))

    def heat_at(height, mask):
        tw, air = hot[height]
        m = np.broadcast_to(mask, tw.shape)
        x, a, ww = tw[m], air[m], np.broadcast_to(wl, tw.shape)[m]
        cells = wl * mask
        mean = lambda f: float((cells * f).sum() / cells.sum())
        return dict(wet_bulb=wstats(x, ww), above_28=share(x, ww, 28.0), above_31=share(x, ww, 31.0),
                    air=wstats(a, ww),
                    tenths={by: {end: {name: mean(v) for name, v in s.items()} for end, s in t.items()}
                            for by, t in sets[height].items()})

    def heat(mask):
        out = heat_at('lowest_level', mask)
        out['at_2_m'] = heat_at('at_2_m', mask)
        return out

    def seas(mask):
        m = np.broadcast_to(sea & mask, sst.shape)
        return wstats(sst[m], np.broadcast_to(ws, sst.shape)[m])

    land_heat = heat(island)
    land_heat['means'] = dict(of_the_run=int(sst.shape[0]), in_the_heat=n)
    land_heat['day'], land_heat['night'] = heat_in(hot, wl, on_land, day), heat_in(hot, wl, on_land, ~day)
    land_heat['day']['czen_above'] = dict(times_cos_latitude=float(half_sunlit(window)),
                                          window_degrees=float(np.degrees(window)))
    # the first readout of 17 October, which took the output's humidity as written and counted every cell with a land
    # fraction above one half alike
    x = (wet_bulb(ta[:, -1], q[:, -1], ps * sigma[-1]) - 273.15)[np.broadcast_to(land > 0.5, sst.shape)]
    land_heat['readout_of_17_october'] = dict(
        wet_bulb=dict(median=float(np.median(x)), p99=float(np.percentile(x, 99))),
        above_28=float((x > 28.0).mean()), above_31=float((x > 31.0).mean()))
    # what the 2 m figures are: the lowest level's air down the dry adiabat, the saturation caps, the ground's pressure
    adiabat = (tas_h - t0 * sigma[-1] ** -AKAP)[on_land]
    land_heat['at_2_m']['air_down_the_dry_adiabat'] = dict(
        mean_k=wmean(adiabat, lw), sd_k=float(np.sqrt(wmean((adiabat - wmean(adiabat, lw)) ** 2, lw))),
        above_the_lowest_level_k=wmean((tas_h - t0)[on_land], lw))
    land_heat['saturated'] = dict(
        capped_at_2_m=share((qm - qsat(tas_h, ps_h))[on_land], lw, 0.0),
        supersaturated_lowest_level=share((qm - qsat(t0, p0))[on_land], lw, 0.0),
        supersaturated_lowest_level_as_written=share((q[:, -1] - qsat(ta[:, -1], ps * sigma[-1]))[on_land_all],
                                                     lw_all, 0.0))
    land_heat['surface_pressure_bar'] = {f'p{round(100 * p):02d}': wquantile(ps_h[on_land] / 1e5, lw, p)
                                         for p in (0.05, 0.5, 0.95)}
    in_water = wl * sea                                  # the boundary file's land in cells the model runs as water
    out = dict(top_of_air_w_m2=dict(incoming=t['rst'] - t['rsut'], reflected=-t['rsut'], absorbed=t['rst'],
                                     absorbed_clear_sky=t['rstcs'], outgoing_infrared=-t['rlut'],
                                     outgoing_infrared_clear_sky=-t['rltcs'], net=t['ntr']),
               ground_k=gmean(ts.mean(axis=0)), rain_mm_day=gmean(rain),
               land_and_water=dict(land_cells=int(island.sum()),
                                   water_share_of_area=float(area[sea].sum() / area.sum()),
                                   water_cells_with_land_fraction_above_half=int((sea & (land > 0.5)).sum()),
                                   land_in_water_cells=float(in_water.sum() / wl.sum()),
                                   its_rain_mm_day=float((in_water * rain).sum() / in_water.sum())),
               land=dict(every_cell=rain_and_deserts(wl), model_land=rain_and_deserts(wl * island, True)),
               heat=dict(land=land_heat, coasts=heat(coast),
                         bands={f'{lo}-{hi}': heat(island & (alat >= lo) & (alat < hi)) for lo, hi in BANDS},
                         below_1_km=heat(island & (elevation < 1000.0)),
                         from_4_km=heat(island & (elevation >= 4000.0))),
               seas=dict(all=seas(np.ones_like(sea)),
                         bands={f'{lo}-{hi}': seas((alat >= lo) & (alat < hi)) for lo, hi in BANDS}))
    # the sky over land: each column's air by height above its ground, means by land area, hypsometric heights
    z = np.zeros_like(ta)
    z[:, -1] = RD * ta[:, -1] / GRAVITY * np.log(1.0 / sigma[-1])
    for lev in range(sigma.size - 2, -1, -1):
        z[:, lev] = (z[:, lev + 1]
                     + RD * 0.5 * (ta[:, lev] + ta[:, lev + 1]) / GRAVITY * np.log(sigma[lev + 1] / sigma[lev]))
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
    # the median of the land's column crossings of freezing, over the columns that freeze at all
    fh = freezing_height(ta, z)[on_land_all]
    freezes = np.isfinite(fh)
    out['freezing_level_over_land_km'] = wquantile(fh[freezes], lw_all[freezes], 0.5) / 1000.0
    out['land_columns_that_never_freeze'] = float(lw_all[~freezes].sum() / lw_all.sum())
    return out, files, island


def record(run):
    """A run's own yearly record and settings (progress.json) and its sha256."""
    path = RUNS_DIR / run / 'progress.json'
    return json.loads(path.read_text()), {str(path.relative_to(ROOT)): sha256(path)}


def progress(p):
    return [dict(year=y['year'], ground_k=y['surface_k'], top_of_air_net_w_m2=y['toa_net_w_m2']) for y in p['years']]


def differs(a, b):
    """The settings in which two runs differ: the model's by name, and every other part of their configuration but
    the run key and the state each starts from."""
    ma, mb = a['configuration']['model'], b['configuration']['model']
    ca, cb = a['configuration'], b['configuration']
    rest = [k for k in set(ca) | set(cb) if k not in ('model', 'key', 'initial_state') and ca.get(k) != cb.get(k)]
    return sorted([k for k in set(ma) | set(mb) if ma.get(k) != mb.get(k)] + rest)


def start(records, run, files):
    """That a run starts from the design case's year-19 restart, with that file's sha256, which the record keeps."""
    state = records[run]['configuration']['initial_state']
    digest = sha256(ROOT / START)
    if state['path'] != START or not digest.startswith(state['sha256']):
        raise ValueError(f'{run} does not start from {START}')
    files[START] = digest
    return dict(path=START, sha256=digest)


# What a glow branch changes from the design case: its sunlight, the CO2 that carries the gases (with the air's gases
# and the sunlight bands the CO2 sets), and the infrared added at the top of the air.
GLOW_SETTINGS = ['co2_ppm', 'gas_bands', 'gases_bar', 'sunlight_scale', 'top_infrared_w_m2']


def glow_pair(records, files):
    """The glow pair's ground and settings over years 10-19 beside the design case's, from each run's yearly record.
    The branches start from the design case's year-19 restart and differ from it in their sunlight, their CO2 and the
    infrared added at the top of the air alone."""
    design = records[DESIGN]
    out = {}
    for run in (DESIGN,) + GLOW_PAIR:
        p = records[run]
        c = p['configuration']
        if run != DESIGN:
            if differs(design, p) != GLOW_SETTINGS:
                raise ValueError(f'{run} is not a branch of {DESIGN} in its sunlight, gases and glow alone')
            start(records, run, files)
        rows = [y for y in p['years'] if YEARS[0] <= y['year'] <= YEARS[1]]
        if len(rows) != YEARS[1] - YEARS[0] + 1:
            raise ValueError(f'{run}: years {YEARS[0]}-{YEARS[1]} are not all in its record')
        ground = [y['surface_k'] for y in rows]
        out[run] = dict(sunlight_scale=c['model']['sunlight_scale'], co2_ppm=c.get('co2_ppm'),
                        top_infrared_w_m2=c.get('top_infrared_w_m2'), key=c['key'], ground_k=float(np.mean(ground)),
                        ground_trend_k_per_decade=float(10.0 * np.polyfit([y['year'] for y in rows], ground, 1)[0]),
                        rain_mm_day=float(np.mean([y['precipitation_mm_day'] for y in rows])))
    return out


def eddies(b, records, island):
    """The design case's years written every 6.5 hours: the land's heat at both heights in those means and in the
    3-day means their elevens make, over the same years, on the design case's own land."""
    import netCDF4
    p, design = records[EDDIES], records[DESIGN]
    if differs(design, p) != ['writes_per_lunar_day']:
        raise ValueError(f'{EDDIES} differs from {DESIGN} in more than how often it writes')
    per, rest = divmod(p['configuration']['model']['writes_per_lunar_day'],
                       design['configuration']['model']['writes_per_lunar_day'])
    if rest:
        raise ValueError(f'{EDDIES}: its means do not fall into whole 3-day means')
    files = {}
    began = start(records, EDDIES, files)
    f = {k: [] for k in ('tas', 'ta', 'hus', 'ps', 'lsm')}
    for year in range(EDDY_YEARS[0], EDDY_YEARS[1] + 1):
        path = output(EDDIES, year)
        files[str(path.relative_to(ROOT))] = sha256(path)
        with netCDF4.Dataset(path) as d:
            lat = np.asarray(d['lat'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            if d['tas'].shape[0] % per:
                raise ValueError(f'{path}: {d["tas"].shape[0]} means do not make whole 3-day means')
            f['tas'].append(np.asarray(d['tas'][:], float))
            f['ta'].append(np.asarray(d['ta'][:, -1], float))           # the lowest level alone
            f['hus'].append(np.asarray(d['hus'][:, -1], float))
            f['ps'].append(np.asarray(d['ps'][:], float) * 100.0)
            f['lsm'].append(np.asarray(d['lsm'][:], float))
    if not (land_and_water(f.pop('lsm')) == island).all():
        raise ValueError(f'{EDDIES} runs on another land than {DESIGN}')
    f = {key: np.concatenate(v) for key, v in f.items()}
    order = [int(np.argmin(np.abs(b['lat'] - x))) for x in lat]
    wl = b['land'][order] * b['area'][order]

    def heat(tas, ta, hus, ps):
        hot = dict(lowest_level=(wet_bulb(ta, hus, ps * sigma[-1]) - 273.15, ta - 273.15),
                   at_2_m=(wet_bulb(tas, np.minimum(hus, qsat(tas, ps)), ps) - 273.15, tas - 273.15))
        land = np.broadcast_to(island, tas.shape)
        out = heat_in(hot, wl, land)
        del out['share_of_land_time']
        out['lowest_level']['supersaturated'] = share((hus - qsat(ta, ps * sigma[-1]))[land],
                                                      np.broadcast_to(wl, tas.shape)[land], 0.0)
        return out

    # each 6.5-hour mean's humidity from its start and its end (window_humidity); the last record has no next one,
    # so the whole 3-day means the rest make are kept, and the 6.5-hour means are the same ones
    written = f['hus']
    blocks = (written.shape[0] - 1) // per
    keep = blocks * per
    fine = dict(tas=f['tas'][:keep], ta=f['ta'][:keep], hus=window_humidity(written)[:keep], ps=f['ps'][:keep])
    m = p['configuration']['model']
    return dict(run=EDDIES, years=list(EDDY_YEARS), start=began,
                writes_per_lunar_day=m['writes_per_lunar_day'], means_to_a_3_day_mean=per,
                means=int(written.shape[0]), means_kept=keep, three_day_means_kept=blocks,
                hours_each=m['steps_per_lunar_day'] * m['timestep_min'] / 60.0 / m['writes_per_lunar_day'],
                at_6_5_hours=heat(**fine),
                in_3_day_means=heat(**{k: v.reshape(-1, per, *v.shape[1:]).mean(axis=1)
                                       for k, v in fine.items()})), files


# ------------------------------------------------------------------------------------- the summit tower's air
# The summit's place and height come from the summit tower study's product; the column is read at the T21 cells
# round it, on heights above the model's sea level, and the runner maps the tower onto it.
SUMMIT_TOWER = ROOT / 'research' / 'studies' / 'summit_tower' / 'results' / 'summit_tower.json'
OLD_DESIGN, OLD_YEARS = 'A28_dim5_moon', (20, 29)    # the earlier design case, whose air the metropolis brief reads
TEST_YEARS = (0, 19)                                  # the test's every year, for the line its own warming draws
LADDER_M = tuple(10000.0 + 250.0 * i for i in range(121))     # heights above sea level, 10 to 40 km
NEAR = 2                                              # the block read: the holding cell and two cells either side
RLAPSE = 0.0065      # the lapse rate pyburn extrapolates the sea-level pressure with (pyburn.py, line 217)


def ground_geopotential(ps_pa, psl_pa, t_lowest, sigma_lowest):
    """The model's own surface geopotential (m2/s2) from its sea-level pressure. The output carries no orography;
    pyburn writes the sea-level pressure from it, extrapolating the lowest level's air to the ground at the standard
    atmosphere's 6.5 K/km and down to sea level along that lapse rate (pyburn.py, lines 1933-1948):
    ln(psl/ps) = a - c a^2/2 + c^2 a^3/3, with a the geopotential over R T* and c = R 0.0065/g. Solved for a by
    Newton's method, it returns the geopotential pyburn used, the same in every mean."""
    c = RD * RLAPSE / GRAVITY
    tstar = (1.0 + c * (1.0 / sigma_lowest - 1.0)) * t_lowest
    tstar = np.where(tstar < 255.0, 0.5 * (255.0 + tstar), tstar)
    target = np.log(psl_pa / ps_pa)
    a = target.copy()
    for _ in range(60):
        a = a - (a - c * a ** 2 / 2.0 + c ** 2 * a ** 3 / 3.0 - target) / (1.0 - c * a + c ** 2 * a ** 2)
    return a * RD * tstar


def geometric(phi):
    """Height (m) of a geopotential under inverse-square gravity; the GCM's own heights take gravity constant."""
    return MOON_RADIUS * phi / (GRAVITY * MOON_RADIUS - phi)


def geopotential(z):
    return GRAVITY * MOON_RADIUS * z / (MOON_RADIUS + z)


def cell_edges(lat):
    """The Gaussian rows' edges (degrees), from the cumulative Gaussian weights, north to south as the rows run."""
    _, w = np.polynomial.legendre.leggauss(lat.size)
    edges = np.degrees(np.arcsin(np.clip(-1.0 + np.concatenate([[0.0], np.cumsum(w)]), -1.0, 1.0)))
    return edges[::-1] if lat[0] > lat[-1] else edges


def site_weights(lat, lon, site_lat, site_lon):
    """The cell holding the site (row, column), and three weightings over the whole grid: bilinear between the four
    cell centres round the site, the holding cell alone, and the holding cell with its eight neighbours by area."""
    edges = cell_edges(lat)
    lo, hi = np.minimum(edges[:-1], edges[1:]), np.maximum(edges[:-1], edges[1:])
    row = int(np.where((site_lat >= lo) & (site_lat < hi))[0][0])
    step = 360.0 / lon.size
    col = int(np.floor(((site_lon - lon[0]) % 360.0) / step + 0.5)) % lon.size
    order = np.argsort(lat)
    la = lat[order]
    k = int(np.clip(np.searchsorted(la, site_lat) - 1, 0, la.size - 2))
    wy = (site_lat - la[k]) / (la[k + 1] - la[k])
    x = ((site_lon - lon[0]) % 360.0) / step
    j = int(np.floor(x)) % lon.size
    wx = x - np.floor(x)
    bilinear = np.zeros((lat.size, lon.size))
    for r, fy in ((order[k], 1.0 - wy), (order[k + 1], wy)):
        for cc, fx in ((j, 1.0 - wx), ((j + 1) % lon.size, wx)):
            bilinear[r, cc] += fy * fx
    holding = np.zeros_like(bilinear)
    holding[row, col] = 1.0
    block = np.zeros_like(bilinear)
    for r in (row - 1, row, row + 1):
        for cc in (col - 1, col, col + 1):
            block[r, cc % lon.size] = np.cos(np.radians(lat[r]))
    return (row, col), dict(bilinear=bilinear, holding=holding, block=block / block.sum())


def level_geopotential(t, q, ps_pa, phis, sigma):
    """Each level's geopotential (m2/s2) over the model's own ground, hypsometric with the virtual temperature: the
    lowest level at its own air above the ground, each step up at the mean of the two levels'. The levels run from the
    top down; t and q are (means, levels, rows, columns)."""
    tv = t * (1.0 + q * (R_VAPOUR / RD - 1.0))
    phi = np.empty_like(t)
    phi[:, -1] = phis + RD * tv[:, -1] * np.log(1.0 / sigma[-1])
    for lev in range(sigma.size - 2, -1, -1):
        phi[:, lev] = phi[:, lev + 1] + RD * 0.5 * (tv[:, lev] + tv[:, lev + 1]) * np.log(sigma[lev + 1] / sigma[lev])
    return phi


def on_ladder(t, phi, ps_pa, sigma, weights):
    """The air (K) and its pressure (Pa) at the ladder's heights above sea level, each cell's column interpolated
    linearly in geopotential between its levels, then weighted over the cells: (means, heights)."""
    target = geopotential(np.asarray(LADDER_M))
    cells = [(r, c, w) for (r, c), w in np.ndenumerate(weights) if w > 0.0]
    air = np.zeros((t.shape[0], target.size))
    pressure = np.zeros_like(air)
    lnsig = np.log(sigma[::-1])
    for r, c, w in cells:
        for m in range(t.shape[0]):
            p = phi[m, ::-1, r, c]                       # from the lowest level up
            if target[0] < p[0]:
                raise ValueError('the ladder starts below a column\'s lowest level')
            air[m] += w * np.interp(target, p, t[m, ::-1, r, c])
            pressure[m] += w * np.exp(np.interp(target, p, lnsig)) * ps_pa[m, r, c]
    return air, pressure


def summit_block(run, years, rows, cols, files):
    """The block's air, humidity, surface and sea-level pressure from a run's yearly output, its files hashed."""
    import netCDF4
    f = {k: [] for k in ('ta', 'hus', 'ps', 'psl')}
    for year in range(years[0], years[1] + 1):
        path = output(run, year)
        key = str(path.relative_to(ROOT))
        if key not in files:
            files[key] = sha256(path)
        with netCDF4.Dataset(path) as d:
            sigma = np.asarray(d['lev'][:], float)
            f['ta'].append(np.asarray(d['ta'][:, :, rows, cols], float))
            f['hus'].append(np.asarray(d['hus'][:, :, rows, cols], float))
            f['ps'].append(np.asarray(d['ps'][:, rows, cols], float) * 100.0)
            f['psl'].append(np.asarray(d['psl'][:, rows, cols], float))     # written in Pa, labelled hPa
    return {k: np.concatenate(v) for k, v in f.items()}, sigma


def ladder_stats(x):
    """A series of means at the ladder's heights: the mean, the warmest and coolest tenth ranked at each height, the
    10th and 90th percentiles and the extremes."""
    order = np.sort(x, axis=0)
    k = max(1, x.shape[0] // 10)
    r = lambda v: [round(float(a), 4) for a in v]
    return dict(mean=r(x.mean(axis=0)), warm_tenth=r(order[-k:].mean(axis=0)), cold_tenth=r(order[:k].mean(axis=0)),
                p10=r(np.percentile(x, 10, axis=0)), p90=r(np.percentile(x, 90, axis=0)), low=r(x.min(axis=0)),
                high=r(x.max(axis=0)), means=int(x.shape[0]))


def summit(records, files):
    """The summit's column in the design case and the test over years 10-19, the test's every year, the earlier
    design case's years 20-29 and the design case's 6.5-hour years: the model's own ground round the summit, and the
    air at heights above sea level in each 3-day mean, under three weightings of the cells round it."""
    import netCDF4
    site = json.loads(SUMMIT_TOWER.read_text())['site']
    files[str(SUMMIT_TOWER.relative_to(ROOT))] = sha256(SUMMIT_TOWER)
    with netCDF4.Dataset(output(DESIGN, YEARS[0])) as d:
        lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
    (row, col), weights = site_weights(lat, lon, site['lat_deg'], site['lon_deg'])
    rows, cols = slice(row - NEAR, row + NEAR + 1), np.arange(col - NEAR, col + NEAR + 1) % lon.size
    local = {k: w[rows][:, cols] for k, w in weights.items()}
    if any(abs(w.sum() - 1.0) > 1e-12 for w in local.values()):
        raise ValueError('a weighting reaches beyond the block read')
    b = boundary()
    order = [int(np.argmin(np.abs(b['lat'] - x))) for x in lat]
    out = dict(site=dict(lat_deg=site['lat_deg'], lon_deg=site['lon_deg'], ground_above_sea_m=site['ground_above_sea_m']),
               holding_cell=dict(row=row, column=col, lat_deg=float(lat[row]), lon_deg=float(lon[col])),
               block=dict(lat_deg=[float(v) for v in lat[rows]], lon_deg=[float(v) for v in lon[cols]],
                          weights={k: [[round(float(v), 6) for v in r] for r in w] for k, w in local.items()},
                          boundary_land_fraction=[[round(float(v), 4) for v in r] for r in b['land'][order][rows][:, cols]],
                          boundary_land_elevation_m=[[round(float(v), 1) for v in r]
                                                     for r in b['elevation'][order][rows][:, cols]]),
               ladder_m=list(LADDER_M), runs={})

    def column(run, years, per=1):
        f, sigma = summit_block(run, years, rows, cols, files)
        phis = ground_geopotential(f['ps'], f['psl'], f['ta'][:, -1], sigma[-1])
        spread = float(np.abs(phis - phis.mean(axis=0)).max() / GRAVITY)
        phis = phis.mean(axis=0)
        # each mean's humidity from its start and its end (window_humidity); the last mean has no next record
        q = window_humidity(f['hus'])
        t, ps = f['ta'][:-1], f['ps'][:-1]
        if per > 1:                                       # the 3-day means the 6.5-hour means' elevens make
            n = (t.shape[0] // per) * per
            t, q, ps = (v[:n].reshape(-1, per, *v.shape[1:]).mean(axis=1) for v in (t, q, ps))
        phi = level_geopotential(t, q, ps, phis, sigma)
        ground = geometric(phis)
        levels = dict(sigma=[float(v) for v in sigma],
                      height_m=[round(float((local['bilinear'] * geometric(phi[:, k]).mean(axis=0)).sum()), 1)
                                for k in range(sigma.size)],
                      air_c=[round(float((local['bilinear'] * t[:, k].mean(axis=0)).sum()) - 273.15, 3)
                             for k in range(sigma.size)])
        ladders, pressure = {}, None
        for name, w in local.items():
            air, p = on_ladder(t, phi, ps, sigma, w)
            ladders[name] = air - 273.15
            if name == 'bilinear':
                pressure = p
        return dict(ground_m=[[round(float(v), 1) for v in r] for r in ground], ground_spread_m=spread,
                    ground_under_m={k: float((w * ground).sum()) for k, w in local.items()},
                    levels=levels, pressure_atm=[round(float(v), 5) for v in pressure.mean(axis=0) / 101325.0],
                    air_c={k: ladder_stats(v) for k, v in ladders.items()}), ladders['bilinear']

    out['runs'][DESIGN], _ = column(DESIGN, YEARS)
    whole, series = column(TEST, TEST_YEARS)
    # years 10-19 of the test's twenty, with each mean's humidity from the next record as in a read of those years
    per_year = records[TEST]['configuration']['model']['writes_per_lunar_day'] * \
        records[TEST]['configuration']['model']['lunar_days_per_year']
    start = (YEARS[0] - TEST_YEARS[0]) * per_year
    out['runs'][TEST], _ = column(TEST, YEARS)
    out['test_years'] = [dict(year=TEST_YEARS[0] + i, air_c=[round(float(v), 4) for v in
                                                            series[i * per_year:(i + 1) * per_year].mean(axis=0)])
                         for i in range(TEST_YEARS[1] - TEST_YEARS[0] + 1)]
    if not np.allclose(series[start:].mean(axis=0), out['runs'][TEST]['air_c']['bilinear']['mean'], atol=1e-3):
        raise ValueError('the test\'s years 10-19 differ between the two reads')
    out['runs'][OLD_DESIGN], _ = column(OLD_DESIGN, OLD_YEARS)
    p = records[EDDIES]['configuration']['model']['writes_per_lunar_day'] // \
        records[DESIGN]['configuration']['model']['writes_per_lunar_day']
    out['eddies'] = dict(at_6_5_hours=column(EDDIES, EDDY_YEARS)[0], in_3_day_means=column(EDDIES, EDDY_YEARS, p)[0])
    return out


def extract():
    b = boundary()
    files = {str(BOUNDARY.relative_to(ROOT)): sha256(BOUNDARY)}
    records = {}
    for run in (DESIGN, TEST) + GLOW_PAIR + (EDDIES, OLD_DESIGN):
        records[run], f = record(run)
        files.update(f)
    runs, lands = {}, {}
    for run in (DESIGN, TEST):
        runs[run], f, lands[run] = read(run, b, records[run]['configuration'])
        files.update(f)
    if not (lands[DESIGN] == lands[TEST]).all():
        raise ValueError(f'{DESIGN} and {TEST} run on different land')
    ed, f = eddies(b, records, lands[DESIGN])
    files.update(f)
    pair = glow_pair(records, files)
    top = summit(records, files)
    return dict(schema=SCHEMA, years=list(YEARS), runs=runs, test_progress=progress(records[TEST]),
                glow_pair=pair, eddies=ed, summit=top,
                producer=dict(files={str(Path(__file__).resolve().relative_to(ROOT)): sha256(__file__)}, read=files,
                              sizes={p: (ROOT / p).stat().st_size for p in files}))


def main(argv=None) -> int:
    out = extract()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
