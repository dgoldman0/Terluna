"""What the study reads from the GCM's own output: the design case and the 20-year test of the longwave refitted by
height with ice optics for the cold cloud, over years 10-19 of each, and the test's year-by-year progress; the glow
pair on the corrected radiation, from its runs' yearly records; and the design case's two years written every
6.5 hours.

    climate/gcm/.venv/bin/python -m research.studies.climate_targets.extract     # about a minute; needs netCDF4

It reads the ExoPlaSim output in climate/gcm/runs (ignored by Git) and the T21 boundary file in climate/gcm/products,
and writes research/runs/climate_targets/extract.json (ignored), with the sha256 and size of every file it read;
run.py binds it.

Land and water are the model's own mask (lsm), for every statistic. The runner turns cells to water, least land
first, until the planned water share is met, so cells with a land fraction a little above one half run as water.
Land cells are weighted by their land area (land fraction times cell area) and water cells by their water area, the
quantiles included. The coasts are land cells with a water cell among their eight neighbours, joined across the
longitude seam and kept apart across the poles. The land's rain and deserts are the land cells' too.

The heat is in the 3-day means, each about a tenth of the lunar day, at two heights. One is the lowest model level,
about 0.9 km above the ground. The other is 2 m, where the output's air is the lowest level's carried down the dry
adiabat to the ground's pressure (PlaSim's mkshfl, with its kappa of 0.286), and the wet-bulb temperature carries the
lowest level's specific humidity down with it, capped at saturation there: the reading of a well-mixed lowest
kilometre, as by day. The wet-bulb temperature is the psychrometric one over water with Bolton's saturation, capped at
the air's own temperature where the air is supersaturated. Each place's warmest and coolest tenth are each cell's own
tenth of its 3-day means ranked by one quantity, with the others averaged over the same means.

The day's 3-day means are those the Sun is up for more than half of. The output's czen is each mean's average cosine
of the Sun's zenith angle, and a mean is the day's when its czen exceeds half_sunlit() of the Sun's travel through the
mean times the cosine of the latitude; the cut is exact where the Sun's declination is zero, and the fixed orbit
swings it through the tilt, 1.54 degrees either way. Of the day's means, those with the Sun past its highest are the
afternoon's. The design case's 6.5-hour years (A28_corrected_970_eddies, two years from the design case's year-19
restart, eleven means to each 3-day mean) give the land's heat in those means and in the 3-day means their elevens
make, over the same two years.
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
SCHEMA = 'terluna.research.climate-targets-extract/4'
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
    t0, q0, p0 = ta[:, -1], q[:, -1], ps * sigma[-1]
    # wet-bulb and air temperatures (C) at the two heights in every 3-day mean
    hot = dict(lowest_level=(wet_bulb(t0, q0, p0) - 273.15, t0 - 273.15),
               at_2_m=(wet_bulb(tas, np.minimum(q0, qsat(tas, ps)), ps) - 273.15, tas - 273.15))
    sst = ts - 273.15
    alat = np.abs(lat)[:, None] * np.ones(land.shape)
    coast = island & beside(sea)
    wl, ws = land * area, (1.0 - land) * area          # each cell's land and water area
    on_land = np.broadcast_to(island, sst.shape)
    lw = np.broadcast_to(wl, sst.shape)[on_land]        # the land means' weights
    # rain, deserts and the ground's temperature on the land, by land area
    rain = pr.mean(axis=0) * 86400.0 * 1000.0
    dry = rain < DESERT_MM_DAY
    ground = sst.mean(axis=0)

    def land_share(mask):
        ww = wl * (island & mask)
        return dict(desert=float((ww * dry).sum() / ww.sum()), rain_mm_day=float((ww * rain).sum() / ww.sum()),
                    ground_c=float((ww * ground).sum() / ww.sum()))

    def elevations(lo, hi):
        return ((elevation >= lo) if lo is not None else True) & ((elevation < hi) if hi is not None else True)

    # the day's 3-day means, the Sun up for more than half of each, and among them the afternoon's, the Sun past its
    # highest: the next mean's czen below the last one's
    window = sun_window(configuration)
    day = czen > half_sunlit(window) * np.cos(np.radians(lat))[:, None]
    falling = np.empty_like(day)
    falling[1:-1] = czen[2:] < czen[:-2]
    falling[0], falling[-1] = czen[1] < czen[0], czen[-1] < czen[-2]
    afternoon = day & falling
    # each cell's warmest and coolest tenth of its 3-day means, ranked by the wet-bulb temperature and by the air, with
    # the other quantities over the same means
    k = max(1, sst.shape[0] // 10)
    sets = {}
    for height, (tw, air) in hot.items():
        others = dict(wet_bulb=tw, air=air, ground=sst, day=day.astype(float), afternoon=afternoon.astype(float))
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
    land_heat['day'], land_heat['night'] = heat_in(hot, wl, on_land, day), heat_in(hot, wl, on_land, ~day)
    land_heat['day']['czen_above'] = dict(times_cos_latitude=float(half_sunlit(window)),
                                          window_degrees=float(np.degrees(window)))
    # the first readout of 17 October, which counted every cell with a land fraction above one half alike
    x = hot['lowest_level'][0][np.broadcast_to(land > 0.5, sst.shape)]
    land_heat['readout_of_17_october'] = dict(
        wet_bulb=dict(median=float(np.median(x)), p99=float(np.percentile(x, 99))),
        above_28=float((x > 28.0).mean()), above_31=float((x > 31.0).mean()))
    # what the 2 m figures are: the lowest level's air down the dry adiabat, the saturation caps, the ground's pressure
    adiabat = (tas - t0 * sigma[-1] ** -AKAP)[on_land]
    land_heat['at_2_m']['air_down_the_dry_adiabat'] = dict(
        mean_k=wmean(adiabat, lw), sd_k=float(np.sqrt(wmean((adiabat - wmean(adiabat, lw)) ** 2, lw))),
        above_the_lowest_level_k=wmean((tas - t0)[on_land], lw))
    land_heat['saturated'] = dict(capped_at_2_m=share((q0 - qsat(tas, ps))[on_land], lw, 0.0),
                                  supersaturated_lowest_level=share((q0 - qsat(t0, p0))[on_land], lw, 0.0))
    land_heat['surface_pressure_bar'] = {f'p{round(100 * p):02d}': wquantile(ps[on_land] / 1e5, lw, p)
                                         for p in (0.05, 0.5, 0.95)}
    out = dict(top_of_air_w_m2=dict(incoming=t['rst'] - t['rsut'], reflected=-t['rsut'], absorbed=t['rst'],
                                     absorbed_clear_sky=t['rstcs'], outgoing_infrared=-t['rlut'],
                                     outgoing_infrared_clear_sky=-t['rltcs'], net=t['ntr']),
               ground_k=gmean(ts.mean(axis=0)), rain_mm_day=gmean(rain),
               land_and_water=dict(land_cells=int(island.sum()),
                                   water_share_of_area=float(area[sea].sum() / area.sum()),
                                   water_cells_with_land_fraction_above_half=int((sea & (land > 0.5)).sum())),
               land=dict(all=land_share(np.ones_like(island)),
                         bands={f'{lo}-{hi}': land_share((alat >= lo) & (alat < hi)) for lo, hi in RAIN_BANDS},
                         heights={name: land_share(elevations(lo, hi)) for name, lo, hi in HEIGHT_CLASSES}),
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
    fh = freezing_height(ta, z)[on_land]
    freezes = np.isfinite(fh)
    out['freezing_level_over_land_km'] = wquantile(fh[freezes], lw[freezes], 0.5) / 1000.0
    out['land_columns_that_never_freeze'] = float(lw[~freezes].sum() / lw.sum())
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

    m = p['configuration']['model']
    return dict(run=EDDIES, years=list(EDDY_YEARS), start=began,
                writes_per_lunar_day=m['writes_per_lunar_day'], means_to_a_3_day_mean=per, means=int(f['tas'].shape[0]),
                hours_each=m['steps_per_lunar_day'] * m['timestep_min'] / 60.0 / m['writes_per_lunar_day'],
                at_6_5_hours=heat(**f),
                in_3_day_means=heat(**{k: v.reshape(-1, per, *v.shape[1:]).mean(axis=1) for k, v in f.items()})), files


def extract():
    b = boundary()
    files = {str(BOUNDARY.relative_to(ROOT)): sha256(BOUNDARY)}
    records = {}
    for run in (DESIGN, TEST) + GLOW_PAIR + (EDDIES,):
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
    return dict(schema=SCHEMA, years=list(YEARS), runs=runs, test_progress=progress(records[TEST]),
                glow_pair=pair, eddies=ed,
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
