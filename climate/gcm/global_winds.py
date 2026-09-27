"""Winds of an ExoPlaSim run over the whole Moon, by height above sea level, written as a data product for other
lanes.

    climate/gcm/.venv/bin/python -m climate.gcm.global_winds A28_dim5:15-24
    # -> results/gcm/global_winds_A28_dim5.json

For the chosen model years it finds, in every cell and output, the height of each model layer above the model's
sea level: the ground's height from the surface and sea-level pressures, and each layer's height above the ground
from the column's own temperature and humidity, as site_winds.py does. It interpolates the horizontal wind speed
and the temperature linearly in height, and the air density and pressure logarithmically, to fixed heights above
sea level. Heights below a cell's ground are left out; between the ground and the lowest layer's midpoint (about 0.9 km up) the
lowest layer's values stand. Everything is weighted by cell area. Speeds and densities are gathered in fine
histograms, which give the percentiles; the highest speed, each model year's highest and a 50-year return level
(Gumbel, as in site_winds.py) are kept at every height, and the speed percentiles by band of latitude. The output
is a 3-day mean at T21, so storms, gusts and turbulence are averaged out.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from climate.gcm.site_winds import RESULTS, gumbel_return, layer_heights

SCHEMA = 'terluna.climate.gcm-global-winds/1'
HEIGHTS_KM = np.arange(0.0, 90.01, 2.5)
LATITUDE_BANDS_DEG = ((0.0, 20.0), (20.0, 40.0), (40.0, 60.0), (60.0, 90.0))
SPEED_BINS = np.arange(0.0, 60.0001, 0.02)
DENSITY_BINS = np.arange(0.0, 2.0001, 0.0005)
PERCENTILES = (10, 50, 90, 99, 99.9)
SEA_TOLERANCE_M = 50.0   # sea cells whose pressures put the ground a few metres off sea level still count at 0 km
EVIDENCE = ('ExoPlaSim 3-day output means at T21 (about 170 km cells) on 10 sigma layers over the whole Moon, '
            'interpolated in height to fixed heights above the model\'s sea level. Averaging over 3 days removes '
            'storms, gusts and turbulence; the speeds are those of the large-scale flow. The model\'s ground is its '
            'cells\' mean elevation, so heights just above high ground are under-sampled.')
READING_RULE = ('Arrays run over height_km, heights above the model\'s sea level. Statistics at a height use only '
                'the area where the ground lies below that height (coverage gives its share of the Moon), weighted '
                'by cell area. speed_* are the horizontal speed sqrt(u^2 + v^2) of the 3-day means in m/s, '
                'density_* the air density in kg/m3, pressure_pa and temperature_k area-weighted means. yearly_max_m_s'
                '[year][height] is each model year\'s highest speed and return_50yr_m_s a Gumbel fit to them. '
                'by_latitude gives the speed percentiles in bands of absolute latitude (both hemispheres).')


def weighted_percentiles(counts, edges, percentiles):
    """Percentiles from a weighted histogram (last axis the bins), interpolating within the bin."""
    counts = np.asarray(counts, dtype=float)
    cum = np.cumsum(counts, axis=-1)
    total = cum[..., -1:]
    out = []
    for p in percentiles:
        target = total * p / 100.0
        k = np.clip((cum < target).sum(axis=-1), 0, counts.shape[-1] - 1)
        below = np.take_along_axis(cum, k[..., None], -1)[..., 0] - np.take_along_axis(counts, k[..., None], -1)[..., 0]
        inside = np.take_along_axis(counts, k[..., None], -1)[..., 0]
        share = np.where(inside > 0, (target[..., 0] - below) / np.where(inside > 0, inside, 1.0), 0.0)
        value = edges[k] + np.clip(share, 0.0, 1.0) * (edges[k + 1] - edges[k])
        out.append(np.where(total[..., 0] > 0, value, np.nan))
    return np.array(out)


def to_heights(z, values, heights_m, log=False):
    """Interpolate each sample's column (samples, layers; z rising along the layers) to fixed heights. Below the
    lowest layer the lowest value stands; above the highest layer the result is NaN."""
    n, layers = z.shape
    out = np.full((heights_m.size, n), np.nan)
    v = np.log(values) if log else values
    for i, h in enumerate(heights_m):
        k = np.clip((z <= h).sum(axis=1) - 1, 0, layers - 2)
        z0 = np.take_along_axis(z, k[:, None], 1)[:, 0]
        z1 = np.take_along_axis(z, (k + 1)[:, None], 1)[:, 0]
        w = np.clip((h - z0) / (z1 - z0), 0.0, 1.0)
        a = np.take_along_axis(v, k[:, None], 1)[:, 0]
        b = np.take_along_axis(v, (k + 1)[:, None], 1)[:, 0]
        r = a + w * (b - a)
        r = np.where(h > z[:, -1], np.nan, r)
        out[i] = np.exp(r) if log else r
    return out


def read_year(path, gas_constant, gm, radius):
    """One model year: speed, density, pressure and temperature at the fixed heights, with the ground mask and the
    latitude of every sample."""
    import netCDF4
    with netCDF4.Dataset(path) as d:
        lat = np.asarray(d['lat'][:], dtype=float)
        sigma = np.asarray(d['lev'][:], dtype=float)
        up = np.argsort(sigma)[::-1]                         # lowest layer first
        field = lambda name: np.moveaxis(np.asarray(d[name][:], dtype=float)[:, up], 1, -1).reshape(-1, sigma.size)
        u, v, t, q = (field(n) for n in ('ua', 'va', 'ta', 'hus'))
        ps = np.asarray(d['ps'][:], dtype=float).reshape(-1) * 100.0      # hPa, as labelled
        psl = np.asarray(d['psl'][:], dtype=float).reshape(-1)            # Pa, although labelled hPa
        times = d['ua'].shape[0]
    sig = sigma[up]
    above_ground, density = layer_heights(sig, ps, t, q, gas_constant, gm, radius)
    eps = 0.01801528 / (8.314462618 / gas_constant)
    tv0 = t[:, 0] * (1.0 + q[:, 0] * (1.0 / eps - 1.0))
    phi = gas_constant * tv0 * np.log(np.maximum(psl, ps) / ps)
    ground = gm / (gm / radius - phi) - radius
    z = above_ground + ground[:, None]
    hm = HEIGHTS_KM * 1e3
    speed = to_heights(z, np.hypot(u, v), hm)
    rho = to_heights(z, density, hm, log=True)
    pressure = to_heights(z, sig[None, :] * ps[:, None], hm, log=True)
    temperature = to_heights(z, t, hm)
    under = hm[:, None] < ground[None, :] - SEA_TOLERANCE_M
    for a in (speed, rho, pressure, temperature):
        a[under] = np.nan
    lat_s = np.broadcast_to(lat[None, :, None], (times, lat.size, ps.size // (times * lat.size))).reshape(-1)
    return speed, rho, pressure, temperature, lat_s


def global_winds(folder: str, first: int, last: int, gas_constant, gm, radius) -> dict:
    from climate.gcm.exoplasim_run import RUNS
    nh = HEIGHTS_KM.size
    speed_hist = np.zeros((nh, SPEED_BINS.size - 1))
    band_hist = np.zeros((len(LATITUDE_BANDS_DEG), nh, SPEED_BINS.size - 1))
    rho_hist = np.zeros((nh, DENSITY_BINS.size - 1))
    sums = {k: np.zeros(nh) for k in ('w', 'p', 't', 'area')}
    yearly, highest, samples = [], np.zeros(nh), 0
    for year in range(first, last + 1):
        speed, rho, pressure, temperature, lat = read_year(RUNS / folder / 'model' / f'MOST.{year:05d}.nc',
                                                           gas_constant, gm, radius)
        weight = np.cos(np.radians(lat))
        samples += lat.size
        top = np.full(nh, 0.0)
        for i in range(nh):
            ok = np.isfinite(speed[i])
            w = weight[ok]
            speed_hist[i] += np.histogram(np.clip(speed[i][ok], 0, SPEED_BINS[-1] - 1e-9), SPEED_BINS, weights=w)[0]
            rho_hist[i] += np.histogram(np.clip(rho[i][ok], 0, DENSITY_BINS[-1] - 1e-9), DENSITY_BINS, weights=w)[0]
            for b, (lo, hi) in enumerate(LATITUDE_BANDS_DEG):
                inb = ok & (np.abs(lat) >= lo) & (np.abs(lat) < hi)
                band_hist[b, i] += np.histogram(np.clip(speed[i][inb], 0, SPEED_BINS[-1] - 1e-9), SPEED_BINS,
                                                weights=weight[inb])[0]
            sums['w'][i] += w.sum()
            sums['area'][i] += w.sum() / weight.sum()
            sums['p'][i] += (pressure[i][ok] * w).sum()
            sums['t'][i] += (temperature[i][ok] * w).sum()
            top[i] = speed[i][ok].max() if ok.any() else np.nan
        yearly.append(top)
        highest = np.fmax(highest, top)
    years = last - first + 1
    yearly = np.array(yearly)
    r = lambda a, n=3: [None if not np.isfinite(x) else round(float(x), n) for x in np.asarray(a, dtype=float)]
    with np.errstate(invalid='ignore', divide='ignore'):
        sp = weighted_percentiles(speed_hist, SPEED_BINS, PERCENTILES)
        dp = weighted_percentiles(rho_hist, DENSITY_BINS, (1, 50, 99))
        mean_rho = (rho_hist * 0.5 * (DENSITY_BINS[1:] + DENSITY_BINS[:-1])).sum(axis=1) / rho_hist.sum(axis=1)
        mean_speed = (speed_hist * 0.5 * (SPEED_BINS[1:] + SPEED_BINS[:-1])).sum(axis=1) / speed_hist.sum(axis=1)
        bands = []
        for b, (lo, hi) in enumerate(LATITUDE_BANDS_DEG):
            bp = weighted_percentiles(band_hist[b], SPEED_BINS, (50, 99))
            bands.append(dict(abs_lat_deg=[lo, hi], speed_p50_m_s=r(bp[0], 2), speed_p99_m_s=r(bp[1], 2)))
    return dict(
        height_km=r(HEIGHTS_KM, 1), coverage=r(sums['area'] / years),
        pressure_pa=r(sums['p'] / sums['w'], 0), temperature_k=r(sums['t'] / sums['w'], 2),
        density_mean_kg_m3=r(mean_rho, 4), density_p1_kg_m3=r(dp[0], 4), density_p50_kg_m3=r(dp[1], 4),
        density_p99_kg_m3=r(dp[2], 4), speed_mean_m_s=r(mean_speed, 2),
        **{f'speed_p{str(p).replace(".", "_")}_m_s': r(sp[k], 2) for k, p in enumerate(PERCENTILES)},
        speed_max_m_s=r(highest, 2), yearly_max_m_s=[r(y, 2) for y in yearly],
        return_50yr_m_s=r([gumbel_return(yearly[:, i], 50.0) if np.all(np.isfinite(yearly[:, i])) else np.nan
                           for i in range(HEIGHTS_KM.size)], 2),
        by_latitude=bands, samples=int(samples))


def main(argv=None) -> int:
    args = list(argv or sys.argv[1:])
    if len(args) != 1 or ':' not in args[0]:
        raise SystemExit('Give a run and year range, for example A28_dim5:15-24')
    folder, span = args[0].split(':')
    first, last = (int(v) for v in span.split('-'))
    from atmosphere.radiative_convective.thermodynamics import MOON, earthlike_air
    from climate.gcm.exoplasim_run import RUNS
    progress = json.loads((RUNS / folder / 'progress.json').read_text())
    air = earthlike_air(progress['configuration']['pressure_pa'], 400.0)
    winds = global_winds(folder, first, last, air.gas_constant, MOON.gm_m3_s2, MOON.radius_m)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    product = dict(schema=SCHEMA, run=folder, years=[first, last], configuration=progress['configuration'],
                   producer=dict(domain='climate', files={'gcm/global_winds.py': digest(__file__),
                                                          'gcm/site_winds.py': digest(Path(__file__).with_name('site_winds.py'))}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, **winds)
    RESULTS.mkdir(parents=True, exist_ok=True)
    out = RESULTS / f'global_winds_{folder}.json'
    out.write_text(json.dumps(product, indent=1) + '\n')
    for k, z in enumerate(winds['height_km']):
        print(f"{z:5.1f} km  cover {winds['coverage'][k]:.2f}  density {winds['density_mean_kg_m3'][k]}  "
              f"median {winds['speed_p50_m_s'][k]}  p99 {winds['speed_p99_m_s'][k]}  max {winds['speed_max_m_s'][k]}  "
              f"50-yr {winds['return_50yr_m_s'][k]}")
    print(out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
