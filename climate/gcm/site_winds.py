"""Winds of an ExoPlaSim run at one site, by height, written as a data product for other lanes.

    climate/gcm/.venv/bin/python -m climate.gcm.site_winds A28_dim5:15-24 summit 5.375 201.375
    # -> results/gcm/site_winds_A28_dim5_summit.json

For the chosen model years it interpolates the eastward and northward wind of every model layer bilinearly
to the site, finds each layer's height above the model's flat surface from the site's own temperature and
humidity (the hypsometric equation with the design air's gas constant and inverse-square gravity), and
reports each layer's wind speed distribution, the mean cube of the speed (for wind power), the mean vector
and the steadiness of direction, the air density, the speed by local hour angle and the yearly highest
speeds. The output is a 3-day mean, a tenth of a lunar day: the distribution is that of 3-day winds, whose
spread is smaller than that of instantaneous winds.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-site-winds/1'
PERCENTILES = (10, 50, 90, 99)
HOUR_BINS = 10
EVIDENCE = ('ExoPlaSim 3-day output means at T21 (about 170 km cells) on 10 sigma layers, interpolated bilinearly '
            'to the site. The model\'s ground is its T21 cells\' mean elevation, so on high ground it stands below '
            'the real peaks; heights are above that model ground, whose own height above the model\'s sea level is '
            'given. Averaging over 3 days removes storms, gusts and turbulence; the speeds are those of the '
            'large-scale flow.')
READING_RULE = ('Arrays run over the model layers from the lowest up. height_m is the mean height of the layer '
                'midpoint above the model ground at the site, and ground_height_m that ground\'s height above the '
                'model\'s sea level (from the surface and sea-level pressures); speed_* are the horizontal speed sqrt(u^2 + v^2) of the 3-day '
                'means in m/s; mean_cube_m3_s3 is the mean of the speed cubed and mean_cube_u_m3_s3 that of the eastward component alone; steadiness is |mean vector| / '
                'mean speed (1 for a wind that never turns); density_kg_m3 is the mean air density. by_hour_angle '
                'gives the mean speed and mean cube by local hour angle (degrees east of the subsolar point, 0 = '
                'noon, negative = morning). yearly_max_m_s[year][layer] is each model year\'s highest speed.')


def bilinear(lat, lon, site_lat, site_lon):
    """Indices and weights for bilinear interpolation to a site on a grid with rows in either order and
    periodic longitudes (degrees east)."""
    lat = np.asarray(lat, dtype=float)
    order = np.argsort(lat)
    la = lat[order]
    k = int(np.clip(np.searchsorted(la, site_lat) - 1, 0, la.size - 2))
    wy = float(np.clip((site_lat - la[k]) / (la[k + 1] - la[k]), 0.0, 1.0))
    rows = (int(order[k]), int(order[k + 1]))
    lon = np.asarray(lon, dtype=float)
    step = 360.0 / lon.size
    x = ((site_lon - lon[0]) % 360.0) / step
    j = int(np.floor(x)) % lon.size
    wx = float(x - np.floor(x))
    cols = (j, (j + 1) % lon.size)
    return rows, cols, (1 - wy, wy), (1 - wx, wx)


def at_site(field, weights):
    """A field (..., lat, lon) interpolated to the site."""
    rows, cols, wy, wx = weights
    return sum(wy[a] * wx[b] * field[..., rows[a], cols[b]] for a in (0, 1) for b in (0, 1))


def layer_heights(sigma, ps_pa, t_k, q, gas_constant, gm, radius):
    """Heights (m) of the layer midpoints above the surface for each sample (time, layer), bottom layer first.

    sigma, t_k and q run from the lowest layer up. The geopotential rises by R T_v ln(p_below/p_above), with the
    surface at the lowest layer's virtual temperature and each step between midpoints at the mean of the two;
    heights follow from inverse-square gravity, z = GM / (GM/R - phi) - R.
    """
    eps = 0.01801528 / (8.314462618 / gas_constant)
    tv = np.asarray(t_k) * (1.0 + np.asarray(q) * (1.0 / eps - 1.0))
    p = np.asarray(sigma)[None, :] * np.asarray(ps_pa)[:, None]
    phi = np.empty_like(tv)
    phi[:, 0] = gas_constant * tv[:, 0] * np.log(np.asarray(ps_pa) / p[:, 0])
    for k in range(1, tv.shape[1]):
        phi[:, k] = phi[:, k - 1] + gas_constant * 0.5 * (tv[:, k - 1] + tv[:, k]) * np.log(p[:, k - 1] / p[:, k])
    return gm / (gm / radius - phi) - radius, p / (gas_constant * tv)


def gumbel_return(maxima, years):
    """Return level for a return period (years) from yearly maxima, by the method of moments."""
    m = np.asarray(maxima, dtype=float)
    beta = np.sqrt(6.0) * m.std(ddof=1) / np.pi
    mu = m.mean() - 0.5772156649 * beta
    return float(mu - beta * np.log(-np.log(1.0 - 1.0 / years)))


def read(folder: str, first: int, last: int, site_lat: float, site_lon: float):
    import netCDF4
    from climate.gcm.climatology import subsolar_longitudes
    from climate.gcm.exoplasim_run import RUNS
    out = {k: [] for k in ('u', 'v', 't', 'q', 'ps', 'psl', 'hour', 'year')}
    weights = None
    for year in range(first, last + 1):
        with netCDF4.Dataset(RUNS / folder / 'model' / f'MOST.{year:05d}.nc') as d:
            lat = np.asarray(d['lat'][:], dtype=float)
            lon = np.asarray(d['lon'][:], dtype=float)
            sigma = np.asarray(d['lev'][:], dtype=float)
            weights = weights or bilinear(lat, lon, site_lat, site_lon)
            rows = sorted(weights[0])
            sub = slice(rows[0], rows[1] + 1)
            local = (tuple(r - rows[0] for r in weights[0]), weights[1], weights[2], weights[3])
            for key, name in (('u', 'ua'), ('v', 'va'), ('t', 'ta'), ('q', 'hus')):
                out[key].append(at_site(np.asarray(d[name][:, :, sub, :], dtype=float), local))
            out['ps'].append(at_site(np.asarray(d['ps'][:, sub, :], dtype=float), local) * 100.0)
            # PlaSim labels psl hPa but writes Pa (about 1.2e5 at sea level); ps is in hPa as labelled.
            out['psl'].append(at_site(np.asarray(d['psl'][:, sub, :], dtype=float), local))
            band = np.abs(lat) < 20.0
            czen = np.asarray(d['czen'][:, band, :], dtype=float)
            solar = subsolar_longitudes(czen, lat[band], lon)
            out['hour'].append((site_lon - solar + 180.0) % 360.0 - 180.0)
            out['year'].append(np.full(czen.shape[0], year))
    data = {k: np.concatenate(v) for k, v in out.items()}
    bottom_up = np.argsort(sigma)[::-1]
    for k in ('u', 'v', 't', 'q'):
        data[k] = data[k][:, bottom_up]
    return data, sigma[bottom_up]


def ground_height(ps_pa, psl_pa, t_k, q, gas_constant, gm, radius):
    """Mean height (m) of the model ground above the model's sea level: the hypsometric thickness between the
    sea-level and surface pressures at the lowest layer's virtual temperature, with inverse-square gravity."""
    eps = 0.01801528 / (8.314462618 / gas_constant)
    tv = np.asarray(t_k) * (1.0 + np.asarray(q) * (1.0 / eps - 1.0))
    phi = gas_constant * tv * np.log(np.asarray(psl_pa) / np.asarray(ps_pa))
    return float((gm / (gm / radius - phi) - radius).mean())


def site_winds(data: dict, sigma, gas_constant, gm, radius) -> dict:
    height, density = layer_heights(sigma, data['ps'], data['t'], data['q'], gas_constant, gm, radius)
    ground = ground_height(data['ps'], data['psl'], data['t'][:, 0], data['q'][:, 0], gas_constant, gm, radius)
    speed = np.hypot(data['u'], data['v'])
    mean_u, mean_v = data['u'].mean(axis=0), data['v'].mean(axis=0)
    edges = np.linspace(-180.0, 180.0, HOUR_BINS + 1)
    index = np.clip(np.digitize(data['hour'], edges) - 1, 0, HOUR_BINS - 1)
    by_speed = [speed[index == b].mean(axis=0) if (index == b).any() else np.full(sigma.size, np.nan)
                for b in range(HOUR_BINS)]
    by_cube = [(speed[index == b] ** 3).mean(axis=0) if (index == b).any() else np.full(sigma.size, np.nan)
               for b in range(HOUR_BINS)]
    years = np.unique(data['year'])
    yearly = np.array([speed[data['year'] == y].max(axis=0) for y in years])
    r = lambda a, n=3: np.round(np.asarray(a, dtype=float), n).tolist()
    return dict(
        ground_height_m=round(ground, 0), sigma=r(sigma, 4), height_m=r(height.mean(axis=0), 0), pressure_pa=r((sigma[None] * data['ps'][:, None]).mean(axis=0), 0),
        temperature_k=r(data['t'].mean(axis=0), 2), density_kg_m3=r(density.mean(axis=0), 4),
        mean_u_m_s=r(mean_u), mean_v_m_s=r(mean_v), speed_mean_m_s=r(speed.mean(axis=0)),
        **{f'speed_p{p}_m_s': r(np.percentile(speed, p, axis=0)) for p in PERCENTILES},
        speed_max_m_s=r(speed.max(axis=0)), mean_cube_m3_s3=r((speed ** 3).mean(axis=0), 2),
        mean_cube_u_m3_s3=r((np.abs(data['u']) ** 3).mean(axis=0), 2),
        steadiness=r(np.hypot(mean_u, mean_v) / speed.mean(axis=0)),
        by_hour_angle=dict(hour_angle_deg=r(0.5 * (edges[:-1] + edges[1:]), 1),
                           speed_mean_m_s=[r(s) for s in by_speed], mean_cube_m3_s3=[r(c, 2) for c in by_cube]),
        yearly_max_m_s=[r(y) for y in yearly],
        return_50yr_m_s=r([gumbel_return(yearly[:, k], 50.0) for k in range(sigma.size)]),
        samples=int(speed.shape[0]))


def main(argv=None) -> int:
    args = list(argv or sys.argv[1:])
    if len(args) != 4 or ':' not in args[0]:
        raise SystemExit('Give a run and year range, a site name, its latitude and its longitude east, for example '
                         'A28_dim5:15-24 summit 5.375 201.375')
    folder, span = args[0].split(':')
    first, last = (int(v) for v in span.split('-'))
    name, site_lat, site_lon = args[1], float(args[2]), float(args[3]) % 360.0
    from atmosphere.radiative_convective.thermodynamics import MOON, earthlike_air
    from climate.gcm.exoplasim_run import RUNS
    progress = json.loads((RUNS / folder / 'progress.json').read_text())
    air = earthlike_air(progress['configuration']['pressure_pa'], 400.0)
    data, sigma = read(folder, first, last, site_lat, site_lon)
    winds = site_winds(data, sigma, air.gas_constant, MOON.gm_m3_s2, MOON.radius_m)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    product = dict(schema=SCHEMA, run=folder, years=[first, last], site=dict(name=name, lat_deg=site_lat, lon_deg=site_lon),
                   configuration=progress['configuration'],
                   producer=dict(domain='climate', files={'gcm/site_winds.py': digest(__file__)}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, **winds)
    RESULTS.mkdir(parents=True, exist_ok=True)
    out = RESULTS / f'site_winds_{folder}_{name}.json'
    out.write_text(json.dumps(product, indent=1) + '\n')
    for k, z in enumerate(winds['height_m']):
        print(f"{z / 1000:6.1f} km  median {winds['speed_p50_m_s'][k]:5.2f}  p99 {winds['speed_p99_m_s'][k]:5.2f}  "
              f"max {winds['speed_max_m_s'][k]:5.2f}  mean vector ({winds['mean_u_m_s'][k]:+.2f}, {winds['mean_v_m_s'][k]:+.2f}) m/s  "
              f"density {winds['density_kg_m3'][k]:.3f}")
    print(out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
