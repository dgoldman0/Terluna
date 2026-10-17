"""The heat budget of a CM1 ring's air by height: what warms and cools it over its second lunar day, how its mean
temperature drifts from the GCM's profile it starts from, and the heat the GCM's circulation carries sideways at the
ring's latitude, which the ring does not get; python -m climate.crm.ring_heat_budget ring_45n_corrected.

From the ring's snapshots, as ring means: radiation (swten + lwten, temperature tendencies), the imposed large-scale
sinking (-w dtheta/dz on each column's potential temperature, with the vertical wind of terluna_wls.txt, as CM1's
wsub applies it), the nudging of the mean potential temperature toward the GCM's reference (lsnudge_0001.dat) over
lsnudge_tau, weighted from var14 to var15 m, and the storage; the rest is CM1's own transport, turbulence and
convection, which carry the ground's heat up. From the GCM's 3-day means at the ring's row and the rows either side:
the zonal-mean meridional circulation's sideways advection, -[v] d[T]/dy, and the convergence of the eddies' heat
flux, -1/(a cos) d([v'T'] cos)/dphi, as far as 3-day means resolve them. Layers are the GCM's lowest four over its sea,
the ring's flat ground standing at sea level.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-heat-budget/1'
RD, CP = 287.04, 1005.7
KAPPA = RD / CP
DAY = 86400.0


def namelist_number(text: str, key: str) -> float:
    m = re.search(rf'^\s*{key}\s*=\s*([-+0-9.eEdD]+)', text, re.M)
    if not m:
        raise KeyError(key)
    return float(m.group(1).replace('d', 'e').replace('D', 'e'))


def read_reference(path: Path):
    """Heights (m) and potential temperatures (K) of CM1's nudging reference."""
    rows = [line.split() for line in path.read_text().splitlines()[4:] if line.strip()]
    z, th = (np.array([float(r[i]) for r in rows]) for i in (0, 1))
    return z, th


def read_wls(path: Path):
    """Heights (m) and vertical wind (m/s, up positive) of terluna_wls.txt."""
    lines = path.read_text().split('\n')
    n = int(lines[0])
    values = np.array([[float(v) for v in lines[1 + k].split()] for k in range(n)])
    return values[:, 0], values[:, 1]


def nudge_weight(z, bottom, top):
    """CM1's nudging weight: nothing below `bottom`, rising linearly to full at `top`."""
    return np.clip((np.asarray(z, float) - bottom) / max(top - bottom, 1e-9), 0.0, 1.0)


def layer_integral(values, rho, zw, z0, z1):
    """The column integral (per m2) of rho cp times a profile on scalar levels between heights z0 and z1, each level
    weighted by its own depth between w levels."""
    dz = np.diff(zw)
    zh = 0.5 * (zw[1:] + zw[:-1])
    inside = (zh >= z0) & (zh < z1)
    return float((np.asarray(values)[inside] * np.asarray(rho)[inside] * CP * dz[inside]).sum())


def layer_mean(values, rho, zw, z0, z1):
    """The mass-weighted mean of a profile on scalar levels between heights z0 and z1."""
    dz = np.diff(zw)
    zh = 0.5 * (zw[1:] + zw[:-1])
    inside = (zh >= z0) & (zh < z1)
    w = np.asarray(rho)[inside] * dz[inside]
    return float((np.asarray(values)[inside] * w).sum() / w.sum())


def ring_budget(name, from_day):
    """The ring's mean profiles and heating rates (K/day of temperature) over the span after from_day, and the
    drift of its mean temperature through the whole run."""
    from climate.crm import ring_analysis as ra
    case = ra.RUNS / name
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    tap = geo['record']['configuration']['output_s']
    text = (case / 'namelist.input').read_text()
    tau = namelist_number(text, 'lsnudge_tau')
    weight = nudge_weight(zh, namelist_number(text, 'var14'), namelist_number(text, 'var15'))
    zr, thr = read_reference(case / 'lsnudge_0001.dat')
    th_ref = np.interp(zh, zr, thr)
    zl, wl = read_wls(case / 'terluna_wls.txt')
    w_ls = np.interp(zh, zl, wl)
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    names = ('t', 'theta', 'radiation', 'sinking', 'nudging', 'rho')
    span = {k: [] for k in names}
    times, drift = [], []
    surface = []
    for n in outputs:
        t = (n - 1) * tap
        d = ra.read_snapshot(case, n)
        exner = (d['prs'] / 1.0e5) ** KAPPA
        temp = d['th'] * exner
        rho = d['prs'] / (RD * temp * (1.0 + 0.61 * d['qv']))
        mean_t = temp.mean(axis=1)
        times.append(t / DAY)
        drift.append(mean_t)
        if t < from_day * DAY:
            continue
        dthdz = np.gradient(d['th'], zh, axis=0)
        mean_th = d['th'].mean(axis=1)
        span['t'].append(mean_t)
        span['theta'].append(mean_th)
        span['radiation'].append((d['swten'] + d['lwten']).mean(axis=1) * DAY)
        span['sinking'].append((-w_ls[:, None] * dthdz * exner).mean(axis=1) * DAY)
        span['nudging'].append(-weight * (mean_th - th_ref) / tau * exner.mean(axis=1) * DAY)
        span['rho'].append(rho.mean(axis=1))
        surface.append(float((d['thflux'] * ra.air_density(d)).mean() * CP))
    prof = {k: np.mean(v, axis=0) for k, v in span.items()}
    span_t = np.array(span['t'])
    span_days = np.array([x for x in times if x >= from_day])
    storage = np.polyfit(span_days, span_t, 1)[0]                     # K/day at each height
    prof['storage'] = storage
    prof['rest'] = storage - prof['radiation'] - prof['sinking'] - prof['nudging']
    return dict(zh=zh, zw=zw, profiles=prof, surface_heat_w_m2=float(np.mean(surface)), times_day=np.array(times),
                drift_t=np.array(drift), snapshots=len(span_days), w_ls=w_ls, tau_s=tau,
                nudge_m=(namelist_number(text, 'var14'), namelist_number(text, 'var15')))


def block_means(a, block):
    """Means over consecutive blocks of `block` outputs along the first axis, a trailing part block left out."""
    if block <= 1:
        return a
    n = a.shape[0] // block * block
    return a[:n].reshape((n // block, block) + a.shape[1:]).mean(axis=1)


def zonal_eddy_sums(v, t, q, block=1):
    """Sums over time of the zonal means of the northward wind, temperature and vapour (time, level, row, longitude)
    and of their zonal-eddy covariances, the deviations from each time's zonal mean; with block > 1 the outputs are
    first averaged in blocks, as coarser means would hold them. Returns the sums and the number of times."""
    v, t, q = (block_means(a, block) for a in (v, t, q))
    vm, tm, qm = (a.mean(axis=3, keepdims=True) for a in (v, t, q))
    eddy = ((v - vm) * (t - tm)).mean(axis=3)                             # (time, lev, rows)
    eddy_q = ((v - vm) * (q - qm)).mean(axis=3)
    return dict(v=vm[..., 0].sum(axis=0), t=tm[..., 0].sum(axis=0), q=qm[..., 0].sum(axis=0),
                eddy=eddy.sum(axis=0), eddy_q=eddy_q.sum(axis=0)), v.shape[0]


def gcm_sideways(run, years, latitude, band, block=1):
    """The GCM's sideways heat transport at the ring's row (K/day of temperature, by level, the top level first): the
    zonal-mean meridional circulation's advection and the resolved eddies' heat flux convergence, from the run's
    means (3-day means for the runner's standard output) at the row and the rows either side, with the zonal-mean
    temperature and northward wind. block > 1 first averages the outputs in blocks of that many."""
    import netCDF4
    from climate.crm.cm1_run import planet
    radius = planet()['radius_m']
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    sums = None
    count = 0
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            lat = np.asarray(d['lat'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            j = int(np.argmin(np.abs(lat - latitude)))
            rows = [j - 1, j, j + 1] if lat[j - 1] > lat[j + 1] else [j + 1, j, j - 1]   # north to south
            v = np.asarray(d['va'][:, :, rows, :], float)                  # (time, lev, 3, lon)
            t = np.asarray(d['ta'][:, :, rows, :], float)
            q = np.asarray(d['hus'][:, :, rows, :], float)
        part, n = zonal_eddy_sums(v, t, q, block)
        sums = part if sums is None else {k: sums[k] + part[k] for k in sums}
        count += n
    mean = {k: v / count for k, v in sums.items()}                      # (lev, 3): rows north, middle, south
    phi = np.radians(lat[rows])
    dphi = phi[0] - phi[2]                                               # north minus south
    dtdy = (mean['t'][:, 0] - mean['t'][:, 2]) / (radius * dphi)
    mmc = -mean['v'][:, 1] * dtdy * DAY
    flux = mean['eddy'] * np.cos(phi)[None, :]
    eddies = -(flux[:, 0] - flux[:, 2]) / (radius * dphi * np.cos(phi[1])) * DAY
    dqdy = (mean['q'][:, 0] - mean['q'][:, 2]) / (radius * dphi)
    flux_q = mean['eddy_q'] * np.cos(phi)[None, :]
    vapour = dict(mean_meridional_g_kg_day=(-mean['v'][:, 1] * dqdy * DAY * 1000.0).tolist(),
                  eddies_g_kg_day=(-(flux_q[:, 0] - flux_q[:, 2]) / (radius * dphi * np.cos(phi[1])) * DAY
                                   * 1000.0).tolist())
    return dict(rows_deg=[float(x) for x in lat[rows]], sigma=sigma.tolist(),
                mean_meridional_k_day=mmc.tolist(), eddies_k_day=eddies.tolist(), vapour=vapour,
                northward_wind_m_s=mean['v'][:, 1].tolist(),
                temperature_k=mean['t'][:, 1].tolist(), eddy_flux_k_m_s=mean['eddy'][:, 1].tolist())


def analyse(name, from_day):
    from climate.crm import gcm_profiles as gp
    from climate.crm import ring_analysis as ra
    record = json.loads((ra.RUNS / name / 'case.json').read_text())
    run, years = record['gcm']['run'], record['gcm']['years']
    latitude, band = record['latitude_deg'], record['reference']['band_deg']
    ring = ring_budget(name, from_day)
    gcm_composite = gp.gcm_composite(run, years, latitude, band)
    spans = np.asarray(gcm_composite['water']['heights_m']).mean(axis=0)
    sideways = gcm_sideways(run, years, latitude, band)
    p, zw, rho = ring['profiles'], ring['zw'], ring['profiles']['rho']
    terms = ('radiation', 'sinking', 'nudging', 'rest', 'storage')
    layers = []
    for m in range(gp.LAYERS):
        z0, z1 = float(spans[m]), float(spans[m + 1])
        level = len(sideways['sigma']) - 1 - m                             # the GCM's level of this layer
        row = dict(layer=m + 1, km=[round(z0 / 1000.0, 2), round(z1 / 1000.0, 2)],
                   ring_k_day={k: layer_mean(p[k], rho, zw, z0, z1) for k in terms},
                   ring_w_m2={k: layer_integral(p[k] / DAY, rho, zw, z0, z1) for k in terms},
                   gcm_sideways_k_day=dict(mean_meridional=sideways['mean_meridional_k_day'][level],
                                           eddies=sideways['eddies_k_day'][level]),
                   gcm_sideways_vapour_g_kg_day=dict(
                       mean_meridional=sideways['vapour']['mean_meridional_g_kg_day'][level],
                       eddies=sideways['vapour']['eddies_g_kg_day'][level]),
                   ring_temperature_c=layer_mean(p['t'], rho, zw, z0, z1) - 273.15,
                   gcm_temperature_c=float(np.asarray(gcm_composite['water']['layers_c']).mean(axis=0)[m]))
        layers.append(row)
    below = dict(km=[0.0, 8.0], ring_w_m2={k: layer_integral(p[k] / DAY, rho, zw, 0.0, 8000.0) for k in terms})
    t = ring['times_day']
    drift = []
    for m in range(gp.LAYERS):
        z0, z1 = float(spans[m]), float(spans[m + 1])
        series = np.array([layer_mean(profile, rho, zw, z0, z1) for profile in ring['drift_t']]) - 273.15
        late = t >= from_day
        drift.append(dict(layer=m + 1, start_c=float(series[0]), at_from_day_c=float(np.interp(from_day, t, series)),
                          end_c=float(series[-1]), span_trend_k_day=float(np.polyfit(t[late], series[late], 1)[0])))
    return dict(schema=SCHEMA, case=name, from_day=from_day, snapshots=ring['snapshots'],
                gcm=dict(run=run, years=years),
                nudging=dict(tau_days=ring['tau_s'] / DAY, from_m=ring['nudge_m'][0], full_m=ring['nudge_m'][1]),
                surface_heat_w_m2=ring['surface_heat_w_m2'], layers=layers, below_8_km=below, drift=drift,
                gcm_sideways=sideways,
                evidence=('Ring means of CM1\'s own tendencies and of the imposed forcing, from its 3-hourly '
                          'snapshots; the rest closes the budget and holds the model\'s transport, turbulence and '
                          'convection. The '
                          'GCM\'s sideways terms come from 3-day means, which miss its faster eddies, so they bound '
                          'the transport from below.'),
                reading_rule=('Heating rates in K/day of temperature, layer means by mass; W/m2 are rho cp times them '
                              'integrated over each layer. Storage is the trend of the ring mean over the span. The '
                              'GCM\'s sideways terms are at the level of each layer, at the ring\'s row.'))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    args = parser.parse_args(argv)
    out = analyse(args.case, args.from_day)
    out['producer'] = dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'gcm_profiles.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f'heat_budget_{args.case}.json').write_text(json.dumps(out, indent=1) + '\n')
    print(f'surface heat {out["surface_heat_w_m2"]:.1f} W/m2; nudging {out["nudging"]}')
    for row in out['layers']:
        r, g = row['ring_k_day'], row['gcm_sideways_k_day']
        print(f'layer {row["layer"]} {row["km"]} km: ring {row["ring_temperature_c"]:.1f} C, '
              f'GCM {row["gcm_temperature_c"]:.1f} C; ' + ', '.join(f'{k} {v:+.3f}' for k, v in r.items())
              + ' K/day; GCM sideways ' + ', '.join(f'{k} {v:+.3f}' for k, v in g.items()))
    print('below 8 km, W/m2:', {k: round(v, 1) for k, v in out['below_8_km']['ring_w_m2'].items()})
    for row in out['layers']:
        print(f'layer {row["layer"]} W/m2', {k: round(v, 1) for k, v in row['ring_w_m2'].items()},
              'GCM sideways vapour g/kg/day', {k: round(v, 4) for k, v in row['gcm_sideways_vapour_g_kg_day'].items()})
    for row in out['drift']:
        print(f'drift layer {row["layer"]}: {row["start_c"]:.2f} -> {row["at_from_day_c"]:.2f} -> '
              f'{row["end_c"]:.2f} C, {row["span_trend_k_day"]:+.4f} K/day over the span')
    return 0


if __name__ == '__main__':
    sys.exit(main())
