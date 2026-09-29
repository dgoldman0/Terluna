"""Sunlight in the air at a CM1 box's site: CM1's radiation (RRTMG) and the GCM's against the repository's line-by-line
model (atmosphere/radiative_convective), each on its own clear column.

    climate/gcm/.venv/bin/python -m climate.crm.radiation_check box_0e_small --from-day 29.5

writes ../results/crm/radiation_check_<case>.json. Needs the line-by-line model's spectroscopic inputs
(python -m atmosphere.radiative_convective.fetch_inputs --download) and takes about a quarter of an hour on 8 cores.

Two comparisons. At two Sun heights, the shares of the arriving sunlight that CM1's cloud-free columns reflect to
space, absorb in the air and pass to the ground, against line-by-line on their mean profile. Over the lunar day, the
same three parts in W/m2, clear sky, for the GCM's mean column on its land cells at the site and for CM1's mean
cloud-free daytime column, each against line-by-line: the daily mean weights each Sun height by the time the Sun
spends there (six Gauss-Legendre nodes in hour angle from noon to sunset). CM1's sunlight has RRTMG's own spectrum
scaled to the shield's total and the GCM's comes through the titania stack, so line-by-line is given both ways.
"""
from __future__ import annotations
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from climate.crm import ring_analysis as ra

RUNS, RESULTS = ra.RUNS, ra.RESULTS
SCHEMA = 'terluna.climate.crm-radiation-check/1'
SURFACE_ALBEDO = 0.20                   # the placeholder land's and the GCM's land albedo
SUN_HEIGHTS = {'near_noon': (0.85, 1.0), 'lower_sun': (0.45, 0.60)}   # cosine of the Sun's zenith angle
KAPPA = 287.04 / 1004.64
EVIDENCE = ('Clear-sky solar fluxes only, with no cloud. CM1\'s fluxes are those RRTMG gave its cloud-free columns in '
            'the box; the GCM\'s are its own clear-sky diagnostics (3-day means). The line-by-line model (HITRAN lines, '
            'MT_CKD water continuum, collision-induced absorption, Rayleigh scattering, a delta-scaled two-stream solver '
            'with a pseudo-spherical direct beam) runs on mean profiles, so it leaves out how the fluxes vary with the '
            'profile from hour to hour; the GCM\'s column is its ten levels, extended isothermal and dry above its top.')
READING_RULE = ('Shares are of the sunlight arriving at the top of the air, for a ground albedo of 0.20. Daily means are '
                'W/m2 over the whole lunar day at the site, of the sunlight each model has arriving there: 373.5 W/m2 in '
                'CM1, and in the GCM its own time mean, which its uneven Sun makes 12% less at 0 E. Line-by-line shares '
                'are applied to each model\'s own arriving sunlight. unfiltered is sunlight with '
                'the Sun\'s own spectrum, as CM1\'s RRTMG takes it; titania is sunlight through the titania stack, as the '
                'GCM takes it; both are scaled to the shield\'s total of 1,173.5 W/m2.')


def line_by_line_column(p, t, qv, surface_pa, surface_k):
    """A line-by-line column from levels above the ground (any order) and the surface state; returns it and its water
    vapour (kg/m2)."""
    from atmosphere.radiative_convective import thermodynamics as th
    eps = th.MOLAR_MASS['H2O'] / 0.028964
    order = np.argsort(-np.asarray(p, float))
    p = np.concatenate([[surface_pa], np.asarray(p, float)[order]])
    t = np.concatenate([[surface_k], np.asarray(t, float)[order]])
    q = np.asarray(qv, float)[order]
    q = np.concatenate([[q[0]], q])
    x = (q / eps) / (1 + q / eps)
    air = th.earthlike_air(surface_pa * (1 - x[0]), 400.0)
    col = th.column_from_levels(th.MOON, air, p, t, x)
    water = float(np.sum(col['layer_mass_kg_m2'] * col['layer_x_h2o'] * th.MOLAR_MASS['H2O']
                         / (th.MOLAR_MASS['H2O'] * col['layer_x_h2o'] + air.molar_mass * (1 - col['layer_x_h2o']))))
    return col, water


def solar_parts(col, cos_zenith, processes=8, surface_albedo=SURFACE_ALBEDO, shields=('unfiltered', 'titania')):
    """For each sunlight (unfiltered, titania) and each Sun height: the shares of the arriving sunlight reflected,
    absorbed in the air and absorbed by the ground, clear sky."""
    from atmosphere.radiative_convective import optics, shortwave as sw, climate as cl
    radius = (col['planet'].radius_m + col['z_m'])[::-1]
    nu, tau_abs = optics.optical_depth(col, replace(optics.SHORTWAVE, processes=processes))
    grids = [(nu, tau_abs), (cl.UV_GRID.nu, cl._uv_optics(col, cl.UV_GRID.nu))]
    rayleigh = [optics.rayleigh_optical_depth(col, g) for g, _ in grids]
    out = {}
    for label in shields:
        spectra = [sw.solar_spectrum(g, shield=cl.load_shield('titania_stack') if label == 'titania' else None)
                   for g, _ in grids]
        rows = []
        for mu0 in np.atleast_1d(cos_zenith):
            incident = reflected = ground = 0.0
            for (g, tau), tau_r, (spectrum, longward) in zip(grids, rayleigh, spectra):
                w = np.full(g.size, g[1] - g[0]); w[0] *= 0.5; w[-1] *= 0.5
                direct, diffuse, up = sw.column_fluxes(tau[::-1], tau_r[::-1], np.zeros_like(tau_r[::-1]),
                                                       surface_albedo, float(mu0), radius)
                f = spectrum * w
                incident += mu0 * float(spectrum @ w)
                reflected += float(up[0] @ f)
                ground += float((direct[-1] + diffuse[-1] - up[-1]) @ f)
            incident += mu0 * spectra[0][1]                  # sunlight longward of the grid: absorbed in the air
            rows.append(dict(reflected=reflected / incident, ground=ground / incident,
                             air=(incident - reflected - ground) / incident))
        out[label] = rows
    return out


def daily_mean(parts, nodes, weights, sunlight_w_m2, lat_deg=0.0):
    """W/m2 over the lunar day at a latitude (the Sun on the equator) from shares at Sun heights cos(lat) cos(h)
    for hour angles h at the nodes."""
    mu0 = np.cos(np.radians(lat_deg)) * np.cos(nodes)
    return {part: float(sum(w * m * row[part] for w, m, row in zip(weights, mu0, parts)) * sunlight_w_m2 / np.pi)
            for part in ('reflected', 'air', 'ground')}


def cm1_clear(case: Path, from_day: float, lo: float, hi: float):
    """CM1's cloud-free columns in snapshots whose Sun height lies in [lo, hi] (all daytime ones for lo = 0): the mean
    profile and the mean shares of the arriving sunlight, and over the whole span, the clear-sky daily mean in W/m2."""
    tap = json.loads((case / 'case.json').read_text())['configuration']['output_s']
    acc = dict(p=0.0, t=0.0, qv=0.0, ps=0.0, t2=0.0, reflected=0.0, air=0.0, ground=0.0, cos_zenith=0.0)
    flux = dict(reflected=0.0, air=0.0, ground=0.0)
    used = total = 0
    for n in sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat')):
        if (n - 1) * tap < from_day * 86400.0:
            continue
        total += 1
        d = ra.read_snapshot(case, n)
        mu0 = float(d['coszen'].mean())
        clear = ~((d['qc'] + d['qi'] + d['qs'] + d['qr'] + d['qg']) >= 1e-6).any(axis=0)
        if not (lo < mu0 <= hi) or clear.sum() < 50:
            continue
        top, up = d['swdnt'][clear], d['swupt'][clear]
        ground = d['swdnb'][clear] - d['swupb'][clear]
        for key, value in (('reflected', up), ('air', top - up - ground), ('ground', ground)):
            acc[key] += float((value / top).mean())
            flux[key] += float(value.mean())
        acc['p'] += d['prs'][:, clear].mean(axis=1)
        acc['qv'] += d['qv'][:, clear].mean(axis=1)
        acc['t'] += (d['th'] * (d['prs'] / 1e5) ** KAPPA)[:, clear].mean(axis=1)
        acc['ps'] += float(d['psfc'][clear].mean())
        acc['t2'] += float(d['t2'][clear].mean())
        acc['cos_zenith'] += mu0
        used += 1
    return {k: v / used for k, v in acc.items()}, {k: v / total for k, v in flux.items()}, used


def gcm_column(lat_deg: float, lon_deg: float, run='A28_dim5', years=(15, 24)):
    """The GCM's mean column on its land cells nearest the site, the sunlight arriving there in the GCM, and the GCM's
    own clear-sky solar fluxes (W/m2)."""
    import netCDF4
    folder = ra.HERE.parent / 'gcm' / 'runs' / run / 'model'
    acc = {k: [] for k in ('ta', 'hus', 'ps', 'tas', 'prw', 'rstcs', 'rsscs', 'rst', 'rsut')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            lsm = np.asarray(d['lsm'][:], float)[0] > 0.5
            rows = np.argsort(np.abs(lat - lat_deg))[:2]
            col = int(np.argmin(np.abs((lon - lon_deg + 180.0) % 360.0 - 180.0)))
            rows = [r for r in rows if lsm[r, col]]
            for k in acc:
                v = np.asarray(d[k][:], float)
                acc[k].append(v[..., rows, col])
            sigma = np.asarray(d['lev'][:], float)
    g = {k: np.concatenate(v).mean(axis=0) for k, v in acc.items()}
    ps = float(g['ps'].mean()) * 100.0
    t, q = g['ta'].mean(axis=-1), g['hus'].mean(axis=-1)
    above = np.exp(np.linspace(np.log(sigma[0] * ps * 0.8), np.log(50.0), 12))
    p = np.concatenate([sigma * ps, above])
    t = np.concatenate([t, np.full(above.size, t[0])])
    q = np.concatenate([q, np.full(above.size, 1e-7)])
    arriving = float((g['rst'] - g['rsut']).mean())                       # rsut is negative, upward
    own = dict(reflected=float(arriving - g['rstcs'].mean()), air=float(g['rstcs'].mean() - g['rsscs'].mean()),
               ground=float(g['rsscs'].mean()))
    return dict(p=p, t=t, qv=q, ps=ps, t2=float(g['tas'].mean()), prw=float(g['prw'].mean()), arriving=arriving,
                cells=[[float(lat[r]), float(lon[col])] for r in rows]), own


def analyse(name: str, from_day: float, processes: int = 8) -> dict:
    from climate.crm.cm1_run import CASES, solar_day_s
    case = RUNS / name
    site = CASES[name]['site']
    record = json.loads((case / 'case.json').read_text())
    sunlight = float(record['build']['air']['sunlight_w_m2'])
    nodes, weights = np.polynomial.legendre.leggauss(6)
    hours, hour_weights = 0.25 * np.pi * (nodes + 1.0), 0.25 * np.pi * weights
    heights = {}
    for label, (lo, hi) in SUN_HEIGHTS.items():
        state, _, used = cm1_clear(case, from_day, lo, hi)
        col, water = line_by_line_column(state['p'], state['t'], state['qv'], state['ps'], state['t2'])
        lbl = solar_parts(col, [state['cos_zenith']], processes)
        heights[label] = dict(snapshots=used, cos_zenith=state['cos_zenith'], water_kg_m2=water,
                              cm1={k: state[k] for k in ('reflected', 'air', 'ground')},
                              line_by_line={k: v[0] for k, v in lbl.items()})
    shares = lambda parts: {k: v / sum(parts.values()) for k, v in parts.items()}
    state, own, used = cm1_clear(case, from_day, 0.0, 1.0)
    col, water = line_by_line_column(state['p'], state['t'], state['qv'], state['ps'], state['t2'])
    lbl = solar_parts(col, np.cos(hours), processes)
    lbl = {k: shares(daily_mean(v, hours, hour_weights, sunlight)) for k, v in lbl.items()}
    arriving = sunlight / np.pi
    cm1_day = dict(snapshots=used, water_kg_m2=water, arriving_w_m2=arriving, own=own, own_shares=shares(own),
                   line_by_line_shares=lbl,
                   line_by_line={k: {p: s * arriving for p, s in v.items()} for k, v in lbl.items()})
    gcm, gcm_own = gcm_column(site['lat_deg'], site['lon_deg'])
    col, water = line_by_line_column(gcm['p'], gcm['t'], gcm['qv'], gcm['ps'], gcm['t2'])
    lbl = solar_parts(col, np.cos(hours), processes)
    lbl = {k: shares(daily_mean(v, hours, hour_weights, sunlight)) for k, v in lbl.items()}
    gcm_day = dict(cells=gcm['cells'], water_kg_m2=water, gcm_prw_kg_m2=gcm['prw'], arriving_w_m2=gcm['arriving'],
                   own=gcm_own, own_shares=shares(gcm_own), line_by_line_shares=lbl,
                   line_by_line={k: {p: s * gcm['arriving'] for p, s in v.items()} for k, v in lbl.items()})
    return dict(schema=SCHEMA, case=name, site=site, evidence=EVIDENCE, reading_rule=READING_RULE,
                surface_albedo=SURFACE_ALBEDO, sunlight_w_m2=sunlight, arriving_daily_w_m2=sunlight / np.pi,
                from_day=from_day, cm1_at_sun_heights=heights, daily_mean_clear_sky=dict(cm1=cm1_day, gcm=gcm_day))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    parser.add_argument('--processes', type=int, default=8)
    args = parser.parse_args(argv)
    result = analyse(args.case, args.from_day, args.processes)
    result['producer'] = dict(domain='climate', files={'crm/radiation_check.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    from climate.crm.ring_crossings import rounded
    path = RESULTS / f'radiation_check_{args.case}.json'
    path.write_text(json.dumps(rounded(result, 4), indent=1) + '\n')
    print(f'{path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
