#!/usr/bin/env python3
"""The GCM's heating by process at a ring's row, and its longwave against the correlated-k model.

    climate/gcm/.venv/bin/python -m climate.crm.gcm_heating A28_corrected_970_heating 0-0 --latitude 45

Reads a run written with the runner's --heating-diagnostics, PlaSim's temperature tendency of each process at every
level (codes 21-29), and averages them over the row nearest the latitude and over the run's years, by level and in
W/m² through each level's mass. Their sum is the storage, near zero in a settled run. On the 3-day means at that row
whose columns hold no cloud, it also sets PlaSim's longwave heating against the atmosphere domain's correlated-k model
on the same columns, as ring_radiation.py does for CM1's. Writes results/crm/gcm_heating_<run>.json.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-gcm-heating/1'
DAY = 86400.0
CP = 1004.64
TERMS = ('dt_diab', 'dt_nudge', 'dt_sfc', 'dt_vdif', 'dt_lw', 'dt_sw', 'dt_lsp', 'dt_conv', 'dt_dca')
OLD_NAMES = {'dt_diab': 'dt_dyn'}   # the first run wrote code 21 under the name the runner first gave it
CLOUDLESS = 0.01               # a 3-day mean column counts as cloud-free when no level's cloud cover reaches this


def gcm_column_levels(sigma, ta, hus, ps, top_pa=1.0, extra=8):
    """Levels for the correlated-k model, ground first, from a GCM column on its full levels (top first): the half
    levels' pressures, temperature and vapour carried to them, then isothermal levels above its top down to top_pa."""
    from climate.crm.gcm_eddies import sigma_depths
    _, half = sigma_depths(np.asarray(sigma))
    p_half = half[::-1] * ps                                              # ground first, ending at 0 at the top
    lp_full, t_full, q_full = np.log(np.asarray(sigma)[::-1] * ps), np.asarray(ta)[::-1], np.asarray(hus)[::-1]
    lp = np.log(np.maximum(p_half[:-1], 1e-9))
    t = np.interp(lp[::-1], lp_full[::-1], t_full[::-1])[::-1]
    q = np.interp(lp[::-1], lp_full[::-1], q_full[::-1])[::-1]
    t[0] = t_full[0] + (lp[0] - lp_full[0]) / (lp_full[1] - lp_full[0]) * (t_full[1] - t_full[0])
    top = np.exp(lp_full[-1])                                             # the top full level's pressure
    above = np.exp(np.linspace(np.log(top), np.log(top_pa), extra + 1))
    p = np.concatenate([p_half[:-1], above])
    return p, np.concatenate([t, np.full(above.size, t_full[-1])]), np.concatenate([q, np.full(above.size, q_full[-1])])


def analyse(run, years, latitude, every=1):
    import netCDF4
    from climate.crm.cm1_run import planet
    from climate.crm.gcm_eddies import sigma_depths
    from climate.crm.ring_radiation import mole_fraction
    from atmosphere.radiative_convective import ck
    from atmosphere.radiative_convective import thermodynamics as td
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    record = json.loads((HERE.parent / 'gcm' / 'runs' / run / 'progress.json').read_text())
    cfg = record['configuration']
    cfg = ast.literal_eval(cfg) if isinstance(cfg, str) else cfg
    gravity = planet()['gravity_m_s2']
    sums, count, ps_sum = None, 0, 0.0
    model = ck.CKLongwave()
    air = td.earthlike_air(float(cfg['pressure_pa']), cfg.get('co2_ppm', 400.0))
    lw_pairs = []
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            lat = np.asarray(d['lat'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            j = int(np.argmin(np.abs(lat - latitude)))
            row = {k: np.asarray(d[k if k in d.variables else OLD_NAMES[k]][:, :, j, :], float)
                   for k in TERMS + ('ta', 'hus', 'cl')}
            ps = np.asarray(d['ps'][:, j, :], float) * 100.0                 # hPa in the output
            ts = np.asarray(d['ts'][:, j, :], float)
        part = {k: row[k].mean(axis=2).sum(axis=0) for k in TERMS}            # (lev,): zonal means summed over time
        sums = part if sums is None else {k: sums[k] + part[k] for k in sums}
        count += row['ta'].shape[0]
        ps_sum += float(ps.mean()) * row['ta'].shape[0]
        clear = np.argwhere(row['cl'].max(axis=1) < CLOUDLESS)                 # (time, lon) of cloud-free columns
        for it, il in clear[::every]:
            hus = row['hus'][it, :, il]
            p, t, q = gcm_column_levels(sigma, row['ta'][it, :, il], hus / (1.0 - hus), ps[it, il])   # mixing ratio
            col = td.column_from_levels(td.MOON, air, p, t, mole_fraction(q, air.molar_mass), constant_gravity=True)
            rate, _ = model.heating_rate(col, CP, t_surface=float(ts[it, il]))
            lw_pairs.append((row['dt_lw'][it, ::-1, il] * DAY, rate[:sigma.size]))   # both ground first
    mean = {k: v / count * DAY for k, v in sums.items()}                         # K/day, top level first
    surface_pa = ps_sum / count
    depth, half = sigma_depths(sigma)
    per_k_day = CP * surface_pa * depth / gravity / DAY
    levels = []
    for k in range(sigma.size)[::-1]:                                             # ground first
        rates = {name[3:]: float(mean[name][k]) for name in TERMS}
        processes = sum(v for n, v in rates.items() if n != 'diab')
        levels.append(dict(sigma_span=[float(half[k + 1]), float(half[k])], w_m2_per_k_day=float(per_k_day[k]),
                           k_day=rates, w_m2={n: v * per_k_day[k] for n, v in rates.items()},
                           processes_sum_k_day=float(processes), diab_less_processes_k_day=float(rates['diab'] - processes)))
    out = dict(schema=SCHEMA, run=run, years=list(years), latitude_deg=float(lat[j]), means=count,
               surface_pa=surface_pa, levels=levels)
    if lw_pairs:
        plasim = np.mean([a for a, _ in lw_pairs], axis=0)
        kd = np.mean([b for _, b in lw_pairs], axis=0)
        out['longwave_cloud_free'] = dict(columns=len(lw_pairs), plasim_k_day=plasim.tolist(), correlated_k_k_day=kd.tolist(),
                                          difference_w_m2=((plasim - kd) * per_k_day[::-1]).tolist())
    out['evidence'] = ('PlaSim\'s own temperature tendency of each process, zonal means at the row over the run\'s 3-day '
                       'means; the correlated-k model follows the line-by-line model and runs on the same cloud-free '
                       '3-day mean columns, so it leaves out how the fluxes vary within the 3 days.')
    out['reading_rule'] = ('K/day of temperature at each GCM level, the ground first, and W/m2 through each level\'s mass; '
                           'positive warms. diab is the physics\' sum as the diabatic step hands it to the dynamics, and '
                           'its difference from the processes\' own sum checks them; sfc is the ground\'s sensible heat '
                           'into the lowest layer, vdif turbulence, lsp large-scale condensation, conv Kuo convection, dca '
                           'dry convective adjustment, nudge the top level\'s nudging. In a settled run the dynamics '
                           'carries away the physics\' sum.')
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('run')
    parser.add_argument('years', help='first-last model years, for example 0-0')
    parser.add_argument('--latitude', type=float, default=45.0)
    parser.add_argument('--every', type=int, default=1, help='take every this many cloud-free columns')
    args = parser.parse_args(argv)
    years = tuple(int(x) for x in args.years.split('-'))
    out = analyse(args.run, years, args.latitude, args.every)
    out['producer'] = dict(domain='climate', files={f'crm/{Path(__file__).name}':
                                                    hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f'gcm_heating_{args.run}.json').write_text(json.dumps(out, indent=1) + '\n')
    print(f'row {out["latitude_deg"]:.2f}, {out["means"]} means, surface {out["surface_pa"]:.0f} Pa')
    for row in out['levels'][:4]:
        print(f'sigma {row["sigma_span"][0]:.3f}-{row["sigma_span"][1]:.3f}: ' +
              ', '.join(f'{k} {v:+.3f}' for k, v in row['k_day'].items())
              + f' K/day; processes {row["processes_sum_k_day"]:+.3f}, diab less them {row["diab_less_processes_k_day"]:+.3f}')
    if 'longwave_cloud_free' in out:
        lw = out['longwave_cloud_free']
        print(f'longwave on {lw["columns"]} cloud-free columns, PlaSim against correlated-k, K/day by level from the '
              'ground: ' + ', '.join(f'{a:+.3f}/{b:+.3f}' for a, b in list(zip(lw['plasim_k_day'], lw['correlated_k_k_day']))[:4]))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
