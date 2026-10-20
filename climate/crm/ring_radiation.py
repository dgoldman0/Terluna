#!/usr/bin/env python3
"""CM1's longwave heating against the atmosphere domain's correlated-k model, for a ring's own cloud-free columns.

    climate/gcm/.venv/bin/python -m climate.crm.ring_radiation ring_45n_corrected_sideways --from-day 29.5

The ring's RRTMG carries the design air's CO2 and O2 and no ozone, methane, nitrous oxide or CFCs. Every column of
the 45° N ring holds some ice 60-150 km up, so the columns compared are those cloud-free below 25 km; the correlated-k
model leaves that ice out, which touches only the layers near it. CM1 stores its radiative tendencies at constant
volume (cm1_output.py), so its lwten is first turned back into RRTMG's own heating rate. The correlated-k model,
built from the line-by-line spectroscopy, takes the same air on the same layers from each sampled column's
temperature, vapour, pressure and ground temperature, with the air at the ground at the ground's temperature and
CM1's cp, as CM1 hands them to RRTMG, plane-parallel at the Moon's surface gravity as RRTMG is. Above the ring's top,
which still has about 5% of the air over it, the column continues isothermally to 1 Pa. Writes
results/crm/radiation_<case>.json: both heating rates by height as means over the sampled columns, with CM1's lwten
as written beside them, and their means over the GCM's lowest layers.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from climate.crm import cm1_output as co

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-ring-radiation/2'
DAY = 86400.0
CP = co.CP                     # J/kg/K, CM1's, which it hands RRTMG
LAYERS_KM = ((0.0, 1.77), (1.77, 6.53), (6.53, 13.74), (13.74, 23.06))   # the GCM's lowest four, as the budget has them
CLEAR_BELOW_M = 25000.0        # columns compared are cloud-free below this; the ring's ice 60-150 km up stays in both


def column_levels(zw, zh, prs, temp, qv, psfc, top_pa=1.0, extra=6, t_ground=None):
    """Levels for the correlated-k model, ground first: CM1's interfaces with the pressure, temperature and vapour of
    its layers carried to them, then isothermal levels above the top down to top_pa. With t_ground the air at the
    ground takes the ground's temperature, as CM1 hands it to RRTMG."""
    lp = np.interp(zw, zh, np.log(prs), left=np.nan, right=np.nan)
    t = np.interp(zw, zh, temp, left=np.nan, right=np.nan)
    q = np.interp(zw, zh, qv, left=np.nan, right=np.nan)
    # Below the first layer centre and above the last, carry the nearest two layers' gradient.
    for i, (a, b) in ((0, (0, 1)), (-1, (-1, -2))):
        f = (zw[i] - zh[a]) / (zh[b] - zh[a])
        lp[i] = np.log(prs[a]) + f * (np.log(prs[b]) - np.log(prs[a]))
        t[i] = temp[a] + f * (temp[b] - temp[a])
        q[i] = max(qv[a] + f * (qv[b] - qv[a]), 0.0)
    lp[0] = np.log(psfc)
    if t_ground is not None:
        t[0] = t_ground
    p = np.exp(lp)
    above = np.exp(np.linspace(np.log(p[-1]), np.log(top_pa), extra + 1)[1:])
    return (np.concatenate([p, above]), np.concatenate([t, np.full(extra, t[-1])]),
            np.concatenate([q, np.full(extra, q[-1])]))


def mole_fraction(qv, dry_molar):
    """Water vapour mole fraction from the mixing ratio (kg per kg of dry air)."""
    n = qv / 0.01801528
    return n / (n + 1.0 / dry_molar)


def compare(case_name, from_day, every_snapshot=8, every_column=40):
    from climate.crm import ring_analysis as ra
    from atmosphere.radiative_convective import ck
    from atmosphere.radiative_convective import thermodynamics as td
    case = ra.RUNS / case_name
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    tap = geo['record']['configuration']['output_s']
    equation_set = co.eqtset(case)
    air = td.earthlike_air(121590.0, 400.0)
    model = ck.CKLongwave()
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    span = [n for n in outputs if (n - 1) * tap >= from_day * DAY][::every_snapshot]
    cm1, written, kd, mass, kept, skipped = [], [], [], [], 0, 0
    for n in span:
        d = ra.read_snapshot(case, n)
        temp = co.temperature(d)
        rate = co.rrtmg_rate(d['lwten'], d, equation_set) * DAY
        low = zh < CLEAR_BELOW_M
        for i in range(0, temp.shape[1], every_column):
            if d['cldfra'][low, i].max() > 0.01 or (d['qc'][low, i] + d['qi'][low, i]).max() > 1e-7:
                skipped += 1
                continue
            p, t, q = column_levels(zw, zh, d['prs'][:, i], temp[:, i], d['qv'][:, i], d['psfc'][i],
                                    t_ground=float(d['tsk'][i]))
            col = td.column_from_levels(td.MOON, air, p, t, mole_fraction(q, air.molar_mass), constant_gravity=True)
            heating, _ = model.heating_rate(col, CP, t_surface=float(d['tsk'][i]))
            kd.append(heating[:zh.size])
            cm1.append(rate[:, i])
            written.append(d['lwten'][:, i] * DAY)
            mass.append(col['layer_mass_kg_m2'][:zh.size])
            kept += 1
    cm1, written, kd, mass = (np.mean(a, axis=0) for a in (cm1, written, kd, mass))
    layers = []
    for z0, z1 in LAYERS_KM:
        m = (zh >= z0 * 1000.0) & (zh < z1 * 1000.0)
        w = mass[m]
        layers.append(dict(km=[z0, z1], cm1_rrtmg_k_day=float((cm1[m] * w).sum() / w.sum()),
                           cm1_lwten_as_written_k_day=float((written[m] * w).sum() / w.sum()),
                           correlated_k_k_day=float((kd[m] * w).sum() / w.sum()),
                           w_m2_per_k_day=float(w.sum() * CP / DAY)))
    return dict(schema=SCHEMA, case=case_name, from_day=from_day, eqtset=equation_set, columns=kept,
                cloudy_columns_left_out=skipped, snapshots=len(span), every_column=every_column, z_m=zh.tolist(),
                cm1_rrtmg_k_day=cm1.tolist(), cm1_lwten_as_written_k_day=written.tolist(),
                correlated_k_k_day=kd.tolist(), layers=layers, ck_table_key=ck.table_key(),
                evidence=('Longwave heating of the same cloud-free columns from CM1\'s RRTMG, its stored lwten turned '
                          'back into RRTMG\'s own rate (cm1_output.py), and from the atmosphere domain\'s correlated-k '
                          'model, which follows the line-by-line model to within the accuracy its validation records, '
                          'on the column CM1 hands RRTMG. Means over the sampled columns.'),
                reading_rule=('K/day of temperature at constant pressure at CM1\'s layer centres (z_m), and '
                              'mass-weighted means over the GCM\'s lowest four layers; negative cools. '
                              'cm1_lwten_as_written_k_day is CM1\'s stored tendency at constant volume, cp/cvm times '
                              'RRTMG\'s rate under eqtset 2. W/m2 per K/day gives each layer\'s heat capacity per unit '
                              'area.'))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    parser.add_argument('--every-snapshot', type=int, default=8)
    parser.add_argument('--every-column', type=int, default=40)
    args = parser.parse_args(argv)
    out = compare(args.case, args.from_day, args.every_snapshot, args.every_column)
    out['producer'] = dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'cm1_output.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f'radiation_{args.case}.json').write_text(json.dumps(out, indent=1) + '\n')
    print(f'{out["columns"]} cloud-free columns from {out["snapshots"]} snapshots ({out["cloudy_columns_left_out"]} '
          'cloudy left out)')
    for row in out['layers']:
        print(f'{row["km"][0]:5.2f}-{row["km"][1]:5.2f} km: RRTMG {row["cm1_rrtmg_k_day"]:+.3f} (written '
              f'{row["cm1_lwten_as_written_k_day"]:+.3f}), correlated-k '
              f'{row["correlated_k_k_day"]:+.3f} K/day; difference {row["cm1_rrtmg_k_day"] - row["correlated_k_k_day"]:+.3f} '
              f'K/day = {(row["cm1_rrtmg_k_day"] - row["correlated_k_k_day"]) * row["w_m2_per_k_day"]:+.1f} W/m2')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
