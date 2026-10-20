#!/usr/bin/env python3
"""RRTMG's tables ending at 1,053.63 hPa, and CM1's call of RRTMG, against the correlated-k model in the lunar air.

    climate/gcm/.venv/bin/python -m climate.crm.rrtmg_edge ring_45n_corrected_sideways --from-day 29.5

RRTMG's absorption tables start at 1,053.63 hPa (pref(1) in module_ra_rrtmg_lw.F), and a layer at higher pressure
takes its absorption by extrapolating linearly in ln p from the first two reference pressures. The lunar air stands at
about 1,216 hPa at the water level, so the lowest ~7 km over the seas are extrapolated. Here CM1's own RRTMG
longwave, compiled from the lunar build's objects into a driver of its own (rrtmg_offline/lw_columns.f90), takes the
45° N ring's cloud-free columns three ways:

- as they are, the layers below the edge extrapolated, as CM1 runs them;
- as twins whose pressures are scaled so that the ground stands at 1,040 hPa, every layer inside the tables, with the
  same temperatures and mole fractions;
- as they are, with the extrapolation held at the edge (fp no lower than 0 in setcoef), the absorption coefficients
  of 1,053.63 hPa carried to the ground.

The correlated-k model takes the columns and their twins. If RRTMG's error against it in the lowest layers is the same
on the twins as on the columns, the extrapolation adds little to it; the change between the first and third ways
shows how far the extrapolated layers' heating depends on the extrapolation at all. Every column takes the same layers
in both models (the correlated-k column's), clear sky, a black ground under air at the ground's temperature, CM1's cp
and the design air's CO2 and O2 alone, as the ring's RRTMG carries them. CM1's stored lwten, turned back into RRTMG's
own rate (cm1_output.py), stands beside them as CM1's call. The comparison of 17 October (ring_radiation.py) read it
as stored, a tendency at constant volume, and so found CM1 cooling these columns more than the correlated-k model at
every height. Writes results/crm/rrtmg_edge_<case>.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from climate.crm import cm1_output as co

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-rrtmg-edge/2'
BUILD = Path('/media/projectspace/terluna-research/cm1/build/moon_omp_rain_limits/src')
WORK = Path('/media/projectspace/terluna-research/crm-runs/rrtmg_edge')
DRIVER = HERE / 'rrtmg_offline' / 'lw_columns.f90'
MODULE = 'module_ra_rrtmg_lw.f90'
OPTS = ['-ffree-form', '-ffree-line-length-none', '-O3', '-march=native', '-finline-functions',
        '--param=max-vartrack-size=0', '-fopenmp']                      # the lunar build's own (its Makefile, OPTS)
TABLE_EDGE_HPA = 1053.63
TWIN_GROUND_HPA = 1040.0
# The edge held: setcoef's interpolation fraction in ln p kept from going negative, which only the layers beyond
# the first reference pressure reach.
EDGE_EDIT = ("         fp = 5._rb *(preflog(jp(lay)) - plog)\n",
             "         fp = max(0._rb, 5._rb *(preflog(jp(lay)) - plog))   ! Terluna: the edge held (rrtmg_edge.py)\n")
LAYERS_KM = ((0.0, 1.77), (1.77, 6.53), (6.53, 13.74), (13.74, 23.06))   # the GCM's lowest four (ring_radiation.py)
DAY = 86400.0


def sha16(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def build(variant: str, work: Path = WORK) -> Path:
    """The driver linked against the lunar build's objects: as built, or with the edge held."""
    folder = work / variant
    folder.mkdir(parents=True, exist_ok=True)
    exe = folder / 'lw_columns'
    objects = [p for p in sorted(BUILD.glob('*.o')) if p.name != 'cm1.o']
    include = ['-I', str(BUILD)]
    if variant == 'edge_held':
        text = (BUILD / MODULE).read_text()
        if text.count(EDGE_EDIT[0]) != 1:
            raise RuntimeError(f'{MODULE} is not the version the edge edit expects')
        (folder / MODULE).write_text(text.replace(EDGE_EDIT[0], EDGE_EDIT[1]))
        subprocess.run(['gfortran', *OPTS, '-J', str(folder), '-I', str(BUILD), '-c', MODULE], cwd=folder, check=True)
        objects = [folder / 'module_ra_rrtmg_lw.o' if p.name == 'module_ra_rrtmg_lw.o' else p for p in objects]
        include = ['-I', str(folder)] + include
    subprocess.run(['gfortran', '-ffree-form', '-ffree-line-length-none', '-O2', '-fopenmp', *include, '-J',
                    str(folder), '-c', str(DRIVER), '-o', str(folder / 'lw_columns.o')], check=True)
    subprocess.run(['gfortran', '-fopenmp', '-o', str(exe), str(folder / 'lw_columns.o'), *map(str, objects)],
                   check=True, capture_output=True)
    return exe


def rrtmg(exe: Path, cols, t_ground, cp, gravity, work: Path = WORK):
    """Clear-sky longwave heating (K/day) of each column's layers from the driver."""
    nlay = cols[0]['layer_p_pa'].size
    fmt = lambda values: ' '.join(f'{float(v):.10e}' for v in values)
    rows = [f'{len(cols)} {nlay} {float(cp):.10e} {float(gravity):.10e}']
    for col, ts in zip(cols, t_ground):
        f = col['dry_fractions']
        xl = col['layer_x_h2o']
        vmr = xl / (1.0 - xl)                                  # RRTMG's mixing ratios are per mole of dry air
        rows += [fmt([ts]), fmt(col['p_pa'] / 100.0), fmt(col['t_k']), fmt(col['layer_p_pa'] / 100.0),
                 fmt(col['layer_t_k']), fmt(vmr), fmt([f['CO2']] * nlay), fmt([f['O2']] * nlay)]
    link = work / 'RRTMG_LW_DATA'                              # rrtmg_lw_ini reads it from the working folder
    if not link.exists():
        link.symlink_to(rrtmg_data())
    out = subprocess.run([str(exe)], input='\n'.join(rows) + '\n', capture_output=True, text=True, check=True,
                         cwd=work).stdout.split('\n')
    return np.array([[float(v) for v in out[3 * i].split()] for i in range(len(cols))])


def rrtmg_data() -> Path:
    """RRTMG's longwave tables as CM1 r22.0 ships them, which every lunar case links into its run folder."""
    return Path('/media/projectspace/terluna-research/cm1/cm1r22.0/run/RRTMG_LW_DATA')


def columns(case_name, from_day, every_snapshot, every_column):
    """The ring's cloud-free columns, sampled as ring_radiation.py samples them, as correlated-k columns with the air
    at the ground at the ground's temperature, as CM1 hands them to RRTMG; with their ground temperatures and CM1's own
    longwave heating at its layer centres, as RRTMG's rate within CM1 and as the lwten CM1 wrote."""
    from climate.crm import ring_analysis as ra
    from climate.crm import ring_radiation as rr
    from atmosphere.radiative_convective import thermodynamics as td
    case = ra.RUNS / case_name
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    tap = geo['record']['configuration']['output_s']
    equation_set = co.eqtset(case)
    air = td.earthlike_air(121590.0, 400.0)
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    span = [n for n in outputs if (n - 1) * tap >= from_day * DAY][::every_snapshot]
    cols, ground, cm1, written, levels = [], [], [], [], []
    for n in span:
        d = ra.read_snapshot(case, n)
        temp = co.temperature(d)
        rate = co.rrtmg_rate(d['lwten'], d, equation_set) * DAY
        low = zh < rr.CLEAR_BELOW_M
        for i in range(0, temp.shape[1], every_column):
            if d['cldfra'][low, i].max() > 0.01 or (d['qc'][low, i] + d['qi'][low, i]).max() > 1e-7:
                continue
            p, t, q = rr.column_levels(zw, zh, d['prs'][:, i], temp[:, i], d['qv'][:, i], d['psfc'][i],
                                       t_ground=float(d['tsk'][i]))
            x = rr.mole_fraction(q, air.molar_mass)
            cols.append(td.column_from_levels(td.MOON, air, p, t, x, constant_gravity=True))
            levels.append((p, t, x))
            ground.append(float(d['tsk'][i]))
            cm1.append(rate[:, i])
            written.append(d['lwten'][:, i] * DAY)
    return dict(cols=cols, levels=levels, ground=ground, cm1=np.array(cm1), cm1_written=np.array(written), zh=zh,
                air=air, snapshots=len(span), eqtset=equation_set)


def twins(found):
    """Each column with its pressures scaled so its ground stands at TWIN_GROUND_HPA, all else kept."""
    from atmosphere.radiative_convective import thermodynamics as td
    out = []
    for p, t, x in found['levels']:
        out.append(td.column_from_levels(td.MOON, found['air'], p * (TWIN_GROUND_HPA * 100.0 / p[0]), t, x,
                                         constant_gravity=True))
    return out


def layer_means(rates, mass, zh, extrapolated):
    """Mass-weighted means over the GCM's lowest four layers and over the layers below the table's edge."""
    rate, w = np.mean(rates, axis=0)[:zh.size], np.mean(mass, axis=0)[:zh.size]
    rows = {}
    for z0, z1 in LAYERS_KM:
        m = (zh >= z0 * 1000.0) & (zh < z1 * 1000.0)
        rows[f'{z0:g}-{z1:g} km'] = float((rate[m] * w[m]).sum() / w[m].sum())
    m = extrapolated[:zh.size]
    rows['below the edge'] = float((rate[m] * w[m]).sum() / w[m].sum())
    rows['above the edge to 13.74 km'] = float((rate[~m & (zh < 13740.0)] * w[~m & (zh < 13740.0)]).sum()
                                               / w[~m & (zh < 13740.0)].sum())
    return rows


def compare(case_name, from_day, every_snapshot=8, every_column=40):
    from atmosphere.radiative_convective import ck
    from atmosphere.radiative_convective import thermodynamics as td
    from climate.crm import ring_radiation as rr
    found = columns(case_name, from_day, every_snapshot, every_column)
    cols, zh = found['cols'], found['zh']
    twin = twins(found)
    g = td.MOON.surface_gravity
    model = ck.CKLongwave()
    ck_cols = np.array([model.heating_rate(c, co.CP, t_surface=ts)[0] for c, ts in zip(cols, found['ground'])])
    ck_twin = np.array([model.heating_rate(c, co.CP, t_surface=ts)[0] for c, ts in zip(twin, found['ground'])])
    built, held = build('as_built'), build('edge_held')
    rr_cols = rrtmg(built, cols, found['ground'], co.CP, g)
    rr_twin = rrtmg(built, twin, found['ground'], co.CP, g)
    rr_held = rrtmg(held, cols, found['ground'], co.CP, g)
    mass = np.array([c['layer_mass_kg_m2'] for c in cols])
    mass_twin = np.array([c['layer_mass_kg_m2'] for c in twin])
    beyond = np.mean([c['layer_p_pa'] / 100.0 > TABLE_EDGE_HPA for c in cols], axis=0) > 0.5
    beyond_twin = np.zeros_like(beyond)
    sets = dict(
        columns=dict(rrtmg=layer_means(rr_cols, mass, zh, beyond), correlated_k=layer_means(ck_cols, mass, zh, beyond),
                     rrtmg_edge_held=layer_means(rr_held, mass, zh, beyond),
                     cm1_as_run=layer_means(found['cm1'], mass, zh, beyond),
                     cm1_as_written=layer_means(found['cm1_written'], mass, zh, beyond)),
        twins=dict(rrtmg=layer_means(rr_twin, mass_twin, zh, beyond), correlated_k=layer_means(ck_twin, mass_twin, zh,
                                                                                                  beyond)))
    for key in sets:
        s = sets[key]
        s['rrtmg_minus_correlated_k'] = {k: s['rrtmg'][k] - s['correlated_k'][k] for k in s['rrtmg']}
    s = sets['columns']
    s['edge_held_minus_as_built'] = {k: s['rrtmg_edge_held'][k] - s['rrtmg'][k] for k in s['rrtmg']}
    s['cm1_minus_rrtmg'] = {k: s['cm1_as_run'][k] - s['rrtmg'][k] for k in s['rrtmg']}
    files = {str(p.relative_to(ROOT)): sha16(p) for p in (Path(__file__), DRIVER, HERE / 'cm1_output.py')}
    return dict(
        schema=SCHEMA, case=case_name, from_day=from_day, eqtset=found['eqtset'], columns=len(cols),
        snapshots=found['snapshots'],
        every_column=every_column, table_edge_hpa=TABLE_EDGE_HPA, twin_ground_hpa=TWIN_GROUND_HPA,
        ground_hpa_mean=float(np.mean([c['p_pa'][0] for c in cols]) / 100.0),
        edge_height_km=float(np.interp(TABLE_EDGE_HPA, np.mean([c['layer_p_pa'] for c in cols], axis=0)[::-1] / 100.0,
                                       np.mean([c['layer_z_m'] for c in cols], axis=0)[::-1]) / 1000.0),
        sets=sets, z_m=zh.tolist(),
        profiles=dict(rrtmg=np.mean(rr_cols, axis=0)[:zh.size].tolist(),
                      correlated_k=np.mean(ck_cols, axis=0)[:zh.size].tolist(),
                      rrtmg_edge_held=np.mean(rr_held, axis=0)[:zh.size].tolist(),
                      rrtmg_twin=np.mean(rr_twin, axis=0)[:zh.size].tolist(),
                      correlated_k_twin=np.mean(ck_twin, axis=0)[:zh.size].tolist(),
                      cm1_as_run=np.mean(found['cm1'], axis=0).tolist(),
                      cm1_as_written=np.mean(found['cm1_written'], axis=0).tolist()),
        producer=dict(files=files, build=str(BUILD), module_sha256=sha16(BUILD / MODULE),
                      rrtmg_lw_data_sha256=sha16(rrtmg_data()), ck_table_key=ck.table_key()),
        evidence=('CM1\'s RRTMG longwave, compiled from the lunar build\'s objects into an offline driver, and the '
                  'atmosphere domain\'s correlated-k model on the same clear-sky layers of the ring\'s cloud-free '
                  'columns, as CM1 hands them to RRTMG, and of twins scaled inside RRTMG\'s tables, with CM1\'s stored '
                  'lwten turned back into RRTMG\'s rate (cm1_output.py); means over the sampled columns. A diagnosis '
                  'of one ring\'s columns over its second lunar day.'),
        reading_rule=('K/day of temperature at constant pressure, negative cools, mass-weighted over the GCM\'s lowest '
                      'four layers, over the layers below the table\'s edge and over those above it to 13.74 km. '
                      'rrtmg_minus_correlated_k on the columns against the twins gives the extrapolation\'s share of '
                      'RRTMG\'s error; edge_held_minus_as_built how far the extrapolated layers depend on the '
                      'extrapolation; cm1_minus_rrtmg CM1\'s call against the same RRTMG offline. cm1_as_written is '
                      'the tendency CM1 stores, at constant volume, cp/cvm times RRTMG\'s rate under eqtset 2.'))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('case')
    ap.add_argument('--from-day', type=float, default=29.5)
    ap.add_argument('--every-snapshot', type=int, default=8)
    ap.add_argument('--every-column', type=int, default=40)
    args = ap.parse_args(argv)
    product = compare(args.case, args.from_day, args.every_snapshot, args.every_column)
    RESULTS.mkdir(parents=True, exist_ok=True)
    out = RESULTS / f'rrtmg_edge_{args.case}.json'
    out.write_text(json.dumps(product, indent=1) + '\n')
    for key, s in product['sets'].items():
        print(key, json.dumps({k: round(v, 3) for k, v in s['rrtmg_minus_correlated_k'].items()}))
    print('edge held - as built', json.dumps({k: round(v, 3) for k, v in
                                              product['sets']['columns']['edge_held_minus_as_built'].items()}))
    print('cm1 - rrtmg', json.dumps({k: round(v, 4) for k, v in product['sets']['columns']['cm1_minus_rrtmg'].items()}))
    print(f'wrote {out.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
