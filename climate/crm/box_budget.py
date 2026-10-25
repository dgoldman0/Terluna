"""The heat and water of a box's air by process, from CM1's domain diagnostics, against a reference box at the same
site.

    climate/gcm/.venv/bin/python -m climate.crm.box_budget box_0e_small_lunar_rain_day29 \\
        --reference box_0e_small_rain_control_day29 --from-day 29.5

CM1 writes the domain means of its potential-temperature budget by process at every diagnostic output (cm1out_diag_*,
with dodomaindiag), the microphysics' share split into condensation, evaporation from cloud and from rain,
deposition, sublimation, melting and freezing. Each term is taken as the temperature tendency at constant pressure,
pi times its potential-temperature tendency, and over the span as its mean by height and its sum over the GCM's lowest
four layers and over three bands of the column (0-16, 16-60 and 60 km to the top) in W/m2 (rho cp times it, level by
level). The box's prescribed day-night heating (terluna_lsadv.txt, var11) is added to the potential temperature
outside CM1's budget, so it is reckoned from its table by local time as the patch applies it; the storage is the
trend of the domain-mean temperature, and what is left closes the budget. At the ground: the heat carried into the
air, the evaporation, the radiation, the air at 2 m, the rain and the cloud fraction of the lowest level, about the
lowest 100 m.

The water: the column's vapour by CM1's own budget terms (rho times each, level by level, in mm/day), the table's
prescribed moistening reckoned as for the heat, against the evaporation, the rain and the change of the column's
vapour and condensate from the span's first output to its last; what is left closes it. CM1 splits its vertical
advection from its implicit diffusion with the imposed sinking taken out twice (adv.F: the slot for the diffusion
holds the sinking too), so only their sum with the horizontal advection, the box's own transport, stands; over a
doubly periodic box it is near nil. The boundary layer's column over the evaporation, and its heating over the
sensible heat, should be one: they stood at 1.32-1.51 in CM1's simple boundary layer before its density was divided
out (cm1_run.PBL_DENSITY_PATCH). Writes results/crm/box_budget_<case>.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np

from climate.crm import cm1_output as co
from climate.crm import ring_analysis as ra

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-box-budget/2'
DAY = 86400.0
XLV = 2501000.0                   # J/kg, CM1's latent heat of vaporization (constants.F)
LAYERS_KM = ((0.0, 1.77), (1.77, 6.53), (6.53, 13.74), (13.74, 23.06))   # the GCM's lowest four, as elsewhere
BANDS_KM = ((0.0, 16.0), (16.0, 60.0), (60.0, 1000.0))                  # the first boxes' forcing top, and above
WATER = {                         # CM1's vapour budget terms (domain diagnostics) by process
    'boundary_layer': ('qvb_pbl', 'qvb_vturb', 'qvb_hturb'),
    'microphysics': ('qvb_mp',),
    'own_transport': ('qvb_hadv', 'qvb_vadv', 'qvb_hidiff', 'qvb_vidiff', 'qvb_hediff', 'qvb_vediff'),
    'imposed_ascent': ('qvb_lsw',),
    'nudging': ('qvb_nudge',),
}
CONDENSATE = ('qc', 'qr', 'qi', 'qs', 'qg')  # the condensate's mixing ratios, kg/kg by level
GROUPS = {                        # CM1's potential-temperature budget terms (domain diagnostics) by process
    'radiation': ('ptb_rad',),
    'microphysics': ('ptb_mp',),
    'turbulence': ('ptb_pbl', 'ptb_diss'),
    'resolved_transport': ('ptb_hadv', 'ptb_vadv', 'ptb_div', 'ptb_hidiff', 'ptb_vidiff'),
    'imposed_sinking': ('ptb_lsw',),
    'nudging': ('ptb_nudge',),
    'damping': ('ptb_rdamp',),
}
MICROPHYSICS = {'condensation': 'tt_cond', 'evaporation_from_cloud': 'tt_evac', 'evaporation_from_rain': 'tt_evar',
                'deposition': 'tt_dep', 'sublimation': 'tt_subl', 'melting': 'tt_melt', 'freezing': 'tt_frz'}


def read_diag(case: Path, path: Path, layout=None) -> dict:
    """One domain-diagnostic output: scalars as floats, profiles as arrays by level."""
    if layout is None:
        layout = ra.grads_layout(case / 'cm1out_diag_s.ctl')[1]
    raw = np.fromfile(path, dtype='<f4')
    out, offset = {}, 0
    for name, levels in layout:
        out[name] = raw[offset:offset + levels].astype(float) if levels > 1 else float(raw[offset])
        offset += levels
    if offset != raw.size:
        raise ValueError(f'{path.name}: {raw.size} values, layout expects {offset}')
    return out


def namelist_value(text: str, key: str) -> float:
    m = re.search(rf'^\s*{key}\s*=\s*([-+0-9.eEdD]+)', text, re.M)
    if not m:
        raise KeyError(key)
    return float(m.group(1).replace('d', 'e').replace('D', 'e'))


def read_day_night(path: Path, field: str = 'theta'):
    """terluna_lsadv.txt: the level heights (m) and, by local-time bin, the potential-temperature tendency (K/s), or
    with field 'qv' the vapour's (kg/kg/s)."""
    rows = path.read_text().split('\n')
    nb, nl = (int(v) for v in rows[0].split())
    z = np.array([float(v) for v in rows[1].split()])
    offset = {'theta': 2, 'qv': 3}[field]
    th = np.array([[float(v) for v in rows[offset + 2 * m].split()] for m in range(nb)])
    if z.size != nl or th.shape != (nb, nl):
        raise ValueError(f'{path.name}: expected {nb} bins of {nl} levels')
    return z, th


def day_night_heating(z_table, th_table, zh, mtime, var18, var19):
    """The patch's potential-temperature tendency at heights zh and model time mtime (s): linear in local time between
    bin centres, evenly spaced in hour angle from -180 degrees, and linear in height, zero outside the table."""
    nb = th_table.shape[0]
    h = ((var18 + 360.0 * mtime / var19 + 180.0) % 360.0) * nb / 360.0 - 0.5
    m = int(np.floor(h))
    f = h - m
    profile = (1.0 - f) * th_table[m % nb] + f * th_table[(m + 1) % nb]
    inside = (zh >= z_table[0]) & (zh <= z_table[-1])
    return np.where(inside, np.interp(zh, z_table, profile), 0.0)


def layer_sums(rate_k_day, rho, zw, layers=LAYERS_KM):
    """W/m2 in each layer and in the whole column: rho cp times a temperature tendency (K/day) on scalar levels, each
    level by its depth."""
    dz = np.diff(zw)
    zh = 0.5 * (zw[1:] + zw[:-1])
    energy = np.asarray(rate_k_day) * rho * co.CP * dz / DAY
    out = {}
    for z0, z1 in layers:
        m = (zh >= z0 * 1000.0) & (zh < z1 * 1000.0)
        out[f'{z0:g}-{z1:g} km'] = float(energy[m].sum())
    out['column'] = float(energy.sum())
    return out


def layer_means(values, rho, zw, layers=LAYERS_KM):
    """The mass-weighted mean of a profile over each of the GCM's lowest four layers."""
    mass = np.asarray(rho, float) * np.diff(np.asarray(zw, float))
    zh = 0.5 * (np.asarray(zw[1:], float) + np.asarray(zw[:-1], float))
    out = {}
    for z0, z1 in layers:
        m = (zh >= z0 * 1000.0) & (zh < z1 * 1000.0)
        out[f'{z0:g}-{z1:g} km'] = float((np.asarray(values, float)[m] * mass[m]).sum() / mass[m].sum())
    return out


def drift(name: str, from_day: float) -> dict:
    """How the box's air moved from its start: the domain-mean temperature (C) and vapour (g/kg) of each of the GCM's
    lowest four layers at every diagnostic output, with their values at the start, at from_day and at the end, and
    their trend per day over the span after from_day."""
    case = ra.RUNS / name
    geo = ra.case_geometry(case)
    layout = ra.grads_layout(case / 'cm1out_diag_s.ctl')[1]
    times, series = [], {}
    for path in sorted(case.glob('cm1out_diag_*_s.dat')):
        d = read_diag(case, path, layout)
        times.append(d['mtime'] / DAY)
        for quantity, values in (('temperature_c', d['t'] - 273.15), ('vapour_g_kg', d['qv'] * 1000.0)):
            for layer, v in layer_means(values, d['rho'], geo['zw']).items():
                series.setdefault(layer, {}).setdefault(quantity, []).append(v)
    times = np.array(times)
    late = times >= from_day
    out = dict(times_day=[round(float(x), 4) for x in times], layers={})
    for layer, quantities in series.items():
        out['layers'][layer] = {q: dict(series=[round(float(x), 4) for x in v], start=float(v[0]),
                                        at_from_day=float(np.interp(from_day, times, v)), end=float(v[-1]),
                                        trend_per_day=float(np.polyfit(times[late], np.array(v)[late], 1)[0])
                                        if late.sum() > 1 else None)
                                for q, v in quantities.items()}
    return out


def budget(name: str, from_day: float) -> dict:
    """A box's heat budget over the span after from_day, as profiles (K/day of temperature) and layer sums (W/m2)."""
    case = ra.RUNS / name
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    text = (case / 'namelist.input').read_text()
    prescribed = namelist_value(text, 'var11') > 0.5
    if prescribed:
        table = read_day_night(case / 'terluna_lsadv.txt')
        var18, var19 = namelist_value(text, 'var18'), namelist_value(text, 'var19')
    layout = ra.grads_layout(case / 'cm1out_diag_s.ctl')[1]
    times, t_mean, rho, day_night, pbl = [], [], [], [], []
    terms = {k: [] for k in GROUPS}
    parts = {k: [] for k in MICROPHYSICS}
    surface = {k: [] for k in ('t2_k', 'q2_g_kg', 'sensible_w_m2', 'latent_w_m2', 'net_sunlight_w_m2',
                               'net_infrared_w_m2', 'rain_mm_day', 'cloud_fraction_lowest_level')}
    dz = np.diff(zw)
    for path in sorted(case.glob('cm1out_diag_*_s.dat')):
        d = read_diag(case, path, layout)
        if d['mtime'] < from_day * DAY:
            continue
        exner = (d['prs'] / co.P00) ** (co.RD / co.CP)
        times.append(d['mtime'] / DAY)
        t_mean.append(d['th'] * exner)
        rho.append(d['rho'])
        pbl.append(float((d['ptb_pbl'] * exner * d['rho'] * dz).sum() * co.CP))
        for k, names in GROUPS.items():
            terms[k].append(sum(d[n] for n in names) * exner * DAY)
        for k, n in MICROPHYSICS.items():
            parts[k].append(d[n] * exner * DAY)
        if prescribed:
            day_night.append(day_night_heating(*table, zh, d['mtime'], var18, var19) * exner * DAY)
        surface['t2_k'].append(d['t2'])
        surface['q2_g_kg'].append(d['q2'] * 1000.0)
        surface['sensible_w_m2'].append(d['hfx'])
        surface['latent_w_m2'].append(d['qfx'] * XLV)
        surface['net_sunlight_w_m2'].append(d['swdnb'] - d['swupb'])
        surface['net_infrared_w_m2'].append(d['lwdnb'] - d['lwupb'])
        surface['rain_mm_day'].append(d['prate'] * DAY)
        surface['cloud_fraction_lowest_level'].append(d['cldfra'][0])
    if len(times) < 2:
        raise RuntimeError(f'{name}: fewer than two diagnostic outputs after day {from_day}')
    times = np.array(times)
    rho = np.mean(rho, axis=0)
    prof = {k: np.mean(v, axis=0) for k, v in terms.items()}
    if day_night:
        prof['day_night_prescribed'] = np.mean(day_night, axis=0)
    prof['storage'] = np.polyfit(times, np.array(t_mean), 1)[0]
    prof['rest'] = prof['storage'] - sum(v for k, v in prof.items() if k != 'storage')
    micro = {k: np.mean(v, axis=0) for k, v in parts.items()}
    top = zh < 23060.0
    sensible = float(np.mean(surface['sensible_w_m2']))
    return dict(case=name, span_days=[float(times[0]), float(times[-1])], outputs=int(times.size),
                z_m=zh[top].tolist(),
                profiles_k_day={k: v[top].tolist() for k, v in prof.items()},
                microphysics_k_day={k: v[top].tolist() for k, v in micro.items()},
                layers_w_m2={k: layer_sums(v, rho, zw) for k, v in prof.items()},
                bands_w_m2={k: layer_sums(v, rho, zw, BANDS_KM) for k, v in prof.items()},
                microphysics_layers_w_m2={k: layer_sums(v, rho, zw) for k, v in micro.items()},
                boundary_layer_heating_over_sensible=float(np.mean(pbl)) / sensible if abs(sensible) > 1.0 else None,
                surface={k: float(np.mean(v)) for k, v in surface.items()})


def water(name: str, from_day: float, to_day: float | None = None) -> dict:
    """A box's vapour budget over the span after from_day, to the case's end or to to_day: each process's column and
    band sums (mm/day of water), the evaporation, the rain, the storage from the first output to the last and what
    closes it, and the boundary layer's column over the evaporation."""
    case = ra.RUNS / name
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    dz = np.diff(zw)
    text = (case / 'namelist.input').read_text()
    prescribed = namelist_value(text, 'var11') > 0.5
    if prescribed:
        z_table, qv_table = read_day_night(case / 'terluna_lsadv.txt', 'qv')
        var18, var19 = namelist_value(text, 'var18'), namelist_value(text, 'var19')
    layout = ra.grads_layout(case / 'cm1out_diag_s.ctl')[1]
    present = {name for name, _ in layout}
    groups = {k: [n for n in names if n in present] for k, names in WATER.items()}
    times, evaporation, rain, held = [], [], [], []
    columns = {k: [] for k in list(groups) + (['day_night_prescribed'] if prescribed else [])}
    bands = {k: [] for k in columns}
    for path in sorted(case.glob('cm1out_diag_*_s.dat')):
        d = read_diag(case, path, layout)
        if d['mtime'] < from_day * DAY or (to_day is not None and d['mtime'] > to_day * DAY + 1.0):
            continue
        mass = d['rho'] * dz
        times.append(d['mtime'] / DAY)
        evaporation.append(d['qfx'] * DAY)
        rain.append(d['prate'] * DAY)
        held.append(float(((d['qv'] + sum(d[n] for n in CONDENSATE if n in present)) * mass).sum()))
        rates = {k: sum((d[n] for n in names), np.zeros_like(zh)) for k, names in groups.items()}
        if prescribed:
            rates['day_night_prescribed'] = day_night_heating(z_table, qv_table, zh, d['mtime'], var18, var19)
        for k, r in rates.items():
            columns[k].append(float((r * mass).sum()) * DAY)
            bands[k].append([float((r * mass)[(zh >= a * 1000.0) & (zh < b * 1000.0)].sum()) * DAY
                             for a, b in BANDS_KM])
    if len(times) < 2:
        raise RuntimeError(f'{name}: fewer than two diagnostic outputs after day {from_day}')
    column = {k: float(np.mean(v)) for k, v in columns.items()}
    storage = (held[-1] - held[0]) / (times[-1] - times[0])
    e, p = float(np.mean(evaporation)), float(np.mean(rain))
    rest = storage - (sum(v for k, v in column.items() if k != 'microphysics') - p)
    return dict(case=name, span_days=[float(times[0]), float(times[-1])], outputs=len(times),
                column_mm_day=column,
                bands_mm_day={k: dict(zip((f'{a:g}-{b:g} km' for a, b in BANDS_KM), np.mean(v, axis=0).tolist()))
                              for k, v in bands.items()},
                evaporation_mm_day=e, rain_mm_day=p, storage_mm_day=storage,
                column_water_kg_m2=dict(start=held[0], end=held[-1]),
                microphysics_less_rain_mm_day=column['microphysics'] + p, rest_mm_day=rest,
                boundary_layer_over_evaporation=column['boundary_layer'] / e if e > 0.0 else None)


def difference(a: dict, b: dict) -> dict:
    """a minus b, layer by layer and at the ground."""
    sub = lambda x, y: {k: (sub(x[k], y[k]) if isinstance(x[k], dict) else x[k] - y[k]) for k in x}
    return dict(layers_w_m2=sub(a['layers_w_m2'], b['layers_w_m2']),
                microphysics_layers_w_m2=sub(a['microphysics_layers_w_m2'], b['microphysics_layers_w_m2']),
                surface=sub(a['surface'], b['surface']))


def analyse(name: str, reference: str | None, from_day: float) -> dict:
    out = dict(schema=SCHEMA, case=name, from_day=from_day, budget=budget(name, from_day), water=water(name, from_day),
               drift=drift(name, from_day))
    if reference:
        out['reference'] = reference
        out['reference_budget'] = budget(reference, from_day)
        out['reference_water'] = water(reference, from_day)
        out['difference'] = difference(out['budget'], out['reference_budget'])
    out['evidence'] = ('Domain means of CM1\'s own budget terms, written every diagnostic output, and the patch\'s '
                       'prescribed day-night heating and moistening reckoned from its table; means over the span. One '
                       'realization of each run.')
    out['reading_rule'] = ('Heat: K/day of temperature at constant pressure (pi times the potential-temperature '
                           'tendencies) by height, and W/m2 summed over the GCM\'s lowest four layers and over the '
                           'bands 0-16, 16-60 and 60 km to the top; positive warms. Turbulence holds the boundary layer '
                           'scheme and the dissipation; resolved transport the box\'s own circulation with its '
                           'advection\'s implicit diffusion, near nil over the column; rest closes the budget against '
                           'the storage and holds what CM1 does not budget. Water: mm/day of the column\'s vapour by '
                           'process and band, positive moistening; rest closes the vapour and condensate against their '
                           'storage and the rain. At the ground the heat into the air and the evaporation are positive '
                           'upward, the radiation positive downward; the difference is the case minus the reference. '
                           'Products before schema 2 counted the imposed sinking\'s heating in turbulence and took it '
                           'from resolved transport.')
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--reference')
    parser.add_argument('--from-day', type=float, default=29.5)
    args = parser.parse_args(argv)
    out = analyse(args.case, args.reference, args.from_day)
    out['producer'] = dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'cm1_output.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f'box_budget_{args.case}.json'
    path.write_text(json.dumps(out, indent=1) + '\n')
    for label, b, w in (('case', out['budget'], out['water']),
                        ('reference', out.get('reference_budget'), out.get('reference_water'))):
        if b is None:
            continue
        print(label, b['case'], 'surface', {k: round(v, 2) for k, v in b['surface'].items()})
        for layer in b['layers_w_m2']['radiation']:
            print(f'  {layer}', {k: round(v[layer], 1) for k, v in b['layers_w_m2'].items()})
        for band in b['bands_w_m2']['radiation']:
            print(f'  {band}', {k: round(v[band], 1) for k, v in b['bands_w_m2'].items()})
        print(f"  water (mm/day): evaporation {w['evaporation_mm_day']:.3f}, rain {w['rain_mm_day']:.3f}, storage "
              f"{w['storage_mm_day']:.3f}, rest {w['rest_mm_day']:.3f};",
              {k: round(v, 3) for k, v in w['column_mm_day'].items()})
        ratio = lambda x: 'n/a' if x is None else f'{x:.3f}'
        print(f"  boundary layer over the surface: water {ratio(w['boundary_layer_over_evaporation'])}, heat "
              f"{ratio(b['boundary_layer_heating_over_sensible'])}")
    if 'difference' in out:
        print('difference surface', {k: round(v, 2) for k, v in out['difference']['surface'].items()})
        for layer in out['difference']['layers_w_m2']['radiation']:
            diff = out['difference']
            print(f'  {layer}', {k: round(v[layer], 1) for k, v in diff['layers_w_m2'].items()})
            print('     microphysics', {k: round(v[layer], 1) for k, v in diff['microphysics_layers_w_m2'].items()})
    print(f'wrote {path.relative_to(HERE.parent.parent)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
