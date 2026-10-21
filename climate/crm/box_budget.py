"""The heat of a box's air by process, from CM1's domain diagnostics, against a reference box at the same site.

    climate/gcm/.venv/bin/python -m climate.crm.box_budget box_0e_small_lunar_rain_day29 \\
        --reference box_0e_small_rain_control_day29 --from-day 29.5

CM1 writes the domain means of its potential-temperature budget by process at every diagnostic output (cm1out_diag_*,
with dodomaindiag), the microphysics' share split into condensation, evaporation from cloud and from rain,
deposition, sublimation, melting and freezing. Each term is taken as the temperature tendency at constant pressure,
pi times its potential-temperature tendency, and over the span as its mean by height and its sum over the GCM's lowest
four layers in W/m2 (rho cp times it, level by level). The box's prescribed day-night heating (terluna_lsadv.txt, var11)
is added to the potential temperature outside CM1's budget, so it is reckoned from its table by local time as the
patch applies it; the storage is the trend of the domain-mean temperature, and what is left closes the budget. At the
ground: the heat carried into the air, the evaporation, the radiation, the air at 2 m, the rain and the cloud fraction
of the lowest level, about the lowest 100 m. Writes results/crm/box_budget_<case>.json.
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
SCHEMA = 'terluna.climate.crm-box-budget/1'
DAY = 86400.0
XLV = 2501000.0                   # J/kg, CM1's latent heat of vaporization (constants.F)
LAYERS_KM = ((0.0, 1.77), (1.77, 6.53), (6.53, 13.74), (13.74, 23.06))   # the GCM's lowest four, as elsewhere
GROUPS = {                        # CM1's potential-temperature budget terms (domain diagnostics) by process
    'radiation': ('ptb_rad',),
    'microphysics': ('ptb_mp',),
    'turbulence': ('ptb_pbl', 'ptb_hidiff', 'ptb_vidiff', 'ptb_diss'),
    'resolved_transport': ('ptb_hadv', 'ptb_vadv', 'ptb_div'),
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


def read_day_night(path: Path):
    """terluna_lsadv.txt: the level heights (m) and, by local-time bin, the potential-temperature tendency (K/s)."""
    rows = path.read_text().split('\n')
    nb, nl = (int(v) for v in rows[0].split())
    z = np.array([float(v) for v in rows[1].split()])
    th = np.array([[float(v) for v in rows[2 + 2 * m].split()] for m in range(nb)])
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
    times, t_mean, rho, day_night = [], [], [], []
    terms = {k: [] for k in GROUPS}
    parts = {k: [] for k in MICROPHYSICS}
    surface = {k: [] for k in ('t2_k', 'q2_g_kg', 'sensible_w_m2', 'latent_w_m2', 'net_sunlight_w_m2',
                               'net_infrared_w_m2', 'rain_mm_day', 'cloud_fraction_lowest_level')}
    for path in sorted(case.glob('cm1out_diag_*_s.dat')):
        d = read_diag(case, path, layout)
        if d['mtime'] < from_day * DAY:
            continue
        exner = (d['prs'] / co.P00) ** (co.RD / co.CP)
        times.append(d['mtime'] / DAY)
        t_mean.append(d['th'] * exner)
        rho.append(d['rho'])
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
    return dict(case=name, span_days=[float(times[0]), float(times[-1])], outputs=int(times.size),
                z_m=zh[top].tolist(),
                profiles_k_day={k: v[top].tolist() for k, v in prof.items()},
                microphysics_k_day={k: v[top].tolist() for k, v in micro.items()},
                layers_w_m2={k: layer_sums(v, rho, zw) for k, v in prof.items()},
                microphysics_layers_w_m2={k: layer_sums(v, rho, zw) for k, v in micro.items()},
                surface={k: float(np.mean(v)) for k, v in surface.items()})


def difference(a: dict, b: dict) -> dict:
    """a minus b, layer by layer and at the ground."""
    sub = lambda x, y: {k: (sub(x[k], y[k]) if isinstance(x[k], dict) else x[k] - y[k]) for k in x}
    return dict(layers_w_m2=sub(a['layers_w_m2'], b['layers_w_m2']),
                microphysics_layers_w_m2=sub(a['microphysics_layers_w_m2'], b['microphysics_layers_w_m2']),
                surface=sub(a['surface'], b['surface']))


def analyse(name: str, reference: str | None, from_day: float) -> dict:
    out = dict(schema=SCHEMA, case=name, from_day=from_day, budget=budget(name, from_day))
    if reference:
        out['reference'] = reference
        out['reference_budget'] = budget(reference, from_day)
        out['difference'] = difference(out['budget'], out['reference_budget'])
    out['evidence'] = ('Domain means of CM1\'s own budget terms, written every diagnostic output, and the patch\'s '
                       'prescribed day-night heating reckoned from its table; means over the span. One realization '
                       'of each run.')
    out['reading_rule'] = ('K/day of temperature at constant pressure (pi times the potential-temperature tendencies) by '
                           'height, and W/m2 summed over the GCM\'s lowest four layers; positive warms. Turbulence '
                           'holds the boundary layer scheme and the advection\'s implicit diffusion; resolved transport '
                           'the box\'s own circulation; rest closes the budget against the storage and holds what CM1 '
                           'does not budget. At the ground the heat into the air and the evaporation are positive '
                           'upward, the radiation positive downward; the difference is the case minus the reference.')
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
    for label, b in (('case', out['budget']), ('reference', out.get('reference_budget'))):
        if b is None:
            continue
        print(label, b['case'], 'surface', {k: round(v, 2) for k, v in b['surface'].items()})
        for layer in b['layers_w_m2']['radiation']:
            print(f'  {layer}', {k: round(v[layer], 1) for k, v in b['layers_w_m2'].items()})
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
