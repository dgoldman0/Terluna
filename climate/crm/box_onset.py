"""Why a forced box's convection starts or stays away, by quarter and by twelfth of the lunar day: its domain-mean air's
buoyancy for parcels from the lowest level and from mixed layers, undiluted and mixing in surrounding air, their
inhibition, condensation and free-convection levels against the boundary layer's depth, the rate at which the
prescribed forcing changes the CAPE (dynamic CAPE, Xie and Zhang 2000), and in the 3-D output the strongest updrafts,
the highest cloud tops, the shallow cloud and CM1's own CAPE and inhibition column by column.

    climate/gcm/.venv/bin/python -m climate.crm.box_onset box_near_side_sea_day3 59 88.5
    climate/gcm/.venv/bin/python -m climate.crm.box_onset box_0e_design 29.5 59
                                                                  # each writes ../results/crm/box_onset_<case>.json

The parcels rise pseudo-adiabatically on CM1's levels through the domain-mean sounding of each 3-hourly output,
keeping their moist static energy h = cp T + g z + L qv and their water (condensate falls out as it forms), at lunar
gravity, with Bolton's saturation over liquid as CM1's own CAPE routine takes it (getcape.F, its pseudo-adiabat without
ice). A mixed-layer parcel starts at the lowest level with the mean potential temperature and vapour of the lowest 0.5
or 1 km by mass. A mixing parcel takes in surrounding air between levels, a share 1 - (p_above/p_below)^eps of its h and
water, eps per e-folding of pressure as the scheme's eps_bm; the column's lowest scale height, R T/g, is about 48 km, so
eps 1 is 0.021/km. If lunar convection is Earth's stretched by the ratio of the gravities, as the gravity pair tests,
mixing per e-folding of pressure carries over unchanged, and Earth's 1/km (Neale and others 2008 used 1e-3/m) is eps 8
here. Buoyancy is the virtual temperature's excess; CAPE sums it from the free-convection level (the first buoyant
level at or above condensation) to the highest buoyant level, the inhibition the negative part below.

The dynamic CAPE is the CAPE's change over an hour under the prescribed forcing alone: the day-night table at that
time (terluna_lsadv.txt, by local time and height) and the mean vertical wind acting on the output's own profile
(terluna_wls.txt), applied above 1 km so the 0.5 km parcel is unchanged (Xie and Zhang's free-troposphere form), and
through the whole column.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from climate.crm import box_budget as bb
from climate.crm import cm1_output as co
from climate.crm import ring_analysis as ra

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.crm-box-onset/1'
DAY = 86400.0
LV = 2.501e6                       # CM1's latent heat of vaporization (constants.F)
EPSV = co.RD / co.RV
SOURCES = {'lowest level': 0.0, 'mixed 0.5 km': 500.0, 'mixed 1 km': 1000.0}
EPS = (0.0, 1.0, 2.0, 4.0, 8.0)
FREE_FROM_M = 1000.0               # the dynamic CAPE's free troposphere starts here
DCAPE_S = 3600.0                   # the dynamic CAPE's interval, an hour
QUARTERS = (('sunrise-noon', -90.0, 0.0), ('noon-sunset', 0.0, 90.0), ('sunset-midnight', 90.0, 180.0),
            ('midnight-sunrise', -180.0, -90.0))
TWELFTHS = 12
SHALLOW_TOP_M, DEEP_TOP_M = 16000.0, 40000.0


def saturation(t, p):
    """Vapour mixing ratio at saturation over liquid (kg/kg) at temperature t (K) and pressure p (Pa), Bolton's."""
    es = np.minimum(611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65)), 0.5 * p)
    return EPSV * es / (p - es)


def virtual(t, q):
    return t * (1.0 + q / EPSV) / (1.0 + q)


def saturate(h, qt, z, p, g):
    """Temperature, vapour and saturation of a parcel of moist static energy h and water qt at height z and pressure
    p: dry where qt stays below saturation at its dry temperature, else on the saturated state of the same h."""
    t_dry = (h - g * z - LV * qt) / co.CP
    if qt <= saturation(t_dry, p):
        return t_dry, qt, False
    t = t_dry
    for _ in range(20):                                      # Newton on cp T + L qs(T) = h - g z
        qs = saturation(t, p)
        es = min(611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65)), 0.5 * p)
        dqs = qs * p / (p - es) * 17.67 * 243.5 / (t - 29.65) ** 2
        step = (co.CP * t + LV * qs - (h - g * z)) / (co.CP + LV * dqs)
        t -= step
        if abs(step) < 1e-6:
            break
    return t, saturation(t, p), True


def source(p, t, q, mass, depth):
    """A parcel's starting temperature and vapour at the lowest level: the lowest level's own for depth 0, or the
    mean potential temperature and vapour by mass of the levels in the lowest `depth` metres (each level's mass the
    share of its layer below depth), carried to the lowest level's pressure."""
    if depth <= 0.0:
        return float(t[0]), float(q[0])
    exner = (p / co.P00) ** (co.RD / co.CP)
    theta = np.sum(t / exner * mass) / np.sum(mass)
    return float(theta * exner[0]), float(np.sum(q * mass) / np.sum(mass))


def layer_share(zw, depth):
    """Each level's share of its layer lying within the lowest `depth` metres."""
    zw = np.asarray(zw, float)
    return np.clip((np.minimum(zw[1:], depth) - zw[:-1]) / np.diff(zw), 0.0, 1.0)


def ascent(z, p, t, q, t0, q0, eps, g):
    """A parcel from the lowest level up through a sounding (levels the ground up), mixing a share
    1 - (p_k/p_(k-1))^eps of its h and water with the surrounding air between levels: its buoyancy (m/s2) and
    whether it is saturated, level by level."""
    h, qt = co.CP * t0 + g * z[0] + LV * q0, q0
    he = co.CP * t + g * z + LV * q
    b, wet = np.zeros(z.size), np.zeros(z.size, bool)
    for k in range(z.size):
        if k and eps > 0.0:
            mix = 1.0 - (p[k] / p[k - 1]) ** eps
            h, qt = (1.0 - mix) * h + mix * he[k], (1.0 - mix) * qt + mix * q[k]
        tk, qk, wet[k] = saturate(h, qt, z[k], p[k], g)
        qt = min(qt, qk)                                     # the condensate falls out
        b[k] = g * (virtual(tk, qk) - virtual(t[k], q[k])) / virtual(t[k], q[k])
    return b, wet


def measures(b, wet, z, dz) -> dict:
    """CAPE and inhibition (J/kg) and the condensation, free-convection and highest buoyant levels (m) from a parcel's
    buoyancy and saturation on levels z of depth dz; NaN for levels the parcel never reaches, and nil CAPE for a parcel
    never buoyant above its condensation."""
    k = np.arange(z.size)
    lcl = int(np.argmax(wet)) if wet.any() else None
    up = (k >= (lcl if lcl is not None else z.size)) & (b > 0.0)
    if not up.any():
        return dict(cape=0.0, cin=None, lcl_m=float(z[lcl]) if lcl is not None else None, lfc_m=None, el_m=None)
    lfc = int(np.argmax(up))
    el = int(np.flatnonzero(b > 0.0).max())
    cape = float(np.sum(np.where((k >= lfc) & (k <= el) & (b > 0.0), b * dz, 0.0)))
    cin = float(np.sum(np.where((k < lfc) & (b < 0.0), -b * dz, 0.0)))
    return dict(cape=cape, cin=cin, lcl_m=float(z[lcl]), lfc_m=float(z[lfc]), el_m=float(z[el]))


def parcel(z, p, t, q, mass, dz, depth, eps, g) -> dict:
    """One parcel's measures: from the lowest level, or from the lowest `depth` metres mixed, at mixing eps."""
    share = layer_share(np.concatenate([[0.0], np.cumsum(dz)]), depth)
    t0, q0 = source(p, t, q, mass * share, depth)
    return measures(*ascent(z, p, t, q, t0, q0, eps, g), z, dz)


def forcing(case: Path, d: dict, zh, zw):
    """The prescribed forcing's temperature (K/s) and vapour (kg/kg/s) tendencies at an output: the day-night table at
    its time and the mean vertical wind acting on its own domain-mean profile."""
    text = (case / 'namelist.input').read_text()
    var18, var19 = bb.namelist_value(text, 'var18'), bb.namelist_value(text, 'var19')
    exner = (d['prs'] / co.P00) ** (co.RD / co.CP)
    th = bb.day_night_heating(*bb.read_day_night(case / 'terluna_lsadv.txt'), zh, d['mtime'], var18, var19)
    qv = bb.day_night_heating(*bb.read_day_night(case / 'terluna_lsadv.txt', 'qv'), zh, d['mtime'], var18, var19)
    wls = np.loadtxt(case / 'terluna_wls.txt', skiprows=1)
    w = np.interp(zh, wls[:, 0], wls[:, 1])
    th = th - w * np.gradient(d['th'], zh)
    qv = qv - w * np.gradient(d['qv'], zh)
    return th * exner, qv


def snapshot(case: Path, n: int, zh) -> dict:
    """From one 3-D output: the strongest updraft and its height, the 99.9th percentile of the vertical wind below and
    above 16 km, the highest cloud top, the shares of columns whose cloud tops lie below 16 km, at 16-40 km and above,
    and CM1's own CAPE, inhibition, condensation and free-convection levels and boundary-layer depth column by
    column."""
    d = ra.read_snapshot(case, n)
    w = d['winterp']
    k = np.unravel_index(int(np.argmax(w)), w.shape)
    cloudy = (d['qc'] + d['qi']) >= ra.CLOUD_KG_KG
    top = np.where(cloudy.any(axis=0), zh[np.where(cloudy, np.arange(zh.size)[:, None], 0).max(axis=0)], np.nan)
    has = np.isfinite(top)
    low, high = zh < SHALLOW_TOP_M, (zh >= SHALLOW_TOP_M) & (zh < 60000.0)
    cape, cin, lfc, hpbl = (np.asarray(d[x], float) for x in ('cape', 'cin', 'lfc', 'hpbl'))
    return dict(w_max_m_s=float(w.max()), w_max_height_m=float(zh[k[0]]),
                w_p999_below_16km_m_s=float(np.percentile(w[low], 99.9)),
                w_p999_16_60km_m_s=float(np.percentile(w[high], 99.9)),
                cloud_top_max_m=float(np.nanmax(top)) if has.any() else None,
                cloud_tops_share=dict(below_16km=float((has & (top < SHALLOW_TOP_M)).mean()),
                                      from_16_to_40km=float((has & (top >= SHALLOW_TOP_M) & (top < DEEP_TOP_M)).mean()),
                                      above_40km=float((has & (top >= DEEP_TOP_M)).mean())),
                cloud_below_16km_share=float(cloudy[low].any(axis=0).mean()),
                columns=dict(cape_mean_j_kg=float(cape.mean()), cape_p90_j_kg=float(np.percentile(cape, 90)),
                             cin_mean_j_kg=float(cin.mean()), cin_below_1_share=float((cin < 1.0).mean()),
                             lcl_mean_m=float(np.mean(d['lcl'])), lfc_mean_m=float(lfc.mean()),
                             lfc_min_m=float(lfc.min()), hpbl_mean_m=float(hpbl.mean()),
                             lfc_within_boundary_layer_share=float((lfc <= hpbl).mean())))


def analyse(case_name: str, span) -> dict:
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    case = ra.RUNS / case_name
    record = json.loads((case / 'case.json').read_text())
    g = record['build']['gravity_m_s2']
    start, lon = record['configuration']['start_hour_angle_deg'], record['configuration']['site']['lon_deg']
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    dz = np.diff(zw)
    tap = record['configuration']['output_s']
    layout = ra.grads_layout(case / 'cm1out_diag_s.ctl')[1]
    rows = []
    for path in sorted(case.glob('cm1out_diag_*_s.dat')):
        d = bb.read_diag(case, path, layout)
        if not span[0] * DAY - 1.0 <= d['mtime'] <= span[1] * DAY + 1.0 or d['mtime'] <= 0.0:
            continue
        t = d['th'] * (d['prs'] / co.P00) ** (co.RD / co.CP)
        p, q, mass = d['prs'], d['qv'], d['rho'] * dz
        row = dict(day=d['mtime'] / DAY, hour=local_hour_angle(start, lon, d['mtime'], solar_day_s()),
                   rain_mm_day=d['prate'] * DAY, hpbl_m=d['hpbl'], cm1_cape_j_kg=d['cape'], cm1_cin_j_kg=d['cin'],
                   parcels={})
        for label, depth in SOURCES.items():
            for eps in EPS:
                row['parcels'][f'{label}, eps {eps:g}'] = parcel(zh, p, t, q, mass, dz, depth, eps, g)
        tt, qq = forcing(case, d, zh, zw)
        free = zh >= FREE_FROM_M
        row['dcape_j_kg_h'] = {}
        for eps in (0.0, 4.0):
            base = row['parcels'][f'mixed 0.5 km, eps {eps:g}']['cape']
            for label, mask in (('free troposphere', free), ('whole column', np.ones_like(free))):
                t2, q2 = t + np.where(mask, tt, 0.0) * DCAPE_S, q + np.where(mask, qq, 0.0) * DCAPE_S
                after = parcel(zh, p, t2, q2, mass, dz, 500.0, eps, g)['cape']
                row['dcape_j_kg_h'][f'{label}, eps {eps:g}'] = after - base
        row['snapshot'] = snapshot(case, int(round(d['mtime'] / tap)) + 1, zh)
        rows.append(row)
    return dict(schema=SCHEMA, case=case_name, span_days=[rows[0]['day'], rows[-1]['day']], outputs=len(rows),
                gravity_m_s2=g, sources_m=SOURCES, eps=list(EPS),
                eps_per_km=[e / (co.RD * 270.0 / g / 1000.0) for e in EPS],
                quarters=summary(rows, [(q, a, b) for q, a, b in QUARTERS]),
                twelfths=summary(rows, [(f'{a:g} to {a + 30.0:g}', a, a + 30.0)
                                        for a in np.linspace(-180.0, 150.0, TWELFTHS)]),
                evidence=('Domain-mean soundings and 3-D output of one CM1 box every 3 hours over one lunar day, one '
                          'realization; parcels computed here on CM1\'s levels, CM1\'s own column CAPE beside them.'),
                reading_rule=('By quarter and by twelfth of the lunar day (hour angle from noon, degrees): the means '
                              'over the outputs in each. parcels: CAPE and inhibition in J/kg and the condensation '
                              '(lcl), free-convection (lfc) and highest buoyant (el) levels in metres for each source '
                              '(the lowest level, or the lowest 0.5 or 1 km mixed by mass) and mixing eps per '
                              'e-folding of pressure (eps_per_km at 270 K); inhibition and levels are means over the '
                              'outputs where the parcel reaches them, with the share that do. dcape: the CAPE\'s '
                              'change over an hour (J/kg/h) under the prescribed forcing, above 1 km or through the '
                              'column, for the 0.5 km parcel undiluted and at eps 4. snapshot: the strongest '
                              'updraft, the 99.9th percentile of the vertical wind, the highest cloud top (cloud: '
                              'liquid and ice at least 0.01 g/kg), the shares of columns by cloud top, and CM1\'s own '
                              'column values '
                              '(its 0.5 km mixed-layer parcel; hpbl from its bulk Richardson number).'))


def summary(rows, bins) -> dict:
    out = {}
    for label, a, b in bins:
        inside = [r for r in rows if a <= r['hour'] < b]
        if not inside:
            continue
        mean = lambda xs: float(np.mean(xs)) if xs else None
        entry = dict(outputs=len(inside), rain_mm_day=mean([r['rain_mm_day'] for r in inside]),
                     hpbl_m=mean([r['hpbl_m'] for r in inside]),
                     cm1_cape_j_kg=mean([r['cm1_cape_j_kg'] for r in inside]),
                     cm1_cin_j_kg=mean([r['cm1_cin_j_kg'] for r in inside]), parcels={}, dcape_j_kg_h={}, snapshot={})
        for key in inside[0]['parcels']:
            vals = [r['parcels'][key] for r in inside]
            entry['parcels'][key] = dict(cape=mean([v['cape'] for v in vals]),
                                         **{m: mean([v[m] for v in vals if v[m] is not None])
                                            for m in ('cin', 'lcl_m', 'lfc_m', 'el_m')},
                                         reaching_lfc_share=float(np.mean([v['lfc_m'] is not None for v in vals])))
        for key in inside[0]['dcape_j_kg_h']:
            entry['dcape_j_kg_h'][key] = mean([r['dcape_j_kg_h'][key] for r in inside])
        snaps = [r['snapshot'] for r in inside]
        for key in ('w_max_m_s', 'w_max_height_m', 'w_p999_below_16km_m_s', 'w_p999_16_60km_m_s',
                    'cloud_below_16km_share'):
            entry['snapshot'][key] = mean([s[key] for s in snaps])
        tops = [s['cloud_top_max_m'] for s in snaps if s['cloud_top_max_m'] is not None]
        entry['snapshot']['cloud_top_max_m'] = dict(mean=mean(tops), lowest=min(tops) if tops else None,
                                                    highest=max(tops) if tops else None)
        entry['snapshot']['cloud_tops_share'] = {k: mean([s['cloud_tops_share'][k] for s in snaps])
                                                 for k in snaps[0]['cloud_tops_share']}
        entry['snapshot']['columns'] = {k: mean([s['columns'][k] for s in snaps]) for k in snaps[0]['columns']}
        out[label] = entry
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('first_day', type=float)
    parser.add_argument('last_day', type=float)
    args = parser.parse_args(argv)
    out = analyse(args.case, (args.first_day, args.last_day))
    out['producer'] = dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                    for p in (Path(__file__), HERE / 'box_budget.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f'box_onset_{args.case}.json'
    path.write_text(json.dumps(out, indent=1) + '\n')
    for label, e in out['quarters'].items():
        s, c = e['snapshot'], e['snapshot']['columns']
        print(f"{label}: rain {e['rain_mm_day']:.2f}, hpbl {e['hpbl_m']:.0f} m, CM1 CAPE {e['cm1_cape_j_kg']:.0f} "
              f"CIN {e['cm1_cin_j_kg']:.1f}; w max {s['w_max_m_s']:.2f} m/s at {s['w_max_height_m'] / 1000:.1f} km; "
              f"cloud top {(s['cloud_top_max_m']['mean'] or float('nan')) / 1000:.1f} km; "
              f"shallow {s['cloud_tops_share']['below_16km']:.2f} deep {s['cloud_tops_share']['above_40km']:.2f}; "
              f"columns lcl {c['lcl_mean_m']:.0f} lfc {c['lfc_mean_m']:.0f} (min {c['lfc_min_m']:.0f}) hpbl "
              f"{c['hpbl_mean_m']:.0f}, lfc within the boundary layer {c['lfc_within_boundary_layer_share']:.2f}")
        for key, v in e['parcels'].items():
            fmt = lambda x: f'{x:7.0f}' if x is not None else '      -'
            print(f"   {key:22s} CAPE {v['cape']:7.0f} CIN {v['cin'] if v['cin'] is not None else float('nan'):6.1f} "
                  f"lcl {fmt(v['lcl_m'])} lfc {fmt(v['lfc_m'])} el {fmt(v['el_m'])} "
                  f"reach {v['reaching_lfc_share']:.2f}")
        print('   dCAPE J/kg/h', {k: round(v, 1) for k, v in e['dcape_j_kg_h'].items()})
    print(f'wrote {path.relative_to(HERE.parents[1])}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
