"""The sea box's third lunar day with Earth's boundary-layer length beside the stretched one (27 October 2026; the record
is the section "The sea box with Earth's boundary-layer length" in README.md beside this file).

    climate/gcm/.venv/bin/python -m climate.crm.box_length     # ../results/crm/box_length_near_side_sea.json, a minute

CM1's simple boundary layer (ipbl 2) mixes with K = l^2 sqrt(max(S^2 - N^2/Pr, small)), its length
l = (l_inf^-2 + (kappa z)^-2)^-1/2, and cm1_run stretches the asymptotic length l_inf by the gravities from Earth's 75 m
to 453 m. Case box_near_side_sea_day3_linf75 runs the sea box's third lunar day again from the same day-59 restart with
l_inf at 75 m; its domain diagnostics come every 3 hours, as box_near_side_sea_day3's do, and its 3-D output every 6.
It writes fewer 3-D fields than its source did, so CM1 puts all its output under a second output configuration, the
files suffixed .02 and numbered from 1, which the readers here follow (box_lid.diag_files, output_times, stats_rows).
By quarter of the lunar day this reads both runs over days 59-88.5:
- from the domain diagnostics: the mixed layer's depth (pbl_mixing's rule: the first height above 0.31 km where the
  mean virtual potential temperature has risen 0.5 K above its value there), the vapour at fixed heights and its drop
  from 0.31 to 6.2 km, the virtual potential temperature's rise from 6.2 to 8.2 km, CM1's own boundary-layer depth
  (hpbl, the mean of its columns' bulk-Richardson heights) and the highest level holding any cloud;
- the rain at the ground, from CM1's running domain total;
- the cloud tops from the 3-D outputs at the times both runs wrote one (every 6 hours): the share of columns holding
  cloud (liquid and ice of at least 0.01 g/kg) below 40 km, the median and 90th percentile of those columns' highest
  cloudy level, and the share whose cloud reaches above 16 km;
- the heat and vapour carried through the lid at 6.45 and 6.99 km by part, resolved, the advection's implicit diffusion
  and the boundary-layer scheme's, with the scheme's diffusivity (box_lid).
The verdict after midnight follows CRITERIA, set before the run. The held night at 6 km with the same length
(box_near_side_sea_night_6km_linf75), when it has run, is read beside box_near_side_sea_night_6km over its second day.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from climate.crm import box_budget as bb
from climate.crm import box_lid as bl
from climate.crm import cm1_output as co
from climate.crm import ring_analysis as ra

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
OUT = RESULTS / 'box_length_near_side_sea.json'
SCHEMA = 'terluna.climate.crm-box-length/1'
DAY, HOUR = 86400.0, 3600.0
RUNS = dict(stretched='box_near_side_sea_day3', earth='box_near_side_sea_day3_linf75')
SPAN_DAYS = (59.0, 88.5)
NIGHTS = dict(stretched='box_near_side_sea_night_6km', earth='box_near_side_sea_night_6km_linf75')
Z_LOW_M, Z_MID_M, Z_LID_M = 314.0, 6201.0, 8238.0   # the GCM's full levels at the sea site (pbl_mixing_l15.json)
PROFILE_M = (314.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6201.0, 7000.0, 8238.0, 10000.0)
DEEP_M, TOP_M = 16000.0, 40000.0
# Set before the run (27 October). After midnight the run with Earth's length is robust if its mixed depth lies within
# 1 km of the stretched run's 6.48 km and its 0.31-6.2 km vapour drop within 1 g/kg of 0.35 g/kg; it depends on the
# length if the depth falls by more than 1.5 km or the drop grows by more than 1 g/kg; otherwise it is ambiguous.
CRITERIA = dict(reference_depth_km=6.48, reference_drop_g_kg=0.35, robust_depth_km=1.0, robust_drop_g_kg=1.0,
                dependent_depth_fall_km=1.5, dependent_drop_rise_g_kg=1.0)


def verdict(depth_km: float, drop_g_kg: float, c=CRITERIA) -> str:
    """The criteria's reading of the run with Earth's length after midnight."""
    if (c['reference_depth_km'] - depth_km > c['dependent_depth_fall_km']
            or drop_g_kg - c['reference_drop_g_kg'] > c['dependent_drop_rise_g_kg']):
        return 'dependent on the length'
    if (abs(depth_km - c['reference_depth_km']) <= c['robust_depth_km']
            and abs(drop_g_kg - c['reference_drop_g_kg']) <= c['robust_drop_g_kg']):
        return 'robust'
    return 'ambiguous'


def configurations(case: Path, stem: str) -> list:
    """The suffixes of a case's output configurations that hold a data file of one kind (stem 'cm1out_stats' or
    'cm1out_metadata'): '' for CM1's first configuration, '.02' and on for later ones (box_lid.diag_files)."""
    import re
    return [m.group(1) or '' for p in sorted(case.glob(f'{stem}*.dat'))
            if (m := re.match(rf'{re.escape(stem)}(\.\d+)?\.dat$', p.name))]


def output_times(case: Path) -> dict:
    """Each 3-D output the case holds, by model time (s, rounded): its file and control file, from the metadata file
    of each output configuration (cm1out_metadata.dat, .02 ...), which numbers the outputs from 1. A run on from a
    restart copies its source's metadata, whose outputs it does not hold, and they are left out."""
    out = {}
    for suffix in configurations(case, 'cm1out_metadata'):
        names = [n for n, _ in ra.grads_layout(case / f'cm1out_metadata{suffix}.ctl')[1]]
        table = np.fromfile(case / f'cm1out_metadata{suffix}.dat', dtype='<f4')
        table = table[:table.size // len(names) * len(names)].reshape(-1, len(names))
        for n, t in zip(table[:, names.index('nwrite')], table[:, names.index('mtime')]):
            path = case / f'cm1out_t{int(round(float(n))):06d}_s{suffix}.dat'
            if n > 0 and path.exists():
                out[round(float(t))] = (path, case / f'cm1out_s{suffix}.ctl')
    return out


def read_snapshot(path: Path, ctl: Path) -> dict:
    """One 3-D output by its control file's layout, as ring_analysis.read_snapshot reads it: sections (nz, columns)."""
    nx, layout = ra.grads_layout(ctl)
    raw = np.fromfile(path, dtype='<f4')
    out, offset = {}, 0
    for name, levels in layout:
        field = raw[offset:offset + levels * nx].reshape(levels, nx)
        out[name] = field[0] if levels == 1 else field
        offset += levels * nx
    if offset != raw.size:
        raise ValueError(f'{path.name}: {raw.size} values, {ctl.name} expects {offset}')
    return out


def stats_rows(case: Path) -> dict:
    """CM1's statistics over every output configuration (cm1out_stats.dat, .02 ...), as arrays in time order, the
    records a restart leaves unwritten (time nil) left out."""
    parts = []
    for suffix in configurations(case, 'cm1out_stats'):
        lines = (case / f'cm1out_stats{suffix}.ctl').read_text().splitlines()
        i = next(k for k, line in enumerate(lines) if line.lower().startswith('vars'))
        names = [line.split()[0] for line in lines[i + 1:i + 1 + int(lines[i].split()[1])]]
        raw = np.fromfile(case / f'cm1out_stats{suffix}.dat', dtype='<f4')
        data = raw[:raw.size // len(names) * len(names)].reshape(-1, len(names)).astype(float)
        parts.append({n: data[:, k] for k, n in enumerate(names)})
    keys = set.intersection(*(set(p) for p in parts))
    out = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    keep = out['mtime'] > 0.0
    order = np.argsort(out['mtime'][keep], kind='stable')
    return {k: v[keep][order] for k, v in out.items()}


def rain_by_quarter(name: str, span_days) -> dict:
    """Rain at the ground (mm/day) by quarter of the lunar day over a span, from CM1's running domain total (train),
    each interval between statistics binned by its midpoint's hour angle, as massflux.cm1_microphysics bins it."""
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    cfg, grid = record['configuration'], record['grid']
    s = stats_rows(case)
    area = grid['nx'] * grid.get('ny', 1) * grid['dx_m'] ** 2
    t = s['mtime']
    idx = np.flatnonzero((t > span_days[0] * DAY) & (t <= span_days[1] * DAY + 1.0))
    acc = {q[0]: [0.0, 0.0] for q in bl.QUARTERS}
    for a, b in zip(idx[:-1], idx[1:]):
        hour = local_hour_angle(cfg['start_hour_angle_deg'], cfg['site']['lon_deg'], 0.5 * (t[a] + t[b]), solar_day_s())
        acc[bl.quarter_name(hour)][0] += float(s['train'][b] - s['train'][a]) / area
        acc[bl.quarter_name(hour)][1] += t[b] - t[a]
    return {q: v[0] / v[1] * DAY if v[1] > 0.0 else None for q, v in acc.items()}


def cloud_tops(qc, qi, zh) -> dict:
    """The cloud of one snapshot column by column ((nz, columns) arrays): the share of columns holding cloud below
    TOP_M, the median and 90th percentile of their highest cloudy level below it, and the share of all columns whose
    cloud below it reaches above DEEP_M."""
    zh = np.asarray(zh, float)
    below = zh < TOP_M
    cloudy = ((np.asarray(qc) + np.asarray(qi)) >= ra.CLOUD_KG_KG) & below[:, None]
    has = cloudy.any(axis=0)
    tops = zh[np.where(cloudy, np.arange(zh.size)[:, None], 0).max(axis=0)][has]
    return dict(cloudy_share=float(has.mean()), top_median_m=float(np.median(tops)) if has.any() else None,
                top_p90_m=float(np.percentile(tops, 90)) if has.any() else None,
                above_16km_share=float((tops > DEEP_M).sum() / has.size))


def diagnostics(name: str, first_s: float, last_s: float) -> list:
    """Every domain diagnostic of a case after first_s up to last_s: its time, hour angle, virtual potential temperature
    and vapour (g/kg) profiles, hpbl and highest cloudy level."""
    from climate.crm.box_resolution import cloud_top_m
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    cfg = record['configuration']
    zh = ra.case_geometry(case)['zh']
    held = record.get('held')
    rows = []
    for path, _, layout, _ in bl.diag_files(case):
        d = bb.read_diag(case, path, layout)
        if not first_s + 1.0 < d['mtime'] <= last_s + 1.0:
            continue
        hour = held['hour_angle_deg'] if held else local_hour_angle(cfg['start_hour_angle_deg'], cfg['site']['lon_deg'],
                                                                     d['mtime'], solar_day_s())
        rows.append(dict(mtime=d['mtime'], hour=hour, theta_v=d['th'] * (1.0 + (co.RV / co.RD - 1.0) * d['qv']),
                         q_g_kg=d['qv'] * 1000.0, hpbl_m=float(d['hpbl']),
                         cloud_top_m=cloud_top_m(d['qcfrac'], d['qifrac'], zh)))
    return sorted(rows, key=lambda r: r['mtime'])


def profile_measures(zh, theta_v, q_g_kg) -> dict:
    """The comparison's measures on mean profiles at CM1's scalar levels (ground first)."""
    from climate.gcm import pbl_mixing as pm
    zh = np.asarray(zh, float)
    at = lambda v, z: float(np.interp(z, zh, v))
    depth = pm.rise_height(zh, theta_v, Z_LOW_M)
    return dict(mixed_depth_km=None if depth is None else depth / 1000.0,
                vapour_drop_g_kg=at(q_g_kg, Z_LOW_M) - at(q_g_kg, Z_MID_M),
                lid_rise_k=at(theta_v, Z_LID_M) - at(theta_v, Z_MID_M),
                vapour_g_kg={f'{z / 1000:g} km': at(q_g_kg, z) for z in PROFILE_M},
                theta_v_k={f'{z / 1000:g} km': at(theta_v, z) for z in PROFILE_M})


def summarize(rows, zh) -> dict:
    """Means over a group of diagnostics: the measures on the mean profiles, hpbl and the highest cloudy level."""
    mean = lambda key: np.mean([r[key] for r in rows], axis=0)
    tops = [r['cloud_top_m'] for r in rows if r['cloud_top_m'] is not None]
    return dict(outputs=len(rows), **profile_measures(zh, mean('theta_v'), mean('q_g_kg')),
                hpbl_km=float(mean('hpbl_m')) / 1000.0,
                highest_cloud_km=float(np.mean(tops)) / 1000.0 if tops else None)


def snapshots(name: str, times_s, zh) -> dict:
    """cloud_tops of a case's 3-D outputs at the given model times, by time (rounded seconds)."""
    held = output_times(ra.RUNS / name)
    out = {}
    for t in times_s:
        if round(t) in held:
            d = read_snapshot(*held[round(t)])
            out[round(t)] = cloud_tops(d['qc'], d['qi'], zh)
    return out


def snapshot_means(shots: list) -> dict:
    if not shots:
        return None
    mean = lambda key: (lambda v: float(np.mean(v)) if v else None)([s[key] for s in shots if s[key] is not None])
    return dict(snapshots=len(shots), **{k: mean(k) for k in ('cloudy_share', 'top_median_m', 'top_p90_m',
                                                               'above_16km_share')})


def lid_mean(levels) -> dict:
    """box_lid's parts averaged over the lid levels."""
    if not levels:
        return None
    out = dict(heights_km=[lev['height_m'] / 1000.0 for lev in levels],
               khv_m2_s=float(np.mean([lev['khv_m2_s'] for lev in levels])))
    for q in ('heat', 'vapour'):
        out[q] = {k: float(np.mean([lev[q][k] for lev in levels])) for k in ('resolved', 'diffusion', 'subgrid', 'total')}
        out[q]['subgrid_share'] = out[q]['subgrid'] / out[q]['total'] if out[q]['total'] else None
    return out


def compare_days() -> dict:
    """Both runs by quarter of the lunar day over SPAN_DAYS."""
    first_s, last_s = SPAN_DAYS[0] * DAY, SPAN_DAYS[1] * DAY
    zh = ra.case_geometry(ra.RUNS / RUNS['stretched'])['zh']
    rows = {k: diagnostics(name, first_s, last_s) for k, name in RUNS.items()}
    common = sorted(set(round(r['mtime']) for r in rows['earth']) & set(round(r['mtime']) for r in rows['stretched']))
    shot_times = sorted(t for t in output_times(ra.RUNS / RUNS['earth']) if first_s + 1.0 < t <= last_s + 1.0)
    shots = {k: snapshots(name, shot_times, zh) for k, name in RUNS.items()}
    lids = {k: bl.lid_split(name, bl.LID_M, first_s, last_s, by_quarter=True) for k, name in RUNS.items()}
    rains = {k: rain_by_quarter(name, SPAN_DAYS) for k, name in RUNS.items()}
    quarters = {}
    for qname, _, _ in bl.QUARTERS:
        entry = {}
        for k in RUNS:
            group = [r for r in rows[k] if round(r['mtime']) in common and bl.quarter_name(r['hour']) == qname]
            hours = {round(r['mtime']): bl.quarter_name(r['hour']) for r in rows[k]}
            entry[k] = dict(summarize(group, zh), rain_mm_day=rains[k][qname],
                            clouds=snapshot_means([s for t, s in shots[k].items()
                                                   if t in shots['earth'] and t in shots['stretched']
                                                   and hours.get(t) == qname]),
                            lid=lid_mean(lids[k]['quarters'][qname]))
        quarters[qname] = entry
    night = quarters['midnight-sunrise']
    reading = dict(depth_km=night['earth']['mixed_depth_km'], drop_g_kg=night['earth']['vapour_drop_g_kg'],
                   stretched_depth_km=night['stretched']['mixed_depth_km'],
                   stretched_drop_g_kg=night['stretched']['vapour_drop_g_kg'])
    reading['verdict'] = verdict(reading['depth_km'], reading['drop_g_kg'])
    return dict(span_days=list(SPAN_DAYS), common_outputs=len(common), snapshots_compared=len(shot_times),
                l_inf_m={k: lids[k]['l_inf_m'] for k in RUNS}, quarters=quarters, after_midnight=reading)


def compare_nights() -> dict | None:
    """The held nights at 6 km, both lengths, over their second day; None until the Earth-length night has run."""
    if not (ra.RUNS / NIGHTS['earth'] / 'progress.json').exists():
        return None
    first_s, last_s = bl.NIGHT_SPAN_H[0] * HOUR, bl.NIGHT_SPAN_H[1] * HOUR
    zh = ra.case_geometry(ra.RUNS / NIGHTS['stretched'])['zh']
    out = {}
    for k, name in NIGHTS.items():
        rows = diagnostics(name, first_s, last_s)
        out[k] = dict(summarize(rows, zh), lid=lid_mean(bl.lid_split(name, bl.LID_M, first_s, last_s)['levels']),
                      rain_mm_day=rain_over_span(name, first_s, last_s))
    return dict(span_hours=list(bl.NIGHT_SPAN_H), runs=out)


def rain_over_span(name: str, first_s: float, last_s: float) -> float:
    """Rain at the ground (mm/day) between two model times, from CM1's running domain total."""
    case = ra.RUNS / name
    s = stats_rows(case)
    grid = json.loads((case / 'case.json').read_text())['grid']
    area = grid['nx'] * grid.get('ny', 1) * grid['dx_m'] ** 2
    t, total = s['mtime'], s['train'] / area
    return float((np.interp(last_s, t, total) - np.interp(first_s, t, total)) / ((last_s - first_s) / DAY))


def run_record(name: str) -> dict | None:
    """What a run cost: its wall time, threads and the size of what it wrote."""
    case = ra.RUNS / name
    if not (case / 'progress.json').exists():
        return None
    progress = json.loads((case / 'progress.json').read_text())
    record = json.loads((case / 'case.json').read_text())
    written = sum(p.stat().st_size for p in case.iterdir() if p.is_file() and not p.is_symlink())
    return dict(case=name, runner=record.get('runner'), executable=record['build'].get('executable_sha256'),
                segments=len(progress['segments']),
                wall_hours=round(sum(s['wall_s'] for s in progress['segments']) / HOUR, 3),
                threads=sorted({s['threads'] for s in progress['segments']}), folder_mb=round(written / 1e6, 1))


def analyse() -> dict:
    return dict(schema=SCHEMA, runs=RUNS, criteria=CRITERIA, days=compare_days(), held_nights=compare_nights(),
                records={k: run_record(n) for k, n in {**RUNS, **{f'night_{k}': v for k, v in NIGHTS.items()}}.items()},
                evidence=('Two CM1 runs of the sea box\'s third lunar day from one restart, one realization each, which '
                          'differ in the boundary-layer scheme\'s asymptotic length alone (453 and 75 m); domain '
                          'diagnostics every 3 hours compared output for output, 3-D outputs every 6 hours, CM1\'s own '
                          'flux diagnostics split by part. A sensitivity of CM1 to its own scheme on a 6 km grid, not a '
                          'test of either length against observation.'),
                reading_rule=('days.quarters: by quarter of the lunar day (hour angle from noon), each run\'s means over '
                              'the domain diagnostics both wrote: mixed_depth_km by pbl_mixing.rise_height from 0.314 km '
                              '(0.5 K of virtual potential temperature), vapour_drop_g_kg the mixing ratio at 0.314 km '
                              'less that at 6.201 km, lid_rise_k the virtual potential temperature at 8.238 km less that '
                              'at 6.201 km (the GCM\'s full levels at the site), vapour_g_kg and theta_v_k at fixed '
                              'heights, hpbl_km CM1\'s own depth, highest_cloud_km the highest level any column holds '
                              'cloud; rain_mm_day at the ground from CM1\'s running totals; clouds: means over the 3-D '
                              'outputs both wrote (6-hourly) of the share of columns holding cloud of 0.01 g/kg below '
                              '40 km, the median and 90th percentile of their highest cloudy level (m) and the share of '
                              'all columns whose cloud reaches above 16 km; lid: box_lid\'s parts (W/m2, upward) averaged '
                              'over 6.45 and 6.99 km, with the scheme\'s share and diffusivity. after_midnight: the '
                              'verdict under criteria, set before the run. held_nights: the 6 km held nights\' second '
                              'day, both lengths.'),
                producer=dict(domain='climate', files={f'crm/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                                                       for p in (Path(__file__), HERE / 'box_lid.py',
                                                                 HERE / 'box_budget.py')}))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args(argv)
    out = analyse()
    args.output.write_text(json.dumps(bl.compact(out), indent=1) + '\n')
    fmt = lambda v, f='.2f': '-' if v is None else format(v, f)
    for qname, e in out['days']['quarters'].items():
        for k in RUNS:
            r = e[k]
            lid = r['lid'] or {}
            print(f"{qname:17s} {k:9s} depth {fmt(r['mixed_depth_km'])} km, drop {fmt(r['vapour_drop_g_kg'])} g/kg, lid "
                  f"{fmt(r['lid_rise_k'])} K, hpbl {fmt(r['hpbl_km'])} km, rain {fmt(r['rain_mm_day'])} mm/day, "
                  f"heat sub {fmt(lid.get('heat', {}).get('subgrid'), '.1f')} of {fmt(lid.get('heat', {}).get('total'), '.1f')}, "
                  f"vapour sub {fmt(lid.get('vapour', {}).get('subgrid'), '.1f')} of "
                  f"{fmt(lid.get('vapour', {}).get('total'), '.1f')} W/m2, khv {fmt(lid.get('khv_m2_s'), '.1f')}")
    print('after midnight:', out['days']['after_midnight'])
    print(f'wrote {args.output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
