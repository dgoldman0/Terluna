"""What an electrified CM1 storm shows: how often and where it flashes, the field it reaches, the charge it carries
and which collisions charged it. Stage 2 of the atmospheric-electricity study (climate/crm/cm1_elec.py).

    python3 -m climate.crm.elec_analysis supercell_elec_msz

writes ../results/crm/elec_<case>.json from:
- the flash log (terluna_flashes.txt): each flash of WRF-ELEC's branched scheme (in cloud, or negative or positive
  charge to ground) with its starting point and the field there, the charge it neutralized, the electrostatic energy
  before and after, the levels and area its channels reached and the nitrogen oxides it made; or each call of the
  cylindrical scheme, with the points above breakdown, the columns and separate regions it took and the charge it
  removed. Each row also carries the domain's largest field over breakdown just before it;
- the field log (terluna_field.txt, every ten steps): the largest field, the domain's positive and negative charge and
  its electrostatic energy, and the largest charging rates;
- the charging totals WRF-ELEC's driver prints every step into CM1's log (C/s over the domain, by collision pair:
  graupel and hail with cloud ice and with snow, snow with ice, and inductive charging of graupel by cloud droplets);
- the charge tracers in the output snapshots (C/kg; times the air density, C/m3): the charge by height and particle
  type, and where the main negative and the positive charge regions sit in height and temperature.

A branched flash is one flash as a lightning mapping array would count it, though its channels run on the model's
grid. A cylindrical discharge acts at once in every cylinder around every point above the breakdown field, so a call
that takes several separate regions counts them in its regions column and stands for several flashes.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path
import sys
import numpy as np
from climate.crm import ring_analysis as ra
from climate.crm.mixed_phase_analysis import significant, stats

RESULTS = ra.RESULTS
SCHEMA = 'terluna.climate.crm-electricity/2'                 # 2: series per bin of bin_s (2026-10-05)
RD, CP, P00 = 287.04, 1005.7, 1.0e5
SPECIES = ('cloud_water', 'rain', 'cloud_ice', 'snow', 'graupel', 'small_ions', 'hail')   # passive tracers 1-7
FLASH_COLUMNS = ('time_s', 'step', 'substep', 'kind', 'x_m', 'y_m', 'z_m', 'e_v_m', 'e_break_v_m', 'e_over_break_max',
                 'extent_m2', 'regions', 'positive_c', 'negative_c', 'energy_before_j', 'energy_after_j', 'z_low_m',
                 'z_high_m', 'points', 'nox_mol')
FIELD_COLUMNS = ('time_s', 'step', 'e_max_v_m', 'x_m', 'y_m', 'z_m', 'positive_c', 'negative_c', 'energy_j',
                 'noninductive_max_c_m3_s', 'inductive_max_c_m3_s')
KINDS = {1: 'in_cloud', 2: 'negative_to_ground', 3: 'positive_to_ground', 9: 'cylinders'}
# the cloud ice the radiation reads with its size capped at 140 microns (terluna_ice_optics.txt, runs with radiation)
ICE_COLUMNS = ('time_s', 'step', 'ice_kg', 'ice_above_140um_kg', 'mean_radius_above_um', 'max_radius_um', 'cells_above')
# the field at the ground and point discharge (terluna_ground.txt): the charge the ground's points gave off since the
# line before (where a case turns point discharge on), the largest field at the ground before discharge since then
# over all land and where no particles reach the lowest level, and the columns discharging, all and wet
GROUND_COLUMNS = ('time_s', 'step', 'positive_c', 'negative_c', 'e_ground_max_v_m', 'e_ground_dry_max_v_m',
                  'columns', 'columns_wet')
# WRF-ELEC's charging totals (module_mp_nssl_2mom.F, its driver): negative and positive parts of the domain's charging
# rate, C/s, by the collisions that separate it
CHARGING = {'ctghi': 'graupel and hail with cloud ice', 'ctghs': 'graupel and hail with snow',
            'ctswi': 'snow with cloud ice', 'ctghw': 'graupel with cloud droplets (inductive)',
            'ctgi': 'graupel with cloud ice', 'ctgs': 'graupel with snow', 'cthi': 'hail with cloud ice',
            'cths': 'hail with snow', 'ctsww': 'snow with cloud droplets'}
TIME_UNITS = {'sec': 1.0, 'min': 60.0, 'hour': 3600.0, 'hrs': 3600.0, 'day': 86400.0}
BIN_S = 300.0
LONG_RUN_S = 2 * 86400.0          # longer runs are binned by their output interval, so the summary stays small
EVIDENCE = ('CM1 r22.0 with WRF-ELEC\'s NSSL two-moment microphysics (MicroTed/wrf4-elec e43041b) and its charging '
            '(non-inductive by collisions of graupel and hail with ice and snow, inductive by droplets rebounding from '
            'graupel in the field), the small ions\' net charge attaching to particles, an FFT field solver, and '
            'WRF-ELEC\'s lightning on its own sub-steps of the sedimentation: the branched scheme (lightmsz, '
            'MacGorman, Straka and Ziegler 2001) or the cylindrical one. Charge is carried by bulk particle '
            'populations; branched channels run on the model grid and neutralize the charge where they reach; the '
            'ground-strike rule, the 50-MV potential and 10-kV/m field a ground flash needs at its start and the '
            'nitrogen oxide yield are WRF-ELEC\'s calibrations for Earth; nothing conducts charge away except where a '
            'case turns on leakage, and the ground gives off ions only where a case turns on point discharge '
            '(Standler and Winn\'s 1 nA/m2 at 8 kV/m over Earth\'s vegetation, its onset scaled by density).')
READING_RULE = ('Times are model seconds from the start. Series are per bin of bin_s seconds from start_s: 5 minutes, or the output interval for runs longer than two days. Kinds: 1 in cloud, 2 negative and 3 positive charge to ground '
                '(branched flashes), 9 a call of the cylindrical scheme, whose regions are the separate areas one call '
                'takes. A flash\'s positive and negative charge are what it neutralized; a ground flash neutralizes '
                'one sign only, the ground supplying the other. Nitrogen oxides are lightmsz\'s: Wang et al.\'s yield '
                'per metre of channel at its pressure, times 0.1 where the flash changes the charge by 1 nC/kg or '
                'more. Charging rates are the domain\'s totals of each sign (C/s). Charge densities are nC/m3; '
                'heights are above the ground; temperatures are the air\'s at the level, weighted by the charge '
                'there.')


def read_log(path: Path, names: tuple) -> dict:
    """A Terluna log as columns. Each run of the model marks the time it takes up from ('# run from'), and the rows an
    earlier run wrote after that time, which the restart replaces, are dropped."""
    rows = []
    if path.exists():
        for line in path.read_text().splitlines():
            if line.startswith('# run from'):
                start = float(line.split()[-1])
                rows = [r for r in rows if r[0] <= start]
            elif line.strip() and not line.startswith('#'):
                rows.append([float(v) for v in line.split()])
    rows = np.array(rows, dtype=float).reshape(-1, len(names))
    return {name: rows[:, k] for k, name in enumerate(names)}


def segment_logs(case: Path) -> list:
    logs = sorted(case.glob('cm1_segment_*.log'))
    return logs if logs else [p for p in (case / 'cm1.log',) if p.exists()]


def charging(case: Path, logs: list | None = None) -> dict:
    """WRF-ELEC's charging totals by step from CM1's logs (all of a case's, or the ones given), each assigned to the
    step whose time line follows it."""
    pattern = re.compile(r'^(ct[a-z]+)n,\1p\s*=\s*(\S+),\s*(\S+)')
    step = re.compile(r'^\s+(\d+)\s+([0-9.Ee+-]+)\s+(sec|min|hour|hrs|day)')
    pending, rows = {}, {}
    for log in segment_logs(case) if logs is None else logs:
        for line in log.read_text(errors='replace').splitlines():
            m = pattern.match(line)
            if m:
                pending[m.group(1)] = (float(m.group(2)), float(m.group(3)))
                continue
            m = step.match(line)
            if m and pending:
                rows[int(m.group(1))] = (float(m.group(2)) * TIME_UNITS[m.group(3)], pending)
                pending = {}
    steps = sorted(rows)
    out = dict(step=np.array(steps), time_s=np.array([rows[s][0] for s in steps]))
    for name in CHARGING:
        out[name] = np.array([rows[s][1].get(name, (np.nan, np.nan)) for s in steps]).reshape(-1, 2)
    return out


def heights_km(ctl: Path) -> np.ndarray:
    """The scalar levels of a GrADS description (km)."""
    lines = ctl.read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith('zdef'))
    parts = lines[i].split()
    n = int(parts[1])
    if parts[2] == 'linear':
        return float(parts[3]) + float(parts[4]) * np.arange(n)
    values = [float(v) for v in parts[3:]]
    j = i + 1
    while len(values) < n:
        values += [float(v) for v in lines[j].split()]
        j += 1
    return np.array(values[:n])


def snapshot_charge(case: Path, n: int, dx: float, dy: float, zw: np.ndarray) -> dict:
    """One output time's charge (C/m3 by tracer, net), air temperature (K) and density (kg/m3), the vertical wind and
    the reflectivity, as (nz, ncol) arrays."""
    snap = ra.read_snapshot(case, n)
    p, th, qv = snap['prs'], snap['th'], snap['qv']
    t = th * (p / P00) ** (RD / CP)
    rho = p / (RD * t * (1.0 + 0.608 * qv))
    names = sorted(k for k in snap if re.fullmatch(r'pt\d+', k))
    charges = {SPECIES[int(k[2:]) - 1]: snap[k] * rho for k in names}
    net = sum(charges.values())
    volume = dx * dy * np.diff(zw)[:, None]
    graupel = sum(snap[q] for q in ('qg', 'qhl') if q in snap) * rho if 'qg' in snap else None
    return dict(charges=charges, net=net, t=t, rho=rho, volume=volume, w=snap.get('winterp'), dbz=snap.get('dbz'),
                graupel=graupel)


def structure(s: dict, zh: np.ndarray, threshold: float = 0.1e-9) -> dict:
    """Where the charge sits: by height the positive and negative charge (C) of the net and of each particle type, and
    the main negative region (the level holding the most negative charge) with the positive regions above and below
    it, each with its height, the temperature weighted by the charge there and the largest density."""
    net, vol = s['net'], s['volume']
    pos = np.where(net > threshold, net, 0.0) * vol
    neg = np.where(net < -threshold, net, 0.0) * vol
    pos_z, neg_z = pos.sum(axis=1), neg.sum(axis=1)
    out = dict(positive_c=float(pos.sum()), negative_c=float(neg.sum()), max_nc_m3=float(net.max() * 1e9),
               min_nc_m3=float(net.min() * 1e9),
               species_max_abs_nc_m3={k: float(np.abs(v).max() * 1e9) for k, v in s['charges'].items()})
    if neg_z.min() >= 0.0:
        return out
    k = int(np.argmin(neg_z))

    def region(level, weights):
        w = np.abs(weights[level])
        return dict(z_km=float(zh[level]), t_c=float(np.sum(s['t'][level] * w) / w.sum() - 273.15),
                    charge_c=float(weights[level].sum()),
                    density_nc_m3=float(net[level][np.argmax(np.abs(net[level]))] * 1e9))

    out['main_negative'] = region(k, neg)
    above = np.arange(k + 1, len(zh))
    below = np.arange(0, k)
    if above.size and pos_z[above].max() > 0.0:
        out['upper_positive'] = region(int(above[np.argmax(pos_z[above])]), pos)
    if below.size and pos_z[below].max() > 0.0:
        out['lower_positive'] = region(int(below[np.argmax(pos_z[below])]), pos)
    # who carries the main negative region's charge
    carried = {name: float((np.where(net[k] < -threshold, q[k], 0.0) * vol[k]).sum())
               for name, q in s['charges'].items()}
    out['main_negative']['carried_c'] = carried
    return out


def binned(times: np.ndarray, values: np.ndarray, edges: np.ndarray, how: str = 'sum') -> list:
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (times > a) & (times <= b)
        if how == 'sum':
            out.append(float(values[sel].sum()))
        elif how == 'max':
            out.append(float(values[sel].max()) if sel.any() else 0.0)
        else:
            out.append(float(values[sel].mean()) if sel.any() else float('nan'))
    return out


def ice_cap_summary(ice: dict) -> dict | None:
    """The ice the 140-micron cap cuts from the radiation's view: the share of the domain's cloud ice above it (largest,
    and over all logged steps by mass), the mass-weighted size of that ice and the largest size where ice is 0.001
    g/kg or more. For the radiation and cloud-optics work: that ice reaches the radiation optically thicker than its
    size gives, by its size over 140 microns."""
    if not ice['time_s'].size or not np.any(ice['ice_kg'] > 0.0):
        return None
    has = ice['ice_kg'] > 0.0
    above = ice['ice_above_140um_kg']
    return dict(rows=int(has.sum()), share_above_max=float(np.max(above[has] / ice['ice_kg'][has])),
                share_above_overall=float(above.sum() / ice['ice_kg'].sum()),
                mean_radius_above_um=float(np.sum(above * ice['mean_radius_above_um']) / above.sum())
                if above.sum() > 0 else None,
                max_radius_um=float(ice['max_radius_um'].max()), optical_depth_overstated=(
                    float(np.sum(above * ice['mean_radius_above_um']) / above.sum() / 140.0) if above.sum() > 0 else None))


def ground_summary(g: dict) -> dict | None:
    """The field at the ground over land, before any point discharge: largest and median over the logged intervals,
    over all land and where no particles reach the lowest level; and the charge the ground's points gave off of each
    sign, in all and at most in one logged interval, with the columns discharging (all and wet) at most."""
    if not g['time_s'].size:
        return None
    return dict(rows=int(g['time_s'].size), e_ground_max_kv_m=float(g['e_ground_max_v_m'].max() / 1e3),
                e_ground_median_kv_m=float(np.median(g['e_ground_max_v_m']) / 1e3),
                e_ground_dry_max_kv_m=float(g['e_ground_dry_max_v_m'].max() / 1e3),
                e_ground_dry_median_kv_m=float(np.median(g['e_ground_dry_max_v_m']) / 1e3),
                positive_c=float(g['positive_c'].sum()), negative_c=float(g['negative_c'].sum()),
                positive_c_max=float(g['positive_c'].max()), negative_c_min=float(g['negative_c'].min()),
                columns_max=int(g['columns'].max()), columns_wet_max=int(g['columns_wet'].max()))


def flash_summary(f: dict, edges: np.ndarray) -> dict:
    """Flashes of one kind: how many and when, the charge each neutralized, the energy it released, where it started
    and how far its channels reached, and the nitrogen oxides it made."""
    n = f['time_s'].size
    dissipated = f['energy_before_j'] - f['energy_after_j']
    out = dict(count=int(n), first_s=float(f['time_s'].min()), per_bin=binned(f['time_s'], np.ones(n), edges),
               positive_c=stats(f['positive_c']), negative_c=stats(f['negative_c']),
               energy_dissipated_j=stats(dissipated), initiation_z_km=stats(f['z_m'] / 1000.0),
               e_at_initiation_kv_m=stats(f['e_v_m'] / 1000.0),
               e_break_at_initiation_kv_m=stats(f['e_break_v_m'] / 1000.0), extent_km2=stats(f['extent_m2'] / 1.0e6),
               z_low_km=stats(f['z_low_m'] / 1000.0), z_high_km=stats(f['z_high_m'] / 1000.0),
               points=stats(f['points']))
    if np.any(f['nox_mol'] > 0.0):
        out.update(nox_mol=stats(f['nox_mol']), nox_total_mol=float(f['nox_mol'].sum()))
    if np.any(f['regions'] != 1):
        out['regions'] = int(f['regions'].sum())
    return out


def analyse(name: str) -> dict:
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    cfg = record['configuration']
    grid = record['grid']
    dx, dy = grid['dx_m'], grid.get('dy_m', grid['dx_m'])
    zh = heights_km(case / 'cm1out_s.ctl')
    zw = (np.loadtxt(case / 'input_grid_z') if (case / 'input_grid_z').exists()
          else np.concatenate([[0.0], 1000.0 * (zh[:-1] + zh[1:]) / 2.0, [1000.0 * (2 * zh[-1] - (zh[-1] + zh[-2]) / 2.0)]]))
    flashes = read_log(case / 'terluna_flashes.txt', FLASH_COLUMNS)
    field = read_log(case / 'terluna_field.txt', FIELD_COLUMNS)
    ice = read_log(case / 'terluna_ice_optics.txt', ICE_COLUMNS)
    ground = read_log(case / 'terluna_ground.txt', GROUND_COLUMNS)
    rates = charging(case)
    end = float(max(field['time_s'].max() if field['time_s'].size else 0.0,
                    rates['time_s'].max() if rates['time_s'].size else 0.0))
    start = float(json.loads((case / 'progress.json').read_text()).get('start_s', 0.0)) \
        if (case / 'progress.json').exists() else 0.0
    bin_s = BIN_S if end - start <= LONG_RUN_S else cfg['output_s'] * record.get('time_scale', 1.0)
    edges = np.arange(start, end + bin_s, bin_s)

    lightning = {}
    first = np.r_[True, (np.diff(flashes['step']) != 0) | (np.diff(flashes['substep']) != 0)] if flashes['step'].size \
        else np.array([], bool)
    for code, label in KINDS.items():
        sel = flashes['kind'] == code
        if sel.any():
            lightning[label] = flash_summary({k: v[sel] for k, v in flashes.items()}, edges)
    if flashes['step'].size:
        groups = np.cumsum(first)
        over = flashes['e_over_break_max']
        lightning['all'] = dict(
            count=int(flashes['step'].size), first_s=float(flashes['time_s'].min()),
            per_bin=binned(flashes['time_s'], np.ones_like(flashes['time_s']), edges),
            # the field over breakdown before a sub-step's first flash (what charging and transport built since the
            # last) and before its later ones (what the flashes before them left)
            over_breakdown_first=stats(over[first]), over_breakdown_later=stats(over[~first]),
            per_substep_max=int(np.bincount(groups).max()))
    fields = dict(
        e_max_kv_m=float(field['e_max_v_m'].max() / 1000.0) if field['time_s'].size else None,
        e_max_per_bin_kv_m=binned(field['time_s'], field['e_max_v_m'] / 1000.0, edges, 'max'),
        positive_c_per_bin=binned(field['time_s'], field['positive_c'], edges, 'max'),
        negative_c_per_bin=binned(field['time_s'], -field['negative_c'], edges, 'max'),
        energy_j_per_bin=binned(field['time_s'], field['energy_j'], edges, 'max'),
        noninductive_max_pc_m3_s=float(field['noninductive_max_c_m3_s'].max() * 1e12) if field['time_s'].size else None,
        inductive_max_pc_m3_s=float(field['inductive_max_c_m3_s'].max() * 1e12) if field['time_s'].size else None,
        noninductive_max_per_bin_pc_m3_s=binned(field['time_s'], field['noninductive_max_c_m3_s'] * 1e12, edges, 'max'),
        inductive_max_per_bin_pc_m3_s=binned(field['time_s'], field['inductive_max_c_m3_s'] * 1e12, edges, 'max'))
    dt = np.diff(np.concatenate([[0.0], rates['time_s']])) if rates['time_s'].size else np.array([])
    separated = {}
    for key, label in CHARGING.items():
        r = rates[key]
        if not r.size or np.all(np.isnan(r)):
            continue
        separated[key] = dict(collisions=label, positive_c=float(np.nansum(r[:, 1] * dt)),
                              negative_c=float(np.nansum(r[:, 0] * dt)),
                              peak_c_s=float(np.nanmax(np.abs(r))),
                              per_bin_c=binned(rates['time_s'], (r[:, 1] - r[:, 0]) * dt, edges))
    snaps = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    structures = []
    for n in snaps:
        s = snapshot_charge(case, n, dx, dy, zw)
        st = structure(s, zh)
        st.update(time_s=(n - 1) * cfg['output_s'], w_max_m_s=float(np.nanmax(s['w'])) if s['w'] is not None else None,
                  dbz_max=float(np.nanmax(s['dbz'])) if s['dbz'] is not None else None)
        structures.append(st)
    return dict(
        schema=SCHEMA, case=name, purpose=cfg.get('purpose'), electricity=record.get('electricity'),
        build=dict(label=record['build']['label'], gravity_m_s2=record['build']['gravity_m_s2'],
                   executable_sha256=record['build']['executable_sha256']),
        producer=dict(domain='climate', files={'crm/elec_analysis.py': hashlib.sha256(
            Path(__file__).read_bytes()).hexdigest()[:16]}),
        evidence=EVIDENCE, reading_rule=READING_RULE, start_s=start, end_s=end, bin_s=bin_s,
        lightning=lightning, field=fields, charging=separated, structure=structures, ice_cap=ice_cap_summary(ice),
        ground=ground_summary(ground))


def against(name: str, other: str) -> dict:
    """The storm beside another run of the same case: the largest updraft, reflectivity, graupel, hail, cloud ice and
    its number, and accumulated rain at each shared snapshot, and the condensate by height as a ratio."""
    a, b = ra.RUNS / name, ra.RUNS / other
    shared = sorted({int(p.name[8:14]) for p in a.glob('cm1out_t*_s.dat')} &
                    {int(p.name[8:14]) for p in b.glob('cm1out_t*_s.dat')})
    output_s = json.loads((a / 'case.json').read_text())['configuration']['output_s']
    fields = dict(w_m_s=('winterp', 1.0), dbz=('dbz', 1.0), graupel_g_kg=('qg', 1e3), hail_g_kg=('qhl', 1e3),
                  cloud_ice_g_kg=('qi', 1e3), cloud_ice_per_kg=('cci', 1.0), rain_cm=('rain', 1.0))
    rows = []
    for n in shared:
        sa, sb = ra.read_snapshot(a, n), ra.read_snapshot(b, n)
        row = dict(time_s=(n - 1) * output_s)
        for key, (var, scale) in fields.items():
            row[key] = [float(np.nanmax(sa[var]) * scale), float(np.nanmax(sb[var]) * scale)]
        condensate = [sum(s[q] for q in ('qc', 'qr', 'qi', 'qs', 'qg', 'qhl')).sum(axis=1) for s in (sa, sb)]
        row['condensate_ratio_by_level'] = [float(x / y) if y > 0 else None for x, y in zip(*condensate)]
        rows.append(row)
    return dict(case=name, other=other, columns='each pair is [this run, the other run]', snapshots=rows)


# Windows of box_0e_elec's storms run again with one setting changed, beside the runs of the same days without the
# change: Takahashi's charging law (WRF-ELEC's isaund 1) for Saunders and Peck's (isaund 12), and leakage through the
# air's conductivity. CM1's threads make every run of a window a different realization of its storms, so each window has
# more than one run without the change to show how far chance alone moves a number; each case's own settings are
# recorded with it.
SAUNDERS_PECK_RUNS = {'first_lunar_day': ('box_0e_elec_uncapped', 'box_0e_elec_uncapped_corona', 'box_0e_elec'),
                      'second_lunar_day': ('box_0e_elec_ground_rule', 'box_0e_elec')}
WINDOW_DAYS = {'first_lunar_day': (10.5, 12.0), 'second_lunar_day': (40.5, 42.0)}
# each comparison's runs by window, those without the change first; the leader's crossing is compared with leakage on,
# beside the leakage window with WRF-ELEC's rule alone
COMPARISONS = {
    'charging_laws': dict(first_lunar_day=SAUNDERS_PECK_RUNS['first_lunar_day'] + ('box_0e_elec_takahashi_first',),
                          second_lunar_day=SAUNDERS_PECK_RUNS['second_lunar_day'] + ('box_0e_elec_takahashi',)),
    'leakage': dict(first_lunar_day=SAUNDERS_PECK_RUNS['first_lunar_day'] + ('box_0e_elec_leakage_first',),
                    second_lunar_day=SAUNDERS_PECK_RUNS['second_lunar_day'] + ('box_0e_elec_leakage',)),
    'leader_crossing': dict(second_lunar_day=('box_0e_elec_leakage', 'box_0e_elec_leader_1kv', 'box_0e_elec_leader_10kv')),
}
COMMANDS = {'laws': 'charging_laws', 'leakage': 'leakage', 'leader': 'leader_crossing'}
LAWS = {1: 'takahashi', 11: 'saunders_peck_wrf', 12: 'saunders_peck'}
# the non-inductive totals over graupel and hail together (ctghi is ctgi and cthi, ctghs ctgs and cths) and snow with
# cloud ice; the inductive ones, droplets rebounding from graupel and from snow in the field
NONINDUCTIVE, INDUCTIVE = ('ctghi', 'ctghs', 'ctswi'), ('ctghw', 'ctsww')


def window_segments(case: Path, t0: float, t1: float) -> list:
    """The runner's records of the segments whose model time overlaps (t0, t1]."""
    progress = json.loads((case / 'progress.json').read_text())
    start, out = progress.get('start_s', 0.0), []
    for seg in progress['segments']:
        if seg['model_s'] > t0 and start < t1:
            out.append(seg)
        start = seg['model_s']
    return out


def window(name: str, t0: float, t1: float) -> dict:
    """One run's storms over model seconds (t0, t1]: its flashes (all, by hour and by kind), the charge each collision
    type separated, the largest field and charge, and at each output time the graupel aloft, the strongest updraft and
    where the charge sits."""
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    cfg = record['configuration']
    grid = record['grid']
    dx, dy = grid['dx_m'], grid.get('dy_m', grid['dx_m'])
    zh = heights_km(case / 'cm1out_s.ctl')
    zw = np.loadtxt(case / 'input_grid_z')
    inside = lambda t: (t > t0) & (t <= t1)

    f = read_log(case / 'terluna_flashes.txt', FLASH_COLUMNS)
    f = {k: v[inside(f['time_s'])] for k, v in f.items()}
    edges = np.array([t0, t1])
    lightning = dict(count=int(f['time_s'].size))
    if f['time_s'].size:
        hours = ((f['time_s'] - t0) // 3600.0).astype(int)
        per_hour = np.bincount(hours)
        lightning.update(first_h=float((f['time_s'].min() - t0) / 3600.0),
                         last_h=float((f['time_s'].max() - t0) / 3600.0),
                         hours_with_flashes=int(np.count_nonzero(per_hour)), most_in_an_hour=int(per_hour.max()),
                         nox_total_mol=float(f['nox_mol'].sum()))
        for code, label in KINDS.items():
            sel = f['kind'] == code
            if sel.any():
                lightning[label] = flash_summary({k: v[sel] for k, v in f.items()}, edges)
                del lightning[label]['per_bin']

    segments = window_segments(case, t0, t1)
    rates = charging(case, [case / seg['log'] for seg in segments if (case / seg['log']).exists()])
    sel = inside(rates['time_s'])
    times = rates['time_s'][sel]
    dt = np.diff(times, prepend=times[0]) if times.size else times
    separated = {}
    for key, label in CHARGING.items():
        r = rates[key][sel]
        if r.size and not np.all(np.isnan(r)):
            separated[key] = dict(collisions=label, positive_c=float(np.nansum(r[:, 1] * dt)),
                                  negative_c=float(np.nansum(r[:, 0] * dt)), peak_c_s=float(np.nanmax(np.abs(r))))
    total = lambda keys: float(sum(separated[k]['positive_c'] - separated[k]['negative_c'] for k in keys
                                   if k in separated))
    charge = dict(by_collisions=separated, noninductive_c=total(NONINDUCTIVE), inductive_c=total(INDUCTIVE),
                  logged_hours=float(dt.sum() / 3600.0))

    fl = read_log(case / 'terluna_field.txt', FIELD_COLUMNS)
    fl = {k: v[inside(fl['time_s'])] for k, v in fl.items()}
    field = dict(e_max_kv_m=float(fl['e_max_v_m'].max() / 1e3), positive_c_max=float(fl['positive_c'].max()),
                 negative_c_max=float(-fl['negative_c'].min()), energy_j_max=float(fl['energy_j'].max()),
                 noninductive_max_pc_m3_s=float(fl['noninductive_max_c_m3_s'].max() * 1e12),
                 hours_above_100_kv_m=float(np.sum(fl['e_max_v_m'] > 1.0e5) * np.median(np.diff(fl['time_s'])) / 3600.0)
                 if fl['time_s'].size > 1 else None) if fl['time_s'].size else {}

    g = read_log(case / 'terluna_ground.txt', GROUND_COLUMNS)
    ground = ground_summary({k: v[inside(g['time_s'])] for k, v in g.items()})

    output_s = cfg['output_s'] * record.get('time_scale', 1.0)
    snapshots = []
    for n in sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat')):
        t = (n - 1) * output_s
        if not t0 < t <= t1:
            continue
        s = snapshot_charge(case, n, dx, dy, zw)
        st = structure(s, zh)
        ions = s['charges']['small_ions'] * s['volume']
        st.update(hour=(t - t0) / 3600.0, graupel_kg=float(np.sum(s['graupel'] * s['volume'])),
                  w_max_m_s=float(np.nanmax(s['w'])), net_c=float(np.sum(s['net'] * s['volume'])),
                  ions_positive_c=float(ions[ions > 0].sum()), ions_negative_c=float(ions[ions < 0].sum()))
        snapshots.append(st)
    return dict(case=name, law=LAWS.get(int(record['electricity']['isaund']), record['electricity']['isaund']),
                electricity={k: record['electricity'].get(k) for k in
                             ('isaund', 'lightning', 'ground_m', 'corona_v_m', 'leakage', 'substep_s', 'leader_v_m')},
                executables=sorted({seg.get('executable') for seg in segments}),
                lightning=lightning, charge=charge, field=field, ground=ground, snapshots=snapshots)


def comparison(name: str) -> dict:
    """Each window's run with the change (COMPARISONS) and its runs without it, each summarized over the window's
    days."""
    out = {}
    for label, runs in COMPARISONS[name].items():
        t0, t1 = (86400.0 * d for d in WINDOW_DAYS[label])
        out[label] = dict(days=list(WINDOW_DAYS[label]),
                          runs=[window(run, t0, t1) for run in runs if (ra.RUNS / run / 'progress.json').exists()])
    return dict(schema=SCHEMA, kind=name.replace('_', '-'), windows=out,
                producer=dict(domain='climate', files={'crm/elec_analysis.py': hashlib.sha256(
                    Path(__file__).read_bytes()).hexdigest()[:16]}),
                evidence=EVIDENCE, reading_rule=READING_RULE + (
                    ' Each window\'s runs start from the same restart of box_0e_elec (or are box_0e_elec itself) and '
                    'differ in the settings listed with each; the runs with the change compared come last. Hours count '
                    'from the window\'s start. Separated charge is the domain\'s '
                    'charging of both signs integrated over the window (C); graupel is graupel and hail together (kg) '
                    'at each three-hourly output, where the net charge and the charge on the small ions are also '
                    'given (C).'))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case', help='a run, or "laws", "leakage" or "leader" for the windows run again under '
                                     'Takahashi\'s law, with leakage or with the leader\'s crossing '
                                     '(elec_charging_laws.json, elec_leakage.json, elec_leader_crossing.json)')
    parser.add_argument('--against', help='another run of the same case to set the storm beside')
    args = parser.parse_args(argv)
    if args.case in COMMANDS:
        RESULTS.mkdir(parents=True, exist_ok=True)
        name = COMMANDS[args.case]
        path = RESULTS / f'elec_{name}.json'
        path.write_text(json.dumps(significant(comparison(name)), indent=1, default=float) + '\n')
        print(path)
        return 0
    result = analyse(args.case)
    if args.against:
        result['against'] = against(args.case, args.against)
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f'elec_{args.case}.json'
    path.write_text(json.dumps(significant(result), indent=1, default=float) + '\n')
    print(path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
