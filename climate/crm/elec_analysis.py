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
SCHEMA = 'terluna.climate.crm-electricity/1'
RD, CP, P00 = 287.04, 1005.7, 1.0e5
SPECIES = ('cloud_water', 'rain', 'cloud_ice', 'snow', 'graupel', 'small_ions', 'hail')   # passive tracers 1-7
FLASH_COLUMNS = ('time_s', 'step', 'substep', 'kind', 'x_m', 'y_m', 'z_m', 'e_v_m', 'e_break_v_m', 'e_over_break_max',
                 'extent_m2', 'regions', 'positive_c', 'negative_c', 'energy_before_j', 'energy_after_j', 'z_low_m',
                 'z_high_m', 'points', 'nox_mol')
FIELD_COLUMNS = ('time_s', 'step', 'e_max_v_m', 'x_m', 'y_m', 'z_m', 'positive_c', 'negative_c', 'energy_j',
                 'noninductive_max_c_m3_s', 'inductive_max_c_m3_s')
KINDS = {1: 'in_cloud', 2: 'negative_to_ground', 3: 'positive_to_ground', 9: 'cylinders'}
# WRF-ELEC's charging totals (module_mp_nssl_2mom.F, its driver): negative and positive parts of the domain's charging
# rate, C/s, by the collisions that separate it
CHARGING = {'ctghi': 'graupel and hail with cloud ice', 'ctghs': 'graupel and hail with snow',
            'ctswi': 'snow with cloud ice', 'ctghw': 'graupel with cloud droplets (inductive)',
            'ctgi': 'graupel with cloud ice', 'ctgs': 'graupel with snow', 'cthi': 'hail with cloud ice',
            'cths': 'hail with snow', 'ctsww': 'snow with cloud droplets'}
TIME_UNITS = {'sec': 1.0, 'min': 60.0, 'hour': 3600.0, 'hrs': 3600.0, 'day': 86400.0}
BIN_S = 300.0
EVIDENCE = ('CM1 r22.0 with WRF-ELEC\'s NSSL two-moment microphysics (MicroTed/wrf4-elec e43041b) and its charging '
            '(non-inductive by collisions of graupel and hail with ice and snow, inductive by droplets rebounding from '
            'graupel in the field), the small ions\' net charge attaching to particles, an FFT field solver, and '
            'WRF-ELEC\'s lightning on its own sub-steps of the sedimentation: the branched scheme (lightmsz, '
            'MacGorman, Straka and Ziegler 2001) or the cylindrical one. Charge is carried by bulk particle '
            'populations; branched channels run on the model grid and neutralize the charge where they reach; the '
            'ground-strike rule, the 50-MV potential and 10-kV/m field a ground flash needs at its start and the '
            'nitrogen oxide yield are WRF-ELEC\'s calibrations for Earth; nothing conducts charge away except where a '
            'case turns on leakage.')
READING_RULE = ('Times are model seconds from the start. Kinds: 1 in cloud, 2 negative and 3 positive charge to ground '
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


def charging(case: Path) -> dict:
    """WRF-ELEC's charging totals by step from CM1's logs, each assigned to the step whose time line follows it."""
    pattern = re.compile(r'^(ct[a-z]+)n,\1p\s*=\s*(\S+),\s*(\S+)')
    step = re.compile(r'^\s+(\d+)\s+([0-9.Ee+-]+)\s+(sec|min|hour|hrs|day)')
    pending, rows = {}, {}
    for log in segment_logs(case):
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
    return dict(charges=charges, net=net, t=t, rho=rho, volume=volume, w=snap.get('winterp'), dbz=snap.get('dbz'))


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


def flash_summary(f: dict, edges: np.ndarray) -> dict:
    """Flashes of one kind: how many and when, the charge each neutralized, the energy it released, where it started
    and how far its channels reached, and the nitrogen oxides it made."""
    n = f['time_s'].size
    dissipated = f['energy_before_j'] - f['energy_after_j']
    out = dict(count=int(n), first_s=float(f['time_s'].min()), per_5_min=binned(f['time_s'], np.ones(n), edges),
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
    rates = charging(case)
    end = float(max(field['time_s'].max() if field['time_s'].size else 0.0,
                    rates['time_s'].max() if rates['time_s'].size else 0.0))
    edges = np.arange(0.0, end + BIN_S, BIN_S)

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
            per_5_min=binned(flashes['time_s'], np.ones_like(flashes['time_s']), edges),
            # the field over breakdown before a sub-step's first flash (what charging and transport built since the
            # last) and before its later ones (what the flashes before them left)
            over_breakdown_first=stats(over[first]), over_breakdown_later=stats(over[~first]),
            per_substep_max=int(np.bincount(groups).max()))
    fields = dict(
        e_max_kv_m=float(field['e_max_v_m'].max() / 1000.0) if field['time_s'].size else None,
        e_max_per_5_min_kv_m=binned(field['time_s'], field['e_max_v_m'] / 1000.0, edges, 'max'),
        positive_c_per_5_min=binned(field['time_s'], field['positive_c'], edges, 'max'),
        negative_c_per_5_min=binned(field['time_s'], -field['negative_c'], edges, 'max'),
        energy_j_per_5_min=binned(field['time_s'], field['energy_j'], edges, 'max'),
        noninductive_max_pc_m3_s=float(field['noninductive_max_c_m3_s'].max() * 1e12) if field['time_s'].size else None,
        inductive_max_pc_m3_s=float(field['inductive_max_c_m3_s'].max() * 1e12) if field['time_s'].size else None,
        noninductive_max_per_5_min_pc_m3_s=binned(field['time_s'], field['noninductive_max_c_m3_s'] * 1e12, edges, 'max'),
        inductive_max_per_5_min_pc_m3_s=binned(field['time_s'], field['inductive_max_c_m3_s'] * 1e12, edges, 'max'))
    dt = np.diff(np.concatenate([[0.0], rates['time_s']])) if rates['time_s'].size else np.array([])
    separated = {}
    for key, label in CHARGING.items():
        r = rates[key]
        if not r.size or np.all(np.isnan(r)):
            continue
        separated[key] = dict(collisions=label, positive_c=float(np.nansum(r[:, 1] * dt)),
                              negative_c=float(np.nansum(r[:, 0] * dt)),
                              peak_c_s=float(np.nanmax(np.abs(r))),
                              per_5_min_c=binned(rates['time_s'], (r[:, 1] - r[:, 0]) * dt, edges))
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
        evidence=EVIDENCE, reading_rule=READING_RULE, end_s=end, bin_s=BIN_S,
        lightning=lightning, field=fields, charging=separated, structure=structures)


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


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--against', help='another run of the same case to set the storm beside')
    args = parser.parse_args(argv)
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
