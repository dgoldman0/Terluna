"""Export resolved dusk microphysics from a corrected CM1 equatorial ring.

The output retains individual columns and species for illumination studies.
Thermodynamic coefficients are read from the run's hash-checked CM1 source.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

import numpy as np

from climate.crm.cm1_run import CM1_HOME
from climate.crm.ring_analysis import case_geometry, grads_layout, hour_angle, read_snapshot
from shared.constants import SYNODIC_MONTH_DAYS

SCHEMA = 'terluna.climate.cm1-dusk-columns/1'
SPECIES = ('qc', 'qr', 'qi', 'qs', 'qg')
NUMBERS = ('nci', 'ncr', 'ncs', 'ncg')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_cm1_constants(record, source_root=CM1_HOME):
    """Read the actual default thermodynamics and Morrison liquid density."""
    directory = Path(source_root) / 'build' / record['build']['label'] / 'src'
    sources = {}
    for name in ('constants.F', 'solve2.F', 'morrison.F'):
        path = directory / name
        value = sha256(path)
        if not value.startswith(record['build']['patched'][name]):
            raise ValueError(f'CM1 source differs from the recorded build: {path}')
        sources[name] = value
    constants = (directory / 'constants.F').read_text()
    default = constants.split('!  cm1 defaults:')[1].split('ENDIF')[0]

    def number(text, name):
        match = re.search(r'(?im)^\s*(?:real,\s*parameter\s*::\s*)?' + name + r'\s*=\s*([0-9.eEdD+\-]+)', text)
        if match is None:
            raise ValueError(f'Missing CM1 coefficient {name}')
        return float(match[1].lower().replace('d', 'e'))

    values = {key: number(default if key in ('rd', 'cp') else constants, key) for key in ('rd', 'cp', 'rv', 'p00')}
    values['liquid_density_kg_m3'] = number((directory / 'morrison.F').read_text(), 'RHOW')
    return values, sources


def dry_density(pressure, potential_temperature, vapor_mixing_ratio, constants):
    p, theta, qv = map(lambda x: np.asarray(x, float), (pressure, potential_temperature, vapor_mixing_ratio))
    if (p.shape != theta.shape or p.shape != qv.shape or
            not all(np.isfinite(x).all() for x in (p, theta, qv)) or
            np.any(p <= 0) or np.any(theta <= 0) or np.any(qv < 0)):
        raise ValueError('Invalid CM1 thermodynamic fields')
    temperature = theta * (p / constants['p00']) ** (constants['rd'] / constants['cp'])
    return p / (temperature * (constants['rd'] + constants['rv'] * qv)), temperature


def layer_overlap(height_edges_m, lower_m, upper_m):
    z = np.asarray(height_edges_m, float)
    if (z.ndim != 1 or len(z) < 2 or not np.isfinite(z).all() or np.any(np.diff(z) <= 0) or
            not np.isfinite([lower_m, upper_m]).all() or lower_m >= upper_m):
        raise ValueError('Invalid layer boundaries')
    return np.maximum(0, np.minimum(z[1:], upper_m) - np.maximum(z[:-1], lower_m))


def water_path(density, mixing_ratio, height_edges_m, lower_m=0., upper_m=None):
    """kg/m2 from dry-air mixing ratio, with exact partial-layer band overlap."""
    z = np.asarray(height_edges_m)
    rho, q = np.asarray(density), np.asarray(mixing_ratio)
    if (rho.shape != q.shape or rho.shape[-1] != len(z) - 1 or
            not np.isfinite(rho).all() or not np.isfinite(q).all() or np.any(rho <= 0) or np.any(q < 0)):
        raise ValueError('Invalid condensate columns')
    return (rho * q) @ layer_overlap(z, lower_m, z[-1] if upper_m is None else upper_m)


def diagnostic_field(case, snapshot, name):
    nx, layout = grads_layout(case / 'cm1out_diag_s.ctl')
    path = case / f'cm1out_diag_{snapshot:06d}_s.dat'
    data = np.fromfile(path, dtype='<f4')
    if nx != 1 or data.size != sum(n for _, n in layout):
        raise ValueError('Expected CM1 domain-mean diagnostic columns')
    offset = 0
    for key, levels in layout:
        if key == name:
            return data[offset:offset + levels].astype(float), path
        offset += levels
    raise ValueError(f'Missing diagnostic {name}')


def export(case, output):
    case, output = Path(case), Path(output)
    geo = case_geometry(case)
    record = geo['record']
    if (record['case'] != 'ring_equator' or record['reference']['run'] != 'A28_dim5_moon' or
            record['path']['tilt_deg'] != 0 or record['time_scale'] != 1 or grads_layout(case / 'cm1out_s.ctl')[0] != record['grid']['nx']):
        raise ValueError('Corrected flat equatorial ring required')
    namelist = (case / 'namelist.input').read_text()
    if not re.search(r'\btestcase\s*=\s*0\s*,', namelist):
        raise ValueError('CM1 default thermodynamics required')
    constants, source_hashes = read_cm1_constants(record)
    fields = ('prs', 'th', 'qv', *SPECIES, *NUMBERS)
    retained = {k: [] for k in (*fields, 'density_kg_m3', 'temperature_k')}
    times, columns, angles, snapshots = [], [], [], []
    raw_inputs, density_errors = [], []
    dt = record['configuration']['output_s']
    expected = range(int(np.ceil(SYNODIC_MONTH_DAYS * 86400 / dt)) + 1,
                     int(np.floor(min(record['configuration']['days'], 2 * SYNODIC_MONTH_DAYS) * 86400 / dt)) + 2)
    for n in expected:
        path = case / f'cm1out_t{n:06d}_s.dat'
        day = (n - 1) * dt / 86400
        if day > 2 * SYNODIC_MONTH_DAYS:
            continue
        d = read_snapshot(case, n)
        if any(not np.isfinite(d[k]).all() or np.any(d[k] < 0) for k in fields):
            raise ValueError(f'Invalid thermodynamics or microphysics in {path.name}')
        rho, temperature = dry_density(d['prs'], d['th'], d['qv'], constants)
        diagnostic, diagnostic_path = diagnostic_field(case, n, 'rho')
        discrepancy = float(np.max(abs(rho.mean(axis=1) / diagnostic - 1)))
        if discrepancy > 2e-5:
            raise ValueError(f'Density reconstruction disagrees with CM1 diagnostics: {discrepancy}')
        density_errors.append(discrepancy)
        h = hour_angle(geo, day * 86400, SYNODIC_MONTH_DAYS * 86400)
        selected = np.flatnonzero((h >= 90) & (h < 110))
        selected = selected[np.argsort(h[selected])]
        times.extend([day] * len(selected)); columns.extend(selected.tolist())
        angles.extend(h[selected].tolist()); snapshots.extend([n] * len(selected))
        for key in fields:
            retained[key].append(d[key][:, selected].T.copy())
        retained['density_kg_m3'].append(rho[:, selected].T.astype('f4'))
        retained['temperature_k'].append(temperature[:, selected].T.astype('f4'))
        raw_inputs.append(dict(snapshot=n, time_day=day, scalar_file=path.name, scalar_sha256=sha256(path),
                               diagnostic_file=diagnostic_path.name, diagnostic_sha256=sha256(diagnostic_path)))
        if len(raw_inputs) % 40 == 0:
            print(f'Recovered {len(raw_inputs)} dusk snapshots', flush=True)
    if not raw_inputs:
        raise ValueError('Second-cycle snapshots absent')
    metadata = dict(schema=SCHEMA,
                    producer=dict(lane='climate', file='climate/crm/cloud_columns.py', sha256=sha256(__file__),
                                  reader_sha256=sha256(Path(__file__).with_name('ring_analysis.py'))),
                    evidence='Resolved microphysics from the corrected CM1 equatorial ring, at 6-km horizontal spacing. '
                             'Only samples after the first complete synodic cycle enter this product. '
                             'The 2-D run represents a selected circulation experiment; particle habits and '
                             'three-dimensional cloud morphology remain model assumptions.',
                    reading_rule='Rows are individual columns, grouped by increasing snapshot time then local hour angle. '
                                 'Dusk is local hour angle [90,110) degrees. Arrays are (row,height); '
                                 'mixing ratios and particle counts use dry-air mass. Raw model levels are preserved. '
                                 'First-cycle data serve spin-up. Coverage ends at the last saved three-hour snapshot.',
                    units=dict(time_day='Earth days since simulation start', hour_angle_deg='degrees after local noon',
                               height_m='metres above sea level', species='kg per kg dry air', numbers='particles per kg dry air',
                               density_kg_m3='dry-air kg m^-3', temperature_k='K'),
                    case=record['case'], gcm_run=record['reference']['run'],
                    requested_span_days=[SYNODIC_MONTH_DAYS, 2 * SYNODIC_MONTH_DAYS],
                    sampled_span_days=[raw_inputs[0]['time_day'], raw_inputs[-1]['time_day']],
                    snapshot_spacing_hours=dt / 3600, snapshots=len(raw_inputs), columns=len(times),
                    grid=record['grid'], constants=constants, source_hashes=source_hashes,
                    source_archive=record['build']['source'], cloud_droplet_number_cm3=record['configuration']['droplets_cm3'],
                    fixed_inputs={name: sha256(case / name) for name in
                                  ('case.json', 'namelist.input', 'cm1out_s.ctl', 'cm1out_diag_s.ctl', 'input_grid_z', 'terluna_surface.txt')},
                    raw_inputs=raw_inputs, checks=dict(maximum_mean_density_relative_error=max(density_errors)))
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, metadata=json.dumps(metadata), time_day=times, column_index=columns,
                        hour_angle_deg=angles, snapshot=snapshots, height_edges_m=geo['zw'],
                        **{key: np.concatenate(value) for key, value in retained.items()})
    return metadata
