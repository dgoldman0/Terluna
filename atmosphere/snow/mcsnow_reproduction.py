"""Reconstruct Welss's archived single-crystal cases with identified native code.

The original figure executable/driver was not found in the inspected public
artifacts. The native 2024 source and
its explicitly documented historical ventilation change are separate candidates.
No archived source/build or preceding runner is modified. External bytes stay in
the existing snow cache, which may be linked to the project drive.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time
import zipfile

import numpy as np
from atmosphere.snow import mcsnow_crystal_growth as native

ROOT = native.ROOT
CACHE = ROOT / 'research/runs/snow_scale/mcsnow_reproduction'
BASELINE = ROOT / 'research/runs/snow_scale/mcsnow_growth/build_2024'
DEPOSITION_SHA = 'a0904ce5e1061b559d9fe3d9c5c14329803dc823c9179fe43a3048f87f402232'
ARCHIVE_MD5 = 'dda30125b9ef8607b290d03649ee1d55'
ARCHIVE_URL = 'https://zenodo.org/api/records/10077724/files/McSnow_paper_plot_data.zip/content'
TARGETS = {
    'HP': '1_IGF_1_depend_0_depo_1.out',
    'IGF1': '1_IGF_1_depend_1_depo_1.out',
    'IGF2plus': '1_IGF_2_depend_1_depo_1.out',
    'IGF2_limit': '1_IGF_2_depend_1_depo_1_limit.out',
}
QUANTITIES = ('mass_kg', 'equatorial_dimension_m', 'polar_dimension_m',
              'ice_volume_m3', 'aspect_ratio', 'terminal_speed_m_s', 'diameter_m')


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def replace_once(body, before, after):
    if body.count(before) != 1:
        raise ValueError('Patch context differs: ' + before[:100])
    return body.replace(before, after, 1)


def historical_deposition(body):
    """Apply only the two deletions identified by the pinned source comment."""
    if hashlib.sha256(body.encode()).hexdigest() != DEPOSITION_SHA:
        raise ValueError('Unrecognized deposition source')
    first = '    Dv = Dv*ventilation_coef_i( sqrtRe*cubrSc, sp%phi )  ! CHRS for old bugged version (reproduce paper) comment these two lines'
    second = '    Ka = Ka*ventilation_coef_i( sqrtRe*cubrPr, sp%phi )'
    pair = first + '\n' + second
    return replace_once(body, pair, '! Historical reconstruction: two operations disabled as directed by the source.\n!' + first + '\n!' + second)


def reconstruction_driver(body):
    """Add a native saturation switch; retain every growth/shape integration call."""
    body = replace_once(body, '  logical :: delete_particle, prolate_interpolation',
                        '  logical :: delete_particle, prolate_interpolation, native_saturation')
    body = replace_once(body, 'initial_polar_m, growth_mode, prolate_interpolation\n',
                        'initial_polar_m, growth_mode, prolate_interpolation, native_saturation\n')
    body = replace_once(body, '  call get_command_argument(1,input_path)',
                        '  native_saturation=.false.\n  call get_command_argument(1,input_path)')
    body = replace_once(body, '  air%psatw=psat_water(air%T); air%psati=psat_ice(air%T)',
                        '  air%psatw=psat_water(air%T); air%psati=psat_ice(air%T)\n'
                        '  if (native_saturation) air%pv=air%psatw')
    return body


def deposition_candidate(body, historical=False, leading_heat=False):
    """Keep each source-supported diagnostic distinct from an adopted repair."""
    if hashlib.sha256(body.encode()).hexdigest() != DEPOSITION_SHA:
        raise ValueError('Unrecognized deposition source')
    if historical:
        body = historical_deposition(body)
    if leading_heat:
        # Equation 4's leading thermal term; all other source choices stay fixed.
        # With ventilation enabled this still retains native unequal fv/fh.
        body = replace_once(body,
            '  Hi = (L/(Rv*atmo%T) - 1.0_wp)*(L*Dv*es)/(Ka*Rv*atmo%T**2)',
            '  Hi = L**2*Dv*es/(Ka*Rv**2*atmo%T**3)')
    return body


def build(output, historical=False, leading_heat=False):
    """Reuse the native dependency selector and flags; preserve a complete diff."""
    output = Path(output).resolve()
    baseline = native.verify_build(BASELINE)
    if baseline['source_identity']['version'] != '2024':
        raise ValueError('Historical comparison requires the pinned 2024 source')
    if output.exists():
        result = verify_build(output)
        if (result['historical_mean_ventilation'] != historical or
                result.get('leading_heat_term', False) != leading_heat):
            raise ValueError('Existing build is a different candidate')
        return result
    output.mkdir(parents=True)
    shutil.copytree(BASELINE / 'source', output / 'source')
    original = (output / 'source/src/mo_depo.f90').read_text()
    if native.sha(output / 'source/src/mo_depo.f90') != DEPOSITION_SHA:
        raise ValueError('Unrecognized baseline')
    changed = deposition_candidate(original, historical, leading_heat)
    (output / 'source/src/mo_depo.f90').write_text(changed)
    old_driver = native.DRIVER.read_text()
    caller = output / native.DRIVER.name
    caller.write_text(reconstruction_driver(old_driver))
    diff = ''.join(difflib.unified_diff(original.splitlines(True), changed.splitlines(True),
                                      fromfile='native/src/mo_depo.f90', tofile='candidate/src/mo_depo.f90'))
    diff += ''.join(difflib.unified_diff(old_driver.splitlines(True), caller.read_text().splitlines(True),
                                       fromfile='existing/' + native.DRIVER.name,
                                       tofile='candidate/' + native.DRIVER.name))
    (output / 'reconstruction.patch').write_text(diff)
    for name in native.TABLE_NAMES:
        shutil.copy2(BASELINE / name, output / name)
    compiler = baseline['compiler_path']
    order = native.selected_sources(output / 'source/src', caller, compiler)
    commands, objects = [], []
    started = time.monotonic()
    with (output / 'build.log').open('w') as log:
        for i, p in enumerate(order):
            obj = output / f'{i:02d}_{p.stem}.o'
            objects.append(str(obj))
            cmd = [compiler, *native.FLAGS, '-I' + str(output / 'source/src'),
                   '-I' + str(output), '-J' + str(output), '-c', str(p), '-o', str(obj)]
            commands.append(cmd)
            subprocess.run(cmd, cwd=output, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
        cmd = [compiler, *native.FLAGS, *objects, '-o', str(output / 'crystal')]
        commands.append(cmd)
        subprocess.run(cmd, cwd=output, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
    record = dict(schema='terluna.atmosphere.mcsnow-reconstruction-build/1',
                  producer_sha256=native.sha(__file__), baseline_record_sha256=native.sha(BASELINE / 'build.json'),
                  source_identity=baseline['source_identity'], historical_mean_ventilation=historical,
                  leading_heat_term=leading_heat,
                  source_files={str(p.relative_to(output / 'source')): native.sha(p)
                                for p in sorted((output / 'source').rglob('*')) if p.is_file()},
                  driver_sha256=native.sha(caller), original_driver_sha256=native.sha(native.DRIVER),
                  table_identity=baseline['table_identity'], patch_sha256=native.sha(output / 'reconstruction.patch'),
                  compiler_path=compiler, compiler_sha256=native.sha(compiler),
                  compiler_version=subprocess.check_output([compiler, '--version'], text=True).splitlines()[0],
                  flags=native.FLAGS, commands=commands, executable_sha256=native.sha(output / 'crystal'),
                  build_log_sha256=native.sha(output / 'build.log'), elapsed_s=time.monotonic()-started,
                  evidence='Reconstruction candidate; original figure driver/executable not recovered from inspected artifacts.')
    save(output / 'build.json', record)
    return record


def verify_build(path):
    path = Path(path)
    record = native.verify_build(path)
    for name, key in [('reconstruction.patch', 'patch_sha256'), ('build.log', 'build_log_sha256')]:
        if native.sha(path / name) != record[key]:
            raise ValueError('Changed build artifact: ' + name)
    expected = deposition_candidate((BASELINE / 'source/src/mo_depo.f90').read_text(),
                                   record['historical_mean_ventilation'],
                                   record.get('leading_heat_term', False))
    if (path / 'source/src/mo_depo.f90').read_text() != expected:
        raise ValueError('Build does not implement its declared candidate')
    if (path / native.DRIVER.name).read_text() != reconstruction_driver(native.DRIVER.read_text()):
        raise ValueError('Build uses a different reconstruction caller')
    return record


def author_data(cache=CACHE, variant='IGF2plus', duration=600):
    """Retain decimal rounding bounds and explicit semiaxis-to-full-axis mapping."""
    archive = Path(cache) / 'McSnow_paper_plot_data.zip'
    if hashlib.md5(archive.read_bytes()).hexdigest() != ARCHIVE_MD5:
        raise ValueError('Author archive checksum mismatch')
    member = f'T_range_-30-0_p1010.0_dt1.0_t_end{duration}_r_start2e-06_mD_' + TARGETS[variant]
    with zipfile.ZipFile(archive) as z:
        raw = z.read(member)
    rows = parse_author_rows(raw, duration)
    return dict(member=member, member_sha256=hashlib.sha256(raw).hexdigest(), rows=rows,
                archive_md5=ARCHIVE_MD5, archive_sha256=native.sha(archive),
                axes='Source a,c are semiaxes: verified by D=2*max(a,c) and V=4*pi*a*a*c/3 within rounding.',
                attribution='Welss, Siewert and Seifert; doi:10.5281/zenodo.10077724; CC-BY-4.0')


def parse_author_rows(raw, duration):
    """Parse nine columns explicitly; never infer a hidden dataframe index."""
    rows = []
    for line in raw.decode().splitlines():
        fields = line.split(',')
        if len(fields) != 9:
            raise ValueError('Expected elapsed time plus eight author columns')
        time_s, T, a, c, phi, m, V, vt, D = map(float, fields)
        if (time_s != duration or not all(math.isfinite(v) for v in (T, a, c, phi, m, V, vt, D))
                or min(T, a, c, phi, m, V, D) <= 0 or vt < 0):
            raise ValueError('Malformed author state')
        values = dict(mass_kg=m, equatorial_dimension_m=2*a, polar_dimension_m=2*c,
                      ice_volume_m3=V, aspect_ratio=phi, terminal_speed_m_s=vt, diameter_m=D)
        # Author fields are fixed decimal, with distinct precision for each column.
        def half_quantum(field):
            from decimal import Decimal
            return float(Decimal(10) ** Decimal(field).as_tuple().exponent / 2)
        errors = dict(zip(QUANTITIES, [half_quantum(fields[5]), 2*half_quantum(fields[2]),
                                      2*half_quantum(fields[3]), half_quantum(fields[6]),
                                      half_quantum(fields[4]), half_quantum(fields[7]), half_quantum(fields[8])]))
        da, dc = errors['equatorial_dimension_m']/2, errors['polar_dimension_m']/2
        for quantity, lower, upper in (
                ('diameter_m', 2*max(a-da, c-dc), 2*max(a+da, c+dc)),
                ('ice_volume_m3', 4*math.pi/3*(a-da)**2*(c-dc), 4*math.pi/3*(a+da)**2*(c+dc)),
                ('aspect_ratio', (c-dc)/(a+da), (c+dc)/(a-da))):
            if values[quantity]+errors[quantity] < lower or values[quantity]-errors[quantity] > upper:
                raise ValueError('Author geometry inconsistent within printed rounding: ' + quantity)
        rows.append(dict(temperature_K=T, elapsed_s=time_s, values=values, rounding_half_width=errors,
                         raw_density_kg_m3=m/V, plotted_density_kg_m3=min(m/V, 919.)))
    if not rows or len({r['temperature_K'] for r in rows}) != len(rows):
        raise ValueError('Empty or duplicate author temperature rows')
    return rows


def configuration(T, dt=1., duration=600., growth=2, native_saturation=True, prolate=True, pressure=101000.):
    for value in (T, dt, duration, pressure):
        if not math.isfinite(value) or value <= 0:
            raise ValueError('Positive finite environmental/time settings required')
    if growth not in (1, 2) or T > 273.15:
        raise ValueError('Unsupported growth package or warm temperature')
    for interval in (duration, 60.):
        if not math.isclose(interval / dt, round(interval / dt), rel_tol=0, abs_tol=1e-8):
            raise ValueError('Timestep must divide duration and output interval')
    return dict(temperature_K=float(T), pressure_Pa=float(pressure),
                vapor_pressure_Pa=float(native.scale.lunar.flakes.liquid_saturation_pressure(T)),
                native_saturation=bool(native_saturation), dt_s=float(dt), duration_s=float(duration),
                output_interval_s=60., initial_mass_kg=-1., initial_equatorial_m=4e-6,
                initial_polar_m=4e-6, growth_mode=growth, prolate_interpolation=bool(prolate))


def verify_run(path, build_dir=None):
    path = Path(path)
    record = json.loads((path / 'run.json').read_text())
    if build_dir is not None:
        verify_build(build_dir)
        if native.sha(Path(build_dir) / 'build.json') != record['build_record_sha256']:
            raise ValueError('Run/build identity mismatch')
    for name, digest in record['file_sha256'].items():
        if native.sha(path / name) != digest:
            raise ValueError('Changed raw run artifact: ' + name)
    if (path / 'crystal.nml').read_text() != native.namelist(record['configuration']):
        raise ValueError('Run configuration differs from saved input')
    rows, air, errors = analyze_output(path, record['configuration'], record['returncode'])
    if (rows != record['snapshots'] or air != record['air'] or errors != record['validation_errors'] or
            record['complete'] != (record['returncode'] == 0 and not errors)):
        raise ValueError('Saved analysis differs from raw output replay')
    return record


def run(build_dir, output, config):
    build_dir, output = Path(build_dir).resolve(), Path(output).resolve()
    if output.exists():
        existing = verify_run(output, build_dir)
        if existing['configuration'] != config:
            raise ValueError('Existing run has different configuration')
        return existing
    identity = verify_build(build_dir)
    output.mkdir(parents=True)
    (output / 'crystal.nml').write_text(native.namelist(config))
    for name in native.TABLE_NAMES:
        shutil.copy2(build_dir / name, output / name)
    started = time.monotonic()
    with (output / 'run.log').open('w') as log:
        completed = subprocess.run([str(build_dir / 'crystal'), 'crystal.nml'], cwd=output,
                                   env=native.ENV, stdout=log, stderr=subprocess.STDOUT)
    rows, air, errors = analyze_output(output, config, completed.returncode)
    record = dict(schema='terluna.atmosphere.mcsnow-reconstruction-run/1',
                  producer_sha256=native.sha(__file__), build_record_sha256=native.sha(build_dir / 'build.json'),
                  historical_mean_ventilation=identity['historical_mean_ventilation'], configuration=config,
                  leading_heat_term=identity.get('leading_heat_term', False),
                  returncode=completed.returncode, complete=completed.returncode == 0 and not errors,
                  validation_errors=errors, snapshots=rows, air=air, elapsed_s=time.monotonic()-started,
                  file_sha256={p.name: native.sha(p) for p in sorted(output.iterdir()) if p.is_file()},
                  evidence='Prescribed-air source reconstruction; no finite-water storm or lunar validation.')
    save(output / 'run.json', record)
    return record


def analyze_output(output, config, returncode):
    """Retain partial numeric prefixes and solver failure as separate evidence."""
    output = Path(output)
    rows, errors = native.read_csv(output / 'trajectory.csv', ('elapsed_s', *QUANTITIES, 'density_kg_m3',
                                                             'ice_mass_balance_residual_kg'))
    air, air_errors = native.read_csv(output / 'air.csv', ('temperature_K', 'pressure_Pa', 'rh_percent',
                                                        'vapor_pressure_Pa', 'native_ice_density_kg_m3'))
    errors += air_errors
    expected = [float(t) for t in np.arange(0, config['duration_s']+1e-8, config['output_interval_s'])]
    if expected[-1] != config['duration_s']:
        expected.append(config['duration_s'])
    if len(air) != 1:
        errors.append('Expected one prescribed-air state')
    else:
        for key in ('temperature_K', 'pressure_Pa'):
            if not math.isclose(air[0][key], config[key], rel_tol=1e-12):
                errors.append('Prescribed-air input mismatch: ' + key)
        if config['native_saturation']:
            if abs(air[0]['rh_percent']-100.) > 1e-12:
                errors.append('Native water-saturation mismatch')
        elif not math.isclose(air[0]['vapor_pressure_Pa'], config['vapor_pressure_Pa'], rel_tol=1e-12):
            errors.append('Prescribed vapour input mismatch')
    if [r['elapsed_s'] for r in rows] != expected[:len(rows)]:
        errors.append('Irregular saved trajectory times')
    if rows and air:
        mass = config['initial_mass_kg']
        if mass < 0:
            mass = math.pi/6*config['initial_equatorial_m']**2*config['initial_polar_m']*air[0]['native_ice_density_kg_m3']
        initial = dict(mass_kg=mass, equatorial_dimension_m=config['initial_equatorial_m'],
                       polar_dimension_m=config['initial_polar_m'])
        if any(not math.isclose(rows[0][k], v, rel_tol=1e-12) for k,v in initial.items()):
            errors.append('Initial mass/axis mapping mismatch')
    for r in rows:
        if min(r[k] for k in QUANTITIES if k != 'terminal_speed_m_s') <= 0 or r['terminal_speed_m_s'] < 0:
            errors.append('Nonpositive crystal state')
            break
        if abs(r['ice_mass_balance_residual_kg']) / r['mass_kg'] > 1e-10:
            errors.append('Ice deposition account differs')
            break
    if returncode == 0 and len(rows) != len(expected):
        errors.append('Successful exit with incomplete trajectory')
    return rows, air, errors


def compare(snapshot, author):
    residual = {k: snapshot[k]-author['values'][k] for k in QUANTITIES}
    return dict(relative_error={k: residual[k]/author['values'][k] for k in QUANTITIES},
                rounding_units={k: abs(residual[k])/author['rounding_half_width'][k] for k in QUANTITIES},
                all_within_printed_rounding=all(abs(residual[k]) <= author['rounding_half_width'][k] for k in QUANTITIES))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['build', 'probe'])
    parser.add_argument('--cache', type=Path, default=CACHE)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    for historical, name in [(False, 'native'), (True, 'historical')]:
        build(args.cache / ('build_' + name), historical)
    if args.action == 'build':
        return
    product = []
    for variant, growth in [('IGF1', 1), ('IGF2plus', 2)]:
        target = author_data(args.cache, variant)
        for T in (258.65, 262.65, 267.9):
            author = next(r for r in target['rows'] if r['temperature_K'] == T)
            for mode in ('native', 'historical'):
                for saturation in (True, False):
                    config = configuration(T, growth=growth, native_saturation=saturation)
                    name = f'{variant}_{mode}_T{T}_sat{int(saturation)}'
                    result = run(args.cache / ('build_' + mode), args.cache / 'runs' / name, config)
                    row = dict(case=name, complete=result['complete'], target=target['member'],
                               **compare(result['snapshots'][-1], author))
                    product.append(row)
                    print(name, result['complete'], {k: round(v, 6) for k, v in row['relative_error'].items()}, flush=True)
    save(args.cache / 'probe.json', product)


if __name__ == '__main__':
    main()
