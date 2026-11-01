"""Trace archived oblate states through unmodified native and author Python drag.

The Python helper comes verbatim from the identified authors' plotting archive.
There are no growth steps, fitted air parameters or changes to the native source.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import zipfile

from atmosphere.snow import mcsnow_reproduction as model

native = model.native
HERE = Path(__file__).resolve().parent
DRIVER = HERE / 'mcsnow_drag_driver.f90'
CACHE = model.CACHE / 'drag_trace'
PRODUCT = HERE / 'results/benchmark_drag_trace.json'
PLOT_SHA = '789e5a505f8680358fb427bdb53f5f814a139e3edc6ed55e15cc16cd7188dfb2'
UTILITY_SHA = 'eedcd82e73bf2d72108c3bf717c40d057c862365a21ae20ad58008598e8efbf1'
CONTROLS = ('archived', 'area_only', 'mass_only', 'solid_envelope')


def author_helper():
    """Execute only the inspected, hash-identified numerical helper, never plots."""
    archive = model.CACHE / 'McSnow_paper_plotscripts.zip'
    if native.sha(archive) != PLOT_SHA:
        raise ValueError('Plot archive checksum mismatch')
    with zipfile.ZipFile(archive) as z:
        raw = z.read('utility_functions.py')
        members = {name: hashlib.sha256(z.read(name)).hexdigest() for name in
                   ('utility_functions.py', 'fig10.py', 'fig13.py', 'fig16+18+20.py')}
    if hashlib.sha256(raw).hexdigest() != UTILITY_SHA:
        raise ValueError('Author helper checksum mismatch')
    namespace = {}
    exec(compile(raw, 'author:utility_functions.py', 'exec'), namespace)
    return namespace['vt_boehm'], dict(archive_sha256=PLOT_SHA, member_sha256=members,
        source='https://doi.org/10.5281/zenodo.10077724',
        attribution='Welss, Siewert and Seifert; CC-BY-4.0',
        scope='The oblate branch is used here. Prolate interpolation differs from the pinned native source.')


def selected_states():
    states, sources = [], {}
    for variant in ('HP', 'IGF1', 'IGF2plus'):
        for duration in ((600,) if variant == 'HP' else (180, 600)):
            target = model.author_data(variant=variant, duration=duration)
            sources[target['member']] = {k: v for k, v in target.items() if k != 'rows'}
            temperatures = (258.65,) if duration == 180 else (258.4, 258.65, 258.9, 262.65)
            for temperature in temperatures:
                row = next(r for r in target['rows'] if r['temperature_K'] == temperature)
                if row['values']['aspect_ratio'] > 1:
                    raise ValueError('Trace only supports oblate states')
                states.append(dict(case_id=len(states), variant=variant, source_member=target['member'], **row))
    return states, sources


def input_text(states):
    lines = []
    for row in states:
        v = row['values']
        values = (row['temperature_K'], 101000., v['mass_kg'], v['ice_volume_m3'], v['aspect_ratio'])
        lines.append(str(row['case_id'])+' '+ ' '.join(format(x, '.17g') for x in values))
    return '\n'.join(lines)+'\n'


def build(path):
    """Fresh compilation of unchanged baseline modules, without growth tables."""
    baseline = native.verify_build(model.BASELINE)
    path = Path(path)
    if path.exists():
        return verify_build(path)
    path.mkdir(parents=True)
    shutil.copyfile(DRIVER, path / DRIVER.name)
    source = model.BASELINE / 'source/src'
    compiler = baseline['compiler_path']
    selected = native.selected_sources(source, path / DRIVER.name, compiler)
    objects, commands = [], []
    with (path / 'build.log').open('w') as log:
        for i, src in enumerate(selected):
            obj = path / f'{i:02d}_{src.stem}.o'
            objects.append(str(obj))
            cmd = [compiler, *native.FLAGS, '-I'+str(source), '-I'+str(path), '-J'+str(path),
                   '-c', str(src), '-o', str(obj)]
            commands.append(cmd)
            subprocess.run(cmd, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
        cmd = [compiler, *native.FLAGS, *objects, '-o', str(path / 'trace')]
        commands.append(cmd)
        subprocess.run(cmd, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
    record = dict(schema='terluna.atmosphere.mcsnow-drag-build/1',
        baseline_record_sha256=native.sha(model.BASELINE/'build.json'), source_identity=baseline['source_identity'],
        compiled_sources={str(p): native.sha(p) for p in selected}, physics_patches=[],
        compiler_path=compiler, compiler_sha256=native.sha(compiler), compiler_version=baseline['compiler_version'],
        flags=native.FLAGS, commands=commands,
        files_sha256={p.name: native.sha(p) for p in [path/DRIVER.name, path/'trace', path/'build.log']})
    model.save(path/'build.json', record)
    return record


def verify_build(path):
    path = Path(path)
    record = json.loads((path/'build.json').read_text())
    native.verify_build(model.BASELINE)
    if native.sha(model.BASELINE/'build.json') != record['baseline_record_sha256']:
        raise ValueError('Baseline record changed')
    for src, digest in record['compiled_sources'].items():
        if native.sha(src) != digest:
            raise ValueError('Compiled source changed: '+src)
    for name, digest in record['files_sha256'].items():
        if native.sha(path/name) != digest:
            raise ValueError('Build artifact changed: '+name)
    if native.sha(DRIVER) != record['files_sha256'][DRIVER.name]:
        raise ValueError('Diagnostic caller changed')
    return record


def run(cache=CACHE):
    build(cache/'build')
    states, _ = selected_states()
    path = cache/'execution'
    if path.exists():
        return verify_run(cache, states)
    path.mkdir()
    (path/'input.txt').write_text(input_text(states))
    with (path/'trace.csv').open('w') as out, (path/'run.log').open('w') as log:
        result = subprocess.run([str(cache/'build/trace'), str(path/'input.txt')], cwd=path,
                                env=native.ENV, stdout=out, stderr=log)
    record = dict(returncode=result.returncode, build_record_sha256=native.sha(cache/'build/build.json'),
        files_sha256={name: native.sha(path/name) for name in ('input.txt', 'trace.csv', 'run.log')})
    model.save(path/'run.json', record)
    return verify_run(cache, states)


def verify_run(cache, states):
    verify_build(cache/'build')
    path = cache/'execution'
    record = json.loads((path/'run.json').read_text())
    if record['returncode'] != 0:
        raise ValueError('Failed trace retained; inspect run.log')
    if record['build_record_sha256'] != native.sha(cache/'build/build.json'):
        raise ValueError('Trace build identity differs')
    for name, digest in record['files_sha256'].items():
        if native.sha(path/name) != digest:
            raise ValueError('Trace artifact changed: '+name)
    if (path/'input.txt').read_text() != input_text(states):
        raise ValueError('Trace inputs differ from selected archive states')
    return record


def validate_trace(rows, states):
    """Reject stale particle inputs and missing/duplicate controls before analysis."""
    expected = {(r['case_id'], i) for r in states for i in range(len(CONTROLS))}
    keys = [(r['case_id'], r['control_id']) for r in rows]
    if len(set(keys)) != len(keys) or set(keys) != expected:
        raise ValueError('Missing, duplicate or unexpected trace rows')
    by_id = {r['case_id']: r for r in states}
    for r in rows:
        if not all(math.isfinite(v) for v in r.values()):
            raise ValueError('Nonfinite trace')
        state = by_id[r['case_id']]
        for key, expected_value in dict(temperature_K=state['temperature_K'], pressure_Pa=101000.,
                input_mass_kg=state['values']['mass_kg'], ice_volume_m3=state['values']['ice_volume_m3'],
                aspect_ratio=state['values']['aspect_ratio']).items():
            if not math.isclose(r[key], expected_value, rel_tol=2e-14):
                raise ValueError('Trace input mismatch: '+key)
        speeds = [r[k] for k in ('speed_from_reynolds_m_s', 'speed_boehm_m_s', 'terminal_speed_m_s')]
        if not all(math.isclose(v, speeds[0], rel_tol=2e-14) for v in speeds):
            raise ValueError('Native drag call paths disagree')
        a, c, volume = r['equatorial_dimension_m'], r['polar_dimension_m'], r['ice_volume_m3']
        if not (math.isclose(math.pi/6*a*a*c, volume, rel_tol=2e-14) and
                math.isclose(c/a, r['aspect_ratio'], rel_tol=2e-14) and
                math.isclose(r['diameter_m'], max(a, c), rel_tol=2e-14)):
            raise ValueError('Native spheroid mapping differs')
        control = int(r['control_id'])
        expected_mass = r['ice_density_kg_m3']*volume if control in (2, 3) else r['input_mass_kg']
        expected_area = math.pi/4*a*a if control in (1, 3) else r['native_projected_area_m2']
        if not (math.isclose(r['drag_mass_kg'], expected_mass, rel_tol=2e-14) and
                math.isclose(r['projected_area_m2'], expected_area, rel_tol=2e-14)):
            raise ValueError('Diagnostic intervention differs')


def figure_fields_summary(previous):
    """Reassess existing growth results without requiring the unplotted speed."""
    plotted = ('mass_kg', 'equatorial_dimension_m', 'polar_dimension_m', 'aspect_ratio')
    groups = {}
    for row in previous['comparisons']:
        if row['group'] != 'grid' or row['temperature_K'] >= 273.15:
            continue
        key = (row['variant'], row['mode'], row['elapsed_s'])
        group = groups.setdefault(key, dict(variant=key[0], mode=key[1], elapsed_s=key[2],
            available_endpoints=0, matches_plotted_mass_axes_and_aspect_ratio=0,
            matches_all_fields_except_speed=0))
        comparison = row['comparison']
        if comparison is None:
            continue
        group['available_endpoints'] += 1
        within = {k: v <= 1 for k, v in comparison['rounding_units'].items()}
        group['matches_plotted_mass_axes_and_aspect_ratio'] += int(all(within[k] for k in plotted))
        group['matches_all_fields_except_speed'] += int(all(within[k] for k in within if k != 'terminal_speed_m_s'))
    return list(groups.values())


def assemble(cache=CACHE):
    states, sources = selected_states()
    record = verify_run(cache, states)
    helper, helper_source = author_helper()
    with (cache/'execution/trace.csv').open() as f:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]
    validate_trace(rows, states)
    by_id = {r['case_id']: r for r in states}
    for r in rows:
        state = by_id[r['case_id']]
        args = (r['drag_mass_kg'], r['projected_area_m2'], r['diameter_m'], r['aspect_ratio'])
        speed = float(helper(*args, r['dynamic_viscosity_Pa_s'], r['air_density_kg_m3']))
        r['control'] = CONTROLS[int(r['control_id'])]
        r['author_helper_speed_m_s'] = speed
        r['helper_native_relative_difference'] = speed/r['terminal_speed_m_s']-1
        # These are fixed surface values in the authors' COLUMN plot, not known single-crystal inputs.
        r['author_helper_column_surface_air_speed_m_s'] = float(helper(*args, 1.715772e-5, 1.28780))
        r['archive_speed_m_s'] = state['values']['terminal_speed_m_s']
        r['native_archive_speed_ratio'] = r['terminal_speed_m_s']/r['archive_speed_m_s']
        r['within_archived_speed_rounding'] = abs(r['terminal_speed_m_s']-r['archive_speed_m_s']) <= state['rounding_half_width']['terminal_speed_m_s']
    previous_path = HERE/'results/benchmark_reproduction.json'
    previous = json.loads(previous_path.read_text())
    return dict(schema='terluna.atmosphere.mcsnow-drag-trace/1',
        producer={p.name: native.sha(p) for p in (Path(__file__), DRIVER, Path(model.__file__), Path(native.__file__))},
        archive_sources=sources, author_helper_source=helper_source,
        build=verify_build(cache/'build'), run_record=record,
        cache_path=str(cache.relative_to(model.ROOT)),
        states=states, traces=rows,
        maximum_helper_native_relative_difference=max(abs(r['helper_native_relative_difference']) for r in rows),
        prior_growth_product_sha256=native.sha(previous_path),
        prior_growth_figure_fields=figure_fields_summary(previous),
        evidence='Fixed-state diagnostic, with no growth integration. Native source is unchanged. The author helper is an independent archived implementation, not an empirical validation. Solid-envelope and one-variable changes are sensitivities, not identified original settings or adopted repairs. The archive lacks projected area and the original generating driver/configuration. Rounded mass and volume are used without clipping; selected compact rows can therefore slightly exceed the native solid-ice density and q=1.',
        figure_scope='Figures 10 and 13 plot mass, semiaxes, aspect ratio and clipped mass/volume density. Their scripts read but do not use vt. The four plotted primitive fields are a necessary subset of figure agreement; the additional all-except-speed count also checks volume and maximum dimension.',
        controls=dict(archived='Native geometry from archived m,V,phi, native moist air at water saturation and 101000 Pa.',
            area_only='Keep archived mass and dimensions; replace projected area by the full enclosing ellipse.',
            mass_only='Keep native projected area and dimensions; replace mass by native ice density times envelope volume.',
            solid_envelope='Replace both mass and area as above. This changes the archived mass.',
            column_surface_air='Authors fig16+18+20.py fixed rho0=1.28780 kg/m3, eta0=1.715772e-5 Pa s. Diagnostic alternative only; not provenance for Figures 10/13.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('run', 'assemble'))
    args = parser.parse_args()
    if args.action == 'run':
        run()
    else:
        model.save(PRODUCT, assemble())
