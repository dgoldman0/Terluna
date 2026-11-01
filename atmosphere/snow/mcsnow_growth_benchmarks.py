"""Compare pinned native McSnow growth packages and initial histories with lab fits.

The laboratory forcing is shared with SCALE. Individual raw runs and native builds
stay in the ignored run cache. Completed trajectories and trapped native failures
are retained separately; a partial trajectory is never an endpoint forecast.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import re
from pathlib import Path
from atmosphere.snow import mcsnow_crystal_growth as native

ROOT = native.ROOT
CACHE = ROOT / 'research/runs/snow_scale/mcsnow_growth'
RESULT = Path(__file__).with_name('results') / 'mcsnow_growth_benchmarks.json'
TEMPERATURES = (-3.7, -5.3, -10.6, -14.4, -20.1)
QUANTITIES = native.QUANTITIES


def case_name(version, mode, temperature, initialization, dt):
    return f'{version}_igf{mode}_{initialization}_T{abs(temperature):g}_dt{dt:g}'.replace('.', 'p')


def matrix():
    for version in ('2024', '2026'):
        for mode in (1, 2):
            for temperature in TEMPERATURES:
                for initialization, steps in (('fit', (1., .2)), ('sphere', (.2, .05))):
                    for dt in steps:
                        yield version, mode, temperature, initialization, dt


def relative_difference(coarse, fine):
    """Absolute relative changes at common output ages; no extrapolated endpoints."""
    by_age = {row['age_min']: row for row in fine['snapshots']}
    differences = {key: [] for key in QUANTITIES}
    common_ages = []
    for row in coarse['snapshots']:
        reference = by_age.get(row['age_min'])
        if reference is None:
            continue
        common_ages.append(row['age_min'])
        for key in QUANTITIES:
            differences[key].append(abs(row[key] / reference[key] - 1.))
    if not common_ages:
        return dict(both_complete=False, common_age_min=[], trajectory_maximum=None, last_common=None, endpoint=None)
    complete = coarse['complete'] and fine['complete']
    result = dict(both_complete=complete, common_age_min=common_ages,
                  trajectory_maximum={key: max(values) for key, values in differences.items()},
                  last_common={key: values[-1] for key, values in differences.items()})
    same_endpoint = (complete and coarse['snapshots'][-1]['age_min'] == fine['snapshots'][-1]['age_min']
                     and common_ages[-1] == coarse['snapshots'][-1]['age_min'])
    result['endpoint'] = result['last_common'] if same_endpoint else None
    return result


def verify_run(path, build_dir):
    """Re-read raw CSVs and verify inputs, source, executable and derived snapshots."""
    path, build_dir = Path(path), Path(build_dir)
    record = json.loads((path / 'run.json').read_text())
    build = native.verify_build(build_dir)
    if record['build_identity'] != build or record['build_record_sha256'] != native.sha(build_dir / 'build.json'):
        raise ValueError('Build identity differs from execution record')
    if record['producer'] != native.producer():
        raise ValueError('Growth producer changed')
    for group in ('inputs', 'outputs'):
        for name, digest in record[group].items():
            if native.sha(path / name) != digest:
                raise ValueError(f'Changed {group} file: {name}')
    config = record['configuration']
    constants = (build_dir / 'source/src/mo_atmo_types.f90').read_text()
    density_match = re.search(r'\brhoi\s*=\s*([0-9.]+)_wp', constants, re.I)
    if not density_match:
        raise ValueError('Native ice density not found')
    density = float(density_match[1])
    regenerated = native.configuration(record['temperature_C'], config['dt_s'], config['growth_mode'],
        record['initialization'], ice_density=density)
    start = native.initial(record['temperature_C'], record['initialization'], density)
    if config != regenerated or record['initial'] != start or (path / 'crystal.nml').read_text() != native.namelist(regenerated):
        raise ValueError('Physics inputs do not reproduce')
    # A native crash can leave absent files or a truncated terminal row. Recompute
    # the diagnostic classification and retain its explicit errors and valid prefix.
    analyzed = native.analyze_output(path, config, start, record['temperature_C'], config['growth_mode'], record['returncode'])
    for key, value in analyzed.items():
        if record[key] != value:
            raise ValueError('Output classification or analysis differs: ' + key)
    rows = []
    if (path / 'trajectory.csv').exists():
        with (path / 'trajectory.csv').open() as stream:
            rows = list(csv.DictReader(stream))
    if len(rows) < len(record['snapshots']):
        raise ValueError('Raw trajectory lacks recorded valid prefix')
    comparisons = 0
    for raw, saved in zip(rows, record['snapshots']):
        for key, text in raw.items():
            value = float(text)
            if not math.isfinite(value) or value != saved[key]:
                raise ValueError('Raw snapshot differs: ' + key)
            comparisons += 1
        for key in QUANTITIES:
            if saved[key] <= 0:
                raise ValueError('Nonpositive ice state')
        volume = math.pi / 6 * saved['equatorial_dimension_m']**2 * saved['polar_dimension_m']
        if not math.isclose(volume, saved['ice_volume_m3'], rel_tol=2e-12):
            raise ValueError('Axis/volume mapping failed')
        if not math.isclose(saved['diameter_m'], max(saved['equatorial_dimension_m'], saved['polar_dimension_m']), rel_tol=2e-12):
            raise ValueError('Maximum dimension differs from axes')
        if not math.isclose(saved['aspect_ratio'], saved['polar_dimension_m'] / saved['equatorial_dimension_m'], rel_tol=2e-12):
            raise ValueError('Aspect ratio differs from axes')
        if not math.isclose(saved['density_kg_m3'], saved['mass_kg'] / volume, rel_tol=2e-12):
            raise ValueError('Density differs from mass and axes')
        if saved['age_min'] >= 3:
            observed = native.lab.observed(record['temperature_C'], saved['age_min'])
            if saved['observed_fit'] != observed:
                raise ValueError('Observation fit differs')
            for key in QUANTITIES:
                if saved['ratio_to_fit'][key] != saved[key] / observed[key]:
                    raise ValueError('Derived ratio differs')
    expected_end = native.lab.FITS[record['temperature_C']][0]
    if record['complete'] != (record['returncode'] == 0 and not record['output_validation_errors'] and bool(record['snapshots']) and record['snapshots'][-1]['age_min'] == expected_end):
        raise ValueError('Completion incorrectly classified')
    return dict(raw_snapshot_values_checked=comparisons,
                maximum_relative_ice_deposition_residual=record['maximum_relative_ice_deposition_residual'])


def summarize(record, path):
    snapshots = record['snapshots']
    complete = record['complete']
    return dict(run=str(Path(path).relative_to(ROOT)), run_sha256=native.sha(Path(path) / 'run.json'),
        version=record['build_identity']['source_identity']['version'], temperature_C=record['temperature_C'],
        growth_mode=record['configuration']['growth_mode'], growth_package=record['growth_package'],
        initialization=record['initialization'], dt_s=record['configuration']['dt_s'],
        complete=complete, returncode=record['returncode'], output_validation_errors=record['output_validation_errors'],
        requested_end_age_min=native.lab.FITS[record['temperature_C']][0],
        last_saved=snapshots[-1] if snapshots else None, endpoint=snapshots[-1] if complete else None,
        age_three=next((row for row in snapshots if row['age_min'] == 3.), None),
        age_ten=next((row for row in snapshots if row['age_min'] == 10.), None),
        branch_diagnostic=record['branch_diagnostic'], step_diagnostic=record['step_diagnostic'],
        air=record['air'], maximum_relative_ice_deposition_residual=record['maximum_relative_ice_deposition_residual'])


def assemble(cache=CACHE):
    cache = Path(cache)
    runs = {}; cases = []; checks = []
    for version, mode, temperature, initialization, dt in matrix():
        name = case_name(version, mode, temperature, initialization, dt)
        path = cache / 'runs' / name
        checks.append(verify_run(path, cache / ('build_' + version)))
        record = json.loads((path / 'run.json').read_text())
        runs[(version, mode, temperature, initialization, dt)] = record
        cases.append(summarize(record, path))
    refinements = []
    for version in ('2024', '2026'):
        for mode in (1, 2):
            for temperature in TEMPERATURES:
                for initialization, steps in (('fit', (1., .2)), ('sphere', (.2, .05))):
                    coarse, fine = (runs[version, mode, temperature, initialization, dt] for dt in steps)
                    refinements.append(dict(version=version, growth_mode=mode, temperature_C=temperature,
                        initialization=initialization, timesteps_s=steps, **relative_difference(coarse, fine)))
    completed = sum(case['complete'] for case in cases)
    return dict(schema='terluna.atmosphere.mcsnow-growth-benchmarks/1',
        evidence='Native Earth process comparisons with laboratory fits and measured timestep sensitivity; no lunar calibration.',
        producer={**native.producer(), Path(__file__).name: native.sha(__file__)},
        builds={v: json.loads((cache / ('build_' + v) / 'build.json').read_text()) for v in ('2024', '2026')},
        matrix=dict(versions=['2024', '2026'], growth_modes=[1, 2], temperatures_C=TEMPERATURES,
                    fit_timesteps_s=[1, .2], sphere_timesteps_s=[.2, .05],
                    fit_start_age_min=3, sphere_start_age_min=0, sphere_radius_m=2e-6,
                    executions=len(cases), complete=completed, failed=len(cases)-completed),
        verification=dict(executions=len(checks), raw_snapshot_values_checked=sum(c['raw_snapshot_values_checked'] for c in checks),
            inputs_regenerated=True, source_and_output_hashes_verified=True,
            maximum_relative_ice_deposition_residual=max((c['maximum_relative_ice_deposition_residual'] for c in checks if c['maximum_relative_ice_deposition_residual'] is not None), default=None)),
        cases=cases, refinements=refinements,
        limits=['All cases use fixed Earth laboratory air; no aggregate collisions, fracture or lunar adaptation.',
                'The same Takahashi evidence informed native growth choices; these are not independent validation cases.',
                'Fit starts impose a spheroid representation of empirical axes and mass; sphere starts test a different history.',
                'Failed runs have no endpoint. Last saved states only describe the failure trajectory.',
                'Timestep differences measure numerical sensitivity; they do not establish an asymptotic convergence order.',
                'The 2026 archive lacks habit tables; both versions use the explicitly pinned 2024 tables.',
                'The presently identified native source is not an exact reproduction of archived paper curves.'])


def refinement_cases():
    """Follow up sensitivity observed in the completed 80-case matrix."""
    yield '2026', 1, -5.3, 'sphere', .0125
    yield '2026', 1, -14.4, 'sphere', .0125
    yield '2026', 2, -14.4, 'sphere', .0125
    yield '2026', 1, -14.4, 'fit', .05
    yield '2026', 2, -14.4, 'fit', .05


def refine(cache=CACHE):
    cache = Path(cache)
    cases = []; comparisons = []; checks = []
    for version, mode, temperature, initialization, dt in refinement_cases():
        name = case_name(version, mode, temperature, initialization, dt)
        path = cache / 'runs' / name
        if path.exists():
            record = json.loads((path / 'run.json').read_text())
        else:
            record = native.run(cache / ('build_' + version), path, temperature, dt, mode, initialization)
        checks.append(verify_run(path, cache / ('build_' + version)))
        cases.append(summarize(record, path))
        coarse_dt = .05 if initialization == 'sphere' else .2
        coarse_path = cache / 'runs' / case_name(version, mode, temperature, initialization, coarse_dt)
        verify_run(coarse_path, cache / ('build_' + version))
        coarse = json.loads((coarse_path / 'run.json').read_text())
        comparisons.append(dict(version=version, growth_mode=mode, temperature_C=temperature,
            initialization=initialization, timesteps_s=[coarse_dt, dt], **relative_difference(coarse, record)))
        print(name, 'complete' if record['complete'] else 'FAILED', flush=True)
    return dict(schema='terluna.atmosphere.mcsnow-growth-refinements/1',
        evidence='Additional timestep checks selected after observed sensitivity in the 80-case matrix.',
        producer={**native.producer(), Path(__file__).name: native.sha(__file__)},
        cases=cases, refinements=comparisons, verification=checks,
        limits=['Selected refinements, not a second complete factorial matrix.',
                'Failed trajectories are compared only at common saved ages; no endpoint or convergence is claimed.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=CACHE)
    parser.add_argument('--output', type=Path, default=RESULT)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--refine', action='store_true')
    args = parser.parse_args()
    if args.run and args.refine:
        parser.error('Choose the main matrix or supplemental refinements')
    if args.refine:
        result = refine(args.cache)
        path = args.output.with_name('mcsnow_growth_refinements.json')
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    if args.run:
        for version, mode, temperature, initialization, dt in matrix():
            name = case_name(version, mode, temperature, initialization, dt)
            path = args.cache / 'runs' / name
            if path.exists():
                raise FileExistsError(path)
            record = native.run(args.cache / ('build_' + version), path, temperature, dt, mode, initialization)
            print(name, 'complete' if record['complete'] else 'FAILED',
                  'last_age', record['snapshots'][-1]['age_min'] if record['snapshots'] else None, flush=True)
    product = assemble(args.cache)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=2, allow_nan=False) + '\n')
    print(json.dumps(product['matrix'], indent=2))


if __name__ == '__main__':
    main()
