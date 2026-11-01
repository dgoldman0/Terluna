"""Guard empirical/native units, source identity, and preservation of failed runs.

All subprocess outcomes below are fixtures. These tests neither compile nor run
McSnow and do not establish numerical convergence or empirical agreement.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest
from atmosphere.snow import mcsnow_crystal_growth as growth


def test_fit_initialization_preserves_full_lab_axes_and_absolute_age():
    state = growth.initial(-14.4, 'fit', 919.)
    config = growth.configuration(-14.4, .2, 2, 'fit', ice_density=919.)
    # Bounds independently catch mm/m, g/kg, radius/full-axis and age-zero errors.
    assert .397e-3 < config['initial_equatorial_m'] < .398e-3
    assert 7.78e-6 < config['initial_polar_m'] < 7.80e-6
    assert 3.83e-10 < config['initial_mass_kg'] < 3.84e-10
    assert state['age_min'] == 3
    assert config['duration_s'] == 27 * 60
    assert config['temperature_K'] == pytest.approx(258.75)
    assert config['pressure_Pa'] == 101000.
    assert config['prolate_interpolation'] is True


def test_author_seed_has_two_micrometre_radius_and_starts_at_zero():
    seed = growth.initial(-14.4, 'sphere', 919.)
    config = growth.configuration(-14.4, .05, 2, 'sphere', ice_density=919.)
    assert seed['age_min'] == 0
    assert seed['equatorial_dimension_m'] == seed['polar_dimension_m'] == 4e-6
    assert 3.07e-14 < seed['mass_kg'] < 3.09e-14
    assert seed['density_kg_m3'] == 919.
    assert config['duration_s'] == 30 * 60


@pytest.mark.parametrize('temperature', [-8.6, -12.2, -16.5, -18.2, -22.])
def test_inadmissible_fitted_spheroid_is_rejected_without_clipping(temperature):
    with pytest.raises(ValueError, match='above bulk ice'):
        growth.initial(temperature, 'fit', 919.)


@pytest.mark.parametrize('step', [0., -1., float('nan'), float('inf'), .7])
def test_invalid_sampling_timestep_is_rejected(step):
    with pytest.raises(ValueError, match='timestep'):
        growth.configuration(-14.4, step, 1, 'fit', ice_density=919.)


def _identity(data, name):
    return dict(path='input/' + name, bytes=len(data),
                sha256=hashlib.sha256(data).hexdigest(),
                git_blob_sha1=hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest())


def _build_fixture(path):
    (path / 'source/src').mkdir(parents=True)
    (path / 'source/src/mo_atmo_types.f90').write_text('REAL(wp), PARAMETER :: rhoi = 919_wp\n')
    (path / 'crystal').write_bytes(b'identified executable fixture')
    (path / growth.DRIVER.name).write_bytes(b'identified caller fixture')
    identities = {}
    for name in growth.TABLE_NAMES:
        data = (name + ': identified table fixture').encode()
        (path / name).write_bytes(data)
        identities[name] = _identity(data, name)
    record = dict(executable_sha256=growth.sha(path / 'crystal'),
                  driver_sha256=growth.sha(path / growth.DRIVER.name),
                  source_files={'src/mo_atmo_types.f90': growth.sha(path / 'source/src/mo_atmo_types.f90')},
                  table_identity=identities)
    (path / 'build.json').write_text(json.dumps(record))
    return record


@pytest.mark.parametrize('relative_path', [
    'crystal', growth.DRIVER.name, 'source/src/mo_atmo_types.f90', growth.TABLE_NAMES[0],
])
def test_changed_executable_caller_source_or_table_cannot_be_reused(tmp_path, relative_path):
    record = _build_fixture(tmp_path)
    assert growth.verify_build(tmp_path) == record
    with (tmp_path / relative_path).open('ab') as stream:
        stream.write(b'changed')
    with pytest.raises(ValueError):
        growth.verify_build(tmp_path)


@pytest.fixture
def recorded_native_outcome(tmp_path, monkeypatch):
    """Install a fake process with well-defined outputs; no numerical model runs."""
    build_dir = tmp_path / 'build'
    build_record = _build_fixture(build_dir)
    build_record['producer'] = {'unit_test_fixture': '1'}
    (build_dir / 'build.json').write_text(json.dumps(build_record))
    monkeypatch.setattr(growth, 'producer', lambda: {'unit_test_fixture': '1'})
    cache = tmp_path / 'cache'
    source = cache / growth.scale.SOURCE_PDF
    source.parent.mkdir(parents=True)
    source.write_bytes(b'identified empirical source fixture')
    monkeypatch.setattr(growth.scale.lunar, 'CACHE', cache)
    monkeypatch.setattr(growth.lab, 'SOURCE_SHA', growth.sha(source))

    def install(returncode, elapsed_times, *, wrong_initial_dimension=False):
        def process(args, *, cwd, **kwargs):
            seed = growth.initial(-14.4, 'fit', 919.)
            if wrong_initial_dimension:
                seed['equatorial_dimension_m'] *= 2
                seed['diameter_m'] *= 2
                seed['density_kg_m3'] /= 4
            trajectory = []
            for elapsed in elapsed_times:
                trajectory.append(dict(elapsed_s=elapsed, mass_kg=seed['mass_kg'],
                    equatorial_dimension_m=seed['equatorial_dimension_m'],
                    polar_dimension_m=seed['polar_dimension_m'], diameter_m=seed['diameter_m'],
                    density_kg_m3=seed['density_kg_m3'],
                    ice_volume_m3=seed['mass_kg']/seed['density_kg_m3'],
                    aspect_ratio=seed['polar_dimension_m']/seed['equatorial_dimension_m'],
                    terminal_speed_m_s=.1, ice_mass_balance_residual_kg=0.))
            for filename, rows in [('trajectory.csv', trajectory), ('air.csv', [{k:v for k,v in growth.configuration(-14.4,.2,1,'fit',ice_density=919.).items() if k in ('temperature_K','pressure_Pa','vapor_pressure_Pa')}])]:
                with (Path(cwd) / filename).open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                    writer.writeheader(); writer.writerows(rows)
            return SimpleNamespace(returncode=returncode)
        monkeypatch.setattr(growth.subprocess, 'run', process)
        return build_dir, tmp_path / 'run'
    return install


def test_native_failure_retains_partial_valid_trajectory_and_returncode(recorded_native_outcome):
    build_dir, output = recorded_native_outcome(-8, [0., 60.])
    result = growth.run(build_dir, output, -14.4, .2, 1, 'fit')
    assert result['complete'] is False
    assert result['returncode'] == -8
    assert [row['age_min'] for row in result['snapshots']] == [3., 4.]
    assert result['step_diagnostic'] is None
    assert 'trajectory.csv' in result['outputs']
    assert 'run.log' in result['outputs']
    assert json.loads((output / 'run.json').read_text())['complete'] is False


def test_zero_exit_with_missing_expected_snapshots_is_not_success(recorded_native_outcome):
    build_dir, output = recorded_native_outcome(0, [0., 60.])
    result = growth.run(build_dir, output, -14.4, .2, 1, 'fit')
    assert result['complete'] is False
    assert result['returncode'] == 0
    assert 'Successful native exit with incomplete trajectory' in result['output_validation_errors']


def test_bad_initial_mapping_is_not_misclassified_as_physical_runaway(recorded_native_outcome):
    build_dir, output = recorded_native_outcome(-8, [0.], wrong_initial_dimension=True)
    result = growth.run(build_dir, output, -14.4, .2, 1, 'fit')
    assert result['complete'] is False
    assert any('Initial mapping differs' in error for error in result['output_validation_errors'])


def test_failed_run_directory_cannot_be_silently_overwritten(recorded_native_outcome):
    build_dir, output = recorded_native_outcome(-8, [0.])
    growth.run(build_dir, output, -14.4, .2, 1, 'fit')
    before = (output / 'run.json').read_bytes()
    with pytest.raises(FileExistsError):
        growth.run(build_dir, output, -14.4, .2, 1, 'fit')
    assert (output / 'run.json').read_bytes() == before


def test_failure_before_first_output_is_retained(tmp_path):
    start = growth.initial(-14.4, 'fit', 919.)
    config = growth.configuration(-14.4, .2, 1, 'fit', ice_density=919.)
    result = growth.analyze_output(tmp_path, config, start, -14.4, 1, -8)
    assert result['complete'] is False
    assert result['returncode'] == -8
    assert result['snapshots'] == []
    assert result['air'] == {}
    assert 'Missing trajectory.csv' in result['output_validation_errors']
    assert 'No valid crystal snapshots' in result['output_validation_errors']


@pytest.mark.parametrize('last_line,error', [
    ('1,nan\n', 'Nonfinite'),
    ('1,2,3\n', 'Truncated or malformed'),
    ('1,not-a-number\n', 'Nonnumeric'),
])
def test_invalid_terminal_csv_row_preserves_valid_prefix(tmp_path, last_line, error):
    path = tmp_path / 'partial.csv'
    path.write_text('elapsed_s,mass_kg\n0,1e-10\n' + last_line)
    rows, errors = growth.read_csv(path, ('elapsed_s', 'mass_kg'))
    assert rows == [dict(elapsed_s=0., mass_kg=1e-10)]
    assert len(errors) == 1
    assert error in errors[0]


def test_duplicate_header_is_rejected(tmp_path):
    path = tmp_path / 'duplicate.csv'
    path.write_text('mass_kg,mass_kg\n1,2\n')
    rows, errors = growth.read_csv(path, ('mass_kg',))
    assert rows == []
    assert errors == ['Invalid header in duplicate.csv']


def _refinement_record(ages, values, complete):
    from atmosphere.snow import mcsnow_growth_benchmarks as benchmark
    return dict(complete=complete, snapshots=[dict(age_min=age, **{key:value for key in benchmark.QUANTITIES})
        for age, value in zip(ages, values)])


def test_refinement_compares_only_common_ages_and_retains_failure_status():
    from atmosphere.snow import mcsnow_growth_benchmarks as benchmark
    coarse = _refinement_record([3.,4.,5.], [1.,4.,100.], True)
    fine = _refinement_record([3.,4.], [2.,2.], False)
    result = benchmark.relative_difference(coarse, fine)
    assert result['common_age_min'] == [3.,4.]
    assert result['both_complete'] is False
    assert result['endpoint'] is None
    assert all(value == 1. for value in result['trajectory_maximum'].values())
    assert all(value == 1. for value in result['last_common'].values())


def test_native_failure_at_final_saved_age_is_still_not_an_endpoint():
    from atmosphere.snow import mcsnow_growth_benchmarks as benchmark
    coarse = _refinement_record([3.,30.], [1.,3.], False)
    fine = _refinement_record([3.,30.], [1.,4.], True)
    result = benchmark.relative_difference(coarse, fine)
    assert result['endpoint'] is None
    assert result['both_complete'] is False


def test_complete_refinement_requires_the_same_terminal_age():
    from atmosphere.snow import mcsnow_growth_benchmarks as benchmark
    coarse = _refinement_record([3.,4.], [1.,3.], True)
    fine = _refinement_record([3.,4.,5.], [1.,4.,10.], True)
    result = benchmark.relative_difference(coarse, fine)
    assert result['both_complete'] is True
    assert result['endpoint'] is None


def test_complete_shared_endpoint_is_reported_separately_from_maximum_change():
    from atmosphere.snow import mcsnow_growth_benchmarks as benchmark
    coarse = _refinement_record([3.,4.], [1.,3.], True)
    fine = _refinement_record([3.,4.], [2.,4.], True)
    result = benchmark.relative_difference(coarse, fine)
    assert all(value == .5 for value in result['trajectory_maximum'].values())
    assert all(value == .25 for value in result['endpoint'].values())


def test_no_common_age_does_not_invent_a_refinement():
    from atmosphere.snow import mcsnow_growth_benchmarks as benchmark
    result = benchmark.relative_difference(_refinement_record([3.], [1.], False),
                                           _refinement_record([4.], [1.], False))
    assert result['common_age_min'] == []
    assert result['trajectory_maximum'] is None
    assert result['last_common'] is None
    assert result['endpoint'] is None


def test_failed_summary_retains_last_saved_state_without_endpoint(tmp_path, monkeypatch):
    from atmosphere.snow import mcsnow_growth_benchmarks as benchmark
    monkeypatch.setattr(benchmark, 'ROOT', tmp_path)
    output = tmp_path / 'failed'; output.mkdir()
    (output / 'run.json').write_text('{}\n')
    record = dict(snapshots=[dict(age_min=30.)], complete=False, returncode=-8,
                  build_identity={'source_identity':{'version':'2026'}}, temperature_C=-14.4,
                  configuration={'growth_mode':2,'dt_s':.2}, growth_package='revised',
                  initialization='fit', output_validation_errors=[], branch_diagnostic=[],
                  step_diagnostic=None, air={}, maximum_relative_ice_deposition_residual=0.)
    summary = benchmark.summarize(record, output)
    assert summary['last_saved'] == {'age_min':30.}
    assert summary['endpoint'] is None
    assert summary['requested_end_age_min'] == 30.
    assert summary['returncode'] == -8
