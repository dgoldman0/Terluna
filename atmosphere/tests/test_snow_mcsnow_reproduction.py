"""Interface checks; synthetic subprocess fixtures are not physics validation."""
import csv
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest
from atmosphere.snow import mcsnow_reproduction as model


def author_line(mass=1e-9):
    # Known spheroid with 100 and 50 micrometre semiaxes.
    volume = 4*math.pi/3*1e-4**2*5e-5
    return f'600.0,258.650,0.00010000,0.00005000,0.50000000,{mass:.15f},{volume:.18f},0.2000,0.00020000\n'.encode()


def test_author_columns_axes_and_decimal_uncertainties():
    row = model.parse_author_rows(author_line(), 600)[0]
    assert row['temperature_K'] == 258.65
    assert row['values']['equatorial_dimension_m'] == 200e-6
    assert row['values']['polar_dimension_m'] == 100e-6
    assert row['values']['mass_kg'] == 1e-9
    assert row['rounding_half_width']['equatorial_dimension_m'] == 1e-8
    assert row['rounding_half_width']['mass_kg'] == 5e-16
    assert row['rounding_half_width']['terminal_speed_m_s'] == 5e-5


def test_plot_density_clipping_does_not_change_raw_mass_or_density():
    row = model.parse_author_rows(author_line(3e-9), 600)[0]
    assert row['raw_density_kg_m3'] > 1400
    assert row['plotted_density_kg_m3'] == 919
    assert row['values']['mass_kg'] == 3e-9


@pytest.mark.parametrize('raw,match', [
    (author_line().replace(b'600.0,', b''), 'nine|eight'),
    (author_line().replace(b'600.0,', b'180.0,'), 'Malformed'),
    (author_line().replace(b'0.00020000', b'0.00040000'), 'geometry'),
    (author_line()+author_line(), 'duplicate'),
    (author_line().replace(b'258.650', b'nan'), 'Malformed'),
])
def test_wrong_column_count_time_geometry_or_duplicate_is_rejected(raw, match):
    with pytest.raises(ValueError, match=match):
        model.parse_author_rows(raw, 600)


def test_corrupt_author_archive_is_rejected_before_parsing(tmp_path):
    (tmp_path/'McSnow_paper_plot_data.zip').write_bytes(b'not the identified archive')
    with pytest.raises(ValueError, match='checksum mismatch'):
        model.author_data(tmp_path)


def test_precision_match_does_not_hide_a_wrong_axis_conversion():
    author = model.parse_author_rows(author_line(), 600)[0]
    state = dict(author['values'])
    assert model.compare(state, author)['all_within_printed_rounding']
    state['equatorial_dimension_m'] /= 2
    assert not model.compare(state, author)['all_within_printed_rounding']


def test_historical_patch_preserves_directional_growth_and_rejects_unknown_source(monkeypatch):
    with pytest.raises(ValueError, match='Unrecognized'):
        model.deposition_candidate('unknown source', historical=True)
    body = ('    Dv = Dv*ventilation_coef_i( sqrtRe*cubrSc, sp%phi )  ! CHRS for old bugged version (reproduce paper) comment these two lines\n'
            '    Ka = Ka*ventilation_coef_i( sqrtRe*cubrPr, sp%phi )\n'
            '  Hi = (L/(Rv*atmo%T) - 1.0_wp)*(L*Dv*es)/(Ka*Rv*atmo%T**2)\n'
            '  call habit_predict(sp, directional_ratio)\n')
    monkeypatch.setattr(model, 'DEPOSITION_SHA', hashlib.sha256(body.encode()).hexdigest())
    changed = model.deposition_candidate(body, historical=True)
    active = '\n'.join(line for line in changed.splitlines() if not line.lstrip().startswith('!'))
    assert 'ventilation_coef_i' not in active
    assert 'Hi = (L/(Rv*atmo%T) - 1.0_wp)' in active
    assert 'call habit_predict(sp, directional_ratio)' in active
    heat = model.deposition_candidate(body, leading_heat=True)
    assert heat.count('ventilation_coef_i') == 2
    assert 'Hi = L**2*Dv*es/(Ka*Rv**2*atmo%T**3)' in heat


def make_outputs(path, config, times, wrong_axis=False):
    volume = math.pi/6*config['initial_equatorial_m']**2*config['initial_polar_m']
    rows = [dict(elapsed_s=t, mass_kg=volume*919,
                 equatorial_dimension_m=config['initial_equatorial_m']*(2 if wrong_axis else 1),
                 polar_dimension_m=config['initial_polar_m'], diameter_m=4e-6,
                 density_kg_m3=919, ice_volume_m3=volume, aspect_ratio=1,
                 terminal_speed_m_s=.01, ice_mass_balance_residual_kg=0) for t in times]
    air = [dict(temperature_K=config['temperature_K'], pressure_Pa=config['pressure_Pa'],
                rh_percent=100, vapor_pressure_Pa=config['vapor_pressure_Pa'], native_ice_density_kg_m3=919)]
    for name, values in [('trajectory.csv', rows), ('air.csv', air)]:
        with (path/name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]))
            writer.writeheader(); writer.writerows(values)


@pytest.fixture
def fake_run(tmp_path, monkeypatch):
    build = tmp_path/'build'; build.mkdir()
    (build/'build.json').write_text('{}')
    monkeypatch.setattr(model, 'verify_build', lambda path: dict(historical_mean_ventilation=False))
    monkeypatch.setattr(model.native, 'TABLE_NAMES', [])
    config = model.configuration(258.65, duration=90)

    def execute(returncode, times, wrong_axis=False):
        def process(args, *, cwd, **kwargs):
            make_outputs(Path(cwd), config, times, wrong_axis)
            return SimpleNamespace(returncode=returncode)
        monkeypatch.setattr(model.subprocess, 'run', process)
        output = tmp_path/'run'
        return output, model.run(build, output, config)
    return execute


def test_final_partial_output_interval_is_retained(fake_run):
    output, record = fake_run(0, [0, 60, 90])
    assert record['complete']
    assert model.verify_run(output)['snapshots'][-1]['elapsed_s'] == 90


def test_failed_execution_preserves_partial_trajectory_and_exit_code(fake_run):
    output, record = fake_run(-8, [0, 60])
    assert not record['complete']
    assert record['returncode'] == -8
    assert len(record['snapshots']) == 2
    assert model.verify_run(output) == record
    assert 'trajectory.csv' in record['file_sha256']


def test_zero_exit_with_missing_snapshot_is_not_complete(fake_run):
    _, record = fake_run(0, [0, 60])
    assert not record['complete']
    assert 'Successful exit with incomplete trajectory' in record['validation_errors']


def test_bad_initial_mapping_is_distinct_from_growth_failure(fake_run):
    _, record = fake_run(0, [0, 60, 90], wrong_axis=True)
    assert not record['complete']
    assert 'Initial mass/axis mapping mismatch' in record['validation_errors']


@pytest.mark.parametrize('tamper', ['raw', 'summary', 'configuration'])
def test_replay_rejects_altered_raw_summary_or_input(fake_run, tamper):
    output, record = fake_run(0, [0, 60, 90])
    if tamper == 'raw':
        with (output/'trajectory.csv').open('a') as stream:
            stream.write('changed\n')
    elif tamper == 'summary':
        record['snapshots'][-1]['mass_kg'] *= 2
        model.save(output/'run.json', record)
    else:
        record['configuration']['pressure_Pa'] *= .5
        model.save(output/'run.json', record)
    with pytest.raises(ValueError):
        model.verify_run(output)


def test_repository_product_records_actual_parity_and_no_claim_of_exact_reproduction():
    product = json.loads((model.ROOT/'atmosphere/snow/results/benchmark_reproduction.json').read_text())
    assert all(v['byte_identical'] for v in product['parity'].values())
    assert product['counts']['executions'] == len(product['manifest'])
    assert all(s['matches_printed_precision'] == 0 for s in product['grid_summary'])
    assert 'unresolved' in product['evidence']
    for extension in product['extensions']:
        for point in extension['points']:
            if point['laboratory_fit'] is not None:
                from atmosphere.snow import takahashi_growth as lab
                assert 3 <= point['age_min'] <= lab.FITS[extension['temperature_C']][0]
