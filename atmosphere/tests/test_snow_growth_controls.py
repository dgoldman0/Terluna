"""Guard source interventions and arbitrary-state growth diagnostics.

Subprocess outputs here are deliberate fixtures; tests do not run native physics,
validate a parameterization, or replace numerical parity/refinement records.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest
from atmosphere.snow import mcsnow_ventilation_controls as controls


def test_patch_identity_agrees_with_registered_2026_source():
    register=json.loads(controls.native.MANIFEST.read_text())['mcsnow_stochastic_geometry']
    selected={item['path']:item for item in register['selected_source_files']}
    assert controls.NATIVE_SHA==selected[controls.TARGET]['sha256']
    with pytest.raises(ValueError,match='Unrecognized 2026'):
        controls.patch_source('unidentified or altered deposition implementation')


def test_missing_or_ambiguous_anchor_is_rejected_before_substitution():
    for text in ('target omitted','anchor anchor'):
        with pytest.raises(ValueError,match='Expected one patch anchor'):
            controls.replace_once(text,'anchor','replacement')


def test_patch_preserves_native_mean_directional_and_capacitance_implementations():
    source=controls.native.ROOT/'research/runs/snow_scale/mcsnow_growth/build_2026/source'/controls.TARGET
    if not source.exists():pytest.skip('Optional pinned external McSnow source has not been restored')
    original=source.read_text();original_bytes=source.read_bytes()
    patched=controls.patch_source(original)
    # These mathematical functions remain byte-for-byte native; the caller
    # selects interventions separately. Their calibration constants are intact.
    for opening,closing in [
        ('PURE ELEMENTAL FUNCTION ventilation_coef_i(', 'END FUNCTION ventilation_coef_i'),
        ('PURE ELEMENTAL FUNCTION ventilation_coef_IGF(', 'END FUNCTION ventilation_coef_IGF'),
        ('FUNCTION capacitance(', 'END FUNCTION capacitance'),
        ('SUBROUTINE vapourDiffusionTemperature(', 'END SUBROUTINE vapourDiffusionTemperature'),
    ]:
        before=original[original.index(opening):original.index(closing)+len(closing)]
        after=patched[patched.index(opening):patched.index(closing)+len(closing)]
        assert after==before
    native_direction=original[original.index('  IGF_vent = IGF *'):original.index('  expo     =')].strip()
    assert native_direction in patched
    assert source.read_bytes()==original_bytes


def test_new_upstream_duplicate_anchor_is_rejected_even_with_updated_pin(monkeypatch):
    source=controls.native.ROOT/'research/runs/snow_scale/mcsnow_growth/build_2026/source'/controls.TARGET
    if not source.exists():pytest.skip('Optional pinned external McSnow source has not been restored')
    original=source.read_text()
    changed=original.replace('  C = capacitance( sp, Re )','  C = capacitance( sp, Re )\n  C = capacitance( sp, Re )')
    monkeypatch.setattr(controls,'NATIVE_SHA',hashlib.sha256(changed.encode()).hexdigest())
    with pytest.raises(ValueError,match='Expected one patch anchor'):
        controls.patch_source(changed)


@pytest.mark.parametrize('mean,directional',[(4,0),(-1,0),(0,3),(0,-1),(1.5,0)])
def test_invalid_modes_are_rejected(mean,directional):
    with pytest.raises(ValueError,match='Unknown ventilation intervention'):
        controls.validate_modes(mean,directional)


def _write_csv(path,rows):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


@pytest.fixture
def fixture_run(tmp_path,monkeypatch):
    """Supply an identified fake process; native model builds are never invoked."""
    build_dir=tmp_path/'build';(build_dir/'source/src').mkdir(parents=True)
    (build_dir/'source/src/mo_atmo_types.f90').write_text('REAL(wp), PARAMETER :: rhoi = 919_wp\n')
    identity={'producer':{'unit_test_fixture':'1'}}
    (build_dir/'build.json').write_text(json.dumps(identity))
    for name in controls.native.TABLE_NAMES:(build_dir/name).write_text('test table '+name)
    monkeypatch.setattr(controls,'producer',lambda:identity['producer'])
    monkeypatch.setattr(controls,'verify_build',lambda path:identity)
    cache=tmp_path/'cache';source=cache/controls.native.scale.SOURCE_PDF
    source.parent.mkdir(parents=True);source.write_bytes(b'empirical source test fixture')
    monkeypatch.setattr(controls.native.scale.lunar,'CACHE',cache)
    monkeypatch.setattr(controls.native.lab,'SOURCE_SHA',controls.native.sha(source))

    def install(config,*,returncode=0,missing_trace=False,negative_ventilation=False,
                wrong_initial_axes=False,wrong_trace_time=False,truncate=False):
        config=dict(config)
        def process(args,*,cwd,**kwargs):
            directory=Path(cwd)
            timestep=config['dt_s'];nsteps=round(config['duration_s']/timestep)
            outsteps=round(config['output_interval_s']/timestep)
            saved=[0,*range(outsteps,nsteps+1,outsteps)]
            if saved[-1]!=nsteps:saved.append(nsteps)
            if truncate:saved=saved[:2]
            a=config['initial_equatorial_m']*(2 if wrong_initial_axes else 1)
            c=config['initial_polar_m'];m=config['initial_mass_kg']
            volume=math.pi*a*a*c/6
            rows=[dict(elapsed_s=step*timestep,mass_kg=m,equatorial_dimension_m=a,
                polar_dimension_m=c,diameter_m=max(a,c),density_kg_m3=m/volume,
                ice_volume_m3=volume,aspect_ratio=c/a,terminal_speed_m_s=.2,
                ice_mass_balance_residual_kg=0.) for step in saved]
            _write_csv(directory/'trajectory.csv',rows)
            _write_csv(directory/'air.csv',[{k:config[k] for k in ('temperature_K','pressure_Pa','vapor_pressure_Pa')}])
            trace_steps=list(range(1,nsteps+1)) if config.get('trace_every_step',False) else sorted(set([1,*saved[1:]]))
            trace=[]
            for step in trace_steps:
                row={k:1. for k in controls.TRACE_FIELDS}
                row.update(evaluation_elapsed_s=(step-1)*timestep+(1 if wrong_trace_time else 0),
                    state_elapsed_s=step*timestep,mass_kg=m,equatorial_dimension_m=a,
                    polar_dimension_m=c,aspect_ratio=c/a,branch_criterion_met=0.,
                    deposition_density_kg_m3=300.,ventilation_vapor=-1. if negative_ventilation else 1.)
                trace.append(row)
            if not missing_trace:_write_csv(directory/'growth_trace.csv',trace)
            if config['growth_mode']==2:
                _write_csv(directory/'branch_diagnostic.csv',[dict(elapsed_s=0.,branch_threshold_ratio=.5,branch_criterion_met=0.)])
            if returncode==0:_write_csv(directory/'step_diagnostic.csv',[dict(maximum_step_mass_fraction=0.,steps=nsteps)])
            return SimpleNamespace(returncode=returncode)
        monkeypatch.setattr(controls.native.subprocess,'run',process)
        return build_dir,tmp_path/'run'
    return install


def _config():
    return controls.native.configuration(-14.4,.2,1,'fit',ice_density=919.)


@pytest.mark.parametrize('time_field,value',[
    ('dt_s',0.),('dt_s',-.1),('dt_s',float('nan')),('dt_s',float('inf')),
    ('duration_s',0.),('output_interval_s',float('nan')),('output_interval_s',.3),
])
def test_invalid_arbitrary_replay_times_fail_before_model_execution(tmp_path,time_field,value):
    config=_config();config[time_field]=value
    with pytest.raises(ValueError,match='Invalid time|does not divide timestep'):
        controls.run_configuration(tmp_path/'missing_build',tmp_path/'not_created',config)
    assert not (tmp_path/'not_created').exists()


def test_arbitrary_observed_age_preserves_full_axes_mass_and_short_step(fixture_run):
    obs=controls.native.lab.observed(-14.4,10.)
    config=_config();config.update(initial_mass_kg=obs['mass_kg'],
        initial_equatorial_m=obs['equatorial_dimension_m'],initial_polar_m=obs['polar_dimension_m'],
        dt_s=.01,duration_s=.01,output_interval_s=.01,trace_every_step=True)
    build_dir,output=fixture_run(config)
    result=controls.run_configuration(build_dir,output,config,metadata={'observed_age_min':10.})
    assert result['complete'] is True
    assert len(result['raw_snapshots'])==2 and len(result['growth_trace'])==1
    first=result['raw_snapshots'][0]
    assert 1.36e-3<first['equatorial_dimension_m']<1.37e-3
    assert first['equatorial_dimension_m']==obs['equatorial_dimension_m']
    assert first['polar_dimension_m']==obs['polar_dimension_m']
    assert first['mass_kg']==obs['mass_kg']
    assert first['ice_volume_m3']==pytest.approx(math.pi/6*obs['equatorial_dimension_m']**2*obs['polar_dimension_m'])
    assert first['aspect_ratio']==pytest.approx(obs['polar_dimension_m']/obs['equatorial_dimension_m'])
    assert result['growth_trace'][0]['evaluation_elapsed_s']==0.
    assert result['growth_trace'][0]['state_elapsed_s']==.01
    assert result['metadata']['observed_age_min']==10.


def test_wrong_full_axis_mapping_cannot_look_like_growth_failure(fixture_run):
    config=_config();config.update(duration_s=.2,output_interval_s=.2)
    build_dir,output=fixture_run(config,wrong_initial_axes=True)
    result=controls.run_configuration(build_dir,output,config)
    assert result['complete'] is False
    assert 'Initial mapping differs: equatorial_dimension_m' in result['output_validation_errors']


@pytest.mark.parametrize('problem,expected',[
    ({'missing_trace':True},'Missing growth_trace.csv'),
    ({'negative_ventilation':True},'Nonpositive source growth diagnostic'),
    ({'wrong_trace_time':True},'Incorrect source evaluation time'),
])
def test_standard_fit_analysis_preserves_generic_trace_errors(fixture_run,problem,expected):
    config=_config();build_dir,output=fixture_run(config,**problem)
    result=controls.run(build_dir,output,-14.4,.2,1,'fit')
    assert result['returncode']==0
    assert len(result['snapshots'])==28
    assert result['complete'] is False
    assert expected in result['output_validation_errors']
    saved=json.loads((output/'run.json').read_text())
    assert saved['complete'] is False and expected in saved['output_validation_errors']


def test_native_failure_retains_trace_and_partial_state_without_final_diagnostic(fixture_run):
    config=_config();build_dir,output=fixture_run(config,returncode=-8,truncate=True)
    result=controls.run(build_dir,output,-14.4,.2,1,'fit')
    assert result['complete'] is False and result['returncode']==-8
    assert [row['age_min'] for row in result['snapshots']]==[3.,4.]
    assert result['growth_trace']
    assert result['step_diagnostic'] is None
    assert 'growth_trace.csv' in result['outputs']


def test_existing_control_run_is_never_overwritten(fixture_run):
    config=_config();config.update(duration_s=.2,output_interval_s=.2)
    build_dir,output=fixture_run(config)
    controls.run_configuration(build_dir,output,config)
    before=(output/'run.json').read_bytes()
    with pytest.raises(FileExistsError):controls.run_configuration(build_dir,output,config)
    assert (output/'run.json').read_bytes()==before


def _parity_pair(tmp_path,returncode):
    directories=[tmp_path/'native',tmp_path/'control']
    for directory in directories:
        directory.mkdir()
        (directory/'run.json').write_text(json.dumps({'returncode':returncode,'complete':returncode==0}))
        for name in ('trajectory.csv','air.csv','branch_diagnostic.csv'):
            (directory/name).write_bytes((name+': fixture output').encode())
        if returncode==0:(directory/'step_diagnostic.csv').write_bytes(b'final step fixture')
    return directories


def test_failure_parity_compares_three_files_and_reports_missing_final_file(tmp_path):
    native,control=_parity_pair(tmp_path,-8)
    result=controls.parity(native,control)
    assert result['all_equal'] is True
    assert result['compared_files']==3
    assert result['files']['step_diagnostic.csv']['equal'] is None
    assert result['files']['step_diagnostic.csv']['matched_failure_absence'] is True


def test_success_parity_requires_four_files(tmp_path):
    native,control=_parity_pair(tmp_path,0)
    result=controls.parity(native,control)
    assert result['all_equal'] is True and result['compared_files']==4
    (control/'step_diagnostic.csv').unlink()
    result=controls.parity(native,control)
    assert result['all_equal'] is False
    assert result['files']['step_diagnostic.csv']['matched_failure_absence'] is False


def test_equal_partial_bytes_do_not_hide_different_failure_status(tmp_path):
    native,control=_parity_pair(tmp_path,-8)
    (control/'run.json').write_text(json.dumps({'returncode':-11,'complete':False}))
    result=controls.parity(native,control)
    assert result['status_equal'] is False and result['all_equal'] is False


def test_different_partial_bytes_are_not_unchanged_mode_parity(tmp_path):
    native,control=_parity_pair(tmp_path,-8)
    (control/'trajectory.csv').write_bytes(b'different trajectory')
    assert controls.parity(native,control)['all_equal'] is False
