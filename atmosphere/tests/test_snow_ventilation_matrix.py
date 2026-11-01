"""Guard factorial coverage, evidence identities, and failed-run interpretation."""
import json
from pathlib import Path
import pytest
from atmosphere.snow import ventilation_benchmarks as matrix


def test_matrix_has_complete_independent_factorial_and_distinct_names():
    cases=list(matrix.matrix())
    assert len(cases)==240==len(set(cases))
    assert len({matrix.case_name(*case) for case in cases})==240
    assert set((case[0],case[1]) for case in cases)=={(0,0),(1,0),(0,1),(1,1),(2,0),(2,1)}
    for combo in matrix.COMBINATIONS:
        subset=[c for c in cases if c[:2]==combo]
        assert len(subset)==40
        assert {c[2] for c in subset}=={1,2}
        assert {c[3] for c in subset}=={-3.7,-5.3,-10.6,-14.4,-20.1}
        assert {c[5] for c in subset if c[4]=='fit'}=={1.,.2}
        assert {c[5] for c in subset if c[4]=='sphere'}=={.05,.0125}


def _trajectory(values,complete=True):
    return dict(complete=complete,snapshots=[dict(age_min=age,**{key:value for key in matrix.QUANTITIES}) for age,value in values])


def test_history_comparison_reports_signed_ratios_at_shared_complete_endpoint():
    fit=_trajectory([(3.,2.),(30.,6.)]);sphere=_trajectory([(0.,1.),(3.,1.),(30.,2.)])
    comparison=matrix.history_comparison(fit,sphere)
    assert comparison['common_age_min']==[3.,30.]
    assert comparison['both_complete'] is True
    assert all(x==3. for x in comparison['fitted_over_spherical_endpoint'].values())
    assert all(x==2. for x in comparison['maximum_relative_difference'].values())


@pytest.mark.parametrize('fit,sphere',[
    (_trajectory([(3.,1.),(30.,2.)],False),_trajectory([(0.,1.),(3.,1.),(30.,3.)])),
    (_trajectory([(3.,1.),(12.,2.)]),_trajectory([(0.,1.),(3.,1.),(30.,3.)])),
    (_trajectory([(3.,1.)]),_trajectory([(0.,1.)])),
])
def test_failed_or_unmatched_history_cannot_produce_an_endpoint_ratio(fit,sphere):
    assert matrix.history_comparison(fit,sphere)['fitted_over_spherical_endpoint'] is None


def test_failed_summary_retains_final_saved_state_but_no_accepted_endpoint(tmp_path,monkeypatch):
    monkeypatch.setattr(matrix,'ROOT',tmp_path)
    path=tmp_path/'failed';path.mkdir();(path/'run.json').write_text('{}')
    record=dict(snapshots=[{'age_min':30.}],configuration=dict(mean_mode=0,directional_mode=0,growth_mode=2,dt_s=.2),
        temperature_C=-14.4,initialization='fit',complete=False,returncode=-8,output_validation_errors=[],
        growth_trace=[{'evaluation_elapsed_s':1619.8}],branch_diagnostic=[],step_diagnostic=None,
        maximum_relative_ice_deposition_residual=0.)
    result=matrix.summarize(record,path)
    assert result['endpoint'] is None
    assert result['last_saved']=={'age_min':30.}
    assert result['last_source_trace']=={'evaluation_elapsed_s':1619.8}


@pytest.fixture
def verification_fixture(tmp_path,monkeypatch):
    identity={'producer':{'test_fixture':'1'}}
    monkeypatch.setattr(matrix.controls,'producer',lambda:identity['producer'])
    monkeypatch.setattr(matrix.controls,'verify_build',lambda path:identity)
    build=tmp_path/'build';(build/'source/src').mkdir(parents=True)
    (build/'build.json').write_text(json.dumps(identity))
    (build/'source/src/mo_atmo_types.f90').write_text('REAL(wp), PARAMETER :: rhoi = 919_wp\n')
    source=tmp_path/'cache'/matrix.native.scale.SOURCE_PDF;source.parent.mkdir(parents=True)
    source.write_bytes(b'empirical test source');monkeypatch.setattr(matrix.native.scale.lunar,'CACHE',tmp_path/'cache')
    monkeypatch.setattr(matrix.native.lab,'SOURCE_SHA',matrix.native.sha(source))
    config=matrix.native.configuration(-14.4,.2,1,'fit',ice_density=919.)
    config.update(mean_mode=0,directional_mode=0,trace_every_step=False)
    path=tmp_path/'run';path.mkdir();(path/'crystal.nml').write_text(matrix.native.namelist(config))
    (path/'trajectory.csv').write_bytes(b'fixture raw trajectory')
    record=dict(producer=identity['producer'],build_record_sha256=matrix.native.sha(build/'build.json'),
        source_pdf_sha256=matrix.native.lab.SOURCE_SHA,temperature_C=-14.4,initialization='fit',configuration=config,
        initial=matrix.native.initial(-14.4,'fit',919.),
        inputs={'crystal.nml':matrix.native.sha(path/'crystal.nml')},
        outputs={'trajectory.csv':matrix.native.sha(path/'trajectory.csv')})
    (path/'run.json').write_text(json.dumps(record))
    return path,build,record


@pytest.mark.parametrize('name',['crystal.nml','trajectory.csv'])
def test_verifier_rejects_changed_input_or_output_before_interpretation(verification_fixture,name):
    path,build,record=verification_fixture
    with (path/name).open('ab') as stream:stream.write(b'tampered')
    with pytest.raises(ValueError,match='Changed .* file'):
        matrix.verify_run(path,build)


def test_verifier_rejects_changed_build_identity(verification_fixture):
    path,build,record=verification_fixture
    (build/'build.json').write_text('{"different_build":true}')
    with pytest.raises(ValueError,match='Control build identity differs'):
        matrix.verify_run(path,build)


def test_input_hash_alone_cannot_accept_a_changed_empirical_initial_state(verification_fixture):
    path,build,record=verification_fixture
    record['configuration']['initial_equatorial_m']*=2
    (path/'crystal.nml').write_text(matrix.native.namelist(record['configuration']))
    record['inputs']['crystal.nml']=matrix.native.sha(path/'crystal.nml')
    (path/'run.json').write_text(json.dumps(record))
    with pytest.raises(ValueError,match='Physics inputs do not reproduce'):
        matrix.verify_run(path,build)


def test_matrix_execution_refuses_to_overwrite_a_completed_or_failed_run(tmp_path,monkeypatch):
    (tmp_path/'build').mkdir();(tmp_path/'build/build.json').write_text('{}')
    args=(0,0,1,-14.4,'fit',.2);path=tmp_path/'runs'/matrix.case_name(*args);path.mkdir(parents=True)
    (path/'run.json').write_text('{"previous_evidence":"preserved"}')
    monkeypatch.setattr(matrix,'matrix',lambda:iter([args]))
    monkeypatch.setattr(matrix.controls,'verify_build',lambda path:{})
    def forbidden(*args,**kwargs):raise AssertionError('Existing evidence must prevent execution')
    monkeypatch.setattr(matrix.controls,'run',forbidden)
    with pytest.raises(FileExistsError):matrix.run_matrix(tmp_path)
    assert json.loads((path/'run.json').read_text())=={'previous_evidence':'preserved'}
