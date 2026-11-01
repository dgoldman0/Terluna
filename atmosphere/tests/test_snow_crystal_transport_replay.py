"""Check conditional transport algebra, geometry, domains and source identities.

These controls verify the calculation, not autonomous crystal growth or lunar
transfer. Synthetic quadrature cases have analytical endpoint integrals.
"""
import csv
import hashlib
import json
import math
from pathlib import Path

import pytest
from atmosphere.rain import recovered
from atmosphere.snow import crystal_transport_replay as replay


def _state():
    return dict(temperature_K=260.,vapor_pressure_Pa=220.,ice_saturation_pressure_Pa=200.,
                conductivity_W_m_K=.025,diffusivity_m2_s=2e-5)


def test_no_latent_resistance_recovers_analytic_maxwell_vapor_flux(monkeypatch):
    state=_state();cap=1e-4;vent=1.7
    monkeypatch.setattr(replay.drops,'latent_heat',lambda temperature:-replay.flakes.LATENT_FUSION)
    expected=4*math.pi*cap*vent*state['diffusivity_m2_s']*(state['vapor_pressure_Pa']-state['ice_saturation_pressure_Pa'])/(replay.flakes.R_VAPOUR*state['temperature_K'])
    assert replay.mass_rate(cap,vent,state)==pytest.approx(expected,rel=2e-15)


def test_latent_heat_reduces_flux_and_common_ventilation_scales_both_resistances():
    state=_state();cap=1e-4
    pure_diffusion=4*math.pi*cap*state['diffusivity_m2_s']*(state['vapor_pressure_Pa']-state['ice_saturation_pressure_Pa'])/(replay.flakes.R_VAPOUR*state['temperature_K'])
    rate=replay.mass_rate(cap,1.,state)
    assert 0<rate<pure_diffusion
    assert replay.mass_rate(2*cap,1.,state)==pytest.approx(2*rate)
    assert replay.mass_rate(cap,3.,state)==pytest.approx(3*rate)


def test_ice_equilibrium_has_zero_growth_and_supersaturation_retains_sign():
    state=_state();state['vapor_pressure_Pa']=state['ice_saturation_pressure_Pa']
    assert replay.mass_rate(1e-4,2.,state)==0.
    state['vapor_pressure_Pa']=180.
    negative=replay.mass_rate(1e-4,2.,state)
    state['vapor_pressure_Pa']=220.
    assert negative<0 and replay.mass_rate(1e-4,2.,state)==pytest.approx(-negative)


def test_full_axis_volumes_equal_independent_regular_hexagon_and_ellipsoid_geometry():
    # Basal full span2 means regular hexagon side1, with six unit triangles.
    triangle_area=math.sqrt(3)/4
    expected_hex=6*triangle_area*3
    assert replay.volume(2.,3.,'hexagonal_prism')==pytest.approx(expected_hex)
    assert replay.volume(2.,3.,'spheroid')==pytest.approx(4*math.pi/3*1*1*1.5)
    assert replay.volume(2.,3.,'hexagonal_prism')/replay.volume(2.,3.,'spheroid')==pytest.approx(replay.HEX_TO_SPHEROID)
    assert replay.HEX_TO_SPHEROID==pytest.approx(1.240490014699032)


@pytest.mark.parametrize('a,c',[(0.,1.),(1.,-1.),(float('nan'),1.),(1.,float('inf'))])
def test_invalid_geometry_is_rejected(a,c):
    with pytest.raises(ValueError):replay.volume(a,c,'spheroid')


@pytest.mark.parametrize('aspect',[.009999,50.0001])
def test_hex_capacitance_rejects_extrapolation_past_sampled_thin_or_column_limits(aspect):
    with pytest.raises(ValueError,match='sampled hexagonal-prism aspect'):
        replay.hex_capacitance(1.,aspect)


@pytest.mark.parametrize('aspect',[.01,1.,50.])
def test_capacitance_is_linear_in_size_at_fixed_shape_with_inclusive_limits(aspect):
    a=2e-4;c=aspect*a
    assert replay.hex_capacitance(2*a,2*c)==pytest.approx(2*replay.hex_capacitance(a,c))


@pytest.mark.parametrize('ratio',[0.,.019999,.350001,float('nan')])
def test_branched_scenarios_reject_unselected_arm_ratios(ratio):
    with pytest.raises(ValueError):replay.branch_capacitance(1e-3,ratio,'dendrite')


def test_secondary_branches_raise_capacitance_without_exceeding_solid_plate_limit():
    for ratio in (.02,.35):
        star=replay.branch_capacitance(1e-3,ratio,'stellar')
        dendrite=replay.branch_capacitance(1e-3,ratio,'dendrite')
        assert 0<star<dendrite<.596*.001/2


@pytest.mark.parametrize('kind,lower,upper',[('column',.2,20.),('plate',1.,120.),('broad_branch',1.,120.)])
def test_jw_source_reynolds_intervals_are_inclusive_and_never_spliced(kind,lower,upper):
    assert replay.jw_ventilation(lower,.63,kind)>0
    assert replay.jw_ventilation(upper,.63,kind)>0
    for value in (lower*.999,upper*1.001):
        with pytest.raises(ValueError,match='published Reynolds interval'):replay.jw_ventilation(value,.63,kind)


@pytest.mark.parametrize('reynolds,schmidt',[(0.,.63),(float('nan'),.63),(1.,0.),(1.,float('inf'))])
def test_jw_rejects_nonphysical_transport_inputs(reynolds,schmidt):
    with pytest.raises(ValueError):replay.jw_ventilation(reynolds,schmidt,'column')


def test_jw_column_uses_basal_diameter_while_hpd_uses_maximum_dimension(monkeypatch):
    calls=[]
    monkeypatch.setattr(replay.lab,'observed',lambda temperature,age:dict(equatorial_dimension_m=1e-4,polar_dimension_m=8e-4))
    monkeypatch.setattr(replay.lab,'fall_speed',lambda temperature,age:.2)
    monkeypatch.setattr(replay,'capacitance_at',lambda *args:1e-4)
    monkeypatch.setattr(replay,'jw_ventilation',lambda re,sc,kind:calls.append((re,sc,kind)) or 2.)
    state={**_state(), 'viscosity_Pa_s':1e-5,'density_kg_m3':1.}
    jw=replay.diagnostics(-5.3,10.,'hexagonal_prism','JW',state)
    assert jw['reynolds_basal']==pytest.approx(2.)
    assert jw['reynolds_maximum_dimension']==pytest.approx(16.)
    assert calls==[(pytest.approx(2.),pytest.approx(.5),'column')]
    assert jw['jw_source_schmidt']==.63
    hpd=replay.diagnostics(-5.3,10.,'hexagonal_prism','HPD',state)
    assert hpd['ventilation']==pytest.approx(.86+.28*.5**(1/3)*math.sqrt(16.))
    assert hpd['jw_source_schmidt'] is None


def test_late_hexagonal_dendrite_has_no_endpoint_outside_thickness_domain():
    result=replay.replay(-14.4,'hexagonal_prism','JW',1e-9)
    assert result['complete'] is False
    assert 'aspect range' in result['domain_error']
    assert result['predicted_end_mass_kg'] is None
    assert result['mass_ratio_to_fit'] is None


def test_piecewise_quadrature_splits_all_axis_mass_speed_and_habit_joins(monkeypatch):
    monkeypatch.setattr(replay.lab,'FITS',{-1.:(12.,'test',[(4.,1.,1.),(12.,1.,1.)],[(7.,1.,1.),(12.,1.,1.)],[(12.,1.,1.)])})
    monkeypatch.setattr(replay.lab,'FALL_SPEED',{-1.:[(6.,1.,1.,0.),(12.,1.,1.,0.)]})
    monkeypatch.setattr(replay.lab,'observed',lambda temperature,age:{'mass_kg':1e-6})
    monkeypatch.setattr(replay,'air_state',lambda temperature:{})
    def diagnostic(temperature,age,*args):
        rate=1e-9 if age<4 else 2e-9 if age<7 else 3e-9 if age<10 else 4e-9
        return {'predicted_mass_rate_kg_s':rate}
    monkeypatch.setattr(replay,'diagnostics',diagnostic)
    original_quad=replay.quad;intervals=[]
    def traced_quad(function,lo,hi,**kwargs):
        intervals.append((lo,hi));return original_quad(function,lo,hi,**kwargs)
    monkeypatch.setattr(replay,'quad',traced_quad)
    result=replay.replay(-1.,'hexagonal_prism','diffusion',1e-9)
    assert intervals==[(3.,4.),(4.,6.),(6.,7.),(7.,10.),(10.,12.)]
    assert result['complete'] is True
    assert result['predicted_end_mass_kg']==pytest.approx(1e-6+60*(1e-9+3*2e-9+3*3e-9+2*4e-9),rel=1e-14)


def test_fitted_derivative_is_rejected_at_an_actual_piecewise_join():
    with pytest.raises(ValueError,match='undefined at a piecewise join'):replay.fitted_exponents(-14.4,4.)
    assert replay.fitted_exponents(-14.4,10.)==[2.37,1.02,.52]


@pytest.fixture
def source_fixture(tmp_path,monkeypatch):
    source=tmp_path/'formula_sources';source.mkdir()
    entries=[]
    for name in replay.REQUIRED_SOURCES:
        data=('identified '+name).encode();(source/name).write_bytes(data)
        entries.append(dict(name=name,status='failed',error='earlier retrieval failed'))
        entries.append(dict(name=name,status='retrieved',url='https://example.test/'+name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    record=tmp_path/'sources.json';record.write_text(json.dumps(entries))
    empirical=tmp_path/'takahashi.pdf';empirical.write_bytes(b'identified empirical source')
    monkeypatch.setattr(replay.lab,'SOURCE_SHA',replay.sha(empirical))
    return source,record,empirical,entries


def test_source_verification_requires_actual_bytes_and_preserves_retrieval_history(source_fixture):
    source,record,empirical,entries=source_fixture
    result=replay.verify_sources(source,record,empirical)
    assert len(result['formula_sources'])==4
    assert result['source_register_sha256']==replay.sha(record)
    assert json.loads(record.read_text())==entries
    assert result['empirical_source']['sha256']==replay.sha(empirical)


@pytest.mark.parametrize('failure',['same_size_tamper','missing','no_success_identity','empirical_tamper'])
def test_failed_or_changed_source_cannot_silently_produce_results(source_fixture,failure):
    source,record,empirical,entries=source_fixture
    target=source/replay.REQUIRED_SOURCES[0]
    if failure=='same_size_tamper':target.write_bytes(b'X'+target.read_bytes()[1:])
    elif failure=='missing':target.unlink()
    elif failure=='no_success_identity':record.write_text(json.dumps([e for e in entries if e.get('status')!='retrieved' or e.get('name')!=target.name]))
    else:empirical.write_bytes(b'wrong measurement source')
    with pytest.raises((ValueError,FileNotFoundError)):replay.verify_sources(source,record,empirical)


def test_calculation_checks_sources_before_starting_any_comparison(monkeypatch):
    def fail():raise ValueError('source identity unavailable')
    monkeypatch.setattr(replay,'verify_sources',fail)
    def forbidden(*args,**kwargs):raise AssertionError('No result should be computed first')
    monkeypatch.setattr(replay,'replay',forbidden)
    with pytest.raises(ValueError,match='source identity unavailable'):replay.calculate()


def test_product_identity_includes_thermodynamics_and_shared_constant_inputs(monkeypatch):
    monkeypatch.setattr(replay,'verify_sources',lambda:{'fixture':'verified'})
    monkeypatch.setattr(replay.lab,'FITS',{})
    monkeypatch.setattr(replay,'air_transport_sensitivity',lambda cases:{'fixture':'separate verified native inputs'})
    result=replay.calculate()
    path='atmosphere/radiative_convective/thermodynamics.py'
    assert result['producer'][path]==replay.sha(replay.ROOT/path)
    assert 'shared/constants.py' not in result['producer']
    assert 'shared/constants.json' not in result['producer']
    assert result['imported_constants']['GAS_CONSTANT']==replay.flakes.GAS_CONSTANT
    assert result['source_verification']=={'fixture':'verified'}


def test_fixed_source_schmidt_changes_only_jw_ventilation_and_not_physical_air(monkeypatch):
    # The replay is a 22 October diagnosis built with the rain solver of that time, whose Chapman-Enskog diffusivity
    # put the air's Schmidt number near 0.73 at -5.3 C; with Massman's measured fit it is near 0.61, below the fixed
    # 0.63. The test checks the product's own physics, the solver it was built with (shared/baselines.py).
    monkeypatch.setattr(replay,'drops',recovered.baseline())
    state=replay.air_state(-5.3)
    nominal=replay.diagnostics(-5.3,10.,'hexagonal_prism','JW',state)
    fixed=replay.diagnostics(-5.3,10.,'hexagonal_prism','JW',state,jw_schmidt=.63)
    assert fixed['ventilation_schmidt']==.63
    assert nominal['ventilation_schmidt']==nominal['schmidt']
    for key in ('schmidt','capacitance_m','reynolds_basal','reynolds_maximum_dimension','observed_fit_speed_m_s'):
        assert fixed[key]==nominal[key]
    assert fixed['predicted_mass_rate_kg_s']/nominal['predicted_mass_rate_kg_s']==pytest.approx(fixed['ventilation']/nominal['ventilation'])
    assert fixed['ventilation']<nominal['ventilation']


def test_schmidt_override_cannot_be_silently_applied_to_another_control():
    with pytest.raises(ValueError,match='only to JW'):
        replay.diagnostics(-5.3,10.,'hexagonal_prism','diffusion',replay.air_state(-5.3),jw_schmidt=.63)


def test_fixed_schmidt_replay_keeps_initial_mass_geometry_and_air(monkeypatch):
    monkeypatch.setattr(replay,'drops',recovered.baseline())   # the solver the replay was built with (above)
    nominal=replay.replay(-5.3,'hexagonal_prism','JW',1e-9)
    fixed=replay.replay(-5.3,'hexagonal_prism','JW',1e-9,jw_schmidt=.63)
    assert nominal['complete'] and fixed['complete']
    assert fixed['air']==nominal['air']
    assert fixed['initial_age_min']==nominal['initial_age_min']==3.
    assert fixed['observed_end_fit_mass_kg']==nominal['observed_end_fit_mass_kg']
    assert fixed['predicted_end_mass_kg']<nominal['predicted_end_mass_kg']
    assert fixed['jw_schmidt_override']==.63 and nominal['jw_schmidt_override'] is None


def test_native_coefficient_control_changes_only_diffusivity_and_conductivity(monkeypatch):
    state=replay.air_state(-5.3)
    source=dict(temperature_C=-5.3,mixture_air=state,native_diffusivity_m2_s=state['diffusivity_m2_s']*1.15,
                native_conductivity_W_m_K=state['conductivity_W_m_K']*1.02,native_density_kg_m3=999.,native_viscosity_Pa_s=999.)
    monkeypatch.setattr(replay,'native_transport_inputs',lambda:[source])
    baseline=replay.replay(-5.3,'hexagonal_prism','JW',1e-9)
    result=replay.air_transport_sensitivity([baseline]);case=result['cases'][0]
    changed=case['native_coefficients']['air']
    assert {k for k in state if state[k]!=changed[k]}=={'diffusivity_m2_s','conductivity_W_m_K'}
    assert changed['diffusivity_m2_s']==source['native_diffusivity_m2_s']
    assert changed['conductivity_W_m_K']==source['native_conductivity_W_m_K']
    assert baseline['air']==state
    assert case['both_complete'] is True and case['native_over_mixture_endpoint_mass']>1.
    initial=replay.lab.observed(-5.3,3.)['mass_kg']
    assert case['native_over_mixture_deposited_mass']==pytest.approx(
        (case['native_coefficients']['predicted_end_mass_kg']-initial)/(baseline['predicted_end_mass_kg']-initial))
    for before,after in zip(baseline['diagnostics'],case['native_coefficients']['diagnostics']):
        assert before['capacitance_m']==after['capacitance_m']
        assert before['reynolds_basal']==after['reynolds_basal']
        assert after['schmidt']==pytest.approx(before['schmidt']/1.15)


@pytest.fixture
def native_input_fixture(tmp_path,monkeypatch):
    from atmosphere.snow import mcsnow_ventilation_controls as controls
    producer={'fixture':'native source identity'};build=tmp_path/'build';build.mkdir()
    (build/'build.json').write_text(json.dumps({'producer':producer}))
    monkeypatch.setattr(controls,'producer',lambda:producer)
    monkeypatch.setattr(controls,'verify_build',lambda path:{'producer':producer})
    monkeypatch.setattr(replay,'TRANSPORT_TEMPERATURES',(-5.3,))
    path=tmp_path/'runs/m0_d0_igf2_sphere_T5p3_dt0p05';path.mkdir(parents=True)
    config={'mean_mode':0,'directional_mode':0}
    (path/'crystal.nml').write_text(controls.native.namelist(config))
    state=replay.air_state(-5.3)
    air={k:state[k] for k in ('temperature_K','pressure_Pa','vapor_pressure_Pa','density_kg_m3')}
    air['dynamic_viscosity_Pa_s']=state['viscosity_Pa_s']
    trace=[dict(diffusivity_raw_m2_s=state['diffusivity_m2_s']*1.15,conductivity_raw_W_m_K=state['conductivity_W_m_K']*1.02)]
    for filename,rows in [('air.csv',[air]),('growth_trace.csv',trace)]:
        with (path/filename).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    record=dict(producer=producer,build_record_sha256=replay.sha(build/'build.json'),complete=True,
        configuration=config,temperature_C=-5.3,inputs={'crystal.nml':replay.sha(path/'crystal.nml')},
        outputs={n:replay.sha(path/n) for n in ('air.csv','growth_trace.csv')},air=air,growth_trace=trace)
    (path/'run.json').write_text(json.dumps(record))
    return tmp_path,path,record


def test_native_transport_input_retains_exact_source_coefficients_and_hashes(native_input_fixture):
    cache,path,record=native_input_fixture
    source=replay.native_transport_inputs(cache)[0]
    assert source['native_diffusivity_m2_s']==record['growth_trace'][0]['diffusivity_raw_m2_s']
    assert source['run_sha256']==replay.sha(path/'run.json')
    assert source['inputs']==record['inputs'] and source['outputs']==record['outputs']
    assert source['producer']==record['producer']


@pytest.mark.parametrize('tamper',['raw_bytes','saved_numeric_value','intervention_mode','wrong_forcing'])
def test_native_transport_control_rejects_unidentified_or_mismatched_inputs(native_input_fixture,tamper):
    cache,path,record=native_input_fixture
    if tamper=='raw_bytes':
        with (path/'growth_trace.csv').open('ab') as stream:stream.write(b'changed')
    elif tamper=='saved_numeric_value':record['growth_trace'][0]['diffusivity_raw_m2_s']*=2
    elif tamper=='intervention_mode':record['configuration']['mean_mode']=1
    else:record['air']['temperature_K']+=1
    (path/'run.json').write_text(json.dumps(record))
    with pytest.raises(ValueError):replay.native_transport_inputs(cache)
