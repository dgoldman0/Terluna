"""Independent forcing, trajectory and native interface checks for lunar probes."""
import json
import math

import numpy as np
import pytest

from atmosphere.snow import mcsnow_weather as weather
from atmosphere.snow import mcsnow_weather_analysis as analysis
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT, CM1_VAPOUR_GAS_CONSTANT


def test_cm1_mixing_ratio_recovers_known_partial_pressure_and_temperature():
    T,p,pv=259.,58000.,200.
    theta=T/(p/weather.history.t.CM1_P00)**(CM1_DRY_AIR_GAS_CONSTANT/weather.history.t.CM1_CP)
    ratio=(CM1_DRY_AIR_GAS_CONSTANT/CM1_VAPOUR_GAS_CONSTANT)*pv/(p-pv)
    raw=np.array([[theta,p,ratio,0,0,0,0,0,0]])
    assert weather.weather_air(raw)[:,0]==pytest.approx([T,p,pv],rel=1e-14)


def test_thermal_exit_fraction_retains_both_boundaries():
    before=np.array([260.,260.,260.]);after=np.array([280.,240.,265.])
    assert weather.thermal_fraction(before,after)==pytest.approx([13.15/20,16.85/20,1.])


def constant_weather():
    data=np.zeros((2,9,2,2,2))
    T,p,pv=260.,60000.,200.
    data[:,0]=T/(p/weather.history.t.CM1_P00)**(CM1_DRY_AIR_GAS_CONSTANT/weather.history.t.CM1_CP)
    data[:,1]=p
    data[:,2]=(CM1_DRY_AIR_GAS_CONSTANT/CM1_VAPOUR_GAS_CONSTANT)*pv/(p-pv)
    data[:,6]=2.
    record=dict(z_m=[20000.,21000.],dx_m=1000.,releases=[[0.,500.,500.,20200.]])
    return data,record


class ConstantParticle:
    """Independent exact transport control: no growth, one metre/second fall."""
    def seed(self,forcing,diameter):
        return np.asfortranarray([[1.],[.001],[1.],[260.],[0.]])

    def advance(self,state,forcing,duration,**options):
        n=state.shape[1];out=np.zeros((14,n),order='F')
        out[0:3]=.1;out[3]=1000.;out[4]=1.;out[7]=.1
        out[8]=np.broadcast_to(duration,n);out[10:12]=100.;out[12]=forcing[2]
        return out,np.zeros(n,np.int32)


def test_growing_probe_transport_reduces_to_exact_constant_settling():
    data,record=constant_weather()
    result=weather.integrate(ConstantParticle(),data,record,step_s=20.)
    assert result['status'].tolist()==[5]
    assert result['age_s'].tolist()==[200.]
    assert result['final_position_m'][0]==pytest.approx([900.,500.,20000.])
    assert np.array_equal(result['states'][0],result['initial_mass_kg'])
    assert result['mass_balance_residual_kg'].tolist()==[0.]


def test_different_particle_times_sample_the_same_linear_weather():
    data,record=constant_weather();data[1,6]=8.
    positions=np.array([[500.,500.,20500.],[500.,500.,20500.]])
    result=weather.variable_sample(data,positions,np.array(record['z_m']),record['dx_m'],np.array([150.,450.]),1)
    assert result[:,6]==pytest.approx([3.5,6.5])


@pytest.fixture(scope='module')
def library():
    if not (weather.CACHE/'build_moon/build.json').exists():
        pytest.skip('Pinned external McSnow build is not installed')
    return weather.Library()


def test_native_batch_matches_saved_lunar_growth(library):
    p=json.loads(weather.lunar.PRODUCT.read_text())
    case=next(r for r in p['fixed_cases'] if r['name']=='fixed_-14.4_moon_lunar_plate_base_dt0.0078125')
    expected=case['points'][-1]['state']
    air=np.array([[258.75],[47344.44821141621],[200.]],order='F')
    state=library.seed(air,4e-6)
    out,code=library.advance(state,air,1800.,saturated=True,max_step=.0078125,max_fraction=1.)
    assert code.tolist()==[0]
    assert state[0,0]==expected['mass_kg']
    assert out[0,0]==expected['diameter_m']
    assert out[1,0]==expected['equatorial_dimension_m']
    assert out[2,0]==expected['polar_dimension_m']


def test_condensation_ceiling_is_a_declared_separate_control(library):
    air=np.array([[260.],[60000.],[1000.]],order='F')
    initial=library.seed(air,4e-6)
    capped=initial.copy(order='F');saturated=initial.copy(order='F');raw=initial.copy(order='F')
    outc,_=library.advance(capped,air,60.,saturated=2)
    outs,_=library.advance(saturated,air,60.,saturated=1)
    outr,_=library.advance(raw,air,60.,saturated=0)
    assert np.array_equal(capped,saturated)
    assert outc[12,0]==outc[10,0]==outs[12,0]
    assert outr[12,0]==1000.
    assert raw[0,0]>capped[0,0]


def test_sublimation_retains_a_small_residual_and_signed_heat(library):
    air=np.array([[260.],[60000.],[100.]],order='F')
    state=library.seed(air,4e-6);mass=state[0,0]
    out,code=library.advance(state,air,600.)
    assert code.tolist()==[1]
    assert 0<state[0,0]<=library.ice_density*math.pi/6*(2e-7)**3
    assert state[0,0]-mass==pytest.approx(out[5,0],rel=1e-12)
    assert out[6,0]<0
    assert out[8,0]<600.


@pytest.mark.parametrize('bad',['time','forcing','shape','step'])
def test_native_interface_rejects_invalid_inputs(library,bad):
    air=np.array([[260.],[60000.],[200.]],order='F');state=library.seed(air,4e-6)
    options={};duration=1.
    if bad=='time':duration=math.nan
    elif bad=='forcing':air[1]=0
    elif bad=='shape':state[2]=-1
    else:options['max_step']=0.
    with pytest.raises(ValueError):library.advance(state,air,duration,**options)


@pytest.mark.parametrize('humidity',[0,1,2])
def test_release_humidity_record_matches_the_selected_forcing(library,humidity):
    data,record=constant_weather()
    result=weather.integrate(library,data,record,step_s=20.,saturated=humidity)
    initial=result['rows'][0]
    assert result['states'][4,0]==pytest.approx(initial[17]-1,rel=1e-14,abs=1e-14)


@pytest.fixture
def saved_arrays(library):
    folder=weather.CACHE/'runs/actual_small_step10'
    if not (folder/'raw.npz').exists():
        pytest.skip('Archived-weather integration has not been run')
    with np.load(folder/'raw.npz',allow_pickle=False) as f:
        arrays={k:f[k] for k in f.files}
    return arrays,json.loads((folder/'run.json').read_text()),json.loads((weather.CACHE/'weather/input.json').read_text())


def test_retained_paths_satisfy_geometry_mass_and_termination_checks(saved_arrays,library):
    arrays,meta,record=saved_arrays
    result=analysis.validate_arrays(arrays,meta,record,library.ice_density)
    assert result['mass_balance_relative']<1e-11


@pytest.mark.parametrize('damage',['geometry','mass_ledger','lost_censor','unfinished'])
def test_path_audit_rejects_inconsistent_physical_records(saved_arrays,library,damage):
    arrays,meta,record=saved_arrays
    if damage=='geometry':arrays['states'][1,0]*=1.01
    elif damage=='mass_ledger':arrays['deposited_mass_kg'][0]+=1e-4
    elif damage=='lost_censor':arrays['peak_diameter_m'][arrays['status']==2]=1e-6
    else:arrays['status'][0]=0
    with pytest.raises(ValueError):analysis.validate_arrays(arrays,meta,record,library.ice_density)


def test_run_file_identity_rejects_a_changed_raw_file(tmp_path):
    path=tmp_path/'raw.npz';path.write_bytes(b'original raw')
    meta={'raw_sha256':weather.native.sha(path)}
    path.write_bytes(b'altered raw')
    with pytest.raises(ValueError,match='Run file changed: raw.npz'):
        analysis.verify_files(tmp_path,meta)
