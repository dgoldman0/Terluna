"""The observing driver must preserve the frozen stochastic calculation."""
import copy
import shutil
import numpy as np
import pytest
from atmosphere.snow import column_diagnosis as d, event_column as e, flakes as sf

pytestmark=pytest.mark.skipif(shutil.which('g++') is None,reason='Diagnosis requires a C++17 compiler.')

@pytest.fixture(scope='module')
def moon():return sf.moon_column()[0]


def scientific(result):
    result=copy.deepcopy(result)
    result.pop('schema');result.pop('process_diagnostics',None)
    result['diagnostics'].pop('build');result['diagnostics'].pop('wall_seconds')
    return result


@pytest.mark.parametrize('changes',[{}, {'growth':False},{'collisions':False},
                                   {'breakup':False},{'settling_scheme':'ballistic'},
                                   {'growth':False,'collisions':False,'settling':False}])
def test_observation_preserves_every_original_field(moon,changes):
    kwargs=dict(top_c=-2.,layers=4,area_m2=.0001,fill=.3,gather=.3,max_step_s=30.,seed=210731)
    kwargs.update(changes)
    original=e.column(moon,**kwargs);observed=d.column(moon,**kwargs)
    assert scientific(original)==scientific(observed)


@pytest.mark.parametrize('changes',[{}, {'growth':False},{'collisions':False},{'breakup':False}])
def test_diagnostics_close_independent_budgets_and_event_counts(moon,changes):
    result=d.column(moon,top_c=-2.,layers=4,area_m2=.0001,fill=.3,gather=.3,seed=812,**changes)
    a=np.array(result['process_diagnostics']['blocks']);s=result['settings'];area=s['area_m2']
    assert np.isfinite(a).all() and (a>=0).all()
    total=s['total_time_s'];gather=s['gathering_time_s'];times=np.r_[total-gather,np.full(4,gather/4)]
    np.testing.assert_allclose(a[:,:,0],np.repeat(times[:,None],s['layers'],axis=1),rtol=2e-14)
    assert (a[:,:,1:4]<=a[:,:,:1]*(1+1e-14)).all()
    np.testing.assert_allclose(a[:,:,5].sum()/area,result['mass_ledger_kg_m2']['deposited'],rtol=2e-12,atol=1e-12)
    supply=(s['snowfall_mm_h']/3600/s['depth_m']*s['layer_volume_m3'] if s['growth'] else 0)
    np.testing.assert_allclose(a[:,:,5]+a[:,:,6]+a[:,:,7],supply*a[:,:,0],rtol=2e-13,atol=1e-18)
    np.testing.assert_allclose(a[:,:,6],supply*a[:,:,1],rtol=2e-13,atol=1e-18)
    assert a[:,:,12].sum()/area==result['number_ledger_m2']['coagulation_loss']
    assert a[:,:,13].sum()==result['diagnostics']['chipping_events']
    assert a[:,:,14].sum()==result['diagnostics']['emitted_fragments']
    assert (a[:,:,10]<=a[:,:,8]*(1+1e-14)).all()
    assert (a[:,:,11]<=a[:,:,0]*(1+1e-14)).all()
    # Integrated N^2 must dominate integrated N and its time-normalized square.
    assert (a[:,:,9]+1e-10>=a[:,:,8]).all()
    assert (a[:,:,9]*a[:,:,0]+1e-6>=a[:,:,8]**2).all()
    if not s['growth']:assert (a[:,:,3:8]==0).all()
    if not s['collisions']:assert (a[:,:,12:]==0).all()
    if not s['breakup']:assert (a[:,:,13:]==0).all()


def test_empty_column_accounts_for_all_unused_supply(moon):
    # Practically zero injection gives a reproducible, explicitly empty realization.
    r=d.column(moon,top_c=-.2,layers=3,area_m2=1e-6,number_per_m3=1e-20,fill=.1,gather=.1,seed=314)
    a=np.array(r['process_diagnostics']['blocks'])
    assert r['number_ledger_m2']['injected']==0
    np.testing.assert_array_equal(a[:,:,0],a[:,:,1])
    assert (a[:,:,2:6]==0).all() and (a[:,:,7:]==0).all()
    assert r['mass_ledger_kg_m2']['relative_residual']==0
    assert (a[:,:,6]>0).all()


def test_source_only_occupancy_matches_poisson_expectation(moon):
    # One layer with no transport or microphysics is a homogeneous Poisson source.
    # Its mean integrated occupation is lambda*T^2/2 under the symmetric source split.
    sample=[];empty=[]
    for seed in range(128):
        r=d.column(moon,top_c=-.1,layers=1,area_m2=1e-7,fill=0,gather=.03,
                   max_step_s=2.,collisions=False,growth=False,settling=False,seed=48100+seed)
        a=np.array(r['process_diagnostics']['blocks'])
        sample.append(a[:,:,8].sum());empty.append(a[:,:,1].sum())
    s=r['settings'];rate=s['injected_number_flux_m2_s']*s['area_m2'];T=s['total_time_s']
    mean=rate*T*T/2;variance=rate*T**3/3
    assert abs(np.mean(sample)-mean)<6*np.sqrt(variance/len(sample))
    # Exact step-sampled empty expectation at source midpoint times, across block boundaries.
    bounds=np.linspace(0,T,5);expected=0.
    for lo,hi in zip(bounds[:-1],bounds[1:]):
        t=lo
        while t<hi:
            h=min(s['max_step_s'],hi-t);expected+=h*np.exp(-rate*(t+h/2));t+=h
    assert abs(np.mean(empty)-expected)<6*np.std(empty,ddof=1)/np.sqrt(len(empty))
