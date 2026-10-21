"""Analytical controls for archived-flow exposure and censored tracer histories."""
import numpy as np
import pytest

from atmosphere.snow import moist_history as m


def constant_fields(frames=7):
    T,p,rh=263.15,70000.,1.05
    e=rh*np.exp(9.550426-5723.265/T+3.53068*np.log(T)-.00728332*T)
    eps=m.CM1_DRY_AIR_GAS_CONSTANT/m.CM1_VAPOUR_GAS_CONSTANT
    values=[T/(p/m.t.CM1_P00)**(m.CM1_DRY_AIR_GAS_CONSTANT/m.t.CM1_CP),
            p,eps*e/(p-e),2e-5,0.,2e-5,2.,1.,0.]
    return np.broadcast_to(np.asarray(values)[None,:,None,None,None],(frames,9,3,4,4)).copy()


def test_primitive_interpolation_and_periodic_scalar_grid():
    z=np.array([20000.,21000.,22000.]);frame=np.zeros((1,3,4,4))
    frame[0]=z[:,None,None]*.1+np.arange(4)[None,None,:]*3+np.arange(4)[None,:,None]*7
    points=np.array([[1500.,2500.,20500.],[1500.+4000.,2500.-4000.,20500.]])
    values=m.spatial(frame,points,z,1000.)[:,0]
    np.testing.assert_allclose(values,2050+3+14)
    with pytest.raises(ValueError,match='outside'):
        m.spatial(frame,np.array([[0.,0.,19999.]]),z,1000.)


def test_thermodynamic_conversion_recovers_prescribed_humidity():
    values=constant_fields()[0,:,0,0,0]
    state=m.air(values[None,:])
    assert state['T'][0]==pytest.approx(263.15,rel=1e-14)
    assert state['rh'][0]==pytest.approx(1.05,rel=1e-14)
    assert all(m.masks(state)[key][0] for key in ('winter','moist','near_saturated','liquid','snow'))
    assert not m.masks(state)['dry'][0]


def test_constant_flow_exposure_and_height_censoring_are_exact():
    data=constant_fields();z=np.array([20000.,21000.,22000.])
    seeds=np.array([[0.,1500.,1500.,21500.,v] for v in m.DRIFTS])
    result=m.integrate(data,z,1000.,seeds,step_s=60.)
    duration=np.minimum(3600.,1500./np.array([np.inf,.3,.6,1.]))
    duration[0]=3600.
    for key in ('tracked_h','winter_h','moist_h','liquid_h','first_moist_episode_h','longest_moist_episode_h'):
        np.testing.assert_allclose(result['per_tracer'][key],duration/3600,rtol=2e-13,atol=1e-13)
    np.testing.assert_allclose(result['per_tracer']['dry_h'],0,atol=0)
    assert result['censoring']==['time_window','time_window','lower_height','lower_height']
    assert result['first_moist_episode_right_censored']==[True]*4
    assert result['final_moist_episode_right_censored']==[True]*4
    np.testing.assert_allclose(result['per_tracer']['ice_supersaturation_h'],duration/3600*.05,rtol=1e-12)
    fine=m.integrate(data,z,1000.,seeds,step_s=30.)
    np.testing.assert_allclose(result['end_height_m'],fine['end_height_m'],rtol=1e-14)


def test_temperature_exit_is_separate_from_height_censoring():
    data=constant_fields();data[:,0]*=280./263.15
    seeds=np.array([[0.,1500.,1500.,21000.,0.]])
    result=m.integrate(data,np.array([20000.,21000.,22000.]),1000.,seeds)
    assert result['per_tracer']['tracked_h']==pytest.approx([1.])
    assert result['per_tracer']['winter_h']==[0.]
    assert result['censoring']==['time_window']
    assert result['first_moist_episode_right_censored']==[False]


def test_releases_are_reproducible_and_shared_between_drift_classes():
    data=constant_fields(211);z=np.array([20000.,21000.,22000.])
    a,manifest=m.releases(data,z,1000.,per_release=2)
    b,_=m.releases(data,z,1000.,per_release=2)
    np.testing.assert_array_equal(a,b)
    assert len(a)==6*2*len(m.DRIFTS)
    for group in a.reshape(-1,4,5):
        np.testing.assert_array_equal(group[:,:4],np.broadcast_to(group[0,:4],(4,4)))
    assert all(row['qualifying_cells']==48 for row in manifest)


def test_temporal_sampling_retains_linear_fields_and_rejects_misaligned_end():
    data=constant_fields();data[:,6]=np.arange(7)[:,None,None,None]
    pos=np.array([[1500.,1500.,21000.]]);z=np.array([20000.,21000.,22000.])
    for stride in (1,2,3):
        assert m.sample(data,pos,z,1000.,450.,stride)[0,6]==pytest.approx(.75)
    seeds=np.array([[0.,1500.,1500.,21000.,0.]])
    with pytest.raises(ValueError,match='share exact'):
        m.integrate(data[:-1],z,1000.,seeds,stride=3)


def test_boundary_fraction_marks_upper_and_lower_crossing_without_clipping_time():
    np.testing.assert_allclose(m.boundary_fraction(np.array([1.,1.,1.]),np.array([-1.,3.,1.5]),0.,2.),[.5,.5,1.])


def test_saved_history_identity_censoring_and_exposure_accounting():
    import json
    from shared.provenance import constants_changed
    product=json.loads(m.OUT.read_text())
    assert product['schema']==m.SCHEMA
    for path,digest in product['producer']['files'].items():
        assert m.sha(m.ROOT/path)==digest
    assert not constants_changed(product['producer']['constants'])
    seeds=np.asarray(product['configuration']['seeds'])
    assert len(seeds)==1416
    np.testing.assert_array_equal(seeds.reshape(-1,4,5)[:,:,:4],
                                  np.repeat(seeds.reshape(-1,4,5)[:,0:1,:4],4,axis=1))
    for case in product['cases'].values():
        a={k:np.asarray(v) for k,v in case['per_tracer'].items()}
        assert all(len(v)==len(seeds) and np.isfinite(v).all() for v in a.values())
        np.testing.assert_allclose(a['winter_h'],a['near_saturated_h']+a['dry_h'],atol=1e-11,rtol=1e-13)
        for lesser,greater in [('first_moist_episode_h','longest_moist_episode_h'),
                               ('longest_moist_episode_h','moist_h'),('moist_h','winter_h'),
                               ('winter_h','tracked_h'),('liquid_h','winter_h')]:
            assert np.all(a[lesser]<=a[greater]+1e-11)
        time_end=np.asarray(case['censoring'])=='time_window'
        np.testing.assert_allclose(a['tracked_h'][time_end],35-seeds[time_end,0]/3600,atol=1e-11)
        still_first=np.asarray(case['first_moist_episode_right_censored'])
        np.testing.assert_allclose(a['first_moist_episode_h'][still_first],a['tracked_h'][still_first],atol=1e-11)
