"""Independent rate, generator, water, source and coupled-column controls."""
import math
import shutil
import numpy as np
import pytest
from scipy.linalg import expm
from atmosphere.snow import event_column as e, aggregation as ag, flakes as sf
from atmosphere.snow import collision_reference as reference


pytestmark=pytest.mark.skipif(shutil.which('g++') is None,reason='Event-column controls require a C++17 compiler (g++).')

@pytest.fixture(scope='module')
def moon(): return sf.moon_column()[0]


@pytest.mark.parametrize('temperature',[-15.,-10.,-1.])
def test_physical_properties_match_frozen_python(moon,temperature):
    diameter=np.geomspace(2e-6,10.,91);mass=e.GEOMETRY.mass(diameter)
    result=e.properties(moon,mass,temperature)
    state=moon.state(moon.height_of(sf.T0+temperature))
    v,re=ag.speeds(diameter,mass,moon.gravity,state,e.GEOMETRY.area_ratio)
    growth=ag.deposition(diameter,re,state,e.GEOMETRY.capacitance_ratio)
    np.testing.assert_allclose(result[:,0],diameter,rtol=3e-15)
    np.testing.assert_allclose(result[:,1],v,rtol=3e-13)
    np.testing.assert_allclose(result[:,2],growth,rtol=1e-13)


@pytest.mark.parametrize('particles',[2,8,25])
def test_constant_kernel_finite_population_generator(moon,particles):
    # Exact finite-volume pure-death chain, not the infinite-volume mean ODE.
    rate=2e-6;volume=1e-4;duration=10.
    Q=np.zeros((particles,particles))
    for n in range(2,particles+1):
        hazard=rate*n*(n-1)/(2*volume)
        Q[n-2,n-1]=hazard;Q[n-1,n-1]=-hazard
    probabilities=expm(duration*Q)[:,-1];values=np.arange(1,particles+1)
    expected=probabilities@values;variance=probabilities@(values-expected)**2
    sample=[]
    for seed in range(512):
        r=e.box(moon,np.full(particles,1e-9),volume,[duration],[particles*1e-9],
                seed=81000+seed,constant_kernel=rate)
        sample.append(r['moments'][0,0]*volume)
        assert abs(r['relative_mass_residual'][0])<4e-15
        assert r['coagulations']==particles-round(sample[-1])
    assert abs(np.mean(sample)-expected)<6*math.sqrt(variance/len(sample))


def test_gravitational_two_particle_hazard(moon):
    masses=e.GEOMETRY.mass(np.array([.0001,.0004]));kernel=reference.SnowKernel(moon)
    volume=.00001;rate=float(kernel(*masses))/volume;end=math.log(2)/rate
    hit=[]
    for seed in range(1024):
        r=e.box(moon,masses,volume,[end],[sum(masses)*1.00001],seed=19123+seed)
        hit.append(r['coagulations'])
        assert r['maximum_acceptance']<=1
        np.testing.assert_allclose(r['moments'][0,1],sum(masses)/volume,rtol=1e-15)
    assert abs(np.mean(hit)-.5)<6*math.sqrt(.25/len(hit))


def test_equal_mass_snow_cannot_collide(moon):
    masses=np.full(64,1e-9)
    r=e.box(moon,masses,.001,[0.,2000.],[1e-9,2e-9])
    np.testing.assert_array_equal(r['moments'][0],r['moments'][1])
    assert r['coagulations']==0


def test_closed_box_matches_independent_direct_sampler(moon):
    initial=np.repeat(reference.ordinary_initial(),8);volume=len(initial)/1e5
    query=e.GEOMETRY.mass(np.array([.00075,.0015,.003]))
    direct=[];candidate=[]
    for seed in range(32):
        candidate.append(e.box(moon,initial,volume,[21600.],query,seed=7100+seed)['moments'][0])
        r=reference.coagulate(initial,volume,[21600.],reference.SnowKernel(moon),seed=8100+seed)
        m=r['snapshots'][0]['masses_kg']
        direct.append([len(m)/volume,sum(m)/volume,sum(m*m)/volume])
    a,b=np.array(candidate),np.array(direct)
    error=np.sqrt(a.var(axis=0,ddof=1)/len(a)+b.var(axis=0,ddof=1)/len(b))
    assert abs(a[:,0].mean()-b[:,0].mean())<6*error[0]
    assert abs(a[:,2].mean()-b[:,2].mean())<6*error[2]
    np.testing.assert_allclose(a[:,1],b[:,1],rtol=2e-15)


@pytest.mark.parametrize('breakup',[False,True])
def test_coupled_shallow_column_water_number_and_supply(moon,breakup):
    r=e.column(moon,top_c=-2.,number_per_m3=1e5,snowfall_mm_h=3.,layers=4,
               area_m2=.0003,fill=.4,gather=.4,max_step_s=30.,breakup=breakup,seed=431)
    ledger=r['mass_ledger_kg_m2'];n=r['number_ledger_m2'];s=r['settings']
    assert abs(ledger['relative_residual'])<1e-12
    assert abs(n['relative_residual'])<1e-12
    assert ledger['deposited']<=3*s['total_time_s']/3600*(1+1e-13)
    assert ledger['fractional_fragment_remaining']<=r['diagnostics']['maximum_fractional_fragment_mass_bound_kg_m2']
    assert r['diagnostics']['final_time_s']==s['total_time_s']
    moments=np.array(r['arriving_cumulative_moments'])
    assert (np.diff(moments,axis=0)>=-1e-14).all()
    np.testing.assert_allclose(np.sum(r['arrival_time_blocks']['moments'],axis=0),moments[-1],rtol=3e-13)
    assert moments[-1,0]>0
    assert all(v>=0 for v in ledger.values() if v not in [ledger['residual'],ledger['relative_residual']])
    if not breakup:
        assert n['chipping_production']==ledger['fractional_fragment_remaining']==0


def test_source_only_keeps_every_particle_and_mass(moon):
    r=e.column(moon,top_c=-2.,layers=2,area_m2=.0001,fill=0,gather=.01,
               collisions=False,growth=False,settling=False,seed=44)
    m=r['mass_ledger_kg_m2'];n=r['number_ledger_m2']
    assert m['injected']==m['remaining']
    assert n['injected']==n['remaining']
    assert m['deposited']==m['freezing_exit']==n['coagulation_loss']==0
    assert r['summary']['mass_d50_m'] is None


def test_poisson_injection_mean_and_variance(moon):
    counts=[]
    for seed in range(256):
        r=e.column(moon,top_c=-.2,layers=1,area_m2=.0001,fill=0,gather=.01,
                   collisions=False,growth=False,settling=False,seed=1000+seed)
        counts.append(r['number_ledger_m2']['injected']*.0001)
    s=r['settings'];expected=s['injected_number_flux_m2_s']*s['area_m2']*s['total_time_s']
    assert abs(np.mean(counts)-expected)<6*math.sqrt(expected/len(counts))
    assert .65<np.var(counts,ddof=1)/expected<1.4


@pytest.mark.parametrize('kwargs',[{'area_m2':0},{'layers':0},{'max_step_s':-1},
                                  {'top_c':1},{'seed':-1},{'gather':0}])
def test_invalid_column_is_rejected(moon,kwargs):
    with pytest.raises(ValueError):e.column(moon,**kwargs)


def test_compiled_chipping_matches_frozen_mean_loss(moon):
    import ctypes as ct
    state=moon.state(moon.height_of(sf.T0-10))
    p,a=e.parameters(moon,[state])
    d=np.geomspace(1e-5,.6,81);other=np.full_like(d,.0017)
    left=np.ascontiguousarray(e.GEOMETRY.mass(d));right=np.ascontiguousarray(e.GEOMETRY.mass(other))
    out=np.zeros_like(left);err=ct.create_string_buffer(1024)
    assert e.library().snow_fragment_loss(p,a,left,right,len(left),out,err)==0
    vl=ag.speeds(d,left,moon.gravity,state,e.GEOMETRY.area_ratio)[0]
    vr=ag.speeds(other,right,moon.gravity,state,e.GEOMETRY.area_ratio)[0]
    small=np.minimum(left,right);energy=left*right/(left+right)*(vl-vr)**2
    expected=np.minimum(ag.fragments(small,energy,mass_to_size=e.GEOMETRY.size)*p[4],.5*small)
    np.testing.assert_allclose(out,expected,rtol=2e-11,atol=1e-25)
    assert (out>=0).all() and (out<=.5*small).all()


def test_upwind_settling_matches_independent_linear_generator(moon):
    # Marked Poisson injection and noninteracting transport admit exact compound
    # Poisson moments. This control may use particle counting variance; coupled
    # aggregation outputs instead need whole-realization uncertainty.
    top=-2.;layers=4;area=.00005;fill=1.;gather=1.
    world=moon;depth=world.height_of(sf.T0+top)-world.z_freeze;dz=depth/layers
    seed=.0002
    speed=float(ag.speeds(np.array([seed]),np.array([e.GEOMETRY.mass(seed)]),world.gravity,
                         world.state(world.height_of(sf.T0+top)),e.GEOMETRY.area_ratio)[0][0])
    crossing=depth/speed;flux=1e5*speed
    nodes,weights=np.polynomial.legendre.leggauss(48)
    d=np.exp((math.log(ag.START_M[0])+math.log(ag.START_M[1]))/2
             +nodes*math.log(ag.START_M[1]/ag.START_M[0])/2)
    masses=e.GEOMETRY.mass(d);expected=np.zeros(3)
    for m,diam,w in zip(masses,d,weights/2):
        matrix=np.zeros((layers+2,layers+2))
        for k in range(layers):
            state=world.state(world.z_freeze+(k+.5)*dz)
            v=float(ag.speeds(np.array([diam]),np.array([m]),world.gravity,state,e.GEOMETRY.area_ratio)[0][0])
            matrix[k,k]=-v/dz
            matrix[k-1 if k else layers,k]=v/dz
        matrix[layers-1,layers+1]=1.
        initial=np.zeros(layers+2);initial[-1]=1.
        before=(expm(fill*crossing*matrix)@initial)[layers]
        after=(expm((fill+gather)*crossing*matrix)@initial)[layers]
        expected+=area*flux*w*(after-before)*np.array([1.,m,m*m])
    r=e.column(world,top_c=top,layers=layers,area_m2=area,fill=fill,gather=gather,
               max_step_s=30.,collisions=False,growth=False,seed=3187)
    observed=np.array(r['arriving_cumulative_moments'])[-1,:3]*area
    assert abs(observed[0]-expected[0])<6*math.sqrt(expected[0])
    assert abs(observed[1]-expected[1])<6*math.sqrt(expected[2])
    assert abs(r['mass_ledger_kg_m2']['relative_residual'])<1e-12


def test_fixed_seed_is_reproducible(moon):
    kwargs=dict(top_c=-.5,layers=2,area_m2=.0001,fill=.2,gather=.2,seed=110)
    a,b=e.column(moon,**kwargs),e.column(moon,**kwargs)
    a['diagnostics'].pop('wall_seconds');b['diagnostics'].pop('wall_seconds')
    assert a==b
