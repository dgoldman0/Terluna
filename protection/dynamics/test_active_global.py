import numpy as np
import pytest
from shared import constants as K

from .active_global import phase_warp, square_pair_separation, Layout, Traffic, RayIndex, sun_source
from .optical import length, unit


def test_registration_derivatives_and_smooth_return_join():
    theta=np.r_[np.linspace(-1.1,1.1,59),[-.95,-.65,.65,.95]]
    ratio=.998
    g,gp,gpp=phase_warp(theta,ratio)
    h=2e-5
    plus=phase_warp(theta+h,ratio)[0];minus=phase_warp(theta-h,ratio)[0]
    np.testing.assert_allclose((plus-minus)/(2*h),gp,atol=1e-9)
    # C2 joins have a third-derivative jump: a centred second difference
    # converges linearly there, quadratically in each smooth segment.
    np.testing.assert_allclose((plus-2*g+minus)/h**2,gpp,atol=3e-5)
    use=np.abs(theta)<=.65
    np.testing.assert_allclose(np.sin(g[use]),ratio*np.sin(theta[use]),atol=1e-14)
    assert np.all(g[np.abs(theta)>.95]==theta[np.abs(theta)>.95])


def test_square_test_detects_real_crossings_not_enclosing_sphere_overlap():
    q=np.array([[0.,0.,0.],[0,0,10],[0,0,0],[20,0,0]])
    a=np.tile([1.,0,0],(4,1));b=np.tile([0,1.,0],(4,1))
    b[2]=[0,0,1.]
    pairs=np.array([[0,1],[0,2],[0,3]])
    d=square_pair_separation(q,a,b,pairs,side=10.)
    assert d[0]>0 and d[1]<=0 and d[2]>0


@pytest.fixture(scope='module')
def environment():
    return make_environment()


def make_environment():
    """Domain-owned ephemeris fixture; independent of the study scenario."""
    from .ephemeris import Ephemeris, DEFAULT_KERNEL
    from .cycling import CyclingEnvironment
    pytest.importorskip('jplephem')
    if not DEFAULT_KERNEL.exists():pytest.skip('Pinned DE440s kernel is not installed')
    eph=Ephemeris()
    samples=eph.sample(np.arange(-3600, K.JULIAN_DAY+7200, 1800.))
    eph.close()
    return CyclingEnvironment(samples,.05,.25,'filter_only')


def test_global_shell_bound_covers_all_normal_clocks(environment):
    f=Traffic(environment,1.)
    bound=f.separation_certificate()
    assert bound['passes']
    rng=np.random.default_rng(938)
    rad=f.layout.first_radius
    tilt=np.deg2rad(f.layout.normal_limit_deg)
    az=rng.uniform(0,2*np.pi,10000)
    n=np.c_[np.full(len(az),np.cos(tilt)),np.sin(tilt)*np.cos(az),np.sin(tilt)*np.sin(az)]
    a=unit(np.cross(n,[0.,0,1.]));b=np.cross(n,a)
    roll=rng.uniform(0,2*np.pi,len(az));u=a*np.cos(roll)[:,None]+b*np.sin(roll)[:,None]
    v=-a*np.sin(roll)[:,None]+b*np.cos(roll)[:,None]
    corners=np.array([rad,0,0])+f.layout.side/2*(u+v)
    assert np.abs(length(corners)-rad).max()<bound['radial_half_extent_bound_m']
    rr=f.layout.first_radius/(f.layout.first_radius+f.layout.radial_step)
    assert phase_warp(np.linspace(-np.pi,np.pi,10001),rr)[1].min()>=bound['minimum_phase_derivative_bound']


def test_inverse_dynamics_acceleration_matches_actual_date_path_derivatives(environment):
    f=Traffic(environment,1.)
    j=np.array([0,101,981,1151,1152,1337,2201,2303])
    k=np.array([17,432,785,1203,990,871,23,689])
    t=12379.;h=.25
    state,acc,_=f.path(t,j,k)
    before=f.path(t-h,j,k)[0];after=f.path(t+h,j,k)[0]
    np.testing.assert_allclose((after[:,:3]-before[:,:3])/(2*h),state[:,3:],atol=2e-6,rtol=0)
    np.testing.assert_allclose((after[:,3:]-before[:,3:])/(2*h),acc,atol=2e-10,rtol=0)


def brute_hits(f,t,source,xy):
    j=np.concatenate([np.full(n,i) for i,n in enumerate(f.count)])
    k=np.concatenate([np.arange(n) for n in f.count])
    d=f.facets(t,j,k);frame=f.frames(t);sample=f.env.at(t)
    q=d['state'][:,:3];tau=np.maximum(q@frame[:,0],0)/K.SPEED_OF_LIGHT
    q=(q-(d['state'][:,3:]+sample['moon_v'])*tau[:,None])@frame
    n=d['normal']@frame;a=d['a']@frame;b=d['b']@frame;s=source@frame
    out=[]
    for x in xy:
        start=np.r_[0.,x];ray=unit(s-start)
        den=n@ray;distance=np.sum((q-start)*n,axis=1)/np.where(np.abs(den)>1e-12,den,1e-100)
        p=start+distance[:,None]*ray-q
        hit=(distance>0)&(distance<length(s-start))&(np.abs(np.sum(p*a,axis=1))<=f.layout.clear_side/2)&(np.abs(np.sum(p*b,axis=1))<=f.layout.clear_side/2)
        out.append(np.any(hit))
    return np.array(out)


def test_index_matches_exhaustive_finite_tile_ray_intersections(environment):
    # Sparse test population, actual hit interiors/edges plus unshielded rays.
    f=Traffic(environment,1.,Layout(planes=8,pitch=250000.))
    t=1379.;frame=f.frames(t);source=sun_source(environment.at(t),frame,[.8,-.5])
    j=np.tile(np.arange(8),4);k=np.rint(-f.omega[j]*t*f.count[j]/(2*np.pi)).astype(int)+np.repeat([-2,-1,1,2],8)
    d=f.facets(t,j,k);q=d['state'][:,:3]
    tau=np.maximum(q@frame[:,0],0)/K.SPEED_OF_LIGHT
    q=(q-(d['state'][:,3:]+environment.at(t)['moon_v'])*tau[:,None])@frame
    s=source@frame;xy=s[1:]+(q[:,1:]-s[1:])*(s[0]/(s[0]-q[:,0]))[:,None]
    xy=np.vstack([xy,xy+[0,4900],xy+[7000,7000],[[1e6,-3e6],[2e6,3e6]]])
    indexed=RayIndex(f,t,source).intercept(xy)
    brute=brute_hits(f,t,source,xy)
    assert np.any(brute) and np.any(~brute)
    np.testing.assert_array_equal(indexed,brute)


def test_shadow_force_bracket_contains_all_illumination_fractions():
    rng=np.random.default_rng(1924);b=rng.normal(size=(100,3));s=rng.normal(size=(100,3))
    lam=np.clip(np.sum(b*s,axis=-1)/np.sum(s*s,axis=-1),0,1)
    lo=length(b-lam[:,None]*s);hi=np.maximum(length(b),length(b-s))
    for light in np.linspace(0,1,100):
        x=length(b-light*s)
        assert np.all(x>=lo-1e-12) and np.all(x<=hi+1e-12)
