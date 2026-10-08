"""Analytic checks of cloud optics and spherical volume transport."""
import math
import numpy as np
import pytest
from .microphysics import effective_radius, extinction


def test_gamma_effective_radius_preserves_mass_area_identity():
    rho=np.array([.4,.8]);q=np.array([1e-4,2e-4])
    re=effective_radius(q,None,rho,'qc')
    fields={s:q if s=='qc' else np.zeros(2) for s in ('qc','qr','qi','qs','qg')}
    fields.update({n:np.ones(2)*1e5 for n in ('nci','ncr','ncs','ncg')})
    beta,_=extinction(fields,rho)
    np.testing.assert_allclose(beta['qc']*re,3*rho*q/(2*997.))
    # Four times the liquid mass per volume gives the cube-root radius scaling.
    assert re[1]/re[0]==pytest.approx(4**(1/3))


def test_exponential_radius_follows_recorded_number_and_mass():
    q=1e-4;n=1e5;rho=.5
    expected=1.5/(math.pi*500*n/q)**(1/3)
    assert effective_radius(q,n,rho,'qi')==pytest.approx(expected)
    with pytest.raises(ValueError):
        effective_radius(q,-1,rho,'qi')


def test_hg_normalization_mean_and_inverse_sampling():
    pytest.importorskip('numba')
    from .volume import phase_hg,sample_hg
    x,w=np.polynomial.legendre.leggauss(400)
    for g in (0.,.5,.85,.95):
        p=np.array([phase_hg(v,g) for v in x])
        assert 2*np.pi*(p@w)==pytest.approx(1.,abs=1e-10)
        assert 2*np.pi*((x*p)@w)==pytest.approx(g,abs=1e-10)
        draws=np.array([sample_hg(v,g) for v in (np.arange(10000)+.5)/10000])
        assert draws.mean()==pytest.approx(g,abs=3e-6)


def test_radial_cloud_segments_match_beer_lambert_depth():
    pytest.importorskip('numba')
    from .volume import cloud_depth
    edges=np.array([1000.,1100.,1300.])
    beta=np.zeros((2,2,2));beta[0,:,0]=.002;beta[1,:,1]=.003
    value=cloud_depth(np.array([0.,0.,1000.01]),np.array([0.,0.,1.]),edges,-.1,.1,beta,50.)
    assert value==pytest.approx(.002*99.99+.003*200,abs=2e-5)
    assert cloud_depth(np.array([0.,100.,1000.01]),np.array([0.,0.,1.]),edges,-.1,.1,beta,50.)==0


def test_volume_absorbing_atmosphere_lambertian_ground():
    pytest.importorskip('numba')
    from .volume import estimate
    edges=np.array([1000.,1100.]);beta=np.zeros((1,1,2))
    values,cuts=estimate(np.array([0.,0.,1000.01]),np.array([[0.,0.,-1.]]),np.array([0.,0.,1.]),
                        edges,np.zeros((1,1)),np.array([[.01]]),np.ones((1,3)),edges,-.1,.2,beta,50.,
                        .2,.85,.8,photons=160,groups=4)
    # The solar beam traverses the full absorbing layer; viewing starts at ground level.
    assert values.mean(axis=1)[0,1]==pytest.approx(.2/math.pi*math.exp(-1),rel=.001)
    assert cuts.sum()==0


def test_direct_finite_sun_and_cloud_shadow():
    pytest.importorskip('numba')
    from .volume import direct_horizontal
    from illumination.sky.solved_transport import SUN_RADIUS_RAD
    edges=np.array([1000.,1100.,1300.])
    beta=np.zeros((2,2,2))
    origin=np.array([0.,0.,1000.01]);sun=np.array([0.,0.,1.])
    ext=np.zeros((2,1));xyz=np.ones((1,3))
    clear=direct_horizontal(origin,sun,edges,ext,xyz,edges,-.1,.1,beta,50.)
    assert clear[1]==pytest.approx(1-SUN_RADIUS_RAD**2/4,abs=1e-8)
    beta[:,:,0]=.005
    cloud=direct_horizontal(origin,sun,edges,ext,xyz,edges,-.1,.1,beta,50.)
    assert cloud[1]/clear[1]==pytest.approx(math.exp(-.005*299.99),rel=2e-5)
    beta*=0
    contact=direct_horizontal(origin,np.array([1.,0.,0.]),edges,ext,xyz,edges,-.1,.1,beta,50.,64,256)
    assert contact[1]==pytest.approx(2*SUN_RADIUS_RAD/(3*math.pi),rel=.001)


def test_earthlight_context_keeps_phase_and_visibility_separate():
    from .earthlight import context
    from shared.constants import EARTH_VISUAL_PHASE_NORMALISATION,EARTH_RADIUS,EARTH_MOON_DISTANCE,SYNODIC_MONTH_DAYS
    sky=[{'sun_deg':0.,'total_horizontal_lux':100.},{'sun_deg':90.,'total_horizontal_lux':100.}]
    earth=context(0.,0.,SYNODIC_MONTH_DAYS*24/4,sky)
    assert earth['earth_illuminated_fraction']==pytest.approx(1.)
    assert earth['clear_surface_proxy_lux']==pytest.approx(100*EARTH_VISUAL_PHASE_NORMALISATION*(EARTH_RADIUS/EARTH_MOON_DISTANCE)**2)
    far=context(0.,180.,0.,sky)
    assert far['earth_elevation_deg']==pytest.approx(-90.)
    assert far['clear_surface_proxy_lux'] is None
