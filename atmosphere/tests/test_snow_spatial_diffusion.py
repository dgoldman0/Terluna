"""Analytic limits and conserved geometry for the independent spatial solver."""
import math
import numpy as np
import pytest
from atmosphere.snow.spatial_diffusion import Attachment, Domain, surface_flux


@pytest.mark.parametrize('law',[Attachment('constant',.3),Attachment('nucleation',.2,.007),Attachment('harrington',1,.03,10),Attachment('harrington',1,.03,1)])
def test_attachment_derivative_and_physical_bounds(law):
    s=np.logspace(-5,0,100)
    value,derivative=law.flux(s)
    eps=1e-5
    numeric=(law.flux(s*(1+eps))[0]-law.flux(s*(1-eps))[0])/(2*eps*s)
    assert np.all(value>=0) and np.all(value<=s)
    assert np.all(derivative>=0)
    assert np.allclose(derivative,numeric,rtol=2e-6,atol=1e-10)


def test_half_cell_constant_attachment_exact_series_resistance():
    s=np.array([0.,1e-8,.01,.1]); d=np.array([.01,.1,.3,1.]); lam=345.
    q,dq,surface=surface_flux(s,d,lam,Attachment('constant',.17))
    factor=lam*.17/(1+d*lam*.17)
    np.testing.assert_allclose(q,factor*s,rtol=1e-12)
    np.testing.assert_allclose(dq,factor,rtol=1e-12)
    np.testing.assert_allclose(surface+d*q,s,rtol=1e-12)


def test_nonlinear_surface_balance_has_no_negative_supply():
    center=np.logspace(-7,-.5,30)
    for law in [Attachment('nucleation',1,.007),Attachment('harrington',1,.173,10)]:
        q,_,s=surface_flux(center,.1,1000,law)
        np.testing.assert_allclose(s+.1*q,center,rtol=2e-9)
        assert np.all(s>=0) and np.all(s<=center)


@pytest.mark.parametrize('radius,depth',[(0.,0.),(.7,.4),(.9,1.5)])
def test_resolved_surface_areas_match_cylinder_and_cavity(radius,depth):
    phi=2.;d=Domain(aspect=phi,resolution=8,cavity_radius=radius,cavity_depth=depth)
    basal=2*np.sum(d.ice[d.ice[:,3]==1,1]);prism=2*np.sum(d.ice[d.ice[:,3]==0,1])
    assert basal==pytest.approx(2*math.pi)
    assert prism==pytest.approx(4*math.pi*(phi+radius*depth))


@pytest.mark.parametrize('lam',[.1,10.,1000.])
def test_refined_sphere_between_analytic_outer_boundary_bounds(lam):
    extent=80.;d=Domain(resolution=64,extent=extent,shape='sphere')
    ans=d.solve(.1,lam,[Attachment('constant')]*2,tolerance=1e-12)
    measured=ans['half_flux']/(2*math.pi*.1)
    assert 1/(1+1/lam-1/(math.sqrt(2)*extent))<=measured<=1/(1+1/lam-1/extent)
    assert ans['relative_mass_residual']<1e-8


def test_uniform_cylinder_corner_matches_original_numerical_table():
    values=[]
    for n in (32,64):
        d=Domain(resolution=n,extent=80)
        row=[]
        for face in (0,1):
            ans=d.imposed_flux(prism=1-face,basal=face)
            ids=np.where(d.ice[:,3]==face)[0]
            k=ids[np.argmin((d.ice[ids,4]-1)**2+(d.ice[ids,5]-1)**2)]
            row.append(ans['surface'][k])
            assert ans['relative_mass_residual']<1e-8
        values.append(row)
    extrap=(2**(2/3)*np.asarray(values[1])-values[0])/(2**(2/3)-1)
    # Nelson & Baker 1996 original Table B1; not the later analytic fit.
    np.testing.assert_allclose(extrap,[.690,.445],rtol=.015)


def test_heat_feedback_linear_attachment_matches_independent_algebra():
    d=Domain(aspect=2,resolution=12,cavity_radius=.7,cavity_depth=1.)
    laws=[Attachment('constant',.1),Attachment('constant',.4)]
    cold=d.solve(.05,1000,laws,tolerance=1e-12)
    h=.5;cap=d.dirichlet()['half_conductance']
    warm=d.solve(.05,1000,laws,heat_parameter=h,tolerance=1e-12)
    ratio=1/(1+h*cold['half_flux']/(.05*cap))
    assert warm['half_flux']/cold['half_flux']==pytest.approx(ratio,rel=1e-8)
    assert abs(warm['thermal_balance_residual'])<1e-12
    assert warm['relative_mass_residual']<1e-8


def test_zero_forcing_and_linear_humidity_response():
    d=Domain(resolution=8)
    laws=[Attachment('constant',.2)]*2
    assert d.solve(0,100,laws)['half_flux']==0
    low=d.solve(.01,100,laws);high=d.solve(.1,100,laws)
    assert high['half_flux']==pytest.approx(10*low['half_flux'],rel=1e-9)


@pytest.mark.parametrize('kwargs',[dict(aspect=0),dict(extent=2),dict(resolution=3),dict(cavity_radius=1,cavity_depth=.5),dict(cavity_radius=.5),dict(cavity_radius=.5,cavity_depth=2),dict(shape='sphere',aspect=2)])
def test_invalid_geometry_rejected(kwargs):
    with pytest.raises(ValueError):Domain(**kwargs)


def test_steep_attachment_root_cannot_cycle_between_bracket_ends():
    law=Attachment('harrington',1.,.004036860788000002,10.)
    center=np.array([.018773136020175496,.02,.001,0.])
    q,_,surface=surface_flux(center,1/48,236.53806137367076,law)
    np.testing.assert_allclose(surface+q/48,center,rtol=2e-12)
    from scipy.optimize import brentq
    expected=brentq(lambda s:s+236.53806137367076/48*law.flux(np.array([s]))[0][0]-center[0],0.,center[0],xtol=1e-16)
    assert surface[0]==pytest.approx(expected,rel=1e-12)
