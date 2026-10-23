import math
import numpy as np
import pytest
from atmosphere.snow.surface_transfer import capture_average, transfer_rates, prism_spiral_rate
from atmosphere.snow.spatial_diffusion import Domain, Attachment
from atmosphere.snow.facet_growth import ICE_DENSITY

AIR=dict(ice_saturation_pressure_Pa=400.,vapor_gas_constant_J_kg_K=460.,
    vapor_molecular_speed_m_s=560.,latent_sublimation_J_kg=2.8e6,
    diffusivity_m2_s=2e-5,conductivity_W_m_K=.024)

@pytest.mark.parametrize('xs',[14e-9,.5e-6,5e-6,1e-3])
def test_capture_integrates_to_the_same_finite_supply_on_unequal_grids(xs):
    height=100e-6
    expected=xs*(-math.expm1(-height/xs))
    for power in (1,2,4):
        x=height*np.linspace(0,1,47)**power
        p=capture_average(x[:-1],x[1:],xs)
        assert np.dot(np.diff(x),p)==pytest.approx(expected,rel=2e-14)
        assert np.all((p>=0)&(p<=1))


def test_transfer_redistributes_only_the_solved_gas_supply():
    d=Domain(3.,resolution=12,extent=40.,cavity_radius=.9,cavity_depth=1.5,mesh='edge')
    laws=[Attachment('nucleation',1.,.006),Attachment('nucleation',1.,.002)]
    off=transfer_rates(d,AIR,40e-6,268.,.05,laws,0.)
    on=transfer_rates(d,AIR,40e-6,268.,.05,laws,5e-6)
    blocked=transfer_rates(d,AIR,40e-6,268.,.05,laws,5e-6,transmission=0.)
    assert blocked['mass_rate_kg_s']==pytest.approx(off['mass_rate_kg_s'],rel=1e-11)
    assert on['captured_fraction_of_total']>0
    assert on['mass_rate_kg_s']>off['mass_rate_kg_s']
    assert on['tip_mean_normal_rate_m_s']>off['tip_mean_normal_rate_m_s']
    assert on['outer_prism_mean_normal_rate_m_s']<off['outer_prism_mean_normal_rate_m_s']
    for row in (off,on,blocked):
        assert row['transfer_mass_residual']<1e-10
        assert row['transfer_recipient_residual']<1e-10
        assert row['gas_water_residual']<1e-8
        assert 0<=row['captured_fraction_of_total']<=row['tip_fraction_of_deposition']<=1


def test_zero_capture_recovers_the_independent_stationary_flux():
    radius=40e-6;T=268.;si=.05
    d=Domain(2.,resolution=12,extent=40.)
    laws=[Attachment('constant',.2),Attachment('constant',1.)]
    expected=d.solve(si,radius*(AIR['vapor_molecular_speed_m_s']/4)/AIR['diffusivity_m2_s'],laws,tolerance=1e-11)
    got=transfer_rates(d,AIR,radius,T,si,laws,0.,heating=False)
    rhoeq=AIR['ice_saturation_pressure_Pa']/(AIR['vapor_gas_constant_J_kg_K']*T)
    assert got['mass_rate_kg_s']==pytest.approx(2*rhoeq*AIR['diffusivity_m2_s']*radius*expected['half_flux'],rel=1e-8)


def test_spiral_normal_growth_is_quadratic_at_small_supersaturation():
    a=prism_spiral_rate(268.,1e-6,2e-10)
    b=prism_spiral_rate(268.,2e-6,2e-10)
    assert b['normal_rate_m_s']/a['normal_rate_m_s']==pytest.approx(4.,rel=1e-6)
    assert a['normal_rate_m_s']==pytest.approx(a['lateral_step_velocity_m_s']*.39e-9/a['spacing_m'])
    assert 1e-5<a['step_coefficient_m_s']<1e-4


def test_capture_approaches_its_geometric_collection_area():
    radius=40e-6;xs=.5e-6;phi=3.
    d=Domain(phi,resolution=12,extent=40.)
    got=transfer_rates(d,AIR,radius,268.,.05,[Attachment('constant',.2)]*2,xs)
    # Two reflected outer prism faces together have circumference 2*pi*a,
    # each with an integrated collection length xs*(1-exp(-c/xs)).
    expected=4*math.pi*radius*xs*(-math.expm1(-phi*radius/xs))
    assert got['capture_area_m2']==pytest.approx(expected,rel=1e-12)


@pytest.mark.parametrize('T',[240.,273.15,float('nan')])
def test_spiral_law_rejects_extrapolation_outside_its_temperature_interval(T):
    with pytest.raises(ValueError):prism_spiral_rate(T,.01,2e-10)


def test_curvature_uses_the_same_recipient_equilibrium_for_direct_and_captured_growth():
    from atmosphere.snow.surface_transfer import CaptureBoundary
    law=CaptureBoundary(Attachment('nucleation',1.,.006),.8,False,recipient_offset=.002)
    s=np.array([.001,.003,.01,.03])
    p,dp=law.flux(s)
    h=1e-8
    numerical=(law.flux(s+h)[0]-law.flux(s-h)[0])/(2*h)
    np.testing.assert_allclose(dp,numerical,rtol=1e-6,atol=1e-9)
    captured=.8*np.maximum(s-.002,0)
    assert np.all(p>=captured)
    assert np.all(p<=s)
    d=Domain(3.,resolution=12,extent=40.,cavity_radius=.9,cavity_depth=1.5,mesh='edge')
    # Deliberately high curvature puts the recipient above even ambient vapour
    # chemical potential. Flat side faces can grow; the tip receives no ice.
    result=transfer_rates(d,AIR,40e-6,268.,.05,[Attachment('constant',.2)]*2,
        5e-6,capillary_length=200e-9)
    assert result['mass_rate_kg_s']>0
    assert result['tip_mean_normal_rate_m_s']==0
    assert result['captured_fraction_of_total']==0
    assert result['transfer_mass_residual']<1e-10


def test_graded_corner_mesh_closes_mass_without_demanding_sub_roundoff_residuals():
    d=Domain(4.75,resolution=64,extent=80.,cavity_radius=.9,cavity_depth=2.375,mesh='edge')
    laws=[Attachment('nucleation',1.,.006),Attachment('nucleation',1.,.002)]
    result=transfer_rates(d,AIR,42e-6,268.,.05,laws,14e-9)
    assert result['nonlinear_residual']<1e-9
    assert result['gas_water_residual']<1e-8
    assert result['transfer_mass_residual']<1e-10
