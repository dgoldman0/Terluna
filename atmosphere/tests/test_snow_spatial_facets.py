"""Independent balances, geometric limits and published facet-law controls."""
import math
import numpy as np
import pytest
from atmosphere.snow.spatial_diffusion import Domain, Attachment
from atmosphere.snow.spatial_facet_growth import FacetGeometry, surface_balances, libbrecht_parameters, libbrecht_laws


def state():
    # Deliberately synthetic physical inputs for algebraic controls, not air data.
    return dict(ice_saturation_pressure_Pa=100.,vapor_gas_constant_J_kg_K=460.,
        vapor_molecular_speed_m_s=500.,latent_sublimation_J_kg=2.8e6,
        diffusivity_m2_s=2e-5,conductivity_W_m_K=.02)


@pytest.mark.parametrize('scales',[(.1,.2),(1.,1.),(.001,.7)])
def test_two_source_linear_laws_match_matrix_inverse(scales):
    matrix=np.array([[30.,5.],[7.,40.]])
    laws=[Attachment('constant',x) for x in scales]
    expected=np.linalg.solve(np.eye(2)+matrix@np.diag(scales),[.05,.05])
    np.testing.assert_allclose(surface_balances(matrix,.05,laws),expected,rtol=1e-10)


def test_common_source_identical_kinetics_force_equal_normal_rates():
    g=FacetGeometry(Domain(3.,resolution=12,mesh='edge',cavity_radius=.8,cavity_depth=1.5))
    ans=g.solve(state(),30e-6,230.,.3,[Attachment('harrington',critical=.1)]*2)
    assert ans['prism_rate_m_s']==pytest.approx(ans['basal_rate_m_s'],rel=1e-10)
    assert ans['source_balance_relative_residual']<1e-9
    assert ans['transport_relative_mass_residual']<1e-8


def test_finite_rim_area_and_integrated_mass_are_independent_of_cavity_depth():
    from atmosphere.snow.facet_growth import ICE_DENSITY
    for depth in (.3,1.,1.8):
        g=FacetGeometry(Domain(2.,resolution=12,mesh='edge',cavity_radius=.7,cavity_depth=depth))
        np.testing.assert_allclose(g.areas,[4*math.pi,.51*math.pi],rtol=1e-12)
        ans=g.solve(state(),20e-6,240.,.05,[Attachment('constant',.01)]*2,sources=('middle','middle'))
        volume=2*(20e-6)**2*np.dot(g.areas,[ans['prism_rate_m_s'],ans['basal_rate_m_s']])
        assert ans['mass_rate_kg_s']==pytest.approx(ICE_DENSITY*volume,rel=1e-12)


def test_zero_supply_gives_zero_growth_and_finite_diagnostics():
    g=FacetGeometry(Domain(resolution=8))
    ans=g.solve(state(),1e-5,240.,0.,libbrecht_laws(-10.))
    assert ans['mass_rate_kg_s']==0
    assert ans['temperature_increment_K']==0
    assert ans['source_balance_relative_residual']==0
    assert ans['passes_uniform_local_uptake_conditions']


def test_constant_attachment_small_particle_limit():
    from atmosphere.snow.facet_growth import ICE_DENSITY
    g=FacetGeometry(Domain(resolution=8))
    a=1e-14;T=240.;si=.01;alpha=.03
    ans=g.solve(state(),a,T,si,[Attachment('constant',alpha)]*2,heating=False)
    expected=100/(460*T)*125*si*alpha/ICE_DENSITY
    assert ans['prism_rate_m_s']==pytest.approx(expected,rel=1e-6)
    assert ans['basal_rate_m_s']==pytest.approx(expected,rel=1e-6)


def test_full_surface_can_reject_an_algebraically_closed_source_balance():
    g=FacetGeometry(Domain(5.,resolution=16,mesh='edge'))
    ans=g.solve(state(),50e-6,240.,.1,[Attachment('constant',1.)]*2)
    assert ans['source_balance_relative_residual']<1e-9
    assert ans['minimum_surface_supersaturation']<0
    assert not ans['passes_uniform_local_uptake_conditions']


def test_2023_peak_dip_strengths_and_percent_conversion():
    b=libbrecht_parameters(-4.5);p=libbrecht_parameters(-14.4)
    assert b['basal_critical']/libbrecht_parameters(-4.5,narrow=False)['basal_critical']==pytest.approx(.13)
    assert p['prism_critical']/libbrecht_parameters(-14.4,narrow=False)['prism_critical']==pytest.approx(.05)
    cold=libbrecht_parameters(-50.)
    assert cold['prism_critical']==pytest.approx(.3743449844,rel=1e-8)
    assert cold['basal_critical']==pytest.approx(.1910301274,rel=1e-8)
    assert cold['prism_critical']>cold['basal_critical']


def test_printed_log_ambiguity_remains_a_distinct_sensitivity():
    decimal=libbrecht_parameters(-10.);natural=libbrecht_parameters(-10.,logarithm='ln')
    assert decimal['prism_critical']==pytest.approx(.0059163,rel=1e-5)
    assert natural['prism_critical']>.014
    for t in (-.5,-5.,-15.,-50.):
        broad=libbrecht_parameters(t,narrow=False)
        narrow=libbrecht_parameters(t)
        assert 0<narrow['prism_critical']<=broad['prism_critical']
        assert 0<narrow['basal_critical']<=broad['basal_critical']


@pytest.mark.parametrize('tc',[-50.01,-.49,0.,float('nan')])
def test_law_does_not_extrapolate_past_the_published_figure(tc):
    with pytest.raises(ValueError):libbrecht_parameters(tc)


@pytest.mark.parametrize('kwargs',[dict(mesh='unknown'),dict(resolution=3)])
def test_bad_grid_rejected(kwargs):
    with pytest.raises(ValueError):Domain(**kwargs)


def test_edge_mesh_refines_thin_rim_and_preserves_exact_interfaces():
    counts=[];widths=[]
    for n in (16,32,64):
        d=Domain(2.,resolution=n,mesh='edge',cavity_radius=.99,cavity_depth=1.)
        counts.append(np.sum((d.r>.99)&(d.r<1)))
        widths.append(min(np.diff(d.rf)))
        assert .99 in d.rf and 1. in d.zf and 2. in d.zf
    assert counts==[4,8,16]
    assert widths[2]<widths[1]/3 and widths[1]<widths[0]/3


@pytest.mark.parametrize('matrix',[np.array([[1.,2.],[2.,1.]]),np.zeros((2,2)),np.array([[-1.,0.],[0.,1.]])])
def test_unsupported_resistance_matrix_is_rejected(matrix):
    with pytest.raises(ValueError):surface_balances(matrix,.1,[Attachment('constant')]*2)
