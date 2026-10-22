"""Analytic limits, independent geometry balances and physical rejection checks."""
import math
import numpy as np
import pytest
from atmosphere.snow import uniform_face_growth as model, uniform_face_assessment as audit
from atmosphere.snow import hollow_shape as shape, facet_attachment as attachment, facet_growth as baseline


@pytest.fixture(scope='module')
def libraries():
    if not (shape.BUILD/'build.json').exists() or not (baseline.BUILD/'build.json').exists():
        pytest.skip('Optional pinned attachment and shape builds required')
    return attachment.AttachmentLibrary(baseline.BUILD), shape.ShapeLibrary()


def test_independent_faces_have_analytic_low_supersaturation_roots():
    si=1e-5;critical=[.1,.2];matrix=np.diag([1e8,2e9])
    result=model.coupled_surfaces(matrix,si,critical,[1.,1.])
    # In this range tanh(critical/surface) rounds to one; each equation is
    # exactly a quadratic, evaluated here in cancellation-free form.
    expected=[2*si/(1+math.sqrt(1+4*length*si/c)) for length,c in zip(np.diag(matrix),critical)]
    assert result['surface_supersaturations']==pytest.approx(expected,rel=1e-10)


def test_common_corner_has_an_analytic_quadratic_limit():
    si=1e-5;critical=.1;la,lc=1e7,3e8
    result=model.coupled_surfaces([[la,lc],[la,lc]],si,[critical,critical],[1.,1.])
    expected=2*si/(1+math.sqrt(1+4*(la+lc)*si/critical))
    assert result['surface_supersaturations']==pytest.approx([expected,expected],rel=1e-10)


def test_exchange_of_face_labels_preserves_solution():
    matrix=np.array([[420.,180.],[90.,620.]])
    first=model.coupled_surfaces(matrix,.07,[.003,.009],[1.,10.])
    second=model.coupled_surfaces(matrix[::-1,::-1],.07,[.009,.003],[10.,1.])
    assert first['surface_supersaturations']==pytest.approx(second['surface_supersaturations'][::-1],rel=1e-9)
    assert first['efficiencies']==pytest.approx(second['efficiencies'][::-1],rel=1e-8)


@pytest.mark.parametrize('matrix', [[[0,1],[1,2]],[[1,3],[3,1]],[[1,-1],[0,1]],[[1,float('nan')],[0,1]],[[1,2,3]]])
def test_rejects_unsupported_transport_systems(matrix):
    with pytest.raises(ValueError):model.coupled_surfaces(matrix,.1,[.01,.01],[10.,10.])


def test_zero_forcing_still_validates_mechanisms():
    with pytest.raises(ValueError):model.coupled_surfaces([[1,0],[0,1]],0.,[.01,.01],[0.,10.])
    result=model.coupled_surfaces([[1,0],[0,1]],0.,[.01,.01],[10.,10.])
    assert result['surface_supersaturations']==[0.,0.]


def test_original_1996_table_independently_constrains_geometry(libraries):
    _,shapes=libraries
    # NB96 Table B1; decimal values are rounded in the original paper.
    reference=[(.01,.0233,62.),(.05,.09,11.7),(.1,.155,5.62),(.2,.258,2.65),
        (.25,.300,2.07),(.5,.469,.959),(1.,.690,.445),(2.,.963,.208),
        (5.,1.38,.0781),(10.,1.70,.0374),(20.,2.05,.0184)]
    # Published analytic fits approximate the tabulated numerical solution;
    # their largest difference is 5.46%, at phi=0.01, not roundoff accuracy.
    errors=[]
    for phi,ha,hc in reference:
        result=shapes.evaluate(phi,'cylinder')
        errors.extend([abs(result['prism_edge']/ha-1),abs(result['basal_edge']/hc-1)])
    assert max(errors)==pytest.approx(.0545907518,rel=1e-6)


@pytest.mark.parametrize('mechanism', [1.,3.,10.])
@pytest.mark.parametrize('rim', [None,.01,.1])
def test_equal_faces_at_common_corner_have_equal_normal_speeds(libraries,mechanism,rim):
    library,shapes=libraries;a,c=30e-6,150e-6
    state=library.evaluate(.48,223.15,97000.,a,c,mechanism)
    assert state['characteristic_supersaturation_a']==state['characteristic_supersaturation_c']
    result=model.solve(state,shapes,a,c,223.15,.48,mechanisms=(mechanism,mechanism),thin_basal_rim=rim)
    assert result['facet_rate_ratio']==pytest.approx(1.,rel=1e-8)


def test_dense_cylinder_flux_closes_independent_volume_derivative(libraries):
    library,shapes=libraries;a,c=30e-6,90e-6
    state=library.evaluate(.003,263.15,70000.,a,c,10.)
    result=model.solve(state,shapes,a,c,263.15,.003)
    da,dc=result['equatorial_rate_m_s'],result['polar_rate_m_s']
    step=.01*min(a/da,c/dc)
    volume=lambda aa,cc:2*math.pi*aa*aa*cc
    # Richardson eliminates the cubic remainder in this exact polynomial volume.
    derivative=lambda dt:(volume(a+da*dt,c+dc*dt)-volume(a-da*dt,c-dc*dt))/(2*dt)
    exact=(4*derivative(step/2)-derivative(step))/3
    assert baseline.ICE_DENSITY*exact==pytest.approx(result['mass_rate_kg_s'],rel=1e-10)


def test_source_thermal_field_responds_to_heat_conductance(libraries):
    library,shapes=libraries;a,c=40e-6,50e-6
    state=library.evaluate(.01,268.15,101000.,a,c,10.)
    first=model.solve(state,shapes,a,c,268.15,.01)
    cooled=model.solve(state,shapes,a,c,268.15,.01,heat_ventilation=2.)
    assert cooled['mass_rate_kg_s']>first['mass_rate_kg_s']
    assert cooled['position_states']['edge']['latent_temperature_increment_K']<first['position_states']['edge']['latent_temperature_increment_K']
    edge=first['position_states']['edge']
    beta=(state['latent_sublimation_J_kg']/state['vapor_gas_constant_J_kg_K']/268.15-1)/268.15
    assert edge['surface_supersaturation']+edge['vapor_depletion']+beta*edge['latent_temperature_increment_K']==pytest.approx(.01)


def test_thin_ring_has_the_published_long_column_asymptote(libraries):
    _,shapes=libraries;a,c=5e-6,250e-6
    solid=model.geometry(shapes,a,c);rim=model.geometry(shapes,a,c,thin_basal_rim=.03)
    # hc approximately proportional to 1/phi: C1 tends to hc * rim_width/a.
    assert rim['source_lengths_m'][0,1]/solid['source_lengths_m'][0,1]==pytest.approx(.03,rel=.02)
    assert rim['active_basal_area_fraction']==pytest.approx(1-.97**2)
    assert rim['source_lengths_m'][0,0]==solid['source_lengths_m'][0,0]


@pytest.mark.parametrize('rim', [0.,-.01,.101,float('nan')])
def test_thick_rims_not_silently_used_as_thin_rings(libraries,rim):
    with pytest.raises(ValueError):model.geometry(libraries[1],30e-6,60e-6,thin_basal_rim=rim)


def test_thin_ring_does_not_invent_off_edge_diffusion_field(libraries):
    with pytest.raises(ValueError):model.geometry(libraries[1],30e-6,60e-6,sources=('middle','edge'),thin_basal_rim=.03)
    state=libraries[0].evaluate(.1,268.15,101000.,30e-6,60e-6,10.)
    result=model.solve(state,libraries[1],30e-6,60e-6,268.15,.1,thin_basal_rim=.03)
    assert set(result['position_states'])=={'edge'}
    assert result['passes_sampled_uniform_growth_conditions'] is None


def test_mass_budget_distinguishes_maximum_dimensions_from_growing_area():
    a,c,da,dc=30e-6,200e-6,1e-8,2e-7;rho=baseline.ICE_DENSITY
    prism=6*math.sqrt(3)*rho*a*c*da;basal=3*math.sqrt(3)*rho*a*a*dc
    full=audit.solid_mass_budget(a,c,da,dc,prism+basal)
    assert full['basal_active_fraction_with_full_prism']==pytest.approx(1.)
    hollow=audit.solid_mass_budget(a,c,da,dc,prism+.2*basal)
    assert hollow['basal_active_fraction_with_full_prism']==pytest.approx(.2)
    impossible=audit.solid_mass_budget(a,c,da,dc,.9*prism)
    assert not impossible['thin_tip_only_mass_budget_feasible']
    assert impossible['equivalent_rim_width_fraction'] is None


def test_equal_axis_invariant_against_actual_measurement_error_bounds(libraries):
    result=audit.cold_constraints(*libraries)
    assert [sum(x['inside_declared_thin_rim_regime'] for x in c['states']) for c in result['cases']]==[0,0]
    windows=[w for c in result['cases'] for w in c['windows']]
    assert [w['equal_axis_rate_excluded_by_reported_error_bounds'] for w in windows]==[True,True,False]
