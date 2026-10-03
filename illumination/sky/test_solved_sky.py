"""Conservation, geometry and spectral tests for the solved-column sky."""
import math

import numpy as np
import pytest

from illumination.cloud_light.model import MolecularColumn
from illumination.sky.colour_matching import cmf
from illumination.sky.solved_optics import adaptive_spectrum, reduce_spectrum
from illumination.sky.solved_transport import (SUN_RADIUS_RAD, layer_lengths, quadrature,
                                               solar_table, weights_path)
from illumination.sky.solver import ray_path, transport


@pytest.mark.parametrize('height,mu',[(0.,1.),(0.,0.),(10.,-.05),(10.,-.8),(90.,-.2),(90.,.6)])
def test_shell_columns_match_independent_nested_spheres(height,mu):
    radius = 1000.
    h = np.array([0.,5.,20.,60.,100.])
    col = MolecularColumn(radius,h,np.ones((4,1)),np.array([550.]))
    expected,ground = col.path(height,mu)
    actual,blocked = layer_lengths(radius+height,mu,radius+h)
    assert blocked == ground
    if not ground:
        np.testing.assert_allclose(actual,expected,atol=1e-10)


def test_vertical_beam_is_beer_lambert_and_scattering_weights_close():
    edges = np.array([1000.,1010.,1040.,1100.])
    sca = np.array([[.01,.02],[.02,.03],[.003,.01]])
    absorb = sca*.3
    p,*_ = ray_path(edges[0],1.,edges[0],edges[-1],edges,1.,1.)
    w,t = weights_path(p,edges,sca,absorb)
    expected = np.exp(-np.diff(edges)@(sca+absorb))
    np.testing.assert_allclose(t,expected,rtol=1e-12)
    np.testing.assert_allclose(w.sum(axis=0)*1.3+t,1.,rtol=1e-12)


def test_finite_sun_occultation_in_vacuum():
    edges = np.array([1000.,1100.])
    a = np.array([-2*SUN_RADIUS_RAD,0.,2*SUN_RADIUS_RAD,math.pi/2])
    b = solar_table(edges[:1],a,edges,np.zeros((1,1)))
    np.testing.assert_allclose(b[0,:,0,0],[0.,.5,1.,1.],atol=2e-15)
    assert b[0,1,0,1] > 0
    assert b[0,3,0,1] == pytest.approx(1.,abs=SUN_RADIUS_RAD**2)


def test_altitude_quadrature_integrates_ground_solid_angle_and_polynomials():
    radius = 1000.
    radii = np.array([radius,1010.,1200.])
    mu,w,*_ = quadrature(radii,radius,5,12)
    for i,r in enumerate(radii):
        horizon = -np.sqrt(1-(radius/r)**2)
        assert w[i,mu[i]<horizon].sum() == pytest.approx(1+horizon,abs=1e-13)
        for power in range(5):
            expected = 2/(power+1) if power%2==0 else 0.
            assert np.dot(w[i],mu[i]**power) == pytest.approx(expected,abs=1e-13)


def test_empty_transport_preserves_vacuum_and_moment_trace():
    # A layer with zero optical thickness and a black ground has zero diffuse light.
    from illumination.sky.solved_transport import build_paths
    r = np.array([1000.,1010.,1100.])
    a = np.radians([-90.,0.,90.])
    mu,w,cp,sp,dphi = quadrature(r,1000.,2,4)
    zero = np.zeros((2,1))
    paths = build_paths(r,mu,r,zero,zero,r)
    beam = solar_table(r,a,r,zero)
    prev = np.zeros((3,3,1,7))
    out = transport(prev,beam,r,a,r,a,mu,w,cp,sp,dphi,*paths,0.,True)
    assert np.max(abs(out)) == 0.
    # With scattering, the second-moment trace equals the angular intensity.
    sca = np.full((2,1),.001)
    paths = build_paths(r,mu,r,sca,zero,r)
    beam = solar_table(r,a,r,sca)
    out = transport(prev,beam,r,a,r,a,mu,w,cp,sp,dphi,*paths,.1,True)
    assert np.min(out[:,:,:,0]) >= 0.
    np.testing.assert_allclose(out[:,:,:,0],out[:,:,:,1:4].sum(axis=3),rtol=1e-12,atol=1e-15)


def synthetic_lines():
    wave = np.linspace(500.,509.99,200)
    energy = np.ones(len(wave))*.05
    # A narrow saturated absorption feature tests averaging across lines.
    tau = .002+20*np.exp(-((wave-505)/.04)**2)
    return dict(height=np.array([0.,1000.,5000.]),wavelength=wave,energy=energy,
                xyz=683*energy[:,None]*cmf(wave),absorption=np.array([.9*tau,.1*tau]),
                scattering=np.array([np.full(len(wave),.3),np.full(len(wave),.1)]))


def test_reduction_preserves_source_colour_and_adapts_to_a_saturated_line():
    fine = synthetic_lines()
    coarse = reduce_spectrum(fine,10.,1)
    reduced,audit = adaptive_spectrum(fine,1e6,tolerance=.005)
    np.testing.assert_allclose(reduced['energy'].sum(),fine['energy'].sum(),rtol=1e-13)
    np.testing.assert_allclose(reduced['xyz'].sum(axis=0),fine['xyz'].sum(axis=0),rtol=1e-13)
    assert len(reduced['energy']) > len(coarse['energy'])
    exact = np.exp(-(fine['absorption']+fine['scattering']).sum(axis=0))@fine['xyz'][:,1]
    naive = np.exp(-(coarse['absorption']+coarse['scattering']).sum(axis=0))@coarse['xyz'][:,1]
    value = np.exp(-(reduced['absorption']+reduced['scattering']).sum(axis=0))@reduced['xyz'][:,1]
    assert abs(value/exact-1) < .005
    assert abs(value-exact) < abs(naive-exact)/10
    assert audit['maximum_training_xyz_error_bound'] <= .005


def test_reduction_fails_if_channel_budget_cannot_meet_tolerance():
    with pytest.raises(ValueError,match='channel cap'):
        adaptive_spectrum(synthetic_lines(),1e6,tolerance=1e-8,max_channels=1)


def test_wavelength_sampled_reference_matches_absorbing_lambertian_ground():
    from illumination.sky.solved_reference import monte_carlo
    edges = np.array([1000.,1100.])
    scattering = np.zeros((2,1))
    absorption = np.array([[0.],[.01]])
    means,truncated = monte_carlo(edges,scattering,absorption,.2,math.pi/2,-math.pi/2,0.,
                                 40,500,921,False,np.array([.25,1.]))
    # The two sampled wavelengths have vertical optical depths zero and one.
    exact = .2/math.pi*(.25+.75*math.exp(-1))
    standard_error = means.std(ddof=1)/math.sqrt(len(means))
    assert abs(means.mean()-exact) < 4*standard_error
    assert truncated.sum()==0
