"""Source radius must affect finite-disk horizon illumination and Monte Carlo light."""
import math
import numpy as np
import pytest
pytest.importorskip('numba')
from illumination.cloud_light.finite_source import direct_horizontal,estimate


def test_large_source_at_horizon_has_correct_projected_flux():
    edges=np.array([1000.,1100.]);beta=np.zeros((1,1,2));origin=np.array([0.,0.,1000.00001])
    radius=math.radians(1.)
    light=direct_horizontal(origin,np.array([1.,0.,0.]),edges,np.zeros((1,1)),np.ones((1,3)),
        edges,-.1,.2,beta,50.,64,256,source_radius_rad=radius)
    assert light[1]==pytest.approx(2*radius/(3*math.pi),rel=.001)


def test_point_source_monte_carlo_recovers_absorbing_lambertian_solution():
    edges=np.array([1000.,1100.]);beta=np.zeros((1,1,2))
    value,cuts=estimate(np.array([0.,0.,1000.01]),np.array([[0.,0.,-1.]]),np.array([0.,0.,1.]),
        edges,np.zeros((1,1)),np.array([[.01]]),np.ones((1,3)),edges,-.1,.2,beta,50.,
        .2,.85,.8,photons=160,groups=4,source_radius_rad=0.)
    assert value.mean(axis=1)[0,1]==pytest.approx(.2/math.pi*math.exp(-1),rel=2e-4)
    assert cuts.sum()==0


def test_new_driver_preserves_the_pinned_solar_baseline():
    from illumination.cloud_light import volume
    edges=np.array([1000.,1100.]);beta=np.zeros((1,1,2));origin=np.array([0.,0.,1000.01])
    source=np.array([0.,0.,1.]);directions=np.array([[0.,0.,-1.],[1.,0.,0.]])
    args=(origin,directions,source,edges,np.array([[.002]]),np.array([[.001]]),
        np.ones((1,3)),edges,-.1,.2,beta,50.,.2,.85,.8)
    before,cuts_before=volume.estimate(*args,photons=64,groups=4)
    after,cuts_after=estimate(*args,photons=64,groups=4)
    np.testing.assert_array_equal(before,after)
    np.testing.assert_array_equal(cuts_before,cuts_after)
    args=(origin,source,edges,np.array([[.003]]),np.ones((1,3)),edges,-.1,.2,beta,50.)
    np.testing.assert_array_equal(volume.direct_horizontal(*args),direct_horizontal(*args))
