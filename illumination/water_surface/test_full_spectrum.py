"""Full-range water spectra must have finite height and retain the short-wave limit."""
import numpy as np
import pytest
from illumination.water_surface.short_waves import curvature
from illumination.water_surface.full_spectrum import full_curvature
from illumination.water_surface.full_spectrum import unified_density
from shared.constants import MOON_SURFACE_GRAVITY as G


def test_height_variance_converges_as_low_wavenumber_limit_is_reduced():
    values=[]
    for low in [1e-2,1e-4,1e-6]:
        k=np.geomspace(low,2000,30000);b,c=full_curvature(k,1.1,.074,.84,G)
        values.append(np.trapezoid((b+c)/k**2,np.log(k)))
    np.testing.assert_allclose(values,values[0],rtol=1e-8)
    assert 0 < 4*np.sqrt(values[0]) < .3


def test_full_spectrum_retains_short_wave_asymptote():
    k=np.array([100.,200.,400.]);_,before=curvature(k,1.1,.074,.84,G)
    _,after=full_curvature(k,1.1,.074,.84,G)
    np.testing.assert_allclose(after,before,rtol=2e-4)


def test_directional_density_integrates_to_omnidirectional_spectrum():
    density=unified_density(1.1,.074,.84,37.,G)
    k=np.array([.1,1.,10.,100.]);phi=(np.arange(720)+.5)*2*np.pi/720
    actual=density(k[:,None]*np.cos(phi),k[:,None]*np.sin(phi)).mean(axis=1)*2*np.pi*k**4
    b,c=full_curvature(k,1.1,.074,.84,G)
    np.testing.assert_allclose(actual,b+c,rtol=1e-12)
    assert density(np.array([0.]),np.array([0.]))[0]==0
