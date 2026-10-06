"""The display filter must preserve the clear control and constant spectra."""
import importlib.util
from pathlib import Path
import numpy as np

spec = importlib.util.spec_from_file_location('sea_cloud_display',Path(__file__).with_name('cloud_display.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_clear_cloud_control_stays_identically_zero():
    e=np.array([.25,10,30,60,89]);a=np.array([0,90,180,240,270,300,330.])
    sky=np.random.default_rng(17).uniform(.1,2,(len(e),len(a),3))
    f=module.delta_interpolator(e,a,sky,sky,4.)
    assert np.array_equal(f([[1,0],[45,359.9],[88,1]]),np.zeros((3,3)))


def test_constant_spectral_change_survives_filter_and_azimuth_seam():
    e=np.array([.25,8,26,55,89]);a=np.array([0,90,180,240,270,300,330.])
    sky=np.ones((len(e),len(a),3));change=np.array([.12,-.08,.03])
    f=module.delta_interpolator(e,a,sky+change,sky,4.)
    assert np.allclose(f([[1,0],[45,359.9],[88,1]]),change,atol=1e-12)
    e2,a2=np.meshgrid(e,a,indexing='ij')
    cloud=sky+np.cos(np.radians(a2))[...,None]*.1
    f=module.delta_interpolator(e,a,cloud,sky,4.)
    assert np.allclose(f([[40,0]]),f([[40,360]]),atol=1e-12)
