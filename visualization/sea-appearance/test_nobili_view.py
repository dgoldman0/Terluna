"""Conserved, nonnegative Earth radiance for a low-elevation out-of-gamut source."""
from pathlib import Path
import importlib.util
import numpy as np

def test_warm_earth_normalization_preserves_xyz_without_negative_pixels():
    path=Path(__file__).with_name('nubium_human_view.py')
    spec=importlib.util.spec_from_file_location('nobili_earth_view_test',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    texture=np.array([[.02,.08,.9],[.8,.85,.9],[.12,.25,.08],[.6,.4,.15]])
    weights=np.array([.001,.002,.003,.001]);target=np.array([.3367,.2875,.03114])
    from research.studies.sea_appearance.scenes import XYZ_TO_SRGB
    assert (target@XYZ_TO_SRGB.T).min()<0  # the measured colour cannot fit sRGB
    direct=mod.normalize_positive_xyz(texture,weights,target)
    assert (direct>=0).all()
    np.testing.assert_allclose((direct*weights[:,None]).sum(axis=0),target,rtol=1e-12)
