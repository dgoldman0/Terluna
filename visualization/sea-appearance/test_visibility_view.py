"""Display checks: monotonic luminance, declared anchor, finite gamut."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest


def load_view():
    path=Path(__file__).with_name('nubium_visibility_view.py')
    spec=importlib.util.spec_from_file_location('visibility_view',path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_global_curve_keeps_bright_detail_and_sky_anchor():
    v=load_view()
    y=np.r_[np.geomspace(1e-5,17000,400),.5291728]
    xyz=y[:,None]*np.array([.95047,1,1.08883])
    rgb,record=v.tone(xyz,18000)
    display_y=rgb@np.array([.2126729,.7151522,.0721750])
    assert np.all(np.diff(display_y[:-1])>0)
    assert display_y[-1] == pytest.approx(.035,abs=1e-8)
    # The old operator clipped all these physically distinct bright values.
    bright=y[:-1]>1000
    assert np.ptp(display_y[:-1][bright])>.2
    assert rgb.max()<1


def test_gamut_mapping_preserves_requested_luminance():
    v=load_view()
    xyz=np.array([[2,1,.01],[.01,1,8],[.1,1,.1],[1e4,8e3,1e2]])
    rgb,record=v.tone(xyz,18000)
    expected=(np.log1p(xyz[:,1]/record['knee_world_cd_m2']) /
              np.log1p(18000/record['knee_world_cd_m2']))**record['power']
    actual=rgb@np.array([.2126729,.7151522,.0721750])
    np.testing.assert_allclose(actual,expected,atol=2e-7,rtol=0)
    assert np.isfinite(rgb).all() and rgb.min()>=0 and rgb.max()<=1
