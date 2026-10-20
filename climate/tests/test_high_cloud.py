"""The high cloud's range: its layer optics follow the GCM's ice optics, and its changes follow their signs."""
import numpy as np
import pytest

from climate.gcm import high_cloud as hc
from climate.gcm import plasim_ice_optics as io


def test_the_layer_optics_are_the_patchs_at_the_tables_radius():
    water, t, visible = 2.0, 192.0, 0.5
    reflection, emissivity = hc.layer_optics(water, float(io.ice_radius(t)), visible, nir_weights=io.NIR_WEIGHTS)
    expected = io.disk_mean(lambda m: visible * io.ice_sw(water, t, m)['r1']
                            + (1 - visible) * io.ice_sw(water, t, m)['r2'])
    assert reflection == pytest.approx(expected, rel=1e-9)
    assert emissivity == pytest.approx(float(io.emissivity(water, io.ice_kappa(t))), rel=1e-12)


def test_more_water_and_smaller_crystals_trap_more_infrared_and_reflect_more():
    cloud = dict(cover=[0.1, 0.3], water_g_m2=[1.0, 2.0], t_k=[165.0, 192.0], sunlight_w_m2=286.0,
                 planetary_albedo=0.3, olr_clear_w_m2=230.0)
    same = hc.effect(cloud, 0.5)
    assert same['net_w_m2'] == pytest.approx(0.0, abs=1e-12)
    for case in (hc.effect(cloud, 0.5, water_scale=2.0), hc.effect(cloud, 0.5, radius_scale=0.5)):
        assert all(x['infrared_w_m2'] > 0 and x['sunlight_w_m2'] < 0 for x in case['layers'])
    less = hc.effect(cloud, 0.5, water_scale=0.5)
    assert less['net_w_m2'] < 0 and np.allclose(less['ground_k'], [less['net_w_m2'] / f for f in hc.FEEDBACK_W_M2_K])
