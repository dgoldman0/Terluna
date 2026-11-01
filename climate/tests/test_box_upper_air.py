"""The upper-air diagnosis (box_upper_air.py) on made-up profiles: the box's levels shared out over PlaSim's layers by
pressure, the curved air's extra escape, an undiluted moist adiabat that keeps its moist static energy, and
zero-buoyancy plumes that fall below it as they entrain drier air."""
import numpy as np
import pytest

from climate.crm import box_onset as bo
from climate.crm import box_upper_air as bu
from climate.crm import cm1_output as co

G = 1.6242


def test_the_levels_share_out_over_the_layers_and_keep_their_total():
    pw = np.array([1000.0, 800.0, 600.0, 400.0, 200.0, 0.0])          # interfaces, ground first
    edges = np.array([0.0, 0.3, 0.7, 1.0])                             # PlaSim's half levels, top first
    per_level = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    out = bu.on_layers(per_level, pw, edges, per_level=True)
    assert out.sum() == pytest.approx(per_level.sum())
    assert out[2] == pytest.approx(1.0 + 0.5 * 2.0)                    # 1000-700 hPa: all of one level, half the next
    rate = np.full(5, 1e-5)                                            # K/s on every level
    assert bu.on_layers(rate, pw, edges, g=G).sum() == pytest.approx(1e-5 * co.CP * 1000.0 / G)
    means = bu.layer_means(np.array([10.0, 10.0, 20.0, 20.0, 20.0]), pw, edges)
    assert means[2] == pytest.approx((10.0 * 300.0) / 300.0) and means[0] == pytest.approx(20.0)
    assert bu.bands(np.arange(10.0))['0-14 km'] == pytest.approx(7.0 + 8.0 + 9.0)


def test_the_curved_air_escape_scales_with_the_fourth_power_and_stays_in_its_shells():
    shells = (np.array([0.0, 1000.0, 2000.0]), np.array([2.0, 4.0]), np.array([250.0, 250.0]))
    zw = np.array([0.0, 500.0, 1500.0, 2500.0])
    out = bu.cm1_escape(np.array([250.0, 250.0, 250.0]), zw, shells)
    assert out.tolist() == pytest.approx([1.0, 1.0 + 2.0, 2.0])
    warm = bu.cm1_escape(np.array([275.0, 275.0, 275.0]), zw, shells)
    assert warm[0] == pytest.approx(1.1 ** 4)


def sounding():
    z = np.arange(50.0, 60000.0, 250.0)
    t = np.maximum(300.0 - 1.2e-3 * z, 200.0)
    lnp = np.log(121500.0) - np.concatenate([[0.0], np.cumsum(np.diff(z) * G / (co.RD * 0.5 * (t[1:] + t[:-1])))])
    return z, np.exp(lnp), t


def test_the_undiluted_adiabat_keeps_its_moist_static_energy():
    z, p, _ = sounding()
    tp = bu.pseudo_adiabat(z, p, 300.0, 0.015, z[0], G)
    k = 120
    qs = bo.saturation(tp[k], p[k])
    assert co.CP * tp[k] + G * z[k] + bo.LV * qs == pytest.approx(co.CP * 300.0 + G * z[0] + bo.LV * 0.015, rel=1e-6)
    assert np.all(np.diff(tp) < 0.0)


def test_zero_buoyancy_plumes_fall_below_the_adiabat_as_they_mix_in_drier_air():
    z, p, _ = sounding()
    tp = bu.pseudo_adiabat(z, p, 300.0, 0.015, z[0], G)
    moist = bu.zero_buoyancy(z, p, 300.0, 0.015, z[0], G, 1e-4, np.ones(z.size))
    assert moist == pytest.approx(tp, abs=1e-3)                        # saturated surroundings dilute nothing
    dry = [bu.zero_buoyancy(z, p, 300.0, 0.015, z[0], G, eps, np.full(z.size, 0.7))[-40] for eps in (0.0, 5e-5, 2e-4)]
    assert dry[0] == pytest.approx(tp[-40], abs=1e-3) and dry[0] > dry[1] > dry[2]
    assert bu.saturated_temperature(co.CP * 280.0 + G * 1000.0 + bo.LV * bo.saturation(280.0, 8e4), 1000.0, 8e4, G,
                                    270.0) == pytest.approx(280.0, abs=1e-6)
