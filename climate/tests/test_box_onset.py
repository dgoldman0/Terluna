"""Why a box convects (box_onset.py), its parcel pieces on made-up soundings: Bolton's saturation, a parcel that keeps
its moist static energy when it condenses, the share of each level inside a mixed layer, and a conditionally unstable
sounding whose CAPE falls as the parcel mixes in more of the surrounding air."""
import numpy as np
import pytest

from climate.crm import box_onset as bo
from climate.crm import cm1_output as co

G = 1.6242


def made_up_sounding():
    """Levels 100 m apart to 60 km: a well-mixed moist layer below 2 km, then air cooling 1.4 K/km, drying with height
    to 40% relative humidity, in hydrostatic balance from 1215 hPa at lunar gravity."""
    z = np.arange(50.0, 60000.0, 100.0)
    theta_bl = 300.0
    t = np.where(z < 2000.0, theta_bl - G / co.CP * z, theta_bl - G / co.CP * 2000.0 - 1.4e-3 * (z - 2000.0))
    t = np.maximum(t, 180.0)
    lnp = np.log(121500.0) - np.concatenate([[0.0], np.cumsum(np.diff(z) * G / (co.RD * 0.5 * (t[1:] + t[:-1])))])
    p = np.exp(lnp)
    q = np.where(z < 2000.0, 0.017, 0.4 * bo.saturation(t, p))
    return z, p, t, q, np.full(z.size, 100.0)


def test_saturation_follows_bolton():
    es = 611.2 * np.exp(17.67 * 26.85 / (300.0 - 29.65))
    assert bo.saturation(300.0, 1.0e5) == pytest.approx(bo.EPSV * es / (1.0e5 - es))


def test_a_parcel_keeps_its_moist_static_energy_when_it_condenses():
    h, z, p = co.CP * 300.0 + G * 5000.0 + bo.LV * 0.03, 5000.0, 9.0e4
    t, q, wet = bo.saturate(h, 0.03, z, p, G)
    assert wet and q == pytest.approx(bo.saturation(t, p)) and q < 0.03
    assert co.CP * t + G * z + bo.LV * q == pytest.approx(h, rel=1e-9)
    t_dry, q_dry, wet_dry = bo.saturate(co.CP * 300.0 + G * z + bo.LV * 1e-4, 1e-4, z, p, G)
    assert not wet_dry and t_dry == pytest.approx(300.0) and q_dry == 1e-4


def test_each_level_counts_its_share_of_the_mixed_layer():
    assert bo.layer_share([0.0, 100.0, 300.0, 600.0], 200.0).tolist() == [1.0, 0.5, 0.0]
    assert bo.layer_share([0.0, 100.0, 300.0, 600.0], 0.0).tolist() == [0.0, 0.0, 0.0]


def test_a_conditionally_unstable_sounding_loses_cape_to_mixing():
    z, p, t, q, dz = made_up_sounding()
    mass = 1.4 * dz
    plain = bo.parcel(z, p, t, q, mass, dz, 500.0, 0.0, G)
    assert plain['cape'] > 100.0 and plain['cin'] >= 0.0
    assert plain['lcl_m'] <= plain['lfc_m'] < plain['el_m']
    mixed = [bo.parcel(z, p, t, q, mass, dz, 500.0, eps, G)['cape'] for eps in (1.0, 4.0, 16.0)]
    assert plain['cape'] > mixed[0] > mixed[1] > mixed[2] >= 0.0
    dry = bo.parcel(z, p, t, np.full(z.size, 1e-5), mass, dz, 0.0, 0.0, G)
    assert dry['cape'] == 0.0 and dry['lfc_m'] is None
