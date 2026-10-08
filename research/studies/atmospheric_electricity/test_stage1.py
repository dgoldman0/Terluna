"""Checks of stage 1's comparison of storm currents with the air's leakage."""
import numpy as np
import pytest

from research.studies.atmospheric_electricity import stage1 as s1
from atmosphere.electricity import conductivity as cd


def test_the_threshold_field_is_284_kv_m_at_standard_number_density():
    rho0 = s1.N0 * s1.AIR_MOLAR_MASS / s1.AVOGADRO
    assert s1.threshold_field(rho0) == pytest.approx(284.0e3)
    assert s1.threshold_field(rho0 / 2.0) == pytest.approx(142.0e3)


def test_the_share_above_a_current_reads_the_log_histogram():
    hist = dict(edges_log10_na_m2=[-2.0, -1.0, 0.0, 1.0], counts=[10, 30, 60])
    assert s1.exceeding(hist, 1.0) == pytest.approx(0.6)
    assert s1.exceeding(hist, 10 ** -0.5) == pytest.approx((15 + 60) / 100.0)
    assert s1.exceeding(dict(edges_log10_na_m2=[0, 1], counts=[0]), 1.0) == 0.0


def test_breakdown_comes_only_above_the_critical_current_and_without_leakage_takes_eps0_e_over_j():
    tau = 1.0e5
    assert s1.time_to_breakdown_h(0.5, 1.0, tau) is None
    j, j_crit, e = 1000.0, 0.01, 1.8e5
    sigma = j_crit * 1e-9 / e
    hours = s1.time_to_breakdown_h(j, j_crit, cd.relaxation_time(sigma))
    assert hours == pytest.approx(cd.EPSILON_0 * e / (j * 1e-9) / 3600.0, rel=1e-4)


def test_the_zone_is_weighted_by_band_depth():
    bands = [dict(depth_m=3.0, height_km=30.0, pressure_hpa=650.0, air_density_kg_m3=0.85, cloud_water_g_m3=0.15),
             dict(depth_m=1.0, height_km=40.0, pressure_hpa=540.0, air_density_kg_m3=0.72, cloud_water_g_m3=0.10),
             dict(depth_m=0.0, warm_c=-30.0)]
    z = s1.zone_mean(bands)
    assert z['height_km'] == pytest.approx(32.5) and z['cloud_water_g_m3'] == pytest.approx(0.1375)
