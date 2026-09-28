"""Checks of the box comparison's pure pieces: the air as the GCM measures it and the variables it adds."""
import numpy as np
import pytest

from climate.crm import box_analysis as ba


def column(theta_of_p, pressures, surface=1.0e5):
    p = np.asarray(pressures, dtype=float)[:, None]
    return dict(prs=p, th=theta_of_p(p), psfc=np.array([surface]))


def test_the_gcm_layer_of_an_isothermal_column_is_carried_down_a_dry_adiabat():
    kappa = 287.04 / 1004.64
    d = column(lambda p: 280.0 * (1.0e5 / p) ** kappa, [99000.0, 97000.0, 95000.0, 90000.0, 80000.0])
    got = ba.gcm_layer_air_c(d, 0.9833)
    assert got[0] == pytest.approx(280.0 * 0.9833 ** -kappa - 273.15, abs=1e-6)


def test_the_gcm_layer_of_a_well_mixed_column_gives_the_air_at_the_ground():
    levels = 1.0e5 - np.arange(0.5, 40.0) * 250.0                         # 250-Pa layers through the bottom 10 kPa
    d = column(lambda p: np.full_like(p, 300.0), levels)                   # constant potential temperature
    assert ba.gcm_layer_air_c(d, 0.9833)[0] == pytest.approx(300.0 - 273.15, abs=0.01)


def test_the_box_adds_the_gcm_measure_and_rain_cover():
    assert ba.KEYS[-2:] == ('air_gcm_layer_c', 'raining') and 'air_c' in ba.KEYS
