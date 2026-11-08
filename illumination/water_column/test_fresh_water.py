"""Checks of fresh pure water against model.py's seawater and the IOCCG absorption composite."""
import numpy as np
import pytest

from illumination.water_column import fresh_water as F
from illumination.water_column import model as M


def test_fresh_water_scatters_a_salt_share_less_than_the_seawater_table():
    assert F.scattering_factor(38.4) == pytest.approx(1.0)
    assert F.scattering_factor(0.0) == pytest.approx(1 / 1.30)
    assert F.scattering_factor(F.SEAS_SALINITY) == pytest.approx(0.770, abs=0.001)
    w = np.arange(400.0, 701.0, 10.0)
    assert np.allclose(F.pure_water(w)[1], M.seawater(w)[1] * F.scattering_factor(F.SEAS_SALINITY))


def test_ultraviolet_absorption_is_pure_waters_and_joins_pope_and_fry():
    assert F.pure_water(np.array([300.0]))[0][0] == pytest.approx(0.0147)              # IOCCG Table 1.1
    assert F.pure_water(np.array([350.0]))[0][0] == pytest.approx(0.0071)
    assert F.pure_water(np.array([300.0]))[0][0] < 0.2 * M.seawater(np.array([300.0]))[0][0]
    w = np.arange(420.0, 701.0, 5.0)
    assert np.allclose(F.pure_water(w)[0], M.seawater(w)[0])
    near = F.pure_water(np.array([415.0, 419.9, 420.0]))[0]
    assert abs(near[1] - near[2]) < 1e-5
    assert np.all(F.pure_water(np.array([300.0, 350.0]), uv_absorption=False)[0] == 0)


def test_water_is_the_models_water_over_fresh_pure_water():
    w = np.arange(420.0, 701.0, 10.0)
    fresh = F.water(w, 0.3, 0.02, 1.0, "developed")
    sea = M.water(w, 0.3, 0.02, 1.0, "developed")
    assert np.allclose(fresh["a"], sea["a"])
    b_w = M.seawater(w)[1]
    assert np.allclose(fresh["bb"], sea["bb"] - 0.5 * b_w * (1 - F.scattering_factor(F.SEAS_SALINITY)))
