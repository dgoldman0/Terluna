"""Checks of the thunder study's own pieces (thunder.py): Earth's reference air and the flash's sources."""
import numpy as np
import pytest

from research.studies.atmospheric_electricity import thunder as tt


def test_earth_air_is_the_us_standard_atmosphere():
    air = tt.earth_air()
    t = lambda z: float(np.interp(z, air.z, air.t))
    assert t(0.0) == pytest.approx(288.15) and t(11000.0) == pytest.approx(216.65, abs=0.01)
    assert t(20000.0) == pytest.approx(216.65, abs=0.01) and t(32000.0) == pytest.approx(228.65, abs=0.01)
    assert float(np.interp(10000.0, air.z, air.p)) == pytest.approx(26436.0, rel=0.005)


def test_a_flash_s_sources_lie_a_kilometre_apart_along_its_channel():
    air = tt.earth_air()
    assert tt.sources(air, 4750.0, 10750.0).tolist() == [5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 11000.0]
    f_peak = tt.peak_frequency(air, 1.2e9, 4750.0, 10750.0, 27.0e6)
    assert 80.0 < f_peak < 150.0                                    # Earth's thunder peaks near 100 Hz
