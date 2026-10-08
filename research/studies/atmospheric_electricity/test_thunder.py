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


def test_a_ground_strike_s_channel_reaches_the_ground(tmp_path, monkeypatch):
    case = tmp_path / 'run'
    case.mkdir()
    row = '{t:.1f} 10 1 {kind} 1000.0 2000.0 35000.0 1.2e5 1.0e5 1.0 3.0e7 1 50.0 50.0 2.0e10 1.0e10 {lo} 40000.0 40 1.0\n'
    (case / 'terluna_flashes.txt').write_text('# header\n# run from 0.0\n' + row.format(t=60.0, kind=1, lo=30000.0) +
                                              row.format(t=66.0, kind=2, lo=31000.0) + row.format(t=70.0, kind=1, lo=-1.0))
    monkeypatch.setattr(tt, 'RUNS', tmp_path)
    monkeypatch.setattr(tt, 'FLASH_SETS', (('run', 0.0, None),))
    f = tt.flashes()
    assert f['ground'].tolist() == [False, True, False]
    assert f['z_low_m'].tolist() == [30000.0, 0.0, 35000.0]          # in cloud as logged; to the ground; at its start
    assert f['energy_j'].tolist() == [1.0e10] * 3
