"""The surface wind product's friction velocity and its daylit percentiles."""
import numpy as np

from climate.gcm import surface_wind as sw


def test_friction_velocity_from_stress_and_the_air_density():
    # 0.12 Pa of stress on air at 1,215.9 hPa and 300 K with R = 288.4 J/kg/K (rho 1.405 kg/m3): u* = 0.292 m/s.
    u = sw.friction_velocity(0.072, -0.096, 1215.9, 300.0, 288.4)
    assert np.isclose(u, np.sqrt(0.12 / (121590.0 / (288.4 * 300.0))))
    assert 0.29 < float(u) < 0.30


def test_daylit_percentiles_keep_only_the_daylit_intervals():
    values = np.arange(10, dtype=float)[:, None, None] * np.ones((10, 1, 2))
    czen = np.zeros((10, 1, 2))
    czen[5:, 0, 0] = 0.6                       # the first cell is daylit in the last five intervals only
    out, count = sw.daylit_percentiles(values, czen, percentiles=(50,))
    assert count[0, 0] == 5 and count[0, 1] == 0
    assert out[0, 0, 0] == 7.0                 # the median of 5..9
    assert np.isnan(out[0, 0, 1])              # no daylit interval


def test_gas_constant_reads_the_run_namelist(tmp_path):
    (tmp_path / 'planet_namelist').write_text(' &PLANET_NML\n GA = 1.624\n GASCON = 288.376720919906\n /\n')
    assert np.isclose(sw.gas_constant(tmp_path), 288.376720919906)
