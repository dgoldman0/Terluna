import numpy as np
import pytest

from research.studies.cold_waters import run as R


def test_air_carried_up_is_colder_by_the_lapse_rate():
    assert R.carried(300.0, 1000.0, 3000.0, 1.0e-3) == pytest.approx(298.0)      # 2 km above the cell's ground
    assert R.carried(300.0, 1000.0, 0.0, 1.0e-3) == pytest.approx(301.0)         # 1 km below it


def test_lapse_rate_of_a_column_with_a_known_lapse():
    """A column cooling 6.5 K a km with height gives 6.5 K a km between the band's heights."""
    sigma = np.array([0.5, 0.7, 0.8, 0.9, 0.95, 0.99])
    rd, g = 287.0, 9.81
    t_surface = 290.0
    # an atmosphere of constant lapse G: T = Ts - G z and p/ps = (T/Ts)^(g/(Rd G))
    gl = 6.5e-3
    t = t_surface * sigma ** (rd * gl / g)
    ta = t[None, :, None, None] * np.ones((2, 1, 1, 1))
    hus = np.zeros_like(ta)
    lapse = R.lapse_rates(ta, hus, sigma, gas_constant=rd, gravity=g, levels_km=(0.9, 2.5))
    assert np.allclose(lapse, gl, rtol=0.02)


def test_nearest_cell_wraps_longitude():
    j, i = R.nearest_cell(np.array([10.0]), np.array([359.0]), np.array([20.0, 10.0, 0.0]), np.array([0.0, 90.0, 180.0, 270.0]))
    assert j[0] == 1 and i[0] == 0


def test_pond_under_a_cold_clear_sky_cools_below_its_air(monkeypatch):
    """With no sunlight, a dry sky's infrared and still air, a shallow pond loses heat below the air's temperature."""
    monkeypatch.setattr(R, 'air_constants', lambda: (288.4, 1008.9))
    n = 400
    f = R.L.Forcing(dt=3575.0, sw_down=np.zeros(n), lw_down=np.full(n, 330.0), air_k=np.full(n, 290.0),
                    q_air=np.full(n, 0.002), wind10=np.full(n, 0.5), pressure=np.full(n, 121590.0), o2_fraction=0.2)
    out = R.run_water(f, 0.5, R.POND_LAYER_M, 17.0)
    assert out['surface_t'][-1] < 290.0 - 273.15
