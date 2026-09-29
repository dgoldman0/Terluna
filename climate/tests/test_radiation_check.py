"""Checks of the radiation comparison's pure pieces: the daily mean over Sun heights and the line-by-line column."""
import numpy as np
import pytest

from climate.crm import radiation_check as rad


def test_the_daily_mean_weights_each_sun_height_by_its_time():
    nodes, weights = np.polynomial.legendre.leggauss(6)
    hours, hour_weights = 0.25 * np.pi * (nodes + 1.0), 0.25 * np.pi * weights
    parts = [dict(reflected=0.25, air=0.30, ground=0.45)] * 6               # the same split at every Sun height
    day = rad.daily_mean(parts, hours, hour_weights, 1173.5)
    assert day['reflected'] == pytest.approx(0.25 * 1173.5 / np.pi, rel=1e-9)
    assert sum(day.values()) == pytest.approx(1173.5 / np.pi, rel=1e-9)   # all the sunlight of the day, 373.5 W/m2


def test_a_line_by_line_column_holds_the_water_it_is_given():
    p = np.array([100000.0, 90000.0, 50000.0, 10000.0])
    t = np.array([295.0, 290.0, 270.0, 230.0])
    col, water = rad.line_by_line_column(p[1:], t[1:], np.full(3, 0.01), p[0], t[0])
    mass = col['layer_mass_kg_m2'].sum()
    assert water == pytest.approx(0.01 * mass, rel=0.02)                   # vapour 1% of the mass above 10 kPa
    assert col['p_pa'][0] == p[0] and np.all(np.diff(col['p_pa']) < 0)
