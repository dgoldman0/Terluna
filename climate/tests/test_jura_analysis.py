import numpy as np
import pytest

from climate.crm import jura_analysis as j


def test_coast_geometry_points_from_the_nearest_water_to_the_land():
    water = np.ones((9, 9), bool)
    water[2:7, 2:7] = False                                               # a square island in the sea
    distance, ex, ey, jw, iw = j.coast_geometry(water, 6000.0)
    assert distance[4, 4] == pytest.approx(18.0) and distance[2, 4] == pytest.approx(6.0) and distance[0, 0] == 0.0
    assert (ex[4, 2], ey[4, 2]) == (pytest.approx(1.0), pytest.approx(0.0))   # the west coast: the sea lies to the west
    assert (ex[6, 4], ey[6, 4]) == (pytest.approx(0.0), pytest.approx(-1.0))  # the north coast: onshore is south
    assert water[jw, iw].all()


def test_great_circle_distances_on_the_moon():
    r = 1737.4
    assert j.great_circle_km(0.0, 0.0, 0.0, 1.0, r) == pytest.approx(2 * np.pi * r / 360.0)
    assert j.great_circle_km(45.0, 10.0, 45.0, 10.0, r) == 0.0
    assert j.great_circle_km(45.0, 0.0, 45.0, 1.0, r) == pytest.approx(2 * np.pi * r / 360.0 * np.cos(np.radians(45.0)),
                                                                        rel=1e-4)
