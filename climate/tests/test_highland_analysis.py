"""Checks of the highland analysis's terrain pieces on synthetic ground."""
import numpy as np
import pytest

from climate.crm import highland_analysis as ha


def test_a_column_s_place_is_its_height_against_the_ground_around_it():
    h = np.full((30, 30), 2000.0)
    h[15, 15] = 3000.0                                                   # a summit
    h[5, 5] = 1000.0                                                     # a hollow
    place = ha.terrain_place(h, 6000.0, radius_m=18000.0, margin_m=150.0)
    assert place[15, 15] == 1 and place[5, 5] == -1 and place[25, 25] == 0
    assert ha.periodic_mean(np.arange(9.0).reshape(3, 3), 1)[1, 1] == pytest.approx(4.0)   # the cross about the centre


def test_slopes_point_uphill_on_the_periodic_box():
    x = np.arange(20) * 6000.0
    h = 1000.0 + 0.1 * np.tile(np.sin(2 * np.pi * x / x.size / 6000.0) * x.size * 6000.0 / (2 * np.pi), (20, 1))
    gx, gy = ha.slopes(h, 6000.0)
    assert np.allclose(gy, 0.0) and gx[0, 0] == pytest.approx(0.1, rel=0.05) and gx[0, 10] == pytest.approx(-0.1, rel=0.05)


def test_the_interior_leaves_out_the_blended_edges_and_stats_allow_an_empty_span():
    keep = ha.interior((64, 64), 6000.0, 50.0e3)
    assert keep.sum() == (64 - 2 * 9) ** 2 and not keep[8].any() and keep[9, 9]
    assert ha.stat(np.median, np.array([])) is None and ha.stat(np.percentile, [1.0, 3.0], 50) == pytest.approx(2.0)
    assert ha.line([0.0, 1.0, 2.0], [5.0, 4.0, 3.0]) == pytest.approx((-1.0, 5.0))


def test_departures_name_real_cases_and_point_to_their_reruns():
    from climate.crm.cm1_run import CASES
    assert set(ha.DEPARTURES) <= set(CASES)
    assert 'box_highland_own_height' in ha.DEPARTURES['box_highland'] and 'box_highland_own_height' in CASES
