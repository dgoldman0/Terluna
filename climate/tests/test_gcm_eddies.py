import numpy as np

from climate.crm import gcm_eddies as ge


def test_layer_depths_come_back_from_their_midpoints():
    half = np.array([0.0, 0.2, 0.5, 0.8, 1.0])
    depth, back = ge.sigma_depths((half[:-1] + half[1:]) / 2.0)
    assert np.allclose(back, half)
    assert np.allclose(depth, np.diff(half))


def test_the_gcm_levels_give_its_lowest_layer_a_thirtieth_of_the_column():
    sigma = np.array([0.0383, 0.1191, 0.2108, 0.3168, 0.4368, 0.5668, 0.6994, 0.8233, 0.9241, 0.9833])
    depth, half = ge.sigma_depths(sigma)
    assert abs(depth.sum() - 1.0) < 1e-3                                 # sigmas rounded to four places
    assert abs(half[5] - 0.5) < 2e-3 and abs(depth[-1] - 0.0334) < 1e-3
