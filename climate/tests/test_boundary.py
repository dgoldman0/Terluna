"""GCM boundary: seas and rain-fed lakes combined on a coarse grid (no geography files needed)."""
import unittest
import numpy as np

from climate.gcm import boundary


class LakeTests(unittest.TestCase):
    def test_lakes_become_fractions_with_their_own_surface(self):
        height = np.array([[-2000.0, 500.0], [800.0, 1200.0]])          # coarse cells, 2 x 2 fine cells each
        sea = np.array([[True, False], [False, False]])
        fine_cols = 4
        # a lake over three fine cells of the top-right coarse cell, level 700 m; one over a sea cell is dropped
        cells = np.array([0 * fine_cols + 2, 0 * fine_cols + 3, 1 * fine_cols + 2, 0])
        levels = np.array([700.0, 700.0, 700.0, 900.0])
        wet, depth, surface = boundary.combine_lakes(height, sea, cells, levels, fine_cols, 2, -1654.0)
        np.testing.assert_allclose(wet, [[1.0, 0.75], [0.0, 0.0]])
        np.testing.assert_allclose(surface, [[-1654.0, 700.0], [-1654.0, -1654.0]])
        np.testing.assert_allclose(depth, [[346.0, 200.0], [0.0, 0.0]])


if __name__ == '__main__':
    unittest.main()
