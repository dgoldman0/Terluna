"""Groundwater estimate: the porosity profile meets its two constraints, and pore columns integrate it."""
import math
import unittest

from geography import groundwater as gw


class GroundwaterTests(unittest.TestCase):
    def test_profile_meets_both_constraints(self):
        phi0, d = gw.fit_profile()
        (mean, over), (deep, depth) = gw.MEAN_POROSITY_UPPER, gw.POROSITY_AT_DEPTH
        self.assertAlmostEqual(gw.pore_column_m(over, phi0, d) / over, mean, places=4)
        self.assertAlmostEqual(phi0 * math.exp(-depth / d), deep, places=4)

    def test_pore_column(self):
        self.assertAlmostEqual(gw.pore_column_m(1e9, 0.1, 5000.0), 500.0)          # phi0 * d when saturated throughout
        self.assertLess(gw.pore_column_m(1000.0, 0.1, 5000.0), 0.1 * 1000.0)         # porosity falls within the layer


if __name__ == '__main__':
    unittest.main()
