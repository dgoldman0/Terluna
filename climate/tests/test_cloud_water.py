"""The mirror of PlaSim's cloud-water formula and the runner's cloud_water choices."""
import unittest
import numpy as np

from climate.gcm import cloud_water as cw
from shared.constants import STANDARD_GRAVITY, MOON_SURFACE_GRAVITY

S = np.arange(1, 11) / 10.0
HALF = 0.75 * S + 1.75 * S ** 3 - 1.5 * S ** 4           # PlaSim's 10-level boundaries
SIGMA = 0.5 * (np.concatenate([[0.0], HALF[:-1]]) + HALF)


def column(ps=121590.0):
    """A warm, moist column: temperature along a rough moist profile with a 200 K floor aloft."""
    ta = np.maximum(200.0, 296.0 * SIGMA ** 0.19)[None, :, None, None]
    q = (0.016 * SIGMA ** 3)[None, :, None, None]
    return ta, q, np.full((1, 1, 1), ps)


class CloudWaterTests(unittest.TestCase):
    def test_boundaries_recovered_from_full_levels(self):
        np.testing.assert_allclose(cw.half_levels(SIGMA), np.concatenate([[0.0], HALF]), atol=1e-12)

    def test_earth_path_on_the_moon_is_plasim_on_earth(self):
        # Counting heights and precipitable water as on Earth, with Earth's water path, the Moon's
        # clouds hold exactly what the same pressures and temperatures give on Earth.
        ta, q, ps = column()
        moon, _ = cw.cloud_water_path(ta, q, ps, SIGMA, MOON_SURFACE_GRAVITY, 'earth_path')
        earth, _ = cw.cloud_water_path(ta, q, ps, SIGMA, STANDARD_GRAVITY, 'plasim')
        np.testing.assert_allclose(moon, earth, rtol=1e-12)

    def test_full_column_is_six_times_earth_path(self):
        ta, q, ps = column(ps=20000.0)                     # a thin column keeps paths below the cap
        full, _ = cw.cloud_water_path(ta, q, ps, SIGMA, MOON_SURFACE_GRAVITY, 'full_column')
        earth, _ = cw.cloud_water_path(ta, q, ps, SIGMA, MOON_SURFACE_GRAVITY, 'earth_path')
        self.assertLess(full.max(), cw.MAX_PATH_G_M2)
        np.testing.assert_allclose(full, earth * STANDARD_GRAVITY / MOON_SURFACE_GRAVITY, rtol=1e-12)

    def test_plasim_formula_empties_the_moons_high_clouds(self):
        ta, q, ps = column()
        plasim, z = cw.cloud_water_path(ta, q, ps, SIGMA, MOON_SURFACE_GRAVITY, 'plasim')
        earth, _ = cw.cloud_water_path(ta, q, ps, SIGMA, MOON_SURFACE_GRAVITY, 'earth_path')
        upper = SIGMA < 0.5
        self.assertTrue(np.all(plasim[0, upper] < 0.01 * earth[0, upper]))
        self.assertGreater(plasim[0, -1], earth[0, -1])    # while the lowest layer, in a heavy column, is wetter
        self.assertGreater(z[0, 4], 25000.0)               # sigma 0.44 stands over 25 km up on the Moon

    def test_optical_depth_fit(self):
        self.assertAlmostEqual(float(cw.optical_depth(np.array(100.0))), 2.0 * np.log10(101.5) ** 3.9)


if __name__ == '__main__':
    unittest.main()
