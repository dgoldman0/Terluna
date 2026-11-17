"""The mapped mare boundaries on the grid: the basalt's share of the Moon and known places inside and outside it."""
from __future__ import annotations
import unittest

import numpy as np

from geography import fetch_surface_inputs as fetch_inputs, mare_map as M


def archive_present():
    try:
        return fetch_inputs.path(M.ARCHIVE).is_file()
    except (OSError, ValueError, KeyError):
        return False


@unittest.skipUnless(archive_present(), 'the mare boundaries are not restored')
class MareMapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = M.basalt(4, left_lon=0.0)

    def at(self, lat, lon):
        return bool(self.grid[int((90.0 - lat) * 4), int((lon % 360.0) * 4)])

    def test_the_maria_cover_a_sixth_of_the_moon(self):
        # Maria cover about 16-17% of the Moon (Head 1976; Nelson et al. 2014's polygons give 16.2% at 16 px/deg).
        self.assertAlmostEqual(M.area_share(self.grid), 0.162, delta=0.004)

    def test_known_places(self):
        self.assertTrue(self.at(32.8, -15.6))       # Mare Imbrium
        self.assertTrue(self.at(8.5, 31.4))         # Mare Tranquillitatis
        self.assertTrue(self.at(27.3, 147.9))       # Mare Moscoviense, on the far side
        self.assertTrue(self.at(-33.7, 163.5))      # Mare Ingenii, in the South Pole-Aitken basin
        self.assertFalse(self.at(-43.3, -11.2))     # Tycho, in the southern highlands
        self.assertFalse(self.at(10.0, 180.0))      # the far side's highlands

    def test_the_longitude_origin_only_rolls_the_grid(self):
        self.assertTrue(np.array_equal(np.roll(M.basalt(4), -180 * 4, axis=1), self.grid))


if __name__ == '__main__':
    unittest.main()
