"""Lunar Prospector's oxides on the grid: the maria's iron and titanium, the highlands' and the basin floor's, and the
highlands' iron with the maria left out of the smoothing."""
from __future__ import annotations
import unittest

import numpy as np

from geography import fetch_surface_inputs as fetch_inputs, mare_map as M, surface_composition as C


def table_present():
    try:
        return fetch_inputs.path(C.TABLE).is_file() and fetch_inputs.path(M.ARCHIVE).is_file()
    except (OSError, ValueError, KeyError):
        return False


@unittest.skipUnless(table_present(), 'the Lunar Prospector table or the mare boundaries are not restored')
class CompositionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.whole = C.grid(4, left_lon=0.0)
        mare = M.basalt(4, left_lon=0.0)
        cls.highland = C.grid(4, left_lon=0.0, exclude=mare)

    def at(self, field, name, lat, lon):
        return float(field[name][int((90.0 - lat) * 4), int((lon % 360.0) * 4)])

    def test_known_places(self):
        self.assertGreater(self.at(self.whole, 'feo', 32.8, -15.6), 15.0)        # Mare Imbrium
        self.assertGreater(self.at(self.whole, 'tio2', 8.5, 31.4), 5.0)          # Tranquillitatis, high in titanium
        self.assertLess(self.at(self.whole, 'feo', 10.0, 175.0), 6.0)            # the far side's highlands
        self.assertTrue(7.0 < self.at(self.whole, 'feo', -50.0, -170.0) < 13.0)  # the South Pole-Aitken floor

    def test_the_highlands_iron_leaves_the_maria_out(self):
        # Leaving nothing out changes nothing; leaving the maria out lowers the iron of the Taurus highlands between
        # Serenitatis and Crisium, whose pixels saw the maria, and keeps the open highlands' and the basin floor's.
        self.assertTrue(np.allclose(C.grid(4, left_lon=0.0, exclude=np.zeros((720, 1440), bool))['feo'], self.whole['feo']))
        self.assertLess(self.at(self.highland, 'feo', 27.0, 30.0), self.at(self.whole, 'feo', 27.0, 30.0) - 2.0)
        for lat, lon in ((10.0, 175.0), (-50.0, -170.0)):
            self.assertAlmostEqual(self.at(self.highland, 'feo', lat, lon), self.at(self.whole, 'feo', lat, lon), delta=0.3)
        self.assertTrue(np.isfinite(self.highland['feo']).all())


if __name__ == '__main__':
    unittest.main()
