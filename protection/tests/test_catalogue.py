"""The module catalogue reproduces the September report's per-cell and whole-screen figures."""
from __future__ import annotations
import unittest

from protection.modules import catalogue


class Catalogue(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.product = catalogue.build()
        cls.modules = {m['id']: m for m in cls.product['modules']}

    def test_cell_count_and_mass(self):
        # protection/report.md, sections 8 and 14: about 976,000 cells of about 5 million kg each.
        cell = self.modules['optical_cell']
        self.assertAlmostEqual(cell['count'] / 976_000, 1.0, places=2)
        self.assertAlmostEqual(cell['unit_mass_kg'] / 5.0e6, 1.0, places=3)

    def test_cell_with_its_share_of_holding(self):
        # Section 14: about 6.3 million kg with power and fuel, 183 MW mean power, 0.29 kg/s of fuel.
        pack, buffer = self.modules['holding_pack'], self.modules['propellant_buffer']
        total = self.modules['optical_cell']['unit_mass_kg'] + pack['unit_mass_kg'] + buffer['unit_mass_kg']
        self.assertAlmostEqual(total / 6.3e6, 1.0, places=1)
        self.assertAlmostEqual(pack['power_mean_w'] / 183e6, 1.0, places=2)
        self.assertAlmostEqual(pack['consumables_kg_s']['propellant'], 0.285, places=3)

    def test_materials_sum_to_unit_mass(self):
        for m in self.product['modules']:
            self.assertAlmostEqual(sum(m['materials_kg'].values()) / m['unit_mass_kg'], 1.0, places=12)

    def test_magnets(self):
        # Section 11: four installations, 2.14e14 kg ideal conductor and support in all.
        magnet = self.modules['regional_magnet_installation']
        self.assertEqual(magnet['count'], 4)
        self.assertAlmostEqual(magnet['count'] * magnet['unit_mass_kg'] / 2.14e14, 1.0, places=2)


if __name__ == '__main__':
    unittest.main()
