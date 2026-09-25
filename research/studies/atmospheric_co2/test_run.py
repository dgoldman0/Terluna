"""Carbon study: the design air's conversions, the leaf model against its sources, and the weathering law."""
import math
import unittest

from research.studies.atmospheric_co2 import run as c


class CarbonTests(unittest.TestCase):
    def setUp(self):
        air = c.design_air()
        self.o2_pa = air.fractions['O2'] * air.surface_pressure_pa

    def test_design_air_conversions(self):
        d = c.results()['design_air']
        self.assertAlmostEqual(d['ppm_per_pa'], 1e6 / d['pressure_pa'], places=3)
        self.assertAlmostEqual(d['co2_pa'], d['co2_ppm'] * 1e-6 * d['pressure_pa'], places=2)
        self.assertTrue(34.0 < d['gt_co2_per_pa'] < 37.0)       # about 0.94 kg per m2 per Pa over 3.79e13 m2

    def test_leaf_kinetics_match_bernacchi_at_25c(self):
        gamma, kc, ko = c.kinetics(25.0)                          # 42.75 and 404.9 umol/mol, 278.4 mmol/mol
        self.assertAlmostEqual(gamma, 4.275, delta=0.05)
        self.assertAlmostEqual(kc, 40.49, delta=0.5)
        self.assertAlmostEqual(ko / 1e3, 27.84, delta=0.3)

    def test_photosynthesis_rises_with_co2_and_halves_at_glacial_levels(self):
        values = [c.relative_photosynthesis(ca, 25.0, self.o2_pa) for ca in (18.0, 28.0, 41.0, 60.0)]
        self.assertEqual(values, sorted(values))
        self.assertAlmostEqual(values[2], 1.0)
        self.assertLess(values[0], 0.45)

    def test_higher_total_pressure_costs_a_few_pascals(self):
        ratio, match = c.pressure_penalty(25.0, self.o2_pa, c.design_air().surface_pressure_pa)
        self.assertTrue(0.9 < ratio < 1.0)
        self.assertTrue(c.EARTH_TODAY_PA < match < 46.0)

    def test_dessert_law(self):
        self.assertAlmostEqual(c.weathering_mol_km2_yr(1000.0, 26.0) / 1e6, 1.717, places=2)
        self.assertAlmostEqual(c.weathering_mol_km2_yr(1.0, 10.0) / c.weathering_mol_km2_yr(1.0, 0.0), math.exp(0.642))


if __name__ == '__main__':
    unittest.main()
