"""Checks on the conservation screens, the site classification and the register."""
from __future__ import annotations
import json
import tempfile
import unittest
from pathlib import Path

from geography import fetch_inputs
from research.studies.conservation import screens as sc
from research.studies.conservation.run import ROOT, PRODUCT, classify, run


def product_present():
    try:
        share = json.loads(PRODUCT.read_text())['share']
        grid = ROOT / 'geography' / 'products' / f'atlas_{round(100 * share)}pct_4ppd.npz'
        return grid.is_file() and fetch_inputs.restore()['status'] == 'PASS'
    except (OSError, ValueError, KeyError):
        return False


class ScreenTests(unittest.TestCase):
    def test_one_atmosphere_per_62_metres_of_lunar_water(self):
        self.assertAlmostEqual(sc.pressure_at_depth(62.4, 0.0) / sc.ATM, 1.0, delta=0.01)

    def test_nitrogen_frost_reproduces_pluto(self):
        # New Horizons measured about 1 Pa over Pluto's nitrogen ice at about 37 K.
        self.assertGreater(sc.frost_pressure('N2', 37.0), 0.5)
        self.assertLess(sc.frost_pressure('N2', 37.0), 2.0)

    def test_saltation_threshold_on_earth(self):
        # Sand on Earth starts moving at a friction speed near 0.2 m/s (Shao and Lu 2000).
        earth = sc.saltation_threshold(101325.0, 288.0, g=9.80665, grain_density=2650.0)
        self.assertGreater(earth['u_star_m_s'], 0.18)
        self.assertLess(earth['u_star_m_s'], 0.30)

    def test_thinner_air_needs_stronger_wind(self):
        winds = [sc.saltation_threshold(p * sc.ATM)['wind_10m_m_s'] for p in (0.01, 0.1, 1.2)]
        self.assertEqual(winds, sorted(winds, reverse=True))

    def test_runoff_scaling_is_consistent(self):
        s = sc.runoff_scaling()
        self.assertAlmostEqual(s['bedload_far_above_threshold'], 1.0, places=12)
        self.assertAlmostEqual(s['shields_number'], s['flow_depth'], places=12)
        self.assertAlmostEqual(s['stream_power'] * s['cohesive_slope_height'], 1.0, places=12)

    def test_dome_shell_takes_the_governing_limit(self):
        buckling = sc.dome_shell_m(1e5, 200.0, 200e9, 200e6)
        self.assertAlmostEqual(buckling, 200.0 * (1e5 / (0.2 * 2 / (3 * 0.91) ** 0.5 * 200e9)) ** 0.5, places=9)
        membrane = sc.dome_shell_m(5e7, 1.0, 200e9, 200e6)
        self.assertAlmostEqual(membrane, 5e7 * 1.0 / (2 * 200e6), places=12)

    def test_kardashev_scale(self):
        self.assertEqual(sc.kardashev_watts(1.0), 1e16)


class ClassificationTests(unittest.TestCase):
    def test_submerged_shore_and_dry(self):
        self.assertEqual(classify(-1700.0, -1654.0), 'submerged')
        self.assertEqual(classify(-1600.0, -1654.0), 'shore')
        self.assertEqual(classify(-1200.0, -1654.0), 'dry')


class ScopeTests(unittest.TestCase):
    def test_immersion_output_rejected(self):
        with self.assertRaises(ValueError):
            run(ROOT / 'immersion' / 'forbidden-conservation-output')
        self.assertFalse((ROOT / 'immersion' / 'forbidden-conservation-output').exists())


@unittest.skipUnless(product_present(), 'geography atlas product or inputs missing')
class RegisterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with tempfile.TemporaryDirectory() as tmp:
            cls.reg, cls.targets = run(Path(tmp))

    def test_every_site_is_placed(self):
        sites = json.loads((Path(__file__).with_name('sites.json')).read_text())['sites']
        self.assertEqual([r['id'] for r in self.reg], [s['id'] for s in sites])

    def test_tranquility_base_lies_under_the_sea(self):
        a11 = next(r for r in self.reg if r['id'] == 'apollo-11')
        self.assertEqual(a11['status'], 'submerged')
        self.assertTrue(400 < a11['depth_m'] < 550)
        self.assertEqual(a11['evaluation_status'], 'decided')

    def test_highland_landing_stays_dry(self):
        a16 = next(r for r in self.reg if r['id'] == 'apollo-16')
        self.assertEqual(a16['status'], 'dry')
        self.assertGreater(a16['above_sea_m'], 1000)

    def test_all_targets_found_in_the_gazetteer(self):
        names = [t['name'] for t in json.loads(Path(__file__).with_name('targets.json').read_text())['targets']]
        self.assertEqual(sorted(t['name'] for t in self.targets), sorted(names))


if __name__ == '__main__':
    unittest.main()
