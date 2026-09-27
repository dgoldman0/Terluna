"""Sky-ship study: the air's equal-density heights, the calibration against LZ 129, the scans and the stored results."""
import json
import unittest

import numpy as np

from research.studies.sky_ships import run as sk


class SkyShipTests(unittest.TestCase):
    def test_equal_density_heights(self):
        self.assertAlmostEqual(sk.earth_height_of_density(1.225), 0.0, places=3)
        self.assertAlmostEqual(sk.earth_height_of_density(float(sk.earth_air(5000.0))), 5.0, places=2)
        self.assertAlmostEqual(sk.earth_height_of_density(1.4), -1.4, delta=0.05)    # below Earth's sea level
        self.assertIsNone(sk.earth_height_of_density(1.6))

    def test_crossing(self):
        x = np.array([1.0, 2.0, 3.0, 4.0])
        self.assertAlmostEqual(sk.crossing(x, [0.1, 0.3, 0.5, 0.7], 0.4), 2.5)
        self.assertIsNone(sk.crossing(x, [0.1, 0.2, 0.3, 0.35], 0.4))

    def test_calibration_reproduces_the_weight_statements(self):
        """Fitted to four rigid airships, the model gives each one's empty weight within 8% and its longitudinals and
        frames within 10%."""
        factors, checks = sk.calibrate()
        self.assertEqual(set(checks), set(sk.SHIPS))
        for name, c in checks.items():
            for key in ('empty_t', 'longitudinals_t', 'frames_t'):
                actual, model = c[key]
                self.assertLess(abs(model / actual - 1.0), 0.08 if key == 'empty_t' else 0.10, (name, key))
        self.assertAlmostEqual(checks['LZ 129 Hindenburg']['gross_lift_t'][1], 206.9, delta=0.1)

    def test_fit_in_relative_error(self):
        x = np.array([[1.0], [2.0], [4.0]])
        self.assertAlmostEqual(float(sk.fit(x, 3.0 * x[:, 0])[0]), 3.0, places=12)

    def test_similarity_in_the_scans(self):
        """The calibrated 1930s hull on the Moon, in air of Earth's sea-level density, has the Earth hull's weight-driven
        and aerodynamic shares at 6.04 times the length."""
        factors, _ = sk.calibrate()
        design = sk.bu.Design(**sk.DESIGN_1930, **factors)
        moon = sk.bu.Air(sk.G_MOON, 1.225, 288.15, sk.EARTH_MOLAR_MASS)
        lift = 0.95 * sk.UNIT_LIFT_KG_M3['hydrogen']
        e = sk.scan(sk.EARTH_SEA, design, lift, [245.0])
        m = sk.scan(moon, design, lift, [245.0 * sk.SIMILAR])
        self.assertAlmostEqual(float(m['weight_driven_share'][0] / e['weight_driven_share'][0]), 1.0, places=9)
        self.assertAlmostEqual(float(m['aerodynamic_share'][0] / e['aerodynamic_share'][0]), 1.0, places=9)
        self.assertGreater(float(m['useful_share'][0]), float(e['useful_share'][0]))

    def test_stored_results(self):
        path = sk.HERE / 'results' / 'sky_ships.json'
        if not path.is_file():
            self.skipTest('results have not been generated')
        r = json.loads(path.read_text())
        self.assertEqual(r['schema'], sk.SCHEMA)
        rho = [row['density_kg_m3'] for row in r['air']]
        self.assertEqual(rho, sorted(rho, reverse=True))
        band = [row for row in r['air'] if 35.0 <= row['height_km'] <= 45.0]
        self.assertTrue(all(0.6 < row['density_kg_m3'] < 0.85 for row in band))
        self.assertTrue(all(row['lift_helium_kg_m3'] < row['lift_hydrogen_kg_m3'] for row in r['air']))
        for name, case in r['buoyant']['cases'].items():
            u = np.array(case['useful_share'])
            best = int(np.argmax(u))
            self.assertTrue(np.all(np.isfinite(u)) and np.all(u < 1.0), name)
            self.assertTrue(np.all(np.diff(u[best:]) <= 1e-9), name)     # past the best size, bigger is worse


if __name__ == '__main__':
    unittest.main()
