"""The shared constants load, derive consistently, and match the values the models used."""
import json
import unittest
from pathlib import Path

from shared import constants as K


class ConstantsTests(unittest.TestCase):
    def test_literals_preserved(self):
        # Values the models used before they moved here; outputs must not shift.
        self.assertEqual(K.MOON_RADIUS, 1_737_400.0)
        self.assertEqual(K.MOON_GM, 4.902800118e12)
        self.assertEqual(K.GAS_CONSTANT, 1.380649e-23 * 6.02214076e23)
        self.assertEqual(K.SYNODIC_MONTH_DAYS, 29.53059)
        self.assertEqual(K.SOLAR_CONSTANT, 1361.0)

    def test_types_are_floats(self):
        for name in ("MOON_RADIUS", "EARTH_RADIUS", "SOLAR_CONSTANT", "EARTH_MOON_DISTANCE"):
            self.assertIsInstance(getattr(K, name), float, name)

    def test_derived_gravity_and_recorded_discrepancy(self):
        self.assertAlmostEqual(K.MOON_SURFACE_GRAVITY, 1.62422, places=5)
        self.assertEqual(K.LEGACY_MOON_GRAVITY, 1.62)
        self.assertAlmostEqual(K.MOON_SURFACE_GRAVITY / K.LEGACY_MOON_GRAVITY - 1, 0.0026, places=4)

    def test_synodic_from_sidereal_and_year(self):
        year = 365.256363  # sidereal year, days; checks the stored months are consistent
        synodic = 1 / (1 / K.SIDEREAL_MONTH_DAYS - 1 / year)
        self.assertAlmostEqual(synodic, K.SYNODIC_MONTH_DAYS, places=2)

    def test_schema(self):
        data = json.loads(Path(K.__file__).with_name("constants.json").read_text())
        self.assertEqual(data["schema"], "terluna.shared.constants/1")


if __name__ == "__main__":
    unittest.main()
