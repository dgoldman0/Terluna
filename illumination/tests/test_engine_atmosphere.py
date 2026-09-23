"""Engine-sky parameters follow from the solver's atmospheres, and their error is stated."""
import json
import math
import unittest

from illumination.sky.atmospheres import EARTH, MOON
from illumination.sky.engine_atmosphere import (SCHEMA, direct_sun_check, engine_parameters,
                                                product, sha256, HERE)


class EngineAtmosphereTests(unittest.TestCase):
    def test_earth_matches_the_engine_defaults(self):
        # Unreal Engine's SkyAtmosphere Earth defaults (per km), from the same Bruneton values.
        p = engine_parameters(EARTH)
        for got, want in zip(p['rayleigh_scattering_per_km'], [0.005802, 0.013558, 0.0331]):
            self.assertAlmostEqual(got / want, 1, delta=1e-3)
        for got, want in zip(p['absorption_per_km'], [0.00065, 0.001881, 0.000085]):
            self.assertAlmostEqual(got / want, 1, delta=2e-3)
        self.assertEqual(p['rayleigh_scale_height_km'], 8.0)
        self.assertEqual(p['absorption_tent'], {'tip_altitude_km': 25.0, 'width_km': 15.0})

    def test_moon_follows_from_its_definition(self):
        e, m = engine_parameters(EARTH), engine_parameters(MOON)
        stretch = 9.80665 / 1.62
        self.assertAlmostEqual(m['bottom_radius_km'], 1737.4)
        self.assertAlmostEqual(m['rayleigh_scale_height_km'], 8 * stretch)
        self.assertAlmostEqual(m['atmosphere_height_km'], 80 * stretch)
        for a, b in zip(m['rayleigh_scattering_per_km'], e['rayleigh_scattering_per_km']):
            self.assertAlmostEqual(a / b, 1.2)
        for a, b in zip(m['absorption_per_km'], e['absorption_per_km']):
            self.assertAlmostEqual(a * stretch / b, 1)
        self.assertAlmostEqual(m['absorption_tent']['tip_altitude_km'], 25 * stretch)
        self.assertEqual(m['mie_scattering_per_km'], [0.0, 0.0, 0.0])

    def test_top_of_atmosphere_sun(self):
        p = engine_parameters(MOON)
        self.assertTrue(125_000 < p['sun_illuminance_lux'] < 140_000, p['sun_illuminance_lux'])
        self.assertEqual(max(p['sun_colour_linear_srgb']), 1.0)

    def test_the_three_wavelength_error_is_reported(self):
        rows = {r['sun_elevation_deg']: r for r in direct_sun_check(MOON)}
        # Sampled coefficients overstate low-Sun light badly in air this thick...
        self.assertGreater(rows[5]['sampled_luminance_ratio'], 1.5)
        # ...fitted ones stay close to the spectral calculation at every elevation.
        for r in rows.values():
            self.assertLess(abs(r['fitted_luminance_ratio'] - 1), 0.06, r['sun_elevation_deg'])

    def test_product_records_its_sources(self):
        p = json.loads(json.dumps(product()))
        self.assertEqual(p['schema'], SCHEMA)
        self.assertEqual(p['producer']['model_sha256'], sha256(HERE / 'atmospheres.py'))
        self.assertEqual(set(p['atmospheres']), {'moon', 'moon_no_ozone', 'earth'})
        self.assertTrue(all(math.isfinite(v) for v in p['atmospheres']['moon']['rayleigh_scattering_fitted_per_km']))


if __name__ == '__main__':
    unittest.main()
