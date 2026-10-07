"""Eddy mixing by breaking gravity waves on Earth and the Moon (atmosphere/loss_response/gravity_waves.py)."""
from __future__ import annotations
import json
import math
import unittest
from unittest import mock

import numpy as np

from shared.provenance import constants_changed
from atmosphere.loss_response import gravity_waves as gw


def isothermal(t=200.0, g=1.62, molar=0.0289, top_efolds=20.0, levels=2001, surface_pa=1e5):
    """An isothermal column with constant gravity: its scale height and buoyancy frequency are exact."""
    h = gw.GAS_CONSTANT * t / (molar * g)
    z = np.linspace(0.0, top_efolds * h, levels)
    return dict(z_m=z, p_pa=surface_pa * np.exp(-z / h), t_k=np.full_like(z, t), g_m_s2=np.full_like(z, g),
                molar_kg_mol=molar), h


class StandardAtmosphere(unittest.TestCase):
    def test_the_column_reproduces_the_standards_tables(self):
        col = gw.us76_column()
        z = col['z_m'] / 1e3
        table = {11.019: (216.65, 22632.1), 32.162: (228.65, 868.019), 71.802: (214.65, 3.95642),
                 86.0: (186.87, 0.37338), 100.0: (195.08, 0.032011), 110.0: (240.0, 0.0071042),
                 120.0: (360.0, 0.0025382), 150.0: (634.39, 4.5422e-4)}
        for km, (t, p) in table.items():
            self.assertAlmostEqual(float(np.interp(km, z, col['t_k'])), t, delta=0.2)
            self.assertAlmostEqual(float(np.exp(np.interp(km, z, np.log(col['p_pa'])))) / p, 1.0, delta=0.02)


class WaveTheory(unittest.TestCase):
    def setUp(self):
        self.still = mock.patch.multiple(gw, viscosity=lambda t: 0.0 * t, conductivity=lambda t: 0.0 * t)

    def test_a_saturated_hydrostatic_wave_needs_lindzens_diffusion(self):
        col, h = isothermal()
        n = col['g_m_s2'][0] / math.sqrt(1005.0 * col['t_k'][0])
        c, k = 20.0, 2 * math.pi / 1e6
        with self.still:
            total, waves = gw.diffusion(col, [(c, k, 1e-4)], col['p_pa'][10])
        broke = waves[0]['breaks_at_pa']
        self.assertIsNotNone(broke)
        above = (col['p_pa'] < broke / 3) & (col['p_pa'] > col['p_pa'][-1] * 10)
        lindzen = k * c ** 4 / (2 * h * n ** 3)
        self.assertTrue(np.allclose(total[above] / lindzen, 1.0, atol=0.01))

    def test_a_wave_holds_its_flux_until_it_breaks(self):
        col, h = isothermal()
        n = col['g_m_s2'][0] / math.sqrt(1005.0 * col['t_k'][0])
        c, k, flux = 20.0, 2 * math.pi / 1e6, 1e-4
        launch = col['p_pa'][10]
        with self.still:
            _, waves = gw.diffusion(col, [(c, k, flux)], launch)
        m = math.sqrt(gw.vertical_wavenumber(k, c, n, 1.0 / h))
        rho_break = 2 * flux * m / (k * c ** 2)
        p_break = rho_break * gw.GAS_CONSTANT * col['t_k'][0] / col['molar_kg_mol']
        self.assertAlmostEqual(waves[0]['breaks_at_pa'] / p_break, 1.0, delta=0.02)

    def test_a_wave_faster_than_the_buoyancy_allows_turns_back(self):
        col, h = isothermal()
        n = col['g_m_s2'][0] / math.sqrt(1005.0 * col['t_k'][0])
        k = 2 * math.pi / 1e5
        total, waves = gw.diffusion(col, [(1.2 * n / k, k, 1e-4)], col['p_pa'][10])
        self.assertIsNotNone(waves[0]['turns_at_pa'])
        self.assertFalse(total.any())

    def test_molecular_damping_ends_the_turbulence(self):
        col, h = isothermal(top_efolds=30.0)
        n = col['g_m_s2'][0] / math.sqrt(1005.0 * col['t_k'][0])
        c, k = 20.0, 2 * math.pi / 1e6
        total, _ = gw.diffusion(col, [(c, k, 1e-4)], col['p_pa'][10])
        q = gw.properties(col)
        quiet = q['nu'] + q['kappa'] > k * c ** 4 / (h * n ** 3) * 1.05
        self.assertTrue(quiet.any())
        self.assertFalse(total[quiet].any())
        self.assertTrue(total[~quiet].any())

    def test_the_calibration_gives_earths_homopause_mixing(self):
        factor, raw = gw.calibration()
        col, total, _ = gw.earth()
        at = float(np.interp(gw.EARTH_HOMOPAUSE_KM * 1e3, col['z_m'], total)) * factor * 1e4
        self.assertAlmostEqual(at / gw.EARTH_HOMOPAUSE_K_CM2_S, 1.0, places=9)
        self.assertTrue(0 < factor < 1)


@unittest.skipUnless(gw.OUT.is_file(), 'gravity_waves.json not written')
class ProductTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.product = json.loads(gw.OUT.read_text())

    def test_product_is_current_with_its_code_and_products(self):
        self.assertEqual(self.product['schema'], gw.SCHEMA)
        for name, value in self.product['producer']['files'].items():
            self.assertEqual(gw.digest(gw.ROOT / name), value, name)
        for name, value in self.product['producer']['products'].items():
            self.assertEqual(gw.digest(gw.ROOT / name), value, name)
        self.assertFalse(constants_changed(self.product['producer']['constants']))

    def test_the_waves_mix_the_upper_air_between_earths_mixing_and_the_moon_scaled(self):
        for row in self.product['lunar'].values():
            waves = row['waves_stretched_homopause_pa']
            self.assertLess(row['moon_scaled_homopause_pa'], row['earth_homopause_pa'])
            self.assertIsNotNone(waves)
            self.assertLess(waves, row['earth_homopause_pa'])


if __name__ == '__main__':
    unittest.main()
