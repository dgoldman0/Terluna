"""Atomic oxygen and hydrogen made in the upper air, and what of them escapes (atmosphere/loss_response/oxygen.py)."""
from __future__ import annotations
import json
import math
import unittest

import numpy as np
from scipy.special import expn

from shared.provenance import constants_changed
from atmosphere.loss_response import oxygen as ox


class ClosedForms(unittest.TestCase):
    def test_the_drift_weights_hold_diffusive_equilibrium(self):
        # With no net flux the exponential fitting keeps f_{k+1} / f_k = exp(-P), the species' own scale height.
        for pe in (-30.0, -1.0, -1e-9, 0.0, 1e-9, 0.5, 40.0):
            up, down = ox.bernoulli(np.array([pe]))[0], ox.bernoulli(np.array([-pe]))[0]
            self.assertAlmostEqual(up / down, math.exp(-pe), delta=1e-12 * max(1.0, math.exp(-pe)))
        self.assertEqual(ox.bernoulli(np.array([0.0]))[0], 1.0)

    def test_the_glow_is_absorbed_where_the_flux_falls(self):
        # Isotropic light from above: what a thick column absorbs is the flux onto it, pi I.
        n_layers = 60
        p = 0.3 * np.exp(-np.linspace(-4.0, 12.0, n_layers + 1))
        col = dict(layer_p_pa=np.sqrt(p[1:] * p[:-1]), layer_t_k=np.full(n_layers, 180.0),
                   layer_x_h2o=np.zeros(n_layers), dry_fractions=dict(O2=0.175, N2=0.825),
                   layer_dz_cm=np.full(n_layers, 5e6))
        j_o2, j_h2o = ox.glow_rates(col, 'quiet')
        n = col['layer_p_pa'] / (1.380649e-23 * 180.0) * 1e-6
        absorbed = np.sum(j_o2 * 0.175 * n * col['layer_dz_cm'])
        intensity = ox.lr.GLOW_RAYLEIGH * 1e10 / (4 * math.pi) * 1e-4
        tau = ox.escape.O2_LYMAN_ALPHA_CM2 * np.sum(0.175 * n * col['layer_dz_cm'])
        self.assertAlmostEqual(absorbed / (math.pi * intensity * (1 - 2 * expn(3, tau))), 1.0, places=9)
        self.assertTrue(np.all(j_h2o > 0) and j_h2o[-1] > j_h2o[0])

    def test_hot_atoms_from_the_top_of_a_thick_column_follow_the_closed_form(self):
        # Atoms made at collision depth s escape on straight paths with chance E2(s) / 2; summed over a column the glow
        # crosses whole, its top well above the exobase, they are the share sigma x / (2 sigma_c) of the atoms made.
        n_layers, t = 400, 180.0
        lnp = np.linspace(-4.0, 14.0, n_layers + 1)
        p = 0.3 * np.exp(-lnp)
        scale_height_cm = 1.380649e-23 * t / (0.0289 / 6.02214076e23 * 1.62) * 100.0
        col = dict(layer_p_pa=np.sqrt(p[1:] * p[:-1]), layer_t_k=np.full(n_layers, t), layer_x_h2o=np.zeros(n_layers),
                   dry_fractions=dict(O2=0.175, N2=0.825),
                   layer_dz_cm=np.full(n_layers, scale_height_cm * (lnp[1] - lnp[0])))
        j_o2, _ = ox.glow_rates(col, 'quiet')
        n = col['layer_p_pa'] / (1.380649e-23 * t) * 1e-6
        made = 2 * j_o2 * 0.175 * n * col['layer_dz_cm']
        column = n * col['layer_dz_cm']
        above = np.cumsum(column[::-1])[::-1] - column / 2
        sigma_c = 1e-17
        straight = np.sum(made * expn(2, sigma_c * above) / 2)
        closed = ox.hot_escape(made.sum(), ox.escape.O2_LYMAN_ALPHA_CM2, 0.175, sigma_c)
        self.assertAlmostEqual(straight / closed, 1.0, delta=0.02)

    def test_hard_sphere_scaling_keeps_oxygens_own_parameter(self):
        self.assertAlmostEqual(ox.diffusion_parameter('O', 200.0), 9.69e16 * 200.0 ** 0.774, delta=1.0)
        self.assertGreater(ox.diffusion_parameter('O', 200.0), ox.diffusion_parameter('O3', 200.0))


@unittest.skipUnless(ox.OUT.is_file(), 'oxygen_escape.json not written')
class ProductTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.product = json.loads(ox.OUT.read_text())

    def test_product_is_current_with_its_code_and_products(self):
        self.assertEqual(self.product['schema'], ox.SCHEMA)
        for name, value in self.product['producer']['files'].items():
            self.assertEqual(ox.digest(ox.ROOT / name), value, name)
        for name, value in self.product['producer']['products'].items():
            self.assertEqual(ox.digest(ox.ROOT / name), value, name)
        for name, value in self.product['producer']['inputs'].items():
            self.assertEqual(ox.digest(ox.ROOT / name), value, name)
        self.assertFalse(constants_changed(self.product['producer']['constants']))

    def test_the_atoms_escape_no_faster_than_they_are_made_and_stay_a_minor_gas(self):
        for run in self.product['runs']:
            r = run['result']
            if r['status'] != 'solved':
                continue
            made = r['oxygen_made_kg_s']['glow'] + r['oxygen_made_kg_s']['traced']
            self.assertLess(r['oxygen_loss_kg_s']['no_magnetosphere'], made)
            self.assertLessEqual(r['oxygen_loss_kg_s']['september_magnetosphere'],
                                 r['oxygen_loss_kg_s']['no_magnetosphere'])
            # The column stays molecular: with the Moon's mixing the atoms are a few percent at most at the exobase;
            # Earth's weaker mixing, the bound, lets them reach about a tenth there.
            self.assertLess(r['oxygen_share_at_exobase'], 0.05 if run['setup'] == 'moon_mixing' else 0.15)

    def test_weaker_mixing_lets_more_oxygen_reach_the_exobase(self):
        runs = {(r['setup'], r['treatment'], r['activity']): r['result'] for r in self.product['runs']
                if r['shield'] == 'titania_stack' and r['result']['status'] == 'solved'}
        for (setup, treatment, activity), r in runs.items():
            if setup != 'earth_mixing' or ('moon_mixing', treatment, activity) not in runs:
                continue
            moon = runs[('moon_mixing', treatment, activity)]
            self.assertGreater(r['oxygen_share_at_exobase'], moon['oxygen_share_at_exobase'])
            titan = runs.get(('titan_analogue', treatment, activity))
            if titan:
                self.assertTrue(moon['oxygen_share_at_exobase'] < titan['oxygen_share_at_exobase']
                                < r['oxygen_share_at_exobase'], (treatment, activity))


if __name__ == '__main__':
    unittest.main()
