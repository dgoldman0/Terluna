"""Annulus films: construction against the stored stack, mass bookkeeping and the product's findings."""
import hashlib
import json
import unittest
from pathlib import Path

import numpy as np

from protection.spectra import annulus_film as af
from protection.spectra import short_wave as sw

PRODUCT_PATH = Path(af.HERE) / 'annulus_film.json'
INPUTS_PRESENT = all((sw.SOURCES / name).is_file() for name in af.INPUTS)
NEEDS_INPUTS = unittest.skipUnless(INPUTS_PRESENT, 'protection inputs not restored (fetch_inputs.py --download)')


class MassTests(unittest.TestCase):
    def test_areal_mass_follows_the_design_densities(self):
        self.assertAlmostEqual(af.areal_mass(1.0, 10.0, False), 3.9 + 22.0)
        # The stored anti-reflection and matching layers add about a gram per square metre.
        extra = af.areal_mass(0.1, 2.0, True) - af.areal_mass(0.1, 2.0, False)
        self.assertTrue(0.5 < extra < 1.5)


@NEEDS_INPUTS
class FilmTests(unittest.TestCase):
    def test_the_thickest_coated_film_is_the_stored_stack_below_210_nm(self):
        w = np.linspace(20.0, 205.0, 400)
        silica, titania = af.indices(w)
        _, t = af.film(w, silica, titania, 1.0, 10.0, True)
        stored = sw.design_transmission(w, silica, titania, sw.stored_ar_layers())
        np.testing.assert_allclose(t, stored, rtol=1e-12, atol=1e-300)

    def test_vanishing_layers_pass_all_the_light(self):
        w = np.linspace(250.0, 2000.0, 50)
        silica, titania = af.indices(w)
        r, t = af.film(w, silica, titania, 0.0, 0.0, False)
        np.testing.assert_allclose(t, 1.0, atol=1e-12)
        np.testing.assert_allclose(r, 0.0, atol=1e-12)


@unittest.skipUnless(PRODUCT_PATH.is_file(), 'annulus_film.json not written')
class ProductTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.product = json.loads(PRODUCT_PATH.read_text())
        cls.films = {(f['titania_um'], f['silica_um'], f['coated']): f for f in cls.product['films']}

    def test_product_is_current_with_its_code(self):
        self.assertEqual(self.product['schema'], af.SCHEMA)
        from atmosphere.middle_atmosphere import escape
        paths = {'annulus_film.py': af.HERE / 'annulus_film.py', 'short_wave.py': af.HERE / 'short_wave.py',
                 'model.py': af.PROTECTION / 'model.py', 'escape.py': Path(escape.__file__)}
        for name, path in paths.items():
            self.assertEqual(self.product['producer']['files'][name],
                             hashlib.sha256(path.read_bytes()).hexdigest()[:16], name)

    def test_coating_passes_the_visible_that_bare_titania_reflects(self):
        bare, coated = self.films[(0.1, 2.0, False)], self.films[(0.1, 2.0, True)]
        self.assertGreater(coated['band_transmission']['visible_400_700'], 0.95)
        self.assertLess(bare['band_transmission']['visible_400_700'], 0.9)
        self.assertLess(coated['pressure_coefficient'], 0.5 * bare['pressure_coefficient'])

    def test_a_tenth_of_a_micrometre_of_titania_stops_the_far_and_extreme_ultraviolet(self):
        film = self.films[(0.1, 1.0, True)]
        self.assertLess(film['uv_transmission']['max_10_to_175'], 1e-6)
        self.assertLess(film['areal_mass_g_m2'], 4.0)

    def test_soft_x_rays_set_the_silica_and_hard_x_rays_pass_every_light_film(self):
        rows = [self.films[(0.1, s, True)] for s in (1.0, 2.0, 4.0, 10.0)]
        soft = [r['uv_transmission']['from_2_5_to_175'] for r in rows]
        self.assertTrue(all(b < a for a, b in zip(soft, soft[1:])))
        window = self.product['window_stack']
        for r in rows[:3]:
            self.assertGreater(r['uv_transmission']['below_175'], 10 * window['uv_transmission']['below_175'])
            self.assertGreater(r['transmitted_W_m2']['xray_0.1_2.5'], r['transmitted_W_m2']['fuv_122.5_175'])


if __name__ == '__main__':
    unittest.main()
