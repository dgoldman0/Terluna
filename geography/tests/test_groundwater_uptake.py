"""The crust's uptake: the solvers against analytic solutions, the product's binding and regeneration, and the order
of its results through time, across the cases and by latitude."""
import hashlib
import json
import math
import unittest

import numpy as np
from scipy.special import erfc

from geography import groundwater_uptake as gu
from shared.provenance import constants_changed

PRODUCT = json.loads(gu.OUT.read_text())


def digest(path):
    return hashlib.sha256((gu.ROOT / path).read_bytes()).hexdigest()


class Solvers(unittest.TestCase):
    def test_conduction_matches_the_error_function(self):
        """A step at the top of a linear profile with constant properties, against the analytic solution."""
        edges, z, dz = gu.grid()
        k, c, gradient = 1.5, 2.0e6, 0.02
        t0 = 250.0 + gradient * z
        t = t0.copy()
        years = np.concatenate([[0.0], 10 ** np.arange(-1.0, 4.0 + 1e-9, 1.0 / gu.STEPS_PER_DECADE)])
        for y0, y1 in zip(years[:-1], years[1:]):
            t = gu.heat_step(t, (y1 - y0) * gu.YEAR_S, np.full(z.size, k), np.full(z.size, c), dz, 290.0,
                             k * gradient)
        exact = t0 + 40.0 * erfc(z / (2.0 * math.sqrt(k / c * 1e4 * gu.YEAR_S)))
        upper = z < 3000.0
        self.assertLess(float(np.max(np.abs(t[upper] - exact[upper]))), 0.5)

    def test_infiltration_matches_green_and_ampt(self):
        """Constant conductivity and porosity: t = (phi/K)(L - H ln(1 + L/H)) ponded, t = phi L / K under gravity."""
        edges, z, dz = gu.grid()
        for cond, head, end in ((1e-10, 500.0, 1e3), (3e-10, 2000.0, 1e4), (1e-11, 0.0, 1e4)):
            years = np.concatenate([[0.0], 10 ** np.arange(-1.0, math.log10(end) + 1e-9, 1.0 / gu.STEPS_PER_DECADE)])
            front = 0.0
            for dt in np.diff(years) * gu.YEAR_S:
                front = gu.advance(max(front, 1e-6), dt, 2e4, edges, np.full(z.size, 0.1), np.full(z.size, cond),
                                   head, None)
            exact = (0.1 / cond) * (front - head * math.log(1.0 + front / head)) if head else 0.1 * front / cond
            self.assertAlmostEqual(exact / gu.YEAR_S / end, 1.0, delta=0.02)


class Product(unittest.TestCase):
    def test_binds_its_code_and_inputs(self):
        self.assertEqual(PRODUCT['schema'], gu.SCHEMA)
        for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
            self.assertEqual(digest(path), h, path)
        self.assertEqual(set(PRODUCT['producer']['ignored_inputs']), set(gu.IGNORED))
        for path, h in PRODUCT['producer']['ignored_inputs'].items():
            if (gu.ROOT / path).exists():
                self.assertEqual(digest(path), h, path)
        self.assertFalse(constants_changed(PRODUCT['producer']['constants']))

    @unittest.skipUnless(all((gu.ROOT / p).exists() for p in gu.IGNORED), 'the grid inputs are not restored')
    def test_regenerates(self):
        fresh = json.loads(json.dumps(gu.rnd(gu.results())))
        self.assertEqual(fresh, {k: v for k, v in PRODUCT.items() if k != 'producer'})

    def test_the_crust_fills_in_order(self):
        years = [f'{y:g}' for y in gu.REPORT_YEARS]
        pore_column = PRODUCT['setting']['porosity_phi0'] * PRODUCT['setting']['porosity_e_folding_m'] * (
            1.0 - math.exp(-gu.PORE_BOTTOM_M / PRODUCT['setting']['porosity_e_folding_m']))
        layers = {}
        for case in gu.CASES:
            held = PRODUCT['cases'][case]['moon']['held']
            series = [held[y]['global_layer_m'] for y in years]
            self.assertTrue(all(a <= b for a, b in zip(series, series[1:])), case)
            self.assertLessEqual(series[-1], pore_column * 1.0001)
            layers[case] = series
            for band in PRODUCT['cases'][case]['bands'].values():
                for col in band['columns'].values():
                    held_m = [col['by_year'][y]['held_m'] for y in years]
                    self.assertTrue(all(a <= b + 1e-9 for a, b in zip(held_m, held_m[1:])))
        for i in range(len(years)):
            self.assertLessEqual(layers['slow'][i], layers['central'][i] * 1.0001)
            self.assertLessEqual(layers['central'][i], layers['fast'][i] * 1.0001)
        warm = PRODUCT['warm_crust_central']['moon']['held']
        self.assertGreater(warm['1000']['global_layer_m'], 1.5 * layers['central'][years.index('1000')])

    def test_the_cold_lasts_longer_toward_the_poles(self):
        bands = PRODUCT['cases']['central']['bands']
        names = [b['band'] for b in PRODUCT['setting']['bands']]
        thawed = [bands[n]['columns']['land']['thawed_through_yr'] for n in names]
        frozen = [bands[n]['columns']['land']['airless_frozen_layer_m'][1] for n in names]
        self.assertTrue(all(a < b for a, b in zip(thawed, thawed[1:])))
        self.assertTrue(all(a < b for a, b in zip(frozen, frozen[1:])))
        self.assertGreater(frozen[0], 100.0)

    def test_the_seas_fall_spreads_the_held_water_over_them(self):
        area = 4.0 * math.pi * gu.MOON_RADIUS ** 2
        sea = sum(b['sea_m2'] for b in PRODUCT['setting']['bands'])
        for row in PRODUCT['cases']['central']['moon']['held'].values():
            self.assertAlmostEqual(row['seas_fall_m'], row['global_layer_m'] * area / sea,
                                   delta=1e-3 * row['seas_fall_m'] + 1e-6)
            self.assertAlmostEqual(row['share_of_seas'], row['kg'] / PRODUCT['setting']['seas_kg'], places=4)


if __name__ == '__main__':
    unittest.main()
