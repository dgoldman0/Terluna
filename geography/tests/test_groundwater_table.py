"""The water table: the solver against the strip's exact solution and with seepage, the product's binding, its water
balance in every state, and the young states' smaller recharge."""
import hashlib
import json
import math
import os
import unittest

import numpy as np

from geography import groundwater_table as gt


def strip(ground):
    rows, cols = 3, 1440
    sea = np.ones((rows, cols), bool)
    sea[:, 100:300] = False
    return dict(lat=np.array([0.25, 0.0, -0.25]), lon=(np.arange(cols) + 0.5) * 0.25, rows=rows, cols=cols,
                ground=np.full((rows, cols), ground), sea=sea, sea_level=0.0, lake=np.zeros((rows, cols), bool),
                lake_level=np.full((rows, cols), np.nan))


CAP, EXPONENT = 1e-12, 3.2
K_DEEP = CAP * (1e9 / 1000.0) ** EXPONENT       # the cap reaches a million km down: transmissivity all but constant


class Solver(unittest.TestCase):
    def test_strip_between_two_seas_is_a_parabola(self):
        s = strip(1000.0)                        # ground well above the mound, so nothing seeps
        rech = 8e-7
        r = gt.solve(s, np.full((s['rows'], s['cols']), rech), K_DEEP, CAP, EXPONENT, h0=np.zeros((3, 1440)))
        t = float(gt.transmissivity(np.array(900.0), K_DEEP, CAP, EXPONENT))
        dx = gt.MOON_RADIUS * math.radians(0.25)
        x = (np.arange(100, 300) - 99) * dx
        exact = rech * x * (201 * dx - x) / (2.0 * t)
        self.assertLess(float(np.max(np.abs(r['h'][1, 100:300] - exact))) / exact.max(), 1e-3)
        self.assertAlmostEqual(r['into_fixed'].sum() / r['recharge'].sum(), 1.0, places=4)

    def test_ground_below_the_mound_seeps_and_the_balance_closes(self):
        s = strip(20.0)
        r = gt.solve(s, np.full((s['rows'], s['cols']), 8e-7), K_DEEP, CAP, EXPONENT)
        land = ~s['sea']
        self.assertLessEqual(float(r['h'][land].max()), 20.0 + 0.5)
        self.assertGreater(float(r['seep'].sum()), 0.0)
        self.assertAlmostEqual((r['seep'].sum() + r['into_fixed'].sum()) / r['recharge'].sum(), 1.0, places=3)


class Product(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.product = json.loads(gt.OUT.read_text())

    def test_binds_its_code_and_inputs(self):
        p = self.product
        self.assertEqual(p['schema'], gt.SCHEMA)
        for path, h in {**p['producer']['files'], **p['producer']['inputs']}.items():
            self.assertEqual(hashlib.sha256((gt.ROOT / path).read_bytes()).hexdigest(), h, path)
        self.assertEqual(set(p['producer']['ignored_inputs']), set(gt.IGNORED))

    def test_every_state_balances_and_converges(self):
        for case, v in self.product['cases'].items():
            for state, row in v['states'].items():
                out = row['seepage_km3_yr'] + row['into_seas_km3_yr'] + row['into_lakes_km3_yr']
                self.assertAlmostEqual(out / row['recharge_km3_yr'], 1.0, delta=0.01, msg=(case, state))
                self.assertLess(row['rms_imbalance'], 10 * gt.NEWTON_TOLERANCE, (case, state))
                shares = row['depth_shares']
                self.assertAlmostEqual(sum(shares.values()), 1.0, delta=1e-4)

    def test_the_young_moon_recharges_less(self):
        for v in self.product['cases'].values():
            mature = v['states']['mature']['recharge_km3_yr']
            for state in ('young_1000_yr', 'young_10000_yr'):
                self.assertLessEqual(v['states'][state]['recharge_km3_yr'], mature * (1 + 1e-9))

    @unittest.skipUnless(os.environ.get('TERLUNA_GROUNDWATER_REGENERATE')
                         and all((gt.ROOT / p).exists() for p in gt.IGNORED),
                         'set TERLUNA_GROUNDWATER_REGENERATE=1 with the grid inputs restored (several minutes)')
    def test_regenerates(self):
        from geography.groundwater_uptake import rnd
        product, _ = gt.results()
        self.assertEqual(json.loads(json.dumps(rnd(product))),
                         {k: v for k, v in self.product.items() if k != 'producer'})


if __name__ == '__main__':
    unittest.main()
