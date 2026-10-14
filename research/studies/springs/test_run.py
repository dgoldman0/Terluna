"""The springs study: its chemistry against the carbon cycle's solver, its transit times against Vogel's exponential
case, its mixing's conservation, and the products' binding and agreement with stage 2."""
import hashlib
import json
import math
import os
import unittest

import numpy as np

from geography import groundwater_table as gt
from research.studies.carbon_cycle import run as cc
from research.studies.springs import carbon as sc
from research.studies.springs import run as sr

SHARES = dict(ca=0.65, mg=0.33, mono=0.02)


def digest(path):
    return hashlib.sha256((sr.ROOT / path).read_bytes()).hexdigest()


class Chemistry(unittest.TestCase):
    def test_constants_and_speciation_match_the_carbon_cycle(self):
        for t in (285.0, 294.3, 300.0):
            a, b = cc.fresh_constants(t), sr.constants(np.array([t]))
            for key in ('kh', 'k1', 'k2', 'ksp', 'kw', 'davies'):
                self.assertAlmostEqual(float(b[key][0]) / a[key], 1.0, places=12)
        box = dict(mass=1.0, ta=2.0e-3, ca=0.65e-3, mg=0.33e-3, mono=0.02e-3, dic=0.0)
        ref = cc.speciate(box, 4.8e-4, cc.fresh_constants(294.3))
        w = {key: np.array([box[key]]) for key in sr.SPECIES}
        out, _, ph = sr.open_water(w, 4.8e-4, sr.constants(np.array([294.3])), 1e9)
        self.assertAlmostEqual(float(ph[0]), ref['ph'], places=8)
        self.assertAlmostEqual(float(out['dic'][0]) / ref['dic'], 1.0, places=8)

    def test_calcite_comes_out_to_saturation_and_the_rock_stops_at_its_ph(self):
        c = sr.constants(np.full(3, 294.0))
        soil, _, _ = sr.open_water(sr.add(sr.zero(3), np.array([6e-3, 1.5e-3, 4e-3]), SHARES),
                                   np.array([0.02, 0.02, 0.005]), c, 1.0)
        omega = sr.closed_water(dict(soil), c, 1e9)[4]
        self.assertTrue(np.all(omega < 1.0 + 1e-3))
        sealed = sr.closed_water(sr.add(soil, np.full(3, 3e-3), SHARES), c, 1.0)
        self.assertTrue(np.all(np.abs(sealed[4] - 1.0) < 1e-3))
        room = sr.reaction_room(soil, c, SHARES, sr.PH_MAX, 1.0)
        ph = sr.closed_water(sr.add(soil, room, SHARES), c, 1.0)[2]
        self.assertTrue(np.all(np.abs(ph - sr.PH_MAX) < 0.05))


class Transit(unittest.TestCase):
    def table(self, k, phi, bottom):
        z = np.arange(0.0, bottom + 1.0, 1.0)
        kk, pp = np.full(z.size, k), np.full(z.size, phi)
        seg = 0.5 * (kk[1:] + kk[:-1]) * np.diff(z)
        ki = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
        mid = 0.5 * (ki[1:] + ki[:-1])
        big = np.concatenate([[0.0], np.cumsum(phi / np.maximum(mid, 1e-300) * np.diff(z))])
        pores = np.concatenate([[0.0], np.cumsum(phi * np.diff(z))])
        return dict(z=z, k=kk, phi=pp, ki=ki, big=big, pores=pores)

    def test_constant_properties_give_vogels_exponential_ages(self):
        k, phi, h, rech = 1e-5, 0.2, 500.0, 0.25           # m/s, porosity, saturated thickness, m/yr
        tab = self.table(k, phi, h)
        depth, r = np.array([0.0]), np.array([rech])
        _, ages, mean = sr.saturated_spread(tab, depth, r)
        tau = phi * h / rech
        self.assertAlmostEqual(float(mean[0]) / tau, 1.0, places=6)
        nodes = 0.5 * (sr.P_EDGES[:-1] + sr.P_EDGES[1:])
        exact = -tau * np.log(nodes)
        # The grid resolves the ages to 0.2% down to the deepest hundredth of the flow; the last thousandth, which runs
        # along the bottom where the share of flow below falls to zero, to 1%.
        keep, tail = nodes > 1e-2, nodes > 1e-3
        self.assertLess(float(np.max(np.abs(ages[0, keep] / exact[keep] - 1.0))), 2e-3)
        self.assertLess(float(np.max(np.abs(ages[0, tail] / exact[tail] - 1.0))), 1e-2)
        years = np.array([100.0, 400.0, 1000.0])
        p = sr.older_than(tab, depth, r, np.zeros(1), years)[0]
        self.assertLess(float(np.max(np.abs(p - np.exp(-years / tau)))), 2e-3)

    def test_unsaturated_rock_holds_the_water_its_conductivity_needs(self):
        k, phi, rech, d = 1e-6, 0.15, 0.05, 200.0
        tab = self.table(k, phi, 2000.0)
        soil, rock = sr.unsaturated_years(tab, np.array([d]), np.array([rech]))
        theta = phi * (rech / (k * sr.YEAR_S)) ** (1.0 / (2.0 * sr.B_ROCK + 3.0))
        self.assertAlmostEqual(float(rock[0]) / (theta * (d - sr.SOIL_DEPTH_M) / rech), 1.0, places=4)


class Mixing(unittest.TestCase):
    def test_the_mixing_conserves_water_and_what_it_carries(self):
        rows, cols = 3, 1440
        sea = np.ones((rows, cols), bool)
        sea[:, 100:300] = False
        s = dict(lat=np.array([0.25, 0.0, -0.25]), lon=(np.arange(cols) + 0.5) * 0.25, rows=rows, cols=cols,
                 ground=np.full((rows, cols), 20.0), sea=sea, sea_level=0.0, lake=np.zeros((rows, cols), bool),
                 lake_level=np.full((rows, cols), np.nan))
        k_deep = 1e-12 * (1e9 / 1000.0) ** 3.2
        r = gt.solve(s, np.full((rows, cols), 8e-7), k_deep, 1e-12, 3.2)
        solve, fx = sr.mixer(s, r, k_deep)
        n = rows * cols
        prop = np.zeros((n, 2))
        prop[:, 0] = 1.0
        prop[:, 1] = np.arange(n) % 7
        out = solve(prop)
        land = ~sea.ravel()
        seep = r['seep'].ravel()
        self.assertTrue(np.all(np.abs(out[land & (seep > 0), 0] - 1.0) < 1e-9))
        carried = (seep[land] * out[land, 1]).sum() + (fx['q'] * out[fx['src'], 1]).sum()
        given = ((fx['rech'] + fx['gift']) * prop[:, 1])[land].sum()
        self.assertAlmostEqual(carried / given, 1.0, places=3)


class Products(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.product = json.loads(sr.OUT.read_text())
        cls.carbon = json.loads(sc.OUT.read_text())

    def test_bind_their_code_and_inputs(self):
        p = self.product
        self.assertEqual(p['schema'], sr.SCHEMA)
        for path, h in {**p['producer']['files'], **p['producer']['inputs']}.items():
            self.assertEqual(digest(path), h, path)
        self.assertEqual(set(p['producer']['ignored_inputs']), set(sr.IGNORED))
        c = self.carbon
        self.assertEqual(c['schema'], sc.SCHEMA)
        for path, h in {**c['producer']['files'], **c['producer']['inputs']}.items():
            self.assertEqual(digest(path), h, path)

    def test_it_re_solves_stage_2(self):
        table = json.loads((sr.ROOT / 'geography/results/groundwater_table.json').read_text())
        for case, v in self.product['cases'].items():
            for state, row in v['states'].items():
                moon = row['variants']['central']['moon']
                ref = table['cases'][case]['states'][state]
                self.assertAlmostEqual(moon['seepage_km3_yr'] / ref['seepage_km3_yr'], 1.0, delta=0.01,
                                       msg=(case, state))
                self.assertAlmostEqual(moon['recharge_km3_yr'] / ref['recharge_km3_yr'], 1.0, delta=0.01)

    def test_the_springs_read_consistently(self):
        for v in self.product['cases'].values():
            for row in v['states'].values():
                for variant in row['variants'].values():
                    for block in [variant['moon'], *variant['by_band'].values()]:
                        if block.get('seepage_km3_yr', 0.0) <= 0.0:
                            continue
                        hard = (block['calcium_mmol_kg'] + block['magnesium_mmol_kg']) * sr.CACO3_G_MOL
                        self.assertAlmostEqual(block['hardness_mg_l_caco3'] / hard, 1.0, delta=0.01)
                        shares = list(block['younger_than'].values())
                        self.assertTrue(all(a <= b + 1e-9 for a, b in zip(shares, shares[1:])))
                        self.assertLessEqual(block['alkalinity_meq_kg'], sr.TA_CEILING * 1e3 * 1.001)
                        b = block['building_share_of_flow']
                        self.assertLessEqual(b['tufa'], b['designed'] + 1e-9)

    def test_the_filling_crust_adds_up(self):
        for rows in self.carbon['filling_crust'].values():
            for row in rows.values():
                c = row['carbon_gt_c_yr']
                self.assertAlmostEqual(c['seas'] + c['lakes'] + c['land'], c['total'], delta=1e-3 * c['total'] + 1e-9)

    @unittest.skipUnless(os.environ.get('TERLUNA_SPRINGS_REGENERATE')
                         and all((sr.ROOT / p).exists() for p in sr.IGNORED),
                         'set TERLUNA_SPRINGS_REGENERATE=1 with the grid inputs restored (about twenty minutes)')
    def test_regenerates(self):
        from geography.groundwater_uptake import rnd
        product, _ = sr.results()
        self.assertEqual(json.loads(json.dumps(rnd(product))),
                         {k: v for k, v in self.product.items() if k != 'producer'})


if __name__ == '__main__':
    unittest.main()
