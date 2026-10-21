"""GCM boundary: seas and rain-fed lakes combined on a coarse grid (no geography files needed), the lakes surface named by
its drainage's content, and its committed record when present."""
import hashlib
import json
import unittest
import numpy as np

from climate.gcm import boundary, exoplasim_run as er


class LakeTests(unittest.TestCase):
    def test_a_lakes_surface_takes_its_drainages_name_and_never_the_superseded_one(self):
        name = boundary.lakes_name('98842fc2a648f889', 'gaussian_T21')
        self.assertEqual(name, 'moon_28pct_water_lakes_98842fc2a648f889_gaussian_T21.nc')
        self.assertNotEqual(name, boundary.SUPERSEDED['files'].replace('<grid>', 'gaussian_T21'))
        # The runner's A28 keeps the surface of 25 September; A28R reads the one built from today's routing.
        self.assertEqual(er.surface_product(28, True), 'moon_28pct_water_lakes_gaussian_T21.nc')
        self.assertEqual(er.surface_product(28, er.EXPERIMENTS['A28R']['lakes']), name)
        self.assertEqual(er.surface_product(25, False), 'moon_25pct_water_gaussian_T21.nc')
        self.assertIn('060219b12f480619', boundary.SUPERSEDED['source_sha256'])

    def test_the_record_binds_its_files_and_counts_the_cells_that_change(self):
        path = boundary.RESULTS / f"boundary_lakes_{er.EXPERIMENTS['A28R']['lakes']}.json"
        if not path.exists():
            self.skipTest('the lakes surface has not been built')
        record = json.loads(path.read_text())
        self.assertEqual(record['schema'], boundary.LAKES_SCHEMA)
        drainage = record['inputs']['geography/products/drainage_28pct_16ppd.npz']
        self.assertTrue(drainage.startswith(er.EXPERIMENTS['A28R']['lakes']))
        shares, t21 = record['shares'], record['t21']
        self.assertAlmostEqual(shares['seas'] + shares['lakes'], shares['water'], places=12)
        self.assertAlmostEqual(t21['water_share'], shares['water'], places=6)          # the regridding keeps the water
        self.assertEqual(t21['cells'], 32 * 64)
        for name, sha in record['files'].items():
            local = boundary.HERE / 'products' / name
            if local.exists():
                self.assertEqual(hashlib.sha256(local.read_bytes()).hexdigest(), sha)
        if 'changed_cells' in record:
            changed = record['changed_cells']
            self.assertEqual(changed['count'], len(changed['cells']))
            self.assertEqual(changed['land_to_water'] + changed['water_to_land'], changed['count'])
            self.assertEqual(sum(not c['land'] for c in changed['cells']), changed['land_to_water'])
    def test_lakes_become_fractions_with_their_own_surface(self):
        height = np.array([[-2000.0, 500.0], [800.0, 1200.0]])          # coarse cells, 2 x 2 fine cells each
        sea = np.array([[True, False], [False, False]])
        fine_cols = 4
        # a lake over three fine cells of the top-right coarse cell, level 700 m; one over a sea cell is dropped
        cells = np.array([0 * fine_cols + 2, 0 * fine_cols + 3, 1 * fine_cols + 2, 0])
        levels = np.array([700.0, 700.0, 700.0, 900.0])
        wet, depth, surface = boundary.combine_lakes(height, sea, cells, levels, fine_cols, 2, -1654.0)
        np.testing.assert_allclose(wet, [[1.0, 0.75], [0.0, 0.0]])
        np.testing.assert_allclose(surface, [[-1654.0, 700.0], [-1654.0, -1654.0]])
        np.testing.assert_allclose(depth, [[346.0, 200.0], [0.0, 0.0]])


if __name__ == '__main__':
    unittest.main()
