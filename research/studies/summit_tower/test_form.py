"""Summit tower form: the hand-sized chosen form balances its loads and grows its members downward."""
import json
import unittest

from research.studies.summit_tower import form


class FormTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = form.results()

    def test_members_grow_downward(self):
        forces = [row['member_force_n'] for row in self.r['frame']['profile']]
        areas = [row['area_m2'] for row in self.r['frame']['profile']]
        self.assertEqual(forces, sorted(forces, reverse=True))   # profile runs from the base up
        self.assertGreater(areas[0], 10 * areas[-1])

    def test_base_carries_the_transfer_load(self):
        b = self.r['base']
        # six nodes on the legs and three on each of six arches: all 24 node loads
        self.assertAlmostEqual((6 + 6 * 3) * b['node_load_n'], b['weight_at_transfer_n'], delta=1.0)
        self.assertGreater(b['arch']['force_n'], b['arch']['thrust_n'])

    def test_rings_and_disks(self):
        for ring in self.r['rings']:
            self.assertAlmostEqual(ring['truss']['root_depth_m'], max(ring['width_m'] / 6, 8.0))
            self.assertLess(ring['width_m'], ring['outer_radius_m'])
        radii = [d['radius_m'] for d in self.r['disks']]
        self.assertEqual(radii, sorted(radii, reverse=True))
        for d in self.r['disks']:
            self.assertAlmostEqual(d['lens_depth_m'], form.LENS * d['radius_m'])
            self.assertAlmostEqual(d['hub_ring_area_m2'], d['horizontal_force_per_radian_n'] / form.SIGMA)

    def test_stored_results(self):
        path = form.HERE / 'results' / 'summit_tower_form.json'
        if not path.is_file():
            self.skipTest('results have not been generated')
        stored = json.loads(path.read_text())
        self.assertEqual(stored['schema'], form.SCHEMA)
        self.assertAlmostEqual(stored['steel_total_kg'], self.r['steel_total_kg'], delta=1e6)


if __name__ == '__main__':
    unittest.main()
