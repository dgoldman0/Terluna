"""Summit tower study: the standard atmosphere, distances, the design wind and the stored results."""
import json
import unittest

from research.studies.summit_tower import run as st


class SummitTowerTests(unittest.TestCase):
    def test_standard_atmosphere(self):
        self.assertAlmostEqual(float(st.earth_air(0.0)), 1.225, places=3)
        self.assertAlmostEqual(float(st.earth_air(11e3)), 0.3639, places=3)
        self.assertAlmostEqual(float(st.earth_air(20e3)), 0.0880, places=3)

    def test_great_circle(self):
        quarter = float(st.great_circle_km(0.0, 0.0, 0.0, 90.0))
        self.assertAlmostEqual(quarter, st.MOON.radius_m / 1e3 * 3.141592653589793 / 2, places=6)

    def test_design_wind_factors(self):
        climate = [dict(height_above_ground_km=0.0, ring_3hourly=dict(max=10.0), gcm_3day=dict(max=6.0)),
                   dict(height_above_ground_km=5.0, ring_3hourly=dict(max=12.0), gcm_3day=dict(max=7.0)),
                   dict(height_above_ground_km=20.0, ring_3hourly=dict(max=30.0), gcm_3day=dict(max=9.0))]
        w = st.design_winds(climate, 10.0)
        self.assertEqual(w['models_highest_m_s'], 12.0)
        self.assertAlmostEqual(w['model_based_design_gust_m_s'], round(12.0 * 1.4 * 1.2 * 1.2, 1))

    def test_stored_results(self):
        path = st.HERE / 'results' / 'summit_tower.json'
        if not path.is_file():
            self.skipTest('results have not been generated')
        r = json.loads(path.read_text())
        self.assertEqual(r['schema'], st.SCHEMA)
        self.assertTrue(11000 < r['site']['ground_above_sea_m'] < 13000)
        gusts = [f'{v:g}' for v in r['design_winds']['scenarios_m_s']]
        for height, row in r['moon'].items():
            masses = [row[g]['frame_mass_mt'] for g in gusts if row[g] is not None]
            self.assertEqual(masses, sorted(masses), height)          # more wind, more steel
        for height, row in r['earth'].items():
            moon = r['moon'].get(height)
            if moon and moon['50'] and row['50']:
                self.assertLess(moon['50']['frame_mass_mt'], row['50']['frame_mass_mt'])
        ref = r['reference']['summary']
        self.assertLessEqual(ref['sway_at_service_m'], ref['height_km'] * 1e3 / 500)
        self.assertGreaterEqual(ref['buckling_factor'], 3.0)
        e = r['energy']
        for dev in e['devices'].values():
            self.assertLess(dev['mean_mw'], e['resource_mw'])
        self.assertLess(e['devices']['bladeless']['mean_mw'], e['devices']['slow_rotor']['mean_mw'])
        loads = r['loads']
        self.assertGreater(loads['closed_tower']['base_moment_ratio_same_outline'], 3.0)
        for key in ('slow_rotor_stopped', 'bladeless_stopped'):
            self.assertGreater(loads[key]['lightest_frame']['frame_mass_mt'], loads['open_frame']['lightest_frame']['frame_mass_mt'])
        p = r['collisions']['probability_per_crossing']
        self.assertGreater(p['slow_rotor']['moth_0.33g'], p['fast_rotor']['moth_0.33g'])
        self.assertGreater(r['night_storage']['water_mm3']['city'], r['night_storage']['water_mm3']['tower'])


if __name__ == '__main__':
    unittest.main()
