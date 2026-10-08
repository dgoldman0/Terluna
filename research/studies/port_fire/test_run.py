"""Port fire study: its helpers, and what the stored results must show."""
import json
import unittest

from research.studies.port_fire import run as pf
from shared.provenance import constants_changed


def stored():
    path = pf.HERE / 'results' / 'port_fire.json'
    return json.loads(path.read_text()) if path.is_file() else None


class Helpers(unittest.TestCase):
    def test_door_limit(self):
        self.assertAlmostEqual(pf.door_limit_pa(), 85.4, delta=0.1)

    def test_width_profile(self):
        port = json.loads(pf.PORT.read_text())['port']
        self.assertAlmostEqual(pf.width_at(port, 0.0), port['design']['base_width_m'], places=6)
        self.assertAlmostEqual(pf.width_at(port, port['design']['height_km'] * 1e3), pf.tower.FRAME.top_width_m, places=6)

    def test_air_cases(self):
        cases = pf.air_cases(json.loads(pf.PORT.read_text())['port'])
        self.assertEqual(len(cases), 5)
        for key, (_, air) in cases.items():
            if key != 'earth':
                self.assertAlmostEqual(air.oxygen_mole_fraction, 0.1746, places=4)
                self.assertLess(air.gravity_m_s2, 1.7)


class StoredResults(unittest.TestCase):
    def setUp(self):
        self.r = stored()
        if self.r is None:
            self.skipTest('results have not been generated')

    def test_schema_and_air(self):
        self.assertEqual(self.r['schema'], pf.SCHEMA)
        for name, h in {**self.r['producer']['files'], **self.r['producer']['inputs']}.items():
            self.assertEqual(pf.digest(pf.ROOT / name), h, name)
        self.assertFalse(constants_changed(self.r['producer']['constants']))
        for key, a in self.r['air'].items():
            if key == 'earth':
                continue
            self.assertAlmostEqual(a['time_factor'], 2.46, places=2)
            self.assertGreater(a['acts_as_earth_fire_times'], 2.4)
            self.assertGreater(a['flame_height_vs_earth'], 1.4)
            self.assertLess(a['room_flashover_vs_earth'], 1.0)
            self.assertLess(a['ventilation_limit_vs_earth'], 0.4)

    def test_sprinklers_answer_to_smaller_fires_and_smoke_fills_slower(self):
        for space in self.r['spaces'].values():
            for rows in space['fires'].values():
                earth = rows['earth']['earth_fire']
                for key, row in rows.items():
                    if key == 'earth':
                        continue
                    moon = row['earth_fire']
                    self.assertLess(moon['fire_at_sprinkler_kw'], earth['fire_at_sprinkler_kw'])
                    if moon['untenable_sprinklered_s'] is not None:
                        self.assertGreater(moon['untenable_sprinklered_s'], earth['untenable_sprinklered_s'])
                    self.assertGreater(row['lunar_estimate']['smoke_detector_s'], moon['smoke_detector_s'])

    def test_rooms_flash_over_smaller_and_burn_longer(self):
        earth = self.r['room_fire']['earth']['earth_fire']
        for key, row in self.r['room_fire'].items():
            if key != 'earth':
                self.assertLess(row['earth_fire']['flashover_kw'], earth['flashover_kw'])
                self.assertGreater(row['earth_fire']['burning_min'], 2 * earth['burning_min'])

    def test_sprinklered_compartments_leave_people_time(self):
        for use, e in self.r['egress'].items():
            for key in ('earth', 'enclosed_3_10km'):
                margin = e[key]['margin_sprinklered_s']
                self.assertTrue(margin is None or margin > 0, (use, key))

    def test_materials_margin(self):
        m = self.r['materials']
        self.assertGreaterEqual(m['test_oxygen_percent'], m['design_oxygen_percent'] + m['margin_points'])

    def test_systems(self):
        s = self.r['systems']
        z = s['water']['pressure_zone_m']['standard_175_psi']
        self.assertGreater(z['moon'], 5 * z['earth'])
        earth_run = s['shafts']['earth_comparison']['unbroken_shaft_within_door_window_m']
        for zone in s['shafts']['zones'].values():
            self.assertGreater(zone['unbroken_shaft_within_door_window_m'], 10 * earth_run)
        self.assertLess(s['firefighting']['arrival_s'], s['firefighting']['nfpa1710_first_engine_travel_s'])
        for use in s['compartment_burnout'].values():
            earth = use['by_air']['earth']['still_air']
            for key, row in use['by_air'].items():
                if key != 'earth':
                    self.assertGreater(row['still_air']['burning_h'], earth['burning_h'])
                    self.assertGreater(row['p90_wind']['fire_mw'], earth['fire_mw'])


if __name__ == '__main__':
    unittest.main()
