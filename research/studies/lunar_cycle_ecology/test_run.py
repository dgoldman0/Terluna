"""Ecology study: flight scaling, the air column, nutrient arithmetic and the stored results."""
import json
import math
import unittest

from research.studies.lunar_cycle_ecology import run as e
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY


class EcologyTests(unittest.TestCase):
    def test_flight_scales_with_weight_and_air_density(self):
        f = e.flight()
        g = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
        rho = f['air_density_kg_m3']
        self.assertTrue(1.40 < rho['moon'] < 1.46 and abs(rho['earth'] - 1.225) < 0.005)
        self.assertAlmostEqual(f['hover_induced_power_ratio'], g ** 1.5 * math.sqrt(rho['earth'] / rho['moon']))
        for flier in f['fliers'].values():
            per_km = flier['sugar_g_per_km']
            self.assertAlmostEqual(per_km['moon'] / per_km['earth'], g)

    def test_air_column_holds_the_design_oxygen(self):
        a = e.air(e.DESIGN_PRESSURE_PA, e.MOON_AIR_K, MOON_SURFACE_GRAVITY)
        self.assertAlmostEqual(a['column_kg_m2'], e.DESIGN_PRESSURE_PA / MOON_SURFACE_GRAVITY)
        self.assertAlmostEqual(a['o2_mole_fraction'], 21227.0 / e.DESIGN_PRESSURE_PA, places=6)
        self.assertTrue(0.19 < a['o2_kg_m2'] / a['column_kg_m2'] < 0.20)

    def test_c4_crossover_moves_warmer_with_more_co2(self):
        c4 = e.c4_crossover_shift()
        self.assertTrue(3.0 < c4['shift_k'] < 4.2)
        self.assertAlmostEqual(c4['moon_co2_pa'] / c4['earth_co2_pa'], 121590.0 / 101325.0)

    def test_rock_stock_arithmetic(self):
        k_per_pct = e.SOIL_BULK_KG_M3 * e.SOIL_DEPTH_M / 100 * e.K_IN_K2O * 1e3      # g K per m2 per % K2O
        self.assertAlmostEqual(e.K_IN_K2O, 0.830, places=3)
        self.assertAlmostEqual(e.P_IN_P2O5, 0.436, places=3)
        self.assertAlmostEqual(0.1 * k_per_pct, 1245.2, places=0)                        # a metre of 0.1% K2O soil

    def test_stored_results(self):
        path = e.HERE / 'results' / 'lunar_cycle_ecology.json'
        if not path.is_file():
            self.skipTest('results have not been generated')
        r = json.loads(path.read_text())
        self.assertEqual(r['schema'], e.SCHEMA)
        d = r['distances']
        self.assertLess(d['any_water']['p50'], d['sea']['p50'])
        within = list(d['sea']['land_share_within'].values())
        self.assertEqual(within, sorted(within))
        n = r['nutrients']
        self.assertAlmostEqual(n['runoff_m_per_yr_over_land'] * n['land_area_m2'] / (365.25 * 86400), 279167, delta=300)
        export = n['harvest_export']['earth_like_0.6']
        self.assertAlmostEqual(export['potassium_g_m2_yr'], export['fruit_kg_m2_yr'] * 1.12, delta=0.1)
        self.assertGreater(export['potassium_g_m2_yr'], 20 * max(n['earth_river_load_at_moon_runoff_g_m2_yr']['potassium']))
        mare = n['weathering_supply']['mare_soil']
        self.assertAlmostEqual(mare['rock_weathered_g_m2_yr'][1] / mare['rock_weathered_g_m2_yr'][0], 10.0, places=3)
        k_share = mare['potassium_g_m2_yr'][0] / mare['rock_weathered_g_m2_yr'][0]      # results keep 4 digits
        self.assertAlmostEqual(k_share / (e.ROCK_WT_PCT['mare_soil']['K2O'] / 100 * e.K_IN_K2O), 1.0, places=3)
        leaves = r['leaves']
        self.assertGreater(leaves['plant_keeps_g_c_per_cycle_equator']['earth_like']['0.6'],
                           leaves['leaf_replacement_g_c_per_cycle']['6_months'])
        self.assertLess(r['night']['night_draw_equator']['0.6']['o2_share_of_column'], 1e-4)


if __name__ == '__main__':
    unittest.main()
