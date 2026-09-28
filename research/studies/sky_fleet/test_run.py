"""The sky-fleet study reproduces the sky-ship study's long-haul ships, fits its reference gliders, lays out berths that
keep their swing, counts calls and flyers consistently, and stores a product whose hashes match its producer."""
from __future__ import annotations
import json
import unittest

import numpy as np

from engineering.flight import winged as wi
from research.studies.sky_fleet import run
from research.studies.sky_ships import run as ships


class Ships(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = run.load()
        cls.factors = cls.p['ships']['buoyant']['factors']
        cls.band = run.air(cls.p['winds'], ships.BAND_KM)
        cls.top = run.air(cls.p['winds'], ships.FLIGHT_BAND_KM[1])

    def test_long_haul_ships_are_the_sky_ship_studys(self):
        """At 30 m/s over the 48-hour route with 400 kg a passenger, the ship is the sky-ship study's."""
        gust = self.p['ships']['band']['updraft_max_m_s']
        stored = {r['length_m']: r for r in self.p['ships']['operations']['ships']}
        for L in (500.0, 1000.0):
            r = run.ship(self.band, self.top, self.factors, L, ships.CRUISE_M_S, run.LONG_HAUL_ROUTE_KM,
                         ships.PASSENGER_KG[1], gust)
            self.assertAlmostEqual(r['passengers'] / stored[L]['passengers'][1], 1.0, delta=0.01)
            self.assertAlmostEqual(r['hydrogen_t'], stored[L]['hydrogen_t'], delta=0.2)

    def test_faster_ships_carry_fewer_for_more_energy(self):
        gust = self.p['ships']['band']['updraft_max_m_s']
        slow, fast = (run.ship(self.band, self.top, self.factors, 750.0, v, run.LONG_HAUL_ROUTE_KM, 400.0, gust)
                      for v in (30.0, 50.0))
        self.assertLess(fast['passengers'], slow['passengers'])
        self.assertGreater(fast['wh_per_passenger_km'], slow['wh_per_passenger_km'])

    def test_berths_hold_the_class_and_its_swing(self):
        row = dict(length_m=500.0, diameter_m=83.3, fin_span_m=20.9)
        across = row['diameter_m'] + 2 * row['fin_span_m']
        for n in (1, 4, 27):
            b = run.crown_berths(n, row, 351.0, 23481.0, 24000.0)
            self.assertGreaterEqual(b['levels'] * 2 * b['ships_per_arm'] if n > 1 else 1, n)
            self.assertGreaterEqual(b['levels_above_summit_m'][0] - 0.5 * across, 23481.0)   # clears the terminal
            self.assertLessEqual(b['levels_above_summit_m'][-1], 24000.0 + 1.0)               # arms on the frame
            self.assertAlmostEqual(np.degrees(np.arccos(across / b['pitch_m'])), run.SWING_DEG, delta=0.3)

    def test_network_calls_carry_the_traffic(self):
        n = run.network([40000.0], 2000.0, 2729.0, 30.0, 2.0)
        self.assertAlmostEqual(n['calls_per_hour'][0] * 2 * 2000.0, 40000.0, delta=1.0)
        self.assertAlmostEqual(n['berths'][0], 20.0, places=6)
        self.assertAlmostEqual(n['destinations_daily'][0], 240.0, places=6)


class People(unittest.TestCase):
    def test_fitted_polar_holds_its_glide_and_speed(self):
        pol = run.fit_polar(10.0, 14.0, 110.0, 12.0, 12.0, 1.3)
        flat = wi.best_glide(pol, 110.0, run.G_EARTH, 1.225, stall_margin=0.0)
        self.assertAlmostEqual(float(flat['lift_to_drag']), 12.0, places=6)
        self.assertAlmostEqual(float(flat['speed_m_s']), 12.0, places=6)

    def test_soaring_hours_add_up(self):
        ring = json.loads(run.RING.read_text())
        s = run.soaring(ring)
        day = run.SYNODIC_MONTH_DAYS * 24.0
        self.assertLessEqual(s['hours_dry_and_over_2_km'], s['hours_mixed_layer_over_2_km'])
        self.assertLess(s['hours_mixed_layer_over_2_km'], day)
        self.assertGreater(s['mixed_layer_max_km'], s['night_mixed_layer_km'])


class City(unittest.TestCase):
    def test_flyers_aloft_scale_with_their_share(self):
        modes = dict(a=dict(share=0.1, km=10.0, occupancy=1.0, flyer='x'), b=dict(share=0.2, km=10.0, occupancy=1.0,
                                                                               flyer='x'))
        per = dict(a=dict(minutes=12.0, kwh=1.0), b=dict(minutes=12.0, kwh=1.0))
        city = dict(population=1e6, area_km2=100.0, trips_per_person=(3.0, 3.0, 3.0), busy_hour_share=(0.1, 0.1, 0.1))
        m = run.metropolis(city, modes, per, dict(passengers=1000.0))['base']
        self.assertAlmostEqual(m['b']['aloft_busy_hour'] / m['a']['aloft_busy_hour'], 2.0, places=2)
        self.assertAlmostEqual(m['a']['aloft_busy_hour'], 3e6 * 0.1 * 0.1 * 12.0 / 60.0, delta=1.0)

    def test_nearest_neighbour_in_a_layer(self):
        bands = dict(x=dict(bottom_m=0.0, top_m=1000.0, layers=1, speed_m_s=10.0, modes=dict(a=1.0)))
        a = run.airspace(dict(a=100.0), bands, 1.0)
        self.assertAlmostEqual(a['x']['nearest_in_layer_m'], 50.0, places=6)     # 0.5 / sqrt(100 per km2)


class Product(unittest.TestCase):
    def test_stored_product(self):
        r = json.loads((run.HERE / 'results' / 'sky_fleet.json').read_text())
        self.assertEqual(r['schema'], run.SCHEMA)
        for name, h in r['producer']['files'].items():
            self.assertEqual(run.digest(run.ROOT / name), h, name)
        for name, h in r['producer']['inputs'].items():
            self.assertEqual(run.digest(run.ROOT / name), h, name)
        sizes = [o['length_m'] for o in r['long_haul']['options']]
        self.assertEqual(sizes, sorted(sizes))
        berths = [o['crown']['berths'] for o in r['long_haul']['options']]
        self.assertEqual(berths, sorted(berths, reverse=True))


if __name__ == '__main__':
    unittest.main()
