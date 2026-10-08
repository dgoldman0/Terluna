"""The fire model reproduces its Earth correlations in Earth air and carries them to other gravity and air by Froude
scaling; the egress and design-fire helpers keep their published values."""
from __future__ import annotations
import math
import unittest
import numpy as np

from atmosphere.radiative_convective import thermodynamics as th
from engineering.fire import design_fires as df
from engineering.fire import egress as eg
from engineering.fire import smoke as sm
from shared.constants import MOON_SURFACE_GRAVITY

EARTH = sm.EARTH_REFERENCE
MOON_SAME_AIR = sm.Air(MOON_SURFACE_GRAVITY, EARTH.pressure_pa, EARTH.temperature_k)
MOON = sm.Air(MOON_SURFACE_GRAVITY, 0.83 * 101325.0, 293.15, th.earthlike_air(121590.0, 400.0))


class EarthCorrelations(unittest.TestCase):
    def test_flame_height_is_heskestads_in_earth_air(self):
        for q, d in ((300.0, 0.5), (1000.0, 1.0), (5000.0, 2.0)):
            heskestad = 0.235 * q ** 0.4 - 1.02 * d
            self.assertAlmostEqual(float(sm.flame_height_m(q, d, EARTH)) / heskestad, 1.0, delta=0.015)

    def test_ceiling_jet_constants_are_nfpa72s(self):
        arrival, offset, slope = sm.hd_constants()
        for got, published in ((arrival, 0.861), (offset, 0.146), (slope, 0.242)):
            self.assertAlmostEqual(got, published, delta=0.001)

    def test_flashover_is_mqhs_in_earth_air(self):
        ao, ho, hk, at = 1.89, 2.1, 0.03, 124.0
        q = float(sm.flashover_kw(ao, ho, hk, at, EARTH))
        self.assertAlmostEqual(q / (610.0 * math.sqrt(hk * at * ao * math.sqrt(ho))), 1.0, delta=0.03)
        self.assertAlmostEqual(float(sm.layer_rise_k(q, ao, ho, hk, at, EARTH)), sm.FLASHOVER_RISE_K, places=6)

    def test_ventilation_limit_is_the_eurocodes_in_earth_air(self):
        self.assertAlmostEqual(float(sm.ventilation_limit_kw(4.0, 2.25, EARTH)), 1400.0 * 4.0 * 1.5, places=6)

    def test_smoke_filling_is_slower_than_nfpa92s_first_smoke(self):
        # NFPA 92's steady-fire equation marks the first indication of smoke; the zone model's layer base comes later.
        q, h, area = 1000.0, 10.0, 400.0
        nfpa = math.exp((1.11 - 0.5) / 0.28) * (area / h ** 2) * h ** (4 / 3) / q ** (1 / 3)
        zone = sm.smoke_filling(lambda t: q, area, h, EARTH, clear_m=0.5 * h)['clear_s']
        self.assertTrue(1.0 < zone / nfpa < 2.5, zone / nfpa)


class OtherGravity(unittest.TestCase):
    def test_ceiling_jet_follows_froude_scaling(self):
        f = sm.froude(MOON)
        alpha, h, r, t = 0.0117, 3.5, 4.0, 400.0
        rise_m, u_m = sm.ceiling_jet(t, alpha, h, r, MOON)
        rise_e, u_e = sm.ceiling_jet(t / f['time'], alpha * f['heat'] * f['time'] ** 2, h, r, EARTH)
        self.assertAlmostEqual(float(rise_m) / float(rise_e * f['temperature']), 1.0, places=9)
        self.assertAlmostEqual(float(u_m) / float(u_e * f['velocity']), 1.0, places=9)

    def test_smoke_filling_follows_froude_scaling(self):
        f = sm.froude(MOON)
        moon = sm.smoke_filling(lambda t: 800.0, 500.0, 4.0, MOON, step_s=0.1)
        earth = sm.smoke_filling(lambda t: 800.0 * f['heat'], 500.0, 4.0, EARTH, step_s=0.1 / f['time'])
        self.assertAlmostEqual(moon['clear_s'] / (earth['clear_s'] * f['time']), 1.0, delta=0.01)

    def test_lunar_plume_entrains_as_g_rho_squared_over_t_to_a_third(self):
        q, z = 50.0, 30.0   # far above a small fire, where the first term dominates
        ratio = float(sm.plume_mass_flow_kg_s(q, z, MOON) / sm.plume_mass_flow_kg_s(q, z, EARTH))
        expected = (MOON.gravity_m_s2 * MOON.density ** 2 * EARTH.temperature_k
                    / (EARTH.gravity_m_s2 * EARTH.density ** 2 * MOON.temperature_k)) ** (1 / 3)
        self.assertAlmostEqual(ratio / expected, 1.0, delta=0.01)

    def test_flames_are_taller_in_weak_gravity_and_poorer_air(self):
        q, d = 1000.0, 1.0
        base = sm.flame_height_m(q, d, EARTH) + 1.02 * d
        weak = sm.flame_height_m(q, d, MOON_SAME_AIR) + 1.02 * d
        self.assertAlmostEqual(float(weak / base), (EARTH.gravity_m_s2 / MOON_SURFACE_GRAVITY) ** 0.2, places=9)
        self.assertGreater(float(sm.flame_height_m(q, d, MOON)), float(sm.flame_height_m(q, d, MOON_SAME_AIR)))

    def test_flashover_comes_at_a_smaller_fire_in_weak_gravity(self):
        args = (1.89, 2.1, 0.03, 124.0)
        ratio = float(sm.flashover_kw(*args, MOON_SAME_AIR) / sm.flashover_kw(*args, EARTH))
        self.assertAlmostEqual(ratio, (MOON_SURFACE_GRAVITY / EARTH.gravity_m_s2) ** 0.25, places=9)

    def test_wind_feeds_a_fire_regardless_of_gravity(self):
        earth = float(sm.wind_ventilation_limit_kw(100.0, 2.5, EARTH))
        self.assertAlmostEqual(earth, 0.6 * 100.0 / math.sqrt(2.0) * EARTH.density * 2.5 * EARTH.heat_per_air_kj_kg, places=6)
        self.assertAlmostEqual(float(sm.wind_ventilation_limit_kw(100.0, 2.5, MOON_SAME_AIR)), earth, places=6)
        self.assertGreater(earth, float(sm.ventilation_limit_kw(100.0, 2.0, EARTH)))

    def test_smoke_filling_names_what_governs(self):
        cool = sm.smoke_filling(lambda t: 200.0, 400.0, 10.0, EARTH)
        self.assertEqual(cool['governed_by'], 'smoke at head height')
        hot = sm.smoke_filling(lambda t: 3000.0, 30.0, 3.0, EARTH, loss=0.0)
        self.assertEqual(hot['governed_by'], 'hot layer overhead')

    def test_stack_pressure(self):
        outside = sm.Air(EARTH.gravity_m_s2, 101325.0, 273.15)
        dp = float(sm.stack_pressure_pa(100.0, outside, 293.15))
        self.assertAlmostEqual(dp, EARTH.gravity_m_s2 * 100.0 * outside.density * (1 - 273.15 / 293.15), places=6)
        self.assertLess(float(sm.stack_pressure_pa(100.0, outside, 263.15)), 0.0)


class Sensors(unittest.TestCase):
    def test_sensor_follows_a_step_exactly(self):
        t = np.linspace(0.0, 120.0, 241)
        out = sm.sensor_rise(t, np.full_like(t, 50.0), np.full_like(t, 2.0), 80.0)
        exact = 50.0 * (1 - np.exp(-math.sqrt(2.0) * t / 80.0))
        self.assertTrue(np.allclose(out, exact, atol=1e-9))

    def test_a_fast_sensor_answers_when_the_gas_does(self):
        args = (0.0469, 3.0, 3.0, EARTH)
        gas = float(sm.gas_rise_time_s(30.0, *args))
        self.assertAlmostEqual(sm.activation_time_s(30.0, 0.01, *args, step_s=0.05), gas, delta=0.2)
        self.assertGreater(sm.activation_time_s(30.0, 50.0, *args), gas)


class PeopleAndFires(unittest.TestCase):
    def test_hydraulic_model(self):
        self.assertAlmostEqual(float(eg.speed_m_s(0.3)), 1.19, places=2)
        self.assertAlmostEqual(float(eg.specific_flow(1.88)), 1.32, places=2)
        m = eg.movement_s(260.0, 30.0, [1.3, 1.3])
        self.assertAlmostEqual(m['queue_s'], 100.0, places=6)

    def test_design_fires(self):
        self.assertAlmostEqual(df.EN1991_ANNEX_E['office'].earth_fire().alpha_kw_s2, 1000 / 300 ** 2)
        self.assertAlmostEqual(df.burning_rate_factor(df.EARTH_OXYGEN), 1.0)
        self.assertEqual(df.burning_rate_factor(0.10), 0.0)
        self.assertAlmostEqual(df.burning_rate_factor(0.1746), 0.649, places=3)
        shifts = df.lunar_oxygen_shifts()
        self.assertEqual(max(max(s.values()) for s in shifts.values()), 5.9)
        self.assertEqual(df.materials_test_oxygen_percent(0.1746, 5.9), 23.5)


if __name__ == '__main__':
    unittest.main()
