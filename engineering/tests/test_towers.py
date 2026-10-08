"""The tower model reproduces beam theory and the constant-stress column; the wind-energy helpers keep their limits."""
from __future__ import annotations
import math
import unittest
import numpy as np

from engineering.towers import lattice as lt
from engineering.towers import wind_energy as we

STILL = lambda z: np.zeros_like(np.asarray(z, dtype=float))


class TowerSizing(unittest.TestCase):
    def test_constant_stress_column_grows_exponentially(self):
        g, top, H = 1.6242, 1e7, 20e3
        frame = lt.Frame(bracing_fraction=0.0, top_mass_kg=top, top_width_m=300.0, min_leg_area_m2=1e-6)
        r = lt.size(H, 300.0, 1.0, 0.0, lt.World('test', g, STILL), frame, floors=None, steps=2000)
        length = frame.allowable_pa / (frame.material.density_kg_m3 * g)
        self.assertAlmostEqual(r['frame_mass_kg'][0] / (top * math.expm1(H / length)), 1.0, places=3)

    def test_cantilever_response_matches_beam_theory(self):
        n, H, EI, m, w, q = 4000, 100.0, 2e9, 50.0, 30.0, 400.0
        dz = np.array([H / n])
        x = (np.arange(n) + 0.5)[:, None] * dz
        sway, omega2, buckling = lt.response(dz, np.full((n, 1), EI), np.full((n, 1), m), np.full((n, 1), w),
                                             w * (H - x) ** 2 / 2, q * (H - x))
        self.assertAlmostEqual(sway[0] / (w * H ** 4 / (8 * EI)), 1.0, places=5)
        # Rayleigh with the uniform-load deflection: 3.530 against the exact 3.516; Rayleigh-Ritz for Greenhill's
        # self-weight column with the same shape: 8.0 against the exact 7.837.
        self.assertAlmostEqual(math.sqrt(omega2[0] * m * H ** 4 / EI), math.sqrt(12.4615), places=3)
        self.assertAlmostEqual(buckling[0] * q * H ** 3 / EI, 8.0, places=3)

    def test_lattice_force_coefficient(self):
        f = lt.Frame(solidity=0.2, round_members=False)
        self.assertAlmostEqual(f.force_coefficient(), (4.0 * 0.04 - 5.9 * 0.2 + 4.0) * 1.15)
        rr = 0.57 - 0.14 * 0.2 + 0.86 * 0.04 - 0.24 * 0.008
        self.assertAlmostEqual(lt.Frame(solidity=0.2).force_coefficient(), f.force_coefficient() * rr)
        self.assertAlmostEqual(lt.Frame(solidity=0.4, round_members=False).force_coefficient(),
                               (4.0 * 0.16 - 5.9 * 0.4 + 4.0) * 1.2)

    def test_wind_adds_mass_and_the_moon_needs_less_than_earth(self):
        air = lambda z: np.full_like(np.asarray(z, dtype=float), 1.2)
        moon, earth = lt.World('moon', 1.6242, air), lt.World('earth', 9.81, air)
        wide = lt.Frame(top_width_m=lt.Floors().full_width_m)
        calm = lt.size(5e3, 1000.0, 2.0, 0.0, moon, wide)
        windy = lt.size(5e3, 1000.0, 2.0, 40.0, moon, wide)
        self.assertGreater(windy['frame_mass_kg'][0], calm['frame_mass_kg'][0])
        self.assertGreater(lt.size(5e3, 1000.0, 2.0, 40.0, earth, wide)['frame_mass_kg'][0], windy['frame_mass_kg'][0])
        self.assertAlmostEqual(windy['floor_area_m2'][0] / (5e3 * lt.Floors().area_per_m(1e3)), 1.0, places=6)
        narrow = lt.size(5e3, 200.0, 1.0, 40.0, moon, lt.Frame(top_width_m=200.0))
        self.assertAlmostEqual(narrow['floor_area_m2'][0], 5e3 * lt.Floors().storeys * 0.5 * 200.0 ** 2 / 200.0, places=0)

    def test_added_drag_raises_the_base_moment(self):
        air = lambda z: np.full_like(np.asarray(z, dtype=float), 1.0)
        world = lt.World('moon', 1.6242, air)
        a = lt.size(3e3, 600.0, 1.0, 30.0, world)
        b = lt.size(3e3, 600.0, 1.0, 30.0, world, added_drag=0.3)
        self.assertGreater(b['base_moment_nm'][0], a['base_moment_nm'][0])
        self.assertGreater(b['frame_mass_kg'][0], a['frame_mass_kg'][0])


class WindEnergy(unittest.TestCase):
    def test_weibull_fit_recovers_its_percentiles(self):
        k, c = we.weibull_from_percentiles(2.5, 4.5)
        self.assertAlmostEqual(we.weibull_quantile(k, c, 0.5), 2.5)
        self.assertAlmostEqual(we.weibull_quantile(k, c, 0.9), 4.5)
        self.assertAlmostEqual(we.weibull_moment(2.0, 1.0, 2.0), 1.0)

    def test_device_output_stays_below_its_share_of_the_resource(self):
        k, c, rho = 2.0, 5.0, 1.0
        dev = we.Device('test', 0.4, 0.0, 1e3, 1e3, 0.8, 0.1)
        self.assertAlmostEqual(we.mean_output_w_m2(dev, k, c, rho) / (0.4 * we.resource_w_m2(k, c, rho)), 1.0, places=3)
        slow = we.Device('slow', 0.3, 3.0, 8.0, 20.0, 0.9, 0.6)
        self.assertLess(we.mean_output_w_m2(slow, k, c, rho), 0.3 * we.resource_w_m2(k, c, rho))
        self.assertLess(we.working_share(slow, k, c), 1.0)

    def test_collision_probability_limits(self):
        rotor = we.Device('rotor', 0.4, 3.0, 11.0, 25.0, 0.8, 0.1, tip_speed_ratio=6.0, blades=3, chord_share=0.06,
                          pitch_deg=5.0)
        self.assertEqual(we.band_collision_probability(we.Device('mast', 0.05, 3, 8, 15, 2, 1.2), 5.0, 5.0, 0.1, 0.2, 5.0), 0.0)
        slow_flier = we.band_collision_probability(rotor, 40.0, 6.0, 0.2, 0.35, 2.0)
        fast_flier = we.band_collision_probability(rotor, 40.0, 6.0, 0.2, 0.35, 10.0)
        self.assertTrue(0.0 < fast_flier < slow_flier <= 1.0)


if __name__ == '__main__':
    unittest.main()
