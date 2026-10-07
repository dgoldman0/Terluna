"""Limb heat: the aperture's zones, the ray tally's energy budget, the self-consistent state, and the product.

Tests that need the CXRO and Leiden absorption tables skip when fetch_limb_inputs has not restored them; the
product's tests skip when results/limb_heat.json has not been written.
"""
from __future__ import annotations
import json
import math
import unittest

import numpy as np

from shared.provenance import constants_changed
from atmosphere.loss_response import traced
from atmosphere.middle_atmosphere import fetch_limb_inputs
from atmosphere.middle_atmosphere import limb_heat as lh

R = lh.MOON_RADIUS


def _present(module):
    try:
        return module.restore()['status'] == 'PASS'
    except (OSError, ValueError):
        return False


NEEDS_TABLES = unittest.skipUnless(_present(fetch_limb_inputs),
                                   'absorption tables not restored (fetch_limb_inputs --download)')


def exponential_air(n0=1e20, scale=50e3, top=2 * R, points=4000):
    """A grey test atmosphere: one species falling off exponentially from the surface."""
    r = R + np.linspace(0., top - R, points)
    zero = np.full(points, 1e-300)
    return dict(radius_m=r, density=dict(N2=n0 * np.exp(-(r - R) / scale), O2=zero, O=zero))


class ZoneTests(unittest.TestCase):
    def test_the_zones_share_out_each_ray(self):
        b = np.linspace(0., 6 * R, 2000)
        zones = lh.zone_weights(b)
        np.testing.assert_allclose(zones['window'] + zones['annulus'] + zones['outside'], 1., atol=1e-12)
        self.assertEqual(zones['window'][0], 1.)
        self.assertAlmostEqual(float(zones['annulus'][np.argmin(abs(b - 2.5 * R))]), 1.)
        self.assertEqual(zones['outside'][-1], 1.)

    def test_the_rays_cover_the_disk_and_every_aperture_they_reach(self):
        b = lh.impact_parameters()
        areas = lh.ring_areas(b)
        self.assertAlmostEqual(areas[:lh.DESIGN['disk_rays']].sum() / (math.pi * R**2), 1., places=12)
        for x in (2, 2.5, 3, 4, 5, 6):
            zones = lh.zone_weights(b, x)
            self.assertAlmostEqual((areas * zones['aperture']).sum() / (math.pi * lh.aperture_edge(x)**2), 1., places=2)


class DepositTests(unittest.TestCase):
    def setUp(self):
        self.air = exponential_air()
        self.sigma = dict(N2=np.array([1e-26, 1e-24, 1e-22]), O2=np.zeros(3), O=np.zeros(3))

    def deposit(self, b):
        """Each wavelength's share left in each profile shell, and on the ground."""
        shares, shell, ground = lh.deposit(b, self.air, self.sigma)
        shells = np.zeros((shares.shape[0], len(self.air['radius_m']) - 1))
        np.add.at(shells.T, shell, shares.T)
        return shells, ground

    def test_a_ray_onto_the_disk_leaves_its_energy_in_the_air_or_on_the_ground(self):
        shells, ground = self.deposit(.3 * R)
        np.testing.assert_allclose(shells.sum(axis=1) + ground, 1., rtol=1e-10)

    def test_a_vertical_ray_sees_the_vertical_column(self):
        _, ground = self.deposit(0.)
        column = 1e20 * 50e3 * -math.expm1(-R / 50e3)
        np.testing.assert_allclose(ground, np.exp(-self.sigma['N2'] * column), rtol=1e-6)

    def test_a_ray_past_the_limb_crosses_the_air_twice(self):
        b = R + 200e3
        shells, ground = self.deposit(b)
        self.assertTrue(np.all(ground == 0.))
        # Half the grazing column, by direct quadrature along the ray from its closest approach outward. Near the
        # tangent point the density falls as a Gaussian in path length, which exponential segments follow to 2e-4.
        s = np.linspace(0., math.sqrt((2 * R)**2 - b**2), 400001)
        half = np.trapezoid(1e20 * np.exp(-(np.hypot(b, s) - R) / 50e3), s)
        np.testing.assert_allclose(shells.sum(axis=1), -np.expm1(-2 * self.sigma['N2'] * half), rtol=5e-4)

    def test_each_shell_is_credited_where_the_ray_crosses_it(self):
        b = R + 200e3
        shells, _ = self.deposit(b)
        r = self.air['radius_m']
        self.assertTrue(np.all(shells[:, r[:-1] < b - 1e3] == 0.))
        self.assertGreater(shells[0, np.searchsorted(r, b) - 1], 0.)


class StateTests(unittest.TestCase):
    def setUp(self):
        self.q = np.linspace(0., 1., 11)
        zero = np.zeros(11)
        self.parts = dict(window=.2 + .5 * self.q, annulus=zero, beyond_aperture=zero, glow=zero,
                          gap_unit=np.ones(11))

    def test_the_state_is_the_first_heat_the_air_returns(self):
        self.assertAlmostEqual(traced.first_state(self.q, self.parts, 0.), .4, places=9)
        self.assertAlmostEqual(traced.first_state(self.q, self.parts, .1), .6, places=9)

    def test_a_heat_that_outgrows_itself_runs_away(self):
        self.parts['window'] = .2 + 1.1 * self.q
        self.assertIsNone(traced.first_state(self.q, self.parts, 0.))

    def test_the_allowed_transmission_inverts_the_state(self):
        self.assertAlmostEqual(traced.allowed(self.q, self.parts, .6), .1, places=6)
        self.assertIsNone(traced.allowed(self.q, self.parts, .3))

    def test_the_share_from_beyond_the_aperture_is_read_at_the_state(self):
        self.parts['beyond_aperture'] = .1 * self.q
        state = traced.first_state(self.q, self.parts, 0.)
        self.assertAlmostEqual(state, .2 / .4, places=9)
        self.assertAlmostEqual(traced.shares(self.q, self.parts, 0., state)['beyond_share'], .1, places=9)


@NEEDS_TABLES
class CrossSectionTests(unittest.TestCase):
    def test_the_atomic_and_molecular_tables_meet_at_the_crossover(self):
        w = np.array([24.9, 25.0])
        sigma = lh.cross_sections(w)
        for species, values in sigma.items():
            self.assertLess(abs(math.log(values[0] / values[1])), math.log(1.6), species)

    def test_soft_x_rays_are_absorbed_more_weakly_than_the_extreme_ultraviolet(self):
        sigma = lh.cross_sections(np.array([1.0, 30.0, 80.0]))
        for species, values in sigma.items():
            self.assertGreater(values[1], 30 * values[0], species)
            self.assertGreater(values[2], 1e-22, species)


@unittest.skipUnless(lh.OUT.is_file(), 'limb_heat.json not written')
class ProductTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.product = json.loads(lh.OUT.read_text())

    def test_product_is_current_with_its_code(self):
        self.assertEqual(self.product['schema'], lh.SCHEMA)
        for name, path in lh.FILES.items():
            self.assertEqual(self.product['producer']['files'][name], lh.digest(path), name)
        for name, value in self.product['producer']['products'].items():
            self.assertEqual(lh.digest(lh.ROOT / name), value, name)
        self.assertFalse(constants_changed(self.product['producer']['constants']))
        self.assertEqual(sorted(self.product['cases']), sorted(lh.CASES))

    def test_the_disk_rays_reproduce_the_escape_models_count(self):
        counts = self.product['escape_count']
        for case in self.product['cases'].values():
            table = case['by_profile_heat'][0]['table']
            for activity in ('quiet', 'solar_maximum'):
                window = table['window']['disk'][activity]['above_base_W_m2']
                self.assertAlmostEqual(window / counts['window_film_heat'][activity], 1., delta=.01)
                # The escape model absorbs Lyman-alpha along a mean slant path; the rays find it within 10 per cent.
                open_band = table['open']['disk'][activity]['above_base_W_m2']
                self.assertAlmostEqual(open_band / counts['band_over_disk'][activity], 1., delta=.1)

    def test_every_zone_accounts_for_the_light_arriving(self):
        for case in self.product['cases'].values():
            for entry in case['by_profile_heat']:
                for source in entry.get('table', {}).values():
                    for zone in source.values():
                        for row in zone.values():
                            self.assertAlmostEqual(sum(row['shares'].values()), 1., places=9)
                            self.assertGreater(row['shares']['passing'], -1e-6)

    def test_light_through_gaps_heats_more_than_the_escape_model_counts(self):
        counts = self.product['escape_count']
        for case in self.product['cases'].values():
            table = case['by_profile_heat'][0]['table']
            ratio = table['open']['aperture']['quiet']['thermosphere_W_m2'] / counts['band_over_disk']['quiet']
            self.assertGreater(ratio, 2.)

    def test_the_bands_weighted_by_activity_give_the_detailed_table(self):
        bands = self.product['design']['bands_nm']
        for case in self.product['cases'].values():
            for entry in case['by_profile_heat']:
                if 'table' not in entry:
                    continue
                for source, zone in lh.COMPACT:
                    for activity, row in entry['table'][source][zone].items():
                        heat = float(np.dot(entry['radii']['4'][source][zone], traced.scales(activity, bands)))
                        self.assertTrue(math.isclose(heat, row['thermosphere_W_m2'], rel_tol=1e-9, abs_tol=1e-300))

    def test_each_parts_heat_lies_between_the_base_and_the_exobase(self):
        step = self.product['design']['log_pressure_step']
        for case in self.product['cases'].values():
            for entry in case['by_profile_heat']:
                for zone in entry.get('shapes', {}).values():
                    for groups in zone.values():
                        for quantiles in groups:
                            if quantiles is None:
                                continue
                            self.assertTrue(all(b >= a for a, b in zip(quantiles, quantiles[1:])))
                            self.assertGreaterEqual(quantiles[0], 0.)
                            self.assertLessEqual(quantiles[-1], entry['exobase_log_pressure'] + 2 * step)

    def test_lyman_alpha_through_gaps_heats_low_and_the_extreme_ultraviolet_high(self):
        groups = self.product['design']['shape_groups_nm']
        euv, lyman = groups.index(10.0), groups.index(121.0)
        middle = self.product['design']['shape_quantiles'].index(0.5)
        for case in self.product['cases'].values():
            entry = case['by_profile_heat'][0]
            gaps = entry['shapes']['aperture_4']['open']
            self.assertLess(gaps[lyman][middle], .5 * gaps[euv][middle])

    def test_a_wider_aperture_never_lets_more_unfiltered_light_reach_the_thermosphere(self):
        for case in self.product['cases'].values():
            for entry in case['by_profile_heat']:
                if 'radii' not in entry:
                    continue
                beyond = [sum(entry['radii'][f'{x:g}']['open']['outside']) for x in lh.DESIGN['summary_radii_R']]
                self.assertTrue(all(b <= a * (1 + 1e-9) + 1e-30 for a, b in zip(beyond, beyond[1:])))

    def test_the_thicker_film_settles_cooler(self):
        for case in self.product['cases'].values():
            for activity in ('quiet', 'solar_maximum'):
                thin = case['scenarios'][f'annulus_2_um_{activity}_gaps_0']['traced']
                thick = case['scenarios'][f'annulus_4_um_{activity}_gaps_0']['traced']
                if 'heat_W_m2' in thin and 'heat_W_m2' in thick:
                    self.assertLess(thick['heat_W_m2'], thin['heat_W_m2'])

    def test_fism2_puts_the_films_x_rays_near_twenty_times_quiet_at_maximum(self):
        sources = self.product['xray_cycle']['sources']
        for period in ('fism2', *lh.MAXIMA):
            for film in ('window', 'annulus_2_um', 'annulus_4_um'):
                self.assertTrue(10. < sources[period][film] < 30., (period, film))
            self.assertTrue(3. < sources[period]['open'] < 8., period)

    def test_measured_x_rays_leave_the_4_um_film_well_under_the_glow(self):
        for name in ('titania_stack_collisional', 'titania_stack_lte', 'titania_stack_all_heats'):
            for activity in ('solar_maximum', *lh.MAXIMA):
                row = self.product['cases'][name]['scenarios'][f'annulus_4_um_{activity}_gaps_0']
                parts = row['traced']['parts_W_m2']
                self.assertLess(parts['annulus'], .3 * parts['glow'])

    def test_tracing_cuts_the_gaps_allowed_to_a_quarter_or_a_third_in_the_cooler_cases(self):
        for name in ('titania_stack_collisional', 'titania_stack_lte'):
            for row in self.product['cases'][name]['budgets']:
                if row['film'] != 'annulus_4_um' or row['activity'] not in ('quiet', 'solar_maximum'):
                    continue
                self.assertTrue(.2 < row['traced'] / row['disk_count'] < .4, (name, row['activity'], row['budget_kg_s']))

    def test_tracing_never_allows_more_transmission_than_the_count(self):
        for case in self.product['cases'].values():
            for row in case['budgets']:
                if row.get('traced') and row.get('disk_count'):
                    self.assertLess(row['traced'], row['disk_count'])


if __name__ == '__main__':
    unittest.main()
