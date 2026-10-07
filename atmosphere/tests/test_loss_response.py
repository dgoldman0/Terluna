"""Loss response: the energy-limited form, the Lyman-alpha glow, the solar-wind cap, and the sweep.

Tests that need the WHI solar spectrum and the stored middle-atmosphere results
skip with that reason when the spectrum has not been restored.
"""
from __future__ import annotations
import math
import unittest

from atmosphere.loss_response import model as lr
from atmosphere.middle_atmosphere import fetch_inputs


def _present(module):
    try:
        return module.restore()['status'] == 'PASS'
    except (OSError, ValueError):
        return False


NEEDS_UV = unittest.skipUnless(_present(fetch_inputs), 'ultraviolet inputs not restored (fetch_inputs --download)')


class ClosedForms(unittest.TestCase):
    def test_energy_limited_matches_feasibility_report(self):
        # The September feasibility report: 12,500 kg/s at 2 lunar radii and eta 0.1 for 4.64e-3 W/m^2.
        self.assertAlmostEqual(lr.energy_limited_loss(4.64e-3, 2.0, 0.1) / 12474.2, 1.0, places=3)

    def test_glow_flux_for_1000_rayleigh(self):
        # 1000 R of Lyman-alpha onto a surface facing open sky: about 4.08e-6 W/m^2 (protection/report.md).
        self.assertAlmostEqual(lr.lyman_glow_flux(1000.0) / 4.08e-6, 1.0, places=2)

    def test_glow_absorption_is_bounded(self):
        heat = lr.lyman_glow_heat(1.3)
        upper = lr.lyman_glow_flux() * 1.3 ** 2 * 0.4
        self.assertGreater(heat, 0.5 * upper)
        self.assertLess(heat, upper)

    def test_mass_loading_cap(self):
        # 5 protons per cm^3 at 400 km/s through a disc three lunar radii across.
        cap = lr.solar_wind_losses(3.0)['mass_loading_cap_kg_s']
        expected = 5e6 * 4e5 * lr.PROTON_KG * math.pi * (3 * lr.MOON_RADIUS) ** 2
        self.assertAlmostEqual(cap / expected, 1.0, places=12)
        self.assertAlmostEqual(cap, 0.286, places=3)

    def test_solar_wind_range_is_ordered(self):
        wind = lr.solar_wind_losses()
        self.assertLess(wind['low']['total_kg_s'], wind['central']['total_kg_s'])
        self.assertLess(wind['central']['total_kg_s'], wind['high']['total_kg_s'])
        self.assertLessEqual(wind['high']['pickup_kg_s'], wind['mass_loading_cap_kg_s'])

    def test_allowed_leak_interpolates_in_log_space(self):
        rows = [dict(leak_fraction=f, status='thermal_column', molecular_loss_kg_s=loss)
                for f, loss in ((1e-4, 0.01), (1e-3, 1.0), (1e-2, 100.0))]
        self.assertAlmostEqual(lr.allowed_leak(rows, 10.0), math.sqrt(1e-3 * 1e-2), places=12)
        self.assertIsNone(lr.allowed_leak(rows, 0.001))
        self.assertIsNone(lr.allowed_leak(rows, -1.0))
        self.assertEqual(lr.allowed_leak(rows, 1e3), 1e-2)


@NEEDS_UV
class Sweep(unittest.TestCase):
    def test_sweep_takes_the_traced_state(self):
        row = lr.euv_sweep('titania_stack', 'lte', 'quiet', glow=False, leaks=(1e-3,))[0]
        # Light through gaps reaches the air above the limb as well as the disk: several times a quarter of the
        # light that reaches the disk, the escape model's count.
        self.assertTrue(2.0 < row['leak_heat_w_m2'] / row['disk_count_heat_w_m2'] < 6.0)
        parts = row['leak_heat_w_m2'] + row['film_heat_w_m2'] + row['beyond_aperture_heat_w_m2'] + row['glow_heat_w_m2']
        self.assertAlmostEqual(parts / row['deposited_heat_w_m2'], 1.0, delta=1e-3)

    def test_the_sweep_marks_where_the_air_runs_away(self):
        rows = lr.euv_sweep('titania_stack', 'all_heats', 'solar_maximum')
        status = [r['status'] for r in rows]
        self.assertIn('runaway', status)
        first = status.index('runaway')
        self.assertTrue(all(s == 'runaway' for s in status[first:]))
        self.assertLess(rows[first - 1]['leak_fraction'], rows[first]['leak_fraction'])
        self.assertEqual(lr.loss_estimate(rows[first]), math.inf)

    def test_loss_rises_with_leak_inside_the_column_domain(self):
        rows = [r for r in lr.euv_sweep('titania_stack', 'lte', 'quiet') if r['status'] == 'thermal_column']
        losses = [r['molecular_loss_kg_s'] for r in rows]
        self.assertGreater(len(losses), 5)
        self.assertTrue(all(b > a for a, b in zip(losses, losses[1:])))

    def test_glow_raises_the_floor(self):
        with_glow = lr.euv_sweep('titania_stack', 'lte', 'quiet', glow=True, leaks=(1e-6,))[0]
        without = lr.euv_sweep('titania_stack', 'lte', 'quiet', glow=False, leaks=(1e-6,))[0]
        self.assertGreater(with_glow['exobase_temperature_k'], without['exobase_temperature_k'] + 20)


if __name__ == '__main__':
    unittest.main()
