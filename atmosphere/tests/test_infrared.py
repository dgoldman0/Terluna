"""What the thermal column radiates with: its coolants above the base and CO2's band escape (loss_response/infrared.py)."""
from __future__ import annotations
import unittest

import numpy as np

from atmosphere.radiative_convective import fetch_inputs as rc_inputs
from atmosphere.loss_response import infrared as ir, model as lr


def _lines_present():
    try:
        return all(rc_inputs.path(n).is_file() for n in ('hitran_co2_0-3500.par', 'hitran_co2_3500-10000.par'))
    except (OSError, ValueError, KeyError):
        return False


class Coolants(unittest.TestCase):
    def test_co2_stays_mixed_where_it_radiates_and_separates_above_its_homopause(self):
        base = lr.base_conditions('titania_stack', 'all_heats')
        x = ir.LOG_PRESSURE
        co2 = ir.co2_profile(base['case'], lr.BASE_PA, base['base_temperature_k'], 0.0289)
        self.assertAlmostEqual(co2[0], ir.CO2_PPM * 1e-6)
        self.assertTrue(np.all(np.diff(co2) <= 0))
        self.assertGreater(co2[x == 4.0][0], .99 * co2[0])
        # Far above the homopause each species keeps its own scale height: d ln x / d ln(1/p) = -(m_CO2/m - 1).
        slope = -np.diff(np.log(co2))[-20:] / np.diff(x)[-20:]
        np.testing.assert_allclose(slope, ir.CO2_MOLAR / 0.0289 - 1, rtol=.02)

    def test_the_bases_atoms_and_no_come_from_the_middle_atmosphere(self):
        titania = lr.base_conditions('titania_stack', 'collisional')['case']
        edge = lr.base_conditions('edge_200nm', 'lte')['case']
        self.assertLess(ir.base_mixing(titania, 'O', lr.BASE_PA), 1e-15)        # no light splits O2 behind titania
        self.assertGreater(ir.base_mixing(edge, 'O', lr.BASE_PA), 1e-5)
        self.assertGreater(ir.base_mixing(edge, 'NO', lr.BASE_PA), 1e-9)


@unittest.skipUnless(_lines_present(), 'HITRAN CO2 lines not restored (radiative_convective.fetch_inputs --download)')
class BandEscape(unittest.TestCase):
    def test_the_band_escape_falls_from_one_and_is_thick_over_the_base(self):
        columns, escape = ir.band_escape(190.0)
        self.assertGreater(escape[0], .999)
        self.assertTrue(np.all(np.diff(escape) <= 0))
        # 400 ppm of the air above the 0.3 Pa base, about 2e17 CO2 cm^-2: the band's strong lines are thick.
        above = 4e-4 * .3 / (.0289 / 6.02214076e23 * 1.4) * 1e-4
        self.assertTrue(.005 < np.interp(np.log(above), np.log(columns), escape) < .1)

    def test_warmer_air_lets_more_of_the_band_out(self):
        cold, warm = ir.band_escape(150.0)[1], ir.band_escape(250.0)[1]
        middle = slice(30, 70)
        self.assertTrue(np.all(warm[middle] >= cold[middle]))

    def test_the_loss_response_radiates_unless_told_not_to(self):
        cooled = lr.column_config('titania_stack', 'lte')
        plain = lr.column_config('titania_stack', 'lte', cooling=False)
        self.assertEqual(bool(cooled.coolant_log_pressure), lr.INFRARED_COOLING)
        self.assertFalse(plain.coolant_log_pressure)
        self.assertEqual(cooled.lower_temperature_k, plain.lower_temperature_k)


if __name__ == '__main__':
    unittest.main()
