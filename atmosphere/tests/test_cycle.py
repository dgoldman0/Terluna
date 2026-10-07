"""The solar cycle of each band of the light below 175 nm, from FISM2's daily record (atmosphere/loss_response/cycle.py)."""
from __future__ import annotations
import datetime as dt
import json
import unittest

import numpy as np

from shared.provenance import constants_changed
from atmosphere.middle_atmosphere import escape, fetch_limb_inputs
from atmosphere.loss_response import cycle, traced


def _inputs_present():
    try:
        return traced.LIMB.is_file() and all(
            fetch_limb_inputs.path(cycle.file_name(lo, hi)).is_file()
            for lo, hi in zip(cycle.band_edges()[:-1], cycle.band_edges()[1:]))
    except (OSError, ValueError, KeyError):
        return False


@unittest.skipUnless(_inputs_present(), 'FISM2 band records not restored (fetch_limb_inputs --download)')
class CycleTests(unittest.TestCase):
    def test_the_quiet_week_is_the_unit(self):
        values, days = cycle.factors(*cycle.QUIET_WEEK)
        self.assertEqual(days, 7)
        self.assertTrue(np.allclose(values, 1.0, rtol=1e-12))

    def test_the_x_ray_bands_from_a_nanometre_up_match_the_escape_models_maxima(self):
        edges = cycle.band_edges()
        table = escape.xray_cycle()
        for name, (start, end) in cycle.MAXIMA.items():
            values, _ = cycle.factors(start, end)
            for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                if lo >= 1.0 and hi <= 10.0:
                    self.assertAlmostEqual(values[i] / table['maxima'][name][i], 1.0, delta=.05, msg=(name, lo))

    def test_the_far_ultraviolet_rises_far_less_than_the_conventions_factor(self):
        for name, (start, end) in cycle.MAXIMA.items():
            a = cycle.activity(name, start, end)
            self.assertLess(a['uv'], escape.SOLAR_MAXIMUM, name)
            far = [f for f, lo in zip(a['bands'], cycle.band_edges()[:-1]) if lo >= 122.0]
            self.assertLess(max(far), 1.5, name)
            self.assertGreater(min(far), 1.0, name)

    def test_a_measured_activity_weights_the_bands_it_names(self):
        a = cycle.activity('2001', dt.date(2001, 1, 1), dt.date(2001, 12, 31))
        self.assertTrue(np.array_equal(traced.scales(a), np.asarray(a['bands'])))
        w = np.array([0.3, 15.0, 121.55, 170.0])
        factors = traced.band_factors(w, a)
        self.assertEqual(factors[2], a['glow'])
        self.assertEqual(factors[3], a['bands'][-1])

    def test_the_years_cover_the_two_complete_cycles(self):
        spans = cycle.periods()
        years = [a['label'] for a in spans['cycles_23_24']['years']]
        self.assertEqual(years[0], '1997')
        self.assertEqual(years[-1], '2019')
        self.assertEqual(sum(a['days'] for a in spans['cycles_23_24']['years']),
                         (dt.date(2019, 12, 31) - dt.date(1997, 1, 1)).days + 1)


@unittest.skipUnless(cycle.OUT.is_file(), 'solar_cycle.json not written')
class ProductTests(unittest.TestCase):
    def test_product_is_current_with_its_code_and_inputs(self):
        product = json.loads(cycle.OUT.read_text())
        self.assertEqual(product['schema'], cycle.SCHEMA)
        for name, value in product['producer']['files'].items():
            self.assertEqual(cycle.digest(cycle.HERE.parents[1] / name), value, name)
        for name, value in product['producer']['inputs'].items():
            path = (fetch_limb_inputs.path(name) if name.startswith('fism2') else
                    __import__('atmosphere.middle_atmosphere.fetch_inputs', fromlist=['path']).path(name))
            self.assertEqual(cycle.digest(path), value, name)
        self.assertEqual(product['producer']['products']['atmosphere/middle_atmosphere/results/limb_heat.json'],
                         cycle.digest(traced.LIMB))
        self.assertFalse(constants_changed(product['producer']['constants']))


if __name__ == '__main__':
    unittest.main()
