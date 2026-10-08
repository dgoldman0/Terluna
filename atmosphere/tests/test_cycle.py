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
        return traced.LIMB.is_file() and all(fetch_limb_inputs.path(name).is_file() for name in cycle.file_names())
    except (OSError, ValueError, KeyError):
        return False


@unittest.skipUnless(_inputs_present(), 'FISM2 band records not restored (fetch_limb_inputs --download)')
class CycleTests(unittest.TestCase):
    def test_the_quiet_week_is_the_unit(self):
        values, days = cycle.factors(*cycle.QUIET_WEEK)
        self.assertEqual(days, 7)
        self.assertTrue(np.allclose(values, 1.0, rtol=1e-12))

    def test_the_record_runs_day_by_day_from_its_start(self):
        dates, values = cycle.record()
        self.assertEqual((dates[0], dates[-1]), (dt.date(1947, 2, 14), dt.date(2025, 12, 31)))
        self.assertEqual(len(dates), (dates[-1] - dates[0]).days + 1)
        # FISM2 gives zero on some days only in the X-rays below 5 nm, most of them in the hardest band; kept as given.
        self.assertTrue(np.all(values >= 0))
        edges = cycle.band_edges()
        self.assertTrue(all(np.all(values[i] > 0) for i, hi in enumerate(edges[1:]) if hi > 5.0))

    def test_each_maximum_is_the_year_around_its_sunspot_maximum(self):
        self.assertEqual(cycle.year_around(2001, 11), (dt.date(2001, 5, 1), dt.date(2002, 4, 30)))
        self.assertEqual(cycle.year_around(2014, 4), (dt.date(2013, 10, 1), dt.date(2014, 9, 30)))
        self.assertEqual(cycle.year_around(1979, 12), (dt.date(1979, 6, 1), dt.date(1980, 5, 31)))
        for start, end in cycle.MAXIMA.values():
            self.assertIn((end - start).days + 1, (365, 366))

    def test_the_x_ray_bands_from_a_nanometre_up_match_the_escape_models_maxima(self):
        edges = cycle.band_edges()
        table = escape.xray_cycle()
        for name in table['maxima']:
            values, _ = cycle.factors(*cycle.MAXIMA[name])
            for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                if lo >= 1.0 and hi <= 10.0:
                    self.assertAlmostEqual(values[i] / table['maxima'][name][i], 1.0, delta=.05, msg=(name, lo))

    def test_the_far_ultraviolet_rises_far_less_than_the_conventions_factor(self):
        for name, (start, end) in cycle.MAXIMA.items():
            a = cycle.activity(name, start, end)
            self.assertLess(a['uv'], escape.SOLAR_MAXIMUM, name)
            far = [f for f, lo in zip(a['bands'], cycle.band_edges()[:-1]) if lo >= 122.0]
            self.assertLess(max(far), 1.65, name)
            self.assertGreater(min(far), 1.0, name)

    def test_the_solar_maximum_is_the_strongest_year_since_1978_taken_whole(self):
        measured = cycle.maxima()
        top, record = cycle.strongest(measured), cycle.strongest(measured, None, 'solar_maximum_cycle_19')
        self.assertEqual((top['cycle'], record['cycle']), ('cycle_21', 'cycle_19'))
        self.assertEqual((top['label'], record['label']), ('solar_maximum', 'solar_maximum_cycle_19'))
        recent = [m for m in measured.values() if dt.date.fromisoformat(m['start']) >= cycle.SATELLITE_ERA]
        self.assertTrue(all(top['uv'] >= m['uv'] for m in recent))
        self.assertTrue(all(record['uv'] >= m['uv'] for m in measured.values()))
        for key in ('bands', 'glow', 'xray', 'uv', 'start', 'end'):
            self.assertEqual(top[key], measured['cycle_21'][key])

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
        from atmosphere.middle_atmosphere import limb_heat
        self.assertEqual(product['bands_nm'], limb_heat.DESIGN['bands_nm'])
        self.assertEqual(product['solar_maximum'], cycle.strongest(product['maxima']))
        self.assertEqual(product['record_maximum'], cycle.strongest(product['maxima'], None, 'solar_maximum_cycle_19'))
        self.assertFalse(constants_changed(product['producer']['constants']))


if __name__ == '__main__':
    unittest.main()
