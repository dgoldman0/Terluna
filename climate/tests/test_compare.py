"""The run comparison's arithmetic on synthetic output (no model files)."""
import unittest
import numpy as np

from climate.gcm import compare


def year(lat, nlon=8, times=12, clear=False):
    shape = (times, len(lat), nlon)
    t = np.arange(times)[:, None, None] * np.ones(shape)
    land = np.zeros(shape); land[:, :, :nlon // 2] = 1.0
    fields = dict(ts=np.full(shape, 290.0) + 10.0 * (1 - land), tas=288.0 + (t % 2) * 6.0 * land,
                  ntr=np.full(shape, 0.5), rst=np.full(shape, 240.0), rsut=np.full(shape, -100.0),
                  rlut=np.full(shape, -239.5), clt=np.full(shape, 0.4), pr=np.full(shape, 3.0 / 86400e3),
                  prw=np.full(shape, 200.0), sic=np.zeros(shape), lsm=land,
                  czen=0.6 * (t % 2))                                    # the warm times are sunlit
    if clear:
        fields.update(rstcs=np.full(shape, 260.0), rltcs=np.full(shape, -249.5))
    return fields


class CompareTests(unittest.TestCase):
    lat = np.degrees(np.arcsin(np.polynomial.legendre.leggauss(16)[0]))[::-1]      # includes a row near 60°

    def test_weights_cover_the_sphere(self):
        w = compare.area_weights(self.lat, 8)
        self.assertAlmostEqual(float(w.sum()), 1.0)
        np.testing.assert_allclose(w, w[::-1])                     # symmetric about the equator

    def test_summary(self):
        out = compare.summarise([year(self.lat), year(self.lat)], self.lat)
        self.assertAlmostEqual(out['surface_k'], 295.0)            # half land at 290 K, half sea at 300 K
        self.assertAlmostEqual(out['sea_surface_k'], 300.0)
        self.assertAlmostEqual(out['planetary_albedo'], 100.0 / 340.0)
        self.assertAlmostEqual(out['precipitation_mm_day'], 3.0)
        self.assertAlmostEqual(out['equator_land_range_k'], 6.0)  # alternating 288 and 294 K over land
        self.assertAlmostEqual(out['equator_land_day_k'], 294.0)
        self.assertAlmostEqual(out['equator_land_night_k'], 288.0)
        self.assertNotIn('cloud_effect_net_w_m2', out)
        self.assertAlmostEqual(out['near_land_air_k'], 291.0)      # land spans the near side and beyond it
        self.assertAlmostEqual(out['far_land_air_k'], 291.0)
        self.assertAlmostEqual(out['land_air_p95_k'], 291.0)
        self.assertAlmostEqual(out['dry_land_share'], 0.0)         # 3 mm/day everywhere

    def test_weighted_percentile(self):
        values, weights = np.array([1.0, 2.0, 3.0, 4.0]), np.array([1.0, 1.0, 1.0, 7.0])
        self.assertEqual(compare.weighted_percentile(values, weights, 50), 4.0)   # the heavy value holds the middle
        self.assertEqual(compare.weighted_percentile(values, np.ones(4), 25), 1.0)

    def test_cloud_effect(self):
        out = compare.summarise([year(self.lat, clear=True)], self.lat)
        self.assertAlmostEqual(out['cloud_effect_sw_w_m2'], -20.0)
        self.assertAlmostEqual(out['cloud_effect_lw_w_m2'], 10.0)
        self.assertAlmostEqual(out['cloud_effect_net_w_m2'], -10.0)


    def test_level_means_weigh_the_globe_the_tropical_land_and_its_day_and_night(self):
        lat = np.array([60.0, 10.0, -10.0, -60.0])
        land = np.zeros((4, 8), bool); land[:, :4] = True
        ta = np.ones((2, 3, 4, 8)) * np.array([250.0, 270.0, 290.0])[None, :, None, None]
        ta[:, :, 1:3, :4] += 2.0                                       # the tropical land runs 2 K warmer at every level
        ta[1, 2, 1:3, :4] += 4.0                                       # and its lowest level 4 K more at the dark time
        sun = np.stack([np.full((4, 8), 0.6), np.zeros((4, 8))])        # the first time sunlit, the second dark
        out = compare.level_means(ta, land, lat, sun)
        np.testing.assert_allclose(out['tropical_land'], [252.0, 272.0, 294.0])
        np.testing.assert_allclose(out['tropical_land_sunlit'], [252.0, 272.0, 292.0])
        np.testing.assert_allclose(out['tropical_land_dark'], [252.0, 272.0, 296.0])
        w = compare.area_weights(lat, 8)
        tropic = 2.0 * w[1:3, :4].sum()
        np.testing.assert_allclose(out['globe'][:2], [250.0 + tropic, 270.0 + tropic])

    def test_settle(self):
        # Imbalance falling 1.2 W/m2 per K from +3 W/m2 at 295 K reaches a +0.3 W/m2 baseline at 297.25 K.
        temps = np.array([295.0, 295.5, 296.0, 296.4])
        t, feedback = compare.settle(temps, 3.0 - 1.2 * (temps - 295.0), 0.3)
        self.assertAlmostEqual(t, 297.25)
        self.assertAlmostEqual(feedback, 1.2)


if __name__ == '__main__':
    unittest.main()
