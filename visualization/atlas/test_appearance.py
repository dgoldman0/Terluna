"""The globe's true-colour pieces: the soil and plant colours, the rivers' sediment-bearing flow, the plumes it makes and
the sky each place sees."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from visualization.atlas import appearance as A


class ColourTests(unittest.TestCase):
    def test_a_munsell_grey_keeps_its_reflectance_in_every_channel(self):
        # Value 5 is 19.77% of magnesium oxide in the renotation, 19.3% of the perfect diffuser.
        grey = A.munsell_srgb((0.31006, 0.31616, 19.77))
        self.assertTrue(np.allclose(grey, 0.193, atol=0.002), grey)

    def test_soils_are_the_anchors_at_the_anchors(self):
        for end, feo in (('highland', A.SOIL_FEO[0]), ('basalt', A.SOIL_FEO[1])):
            for aridity, _, _, xyY in A.SOIL_MUNSELL[end]:
                self.assertTrue(np.allclose(A.soil_colour(aridity, feo, 2.0), A.munsell_srgb(xyY), atol=1e-5))

    def test_basalt_soil_is_darker_and_titanium_darkens_it(self):
        lum = lambda c: float(A.LUMINANCE @ c)
        for aridity in (0.02, 0.3, 1.0):
            self.assertLess(lum(A.soil_colour(aridity, 20.0, 1.0)), lum(A.soil_colour(aridity, 4.5, 0.3)))
        self.assertLess(lum(A.soil_colour(0.1, 20.0, 9.0)), lum(A.soil_colour(0.1, 20.0, 1.0)))

    def test_the_cover_follows_the_stated_curve(self):
        for aridity, share in ((0.1, 0.17), (0.2, 0.46), (1.0, 0.97)):
            self.assertAlmostEqual(float(A.plant_cover(aridity)), share, delta=0.01)

    def test_the_ground_greens_and_darkens_with_the_rain(self):
        soil = lambda a: A.soil_colour(a, 4.5, 0.3)
        dry, wet = A.ground_colour(0.05, soil(0.05)), A.ground_colour(2.0, soil(2.0))
        self.assertGreater(wet[1] / wet[0], dry[1] / dry[0])
        self.assertLess(float(A.LUMINANCE @ wet), float(A.LUMINANCE @ dry))


def routing(path, receiver, discharge, sea, lake_cells, neck=None):
    n = len(receiver)
    np.savez(path, metadata=json.dumps(dict(neck=neck)), receiver=np.asarray(receiver, np.int32),
             discharge_m3s=np.asarray(discharge, float), sea=np.asarray(sea, np.int8),
             routing_lake=np.zeros(n, bool), outlet=np.asarray(sea) >= 0, lake_cell=np.asarray(lake_cells, np.int32))


class SedimentTests(unittest.TestCase):
    def test_lakes_and_seas_pass_no_sediment_on_and_the_neck_spill_is_clear(self):
        # 0 -> 1 -> 2 (a lake) -> 3 -> 4 -> 5 (the sea); 1 m3/s of runoff on each dry cell; 10 m3/s of spill put in at 3.
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'routing.npz'
            routing(path, [1, 2, 3, 4, 5, 5], [1, 2, 2, 13, 14, 14], [-1, -1, -1, -1, -1, 0], [2],
                    neck=dict(spill_cell=3, flow_m3s=10.0))
            carried = A.land_discharge(path)
        self.assertTrue(np.allclose(carried, [1, 2, 0, 1, 2, 0]), carried)

    def test_no_mouth_no_plume_and_one_mouth_one_round_plume(self):
        water = np.zeros((40, 80), bool)
        water[:, 40:] = True
        lat = 1.0 - (np.arange(40) + 0.5) * 0.05             # near the equator, 2 km pixels
        self.assertEqual(float(A.plume_share(np.zeros(water.shape, np.float32), water, lat, 2.0).max()), 0.0)
        flow = np.zeros(water.shape, np.float32)
        flow[20, 39] = 5000.0                               # a plume about 17 km, 8 pixels, across
        plume = A.plume_share(flow, water, lat, 2.0)
        self.assertGreater(plume[20, 40], 0.5)
        self.assertEqual(float(plume[:, :40].max()), 0.0)
        self.assertLess(plume[20, 60], plume[20, 41])


class SkyTests(unittest.TestCase):
    def test_level_ground_sees_the_whole_sky_and_a_pit_less(self):
        flat = np.zeros((360, 720), np.float32)
        self.assertTrue(np.allclose(A.sky_view(flat, 2), 1.0, atol=1e-5))
        pit = flat.copy()
        pit[180, 360] = -10000.0
        seen = A.sky_view(pit, 2)
        self.assertLess(seen[180, 360], 0.9)
        self.assertTrue(np.allclose(seen[:100], 1.0, atol=1e-5))


class SeaTests(unittest.TestCase):
    def test_open_water_keeps_its_stress_and_land_carries_a_water_wind_back_to_itself(self):
        from illumination.water_surface import short_waves as S
        gravity, closures = A._sea_closures()
        ustar = np.array([0.094, 0.06])
        u10 = np.array([2.0, float(S.neutral_wind(0.06, gravity, *closures))])    # the land cell's air acts as water's
        wind = dict(ustar_day_ms=np.stack([ustar] * 3), ustar_all_ms=ustar, wind10_day_ms=np.stack([u10] * 3),
                    wind10_all_ms=u10, lsm=np.array([0.0, 1.0]))
        out = A.water_friction_velocity(wind)
        self.assertAlmostEqual(out[0], 0.094)                     # open water: the stress's own
        self.assertAlmostEqual(out[1], 0.06, delta=1e-4)          # the blend over water gives the same air back

    def test_the_sea_steepens_with_the_wind_and_rough_water_mirrors_less_sky_at_grazing(self):
        table = A.sea_slope_table()
        self.assertTrue(np.all(np.diff(table) >= 0))
        self.assertTrue(0.015 < np.interp(np.log(0.094), np.log(A.SEA_USTAR), table) < 0.035)
        sky, index = A.sky_reflectance_table()
        self.assertTrue(np.allclose(sky[:, -1], 0.0212, atol=5e-4))          # straight down, level water's Fresnel
        self.assertTrue(np.all(np.diff(sky[:, 0]) < 0))                      # at grazing, less as the slopes grow
        # and more toward grazing, apart from the quadrature's wiggles of a few thousandths there
        self.assertTrue(np.all(sky[:, 0] > sky[:, 10]) and np.all(sky[:, 10] > sky[:, 20]) and np.all(sky[:, 20] > sky[:, -1]))


if __name__ == '__main__':
    unittest.main()
