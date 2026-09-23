"""Geometry checks for the site ephemeris (mean-orbit model)."""
import math
import unittest

import numpy as np

from illumination.ephemeris import (DAY, EPSILON, Site, celestial_to_local, earth_body, enu_basis,
                                    lambert_phase, state, sun_body)
from shared import constants as K

SYNODIC = K.SYNODIC_MONTH_DAYS * DAY
SITE = Site(0.0, -65.0, 90.0)
FIXED = Site(0.0, -65.0, 90.0, libration=False)


class EphemerisTests(unittest.TestCase):
    def test_noon_zenith_and_synodic_return(self):
        np.testing.assert_allclose(state(SITE, 0)["sun_enu"], [0, 0, 1], atol=1e-12)
        np.testing.assert_allclose(state(SITE, SYNODIC)["sun_enu"], [0, 0, 1], atol=1e-9)

    def test_sun_sets_in_the_west(self):
        np.testing.assert_allclose(state(SITE, SYNODIC / 4)["sun_enu"], [-1, 0, 0], atol=1e-9)

    def test_earth_fixed_without_libration(self):
        for t in (0, 3 * DAY, SYNODIC / 2):
            e = state(FIXED, t)["earth_enu"]
            self.assertAlmostEqual(math.degrees(math.asin(e[2])), 25.0, places=9)
            self.assertAlmostEqual(e[1], 0.0, places=12)
            self.assertGreater(e[0], 0)  # east

    def test_libration_is_bounded(self):
        limit = math.radians(math.hypot(K.LIBRATION_LONGITUDE_DEG, K.LIBRATION_LATITUDE_DEG)) + 1e-9
        for t in np.linspace(0, 3 * SYNODIC, 400):
            self.assertLessEqual(math.acos(np.clip(earth_body(SITE, t) @ [1, 0, 0], -1, 1)), limit)

    def test_phase_and_earthlight(self):
        st = state(FIXED, SYNODIC * 115 / 360)  # sub-solar longitude 180: Sun opposite the sub-Earth point
        self.assertAlmostEqual(st["earth_illuminated_fraction"], 1.0, places=9)
        full = K.EARTH_GEOMETRIC_ALBEDO * (K.EARTH_RADIUS / K.EARTH_MOON_DISTANCE) ** 2
        self.assertAlmostEqual(st["earthlight_ratio"] / full, 1.0, places=9)
        self.assertAlmostEqual(full, 1.008e-4, delta=0.002e-4)
        self.assertAlmostEqual(lambert_phase(math.pi), 0.0, places=12)
        for t in np.linspace(0, SYNODIC, 50):
            st = state(SITE, t)
            self.assertAlmostEqual(st["earth_illuminated_fraction"], (1 + math.cos(st["earth_phase_angle_rad"])) / 2, places=12)

    def test_celestial_frame_is_rotation_with_sidereal_period(self):
        for t in (0, 5 * DAY, 20 * DAY):
            m = celestial_to_local(SITE, t)
            np.testing.assert_allclose(m @ m.T, np.eye(3), atol=1e-12)
            self.assertAlmostEqual(np.linalg.det(m), 1.0, places=12)
            np.testing.assert_allclose(celestial_to_local(SITE, t + K.SIDEREAL_MONTH_DAYS * DAY), m, atol=1e-9)

    def test_ecliptic_pole_is_north_at_an_equatorial_site(self):
        pole = np.array([0, -math.sin(EPSILON), math.cos(EPSILON)])  # ecliptic pole, J2000 equatorial
        np.testing.assert_allclose(celestial_to_local(SITE, 7 * DAY) @ pole, [0, 1, 0], atol=1e-12)

    def test_sun_agrees_with_the_star_frame(self):
        # The Sun's ecliptic longitude advances at (sidereal - synodic) rate from its noon value.
        omega = 2 * math.pi * (1 / K.SIDEREAL_MONTH_DAYS - 1 / K.SYNODIC_MONTH_DAYS) / DAY
        for t in np.linspace(0, 2 * SYNODIC, 17):
            lam = math.radians(SITE.sun_longitude_at_noon_deg) + omega * t
            ecl = np.array([math.cos(lam), math.sin(lam), 0.0])
            c, s = math.cos(EPSILON), math.sin(EPSILON)
            eq = np.array([ecl[0], c * ecl[1] - s * ecl[2], s * ecl[1] + c * ecl[2]])  # ecliptic -> equatorial
            np.testing.assert_allclose(celestial_to_local(SITE, t) @ eq, enu_basis(SITE) @ sun_body(SITE, t), atol=1e-9)


if __name__ == "__main__":
    unittest.main()
