"""Checks of the stars' places in the Moon's sky."""
import math

import numpy as np
import pytest

from geography import lunar_ephemeris as ephemeris
from illumination.stars import sky

JD = 2465100.5                                         # February 2038


def test_ecliptic_pole_stands_by_the_moons_spin_axis():
    pole = sky.body_directions(math.radians(270.0), math.radians(66.5607), JD)[0]
    assert math.degrees(math.acos(pole[2])) == pytest.approx(1.54, abs=0.05)


def test_the_sun_placed_as_a_star_matches_the_ephemeris():
    longitude, _ = ephemeris.sun_geocentric(np.array([JD]))
    # Back from the ecliptic of date to J2000 equatorial coordinates.
    lon = longitude[0] - math.radians(sky.PRECESSION_ARCSEC_PER_CENTURY * (JD - ephemeris.J2000) / 36525 / 3600)
    c, s = math.cos(sky.OBLIQUITY_J2000), math.sin(sky.OBLIQUITY_J2000)
    x, y = math.cos(lon), math.sin(lon)
    ra, dec = math.atan2(c * y, x), math.asin(s * y)
    placed = sky.body_directions(ra, dec, JD)[0]
    expected = ephemeris.earth_and_sun(np.array([JD]))["sun"][0]
    assert math.degrees(math.acos(min(1.0, placed @ expected))) < 0.3     # the Moon's parallax on the Sun is 0.15


def test_limiting_magnitude_brightens_with_the_sky():
    dark, twilight, bright = sky.naked_eye_limit(np.array([2.7e-4, 0.05, 5.0]))
    assert dark == pytest.approx(6.4, abs=0.2) and twilight < dark and bright < twilight
