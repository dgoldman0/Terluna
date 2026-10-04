"""Where the catalogue's stars stand in the Moon's sky, and which an eye picks out against a bright sky.

Positions: the J2000 places of bright_stars.json, turned to the ecliptic, carried to the equinox of date by
general precession in longitude (5029.0966 arcseconds per Julian century), and placed in the Moon's body frame
by the optical-libration transform that places the Earth and the Sun (geography/lunar_ephemeris.py, Meeus
chapter 53). Proper motion is not applied, as the catalogue notes; parallax and aberration are left out.

The Moon spins about an axis 1.5 degrees from the ecliptic pole, so its sky turns once a sidereal month about a
point in Draco in the north and in Dorado in the south.

Visibility: the naked-eye limiting magnitude against a sky of luminance L (cd/m2) follows Schaefer (1990,
PASP 102, 212), m = 7.93 - 5 log10(10^(4.316 - b/5) + 1), with the sky's surface brightness
b = 12.58 - 2.5 log10(L) in V magnitudes per square arcsecond.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from geography import lunar_ephemeris as ephemeris

CATALOGUE = Path(__file__).with_name("bright_stars.json")
OBLIQUITY_J2000 = math.radians(23.4392911)
PRECESSION_ARCSEC_PER_CENTURY = 5029.0966
V0_LUX = 2.54e-6                      # illuminance of a V = 0 star above the air


def catalogue(faintest=6.5):
    """Rows of J2000 right ascension and declination (rad), V and B-V, down to a magnitude."""
    stars = np.array(json.loads(CATALOGUE.read_text())["stars"], float)
    return stars[stars[:, 2] <= faintest]


def body_directions(ra, dec, jd_tt):
    """Unit vectors (N, 3) toward stars at J2000 places, in the Moon's body frame at a date."""
    ra, dec = np.atleast_1d(np.asarray(ra, float)), np.atleast_1d(np.asarray(dec, float))
    x, y, z = np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)
    c, s = math.cos(OBLIQUITY_J2000), math.sin(OBLIQUITY_J2000)
    ey, ez = c * y + s * z, -s * y + c * z
    lon = np.arctan2(ey, x) + math.radians(PRECESSION_ARCSEC_PER_CENTURY * (jd_tt - ephemeris.J2000) / 36525 / 3600)
    lat = np.arcsin(np.clip(ez, -1, 1))
    _, _, _, elements = ephemeris.moon_geocentric(np.array([jd_tt]))
    # selenographic() takes the direction from a body to the Moon: the opposite of the star's direction.
    sel_lon, sel_lat = ephemeris.selenographic(lon + math.pi, -lat, elements["node"][0], elements["F"][0])
    return ephemeris.unit(sel_lon, sel_lat)


def naked_eye_limit(sky_luminance_cd_m2):
    """Faintest V magnitude an eye picks out against a sky of this luminance (Schaefer 1990)."""
    b = 12.58 - 2.5 * np.log10(np.maximum(sky_luminance_cd_m2, 1e-12))
    return 7.93 - 5 * np.log10(10 ** (4.316 - b / 5) + 1)
