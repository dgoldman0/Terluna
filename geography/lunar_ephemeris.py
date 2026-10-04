"""Earth and Sun seen from the Moon: direction in the Moon's body frame and distance.

The tide product needs where the Earth and the Sun stand over the Moon through
the months and the 18.6-year nodal cycle. The geocentric Moon comes from the
truncated ELP-2000/82 lunar theory (Chapront-Touze and Chapront), as tabulated
in Meeus, Astronomical Algorithms (2nd ed., 1998), chapter 47: mean elements,
the leading periodic terms in longitude and distance (Table 47.A) and in
latitude (Table 47.B). The Earth's selenographic position is the optical
libration of chapter 53, with the lunar equator inclined 1.54242 degrees to
the ecliptic (Cassini's laws); the physical librations, under 0.05 degrees,
are left out. The Sun's geocentric position is the low-precision solar theory
of chapter 25.

Frames: the Moon's mean-Earth/polar-axis frame, x toward longitude 0 and
latitude 0, z along the spin axis, east longitudes positive. The check against
JPL Horizons (tides_ephemeris_check.csv) bounds the errors; tides.py records it.
"""
from __future__ import annotations

import math

import numpy as np

J2000 = 2451545.0
EQUATOR_TO_ECLIPTIC = math.radians(1.54242)  # Meeus 53, the IAU value used with these librations
MEAN_DISTANCE_KM = 385000.56
AU_KM = 149597870.7

# Table 47.A: multiples of D, M, M', F; longitude (1e-6 degree); distance (1e-3 km).
LONGITUDE_DISTANCE = np.array([
    (0, 0, 1, 0, 6288774, -20905355), (2, 0, -1, 0, 1274027, -3699111), (2, 0, 0, 0, 658314, -2955968),
    (0, 0, 2, 0, 213618, -569925), (0, 1, 0, 0, -185116, 48888), (0, 0, 0, 2, -114332, -3149),
    (2, 0, -2, 0, 58793, 246158), (2, -1, -1, 0, 57066, -152138), (2, 0, 1, 0, 53322, -170733),
    (2, -1, 0, 0, 45758, -204586), (0, 1, -1, 0, -40923, -129620), (1, 0, 0, 0, -34720, 108743),
    (0, 1, 1, 0, -30383, 104755), (2, 0, 0, -2, 15327, 10321), (0, 0, 1, 2, -12528, 0),
    (0, 0, 1, -2, 10980, 79661), (4, 0, -1, 0, 10675, -34782), (0, 0, 3, 0, 10034, -23210),
    (4, 0, -2, 0, 8548, -21636), (2, 1, -1, 0, -7888, 24208), (2, 1, 0, 0, -6766, 30824),
    (1, 0, -1, 0, -5163, -8379), (1, 1, 0, 0, 4987, -16675), (2, -1, 1, 0, 4036, -12831),
    (2, 0, 2, 0, 3994, -10445), (4, 0, 0, 0, 3861, -11650), (2, 0, -3, 0, 3665, 14403),
    (0, 1, -2, 0, -2689, -7003), (2, 0, -1, 2, -2602, 0), (2, -1, -2, 0, 2390, 10056),
    (1, 0, 1, 0, -2348, 6322), (2, -2, 0, 0, 2236, -9884), (0, 1, 2, 0, -2120, 5751),
    (0, 2, 0, 0, -2069, 0), (2, -2, -1, 0, 2048, -4950), (2, 0, 1, -2, -1773, 4130),
    (2, 0, 0, 2, -1595, 0), (4, -1, -1, 0, 1215, -3958), (0, 0, 2, 2, -1110, 0),
    (3, 0, -1, 0, -892, 3258), (2, 1, 1, 0, -810, 2616), (4, -1, -2, 0, 759, -1897),
    (0, 2, -1, 0, -713, -2117), (2, 2, -1, 0, -700, 2354), (2, 1, -2, 0, 691, 0),
    (2, -1, 0, -2, 596, 0), (4, 0, 1, 0, 549, -1423), (0, 0, 4, 0, 537, -1117),
    (4, -1, 0, 0, 520, -1571), (1, 0, -2, 0, -487, -1739), (2, 1, 0, -2, -399, 0),
    (0, 0, 2, -2, -381, -4421), (1, 1, 1, 0, 351, 0), (3, 0, -2, 0, -340, 0),
    (4, 0, -3, 0, 330, 0), (2, -1, 2, 0, 327, 0), (0, 2, 1, 0, -323, 1165),
    (1, 1, -1, 0, 299, 0), (2, 0, 3, 0, 294, 0), (2, 0, -1, -2, 0, 8752),
], dtype=float)

# Table 47.B: multiples of D, M, M', F; latitude (1e-6 degree).
LATITUDE = np.array([
    (0, 0, 0, 1, 5128122), (0, 0, 1, 1, 280602), (0, 0, 1, -1, 277693), (2, 0, 0, -1, 173237),
    (2, 0, -1, 1, 55413), (2, 0, -1, -1, 46271), (2, 0, 0, 1, 32573), (0, 0, 2, 1, 17198),
    (2, 0, 1, -1, 9266), (0, 0, 2, -1, 8822), (2, -1, 0, -1, 8216), (2, 0, -2, -1, 4324),
    (2, 0, 1, 1, 4200), (2, 1, 0, -1, -3359), (2, -1, -1, 1, 2463), (2, -1, 0, 1, 2211),
    (2, -1, -1, -1, 2065), (0, 1, -1, -1, -1870), (4, 0, -1, -1, 1828), (0, 1, 0, 1, -1794),
    (0, 0, 0, 3, -1749), (0, 1, -1, 1, -1565), (1, 0, 0, 1, -1491), (0, 1, 1, 1, -1475),
    (0, 1, 1, -1, -1410), (0, 1, 0, -1, -1344), (1, 0, 0, -1, -1335), (0, 0, 3, 1, 1107),
    (4, 0, 0, -1, 1021), (4, 0, -1, 1, 833),
], dtype=float)


def centuries(jd_tt):
    return (np.asarray(jd_tt, dtype=float) - J2000) / 36525.0


def mean_elements(t):
    """Mean longitude L', elongation D, solar anomaly M, lunar anomaly M', argument F, node, in radians."""
    deg = lambda *c: np.radians(c[0] + c[1] * t + c[2] * t**2 + c[3] * t**3 + c[4] * t**4)
    return dict(
        L=deg(218.3164477, 481267.88123421, -0.0015786, 1 / 538841, -1 / 65194000),
        D=deg(297.8501921, 445267.1114034, -0.0018819, 1 / 545868, -1 / 113065000),
        M=deg(357.5291092, 35999.0502909, -0.0001536, 1 / 24490000, 0.0),
        Mp=deg(134.9633964, 477198.8675055, 0.0087414, 1 / 69699, -1 / 14712000),
        F=deg(93.2720950, 483202.0175233, -0.0036539, -1 / 3526000, 1 / 863310000),
        node=deg(125.0445479, -1934.1362891, 0.0020754, 1 / 467441, -1 / 60616000))


def moon_geocentric(jd_tt):
    """Geocentric ecliptic longitude and latitude (radians, mean equinox of date) and distance (km)."""
    t = centuries(jd_tt)
    e = mean_elements(t)
    eccentricity_factor = 1 - 0.002516 * t - 0.0000074 * t**2
    def series(table, value_column, trig):
        total = 0.0
        for row in table:
            d, m, mp, f = row[:4]
            arg = d * e["D"] + m * e["M"] + mp * e["Mp"] + f * e["F"]
            total = total + row[value_column] * eccentricity_factor**abs(m) * trig(arg)
        return total
    a1 = np.radians(119.75 + 131.849 * t)
    a2 = np.radians(53.09 + 479264.290 * t)
    a3 = np.radians(313.45 + 481266.484 * t)
    sum_l = series(LONGITUDE_DISTANCE, 4, np.sin) + 3958 * np.sin(a1) + 1962 * np.sin(e["L"] - e["F"]) + 318 * np.sin(a2)
    sum_r = series(LONGITUDE_DISTANCE, 5, np.cos)
    sum_b = (series(LATITUDE, 4, np.sin) - 2235 * np.sin(e["L"]) + 382 * np.sin(a3) + 175 * np.sin(a1 - e["F"])
             + 175 * np.sin(a1 + e["F"]) + 127 * np.sin(e["L"] - e["Mp"]) - 115 * np.sin(e["L"] + e["Mp"]))
    longitude = e["L"] + np.radians(sum_l * 1e-6)
    latitude = np.radians(sum_b * 1e-6)
    distance = MEAN_DISTANCE_KM + sum_r / 1000.0
    return longitude, latitude, distance, e


def sun_geocentric(jd_tt):
    """Geocentric ecliptic longitude (radians, mean equinox of date) and distance (km) of the Sun."""
    t = centuries(jd_tt)
    mean_longitude = np.radians(280.46646 + 36000.76983 * t + 0.0003032 * t**2)
    anomaly = np.radians(357.52911 + 35999.05029 * t - 0.0001537 * t**2)
    eccentricity = 0.016708634 - 0.000042037 * t - 0.0000001267 * t**2
    centre = np.radians((1.914602 - 0.004817 * t - 0.000014 * t**2) * np.sin(anomaly)
                        + (0.019993 - 0.000101 * t) * np.sin(2 * anomaly) + 0.000289 * np.sin(3 * anomaly))
    true_anomaly = anomaly + centre
    distance_au = 1.000001018 * (1 - eccentricity**2) / (1 + eccentricity * np.cos(true_anomaly))
    return mean_longitude + centre, distance_au * AU_KM


def selenographic(longitude, latitude, node, argument_f):
    """Selenographic longitude and latitude of the point beneath a body (Meeus 53.1).

    `longitude` and `latitude` give the ecliptic direction from that body to the
    Moon: the geocentric Moon for the sub-Earth point, the heliocentric Moon
    for the sub-solar point.
    """
    w = longitude - node
    i = EQUATOR_TO_ECLIPTIC
    a = np.arctan2(np.sin(w) * np.cos(latitude) * np.cos(i) - np.sin(latitude) * np.sin(i),
                   np.cos(w) * np.cos(latitude))
    lon = np.mod(a - argument_f + np.pi, 2 * np.pi) - np.pi
    lat = np.arcsin(-np.sin(w) * np.cos(latitude) * np.sin(i) - np.sin(latitude) * np.cos(i))
    return lon, lat


def unit(lon, lat):
    return np.stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)), axis=-1)


def earth_and_sun(jd_tt):
    """Unit vectors to the Earth and the Sun in the Moon's body frame, and their distances (m)."""
    moon_lon, moon_lat, moon_km, e = moon_geocentric(jd_tt)
    earth_lon, earth_lat = selenographic(moon_lon, moon_lat, e["node"], e["F"])
    sun_lon, sun_km = sun_geocentric(jd_tt)
    # The heliocentric Moon: the geocentric Moon minus the geocentric Sun.
    moon_vec = unit(moon_lon, moon_lat) * moon_km[..., None]
    sun_vec = unit(sun_lon, np.zeros_like(sun_lon)) * sun_km[..., None]
    rel = moon_vec - sun_vec
    rel_km = np.linalg.norm(rel, axis=-1)
    rel_lon, rel_lat = np.arctan2(rel[..., 1], rel[..., 0]), np.arcsin(rel[..., 2] / rel_km)
    sun_sel_lon, sun_sel_lat = selenographic(rel_lon, rel_lat, e["node"], e["F"])
    return dict(earth=unit(earth_lon, earth_lat), earth_distance_m=moon_km * 1000.0,
                earth_lon_lat_deg=(np.degrees(earth_lon), np.degrees(earth_lat)),
                sun=unit(sun_sel_lon, sun_sel_lat), sun_distance_m=rel_km * 1000.0,
                sun_lon_lat_deg=(np.degrees(sun_sel_lon), np.degrees(sun_sel_lat)))
