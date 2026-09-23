"""Sun, Earth and star geometry above a site on the Moon, and the earthlight it sets.

A mean-orbit model, good to a few degrees and not an ephemeris for dates:

- The Moon rotates synchronously; its equator is taken in the ecliptic (the 1.54 degree
  tilt is neglected), so the Sun crosses the zenith of an equatorial site.
- The Earth's direction in the Moon's frame wobbles by optical libration, modelled as
  sinusoids with the anomalistic and draconic months.
- The Sun's ecliptic longitude advances with the tropical-free "year" implied by the
  sidereal and synodic months, so the Sun returns to local noon every synodic month.
- Earthlight: the Earth is a Lambert-phase sphere with the visible geometric albedo in
  shared/constants.json at the mean Earth-Moon distance.

Frames. Body: z along the Moon's spin axis, x toward the mean Earth. Local: east, north,
up at the site. Celestial: J2000 equatorial. Time t is seconds after local noon.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from shared import constants as K

DAY = 86400.0
TAU = 2 * math.pi
OMEGA_SYNODIC = TAU / (K.SYNODIC_MONTH_DAYS * DAY)
OMEGA_SIDEREAL = TAU / (K.SIDEREAL_MONTH_DAYS * DAY)
OMEGA_ANOMALISTIC = TAU / (K.ANOMALISTIC_MONTH_DAYS * DAY)
OMEGA_DRACONIC = TAU / (K.DRACONIC_MONTH_DAYS * DAY)
OMEGA_EARTH = TAU / K.EARTH_SIDEREAL_DAY_S
EPSILON = math.radians(K.EARTH_OBLIQUITY_DEG)
EARTH_ANGULAR_RADIUS = math.asin(K.EARTH_RADIUS / K.EARTH_MOON_DISTANCE)


@dataclass(frozen=True)
class Site:
    latitude_deg: float
    longitude_deg: float
    sun_longitude_at_noon_deg: float = 90.0
    earth_rotation_at_noon_deg: float = 0.0
    libration_phase_longitude_deg: float = 0.0
    libration_phase_latitude_deg: float = 0.0
    libration: bool = True

    @classmethod
    def from_scenario(cls, spec: dict, libration: bool = True) -> "Site":
        return cls(spec["selenographic_latitude_deg"], spec["selenographic_longitude_deg"],
                   spec["sun_ecliptic_longitude_at_local_noon_deg"],
                   spec["earth_rotation_angle_at_local_noon_deg"],
                   spec["libration_phase_longitude_deg"], spec["libration_phase_latitude_deg"],
                   libration)


def unit(lat: float, lon: float) -> np.ndarray:
    return np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)])


def enu_basis(site: Site) -> np.ndarray:
    """Rows are east, north and up at the site, in the body frame."""
    lat, lon = math.radians(site.latitude_deg), math.radians(site.longitude_deg)
    east = np.array([-math.sin(lon), math.cos(lon), 0.0])
    north = np.array([-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat)])
    return np.vstack([east, north, unit(lat, lon)])


def sun_body(site: Site, t: float) -> np.ndarray:
    """Sun direction in the body frame; the sub-solar point moves west at the synodic rate."""
    return unit(0.0, math.radians(site.longitude_deg) - OMEGA_SYNODIC * t)


def earth_body(site: Site, t: float) -> np.ndarray:
    """Earth direction in the body frame: the sub-Earth point, displaced by libration."""
    if not site.libration:
        return unit(0.0, 0.0)
    lon = math.radians(K.LIBRATION_LONGITUDE_DEG) * math.sin(
        OMEGA_ANOMALISTIC * t + math.radians(site.libration_phase_longitude_deg))
    lat = math.radians(K.LIBRATION_LATITUDE_DEG) * math.sin(
        OMEGA_DRACONIC * t + math.radians(site.libration_phase_latitude_deg))
    return unit(lat, lon)


def body_rotation(site: Site, t: float) -> float:
    """Ecliptic longitude of the body x axis (the mean Earth direction seen from the Moon)."""
    moon_longitude_at_noon = math.radians(site.sun_longitude_at_noon_deg) - math.pi - math.radians(site.longitude_deg)
    return moon_longitude_at_noon + math.pi + OMEGA_SIDEREAL * t


def celestial_to_local(site: Site, t: float) -> np.ndarray:
    """Matrix taking J2000 equatorial unit vectors to local east-north-up."""
    c, s = math.cos(EPSILON), math.sin(EPSILON)
    eq_to_ecl = np.array([[1, 0, 0], [0, c, s], [0, -s, c]])
    theta = body_rotation(site, t)
    ecl_to_body = np.array([[math.cos(theta), math.sin(theta), 0], [-math.sin(theta), math.cos(theta), 0], [0, 0, 1]])
    return enu_basis(site) @ ecl_to_body @ eq_to_ecl


def earth_rotation_angle(site: Site, t: float) -> float:
    """Angle of the Earth's prime meridian about the celestial pole."""
    return math.radians(site.earth_rotation_at_noon_deg) + OMEGA_EARTH * t


def lambert_phase(alpha: float) -> float:
    return (math.sin(alpha) + (math.pi - alpha) * math.cos(alpha)) / math.pi


def state(site: Site, t: float) -> dict:
    """Directions (local east-north-up), Earth's phase and the earthlight ratio at time t."""
    basis = enu_basis(site)
    sun, earth = sun_body(site, t), earth_body(site, t)
    cos_se = float(np.clip(sun @ earth, -1.0, 1.0))
    alpha = math.acos(-cos_se)  # Sun-Earth-Moon angle: 0 at full Earth
    ratio = K.EARTH_GEOMETRIC_ALBEDO * (K.EARTH_RADIUS / K.EARTH_MOON_DISTANCE) ** 2 * lambert_phase(alpha)
    return {
        "t_s": t,
        "sun_enu": (basis @ sun).tolist(),
        "earth_enu": (basis @ earth).tolist(),
        "earth_illuminated_fraction": (1.0 - cos_se) / 2.0,
        "earth_phase_angle_rad": alpha,
        "earthlight_ratio": ratio,
        "celestial_to_enu": celestial_to_local(site, t).tolist(),
        "earth_rotation_rad": earth_rotation_angle(site, t),
    }


def product(site_id: str, site: Site, samples: int = 97) -> dict:
    """The site-sky product: model parameters and golden samples over one synodic month."""
    period = K.SYNODIC_MONTH_DAYS * DAY
    return {
        "schema": "terluna.illumination.site-sky/1",
        "producer": {"domain": "illumination", "model": "illumination/ephemeris.py"},
        "evidence": ("Mean-orbit geometry: synchronous rotation, lunar equator in the ecliptic, "
                     "sinusoidal optical libration, Lambert-phase Earth at mean distance. Good to a few "
                     "degrees; not an ephemeris for dates."),
        "site_id": site_id,
        "site": site.__dict__,
        "constants": {
            "synodic_month_s": period,
            "sidereal_month_s": K.SIDEREAL_MONTH_DAYS * DAY,
            "anomalistic_month_s": K.ANOMALISTIC_MONTH_DAYS * DAY,
            "draconic_month_s": K.DRACONIC_MONTH_DAYS * DAY,
            "earth_sidereal_day_s": K.EARTH_SIDEREAL_DAY_S,
            "obliquity_rad": EPSILON,
            "libration_longitude_rad": math.radians(K.LIBRATION_LONGITUDE_DEG),
            "libration_latitude_rad": math.radians(K.LIBRATION_LATITUDE_DEG),
            "earth_geometric_albedo": K.EARTH_GEOMETRIC_ALBEDO,
            "earth_radius_m": K.EARTH_RADIUS,
            "earth_moon_distance_m": K.EARTH_MOON_DISTANCE,
            "earth_angular_radius_rad": EARTH_ANGULAR_RADIUS,
        },
        # Golden samples span 1.37 synodic months, so consumers are checked across the wrap.
        "samples": [state(site, period * 1.37 * i / (samples - 1)) for i in range(samples)],
    }
