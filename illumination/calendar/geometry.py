"""Topocentric Sun/Earth geometry in the lunar mean-Earth frame; input dates are TT."""
from __future__ import annotations

import numpy as np
from geography.lunar_ephemeris import earth_and_sun
from shared.constants import AU, EARTH_RADIUS, MOON_RADIUS, SUN_RADIUS


def local_frame(longitude_deg, latitude_deg):
    lon, lat = np.radians([longitude_deg, latitude_deg])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    north = np.array([-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)])
    up = np.array([np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)])
    return np.stack([east, north, up])


def geometry(jd_tt, longitude_deg, latitude_deg):
    if not -90 <= latitude_deg <= 90 or not -180 <= longitude_deg <= 180:
        raise ValueError('Latitude/longitude outside their supported ranges')
    if not np.all(np.isfinite(jd_tt)):
        raise ValueError('Date must be finite')
    g = earth_and_sun(np.asarray(jd_tt, dtype=float))
    return geometry_from_vectors(g,longitude_deg,latitude_deg)


def geometry_from_vectors(g,longitude_deg,latitude_deg):
    """Apply the same observer geometry to an independently supplied ephemeris."""
    frame = local_frame(longitude_deg, latitude_deg)
    observer = MOON_RADIUS*frame[2]
    earth_vector = g['earth']*g['earth_distance_m'][..., None]
    sun_vector = g['sun']*g['sun_distance_m'][..., None]
    earth = (earth_vector-observer)@frame.T
    sun = (sun_vector-observer)@frame.T
    de, ds = np.linalg.norm(earth, axis=-1), np.linalg.norm(sun, axis=-1)
    eu, su = earth/de[..., None], sun/ds[..., None]
    earth_to_sun = sun_vector-earth_vector
    earth_sun_distance = np.linalg.norm(earth_to_sun, axis=-1)
    sun_from_earth = earth_to_sun/earth_sun_distance[..., None]@frame.T
    phase = np.arccos(np.clip(np.sum(sun_from_earth*(-eu), axis=-1), -1, 1))
    # Upward tangent to Earth's disk. At the zenith choose local north.
    vertical = np.array([0.0, 0.0, 1.0])-eu*eu[..., 2, None]
    length = np.linalg.norm(vertical, axis=-1)
    vertical = np.where((length > 1e-12)[..., None], vertical/np.maximum(length[..., None],1e-12),
                        np.array([0.0,1.0,0.0]))
    horizontal = np.cross(vertical, eu)
    beta = np.arctan2(np.sum(sun_from_earth*horizontal, axis=-1),
                      np.sum(sun_from_earth*vertical, axis=-1))
    out = dict(earth_distance_m=de, sun_distance_m=ds,
               earth_phase_deg=np.degrees(phase), earth_lit_fraction=(1+np.cos(phase))/2,
               earth_bright_limb_rad=beta, earth_radius_deg=np.degrees(np.arcsin(EARTH_RADIUS/de)),
               sun_radius_deg=np.degrees(np.arcsin(SUN_RADIUS/ds)),
               sun_distance_factor=(AU/ds)**2, earth_sunlight_factor=(AU/earth_sun_distance)**2,
               source_separation_deg=np.degrees(np.arccos(np.clip(np.sum(eu*su,axis=-1),-1,1))))
    for name, v in [('sun',su),('earth',eu)]:
        out[name+'_elevation_deg']=np.degrees(np.arcsin(np.clip(v[...,2],-1,1)))
        out[name+'_azimuth_deg']=np.degrees(np.arctan2(v[...,0],v[...,1]))%360
    return out
