"""Mean-orbit earthlight context using the existing solar-coloured approximation."""
from __future__ import annotations
import numpy as np
from illumination.ephemeris import Site,state
from shared.constants import SYNODIC_MONTH_DAYS


def context(latitude_deg,longitude_deg,hours_after_sunset,clear_sky_samples):
    """Screen Earth phase and clear-ground light separately from the solar scene.

    Earth is a Lambert-phase sphere at mean distance, with libration disabled.
    Scaling the solar sky assumes the same source colour and uses the smaller
    solar disk as a point-source approximation; the returned light is a proxy.
    """
    site=Site(latitude_deg,longitude_deg,libration=False)
    seconds=(SYNODIC_MONTH_DAYS*24/4+hours_after_sunset)*3600
    earth=state(site,seconds)
    elevation=float(np.degrees(np.arcsin(earth['earth_enu'][2])))
    proxy=None
    if elevation>0:
        angles=np.array([r['sun_deg'] for r in clear_sky_samples])
        light=np.array([r['total_horizontal_lux'] for r in clear_sky_samples])
        proxy=float(np.exp(np.interp(elevation,angles,np.log(light)))*earth['earthlight_ratio'])
    return dict(earth_elevation_deg=elevation,earth_illuminated_fraction=earth['earth_illuminated_fraction'],
                source_ratio=earth['earthlight_ratio'],clear_surface_proxy_lux=proxy,
                reading_rule='Separate clear-atmosphere, solar-coloured Earth proxy. Cloud attenuation and reflection, '
                             'Earth spectrum, larger angular disk and libration require a separate Earth-source render. '
                             'The surface proxy is evaluated where Earth is above the horizon.')
