"""The bright stars over a frame, as a separate layer of point images.

Positions: the Yale Bright Star Catalogue's places in the Moon's sky at the frame's date (illumination/stars/sky.py).
Brightness: V magnitude, with V = 0 giving 2.54e-6 lux above the air; colour from B-V through a blackbody
(Ballesteros 2012), per channel of the solved sky relative to the shielded sunlight it carries (the shield does
not dim the stars), dimmed by the solved air's direct transmittance at the star's elevation. Each star is drawn
as a small Gaussian spot holding its illuminance; terrain and the sea horizon hide it. Reflections of stars in
calm water are left out.
"""
from __future__ import annotations

import json
import math

import numpy as np

import scene
from illumination.stars.sky import V0_LUX, body_directions
from research.studies.sea_appearance import lighting

PLANCK_C2 = 1.4388e-2          # m K
SPOT_SIGMA_PX = 0.6
FAINTEST = 6.5


def blackbody(wavelength_nm, temperature):
    lam = wavelength_nm * 1e-9
    return lam ** -5 / np.expm1(PLANCK_C2 / (lam * temperature))


def colour_temperature(b_minus_v):
    """Ballesteros (2012, EPL 97, 34008)."""
    return 4600 * (1 / (0.92 * b_minus_v + 1.7) + 1 / (0.92 * b_minus_v + 0.62))


def horizon_profile(arrays, step_deg=0.1):
    """The terrain's highest elevation angle (degrees) seen from the camera, by azimuth."""
    p = arrays["params"]
    h, x0, y0, dx, offset = arrays["ter_h"], p[scene.k.P_TER_X0], p[scene.k.P_TER_Y0], p[scene.k.P_TER_DX], \
        p[scene.k.P_TER_OFFSET]
    azimuths = np.arange(0.0, 360.0, step_deg)
    dist = np.arange(50.0, 1.2e5, 50.0)
    out = np.full(len(azimuths), -90.0)
    for i, az in enumerate(np.radians(azimuths)):
        x, y = dist * math.sin(az), dist * math.cos(az)
        u, v = (x - x0) / dx, (y - y0) / dx
        inside = (u >= 0) & (v >= 0) & (u < h.shape[1] - 1) & (v < h.shape[0] - 1)
        if not inside.any():
            continue
        ui, vi = u[inside].astype(int), v[inside].astype(int)
        z = h[vi, ui] + offset - dist[inside] ** 2 / (2 * scene.MOON_RADIUS) - scene.EYE_HEIGHT_M
        out[i] = np.degrees(np.max(np.arctan2(z, dist[inside])))
    return azimuths, out


def layer(common, record, geometry, camera, arrays, width, height, fov_deg):
    catalogue = json.loads(scene.STARS.read_text())
    stars = np.array(catalogue["stars"], float)
    stars = stars[stars[:, 2] <= FAINTEST]
    body = body_directions(stars[:, 0], stars[:, 1], record["date_tt"])
    lon, lat = record["camera_lon_lat"]
    east, north, up = lighting.local_frame(lon, lat)
    local = np.stack([body @ east, body @ north, body @ up], -1)
    elevation = np.degrees(np.arcsin(np.clip(local[:, 2], -1, 1)))
    azimuth = np.degrees(np.arctan2(local[:, 0], local[:, 1])) % 360
    profile_az, profile = horizon_profile(arrays)
    blocked = elevation <= np.interp(azimuth, profile_az, profile, period=360)
    dip = math.degrees(math.sqrt(2 * scene.EYE_HEIGHT_M / scene.MOON_RADIUS))
    visible = (~blocked) & (elevation > -dip + 0.05)
    # Colour per channel relative to the shielded sunlight, scaled so that V gives the illuminance above the air.
    sky = common.sky
    bands = np.unique(sky.channel_band, axis=0)
    band_index = np.array([int(np.flatnonzero((bands == b).all(1))[0]) for b in sky.channel_band])
    sun_band = np.array([sky.channel_energy[band_index == i].sum() for i in range(len(bands))])
    centre = bands.mean(axis=1)
    width_nm = np.diff(bands, axis=1)[:, 0]
    xyz = common.xyz
    image = np.zeros((height, width, 3))
    f, right, upv = camera[3:6], camera[6:9], camera[9:12]
    th = math.tan(math.radians(fov_deg) / 2)
    aspect = width / height
    focal = width / 2 / th                     # a rectilinear lens: a pixel at angle a off the axis spans cos^3 a / f^2
    drawn, brightest = 0, []
    for i in np.flatnonzero(visible):
        d = local[i]
        z = d @ f
        if z <= 0:
            continue
        u = (d @ right / z / th + 1) / 2 * width
        v = (1 - d @ upv / z / th * aspect) / 2 * height
        if not (-3 <= u < width + 3 and -3 <= v < height + 3):
            continue
        ratio = blackbody(centre, colour_temperature(stars[i, 3])) * width_nm / sun_band
        per_channel = ratio[band_index]
        above = per_channel @ xyz
        transmit = np.array([np.interp(elevation[i], sky.beam_angles, sky.direct_channels[:, c])
                             for c in range(sky.direct_channels.shape[1])]) / math.sin(math.radians(max(elevation[i], 0.5)))
        seen = (per_channel * np.clip(transmit, 0, 1)) @ xyz
        e = V0_LUX * 10 ** (-0.4 * stars[i, 2]) * seen / above[1]
        # A Gaussian spot holding the illuminance: radiance = E * spot / the solid angle of a pixel there.
        omega = z ** 3 / focal ** 2
        r = 3
        cu, cv = int(math.floor(u)), int(math.floor(v))
        ys, xs = np.mgrid[cv - r:cv + r + 1, cu - r:cu + r + 1]
        w = np.exp(-(((xs + 0.5) - u) ** 2 + ((ys + 0.5) - v) ** 2) / (2 * SPOT_SIGMA_PX ** 2))
        w /= w.sum()
        ok = (xs >= 0) & (ys >= 0) & (xs < width) & (ys < height)
        image[ys[ok], xs[ok]] += w[ok][:, None] * e / omega
        drawn += 1
        brightest.append((float(stars[i, 2]), round(float(e[1]), 9)))
    brightest.sort()
    return image, dict(catalogue=catalogue["source"]["catalogue"], faintest_v=FAINTEST, drawn=drawn,
                       above_terrain_and_sea=int(visible.sum()), brightest_v_and_lux=brightest[:5],
                       note="point images of 0.6-pixel Gaussian spots; no reflections in the water")
