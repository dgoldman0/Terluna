"""How the close views are lit: the page's light worked on a tile outside the browser, and the light proposed for the
camera close up.

The page lights the ground with the one Sun through the air, or, fully lit, with a light every place has alike: 45
degrees up from the screen's upper left, the direct sunlight and the diffuse skylight of that Sun height from the
page's light table, the skylight cut to the share of the sky each place sees. shade() works the same sum on a tile's
layers (colour as reflectance, shading normals, water and sky) with north up, without the air between the ground and
the eye, which the page peels away as the camera comes close. It is the page's light re-worked in Python for
comparing tiles; the browser remains the test of the page.

Two tone rules turn that light into a picture:

- tone_page: the page's own, the exposure in stops (one by default), a soft shoulder above 0.8 and the sRGB curve;
- tone_close: the close view's as proposed (CLOSE): the light divided, colour by colour, by what level white ground
  would send back, so that a pixel holds the ground's reflectance times its relief shading whatever the Sun's height
  and whatever colour the air has given the light, times one fixed gain, with the same shoulder. Under this air the
  light on level ground at 45 degrees is orange (red 0.53, green 0.35, blue 0.24 of the sunlight above the air), and
  the page's own tone shows the ground in that light; a satellite mosaic is corrected to the ground's reflectance,
  which is what this division does. A display rule.

The close view as proposed also leaves the sky's mirror image off the water (CLOSE['sky_mirror']), so that water
shows its own colour: under this air the sky is bright, and its image on the water (2% of the skylight) is several
times the deep water's own light, which makes every lake one flat navy. The relief is drawn at its true strength
(CLOSE['relief'] is 1; shade's `relief` can exaggerate it). Under this air about half the light on level ground is
skylight, which no slope turns away from, so the true relief shades gently; drawn twice as steep, as first proposed,
the craters' shading was the first thing the eye found in dry country, where a mosaic is read by its ground.

night_blend() gives the weight of the even light where the Sun is low, as the patch to the page proposes.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

EVEN_ELEVATION, EVEN_AZIMUTH = 45.0, 315.0      # the page's fully lit light: 45 degrees up, from the upper left
PAGE_STOPS = 1.0                                # the page's opening exposure
CLOSE = dict(gain=1.6,                          # times the ground's reflectance, before the shoulder: bare highland soil
             #                                    then stands near 190 of 255, as desert sand does in the reference mosaics
             shoulder=0.8,                      # the page's shoulder, kept
             relief=1.0,                        # the relief's slopes drawn this many times as steep: 1, their true strength
             sky_mirror=0.0,                    # the share of the sky's mirror image left on the water: none, so that
             #                                    water shows its own colour, as the references' near-black water does
             sun_low=(6.0, 24.0),               # the even light comes in as the Sun sinks from 24 to 6 degrees up
             zoom=(3.0, 2.2))                   # and only as the camera comes from 3.0 to 2.2 Moon radii of the centre
WATER_SKY = 0.021                               # still water's reflectance of the sky seen from above (Fresnel, 550 nm)


def page_parameters(page_html):
    """The built page's light: its table of direct and diffuse light on level ground against the Sun's height (both
    relative to sunlight above the air), the Sun's colour and the relief's strength in the shading normals."""
    text = Path(page_html).read_text(encoding='utf-8')
    data = json.loads(re.search(r'<script[^>]*id="atlas-data"[^>]*>(.*?)</script>', text, re.S).group(1))
    p = data['appearance']
    return dict(light=np.asarray(p['light'], float), lut=tuple(p['lutElevations']), sun_colour=np.asarray(p['sunColour'], float),
                normal_strength=float(p['normalStrength']), tau0=np.asarray(p['tau0'], float))


def table(p, row, elevation):
    """The page's lightTable: a row of the light table at a Sun height (degrees), linear between its entries."""
    light = p['light'][row]
    lo, hi = p['lut'][0], p['lut'][-1]
    x = np.clip((elevation - lo) / (hi - lo), 0.0, 1.0) * (len(light) - 1)
    i = int(np.floor(x))
    j = min(i + 1, len(light) - 1)
    return light[i] * (1.0 - (x - i)) + light[j] * (x - i)


def level_light(p, elevation=EVEN_ELEVATION):
    """Direct sunlight per unit of the cosine at the ground, and the skylight on level ground, at a Sun height."""
    direct = table(p, 0, elevation) / max(np.sin(np.radians(elevation)), 1e-3) if elevation > 0 else np.zeros(3)
    return direct, table(p, 1, elevation)


def decode_srgb(c):
    c = np.asarray(c, np.float32) / 255.0
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4).astype(np.float32)


def encode_srgb(c):
    c = np.clip(np.asarray(c, np.float32), 0.0, 1.0)
    return np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1.0 / 2.4) - 0.055)


def shade(colour, normal, water, sky_seen, p, elevation=EVEN_ELEVATION, azimuth=EVEN_AZIMUTH, relief=1.0, mirror=1.0):
    """The light a tile sends up, relative to sunlight above the air (linear, a trailing axis of 3): its reflectance
    (colour, uint8 sRGB) under direct light from `azimuth` (degrees east of north) `elevation` up on the relief its
    shading normals hold (uint8, x east and y north, at the page's strength; `relief` exaggerates it), plus the
    skylight in the share of the sky each place sees; water lies level and mirrors the sky (`mirror` the share of
    that mirror image kept). As the page's shader sums it with the air between ground and eye taken away."""
    albedo = decode_srgb(colour)
    nm = np.asarray(normal, np.float32) / 255.0 * 2.0 - 1.0
    w = np.asarray(water, np.float32)[..., None]
    k = relief / p['normal_strength'] * (1.0 - w[..., 0]) / np.maximum(nm[..., 2], 0.2)
    east, north = nm[..., 0] * k, nm[..., 1] * k
    norm = np.sqrt(east * east + north * north + 1.0)
    el, az = np.radians(elevation), np.radians(azimuth)
    mus = np.clip((east * np.cos(el) * np.sin(az) + north * np.cos(el) * np.cos(az) + np.sin(el)) / norm, 0.0, None)
    direct, sky = level_light(p, elevation)
    seen = (np.asarray(sky_seen, np.float32) * (1.0 - w[..., 0]) + w[..., 0])[..., None]
    light = direct * mus[..., None] + sky * seen
    return ((albedo * light + w * WATER_SKY * mirror * sky * seen) * p['sun_colour']).astype(np.float32)


def white_level(p, elevation=EVEN_ELEVATION):
    """What level white ground under an open sky sends up in that light."""
    direct, sky = level_light(p, elevation)
    return (direct * np.sin(np.radians(max(elevation, 0.0))) + sky) * p['sun_colour']


def shoulder(c, knee=0.8):
    c = np.asarray(c, np.float32)
    return np.where(c > knee, knee + (1.0 - knee) * (1.0 - np.exp(-(c - knee) / (1.0 - knee))), c)


def tone_page(c, stops=PAGE_STOPS):
    """The page's picture of a light: the exposure, its shoulder above 0.8 and the sRGB curve, as uint8."""
    return np.round(encode_srgb(shoulder(c * 2.0 ** stops)) * 255).astype(np.uint8)


def tone_close(c, p, elevation=EVEN_ELEVATION, gain=None):
    """The close view's picture as proposed: the light over level white ground's, colour by colour, times
    CLOSE['gain'], the shoulder and the sRGB curve."""
    white = np.maximum(white_level(p, elevation), 1e-4).astype(np.float32)
    return np.round(encode_srgb(shoulder(c / white * (CLOSE['gain'] if gain is None else gain), CLOSE['shoulder'])) * 255).astype(np.uint8)


def smoothstep(a, b, x):
    t = np.clip((np.asarray(x, float) - a) / (b - a), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def close_weight(distance):
    """How far the close view's light and tone have come in at a camera distance (Moon radii from the centre)."""
    return float(smoothstep(CLOSE['zoom'][0], CLOSE['zoom'][1], distance))


def night_blend(sun_elevation, distance):
    """The even light's weight in the ground's light: none far out or under a high Sun, all of it close up where the
    Sun stands under 6 degrees, so that a region in the Moon's night still shows its ground."""
    return close_weight(distance) * float(smoothstep(CLOSE['sun_low'][1], CLOSE['sun_low'][0], sun_elevation))


def lit(colour, normal, water, sky_seen, p, sun_elevation, distance, sun_azimuth=90.0):
    """The ground's light as the patched page would give it at a Sun height and camera distance: the Sun's own light
    and the even light mixed by night_blend, and the page's tone and the close view's mixed by close_weight. Returns
    the picture as uint8."""
    k, n = close_weight(distance), night_blend(sun_elevation, distance)
    steep = 1.0 + (CLOSE['relief'] - 1.0) * k
    mirror = 1.0 - k * (1.0 - CLOSE['sky_mirror'])
    own = shade(colour, normal, water, sky_seen, p, sun_elevation, sun_azimuth, relief=steep, mirror=mirror) if sun_elevation > -18 else 0.0
    even = shade(colour, normal, water, sky_seen, p, relief=steep, mirror=mirror)
    c = own * (1.0 - n) + even * n
    white = np.maximum(white_level(p, sun_elevation) * (1.0 - n) + white_level(p) * n, 1e-4).astype(np.float32)
    return np.round(encode_srgb(shoulder(c * (2.0 ** PAGE_STOPS * (1.0 - k) + CLOSE['gain'] / white * k))) * 255).astype(np.uint8)
