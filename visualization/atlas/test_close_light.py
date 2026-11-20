"""Checks of the close view's light: the page's light table read between its entries, level ground and slopes under
the even light, level water, the close view's tone holding a reflectance whatever the Sun's height, and the even light
coming in close up where the Sun is low."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import close_light as CL  # noqa: E402

ELEVATIONS = np.linspace(-90.0, 90.0, 121)
P = dict(light=np.stack([np.clip(np.sin(np.radians(ELEVATIONS)), 0, None)[:, None] * np.array([0.5, 0.55, 0.6]),
                         np.clip(np.sin(np.radians(ELEVATIONS + 10.0)), 0, None)[:, None] * np.array([0.2, 0.25, 0.3]),
                         np.ones((121, 3))]),
         lut=(-90.0, 90.0), sun_colour=np.array([1.0, 0.98, 0.95]), normal_strength=10.0, tau0=np.zeros(3))


def tile(grey=128, east=0.0, north=0.0, water=0.0, sky=1.0):
    """One pixel: a grey ground whose slope rises by `east` and `north` per metre, as the tiles hold it."""
    n = np.array([-east * 10.0, -north * 10.0, 1.0])
    n /= np.linalg.norm(n)
    return (np.full((1, 1, 3), grey, np.uint8), np.round((n + 1) / 2 * 255).astype(np.uint8)[None, None],
            np.full((1, 1), water, np.float32), np.full((1, 1), sky, np.float32))


def test_the_light_table_is_read_between_its_entries():
    assert np.allclose(CL.table(P, 0, 45.0), np.sin(np.radians(45.0)) * np.array([0.5, 0.55, 0.6]), atol=2e-4)
    direct, sky = CL.level_light(P, 45.0)
    assert np.allclose(direct, [0.5, 0.55, 0.6], atol=5e-4) and np.allclose(sky, CL.table(P, 1, 45.0))


def test_level_ground_sends_back_its_reflectance_of_the_level_light():
    c = CL.shade(*tile(), P)[0, 0]
    assert np.allclose(c, CL.decode_srgb(128) * CL.white_level(P), rtol=2e-3)
    grey = np.round(CL.encode_srgb(CL.shoulder(CL.decode_srgb(128) * CL.CLOSE['gain'])) * 255)
    assert np.allclose(CL.tone_close(CL.shade(*tile(), P), P)[0, 0], grey, atol=1)        # grey ground shows grey: the light's colour is taken out


def test_a_slope_towards_the_light_is_brighter_and_a_pit_darker():
    level = CL.shade(*tile(), P)[0, 0, 1]
    towards = CL.shade(*tile(east=0.2, north=-0.2), P)[0, 0, 1]          # it falls to the north-west, where the light stands
    away = CL.shade(*tile(east=-0.2, north=0.2), P)[0, 0, 1]
    assert towards > level > away
    assert CL.shade(*tile(sky=0.5), P)[0, 0, 1] < level
    assert CL.shade(*tile(east=0.2, north=-0.2), P, relief=2.0)[0, 0, 1] > towards


def test_water_lies_level_and_mirrors_the_sky():
    still = CL.shade(*tile(east=-0.5, north=0.5, water=1.0), P)[0, 0]
    assert np.allclose(still, (CL.decode_srgb(128) * CL.white_level(P) / P['sun_colour'] + CL.WATER_SKY * CL.level_light(P)[1])
                       * P['sun_colour'], rtol=2e-3)
    own = CL.shade(*tile(water=1.0), P, mirror=0.0)[0, 0]        # with the mirror image left off: the water's own colour
    assert np.allclose(own, CL.decode_srgb(128) * CL.white_level(P), rtol=2e-3)


def test_the_close_tone_holds_a_reflectance_whatever_the_suns_height():
    high = CL.tone_close(CL.shade(*tile(), P, elevation=60.0, azimuth=90.0), P, elevation=60.0)[0, 0]
    low = CL.tone_close(CL.shade(*tile(), P, elevation=30.0, azimuth=90.0), P, elevation=30.0)[0, 0]
    assert np.abs(high.astype(int) - low.astype(int)).max() <= 1          # neither the light's strength nor its colour shows
    assert CL.tone_page(CL.shade(*tile(), P, elevation=30.0))[0, 0, 1] < CL.tone_page(CL.shade(*tile(), P, elevation=60.0))[0, 0, 1]


def test_the_even_light_comes_in_close_up_where_the_sun_is_low():
    assert CL.night_blend(-30.0, 5.8) == 0.0 and CL.night_blend(-30.0, 1.5) == 1.0 and CL.night_blend(60.0, 1.5) == 0.0
    assert 0.0 < CL.night_blend(15.0, 1.5) < 1.0 and 0.0 < CL.close_weight(2.6) < 1.0
    night_far = CL.lit(*tile(), P, sun_elevation=-30.0, distance=5.8)[0, 0]
    night_close = CL.lit(*tile(), P, sun_elevation=-30.0, distance=1.5)[0, 0]
    assert night_far.max() == 0 and np.array_equal(night_close, CL.tone_close(CL.shade(*tile(), P), P)[0, 0])   # level ground: no relief to draw steeper
    day_close = CL.lit(*tile(), P, sun_elevation=60.0, distance=1.5, sun_azimuth=90.0)[0, 0]
    assert np.abs(day_close.astype(int) - night_close.astype(int)).max() <= 1
