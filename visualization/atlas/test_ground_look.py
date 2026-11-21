"""Checks of the ground's look: the gathered share read from the moisture map's table, water run down the ground and
conserved, a hollow wetted over the area its water can keep wet, the ground beside water with its marsh, ground that
stands out, and the way a slope faces. The kinds of ground have their own checks in test_ground_kinds.py."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ground_look as GL  # noqa: E402

TABLE = {'biological crust': {'arid': [1.5, 1.8], 'semi-arid': [2.4, 2.8], 'dry sub-humid': [3.5, 3.9]}}


def valley():
    """Five rows by nine columns: ground falling to the middle column from both sides and down it to the top row."""
    r, c = np.mgrid[0:5, 0:9]
    return (np.abs(c - 4) + 0.01 * r).astype(np.float64)


def test_the_gathered_share_is_the_tables():
    s = GL.gathered_share(np.array([0.01, 0.1, 0.316, 0.57, 2.0]), TABLE)
    assert np.allclose(s, [0.0165, 0.0165, 0.026, 0.037, 0.037], atol=1e-4)


def test_water_runs_down_and_none_is_lost():
    s = valley()
    wet = np.zeros(s.shape, bool)
    rcv = GL.receivers(s, wet, 1.0, np.ones(5))
    assert rcv[0, 4] == -1 and rcv[3, 4] == 2 * 9 + 4 and rcv[2, 0] in (1 * 9 + 1, 2 * 9 + 1)
    own, room = np.ones(s.shape), np.full(s.shape, 1e9)
    kept, held = GL.gathered(s, wet, rcv, own, room, 0.0)
    assert kept.sum() == 0.0 and held[0, 4] == 45.0
    kept, held = GL.gathered(s, wet, rcv, own, room, 0.5)
    assert np.isclose(kept.sum() + held.sum(), 45.0) and kept.sum() > 0.0 and held[0, 4] < 45.0
    kept, held = GL.gathered(s, wet, rcv, own, np.zeros(s.shape), 0.5)    # ground with no room takes none
    assert kept.sum() == 0.0 and held[0, 4] == 45.0
    assert (GL.roots(rcv) == 4).all()


def test_water_that_reaches_open_water_is_the_waters():
    s = valley()
    wet = np.zeros(s.shape, bool)
    wet[0] = True
    rcv = GL.receivers(s, wet, 1.0, np.ones(5))
    kept, held = GL.gathered(s, wet, rcv, np.ones(s.shape), np.full(s.shape, 1e9), 0.0)
    assert held.sum() == 0.0 and kept.sum() == 0.0


def test_a_hollow_is_wet_over_the_area_its_water_can_keep_wet():
    s = valley()
    root = np.full(s.size, 4)
    deficit, area = np.ones(s.shape), np.ones(s.shape)
    held = np.zeros(s.shape)
    held[0, 4] = 2.5                                      # enough for two and a half cells
    wetted = GL.floors(s, root, held, deficit, area)
    assert np.isclose(wetted.sum(), 2.5) and wetted[0, 4] == 1.0 and wetted[1, 4] == 1.0 and np.isclose(wetted[2, 4], 0.5)
    held[0, 4] = 1e6                                      # however much water, no further than the hollow's floor
    wetted = GL.floors(s, root, held, deficit, area)
    assert wetted.sum() == (s - s[0, 4] <= GL.LOOK['floor_m']).sum()
    held[0, 4] = 2.5
    assert np.isclose(GL.floors(s, root, held, 5.0 * deficit, area).sum(), 0.5)    # five times the deficit, a fifth the area


def test_ground_beside_water_is_wet_while_low_and_near():
    h = np.repeat(np.arange(12, dtype=np.float32)[None, :], 5, axis=0) * 0.5   # rising half a metre a cell eastward
    wet = np.zeros(h.shape, bool)
    wet[:, 0] = True
    level = np.zeros(h.shape, np.float32)
    none, lat = np.zeros(h.shape, np.float32), np.zeros(5)
    slope = np.full(h.shape, 0.5 / 474.0, np.float32)
    lack = np.full(h.shape, 0.5, np.float32)                 # the ground's deficit, half a metre a year
    w = GL.beside_water(h, wet, level, np.ones(h.shape, bool), none, lat, 474.0, slope, lack)
    assert w[2, 0] == 0.0 and w[2, 1] > 0.99                 # the cell at the water's edge, half a metre up: all of it
    assert 0.4 < w[2, 2] < 0.7 and w[2, 4] == 0.0            # a metre up, two cells off: thinning with the reach; beyond the reach: dry
    steep = GL.beside_water(200.0 * h, wet, level, np.ones(h.shape, bool), none, lat, 474.0, 200.0 * slope, lack)
    assert np.isclose(steep[2, 1], GL.LOOK['shore_m'] / 200.0, atol=0.01)       # a hundred metres up at the edge: a narrow margin
    marsh = GL.beside_water(h, wet, level, np.ones(h.shape, bool), none, lat, 474.0, slope, lack, rise=GL.LOOK['marsh_m'])
    assert marsh[2, 1] > 0.99 and 0.2 < marsh[2, 2] < 0.35 and (marsh <= w + 1e-6).all()     # the marsh: the lowest of the belt
    flat = GL.beside_water(np.zeros_like(h), wet, level, np.ones(h.shape, bool), none, lat, 474.0, 0.0 * slope, lack)
    assert flat[2, 1] > 0.99 and flat[2, 6] == 0.0           # level ground: as far as the reach, 1.5 km
    assert GL.beside_water(h, wet, level, np.zeros(h.shape, bool), none, lat, 474.0, slope, lack).max() == 0.0   # salt water keeps none
    # a stream of a quarter of a m3/s on dry level ground: a hundredth of its flow per km meets the deficit over 79 m to
    # either side, a third of its own cell; where the ground lacks a tenth as much, ten times as far
    flow = np.zeros(h.shape, np.float32)
    flow[2, 6] = 0.25
    dry = np.zeros(h.shape, bool)
    by = GL.beside_water(np.zeros_like(h), dry, level, np.ones(h.shape, bool), flow, lat, 474.0, 0.0 * slope, lack)
    assert np.isclose(by[2, 6], 2.0 * 78.9 / 474.0, atol=0.01) and by[2, 7] == 0.0
    wide = GL.beside_water(np.zeros_like(h), dry, level, np.ones(h.shape, bool), flow, lat, 474.0, 0.0 * slope, 0.1 * lack)
    assert wide[2, 6] == 1.0 and 0.0 < wide[2, 7] < 1.0 and wide[2, 9] == 0.0


def test_ground_that_stands_out_has_less_water_and_a_hollow_no_more():
    y, x = np.mgrid[0:80, 0:80].astype(np.float32)
    ridge = 400.0 * np.exp(-((x - 40.0) / 3.0) ** 2)
    steep = np.abs(np.gradient(ridge, axis=1)) / 474.0
    s = GL.standing(ridge, steep)
    assert s[40, 5] == 0.0 and s[40, 40] > 0.14 and s[40, 40] >= s[40, 42] > s[40, 5] and s.max() <= GL.LOOK['stand_most']
    assert np.exp(-GL.LOOK['stand_most']) > 1.0 - 0.199      # the plants' water cut by no more than a bare loam sheds of its rain
    worn = GL.standing(ridge, 0.2 * steep)                  # the same rise with faces a fifth as steep: it carries cover
    assert worn[40, 40] < 0.6 * s[40, 40]
    assert GL.standing(0.02 * ridge, 0.02 * steep).max() < 0.002             # a swell of eight metres: no cut at all
    bowl = -400.0 * np.exp(-((x - 40.0) ** 2 + (y - 40.0) ** 2) / 36.0)
    s = GL.standing(bowl, np.hypot(*np.gradient(bowl)) / 474.0)
    assert s[40, 40] < 0.004 and s.min() >= 0.0              # a floor keeps the climate's water; nothing gains


def test_the_covers_colours_are_the_references_over_the_tones_gain():
    r = GL.cover_reflectances()
    back = np.round(GL.appearance.srgb_encode(r['dry'] * GL.REFERENCE_GAIN) * 255)
    assert np.allclose(back, GL.COVERS['dry'], atol=1) and r['swamp'].sum() < r['forest'].sum() < r['green'].sum() < r['dry'].sum()
    assert r['marsh'].sum() < r['oasis'].sum()               # the marsh at a water's edge is the darker of the two wet greens


def test_the_mottle_is_noise_of_the_stated_spread_and_can_be_left_out():
    lat, lon = np.linspace(-10.0, -8.0, 128), np.linspace(20.0, 22.0, 128)
    m = GL.mottle(lat, lon)
    assert abs(float(m.mean()) - 1.0) < 0.01 and 0.8 * GL.LOOK['noise'] < float(m.std()) < 1.25 * GL.LOOK['noise']
    kept = GL.LOOK['noise']
    GL.LOOK['noise'] = 0.0
    try:
        assert float(GL.mottle(lat, lon)) == 1.0
    finally:
        GL.LOOK['noise'] = kept


def test_a_slope_facing_the_noon_sun_has_more_light():
    lat, lon = np.array([-50.0]), np.array([10.0])
    flat = GL.facing(np.zeros((1, 1)), np.zeros((1, 1)), lat, lon, np.eye(3))
    to_sun = GL.facing(np.zeros((1, 1)), np.full((1, 1), 0.3), lat, lon, np.eye(3))      # rising southward: it faces north
    away = GL.facing(np.zeros((1, 1)), np.full((1, 1), -0.3), lat, lon, np.eye(3))
    assert np.isclose(flat[0, 0], 1.0) and to_sun[0, 0] > 1.1 and away[0, 0] < 0.9
    east = GL.facing(np.full((1, 1), 0.3), np.zeros((1, 1)), lat, lon, np.eye(3))        # a slope facing west or east: as level
    assert abs(east[0, 0] - 1.0) < 0.05


def test_the_grounds_dryness_is_held_to_what_the_rain_can_cover():
    """cover_dryness undoes the page's cover curve, so that 40 mm of rain a year over a small demand classes as desert."""
    import appearance
    for a in (0.03, 0.1, 0.22, 0.5, 1.5):
        assert abs(float(GL.cover_dryness(appearance.plant_cover(a))) - a) < 1e-3 * a + 1e-6
    assert np.isinf(GL.cover_dryness(1.0)) and np.isinf(GL.cover_dryness(1.3))      # full cover: no bound
    polar = float(GL.cover_dryness(40.0 / GL.LOOK['rain_full_mm']))                 # 40 mm a year
    assert 0.05 < polar < 0.2                                                       # the arid class, whatever the demand
    assert float(GL.cover_dryness(0.0)) < 0.01
