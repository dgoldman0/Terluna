"""Tests of far_from_look.compose on a small made-up Moon: where the look's average is laid into the far textures."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import far_from_look as FF  # noqa: E402

ROWS, COLS = 180, 360                  # one pixel a degree


def made_up():
    albedo = np.full((ROWS, COLS, 3), 200, np.uint8)
    water = np.zeros((ROWS, COLS), np.uint8)
    average = dict(colour=np.full((ROWS, COLS, 3), 0.05, np.float32), water=np.full((ROWS, COLS), 0.5, np.float32),
                   covered=np.ones((ROWS, COLS), np.float32))
    return albedo, water, average


def test_the_look_is_taken_inside_the_tiles_and_the_texture_kept_beyond_them():
    albedo, water, average = made_up()
    colour, wet, share = FF.compose(albedo, water, average, np.eye(3), chunk=60)
    assert (colour[:9] == 200).all() and (colour[-9:] == 200).all()         # poleward of 81 degrees: as it was
    assert (wet[:9] == 0).all()
    middle = colour[20:160]
    assert middle.max() < 70 and middle.min() > 50                          # linear 0.05 is about 63 in sRGB
    assert abs(int(wet[90, 180]) - 128) <= 1
    assert 0.80 < share < 0.90                                              # the 156 to 160 degrees of latitude cut


def test_ground_the_look_left_empty_keeps_the_texture():
    albedo, water, average = made_up()
    average['covered'][80:100, 20:60] = 0.0          # 10 N to 10 S, 20 E to 60 E on today's Moon
    colour, _, _ = FF.compose(albedo, water, average, np.eye(3), chunk=60)
    # with no turn, today's 20 to 60 E is the texture's columns 200 to 240 (its left edge stands at 180 W)
    assert (colour[82:98, 202:238] == 200).all()
    assert colour[90, 100, 0] < 70


def test_a_turned_frame_moves_the_look_with_the_ground():
    albedo, water, average = made_up()
    average['colour'][:, :180] = 0.5                 # today's eastern half bright
    quarter = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])   # today's 90 E becomes the turned frame's 0
    colour, _, _ = FF.compose(albedo, water, average, quarter, chunk=60)
    # turned longitude L is today's L + 90: turned -90..90 (columns 90 to 270) is today's 0..180, the bright half
    assert colour[90, 180, 0] > 150 and colour[90, 20, 0] < 70 and colour[90, 340, 0] < 70


def test_the_caps_take_the_last_rings_colour_by_longitude_and_keep_their_own_detail():
    albedo, water, average = made_up()
    average['water'][:] = 0.0                          # all land
    average['colour'][10:14, :180] = 0.30              # the northern ring (80 to 76 N): bright over today's eastern half
    albedo[2:4, 20:40, :] = 120                        # a darker patch inside the northern cap
    colour, _, _ = FF.compose(albedo, water, average, np.eye(3), chunk=60, caps=True)
    # with no turn, today's 0 to 180 E is the texture's columns 180 to 360
    east, west = int(colour[5, 270, 0]), int(colour[5, 90, 0])
    assert abs(east - 148) <= 6 and abs(west - 63) <= 6            # linear 0.30 and 0.05 in sRGB
    assert int(colour[3, 30, 0]) < west - 5                         # the darker patch is still darker than its surroundings
    assert (colour[175:, :, 0] > 55).all() and (colour[175:, :, 0] < 72).all()   # the southern cap: its own ring's 0.05
