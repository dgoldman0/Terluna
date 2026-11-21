"""Checks of the rivers' courses: the far lines' smoothing, parting and width rule, the way over a lake, and the
close views' rivers on a small made ground (their place, their joining and the area their water covers)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import appearance  # noqa: E402
import river_courses as RC  # noqa: E402


def test_a_course_parts_into_channel_and_lake():
    assert RC.pieces(np.array([0, 0, 1, 1, 1, 0], bool)) == [('channel', 0, 1), ('lake', 2, 4), ('channel', 5, 5)]


def test_smoothing_takes_out_the_staircase_and_keeps_the_ends():
    lon = np.arange(0.0, 3.0, 0.0625)
    lat = np.where(np.arange(lon.size) % 2, 0.0625, 0.0)                 # a step up and down at every cell
    u = RC.unit(lat, lon)
    p, km, flow = RC.smooth_line(u, values=np.arange(lon.size, dtype=float))
    la, _ = RC.lat_lon(p)
    assert np.allclose(p[0], u[0]) and np.allclose(p[-1], u[-1])
    middle = slice(len(la) // 4, 3 * len(la) // 4)
    assert np.ptp(la[middle]) < 0.1 * 0.0625
    assert np.all(np.diff(flow) >= 0) and np.isclose(flow[-1], lon.size - 1)
    assert np.isclose(km[-1], np.arccos(np.clip((u[1:] * u[:-1]).sum(axis=1), -1, 1)).sum() * RC.MOON_KM)


def test_the_way_over_a_lake_takes_its_deep_water():
    rows, cols = 40, 200                                                 # the lake spans under half the way round
    water = np.zeros((rows, cols), bool)
    water[10:30, 5:55] = True
    depth = np.full((rows, cols), 2.0)
    depth[25:28, 5:55] = 200.0                                           # a trough along the south
    depth[10:28, 5:8] = depth[10:28, 52:55] = 200.0                      # reached down both ends
    a, b = 10 * cols + 6, 10 * cols + 53
    path = RC.lake_path(a, b, water.ravel(), depth.ravel(), rows, cols, margin=20)
    assert path[0] == a and path[-1] == b and water.ravel()[path].all()
    assert (path // cols).max() >= 25
    apart = RC.lake_path(a, 35 * cols + 30, water.ravel(), depth.ravel(), rows, cols, margin=5)   # no water joins them
    assert list(apart) == [a, 35 * cols + 30]


def test_the_width_drawn_has_a_floor_far_out_and_is_the_waters_close_up():
    half, factor = RC.display_half_width_px(np.array([RC.HEAD_M3S, RC.FLOOR_FLOW]), km_per_px=3.85)
    assert np.allclose(half, RC.FLOOR_HALF_PX) and factor[1] > 4.0
    half, factor = RC.display_half_width_px(np.array([RC.FLOOR_FLOW]), km_per_px=0.25)
    assert factor[0] == 1.0 and np.isclose(half[0], 0.5 * appearance.RIVER_WIDTH * RC.FLOOR_FLOW ** 0.5 / 250.0)


def made_ground():
    """A small globe of 16 x 32 routed cells; the window is its western 16 x 16. A river of 1,000 m3/s runs south down
    routed column 8 into a lake over the four southern rows, in a valley along fine column 33; a tributary of 100 m3/s
    comes in from the west along routed row 5."""
    rows16, cols16 = 16, 32
    q = np.zeros((rows16, cols16), np.float32)
    water = np.zeros((rows16, cols16), bool)
    cells = np.arange(rows16 * cols16).reshape(rows16, cols16)
    receiver = cells.copy()
    water[12:] = True
    q[:12, 8] = 1000.0
    receiver[:12, 8] = cells[1:13, 8]
    q[5, 5:8] = 100.0
    receiver[5, 5:8] = cells[5, 6:9]
    routing = dict(rows=rows16, cols=cols16, q=q.ravel(), water=water.ravel(), receiver=receiver.ravel().astype(np.int64))
    y, x = np.mgrid[0:64, 0:64].astype(np.float32)
    h = 100.0 - 0.5 * y + 30.0 * np.abs(x - 33.0)
    wet = y >= 48
    return h, wet, routing


def test_fine_rivers_lie_in_the_valley_at_their_width():
    h, wet, routing = made_ground()
    out = RC.fine_rivers(h, wet, np.zeros(64), 0, 0, routing)
    main, side = out['courses'][0], out['courses'][1]
    assert np.allclose(main['x'][3:-3], 33.0, atol=0.3) and main['y'][-1] >= 47
    dy = RC.MOON_KM * 1000.0 * np.radians(1.0 / 64)
    width = appearance.RIVER_WIDTH * 1000.0 ** 0.5 / dy                  # in fine cells: 0.77
    assert np.allclose(out['share'][10:18, 32:35].sum(axis=1), width, atol=0.05)        # under a pixel wide: the share of the pixel
    assert out['share'][10:18, 30].max() == 0.0
    assert out['flow'][15, 33] == 1000.0 and out['sediment'].max() == 0.0
    # the tributary ends on the river it joins, and its own water is a seventh of a pixel wide
    assert np.isclose(side['x'][-1], 33.0, atol=0.3) and 20 <= side['y'][-1] <= 23
    column = int(round(side['x'][4]))
    assert np.isclose(out['share'][18:26, column].sum(), appearance.RIVER_WIDTH * 10.0 / dy, atol=0.05)
    # the water's area is each course's length times its width
    area = sum(float((appearance.RIVER_WIDTH * np.sqrt(0.5 * (k['flow'][1:] + k['flow'][:-1])) / dy
                      * np.hypot(np.diff(k['y']), np.diff(k['x']))).sum()) for k in out['courses'])
    assert abs(float(out['share'].sum()) - area) < 0.06 * area


def test_no_river_no_water():
    h, wet, routing = made_ground()
    routing['q'] = np.zeros_like(routing['q'])
    out = RC.fine_rivers(h, wet, np.zeros(64), 0, 0, routing)
    assert out['share'].max() == 0.0 and out['courses'] == []


def test_a_river_between_two_lakes_runs_from_one_to_the_other():
    """A lake that takes a stream in and sends the river out is the end of the stream and the head of the river: the
    river's course must still run its whole way down to the lower lake."""
    h, wet, routing = made_ground()
    rows16, cols16 = routing['rows'], routing['cols']
    q, water, receiver = (routing[k].reshape(rows16, cols16).copy() for k in ('q', 'water', 'receiver'))
    cells = np.arange(rows16 * cols16).reshape(rows16, cols16)
    water[:2] = True                                      # an upper lake over the two northern rows
    q[:2, 8] = 0.0
    receiver[1, 8] = cells[2, 8]                          # it passes its water to the river
    q[1, 8] = 1000.0
    q[2:5, 3] = 50.0                                      # and a stream runs north into it
    receiver[4, 3], receiver[3, 3], receiver[2, 3] = cells[3, 3], cells[2, 3], cells[1, 3]
    routing = dict(routing, q=q.ravel(), water=water.ravel(), receiver=receiver.ravel())
    wet = wet | (np.mgrid[0:64, 0:64][0] < 8)
    out = RC.fine_rivers(h, wet, np.zeros(64), 0, 0, routing)
    main = out['courses'][0]
    assert out['unplaced'] == 0 and main['y'][0] <= 8 and main['y'][-1] >= 47
    drawn = out['flow'].reshape(16, 4, 16, 4).max(axis=(1, 3)) >= 1000.0
    assert drawn[2:12, 8].all()                           # every routed cell of the river carries its course
    stream = out['flow'].reshape(16, 4, 16, 4).max(axis=(1, 3)) == 50.0
    assert stream[2:5, 3].all()
