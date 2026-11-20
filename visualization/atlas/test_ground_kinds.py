"""Checks of the kinds of ground: storm water's load kept on an even slope and laid where it slackens, a hollow's floor
and margin laid from its lowest cell, a hollow's depth, the alluvial ground beside a line, a plain's swales and swells
with its craters left out, the rock on each side of a mare's edge, the kinds sharing out a small piece of ground by
dryness class, and the sand layer."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ground_kinds as K  # noqa: E402
import ground_look as GL  # noqa: E402

INSIDE = (slice(1, -1), slice(1, -1))


def ramp(steep=0.2, gentle=0.02, cell=100.0):
    """Three rows by twenty columns falling eastward: steep over the first ten cells, gentle over the rest."""
    slope = np.where(np.arange(20) < 10, steep, gentle)
    s = np.repeat((slope[::-1].cumsum()[::-1] * cell)[None, :], 3, axis=0) + 0.001 * np.abs(np.arange(3) - 1)[:, None]
    return s.astype(np.float64), np.repeat(slope[None, :], 3, axis=0)


def test_storm_water_keeps_its_load_on_an_even_slope_and_lays_it_where_the_slope_slackens():
    s, slope = ramp()
    wet = np.zeros(s.shape, bool)
    rcv = GL.receivers(s, wet, 100.0, np.full(3, 100.0))
    own = np.ones(s.shape)
    basalt = (np.arange(20) < 5)[None, :] * np.ones((3, 1))
    flow, laid, dark = K.through_flow(s, wet, rcv, own, slope, basalt, settle=0.3)
    assert np.isclose(flow.max(), own.sum()) and np.isclose(laid.sum(), (own * slope).sum())      # none lost, all the load laid
    assert np.allclose(laid[:, :10], 0.0) and laid[1, 10] > 0.0                                   # none on the even slope
    line = laid[:, 10:19].sum(axis=0)
    assert (np.diff(line) < 0.0).all()                                                            # and less with each cell
    assert np.isclose(dark[1, 4], 1.0) and np.isclose(dark.ravel()[np.argmax(flow)], 0.25)        # a quarter of it fell on basalt


def bowl():
    """Nine by nine: a floor of nine cells with a 0.2 m step, and rings 10, 15 and 20 m above its middle."""
    y, x = np.mgrid[0:9, 0:9]
    ring = np.maximum(np.abs(y - 4), np.abs(x - 4))
    return np.where(ring <= 1, 0.2 * ring, 5.0 * ring).astype(np.float64)


def test_a_floor_and_its_margin_are_laid_from_the_lowest_cell_upward():
    s = bowl()
    root = np.full(s.size, 4 * 9 + 4)
    area = np.ones(s.shape)
    first, second = np.zeros(s.shape), np.zeros(s.shape)
    first[4, 4], second[4, 4] = 3.0, 4.0
    floor, margin = K.laid(s, root, first, second, area, floor_m=1.5, margin_m=12.0)
    assert np.isclose(floor.sum(), 3.0) and floor[4, 4] == 1.0 and np.isclose(margin.sum(), 4.0)
    assert np.isclose((floor + margin)[3:6, 3:6].sum(), 7.0) and (floor + margin).max() <= 1.0     # the margin rings the floor
    first[4, 4] = 100.0                                    # however much stands, no more than the floor holds; the margin beyond
    floor, margin = K.laid(s, root, first, second, area, floor_m=1.5, margin_m=12.0)
    assert np.isclose(floor.sum(), 9.0) and np.isclose(margin.sum(), 4.0) and margin[3:6, 3:6].sum() == 0.0
    second[4, 4] = 1000.0                                  # and the margin no higher than the roots reach
    assert np.isclose(K.laid(s, root, first, second, area, 1.5, 12.0)[1].sum(), 16.0)              # the ring 10 m up; none 15 m up
    assert np.allclose(GL.floors(s, root, 2.5 * (first > 0), np.ones(s.shape), area), K.laid(s, root, 0 * first, 2.5 * (first > 0), area, 1.5, 1.5)[1])


def test_a_hollows_depth_is_to_its_lowest_pass():
    s = np.full((5, 11), 50.0)
    s[2, 1:10] = [0.0, 3.0, 9.0, 3.0, 2.0, 6.0, 20.0, 12.0, 10.0]      # three pits along a trough, passes at 9 and 20
    wet = np.zeros(s.shape, bool)
    root = GL.roots(GL.receivers(s, wet, 1.0, np.ones(5)))
    depth = K.hollow_depths(s, root)
    assert depth[2, 1] == 9.0 and depth[2, 5] == 7.0 and depth[2, 9] == 10.0 and depth[2, 2] == depth[2, 1]


def test_alluvial_ground_lies_beside_a_line_that_crosses_the_country_and_stays_low():
    surf = np.repeat((np.abs(np.arange(15) - 7) * 0.5)[None, :], 9, axis=0).astype(np.float32)     # a valley half a metre a cell
    wet = np.zeros(surf.shape, bool)
    slope = np.full(surf.shape, 0.5 / 400.0, np.float32)
    upstream = np.zeros(surf.shape)
    upstream[:, 7] = 16.0 * K.KINDS['line_km2'][1]         # a line draining 2,560 km2: 1.3 km to either side, within 4.2 m of its bed
    floor = K.valley_floors(surf, wet, upstream, 400.0, slope)
    assert floor[4, 7] == 1.0 and floor[4, 8] > 0.9 and 0.0 < floor[4, 10] < floor[4, 9] and floor[4, 13] == 0.0
    high = K.valley_floors(10.0 * surf, wet, upstream, 400.0, 10.0 * slope)                        # five metres a cell: little beside the bed
    assert high[4, 8] < 0.4 and high[4, 9] == 0.0
    upstream[:, 7] = 0.5 * K.KINDS['line_km2'][0]          # a gully of one slope leaves none
    assert K.valley_floors(surf, wet, upstream, 400.0, slope).max() == 0.0
    marked = np.zeros(surf.shape, np.float32)              # unless it is marked: a line of a plain's network, half-marked
    marked[:, 7] = 0.5
    assert np.isclose(K.valley_floors(surf, wet, upstream, 400.0, slope, marked)[4, 7], 0.5 * 2.0 * K.KINDS['valley_m'] * np.sqrt(20.0) / 400.0)


def test_a_plains_swales_and_swells_are_read_and_its_craters_left_out():
    y, x = np.mgrid[0:240, 0:240].astype(np.float32)
    swale = -50.0 * np.exp(-((x - 60.0) / 8.0) ** 2)                               # a trough 50 m deep, 8 km wide, north to south
    swell = 50.0 * np.exp(-((x - 120.0) / 8.0) ** 2)
    r = np.hypot(y - 120.0, x - 190.0)
    crater = -1000.0 * np.clip((16.0 - r) / 4.0, 0.0, 1.0)                         # a crater 15 km across, 1 km deep, flat-floored
    surf = swale + swell + crater
    wet = np.zeros(surf.shape, bool)
    lat = -10.0 - (np.arange(240) + 0.5) / 64.0
    rel = K.plain_ground(surf, wet, lat, 474.0)
    assert rel[120, 60] < -25.0 and rel[120, 120] > 25.0 and abs(rel[120, 10]) < 8.0   # the swale, the swell, the plain between
    assert abs(rel[120, 190]) < 8.0 and rel[120, 176] == 0.0                       # the crater's floor is its own plain; its wall is none
    low, high = K.provinces(rel)
    assert low[120, 60] > 0.3 and high[120, 120] > 0.3 and low[120, 190] == 0.0 and high[120, 60] == 0.0
    part = K.plain_ground(surf[40:, 30:], wet[40:, 30:], lat[40:], 474.0)          # another window reads the same swale
    assert abs(part[80, 30] - rel[120, 60]) < 2.0


def test_the_rock_is_read_on_each_side_of_a_mares_edge():
    import cv2
    mafic = np.zeros((40, 80), np.float32)
    mafic[:, 40:] = 1.0
    feo = cv2.GaussianBlur(np.where(mafic > 0, 18.0, 5.0).astype(np.float32), (0, 0), 6.0)         # a footprint blurs the edge
    sides = K.rock_sides(feo, 0.2 * feo, mafic, cells=4.0)
    assert sides['feo_mare'][20, 40] >= 15.0 and sides['feo_highland'][20, 39] < 9.0               # at the edge itself
    assert abs(sides['feo_mare'][20, 60] - 18.0) < 0.5 and abs(sides['feo_highland'][20, 20] - 5.0) < 0.5


def piece(aridity, basalt=False, sand=None, broad=None, n=24):
    """A bowl 2 km across and 40 m deep on ground that falls gently away from it to open water all round, with the
    fields ground() reads, at one aridity."""
    y, x = np.mgrid[0:n, 0:n].astype(np.float64)
    r = np.hypot(y - n / 2, x - n / 2)
    surf = 40.0 * np.clip(r / 2.0, 0.0, 1.0) ** 2 - 0.5 * np.clip(r - 2.0, 0.0, None)
    wet = np.ones(surf.shape, bool)
    wet[1:-1, 1:-1] = False
    dy = 474.0
    rcv = GL.receivers(surf, wet, dy, np.full(n, dy))
    area = np.full(surf.shape, dy * dy)
    rain, demand = 0.5 * aridity, 0.5
    own = np.full(surf.shape, 0.025 * rain * dy * dy)
    deficit = np.full(surf.shape, max(demand - rain, 0.0))
    kept, held = GL.gathered(surf, wet, rcv, own, deficit * area, 0.0)
    gy, gx = np.gradient(surf, dy)
    a = np.full(surf.shape, aridity, np.float32)
    mare = np.full(surf.shape, basalt)
    rock = dict(feo_highland=np.full(surf.shape, 5.0, np.float32), tio2_highland=np.full(surf.shape, 0.3, np.float32),
                feo_mare=np.full(surf.shape, 18.0, np.float32), tio2_mare=np.full(surf.shape, 3.0, np.float32))
    soil = K.appearance.soil_colour(a, np.where(mare, 18.0, 5.0), np.where(mare, 3.0, 0.3))
    lat = -10.0 - (np.arange(n) + 0.5) / 64.0
    lon = 20.0 + (np.arange(n) + 0.5) / 64.0
    out = K.ground(surf, surf, wet, rcv, GL.roots(rcv), own, held, deficit, 0.05, area, dy, np.hypot(gy, gx).astype(np.float32), a, mare, rock, soil,
                   broad=None if broad is None else np.full(surf.shape, broad, np.float32), sand=sand, lat=lat, lon=lon, dy=dy)
    return out, held


def lum(colour):
    return float((colour[INSIDE] @ K.LUMINANCE).mean())


def test_the_kinds_share_out_the_ground_by_dryness_class():
    arid, held = piece(0.12)
    shares = arid['shares']
    assert all(np.isfinite(v).all() and v.min() >= 0.0 and v.max() <= 1.0 for v in shares.values()) and np.isfinite(arid['open']).all()
    floor_m2 = K.KINDS['stand_ratio'] * held.max() / K.KINDS['pool_evaporation_m']
    assert np.isclose(shares['pan floor, cemented clay'].sum() * 474.0 ** 2, floor_m2, rtol=1e-3)     # the pool's area, bare
    assert arid['bare'].max() > 0.0 and shares['desert soil with lime under a crust'][INSIDE].mean() > 0.9
    assert shares['ground under a biological crust'][INSIDE].mean() > 0.5 and shares['steppe soil, dark humus under a crust'].max() == 0.0
    again, _ = piece(0.12)
    assert all(np.array_equal(shares[k], again['shares'][k]) for k in shares) and np.array_equal(arid['open'], again['open'])   # no chance in it
    steppe, _ = piece(0.4)                                 # a planting pit: the floor planted; the ground the steppe reference's colour
    assert steppe['shares']['pan floor, cemented clay'].max() == 0.0 and steppe['margin'].sum() > 0.0
    assert steppe['shares']['steppe soil, dark humus under a crust'][INSIDE].mean() > 0.9
    assert np.allclose(steppe['open'][5, 5], K.measured('steppe'), atol=1e-4) and lum(steppe['open']) < 0.5 * lum(arid['open'])
    meadow, _ = piece(0.6)                                 # the dry sub-humid land: darker still
    assert meadow['shares']['meadow soil, deep dark humus'][INSIDE].mean() > 0.9 and lum(meadow['open']) < lum(steppe['open'])
    basalt, _ = piece(0.4, basalt=True)                    # cracking clay on basalt: darker, and grass without trees
    assert basalt['shares']['cracking clay on basalt'][INSIDE].mean() > 0.5 and basalt['grassland'][INSIDE].mean() > 0.5
    assert lum(basalt['open']) < 0.75 * lum(steppe['open'])
    desert, _ = piece(0.03)                                # a pavement, no crust to speak of
    assert desert['shares']['stony ground: pavement and steep slopes'][INSIDE].mean() > 0.6
    assert desert['shares']['ground under a biological crust'].max() < 0.01 and lum(desert['open']) < lum(arid['open'])
    wet, _ = piece(1.5)                                    # wet country: no bare floor, pale soil where highland ground is open
    assert wet['bare'].max() == 0.0 and wet['shares']['pale leached soil on highland rock'][INSIDE].mean() > 0.4


def test_a_plains_low_ground_and_swells_take_their_classs_tones():
    level, _ = piece(0.12)
    low, _ = piece(0.12, broad=-200.0)                     # arid: the low ground sandy and cemented, pale; the swell stony, dark
    swell, _ = piece(0.12, broad=200.0)
    assert lum(low['open']) > 1.15 * lum(level['open']) and lum(swell['open']) < 0.85 * lum(level['open'])
    assert low['shares']['low ground of the plains: sandy and cemented where arid, dark and fine beyond'][INSIDE].mean() > 0.9
    level, _ = piece(0.4)
    low, _ = piece(0.4, broad=-200.0)                      # steppe: the low ground dark and fine; the swell the paler, thinner soil
    swell, _ = piece(0.4, broad=200.0)
    assert lum(low['open']) < 0.85 * lum(level['open']) and lum(swell['open']) > 1.15 * lum(level['open'])
    assert swell['shares']['swells of the plains: stony where arid, paler soil beyond'][INSIDE].mean() > 0.9 and swell['swell'].min() == 1.0


def test_a_soda_pan_where_a_hollow_gathers_the_most_ground(monkeypatch):
    plain, _ = piece(0.12)
    assert plain['shares']['soda pan'].max() == 0.0        # a bowl 2 km across gathers less ground than a soda pan needs
    monkeypatch.setitem(K.KINDS, 'soda_km2', 1.0)
    soda, _ = piece(0.12)
    assert soda['shares']['soda pan'].sum() > 0.0 and soda['shares']['pan floor, cemented clay'].max() == 0.0
    steppe, _ = piece(0.4)                                 # and it stays bare where other floors are planted
    assert steppe['shares']['soda pan'].sum() > 0.0
    assert (K.measured('soda') @ K.LUMINANCE) > (K.measured('pan_floor') @ K.LUMINANCE) and K.REFERENCE_GAIN == GL.REFERENCE_GAIN


def made_up_sand(shape, share, active, form, azimuth=30.0, spacing=1500.0):
    """A made-up sand field at a window's fine cells, one value everywhere."""
    full = lambda v, kind=np.float32: np.full(shape, v, kind)   # noqa: E731
    return dict(share=full(share), active=full(active), form=full(form, np.uint8), azimuth_deg=full(azimuth), spacing_m=full(spacing))


def test_without_a_product_or_with_an_empty_one_no_sand_is_drawn():
    rows, cols = np.arange(4000, 4064), np.arange(23000, 23100) % 23040        # across the grid's seam at 0 E
    assert K.sand_window(None, rows, cols) is None and K.read_sand('/nowhere/sand_16ppd.npz', log=lambda s: None) is None
    product = {k: np.zeros(K.SAND_SHAPE, np.uint8 if k == 'form' else np.float32) for k in K.SAND_FIELDS}
    assert K.sand_window(product, rows, cols) is None
    product['sand_share'][1000:1016, :] = 0.5              # a band of sand 16 product cells deep, round the whole Moon
    product['crest_azimuth_deg'][:] = np.nan
    product['crest_azimuth_deg'][1000:1016, :], product['crest_spacing_m'][1000:1016, :], product['form'][1000:1016, :] = 170.0, 1200.0, 3
    window = K.sand_window(product, rows, cols)
    assert window['share'].shape == (64, 100) and np.isclose(window['share'][32, 50], 0.5) and window['share'][0, 0] < 0.5
    assert abs(abs(window['azimuth_deg'][32, 50]) - 10.0) < 0.5 and np.isclose(window['spacing_m'][32, 50], 1200.0) and window['form'][32, 50] == 3
    bare, _ = piece(0.12)
    assert bare['shares']['moving sand'].max() == 0.0 and bare['shares']['sand held under plants or crust'].max() == 0.0


def test_sand_covers_its_share_and_only_the_moving_part_is_bare():
    without, _ = piece(0.12)
    sheet, _ = piece(0.12, sand=made_up_sand((24, 24), 0.4, 0.5, 1, spacing=0.0))
    shares = sheet['shares']
    assert np.allclose(shares['moving sand'][INSIDE], 0.2) and np.allclose(shares['sand held under plants or crust'][INSIDE], 0.2)
    assert shares['moving sand'][0, 0] == 0.0                                      # none over water
    assert np.allclose(sheet['bare'][INSIDE], 0.2 + 0.6 * without['bare'][INSIDE], atol=1e-5)     # the sand lies over the other kinds
    stony = 'stony ground: pavement and steep slopes'
    assert np.allclose(shares[stony], 0.6 * without['shares'][stony], atol=1e-5)


def test_crests_are_a_fixed_function_of_position_and_keep_the_sands_share():
    shape = (200, 240)
    lat, lon = -50.0 - (np.arange(200) + 0.5) / 64.0, 30.0 + (np.arange(240) + 0.5) / 64.0
    whole = K.crests(lat, lon, np.full(shape, 30.0), np.full(shape, 1500.0))
    part = K.crests(lat[50:150], lon[100:200], np.full((100, 100), 30.0), np.full((100, 100), 1500.0))
    assert np.allclose(whole[50:150, 100:200], part, atol=1e-5)                    # another window draws the same crests
    dy = 474.0
    dunes = K.sand_cover(made_up_sand(shape, 0.4, 1.0, 3), lat, lon, dy)
    assert abs(float(dunes.mean()) - 0.4) < 0.05 and dunes.min() < 0.05 and dunes.max() > 0.95      # crests and the ground between
    south, east = np.gradient(dunes.astype(np.float64), dy, dy * np.cos(np.radians(51.5)))         # per metre southward and eastward
    theta = np.radians(30.0)
    along = np.abs(east * np.sin(theta) - south * np.cos(theta)).mean()            # the crests run 30 degrees east of north
    across = np.abs(east * np.cos(theta) + south * np.sin(theta)).mean()
    assert along < 0.3 * across
    star = K.sand_cover(made_up_sand(shape, 0.4, 1.0, 4, spacing=2000.0), lat, lon, dy)             # star dunes: two directions
    assert abs(float(star.mean()) - 0.4) < 0.07
    close = K.sand_cover(made_up_sand(shape, 0.4, 1.0, 3, spacing=600.0), lat, lon, dy)
    assert np.allclose(close, 0.4)                                                 # crests under two pixels apart: an even cover
