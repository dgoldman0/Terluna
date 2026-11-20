"""The land's close-up look made from the ground itself, in windows of 2 x 2 tiles: a second trial beside fine_look.py.

    ulimit -v 10000000; OPENBLAS_NUM_THREADS=1 python3 visualization/atlas/ground_look.py \\
        --surface <globe_surface_16ppd.npz> --moisture research/studies/atlas_round/results/moisture_map_final.json \\
        --window Paqarinwayqu:-17,174 --window "Erg Ahmar:-50,12" --out <folder> [--sand <sand_16ppd.npz>]
                                                                                   # 15 to 40 seconds a window

fine_look.py multiplied the moisture map's dryness by terms for hollows, drainage lines and noise, with nothing to
bound how green a hollow became. Here a 474 m pixel's colour is the area-weighted mix of what stands in it: open
water, ground kept wet, rain-fed plants and the open ground between them, which is one of the soils study's kinds of
ground. The moisture map gives each place its yearly rain and demand (the regional base, from the climate's 170 km
cells); everything finer comes from the fine heights, the water and the mare map's outlines. Every rule below is a
display guess unless it names a product or a study, and LOOK and ground_kinds.KINDS mark each setting that is one.

- Gathered water. Of the rain on each cell the share that storms gather (the moisture map's storm_water_surfaces, the
  biological crust's row, the pair's mean, by dryness class) runs down the steepest way over the fine ground. Each
  cell it crosses takes up LOOK['loss'] of what arrives, as far as its ground's deficit (demand less its own rain)
  allows; what reaches a cell with no lower neighbour stays there. Water kept in a cell keeps wet the area whose
  deficit it meets. A hollow's water keeps a planted margin round the floor where its storm water stands and dries
  (ground_kinds.py), over the area that water can keep wet and no further, and the yearly rain bounds it.
- Corridors and shores. A river, or a routed stream of 1 m3/s or more that drains 100 km2, keeps its valley floor
  wet: the ground within a few metres of its surface (more for a larger river), as far to either side as
  LOOK['corridor_share'] of its yearly flow, per kilometre of course, meets the ground's deficit, so a small stream in
  dry country has a narrow belt and any stream in wet country a wide one. Open fresh water keeps a belt over its
  groundwater, as the soils study describes the shores (2.6): the ground within LOOK['shore_m'] of its surface, to
  LOOK['shore_reach_m'] inland, with marsh on the ground within LOOK['marsh_m'] of the water's level. A pixel beside
  water is taken to rise from the water's level, so a steep shore keeps a narrow margin.
- The rest lives on its own rain: the page's cover curve (appearance.plant_cover) on rain over demand, never more
  cover than LOOK['rain_full_mm'] of yearly rain supports, with what the ground does to it: the demand raised or
  lowered by the way a slope faces the noon Sun; and a little less water for plants where the ground stands well out
  of its surroundings (knolls, rims and ridges at 0.7, 1.9 and 5.7 km), stands in the upper half of the relief round
  it, is steep (standing, and the page's rule for slopes) or lies on a swell of its plain (ground_kinds.provinces). Ground
  that stands out sheds its storm water, and the cut is held to the share of its rain that a bare loam sheds
  (LOOK['stand_most']); in mosaics of Earth's dry country rough ground is darker than the plains beside it, so a
  rim is drawn close to its surroundings' tone.
  Ground that stands out never gains water, so no cell has more than its rain gives.
- The open ground between the plants is placed kind by kind by ground_kinds.py: each dryness class's own soil, the
  deeper soil of a plain's low ground and the thinner of its swells, stony ground, the beds, fans and valley floors
  that storm water leaves, the floors and soda pans of hollows, cracking clay on the mare remnants (their outlines
  read on the fine grid), silt beside laden rivers, and loose sand where the sand product (--sand) lays it. Every
  patch is placed by the heights, the water, the rock maps, the routing or the climate's cells. The kinds and the
  soil's colour follow the climate's dryness, held no wetter than the dryness at which the cover curve gives the
  cover the yearly rain can support (cover_dryness): towards the poles the demand falls to 90 to 140 mm a year, and
  30 to 60 mm of rain over it would otherwise class a desert's ground as steppe.
- Each kind of plant cover has a colour measured in the reference mosaics (COVERS: wet forest, forest, green herbs,
  dry herbs, the green of a corridor in dry country and the marsh at a water's edge), as the reflectance that the
  close view's tone shows as that colour; Earth's colours guide here and do not bind.
- Rivers of 30 m3/s or more lie in the valleys of the fine ground at their computed widths (river_courses.fine_rivers).
  Water takes the designed waters' colours by the page's rule for depth (the productive coast's colour from 20 m,
  fully at 3 m) and by the plumes the rivers bring; a cell at the water's edge is taken to shelve from the shore
  down to twice its mean depth, so that its share of shallows shows, and the bed shows through the last metre.
- The finished picture is blurred just under half a pixel, as an instrument's optics blur what they record.
- Noise, and the only noise: the plants' brightness is mottled by LOOK['noise'] (standard deviation) at 0.7 and
  2.2 km, as a canopy's is. With `noise=0` every difference between two pixels comes from the heights, the water, the
  rock maps or the climate's cells; the record gives the picture's measures both ways.

It says nothing of which plants grow.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import appearance  # noqa: E402
import fine_tiles as FT  # noqa: E402
import ground_kinds as K  # noqa: E402
import river_courses as RC  # noqa: E402
from research.studies.atlas_round import fine_surface as FS  # noqa: E402

PPD, TILE = FT.PPD, FT.TILE
MARGIN = 128                        # fine cells kept round a window (61 km), so that water gathers across its edges
SALT = 3                            # the water code of the salt sea, whose shores keep no fresh belt
CLASS_MIDDLES = (('arid', 0.1), ('semi-arid', 0.316), ('dry sub-humid', 0.57))   # the classes' middles, in aridity
# The covers' colours, sRGB, measured in the reference mosaics' colour groups (NASA's Blue Marble: Next Generation at
# about 550 m a pixel; the agents' reference_forms.json of 19 November 2026), and the gain of the close view's tone
# (close_light.CLOSE['gain']): a cover's reflectance here is its reference colour over that gain, so that level ground
# of it shows as the reference does. A display guess: the mosaics are stretched pictures, not reflectances.
REFERENCE_GAIN = 1.6
COVERS = dict(swamp=(24, 39, 11),    # closed wet forest: the Amazon scene's largest group
              forest=(38, 51, 17),   # forest on drained ground: the Uganda scene's forest group
              green=(56, 62, 28),    # green herbs and savanna in leaf: the Uganda scene's green cover (Turkana's 60, 63, 30)
              dry=(110, 88, 57),     # dry herbs: the savanna and steppe groups of the Okavango, Balkhash and Turkana scenes
              oasis=(45, 58, 28),    # a green corridor in dry country: between the Nile's (60, 66, 37) and the Ili delta's (30, 50, 23)
              marsh=(30, 50, 23))    # reed and sedge marsh at a fresh water's edge: the Ili delta's, in the Balkhash scene
LOOK = dict(
    gathered_surface='biological crust',   # the soil surface whose gathered share is read (the drylands' kept surface)
    loss=0.15,                  # guess: the share of the water running onto a cell that its ground takes up (per 474 m)
    floor_m=1.5,                # guess: a hollow's floor is the ground within this of its lowest cell
    deficit_least_m=0.05,       # the least deficit gathered or river water is spread against (m a year)
    rain_full_mm=250.0,         # guess: the yearly rain under which no ground is fully covered, whatever the demand
    direct=0.5,                 # guess: the direct beam's share of a day's light, for the way a slope faces
    shore_m=20.0,               # the soils study, 2.6: the belt over fresh groundwater reaches 15 to 25 m above the water
    shore_reach_m=1500.0,       # on steep shores and 0.6 to 2.5 km inland on flat ones; the middle of each
    marsh_m=1.0,                # guess: ground within this of open fresh water's level is marsh
    corridor_share=0.01,        # guess: the share of a river's yearly flow that each km of its course gives its valley floor
    corridor_reach_m=6000.0,    # guess: the farthest a corridor reaches to either side
    corridor_base_m=2.0,        # guess: a stream's valley floor is the ground within this of its surface,
    corridor_rise_m=1.0,        # guess: and this much higher for each tenfold of its flow
    stream_m3s=1.0,             # the routed streams that keep a corridor, from this flow up (water is drawn from river_m3s),
    stream_km2=100.0,           # guess: where they drain at least this much ground: in very wet country every 1.9 km cell
    #                             carries 1 m3/s, and corridors on all of them comb the slopes at the routing's spacing
    stand_out_px=(1.5, 4.0, 12.0),     # guess: ground standing above the mean of its surroundings at these scales (cells)
    stand_out_m=(40.0, 100.0, 300.0),  # by more than half of these heights has its plants' water cut: by a factor
    stand_out_ln=(0.04, 0.07, 0.07),   # exp(-stand_out_ln) at one and a half times the height, in proportion to twice it
    height_px=6, height_ln=0.1,        # guess: and by up to exp(-height_ln) in the upper half of the relief within these cells,
    relief_m=(40.0, 400.0),            # where that relief is real: none of it under the first height, all over the second
    sharp=(0.08, 0.30), worn=0.3,      # guess: the whole cut scaled by the steepest slope within two cells: by `worn` under
    #                                    the first slope (5 degrees), in full over the second (17 degrees)
    stand_most=0.22,                   # all of it together cuts the water by no more than the share of its rain that a bare
    #                                    loam sheds, 19.8% (the soils study, 3.1), and the three settings above are scaled to
    #                                    that; in eight mosaics of Earth's dry country rough ground is 16 to 35% darker than
    #                                    smooth ground and nowhere paler
    facing_px=1.2,              # the heights smoothed over this before the slopes' steepness and facing are read for the
    #                             cover and the rock: their track noise is finer
    wet_between=(0.3, 1.0),     # guess: wet ground's cover turns from the dry country's corridor green to wet forest
    #                             between these aridities
    psf_px=0.45,                # display: the picture blurred as an instrument's optics blur it (a Gaussian, in pixels)
    noise=0.05, noise_km=(0.7, 2.2),   # NOISE: the plants' brightness mottled by this standard deviation at these scales
    bed_m=1.0, bed_share=0.4, bed_dark=0.4,   # guess: the bed shows through the last metre of water, as dark wet soil
    river_m3s=RC.FINE_FLOW)


def cover_reflectances():
    """COVERS as reflectances (linear sRGB): each reference colour decoded and divided by REFERENCE_GAIN."""
    out = {}
    for k, c in COVERS.items():
        v = np.asarray(c, np.float32) / 255.0
        out[k] = np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4).astype(np.float32) / REFERENCE_GAIN
    return out


def carried(record, name, rows=720):
    """The logarithm of a column of the moisture map's cells on a grid of 4 px/deg in the settled frame (180 W at
    the left edge), as appearance.moisture_field carries the aridity: each land cell at its centre, every other cell
    its nearest land cell's, regridded and smoothed as the climate's fields are."""
    from scipy.spatial import cKDTree
    lat = np.asarray(record['grid']['lat'], float)
    lon = np.asarray(record['grid']['lon'], float)
    cells = record['cells']
    f = np.full((lat.size, lon.size), np.nan)
    f[np.asarray(cells['row']), np.asarray(cells['col'])] = np.log(np.asarray(cells[name], float))
    la, lo = np.meshgrid(np.radians(lat), np.radians(lon), indexing='ij')
    points = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)
    land = np.isfinite(f)
    _, nearest = cKDTree(points[land]).query(points[~land])
    f[~land] = f[land][nearest]
    return appearance.grid_field(f, lat, lon, (rows, 2 * rows)).astype(np.float32)


def sample4(field, lat_t, lon_t):
    """A 4 px/deg field in the settled frame read bilinear at settled latitudes and longitudes."""
    import cv2
    x = ((lon_t + 180.0) * 4.0 - 0.5).astype(np.float32)
    y = ((90.0 - lat_t) * 4.0 - 0.5).astype(np.float32)
    return cv2.remap(field, x, y, cv2.INTER_LINEAR, borderMode=cv2.BORDER_WRAP)


def gathered_share(aridity, table, surface=None):
    """The share of the rain that storms gather into the low ground, from the moisture map's storm_water_surfaces
    (per cent of the rain, by soil surface and dryness class; the mean of each pair): interpolated in the logarithm
    of the aridity between the classes' middles and held beyond them."""
    row = table[surface or LOOK['gathered_surface']]
    at = np.log([m for _, m in CLASS_MIDDLES])
    values = [float(np.mean(row[k])) for k, _ in CLASS_MIDDLES]
    return (np.interp(np.log(np.asarray(aridity, float)), at, values) / 100.0).astype(np.float32)


def receivers(s, wet, dy, dx):
    """The cell each cell's water runs to: the steepest way down among its eight neighbours (dy the rows' spacing in
    metres, dx each row's column spacing), or -1 where none lies lower. Water takes what reaches it and passes none on."""
    rows, cols = s.shape
    pad = np.pad(s, 1, mode='edge')
    best = np.zeros(s.shape, np.float64)
    rcv = np.full(s.shape, -1, np.int64)
    idx = np.arange(rows * cols).reshape(rows, cols)
    rr, cc = np.arange(rows)[:, None], np.arange(cols)[None, :]
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            drop = (s - pad[1 + dr:1 + dr + rows, 1 + dc:1 + dc + cols]) / np.hypot(dr * dy, dc * np.asarray(dx))[:, None]
            better = (drop > best) & (rr + dr >= 0) & (rr + dr < rows) & (cc + dc >= 0) & (cc + dc < cols)
            rcv = np.where(better, idx + dr * cols + dc, rcv)
            best = np.where(better, drop, best)
    rcv[wet] = -1
    return rcv


def gathered(s, wet, rcv, own, room, loss):
    """Water run down the receivers from the highest cell to the lowest. own is the water each cell sheds (m3 a year)
    and room what its ground could still take up. A cell on the way takes `loss` of what arrives from above, as far as
    its room allows, and passes the rest and its own on; a cell with no receiver holds all of it. Returns what each
    cell on the way keeps and what each end holds. Water that reaches open water is the water's."""
    n = s.size
    order = np.argsort(s, axis=None)[::-1]
    order = order[~wet.ravel()[order]]
    r, o, c = rcv.ravel().tolist(), own.ravel().tolist(), room.ravel().tolist()
    arriving, kept, held = [0.0] * n, [0.0] * n, [0.0] * n
    for i in order.tolist():
        a = arriving[i]
        j = r[i]
        if j < 0:
            held[i] = a + o[i]
            continue
        u = loss * a
        if u > c[i]:
            u = c[i]
        kept[i] = u
        arriving[j] += a - u + o[i]
    return np.asarray(kept).reshape(s.shape), np.asarray(held).reshape(s.shape)


def roots(rcv):
    """The cell each cell's water ends in, following the receivers."""
    flat = rcv.ravel()
    root = np.where(flat >= 0, flat, np.arange(flat.size))
    for _ in range(64):
        further = root[root]
        if np.array_equal(further, root):
            break
        root = further
    return root


def floors(s, root, held, deficit, area):
    """The share of each cell a hollow's held water keeps wet: the water over the deficit is an area, laid over the
    hollow's cells from its lowest upward, within LOOK['floor_m'] of the lowest."""
    budget = held / np.maximum(deficit, LOOK['deficit_least_m'])
    return K.laid(s, root, np.zeros(s.shape), budget, area, LOOK['floor_m'], LOOK['floor_m'])[1]
    order = cand[np.lexsort((z[cand], root[cand]))]
    a = area.ravel()[order].astype(np.float64)
    upto = np.cumsum(a)
    group = root[order]
    first = np.r_[True, group[1:] != group[:-1]]
    before = upto - a - np.maximum.accumulate(np.where(first, upto - a, 0.0))
    out[order] = np.clip((budget[order] - before) / a, 0.0, 1.0)
    return out.reshape(s.shape)


def beside_water(h, wet, level, fresh, river_flow, lat, dy, slope, deficit, rise=None):
    """The share of each cell kept wet by the fresh water, river or stream nearest it: the part of the cell's ground
    within reach of that water's surface, in height and in distance. In height: LOOK['shore_m'] above open water (or
    `rise`), and above a river LOOK['corridor_base_m'] and LOOK['corridor_rise_m'] more for each tenfold of its flow;
    a cell's ground spans its slope times its width
    about its mean height, and a cell touching the water is taken to rise from the water's level to twice its mean
    height above it. In distance: LOOK['shore_reach_m'] from open water; from a river, as far to either side as
    LOOK['corridor_share'] of its yearly flow per kilometre of course meets the cell's deficit (m a year), at most
    LOOK['corridor_reach_m'], and where that is under half a cell, that share of the river's own cell."""
    import cv2
    source = (wet & fresh) | (river_flow > 0.0)
    if not source.any():
        return np.zeros(h.shape, np.float32)
    dist, label = cv2.distanceTransformWithLabels((~source).astype(np.uint8), cv2.DIST_L2, 5, labelType=cv2.DIST_LABEL_PIXEL)
    index = np.zeros(int(label.max()) + 1, np.int64)
    index[label[source]] = np.flatnonzero(source)
    near = index[label]
    surface = np.where(wet, level, h).ravel()[near]
    flow = np.where(wet, 0.0, river_flow).ravel()[near]
    z = h - surface
    cell = dy * np.sqrt(np.cos(np.radians(lat)))[:, None]              # a cell's width, the mean of its two sides
    corridor = 0.5 * LOOK['corridor_share'] * flow * 3.156e7 / (np.maximum(deficit, LOOK['deficit_least_m']) * 1000.0)
    reach = np.where(flow > 0.0, np.minimum(corridor, LOOK['corridor_reach_m']), LOOK['shore_reach_m'])
    rise = np.where(flow > 0.0, LOOK['corridor_base_m'] + LOOK['corridor_rise_m'] * np.log10(np.maximum(flow, 1.0)),
                    LOOK['shore_m'] if rise is None else rise)
    span = np.maximum(slope * cell, 0.5)
    touching = dist <= 1.5
    low = np.where(touching, 0.0, z - 0.5 * span)
    high = np.where(touching, np.maximum(2.0 * z, 0.5), z + 0.5 * span)
    part = np.clip((rise - low) / (high - low), 0.0, 1.0)
    within = np.where(dist < 0.5, np.clip(2.0 * reach / cell, 0.0, 1.0), appearance.smoothstep(reach, 0.3 * reach, dist * cell))
    return np.where(wet, 0.0, part * within).astype(np.float32)


def facing(gx, gy, lat, lon, rotation, direct=None):
    """The day's light on a slope over that on level ground, by the noon Sun's angle on it: the Sun stands over the
    settled frame's equator, on the place's own meridian. gx and gy are appearance.gradients' (east, and south)."""
    direct = LOOK['direct'] if direct is None else direct
    z = np.asarray(rotation, float)[2]                     # the settled pole, in today's frame
    la, lo = np.radians(lat)[:, None], np.radians(lon)[None, :]
    east = np.stack(np.broadcast_arrays(-np.sin(lo), np.cos(lo), np.zeros_like(lo * la)), axis=-1)
    north = np.stack(np.broadcast_arrays(-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la) + 0.0 * lo), axis=-1)
    up = np.stack(np.broadcast_arrays(np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la) + 0.0 * lo), axis=-1)
    pz = up @ z
    level = np.sqrt(np.clip(1.0 - pz * pz, 1e-4, None))    # the noon Sun's height on level ground, as a cosine
    se, sn = -pz * (east @ z) / level, -pz * (north @ z) / level
    on_slope = (-gx * se + gy * sn + level) / np.sqrt(1.0 + gx * gx + gy * gy)
    return np.clip((1.0 - direct) + direct * np.clip(on_slope, 0.0, None) / level, 0.4, 1.8).astype(np.float32)


def standing(surf, slope):
    """How far the ground at each cell stands out of or above its surroundings, as the logarithm of the factor that
    cuts its plants' water: 0 on floors, slopes and level ground, and on any ground that stands less than half a
    scale's height above the mean of its surroundings; rising to LOOK['stand_out_ln'] at one and a half times that
    height and on in proportion to twice it, at each of three scales; up to LOOK['height_ln'] more in the upper half
    of the relief within LOOK['height_px'] cells, where that relief is real (none under LOOK['relief_m'][0], all of
    it over LOOK['relief_m'][1]). The whole is then scaled by how sharp the ground is: by LOOK['worn'] where the
    steepest slope within two cells is under LOOK['sharp'][0], in full where it is over LOOK['sharp'][1], so that a
    sharp, high rim carries a little less cover than a worn one. LOOK['stand_most'] at most. Heights in metres and
    fixed scales, so that neighbouring windows agree."""
    import cv2
    surf = np.asarray(surf, np.float32)
    out = np.zeros(surf.shape, np.float32)
    for px, metres, ln in zip(LOOK['stand_out_px'], LOOK['stand_out_m'], LOOK['stand_out_ln']):
        out += ln * np.clip((surf - cv2.GaussianBlur(surf, (0, 0), px)) / metres - 0.5, 0.0, 1.5)
    round_it = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * LOOK['height_px'] + 1,) * 2)
    low, high = cv2.erode(surf, round_it), cv2.dilate(surf, round_it)
    out += (LOOK['height_ln'] * np.clip(2.0 * (surf - low) / (high - low + 5.0) - 1.0, 0.0, 1.0)
            * appearance.smoothstep(LOOK['relief_m'][0], LOOK['relief_m'][1], high - low))
    steepest = cv2.dilate(np.asarray(slope, np.float32), cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
    out *= LOOK['worn'] + (1.0 - LOOK['worn']) * appearance.smoothstep(LOOK['sharp'][0], LOOK['sharp'][1], steepest)
    return np.minimum(out, LOOK['stand_most'])


def mottle(lat, lon):
    """NOISE: a factor on the plants' brightness, mean 1 and standard deviation LOOK['noise'], from value noise on the
    sphere at LOOK['noise_km'] (the page's own noise function); nothing in the ground places it."""
    if not LOOK['noise']:
        return np.float32(1.0)
    la, lo = np.radians(lat)[:, None], np.radians(lon)[None, :]
    p = np.stack(np.broadcast_arrays(np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la) + 0.0 * lo), axis=-1)
    n = sum(appearance.value_noise(p * (appearance.MOON_KM / km), 97 + 13 * i) for i, km in enumerate(LOOK['noise_km']))
    n = n / len(LOOK['noise_km']) - 0.5
    return (1.0 + LOOK['noise'] * n / 0.142).astype(np.float32)       # 0.142: the spread of two such fields' mean (measured)


def cover_dryness(cover):
    """The aridity (rain over demand) at which the page's cover curve (appearance.plant_cover) gives a cover: without
    bound from the curve's fullest cover on, and next to nothing at its least."""
    c = appearance.COVER
    x = (c['most'] - c['least']) / np.clip(np.asarray(cover, np.float64) - c['least'], 1e-6, None) - 1.0
    return np.where(x > 1e-6, c['half'] / np.maximum(x, 1e-6) ** (1.0 / c['steep']), np.inf)


def look_window(tr, tc, g, size=2 * TILE, margin=MARGIN, origin=None):
    """The look of the window whose top left tile is (tr, tc), or of the piece of ground whose first fine row and
    column are `origin` (multiples of 4) and whose rows and columns are `size`: its colour (uint8 sRGB), the share of
    each pixel under water, and the fields the colour was mixed from."""
    import cv2
    from climate.gcm.turned import to_turned
    r0, c0 = origin if origin is not None else ((90 - FT.LAT_TOP + FT.TILE_DEG * tr) * PPD, FT.TILE_DEG * tc * PPD)
    size_r, size_c = (size, size) if np.isscalar(size) else size
    rows = np.arange(r0 - margin, r0 + size_r + margin)
    cols = np.arange(c0 - margin, c0 + size_c + margin) % FS.COLS
    h = g['height'][rows[0]:rows[-1] + 1][:, cols].astype(np.float32) * 0.5
    level = g['level'][rows[0]:rows[-1] + 1][:, cols].astype(np.float32) * 0.5
    code = g['water'][rows[0]:rows[-1] + 1][:, cols]
    wet = code > 0
    lat, lon = 90.0 - (rows + 0.5) / PPD, (cols + 0.5) / PPD
    dy = appearance.MOON_KM * 1000.0 * np.radians(1.0 / PPD)
    dx = dy * np.cos(np.radians(lat))
    area = np.repeat((dy * dx)[:, None], len(cols), axis=1)
    surf = np.where(wet, level, h)

    # the climate's cells: the yearly rain and demand (m), and the rock
    lat_t, lon_t = to_turned(np.repeat(lat[:, None], len(cols), axis=1), np.repeat(lon[None, :], len(rows), axis=0), g['rotation'])
    rain = np.exp(sample4(g['ln_rain'], lat_t, lon_t)) / 1000.0
    demand = np.exp(sample4(g['ln_demand'], lat_t, lon_t)) / 1000.0
    aridity = rain / demand
    share = gathered_share(aridity, g['storm'])
    gy, gx = appearance.gradients(surf, dy, lat)
    steepness = np.hypot(gx, gy)                           # as the heights stand, for the ground a pixel spans at a shore
    gy, gx = appearance.gradients(cv2.GaussianBlur(surf, (0, 0), LOOK['facing_px']), dy, lat)
    slope = np.hypot(gx, gy)                               # with the heights' track noise smoothed out, for the cover and the rock
    light = facing(gx, gy, lat, lon, g['rotation'])
    own_rain = rain * (1.0 - share)
    deficit = np.maximum(demand * light - own_rain, 0.0)

    # where water gathers
    s = surf.astype(np.float64) + 1e-3 * cv2.GaussianBlur(surf, (0, 0), 4.0)      # level ground drains the way the land round it leans
    rcv = receivers(s, wet, dy, dx)
    kept, held = gathered(s, wet, rcv, share * rain * area, deficit * area, LOOK['loss'])
    on_lines = np.where(deficit > 0.0, np.clip(kept / np.maximum(deficit * area, 1e-9), 0.0, 1.0), 0.0)
    root = roots(rcv)
    rivers = RC.fine_rivers(h, wet, lat, int(rows[0]), int(cols[0]), g['streams'], g.get('sediment'), LOOK['stream_m3s'])
    rivers['share'] = np.where(rivers['flow'] >= LOOK['river_m3s'], rivers['share'], 0.0).astype(np.float32)   # water from 30 m3/s up
    fresh = code != SALT
    shore = beside_water(h, wet, level, fresh, rivers['flow'], lat, dy, steepness, deficit)
    marsh = beside_water(h, wet, level, fresh, np.zeros_like(rivers['flow']), lat, dy, steepness, deficit, rise=LOOK['marsh_m'])
    shore = np.maximum(shore, marsh)

    # the plants' own rain, and the kinds of ground between them (ground_kinds.py)
    broad = K.plain_ground(surf, wet, lat, dy)             # each cell's place on its plain (m below or above the plain round it)
    stands = np.minimum(standing(surf, steepness) + K.KINDS['swell_ln'] * K.provinces(broad)[1], LOOK['stand_most'])
    fed = own_rain / (demand * light) * np.exp(-stands - appearance.SLOPE_LN * appearance.smoothstep(0.2, 0.6, slope))
    cover = np.minimum(appearance.plant_cover(fed), rain * 1000.0 / LOOK['rain_full_mm'])
    # The dryness the ground's own kinds and soil colours follow: the climate's rain over demand, and no wetter than the
    # dryness at which the cover curve gives the cover the yearly rain can support. Where the demand is small, as in the
    # polar dimness, 30 to 60 mm of rain a year over 90 to 140 mm of demand would otherwise class a desert's ground as
    # steppe and give it a steppe's dark humus soil. A display guess, the same bound as the cover's.
    ground_aridity = np.minimum(aridity, cover_dryness(rain * 1000.0 / LOOK['rain_full_mm'])).astype(np.float32)
    if 'mare' in g:                                        # the mare map on the fine grid, the rock on each side of its edge
        mare = g['mare'][rows[0]:rows[-1] + 1][:, cols]
        rock = {k: FT.sample16(g[k], lat_t, lon_t) for k in ('feo_highland', 'tio2_highland', 'feo_mare', 'tio2_mare')}
    else:
        mare = FT.sample16(g['mafic'], lat_t, lon_t) > 0.5
        feo, tio2 = FT.sample16(g['feo'], lat_t, lon_t), FT.sample16(g['tio2'], lat_t, lon_t)
        rock = dict(feo_highland=feo, tio2_highland=tio2, feo_mare=feo, tio2_mare=tio2)
    soil = appearance.soil_colour(ground_aridity, np.where(mare, rock['feo_mare'], rock['feo_highland']),
                                  np.where(mare, rock['tio2_mare'], rock['tio2_highland']))
    plume = np.clip(FT.sample16(g['plume'], lat_t, lon_t), 0.0, 1.0)
    reach = K.KINDS['mouth_px']                            # the plume of the water within a few pixels, for the silt at a mouth
    beside = cv2.GaussianBlur(plume * wet, (0, 0), reach) / np.maximum(cv2.GaussianBlur(wet.astype(np.float32), (0, 0), reach), 0.05)
    cell = dy * np.sqrt(np.cos(np.radians(lat)))[:, None]
    kinds = K.ground(surf, s, wet, rcv, root, share * rain * area, held, deficit, LOOK['deficit_least_m'], area, cell, slope,
                     ground_aridity, mare, rock, soil, broad=broad, laden=rivers['share'] * rivers['sediment'], mouth=marsh * beside,
                     sand=K.sand_window(g.get('sand'), rows, cols), lat=lat, lon=lon, dy=dy)
    on_floors = kinds['margin']
    wetted = 1.0 - (1.0 - np.clip(on_lines + on_floors, 0.0, 1.0)) * (1.0 - shore)

    # the mix: the kinds of ground that carry no plants; on the rest wet ground, rain-fed plants and the open ground
    # between them; then the rivers and the open water
    ref = cover_reflectances()
    trees = appearance.share_of(fed, appearance.TREES)[..., None] * (1.0 - kinds['grassland'][..., None])
    green = appearance.share_of(fed, appearance.GREEN)[..., None]
    grain = np.asarray(mottle(lat, lon))[..., None]
    canopy = (ref['forest'] * trees + (ref['dry'] * (1.0 - green) + ref['green'] * green) * (1.0 - trees)) * grain
    turn = appearance.smoothstep(np.log(LOOK['wet_between'][0]), np.log(LOOK['wet_between'][1]), np.log(aridity))[..., None]
    reeds = (np.minimum(marsh, wetted) / np.maximum(wetted, 1e-6))[..., None]      # the marsh's share of the wet ground
    wet_canopy = ((ref['oasis'] * (1.0 - turn) + ref['swamp'] * turn) * (1.0 - reeds) + ref['marsh'] * reeds) * grain
    wet_cover = appearance.COVER['most'] * (1.0 - 0.6 * appearance.smoothstep(0.35, 0.7, slope))
    bare = kinds['bare'][..., None]
    plants = (1.0 - kinds['bare']) * (wetted * wet_cover + (1.0 - wetted) * cover)     # the share of the ground plants cover
    w3, c3 = wetted[..., None], cover[..., None]
    land = (w3 * (wet_canopy * wet_cover[..., None] + kinds['open'] * (1.0 - wet_cover[..., None]))
            + (1.0 - w3) * (canopy * c3 + kinds['open'] * (1.0 - c3)))
    land = kinds['bare_colour'] * bare + land * (1.0 - bare)
    waters = g['waters']
    silt = rivers['sediment'][..., None]
    clear = 0.5 * (waters['open_sea'] + waters['productive_coast'])
    river = rivers['share'][..., None]
    land = land * (1.0 - river) + (clear * (1.0 - silt) + waters['river_mouth'] * silt) * river
    # the water: a cell at the water's edge shelves from the shore down to twice its mean depth
    depth = np.where(wet, level - h, 0.0)
    edge = wet & (cv2.distanceTransform(wet.astype(np.uint8), cv2.DIST_L2, 3) <= 1.5)
    spread = 2.0 * np.maximum(depth, 0.0)[..., None] * ((np.arange(8, dtype=np.float32) + 0.5) / 8.0)
    shelf = np.where(edge, appearance.smoothstep(appearance.SHELF_M[0], appearance.SHELF_M[1], spread).mean(axis=-1),
                     appearance.smoothstep(appearance.SHELF_M[0], appearance.SHELF_M[1], depth))[..., None]
    bed = (LOOK['bed_share'] * np.where(edge, np.exp(-spread / LOOK['bed_m']).mean(axis=-1),
                                        np.exp(-np.maximum(depth, 0.0) / LOOK['bed_m'])))[..., None]
    plume = plume[..., None]
    sea = (waters['open_sea'] * (1.0 - shelf) + waters['productive_coast'] * shelf) * (1.0 - plume) + waters['river_mouth'] * plume
    sea = sea * (1.0 - bed) + soil * LOOK['bed_dark'] * bed
    colour = cv2.GaussianBlur(np.where(wet[..., None], sea, land).astype(np.float32), (0, 0), LOOK['psf_px'])
    under_water = cv2.GaussianBlur(np.maximum(wet.astype(np.float32), rivers['share']), (0, 0), LOOK['psf_px'])

    inner = (slice(margin, margin + size_r), slice(margin, margin + size_c))
    dry, wi = ~wet[inner], wet[inner]
    mean = lambda a: float(a[inner][dry].mean()) if dry.any() else 0.0       # noqa: E731
    km = stream_km = 0.0                                  # the rivers' and all the streams' length inside the window
    for k in rivers['courses']:
        y, x = k['y'], k['x']
        inside = (np.minimum(y[1:], x[1:]) >= margin) & (y[1:] < margin + size_r) & (x[1:] < margin + size_c)
        length = np.hypot(np.diff(y) * dy, np.diff(x) * np.interp(y[1:], np.arange(len(lat)), dx)) * inside / 1000.0
        km += float(length[k['flow'][1:] >= LOOK['river_m3s']].sum())
        stream_km += float(length.sum())
    ends = (held > 0.0)
    ends[:margin] = ends[-margin:] = False
    ends[:, :margin] = ends[:, -margin:] = False
    patches = held[ends] / np.maximum(deficit[ends], LOOK['deficit_least_m']) / 1e6     # km2 each hollow's water can keep wet
    q = lambda a, k: round(float(np.quantile(a, k)), 4) if a.size else 0.0   # noqa: E731
    stats = dict(tiles=[[tr + i, (tc + j) % FT.TILE_COLS] for i in (0, 1) for j in (0, 1)] if origin is None else None,
                 land_share=round(float(dry.mean()), 4),
                 rain_mm=round(mean(rain) * 1000.0, 1), demand_mm=round(mean(demand) * 1000.0, 1),
                 aridity=[round(float(np.quantile(aridity[inner], k)), 4) for k in (0.02, 0.5, 0.98)],
                 gathered_share=round(mean(share), 4), plant_cover=round(mean(plants), 4),
                 rain_fed_cover=round(mean(cover), 4), cover_on_level_deep_ground=round(mean(appearance.plant_cover(own_rain / demand)), 4),
                 water_cut_by_standing=dict(mean=round(mean(np.exp(-stands)), 3), ground_cut_at_all=round(mean((stands > 0.02).astype(np.float32)), 3),
                                            ground_cut_to_half_or_less=round(mean((stands > 0.69).astype(np.float32)), 3)),
                 wetted=round(mean(wetted), 5),
                 wetted_by_lines=round(mean(on_lines), 5), wetted_by_floors=round(mean(on_floors), 5),
                 wetted_by_water_beside=round(mean(shore), 5),
                 hollows=dict(count=int(patches.size), km2_each_median=q(patches, 0.5), km2_each_90th=q(patches, 0.9),
                              km2_largest=q(patches, 1.0), km2_together=round(float(patches.sum()), 1),
                              over_one_pixel=int((patches > dy * dy / 1e6).sum()), pixel_km2=round(dy * dy / 1e6, 3)),
                 water=dict(under_20_m=round(float((depth[inner][wi] < 20.0).mean()), 4) if wi.any() else 0.0,
                            under_3_m=round(float((depth[inner][wi] < 3.0).mean()), 4) if wi.any() else 0.0,
                            mean_shallows_share=round(float(shelf[inner][wi].mean()), 4) if wi.any() else 0.0,
                            plume_over_5_percent=round(float((plume[inner][wi] > 0.05).mean()), 5) if wi.any() else 0.0,
                            plume_most=round(float(plume[inner][wi].max()), 3) if wi.any() else 0.0),
                 river_courses=len(rivers['courses']), river_water_km2=round(float((rivers['share'] * area)[inner].sum() / 1e6), 1),
                 river_km=round(km, 1), stream_km=round(stream_km, 1))
    showing = (1.0 - kinds['bare']) * (1.0 - wetted) * (1.0 - cover)       # open ground between plants, as a share of the cell
    stats['kinds'] = {name: round(mean(v), 5) for name, v in kinds['shares'].items()}      # each kind's share of the land
    stats['kinds'].update({'planted margins and marshes of hollows': round(mean(on_floors), 5),
                           'ground kept wet along lines of gathered water': round(mean(on_lines), 5),
                           'marsh at the water\'s edge': round(mean(marsh), 5),
                           'belt over fresh groundwater, and river corridors': round(mean(shore - marsh), 5)})
    stats['ground_showing'] = dict(bare_whatever_the_climate=round(mean(kinds['bare']), 5), open_between_plants=round(mean(showing), 5))
    drawn = (aridity[ends] < K.KINDS['pan_between'][0]) & (kinds['hollow_depth'][ends] >= K.KINDS['pan_depth_m'])
    gathers = kinds['gathering_km2'][ends][drawn]                        # km2, the hollows of dry land with a floor drawn
    stats['pans'] = dict(hollows=int(ends.sum()), with_a_floor_drawn=int(gathers.size), soda=int((gathers >= K.KINDS['soda_km2']).sum()),
                         soda_share_of_floors_drawn=round(float((gathers >= K.KINDS['soda_km2']).mean()), 4) if gathers.size else 0.0,
                         depth_m=[q(np.minimum(kinds['hollow_depth'][ends], 9999.0), x) for x in (0.25, 0.5, 0.75, 0.9)],
                         gathering_km2=[q(gathers, x) for x in (0.5, 0.9, 0.916, 0.936, 0.95, 0.99, 1.0)],
                         floor_km2_largest=round(float(held[ends].max()) * K.KINDS['stand_ratio'] / K.KINDS['pool_evaporation_m'] / 1e6, 3)
                         if ends.any() else 0.0)
    return dict(colour=np.round(appearance.srgb_encode(np.clip(colour[inner], 0.0, 1.0)) * 255).astype(np.uint8),
                water=under_water[inner],
                plants=plants[inner].astype(np.float32), wetted=wetted[inner].astype(np.float32),
                river=rivers['share'][inner], stats=stats,
                kinds=dict({name: v[inner] for name, v in kinds['shares'].items()}, bare=kinds['bare'][inner],
                           margin=on_floors[inner], marsh=marsh[inner], cover=cover[inner].astype(np.float32)),
                dryness_class=np.digitize(aridity[inner], [0.05, 0.2, 0.5, 0.65]).astype(np.uint8))


def read_inputs(surface, moisture, cache=None, log=print, sand=None):
    """What the windows read: the fine surface, the moisture map's rain and demand, the rock and plumes the page's
    look keeps (fine_look.py's fields, checked against their inputs), the rock on each side of the mare map's edge and
    the map itself on the fine grid, the routing, the waters' colours, and the sand product where one is given."""
    from research.studies.atlas_round import design_routing as DRt
    started = time.time()
    record = appearance.moisture_record(moisture)
    key = FT.sha256(surface)[:16] + FT.sha256(moisture)[:16] + FT.sha256(HERE / 'appearance.py')[:16]
    with np.load(FS.BULK / 'look_fields.npz', allow_pickle=False) as z:
        if str(z['key']) != key:
            raise SystemExit('look_fields.npz was made from other inputs; run fine_look.py to make it again')
        g = {k: z[k] for k in ('plume', 'feo', 'tio2', 'mafic')}
    g.update(K.rock_sides(g['feo'], g['tio2'], g['mafic']))
    with np.load(surface, allow_pickle=False) as z:
        g['rotation'] = np.asarray(json.loads(str(z['metadata']))['frame']['rotation_today_to_turned'], float)
    g.update(ln_rain=carried(record, 'rain_mm_yr'), ln_demand=carried(record, 'demand_mm_yr'),
             storm=record['settings']['storm_water_surfaces'], waters=appearance.water_albedos())
    routing = json.loads(DRt.OUT.read_text())['product']
    if FT.sha256(routing['path']) != routing['sha256']:
        raise SystemExit(f"{routing['path']} no longer matches its record")
    g['routing'] = RC.routing_arrays(routing['path'])
    q, drained = g['routing']['q'], g['routing']['drained']
    g['streams'] = dict(g['routing'], q=np.where((q >= LOOK['river_m3s']) | (drained >= LOOK['stream_km2']), q, 0.0).astype(np.float32))
    kept = Path(cache) / f"sediment_{routing['sha256'][:16]}.npy" if cache else None
    if kept is not None and kept.exists():
        g['sediment'] = np.load(kept)
    else:
        g['sediment'] = appearance.land_discharge(routing['path']).astype(np.float32)
        if kept is not None:
            np.save(kept, g['sediment'])
    log(f'the climate, rock, routing and sediment read, {time.time() - started:.0f} s')
    with np.load(FS.PRODUCT, allow_pickle=False) as z:
        g.update(height=z['height'], water=z['water'], level=z['level'])
    log(f'the fine surface read, {time.time() - started:.0f} s')
    g['sources'] = {str(surface): FT.sha256(surface), str(moisture): FT.sha256(moisture), routing['path']: routing['sha256'],
                    str(FS.PRODUCT): FT.sha256(FS.PRODUCT)}
    from geography import fetch_surface_inputs, mare_map
    g['mare'] = mare_map.basalt(PPD, left_lon=0.0)         # the mapped outlines on the fine grid: today's Moon, from 0 E
    outlines = fetch_surface_inputs.path(mare_map.ARCHIVE)
    g['sources'][str(outlines)] = FT.sha256(outlines)
    g['sand'] = K.read_sand(sand, log) if sand else None
    if g['sand'] is not None:
        g['sources'][str(sand)] = FT.sha256(sand)
    log(f'the mare map on the fine grid, {time.time() - started:.0f} s')
    return g


def main(argv=None) -> int:
    import fine_look as FL
    from PIL import Image
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--surface', type=Path, required=True)
    ap.add_argument('--moisture', type=Path, required=True)
    ap.add_argument('--window', action='append', required=True, help='name:latitude,longitude on today\'s Moon')
    ap.add_argument('--out', type=Path, required=True, help='a folder outside the page\'s sources')
    ap.add_argument('--noise', type=float, help='the plants\' mottling (LOOK[\'noise\']); 0 leaves the ground alone to place every difference')
    ap.add_argument('--sand', type=Path, help='the sand product (sand_16ppd.npz); without it, or where it holds none, no sand is drawn')
    a = ap.parse_args(argv)
    if a.noise is not None:
        LOOK['noise'] = a.noise
    a.out.mkdir(parents=True, exist_ok=True)
    g = read_inputs(a.surface.absolute(), a.moisture.absolute(), cache=a.out, log=lambda s: print(s, flush=True), sand=a.sand)
    for w in a.window:
        started = time.time()
        name, at = w.rsplit(':', 1)
        lat, lon = (float(v) for v in at.split(','))
        tr, tc = FL.window_tiles(lat, lon)[0]
        out = look_window(tr, tc, g)
        folder = a.out / name.lower().replace(' ', '_')
        folder.mkdir(exist_ok=True)
        Image.fromarray(out['colour']).save(folder / 'colour.png')
        for k in ('water', 'plants', 'wetted', 'river'):
            Image.fromarray(np.round(np.clip(out[k], 0.0, 1.0) * 255).astype(np.uint8)).save(folder / f'{k}.png')
        np.savez_compressed(folder / 'kinds.npz', dryness_class=out['dryness_class'],     # 0 hyper-arid to 4 humid; the shares as 0-255
                            **{k: np.round(np.clip(v, 0.0, 1.0) * 255).astype(np.uint8) for k, v in out['kinds'].items()})
        out['stats'].update(name=name, today=[lat, lon], seconds=round(time.time() - started, 1), settings=LOOK, kinds_settings=K.KINDS, covers=COVERS,
                            reference_gain=REFERENCE_GAIN, sources=g['sources'])
        (folder / 'record.json').write_text(json.dumps(out['stats'], indent=1, default=list) + '\n')
        print(name, json.dumps({k: v for k, v in out['stats'].items() if k not in ('settings', 'kinds_settings', 'sources', 'tiles', 'covers')}), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
