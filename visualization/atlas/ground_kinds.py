"""The kinds of ground under the close-up look, and its layer for loose sand.

ground_look.py mixes a 474 m pixel from open water, ground kept wet, rain-fed plants and open ground. This module says
what the open ground is. The kinds are the soils study's (the research drive's session-tools/agents/soils/
soils_maintenance.md of 3 November 2026, sections 1.2, 1.3, 2.1 to 2.6 and 3.1 to 3.5), and each is placed by what the
ground gives at that scale: the fine heights, the water, the mare map's outlines, the way storm water runs over the
ground, and the climate's cells. No kind is scattered by noise, and a kind smaller than a pixel enters as its share of
the pixel. KINDS holds every setting and says of each whether it comes from the study, from a product or is a display
guess. Each colour starts from the page's soil colours (appearance.soil_colour) and departs from them by a Munsell
neutral, the page's basalt rock, or a colour or a tone measured in the reference mosaics (MEASURED).

Where each kind lies:

- Soil by dryness class, on level ground (study 2.1 to 2.5). In the arid land the page's desert soil with lime in it,
  under a biological crust that darkens it (study 2.2). From the semi-arid land on a steppe soil with dark humus
  under its crust (Kastanozems), drawn in the colour of Earth's steppe at that dryness; in the dry sub-humid land the
  deep humus soils of meadow steppe (Chernozems, Phaeozems), darker still; in humid land the page's forest soils, pale
  where highland rock is leached (study 2.5) and red on basalt.
- Place on a plain (plain_ground, provinces): how far a cell of a plain lies below or above the plain round it,
  within about 15 km, a crater's floor weighed against itself and its rim and walls left out, so that the swales and
  swells between the craters are drawn and single craters are left alone. Low ground gathers the finer, deeper soil:
  sandy and cemented where arid (Calcisols, Durisols and Arenosols where water gathers), dark and fine from the
  semi-arid land on (Vertisols in lowlands; study 2.4). A swell keeps the thinner soil: towards stony desert where
  arid, the page's paler soil beyond, with a little less water for its plants. The study says that what slopes lose
  settles in the low ground close by (2.4, 2.6); the scales and strengths are display guesses.
- Stony ground under a lag of rock fragments (Leptosols, Regosols; study 2.1, 2.6): all hyper-arid land, in the colour
  of Earth's stony deserts, fading out through the dry end of the arid class; and steep ground in every class, in
  its own soil's colour, darker by the measure that rough ground is darker than smooth ground in mosaics of Earth's
  dry country.
- Storm water's ground. The gathered water of ground_look.py is followed down the fine ground with the sediment it
  can carry (through_flow). A line leaves its marks once it drains KINDS['line_km2'], which leaves out the gullies
  of a single slope, or KINDS['plain_km2'] on gentle ground off the swells where it joins such a line: a bare
  sandy bed, wider with its catchment; fans where the flow slackens and its load settles; and a floor of alluvial
  ground beside it as far as the ground stays within a metre or two of its bed (valley_floors). That ground is the
  low ground's: sandy and cemented where arid (study 2.1, 2.2, 3.4), dark and fine beyond. Its sand is its catchment's: pale
  feldspar and glass off highland rock, dark off basalt.
- Hollows (study 2.6, 3.3 to 3.5). The water a hollow holds is the seepage through its floor, which the plants round
  it draw on: a planted margin over the area that water meets the ground's deficit, as before. The water that stands
  and dries on the floor is KINDS['stand_ratio'] times that, and it covers, as a pool that evaporates
  KINDS['pool_evaporation_m'] a year, the floor laid from the hollow's lowest cell upward, and the margin rings it.
  In the arid and hyper-arid land the floor is bare: pale clay cemented by calcite and opal. From the semi-arid land
  on the hollows are planting pits (study 2.3) and the floor is planted. Where a hollow gathers the most ground's
  water the floor is a soda pan, bare in every dry class. A floor is drawn in hollows at least KINDS['pan_depth_m']
  deep (hollow_depths).
- Cracking clay on basalt (Vertisols on the mare remnants; study 1.2, 2.3, 2.4): gentle ground inside the mare map's
  outlines in the semi-arid and dry sub-humid land, dark and grey, under grass without trees. The outlines are
  rasterised on the fine grid, so the soil boundary at a mare's edge is as sharp as the map draws it, and the rock's
  iron and titanium are read on each side of it separately (rock_sides).
- Fresh silt (Fluvisols; study 2.5): bare bars beside a river in proportion to its load, and bare flats on the low
  ground beside the plume at its mouth.
- Loose sand, from the sand product when one is given (sand_window): its share of each cell covers everything else
  on land; the active part is bare, and the held part carries the cover its climate gives over sand-coloured ground.
  Where the product's crest spacing is two pixels or more the sand is drawn as crests with the ground between them
  (crests): texture laid along the product's direction at its spacing, the place of each crest a fixed function of
  position on the Moon and otherwise made up.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import appearance  # noqa: E402

LUMINANCE = np.array([0.2126, 0.7152, 0.0722], np.float32)
NEUTRAL_8 = (0.31006, 0.31616, 59.10)    # Munsell N 8/ as xyY under Illuminant C (the renotation's value)
REFERENCE_GAIN = 1.6                     # as ground_look's: a measured colour's reflectance is the colour over the close tone's gain
MEASURED = dict(                         # sRGB, colour groups of NASA's Blue Marble: Next Generation at about 550 m a pixel
    steppe=(108, 92, 62),                # semi-arid steppe on dark-humus soil: eastern Mongolia's largest group, its median
    meadow=(86, 75, 47),                 # its darker steppe, towards the wetter north of the scene
    wet_meadow=(69, 65, 36),             # and its darkest ground short of forest
    pavement=(135, 106, 74),             # stony desert: the mean of Lake Eyre's stony plains (141, 100, 62), the Messak plateau
    #                                      (130, 101, 73) and the Sistan basin's dark gravel plain (135, 116, 86)
    pavement_basalt=(124, 95, 62),       # the Harrat ash Shaam's basalt plain
    pan_floor=(190, 175, 149),           # a dry lake's clay floor: the Sistan basin's (the Helmand's end)
    soda=(204, 191, 174))                # a salt crust: Lake Eyre's and its neighbours'; soda's trona and natron are as white
SAND_SHAPE = (2880, 5760)                # the sand product's grid: 16 cells per degree from 90 N and 0 E, today's Moon
SAND_FIELDS = ('sand_share', 'active_share', 'form', 'crest_azimuth_deg', 'crest_spacing_m', 'drift_azimuth_deg')
SAND_MADE = ('made_sea', 'dune_height_m', 'made_sand_part', 'made_sand_basalt')   # what a product with made seas adds (sand_seas.py)
DUNES = dict(                            # display: the relief drawn for a made sea's dunes; the product gives form, direction, spacing and height
    pad_cells=8,                         # product cells (15 km) kept round a sea's outline, where its fringe may carry dunes
    star_ring=0.18,                      # star dunes: the spread of the pattern's wavelengths about the spacing, as a share of the wavenumber
    star_floor=0.35,                     # and how far the ground between the dunes dips, against a dune's height
    star_peak=2.4,                       # a star dune of the product's height stands where the pattern is this many of its standard deviations
    ridge_along=(0.16, 0.11),            # transverse ridges, linear dunes: the wavenumber's spread along the pattern's own direction
    ridge_across=(0.38, 0.10),           # and across it: wide for sinuous, broken crests, narrow for long straight ones
    lee=0.25,                            # transverse ridges: the second harmonic that steepens the lee face, where the pixels resolve it
    crest=0.15,                          # linear dunes: the second harmonic that sharpens the crest
    amplitude=(0.35, 1.6),               # a dune's height against its field's mean, held between these
    streak_km=(12.0, 2.5, 35.0),         # darker and paler sand: streaks down the drift of this length and width, on drifts of tone this broad
    streak_parts=(0.8, 0.6),             # the streaks' and the broad drifts' parts in the pattern
    streak=0.08,                         # how far the pattern moves the basalt's share of a made sea's sand, as its standard deviation
    margin_soft=0.12,                    # where a made sea thins out its sand is drawn on the highest dunes first; how soft their edges are, in rank
    tone=0.04)                           # crests are drawn paler and the ground between darker by this share, at a dune of the usual height
KINDS = dict(
    # the rock on each side of a mare's mapped edge
    rock_cells=8.0,             # the page's iron and titanium (16 px/deg) averaged over this many cells on each side (15 km)
    mare_feo_least=15.0,        # the page's rule (appearance.rock_fields): mare basalt soils hold at least 15% FeO
    # stony ground
    pavement=(0.05, 0.10),      # study 2.1: the hyper-arid class (aridity under 0.05) under a pavement; guess: gone by 0.10
    pavement_cover=0.7,         # guess: the share of such ground that stones cover
    steep=(0.2, 0.5),           # guess: ground steepening between these slopes (11 to 27 degrees) turns stony,
    steep_cover=0.8,            # to this share
    rough_tone=0.8,             # measured: steep stony ground is its soil's colour at this share of its brightness. In eight mosaics
    #                             of Earth's dry country the roughest fifth of the land is 0.65 to 0.84 as bright as the smoothest
    # a place on its plain: swales and swells between the craters
    plain_px=5, plain_relief_m=300.0,   # guess: a plain is ground whose heights span less than this within 2.4 km (five cells):
    #                             a little over half of the two dry highland windows
    plain_km=15.0,              # guess: a cell of a plain is set against the plain within about this distance,
    plain_band_m=150.0,         # that lies within about this of its own height: a crater's floor is its own plain, and its rim
    #                             and walls are none, so that swales and swells are drawn and single craters are left alone
    plain_m=(8.0, 35.0),        # guess: ground that far below the plain round it is low ground, from none to all; as far above
    #                             it, a swell (nine tenths of those windows' plains lie within 40 m of the plain round them)
    swell_tone=0.6,             # guess: how far a full swell's ground takes the stony or thin soil's colour
    swell_ln=0.22,              # guess: a swell's thin soil cuts its plants' water by up to this (the bound of ground_look's
    #                             rule for ground that stands out, with which it shares LOOK['stand_most'])
    # the biological crust
    crust_from=(0.03, 0.10),    # study 2.1, 2.2: little crust in the hyper-arid land, a full one in the arid;
    crust_to=(0.5, 0.65),       # study 2.4: a closed cover with litter takes over in the dry sub-humid land
    crust_cover=0.85,           # guess: its share of the stable open ground (study 2.2: it holds the wind above two thirds)
    crust_dark=0.6,             # guess: a crusted soil's reflectance over the bare soil's (dark cyanobacterial and lichen
    #                             crusts lower a soil's reflectance by a third to a half on Earth; recalled, not re-read)
    crust_on_alluvium=0.5,      # guess: alluvial ground, worked by each storm, keeps this much of that crust
    # storm water through the ground
    settle=0.3,                 # guess: of the load beyond what the flow can carry, the share that settles in each cell
    fan_index=3.0,              # guess: the settled load, in cells of upstream ground times slope lost, that covers 63% of a cell
    line_km2=(40.0, 160.0),     # guess: a line of storm water leaves its marks (a bed, fans, alluvial ground beside it) once it
    #                             drains the first of these areas, and in full from the second: lines that cross the country,
    #                             and none of the gullies of a single slope;
    plain_km2=(5.0, 20.0),      # guess: and from these areas on gentle ground off the swells, where the line belongs to a
    #                             network: the water it joins drains line_km2 by its end
    plain_slope=(0.03, 0.08),   # guess: gentle is under the first of these slopes, and none of it over the second
    wash_m=12.0,                # guess: its bare bed and bars, metres wide per square root of a km2 drained
    valley_m=25.0,              # guess: the alluvial ground beside it, this many metres to either side per square root of a km2 drained,
    valley_most_m=3000.0,       # at most this far,
    valley_rise_m=1.5,          # on ground within this of its bed, and as much again for each tenfold of its catchment
    wet_fade=(0.5, 1.0),        # guess: bare beds and open alluvium give way to covered ground between these aridities
    # hollows
    stand_ratio=1.48,           # study 3.1, 3.3: the crusted surface sheds 7.2% of the drylands' rain and 2.9% (2.7 to
    #                             3.1) seeps through the hollows' floors, so (7.2 - 2.9) / 2.9 times the seepage stands
    #                             and dries; guess: the same ratio in every class and hollow
    pool_evaporation_m=0.84,    # study 3.5: a shallow pool evaporates 68 mm a lunar day, 0.84 m a year
    pan_between=(0.5, 0.8),     # guess: standing water dries to a bare floor under these aridities and keeps a marsh above
    pit_between=(0.15, 0.3),    # study 2.3 and the moisture map's soil zones: the semi-arid land's hollows are planting pits,
    #                             their floors planted; guess: bare floors give way to planted ones between these aridities.
    #                             A floor that turns to soda stays bare in every dry class
    floor_m=1.5,                # guess (ground_look's): a hollow's floor is the ground within this of its lowest cell
    pan_depth_m=5.0,            # guess: a bare floor is drawn in hollows at least this deep; shallower pits in the fine heights
    #                             lie within the marks its instrument's tracks leave on level ground
    margin_m=12.0,              # study 1.3: roots draw on a water table from about 12 m above it in a sandy loam
    soda_km2=8.5,               # guess, set against study 3.5 (6.4-8.4% of the depressions turn saline within ten thousand years
    #                             on clay floors): a floor turns to soda where its hollow holds the storm water of more ground
    #                             than this, the solutes of the rain and dust on that ground with it
    # clays, humus and leached soil
    vertisol=(0.12, 0.25),      # study 2.3, 2.4: cracking clay on mare ground from the semi-arid land (guess: from these
    vertisol_to=(0.65, 1.0),    # aridities) to the dry sub-humid; red clay beyond (the page's 2.5YR 3/6)
    vertisol_slope=(0.03, 0.08),   # guess: on ground gentler than these slopes
    clay_value=0.0639,          # the luminance of Munsell value 3 (Y 6.555%): a very dark grey clay (typically 10YR 3/1)
    clay_chroma=0.5,            # guess: with half the chroma of the page's soil colour there
    steppe_soil=(0.15, 0.25),   # study 2.3: Kastanozems, dark humus under a crust, across the semi-arid class (from 0.2); guess:
    #                             coming in between these aridities. The ground takes the steppe reference's colour
    meadow_soil=(0.4, 0.57),    # study 2.4: Chernozems and Phaeozems in the dry sub-humid class (from 0.5); darker still
    humus_to=(1.0, 2.0),        # guess: giving way to the forest soils of wet land
    leached=(0.65, 1.5),        # study 2.5: highland soils of humid land stay pale; guess: over these aridities, taking the
    #                             page's palest highland colour (7.5YR 6/4)
    # fresh silt
    bar=0.5, bar_px=1.0,        # guess: bare bars cover half as much ground as a fully laden river's water, within a pixel
    mouth=0.5, mouth_px=3.0,    # guess: half of the low wet ground beside a plume is bare silt; the plume read within 3 pixels
    # colours
    sand_white=0.1,             # guess: highland sand is the rock's least weathered soil colour with this share of N 8/
    #                             (white feldspar); the rest of its tone is glass and iron stain
    sand_rock=0.5,              # guess: basalt sand is that soil colour with this share of the page's basalt rock
    hauled_highland=(5.0, 0.3),  # FeO and TiO2 (weight per cent) of the highland rock a made sea's pale sand is hauled from: the highlands' usual 4-6% FeO
    # loose sand
    dune_px=2.0,                # crests are drawn where their spacing is at least this many pixels
    dune_soft=0.25,             # display: how soft a crest's edge is drawn, in the cosine of the phase
    anchor_deg=0.25,            # the crests' phase is anchored on a lattice of this spacing fixed on the Moon
    sand_rock_px=2.0)           # guess: wind-laid sand takes the colour of the rock within this many pixels


def S(a, b, x):
    return appearance.smoothstep(a, b, x)


def between(lo_hi, aridity_ln):
    """A share rising from 0 to 1 as the aridity crosses a pair of bounds, in its logarithm."""
    return S(np.log(lo_hi[0]), np.log(lo_hi[1]), aridity_ln)


def rock_sides(feo, tio2, mafic, cells=None):
    """The rock's iron and titanium on each side of the mare map's edge, from the page's 16 px/deg fields (FeO and TiO2
    in weight per cent, and the basalt's share): each side's own cells averaged over KINDS['rock_cells'], so that ground
    just inside a mare reads the basalt's oxides and ground just outside it the highland's. The basalt's iron is at
    least KINDS['mare_feo_least'], as on the page."""
    import cv2
    k = KINDS['rock_cells'] if cells is None else cells
    pad = int(np.ceil(3 * k))
    blur = lambda a: cv2.GaussianBlur(np.concatenate([a[:, -pad:], a, a[:, :pad]], axis=1), (0, 0), k)[:, pad:-pad]   # noqa: E731
    m = np.clip(np.asarray(mafic, np.float32), 0.0, 1.0)
    out = {}
    for name, f in (('feo', feo), ('tio2', tio2)):
        f = np.asarray(f, np.float32)
        for side, w in (('mare', m), ('highland', 1.0 - m)):
            weight = blur(w)
            out[f'{name}_{side}'] = np.where(weight > 1e-3, blur(f * w) / np.maximum(weight, 1e-3), f).astype(np.float32)
    out['feo_mare'] = np.maximum(out['feo_mare'], KINDS['mare_feo_least'])
    return out


def through_flow(s, wet, rcv, own, slope, basalt, settle=None):
    """Storm water followed down the receivers with the sediment it carries. own is the water each cell sheds (m3 a
    year), slope its ground's and basalt the share of its ground on basalt. A flow can carry a load of its water
    times its slope; a cell adds its own water times its slope, so on an even slope the load keeps to what the flow
    can carry and nothing settles. Where the load exceeds it (the slope slackens), `settle` of the excess is laid
    down in the cell and the rest passes on; an end with no receiver takes all that reaches it. Returns the water
    through each cell (its own included), the load laid down in it (m3 a year times slope) and the share of its water
    that fell on basalt."""
    settle = KINDS['settle'] if settle is None else settle
    n = s.size
    order = np.argsort(s, axis=None)[::-1]
    order = order[~wet.ravel()[order]]
    r, o, g = rcv.ravel().tolist(), np.asarray(own, np.float64).ravel().tolist(), np.asarray(slope, np.float64).ravel().tolist()
    d = (np.asarray(own, np.float64) * basalt).ravel().tolist()
    water, load, dark, laid = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
    for i in order.tolist():
        q = water[i] + o[i]
        b = dark[i] + d[i]
        carried = load[i] + o[i] * g[i]
        water[i], dark[i] = q, b
        over = carried - q * g[i]
        if over > 0.0:
            over *= settle
            laid[i] = over
            carried -= over
        j = r[i]
        if j >= 0:
            water[j] += q
            dark[j] += b
            load[j] += carried
        else:
            laid[i] += carried
    water = np.asarray(water).reshape(s.shape)
    return water, np.asarray(laid).reshape(s.shape), np.asarray(dark).reshape(s.shape) / np.maximum(water, 1e-12)


def laid(s, root, first, second, area, floor_m, margin_m):
    """Two areas laid over each hollow from its lowest cell upward, as shares of each cell. `first` and `second` are
    areas (m2) read at each hollow's lowest cell (root gives every cell's): the first covers the lowest ground within
    floor_m of that cell and no more than that floor holds; the second covers the next lowest ground round it, within
    margin_m."""
    flat = s.ravel()
    z = flat - flat[root]
    b1, b2 = np.asarray(first, np.float64).ravel()[root], np.asarray(second, np.float64).ravel()[root]
    one, two = np.zeros(s.size, np.float32), np.zeros(s.size, np.float32)
    cand = np.flatnonzero(((b1 > 0.0) | (b2 > 0.0)) & (z <= max(floor_m, margin_m)))
    if cand.size:
        order = cand[np.lexsort((z[cand], root[cand]))]
        a = area.ravel()[order].astype(np.float64)
        upto = np.cumsum(a)
        group = root[order]
        start = np.r_[True, group[1:] != group[:-1]]
        before = upto - a - np.maximum.accumulate(np.where(start, upto - a, 0.0))
        on_floor = z[order] <= floor_m
        floor = np.bincount(group, weights=a * on_floor, minlength=s.size)[group]
        f = np.minimum(b1[order], floor)
        inner = np.clip((f - before) / a, 0.0, 1.0)
        one[order] = inner * on_floor
        two[order] = np.where(z[order] <= margin_m, np.clip((f + b2[order] - before) / a, 0.0, 1.0) - inner, 0.0)
    return one.reshape(s.shape), two.reshape(s.shape)


def valley_floors(surf, wet, upstream_km2, cell, slope, marked=None):
    """The share of each cell that is the alluvial floor of the nearest marked line (`marked`, how far each cell is one,
    0 to 1; without it, the lines draining KINDS['line_km2']): the
    part of its ground within reach of that line's bed in height (KINDS['valley_rise_m'], and as much again for each
    tenfold of the line's catchment; a cell's ground spans its slope times its width about its mean height) and in
    distance (KINDS['valley_m'] to either side per square root of a km2 drained, at most KINDS['valley_most_m'])."""
    import cv2
    k = KINDS
    if marked is None:
        marked = S(k['line_km2'][0], k['line_km2'][1], upstream_km2)
    source = (marked > 0.0) & ~wet
    if not source.any():
        return np.zeros(surf.shape, np.float32)
    dist, label = cv2.distanceTransformWithLabels((~source).astype(np.uint8), cv2.DIST_L2, 5, labelType=cv2.DIST_LABEL_PIXEL)
    index = np.zeros(int(label.max()) + 1, np.int64)
    index[label[source]] = np.flatnonzero(source)
    near = index[label]
    z = surf - surf.ravel()[near]
    drained = np.maximum(upstream_km2.ravel()[near], 1.0)
    reach = np.minimum(k['valley_m'] * np.sqrt(drained), k['valley_most_m'])
    rise = k['valley_rise_m'] * (1.0 + np.clip(np.log10(drained / k['line_km2'][0]), 0.0, None))
    span = np.maximum(slope * cell, 0.5)
    part = np.clip((rise - (z - 0.5 * span)) / span, 0.0, 1.0)
    within = np.where(dist < 0.5, np.clip(2.0 * reach / cell, 0.0, 1.0), S(reach, 0.3 * reach, dist * cell))
    return np.where(wet, 0.0, part * within * np.asarray(marked, np.float32).ravel()[near]).astype(np.float32)


def ground_blur(a, km, lat, dy):
    """A Gaussian blur over a distance on the ground: km north-south and, in each whole degree of latitude, km east-west
    by that degree's own column spacing, so that every window blurs a place alike."""
    import cv2
    a = np.ascontiguousarray(a, np.float32)
    sigma = km * 1000.0 / dy
    out = cv2.GaussianBlur(a, (1, 2 * int(np.ceil(3.0 * sigma)) + 1), sigma)
    band = np.floor(lat)
    for degree in np.unique(band):
        rows = np.flatnonzero(band == degree)
        across = sigma / max(float(np.cos(np.radians(degree + 0.5))), 1e-3)
        out[rows[0]:rows[-1] + 1] = cv2.GaussianBlur(out[rows[0]:rows[-1] + 1], (2 * int(np.ceil(3.0 * across)) + 1, 1), across)
    return out


def plain_ground(surf, wet, lat, dy):
    """How far each cell of a plain lies below or above the plain round it (m); 0 off the plains. A plain is ground
    whose heights (smoothed over a cell) span less than about KINDS['plain_relief_m'] within KINDS['plain_px'] cells,
    which leaves out rims, walls and the small craters. A cell is set against the mean height of the plain within about
    KINDS['plain_km'] that lies within about KINDS['plain_band_m'] of its own height, so that the floor of a crater
    is weighed against itself and the country round the crater against the country. Worked on every second cell in
    bands of height fixed in metres, with distances on the ground: a fixed function of the place."""
    import cv2
    k = KINDS
    z = cv2.GaussianBlur(np.asarray(surf, np.float32), (0, 0), 1.0)
    reach = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * k['plain_px'] + 1,) * 2)
    span = cv2.dilate(z, reach) - cv2.erode(z, reach)
    plain = (1.0 - S(0.7 * k['plain_relief_m'], 1.3 * k['plain_relief_m'], span)) * ~np.asarray(wet, bool)     # how far a cell is plain
    zs, ps, lats = z[::2, ::2], plain[::2, ::2].astype(np.float32), np.asarray(lat)[::2]
    band = k['plain_band_m']
    total, weight = np.zeros(zs.shape, np.float32), np.zeros(zs.shape, np.float32)
    for level in np.arange(np.floor(float(zs.min()) / band) * band, float(zs.max()) + band, band):
        near = np.clip(1.0 - np.abs(zs - level) / band, 0.0, 1.0)
        if float((near * ps).max()) > 0.0:
            total += near * ground_blur(near * ps * zs, k['plain_km'], lats, 2.0 * dy)
            weight += near * ground_blur(near * ps, k['plain_km'], lats, 2.0 * dy)
    rel = np.where(weight > 0.02, zs - total / np.maximum(weight, 1e-6), 0.0).astype(np.float32) * ps
    rel = cv2.resize(cv2.GaussianBlur(rel, (0, 0), 1.0), (z.shape[1], z.shape[0]), interpolation=cv2.INTER_LINEAR)
    return (rel * plain).astype(np.float32)


def at_fine(field, rows, cols, step=4, nearest=False):
    """A field kept on every `step`-th cell of the fine grid (the sand product) read at a window's fine rows and
    columns: bilinear, or its nearest cell's, wrapping in longitude."""
    import cv2
    rows, cols = np.asarray(rows), np.asarray(cols)
    y = (rows + 0.5) / step - 0.5
    x = (cols[0] + np.arange(len(cols)) + 0.5) / step - 0.5
    r = np.clip(np.arange(int(np.floor(y[0])), int(np.floor(y[-1])) + 2), 0, field.shape[0] - 1)
    c = np.arange(int(np.floor(x[0])), int(np.floor(x[-1])) + 2)
    crop = np.nan_to_num(np.asarray(field[r[0]:r[-1] + 1][:, c % field.shape[1]], np.float32))
    mx = np.repeat((x - c[0]).astype(np.float32)[None, :], len(rows), axis=0)
    my = np.repeat(np.clip(y - r[0], 0, len(r) - 1).astype(np.float32)[:, None], len(cols), axis=1)
    return cv2.remap(crop, mx, my, cv2.INTER_NEAREST if nearest else cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)


def provinces(broad):
    """A place on its plain as two shares, low ground and swell, from how far it lies below or above the plain round
    it (plain_ground's field, m) by KINDS['plain_m']. Fixed scales in metres."""
    b = KINDS['plain_m']
    broad = np.asarray(broad, np.float32)
    return S(-b[0], -b[1], broad), S(b[0], b[1], broad)


def hollow_depths(s, root):
    """How deep the hollow each cell drains to is (m): from its lowest cell up to the lowest pass on its edge, where its
    water would spill into a neighbouring hollow or reach open water. A hollow that meets no other within the heights
    given has no bound."""
    lab = root.reshape(s.shape)
    spill = np.full(s.size, np.inf)
    whole, less, more = slice(None), slice(None, -1), slice(1, None)
    for a, b in (((whole, less), (whole, more)), ((less, whole), (more, whole)), ((less, less), (more, more)), ((less, more), (more, less))):
        differ = lab[a] != lab[b]
        over = np.maximum(s[a], s[b])[differ]
        np.minimum.at(spill, lab[a][differ], over)
        np.minimum.at(spill, lab[b][differ], over)
    return (spill - s.ravel())[root].reshape(s.shape)


def measured(name):
    """A colour measured in the reference mosaics, as the reflectance (linear sRGB) the close tone shows as that colour."""
    v = np.asarray(MEASURED[name], np.float32) / 255.0
    return np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4).astype(np.float32) / REFERENCE_GAIN


def greyed(colour, value, chroma):
    """A colour brought to a luminance, keeping `chroma` of its difference from the grey of that luminance."""
    lum = np.maximum(colour @ LUMINANCE, 1e-4)[..., None]
    return (chroma * colour * (value / lum) + (1.0 - chroma) * value).astype(np.float32)


def read_sand(path, log=print):
    """The sand product's fields, or None where the file is absent or holds no sand."""
    path = Path(path)
    if not path.exists():
        log(f'no sand product at {path}: the sand layer draws nothing')
        return None
    with np.load(path, allow_pickle=False) as z:
        missing = [k for k in SAND_FIELDS if k not in z.files]
        if missing or z['sand_share'].shape != SAND_SHAPE:
            raise SystemExit(f'{path} is not the sand product: missing {missing}, sand_share {z["sand_share"].shape}')
        product = {k: z[k] for k in SAND_FIELDS}
        if all(k in z.files for k in SAND_MADE):
            product.update({k: z[k] for k in SAND_MADE})
            product['seas'] = json.loads(str(z['metadata'])).get('seas', [])
    return product if float(np.nanmax(product['sand_share'])) > 0.0 else None


def sand_window(product, rows, cols):
    """The sand product at a window's fine cells (rows and columns of the 64 px/deg grid): the shares and the crest
    spacing bilinear, the crests' direction bilinear as an axis (a crest running north runs south as well), the form
    its nearest cell's. None where the product is absent or the window holds no sand."""
    if product is None:
        return None
    read = lambda name, nearest=False: at_fine(np.asarray(product[name]), rows, cols, 4, nearest)   # noqa: E731
    share = read('sand_share')
    if float(share.max()) <= 0.0:
        return None
    dunes = np.nan_to_num(np.asarray(product['crest_spacing_m'], np.float32)) > 0.0
    azimuth = np.radians(np.nan_to_num(np.asarray(product['crest_azimuth_deg'], np.float32)))
    turn = lambda f: at_fine(f(2.0 * azimuth) * dunes, rows, cols, 4)   # noqa: E731
    out = dict(share=np.clip(share, 0.0, 1.0), active=np.clip(read('active_share'), 0.0, 1.0), form=read('form', True).astype(np.uint8),
               azimuth_deg=np.degrees(0.5 * np.arctan2(turn(np.sin), turn(np.cos))),
               spacing_m=read('crest_spacing_m') * (at_fine(dunes.astype(np.float32), rows, cols, 4) > 0.5))
    if 'made_sand_part' in product:                          # a made sea's sand, and the basalt's share of it by where it was hauled from
        part = np.clip(read('made_sand_part'), 0.0, 1.0)
        dark = at_fine(np.asarray(product['made_sand_part'], np.float32) * np.asarray(product['made_sand_basalt'], np.float32), rows, cols, 4)
        out.update(made=part, made_basalt=np.clip(dark / np.maximum(part, 1e-4), 0.0, 1.0))
    return out


def crests(lat, lon, azimuth_deg, spacing_m):
    """The cosine of the dunes' phase at each cell: 1 on a crest, -1 midway between two. Crests run along azimuth_deg
    (clockwise from north) spacing_m apart. The phase is anchored on a lattice fixed on the Moon, a point every
    KINDS['anchor_deg'] of latitude and longitude: each of the four anchors round a cell proposes the plane wave of
    the cell's direction and spacing that has a crest through the anchor, and the four are summed as phasors with
    bilinear weights. The result depends on the position and on the two fields alone, so every window draws the same
    crests; where neighbouring anchors disagree the crests fork, as dunes do. Texture: the product gives a direction
    and a spacing, and the place of each crest is made up."""
    a = KINDS['anchor_deg']
    la = np.asarray(lat, np.float64)[:, None] + np.zeros((1, len(lon)))
    lo = np.asarray(lon, np.float64)[None, :] + np.zeros((len(lat), 1))
    u, v = la / a - 0.5, lo / a - 0.5
    i0, j0 = np.floor(u), np.floor(v)
    fu, fv = u - i0, v - j0
    theta = np.radians(np.asarray(azimuth_deg, np.float64))
    across_east, across_north = np.cos(theta), -np.sin(theta)
    wave = 2.0 * np.pi / np.maximum(np.asarray(spacing_m, np.float64), 1.0)
    radius = appearance.MOON_KM * 1000.0
    total = np.zeros(la.shape, np.complex128)
    for di in (0, 1):
        for dj in (0, 1):
            weight = (fu if di else 1.0 - fu) * (fv if dj else 1.0 - fv)
            east = radius * np.cos(np.radians(la)) * np.radians(lo - (j0 + dj + 0.5) * a)
            north = radius * np.radians(la - (i0 + di + 0.5) * a)
            total += weight * np.exp(1j * wave * (east * across_east + north * across_north))
    return np.cos(np.angle(total)).astype(np.float32)


def dune_relief(product, log=print):
    """The relief of each made sea's dunes on the fine grid, as slopes to add to the ground's: a list of patches, each
    with its first fine row and column and the dunes' height's slope eastward and northward (metres a metre). The
    product gives each sea's form, the direction of its net drift, its crest spacing and the dunes' height cell by
    cell; where each dune stands is made up here, as texture. A field of random numbers, the same at every run for the
    same outline, is passed through a filter that keeps the wavelengths near the spacing: every direction alike for
    star dunes, which leaves mounds a spacing apart; a band across the drift for transverse ridges, with a steeper lee
    face; a narrow band along it for linear dunes. Distances are measured on the ground at the sea's own latitude."""
    if product is None or not product.get('seas') or 'made_sea' not in product:
        return []
    k, made = DUNES, np.asarray(product['made_sea'])
    dy = appearance.MOON_KM * 1000.0 * np.radians(1.0 / (4 * 16))
    patches = []
    for sea in product['seas']:
        cells = np.argwhere(made == sea['id'])
        if not len(cells):
            continue
        (ra, ca), (rb, cb) = cells.min(axis=0), cells.max(axis=0)
        if cb - ca > made.shape[1] // 2:
            log(f"the sea {sea['name']} lies across 0 E: its dunes are not drawn")
            continue
        r0, r1 = max(ra - k['pad_cells'], 0), min(rb + k['pad_cells'] + 1, made.shape[0])
        c0, c1 = max(ca - k['pad_cells'], 0), min(cb + k['pad_cells'] + 1, made.shape[1])
        rows, cols = np.arange(4 * r0, 4 * r1), np.arange(4 * c0, 4 * c1)
        height = at_fine(np.asarray(product['dune_height_m'], np.float32), rows, cols, 4)
        if float(height.max()) <= 0.0:
            continue
        lat = 90.0 - (rows + 0.5) / 64.0
        dx = dy * math.cos(math.radians(sea['centre_today'][0]))
        spacing, towards = float(sea['crest_spacing_m']), math.radians(sea['net_towards_deg'])
        form = 'star' if sea['form'].startswith('star') else 'linear' if sea['form'].startswith('linear') else 'transverse'
        rng = np.random.default_rng([int(sea['id']), int(rows[0]), int(cols[0])])
        white = np.fft.fft2(rng.standard_normal((len(rows), len(cols))))
        k_north = -2.0 * np.pi * np.fft.fftfreq(len(rows), d=dy)[:, None]      # rows run south
        k_east = 2.0 * np.pi * np.fft.fftfreq(len(cols), d=dx)[None, :]
        along = k_east * math.sin(towards) + k_north * math.cos(towards)       # the wavenumber along the drift, and across it
        across = -k_east * math.cos(towards) + k_north * math.sin(towards)
        k0 = 2.0 * np.pi / spacing
        if form == 'star':
            field = np.fft.ifft2(white * np.exp(-(np.hypot(along, across) - k0) ** 2 / (2.0 * (k['star_ring'] * k0) ** 2))).real
            field /= k['star_peak'] * field.std()                              # a dune of the usual size stands a little under 1
            pattern = np.clip(np.where(field > 0.0, np.abs(field) ** 1.25, k['star_floor'] * field), -1.0, k['amplitude'][1] * 0.75)
        else:
            i = 0 if form == 'transverse' else 1
            own, other = (along, across) if form == 'transverse' else (across, along)
            wave = np.fft.ifft2(white * np.exp(-((own - k0) ** 2 / (2.0 * (k['ridge_along'][i] * k0) ** 2) + other ** 2 / (2.0 * (k['ridge_across'][i] * k0) ** 2))))
            size = np.clip(np.abs(wave) / np.abs(wave).mean(), *k['amplitude'])
            phase = np.angle(wave)                                             # it grows downwind for transverse ridges
            if form == 'transverse':
                pixel = math.hypot(dx * math.sin(towards), dy * math.cos(towards))
                shape = np.cos(phase) + k['lee'] * float(np.clip((spacing / pixel - 3.0) / 3.0, 0.0, 1.0)) * np.sin(2.0 * phase)
            else:
                shape = np.cos(phase) + k['crest'] * np.cos(2.0 * phase)
            pattern = 0.5 * size * shape                                       # from trough to crest, a dune of the usual size is the product's height
        relief = (height * pattern).astype(np.float32)
        north = -np.gradient(relief, axis=0) / dy
        east = np.gradient(relief, axis=1) / (dy * np.cos(np.radians(lat)))[:, None]
        long, wide, broad = (1000.0 * v for v in k['streak_km'])     # sand from different sources: streaks drawn out down the drift, on broad drifts of tone
        smooth = lambda shape: np.fft.ifft2(np.fft.fft2(rng.standard_normal(white.shape)) * shape).real   # noqa: E731
        streak = smooth(np.exp(-0.5 * ((along * long) ** 2 + (across * wide) ** 2)))
        drift = smooth(np.exp(-0.5 * (along ** 2 + across ** 2) * broad ** 2))
        streak = k['streak_parts'][0] * streak / streak.std() + k['streak_parts'][1] * drift / drift.std()
        rank = np.empty(pattern.size, np.float32)                    # each pixel's place among the patch's, lowest ground to highest dune
        rank[np.argsort(pattern, axis=None)] = (np.arange(pattern.size) + 0.5) / pattern.size
        patches.append(dict(row=int(rows[0]), col=int(cols[0]), east=east.astype(np.float32), north=north.astype(np.float32), relief=relief,
                            streak=(streak / streak.std()).astype(np.float32), grain=pattern.astype(np.float32) * (height > 0),
                            rank=rank.reshape(pattern.shape), covered=np.ones(pattern.shape, np.float32), name=sea['name'], form=form))
        log(f"the dunes of {sea['name']}: {form}, {spacing:.0f} m apart, up to {float(np.abs(relief).max()):.0f} m of relief on {len(rows)} x {len(cols)} pixels")
    return patches


def relief_window(patches, rows, cols):
    """The made seas' dunes at a window's fine rows and columns: their slopes eastward and northward (`east`, `north`),
    their relief against their height (`grain`), each pixel's place from the lowest ground to the highest dune (`rank`),
    the streaks and drifts of tone (`streak`) and where a sea's patch lies (`covered`), or None where the window
    holds none."""
    out = None
    for p in patches or ():
        nr, nc = p['east'].shape
        ra, rb = max(rows[0], p['row']), min(rows[-1] + 1, p['row'] + nr)
        if ra >= rb:
            continue
        inside = (cols >= p['col']) & (cols < p['col'] + nc)
        if not inside.any():
            continue
        if out is None:
            out = {name: np.zeros((len(rows), len(cols)), np.float32) for name in ('east', 'north', 'grain', 'streak', 'rank', 'covered')}
        take = (slice(ra - p['row'], rb - p['row']), cols[inside] - p['col'])
        for name, field in out.items():
            field[ra - rows[0]:rb - rows[0], inside] += p[name][take]
    return out


def sand_cover(sand, lat, lon, dy, dunes=None):
    """The share of each cell under loose sand as drawn: the product's share, gathered onto crests where the crests
    stand KINDS['dune_px'] pixels or more apart (crescent, transverse, linear and held dunes along one direction; star
    dunes along two at right angles), and laid evenly elsewhere. Over a field of crests the mean stays the product's.
    Where a made sea's dunes are given (relief_window) its sand lies on them: all of it where the sea is whole, and
    on the highest dunes first where it thins out, so that its margin breaks up into dunes on bare ground."""
    share = sand['share']
    made = None
    if dunes is not None and sand.get('made') is not None:
        soft = DUNES['margin_soft']
        at = (1.0 - share) * (1.0 + soft) - 0.5 * soft               # the rank above which sand lies: under every rank where the sea is whole
        made = (dunes['covered'] > 0.5) & (sand['made'] > 0.5)
        on_dunes = S(at - 0.5 * soft, at + 0.5 * soft, dunes['rank'])
    drawn = (sand['spacing_m'] >= KINDS['dune_px'] * dy) & (sand['form'] >= 2)
    if made is not None:
        drawn &= ~made
    if not drawn.any():
        return share if made is None else np.where(made, on_dunes, share).astype(np.float32)

    def gathered(part, azimuth):
        soft = KINDS['dune_soft'] * np.clip(4.0 * part * (1.0 - part), 1e-3, 1.0)
        return np.clip(0.5 + (crests(lat, lon, azimuth, sand['spacing_m']) - np.cos(np.pi * part)) / (2.0 * soft), 0.0, 1.0)

    one = gathered(share, sand['azimuth_deg'])
    star = sand['form'] == 4
    if star.any():
        part = 1.0 - np.sqrt(1.0 - share)
        two = 1.0 - (1.0 - gathered(part, sand['azimuth_deg'])) * (1.0 - gathered(part, sand['azimuth_deg'] + 90.0))
        one = np.where(star, two, one)
    out = np.where(drawn, one, share)
    return (out if made is None else np.where(made, on_dunes, out)).astype(np.float32)


def ground(surf, s, wet, rcv, root, own, held, deficit, least, area, cell, slope, aridity, mare, rock, soil,
           broad=None, laden=None, mouth=None, sand=None, lat=None, lon=None, dy=None, dunes=None):
    """The open ground of a window, kind by kind.

    surf and s are the ground's heights (s with ground_look's lean for level ground), wet the open water, rcv and
    root the receivers and the ends they lead to, own the storm water each cell sheds and held what each end holds (m3
    a year), deficit the ground's yearly lack of water (m) and least the least it is spread against, area and cell
    the cells' areas and widths, slope and aridity the ground's and the climate's, mare the mare map on the fine grid,
    rock the iron and titanium on each side of its edge (rock_sides' four fields at the fine cells), soil the page's
    soil colour there, broad the place of each cell on its plain (plain_ground's field), laden the share of each cell
    under sediment-laden river water, mouth the share of low wet ground beside a plume, sand the sand product at
    the fine cells (sand_window) with lat, lon and dy to lay it, and dunes a made sea's dunes there (relief_window).

    Returns the hollows' planted margins and marshes as shares of each cell (`margin`), the share of each cell that is
    bare whatever its climate (`bare`: beds, pan floors, bars and moving sand) and its colour (`bare_colour`), the
    colour of the open ground between the plants on the rest (`open`), the share of that rest on which grass grows
    without trees (`grassland`), the plains' swells (`swell`), and each kind's share of the cell (`shares`)."""
    import cv2
    k = KINDS
    ln = np.log(aridity)
    mare = np.asarray(mare, bool)
    dry = 1.0 - between(k['wet_fade'], ln)
    own = np.maximum(np.asarray(own, np.float64), 1e-9)
    low, swell = provinces(broad) if broad is not None else (np.zeros(aridity.shape, np.float32),) * 2

    # the storm water's ground: beds, fans and the floors of the lines that cross the country or net a plain
    flow, settled, basalt = through_flow(s, wet, rcv, own, slope, mare)
    upstream_km2 = flow / own * area / 1e6
    gentle = (1.0 - S(k['plain_slope'][0], k['plain_slope'][1], slope)) * (1.0 - swell)
    network = S(k['line_km2'][0], k['line_km2'][1], upstream_km2.ravel()[root].reshape(aridity.shape))     # by the end of its way
    marked = np.maximum(S(k['line_km2'][0], k['line_km2'][1], upstream_km2),
                        S(k['plain_km2'][0], k['plain_km2'][1], upstream_km2) * gentle * network)
    wash = np.clip(k['wash_m'] * np.sqrt(upstream_km2) / cell, 0.0, 1.0) * marked * dry
    fan = (1.0 - np.exp(-settled / own / k['fan_index'])) * marked
    alluvial = (1.0 - (1.0 - fan) * (1.0 - valley_floors(surf, wet, upstream_km2, cell, slope, marked))).astype(np.float32)

    # the hollows: a floor where the storm water stands and dries, and the planted margin round it
    depth = hollow_depths(s, root)
    standing = k['stand_ratio'] * held / k['pool_evaporation_m']
    floor, margin = laid(s, root, standing * (depth >= k['pan_depth_m']), held / np.maximum(deficit, least), area, k['floor_m'], k['margin_m'])
    dries = (1.0 - S(k['pan_between'][0], k['pan_between'][1], aridity)).ravel()[root].reshape(aridity.shape)
    gathering_km2 = held / (own / area) / 1e6                                     # the ground whose storm water a hollow holds
    soda = (gathering_km2.ravel()[root] >= k['soda_km2']).reshape(aridity.shape)
    planted = np.where(soda, 0.0, between(k['pit_between'], ln).ravel()[root].reshape(aridity.shape))
    pan = floor * dries * (1.0 - planted)
    margin = margin + floor - pan

    # fresh silt beside laden rivers and at their mouths
    silt = np.zeros(aridity.shape, np.float32)
    if laden is not None:
        silt += np.clip(k['bar'] * cv2.GaussianBlur(np.asarray(laden, np.float32), (0, 0), k['bar_px']), 0.0, 1.0)
    if mouth is not None:
        silt += k['mouth'] * np.asarray(mouth, np.float32)
    silt = np.where(wet, 0.0, np.clip(silt, 0.0, 1.0))

    # loose sand, from the product
    loose = sand_cover(sand, lat, lon, dy, dunes) if sand is not None else np.zeros(aridity.shape, np.float32)
    loose = np.where(wet, 0.0, loose)
    moving = loose * (sand['active'] if sand is not None else 0.0)
    held_sand = loose - moving

    # the bare ground, in turn: moving sand over everything, then the pan floors, the bars and the beds
    rest = 1.0 - loose
    pan, rest = pan * rest, rest * (1.0 - pan)
    silt, rest = silt * rest, rest * (1.0 - silt)
    wash, rest = wash * rest, rest * (1.0 - wash)
    bare = moving + pan + silt + wash

    # colours of the bare kinds
    least_weathered = np.full(aridity.shape, 0.05, np.float32)
    fresh_highland = appearance.soil_colour(least_weathered, rock['feo_highland'], rock['tio2_highland'])
    fresh_mare = appearance.soil_colour(least_weathered, rock['feo_mare'], rock['tio2_mare'])
    sand_highland = (1.0 - k['sand_white']) * fresh_highland + k['sand_white'] * appearance.munsell_srgb(NEUTRAL_8)
    sand_mare = (1.0 - k['sand_rock']) * fresh_mare + k['sand_rock'] * appearance.BASALT_ROCK
    carried = basalt[..., None].astype(np.float32)
    water_sand = sand_highland * (1.0 - carried) + sand_mare * carried            # sand that storm water brought
    near = cv2.GaussianBlur(mare.astype(np.float32), (0, 0), k['sand_rock_px'])[..., None]
    wind_sand = sand_highland * (1.0 - near) + sand_mare * near                   # sand that the wind laid
    if sand is not None and sand.get('made') is not None:                         # a made sea's sand keeps the colour of the rock it was hauled from:
        dark, part = sand['made_basalt'][..., None], sand['made'][..., None]      # highland sand of the highlands' usual make-up, basalt sand of the rock at hand
        if dunes is not None:                                                     # display: streaks down the drift move the basalt's share
            dark = np.clip(dark + DUNES['streak'] * dunes['streak'][..., None], 0.0, 1.0)
        far = appearance.soil_colour(least_weathered, np.full(aridity.shape, k['hauled_highland'][0], np.float32), np.full(aridity.shape, k['hauled_highland'][1], np.float32))
        hauled = (1.0 - k['sand_white']) * far + k['sand_white'] * appearance.munsell_srgb(NEUTRAL_8)
        wind_sand = wind_sand * (1.0 - part) + (hauled * (1.0 - dark) + sand_mare * dark) * part
        if dunes is not None:                                                     # display: crests a little paler than the ground between them
            wind_sand = wind_sand * (1.0 + DUNES['tone'] * dunes['grain'][..., None] * part)
    m3 = mare[..., None]
    fresh_here = np.where(m3, fresh_mare, fresh_highland)
    pan_colour = np.where(soda[..., None], measured('soda'), measured('pan_floor'))
    bare_colour = ((moving[..., None] * wind_sand + pan[..., None] * pan_colour + silt[..., None] * fresh_here
                    + wash[..., None] * water_sand) / np.maximum(bare, 1e-6)[..., None]).astype(np.float32)

    # the soil between the plants on the rest of the ground, by dryness class and by its place on the broad ground
    mix = lambda a, b, share: a * (1.0 - share[..., None]) + b * share[..., None]   # noqa: E731
    crusted = k['crust_cover'] * between(k['crust_from'], ln) * (1.0 - between(k['crust_to'], ln))
    shade = (1.0 - crusted * (1.0 - k['crust_dark']))[..., None]
    steppe, meadow, woods = between(k['steppe_soil'], ln), between(k['meadow_soil'], ln), between(k['humus_to'], ln)
    leached = between(k['leached'], ln) * ~mare
    plain = mix(soil, fresh_highland, leached) * shade                            # the page's soil, pale where leached, under its crust
    level = mix(mix(mix(plain, measured('steppe'), steppe), measured('meadow'), meadow), plain, woods)
    vertisol = mare * between(k['vertisol'], ln) * (1.0 - between(k['vertisol_to'], ln)) * (1.0 - S(k['vertisol_slope'][0], k['vertisol_slope'][1], slope))
    level = mix(level, greyed(soil, k['clay_value'], k['clay_chroma']) * shade, vertisol)
    deep = mix(mix(measured('meadow'), measured('wet_meadow'), meadow), plain, woods)        # the low ground's deeper, finer, darker soil
    deep = mix(deep, greyed(soil, k['clay_value'], k['clay_chroma']) * shade, vertisol)
    lowland = mix(water_sand * (1.0 - k['crust_on_alluvium'] * (1.0 - shade)), deep, steppe)  # sandy and cemented where arid, dark and fine beyond
    stones = np.where(m3, measured('pavement_basalt'), measured('pavement'))
    upland = mix(level, mix(stones, plain, steppe), np.full(aridity.shape, k['swell_tone'], np.float32))   # a swell: towards stony
    #                                                                               desert where arid, the page's paler soil beyond
    open_ground = level * (1.0 - low - swell)[..., None] + lowland * low[..., None] + upland * swell[..., None]
    open_ground = mix(open_ground, lowland, alluvial)
    steep = (k['steep_cover'] * S(k['steep'][0], k['steep'][1], slope)).astype(np.float32)
    open_ground = mix(open_ground, level * k['rough_tone'], steep)
    pavement = (k['pavement_cover'] * (1.0 - S(k['pavement'][0], k['pavement'][1], aridity))).astype(np.float32)
    open_ground = mix(open_ground, stones, pavement)
    under = held_sand / np.maximum(1.0 - bare, 1e-6)                           # the held sand's share of the ground left
    open_ground = mix(open_ground, wind_sand * shade, under)

    left = (1.0 - bare) * (1.0 - under)
    here = left * (1.0 - pavement) * (1.0 - steep) * (1.0 - alluvial)             # ground that keeps its own place's soil
    on_level = here * (1.0 - low - swell)
    humus = steppe * (1.0 - woods) * (1.0 - vertisol)
    shares = {'moving sand': moving, 'sand held under plants or crust': held_sand, 'pan floor, cemented clay': pan * ~soda,
              'soda pan': pan * soda, 'fresh silt': silt, 'bare wash bed': wash,
              'stony ground: pavement and steep slopes': left * (1.0 - (1.0 - pavement) * (1.0 - steep)),
              'swells of the plains: stony where arid, paler soil beyond': here * swell,
              'low ground of the plains: sandy and cemented where arid, dark and fine beyond': here * low,
              'alluvial ground beside lines, sandy and cemented': left * alluvial * (1.0 - steppe),
              'alluvial ground beside lines, dark and fine': left * alluvial * steppe,
              'cracking clay on basalt': on_level * vertisol,
              'steppe soil, dark humus under a crust': on_level * humus * (1.0 - meadow),
              'meadow soil, deep dark humus': on_level * humus * meadow,
              'desert soil with lime under a crust': on_level * (1.0 - steppe) * (1.0 - vertisol),
              'pale leached soil on highland rock': on_level * leached * woods,
              'ground under a biological crust': left * (1.0 - pavement) * crusted}
    shares = {name: np.where(wet, 0.0, v).astype(np.float32) for name, v in shares.items()}
    return dict(margin=margin.astype(np.float32), bare=np.where(wet, 0.0, bare).astype(np.float32), bare_colour=bare_colour,
                open=open_ground.astype(np.float32), grassland=(vertisol * (1.0 - alluvial)).astype(np.float32), swell=swell, low=low,
                shares=shares, soda=soda, hollow_depth=depth, gathering_km2=gathering_km2, upstream_km2=upstream_km2.astype(np.float32),
                basalt_water=basalt.astype(np.float32))
