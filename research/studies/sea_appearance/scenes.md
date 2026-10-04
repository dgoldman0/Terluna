# The four coasts, scene by scene

Twenty-four views of the Open Moon's seas: four coasts at the six moments of the
[results by regime](README.md#the-seas-by-regime), from the Sun at its highest to the
darkest hour of the night. Each scene is the frame the
[sea-scene renderer](../../../visualization/reference-renderer/seas/) will draw. They are described here
in words and numbers, so that they can be pictured, or illustrated with an image
generator, before the full renderings are made.

Every number comes from [scenes.json](results/scenes.json), written by
[scenes.py](scenes.py) from the study's products:
- the colours and brightness of sky and sea from the regime panoramas (the solved sky, the
  rough-water reflection, both disks' glitter, the water's own light);
- the coast from LOLA heights at the moment's tide;
- the waves from that hour's wave spectrum;
- the stars from the bright-star catalogue.

These are calculated scenes with stated gaps; an image made from them illustrates the
calculation. Two kinds of number are approximate:
- **Land colours:** a placeholder ground, lit by the clear sky and hazed by the air.
- **Sea colours:** the panoramas assume open water in every direction, so where land rises
  behind the sea, the water just below the horizon also mirrors the land, as each scene notes.

## Reading the scenes

- **The frame** is a photograph 1920 by 1080 pixels with a 65-degree horizontal field of view (a 28 mm
  lens on a full-frame camera; 39.4 degrees vertically), taken from a small boat with the eye 2 m above
  the water, the camera level and pitched as stated. Pixel positions count from the top left. The sea's
  horizon lies 2.6 km away and 0.087 degrees below level; its row is given for the frame's centre.
- **Colours** are sRGB hex values at one exposure per frame: the luminance shown white is stated (it is
  1.25 times the brightest twentieth of the sky in the frame), brighter light rolls off smoothly toward
  white from 80% of it, and the display's white is Earth daylight (D65), without adapting to the light of
  the moment. That is a photograph's rendering: in the night and late-evening scenes it is a long
  exposure, and the eye sees less colour than the values show. Luminances are in cd/m².
- **Directions** are compass bearings. Wave directions say where the waves travel toward, and how that
  lies against the view.

## What every scene shares

- **Clear air, seven times Earth's column.** The air is 1.2 atmospheres at the surface and very tall.
  The sky is cloudless in every scene. The calculated air holds no dust or haze droplets, so the Sun
  stands in a broad, smooth glow without a bright halo round it.
- **The Sun's colour.** Through this much air the Sun is pale gold even high in the sky. It is orange at
  15 degrees up and deep orange near the horizon. Its disk is 0.53 degrees across, about 16 pixels.
- **The daytime sky** is a very pale cyan-blue, almost white, becoming a pale grey-green-cream toward the
  horizon. The sea under it is a grey-teal slate.
- **The evening lasts days.** The Sun sinks half a degree an hour. With it 4 degrees below the horizon
  the western sky is amber and orange. At 15 degrees below the western sky is amber and coral and the
  eastern a soft golden cream. At 35 degrees below, three to four days after sunset, a deep orange arch
  still glows over the western horizon.
- **The horizon is close and the air hazy.**
  - Distant hills keep two-thirds of their contrast at 30 km in green light and 40% in blue.
  - At 65 km they keep 40% and 15%, so far ranges fade toward the colour of the sky at the horizon.
  - The Moon's curvature hides the lowest metres of coasts beyond 2.6 km.
- **The Earth.** It is 1.9 to 2.0 degrees across (54 to 60 pixels), nearly four times the Sun's width. It
  hangs nearly still in the nearside sky, shows phases, and never rises on the far side.
  - Above the air it is bluish white.
  - Through the air it is amber to orange, deeper the lower it stands.
  - The calculation gives it no visible continents or clouds; it is a smooth lit shape. The unlit part
    of its disk looks like the surrounding sky.
- **The sea.** Deep water carries long swell from the wave runs, covered with wind ripples where the wind
  blows. At lunar gravity the smallest wind ripples are 4.2 cm long, coarser than Earth's 1.7 cm, and the
  wind raises them harder. In light winds the sea goes glassy.
  - The water is the biosphere's productive coast: looking straight down into it, it is olive green.
    Seen at low angles the sea is mostly reflected sky.
  - The water is unbroken: no foam or whitecaps, which the study does not yet model.
- **The land is a placeholder:** bare ground with the colour of Apollo soils, dark grey-brown mare soil,
  and a brighter tan highland soil at the Ingenii coast. There are no plants, buildings, roads, boats or
  people.
- **Few stars.** Even the darkest hours are bright. The earthlit or twilit sky is 0.03 to 0.8 cd/m²,
  100 to 3,000 times Earth's moonless night, and the thick air dims stars near the horizon by two to four
  magnitudes. Only one frame records stars.

## Western Oceanus Procellarum, by Russell

The most exposed coast of the study. The camera floats at the wave runs' shore point (73.875° W,
28.875° N) in 78 m of water. Low islets 20 to 60 m high lie 1.8 to 2.5 km to the east-northeast. The
mainland shore runs 3 to 8 km away to the south and west, behind it a belt of hills and crater rims
rising to 650 to 1,360 m at 15 to 42 km. Long ocean swell from the open sea reaches this coast. The Earth
stands in the east all month, 6 to 21 degrees up.

### WP-1 · The Sun at its highest

*Day 18.3 of the month (from 7 February 2038), seven and a half days before sunset. Sun 59.9° up in
the south. 72,500 lux on the ground; day vision.*

**Image prompt.** Photograph from a small boat on a wide sea at midday under a cloudless sky, eye 2 m
above the water, 28 mm lens, looking south-southwest, camera tilted 3 degrees down so the horizon sits
43% from the top of the frame. The high Sun is out of frame, ahead and to the left, so the light falls
from high in front: the water shows no glitter in view. The sky is a very pale cyan-blue, almost white
(#cdece7 at the top), a little brighter than near the horizon, where it turns a pale grey-green-cream
(#ced8c5) and is evenly bright from side to side. Along the whole horizon runs a low, faint range of
hills and crater rims 16 to 42 km away, 25 to 63 pixels tall, highest at three-quarters of the way across. The hills are
a hazy, desaturated blue-grey (#8a96a1), only slightly darker than the sky above them, sunlit but softened
by the long path through the air. The sea is a grey-teal slate (#859895 just below the horizon), much
darker than the sky, darkening smoothly to #565f57 at the bottom of the frame with a faint olive tint in
the nearest water. A long ocean swell of 2.2 m significant height (crest to trough) and 71 m between crests,
runs away from the camera toward the hills. Its crests lie parallel to the horizon as broad, low, rounded
ridges, so the near water rises and falls in long gentle slopes. Everything is covered with lively wind
ripples a few centimetres to a few decimetres long, coarser than ripples on Earth. A 5 m/s breeze blows
from behind the camera. The ripples break the reflected sky into a fine, high-contrast mosaic of
pale and slate facets; the waves show strongly. Unbroken water, no foam. No boats, no people, no plants,
no buildings.

**Frame data.**
- **View:** 205° (south-southwest), pitched 3° down; horizon at row 463. The frame is 39% sky,
  4% land and 57% sea.
- **Exposure:** white is 14,100 cd/m².
- **Sky:**
  - top #cdece7 (11,100 cd/m²);
  - 8° up #d2e7d9;
  - 3° up #d1dfcd;
  - at the horizon #ced8c5 (9,330);
  - the left edge 7 to 9% brighter than the right.
- **Sun and Earth:**
  - The Sun is 59.9° up at 178.6°, out of frame above the top left. Its disk is pale gold #ffd979,
    750 million cd/m².
  - The Earth is 6.5° up in the east-northeast, behind the camera, 40% lit.
- **Coast:** a continuous band of hills across the full width, 25 to 63 px above the horizon line
  (1.0 to 2.2°):
  - left third 25 to 36 px, at 34 to 42 km;
  - the centre rising from 30 to 55 px;
  - highest near x = 1500 (63 px, 1,400 m at 30 km);
  - falling to 44 to 45 px at the right edge, where a nearer ridge stands at 16 km.
  - The low shore 3 to 8 km away forms the foot of the band just above the sea horizon, less hazed than
    the hills; the Moon's curvature hides its lowest 0.1 to 8 m.
  - Sunlit faces; haze transmittance 0.67 (green) and 0.41 (blue); colour about #8a96a1.
- **Sea:**
  - below the horizon #859895 (4,180 cd/m², 45% of the sky above);
  - 3° down #7c8d8a;
  - 8° down #70807b;
  - bottom #565f57 (1,530).
- **Waves:**
  - swell: significant height 2.2 m, period 16.6 s, 71 m long, travelling toward 207°, straight away
    from the camera;
  - wind 4.9 m/s toward 207°;
  - mean square slope 0.037 (rms tilt 10.9°), of which 0.025 from waves shorter than 1 m;
  - from the nearside wave run, 4 h after the moment.
- **Water:** its own light looking down is olive (chromaticity 0.34, 0.41), 204 cd/m².
- **Tide:** 0.9 m below its mean level.

### WP-2 · The Sun 4° up

*Day 25.3, nine hours before sunset. Sun 3.7° up in the west. 6,700 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the late afternoon, eye 2 m above the water,
28 mm lens, looking west into the low Sun, horizon 53% from the top. The Sun is a small, intensely
bright orange disk, 16 pixels wide (#ff9b00 at its own exposure, here blown out to white with a warm
bloom), standing 3.7 degrees above the horizon just right of centre (x = 1078, y = 470). It hangs above
a dark, hazy line of coastal hills. Beneath it a broad, dazzling orange glitter path, 100 to 270 pixels
wide, runs from the horizon down to the bottom of the frame. It is made of countless sparkling highlights
on a rough, rippled sea, brightest at the horizon. The whole sky is warm:
- pale amber-cream at the top (#ffe59c), slightly brighter above the Sun;
- deepening through #ffd484 at 8 degrees and #ffc17a at 3 degrees;
- to warm orange at the horizon (#f4b474).

There is no halo: the glow is broad and smooth. The coast across the whole width of the frame, 4 to 8 km
away with hills 15 to 30 km behind it, is a low backlit silhouette 17 to 61 pixels tall, highest left of
centre. It is dull red-brown (#7d5d51), softened by haze. Away from the glitter the sea is a warm
tan-bronze (#b49661 below the horizon), darkening to dark olive-brown (#6d5e40) in the lower corners.
A long, low swell 2.9 m high and 150 m from crest to crest travels away to the left. Its crests run
diagonally, and the glitter path wobbles and breaks across them. Wind ripples, coarser than on Earth,
cover everything under a 4 m/s breeze. Unbroken water, no foam, no boats, no people.

**Frame data.**
- **View:** 262°, pitched 1° up; horizon at row 569. The frame is 49% sky, 4% land and 47% sea.
- **Exposure:** white is 4,290 cd/m².
- **Sun:**
  - at (1078, 470), 3.67° up at 266.5°;
  - its disk is 37 million cd/m² (8,600 times white), deep orange #ff9b00.
- **Glitter:** from the horizon to the bottom, across x = 965 to 1235. Its peak at the horizon is
  95,000 cd/m² (22 times white), vivid orange; it covers 10% of the sea in the frame.
- **Sky:** at the top centre #ffe59c (3,440 cd/m²) and at the edges #f6d895 (3,050), then #ffd484 at 8°,
  #ffc17a at 3° and #f4b474 at the horizon (2,270).
- **Coast:**
  - from 52 px at the left edge (29.8 km) to the highest ridge near x = 480 (61 px, 647 m at 15.4 km);
  - 31 to 59 px across the centre (hills 320 to 650 m at 14 to 18 km);
  - lower on the right (17 to 41 px, with an islet ridge 93 m high at 8 km).
  - The shore is 4.3 to 10.7 km away.
  - Backlit; haze transmittance 0.82 (green) and 0.64 (blue); colour #7d5d51.
- **Sea beside the glitter:** #b49661 (1,390 cd/m²) below the horizon, #aa8c5b at 3° down, and #6d5e40
  to #756342 at the bottom corners.
- **Waves:**
  - swell 2.9 m, period 24.4 s, 153 m long, travelling toward 222°, 40° left of straight away;
  - wind 4.3 m/s;
  - mean square slope 0.031 (rms tilt 9.9°);
  - spectrum 21 h before the moment.
- **Earth:** 12.4° up in the east behind the camera, 98% lit.
- **Tide:** +0.25 m.

### WP-3 · The Sun 4° below the horizon

*Day 26.1, nine hours after sunset. 2,770 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the early evening twilight, eye 2 m above the
water, 28 mm lens, looking west toward where the Sun has set, horizon 57.5% from the top. The Sun itself
is 4.35 degrees below the horizon. The sky is a luminous warm gradient, brightest above the place of
sunset just right of centre:
- soft amber at the top (#ffe194);
- light orange (#ffce78) at 8 degrees and #ffb76f at 3 degrees;
- coral-orange at the horizon (#f0a76b).

The amber fills the frame from the horizon to its top, 980 to 1,650 cd/m². A low, hazy silhouette of coastal hills, dull red-brown (#79554a), lies along the horizon across
the frame. It is 14 to 62 pixels tall, highest in the left half (630 m hills at 15 km), lower and fainter
on the right. The sea is a choppy, wind-roughened surface reflecting the amber sky diffusely, at about
60% of its brightness: warm tan-orange (#af8a56 just below the horizon), darkening to deep bronze
(#745c39) at the bottom. A long swell 3.3 m high and 150 m from crest to crest travels away to the left
under a 5 m/s wind. Ripples cover every slope and break the reflection into a fine texture of orange
highlights and darker brown facets, with no distinct glitter path. Unbroken water, no foam, no boats,
no people.

**Frame data.**
- **View:** 268°, pitched 3° up; horizon at row 621. The frame is 54% sky, 4% land and 43% sea.
- **Exposure:** white is 2,080 cd/m².
- **Sky:**
  - top #ffe194 (1,650 cd/m²);
  - 15° up #ffde87;
  - 8° up #ffce78;
  - 3° up #ffb76f;
  - at the horizon #f0a76b (984);
  - the edges 8 to 12% darker than the centre.
- **Coast:**
  - left half 45 to 62 px (hills 490 to 710 m at 15 to 23 km);
  - right half 14 to 39 px (300 to 900 m at 17 to 42 km).
  - The shore is 4.3 to 10.6 km away.
  - Backlit; colour #79554a.
- **Sea:**
  - below the horizon #af8a56 (581 cd/m²);
  - 3° down #a58150;
  - 8° down #917247;
  - bottom #745c39 (244).
- **Waves:**
  - swell 3.3 m, period 24.4 s, 153 m long, travelling toward 223°, 45° left of straight away;
  - wind 5.0 m/s;
  - mean square slope 0.042 (rms tilt 11.6°), the roughest moment at this coast.
  - At this wind a few whitecaps would form on Earth; the study has none.
- **Earth:** 13.4° up in the east behind the camera, full.
- **Tide:** +0.39 m.

### WP-4 · The Sun 15° below the horizon

*Day 27.1, a day and a half after sunset. 479 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the long evening, eye 2 m above the water,
28 mm lens, looking east, away from the sunset, horizon 57.5% from the top. A full Earth hangs in the
sky just left of centre, 15 degrees above the horizon (centre x = 952, y = 226). It is a round disk
54 pixels across, bright amber-gold in colour (#ffc600 at its own exposure), here overexposed to a warm
white with a soft bloom. The sky is a soft, even, pale golden cream (#e7e7cc at the top, #d8cdac just
above the horizon), a fifth duller and greyer toward the horizon, where the Moon's own shadow
rises into the air. A row of low, flat islets 1.8 to 2.5 km away runs from the left edge to about
nine-tenths of the way across. They are 12 to 60 m high, 8 to 43 pixels tall, highest on the left, and
their bases sit on the water just under the horizon. They are dark muted brown silhouettes (#695847),
lit from behind by the Earth. The sea is a muted grey-gold (#9a9580 to #b6a080 just below the horizon),
slightly warmer and brighter in a broad soft column below the Earth. It darkens to olive-grey (#726b58)
at the bottom. There is no sharp glitter path: the earthlight is faint against the reflected sky. A long
swell 3.3 m high, 150 m from crest to crest, comes toward the camera from the front left. Its crests run
diagonally. Wind ripples under a 4.7 m/s breeze give the water a fine, matte texture. Unbroken water,
no foam, no boats, no people.

**Frame data.**
- **View:** 92° (east), pitched 3° up; horizon at row 621. The frame is 55% sky, 2% land and 43% sea.
- **Exposure:** white is 99 cd/m².
- **Earth:**
  - at (952, 226), 14.8° up at 91.7°;
  - 53.5 px across, 99.8% lit;
  - 905 cd/m² (9 times white), amber-gold #ffc600; above the air it would be bluish white.
- **Sky:**
  - top #e7e7cc (78 cd/m²);
  - 15° up #e8e3c4;
  - 8° up #e5dcbb;
  - 3° up #ded3b2;
  - at the horizon #d8cdac (61);
  - the left edge 3% brighter than the right.
- **Islets:**
  - 16 px at the left edge, 43 px near x = 480, 30 to 32 px across the centre, 21 px at x = 1560,
    8 px at x = 1680, sea beyond x = 1740;
  - 1.75 to 2.5 km away, 12 to 60 m high; colour #695847.
- **Sea:**
  - below the horizon #b6a080 at the centre (36 cd/m²) and #9a9580 at the edges (30);
  - 8° down #90846c;
  - bottom #726b58 (15).
  - The Earth's glitter stays below half the reflected sky everywhere in the frame.
- **Waves:**
  - swell 3.3 m, period 24.4 s, 153 m long, travelling toward 223° (131° from the view, so coming from
    the front left);
  - wind 4.7 m/s;
  - mean square slope 0.039.
- **Tide:** +0.58 m.

### WP-5 · The Sun 35° below the horizon

*Day 29.0, three and a third days after sunset. 9.6 lux; between day and night vision.*

**Image prompt.** Photograph from a small boat at sea late in the long evening, eye 2 m above the water,
28 mm lens, looking west-northwest, horizon 57.5% from the top; a long exposure. A great twilight arch
glows over the western horizon. It is a broad dome of deep gold and orange centred a little right of the
middle of the frame:
- brightest amber-gold about 8 degrees up (#ffcc42);
- vivid orange at the horizon (#ff9d32);
- duller ochre-orange at the top of the frame (#ecbc6b);
- redder and darker toward the left edge (#f8873f low down).

The Sun is 35 degrees below the horizon, and this light comes from the very high air it still reaches.
A hazy silhouette of coastal hills, dark red-brown (#7c422f), runs along the horizon:
- tallest at the far left, 66 pixels (640 m hills at 16 km);
- sloping down to a thin, faint line of distant ranges only 6 to 15 pixels tall on the right;
- with a few low islands near the right edge.

The sea mirrors the arch: glowing orange (#ff9f36) just below the horizon at the centre, darker toward
the edges (#cb833c), and deep brown-orange (#8a5e28) at the bottom. A smooth swell 1.8 m high and 105 m
between crests runs away to the left. Light wind ripples under a 2.6 m/s breeze soften the reflection
into a warm, slightly grainy sheen. To the eye at this light level colours are muted; the photograph
records them. No stars, no foam, no boats, no people.

**Frame data.**
- **View:** 285°, pitched 3° up; horizon at row 621. The frame is 55% sky, 3% land and 43% sea.
- **Exposure:** white is 16.5 cd/m².
- **Sky:**
  - top #ecbc6b (9.0 cd/m²);
  - 15° up #ffcc5d;
  - 8° up #ffcc42 (13.4, the brightest);
  - 3° up #ffb433;
  - at the horizon #ff9d32 (9.8) at the centre;
  - the left edge #f8873f (6.2) and the right edge #ff9239 (7.9).
  - The arch is centred toward 291°, at x = 1124.
- **Coast:**
  - 66 px at the left edge (643 m at 15.8 km), 55 to 60 px to x = 240;
  - 29 to 41 px to x = 960;
  - 6 to 19 px beyond (ranges at 18 to 32 km and islets at 9 km).
  - Backlit; colour #7c422f.
- **Sea:**
  - below the horizon #ff9f36 (7.9 cd/m²) at the centre and #cb833c at the left edge;
  - 8° down #cc7f2e;
  - bottom #8a5e28 (2.2).
- **Waves:**
  - swell 1.8 m, period 20.1 s, 105 m long, travelling toward 221°, 64° left of the view;
  - wind 2.6 m/s;
  - mean square slope 0.034 (rms tilt 10.4°), mostly from the swell.
- **Earth:** 17.2° up in the east behind the camera, 96% lit; its glitter path is on the sea behind.
- **Stars:** none recorded.
- **Tide:** +0.97 m.

### WP-6 · The darkest hour

*Day 4.5, eight days after sunset; the Sun 60° below the horizon. 0.48 lux, from the Earth and the
twilight; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat on a glassy night sea, eye 2 m above the
water, 28 mm lens, looking east toward the Earth, horizon 60% from the top. The Earth hangs near the top
centre, 20 degrees up (centre x = 970, y = 100). It is a gibbous disk 54 pixels across, 58% lit, its
bright side toward the lower left. It is golden-amber in colour (#ffcd2d at its own exposure) but here
overexposed to near white with a bloom, the unlit part invisible against the sky. Straight down from it
a vivid golden-yellow glitter path, up to 260 pixels wide, runs from the horizon to the bottom of the
frame (#ffe93c at its brightest). Its edge is soft and its interior a ladder of long smooth highlights,
because the sea is glassy. The sky is pale and warm:
- cream-yellow overhead (#dae0c6);
- a soft peach-orange glow along the left horizon (#ffc596), on the side toward the north, where the
  twilight of the deep-set Sun lingers;
- a greyer, duller yellow-beige toward the right (#c9b796).

A row of low flat islets 2 km away, 20 to 60 m high and 15 to 43 pixels tall, stands as near-black
silhouettes (#443c3a) from the left edge to four-fifths of the way across. The calm sea mirrors the sky: peach-tan
(#e3ba8c) on the left, grey-beige (#b0a78a) on the right, darkening downward. It moves in a long, low,
glassy swell 1.1 m high and 71 m from crest to crest, coming from the front left, with no wind ripples.
To the eye this is a dim, almost colourless night with the bright Earth and its path. No stars show
against this bright sky. No foam, no boats, no people.

**Frame data.**
- **View:** 94° (east), pitched 4° up; horizon at row 648. The frame is 58% sky, 3% land and 40% sea.
- **Exposure:** white is 0.237 cd/m².
- **Earth:**
  - at (970, 100), 20.3° up at 94.4°;
  - 53.6 px across, 58% lit, lit side toward the lower left;
  - 729 cd/m², amber #ffcd2d.
- **Glitter:** from the horizon to the bottom, across x = 845 to 1105. Its peak at the bottom is
  2.4 cd/m² (10 times white), golden yellow; it covers 9% of the sea.
- **Sky:**
  - top #dae0c6 (0.17 cd/m²);
  - 8° up #eed8ab;
  - at the horizon #e2c197 (0.13);
  - the left edge low down #ffc596 (0.15) and the right edge #c9b796 (0.11).
- **Islets:**
  - 18 px at the left edge, 43 px near x = 480, about 30 px to x = 1440, 15 px at x = 1560, sea beyond;
  - 1.75 to 2.5 km away; colour #443c3a.
- **Sea beside the glitter:** #e3ba8c on the left and #b0a78a on the right (0.09 to 0.13 cd/m²); bottom
  corners #8d7e62 and #727260.
- **Waves:**
  - swell 1.1 m, period 16.6 s, 71 m long, travelling toward 219°, 125° from the view;
  - wind 1.6 m/s, too light for short waves;
  - mean square slope 0.0046 (rms tilt 3.9°).
- **Stars:** none. The naked-eye limit here is about magnitude 0.7; the brightest stars in the frame
  appear at 4.6 to 5.5 through the air.
- **Tide:** +1.35 m.

## The eastern Smythii headland

The camera floats about 1 km off the west face of a narrow headland in the eastern Smythii basin
(93.22° E, 2.36° N), over water about 270 m deep. Its waves are those of the shore history's west-face
station, 1 km inshore in 22 m of water. It looks west across the open sea; the headland is behind the camera. No
land enters these frames. The Earth stands low in the west, rising and setting with the libration: below
the horizon for the first five scenes, just risen at the darkest hour.

### SH-1 · The Sun at its highest

*Day 4.6, seven and a half days before sunset. Sun 86.8° up, almost overhead. 87,200 lux; day vision.*

**Image prompt.** Photograph from a small boat on the open sea at midday, eye 2 m above the water, 28 mm
lens, looking west across empty water to a perfectly straight horizon, camera tilted 3 degrees down,
horizon 43% from the top. The Sun is almost directly overhead, out of frame, its glitter beneath the boat
out of view. The sky is a very pale, luminous cyan-blue (#c6eef6 at the top), almost white, fading to a pale
grey-green-cream at the horizon (#cfe0d6). It is evenly bright across the frame: a clean, washed-out
noon sky. The sea is a cool grey-teal (#8ca5a9 just below the horizon), about half as bright as the sky.
It darkens steadily to deep slate-teal (#4e5f5f) at the bottom of the frame, with a faint olive cast in
the nearest water. A moderate swell 0.6 m high and 23 m from crest to crest rolls straight toward the
camera from the horizon. Its crests run parallel to the horizon and the camera looks into their faces.
Lively wind ripples under a 3 m/s breeze blowing toward the camera cover the surface in a fine, crisp
texture of pale and dark facets. Empty sea, unbroken water, no foam, no land, no boats, no people.

**Frame data.**
- **View:** 270° (west), pitched 3° down; horizon at row 463. The frame is 43% sky and 57% sea.
- **Exposure:** white is 13,100 cd/m².
- **Sky:**
  - top #c6eef6 (10,400 cd/m²);
  - 8° up #ceecea;
  - at the horizon #cfe0d6 (9,360), uniform across.
- **Sea:**
  - below the horizon #8ca5a9 (4,610 cd/m²);
  - 3° down #859c9f;
  - 8° down #768b8c;
  - bottom #4e5f5f (1,380), with the corners #576869.
- **Waves:**
  - swell 0.6 m, period 9.4 s, 23 m long, travelling toward 94°, straight at the camera;
  - wind 3.1 m/s toward 97°;
  - mean square slope 0.025 (rms tilt 9.0°).
- **Sun and Earth:** the Sun 86.8° up. The Earth is 7.3° below the horizon straight ahead, so not
  visible.
- **Tide:** +0.21 m.

### SH-2 · The Sun 4° up

*Day 11.7, eight hours before sunset. Sun 4.0° up in the west. 6,800 lux; day vision.*

**Image prompt.** Photograph from a small boat on a glassy sea in the late afternoon, eye 2 m above the
water, 28 mm lens, looking west into the low Sun over an empty, perfectly flat horizon, horizon 53% from
the top. The Sun is a small, dazzling orange disk, 16 pixels across (#ff9d00 at its own exposure, blown
out here), 4 degrees above the horizon at the centre of the frame (x = 980, y = 461). Below it lies an
almost mirror-calm sea. The Sun's reflection is a narrow, intensely brilliant column of molten orange,
only about 50 pixels wide, running from the horizon straight down to the bottom of the frame. It is
broken into long, smooth, stretched highlights by a low swell. Its core is blinding white-orange and its
edges sharp. The sky is a warm amber-cream at the top (#ffe59d), deepening to #ffd584 at 8 degrees and
to soft orange at the horizon (#f4b475). It is brightest above the Sun, with no halo. Beside the
glitter column the calm water mirrors the sky almost exactly, slightly darker: warm apricot (#daaa6e)
just below the horizon, darkening to bronze-olive (#7d6b46) at the bottom. A gentle, glassy swell
0.9 m high and 48 m from crest to crest rolls toward the camera from straight ahead. Its smooth
undulations bend the reflections into long wavy bands. There are no wind ripples. Empty sea, no land, no
foam, no boats, no people.

**Frame data.**
- **View:** 268°, pitched 1° up; horizon at row 569. The frame is 53% sky and 47% sea.
- **Exposure:** white is 4,350 cd/m².
- **Sun:** at (980, 461), 4.02° up at 268.8°; 40 million cd/m², deep orange.
- **Glitter:** a narrow column across x = 975 to 1025, from the horizon to the bottom. Its peak is
  1.2 million cd/m² (290 times white); it covers 2% of the sea.
- **Sky:**
  - top #ffe59d (3,490 cd/m²);
  - 8° up #ffd584;
  - 3° up #ffc17a;
  - at the horizon #f4b475 (2,310);
  - the edges about 9% darker.
- **Sea beside the glitter:**
  - below the horizon #daaa6e (1,940 cd/m², 92% of the sky above);
  - 8° down #ab8a58;
  - bottom #7d6b46.
- **Waves:**
  - swell 0.9 m, period 13.7 s, 48 m long, travelling toward 75°, coming at the camera;
  - wind 1.2 m/s, too light for short waves;
  - mean square slope 0.0027 (rms tilt 3.0°).
- **Earth:** 8.5° below the horizon straight ahead, a new Earth.
- **Tide:** +0.31 m.

### SH-3 · The Sun 4° below the horizon

*Day 12.3, eight hours after sunset. 2,800 lux; day vision.*

**Image prompt.** Photograph from a small boat on a calm sea in the early evening twilight, eye 2 m
above the water, 28 mm lens, looking west to the place of sunset over a flat, empty horizon, horizon
57.5% from the top. The Sun is 4 degrees below the horizon at the centre of the frame. The sky glows
amber and orange:
- soft amber at the top (#ffe194), brightest above the centre;
- light orange (#ffce78) at 8 degrees;
- coral-orange just above the horizon (#f1a86b).

The sea is nearly calm and mirrors this sky as a slightly darker, softer copy: apricot-orange
(#d9a260) just below the horizon, about 90% as bright as the sky above, fading to warm bronze
(#83683f) at the bottom. The horizon line is barely visible as a slight darkening between two glowing
halves of the frame. A low glassy swell 0.8 m high and 58 m between crests rolls toward the camera.
Faint patches of the smallest ripples under a 1.9 m/s breeze dull the mirror slightly, and its long
smooth undulations stretch the reflected glow into gentle horizontal bands. Empty sea, no land, no foam,
no boats, no people.

**Frame data.**
- **View:** 268°, pitched 3° up; horizon at row 621. The frame is 57% sky and 43% sea.
- **Exposure:** white is 2,100 cd/m².
- **Sky:**
  - top #ffe194 (1,670 cd/m²);
  - 15° up #ffde87;
  - 8° up #ffce78;
  - 3° up #ffb76f;
  - at the horizon #f1a86b (997);
  - the edges 10% darker.
- **Sea:**
  - below the horizon #d9a260 (873 cd/m²);
  - 3° down #cc985a;
  - 8° down #af844f;
  - bottom #83683f (316).
- **Waves:**
  - swell 0.8 m, period 15.1 s, 58 m long, coming at the camera;
  - wind 1.9 m/s;
  - mean square slope 0.0028 (rms tilt 3.0°).
- **Sun and Earth:** the Sun is 4.1° below the horizon at 269°. The Earth is 7.6° below and new.
- **Tide:** +0.25 m.

### SH-4 · The Sun 15° below the horizon

*Day 13.3, 30 hours after sunset. 450 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the long evening, eye 2 m above the water,
28 mm lens, looking west over an empty horizon, horizon 57.5% from the top. A deep amber-orange twilight
fills the western sky:
- warm amber at the top (#ffd980);
- golden orange (#ffcb60) at 8 degrees;
- a rich coral-orange band along the horizon (#ff9e55), brightest at the centre and duller toward
  both edges.

The sea, gently rippled by a 2.8 m/s breeze, reflects it at about 80% of its brightness: amber-orange
(#d6964c) below the horizon, softening into a smooth, slightly grainy sheen, then darkening to brown-gold
(#856032) at the bottom. A swell 1.1 m high and 69 m from crest to crest comes from the front left; its
crests run diagonally. Empty sea, no land, no foam, no boats, no people.

**Frame data.**
- **View:** 268°, pitched 3° up; horizon at row 621. The frame is 57% sky and 43% sea.
- **Exposure:** white is 513 cd/m².
- **Sky:**
  - top #ffd980 (390 cd/m²);
  - 15° up #ffdb72;
  - 8° up #ffcb60;
  - 3° up #ffb058;
  - at the horizon #ff9e55 (239);
  - the edges 15% darker.
- **Sea:**
  - below the horizon #d6964c (188 cd/m²);
  - 3° down #c88c47;
  - 8° down #ad7a3e;
  - bottom #856032 (70).
- **Waves:**
  - swell 1.1 m, period 16.6 s, 69 m long, travelling toward 44°, 136° from the view;
  - wind 2.8 m/s;
  - mean square slope 0.027.
- **Sun and Earth:** the Sun is 15.2° below the horizon. The Earth is 6.1° below and 1% lit.
- **Tide:** +0.16 m.

### SH-5 · The Sun 35° below the horizon

*Day 14.9, nearly three days after sunset. 9.1 lux; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat at sea late in the long evening, eye 2 m
above the water, 28 mm lens, looking west over an empty horizon, horizon 57.5% from the top. A deep
twilight arch spans the western sky:
- glowing amber-gold about 8 degrees up at the centre (#ffcc3d);
- vivid orange at the horizon (#ff9d2c);
- duller ochre at the top of the frame (#eebb68);
- redder and darker toward both edges (#fe8739 low on the left).

It is the light of the very high air, which still sees the Sun 35 degrees below the horizon. The
sea, roughened by a 3.3 m/s breeze and a confused low swell, reflects the arch at half its brightness: orange (#c98436) just below
the horizon, deepening to dark brown-orange (#805423) at the bottom. Its surface is a soft-grained
texture of orange and brown facets. The Earth, a thin crescent, is 3 degrees below the horizon and not
visible. Muted to the eye; the photograph records the colour. Empty sea, no land, no stars, no foam,
no boats, no people.

**Frame data.**
- **View:** 266°, pitched 3° up; horizon at row 621. The frame is 57% sky and 43% sea.
- **Exposure:** white is 17.0 cd/m².
- **Sky:**
  - top #eebb68 (9.3 cd/m²);
  - 15° up #ffcc5a;
  - 8° up #ffcc3d (13.9);
  - 3° up #ffb42b;
  - at the horizon #ff9d2c (10.1);
  - the left edge low down #fe8739 (6.6).
- **Sea:**
  - below the horizon #c98436 (5.0 cd/m²);
  - 8° down #a86b2b;
  - bottom #805423 (1.9).
- **Waves:**
  - a confused swell 0.5 m high, period 15.1 s, 58 m long, spread 60°;
  - wind 3.3 m/s toward 245°;
  - mean square slope 0.038 (rms tilt 11.0°).
- **Earth:** 3.2° below the horizon at 263.6°, 8% lit.
- **Tide:** −0.03 m.

### SH-6 · The darkest hour: the Earth rising

*Day 18.1, six days after sunset; the Sun 74° below the horizon. 0.085 lux; night vision.*

**Image prompt.** Long-exposure night photograph from a small boat on a mirror-calm sea, eye 2 m above
the water, 28 mm lens, looking west, horizon 55% from the top. The Earth has just risen. It sits 1.5
degrees above the horizon at the centre of the frame (centre x = 958, y = 552), a disk 58 pixels across,
38% lit. Its lit part is a thick crescent along the lower edge of the disk, like a bowl, deep orange
(#ff8f00) from the long path through the air, overexposed with a bloom in this exposure. The rest of the
disk is invisible against the sky. Its reflection runs from the horizon straight down toward the camera
as a narrow, vivid orange column about 60 pixels wide, broken into a few long smooth bands by a very low
swell. It is an almost perfect mirror image, reaching halfway from the horizon to the bottom of the frame. The
sky is a faint warm glow:
- pale apricot-cream overhead (#fce4ae);
- golden-orange at 8 degrees (#ffd592);
- soft salmon-orange at the horizon (#fcb381), lit by the low Earth.

The bright star Aldebaran shows as a small point near the top of the frame, just left of centre
(x = 906, y = 84). Seven fainter stars are scattered across the sky:
- (899, 136);
- (623, 37);
- (1375, 250);
- (1612, 167);
- (1891, 180);
- (1713, 428);
- one more, fainter still.

The water away from the Earth's reflection is glassy and mirrors the sky: salmon-tan (#e2aa7c) at the
horizon, darkening to dark olive-brown (#816f50) at the bottom. A very low swell 0.3 m high and 48 m
between crests barely moves it, without ripples. To the eye this scene is a dark, almost colourless
night with the orange Earth and its path. Empty sea, no land, no foam, no boats, no people.

**Frame data.**
- **View:** 264.5°, pitched 2° up; horizon at row 595. The frame is 55% sky and 45% sea.
- **Exposure:** white is 0.059 cd/m².
- **Earth:**
  - at (958, 552), 1.55° up at 264.4°;
  - 58 px across, 38% lit, lit side toward the bottom;
  - 59 cd/m² (1,000 times white), deep orange #ff8f00; above the air it would be bluish white.
- **Glitter:** across x = 935 to 995, from the horizon to y = 835. Its peak is 10.4 cd/m²
  (180 times white), vivid orange.
- **Sky:**
  - top #fce4ae (0.047 cd/m²);
  - 15° up #ffe3a3;
  - 8° up #ffd592;
  - 3° up #ffc187;
  - at the horizon #fcb381 (0.032);
  - the edges 8 to 10% darker.
- **Sea beside the glitter:**
  - below the horizon #e2aa7c (0.027 cd/m²);
  - 8° down #ad8b62;
  - bottom #816f50 (0.0097).
- **Waves:**
  - swell 0.3 m, period 13.7 s, 48 m long, coming at the camera;
  - wind 0.6 m/s, glassy;
  - mean square slope 0.0014.
- **Stars recorded:**
  - Aldebaran (magnitude 0.85, seen at 3.0 through the air) at (906, 84);
  - six others seen at magnitudes 5.2 to 5.9: (623, 37), (899, 136), (1612, 167), (1891, 180),
    (1375, 250) and (1713, 428);
  - one fainter.
  - The naked-eye limit here is about 2.1, so the eye sees none.
- **Tide:** −0.30 m.

## Southern Mare Nubium, by Pitatus

The coast with the largest tide on the exposed tenth of the nearside's shores: this month's range is
7.7 m. The wave grid's shore point falls on the west shore of an island 4 by 6 km, rising to about
290 m. The camera floats 1.5 km off that shore (11.93° W, 28.11° S) in 65 m of water, with open sea to the
west and north. The mainland coast, with hills and crater rims of 400 to 1,100 m, lies 12 to 30 km to the south
and southwest. The Earth stands high in the northeast, 53 to 65 degrees up, out of every frame.

### SN-1 · The Sun at its highest

*Day 13.2, seven and a half days before sunset. Sun 63° up in the north, behind the camera. 75,000 lux;
day vision.*

**Image prompt.** Photograph from a small boat on a choppy sea at midday, eye 2 m above the water,
28 mm lens, looking south with the high Sun behind the camera, tilted 3 degrees down, horizon 43% from
the top. On the left rises the sunlit island, 1.7 to 5 km away: a smooth-shouldered hill of bare, dark
grey-brown ground (#847465) rising to 120 pixels above the horizon at the left edge (280 m high). It
slopes down to the water about a fifth of the way across the frame. It is front-lit and only lightly hazed.
It is high tide, 4 m above the mean level, and the sea laps high on the island's lowest slopes. Beyond
it, across the rest of the horizon, runs a faint low line of mainland hills 17 to 30 km away, 38 to 65
pixels tall, pale hazy blue-grey. The sky is very pale cyan-blue, almost white (#c6ecef at the top),
fading to pale grey-green-cream at the horizon (#d8e3d1). The sea is the roughest of the 24 scenes:
- a 6 m/s wind blows from behind the camera, over the right shoulder, raising steep, crowded wind waves
  and ripples over a 1.2 m swell 33 m from crest to crest that runs away to the left;
- the surface is a busy, high-contrast mosaic of pale sky-reflecting facets and dark slate troughs;
- grey-teal (#7b9295) below the horizon, darkening to deep slate (#4b5a58) at the bottom;
- faint olive in the nearest water.

Unbroken water, no foam. No boats, no people, no plants.

**Frame data.**
- **View:** 180° (south), pitched 3° down; horizon at row 463. The frame is 38% sky, 5% land and 57% sea.
- **Exposure:** white is 12,400 cd/m².
- **Sky:**
  - top #c6ecef (9,700 cd/m²);
  - 8° up #d4ede4;
  - at the horizon #d8e3d1 (9,250).
- **Island:**
  - 120 px at the left edge (284 m at 4.2 km), 104 px at x = 120, 82 px at x = 240, 54 px at x = 360;
  - the shore 1.7 to 2.1 km away.
  - Mainland beyond: 38 to 65 px (550 to 1,080 m at 17 to 30 km), the shore 11 to 21 km away.
  - Sunlit; haze transmittance 0.95 (green) and 0.89 (blue) on the island; island colour #847465.
- **Sea:**
  - below the horizon #7b9295 (3,330 cd/m²);
  - 8° down #6b7e7d;
  - bottom #4b5a58 (1,180).
- **Waves:**
  - swell 1.2 m, period 11.3 s, 33 m long, travelling toward 138°, 42° left of straight away;
  - wind 6.05 m/s;
  - mean square slope 0.053 (rms tilt 12.9°), the steepest of the 24 scenes;
  - spectrum 18 h before the moment.
  - At 6 m/s the first whitecaps would form; the study does not model them.
- **Earth:** 65° up in the north-northeast, a new Earth, not visible.
- **Tide:** +4.08 m.

### SN-2 · The Sun 4° up

*Day 20.3, nine hours before sunset. Sun 4.0° up in the west. 6,900 lux; day vision.*

**Image prompt.** Photograph from a small boat on a mirror-calm open sea in the late afternoon, eye 2 m
above the water, 28 mm lens, looking west into the low Sun, horizon 53% from the top. The Sun is a
small, blinding orange disk 16 pixels wide, 4 degrees above the horizon at the centre (x = 980,
y = 460). Below it a narrow, brilliant column of reflected sunlight runs from the horizon nearly to the
bottom of the frame. It is about 60 pixels wide, molten white-orange at its core, and stretched into
long smooth bands by a low swell. The sky is warm amber-cream at the top (#ffe59d), deepening to soft
orange at the horizon (#f4b475), brightest above the Sun. On the far left, a thin, faint line of distant
coast 12 to 23 km away rises only 2 to 20 pixels above the horizon, a hazy dull red-brown (#8b695b);
the rest of the horizon is open sea. The calm water beside the glitter mirrors the warm sky: apricot
(#dbaa6f) at the horizon, darkening to bronze-olive (#7d6b46) at the bottom. A smooth swell 1.1 m high
and 71 m from crest to crest comes from the front right; its long glassy undulations bend the
reflections. There are no ripples: the wind is under 1 m/s. No foam, no boats, no people.

**Frame data.**
- **View:** 270° (west), pitched 1° up; horizon at row 569. The frame is 52% sky, 0.4% land and 47% sea.
- **Exposure:** white is 4,390 cd/m².
- **Sun:** at (980, 460), 4.03° up at 270.8°; 40 million cd/m².
- **Glitter:** across x = 945 to 1005, from the horizon to y = 1005. Its peak is 1.7 million cd/m²
  (390 times white).
- **Sky:** top #ffe59d (3,520 cd/m²), 8° up #ffd585, at the horizon #f4b475 (2,340), the edges 9% darker.
- **Coast:** across x = 0 to 600, 2 to 20 px tall, 47 to 390 m high at 12 to 23 km; the shore 11.7 to
  21 km away. Backlit; colour #8b695b.
- **Sea beside the glitter:** below the horizon #dbaa6f (1,980 cd/m²), 8° down #ab8a59, bottom #7d6b46.
- **Waves:**
  - swell 1.1 m, period 16.6 s, 71 m long, travelling toward 134°, coming from the front right;
  - wind 0.6 m/s;
  - mean square slope 0.0024.
- **Earth:** 59° up in the northeast behind the camera, 62% lit.
- **Tide:** +0.79 m.

### SN-3 · The Sun 4° below the horizon

*Day 21.1, nine hours after sunset. 2,870 lux; day vision.*

**Image prompt.** Photograph from a small boat on a mirror-calm sea in the early evening twilight, eye
2 m above the water, 28 mm lens, looking west to the place of sunset, horizon 57.5% from the top. The
sky is a luminous gradient:
- soft amber at the top (#ffe194);
- orange through the middle;
- coral (#f1a86b) at the horizon, brightest above the centre where the Sun has gone down.

The sea is almost perfectly calm and mirrors the sky at nearly its full brightness. Just below the
horizon it is apricot-orange (#eba865), 98% as bright as the sky just above, so the horizon line all but
disappears. It darkens to warm bronze (#876c40) at the bottom of the frame. A long, low glassy swell
1.1 m high and 71 m between crests rolls in from the front right, drawing slow, broad bends into the
mirrored glow. On the far left a hair-thin, hazy line of distant coast rises 1 to 18 pixels, dull
red-brown (#8a6153). No ripples, no foam, no boats, no people.

**Frame data.**
- **View:** 266°, pitched 3° up; horizon at row 621. The frame is 57% sky, 0.3% land and 43% sea.
- **Exposure:** white is 2,140 cd/m².
- **Sky:**
  - top #ffe194 (1,700 cd/m²);
  - 8° up #ffce79;
  - 3° up #ffb76f;
  - at the horizon #f1a86b (1,020);
  - the edges 10% darker.
- **Sea:**
  - below the horizon #eba865 (997 cd/m²);
  - 3° down #db9e5e;
  - 8° down #bb8b51;
  - bottom #876c40 (346).
- **Coast:** across x = 0 to 720, 1 to 18 px (60 to 380 m at 12 to 23 km); colour #8a6153.
- **Waves:**
  - swell 1.1 m, period 16.6 s, 71 m long, coming from the front right;
  - wind 1.4 m/s, below the threshold for short waves;
  - mean square slope 0.0024.
- **Sun and Earth:** the Sun is 4.0° below the horizon at 266.4°. The Earth is 58° up behind the camera,
  70% lit.
- **Tide:** +0.18 m.

### SN-4 · The Sun 15° below the horizon

*Day 22.1, a day and a half after sunset. 466 lux; day vision.*

**Image prompt.** Photograph from a small boat on a wind-roughened sea in the long evening, eye 2 m
above the water, 28 mm lens, looking west, horizon 57.5% from the top. The western sky is a broad
twilight of warm amber (#ffd980 at the top) deepening to a coral-orange band along the horizon
(#ff9e56), brightest at the centre. A thin, hazy line of distant coast 12 to 23 km away runs along the
left half of the horizon, 5 to 21 pixels tall, dull red-brown (#945a44). The sea is choppy under a
5 m/s wind. It reflects the sky at about 60% of its brightness as a busy texture of orange highlights
and brown troughs: bronze-orange (#b48447) below the horizon, dark brown-gold (#78572f) at the bottom.
A swell 1.6 m high and 49 m from crest to crest runs away from the camera, slightly to the right; short
steep wind waves ride on it. Unbroken water, no foam, no boats, no people.

**Frame data.**
- **View:** 260°, pitched 3° up; horizon at row 621. The frame is 57% sky, 0.4% land and 43% sea.
- **Exposure:** white is 528 cd/m².
- **Sky:**
  - top #ffd980 (402 cd/m²);
  - 15° up #ffdb72;
  - 8° up #ffcb61;
  - 3° up #ffb058;
  - at the horizon #ff9e56 (246);
  - the edges 15% darker.
- **Coast:** across x = 0 to 840, 5 to 21 px (70 to 400 m at 12 to 23 km); colour #945a44.
- **Sea:**
  - below the horizon #b48447 (140 cd/m²);
  - 8° down #966c3a;
  - bottom #78572f (58).
- **Waves:**
  - swell 1.6 m, period 13.7 s, 49 m long, travelling toward 286°, 26° right of straight away;
  - wind 5.1 m/s;
  - mean square slope 0.046 (rms tilt 12.2°).
  - At this wind the first whitecaps would form; the study does not model them.
- **Tide:** −0.59 m.

### SN-5 · The Sun 35° below the horizon

*Day 24.1, three and a third days after sunset. 11.7 lux; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat at sea late in the long evening, eye 2 m
above the water, 28 mm lens, looking west-southwest, horizon 57.5% from the top. A deep golden-orange
twilight arch fills the western sky:
- brightest about 8 degrees up at the centre (#ffce4e);
- vivid orange along the horizon (#ffa03f);
- duller ochre at the top (#ecbd73);
- redder toward the edges.

On the left, coastal hills 18 km away rise as a dark red-brown silhouette (#8f4c38), 63 pixels tall at the left edge, falling
toward the centre into a thin, hazy line of distant coast that ends two-thirds of the way across. The
sea, rippled by a 4 m/s breeze over a 1.1 m swell running away to the right, reflects the arch at half
its brightness: orange (#c9863f) below the horizon, dark brown-orange (#805529) at the bottom. Muted to
the eye; the photograph records the colour. No stars, no foam, no boats, no people.

**Frame data.**
- **View:** 246°, pitched 3° up; horizon at row 621. The frame is 56% sky, 1% land and 43% sea.
- **Exposure:** white is 16.8 cd/m².
- **Sky:**
  - top #ecbd73 (9.3 cd/m²);
  - 15° up #ffcd66;
  - 8° up #ffce4e (13.8);
  - 3° up #ffb640;
  - at the horizon #ffa03f (10.2);
  - the edges low down #ff9047 (7.3).
- **Coast:**
  - 63 px at the left edge (719 m at 18.2 km), 42 px at x = 120, 27 px at x = 240;
  - 7 to 16 px to x = 720, 20 px at x = 960;
  - sea beyond x = 1260.
  - Colour #8f4c38.
- **Sea:**
  - below the horizon #c9863f (5.0 cd/m²);
  - 8° down #a86d33;
  - bottom #805529 (1.9).
- **Waves:**
  - swell 1.1 m, period 11.3 s, 33 m long, travelling toward 302°, 56° right of the view;
  - wind 4.2 m/s;
  - mean square slope 0.037.
- **Earth:** 55° up behind the camera, 93% lit.
- **Tide:** −1.79 m.

### SN-6 · The darkest hour

*Day 0.5, nine days after sunset; the Sun 52.8° below the horizon, straight ahead. 2.9 lux, mostly
from a high, nearly full Earth behind the camera; between day and night vision.*

**Image prompt.** Long-exposure night photograph from a small boat, eye 2 m above the water, 28 mm lens,
looking southeast toward the island, horizon 57.5% from the top. The island fills the lower part of the
view from edge to edge, 1.5 to 5 km away. It is a long, rounded ridge of bare ground rising 51 to 110
pixels above the horizon, highest at the centre (250 to 290 m). Beyond its right end, far hazy hills of
the mainland 23 km away stand 51 to 58 pixels tall. The island's slopes are lit softly from behind the
camera by a high, nearly full Earth (out of frame), so they show as a dim muted red-brown (#6c544f), their
shapes visible. It is low tide, 3 m below the mean level: a band of wet, darker exposed shore
2 to 4 pixels high runs along the island's foot at the waterline. The sky is a pale, nearly neutral
light grey at the top (#d6d5d3), warming downward into a soft pinkish peach glow toward the horizon
(#ffdebd at 8 degrees, #ffc6ab at the horizon) above the hidden Sun. The sea is lightly rippled by a
4 m/s breeze over a short swell 0.9 m high and 23 m from crest to crest, coming toward the camera. It
reflects the glow as a dim greyish rose-taupe (#a18d83) below the horizon, darkening to dark grey-brown
(#6d6059) at the bottom. No stars show in this bright earthlit sky. No foam, no boats, no people.

**Frame data.**
- **View:** 137° (southeast), pitched 3° up; horizon at row 621. The frame is 50% sky, 8% land and
  43% sea.
- **Exposure:** white is 0.97 cd/m².
- **Sky:**
  - top #d6d5d3 (0.64 cd/m²);
  - 15° up #f2ddcb;
  - 8° up #ffdebd (0.78);
  - 3° up #ffd2b2;
  - at the horizon #ffc6ab (0.68) at the centre and #f7bfac at the edges (0.58).
- **Island:**
  - 66 px at the left edge, rising through 73, 78, 81 and 88 px to 110 px at x = 840 to 1080, 97 px at
    x = 1320, 54 px at x = 1560;
  - 2.8 to 4.9 km away, 105 to 288 m high; the shore 1.45 to 2.1 km away.
  - Mainland beyond: 51 to 58 px (850 to 925 m at 23 to 25 km).
  - Earthlit; colour #6c544f.
- **Sea:**
  - below the horizon #a18d83 (0.27 cd/m²);
  - 3° down #9a857b;
  - 8° down #89766d;
  - bottom #6d6059 (0.12).
- **Waves:**
  - swell 0.9 m, period 9.3 s, 23 m long, travelling toward 303°, coming at the camera;
  - wind 4.1 m/s;
  - mean square slope 0.035.
- **Earth:** 53.3° up at 22° behind the camera, 90% lit. The sky overhead and behind is the blue
  earthlit sky of the regime results (about 0.5 cd/m², chromaticity 0.27, 0.30).
- **Stars:** none. The naked-eye limit here is about magnitude −0.7.
- **Tide:** −3.06 m. This month's range is 7.7 m, from +4.1 m at noon.

## The South Pole–Aitken coast, by Mare Ingenii

The far side, where the Earth never rises. The camera floats 5.4 km off a coast that runs east to west
(162.875° E, 37.125° S), over water 1.4 km deep. Along the coast to the south stands a range of bare
highland mountains 700 to 1,330 m high, 10 to 21 km away. Behind it, 50 to 72 km to the southwest, rises
a higher massif of 2,300 to 3,500 m. The open sea lies to the north.
No wave run reaches this sea: its waves borrow a nearside spectrum of median steepness, turned to the
local wind. The land here is the brighter, tan highland soil.

### IC-1 · The Sun at its highest

*Day 28.4, seven and a half days before sunset. Sun 54° up in the north, behind the camera. 66,000 lux;
day vision.*

**Image prompt.** Photograph from a small boat at sea at midday, eye 2 m above the water, 28 mm lens,
looking south across 5 km of open water to a mountainous coast, with the high Sun behind the camera,
tilted 2 degrees down, horizon 45% from the top. A range of bare, sunlit mountains spans the whole
width of the frame, 10 to 21 km away. It rises 77 pixels above the horizon at the left edge, climbing
in rounded peaks and shoulders to 174 pixels at two-thirds of the way across (1,270 m), then falls to
115 pixels at the right edge. The mountains are pale tan-beige highland ground (#b9a48d), softened and
cooled by haze. Their front faces are lit, with gentle shading in the folds and no vegetation. The sky
is a very pale cyan-blue, almost white (#c8eae9 at the top), fading to pale cream-grey at the horizon
behind the peaks (#dce2ca). The sea is a grey-teal slate (#8ea29d below the horizon), darkening to deep
slate-teal (#505e5b) at the bottom. The band of water nearest the horizon mirrors the pale mountains in
blurred streaks. A swell 1.2 m high and 33 m from crest to crest comes from the front left; wind ripples
under a 3.7 m/s breeze give it a lively, fine texture. No foam, no boats, no people, no plants, no
buildings.

**Frame data.**
- **View:** 180° (south), pitched 2° down; horizon at row 490. The frame is 32% sky, 13% land and 55% sea.
- **Exposure:** white is 12,200 cd/m².
- **Sky:**
  - top #c8eae9 (9,400 cd/m²);
  - 8° up #d8eddc;
  - at the horizon #dce2ca (9,020).
- **Mountains:**
  - 77 px at the left edge, 97 px at x = 120, 112 to 132 px to x = 600;
  - 140 to 153 px to x = 1080, 167 to 174 px at x = 1200 to 1320, 163 to 168 px to x = 1680;
  - 139 px at x = 1800, 115 px at the right edge.
  - 10.5 to 21 km away, 715 to 1,325 m high; the shore 5.5 to 7.1 km away.
  - Sunlit; haze transmittance 0.86 (green) and 0.72 (blue); colour #b9a48d.
- **Sea:**
  - below the horizon #8ea29d (4,140 cd/m²);
  - 3° down #879994;
  - 8° down #788883;
  - bottom #505e5b (1,275).
- **Waves:**
  - swell 1.2 m, period 11.3 s, 33 m long, borrowed, travelling toward 325°, 145° from the view;
  - wind 3.7 m/s;
  - mean square slope 0.028.
- **Earth:** 55° below the horizon.
- **Tide:** −0.53 m.

### IC-2 · The Sun 4° up

*Day 5.9, ten hours before sunset. Sun 3.9° up in the west. 6,900 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the late afternoon, eye 2 m above the water,
28 mm lens, looking west along a mountainous coast into the low Sun, horizon 53% from the top. The Sun
is a small, blinding orange disk 16 pixels wide, 3.9 degrees up just right of centre (x = 1061,
y = 463). Below it a broad, sparkling orange glitter path, 100 to 260 pixels wide, runs from the horizon
to the bottom of the frame over a rippled sea. On the left the coast runs out toward the horizon:
- a near coastal ridge 20 km away stands as a dark, hazy red-brown silhouette, 99 pixels tall at the
  left edge;
- beyond it, from an eighth to half of the way across, a far, much higher mountain massif 50 to 72 km
  away rises 22 to 63 pixels above the horizon;
- the massif is so hazed that it is only a faint, pale, warm silhouette barely darker than the sky
  behind it.

The right half of the horizon is open sea. The sky is warm amber-cream at the top (#ffe59d), deepening
to soft orange at the horizon (#f4b475). The sea beside the glitter is warm tan-gold (#bd9c64), darkening
to dark olive-bronze (#726242) in the lower corners. A swell 1.2 m high and 33 m between crests comes
from the front left; a 3.3 m/s breeze ripples it. No foam, no boats, no people.

**Frame data.**
- **View:** 268°, pitched 1° up; horizon at row 569. The frame is 50% sky, 3% land and 47% sea.
- **Exposure:** white is 4,380 cd/m².
- **Sun:** at (1061, 463), 3.93° up at 271.8°; 39 million cd/m².
- **Glitter:** across x = 975 to 1235. Its peak at the horizon is 124,000 cd/m² (28 times white).
- **Sky:** top #ffe59d (3,510 cd/m²), 8° up #ffd585, at the horizon #f4b475 (2,330).
- **Coast:**
  - near ridge: 99 px at the left edge and 71 px at x = 120 (930 to 1,230 m at 20 km);
  - far massif: 55 to 63 px at x = 240 to 600, 35 px at x = 720, 22 to 25 px at x = 840 to 960
    (2,300 to 3,440 m at 50 to 72 km);
  - sea beyond x = 1060.
  - Backlit. The near ridge has haze transmittance 0.76 (green) and 0.54 (blue), colour #8b6d5d. The
    far massif has about 0.42 and 0.15, nearly the horizon sky's colour.
- **Sea beside the glitter:**
  - below the horizon #bd9c64 (1,550 cd/m²);
  - 8° down #907951 on the left;
  - bottom corners #726242.
- **Waves:**
  - swell 1.2 m, period 11.3 s, 33 m long, borrowed, travelling toward 43°, 135° from the view;
  - wind 3.3 m/s;
  - mean square slope 0.026.
- **Tide:** +0.56 m.

### IC-3 · The Sun 4° below the horizon

*Day 6.8, ten hours after sunset. 2,850 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the early evening twilight, eye 2 m above the
water, 28 mm lens, looking west along a mountainous coast, horizon 57.5% from the top. The sky glows:
- soft amber at the top (#ffe194);
- orange through the middle;
- coral (#f1a86b) along the horizon, brightest just right of centre where the Sun has set.

On the left the coast's near ridge, 15 to 21 km away, stands as a dark, hazy red-brown silhouette
(#7d5d4e), 60 to 105 pixels tall. Beyond it, across the middle of the frame, the far high massif 50 to
72 km away shows as a pale, ghostly warm silhouette 25 to 58 pixels tall, little darker than the sky.
The right third of the horizon is open sea. The sea, gently rippled by a 2 m/s breeze over a 1.2 m swell,
reflects the glow at about 80% of its brightness: apricot-orange (#cf9d5e) below the horizon, warm
bronze (#80653e) at the bottom. The water just below the horizon on the left mirrors the dark ridge in
soft, broken streaks. No foam, no boats, no people.

**Frame data.**
- **View:** 262°, pitched 3° up; horizon at row 621. The frame is 54% sky, 3% land and 43% sea.
- **Exposure:** white is 2,120 cd/m².
- **Sky:**
  - top #ffe194 (1,690 cd/m²);
  - 8° up #ffce78;
  - at the horizon #f1a86b (1,010);
  - the left edge 12% darker.
- **Coast:**
  - near ridge: 105 px at the left edge, 97 px at x = 120, 87 px at x = 240, 60 px at x = 360
    (580 to 1,250 m at 15 to 21 km);
  - far massif: 57 to 59 px at x = 480 to 720, 41 px at x = 840, 25 to 29 px to x = 1080, 7 px at
    x = 1200 (2,400 to 3,500 m at 49 to 72 km);
  - sea beyond x = 1230.
  - Backlit; colour #7d5d4e.
- **Sea:**
  - below the horizon #cf9d5e (815 cd/m²);
  - 8° down #a7804d;
  - bottom #80653e (304).
- **Waves:**
  - swell 1.2 m, period 11.3 s, 33 m long, borrowed, travelling toward 19°;
  - wind 2.1 m/s;
  - mean square slope 0.019.
- **Tide:** +0.79 m.

### IC-4 · The Sun 15° below the horizon

*Day 7.9, a day and a half after sunset. 456 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the long evening, eye 2 m above the water,
28 mm lens, looking west-southwest along a mountainous coast, horizon 57.5% from the top. A broad amber
twilight fills the sky: warm amber at the top (#ffd980), deepening to a coral-orange band along the
horizon (#ff9e55). The near coastal mountains, 15 to 21 km away, stand on the left as a dark red-brown
silhouette (#7b513f), 77 to 114 pixels tall. They lead into the pale, hazed far massif across the
middle, 50 to 72 km away and 24 to 61 pixels tall, its outline soft against the glow. The right quarter
is open sea. The sea, rippled by a 3.6 m/s breeze over a 1.2 m swell crossing from left to right,
reflects the sky at three-quarters of its brightness: amber-orange (#d1944b) below the horizon, brown-gold
(#835f31) at the bottom, with the mountains' dark reflection blurred into the water below them on the
left. No foam, no boats, no people.

**Frame data.**
- **View:** 255°, pitched 3° up; horizon at row 621. The frame is 53% sky, 4% land and 43% sea.
- **Exposure:** white is 520 cd/m².
- **Sky:**
  - top #ffd980 (395 cd/m²);
  - 8° up #ffcb60;
  - at the horizon #ff9e55 (242);
  - the left edge 17% darker.
- **Coast:**
  - near mountains: 114 px at the left edge, 106, 98 and 92 px at x = 120 to 360, 77 px at x = 480
    (975 to 1,270 m at 15 to 21 km);
  - far massif: 51 to 61 px at x = 600 to 960, 24 to 33 px at x = 1080 to 1320 (2,500 to 3,400 m at
    52 to 72 km);
  - sea beyond x = 1420.
  - Backlit; colour #7b513f.
- **Sea:**
  - below the horizon #d1944b (183 cd/m²);
  - 8° down #aa783d;
  - bottom #835f31 (69).
- **Waves:**
  - swell 1.2 m, period 11.3 s, 33 m long, borrowed, travelling toward 357°, crossing from left to right;
  - wind 3.6 m/s;
  - mean square slope 0.027.
- **Tide:** +1.10 m.

### IC-5 · The Sun 35° below the horizon

*Day 10.2, nearly four days after sunset. 8.9 lux; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat at sea late in the long evening, eye 2 m
above the water, 28 mm lens, looking southwest along the mountainous coast, horizon 57.5% from the top.
A deep twilight arch glows over the mountains:
- brightest amber-gold about 8 degrees up at the centre (#ffcc3c);
- vivid orange along the horizon (#ff9d2b);
- ochre at the top (#eebb68);
- redder toward the edges.

The coastal mountains fill the whole width of the horizon as a dark red-brown silhouette (#7e432c):
- tallest at the left edge, 178 pixels (1,240 m at 12 km);
- stepping down through rounded shoulders to about 90 pixels across the centre;
- passing on the right into the far high massif, 51 to 72 km away and 27 to 61 pixels tall, paler
  with distance.

The sea is nearly calm: no wind ripples under a 1.7 m/s breeze, only a smooth 1.2 m swell coming toward
the camera. It mirrors the glow: orange (#f49834) below the horizon, with the mountains' dark shapes
reflected in the nearest band of water, darkening to brown-orange (#895b25) at the bottom. Muted to the
eye; the photograph records the colour. No stars, no foam, no boats, no people.

**Frame data.**
- **View:** 236°, pitched 3° up; horizon at row 621. The frame is 50% sky, 7% land and 43% sea.
- **Exposure:** white is 16.7 cd/m².
- **Sky:**
  - top #eebb68 (9.1 cd/m²);
  - 15° up #ffcc59;
  - 8° up #ffcc3c (13.7);
  - 3° up #ffb42b;
  - at the horizon #ff9d2b (10.0);
  - the edges low down #ff8b37 (7.1).
- **Coast:**
  - 178 px at the left edge, 156 px at x = 120, 121 px at x = 240;
  - 82 to 99 px from x = 360 to 960, 55 to 56 px at x = 1080 to 1200;
  - 58 to 61 px at x = 1320 to 1440 (the far massif at 61 to 65 km);
  - 27 to 43 px to the right edge.
  - Backlit; colour #7e432c.
- **Sea:**
  - below the horizon #f49834 (7.0 cd/m²);
  - 8° down #c3792b;
  - bottom #895b25 (2.2).
- **Waves:**
  - swell 1.2 m, period 11.3 s, 33 m long, borrowed, travelling toward 32°, coming at the camera;
  - wind 1.7 m/s, too light for short waves;
  - mean square slope 0.015.
- **Stars:** Achernar is 18.8° up at the top left. It appears at magnitude 2.6 through the air, below
  what the frame records against this sky.
- **Tide:** +1.49 m.

### IC-6 · The darkest hour: the far side's twilit night

*Day 13.7, seven and a third days after sunset; the Sun 51.8° below the southern horizon, straight
ahead. 0.46 lux, all from the twilight; between day and night vision.*

**Image prompt.** Long-exposure night photograph from a small boat at sea, eye 2 m above the water,
28 mm lens, looking south toward a mountainous coast, horizon 57.5% from the top. All night a golden
twilight arch stands over the southern mountains, lit by a Sun more than 50 degrees below the horizon:
- warm amber-gold, brightest about 8 degrees up (#ffcf52);
- orange near the mountain tops;
- ochre-orange at the top of the frame (#f0bf74).

The mountain range, 10 to 21 km away, spans the whole frame as a dark silhouette, 76 pixels tall at the
left edge, rising to 172 pixels two-thirds of the way across (1,270 m), 115 pixels at the right edge.
The haze tints it a deep rust-brown (#9a4f34) in this long exposure. It hides the arch's lowest, reddest
part. The sea is lightly rippled by a 2 m/s breeze over a 1.2 m swell running away to the right. It
reflects the arch as a warm orange sheen (#e29443 below the horizon at the centre), with the mountains'
dark outline mirrored, blurred, in the band of water just below the horizon. It darkens to brown-orange
(#865a2b) at the bottom. The Earth never rises here. No stars show against this glowing sky; the
brightest in view, Achernar, is too dim through the air. To the eye: a dim, nearly colourless twilight
with a pale glow over black mountains. No foam, no boats, no people.

**Frame data.**
- **View:** 180° (south), pitched 3° up; horizon at row 621. The frame is 45% sky, 13% land and 43% sea.
- **Exposure:** white is 0.62 cd/m².
- **Sky:**
  - top #f0bf74 (0.35 cd/m²);
  - 15° up #ffcc69;
  - 8° up #ffcf52 (0.51);
  - 3° up #ffbd41;
  - at the horizon, hidden behind the range, #ffa93e (0.43);
  - the edges 11 to 29% darker.
  - The sky just above the skyline is about #ffbd41 at the centre.
- **Mountains:** as at noon. 76 px at the left edge, 111 to 131 px to x = 600, 139 to 152 px to
  x = 1080, 165 to 172 px at x = 1200 to 1440, 161 to 164 px to x = 1680, 138 px at x = 1800 and 115 px
  at the right edge. Backlit by the arch; colour #9a4f34.
- **Sea:**
  - below the horizon #e29443 (0.24 cd/m²) for open sky; darker where it mirrors the range;
  - 8° down #b77736;
  - bottom #865a2b (0.08).
- **Waves:**
  - swell 1.2 m, period 11.3 s, 33 m long, borrowed, travelling toward 232°, 52° right of straight away;
  - wind 2.1 m/s;
  - mean square slope 0.019.
- **Earth:** 45° below the horizon.
- **Stars:** none recorded. Achernar is 8.3° up at (1119, 398), appearing at magnitude 3.75; the
  naked-eye limit here is about −0.4.
- **Tide:** +0.93 m.

## Boundaries

- **Sky and sea colours** come from the regime panoramas, at 1 by 0.5 degrees. They are statistical: no
  individual waves, open sea in every direction. Where land rises behind the sea, the scene text says
  how the nearest water mirrors it, without numbers.
- **The coast** is the LOLA altimetry grid, about 118 m between points, at the moment's tide. Its colour
  is a placeholder: bare Apollo soil, lit by the clear sky and the disks on a slope facing the camera,
  and hazed by extinction at the ground toward the colour of the horizon sky. Beaches, cliffs, surf, wet
  sand and plants are not modelled. The Smythii frames show no land: the headland is behind the camera.
- **Waves:** significant heights, periods and directions come from the wave spectra.
  - At the nearside shores these are the restart files up to 21 hours from the moment.
  - At Smythii they are the hourly west-face spectra.
  - At Ingenii they are a borrowed nearside spectrum.
  - Short waves follow the GCM's stress at the hour. Whitecaps and foam are not modelled.
- **Stars:** the eye's limit follows Schaefer (1990). A star counts as recorded when its point image,
  spread over about four pixels of this frame, reaches a fifth of the sky's brightness around it.
  Planets are not placed.
- **Approximate reading:** the Earth and the Sun are uniform disks, without the Earth's clouds and
  continents and without limb darkening. Refraction is left out, about 0.14 degrees at the horizon in
  this air.

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.sea_appearance.scenes
```

The runner reads the regime panoramas and the coast grids on the research drive
(`research/runs/sea_appearance/`); `python -m geography.coast_terrain --download` restores the grids'
LOLA rows.
