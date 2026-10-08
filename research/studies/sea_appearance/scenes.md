<!-- The words of scenes.md. Names in braces, {SCENE:name}, are filled from results/scenes.json by
     scenes.py, which writes scenes.md; edit the words here and the numbers in the product. -->
# The four coasts, scene by scene

Twenty-four views of the Open Moon's seas: four coasts at the six moments of the
[results by regime](README.md#the-seas-by-regime), from the Sun at its highest to the
darkest hour of the night. Each scene is the frame the
[sea-scene renderer](../../../visualization/reference-renderer/seas/) will draw. They are described here
in words and numbers, so that they can be pictured, or illustrated with an image
generator, before the full renderings are made.

These frames isolate the sea, sky and coastal lighting. **Land vegetation is omitted for this
optical study; the bare shorelines do not represent the intended mature landscape.** The coasts
have not been established as naturally barren. The water is modelled as a living, productive
coast, with phytoplankton and dissolved organic matter supplied as biosphere design guesses.
Coastal vegetation and habitats remain to be supplied by the domain research.

Every number is filled from [scenes.json](results/scenes.json), which
[scenes.py](scenes.py) writes from the study's products:
- **the colours and brightness of sky and sea:** the model of the results by regime (the solved
  sky, the rough-water reflection, both disks' glitter, the water's own light), evaluated for
  each frame's directions with that scene's own wave slopes, and the Sun and the Earth as its
  camera sees them;
- **the coast:** LOLA heights at the moment's tide;
- **the waves:** that hour's wave spectrum;
- **the stars:** the bright-star catalogue.

These are calculated scenes with stated gaps; an image made from them illustrates the
calculation. Two kinds of number are approximate:
- **Land colours:** a placeholder ground, lit by the clear sky and hazed by the air.
- **Sea colours:** they assume open water in every direction. Where land rises behind the sea,
  the water just below the horizon also mirrors the land, as each scene notes.

The supplemental [Nubium local-midnight view](README.md#nubium-at-local-midnight)
uses a taller frame to include a nearly full Earth high above the sea. Its
[human-vision display study](../../../visualization/sea-appearance/README.md#nubium-midnight-and-human-vision)
is separate from the photographic colour convention used by the 24 frames below.

## Reading the scenes

- **Generated illustrations** use an immersive, high-fidelity photographic treatment. Keep scene IDs,
  captions and modelling notes outside the image; the image itself has no text or labels. Use
  calculated renders or geometry and colour guides as inputs, refining components and composing them
  across passes when needed. Review framing, lighting, shadows, atmospheric depth and wave structure
  against the study inputs before selecting a final composite.
- **The frame** is a photograph 1920 by 1080 pixels through a rectilinear lens with a 65-degree
  horizontal field of view (a 28 mm lens on a full-frame camera; 39.4 degrees vertically), taken from
  a small boat with the eye 2 m above the water, the camera level and pitched as stated. Pixel
  positions count from the top left.
  - A pixel at the centre of the frame spans 0.038 degrees, and objects grow slightly toward the edges.
  - The sea's horizon lies 2.6 km away and 0.087 degrees below level; its row is given for the
    frame's centre.
- **Colours** are sRGB hex values at one exposure per frame, the luminance shown white being stated.
  - White is 1.25 times the brightest twentieth of the sky in the frame.
  - Brighter light rolls off smoothly toward white from 80% of it.
  - The display's white is Earth daylight (D65), without adapting to the light of the moment.
  - That is a photograph's rendering: in the night and late-evening scenes it is a long exposure,
    and the eye sees less colour than the values show. Luminances are in cd/m².
- **Directions** are compass bearings. Wave directions say where the waves travel toward, and how that
  lies against the view.

## What every scene shares

- **Clear air, seven times Earth's column.** The air is 1.2 atmospheres at the surface and very tall.
  The sky is cloudless in every scene. The calculated air holds no dust or haze droplets, so the Sun
  stands in a broad, smooth glow without a bright halo round it.
- **The Sun's colour.** Through this much air the Sun is pale gold even high in the sky. It is orange at
  15 degrees up and deep orange near the horizon. Its disk is 0.53 degrees across, about 14
  pixels.
- **The daytime sky** is a very pale cyan-blue, almost white, becoming a pale grey-green-cream toward the
  horizon. The sea under it is a grey-teal slate.
- **The evening lasts days.** The Sun sinks half a degree an hour.
  - With it 4 degrees below the horizon the western sky is amber and orange.
  - At 15 degrees below the western sky is amber and coral, and the eastern a soft golden cream.
  - At 35 degrees below, three to four days after sunset, a deep orange arch still glows over the
    western horizon.
- **The horizon is close and the air hazy.**
  - Distant hills keep two-thirds of their contrast at 30 km in green light and 40% in blue.
  - At 65 km they keep 40% and 15%, so far ranges fade toward the colour of the sky at the horizon.
  - The Moon's curvature hides the lowest metres of coasts beyond 2.6 km.
- **The Earth** is 1.80 to 2.02 degrees across as its distance changes through
  the month, 47 to 53 pixels near the centre of the frame, nearly four times the
  Sun's width. It hangs nearly still in the nearside sky,
  shows phases, and never rises on the far side.
  - Above the air it is bluish white; through the air it is amber to orange, deeper the lower it stands.
  - The calculation gives it no visible continents or clouds: it is a smooth lit shape. The unlit part
    of its disk looks like the surrounding sky.
- **The sea.** Deep water carries long swell from the wave runs, covered with wind ripples where the
  wind blows. In light winds the sea goes glassy.
  - At lunar gravity the waves in which gravity and surface tension balance are 4.2 cm long, against
    1.7 cm on Earth. These are the slowest of all water waves, and much of the short waves' slope lies
    near that length.
  - Capillary ripples continue below it to a few millimetres. The wind raises them harder than on
    Earth.
  - The water is the biosphere's productive coast: looking straight down into it, it is olive green.
    Seen at low angles the sea is mostly reflected sky.
  - The water is unbroken: no foam or whitecaps, which the study does not yet model.
- **The land is a placeholder:** bare ground with the colour of Apollo soils, dark grey-brown mare soil,
  and a brighter tan highland soil at the Ingenii coast. Vegetation is omitted to isolate coastal
  optics. Buildings, roads, boats and people are also omitted from these study frames.
- **Few stars.** Even the darkest hours are bright. The earthlit or twilit sky is 0.03 to 0.8 cd/m²,
  100 to 3,000 times Earth's moonless night, and the thick air dims stars near the horizon by two to four
  magnitudes. Only one frame records stars.

## Western Oceanus Procellarum, by Russell

The most exposed coast of the study. The camera floats at the wave runs' shore point (73.875° W,
28.875° N) in 78 m of water. Low islets 20 to 60 m high lie 1.8 to 2.5 km to the east-northeast. The
mainland shore runs 3 to 8 km away to the south and west, and behind it a belt of hills and crater rims
rises to 650 to 1,360 m at 15 to 42 km. Long ocean swell from the open sea reaches this coast. The Earth
stands in the east all month, 6 to 21 degrees up.

### WP-1 · The Sun at its highest

*Day 18.3 of the month (from 7 February 2038), seven and a half days before sunset. Sun 59.9°
up in the south. 72,500 lux on the ground; day vision.*

**Image prompt.** Photograph from a small boat on a wide sea at midday under a cloudless sky, eye 2 m
above the water, 28 mm lens, looking south-southwest, camera tilted 3 degrees down so the horizon sits
43% from the top of the frame. The high Sun is out of frame, ahead and to the left, so
the light falls from high in front: the water shows no glitter in view.

The sky is a very pale cyan-blue, almost white (#cdece7 at the top). It is a little brighter there
than near the horizon, where it turns a pale grey-green-cream (#ced8c5), and it is evenly bright
from side to side.

Along the whole horizon runs a low, faint range of hills and crater rims 16 to 42 km away, 25 to 63
pixels tall, highest at three-quarters of the way across. The hills are a hazy, desaturated blue-grey
(#8a96a1), only slightly darker than the sky above them, sunlit but softened by the long path
through the air.

The sea is a grey-teal slate (#859895 just below the horizon), much darker than the sky. It
darkens smoothly to #565f57 at the bottom of the frame, with a faint olive tint in the nearest
water. A long ocean swell of 2.2 m significant height (crest to trough), 71 m between
crests, runs away from the camera toward the hills. Its crests lie parallel to the horizon as broad,
low, rounded ridges, so the near water rises and falls in long gentle slopes. Wind ripples from a few
centimetres to a few decimetres long, coarser than ripples on Earth, cover everything under a
4.9 m/s breeze from behind the camera. They break the reflected sky into a fine, high-contrast
mosaic of pale and slate facets; the waves show strongly. Unbroken water, no foam. No boats, no people,
no buildings. Vegetation is omitted for this optical study.

**Frame data.**
- **View:** 205° (south-southwest), pitched 3° down; horizon at row 463. The frame is 39% sky, 4% land and 57% sea.
- **Exposure:** white is 14,100 cd/m².
- **Sky:**
  - top #cdece7 (11,100 cd/m²); edges #d2eee4 and #c6e5e1 (11,300 and 10,400);
  - 15° up #ceece6 (11,100 cd/m²);
  - 8° up #d2e7d9 (10,700 cd/m²); edges #d5e9da and #cae1d7 (10,900 and 10,100);
  - 3° up #d1dfcd (9,960 cd/m²); edges #d3e0ce and #cadacb (10,100 and 9,440);
  - 1° up #cfdac7 (9,520 cd/m²); edges #d1dbc8 and #c8d5c6 (9,690 and 9,050);
  - at the horizon #ced8c5 (9,330 cd/m²); edges #d0d9c6 and #c8d3c3 (9,490 and 8,870).
- **Sun:** 59.9° up at 179°, out of frame; 7.5e+02 million cd/m² (53,000 times white), #ffd979.
- **Earth:** 6.5° up at 100°, out of frame, 40% lit; 170 cd/m² (0 times white), #ffad00.
- **Coast:** skyline above the horizon line, x 0: 27 px (41.5 km, 1,073 m); x 120: 25 px (34.4 km, 802 m); x 240: 33 px (37.5 km, 1,104 m); x 360: 36 px (37.6 km, 1,182 m); x 480: 32 px (30.8 km, 844 m); x 600: 30 px (37.8 km, 1,090 m); x 720: 41 px (29.5 km, 997 m); x 840: 45 px (31.4 km, 1,168 m); x 960: 52 px (31.5 km, 1,321 m); x 1079: 55 px (30.8 km, 1,335 m); x 1199: 54 px (30.1 km, 1,270 m); x 1319: 54 px (30.6 km, 1,282 m); x 1439: 60 px (29.9 km, 1,357 m); x 1559: 55 px (29.9 km, 1,231 m); x 1679: 47 px (29.6 km, 1,050 m); x 1799: 44 px (22.5 km, 687 m); x 1919: 45 px (16.2 km, 469 m). The shore is 3.2 to 20.1 km away. Sunlit faces; haze transmittance 0.672 (green) and 0.413 (blue) at the highest point; colour #8a96a1 (4,190 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #859895 (4,180 cd/m²); edges #899c97 and #7f9393;
  - 1° down #829592 (4,010 cd/m²); edges #879a95 and #7e9291;
  - 3° down #7b8d8a (3,560 cd/m²); edges #83948e and #798c8a;
  - 8° down #70807b (2,850 cd/m²); edges #76857e and #6c7d7b;
  - bottom #565f57 (1,530 cd/m²); edges #636a61 and #515f5d.
- **Waves:** swell 2.2 m, period 16.6 s, 71 m long, travelling toward 207° (2° right of the view's direction); wind 4.9 m/s; mean square slope 0.0373 (rms tilt 10.9°), 0.0246 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1152 (+4 h).
- **Water:** its own light looking straight down, pale yellow-green (204 cd/m²).
- **Tide:** −0.90 m from its mean level.

### WP-2 · The Sun 4° up

*Day 25.3, nine hours before sunset. Sun 3.7° up in the west. 6,700 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the late afternoon, eye 2 m above the water,
28 mm lens, looking west into the low Sun, horizon 53% from the top.

The Sun is a small, intensely bright orange disk, 14 pixels wide (#ff9b00 at its own
exposure, here blown out to white with a warm bloom). It stands 3.7 degrees above the horizon
just right of centre (x = 1078, y = 470), above a dark, hazy line of coastal hills.
Beneath it a broad, dazzling orange glitter path, up to 270 pixels wide, runs from the horizon
down to the bottom of the frame. It is made of countless sparkling highlights on a rough, rippled sea,
brightest at the horizon.

The whole sky is warm:
- pale amber-cream at the top (#ffe59c), slightly brighter above the Sun;
- deepening through #ffd484 at 8 degrees and #ffc17a at 3 degrees;
- to warm orange at the horizon (#f4b474).

There is no halo: the glow is broad and smooth.

The coast across the whole width of the frame, 4 to 8 km away with hills 15 to 30 km behind it, is a
low backlit silhouette 17 to 61 pixels tall, highest left of centre. It is dull red-brown
(#7d5d51), softened by haze.

Away from the glitter the sea is a warm tan-bronze (#a48b5f below the horizon), darkening to
dark olive-brown (#6d5e40) in the lower corners. A long, low swell 2.9 m high and
153 m from crest to crest travels away to the left. Its crests run diagonally, and the glitter path
wobbles and breaks across them. Wind ripples, coarser than on Earth, cover everything under a
4.3 m/s breeze. Unbroken water, no foam, no boats, no people.

**Frame data.**
- **View:** 262° (west), pitched 1° up; horizon at row 569. The frame is 49% sky, 4% land and 47% sea.
- **Exposure:** white is 4,290 cd/m².
- **Sky:**
  - top #ffe59c (3,440 cd/m²); edges #f6d895 and #fddd96 (3,050 and 3,220);
  - 15° up #ffe293 (3,430 cd/m²); edges #f9d791 and #ffdc92 (3,030 and 3,200);
  - 8° up #ffd484 (3,120 cd/m²); edges #faca83 and #ffce84 (2,760 and 2,900);
  - 3° up #ffc17a (2,620 cd/m²); edges #efb87a and #f6bc7a (2,330 and 2,450);
  - 1° up #f8b776 (2,370 cd/m²); edges #e7b076 and #eeb376 (2,110 and 2,220);
  - at the horizon #f4b474 (2,270 cd/m²); edges #e3ac74 and #eaaf74 (2,030 and 2,130).
- **Sun:** at (1078, 470), 3.67° up at 266.5°, 14 by 14 px; 37 million cd/m² (8,600 times white), #ff9b00.
- **Earth:** 12.4° up at 92°, out of frame, 98% lit; 670 cd/m² (0.20 times white), #ffc100.
- **Glitter of the Sun:** across x = 965 to 1235, y = 575 to 1075; peak 95,000 cd/m² (22 times white) at (1075, 575); 10% of the sea in the frame.
- **Coast:** skyline above the horizon line, x 0: 52 px (29.8 km, 1,081 m); x 120: 46 px (20.7 km, 652 m); x 240: 39 px (16.3 km, 441 m); x 360: 52 px (15.6 km, 545 m); x 480: 61 px (15.4 km, 647 m); x 600: 59 px (15.6 km, 641 m); x 720: 57 px (16.0 km, 647 m); x 840: 51 px (17.3 km, 642 m); x 960: 39 px (18.1 km, 535 m); x 1079: 36 px (14.6 km, 383 m); x 1199: 40 px (14.0 km, 407 m); x 1319: 31 px (14.0 km, 321 m); x 1439: 36 px (19.6 km, 535 m); x 1559: 41 px (19.9 km, 596 m); x 1679: 34 px (18.2 km, 450 m); x 1799: 17 px (8.2 km, 93 m); x 1919: 19 px (31.9 km, 590 m). The shore is 4.3 to 10.7 km away. Earthlit faces; haze transmittance 0.817 (green) and 0.637 (blue) at the highest point; colour #7d5d51 (548 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #b49661 (1,380 cd/m²); edges #a48b5f and #b49663;
  - 1° down #b29460 (1,350 cd/m²); edges #a1895d and #b29461;
  - 3° down #a98c5b (1,210 cd/m²); edges #998258 and #aa8d5d;
  - 8° down #a97f4f (1,040 cd/m²); edges #87734e and #947b51;
  - bottom #c67134 (1,030 cd/m²); edges #6d5e40 and #756342.
- **Waves:** swell 2.9 m, period 24.4 s, 153 m long, travelling toward 222° (40° left of the view's direction); wind 4.3 m/s; mean square slope 0.0306 (rms tilt 9.9°), 0.0182 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1296 (-21 h).
- **Water:** its own light looking straight down, pale yellow-green (18.9 cd/m²).
- **Tide:** +0.25 m from its mean level.

### WP-3 · The Sun 4° below the horizon

*Day 26.1, nine hours after sunset. 2,770 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the early evening twilight, eye 2 m above the
water, 28 mm lens, looking west toward where the Sun has set, horizon 58% from the top.
The Sun itself is 4.3 degrees below the horizon.

The sky is a luminous warm gradient, brightest above the place of sunset just right of centre:
- soft amber at the top (#ffe194);
- light orange (#ffce78) at 8 degrees and #ffb76f at 3 degrees;
- coral-orange at the horizon (#f0a76b).

The amber fills the frame from the horizon to its top, 984 to 1,650 cd/m².

A low, hazy silhouette of coastal hills, dull red-brown (#79554a), lies along the horizon across
the frame. It is 14 to 62 pixels tall, highest in the left half (630 m hills at 15 km), lower and fainter
on the right.

The sea is a choppy, wind-roughened surface that reflects the amber sky diffusely, at
59% of its brightness. It is warm tan-orange (#af8a56 just below the horizon),
darkening to deep bronze (#745c39) at the bottom. A long swell 3.3 m high and 153 m
from crest to crest travels away to the left under a 5.0 m/s wind. Ripples cover every slope and
break the reflection into a fine texture of orange highlights and darker brown facets, with no distinct
glitter path. Unbroken water, no foam, no boats, no people.

**Frame data.**
- **View:** 268° (west), pitched 3° up; horizon at row 621. The frame is 54% sky, 4% land and 43% sea.
- **Exposure:** white is 2,080 cd/m².
- **Sky:**
  - top #ffe194 (1,650 cd/m²); edges #fed48b and #ffd88c (1,450 and 1,510);
  - 15° up #ffde87 (1,660 cd/m²); edges #ffd185 and #ffd585 (1,440 and 1,500);
  - 8° up #ffce78 (1,480 cd/m²); edges #ffc277 and #ffc577 (1,270 and 1,330);
  - 3° up #ffb76f (1,180 cd/m²); edges #f0ad6f and #f5b06f (1,030 and 1,070);
  - 1° up #f6ac6c (1,040 cd/m²); edges #e4a46c and #e9a66c (912 and 949);
  - at the horizon #f0a76b (984 cd/m²); edges #dea06a and #e3a26b (865 and 899).
- **Sun:** 4.3° below the horizon at 271°.
- **Earth:** 13.4° up at 92°, out of frame, 99% lit; 790 cd/m² (0.40 times white), #ffc400.
- **Coast:** skyline above the horizon line, x 0: 45 px (23.4 km, 711 m); x 120: 48 px (15.8 km, 489 m); x 240: 62 px (15.4 km, 625 m); x 360: 62 px (15.4 km, 636 m); x 480: 61 px (15.7 km, 656 m); x 600: 55 px (16.1 km, 624 m); x 720: 50 px (17.9 km, 659 m); x 840: 32 px (17.2 km, 426 m); x 960: 39 px (14.5 km, 416 m); x 1079: 39 px (14.0 km, 402 m); x 1199: 29 px (19.6 km, 451 m); x 1319: 38 px (19.9 km, 570 m); x 1439: 39 px (18.9 km, 547 m); x 1559: 20 px (18.1 km, 294 m); x 1679: 17 px (42.4 km, 901 m); x 1799: 18 px (30.8 km, 562 m); x 1919: 14 px (31.6 km, 494 m). The shore is 4.3 to 10.6 km away. Earthlit faces; haze transmittance 0.817 (green) and 0.637 (blue) at the highest point; colour #79554a (230 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #af8a56 (581 cd/m²); edges #a08153 and #ae8957;
  - 1° down #ac8854 (561 cd/m²); edges #9d7e52 and #aa8755;
  - 3° down #a58150 (506 cd/m²); edges #94774d and #a38051;
  - 8° down #917247 (385 cd/m²); edges #836944 and #8f7147;
  - bottom #745c39 (244 cd/m²); edges #725b3b and #7a613d.
- **Waves:** swell 3.3 m, period 24.4 s, 153 m long, travelling toward 223° (45° left of the view's direction); wind 5.0 m/s; mean square slope 0.0422 (rms tilt 11.6°), 0.0264 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1344 (+9 h).
- **Water:** its own light looking straight down, pale yellow-green (7.78 cd/m²).
- **Tide:** +0.39 m from its mean level.
- At this wind a few whitecaps would form on Earth; the study does not model them.

### WP-4 · The Sun 15° below the horizon

*Day 27.1, a day and a half after sunset. 479 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the long evening, eye 2 m above the water,
28 mm lens, looking east, away from the sunset, horizon 58% from the top.

A full Earth hangs in the sky just left of centre, 14.8 degrees above the horizon (centre
x = 952, y = 226). It is a round disk 49 pixels across, bright
amber-gold in colour (#ffc600 at its own exposure), here overexposed to a warm white with a soft
bloom.

The sky is a soft, even, pale golden cream: #e7e7cc at the top and #d8cdac just above the
horizon. It is a fifth duller and greyer toward the horizon, where the Moon's own shadow rises into the
air.

A row of low, flat islets 1.8 to 2.5 km away runs from the left edge to about nine-tenths of the way
across. They are 12 to 60 m high and 8 to 43 pixels tall, highest on the left, and their bases sit on
the water just under the horizon. They are dark muted brown silhouettes (#695847), lit from
behind by the Earth.

The sea is a muted grey-gold, #9a9580 to #b7a080 just below the horizon, slightly
warmer and brighter in a broad soft column below the Earth. It darkens to olive-grey
(#726b58) at the bottom. There is no sharp glitter path: the earthlight is faint against the
reflected sky. A long swell 3.3 m high, 153 m from crest to crest, comes toward the camera
from the front left, its crests running diagonally. Wind ripples under a 4.7 m/s breeze give
the water a fine, matte texture. Unbroken water, no foam, no boats, no people.

**Frame data.**
- **View:** 92° (east), pitched 3° up; horizon at row 621. The frame is 55% sky, 2% land and 43% sea.
- **Exposure:** white is 99.0 cd/m².
- **Sky:**
  - top #e7e7cc (77.5 cd/m²); edges #eeebcc and #ebe8cb (81.4 and 79.1);
  - 15° up #e8e3c4 (75.4 cd/m²); edges #efe9c8 and #ebe6c6 (79.9 and 77.6);
  - 8° up #e5dcbb (70.5 cd/m²); edges #ebe1be and #e8debc (74.6 and 72.5);
  - 3° up #ded3b2 (64.9 cd/m²); edges #e4d8b4 and #e1d6b3 (68.3 and 66.6);
  - 1° up #dacfad (61.9 cd/m²); edges #e0d3b0 and #ddd1af (65.1 and 63.5);
  - at the horizon #d8cdac (60.6 cd/m²); edges #ded1ae and #dbcfad (63.6 and 62.1).
- **Sun:** 15.0° below the horizon at 277°.
- **Earth:** at (952, 226), 14.76° up at 91.7°, 49 by 50 px, 100% lit, lit side toward the left; 900 cd/m² (9.1 times white), #ffc600.
- **Coast:** skyline above the horizon line, x 0: 15 px (2.5 km, 22 m); x 120: 21 px (2.4 km, 29 m); x 240: 30 px (2.3 km, 42 m); x 360: 39 px (2.2 km, 54 m); x 480: 43 px (2.2 km, 60 m); x 600: 41 px (2.1 km, 58 m); x 720: 36 px (2.1 km, 50 m); x 840: 32 px (2.1 km, 45 m); x 960: 30 px (2.2 km, 44 m); x 1079: 31 px (2.2 km, 46 m); x 1199: 32 px (2.4 km, 49 m); x 1319: 31 px (2.4 km, 48 m); x 1439: 28 px (2.5 km, 43 m); x 1559: 21 px (2.4 km, 32 m); x 1679: 8 px (2.4 km, 12 m). The shore is 1.8 to 2.2 km away. Backlit; haze transmittance 0.972 (green) and 0.939 (blue) at the highest point; colour #695847 (10.4 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #b7a080 (36.4 cd/m²); edges #9a9580 and #a09c86;
  - 1° down #b39e7f (35.4 cd/m²); edges #98937e and #9e9a84;
  - 3° down #a89679 (31.1 cd/m²); edges #918d79 and #98937e;
  - 8° down #90846c (23.2 cd/m²); edges #827e6b and #87836f;
  - bottom #726b58 (14.6 cd/m²); edges #726f5d and #747160.
- **Waves:** swell 3.3 m, period 24.4 s, 153 m long, travelling toward 223° (131° right of the view's direction); wind 4.7 m/s; mean square slope 0.039 (rms tilt 11.2°), 0.0232 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1344 (-15 h).
- **Water:** its own light looking straight down, pale yellow-green (1.34 cd/m²).
- **Tide:** +0.58 m from its mean level.

### WP-5 · The Sun 35° below the horizon

*Day 29.0, three and a third days after sunset. 9.6 lux; between day and night vision.*

**Image prompt.** Photograph from a small boat at sea late in the long evening, eye 2 m above the water,
28 mm lens, looking west-northwest, horizon 58% from the top; a long exposure.

A great twilight arch glows over the western horizon. It is a broad dome of deep gold and orange
centred a little right of the middle of the frame:
- brightest amber-gold about 8 degrees up (#ffcc42);
- vivid orange at the horizon (#ff9d32);
- duller ochre-orange at the top of the frame (#ecbc6b);
- redder and darker toward the left edge (#f8873f low down).

The Sun is 35 degrees below the horizon, and this light comes from the very high air it still reaches.

A hazy silhouette of coastal hills, dark red-brown (#7c422f), runs along the horizon:
- tallest at the far left, 66 pixels (640 m hills at 16 km);
- sloping down to a thin, faint line of distant ranges only 6 to 15 pixels tall on the right;
- with a few low islands near the right edge.

The sea is roughened by a swell 1.8 m high and 105 m between crests, running away to the
left, and by wind ripples under a 2.6 m/s breeze. It reflects the arch at about half its
brightness (54% just below the horizon): warm orange (#d38a39) at the centre,
duller toward the edges (#a07039 on the left), and deep brown-orange (#825625) at
the bottom. The tilted facets break the reflection into a grainy orange texture. To the eye at this
light level colours are muted; the photograph records them. No stars, no foam, no boats, no people.

**Frame data.**
- **View:** 285° (west-northwest), pitched 3° up; horizon at row 621. The frame is 55% sky, 3% land and 43% sea.
- **Exposure:** white is 16.5 cd/m².
- **Sky:**
  - top #ecbc6b (9.02 cd/m²); edges #d3a762 and #eab564 (6.95 and 8.45);
  - 15° up #ffcc5e (11.6 cd/m²); edges #e4ab5b and #ffbb5d (7.64 and 9.50);
  - 8° up #ffcc42 (13.4 cd/m²); edges #ffa94b and #ffba47 (8.31 and 10.7);
  - 3° up #ffb432 (11.9 cd/m²); edges #ff9741 and #ffa63b (7.37 and 9.51);
  - 1° up #ffa432 (10.4 cd/m²); edges #fc8c3f and #ff983a (6.55 and 8.40);
  - at the horizon #ff9d32 (9.81 cd/m²); edges #f8873f and #ff9239 (6.20 and 7.91).
- **Sun:** 35.3° below the horizon at 291°.
- **Earth:** 17.2° up at 92°, out of frame, 96% lit; 970 cd/m² (59 times white), #ffca08.
- **Coast:** skyline above the horizon line, x 0: 66 px (15.8 km, 643 m); x 120: 60 px (16.7 km, 642 m); x 240: 55 px (18.2 km, 670 m); x 360: 35 px (17.2 km, 431 m); x 480: 39 px (14.6 km, 404 m); x 600: 41 px (14.0 km, 408 m); x 720: 29 px (14.0 km, 307 m); x 840: 35 px (19.6 km, 535 m); x 960: 38 px (18.9 km, 559 m); x 1079: 19 px (18.1 km, 301 m); x 1199: 15 px (32.0 km, 562 m); x 1319: 13 px (31.2 km, 496 m); x 1439: 6 px (32.4 km, 386 m); x 1559: 9 px (9.1 km, 66 m); x 1679: 13 px (9.4 km, 90 m); x 1799: 10 px (9.7 km, 74 m); x 1919: 8 px (4.8 km, 23 m). The shore is 4.4 to 7.3 km away. Earthlit faces; haze transmittance 0.813 (green) and 0.631 (blue) at the highest point; colour #7c422f (1.38 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #d38a39 (5.29 cd/m²); edges #a07039 and #b97c39;
  - 1° down #cf8738 (5.07 cd/m²); edges #9f6f38 and #b67a38;
  - 3° down #c58035 (4.53 cd/m²); edges #976835 and #af7435;
  - 8° down #ad6f2e (3.35 cd/m²); edges #865b2e and #9a652e;
  - bottom #825625 (1.89 cd/m²); edges #6f4d27 and #7f5527.
- **Waves:** swell 1.8 m, period 20.1 s, 105 m long, travelling toward 221° (64° left of the view's direction); wind 2.6 m/s; mean square slope 0.0336 (rms tilt 10.4°), 0.0067 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1419 (+13 h).
- **Water:** its own light looking straight down, pale yellow-green (0.0268 cd/m²).
- **Stars:** none recorded; the eye's limit here is magnitude -3.9.
- **Tide:** +0.97 m from its mean level.

### WP-6 · The darkest hour

*Day 4.5, eight days after sunset; the Sun 60° below the horizon. 0.48 lux, from the Earth and the
twilight; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat on a glassy night sea, eye 2 m above the
water, 28 mm lens, looking east toward the Earth, horizon 60% from the top.

The Earth hangs near the top centre, 20.3 degrees up (centre x = 970,
y = 100). It is a gibbous disk 50 pixels across, 58% lit, its
bright side toward the lower left. It is golden-amber in colour (#ffcd2d at its own
exposure) but here overexposed to near white with a bloom, the unlit part invisible against the sky.

Straight down from it a vivid golden-yellow glitter path, up to 260 pixels wide, runs from
the horizon to the bottom of the frame (#ffe93c at its brightest). Its edge is soft and its
interior a ladder of long smooth highlights, because the sea is glassy.

The sky is pale and warm:
- cream-yellow overhead (#dae0c6);
- a soft peach-orange glow along the left horizon (#ffc596), on the side toward the north,
  where the twilight of the deep-set Sun lingers;
- a greyer, duller yellow-beige toward the right (#c9b796).

A row of low flat islets 2 km away, 20 to 60 m high and 15 to 43 pixels tall, stands as near-black
silhouettes (#443c3a) from the left edge to four-fifths of the way across.

The calm sea mirrors the sky, peach-tan (#e4ba8c) on the left and grey-beige
(#b0a78a) on the right, darkening downward. It moves in a long, low, glassy swell
1.1 m high and 71 m from crest to crest, coming from the front left, with no wind ripples.

To the eye this is a dim, almost colourless night with the bright Earth and its path. No stars show
against this bright sky. No foam, no boats, no people.

**Frame data.**
- **View:** 94° (east), pitched 4° up; horizon at row 648. The frame is 58% sky, 2% land and 40% sea.
- **Exposure:** white is 0.237 cd/m².
- **Sky:**
  - top #dae0c6 (0.171 cd/m²); edges #eee6c1 and #cbd4bd (0.186 and 0.150);
  - 15° up #e7e0ba (0.176 cd/m²); edges #fbe7ba and #d0d3b6 (0.192 and 0.150);
  - 8° up #eed8ab (0.167 cd/m²); edges #ffe0aa and #d3cba8 (0.190 and 0.141);
  - 3° up #eacb9f (0.148 cd/m²); edges #ffd29d and #cfc09d (0.171 and 0.126);
  - 1° up #e5c49a (0.138 cd/m²); edges #ffc998 and #cbb998 (0.158 and 0.118);
  - at the horizon #e2c197 (0.134 cd/m²); edges #ffc596 and #c9b796 (0.152 and 0.114).
- **Sun:** 60.0° below the horizon at 23°.
- **Earth:** at (970, 100), 20.30° up at 94.4°, 50 by 52 px, 58% lit, lit side toward the lower left; 730 cd/m² (3,100 times white), #ffcd2d.
- **Glitter of the Earth:** across x = 845 to 1105, y = 655 to 1075; peak 2.4 cd/m² (10 times white) at (975, 1075); 9% of the sea in the frame.
- **Coast:** skyline above the horizon line, x 0: 18 px (2.4 km, 25 m); x 120: 25 px (2.4 km, 35 m); x 240: 36 px (2.2 km, 48 m); x 360: 42 px (2.2 km, 57 m); x 480: 43 px (2.1 km, 58 m); x 600: 39 px (2.1 km, 54 m); x 720: 34 px (2.1 km, 46 m); x 840: 31 px (2.1 km, 43 m); x 960: 30 px (2.2 km, 44 m); x 1079: 31 px (2.3 km, 47 m); x 1199: 31 px (2.5 km, 50 m); x 1319: 29 px (2.4 km, 45 m); x 1439: 25 px (2.5 km, 38 m); x 1559: 15 px (2.4 km, 21 m). The shore is 1.8 to 2.1 km away. Backlit; haze transmittance 0.972 (green) and 0.939 (blue) at the highest point; colour #443c3a (0.0112 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #ffe88d (0.225 cd/m²); edges #e4ba8c and #b0a78a;
  - 1° down #ffea8a (0.237 cd/m²); edges #e0b78a and #aba386;
  - 3° down #ffed7a (0.308 cd/m²); edges #d4ad83 and #a39c81;
  - 8° down #ffeb50 (0.834 cd/m²); edges #b09572 and #8a8670;
  - bottom #ffe93c (2.21 cd/m²); edges #8d7e62 and #727260.
- **Waves:** swell 1.1 m, period 16.6 s, 71 m long, travelling toward 219° (125° right of the view's direction); wind 1.6 m/s; mean square slope 0.0046 (rms tilt 3.9°), 0.0000 of it from waves shorter than 1 m; nearside wave run, restart file at hour 816 (+0 h).
- **Water:** its own light looking straight down, pale yellow-green (0.00140 cd/m²).
- **Stars:** none recorded; the eye's limit here is magnitude 0.7.
- **Tide:** +1.35 m from its mean level.

## The eastern Smythii headland

The camera floats about 1 km off the west face of a narrow headland in the eastern Smythii basin
(93.22° E, 2.36° N), over water about 270 m deep. Its waves are those of the shore history's west-face
station, 1 km inshore in 22 m of water. It looks west across the open sea; the headland is behind the
camera. No land enters these frames. The Earth stands low in the west, rising and setting with the
libration: below the horizon for the first five scenes, just risen at the darkest hour.

### SH-1 · The Sun at its highest

*Day 4.6, seven and a half days before sunset. Sun 86.8° up, almost overhead. 87,200 lux; day
vision.*

**Image prompt.** Photograph from a small boat on the open sea at midday, eye 2 m above the water, 28 mm
lens, looking west across empty water to a perfectly straight horizon, camera tilted 3 degrees down,
horizon 43% from the top. The Sun is almost directly overhead, out of frame, its glitter
beneath the boat out of view.

The sky is a very pale, luminous cyan-blue (#c6eef6 at the top), almost white, fading to a pale
grey-green-cream at the horizon (#cfe0d7). It is evenly bright across the frame: a clean,
washed-out noon sky.

The sea is a cool grey-teal (#89a3a7 just below the horizon), 48% as bright as
the sky. It darkens steadily to deep slate-teal (#4e5e5f) at the bottom of the frame, with a
faint olive cast in the nearest water. A moderate swell 0.6 m high and 23 m from crest to
crest rolls straight toward the camera from the horizon. Its crests run parallel to the horizon, and
the camera looks into their faces. Lively wind ripples under a 3.1 m/s breeze blowing toward the
camera cover the surface in a fine, crisp texture of pale and dark facets. Empty sea, unbroken water,
no foam, no land, no boats, no people.

**Frame data.**
- **View:** 270° (west), pitched 3° down; horizon at row 463. The frame is 43% sky, 0% land and 57% sea.
- **Exposure:** white is 13,100 cd/m².
- **Sky:**
  - top #c6eef6 (10,400 cd/m²); edges #c9eef4 and #c8edf3 (10,500 and 10,400);
  - 15° up #c6eef5 (10,400 cd/m²);
  - 8° up #ceecea (10,300 cd/m²); edges #ceecea and #ceebea (10,400 and 10,300);
  - 3° up #d0e6df (9,840 cd/m²); edges #d0e6df and #d0e6df (9,860 and 9,830);
  - 1° up #cfe2d9 (9,510 cd/m²); edges #d0e2d9 and #cfe2d9 (9,520 and 9,510);
  - at the horizon #cfe0d7 (9,360 cd/m²); edges #cfe0d7 and #cfe0d6 (9,370 and 9,360).
- **Sun:** 86.8° up at 178°, out of frame; 8.3e+02 million cd/m² (63,000 times white), #ffda84.
- **Earth:** 7.3° below the horizon at 276°, 57% lit, lit side toward the top.
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #89a3a7 (4,470 cd/m²); edges #8ca5a9 and #8aa3a7;
  - 1° down #879fa4 (4,280 cd/m²); edges #89a2a6 and #88a1a5;
  - 3° down #80989b (3,840 cd/m²); edges #859c9f and #839b9e;
  - 8° down #738789 (3,000 cd/m²); edges #768b8c and #758a8c;
  - bottom #4e5e5f (1,370 cd/m²); edges #566768 and #556768.
- **Waves:** swell 0.6 m, period 9.3 s, 23 m long, travelling toward 94° (176° left of the view's direction); wind 3.1 m/s; mean square slope 0.0253 (rms tilt 9.0°), 0.0098 of it from waves shorter than 1 m; Smythii shore history, west-face station, hour 819.
- **Water:** its own light looking straight down, pale yellow-green (246 cd/m²).
- **Tide:** +0.21 m from its mean level.

### SH-2 · The Sun 4° up

*Day 11.7, eight hours before sunset. Sun 4.0° up in the west. 6,800 lux; day vision.*

**Image prompt.** Photograph from a small boat on a glassy sea in the late afternoon, eye 2 m above the
water, 28 mm lens, looking west into the low Sun over an empty, perfectly flat horizon, horizon
53% from the top.

The Sun is a small, dazzling orange disk, 14 pixels across (#ff9d00 at its own
exposure, blown out here), 4.0 degrees above the horizon at the centre of the frame
(x = 980, y = 460). Below it lies an almost mirror-calm sea. The Sun's reflection is
a narrow, intensely brilliant column of molten orange, no more than 70 pixels wide, running
from the horizon straight down to the bottom of the frame. It is broken into long, smooth, stretched
highlights by a low swell; its core is blinding white-orange and its edges sharp.

The sky is a warm amber-cream at the top (#ffe59d), deepening to #ffd585 at 8 degrees and to
soft orange at the horizon (#f4b475). It is brightest above the Sun, with no halo.

Beside the glitter column the calm water mirrors the sky at 89% of its brightness: warm
apricot (#d4a86d) just below the horizon, darkening to bronze-olive (#7c6a46) at
the bottom. A gentle, glassy swell 0.9 m high and 48 m from crest to crest rolls toward the
camera from straight ahead. Its smooth undulations bend the reflections into long wavy bands. There are
no wind ripples. Empty sea, no land, no foam, no boats, no people.

**Frame data.**
- **View:** 268° (west), pitched 1° up; horizon at row 569. The frame is 53% sky, 0% land and 47% sea.
- **Exposure:** white is 4,400 cd/m².
- **Sky:**
  - top #ffe59d (3,530 cd/m²); edges #f9db96 and #fadc96 (3,210 and 3,240);
  - 15° up #ffe394 (3,530 cd/m²); edges #fbd992 and #fdda92 (3,190 and 3,220);
  - 8° up #ffd585 (3,210 cd/m²); edges #fdcc84 and #fecd84 (2,900 and 2,930);
  - 3° up #ffc27a (2,700 cd/m²); edges #f2ba7a and #f4bb7a (2,450 and 2,470);
  - 1° up #f8b877 (2,450 cd/m²); edges #eab276 and #ebb276 (2,230 and 2,240);
  - at the horizon #f4b475 (2,350 cd/m²); edges #e7ae75 and #e8af75 (2,140 and 2,150).
- **Sun:** at (980, 460), 4.02° up at 268.8°, 14 by 14 px; 40 million cd/m² (9,100 times white), #ff9d00.
- **Earth:** 8.5° below the horizon at 267°, 1% lit, lit side toward the top.
- **Glitter of the Sun:** across x = 955 to 1025, y = 575 to 1075; peak 1.2 million cd/m² (280 times white) at (975, 575); 3% of the sea in the frame.
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #e4b06d (2,140 cd/m²); edges #d4a86d and #daaa6e;
  - 1° down #efb06a (2,220 cd/m²); edges #d0a56b and #d7a86c;
  - 3° down #ffc938 (5,470 cd/m²); edges #c69d66 and #cb9f67;
  - 8° down #ffbb00 (12,500 cd/m²); edges #a78858 and #ab8a58;
  - bottom #8f6a40 (726 cd/m²); edges #7c6a46 and #7d6b46.
- **Waves:** swell 0.9 m, period 13.7 s, 48 m long, travelling toward 75° (167° right of the view's direction); wind 1.2 m/s; mean square slope 0.0027 (rms tilt 3.0°), 0.0000 of it from waves shorter than 1 m; Smythii shore history, west-face station, hour 989.
- **Water:** its own light looking straight down, pale yellow-green (19.6 cd/m²).
- **Tide:** +0.31 m from its mean level.

### SH-3 · The Sun 4° below the horizon

*Day 12.3, eight hours after sunset. 2,800 lux; day vision.*

**Image prompt.** Photograph from a small boat on a calm sea in the early evening twilight, eye 2 m
above the water, 28 mm lens, looking west to the place of sunset over a flat, empty horizon, horizon
58% from the top. The Sun is 4 degrees below the horizon at the centre of the frame.

The sky glows amber and orange:
- soft amber at the top (#ffe194), brightest above the centre;
- light orange (#ffce78) at 8 degrees;
- coral-orange just above the horizon (#f1a86b).

The sea is nearly calm and mirrors this sky almost exactly: apricot-orange (#e8a864) just below
the horizon, 96% as bright as the sky above, fading to warm bronze (#866b40) at
the bottom. The horizon line is barely visible, a faint darkening between two glowing halves of the
frame. A low glassy swell 0.8 m high and 58 m between crests rolls toward the camera.
Faint patches of capillary ripples, the first a 1.9 m/s breeze raises, dull the mirror slightly.
The swell's long smooth undulations stretch the reflected glow into gentle horizontal bands. Empty sea,
no land, no foam, no boats, no people.

**Frame data.**
- **View:** 268° (west), pitched 3° up; horizon at row 621. The frame is 57% sky, 0% land and 43% sea.
- **Exposure:** white is 2,140 cd/m².
- **Sky:**
  - top #ffe194 (1,700 cd/m²); edges #ffd58c and #ffd78c (1,520 and 1,540);
  - 15° up #ffdf88 (1,710 cd/m²); edges #ffd385 and #ffd486 (1,500 and 1,530);
  - 8° up #ffce78 (1,520 cd/m²); edges #ffc378 and #ffc578 (1,340 and 1,360);
  - 3° up #ffb76f (1,220 cd/m²); edges #f2af6f and #f4b06f (1,080 and 1,100);
  - 1° up #f6ac6c (1,070 cd/m²); edges #e6a56c and #e8a66c (954 and 969);
  - at the horizon #f1a86b (1,020 cd/m²); edges #e0a16b and #e2a26b (906 and 920).
- **Sun:** 4.1° below the horizon at 269°.
- **Earth:** 7.6° below the horizon at 266°, 0% lit, lit side toward the upper right.
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #e8a864 (981 cd/m²); edges #d89f63 and #dba164;
  - 1° down #e4a562 (946 cd/m²); edges #d59d62 and #d79e62;
  - 3° down #d89d5d (844 cd/m²); edges #ca955d and #cc965d;
  - 8° down #b88a50 (617 cd/m²); edges #ac8350 and #ae8450;
  - bottom #866b40 (342 cd/m²); edges #8a6d43 and #8c6e43.
- **Waves:** swell 0.8 m, period 15.1 s, 58 m long, travelling toward 71° (163° right of the view's direction); wind 1.9 m/s; mean square slope 0.0028 (rms tilt 3.0°), 0.0015 of it from waves shorter than 1 m; Smythii shore history, west-face station, hour 1005.
- **Water:** its own light looking straight down, pale yellow-green (8.07 cd/m²).
- **Tide:** +0.25 m from its mean level.

### SH-4 · The Sun 15° below the horizon

*Day 13.3, 30 hours after sunset. 450 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the long evening, eye 2 m above the water,
28 mm lens, looking west over an empty horizon, horizon 58% from the top.

A deep amber-orange twilight fills the western sky:
- warm amber at the top (#ffd980);
- golden orange (#ffcb61) at 8 degrees;
- a rich coral-orange band along the horizon (#ff9e55), brightest at the centre and duller toward
  both edges.

The sea, gently rippled by a 2.8 m/s breeze, reflects the sky at 80% of its
brightness. It is amber-orange (#d6974c) below the horizon, softening into a smooth, slightly
grainy sheen, then darkening to brown-gold (#856032) at the bottom. A swell 1.1 m high
and 69 m from crest to crest comes from the front left, its crests running diagonally. Empty sea,
no land, no foam, no boats, no people.

**Frame data.**
- **View:** 268° (west), pitched 3° up; horizon at row 621. The frame is 57% sky, 0% land and 43% sea.
- **Exposure:** white is 527 cd/m².
- **Sky:**
  - top #ffd980 (401 cd/m²); edges #ffc976 and #ffcc77 (340 and 352);
  - 15° up #ffdb72 (427 cd/m²); edges #ffc870 and #ffcb70 (346 and 358);
  - 8° up #ffcb61 (390 cd/m²); edges #ffba61 and #ffbd61 (315 and 326);
  - 3° up #ffb058 (305 cd/m²); edges #faa359 and #ffa559 (249 and 258);
  - 1° up #ffa356 (261 cd/m²); edges #ec9857 and #f09a57 (216 and 223);
  - at the horizon #ff9e55 (245 cd/m²); edges #e69456 and #ea9656 (203 and 210).
- **Sun:** 15.2° below the horizon at 270°.
- **Earth:** 6.1° below the horizon at 265°, 1% lit, lit side toward the lower right.
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #d6974c (195 cd/m²); edges #ba864a and #bf8a4a;
  - 1° down #d2944b (186 cd/m²); edges #b88549 and #bd8849;
  - 3° down #c88c47 (166 cd/m²); edges #ae7d45 and #b48146;
  - 8° down #ad7a3e (122 cd/m²); edges #986e3c and #9c703d;
  - bottom #856032 (71.4 cd/m²); edges #7f5d33 and #825f34.
- **Waves:** swell 1.1 m, period 16.6 s, 69 m long, travelling toward 44° (136° right of the view's direction); wind 2.8 m/s; mean square slope 0.0268 (rms tilt 9.3°), 0.0080 of it from waves shorter than 1 m; Smythii shore history, west-face station, hour 1027.
- **Water:** its own light looking straight down, pale yellow-green (1.30 cd/m²).
- **Tide:** +0.16 m from its mean level.

### SH-5 · The Sun 35° below the horizon

*Day 14.9, nearly three days after sunset. 9.1 lux; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat at sea late in the long evening, eye 2 m
above the water, 28 mm lens, looking west over an empty horizon, horizon 58% from the
top.

A deep twilight arch spans the western sky:
- glowing amber-gold about 8 degrees up at the centre (#ffcc3d);
- vivid orange at the horizon (#ff9c2c);
- duller ochre at the top of the frame (#eebc68);
- redder and darker toward both edges (#fd8739 low on the left).

It is the light of the very high air, which still sees the Sun 35 degrees below the horizon.

A confused low swell and a 3.2 m/s breeze roughen the sea. It reflects the arch at
45% of its brightness: orange (#c18035) just below the horizon, deepening to dark
brown-orange (#7e5222) at the bottom. Its surface is a soft-grained texture of orange and
brown facets. The Earth, a thin crescent, is 3 degrees below the horizon and not visible. Muted to the
eye; the photograph records the colour. Empty sea, no land, no stars, no foam, no boats, no people.

**Frame data.**
- **View:** 266° (west), pitched 3° up; horizon at row 621. The frame is 57% sky, 0% land and 43% sea.
- **Exposure:** white is 17.6 cd/m².
- **Sky:**
  - top #eebc68 (9.66 cd/m²); edges #d7a85e and #e8b260 (7.58 and 8.77);
  - 15° up #ffcc5a (12.5 cd/m²); edges #e9ac58 and #fdb859 (8.38 and 9.86);
  - 8° up #ffcc3d (14.4 cd/m²); edges #ffaa46 and #ffb743 (9.19 and 11.1);
  - 3° up #ffb42c (12.8 cd/m²); edges #ff983b and #ffa337 (8.13 and 9.80);
  - 1° up #ffa32c (11.2 cd/m²); edges #ff8c3a and #ff9535 (7.20 and 8.63);
  - at the horizon #ff9c2c (10.5 cd/m²); edges #fd8739 and #ff9035 (6.80 and 8.13).
- **Sun:** 34.9° below the horizon at 270°.
- **Earth:** 3.2° below the horizon at 264°, 8% lit, lit side toward the bottom.
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #c18035 (4.74 cd/m²); edges #9e6d36 and #ad7536;
  - 1° down #bd7d34 (4.55 cd/m²); edges #9b6b34 and #ab7335;
  - 3° down #b77731 (4.13 cd/m²); edges #916431 and #a46e32;
  - 8° down #a2682a (3.12 cd/m²); edges #81572a and #91602b;
  - bottom #7e5222 (1.86 cd/m²); edges #6f4b24 and #7a5124.
- **Waves:** swell 0.5 m, period 15.1 s, 58 m long, travelling toward 41° (135° right of the view's direction); wind 3.2 m/s; mean square slope 0.0378 (rms tilt 11.0°), 0.0106 of it from waves shorter than 1 m; Smythii shore history, west-face station, hour 1066.
- **Water:** its own light looking straight down, pale yellow-green (0.0264 cd/m²).
- **Stars:** none recorded; the eye's limit here is magnitude -4.0.
- **Tide:** −0.03 m from its mean level.

### SH-6 · The darkest hour: the Earth rising

*Day 18.1, six days after sunset; the Sun 74° below the horizon. 0.085 lux; night vision.*

**Image prompt.** Long-exposure night photograph from a small boat on a mirror-calm sea, eye 2 m above
the water, 28 mm lens, looking west, horizon 55% from the top.

The Earth has just risen. It sits 1.6 degrees above the horizon at the centre of the frame
(centre x = 958, y = 552), a disk 52 pixels across,
38% lit. Its lit part is a thick crescent along the bottom edge of the disk,
like a bowl. It is deep orange (#ff8f00) from the long path through the air, overexposed with a
bloom in this exposure. The rest of the disk is invisible against the sky.

Its reflection runs from the horizon straight down toward the camera as a narrow, vivid orange column
about 60 pixels wide, broken into a few long smooth bands by a very low swell. It is an
almost perfect mirror image, reaching to row 825, about halfway from the horizon to the
bottom of the frame.

The sky is a faint warm glow:
- pale apricot-cream overhead (#fce4ae);
- golden-orange at 8 degrees (#ffd592);
- soft salmon-orange at the horizon (#fcb381), lit by the low Earth.

8 stars are recorded: Aldebaran at (906, 84), a small point near the top of the frame,
and fainter ones at (1713, 428); (1375, 250); (1612, 167); (1891, 180); (623, 37); (582, 35); (485, 29).

The water away from the Earth's reflection is glassy and mirrors the sky at 94% of its
brightness: salmon-tan (#dea97a) at the horizon, darkening to dark olive-brown
(#806e50) at the bottom. A very low swell 0.3 m high and 48 m between crests
barely moves it, without ripples.

To the eye this scene is a dark, almost colourless night with the orange Earth and its path. Empty sea,
no land, no foam, no boats, no people.

**Frame data.**
- **View:** 264.5° (west), pitched 2° up; horizon at row 595. The frame is 55% sky, 0% land and 45% sea.
- **Exposure:** white is 0.0596 cd/m².
- **Sky:**
  - top #fce4ae (0.0473 cd/m²); edges #f3daa6 and #f6dba6 (0.0429 and 0.0436);
  - 15° up #ffe3a3 (0.0480 cd/m²); edges #f8d8a1 and #fcdaa1 (0.0428 and 0.0436);
  - 8° up #ffd592 (0.0445 cd/m²); edges #fecb91 and #ffcd91 (0.0394 and 0.0404);
  - 3° up #ffc187 (0.0376 cd/m²); edges #f5b987 and #fbba87 (0.0332 and 0.0343);
  - 1° up #ffb783 (0.0340 cd/m²); edges #edaf83 and #f3b182 (0.0301 and 0.0310);
  - at the horizon #fcb381 (0.0325 cd/m²); edges #e9ac81 and #efad81 (0.0288 and 0.0297).
- **Sun:** 74.3° below the horizon at 274°.
- **Earth:** at (958, 552), 1.55° up at 264.4°, 52 by 52 px, 38% lit, lit side toward the bottom; 59 cd/m² (980 times white), #ff8f00.
- **Glitter of the Earth:** across x = 935 to 995, y = 595 to 825; peak 10 cd/m² (170 times white) at (955, 595); 2% of the sea in the frame.
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #ffa400 (6.19 cd/m²); edges #dea97a and #e8ab7c;
  - 1° down #ffa400 (6.04 cd/m²); edges #daa678 and #e4a97a;
  - 3° down #ffa400 (3.72 cd/m²); edges #ce9e72 and #d7a173;
  - 8° down #ffb919 (0.0835 cd/m²); edges #ac8a62 and #b18c62;
  - bottom #7a6a4c (0.00896 cd/m²); edges #806e50 and #837050.
- **Waves:** swell 0.3 m, period 13.7 s, 48 m long, travelling toward 61° (156° right of the view's direction); wind 0.6 m/s; mean square slope 0.0014 (rms tilt 2.1°), 0.0000 of it from waves shorter than 1 m; Smythii shore history, west-face station, hour 1144.
- **Water:** its own light looking straight down, pale green (0.000200 cd/m²).
- **Stars:** 8 recorded (Aldebaran at (906, 84); others at (1713, 428); (1375, 250); (1612, 167); (1891, 180); (623, 37); (582, 35); (485, 29)); the eye sees none, its limit here magnitude 2.1.
- **Tide:** −0.30 m from its mean level.

## Southern Mare Nubium, by Pitatus

The coast with the largest tide on the exposed tenth of the nearside's shores: this month's range is
7.7 m. The wave grid's shore point falls on the west shore of an island 4 by 6 km, rising to about
290 m. The camera floats 1.5 km off that shore (11.93° W, 28.11° S) in 65 m of water, with open sea to
the west and north. The mainland coast, with hills and crater rims of 400 to 1,100 m, lies 12 to 30 km
to the south and southwest. The Earth stands high in the northeast, 53 to 65 degrees up, out of every
frame.

### SN-1 · The Sun at its highest

*Day 13.2, seven and a half days before sunset. Sun 63.0° up in the north, behind the camera.
75,000 lux; day vision.*

**Image prompt.** Photograph from a small boat on a choppy sea at midday, eye 2 m above the water,
28 mm lens, looking south with the high Sun behind the camera, tilted 3 degrees down, horizon
43% from the top.

On the left rises the sunlit island, 1.7 to 5 km away: a smooth-shouldered hill of bare, dark
grey-brown ground (#847465). It rises to 120 pixels above the horizon at the left edge (280 m
high) and slopes down to the water about a fifth of the way across the frame. It is front-lit and only
lightly hazed. It is high tide, +4.08 m from the mean level, and the sea laps high on the island's
lowest slopes. Beyond it, across the rest of the horizon, runs a faint low line of mainland hills 17 to
30 km away, 38 to 65 pixels tall, pale hazy blue-grey.

The sky is very pale cyan-blue, almost white (#c6ecef at the top), fading to pale
grey-green-cream at the horizon (#d8e3d1).

The sea is the roughest of the 24 scenes:
- a 6.0 m/s wind blows from behind the camera, over the right shoulder, raising steep, crowded
  wind waves and ripples over a 1.2 m swell 33 m from crest to crest that runs away to the
  left;
- the surface is a busy, high-contrast mosaic of pale sky-reflecting facets and dark slate troughs;
- grey-teal (#7b9295) below the horizon, darkening to deep slate (#4b5a58) at the
  bottom, with a faint olive cast in the nearest water.

Unbroken water, no foam. No boats, no people. Vegetation is omitted for this optical study.

**Frame data.**
- **View:** 180° (south), pitched 3° down; horizon at row 463. The frame is 38% sky, 5% land and 57% sea.
- **Exposure:** white is 12,400 cd/m².
- **Sky:**
  - top #c6ecef (9,700 cd/m²); edges #c9ecec and #c9ecec (9,740 and 9,750);
  - 15° up #c7ecee (9,740 cd/m²);
  - 8° up #d4ede4 (9,950 cd/m²); edges #d2ece3 and #d2ece3 (9,830 and 9,840);
  - 3° up #d9e8d9 (9,650 cd/m²); edges #d6e7d9 and #d6e7d9 (9,500 and 9,520);
  - 1° up #d9e5d4 (9,380 cd/m²); edges #d6e3d3 and #d7e4d3 (9,220 and 9,240);
  - at the horizon #d8e3d1 (9,250 cd/m²); edges #d6e2d1 and #d6e2d1 (9,090 and 9,110).
- **Sun:** 63.0° up at 1°, out of frame; 7.6e+02 million cd/m² (61,000 times white), #ffd97c.
- **Earth:** 65.4° up at 23°, out of frame, 1% lit; 1,600 cd/m² (0.10 times white), #ffdc94.
- **Coast:** skyline above the horizon line, x 0: 120 px (4.2 km, 284 m); x 120: 104 px (4.5 km, 271 m); x 240: 82 px (4.8 km, 233 m); x 360: 54 px (4.9 km, 166 m); x 480: 48 px (23.4 km, 832 m); x 600: 53 px (23.4 km, 923 m); x 720: 44 px (29.9 km, 1,077 m); x 840: 40 px (26.9 km, 870 m); x 960: 46 px (24.2 km, 863 m); x 1079: 46 px (22.9 km, 819 m); x 1199: 38 px (19.1 km, 555 m); x 1319: 44 px (20.1 km, 663 m); x 1439: 46 px (19.5 km, 643 m); x 1559: 48 px (17.3 km, 575 m); x 1679: 54 px (17.4 km, 622 m); x 1799: 62 px (17.9 km, 717 m); x 1919: 65 px (17.7 km, 707 m). The shore is 1.8 to 21.9 km away. Sunlit faces; haze transmittance 0.948 (green) and 0.888 (blue) at the highest point; colour #847465 (2,270 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #7b9295 (3,330 cd/m²); edges #778e92 and #81989a;
  - 1° down #7a9093 (3,250 cd/m²); edges #758c8f and #809698;
  - 3° down #768b8c (3,000 cd/m²); edges #708588 and #7b9191;
  - 8° down #6b7e7d (2,420 cd/m²); edges #687a7b and #6f8281;
  - bottom #4b5a58 (1,180 cd/m²); edges #505f5f and #536362.
- **Waves:** swell 1.2 m, period 11.3 s, 33 m long, travelling toward 138° (42° left of the view's direction); wind 6.0 m/s; mean square slope 0.0525 (rms tilt 12.9°), 0.0359 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1008 (-18 h).
- **Water:** its own light looking straight down, pale yellow-green (212 cd/m²).
- **Tide:** +4.08 m from its mean level.
- At 6.0 m/s the first whitecaps would form; the study does not model them.

### SN-2 · The Sun 4° up

*Day 20.3, nine hours before sunset. Sun 4.0° up in the west. 6,900 lux; day vision.*

**Image prompt.** Photograph from a small boat on a mirror-calm open sea in the late afternoon, eye 2 m
above the water, 28 mm lens, looking west into the low Sun, horizon 53% from the top.

The Sun is a small, blinding orange disk 14 pixels wide, 4.0 degrees above the
horizon at the centre (x = 980, y = 460). Below it a narrow, brilliant column of
reflected sunlight, no more than 70 pixels wide, runs from the horizon nearly to the bottom of
the frame. It is molten white-orange at its core and stretched into long smooth bands by a low swell.

The sky is warm amber-cream at the top (#ffe59d), deepening to soft orange at the horizon
(#f4b475), brightest above the Sun.

On the far left, a thin, faint line of distant coast 12 to 23 km away rises only 2 to 20 pixels above
the horizon, a hazy dull red-brown (#8b695b); the rest of the horizon is open sea.

The calm water beside the glitter mirrors the warm sky at 93% of its brightness:
apricot (#dbaa6f) at the horizon, darkening to bronze-olive (#7d6b46) at the
bottom. A smooth swell 1.1 m high and 71 m from crest to crest comes from the front right;
its long glassy undulations bend the reflections. There are no ripples: the wind is under 1 m/s. No
foam, no boats, no people.

**Frame data.**
- **View:** 270° (west), pitched 1° up; horizon at row 569. The frame is 52% sky, 0% land and 47% sea.
- **Exposure:** white is 4,410 cd/m².
- **Sky:**
  - top #ffe59d (3,540 cd/m²); edges #f9db96 and #fadc96 (3,210 and 3,240);
  - 15° up #ffe394 (3,530 cd/m²); edges #fbd992 and #fdda92 (3,190 and 3,220);
  - 8° up #ffd585 (3,210 cd/m²); edges #fdcc84 and #fecd84 (2,900 and 2,930);
  - 3° up #ffc27a (2,710 cd/m²); edges #f2ba7a and #f4bb7a (2,450 and 2,480);
  - 1° up #f8b877 (2,450 cd/m²); edges #eab276 and #ebb276 (2,230 and 2,250);
  - at the horizon #f4b475 (2,350 cd/m²); edges #e7ae75 and #e8af75 (2,140 and 2,160).
- **Sun:** at (980, 460), 4.03° up at 270.8°, 14 by 14 px; 40 million cd/m² (9,100 times white), #ff9d00.
- **Earth:** 59.1° up at 38°, out of frame, 62% lit; 2,000 cd/m² (0.40 times white), #ffdb8e.
- **Glitter of the Sun:** across x = 935 to 1005, y = 575 to 1005; peak 1.7 million cd/m² (390 times white) at (985, 575); 3% of the sea in the frame.
- **Coast:** skyline above the horizon line, x 0: 6 px (12.4 km, 75 m); x 120: 2 px (12.9 km, 47 m); x 240: 20 px (22.8 km, 387 m); x 360: 17 px (21.5 km, 330 m); x 480: 11 px (19.6 km, 224 m); x 600: 3 px (20.4 km, 125 m). The shore is 11.7 to 21.1 km away. Earthlit faces; haze transmittance 0.75 (green) and 0.527 (blue) at the highest point; colour #8b695b (716 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #f8b46e (2,370 cd/m²); edges #dbaa6f and #d8aa6e;
  - 1° down #ffbf65 (3,050 cd/m²); edges #d8a86d and #d5a76c;
  - 3° down #ffbb00 (18,900 cd/m²); edges #cc9f67 and #ca9f66;
  - 8° down #ffb500 (71,000 cd/m²); edges #ab8a59 and #aa8a58;
  - bottom #956c40 (766 cd/m²); edges #7d6b46 and #7d6b46.
- **Waves:** swell 1.1 m, period 16.6 s, 71 m long, travelling toward 134° (136° left of the view's direction); wind 0.6 m/s; mean square slope 0.0024 (rms tilt 2.8°), 0.0000 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1200 (+3 h).
- **Water:** its own light looking straight down, pale yellow-green (19.6 cd/m²).
- **Tide:** +0.79 m from its mean level.

### SN-3 · The Sun 4° below the horizon

*Day 21.1, nine hours after sunset. 2,870 lux; day vision.*

**Image prompt.** Photograph from a small boat on a mirror-calm sea in the early evening twilight, eye
2 m above the water, 28 mm lens, looking west to the place of sunset, horizon 58% from
the top.

The sky is a luminous gradient:
- soft amber at the top (#ffe194);
- orange through the middle;
- coral (#f1a86b) at the horizon, brightest above the centre where the Sun has gone down.

The sea is almost perfectly calm and mirrors the sky at nearly its full brightness. Just below the
horizon it is apricot-orange (#eba865), 98% as bright as the sky just above, so
the horizon line all but disappears. It darkens to warm bronze (#876c40) at the bottom of the
frame. A long, low glassy swell 1.1 m high and 71 m between crests rolls in from the front
right, drawing slow, broad bends into the mirrored glow.

On the far left a hair-thin, hazy line of distant coast rises 1 to 18 pixels, dull red-brown
(#8a6153). No ripples, no foam, no boats, no people.

**Frame data.**
- **View:** 266° (west), pitched 3° up; horizon at row 621. The frame is 57% sky, 0% land and 43% sea.
- **Exposure:** white is 2,150 cd/m².
- **Sky:**
  - top #ffe194 (1,710 cd/m²); edges #ffd68c and #ffd68c (1,530 and 1,540);
  - 15° up #ffdf88 (1,720 cd/m²); edges #ffd386 and #ffd486 (1,520 and 1,530);
  - 8° up #ffce79 (1,530 cd/m²); edges #ffc478 and #ffc478 (1,350 and 1,360);
  - 3° up #ffb770 (1,230 cd/m²); edges #f2af6f and #f3af6f (1,090 and 1,100);
  - 1° up #f6ac6c (1,080 cd/m²); edges #e6a56c and #e7a56c (964 and 970);
  - at the horizon #f1a86b (1,020 cd/m²); edges #e1a16b and #e2a26b (916 and 921).
- **Sun:** 4.0° below the horizon at 266°.
- **Earth:** 58.2° up at 37°, out of frame, 70% lit; 2,100 cd/m² (1.0 times white), #ffdb8d.
- **Coast:** skyline above the horizon line, x 0: 10 px (12.7 km, 102 m); x 120: 8 px (12.4 km, 84 m); x 240: 3 px (12.9 km, 58 m); x 360: 18 px (22.8 km, 376 m); x 480: 16 px (21.5 km, 328 m); x 600: 12 px (19.5 km, 227 m); x 720: 1 px (20.5 km, 110 m). The shore is 11.6 to 21.1 km away. Earthlit faces; haze transmittance 0.75 (green) and 0.527 (blue) at the highest point; colour #8a6153 (315 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #eba865 (1,010 cd/m²); edges #dca064 and #daa064;
  - 1° down #e7a563 (966 cd/m²); edges #da9f64 and #d79e62;
  - 3° down #db9e5e (865 cd/m²); edges #ce975e and #cc965d;
  - 8° down #bb8b51 (634 cd/m²); edges #af8450 and #ae8450;
  - bottom #876c40 (347 cd/m²); edges #8c6f43 and #8b6e43.
- **Waves:** swell 1.1 m, period 16.6 s, 71 m long, travelling toward 134° (132° left of the view's direction); wind 1.4 m/s; mean square slope 0.0024 (rms tilt 2.8°), 0.0000 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1200 (-15 h).
- **Water:** its own light looking straight down, pale yellow-green (8.12 cd/m²).
- **Tide:** +0.18 m from its mean level.

### SN-4 · The Sun 15° below the horizon

*Day 22.1, a day and a half after sunset. 466 lux; day vision.*

**Image prompt.** Photograph from a small boat on a wind-roughened sea in the long evening, eye 2 m
above the water, 28 mm lens, looking west, horizon 58% from the top.

The western sky is a broad twilight of warm amber (#ffd980 at the top), deepening to a
coral-orange band along the horizon (#ff9e56), brightest at the centre.

A thin, hazy line of distant coast 12 to 23 km away runs along the left half of the horizon, 5 to 21
pixels tall, dull red-brown (#945a44).

The sea is choppy under a 5.1 m/s wind. It reflects the sky at 57% of its
brightness as a busy texture of orange highlights and brown troughs: bronze-orange (#b48447)
below the horizon, dark brown-gold (#78572f) at the bottom. A swell 1.6 m high and
49 m from crest to crest runs away from the camera, slightly to the right; short steep wind
waves ride on it. Unbroken water, no foam, no boats, no people.

**Frame data.**
- **View:** 260° (west), pitched 3° up; horizon at row 621. The frame is 57% sky, 0% land and 43% sea.
- **Exposure:** white is 531 cd/m².
- **Sky:**
  - top #ffd980 (405 cd/m²); edges #ffca77 and #ffcb77 (348 and 350);
  - 15° up #ffdb73 (431 cd/m²); edges #ffc970 and #ffca70 (354 and 356);
  - 8° up #ffcb61 (394 cd/m²); edges #ffbb61 and #ffbb61 (323 and 324);
  - 3° up #ffb059 (307 cd/m²); edges #fca45a and #fda459 (255 and 256);
  - 1° up #ffa357 (264 cd/m²); edges #ee9957 and #ee9957 (221 and 222);
  - at the horizon #ff9e56 (247 cd/m²); edges #e89556 and #e89556 (208 and 209).
- **Sun:** 15.1° below the horizon at 260°.
- **Earth:** 57.0° up at 35°, out of frame, 80% lit; 2,300 cd/m² (4.3 times white), #ffdb8c.
- **Coast:** skyline above the horizon line, x 0: 16 px (13.3 km, 154 m); x 120: 11 px (13.2 km, 118 m); x 240: 10 px (12.7 km, 105 m); x 360: 5 px (12.6 km, 67 m); x 480: 7 px (23.2 km, 220 m); x 600: 21 px (21.9 km, 401 m); x 720: 10 px (20.2 km, 223 m); x 840: 5 px (19.9 km, 158 m). The shore is 11.6 to 20.0 km away. Earthlit faces; haze transmittance 0.75 (green) and 0.527 (blue) at the highest point; colour #945a44 (74.1 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #b48447 (141 cd/m²); edges #ac7e47 and #a27845;
  - 1° down #b18146 (136 cd/m²); edges #a97c46 and #9f7643;
  - 3° down #a97b42 (122 cd/m²); edges #a17642 and #946f3f;
  - 8° down #966c3a (93.1 cd/m²); edges #8e673a and #856238;
  - bottom #78572f (58.5 cd/m²); edges #795932 and #735530.
- **Waves:** swell 1.6 m, period 13.7 s, 49 m long, travelling toward 286° (26° right of the view's direction); wind 5.1 m/s; mean square slope 0.0463 (rms tilt 12.2°), 0.0272 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1248 (+8 h).
- **Water:** its own light looking straight down, pale yellow-green (1.32 cd/m²).
- **Tide:** −0.59 m from its mean level.
- At this wind the first whitecaps would form; the study does not model them.

### SN-5 · The Sun 35° below the horizon

*Day 24.1, three and a third days after sunset. 11.7 lux; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat at sea late in the long evening, eye 2 m
above the water, 28 mm lens, looking west-southwest, horizon 58% from the top.

A deep golden-orange twilight arch fills the western sky:
- brightest about 8 degrees up at the centre (#ffce4e);
- vivid orange along the horizon (#ffa03f);
- duller ochre at the top (#ecbd73);
- redder toward the edges.

On the left, coastal hills 18 km away rise as a dark red-brown silhouette (#8f4c38), 63 pixels
tall at the left edge. They fall toward the centre into a thin, hazy line of distant coast that ends
two-thirds of the way across.

The sea, rippled by a 4.2 m/s breeze over a 1.1 m swell running away to the right,
reflects the arch at 49% of its brightness: orange (#c9863f) below the horizon,
dark brown-orange (#805529) at the bottom. Muted to the eye; the photograph records the
colour. No stars, no foam, no boats, no people.

**Frame data.**
- **View:** 246° (west-southwest), pitched 3° up; horizon at row 621. The frame is 56% sky, 1% land and 43% sea.
- **Exposure:** white is 16.9 cd/m².
- **Sky:**
  - top #ecbd73 (9.38 cd/m²); edges #deaf6b and #deaf6b (8.01 and 7.97);
  - 15° up #ffcd66 (12.1 cd/m²); edges #f2b565 and #f1b465 (8.90 and 8.85);
  - 8° up #ffce4e (13.9 cd/m²); edges #ffb454 and #ffb354 (9.87 and 9.80);
  - 3° up #ffb640 (12.4 cd/m²); edges #ffa14a and #ffa049 (8.80 and 8.73);
  - 1° up #ffa63f (10.9 cd/m²); edges #ff9548 and #ff9447 (7.82 and 7.75);
  - at the horizon #ffa03f (10.3 cd/m²); edges #ff9047 and #ff8f47 (7.39 and 7.33).
- **Sun:** 35.2° below the horizon at 246°.
- **Earth:** 55.0° up at 31°, out of frame, 93% lit; 2,900 cd/m² (170 times white), #ffdb8b.
- **Coast:** skyline above the horizon line, x 0: 63 px (18.2 km, 719 m); x 120: 42 px (16.1 km, 448 m); x 240: 27 px (21.9 km, 461 m); x 360: 16 px (13.8 km, 171 m); x 480: 14 px (13.2 km, 146 m); x 600: 9 px (12.8 km, 102 m); x 720: 7 px (12.4 km, 85 m); x 840: 2 px (23.4 km, 157 m); x 960: 20 px (22.7 km, 413 m); x 1079: 11 px (20.6 km, 248 m); x 1199: 7 px (19.8 km, 179 m). The shore is 11.6 to 20.4 km away. Earthlit faces; haze transmittance 0.786 (green) and 0.584 (blue) at the highest point; colour #8f4c38 (1.89 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #c9863f (5.03 cd/m²); edges #b17a40 and #a1713e;
  - 1° down #c7843e (4.89 cd/m²); edges #ae773f and #9e6f3d;
  - 3° down #bf7d3a (4.40 cd/m²); edges #a5713b and #996a3a;
  - 8° down #a86d33 (3.29 cd/m²); edges #916234 and #885d32;
  - bottom #805529 (1.90 cd/m²); edges #79522c and #724f2b.
- **Waves:** swell 1.1 m, period 11.3 s, 33 m long, travelling toward 302° (56° right of the view's direction); wind 4.2 m/s; mean square slope 0.0366 (rms tilt 10.8°), 0.0170 of it from waves shorter than 1 m; nearside wave run, restart file at hour 1296 (+9 h).
- **Water:** its own light looking straight down, pale yellow-green (0.0331 cd/m²).
- **Stars:** none recorded; the eye's limit here is magnitude -3.9.
- **Tide:** −1.79 m from its mean level.

### SN-6 · The darkest hour

*Day 0.5, nine days after sunset; the Sun 52.8° below the horizon, straight ahead. 2.9 lux, mostly
from a high, nearly full Earth behind the camera; between day and night vision.*

**Image prompt.** Long-exposure night photograph from a small boat, eye 2 m above the water, 28 mm lens,
looking southeast toward the island, horizon 58% from the top.

The island fills the lower part of the view from edge to edge, 1.5 to 5 km away. It is a long, rounded
ridge of bare ground rising 51 to 110 pixels above the horizon, highest at the centre (250 to 290 m).
Beyond its right end, far hazy hills of the mainland 23 km away stand 51 to 58 pixels tall. The
island's slopes are lit softly from behind the camera by a high, nearly full Earth (out of frame), so
they show as a dim muted red-brown (#6c544f), their shapes visible.

It is low tide, −3.06 m from the mean level: a band of shore 3 m high, 2 to 3 pixels, lies
exposed along the island's foot at the waterline. The study does not calculate how wet ground looks.
Drawn darker and glistening, the band follows Earth's shores, which is an assumption here; drawn in
the island's soil colour, it follows the calculation.

The sky is a pale, nearly neutral light grey at the top (#d6d5d3), warming downward into a soft
pinkish peach glow above the hidden Sun: #ffdebd at 8 degrees and #ffc6ab at the horizon.

The sea is lightly rippled by a 4.1 m/s breeze over a short swell 0.9 m high and
23 m from crest to crest, coming toward the camera. It reflects the glow as a dim greyish
rose-taupe (#a18d83) below the horizon, darkening to dark grey-brown (#6c6059) at the
bottom. No stars show in this bright earthlit sky. No foam, no boats, no people.

**Frame data.**
- **View:** 137° (southeast), pitched 3° up; horizon at row 621. The frame is 50% sky, 8% land and 43% sea.
- **Exposure:** white is 0.963 cd/m².
- **Sky:**
  - top #d6d5d3 (0.642 cd/m²); edges #d5d4d1 and #d3d2cf (0.635 and 0.623);
  - 15° up #f2ddcb (0.723 cd/m²); edges #e1d7cc and #e0d6cb (0.662 and 0.656);
  - 8° up #ffdebd (0.780 cd/m²); edges #f7d5bf and #f7d5bf (0.682 and 0.688);
  - 3° up #ffd2b2 (0.743 cd/m²); edges #fcc9b4 and #feccb4 (0.635 and 0.649);
  - 1° up #ffcaad (0.697 cd/m²); edges #f9c2af and #fbc5af (0.595 and 0.611);
  - at the horizon #ffc6ab (0.675 cd/m²); edges #f6bfad and #f9c2ad (0.576 and 0.594).
- **Sun:** 52.8° below the horizon at 137°.
- **Earth:** 53.3° up at 22°, out of frame, 90% lit; 2,700 cd/m² (2,800 times white), #ffdb89.
- **Coast:** skyline above the horizon line, x 0: 66 px (2.8 km, 105 m); x 120: 73 px (2.8 km, 117 m); x 240: 78 px (2.8 km, 131 m); x 360: 81 px (2.9 km, 143 m); x 480: 88 px (3.8 km, 213 m); x 600: 94 px (3.8 km, 231 m); x 720: 105 px (3.3 km, 227 m); x 840: 110 px (3.4 km, 248 m); x 960: 110 px (3.6 km, 267 m); x 1079: 110 px (3.9 km, 282 m); x 1199: 107 px (4.1 km, 288 m); x 1319: 97 px (4.5 km, 281 m); x 1439: 79 px (4.8 km, 238 m); x 1559: 54 px (4.9 km, 166 m); x 1679: 51 px (23.6 km, 847 m); x 1799: 58 px (23.5 km, 921 m); x 1919: 55 px (25.2 km, 925 m). The shore is 1.4 to 20.6 km away. Earthlit faces; haze transmittance 0.955 (green) and 0.902 (blue) at the highest point; colour #6c544f (0.0982 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #a18d83 (0.274 cd/m²); edges #948984 and #9a8d86;
  - 1° down #9f8c81 (0.266 cd/m²); edges #928782 and #9a8d85;
  - 3° down #9a857b (0.242 cd/m²); edges #8b817b and #93867e;
  - 8° down #89776d (0.188 cd/m²); edges #7d736e and #837770;
  - bottom #6c6059 (0.118 cd/m²); edges #6b635f and #6e6660.
- **Waves:** swell 0.9 m, period 9.3 s, 23 m long, travelling toward 303° (166° right of the view's direction); wind 4.1 m/s; mean square slope 0.0346 (rms tilt 10.5°), 0.0157 of it from waves shorter than 1 m; nearside wave run, restart file at hour 720 (+0 h).
- **Water:** its own light looking straight down, pale green (0.00830 cd/m²).
- **Stars:** none recorded; the eye's limit here is magnitude -0.8.
- **Tide:** −3.06 m from its mean level.
- **Earth's sky:** the sky overhead and behind the camera is the blue earthlit sky of the results by
  regime (about 0.5 cd/m², chromaticity 0.27, 0.30).
- **Tide's range:** this month the water moves through 7.7 m, from +4.1 m at the noon scene.

## The South Pole–Aitken coast, by Mare Ingenii

The far side, where the Earth never rises. The camera floats 5.4 km off a coast that runs east to west
(162.875° E, 37.125° S), over water 1.4 km deep. Along the coast to the south stands a range of bare
highland mountains 700 to 1,330 m high, 10 to 21 km away. Behind it, 50 to 72 km to the southwest, rises
a higher massif of 2,300 to 3,500 m. The open sea lies to the north. No wave run reaches this sea: its
waves borrow a nearside spectrum of median steepness, turned to the local wind. The land here is the
brighter, tan highland soil.

### IC-1 · The Sun at its highest

*Day 28.4, seven and a half days before sunset. Sun 54.2° up in the north, behind the camera.
66,000 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea at midday, eye 2 m above the water, 28 mm lens,
looking south across 5 km of open water to a mountainous coast, with the high Sun behind the camera,
tilted 2 degrees down, horizon 45% from the top.

A range of bare, sunlit mountains spans the whole width of the frame, 10 to 21 km away. It rises 77
pixels above the horizon at the left edge, climbs in rounded peaks and shoulders to 174 pixels at
two-thirds of the way across (1,270 m), then falls to 115 pixels at the right edge. The mountains are
pale tan-beige highland ground (#b9a48d), softened and cooled by haze. Their front faces are
lit, with gentle shading in the folds. Vegetation is omitted for this optical study; the bare ground
is a rendering placeholder.

The sky is a very pale cyan-blue, almost white (#c8eae9 at the top), fading to pale cream-grey at
the horizon behind the peaks (#dce2ca).

The sea is a grey-teal slate (#8ea29d below the horizon), darkening to deep slate-teal
(#505e5b) at the bottom. The band of water nearest the horizon mirrors the pale mountains in
blurred streaks. A swell 1.2 m high and 33 m from crest to crest comes from the front
left; wind ripples under a 3.7 m/s breeze give it a lively, fine texture. No foam, no boats,
no people, no buildings.

**Frame data.**
- **View:** 180° (south), pitched 2° down; horizon at row 490. The frame is 32% sky, 13% land and 55% sea.
- **Exposure:** white is 12,200 cd/m².
- **Sky:**
  - top #c8eae9 (9,400 cd/m²); edges #caeae6 and #caeae6 (9,400 and 9,410);
  - 15° up #cbece7 (9,530 cd/m²); edges #c9eae7 and #c9eae7 (9,360 and 9,380);
  - 8° up #d8eddc (9,750 cd/m²); edges #d5eadc and #d5ebdc (9,540 and 9,560);
  - 3° up #dde8d2 (9,440 cd/m²); edges #d9e5d1 and #d9e6d1 (9,200 and 9,220);
  - 1° up #dde4cc (9,160 cd/m²); edges #d9e2cc and #d9e2cc (8,910 and 8,940);
  - at the horizon #dce2ca (9,020 cd/m²); edges #d8e0c9 and #d9e0c9 (8,780 and 8,800).
- **Sun:** 54.2° up at 1°, out of frame; 7.1e+02 million cd/m² (58,000 times white), #ffd874.
- **Earth:** 55.5° below the horizon at 212°, 98% lit.
- **Coast:** skyline above the horizon line, x 0: 77 px (21.1 km, 1,004 m); x 120: 97 px (17.2 km, 1,028 m); x 240: 112 px (12.6 km, 867 m); x 360: 123 px (12.9 km, 1,009 m); x 480: 132 px (13.2 km, 1,132 m); x 600: 132 px (13.6 km, 1,190 m); x 720: 140 px (13.4 km, 1,254 m); x 840: 146 px (10.5 km, 1,026 m); x 960: 145 px (10.9 km, 1,066 m); x 1079: 153 px (11.7 km, 1,201 m); x 1199: 167 px (11.9 km, 1,325 m); x 1319: 174 px (11.2 km, 1,273 m); x 1439: 168 px (11.1 km, 1,190 m); x 1559: 165 px (12.1 km, 1,247 m); x 1679: 162 px (11.8 km, 1,174 m); x 1799: 139 px (11.1 km, 913 m); x 1919: 115 px (10.8 km, 716 m). The shore is 5.5 to 7.1 km away. Sunlit faces; haze transmittance 0.863 (green) and 0.72 (blue) at the highest point; colour #b9a48d (4,700 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #8ea29d (4,130 cd/m²); edges #899d9a and #94a8a2;
  - 1° down #8da09b (4,030 cd/m²); edges #869a97 and #92a5a0;
  - 3° down #879994 (3,670 cd/m²); edges #809390 and #8c9f99;
  - 8° down #788883 (2,840 cd/m²); edges #738380 and #7c8c86;
  - bottom #505e5b (1,280 cd/m²); edges #566462 and #596865.
- **Waves:** swell 1.2 m, period 11.3 s, 33 m long, travelling toward 325° (145° right of the view's direction); wind 3.7 m/s; mean square slope 0.0282 (rms tilt 9.5°), 0.0127 of it from waves shorter than 1 m; borrowed: the Nectaris restart spectrum at hour 816, the month's median resolved slope, turned +230 degrees to the local wind, as deep water.
- **Water:** its own light looking straight down, pale yellow-green (187 cd/m²).
- **Tide:** −0.53 m from its mean level.

### IC-2 · The Sun 4° up

*Day 5.9, ten hours before sunset. Sun 3.9° up in the west. 6,900 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the late afternoon, eye 2 m above the water,
28 mm lens, looking west along a mountainous coast into the low Sun, horizon 53% from
the top.

The Sun is a small, blinding orange disk 14 pixels wide, 3.9 degrees up just right
of centre (x = 1061, y = 463). Below it a broad, sparkling orange glitter path, up to
280 pixels wide, runs from the horizon to the bottom of the frame over a rippled sea.

On the left the coast runs out toward the horizon:
- a near coastal ridge 20 km away stands as a dark, hazy red-brown silhouette (#8b6d5d),
  99 pixels tall at the left edge;
- beyond it, from an eighth to half of the way across, a far, much higher mountain massif 50 to 72 km
  away rises 22 to 63 pixels above the horizon;
- the massif is so hazed that it is only a faint, pale, warm silhouette barely darker than the sky
  behind it.

The right half of the horizon is open sea.

The sky is warm amber-cream at the top (#ffe59d), deepening to soft orange at the horizon
(#f4b475). The sea beside the glitter is warm tan-gold (#a88e60), darkening to dark
olive-bronze (#6f5f41) in the lower corners. A swell 1.2 m high and 33 m
between crests comes from the front left, and a 3.3 m/s breeze ripples it. No foam, no boats,
no people.

**Frame data.**
- **View:** 268° (west), pitched 1° up; horizon at row 569. The frame is 50% sky, 3% land and 47% sea.
- **Exposure:** white is 4,380 cd/m².
- **Sky:**
  - top #ffe59d (3,510 cd/m²); edges #f7d995 and #fddd96 (3,130 and 3,270);
  - 15° up #ffe294 (3,500 cd/m²); edges #f9d792 and #ffdb93 (3,110 and 3,250);
  - 8° up #ffd585 (3,180 cd/m²); edges #facb84 and #ffce84 (2,830 and 2,960);
  - 3° up #ffc17a (2,680 cd/m²); edges #f0b97a and #f6bc7a (2,390 and 2,490);
  - 1° up #f8b876 (2,430 cd/m²); edges #e8b076 and #eeb376 (2,170 and 2,260);
  - at the horizon #f4b475 (2,330 cd/m²); edges #e4ad75 and #eab075 (2,080 and 2,170).
- **Sun:** at (1061, 463), 3.93° up at 271.8°, 14 by 14 px; 39 million cd/m² (8,900 times white), #ff9d00.
- **Earth:** 56.3° below the horizon at 201°, 44% lit.
- **Glitter of the Sun:** across x = 945 to 1225, y = 575 to 1075; peak 120,000 cd/m² (28 times white) at (1065, 575); 10% of the sea in the frame.
- **Coast:** skyline above the horizon line, x 0: 99 px (20.4 km, 1,226 m); x 120: 71 px (20.3 km, 928 m); x 240: 55 px (50.2 km, 2,308 m); x 360: 59 px (52.6 km, 2,628 m); x 480: 63 px (60.3 km, 3,371 m); x 600: 55 px (65.0 km, 3,438 m); x 720: 35 px (65.2 km, 2,632 m); x 840: 25 px (71.5 km, 2,533 m); x 960: 22 px (70.5 km, 2,368 m). The shore is 11.8 to 35.2 km away. Backlit; haze transmittance 0.76 (green) and 0.543 (blue) at the highest point; colour #8b6d5d (756 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #ba9a63 (1,500 cd/m²); edges #a88e60 and #ba9a65;
  - 1° down #b89863 (1,470 cd/m²); edges #a68c5f and #b79763;
  - 3° down #af905d (1,310 cd/m²); edges #9e865a and #af905f;
  - 8° down #d18a4d (1,410 cd/m²); edges #8b764f and #987e53;
  - bottom #e17a30 (1,320 cd/m²); edges #6f5f41 and #776543.
- **Waves:** swell 1.2 m, period 11.3 s, 33 m long, travelling toward 43° (135° right of the view's direction); wind 3.3 m/s; mean square slope 0.026 (rms tilt 9.2°), 0.0106 of it from waves shorter than 1 m; borrowed: the Nectaris restart spectrum at hour 816, the month's median resolved slope, turned +153 degrees to the local wind, as deep water.
- **Water:** its own light looking straight down, pale yellow-green (19.4 cd/m²).
- **Tide:** +0.56 m from its mean level.

### IC-3 · The Sun 4° below the horizon

*Day 6.8, ten hours after sunset. 2,850 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the early evening twilight, eye 2 m above the
water, 28 mm lens, looking west along a mountainous coast, horizon 58% from the top.

The sky glows:
- soft amber at the top (#ffe194);
- orange through the middle;
- coral (#f1a86b) along the horizon, brightest just right of centre where the Sun has set.

On the left the coast's near ridge, 15 to 21 km away, stands as a dark, hazy red-brown silhouette
(#7d5d4e), 60 to 105 pixels tall. Beyond it, across the middle of the frame, the far high massif
50 to 72 km away shows as a pale, ghostly warm silhouette 25 to 58 pixels tall, little darker than the
sky. The right third of the horizon is open sea.

The sea, gently rippled by a 2.1 m/s breeze over a 1.2 m swell, reflects the glow at
81% of its brightness: apricot-orange (#d09e5e) below the horizon, warm bronze
(#81663e) at the bottom. The water just below the horizon on the left mirrors the dark ridge
in soft, broken streaks. No foam, no boats, no people.

**Frame data.**
- **View:** 262° (west), pitched 3° up; horizon at row 621. The frame is 54% sky, 3% land and 43% sea.
- **Exposure:** white is 2,120 cd/m².
- **Sky:**
  - top #ffe194 (1,690 cd/m²); edges #fdd38b and #ffd88d (1,480 and 1,560);
  - 15° up #ffde88 (1,700 cd/m²); edges #ffd185 and #ffd686 (1,460 and 1,540);
  - 8° up #ffce78 (1,510 cd/m²); edges #ffc177 and #ffc678 (1,300 and 1,370);
  - 3° up #ffb76f (1,210 cd/m²); edges #efad6f and #f6b16f (1,050 and 1,110);
  - 1° up #f6ac6c (1,060 cd/m²); edges #e3a36c and #eaa66c (930 and 977);
  - at the horizon #f1a86b (1,010 cd/m²); edges #dda06b and #e4a36b (882 and 927).
- **Sun:** 4.1° below the horizon at 266°.
- **Earth:** 55.6° below the horizon at 199°, 36% lit.
- **Coast:** skyline above the horizon line, x 0: 105 px (17.9 km, 1,116 m); x 120: 97 px (20.6 km, 1,252 m); x 240: 87 px (20.1 km, 1,137 m); x 360: 60 px (14.6 km, 580 m); x 480: 59 px (48.8 km, 2,431 m); x 600: 57 px (61.0 km, 3,216 m); x 720: 58 px (64.2 km, 3,523 m); x 840: 41 px (65.0 km, 2,889 m); x 960: 29 px (68.0 km, 2,523 m); x 1079: 25 px (69.2 km, 2,430 m); x 1199: 7 px (72.0 km, 1,711 m). The shore is 11.8 to 36.2 km away. Backlit; haze transmittance 0.821 (green) and 0.645 (blue) at the highest point; colour #7d5d4e (271 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #d09e5e (820 cd/m²); edges #b8905b and #c8995e;
  - 1° down #cd9b5d (792 cd/m²); edges #b58e59 and #c5975d;
  - 3° down #c39358 (707 cd/m²); edges #ad8755 and #bb8f58;
  - 8° down #a8814d (521 cd/m²); edges #96764b and #a27d4d;
  - bottom #81663e (306 cd/m²); edges #7e6440 and #866a41.
- **Waves:** swell 1.2 m, period 11.3 s, 33 m long, travelling toward 19° (117° right of the view's direction); wind 2.1 m/s; mean square slope 0.0191 (rms tilt 7.9°), 0.0037 of it from waves shorter than 1 m; borrowed: the Nectaris restart spectrum at hour 816, the month's median resolved slope, turned +176 degrees to the local wind, as deep water.
- **Water:** its own light looking straight down, pale yellow-green (8.00 cd/m²).
- **Tide:** +0.79 m from its mean level.

### IC-4 · The Sun 15° below the horizon

*Day 7.9, a day and a half after sunset. 456 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the long evening, eye 2 m above the water,
28 mm lens, looking west-southwest along a mountainous coast, horizon 58% from the top.

A broad amber twilight fills the sky: warm amber at the top (#ffd980), deepening to a
coral-orange band along the horizon (#ff9e55).

The near coastal mountains, 15 to 21 km away, stand on the left as a dark red-brown silhouette
(#7b513f), 77 to 114 pixels tall. They lead into the pale, hazed far massif across the middle,
50 to 72 km away and 24 to 61 pixels tall, its outline soft against the glow. The right quarter is
open sea.

The sea, rippled by a 3.6 m/s breeze over a 1.2 m swell crossing from left to right,
reflects the sky at 78% of its brightness: amber-orange (#d5964c) below the
horizon, brown-gold (#846032) at the bottom. The mountains' dark reflection blurs into the
water below them on the left. No foam, no boats, no people.

**Frame data.**
- **View:** 255° (west-southwest), pitched 3° up; horizon at row 621. The frame is 53% sky, 4% land and 43% sea.
- **Exposure:** white is 520 cd/m².
- **Sky:**
  - top #ffd980 (395 cd/m²); edges #ffc876 and #ffcc77 (334 and 348);
  - 15° up #ffdb72 (421 cd/m²); edges #ffc76f and #ffcb70 (339 and 354);
  - 8° up #ffcb60 (385 cd/m²); edges #ffb961 and #ffbd61 (309 and 323);
  - 3° up #ffb058 (301 cd/m²); edges #faa359 and #ffa559 (245 and 255);
  - 1° up #ffa356 (258 cd/m²); edges #ec9857 and #f19a57 (212 and 221);
  - at the horizon #ff9e55 (242 cd/m²); edges #e59456 and #eb9656 (200 and 207).
- **Sun:** 15.3° below the horizon at 257°.
- **Earth:** 54.3° below the horizon at 198°, 25% lit.
- **Coast:** skyline above the horizon line, x 0: 114 px (14.6 km, 975 m); x 120: 106 px (15.4 km, 997 m); x 240: 98 px (17.7 km, 1,104 m); x 360: 92 px (20.6 km, 1,267 m); x 480: 77 px (20.2 km, 1,075 m); x 600: 51 px (60.8 km, 2,957 m); x 720: 55 px (52.2 km, 2,590 m); x 840: 61 px (60.3 km, 3,371 m); x 960: 52 px (65.4 km, 3,391 m); x 1079: 33 px (66.6 km, 2,635 m); x 1199: 27 px (71.5 km, 2,625 m); x 1319: 24 px (72.0 km, 2,483 m). The shore is 7.8 to 35.5 km away. Backlit; haze transmittance 0.827 (green) and 0.654 (blue) at the highest point; colour #7b513f (54.2 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #d5964c (189 cd/m²); edges #b88549 and #c08a4a;
  - 1° down #d1934a (182 cd/m²); edges #b68449 and #bd8849;
  - 3° down #c78c47 (162 cd/m²); edges #ac7c44 and #b48145;
  - 8° down #ad7a3e (119 cd/m²); edges #976d3c and #9d703d;
  - bottom #846032 (70.1 cd/m²); edges #7e5c33 and #835f34.
- **Waves:** swell 1.2 m, period 11.3 s, 33 m long, travelling toward 356° (102° right of the view's direction); wind 3.6 m/s; mean square slope 0.0275 (rms tilt 9.4°), 0.0120 of it from waves shorter than 1 m; borrowed: the Nectaris restart spectrum at hour 816, the month's median resolved slope, turned +199 degrees to the local wind, as deep water.
- **Water:** its own light looking straight down, pale yellow-green (1.28 cd/m²).
- **Tide:** +1.10 m from its mean level.

### IC-5 · The Sun 35° below the horizon

*Day 10.2, nearly four days after sunset. 8.9 lux; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat at sea late in the long evening, eye 2 m
above the water, 28 mm lens, looking southwest along the mountainous coast, horizon 58%
from the top.

A deep twilight arch glows over the mountains:
- brightest amber-gold about 8 degrees up at the centre (#ffcc3c);
- vivid orange along the horizon (#ff9d2b);
- ochre at the top (#eebb68);
- redder toward the edges.

The coastal mountains fill the whole width of the horizon as a dark red-brown silhouette
(#7e432c):
- tallest at the left edge, 178 pixels (1,240 m at 12 km);
- stepping down through rounded shoulders to about 90 pixels across the centre;
- passing on the right into the far high massif, 51 to 72 km away and 27 to 61 pixels tall, paler with
  distance.

The sea is nearly calm, with no wind ripples under a 1.7 m/s breeze, only a smooth 1.2 m
swell coming toward the camera. It mirrors the glow at 62% of its brightness: orange
(#e59235) below the horizon, with the mountains' dark shapes reflected in the nearest band of
water, darkening to brown-orange (#875924) at the bottom. Muted to the eye; the photograph
records the colour. No stars, no foam, no boats, no people.

**Frame data.**
- **View:** 236° (southwest), pitched 3° up; horizon at row 621. The frame is 50% sky, 7% land and 43% sea.
- **Exposure:** white is 16.7 cd/m².
- **Sky:**
  - top #eebb68 (9.13 cd/m²); edges #dfad5f and #dfac5f (7.71 and 7.69);
  - 15° up #ffcc59 (11.9 cd/m²); edges #f3b258 and #f2b258 (8.61 and 8.58);
  - 8° up #ffcc3c (13.7 cd/m²); edges #ffb144 and #ffb144 (9.58 and 9.55);
  - 3° up #ffb42b (12.2 cd/m²); edges #ff9d39 and #ff9d39 (8.51 and 8.47);
  - 1° up #ffa32a (10.7 cd/m²); edges #ff9137 and #ff9037 (7.52 and 7.49);
  - at the horizon #ff9d2b (10.0 cd/m²); edges #ff8b37 and #ff8b37 (7.09 and 7.06).
- **Sun:** 35.2° below the horizon at 236°.
- **Earth:** 50.8° below the horizon at 197°, 8% lit.
- **Coast:** skyline above the horizon line, x 0: 178 px (12.1 km, 1,238 m); x 120: 156 px (11.5 km, 1,067 m); x 240: 121 px (11.0 km, 819 m); x 360: 97 px (10.7 km, 654 m); x 480: 89 px (14.0 km, 824 m); x 600: 99 px (14.6 km, 974 m); x 720: 93 px (14.9 km, 956 m); x 840: 89 px (18.1 km, 1,133 m); x 960: 81 px (20.2 km, 1,181 m); x 1079: 55 px (14.6 km, 570 m); x 1199: 56 px (51.1 km, 2,544 m); x 1319: 61 px (61.5 km, 3,420 m); x 1439: 58 px (64.5 km, 3,482 m); x 1559: 43 px (64.8 km, 2,846 m); x 1679: 33 px (68.0 km, 2,583 m); x 1799: 33 px (71.8 km, 2,750 m); x 1919: 27 px (72.0 km, 2,483 m). The shore is 6.0 to 35.5 km away. Backlit; haze transmittance 0.853 (green) and 0.702 (blue) at the highest point; colour #7e432c (1.44 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #e59235 (6.23 cd/m²); edges #be7d37 and #ce8437;
  - 1° down #e28f34 (6.02 cd/m²); edges #ba7b36 and #cc8236;
  - 3° down #d88731 (5.38 cd/m²); edges #b27533 and #c27b33;
  - 8° down #ba752b (3.87 cd/m²); edges #99652c and #a46a2d;
  - bottom #875924 (2.06 cd/m²); edges #7c5326 and #825727.
- **Waves:** swell 1.2 m, period 11.3 s, 33 m long, travelling toward 32° (156° right of the view's direction); wind 1.7 m/s; mean square slope 0.0154 (rms tilt 7.1°), 0.0000 of it from waves shorter than 1 m; borrowed: the Nectaris restart spectrum at hour 816, the month's median resolved slope, turned +163 degrees to the local wind, as deep water.
- **Water:** its own light looking straight down, pale yellow-green (0.0250 cd/m²).
- **Stars:** none recorded; the eye's limit here is magnitude -3.8.
- **Tide:** +1.49 m from its mean level.

### IC-6 · The darkest hour: the far side's twilit night

*Day 13.7, seven and a third days after sunset; the Sun 51.8° below the southern horizon, straight
ahead. 0.46 lux, all from the twilight; between day and night vision.*

**Image prompt.** Long-exposure night photograph from a small boat at sea, eye 2 m above the water,
28 mm lens, looking south toward a mountainous coast, horizon 58% from the top.

All night a golden twilight arch stands over the southern mountains, lit by a Sun more than 50 degrees
below the horizon:
- warm amber-gold, brightest about 8 degrees up (#ffcf52);
- orange near the mountain tops (#ffbd41 at 3 degrees);
- ochre-orange at the top of the frame (#f0bf74).

The mountain range, 10 to 21 km away, spans the whole frame as a dark silhouette: 76 pixels tall at the
left edge, rising to 172 pixels two-thirds of the way across (1,270 m), and 115 pixels at the right
edge. The haze tints it a deep rust-brown (#9b4f34) in this long exposure. It hides the arch's
lowest, reddest part.

The sea is lightly rippled by a 2.1 m/s breeze over a 1.2 m swell running away to the
right. It reflects the arch as a warm orange sheen (#e19442 below the horizon at the centre),
with the mountains' dark outline mirrored, blurred, in the band of water just below the horizon. It
darkens to brown-orange (#865b2c) at the bottom.

The Earth never rises here. No stars show against this glowing sky; the brightest in view, Achernar, is
too dim through the air. To the eye: a dim, nearly colourless twilight with a pale glow over black
mountains. No foam, no boats, no people.

**Frame data.**
- **View:** 180° (south), pitched 3° up; horizon at row 621. The frame is 45% sky, 13% land and 43% sea.
- **Exposure:** white is 0.623 cd/m².
- **Sky:**
  - top #f0bf74 (0.355 cd/m²); edges #e7b46d and #e8b56d (0.317 and 0.319);
  - 15° up #ffcc69 (0.439 cd/m²); edges #f8b967 and #f9b967 (0.346 and 0.348);
  - 8° up #ffcf52 (0.510 cd/m²); edges #ffb856 and #ffb856 (0.379 and 0.382);
  - 3° up #ffbd41 (0.486 cd/m²); edges #ffa74a and #ffa84a (0.349 and 0.352);
  - 1° up #ffaf3e (0.445 cd/m²); edges #ff9c47 and #ff9d47 (0.318 and 0.321);
  - at the horizon #ffa93e (0.425 cd/m²); edges #ff9746 and #ff9846 (0.304 and 0.306).
- **Sun:** 51.8° below the horizon at 180°.
- **Earth:** 45.0° below the horizon at 202°, 2% lit.
- **Coast:** skyline above the horizon line, x 0: 76 px (21.1 km, 1,001 m); x 120: 96 px (17.2 km, 1,028 m); x 240: 111 px (12.6 km, 867 m); x 360: 123 px (12.9 km, 1,008 m); x 480: 131 px (13.2 km, 1,131 m); x 600: 131 px (13.6 km, 1,188 m); x 720: 139 px (13.4 km, 1,253 m); x 840: 145 px (10.5 km, 1,024 m); x 960: 144 px (10.9 km, 1,064 m); x 1079: 151 px (11.7 km, 1,199 m); x 1199: 165 px (11.9 km, 1,323 m); x 1319: 172 px (11.2 km, 1,272 m); x 1439: 167 px (11.1 km, 1,189 m); x 1559: 164 px (12.1 km, 1,246 m); x 1679: 161 px (11.8 km, 1,173 m); x 1799: 138 px (11.1 km, 913 m); x 1919: 114 px (10.8 km, 715 m). The shore is 5.5 to 7.1 km away. Backlit; haze transmittance 0.863 (green) and 0.72 (blue) at the highest point; colour #9b4f34 (0.0792 cd/m²).
- **Sea** (reflecting the sky through this scene's slopes; open water assumed behind):
  - below the horizon #e19442 (0.233 cd/m²); edges #cf8a44 and #bb8042;
  - 1° down #df9141 (0.227 cd/m²); edges #cc8743 and #b87e41;
  - 3° down #d58a3d (0.204 cd/m²); edges #c3803f and #b1783d;
  - 8° down #b87736 (0.148 cd/m²); edges #a76f37 and #9a6835;
  - bottom #865b2c (0.0790 cd/m²); edges #855b2f and #7e572d.
- **Waves:** swell 1.2 m, period 11.3 s, 33 m long, travelling toward 232° (52° right of the view's direction); wind 2.1 m/s; mean square slope 0.0189 (rms tilt 7.8°), 0.0035 of it from waves shorter than 1 m; borrowed: the Nectaris restart spectrum at hour 816, the month's median resolved slope, turned +323 degrees to the local wind, as deep water.
- **Water:** its own light looking straight down, pale yellow (0.00130 cd/m²).
- **Stars:** none recorded; the eye's limit here is magnitude -0.4.
- **Tide:** +0.93 m from its mean level.

## Boundaries

- **Sky and sea colours** come from the model of the results by regime, evaluated for each frame's
  directions with its own slopes and camera. The sea in it is statistical: no individual waves, and open
  water in every direction. Where land rises behind the sea, the scene text says how the nearest water
  mirrors it, without numbers.
- **The coast** is the LOLA altimetry grid, about 118 m between points, at the moment's tide. Its colour
  is a placeholder: bare Apollo soil, lit by the clear sky and the disks on a slope facing the camera,
  and hazed by extinction at the ground toward the colour of the horizon sky. Beaches, cliffs, surf, wet
  ground and plants are not modelled. The vegetation omission is a boundary of this optical study,
  not a conclusion about which coasts support life. The Smythii frames show no land: the headland is
  behind the camera.
- **Waves:** significant heights, periods and directions come from the wave spectra.
  - At the nearside shores these are the restart files nearest the moment, up to 21 hours away.
  - At Smythii they are the hourly west-face spectra.
  - At Ingenii they are a borrowed nearside spectrum.
  - Short waves follow the GCM's stress at the hour. Whitecaps and foam are not modelled.
- **Stars:** the eye's limit follows Schaefer (1990). A star counts as recorded when its point image,
  spread over about four pixels where it falls, reaches a fifth of the sky's brightness around it.
  Planets are not placed.
- **Sizes in pixels** follow the rectilinear lens's projection where each object falls.
- **Approximate reading:** the Earth and the Sun are uniform disks, without the Earth's clouds and
  continents and without limb darkening. Refraction is left out, about 0.14 degrees at the horizon in
  this air.

```sh
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=2 <numba environment>/bin/python -m research.studies.sea_appearance.scenes
python -m research.studies.sea_appearance.scenes --text     # rewrite scenes.md from the product only
```

The runner reads the solved sky and the coast grids on the research drive
(`research/runs/optical_comfort/spherical/`, `research/runs/sea_appearance/terrain/`);
`python -m geography.coast_terrain --download` restores the grids' LOLA rows.
