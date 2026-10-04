<!-- The words of scenes.md. Names in braces, {SCENE:name}, are filled from results/scenes.json by
     scenes.py, which writes scenes.md; edit the words here and the numbers in the product. -->
# The four coasts, scene by scene

Twenty-four views of the Open Moon's seas: four coasts at the six moments of the
[results by regime](README.md#the-seas-by-regime), from the Sun at its highest to the
darkest hour of the night. Each scene is the frame the
[sea-scene renderer](../../../visualization/reference-renderer/seas/) will draw. They are described here
in words and numbers, so that they can be pictured, or illustrated with an image
generator, before the full renderings are made.

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

## Reading the scenes

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
  15 degrees up and deep orange near the horizon. Its disk is 0.53 degrees across, about {ALL:sun_px}
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
- **The Earth** is {ALL:earth_deg_min} to {ALL:earth_deg_max} degrees across as its distance changes through
  the month, {ALL:earth_px_min} to {ALL:earth_px_max} pixels near the centre of the frame, nearly four times the
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
  and a brighter tan highland soil at the Ingenii coast. There are no plants, buildings, roads, boats or
  people.
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

*Day 18.3 of the month (from 7 February 2038), seven and a half days before sunset. Sun {WP-1:sun_el}°
up in the south. 72,500 lux on the ground; day vision.*

**Image prompt.** Photograph from a small boat on a wide sea at midday under a cloudless sky, eye 2 m
above the water, 28 mm lens, looking south-southwest, camera tilted 3 degrees down so the horizon sits
{WP-1:horizon_pct}% from the top of the frame. The high Sun is out of frame, ahead and to the left, so
the light falls from high in front: the water shows no glitter in view.

The sky is a very pale cyan-blue, almost white ({WP-1:sky_top} at the top). It is a little brighter there
than near the horizon, where it turns a pale grey-green-cream ({WP-1:sky_hz}), and it is evenly bright
from side to side.

Along the whole horizon runs a low, faint range of hills and crater rims 16 to 42 km away, 25 to 63
pixels tall, highest at three-quarters of the way across. The hills are a hazy, desaturated blue-grey
({WP-1:land_hex}), only slightly darker than the sky above them, sunlit but softened by the long path
through the air.

The sea is a grey-teal slate ({WP-1:sea_below} just below the horizon), much darker than the sky. It
darkens smoothly to {WP-1:sea_bottom} at the bottom of the frame, with a faint olive tint in the nearest
water. A long ocean swell of {WP-1:hs} m significant height (crest to trough), {WP-1:wl} m between
crests, runs away from the camera toward the hills. Its crests lie parallel to the horizon as broad,
low, rounded ridges, so the near water rises and falls in long gentle slopes. Wind ripples from a few
centimetres to a few decimetres long, coarser than ripples on Earth, cover everything under a
{WP-1:wind} m/s breeze from behind the camera. They break the reflected sky into a fine, high-contrast
mosaic of pale and slate facets; the waves show strongly. Unbroken water, no foam. No boats, no people,
no plants, no buildings.

**Frame data.**
{WP-1:frame_data}

### WP-2 · The Sun 4° up

*Day 25.3, nine hours before sunset. Sun {WP-2:sun_el}° up in the west. 6,700 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the late afternoon, eye 2 m above the water,
28 mm lens, looking west into the low Sun, horizon {WP-2:horizon_pct}% from the top.

The Sun is a small, intensely bright orange disk, {WP-2:sun_px} pixels wide ({WP-2:sun_hex} at its own
exposure, here blown out to white with a warm bloom). It stands {WP-2:sun_el} degrees above the horizon
just right of centre (x = {WP-2:sun_x}, y = {WP-2:sun_y}), above a dark, hazy line of coastal hills.
Beneath it a broad, dazzling orange glitter path, up to {WP-2:glit_w} pixels wide, runs from the horizon
down to the bottom of the frame. It is made of countless sparkling highlights on a rough, rippled sea,
brightest at the horizon.

The whole sky is warm:
- pale amber-cream at the top ({WP-2:sky_top}), slightly brighter above the Sun;
- deepening through {WP-2:sky_8} at 8 degrees and {WP-2:sky_3} at 3 degrees;
- to warm orange at the horizon ({WP-2:sky_hz}).

There is no halo: the glow is broad and smooth.

The coast across the whole width of the frame, 4 to 8 km away with hills 15 to 30 km behind it, is a
low backlit silhouette 17 to 61 pixels tall, highest left of centre. It is dull red-brown
({WP-2:land_hex}), softened by haze.

Away from the glitter the sea is a warm tan-bronze ({WP-2:sea_below_l} below the horizon), darkening to
dark olive-brown ({WP-2:sea_bottom_l}) in the lower corners. A long, low swell {WP-2:hs} m high and
{WP-2:wl} m from crest to crest travels away to the left. Its crests run diagonally, and the glitter path
wobbles and breaks across them. Wind ripples, coarser than on Earth, cover everything under a
{WP-2:wind} m/s breeze. Unbroken water, no foam, no boats, no people.

**Frame data.**
{WP-2:frame_data}

### WP-3 · The Sun 4° below the horizon

*Day 26.1, nine hours after sunset. 2,770 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the early evening twilight, eye 2 m above the
water, 28 mm lens, looking west toward where the Sun has set, horizon {WP-3:horizon_pct}% from the top.
The Sun itself is {WP-3:sun_depth} degrees below the horizon.

The sky is a luminous warm gradient, brightest above the place of sunset just right of centre:
- soft amber at the top ({WP-3:sky_top});
- light orange ({WP-3:sky_8}) at 8 degrees and {WP-3:sky_3} at 3 degrees;
- coral-orange at the horizon ({WP-3:sky_hz}).

The amber fills the frame from the horizon to its top, {WP-3:sky_hz_lum} to {WP-3:sky_top_lum} cd/m².

A low, hazy silhouette of coastal hills, dull red-brown ({WP-3:land_hex}), lies along the horizon across
the frame. It is 14 to 62 pixels tall, highest in the left half (630 m hills at 15 km), lower and fainter
on the right.

The sea is a choppy, wind-roughened surface that reflects the amber sky diffusely, at
{WP-3:sea_ratio}% of its brightness. It is warm tan-orange ({WP-3:sea_below} just below the horizon),
darkening to deep bronze ({WP-3:sea_bottom}) at the bottom. A long swell {WP-3:hs} m high and {WP-3:wl} m
from crest to crest travels away to the left under a {WP-3:wind} m/s wind. Ripples cover every slope and
break the reflection into a fine texture of orange highlights and darker brown facets, with no distinct
glitter path. Unbroken water, no foam, no boats, no people.

**Frame data.**
{WP-3:frame_data}
- At this wind a few whitecaps would form on Earth; the study does not model them.

### WP-4 · The Sun 15° below the horizon

*Day 27.1, a day and a half after sunset. 479 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the long evening, eye 2 m above the water,
28 mm lens, looking east, away from the sunset, horizon {WP-4:horizon_pct}% from the top.

A full Earth hangs in the sky just left of centre, {WP-4:earth_el} degrees above the horizon (centre
x = {WP-4:earth_x}, y = {WP-4:earth_y}). It is a round disk {WP-4:earth_px} pixels across, bright
amber-gold in colour ({WP-4:earth_hex} at its own exposure), here overexposed to a warm white with a soft
bloom.

The sky is a soft, even, pale golden cream: {WP-4:sky_top} at the top and {WP-4:sky_hz} just above the
horizon. It is a fifth duller and greyer toward the horizon, where the Moon's own shadow rises into the
air.

A row of low, flat islets 1.8 to 2.5 km away runs from the left edge to about nine-tenths of the way
across. They are 12 to 60 m high and 8 to 43 pixels tall, highest on the left, and their bases sit on
the water just under the horizon. They are dark muted brown silhouettes ({WP-4:land_hex}), lit from
behind by the Earth.

The sea is a muted grey-gold, {WP-4:sea_below_l} to {WP-4:sea_below} just below the horizon, slightly
warmer and brighter in a broad soft column below the Earth. It darkens to olive-grey
({WP-4:sea_bottom}) at the bottom. There is no sharp glitter path: the earthlight is faint against the
reflected sky. A long swell {WP-4:hs} m high, {WP-4:wl} m from crest to crest, comes toward the camera
from the front left, its crests running diagonally. Wind ripples under a {WP-4:wind} m/s breeze give
the water a fine, matte texture. Unbroken water, no foam, no boats, no people.

**Frame data.**
{WP-4:frame_data}

### WP-5 · The Sun 35° below the horizon

*Day 29.0, three and a third days after sunset. 9.6 lux; between day and night vision.*

**Image prompt.** Photograph from a small boat at sea late in the long evening, eye 2 m above the water,
28 mm lens, looking west-northwest, horizon {WP-5:horizon_pct}% from the top; a long exposure.

A great twilight arch glows over the western horizon. It is a broad dome of deep gold and orange
centred a little right of the middle of the frame:
- brightest amber-gold about 8 degrees up ({WP-5:sky_8});
- vivid orange at the horizon ({WP-5:sky_hz});
- duller ochre-orange at the top of the frame ({WP-5:sky_top});
- redder and darker toward the left edge ({WP-5:sky_hz_l} low down).

The Sun is 35 degrees below the horizon, and this light comes from the very high air it still reaches.

A hazy silhouette of coastal hills, dark red-brown ({WP-5:land_hex}), runs along the horizon:
- tallest at the far left, 66 pixels (640 m hills at 16 km);
- sloping down to a thin, faint line of distant ranges only 6 to 15 pixels tall on the right;
- with a few low islands near the right edge.

The sea is roughened by a swell {WP-5:hs} m high and {WP-5:wl} m between crests, running away to the
left, and by wind ripples under a {WP-5:wind} m/s breeze. It reflects the arch at about half its
brightness ({WP-5:sea_ratio}% just below the horizon): warm orange ({WP-5:sea_below}) at the centre,
duller toward the edges ({WP-5:sea_below_l} on the left), and deep brown-orange ({WP-5:sea_bottom}) at
the bottom. The tilted facets break the reflection into a grainy orange texture. To the eye at this
light level colours are muted; the photograph records them. No stars, no foam, no boats, no people.

**Frame data.**
{WP-5:frame_data}

### WP-6 · The darkest hour

*Day 4.5, eight days after sunset; the Sun 60° below the horizon. 0.48 lux, from the Earth and the
twilight; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat on a glassy night sea, eye 2 m above the
water, 28 mm lens, looking east toward the Earth, horizon {WP-6:horizon_pct}% from the top.

The Earth hangs near the top centre, {WP-6:earth_el} degrees up (centre x = {WP-6:earth_x},
y = {WP-6:earth_y}). It is a gibbous disk {WP-6:earth_px} pixels across, {WP-6:earth_lit}% lit, its
bright side toward the {WP-6:earth_side}. It is golden-amber in colour ({WP-6:earth_hex} at its own
exposure) but here overexposed to near white with a bloom, the unlit part invisible against the sky.

Straight down from it a vivid golden-yellow glitter path, up to {WP-6:glit_w} pixels wide, runs from
the horizon to the bottom of the frame ({WP-6:sea_bottom} at its brightest). Its edge is soft and its
interior a ladder of long smooth highlights, because the sea is glassy.

The sky is pale and warm:
- cream-yellow overhead ({WP-6:sky_top});
- a soft peach-orange glow along the left horizon ({WP-6:sky_hz_l}), on the side toward the north,
  where the twilight of the deep-set Sun lingers;
- a greyer, duller yellow-beige toward the right ({WP-6:sky_hz_r}).

A row of low flat islets 2 km away, 20 to 60 m high and 15 to 43 pixels tall, stands as near-black
silhouettes ({WP-6:land_hex}) from the left edge to four-fifths of the way across.

The calm sea mirrors the sky, peach-tan ({WP-6:sea_below_l}) on the left and grey-beige
({WP-6:sea_below_r}) on the right, darkening downward. It moves in a long, low, glassy swell
{WP-6:hs} m high and {WP-6:wl} m from crest to crest, coming from the front left, with no wind ripples.

To the eye this is a dim, almost colourless night with the bright Earth and its path. No stars show
against this bright sky. No foam, no boats, no people.

**Frame data.**
{WP-6:frame_data}

## The eastern Smythii headland

The camera floats about 1 km off the west face of a narrow headland in the eastern Smythii basin
(93.22° E, 2.36° N), over water about 270 m deep. Its waves are those of the shore history's west-face
station, 1 km inshore in 22 m of water. It looks west across the open sea; the headland is behind the
camera. No land enters these frames. The Earth stands low in the west, rising and setting with the
libration: below the horizon for the first five scenes, just risen at the darkest hour.

### SH-1 · The Sun at its highest

*Day 4.6, seven and a half days before sunset. Sun {SH-1:sun_el}° up, almost overhead. 87,200 lux; day
vision.*

**Image prompt.** Photograph from a small boat on the open sea at midday, eye 2 m above the water, 28 mm
lens, looking west across empty water to a perfectly straight horizon, camera tilted 3 degrees down,
horizon {SH-1:horizon_pct}% from the top. The Sun is almost directly overhead, out of frame, its glitter
beneath the boat out of view.

The sky is a very pale, luminous cyan-blue ({SH-1:sky_top} at the top), almost white, fading to a pale
grey-green-cream at the horizon ({SH-1:sky_hz}). It is evenly bright across the frame: a clean,
washed-out noon sky.

The sea is a cool grey-teal ({SH-1:sea_below} just below the horizon), {SH-1:sea_ratio}% as bright as
the sky. It darkens steadily to deep slate-teal ({SH-1:sea_bottom}) at the bottom of the frame, with a
faint olive cast in the nearest water. A moderate swell {SH-1:hs} m high and {SH-1:wl} m from crest to
crest rolls straight toward the camera from the horizon. Its crests run parallel to the horizon, and
the camera looks into their faces. Lively wind ripples under a {SH-1:wind} m/s breeze blowing toward the
camera cover the surface in a fine, crisp texture of pale and dark facets. Empty sea, unbroken water,
no foam, no land, no boats, no people.

**Frame data.**
{SH-1:frame_data}

### SH-2 · The Sun 4° up

*Day 11.7, eight hours before sunset. Sun {SH-2:sun_el}° up in the west. 6,800 lux; day vision.*

**Image prompt.** Photograph from a small boat on a glassy sea in the late afternoon, eye 2 m above the
water, 28 mm lens, looking west into the low Sun over an empty, perfectly flat horizon, horizon
{SH-2:horizon_pct}% from the top.

The Sun is a small, dazzling orange disk, {SH-2:sun_px} pixels across ({SH-2:sun_hex} at its own
exposure, blown out here), {SH-2:sun_el} degrees above the horizon at the centre of the frame
(x = {SH-2:sun_x}, y = {SH-2:sun_y}). Below it lies an almost mirror-calm sea. The Sun's reflection is
a narrow, intensely brilliant column of molten orange, no more than {SH-2:glit_w} pixels wide, running
from the horizon straight down to the bottom of the frame. It is broken into long, smooth, stretched
highlights by a low swell; its core is blinding white-orange and its edges sharp.

The sky is a warm amber-cream at the top ({SH-2:sky_top}), deepening to {SH-2:sky_8} at 8 degrees and to
soft orange at the horizon ({SH-2:sky_hz}). It is brightest above the Sun, with no halo.

Beside the glitter column the calm water mirrors the sky at {SH-2:sea_ratio_l}% of its brightness: warm
apricot ({SH-2:sea_below_l}) just below the horizon, darkening to bronze-olive ({SH-2:sea_bottom_l}) at
the bottom. A gentle, glassy swell {SH-2:hs} m high and {SH-2:wl} m from crest to crest rolls toward the
camera from straight ahead. Its smooth undulations bend the reflections into long wavy bands. There are
no wind ripples. Empty sea, no land, no foam, no boats, no people.

**Frame data.**
{SH-2:frame_data}

### SH-3 · The Sun 4° below the horizon

*Day 12.3, eight hours after sunset. 2,800 lux; day vision.*

**Image prompt.** Photograph from a small boat on a calm sea in the early evening twilight, eye 2 m
above the water, 28 mm lens, looking west to the place of sunset over a flat, empty horizon, horizon
{SH-3:horizon_pct}% from the top. The Sun is 4 degrees below the horizon at the centre of the frame.

The sky glows amber and orange:
- soft amber at the top ({SH-3:sky_top}), brightest above the centre;
- light orange ({SH-3:sky_8}) at 8 degrees;
- coral-orange just above the horizon ({SH-3:sky_hz}).

The sea is nearly calm and mirrors this sky almost exactly: apricot-orange ({SH-3:sea_below}) just below
the horizon, {SH-3:sea_ratio}% as bright as the sky above, fading to warm bronze ({SH-3:sea_bottom}) at
the bottom. The horizon line is barely visible, a faint darkening between two glowing halves of the
frame. A low glassy swell {SH-3:hs} m high and {SH-3:wl} m between crests rolls toward the camera.
Faint patches of capillary ripples, the first a {SH-3:wind} m/s breeze raises, dull the mirror slightly.
The swell's long smooth undulations stretch the reflected glow into gentle horizontal bands. Empty sea,
no land, no foam, no boats, no people.

**Frame data.**
{SH-3:frame_data}

### SH-4 · The Sun 15° below the horizon

*Day 13.3, 30 hours after sunset. 450 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the long evening, eye 2 m above the water,
28 mm lens, looking west over an empty horizon, horizon {SH-4:horizon_pct}% from the top.

A deep amber-orange twilight fills the western sky:
- warm amber at the top ({SH-4:sky_top});
- golden orange ({SH-4:sky_8}) at 8 degrees;
- a rich coral-orange band along the horizon ({SH-4:sky_hz}), brightest at the centre and duller toward
  both edges.

The sea, gently rippled by a {SH-4:wind} m/s breeze, reflects the sky at {SH-4:sea_ratio}% of its
brightness. It is amber-orange ({SH-4:sea_below}) below the horizon, softening into a smooth, slightly
grainy sheen, then darkening to brown-gold ({SH-4:sea_bottom}) at the bottom. A swell {SH-4:hs} m high
and {SH-4:wl} m from crest to crest comes from the front left, its crests running diagonally. Empty sea,
no land, no foam, no boats, no people.

**Frame data.**
{SH-4:frame_data}

### SH-5 · The Sun 35° below the horizon

*Day 14.9, nearly three days after sunset. 9.1 lux; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat at sea late in the long evening, eye 2 m
above the water, 28 mm lens, looking west over an empty horizon, horizon {SH-5:horizon_pct}% from the
top.

A deep twilight arch spans the western sky:
- glowing amber-gold about 8 degrees up at the centre ({SH-5:sky_8});
- vivid orange at the horizon ({SH-5:sky_hz});
- duller ochre at the top of the frame ({SH-5:sky_top});
- redder and darker toward both edges ({SH-5:sky_hz_l} low on the left).

It is the light of the very high air, which still sees the Sun 35 degrees below the horizon.

A confused low swell and a {SH-5:wind} m/s breeze roughen the sea. It reflects the arch at
{SH-5:sea_ratio}% of its brightness: orange ({SH-5:sea_below}) just below the horizon, deepening to dark
brown-orange ({SH-5:sea_bottom}) at the bottom. Its surface is a soft-grained texture of orange and
brown facets. The Earth, a thin crescent, is 3 degrees below the horizon and not visible. Muted to the
eye; the photograph records the colour. Empty sea, no land, no stars, no foam, no boats, no people.

**Frame data.**
{SH-5:frame_data}

### SH-6 · The darkest hour: the Earth rising

*Day 18.1, six days after sunset; the Sun 74° below the horizon. 0.085 lux; night vision.*

**Image prompt.** Long-exposure night photograph from a small boat on a mirror-calm sea, eye 2 m above
the water, 28 mm lens, looking west, horizon {SH-6:horizon_pct}% from the top.

The Earth has just risen. It sits {SH-6:earth_el} degrees above the horizon at the centre of the frame
(centre x = {SH-6:earth_x}, y = {SH-6:earth_y}), a disk {SH-6:earth_px} pixels across,
{SH-6:earth_lit}% lit. Its lit part is a thick crescent along the {SH-6:earth_side} edge of the disk,
like a bowl. It is deep orange ({SH-6:earth_hex}) from the long path through the air, overexposed with a
bloom in this exposure. The rest of the disk is invisible against the sky.

Its reflection runs from the horizon straight down toward the camera as a narrow, vivid orange column
about {SH-6:glit_w} pixels wide, broken into a few long smooth bands by a very low swell. It is an
almost perfect mirror image, reaching to row {SH-6:glit_y1}, about halfway from the horizon to the
bottom of the frame.

The sky is a faint warm glow:
- pale apricot-cream overhead ({SH-6:sky_top});
- golden-orange at 8 degrees ({SH-6:sky_8});
- soft salmon-orange at the horizon ({SH-6:sky_hz}), lit by the low Earth.

{SH-6:stars_recorded} stars are recorded: {SH-6:stars_named}, a small point near the top of the frame,
and fainter ones at {SH-6:stars_list}.

The water away from the Earth's reflection is glassy and mirrors the sky at {SH-6:sea_ratio_l}% of its
brightness: salmon-tan ({SH-6:sea_below_l}) at the horizon, darkening to dark olive-brown
({SH-6:sea_bottom_l}) at the bottom. A very low swell {SH-6:hs} m high and {SH-6:wl} m between crests
barely moves it, without ripples.

To the eye this scene is a dark, almost colourless night with the orange Earth and its path. Empty sea,
no land, no foam, no boats, no people.

**Frame data.**
{SH-6:frame_data}

## Southern Mare Nubium, by Pitatus

The coast with the largest tide on the exposed tenth of the nearside's shores: this month's range is
7.7 m. The wave grid's shore point falls on the west shore of an island 4 by 6 km, rising to about
290 m. The camera floats 1.5 km off that shore (11.93° W, 28.11° S) in 65 m of water, with open sea to
the west and north. The mainland coast, with hills and crater rims of 400 to 1,100 m, lies 12 to 30 km
to the south and southwest. The Earth stands high in the northeast, 53 to 65 degrees up, out of every
frame.

### SN-1 · The Sun at its highest

*Day 13.2, seven and a half days before sunset. Sun {SN-1:sun_el}° up in the north, behind the camera.
75,000 lux; day vision.*

**Image prompt.** Photograph from a small boat on a choppy sea at midday, eye 2 m above the water,
28 mm lens, looking south with the high Sun behind the camera, tilted 3 degrees down, horizon
{SN-1:horizon_pct}% from the top.

On the left rises the sunlit island, 1.7 to 5 km away: a smooth-shouldered hill of bare, dark
grey-brown ground ({SN-1:land_hex}). It rises to 120 pixels above the horizon at the left edge (280 m
high) and slopes down to the water about a fifth of the way across the frame. It is front-lit and only
lightly hazed. It is high tide, {SN-1:tide} m from the mean level, and the sea laps high on the island's
lowest slopes. Beyond it, across the rest of the horizon, runs a faint low line of mainland hills 17 to
30 km away, 38 to 65 pixels tall, pale hazy blue-grey.

The sky is very pale cyan-blue, almost white ({SN-1:sky_top} at the top), fading to pale
grey-green-cream at the horizon ({SN-1:sky_hz}).

The sea is the roughest of the 24 scenes:
- a {SN-1:wind} m/s wind blows from behind the camera, over the right shoulder, raising steep, crowded
  wind waves and ripples over a {SN-1:hs} m swell {SN-1:wl} m from crest to crest that runs away to the
  left;
- the surface is a busy, high-contrast mosaic of pale sky-reflecting facets and dark slate troughs;
- grey-teal ({SN-1:sea_below}) below the horizon, darkening to deep slate ({SN-1:sea_bottom}) at the
  bottom, with a faint olive cast in the nearest water.

Unbroken water, no foam. No boats, no people, no plants.

**Frame data.**
{SN-1:frame_data}
- At {SN-1:wind} m/s the first whitecaps would form; the study does not model them.

### SN-2 · The Sun 4° up

*Day 20.3, nine hours before sunset. Sun {SN-2:sun_el}° up in the west. 6,900 lux; day vision.*

**Image prompt.** Photograph from a small boat on a mirror-calm open sea in the late afternoon, eye 2 m
above the water, 28 mm lens, looking west into the low Sun, horizon {SN-2:horizon_pct}% from the top.

The Sun is a small, blinding orange disk {SN-2:sun_px} pixels wide, {SN-2:sun_el} degrees above the
horizon at the centre (x = {SN-2:sun_x}, y = {SN-2:sun_y}). Below it a narrow, brilliant column of
reflected sunlight, no more than {SN-2:glit_w} pixels wide, runs from the horizon nearly to the bottom of
the frame. It is molten white-orange at its core and stretched into long smooth bands by a low swell.

The sky is warm amber-cream at the top ({SN-2:sky_top}), deepening to soft orange at the horizon
({SN-2:sky_hz}), brightest above the Sun.

On the far left, a thin, faint line of distant coast 12 to 23 km away rises only 2 to 20 pixels above
the horizon, a hazy dull red-brown ({SN-2:land_hex}); the rest of the horizon is open sea.

The calm water beside the glitter mirrors the warm sky at {SN-2:sea_ratio_l}% of its brightness:
apricot ({SN-2:sea_below_l}) at the horizon, darkening to bronze-olive ({SN-2:sea_bottom_l}) at the
bottom. A smooth swell {SN-2:hs} m high and {SN-2:wl} m from crest to crest comes from the front right;
its long glassy undulations bend the reflections. There are no ripples: the wind is under 1 m/s. No
foam, no boats, no people.

**Frame data.**
{SN-2:frame_data}

### SN-3 · The Sun 4° below the horizon

*Day 21.1, nine hours after sunset. 2,870 lux; day vision.*

**Image prompt.** Photograph from a small boat on a mirror-calm sea in the early evening twilight, eye
2 m above the water, 28 mm lens, looking west to the place of sunset, horizon {SN-3:horizon_pct}% from
the top.

The sky is a luminous gradient:
- soft amber at the top ({SN-3:sky_top});
- orange through the middle;
- coral ({SN-3:sky_hz}) at the horizon, brightest above the centre where the Sun has gone down.

The sea is almost perfectly calm and mirrors the sky at nearly its full brightness. Just below the
horizon it is apricot-orange ({SN-3:sea_below}), {SN-3:sea_ratio}% as bright as the sky just above, so
the horizon line all but disappears. It darkens to warm bronze ({SN-3:sea_bottom}) at the bottom of the
frame. A long, low glassy swell {SN-3:hs} m high and {SN-3:wl} m between crests rolls in from the front
right, drawing slow, broad bends into the mirrored glow.

On the far left a hair-thin, hazy line of distant coast rises 1 to 18 pixels, dull red-brown
({SN-3:land_hex}). No ripples, no foam, no boats, no people.

**Frame data.**
{SN-3:frame_data}

### SN-4 · The Sun 15° below the horizon

*Day 22.1, a day and a half after sunset. 466 lux; day vision.*

**Image prompt.** Photograph from a small boat on a wind-roughened sea in the long evening, eye 2 m
above the water, 28 mm lens, looking west, horizon {SN-4:horizon_pct}% from the top.

The western sky is a broad twilight of warm amber ({SN-4:sky_top} at the top), deepening to a
coral-orange band along the horizon ({SN-4:sky_hz}), brightest at the centre.

A thin, hazy line of distant coast 12 to 23 km away runs along the left half of the horizon, 5 to 21
pixels tall, dull red-brown ({SN-4:land_hex}).

The sea is choppy under a {SN-4:wind} m/s wind. It reflects the sky at {SN-4:sea_ratio}% of its
brightness as a busy texture of orange highlights and brown troughs: bronze-orange ({SN-4:sea_below})
below the horizon, dark brown-gold ({SN-4:sea_bottom}) at the bottom. A swell {SN-4:hs} m high and
{SN-4:wl} m from crest to crest runs away from the camera, slightly to the right; short steep wind
waves ride on it. Unbroken water, no foam, no boats, no people.

**Frame data.**
{SN-4:frame_data}
- At this wind the first whitecaps would form; the study does not model them.

### SN-5 · The Sun 35° below the horizon

*Day 24.1, three and a third days after sunset. 11.7 lux; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat at sea late in the long evening, eye 2 m
above the water, 28 mm lens, looking west-southwest, horizon {SN-5:horizon_pct}% from the top.

A deep golden-orange twilight arch fills the western sky:
- brightest about 8 degrees up at the centre ({SN-5:sky_8});
- vivid orange along the horizon ({SN-5:sky_hz});
- duller ochre at the top ({SN-5:sky_top});
- redder toward the edges.

On the left, coastal hills 18 km away rise as a dark red-brown silhouette ({SN-5:land_hex}), 63 pixels
tall at the left edge. They fall toward the centre into a thin, hazy line of distant coast that ends
two-thirds of the way across.

The sea, rippled by a {SN-5:wind} m/s breeze over a {SN-5:hs} m swell running away to the right,
reflects the arch at {SN-5:sea_ratio}% of its brightness: orange ({SN-5:sea_below}) below the horizon,
dark brown-orange ({SN-5:sea_bottom}) at the bottom. Muted to the eye; the photograph records the
colour. No stars, no foam, no boats, no people.

**Frame data.**
{SN-5:frame_data}

### SN-6 · The darkest hour

*Day 0.5, nine days after sunset; the Sun 52.8° below the horizon, straight ahead. 2.9 lux, mostly
from a high, nearly full Earth behind the camera; between day and night vision.*

**Image prompt.** Long-exposure night photograph from a small boat, eye 2 m above the water, 28 mm lens,
looking southeast toward the island, horizon {SN-6:horizon_pct}% from the top.

The island fills the lower part of the view from edge to edge, 1.5 to 5 km away. It is a long, rounded
ridge of bare ground rising 51 to 110 pixels above the horizon, highest at the centre (250 to 290 m).
Beyond its right end, far hazy hills of the mainland 23 km away stand 51 to 58 pixels tall. The
island's slopes are lit softly from behind the camera by a high, nearly full Earth (out of frame), so
they show as a dim muted red-brown ({SN-6:land_hex}), their shapes visible.

It is low tide, {SN-6:tide} m from the mean level: a band of shore 3 m high, 2 to 3 pixels, lies
exposed along the island's foot at the waterline. The study does not calculate how wet ground looks.
Drawn darker and glistening, the band follows Earth's shores, which is an assumption here; drawn in
the island's soil colour, it follows the calculation.

The sky is a pale, nearly neutral light grey at the top ({SN-6:sky_top}), warming downward into a soft
pinkish peach glow above the hidden Sun: {SN-6:sky_8} at 8 degrees and {SN-6:sky_hz} at the horizon.

The sea is lightly rippled by a {SN-6:wind} m/s breeze over a short swell {SN-6:hs} m high and
{SN-6:wl} m from crest to crest, coming toward the camera. It reflects the glow as a dim greyish
rose-taupe ({SN-6:sea_below}) below the horizon, darkening to dark grey-brown ({SN-6:sea_bottom}) at the
bottom. No stars show in this bright earthlit sky. No foam, no boats, no people.

**Frame data.**
{SN-6:frame_data}
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

*Day 28.4, seven and a half days before sunset. Sun {IC-1:sun_el}° up in the north, behind the camera.
66,000 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea at midday, eye 2 m above the water, 28 mm lens,
looking south across 5 km of open water to a mountainous coast, with the high Sun behind the camera,
tilted 2 degrees down, horizon {IC-1:horizon_pct}% from the top.

A range of bare, sunlit mountains spans the whole width of the frame, 10 to 21 km away. It rises 77
pixels above the horizon at the left edge, climbs in rounded peaks and shoulders to 174 pixels at
two-thirds of the way across (1,270 m), then falls to 115 pixels at the right edge. The mountains are
pale tan-beige highland ground ({IC-1:land_hex}), softened and cooled by haze. Their front faces are
lit, with gentle shading in the folds and no vegetation.

The sky is a very pale cyan-blue, almost white ({IC-1:sky_top} at the top), fading to pale cream-grey at
the horizon behind the peaks ({IC-1:sky_hz}).

The sea is a grey-teal slate ({IC-1:sea_below} below the horizon), darkening to deep slate-teal
({IC-1:sea_bottom}) at the bottom. The band of water nearest the horizon mirrors the pale mountains in
blurred streaks. A swell {IC-1:hs} m high and {IC-1:wl} m from crest to crest comes from the front
left; wind ripples under a {IC-1:wind} m/s breeze give it a lively, fine texture. No foam, no boats,
no people, no plants, no buildings.

**Frame data.**
{IC-1:frame_data}

### IC-2 · The Sun 4° up

*Day 5.9, ten hours before sunset. Sun {IC-2:sun_el}° up in the west. 6,900 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the late afternoon, eye 2 m above the water,
28 mm lens, looking west along a mountainous coast into the low Sun, horizon {IC-2:horizon_pct}% from
the top.

The Sun is a small, blinding orange disk {IC-2:sun_px} pixels wide, {IC-2:sun_el} degrees up just right
of centre (x = {IC-2:sun_x}, y = {IC-2:sun_y}). Below it a broad, sparkling orange glitter path, up to
{IC-2:glit_w} pixels wide, runs from the horizon to the bottom of the frame over a rippled sea.

On the left the coast runs out toward the horizon:
- a near coastal ridge 20 km away stands as a dark, hazy red-brown silhouette ({IC-2:land_hex}),
  99 pixels tall at the left edge;
- beyond it, from an eighth to half of the way across, a far, much higher mountain massif 50 to 72 km
  away rises 22 to 63 pixels above the horizon;
- the massif is so hazed that it is only a faint, pale, warm silhouette barely darker than the sky
  behind it.

The right half of the horizon is open sea.

The sky is warm amber-cream at the top ({IC-2:sky_top}), deepening to soft orange at the horizon
({IC-2:sky_hz}). The sea beside the glitter is warm tan-gold ({IC-2:sea_below_l}), darkening to dark
olive-bronze ({IC-2:sea_bottom_l}) in the lower corners. A swell {IC-2:hs} m high and {IC-2:wl} m
between crests comes from the front left, and a {IC-2:wind} m/s breeze ripples it. No foam, no boats,
no people.

**Frame data.**
{IC-2:frame_data}

### IC-3 · The Sun 4° below the horizon

*Day 6.8, ten hours after sunset. 2,850 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the early evening twilight, eye 2 m above the
water, 28 mm lens, looking west along a mountainous coast, horizon {IC-3:horizon_pct}% from the top.

The sky glows:
- soft amber at the top ({IC-3:sky_top});
- orange through the middle;
- coral ({IC-3:sky_hz}) along the horizon, brightest just right of centre where the Sun has set.

On the left the coast's near ridge, 15 to 21 km away, stands as a dark, hazy red-brown silhouette
({IC-3:land_hex}), 60 to 105 pixels tall. Beyond it, across the middle of the frame, the far high massif
50 to 72 km away shows as a pale, ghostly warm silhouette 25 to 58 pixels tall, little darker than the
sky. The right third of the horizon is open sea.

The sea, gently rippled by a {IC-3:wind} m/s breeze over a {IC-3:hs} m swell, reflects the glow at
{IC-3:sea_ratio}% of its brightness: apricot-orange ({IC-3:sea_below}) below the horizon, warm bronze
({IC-3:sea_bottom}) at the bottom. The water just below the horizon on the left mirrors the dark ridge
in soft, broken streaks. No foam, no boats, no people.

**Frame data.**
{IC-3:frame_data}

### IC-4 · The Sun 15° below the horizon

*Day 7.9, a day and a half after sunset. 456 lux; day vision.*

**Image prompt.** Photograph from a small boat at sea in the long evening, eye 2 m above the water,
28 mm lens, looking west-southwest along a mountainous coast, horizon {IC-4:horizon_pct}% from the top.

A broad amber twilight fills the sky: warm amber at the top ({IC-4:sky_top}), deepening to a
coral-orange band along the horizon ({IC-4:sky_hz}).

The near coastal mountains, 15 to 21 km away, stand on the left as a dark red-brown silhouette
({IC-4:land_hex}), 77 to 114 pixels tall. They lead into the pale, hazed far massif across the middle,
50 to 72 km away and 24 to 61 pixels tall, its outline soft against the glow. The right quarter is
open sea.

The sea, rippled by a {IC-4:wind} m/s breeze over a {IC-4:hs} m swell crossing from left to right,
reflects the sky at {IC-4:sea_ratio}% of its brightness: amber-orange ({IC-4:sea_below}) below the
horizon, brown-gold ({IC-4:sea_bottom}) at the bottom. The mountains' dark reflection blurs into the
water below them on the left. No foam, no boats, no people.

**Frame data.**
{IC-4:frame_data}

### IC-5 · The Sun 35° below the horizon

*Day 10.2, nearly four days after sunset. 8.9 lux; between day and night vision.*

**Image prompt.** Long-exposure photograph from a small boat at sea late in the long evening, eye 2 m
above the water, 28 mm lens, looking southwest along the mountainous coast, horizon {IC-5:horizon_pct}%
from the top.

A deep twilight arch glows over the mountains:
- brightest amber-gold about 8 degrees up at the centre ({IC-5:sky_8});
- vivid orange along the horizon ({IC-5:sky_hz});
- ochre at the top ({IC-5:sky_top});
- redder toward the edges.

The coastal mountains fill the whole width of the horizon as a dark red-brown silhouette
({IC-5:land_hex}):
- tallest at the left edge, 178 pixels (1,240 m at 12 km);
- stepping down through rounded shoulders to about 90 pixels across the centre;
- passing on the right into the far high massif, 51 to 72 km away and 27 to 61 pixels tall, paler with
  distance.

The sea is nearly calm, with no wind ripples under a {IC-5:wind} m/s breeze, only a smooth {IC-5:hs} m
swell coming toward the camera. It mirrors the glow at {IC-5:sea_ratio}% of its brightness: orange
({IC-5:sea_below}) below the horizon, with the mountains' dark shapes reflected in the nearest band of
water, darkening to brown-orange ({IC-5:sea_bottom}) at the bottom. Muted to the eye; the photograph
records the colour. No stars, no foam, no boats, no people.

**Frame data.**
{IC-5:frame_data}

### IC-6 · The darkest hour: the far side's twilit night

*Day 13.7, seven and a third days after sunset; the Sun 51.8° below the southern horizon, straight
ahead. 0.46 lux, all from the twilight; between day and night vision.*

**Image prompt.** Long-exposure night photograph from a small boat at sea, eye 2 m above the water,
28 mm lens, looking south toward a mountainous coast, horizon {IC-6:horizon_pct}% from the top.

All night a golden twilight arch stands over the southern mountains, lit by a Sun more than 50 degrees
below the horizon:
- warm amber-gold, brightest about 8 degrees up ({IC-6:sky_8});
- orange near the mountain tops ({IC-6:sky_3} at 3 degrees);
- ochre-orange at the top of the frame ({IC-6:sky_top}).

The mountain range, 10 to 21 km away, spans the whole frame as a dark silhouette: 76 pixels tall at the
left edge, rising to 172 pixels two-thirds of the way across (1,270 m), and 115 pixels at the right
edge. The haze tints it a deep rust-brown ({IC-6:land_hex}) in this long exposure. It hides the arch's
lowest, reddest part.

The sea is lightly rippled by a {IC-6:wind} m/s breeze over a {IC-6:hs} m swell running away to the
right. It reflects the arch as a warm orange sheen ({IC-6:sea_below} below the horizon at the centre),
with the mountains' dark outline mirrored, blurred, in the band of water just below the horizon. It
darkens to brown-orange ({IC-6:sea_bottom}) at the bottom.

The Earth never rises here. No stars show against this glowing sky; the brightest in view, Achernar, is
too dim through the air. To the eye: a dim, nearly colourless twilight with a pale glow over black
mountains. No foam, no boats, no people.

**Frame data.**
{IC-6:frame_data}

## Boundaries

- **Sky and sea colours** come from the model of the results by regime, evaluated for each frame's
  directions with its own slopes and camera. The sea in it is statistical: no individual waves, and open
  water in every direction. Where land rises behind the sea, the scene text says how the nearest water
  mirrors it, without numbers.
- **The coast** is the LOLA altimetry grid, about 118 m between points, at the moment's tide. Its colour
  is a placeholder: bare Apollo soil, lit by the clear sky and the disks on a slope facing the camera,
  and hazed by extinction at the ground toward the colour of the horizon sky. Beaches, cliffs, surf, wet
  ground and plants are not modelled. The Smythii frames show no land: the headland is behind the camera.
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
