# Atlas sheets

Map sheets of the Open Moon at its selected standing-water share. They draw the geography [atlas](../../geography/README.md#atlas-at-the-selected-water-share) product and the heritage sites and scientific targets of the [conservation study](../../research/studies/conservation/README.md).

```sh
python -m geography.atlas
python -m research.studies.conservation.run
climate/gcm/.venv/bin/python -m climate.gcm.climatology A28_corrected_970:10-19   # the design case's clouds and rain
python visualization/atlas/render.py      # out/: three sheets, globe.html and manifest.json (kept out of Git)
```

The globe on the turned Moon's design climate (9 November 2026), from the climate domain's corrected run:

```sh
V=climate/gcm/.venv/bin/python; O=<domain/geology checkout>/geology/results/moon_orientation.json
C=A28R_301k_design_l15_snow_s938_turned_refit_precession
$V -m climate.gcm.climatology $C:10-19
$V -m climate.gcm.turned_climatology $O climate/gcm/products/climatology_$C.npz \
    --source domain/geology@45c3d817:geology/results/moon_orientation.json
$V -m geography.drainage --out geography/turned_corrected --climatology climate/gcm/products/climatology_${C}_today.npz
python3 visualization/atlas/render.py --frame $O --source domain/geology@45c3d817:geology/results/moon_orientation.json \
    --climatology climate/gcm/products/climatology_$C.npz --drainage geography/turned_corrected
```

The 7 November globe ran the same chain on `A28R_301k_design_l15_snow_s924_turned_refit` into `geography/turned`. The 16 November globe adds the years control: `research.studies.drainage_network.stops` writes the drainage network's water at eight years, and the last command takes `--stops research/studies/drainage_network/results/globe_stops.json`.

The globe on the atlas round's design surface (17 November 2026), from the round's surface product on the research drive, since 18 November on the design climate as estimated by pattern scaling (the run on the design surface at sunlight 0.936, years 45–59, carried to 300.31 K by `climate.gcm.pattern_scaling`; until then the run at 0.930 on the final surface, years 35–49):

```sh
C=climate/gcm/products/climatology_A28D_301k_design_l15_snow_s936_settled_y45-59_scaled300.31K.npz
W=climate/gcm/products/surface_wind_A28D_301k_design_l15_snow_s936_settled_y45-59_scaled300.31K.npz
python3 visualization/atlas/render.py --surface /media/projectspace/terluna-research/atlas-round/globe_surface/globe_surface_16ppd.npz \
    --surface-record research/studies/atlas_round/results/globe_surface.json --climatology $C \
    --moisture research/studies/atlas_round/results/moisture_map_final.json \
    --falls research/studies/atlas_round/results/falls.json --wind $W \
    --districts research/studies/atlas_round/results/districts.json
```

With `--moisture` (since 17 November) the ground is drawn in true colour (below, The ground in true colour) and the page's textures are written as files beside it, in `out/textures/`, which a publish carries with the page. The moisture map must be of the climatology's run and years; the build stops otherwise. It takes about 8 minutes, half of it the sky each place sees. `--wind` takes the climate domain's surface wind of the same run and years, made with `climate/gcm/.venv/bin/python -m climate.gcm.surface_wind <run>:<first>-<last>` (about a minute), which sets the sea's slopes (below, The sea's slopes); the build stops on a different run or years. An estimate's climatology and surface wind come from `climate.gcm.pattern_scaling` together, and the page's climate sentence reads its ground, rain and drylands from the climatology and the moisture map it draws (`climate_figures` in render.py) where the run's entry in `TURNED_CLIMATES` leaves them open.

| Sheet | Shows |
|---|---|
| `atlas_<pct>pct_near_far.png` | The near side as seen from Earth and the far side, in orthographic projection |
| `atlas_<pct>pct_global.png` | The whole Moon in longitude and latitude, with a 30° grid |
| `atlas_<pct>pct_poles.png` | Both polar regions poleward of 80°, from the product's 16 px/deg polar caps, with the permanently shadowed craters |
| `globe.html` | An interactive 3D globe in one self-contained page (three.js from cdnjs, textures embedded), with an appearance mode and a topography mode, which can also lie over the appearance: views of the near side, far side, poles and the Moon from Earth, the Sun set by the Moon's age or behind the viewer, labels, a grid at chosen spacings, going to a place, the height or depth of a clicked place, and a record card for each heritage site and target |

## How the sheets are drawn

- **Water.** Six depth bands of one blue hue: 0–250 m, 250–500 m, 0.5–1 km, 1–2 km, 2–4 km and over 4 km.
- **Land.** Neutral shaded relief, lit from the north-west.
- **Symbols.** Orange marks heritage sites: circles for submerged sites and triangles for sites on land. Dark squares mark scientific targets.
- **Labels.**
  - Sea labels are working names at hand-placed positions.
  - Island labels come from the IAU features nearest each summit.
  - Site labels follow the register's identifiers.
- **Polar sites.** Sites nearer the poles than 78° are labelled on the polar sheet.

## The globe

[globe.template.html](globe.template.html) holds the viewer. `render.py` fills it with the appearance data from
[appearance.py](appearance.py) and with the relief of the atlas grid. The appearance mode's surface is 16 pixels
per degree (1.9 km), with relief shading at 8 (at 16 on the design surface's true-colour page). The topography
takes its heights from the atlas's 4 and its shading from the same relief. JPEG keeps the single-file page near
12 MB; the true-colour page keeps its textures as files beside it.

On 16 October the globe took its clouds and ground cover from the corrected design case in place of run A, and
the lakes as the drainage product holds them, 10.3% of the Moon in place of the first estimate's 12%, and was
republished to its private claude.ai page. The lakes follow the design case once it settles, after the glow pair
([register](../../research/decisions.md#the-globe-and-the-glow-pair-16-october-2026)). The same day the
topography replaced the map mode, as a mode of its own and as an overlay on the appearance, and the globe gained
the fully lit option, the grid's spacings, going to a place and the readout of a clicked place.

On 7 November the globe was redrawn on the turned Moon and republished to the same page. With `--frame`, render.py carries the 4 and 16 pixel/degree heights, the rivers, the lakes and every label into the frame of domain/geology's orientation product with climate.gcm.turned's rotation; each turned pixel takes the height at its place today, the nearest cell's, which keeps the atlas's water share (appearance.turn_texture), and each channel cell lands in the turned pixel holding its centre as well, so channels stay continuous (turn_channels). The clouds and the ground cover come from the design run on the turned Moon at the target (`A28R_301k_design_l15_snow_s924_turned_refit`, years 10–19), whose climatology is in the turned frame already. The rivers and lakes come from the drainage routed on that run's rain and evaporation carried back to today's ground (climate.gcm.turned_climatology, `geography/turned/`): 33,262 rain-fed lakes over 10.28% of the Moon holding 2.84 million km³ and 199,217 m³/s reaching the seas, against the September climate's 27,865 lakes, 3.11 million km³ and 229,410 m³/s. The water level, the lakes' levels and every height stay on the dry Moon's geoid, as on the turned run's own surface. The sheets stay on today's grid.

On 9 November the globe was redrawn on the corrected design climate and republished to the same page (version 9). The clouds and the ground cover come from `A28R_301k_design_l15_snow_s938_turned_refit_precession`, years 10–19, the design climate on the turned Moon, with PlaSim's night-side twilight counted once and the equinox regressing with the node (climate/gcm/README.md, The corrected design climate on the turned Moon). The rivers and lakes come from the drainage routed on its rain and evaporation carried back to today's ground (`geography/turned_corrected/`): 32,636 rain-fed lakes over 9.80% of the Moon holding 2.65 million km³ and 199,973 m³/s reaching the seas, against the target run's 33,262 lakes, 10.28%, 2.84 million km³ and 199,217 m³/s. The page's paragraph on the turned Moon takes each run's ground, rain and drylands from `TURNED_CLIMATES` in render.py, and the build stops for a run the table does not hold. The build's own checks hold: the relief within 0.5 m, and the seas' and lakes' shares on the page as the products give them. The sheets, on today's grid, are unchanged.

On 16 November the globe gained a years control and was republished to the same page (version 10). The control steps the water through the drainage network study's history of the six large seas with the routing's water budget closed, at 0, 250, 500, 750, 1,000, 2,000, 5,000 and 10,000 years ([the globe's years](../../research/studies/drainage_network/README.md#the-globes-years)), and the page opens at year 0, the water as the atlas draws it. With `--stops`, render.py reads the study's `globe_stops.json` and its maps, carries the maps into the turned frame as it carries the surface grid, and embeds two more textures. The first holds, for each pixel, the years in which one of the study's bodies of water covers it and the first and second body to cover it. The second holds the years in which that second body is the pixel's surface, and the overflow rivers on dry ground: their years and the share of the pixel that a river of their discharge covers. The shaders take each body's level in the year shown and colour its water by its depth over the 0.25° heights with the surface colour's own depth ramp, draw an overflow river at the width its discharge gives, and give floor a sea has left the basalt soil colour, a placeholder the legend names. Every other water stays as version 9 draws it. A clicked place reads the year's water and its depth there, and `y<years>` in the link opens the page at the recorded year nearest it. The build checks each year's drawn water share against the study's within 0.002 of the Moon (all eight agree within 0.0001) and counts the pixels where a body at year 0 lies off the atlas's sea, which it allows up to 0.05% of the grid; this build counts none.

The same build smooths the climate run's rain and soil moisture before they set the ground cover, with a Gaussian round on the sphere whose width is 0.7 of the run's 5.5° grid spacing, about 3.9° (`CLIMATE_SMOOTHING_CELLS` in appearance.py). Interpolated straight from the run's 32 × 64 grid, the fields carried kinks along the grid's rows and columns and ripples at its spacing, finer than the spectral model resolves; the cover's thresholds turned them into straight east–west streaks along the turned Moon's latitude circles, which are the run's rows. At 0.7 of a cell, ripples two cells long keep a tenth of their contrast and waves 30° long about 70%.

**Appearance mode** estimates how the Open Moon looks from outside its air. One full-screen pass traces each
pixel's ray through the atmosphere, two cloud decks and the surface:

- **Computed:**
  - land, sea and relief: LOLA above the GRAIL geoid at 16 pixels per degree through the geography domain,
    filled to the atlas product's sea level;
  - rivers and rain-fed lakes above sea level, from the geography domain's drainage product (10.3% of the Moon
    in lakes, routed with the dimmer design case's climate before its radiation was corrected); river width
    follows discharge;
  - the air, from the illumination domain's engine atmosphere: optical depth and scale height;
  - the halo around the limb and the sunlit air beyond the terminator, which come out of the same ray-march
    through the air's full 484 km;
  - the ground's direct and diffuse light at every Sun height, from the Open Moon sky atlas.
- **Informed:**
  - the low and high cloud decks, from the climate domain's GCM climatology of the design case
    (`A28_corrected_970`, years 10–19; run A's until 16 October). Each moves with the Sun
    through the run's composite by latitude and local time and is carried east by its mean wind: 1.3 m/s for
    the low deck and 6.3 m/s for the high deck, whose clouds are drawn out east-west as that wind is;
  - light scattered more than once in the air, an isotropic source proportional to the sky atlas's diffuse
    light, scaled so that the air seen straight down matches Eddington's reflectance.
- **Guesstimate:**
  - surface cover, from rainfall and soil moisture (smoothed over 0.7 of the climate run's grid spacing),
    nearness to water, floodplains, height, hollows and slope;
  - water colour, from depth, with sediment-laden rivers;
  - river width, from Earth's width–discharge relation widened by 1.6 for the slower flow of lunar gravity;
  - cloud shapes and opacity. The shapes are noise that the page evaluates in 3D, so they stay sharp at any
    zoom. `appearance.py` replicates the page's hash and measures the noise's distribution, so that the
    share of sky each deck covers equals the run's cover.
- **Display:** exposure in stops, D65 white with no white balance, a soft shoulder above 80%, dithering
  against banding, and a view in red light only.

Towards the poles, where the grid's columns crowd together, the surface is worked out round on the sphere. Each
east-west slope spans as many columns as make up the north-south spacing, the cover's noise is drawn in 3D, its
hollows and floodplains are blurred alike in every direction on the ground, and nearness to water is measured on
the sphere. The page takes each sample's detail from the ground a pixel covers, a column counting for
cos(latitude) of a row, so the poles show neither spokes nor radial smears.

Refraction, calibrated multiple scattering along the limb and in twilight, Earthshine and clouds with depth
are the next refinements.

**Topography** colours the ground by its height above sea level and the sea floor by its depth below it, with
the water taken away, on a scale of 18 stops interpolated in sRGB:

- below sea level each of the sheets' six depth bands keeps its colour at the band's middle depth, so that the
  globe and the sheets read alike;
- contours mark a chosen interval (1 km at first) and draw sea level darker;
- light that every place has alike, from 60° up in its own north-west, shades the full relief, so that level
  ground keeps the legend's colours anywhere on the globe;
- as a mode, it raises the relief by the slider's exaggeration, down into the sea floors as well as up into the
  highlands, and carries the labels and markers the map mode carried (`#map` links open it);
- over the appearance, a checkbox lays the same colours and contours on the rendered globe at a chosen opacity,
  unlit by the Sun so that they read on the night side and through the air;
- rain-fed lakes above sea level are left out, so their beds show.

It reads two textures. The relief texture holds the atlas grid's height relative to sea level in whole metres,
its high byte in red and its low byte in green, which the page decodes. The tangent-space normal map holds the
full relief, which the appearance lays flat under water.

**In both modes:**

- the grid draws latitude and longitude at a major and a minor spacing inside the shaders, so its lines lie on
  the relief;
- Go to takes degrees (45 N, 30 W, or 45, −30 with east positive) or a name the globe labels, turns the globe to
  the place and leaves a pin;
- a click on the globe pins that place and reads the height of its 0.25° atlas cell;
- in the appearance mode, Fully lit gives every place the same Sun, 60° up in the north-west of its own sky and
  overhead near the poles, so that the globe reads evenly from limb to limb, with no sunglint; it keeps the clouds
  of the slider's day. The dense air veils the ground in blue there as in sunlight (its vertical optical depth is
  2.0 in blue); turning off the air, or red light only, shows the ground beneath.

The page opens in the view named by its link: `#near`, `#far`, `#south`, `#north`, `#earth` or `#topo` (`#map`
opens the topography too). It can add:

- `e<degrees>` for the Sun's elongation from Earth (`#earth.e60` is a crescent), and `lit` for the Sun behind the
  viewer;
- `x<stops>` for exposure;
- `z<tenths>` for the camera's distance in Moon radii (from Earth, the field of view in tenths of a degree);
- `lat<degrees>n|s` and `lon<degrees>e|w` to centre a place, and `pin` to mark it;
- `grid`, `g<degrees>` for the major spacing and `gm<degrees>` for the minor, with `p` for a decimal point
  (`gm0p5`);
- `c<metres>` for the contour interval (`c0` for none);
- `overlay` for the topography over the appearance, and `o<percent>` for its opacity;
- `y<years>` for the drainage network's water at the recorded year nearest it (`#far.y1000`);
- `red`, `noclouds` or `noair`;
- `peel<km>` for the air in view: only the air below that height, at the control's nearest stop below it (`#far.peel30`; `peel0` for none of it).

The Air in view control peels the air away from the top down, as a descent through it would show the ground. Its stops are the whole air (484 km deep, the view from space), 250, 150, 100, 75, 50, 30, 20, 10, 5, 2 and 1 km, and the ground. The page ray-marches only the air below the chosen height: the column between each point and the eye ends where the ray crosses that height, the haze and the halo come from the air below it, and a cloud deck above it leaves the view. Every place keeps the light the whole air gives it: the sunlight through the full column, the sky's diffuse light and the decks' shade on the ground. With the air's 48.4 km scale height, the air below 100 km holds 87% of the column, below 30 km 46% and below 10 km 19%. The camera stays outside the air (the nearest zoom is 608 km up, above its 484 km top); a camera inside it goes with the close views.

The page was also published privately to claude.ai.

On 17 November the globe was redrawn on the atlas round's design surface and republished to the same page (version 11), without a years control. With `--surface`, render.py draws the globe alone from the round's surface product ([globe_surface.py](../../research/studies/atlas_round/globe_surface.py)): the floor as the water's weight redraws it at the settled state, turned 6.97° into the settled frame the product names, every sea and lake at its own level, and the rivers routed on the design climate on the final surface with the near-side sea parted at its neck. The heights are measured from the common level. The appearance takes every water body as a lake at its own level (appearance.surface_grids) in place of the sea below one level. The topography gets one more texture, each 0.25° cell's mean water depth and the body over it, in which a sea takes the depth band's colour for its depth under its own surface, the lakes are left out so their beds show, and ground below the common level that no sea covers takes the lowland's colour; a click reads the water over a place, its depth and the sea's level. The heritage sites and targets are placed on the surface with the water over them there, and the five whose water differs from the register's say so on their card. The figures, the page's notes and the provenance come from the surface's record, and the climate sentence from `TURNED_CLIMATES`. The build checks the surface's water share against the product's within 0.002 of the Moon (0.3755 against 0.3755) and the map texture's within 0.005 (0.3736), and the relief within 0.5 m. The sheets are not drawn in this mode. A copy of the build and headless views are on the research drive in `atlas-round/globe_surface/globe_v11_2026-11-17/`.

### The ground in true colour

On 17 November the design surface's globe took its ground in true colour, as a satellite mosaic shows Earth's, in place of six cover types with one colour each, and more shading from the relief. The colours are reflectances in linear sRGB that the page lights with its computed sunlight and skylight; with the atmosphere and clouds off the page shows them as a mosaic does, the air's haze and the clouds taken away.

- **Plants** follow the design climate's dryness. The moisture map of the climatology's run and years gives each of the run's 170 km land cells its aridity, rain over the water the air could draw off; the run's water cells take their nearest land cell's, as the map assigns the design land there, and the field is smoothed as the clouds' fields are (appearance.moisture_field). Within a cell, ground along larger rivers and in hollows counts as wetter and steep slopes drier, and noise draws patches. Earth's land under the same dryness sets the mix, a guesstimate: plants cover 2.5% of the ground at an aridity of 0.02, 17% at 0.1, 46% at 0.2, 86% at 0.5 and 97% at 1 (`COVER`; Earth's drylands lose their cover fast below 0.2, Berdugo et al. 2020), trees take a growing share above 0.3 (`TREES`), and the grass is straw where it is dry and green where it is wet (`GREEN`).
- **Bare ground** is weathered soil. Under the air space weathering ends, and oxygen and water rust the rock's iron, so the bare ground takes the dry colours of Earth's soils on like rock under like dryness. The terrain's upkeep (register, 13 October) carries rock dust back uphill at 1–42 mm a century and renews the top metre every 2,400–100,000 years, so the topsoil stays moderately weathered. Each end of the rock's iron has a soil colour in Munsell notation at three aridities (`SOIL_MUNSELL`), converted exactly from the 1943 renotation's xyY (Bradford-adapted from Illuminant C to D65): on the feldspathic highland rock, whose 4–6% FeO matches Earth's upper continental crust, light brown 7.5YR 6/4 in the deserts, brown 7.5YR 5/4 in grassland and yellowish red 5YR 5/6 where it is humid; on basalt, reddish brown 5YR 4/4, dark reddish brown 5YR 3/4 and the dark red of terra roxa, 2.5YR 3/6, darkened 3% for each per cent of TiO2 above 2. The rock comes from two geography inputs, restored with `python -m geography.fetch_surface_inputs --download` (about 19 MB, hash-checked): the LROC team's mapped mare boundaries (geography/mare_map.py; Nelson et al. 2014), where basalt holds 7.6% of the design surface's land, most of it the dry plain of Oceanus Procellarum above the near-side sea, and Lunar Prospector's gamma-ray FeO and TiO2 (geography/surface_composition.py; Prettyman et al. 2006), whose 5° footprint gives the large-scale iron: the South Pole–Aitken floor at about 9.5% FeO, the open highlands at 3–5%, the near side's highland margins at 10–13%, partly the instrument's smear and partly the ejecta and mixing the Apollo 15 soils show. The maria take at least 15% FeO. The drylands' soils and living crusts are chosen with the biome design (register, 3 November), so these are the soils the rock and climate give unaided.
- **Water** takes the sea appearance study's colours for the designed waters, pi times their remote-sensing reflectance weighted by sunlight in each channel (appearance.water_albedos): the open sea's deep blue, the productive coast's green over water shallower than 150–20 m and in shallow lakes, and the river mouth's ochre in rivers and their plumes. A plume reaches about 0.3 km per (m³/s)^½ of the flow a river has gathered off land since the last lake or sea (appearance.land_discharge, on the design routing's receivers), since lakes and seas trap their rivers' sediment and send clear water on; the neck's 11,032 m³/s of spill is sea water. Of the 268,824 m³/s the rivers carry into the seas, 3,382 m³/s has gathered since the last lake, and no mouth into a sea carries more than 54 m³/s of it, so the seas' coasts show almost no plumes and the lakes the rivers enter hold them. The water also mirrors the sky: the shader adds Fresnel's reflectance of the diffuse skylight at the view's angle, 2% straight down.
- **Relief** shades the ground from the heights at the full 16 pixels per degree in place of 8, where the land's median slope is 4.0° against 3.4°. The diffuse skylight reaches each place in the share of an even sky its surroundings leave open (appearance.sky_view: the horizon out to 50 km in 16 directions, the Moon's curvature taken off), half the land seeing more than 99.7% of it and the steepest one per cent less than 94%; under this air the sky gives 43% of the green light on level ground with the Sun 60° up, so hollows and crater walls darken where they see less of it. The fully lit view's Sun stands 45° up, the cartographers' hill-shading light, in place of 60°, which roughly doubles the relief's light and shade on gentle slopes; the topography's shading takes the same light.

**Falls and cuts.** With `--falls` the page marks the twelve largest falls of the atlas round's falls study (research/studies/atlas_round/falls.py), each with a card giving its flow, drop, length, steepest 1.9 km step and the power of its falling water and what LOLA's 59 m heights show of the largest (long cascades, no step steeper than 45°), under a Falls and cascades checkbox shown by default in both modes; rivers turn to white water as their steps steepen from 30 to 200 m a km, and spray whitens up to about 5 km round the falls of 2 GW or more (`FOAM_SLOPE`, `MIST_GW`, `MIST_KM`, `MIST_MAX`), display guesstimates. The design surface now carries the works' cuts (design_cuts.py): the trench through Moscoviense's northern ridge is cut into the heights the topography and the relief draw, with the branches' 14,176 m³/s drawn through it as a river. At 1.9 km it is one cell wide, and the topography mode does not show it: its colours and contours average the heights over 7.6 km cells, where the trench lowers the deepest by 61 m beside neighbours 350 m apart on average, its raised relief takes points 21 km apart, and its 1.9 km shading draws a line no steeper than the ridge's own slopes. So a Works checkbox marks each piece of trench at its deepest point, with a card of its depth, length, water level and flow, and draws each cut's course as atlases draw canals: an amber line with a dark edge a few pixels wide at any zoom, in both modes, along the cells the cuts product lists in the order the water takes them (`course_cell`, from the edge of the water a stretch leaves over its passes to the edge of the water it falls into; render.work_lines), far wider than the trench at most zooms.

### The sea's slopes

With `--wind` the water mirrors the Sun and the sky through slopes that follow the design climate's daytime wind. The climate domain's surface wind (climate/gcm/surface_wind.py) gives each 170 km cell's medians over the run's daylit 3-day intervals of the friction velocity of the surface stress, from the run's own air density, and of the 10 m wind. On the run's open water the friction velocity is the stress's; over its land cells, which hold the lakes and the seas' edges at that grid, the land's roughness has slowed the 10 m wind, so the wind at a blending height of 100 m, carried up from the land's 10 m wind with its own friction velocity by the log law, comes down over the water through the neutral surface layer with Charnock's roughness (appearance.water_friction_velocity; on open water that carry gives about a tenth less than the stress). The illumination domain's model of a wind-roughened sea at lunar gravity turns the friction velocity into a mean square slope (illumination/water_surface/open_sea.py, open_sea_slopes): the short-wave regime of Elfouhaily et al. (1997) above 1 m waves, as the sea appearance study joins it to its wave runs, and in place of a wave run the long-wave regime of a fully developed sea, about a quarter steeper than the study's wave runs at its coasts' median winds. On the design climate's estimate (the run on the design surface, years 45–59, carried to 300.31 K) the open seas' median daytime friction velocity is 0.093 m/s and their mean square slope 0.0242, 0.017–0.030 across them, an rms tilt of 9°, as on the run at 0.930 (0.094 m/s, 0.0245); the page drew one slope of 0.04 for every water before. The map sits on the clouds' 2° grid (`sea_slope.png`, full scale 0.08). The shader takes Gaussian facet slopes, as Cox and Munk measured Earth's, for the glint: Fresnel's law for seawater (index 1.341 at 550 nm), the slopes' density at the facet that mirrors the Sun, and Smith's shadowing for both rays (Walter and others' rational fit). The sky the water mirrors comes from a table of rough water's reflectance of an even sky (illumination/water_surface/reflection.py, sky_reflection; `sky_reflectance_table`, 32 view angles crowded toward grazing by 8 rms slopes), in place of level water's reflectance: seen 3° above the horizon, level water mirrors 0.73 of the sky (the page's former formula, Schlick's, gave 0.78) and water of the open seas' median slopes 0.42. The run's 3-day means leave out gusts and calm hours, which would break the glitter into patches, and the slopes are isotropic, with no direction from the wind.

### Regions and districts

With `--districts` the page shows the atlas round's first draft of regions and districts (research/studies/atlas_round/districts.md). A texture in the globe's frame, `districts.png` at 8 pixels per degree, holds each cell's district (red and green, 0 on open water) and region (blue), read cell by cell. Both modes draw the regions' edges, and zoomed in past 2.6 Moon radii from the centre, the districts' edges as well, where a neighbour a pixel or so away differs. Labels follow the same switch: the regions' names when zoomed out, the districts' and the coastal seas' names when zoomed in. A Regions and districts box shows or hides them. A label's card gives the name, code, region, tongue, area and short description. Selecting a coastal sea tints it with the land districts it faces and whose rivers reach it, and selecting a district tints it with its coastal seas: the hybrid the author chose so that a coast and its land are planned together. The old sea labels, Latin working names, are now off by default in the topography mode too. A page section explains the draft and lists its tongues and the Terlunan words coined for it.

### Peeling the air with the zoom

The Air in view control has a With the zoom box, ticked by default: the zoom stands in for a descent through the air, and from the appearance's opening view (5.8 Moon radii from the centre) to the nearest zoom (1.35) the share of the air's column left in view falls from all of it to none, evenly with the logarithm of the camera's height above the ground; the slider's thumb follows to the nearest stop and the readout gives the height. Moving the slider takes the peel off the zoom, and `peel<km>` in a link opens at a stop with the box clear. The camera itself stays above the air; one that truly descended would keep most of the haze until the last few tens of kilometres, where the ground this page holds is a blur, so the zoom stands in.

The tests (visualization/atlas/test_appearance.py, geography/tests/test_mare_map.py, test_surface_composition.py) check the Munsell conversion and the soils at their anchors, the cover curve, the sediment flow on a small graph with a lake and a spill put in, the plumes and the sky on level ground and in a pit, the mare map's share and known places, and Lunar Prospector's oxides at known places.

## Faithfulness

`manifest.json` records the SHA-256 of every product the sheets read, and a round-trip check. The water share of the
rendered near-side disk must match the product's Earth-facing disk share within 0.01; at 28% both are 0.415. The relief
texture, decoded as the page decodes it, must give back every cell's height within half a metre and on its own side of
sea level. Its water and the 16 pixels-per-degree surface must each cover the product's sea share within 0.002 by area
(at 28% both are 0.280). The surface's rain-fed lakes must cover the drainage product's lake share within 0.002
(0.1034), and the water mask, which adds river channels, may exceed sea and lakes by at most 0.004 (0.3836).

The sheets display the product's hydrostatic geometry, so every depression below sea level appears as water. Shorelines and lake levels change once the hydrology is added.
