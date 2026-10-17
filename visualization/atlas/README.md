# Atlas sheets

Map sheets of the Open Moon at its selected standing-water share. They draw the geography [atlas](../../geography/README.md#atlas-at-the-selected-water-share) product and the heritage sites and scientific targets of the [conservation study](../../research/studies/conservation/README.md).

```sh
python -m geography.atlas
python -m research.studies.conservation.run
climate/gcm/.venv/bin/python -m climate.gcm.climatology A28_corrected_970:10-19   # the design case's clouds and rain
python visualization/atlas/render.py      # out/: three sheets, globe.html and manifest.json (kept out of Git)
```

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
per degree (1.9 km), with relief shading at 8. The topography takes its heights from the atlas's 4 and its
shading from the same relief at 8. JPEG keeps the page near 12 MB.

On 16 October the globe took its clouds and ground cover from the corrected design case in place of run A, and
the lakes as the drainage product holds them, 10.3% of the Moon in place of the first estimate's 12%, and was
republished to its private claude.ai page. The lakes follow the design case once it settles, after the glow pair
([register](../../research/decisions.md#the-globe-and-the-glow-pair-16-october-2026)). The same day the
topography replaced the map mode, as a mode of its own and as an overlay on the appearance, and the globe gained
the fully lit option, the grid's spacings, going to a place and the readout of a clicked place.

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
  - surface cover, from rainfall, soil moisture, nearness to water, floodplains, height, hollows and slope;
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
- `red`, `noclouds` or `noair`.

The page was also published privately to claude.ai.

## Faithfulness

`manifest.json` records the SHA-256 of every product the sheets read, and a round-trip check. The water share of the
rendered near-side disk must match the product's Earth-facing disk share within 0.01; at 28% both are 0.415. The relief
texture, decoded as the page decodes it, must give back every cell's height within half a metre and on its own side of
sea level. Its water and the 16 pixels-per-degree surface must each cover the product's sea share within 0.002 by area
(at 28% both are 0.280). The surface's rain-fed lakes must cover the drainage product's lake share within 0.002
(0.1034), and the water mask, which adds river channels, may exceed sea and lakes by at most 0.004 (0.3836).

The sheets display the product's hydrostatic geometry, so every depression below sea level appears as water. Shorelines and lake levels change once the hydrology is added.
