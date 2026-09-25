# Atlas sheets

Map sheets of the Open Moon at its selected standing-water share. They draw the geography [atlas](../../geography/README.md#atlas-at-the-selected-water-share) product and the heritage sites and scientific targets of the [conservation study](../../research/studies/conservation/README.md).

```sh
python -m geography.atlas
python -m research.studies.conservation.run
climate/gcm/.venv/bin/python -m climate.gcm.climatology A:30-39   # the cloud and rainfall climatology
python visualization/atlas/render.py      # out/: three sheets, globe.html and manifest.json (kept out of Git)
```

| Sheet | Shows |
|---|---|
| `atlas_<pct>pct_near_far.png` | The near side as seen from Earth and the far side, in orthographic projection |
| `atlas_<pct>pct_global.png` | The whole Moon in longitude and latitude, with a 30° grid |
| `atlas_<pct>pct_poles.png` | Both polar regions poleward of 80°, from the product's 16 px/deg polar caps, with the permanently shadowed craters |
| `globe.html` | An interactive 3D globe in one self-contained page (three.js from cdnjs, textures embedded), with an appearance mode and a map mode: views of the near side, far side, poles and the Moon from Earth, the Sun set by the Moon's age, labels, and a record card for each heritage site and target |

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
[appearance.py](appearance.py) and with map textures built from the atlas grid. The appearance mode's surface is
16 pixels per degree (1.9 km), with relief at 8. The map mode stays at the atlas's 4. JPEG keeps the page
near 9 MB.

**Appearance mode** estimates how the Open Moon looks from outside its air. One full-screen pass traces each
pixel's ray through the atmosphere, two cloud decks and the surface:

- **Computed:**
  - land, sea and relief: LOLA above the GRAIL geoid at 16 pixels per degree through the geography domain,
    filled to the atlas product's sea level;
  - rivers and rain-fed lakes above sea level, from the geography domain's drainage product (12% of the Moon
    in lakes); river width follows discharge;
  - the air, from the illumination domain's engine atmosphere: optical depth and scale height;
  - the halo around the limb and the sunlit air beyond the terminator, which come out of the same ray-march
    through the air's full 484 km;
  - the ground's direct and diffuse light at every Sun height, from the Open Moon sky atlas.
- **Informed:**
  - the low and high cloud decks, from the climate domain's GCM climatology. Each moves with the Sun
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

Refraction, calibrated multiple scattering along the limb and in twilight, Earthshine and clouds with depth
are the next refinements.

**Map mode** colours water in depth bands, shows land as a neutral placeholder grey, and lets the light follow
the viewer so that every face reads as a map. It uses three textures: the colour map, a tangent-space normal
map of the relief and the sea-clamped surface as displacement.

The page opens in the view named by its link: `#near`, `#far`, `#south`, `#north`, `#earth` or `#map`. It
can add `e<degrees>` for the Sun's elongation from Earth (`#earth.e60` is a crescent), `x<stops>` for
exposure, `z<tenths>` for the camera's distance in Moon radii (from Earth, the field of view in tenths of a degree), `lat<degrees>n|s` and `lon<degrees>e|w` to centre a place, `red`, `noclouds` or `noair`. The page was also published privately to claude.ai.

## Faithfulness

`manifest.json` records the SHA-256 of every product the sheets read, and a round-trip check. The water share of the rendered near-side disk must match the product's Earth-facing disk share within 0.01; at 28% both are 0.415. The map's colour texture and the 16 pixels-per-degree surface must each cover the product's sea share within 0.002 by area (at 28% both are 0.280); the surface's rain-fed lakes must cover the drainage product's lake share within 0.002 (0.120), and the water mask, which adds river channels, may exceed sea and lakes by at most 0.004 (0.4005).

The sheets display the product's hydrostatic geometry, so every depression below sea level appears as water. Shorelines and lake levels change once the hydrology is added.
