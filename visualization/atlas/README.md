# Atlas sheets

Map sheets of the Open Moon at its selected standing-water share. They draw the geography [atlas](../../geography/README.md#atlas-at-the-selected-water-share) product and the heritage sites and scientific targets of the [conservation study](../../research/studies/conservation/README.md).

```sh
python -m geography.atlas
python -m research.studies.conservation.run
python visualization/atlas/render.py      # out/: three sheets, globe.html and manifest.json (kept out of Git)
```

| Sheet | Shows |
|---|---|
| `atlas_<pct>pct_near_far.png` | The near side as seen from Earth and the far side, in orthographic projection |
| `atlas_<pct>pct_global.png` | The whole Moon in longitude and latitude, with a 30° grid |
| `atlas_<pct>pct_poles.png` | Both polar regions poleward of 80°, from the product's 16 px/deg polar caps, with the permanently shadowed craters |
| `globe.html` | An interactive 3D globe in one self-contained page (three.js from cdnjs, textures embedded): views of the near side, far side and poles, adjustable relief with the seas flat at sea level, labels, and a record card for each heritage site and target |

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

[globe.template.html](globe.template.html) holds the viewer, and `render.py` fills it with three textures built from the same grid:

- the colour map, with water in depth bands and land in a neutral placeholder grey;
- a tangent-space normal map of the relief;
- the sea-clamped surface as displacement.

The light follows the viewer, so every face reads as a map. The page works from the file itself and was also published privately to claude.ai.

## Faithfulness

`manifest.json` records the SHA-256 of every product the sheets read, and a round-trip check. The water share of the rendered near-side disk must match the product's Earth-facing disk share within 0.01; at 28% both are 0.415. The globe's colour texture must cover the product's water share within 0.002 by area; at 28% both are 0.280.

The sheets display the product's hydrostatic geometry, so every depression below sea level appears as water. Shorelines and lake levels change once the hydrology is added.
