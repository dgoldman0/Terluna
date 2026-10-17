# Jura Island flyover

A 44-second flight over Jura Island and a map of its topography and water depth, both a common-water-level view:
LOLA terrain at 16 px/deg with every basin below the atlas's sea level filled to that one level. The flight crosses
the Sea of Rains into the Bay of Rainbows, runs along the Jura front and over the western uplands, and climbs to see
the whole island. Each of its five chapters holds one fixed panel naming two reference points, the summit's height or
a sounding's depth, with their places marked A and B on a route map. The map adds depth contours, the soundings, a
profile along 44.5° N and a 3-D inset of the land and seafloor.

The view reads the [Jura Island study](../../research/studies/jura_island/README.md)'s product and the geography
[atlas](../../geography/README.md#atlas-at-the-selected-water-share) it was measured against:

- `research/studies/jura_island/results/shores.json` (`terluna.research.jura-shores/1`) gives the common water level,
  −1,654.2 m relative to the geoid, the island's summit and the input hashes the study measured;
- `geography/results/atlas.json` (`terluna.geography.atlas/1`) gives the water share, 28% of the Moon, and must be the
  atlas whose sha256 the island study recorded;
- LOLA's LDEM_16 and LDEM_4 and the GRAIL GL0420A coefficients, through
  [geography's inputs](../../geography/fetch_inputs.py), must match `geography/inputs.json` and the hashes in the
  island study's product.

[terrain.py](terrain.py) regenerates the terrain from those inputs in about a second, as the island study does: LOLA's
radii at their native pixel centres less the degree-200 geoid evaluated on the same centres, less the common water
level. The island is the connected land holding the product's summit. Its area, 57,824 km², its summit, 2,717 m at
44.41° N, 39.28° W, and its extent equal the island study's. Any difference in a product or an input stops the build,
and no array is kept in Git.

```sh
python3 -m pytest visualization/jura-island/tests                                    # no OpenGL needed
python3 visualization/jura-island/render_map.py                                      # out/jura-island-topography.png
python3 visualization/jura-island/render_flyover.py --frames 0 22 43.9 --scale 0.5   # stills at 800 x 450
python3 visualization/jura-island/render_flyover.py --preview                        # eight stills and a contact sheet
python3 visualization/jura-island/render_flyover.py                                  # the film, 2,640 frames
```

Outputs go to `out/` (kept out of Git) or to `--out`, each with a record (`render-record.json`, `map-record.json`)
listing the products, inputs and code read, by sha256, with the camera's checks and the soundings. The flyover
needs moderngl and glcontext ([requirements.txt](requirements.txt)) and an EGL driver; the film also needs ffmpeg.
The map needs only the repository's packages. Both use the DejaVu Sans fonts.

## What the view shows

- **Computed:** the land and seafloor heights relative to the common water level, the water's extent and depth, the
  island's outline and summit, and the soundings, sampled bilinearly and rounded to 10 m (the terrain's spacing is
  about 1.9 km north–south).
- **Display:** the land stands at true relief on the 1,737.4 km reference sphere and the water surface is the sphere
  itself, a mean surface with the geoid-relative relief draped on it. Depth bands and 250 m contours overlay the
  water. The light comes from a fixed direction in the north-west, about 21° up, a cartographer's convention; the haze
  is a display formula; the land's colours follow height alone. The map's hillshade is lit from the north-west too.
  The map's 3-D inset exaggerates relief 20 times and its profile's axes have independent scales.
- **Left out:** lake levels apart from the common one, the monthly tide, erosion and crustal loading. Travel time is
  compressed; altitudes and distances are true, the camera never lower than 16.9 km above the water.

## The original

The flyover was made in ChatGPT on 15 October 2026 and recovered on 17 October. Its film, its map and its source zip
are catalogued with their sha256 on the research drive, in
`/media/projectspace/terluna-research/jura-island/flyover/` (`catalogue.json` there), and stay out of Git:

| File | Bytes | sha256 |
|---|---|---|
| `jura-island-flyover.mp4`, 1600 × 900, 60 fps, 2,640 frames | 18,501,007 | `c5dc1a9c1ca35acb4f0ba4ff2b30a8885f2b71dfdfa7aee2044f88ad9be7a952` |
| `jura-island-topography.png`, 3240 × 2160 | 2,564,705 | `8fd61ddb7198f79b7d1c775bd5a34056be7fb5b19365181cad1338361301bfa8` |
| `jura-island-flyover-source.zip`, ten members | 12,120,259 | `981b6e6f2a9660ebb24372be7f76a4ae34dced69333e4751040e1a89201d567b` |

[repository_import.json](repository_import.json) records each member's hash and what became of it, and keeps the
original's settings, checks, chapters and soundings, without its ChatGPT Library identifier. The tests hold the port
to that record: the same lowest altitude (16.924 km) and height above the ground (16.768 km), the same largest turn
of the camera (0.435° a frame), four panel changes with 0.75 s fades and holds of at least 6.5 s, and the seven
soundings within 0.5 m and identical when rounded.

### How this port departs from the original

- **Terrain from the inputs.** The saved `terrain.npz` gives way to terrain.py. The geoid is evaluated at the 16 px/deg
  pixel centres, as the island study evaluates it, where the original interpolated it from 4 px/deg; heights move by
  at most 0.49 m, and the island keeps every one of its 22,869 pixels.
- **The water level from the island study.** The common water level is the product's −1,654.2 m, the atlas's value,
  in place of the −1,654.16 m saved with the original arrays and the hard-coded text; heights move by 0.037 m.
- **The water share from the atlas,** checked against the hash the island study recorded, in place of the hard-coded
  28%; the summit's place and height, on the panels and the map, come from the island study's product.
- **The view named.** The frame and the map call it a common-water-level view with illustrative lighting and haze.
- **Depth colours on their pixels.** The depth texture's texel centres sit on the terrain's pixel centres. The
  original mapped the crop's first and last pixel centres to the texture's edges, which misplaced the depth colours
  by up to half a pixel towards the crop's edges and by almost nothing near the island.
- **Smaller changes.** The 2 px/deg surround closes its seam at 0° E; the band edges and the contour interval reach
  the shader from [water_depth.py](water_depth.py); frames render at a fraction of full size with `--scale`; each
  output carries its record.

The camera's keys and splines, the unwrapped heading and pitch, the chapters, their panels and fades, the route map
and its position marker are the original's. Frames rendered at full size at 0, 22 and 43.9 s differ from the original
film's by a mean of 1.4–2.9 of 255 levels over the scene, the larger differences along band edges, contours and
coasts, and the captions by design. Test frames at 800 × 450 and the rendered map are kept in the catalogue folder.

## Files

| File | Holds |
|---|---|
| [film.py](film.py) | The frame size and rate, the camera path, the chapters and fades, and the frame's texts; no OpenGL |
| [terrain.py](terrain.py) | The products' and inputs' checks, the regenerated terrain, the island and the soundings |
| [water_depth.py](water_depth.py) | The depth bands, their colours, the contour interval and the named points |
| [render_flyover.py](render_flyover.py) | The scene, its overlays, stills, contact sheets and the film |
| [render_map.py](render_map.py) | The topography and water-depth map |
| [repository_import.json](repository_import.json) | The recovered source, its outputs and its record |
| [tests/test_flyover.py](tests/test_flyover.py) | The camera, the panels, the inputs' hashes, the water level and the terrain against the products |
