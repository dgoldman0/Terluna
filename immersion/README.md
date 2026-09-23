# Immersion

Interactive experiences of what an Open Moon could look and feel like. They are
illustrations built on the project's research, not evidence for it: the
research domains own the science, and the experience reads their results.

## Layout

| Folder | Holds |
|---|---|
| [engine/](engine/) | Runtime systems any world can use: sky, Earth, stars and clouds drawn from baked products, clipmap terrain and surface materials, sea, rain and surface water, procedural vegetation geometry, the far forest drawn as impostors, audio. The engine never imports a world; it reads the active one from `OM.world`. |
| [world/](world/) | World definitions. [cove.js](world/cove.js) is the development cove: its landscape field, viewpoints, rain shelter and path, authored weather episode, and planting rules, all placeholders. `cove-forest.worker.js` places its distant trees off the main thread. |
| [experiences/](experiences/) | Entry points: [shoreline](experiences/shoreline/) (the walkable development scene) and [renderer-lab](experiences/renderer-lab/) (the WebGPU/TSL experiment) |
| [bake/](bake/) | Turns domain products into runtime assets: `sky_atlas.py` packs the [illumination/sky](../illumination/sky/) atlases into `assets/sky/`; `columns.mjs` packs the atmosphere domain's column set; `sky_bodies.py` bakes the site ephemeris, the bright stars and the Earth-disk calibration; `surfaces.py` generates the surface textures |
| [assets/](assets/) | Baked runtime assets and their manifests |
| [tests/unit/](tests/unit/) | Node tests of engine and world modules |
| [docs/](docs/) | The [roadmap](docs/roadmap.md), rendering notes, the renderer roadmap, and the [checkpoint history](docs/history/) |

## Run

```sh
cd immersion
npm install
python -m pip install -r bake/requirements.txt
npm run bake                # column set, site sky and stars, surface textures
npm run dev                 # then open /experiences/shoreline/
npm test                    # unit tests
npm run build               # static site in dist/; `npm run preview` serves it
```

`dev`, `build` and `test` re-run the fast bake steps first. The surface bake
skips work when the textures already match their manifest.

The baked sky (`assets/sky/atmosphere.json`) is committed, because regenerating it
needs the numba sky solve. After the illumination atlases change, re-bake it with
`python bake/sky_atlas.py`. The experiences need a WebGL2 browser.

## How the experience relates to the research

- **Science stays in its domain.** The clear-sky solver lives in
  [illumination/sky](../illumination/sky/), the column model in
  [atmosphere/column](../atmosphere/column/), and the light-transport references
  that measure this renderer in [illumination/references](../illumination/references/).
  Tools that render science faithfully are in [visualization](../visualization/).
- **The experience reads results.** `bake/` turns domain products into assets:
  the sky from the illumination atlases, and the clouds' atmospheric columns from
  the atmosphere domain's column set. At runtime the experience interpolates
  those data and runs no domain model.
- **The scene is a placeholder.** The shoreline cove, its trees and its weather
  episode were authored for development. They make no claim about what Open Moon
  environments will be. Specific environments should come from the biology,
  ecology, geology and hydrology research as it matures.

## The land to the horizon

- **Trees** stand on an 8 m lattice by the world's planting rule (for the cove, its
  authored habitat rule: artistic placeholder). The cells around the viewpoints are
  built as full tree models; every other cell within range is drawn by
  `engine/forest.js` as a camera-facing impostor card. The cards are images of the
  same tree models from eight directions, with a random yaw per tree, baked on the
  GPU at start-up (`engine/impostors.js`). A worker places the trees, nearest
  first, so distant woodland fills in over a few seconds without stalling frames.
- **Range** follows the quality setting: 1.8 km (economy), 3 km (balanced) or
  4.5 km (high). Cards thin out toward the limit. Beyond it the terrain shader
  carries the woodland: the world's expected canopy cover per terrain vertex, seen
  along the view ray (more ground hidden at grazing angles). Its strength and
  colour were calibrated against the cards' mean colour at the same distances.
- **Perceptual conventions:** the cards and the far woodland shading stand in for
  trees too small to model individually. Crown occlusion and rounded crown normals
  are baked into the cards; distant trees cast no shadows. On a small world the
  horizon is close (about 2.4 km for a standing person, 8 km from the overlook),
  so higher ground beyond it still shows.

## The sky

- **The Sun, Earth and stars** come from the illumination domain's site ephemeris
  for the placeholder development site (0° N, 65° W). Earth hangs in the east,
  17–33° up as it wobbles with libration, and goes through phases opposite to the
  lunar day: a thin crescent near dawn, nearly full near sunset, 72% lit at midnight.
  Its disk is drawn from NASA Blue Marble imagery. Seen through the thick Open Moon
  air, its light is reddened.
- **Earthlight** lights the land at night: it becomes the shadow-casting key light
  once the Sun is down, and its scattered glow is added to the sky. The glow
  reuses the solar clear-sky atlas at Earth's elevation, so it assumes earthlight
  is sunlight-coloured.
- **Stars** come from the Yale Bright Star Catalogue, dimmed by the atmosphere's
  transmission, cloud and fog.
- **Colour** is shown as calculated: the atlas's spectral radiance converted to
  sRGB, whose white is average Earth daylight (D65). This is how the place would
  look to a daylight-balanced camera or to eyes that just arrived from Earth. The
  sky solver's proxy atmosphere (an exponential profile without haze, not yet the
  atmosphere domain's solved columns) holds about seven times Earth's column of
  air. That reddens the noon Sun to about 3,500 K and turns the sky pale and the
  horizon cream. **Adapted to local daylight** instead scales
  each channel so local noon light renders neutral, roughly as a resident's eyes
  would adapt; that makes the sky look more like Earth's than it is.
- **Perceptual conventions:** the adapted exposure follows light as L^−0.85 (eye-like,
  so earthlit nights are dim but visible). It meters the rendered sky itself,
  cloud included, plus the direct Sun and Earth. It brightens within about a
  second, darkens over a few seconds, and resets at once when you jump to another
  place, time or weather. The Earth disk is compressed for display so its surface
  stays readable, and stars are drawn at the eye's resolution (about one
  arcminute) rather than the screen's. The fixed exposure and the M1 noon
  reference (which keeps the local-daylight colour balance) are unchanged.

Past checkpoint methods, handoffs and validation records are in
[docs/history](docs/history/).
