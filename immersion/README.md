# Immersion

Interactive experiences of what an Open Moon could look and feel like. They are
illustrations built on the project's research, not evidence for it: the
research domains own the science, and the experience reads their results.

## Layout

| Folder | Holds |
|---|---|
| [engine/](engine/) | Runtime systems any experience can use: sky and clouds drawn from baked atlases, clipmap terrain and surface materials, sea, rain and surface water, procedural vegetation geometry, audio |
| [world/](world/) | World content. At present this is only the authored development landscape. |
| [experiences/](experiences/) | Entry points: [shoreline](experiences/shoreline/) (the walkable development scene) and [renderer-lab](experiences/renderer-lab/) (the WebGPU/TSL experiment) |
| [bake/](bake/) | Turns domain products into runtime assets: `sky_atlas.py` packs the [illumination/sky](../illumination/sky/) atlases into `assets/sky/`; `columns.mjs` packs the atmosphere domain's column set; `sky_bodies.py` bakes the site ephemeris, the bright stars and the Earth-disk calibration; `surfaces.py` generates the surface textures |
| [assets/](assets/) | Baked runtime assets and their manifests |
| [tests/unit/](tests/unit/) | Node tests of engine and world modules |
| [docs/](docs/) | Rendering notes, the renderer roadmap, and the [checkpoint history](docs/history/) |

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

## The sky

- **The Sun, Earth and stars** come from the illumination domain's site ephemeris
  for the placeholder development site (0° N, 65° W). Earth hangs about 25° up in
  the east, wobbles with libration, and goes through phases opposite to the lunar
  day: a thin crescent near dawn, nearly full near sunset, 72% lit at midnight.
  Its disk is drawn from NASA Blue Marble imagery. Seen through the thick Open Moon
  air, its light is reddened.
- **Earthlight** lights the land at night: it becomes the shadow-casting key light
  once the Sun is down, and its scattered glow is added to the sky. The glow
  reuses the solar clear-sky atlas at Earth's elevation, so it assumes earthlight
  is sunlight-coloured.
- **Stars** come from the Yale Bright Star Catalogue, dimmed by the atmosphere's
  transmission, cloud and fog.
- **Perceptual conventions:** the adapted exposure follows light as L^−0.85 (eye-like,
  so earthlit nights are dim but visible), the Earth disk is compressed for display
  so its surface stays readable, and stars are drawn at the eye's resolution (about
  one arcminute) rather than the screen's. The fixed exposure and the M1 noon
  reference are unchanged.

Past checkpoint methods, handoffs and validation records are in
[docs/history](docs/history/).
