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
| [bake/](bake/) | Turns domain products into runtime assets: `sky_atlas.py` packs the [illumination/sky](../illumination/sky/) atlases into `assets/sky/`; `surfaces.py` generates the surface textures |
| [assets/](assets/) | Baked runtime assets and their manifests |
| [tests/unit/](tests/unit/) | Node tests of engine and world modules |
| [docs/](docs/) | Rendering notes, the renderer roadmap, and the [checkpoint history](docs/history/) |

## Run

```sh
cd immersion
npm install
python -m pip install -r bake/requirements.txt
python bake/surfaces.py     # surface textures, checked against assets/surfaces/manifest.json
npm run dev                 # then open /experiences/shoreline/
npm test                    # unit tests
npm run build               # static site in dist/; `npm run preview` serves it
```

The baked sky (`assets/sky/atmosphere.json`) is committed, because regenerating it
needs the numba sky solve. After the illumination atlases change, re-bake it with
`python bake/sky_atlas.py`. The experiences need a WebGL2 browser.

## How the experience relates to the research

- **Science stays in its domain.** The clear-sky solver lives in
  [illumination/sky](../illumination/sky/), the column model in
  [atmosphere/column](../atmosphere/column/), and the light-transport references
  that measure this renderer in [illumination/references](../illumination/references/).
  Tools that render science faithfully are in [visualization](../visualization/).
- **The experience reads results.** The sky comes from a baked atlas. One
  exception remains: the cloud renderer still imports the atmosphere domain's
  column model directly. Baked column products will replace that import.
- **The scene is a placeholder.** The shoreline cove, its trees and its weather
  episode were authored for development. They make no claim about what Open Moon
  environments will be. Specific environments should come from the biology,
  ecology, geology and hydrology research as it matures.

Past checkpoint methods, handoffs and validation records are in
[docs/history](docs/history/).
