# Working on the immersion

The immersion is the explorable world: what an Open Moon could look and feel
like to be in. It is judged as an experience, by whether it is believable,
coherent with the research and pleasant to explore. Its images are illustrations.
Vividness never upgrades the evidence behind them.

## What belongs here

- `engine/`: rendering and runtime systems that any environment can use. The
  engine never imports a world; it reads the active one from `OM.world`.
- `world/`: world definitions, meaning landscapes, communities, structures, sites
  and weather. Each world is one object, like `world/cove.js`, that an experience
  registers on `OM.world`.
- `experiences/`: entry points that combine engine and world.
- `bake/`: turning domain products into runtime assets.

What does not belong here: radiative-transfer benchmarks, accuracy programs,
reference renderers, validation campaigns, or any model whose outputs someone
might cite. Those go to their domain (for example `illumination/references`) or
to `visualization/`. If the experience needs a quantity no domain provides, add it
to the domain with tests, then read its product.

## Where content comes from

Environments, biomes, species and scenes come from the domain research as it
matures: terrain and water from geography, climate from climate, living
communities from the biosphere, structures from habitation. Do not design
environments from the manuscripts' narrative passages; they presuppose what the
environment is like. Until the research can supply an element, use a neutral
placeholder and label it as one. The shoreline cove is a development scene.

## Reading the science

- Runtime code imports only `immersion/` and `shared/`, and engine code never
  imports `world/`; `bake/` may use domain code and products.
  `shared/check_layers.py` enforces both.
- Read products by their stated rule (interpolation, units) and record which
  product version a baked asset came from.
- Take constants from `shared/constants.js`.

## Saying what is what

Every element falls into one of these kinds. Say which in code comments, the
README and the in-app notes:

- **computed**: taken straight from a domain product (the clear-sky atlas,
  earthlight, column clouds);
- **informed**: shaped by research but adapted for display (Earth-disk
  calibration, tree forms scaled from biosphere allometry);
- **artistic**: no research basis yet (the cove's layout, placeholder weather);
- **perceptual convention**: a choice about how things are shown (adapted
  exposure, star brightness at eye resolution).

## Practice

- Write readable ES modules. Prettier formats them (`npm run format`).
- `npm test` runs the unit tests; `npm run build` must succeed.
- For refactors, compare deterministic renders before and after. A refactor
  should be pixel-identical; say so when a visual change is intended.
- Show progress with a playable build and screenshots, and measure frame times
  on real GPUs when available. Keep documentation to one README per experience
  plus the roadmap. Do not accumulate handoff and validation-record layers; Git
  history keeps the past.
