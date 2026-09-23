# Immersion roadmap

The explorable world is moving to Unreal Engine 5, in its own repository
([dgoldman0/terluna-game](https://github.com/dgoldman0/terluna-game)). This lane
keeps the research side: the world definitions and the exports the game reads
(`bake/game/`). The web engine and the shoreline viewer are frozen: they keep
working as a preview and get fixes only if they break. Specific environments still
wait for the domain research to supply them.

## Next

1. **Unreal test level** (in the game repository): the exported landscape with a
   curved horizon; the sky drawn from the atlas, with Unreal's atmosphere set from
   `engine_atmosphere.json` for haze and cloud light and its mismatch measured; the
   Sun and Earth from the ephemeris; the stars; one volumetric cloud layer from a
   column preset; lunar gravity.
2. **Retire the web engine** once that level shows the sky, Earth and stars,
   terrain, trees and clouds at least as well: move `engine/`, `experiences/`, the
   Vite build and their tests to `archive/`, rename this lane (for example to
   `world/`), and update the lane rules and guidance.
3. **Keep the exports the only path** from research to game, adding what the game
   needs as the domains provide it.

## Done (the web engine, now frozen)

- **Foundation:** readable ES modules built with Vite; Three.js from npm; unit
  tests; deterministic render comparisons for refactors.
- **Reading science as products:** the clear-sky atlas (illumination) and the
  atmospheric column set (atmosphere) are baked; no domain code runs in the
  experience.
- **The lunar sky:** Sun, Earth and stars from the illumination domain's site
  ephemeris, Earth's phases and libration, earthlight on the land and in the sky,
  and 9,096 catalogue stars. Adapted exposure makes earthlit nights visible.
- **Vegetation to the horizon:** a world gives its planting rule per lattice cell;
  the engine draws cells beyond the detailed planting as eight-view impostors
  (streamed from a worker) and carries woodland on the distant terrain beyond
  that. The column study moved out to `visualization/atmospheric-columns`.
- **Engine and world apart:** the cove's landscape, sites, shelter, path, weather
  episode and planting rules form one world object (`world/cove.js`) that the
  experience registers; the engine reads it and never imports it.

Engine capability (time and night, low-gravity motion, sound, performance) now
belongs to the game's own roadmap.

## World system (fed by the domains)

A world is an object like `world/cove.js`: a landscape field, sites, weather and
flora (and later structures). Future worlds draw these from domain products:

- terrain and water from geography (elevation, basin filling, drainage);
- climate from climate (temperature, wind, cloud and rain over the lunar day);
- living communities from the biosphere (which organisms, where, at what size,
  including tree forms from its mechanics results);
- structures and settlement from habitation.

Where a product does not exist yet, the world uses a neutral placeholder and
labels it as one. The shoreline cove remains the development scene.

## Environments (when the research supports them)

Detailed environments and any opening scene are chosen as the biology, ecology,
geology and hydrology research matures. They are not designed from the
manuscripts' narrative passages, which presuppose the environment.

## Open decision

Lunar gravity: the sky solver, the column model and the experience use a rounded
1.62 m/s², while other models use GM/R² = 1.6242 m/s² (`shared/constants.json`).
Adopting GM/R² everywhere means regenerating the sky atlas, the column set and the
A1–A3 references.
