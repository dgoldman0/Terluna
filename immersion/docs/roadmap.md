# Immersion roadmap

The experience grows in three layers. Work on engine capability that holds for
any environment comes first. Specific environments wait for the domain research
to supply them.

## Done

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

## Engine capability (any environment)

1. **An explicit renderer interface.** Replace the transitional `OM` registry with
   direct imports and a material-factory interface that the WebGL path and the
   WebGPU/TSL lab both implement.
2. **Scale.** Stream terrain tiles from elevation products (geography has none
   yet). Add a middle level of 3D trees so that walking out of the detailed
   planting does not bring impostor cards close; let distant trees darken the
   ground with their shade; take the aerial-perspective extinction from the
   illumination domain instead of the engine's band constants.
3. **Night and time.** Mesopic colour at low light, a continuous multi-month
   clock (so libration and the stars do not jump at the month boundary), and
   earthlit sky glow with earthlight's own spectrum once illumination computes it.
4. **Low gravity.** Motion and physics that feel lunar: long jumps and falls,
   gliding, slow rain and waves.
5. **Sound.** Ambience that follows wind, water, rain and distance.
6. **Performance.** Frame times measured on real GPUs, with budgets per system. On
   an Intel UHD (TGL GT1) iGPU the shoreline takes roughly 85–150 ms a frame at
   1280×720 before the far forest (which adds 0–17%); terrain shading is the
   largest share, then vegetation, then the sky.

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
