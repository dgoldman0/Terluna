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

## Engine capability (any environment)

1. **An explicit renderer interface.** Replace the transitional `OM` registry with
   direct imports and a material-factory interface that the WebGL path and the
   WebGPU/TSL lab both implement.
2. **Separate engine from world in the shoreline code.** Waypoints, the shelter,
   the authored weather episode and the community placement rules move to
   `world/`; geometry, rendering and simulation stay in `engine/`.
3. **Scale.** Stream terrain tiles from elevation products; instance vegetation
   with levels of detail and impostors, so forests of tall trees reach the
   horizon; add atmospheric perspective over tens of kilometres.
4. **Night and time.** Mesopic colour at low light, a continuous multi-month
   clock (so libration and the stars do not jump at the month boundary), and
   earthlit sky glow with earthlight's own spectrum once illumination computes it.
5. **Low gravity.** Motion and physics that feel lunar: long jumps and falls,
   gliding, slow rain and waves.
6. **Sound.** Ambience that follows wind, water, rain and distance.
7. **Performance.** Frame times measured on real GPUs, with budgets per system.

## World system (fed by the domains)

`world/` describes an environment as data drawn from domain products:

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
