# Shoreline immersion prototype — WIP checkpoint

This directory preserves the current first-person shoreline prototype before a visual-quality rework.

The prototype adds traversable shoreline geometry, a woodland path, rain shelter and overlook, shared weather state, persistent wetness, spatial audio, water, and Three.js integration around the existing Open Moon light-cycle data.

## Current status

This is a preservation checkpoint, **not a visual-quality milestone**. The current rendered scene has serious known weaknesses observed in-browser:

- distant terrain can read as floating or disconnected ribbons instead of coherent shoreline relief;
- water reflections are coarse and visibly aliased;
- vegetation, rocks, ground materials and their distributions remain placeholder-like;
- lighting, exposure and atmospheric display can become muddy enough to obscure scene structure;
- the scene has not yet been calibrated to a lifelike visual standard.

The next milestone should focus on one convincing shoreline view at noon, then the same view under low Sun and twilight, before adding more environmental complexity.

## Preserved source

`shoreline-code-wip.tgz` contains the editable implementation source, template, scenario, build scripts, and atmosphere repacking tool from this checkpoint. It excludes the large generated standalone HTML and packed runtime atmosphere data.

Reference hashes for the complete generated artifacts produced from this checkpoint:

- `Open_Moon_Shoreline.html` — SHA-256 `8cecea048a67040c532f1660b4b071dd0fe1db1927dfc063aba862e45ca3b218`
- complete `Open_Moon_Shoreline_Source.zip` — SHA-256 `607d8917aa584b1c49504b7ce466cb4ad66729f3bf89eeaa4fad550470890219`

The preceding full-cycle viewer remains unchanged.
