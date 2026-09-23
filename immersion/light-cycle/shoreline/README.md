# Shoreline (development scene)

A walkable Three.js r186 scene on an authored cove: clipmap terrain from one
landscape field, habitat-placed vegetation and rocks, a lunar-gravity sea, surface
water after rain, the clear-sky atlas, and column-based clouds. It is the current
immersion development scene, not a claim about what an Open Moon shore would be
like; its landscape and ecology are placeholders until the geography and
biosphere research can supply them.

This folder is being reorganized into the new immersion layout (engine, world,
bake, experiences). The science developed here has already moved to its domains:

| Moved work | Now in |
|---|---|
| Clear-sky solver that produced `data/atmosphere.json` | [illumination/sky](../../../illumination/sky/) |
| Column thermodynamics (`weather-column.js`, `atmospheric-profile.js`) | [atmosphere/column](../../../atmosphere/column/) |
| A1–A3 light-transport references and probes | [illumination/references](../../../illumination/references/) |
| A1–A3 lab pages | [visualization/labs](../../../visualization/labs/) |
| B1 offline scene renderer | [visualization/reference-renderer](../../../visualization/reference-renderer/) |
| Checkpoint methods, handoffs and validation records | [immersion/docs/history](../../docs/history/) |

## Build and check

Place the four pinned Three.js r186 build files (`three.core.js`,
`three.module.js`, `three.tsl.js`, `three.webgpu.js`, hashes in
`vendor/three-r186/manifest.json`) in `vendor/three-r186/`, either with
`python tools/import_three.py /path/to/three.js-r186.zip` or from the `three@0.186.0`
npm package. Then:

```sh
python -m pip install -r requirements-surfaces.txt
python tools/generate_surfaces.py     # surface textures (checked against the manifest)
python build_landscape.py             # Open_Moon_Shoreline.html
node --test --test-concurrency=1 tests/*.test.cjs
```

`data/atmosphere.json` is the packed clear-sky atlas (`pack_atmosphere.py`). It is
committed because regenerating it needs the numba sky solve; regenerating it from
`illumination/sky` reproduces the committed data exactly.
