# Open Moon — Shoreline M1 progress checkpoint

This directory contains the editable source for the current M1 clear-noon shoreline pass, based on Terluna checkpoint `4877d8e6c020de51e237deb377955442778917da`.

## Status

This is a **progress checkpoint**, not M1 visual acceptance. It is substantially improved over the preceding shoreline prototype: one indexed terrain field now continues from foreground land through the seabed and regional relief; water reads rendered seabed depth for optical attenuation; reflections scale with the viewport; foreground rocks, pebbles, sand detail, vegetation geometry, exposure calibration, and daylight white balance have all been reworked.

The scene still reads as visibly procedural. Foliage silhouettes and thin-leaf response, foreground material structure, regional landform/material detail, and shoreline contact remain the main visual weaknesses. M1's convincing/lifelike appearance criterion therefore remains open.

## Build

`src/` contains ordinary editable JavaScript and `index.template.html` is the retained checkpoint page shell; `build.py` applies the M1 shell changes explicitly while assembling the current viewer. `data/atmosphere.json` retains the inherited packed sky data unchanged. Run:

```sh
python build.py
```

If `vendor/three.cjs` is present and matches the pinned Three.js r180 Git blob `ca4833532c363b72477b2e8a6f47cc0e2fc7b09a`, the build embeds it. Otherwise the generated HTML requests that pinned library on first opening and retains the existing offline-export path. The generated `Open_Moon_Shoreline.html` is intentionally not retained in this progress commit so an older packaged build cannot be mistaken for the current source state.

## Test

For the complete 54-test M1 validation path, first place the verified r180 bundle at `vendor/three.cjs`, run `python build.py`, then run:

```sh
node --test tests/core.test.cjs tests/m1.test.cjs
```

The M1 package used for this checkpoint passed 54/54 unit tests and was rendered in Chromium/SwiftShader at 1440×900 and 390×844. `VALIDATION.json` records the executed checks and their evidence boundary. Browser screenshots and the self-contained review HTML were retained as review artifacts outside this source checkpoint; consumer hardware-GPU performance remains untested.

See `M1_HANDOFF.md` for the next visual work and `METHODS.md` for physical/model scope.
