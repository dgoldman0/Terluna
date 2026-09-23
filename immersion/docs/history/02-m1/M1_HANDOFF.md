# M1 progress handoff

## Working state

This source checkpoint advances the shoreline prototype from base commit `4877d8e6c020de51e237deb377955442778917da`. The editable source, build scripts, model inputs, tests, and compact validation records are preserved directly in the repository. The previous generated HTML remains recoverable from the base commit; it is removed from the current tree to prevent that older render from being confused with the revised M1 source.

The validated local review package used the exact Three.js r180 bundle with Git blob `ca4833532c363b72477b2e8a6f47cc0e2fc7b09a` and embedded it into a self-contained HTML build. The repository build keeps the dependency optional: placing that verified bundle at `vendor/three.cjs` embeds it, while its absence leaves the pinned first-open dependency path intact.

## Current result

The noon scene now has continuous land/seabed, explicit optical water depth, viewport-scale reflections, reworked foreground geometry/materials, calibrated fixed exposure, and a separately reversible daylight white balance. The final 1440×900 view/reset sequence in the review package passed and restored a pixel-identical reference PNG. A current-build GPU probe measured gray-card output and submerged scene depth. The 390×844 layout was also rendered and checked.

A prior full-resolution debug-readback stress run lost its software WebGL context. That diagnostic remains recorded. Clear-sky lighting is cached across camera/reset operations; final full-resolution ordinary view/reset and separate GPU readback tests passed. Consumer-GPU behaviour remains to be measured.

## Keep M1 open

The scene remains visibly procedural. The next pass should improve foliage silhouettes and thin-leaf response, foreground material structure, regional landform/material detail, and shoreline contact while holding the same clear-noon reference camera and display calibration fixed. Asset or procedural upgrades should preserve provenance and the existing physical assumptions.

The water still uses selected wave modes and a fixed 12 m dispersion depth; its optical path now comes from the actual rendered seabed. Run-up, breaking, wind-wave growth, spatial hydrology, and predicted ecological placement remain later work.
