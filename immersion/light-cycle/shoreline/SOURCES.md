# Source and dependency record

## Existing Terluna material

The optical data come from the already-produced `Open_Moon_Full_Cycle_Source.zip`
(SHA-256 in `METHODS.md` and `scenario.json`). Its corresponding source/method
workspace was committed at:

- https://github.com/dgoldman0/Terluna/tree/71831e89f154cbf1d7f0b2eb3f3fe899847f5d1b/immersion/light-cycle
- https://github.com/dgoldman0/Terluna/blob/71831e89f154cbf1d7f0b2eb3f3fe899847f5d1b/immersion/light-cycle/METHODS.md

This release reads and repacks the archived NPZ outputs. It does not claim new
validation from those preceding source records. Original atmosphere/array hashes
are in `data/atmosphere.json`. No theoretical or engineering model was replaced.

## Three.js

Pinned release: r180, package version 0.180.0. This is a deliberate reproducible
version pin, not a claim to use the latest available release.

- Package metadata: https://github.com/mrdoob/three.js/blob/r180/package.json
- CommonJS bundle: https://github.com/mrdoob/three.js/blob/r180/build/three.cjs
- Distribution URL: https://cdn.jsdelivr.net/npm/three@0.180.0/build/three.cjs
- Expected Git blob: `ca4833532c363b72477b2e8a6f47cc0e2fc7b09a`
- License: https://github.com/mrdoob/three.js/blob/r180/LICENSE
- Material extension hook inspected in the pinned source:
  https://github.com/mrdoob/three.js/blob/r180/src/renderers/shaders/ShaderChunk/lights_fragment_begin.glsl.js
- WebGL capabilities inspected in the pinned source:
  https://github.com/mrdoob/three.js/blob/r180/src/renderers/webgl/WebGLCapabilities.js

The loader computes the Git blob hash over `"blob " + byte_count + NUL + bytes`
and rejects an altered file before evaluating code. Metadata, license and the
small interface sources were accessible through the connected GitHub reader.
The bundle bytes could not be retrieved into the preparation container. A failed
attempt to materialize its connector reference produced no file. No placeholder
library is shipped. The first-open runtime request is explicit in the UI.

## Inherited optical method references

- Bruneton, E. (2017), *Precomputed Atmospheric Scattering: a New Implementation*.
  https://ebruneton.github.io/precomputed_atmospheric_scattering/
  The preceding optical implementation used the reference solar and ozone arrays;
  its third-party notice is preserved in `LICENSES.txt` and in the HTML.
- Wyman, C., Sloan, P.-P., and Shirley, P. (2013), *Simple Analytic Approximations
  to the CIE XYZ Color Matching Functions*. Journal of Computer Graphics
  Techniques 2(2). https://jcgt.org/published/0002/02/01/
  This is inherited provenance for conversion of the archived spectra to RGB.

These links identify methods and input provenance; this build does not claim a
new full-paper review or independent reproduction of their numerical results.

## Browser audio

The spatial listener and sources use the Web Audio API PannerNode interface:
https://developer.mozilla.org/en-US/docs/Web/API/PannerNode

Sound textures, local scene geometry and material texture patterns are procedural
project code. No photographs, scanned commercial assets, recorded sound libraries,
font files, or external papers are redistributed.

## Model equations

The water reservoir and solar/curvature/wave relations are specified explicitly
in `METHODS.md` and `src/core.js`. The selected spherical drag correlation is
implemented as `Cd = 24/Re × (1+0.15 Re^0.687)` below Re=1000 and 0.44 above.
Using this selected correlation with a single rigid spherical drop is an input
assumption; no lunar-drop measurement or complete precipitation model is claimed.
