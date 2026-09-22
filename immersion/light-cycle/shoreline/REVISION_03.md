# Shoreline revision 03 — r186, connected ponding and the renderer lab

Base repository commit: `caa4caa2b2855e9ed8f03597317a05d8224949e6`.
This revision is prepared locally. Its application utility and patch operate on
ordinary source files; neither creates a commit nor moves a remote branch.

## The controlled dependency migration

The supplied `three.js-r186.zip` is the dependency input. Its version, archive
hash and selected module hashes are recorded in `vendor/three-r186/manifest.json`.
Four original modules, their MIT license and the shared core are sufficient for
the two viewers. Examples, manual assets and fonts are excluded.

The previous browser loader executed the r180 CommonJS bundle. The r186 build
uses ES modules, linked through local Blob URLs after verifying each original
module. The WebGL viewer embeds core and WebGL modules; the lab additionally
embeds WebGPU and TSL. Both outputs carry their atmosphere, material images and
required dependency bytes. No dependency download occurs during building or
opening. The loader verifies the original bytes before rewriting their declared
import edges. Altered-byte rejection is covered by a test.

The unchanged-scene r186 comparison retained the original camera, Sun, weather,
terrain, plants, material data and exposure. PCFShadowMap replaces the deprecated
soft-PCF selection. At the original exposure of 0.236657345, the GPU gray-card
probe returned 0.17222686, compared with approximately 0.18 under r180. The original
failed 18% assertion is retained in `validation/r186/migration-optics.json`.

Only after recording that comparison was the new default exposure explicitly
recalibrated to **0.244647046**, a 3.376% increase. The final gray-card check returns
approximately **0.17999983**. The main viewer and renderer lab use that same new
exposure. Atmospheric input bytes and input radiances were retained. This is a
recorded camera calibration change, separate from the landscape-model changes.

## Model state now reaches additional visible surfaces

### Connected pond reconstruction

`ponds.js` reconstructs visible water from the existing conservative surface-water
columns. A one-metre triangulated sample of the authoritative terrain supplies
spill basins inside the local four-metre water-grid footprint. Connected samples
at the same spill elevation form a basin. Column pond storage is apportioned to
intersecting basins by their sub-grid storage capacities.

For a basin, the model integrates the depth above each terrain triangle and solves
for one water level whose integrated volume equals its assigned storage, capped
at the basin's spill capacity. Partially flooded triangles are clipped at that
level. This yields an actual free-surface mesh and a terrain-following shoreline,
with the same curvature datum as the ground. Each basin has a single level;
independent columns no longer become independent flat water patches.

The allocation is an explicit sub-grid reconstruction assumption. It leaves the
column water ledger unchanged. Water assigned outside resolved basins remains
reported as mobile water; water exceeding a basin's capacity remains reported as
above-spill storage. Neither term is silently rendered as a larger pond or
removed from the hydrologic budget. The operation is not a shallow-water dynamics
solver. Its fixed terrain assumes unchanged landforms during the weather replay.

In the four-hour after-rain scenario, the current model has 4.565945 m³ of column
pond storage. Approximately 3.408699 m³ is reconstructed in two wet basins,
0.742043 m³ remains mobile outside the basins and 0.415204 m³ is above-spill storage.
The reconstruction budget residual is approximately 2.6e-11 m³. These are internal
model results from selected scenario parameters, rather than measurements of a
real landscape. The triangulated water meshes account for the reconstructed
volume within the tested floating-point tolerances.

The material system reads the same reconstructed head texture. Wet-bank shading
uses the difference between ground elevation and pond level, alongside distinct
film and soil stores. This avoids treating an entire coarse cell as a visible
pond. The viewer's **Inspect ponds** control moves to the largest wet basin after
a rain replay; **M1 noon** restores the original reference camera and dry surface
state with the new, documented r186 camera calibration.

### Shared bathymetry controls the wave envelope

The sea now reads the same elevation raster used by the landscape fields. The
prescribed wave amplitude envelope tends to zero as water depth tends to zero.
For sampled depth h, deep-mode amplitude sum A and limit q = 0.45, the envelope is

    e = h / (h + A/q),     A e <= q h.

The CPU diagnostic, WebGL shader and node renderer use this same construction.
The normal calculation includes the spatial derivative of the envelope, so a
changing amplitude does not leave the old deep-water normal behind. Tests cover
land suppression, boundedness, the depth constraint and finite-difference normals.

This is a depth-limited prescribed-wave rendering model. Frequency/phase modes
remain the existing fixed-depth scenario. Refraction, energy-conserving shoaling,
breaking, run-up, variable-depth dispersion and sediment transport remain open.
The local bathymetry sample is two metres; the raster-to-fine-terrain discrepancy
and the fallback beyond its footprint remain explicit resolution boundaries.

### Moving terrain retains its existing topology

The terrain update now reuses a parent-sample scratch object, avoids redundant
height evaluation and normal normalization, accumulates bounds in the existing
vertex pass and marks only changed attributes for upload. The geometry still
uses the same eleven nested levels, stitched edge topology and camera-following
resolution. Existing continuity, collision and moving-grid tests remain active.
This is CPU work reduction; no GPU terrain generation or compute culling is claimed.
Across three alternating process pairs, the median of per-run median update times
fell from 22.21 ms to 16.40 ms (approximately 26%). Individual tail times remain
material; this CPU-only result makes no GPU frame-rate claim.

The regional landforms, tree architectures, planting identities, canopy suitability
and soil parameterizations remain those of the recovered checkpoint. More natural
vegetation, detailed material-height displacement and regional vegetation LOD are
still substantive development work.

## The optional TSL / WebGPU renderer

`Open_Moon_Renderer_Lab.html` builds the production landscape, ecology and pond
geometry through the shared modules. Renderer-specific work lives in
`node-materials.js` and `renderer-lab.js`. It ports coordinated surface textures,
field-driven materials, clear-sky LUT lookup, local aerial perspective, direct
shadows, depth-based water absorption, planar reflection and pond shading into
Three.js node materials. It does not create a simplified replacement landscape.

The lab selects a real WebGPU adapter when the browser exposes one, or labels its
TSL/WebGL2 fallback explicitly. An explicit native-WebGPU request fails visibly
when it cannot initialize. The present native path requires float32-filterable
field textures; unsupported adapters receive an explicit explanation. Adapter
information, backend, software-renderer detection, resolution, scene state,
feature coverage and timing method are included in exported measurements.

The runner permitted the in-memory browser tests but blocked local-file navigation;
its in-memory pages did not expose WebGPU. Therefore, the node implementation was
executed through **WebGPURenderer with a WebGL2 backend on SwiftShader**. This proves
that the node implementation renders and consumes the shared fields through that
backend. Native WebGPU shader execution, device stability, consumer-GPU speed and
visual headroom remain unmeasured.

Clear-noon image comparison is appropriate for the implemented overlap. Exact
foliage transmission, leaf sway, cloud transport and animated rain remain parity
gaps. Wetness replay in the lab changes water stores while retaining clear-sky
lighting. Compare its hydrologic state with the reference independently of the
reference viewer's weather-dependent cloud appearance.

The lab includes a repeatable route and measurement export. Frame measurements
include CPU submission and all scene passes. WebGL2 completion is enforced by a
synchronous one-pixel readback followed by finish; native WebGPU uses the device
queue's completion promise. These diagnostic latencies include synchronization
cost and are not unsynchronized display frame rates. The preliminary finish-only
records are preserved as diagnostics and are excluded from the final comparison.
The prototype's incomplete parity and software adapter prevent a promotion claim.
In the synchronized 12-frame software sample, the WebGL viewer had a median of
3.10 seconds and the TSL/WebGL2 prototype 3.75 seconds. Shadow update scheduling
and exact leaf lighting still differ. These end-to-end diagnostic samples neither
measure native WebGPU nor settle a hardware renderer choice.

## Evidence and remaining gates

`REVISION_VALIDATION.json` indexes the current unit results, browser captures,
actual GPU field readback, gray-card/depth checks, reset, mobile controls, route
checks and diagnostic failures. Current short walking checks succeeded; they do
not isolate the earlier context-loss cause or establish long-duration stability.
The initial shader-port mistakes and the readback-harness RAF scheduling failure
are preserved separately from accepted final results.

The next visual priorities are material-height geometry, coherent exposed-rock
and sediment structure, ground-cover/crown realism and vegetation LOD. The next
renderer gate is a consumer-machine native-WebGPU run against the same reference
state, including missing-feature parity and longer traversal. The node path stays
optional while those gates are open. Browser/library changes alone do not settle
the scene's visual acceptance.

## Implementation references

The pinned release source is the API/shader authority for this revision. The
Three.js migration guide and renderer documentation informed the dependency
transition:

- https://github.com/mrdoob/three.js/wiki/Migration-Guide
- https://threejs.org/docs/pages/WebGPURenderer.html

The inherited landscape methods and their design references remain in
`LANDSCAPE_METHODS.md`. Third-party licensing is retained in
`THIRD_PARTY_NOTICES.txt` and embedded in both viewers.
