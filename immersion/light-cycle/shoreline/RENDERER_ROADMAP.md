# Renderer and model roadmap

Decision recorded 2026-09-22. The recovered landscape checkpoint remains on
Three.js r180 / WebGL. The supplied r186 archive is the dependency candidate for
the next revision. Keep the recovery and dependency migration in separate commits.

## 1. Controlled r180 to r186 WebGL migration

Hold the reference camera, Sun, atmosphere bytes, exposure, white balance, terrain,
plant identities and surface-water state fixed. Inspect the supplied package and
migration guidance; pin the exact release files and hashes. Adapt the loader,
build process, shader hooks and tests as required, including module dependencies.
Retain all third-party license notices and exclude unrelated examples/manual
assets from the project dependency bundle.

Compare shader compilation, gray-card output, actual seabed depth/refraction,
reflections, shadows, wetness upload, reset, controls and continuous traversal.
Record image differences and regressions before changing the world model.
A version upgrade alone is not a visual-quality acceptance result.

## 2. Better model-to-render integration

Prioritize connected drainage morphology and sediment/substrate structure,
spatial moisture and ponding, visible connected pond surfaces, meaningful near
surface displacement, vegetation-community/canopy/light relationships, shoreline
contact and regional terrain/vegetation LOD. Keep the same world fields
responsible for geometry, collision, material identity, placement and water.
Preserve explicit distinctions between authored scenarios, numerical checks and
empirical ecological or hydraulic validation.

## 3. Separate r186 WebGPU/TSL evaluation

Keep WebGL available as the reference path. Build an experimental node/TSL
renderer from the same landscape/ecology/surface-water state, with explicit
feature-parity checks for atmosphere, materials, shadows, water and camera
calibration. Isolate renderer-specific code; avoid adding unnecessary new
onBeforeCompile coupling. Porting existing GLSL hooks is real implementation work.

Record the actual backend and adapter, including software adapters and WebGL
fallbacks. A renderer class name alone does not establish hardware WebGPU use.
Compare identical scene state, camera, resolution, quality and visible content.
Measure startup, shader compilation, median and tail frame times, memory,
terrain-update hitches, culling/LOD cost and stability along a repeatable route.
Measure consumer hardware separately from software rendering. Only make compute,
placement or culling claims for paths actually implemented and exercised.

Promotion requires preserved model/optical behavior and demonstrated useful
performance, stability or visual headroom. Keep the experiment optional until
that evidence exists. Evaluate Babylon.js or a purpose-built terrain renderer
only when measured needs justify the additional migration and maintenance cost.

## Exit evidence

Produce separate dependency-regression, model-improvement and backend-comparison
records. A useful checkpoint should show both whether the landscape is visibly
more coherent and whether WebGPU renders that same richer model more effectively.
Keep unresolved failures visible, including the present SwiftShader context loss.
