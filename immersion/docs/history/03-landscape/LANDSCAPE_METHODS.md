# Landscape, materials and ecology checkpoint

## Evidence and scope

This checkpoint extends the eight verified M1 runtime files from commit
`4c7affae780d3ff7ef5581e88aac7ae32ada7505`. The inherited atmosphere dataset,
atmospheric runtime, audio runtime, pinned Three.js r180 bytes, noon exposure,
white balance, and reference camera are retained. The new landscape is an
explicit authored scenario for an introduced biosphere. Numerical conservation,
rendered appearance, and empirical physical validation are separate evidence
states. The current build provides the first two to the extent recorded in
`LANDSCAPE_VALIDATION.json`; ecological and hydraulic parameters remain illustrative.

## Terrain authority and drainage

`landscape.js` owns a deterministic two-dimensional elevation field relative to
mean water level. Authored regional ridges, valley polylines and coastal shelf
forms establish composition. Warped ridged structure and substrate-dependent
small relief add spatial scales. No independent sinusoidal shoreline guides the
renderer: the coastal query finds a bracketed zero of the same elevation field.
That query selects one local crossing; it is not a global nearest-coast solver.

A 257 by 257 raster at 2 m spacing covers x = [-256, 256] m and
z = [-128, 384] m. Priority flooding establishes spill paths for depressions;
receiver routing and D8 contributing-area accumulation provide connected
catchment information. An initial routing pass guides a modest authored incision,
and final drainage is recalculated from the resulting surface. This produces a
coherent drainage scenario, without solving geological evolution or calibrating
lunar erosion rates. Raster values describe node-centred 2 m cells, so the
summed contributing area includes the outer half-cell margin.

Substrate affinity, sediment/rooting depth, slope, curvature, exposure, moisture
suitability, material weights and plant-community suitability are shared fields.
Some are geometry-derived; others are declared scenario parameterizations. Two
managed shore plantings have corresponding improved rooting zones. A prescribed
moisture-suitability field supplies the initial effective soil storage. It does
not change into an empirically validated climate prediction through this coupling.

The benchmark elevation is anchored to the M1 surface value at (0,9) m, preserving
the original camera height. Body curvature is applied separately, consistently
for Moon/Earth comparison, terrain, placed objects and collision.

## Camera-following terrain

`terrain.js` constructs eleven nested 128-cell grids, from 0.25 m spacing near
the observer to a 16,384 m outer half-extent. Grid origins snap to aligned parent
spacing. Outer bands morph towards the coarser triangle surface; explicit edge
collapse and a split corner triangle resolve inter-level topology. Index
patterns and sampled elevations are cached with bounded storage. A child-grid
origin change can update its parent's hole topology without recalculating all
parent vertices.

All passes consume the same displaced geometry: colour, depth, reflection and
shadow. Terrain collision samples the authoritative elevation field; the
reported near-camera triangle-versus-height discrepancies are finite-resolution
errors, measured at selected camera poses rather than globally bounded.
Centimetric procedural relief exists in the terrain geometry. The packed
material height channel generates coordinated normals and remains available for
future surface displacement; it does not currently displace terrain vertices.

Tests weld active boundary vertices and count triangle-edge incidence across all
levels, including three moved poses. The only unpaired edges should be the outer
terrain perimeter. These tests catch cracks and the previously discovered
upper-right transition T-junction. They do not establish perceptually invisible
LOD morphing on every walking route. CPU grid rebuild cost is recorded separately
from GPU rendering and is still an optimization target.

## Material channels and scale

Six original deterministic materials cover sand, gravel, basalt, soil, litter
and bark. `tools/generate_surfaces.py` produces coordinated albedo/AO and
normal/roughness/height images using authored grain, pebble, fissure, aggregate,
leaf-litter and bark-ridge morphology. `assets/surfaces/manifest.json` records
seeds, physical tile scales, channel encodings, provenance and file hashes.
No external scenery image assets are included; licensing remains the author's
choice under the project's existing license practice.

The PNGs pack layers for transport. At load time, each of their six tiles becomes
a separate repeating, mipmapped 512-pixel texture-array layer, avoiding colour
bleeding between materials at coarse mip levels. Albedo uses sRGB decoding;
normals, roughness, height and AO remain linear data. Ground blending and object
placement sample the same material field. Rocks use triplanar projection and
branches use length/circumference-scaled bark coordinates. Foliage has an explicit
lamina, midrib/vein detail and an approximate shadowed thin-sheet light term.
Wetness modifies albedo, roughness and dielectric specular response. This is an
approximate layered response rather than a measured optical transport solution.

## Vegetation and objects

`ecology.js` selects stable world-cell identities using suitability, rooting
depth, slope, exposure, moisture, canopy, spacing and explicit path/shelter
exclusions. Mature plants retain their identities across observer movement,
quality changes, weather replay and instantaneous Sun direction. Three branching
families and two age classes vary crown structure and leaf form. Main limbs
start at explicit trunk nodes, child limbs connect to parent endpoints, and
leaves attach to terminal twig segments through geometric petioles. The
architecture is procedural, without species-specific growth validation.

Canopy coverage is a planted-crown proxy used for litter, interception and
understory placement. Actual geometry supplies direct shadows. Habitat light
availability is a prescribed exposure/canopy proxy; monthly radiative ecology,
plant growth, succession and viability predictions remain future work. Grass,
moist-margin tufts, understory and substrate-associated rocks share the fields.
Instance batches support spatial culling. Large-scale vegetation streaming and
hierarchical distant vegetation representations remain open.

## Spatial water bookkeeping

`surface-water.js` maintains a 65 by 65 grid at 4 m spacing over
x = [-128,128] m and z = [-64,192] m. Each cell has leaf interception,
surface film, effective soil storage and ponded water. Initial soil storage is
65% of the capacity-scaled static moisture-suitability value on land; film, leaf
and ponded stores start dry. The effective capacity and permeability are selected
material-dependent proxies, not measured volumetric hydraulic properties.

Rain fills interception and surface stores, then infiltration and drainage act.
D8 routing exports water at sea/domain boundaries. Downstream-first processing
moves water at most one cell per substep. The shelter bypasses incident rain as
an explicit roof export. Evaporation, deep drainage, roof loss and boundary loss
are individually accounted. The ledger is

    initial storage + rainfall = current storage + evaporation
                               + deep drainage + boundary export + roof export

All terms are accumulated in mm-cell units and converted to cubic metres by the
16 square metre cell area divided by 1000. Substeps are at most 10 seconds.
Backward seeking reconstructs the same antecedent state and replays forcing;
drying continues beyond the four-hour prescribed weather episode.

Texture channels expose film fraction, effective soil-storage fraction, ponded
millimetres and leaf wetness to the renderer. A GPU readback test compares their
actual uploaded values with CPU state. Ground-level pond storage modifies surface
appearance, while connected free-surface puddle geometry, shallow-water dynamics,
beach run-up, breaking and variable-depth wave dispersion remain unimplemented.
The sea retains the M1 optical depth path and selected fixed-depth wave modes.

## Rendering and review boundaries

Clear-noon exposure and camera are held fixed for comparison. Software WebGL is
detected explicitly; SwiftShader/llvmpipe start in Economy and disable
multisampling after higher-setting software paths failed during review. Hardware
keeps Balanced at desktop width and the multisample path. This changes anti-aliasing, not the inherited exposure or atmospheric data.
The exact backend is recorded in saved state. Hardware-GPU frame rate, high-DPI
quality and mobile-device thermal performance have not been measured.

The browser harness queues actual application animation callbacks for repeatable
stills. Restricted environments can inject the exact HTML into about:blank and
bridge its SHA-1 check to Python; that bridge exists solely in the test harness.
Explicit GPU completion between walking test callbacks bounds the queued
software workload. Sustained SwiftShader traversal still fails, including an
Economy run at 480 by 300. These failures remain in the validation record
alongside successful reference, woodland, weather, reset and mobile checks.
A Mesa/ANGLE browser attempt did not expose WebGL2 in this runner.

## Design references

The nested-grid design follows the general geometry-clipmap family described by
Asirvatham and Hoppe in *GPU Gems 2*, chapter 2, "Terrain Rendering Using
GPU-Based Geometry Clipmaps" (2005). This implementation uses CPU geometry updates
and its own transition topology rather than their GPU elevation-update pipeline.
The shipped Three.js r180 source is the authoritative shader/API reference.
Inherited atmospheric and wave references remain in `SOURCES.md` and `METHODS.md`.
