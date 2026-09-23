# Landscape/ecology progress handoff

Base shoreline commit: `4c7affae780d3ff7ef5581e88aac7ae32ada7505`.
Recovery preserves the exact reviewed runtime and page-template bytes.
The original review record remains unchanged; see validation/recovery.json.

## Preserved landscape work

The authoritative two-dimensional terrain feeds rendering, collision, material
classification, water depth and ecological placement. Its local 2 m raster owns
final routed drainage, substrate/rooting depth and habitat fields. Eleven moving
nested grids provide 0.25 m near spacing and 16.4 km outer half-extent. Active
triangle-edge tests verify stitched transitions over four tested observer poses.

Six physically scaled material sets include albedo, normal, roughness, height and
micro-occlusion data. Each layer has independent repeating mipmaps. The generator,
manifest and original source are included. The transported height channel is
available for future geometric/parallax refinement; the terrain already includes
separate small geometric relief.

Habitat-selected trees, grasses, understory, rocks and gravel retain stable
identities. Tree limbs connect at explicit structural nodes and leaf petioles
start on terminal twigs. Crown proxies supply canopy interception and litter;
actual mesh geometry supplies direct shadows. Two managed shore trees and their
rooting zones are explicit scene-authoring choices.

A local 4 m grid tracks interception, surface film, soil and ponded stores, with
antecedent soil moisture tied to the static habitat field. Input, evaporation,
deep drainage, roof loss and boundary exports close a conservative ledger.
The renderer samples the uploaded state. Existing scene controls, atmospheric
inputs, noon exposure, initial camera and source-first build remain available.

## Review targets

Use the exact opening clear-noon reference, Water's edge, Woodland path and
Overlook. Turn inland on the woodland path to inspect density and connected
branching. Conditions → Landscape fields exposes elevation, material weights,
drainage, spatial water and canopy/moisture. Seek rain and clearing while watching
the local film, ponded water, soil-storage fraction and leaf wetness readout.
“M1 noon” should restore both model state and the reference image.

## Known limitations and first follow-on work

The scene still reads as procedurally constructed. The improvements establish a
shared landscape and material/ecology system; visual acceptance remains open.
High-value next work is foreground ground cover and foliage naturalism, more
convincing exposed-rock and sediment structure, material-height displacement,
shoreline contact/run-up, and connected visible pond surfaces. Habitat exposure
and moisture are scenario proxies, with no monthly ecological viability solve.

Performance needs a consumer-GPU profile before increasing geometry. Current CPU
near-grid rebuilds can consume a substantial part of a frame. Vegetation uses
spatial instance batches but still needs hierarchical LOD for distant crowns and
regional coverage. Software WebGL starts in Economy, and sustained SwiftShader movement still
exposes a graphics-context failure. Short, separated review cases pass; hardware
retains its separate default and remains unmeasured.
Consult the exact-build validation for successes and failures rather than treating
an earlier screenshot or original 54-test M1 record as current evidence.

## Repository recovery and next revision

Source comes from the saved `Open_Moon_Shoreline_Landscape_Source.zip`, verified
against its manifest and the previously delivered HTML. The previous upload
attempt left unreferenced Git objects and never updated main. Recovery uses only
objects matching the saved source; incomplete or differing objects are excluded.

The M1 repository builder is preserved with its landscape delegation. Surface
maps are reproduced from the ordinary generator and reference manifest; README
records the build steps. Original browser results are historical evidence,
including the failed walking run. Recovery reruns the build and numerical tests;
it does not imply a new browser or consumer-GPU test pass.

Follow `RENDERER_ROADMAP.md` for the agreed r186 migration, improved shared-field
integration and parallel WebGPU/TSL experiment. Keep those changes separate from
this source recovery. Compare equal scene state and record the actual backend;
retain WebGL until the alternative demonstrates its value.
