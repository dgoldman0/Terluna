# Rendering the column clouds

*From shoreline revision 04. The column thermodynamics these passes draw are
documented with the model in [atmosphere/column/METHODS.md](../../atmosphere/column/METHODS.md);
the reference light transport that measures this renderer is in
[illumination/references](../../illumination/references/).*

## The shared fields

`src/cloud-renderer.js` uploads a 256-sample vertical texture containing
extinction, ice fraction and two horizontal wind components. Its altitude bounds
come from the diagnosed layer. The fog is locally stratified; fair clouds use a
continuous multiscale density field; convection uses unequal overlapping lobes
and spreading upper heads within the diagnosed buoyant envelope. Horizontal
scales, occupancy and the shape of the lobes are selected morphology inputs.
They do not emerge from a cloud-resolving fluid simulation.

Wind advects each altitude independently, producing shear from the sounding.
A frozen study retains one deterministic state. Resuming motion advects that
state while holding its thermodynamics fixed. Cloud formation, dissipation,
vertical velocities and precipitation evolution remain further modeling work.

`src/cloud-optics.js` integrates spherical paths through the inherited molecular
profile to obtain solar transmittance at each cloud height and local solar angle.
A 128-by-64 RGB lookup stores those results. The ray/planet intersection is also
checked directly in the shader, preventing interpolation from lighting a cloud
through the solid Moon. This permits elevated cloud to remain sunlit after the
local surface enters night. The source irradiance is anchored to the inherited
noon direct irradiance; no camera-exposure adjustment accompanies this revision.

The optical proxy uses RGB molecular extinction and an explicitly selected
10–40-km triangular ozone layer, with effective RGB optical depths
[0.015, 0.035, 0.003]. The no-ozone comparison disables this term. The solar path
includes these molecular and ozone terms; the eye-to-cloud path currently uses
molecular extinction. Diffuse cloud illumination combines downwelling sky with upwelling ground
reflection using a selected regional Lambertian albedo of 0.18, vertical molecular
transmittance and a height-dependent boundary weighting. This replaces a sky-only
closure that gave excessively blue undersides. Regional albedo, internal multiple
scattering and foreground molecular in-scattering remain approximations. The original
spectral-sky calculation and its previously documented energy-accounting and
twilight uncertainties remain unchanged. Bruneton's rendering work [4] is a
reference for a fuller future transport solution, not a claim that the present
cloud closure reproduces that solver.

## Bounded render passes

An initial direct integration path lost its SwiftShader context. That failed
record is retained. The accepted implementation integrates cloud radiance and
transmittance into a cached angular sky and renders a separate cloud-shadow map.
The main view, sea reflection and environment lighting reuse these fields.

A deterministic 64³ R8 noise lattice supplies cubic-interpolated, periodic
world-space noise. This removes repeated hash/trigonometric evaluation from the
view and shadow steps. Its bytes and actual GPU samples are checked. The complex
raymarching shaders are confined to the bake passes; surface materials sample the
cloud-shadow cache.

The sky cache has 768×384, 1536×768 or 3072×1536 pixels according to the existing
Economy, Balanced or High setting. Its elevation mapping allocates extra samples
near the horizon. Integration uses 128 primary samples and eight secondary
cloud-shadow samples. Expensive baking is tiled into at most 256×128-pixel draws;
this bounds individual work submissions while preserving global sample jitter.
The solar disk is rendered separately at full output resolution.

Wind updates are quantized to two seconds. Observer translation invalidates the
cache at 32-m horizontal cells and one-metre rounded eye-height changes. Solar
angle, atmospheric regime, world and render quality invalidate it immediately.
The ground shadow map spans 32.768 km on a curved sea-level surface at 128×128
resolution. An approximate coverage-weighted shadow value supplies the far
boundary transition. Building-height differences in cloud shadow, sharp small
cloud-shadow features and finite cache angular resolution remain approximations.
These budgets are exported with the saved scene state. Submission duration is
explicitly labeled CPU submission time, not completed GPU time.

