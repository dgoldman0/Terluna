# M1 methods and reference protocol

## Scope and provenance

This build revises the immersion renderer at Terluna checkpoint
`4877d8e6c020de51e237deb377955442778917da`. The seven original runtime files
were checked against their repository Git blob identities before editing.
The recovered source archive also contained earlier versions of some
packaging documents; new documentation here describes the actual revised build.
The previous method text is preserved in `docs/checkpoint-methods.md`.

The Moon, Earth and zero-ozone Moon atlas arrays are byte-for-byte unchanged.
They remain outputs of the previous optical proxy, including its prescribed
exponential density profiles, fixed-observer calculation and unresolved global
energy-accounting discrepancy of up to about 4.6%. No atmospheric transport,
climate, biosphere or orbital solver was rerun. M1 validates a rendered
clear-noon scenario and selected numerical/rendering properties.

## Reference camera and forcing

The reference uses Open Moon, zero-declination equatorial solar noon, phase 0,
weather kind `clear`, motion time 0 s and a frozen environment. The camera is
at local x=0, z=9 m, 1.7 m above the authored surface (absolute y approximately
3.045873320927338 m), yaw 0, pitch -0.04 rad, vertical field of view 58 degrees.
Reference captures use 1440 by 900 pixels and pixel ratio 1. The 960 by 600
integration capture uses the same aspect ratio and physical camera.

Clear forcing has zero cloud coverage, liquid water and precipitation, with
selected wind 1.5 m/s, air temperature 294 K and relative humidity 0.48. The
extra visibility parameter is 180 km; this is a selected scene parameter,
not an inferred lunar meteorological measurement. The inherited clear-noon
horizontal illuminance is approximately 94,158.49 lx for the Moon profile.

## Continuous land and seabed

`surfaceHeight(x,z)` defines a single authored elevation field relative to a
local sea datum. A smooth coastal transition, regional relief, local slopes
and nested noise scales extend continuously above and below that datum.
`waterDepth=max(0,-surfaceHeight)` is the corresponding vertical depth.

A single radial indexed mesh has 384 angular divisions and 512 nested radial
rings out to 16,384 m. It contains 196,609 vertices. Adjacent rings share
indices; the only open topological boundary is the outer circumference.
Geometric height includes spherical sag evaluated by the stable expression

    sag = (x²+z²) / (R + sqrt(R²-x²-z²)).

The finite mesh approximates the continuous field between vertices. Normals
use sampled field gradients and a small-angle curvature-gradient term. This
is authored terrain, with no claim to measured lunar coordinates or erosional
history. The original coast guide remains an input to the height field;
water intersection comes from the resulting geometry and depth buffer.

Rocks use eight seeded bevelled fracture-shape families. Substrate and slope
control material masks; deposition bands control pebble placement. Trees and
grass use hierarchical geometry and spatial instance batches. Materials use
linear-colour world-space detail and derivative filtering. All scenery remains
procedural; there are no scanned terrain or vegetation assets in this build.

## Water rendering

Water uses a shared-ring curved surface. Six selected long-wave modes displace
vertices; twelve modes contribute to shading normals. For each mode,

    omega² = g k tanh(k h),       h = 12 m,

with lunar gravity 1.62 m/s². The shader evaluates tanh through
`(1-exp(-2kh))/(1+exp(-2kh))`, avoiding an overflow observed in the software
renderer. Wave amplitudes and directions are prescribed. They provide a
finite-depth dispersion-based render model; fetch, spectral energy growth,
shoaling, breaking, tidal flow and shoreline run-up remain future systems.
The h=12 m dispersion depth is still fixed, even though optical depth now
uses the actual seabed.

The opaque pass writes half-float colour plus scene depth. Reconstructed
camera-space depth supplies the underwater optical path for refraction.
Foreground-depth rejection reduces refraction leaking through nearer objects.
The selected band absorption coefficients are [0.29, 0.071, 0.036] m⁻¹, with
Beer–Lambert transmittance exp(-coefficient × path length). These RGB proxies
are explicit artistic/optical inputs, not sampled spectral water measurements.
A selected in-scattering term supplies the remaining water-body contribution.

Water has F0=0.0204 Fresnel reflection. The near field uses a viewport-sized
half-float planar reflection target with four multisamples; the economy preset
halves reflection dimensions. Distant/off-screen reflection transitions to the
inherited atmosphere. Curvature and wave displacement make the planar mirror
an approximation. A finite solar-sized highlight is added explicitly while
the solar disk is disabled in the mirror capture. A small shallow-contact
ripple replaces the previous decorative white foam bands.

Each active frame refreshes the required opaque and reflection passes. Renderer
target, exposure, tone mapping, clipping, disk visibility and shadow state are
restored after those passes. Frozen idle frames issue no new rendering work. Clear-sky filtered environmental lighting is cached by optical state, so camera translation and restoring the same noon state reuse it. Cloudy states retain timed updates.

## Lighting and the display transform

The inherited radiance and irradiance retain the common scale factor 1/8500.
Environmental lighting is computed from a 256-pixel cubemap and filtered with
Three.js PMREM. Near shadows use 4096² maps at balanced/high quality and 1024²
in economy. Spatial batching culls inland geometry from the reference camera.

A selected local aerial-perspective proxy uses RGB molecular coefficients
[5.8e-6, 13.5e-6, 33.1e-6] m⁻¹, scaled by 1.2 for the lunar surface-density
proxy and 1 for Earth. Additional selected extinction depends on visibility.
In-scattering samples the inherited sky in the view direction. This local RGB
approximation supplements the atlas; it is not a new reciprocal spectral
terrain–atmosphere transport solution.

Daylight white balance uses Y(E)/E for each RGB channel of the inherited noon
direct-plus-diffuse irradiance. Moon gains are approximately
[0.8055076246, 1.0428335463, 1.4370070266]. This transform is applied once at
final tone-mapped output. It is disabled for linear intermediate captures and
can be disabled in the Camera colour control. It is a photographic display
choice, not a physiological eye-adaptation result or a modification of sky data.

An actual horizontal MeshStandardMaterial gray card with linear reflectance
0.18 is lit by the same noon direct light and environment. Float-target GPU
readback and the actual r180 ACES tone-mapping shader give a fixed exposure
of 0.236657345 for output linear luminance approximately 0.18. The camera
compensation control acts relative to that reference. The gray-card probe
validates that selected display convention, not absolute atmospheric accuracy.

## Browser validation protocol

The supplied dependency is the exact r180 CommonJS bundle, Git blob
`ca4833532c363b72477b2e8a6f47cc0e2fc7b09a`. The HTML embeds its original bytes
and the loader verifies the pin before evaluating code.

The preparation runner blocks URL navigation and provides software WebGL only
in headed Chromium on Xvfb. The recorded runs therefore inject the complete
HTML into about:blank using Playwright. This origin lacks native Web Crypto,
so the harness bridges the single SHA-1 operation to Python hashlib; the
application's own expected-hash comparison remains active. This bridge exists
only in `tests/capture_m1.py`. No browser policy is changed. The normal harness
path serves localhost and uses native browser cryptography; that path and
native file opening were not executed in this runner.

requestAnimationFrame callbacks are queued by the harness so each capture has
a defined simulation time. The real application rendering and shader code are
used. Integration checks explicitly advance those callbacks to inspect paused
redraw behaviour and exercise actual view/control operations. GPU depth samples
and float gray-card readbacks supplement screenshots. Software-run timings are
recorded as test duration, not consumer-device frame-rate estimates.

## Remaining quality work

The reference has coherent land/water contact and improved exposure/detail.
Foliage remains visibly procedural, distant terrain has limited material
structure, and shoreline wetting/contact is simplified. The lifelike M1 image
criterion remains open. M2–M4 require separate motion, weather/hydrology and
full-cycle validation. Earth mode inherits a Moon-created local mesh and has
not been certified as a geometrically complete Earth comparison.
