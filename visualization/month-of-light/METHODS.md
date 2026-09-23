# Full-cycle methods and interpretation

## Clear-sky input

The sky radiance, irradiance and cloud-altitude light come from the illumination
domain's clear-sky solver. Its method, optical profiles and numerical checks are
documented in [illumination/sky/METHODS.md](../../illumination/sky/METHODS.md); this
viewer displays those atlases and adds the display models below.

## Full-cycle geometry and clocks

Let p range from 0 to 1 and H = 2πp. The idealized equatorial, zero-declination
solar direction in the local scene is

    s = (sin H, cos H, 0).

Local +x is west, +y is up and the camera can rotate through 360°.
The Sun's elevation is asin(cos H). The sequence is noon (+90°), sunset (0°),
midnight (−90°), sunrise (0°), and noon (+90°). The direction changes continuously;
the horizontal direction reverses at the zenith/nadir where its length vanishes.
The clear-sky field is sampled at the matching elevation and rotated to the
appropriate west/east half of the trajectory. The stationary optical atmosphere
is identical on the morning and evening branches.

Elapsed time is p × 29.53059 × 86400 seconds for the Moon or p × 86400 for Earth.
The clocks change the timing, not the optical data. Longitude, latitude, season,
orbital eccentricity, libration and actual ephemerides are outside this control.
Atmospheric refraction and local apparent solar contact times remain uncomputed.

## Prescribed weather display

The added conditions are reproducible optical/graphic scenarios. No convection,
moisture budget, condensation, precipitation microphysics, atmospheric momentum,
weather probability, or ecosystem response has been solved.

The WebGL cloud layer occupies a selected curved shell from 1.5 to 3.6 km. A
procedural three-dimensional density field is integrated using 16, 32 or 64
view-ray samples according to the graphics setting. The layer advances with
selected drift speeds on its own clock. Lighting uses the cloud-altitude optical
inputs, an approximate forward/backward phase mixture and a small number of
light-path density samples. The renderer integrates emission-like scattered
source contributions with exponential transmittance. Cloud multiple scattering
and reciprocal interaction with the molecular atmosphere are not solved.
The density controls and top lighting produce illustrative cloud structure;
small-scale microphysical or weather predictions should not be read from it.

A separate approximate flux modifier is used to light local surfaces. With
coverage C and optical-depth control τ, its conservative slab proxy is

    T = 1 / (1 + 0.75 (1 − g) τ),       g = 0.85,
    B = exp[−τ / max(sin a, 0.035)],
    F_diffuse' = F_diffuse (1 − C + C T)
                 + F_direct,horizontal C max(T − B, 0).

This expression is an explicit low-order graphics approximation. Direct ground
shading also samples a projected cloud shadow. It is not the diffuse hemisphere
integral of the rendered three-dimensional cloud field. The software renderer
uses a simpler spherical cloud sheet and coverage-averaged direct attenuation;
its detailed cloud shapes and local shadows are therefore different from the
WebGL path. Clear-sky maps and their numerical readouts are common to both.

Fog and haze add Beer–Lambert attenuation with extinction k = 3.912 / V, where V
is a chosen 2% contrast visibility distance. The in-scattering colour is derived
from the current calculated light with a display-only neutralization term.
The layer height is prescribed. Rayleigh aerial perspective on local objects is
also an approximation using surface scattering coefficients. Neither term is a
new fully coupled aerosol calculation or an empirical fog simulation.

Rain streaks and surface wetness are exposure and material effects. Their speed,
amount and wetness control are illustrative; they do not report rain rates or
terminal drop velocities. Earth and Moon use the same chosen weather controls
in comparison mode. The application labels the lux badges as the clear-sky
baseline, never as calibrated weather-altered ground illumination.

## Landscape and water

A deterministic procedural wooded cove supplies visual scale references.
Terrain combines smooth coastline and ridged height functions; 56 trees and
associated branches/crown clusters, rocks and shrubs supply 8,558 analytic
primitives. A BVH traces the fixed viewing point once into a 3072 × 1536 full
panorama of surface albedo, normals, distance/material, approximate horizon
shadows and local ambient occlusion. The visible rays then receive the changing
sky and solar light at interactive rates. This is a fixed-viewpoint scene, not
free navigation through a full three-dimensional world.

The scene is local tangent geometry, chosen for visual comparison. It does not
establish measured basin locations, full terrain-curvature effects, tree health
or organism dimensions. The precomputed sky and its ground boundary use a
smooth spherical planet independently of the local illustrative terrain.

Diffuse surface light is reconstructed from the sky's nine spherical-harmonic
coefficients. Direct light uses a selected soft-shadow horizon approximation;
canopy obstacles do not capture every open gap under a tree. Foliage includes a
small backlighting term. Spectral sky results are applied to RGB surface
reflectances. Local albedos do not feed back into the global sky calculation.

Water reflects the calculated sky and visual cloud layer, with Fresnel weighting,
a small body-colour term and a microfacet-like solar highlight. Several selected
wave modes use ω = sqrt(g |k|), with g = 1.62 or 9.80665 m/s², holding the chosen
wavenumbers fixed. Unresolved short waves are filtered by the projected sampling
footprint to reduce aliasing. Wave amplitudes and drift are selected display
parameters; water depth, wind-wave generation, flow and shoreline hydrodynamics
have not been solved.

## Exposure, night and precision

The fixed daylight display exposure is 1/8500 before user compensation. The
luminance-based tone curve is x/(1 + luminance(x)), followed by gamut compression
and sRGB encoding. Automatic exposure uses a smoothed clear-sky horizontal-light
reference, shared across the Earth/Moon comparison. It is a camera aid, not a
physiological eye model or absolute reproduction of outdoor brightness.

The Sun is the sole illumination source. Earthlight, stars, airglow and artificial
lights are omitted. Deep-night outputs inherit the coarser source grid and are
exploratory. Neither the viewer's display visibility nor an automatically exposed
night image establishes human visual acuity under those conditions.

The FP16 GPU sky packing retains a per-frame scale and records quantization
checks (`validation/packing_checks.json`). `test_viewer_data.py` checks the
display-model identities, clocks and landscape assets; browser controls and exact
Mesa shader renders are recorded in `validation/`. The solver's own tests and
noon Monte Carlo checks live with it in `illumination/sky/`. These address
selected numerical and implementation properties, not empirical lunar weather
validity.
