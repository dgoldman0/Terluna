# Atmospheric columns and clouds — revision 04

Baseline: `032765c51b1aeb140c732e8f484e47efda15cf4d`, shoreline revision 03,
Three.js r186. This revision is a local development checkpoint. It changes the
WebGL shoreline viewer and introduces three selectable atmospheric studies.
The TSL renderer lab retains its clear-sky implementation.

## Scope and evidence

Pressure, temperature, humidity and wind profiles now determine a parcel's
condensation level, accessible buoyant layer, optical condensate, cloud ceiling
and vertical wind shear. Those same data are uploaded to the renderer. Cloud
illumination follows physical altitude and spherical Sun paths. These are
specified sounding experiments. The soundings are inputs; the numerical cloud
boundaries and optical integrals are outputs. A general circulation model,
monthly radiative–convective equilibrium and a precipitation microphysics model
have yet to supply those inputs dynamically.

The existing molecular-sky atlas, surface assets, terrain, ecological placement
and hydrologic equations are preserved. **The thermodynamic sounding and the
inherited exponential molecular optical background remain separate models.**
The new altitude-dependent cloud-light calculation uses that fixed optical
background so the reference comparison retains its original sky. Regenerating
the spectral sky from the diagnosed density, water-vapour and ozone profiles is
a subsequent coupling step. The present viewer should not be interpreted as
having completed that coupling.

The authored four-hour rain replay is still available and is labeled separately.
Its older slab-cloud geometry remains a historical weather-rendering path. The
new column studies hold the surface rain input at zero and do not manufacture
rainfall from a diagnostic cloud-water content. Existing surface-water stores
therefore retain their established accounting.

## Sounding coordinates and hydrostatic height

`src/weather-column.js` uses metres, seconds, kelvin and pascals. Soundings are
specified against x = ln(p_surface / p), allowing the same temperature structure
in pressure coordinates to be compared between gravitational environments.
The Earth reference has p_surface = 101325 Pa, R = 6371000 m and g_surface =
9.80665 m/s². The lunar case has p_surface = 121590 Pa, R = 1737400 m and
g_surface = 1.62 m/s². These comparisons include the selected surface-pressure
difference as well as gravity; they are not gravity-only controlled experiments.

Virtual temperature includes water vapour. Hydrostatic geopotential is integrated
as dPhi/dx = R_d T_virtual. For spherical gravity,

    Phi(z) = g_surface R z / (R + z)
    z(Phi) = R Phi / (g_surface R - Phi)
    g(z) = g_surface [R / (R + z)]².

The local scale height is R_d T_virtual / g. The lunar surface values in these
soundings are about 51–53 km. The dry parcel lapse diagnostic is g/c_pd, about
1.61 K/km at the lunar surface. The inherited optical scale height remains
48.4 km, as explicitly distinguished above.

The default x step is 0.001, refined to 0.00005 near the surface to resolve the
fog layer. The exact condensation coordinate is inserted into the grid. Pressure
falls to exp(-2.2) of its surface value. The grid is ordered and deduplicated.
Pressure-drop versus integrated rho g dz, step refinement and the spherical
height inversion have independent numerical checks.

## Parcel thermodynamics and cloud support

Saturation vapour pressure follows the temperature-dependent latent-heat form
presented by Ambaum and documented in MetPy [1]. Liquid and ice saturation share
a triple-point reference of 273.16 K and 611.657 Pa. Constant heat capacities and
rounded gas constants are recorded in the source. MetPy's published example
uses a slightly different reference vapour pressure, so the comparison test
allows that explicit difference.

For unsaturated ascent the parcel retains its vapour mixing ratio and follows
T = T_initial exp(-kappa x), with moist-gas kappa. The lifting-condensation level
is the numerical root where the parcel reaches saturation. The root is checked
against the independently documented Romps/MetPy example [2].

Above condensation, a fourth-order Runge–Kutta integration follows the dilute
liquid pseudoadiabatic equation documented by MetPy [3]:

    dT/dx = -(R_d T + L_v r_s)
             / [c_pd + L_v² r_s epsilon / (R_d T²)].

Here L_v is held constant in the lapse equation, while saturation pressure uses
the temperature-dependent expression. An independent pressure/temperature
example is checked down to 200 hPa. The parcel is undiluted. Buoyancy uses virtual
temperature, B = g (Tv_parcel - Tv_environment) / Tv_environment. A selected small
launch kinetic energy plus integrated buoyancy determines accessibility. The
first neutral-buoyancy crossing after cloud initiation supplies the convective
ceiling. An LCL alone consequently does not guarantee a rendered cloud.

The reported positive-buoyancy integral is a parcel diagnostic. Entrainment,
condensate loading, turbulence, overshoot and environmental response can change
actual cloud depth. The three selected soundings are fully recorded in
`validation/clouds/scenarios.json` and exported by the viewer.

## Selected regimes and their calculated boundaries

| Study | Lunar cloud base–ceiling | Earth cloud base–ceiling | Specified origin |
| --- | --- | --- | --- |
| Coastal inversion fog | 0–0.123 km | 0–0.020 km | Saturated air cooled 0.8 K beneath an inversion |
| Elevated fair-weather cloud | 9.289–12.500 km | about 1.53–2.06 km | 294 K, 48% RH, 0.8 K parcel heating, capping inversion |
| Deep convective column | 3.982–58.697 km | about 0.66–9.83 km | 299 K, 78% RH, 1 K parcel heating, unstable layer and upper cap |

These heights describe the prescribed experiments, not the expected everyday
weather of a terraformed Moon. In particular, the approximately 59-km ceiling
is an undiluted parcel result for the selected deep sounding. It is not a
published prediction or a validated climate result. With parcel heating removed,
the 294-K, 48%-RH reference gives the previously discussed roughly 8.6-km lunar
condensation level; the fair-weather study explicitly adds heating.

The optical condensate is a selected retained fraction of cumulative parcel
condensation, reduced with ascent. Retention factors are 0.12 for the fair case
and 0.055 for convection. The deep case gives an unmasked liquid water path of
8.35 kg/m² and ice water path of 2.20 kg/m². The fair case gives 0.190 kg/m² of
liquid. These are diagnostic optical columns, not a conserved prognostic cloud
water budget. Horizontal morphology reduces their actual scene occupancy.

A continuous 273.15–233.15 K partition provides the optical liquid/ice fraction.
The buoyancy calculation continues using the liquid pseudoadiabat; ice latent
heat and deposition do not feed back on the parcel in this version. Effective
radii are specified per regime. Extinction uses the large-particle Q_ext ≈ 2
approximation, beta = 3 liquid_water_content/(2 rho_water r_effective), with an
analogous equivalent-sphere ice term. Ice crystal habits and halos are outside
this optical closure.

## Rendering the shared fields

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

## Verification and remaining work

The evidence manifest is `CLOUD_VALIDATION.json`. It separates the numerical
suite, actual GPU field readback, reviewed images, reference/reset comparisons,
intermediate failures and untested capabilities. A synthetic renderer checks
cache invalidation and state restoration; that check does not establish pixel
correctness. Pixel and float-target checks run in the real WebGL viewer.

The new capture harness treats any graphics-context-loss message as a failure,
including late console events. This is stricter than the historical harness.
Short successful scenarios do not establish sustained traversal stability.
Higher-resolution Balanced SwiftShader tests, including isolated 1536×768 sky
caches with Economy scene quality, still lost their contexts after caching.
Those records are retained as failures. Higher cloud-cache resolution can also be tested
independently while holding scene shadows and reflection resolution at Economy.
Consumer-GPU performance and native WebGPU cloud rendering remain unmeasured.
The TSL lab has an explicit notice directing column-cloud studies to the WebGL
viewer. No new node-cloud parity claim is made.

Priority follow-through is density-consistent spectral background regeneration,
entraining/mixed-phase parcel or cloud-resolving calculations, precipitation and
virga linked to surface input, calibrated cloud optical closures, longer GPU
stability tests and a shared TSL implementation.

## References used for equations and implementation checks

[1] Unidata, MetPy `saturation_vapor_pressure`, documenting the Ambaum (2020)
liquid and solid saturation expressions and an independent numerical example.
https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_vapor_pressure.html

[2] Unidata, MetPy `lcl`, documenting Romps (2017) and a pressure/temperature
example used here as an independent check of the root calculation.
https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.lcl.html

[3] Unidata, MetPy `moist_lapse`, including its governing equation and numerical
pressure-profile example. The present implementation was independently written.
https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.moist_lapse.html

[4] Eric Bruneton, Precomputed Atmospheric Scattering: reference implementation,
testing and radiance/colour conversion discussion.
https://ebruneton.github.io/precomputed_atmospheric_scattering/

All external reference material remains outside the redistributed source. The
existing Third-Party Notices and Three.js MIT license are retained.
