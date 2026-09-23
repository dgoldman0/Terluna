# Shoreline methods and evidence boundary

## Purpose and scope

This implementation develops a navigable local experience from the preceding
full-cycle viewer. The scene is hypothetical: no surveyed basin, planted
community, validated architecture, climate prediction or physiological
low-gravity gait is asserted. Input values and optical results remain separable
from geometric detail and artistic choices.

## Inherited atmospheric data

The source is `Open_Moon_Full_Cycle_Source.zip`, SHA-256
`e753eaf67d3392fa48617fddad12376c3564ffb7cd10d2bd2705d4640e3a614b`.
`data/atmosphere.json` retains each source NPZ filename and SHA-256. The three
profiles cover Earth, Open Moon with 300 DU ozone, and Open Moon with zero ozone.
They contain 174 solar elevations from −90° to +90°, direct and diffuse ground
illumination, and cloud-altitude reference illumination. This release uses those
already-calculated results without modifying or rerunning the optical solver.

The original assumptions include radii 6,371 km / 1,737.4 km, exponential molecular
scale heights 8 km / about 48.428 km, and surface-density multipliers 1 / 1.2.
The ozone-bearing lunar case stretches a 300 DU profile. The incident solar
spectrum is the original unfiltered reference. These profiles are optical proxies;
they do not close the project's thermal structure, composition or escape model.
The prior up-to-4.6% global energy-accounting discrepancy remains unresolved.
Deep twilight inherits a coarser source field and its lower verification level.

The import selects every other view angle, reducing each 81×65 map to 41×33.
Signed linear RGB outside the positive gamut is desaturated toward luminance,
then stored as positive FP16 values normalized by each frame's peak. The largest
quantization change is below 0.00025 of the frame maximum for each profile. That
quantity does not bound angular downsampling, color-gamut projection, spectral
or underlying atmospheric errors. The original ground irradiance arrays are
retained separately for the clear-sky numerical badge.

The renderer interpolates between solar elevations and the quadratically spaced
view-elevation/azimuth atlas. Radiance and irradiance are scaled by 1/8500 within
the scene. Three.js performs material shading, environment-map filtering, ACES
camera tone mapping and sRGB output. Exposure settings therefore affect display,
not the physical-unit clear-sky badge.

The atlas remains the original fixed-observer clear atmosphere even when the
viewer walks to the overlook. Local parallax, terrain occlusion and selected
cloud-shell ray origins follow the actual camera. Altitude-dependent molecular
transport is not newly solved. Refraction, full gas/aerosol optical properties,
Earthlight, stars, airglow and actual Sun–Earth–Moon ephemerides remain absent.

## Geometry and solar time

The synthetic scene has a navigable area within approximately ±130 m in x,
with an irregular coast, inland wooded trail, a shelter at (20,33) m and an
approximately 20 m-scale overlook near (−56,47) m. Instanced rocks, branches,
leaves and grasses have real geometric depth. Nearby water has finer tessellation
than the distant curved water surface. Landscape/asset distribution is seeded.

The solar vector is `(sin(2πp), cos(2πp), 0)` for cycle phase p. It represents an
idealized equatorial zero-declination day. The selected solar period is 29.53059
Earth days or 24 hours. The geometric Sun radius is 0.00465421 radians. Glare is
an additional display halo. Horizon/refraction corrections and actual season or
latitude are not solved. A simple curvature sag is used in local surfaces, with
the same built foreground retained for Earth comparison.

Walking uses a ground-following camera and collision checks for shoreline,
large rocks and trunks. The deck/ramp is included in the walking surface. The
Space key launches a selected 1.2 m/s vertical hop with the chosen world's gravity;
no biological gait, balance or comfort result is inferred from it. The observer
controls operate in real time even while the environment is surveyed faster.

## Authored environmental forcing

The scenario has a four-hour prescribed passage through gathering cloud, rain,
clearing and retained wetness. `scenario.json` contains the keyframes. Temperature,
humidity, wind, cloud water content, cloud base/thickness, visibility and rainfall
are smoothly interpolated. Rainfall peaks at 10 mm/h; effective droplet radius
inside the cloud is chosen as 12 μm. These are selected inputs, not predictions
of lunar meteorology or event frequency. Moisture/energy/momentum closure for the
air column is not solved, and the episode does not change with the solar phase.
Consequently the user can intentionally inspect the same forcing at noon or night.

For a liquid-water cloud, the optical-depth proxy is

    tau = 3 LWC H / (2 rho_water r_eff),

with LWC in kg/m³, thickness H in m, water density 1000 kg/m³, and effective radius
in m. This uses a geometric-optics extinction efficiency of approximately two.
A procedural density field shapes the cloud; tau is its reference normalization,
not a guaranteed column-integrated water amount for every rendered pixel.

The cloud field is advected with integrated wind displacement. The same GLSL
function is sampled by the visible sky ray march and the surface direct-shadow
calculation. Cloud ray origins use the camera height. Cloud source illumination
uses inherited 2.5 km reference inputs although this episode places cloud layers
at other heights. Its finite-step lighting and multiple-scattering treatment are
approximate. Environment-map filtering conveys the rendered sky/cloud field to
materials. This makes inputs and selected geometry common without establishing a
fully reciprocal cloud–molecular–terrain transport solution or an energy bound.

Prescribed visibility supplies `k = 3.912 / V` for a 2% contrast extinction proxy.
A selected ambient in-scattering color is used in both the sky and objects.
Even the clearest preset retains a chosen 35 km visibility. This is a display
haze/fog envelope, not new aerosol microphysics; the clear lux badge is unchanged.

## Persistent water accounting

Representative surface reservoirs use the exact bounded constant-forcing equation

    dS/dt = I − (ke + kd) S,        0 ≤ S ≤ capacity.

The solver explicitly integrates evaporation `ke S`, drainage `kd S`, and overflow
when the cap is reached. The exposed surface holds at most 2 mm. A canopy-ground
store holds 1.4 mm, with a 0.35 mm interception store on leaves. Forty-five percent
of rain enters interception; the rest reaches canopy ground directly. Drainage
and overflow from leaves become throughfall. Sheltered ground receives no direct
rain. Evaporation coefficients respond to the chosen humidity, temperature and
wind, using an explicit empirical-style proxy rather than a surface energy balance.
Drainage coefficients and capacities are selected scenario inputs.

The ledger represents an exposed one-m² reference column and a second one-m²
canopy column. Thus adding their water depths is equivalent to summing litres
from the two reference columns; it is not an area-weighted terrain rainfall total.
Interception is internal to the canopy column. A zero-input sheltered store is
reported as a separate reference. The numerical residual compares cumulative
input with remaining stores, evaporation and external drainage/overflow.

Forcing is advanced in substeps of at most ten seconds. Within each substep,
throughfall from the leaf reservoir is supplied at its interval-mean rate to the
ground reservoir. This split is conservative; timestep refinement tests its
storage approximation. The code is exact for each scalar constant-forcing
reservoir but not an exact solution to time-varying coupled canopy dynamics.
The surface shaders interpolate representative exposed/canopy wetness spatially.
Wet-film patches are decorative receivers of that state; the terrain does not
solve lateral runoff or actual puddle depth/area.

In optical-study mode, selected stores remain fixed for comparison. In live
mode they evolve, including drying after the four-hour episode reaches its last
forcing state. Selecting an episode timestamp reconstructs its preceding history
from dry initial conditions; it is a repeatable scenario snapshot.

## Rain, waves and sound

The chosen visible raindrop radius is 0.7 mm. Terminal speed balances spherical
buoyancy-corrected gravity and drag, using Schiller–Naumann below Reynolds number
1000 and Cd=0.44 above. Selected air density is 1.45 kg/m³ for the Moon and
1.2 kg/m³ for Earth; dynamic viscosity is 1.8e−5 Pa s. Deformation, drop-size
distribution, acceleration, evaporation and breakup are omitted. Rain streaks
subsample this speed and the chosen rain amount; their particle number is not a
calibrated volumetric rain distribution. Roof masking suppresses drops beneath
its footprint. Trunk/leaf interception does not individually collide with particles.

Four selected water-wave modes use `omega² = g k tanh(kh)` at a prescribed depth
h=12 m. Wave amplitudes respond to wind with a selected scaling; there is no
wind-wave growth, breaking, fetch, tidal or shoreline-flow solution. Smaller
normal detail and foam are graphics approximations. A planar reflected camera
supplies reflection, with Fresnel weighting and a microfacet-style solar glint.
The sea curvature and wave perturbations make planar reflection approximate.

Web Audio synthesizes positioned surf, foliage, rain, roof and drip textures.
Camera translation/orientation updates the listener, while rain/wind/water stores
control amplitudes. Shelter exposure changes the mix. Steps vary between wood and
ground and respond illustratively to wetness. These are procedural audio cues,
not recorded or predicted lunar acoustics. Accelerated live modes mute sound.

## Validation level

`VALIDATION.json` records the checks actually performed. Model unit tests,
source/packing integrity and isolated GLSL compilation are separate from
integrated Three.js browser execution and physical validation. The latter two
are unperformed. The previous scientific solver was not changed in this work.
