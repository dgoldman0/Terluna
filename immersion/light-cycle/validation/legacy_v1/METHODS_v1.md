# Methods and interpretation

## Controlled comparison

The experiment compares directional visible radiance at the same geometric Sun
elevation in two prescribed spherical atmospheres. Timing is applied afterward.
A constant equatorial path uses 15 degrees/hour for Earth and
360/(29.53059 × 24) degrees/hour for the lunar case. This isolates optical
appearance from elapsed time. It excludes latitude, season, ephemerides, terrain
horizons and libration. Sunrise reverses the optical sequence without predicting
weather differences between morning and evening.

Earth: radius 6,371 km; effective molecular density scale height 8 km; surface
optical density multiplier 1. Moon: radius 1,737.4 km; scale height
8 × 9.80665/1.62 = 48.428 km; density multiplier 1.2. The ratio of the idealized
vertical molecular columns is 7.264. Density is prescribed as exp(-h/H), and
scattering coefficients use the terrestrial molecular proxy at the stated
multiplier. These are convenient optical profiles. They are neither a solved
hydrostatic temperature profile with varying gravity nor the earlier thermal
column results in the repository. The 1.2 multiplier approximates the number
density change at comparable surface temperature; composition-specific refractive
indices and depolarization factors have not been substituted for the reference
air proxy.

The top is ten scale heights above the surface. Discarded overhead exponential
column fraction is exp(-10). Global optical depth accuracy is a separate issue
from this small vertical truncation. Ground reflectance is 0.1, Lambertian and
spectrally neutral. Observer height is 1.7 m. Visible shield transmission is 1:
an unfiltered control, not a claimed protection-system spectrum.

The Earth reference and first lunar scenario each contain a total 300 DU ozone
column. A triangle rises from zero at 10 km to one at 25 km and returns to zero
at 40 km on Earth. Lunar heights are stretched by H_M/H_E, with peak density
reduced to preserve the total column. The second lunar case has zero ozone.
These choices expose an important optical uncertainty. Neither establishes
lunar photochemical equilibrium or a biological UV environment.

## Scalar transport

At each wavelength the radiance equation is

    dI/ds = -chi I + beta_R(r) integral[p(d dot d') I(d') dOmega']

with chi = beta_R + beta_O3 and the scalar Rayleigh phase

    p(nu) = 3/(16 pi) (1 + nu^2).

The direct solar beam is attenuated by exp(-integral chi ds) along the actual
straight spherical ray. A planet-intersecting ray is occulted. Solar strips
span the visible part of a 32 arcminute uniform disk. Their integrated weights
are normalized to the analytic visible circular-segment area. Twelve strips
are used. The source phase and subsequent sky scattering use the disk-centre
direction; this is a small-angle approximation. Limb darkening and refractive
ray curvature are absent.

Successive scattering orders include the direct beam once and the preceding
diffuse order subsequently. Lambertian ground emission receives the direct
horizontal solar irradiance and diffuse downward irradiance of the appropriate
order. Each order's angular field is summarized by J = integral I dOmega and
K_ij = integral I d_i d_j dOmega. The Rayleigh source is exactly
3/(16 pi) [J + d_i K_ij d_j] for an unpolarized scalar field. This second-moment
representation is exact for this phase function; it does not isotropize the
scattering or make a diffusion approximation. Spatial interpolation and angular
quadrature do remain approximations.

Axial symmetry about the Sun permits a two-dimensional field indexed by radius
and local solar elevation. Symmetry removes K_rp and K_tp. Stored quantities
are J, K_rr, K_tt, K_pp, K_rt and downward horizontal diffuse flux. Rays are split
at concentric radial shells. Within each segment, extinction and the scattering
source are sampled at its midpoint; the exponential formal solution is exact
for those constant values. Multiple radial crossings are retained for tangent
rays. The solar beam uses a separate, finer lookup grid and direct slant-column
quadrature. It includes the planet shadow and ozone profile.

The standard source grid has 35 radial samples, 0.75-degree solar-angle spacing
through the twilight region, a split Gauss–Legendre zenith quadrature with 20
nodes, 16 azimuth samples, and 56 transport shells. Finite solar attenuation has
101 radial nodes, a 0.125-degree twilight grid and 240 slant-column integration
segments. The model includes 22 wavelengths at 20 nm spacing from 380 to 800 nm.
Iterations stop once the last order adds less than 0.0005 of downward diffuse
flux at every sampled surface wavelength and Sun angle from -24 to +12 degrees,
or at the requested order limit. This is an increment test, not a rigorous
bound on the sum of all omitted orders. Orders and residual increments are
stored per atmosphere. Deep-blue lunar paths need many more scattering orders
than the terrestrial reference.

Panoramas cover 126 Sun positions, -90 to +12 degrees. The -35 to +12 degree
interval has 0.5-degree spacing with extra 0.125-degree samples near the horizon;
the deeper-night extension has 3-degree output spacing. Below -35 degrees the
underlying diffuse source field is sampled every 5 degrees. The stopping and
spot-check criteria focused on earlier twilight, so this extension is explicitly
exploratory. Other night-light sources are absent. Each contains 65 elevations and
49 azimuths across the symmetry half-sky. Elevations are quadratically spaced
toward the horizon. A separate 64-shell ray integration evaluates spectra from
the solved source field. Smooth interpolation supplies viewing directions and
times between these samples. The Sun disk is drawn separately at its actual
angular size. These maps are appropriate for smooth molecular skies; they
would be insufficient to resolve detailed cloud structure.

## Spectra and display

Solar irradiance and ozone absorption coefficients are binned arrays from the
Bruneton reference implementation; see SOURCES.md and LICENSES.txt. Spectra
are integrated against analytic approximations to the CIE 1931 2-degree colour
matching functions (Wyman, Sloan & Shirley, 2013), including the 683 lm/W factor.
XYZ is converted to linear sRGB. Raw atlas RGB may contain negative components
for colours outside the sRGB gamut; physical luminance can still be recovered.
A display-only gamut mapping and common luminance tone curve map the result to
the monitor. FP16 atlas packing uses one scale per Sun frame to retain dynamic
range; quantization is checked against the full-precision display atlas.

The readout is horizontal photopic illuminance from direct Sun plus hemispheric
diffuse light. Its diffuse component is integrated over the panorama angular
grid. It is a model output, not a measured lux reading. The fixed-exposure mode
holds the exposure constant across worlds and time. Automatic exposure is a
shared camera-style aid. It is not a physiological dark-adaptation model, and
can make very dark conditions look far brighter than their lux readout.

## Foreground is a reference, with simpler shading

The distant ground and sky come from the spherical calculation. A few nearby
ray-traced reference objects establish scale and provide surfaces and shadows.
Their diffuse environment lighting is reconstructed using nine spherical-
harmonic coefficients; this deliberately smooths the angular illumination.
Local object materials are RGB approximations. Shadows use a soft visibility
approximation. A short-distance molecular aerial-perspective approximation and
fade into the physical ground are applied. These objects are visual aids;
the quantitative sky field and irradiance records are the experimental outputs.
There is no simulated terrain, vegetation, ocean, city, weather or acoustic scene.

## What the comparison establishes

It demonstrates consequences of these specified optical profiles under the
stated scalar straight-ray model. Matched Sun elevations isolate optical
changes, while the two clocks expose the duration change. The zero-ozone switch
shows a scenario dependency rather than an error bar around a preferred result.

External validation against measured Earth twilight, vector/polarized transport,
refraction, denser spectral sampling, cloud and aerosol scenarios, physically
solved ozone and temperature profiles, and Earth illumination are future steps.
The present Monte Carlo checks test the same scalar physics with independent
spatial sampling. They cannot establish whether the chosen atmosphere is
physically achievable or whether the model includes every visually important
process. A head-mounted display would also require separate display calibration
and comfort testing; this version is a mouse/touch 360-degree browser viewer.

## Actual validation outcome

Twenty-one unit tests pass. Twenty-four selected direction/wavelength tests
used 500,000 backward Monte Carlo paths each. Differences from the deterministic
calculation are below 0.7% in the selected lunar cases; Earth differences reach
about 7%. Some terrestrial differences exceed the Monte Carlo sampling errors.
These observations are not whole-sky error bounds. Doubling output-ray shell
resolution from 64 to 128 changes selected radiances by at most 0.65% above the
stated radiance threshold. An earlier combined fine-grid/solar-quadrature/order-
stopping sensitivity run is archived with its distinct conditions.

Whole-sphere energy accounting at 440, 560 and 680 nm leaves residuals of up to
0.27% for Earth and 4.6% for the lunar profiles. The present model and its audit
quadrature therefore retain a few-percent global numerical closure issue.
No empirical Earth-twilight benchmark has been completed. Refer to the actual
JSON records rather than interpreting a test count as physical validation.

The browser was exercised in Chromium with WebGL unavailable. The software
fallback, loaded from the same self-contained HTML, performed the 360-degree
rendering and control checks without external network requests or JavaScript
errors. The exact GPU vertex and fragment shaders separately compiled and linked
in Mesa OpenGL ES 3.2. Browser WebGL execution itself was unavailable in this
environment. The software foreground omits the thin pole's cast shadow; sky
sampling, clock, colour conversion, and exposure use the same numerical data.
