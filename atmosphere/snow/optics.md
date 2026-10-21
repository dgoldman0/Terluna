# Snow light in the aerial winter

The present snow-population outputs are exploratory numerical distributions. Their uncapped giant tails and
representation-dependent particle arithmetic prevent using them as reliable storm forecasts; the
[numerical audit](results/numerical_audit.json) separates those failures from the physical questions about
breakup. The optical calculation exposes which size ranges contribute to the result, retaining the complete
area and mass ledger.

A snow population determines how much projected ice area a ray encounters. That area, together with the
light reaching the winter layer, determines how a storm obscures a distant view and adds scattered light along
the sightline. [optics.py](optics.py) calculates those quantities from the population distributions and the existing
solved molecular sky. The results are conditional optical calculations; they do not establish an observed lunar
visibility or a photographically resolved appearance of individual flakes.

An important consequence of the aggregate mass law is already clear. With mass proportional to diameter to the
power 1.9 and projected area proportional to diameter squared, area per unit ice mass changes only as diameter
to the power 0.1. Increasing diameter fourfold raises that ratio by about 15%. Giant loose flakes can therefore
carry substantial total projected area even when their number decreases. Their slow settling increases the ice
concentration needed to carry a given snowfall flux. Large flakes alone do not imply a transparent storm.

## Distribution and extinction

The [population product](results/population.json) and the separate [growth-depth sensitivity](results/growth_depth.json)
supply both number concentration and number-weighted particle
area in each size bin. The optical calculation retains their sum exactly:

\[
\beta_{\rm snow}=Q_{\rm ext}\sum_i n_i\overline{A}_i,
\qquad T(x)=\exp(-\beta_{\rm snow}x).
\]

Here each \(n_i\) is the whole bin's number per cubic metre, rather than a density per unit diameter. At the
freezing level the population producer converts the arriving number flux to concentration using each particle's
own settling speed; inside the cloud it samples the actual particles occupying each height layer. Reconstructing
the area from a representative diameter would lose the separate area moment, particularly across broad bins.
The reported attenuation treats the time-mean population as a homogeneous medium. Because transmission is
nonlinear in optical depth, exp(−mean optical depth) is not the time-mean transmission through intermittent snow;
the results do not give the frequency or duration of poor visibility.

The uncapped population law can also produce extreme numerical tails far outside physically supported
aggregate sizes. The optical sum retains all their area, with separate 1–10 m and at-least-10 m contribution
bins. Where those tails contribute appreciably, the total describes the runaway population model; it cannot
be interpreted as a credible storm prediction. No tail is silently clipped or recast as a physical snowflake.
A separate diagnostic reports the existing below-1-m component, leaving all its numbers and areas unchanged.
Excluded tail mass is not redistributed; this component calculation cannot stand in for a physical breakup
model, and its smaller bins still inherit the particle algorithm's numerical limitations. Mass and optical shares are both reported by size to distinguish a rare-count tail from a negligible
mass or area contribution.

For the central corrected-mass, dendrite-area cases, the freezing-level optical diagnostic is:

| Imposed vapour supply, mm/h | Arriving snowfall, mm/h | Initial crystals/L | Full 5% beam distance, m | Extinction from sizes ≥1 m | Below-1-m component distance, m |
|---:|---:|---:|---:|---:|---:|
| 0.3 | 0.298 | 10 | 813 | 0% | 813 |
| 0.3 | 0.347 | 100 | 698 | 0% | 698 |
| 1 | 0.945 | 10 | 243 | 0% | 243 |
| 1 | 0.995 | 100 | 229 | 0% | 229 |
| 3 | 2.695 | 10 | 79 | 2.49% | 81 |
| 3 | 2.775 | 100 | 44 | 45.79% | 82 |

These are homogeneous optical diagnostics of the inherited numerical populations, not plausible-weather
forecasts. In the last case, the at-least-10-m tail alone supplies **43.37% of extinction**, although its
represented particle count is tiny. Its full 44 m optical distance is therefore materially contaminated by
nonphysical aggregates. Two other random seeds give full distances of 76 m and 76 m, while their existing
sub-metre components give 80–81 m. The apparent agreement of those smaller components does not validate
them: the numerical audit identifies representation-dependent particle arithmetic and batch-collision
sampling limitations that also affect ordinary sizes. Changing only the representative area closure can
make the extreme tail dominate almost all extinction.

A coupled numerical refinement doubles the vertical cells from 30 to 60, raises injection sampling from
8 to 12 packets per step, and halves the maximum step from 30 to 15 seconds. Its full distance is 76.8 m,
with 1.32% of extinction from sizes at least 10 m; its unchanged sub-metre component gives 78.4 m. The
baseline values are 44.2 m, 43.37% and 81.6 m, respectively. This reduces the tail's optical contribution
in that realization, but the refined population still contains a nonphysical kilometre-scale outlier. One
coupled refinement and the seed spread establish sensitivity, not numerical convergence or a validated
storm prediction.

The separate shallow-cloud experiments hold the 100/L input and the same representative geometry, but change
the saturated growth path. Their model outputs are:

| Cloud-top temperature | Imposed vapour supply, mm/h | Arriving snowfall, mm/h | Full 5% beam distance, m |
|---:|---:|---:|---:|
| −5 °C | 1 | 0.801 | 302 |
| −5 °C | 3 | 1.660 | 138 |
| −10 °C | 1 | 0.943 | 248 |
| −10 °C | 3 | 2.402 | 90 |

The last row has about 0.11% of optical extinction above 1 m; the other three have none. These calculations
show the consequence of the specified growth path within this model, including its changed uptake of the
prescribed vapour supply. They do not establish cloud-depth frequencies, storm optical ranges or a validated
repair of the deep-column model.

The population's settling projected area is used for every optical ray direction. This direction-independent
area closure does not describe horizontally aligned flattened flakes, whose side-view and downward-view
cross sections can differ. The present population supplies no orientation distribution to calculate that
difference.

The geometric-optics closure uses \(Q_{\rm ext}=2\) and conservative, wavelength-independent scattering. This
counts the forward diffraction contribution as extinction from the original collimated beam. A finite angular
receiver can recover some of that light. The present calculation does not resolve electromagnetic scattering
from a porous aggregate or the observer's angular response. It also excludes absorption by impurities and weak
visible ice absorption. A separate \(Q_{\rm ext}=1\) sensitivity counts geometric interception alone and doubles
the corresponding 5% beam distance. This is a change in optical closure, not a confidence interval or a
prediction of how much forward-scattered light a particular observer recovers.

The product reports the distance at which the unscattered beam falls to 5%, \(-\ln(0.05)/\beta\), alongside
transmission over 10, 100 and 1,000 m. This follows the optical-range definition in the
[WMO instrument guide, chapter 9](https://www.weather.gov/media/epz/mesonet/CWOP-WMO8.pdf).
Interpreting it as a person's visual range requires additional conditions: a specified target and background,
their illumination, and the observer's adaptation. Those conditions are especially consequential at night and
in a spatially uneven snowstorm.

The snow cloud also contains whatever liquid water the atmospheric state supplies. The snow population assumes
liquid saturation while prescribing its vapour supply; it does not predict the liquid-droplet population. As a
separate sensitivity, droplets with 10 µm effective radius give:

| Imposed liquid water | Droplet-only 5% beam distance |
|---|---:|
| 0.01 g/m³ | 1,993 m |
| 0.1 g/m³ | 199 m |
| 0.5 g/m³ | 40 m |

These values follow \(\beta=3\mathrm{LWC}/(2\rho_w r_{\rm eff})\), with the existing atmosphere model's water density.
They are imposed cloud cases, not inferred frequencies or measured liquid-water contents. Their extinction adds
to the snow's. A community beneath the liquid cloud and one embedded in it can consequently have very different
views of the same falling aggregates.

## Light at snow altitude

The light calculation reads the existing versioned solved-sky product and checks the hash chain from its
directional atlas to the stored spectral radiation moments. It interpolates those moments in physical radius
and solar elevation, then independently traces the finite Sun through the product's spherical optical shells.
It does not carry surface illuminance upward. At 35 km, the sampled clear molecular atmosphere gives:

| Sun elevation relative to local horizontal | Direct normal illuminance | Diffuse light onto a horizontal upper-facing surface |
|---|---:|---:|
| +60° | 80,399 lux | 23,002 lux |
| +10° | 28,504 lux | 11,267 lux |
| 0° | 8,119 lux | 5,855 lux |
| −5° | 2,089 lux | 3,327 lux |
| −10° | 164 lux | 1,533 lux |
| −20° | 0 lux | 216 lux |
| −40° | 0 lux | 4.26 lux |

The direct column is normal to the solar beam; it should not be added directly to the horizontal diffuse
column. Above the ground, a source below the local horizontal can remain above the depressed solid-body
horizon. At −5° this still leaves a strong, spectrally reddened direct beam in the reference atmosphere; at −20°
the Sun is blocked from this height while molecular multiple scattering still supplies diffuse light.

These are incident clear-sky values before adding the snow and cloud. The optical atmosphere has a 294.9 K
surface and the admitted titania-stack transmission multiplied by 0.95, whereas the snow column uses the 301 K
resort interpolation. The product records this profile mismatch. Both source and cloud extinction must be
coupled before these numbers can describe the illumination deep inside an actual storm.

A separately labelled full-Earth proxy uses the Earthlight domain's rounded 7.8 lux above-air normalization,
with the source at +25°. At 35 km it yields 3.40 lux direct normal and 1.06 lux diffuse downward. Its spectrum is
the shielded-solar reference scaled to that brightness; the Earth's bluer spectrum, phase, changing continents
and larger disk remain outside the proxy. It is a separate source contribution, not an Earthlight colour prediction or a complete night sky. Solar
twilight must be added at a compatible Sun–Earth geometry; at full Earth the Sun lies roughly opposite it.

## Direction, colour and repeated scattering

For a short sightline, the added snow radiance is calculated from the direct beam and the diffuse field's
all-direction integral. The diffuse angular distribution is replaced by its isotropic equivalent. The direct
beam uses a normalized Henyey–Greenstein phase function, with asymmetry parameters 0.7, 0.8 and 0.9 as explicit
sensitivities. The calculation accepts optical depths only through 0.1; the exported examples use 0.05.
Incoming light is held constant along that short path. This makes the assumptions about the surrounding
storm visible and prevents the thin-path formula from being used as a thick-cloud renderer. The examples
describe exposed thin snow with approximately clear incident illumination; a short sightline deep inside an
opaque cloud would receive a different incident field.

At 35 km with the Sun on the local horizontal, the central asymmetry case gives 305 cd/m² added radiance when
the scattering angle is 20°, 85 cd/m² at 90°, and 82 cd/m² at 160°. The corresponding 20° light is much warmer:
its CIE chromaticity is about (0.495, 0.404), compared with (0.356, 0.362) at 90°. In this approximation,
the reddened direct beam dominates the forward view while the cooler diffuse field contributes more strongly
to the side view. These are conditional radiances of a thin population, not luminances of individual flakes or
a prediction of their specular flashes. The product also computes the contrast of an explicitly imposed black
target against a background whose unobstructed luminance equals the diffuse mean. Actual background radiance
and visual adaptation need their own scene.

Repeated scattering changes the interpretation at larger optical depth. A separate forward Monte Carlo follows
photons through a homogeneous, horizontally infinite snow slab until they leave the top or bottom. Every
scattering order is retained. With asymmetry 0.8 and a normal incident beam, a slab of optical depth 10 transmits
about half the total flux, although only \(\exp(-10)\), about 0.0045%, remains in the original unscattered beam.
Diffuse transmitted light can therefore illuminate a snowy space even when very little of the original
collimated beam survives. This flux calculation does not determine how much image detail survives, because
that also depends on the angular redistribution and the target and receiver geometry. This slab has no reflecting ground, surrounding molecular air, liquid cloud or finite storm
edge, so its flux fractions do not supply the missing three-dimensional storm radiance. The product records
sampling errors and intervals, including finite upper intervals when a rare event has no sampled photons.

The HG sensitivities cannot make ice halos or resolve facets. The primary
[Yang et al. ice-optics library description](https://www.giss.nasa.gov/pubs/abs/ya07100h.html) distinguishes habit,
roughness and orientation, and covers sizes through 10 mm. Applying a phase matrix to decimetre aggregates
requires a justified extension. Crystal form and orientation can eventually identify which halo or glint
mechanisms deserve calculation, but an aggregate's large span alone does not determine their strength.

## Reproduction and evidence

```sh
OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.optics
OPENBLAS_NUM_THREADS=1 python -m pytest atmosphere/tests/test_snow_optics.py
```

The runner needs both population products with their numerical assessments attached, the completed numerical audit,
and the admitted ignored sky archives. The assessments must bind the supplied audit. Missing files or
mismatched hashes cause failure; surface light or an invented spectrum is never substituted. The compact
[optics product](results/optics.json) carries model and input hashes, units, reading rules and assumptions.
[Source notes](optics_sources.json) distinguish inspected primary material from inaccessible full-paper links.

Tests check conserved area under rebinning, normalized phase functions, exact empty and absorbing limits,
shell geometry, isotropic-light colour, Monte Carlo energy conservation and unscattered Beer–Lambert transport.
The independent direct-shell readout agrees with the existing surface product to 0.012% or better at the
four sampled elevations. Reading the independently pinned coarse spatial grid at the same winter heights changes
diffuse scalar and horizontal light by at most 1.68% across the sampled conditions. Refining the solar-disk
quadrature from 32 to 128 strips changes direct light by at most 0.029% in four selected grazing or low-Sun
cases at winter heights. These checks establish implementation behaviour and consistency with the saved
optical field. The remaining physical work includes jointly evolving liquid and frozen cloud populations,
breakup-limited aggregate geometry, applicable ice phase matrices, and transport through a finite storm.
