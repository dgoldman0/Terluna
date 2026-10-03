# High clouds through the long twilight

The remote baseline used committed cloud summaries and packed sky data to calculate sunlight reaching elevated features after surface sunset. The current evidence favours studying illuminated sides, rims and thin upper structures against a luminous twilight sky. Their appearance depends on cloud scattering and the surrounding light field.

The [local follow-up](#local-dusk-columns) now recovers the raw dusk condensate and particle fields, with explicit second-cycle sampling. It adds liquid-extinction scenarios while preserving the original mechanism screen below. Both stages reuse completed climate simulations.

## What the corrected cloud run supplies

The committed `ring_ring_equator` products use corrected `A28_dim5_moon`
forcing: model days 29.5–59, 237 snapshots, 6-km horizontal columns, and a
two-dimensional equatorial ring with prescribed land moisture and upper forcing.
Optical paths use the global reference column. Simultaneous storm soundings
remain a further input to local transport.

Cloud-top statistics refer to the highest cloud-liquid/ice cell over each
detected rain footprint in a snapshot. Cloud is `qc + qi >= 1e-5 kg/kg`; a rain
footprint is a contiguous region above 1 mm/h. The resulting percentiles describe these detected footprints at sampled times.

| Storm-footprint statistic | Cloud-top height | Overhead solar-centre sunset after surface sunset | Direct illumination falls below 1,000 lux |
|---|---:|---:|---:|
| Median | 22.1 km | 17.9 h | 10.5 h |
| 90th percentile | 37.8 km | 23.3 h | 15.5 h |
| Largest sampled value | 83.8 km | 34.3 h | 27.2 h |

The last column is a molecular-only, design-shield calculation of illumination
normal to the beam. **1,000 lux is an arbitrary reference level**. Cloud visibility requires a
contrast and observer model. Geometric durations assume fixed height and position at
equatorial equinox. The finite solar disk adds about 0.52 h between centre and
last-limb geometric sunset; refraction and terrain are omitted.

The corrected heights are lower than the old uncorrected ring's often-quoted
28-km median and 68-km 90th percentile. The newer statistics are used throughout.

## What stands out

**There is a long opportunity for elevated sunlight, but the brightness evolves
through it.** At surface sunset an overhead 40-km feature receives about
10,800 lux of direct normal illumination. At a Sun depression of 5 degrees
(9.84 h later), it receives 3,290 lux; at 10 degrees (19.69 h), 387 lux. Its
solar centre is geometrically visible until 24.0 h. These numbers describe
incident light. For a Sun below the cloud’s horizontal plane, direct light reaches exposed
sides or undersides; a horizontal upper face faces away from that beam.

**The sunsetward view deserves the first detailed scene.** With the observer's
Sun 5 degrees down, a 40-km feature 100 km west is 19.9 degrees above the visible
horizontal. The Sun at that feature is only 1.70 degrees down. Its incident
direct normal illumination is 7,760 lux. The cloud-to-observer ray then removes
more light: the spectrally combined Sun-to-cloud-to-eye source factor is 2.70%
of the above-atmosphere design source, compared with 1.97% overhead. Among the
sampled westward ranges, this factor peaks near 75 km at that time. At later
twilight the favourable range moves farther west. This maximum describes the two-path transmission. Cloud brightness also
depends on scattering and viewing geometry.

The same 100-km example has a 24.9-degree scattering angle, where zero denotes
forward scattering. A resolved cloud calculation must use the appropriate
liquid/ice phase functions and internal optical depth. The angle is retained as an input to that further scattering calculation. A 40-km feature 400 km away is hidden by
the solid body from the assumed 1.6-m observer, even though the feature itself
can receive sunlight.

**The direct component becomes very red.** In that example, solar transmission
to the feature is 0.074%, 4.21% and 20.1% at 450, 550 and 650 nm. Including the
viewing path leaves 0.0061%, 1.41% and 11.5%. This spectral filtering supports
reddening of the direct component; perceived cloud colour also depends on
diffuse illumination, cloud optics and atmospheric light along the sightline.

**The existing sky model retains a luminous twilight.** The
unfiltered, zero-ozone proxy retains about 3,220 lux on a horizontal surface
when the Sun is 5 degrees down, 1,650 lux at 10 degrees, and 772 lux at 15 degrees.
Those are a separate atmospheric calculation. The study also computes direct
cloud beams in that same proxy: at 40 km overhead, about 3,540 lux at 5 degrees
and 296 lux at 10 degrees. Even this internally matched comparison is of two
different illumination planes and components. Cloud-to-sky luminance contrast
requires the cloud scattering source and a background sightline.
The current shielded column still needs a spherical diffuse-sky solution.

**Very high dusk clouds are a special case in this run.** In the 90–100-degree
local-hour-angle bin, the largest cloud occurrence at any sampled height from
20 to 40 km is 5.54%; from 40 to 60 km it is 0.11%. These fractions describe sampled columns at individual heights. Observer
sky coverage and visibility depend on the three-dimensional geometry. The
60–80 km scenarios measure height sensitivity within an intermittent population.

The rain-footprint tracking diagnostic has a median lifetime of 6 h and a
90th percentile of 15 h, at 3-h output spacing; duration is the number of
detected samples times that spacing. Anvil lifetimes require separate tracking.
Cloud-shape persistence over an 18–34 h illumination window remains open.

## Calculation and evidence boundary

For radius `R`, height `h`, westward surface arc distance `s`, and observer Sun
depression `d`, the overhead horizon depression is `acos(R/(R+h))` and the
cloud's local solar elevation is `s/R - d` (radians). Its solar-centre sunset
occurs at observer depression `s/R + acos(R/(R+h))`. Time follows the shared
29.53059-day synodic period at equatorial equinox.

The domain model intersects straight rays with concentric shells analytically
and applies wavelength-dependent Beer–Lambert extinction. The current profile
is rebuilt by `illumination.surface_light.model.build` from committed solved
columns, preserving its per-layer Rayleigh depth exactly; total vertical
Rayleigh depth at 550 nm is 0.7567456. Extinction is constant inside each stored
layer, which is a vertical-resolution approximation. Ozone is included through
the existing coarse cross sections and is negligible in this Moon column.

The 48-point, 360–830 nm solar spectrum and photopic matching fit are inherited
from the sky domain’s binned ASTM source. The full surface-light calculation
uses WHI, which was absent from the remote workspace and is now restored locally. The design source applies the stored titania-film
transmission and the additional 0.95 factor. An unfiltered current-column case
isolates the source change; a separate unfiltered exponential case isolates
the atmosphere used by the packed sky atlas.

The solar disk uses 32 vertical strips weighted for a uniform circular disk.
This small-angle approximation omits transverse curvature, limb darkening and
refraction. The calculation includes solid-body solar and viewing occultation.
The two-path factor is
`sum(wavelength_weight * solar_transmission * viewing_transmission) / sum(weight)`.
Multiplying separately integrated photopic transmissions would give 2.22%
instead of 2.70% for the 100-km example and would lose the shared spectral
filtering. The factor describes molecular transmission along the two paths;
cloud scattering is a further operation.

Omitted: water/oxygen and other absorption bands, aerosol and cloud extinction,
cloud self-shadowing and multiple scattering, atmospheric in-scattering along
the viewing ray, refraction, polarization, terrain and a simultaneous local
atmospheric profile. Additional extinction reduces the corresponding ballistic
component; diffuse light and refraction require their own treatment. Total visible cloud radiance requires all of these transport components.

Numerical checks cover vertical/tangent analytic limits, reciprocal finite
paths, path composition, ground blockage, the half-visible finite Sun in
vacuum, and independent continuous integration of the exponential proxy.
For the sampled checks, 32 versus 64 solar strips differ by at most 1.8 lux;
1-km versus 0.5-km proxy shells differ by at most 0.12% where the beam exceeds
1 lux. These checks concern numerical integration. Coarse profiles, omitted optics
and the two-dimensional climate model contribute additional physical uncertainty.

## Local inputs identified by the remote screen

The remote screen identified a **specific dusk scene** as the next physical
calculation, with upper structures around 20–40 km and observers tens to a few
hundred kilometres east. Predicting illuminated surfaces, colour contrast and
duration required raw fields beyond the compact products.

The identified local input was a short sequence of corrected CM1 dusk snapshots
with geometry/time, density and temperature, separate `qc/qi/qr/qs/qg`, and
particle-number or size information where available. Missing microphysical
sizes and ice habits would still require stated optical assumptions. Begin
with a 2-D diagnostic if only a ring is available; a credible scene needs
3-D structure or a clearly labelled extrusion sensitivity. The ignored
`climate/crm/products/ring_ring_equator.npz` may help select events, but raw
CM1 files are needed for fields that the reduction discarded.

The saved maximum-rain storm slice peaks at local hour angle 21.9 degrees,
well before sunset. Its condensate combines liquid, ice, snow and graupel.
Its time and combined condensate limit its use for dusk extinction.
The atlas’s stored cloud illumination applies at 2.5 km; the elevated-feature
calculation evaluates 20–80 km explicitly. Native sky archives improve angular
background sampling. Cloud optics need the separate microphysical fields
recovered in the local follow-up.

## Reproduce and inspect

```bash
OPENBLAS_NUM_THREADS=1 python -m research.studies.cloud_twilight.run
python -m pytest illumination/cloud_light research/studies/cloud_twilight -q
MPLCONFIGDIR=/tmp/terluna-mpl python visualization/cloud-twilight/plot.py
```

The versioned [result](results/cloud_twilight.json) carries source/input SHA-256
hashes, units, evidence and reading rules. The plot generator writes an ignored
four-panel PNG and hash manifest. [sources.json](sources.json) records inspected
sources; [checks.json](checks.json) records tests and known environment blocks.

## Local dusk columns

The original study through commit `9ad3a7e` ran in the remote workspace, using committed cloud summaries and packed sky data. The follow-up here reads the ignored raw CM1 fields from the local research drive. It preserves the earlier screen and adds [local_clouds.json](results/local_clouds.json), schema `terluna.research.cloud-twilight-local/1`. The runner is [local.py](local.py); the domain exporter is [cloud_columns.py](../../../climate/crm/cloud_columns.py).

The exporter selects the corrected `ring_equator` case forced by `A28_dim5_moon`, with 6-km columns and 111 layers. The requested interval is the second complete synodic cycle, days 29.53059–59.06118. Saved three-hour snapshots cover days 29.625–59.0: 236 snapshots and 23,810 columns at local hour angles 90–110°. At equinox on the equator, that sector extends from geometric solar-centre sunset to about 39.4 hours afterward. Each sampled column has equal weight. The results describe this two-dimensional ring experiment.

Raw pressure, potential temperature, vapour, five condensate species and four particle-number fields are retained in `research/runs/optical_comfort/dusk_columns.npz`. Dry density is reconstructed with the coefficients from the run’s hash-checked CM1 source. Its domain means agree with CM1’s saved density diagnostics to within 0.0000122%. The companion source manifest records every consumed snapshot hash. These bulk products live on the research drive.

Cloud liquid plus cloud ice reaches the existing threshold of 0.01 g/kg in 22.53% of sampled columns. Among cloudy columns, cloud-top heights have a median of 23.86 km, a 90th percentile of 35.79 km and a maximum of 79.79 km. Cloud occurrence declines from 28.4% at 0–5° solar depression to 19.4–19.7% at 10–20°.

| Height band | Columns containing cloud in the band |
|---|---:|
| 0–20 km | 13.29% |
| 20–40 km | 13.07% |
| 40–60 km | 0.265% |
| 60–150 km | 0.659% |

Bands can overlap in the same column. These fractions measure cloud occurrence somewhere within a band; the original report’s occurrence at each individual level measures a different quantity. The high-altitude population is intermittent. Its recovered particle numbers and condensate layers now provide inputs for a cloud-optics calculation.

For a prescribed uniform liquid-droplet effective radius, the visible geometric-optics approximation gives `tau = 3 LWP / (2 rho_liquid r_effective)`. The calculation uses Morrison’s recorded liquid density and the actual cloud-liquid mass path. Height-band integration preserves partial layers and closes the total species mass to within 1.5 × 10⁻¹⁴ kg/m². [Method notes](local_sources.json) identify the equation and the inspected CM1 source.

| Prescribed effective radius | All dusk columns with liquid tau ≥ 1 | With liquid tau ≥ 10 |
|---|---:|---:|
| 5 μm | 23.0% | 17.5% |
| 10 μm | 21.9% | 15.2% |
| 20 μm | 20.3% | 12.1% |

These radius values are sensitivity scenarios. The saved cloud-droplet concentration is fixed at 100 cm⁻³ in the simulation; recovering its effective-radius distribution requires the microphysics size-distribution closure. Cloud ice, snow, graupel and rain remain separately reported mass and number fields. Their phase functions and extinction per unit mass need their own treatment. Mass alone gives limited guidance about the relative optical importance of the species.

Selected cloud-liquid and cloud-ice mass maxima have incident molecular-path illumination of roughly 9,100–18,600 lux at their diagnosed tops. Those examples use actual cloud heights and local solar angles with the earlier solved mean molecular column. The illumination uses the original binned ASTM solar source. Cloud extinction along the solar path, self-shadowing, multiple scattering and the surrounding diffuse sky remain to be coupled. Incident top-of-cloud lux describes the available source; a visible cloud’s luminance also depends on these transport processes.

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.cloud_twilight.local
python -m pytest climate/crm/test_cloud_columns.py illumination/cloud_light research/studies/cloud_twilight -q
```

Local input preparation and the full check record are linked from the [optical-comfort study](../optical_comfort/README.md#remote-and-local-data-access). Three-hour sampling and a single two-dimensional cycle support this conditional sample; cloud lifetimes, sky coverage from a fixed observer and scene appearance require further spatial and temporal analysis.
