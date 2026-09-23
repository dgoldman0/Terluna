# Cloud accuracy assessment — revision 04

Assessment date: 2026-09-22. Audited viewer SHA-256:
`c7ead1fd83e0d55518d2951e08cd46d31f323280885125382a8013b5d74636d9`.
The atmospheric studies were prepared against revision-03 commit
`032765c51b1aeb140c732e8f484e47efda15cf4d`.

This checkpoint visualizes prescribed three-dimensional cloud structures inside
an altitude envelope calculated from an idealized low-gravity sounding, with
approximate atmospheric and cloud light transport. The numerical atmospheric
column has stronger support than the specific silhouettes, brightness and colour.
An overall physical error bound for the image remains unavailable.

[Methods](../../atmosphere/column/METHODS.md), [checkpoint validation](../../immersion/docs/history/05-cloud-columns/CLOUD_VALIDATION.json),
[accuracy-audit summary](../../immersion/docs/history/05-cloud-columns/validation/clouds/accuracy-audit-summary.json), and
[next accuracy milestone](CLOUD_NEXT_ACCURACY_MILESTONE.md) distinguish the model,
its implementation and its remaining validation requirements. The milestone is
planned work. This publication preserves the reviewed rendering implementation.

## Supported calculations and conditional cloud heights

`src/weather-column.js` integrates spherical hydrostatic balance using virtual
temperature, unsaturated parcel ascent to the lifting-condensation level, and a
liquid pseudoadiabat above condensation. Existing tests compare saturation,
condensation and moist-ascent calculations against independent documented MetPy
examples and check numerical refinement. Passing these tests supports those
idealized calculations; it supplies no general validation of lunar weather.
Equation references and their implementation boundaries are in `CLOUD_METHODS.md`.

For the selected lunar convective sounding, the diagnosed base is approximately
3.98 km and the first neutral-buoyancy ceiling is approximately 58.7 km. At that
ceiling the sampled environment has pressure 36,004.9 Pa, temperature 238.34 K,
and density 0.526 kg/m³. Cloud formation there is consistent with substantial
atmosphere in the supplied low-gravity profile. These are outputs of this
specified experiment, rather than measurements of an established lunar climate.

The temperature, humidity, wind, initial parcel heating and launch speed are
prescribed inputs. Circulation, radiation and the monthly illumination cycle have
yet to demonstrate that they produce or maintain those inputs. The Earth
comparison reuses pressure-coordinate soundings and also changes surface pressure
from 1.2 to 1 atmosphere. It therefore combines gravity and pressure effects.
Printed hundredths of a kilometre record numerical precision; physical uncertainty
in the ceiling remains unquantified. The approximately 59-km value is conditional.

## Limitation register

| ID | Current implementation and consequence | Accuracy work |
| --- | --- | --- |
| ATM-01 | Thermodynamic density follows the sounding; molecular optical density retains the separate 48.4-km lunar exponential profile. Clouds and light therefore describe different vertical air columns. | Generate all optical paths and spectral background tables from one composition-resolved hydrostatic column. |
| PAR-01 | The parcel is undiluted and follows a liquid pseudoadiabat. Entrainment, condensate loading, ice latent-heat feedback, momentum and overshoot remain absent. The neutral crossing is an idealized ceiling. | Compare entraining and mixed-phase parcel variants, loading and launch assumptions; publish sensitivity ranges with explicit thermodynamic closures. |
| GEO-01 | The renderer explicitly constructs a body, middle lobe, head and anvil. Horizontal scale, aspect and occupancy are selected; scale is tied to layer thickness. Buoyancy constrains vertical support without solving the shape. | Keep morphology fixed during optical validation, label it as prescribed, then test entraining plume/cloud-resolving structures separately. |
| WAT-01 | A selected retained fraction converts cumulative condensation into optical cloud water. Cloud-water conservation, evaporation, fallout and replenishment are unresolved; new studies provide zero surface rain. | Introduce an accounted condensate budget before claiming predicted opacity or coupling cloud precipitation to surface water. |
| RAD-01 | Beer attenuation and spherical occultation are implemented, while diffuse and repeated scattering use empirical source terms and a selected regional albedo. Brightness and underside colours have unmeasured closure errors. | Benchmark identical prescribed media against an independent multiple-scattering reference with energy accounting. |
| RAD-02 | Foreground in-scattering blends toward the inherited clear sky. RGB extinction, an effective ozone prescription, and differing eye/Sun path treatments remain approximations. | Solve the finite eye-to-cloud path and altitude-dependent incident field from the common spectral medium; quantify each reduced-band approximation. |
| NUM-01 | The normal 128/8 integration budget has measurable raywise error against higher sampling of the same model. | Reproduce and extend the fixed-ray comparison, tune primary and secondary sampling separately, and enforce numerical error gates. |
| CACHE-01 | Angular maps, finite profile/light/noise textures, coarse ground-shadow maps, and quantized observer/time updates approximate the continuous scene. The fixed-ray audit excludes their error. | Measure cache resolution, interpolation, parallax, shadow displacement and temporal transitions independently, then end to end. |
| VAL-01 | Current evidence consists of equation examples, conservation/refinement checks, GPU field readback, and selected software-renderer scenarios. Higher-resolution software tests lose their contexts. Native cloud WebGPU, consumer-GPU performance and long traversal remain unmeasured. | Maintain correctness, numerical convergence, physical validation and performance as separate gates; count every context loss as failure. |

### Morphology and condensate

The convection shader's unequal overlapping lobes and flattened head are
physically motivated design choices. Their dimensions and anvil form have yet to
be obtained from a fluid calculation. Individual towers can end below the
column's common diagnosed ceiling. A familiar angular silhouette cannot establish
physical size: proportionately enlarging a cloud and increasing its distance
preserves angular size. Perspective checks should use known distances, heights,
curvature and occlusion, independently of visual familiarity.

The selected convective retention factor is 0.055. The unmasked vertical column
has liquid water path 8.3475 kg/m², ice water path 2.2006 kg/m² and optical depth
974.369. Horizontal morphology changes the occupied paths. These are diagnostic
optical quantities, and their large values are neither a validated prediction nor
a scene-wide average. Optical liquid/ice partitioning currently feeds no ice
latent heat or condensate load back into parcel buoyancy. Equivalent-sphere ice
and selected effective radii omit crystal-habit optical effects.

### Light transport and colour

`src/cloud-renderer.js` combines a prescribed diffuse sky field, surface-reflected
illumination with regional Lambertian albedo 0.18, a directly attenuated sunlight
term, and an empirical softened term proportional to
`0.09 * exp(-0.16 * sunDepth)`. These terms have yet to be calibrated against a
reference transport solution. Repeated scattering is particularly consequential
for interpreting the optically thick convective column. Raising the sample count
only converges the equations already used by the shader.

Foreground blending of each contribution toward the old clear sky approximates
molecular in-scattering. Blue-grey cloud faces and reduced contrast therefore
contain both intended intervening-atmosphere effects and unquantified model error.
The inherited spectral background also retains its documented energy-accounting
and deep-twilight uncertainties. Default daylight white balance and tone mapping
are display transformations; exposure tuning cannot validate the radiance field.

## Recorded numerical rendering audit

The saved audit evaluated 160 fixed perspective directions (16×10) through the
production cloud-transfer GLSL for lunar convection at 45° Sun elevation. The
camera was at [0, 3.045873320927338, 9] metres, yaw 0, pitch 0.46 radians, vertical
field of view 58°, aspect ratio 1.6. It used WebGL2 on ANGLE SwiftShader. Density,
lighting closures and rays were held fixed while primary-view and secondary-Sun
sample counts changed. The generalized shader matched the unmodified production
shader at 128/8 with maximum component difference zero.

The numerical reference used 4,096 primary and 512 secondary samples. The table
reports differences from that reference in linear luminance and cloud transmission.

| Primary / secondary samples | Aggregate normalized absolute luminance difference | Per-ray p95 | Per-ray maximum | Maximum absolute transmission difference |
| --- | ---: | ---: | ---: | ---: |
| 128 / 8 (production) | 1.117% | 6.361% | 12.349% | 0.127458 |
| 512 / 8 | 0.703% | 3.695% | 11.936% | 0.012718 |
| 512 / 32 | 0.131% | 0.690% | 2.534% | 0.012718 |
| 1,024 / 128 | 0.0315% | 0.1451% | 0.4317% | 0.001475 |
| 2,048 / 256 | 0.00962% | 0.04933% | 0.13617% | 0.001394 |

Aggregate error is `sum(abs(Y - Yref)) / sum(Yref)`. Per-ray figures describe the
sampled directions, not whole-image distributions or statistical confidence
intervals. Transmission error is absolute: 0.127458 is approximately 12.75
percentage points. The 2,048/256 row measures disagreement between the two highest
budgets; it supports convergence for these directions without establishing exact
truth. Some intermediate per-ray errors are nonmonotonic. Primary-only sampling
increases leave secondary illumination error, so a single larger loop count is
insufficient as an accuracy argument.

This audit excludes angular-cache interpolation, scene geometry, display mapping,
physical morphology validity, climate plausibility and the true multiple-scattering
solution. Its 1.117% aggregate figure is not a percentage of physical realism.
Neither 137 passing regression tests nor matching CPU/GPU field values supply an
overall physical error bound.

## Evidence provenance

`validation/clouds/accuracy-audit-summary.json` retains the six model samples,
temperature sensitivity experiment, all recorded GPU comparison metrics, camera,
backend and audited source hashes. It is derived from the conversation artifact
`Open_Moon_Clouds_Accuracy_Audit.json`; its source-artifact record identifies the
original byte count and SHA-256. Per-ray raw readbacks are retained in that
original artifact and are explicitly omitted from the repository summary. The
prior GPU experiment was not rerun while writing this documentation. A runnable,
pinned integration-comparison harness with repository-local raw fixtures is an
initial deliverable of the next milestone.

`CLOUD_VALIDATION.json` and `CLOUD_SOURCE_MANIFEST.json` remain historical records
of the delivered local checkpoint; their original publication fields describe
that preparation event. `PUBLICATION_04.json` separately records the files and
fresh build/regression checks used for publication. The source manifest's original
README hash refers to the pre-publication README; the current publication manifest
covers its added accuracy links. No historical test record is silently rewritten.
