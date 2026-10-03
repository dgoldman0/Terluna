# High clouds through the long twilight

The repository supports a useful first study without additional local data.
Elevated cloud material can receive a strongly reddened direct beam well after
surface sunset. The current evidence favours investigating illuminated sides,
rims and thin upper structures against a still-luminous twilight sky. It does
not yet establish bright cloud towers over a dark landscape or their appearance.

The study remains on `study/optical-comfort`, alongside the preceding optical
work. It reuses the current surface-light column, optical constants, protection
transmission and packed-atlas reader. A small new illumination-domain module
adds spherical paths to elevated features. No climate simulation is rerun.

## What the corrected cloud run supplies

The committed `ring_ring_equator` products use corrected `A28_dim5_moon`
forcing: model days 29.5–59, 237 snapshots, 6-km horizontal columns, and a
two-dimensional equatorial ring with prescribed land moisture and upper forcing.
The atmosphere used for optical paths is the global reference column, not a
simultaneous sounding within those storms.

Cloud-top statistics refer to the highest cloud-liquid/ice cell over each
detected rain footprint in a snapshot. Cloud is `qc + qi >= 1e-5 kg/kg`; a rain
footprint is a contiguous region above 1 mm/h. They are not percentiles of all
clouds or independent cloud events.

| Storm-footprint statistic | Cloud-top height | Overhead solar-centre sunset after surface sunset | Direct illumination falls below 1,000 lux |
|---|---:|---:|---:|
| Median | 22.1 km | 17.9 h | 10.5 h |
| 90th percentile | 37.8 km | 23.3 h | 15.5 h |
| Largest sampled value | 83.8 km | 34.3 h | 27.2 h |

The last column is a molecular-only, design-shield calculation of illumination
normal to the beam. **1,000 lux is an arbitrary reference level**, not a cloud
visibility threshold. Geometric durations assume fixed height and position at
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
incident light. An approximately horizontal upper cloud face is not illuminated
by a Sun below its own horizontal plane; exposed sides or undersides may be.

**The sunsetward view deserves the first detailed scene.** With the observer's
Sun 5 degrees down, a 40-km feature 100 km west is 19.9 degrees above the visible
horizontal. The Sun at that feature is only 1.70 degrees down. Its incident
direct normal illumination is 7,760 lux. The cloud-to-observer ray then removes
more light: the spectrally combined Sun-to-cloud-to-eye source factor is 2.70%
of the above-atmosphere design source, compared with 1.97% overhead. Among the
sampled westward ranges, this factor peaks near 75 km at that time. At later
twilight the favourable range moves farther west. This is a transmission
diagnostic, not an optimum for observed cloud brightness.

The same 100-km example has a 24.9-degree scattering angle, where zero denotes
forward scattering. A resolved cloud calculation must use the appropriate
liquid/ice phase functions and internal optical depth. The model does not assign
a brightness boost from this angle. A 40-km feature 400 km away is hidden by
the solid body from the assumed 1.6-m observer, even though the feature itself
can receive sunlight.

**The direct component becomes very red.** In that example, solar transmission
to the feature is 0.074%, 4.21% and 20.1% at 450, 550 and 650 nm. Including the
viewing path leaves 0.0061%, 1.41% and 11.5%. This spectral filtering supports
reddening of the direct component; perceived cloud colour also depends on
diffuse illumination, cloud optics and atmospheric light along the sightline.

**The surface does not quickly become dark in the existing sky model.** The
unfiltered, zero-ozone proxy retains about 3,220 lux on a horizontal surface
when the Sun is 5 degrees down, 1,650 lux at 10 degrees, and 772 lux at 15 degrees.
Those are a separate atmospheric calculation. The study also computes direct
cloud beams in that same proxy: at 40 km overhead, about 3,540 lux at 5 degrees
and 296 lux at 10 degrees. Even this internally matched comparison is of two
different illumination planes/components, not a cloud-to-sky luminance contrast.
The current shielded column still needs a spherical diffuse-sky solution.

**Very high dusk clouds are a special case in this run.** In the 90–100-degree
local-hour-angle bin, the largest cloud occurrence at any sampled height from
20 to 40 km is 5.54%; from 40 to 60 km it is 0.11%. These are fractions of
sampled columns at individual heights, not the fraction of the sky covered or
the chance that an observer sees a cloud. The 60–80 km scenarios expose the
height dependence but should not be presented as routine dusk weather.

The rain-footprint tracking diagnostic has a median lifetime of 6 h and a
90th percentile of 15 h, at 3-h output spacing; duration is the number of
detected samples times that spacing. It is not an anvil-lifetime measurement.
The geometry therefore supplies no evidence that a single cloud maintains the
same shape throughout an 18–34 h illumination window.

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
from the sky domain. They are not the missing WHI spectrum used by the full
surface-light calculation. The design source applies the stored titania-film
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
filtering. No cloud scattering coefficient or phase law is applied to this factor.

Omitted: water/oxygen and other absorption bands, aerosol and cloud extinction,
cloud self-shadowing and multiple scattering, atmospheric in-scattering along
the viewing ray, refraction, polarization, terrain and a simultaneous local
atmospheric profile. Additional extinction reduces the corresponding ballistic
component; diffuse light and refraction require their own treatment. This is
not a rigorous bound on total visible cloud radiance.

Numerical checks cover vertical/tangent analytic limits, reciprocal finite
paths, path composition, ground blockage, the half-visible finite Sun in
vacuum, and independent continuous integration of the exponential proxy.
For the sampled checks, 32 versus 64 solar strips differ by at most 1.8 lux;
1-km versus 0.5-km proxy shells differ by at most 0.12% where the beam exceeds
1 lux. These checks do not constrain the physical error of coarse profiles,
missing optics, or the two-dimensional climate model.

## Where local data becomes necessary

The next useful step is a **specific dusk scene**, initially with upper
structures around 20–40 km and observers tens to a few hundred kilometres
east of them. The compact products suffice for the present screen; they do
not contain what is needed to determine which cloud surfaces glow, the colour
contrast, or how long the scene lasts.

The needed local input is a short sequence of corrected CM1 dusk snapshots
with geometry/time, density and temperature, separate `qc/qi/qr/qs/qg`, and
particle-number or size information where available. Missing microphysical
sizes and ice habits would still require stated optical assumptions. Begin
with a 2-D diagnostic if only a ring is available; a credible scene needs
3-D structure or a clearly labelled extrusion sensitivity. The ignored
`climate/crm/products/ring_ring_equator.npz` may help select events, but raw
CM1 files are needed for fields that the reduction discarded.

The saved maximum-rain storm slice peaks at local hour angle 21.9 degrees,
well before sunset. Its condensate combines liquid, ice, snow and graupel.
It is neither a dusk snapshot nor enough information to assign extinction.
The atlas's stored cloud illumination is only at 2.5 km; it is not reused at
20–80 km. High-resolution clear-sky atlases alone would improve angular
background sampling but would not supply the missing cloud optics.

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
