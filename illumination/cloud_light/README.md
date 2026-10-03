# Cloud light through a lunar evening

The current calculation follows spectral light through the solved molecular
atmosphere and a saved CM1 cloud section. It produces absolute cloud radiance,
colour, contrast against a clear sightline, open-ground illuminance and the
light admitted by a specified viewing recess. The coupled
[cloud-twilight study](../../research/studies/cloud_twilight/) supplies regional
selection and evolving weather histories.

`microphysics.py` derives effective radii from the producing Morrison scheme's
mass and number distributions. Liquid uses its gamma distribution; cloud ice,
rain, snow and graupel use the scheme's exponential distributions and slope
bounds. The visible geometric-optics area is `Q_ext = 2`, with equivalent
spheres at the recorded particle densities. Particle absorption is set to zero.
Liquid has Henyey–Greenstein asymmetry `g = 0.85`; frozen particles have `g = 0.8`,
with `0.7–0.9` sensitivity scenes. These are optical scenarios. Ice habit,
polarization, halos and rainbows need more detailed particle optics.

`volume.py` samples backward photon histories with exact piecewise-constant
optical-depth segments. Molecular and particle scattering, repeated scattering
within clouds, solar shadows, intervening atmospheric light and a Lambertian
ground all contribute to the same transport. Each scattering event samples the
finite solar disk. Spectral source weights come from the admitted shielded
solved-column product. Independent photon blocks retain sampling errors.
Solar-path integration stops beyond optical depth 40, where the remaining
transmission is below 5 × 10⁻¹⁸. Segment and event limits cause a failed render if reached; a saved diagnostic
identifies that state.

Cloud cells are bounded by concentric shells and along-ring angular planes.
The two-dimensional section is extruded across a stated width, bounded by two
parallel planes. Widths of 100, 200 and 400 km test this missing spatial
information. The ground is spherical, with albedo 0.1. A raised observer sees
through the same atmosphere and cloud field; this is a viewpoint-height
sensitivity. Local mountain surfaces need their own geometry.

`render.py` writes XYZ radiance (Y in cd/m²), clear reference radiance, block
estimates, statistical errors and scene geometry to the research drive.
Surface illuminance uses an independent cosine-weighted hemispheric estimate,
including the attenuated finite direct Sun. The viewing-recess probe admits a
20° wide opening at elevations 5–25°; every other surrounding direction has
zero radiance. Reflecting shelter walls would increase the local illumination.
The recess is an explicit viewing geometry, suitable for exploring a dark
foreground with an open view of bright clouds.

`evening_context.py` extends the clear solved sky through Sun elevations of
−60°, compares standard and draft numerical grids, and reports straight-ray
geometric illumination bounds. Clouds, terrain and Earthlight contribute
separate scene-dependent changes. `detail.py` refines the low-horizon angular field and
its 0.5–5° viewing opening. `earthlight.py` supplies the separate mean-orbit,
solar-coloured Earthlight proxy. Faint twilight is also checked against the
independent clear photon histories accompanying the cloud scenes.

The tests cover distribution moments and mass/area consistency, HG
normalization and inverse sampling, analytic cloud-segment optical depths,
finite-Sun illumination and cloud shadow, and an absorbing atmosphere above a
Lambertian surface. Numerical tests establish implementation behaviour;
regional three-dimensional climate and physical cloud appearance require
further evidence. [Source notes](../../research/studies/cloud_twilight/evening_sources.json)
record the inspected equations and producing CM1 source hash.

## Earlier beam-only screen

`model.py` adds a small spherical beam calculation to the existing optics:
analytic distances through concentric shells, solid-body occultation, finite
Sun averaging, and surface-to-feature geometry. It imports the current
surface-light column builder, Rayleigh cross sections, sky spectral data,
photopic fit and protection transmission. No historical solver is modified.

The current column preserves each layer's vertical optical depth exactly and
uses constant extinction within that layer. The comparison atmosphere is the
sky atlas's zero-ozone exponential proxy, conservatively averaged into 1-km
shells. Layer resolution is a model approximation; analytic shell geometry
does not remove it. The uniform finite Sun uses small-angle strip quadrature.

Only molecular extinction of the direct beam is calculated. Water/oxygen
absorption bands, aerosol and cloud extinction, refraction, in-scattering,
cloud phase functions and cloud radiance remain outside this model. Ozone
uses the existing coarse visible cross sections; the design Moon's solved
ozone is negligible. These outputs are therefore a mechanism screen.

The [cloud-twilight study](../../research/studies/cloud_twilight/) couples these
paths to committed CM1 cloud statistics and a separately labelled twilight
sky proxy. It owns the versioned result and full input/source hashes.

Tests cover analytic vertical and tangent paths, ground blockage, finite
segments, reciprocity, path composition, the finite Sun in vacuum, an
independent continuous exponential integral, and the stored column depth.

Method references: [PBRT, Transmittance](https://pbr-book.org/4ed/Volume_Scattering/Transmittance),
and the existing [sky sources](../sky/SOURCES.md). The shell-intersection
implementation is independent, and numerical checks are not physical validation.

The original full-field renderer’s scalar pixel-error diagnostic uses a 1 cd/m²
floor in its denominator. The coupled appearance product computes the true
relative error from each positive pixel’s stored photon blocks. Interpretation
of faint scenes uses those blocks and the reported angular-patch errors.
