# Molecular light paths to elevated clouds

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
