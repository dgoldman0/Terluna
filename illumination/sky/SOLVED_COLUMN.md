# Sky light on the solved atmospheric columns

The follow-up calculation combines the current atmospheric profiles, the chosen
shield and spherical multiple scattering in one radiation field. Its absolute
sky radiance and surface flux feed the optical-comfort study through
[`results/solved_sky.json`](results/solved_sky.json), schema
`terluna.illumination.solved-spherical-sky/1`.

The design Moon has 1.2 atm of dry air and a 294.9 K surface. Its upper-air
temperature and trace gases use the same interpolation between stored solved
columns as `surface_light/model.py`. Titania-stack transmission, multiplied by
0.95, filters the incident TSIS-1 spectrum. The Earth control has 1 atm of dry
air, a 288 K surface and unfiltered sunlight. Both use a Lambertian ground with
albedo 0.1. These are horizontally uniform, clear molecular atmospheres.

## Calculation

`solved_optics.py` retains all 100 atmospheric layers. Each layer conserves the
vertical molecular optical depths supplied by the atmosphere domain and treats
its extinction as homogeneous inside its spherical shell. The visible interval
is 360–830 nm:

- Below 500 nm, the photolysis calculation supplies 1-nm molecular absorption and
  Rayleigh scattering, including the solved trace gases and oxygen collision
  bands.
- Above 500 nm, H2O and O2 lines use the atmosphere model's 0.05 cm⁻¹ grid, with
  the water continuum and collision-induced absorption. The temperature-dependent
  photolysis tables supply the solved trace gases across that grid.
- Incident energy and CIE XYZ weights follow the admitted TSIS-1 HSRS table and
  the shield spectrum. The inherited analytic approximation to the CIE 1931
  functions supplies colour and photopic illuminance.

Transport uses adaptively reduced spectral channels. Within each 10-nm band,
wavelengths are ordered by their column extinction. Groups retain solar-weighted
layer scattering and absorption, plus their exact integrated incident energy and
XYZ weights. Groups are divided until a sampled brightness/colour error budget
is met. Independent slant rays test the reduction; broadband Monte Carlo paths
also draw wavelengths directly from the original 159,177-point spectral input.
This second comparison includes the narrow gas bands throughout multiple scattering.

`solved_transport.py` reuses the existing Rayleigh angular-moment operator in
`solver.py`. Extinction is integrated analytically along each homogeneous shell
segment. The diffuse source is interpolated between radius and solar-elevation
nodes. Angular quadrature splits at the local geometric ground horizon at every
altitude. Solar attenuation uses a finer lookup grid and 12 strips across the
visible segment of a uniform solar disk. Its angular radius comes from the shared
solar radius and astronomical unit. Scattering directions use the disk centre.

Each spectral channel accumulates successive scattering orders until its newest
contribution meets the recorded tolerance. Saved fields retain the last order so
formal ray evaluation reconstructs the same retained series. `solved_checks.py
--finish` continues channels that reach the initial order cap. The reported
increment is a convergence diagnostic; the stored comparisons provide the
additional evidence for numerical accuracy.

The formal-transfer geometry and division into direct and successive scattered
light are consistent with the methods described in
[Bruneton's reference implementation](https://ebruneton.github.io/precomputed_atmospheric_scattering/atmosphere/functions.glsl.html).
The follow-up uses the repository's scalar Rayleigh moment formulation. The
existing [method notes](METHODS.md), [sources](SOURCES.md), atmosphere input
manifests and protection product supply the remaining method and input provenance.
Reading of the external implementation focused on ray geometry, transmittance
and scattering-source integration.

## Numerical evidence

The product keeps several checks separate:

| Check | Question it addresses |
|---|---|
| Analytic tests | Beer–Lambert attenuation, shell geometry, finite-disk occultation, angular integration, moment trace and a wavelength-sampled Lambertian reference |
| Halved absorption-line spacing | Sensitivity to the underlying line spectrum |
| Fine-grid slant beams | Error introduced by spectral reduction |
| Coarse and standard full spectra | Sensitivity of photopic flux to spatial and angular resolution |
| Doubled atlas angular sampling | Convergence of the surface-flux integral of the directional radiance |
| Three-channel spatial refinement | Additional resolution sensitivity in radiance and whole-sphere energy |
| Single-channel Monte Carlo | Independently sampled noon, low-Sun and twilight radiance and diffuse flux |
| Fine-spectrum Monte Carlo | Combined spectral and transport error in broadband diffuse illuminance |
| Whole-sphere energy | Direct escape, diffuse escape, absorption in the air and absorption at the ground |

Monte Carlo rows include photon counts and sampling standard errors. Very faint
Earth twilight can yield few positive contributions. A sample with zero contributions
is recorded as unresolved, and the percentage comparison requires a relative
standard error of at most 2%. The deeper Earth twilight cases need importance
sampling. The Moon's much brighter twilight has separate usable comparisons.

The global ledger integrates absorbed energy through every layer and outgoing
light around the entire sphere using the moment-grid quadrature. Surface light
for the study comes from integrating the finer directional atlas. This resolves
the sharper Earth horizon contribution. Both flux estimates and the energy
residual remain in the product alongside local radiance checks. The historical atlas retains its own energy
record; the new calculation has an independent ledger.

## Products and reproduction

Bulk inputs, optical tables, moment fields and sky atlases live under
`research/runs/optical_comfort/spherical`, linked to the research drive. The
compact committed product records their hashes and the producer files. Existing
atlases continue to serve their established consumers.

```sh
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=3 python -m illumination.sky.solved_sky --compute
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=3 python -m illumination.sky.solved_checks \
  --grid --line-grid --monte-carlo --finish
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=3 python -m illumination.sky.solved_sky --export
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=3 python -m illumination.sky.solved_checks --angular --summarize
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=3 python -m illumination.sky.solved_sky --export --reuse-atlases
OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.spherical
MPLCONFIGDIR=/tmp/terluna-mpl python visualization/optical-comfort/spherical.py
```

The atlas stores XYZ radiance on 41 upper-hemisphere elevations and 65 symmetric
azimuths, at 30 solar elevations. Y has units of cd m⁻²; flux Y has units of lux.
Interpolation is bilinear in sin(elevation) and azimuth. The exporter compares
the angular interpolant's horizontal integral with the internal moment-grid
flux, and uses the atlas integral for the surface product. Doubling the atlas's
angular resolution supplies a separate check. The scene integrator uses that
same integral to normalize its angular field. Both amplitude and angular pattern
therefore describe the same air and incident spectrum.

Regional temperature and humidity profiles, clouds, aerosols, polarization,
refraction, terrain and measured surface reflectance remain additional inputs.
The present output describes a conditional molecular sky; human visual comfort
also requires the viewing task, adaptation and human-factors evidence.
