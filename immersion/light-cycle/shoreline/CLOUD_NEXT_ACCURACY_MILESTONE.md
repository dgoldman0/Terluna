# Next accuracy milestone — A1

> **Status update — 23 September 2026.** A1 through A3 now implement substantial
> portions of the reference-accuracy program described below. `A1_METHODS.md`
> establishes the shared atmospheric state, first-order spectral reference and
> frozen-cloud benchmark. `A2_METHODS.md` adds full-spectrum molecular multiple
> scattering plus cloud-grid/crop convergence. `A3_METHODS.md` adds matched
> gas-plus-cloud transport, records the remaining blue-band molecular convergence
> limitation, and benchmarks the production-style cloud-lighting approximation.
> The production shoreline renderer has not yet adopted these reference results.
> Entraining/mixed-phase parcel physics, conserved condensate/precipitation,
> full-spectrum molecular convergence, coupled cloud-grid/crop convergence, and
> production/cache validation remain open. The original pre-implementation plan
> is retained below as a historical specification; individual claims that every
> task is merely planned should be read in that original context.


## Objective and completion claim

**A1: a density-consistent atmospheric column and quantitatively benchmarked cloud
light transport.** Complete this accuracy work before treating additional visual
polish or a renderer-backend switch as validation of the sky.

The intended result is a scene whose gas, cloud illumination and foreground
scattering share one atmospheric state, with numerical rendering error measured
against a converged fixed-medium reference and cloud-height sensitivity exposed.
Climate prediction and dynamically generated three-dimensional cloud morphology
remain separate future claims. Every task below is planned; publication of this
plan marks none of its gates complete.

[Current limitations](CLOUD_ACCURACY.md) and the
[recorded audit](validation/clouds/accuracy-audit-summary.json) are the starting
point. Freeze the delivered revision-04 hashes before changing a model, numerical
budget or display transform. Separate implementation correctness, physical-model
adequacy, numerical approximation and hardware performance throughout the work.

## A1.0 — Reproducible baseline and error harness

Create a versioned fixture for the audited 160 rays, camera, Sun, sounding, density
seed, noise/profile/light textures and display settings. Restore a runnable
high-sample comparison harness in the repository, preserve raw linear RGB and
transmission readbacks, and verify that its 128/8 variant reproduces production.
The earlier artifact supplies source hashes and raw results; publication of its
summary alone does not make that experiment fully reproducible from the checkout.

Expand the fixtures to all three regimes, both gravitational environments, several
cloud-edge/interior/clear paths, horizon and zenith views, and a Sun-angle sweep
including surface sunset and elevated-cloud illumination after sunset. Include a
matched-surface-pressure Earth/Moon control alongside the existing 1/1.2-atmosphere
comparison. Keep authored morphology identical during renderer comparisons.

Output one machine-readable error report with exact source/data hashes, sampling
settings, ray directions, backend and convergence status. Preserve failed cases.
Brightness, spectral/channel colour, transmission, shadow position and energy
residual require separate reported quantities. Record measurements before display
mapping and repeat end-to-end checks after the documented display transform.

## A1.1 — One gas column for thermodynamics and optics

Define a common column carrying pressure, temperature, vapour, composition,
molecular number density and gravity in SI units. Distinguish moist-air mass
density from dry-air molecular density and from each absorbing species. Ozone,
aerosol and water-vapour prescriptions require explicit profiles and provenance.
Use the same spherical coordinate/datum conventions for cloud positions, observer
height, solar paths and planetary intersections.

Regenerate direct solar transmission, clear-sky radiance, diffuse illumination and
finite eye-to-cloud transport from that common column. Independently integrate
selected vertical, grazing and occulted paths. Run the documented no-ozone case.
Retain the old exponential sky as an explicitly named legacy control so its
pixel-identical regression remains meaningful. Density-coupled skies receive new
reference outputs; their accuracy must not be forced to preserve old colours.

Required evidence: column mass/weight closure, species optical-depth integrals,
pressure-coordinate and height-coordinate consistency, direct-path refinement,
finite-path foreground tests, and revised spectral energy accounting. Resolve or
bound the inherited energy-accounting discrepancy and report twilight separately.
Freeze exposure and white balance while comparing physical transport; document
any subsequent camera recalibration as a distinct change.

## A1.2 — Cloud-height sensitivity with improved parcel physics

Add an entraining parcel/plume option and mixed-phase thermodynamics. Specify
whether condensate is retained, reversibly carried, precipitated or removed
pseudoadiabatically. Apply the chosen latent heating, water budget and condensate
loading consistently to buoyancy; preserve the existing undiluted liquid case as
a named control. Separate the condensation level, free-convection level,
neutral-buoyancy level and any momentum/overshoot estimate.

Sweep initial humidity, parcel heating/launch speed, environmental temperature,
entrainment and condensate treatment across declared ranges. Report changes in
cloud accessibility, base, ceiling and positive/negative buoyancy integrals.
Distinguish numerical convergence from sensitivity to selected inputs. Sensitivity
ranges are scenario results, not observational confidence intervals.

A1 must provide this uncertainty comparison for the approximately 59-km case.
A conserved evolving three-dimensional condensate/precipitation model can follow
separately. Until that model exists, the 5.5% optical-retention prescription,
selected effective particle radii and zero predicted surface rainfall stay
explicit in exports and documentation. State the diagnosed ceiling approximately
in narrative/UI descriptions; retain exact values in machine-readable diagnostics.

## A1.3 — Independent reference cloud transport

Build a slower reference integrator for fixed three-dimensional media with the
same extinction, single-scattering albedo, phase function, solar spectrum,
planetary geometry and surface boundary as the interactive renderer. Start with
analytic or independently calculable limiting cases: vacuum, pure absorption,
homogeneous slabs, optically thin scattering, blocked Sun and Lambertian
boundaries. Confirm reference convergence and uncertainty before using it to
judge the interactive path.

Then solve repeated scattering for thin and thick cloud fixtures, including the
prescribed deep column. An offline Monte Carlo or deterministic transport solver
is acceptable; its numerical uncertainty and domain/boundary assumptions must be
reported. Sufficiently high primary raymarch sampling of the existing empirical
source function serves only as a quadrature reference. It cannot serve as the
independent multiple-scattering reference.

Compare directly lit faces, shaded interiors/undersides, edges, ground/cloud
interaction and finite foreground paths. Quantify the empirical softened-sun term,
regional albedo dependence, equivalent-sphere ice optics and RGB/spectral
reduction separately. Replace or calibrate the interactive diffuse/repeated-
scattering approximation on a training set of fixtures and evaluate it on held-out
fixtures. Preserve a label for optical phenomena excluded by the chosen particle
model, including ice-crystal-habit effects.

## A1.4 — Interactive error budgets and cache fidelity

Select primary and secondary sampling independently using the same-medium
reference. Evaluate adaptive sampling or empty-space skipping by their measured
errors and stability, with fixed-budget controls. Add stratified edge/interior and
horizon sampling so clear rays cannot conceal concentrated cloud-boundary error.

Measure each texture/profile/light lookup and angular sky cache against direct
rays. Vary sky-map resolution independently from scene/reflection quality.
Measure the 32-m observer-cell and one-metre eye-height quantization, two-second
wind updates, ground-shadow-map resolution and altitude assumptions. Test known
cloud ranges/heights, curved-horizon occultation, solar-shadow transitions,
reflections, and observer translations. Report angular or world-space shadow and
edge offsets alongside radiometric error; temporal/cache jumps need their own
regressions rather than still-image approval.

## Provisional acceptance gates

These are proposed engineering tolerances to agree and version before
implementation tuning. They are neither current achievements nor guarantees of
physical truth. A target change requires an explicit rationale and preserved old
results; acceptance must not follow from simply loosening a failed threshold.

| Gate | Required result |
| --- | --- |
| Reference reproducibility | Pinned inputs, runnable commands and repository-local raw outputs; production-equivalent fixture reproduced before budget changes. |
| Common atmosphere | All new sky/cloud/foreground transport reads the same species-resolved column; legacy exponential mode remains separately identified. Hydrostatic column-weight residual below 1e-4 under grid refinement. |
| Direct optics | Independent path quadrature agrees within 0.5% where reference transmission is at least 0.01, and absolute transmission error is below 1e-4 on darker/occulted paths; shadow geometry checked separately. |
| Parcel uncertainty | Independent thermodynamic examples and energy/water accounting pass; entrainment/mixed-phase/input-sensitivity results published for each regime. No universal altitude-error bound inferred from these scenarios. |
| Interactive quadrature | Across each fixed-medium fixture: aggregate normalized absolute luminance difference at most 1%, per-ray p95 at most 2%, maximum at most 5%, and maximum absolute transmission difference at most 0.01. Compare both primary and secondary refinement. |
| Physical transport closure | Converged independent multiple-scattering reference plus held-out fixtures. Initial target: per-ray p95 linear-luminance error at most 5% and maximum at most 10%, with RGB/spectral errors, boundary effects and energy residual also published. Unsupported regimes remain labeled outside validated scope. |
| End-to-end cache | Direct-ray and cached rendered outputs measured on the expanded fixture set; interpolation, shadows, perspective and temporal discontinuities independently reported. Set/version angular and temporal tolerances from these fixtures before approving cache settings. |
| Stability and performance | Zero graphics-context-loss events in accepted runs, including late messages; sustained fixed-route test and explicit backend/adapter record. Consumer hardware and native WebGPU are separate claims requiring their own execution. |

For every relative radiometric gate, report the unfloored metric and an additional
metric using a declared dark-ray denominator floor, initially 0.1% of each
fixture's reference 99th-percentile luminance. Report absolute dark-ray error and
per-regime distributions as well. Keep the historical audit's original metric
unchanged. Reference uncertainty should consume at most one fifth of the relevant
acceptance budget; otherwise classify the comparison as reference-limited.

The historical 128/8 production setting fails the proposed per-ray quadrature and
transmission gates on the recorded deep-convection sample. The high-sample
quadrature test provides no result for the independent physical-closure gate.

## Deliverables, sequence and explicit boundaries

Work in order: A1.0 baseline, A1.1 common atmosphere, A1.2 parcel sensitivity,
A1.3 reference transport, then A1.4 interactive/cache validation. Reference limiting
cases may be developed alongside A1.1. Preserve separate change sets for physical
inputs, solver changes, approximation tuning, and display calibration.

Deliver a shared-column schema and adapters; regenerated optical tables with
input hashes; parcel sensitivity reports; independent transport reference and
limiting-case fixtures; reproducible numerical and end-to-end error reports;
updated model-boundary UI/export labels; and an evidence index linking accepted
and failed runs. Mark each gate pass, fail, reference-limited or untested.

More detailed terrain/vegetation, artistic cloud embellishment, and a TSL/WebGPU
cloud port are deferred from this accuracy milestone. Three-dimensional plume
morphology, evolving cloud microphysics/virga, precipitation coupling and a monthly
climate/circulation solution remain subsequent physical-model extensions.

A1 completion would establish bounded fidelity to specified, tested atmospheric
and optical scenarios. It would still leave the frequency, persistence and full
appearance of real terraformed-Moon weather as open modeling questions.
