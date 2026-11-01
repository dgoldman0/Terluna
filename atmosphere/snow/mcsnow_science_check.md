# Recent-science check of the tuned crystal-growth model

The frozen two-temperature model fails the humidity-dependent morphology check. It reduces mass uptake when vapour availability falls, but continues to prescribe porous, strongly directional growth where the experimental evidence favours compact growth. This confirms a previously recorded limitation in the actual tuned component. The author's condition for proceeding to the proposed freeze-and-aggregation work is therefore unmet. The coefficients, native executable and earlier products remain unchanged.

[Harrington and Pokrifka (2024)](https://doi.org/10.1175/JAS-D-23-0131.1) compile experimental evidence that low vapour excess favours compact crystals, while increasing excess promotes hollowing and branching. Their approximate transition uses relative excess, σ = (e − eᵢ)/(e𝓌 − eᵢ), with a broad range around 0.2–0.4. The range combines different experiments and morphologies; it supplies a qualitative response criterion, not a universal threshold or an exact density target.

The [runner](mcsnow_science_check.py) holds the frozen controls fixed and probes solid spherical initial states at −5.3 and −14.4 °C, 1010 hPa, diameters of 20, 60 and 250 µm, and σ = 0.05, 0.1, 0.2, 0.4 and 1. Each state advances for 1/32 s at two timesteps. Comparing matched initial states isolates the humidity response from changes accumulated through a longer growth history. The [product](results/mcsnow_science_check.json) records all 60 probes, exact inputs, source and executable identities, and raw-output hashes.

For the 60 µm warm crystal, the refined results are:

| Relative excess σ | Ice supersaturation | Mass uptake, kg/s | Added mass / added enclosing volume, kg/m³ | Polar/equatorial logarithmic growth ratio |
|---:|---:|---:|---:|---:|
| 0.05 | 0.262% | 3.5731 × 10⁻¹⁴ | 247.980 | 4.48615 |
| 1.00 | 5.237% | 7.1470 × 10⁻¹³ | 247.983 | 4.48613 |

Reducing excess vapour twentyfold reduces uptake about twentyfold, while the density of the added growth and its axis allocation remain essentially unchanged. The added-volume density is only 27% of the native solid-ice density, 919 kg/m³. It is a model diagnostic of sparse growth, not a measured density or the bulk density of the whole crystal. This check therefore identifies a missing physical response without assigning a new percentage error against an unmatched laboratory trajectory.

The native source explains the result. Its habit factor comes from temperature; warm columns receive hollow-growth density immediately, modified by the tuned size-dependent factors. The fitted compact-stage override applies only below 51.1 µm. At a given crystal state these rules contain no explicit humidity dependence. Humidity can alter subsequent size and ventilation indirectly, but cannot supply the missing compact-growth response at the tested state. The cold 250 µm probes likewise retain approximately 357.46 kg/m³ across the humidity sweep. The source inspected is the hashed `column_density/build/source/src/mo_depo.f90` in the retained run cache.

The other recent sources qualify the scope of this conclusion:

| Source | Consequence for this audit |
|---|---|
| [Takahashi (2014)](https://doi.org/10.1175/JAS-D-14-0043.1) | At maintained water saturation, mass changes with cloud liquid-water content while diameter changes much less. The deposition-only driver sets liquid water to zero, and the inspected growth and drag components do not use that quantity. This experimental dependence is omitted; a numerical error requires a matched reconstruction. The older fitted references remain intact. |
| [Libbrecht (2023)](https://arxiv.org/abs/2306.13087) | Controlled needle-tip observations constrain morphology under varied temperature and humidity. Different initial geometry and undigitized image panels prevent a pointwise validation here. |
| [Pokrifka et al. (2025)](https://doi.org/10.1175/JAS-D-25-0030.1) | The −67 to −46 °C substrate experiments lie outside the two calibrated temperatures. Extrapolating those controls would not provide a valid direct test. |
| [Fuchs et al. (2025)](https://acp.copernicus.org/articles/25/12177/2025/) | Field growth depends on particle competition and moisture supply. Comparing those experiments requires population and atmospheric feedback, beyond these prescribed-air probes. |
| [Harrington and Pokrifka (2026)](https://doi.org/10.1175/JAS-D-26-0016.1) | Rim width and step-source location remain unresolved. Their omission alone does not falsify an effective bulk fit. This source was already recorded in the hollow-growth investigation; its warm comparison reuses the 1991 observations. |
| [Lamb et al. (2026)](https://doi.org/10.1126/sciadv.aee6174) | Learning mass-dependent growth from measured histories does not itself contradict a model that already depends on size and volume. The study does not jointly validate the warm-column dimensions and fall speed. |

The [source register](mcsnow_science_check_sources.json) distinguishes published text, indexed excerpts and preprints. Complete final texts for the two 2026 papers remain unavailable in this review; no missing equations or experimental tables were substituted.

Numerical verification reproduces all 24 sampled frozen-baseline states exactly. Halving the probe timestep changes the derived rates by at most 0.0567%; the maximum mass-ledger residual is zero and the maximum relative volume residual is 2.11 × 10⁻¹⁵. Raw inputs and outputs occupy 38,192 bytes in the existing project-drive cache. The audit runs with one worker:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONDONTWRITEBYTECODE=1 python -B -m atmosphere.snow.mcsnow_science_check
```

Repository verification: `make check` passed with 3,801 Python tests passed and 20 skipped, 191 JavaScript tests passed, 1,338 files checked without layer violations, and passing provenance and ensemble checks. The [complete log](../../research/runs/snow_scale/mcsnow_reproduction/growth_diagnosis/science_check/repository_check.log) is retained with SHA-256 `120d010abb946172b29a0172da53b2102afc4e3f26e414fd2aa5da0fa3ff295e`.

The existing [10% and non-regression requirements](mcsnow_joint_acceptance.md) remain in force, including the unresolved fall-speed errors. This humidity check does not explain those saturated-air errors. A physical humidity response would need independent constraints before reduced-moisture transfer could be accepted. This bounded audit introduces no replacement law or further parameter search; aggregation continuation remains withheld under the author's stated condition.
