# Numerical audit of the snow population

The inherited heavy-snow calculation fails as a quantitative snowfall prediction. Its corrected-dendrite, 3 mm/h, 100 crystals/L, seed-0 case produces a maximum dimension of **2.57893 × 10⁹ m**, also its arriving-mass 99th percentile. The component above 10 m carries 3.0439% of arriving mass and 43.365% of the inferred projected-area concentration at the freezing level. Its relative mass-ledger residual is only −7.69 × 10⁻¹³: conservation of weighted mass does not establish physically meaningful particles or optics. The ordinary centimetre/decimetre range also needs numerical assessment; deleting the extreme bins cannot validate it.

The [audit product](results/numerical_audit.json) separates arithmetic counterexamples, six short area controls, and one instrumented replay of the full failing case. [The runner](numerical_audit.py) preserves the original simulator. A private copy of the column function's global bindings supplies a collision wrapper; a copied random generator diagnoses its draws without advancing the actual generator. The full replay's summary, complete mass ledger and maximum arriving dimension exactly equal the existing case output.

## Representation-area controls

The column's implicit horizontal area is 1 m². A representation of area A uses A times the injected particle flux and condensation supply, and A times each collision volume. The audit divides output weights, both observers and every mass-ledger term by A. Particle masses, temperatures, positions, timestep and random seed remain fixed. Apart from count-dependent operations and floating arithmetic, the inherited equations are homogeneous under this transformation.

Two operations violate that homogeneity:

- Exhausted collision pairs use `floor(xi/2)`, although injection and fragments have real-valued weights. An actual collision of equal weights 1.5 produces normalized weights `(0, 1.5)` at A=1, but `(0.5, 1)` at A=2. Both conserve mass; they leave different numbers of packets at different positions.
- A fragment packet enters only if its absolute count is at least one. A production of 0.5 fragments/m² is discarded at A=1 and retained at A=2. This absolute threshold depends on the representation area.

Shima et al. define positive **integer** multiplicities and use the floor split to preserve that assumption. Their formulation includes coalescence efficiency in collision probability and gives `p ≪ 1` as the rigorous timestep criterion; multiple-collision handling extends the scheme beyond that premise. This audit read the relevant sections, not the whole paper. [Shima et al. (2009), §§4.1.1, 4.1.4 and 5.1.3](https://doi.org/10.1002/qj.441), [peer-reviewed author text](https://arxiv.org/html/physics/0701103v2).

The short controls prescribe a liquid-saturated lunar column from −2°C to freezing, depth 2.171 km, supply 3 mm/h and 10 crystals/L. They use four layers, four injected packets per 30 s step, seed 0, and one initial-crystal crossing time each for filling and gathering: 12.42 model hours total. Geometry is `bf95_dmax_dendrites_area`. Every control has nonempty arrivals; they are millimetre snow populations, not a shortened forecast of the deep cloud.

| Area comparison | Breakup | Freezing-level number change | Mass change | Projected-area change |
|---|---|---:|---:|---:|
| 1 → 2 m² | Off or on | 0 | 0 | 0 |
| 1 → 256 m² | Off | −0.0000253% | −0.00000801% | −0.00000589% |
| 1 → 256 m² | On | −0.490% | −0.551% | −0.558% |

With breakup on, the maximum difference between binned area-weighted cumulative distributions is 0.00822; the corresponding whole-cloud mean area concentration changes −0.260%. The largest collision probability is below 0.050 and every accepted batch has at most one encounter. These paired-seed comparisons establish representation sensitivity in this small control. They do not estimate its ensemble bias or predict sensitivity of the deep case. Once branches differ, subsequent random pairings also diverge. Full number, mass and area distributions are retained separately in the product.

## Full heavy-case replay

The replay uses the production case's −15°C cloud top, 3 mm/h supply, 100 crystals/L, 30 layers, eight injections per step and seed 0. It runs 132.72 model hours and took 296.2 wall seconds on one numerical worker. Times below are ages of the prescribed cloud. Events mark the first positive-weight packet crossing each threshold; successive rows can describe different packets and are not a complete particle genealogy.

| First threshold | Model hour | Packet weight, m⁻² | Accepted encounters in that update |
|---|---:|---:|---:|
| 15 cm | 6.067 | 13.378 | 13 |
| 60 cm | 7.325 | 0.49691 | None recorded for this packet |
| 10 m | 12.375 | 0.023542 | 28 |
| 1 km | 38.125 | 5.1050 × 10⁻⁷ | 2,256,664 |
| 10⁶ m | 65.092 | 5.4012 × 10⁻¹² | 4 |
| 10⁹ m | 65.858 | 2.2204 × 10⁻¹⁶ | 1 |

The first 15 cm crossing already involves centimetre-scale collisional growth: a roughly 14.4 cm collector receives thirteen roughly 1.64 cm donors in one accepted sticking batch, with `p=12.4128`. The simulator decides sticking once for the whole batch. Conditional on k encounters, this gives sticking-count variance `k²E(1−E)`, compared with `kE(1−E)` for independent encounter decisions. The early event therefore exposes the decimetre population to the batching approximation. It does **not** quantify the direction or magnitude of the resulting distribution error. At the first 60 cm crossing, no accepted collision was recorded for the packet in that layer call; crossing during the preceding deposition update is consistent with this record.

Near-roundoff residuals appear substantially later. The first flagged positive subtraction remainder occurs at 65.108 hours, after kilometre growth and just after the first million-metre dimension. The diagnostic finds 684 positive remainders no larger than 32 machine epsilons times their donor weights. At 65.858 hours, two packets each weighted 2.220446049 × 10⁻¹⁶ m⁻² merge into a raw 1.10 × 10⁹ m particle, carrying **0.403 kg/m²** of represented mass. Tiny number weight does not imply negligible mass or optical contribution.

Across 207,504,491 sampled superparticle pairs, 2,809,205 have `p>1`; 1,687,613 use multiple encounters and 87,820 have their proposed count capped. Maximum `p` is 2.7592 × 10²⁰ and maximum accepted batch size is 3.9210 × 10¹⁸. These floating-point counters exceed the range in which every integer is representable. That identifies an extreme arithmetic regime, not a separately measured error magnitude. Counts weight sampled pairs equally; they are not real-particle encounter probabilities.

The evidence separates three issues. Absolute count rules demonstrably depend on an arbitrary representation area. Large all-or-none batches already operate during centimetre/decimetre growth. The later tiny residual weights support catastrophic per-particle sizes while retaining material weighted mass. The early ten-metre and kilometre particles precede the flagged near-roundoff residuals, so rounding alone cannot explain their onset. Unlimited extrapolation of aggregate geometry and incomplete large-aggregate breakup remain independent physical-model deficiencies. This pass does not assign every part of the tail to one cause or introduce a corrected solver. A coherent multiplicity formulation, collision-timestep controls and suitable aggregate physics must precede quantitative use of these populations; clipping cannot undo their prior competition for water and collisions.

## Identity and reproduction

The product records full producer-file hashes, climate-column input hashes, shared constants by name/value, and exact hashes of the existing population outputs read for context. The unchanged replay is bound specifically to `research/runs/snow_population/moon_3_100_bf95_dmax_dendrites_area_s0.json`, SHA-256 `0f61a699668bd50592fa446ff204bc45252e96a19a63ff26758a57fe97cf355a`. This ignored case cache may be absent from another checkout. No shared `partial.json` was used. The model/data products remain exploratory and empirically unvalidated.

From the `domain/atmosphere` worktree:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m atmosphere.snow.numerical_audit
python -m pytest -q atmosphere/tests/test_snow_numerical_audit.py
```

`--skip-heavy` runs only the short controls. The focused checks cover deterministic arithmetic counterexamples, collision-free area normalization of outputs and observers, unchanged random sequencing, product identities and extensive moments, and exact heavy-case replay. No aggregation, population or geometry implementation was modified for this audit.

**Provenance, 25 October 2026.** [results/numerical_audit.json](results/numerical_audit.json) was built with the rain
solver as it stood before 25 October, whose Chapman–Enskog vapour diffusivity stood 17.6% below Massman's (1998)
measured value. It keeps its values; its record names `atmosphere/rain/drops.py` at that solver's hash, kept byte for
byte as [atmosphere/rain/baseline/drops_2026-10-24.py](../rain/baseline/drops_2026-10-24.py) and listed in
`research/provenance.json`, and its identity checks read that copy ([shared/baselines.py](../../shared/baselines.py)).
