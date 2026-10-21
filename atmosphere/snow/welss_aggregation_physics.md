# Earth aggregation component equations

The [aggregation adapter](welss_aggregation_physics.py) reproduces the selected McSnow component equations needed to extend the published Earth comparison. It adds the distinction between monomers and aggregates during vapour growth, and combines hydrodynamic interception with the author's sticking law and collision area. It supplies no new population prediction and changes neither the frozen lunar column nor the independent Earth growth reference.

The selected conditions are pure ice, diagnostic geometry (`habit0`, `agggeo1`, `snowgeo2`), aspect ratio one, temperature-only sticking option 1 and zero outer random-speed perturbation. Under these conditions, diameter and projected area depend on mass alone. Monomer count matters through the condition `mm < 1.01`: the source's porous monomers use capacitance `D/π`, while aggregates use `D/4`. A nonpositive calculated void volume selects `D/2` in both cases. The adapter preserves that rounding-sensitive source decision; it does not adopt the separate compact-sphere sensitivity. The same mass–area geometry gives the same terminal velocity for monomers and aggregates in this particular aspect-ratio-one case.

The collision kernel has two distinct speed treatments. Rough-ice interception uses the source's internal width `SIG = 0.25`; this is a parameter in its encounter calculation, not measured cloud turbulence. The outer rate uses the absolute difference between deterministic terminal velocities. For aspect ratio one, the author's `ckt0` radii are `D/2`, so their combined area is `π(D₁ + D₂)²/4`. Projected-area ratios enter drag and interception, but do not replace this outer area. Raw interception ratios above one remain intact. Equal-speed particles have a zero outer collision rate.

The temperature sticking law is the source's piecewise option 1, including its interval boundaries. This differs from the exponential sticking law in the lunar control. Matching a published Earth experiment therefore requires retaining its complete set of declared equations rather than changing gravity alone.

## Independent component checks

The [fixture builder](references/build_welss_aggregation_reference.py) compiles the retained author geometry, growth, drag, interception and sticking routines without importing the Python physics adapter. Its driver composes their aspect-ratio-one collision kernel. The [182-case fixture](references/welss_aggregation_reference.json) covers monomer/aggregate pairs, compact and porous masses, identical and nearly equal masses, and both sides of the sticking-temperature boundaries. This is source-component parity, distinct from running the author's population algorithm or comparing with observations.

Diameter, capacitance and sticking agree exactly in these cases. Maximum relative differences are 6.46 × 10⁻¹⁶ for deposition, 2.17 × 10⁻¹⁴ for speed and 8.47 × 10⁻¹⁵ for interception. The composed kernel's largest relative difference is 4.89 × 10⁻⁷ for almost equal terminal speeds, where the kernel is only 2.295 × 10⁻¹⁶ m³/s and the absolute difference is about 1.12 × 10⁻²² m³/s. The test bounds this subtraction error against the sum of the two speeds as well as checking ordinary relative agreement. It also checks reciprocal kernels, zero equal-speed encounters, monomer parity with the frozen growth reference, and unchanged diagnostic geometry when the monomer flag changes.

```sh
PYTHONDONTWRITEBYTECODE=1 python atmosphere/snow/references/build_welss_aggregation_reference.py \
  --output /tmp/welss_aggregation_reference.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python -m pytest -p no:cacheprovider \
  atmosphere/tests/test_snow_welss_aggregation_physics.py -q
```

Five tests passed. The fixture records compiler settings and hashes for its generation code and retained source slices. [The source record](welss_aggregation_physics_sources.json) identifies the pinned public revision, exact sticking routine and BSD license; the existing [Böhm source record](boehm_sources.json) and [Earth source record](references/welss_reference_sources.json) cover the reused routines. The selected equations come from [McSnow revision 1e428d7e901cbc5b91f2040e575f736ec5963a74](https://gitlab.dkrz.de/mcsnow/mcsnow/-/tree/1e428d7e901cbc5b91f2040e575f736ec5963a74). That revision has not been identified as the executable behind the older archived population profiles. A coupled independent aggregation column remains a separate deliverable.
