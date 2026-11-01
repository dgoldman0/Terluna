# Final export of growth-repair controls

`growth_products.py` writes the public `results/growth_repair.json` after the
frozen simulations finish. It adds evidence qualifications and a reporting
diagnostic without changing the simulations, their source hashes, distributions,
quantiles or original threshold summaries.

Run the desired cases with the existing runner, then use the final exporter:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.growth_repair --case earth_limited_m144
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.growth_repair --case moon_limited_m96
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.growth_repair --case moon_limited_m144
OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.growth_products
```

Completed, source-matching caches can be reused; the export does not rerun a
column. Its first stage calls the frozen `growth_repair.assemble()` admission
checks with only the intermediate output destination redirected to a temporary
file. Reference qualifications are recovered from the same verified raw caches.
The final public file is replaced atomically only after qualification and
identity checks pass. A failed export leaves the previous public file intact.
Running the older `growth_repair --assemble` command writes an intermediate
product without these added qualifications; finish with `growth_products` when
publishing the current comparison.

Every case and reference carries `optical_prediction_eligible: false`, including
the generic Earth controls. Their finite-domain checks and numerical transfer
tests do not establish calibrated storm optics. The reference records retain
their original `reading_rule`, `evidence` and `numerical_status`. Any inherited
extinction or 5% beam-distance fields remain conditional `Qext=2`, homogeneous
direct-beam diagnostics; they do not predict human visibility.

The nominal 10 m endpoint reconstructs as `9.999999999999998 m` under the
baseline mass–size inverse. Literal `diameter >= 10` therefore excludes its
mass, even though that pivot represents the configured endpoint. The original
threshold fields remain unchanged. A separate `nominal_endpoint_diagnostic`
reports the top resolved pivot's arriving mass share and classification at the
nominal endpoint after matching **only the final pivot** within 64 binary64
machine epsilons in relative diameter, with zero absolute tolerance. This is a
roundoff tolerance, not a bin width or physical uncertainty. A genuinely
different upper diameter is not matched. Both shares normalize surviving
freezing-level arrivals during gathering and exclude numerical-domain exits;
null shares mean there were no arrivals. Endpoint occupancy remains a numerical
diagnostic, especially for a case that fails its domain check.

The product's `assembly` record identifies the two stages and hashes the final
exporter, frozen assemblers and exact raw cache bytes. Missing caches are recorded
with null hashes. Those identities are checked before and after assembly;
each cache is also checked when read. The original per-run `producer` records
remain intact. Cheap tests check qualification recovery, unchanged original
summaries, endpoint classification, source/cache mutation rejection and
preservation of the previous public product on failure.
