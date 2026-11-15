# Checks for enclosed-industry development

The focused stock/flow and frozen-product suites passed 18 tests on 13 November 2026:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest provisioning/enclosed_industry/tests -q
```

Three numerical ferry controls reproduce the original pinned chemical solver, including tank-limited loading. The other model checks cover zero-demand trips, discrete annual mean lifetimes, dated two-year orders, idle retirement, gross seed procurement versus surviving stock, coproduct feed, feed-bound energy, recycled-material balance, independent computing loads, initial external route fuel charged once, and deterministic gzip packaging. Frozen-product tests verify producer/input identities, all 81 base cases, the population cap, accounting identities, adverse-case identities and finite numbers. These checks establish implementation properties, not physical validation.

The required affected check command also passed:

```sh
make check
```

Its selected scope comprised the layer check (1,657 files; zero violations), historical provenance and numerical reproduction, the new enclosed-industry tests, the six existing provisioning suites, six climate dependency controls and the check-selector/layer controls. Pytest reported 92 passed and one skipped in 13.45 seconds. The skipped optional case would reread 2.3 GB of GCM output and requires `TERLUNA_RERUN_EXTRACT`; it was not run. Full repository verification was not requested or run.

All 81 base cases and six physical sensitivity cases converged and passed every reported annual equipment and resident-workforce screen. Capital and operating subtotals reconcile with each directly computed physical cash bill to within USD2041 1 per annual row. All four compressed products retain their producer, input and postprocessor hashes, with separate compressed and decompressed hashes in their manifests. Prose was read against those outputs, including the seed scope, the source-water/ferry feedback and the deterministic-life timing comparison.

The results do not establish qualified source access, a mission-by-mission transport schedule, a fully priced bootstrap or installation programme, empirical component lifetimes, or realised commercial demand. Those evidence limits remain beside the results and in each product's reading rules.
