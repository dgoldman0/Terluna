# Unified branch merge plan

**Proposed, not executed.** This plan joins the five pinned states reviewed on 8 October 2026. It preserves source histories, author decisions, reproducible historical results and unresolved science while producing one usable branch. The current analysis workspace stays the integration workspace after the author approves execution. No step below has been run merely because it appears here.

## Scope and the review decision

The proposed implementation is repository integration and necessary reconciliation: merges, documentation/status repairs, producer/consumer provenance corrections, bounded deterministic postprocessing, the infrastructure input refresh proposed by its cap for joint integration, relevant builds and checks. It excludes new climate or storm simulations, full ring-fleet searches, stage 3 electricity, array-industry research, new illustrations, manuscript revision, remote pushes and publication. A required missing input is reported; no substitute dataset or optimistic PASS is allowed.

The author can adopt this plan as written, narrow it or change the order before any execution. Three recommendations deserve explicit attention in that review:

- Preserve the magnetic architecture result as a historical comparison tied to its actual main-era design-point input, instead of reopening the capped orbital campaign just to make its freshness test pass.
- Preserve electricity's traceable historical diagnostic producer; add honest dependency/reading metadata without claiming a simulation rerun. Recompute cheap derivatives only when changed dependencies require it.
- Use the existing infrastructure form as the current concept geometry implementing the adopted form, with the square lattice and old fire arrangement retained as named comparisons. Correct mixed-geometry consumers; do not choose a new tower shape or certify the old fire design for it.

These are concrete proposed resolutions. Leaving any unresolved does not invalidate the analysis, but it prevents the corresponding implementation gate from being called complete.

## Source identity and isolation

| Ref | Exact source commit |
|---|---|
| `main` | `2f2e0a104a3a1ea962a6dead4179f522762f997e` |
| `study/sea-appearance` | `46681f458a8c21eb9e4589e784965fcc1acc2b48` |
| `research/solar-shield-habitat-array` | `eeaecbc4c014f99366ae965286c63191da9b3183` |
| `research/atmospheric-electricity-plan` | `24d236b03c760171de8bbd98282c8902cd72b158` |
| `habitation/infrastructure` | `625f60e83313cadd21a7999a772d93dff07cee17` |

Working branch: `research/branch-integration-analysis`. Working directory: `/home/kir/Documents/Projects/terluna-integration-analysis`. Keep this branch as the eventual unified working line; no additional source branches or workspace moves are required. Its name may be changed later only if useful to the author.

Before execution, verify each ref still points to the recorded SHA, source worktrees remain untouched, and only the approved analysis is pending in this checkout. Use `python research/integration/audit.py` to repeat the pinned-source and pairwise checks without touching repository refs/index/worktrees. Ref drift requires a delta review, not silently replacing the source set.

Do not fetch/reset to remote tips. The shield cap is local-only relative to its stored origin ref; merging `origin/research/solar-shield-habitat-array` would omit it. Do not delete or rename old source branches or worktrees. Preserve analysis and original histories as the recovery points.

Restore required input bytes into isolated ignored locations with size/hash verification. A read-only mount or verified copy is appropriate for immutable source data; a normal writable symlink is not read-only. Use private output/cache/build directories. Do not recursively copy whole research drives, run builds against shared CM1 homes, or run producer commands in source worktrees. Archived CM1 cases point at shared executables that a routine rebuild could overwrite. Any required build must set `TERLUNA_CM1_HOME` to a private source/build tree, and writing producers must use private `climate/crm/runs` and output paths. Temporary test harnesses may read verified archived sources; no new full build is required merely for constants cleanup.

## Measured Git conflicts

The review ran Git 2.43.0 `merge-tree --write-tree --name-only` for each pair, with `GIT_OBJECT_DIRECTORY` in a temporary directory and the original object store used only as an alternate. This computes the merge tree for inspection without a checkout, index change, commit or ref mutation. The trees are disposable; [snapshot.json](evidence/snapshot.json) retains messages and names. It is not a cumulative four-branch rehearsal.

Abbreviations below: R = `research/README.md`; D = `research/decisions.md`; S = `research/status.json`; G = `climate/gcm/README.md`; C = `climate/crm/README.md`; H = `habitation/README.md`; A = `atmosphere/README.md`; V = `visualization/README.md`.

| Pair | Paths changed by both since their common ancestor | Text conflicts |
|---|---:|---|
| Main + infrastructure | 6 | G, H, R, S |
| Main + electricity | 4 | C, R, D, S |
| Main + shield | 0 | None; main is ancestor |
| Main + sea | 0 | None; main is ancestor |
| Infrastructure + electricity | 4 | G, R, S |
| Infrastructure + shield | 7 | G, H, R, S |
| Infrastructure + sea | 6 | G, H, R, D, S, V |
| Electricity + shield | 5 | A, C, R, D, S |
| Electricity + sea | 4 | C, R, D, S |
| Shield + sea | 4 | None; overlapping R, D, `shared/constants.json`, `shared/constants.py` auto-combine |

Overlaps differ with the pair's ancestor; comparing `git diff main branch` can misleadingly count inherited missing work as branch-authored edits. Use the snapshot's three-way bases. A clean pairwise result says nothing about all scientific interfaces or the final cumulative tree. Check again after each real integration stage; never apply blanket `ours` or `theirs` to a hub file.

## Staged execution after approval

Before any approved commit, restore a usable check environment and run `make check` as required by the root instructions. If a check fails, retain the exact failure and either resolve it within approved scope or leave the stage uncommitted and report the limitation. A permissive protection-input status is not complete spectral restoration.

First, make the approved analysis checkpoint only if the author has requested that commit. Use an informative multi-phrase message such as “Review five research states; record evidence boundaries and the unified integration plan”. Its check record must describe this main-based tree, not pretend the future combined code has passed.

Then use ordinary, non-squashed merges by **exact SHA**, with `--no-ff --no-commit` to inspect each result before making its merge commit. The sketch below describes individual stages; it is not a script to run all merges through unresolved state.

| Stage | Merge target | Reconciliation before its commit | Exit condition |
|---|---|---|---|
| 1 | Sea `46681f4…` | Preserve main constants-by-value policy, sea's photometric decisions, portable calendar/Make target, historical images and frozen-experience brightness correction; add missing sea/calendar current status | Source/product tests, required Make checks; almanac versus complete sky assets distinguished |
| 2 | Shield `eeaecbc…` | Union constants by key; correct stale shield current summaries; explicitly preserve magnetic comparison against its historical parent; audit changed producer identities | Magnetic test has truthful meaning and passes or stage remains unresolved; affected products/current read rules agree |
| 3 | Electricity `24d236b…` | Retain corrected climate and all optical/wave additions; present the branch’s corrected 530-flash result as current; preserve original/race/window evidence; reconcile constants/contracts and historical diagnostic identity | Targeted checks with real skips; current stage-2 state and pending stage 3 visible; no shared archive executable altered |
| 4 | Infrastructure `625f60e…` | Retain all newer environment/status work; point studies to accepted corrected products; fix mixed geometry inputs; run the five ordered postprocessors; refresh prose and relevant diagrams | Corrected results and provenance agree; conceptual designs stay conditional; all source histories reachable |
| 5 | Final integration reconciliation | Finish topic/decision/plan index, dependencies and acceptance record; inspect aggregate diff and check final outputs | All acceptance criteria below met; clean working tree after approved commits; no unapproved research or publication |

Every stage records checks at the tree actually committed. If repairs are sizeable, use informative reviewable repair commits after the merge once that stage's required checks pass. Do not flatten discarded evidence or overwrite prior source refs. The final outcome contains all source ancestry, not merely a copied union of files.

A normal unresolved merge can be aborted in this dedicated checkout after preserving any needed resolution notes; this leaves source branches untouched. Do not use `reset --hard`, recursive cleans or deleting run storage as routine recovery. If an approved stage is already committed, use its named checkpoint and a reviewed follow-up/revert, not destructive history rewriting. No rollback actions are needed in the current analysis-only state.

## File resolution policy

| File or family | Retain and reconcile | Specific trap |
|---|---|---|
| `research/README.md` | All domain capabilities/studies with current boundaries and cap links | One side's entire index hides other work |
| `research/status.json` | One coherent current record per topic; add sea/calendar, electricity/aerosol, shield/industry status; updated date reflects actual edit | Do not concatenate condition strings into mutually contradictory present claims; don't silently reinterpret the inherited `base_commit` field as this audit snapshot |
| `research/decisions.md` | Author choices from all branches, dates and explicit replacements | Five-hour twilight, full-ring tower and older shield gates cannot remain simultaneously current; cap recommendations are not decisions |
| `research/plan.md` | Replace stale active work order with the approved integrated queue and completed/current/pending distinctions | Existing installation/topography/running-box text predates completed work; keep climate pause and no rented compute |
| `climate/gcm/README.md` | Main corrected climate/solar clock/comfort limits plus infrastructure wind exports | Infrastructure started before the corrections |
| `climate/crm/README.md` | Main highland/ring history plus electricity builds, corrected run, race bounds and diagnostics | Do not overwrite corrected climate or count storm windows twice |
| `atmosphere/README.md` | Existing column/escape context, final shield loss/cooling and electricity domain | Old prose says explicit cooling is absent; stage-1 or pre-IR values need historical labels |
| `habitation/README.md` | Main comfort range and author-adopted port/flyer possibilities; shield interfaces as future work | “No sizing” is stale, “safe/complete design” would overstate the new work |
| `engineering/README.md` | Independent tower/fire/flight and orbital-control/heat tools | Network reference ledger is not a current full ring-fleet industry budget |
| `visualization/README.md` | All viewers; current versus historical port forms; light calendar and evidence limits | Generated illustrations, reference renderings and measured outcomes are different |
| `shared/constants.*` | Sea albedo/phase and ExoPlaSim exports; shield proton mass, DE440 GM/system/time additions; main convention | Two `SUN_GM` aliases, real value changes, equivalent float formatting and literal duplicates need different treatment |
| `Makefile`, `.gitignore` | Sea calendar target and run/dist ignore patterns, all prior checks | `make check` must not lose any existing target; its normal calendar target may tolerate missing public-sky assets |
| `immersion/` | Existing branch-authored sky brightness fixes and labelled ring-0 artistic proposal | Merge does not reopen the frozen web engine or export new “science” from the experience |
| Locked/historical material | Exact ensemble trees and original provenance-pinned bytes | No paper role, UUID, locked seed, reference model or source-license cleanup during this integration |

Keep externally sourced bytes ignored and preserve manifests, credit and hashes. Consolidate metadata without vendoring papers, specialist datasets or private local configuration.

## Product and constants reconciliation

For each touched product, classify it before changing it:

1. **Current and unaffected:** its actual consumed model/input bytes and constants agree. Preserve it. A branch being older or the constants file gaining keys is not sufficient to regenerate it.
2. **Current but affected:** an adopted input or meaningful producer changed. Recompute with verified inputs into isolated outputs, compare scientifically meaningful values, then publish the product with its true producer/input identity.
3. **Historical by design:** preserve exact bytes and link the producing commit and actual inputs; current consumers must opt into its historical/scenario meaning rather than assume freshness.
4. **Unreproducible with present storage:** record missing artifacts and stop the affected claim. Do not suppress its check or replace a hash to imply reproduction.

| Product family | Required treatment |
|---|---|
| Main optical/sky/cloud products | Keep names-and-values provenance; sea already changes relevant Earthlight consumers. Verify exact inputs, don't globally rebuild because shield adds constants |
| Sea calendar | Sea+shield temporary tree matched 45 recorded producer-file entries and its recorded constants exclude changed `SUN_GM`; recheck full integrated tree, external assets and build. Preserve historical sea-calendar differences |
| Tides/waves/coasts | Original whole-file constant hashes remain honest run records. Solar GM changes by about 1.75×10⁻¹⁰ relatively and tides consume it. Pin historical values/reading rule or regenerate if claiming current-value outputs; never relabel an old run by updating hash alone |
| Shield magnetic architecture | Preserve original result bytes. Add a separate historical-input manifest mapping the original design-point path to its exact old fixture and producing commit; validate that chain and reading contract, identify it as historical at consumers, and retain both unavailable NPZ replay and modern rerun as unperformed. Alternative is genuine rerun from exact recovered closure/return histories. Simply accepting new hash, xfail, or deleting test is not the proposed repair |
| Shield ring/held eccentricity/heat | Preserve dependence chain and non-full-fleet interpretation; regenerate only if actual code/inputs changed. Even cosmetic changes to whole-file-hashed producers can invalidate freshness |
| Shield protection supply ledger | Keep its September catalogue scope; final loss hash does not make it the ring industry's ledger |
| Electricity corrected diagnostic | Record `94082be` producer identity (`6f158af5afeaa644`) versus later analysis script (`7eef6eaf777860d2` at `496ab7f`); preserve historical output or independently re-postprocess exact snapshots |
| Electricity conductivity/stage 1/aerosol/thunder/nitrogen | Add missing source/dependency contracts and schema checks at touched interfaces; preserve scheme/scenario distinctions. Historical identity supplements must be labelled as later curation, not original-run metadata. Use sidecars when preserving exact product bytes; attach only evidenced historical inputs, and mark unproved or inferred dependencies explicitly. A hash measured now alone does not prove historical consumption |
| Infrastructure studies | Recompute after explicit corrected paths and dependency/provenance repair; refresh sky-ships producer mismatch through real execution; include used constants and indirect source/product dependencies |
| Infrastructure geometry | The current concept form implementing the adopted architecture is the input for form-specific consumers, with its numerical dimensions still proposed. Preserve `form.py`’s fixed 8,193 m base; take crown dimensions from that form product rather than the corrected square lattice’s 9.3 km base. Keep square-lattice comparison separate and old fire outputs explicitly tied to old compartments; record unsolved chosen-form fire work |

Use shared physical constants without silently changing calibration, scenario durations, CM1 gas closures or legacy gravity. Electricity duplicates Boltzmann, charge, light speed and a rounded gas constant; infrastructure has shared physical literals alongside design choices. A migration should preserve each product's actual historical values and then quantify any change before publishing a current result. Record the common domain/model/input/units/evidence/reading-rule contract; tests should exercise meaningful interface failures, not mirror numerical implementations.

## Runtime inputs and bounded regeneration

The fresh worktree has no ignored model inputs or node environment. [Local state](evidence/local-state.json) records actual availability, not a claim that every archive byte has been verified. Choose interpreters by demonstrated dependencies; `climate/gcm/.venv/bin/python` is absent in the electricity worktree despite some README commands.

| Need | Restore or inspect | Avoid |
|---|---|---|
| Baseline repository checks | Research Python dependencies, optional numerical packages as required, immersion node dependencies and bake inputs; protection fetch manifest | Copying a virtualenv blindly; interpreting skip-heavy checks as complete coverage |
| Infrastructure input refresh | Two ignored main geography grids, committed corrected wind/ring/light/drainage products | Old infrastructure drainage under the same filename; new GCM/CM1 runs |
| Electricity focused checks/postprocessing | Pinned three electrical input tables, exact canonical case/snapshot/flash files, source archives for optional Fortran harnesses | Rebuilding shared `moon_omp_elec`, writing into archived runs, rerunning a 27-hour storm |
| Shield historical magnet treatment | Main's exact old `design_point.json` (expected SHA in shield report) and a documented historical fixture/consumer contract | Claiming old magnets consumed modern atmosphere |
| Shield modern magnetic rerun if separately chosen | Exact two closure NPZs, refined return NPZ, DE440, current design point | Guessing missing trajectories or substituting the unrelated expanded-cycle package |
| Sea calendar/viewer | Compact products plus pinned sky cache/atlas, Earth texture assets, star product and water/site inputs | Declaring a partial almanac build a complete public sky |
| Sea scientific rerender if later required | Wave spectra, coastal terrain rows, cloud arrays, Earthlight/water spectra, HDR-VDP/Octave | Treating a merge as a reason to render all 24 expensive scenes |

Exact infrastructure copies from main, at the same relative paths:

- `geography/products/atlas_28pct_4ppd.npz`: 9,875,080 bytes, SHA-256 `dc283827f349fdc134816e1a55deef3d608047fba9fbb7521a6005c171d81693`.
- `geography/products/drainage_28pct_16ppd.npz`: 2,697,170 bytes, SHA-256 `98842fc2a648f8899bd4f5d712445f78c8882c85a1bc30bf631ac9dd14ab1c03`.

After adopting corrected input paths (`site_winds_A28_dim5_moon_summit.json`, `global_winds_A28_dim5_moon.json`, `ring_ring_equator.json`) and repairing current geometry consumers, run sequentially:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.summit_tower.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.summit_tower.form
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.port_fire.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.sky_ships.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.sky_fleet.run
```

Compare the same-baseline cap regression values: design gust 27.3→32.4 m/s, square-lattice steel 32.8→46.0 Mt, chosen-form steel 90→98 Mt, and service sway 30→37 m. Those are recorded cap reruns, not numbers generated in this analysis. Any later geometry/provenance correction that changes results requires an explained delta; don't force exact agreement across changed assumptions. Preserve the chosen 1.5-million/24-km programme separately from optimized comparison outputs. Update READMEs/metropolis brief and relevant diagrams from the final products.

The old fire study may be rerun for corrected wind as a historical-geometry comparison. It does not become a calculation of the chosen wide-ring building. The new common geometry can expose its mismatch without forcing the integration to solve a new fire design.

Default to one numerical thread and serialized checks. Use existing caches only when identities match. If a necessary step turns out to require an expensive physical campaign, pause that step for a scoped decision with its concrete cost and missing input; continue unrelated integration work. No new long simulation is implicit in this plan.

## Verification and acceptance

Before each approved commit, run `make check` with suitable environment and thread limits. Once the final content is settled, require a final complete run against that exact tree; previous branch logs are comparison evidence. If Make stops early, record which targets were not reached and run unaffected downstream checks separately where useful. Do not hide failures by changing discovery, loosening numerical tolerances or removing provenance assertions.

Focused suites and exact historical/actual check distinctions are in the branch reports. Add only necessary meaningful checks for changed product boundaries or geometry consumers. For complete sky delivery run `python visualization/light-calendar/build.py --require-assets`, then its domain and visualization JavaScript tests. Review desktop/phone interaction at an HTTP origin if the build is included as a usable deliverable; this analysis has not provided browser review.

Report `research/check.py`'s protection state verbatim. The strict command `python research/check.py --require-inputs` must pass before claiming restoration/reproduction of the historical protection-input set it checks. Other spectral, DE440, electrical and terrain inputs require their own manifests; this command does not verify all branch inputs. A normal `make check` can return successfully with `PASS_WITH_BLOCKED_INPUTS`; if so, the final handoff must explicitly retain that block. Missing optional modules and skipped Fortran/input tests likewise limit coverage.

The final unified branch is ready for further work only when:

- All four pinned source heads are ancestors of it; main's corrections and branch caps are present; no extra source branch has been silently selected.
- No unresolved conflicts, unintended deletions, untracked deliverables outside the agreed analysis/implementation, or uncommitted intended changes remain after approved commits.
- Every topic's current state, adopted decisions and next-work pointer are coherent; superseded results are labelled, including the old twilight, tower, lightning and shield values.
- Product consumers either use verified current inputs or explicitly documented historical/scenario contracts. The known magnetic test cannot remain an unexplained failure; unresolved missing reproduction is separately declared.
- Corrected infrastructure products are real reruns, with geometry versions explicit and no false claim of a solved current fire design.
- Selected runtime inputs/build outputs are reproducibly described and writes are isolated; all five source worktrees and shared archived executables remain unchanged.
- Required checks ran on the final tree with exact pass/fail/skip/block states. A conditional check pass is described as such, never “everything validated”.
- Locked ensemble/provenance trees retain their recorded identity, with exactly five paper roles; no source-admission, manuscript or PDF approval was fabricated.
- The handoff states current branch/commit, completed integrations, historical datasets retained, missing runtime inputs and ranked next research choices. Remote publication is not part of this plan.

## Review package checks and limits

The analysis has already established branch identity, pairwise conflicts and several concrete provenance defects. It has not performed this cumulative integration or its final tests. [checks.json](evidence/checks.json) separates tests actually run during analysis from historical logs and deferred work. [audit.py](audit.py) repeats source identity, locked-tree, pairwise and artifact-structure checks with temporary Git objects, without creating commits or a merged workspace. The snapshot's temporary tree OIDs are diagnostic outputs; reproduction relies on pinned source commits, not retaining those temporary objects.
