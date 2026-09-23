# Reviewed forest-patch checkpoint

Publication record: 23 September 2026. This directory preserves the results
reviewed in the lunar-megaforest conversation, including the perimeter-root
experiment and the correction about prescribed canopy shapes.

## Two independent experimental lineages

The remote repository already contained a different implementation at
[`879bbb9`](https://github.com/dgoldman0/Terluna/commit/879bbb98257071fe73db3c9d988b4f3f811a9904).
The saved conversation checkpoint was based on the same earlier individual-tree
commit, but uses incompatible definitions under some of the same filenames.
Publication preserves both implementations. The existing working files on `main`
are unchanged; the reviewed source is an additional ancestor of the publication
merge and is also retained on `checkpoint/forest-patch-review-57-tests`.

| Property | Reviewed conversation checkpoint | Previously published channel-flow checkpoint |
|---|---|---|
| Exact source commit | `363b161e7c92be54e7f4907b2186bfc1bace10f1` | `879bbb98257071fe73db3c9d988b4f3f811a9904` |
| Main patch spacing | 90 m | 105 m |
| Intact comparison wind | 20 m/s reservoir target | 15 m/s imposed upstream wind |
| Flow boundary treatment | Periodic horizontal domain with relaxation reservoir | Confined channel with fixed normal inflow/outflow |
| Aboveground mass matching | Specified tissues/substrate; excludes retained water | Includes selected retained water |
| Focused tests in checkpoint | 57 | 62 |
| Structural comparisons | 17 | 39 |

Read [methods.md](methods.md), [checkpoint.json](checkpoint.json),
[shape_comparison.csv](shape_comparison.csv) and [sources.json](sources.json) as
records of the **reviewed** source commit. The other implementation's methods
remain at [research/forest_patch.md](../README.md). Parameters, thresholds,
test counts and result tables must stay with their originating source revision.
The reviewed files' historical `NOT_PUBLISHED` fields describe the earlier
blocked delivery and are superseded by this publication record.

## Canopy-shape interpretation

The dome was specified as one candidate after the author's suggestion, alongside
flat, tapered, ramped, irregular and gapped layouts. Its geometry was an input.
The analysis computes conditional responses; growth, reproduction and evolutionary
selection are unimplemented.

At the same 20-m/s reference and matched aboveground tissue/substrate mass,
maximum anchorage demand/capacity is 0.642 for the flat independent stand, 0.701
for tapered margins and 0.726 for the dome. The ramp gives 0.719 with wind toward
the rising canopy and 0.620 from the reversed direction. These selected cases
support retaining a flat canopy as a serious baseline; they establish no general
optimum or evolved forest shape.

Keeping the dome's perimeter root plates at least 30 m in radius and 4 m deep
reduces its maximum ratio to 0.599 without altering the aboveground geometry.
The summed idealized root-plate volume increases from about 549,851 to 611,513 m3.
That is additional belowground investment whose tissue construction cost remains
unmodeled. An equivalent reinforcement/resource-budget comparison for the other
shapes is required before assigning an advantage to the dome itself. Earlier
emphasis on this improved dome should not be read as model selection of a dome.

The transferable question is how height, crown porosity, stem taper and rooting
should vary with position. Gust forcing, broader domains, nonlinear mechanics,
complete resource budgets and regeneration remain discriminating next steps.
Damage-sequence stopping counts in this checkpoint are unresolved because of
load-increment dependence and the small-displacement guard.

## Reproduce the exact reviewed implementation

Use the saved commit in a separate worktree so the differing `main` modules and
runner are never mixed into this numerical experiment:

```sh
git fetch origin
git worktree add --detach ../Terluna-forest-patch-reviewed 363b161e7c92be54e7f4907b2186bfc1bace10f1
cd ../Terluna-forest-patch-reviewed
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test*forest*.py' -v
OPENBLAS_NUM_THREADS=1 python research/run_forest_patch.py
```

The runner reproduces the case, sensitivity, damage and extra-study stages,
including individual-tree tables, flow arrays, the 125-tree comparison and the
perimeter-root experiment. Follow its own methods for optional stage selection.
Large generated numerical arrays and diagnostic images remain in the original
checkpoint bundle or can be regenerated; they are not added as opaque archives
to Git. This follows the repository's compact-results practice.

## Publication checks

All 60 entries of the saved bundle match its manifest hashes. Its ten-file source
patch is preserved byte-for-byte: the independently reconstructed Git tree is
`1b9aa73da4950fa23287933231a1578bce7ada7a`. The 57 focused tests were rerun during
publication and passed. The full simulation study was not rerun for publication;
its original numerical records remain unchanged. Full-repository environmental
and protection regressions remain unrun.

See [publication.json](publication.json) for exact source hashes and test/runtime
scope. Existing `main` solver files, its tests/results, historical baselines,
manuscript seeds and `immersion/` are preserved. No force push is used.
