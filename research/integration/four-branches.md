# The joint integration of 1 November 2026: four branches after the design run

Five states of the work became one on 1 November 2026: main at `25a98dd`, climate at `a5de646`, atmosphere at `188e03a`, ecology at `e953c19` and provisioning at `b26125c`. The register scheduled the reconciliation for when the design run was complete and analysed ([register](../decisions.md#the-cm1-boxes-rain-and-plasims-convection-21-october-2026), 30 October); since the gap between the GCM and CM1 stayed significant after the run, the branches are integrated and aligned before the work returns to that gap ([register](../decisions.md#the-design-run-31-october-2026), 1 November). Array-industry, Jura Island and resources stood at `98dfcc6` with no work of their own. The merges ran on `research/joint-integration-5` in the worktree `../terluna-joint-integration-5`, climate first, then atmosphere, ecology and provisioning, each an ordinary merge commit by exact head, so every source head is an ancestor of the result. Main moves up to it by fast-forward once it is checked and reviewed, and every branch is then brought level with main.

This round has no synthesis study. The studies that read the design climate, the biome maps, local climates and the Jura box among them, wait until the GCM and CM1 agree, and the one joint calculation is the cold cloud's sensitivity on the design years (below), which the register asked of the reconciliation.

## The commits

| Commit | What it does |
|---|---|
| `f91d450` | Merges climate: the GCM's lunar physics since 20 October and the design run, whose years 10–19 are the design climate; the fall law for CM1's rain rebuilt on main's rebound products, its values unchanged |
| `025ff80` | Merges atmosphere: the snow work, the rain solver's Massman diffusivity and the sub-cloud closure; the two shared READMEs joined, the shared rain and snow products rebuilt, the frozen snow records re-pinned |
| `65d1ce9` | Merges ecology: the biomass comparison and the deep sky's production, rebuilt on the merged products with their values unchanged |
| `c5aa9ac` | Merges provisioning: the aerophytes' uses, the world economy and the monetary system; the status file merged field by field |
| `488d627` | The cold cloud's sensitivity on the design years and the records: this account, the register, the plan, the status, the climate README and this folder's index |
| the answering commit | The independent review's findings answered (below) |

## How the shared documents were joined

Git combined most of the tree. Where both sides changed a document, the resolution kept every branch's entries, as atmosphere's merge notes asked, except climate's descriptions of snow work that atmosphere's later work had superseded:
- **atmosphere/README.md.** The table takes atmosphere's Rain row, with Massman's diffusivity, the breakup scale and CM1's rain limits, and its Snow row joined with climate's slabs of cold air from cirrus.py and the high cloud's ice for the GCM; climate's description of the superseded super-particle column (Phillips, Yano and Khain 2017) and its flakes 15–90 times Earth's gives way to McSnow's.
- **atmosphere/snow/README.md.** The file table takes atmosphere's rows with climate's cirrus.py row beside them, and the investigation's index keeps atmosphere's table. Climate's open items for the cold cloud join the snow work's open list as one item; its open items for the snow itself, the spread of like flakes' fall speeds and a turbulence limit on large flakes, were since taken up by McSnow's aggregation and the turbulence analogy.
- **research/status.json.** Every branch's topics, provisioning's three beside ecology's, merged field by field (`session-tools/integration3/status_merge.py`).
- **The shared rain and snow products.** The tree takes atmosphere's rain drops, arrivals and limits and its snow flakes and aggregation, built on the measured diffusivity, with main's rain rates; settle rebuilt them against main's climate targets with every value unchanged. Climate's fall law for CM1, rebuilt on the new solver, moves only in its two diagnostic ratios to the solver, of evaporation from 0.9613 to 0.9616 and of loss per metre from 0.9085 to 0.9087.
- **Climate's files in atmosphere/.** Methane's lifetime below the 0.3 Pa base (loss_response, the hydrogen escape and the middle atmosphere's README) merged as climate wrote it, replacing atmosphere's earlier call for methane's chemistry in the column.

## The snow records' pins

Forty-seven files on the atmosphere branch pinned what the branch held: thirty-six frozen snow records and one note the climate targets (`80a85def`), and twelve records the shared constants (`dc8713f4`), two of them both. The sub-cloud closure keeps its own `80a85def` binding byte for byte, as the branch's merge notes asked, since PlaSim's snow build records the closure's hash. Main's cascades had since rebound the targets with every value the same, and climate had added constants beside the old ones. Several of the records take hours to rebuild (McSnow, SCALE-SDM, the event columns) or need run folders settle does not hold, so their pins were rewritten instead, after two checks: the targets at `80a85def` and in the merged tree hold the same values with hashes masked, and every entry of the branch's constants stands unchanged in the merged file.

The records pin one another, two snow modules carry their baselines' hashes in code, and the snow assessment binds two products' runs and original content by digest. So the new hashes were carried on through every file that recorded the old ones, in eight rounds until none changed, and the assessment's two digests were recomputed with the assessor's own functions, its record of the bytes first assessed kept. Only hash strings changed. The proof, with every file and replacement, is on the research drive (`session-tools/design_run/pin_rewrite_2026-11-01.json`). Twenty dated check snapshots, eighteen of atmosphere's and two of provisioning's, record the hashes of what was checked when they were taken; the rewrite first carried them forward too, and they are restored to their branches' versions, which no file outside the snapshots pins.

## The world economy's inputs

Provisioning's world economy binds the products it reads by hash in its manifests (`inputs.json`, `amenities_inputs.json` and those of its later modules) and refuses a changed input until it is reviewed and rebound explicitly. Six of main's products (the distributions, the power screen, the food blend, the terrain through time, the carbon cycle and the sea phosphorus) had been rebound by main's cascades since the branch began, each keeping every value once hashes are masked, the fields the study selects among them. They were rebound explicitly: the new hashes were carried through the study's manifests and products, and through the monetary system and provisioning's files that pin them, in nine rounds until none changed, and the construction product, which records a digest of its whole loaded inputs, was regenerated by its own producer with its values unchanged. The proof is on the research drive (`session-tools/design_run/world_economy_rebind_2026-11-01.json`).

## The settle tool

The settle tool on the research drive (`session-tools/settle.py`, its earlier states kept as `settle.py.bak_2026-11-01` and its `…b` to `…g` successors) takes four new producers, each after what it binds: climate's fall law for CM1 after the rain limits, ecology's biomass comparison after the light limitation, its aerial production after the gases-glow study, and provisioning's aerophyte products after the sky's products.

It also gains follow-only files, listed in `follow_only.txt` beside it: the fifty-four records, modules and note of the snow chain and the thirty-one files of the world economy and the monetary system. Settle never rebuilds the snow records, whose producers take hours. It rebinds the follow-only files' JSON pins to products a cascade keeps, checking each rebind as it checks its own products'. In every follow-only file it then rewrites each earlier digest of such a product as a whole hex token, however the pin quotes it, which reaches records that name their inputs relative to their own folder and the note that quotes its pin in backticks. It recomputes the snow assessment's digests (`reassess_snow.py`) only on products whose values held through the run and whose runs equal the committed version's with hashes masked, and regenerates the world economy's construction and comparison, which record a digest of their whole inputs, with their own producers. It repeats these steps until nothing changes. A follow-only file still pinning a product whose values changed in the run, a refused assessment, a regeneration its manifest refuses, and an earlier digest left in any follow-only file each stop it with exit 2 at its end, for a rebuild or a decision by hand; pins to products changed outside the run are listed as notes.

Provisioning's branch asked for the world economy's seventeen producers and the monetary system's run in settle's order, each rebound in its readers' manifests as it rebuilds. Settle does not read pins of the manifests' form, and the producers refuse a changed input until it is rebound by hand, so the integration follows the manifests and products instead and regenerates the two products that digest their whole inputs. The chain keeps its explicit guard: a cascade that changes a value it reads stops at it, and one that keeps every value rebinds it without rebuilding seventeen producers.

Copies of the tree tested the tool, the logs beside them on the research drive (`session-tools/design_run/t2`, `t3` and `cascade_test`). With a one-character register edit, the cascade rebuilt 6 products, four of main's chain and the world economy's construction and comparison, and rebound 197 files in about fourteen minutes while a second test shared the machine; the 1,851 tests of the atmosphere, the world economy, the monetary system and provisioning then passed, and a second run found nothing stale. With the terrain study's rail energy raised from 107 to 108 kJ a tonne-kilometre, settle rebuilt the terrain product with its values changed, the world economy's manifest refused the new product, and settle reported the two regenerations for a rebind by hand and exited 2. In the tree itself, a snow product with one number of its runs changed, restored afterwards, was refused by the assessment helper, exit 2. A first test, on a copy at `025ff80` with the tool as it then stood, had rebuilt 4 products and rebound 110 files in about seven minutes. The records commit's own register edit ran the cascade on the merged tree, rebuilding 6 products and rebinding 196 files in about nine minutes, every value unchanged, and the answering commit's edits ran it again, rebuilding 6 and rebinding 197.

## The cold cloud's factors on the design years

When the design run's diffusivity limit was set, the register asked that the cold cloud's factors be redone at the reconciliation with Massman's measured diffusivity, as a sensitivity on the design years ([register](../decisions.md#the-design-run-31-october-2026)). They were computed twice on the design run's years 10–19: from the merged tree, with Massman's diffusivity, and from the climate merge's tree, whose rain solver still diffuses at the Chapman–Enskog rate (`high_cloud.py` and `cold_cloud.py`, four workers, about 17 minutes each). Each set's first-order effect is read on the design run's own cloud:

| Factors, layers 1–3 | Water | Radius | First-order effect |
|---|---|---|---|
| Adopted (the 20 October years, Chapman–Enskog) | 6.09, 36.61, 42.07 | 1.00, 1.68, 1.60 | +9.86 W/m², 6.3–8.2 K at the ground |
| The design years, Chapman–Enskog | 6.10, 32.40, 43.10 | 1.00, 1.68, 1.61 | +10.13 W/m² |
| The design years, Massman | 6.09, 32.06, 43.13 | 1.00, 1.72, 1.59 | +10.18 W/m² |

On the design years the diffusivity moves the factors by under 2.5%, where it moved the second layer's water 7.7% on the 20 October years; the change of climate moves that water 11.5%. Taken together, the factors recomputed on the design years would warm the design climate by 0.33 W/m² at first order, 0.21–0.27 K at the ground, against the 0.41 K by which the design years lie below the target. The GCM's own response can differ from a fixed-cloud estimate, as the first design run's clouds showed. The products are [cold_cloud_A28R_301k_design_l15_snow_s934.json](../../climate/results/gcm/cold_cloud_A28R_301k_design_l15_snow_s934.json) and [high_cloud_A28R_301k_design_l15_snow_s934.json](../../climate/results/gcm/high_cloud_A28R_301k_design_l15_snow_s934.json); the earlier diffusivity's twin and the comparison (`factor_effects.py`) are on the research drive in `session-tools/cold_cloud_design_years`. The design runs carry the adopted factors; whether the design case takes the recomputed ones stays open.

## Outside Git

- The worktree links main's ignored inputs, the research drive's run folders, atmosphere's twelve snow run folders and the immersion's `node_modules`. Only atmosphere's worktree held the twelve snow run links before; main's checkout and the other branches' worktrees take them when they move up, so the snow tests that read them run everywhere, and every checkout's `make check` then counts 24 skipped.
- The settle tool's changes above, `follow_only.txt` and `reassess_snow.py`, with their earlier states beside them.

## Checks

Each merge, the records commit and the answering commit passed `make check` before its commit:

| Merge | Python tests | The rest |
|---|---|---|
| Climate | 3,098 passed, 24 skipped | JavaScript, calendar, research and ensemble checks passed, once the fall law was rebuilt and the immersion's `node_modules` linked |
| Atmosphere | 4,132 passed, 24 skipped | all passed, after the snow chain's re-pinning |
| Ecology | 4,149 passed, 24 skipped | all passed |
| Provisioning | 4,633 passed, 24 skipped | all passed |
| The records commit, after its register edit's cascade | 4,633 passed, 24 skipped | all passed |
| The answering commit, after its cascade | 4,633 passed, 24 skipped | all passed |

## The review

An independent reviewer read the merged tree at `488d627` without changing it, before main moved (its findings on the research drive, `session-tools/integration5_review/findings.md`). It found nothing blocking, twelve points to fix and fifteen minor ones, and confirmed that every branch's work is in the tree, that both rebinds change only hash strings and that every branch head is an ancestor of the result. The answering commit takes them up:
- **The settle tool**: digests followed as whole hex tokens in any quoting, the snow note's pin set to the merged targets, the assessment recomputed only where nothing but hashes moved, a refused regeneration reported rather than crashing the run, follow-only rebinds checked, and a closing search for any earlier digest left behind (above, with the tests).
- **The records' numbers**: the snow pins counted by what each file pinned; the settling fits' feedbacks given their right signs; the run's stop placed in its year 5, just after year 4 was written; the ground's range through year 11 (300.87–301.05 K) and the soil water's (0.151–0.159 m); the step's sea clouds by year (1.25 W/m² more by year 15, about 0.65 more by years 17–19); the air temperatures as land and all-air means; the sunlight run's cloud changes as means over the land and over the sea with their shares of the global 0.99 W/m²; the adopted third-layer radius 1.60; the cascade test's time.
- **The register's economy rows**: the per-area ratio's waters, the tenant clusters as payload for the mass of a home and of a community, the build's power as purchased process power, mostly at the source sites, flame-spread resistance as credible but unsized, and the biomass study's rerun after the models agree.
- **The departures from provisioning's merge notes**, recorded above, with its citation of the aerial-productivity study given by path and its deferred items carried below.
- **The status file**: the merged topics' next tasks after the merge, atmosphere's merge list marked done, the climate topic's cover and dry-land figures for the design years.
- **Editorial points** in this account, the register, the plan and the climate README.

## Afterwards

Main moves up to the integration branch's head by fast-forward, and every branch then moves to main by fast-forward, each one's head being an ancestor of the result: climate, atmosphere, ecology and provisioning, merged here, and array-industry, Jura Island and resources from `98dfcc6`.

## Still to do from the integration

- **The gap with CM1**, from CM1's sea and land boxes forced by the design years, before the design climate's readers run ([register](../decisions.md#the-design-run-31-october-2026)).
- **Storm water running off dry ground** into channels and hollows, which PlaSim's bucket soil cannot shed below its capacity, soon after this integration ([register](../decisions.md#the-design-run-31-october-2026)).
- **Whether the design case takes the cold cloud's factors** recomputed on the design years (above).
- **The cost of McSnow's cap**, up to 3% of a run's time, trimmed before the next long McSnow batch ([snow README](../../atmosphere/snow/README.md)).
- **The world economy's deferred items**, from provisioning's merge notes: moving its engineering models to `engineering/`, which redesigns the modules and their products since each mixes physical and cost calculations and binds the next by hash; splitting the build's energy out of `construction.py` into resources; sizing interworld travel; and the editorial pass on the industrial write-ups' bodies.
- **Settle's dry run**, which still repins the in-tree pins before it lists what is stale.
