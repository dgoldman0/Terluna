# The joint integration of 1 November 2026: four branches after the design run

Five states of the work became one on 1 November 2026: main at `25a98dd`, climate at `a5de646`, atmosphere at `188e03a`, ecology at `e953c19` and provisioning at `b26125c`. The register set the reconciliation for once the design run was complete and analysed ([register](../decisions.md#the-cm1-boxes-rain-and-plasims-convection-21-october-2026), 30 October), and with the gap between the GCM and CM1 still significant after it, the branches are integrated and aligned before the work returns to that gap ([register](../decisions.md#the-design-run-31-october-2026), 1 November). Array-industry, Jura Island and resources stood at `98dfcc6` with no work of their own. The merges ran on `research/joint-integration-5` in the worktree `../terluna-joint-integration-5`, climate first, then atmosphere, ecology and provisioning, each an ordinary merge commit by exact head, so every source head is an ancestor of the result. Main moves up to it by fast-forward once it is checked and reviewed, and every branch is then brought level with main.

This round has no synthesis study. The studies that read the design climate, the biome maps, local climates and the Jura box among them, wait until the GCM and CM1 agree, and the one joint calculation is the cold cloud's sensitivity on the design years (below), which the register asked of the reconciliation.

## The commits

| Commit | What it does |
|---|---|
| `f91d450` | Merges climate: the GCM's lunar physics since 20 October and the design run, whose years 10–19 are the design climate; the fall law for CM1's rain rebuilt on main's rebound products, its values unchanged |
| `025ff80` | Merges atmosphere: the snow work, the rain solver's Massman diffusivity and the sub-cloud closure; the two shared READMEs joined, the shared rain and snow products rebuilt, the frozen snow records re-pinned |
| `65d1ce9` | Merges ecology: the biomass comparison and the deep sky's production, rebuilt on the merged products with their values unchanged |
| `c5aa9ac` | Merges provisioning: the aerophytes' uses, the world economy and the monetary system; the status file merged field by field |
| the records commit | The cold cloud's sensitivity on the design years and the records: this account, the register, the plan, the status, the climate README and this folder's index |

## How the shared documents were joined

Git combined most of the tree. Where both sides changed a document, the resolution kept every branch's entries, as atmosphere's merge notes asked:
- **atmosphere/README.md.** The table takes atmosphere's Rain row, with Massman's diffusivity, the breakup scale and CM1's rain limits, and its Snow row joined with climate's slabs of cold air from cirrus.py and the high cloud's ice for the GCM.
- **atmosphere/snow/README.md.** The file table takes atmosphere's rows with climate's cirrus.py row beside them, and the investigation's index keeps atmosphere's table. Climate's open items for the cold cloud join the snow work's open list as one item; its open items for the snow itself, the spread of like flakes' fall speeds and a turbulence limit on large flakes, were since taken up by McSnow's aggregation and the turbulence analogy.
- **research/status.json.** Every branch's topics, provisioning's three beside ecology's, merged field by field (`session-tools/integration3/status_merge.py`).
- **The shared rain and snow products.** Atmosphere's rain drops, arrivals and limits and its snow flakes and aggregation, built on the measured diffusivity, and main's rain rates; settle rebuilt them against main's climate targets with every value unchanged. Climate's fall law for CM1, rebuilt on the new solver, moves only in its two diagnostic ratios to the solver's evaporation, 0.9613 to 0.9616 and 0.9085 to 0.9087.
- **Climate's files in atmosphere/.** Methane's lifetime below the 0.3 Pa base (loss_response, the hydrogen escape and the middle atmosphere's README) merged as climate wrote it, replacing atmosphere's earlier call for methane's chemistry in the column.

## The snow records' pins

Forty-six frozen snow records and one note on the atmosphere branch pinned the climate targets as the branch held them (`80a85def`), and twelve of them the shared constants (`dc8713f4`). Main's cascades had since rebound the targets with every value the same, and climate had added constants beside the old ones. Several of the records take hours to rebuild (McSnow, SCALE-SDM, the event columns) or need run folders settle does not hold, so their pins were rewritten instead, after two checks: the targets at `80a85def` and in the merged tree hold the same values with hashes masked, and every entry of the branch's constants stands unchanged in the merged file.

The records pin one another, two snow modules carry their baselines' hashes in code, and the snow assessment binds two products' runs and original content by digest. So the new hashes were carried on through every file that recorded the old ones, in eight rounds until none changed, and the assessment's two digests were recomputed with the assessor's own functions, its record of the bytes first assessed kept. Only hash strings changed. The proof, with every file and replacement, is on the research drive (`session-tools/design_run/pin_rewrite_2026-11-01.json`). Twenty dated check snapshots, eighteen of atmosphere's and two of provisioning's, record the hashes of what was checked when they were taken; the rewrite first carried them forward too, and they are restored to their branches' versions, which pin only one another.

## The world economy's inputs

Provisioning's world economy binds the products it reads by hash in its manifests (`inputs.json`, `amenities_inputs.json` and those of its later modules) and refuses a changed input until it is reviewed and rebound explicitly. Six of main's products (the distributions, the power screen, the food blend, the terrain through time, the carbon cycle and the sea phosphorus) had been rebound by main's cascades since the branch began, each with every value the same with hashes masked and the fields the study selects equal. They were rebound explicitly: the new hashes were carried through the study's manifests and products, and through the monetary system and provisioning's files that pin them, in nine rounds until none changed, and the construction product, which records a digest of its whole loaded inputs, was regenerated by its own producer with its values unchanged. The proof is on the research drive (`session-tools/design_run/world_economy_rebind_2026-11-01.json`).

## The settle tool

The settle tool on the research drive (`session-tools/settle.py`, its earlier states kept as `settle.py.bak_2026-11-01` and its `…b` to `…e` successors) takes four new producers, each after what it binds: climate's fall law for CM1 after the rain limits, ecology's biomass comparison after the light limitation, its aerial production after the gases-glow study, and provisioning's aerophyte products after the sky's products.

It also gains follow-only files, listed in `follow_only.txt` beside it: the fifty-four records, modules and note of the snow chain and the thirty-one files of the world economy and the monetary system. Settle never rebuilds the snow records, whose producers take hours. It rebinds the follow-only files' JSON pins to products a cascade keeps, rewrites their quoted hashes in whatever form a pin takes (some records name their inputs relative to their own folder), recomputes the snow assessment's digests (`reassess_snow.py`), regenerates the world economy's construction and comparison, which record a digest of their whole inputs, with their own producers once those inputs are followed, and repeats until nothing changes. A follow-only file still pinning a product whose values changed in the run stops it with exit 2, for a rebuild or a decision by hand. A copy of the merged tree at `025ff80`, given a one-character register edit, tested it: the cascade rebuilt 4 products and rebound 110 files in under five minutes, all 1,326 atmosphere tests passed after it, and a second run found nothing stale. The records commit's own register edit then ran the cascade on the merged tree: it rebuilt 6 products, four of main's chain and the world economy's construction and comparison, and rebound 196 files in about nine minutes, every value unchanged.

## The cold cloud's factors on the design years

When the design run's diffusivity limit was set, the register asked that the cold cloud's factors be redone at the reconciliation with Massman's measured diffusivity, as a sensitivity on the design years ([register](../decisions.md#the-design-run-31-october-2026)). They were computed on the design run's years 10–19 from the merged tree, and beside it from the climate merge's tree, whose rain solver still diffuses at the Chapman–Enskog rate (`high_cloud.py` and `cold_cloud.py`, four workers, about 17 minutes each). Each set's first-order effect is read on the design run's own cloud:

| Factors, layers 1–3 | Water | Radius | First-order effect |
|---|---|---|---|
| Adopted (the 20 October years, Chapman–Enskog) | 6.09, 36.61, 42.07 | 1.00, 1.68, 1.61 | +9.86 W/m², 6.3–8.2 K at the ground |
| The design years, Chapman–Enskog | 6.10, 32.40, 43.10 | 1.00, 1.68, 1.61 | +10.13 W/m² |
| The design years, Massman | 6.09, 32.06, 43.13 | 1.00, 1.72, 1.59 | +10.18 W/m² |

On the design years the diffusivity moves the factors by under 2.5%, where it moved the second layer's water 7.7% on the 20 October years; the change of climate moves that water 11.5%. Taken together, the factors recomputed on the design years would warm the design climate by 0.33 W/m² at first order, 0.21–0.27 K at the ground, against the 0.41 K by which the design years lie below the target. The GCM's own response can differ from a fixed-cloud estimate, as the first design run's clouds showed. The products are [cold_cloud_A28R_301k_design_l15_snow_s934.json](../../climate/results/gcm/cold_cloud_A28R_301k_design_l15_snow_s934.json) and [high_cloud_A28R_301k_design_l15_snow_s934.json](../../climate/results/gcm/high_cloud_A28R_301k_design_l15_snow_s934.json); the earlier diffusivity's twin and the comparison (`factor_effects.py`) are on the research drive in `session-tools/cold_cloud_design_years`. The design runs carry the adopted factors; whether the design case takes the recomputed ones stays open.

## Outside Git

- The worktree links main's ignored inputs, the research drive's run folders, atmosphere's twelve snow run folders and the immersion's `node_modules`, as the branches' worktrees do. Main's checkout takes the twelve snow run links when it moves up, so the snow tests that read them run there too.
- The settle tool's changes above, `follow_only.txt` and `reassess_snow.py`.

## Checks

Each merge, and the records commit, passed `make check` before its commit:

| Merge | Python tests | The rest |
|---|---|---|
| Climate | 3,098 passed, 24 skipped | JavaScript, calendar, research and ensemble checks passed, once the fall law was rebuilt and the immersion's `node_modules` linked |
| Atmosphere | 4,132 passed, 24 skipped | all passed, after the snow chain's re-pinning |
| Ecology | 4,149 passed, 24 skipped | all passed |
| Provisioning | 4,633 passed, 24 skipped | all passed |
| The records commit, after its register edit's cascade | 4,633 passed, 24 skipped | all passed |

## The review

An independent reviewer reads the merged tree at the records commit, without changing it, before main moves; its findings and their answers follow in a commit of their own.

## Afterwards

Main moves up to the integration branch's head by fast-forward, and every branch then moves to main by fast-forward, each one's head being an ancestor of the result: climate, atmosphere, ecology and provisioning, merged here, and array-industry, Jura Island and resources from `98dfcc6`.

## Still to do from the integration

- **The gap with CM1**, from CM1's sea and land boxes forced by the design years, before the design climate's readers run ([register](../decisions.md#the-design-run-31-october-2026)).
- **Storm water running off dry ground** into channels and hollows, which PlaSim's bucket soil cannot shed below its capacity, soon after this integration ([register](../decisions.md#the-design-run-31-october-2026)).
- **Whether the design case takes the cold cloud's factors** recomputed on the design years (above).
- **The cost of McSnow's cap**, up to 3% of a run's time, trimmed before the next long McSnow batch ([snow README](../../atmosphere/snow/README.md)).
