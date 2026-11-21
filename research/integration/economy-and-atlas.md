# The joint integration of 21 November 2026: the vacuum buildout's economy and the atlas round

Five states of the work became one on 21 November 2026: main, at `58286a18` when the round was called and at `2257668e` with the register's row on it, the vacuum economy's branch at `97177c89`, resources at `346e17f0`, hydrosphere at `1baa2658` and geography at `f1b27235`. The merges ran on the standing integration branch `maintenance/integration` in the worktree `../terluna-integration`, fast-forwarded first from the eighth round's head, `4b31673f`, to main. Each branch's merge is an ordinary merge commit by exact head, in the order economy, resources, hydrosphere, geography: resources' salt study reads the economy's source-access product, and the other two read neither.

The round takes in all the work beyond main, the economy's with it. The register's condition of 15 November, that main takes none of the refactored economy before its report, was met by the report written on 17 and 18 November, and the register's row of 21 November sets the round ([register](../decisions.md#the-atlas-round-and-the-vacuum-economy-closed-for-now-and-the-ninth-integration-21-november-2026)). Both studies were closed for now before the round, each with a section in its README that says where it stands and what is open: [the atlas round](../studies/atlas_round/README.md#where-the-round-stands-21-november-2026) and [the period model of the economy](../studies/vacuum_economy_periods/README.md#where-the-study-stands-21-november-2026).

`domain/provisioning` (`adeaa8ea`) and `domain/humanities` (`c077b606`) come in inside the economy's branch, which merged them on 15 November with their history kept, as do the eight research branches made during the returns work, consolidated-returns, credit-supported-expansion, export-strategies, price-pressure, staged-complex-returns, u-monetary-pressure, u-policy-risk and u-trade-transition, whose commits all lie inside `domain/humanities`. Atmosphere, climate, ecology, food-systems, geology, array-industry and Jura Island stood at main's head of 16 November with no work of their own.

The round has no synthesis study and no joint computation. The economy's settlements are those of the vacuum period, before the water, and no product of one study reads the other's; the one place they touch is named under what is still to do.

## The commits

| Commit | What it does |
|---|---|
| `2862f0dc` | Merges the vacuum economy's branch, 86 commits and 487 files: the studies of 12 to 14 November on the crater towns and tube cities, the staged industry and its power, the sources, sport and entertainment, the returns, U's money, wealth and income and the product routes, as they were written; and the three models built beside them, the coupled study, the compact model and the period model with its report and its fund proposal |
| `e4fb7660` | Merges resources: the salt for an ocean-like South Pole–Aitken Sea, seven files, with the salt product regenerated on the merged tree |
| `12f1dd03` | Merges hydrosphere: two notes in the south sea study's README on its paused oxygen run |
| `e6759adb` | Merges geography, 75 commits of its own among 118 and 160 files: the atlas round, from the surface redrawn under its water to regions, districts, names and the globe on the design surface, with the climate domain's pattern scaling and surface-change estimate, geography's mare map and surface composition and the water-surface model's open sea |
| the records | This account, the register's section of 12 to 14 November and its statuses brought up, the plan, the research index, the status file and this folder's index |

## The conflicts

Three merges had none. The economy's branch and main had both changed the register and the plan since 15 November, and Git joined them by lines: the branch's sections of 12 November beside main's rows since.

Resources' merge conflicted in ten files, all in `resources/source_access`: two input records, the start-up deliveries' inputs and seven products. On the branch they differ from the common base only by the fix that binds the protection gates to the register lines they cite (`770dea09`), picked across from the economy's branch (`b3e0010d`) on 15 November so the branch's tests could run. The economy's branch carries that fix and regenerated the same products twice more, when the copies between studies became direct reads. The merge keeps the economy's versions of all ten. The salt product had recorded the potential product's digest as resources held it; regenerated on the merged tree, it has every value the branch computed and only that one recorded digest differs (`sea_salt_as_merged.json` on the research drive holds the branch's product for the comparison).

## Attribution

Thirty-five commits of 12 to 14 November were written with OpenAI Codex (GPT-6 Astra) and carry no co-author line of its own (register, 15 November). No history is rewritten. The three merges into the economy's branch name them and carry the co-author line, the economy's merge here repeats it, and this is the list:

- on `domain/provisioning`: `3e045911`, `a6df1bec`, `c9e37d59`, `210a20bd`, `a8af8761`, `fd431684` and `adeaa8ea`;
- on `domain/resources`: `38cd446f`, `0d4aee38`, `b73457e7`, `1d652941`, `22b2d1c5` and `84886850`;
- on `domain/humanities`: `1a110763`, `22e0bca3`, `a37de772`, `d1cdbbd7`, `b7d2004e`, `84f01464`, `93b2b25c`, `9140c792`, `ec3aca47`, `a1b2b4d2`, `fb690e64`, `4debbd75`, `71cbcbe3`, `6f3ab3e6`, `e2d077a8`, `65804012`, `487d60f4`, `d06172c0`, `2b7df218`, `a63cd71b` and `c077b606`;
- on main: `4d8931d3`, the check selector.

## The register's section of 12 to 14 November

The author's decisions in the Codex session that began with sport's economics were drafted as register rows on 15 November and had not been entered. The records enter them as a section of 28 rows before the section of 15 November, each Status naming the commit and the place where the work on it stands. Their 43 quotations are exact against the session's messages kept on the research drive (`check_register_quotes.py`, `vacuum-economy/intake/2026-11-12_codex/`).

## Completeness and pins

The completeness script of the fourth round finds every head merged: the four branches', and provisioning's and humanities' inside the economy's. It flags nine prose pointers and two pins. Eight pointers and both pins were on main before the round, in files it did not change: the economy's register rows and plan sections of 12 November, a run note in geology's README, and the pins in the optical-comfort and tile-design studies. The ninth is in the returns study's revised note, which names an export "on `domain/provisioning`" at `a8af8761`. The register's and the plan's sections of 12 November name the same branch for six more pieces of work. Those commits are now in main's history, and the studies they name are in the tree.

The 211 JSON products the round changed hold 1,594 SHA-256 values, 1,451 of them naming files. 1,344 match the files they name, 73 differ and 34 name files the check did not find (`pins.py` and `pins.log` on the research drive, `session-tools/integration9/`). The same check run on each branch's own head flags the same 107 and no other (`pins_geography_branch.log`, `pins_economy_branch.log`, `pins_resources_branch.log`), so the merges changed none of them. Of the 73 that differ, 39 are in nine of the atlas round's products: the districts' record, the dryland rivers, the first lake rule's measure and its balance, the works, the settled turn, the turned surface, the surface's climate and the moisture map. Each records its inputs and code as they stood when it was built, and later steps changed some of those; the round's README names the districts' record and the dryland rivers among what still stands on the earlier design. The other 34 are in the staged studies' input records and manifests of 12 to 14 November. The 34 not found are files kept outside Git, 15 responses from JPL's small-body database and 18 review records and logs of the settlement studies' checks, and one file that no longer exists: the mining study's sources name a retained register at `research/studies/vacuum_buildout/sources.json`, a folder that has held only its README since the proportional appraisal was replaced (`e2d077a8`).

## Outside Git

- **The economy's run folders** are on the research drive, `vacuum-economy/runs/`, linked from `research/runs/` as `complex_returns`, `enclosed_feasibility`, `enclosed_industry`, `enclosed_power`, `enclosed_returns`, `monetary_wealth`, `production_income` and `source_access`. The industry runs' cache of cases is the largest, 890 MB. A checkout without the links runs the studies' tests on the committed products and reruns nothing.
- **The atlas round's grids** are on the research drive, `atlas-round/`, named by absolute path in the products that read them, 161 recorded digests in all; the design climate's runs are in the climate domain's run folders, linked from `climate/gcm/runs`. The ocean models' runs and products (`climate/ocean/runs`, `climate/ocean/products`) and the south sea's run folder (`research/runs/south_sea_ocean`) are links to the drive that the geography and hydrosphere worktrees hold.
- **The naming material** of 19 and 20 November, the map of parent cultures and the four samples of names, is on the research drive, `session-tools/atlas_places/`, and enters the repository when names are adopted.
- **The economy's documents** as PDFs, with their earlier builds, and the lookups behind the inputs are on the research drive, `vacuum-economy/report/` and `vacuum-economy/research/`.
- **The built globe pages** are in the geography worktree's `visualization/atlas/out/`, ignored by Git.

## Checks

The complete suite ran on the merged tree at `e6759adb`, as the project's checks ask of an integration (`make check-full`, 18:10 to 18:35; the log is on the research drive, `session-tools/integration9/make_check_full.log`):

| Check | Result |
|---|---|
| Layers | 1,966 files, no violation |
| Python tests | 7,104 passed, 36 skipped |
| JavaScript tests, the calendar's 14 among them | 191 passed |
| The calendar's build, the provenance reproduction and the ensemble's validation | passed |

It passed in one run. Eleven links to the research drive were made in the worktree before it, as the branches' own worktrees hold them: the economy's eight run folders, the ocean models' runs and products, and the south sea's run folder. The settle pass found no product stale on the merged tree (`settle_merged_dry.log`). The register's section of 12 to 14 November was entered in the worktree two minutes into the run; the register's lines before it, which the source-access gates and the population test read, did not move. The records change records and prose only, and the checks they select ran after them and passed (`make_check_records.log`).

## The review

An independent reviewer reads the merged tree at the records commit, without changing it, before main moves (its brief on the research drive, `session-tools/integration9_review/`). Its findings and the answer to them are entered here by the answering commit.

## Still to do from the integration

- **The two studies** are closed for now, and each README lists what is open: [the atlas round](../studies/atlas_round/README.md#what-is-open), its names in the middle of a rework, and [the economy](../studies/vacuum_economy_periods/README.md#where-the-study-stands-21-november-2026), with trade within the zone designed and unbuilt.
- **Where the two touch.** The four tube cities' siting reads the atlas's waters at 28% (`geography/products/atlas_28pct_4ppd.npz`, in `research/studies/lava_tube_cities/siting_inputs.json`). The atlas round's design surface holds different waters: the common level 39 m higher, the three upper seas filled to their rims, and the lakes of the corrected rule. The search regions have not been read against it.
- **The heritage register** (`research/studies/conservation/results/heritage_register.json`) lists the Shoemaker scattering site and the Cabeus impact as under water; on the design surface both are dry bowls.
- **The economy's layered studies** of 12 to 14 November, and the register's and the plan's sections of 12 November, keep the prose they were written in; they have had no editorial pass.
- **The research index** (`research/README.md`) gains the atlas round's row with this round and still lacks rows for ten studies of earlier rounds: cold_waters, drainage_network, littoral_community, polar_light, river_deltas, seawater_lakes, south_sea_ocean, storm_water, underwater_light and vacuum_buildout.
- **The south sea's oxygen run** stays paused on `domain/hydrosphere` at model year 53; its basin is redrawn in the atlas round's first step, and no ocean run has been made on it.
- **The shared population scenario** (`shared/scenarios/population.json`) still lists comparison cases only, while the register holds the distribution of 11 October (C′: the Moon 8 billion, 4 on the ground and 4 aloft; the array 5; Earth 9; the colonies 3).
- **Carried from earlier rounds:** the soils study, the lakes' distribution among fresh and soda, and how the moisture and the terrain combine into the biome map; the year files' first humidity records still to be carried into four studies; the cold-waters study's forcing on the turned Moon; geology's fallback to the hydrosphere worktree in `geology/interior_response.py`; climate's stale code pins in the waves', the snow's and the cloud twilight's products; and the CM1 boxes on the turned years, held with the gap carried as a range.
