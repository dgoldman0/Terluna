# The joint integration of 16 November 2026: hydrosphere and ecology, without the economy

Three states of the work became one on 16 November 2026: main, at `ab8d278c` when the round began and at `296286a1` by its end, ecology at `0dd5aa0b` and hydrosphere at `2a24c4f8`. The merges ran on the standing integration branch `maintenance/integration` in the worktree `../terluna-integration`, fast-forwarded first from the seventh round's head, `595694ff`, to main. Each branch's merge is an ordinary merge commit by exact head, ecology and then hydrosphere; neither branch's work reads the other's. Main gained one commit while the round ran, a register row of the vacuum buildout's economy, and was merged in before the round's review was answered.

The round leaves out the vacuum buildout's economy, as the register's section of 16 November sets it ([the atlas round's frame](../decisions.md#the-atlas-rounds-frame-16-november-2026)): `research/vacuum-economy`, and the economy's commits of 12 to 14 November on `domain/provisioning` (23 beyond main) and `domain/humanities` (21). Eight research branches made during the economy's returns work, consolidated-returns, credit-supported-expansion, export-strategies, price-pressure, staged-complex-returns, u-monetary-pressure, u-policy-risk and u-trade-transition, hold only commits inside `domain/humanities` and `research/vacuum-economy`, and stay out with them. `domain/resources` stays out whole. Six of its nine commits beyond main are the economy's source studies of 13 November, one merges main and one binds their gates to the register. The ninth, the salt for an ocean-like south sea (`346e17f0`, 15 November), reads the economy's source-access product, `resources/source_access/results/potential.json.gz` from `1d652941`, which main does not hold, so it comes in with the economy. The south sea study names that work as done "on `domain/resources` (`resources/sea_salt.md` there)", which reads correctly on main without it. Atmosphere, climate, food-systems, geology, array-industry and Jura Island stood at main's earlier head with no work of their own.

The register's row on the round had been written that morning as taking in "every branch with work beyond main except the vacuum economy's". Three of those branches carry the economy, so the records commit corrected its Decision to the work beyond main except the economy's, and its Status gives the salt study's place as the round's own judgement of what depends on the economy.

The round has no synthesis study and no joint computation, since the two branches' work meets in no product.

## The commits

| Commit | What it does |
|---|---|
| `7f18df4e` | Merges ecology: the polar light rebuilt on the corrected design climate, with the night side's twilight counted once, and the moisture map, the terrain axis and the littoral water map rebuilt after it; the seventh round's open item |
| `226b44d0` | Merges hydrosphere: the south sea as an ocean (its salt and water balance, what a salt sea does to its shores and air, its 50-year circulation with the monthly tide and a fine patch at the Moscoviense mouth, the Moon's own interior mixing and an oxygen run paused at model year 53), the drainage network as the basins fill, the globe's years control with the climate fields smoothed under its ground cover, and the atlas round's choices of 16 November carried into the studies |
| `a95cb729` | The records: this account, the register's rows of 16 November on the round brought up to date, the plan and this folder's index |
| `9216b3cf` | Merges main at `296286a1`, the register row on the economy's freight that main gained at 10:40 during the round |
| the answering commits | Answer the independent review; the second corrects three phrases the first introduced, which the reviewer found on checking the answer |

No merge had a conflict. The branches' merges touched no file main had changed, and in main's merge only the register changed on both sides: main's row in the economy's section and the records commit's rows in the round's. The settle pass found no product stale after each branch's merge (`settle_ecology_dry.log`, `settle_hydrosphere_dry.log`). No merged file reads another branch by a pinned commit; the atlas README's build command names the orientation product's source as `domain/geology@45c3d817`, a provenance string the builder records.

## Completeness and pins

The completeness script of the fourth round finds both heads merged. It flags five prose pointers and two pins, each on main unchanged at `ab8d278c` and in a file the round did not touch: the economy's register rows and plan sections of 12 November, which name commits on `domain/provisioning`, a run note in geology's README, and the pins in the optical-comfort and tile-design studies.

The 16 JSON products the round changed hold 144 SHA-256 values. 133 name files: 123 by their key, 9 by a sibling field and one, the eddies run's 2.4 GB year file, by its key with its size. Every one matches its file: 114 in the tree, 9 on the research drive by absolute path and the corrected climate run's 10 year files, `MOST.00010` to `MOST.00019`, in main's run folder (`pins.py` and `pins.log` on the research drive, `session-tools/integration8/`, the second version after the review). The other 11 are hashes of code or of identities: the year maps' and the mixing inputs' code, the polar light's frames and heights, and one compare check of the moisture map's.

## Outside Git

- The south sea's mixing inputs (`research/studies/south_sea_ocean/results/lunar_mixing_inputs_1deg.json`) record two things that the study's next rebuild of its inputs rebinds, after the paused oxygen run, whose setup they feed. They name their maps by one absolute path into the hydrosphere worktree, `/home/kir/Documents/Projects/terluna-hydrosphere/research/runs/south_sea_ocean/lunar_mixing_inputs_1deg.npz`. And their producer's code hash, `d91eab2b…`, is given by no committed version of `lunar_mixing.py` and `climate/ocean/stratified.py`, which give `20fa9d99…`. The inputs were computed at 06:23 on 16 November, and three edits followed before the commit: the GGL90 namelist's file names and the tracers' size header in `stratified.py`, both for the one-year trial, and the command line of `lunar_mixing.py`. None of them is in the functions that compute the inputs, and the oxygen run began at 06:47, after the last of them; the edits and their times come from the session's own record of its edits, with which the run folders' times and namelists agree. No settle registry names this product.
- The drainage network's maps and the globe's year maps lie in the hydrosphere worktree's `research/runs/drainage_network/`, and the ocean runs in its `climate/ocean/runs/`, all ignored by Git. The year maps' product names its maps by their repository path. A globe build from another checkout needs `globe_stops.npz` copied into its `research/runs/drainage_network/`, or `stops.py` rerun, which the product records as taking 551 seconds.

## Checks

The complete suite ran on the merged tree at `226b44d0`, as the project's checks ask of an integration (`make check-full`, 10:19 to 10:38; the log is on the research drive, `session-tools/integration8/make_check_full.log`):

| Check | Result |
|---|---|
| Layers | 1,676 files, no violation |
| Python tests | 5,540 passed, 36 skipped |
| JavaScript tests, the calendar's 14 among them | 191 passed |
| The calendar's build, the provenance reproduction and the ensemble's validation | passed |

It passed in one run, with no ignored input copied into the worktree. The south sea study's three regeneration tests, of its inflows, shores and spray, skip unless `TERLUNA_SOUTH_SEA_REGENERATE=1` is set, as they do in any checkout; they take 3.5, 1 and 2 minutes and up to 5 GB. The records, main's merge and the answering commits change records and prose only, and the checks they select passed.

## The review

An independent reviewer read the merged tree at the records commit, `a95cb729`, without changing it, before main moved (its brief, findings and progress on the research drive, `session-tools/integration8_review/`). It confirmed that the two merges are ordinary merges by exact head, each changing exactly its branch's 8 and 39 files, every one equal to the branch's head; that nothing of the economy entered, from any of its branches; that the salt study's exclusion holds as stated; that every test figure agrees with the log; that the register's quotations are the author's exact words, one of them attached to the wrong reply; and, in its sample, the drainage network's table and the globe's version 10 against their products and code. Its one blocking finding was main's move during the round, which left main unable to fast-forward to the records; main is merged in at `9216b3cf`. The answer, which the reviewer then checked finding by finding:

- **The register.** One row's Source gave the reply to another recommendation and now gives the reply to its own; the upkeep budget's row keeps the bound "a little"; the row on the seas' levels speaks of active transport, as its Source does; and the row on the round names its corrected Decision and gives the economy's branch without a commit count, which moved during the round.
- **This account.** Main's third state and its merge; the eight research branches beside the economy; the pin check over all 144 values, with the year files hashed in the logged run; the mixing inputs' code hash beside their absolute path; and the year maps a globe build elsewhere needs.
- **What remains.** The economy's nine commits already on main, reviewed alongside it, and the biome map's axes return to the lists here and in the plan, where the records had dropped them.
- **Smaller points.** The globe's year-0 check is described as the build makes it, a count of pixels allowed up to 0.05% of the grid and zero in this build; the plan's item on the atlas round adds the heritage review and the seas' 28% with the upper seas on top; and the south sea README names the oxygen run's product without linking it before it exists, and says which edits followed its mixing inputs.

## Still to do from the integration

- **The vacuum buildout's economy**, with `domain/provisioning`, `domain/humanities` and `domain/resources`, the south sea's salt among them, in the round after its refactoring's report, with the nine commits already on main reviewed alongside.
- **The atlas round** on `domain/geography`, its frame set on 16 November ([register](../decisions.md#the-atlas-rounds-frame-16-november-2026)).
- **The south sea's oxygen run** resumes on `domain/hydrosphere` from model year 53, about 13–18 hours for its remaining 47 model years.
- **The shared population scenario** (`shared/scenarios/population.json`) still lists comparison cases only, while the register holds the distribution of 11 October (C′: the Moon 8 billion, 4 on the ground and 4 aloft; the array 5; Earth 9; the colonies 3).
- **Carried from earlier rounds:** the soils study, the lakes' distribution, now fresh and soda only (register, 16 November), and how the moisture and the terrain combine into the biome map; the year files' first humidity records still to be carried into four studies; the cold-waters study's forcing on the turned Moon; geology's fallback to the hydrosphere worktree in `geology/interior_response.py`; climate's stale code pins in the waves', the snow's and the cloud twilight's products; and the CM1 boxes on the turned years, held with the gap carried as a range.
