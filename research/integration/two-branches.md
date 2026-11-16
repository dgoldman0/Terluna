# The joint integration of 16 November 2026: hydrosphere and ecology, without the economy

Three states of the work became one on 16 November 2026: main at `ab8d278c`, ecology at `0dd5aa0b` and hydrosphere at `2a24c4f8`. The merges ran on the standing integration branch `maintenance/integration` in the worktree `../terluna-integration`, fast-forwarded first from the seventh round's head, `595694ff`, to main. Each merge is an ordinary merge commit by exact head, ecology and then hydrosphere; neither branch's work reads the other's.

The round leaves out the vacuum buildout's economy, as the register's section of 16 November sets it ([the atlas round's frame](../decisions.md#the-atlas-rounds-frame-16-november-2026)): `research/vacuum-economy`, and the economy's commits of 12 to 14 November on `domain/provisioning` (23 beyond main) and `domain/humanities` (21). `domain/resources` stays out whole. Six of its nine commits beyond main are the economy's source studies of 13 November, one merges main and one binds their gates to the register. The ninth, the salt for an ocean-like south sea (`346e17f0`, 15 November), reads the economy's source-access product, `resources/source_access/results/potential.json.gz` from `1d652941`, which main does not hold, so it comes in with the economy. The south sea study names that work as done "on `domain/resources` (`resources/sea_salt.md` there)", which reads correctly on main without it. Atmosphere, climate, food-systems, geology, array-industry and Jura Island stood at main's earlier head with no work of their own.

The round has no synthesis study and no joint computation, since the two branches' work meets in no product.

## The commits

| Commit | What it does |
|---|---|
| `7f18df4e` | Merges ecology: the polar light rebuilt on the corrected design climate, with the night side's twilight counted once, and the moisture map, the terrain axis and the littoral water map rebuilt after it; the seventh round's open item |
| `226b44d0` | Merges hydrosphere: the south sea as an ocean (its salt and water balance, what a salt sea does to its shores and air, its 50-year circulation with the monthly tide and a fine patch at the Moscoviense mouth, the Moon's own interior mixing and an oxygen run paused at model year 53), the drainage network as the basins fill, the globe's years control with the climate fields smoothed under its ground cover, and the atlas round's choices of 16 November carried into the studies |
| the records commit | This account, the register's rows of 16 November brought up to date, the plan and this folder's index |

No merge had a conflict, and no file was changed on both sides since the merge bases. The settle pass found no product stale after each merge (`settle_ecology_dry.log`, `settle_hydrosphere_dry.log`). No merged file reads another branch by a pinned commit; the atlas README's build command names the orientation product's source as `domain/geology@45c3d817`, a provenance string the builder records.

## Completeness and pins

The completeness script of the fourth round finds both heads merged. It flags five prose pointers and two pins, each on main unchanged at `ab8d278c` and in a file the round did not touch: the economy's register rows and plan sections of 12 November, which name commits on `domain/provisioning`, a run note in geology's README, and the pins in the optical-comfort and tile-design studies.

The 16 JSON products the round changed record 126 digests. 116 name files in the tree and match them; the other 10 name the corrected climate run's year files `MOST.00010` to `MOST.00019`, which match in main's run folder. None differs.

## Outside Git

- One merged product records an absolute path into a branch's worktree: the south sea's mixing inputs (`research/studies/south_sea_ocean/results/lunar_mixing_inputs_1deg.json`) name their maps in `../terluna-hydrosphere/research/runs/`. The paused oxygen run reads those inputs through its own setup, which must not change before it resumes, so the record stays as it is until the study next rebuilds its inputs.
- The drainage network's maps and the globe's year maps lie in the hydrosphere worktree's `research/runs/drainage_network/`, and the ocean runs in its `climate/ocean/runs/`, all ignored by Git. The products record the year maps by their repository path.

## Checks

The complete suite ran on the merged tree at `226b44d0`, as the project's checks ask of an integration (`make check-full`, 10:19 to 10:38; the log is on the research drive, `session-tools/integration8/make_check_full.log`):

| Check | Result |
|---|---|
| Layers | 1,676 files, no violation |
| Python tests | 5,540 passed, 36 skipped |
| JavaScript tests, the calendar's 14 among them | 191 passed |
| The calendar's build, the provenance reproduction and the ensemble's validation | passed |

It passed in one run, with no ignored input copied into the worktree. The south sea study's three regeneration tests, of its inflows, shores and spray, skip unless `TERLUNA_SOUTH_SEA_REGENERATE=1` is set, as they do in any checkout; they take 3.5, 1 and 2 minutes and up to 5 GB.

## The review

An independent reviewer reads the merged tree at the records commit before main moves; its findings and their answer follow here.

## Still to do from the integration

- **The vacuum buildout's economy**, with `domain/provisioning`, `domain/humanities` and `domain/resources`, the south sea's salt among it, in the round after its refactoring's report.
- **The atlas round** on `domain/geography`, its frame set on 16 November ([register](../decisions.md#the-atlas-rounds-frame-16-november-2026)).
- **The south sea's oxygen run** resumes on `domain/hydrosphere` from model year 53, about 13–18 hours for its remaining 47 model years.
- **The shared population scenario** (`shared/scenarios/population.json`) still lists comparison cases only, while the register holds the author's distribution of 11 October (C′: the Moon 8 billion, 4 on the ground and 4 aloft; the array 5; Earth 9; the colonies 3).
- **Carried from earlier rounds:** the soils study and the lakes' distribution, now fresh and soda only (register, 16 November); the year files' first humidity records still to be carried into four studies; the cold-waters study's forcing on the turned Moon; geology's fallback to the hydrosphere worktree in `geology/interior_response.py`; climate's stale code pins in the waves', the snow's and the cloud twilight's products; and the CM1 boxes on the turned years, held with the gap carried as a range.
