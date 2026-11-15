# The joint integration of 15 November 2026: three branches on the corrected design climate

Four states of the work became one on 15 November 2026: main at `4059f8bb`, climate at `39c5ee04`, ecology at `d4b4a37b` and humanities at `91b1a5ac`. The merges ran on the standing integration branch `maintenance/integration` in the worktree `../terluna-integration`, fast-forwarded first from the sixth round's head to main. Each merge is an ordinary merge commit by exact head, in the order the work depends on itself: climate, whose corrected run the ecology studies read, then ecology, then humanities.

The round takes the work of 8 to 10 November, finished before the vacuum buildout's economy was opened on 12 November. That economy, the settlements, industry, sources and returns studied on `domain/provisioning`, `domain/resources` and `domain/humanities` beyond `91b1a5ac`, is to be refactored on a branch of its own, `research/vacuum-economy`, which is yet to be made, and comes in at the next round ([register](../decisions.md#the-seventh-integration-and-the-refactoring-of-the-vacuum-buildouts-economy-15-november-2026)). Main already held nine commits from those sessions when the round began: the register's sections of 12 November, the plan's settlement sections, `habitation/crater_gardens` and the selection of checks by what a change affects. They are part of the round's base and are reviewed with the economy. Atmosphere, food-systems, geology, hydrosphere, array-industry and Jura Island stood at the sixth round's head with no work of their own.

The round has no synthesis study and no joint computation. The products in which these branches meet were built on the ecology branch from 8 to 10 November, from the climate branch's runs and committed products, and the merged tree holds them beside their inputs.

## The commits

| Commit | What it does |
|---|---|
| `ee633090` | Merges climate: the night side's twilight counted once in PlaSim's shortwave, the equinox regressing with the orbit's node, the corrected design run on the turned Moon over forty years with the gap to CM1 carried as a range, the turned-years gap boxes, the weather cycles' rescreen, the globe and the drainage on the corrected climate, and year 10 written every 6.5 hours |
| `7393945f` | Merges ecology: the moisture map and the terrain sheets on the corrected climate's forty years, the polar light of the turned Moon, the terrain tried beside the moisture on the exact data, and the freshwater littoral study (the shore's water, the founders' review, the waters' warmth map, the wind inside a 3-day mean, the shelved simulator) |
| `f58a1e22` | Merges humanities at `91b1a5ac`: the sport study's second round, aimed projectiles and moves off obstacles, the landing loads, the hanging salto's landing and the motion engine's memory and integrator |
| the records commit | This account, the register's section of 15 November, the plan and this folder's index |

No merge had a conflict. Humanities had added six rows of 8 and 9 November to the register and updated two earlier ones, for the sport and the dance; Git's merge set them in the two sections where the branch had placed them, with every line of both sides kept once and in order. Climate and ecology had each changed `geography/.gitignore` and the corrected drainage record, to the same content. No other file was changed on two sides. The settle pass found no product stale after each merge and again once the two drainage grids were in the worktree, and no file in the merged work reads another branch by a pinned commit.

## Where the branches meet

The corrected design climate is the base of the ecology work rebuilt on 9 and 10 November. The moisture map and the terrain sheets read the corrected run's forty years and the drainage routed on them; the littoral forcing and the waters' warmth map read its years 10 to 19; and the littoral wind record reads its year written every 6.5 hours, beside two years of an earlier run written the same way.

The polar light study was built on 8 November, before the corrected run existed, and still reads the turned run at sunlight 0.924, in which the night side's twilight was counted twice. Its cloud cover, its rain-fed lakes and its sunlight diagnostic come from that run; the Sun's geometry, which sets its clear-sky light, is the same for both runs. The moisture map takes its plant darkness poleward of 60° from that study, so those dark hours rest on the clear-sky geometry, with the earlier run's clouds in the cloudy column beside it. Whether the polar light is rebuilt on the corrected run is open (below).

These nine products record 129 input digests. The 94 that name a file within reach all match it: 63 tracked files, 16 ignored or linked ones, the two drainage grids among them, 13 on the research drive and 2 in main's checkout. Of the 26 that name the climate runs' year files, the reviewer hashed three, which match. The other nine are recorded without a path. The humanities work reads the design air and no product of the other two branches.

## Outside Git

- The drainage grids routed on the corrected climate, `geography/turned_corrected/products/drainage_28pct_16ppd.npz` and `geography/turned_corrected_y10-49/products/drainage_28pct_16ppd.npz`, are ignored by Git and lay only in the climate and ecology worktrees. The moisture map and the terrain sheets pin them, and their tests read them. Both are copied into the integration worktree, matching the recorded digests, and go into main's checkout and every worktree when the branches come up.
- The corrected run's climatologies lie in main's products folder, which every worktree links.
- Two kinds of absolute path are recorded in the merged products, neither pinning a file inside a branch's worktree: the corrected drainage records name their climatology in main's checkout, as the earlier drainage records do, and the weather cycles' record names the climate worktree as the folder and command of its run.

## Checks

The complete suite ran on the merged tree at `f58a1e22`, as the project's checks ask of an integration (`make check-full`; the logs are on the research drive, `session-tools/integration7/`):

| Check | Result |
|---|---|
| Layers | 1,658 files, no violation |
| Python tests, first run under `make`, before the two drainage grids were in the worktree | 5,450 passed, 34 skipped, 2 failed: the moisture map's two tests that read the grids |
| Python tests, second run, with the grids, started outside `make` | 5,452 passed, 33 skipped, 1 failed: the sky ecology's exact-regeneration test, which needs the single numerical thread that `make` sets and had passed under it in the first run |
| JavaScript tests, the calendar's 14 among them | 191 passed |
| The calendar's build, the provenance reproduction and the ensemble's validation | passed |
| The records commit and the answering commit | the one test selected for each passed |

Each of the 5,453 tests that run in this worktree passed in the configuration it needs, over the two runs together; neither run alone was clean. The second run came from a test selection that turned out empty and so ran everything; it stands as the record of the suite with the grids in place. The moisture map's 11 tests and the sky ecology's 18 were also rerun on their own under `make` once the grids were copied in, and passed; those two reruns were not logged. One clean run of the complete suite follows in main's checkout as soon as main moves up, before any branch is brought up to it.

## The review

An independent reviewer read the merged tree at the records commit, `b9285b9c`, without changing it, before main moved (its brief and report on the research drive, `session-tools/integration7_review/`). It found nothing blocking, seven points to fix in the records and eight minor ones. It confirmed that the tree is main plus the three branches' changes exactly, the 97 changed paths equal to the difference between main and the merged tree; that Git's merge of the register reproduces byte for byte, with every line of both sides once and in order; that nothing of the vacuum buildout's economy entered; that every recorded digest naming a file within reach matches it; that the test figures agree with the logs; and that the register's quotations are the author's exact words. The answering commit:
- **The polar light.** The account said the polar light study reads the corrected run. It reads the turned run before the correction, as set out above, and its rebuild is listed below as open.
- **The checks in the register.** The row said the complete suite had passed; it now says what the two runs showed, and the clean run follows in main's checkout.
- **Files changed on two branches.** Climate and ecology changed two files to the same content, which the account now says, with where the humanities rows sit in the register.
- **The pins' count** is given in full, 129 digests in nine products, in place of the 85 the lead had checked.
- **The sixth round's leftovers**, geology's fallback and climate's stale code pins, are carried below, with the state of the CM1 boxes on the turned years.
- **The attribution row** gives what its counts are, the scripts' place on the drive and the condition attached to the decision, and is set at the merge.
- **Smaller points.** The branch yet to be made is described as such; the absolute paths recorded in three products are noted; the JavaScript count is said to include the calendar's tests; the settle pass was repeated with the grids in place; and the agents row states only what was accepted.

## Still to do from the integration

- **Main moves up** to the round's head now that the review is answered, with one clean run of the complete suite in its checkout ([register](../decisions.md#the-seventh-integration-and-the-refactoring-of-the-vacuum-buildouts-economy-15-november-2026)). Climate, ecology and the six branches without work of their own then come up to main by fast-forward, each worktree given both corrected drainage grids. Provisioning, resources and humanities keep their heads until the economy's round.
- **The polar light on the corrected run.** The study's clouds, lakes and sunlight diagnostic come from the turned run before the correction, and the moisture map reads its polar darkness. Rebuilding it on the corrected run takes about twelve minutes on four workers, followed by the moisture map's and the terrain sheets' rebuilds; whether to rebuild or to state the limit is open.
- **The vacuum buildout's economy**, refactored on `research/vacuum-economy` and merged at the eighth round, with the nine commits already on main reviewed alongside.
- **The soils study and the lakes' distribution**, and how the moisture and the terrain combine into the biome map, carried over from the sixth round.
- **The year files' first humidity records**, whose corrected reading in the littoral forcing and the map is still to be carried into four other studies ([register](../decisions.md), 10 November).
- **The cold-waters study's forcing**, redone with the ocean studies on the turned Moon ([register](../decisions.md), 10 November).
- **Geology's fallback** to hydrosphere's worktree in `geology/interior_response.py`, carried from the sixth round and dropped at geology's next rebuild.
- **Climate's stale code pins** in the waves', the snow's and the cloud twilight's products, carried from the sixth round and rebound or rebuilt at their next build.
- **The CM1 boxes on the turned years**, the sixth round's third leftover, arrived with climate: the sea box ran, and the land and dry boxes are held with the gap carried as a range ([register](../decisions.md), 8 November).
