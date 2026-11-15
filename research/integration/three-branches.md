# The joint integration of 15 November 2026: three branches on the corrected design climate

Four states of the work became one on 15 November 2026: main at `4059f8bb`, climate at `39c5ee04`, ecology at `d4b4a37b` and humanities at `91b1a5ac`. The merges ran on the standing integration branch `maintenance/integration` in the worktree `../terluna-integration`, fast-forwarded first from the sixth round's head to main. Each merge is an ordinary merge commit by exact head, in the order the work depends on itself: climate, whose corrected run the ecology studies read, then ecology, then humanities.

The round takes the work of 8 to 10 November, finished before the vacuum buildout's economy was opened on 12 November. That economy, the settlements, industry, sources and returns studied on `domain/provisioning`, `domain/resources` and `domain/humanities` beyond `91b1a5ac`, is refactored on a branch of its own, `research/vacuum-economy`, and comes in at the next round ([register](../decisions.md#the-seventh-integration-and-the-refactoring-of-the-vacuum-buildouts-economy-15-november-2026)). Main already held nine commits from those sessions when the round began: the register's sections of 12 November, the plan's settlement sections, `habitation/crater_gardens` and the selection of checks by what a change affects. They are part of the round's base and are reviewed with the economy. Atmosphere, food-systems, geology, hydrosphere, array-industry and Jura Island stood at the sixth round's head with no work of their own.

The round has no synthesis study and no joint computation. The products in which these branches meet were built on the ecology branch from 8 to 10 November, from the climate branch's runs and committed products, and the merged tree holds them beside their inputs.

## The commits

| Commit | What it does |
|---|---|
| `ee633090` | Merges climate: the night side's twilight counted once in PlaSim's shortwave, the equinox regressing with the orbit's node, the corrected design run on the turned Moon over forty years with the gap to CM1 carried as a range, the turned-years gap boxes, the weather cycles' rescreen, the globe and the drainage on the corrected climate, and year 10 written every 6.5 hours |
| `7393945f` | Merges ecology: the moisture map and the terrain sheets on the corrected climate's forty years, the polar light of the turned Moon, the terrain tried beside the moisture on the exact data, and the freshwater littoral study (the shore's water, the founders' review, the waters' warmth map, the wind inside a 3-day mean, the shelved simulator) |
| `f58a1e22` | Merges humanities at `91b1a5ac`: the sport study's second round, aimed projectiles and moves off obstacles, the landing loads, the hanging salto's landing and the motion engine's memory and integrator |
| the records commit | This account, the register's section of 15 November, the plan and this folder's index |

No merge had a conflict. Git joined the register's rows for the sport work of 8 and 9 November into their dated sections, and no other file was changed on both sides. The settle pass on each merged tree found no product stale, and no file in the merged work reads another branch by a pinned commit.

## Where the branches meet

The corrected design climate is the base of everything ecology merged. The moisture map, the terrain sheets and the polar light read the corrected run's forty years and the drainage routed on them; the littoral forcing and the waters' warmth map read its years 10 to 19, and the littoral wind record reads the year written every 6.5 hours. Each of these products records the digests of the files it read: the 85 that name files in the tree match them, and the others name the corrected run's year files on the research drive. The humanities work reads the design air and no product of the other two branches.

## Outside Git

- The drainage grids routed on the corrected climate, `geography/turned_corrected/products/drainage_28pct_16ppd.npz` and `geography/turned_corrected_y10-49/products/drainage_28pct_16ppd.npz`, are ignored by Git and lay only in the climate and ecology worktrees. The moisture map and the terrain sheets pin them, and their tests read them. Both are copied into the integration worktree, matching the recorded digests, and go into main's checkout and every worktree when the branches come up.
- The corrected run's climatologies lie in main's products folder, which every worktree links.

## Checks

The complete suite ran on the merged tree at `f58a1e22`, as the project's checks ask of an integration (`make check-full`; the logs are on the research drive, `session-tools/integration7/`):

| Check | Result |
|---|---|
| Layers | 1,658 files, no violation |
| Python tests, first run, before the two drainage grids were in the worktree | 5,450 passed, 34 skipped, 2 failed: the moisture map's two tests that read the grids |
| The moisture map's tests with the grids copied in | 11 passed |
| Python tests, second run, with the grids, started outside `make` | 5,452 passed, 33 skipped, 1 failed: the sky ecology's exact-regeneration test, which needs the single numerical thread that `make` sets |
| The sky ecology's tests under `make` | 18 passed |
| JavaScript tests | 191 passed |
| The calendar's build and tests, the provenance reproduction and the ensemble's validation | passed |

Each of the 5,453 tests that run in this worktree passed in the configuration it needs, over the two runs together; neither run alone was clean. The second run came from a test selection that turned out empty and so ran everything; it stands as the record of the suite with the grids in place. A single clean run follows in main's checkout when main moves up. The settle pass found nothing stale after each merge, and the 85 input digests that the merged ecology and drainage products record for files in the tree all match those files. The records commit changes prose alone and passed the checks selected for it.

## The review

An independent reviewer reads the merged tree at the records commit before main moves (its brief and findings on the research drive, `session-tools/integration7_review/`). Its findings and the answering commit are recorded here.

## Still to do from the integration

- **Main moves up** to the round's head at the author's word. Climate, ecology and the six branches without work of their own then come up to main by fast-forward, each worktree given both corrected drainage grids. Provisioning, resources and humanities keep their heads until the economy's round.
- **The vacuum buildout's economy**, refactored on `research/vacuum-economy` and merged at the eighth round, with the nine commits already on main reviewed alongside.
- **The soils study and the lakes' distribution**, and how the moisture and the terrain combine into the biome map, carried over from the sixth round.
- **The year files' first humidity records**, whose corrected reading in the littoral forcing and the map is still to be carried into four other studies ([register](../decisions.md), 10 November).
- **The cold-waters study's forcing**, redone with the ocean studies on the turned Moon ([register](../decisions.md), 10 November).
