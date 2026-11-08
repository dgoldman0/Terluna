# The joint integration of 7–8 November 2026: eight branches on the turned Moon

Nine states of the work became one on 7 and 8 November 2026: main at `a6ade546`, geology at `bcc7e753`, climate at `5145d80f`, hydrosphere at `f6fe4026`, ecology at `e95bdfcf`, provisioning at `afe38eda`, array-industry at `1b2358a1`, humanities at `776473f2` and food-systems at `ab9854c5`, after its rework. Atmosphere, resources and Jura Island stood at main's state with no work of their own. The merges ran on the standing integration branch `maintenance/integration` in the worktree `../terluna-integration`, its first round: the branch was fast-forwarded to main, which had just taken the round's first record, the ten rows of 7 November ([register](../decisions.md#the-standing-integration-branch-the-turned-runs-forty-years-the-biome-maps-axes-and-the-food-systems-rework-7-november-2026)). Each merge is an ordinary merge commit by exact head, in the order the work depends on itself: geology, whose orientation frames the climate's turned run, then climate, hydrosphere, the three small branches together, humanities and food-systems.

The round's joint computation is the moisture map, the first biome map renamed and remade on the turned run at the target over forty years and compared with years 10–19 and with the four decades classed alone (below).

## The commits

| Commit | What it does |
|---|---|
| `a6ade546` (main) | The round's first record: the standing branch, the turned run extended to fifty years, the biome map's axes and answers, and the food-systems rework begun |
| `2aa591fa`, `ab9854c5` (food-systems) | Main merged into the food-systems branch, then the branch's rework on its review by two bounded agents (below) |
| `327f8274` | Merges geology: the geology lane, the present Moon's viscoelastic response and the turned Moon's frame |
| `c8f4a462` | Merges climate: the CM1 boxes on the design years, the convection's refinements, the hollows' store in PlaSim, the turned design run at the target and the globe on the turned Moon |
| `6c753ee6` | Merges hydrosphere: the storm water, the strait and the tides, the cold waters, the water's own gravity, the underwater light, the deltas and the seawater lakes |
| `612d073f` | Merges ecology, provisioning and array-industry in one merge: the first biome map, the energy options and the fault across the ring |
| `bf0ee2fe` | Merges humanities: dance at lunar gravity, its garments and renderer, and the sport study from the body |
| `bf55f94d` | Merges food-systems after its rework, with the integration's follow-through |
| `40e28efa` | The moisture map on the turned run at the target over forty years, with its windows |
| the records commit | This account, the register's statuses, the plan and this folder's index |

The three small branches, of one commit or two each and touching separate files, joined in one merge so that the night's test runs stayed few; the earlier rounds merged each branch alone.

## How the conflicts were resolved

Under the register's old pin, every register edit on a branch rebuilt the products that bound the register and rebound their readers, so hydrosphere and humanities each carried rebinds of the same files with different digests. Hydrosphere's merge met 108 conflicts and humanities' 115. Each conflicted file was compared at the merge base, on the integration's side and on the branch's, with every hex token of 12 characters or more masked (`session-tools/integration6/resolve_hash_conflicts.py` on the research drive): where one side had changed only digests, the other side was kept. Hydrosphere kept the integration's side in 105 files and its own in three, the lake census, the lake kinds and the lake methane, whose values it had changed; humanities kept the integration's side in 110, among them the seas' phosphorus code, which on main no longer binds the register. Settle then ran on each merged tree: after hydrosphere it rebuilt 13 products and rebound 189 files, the three lake products differing from hydrosphere's only in their run time; after humanities it found nothing stale.

Some products bind inputs that settle does not follow. Where such a product's pinned input had been rebound or changed by the integration, the input was checked to hold the values the branch read, with digests and run times masked, and its digest was then followed through every file to a fixed point (`round_1107_pin/rebind_amenities.py` and `integration6/rebind_stale_pins.py`; the logs `rebinds_*.log`):
- the amenities study to the rebuilt sea-phosphorus product, of which it reads only the soil phosphorus grade, unchanged;
- the coast tide and the seawater lakes to the integration's sea chemistry, coast-water and climate-target products;
- the storm-water products and the biome map to climate's CM1 site code, of which they call only `subsolar_longitude`, unchanged;
- the energy options and the ring fault to the integration's power product, the domain synthesis's ledgers and the heat levers;
- the eel calculator to the rebound food-systems product;
- the cold waters, and the river deltas that read them, to the integration's lake oxygen, which the independent review found.

Five files were joined by hand. Geology and humanities each added a lane, and the layer check, the shared README, AGENTS.md and pytest's paths take both. Humanities had kept its decisions in its own copy of the register, in four sections dated 1, 2–3, 6 and 7 November; they stand among main's sections at their dates, every line of both kept (`integration6/join_register_humanities.py`).

## The food-systems rework

The branch's Codex work of 4–5 November was reworked on its review of 6 November before it merged, in the worktree `../terluna-food-systems`, by two agents with fixed files, bounds and about two hours each (briefs and reports in `session-tools/agents/food_systems_rework/`). The aquatic study stands on the turned design years; the branch's ideation register became ideas, labelled and set beside what the register holds; the eel chapter leads with a migration through a held seawater lake, with per-person cases (fed eel at 1 kg a person for 8 billion people leave 167 kt of phosphorus a year unretained, three to eight times what the waters' life can spare, so fed eel stay in closed systems); and the sturgeon chapter is rebuilt on the cold-waters study, with reproduction through hatcheries. The lead read the study's frame chapter in full and the rest in part, with a mechanical scan of the writing standards and checks of the figures against their sources.

## The follow-through

- The seawater lakes' README gives the eel chapter's new ranking, and its source register keeps the first state it read.
- The seawater sites read the turned surface in this tree by repository path; their product's values are unchanged.
- Geology's interior response keeps its fallback to hydrosphere's worktree, read only when the in-tree crustal-loading product is missing: editing that code would change the orientation product's digest, which names the climate's turned surface. The fallback goes at geology's next rebuild.
- Climate's commit `12673b97` let four CM1 readers (`wave_forcing.py`, `ring_analysis.py`, `gcm_profiles.py`, `cloud_scene.py`) read output kept raw or compressed, bit for bit, without rebinding the products that pin them as producer files: the waves' forcing, ring coverage and coastal verification, the snow's moist history and the cloud twilight's scenes. Their values stand; the pins are rebound or the products rebuilt at their next build, since a rebind now would also rewrite dated check snapshots.
- An output option for the climatology tool was tried and withdrawn: the zonal winds product pins the tool, and the rebind spread through 126 files. The forty-year climatology was written by copying the ten-year product aside, running the tool and restoring the canonical file, checked identical.

## The moisture map

The first biome map classed the design run of 1 November at today's orientation over ten years. Renamed the moisture map, since it classes the land by its water balance alone, it was remade on the turned run at the target over years 10–49 ([ecology README](../../biosphere/ecology/README.md#the-moisture-map-4-and-8-november-2026), [moisture_map.py](../../biosphere/ecology/moisture_map.py)). Its cells and night are the turned Moon's; geography's atlas, tides and lakes are counted in the cell that holds each pixel on the turned Moon through geology's orientation product; the lakes come from the drainage routed on the same forty years' climatology (`geography/turned_y10-49`).

| Class | Years 10–49 | Robust | Years 10–19 |
|---|---:|---:|---:|
| Hyper-arid | 1.9% | 1.0% | 2.1% |
| Arid | 21.1% | 12.7% | 23.0% |
| Semi-arid | 26.9% | 9.3% | 27.6% |
| Dry sub-humid | 8.0% | 0.0% | 6.8% |
| Humid and wet | 42.1% | 23.0% | 40.5% |

The drylands are 57.9% of the land over the forty years and 59.5% over years 10–19. Forty years leave 46% of the land robust across the gap's scale and the years' band together, against 37% for the turned run's own years 10–19 and 36% in the first map. Classed alone, two independent decades of the one run differ in class on 30.5% of the land (28–34%), or 25.8% in UNEP's five classes (24–29%): that is how far a ten-year map moves by chance. The 45% of ground whose five-class aridity changed between the first map and the turned run's ten years therefore exceeds chance by about 19 points, which bounds from above what the turn and the new climate changed, since matching each cell across the two frames adds to it; the review's rough estimate puts that matching near 14 points, a bound and not a measure. Years 10–19 keep the forty years' class on 75.8% of the land, and their rain-fed lakes on land cover 2.458 million km² against the forty years' 2.457 million; the forty-year drainage holds 33,686 lakes over 10.30% of the Moon against 33,262 over 10.28%. The turned design years' products, the globe and its rivers and lakes among them, therefore stand, and the forty-year versions confirm them. The 10–19 shares reproduce the turned years' land readings to within 3×10⁻¹⁷. The gap's scale and the night ratio still come from the CM1 boxes of the design years of 1 November, and the hollows' store runs, which belong to that run's land, do not enter.

## Outside Git

- The forty-year climatology and its companion on today's ground (`climate/gcm/products/climatology_A28R_301k_design_l15_snow_s924_turned_refit_y10-49.npz` and `_y10-49_today.npz`) and the ten-year product's named copy (`_y10-19.npz`) lie in main's products folder, which every worktree links.
- The drainage grids of both windows (`geography/turned/products/` and `geography/turned_y10-49/products/`), ignored by Git, are copied into main's checkout when it moves up, the ten-year one byte-identical to the grid the globe pinned.
- The moisture map's year sums are cached on the research drive (`session-tools/moisture_map/cache`), and its sheet is drawn there.
- The moisture map's tests bind both windows' drainage grids, so `make check` fails in main's checkout until they are copied in.

## Checks

Each merge and each commit of the round passed `make check` before it was committed:

| Commit | Python tests |
|---|---|
| Main's first record | 4,636 passed, 22 skipped |
| Main into food-systems; the rework | 4,656 passed, 24 skipped; 4,666 passed, 24 skipped |
| Geology | 4,657 passed, 24 skipped |
| Climate | 4,861 passed, 27 skipped |
| Hydrosphere, after the rebinds | 5,048 passed, 30 skipped |
| The three small branches, after the rebinds | 5,089 passed, 31 skipped |
| Humanities | 5,177 passed, 33 skipped |
| Food-systems with the follow-through | 5,209 passed, 33 skipped |
| The moisture map | 5,210 passed, 33 skipped |
| The records | 5,210 passed, 33 skipped |
| The answering commit | 5,210 passed, 33 skipped |

From 21:00 the machine rested thirty minutes between full test runs.

## The review

An independent reviewer read the merged tree at `d291b2f2` without changing it, before main moved (its findings on the research drive, `session-tools/integration6_review/findings.md`). It found nothing blocking, three points to fix and eleven minor ones. It confirmed that every branch's work is in the tree, that both merges' conflicts reproduce and hold under a stricter mask, that every rebind's input holds the values its reader read, that the hand joins keep both lanes and every register line, and that the moisture map's frame, lakes, windows and decades recompute from its inputs. The answering commit:
- **The cold waters' pin** to the lake oxygen, which settle had rebound in the hydrosphere merge, follows it with its values unchanged, and the river deltas after it.
- **The river deltas' inputs.** The deltas recorded their inputs by absolute path in hydrosphere's worktree, and their test read them there. They now record each input by its path in the repository; all eleven in the tree match their recorded digests, and the product's values are unchanged.
- **The 45%.** The decades show that two ten-year maps of one climate differ on about 26% of the land in five classes, so the 45% exceeds chance by about 19 points, an upper bound on the turn's and the climate's share, since matching cells across the two frames adds to it. The account, the register and the ecology README say so.
- **The compare window's lakes** are recorded in the product with their digests, beside its robust share (37% for the turned run's own years 10–19) and the five-class figures, and the lake totals are summed at full precision.
- **The records**: the robust share of the humid and wet land (23.0%), the README's comparison of robust shares, the status file's next climate task and date, the soils row's map, the missing log line, the moisture map's docstring, and the note that main's checkout needs both windows' grids before its tests pass.
- **Climate's code pins** left stale by its commit `12673b97` are recorded above, for the products' next builds.
- **The resolver's masks**, which the review found looser than intended though they hid nothing this round, are tightened in the drive's scripts for reuse.
- **The quotation in ecology finding 15**, which predates this round, stays as it stands.

## Still to do from the integration

- **The independent review** and its answers, then main's fast-forward at the author's word and the branches brought up to main.
- **A further axis beside the moisture**, among the conditions the builders will not set, with the Codex mockups of 7 November tried on the exact data where it exists ([register](../decisions.md), 7 November).
- **The soils study and the lakes' distribution** among fresh, soda and seawater lakes, with the water-balance screen of the dry belts' closed basins, after the integration.
- **Geology's fallback** to hydrosphere's worktree, dropped at geology's next rebuild.
- **Climate's stale code pins** in the waves', the snow's and the cloud twilight's products, rebound or rebuilt at their next build.
- **The CM1 boxes on the turned years**, which would give the moisture map a gap scale and night ratio of the turned climate's own.
