# Joint branch analysis and integration proposal

**Draft for author review, 8 October 2026.** The five states can be brought together without choosing a new climate answer, reopening capped research, or revising the manuscript ensemble. The proposed integration preserves their histories, reconciles current decisions and evidence, repairs product provenance, and updates the infrastructure calculations to the accepted climate inputs. It does not yet establish a working shield, a safe tower, a closed electrical circuit or a validated ecosystem.

This review is on `research/branch-integration-analysis`, in `/home/kir/Documents/Projects/terluna-integration-analysis`, created from current main. Only this analysis package has been written. There are **no commits, branch merges, rebases, source-worktree edits or new scientific simulations** from this task. The source heads are pinned below. Temporary Git trees were computed outside the repository to inspect pairwise conflicts; no merged checkout, index or ref was created. This package remains uncommitted until the author decides whether to adopt it and authorize execution.

## Reading the evaluation

- [Main and sea appearance](main-and-sea.md): shared baseline, light and water findings, dated sky calendar, evidence and downstream effects.
- [Infrastructure](infrastructure.md): author choices, corrected-input consequences, port geometry, aircraft, fire and integration inputs.
- [Atmospheric electricity](electricity.md): corrected storm results, aerosol, adopted rules, model limits and archive/build hazards.
- [Solar shield](solar-shield.md): architecture history, final retention screen, adopted choices, open system gates and industry cap.
- [Merge plan](merge-plan.md): proposed order, file-level resolutions, product treatment, input restoration and acceptance gates.
- [Snapshot](evidence/snapshot.json): full SHAs, histories, changed paths and all ten pairwise conflict analyses.
- [Local state](evidence/local-state.json), [calendar compatibility](evidence/calendar-compatibility.json), and [checks](evidence/checks.json): bounded evidence supporting this review.

The branch reports refer to source paths at their stated commits. Some of those paths do not exist in this main-based analysis worktree yet. Read them from the existing source worktrees or with `git show <pinned-sha>:<path>`. This is intentionally an analysis branch, not a partial implementation.

## Where each state stands

Counts are relative to each branch's merge base with main. “Main only” counts work absent from that branch's ancestry, not a count of fixes that all require manual application.

| State | Pinned head | Main only / branch only commits | Changed paths | Current position |
|---|---|---:|---:|---|
| Main | `2f2e0a1` | 0 / 0 | — | Corrected climate, waves/tides, optical comfort/cloud twilight, shared constants by value; baseline for integration |
| Infrastructure | `625f60e` | 64 / 15 | 71 | Capped; adopted port/metropolis/flight planning; cap recommendations pending; original study products await corrected-input integration |
| Atmospheric electricity | `24d236b` | 27 / 39 | 78 | Capped stages 1–2; corrected 530-flash result; global circuit and upper boundary remain stage 3 |
| Solar shield and habitat array | `eeaecbc` | 0 / 62 | 500 | Capped; lead ring architecture and revised loss screen; complete-system gates open; array-industry branch planned |
| Sea appearance | `46681f4` | 0 / 25 | 255 | Uncapped by the author's choice; usable optics/scene/calendar products plus acknowledged unfinished work |

All five source worktrees were clean at the start of the review. The shield's local cap is one commit ahead of its remote-tracking ref; the review includes that **local** commit. No fetch or push was performed, and remote server freshness is not asserted. `restructure/lanes` (`d73ce9f`) is another local branch name, but is already an ancestor of main and contributes no unmerged work. It is excluded from the four active branches, without deleting it.

The older infrastructure and electricity branches need main's later corrections. Shield and sea already contain main. These relationships favour ordinary history-preserving merges over cherry-picking selected “good” files, which would obscure the rejected experiments, design decisions and evidence progression.

## What the branches give together

Main establishes the common environment. Electricity adds storms, charge, acoustic consequences, aerosol and an uncertain nitrogen source. Infrastructure asks what people can build and how they move through that air. Shield asks how the atmosphere and industrial support can persist. Sea appearance makes light, water and location usable as quantitative environmental conditions as well as views.

Together they expose a more concrete world than any branch alone, but they do not form one fully coupled model. The links fall into three groups: calculations already consuming compatible products, interpretation supported by comparing the branches, and future model interfaces. The distinction matters: combining repositories does not cause aerosol to enter the sky solver, storm charge to interact with the summit tower, or the ring fleet to complete its operational control gates.

| Intersection | What the combined evidence resolves | What it now asks |
|---|---|---|
| Main climate × infrastructure | Accepted corrected inputs preserve much of the air-by-height/flight argument but raise wind-sensitive loads; cap already reports temporary rerun targets | Commit truthful corrected products; then study height-dependent gusts, summit terrain and wakes separately |
| Electricity × infrastructure | The 35–45 km flight band and 35.6 km crown overlap storm charge and median flash initiation near 35.8 km | Grounded-conductor/leader interaction, hydrogen berths, aircraft triggers, warning, closure and rescue |
| Sea light × infrastructure | Long usable twilight replaces five-hour/whole-night-dark descriptions; visual light and energy supply must be separated | Summit-altitude/horizon light, night operations, energy storage and human schedules |
| Shield aperture × atmospheric retention | Slant heating, gaps, solar-history bands and infrared cooling revise earlier loss estimates and the maximum-year/latch story | Full composition/transport feedback, failed-cell coverage and real plasma/magnetic protection |
| Shield geometry × industry | Ring layout, film mass and heat placement make replacement industry and actuator power concrete | Qualified trim/store devices, mass feedback, recurring operations and full source-to-use logistics |
| Electricity × shield | Stage 3's upper boundary and magnetic scenarios can now be specified against more developed UV/loss cases | Conducting ionosphere, global circuit, particle transport and night persistence; current storm statistics do not answer these |
| Electricity × sea optics | Native ice versus radiative cap is explicitly distinguished; fog/aerosol is an identified uncertainty shared by both | Aerosol optical transport and scheme-sensitive fog; a view cannot repair the cloud field's radiative history |
| Waves × sea appearance | Existing wave/tide products now inform slopes, coast views and source reflection | Whitecaps, living-sea films, sediment, remaining basins and missing swell histories |
| Water optics × biology × lightning | Underwater light and nitrogen-source bounds become usable conditional inputs | Full nutrient/oxygen/carbon cycles; no global fertility budget follows from one lightning box or guessed plankton |
| Surface city × orbital array | Demand, heat, freight and exhaust can be assigned to explicit endpoints | Surface–orbit energy/material/outflow ledger; no adopted replacement of fusion or extension of S7 to all launches |

The strongest new design tension is **height**. Low gravity makes tall construction and efficient flight attractive, but the same tall atmosphere carries deep storms and charge into the selected flight band. The strongest scheduling correction is **light versus power**: twilight supports outdoor activity for days without supplying the energy assumed in a solar-day calculation. The strongest industrial correction is **heat location**: energy absorbed by thin orbital film is not automatically useful work, and energy delivered to the city ultimately heats the city.

These are interpretations of the repository results, not new numerical or empirical findings.

## What must remain locked and what must remain open

| Category | Preserve in the integrated record | Avoid promoting into a new decision |
|---|---|---|
| Project-wide author choices | Open Moon; 1.2 atm; 28% water scenario; 5% dimmer; corrected drier land; named flight band; safety rules; climate pause and laptop limits | A single resolved climate/comfort answer, a new loss-budget choice, rented compute |
| Infrastructure | Summit location, approximate population/height, six-leg round form, lower rings/upper disks, metropolis scale, conditional flyer planning basis | A selected crown liner, a new city power source, complete structural/fire safety, cap recommendations |
| Electricity | Three stages; adopted breakdown, leakage and ground-leader rules; retained radiation cap; aerosol held for stage 3 | Stage 3 completion, optional corona as the default, global flash rate, confident polarity or nitrogen yield |
| Shield | Ring lead, held-screen fallback, integrated comparison gates, S7, traced heat, primary cycle-21 maximum, radius-layout direction, cap | Gate pass, all-ring operational safety, qualified hardware, primary regional magnet selection, initiated array-industry campaign |
| Sea appearance | Author-approved photometric twilight, Earth brightness/phase curve, living-water guesses, chosen study/scenes and native ice | A new cap, guessed ecology as observed biology, photographic illustration as naked-eye validation |
| Manuscripts and baselines | All five roles/UUIDs, planning snapshot, selected seeds and imported byte pins | Manuscript source admission, substantive editorial approval or exact-build PDF approval |

Tree identities confirm `ensemble/planning`, `ensemble/papers`, `ensemble/seed_records` and `research/provenance.json` are identical across all five heads. No manuscript rewriting is needed for this integration. Later writing should use the established roles: root `cef466e6-d9f8`; physical `7c83e7a3-5955`; biosphere `74a946e2-8f98`; human `babdd5f6-4920`; engineering `155b039f-a0e8`. Biology and ecology stay together; human society remains the human paper's scope.

The master decision register should retain chronology and explicit supersession. “Author adopted the findings as a planning basis” remains different from “empirically validated”, and a cap is different from approval of every recommendation inside it.

## What prevents a simple merge and declaration of readiness

Textual conflicts are concentrated in eight documentation/status files across all pairwise comparisons: atmosphere README, the two climate READMEs, habitation README, research README/decisions/status, and visualization README. Actual execution may reveal different conflicts as earlier ones are resolved. No branch-native scientific code file is a textual conflict in these pairwise tests. That does **not** mean the results are ready to use together.

Four material issues require explicit treatment:

1. **Stale or contradictory current summaries.** Sea's status is still main's old status. Electricity has completed work still described as running or held. Shield retains earlier tile counts and pre-infrared maximum-year conclusions. Infrastructure simultaneously describes a sized concept and components “not yet sized”. Main's research plan also proposes work already done.
2. **Historical outputs beside newer producers.** The shield magnetic comparison fails its current-parent hash test; sky-ships pins an older winged-flight source; electricity's corrected diagnostic pins an older, identified analysis script. These have different repair paths and cannot all be fixed by replacing hashes.
3. **Competing infrastructure geometry.** Square-lattice sizing, chosen round form, old fire compartments and mixed crown drawings analyse different buildings. Integration must label historical studies and stop current consumers from quietly mixing geometries. A full new fire/structural design remains future work.
4. **Isolation beyond Git.** Data archives, model executables, caches, Python environments and node dependencies live outside the tracked tree. A new worktree does not supply them. Writing through symlinked run/build paths could alter other branches' working evidence.

A targeted test already confirms the shield provenance failure. The preferred plan preserves that magnetic comparison with its exact historical parent and a reading rule, keeping a modern rerun pending. It preserves traceable electricity diagnostics similarly unless isolated postprocessing is justified. Infrastructure's already-planned inexpensive reruns provide the natural way to repair its stale product and update inputs. These are proposals for the implementation stage, not repairs applied in this review.

## Proposed integration outcome

Use the existing new branch and workspace, beginning from the approved analysis checkpoint. Merge **sea appearance → solar shield → atmospheric electricity → infrastructure**, with exact pinned heads, semantic reconciliation and required checks at each commit. Sea and shield combine cleanly in the pairwise Git test. They establish the current light, constants and protection record before the older electricity and infrastructure lines are reconciled. Infrastructure is last because its cap explicitly consumes the other four states and its recalculations need their final inputs.

The outcome should be one branch with all source histories reachable, a clean working tree after approved commits, current per-topic indexes, preserved author decisions, explicit historical products, corrected infrastructure outputs, isolated runtime writes and meaningful validation records. Scientific open questions remain open. The detailed [merge plan](merge-plan.md) defines exactly what this requires and where execution must stop rather than falsely declare success.

## Questions the unified branch should make easier to pursue

Priority here is a recommendation for review, not authorization to run:

1. **Crown and storms:** begin with bounded electrostatics on existing charge fields and chosen tower geometry; establish what the discharge model can and cannot say before terrain-specific storm simulations.
2. **One port geometry:** supply the same adopted form to loading, floor programme, berths and future fire/refuge calculations; keep alternative structures as comparisons.
3. **Array industry:** pursue the capped plan's hardware/heat/power questions with explicit surface and orbital ledgers; close actuator mass feedback before ranking a full fleet.
4. **Atmospheric circuit:** define scenario-tagged upper-boundary inputs from shield work, then decide the scope of stage 3; preserve muon-dominated lower-air reasoning and separate dose from conductivity.
5. **Aerosol and optical life:** couple aerosol extinction and weather variability to sky/water light and conditional biological requirements; resolve model-scheme and site differences before interpreting new pictures.
6. **Water and material cycles:** connect waves, spray, light, sediment, nitrogen and biological return to a budget that distinguishes missing transport from missing empirical traits.

The climate discrepancy cuts across all six. It remains a shared uncertainty, not an invitation to restart the paused programme during integration.

## Verification actually completed

The review recorded source ancestry and cleanliness, all ten pairwise Git conflict analyses, unchanged locked trees, selected product/source hashes and ignored-data locations. In the fresh analysis checkout, layer checking passed across 396 files and ensemble integrity passed for 14 snapshot files/five paper workspaces. Focused source-branch checks returned 89 infrastructure tests passed; 34 shield tests passed and one provenance failure; 32 sea tests passed and one optional skip; and two sea-calendar JavaScript test files passed. Shield's own layer check passed across 676 files.

Electricity's previously saved final `make check` log was inspected: 486 Python passed/29 skipped, 176 JavaScript passed, integrity checks passed, protection verification blocked by input state. That is historical evidence, not a test run in this review. No fresh full `make check`, expensive model campaign, literature admission, browser/image inspection, manuscript edit or PDF review was performed. The package's link/JSON/snapshot checks are recorded separately in [checks](evidence/checks.json).

Approval of this draft would adopt a review and plan. Starting merges, making commits, initiating the new research programmes and publishing anything remain separate actions; none occurred here.
