# The joint integration of 8 October 2026

Five states of the work became one on 8 October 2026: main at `2f2e0a1`, sea appearance at `46681f4`, the solar
shield and habitat array at `eeaecbc`, atmospheric electricity at `24d236b` and infrastructure at `625f60e`. The author
asked for "a proper integration and joint merger" of the branches, with a synthesis of what they say together. The
merges ran on `research/branch-integration-analysis` in the worktree `../terluna-integration-analysis`, in the order
Codex's review proposed and the author accepted. Every source branch is an ancestor of the result, and every merge is
an ordinary merge commit by exact head. The synthesis is in
[research/studies/joint_synthesis](../studies/joint_synthesis/README.md).

Codex's review, written before the merges, stays here as the plan the integration followed: its
[merge plan](merge-plan.md), its reports on [main and sea appearance](main-and-sea.md),
[infrastructure](infrastructure.md), [atmospheric electricity](electricity.md) and the [solar shield](solar-shield.md),
and its [evidence](evidence/) with the [audit script](audit.py), which checks the pre-merge heads and pairwise trees.

## The commits

| Commit | What it does |
|---|---|
| `05810ce` | Codex's review, kept as written |
| `6498da1` | Merges sea appearance (main is its ancestor; no conflicts) |
| `979c5b1` | Merges the solar shield (combines cleanly with sea; the constants hold both branches' keys) |
| `eff585a` | Keeps the magnetic comparison's 26 September design point in solar_shield_array/historical with its reading rule; compares the natural-return products' constants by value; removes a duplicate SUN_GM |
| `f31d8b3` | Merges atmospheric electricity (five shared documents resolved by meaning) |
| `cf371ae` | Merges infrastructure (six shared documents resolved by meaning) |
| `5b597ee` | Reproduces the corrected storm's diagnostic with the current analysis script (every value the same) |
| `a14344d` | Reruns the infrastructure studies on main's corrected climate, with one port geometry and shared constants by value |
| `d1c9f34` | Sizes the sky ships' wings for the corrected storms' gust (21.7 m/s; the wings unchanged) |
| `9a51250` | Brings the infrastructure prose to the corrected products |
| `e219b82`, `741d713`, `0a489a9` | The synthesis: Moon-wide lightning and the tower screen, the fleet's night light and haze, the synthesis study |
| this commit | The shared records brought together: status, research index, decision register, plan and this record |

## How the shared documents were joined

Git combined most of the tree; where both sides changed a shared document, the resolution kept every branch's
entries and took each topic's newest statement.
- **Atmosphere README:** the shield's middle-atmosphere row and the electricity row.
- **CM1 README:** main's wave exports and the electricity study's pointer.
- **GCM README:** main's Sun and radiation checks beside infrastructure's wind exporters.
- **Habitation README:** one condition naming the optical-comfort screens and the flyer studies.
- **Visualization README:** main's newer reference-renderer row with the sea scenes, and infrastructure's port and
  fleet rows.
- **Research README:** every study and command, main's climate, geography and habitation rows, infrastructure's
  engineering row.
- **Decision register:** electricity's programme of 2 October before main's wave decisions of 3 October; main's
  living-seas decision before infrastructure's three sections.
- **Status record:** main's seas topic, the shield's protection task, the electricity topic, and infrastructure's
  lists and condition tokens joined without duplicates.

The shared constants hold sea's Earth albedo (0.242) and phase curve with the shield's DE440 solar GM, proton mass,
planetary GM and time units. Products that recorded the whole constants table by hash keep that hash as their run
record; products that read named constants are checked by value, as shared/provenance.py prescribes.

## Repairs

- **The magnetic comparison** read the design point of 26 September, which the shield rederived six times on 7–8
  October. That design point is kept byte for byte with a manifest naming its producer, successors and reader; the
  comparison's test checks it against the kept copy, and the rerun on the current design point waits for the closure
  trajectories this machine lacks.
- **The storm diagnostic** recorded its analysis script as of `94082be`; rerun on the archived run with the current
  script (82 seconds) it reproduces every value.
- **The infrastructure studies** read the corrected design run (`A28_dim5_moon`) and equatorial ring by name, main's
  hash-checked geography grids, and Earth's gravity and the half-month night from the shared constants. The chosen
  form records its width profile, and the sky fleet and its drawings take the crown from it. Each product records the
  constants it reads and the modules it imports, and its test checks them. The reruns reproduce the infrastructure
  review's corrected-input values: design gust 32.4 m/s, chosen form 98 Mt, screened rotors 140% of the port's use,
  the 750 m liner's nose load 0.43 MN, the high platforms 115 and 460 kW.
- **Prose** in the four infrastructure studies, the metropolis brief and the register's status column follows the
  corrected products. Corrected while reading: the water column's pressure at lunar gravity, the first disk's
  structure, the brace rings' spacing, the lakes round the summit and the crop stand's carbon.
- **Stale statements** in the shared records now say what stands: the thermal column's infrared cooling, the
  electricity topic's canonical run and adopted rules, the photometric dusk, the light calendar and Earthlight, the
  old port drawing's place, and the plan's work order.

## Inputs outside Git

Main's checkout holds every input the branches keep outside Git, copied from the integration worktree once main had
moved, each matching its manifest: the sea branch's Earthlight and water-column inputs with the two Earth textures the
light calendar reads, the shield's FISM2, X-ray, cross-section and optical-constant inputs and its limb-heating cache,
and the electricity branch's three input tables. The CM1 and GCM run archives stay on the research drive, linked from
climate/crm/runs and climate/gcm/runs; the synthesis's run logs are there under research/runs/joint_synthesis.

The five worktrees were then removed and their branches deleted locally. The outputs only they held moved to the same
paths in main's checkout: the sea-appearance renders (visualization/sea-appearance/results), the sky-fleet and
summit-port drawing builds with the ring0 drawings, the integration's run logs (research/runs/integration_checks) and
the shield's check logs. Products that recorded absolute input paths inside a worktree keep them as run records.

## Checks

`make check` ran before every commit except `5b597ee`, whose one changed field was checked by the suites that read it, with the next full check passing. On the final tree: the layer check over 859 files with no violations, 1,345
Python tests passed and 14 skipped (seven want the SWAN or SWASH executables, seven numba, neither installed here), 191 JavaScript tests passed with the
light calendar's, the provenance check PASS with the protection inputs restored, and the ensemble validator PASS (14
snapshot files, five paper workspaces, integrity only). The light calendar's build in the check fetches its two public
NASA Earth textures. The ensemble's planning snapshot, papers and seed records, and research/provenance.json, are the
same as at every source head.

## Still to do from the integration

- Electricity's code keeps local copies of Boltzmann's constant, the elementary charge, the speed of light and a
  rounded gas constant; moving them to the shared constants changes its producers' hashes, so it waits for the
  products' next regeneration.
- The fire study screens the earlier narrow bands; the chosen form has no fire study yet.
- The sky-fleet drawings are built outputs; their re-rendering with the nine-berth crown waits for the author's look.
- The port's programme follows the square lattice until the chosen form has its own.
