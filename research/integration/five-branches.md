# The joint integration of 17 October 2026: five branches

Six states of the work became one on 17 October 2026: main at `fadb0ab`, provisioning at `f134512`, resources at
`c4ef0ed`, array-industry at `684f329`, ecology at `e790666` and Jura Island at `6424aa6`. The author asked for the
branches to be brought together soon, as they had drifted apart while open, on a separate integration branch, with
main moved up to it by fast-forward once the merged tree passed its checks and the work was finished
([register](../decisions.md#the-third-joint-integration-17-october-2026)). Ecology's 20-year GCM test finished during
the merges, and its record was committed on `domain/ecology` (`e790666`) before that branch merged. The merges ran on
`research/joint-integration-3` in the worktree `../terluna-joint-integration-3`, provisioning first, then resources,
array-industry, ecology and Jura Island, each an ordinary merge commit by exact head, so every source head is an
ancestor of the result. The synthesis is [research/studies/climate_targets](../studies/climate_targets/README.md).

## The commits

| Commit | What it does |
|---|---|
| `3e5a6e6` | Merges provisioning: beamed power chosen, C′ placed with the food blend, the power screen |
| `62f6a81` | Merges resources: the carbon cycle, the seas' chemistry, the terrain and its deltas, the three groundwater stages carried into the carbon cycle, the fresh shores and the soda lakes |
| `d8a82ba` | Merges array-industry: the fleet's glow and its roll, the night roll's attitude, the re-flown layout, the tile hardware, the materials and supports, the geometry screen, the coupled flexible model and the constraints review |
| `e790666` | On `domain/ecology`: where the corrected Moon settles, and the target climate opened in the register |
| `e0a2fa5` | Merges ecology: the gases and their pair, where light limits growth, the GCM with the gases, the glow and the dimming, the curved air's light and the corrected radiation, the 45° rings, the longwave refit, ice optics and the settled test |
| `43c6dec` | Merges Jura Island: the 3-D CM1 box, its shores at LOLA's resolution and its analysis |
| the synthesis commit | What each target climate takes, the glow's efficacy carried into the glow study, and the records brought into the merged tree |

## How the shared documents were joined

Git combined most of the tree. Where more than one branch changed a shared document, the resolution kept every side's
entries.
- **Decision register:** every branch's rows stand, with one Decision cell reworded (below, after the review). The dated
  sections follow the order in which each first appeared, from the gas pair of 13 October to the target climate of 17
  October, with Jura Island's among them. In the domain synthesis's decisions, ecology's methane rows and provisioning's
  power rows stand together.
- **Plan:** each branch's paragraphs, provisioning's numbering of the next work with ecology's finished first item and
  array-industry's glow, then brought up to this integration (below).
- **Status record:** merged topic by topic at the JSON level, condition tokens as tag lists and lists by a three-way
  diff, with the domain synthesis's next task written anew for both sides.
- **The domain synthesis's README:** ecology's finished gas item and its paragraph on the gas pair beside provisioning's
  numbering and its power decision.
- **The CRM README and the CM1 runner's tests:** ecology's later additions, the 45° rings and the sounding above the
  GCM's top, with the Jura box's section and tests after them.
- **Products and their pinned scripts:** where two branches had rebuilt the same product, every conflict lay in its
  recorded hashes, so the owning branch's version stood and the settle cascade rebuilt what was stale. Every rebuilt
  product matches the branch that last changed it value for value, compared field by field with hashes set aside.

## Products rebuilt outside the cascade

- **Provisioning's power screen** binds the synthesis's ledgers and the array's heat. The settle tool did not list it;
  it does now.
- **The symbioses product** records its inputs relative to its study's folder, which the settle tool did not resolve.
  It does now, which retires the manual rebuild the sky ecology needed whenever its parent moved.
- **The reduced-gases product**, ecology's line-by-line forcing of methane and N2O, pinned provisioning's population
  product, whose recorded hashes moved in the array-industry merge while its values stayed the same. The line-by-line
  columns depend on neither branch's changes, so the product was re-derived through its own `results()` from its
  stored columns, with its run's seconds kept, in place of about 38 minutes on six workers for the same numbers
  (`session-tools/rederive_reduced_gases.py` on the research drive). Every value is identical.

## The glow's efficacy

- **Carried at the author's direction of 13 October.** The glow study's infrared efficacy is the GCM's 0.71, from the
  held glow branch of the gases, glow and dimming study. It is a constant in
  [glow.py](../studies/solar_shield_array/glow.py) with a test against that study's product, which takes the glow's
  answer from the glow study and so cannot be bound back. The recommended answer's dimming is 2.0% of the design
  sunlight (2.7% at an efficacy of 1), about 1% below the domain synthesis's convention, and gives back 45–67 TW. Every
  answer in [glow.md](../studies/solar_shield_array/glow.md) follows; the absorbing dimmer returns 25% of its dimming as
  glow as built and 14% with the night roll.
- **Downstream:** the tile design, the support comparison, the flexible benchmark and the geometry screen were rebuilt
  on the new glow product, their values unchanged but for timings.
- **The glow note** follows the author's setting aside of the annulus absorber and the ultraviolet mirror, as the
  constraints review's finding 14 asked, and the gases, glow and dimming study reads the glow's answer from the
  in-tree glow product, with a test holding the two together.

## Into the merged tree

- **Pointers to other branches' files** became links: the register's rows on question 10, the soda lakes and the sea
  chemistry; Jura Island's study, ecology's README and questions, geography's Jura bullet, the shield study's README
  and provisioning's opening list.
- **Statuses:** the curved air's light and the corrected radiation's rings done; Jura Island done, with its rerun
  waiting; the glow's dimming at the GCM's efficacy, and the 0.71 carried.
- **Records the glow screen corrects,** annotated where they stand: provisioning's finding 1 (5.7–9.9 W/m² of glow,
  against the screen's 9.4 W/m² as built with 4.7 W/m² of returns, and 6.0 W/m² rolled), the joint synthesis's night
  light (the films' specular light, 1.6 lux as built), and the domain synthesis's heat and next work. Two need model
  work and are on the plan's next work: illumination's fleet light on the layout's radii, and the tile inventory
  against the eccentric rings' length.
- **The shared records:** the research index and the integration index, the plan (the branches' work since 13
  October, the waiting list and the next work ranked anew), the status record and a note at the head of the
  constraints review, which otherwise stands as written on 14 October.

## The synthesis

[What each target climate takes](../studies/climate_targets/README.md) answers the author's open question across the
domains, from the GCM's two runs and the domains' products. Keeping the warm Moon, near 307 K, still takes the 3.9–5.3%
of dimming that answers the gases, the glow and the beamed power; left unanswered they would warm it about 4–7 K more,
and at full watts they would ask its clear sky to send out 7–11 W/m² more than a saturated column can. The GCM holds
the warm Moon stable by its restoring gain of 1.65–2.16 W/m² a kelvin, on air that stays unsaturated, and its cloud
warms it 6.0 W/m² net. 300 K takes 5.2–6.5% more of the sunlight and the design temperature 8.9–11.3% more. At 307 K
the wet-bulb temperature passes 28 °C for a third of the time near the equator and in the lowlands, and the
comfortable air lies 8–12 km up; at 300 K the ground is warm and tropical, with seas near 28 °C. The water the author
keeps delivering to balance the crust brings 1.5–9.5 W/m² of heat in its fall through the first thousand years unless
it is lowered with energy recovery, as the build's deliveries already must be.

## The review

An independent review of the synthesis commit (`8e60aea`) reproduced the product and confirmed the register's rows,
the status merge, the links and the heat items' arithmetic, and found fourteen faults, corrected in the commit after
it:
- **The feedback** paired one fit's slope with the other fit's balance; each fit now keeps its own, which narrows the
  dimming at 300 K from 5.0–6.8% to 5.2–6.5% and at 295 K from 8.7–11.6% to 8.9–11.3%.
- **The stability reading** compared all-sky absorbed sunlight with the saturated column's clear-sky ceiling and gave
  the cloud a share in holding the balance, though it warms the Moon on net. The comparison is now clear sky against
  clear sky, with the added heat at full watts, and the restoring gain and the Manabe–Wetherald columns lead.
- **The costs** now name each target's own dimming and give the shared dimming's cost beside them; the warm Moon's land
  rain falls to 3.10–3.27 mm a day with it.
- **Two sunlight scales**, the target's own design case and the plan's whole dimming, in place of one.
- **Figures:** the land's rain rises 48% (the 44% was the whole Moon's, and 43% when unrounded); the absorbing dimmer
  returns 10% of its dimming as glow with both halves rolled; the glow note's as-built and heat-mirror dimmings
  follow the efficacy; the convention matches to about 1% (1.2%); the beamed power is the beams arriving at the
  receivers and crowns; the GCM README's balance gives each fit its own window.
- **Pointers:** the domain synthesis's plan item and its traffic item, two status tasks still waiting on the 0.9685
  run, and provisioning's finding 1, whose note now follows the decision.

The register's row of the target climate, written on ecology with some of the author's own phrases in its Decision
column, was put in my words at the synthesis commit; the phrases stay in its Source column.

## Outside Git

- **The worktree** linked main's ignored inputs as the domain worktrees do, and held hash-identical copies of the run
  folders only the ecology and array-industry worktrees held: `research/runs/gases_glow_dimming`, `sky_ecology`,
  `solar_shield_array_flexible_benchmark` and `solar_shield_array_geometry_screen`, the last two rewritten by their
  rebuilds there. `research/runs/climate_targets` links `terluna-research/climate-targets-runs`, which holds the
  synthesis's extract. Main's checkout took the same folders and link with the fast-forward.
- **The session tools** on the research drive: the settle tool gained the power screen, study-relative inputs and the
  synthesis (its earlier copy is `settle.py.2026-10-17.bak`); the settled test's analysis scripts are in
  `climate_settled_2026-10-17`; the check logs are `make_check_integration3_*.log`; the merge helpers, among them the
  field-by-field comparison of products with hashes set aside and the status record's merge at the JSON level, are in
  `integration3`.

## Checks

`make check` ran before every commit: the layer check, the Python and JavaScript tests, the provenance check and the
ensemble validator. Python tests passed at each merge: 2,000 with 16 skipped after provisioning, 2,066 with 19 skipped
after resources, 2,334 after array-industry, 2,412 after ecology and 2,421 after Jura Island, each with 19 skipped, and
the synthesis commit 2,430 with 19 skipped. Ecology's record of the settled test passed 2,063 with 16 skipped on its own
branch. Two merges first failed on a product outside the settle tool's reach (the power screen, then the symbioses), and
one on the reduced gases; each was rebuilt or re-derived as above before its commit.

## Afterwards

Main's own checkout passed `make check` after the fast-forward, with 2,430 Python tests passed and 19 skipped. At the
author's word the same day, the integration's worktree was removed and its branch deleted locally, once every folder
it held outside Git was found to match main's checkout or to link the research drive; GitHub keeps the branch for the
author. The five branches then moved up to main by fast-forward. Each was an ancestor of main, its worktree held
nothing uncommitted, and none of main's new files would land on a file already there. Each worktree took the run
folders it lacked from main's checkout, the links to the research drive's run folders and copies of the four folders
above, so its checks read what main's do, and Jura Island's links into ecology's worktree became copies. Each worktree
then passed `make check` with the same 2,430 Python tests passed and 19 skipped; the logs are
`make_check_moved_up_*.log` with the session tools.

## Still to do from the integration

- **The target climate,** the author's choice; then the design case for it, the glow pair from it, and geography's
  lakes, the atlas globe and the studies that read the design climate, which stay on the design case until then.
- **The heat beyond the sunlight:** the power screen's heat in the ledgers, and the water top-up's energy recovery.
- **The fleet's records:** illumination's fleet light on the layout's radii and the tile inventory against the
  eccentric rings' length.
