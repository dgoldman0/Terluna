# The joint integration of 12 October 2026: the three domains

Four states of the work became one on 12 October 2026: main at `2e98de3`, resources at `98d5ae9`, ecology at
`d81f833` and provisioning at `28a7bec`. The author asked for it in one message: "Integrate now. Yuo may create an
integration branch. Afterwards let's bring the three branches up to the new main. Run the food blend before merge."
The blend ran on provisioning first (`28a7bec`). The merges ran on `research/joint-integration-2` in the worktree
`../terluna-joint-integration-2`, resources first, then ecology and provisioning, each an ordinary merge commit by
exact head, so every source head is an ancestor of the result. The synthesis is in
[research/studies/domain_synthesis](../studies/domain_synthesis/README.md).

## The commits

| Commit | What it does |
|---|---|
| `83cd79d` | Merges resources: its evaluation, industry, nitrogen extraction and phosphorus from the waters |
| `428231c` | Merges ecology: the living sky, aerophytes and sky reefs, taxonomy, living light and machines in the wild |
| `9b0f1a8` | Merges provisioning: people on the Moon, food and the blend, cities, distributions, parks and the network |
| `7c230dd` | The four refreshes of main's studies, every former cross-branch input read in-tree, the cascade rebuilt |
| `b0537b8` | The prose brought to the refreshed products, the refreshes followed into the readers they missed, the register's status cells |
| `6964105` | The three domains together (heat, land, gases, phosphorus and the sky's layers at the author's choices), with this record, the research index, the plan and the status |
| `8b401d3` | The synthesis corrected after an independent review: town size, reef passes, the glow's dimming, heights, the sky's phosphorus need and the stated assumptions |
| this commit | The author's decisions on the synthesis's list, with the recommendation for distributing the array's power |

## How the shared documents were joined

Git combined most of the tree. Where more than one side changed a shared document, the resolution kept every
side's entries.
- **Decision register:** main and each branch had recorded the 9 October cross-domain follow-up; one copy remains,
  in ecology's wording, with a note naming main's summaries and the roaming row's copy under Life and people. Every
  section is in date order, from the follow-up and the aerophytes' light history to the machines in the living world.
- **Lane configuration:** both new domain lanes, provisioning and resources, in AGENTS.md, pyproject.toml,
  shared/README.md, the layer check and its test.
- **Research README:** main's follow-up paragraphs, with the current name for the aerophytes, beside ecology's.
- **Plan:** main's text, whose summary already held ecology's sections.
- **Status record:** merged topic by topic at the JSON level: lists by a three-way diff, condition tokens as tag lists.

## Reading the former branches' inputs in-tree

The branches' studies read each other's products at pinned commits through Git. With all of them in one tree, the 28
pins in seven readers (food systems, distributions, ways of living, land cities, daily life, sky ecology and living
light) read their files in-tree, each checked against its sha256; their products record `read_from: in-tree` for each
pin, except land cities, whose one pin records no commit. The tests that compared a pin's Git blob compute it in Python.
Main corrected the joint ledger on 9 October after resources and provisioning had read it, so their products were
rebuilt on the corrected ledger at the merges; only their recorded input hashes changed.

## The four refreshes

The author ordered all four on 12 October (ecology, question 14). They are main's studies, refreshed on the merged
tree, with everything downstream rebuilt until every product's recorded hashes matched its files.
- **Plant twilight from the solved sky.** The canopy model takes its low-Sun and twilight light from the
  solved-column sky under the chosen shield, in place of the exponential-column atlas scaled by the shield. Its
  twilight is dimmer: dark hours rise from 238 to 254 at the equator and from 34 to 120 at 60°, gross photosynthesis
  falls 1.3% and 4.6%, and the never-dark line on the clear ground moves from 60.4° to 64.5°.
- **The forests' leaf area and air.** The megaforest and forest-patch models take a leaf-area index of 6, the canopy
  model's optimum, in place of 22, and the design case's 294.9 K air in place of 288 K. The reference root threshold
  rises from 23.5 to 23.9 m/s, right-censored cases fall from 37 to 20, and the patch takes more waves to fail.
- **Developed fines in the living waters.** The fines are developed soils' particles with Earth's coastal
  non-algal particle optics (Babin et al. 2003a, 2003b; Petzold 1972) in place of raw Apollo soils. The productive
  coast leaves 547 cd/m² of its own light (156 before) with its 1% depth at 12.4 m (9.5), and the river mouth 1,334
  cd/m² (186) at 2.3 m (1.6). In Nubium at midnight the Earth's glitter keeps its 0.069 cd/m² but the brighter sea
  around it takes it below the scene's rule for recording a band.
- **The lunar-cycle register** reads the corrected plant, fruit and drainage products and takes the design case's
  294.9 K air (294 K before): rivers carry 229,418 m³/s (279,167 before), and the coldest hour in the plant model's
  latitude bands is 19.0 °C, still above the GCM's coldest 3-day mean of 14.3 °C.

Downstream, the people screen, the aerophytes' photoperiod and nursery ages, the food study and the blend, the
distributions, ways of living, living light, the resources screens and industry moved with them: the blend's fields take
0.68–1.86 million km² at 8 billion (0.67–1.81 before), and resources' industry, with its runoff scaling corrected
(below), puts the weathering cycle at 706–7,060 years and the kilns at 35–350 GW (652–6,520 years and 38–380 GW before).
C′'s fits are unchanged. Value-by-value lists of the refreshes' changes, `9b0f1a8` to `7c230dd`, are on the research
drive (`terluna-research/integration2/diffs`); industry's corrected figures came after, in `b0537b8`.

## Follow-through found while bringing the prose up

Bringing the prose to the products found five more places that the refreshes reached:
- **Resources' sea phosphorus** read its six inputs from the ecology and provisioning branches through Git at their
  pinned commits; it now reads them in-tree like the other 28 pins.
- **Resources' industry** scaled the CO2 study's weathering by the first screen's runoff rescale, which compares the
  drainage product's runoff with the register's; the rerun register brought it to 1.0. Industry now scales by the
  present runoff over the CO2 study's own 349 mm, the correction its comment named.
- **The register's air** was still the earlier 294 K, one of the figures ecology's state listed as stale; it is the
  design case's 294.9 K.
- **The taxa register's module texts** quoted the earlier plant figures; they follow the refreshed products, as the
  taxonomy does.
- **Names of the old method** in code: the canopy product's scaling key, the aerophytes' photoperiod notes and the
  food study's evidence now name the solved sky and the domains, and the water-colour figure draws the
  developed-fines waters with the raw fines for comparison.

## Outside Git

- **The sea-appearance regimes** write a 12.5 MB `regimes.npz` to the study's run folder. The integration's run went
  to `terluna-research/sea-appearance-runs-integration2`, which links every other entry of the main run folder.
- **Numba.** The regimes need numba, which no environment here had; `terluna-research/venvs/optics` holds numba
  0.68.0 over the system packages, made for this run.
- **The refresh runs** of the megaforest and forest-patch models, from whose summaries their results were rebuilt.
- **Check logs** are in `terluna-research/check-logs/integration2_*.log`, and the prose pass's reports in
  `terluna-research/integration2/prose`.

## Checks

`make check` ran before every commit: the layer check, the Python and JavaScript tests, the provenance check and the
ensemble validator. Python tests passed at each: 1,387 with 15 skipped after resources, 1,778 with 15 skipped after
ecology, 1,974 with 16 skipped after provisioning, 1,975 with 16 skipped after the refreshes, 1,982 with 16 skipped at
the prose and synthesis commits, whose working tree held the synthesis's seven tests each time, 1,984 with 16 skipped
after the review's corrections added two, and 1,985 with the array-power test. The ensemble's planning snapshot, papers
and seed records, and research/provenance.json, are the same as at every source head.

## Still to do from the integration

- **The sea-appearance illustrations** show the seas with the raw-soil fines; the products now give brighter seas,
  and the illustrations' READMEs say so. Their reference renders can be remade from the developed-fines waters; the
  generated images wait for the author's look.
- **The sea scenes' land** is still bare Apollo soil, a placeholder beside the developed-soil waters.
- **The Earth-visibility result** records Nubium at midnight before the refresh; it runs on an external eye model.
- **C′ with the blend:** the distributions study still places C′ with the power-heavy mix; the synthesis gives the
  totals with the blend, and placing it band by band is next.
- **Ecology's land phosphorus loss** takes Earth's river concentrations at the register's former runoff, in the
  people screen (0.063–0.178 g/m² a year) and the sky-ecology screen's dissolved part; at the rerun register's runoff
  it would be about a fifth lower. Resources' notes keep the bound until ecology moves it.
- **The aerosol study's night** is the 238 dark hours it took from the register before the refresh; the register now
  gives 254.
- **Resources' sea phosphorus** keeps its own food-system stand-ins. Binding the food study's per-system phosphorus,
  as it planned for the integration, would close a loop with the food study, which reads sea phosphorus.
- **The aerophytes' summary of 9 October** ([floater_viability.md](floater_viability.md)) keeps its old name as
  written; the study it summarized is now in this tree as
  [research/studies/aerophytes](../studies/aerophytes/README.md).
