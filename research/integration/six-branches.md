# The joint integration of 19 October 2026: six branches

Seven states of the work became one on 19 October 2026: main at `f03b839`, climate at `82b6ed4`, provisioning at
`a66c0b1`, Jura Island at `c28959d`, resources at `8d31e86`, atmosphere at `a90c23b` and ecology at `b5eff1b`. The
author asked for an integration round like the last, with more branches, main moved up to it by fast-forward once it was
checked complete and clean, and every other branch then brought to match main as the common baseline for new work
([register](../decisions.md#the-fourth-joint-integration-19-october-2026)). Array-industry stood level with main and
held no work of its own. The merges ran on `research/joint-integration-4` in the worktree
`../terluna-joint-integration-4`, climate first, then provisioning, Jura Island, resources, atmosphere and ecology, each
an ordinary merge commit by exact head, so every source head is an ancestor of the result. The synthesis is
[the resort at 301 K across the domains](../studies/resort_301k/README.md).

## The commits

| Commit | What it does |
|---|---|
| `880ab45` | Merges climate: the target climates compared and the resort set at 301 K, the dimming curve and the glow's range, the model's own land mask, the humidity from consecutive records, winter on the tower, the coasts' and islands' water, the deserts as a range |
| `daeef60` | Merges provisioning: the low-gravity design note, the ChatGPT-era habitation decisions, the summit metropolis for the district wind study |
| `1ba66e3` | Merges Jura Island: the flyover in visualization, drawn from the island study's product |
| `4e9749f` | Merges resources: the first seas' finding and the heat levers on the corrected targets |
| `11955b3` | Merges atmosphere: the rain's drop physics, its rates and its sound, rebuilt at the corrected targets |
| `94d17d1` | Merges ecology: the corrected lake model with the wind's deep mixing, the waters' oxygen and anoxia, the lakes by kind, the seas' basins and the waters' methane at 301 K, kept on the targets they ran on |
| the synthesis commit | The synthesis, the comparison brought to the corrected lakes, the targets reading the waters' oxygen in-tree, and the records |

## Before the merges, on the branches

Domain/climate's corrected climate targets (schema 5) superseded the product three branches had built on: ecology's
lake studies and atmosphere's rain (schema 2, ecology's version at `e9c12eb`) and resources' heat levers (schema 1,
main's). Each branch brought its own work up on its own branch, with its own `make check`, before it merged:
- **Resources (`8d31e86`).** The heat levers read the corrected targets and their extract, pinned by commit and hash
  until the branches met, carry the resort at 301 K, and take all nine fixes of the screen's review of 17 October
  (an agent's work, reviewed; two contrasts restated).
- **Atmosphere (`a90c23b`).** The rain's drops, rates and arrivals read the corrected targets the same way and take the
  resort from its 301 K row; the rain column's dry adiabat takes the repository's heat capacity, as the author inclined
  for this rebuild; the rates carry each target's own dimming cut; the sound lab's renders of 19 October follow on the
  research drive. The agent doing it was stopped near the end by auto mode's safety check, which blocked its last
  edits and then its shell, citing content in its own conversation; the author chose that the coordinator review its
  whole change and finish it, which the coordinator did: the blocked test and three README edits, the register's rain
  rows without their contrasts and with the author's hedges, and the plan's two rain items.
- **Ecology (`b5eff1b`).** The lake studies keep the targets they ran on (below), their caches are keyed on the values
  each run reads, and the register records the lake decisions of 17 and 18 October, which had waited for this
  integration so that one rebind would serve.

## The lake studies keep the targets they ran on

The merged targets move the lakes' forcing very little: the land's lapse rates fall 1.2–1.8% under the model's own land
mask, which would cool the lakes' air by 0.026 K at most (0.015 K typically), and the lock's share of the way goes from
−0.002 to 0, the design case the lakes' lock already is; the 301 K share is unchanged. Rerunning the whole lake set
would take six to nine hours on the machine, and the lakes rerun anyway on the 301 K design case. So every lake study
reads a byte-identical copy of the targets it ran on, `research/studies/lake_oxygen/climate_targets_as_run.json`
(ecology's product at `e9c12eb`), checked by hash; the copy also breaks the binding cycle that reading the waters'
oxygen in-tree would otherwise make with the climate targets, which read it.

The lake studies' caches had bound the settle chain's products by file hash, the sea chemistry in every key, so every
edit of the register staled all six lake products and missed every cache. They are now keyed on a fingerprint of the
values each run reads: a cascade that rewrites those products with the same values leaves every key standing, and
every product step rebuilds from cache in about fourteen minutes. Each step has a cache-only mode, which the settle tool
always sets; 1,744 cache entries were re-keyed by hard links, each family with a written proof that only keys, paths
and flags changed, and all six products match their earlier values.

## How the shared documents were joined

Git combined most of the tree. Where more than one branch changed a shared document, the resolution kept every side's
entries.
- **Decision register.** Every branch's rows stand. The dated sections follow the order in which each first appeared:
  after the third integration's, the ChatGPT work confirmed (provisioning), the rain's physics and sound (atmosphere),
  the Jura Island flyover, the lakes of 17–18 October (ecology) and this integration's. The flyover's and the sound
  lab's tests rows, one answer of the author's to one question about both, became one row in the rain's section.
- **Plan.** Each branch's paragraphs, the next work's items from four branches numbered in turn, then ranked anew at the
  synthesis (below).
- **Status record, research index and test paths.** Merged entry by entry; the test paths carry habitation, the
  flyover and the sound lab together.
- **Products and their pinned scripts.** Every conflict lay in recorded hashes, so the owning branch's version stood and
  the settle cascade rebuilt what was stale; every rebuilt product matches its owning branch value for value, compared
  field by field with hashes set aside. The climate targets took domain/climate's version at every merge.

## The settle tool

The settle tool on the research drive (`session-tools/settle.py`, its earlier copy `settle.py.2026-10-19.bak`) now
rebuilds the heat levers after coast water; the lake studies' eight product steps after the resources chain and before
the climate targets, which read the waters' oxygen, always from their caches alone; and the three rain products after
the heat levers, the rates with the climate interpreter that has netCDF4; and the synthesis last. A cascade after an
edit of the register takes about 22–25 minutes and recomputes nothing. The lake products' `seconds` now record each
cache-only product step (lake_oxygen 1.2 s), and their runs' own times stay in the caches and the READMEs. The lake
kinds' `setting.py` keeps its docstring's note of the register "on domain/climate": the file's hash keys the lake
kinds' caches, so even a docstring edit would miss them.

## Into the merged tree

- **Reads in-tree.** The heat levers and the rain read the climate targets in the tree, their values unchanged; the
  climate targets read the corrected waters' oxygen in-tree, which retires the pinned earlier product (`bd0097c`) and
  moves the resort band's oxygen to the corrected lake model.
- **Pointers to other branches' files** became links: the comparison's sources, the targets' README, the heat levers'
  sources and lakes, and the lake studies' notes on the resort's register row.
- **The comparison** of the three targets carries the corrected lakes: the deep lakes' overturns and oxygen by both
  mixing forms, the shallows, the lakes by kind, the seas and the methane, the levers on the corrected targets and the
  rain at 301 K.

## The synthesis

[The resort at 301 K across the domains](../studies/resort_301k/README.md) works out where the branches' pieces meet.
The waters' trophic state sets the air's methane more than the gas pair does: nutrient-poor lakes and seas hold the held
pair with no removal at the methane share's centre and with up to 2.6 Tg a year across its brackets, mesotrophic ones
need 4.7–6.9 Tg. The lakes stay nutrient-poor while people's dissolved phosphorus stays within their room, which turns
on the giants' floors: with 8 billion people at the middle settling, 57–71 g a person a year while those floors stay
oxic and 5–10 g if their depths go anoxic, against a design loss of 44–180 g. On the rivers the seas get almost none and
run ultra-oligotrophic, a question of production for the biodiverse communities the register asks of them, with the
dredged deltas' phosphorus as its lever. The levers' first four, each with the case it belongs to, leave the resort
6.1–8.3% of dimming in all. Watering the deserts, most of them poleward of 45°, takes 35–44 km³ a year and 2.3–3.8 GW
for each 1%, its albedo the larger cost, twice as much of it for each km² at 30–45° as at 60–75°.

## The review

An independent review of the synthesis commit (`a0713f8`) rebuilt the synthesis byte for byte, confirmed the register's
rows and order, the merges' exact heads, the frozen copy, the cache keys, the products' identity across the branches and
the rain work, and found faults, corrected in the commit after it:
- **The levers' pairing.** The resort's dimming with the first four set the dimming's low end against the levers' high
  end, which needs the beams' high end; each case now takes its own glow efficacy and beams' end, in the heat levers
  and the synthesis: the resort is left 6.1–8.3% (5.2–9.3% before), and the four save 43–57%, 17–31% and 11–21% of the
  three targets' dimming (30–81%, 15–37% and 10–24% before).
- **The phosphorus room** took the 352 m lakes' line for every lake; it now takes each kind's line from the lake kinds,
  with the giants' floors oxic and anoxic: 57–71 g a person a year and 5–10 g with 8 billion at the middle settling
  (15–29 g before).
- **The seas.** The author's condition is biodiversity, as the register states it, with the deltas' phosphorus for
  production if needed; the synthesis had read it as a call for nutrient-rich seas. The dust, which can lift the mixed
  seas to 18–28 mg/m³ at the inherited envelope's top, is now given.
- **The methane brackets**, the share's 0.8 and Shih and Bouffard's form, are carried beside the centre; the
  comparison's overturns under the energetic regime take one convention, under which the resort's lakes turn over least;
  the comparison's CO2 growth cost, shallows and food carbon are corrected; the watering pairs each band's lift with its
  own sunlight; the targets' reading rule describes the in-tree oxygen; and the records, the register's status of this
  integration among them, follow.

## Outside Git

- **The worktree** linked main's ignored inputs as the domain worktrees do; `research/runs/climate_targets` links
  `terluna-research/climate-targets-runs-climate`, the extract the merged product binds, and the five lake folders,
  coast water and the drive's other run folders are linked as on the branches; the sound lab's recordings were copied,
  hash-checked.
- **The session tools** on the research drive: the merge helpers (`integration4/hash_conflicts.py`, `reflow.py`,
  `completeness.py`), the record of the day (`integration4/record.md`), the settle logs and the check logs
  (`make_check_integration4_*.log`).

## Checks

`make check` ran before every commit. Python tests passed at each merge: 2,500 with 20 skipped after climate, 2,518
after provisioning, 2,535 after Jura Island, 2,549 after resources, 2,632 after atmosphere and 2,755 after ecology, each
with 20 skipped, and 2,761 at the synthesis commit; on their own branches before merging, resources passed 2,444,
atmosphere 2,517 and ecology 2,557, each with 19 skipped.

## Afterwards

Main's checkout moved up to the reviewed head, `0e6a35e`, by fast-forward on 19 October, as the author asked once the
integration was checked complete and clean, and the seven branches, the five domains, `research/jura-island` and
`research/array-industry`, moved up to it the same way, each clean and an ancestor of main. Every checkout's
`research/runs` now links the research drive's `climate-targets-runs-climate`, the extract the merged targets bind,
with `coast-water` and the five lake folders, and holds the sound lab's three recordings, copied and hash-checked;
atmosphere's local copy of the old extract, identical to the drive's ecology extract, moved to the drive
(`session-tools/integration4/`). Main's checkout and every worktree then passed `make check`, 2,762 Python tests
with 20 skipped (`make_check_after_integration4_*.log`). The review's confirmation found every high and medium fault
fixed and one new error, two tables the reflow had flattened, restored in `0e6a35e`, whose check passed 2,762. The
integration's worktree and its local branch, which equals main, stay until the author says.

## Still to do from the integration

- **The design case at 301 K**, at the author's word, then the glow pair and the lakes' rerun on its weather.
- **The author's choices:** the waters' trophic design, the deserts' watering, food in orbit's carbon, the sea
  chemistry's floor ordering and the summit site product's rebuild.
