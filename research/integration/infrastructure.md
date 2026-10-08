# Infrastructure branch review for unified branch planning

Audit date: 2026-10-08. Review only; no source edits, merges, commits, climate simulations, downloads, or external scholarly verification were performed.

## Snapshot and contribution

`habitation/infrastructure` is at `625f60e83313cadd21a7999a772d93dff07cee17`, based on common ancestor `70eed7cba7c101f6bb802439483c968f46b5fdb4` with current main `2f2e0a104a3a1ea962a6dead4179f522762f997e`. The branch has 15 commits beyond that ancestor and changes 71 files, adding about 50,576 lines and removing 14. A large fraction is structured results and SVG drawings, not new scientific equations. Its worktree was clean before and after this audit.

The branch adds reusable engineering models for towers, wind devices, fires, egress, buoyant ships and winged/rotor flight; two GCM wind exporters; four coupled studies (summit tower, port fire, sky ships, sky fleet); a habitation brief for the summit metropolis; scientific/engineering visualization tools; and an artistic ring-0 world proposal. It does not alter `ensemble/`, `research/provenance.json`, shared constants or the byte-pinned historical imports.

The cap, `research/studies/infrastructure_review/README.md`, was committed on 8 October. It reads exactly main `2f2e0a1`, electricity `24d236b`, sea appearance `46681f4`, and shield/array `eeaecbc`. It is already a substantial cross-branch review, not merely a closing note. Its recommendations explicitly await the author; closing the branch did not adopt them.

Key commits:

| Commit | Contribution |
|---|---|
| `e70b745` | Summit tower and built-in wind-device sizing |
| `23d58cc`, `96321d7` | Author's port scale and programme by height |
| `9f82804` | Fire and evacuation screening |
| `24a6df4`, `faa3941`, `b3c8277` | Author's round six-leg form, population scale, lower rings / upper disks |
| `35ea2be`, `1432efa` | Ship-size and flyer/fleet models and products |
| `aac8c64` | Hand-sized chosen port form |
| `7cbce73` | Ring-0 drawings and world definition |
| `b79f705` | Author adopts flyer findings as planning basis |
| `625f60e` | Cross-branch cap; corrected wind products; deferred integration reruns |

## What is locked, and at what level

The authoritative branch register is `research/decisions.md`, sections “The summit port and metropolis”, “Sky ships and flyers”, and “The infrastructure branch”. Preserve the source dates and epistemic levels when reconciling it with the newer registers.

| State | Content |
|---|---|
| Adopted author design choices | Largest commercial sky port on far-side summit at 5.4° N, 158.6° W; about 1.5 million regular occupants; 24 km tall, top about 35.6 km above sea level; trade/travel/hospitality/short stays, no homes or industry; round diagrid gathered into six splayed legs; lower rings, upper disks around a core from about 12 km; wide cantilevered lower parks; parks limited by air and sealed interiors above 15 km; regional docks at 10–15 km; berths throughout. |
| Adopted planning scale and character | Metropolis about 100 million people, about 40,000/km², radius about 32 km; metro/rail/large sky-ferry backbone with small sky boats and gliders; designated gliding and parkour spaces. These are choices, not demonstrated capacity or social uptake. |
| Explicit conditional assumptions | Fusion supplies most local power with solar/wind supplements; regional agriculture and foraging feed the metropolis. No fusion plant or food-area closure is supplied. |
| Author-adopted conditional findings | Similarity scaling; rigid ships reasonable to roughly 3 km, with 750–900 m favourable in model; small flyers winged/rotor and ferry-scale upward buoyant; personal wings/canopies feasible within the model; 700 kg four-seat sky boat; about 250 m/1,200-person ferries; urban traffic/berth constraints; first height stack; efficient fast winged liners. These remain model findings with stated calibration, atmosphere and design assumptions. |
| Explicitly open author choices | Crown long-haul class; whether winged liners also serve it; sky boats' trip share, ownership and traffic system. The 750 m liner recommendation is not a locked decision. |
| Cap status | Review completed; author has not adopted its recommendations. The branch awaits joint integration. |

The floor arrangement supersedes the prior full-ring proposal; preserve the 28 September decision rather than reviving the 27 September full-ring version. The new sea-appearance decision supersedes the inherited five-hour twilight canon. Neither update authorizes redesigning the metropolis, selecting a liner, abandoning fusion, imposing helium or fixing the magnetic routes during a merge.

## Scientific and practical value

1. **Engineering models are reusable outside the summit.** `engineering/towers/` provides closed-form lattice statics, Rayleigh frequency/buckling and wind devices. `engineering/fire/` provides Froude-scaled Earth fire correlations plus bracketed fires and evacuation. `engineering/flight/` provides calibrated airship and wing sizing, drag polars, canopies and vertical takeoff. These should survive intact as domain tools, with their limits.
2. **The tower separates geometry from feasibility claims.** The old square-lattice study supplies a size/shape trade; `summit_tower/form.py` supplies a hand-sized round diagrid with six legs, arches, lower cantilever rings and upper spoked disks. The latter has 39.3 km² of floor, 90 Mt steel and 57.8 Mt floors/parks in the September case. This is concept sizing, not a joint-, foundation-, construction- or dynamical-response design.
3. **The fire study materially limits architectural freedom.** Modelled flames grow taller; smoke entrainment falls; sprinklers answer to smaller fires; oxygen, pressure and wind matter together. It supplies a first compartment/refuge/water/shaft vocabulary. Its Earth correlations and small-sample partial-gravity extrapolation do not certify room-scale fire behaviour, crowd egress or aircraft fire protection.
4. **Flight scaling becomes a usable network question.** Fitting four historical rigid-airship weight statements (empty masses within roughly 7%) and ten aircraft wing groups (nine within roughly 12%) produces conditional size/energy envelopes. Fleet planning shows that cheap lift does not remove routing, berth, frequency, wake, hydrogen or rescue constraints. A 10% sky-boat trip share outruns the adopted Earth UAM spacing; the approximate 7–8% capacity and owned-versus-shared parking trade are useful planning constraints, not a traffic simulation.
5. **A well-developed world proposal is preserved without becoming evidence.** `immersion/world/summit/ring0/` is labelled informed/artistic, with vegetation placeholders. Its SVGs are intentionally committed design records. The separate visualization build outputs remain ignored.

## What the other branches change

### Main: corrected forcing, water and optics

The September studies read `A28_dim5` and the original `ring_ring.json`. Main corrected the Moon's obliquity and the Sun clock, then adopted `A28_dim5_moon` (years 20–29). The cap exports the corrected GCM winds but deliberately keeps September study products unchanged. It reports temporary-copy reruns with corrected wind products, `ring_ring_equator.json`, main's drainage and sunlight.

| Quantity | September committed study | Cap's temporary corrected-input rerun |
|---|---:|---:|
| Tower design gust | 27.3 m/s | 32.4 m/s |
| Square-lattice port steel; base | 32.8 Mt; 8.2 km | 46.0 Mt; 9.3 km |
| Chosen round-form steel; service sway | 90.0 Mt; 30 m | 98.0 Mt; 37 m |
| Whole-face screened slow rotors | 295 MW; 74% of port use | 575 MW; 140% |
| Moored ship design / service wind | 34.6 / 13.6 m/s | 38.6 / 24.4 m/s |
| 750 m liner nose mooring load | 0.34 MN | 0.43 MN |
| Pressure hull at 25 m/s | 350 m | 330 m |
| 100 / 200 m high platform station keeping | 15 / 62 kW | 115 / 460 kW |
| 200 m high-platform fabric load / reference strength | 0.39 | 1.2 |
| Sealed-zone stores-compartment burnout at median wind | 484 MW for 2.1 h | 818 MW for 1.2 h |
| Level-panel cycle mean | 33.0 W/m² | 32.8 W/m² |

These corrected reruns are reported historical measurements from the cap, not reruns conducted during this audit. Their numbers are valuable expected-regression targets, with rounding tolerance. The chosen form keeps its imposed 8,193 m base constant in `form.py`; it must not silently inherit the newly optimized square lattice's 9.3 km base.

The corrected upper air remains close in density and temperature, preserving much of the flight and park-by-height argument, while winds strengthen. The 32.4 m/s gust still comes from the first 10 km above the summit: corrected upper-port winds could imply 36–39 m/s under the same factors. Resolving that is a new height-dependent loading study, not a prerequisite for merging the corrected baseline honestly.

Main's waves and monthly tides add coastal design boundaries: typical nearside ranges about 3.7 m, 5.5–6 m in southern maria; mean nearside significant wave height about 1.39 m in the sampled forcing. Those are not Korolev-lake wave/seiche predictions. Main's paused CM1/GCM comfort disagreement remains a range and must not be erased by the metropolis brief's confident climate prose. No new GCM/CM1 runs should be hidden inside integration.

Main's optical studies supply improved low-Sun and shaded-scene illumination, useful for panels, urban daylight and parks. The study's fixed 22% panel assumption does not become a complete all-sky/high-altitude PV calculation.

### Atmospheric electricity: port and transport operating envelopes

The electricity cap's representative model is one equatorial, flat, sea-level 385 km box on 6 km cells through two lunar days. The infrastructure cap transfers its interpretation to the summit; it does not simulate the summit's terrain.

The striking intersection is height: the port crown is about 35.0–35.6 km above sea level and modelled flashes start near 35.8 km. A 24 km grounded structure can alter discharge initiation; it cannot be treated as a passive sample of background lightning. The cap reports 0.022 flashes/km²/year, roughly 20% to ground, with median 149 C and maximum 1,025 C ground transfers (larger under an alternative charging law). These are conditional model outputs, not site return periods or peak-current design loads.

The flight band is inside the strongest storms, whose reported tops reach 66–98 km and updrafts 12–36 m/s. Statements that the flight band simply avoids storms need qualification. The existing 14.3 m/s gust envelope and proposed crown hydrogen inventory near 2,100 t make storm avoidance, release/docking time, bonding, ignition control and wake behaviour urgent research interfaces. Ground ferries also need a storm service policy: a hydrogen sky ferry cannot be assumed to operate like metro/rail through the same events.

Infrastructure feeds back concrete test cases to electricity: a grounded 24 km conductor, long conducting aircraft, point/corona discharge from frame geometry, cables/pipes/rails, and hydrogen handling. Peak currents, leader onset and summit-specific storms remain open. The cap's Earth lightning-reference comparisons have not been independently source-verified in this audit.

### Sea appearance: usable twilight, darkness and night power

The adopted photometric definition places practical dusk near 2.98 lux, about 82 h after equatorial sunset and similarly before sunrise under the clear-sky model. At the summit the cap estimates roughly 190 h requiring lamps and about 116 h below 0.1 lux. Its elevation gives direct-Sun extensions, roughly 13 h at summit ground and 22–23 h at the crown, before terrain/refraction corrections.

This replaces “354 h of darkness” as an activity/lighting description, while a long low-solar-power interval remains. The cap explicitly retains full-night pumped-storage sizing: visual twilight is not usable panel energy. Sky boats need lit operations and gliders lose convective lift at night. The light calendar supplies phase/time capability; it does not settle culture, circadian health or a metropolitan schedule.

The reverse interface is also important: tall structures, urban surfaces and the lit shield stack add observers, horizons, reflectors and sources absent from a sea-level clear-sky calendar. Sea-level lux values must not silently become a summit photometry product.

### Solar shield / habitat array: surface–orbit interfaces

The shield cap proposes an enormous service industry in natural orbits and keeps energy-intensive industry/computing in orbit to control lunar heat. Its preliminary power links are orbit-to-orbit. They do not yet replace the author's fusion-led surface plan. A future ground power path must include atmospheric propagation, storms, receivers and local heat.

Infrastructure supplies a concrete 200 GW metropolitan demand scenario, roughly 80 W/m² locally, and high-altitude transport endpoints. It exposes the missing source-to-use network: surface/orbit freight, launch exhaust, capture, orbital receiving capability and allocation between local and orbital industry. S7 applies to protection-system outflows; extending it to ordinary launches is a recommendation for author review.

Regional magnets add a surface grid, long routes, fields, mechanical forces and night refrigeration. The cap reports preliminary four-loop cases with 9–15 GW refrigeration and 12,600–25,100 km route totals; these are screening designs. Its straight-wire estimates suggest land-use and instrument corridors, but are not a settled corridor map or a medical safety assessment. The summit is outside tested routes, while potentially comfortable polar settlement regions intersect them.

The shield's lit ring fleet may change far-side darkness and astronomy-platform usefulness. Hydrogen inventory estimates show the sky fleet is small against the assumed atmospheric H2 inventory, but do not close the leakage/source/escape budget. The shield's ultraviolet edge and dimming spectrum can alter surface panels; a visible appearance study cannot settle that engineering choice.

## Intersections actually resolved, partially resolved, and open

| Intersection | What is resolved now | What remains |
|---|---|---|
| New climate versus adopted port | Most air-by-height and rigid-hull scaling survive; wind-sensitive terms are identified and temporarily recalculated | Updated committed products, upper-height gusts, terrain/wake/dynamics |
| Night versus solar supply | Longer usable twilight is compatible with long energy storage need | Summit terrain/cloud/elevation photometry, lamps, heat and operations |
| Flight band versus storms | The band is not wholly above storms; electrical and wind hazards overlap it | Summit storms, conductor/aircraft triggering, peak currents and operational rules |
| Port form versus fire study | Both can be retained as separately labelled concept studies | Fire geometry still belongs to old narrow bands, not chosen wide rings/disks |
| Cheap personal flight versus metropolis scale | The adopted planning case exposes airspace and parking constraints | Network simulation, access, ownership, fail-safe procedures and service in storms/night |
| Orbital power versus surface fusion | Current calculations address different receivers and can coexist | Surface power link, allocation, heat rejection and author selection |
| Ring/array industry versus lunar city | Service chain and demand endpoints now identifiable | Transport, emissions, materials/energy ledger and public artificial gravity provision |
| Seas versus highland metropolis | Main supplies separate coastal wave/tide boundary conditions | Korolev storage interaction, lake waves/seiches, water/food/heat closure |

## Concrete integration hazards

1. **Textual hubs are cumulative ledgers.** Merge additions in `research/README.md`, `research/status.json`, `research/decisions.md`, `climate/gcm/README.md`, `habitation/README.md`, `engineering/README.md` and `visualization/README.md` semantically. Retain main's corrected climate, optical comfort, waves and pauses; add infrastructure capabilities without restoring the old global status. Parent's dry-tree conflict survey is the authoritative textual-conflict inventory.
2. **New corrected wind files are not yet consumers' defaults.** Update `summit_tower/run.py`, `summit_tower/form.py`, `sky_ships/run.py` (and inherited aliases in `sky_fleet/run.py`) to explicit `_moon` products and `ring_ring_equator.json`. Also update `visualization/sky-fleet/figures.py`, which directly names `ring_ring.json`. Do not repeat the cap experiment's temporary renaming trick in the final producer: current paths and recorded hashes must tell the truth.
3. **Different port geometries currently coexist.** `summit_tower/run.py` carries a square lattice and dense light bands; `form.py` fixes a round-form geometry with 39.3 km²; `port_fire/run.py` explicitly assumes 40 m-deep strips, four-storey bands every 200 m and its own compartment layouts. `visualization/summit-port/` still displays an older hexagonal/wing scheme and knowingly noncredible slender crown arms. `visualization/sky-fleet/crown.py` uses the newer form. Preserve historical studies/proposals with clear current-versus-superseded roles; do not label the old fire study certification for the chosen design or silently redesign it during integration.
4. **A specific geometry propagation trap.** The corrected square lattice's base becomes about 9.3 km; `form.py` fixes `B0=8193.0`. `sky_fleet/run.py` reconstructs crown width from the square-lattice product while also reading terminal/disks from the form product. Audit that mixed-geometry crown layout before promoting regenerated diagrams. The branch's ring-0 artistic data fixes a radius of 3,334 m; it should remain an explicit proposal tied to a form version rather than mutate opportunistically.
5. **Stale status claims inside the branch.** The decisions table and metropolis brief say legs/rings/disks/core have yet to be sized even though `aac8c64` adds their hand sizing. Sky-ships' “next distribution” text predates sky-fleet. Main status also remains unaware of all this. Reconcile status to “concept sized; structural analysis and programme reconciliation remain” without declaring the structures designed.
6. **One actual stale producer hash at branch head.** `sky_ships/results/sky_ships.json` records an `engineering/flight/winged.py` hash from before `1432efa` appended glider/VTOL functionality. This audit found the mismatch. The diff appends functionality rather than changing the previously used equations, but that observation does not justify rewriting the stored hash without a rerun. All 32 other checked available files/input hashes across the five study products matched.
7. **Provenance is incomplete even where listed hashes match.** These products have no named-constant records. `sky_fleet/run.py` directly reads `MOON_RADIUS` and `SYNODIC_MONTH_DAYS`; use main's `shared.provenance` practice. Indirect thermodynamics and cross-study imports are not always pinned (`port_fire` imports tower assumptions; `form` imports lattice); ensure actual dependencies are recorded. Wind exporters record code/configuration/run/year range but not source NetCDF hashes. Consumer schema checks are partial; strengthen relevant product boundaries while touching them. Do not hash the entire constants table and force unrelated recomputation.
8. **Literal constants versus scenario assumptions.** The branch still has `G_EARTH=9.80665`, `NIGHT_HOURS=354.4` and assorted unit/reference values. Review true shared physical constants against main's named constants; keep stated design choices/calibration data distinct. Do not change established GM/R² versus legacy 1.62 conventions in passing or rewrite pinned historical baselines.
9. **Ignored data are essential but not a reason to rerun climate.** The lightweight tower runner needs two ignored geography grids. Main's corrected drainage grid differs from the infrastructure worktree's older grid despite the same filename. Restore from the selected baseline with hashes and metadata; never copy the infrastructure drainage accidentally or rename old products to masquerade as corrected ones.
10. **Generated artifacts require separate handling.** Versioned study JSON must be regenerated after producer changes; ignored visualization builds should be rebuilt only if included in acceptance. Ring-0 SVGs are intentional world design records; preserve them pending explicit revisions. Unreal work mentioned under `/media/projectspace/tmp` is outside this Git branch and is not delivered by a merge.
11. **Locked and historical files remain untouched.** Infrastructure has no manuscript changes; no full editorial review is needed to merge its engineering lane, and no manuscript evidence admission follows from successful tests. Preserve five paper roles/UUIDs, locked planning/seeds and `research/provenance.json` byte pins.

## Exact integration inputs and order

No expensive model rerun is needed to reproduce the cap's corrected integration baseline. The corrected wind JSONs are committed on infrastructure; the corrected ring, sunlight and drainage summary are committed on main. The atlas and drainage binary products exist in main's worktree and are ignored by Git.

| Offline input to copy into integration at the same relative path | Current main size | SHA256 measured in this audit |
|---|---:|---|
| `geography/products/atlas_28pct_4ppd.npz` | 9,875,080 bytes | `dc283827f349fdc134816e1a55deef3d608047fba9fbb7521a6005c171d81693` |
| `geography/products/drainage_28pct_16ppd.npz` | 2,697,170 bytes | `98842fc2a648f8899bd4f5d712445f78c8882c85a1bc30bf631ac9dd14ab1c03` |

Source root: `/home/kir/Documents/Projects/terluna/`. Copy only these required artifacts, preserving metadata and checking hashes; a writable symlink to live source outputs would undermine isolation. The atlas metadata is schema `terluna.geography.atlas-grid/1`; drainage metadata is `terluna.geography.drainage/1`. The branch's old drainage is 3,066,002 bytes, SHA256 `c816483f7ce97f833f88b301705a3d15d36caed07d4726a85d21f1beace12251`, and is the wrong input for the corrected run.

Required tracked inputs after source merges:

- `climate/results/gcm/site_winds_A28_dim5_moon_summit.json`
- `climate/results/gcm/global_winds_A28_dim5_moon.json`
- `climate/results/crm/ring_ring_equator.json`
- `geography/results/drainage.json` from main, matching the selected drainage grid
- `illumination/surface_light/results/surface_light.json` from main, retaining its declared solver/spectral limits
- `research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json`

After correcting explicit input paths, provenance and labels, execute in this dependency order from the unified candidate root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.summit_tower.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.summit_tower.form
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.port_fire.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.sky_ships.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m research.studies.sky_fleet.run
```

These are four studies but five commands, because the tower produces both the square-lattice comparison and chosen-form products. The cap reports seconds per study on one core; this audit did not rerun them. `port_fire` and `sky_ships` can follow tower independently, but serial execution is simpler and considerate of shared resources. Fleet needs both ships and form. Update prose from regenerated products and compare against the cap's expected shifts; do not blindly replace the author-chosen 1.5 million with the square-lattice programme's numerical drift.

If correcting/exporting wind-product provenance requires exporter regeneration, the required ignored inputs are `climate/gcm/runs/A28_dim5_moon/progress.json` and `climate/gcm/runs/A28_dim5_moon/model/MOST.00020.nc` through `MOST.00029.nc`. These exist in main's GCM-run storage (this audit checked the progress file, not all ten NetCDF files). Commands, after making those inputs available read-only in the isolated candidate and ensuring the environment:

```sh
/path/to/main/climate/gcm/.venv/bin/python -m climate.gcm.site_winds A28_dim5_moon:20-29 summit 5.375 201.375
/path/to/main/climate/gcm/.venv/bin/python -m climate.gcm.global_winds A28_dim5_moon:20-29
```

Those commands postprocess existing results; they are not climate simulations. Raw atlas/drainage regeneration would need geography inputs/cache and corrected climatology; avoid that unrelated resource-heavy route when verified products already exist. Preserve failed-restoration states explicitly.

## Verification and acceptance

A scoped check was actually run in the clean infrastructure worktree with bytecode and pytest-cache writes disabled and numerical thread pools limited:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -p no:cacheprovider \
  engineering/tests/test_towers.py engineering/tests/test_fire.py engineering/tests/test_flight.py \
  climate/tests/test_global_winds.py climate/tests/test_site_winds.py \
  research/studies/summit_tower/test_run.py research/studies/summit_tower/test_form.py \
  research/studies/port_fire/test_run.py research/studies/sky_ships/test_run.py \
  research/studies/sky_fleet/test_run.py
```

Result: **89 passed in 0.66 s**. Tests include algebraic/benchmark checks, calibration and stored-product relationships. They do not discover the sky-ships stale producer hash, certify current inputs globally, validate lunar physics, or prove that the fire compartment geometry matches the chosen port form.

Also actually performed: Git ancestry/diff/log inspection; source/README/decision comparison; a 33-dependency hash audit of the five study JSONs (one mismatch); read-only main/infrastructure ignored-grid presence/size/SHA256 comparison; metadata inspection; final source-worktree cleanliness check.

The cap separately reports baseline temporary-copy reproduction and corrected-input study reruns on 8 October; distinguish those historical checks from this audit. No `make check`, full-repository pytest, graphics build, numerical model run, external source reading, scholarly admission, manuscript edit or exact-build PDF review was performed by this sub-review.

Before accepting the unified implementation, run the same focused tests against regenerated products, add meaningful schema/dependency freshness coverage for touched products, inspect geometric consistency of the crown and its diagrams, and run `make check` as the repository requires before any eventual commit. Record the protection-input state separately (`research/check.py --require-inputs` only if declaring complete spectral-input restoration). Keep the unimplemented summit CFD, leader/charge interaction, height-dependent gust model, chosen-form fire model, hydrogen safety and traffic/heat/food closure as named future research rather than merge prerequisites or completed work.

## Recommended follow-on questions, after integration

The first coupled research target is the crown as a grounded conductor with hydrogen berths: it joins the most consequential new evidence from electricity, tower geometry, wind, flight and fire. A bounded electrostatic calculation using existing charge fields can precede a terrain-specific electrified storm simulation. The next engineering target is a common chosen-form geometry product for loads, programme, fire/refuges and berths; it prevents independent models from quietly analysing different buildings.

In parallel, the surface–orbit interface deserves its own material/energy/outflow ledger, while the night calendar should drive explicit operations and energy cases. High-platform height/drift, Korolev lake storage/waves, city heat and food area then become concrete design trades. None requires reopening the paused climate programme merely to assemble a truthful unified branch.
