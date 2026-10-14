# What binds the array: a review of the tile, support and flexible work

14 October 2026 · `research/array-industry` · a review of commits `e1196f5`, `669ddd3` and `25df091` and of the branch
since `fadb0ab`

This review reads the branch's three newest commits, the [tile hardware](tile_design.md),
[materials and supports](support_comparison.md) and [flexible coupling](flexible_coupling.md), against the ring fleet as
designed: about 28.2 million tiles 10 km square in 1,859 nested rings near 20,000 km, a shadow of four lunar radii,
station keeping by photons and the night roll. It sets beside them the climate's corrected radiation on
`domain/ecology`, which reaches this branch at the next integration. Three kinds of number appear below: *product*
values from committed products, *ecology branch* values from `domain/ecology`'s products and run logs, and
*estimates*, this review's arithmetic on those, each with its basis. The estimates are reproduced by `estimates.py`
and `drive_check.py` in `/media/projectspace/terluna-research/array-industry-review-2026-10-14/`.

The three commits are careful work. Every headline number in their notes matches its product, the physics checks out
where this review tested it (the membrane's sag series, the frame and truss sizing, the conserved momentum, the
Hill–Clohessy–Wiltshire faults, the service queues, the curved-edge geometry and corner thrust, the resolved normal
acceleration and the floating-frame mechanics), and two of their findings about the inherited design are right and
consequential: the attitude study's momentum store leaves out the transport of the stored momentum, and the window
tile's allocation is its bare optical stack. The findings below are ranked by what they change in the fleet.

## What binds the array now

The tiles' areal mass and the photons' authority bind together. One sail loading ties every ring's photon acceleration
to the window tiles', which the dynamics take as 0.276 over 26 g/m²: all of the window's turned-away light carried as
specular reflection along the normal, and the window tile at the mass of its bare stack. Against that loading the
annulus tile is allowed 13.765 g/m², 3.55 g/m² above its film. No tile architecture tested fits that, and the photon
authority that would let the tiles grow heavier is closer to its limit than the records state, because the films
reflect only part of the light they take from the beam and only reflected light steers.

| Limit | Where it stands (*product*) | *Estimate* | Product |
|---|---|---|---|
| Annulus tile mass | Allowance 13.765 g/m², film 10.211 g/m². Rigid frame 2,684 g/m². Lightest film-plus-support subtotal passing the 150 m gate: 139–345 g/m² at the flown 1.5° shingle, 29–37 at 2°, 21.5–23.6 at 3° with 5 km tiles in 4.3 times the count | Matched on the eccentricity drive, 16.6 g/m²; 18.9 with the corrected climate window | glow.json; tile_design.json; geometry_screen.json and its archive |
| Window tile mass | Allocation 26 g/m², equal to its stack's 26.921 g/m² | Window structure lowers the common loading in proportion: 3.6 g/m² takes 12% | annulus_film.json; tile_design.json |
| Plane steering of the edge rings, rolled layout | At 1 km between strips, 0.86 of the photon roll's modelled reach in the layout and 0.95 on its sampled rings flown; at 0.6 km, 0.67 | 1.1–1.45 of the reach the annulus films give; 1.23–1.65 with the corrected window; 0.85–1.17 at 0.6 km | ring_layout_rolled.json |
| Translation keeping | 1.58 times the kept run's interior demand at 15,000 km with 62.7 g/m² tiles (0.89 at the edges), at 2° tilts and 25% trim | About 3.8 (2.1 at the edges) at the design loading; about 2.5 (1.4) with only reflected light steering; the 20,000 km demand is unflown | photon_control.json |
| Momentum store for the roll | 1.0–2.3×10¹⁰ N·m·s, 5–12 times the chosen store, from an integral without transport | Rotors of 1.1–2.6 g/m² at 400 m radius, 4.5–10.4 at 100 m (three axes, 50% supports) | attitude_schemes.json; tile_design.json |
| Shingle clearance | 222.5 m between parallel shingles against the 150 m gate, so 72.5 m for both envelopes | | tile_design.json |
| Service | 5,153 planned exchanges and 773 failures a day against 6,184 dispatches, 4.35% spare | | tile_design.json |

## Findings, ranked by consequence

### 1. The force model carries all of the window's turned-away light as reflected, and only reflected light steers

*Where:* attitude_schemes.py:238–253 (`load()`, the pressure at :247) and cycling_search.py:23 (the redirected fraction,
1 − 1,173.25/1,361 = 13.795%); ring_layout.py:187–192 (`reach()`, the photon roll's steering reach for every ring);
frozen_rings.py:101 and ring_screen.py:121 (the same model for the forced eccentricity).

*What:* the dynamics model every tile as a specular reflector of reflectance 13.795%, a force coefficient of 0.276 at
normal incidence, at the window's 26 g/m². The window stack reflects 4.69% of the sunlight and absorbs 4.49%, and the
5% dimmer takes 4.54% more (annulus_film.json `window_stack`); the annulus film reflects 5.79% and absorbs 3.03%
(glow.json). Absorbed light pushes along the sunlight and reflected light along the tile's normal.

*Why it matters:* for a radial-facing tile the orbit-averaged drive on the eccentricity across the Sun line is
proportional to 2(R + A), the light taken out of the beam, because the absorbed light's push along the track adds to
its radial push (*estimate*, Gauss's equations averaged over the lit orbit and checked numerically in
`drive_check.py`). The convention therefore gets the window's drive right, 0.274 for the stack with a reflecting
dimmer against 0.276. A tilt or a roll turns only the reflected force, so the side force that steers a plane goes as
2R: 0.185 for the window and 0.116 for the annulus film. Per unit of drive, the steering and keeping authority is 0.67
of the modelled for window tiles and 0.66 for annulus tiles matched on the drive (finding 2); 0.79 for the annulus at
its present 13.765 g/m².

*Consequence:* the rings that need the most steering are annulus rings near the stack's edges, beyond the window's
middle 26%. They need 0.86–0.95 of the modelled reach at 1 km between strips, so about **1.1–1.45 of the reach their
films give** (*estimate*), and the rolled layout's statement that 1 km stays inside the photon roll's reach rests on
the convention. At 0.6 km between strips the coupled matched layout needs 0.67 of the modelled reach, 0.85–1.03 of the
actual. The translation keeping keeps a margin (finding 6).

*Fix:* fly the steering roll and the frozen rings with each film's own reflectance and absorptance, the window with its
dimmer in the form chosen for it, and lay the stack out again (recommendation 3).

### 2. One sail loading was matched on the normal-incidence coefficient

*Where:* glow.py:1155–1174 (`SAIL`, `sail_loading()`); integrated_comparison.md:384–389.

*What:* the match sets the annulus tile so that its (2R + A)/σ equals the window's 0.276/26 g/m², which gives
13.765 g/m². Nested rings stay nested when their eccentricity drives match, and for radial-facing tiles the drive goes
as 2(R + A)/σ (finding 1). With the annulus film's 2(R + A) = 0.1765 the drive-matched annulus tile is **16.6 g/m²**
(*estimate*); at 13.765 g/m² its rings carry about 21% more drive than the window rings and settle on forced
eccentricities about 21% larger, which the integrated comparison found makes nested rings cross.

*Fix:* match one sail loading on the orbit-averaged drive under the law flown, the 80° roll, with each film's optics,
in frozen_rings (recommendation 3). The room for structure above the annulus film would then be about 6.4 g/m² in place
of 3.55.

### 3. At the flown 1.5° shingle the lightest passing support weighs 139–345 g/m²

*Where:* the geometry screen's archive
`research/runs/solar_shield_array_geometry_screen/configurations.jsonl`, hashed in geometry_screen.json, and
flexible_coupling.md's first table.

*What:* the lightest film-plus-support subtotals passing the member and 150 m clearance screens, 21.5–23.6
g/reference m², are 5 km tiles at a 3° shingle in 4.3 times the count. By shingle angle the lightest passing subtotal
is 138.6–345.0 g/m² at 1.5° (24 of 5,400 configurations pass; 10 km, 2.5% retreat, 0.3% boom depth), 29.0–37.2 at 2°
(518 pass) and 21.5–23.6 at 3° (1,862 pass), from a quarter of the inherited load to all of it (*product*).

*Why:* within a ring the parallel shingles stand 8,500 m × sin 1.5° = 222.5 m apart, so the 150 m gate leaves two
supports 72.5 m between them, and booms 7 km long held to tens of metres of depth are heavy. A 3° shingle with half-size
tiles keeps 222 m with half the span. A larger shingle departs from the flown design: the attitude study's clearance
envelope, shadow torques, store and push are all for 1.5°.

*Fix:* flexible_coupling.md now gives the table by shingle angle; the shingle angle and tile size belong in the next
hardware screen together with the attitude study (recommendation 6).

### 4. The roll's momentum store comes from a body-axis integral without transport

*Where:* attitude_schemes.py:264–274 (`store()`, the update at :272); the tile study found it (tile_design.md).

*What:* `store()` integrates the required torque component by component in the rotating body axes and leaves out the
term −ω × h. The body turns once an orbit about the orbit normal, and the roll turns its axes further; normal-axis
momentum, which no trim unloads in `store()`, would in fact turn into the in-plane axes, where trim can unload it, and
in-plane momentum would turn too. The roll's 1.0–2.3×10¹⁰ N·m·s (5–12 times the chosen store), the night roll's
2.7–9.0×10¹⁰, the node tilt's 1.2–1.7×10¹⁰ and the chosen 1.9×10⁹ all come from this integral, and the glow note, the
README and the register's glow row carry them. The tile study's conserved calculation (array_momentum.py) adds the
transport and leaves out the radiation torque, so neither gives the roll's store.

The same product gives smaller rolls larger stores: 5.8×10¹⁰ N·m·s at 30°, 6.2×10¹⁰ at 60° and 2.3×10¹⁰ at 80°
(62.7 g/m², 25% trim). Whether that survives the transport term is open.

*Fix:* add the transport term to `store()` and fly attitude_schemes again, about an hour of CPU along the kept orbits
(recommendation 2). The register's row reports the product correctly and stands until then.

### 5. The window tile's allocation is its bare stack

*Where:* annulus_film.json `window_stack.areal_mass_g_m2` (26.921); annulus_film.py:7; tile_design.json
`window_bare_film_excess_over_complete_allocation_g_m2` (0.921). The tile study found it.

*What:* the one-sail-loading match took the window tile at the climate stack's mass, so the window tiles carry no
structure, store, trim or power in it. Any window structure lowers the common loading in proportion, 12% for the
3.6 g/m² the annulus has, and with it the steering margins of finding 1.

*Fix:* add the window's structure to its allocation in the match (recommendation 3).

### 6. The photon keeping margin of 1.6 belongs to the 15,000 km kept run with 62.7 g/m² tiles

*Where:* photon_control.json (`tile_areal_mass_kg_m2` 0.0627, `period_h` 45.79, `translation`); photon_control.py:21;
README:68–72 and integrated_comparison.md:65 and :197 (before this review the README gave the 1.6 without its
conditions).

*What:* at the design loading the reach is 2.41 times that run's, the product's own `sail_reach_factor`, so about 3.8
times its interior demand and 2.1 at the edges with 2° tilts, 2.7 with 5° (*estimate*); with only reflected light
steering (finding 1) about 2.5 and 1.4. The demand at 20,000 km, with its faster tidal regression and the two-sided
layout, has not been flown; the kept run's demand is inflated by its fixed keeping frame.

*Fix:* the README now states the run's radius and tile mass; fly the keeping at 20,000 km on the rolled layout
(recommendation 4).

### 7. The roll loads the membranes more than sunlight does

*Where:* tile_design.json `annulus_reference.membrane`; support_comparison.json `loading`.

*What:* on the rolled command the rotation is 79% of the annulus tile's membrane load envelope, and the resolved normal
acceleration at the corners, 7.91×10⁻⁵ m/s², is 1.64 times the sail's own acceleration on the annulus tile
(*product*, and *estimate* for the ratio). The 30° cosine ramps of the roll set it.

*Consequence:* a slower roll (wider ramps) lowers the membrane tension, and with it the frame or support mass in
proportion, at a cost in the light the ramps let onto the night's edge; it belongs in the roll's trade
(recommendation 2).

### 8. The rotors' mass rests on the inherited store and the rotor proportions

*Where:* tile_design.json `hardware.variants` (`annulus_side_1000m`, `annulus_side_3000m`); glow.md:392–393.

*What:* the tile study sizes each axis with the full capacity, three axes, 50% supports and rims at 1,000 m/s. At 3 km
its rotors and supports are 10.3 g/m², sized for the heavy rigid tile's inertia; at 1 km they are 6.7 g/m², the
inherited 1.498×10¹⁰ N·m·s scaled as L³, about eight times what the conserved calculation asks (*product*). With rotors
of L/100 in radius the inherited allocation costs 6.7 g/m² at every tile size.

*Estimate:* on the same model the roll's 1.0–2.3×10¹⁰ N·m·s needs 1.1–2.6 g/m² of rotors and supports at 400 m radius,
a third to three-quarters of the annulus's 3.55 g/m² of room, or 4.5–10.4 g/m² at 100 m; the chosen 1.9×10⁹ needs 0.21
or 0.85 g/m². The glow note's 0.7–1.6 Gt is one axis of rims; three full axes with supports give 3.2–7.3 Gt, 7–16% of
the fleet's 46 Gt.

### 9. The flexible model products carried no reading rule

*Where:* flexible_benchmark.py:155–178 (`export_model()`) and :185–216 (`build_models()`); the products
results/flexible_models.json and its embedded `terluna.engineering.array-flexible/1` models.

*What:* AGENTS.md asks a product that crosses lanes for a versioned schema, the producer and model hash, an evidence
statement and the rule for reading it; these had the first three.

*Fix (done):* the producer writes a reading rule for the collection and for every model, test_flexible_benchmark.py
asserts it, and the products were regenerated (below). Two schema names remain for related objects: the model's
metadata is `terluna.engineering.flexible-array-model/1` (array_flexible.py:262) and the product
`terluna.engineering.array-flexible/1` (flexible_benchmark.py:173); the saved-state audit's script sits in the ignored
`research/runs/` with absolute paths in its record.

### 10. The source registers: every number right, five statements misattributed

The registers' 36 sources were each fetched and checked claim by claim. Every density, modulus, strength, dimension and
date matches its source, and the derived densities and admissible stresses follow from them. Corrected in the
registers (text only, every value unchanged):

- support_materials.json, `beryllium_s200fh`: the data sheet points users to the safety data sheet and the exposure
  limits; the recycling requirement attributed to the supplier is in none of the cited sources.
- support_materials.json, `cnt_ds_2026`: the 14 GPa high-strain-rate figure belongs to a different paper; the
  limitation now covers high-strain-rate strengths in general, and the density derivation states that the paper's
  7.03×10⁶ N·m/kg was measured before rolling, which it reports changed it little.
- support_sources.json, `solar_cruiser_ground_2024`: NASA funded the technology to reach TRL 6, and the article quotes
  the project saying it met the goal; the register said NASA reported reaching it.
- support_sources.json, `heliogyro_ground_2015`: the related 2014 paper is now named (NTRS 20140003899), with the 2015
  paper's change to a 0.46 m rope-ladder blade.
- tile_sources.json, `nasa_flywheel_viscoelasticity_2004`: the report says the relaxation may limit rotor durability,
  with its measured 16.4% and 13% relaxation; the register said it shows that it does.
- Conditions restored: Toray T700S's 90° tensile strength (81 MPa) is a raw value, outside the sheet's 60%
  fibre-volume normalization; Hexcel 8552's values are representative at a nominal 57.42% fibre volume; AGY's sheet
  states no test basis; the material catalogue's general limit names specified minimums as well as typical values.
  Bibliographic notes added for CABLESSail (its version and authors), the 2005 sail test (authors, AIAA 2005-2121),
  the Kevlar guide's editions and the solar-sail thrust model's access.

### 11. Some tests are weaker than their names

*Where:* engineering/test_array_flexible.py:104–113 (the internal force is defined as the energy's gradient, so the
check holds for any symmetric stiffness), :122–134 (the ordering of modal work holds by construction for nested
modes), :76–78 (the inertia check recomputes the producer's formula); protection/dynamics/test_flexible_fleet.py:61–67
and :212–220 (energy minus work closes for any force law, the photon and tidal signs being tested in
test_flexible_optics.py and test_flexible_environment.py). The coupling block G of the mass matrix vanishes for every
benchmark's normal-only modes, so only the random three-dimensional modes of test_flexible_fleet_energy.py exercise it;
flexible_coupling.md now says so. The modules themselves are correct where checked (flexible_fleet.py:180–204 follows
Kane's equations; flexible_optics.py:339–345 has the right signs and angular dependence).

### 12. Design values are repeated as literals in the new producers

*Where:* tile_design.py:78, :175, :257, :271; support_comparison.py:54–64, :115–116, :209; flexible_benchmark.py:451.

*What:* the 8,500 m pitch, the 1.5° shingle, the 8,900 × 9,890 m optical rectangle, the 9,286.456 m height step and the
150 m gate are written out instead of read from `ring_bundle_run.DESIGN`, `ring_bundle.SIDE` and `CLEAR`,
`attitude_schemes.DESIGN` and the layout product. These are design values, outside shared/constants.json, so AGENTS.md's
rule on constants is met, but a change of design would leave these producers behind. *Fix:* read them from their
designs at the next change of these producers.

### 13. The status records an author direction the register lacks

*Where:* status.json, the engineering and protection `next_session_task`.

The tile work's session records that the author asked on 14 October for the cheap screen and the coupled model's
implementation, verification and benchmarks, with the coupled search after them. The register holds no such row; the
status now says so, for the author to confirm.

### 14. The glow note's next steps and efficacy, for the next integration

*Where:* glow.md:104–105, :409–414 and :439.

The glow note gives the annulus absorber and the ultraviolet mirror as the next steps and answers infrared at an
efficacy of 1. On 13 October, as `domain/ecology` records, the author set those two answers aside, with the plants'
water design making up the rain the glow takes, and gave the glow's dimming the GCM's 0.71 of its infrared's watts. The
note follows those decisions when the branches are integrated.

## The climate side, on domain/ecology

These results are on `domain/ecology` (commit `55e5020` and the GCM's runs), not yet integrated here.

- On 14 October the author made the GCM's corrected radiation its standard: the Moon's curved air for sunlight and
  infrared, CO2 and O2 in PlaSim's sunlight, and the water and Rayleigh multipliers refitted with them. The cooling it
  finds is answered with less dimming, holding the design temperature (*ecology branch*, register).
- At the old design dimming, a sunlight scale of 0.95 on the titania stack's 1,235 W/m² (1,173.25 W/m² at the top of the
  air), the corrected model settles near 292.0 K over years 10–19, still cooling 0.3 K a decade, 2.9 K under the design
  case's 294.90 K, with 7.6% less rain (*ecology branch*, climate/gcm/README.md).
- At a scale of 0.9655 it settles at 294.40 K over years 10–19, with 2.25 mm/day of rain against the design case's 2.24
  (*ecology branch*, run log `gcm-runs/A28_corrected_9655.log` on the research drive).
- A run at 0.9685, aimed at 294.9 K, is in progress. From the 0.9655 run and the dimming runs' 1.4–1.9 K per 1% of
  sunlight, the scale that lands on 294.9 K is 0.968–0.969, a dimming of about 3.2% below the shield product in place
  of 5% (*estimate*, pending). The run's years 10–14 average 294.75 K (*ecology branch*, run log
  `gcm-runs/A28_corrected_9685.log`, at 19:16 on 14 October), which at the two runs' 1.1 K per 1% would put the scale
  near 0.970 and the dimming nearer 3.0%; its years 10–19 decide it.
- So the climate window would pass about 1,196 W/m² in place of 1,173.25, and the window tiles would turn away about
  12.1% of the 1,361 W/m² in place of 13.795% (*estimate*, pending; at a scale of 0.970, 1,198 W/m² and 12.0%). The
  trades below take the scale of 0.9685.
- The added dimming for the gases and the glow, 1.69% and 1.80% plus 0.21% of the design sunlight, was sized on the old
  radiation in `domain/ecology`'s gases_glow_dimming study and has not been measured again.

## The trades

### Light pressure and one sail loading (*estimates*)

With the window turning away 12.1%, the convention's coefficient falls from 0.276 to 0.242. With the window at 26 g/m²,
one sail loading then lets the annulus tile weigh 15.7 g/m² in place of 13.765 (5.5 g/m² above its film in place of
3.55), or 18.9 g/m² matched on the drive (8.7 above its film). The common photon acceleration falls 12%, and with it the
sail-force disturbances, the forced eccentricities and the photon authority: the edge rings' steering need rises to
1.23–1.65 of the reach their films give at 1 km between strips, and 0.97–1.17 at 0.6 km. The corrected window gives the
annulus room for structure and takes photon authority, and the edge rings' steering decides how much of the room can be
used.

### The roll (*estimates*)

As built the stack puts 14.12 W/m² on the Moon (glow 9.37, the Moon's infrared returned 4.04, its sunlight returned
0.71); rolled 80° outside the service arc on both halves, on the layout re-flown for the roll's push, 5.95 W/m² (3.60,
1.75 and 0.61) (*product*, glow.json). Dropping the roll adds 8.07 W/m² of infrared and 0.10 W/m² of sunlight. On the
Moon's absorbed design sunlight, 0.72 × 293.3 = 211.2 W/m², that needs 2.75% more dimming at the ecology GCM's 0.71
watts of dimming per watt of the glow's infrared, or 3.85% at this branch's efficacy of 1.

The corrected radiation frees 0.9685/0.95 − 1 = 1.95% of the design sunlight, 71% of the roll's heat at 0.71. In total,
the fleet without the roll under the corrected radiation needs 1,139 W/m² at the top of the air, against 1,150 W/m² for
the old plan with the roll: 0.9% more dimming, so a much smaller roll could keep the planned temperature on heat alone.

- **The rain.** The extra glow takes about 3.7% more of the rain, at 0.46% for each W/m² of the glow's infrared, beside
  the rolled glow's 2.4%. Under the author's decision of 13 October the plants' water design makes up the glow's rain:
  3.2% more carbon fixed per litre for the rolled glow, about 8% in all without the roll.
- **The night.** The author's night decision keeps the fleet's light off the night zones. As built, the films' specular
  return of the day side lights the night hemisphere at 1.6 lux (1.0 lux at midnight, 2.6 within 10° of the
  terminator), a thousand times a far-side midnight, beside the joint synthesis's 25–66 lux of diffuse glow at O5's
  0.1% scatter: an unrolled night half breaks the decision. Rolls of 30° and 60° keep the specular return at 0.08 lux
  on the as-built eccentricity, but their scatter has not been computed below 80°, where the joint synthesis brought a
  far-side midnight to 0.02 lux.
- **The store.** The attitude product gives the smaller rolls larger stores (finding 4), so a smaller roll saves heat
  without saving store until the store is recomputed with its transport term.

### Heavier tiles (*estimates*)

Tiles made uniformly heavier on every ring keep one sail loading and lose photon authority in proportion. The
translation keeping would allow about 3.8 times the design loading's masses on the 15,000 km run's demand (2.5 with
only reflected light steering), but the edge rings' steering binds first: at 1 km between strips it allows 1.05–1.16
times as modelled and already falls short with the films' reflected share, while at 0.6 km it allows about 1.5 times as
modelled and 1.0–1.2 with the reflected share. The flexible screen's 21.5–23.6 g/m² is for 5 km tiles at a 3°
shingle, and at the flown 1.5° shingle no passing support is lighter than 139 g/m², so a heavier tile that still flies
needs a different shingle, closer strips or more reflective films.

### Rotors (*estimates*)

At the roll's store as computed, rotors and supports of 400 m radius take a third to three-quarters of the annulus's
room, and those of 100 m all of it (finding 8). The chosen radial-facing store takes 0.2–0.9 g/m². The roll's amplitude
therefore weighs in rotor mass as well as in heat and night light, once its store is computed with the transport term.

## Recommendations

For the author; each says what it would change from the known design.

1. **On domain/ecology, fly the glow pair again under the corrected radiation**, to measure the 0.71 factor and the
   glow's rain cost anew: the CO2 and O2 bands already thin the high cirrus that the glow was thinning. It would change
   the glow's dimming (2.0% of the design sunlight for the rolled law at 0.71) and the plants' water design.
2. **Recompute the roll's momentum store with the transport term and the radiation torques**, then price the roll's
   amplitude and ramps against rotor mass, membrane load, the Moon's heat and the night's scatter. It would replace the
   5–12 and 14–48 times the chosen store in the register's glow row, the glow note and the README, and could move the
   recommended roll away from 80° on both halves.
3. **Fly the steering roll, the frozen rings and one sail loading with each film's own reflectance and absorptance**,
   the window tile with its structure and the climate window's corrected redirected fraction (0.121 at a scale of
   0.9685), matching the loading on the eccentricity drive under the rolled law, and lay the stack out again. It would
   change the annulus allowance (from 13.765 g/m² toward 16.6–18.9), the forced eccentricities and radii, and the edge
   rings' steering margin; if the edge rings fall beyond the roll's reach, the strips' clearance would narrow toward
   0.6 km or the edge strips be oversized (33 million tiles and 53 Gt in the comparison's count), departing from the
   design's 1 km and 28 million tiles.
4. **Treat the tiles' areal mass as a design variable bounded by the photon authority**, the edge rings' steering
   first and then the keeping, with the keeping flown at 20,000 km on the rolled layout. It would replace the one
   allocation of 13.765 g/m² with a bound that the hardware work designs to.
5. **Carry the climate window's new transmission when the 0.9685 run lands**, through the redirected fraction in the
   dynamics' scenario and the glow's design sunlight. It would change the common photon acceleration (about 12% less)
   and the glow's dimming.
6. **Screen the shingle angle and tile size together** with the supports, the clearance gate and the attitude study.
   It would depart from the flown 1.5° shingle and its clearance envelope, shadow torques and store.
7. **Record the 14 October direction for the coupled model in the register** if the author confirms it (finding 13).

## Cleanup made in this review

- **The Codex notes, rewritten in the author's register**, every number and limit kept: tile_design.md,
  support_comparison.md, flexible_coupling.md, the README's three new paragraphs, array_industry.md's new section, and
  the new sections of engineering/README.md and protection/dynamics/README.md. Denials and repeated disclaimers became
  one positive statement of scope near the outputs they bound; American spellings became British (aluminium, fibre,
  neighbour, modelled, disc); long lines were wrapped at 120 characters.
- **Added from the review:** flexible_coupling.md's table of the lightest passing support by shingle angle, the
  vanishing of G for this product's modes, and the initial projection error in the modal truncation; tile_design.md's
  share of the roll in the membrane load and the side members' overlap behind the clearance result; the README's link
  to this review and the conditions of the photon keeping margin; two British spellings in older README lines.
- **status.json:** the new engineering and protection entries rewritten with their numbers, this review recorded, and
  the author's 14 October direction marked for the register. **plan.md:** item 3 points to this review.
- **Source registers:** the corrections of finding 10, text only.
- **Code:** flexible_benchmark.py writes the reading rules of finding 9, and test_flexible_benchmark.py checks them.
- **Products regenerated** from their producers, single-threaded, after the register and producer edits: tile_design
  (6.1 CPU s), support_comparison (2.3 s), geometry_screen (310 s) and flexible_benchmark with flexible_models (125 s),
  with a climate model on four other cores. Every value is unchanged apart from the input and producer hashes, the run's
  resources and the new reading rules; the benchmark's timings, 15–25% longer on the loaded machine, are in the
  regenerated product and flexible_coupling.md. The check records (tile_checks.json, support_checks.json,
  flexible_checks.json) describe the checks of the original commits and are kept as they were.

## Checks

`OPENBLAS_NUM_THREADS=1 nice -n 10 make check` ran on the final code and products before the commit, in 8 min 36 s
with a climate model on four other cores: the layer check passed 988 files with no violations, **2,253 Python tests
passed and 16 skipped** (475 s), all **191 JavaScript tests passed** (129 immersion, 48 column and transport, 14
calendar), and the provenance check and the ensemble validator passed. The focused tests of the regenerated products
and the changed producer, 104 of them, passed first. The check's log is `make_check.log` in the review's folder on the
research drive. These are numerical checks: the estimates above are arithmetic on products, and flying them is the
work of the recommendations.
