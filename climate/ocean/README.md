# Ocean

Two ocean models for the Open Moon's seas, set up from pinned sources with the Moon's gravity, radius and rotation, and checked against geography's linear tide solver. The first study run on them is Mare Fecunditatis's monthly tide through its straits. Climate-forced runs of the seas' circulation wait on the base climate work (the gap with CM1); the models are ready for them.

| File | Holds |
|---|---|
| [sources.json](sources.json) | The upstream sources, each pinned by tag and commit: MITgcm checkpoint69s and SCHISM v5.13.0, with the NetCDF packages SCHISM links against and their hashes |
| [builds.py](builds.py) | Clones, checks the pins, patches and builds both models on the research drive |
| [schism/lunar.patch](schism/lunar.patch) | SCHISM's lunar build (below) |
| [lunar.py](lunar.py) | The models' lunar settings from `shared/constants.json`, the seas' fresh water and the bed drag |
| [terrain.py](terrain.py) | The Near-side Sea's floor at 16 px/deg above the GRAIL geoid, and its connected water |
| [forcing.py](forcing.py) | The Moon's tide as forcing: geography's degree-2 and degree-3 monomials, their self-attraction, and the files each model reads |
| [mitgcm.py](mitgcm.py) | MITgcm experiments: grids with a fine region, namelists, the tide's geopotential, chunked runs with pickups, output readers |
| [mesh.py](mesh.py) | Triangle meshes for SCHISM from the terrain's pixels (a quadtree, no hanging nodes) |
| [schism.py](schism.py) | SCHISM experiments: mesh files, namelists, the lunar tide files, chunked runs with hotstarts |
| [checks.py](checks.py) | Both models against geography/tide_dynamics.py and an analytic channel ([results](../results/ocean/)) |
| [fecunditatis.py](fecunditatis.py) | The strait study (below) |
| [strait_numerics.py](strait_numerics.py) | The strait study's first day under each model's step and SCHISM's filters ([results](../results/ocean/strait_numerics.json)) |

Builds, runs and products live on the research drive (`TERLUNA_OCEAN_DRIVE`, default `/media/projectspace/terluna-research/ocean`); `products/` and `runs/` here link to it and stay out of Git. The Python tools run in the `ocean` virtual environment beside the project's others (`venvs/ocean`: the system's numpy and scipy, with `netCDF4` and CMake 3 added; SCHISM's CMakeLists sets a policy CMake 4 refuses). The drive's `run.sh` runs a module with the numerical libraries held to two threads, niced and under a memory cap:

```sh
/media/projectspace/terluna-research/ocean/run.sh -m climate.ocean.builds schism      # about two minutes
/media/projectspace/terluna-research/ocean/run.sh -m climate.ocean.checks mitgcm      # channel, seiches, forced line
/media/projectspace/terluna-research/ocean/run.sh -m climate.ocean.checks schism
/media/projectspace/terluna-research/ocean/run.sh -m climate.ocean.fecunditatis grids
/media/projectspace/terluna-research/ocean/run.sh -m climate.ocean.fecunditatis run mitgcm
/media/projectspace/terluna-research/ocean/run.sh -m climate.ocean.fecunditatis run schism
/media/projectspace/terluna-research/ocean/run.sh -m climate.ocean.fecunditatis analyse
```

## The models on the Moon

**MITgcm** builds from unchanged source, with each experiment's compile-time options (its grid size, packages, the forcing package's tide reader and implicit bed drag). Its gravity (`gravity`, `gBaro`), reference sphere (`rSphere`, 1,737.4 km) and rotation (`rotationPeriod`, the sidereal month of 27.32 days, with `selectCoriMap = 2` for 2 Omega sin(latitude)) are namelist values, and its forcing package reads a gridded tidal geopotential, which the momentum equations take as -grad(phi): the file holds the equilibrium height F and the namelist scales it by -g. Binary inputs are big-endian, as the gfortran build reads them. The seas run as one depth-averaged layer. The checks hold the free surface linear, as geography's solver does, on the solver's own cells down to its 1-m floor. The strait study carries each column's full thickness (the z-star free surface, `nonlinFreeSurf = 4` with `select_rStar = 2`), as SCHISM always does: the equilibrium monthly range round Fecunditatis is 4.8-6.0 m, as deep as the strait grid's shallowest cells (5 m; shallower cells count as land). MITgcm cannot dry a cell and stops if a column thins below 5% of its depth (`hFacInf`); SCHISM dries and rewets them.

**SCHISM** holds gravity and rotation as compiled constants and builds Earth's tidal potential from Earth's species patterns, with Earth's Love-number factor 0.69 written into the code. [lunar.patch](schism/lunar.patch) changes four files:

- `grav` = MOON_GM / MOON_RADIUS² = 1.6242 m/s² and `omega_e` = 2π over the sidereal month = 2.6617 × 10⁻⁶ rad/s;
- the reference sphere's default radius is the Moon's;
- `ilunartide = 1` reads the Moon's tide as a basis of equilibrium-height fields at the nodes (`lunar_tide_basis.in`) with their coefficients through time (`lunar_tide.th`), interpolated each step;
- Earth's constituents (`ntip > 0`) and the `sflux` atmospheric forcing (`nws = 2`, whose bulk formulas still carry Earth's 9.81) stop the run with a message;
- the output files' metadata names the Moon's sphere.

The build's ledger (`terluna_ledger.json` beside the binary) records the upstream commit, the patch's hash and each patched file's hashes before and after. SCHISM runs in its two-dimensional mode (two vertical levels), with elevation and depth-averaged velocity written by one writer rank beside the compute ranks. Two quirks of v5.13.0's sample `param.nml` are handled in [schism.py](schism.py): it gives the vegetation settings `veg_lai` and `veg_cw` real values that the code declares as integers, which gfortran cannot read (vegetation is off, so they keep the code's defaults), and the code switches on three-dimensional layer heights by default, which a 2-D run turns off.

Bed friction is quadratic, tau = rho Cd |u| u, with Cd = 0.0025 and a bracket of 0.001-0.005. A drag coefficient from a logarithmic boundary layer holds no gravity, so it carries to the Moon unchanged; Manning's n, which hides g, is not used. MITgcm's strait run adds Smagorinsky viscosity (`viscC2smag = 2.2`), capped at the explicit scheme's stability limit (`viscAhGridMax = 0.5`, L²/8Δt). SCHISM's strait run takes its Smagorinsky-like Shapiro filter at the same strength (below). Both models solve the drag implicitly. SCHISM does by design; MITgcm does with `selectImplicitDrag = 2`, which solves it together with the surface pressure and needs `ALLOW_SOLVE4_PS_AND_DRAG` at compile time. MITgcm's default steps the drag explicitly within its Adams-Bashforth scheme, which holds only while r dt / H stays below about 1: the forced check at a hundred times the base friction reaches 4.5 on the solver's 1-m cells, and with explicit drag it ran away to NaN.

A run proceeds in chunks and resumes from its last pickup or hotstart, and only under the setup it started with: each run folder stamps its namelists or settings, fields, mesh, tide files and build (`setup.sha256`), and a changed setup must start again from a cleared folder. Each check starts from a cleared folder.

## Checks

[checks.py](checks.py) runs each model on cases with known answers and writes [checks_mitgcm.json](../results/ocean/checks_mitgcm.json) and [checks_schism.json](../results/ocean/checks_schism.json). A period comes from the least-squares fit of a mean and one sinusoid with the smallest residual; that recovers a test signal's frequency to one part in 10⁸, where the fitted amplitude's peak would sit about 0.1% off over a dozen cycles.

| Check | Known answer | MITgcm | SCHISM |
|---|---|---|---|
| A channel 10° long and 2° wide on the equator, 1,000 m deep: its slowest seiche | 2L/√(gH) = 4.1801 h; solver 4.1809 h | 4.1810 h (+0.022%) | 4.1799 h (−0.004%) |
| Smythii–Marginis's three slowest seiches on the atlas's quarter-degree cells | solver: 40.18, 20.65, 16.73 h | 40.18, 20.65, 16.73 h (within 0.0013%) | 33.47, 17.72, 13.90 h |
| The same coast with every cell 1,000 m deep | solver: 24.34, 9.61, 8.20 h | | 23.74, 9.30, 7.77 h |
| Fecunditatis under the 9.61-day line 2D + M′, with ten and a hundred times the base friction, at an 1,800 s step with the forcing every 6 hours | solver: 5.434 and 0.614 times the equilibrium height | 5.398 (−0.67%) and 0.613 (−0.21%), 0.33–0.34° behind; over the Near-side Sea's cells 0.8% and 0.3% rms apart | |
| The hundredfold case at a 900 s step with the forcing every 3 hours | solver: 0.614 | 0.614 (−0.06%), 0.17° behind; 0.13% rms apart | |
| Smythii–Marginis under the monthly line M′ through the patch's tide reader, without friction | solver, base friction | | gain 0.9998 on the solver, phase +0.04° |

The channel fixes each model's lunar gravity and radius. MITgcm shares the solver's discretization, an Arakawa C grid whose faces take the shallower cell's depth, with the same four-point Coriolis average, and its seiches agree with the solver's to about one part in 10⁵. Under the forcing, what remains comes from the time-stepping. The forcing's records, interpolated linearly in time, carry a line at 1 − (ωT)²/12 of its amplitude, T the records' interval: the 9.61-day line falls 0.22% short with records every six hours and 0.06% with records every three, against the hundredfold case's deficits of 0.21% and 0.06%. The response lags the solver's by a phase that halves with the step, 0.34° at 1,800 s and 0.17° at 900 s. At ten times the base friction, near resonance, the gain falls 0.67% short; that case was not repeated at the shorter step. The strait study steps 300 s with the forcing every 3 hours, where the interpolation takes 0.007% from the monthly line and 0.06% from the 9.61-day line, and the lag is under 0.06°.

SCHISM puts depths at the corners and interpolates them across each triangle, and its waves run faster. On Smythii's own floor the difference is 14–17%, which the solver's face-depth rule accounts for: with the mean or the deeper cell's depth at each face, its slowest seiche takes 34.9 or 32.5 hours, bracketing SCHISM's 33.5. On a uniform floor 2.5–5% remains, from the two methods' treatments of a staircase coast. These are differences of numerical representation, which the strait study measures directly by running both models on one floor. Through the patch's tide reader, SCHISM's monthly response matches the solver's to 0.02% in gain.

## Mare Fecunditatis's tide through its straits

Geography's linear solver found Mare Fecunditatis joined to the Near-side Sea through narrow, shallow passages that set a 9.5-day exchange, close to the 9.61-day tidal line, so that bed friction decides its tide ([geography/README.md](../../geography/README.md), The monthly tide). [fecunditatis.py](fecunditatis.py) runs the whole Near-side Sea in both models under the Moon's tide from the ephemeris, from 2026-01-01 and ramped over five days, with the passages at 16 pixels per degree and quadratic bed friction (Cd 0.0025):

- MITgcm on a latitude-longitude grid of 1/16° over Fecunditatis, its western channels and its passage to Crisium (33-66° E, 15° S-13° N), widening to 1/2° elsewhere: 814 × 634 cells, 183,698 of them water, 5.23 million km²;
- SCHISM on 263,723 triangles (138,353 nodes) over the same terrain, 1/16° in the same box and up to 1/2° beyond, 5.22 million km².

Both step 300 s. The first 30 days spin the seas up and the next seven anomalistic months (192.9 days) are analysed. The Fecunditatis basin, 244,000 km², is the water joined to its deep core below 160 m, with shallower water given to the nearest core, so its boundary runs through the sills and its volume's change is the exchange through the straits. Each chunk's snapshots (every three hours) are reduced as it finishes: the basin's mean level against its equilibrium in the monthly, fortnightly and nine-day bands, each month's range at every water point, the exchange, the largest currents in the western channels and the Crisium passage, and the bed's dissipation from the mean cube of speed. [visualization/ocean-tides/render.py](../../visualization/ocean-tides/render.py) draws the basin's level, the monthly ranges and the largest currents from the reduced runs.

Before the long runs, [strait_numerics.py](strait_numerics.py) ran the study's first day, as the tide spins up from rest, under each model's step and SCHISM's filters ([strait_numerics.json](../results/ocean/strait_numerics.json)). Each case is read on the largest speeds in the western channels and the Crisium passage at hour 21, the last snapshot both models share, and against its reference over the whole day on those speeds and each water point's largest speed:

| Case | Western channels, Crisium passage at hour 21 | Against its reference over the day |
|---|---|---|
| MITgcm, 300 s | 0.60, 0.72 m/s | |
| MITgcm, 150 s | 0.60, 0.72 m/s | within 0.1% |
| SCHISM's sample filter (`ishapiro = 1`, 0.5), 300 s | 0.43, 0.36 m/s | |
| the same at 150 s and at 600 s | 0.35, 0.21 and 0.49, 0.52 m/s | straits 13-46% apart; points 26% rms |
| SCHISM's Smagorinsky-like filter at MITgcm's strength, 300 s | 0.56, 0.69 m/s | |
| the same at 150 s and at 600 s | 0.56, 0.68 and 0.56, 0.70 m/s | straits within 3.3%; points 7% rms |
| the same at 300 s with backtracking sub-steps chosen by the flow | 0.56, 0.69 m/s | within 0.5% throughout, at half the cost |

MITgcm does not depend on its step. SCHISM's sample filter, the setting its documentation suggests for coasts, estuaries and rivers, smooths the velocities once a step at a fixed strength, so its smoothing per unit time falls as the step grows and the straits' currents moved by up to 46% with the step. Its Smagorinsky-like filter (`ishapiro = 2`) sets each step's strength to 0.5 tanh(γ₀ Δt μ), μ the local deformation rate, which is proportional to the step while the tanh is unsaturated. A filter of strength s acting once a step on cells of size L diffuses momentum like a viscosity s L²/(4Δt), so γ₀ = 8√2 (2.2/π)² = 5.55 matches MITgcm's Smagorinsky viscosity; on triangles the conversion holds within a factor of order one. At that strength the straits' speeds hold within 3.3% from a 150 s to a 600 s step and each point's largest speed within 7% rms, against 13-46% and 26% under the sample filter, and the straits' currents at hour 21 come within 4-6% of MITgcm's. SCHISM's sample also forces at least ten backtracking sub-steps a step (`dtb_max = 30 s`); letting the flow choose them (`dtb_max` = the step) changed nothing that matters and halved its cost. The strait run takes both: a 300 s step, the Smagorinsky-like filter at γ₀ = 5.55 and the flow's sub-steps.

### What the runs find

Both models ran the 222.9 days on four processes, MITgcm in 3.6 hours and SCHISM in about 11, and they agree on Fecunditatis's tide within 5% ([results/ocean/fecunditatis.json](../results/ocean/fecunditatis.json)).

- **The straits halve Fecunditatis's monthly tide and delay it by about two and a half days.** The basin's level swings 2.84–2.99 m in the median month (MITgcm, SCHISM) against 5.21 m in equilibrium. In the monthly band it carries 0.58–0.60 of the equilibrium's amplitude, 31–34° (2.4–2.6 days) behind, and in the fortnightly band 0.40–0.47, 57–69° (2.2–2.6 days) behind. The 9.61-day line, which geography's linear solver found resonating at 39 times its equilibrium height under weak friction, keeps 0.73–0.84 of it, within 16° of its phase, close to the solver's hundredfold-friction case (0.6): the passages' quadratic friction damps the resonance away.
- **The western channels carry tidal races.** The basin's exchange peaks at 1.4–1.5 million m³/s against 2.8 million in equilibrium. The currents in the western channels reach 2.2–2.3 m/s, with a daily peak of 1.5–1.6 m/s on the median day and 2.1–2.2 m/s on the strongest twentieth of days; the Crisium passage reaches 1.8 m/s, 1.2–1.3 m/s on the median day. The bed in the western channels dissipates 0.5–0.7 GW of the Near-side Sea's 6.5–7.4 GW.
- **Crisium, behind its passage, is choked further.** Its range is 0.48–0.63 of the equilibrium's over the same months, the two models' widest difference. At its centre the monthly tide arrives 85° (6.5 days) late and the fortnightly at a quarter of its equilibrium height (SCHISM's station).
- **The open sea follows its equilibrium tide.** Over the same months Procellarum, Imbrium, Frigoris, Serenitatis, Tranquillitatis and Nectaris hold 0.98–1.06 of the equilibrium's median monthly range in both models, and their stations' monthly bands lie within 6% and 5° of it. The maps' low ranges over Imbrium and high ones in the south-west belong to the equilibrium tide of these months. Its two monthly lines, 27.55 and 27.21 days, beat over about six years and redraw the pattern, so the analysed months' equilibrium ranges differ from geography's 2026–2045 medians by −42% to +31% across the sea; each mare is compared with the equilibrium over the same months.
- **The south-western maria swing 8–14% beyond their equilibrium.** Nubium, Humorum and Cognitum, joined to Imbrium by a neck across western Procellarum 120 m deep, take 1.08–1.14 of the equilibrium's range in both models, and their stations lag it by 11–15° in the monthly band and 29–43° in the fortnightly.
- **Seiches and other motions under five days** hold 0.5–4 cm rms at the stations.

| | MITgcm | SCHISM | Equilibrium, same months |
|---|---|---|---|
| Fecunditatis basin, median monthly range | 2.84 m | 2.99 m | 5.21 m |
| Monthly band: gain, lag | 0.58, 34° | 0.60, 31° | |
| Fortnightly band | 0.40, 69° | 0.47, 57° | |
| Nine-day band | 0.84, 1° | 0.73, −16° | |
| Exchange through the straits, largest | 1.41 million m³/s | 1.52 million m³/s | 2.77 million m³/s |
| Western channels: largest speed; median daily peak | 2.23; 1.48 m/s | 2.34; 1.57 m/s | |
| Crisium passage | 1.81; 1.31 m/s | 1.82; 1.22 m/s | |
| Bed dissipation: Near-side Sea; western channels | 7.4; 0.50 GW | 6.5; 0.73 GW | |
| Crisium, median monthly range | 2.11 m | 2.80 m | 4.42 m |

The models differ most where the passages' staircase coasts and the cube of the speed weigh most: by 46% in the western channels' dissipation and 30% in Crisium's range. The stations come from SCHISM's run alone; MITgcm's run predates them. Each run keeps the level through time at the centres of eleven maria and the two strait zones.

## Limits

- The seas are fresh water at one density; no temperature, salinity or stratification enters the tide runs.
- The floor is LOLA's at 16 px/deg and the water level the atlas's, on the dry Moon's geoid. Geography's sea-level equation ([geography/README.md](../../geography/README.md), The crust's response to the water) finds that with 28% cover held, the water's own gravity and the crust's elastic sag leave the western channels within 2 m of their 109 m, deepen the Crisium passage from 155 m to 197 m and shoal the south-western maria's neck from 120 m to 86 m. The water also turns the Moon by about 8°, which moves the tide's pattern across the seas (the Moon's orientation study, research/studies/moon_orientation on main). The atlas work redraws the floor and the level; a screen with geography's linear solver on that geography decides whether these runs repeat.
- Both models resolve the straits with 1/16-degree cells (about 1.9 km); channels narrower than a few cells are rendered as staircases, and the two models' different treatments of those staircases measure part of the uncertainty.
- The water's self-attraction is the equilibrium tide's (9% of the forcing on the largest seas), held fixed through each run.

## Next work

- The seas' circulation under the GCM's forcing (winds, heat and water through the lunar day and night), once the base climate settles: deep-water renewal past the sills against the sea-oxygen study's columns, the depth the night's convection mixes against the GCM's 50 m slab, and the heat the seas carry, returned to ExoPlaSim through its prescribed ocean heat flux (`nfluko = 1`) and its slab depth (`mldepth`).
- MITgcm's forcing package's bulk formulas audited for Earth terms before any atmosphere-forced run.
- One lunar-night convection column against the sea-oxygen study's Procellarum column.
- The South Pole-Aitken sea's Tenebrarum sub-basin, geography's other strait, on the same footing.
- A screen of the strait study with geography's linear solver on the deformed, turned geography of the atlas work, which decides whether its two runs repeat.
