# A mass-flux convection scheme for PlaSim, tested offline on the CM1 boxes (26 October 2026)

The author approved this work on 26 October on the condition that the scheme stand a good chance of bringing PlaSim's
lunar climate and weather nearer CM1's physics. This record describes the scheme in [massflux.py](massflux.py) and its
offline tests, results in [massflux_offline.json](../results/gcm/massflux_offline.json). Nothing was put into PlaSim.

**What the tests show.** The scheme is the IFS's form of Tiedtke's bulk mass-flux scheme. It mixes at the bulk rates
CM1's own clouds mix at, holds the sea box's diluted upper air in a running column, and carries the GCM's warmer column
most of the way to it within three lunar days. It also gives the sea's capped night cloud, the sea's afternoon and
evening precipitation and mass flux, and the land's quiet nights. It misses four things:
- **The first hours of the sea's night.** It convects deep then, forming four times CM1's night precipitation.
- **Most of the land's afternoon on CM1's own columns.**
- **Surface-rooted convection on PlaSim's grids.** Their shared lowest layer, 0–1.73 km, averages away the
  surface-layer air CM1's updraughts carry.
- **CM1's deep cloud tops.** The bulk plume stops 15–25 km below them.

## Why a mass-flux scheme

The simplified Betts–Miller scheme now in PlaSim tells deep from shallow convection by moisture alone, and it could not
carry mixing: a diluted reference moved its heat lower and thinned its cloud ([README](README.md), the scheme
sections). CM1 showed that mixing does most of the work at the sea box. Its night cloud stops at 7–8 km because its
updraughts mix 1.2–1.6×10⁻³ of their mass a metre under a stable layer, while its afternoon cloud mixes 0.2–0.3×10⁻³ and
reaches 22–28 km. A zero-buoyancy plume mixing 0.5–1×10⁻⁴ m⁻¹ reproduces the box's upper air, which runs 7.7–10 K below
the GCM's at 42–56 km ([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)).
A plume scheme carries mixing directly.

## The scheme

**Structure** (IFS Cy36r1 and Cy40r1, Part IV, chapter 6; [massflux.py](massflux.py)):
- **The updraught.** One bulk updraught carries liquid-ice static energy and total water, the variables its
  condensation and freezing keep. Its mixing is solved exactly over each half layer with the surrounding air linear in
  height, so it holds on PlaSim's 5–30 km layers. Its condensate freezes between 273 and 250 K by the IFS's
  quadratic share, and rains at the IFS's rates per second (c0 1.4×10⁻³ s⁻¹, l_crit 0.5 g/kg, nothing below
  0.3 g/kg, with the Bergeron–Findeisen enhancement).
- **Its top.** The updraught's kinetic energy (IFS (6.10)) sets the top where it runs out. Organized detrainment
  shrinks the mass flux as the square root of that energy where the updraught is negatively buoyant.
- **The type.** Test parcels decide it in the IFS's order. A cloud thicker than 200 hPa from a departure level
  within 350 hPa of the ground convects deep. A thinner cloud from the ground convects shallow. A moist level
  (relative humidity at least 0.8) under large-scale ascent gives mid-level convection.
- **The downdraught.** It is saturated and starts at the level of free sinking with a third of the cloud-base flux.
  It evaporates rain to stay saturated, at most what falls through it, and spreads over the lowest 60 hPa.
- **The tendencies.** They are the convergence of the draughts' fluxes against the surrounding air's compensating
  motion, taken upwind. The fluxes fall linearly in pressure below the cloud base. The rain's formation, evaporation
  and melting are the only sources, and detrained condensate evaporates where it leaves the cloud, since PlaSim
  carries no cloud water.
- **Conservation.** The column gains exactly the latent heat of the precipitation that leaves it and loses exactly that
  water, to round-off ([test_massflux.py](../tests/test_massflux.py)).
- **PlaSim's evaporation of falling rain.** For PlaSim the scheme hands its rain and snow by layer, net of the
  downdraught's evaporation, to `mkrain`. A copy of `mkrain` here may evaporate only what the downdraught does not
  need lower down, which keeps the water.

**The closures:**
- **Deep** (Bechtold and others 2014). Convection removes its entraining, density-weighted CAPE (PCAPE), less the part
  the sub-cloud layer's own heating supplies (PCAPEbl), over τ. τ is 5 times the cloud's turnover time H/w̄, with no
  Earth bound on it.
- **Shallow.** It carries out of the sub-cloud layer the moist static energy its other processes bring in, over the
  cloud base's excess.
- **Mid-level.** It carries the large-scale mass flux at its base.
- **The limit.** The mass flux may cross no layer faster than that layer's mass a step.

**Lunar parameters.** Lengths per metre come from CM1 where it measures them. Depths and departure ranges stay in
pressure, which a column stretched by the ratio of the gravities keeps. Earth's rates per metre are divided by that
ratio, 6.04, where CM1 measures none.

| Parameter | Value | Source; evidence |
|---|---|---|
| Entrainment, deep and mid-level | ε = 0.8×10⁻⁴ + 2.0×10⁻⁴ exp(−(z − z_base)/2 km) m⁻¹ | CM1's bulk rates (box_entrainment's sampling, both boxes): 0.06–0.1×10⁻³ m⁻¹ through the middle and upper troposphere, 0.2–0.6×10⁻³ within about 2 km of the base wherever it lies (4 km on land, 10 km over the afternoon sea); the zero-buoyancy plume's 0.5–1×10⁻⁴ aloft; measured |
| The same per e-folding of pressure | ε R T_v/g: 4.1 aloft (3.6 at 40 km), 14.5 at the base | the scale height R T_v/g, 51.8 km at 292 K |
| Shallow entrainment | twice the deep | the IFS's factor; CM1's night cores 0.18–0.57×10⁻³ m⁻¹; consistent |
| Detrainment | δ = ε (1.6 − RH), plus organized | the IFS's shallow form, taken for every type because CM1's deep updraughts detrain nearly as fast as they entrain above 12 km |
| Kinetic energy | the IFS's virtual mass 0.5, factor 2, drag 1.875 × 0.506; 1 m/s at a deep test's departure | Earth values, dimensionless or a speed, which the stretch keeps |
| Test parcels | the lowest level's air plus 0.2 K and 0.1 g/kg, or 1.5 times each surface flux over ρ w* if larger; sub-cloud mixing at 0.4 times the deep base rate | the IFS (6.19)–(6.21), with w* the mixed layer's convective velocity so the formula scales in a 4–6 km boundary layer; elevated departures only under large-scale ascent, the IFS's mid-level condition, because CM1's land evenings convect weakly from a moist layer under sinking |
| Downdraught | η 0.3; mixing 0.33×10⁻⁴ m⁻¹ | the IFS; the mixing stretched |
| τ | 5 H/w̄ | CM1's convection removes the scheme's PCAPE over 1.6–6.9 turnover times in its four quarters of deep convection, median 5; the IFS's 1 + 264/n at the lunar grid's Earth-equivalent truncation, 77, gives 4.4; measured, partly on the sea quarters it is tested on |
| PCAPEbl's T* | g H/c_p, H the cloud's depth | Bechtold and others' definition, which they set to 1 K "and cast the scaling into α"; an interpretation (below) |
| τ_bl | the turnover time over land; cloud-base height over the sub-cloud wind, at least 2 m/s, over water | the IFS |

The (1.3 − RH) factor on entrainment (Bechtold and others 2008) is an option left off: CM1's night and afternoon rates
differ by up to four times, against the factor's 1.2.

## What CM1 measures for it

**The source air.** In convecting quarters CM1's cloudy updraughts reach their base with the moist static energy of the
lowest level's air, 50 m up: 346.8 against 346.5 kJ/kg on land in the afternoon. At night over the sea they carry the
well-mixed 6 km boundary layer's air with 0.5 kJ/kg to spare, and keep a quarter of the surface air's excess. So the
air at the base depends on the regime: cold pools leave a warm sector in the afternoon, and the night's boundary layer
is mixed. No single rate of sub-cloud mixing gives both.

**The adjustment time.** On each output where the scheme convects deep, CM1's own convective tendencies (microphysics,
resolved transport, boundary layer) would remove the scheme's PCAPE, its updraught held, at a rate set against the PCAPE
itself:

| Quarter (deep outputs) | PCAPE, Pa | Turnover, h | CM1's time, h | Ratio |
|---|---|---|---|---|
| Sea afternoon (36) | 95 | 4.4 | 24.5 | 5.7 |
| Sea evening (21) | 144 | 3.6 | 26.3 | 6.9 |
| Land morning (10) | 92 | 3.3 | 11.8 | 3.6 |
| Land afternoon (17) | 63 | 5.4 | 9.5 | 1.6 |
| Sea night (19) | 78 | 5.3 | 926 | 175 |

CM1 removes PCAPE three times faster, per turnover, over land in the day than over the sea. The Bechtold term slows
convection while the boundary layer heats, the other way. Taken with T* = 1 K it would exceed PCAPE 18–37 times on land
and stop the lunar day's strongest convection, so the scheme takes T* as g H/c_p, 40–60 K for these clouds. CM1 hardly
removes the scheme's PCAPE in the sea's early night.

## The tests

**1. CM1's own columns** (every 3-hourly output, 59 a quarter, on CM1's levels; precipitation formed in mm/day against
the falling water CM1 forms, its evaporation on the way down plus its rain at the ground):

| Quarter | Types (deep, shallow) | Cloud: base; top median and 90th percentile | Base mass flux | Formed: scheme; CM1 |
|---|---|---|---|---|
| Sea night | 32%, 68% | shallow 5.9; 8.2, 8.6 km (CM1's tops 7.3, 7.9) | shallow nil (supply below nil); CM1 0.014 at 6 km | 2.21; 0.53 |
| Sea morning | 0, 27% | 5.9; 6.8, 7.1 km (CM1 6.7, 7.3) | nil | 0.00; 0.01 |
| Sea afternoon | 61% elevated, 2% | deep 10.4; 46, 54 km | 0.012; CM1 0.014 at 12 km | 5.22; 6.47 |
| Sea evening | 36% from the ground, 63% | deep 3.9; 35, 41; shallow 5.5; 6.6, 7.7 km | deep 0.026; CM1 0.021 at 6 km | 5.97; 6.62 |
| Land night | none | none | nil | 0.00; 0.33 |
| Land morning | 17%, 5% | deep 6.7; 34, 36 km | 0.012 | 2.23; 7.56 |
| Land afternoon | 29%, 8% | deep 6.4; 40, 43 km | 0.003 | 1.00; 22.63 |
| Land evening | none (mid-level 10%) | none | nil | 0.00; 3.42 |

- **The sea.** The night's cloud is capped where CM1's is, half a kilometre high, and its shallow closure finds the
  sub-cloud layer losing energy without convection. So it gives no flux, where CM1's boundary layer vents the evening's
  stored energy at 0.014 kg m⁻² s⁻¹. The early night's deep outputs come while the boundary layer still holds the
  evening's moisture, its parcel clearing the 6–12 km barrier by 0.3 kJ/kg. The afternoon convects from 10 km, where
  CM1's mass flux grows tenfold, and the evening from the ground at 3.9 km.
- **The mass flux.** Its profiles follow CM1's sampled cloudy updraughts to 30 km: 0.011–0.014 kg m⁻² s⁻¹ against
  0.014–0.022 in the afternoon, 0.010–0.020 against 0.010–0.021 in the evening.
- **The tops.** The plume ends at 35–54 km, where CM1's mass flux runs on to 60–66 km. CM1's cores, buoyant by
  0.4–0.8 K, carry it there, while the bulk plume, mixed at the bulk rates, is near neutral aloft and its kinetic energy
  drains. So the afternoon's heating at 16–35, 35–48 and 48–63 km is 52, 8 and −2 W/m² against CM1's 27, 26 and 44.
- **The land.** The afternoon convects deep on 29% of outputs because the domain-mean boundary layer is stable under
  cold pools, and PCAPEbl takes two-thirds of the PCAPE.
- **The variants.** Tiedtke's undiluted test parcel (sub-cloud mixing nil) forms 8.0 and 7.2 mm/day on land in the
  morning and afternoon, but 5.9 over the sea at night and 9.3–10.0 in the sea's afternoon and evening. Dropping the
  boundary-layer term moves the land's afternoon only to 1.8.

**2. The GCM's site columns** (the rh 0.7 branch, years 10–14, its six sea and four land cells). The scheme convects on
none of the sea's columns at any hour and on 7% of the land's afternoon columns, forming 1.4 mm/day where the branch's
convection rains 13.3. The branch's upper air, built by an undiluted reference, is 5.5–11 K warmer than CM1's, too
warm for a plume mixing at CM1's rates. So the first response in PlaSim would be a pause while radiation cools the air
aloft, and test 4 shows where the column then settles.

**3. Ten levels against fourteen** (each CM1 output laid on PlaSim's ten levels and on the levels agent's `l14`, from
[levels.py](levels.py), by layer means):
- **The lowest layer.** Both grids share it, 0–1.73 km. Its mean air holds 1.4–2.3 kJ/kg less moist static energy
  than CM1's 50 m air over the sea, and 1.2–2.3 kJ/kg less over land by day (medians by quarter), the energy that lets
  surface-rooted cloud form. On both grids the scheme then convects from the ground on at most 2% of the boxes'
  columns: the sea's evening and the land's afternoon go quiet.
- **The night cap.** On ten levels the night cloud (base 5.9 km, tops 7–8 km) and its 6–12 km cap lie inside one
  layer, between half levels at 6.4 and 13.7 km, so no level resolves either. On fourteen, half levels at 5.0, 7.0
  and 9.0 km cut through both.
- **Elevated convection.** Fourteen levels keep the sea afternoon's elevated deep convection on 49% of outputs,
  forming 4.9 mm/day against CM1's 6.5, where ten keep 7% and 1.1.

**4. A single column** driven by each box's own radiation, day-night table, mean ascent and surface fluxes, with the
scheme, `mkrain` at the box's constant, a saturation adjustment and a mixed boundary layer. Each box's column runs one
lunar day from the box's own air, and the GCM's site mean (5–11 K warmer aloft) runs three. Temperature less CM1's
lunar-day mean, K:

| Sea | 4 km | 8 km | 13.5 km | 19 km | 29 km | 42 km | 56 km |
|---|---|---|---|---|---|---|---|
| CM1's levels, from CM1, day 1 | +1.1 | 0.0 | −0.1 | +0.3 | +0.1 | −0.7 | −1.9 |
| 14 levels, from CM1, day 1 | −0.1 | 0.0 | −0.2 | −0.1 | −0.2 | −1.0 | −2.4 |
| 10 levels, from CM1, day 1 | −0.1 | −0.1 | −0.1 | −0.1 | −0.1 | −0.8 | −2.4 |
| 14 levels, from the GCM: start; day 3 | +3.4; +5.6 | +3.6; +3.8 | +4.1; +3.8 | +5.8; +3.9 | +8.1; +3.5 | +9.7; +2.8 | +11.3; +1.0 |
| 10 levels, from the GCM: start; day 3 | +3.4; +1.3 | +3.7; +1.9 | +3.7; +1.9 | +5.5; +3.1 | +8.1; +4.6 | +9.7; +4.2 | +11.3; +3.5 |
| CM1's levels, from the GCM: start; day 3 | +3.2; +7.8 | +3.8; +6.0 | +4.0; +3.2 | +5.6; +3.9 | +7.9; +3.3 | +9.5; +3.0 | +11.2; +1.6 |

- **Aloft.** At 42–56 km the GCM's excess falls by 5.5–10 K in three lunar days, 57–91% of it, and in the column
  started from CM1's air the 42–56 km air ends 0.7–2.4 K cool. That follows the low tops.
- **Near the ground.** Below 8 km the column from the GCM warms to 4–8 K above CM1's on CM1's levels and on fourteen,
  as its boundary layer holds the GCM's extra vapour; ten levels keep it within 2 K.
- **Rain at the ground** (mm/day by quarter, night first; CM1 0.16, 0.00, 2.27, 2.96):

| Single column | Night | Morning | Afternoon | Evening |
|---|---|---|---|---|
| Sea, CM1's levels, from CM1 | 0.48 | 0.06 | 2.12 | 1.66 |
| Sea, 14 levels, from CM1 | 0.74 | 0.00 | 0.46 | 2.58 |
| Sea, 10 levels, from CM1 | 0.95 | 0.00 | 0.85 | 2.65 |
| Land, 14 levels, from CM1 (CM1 0.07, 2.92, 8.71, 1.57) | 0.00 | 0.79 | 8.41 | 0.01 |
| Land, 10 levels, from CM1 | 0.03 | 0.51 | 9.61 | 0.24 |

- **The land on PlaSim's grids.** A running column warms and moistens its lowest layer until the scheme fires from
  the ground, so the afternoon rains at CM1's rate on both grids, though the morning and evening fall short.
- **The land on CM1's levels.** The column drifts 10 K warm below 10 km within a lunar day, and the result is set
  aside as the test column's: its mixed layer and prescribed fluxes stand in for PlaSim's vertical diffusion and do
  not hold the fine grid's land.

**5. Cost.** Counted as exponentials (each saturation, mixing and decay takes one):
- **The scheme.** It evaluates 77–99 a column on the GCM's mostly quiet columns and 148–571 on the boxes' convecting
  columns, laid on ten and fourteen levels.
- **The Betts–Miller scheme.** It evaluates 21 with one step a layer and 291 with sixteen.
- **The estimate.** The sixteen steps cost a quarter more model time a year, about 65 s on 250. Scaled by that, the
  mass-flux scheme would add 20–130 s a year, a tenth to a half of the model's time, more on fourteen levels than on
  ten. This is a count, not a Fortran timing.

## Whether to implement it, and what that would take

**The case for it.** The tests answer the approval's condition for the upper air. Mixing at CM1's measured rates holds
CM1's diluted upper troposphere and moves the GCM's toward it by 5.5–10 K at 42–56 km in three lunar days. That is the
mixing part, 5.3–6.3 K, of the gap the upper-air diagnosis traced, and the Betts–Miller scheme's mixing options could
not carry it.

**Three costs come with that change:**
- **The nights.** The sea's nights stay wetter than CM1's (0.5–1.4 mm/day at the ground in the single column, against
  0.16).
- **The deep heating.** It sits 15–25 km too low, which the column's 1–2 K cooling at 42–56 km already shows.
- **The source air.** PlaSim's lowest layer cannot hold the surface-layer air that roots CM1's afternoon and evening
  convection, so a running model fires only after its lowest layers warm, by 1–8 K in the single column.

**Before PlaSim.** One change belongs first: a surface parcel that carries the near-surface air, measured in CM1 at
1.4–2.3 kJ/kg above PlaSim's lowest-layer mean by day. A neutral surface-layer similarity estimate from the fluxes and the
friction velocity gave 6 kJ/kg over the sea, anti-correlated with CM1's, so the estimate itself needs work. The
fourteen levels help: they keep the elevated afternoon convection and resolve the night cloud's layer.

**The implementation** would be a Fortran routine in `rainmod.f90` behind `nprc` 3, with Kuo (1) and the Betts–Miller
scheme (2) kept bit for bit:
- **Inputs.** PlaSim already holds them at the convection call: the state, `dtdt` and `dqdt` from vertical diffusion
  and radiation, `dw` (dp/dt), the dynamics' moistening `dqt`, and the surface fluxes `dshfl`, `dlhfl` and `devap`.
  The dynamics' heating would need a transform like `mkdqtgp`.
- **Outputs.** The routine would hand `mkrain` its rain and snow by level under the floor rule above, and `mkclouds`
  its cloudy layers.
- **Its tests.** As for the Betts–Miller scheme: the Fortran against this Python copy to round-off, conservation, and
  the other schemes unchanged.
- **The runs.** One-year smoke runs, two-year branches, and a longwave refit, since the column would cool aloft by
  several kelvin; any fifteen-year branch is the author's call.
- **The size.** About 800–1,000 lines, roughly five times the Betts–Miller patch.

## Limits

- **The data.** CM1's columns are domain means of one realization of each box, which its own convection shaped, so the
  scheme's response at once is set against CM1's realised convection. The entrainment and τ were measured on these
  boxes, so the sea quarters are not independent tests of them.
- **The GCM's columns.** They are 3-day means of a branch built by another scheme. Their boundary-layer forcing lacks
  the dynamics' heating and takes radiation and the dynamics' moistening from one-step snapshots.
- **The single column.** It takes CM1's radiation and surface fluxes as given, so they do not answer its own state, and
  its boundary layer is a mixed layer.
- **Rain at the ground.** Applied on CM1's 111 levels, `mkrain`'s per-layer share evaporates nearly all the scheme's
  rain beneath the box's 6 km sub-cloud layer, which holds a deficit enough to evaporate about 13 mm/day at the sea's
  0.027. So rain at the ground is compared on PlaSim's grids and in the single column only.
- **PCAPEbl's T*.** It is an interpretation of Bechtold and others' scaling.
- **The downdraught.** Its parameters are Earth's.
- **The cost.** It is an estimate from counts.

## Sources

- **IFS Documentation Cy36r1, Part IV, chapter 6.** Read, from the text on the research drive
  (`session-tools/agents/climate-explain/ifs_cy36r1_part4.txt`), sections 6.1–6.11.
- **IFS Documentation Cy40r1, Part IV, chapter 6.** Read, sections 6.3–6.4 (ecmwf.int elibrary 9204,
  sha256 f9fec6f0246be4f3…).
- **Bechtold and others (2014), *J. Atmos. Sci.* 71, 734–753.** Read as ECMWF Technical Memorandum 705 (2013), the
  version submitted to the journal (elibrary 8018, sha256 2e057dd36bfae7bd…): sections 2 and 5 and the appendix.
- **Bechtold and others (2008), *QJRMS* 134, 1337–1351.** Read as ECMWF Technical Memorandum 556 (elibrary 8011,
  sha256 431b13baaf05e12e…), section 2.1.
- **Tiedtke (1989), *Mon. Wea. Rev.* 117, 1779–1800.** Not read: the AMS site refused the download, and the scheme
  follows it as the IFS documentation gives it.
- The PDFs and their text are in `session-tools/agents/climate-massflux/literature` on the research drive.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.massflux        # about 5 minutes on one core: ../results/gcm/massflux_offline.json
python3 -m pytest climate/tests/test_massflux.py              # the scheme's tests, under a second
```
