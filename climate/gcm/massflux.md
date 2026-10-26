# A mass-flux convection scheme for PlaSim: offline tests on the CM1 boxes and a Fortran routine (26 October 2026)

The author approved this work on 26 October on the condition that the scheme stand a good chance of bringing PlaSim's
lunar climate and weather nearer CM1's physics, and then approved its Fortran routine for PlaSim. This record
describes the scheme in [massflux.py](massflux.py), its offline tests (results in
[massflux_offline.json](../results/gcm/massflux_offline.json)), and its Fortran routine in
[massflux_fortran.py](massflux_fortran.py), which waits outside the runner's patches until PlaSim runs on fifteen
levels with it.

**What the tests show.** The scheme is the IFS's form of Tiedtke's bulk mass-flux scheme. It mixes at the bulk rates
CM1's own clouds mix at, holds the sea box's diluted upper air in a running column, and carries the GCM's warmer column
most of the way to it within three lunar days. It gives the sea's capped night cloud, the sea's afternoon and evening
precipitation and mass flux, and the land's quiet nights. The second round gave its surface parcel the air near the
ground on PlaSim's fifteen levels, which restores the surface-rooted convection of the sea's evening and the land's
day on PlaSim's grid, and put the scheme into Fortran, which gives the Python copy's tendencies to round-off. Three
things it still misses:
- **The first hours of the sea's night.** It convects deep then, forming four to five times CM1's night
  precipitation on CM1's columns; a running column rains 0.56 mm/day at the ground against CM1's 0.16.
- **The land's afternoon on CM1's own mean columns**, whose boundary layer is stable above its lowest 250 m. A
  running column, whose boundary layer mixes, rains at CM1's rate.
- **CM1's deep cloud tops.** The bulk plume stops 10–25 km below them.

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
- **The surface parcel.** It carries the air the IFS's surface parcel carries, in a column stretched by the ratio of
  the gravities: the air of the IFS's lowest model level, 10 m above Earth's ground and 60 m above the Moon's, with the
  IFS's excess from the surface fluxes at that height (the second round, below).
- **The downdraught.** It is saturated and starts at the level of free sinking with a third of the cloud-base flux.
  It evaporates rain to stay saturated, at most what falls through it, and spreads over the lowest 60 hPa.
- **The tendencies.** They are the convergence of the draughts' fluxes against the surrounding air's compensating
  motion, taken upwind. The fluxes fall linearly in pressure below the cloud base. The rain's formation, evaporation
  and melting are the only sources, and detrained condensate evaporates where it leaves the cloud, since PlaSim
  carries no cloud water. For PlaSim the heating goes over its moist heat capacity, c_p (1 + 0.85 q).
- **Conservation.** The column gains exactly the latent heat of the precipitation that leaves it and loses exactly that
  water, to round-off ([test_massflux.py](../tests/test_massflux.py)).
- **PlaSim's evaporation of falling rain.** For PlaSim the scheme hands its rain and snow by layer, net of the
  downdraught's evaporation, to `mkrain`, in the phase `mkrain` will give them: all of it rain in a layer above
  freezing and snow at or below, the downdraught evaporating that phase. `mkrain` then takes the heat of each change
  of phase as it does for Kuo's rain, and under the mass-flux scheme it may evaporate in a layer only what the water
  handed on below still needs, so no running total falls below nil and the water is kept.

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
ratio, 6.04, where CM1 measures none, and Earth's heights multiplied by it.

| Parameter | Value | Source; evidence |
|---|---|---|
| Entrainment, deep and mid-level | ε = 0.8×10⁻⁴ + 2.0×10⁻⁴ exp(−(z − z_base)/2 km) m⁻¹ | CM1's bulk rates (box_entrainment's sampling, both boxes): 0.06–0.1×10⁻³ m⁻¹ through the middle and upper troposphere, 0.2–0.6×10⁻³ within about 2 km of the base wherever it lies (4 km on land, 10 km over the afternoon sea); the zero-buoyancy plume's 0.5–1×10⁻⁴ aloft; measured |
| The same per e-folding of pressure | ε R T_v/g: 4.1 aloft (3.6 at 40 km), 14.5 at the base | the scale height R T_v/g, 51.8 km at 292 K |
| Shallow entrainment | twice the deep | the IFS's factor; CM1's night cores 0.18–0.57×10⁻³ m⁻¹; consistent |
| Detrainment | δ = ε (1.6 − RH), plus organized | the IFS's shallow form, taken for every type because CM1's deep updraughts detrain nearly as fast as they entrain above 12 km |
| Kinetic energy | the IFS's virtual mass 0.5, factor 2, drag 1.875 × 0.506; 1 m/s at a deep test's departure | Earth values, dimensionless or a speed, which the stretch keeps |
| Surface parcel's air | the air 60 m up, from the lowest level by Monin–Obukhov similarity in a heated surface layer: θ and q differ by (w′θ′)ₛ/u* (or (w′q′)ₛ/u*) times (ln(z_n/zₛ) − ψ(z_n/L) + ψ(zₛ/L))/κ, with φ_h = (1 − 16ζ)^(−1/3) and u* 0.1 m/s | the IFS's lowest level, 10 m, stretched; the IFS's u* of (6.20); the free-convection form of the stability function, since the lunar surface layer at 0.3 km lies near ζ = −10 (an interpretation, below) |
| Surface excess | 1.5 times each surface flux over ρ σ_w, at least 0.2 K and 0.1 g/kg; σ_w = (1.2 u*³ + 1.5 κ g zₛ/T (w′θ_v′)ₛ)^(1/3) at the source height zₛ | the IFS (6.19)–(6.21), with (6.20) evaluated where the parcel's air comes from |
| Test parcels' sub-cloud mixing | 0.4 times the deep base rate, 1.12×10⁻⁴ m⁻¹ | the IFS Cy40r1 deep test's 0.4 times CM1's rate at the base; CM1 mixes its sub-cloud layer with its boundary-layer scheme, so it measures no plume rate there |
| Elevated departures | only under large-scale ascent | the IFS's mid-level condition, because CM1's land evenings convect weakly from a moist layer under sinking |
| Downdraught | η 0.3; mixing 0.33×10⁻⁴ m⁻¹ | the IFS; the mixing stretched |
| τ | 5 H/w̄ | CM1's convection removes the scheme's PCAPE over 2.6–7.4 turnover times in its four quarters of deep convection, median 5; the IFS's 1 + 264/n at the lunar grid's Earth-equivalent truncation, 77, gives 4.4; measured, partly on the sea quarters it is tested on |
| PCAPEbl's T* | g H/c_p, H the cloud's depth | Bechtold and others' definition, which they set to 1 K "and cast the scaling into α"; an interpretation (below) |
| τ_bl | the turnover time over land; cloud-base height over the sub-cloud wind, at least 2 m/s, over water | the IFS; the sub-cloud wind is the mean by mass from the ground to the cloud base |

The (1.3 − RH) factor on entrainment (Bechtold and others 2008) is an option left off: CM1's night and afternoon rates
differ by up to four times, against the factor's 1.2.

## What CM1 measures for it

**The source air.** CM1's cloudy updraughts at their lowest bases, 4.1 km up, carry the moist static energy of the
lowest level's air, 50 m up, and a little more: 0.45 kJ/kg more on land in the afternoon and 0.60 over the sea in the
evening, while the mean air around them at that height holds 4.6 and 3.5 kJ/kg less than the 50 m air. At night over
the sea their bases lie at the top of the well-mixed 6 km boundary layer, whose air they carry with 0.44 kJ/kg to
spare (their buoyant cores 1.15), a quarter of the 50 m air's excess over that layer. So the updraughts that root
deep convection by day lift surface air through a boundary layer they hardly mix with, and the night's shallow cloud
draws on the mixed layer.

**The adjustment time.** On each output where the scheme convects deep, CM1's own convective tendencies (microphysics,
resolved transport, boundary layer) would remove the scheme's PCAPE, its updraught held, at a rate set against the PCAPE
itself:

| Quarter (deep outputs) | PCAPE, Pa | Turnover, h | CM1's time, h | Ratio |
|---|---|---|---|---|
| Sea afternoon (36) | 95 | 4.4 | 24.5 | 5.7 |
| Sea evening (23) | 150 | 3.6 | 27.0 | 7.4 |
| Land morning (13) | 103 | 3.1 | 13.5 | 4.4 |
| Land afternoon (25) | 77 | 4.5 | 11.0 | 2.6 |
| Sea night (25) | 84 | 5.1 | 759 | 150 |

CM1 removes PCAPE twice as fast, per turnover, over land in the day as over the sea. The Bechtold term slows
convection while the boundary layer heats, the other way. Taken with T* = 1 K it would exceed PCAPE 18–37 times on land
and stop the lunar day's strongest convection, so the scheme takes T* as g H/c_p, 40–60 K for these clouds. CM1 hardly
removes the scheme's PCAPE in the sea's early night.

## The second round: the surface parcel on fifteen levels

PlaSim's fifteen levels (levels.py's `l15`) split the shared lowest layer at 0.6 km, which puts the lowest full level
0.3 km up. The question was whether that layer's air carries the energy of the surface air CM1's updraughts draw on,
and, where it does not, what the parcel should take instead.

**The thin layer's air.** Its deficit is a tenth to a third smaller than the shared layer's, and most of it remains
(CM1's 50 m air less the lowest layer's mean, kJ/kg, median by quarter; `surface_layer_excess`):

| Quarter | Ten or fourteen levels (0–1.73 km) | Fifteen levels (0–0.6 km) |
|---|---|---|
| Sea night, morning | 1.38, 1.43 | 1.05, 1.09 |
| Sea afternoon, evening | 2.27, 1.61 | 1.50, 1.16 |
| Land morning, afternoon | 2.22, 2.33 | 1.98, 1.72 |
| Land evening | 1.17 | 0.77 |

Within the lowest 550 m the moist static energy falls 1.5–2.7 kJ/kg below its 50 m value, most of that in the lowest
250 m, so a layer mean at any of these depths lies well below the air the updraughts lift.

**The rule.** The IFS's surface parcel starts at its lowest model level, 10 m up on Earth, with an excess from the
surface fluxes over the velocity spread of the boundary layer's thermals at that height (Cy40r1 (6.19)–(6.20), the
profile of Holtslag and Moeng 1991 and Troen and Mahrt 1986). In a column stretched by the ratio of the gravities that
level stands 60 m up, and the scheme's surface parcel takes the air there:
- **The air at 60 m** follows from PlaSim's lowest level by surface-layer similarity with the surface fluxes, the
  profile models use for their 2 m and 10 m diagnostics. The lunar surface layer is strongly convective (the
  Obukhov length is −26 m over the sea), so the stability function takes the free-convection form φ_h = (1 − 16ζ)^(−1/3)
  (Carl and others 1973; Kader and Yaglom 1990), whose −1/3 exponent free-convection scaling requires; the IFS's own
  surface layer uses Dyer and Hicks' −1/2, fitted at ζ above about −2. In free convection the result hardly depends on
  u*: CM1's own friction velocity gives the same values as the IFS's 0.1 m/s to 0.01 kJ/kg. Where the surface's
  buoyancy flux does not rise there is no surface-rooted convection, and the parcel keeps the lowest level's air.
- **The excess** is the IFS's at the parcel's source height, where σ_w is smaller than at 0.3 km.
- **What it gives** on CM1's columns laid on fifteen levels (kJ/kg, median by quarter; the similarity's correlation
  with CM1's deficit across the outputs in brackets):

| Quarter | CM1's 50 m air less the layer | Similarity | Excess | CM1's 50 m air less the parcel |
|---|---|---|---|---|
| Sea night | 1.05 | 0.78 (−0.06) | 0.85 | −0.58 |
| Sea afternoon | 1.50 | 0.86 (0.52) | 0.88 | −0.17 |
| Sea evening | 1.16 | 0.76 (−0.73) | 0.82 | −0.38 |
| Land morning | 1.98 | 1.26 (−0.38) | 1.26 | −0.95 |
| Land afternoon | 1.72 | 0.97 (0.96) | 0.99 | −0.40 |

The similarity carries 55–75% of the 50 m air's lead over the layer, and the excess puts the parcel 0.2–1.0 kJ/kg
above CM1's 50 m air, where CM1's lowest updraught bases hold 0.45–0.60 above it. With CM1's friction velocity the
Dyer–Hicks form gives about half as much, 0.36–0.46 kJ/kg over the sea and 0.53–0.56 over land by day, and the neutral
profile two to four times CM1's deficit, 2.9–3.7 over the sea and 3.6–5.7 on land (medians by quarter).

**Tested and set aside** (CM1's columns on fifteen levels, `variants_on_15_levels`):
- **The excess at the lowest level alone**, the IFS's formula at 0.3 km without the similarity, gives the first
  round's types and rain in every quarter: the sea's evening and the land's day stay quiet.
- **The IFS's own scaling of the test parcel's sub-cloud mixing**, (q_s(T)/q_s(T_surface))³, with the excess at the
  lowest level, turned 19% of the sea's nights deep while the evening reached only 2.2 mm/day of CM1's 6.6 (a scratch
  test, `session-tools/agents/climate-massflux2` on the research drive; the option is not kept).
- **An undiluted sub-cloud ascent** turns 75% of the sea's nights deep (5.3 mm/day formed).

**On CM1's columns laid on fifteen levels** (59 outputs a quarter; precipitation formed, mm/day; `on_plasim_levels`):

| Quarter | Types (deep, shallow, mid) | Deep base; top median, km | Formed: second round; first round's parcel; CM1 |
|---|---|---|---|
| Sea night | 42%, 37%, 20% | 7.2; 45 | 2.31; 0.03; 0.53 |
| Sea afternoon | 49% elevated, 0, 20% | 11.8; 51 | 4.63; 4.63; 6.47 |
| Sea evening | 31%, 14%, 20% | 5.2; 39 | 3.97; 0.00; 6.62 |
| Land morning | 20%, 0, 34% | 7.3; 39 | 1.77; 0.09; 7.56 |
| Land afternoon | 31%, 0, 0 | 7.3; 48 | 0.98; 0.00; 22.63 |

The fifteen levels with the new parcel behave as CM1's own levels did with the first round's (sea night 32% deep,
evening 36%, land afternoon 29%): the thin layer and the similarity give the parcel the surface air.

**A running column on fifteen levels** (test 4 below; rain at the ground by quarter, mm/day):

| Column | Night | Morning | Afternoon | Evening |
|---|---|---|---|---|
| Sea, from CM1 (CM1 0.16, 0.00, 2.27, 2.96) | 0.56 | 0.00 | 1.47 | 2.50 |
| Sea, from CM1, the first round's parcel | 0.32 | 0.00 | 1.94 | 1.56 |
| Land, from CM1 (CM1 0.07, 2.92, 8.71, 1.57) | 0.11 | 0.99 | 9.23 | 0.22 |
| Land, from CM1, the first round's parcel | 0.00 | 0.73 | 7.88 | 0.06 |

The sea's evening rises from half of CM1's to 85% of it, its afternoon falls to two-thirds, and the land's morning
reaches a third of CM1's. In the column started from the GCM's air the sea's upper air at 41 and 55 km falls from
9.7 and 11.3 K above CM1's to 3.1 and 1.3 K in three lunar days.

**How the shortfalls moved:**
- **The sea's early nights** go deep on CM1's columns as they did on CM1's levels: 42% of outputs, 2.31 mm/day formed
  against CM1's 0.53. The first round's quiet nights on PlaSim's grids came from the averaged-away surface air, which
  also stopped the evening and the land.
- **The land's afternoon on CM1's mean columns** convects on 31% of outputs, forming 0.98 of CM1's 22.6 mm/day; the
  next paragraphs give the reason.
- **The tops** stay at 39–51 km (medians), where CM1's mass flux reaches 60–66 km.
- **The surface air**, 1.05–1.98 kJ/kg short in the fifteen-level layer by day, reaches the parcel to within 0.2–1.0.

**The land's afternoon trigger.** CM1's mean land boundary layer is superadiabatic only in its lowest 250 m. Above that
its density potential temperature rises, by 1.2 K from 0.25 to 6.2 km in the afternoon and 2.6 K in the morning, and
CM1's own diagnostic gives the air of its lowest level 8.7 J/kg of inhibition (median). Cold pools do not make that stratification:
they cover 7% of the afternoon's columns (16% at half a kelvin), and outside them the air above 550 m has the domain's
density potential temperature to within 0.02 K (`cold_pools`). CM1's convection rises through it by resolved ascent:
air rising at more than 0.5 m/s below its cloud carries 0.6–2.0 kJ/kg more than the mean around it, and the lowest
updraught bases carry the surface air. The scheme's surface parcel meets 1.5 J/kg of inhibition (median) on the
afternoon's columns and starts with 0.5 J/kg, the IFS's 1 m/s, of which the kinetic-energy equation counts a third
of the buoyancy (its virtual mass and turbulence factors), so it crosses at most about 1.5 J/kg and reaches its cloud
on 31% of afternoons (`trigger_on_15_levels`). PlaSim's own boundary layer stays neutral under heating: its vertical
diffusion and its dry adjustment (`mkdca`) mix away any layer whose potential temperature falls with height. The
running column, whose mixed layer stands in for both, rains 9.2 mm/day on land afternoons against CM1's 8.7. So the
trigger fails on CM1's stratified mean, which the running column's boundary layer does not hold, and the explanation
points to no change in the trigger; the two-year runs show whether PlaSim's land afternoons keep any inhibition like
CM1's.

## The tests

**1. CM1's own columns** (every 3-hourly output, 59 a quarter, on CM1's levels; precipitation formed in mm/day against
the falling water CM1 forms, its evaporation on the way down plus its rain at the ground):

| Quarter | Types (deep, shallow) | Cloud: base; top median and 90th percentile | Base mass flux | Formed: scheme; CM1 |
|---|---|---|---|---|
| Sea night | 42%, 58% | shallow 5.9; 8.3, 8.7 km (CM1's tops 7.3, 7.9) | shallow nil (supply below nil); CM1 0.014 at 6 km | 2.69; 0.53 |
| Sea morning | 0, 36% | 5.9; 6.8, 7.2 km (CM1 6.7, 7.3) | nil | 0.00; 0.01 |
| Sea afternoon | 61% elevated, 2% | deep 10.4; 46, 54 km | 0.012; CM1 0.014 at 12 km | 4.95; 6.47 |
| Sea evening | 39% from the ground, 61% | deep 3.9; 35, 41; shallow 5.5; 6.7, 7.9 km | deep 0.026; CM1 0.021 at 6 km | 6.56; 6.62 |
| Land night | none | none | nil | 0.00; 0.33 |
| Land morning | 22%, 7% | deep 6.4; 34, 36 km | 0.015 | 3.03; 7.56 |
| Land afternoon | 42%, 2% | deep 5.9; 40, 43 km | 0.003 | 2.63; 22.63 |
| Land evening | none (mid-level 10%) | none | nil | 0.00; 3.42 |

- **The sea.** The night's cloud is capped where CM1's is, half a kilometre high, and its shallow closure finds the
  sub-cloud layer losing energy without convection. So it gives no flux, where CM1's boundary layer vents the evening's
  stored energy at 0.014 kg m⁻² s⁻¹. The early night's deep outputs come while the boundary layer still holds the
  evening's moisture, its parcel clearing the 6–12 km barrier. The afternoon convects from
  10 km, where CM1's mass flux grows tenfold, and the evening from the ground at 3.9 km.
- **The mass flux.** Its profiles follow CM1's sampled cloudy updraughts to 30 km: 0.011–0.014 kg m⁻² s⁻¹ against
  0.014–0.022 in the afternoon, 0.010–0.020 against 0.010–0.021 in the evening.
- **The tops.** The plume ends at 35–54 km, where CM1's mass flux runs on to 60–66 km. CM1's cores, buoyant by
  0.4–0.8 K, carry it there, while the bulk plume, mixed at the bulk rates, is near neutral aloft and its kinetic energy
  drains, so the scheme heats 16–35 km more than CM1 and 35–63 km less.
- **The land.** The afternoon convects deep on 42% of outputs, the trigger stopping the rest on the stratified mean
  boundary layer (above), and PCAPEbl takes about half of those outputs' PCAPE (median).
- **The variants** (`variants_on_cm1_levels`). Tiedtke's undiluted test parcel (sub-cloud mixing nil) forms 9.7 and
  8.9 mm/day on land in the morning and afternoon, but 7.2 over the sea at night and 9.4–11.3 in the sea's afternoon and
  evening. Dropping the boundary-layer term moves the land's afternoon only to 3.7. The first round's surface parcel
  forms 2.14 at the sea's night, 5.81 in its evening and 0.97 on land in the afternoon.

**2. The GCM's site columns** (the rh 0.7 branch, years 10–14, its six sea and four land cells). The scheme convects on
none of the sea's columns at any hour and on 27% of the land's afternoon columns (16% deep), forming 3.2 mm/day where
the branch's convection rains 13.3. The branch's upper air, built by an undiluted reference, is 5.5–11 K warmer than
CM1's, too warm for a plume mixing at CM1's rates. So the first response in PlaSim would be a pause while radiation
cools the air aloft, and test 4 shows where the column then settles.

**3. Ten, fourteen and fifteen levels** (each CM1 output laid on PlaSim's ten levels and on levels.py's `l14` and
`l15`, by layer means):
- **The lowest layer** is the second round's subject (above).
- **The night cap.** On ten levels the night cloud (base 5.9 km, tops 7–8 km) and its 6–12 km cap lie inside one
  layer, between half levels at 6.4 and 13.7 km, so no level resolves either. On fourteen and fifteen, half levels at
  5.0, 7.0 and 9.0 km cut through both.
- **Elevated convection.** Fourteen and fifteen levels keep the sea afternoon's elevated deep convection on 49% of
  outputs, forming 4.6 mm/day against CM1's 6.5, where ten keep 7% and 1.1.

**4. A single column** driven by each box's own radiation, day-night table, mean ascent and surface fluxes, with the
scheme, `mkrain` at the box's constant, a saturation adjustment and a mixed boundary layer. Each box's column runs one
lunar day from the box's own air, and the GCM's site mean (5–11 K warmer aloft) runs three. Temperature less CM1's
lunar-day mean, K:

| Sea | 4 km | 8 km | 13.5 km | 19 km | 29 km | 42 km | 56 km |
|---|---|---|---|---|---|---|---|
| CM1's levels, from CM1, day 1 | +1.5 | +0.4 | −0.1 | +0.2 | 0.0 | −0.7 | −1.7 |
| 15 levels, from CM1, day 1 | +0.9 | −0.3 | −0.2 | 0.0 | −0.1 | −0.8 | −2.4 |
| 10 levels, from CM1, day 1 | +0.1 | 0.0 | 0.0 | +0.2 | 0.0 | −0.7 | −2.2 |
| 15 levels, from the GCM: start; day 3 | +3.4; +2.2 | +3.6; +2.4 | +4.1; +1.9 | +5.8; +2.7 | +8.1; +3.6 | +9.7; +3.1 | +11.3; +1.3 |
| 14 levels, from the GCM: start; day 3 | +3.4; +3.2 | +3.6; +1.7 | +4.1; +2.4 | +5.8; +3.4 | +8.1; +4.5 | +9.7; +4.2 | +11.3; +3.6 |
| 10 levels, from the GCM: start; day 3 | +3.4; +1.4 | +3.7; +2.0 | +3.7; +2.0 | +5.5; +3.1 | +8.1; +4.6 | +9.7; +4.3 | +11.3; +3.4 |
| CM1's levels, from the GCM: start; day 3 | +3.2; +8.1 | +3.8; +6.3 | +4.0; +3.5 | +5.6; +3.7 | +7.9; +4.1 | +9.5; +3.0 | +11.2; +1.7 |

- **Aloft.** At 42–56 km the GCM's excess falls by 5.4–10 K in three lunar days, 56–88% of it, and in the column
  started from CM1's air the 42–56 km air ends 0.7–2.4 K cool. That follows the low tops.
- **Near the ground.** Below 8 km the column from the GCM settles 2–3 K above CM1's on fifteen levels and 6–8 K above
  on CM1's levels, as its boundary layer holds the GCM's extra vapour.
- **Rain at the ground** (mm/day by quarter, night first; CM1 0.16, 0.00, 2.27, 2.96 at sea, 0.07, 2.92, 8.71, 1.57 on
  land):

| Single column | Night | Morning | Afternoon | Evening |
|---|---|---|---|---|
| Sea, CM1's levels, from CM1 | 0.58 | 0.06 | 2.87 | 1.38 |
| Sea, 15 levels, from CM1 | 0.56 | 0.00 | 1.47 | 2.50 |
| Sea, 14 levels, from CM1 | 0.85 | 0.00 | 2.92 | 0.97 |
| Sea, 10 levels, from CM1 | 0.08 | 0.00 | 1.74 | 2.59 |
| Land, 15 levels, from CM1 | 0.11 | 0.99 | 9.23 | 0.22 |
| Land, 14 levels, from CM1 | 0.11 | 1.10 | 9.15 | 0.21 |
| Land, 10 levels, from CM1 | 0.01 | 0.50 | 9.62 | 0.20 |

- **The land on PlaSim's grids.** A running column keeps its boundary layer mixed, and the scheme fires from the
  ground, so the afternoon rains at CM1's rate on all three grids, though the morning and evening fall short.
- **The land on CM1's levels.** The column drifts 10 K warm below 10 km within a lunar day, and the result is set
  aside as the test column's: its mixed layer and prescribed fluxes stand in for PlaSim's vertical diffusion and do
  not hold the fine grid's land.

**5. Cost** of the Python copy, counted as exponentials (each saturation, mixing and decay takes one): 77–119 a column
on the GCM's mostly quiet columns, 212–977 on the boxes' convecting columns laid on ten to fifteen levels, against the
Betts–Miller scheme's 21. The Fortran routine's timing is below.

## The Fortran routine

[massflux_fortran.py](massflux_fortran.py) holds the scheme as Fortran for PlaSim's `rainmod.f90` and the edits that put
it there, behind PlaSim's own switch: `nprc` 3 calls it in place of Kuo's scheme (1) and the simplified Betts–Miller
scheme (2), which stay as they were. The patch goes over the Betts–Miller patch, which every run's source has, with or
without that scheme's optional patches and in either order. It is not yet among the runner's patches.
- **`mkmassflux`** takes each column's state after the other processes' tendencies so far, as Kuo's and the
  Betts–Miller scheme take it. The surface fluxes are PlaSim's `dshfl` and `dlhfl` with their sign turned to upward;
  the boundary layer's forcing is the tendencies so far (radiation, vertical diffusion and the surface fluxes) with the
  dynamics' heating and moistening; the large-scale vertical velocity is `dw`; the wind speed comes from `du` and
  `dv`; land is PlaSim's mask. The scheme's tendencies join `dtdt` and `dqdt`, its rain and snow by layer go to `mkrain`
  as `dprcl` and `dprscl`, and the layers its updraught crosses become `mkclouds`' convective cloud, as Kuo's raining
  levels do. PlaSim's diffusive shallow scheme is not called: the scheme convects shallow itself.
- **`mfcolumn`** is the scheme on one column, `massflux.column` with PlaSim's phase rule step for step.
- **`mkdttgp`** gives the dynamics' temperature tendency at the grid points, transforming the spectral tendency `stt`
  as PlaSim's `mkdqtgp` transforms `sqt`, with the scaling PlaSim's own energy diagnostics give it. The spectral part
  of each step has filled `stt` with the dynamics' tendency before the physics runs.
- **`mkrain`** under `nprc` 3 evaporates in a layer at most the least running total of the handed water from that
  layer down, less what it has evaporated above (massflux.evaporate is its copy). Under `nprc` 1 and 2 that cap is
  10³⁰, so its arithmetic is as before.
- **The parameters** are `massflux.PARAMETERS` under names ending in `_mf`, each in `rainmod`'s namelist; five options
  of the Python copy (the relative-humidity factor, the IFS's deep detrainment, a 50 hPa source layer, the first
  round's excess and the other stability functions) stay in Python.
- **Its size** is about 950 lines of Fortran.

**Its tests** ([test_massflux_fortran.py](../tests/test_massflux_fortran.py), compiled with ExoPlaSim's own flags):
- **Against the Python copy.** On 16 columns, lunar and Earth-like, on PlaSim's ten and fifteen levels and on 30–40
  equal layers, deep, shallow, mid-level and quiet, the Fortran gives the type, the cloud's base and top, the base mass
  flux, PCAPE, the tendencies, the handed rain and snow and both mass fluxes to 1e-10; they agree to about 1e-15.
- **Conservation.** Its column keeps its water and moist enthalpy to 1e-9 with the latent heat of the water it hands
  on, and through the patched `mkrain` the water and heat reaching the ground equal what the scheme and `mkrain` took
  from the air, with no running total below nil but round-off.
- **The types.** A heated moist column convects deep from the ground through more than 200 hPa; a mixed layer under a
  stable layer convects shallow, its top within that layer; a saturated layer under large-scale ascent over a dry
  ground convects from there with the ascent's mass flux; a dry stable column stays quiet.
- **The other schemes.** With `nprc` 1 or 2 the patched `mkrain` prints every digit the old one prints, and every
  other edit adds text after what it finds. The patched source compiles with PlaSim's modules on fifteen levels.

**Its cost.** Time per column, the same columns on fifteen levels, both routines compiled with ExoPlaSim's flags (the
Betts–Miller scheme as the baseline runs take it: the depth test at 200 hPa, the shallow cap, and PlaSim's shallow
scheme for the columns it hands on):

| Columns | Mass-flux scheme | Betts–Miller scheme | Ratio |
|---|---|---|---|
| The GCM's sea cells, 1,024 (4 convect) | 4.7 µs | 2.4 µs | 2.0 |
| The GCM's land cells, 1,024 (95 convect) | 6.3 µs | 1.7 µs | 3.8 |
| CM1's sea columns, 236 (182 convect) | 21.7 µs | 2.3 µs | 9.3 |
| CM1's land columns, 236 (61 convect) | 11.7 µs | 2.9 µs | 4.1 |

At 2,048 columns and 17,160 steps a model year, 2–19 µs a column more than the Betts–Miller scheme's comes to 1.3–11
minutes of one core's time a year, the low end for a model as quiet as the GCM's sea cells and the high end for one
convecting as often as CM1's sea box; on a run's four processes, 20 seconds to 3 minutes. This is a timing of the
routines alone; the runs will give the model's.

## Whether to implement it, and what that would take

**The case for it.** The tests answer the approval's condition for the upper air. Mixing at CM1's measured rates holds
CM1's diluted upper troposphere and moves the GCM's toward it by 5.4–10 K at 42–56 km in three lunar days. That is the
mixing part, 5.3–6.3 K, of the gap the upper-air diagnosis traced, and the Betts–Miller scheme's mixing options could
not carry it. On fifteen levels with its surface parcel the scheme also convects from the ground where CM1 does.

**Two costs come with that change:**
- **The nights.** The sea's nights stay wetter than CM1's (0.56 mm/day at the ground in the single column on fifteen
  levels, against 0.16).
- **The deep heating.** It sits 10–25 km too low, which the column's 1–2 K cooling at 42–56 km already shows.

**The runs** follow: the routine among the runner's patches for runs that set it, ten-, fourteen- and fifteen-level
runs without it unchanged bit for bit, a smoke run, and two-year tests on fifteen levels from the fifteen-level
restart with the provisional longwave fit, judged against CM1. Any fifteen-year branch is the author's call.

## Limits

- **The data.** CM1's columns are domain means of one realization of each box, which its own convection shaped, so the
  scheme's response at once is set against CM1's realised convection. The entrainment and τ were measured on these
  boxes, so the sea quarters are not independent tests of them.
- **CM1's surface layer.** Its 50 m air and the profile above it come from CM1's simple boundary-layer scheme at 6 km
  spacing; the similarity's agreement with them tests the parcel against that model's surface layer.
- **The GCM's columns.** They are 3-day means of a branch built by another scheme. Their boundary-layer forcing lacks
  the dynamics' heating and takes radiation and the dynamics' moistening from one-step snapshots.
- **The single column.** It takes CM1's radiation and surface fluxes as given, so they do not answer its own state, and
  its boundary layer is a mixed layer.
- **Rain at the ground.** Applied on CM1's 111 levels, `mkrain`'s per-layer share evaporates nearly all the scheme's
  rain beneath the box's 6 km sub-cloud layer, which holds a deficit enough to evaporate about 13 mm/day at the sea's
  0.027. So rain at the ground is compared on PlaSim's grids and in the single column.
- **Interpretations.** PCAPEbl's T* as g H/c_p, the surface parcel's source height as the IFS's lowest level
  stretched, and the free-convection form of the stability function in the lunar surface layer.
- **The downdraught.** Its parameters are Earth's.
- **The cost.** The Fortran timing is of the routines alone, with ExoPlaSim's flags, on columns from the boxes and the
  GCM's site cells.
- **A correction.** The first round's cloud-top pressure in the depth test took the pressure of the plume's last point,
  which lies up to a quarter layer above its top; it now takes the pressure at the top's height.

## Sources

- **IFS Documentation Cy36r1, Part IV, chapter 6.** Read, from the text on the research drive
  (`session-tools/agents/climate-explain/ifs_cy36r1_part4.txt`), sections 6.1–6.11.
- **IFS Documentation Cy40r1, Part IV, chapters 3 and 6.** Read: sections 3.2 (the surface layer's stability
  functions) and 3.3 (the boundary layer's updraught, with Holtslag and Moeng's σ_w and Troen and Mahrt's excess) and
  6.3–6.4 (ecmwf.int elibrary 9204, sha256 f9fec6f0246be4f3…).
- **Bechtold and others (2014), *J. Atmos. Sci.* 71, 734–753.** Read as ECMWF Technical Memorandum 705 (2013), the
  version submitted to the journal (elibrary 8018, sha256 2e057dd36bfae7bd…): sections 2 and 5 and the appendix.
- **Bechtold and others (2008), *QJRMS* 134, 1337–1351.** Read as ECMWF Technical Memorandum 556 (elibrary 8011,
  sha256 431b13baaf05e12e…), section 2.1.
- **Tiedtke (1989), *Mon. Wea. Rev.* 117, 1779–1800.** Not read: the AMS site refused the download, and the scheme
  follows it as the IFS documentation gives it.
- **Holtslag and Moeng (1991), *J. Atmos. Sci.* 48, 1690–1698; Troen and Mahrt (1986), *Boundary-Layer Meteorol.* 37,
  129–148.** Not read; taken as the IFS documentation (Cy40r1, section 3.3) gives them.
- **Carl, Tarbell and Panofsky (1973), *J. Atmos. Sci.* 30, 788–794; Kader and Yaglom (1990), *J. Fluid Mech.* 212,
  637–662.** Not read here; the free-convection form of the stability function and its −1/3 exponent as standard
  surface-layer theory gives them.
- The PDFs and their text are in `session-tools/agents/climate-massflux/literature` on the research drive.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.massflux        # about 6 minutes on one core: ../results/gcm/massflux_offline.json
python3 -m pytest climate/tests/test_massflux.py              # the scheme's tests, under a second
python3 -m pytest climate/tests/test_massflux_fortran.py      # the Fortran routine's, about half a minute with gfortran
```
