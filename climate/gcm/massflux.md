# A mass-flux convection scheme for PlaSim: offline tests, a Fortran routine and two years in PlaSim (26 October 2026)

The author approved this work on 26 October on the condition that the scheme stand a good chance of bringing PlaSim's
lunar climate and weather nearer CM1's physics, and then approved its Fortran routine for PlaSim and two-year tests on
fifteen levels, and a starting speed for its surface parcel from the boundary layer. This record describes the scheme in
[massflux.py](massflux.py), its offline tests (results in
[massflux_offline.json](../results/gcm/massflux_offline.json)), its Fortran routine in
[massflux_fortran.py](massflux_fortran.py), which the runner applies for runs that set `--convection massflux`, those
runs (results in [massflux_l15_test.json](../results/gcm/massflux_l15_test.json)), and the starting speed's tests and
run (the third round; results in [trigger_l15_test.json](../results/gcm/trigger_l15_test.json)).

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

**In PlaSim**, over two years on fifteen levels, the scheme convects too seldom: its test parcels meet 1.4–2.4 J/kg of
inhibition in PlaSim's own columns and start with 0.5 J/kg. Its share of the rain falls from 86–87% to 42–45%, the
air near the ground warms about 2 K and gains 3–4 g/kg over the sea, and the land's afternoon rains 24–25 mm/day,
half of it large-scale. Aloft it cools the air at 42–56 km 2.3–2.6 K, a third of the way to CM1's. Of the brief's
criteria it meets only the sea's night and evening, and only with the carried longwave correction.

**The third round** started the surface parcel at the boundary layer's own convective velocity, 2.2 w*, the speed of a
grid cell's strongest thermals as CM1's 2 km box measures it. In PlaSim that speed is 0.7–1.3 m/s over the sea and nil
to 0.8 m/s over the wet land, where the IFS started at 1 m/s, and the share of rain from convection stays at 42%, which
tripped the round's stop rule. The parcel's start hardly decides the trigger in PlaSim's columns: their dry convective
layer is 0.5–4.2 km deep under a condensation level 4–6 km up, and the parcel's sub-cloud mixing spends any start across
the stable air between.

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
afternoon's columns and starts with 0.5 J/kg, the IFS's 1 m/s, and its kinetic-energy equation counts a third of
the buoyancy (its virtual mass and turbulence factors), so it crosses at most about 1.5 J/kg and reaches its cloud
on 31% of afternoons (`trigger_on_15_levels`). PlaSim's own boundary layer stays neutral under heating: its vertical
diffusion and its dry adjustment (`mkdca`) mix away any layer whose potential temperature falls with height. The
running column, whose mixed layer stands in for both, rains 9.2 mm/day on land afternoons against CM1's 8.7. So the
trigger fails on CM1's stratified mean, which the running column's boundary layer does not hold, and the explanation
pointed to no change in the trigger. The two-year runs answered otherwise: PlaSim's own columns hold 1.4–2.4 J/kg of
inhibition against the parcel too (below), so the column's mixed layer had stood in for PlaSim's boundary layer
better than PlaSim's own does.

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

## Two years in PlaSim on fifteen levels

**The patch in the runner.** `--convection massflux` ([exoplasim_run.py](exoplasim_run.py)) records the scheme with
its 33 parameters, part of the run key, writes `nprc` 3 and each parameter into `rainmod`'s namelist, and applies
massflux_fortran.py's patch to the twilight copy's source for a run that sets it, refusing the Betts–Miller scheme's
options beside it. The patch went into the shared source on 26 October, which changed `rainmod.f90`'s hash from
9db0dc61… to 7bf3ebea….

**Runs without the scheme stay as they were.** The month trials on ten, fourteen and fifteen levels, rebuilt on the
patched source, repeat the levels agent's earlier months bit for bit in all 61 output variables and their restarts.

**A smoke month** at fifteen levels with the carried correction keeps every field finite and the humidity nowhere
negative, its strongest wind 26.6 m/s. It rains 1.78 mm/day against the Betts–Miller month's 2.42, 36% of it
convective against 88%, and takes 45.9 s against 31.9 s.

**Two years** from the fifteen-level restart (`gcm-runs/A28R_301k_cold_l15_restart`), each beside the Betts–Miller run
with the same longwave correction, in turn the provisional fifteen-level refit (`A28R_301k_cold_massflux_l15` beside
`A28R_301k_cold_sbm_rh07_deep200_l15cap`) and the carried correction (`…_massflux_l15_carried` beside
`…_l15cap_carried`). Each pair differs only in the convection and the source's hash. Years 0–1, against the brief's
criteria:

| Criterion | Refit: mass-flux; baseline | Carried: mass-flux; baseline | CM1, and the test | Met (refit, carried) |
|---|---|---|---|---|
| (a) Sea's air at 41.8, 55.8 km | 256.2, 237.5; 258.5, 240.1 K | 256.9, 238.6; 259.3, 241.2 K | 251.5, 231.8; half the gap closed | no (31–33%), no (28–31%) |
| (b) Sea, midnight–sunrise | 0.68; 1.98 mm/day | 0.18; 1.48 | 0.14; at most 0.5 | no, yes |
| (b) Sea, sunset–midnight | 1.60; 2.68 | 4.17; 0.87 | 2.97; at least 1.98 | no, yes |
| (b) Sea, noon–sunset | 0.45; 1.45 | 0.97; 2.69 | 2.23; 1.1–3.3 | no, no |
| (c) Land, noon–sunset | 24.3; 17.8 | 24.7; 21.1 | 8.59; nearer than the baseline | no, no |
| (c) Land, sunrise–noon | 0.13; 0.50 | 0.04; 1.19 | 2.87; at least 1 | no, no |
| (d) Sea's lowest 0.6 km, vapour | 15.2; 12.3 g/kg | 17.2; 13.4 | 11.8 (lowest 0.5 km); no further | no, no |
| (e) Imbalance at the top | −0.45; −3.70 W/m² | +3.84; +0.63 | within 3 of the baseline | no (+3.25), no (+3.21) |
| Ground, global | 299.70; 297.81 K | 300.94; 299.03 | shift reported | +1.88, +1.91 K |
| (f) Top two levels, air | 164.2, 190.0; 165.8, 191.6 K | 162.8, 190.9; 164.1, 192.6 | reported | |
| (f) Top two levels, water | 0.021, 1.63; 0.035, 2.33 ppmv | 0.014, 2.01; 0.020, 2.93 | reported | |

The vapour is a mixing ratio and the sea's lowest layer is the fifteen-level bottom layer, compared with CM1's lowest
0.5 km at matching heights; PlaSim's 2 m air reads from the 0.3 km level and is not compared. PlaSim's free-convection
constant in its sea fluxes stays at its Earth value, 0.0016. Year 1 alone gives the same answers.

**What happens.** The scheme convects too seldom. Its share of the global rain falls from 86–87% to 42–45% (1.97 and
1.93 mm/day of rain against 2.18 and 2.15), and large-scale condensation carries the rest: half of the land
afternoon's 24–25 mm/day. The sea site's air below 10 km runs 1.5–2.3 K warmer than the baseline's and its lowest layer
holds 3–4 g/kg more vapour; the cloud cover rises from 0.35–0.36 to 0.39–0.41, and the balance at the top gains
3.2 W/m² while the global ground warms 1.9 K. Aloft the scheme heats 8–30 km where the baseline heats 19–59 km, so the
air at 42–56 km runs 2.3–2.6 K below the baseline's after two years, a third of the way to CM1's box, and the top two
levels cool 1.3–1.8 K and lose a third of their water.

**Why it convects so seldom** (`gcm_trigger`). On the runs' own 3-day mean site columns the surface test parcel meets
1.4–2.4 J/kg of inhibition (medians by quarter) below its level of free convection, and the baselines' columns give
the same. It starts with 0.5 J/kg, the IFS's 1 m/s, and its kinetic-energy equation counts a third of the
buoyancy, so it crosses about 1.5 J/kg at most and reaches its cloud on 5–49% of the sea's columns and 0–32% of the
land's, mostly as cloud 12–100 hPa
thick. PlaSim's boundary layer thus holds inhibition like CM1's mean columns, and the test parcel's starting energy
limits the scheme in PlaSim as it did on CM1's columns; what lifting energy a lunar parcel should start with, through
a boundary layer 4–6 km deep, is the open question these runs leave. The land's site is wet in the GCM, its sensible
heat running downward by 5–17 W/m² in the day's means under 170–260 W/m² of evaporation, so its surface excess is small.

**The longwave check** (`plasim_longwave check` on year 1 of both runs,
[longwave_check_massflux_l15.json](../results/gcm/longwave_check_massflux_l15.json)): on the mass-flux columns the
provisional refit holds the clear-sky longwave within 3.03 W/m² of correlated-k in its worst layer, the carried
correction within 8.76, as on the baseline's columns (3.0 and 8.2).

**The cost.** The runs took 554–586 and 536–550 s a model year on four processes against their baselines' 448–464 and
458–463 s, a fifth to a quarter more.

## The third round: the surface parcel's starting speed

The two-year runs traced the scheme's rare convection to its surface parcel's start, 0.5 J/kg of kinetic energy at the
IFS's fixed 1 m/s, against 1.4–2.4 J/kg of inhibition in PlaSim's columns. On 26 October the author approved replacing
that speed with the boundary layer's own velocity scale, tested offline and then in one two-year run, under a stop rule:
should convection's share of the rain stay below about 75%, the work would stop with a report. It stopped there.
Results are in [trigger_l15_test.json](../results/gcm/trigger_l15_test.json).

**The rule** (`start` 'boundary_layer' in [massflux.py](massflux.py); the runner's `--massflux-start boundary_layer`).
The surface deep test parcel starts at c_start w*, where w* = (g/θ_v (w′θ_v′)ₛ zᵢ)^1/3 is the boundary layer's
convective velocity (Deardorff's, as Lenschow and others 1980 define it, their eq. 12), from PlaSim's own surface
sensible and latent heat. zᵢ is the depth of the dry convective layer the parcel's source air rises through: the height
where the surroundings' virtual potential temperature first exceeds that air's (the parcel method), or where the air
saturates lifted dry, whichever is lower, each linear in height between levels. Where the surface buoyancy flux does
not rise, w* and the starting speed are nil; the rule has no floor. Elevated and mid-level departures keep the IFS's
1 m/s, and the shallow test keeps its σ_w.

**The statistic.** The test parcel decides whether any of a grid cell's thermals reach their level of free convection,
so it stands for the strongest of them. A T21 cell on the Moon is about 170 km across and holds several hundred thermals
if they lie about a mixed-layer depth apart, which puts its strongest near the 99.9th percentile of the mixed layer's
vertical wind (an interpretation). Lenschow and others' fit to their aircraft data over the sea, which also fits the
Minnesota data over land and Deardorff's simulations, σ_w²/w*² = 1.8 (z/zᵢ)^2/3 (1 − 0.8 z/zᵢ)² (their eq. 13), peaks at
0.68 w* at 0.31 zᵢ, where a Gaussian's 99.9th percentile, 3.09 σ_w, is 2.1 w*; the convective boundary layer's upward
skew lengthens that tail. The IFS's own boundary-layer plume takes the mean of a Gaussian's top 5%, 2.05 σ_w (Cy40r1,
section 3.3), for the thermals whose flux it carries; the trigger asks about the strongest.

**What CM1 measures** (`cm1_thermals.py` on every 3-D output; w* from CM1's surface fluxes and the parcel method's depth
on its domain-mean column; means over each period's outputs):

| CM1 box and period | w*, m/s | zᵢ; condensation level, km | Peak σ_w / w* | 99.9th percentile at 0.3–0.9 zᵢ, / w* | 99.9th percentile below 16 km: all air; clear air, m/s | 2.2 w*, m/s |
|---|---|---|---|---|---|---|
| Land, 30° before noon to noon | 1.15 | 6.1; 5.7 | 0.36 | 1.36 | 2.63; 1.84 | 2.52 |
| Land, noon to 30° after | 1.17 | 5.6; 5.3 | 0.37 | 1.76 | 2.91; 1.94 | 2.57 |
| Sea, midnight–sunrise | 0.70 | 5.9; 5.3 | 0.46 | 1.65 | 1.10; 1.10 | 1.55 |
| Sea night, 6 km grid | 0.70 | 6.0; 5.7 | 0.48 | 1.73 | 1.11; 1.11 | 1.55 |
| Sea night, 3 km grid | 0.70 | 5.9; 5.7 | 0.62 | 2.11 | 1.35; 1.35 | 1.54 |
| Sea night, 2 km grid | 0.70 | 5.8; 5.7 | 0.68 | 2.28 | 1.46; 1.45 | 1.54 |

At 2 km spacing CM1's sea night has Lenschow's spread, 0.68 w* at its peak, and the 99.9th percentile of its vertical
wind in the mixed layer's upper part is 2.28 w* (2.11 at 3 km). The 6 km boxes, the land and sea boxes among them,
resolve thermals a mixed layer wide only coarsely, and their spread peaks at 0.36–0.48 w*. c_start is therefore 2.2,
from the 2 and 3 km boxes, within the range the Gaussian and the skew give.

**The check against CM1.** On CM1's own w* the rule starts the parcel at 2.5–2.6 m/s at the land's onset by day, where
CM1's 99.9th percentile below 16 km is 2.6–2.9 m/s in all air and 1.8–1.9 m/s in clear air, and at 1.55 m/s in the
sea's night, where CM1 measures 1.10 m/s at 6 km spacing, 1.35 at 3 km and 1.46 at 2 km.

**Offline on the GCM's columns** (`gcm_trigger_start.py`; the runs' 3-day mean site columns over years 0–1; each cell
gives the new start's median in m/s, then the share of surface test parcels reaching their cloud under the old rule and
the new):

| Columns | Night | Morning | Afternoon | Evening |
|---|---|---|---|---|
| Mass-flux run, sea | 1.01; 49 → 49% | 0.89; 7 → 7% | 0.79; 7 → 7% | 0.96; 24 → 24% |
| Mass-flux run, land | 0.00; 1 → 0% | 0.73; 17 → 17% | 0.36; 32 → 30% | 0.00; 0 → 0% |
| Mass-flux run with the carried correction, sea | 0.80; 23 → 23% | 0.73; 6 → 6% | 0.69; 5 → 5% | 0.76; 16 → 15% |
| Mass-flux run with the carried correction, land | 0.00; 0 → 0% | 0.66; 9 → 9% | 0.00; 18 → 18% | 0.00; 0 → 0% |
| Baseline, sea | 1.20; 36 → 37% | 0.99; 4 → 5% | 0.95; 6 → 6% | 1.27; 17 → 17% |
| Baseline, land | 0.00; 0 → 0% | 0.82; 9 → 8% | 0.00; 34 → 32% | 0.00; 0 → 0% |
| Baseline with the carried correction, sea | 0.86; 18 → 18% | 0.74; 3 → 2% | 0.77; 2 → 2% | 0.85; 3 → 3% |
| Baseline with the carried correction, land | 0.00; 0 → 0% | 0.61; 4 → 4% | 0.00; 34 → 32% | 0.00; 0 → 0% |

PlaSim's w* over the sea is 0.32–0.58 m/s (medians), against CM1's 0.70, because its dry convective layer is 0.5–4.2 km
deep under a condensation level 4–6 km up, where CM1's sea night is mixed to 5.9 km. Over the land site PlaSim's surface
buoyancy flux rises only by day and weakly, its sensible heat running downward under 170–260 W/m² of evaporation, so w*
is nil in about half the afternoon means. The new start ranges from nil to 1.3 m/s, and the share reaching cloud moves
by at most two points in any quarter. The starting-speed run's own columns give the same.

**Why the start matters so little** (`start_sensitivity.py`; the share of surface test parcels reaching their cloud, %,
at fixed starting speeds with the parcel's sub-cloud mixing at its 0.4 of the base rate, and at 1 m/s with none):

| Columns, quarter | 0 | 0.5 | 1 | 1.5 | 2 | 3 | 5 m/s | 1 m/s, no sub-cloud mixing |
|---|---|---|---|---|---|---|---|---|
| Mass-flux run, sea night | 47 | 48 | 49 | 51 | 53 | 62 | 77 | 71 |
| Mass-flux run, sea evening | 22 | 22 | 24 | 26 | 27 | 33 | 44 | 65 |
| Mass-flux run, sea afternoon | 6 | 6 | 7 | 7 | 8 | 10 | 13 | 37 |
| Mass-flux run, land afternoon | 30 | 30 | 32 | 35 | 39 | 48 | 64 | 58 |
| Mass-flux run, land morning | 17 | 17 | 17 | 19 | 21 | 25 | 30 | 51 |
| Baseline, sea evening | 17 | 17 | 17 | 17 | 18 | 22 | 32 | 83 |

Between rest and 2 m/s the share moves by at most nine points. The parcel's sub-cloud mixing, 1.12×10⁻⁴ m⁻¹, also drags
on its kinetic energy, which then decays by a factor e every 2.3 km, so whatever start a thermal brings is mostly spent
across the 2–5 km of stable air PlaSim holds below the condensation level, and the trigger turns on the parcel's
buoyancy near that level. The same columns without the sub-cloud mixing reach cloud 1.4–5 times as often.

**The patch.** [massflux_fortran.py](massflux_fortran.py)'s `START_PATCHES` go over the mass-flux patch for a run that
sets the rule: `mfstart` gives w*, and `nstart_mf` and `cstart_mf` join `rainmod`'s namelist, off by default. A parcel
may now start at rest, so the organized detrainment's ratio of kinetic energies is taken only where the energy below is
positive, as it always is under the IFS's start. The patch went into the shared source on 26 October, changing
`rainmod.f90`'s hash from 7bf3ebea… to 42122b3c…. Its tests
([test_massflux_fortran.py](../tests/test_massflux_fortran.py)):
- with the rule, the Fortran gives the Python copy's type, cloud, base mass flux, tendencies, rain, snow and mass fluxes
  to 1e-10 on 40 columns (the earlier 16, and their 12 lunar ones again under the GCM's sea fluxes and its wet land's),
  and keeps each column's water and moist enthalpy to 1e-9;
- without the rule, the patched routine prints every digit the old one prints on the same columns;
- w* and its depth follow the parcel method and stop at the condensation level, and a falling flux starts the parcel at
  rest, with no division by its energy on a saturated lowest level.

The month trials without the rule, run on the patched source, repeat the earlier months bit for bit in all 61 output
variables and their restarts, with the Betts–Miller scheme and with the mass-flux scheme on fifteen levels. A smoke month
with the rule keeps every field finite and the humidity nowhere negative, and rains 1.69 mm/day, 33% of it convective.

**Two years with the rule** (`A28R_301k_cold_massflux_wstart_l15`, the earlier mass-flux run's flags with the rule and
the provisional refit; years 0–1, criterion (e) against the baseline):

| Criterion | Starting-speed run | Earlier mass-flux run | Baseline | CM1, and the test | Met |
|---|---|---|---|---|---|
| (a) Sea's air at 41.8, 55.8 km | 256.22, 237.55 K | 256.21, 237.51 | 258.53, 240.12 | 251.49, 231.76; half the gap closed | no (33%, 31%) |
| (b) Sea, midnight–sunrise | 0.99 mm/day | 0.68 | 1.98 | 0.14; at most 0.5 | no |
| (b) Sea, sunset–midnight | 2.99 | 1.60 | 2.68 | 2.97; at least 1.98 | yes |
| (b) Sea, noon–sunset | 0.86 | 0.45 | 1.45 | 2.23; 1.1–3.3 | no |
| (c) Land, noon–sunset | 24.05 | 24.29 | 17.84 | 8.59; nearer than the baseline | no |
| (c) Land, sunrise–noon | 0.05 | 0.13 | 0.50 | 2.87; at least 1 | no |
| (d) Sea's lowest 0.6 km, vapour | 15.09 g/kg | 15.23 | 12.34 | 11.8 (lowest 0.5 km); no further | no |
| (e) Imbalance at the top | −0.38 W/m² | −0.45 | −3.70 | within 3 of the baseline | no (+3.32) |
| Ground, global | 299.72 K | 299.70 | 297.81 | shift reported | +1.91 K |
| Rain; its convective share | 1.96 mm/day; 42.4% | 1.97; 42.2% | 2.18; 86.0% | the stop rule: below 75% | tripped |
| (f) Top two levels, air | 164.17, 190.00 K | 164.15, 189.97 | 165.82, 191.56 | reported | |
| (f) Top two levels, water | 0.0210, 1.64 ppmv | 0.0210, 1.63 | 0.0349, 2.33 | reported | |

Year 1 alone gives a convective share of 43.7%, against 43.0% in the earlier run. The run repeats the earlier mass-flux
run's global state to within 0.1 W/m² at the top of the air, 0.02 K at the ground and 0.2 points of convective share. At
the sea site its evening rains 2.99 mm/day, CM1's rate, against 1.60 before (2.8 and 3.2 in years 0 and 1, against 1.3
and 1.9), and its afternoon 0.86 against 0.45; the land site's afternoon still rains 24 mm/day, three-fifths of it
large-scale. The longwave check (`plasim_longwave check` on year 1,
[longwave_check_trigger_l15.json](../results/gcm/longwave_check_trigger_l15.json)) holds the provisional refit within
3.06 W/m² of correlated-k in its worst layer, as on the earlier runs (3.0–3.03). The run took 558 and 561 s a model year
on four processes, against the earlier mass-flux run's 586 and 554 and the baseline's 464 and 448: the rule adds nothing
measurable to the scheme's fifth to a quarter more than the baseline.

**Where the trigger's limit lies.** The stop rule tripped, and this round's evidence places the scheme's rare convection
in PlaSim's boundary layer. PlaSim keeps a dry convective layer 0.5–4.2 km deep over the sea under a condensation level
4–6 km up, so a surface parcel must cross 2–5 km of stable air while it mixes, where CM1's sea night is mixed to 5.9 km
and its condensation level lies within that layer; the same shallow layer holds PlaSim's w* to 0.3–0.6 m/s against CM1's
0.70. A candidate cause, untested: PlaSim's vertical diffusion (`vdiff`) keeps Earth's asymptotic mixing length, 160 m,
and counts heights at a fixed 250 K, while its coefficient in σ goes as g³. For the same wind shear and Richardson
number it mixes PlaSim's interfaces 0.6 and 3.1 km up 7–36 times more slowly at the Moon's gravity than at Earth's (a
reading of its arithmetic), and its dry adjustment mixes only where the air turns superadiabatic. A second part lies in
the scheme: its test parcel mixes in sub-cloud air at 0.4 of the deep base rate, a rate CM1 does not measure (the
parameter table above), and that mixing decides more of the trigger than the start does.

## Whether to implement it, and what that would take

**The case for it.** The tests answer the approval's condition for the upper air. Mixing at CM1's measured rates holds
CM1's diluted upper troposphere and moves the GCM's toward it by 5.4–10 K at 42–56 km in three lunar days. That is the
mixing part, 5.3–6.3 K, of the gap the upper-air diagnosis traced, and the Betts–Miller scheme's mixing options could
not carry it. On fifteen levels with its surface parcel the scheme also convects from the ground where CM1 does.

**Two costs come with that change:**
- **The nights.** The sea's nights stay wetter than CM1's (0.56 mm/day at the ground in the single column on fifteen
  levels, against 0.16).
- **The deep heating.** It sits 10–25 km too low, which the column's 1–2 K cooling at 42–56 km already shows.

**The two-year runs** show a third cost: in PlaSim the scheme convects too seldom, its test parcels held back by the
inhibition PlaSim's columns carry, so large-scale condensation takes half the rain, the near-ground air warms and
moistens, and the balance at the top gains 3.2 W/m². Its cooling of the air aloft goes a third of the way to CM1's in
two years. The third round's starting speed leaves all of this as it was, and places the trigger's limit in PlaSim's
shallow boundary layer and the test parcel's sub-cloud mixing. Any fifteen-year branch is the author's call.

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
  GCM's site cells; the runs' seconds a year include whatever else the machine ran.
- **The two-year runs.** Each is one realization of the first two years after the remap, so they show the direction
  and early size of the scheme's effect; the inhibition its parcels meet is measured on 3-day means.
- **A correction.** The first round's cloud-top pressure in the depth test took the pressure of the plume's last point,
  which lies up to a quarter layer above its top; it now takes the pressure at the top's height.
- **The third round.** The statistic the parcel stands for, the grid cell's strongest thermals at the 99.9th percentile,
  is an interpretation, and c_start comes from one CM1 night at sea on 2 and 3 km grids. The depth zᵢ is read from 15
  levels, linear between them. The offline tests take 3-day means, and the two-year run is one realization; the reading
  of `vdiff` as the shallow boundary layer's cause is untested.

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
- **Lenschow, Wyngaard and Pennell (1980), *J. Atmos. Sci.* 37, 1313–1326.** Read (journals.ametsoc.org, sha256
  21660f4f70942b12…): section 3, its mixed-layer scales (eq. 12, w*) and the vertical wind's variance (eq. 13). Its PDF
  is in `session-tools/agents/climate-trigger/literature` on the research drive.
- **Deardorff's convective velocity scale.** Taken as Lenschow and others define it; Deardorff's own papers not read.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.massflux        # about 6 minutes on one core: ../results/gcm/massflux_offline.json
python3 -m pytest climate/tests/test_massflux.py              # the scheme's tests, under a second
python3 -m pytest climate/tests/test_massflux_fortran.py      # the Fortran routine's, about half a minute with gfortran
M=/media/projectspace/terluna-research/session-tools/agents/climate-massflux2
$M/run_massflux_tests.sh month10                 # and month14, month15: the months without the scheme, for the bit check
climate/gcm/.venv/bin/python $M/bit_check.py A28R_301k_cold_sbm_rh07_deep200_month_mfsource \
    A28R_301k_cold_sbm_rh07_deep200_month_recheck …                               # every variable and the restarts
$M/run_massflux_tests.sh smoke15                 # one model month with the scheme; then smoke_check.py
$M/run_massflux_tests.sh test15                  # two years with the provisional refit, about 20 minutes on 4 workers
$M/run_massflux_tests.sh carried15               # two years with the carried correction
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave check --case A28R_301k_cold_massflux_l15:1 \
    --case A28R_301k_cold_massflux_l15_carried:1 \
    --correction provisional_l15=climate/results/gcm/longwave_A28R_301k_cold_sbm_rh07_deep200_l15cap_month_refit.json \
    --correction carried_301k_l15=climate/results/gcm/longwave_A28_corrected_970_lw_ice_l15.json \
    --output climate/results/gcm/longwave_check_massflux_l15.json                # about a minute and a half
climate/gcm/.venv/bin/python $M/gcm_trigger.py A28R_301k_cold_massflux_l15 …     # the parcel on the runs' columns
climate/gcm/.venv/bin/python $M/massflux_l15_analysis.py > $M/massflux_l15_analysis.txt   # the results file
T=/media/projectspace/terluna-research/session-tools/agents/climate-trigger       # the third round
climate/gcm/.venv/bin/python $T/cm1_thermals.py                # CM1's vertical wind against its w*, minutes on one core
climate/gcm/.venv/bin/python $T/cm1_thermals_summary.py        # the table above
climate/gcm/.venv/bin/python $T/gcm_trigger_start.py A28R_301k_cold_massflux_l15 …   # the offline test, about a minute
climate/gcm/.venv/bin/python $T/start_sensitivity.py A28R_301k_cold_massflux_l15:sea:night …
$T/run_trigger_tests.sh ws_smoke15               # one model month with the rule; it applies the patch
$T/run_trigger_tests.sh ws_month15               # and ws_mfmonth15: the months without it, for the bit check
climate/gcm/.venv/bin/python $T/bit_check.py A28R_301k_cold_sbm_rh07_deep200_l15cap_month_wssource \
    A28R_301k_cold_sbm_rh07_deep200_l15cap_month_fcsource A28R_301k_cold_massflux_l15_month_wssource \
    A28R_301k_cold_massflux_l15_month_fcsource                                     # every variable and the restarts
$T/run_trigger_tests.sh ws15                     # two years with the rule, about 19 minutes on 4 workers
climate/gcm/.venv/bin/python $T/trigger_analysis.py > $T/trigger_analysis.txt   # trigger_l15_test.json
```
