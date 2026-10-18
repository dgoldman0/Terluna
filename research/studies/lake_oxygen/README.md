# The waters' oxygen in the three target climates

17 October 2026. Ecology's [question 10](../../../biosphere/ecology/questions.md) plans a one-dimensional lake model
through the lunar month for the stratification and oxygen of the deep crater lakes, and the author asked for it in each
of the three targets of [what each target climate takes](../climate_targets/README.md): the hot world near 307 K, the
resort at 300 K and the lock at 295 K. This study gives the oxygen the held fresh seas and the lakes hold at saturation
in each, and runs the [lake model](../../../biosphere/ecology/lake.py) through ten years of each target's lunar months,
from the GCM's own output, for a crater lake 352 m deep in each of three latitude bands and a productive shallow lake
3 m deep beside it, under oxygen demands from oligotrophic to eutrophic water. [run.py](run.py) writes
[results/lake_oxygen.json](results/lake_oxygen.json) (schema `terluna.research.lake-oxygen/1`) from
[lakes.py](lakes.py)'s choice of lakes ([results/lakes.json](results/lakes.json)) and [extract.py](extract.py)'s
forcing; [sources.json](sources.json) records each published value and how far it was read.

## What it finds

**Fresh seas hold a fifth more oxygen than seawater would, and each target's warmth sets how much.** At saturation under
the Moon's air, 17.46% oxygen at 1.2 atm, the held seas (0.14 g/kg) hold:

| | The hot world | The resort | The lock |
|---|---|---|---|
| Seas' surface, median | 34.6 °C | 28.7 °C | 24.2 °C |
| Oxygen at saturation | 7.06 mg/L | 7.78 mg/L | 8.43 mg/L |
| Band medians, 0–15° to 60–90° | 6.96–7.21 mg/L | 7.59–8.16 mg/L | 8.15–9.05 mg/L |
| At the seas' 99th-percentile warmth | 6.82 mg/L (36.8 °C) | 7.41 mg/L (31.6 °C) | 7.92 mg/L (27.7 °C) |
| Every water cell's 3-day means at its own pressure: median; 1st–99th percentile | 6.96; 6.53–7.25 mg/L | 7.69; 7.15–8.34 mg/L | 8.25; 7.62–9.19 mg/L |
| The GCM's lakes 3–8 km up, about 1.08 atm, median | 6.64 mg/L | 7.27 mg/L | 7.75 mg/L |
| Earth's seawater (35 g/kg) at the seas' median warmth | 5.86 mg/L | 6.41 mg/L | 6.90 mg/L |

The hot world's seas hold 16% less than the lock's and 9% less than the resort's. The air keeps Earth's oxygen partial
pressure, and at 1.2 atm the vapour takes a smaller share of it than at sea level on Earth, so fresh water here holds
0.4–0.9% more than under Earth's air at the same warmth; the held seas' salts take 0.1%. Fresh water holds 20–22% more
than seawater at the same warmth, so the hot world's seas at 34.6 °C hold 10% more than Earth's tropical seawater at
28 °C (6.44 mg/L). The fit gives 9.09 mg/L at 20 °C and 6.95 mg/L at 35 °C for fresh water under Earth's air, the
standard values. At the earlier sea medians of 34.0, 27.7 and 23.0 °C it gives 7.13, 7.91 and 8.61 mg/L, within
0.02 mg/L of the reviewer's 7.14, 7.93 and 8.63; the held seas' salts take 0.006 of that.

**The lunar night mixes only the top of a deep crater lake.** Each night's cooling and wind mix the 352 m lakes to a
median 12–15 m, and in one night in ten to 25–138 m. The deep water settles at the temperature of the coldest spell in
the record, which an ordinary night's mixed layer never reaches, so the whole column overturns only in that spell: three
to six lunar cycles in ten years, clustered within months of one another, in every lake and target. Its water is then a
mean 2.5–3.8 years from its last contact with the air, and up to 4–8 years. The lakes are oligomictic, as Earth's deep
tropical lakes are where their seasons' cooling seldom reaches the bottom (Fukushima et al. 2017): one lunar month here
differs little from the next, and only the rarest cold mixes the whole column.

| Overturns in ten years; median and 90th-percentile mixing; deep water's mean age | Mandel'shtam, 4° N | Gagarin, 18° S | Rabbi Levi, 35° S |
|---|---|---|---|
| The lock | 3; 15 and 55 m; 3.7 years | 4; 14 and 138 m; 3.3 years | 4; 14 and 71 m; 2.9 years |
| The resort | 3; 15 and 45 m; 3.5 years | 4; 14 and 131 m; 3.6 years | 4; 14 and 81 m; 2.5 years |
| The hot world | 4; 14 and 25 m; 3.5 years | 4; 13 and 30 m; 3.8 years | 6; 12 and 26 m; 2.5 years |

**Their deep water goes anoxic.** In the deepest tenth of each lake's volume:

| Deepest tenth's oxygen, mean; share of time under 0.5 mg/L | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic | 0.6–1.5 mg/L; 34–69% | 0.2–0.8 mg/L; 61–88% | 0.0–0.06 mg/L; 95–100% |
| Mesotrophic | 0.0–0.2 mg/L; 86–100% | 0.0–0.07 mg/L; 94–100% | 0.0 mg/L; 100% |
| Eutrophic, and the new reservoir | anoxic | anoxic | anoxic |

The oxygen the floor and the water below the euphotic zone take, 0.32–1.31 g per m² a day across the states, empties the
deep water within months to a few years of an overturn, and the warmer targets leave it less to start from.

**An overturn brings the anoxic water up.** When the column turns over, the whole lake takes the mix of its oxygenated
top and its anoxic deep water, and the reduced substances the anoxic water has gathered take up oxygen as they meet it.
Just after an overturn the whole lake, its surface with it, holds:

| Whole lake just after an overturn, median over the overturns | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic | 2.4–4.9 mg/L | 1.3–3.8 mg/L | 0.1–0.9 mg/L |
| Mesotrophic | 0.04–2.0 mg/L | 0.01–0.8 mg/L | 0.00–0.02 mg/L |
| Eutrophic, and the new reservoir | about 0.01 mg/L | about 0.01 mg/L | about 0.01 mg/L |

Where the overturned lake falls below 2 mg/L, its surface stays there while the overturns last, and once they stop it
climbs back above 2 mg/L within half a day to three days, as the day's warm layer forms and takes oxygen from the air;
below it the column stays anoxic. Earth's deep tropical lakes show the same: fish die when strong winds bring their
oxygen-poor deep water to the surface (Fukushima et al. 2017).

**Very deep lakes keep their deep water's oxygen for years.** Half the lake area lies in lakes deeper than 1,925 m.
Mandel'shtam's lake taken to that depth, oxygenated at the start as if just overturned, did not overturn in ten years in
the lock or the resort and turned over twice in the hot world. In the lock and the resort its deepest tenth lost
2.2–4.1 mg/L between its first two years and its last two, 0.3–0.5 mg/L a year from oligotrophic to new-reservoir water:
its floor is small for its volume, and the sediment's demand, the only one that deep, takes 11–25 years to bring it from
saturation to 2 mg/L.

**What the deep water's oxygen rests on.** At the 0–15° lake:

| Variant | Overturns in ten years | Median and 90th-percentile mixing | Deepest tenth: oligotrophic; mesotrophic |
|---|---|---|---|
| The lock, central: dissipation 10⁻¹⁰ W/kg, CO2 solutes | 3 | 15 and 55 m | 0.55; 0.00 mg/L |
| Dissipation 10⁻¹¹ W/kg | 2 | 15 and 22 m | 0.00; 0.00 mg/L |
| Dissipation 10⁻⁹ W/kg | 3 | 51 and 143 m | 3.60; 1.39 mg/L |
| No solutes | 4 | 16 and 64 m | 0.95; 0.03 mg/L |
| CO2 dissolving calcite to bicarbonate | 3 | 14 and 20 m | 0.03; 0.00 mg/L |
| Physics of eutrophic water | 4 | 12 and 122 m | 1.62; 0.04 mg/L |
| The hot world, central | 4 | 14 and 25 m | 0.00; 0.00 mg/L |
| Dissipation 10⁻¹¹ W/kg | 4 | 14 and 20 m | 0.00; 0.00 mg/L |
| Dissipation 10⁻⁹ W/kg | 12 | 39 and 260 m | 2.71; 1.18 mg/L |
| The resort, central, from the lock's weather | 3 | 15 and 45 m | 0.16; 0.00 mg/L |
| From the hot world's weather | 8 | 15 and 42 m | 0.21; 0.00 mg/L |

The dissipation below the mixed layer, which sets how fast heat and oxygen cross the thermocline, carries the most: at
ten times the central value the nights mix a median 39–51 m and oligotrophic deep water keeps 2.7–3.6 mg/L, while at a
tenth of it the deep water was still cooling by 0.006 K a pass after eight passes over the record. The solutes
respiration leaves matter little while only its CO2 is counted (0.003 kg/m³ at the bottom, worth a hundredth of a
kelvin); dissolving calcite they would add ten times as much and shorten the deeper nights' mixing from 55 to 20 m. The
resort's deep water hardly changes with the weather it is built from, although its overturns do.

**Productive shallows run hypoxic in the night, most of all in the warm targets.** In the 3 m lake, each night's lowest
surface oxygen (the median night, the worst tenth of nights in brackets):

| Lowest at night: median (worst tenth) | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic | 7.1–7.3 mg/L (6.8–7.1) | 6.4–6.6 mg/L (6.2–6.4) | 5.7–5.9 mg/L (5.6–5.7) |
| Mesotrophic | 6.0–6.6 mg/L (5.4–6.1) | 5.2–5.8 mg/L (4.8–5.3) | 4.6–5.0 mg/L (4.1–4.5) |
| Eutrophic | 1.4–2.7 mg/L (0.5–1.7) | 0.7–1.9 mg/L (0.3–0.8) | 0.4–1.1 mg/L (0.2–0.4) |
| New reservoir | 0.7–2.1 mg/L (0.2–0.9) | 0.3–1.2 mg/L (0.2–0.4) | 0.2–0.6 mg/L (0.1–0.2) |

Eutrophic shallows fall below 2 mg/L on 20–77% of nights in the lock, 57–93% in the resort and 81–99% in the hot world,
and below 0.5 mg/L on up to 12%, 35% and 58%; Gagarin's, with the weakest winds, are the lowest. The lowest point comes
from a third to four fifths of the way through the night, and by dawn the surface has climbed back to 2.8–3.9 mg/L in
the lock, 2.2–3.3 in the resort and 1.3–2.8 in the hot world, as the cooling water holds more oxygen and takes it from
the air. Oligotrophic shallows keep 82–90% of saturation at their lowest, mesotrophic ones 64–75%. Through the 354-hour
day the productive shallows' bottom water is anoxic for 294–341 hours under a surface 1.8–2.2 times saturated, so the
floor of a productive pond is anoxic most of the month in every target.

The bracket matters most in the warm targets. With Cole & Caraco's transfer and respiration's warm kinetics, the
eutrophic shallows' median night falls to 0.4 mg/L in the lock and below 0.1 mg/L in the resort and the hot world, and
mesotrophic shallows in the hot world fall to 0.25–0.42 mg/L, hypoxic in the median night. The warm kinetics raise
respiration 1.4 times at 24 °C, 2.0 at 28.7 °C and 3.2 at 34.6 °C over 20 °C with production held, an upper bound that a
community in balance with its supply would not reach.

**What it means for the targets.** Fresh water gives the hot world's seas a fifth more oxygen than seawater at their
warmth would hold, 10% more than Earth's tropical seawater holds, and the warmth still leaves them 16% below the lock's.
The lakes' oxygen turns on demand and mixing as much as on saturation, and there the warmth tells: the hot world's deep
lakes go anoxic even in oligotrophic water, an overturn leaves even an oligotrophic lake at 0.1–0.9 mg/L, and its
productive shallows run hypoxic on 81–99% of nights. The lock keeps the most: oligotrophic deep water at 0.6–1.5 mg/L,
overturns that leave 2.4–4.9 mg/L, and productive shallows hypoxic on a fifth to three quarters of nights. The resort
lies between them in each. In every target the deep crater lakes turn over rarely and their deep water runs hypoxic or
anoxic even when oligotrophic, so oxygen-breathing life in their depths needs deep water renewed by design, by oxygen
piped in, as some of Earth's eutrophic lakes and reservoirs are kept, or by driven mixing, or deep communities built for
little oxygen; and productive shallows need clear, nutrient-poor water or aeration to keep fish through the night.

## How it is set up

**The lakes.** [lakes.py](lakes.py) labels geography's rain-fed lakes as connected bodies on LOLA's 16 pixel/degree
grid, as the drainage summary does, and reproduces its 27,865 lakes and their median local depth of 352 m. Their area
lies 62% within 15° of the equator, 28% at 15–30° and 9% at 30–45°; 84% of it lies in lakes whose deepest point passes
352 m, and half of it in lakes deeper than 1,925 m. In each band the study takes the lake of at least 300 km² inside a
named crater whose deepest point lies nearest the median ([results/lakes.json](results/lakes.json)):

| Band | Crater | Centre | Area | Mean depth | Deepest | Level above the seas |
|---|---|---|---|---|---|---|
| 0–15° | Mandel'shtam | 3.97° N, 162.95° E | 581 km² | 143 m | 354 m | 1,047 m |
| 15–30° | Gagarin | 18.00° S, 151.81° E | 1,599 km² | 67 m | 363 m | 1,495 m |
| 30–45° | Rabbi Levi | 34.90° S, 23.57° E | 1,096 km² | 126 m | 351 m | 1,735 m |

The deep lake is each one's own hypsography scaled to 352 m; the Gagarin lake is mostly shallow round a deep pit, the
other two are bowls. The productive shallow lake is a pan 3 m deep beside each, at the crater lake's level, standing for
the small lakes (15,000 of the bodies are under 10 km² with median deepest points of 2–14 m) and for the shallows round
the large ones.

**The forcing.** [extract.py](extract.py) reads the two GCM runs' 3-day means at the T21 cell holding each lake (all
three cells are land in PlaSim's mask): the air at 2 m; the humidity, wind and air at the lowest level, about 0.9 km up;
the surface pressure; and the surface's net sunlight and infrared, its albedo and temperature. These are means over each
window, all but the humidity: PlaSim writes the humidity of the step that opens the window, held through it (plasim.f90
computes it only at output steps and outmod.f90 accumulates the held value), so each window takes the mean of its own
record and the next, its opening and closing steps. Downward sunlight is the net over one less the albedo, downward
infrared the net over the surface's emissivity (1 on land, 0.98 on water, as PlaSim's radiation takes them) plus the
surface's own emission. Within each mean the sunlight follows the Sun's height at the lake on the synodic clock, scaled
to the mean; the subsolar longitude fitted to every mean's Sun reproduces the GCM's mean cosine of the zenith angle to
0.0006–0.008. The other fields run linearly between the means' centres. The lowest level's wind comes down to 10 m on
the neutral profile over water with Charnock's roughness at lunar gravity, as PlaSim's surface layer assumes, after
PlaSim's 1 m/s floor. In the two years of the design case written every 6.4 hours, the 6.4-hour means' mean speed runs
19–22% above the speed of their 3-day mean vector (21–24% in the cube's ratio), so the wind is raised by the cube's
ratio at each cell. The lakes lie 0.6–1.9 km below their cells' mean ground: the air is carried down along each target's
own lapse rate over land (0.75 K/km in the lock, 0.92 in the hot world, from climate_targets' ladder), the sky's
infrared with it as the fourth power of the air's temperature, and the pressure hydrostatically; the humidity is kept
and capped at saturation. The resort takes the lock's weather and moves each mean by 0.433 of the difference between the
two runs' mean cycles by solar hour, the share of the way climate_targets uses; a second resort built from the hot
world's weather is a check.

**The model.** [biosphere/ecology/lake.py](../../../biosphere/ecology/lake.py) steps the column every 2 hours for the
deep lakes and every hour for the shallows, on fixed layers (0.5 m to 20 m depth, then growing to 8 m; 0.25 m in the
shallows) over the hypsography. Sunlight's near-infrared 58% is absorbed over 0.35 m (Paulson & Simpson 1977) and the
rest at the water's attenuation; light reaching the sloping floor heats the water over it. The surface loses infrared,
sensible and latent heat by PlaSim's own bulk scheme at 10 m: the neutral coefficient over Charnock's roughness, Louis's
stability factors, and over warm water Miller et al.'s free convection, carried to lunar gravity. The mixed layer
deepens by GLM's energy balance (Hipsey et al. 2019): convective overturn and wind stirring (C_K 0.2, C_W 0.23) against
lifting the water below and accelerating it (C_T 0.51), with lunar gravity in every buoyancy term, so the same wind
mixes deeper through the same warmth than on Earth. Below it heat and solutes diffuse at 0.2 times the dissipation over
the squared buoyancy frequency, 10⁻¹⁰ W/kg central, on top of molecular diffusion and up to 10⁻⁴ m²/s. The floor gives
the Moon's heat flow, geography's 0.015 W/m². Ice forms when the surface reaches freezing; no lake freezes in any
target, and the tests check it.

Oxygen crosses the surface at saturation for the lake's own pressure and warmth, at the surface-renewal transfer
velocity (Tedford et al. 2014, central) or Cole & Caraco's (1998); photosynthesis follows the light (tanh,
150 µmol/m²/s), scaled so that the cycle's mean production is the scenario's; respiration slows as oxygen falls toward
0.5 mg/L, its half-saturation, and what anoxic water cannot respire becomes an oxygen debt of reduced substances that
take oxygen back at a day's time scale where it returns. Respiration's CO2 makes the water denser at Lake Kivu's
coefficient, 0.284 × 10⁻³ per g/kg, 0.39 × 10⁻³ kg/m³ for each g/m³ of oxygen respired, leaving at the surface with the
gas. Each run carries sixteen oxygen columns, four trophic states by two transfer velocities by two kinetics; its
physics, the light's attenuation and the solutes follow one of them, mesotrophic water in the deep lakes and eutrophic
in the shallows; the variant runs keep the four central columns.

**The demand.** Bounded scenarios from Earth's lakes, at 20 °C:

| Trophic state | Gross production | Secchi depth, attenuation | Sediment | Water below the euphotic zone |
|---|---|---|---|---|
| Oligotrophic | 200 mg C, 0.53 g O2 per m² a day | 6 m, 0.28/m | 0.25 g O2 per m² a day | 0.07, over 25 m |
| Mesotrophic | 600 mg C, 1.6 g O2 | 3 m, 0.57/m | 0.30 | 0.30 |
| Eutrophic | 2,000 mg C, 5.3 g O2 | 1 m, 1.7/m | 0.37 | 0.53 |
| New reservoir | 2,000 mg C, 5.3 g O2 | 1 m, 1.7/m | 0.49 | 0.82 |

Production is Wetzel's (2001) ranges (oligotrophic 50–300, mesotrophic 250–1,000, eutrophic over 1,000 mg C per m² a
day); the attenuation is 1.7 over Carlson's (1977) Secchi depths. The sediment and deep water take deep European lakes
(Steinsberger et al. 2020: sediment uptake 7.9–10.6 and water-column mineralization 2.1–27.4 mmol O2 per m² a day, the
sediment's independent of trophic state) and eutrophic ones (Müller et al. 2012: 0.47–1.31 g O2 per m² a day, mean 0.90,
of which 0.37 ± 0.12 a flux of reduced substances, the water's share growing with the hypolimnion's thickness to about
25 m); the new reservoir, eutrophic for its first six to fourteen years (Paterson et al. 2019), takes Müller's highest.
The euphotic water respires four fifths of its production where it is made (a design guess), so that eutrophic
production exports what Müller's eutrophic lakes consume below. Where a lake is shallower than its euphotic zone and
25 m, the deep water's demand falls on the floor. Central demand is fixed by the trophic state's supply of organic
matter; the warm-kinetics bound gives respiration ecosystem respiration's 0.62 eV (Yvon-Durocher et al. 2012) from 20 °C
with production held.

**Spin-up.** Deep water settles at the temperature of the densest water the surface makes, so each deep lake starts 1 K
above the air's mean and saturated, and passes over the whole ten-year record in 8-hour steps until its deep water's
mean moves less than 0.005 K and 0.02 g/m³ from one pass to the next (4–7 passes), before the record is analysed. The
1,925 m lake starts from the 352 m lake's column carried down and oxygenated, as if just overturned.

## The checks

- **The model.** Its tests ([biosphere/tests/test_lake.py](../../../biosphere/tests/test_lake.py)) reproduce the
  standard saturation (9.09, 6.95 and 14.62 mg/L at 20, 35 and 0 °C) and Weiss's independent salinity factor; keep heat
  and oxygen through a freeze and a thaw to 10⁻⁹; decay a cosine mode exactly as the implicit scheme should, to 10⁻¹⁰,
  and within its first-order error of the analytic decay; deepen a convecting mixed layer into linear stratification as
  h = (2Qt/ρcΓ)^½ within 2%, and a wind-stirred one as h³ = 6 C_W u*³ t/N² within 5%, deeper at lunar gravity by the
  cube root of Earth's over it; relax a mixed slab to saturation at its transfer velocity; grow ice by conduction as
  Stefan's law has it; integrate each scenario's demand; and hold dense, CO2-laden deep water against an overturn its
  warmth alone would allow.
- **The runs.** Every run keeps its heat to 2 × 10⁻⁹ and its oxygen to 9 × 10⁻⁹ of their budgets, and no lake freezes.
  The fitted Sun reproduces the GCM's mean cosine of the zenith angle at the lake cells to 0.0006–0.008.
- **Against the GCM's own water.** A 50 m mixed slab forced as the lakes are, at the water cell nearest each lake,
  follows the GCM's slab ocean: with PlaSim's free convection at Earth's constant, 0.12–0.77 K cooler on average with
  correlations of 0.92–0.98; with it carried to lunar gravity, as the lakes take it, between 0.17 K cooler and 1.08 K
  warmer.
- **The 3-day means.** The design case's two years written every 6.4 hours drive the lakes once from those means and
  once from their own 3-day means with the wind raised as the main runs raise it: the 3-day means put the eutrophic
  shallows' median night at 2.57 mg/L against 2.72, their worst tenth at 1.46 against 1.86, and the deep lake's deepest
  tenth at 3.35 and 0.48 mg/L against 3.61 and 0.84 in oligotrophic and mesotrophic water, with six overturns in the two
  years in each. The 3-day means run the shallows' nights and the deep water a little short of oxygen.
- **The GCM's free convection.** PlaSim's fluxmod.f90 carries Miller et al.'s (1992) free-convection transfer over
  water, 0.0016 Δθ^⅓ m/s, at Earth's gravity with no gravity factor; turbulent free convection scales as (gΔθ)^⅓, so at
  lunar gravity it would be 0.55 of it. In the slab the Earth constant cools the GCM's water cells by 0.5–1.2 K against
  the constant carried to lunar gravity, a question for the climate domain.

## Assumptions and limits

- **Earth's ranges extended.** The Benson and Krause fits are stated for 0.5–1.1 atm and 0–40 °C. The Moon's 1.2 atm
  extends the pressure factor, which is thermodynamic, and the next term it omits (the Poynting correction for dissolved
  O2) is 3 × 10⁻⁴. The seas stay under 38 °C; only the lakes' sunlit surface layers pass 40 °C, to 46 °C in the hot
  world's shallows, and there the fit is extrapolated in its own variable. The Schmidt number is Wanninkhof's at 20 °C
  carried by the Stokes–Einstein scaling past his fit's 30 °C. Cole & Caraco's transfer, the GLM coefficients and the
  trophic states' rates are Earth's; Miller et al.'s free convection is carried to lunar gravity by its cube root.
- **One record, repeated.** The deep water settles to the coldest spell in the GCM's ten years. Longer records hold
  colder spells, which would set the deep water colder and the full overturns rarer, and longer isolation lets
  respiration's solutes build further.
- **Solutes.** Respiration's CO2 alone is counted. CO2 dissolving sediment calcite into bicarbonate adds about ten times
  the density (run as a bracket); methane from anoxic decomposition lightens the water and is left out, as are inflows,
  groundwater seeping through the floor and the lakes' own salts.
- **The physics left out.** Inflows and outflows, wind sheltering by the crater rim, shear and billowing at the
  thermocline, internal seiches and the sediment's own heat storage; one column stands for the whole lake.
- **Forcing.** 3-day means at T21 made hourly; the cells' air carried to the lakes' levels; the wind raised by the
  6.4-hour runs' ratio for the design case; the resort mixes two physics.
- **Biology.** Production and demand from Earth's lakes stand for engineered communities that have not been designed;
  the 2 and 0.5 mg/L lines are the usual hypoxic and anoxic marks.

## Running it

```sh
python -m research.studies.lake_oxygen.lakes                             # a few seconds: LOLA and the drainage product
climate/gcm/.venv/bin/python -m research.studies.lake_oxygen.extract      # about a minute; needs netCDF4
python -m research.studies.lake_oxygen.run                               # about 18 minutes on one core
python -m pytest research/studies/lake_oxygen biosphere/tests/test_lake.py
```

The extract reads the GCM's output in `climate/gcm/runs` and writes `research/runs/lake_oxygen/forcing.npz` with
`forcing.json` (ignored, on the research drive), recording the sha256 of every file it read; the run caches each lake in
`research/runs/lake_oxygen/cases`, so a run stopped part way resumes, and `--quick` runs two years to
`research/runs/lake_oxygen/lake_oxygen_quick.json`.
