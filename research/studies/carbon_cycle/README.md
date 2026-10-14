# The Open Moon's carbon cycle

12 October 2026. On the author's request the carbon cycle is worked out again with the project's present parameters.
This box model answers resources' [question 7](../../../resources/questions.md) and the
[CO2 study](../atmospheric_co2/README.md)'s next work: the air, the land's plants and soils with their carbonate, the
sky's aerophytes and sky reefs, the lakes' and seas' dissolved carbon with their alkalinity, the sediments, people's
stocks, and the calcium and magnesium that weathering frees, the mineral strand that carbon shares. Since 14 October it
takes the groundwater too ([springs](../springs/README.md)): the share of the weathering that reaches the waters, the
calcite the rock keeps and the carbon the filling crust takes down. It runs from the seeding of the biosphere to year
10,000, nineteen runs and a step check in about six seconds.
`python -m research.studies.carbon_cycle.run` writes [results/carbon_cycle.json](results/carbon_cycle.json) (schema
`terluna.research.carbon-cycle/1`, with its producer, the sha256 of its fifteen inputs, its evidence statement,
reading rule and assumptions); [sources.json](sources.json) records each source and how far it was read.

The parameters are the project's present ones: the corrected design case (294.9 K, the air over land at 294.3 K, the
sea surface at 297.0 K), the drainage's 265 mm a year of runoff over the land, C′'s 23.4 million km² of dry land with
the food blend, the seas and lakes at 28% water with the rain-fed lakes on 10.3% of the Moon, 8 billion people, CO2
held at 400 ppm (48.6 Pa) and methane held at 1.9 ppm, the author's choice of 13 October, or left at the soils' steady
state. The corners pair the low ends of every input, their middles and their high ends (the reading rule names them).
The author takes carbon imports as part of the plan, so the answers measure the imports and the local return paths
that keep the steady state.

## What it finds

- **400 ppm needs holding.** Left alone once the biosphere is grown (year 922–940), the air falls below the
  plants' comfortable 43 Pa in 207–447 years and to the 29 Pa that diverse ecosystems need in 1,094–2,296 years.
- **Holding it takes 0.14–0.39 Gt C a year of return or import at year 10,000.** Kilns on the carbonate as it forms
  return 0.012–0.017 Gt C of it for 5.5–10 GW, and 0.013–0.031 Gt C for 3.9–14 GW with the magnesite route; dredging
  for phosphorus returns 1–45 Mt C. The lakes' and wet margins' organic burial is the largest sink, and most of it has
  no return yet.
- **The groundwater moves the mineral strand underground.** A fifth to a half of the catchment law's weathering
  reaches the waters in the first 10,000 years (36–42% in the central corner), since nearly half of the rivers' flow is
  rain on the lakes and the filling crust takes much of the recharge. The filling crust takes 0.02–0.07 Gt C a year
  down with its water at year 10,000 and the rock keeps 0.008–0.034 as calcite, sinks no kiln can reach, so the seas
  reach three times calcite saturation only in year 3,975–6,555 and the kilns have less to return.
- **The build brings 951–1,251 Gt C by the year the biosphere is grown** (methane held, magnesium dissolved), the air's
  526 among it, and the steady state then imports 0.12–0.33 Gt C a year, 3,600–10,300 kg C/s, unless more of the
  burial comes back.
- **Burial leaves 0.36–0.86% more oxygen in the air by year 10,000**, 0.08–0.18 kPa, about 0.04–0.09% a thousand years.

### Left alone, the air falls below the plants' comfort within centuries

The free runs start from the held runs when the plants and soils reach 95% of their stocks. The kilns, the magnesite
route and the import stop; dredging and the methane hold go on. Weathering slows as the air cools under the CO2 sweep's
forcing over the GCM's 1.0–1.6 W m⁻² K⁻¹ and as the plants' help fades (GEOCARB III's law), and the land's and the
sky's stocks follow the CO2 study's leaf rate down to its compensation point, releasing their carbon.

| Corner | Grown in year | Years to 43 Pa (354 ppm) | 35 Pa (288 ppm) | 29 Pa (239 ppm) | 20 Pa (165 ppm) | CO2 in year 10,000 | Land carbon given back |
|---|---|---|---|---|---|---|---|
| Low | 940 | 447 | 1,383 | 2,296 | 4,187 | 78 ppm | 235 Gt C |
| Central | 932 | 293 | 908 | 1,530 | 2,858 | 46 ppm | 265 Gt C |
| High | 922 | 207 | 645 | 1,094 | 2,089 | 14 ppm | 268 Gt C |

The thresholds are the CO2 study's plant table: comfortable from 43 Pa, productive C3 farming at 35 Pa, diverse
ecosystems at 28–30 Pa (29 here), the glacial hard floor at 20 Pa. With methane at the soils' steady state the central
times move by under 3% (284, 900, 1,538 and 2,913 years). The sinks that drain the air are the ones a held run must
answer; the next section sizes them.

### The managed return

To hold 400 ppm with no steady import, the Moon must return its whole net sink, the carbon its reservoirs take from
the air each year. In the held runs the kilns calcine the carbonate as it forms and the deltas are dredged for the
people's phosphorus; the import is what they leave (methane held at 1.9 ppm):

| Corner | Magnesium | Year | Net sink | Kilns | Dredged | Import | Kilns' power, minimum–best kilns |
|---|---|---|---|---|---|---|---|
| Low | stays dissolved | 2,000 | 0.168 | 0 | 0.001 | 0.167 | 0 |
| | | 10,000 | 0.146 | 0.017 | 0.001 | 0.129 | 7.9–10.3 GW |
| | magnesite route | 2,000 | 0.157 | 0.011 | 0.001 | 0.144 | 3.4–4.5 GW |
| | | 10,000 | 0.140 | 0.024 | 0.001 | 0.115 | 9.3–12.2 GW |
| Central | stays dissolved | 2,000 | 0.256 | 0 | 0.009 | 0.248 | 0 |
| | | 10,000 | 0.241 | 0.012 | 0.009 | 0.221 | 5.5–7.3 GW |
| | magnesite route | 2,000 | 0.243 | 0.013 | 0.009 | 0.222 | 3.9–5.1 GW |
| | | 10,000 | 0.239 | 0.013 | 0.009 | 0.218 | 3.9–5.1 GW |
| High | stays dissolved | 2,000 | 0.395 | 0 | 0.045 | 0.350 | 0 |
| | | 10,000 | 0.385 | 0.014 | 0.045 | 0.326 | 6.8–8.9 GW |
| | magnesite route | 2,000 | 0.375 | 0.020 | 0.045 | 0.310 | 5.9–7.8 GW |
| | | 10,000 | 0.369 | 0.031 | 0.045 | 0.293 | 10.4–13.7 GW |

Gt C a year. The kilns' minimum is 4.08 GJ a tonne of CO2 from calcite and 2.57 GJ from magnesite; the best lime kilns
use up to 4.2 GJ a tonne of lime (Schorcht et al. 2013), and magnesite's kilns are scaled by the same ratio (a design
guess). The return follows two strands.

**The mineral strand.** Weathering moves the air's carbon into the waters with calcium, magnesium and sodium and
potassium: the highland soil gives them as 65%, 33% and 2% of its alkalinity, and what reaches the waters arrives at
the springs study's 45–60%, 37–61% and 2–4%, since the groundwater lays part of its calcium in the rock. Only a share of
the catchment law's weathering arrives: 36% at 1,000 years and 42% at 10,000 in the central corner (49% and 47% in the
low, 20% and 39% in the high), because nearly half of the rivers' flow is rain on the lakes and the filling crust takes
much of the recharge with its load. The rock's calcite keeps one carbon for each calcium, 0.008–0.034 Gt C a year at
year 10,000, and the water the filling crust draws carries 0.02–0.07 Gt C a year down, the land's with its soil CO2
and the floors' with the lakes' and seas' dissolved load. Calcite takes one carbon of each two and returns the other to
the air; the kilns return the one it keeps. Fresh seas start at pH 5.6 and hold weathering's carbon as bicarbonate until
they reach three times calcite saturation, in year 3,975–6,555 (central 6,555); the lakes stay below it in every
corner. Until then the kilns have only what the magnesite route gives them. Magnesium left to itself stays dissolved
with two carbons, which no kiln can reach; the magnesite route makes the lakes' new magnesium into MgCO3 as the rivers
enter, which the kilns return at 37% less energy a tonne. It lowers the steady import by 1–11% (central, 0.221 to
0.218 Gt C a year) and holds the seas' alkalinity to 1.6 meq/kg by year 10,000, against 2.5 with the magnesium
dissolved. The route's own energy, concentrating magnesium until its carbonate precipitates, lies outside the model.

- **Power.** 3.9–14 GW at year 10,000 over both routes and every corner, and at most 0.038 Gt C a year (0.14 Gt CO2)
  for 11–15 GW, the magnesite route's first years in the high corner. Brought in from outside the Moon's sunlight,
  that is under 0.001 W/m², a few parts in a million of the design sunlight; it sits beside the heat ledger's
  30–284 TW.
- **Lime.** Up to 0.08 Gt of CaO a year at year 10,000, with 0.04–0.08 Gt of MgO by the magnesite route, at most 0.13 Gt
  a year in all: 0.1–3% of the 4.2–41 Gt a year that 8 billion people's built stocks renew at resources' stocks and
  lifetimes. The kilns' return holds while the lime is kept from the air, sealed in stores or in built stocks closed to
  CO2; a quarter taking CO2 back raises the kiln load by a third, half doubles it, and lime spread on soils takes back
  all it released.

**The organic strand.** The rain-fed lakes bury 0.075–0.19 Gt C a year at Mendonça et al.'s lake rates and the wet
margins grow peat at 0.009–0.056 Gt C a year (Loisel et al. 2014): together the largest steady sink, larger than the
carbonate the kilns return. Dredging the deltas for the people's phosphorus returns 1.1–45 Mt C a year (central 8.5) at
Earth's deltaic organic content of 0.7–1.5%. Returning all of the lakes' burial that way would move 10.6–17.5 Gt of
delta sediment a year, at 4–39 GW of fuel by resources' dredging energy. The seas bury 0.003 Gt C a year at Earth's
rates for their shelves and deep floors, an upper bound for waters the rivers leave short of phosphorus.

The remainder is the import: 0.12–0.33 Gt C a year at year 10,000. Methane held at 1.9 ppm costs no carbon once there:
biology at the sources oxidizes 57–76 Tg a year to CO2, the soils take 5.4–5.5 Tg and the contactors 0.8–3.0 Tg,
passing 0.03–0.09% of the air a year.

### The carbon the build imports

By the year the biosphere is grown the Moon holds, in Gt C (held methane, magnesium dissolved):

| Corner | Grown in year | All | The air's CO2 | Plants and soils | Sky | People | Methane | Lakes' and seas' DIC | Organic sediments | Peat | Dry-soil carbonate | Rock's calcite | Filling crust |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Low | 940 | 951 | 526 | 272 | 0.3 | 2.2 | 2.2 | 52 | 69 | 8.0 | 1.1 | 7.9 | 10 |
| Central | 932 | 1,081 | 526 | 275 | 4.3 | 8.6 | 2.5 | 63 | 113 | 20 | 2.5 | 12 | 54 |
| High | 922 | 1,251 | 526 | 268 | 10.3 | 41 | 2.5 | 83 | 129 | 50 | 3.8 | 2.1 | 136 |

- **The stocks the build fills** (the air, the plants and soils, the sky, people and methane) come to 802–847 Gt C; the
  sinks take another 149–404 Gt C meanwhile, most of it into the lakes' sediments, the filling crust and the fresh
  waters' bicarbonate. Beyond the air the build brings 425–726 Gt C (central 555), at up to 3.8–4.7 Gt C a year around
  year 34, when the forests grow fastest.
- **The steady import** after the grown year runs at 0.13–0.33 Gt C a year at year 10,000 with magnesium dissolved
  (central 0.221, 7,000 kg C/s) and 0.12–0.29 with the magnesite route (central 0.218, 6,900 kg C/s): 1,253–2,999 Gt C
  between the grown year and year 10,000, two to six times the air's own CO2.
- **Variants.** Methane at the soils' steady state (25 ppm, central) puts 19–27 Gt C more into the build. Seas salted
  at delivery arrive with 261 Gt C of dissolved carbon, start at Earth's saturation and calcify from the first year,
  and then need 0.207–0.209 Gt C a year (central), their top-up water seawater as they were. Soils that build five
  times slower are grown only in year 4,655, and the build then holds 1,995 Gt C, most of the difference in the
  sinks.

### Oxygen

One O2 stays in the air for each organic carbon buried and kept (lakes, seas and peat, less the dredged return, which
is respired), and the weathering rock's FeO takes one O2 for four FeO:

| Corner | Year | Burial's O2 | Rock's O2 | Net, share of the air's O2 | Net, kPa on 21.2 |
|---|---|---|---|---|---|
| Low | 1,000 / 5,000 / 10,000 | 2.2×10¹⁴ / 1.1×10¹⁵ / 2.3×10¹⁵ kg | −3.9×10¹² / −2.0×10¹³ / −3.9×10¹³ kg | 0.035% / 0.18% / 0.36% | 0.007 / 0.038 / 0.077 |
| Central | the same | 3.8×10¹⁴ / 2.0×10¹⁵ / 3.9×10¹⁵ kg | −5.4×10¹² / −2.3×10¹³ / −4.6×10¹³ kg | 0.061% / 0.32% / 0.64% | 0.013 / 0.067 / 0.135 |
| High | the same | 5.2×10¹⁴ / 2.7×10¹⁵ / 5.4×10¹⁵ kg | −7.3×10¹² / −3.6×10¹³ / −7.7×10¹³ kg | 0.083% / 0.43% / 0.86% | 0.018 / 0.091 / 0.183 |

The rock takes back 1–2% of what burial leaves. If the build's carbon arrives as CO2, the biosphere's growth releases
another 7.4–7.5×10¹⁴ kg of O2 once, 0.12% of the air's. At 0.04–0.09% a thousand years the air's oxygen would grow by
36–86% in a million years: steady O2 needs the buried carbon returned and respired, or the steady import delivered as
reduced carbon whose burning takes up the oxygen that burial leaves (resources' redox account, questions 4 and 8).

## The model

| Reservoir or flux | Basis | Evidence | Range |
|---|---|---|---|
| Air CO2 | 400 ppm of the reference air's 1.094×10²⁰ mol: 526 Gt C | project product | held, or free from the grown year |
| Methane | sources from the gas ledger (wild 58.9 Tg a year with the land's growth, people 4.1–25.7), the soils' central uptake, the held case's 90% cut at the sources and removal to 1.9 ppm | project product and design choice | 21.6–29 ppm at the soils' steady state |
| Silicate weathering | Dessert's basalt law at 294.3 K over the land's discharge, 0.109 Gt C a year at the catchment rate; developed soils start at 1, √10 or 10 times it and age as (1 + t/age)^−0.61; the springs study's share of it reaches the waters | Earth measurement, design guess and project product | 0.05–0.21 Gt C a year reaching the waters at the seeding, 0.05–0.08 at year 10,000 |
| Cations | the springs study's export, the highland soil's Ca, Mg, Na and K less the calcium the groundwater lays in the rock | project product on a lunar measurement | Ca 45–60%, Mg 37–61%, Na and K 2–4% of the alkalinity |
| Groundwater | the springs study's central chemistry for stage 1's slow, central and fast crust (the low, central and high corners): the rock's calcite; the land's water into the filling crust with its carbon at the water table; the floors' water with the lakes' and seas' load, replaced by delivered water of the seas' own make-up | project products | rock 0.008–0.034 and filling crust 0.02–0.07 Gt C a year at year 10,000 |
| Land plants and soils | Earth's biome stocks (IPCC TAR Table 3.2, WBGU) on the rain regions: forest where the wet bands get 2 mm a day, savanna on the rest of the wet land, desert on the fog desert and polar dry land; C′'s fields and towns at cropland stocks, the pasture at savanna stocks | Earth inventory with a design mapping | 282–290 Gt C at maturity |
| Growth | plants logistic to 90% in 50 years (5 for managed land); soils toward their stock at Post & Kwon's 33 g C m⁻² a year | Earth measurement | soils approach in 372 years |
| Peat | the wet margins, 1.3–5.3% of the Moon, at 18–28 g C m⁻² a year | Earth measurement | 0.009–0.056 Gt C a year |
| Dry-soil carbonate | 0.12–0.42 g C m⁻² a year over 9.9 million km² (Schlesinger 1985) | Earth measurement | 1.2–4.1 Mt C a year |
| Sky | the sky ecology's cover (0.1–2%), class shares and standing host carbon (9.2–13.5 kg C m⁻²) | project product | 0.35–10.3 Gt C |
| People | 8 billion: food grown 1.0–1.2 Gt C a year, half of it (0.5–0.6) fixed from CO2 in the blend's closed plants and all of it returned by eating, respiration and decay within the year; half a year of it held, bodies, built stocks at 0.2–5 t C each | design guess | 2.2–41 Gt C |
| Lakes and seas | well-mixed boxes in CO2 equilibrium with the air; lakes fresh (294.3 K), seas fresh or salted (297.0 K); 96.5% of the runoff passes the lakes, which spill to the seas | textbook constants and project products | exchange within 4–19 years |
| Calcite | precipitates where fresh water stands above three times saturation (salted seas: their delivered saturation) | design guess | the seas from year 3,975 (high) to 6,555 (central); the lakes stay below it |
| Organic burial | lakes 19–48 g C m⁻² a year (Mendonça et al. 2017); seas at Earth's shelf and deep-sea rates (Duarte et al. 2005); both with the land's growth | Earth measurement | lakes 0.075–0.19, seas 0.003 Gt C a year |
| Dredging | the people's phosphorus losses (0.35–1.44 Tg P) at 449–2,083 kg of delta sediment a kilogram, 0.7–1.5% organic carbon, respired on the farms | project product and Earth measurement | 1.1–45 Mt C a year |
| Kilns and lime | all carbonate calcined as it forms; lime and magnesia kept from the air | requirement | see the return table |

The scheme steps 0.5 years to year 200, 2 years to 2,000 and 5 years to 10,000, with exact updates for the growth,
relaxation and methane, equilibrium chemistry for the waters and a root solve for the air–water split in free runs.
Carbon, calcium, magnesium, sodium with potassium, alkalinity and lime balance to 2×10⁻¹² of their stocks in every
run. Halving the step moves the answers by 0.2% and the free runs' threshold years by under 1.5 years. The waters reach
CO2 equilibrium within 4 years (lakes) and 19 years (seas) at a piston velocity of 1 m a day, short against everything
else here. The weathering law reproduces the CO2 study's twelve cases to 0.14%, and the carbonate constants their
published 25 °C values.

## Design guesses beside the outputs

| Guess | Value | What it moves |
|---|---|---|
| Developed soils' factor over catchment weathering at the seeding, and their age then | 1, √10, 10; 100–250 years | the kilns' load and the dissolved sink, most in the first millennia |
| Calcite saturation the calcifiers hold | three times (fresh); as delivered (salted) | when the kilns get their feed; where it settles |
| Rain mapping of biomes and land uses | forest from 2 mm a day | the land's 282–290 Gt C |
| Starting shares | plants and sky 1%, soils 5% | the timing of the build's import |
| People's stocks | half a year of food, 0.2–5 t C of built stock each | 2.2–41 Gt C |
| Magnesite route | all the lakes' new magnesium | 18–24% of the steady import |
| Lime kept from the air | all of it | the kiln load, by 1/(1 − f) |
| Free runs' stocks follow the leaf rate | C3 at 25 °C | how much land carbon comes back as CO2 falls |

## What the results revise in the records

- **Burial in the seas.** Resources' finding 6 and `resources/industry.py` buried carbon in the seas at lake rates,
  0.74–1.87 Gt CO2 a year, and imported 6,390–16,200 kg C/s to answer it. Earth's rates for the seas' settings give
  0.010 Gt CO2 a year; the rain-fed lakes at lake rates bury 0.27–0.69 Gt CO2 and the wet margins' peat 0.03–0.21. Both
  records now take the lakes at lake rates and the seas at Earth's marine rates, so the industry screen's burial asks
  2,370–5,970 kg C/s; the whole steady import this model leaves is 3,600–10,300 kg C/s at year 10,000.
- **The kilns.** The industry screen sizes them at 35–350 GW from the CO2 study's weathering at 26 °C; here they take
  3.9–14 GW at year 10,000 and at most 15 GW, the figures resources' findings 6 and 18 now cite.
- **Weathering.** At the design case's 294.3 K and 265 mm a year the catchment rate gives 0.40 Gt CO2 a year gross,
  56% of the CO2 study's 0.72 at 26 °C and 349 mm; resources' state notes the runoff's quarter, and the temperature
  gives the rest.
- **The land biosphere.** The CO2 study's 15–16 kg C m⁻² over 72% of the Moon, 410–437 Gt C, becomes 282–290 Gt C
  (12.1–12.4 kg C m⁻²) over the 23.4 million km² of dry land, because the corrected design case leaves 42% of the land
  under 0.5 mm of rain a day; that is 54–55% of the air's CO2 at the variable-gravity inventory of 1.93×10¹⁵ kg, which
  the CO2 study's 1.73×10¹⁵ kg understates by a tenth.
- **Peat and sediments.** The CO2 study counts them small beside weathering, under the word the author keeps for air
  lost to space; organic burial in the lakes and wet margins is the largest steady sink here, 0.08–0.24 Gt C a year.
- **The waters' carbon.** Resources' screen gives the seas and lakes at Earth's dissolved carbon 0.63–0.72 of the air's
  CO2. Fresh seas start with 3 Gt C and gain weathering's: 52–83 Gt C by the grown year and 339–416 Gt C by year
  10,000 with magnesium dissolved, while the filling crust holds 191–676 Gt C and the rock's calcite 91–240; salted at
  delivery they hold 261 Gt C at once.
- **Consistent records.** The CO2 study's 200–2,100 years to 28 Pa overlaps the free runs' 1,094–2,296 years to
  29 Pa, and ecology's finding 17 (fresh seas undersaturated until weathering supplies alkalinity) gets its timing: the
  seas pass calcite saturation by year 3,000–5,000.

- **Weathering that reaches the waters** (14 October, [springs](../springs/README.md)). The catchment law over the whole
  discharge counted rain on the lakes, nearly half of the rivers' flow, as if it had run off rock. With the springs
  study's share and make-up, the rock's calcite and the filling crust's carbon, the seas reach three times saturation
  two to four thousand years later, the kilns return a fifth to two fifths of what they did at year 10,000, for as
  little of the power, and the import at year 10,000 rises from 0.09–0.27 to 0.12–0.33 Gt C a year; left alone, the air
  takes longer to reach 29 Pa (1,094–2,296 years against 749–2,070).

## Limits

- Every rate is Earth's. Lake burial at Earth's mean lake rates falls on large, deep crater lakes (median 352 m), and
  Earth's large lakes bury less per area than its small lakes and reservoirs, so the organic strand may run high; the
  small common-level waters (5.6% of the Moon) count as seas, and at lake rates they would bury another 0.04–0.10 Gt C a
  year.
- The waters are well mixed at their surface temperatures. Deep basins mixing over centuries, the pressure that makes
  deep water more corrosive to calcite, ion pairs, magnesium carbonates without the route, and the magnesium clays
  that return CO2 are outside the model.
- Calcite precipitation at three times saturation is a design guess, and the kilns' timing rests on it; where the
  calcite settles (86% of the sea floor lies deeper than 200 m) decides how the kilns reach it.
- The land carries Earth's biome stocks on rain bands, with 8 billion people from the seeding and one highland soil
  for all dry land; the free runs' biosphere follows a C3 leaf, without C4 or CO2-concentrating organisms, and their
  climate feedback acts through temperature alone.
- The imports' form is outside the model, and it decides the oxygen; so do the soils' metallic iron and sulfide.
- The groundwater's terms are the springs study's central chemistry at three of its states, carried to the corners'
  soils and to the free runs' cooling in proportion; its own design guesses (the rock's reaction rate and surface, the
  pH at which it stops) stand behind them.

## Next work

1. **The imports' form and the redox account** (resources' questions 4 and 8): CO2, carbonate or reduced carbon for
   the build and the steady import, with the oxygen each leaves, and the soils' oxidation front by terrane.
2. **Burial by lake:** the drainage product's lakes by depth and area with Earth's depth–burial relations, and the
   small common-level waters. Minutes.
3. **Calcite where the kilns can reach it:** calcifiers in the lakes and deltas working at lower saturation, and the
   magnesite route's own energy. Minutes on this model.
4. **The organic returns:** dredging lake sediment beyond the phosphorus need, peat harvest and their costs to the
   wild, with conservation.
5. **The seas by basin** (resources' question 9): their mixing, pressure and delivered composition as scenarios.
6. **Weathering by terrane** on the atlas (resources' question 10), with the dry mare remnants.

## Running it

```
python -m research.studies.carbon_cycle.run      # writes results/carbon_cycle.json, about six seconds
python -m pytest research/studies/carbon_cycle
```

The product binds its producer and its fifteen inputs by sha256: resources' first screen, sea phosphorus and
industry; the CO2 study; the synthesis ledgers; the food study; sky ecology; the lunar-cycle register; ecology's
first screen and people products; the atlas and the drainage; the radiative-convective CO2 sweep; and the springs
study and stage 1 of groundwater. It records the
shared constants it reads by name and value. The tests regenerate it exactly and check its conservation, convergence,
constants and the sums this README states.
