# The lakes' and seas' methane at the resort's 301 K, and the biology that eats it

18 October 2026. The author set the resort at 301 K on 18 October, climate_targets' recommendation within the 300–302 K
band the summit tower's winter allows ([the register](../../decisions.md#the-target-climate-17-october-2026)), and
approved this study after [the lakes' anoxic depths](../lake_anoxia/README.md) reckoned the air's methane with every
lake at the flux of its 352 m lakes, while lakes over 1,000 m deep hold 70% of the lake area, and left the biology of
methanotrophs designed for sudden overturns undone. Here methane is a tracer of its own: [tracer.py](tracer.py) carries
it through the columns the lake model (`biosphere/ecology/lake.py`) saved, with [columns.py](columns.py) reading them
and [biology.py](biology.py) holding the methane-eating biology's ranges; [run.py](run.py) writes
[results/lake_methane.json](results/lake_methane.json) (schema `terluna.research.lake-methane/1`); the published values
and how far each was read are in [sources.json](sources.json).

> **Status.** Traced on the final lake model's columns: lake_anoxia's regenerated 352 m and 1,925 m lakes in the three
> targets, the lake kinds' 301 K ladder, giants and meromictic lakes, lake_oxygen's mixing-forms columns and the seas'
> 90 columns. The figures below are this tracer's arithmetic on them.

## What it finds

At the resort's 301 K (0.52 of the way from the design case to the warm Moon, absorbing 218.6 W/m²), with Earth's lake
methanotrophs at Earth's central rates and the methane share's central 0.5, the Moon's lakes (3.92 million km²) give
0.15 Tg of methane a year nutrient-poor, 3.5 mesotrophic, 20 eutrophic and 21 as new reservoirs, and the held fresh seas
(10.6 million km²) 0.42, 4.3, 14 and 37, their shelves to 200 m at the lakes' rungs of the same depth. Nutrient-poor
lakes and seas together, 0.57 Tg a year, hold the air at 1.01–1.75 ppm, under the held 1.87 ppm whichever community the
water keeps, with 0.35 Tg a year of room at the ledger's high-cost end (0.23 with the poor community) and 2.5 at its
low-cost end. That room is thin. At the methane share's 0.8 the same waters give 3.6 Tg a year and put the air at
2.04–2.78 ppm, needing 0.5–2.6 Tg a year of removal; under Shih and Bouffard's form of the deep mixing Earth's community
gives 1.7 Tg a year and 1.39–2.13 ppm, needing 0.77 Tg a year at the high-cost end; the 3,500 m nutrient-poor giant,
still losing oxygen, would add about 0.29 Tg a year once anoxic; against these, the seas' shelves at the lakes' rungs
lean high by up to 0.15 Tg a year. Mesotrophic waters, 7.8 Tg a year, put the air at 3.5–4.2 ppm, needing 4.7–6.9 Tg a
year of removal, and eutrophic waters and new reservoirs at 12.5–21 ppm. In nutrient-poor lakes the floors' bubbles in
shallow water carry most of the methane; in productive ones the giants' anoxic deep water does, through oxyclines whose
methanotrophs oxidise 88–99% of what reaches them. These are arithmetic on one lake model's saved columns with Earth's
measured rates carried to engineered communities; no lunar water has been measured.

**By depth.** Emission per m² of lake surface, mol a year, under Earth's community, each rung's bands weighted by its
lakes' area (the 1,925 m giants run in all three bands nutrient-poor and mesotrophic, the other giants in the 0–15°
band):

| Rung (the lakes whose deepest point lies nearest) | Area share | Nutrient-poor | Mesotrophic | Eutrophic | New reservoir |
|---|---|---|---|---|---|
| 5–50 m (to 71 m) | 4.4% | 0.016–0.076 | 0.10–0.17 | 0.37–0.61 | 0.72–1.10 |
| 100 m (71–188 m) | 5.9% | 0.0067 | 0.063 | 0.19 | 0.39 |
| 352 m (188–593 m) | 12.1% | 0.0021 | 0.011 | 0.044 | 0.30 |
| 1,000 m (593–1,388 m) | 14.6% | 0.00055 | 0.0040 | 0.36 (the 1,925 m giant's) | 0.30 (the 352 m lake's) |
| 1,925 m (1,388–2,596 m) | 25.5% | 0.00042 | 0.061 | 0.36 | 0.30 (the 352 m lake's) |
| 3,500 m (beyond 2,596 m) | 37.5% | 0.00016 | 0.074 | 0.36 (the 1,925 m giant's) | 0.30 (the 352 m lake's) |

To 100 m the nights mix the water through and every floor's bubbles carry almost all the emission; the methane share's
0.3 leaves the floors below their bubble point and cuts it ten to a thousand times, and its 0.8 raises it two to seven
times. At 352 m the oxic deep water lets only the margins' bubbles out unless the lake is a new reservoir, whose anoxic
deep water makes 2.6–3.7 mol/m² a year and sends 0.075–0.26 across the surface after its oxycline has oxidised 95–97%.
The mesotrophic and eutrophic giants turn anoxic below their seasonal mixing (1.3–2.4 mol/m² a year made) and emit what
their oxyclines pass, 0.022–0.36 mol/m² a year. The nutrient-poor giants emit only their margins' bubbles: the 1,000 m
and 1,925 m ones keep oxic deep water, and the 3,500 m one, still losing oxygen, makes 0.13 mol/m² a year in its
anoxic deepest layers and oxidises all of it. The 1,000 m
mesotrophic giant keeps oxic deep water too and emits about a fifteenth of the deeper giants'. No giant was run
eutrophic at 1,000 or 3,500 m or as a new reservoir at any depth; those rungs take the deepest lake run in that state,
which for the new reservoirs' 77% of the area is the 352 m lake's flux, a bound for every newly flooded lake whose
direction is open (below).

**The lakes Moon-wide**, Tg a year (each rung's flux over its lakes' area; the lakes between the 100 m and 352 m rungs
split at the lake kinds census's ventilated depth):

| Lakes, Tg a year | Earth's community | Designed–poor | No methanotrophs | Ventilated depth at 100–352 m | Shih and Bouffard's form | lake_anoxia's uniform scale-up |
|---|---|---|---|---|---|---|
| Nutrient-poor | 0.147 | 0.147–0.254 | 0.88 | 0.134–0.164 | 0.75 | 0.135 |
| Mesotrophic | 3.47 | 1.76–16.2 | 46.6 | 3.27–3.61 | 8.12 | 0.70 |
| Eutrophic | 20.2 | 10.2–65.0 | 114 | 19.6–20.5 | 27.1 | 2.79 |
| New reservoir | 20.6 | 13.7–125 | 185 | 20.2–20.8 | 42.6 | 18.5 |

The rungs of 1,000 m and deeper carry 78% of the area and 79–87% of the mesotrophic and eutrophic methane; the shallow
lakes to 50 m, 4.4% of the area, carry 62% of the nutrient-poor lakes' methane. lake_anoxia's uniform scale-up, every
lake at its 352 m lakes' flux, comes out close for nutrient-poor lakes and new reservoirs and five to seven times low
for mesotrophic and eutrophic ones, whose giants turn anoxic where the 352 m lake stays oxic. The lane without
methanotrophs is its traced transient after four to ten decades; all it makes would leave once its store filled.
Splitting the 100–352 m lakes at either rung, the census's bracket on its self-ventilated lakes (226,500–629,700 km²,
central 403,700), moves the nutrient-poor total by 9–12% and the others by 1–6%. Under Shih and Bouffard's form the
1,925 m nutrient-poor giant turns anoxic below its seasonal mixing (0.74 mol/m² a year made) and the nutrient-poor total
rises fivefold; the productive states rise 1.3–2.3 times. The borrowed fluxes lean as follows. The eutrophic 1,000 m and
3,500 m rungs take the 1,925 m giant's 0.36, run only at 0–15° and still settling, and carry 11.9 of the eutrophic
lakes' 20.2 Tg a year; the 1,000 m borrowing is probably an upper value, since the mesotrophic giant at that depth keeps
oxic deep water, and at the 352 m lake's 0.044 the eutrophic total would be about 2.9 Tg lower. The giants run at 0–15°
alone stand for every band, and the 0–15° runs sit at the high end of the mesotrophic giants' spread (1,925 m: 0.082
against 0.022 at 15–30° and 0.072 at 30–45°); the 3,500 m mesotrophic giant carries 1.73 of 3.47 Tg a year, about 9%
less if its 15–30° part behaved like the 1,925 m giant's. The new reservoirs' giant rungs take the 352 m lake's 0.30, a
bound for every newly flooded lake, since the census counts no giant reservoirs; lake_anoxia's 1,925 m new reservoirs
emit 0.07–0.10 by their fifth decade, so the bound's direction is open. In the kinds study's record the 3,500 m
nutrient-poor giant, 37.5% of the lake area, is losing oxygen (its deepest tenth from 0.72 to 0.53 mg/L over five
decades) and turns anoxic within 12.6 years of its record's start, so its methane is a lower value: settled anoxic,
as the 1,925 m nutrient-poor giant does under Shih and Bouffard's form (0.012 mol/m² a year more), it would add about
0.29 Tg a year to the nutrient-poor lakes. The 1,000 m giants and the 1,925 m nutrient-poor giants gain oxygen
slowly, their methane as traced or slightly high; the mesotrophic and eutrophic giants are anoxic throughout
(kinds.giants).

**The seas**, Tg a year, each basin weighted by what its columns stand for and its water taken by local depth from its
own shelf shares (seas.json; the water it stands for beyond its own at the shares the seas' whole leaves once every
basin's own shelves are counted): to 5 m at the 3 m pan, 5–30 m at the 10 m pan, 30–200 m at the lakes' ladder rungs
nearest that depth in the basin's band, and deeper at its typical deep column:

| Seas, Tg a year | Earth's community | Designed–poor | Coldest columns | Shih and Bouffard's form | Wind share 0.2%; 0.4% | Deep columns from 30 m |
|---|---|---|---|---|---|
| Nutrient-poor | 0.424 | 0.423–0.440 | 0.425 | 0.94 (0.71–1.29) | 0.445; 0.424 | 0.278 |
| Mesotrophic | 4.33 | 3.18–17.3 | 6.90 | 8.8 (7.1–11.5) | 8.03; 3.30 | 3.60 |
| Eutrophic | 13.8 | 8.87–41.2 | 18.8 | 20.8 (20.0–22.1) | 17.9; 11.9 | 11.6 |
| New reservoir | 36.6 | 22.6–119 | 59.5 | 53.3 (49.1–56.1) | 43.1; 31.7 | 33.6 |

The water to 30 m, 2.2% of the seas, gives 0.25 of the nutrient-poor seas' 0.42 Tg a year by its floors' bubbles, and
the shelves from 30 m to 200 m, 12.2%, 0.15 at the lakes' rungs. Those rungs carry their own shallow margins' bubbles,
so the shelves there lean high; with the deep columns' flux there instead the seas give 0.28, 3.6, 11.6 and 33.6 Tg a
year. The deep basins' mesotrophic and eutrophic methane comes mostly from Jules Verne, standing for the crater seas
(24% of the seas: 2.2 and 7.0 Tg a year), and Ingenii (20%: 0.21 and 2.5), whose deep water turns anoxic; Procellarum
(30%) overturns each year and gives 0.31 and 0.83. Under Shih and Bouffard's form the two basins run under it take their
own rows for their deep water and every other basin its own emission plus the change per m² theirs show; the bracket is
each basin's change alone. The seas' study reads several of these columns as unsettled at their analysis step, and the
methane follows: Procellarum's Shih and Bouffard column holds the decade after its one overturn in a century, its deep
water just renewed, so its methane is a lower value; Ingenii's coldest column and Crisium's at the wind share's 0.4% are
still losing oxygen (methane a lower value); Fecunditatis's at 0.2% is still gaining it (an upper value); the four
columns under the lake number's mixing and the fixed dissipation are cut off from the surface and still ageing, their
methane a lower value. In 13 of the 20 deep columns the methane store of Earth's community in some trophic state still
moved by 0.5–12% in the fourth pass of its decade, nearly always shrinking, so the emission there is an upper value.

**The air**, as lake_anoxia reckons it: the lakes and seas at these fluxes, every other source cut 90% by biology, the
soils' 2.92 Tg a year per ppm the only sink; each pair is the domain synthesis's low- and high-cost ends, the room is
the Tg a year the waters could add before the air passes the held 1.87 ppm (a negative room is the removal that holds
it), and the dimming is the forcing beyond the held pair's 3.58 W/m² over the resort's 218.6 W/m²:

| Lakes and seas | Tg a year | Steady methane, ppm | Room under 1.87 ppm, Tg a year | Added dimming |
|---|---|---|---|---|
| Nutrient-poor, Earth's community at the methane share's 0.5 (designed–poor) | 0.57 (0.57–0.69) | 1.01–1.75 (1.01–1.79) | 2.51; 0.35 (2.39; 0.23) | −0.20 to −0.02% |
| Nutrient-poor at the share's 0.3; 0.8 | 0.001; 3.57 | 0.81–1.55; 2.04–2.78 | 3.08; 0.92 and −0.49; −2.65 | −0.27 to −0.06%; 0.04 to 0.17% |
| Mesotrophic at the share's 0.5 (0.3–0.8; designed–poor) | 7.8 (2.8–19.6; 4.9–33.5) | 3.49–4.23 (1.79–8.25; 2.50–13.0) | −4.7; −6.9 | 0.27 to 0.38% |
| Eutrophic | 34.0 (18.2–60.6; 19.1–106) | 12.5–13.2 (7.0–22.3; 7.3–37.9) | −31; −33 | 1.03 to 1.07% |
| New reservoir | 57.2 (30.5–103; 36.3–244) | 20.4–21.2 (11.3–36.8; 13.2–85.2) | −54; −56 | 1.37 to 1.40% |
| Under Shih and Bouffard's form, nutrient-poor (designed–poor) | 1.69 (0.72–12.8) | 1.39–2.13 (1.06–5.95) | 1.39; −0.77 (2.37; 0.21 designed) | −0.09 to 0.06% |
| Shih and Bouffard's, mesotrophic | 17.0 (7.3–74) | 6.6–7.4 (3.3–27) | −14; −16 | 0.64 to 0.70% |

With every anaerobic pathway at the high-erosion end (the land's earliest erosion, 8.3 Gt a year, the soil loss the
phosphorus particles imply, with highland soils) the nutrient-poor waters give 0.50 Tg a year (0.99–1.73 ppm, 0.42 Tg a
year of room at the high-cost end) and the mesotrophic 6.9 (3.2–3.9 ppm); with mare soils 0.43 and 5.1. At the
conservation study's 10 m a million years (0.71 Gt a year) the floors' iron cuts the lakes' emission by about 5%
nutrient-poor and 0.5–2.7% in the other states (on the rungs run with it), which leaves the air as above. Which erosion
the 301 K world has is open: the conservation study keeps 10–100 m a million years (0.71–7.1 Gt a year), the terrain
study's central rate is 50 m (3.55 Gt), and 8.3 Gt is sea chemistry's earliest, the young soils' loss. Splitting the
lakes' 100–352 m rungs at either end moves the steady methane by 0.01–0.23 ppm.

**Meromictic lakes 20–50 m deep.** The kinds study's lakes held by a chemocline at 20 m of 0.1 or 1 kg/m³ do not keep it
at 301 K: the nights mix them to the floor and wear the step to 1–4% of itself, their deep water never turns anoxic in
fifty years, and their methane is that of any nightly mixed lake of their depth, the floors' bubbles, in t of CH4 a year
per km² under Earth's community (designed–poor; all anaerobic pathways at the high-erosion end):

| Lake | Nutrient-poor | Mesotrophic | Eutrophic | New reservoir |
|---|---|---|---|---|
| 30 m | 0.337 (0.337–0.342; 0.175) | 1.96 (1.95–2.04; 1.50) | 6.9 (6.9–7.2; 6.3) | 13.6 (13.4–14.1; 12.8) |
| 50 m | 0.204 (0.204–0.207; 0.097) | 1.21 (1.21–1.27; 0.90) | 4.9 (4.9–5.1; 4.5) | 9.7 (9.7–10.2; 9.2) |

The same per km² holds every decade of the fifty years and for either step. A soda or spring-fed lake whose chemocline
does hold, by a stronger step, was not run; its anoxic deep water would add methane made below the step to these floors'
bubbles.

**The three targets beside lake_anoxia's two rules.** On lake_anoxia's 352 m lakes, regenerated on the final lake model,
the mean of the three bands at the central wind share, mol/m² a year. Those lakes carry lake_anoxia's mesotrophic
physics for every trophic state, where the ladder runs each state with its own, hence a new reservoir's 1.04 mol/m² a
year made in the resort at 300 K against 3.5 in the ladder's 352 m lake at 301 K:

| Lake | Tracer across the surface: Earth's (designed; poor) | lake_anoxia's pelagic flux: net rise; stock | Tracer's bubbles, own floors; lake_anoxia's, the Moon's floors | Made: tracer; lake_anoxia |
|---|---|---|---|---|
| Lock, new reservoir | 0.0057 (0.0002; 0.20) | 0.010; 0.012 | 0.100; 0.125 | 0.279; 0.275 |
| Resort at 300 K, new reservoir | 0.029 (0.0036; 0.68) | 0.059; 0.076 | 0.096; 0.122 | 1.04; 1.02 |
| Hot world, mesotrophic | 0.0008 (0.00002; 0.069) | 0.0049; 0.0052 | 0.010; 0.034 | 0.156; 0.156 |
| Hot world, eutrophic | 0.036 (0.011; 0.55) | 0.040; 0.048 | 0.043; 0.065 | 1.01; 0.99 |
| Hot world, new reservoir | 0.20 (0.088; 1.67) | 0.32; 0.62 | 0.093; 0.119 | 2.85; 2.77 |

The tracer makes what lake_anoxia makes, within 4%. With Earth's community its surface flux lies below both rules in
each lake with much anoxic water, 12–83% under the net rise, its methanotrophs oxidising more of what the overturns
bring up than either rule lets pass; in the resort's eutrophic lakes, whose anoxic water is slight, it lies between
them. With the designed community it lies lower still, and with the poor one well above both. The other lakes keep oxic
deep water and emit only bubbles, which the tracer's own floors give at 30–80% of lake_anoxia's reckoning from the
Moon's floors. lake_anoxia's 1,925 m lake, started oxygenated, stays oxic through fifty years nutrient-poor and
mesotrophic, emitting 0.0003–0.002 mol/m² a year (the hot world's mesotrophic lake begins to make methane in its fifth
decade); as a new reservoir it turns anoxic in its second to fourth decade in every target (2.6–2.8 mol/m² a year made
by the fifth) and its emission rises to 0.07–0.10, and eutrophic in the hot world from its second decade (1.6 made by
the fifth, 0.02 emitted).

## The methane-eating biology, sized

The ranges below are inputs, Earth's measured values carried to engineered communities as design variables; the cuts
they achieve, each set beside Earth's community at its central rates on the same lakes, are the results.

- **At a lasting oxycline (a).** Aerobic methanotrophs oxidise methane by Michaelis–Menten kinetics in methane and
  oxygen, their stock growing on what they eat and lost to grazing, lysis and settling. Each rate constant runs over
  Earth's measured range:
  - half-saturation in methane, 1.7–15 µM, central 5 (Mayr et al. 2020a; Zimmermann et al. 2021; Martinez-Cruz et al.
    2015);
  - in oxygen, 1–20 µM, central 5 (Ren et al. 1997; Lidstrom & Somers 1984; Blees et al. 2014 find it sub-micromolar);
  - growth, 0.16–2 per day at 20 °C, central 0.5: 0.16 per day is the maximum Zimmermann et al. (2021) model for
    Rotsee's overturn, and the oxycline doubling times of 1.8–3.2 days give 0.22–0.39 per day in situ (Oswald et al.
    2015, 2016);
  - carbon conversion efficiency, 0.13–0.6, central 0.3 (Zimmermann et al. 2021);
  - uptake's activation energy, 0.43–0.55 eV, central 0.49 (Thottathil et al. 2019);
  - loss rate, 0.03–0.3 per day, a design guess.

  In the 352 m lakes at 301 K only a new reservoir's oxycline carries enough methane for the kinetics to matter. It
  emits 0.33 mol/m² a year at Earth's central rates, 0.22–0.51 over the growth rate's range (2 to 0.16 a day), 0.28–0.42
  over the half-saturation in methane, 0.31–0.36 over that in oxygen, 0.28–0.38 over the loss rate and 0.29–0.35 over
  the yield; the activation energy moves it by 3%. The designed community, every constant at its favourable end, gives
  0.18 and the poor one 2.2. Eutrophic lakes move within 0.034–0.035, and the other states not at all, their emission
  being the floors' bubbles. The stocks' first-order constants come out at 0.005–0.02 per hour for Earth's community,
  inside Earth's measured 0.002–0.46 (Hudspeth et al. 2025; Thottathil et al. 2019), 0.57–0.81 for the designed one and
  0.0001 for the poor one.
- **At an overturn (b).** The stock the water keeps between overturns (10³–10⁶ cells per mL) and the community's
  response time (growth of 0.16–2 per day) set how much of what an overturn brings up is oxidised, within the oxygen the
  overturned lake holds. At 301 K the 30 m and 100 m lakes mix to their floors about ten times a year and the 352 m
  lakes three to six times, by band. Only a new reservoir's overturn brings up methane worth its oxygen: the 352 m lake
  has stored 1.8 mol/m² when one begins and then holds 3.7 g of oxygen per m², a thirtieth of what that methane needs,
  so the oxidation waits on oxygen from the air through the overturn's window; 0.23 mol/m² a year leaves within the
  windows and 2.2 is oxidised, 91% of what leaves the water there, within Rotsee's 91–98% (Zimmermann et al. 2021; Mayr
  et al. 2020b) and Earth's 54–98%. The stock and response time span its emission from 0.73 mol/m² a year with 10³ cells
  per mL growing 0.16 a day to 0.18 with 10⁶ growing 2 a day (Earth's central 10⁴ and 0.5 give 0.33); either a stock of
  10⁵ or growth of 2 a day cuts the slow, sparse end by 58–67%. A eutrophic 352 m lake stores 0.003 mol/m² at an
  overturn against 177 g of oxygen, and the grid moves it within 0.034–0.049; at 30 m and 100 m the grid moves the
  emission by at most 4%.
- **Anaerobic oxidation (c).** Each pathway is held to its electron acceptor's supply and to rates measured where it
  works:
  - iron(III) and manganese(IV) oxides settled from the rivers onto the floors (8 iron or 4 manganese a methane), at the
    land's erosion (sea chemistry) and the reactive share of Earth's river particles (Poulton & Raiswell 2002), with
    ferruginous Lake Matano's arithmetic (Crowe et al. 2011: iron burial five times the methane flux, iron carrying at
    most 66% of the oxidation, each atom recycled eighty times);
  - nitrate and nitrite from the oxic water, at Lake Constance's measured rates (Deutzmann et al. 2014);
  - sulphate, 1.4 µM in the lakes and none in the held seas, at its renewal.

  At the high-erosion end (the earliest erosion, 8.3 Gt a year, highland soils) all the pathways together cut the lakes'
  emission by 37–50% nutrient-poor, 20–39% mesotrophic, 7–11% eutrophic and 4–16% as new reservoirs, rung by rung; with
  mare soils and nitrate at its highest by 93–97%, 54–95%, 20–33% and 13–46%. Iron oxides carry nearly all of it, taking
  the floors' methane before it bubbles, which is why the shallow and nutrient-poor lakes gain most. At the conservation
  study's 10 m a million years (0.71 Gt a year, highland soils) the iron cuts only 5% of the nutrient-poor lakes'
  emission, 1.7–2.7% of the mesotrophic and under 1% of the eutrophic and the new reservoirs'; at the earliest erosion
  49% and with mare soils 97% of the nutrient-poor. Which erosion the 301 K world has is open (the air, above). Nitrate
  and nitrite matter only where the deep water is anoxic, 3–40% in the mesotrophic giant across their supply; sulphate
  cuts under 1%. Moon-wide the high-erosion end takes the nutrient-poor lakes from 0.147 to 0.079 Tg a year and the
  mesotrophic from 3.47 to 2.52; the conservation rate leaves them at about 0.14 and 3.4.
- **Whether nutrient-poor lakes at 301 K fit 1.87 ppm.** At the methane share's central 0.5 and under the lake model's
  Osborn form they do, with the seas, by a thin margin: 0.57 Tg a year leaves the air at 1.01–1.75 ppm with Earth's
  community and 1.05–1.79 with the poor one, with 0.35 Tg a year of room at the ledger's high-cost end (0.23 with the
  poor community), and the anaerobic pathways at the high-erosion end widen it to 0.42. The room would go to the 3,500 m
  nutrient-poor giant's likely 0.29 Tg a year once it turns anoxic; the seas' shelves at the lakes' rungs lean high and
  may understate it by up to 0.15. At the methane share's 0.8 the waters give 3.6 Tg a year and the air 2.04–2.78 ppm,
  needing 0.5–2.6 Tg a year of removal. Under Shih and Bouffard's form, whose weaker deep mixing turns the 1,925 m
  nutrient-poor giant anoxic, they fit with the designed community (1.06–1.80 ppm, 0.21 Tg a year of room) and with
  Earth's only at the ledger's low-cost end (1.39 ppm); Earth's at the high-cost end needs 0.77 Tg a year of removal,
  and the poor community's 12.8 Tg a year puts the air at 5.2–6.0 ppm. A designed community keeps the nutrient-poor
  waters within the held methane under either form at the central share.

## How it is set up

- **The columns.** The lake model is not run here. Each saved run gives its column about every day and the mixed layer's
  depth, the surface's temperature and oxygen and the dissipation below the mixed layer at every step. The tracer reads
  lake_anoxia's 352 m and 1,925 m lakes in the three targets, the lake kinds' 301 K ladder, giants and meromictic lakes,
  lake_oxygen's mixing-forms columns, and the seas' basins and their 3 m and 10 m pans. Each run's scenarios, physics
  and pairing of the bulk exchange follow the model version that made it, read from its recorded lake.py.
- **Transport.** Every step the recorded mixed layer mixes the tracer over its layers and the share of the next that the
  model entrained. Below it the model's own diffusivity applies: Osborn's rule (`lake.diffusivity`), the mixing
  efficiency times the step's dissipation over the squared buoyancy frequency under its cap, or Shih and Bouffard's
  three regimes where the run used them, or the lake number's mixing. The density comes from the day's temperature and
  solute and the debt's anaerobic products, mixed each step as the model mixed them, with the mixed layer at its saved
  temperature.
- **Methane made.** In water under 0.5 g/m³ of oxygen, the share of the demand respiration cannot meet, K/(K + O2),
  times the methane share of the anoxic mineralisation, over 64 g of oxygen a mole. The share is 0.5 central and 0.3–0.8
  bracket, as lake_anoxia reads it from Lake Matano (Crowe et al. 2011). Each floor's sediment makes the methane share
  of its own demand whatever the water above. In order, that methane goes to anaerobic oxidation in the mud, then to
  bubbles at lake_anoxia's porewater balance for the floor's depth, temperature and the methane in the water above, then
  into anoxic water; what diffuses out under oxic water is oxidised at the mud's surface. A bubble keeps to the surface
  what lake_anoxia's rising bubble keeps and leaves the rest in the layers it crosses.
- **Oxidation and oxygen.** Methanotrophs oxidise methane by the kinetics above, their stock mixed and diffused with the
  water. The model's oxygen already carries the reoxidation of its debt, which holds the methane. Each case therefore
  carries that methane as the model holds it, and the difference between the oxygen the model spends on it and the
  oxygen the methanotrophs use. The difference moves with the water and returns to the air at the surface's oxygen
  transfer. Methane made in anoxic water that escapes, bubbles away or is oxidised anaerobically gives back the oxygen
  the model spent on it, and a floor under oxic water keeps the oxygen the model spent on what its mud bubbles off or
  oxidises with iron, which leaves the water above slightly short of oxygen and its methane oxidation, limited by
  methane, all but unchanged; and where oxygen is short the methane takes its part of what the water holds beside the
  debt's other substances.
- **The air.** The mixed layer exchanges methane with the air at the model's own surface-renewal transfer velocity,
  rebuilt from the run's forcing, carried to methane by the Schmidt numbers' ratio (0.928), toward the held air's 1.87
  ppm.
- **Passes.** A record that repeats as a settled cycle is passed over until the inventories of the lanes whose
  methanotrophs oxidise move less than 0.5% (or 0.001 mol/m²) from one pass to the next; a lane without methanotrophs
  gathers methane in anoxic water until its store is full, when all it makes leaves, and is read so, beside its traced
  transient. The seas' central columns were traced first with every lane held to 0.5%. A giant, spun up by the kinds
  study (the 3,500 m nutrient-poor giant's deep water still losing oxygen when it stopped), is traced over its last
  saved decade as a settled cycle, up to ten passes. A lake
  held meromictic, which starts oxygenated with its chemocline, and lake_anoxia's 1,925 m lake, started oxygenated, run
  once through their fifty years by decades.
- **Scale-up.** The lakes by depth class from geography's rain-fed lakes on LOLA (lake_oxygen's lakes), each lake taking
  the rung nearest its deepest point (between geometric midpoints) and each rung its bands' fluxes weighted by its
  lakes' area in each band; a rung lacking a trophic state (no giants run as new reservoirs) takes the deepest rung run
  in it. The lakes whose deepest point lies between the 100 m and 352 m rungs are split at the lake kinds census's
  ventilated depth, the floor the nights reach on half of them, for each trophic state and band (168–239 m for
  nutrient-poor lakes): those within it take the 100 m rung's band-weighted flux and those beyond it the 352 m rung's.
  Splitting at either rung instead gives the bracket the census gives its self-ventilated lakes, 226,500–629,700 km²
  against its central 403,700 km². Under Shih and Bouffard's form the 352 m rung and the giants take their own emission
  plus the change the form makes per m² in the same trophic state and community where it was run (the 352 m lake and the
  1,925 m giant), and the rungs to 100 m, within the ventilated depth, stay as under Osborn's. The seas by basin, each
  basin's water by local depth from its own shelf shares in seas.json (the water it stands for beyond its own at the
  shares the seas' whole leaves once every basin's own shelves are counted, so that the seas' total keeps seas.json's),
  interpolated in the logarithm of depth: to 5 m at the 3 m pan, 5–30 m at the 10 m pan, 30–200 m at the lakes' ladder
  rungs nearest that depth (between the geometric midpoints) in the basin's band, and deeper at its typical deep column,
  with the coldest columns as the bound. The pans were not run at the methane share's ends; there each takes its
  Earth's-community flux times the ratio the lakes' rung of its depth shows. Under Shih and Bouffard's form the two
  basins run under both forms, Procellarum and Ingenii, take their own traced rows for their deep water, and every other
  basin's typical column takes its own emission plus the change per m² theirs show in the same trophic state and
  community (their mean weighted by what they stand for, each alone for the bracket), its shelves to 200 m as under
  Osborn's. The air through lake_anoxia's own `air()`.

## Assumptions and limits

- **Offline on saved columns.** The tracer takes the model's temperature, solute and oxygen as the model saved them,
  each day, and the mixed layer's oxygen each step. Its own oxygen use enters only through the oxygen difference, so the
  model's other demands keep the oxygen the model gave them.
- **Earth's rates for engineered communities.** The kinetics, the methane share, the bubble sizes and the anaerobic
  rates are Earth's, carried to the Moon as design ranges. No lunar lake or sea has been measured.
- **The lakes the ladder did not run.** No giant was run as a new reservoir, nor eutrophic at 1,000 m or 3,500 m; those
  rungs take the deepest lake run in that state (the 352 m lake for new reservoirs, the 1,925 m giant for the
  eutrophic). The 1,000 m eutrophic borrowing is probably an upper value (about 2.9 Tg a year at the 352 m lake's flux);
  the new reservoirs' giants are a bound whose direction is open. The ladder runs every rung in all three bands; the
  1,000 m and 3,500 m giants and the eutrophic 1,925 m giant run in the 0–15° band only and stand for theirs, at the
  high end of the mesotrophic giants' spread (about 9% on the 3,500 m mesotrophic giant).
- **How far the passes settled.** Every ladder rung and lake_anoxia's 352 m lakes settled within 0.5% in two to four
  passes of their decade, and eleven of the fifteen giant cases (the wind share's brackets among them) in two to eight.
  The other four (the mesotrophic 1,925 m in two bands and at 3,500 m, the eutrophic 1,925 m) ended their ten passes, a
  century of their settled decade, with the stores still shrinking 0.2–0.4% a pass and their emission 0.2–0.5% a pass,
  slight upper values. The giants' own oxygen settles less evenly (kinds.giants): the 3,500 m nutrient-poor giant is
  still losing it and turns anoxic within 12.6 years of its record, its methane a lower value by about 0.29 Tg a year;
  the 1,000 m giants and the 1,925 m nutrient-poor giants gain it slowly, their methane as traced or slightly high. The
  seas' deep columns ran four passes, and 13 of the 20 still moved by 0.5–12% in the last (the seas, above); their
  passes' own emission was not recorded, as they were traced before the tracer recorded it. A lane without methanotrophs
  never settles: its store fills over centuries, and all it makes would then leave, so its traced figure is a transient
  below that limit.
- **The ventilated depth.** The lake kinds census's ventilated depth, interpolated between the 100 m and 352 m rungs,
  splits those lakes; each side takes its rung's band-weighted flux, and the split at either rung is the bracket. The
  census's depths are Osborn's; under Shih and Bouffard's form the same split is kept.
- **Shih and Bouffard's form, carried by addition.** The form was run on the 352 m lake and the 1,925 m giant at 301 K
  and on two sea basins. Their change per m² is added to the rungs and basins not run, which stays bounded where
  Osborn's flux is near nothing, as in a ventilated column, where a ratio would not; the rungs within the ventilated
  depth and the seas' pans stay as under Osborn's. Procellarum's column under the form is the decade after its one
  overturn in a century, its methane a lower value.
- **The seas' shelves.** The shelves from 30 m to 200 m take the lakes' rungs of the same local depth, whose fluxes
  carry each lake's own shallow margins' bubbles; the shelf water at that depth has its floor at that depth, so the
  figure leans high, by up to 0.15 Tg a year nutrient-poor, with the deep columns' flux there as the lower end.
- **Meromixis.** The kinds study's chemoclines of 0.1 and 1 kg/m³ at 20 m do not hold at 301 K, so the meromictic lakes'
  flux is a nightly mixed lake's; a lake held by a stronger step was not run.
- **The air.** lake_anoxia's reckoning: the domain synthesis's ledger with its lakes replaced by these waters, every
  other source cut 90% by biology, the soils the only sink in the bare column, the steady methane the sources over that
  sink.
- **What a column cannot hold.** Sideways exchange, sills, the slopes' mixing and the cold end's sinking water are
  beyond a one-dimensional column, as the lake and sea studies set out.
- **Open in the literature.** A culture's maximum growth rate for the designed end (sources give "a few hours"), the
  uptake's activation energy (an abstract's summary), freshwater iron-dependent rates (Lake Matano's and Kinneret's
  arithmetic stand in), and the methanotrophs' loss rate, a design guess.
- **The targets as run.** The columns ran, and their transfer velocities are rebuilt, on the climate targets the lakes
  ran on, domain/ecology's at e9c12eb (schema /2), which lake_oxygen keeps byte for byte as
  climate_targets_as_run.json until the lakes rerun on the resort's 301 K design case. The merged targets would cool
  the lakes' air by 0.026 K at most and 0.015 K typically and leave the seas' as it is.

## The checks

[test_run.py](test_run.py) checks:

- the methane budget closing by every road;
- the oxygen difference's own budget;
- steady diffusion against first-order oxidation matching its analytic profile and escape;
- the bubbles against lake_anoxia's `ebullition_share` and Lake Soppen's measured share;
- the floors' and water's demands giving back the model's;
- the mixed layers recovered from the recorded depths;
- the methanotrophs growing on methane and falling back to their seed;
- the acceptors' supply;
- on lake_anoxia's saved columns, the methane made where its rule coincides, its overturns, and the tracer's transport
  against the model's own age and debt;
- the product's binding to its code and inputs.

The rebuilt transfer velocity's friction velocity and the wind's dissipation are checked against the values the runs
saved. Over the whole of lake_anoxia's regenerated records (the product's `reproduction`), the tracer makes the methane
lake_anoxia makes within 4% in all nine 352 m lakes and counts the same overturns; the debt's clone keeps the model's
own debt within 7.1% where its store passes 0.01 mol/m² (11% below that), and the water's age within 0.08–0.78 days on
average.

## Running it

```sh
python -m research.studies.lake_methane.run --checks    # lake_anoxia's columns: the reproduction
python -m research.studies.lake_methane.run --seas      # the seas' columns
python -m research.studies.lake_methane.run --targets   # lake_anoxia's regenerated columns, the three targets
python -m research.studies.lake_methane.run --work      # one of up to four processes tracing into the cache
python -m research.studies.lake_methane.run             # everything, then results/lake_methane.json
python -m research.studies.lake_methane.run --cache-only   # the product from the cache, tracing no case
python -m pytest research/studies/lake_methane
```

Each traced case is cached in research/runs/lake_methane/cases (ignored, on the research drive) under a hash of the code
that traces it, the saved files' contents (their sha256, so a producer's renaming leaves it standing), the column's
diffusivity form, its lanes and its options, so a stopped run resumes at the case it stopped in. The transfer series
rebuilt from each run's forcing are cached in research/runs/lake_methane/transfer under a hash of the readers' code,
what the series read from the setting's inputs (lake_oxygen's `Setting.values_fingerprint()`) and the saved files'
contents; with `--cache-only` the run traces no case and builds no series, and stops at the first one not cached. The
--work processes claim each case by a file in research/runs/lake_methane/claims before tracing it, so that up to four
can share the work; the folder is cleared before they are started again.
