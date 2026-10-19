# The waters' oxygen in the three target climates

17 October 2026; rerun on 18 October on the lake model as corrected after that day's independent review and its
re-check (the findings, the re-check and their scripts are on the research drive, `review-lakes-2026-10-18`), and with
the mixing below the mixed layer driven by the wind, the treatment the author chose that day. Ecology's
[question 10](../../../biosphere/ecology/questions.md) plans a one-dimensional lake model through the lunar month for
the stratification and oxygen of the deep crater lakes, and the author asked for it in each of the three targets of
[what each target climate takes](../climate_targets/README.md): the hot world near 307 K, the resort at 300 K and the
lock at 295 K. This study gives the oxygen the held fresh seas and the lakes hold at saturation in each, and runs the
[lake model](../../../biosphere/ecology/lake.py) through ten years of each target's lunar months, from the GCM's own
output, for a crater lake 352 m deep in each of three latitude bands and a productive shallow lake 3 m deep beside it,
under oxygen demands from oligotrophic to eutrophic water. The corrections of 18 October pair the surface exchange as
PlaSim pairs it, at the lowest level's height; carry the anaerobic products of the oxygen debt in the density; draw the
euphotic water's night respiration from a labile pool; drive the mixing below the mixed layer by the share of the
wind's energy that Earth's lakes pass to their stratified water, with the share's measured range as the bracket, the
earlier fixed dissipations and a Lake-number scheme as sensitivity rows, under a cap on the diffusivity that holds
nothing back; fix the demand by supply, the warm kinetics from the temperature it was measured at a transient beside
it; settle each deep lake on its water's warmth, oxygen, age and solute in coarse steps and then at the analysis step;
and read the Lake number at the lunar month's thermocline as well as the day's warm layer. After the review's last
check of that day, [mixing_forms.py](mixing_forms.py) runs the deep lakes with the form that turns the dissipation into
mixing taken from Shih et al. (2005) in place of Osborn's rule, and the 3 m pans at their share's bracket, writing
[results/mixing_forms.json](results/mixing_forms.json) (schema `terluna.research.lake-oxygen-mixing-forms/1`).
[run.py](run.py) writes [results/lake_oxygen.json](results/lake_oxygen.json) (schema `terluna.research.lake-oxygen/2`)
from [lakes.py](lakes.py)'s choice of lakes ([results/lakes.json](results/lakes.json)) and [extract.py](extract.py)'s
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
standard values. At the earlier sea medians of 34.0, 27.7 and 23.0 °C it gives 7.13, 7.91 and 8.63 mg/L, within
0.02 mg/L of the 7.14, 7.93 and 8.63 of the integration history review of 17 October (its finding F04, from the same
fit); the held seas' salts take 0.006 of that.

**The wind's energy below the night's mixing turns the deep crater lakes over every second to fourth lunar month by
Osborn's rule, about once a year by the energetic regime's.** Each
night's cooling and wind mix the top of the 352 m lakes, to 9–14 m even in the calmest tenth of lunar cycles. Below the
mixed layer the share of the wind's energy that Earth's lakes pass to their stratified water dissipates a mean 4.1 ×
10⁻¹⁰ to 1.3 × 10⁻⁹ W/kg, and it keeps the deep water close to neutral: it carries the warmth that the thermocline and
the floor pass down through the whole deep water, so the deep water stays near the warmth that the ordinary cold nights
of many lunar months reach. With Osborn's rule turning that dissipation into mixing, in the lock and the resort the
whole column turns over in 55–74 of the 120 lunar cycles of ten years, in the hot world in 29–30, and its water is then
a mean 0.1–0.3 years from its last contact with the air in the lock and the resort and 0.4–1.4 years in the hot world,
3.8 years at the most. With Shih and Bouffard's energetic regime (below), the 0–15° lake turns over 11 times in the
lock, 8 in the resort and 12 in the hot world, its deep water a mean 1.7–2.8 years old. Earth's deep tropical lakes
overturn only in the coldest of their seasons (Fukushima et al. 2017); in these lakes the near-neutral deep water lets
the coldest nights of many lunar months through.

| Overturns in ten years, of 120 lunar cycles; tenth-percentile and median of each cycle's deepest mixing; deep water's mean age | Mandel'shtam, 4° N | Gagarin, 18° S | Rabbi Levi, 35° S |
|---|---|---|---|
| The lock, 0.3% (central) | 66; 13 and 352 m; 0.2 years | 74; 14 and 352 m; 0.1 years | 62; 11 and 352 m; 0.3 years |
| The lock, 0.2% | 53; 12 and 21 m; 0.2 years | 61; 13 and 352 m; 0.2 years | 47; 11 and 17 m; 0.3 years |
| The lock, 0.4% | 80; 13 and 352 m; 0.1 years | 84; 138 and 352 m; 0.1 years | 73; 12 and 352 m; 0.1 years |
| The resort, 0.3% | 58; 12 and 29 m; 0.2 years | 68; 13 and 352 m; 0.1 years | 55; 11 and 17 m; 0.3 years |
| The resort, 0.2% | 40; 12 and 18 m; 0.3 years | 45; 12 and 146 m; 0.2 years | 42; 11 and 16 m; 0.5 years |
| The resort, 0.4% | 67; 12 and 352 m; 0.2 years | 75; 14 and 352 m; 0.1 years | 62; 11 and 352 m; 0.3 years |
| The hot world, 0.3% | 29; 11 and 16 m; 1.4 years | 29; 11 and 15 m; 0.4 years | 30; 9 and 14 m; 0.7 years |
| The hot world, 0.2% | 22; 11 and 15 m; 2.0 years | 18; 10 and 14 m; 1.0 years | 22; 9 and 14 m; 1.0 years |
| The hot world, 0.4% | 39; 11 and 16 m; 0.8 years | 38; 11 and 137 m; 0.2 years | 38; 9 and 14 m; 0.4 years |

The share's bracket moves the overturns by 10–34% in the lock and the resort and by 24–38% in the hot world, most at
its low end. A median at the floor means more than half the cycles overturn.

**The stirring one column cannot hold.** The Lake number (Imberger & Patterson 1990) weighs the stratification's
resistance against the wind's stress, at the depth of the thermocline. With it at the lunar month's thermocline, the
strongest stratification below 10 m, the 352 m lakes' Lake number is a median 0.30–0.93, under 1 on 51–90% of the days;
with it at the strongest stratification below 2 m, which the review found to be the day's shallow warm layer on two days
in three under the first round's mixing, a median 0.10–0.37, under 1 on 84–100%. On those days the wind tilts their
thermocline to the surface and stirs their deep water by basin-scale motions, upwelling at the shores and internal
waves, which a one-dimensional column cannot represent. At lunar gravity the same warmth resists the wind six times more
weakly than on Earth, where the review found the same water and wind give a Lake number about eight times higher. The
column stands for that stirring by the share of the wind's energy measured in Earth's lakes, with its bracket in every
table of the deep water below; a basin-scale model of one crater lake would test it.

**Their deep water keeps its oxygen in nutrient-poor water; in richer water it turns on the mixing's form and the
warmth.** In the deepest tenth of each lake's volume, by Osborn's rule:

| Deepest tenth's oxygen, mean (at the share's bracket); share of time under 0.5 mg/L | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic | 6.8–7.4 mg/L (6.3–7.5); 0% | 5.9–6.5 mg/L (5.3–6.6); 0% | 4.2–4.9 mg/L (3.1–5.4); 0% |
| Mesotrophic | 5.0–6.1 mg/L (4.2–6.3); 0% | 4.0–5.0 mg/L (3.1–5.4); 0% | 1.8–2.8 mg/L (0.7–3.7); 0–29% |
| Eutrophic | 2.9–4.2 mg/L (2.0–4.5); 0% | 1.9–3.2 mg/L (0.9–3.6); 0–2% | 0.39–0.89 mg/L (0.0–1.7); 51–63% |
| New reservoir | 0.8–2.7 mg/L (0.2–3.1); 0–37% | 0.2–1.5 mg/L (0.0–2.1); 10–79% | 0.00–0.20 mg/L (0.00–0.40); 86–100% |

The ranges run over the three bands, the brackets over the share's 0.2% and 0.4%. The oxygen the floor and the water
below the euphotic zone take, 0.32–1.31 g per m² a day across the states, empties the deep water between overturns,
and each overturn refills it. Oligotrophic and mesotrophic deep water stays oxic in the lock and the resort through the
record; eutrophic water keeps 1.9–4.2 mg/L there and runs anoxic for half the time in the hot world, whose lakes turn
over half as often; new reservoirs run low in every target. By the energetic regime the 0–15° lake's deepest tenth keeps
4.8, 3.4 and 1.6 mg/L oligotrophic in the lock, the resort and the hot world (the hot world's anoxic 35% of the time),
2.9, 1.2 and 0.2 mesotrophic (anoxic 0%, 45% and 86%), and 1.0 mg/L eutrophic in the lock (45%), its eutrophic water
and new reservoirs elsewhere anoxic most or all of the time. The deep-oxygen conclusions run between the two forms:
oligotrophic deep water in the lock and the resort stays oxic under both, and everything richer, or warmer, rests on the
form.

**An overturn brings the deep water up.** When the column turns over, the whole lake takes the mix of its oxygenated
top and its deep water, and the reduced substances the deep water has gathered take up oxygen as they meet it. Just
after an overturn the whole lake, its surface with it, holds:

| Whole lake just after an overturn, median over the overturns (at the share's bracket) | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic | 6.9–7.5 mg/L (6.4–7.6) | 6.0–6.6 mg/L (5.5–6.8) | 4.9–5.5 mg/L (4.0–5.7) |
| Mesotrophic | 5.2–6.3 mg/L (4.4–6.5) | 4.3–5.3 mg/L (3.4–5.6) | 2.9–4.0 mg/L (1.4–4.3) |
| Eutrophic | 3.2–4.5 mg/L (2.2–4.7) | 2.3–3.6 mg/L (1.2–3.9) | 0.7–2.2 mg/L (0.1–2.5) |
| New reservoir | 1.2–3.1 mg/L (0.3–3.4) | 0.4–2.0 mg/L (0.1–2.5) | 0.02–0.41 mg/L (0.01–0.78) |

At the central share an overturn takes the median surface below 2 mg/L in two of the hot world's three eutrophic lakes
and in six of the nine new reservoirs, all three of the hot world's among them; there the surface climbs back above
2 mg/L within 0.4–1.1 days (the median after each), as the day's warm layer forms and takes oxygen from the air, while
in the hot world's new reservoirs the column below stays anoxic. Earth's deep tropical lakes show the same: fish die
when strong winds bring their oxygen-poor deep water to the surface (Fukushima et al. 2017). Oligotrophic and
mesotrophic overturns leave every lake's surface above 2 mg/L.

**The very deep lakes turn over about once a year, and their richer depths go anoxic.** Half the lake area lies in lakes
deeper than 1,925 m. Mandel'shtam's lake taken to that depth and oxygenated at the start, as if just overturned, turns
over in 33–39 of its first 120 lunar cycles while its deep water settles from the 352 m lake's warmth. Passed over the
record until it settles, 130–220 years in coarse steps and then at the analysis step, it turns over in 10–11 of 120
cycles at the central share (7–9 at 0.2%, 12–16 at 0.4%), and its deepest tenth keeps a mean 4.0 mg/L of oxygen in
oligotrophic water in the lock, 2.0 in the resort and 0.8 in the hot world (2.2, 1.2 and 0.03 at 0.2%; 5.0, 2.9 and 1.3
at 0.4%); mesotrophic water is anoxic at depth in the resort and the hot world and keeps 0.7 mg/L in the lock, and
eutrophic water and new reservoirs are anoxic at depth in every target. An overturn mixes the whole lake, but its
oxygenated top is a small share of its volume, so the mix holds 1.1–4.1 mg/L in oligotrophic water and the demand of the
floor and the deep water empties it again before the next; the deep water is a mean 1.1–2.8 years and up to 7.2 years
from its last contact with the air. In the lock at the central share and at 0.4% the giant had not settled at the
analysis step in a hundred years, its mesotrophic deepest tenth moving 0.05 mg/L between the analysed decade's first two
years and its last two. By the energetic regime the same giant at the resort's 301 K turns over twice in fifty years
and goes anoxic at depth even oligotrophic, its deep water a mean 29 years old (mixing_forms).

**What rides on the stirring.** At the 0–15° lake, overturns in ten years, the tenth-percentile and median of each
cycle's deepest mixing, and the deepest tenth's mean oxygen in oligotrophic and mesotrophic water:

| Below the mixed layer | The lock | The resort | The hot world |
|---|---|---|---|
| The wind's 0.3% (central), cap 0.1 m²/s | 66; 13 and 352 m; 7.0 and 5.4 mg/L | 58; 12 and 29 m; 6.1 and 4.4 mg/L | 29; 11 and 16 m; 4.2 and 1.8 mg/L |
| The wind's 0.2% | 53; 12 and 21 m; 6.6 and 4.9 mg/L | 40; 12 and 18 m; 5.6 and 3.6 mg/L | 22; 11 and 15 m; 3.1 and 0.7 mg/L |
| The wind's 0.4% | 80; 13 and 352 m; 7.2 and 5.8 mg/L | 67; 12 and 352 m; 6.3 and 4.8 mg/L | 39; 11 and 16 m; 4.8 and 2.7 mg/L |
| The wind's 0.3%, a cap of 1 m²/s | 66; 13 and 352 m; 7.0 and 5.4 mg/L | 59; 12 and 29 m; 6.1 and 4.4 mg/L | 29; 11 and 16 m; 4.2 and 1.8 mg/L |
| The wind's 0.3%, the first round's cap of 10⁻⁴ m²/s | 6; 13 and 39 m; 4.1 and 2.2 mg/L | 5; 12 and 31 m; 2.7 and 1.0 mg/L | 14; 11 and 16 m; 2.3 and 1.0 mg/L |
| A fixed 10⁻⁹ W/kg, cap 10⁻⁴ m²/s | 6; 19 and 54 m; 4.6 and 3.0 mg/L | 5; 17 and 50 m; 3.5 and 1.8 mg/L | 16; 11 and 46 m; 3.2 and 1.9 mg/L |
| A fixed 10⁻¹⁰ W/kg, cap 10⁻⁴ m²/s | 5; 12 and 16 m; 1.5 and 0.2 mg/L | 5; 12 and 15 m; 0.7 and 0.0 mg/L | 11; 11 and 15 m; 0.5 and 0.0 mg/L |
| DYRESM's Lake-number diffusivity (Yeates & Imberger 2003) | 0; 26 and 35 m; 0.0 and 0.0 mg/L | 0; 24 and 32 m; 0.0 and 0.0 mg/L | 0; 21 and 28 m; 0.0 and 0.0 mg/L |

The cap carries most of the change from the first round. Under that round's cap of 10⁻⁴ m²/s the same wind turns the
lakes over 5–14 times in ten years, where the raised cap gives 29–66, and leaves their oligotrophic deep water
2.3–4.1 mg/L: the old cap held back 45–54% of the interfaces below the mixed layer, near-neutral water whose mixing it
slowed, and so let the deep water's warmth and solutes settle into a weak stratification of their own. A cap ten times
the raised one leaves the lakes as they are, to one overturn in ten years and 0.003 mg/L. The first round's fixed
dissipations under their cap give 5–16 overturns. DYRESM's diffusivity, the one published one-dimensional form found
that scales the deep mixing with the Lake number, sets its mixing by the stratification's own profile, N²/N²_max, so it
mixes the thermocline hard and the near-neutral deep water hardly at all: it turns no lake over in ten years and leaves
the deep water anoxic in every target and state. It stirs these lakes' depths less than Earth's measured share of the
wind's energy does, so it stands as a row of its own beside the share's bracket.

**The form of the mixing below the mixed layer carries as much as the share.** The model turns the dissipation into
mixing by Osborn's rule, K = 0.2 ε/N², at every buoyancy Reynolds number Re_b = ε/(νN²). In the near-neutral deep
water the wind's share keeps, Re_b reaches 10³–10⁵, where Shih et al. (2005) find the mixing efficiency falling:
K = 2ν Re_b^½ above Re_b = 100, Osborn's 0.2 between 7 and 100 and molecular diffusion alone below 7. Bouffard &
Boegman (2013) carry Shih's three regimes to lakes with a fourth for strong thermoclines, left out here, which would add
mixing near the molecular value there. In the energetic regime the diffusivity is Osborn's rule divided by 0.1 Re_b^½,
3–30 times less at 10³–10⁵. [mixing_forms.py](mixing_forms.py) runs the 0–15° lake under those regimes with
everything else as here (ν 10⁻⁶ m²/s, the same cap, unstable rule, spin-up and forcing):

| At the 0–15° lake: overturns in ten years; deep water's mean age; deepest tenth oligotrophic, mesotrophic, eutrophic, new reservoir (share of time under 0.5 mg/L) | Osborn's rule | Shih and Bouffard's regimes |
|---|---|---|
| The lock | 66; 0.2 years; 7.0, 5.4, 3.7, 1.9 mg/L | 11; 1.7 years; 4.8, 2.9, 1.0 (45%), 0.07 mg/L (94%) |
| The resort at 300 K | 58; 0.2 years; 6.1, 4.4, 2.7, 0.7 mg/L (35%) | 8; 1.9 years; 3.4, 1.2 (45%), 0.18 (83%), 0.00 mg/L (100%) |
| The hot world | 29; 1.4 years; 4.2, 1.8 (29%), 0.54 (61%), 0.00 mg/L (100%) | 12; 2.8 years; 1.6 (35%), 0.17 (86%), 0.00, 0.00 mg/L (100%) |
| The resort at 301 K, each state with its own physics | 44–55; 0.2–0.3 years; 6.0, 4.2, 1.0 (7%), 0.00 mg/L (100%) | 5–14; 0.8–2.1 years; 3.5, 1.0 (48%), 0.00, 0.00 mg/L (100%) |
| The 1,925 m giant at 301 K, oligotrophic; mesotrophic: overturns in fifty years; deepest tenth | 40, 2.7 mg/L; 50, anoxic | 2, anoxic; 15, anoxic |

By the energetic regime the lakes turn over a sixth to a seventh as often in the cooler targets and less than half as
often in the hot world. Oligotrophic deep water keeps 3.4–4.8 mg/L in the cooler targets and the lock's mesotrophic
water 2.9 mg/L, while the resort's mesotrophic and richer water and the hot world's run low or anoxic. Holding the
dissipation at most at the share spread over the water below 20 m, deeper than the ordinary nights mix the 352 m lakes,
lowers the lock's mean dissipation by 0.8% and leaves its overturns and deep oxygen as they are: the inverse scaling
with the water left below the mixed layer, which raises ε as an overturn nears, does not make the overturns. The
lake kinds' 301 K ladder and census, the seas' and the methane runs take Osborn's rule, so their figures lean towards
more mixing.

**What respiration's products do.** The same lake, with what its water carries in the density:

| At the 0–15° lake: overturns; tenth-percentile and median mixing; deepest tenth, oligotrophic and mesotrophic | The lock | The hot world |
|---|---|---|
| Respiration's CO2 and the debt's anaerobic CO2 and methane (central) | 66; 13 and 352 m; 6.98 and 5.45 mg/L | 29; 11 and 16 m; 4.18 and 1.79 mg/L |
| No solutes | 66; 13 and 352 m; 7.00 and 5.48 mg/L | 31; 11 and 16 m; 4.25 and 1.90 mg/L |
| The CO2 dissolving calcite into bicarbonate | 53; 12 and 22 m; 6.77 and 5.05 mg/L | 21; 11 and 15 m; 3.38 and 0.84 mg/L |
| The central solutes and iron reduction's iron and bicarbonate | 66; 13 and 352 m; 6.98 and 5.45 mg/L | 29; 11 and 16 m; 4.17 and 1.79 mg/L |
| The physics of eutrophic water | 75; 9 and 352 m; 7.39 and 5.94 mg/L | 22; 5 and 7 m; 3.89 and 0.90 mg/L |

With the deep water renewed every few months, respiration's products hold it apart little: their own density at the
bottom is 0.001–0.002 kg/m³, and without them the lakes overturn as often and keep as much oxygen. CO2 dissolving
sediment calcite into bicarbonate, ten times the CO2's density, takes the lock's overturns from 66 to 53 and the hot
world's from 29 to 21, and the hot world's mesotrophic deep water from 1.79 to 0.84 mg/L; iron reduction's dissolved
iron and bicarbonate, up to Lake Matano's 120 µM of iron and net of the CO2 the reduction consumes, change nothing
where the water stays oxic. The physics of eutrophic water, its sunlight absorbed within the top 2.7 m, lets the
lock's lake overturn more and mixes the hot world's lake less. The resort built from the hot world's weather overturns
46 times in ten years where the lock's weather gives 58, and keeps 5.68 and 3.73 mg/L of oligotrophic and mesotrophic
deep water against 6.09 and 4.42.

**Productive shallows run low in the night, most of all in the warm targets.** In the 3 m lake, each night's lowest
surface oxygen (the median night, the worst tenth of nights in brackets):

| Lowest at night: median (worst tenth) | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic | 7.0–7.2 mg/L (6.9–7.0) | 6.6 mg/L (6.4–6.5) | 6.0 mg/L (5.8–5.9) |
| Mesotrophic | 6.3–6.7 mg/L (6.0–6.3) | 5.8–6.1 mg/L (5.5–5.7) | 5.2–5.4 mg/L (5.0–5.1) |
| Eutrophic | 3.8–4.8 mg/L (3.2–4.2) | 3.3–4.1 mg/L (2.8–3.4) | 2.8–3.3 mg/L (2.2–2.5) |
| New reservoir | 3.4–4.5 mg/L (2.7–3.8) | 2.9–3.8 mg/L (2.4–3.1) | 2.4–3.0 mg/L (1.7–2.2) |

The euphotic water respires four fifths of the day's production from a labile pool with an e-folding of 4 days, so
through the 354-hour night its respiration runs down. The wind's share of its energy, spread over the little water below
a pan's warm layer, mixes the pan through the lunar day, so its floor stays oxic through the day: an extension of a
share measured below deep lakes' thermoclines into a 3 m pan, where it gives 1.1–2.1 × 10⁻⁸ W/kg, thirty times the deep
lakes'. Its bracket (mixing_forms) runs from the share at 0.15%, which keeps the same oxic floors and eutrophic nights
at 2.7–4.8 mg/L, to the first round's fixed dissipation of 10⁻¹⁰ W/kg, under which the eutrophic pans' floors are anoxic
27–45% of the time, 176–242 hours of the day, and their median night falls to 2.7–4.0 mg/L in the lock, 2.2–3.1 in the
resort and 1.7–2.3 in the hot world, below 2 mg/L on up to 12%, 38% and 70% of nights. With the central share, eutrophic
shallows fall below 2 mg/L on at most 1% of the hot world's nights, new reservoirs on up to 3% of the resort's and 4–19%
of the hot world's, and no shallows fall below 0.5 mg/L; Gagarin's, with the weakest winds, are the lowest. The lowest
point comes a sixth to a quarter of the way through the night in productive water, and by dawn the surface has climbed
back to 7.1–7.5 mg/L in the eutrophic lock, 6.4–6.7 in the resort and 5.3–5.7 in the hot world, as the pool runs down,
the cooling water holds more oxygen and takes it from the air. Oligotrophic shallows keep 90–94% of saturation at their
lowest, mesotrophic ones 78–85% and eutrophic ones 42–60%.

The brackets carry the nights. With the pool's e-folding at 2 days the eutrophic shallows' median night keeps 3.7–5.6
mg/L, at 7 days 2.3–4.1; Cole & Caraco's transfer, which takes less from the air in light wind, puts it at 2.6–3.7. The
warm columns let the floor's demand rise from the 6 °C its rates were measured at by ecosystem respiration's 0.62 eV,
5–13 times at the shallows' mean 25–37 °C, with the pool and production held. No settling supply keeps a floor
respiring at that rate, so they show the weeks to months after a warming while the floor's stock of settled matter
lasts: in such a spell even oligotrophic shallows fall to 1.4–2.5 mg/L in the hot world's median night and 4.2–4.6 in
the resort's. The bounds above are the demand the supply fixes.

**What it means for the targets.** Fresh water gives the hot world's seas a fifth more oxygen than seawater at their
warmth would hold, 10% more than Earth's tropical seawater holds, and the warmth still leaves them 16% below the lock's.
The deep crater lakes' oxygen runs between two forms of the same wind-driven mixing. By Osborn's rule they turn over
every second lunar month in the lock and the resort and every fourth in the hot world, their oligotrophic and
mesotrophic deep water stays oxic in the lock and the resort, eutrophic water keeps 1.9–4.2 mg/L there, and in the hot
world mesotrophic deep water runs at 1.8–2.8 mg/L and eutrophic water anoxic half the time. By Shih and Bouffard's
energetic regime they turn over about once a year, oligotrophic deep water keeps 3.4–4.8 mg/L in the cooler targets and
1.6 in the hot world, mesotrophic water runs at 1.2–2.9 mg/L in the cooler targets and anoxic most of the time in the
hot world, and richer water goes anoxic 45–100% of the time. The very deep lakes, half the lake area, turn over
about once a year by Osborn's rule with their richer depths anoxic and oligotrophic water keeping 0.8–4.0 mg/L; by the
energetic regime the oligotrophic giant turns over once in decades and every state's depths go anoxic. The warmth tells
at every depth: the hot world's lakes keep less than the resort's, and the resort's less than the lock's. The share is
Earth's, its drag and its reach in lakes the wind tilts to the surface open, and the first round's cap, which held back
half the deep interfaces, takes the overturns to a tenth in the lock and the resort and to half in the hot world; a
basin-scale model of one crater lake would test the share and the form. Oxygen-breathing life in the deep crater lakes
needs nutrient-poor water, which keeps the 352 m lakes' depths oxic in the cooler targets under both forms; the giants'
depths need oligotrophic water and, below that, renewal by design (oxygen piped in, as some of Earth's eutrophic lakes
and reservoirs are kept, or driven mixing) or deep communities built for little oxygen. Productive shallows stay above 2
mg/L on almost every night by the wind's share, and run down to 1.7–4.0 mg/L in the eutrophic median night, below 2 mg/L
on up to 70% of nights in the hot world, at the bracket's fixed dissipation.

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
three cells are land in PlaSim's mask): PlaSim's tas, which is the lowest level's temperature brought down its dry
adiabat to the surface pressure, that level's potential temperature at the surface (fluxmod.f90's dtsa, 1.42–1.48 K
above the level's own temperature at every cell); the lowest level's humidity and wind, about 0.86 km up; the surface
pressure; and the surface's net sunlight and infrared, its albedo and temperature. These are means over each window,
all but the humidity: PlaSim writes the humidity of the step that opens the window, held through it (plasim.f90
computes it only at output steps and outmod.f90 accumulates the held value), so each window takes the mean of its own
record and the next, its opening and closing steps. Downward sunlight is the net over one less the albedo, downward
infrared the net over the surface's emissivity (1 on land, 0.98 on water, as PlaSim's radiation takes them) plus the
surface's own emission. Within each mean the sunlight follows the Sun's height at the lake on the synodic clock, scaled
to the mean; the subsolar longitude fitted to every mean's Sun reproduces the GCM's mean cosine of the zenith angle to
0.0006–0.008. The other fields run linearly between the means' centres. In the two years of the design case written
every 6.4 hours, the 6.4-hour means' mean speed runs 19–22% above the speed of their 3-day mean vector (21–24% in the
cube's ratio), so the wind is raised by the cube's ratio at each cell, after PlaSim's 1 m/s floor, in every target.

The bulk exchange is paired as PlaSim pairs it: the transfer coefficient at the lowest level's own height (its log of
the height over Charnock's roughness at lunar gravity), with that level's wind, potential temperature and humidity.
Until 18 October the lakes took the coefficient at 10 m with a 10 m wind from the neutral profile and the same level's
air, which exchanged 1.33–1.39 times PlaSim's heat for the same air. The 10 m wind is kept for Cole & Caraco's transfer
velocity, which is fitted to it, and for the wind's energy flux below.

The lakes lie 0.6–1.9 km below their cells' mean ground: the air is carried down along each target's own lapse rate
over land (0.75 K/km in the lock, 0.92 in the hot world, from climate_targets' ladder), the sky's infrared with it as
the fourth power of the air's temperature, and the pressure hydrostatically; the humidity is kept and capped at
saturation. The air the lake's exchange takes is the free air at the lowest level's height above the lake, warmer than
the cell's lowest level by the lapse rate times the depth below the cell; its potential temperature at the lake's own
surface pressure is its temperature times σ^−κ, as tas is the cell's, so carrying tas by the lapse rate is that to
within the factor σ^−κ = 1.005 on the carry, a few thousandths of a kelvin. A dry adiabat would carry it at 1.6 K/km;
the GCM's own column over land has the free-air lapse rate. The resort takes the lock's weather and moves each mean by
0.433 of the difference between the two runs' mean cycles by solar hour, the share of the way climate_targets uses; a
second resort built from the hot world's weather is a check.

**The model.** [biosphere/ecology/lake.py](../../../biosphere/ecology/lake.py) steps the column every 2 hours for the
deep lakes and every hour for the shallows, on fixed layers (0.5 m to 20 m depth, then growing to 8 m; 0.25 m in the
shallows) over the hypsography. Sunlight's near-infrared 58% is absorbed over 0.35 m (Paulson & Simpson 1977) and the
rest at the water's attenuation; light reaching the sloping floor heats the water over it. The surface loses infrared,
sensible and latent heat by PlaSim's own bulk scheme, paired as above: the neutral coefficient over Charnock's
roughness, Louis's stability factors and, over warm water, Miller et al.'s free convection, carried to lunar gravity.
The mixed layer deepens by GLM's energy balance (Hipsey et al. 2019): convective overturn and wind stirring (C_K 0.2,
C_W 0.23) against lifting the water below and accelerating it (C_T 0.51), with lunar gravity in every buoyancy term, so
the same wind mixes deeper through the same warmth than on Earth. The floor gives the Moon's heat flow, geography's
0.015 W/m². Ice forms when the surface reaches freezing; no lake freezes in any target, and the tests check it.

**The stirring below the mixed layer.** Below the mixed layer heat and solutes diffuse at Osborn's rule, 0.2 times the
dissipation over the squared buoyancy frequency at lunar gravity, on top of molecular diffusion. On 18 October the
author chose to drive that dissipation by the wind: a share η of the wind's energy flux 10 m above the lake,
P10 = ρ_air C10 U10³, spread over the mass of water below the mixed layer, M, as ε = η P10 / M. Across Earth's lakes
0.3 ± 0.1% of P10 reaches the stratified water as turbulence, 0.42% in Lake Alpnach, nine tenths of it spent in the
bottom boundary layer (Wüest, Piepke & Van Senden 2000), so η is 0.3% with 0.2% and 0.4% as its bracket; their mixing
efficiency is 0.15, against the model's 0.2. Under ice P10 is zero. C10 is Wüest & Lorke's (2003) drag over lakes, as
Simstrat pairs it with its own energy share (Gaudard et al. 2017). The lakes' winds are light, a mean 1.5–2.0 m/s at
10 m and 2.3–3.0 m/s in the windiest tenth, where that drag is high: the lake model's own surface stress at lunar
gravity, τ U10 with PlaSim's Charnock roughness at g/6, gives 0.51–0.78 of the same P10, the same as a share of
0.15–0.23% against Wüest & Lorke's drag, so the bracket's 0.2% covers only the upper half of that pairing; a constant
drag of 1.3 × 10⁻³ gives 0.61–0.85 of it; and at 4–10 m/s the model's own drag runs 1.5–1.8 times Wüest & Lorke's.
Which drag Wüest et al.'s 0.3% was measured against is open: a short search did not find it. The lakes' P10 averages
14–25 mW/m², against Lake Alpnach's 39. The mean dissipation runs 2.7 × 10⁻¹⁰ to 1.8 × 10⁻⁹ W/kg in the 352 m lakes
across the bracket, its median in proportion to the share; it runs higher in the hours when a deep mixed layer leaves
little water below it. The cap on the diffusivity is 0.1 m²/s. Osborn's rule passes it only where the water is close
to neutral, 43–62% of the interface-steps below the mixed layer, and there a step at 0.1 m²/s already mixes each layer
through with its neighbour (the cap times the step over the layers' spacing squared passes 1), so the cap holds back
none of the interfaces in any run, and a cap ten times higher leaves the 0–15° lakes as they are (above). Statically
unstable interfaces, 0.1–1.6% of them, mix at the cap.

Three treatments run beside the wind as sensitivity rows on the 0–15° lake: the wind's share under the first round's
cap of 10⁻⁴ m²/s; the first round's fixed dissipations, 10⁻¹⁰ and 10⁻⁹ W/kg under that cap; and DYRESM's diffusivity
(Yeates & Imberger 2003, as the Mono Lake boundary-layer study gives it), 200 times the molecular diffusivity of heat
times N²/N²_max over the Lake number, calibrated on Lake Kinneret and two reservoirs, under the wind's cap.

Oxygen crosses the surface at saturation for the lake's own pressure and warmth, at the surface-renewal transfer
velocity (Tedford et al. 2014, central) or Cole & Caraco's (1998); photosynthesis follows the light (tanh,
150 µmol/m²/s), scaled so that the cycle's mean production is the scenario's. The euphotic water respires four fifths
of the day's production from a labile pool that the production feeds and that decays with an e-folding of 4 days, with
2 and 7 days as a bracket: a plankton community's dark respiration falls through the night at 0.04 an hour towards a
steady level (Markager et al. 1992), and the pool's 2–7 days bracket the run-down of that steady part over a lunar
night. Respiration slows as oxygen falls toward 0.5 mg/L, its half-saturation, and what anoxic water cannot respire
becomes an oxygen debt of reduced substances that take oxygen back at a day's time scale where it returns.
Respiration's products make the water denser: aerobic respiration's CO2 at Lake Kivu's coefficient, 0.39 × 10⁻³ kg/m³
for each g/m³ of oxygen respired, leaving at the surface with the gas; and the debt's anaerobic products, half its
oxygen equivalents ending as methane (lake_anoxia's central share), its CO2 and methane's lightening at Kivu's
coefficients, a net 1.36 × 10⁻⁴ kg/m³ for each g/m³ of debt; when the debt is oxidised its products become aerobic
respiration's CO2. Each run carries 24 oxygen columns, four trophic states by two transfer velocities by two kinetics
and the labile pool's two brackets; its physics, the light's attenuation and the solutes follow one of them,
mesotrophic water in the deep lakes and eutrophic in the shallows; the variant runs keep the four central columns.

**The demand.** Bounded scenarios from Earth's lakes:

| Trophic state | Gross production | Secchi depth, attenuation | Sediment | Water below the euphotic zone |
|---|---|---|---|---|
| Oligotrophic | 200 mg C, 0.53 g O2 per m² a day | 6 m, 0.28/m | 0.25 g O2 per m² a day | 0.07, over 25 m |
| Mesotrophic | 600 mg C, 1.6 g O2 | 3 m, 0.57/m | 0.30 | 0.30 |
| Eutrophic | 2,000 mg C, 5.3 g O2 | 1 m, 1.7/m | 0.37 | 0.53 |
| New reservoir | 2,000 mg C, 5.3 g O2 | 1 m, 1.7/m | 0.49 | 0.82 |

Production is Wetzel's (2001) ranges (oligotrophic 50–300, mesotrophic 250–1,000, eutrophic over 1,000 mg C per m² a
day); the attenuation is 1.7 over Carlson's (1977) Secchi depths. The sediment and deep water take deep European lakes
(Steinsberger et al. 2020: sediment uptake 7.9–10.6 and water-column mineralisation 2.1–27.4 mmol O2 per m² a day, the
sediment's independent of trophic state) and eutrophic ones (Müller et al. 2012: 0.47–1.31 g O2 per m² a day, mean 0.90,
of which 0.37 ± 0.12 a flux of reduced substances, the water's share growing with the hypolimnion's thickness to about
25 m); the new reservoir, eutrophic for its first six to fourteen years (Paterson et al. 2019), takes Müller's highest.
The euphotic water respires four fifths of its production (a design value), so that eutrophic production exports what
Müller's eutrophic lakes consume below. Where a lake is shallower than its euphotic zone and 25 m, the deep water's
demand falls on the floor. The demand is fixed by the trophic state's supply of organic matter, whatever the water's
warmth, and the bounds rest on it. The deep and floor rates were measured in hypolimnia at 4–8 °C; the warm columns give
them ecosystem respiration's 0.62 eV (Yvon-Durocher et al. 2012) from 6 °C, with production and the pool held, 5.2
times the supply's rates at 25 °C. No settling supply sustains that, so the warm columns stand for the weeks to months
after a warming while a stock of settled matter lasts, a transient kept beside the bounds.

**Spin-up.** Deep water settles at the temperature of the densest water the surface makes, so each deep lake starts
1 K above the air's mean and saturated, and passes over the whole ten-year record in 8-hour steps until its deep
water's mean temperature, the coupled column's deepest-tenth oxygen, the bottom water's age and its solute each move
less than 0.005 K, 0.02 g/m³, 2% and 0.02 g/m³ from one pass to the next, or for twenty passes: the 352 m lakes settle
in 3–6 passes. The 8-hour steps settle a slightly different deep water from the 2-hour steps of the analysis, so each
deep lake then passes over the record at the analysis step until the same test passes there, or for ten passes, before
the record is analysed: the 352 m lakes settle there in 2–4 more. The review's re-check had found the first 2-hour
decade's bottom water 4–53% older than the 8-hour spin-up's last pass under the first round's mixing. A column whose
deep water has stopped exchanging keeps ageing and so does not pass the test. The 1,925 m lake starts from the 352 m
lake's column carried down and oxygenated, as if just overturned, and either runs the record from there or first passes
over it for up to 300 years (thirty passes) and then for up to 100 at the analysis step, to reach its own regime.

## The checks

- **The model.** Its tests ([biosphere/tests/test_lake.py](../../../biosphere/tests/test_lake.py)) reproduce the
  standard saturation (9.09, 6.95 and 14.62 mg/L at 20, 35 and 0 °C) and Weiss's independent salinity factor; keep heat
  and oxygen through a freeze and a thaw to 10⁻⁹, and with the exchange paired aloft, the labile pool, the debt's
  density and a raised cap; decay a cosine mode exactly as the implicit scheme should, to 10⁻¹⁰, and within its
  first-order error of the analytic decay; deepen a convecting mixed layer into linear stratification as
  h = (2Qt/ρcΓ)^½ within 2%, and a wind-stirred one as h³ = 6 C_W u*³ t/N² within 5%, deeper at lunar gravity by the
  cube root of Earth's over it; relax a mixed slab to saturation at its transfer velocity; grow ice by conduction as
  Stefan's law has it; integrate each scenario's demand; hold dense, CO2-laden deep water and water carrying the debt's
  anaerobic products against an overturn its warmth alone would allow; give the same stress from the exchange paired
  at 10 m and at the lowest level, with 1.25–1.45 times the heat exchange at 10 m; run the labile pool down as the
  implicit step has it and respire all it loses; cap the diffusivity, with unstable interfaces mixing at least as fast
  as the cap; carry and mix the solute where it adds no density; and leave the rates as given at the Arrhenius
  reference. For the wind's deep mixing they check that the dissipation times the water's mass below the mixed layer is
  the share times P10, in W/kg from W/m² over kg/m², with P10 on Wüest & Lorke's drag; that the scheme gives the fixed
  dissipation's column when its share gives the same dissipation; that under the wind a cap of 0.1 m²/s holds back no
  interface, where 10⁻⁴ m²/s holds back a tenth or more, and that a cap ten times higher leaves the deep water within
  0.01 K and 0.01 g/m³; that heat and oxygen close under it; and DYRESM's diffusivity's form. The study's own tests
  check the iron bracket's coefficients, net of the CO2 iron reduction consumes; that every deep lake passed over the
  record at the analysis step before its record was analysed (the 352 m lakes under the wind settle there in 2–4
  passes, while the Lake-number runs and the lock's giant at the central share and at 0.4% reached the ten-pass limit
  unsettled); that the wind's runs carry their share, their median dissipation in proportion to it, with no interface
  held back; that the raised cap's runs match the central ones; that the energetic regime's diffusivity takes Shih's
  three regimes, continuous at Re_b = 100 and Osborn's divided by 0.1 Re_b^½ above it, with lake.py's own functions
  back in place after the wrapper; and that mixing_forms' product binds its code and inputs, its runs keep their
  budgets and its energetic regime turns the 0–15° lakes over no more often than Osborn's rule and leaves their deep
  water no richer.
- **The runs.** Every run keeps its heat to 3.5 × 10⁻⁹ of the column's heat content, 10⁻⁷ W/m² over the run, and its
  oxygen to 2 × 10⁻⁷ of the column's oxygen, 1.3 × 10⁻¹⁰ g/m² a day, and no lake freezes: the implicit diffusion's
  round-off, which the wind's large diffusivities across near-neutral water raise. The fitted Sun reproduces the GCM's
  mean cosine of the zenith angle at the lake cells to 0.0006–0.008.
- **Against the GCM's own water.** A 50 m mixed slab forced as the lakes are, at the water cell nearest each lake,
  against the GCM's slab ocean, mean bias in K (correlation):

| Slab against the GCM's water | Paired at the lowest level, Earth's free convection | Lowest level, lunar free
convection | At 10 m, Earth's | At 10 m, lunar | |---|---|---|---|---| | The lock, 0–15° | +0.03 (0.985) | +1.49 (0.979)
| −0.54 (0.978) | +0.58 (0.966) | | The lock, 15–30° | +0.07 (0.983) | +1.51 (0.978) | −0.46 (0.974) | +0.66 (0.963) | |
The lock, 30–45° | +0.04 (0.991) | +0.87 (0.977) | −0.77 (0.975) | −0.10 (0.958) | | The hot world, 0–15° | +0.28
(0.973) | +1.77 (0.952) | −0.12 (0.959) | +1.08 (0.927) | | The hot world, 15–30° | +0.16 (0.977) | +1.45 (0.958) |
−0.32 (0.962) | +0.65 (0.933) | | The hot world, 30–45° | +0.03 (0.961) | +0.69 (0.936) | −0.64 (0.922) | −0.17 (0.890)
|

  Paired as PlaSim pairs it and with PlaSim's own Earth constant for the free convection, the slab follows the GCM's
  water to 0.03–0.28 K on average with correlations of 0.96–0.99: the lake model's exchange reproduces the GCM's own.
  The lakes take the constant carried to lunar gravity, which runs the slab 0.69–1.77 K warmer. The 10 m pairing
  exchanged a third more heat and ran the slab 0.12–0.77 K cooler with Earth's constant, which offset the lunar
  constant's warming before 18 October.
- **The 3-day means and the step.** The design case's two years written every 6.4 hours, its years 20–21, drive the
  0–15° lakes from the lock's settled columns once from those means and once from their own 3-day means with the wind
  raised as the main runs raise it, both at a half-hour step, and from the 3-day means at the main runs' steps:

| The design case's years 20–21 at the 0–15° lake | 6.4-hour means, half-hour step | 3-day means, half-hour step | 3-day
means, 2-hour and 1-hour steps | |---|---|---|---| | 352 m lake: complete overturns in 23 cycles | 13 | 14 | 13 | | Its
deepest tenth: oligotrophic; mesotrophic | 7.34; 6.11 mg/L | 7.22; 5.86 mg/L | 7.25; 5.91 mg/L | | 3 m eutrophic
shallows: median night's lowest; worst tenth | 4.53; 4.12 mg/L | 4.87; 4.25 mg/L | 4.91; 4.27 mg/L |

  The 3-day means and the main runs' steps leave the deep lake's overturns within one and its deepest tenth within
  0.25 mg/L, and run the eutrophic shallows' nights 0.13–0.39 mg/L richer than the 6.4-hour means. Those two years
  overturn in 13–14 of 23 cycles, where years 10–19 overturn in 66 of 120 at the same lake: the overturns' rate holds
  across the two records.
- **The GCM's free convection.** PlaSim's fluxmod.f90 carries Miller et al.'s (1992) free-convection transfer over
  water, 0.0016 Δθ^⅓ m/s, at Earth's gravity with no gravity factor; turbulent free convection scales as (gΔθ)^⅓, so at
  lunar gravity it would be 0.55 of it. In the slab, paired as PlaSim pairs it, Earth's constant cools the GCM's water
  cells by 0.66–1.49 K against the constant carried to lunar gravity, a question for the climate domain.

## Assumptions and limits

- **The lakes stir beyond what one column holds.** The Lake numbers above, under 1 on many days whichever thermocline
  they are read at, say the wind's stress tilts the thermocline to the surface and stirs the deep water by basin-scale
  motions, upwelling at the shores and internal waves, which a one-dimensional column cannot represent. The column
  stands for that stirring by the share of the wind's energy that Earth's lakes pass to their stratified water, spread
  evenly over the water below the mixed layer, where most of it is spent at the floor. The share was measured in Earth's
  lakes under Earth's gravity, with no lunar measurement behind it; whether a lake whose thermocline the wind tilts to
  the surface passes more of the wind's energy below is open, and so is the drag the share was measured against. The
  form that turns the dissipation into mixing is Osborn's rule at every buoyancy Reynolds number; Shih and Bouffard's
  energetic regime, which mixes the near-neutral deep water 3–30 times less, runs as the other end of the range
  (mixing_forms), and the share is spread evenly where Wüest et al. spend nine tenths of it in the bottom boundary
  layer, and taken at an efficiency of 0.2 against their 0.15: each leans towards more mixing. The lake kinds' ladder
  and census at 301 K, the seas' and the methane runs take Osborn's rule. The 3 m pans take the same share, measured
  below the thermoclines of deep lakes, over the little water below their warm layer, bracketed by half the share and
  the first round's fixed dissipation. Shear and billowing at the thermocline are left out, and their interface's
  reduced gravity is six times smaller at g/6.
- **Earth's ranges extended.** The Benson and Krause fits are stated for 0.5–1.1 atm and 0–40 °C. The Moon's 1.2 atm
  extends the pressure factor, which is thermodynamic, and the next term it omits (the Poynting correction for dissolved
  O2) is 3 × 10⁻⁴. The seas stay under 38 °C; only the lakes' sunlit surface layers pass 40 °C, their warmest
  hundredth of the time reaching 47 °C in the hot world's shallows, and there the fit is extrapolated in its own
  variable. The Schmidt number is Wanninkhof's at 20 °C carried by the Stokes–Einstein scaling past his fit's 30 °C.
  Cole & Caraco's transfer, the GLM coefficients, the wind's share and drag and the trophic states' rates are Earth's;
  Miller et al.'s free convection is carried to lunar gravity by its cube root.
- **One record, repeated.** The giants' deep water settles to the colder spells in the GCM's ten years. Longer records
  hold colder spells, which would set their deep water colder and their overturns rarer; the 352 m lakes, renewed every
  few months, depend on the record less, and the design case's years 20–21 overturn them at the rate of years 10–19.
- **Solutes.** Aerobic respiration's CO2 and the debt's anaerobic CO2 and methane reach the density; CO2 dissolving
  sediment calcite into bicarbonate and iron reduction's dissolved iron and bicarbonate are brackets, the iron carried
  with the debt and capped on it, so that a mixture of capped and uncapped water is up to twice as dense in iron as
  its parts; inflows, groundwater seeping through the floor and the lakes' own salts are left out.
- **The physics left out.** Inflows and outflows, wind sheltering by the crater rim, internal seiches and the
  sediment's own heat storage; one column stands for the whole lake.
- **Forcing.** 3-day means at T21 made hourly; the cells' air carried to the lakes' levels by the free-air lapse rate;
  the design case's wind factor raising every target's wind; the resort mixes two physics.
- **Biology.** Production and demand from Earth's lakes stand for engineered communities that have not been designed;
  the labile pool's e-folding is a design bracket; the 2 and 0.5 mg/L lines are the usual hypoxic and anoxic marks.
- **The targets as run.** The lake studies read the climate targets they ran on, domain/ecology's product at e9c12eb
  (schema /2), copied byte for byte to [climate_targets_as_run.json](climate_targets_as_run.json); the lakes' rerun on
  the resort's 301 K design case replaces it. The merged targets of the fourth joint integration take the lapse rates
  over the GCM's own land mask, 1.2–1.8% lower, which would cool these lakes' air by 0.026 K at most and 0.015 K
  typically; their lock's share of the way is 0 where the copy's is −0.002, the design case the lakes' lock already is;
  and the seas' medians the saturation reads fall 0.05–0.12 K (the lock 24.17 to 24.12 °C, the resort at 300 K 28.71 to
  28.59 °C), which would raise it by 0.006–0.013 mg/L.

## Running it

```sh
python -m research.studies.lake_oxygen.lakes                             # a few seconds: LOLA and the drainage product
climate/gcm/.venv/bin/python -m research.studies.lake_oxygen.extract      # about a minute; needs netCDF4
python -m research.studies.lake_oxygen.run                               # about 75 minutes on two cores
python -m research.studies.lake_oxygen.run --cache-only                  # the product from the cached runs, seconds
python -m research.studies.lake_oxygen.mixing_forms                      # about an hour on one core, after run
python -m research.studies.lake_oxygen.mixing_forms --cache-only         # its product from its cache, a minute
python -m pytest research/studies/lake_oxygen biosphere/tests/test_lake.py
```

The extract reads the GCM's output in `climate/gcm/runs` and writes `research/runs/lake_oxygen/forcing.npz` with
`forcing.json` (ignored, on the research drive), recording the sha256 of every file it read; the run caches each lake in
`research/runs/lake_oxygen/cases`, so a run stopped part way resumes, and `--quick` runs two years to
`research/runs/lake_oxygen/lake_oxygen_quick.json`. The cache key holds run.py and lake.py whole, the case and what the
runs read from their inputs (`Setting.values_fingerprint()`: this study's extracts by their sha256, and the targets'
shares and lapse rates, the held seas' salinity and pH and the floor's heat flow by their values), so a product
written again with the same values, as the settle cascade writes them, leaves the cache standing; `--cache-only`
builds the product from the cache and stops at the first run it would have to compute, here and in every lake study.
mixing_forms.py saves its deep columns daily in the lake kinds' format to `research/runs/lake_oxygen/mixing_forms`,
indexed with each file's sha256 in `mixing_forms_index.json`, and resumes pass by pass.
