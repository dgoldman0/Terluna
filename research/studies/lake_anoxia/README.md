# What the lakes' anoxic depths do, and what design can do about them

17 October 2026; rerun on 18 October on the lake model as corrected after that day's independent review and its re-check
(the findings, the re-check and their scripts are on the research drive, `review-lakes-2026-10-18`), and with the mixing
below the mixed layer driven by the wind, the treatment the author chose that day. [The waters'
oxygen](../lake_oxygen/README.md) finds the wind's energy below the night's mixing turning the Moon's 352 m crater lakes
over every second lunar month in the lock (295 K) and the resort (300 K) and every fourth in the hot world (near 307 K)
by Osborn's rule for the mixing: their oligotrophic and mesotrophic deep water stays oxic in the two cooler targets,
eutrophic water and new reservoirs run low, most of all in the hot world, and the 1,925 m lake that half the lake area
exceeds turns over about once a year with its richer depths anoxic. By Shih and Bouffard's energetic regime, the other
end of its range, they turn over about once a year, oligotrophic deep water keeps its oxygen in the cooler targets and
the lock's mesotrophic water keeps 2.9 mg/L, and richer water runs low or anoxic. This study's columns take Osborn's
rule; the energetic regime's methane comes from the oxygen study's [mixing_forms.py](../lake_oxygen/mixing_forms.py),
which runs this study's arithmetic on its 0–15° lakes. The author
asked what anoxic water does and what design can do about it: the lakes' methane and what it does to the air and its
dimming; whether CO2 and methane gather towards the hydrostatic bubble point behind Lake Nyos's eruption of 1986 and the
gas held in Lake Kivu; phosphorus released from anoxic floors; and the remedies, sized. [profiles.py](profiles.py)
reruns the lake model (`biosphere/ecology/lake.py`) from the oxygen study's own settled columns, at the wind's central
share and at its bracket, and saves each column daily; [run.py](run.py), with the gas physics in [gases.py](gases.py),
writes [results/lake_anoxia.json](results/lake_anoxia.json) (schema `terluna.research.lake-anoxia/2`); sources and how
far each was read are in [sources.json](sources.json). The corrections of 18 October bound what each overturn releases
by two rules, the inventory's net rise since the one before and the stock it brings up within the methane made since
then; count methane made in anoxic water alone; take the lasting oxycline's escape from Lake Kivu; bracket the giants by
the oxycline's flux alone; bracket the phosphorus retention's settling for lunar gravity and bound the internal load by
the steady state; time the bubble point from the mineralisation rate; compare with Earth under Earth's own air; and keep
the warm kinetics, which no settling supply sustains, beside the oxygenation's bounds as a transient. The seas' basins
are outside this study.

## What it finds

**The lakes' methane turns on their trophic state, the warmth and the mixing's form.** Methane is the share of the
anoxic water's and floors' mineralisation that respiration cannot meet with oxygen; the lake model carries it as its
debt of reduced substances, and in water as poor in sulphate as the Moon's lakes (about 1.4 µM) methanogenesis takes
about half of it (central 0.5, bracket 0.3–0.8), as in ferruginous Lake Matano. Moon-wide, over geography's 3.92 million
km² of rain-fed lakes, in Tg of methane a year, at the wind's central share of its energy below the mixed layer (0.3%)
and at the share's bracket:

| Tg a year, central (bracket), each overturn releasing its net rise; its stock; the lakes over 1,000 m on the oxycline's flux alone, by each | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic, 0.3% (central), 0.2% and 0.4% alike | 0.37 (0.27–0.38); 0.37 (0.27–0.38); 0.37, 0.37 | 0.34 (0.25–0.34); 0.34 (0.25–0.34); 0.34, 0.34 | 0.29 (0.23–0.30); 0.29 (0.23–0.30); 0.29, 0.29 |
| Mesotrophic, 0.3% | 2.4 (1.3–2.5); 2.4 (1.3–2.5); 2.4, 2.4 | 2.3 (1.3–2.4); 2.3 (1.3–2.4); 2.3, 2.3 | 2.4 (1.2–5.5); 2.5 (1.2–5.6); 2.3, 2.3 |
| Mesotrophic, 0.2% | 2.4 (1.3–2.5); 2.4 (1.3–2.5); 2.4, 2.4 | 2.3 (1.3–2.4); 2.3 (1.3–2.4); 2.3, 2.3 | 3.8 (1.3–16); 3.9 (1.3–16); 2.8, 2.8 |
| Mesotrophic, 0.4% | 2.4 (1.3–2.5); 2.4 (1.3–2.5); 2.4, 2.4 | 2.3 (1.3–2.4); 2.3 (1.3–2.4); 2.3, 2.3 | 2.1 (1.2–2.2); 2.1 (1.2–2.2); 2.1, 2.1 |
| Eutrophic, 0.3% | 4.3 (2.6–4.4); 4.3 (2.6–4.4); 4.3, 4.3 | 4.2 (2.5–4.4); 4.2 (2.5–4.5); 4.2, 4.2 | 6.6 (2.7–28); 7.1 (2.7–29); 5.1, 5.2 |
| Eutrophic, 0.2% | 4.3 (2.6–4.5); 4.3 (2.6–4.5); 4.3, 4.3 | 4.8 (2.6–8.9); 4.8 (2.6–9.1); 4.4, 4.4 | 15 (7.5–50); 23 (12–61); 7.9, 10 |
| Eutrophic, 0.4% | 4.3 (2.6–4.4); 4.3 (2.6–4.4); 4.3, 4.3 | 4.2 (2.5–4.3); 4.2 (2.5–4.3); 4.2, 4.2 | 5.7 (2.6–18); 5.8 (2.6–18); 4.7, 4.7 |
| New reservoir, 0.3% | 8.5 (4.6–14); 8.6 (4.6–15); 8.1, 8.2 | 11 (4.8–36); 13 (4.9–41); 9.1, 9.4 | 27 (14–81); 46 (24–110); 14, 20 |
| New reservoir, 0.2% | 12 (4.9–37); 13 (4.9–42); 9.3, 9.7 | 20 (8.0–72); 32 (14–93); 12, 16 | 25 (13–76); 52 (29–112); 14, 22 |
| New reservoir, 0.4% | 7.9 (4.5–8.3); 7.9 (4.6–8.5); 7.9, 7.9 | 8.7 (4.6–17); 9.0 (4.6–18); 8.1, 8.2 | 18 (7.0–72); 26 (11–90); 11, 14 |

Every lake takes the 352 m lakes' flux from the water per area by trophic state, with each floor's own bubbles by its
local depth. What an overturn releases is bounded by two rules until a methane tracer settles it. The net rise counts
the inventory's rise since the previous overturn, which in a repeating record is what the model's own oxygen removed
during the overturns: a lower bound. The stock counts the standing stock the overturn brings up, no more than the
methane made since the previous overturn: the study's own assumption taken whole, an upper bound. The share of the
stock the restratifying mixed layer degasses lies between them. The lakes deeper than 1,000 m hold 70% of the area, and
the oxygen study's 1,925 m lake turns over about once a year with its richer depths anoxic between, so the last figures
give them only what crosses a lasting oxycline. The bracket runs the methane share, the bubbles' size and surface, and
the oxidation and escape shares below from their low ends to their high.

The domain synthesis's ledger carries 39 Tg for the lakes (7.8–118), Earth's 10 g m⁻² a year over their area. The Moon's
lakes emit, central, 0.075–0.095 g m⁻² a year oligotrophic, 0.58–0.63 mesotrophic, 1.1–1.8 eutrophic and 2.2–12 in new
reservoirs, from the net rise to the stock. In the lock and the resort, by Osborn's rule, the overturns keep the 352 m
lakes' oligotrophic, mesotrophic and eutrophic depths oxic, so their deep water makes almost no methane and the lakes'
methane is the floors'
own bubbles; their anoxic depths make 0.16 mol of methane per m² a year mesotrophic in the hot world, 0.99 eutrophic,
and 0.28, 1.0 and 2.8 in new reservoirs from the lock to the hot world. Lake Tanganyika's floors give its anoxic water
0.5 mol m⁻² a year (Durisch-Kaiser et al. 2011). At the share's 0.2% the hot world's mesotrophic depths make 0.52 mol
and its eutrophic ones 1.7, and the resort's eutrophic ones 0.18; at 0.4% the hot world's mesotrophic depths stay oxic.
The warmth shows in the hot world's richer lakes, whose floors lie under anoxic water for 222–306 days of the year
eutrophic and 307–355 in new reservoirs.

By Shih and Bouffard's energetic regime the 0–15° lake's anoxic depths make 0.16 mol of methane per m² a year eutrophic
and 1.4 in new reservoirs in the lock, 0.64 and 2.6 in the resort, and 1.8 and 2.8 in the hot world, where its
mesotrophic depths make 1.1. With that lake standing for every lake, the Moon's lakes would send 4.4–5.1 Tg a year
eutrophic in the cooler targets and 20–37 in the hot world, against 4.2–4.3 and 7.8–8.4 by Osborn's rule on the same
lake; 9.7–34 and 31–69 from new reservoirs, against 7.8–9.6 and 38–58; mesotrophic 2.3–2.4 and 4.7–5.5, against 2.3–2.4
and 2.9; oligotrophic 0.33–0.58, against 0.28–0.36 (centrals, from the net rise to the stock). The methane conclusions
run between the two forms: oligotrophic and mesotrophic lakes stay within a few tenths of a Tg in the cooler targets,
the hot world's eutrophic lakes carry two to four times as much by the energetic regime and its mesotrophic ones nearly
twice, and the resort's new reservoirs up to three and a half times, most of it at the stock rule's upper end.

**What the overturns release lies between two bounds.** Between overturns at the central share the anoxic water stores
under 0.01 mol of methane per m² in the lock's and the resort's eutrophic lakes, 0.14–1.0 in the hot world's, and up to
0.58 in the resort's new reservoirs and 1.0–3.9 in the hot world's, with 17–28 µM at the floor in the hot world's
eutrophic lakes and 46–63 µM in its new reservoirs, below the 170–220 µM of Tanganyika's depths. By the net rise the
overturns release 0.24 mol per m² a year in the hot world's eutrophic lakes and 0.059, 0.37 and 0.89 in new reservoirs
from the lock to the hot world; by the stock 0.29, and 0.072, 0.49 and 1.7. Each overturn's oxygen meets all of it
except in the hot world's new reservoirs, where it meets 40–45%. Earth's lakes oxidise 54–98% of their stored methane at
the autumn overturn, while oxygen enters as the mixed layer deepens over weeks; the Moon's deep lakes turn over in one
cold night. What the overturned lake's oxygen meets is oxidised at Earth's central 85%, and half of what it cannot meet
escapes. The rest of the methane made crosses a lasting oxycline, through which Lake Kivu lets 0.2 of the 35 g C m⁻² a
year reaching it escape (Pasche et al. 2011): 0.57% central, from none to the 15% Earth's overturns let through.

**Bubbles reach six times deeper on the Moon, and carry most of the methane.** A floor bubbles where its porewater
methane reaches the local pressure less its nitrogen; at a sixth of Earth's gravity the same pressure lies six times
deeper, and under the Moon's 1.2 atm, whose 0.95 atm of nitrogen nearly offsets the extra surface pressure, a floor
bubbles as Earth's does at about six times the depth; Earth's floors bubble mostly in their top 10–20 m. The steady
porewater balance behind this (Bazhin 2010, as Langenegger et al. 2019 apply it) gives 11% of Lake Soppen's production
as bubbles at 10 m, against the 10–17% measured. Bubbles rise more slowly at lunar gravity (a 5 mm bubble at 0.15 m/s
against 0.23) and lose more of their methane on the way: from 30 m a 10 mm bubble keeps 87% with a contaminated surface
and 44% clean, and a 5 mm one 67% and 5%. 94% of the lake area lies in lakes deeper than 100 m and only 4.2% of the
floor is shallower than 10 m. With the depths oxic, the bubbles carry all of the oligotrophic and mesotrophic lakes'
methane in the lock and the resort and of the eutrophic lakes' there, 87% of the hot world's mesotrophic lakes', 58–62%
of its eutrophic lakes' and 16–93% of new reservoirs', from the hot world to the lock. On Earth, whose lakes are small
and shallow, they carry about half (Bastviken et al. 2011).

**What it does to the air and the dimming.** In the bare column under the stack the air holds no OH (the joint
synthesis), so the soils' 2.92 Tg a year for each ppm (1.07–4.77) is its only sink; a living column may make some OH
(ecology, finding 9). With the lakes in the ledger at these rates and every other source cut 90% by biology, as the
held pair assumes, at the wind's central share:

| Steady methane, ppm; removal that holds 1.87 ppm, Tg a year (share of the air passed a year); added dimming if left | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic, net rise or stock | 0.94–1.7; 0; −0.23 to −0.04% | 0.93–1.7; 0; −0.22 to −0.04% | 0.91–1.7; 0; −0.22 to −0.04% |
| Mesotrophic, net rise or stock | 1.6–2.4; 0–1.5 (0–0.05%); −0.04 to 0.10% | 1.6–2.3; 0–1.4 (0–0.04%); −0.05 to 0.09% | 1.7–2.4; 0–1.5 (0–0.05%); −0.04 to 0.10% |
| Eutrophic, net rise | 2.3–3.0; 1.2–3.4 (0.04–0.10%); 0.09 to 0.21% | 2.3–3.0; 1.1–3.3 (0.03–0.10%); 0.08 to 0.20% | 3.1–3.8; 3.5–5.7 (0.11–0.17%); 0.20 to 0.31% |
| Eutrophic, stock | 2.3–3.0; 1.2–3.4 (0.04–0.10%); 0.09 to 0.21% | 2.3–3.0; 1.2–3.3 (0.04–0.10%); 0.08 to 0.20% | 3.2–4.0; 4.0–6.2 (0.12–0.19%); 0.22 to 0.33% |
| New reservoir, net rise | 3.7–4.5; 5.4–7.6 (0.16–0.23%); 0.32 to 0.43% | 4.7–5.5; 8.4–11 (0.25–0.32%); 0.45 to 0.53% | 10–11; 24–26 (0.74–0.80%); 0.87 to 0.91% |
| New reservoir, stock | 3.8–4.5; 5.5–7.7 (0.17–0.23%); 0.32 to 0.43% | 5.1–5.8; 9.4–12 (0.29–0.35%); 0.49 to 0.57% | 17; 43–45 (1.32–1.38%); 1.20 to 1.22% |

Each pair is the domain synthesis's low- and high-cost ends; the dimming is the forcing beyond the held pair's
3.58 W/m² in the deep column, with N2O at 0.1 ppm, over each target's absorbed sunlight (211–225 W/m²), and comes on
top of the 1.69% the pair already takes; a negative figure is dimming the air's lower methane spares. Oligotrophic lakes
left as they are hold the air at 0.91–1.7 ppm, under the held 1.87 ppm, at every share; mesotrophic lakes at 1.6–2.4
(2.9 in the hot world at the share's 0.2%), which a removal of up to 1.5 Tg a year holds at the held level (3.0 at
0.2%). Eutrophic lakes left alone put the air at 2.3–4.0 ppm, and at the share's 0.2% the hot world's at 6.1–9.4, with
up to 22 Tg a year of removal to hold it; new reservoirs left alone put it at 3.7–17 ppm, 19 ppm at the share's 0.2%.
By Shih and Bouffard's energetic regime, with the 0–15° lake standing for every lake, uncut lakes in the hot world hold
the air at 7.6–14.2 ppm eutrophic and 11.4–25.3 ppm as new reservoirs, against 3.5–4.4 and 13.7–21.4 ppm by Osborn's
rule on the same lake; the cooler targets barely change. With the lakes' methane cut 90% too, every state holds the
air at 0.82–3.3 ppm across the shares, with at most 4.3 Tg a year of removal.

**No lake in the targets approaches a limnic eruption, but lunar gravity leaves little headroom.** The gas saturation
ratio is the dissolved gases' partial pressures with the vapour over the local pressure; bubbles can form where it
reaches 1. Lake Kivu's deep water reaches about 0.5 at 320 m (Bärenbold et al. 2020). In the 352 m lakes the floor's
water reaches 0.15–0.17 at its 6.8 atm, its gases at 1.05–1.17 atm, the air's nitrogen and argon about 0.95 atm of it.
The highest ratio, 0.70–0.77, sits at the hypolimnion's top, in the layer at 20 m where the ratio is first read: the
local pressure there is 1.5 atm and the water still carries the gases it took at the surface, so the ratio comes from
the air's gases more than the lakes' own. The 1,925 m lake turns over 86–105 times in its fifty years from an
oxygenated start, so its deep water gathers no gas toward its 32 atm bubble point; isolated, at its trophic state's
supply-fixed demand over a 50 m hypolimnion (the hazard map's rule), it would reach it in 1,880–2,290 years oligotrophic
and 670–810 eutrophic. An isolated hypolimnion reaches its bubble point 3.5–6.7 times sooner than on Earth at
the same depth under Earth's own air, because its pressure grows six times more slowly:

| Years for an isolated hypolimnion to reach its bubble point at the floor: lock, hot world (Earth at 25 °C) | Oligotrophic | Eutrophic |
|---|---|---|
| Floor at 20 m, 5 m of hypolimnion | 4.3, 3.5 (15) | 1.5, 1.2 (5.3) |
| 30 m, 15 m | 16, 13 (65) | 5.9, 4.8 (23) |
| 50 m, 35 m | 55, 45 (244) | 20, 16 (87) |
| 100 m, 50 m | 139, 113 (679) | 49, 40 (241) |
| 352 m, 50 m | 442, 361 (2,348) | 157, 128 (835) |

The hypolimnion takes the trophic state's floor and deep-water demand over its thickness and exchanges nothing. So the
hazard lies where a lake tens of metres deep never turns over: a meromictic lake held by a salt or soda layer, as closed
basins become soda lakes over millennia (sea chemistry), a lake fed at depth by mineral springs as Kivu is, or a deep
new reservoir in its eutrophic first decade. There a eutrophic hypolimnion at 20–50 m reaches its bubble point within
one to twenty years. Its gas would be mostly methane, which rises and disperses, and burns if lit, where Nyos's CO2
pooled. Nyos's magma fed it 80 mol of CO2 per m² a year (Kling et al. 2005); a eutrophic lunar lake's depths, wholly
anoxic, would make about 2.6 mol of methane and as much CO2.

**Phosphorus: anoxic floors can at most cancel the lakes' retention.** On the land's dissolved phosphorus alone,
0.018–0.048 g m⁻² a year, Vollenweider's relation puts the 352 m lakes at 1.1–3.0 mg/m³ of total phosphorus with Earth's
settling, ultra-oligotrophic. Its retention stands for settling, which slows at lunar gravity: with settling at the
square root of g/6 or at g/6 the same load gives 2.4–6.4 or 4.4–12 mg/m³. Nürnberg's anoxic release, 6.1 mg m⁻² a day
mesotrophic and 10.8 eutrophic, the means of her trophic classes taken as inputs that switch release on at the
mesotrophic boundary by construction, acts only where the floors go anoxic: at the central share no mesotrophic or
eutrophic floor in the lock does, the resort's eutrophic floors for 8 days a year and the hot world's mesotrophic and
eutrophic ones for 71 and 256, and new reservoirs' for 76, 225 and 338 from the lock to the hot world. There it gives up
to 2.8 g m⁻² a year eutrophic in the hot world and 3.6 in its new reservoirs, up to 200 times the dissolved load at
the release's top (new reservoirs' flooded soils, at her hypereutrophic 25.6, up to 8.6 g). Added outside the
retention, as her mass balance has it for lakes drawing on a floor's stock (Nürnberg 1984), that would put those lakes
at 6–2,050 mg/m³. At steady state a floor releases at most what settled on it, so whatever the floors do the lake
keeps at most what flows in, the load over the throughflow: 10–27 mg/m³ on the dissolved load, inside the mesotrophic
range. Anoxic floors can thus at most cancel the retention, which holds an oxic lake 2.3–8.9 times lower by the
settling's bracket, and burial, the vivianite below among it, keeps part of it; Nürnberg's figures hold while a floor
gives up a stock beyond its steady state, as a new reservoir's flooded soils do through its first six to fourteen
eutrophic years (Paterson et al. 2019). The Moon's low sulphate leaves iron to hold phosphate as vivianite (Caraco et
al. 1989), which the bracket takes as a tenth of Earth's release; oligotrophic floors release almost none whatever
their oxygen (Nürnberg 1997). The trophic state is a tipping point: while their floors stay oxic, as the overturns keep
them in oligotrophic and mesotrophic water in the cooler targets, the lakes stay oligotrophic on their dissolved load
with Earth's settling or its square root at g/6, and reach 12 mg/m³ at the high load with Stokes settling at g/6, the
load 1.2 times its 0.041 g permissible; where their floors go anoxic, in the hot world's richer lakes and in new
reservoirs, the same load can bring them into the mesotrophic range; and anything that makes them mesotrophic (the
particulate phosphorus becoming available, people's losses, a new reservoir's flush) starts a release that carries them
further in the hot world.

## The remedies, sized

| Per km² of lake, and for all 3.92 million km², at the wind's central share | The lock | The resort | The hot world |
|---|---|---|---|
| Oxygenating the depths: oxygen, t a day; up to, as the devices' stirring raises it | 0.32–1.31; up to 0.96–3.9 | the same | the same |
| Its power, kW (all lakes, GW) | 3.5–14 (14–57) to 14–58 (56–227) | 3.5–14 (14–57) to 14–58 (56–229) | 3.6–15 (14–57) to 14–59 (56–230) |
| Circulation, Earth's guideline air flow compressed to 352 m; to 1,925 m, kW (GW) | 54 (212); 102 (400) | the same | the same |
| Cold inflows routed to depth: oxygen, share of the deep demand | 0.044 g m⁻² a day, 3–14% | 0.040 g, 3–13% | 0.036 g, 3–11% |
| Degassing pipes where the depths hold methane: water, m³/s; pumping against the methane's energy, kW | new reservoirs: 3.2; 229 against 8 | eutrophic and new: 3.0–3.5; 217–249 against 0.15–29 | mesotrophic to new: 1.5–2.1; 108–152 against 4–78 |

Respiration's warm kinetics from the 6 °C its rates were measured at would raise the stirred demand 3.5 times in the
lock, 5.3 in the resort and 8.4 in the hot world, to 50–495 kW per km²; no settling supply sustains that, so it stands
for the weeks to months after a warming, a transient beside these bounds.

- **Nutrient targets.** The deep lakes stay oligotrophic under 0.16 g of bioavailable phosphorus per m² a year (10 mg/m³
  at their 1.8 m a year of throughflow) with Earth's settling, 0.075 with settling at the square root of g/6 and 0.041
  at g/6, and ultra-oligotrophic under 0.063, 0.030 and 0.016; where their floors go anoxic the retention can at most be
  cancelled, the oligotrophic line then falling towards 0.018 g, the dissolved load's lower end. Their dissolved load is
  a ninth to a third of the 0.16 g and 0.4–1.2 times the 0.041 g; the particulate load, kept in the deltas, must stay
  unavailable, which dredging the deltas for people's phosphorus (resources' plan, 4–51% of what the lakes receive)
  helps. Oligotrophic lakes are also the ones whose methane leaves the air under the held level by either rule, with no
  removal at all.
- **Oxygenating the depths.** The oxygen is the deep demand by trophic state (0.32–1.31 t per km² a day), fixed by the
  supply of settling matter, and up to three times it as the devices' stirring raises it (Newman Lake's design took
  three). Earth's installations deliver 0.29–3.6 t per km² a day (Moore et al. 2015). Oxygen by vacuum swing adsorption
  from the air costs 0.20–0.24 kWh a kilogram, and pressing it down to 352 m or 1,925 m at lunar gravity 0.06–0.12 more,
  half to two thirds of Earth's at the same depth. Line diffusers suit the deep lakes: their bubbles dissolve over the
  long lunar rise, as the methane bubbles' survival above shows. An oxic floor stops its methane. Lakes Baldegg and
  Sempach released as much phosphorus from their floors through ten years of oxygenation and mixing as before (Gächter
  & Wehrli 1998), so the nutrient targets come first and oxygenation backs them up. With the overturns keeping the
  352 m lakes' depths oxic in the lock and the resort, oxygenation is for the hot world's richer lakes, new reservoirs
  and the giants' depths.
- **Circulation.** Lorenzen & Fast's 9.2 m³ of air a minute per km², Earth's guideline for lakes tens of metres deep,
  compressed to 352 m costs 54 kW per km² on the Moon against 92 on Earth, and to 1,925 m 102 against 136, at 60%
  compressor efficiency; no plume is computed, and whether that air mixes a lunar lake hundreds of metres deep is open.
  The hot world's columns hold a quarter more Schmidt stability than the lock's (732 against 579 J/m²) and warm their
  deep water, which raises the demand.
- **Routing cold inflows.** The lakes' throughflow, 1.8 m a year, carries 3–14% of the deep demand's oxygen if all of
  it is taken to depth; the work against its buoyancy is a few watts per km², so it is cheap help for oligotrophic
  lakes and little for productive ones.
- **Lake depth by design.** A lunar night mixes a lake to 9–14 m even in the calmest tenth of lunar cycles, so floors
  within about 12 m are ventilated through most nights, and the overturns reach deeper floors every second to fourth
  lunar month. A reservoir or park lake shallower than that keeps an oxic floor with no running energy when its water is
  oligotrophic, and the 3 m pans' floors stay oxic in every state under the wind's mixing (lake_oxygen). One 20–50 m
  deep must not be held stratified for years: it reaches its bubble point in 1–20 years if eutrophic. Making the deep
  crater lakes shallow is out of reach: bringing the 352 m lakes' floors up to 12 m takes a mean 53–129 m of fill under
  each square metre of lake, 340 m at the deepest point.
- **Degassing.** Pipes as at Nyos lift deep water through a degasser. The lunar lakes' deep water holds 1.05–1.16 atm of
  gas, which lifts itself in no lake; pumped, the pipes pay only where the depths hold methane: in the hot world's new
  reservoirs they move 1.7 m³/s per km², pumped at 123 kW against 78 kW of methane fuel, Kivu's way, and where the
  overturns keep the depths oxic the water holds too little methane to pay for its pumping. Methane burnt to CO2 warms
  the air far less. Degassing keeps meromictic and spring-fed lakes from their bubble point; in the overturning deep
  lakes oxygenation does more.

## How it is set up

- **The columns.** The nine 352 m lakes (three bands, three targets) start from the oxygen study's settled columns
  under the wind's deep mixing, at its central share of P10 (0.3%) and at its bracket (0.2% and 0.4%), and run their
  ten years again, saved daily with the solutes, the labile pool and every column's surface oxygen at every step; the
  1,925 m lake runs five passes, fifty years, from its oxygenated start, at each share. The oxygen study's columns
  settle at the analysis step after their coarse spin-up, so these passes repeat a settled cycle. Each trophic state is
  the central column, surface-renewal transfer with demand fixed by the trophic state; the warm kinetics, from the
  source hypolimnia's 6 °C, give a transient beside the oxygenation's bounds. The lake model carries the debt's
  anaerobic CO2 and methane in the density since the review of 18 October, and the euphotic water's respiration from a
  labile pool.
- **Methane and CO2 from the oxygen.** Methane is the model's debt times the methane share, over 64 g of oxygen a
  mole; CO2 is the water's oxygen deficit below the euphotic zone and its debt's carbon (a mole of CO2 to every mole of
  methane, and one to each 32 g of the rest), added to lake water in equilibrium with the design air at the sea
  chemistry's lake pH, 7.72, an alkalinity of 0.42 meq/kg. N2 and Ar keep the surface's saturation. The methane made is
  the model's demand in each layer times the share respiration cannot meet, K/(K + O2) with its 0.5 g/m³, in anoxic
  water under 0.5 g/m³; between 0.5 and 2 g/m³ the model owes part of its demand and repays it in place within a day,
  which is aerobic respiration in a lake.
- **What an overturn releases.** The model keeps its debt through an overturn: the stock is mixed through the lake and
  waits for oxygen. An overturn here is a run of steps with the mixed layer at the floor, runs under fifteen days apart
  counted as one. By the net rise each overturn releases the inventory's rise since the one before
  (the first since the record's last, the record repeating); by the stock, the stock it brings up; by both, no more than
  the methane made in between. What it releases meets the overturned lake's oxygen, of which the share Earth's autumn
  overturns oxidise (85%, 50–98%) acts as far as it reaches, and what the oxygen cannot meet leaves or is oxidised
  later, half each. The methane made between overturns crosses a lasting oxycline, through which Lake Kivu lets 0.2 of
  the 35 g C m⁻² a year reaching it escape, 0.57% (Pasche et al. 2011), the central share; the bracket runs from all of
  it oxidised to the 15% Earth's overturns let through. The giants' figure gives the lakes deeper than 1,000 m only what
  crosses the oxycline, as if they never turned over.
- **Ebullition.** Production falls off exponentially in the sediment (30 per m) at the trophic state's floor demand
  (and the deep water's, within 25 m of the euphotic zone), the methane share of it; bubbles 80% methane, 10 mm across
  (5–20 mm, since weaker buoyancy lets bubbles grow larger before they break free of the mud, a guess on Earth's
  2–8 mm), rise with Tomiyama's drag for contaminated bubbles and lose gas by Frössling's transfer (contaminated) or
  Higbie's (clean); the bracket's low end takes small clean bubbles, its high end large contaminated ones. Each floor
  takes the mixed layer's temperature above 15 m and the deep water's below.
- **Scale-up.** Geography's lakes on LOLA give the floor by local depth (1.1% under 2 m, 4.2% under 10 m, 23% under
  100 m, 49% under 400 m); lakes over 1,000 m deep hold 70% of the area and those over 100 m 94%. Every lake takes the
  352 m lakes' flux from the water per area, with each floor's own bubbles, and the giants' figure replaces the deepest
  lakes' flux from the water with the oxycline's.
- **The bubble point.** Water isolated at a depth loses its oxygen at its mineralisation rate and then goes on
  mineralising anaerobically, half its oxygen equivalents ending as methane and the rest as CO2 in the lake's carbonate
  system, until its gases with the vapour reach the local pressure. The hazard map takes the trophic state's floor and
  deep-water demand over a hypolimnion from the floor to 15 m below the surface, at most 50 m thick; the 1,925 m lake,
  which turns over, takes the same rule at its own floor and deep warmth for the years it would take isolated. Earth's
  comparison takes Earth's gravity at 25 °C under Earth's own dry air at 1 atm.
- **Phosphorus.** Vollenweider's TP = L/(q_s(1 + √τ)), with the lakes' 793 m mean depth and 445 years of residence
  (1.8 m a year); the 352 m lakes take that throughflow over their 112 m mean depth. Its retention term stands for
  settling, an apparent settling speed of q_s√τ, 14 m a year in the 352 m lakes on Earth; settling slows with gravity,
  as Stokes settling does (g/6, 2.3 m a year) or as aggregates do (about √(g/6), 5.7), so the retention is bracketed
  by the three. The internal load, Nürnberg's release rate times the anoxic factor (the days a year a floor area equal
  to the lake's surface lies under water below 1 g/m³), enters outside the retention, as her mass balance has it
  (Nürnberg 1984), and at steady state the lake keeps at most its load over its throughflow. Her release rates are the
  means of her trophic classes, inputs that switch release on at the mesotrophic boundary by construction; new
  reservoirs take the eutrophic class with their production, and the hypereutrophic 25.6 mg m⁻² d⁻¹ as the bracket for
  their flooded soils.

## Assumptions and limits

- The lake model carries the methane and CO2 of its debt only in the density, as their net effect, and no phosphorus:
  the methane, CO2 and phosphorus here follow from its oxygen and debt offline, so methane leaves the water only by the
  bracketed shares.
- The stirring below the night's mixing is the oxygen study's: a share of the wind's energy measured in Earth's lakes,
  0.3% with 0.2% and 0.4% as its bracket, and every figure here stands at all three; the lakes' Lake numbers, under 1
  on many days whichever thermocline they are read at, say basin-scale motions a column cannot hold stir them, and
  whether they pass more of the wind's energy below than Earth's lakes do is open. Its form is Osborn's rule, which
  leans towards more mixing; the energetic regime's figures above come from the 0–15° lakes alone.
- Every lake takes the 352 m lakes' flux from the water. The giants over 1,000 m, 70% of the area, turn over about once
  a year with their mesotrophic and richer depths anoxic between, so they make more methane at depth than the 352 m
  lakes and store it behind an oxycline for longer; the giants' figure counts only what crosses that oxycline, and what
  their own overturns release is the methane tracer's to settle.
- What an overturn releases lies between the net rise and the stock until a methane tracer settles how much of the stock
  the restratifying mixed layer degasses. The oxidation shares and the methane share are Earth's ranges carried to the
  Moon; the half of the unmet methane that escapes is a guess. Bubble size at lunar gravity is a guess. Earth's anoxic
  release rates come from sulphate-bearing lakes.
- The isolated-hypolimnion years take the trophic state's demand over up to 50 m and no exchange; real lakes that
  never turn over still lose gas upward by diffusion, as Tanganyika and Kivu do.
- The seas' basins, the remedies run in the model and the biology of methanotrophs designed for sudden overturns are
  outside this study.
- The climate targets are the ones the lakes ran on, domain/ecology's at e9c12eb (schema /2), which lake_oxygen keeps
  byte for byte as climate_targets_as_run.json until the lakes rerun on the resort's 301 K design case. The merged
  targets' lapse rates over the GCM's own land mask would cool these lakes' air by 0.026 K at most and 0.015 K
  typically, and their lock's share of 0 for −0.002 is the design case the lakes' lock already is.

## The checks

[test_run.py](test_run.py) checks the product's binding to its code and inputs; the solubilities against standard
values (CH4 1.55 mmol/L/atm and N2 0.70 at 20 °C, CO2 33.9 at 25 °C, and Weiss's O2 against the lake model's Benson and
Krause to 1%); the ebullition share against Lake Soppen, and a lunar floor under the Moon's own air bubbling as Earth's
does at 5.5–7 times the depth; bubbles losing more the deeper they start and the larger keeping more; the carbonate
system's round trip; an isolated hypolimnion's bubble point coming sooner on the Moon at every depth; the deep lakes
below their bubble point; methane rising with the trophic state, scaling with the area and lying within its brackets at
each of the wind's shares by each release rule, the stock's at least the net rise's; uncut lakes holding more methane
than cut ones; the dissolved load's phosphorus under each settling, the internal load as Nürnberg's release over the
anoxic days, added outside the retention, the warmer targets' richer floors anoxic for longer, and the steady state's
ceiling; each overturn releasing within the two rules, neither more than the stock there is or the methane made since
the last; the oxycline's share at Kivu's and the giants' figure lowering the flux; the oxygenation's bounds at the
demand the supply fixes, the warm transient beside them; and the three shares standing side by side.

## Running it

```sh
python -m research.studies.lake_anoxia.profiles      # fifteen minutes on two cores; needs lake_oxygen's cached runs
python -m research.studies.lake_anoxia.run           # about five minutes; reads LOLA through lake_oxygen's lakes
python -m research.studies.lake_anoxia.run --cache-only   # the same, stopping unless the saved columns start from
                                                          # the runs lake_oxygen's cache holds under its current keys
python -m pytest research/studies/lake_anoxia
```

profiles.py writes research/runs/lake_anoxia (ignored, on the research drive), with the sha256 of each saved column
and the lake model's in profiles.json, which run.py checks.
