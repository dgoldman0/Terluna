# What the lakes' anoxic depths do, and what design can do about them

17 October 2026. [The waters' oxygen](../lake_oxygen/README.md) found the Moon's deep crater lakes turning over only in
the rarest cold spells, their deepest tenth anoxic in mesotrophic and richer water in all three targets (the hot world
near 307 K, the resort at 300 K and the lock at 295 K), and overturns bringing the anoxic water up. The author asked
what that anoxic water does and what design can do about it: the lakes' methane and what it does to the air and its
dimming; whether CO2 and methane gather towards the hydrostatic bubble point behind Lake Nyos's eruption of 1986 and
the gas held in Lake Kivu; phosphorus released from anoxic floors; and the remedies, sized. [profiles.py](profiles.py)
reruns the lake model (`biosphere/ecology/lake.py`, unchanged) from the oxygen study's own settled columns and saves
each column daily; [run.py](run.py), with the gas physics in [gases.py](gases.py), writes
[results/lake_anoxia.json](results/lake_anoxia.json) (schema `terluna.research.lake-anoxia/1`); sources and how far
each was read are in [sources.json](sources.json). The seas' own basin run is not done.

## What it finds

**The lakes' methane depends on their trophic state far more than on the target, and only oligotrophic lakes fit the
air's held methane.** Methane is the share of the anoxic water's and floors' mineralization that respiration cannot
meet with oxygen; the lake model carries it as its debt of reduced substances, and in water as poor in sulphate as the
Moon's lakes (about 1.4 µM) methanogenesis takes about half of it (central 0.5, bracket 0.3–0.8), as in ferruginous
Lake Matano. Moon-wide, over geography's 3.92 million km² of rain-fed lakes, in Tg of methane a year:

| Central (bracket) | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic | 2.3 (0.5–11) | 4.5 (0.7–23) | 13 (1.0–50) |
| Mesotrophic | 20 (3.5–70) | 25 (3.5–77) | 53 (3.5–84) |
| Eutrophic | 69 (6.3–127) | 73 (6.4–142) | 110 (6.6–173) |
| New reservoir | 124 (11–204) | 125 (11–223) | 185 (11–288) |

The domain synthesis's ledger carries 39 Tg for the lakes (7.8–118), Earth's 10 g m⁻² a year over their area. The
Moon's lakes emit 0.6–3.3 g m⁻² a year oligotrophic and 18–28 g eutrophic. Their hypoxic depths make 0.2–1.0 mol of
methane per m² a year oligotrophic, 2.5–3.2 mol eutrophic and 4.0–4.9 in new reservoirs, the range of what Lake
Tanganyika's floors give its anoxic water (0.5 mol m⁻² a year; Durisch-Kaiser et al. 2011). The hot world raises the
oligotrophic lakes' emission five- to sixfold over the lock's, because their floors lie under anoxic water for 340 days
of the year against 107.

**Most of it rises with the overturns, which bring up far more methane than their oxygen can oxidise.** Between
overturns the anoxic water stores 0.2–0.8 mol per m² oligotrophic in the lock and the resort, 3–4 in the hot world, and
5–18 eutrophic, with 35–100 µM at the floor oligotrophic and 130–390 µM eutrophic and in new reservoirs, Tanganyika's
range (170–220 µM). When it turns over, an eutrophic lake or new reservoir holds 0–5 g of oxygen per m² against the
80–1,480 g its stored methane would take; an oligotrophic lake in the lock holds 240–430 g against 0–14. Earth's lakes
oxidise 54–98% of their stored methane at the autumn overturn, while oxygen enters as the mixed layer deepens over
weeks; the Moon's deep lakes turn over in one cold spell into water with almost no oxygen. The central flux lets half of
what the overturned lake's oxygen cannot meet escape, and 15% of the rest, the share Earth's lakes let through; the
bracket runs from 98% oxidised at the low methane share to 50% at the high.

**Bubbles reach six times deeper on the Moon but carry little of the lakes' methane.** A floor bubbles where its
porewater methane reaches the local pressure less its nitrogen; at a sixth of Earth's gravity the same pressure lies
six times deeper, so the floors bubble down to 50–60 m where Earth's stop near 10 m. The steady porewater balance
behind this (Bazhin 2010, as Langenegger et al. 2019 apply it) gives 11% of Lake Soppen's production as bubbles at
10 m, against the 10–17% measured. Bubbles rise more slowly at lunar gravity (a 5 mm bubble at 0.15 m/s against 0.23)
and lose more of their methane on the way: from 30 m a 10 mm bubble keeps 87% with a contaminated surface and 44%
clean, and a 5 mm one 67% and 5%. But 94% of the lake area lies in lakes deeper than 100 m, and only 4.2% of the floor
is shallower than 10 m, so bubbles carry 2–16% of the lakes' methane, where on Earth, whose lakes are small and
shallow, they carry about half (Bastviken et al. 2011).

**What it does to the air and the dimming.** The air holds no OH, so the soils' 2.92 Tg a year for each ppm
(1.07–4.77) is its only sink. With the lakes in the ledger at these rates and every other source cut 90% by biology, as
the held pair assumes:

| Steady methane, ppm; removal that holds 1.87 ppm, Tg a year (share of the air passed a year); added dimming if left | The lock | The resort | The hot world |
|---|---|---|---|
| Oligotrophic lakes left as they are | 1.6–2.4; 0–1.4 (0–0.04%); up to 0.10% | 2.4–3.1; 1.4–3.6 (0.04–0.11%); 0.10–0.21% | 5.3–6.0; 10–12 (0.30–0.37%); 0.49–0.57% |
| Mesotrophic, left | 7.5–8.3; 16–19 (0.5–0.6%); 0.74–0.79% | 9.3–10; 22–24 (0.7%); 0.84–0.89% | 19–20; 50–52 (1.5–1.6%); 1.3% |
| Eutrophic, left | 25; 66–68 (2.0–2.1%); 1.56–1.58% | 26–27; 70–72 (2.1–2.2%); 1.55–1.57% | 39; 107–110 (3.3%); 1.79–1.80% |
| Eutrophic, its methane cut 90% too | 3.2–3.9; 3.8–6.0; 0.23–0.35% | 3.3–4.1; 4.2–6.4; 0.24–0.36% | 4.6–5.3; 8.0–10; 0.41–0.50% |

The pairs are the domain synthesis's low- and high-cost ends; the dimming is the forcing beyond the held pair's
3.58 W/m² in the deep column, with N2O at 0.1 ppm, over each target's absorbed sunlight (211–225 W/m²), and comes on
top of the 1.69% the pair already takes. Oligotrophic lakes cut 90% hold the air at 0.9–2.0 ppm, under the held level.
Eutrophic lakes left alone would put the air at 25–39 ppm, past the domain synthesis's central steady state.

**No lake in the targets approaches a limnic eruption, but lunar gravity leaves little headroom.** The gas saturation
ratio is the dissolved gases' partial pressures with the vapour over the local pressure; bubbles can form where it
reaches 1. Lake Kivu's deep water reaches about 0.5 at 320 m (Bärenbold et al. 2020). In the 352 m lakes, between
overturns, the floor's water reaches 0.16–0.20 at its 6.8 atm, its gases at 1.04–1.34 atm, the air's nitrogen and
argon about 0.95 atm of it and methane most of the rest; in the
hypolimnion's top at 20 m, where the water still carries the air's nitrogen at 1.5 atm, it reaches 0.70–0.86, from
the air's gases more than the lakes' own. The 1,925 m lake left unmixed for fifty years goes anoxic in its deepest
tenth (losing 0.09–0.18 g m⁻³ a year) and would take 28,000–48,000 years to reach its 32 atm bubble point. An isolated
hypolimnion reaches its bubble point four to six times sooner than on Earth at the same depth, because its pressure
grows six times more slowly:

| Years for an isolated hypolimnion to reach its bubble point at the floor: lock, hot world (Earth at 25 °C) | Oligotrophic | Eutrophic |
|---|---|---|
| Floor at 20 m, 5 m of hypolimnion | 4.4, 3.5 (15) | 1.6, 1.3 (5.3) |
| 30 m, 15 m | 17, 14 (65) | 6.0, 4.8 (23) |
| 50 m, 35 m | 56, 46 (243) | 20, 16 (86) |
| 100 m, 50 m | 142, 115 (679) | 50, 41 (241) |
| 352 m, 50 m | 451, 367 (2,350) | 160, 130 (835) |

So the hazard lies where a lake tens of metres deep never turns over: a meromictic lake held by a salt or soda layer,
as closed basins become soda lakes over millennia (sea chemistry), a lake fed at depth by mineral springs as Kivu is,
or a deep new reservoir in its eutrophic first decade. There a eutrophic hypolimnion at 20–50 m reaches its bubble
point within one to twenty years. Its gas would be mostly methane, which rises and burns where Nyos's CO2 pooled.
Nyos's magma fed it 80 mol of CO2 per m² a year (Kling et al. 2005); a eutrophic lunar lake's depths make 2.5–3.2 mol
of methane and as much CO2.

**Phosphorus: anoxic floors would turn the lakes eutrophic once they pass mesotrophic.** On the land's dissolved
phosphorus alone, 0.018–0.048 g m⁻² a year, Vollenweider's relation puts the 352 m lakes at 1.1–3.0 mg/m³ of total
phosphorus, ultra-oligotrophic: the load is a ninth to a third of the 0.16 g m⁻² a year that would bring them to the
oligotrophic limit of 10 mg/m³; counting all the particulate load, 0.36–1.0 g, puts them at 24–67 mg/m³. Nürnberg's
anoxic release, 6.1 mg m⁻² a day mesotrophic and 10.8 eutrophic, over floors anoxic 330–350 days of the year, gives
0.2–3.8 g m⁻² a year (new reservoirs up to 9 g), 4–210 times the dissolved load, and 14–240 mg/m³: eutrophic to
hypertrophic. The Moon's low sulphate leaves iron to hold phosphate as vivianite (Caraco et al. 1989), which the bracket
takes as a tenth of Earth's release; oligotrophic floors release almost none whatever their oxygen (Nürnberg 1997). The
trophic state is a tipping point: the lakes stay oligotrophic on their dissolved load, and anything that makes them
mesotrophic (the particulate phosphorus becoming available, people's losses, a new reservoir's flush) starts a release
that carries them further.

## The remedies, sized

| Per km² of lake, and for all 3.92 million km² | The lock | The resort | The hot world |
|---|---|---|---|
| Oxygenating the depths: oxygen, t a day; up to, stirred and at warm kinetics | 0.32–1.31; up to 0.92–3.7 | 0.32–1.31; up to 1.4–5.7 | 0.32–1.31; up to 2.2–9.0 |
| Its power, kW (all lakes, GW) | 3.5–14 (14–57) to 13–55 (53–217) | 3.5–14 to 20–84 (80–329) | 3.6–15 to 33–134 (128–525) |
| Circulation by bubble plumes from 352 m; from 1,925 m, kW (GW) | 54 (212); 102 (400) | the same | the same |
| Cold inflows routed to depth: oxygen, share of the deep demand | 0.045 g m⁻² a day, 3–14% | 0.041 g, 3–13% | 0.037 g, 3–12% |
| Degassing pipes: water, m³/s; pumping against the methane's energy, kW | 0.2–0.7; 13–47 against 6–110 | 0.3–0.7; 23–48 against 12–120 | 0.4–0.5; 29–38 against 28–140 |

- **Nutrient targets.** The deep lakes stay oligotrophic under 0.16 g of bioavailable phosphorus per m² a year (10 mg/m³
  at their 1.8 m a year of throughflow), and ultra-oligotrophic under 0.063. Their dissolved load is a ninth to a third
  of it; the particulate load, kept in the deltas, must stay unavailable, which dredging the deltas for people's
  phosphorus (resources' plan, 4–51% of what the lakes receive) helps. Oligotrophic lakes are also the only ones whose
  methane stays near the held air: left uncut they need 0–1.4 Tg a year of removal in the lock, 1.4–3.6 in the resort
  and 10–12 in the hot world.
- **Oxygenating the depths.** The oxygen is the deep demand by trophic state (0.32–1.31 t per km² a day), raised by
  respiration's warm kinetics at each target's deep water (0.95 times in the lock, 1.44 in the resort, 2.3 in the hot
  world) and by the stirring the devices cause (1–3 times; Newman Lake's design took three). Earth's installations
  deliver 0.29–3.6 t per km² a day (Moore et al. 2015). Oxygen by vacuum swing adsorption from the air costs 0.20–0.24
  kWh a kilogram, and pressing it down to 352 m or 1,925 m at lunar gravity 0.06–0.12 more, half to two thirds of
  Earth's at the same depth. Line diffusers suit the deep lakes: their bubbles dissolve over the long lunar rise.
  Oxygenation stops the methane and the phosphorus release where it keeps the floor oxic; Lakes Baldegg and Sempach kept
  releasing phosphorus after ten years of it, so it backs up the nutrient targets (Gächter & Wehrli 1998).
- **Circulation.** Lorenzen & Fast's 9.2 m³ of air a minute per km² from 352 m costs 54 kW per km² on the Moon against
  92 on Earth, and from 1,925 m 102 against 136, at 60% compressor efficiency. A plume's rise through stratification
  hardly changes with gravity. The hot world's columns hold 40% more Schmidt stability than the lock's (930 against
  660 J/m²) and warm their deep water, which raises the demand.
- **Routing cold inflows.** The lakes' throughflow, 1.8 m a year, carries 3–14% of the deep demand's oxygen if all of
  it is taken to depth; the work against its buoyancy is a few watts per km², so it is cheap help for oligotrophic
  lakes and little for productive ones.
- **Lake depth by design.** A lunar night mixes a lake to a median 12–15 m and one night in ten to 25–140 m, so floors
  within about 12 m are ventilated through most nights. A reservoir or park lake shallower than that keeps an oxic
  floor with no running energy, and one 20–50 m deep must not be left to stratify for years: it reaches its bubble
  point in 1–20 years if eutrophic. Making the deep crater lakes shallow is out of reach: 330 m of fill under each km².
- **Degassing.** Pipes as at Nyos lift deep water through a degasser. The lunar lakes' deep water holds 1.04–1.34 atm
  of gas, so it lifts itself only in the hot world's new reservoir; elsewhere it is pumped, 13–48 kW per km², while
  the methane it carries holds 6–140 kW of fuel, Kivu's way. Methane burnt to CO2 warms the air far less. Degassing
  keeps meromictic and spring-fed lakes from their bubble point; in the overturning deep lakes oxygenation does more.

## How it is set up

- **The columns.** The nine 352 m lakes (three bands, three targets) start from the oxygen study's settled columns and
  run their ten years again, saved daily; the 1,925 m lake runs five passes, fifty years, from its oxygenated start.
  Each trophic state is the central column, surface-renewal transfer with demand fixed by the trophic state; the warm
  kinetics give the oxygenation's upper bracket.
- **Methane and CO2 from the oxygen.** Methane is the model's debt times the methane share, over 64 g of oxygen a
  mole; CO2 is the water's oxygen deficit below the euphotic zone and its debt's carbon (a mole of CO2 to every mole of
  methane, and one to each 32 g of the rest), added to lake water in equilibrium with the design air at the sea
  chemistry's lake pH, 7.72, an alkalinity of 0.42 meq/kg. N2 and Ar keep the surface's saturation. The anoxic methane
  is the model's demand in each layer times the share respiration cannot meet, K/(K + O2) with its 0.5 g/m³, in water
  under 2 mg/L.
- **Ebullition.** Production falls off exponentially in the sediment (30 per m) at the trophic state's floor demand
  (and the deep water's, within 25 m of the euphotic zone), the methane share of it; bubbles 80% methane, 10 mm across
  (5–20 mm, since weaker buoyancy lets bubbles grow larger before they break free of the mud, a guess on Earth's 2–8
  mm), rise with Tomiyama's drag for contaminated bubbles and lose gas by Frössling's transfer (contaminated) or
  Higbie's (clean). Each floor takes the mixed layer's temperature above 15 m and the deep water's below.
- **Scale-up.** Geography's lakes on LOLA give the floor by local depth (1.1% under 2 m, 4.2% under 10 m, 23% under
  100 m, 49% under 400 m); lakes over 1,000 m deep hold 70% of the area and those over 100 m 94%. Every lake takes the
  352 m lakes' pelagic flux per area, with each floor's own bubbles.
- **Phosphorus.** Vollenweider's TP = L/(q_s(1 + √τ)), with the lakes' 793 m mean depth and 445 years of residence
  (1.8 m a year); the 352 m lakes take that throughflow over their 112 m mean depth.

## Assumptions and limits

- The lake model carries no methane, CO2 or phosphorus: they follow from its oxygen and debt offline, so methane
  leaves the water only by the bracketed shares, and its own density (it lightens the water) and CO2's are not fed back.
- Every lake takes the 352 m lakes' pelagic flux. The giants over 1,000 m, 70% of the area, overturn more rarely and
  store more for longer behind a lasting oxycline that oxidises most of what crosses it, which would lower their flux;
  the 1,925 m runs show their depths going anoxic within fifty years.
- The overturn release, the oxidation shares and the methane share are Earth's ranges carried to the Moon. Bubble size
  at lunar gravity is a guess. Earth's anoxic release rates come from sulphate-bearing lakes.
- The isolated-hypolimnion years take the trophic state's demand over up to 50 m and no exchange; real lakes that
  never turn over still lose gas upward by diffusion, as Tanganyika and Kivu do.
- The seas' basin run, the remedies run in the model and the biology of methanotrophs designed for sudden overturns
  are not done.

## The checks

[test_run.py](test_run.py) checks the product's binding to its code and inputs; the solubilities against standard
values (CH4 1.55 mmol/L/atm and N2 0.70 at 20 °C, CO2 33.9 at 25 °C, and Weiss's O2 against the lake model's Benson and
Krause to 1%); the ebullition share against Lake Soppen and its exact gravity scaling; bubbles losing more the deeper
they start and the larger keeping more; the carbonate system's round trip; an isolated hypolimnion's bubble point
coming sooner on the Moon at every depth; the deep lakes below their bubble point; methane rising with the trophic
state and scaling with the area; and the dissolved load keeping the lakes oligotrophic.

## Running it

```sh
python -m research.studies.lake_anoxia.profiles      # about five minutes; needs lake_oxygen's cached runs
python -m research.studies.lake_anoxia.run           # about a minute and a half; reads LOLA through lake_oxygen's lakes
python -m pytest research/studies/lake_anoxia
```

profiles.py writes research/runs/lake_anoxia (ignored, on the research drive), with the sha256 of each saved column
in profiles.json, which run.py checks.
