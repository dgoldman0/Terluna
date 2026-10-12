# Food systems for the Open Moon

This study prices eight ways of growing the lunar design's plant-rich diet for each billion residents, fits each to the
29.5-day light cycle, and combines them into seven mixes for 4, 6 and 8 billion people. It reads ecology's diets,
crop-model yields and nutrient targets and resources' phosphorus accounts for the sky and the waters (from their
branches, read in the tree at pinned hashes since the joint integration), the canopy, climate, drainage, aerosol, heat
and power products, and published measurements for every food system (`sources.json`, 100 sources).
Numbers computed here are marked *screen*; arithmetic on cited values is marked *derived*. The product is
`results/food_systems.json` (schema `terluna.research.food-systems/1`).

## Findings

1. **Rain-fed fields feed the plant-rich diet at the least cost in every currency, and the rain alone feeds 9.1–27.0
   billion people outside a kept 30%, or 6.5–19.3 billion outside a kept half** (*screen*). A billion residents need
   141,000–401,000 km² of crops at the crop model's equatorial yields, 248,000–731,000 km² on the wettest band: 2.6–7.7%
   of what the rain can grow in the bands where most land is wet, 5.22–5.45 million km² of equivalent once the crops'
   water need is held to what the design sunlight can evaporate (Departures). Their only power is nitrogen fixation,
   6.3–14 W a resident. They lose 0.044–0.18 Mt of phosphorus a year and emit 17–65 kt N of N2O and 0.3–3 Mt of methane
   at the lunar design's targets. Farms that reach the 70–80% of potential yield that Earth's best farms reach (Cassman
   1999; Lobell et al. 2009) need 25–43% more land, and the rain then feeds 6.4–22 billion outside a kept 30%.

2. **Water from the seas grows what the rain cannot for 0.92–1.44 kWh a cubic metre** (*screen*). Lifting it to the
   fields' mean height (1.4–1.5 km at a sixth of Earth's gravity), piping it 69–95 km and pressing it through the
   emitters costs about half of what California's State Water Project spends (2.4 kWh/m³; Cohen et al. 2004). Topping up
   a billion residents' rain-belt fields takes 115–540 km³ a year, 1.6–7.5% of the rivers' flow, and 13–89 W a resident.
   Watered fields in the fog desert take 247,000–751,000 km², 263–1,180 km³ and 28–176 W a resident; the sunlight at
   45–60° lets them transpire enough for 73–78% of the crop model's production there.

3. **Sunlit glass needs 1.2–2.5 times the equatorial field area, 4–23 Gt of glass and steel a billion and 0.16–1.2 kW a
   resident** (*screen*). The glass passes 66–73% of the light (Li et al. 2014; Kempkes et al. 2022), so a house yields
   66–89% of an open field at the crop model's potential, CO2 enrichment giving back 0–10% (*derived*, Nederhoff 1994),
   and the fog desert's sun at 47° gives 72% of the equator's yield. The sunlight absorbed under the glass lets the
   crops transpire 2.0–2.3 mm a day there, so a fog-desert house grows 85–90% of what its light allows. Vented houses
   lose the crops' transpiration, 204–736 km³ a year a billion: in the fog desert the rain on the glass covers 7% of it
   and the seas supply 189–685 km³; on the highlands the rain on the glass covers it all. Closing a house saves that
   water and costs 2.4–13 kW a resident to pump out the absorbed sunlight. Renewing the glass every 15 years takes
   133–788 Mt a year a billion. Fields at Earth's 70–80% of potential yield about as much a square metre as glass at the
   crop model's potential.

4. **Fog nets catch little in the Moon's fogs, which form in nearly still air, so the fog desert's water comes by pipe
   from the seas** (*screen*). Calm radiation fogs make poor nets on Earth too: collection stops near 1 m/s (Montecinos
   et al. 2018), Raschel nets work in 2–6.5 m/s winds and catch most of their water above 4 m/s (Rivera 2011), and thin
   radiation fog holds a median 0.010 g/m³ of liquid water (Mazoyer et al. 2022). The design run's fogs on the circles
   of latitude at 45° N and 45° S, the CM1 rings in the fog desert's latitudes, come in the night's calmest hours, in
   hour-angle bins whose mean wind in fog is 0.004–0.12 m/s against whole-cycle medians of 0.9–1.1 m/s. The ring
   products give fog and mean wind only for each 10° bin of hour angle, over every land column and snapshot and with the
   winds to 0.1 m/s, so the bracket spreads each bin's winds as a Weibull distribution about its mean, from a Rayleigh
   spread with fog and wind independent to each ring's whole-cycle spread with the fog in the windiest share of the bin.
   Nets then catch 0–98 L/m² a year above the 1 m/s threshold (0.007–145 L/m² counting all wind), against 1,100–3,650
   L/m² at Earth's working fog-collection sites (Klemm et al. 2012); the two rings tilted 70° across the equator give
   0.004–68 L/m². At the bracket's top, nets for a billion residents' fog-desert glass would cover 1.9 million km²,
   two-fifths of the fog desert. Rain on glass gives 51–61 L/m² a year. Water drawn from the air costs 238–500 kWh/m³
   (Peeters et al. 2021; Lord et al. 2021), 5.1–39 TW a billion for fog-desert glass, against 0.019–0.09 TW to lift and
   pipe it from the seas at 0.89–1.15 kWh/m³.

5. **Lamp-lit staples take 5.6–61 kW a resident, and the fully lamp-lit diet for 8 billion takes 41.5–425 TW, whose heat
   needs 0.50–8.1% more dimming** (*screen*). Under LEDs a kilowatt-hour grows 7.4–39 kcal of potato or 3.7–16 kcal of
   wheat, and the varied diet with its losses needs 2.1–4.3 times its 2,500 kcal grown. Lit by ground solar with night
   storage, lamp farms need 191,000–2,090,000 km² of panels a billion, more than the 141,000–401,000 km² of fields they
   replace; lit by the array, 114,000–1,240,000 km² of receivers. The array intercepts 1,900–2,400 PW; the bound is
   heat. Each 1% of added dimming answers 72–115 TW of heat on the Moon (population study), and the dimming takes the
   same share of sunlight from every wild plant. Fusion's heat raises the dimming to 0.7–18% at 8 billion.

6. **Microbial protein is the cheapest food made on power: 40–133 W a resident replaces the diet's animal part**
   (*screen*). Hydrogen-oxidising bacteria take 30–87 kWh per kg of protein across three process accounts (*derived*;
   Järviö et al. 2021, Science of the Total Environment; García Martínez et al. 2021; Leger et al. 2021), grown with the
   diet's losses of a quarter to a third. Their nucleic acids hold processed biomass to 50–80 g of protein a day
   (*derived*; Ritala et al. 2017; Coelho et al. 2020), so half the reference diet's protein, 47 g, takes 77–257 W. A
   whole diet made on power, with the other half of its protein fermented by a fungus on sugar made from power and
   electrochemical sugar for the food energy the protein and produce leave, takes 0.70–2.8 kW a resident.

7. **Aquaponic fish take four to six times the cropland of the animal part they replace and emit three to seven times
   the whole field diet's N2O** (*screen*). Replacing the diet's 24.3 g of animal protein a day, with its losses, takes
   128–190 Mt of live fish a billion, fed on 141–425 Mt of feed grown on 45,100–315,000 km² of crops; ecology's animal
   part takes 10,000–53,000 km² of crops and up to a million km² of grazing. The fish emit 52–446 kt N of N2O against
   the field diet's 17–65, and send 0.28–4.8 Mt of phosphorus to their sludge, which must go back to the fields.

8. **The waters' phosphorus holds sea farms to the animal part of 0.04–0.26 billion people, so the sea-supported mix
   cannot be built at 4–8 billion** (*screen*). Per billion, seaweed would take 20,500–83,100 km² of fresh sea and
   shellfish 39,900–407,000 km² for the whole animal part; half of each draws 0.22–0.51 Mt of phosphorus a year, three
   to five times what the people lose. The rain-fed lakes keep 96% of the land's phosphorus, the common-level waters
   receive 4.6–12 kt a year of it dissolved, and the waters' life can spare 21–56 kt a year in all to seaweed, fish and
   sorption, 1.4–3.7 kt of it from the common-level waters (resources' sea-phosphorus study). At 4–8 billion the
   sea-supported mix draws 0.86–4.1 Mt a year, 15–195 times that; if particles release phosphorus as in Earth's
   estuaries (44–377 kt a year), 2.3–94 times. At the high-cost end the farms would also cover 2.9–5.9 times the seas
   shallower than 50 m. Fed fish pens in a Chinese lake emit 40 g of methane a square metre a year above open water (Pu
   et al. 2023); taken as the farms' upper bound, it brings the sea-supported mix's methane to 102 Mt a year at 8
   billion.

9. **N2O accumulates in every mix, and glass and fish multiply it** (*screen*). N2O has no sink in the design air, so
   every 3.07 Tg N stays as 1 ppb (ecology's first screen), and the wild biosphere alone adds 0.31–0.42 ppb a year.
   Fields at the lunar design's 0.2–0.5% emission add 0.045–0.17 ppb a year at 8 billion. Soilless culture emits
   0.72–1.19% of its nitrogen (Llorach-Massana et al. 2017; Daum & Schenk 1996) and aquaponic fish 0.8–1.9% of their
   feed's (Hu et al. 2015; Yuan et al. 2019). The fish feed's crops take their own nitrogen, so the fog-desert glass mix
   adds 0.50–3.4 ppb a year at 8 billion (0.50–0.71 in its feasible range), 500–710 ppb in a thousand years.

10. **At 4–8 billion, pasture cannot carry the diet's animal part** (*screen*). Ecology's design allows 0–1,000 m² of
    grazing a person on the fog desert. At the upper allowance the fog desert outside a kept 30% holds pasture for 3.4
    billion people, and 2.4 billion outside a kept half, so the rain-fed-first mix's high-cost end takes 119–237% of the
    fog desert outside a kept 30% at 4–8 billion and cannot be built; its feasible range ends 42–84% of the way to that
    end. At 4–8 billion the diet's animal part comes from microbial protein, aquaponics, sea farms or animals fed on
    crops.

11. **At 8 billion, the feasible cases of every buildable mix with fields leave 35–89% of the dry land wild outside a
    kept 30%; lamps and power-made food leave 97% or more** (*screen*). Rain-fed first, with grazing held to the fog
    desert, uses 2.08–3.8 million km², 0.05–0.077 TW and no sea water. Greenhouse-heavy fits to 95% of the way to its
    high-cost end, where its glass fills the highlands outside a kept 30%. Fog-desert glass spares the rain belt whole
    and fits the fog desert only to 16% of the way (2.6–3.5 million km², 2.8–3.9 TW, 1,970–2,630 km³ of sea water a
    year); outside a kept half it fits at neither end. Power-made food takes 5.6–22.4 TW (0.067–0.43% more dimming), the
    fully lamp-lit diet 41.5–425 TW, against 12–36 TW for the residents' other energy.

12. **Under the author's wild share, at least half wild or 30% wild with another 20% hybrid, the power-heavy mix fits
    at 4–8 billion at both ends** (*screen*). The hybrid fifth is the author's megaparks and great gardens, closed to
    fields, glass, towns and grazing, so the compromise closes the same half of each zone as a kept half and every mix
    matches the half's results: the rain alone feeds 6.5–19.3 billion. Half wild with another 20% hybrid, run beside
    them for comparison at the author's request, leaves fields, glass, towns and grazing 30% of each zone: the rain
    alone feeds 3.9–11.6 billion, the rain-fed-first mix fits 18–36% of the way to its high-cost end, and at 8 billion
    the power-heavy mix's high-cost end grows 1.12 million km² of its crops' equivalent on watered fields, with 3,430
    km³ a year from the seas and 0.56 TW of lift.

13. **The author's blend, power-heavy 40%, power-made 50% and fully lamp-lit 10%, takes the least sunlight from wild
    life at the high-cost end within the limits set for it** (*screen*; [blend.py](blend.py),
    [results/food_blend.json](results/food_blend.json)). With each mix between a tenth and a half of the diet,
    power-made takes the cap and fully lamp-lit the floor, so the limits set the weights. At 8 billion under the
    compromise its fields take 0.68–1.86 million km², none of them in the middle latitudes, and its heat is 15–171 TW
    with fusion, 0.13–2.4% more dimming. Heritage meat at a kilogram a person a year takes 0.30 million km², a fifth of
    the homestead region, and its cattle add 1.6–1.7 Mt of methane a year.

## The systems, per billion residents

Each row grows one part of the diet for a billion residents: the whole plant-rich diet for fields, glass and lamp farms,
and the diet's 24.3 g of animal-source protein a day for the protein systems (*screen*). TW per billion equals kW a
resident.

| System | Land or sea, thousand km² | Water, km³ a year | Power a resident | Phosphorus, Mt a year | N2O, kt N a year | Structure |
|---|---|---|---|---|---|---|
| Rain-fed fields | 141–401 eq; 248–731 on the wettest band | 247–736 rain | 6.25–14.4 W | 0.044–0.18 | 17.3–65 | – |
| Irrigated top-up, rain belt | 144–409 | 115–540 from the seas | 13.2–88.8 W | as fields | as fields | – |
| Irrigated fields, fog desert | 247–751 | 263–1,180 from the seas | 27.6–176 W | as fields | as fields | – |
| Vented glass, fog desert | 246–993 | 189–685 from the seas | 244–1,220 W | 0–0.09 discharged | 62.4–155 | 2.93–11.8 Gt glass, 2.7–10.9 Gt steel |
| Closed glass, fog desert | 207–841 | 0 from the seas | 2,410–13,000 W | 0–0.09 discharged | 62.4–155 | 2.46–10 Gt glass, 2.28–9.24 Gt steel |
| Vented glass, highlands | 177–709 | 0 from the seas | 162–809 W | 0–0.09 discharged | 62.4–155 | 2.11–8.44 Gt glass, 1.95–7.79 Gt steel |
| Vented glass, rain belt | 167–658 | 21.8–133 from the seas | 155–770 W | 0–0.09 discharged | 62.4–155 | 1.99–7.82 Gt glass, 1.84–7.22 Gt steel |
| Vertical farms, potato | 5.01–64.7 floor | recycled | 5.56–60.9 kW | – | – | 2.5–32.3 Gt building |
| Vertical farms, wheat | 2.83–28.7 floor | recycled | 14–120 kW | – | – | 1.42–14.3 Gt building |
| Lamp-lit vegetables and fruit | – | – | 90.4–704 W | – | – | – |
| Microbial protein, animal part | – | – | 39.7–133 W | 0.0913–0.295 in biomass | – | – |
| Microbial protein, half the protein | – | – | 76.9–257 W | 0.177–0.572 in biomass | – | – |
| Precision fermentation, animal part | 7.15–11.4 eq of sugar crops | – | 14.7–17 W (+46.8–167 W with electro-sugar) | – | – | – |
| Power-made food | – | – | 0.695–2.79 kW | 0.177–0.572 taken up | – | – |
| Aquaponics, animal part | 45.1–315 eq of feed crops | 9.72–14.4 | 45–66.8 W | 0.283–4.79 to sludge | 52.2–446 | – |
| Seaweed, animal part | 20.5–83.1 sea | – | – | 0.289–0.627 drawn from the seas | – | – |
| Shellfish, animal part | 39.9–407 sea | – | – | 0.144–0.399 drawn from the seas | – | – |

Every system feeds people who eat 0.44–0.60 kg of phosphorus a year and lose 0.044–0.18 kg of it at the lunar design's
recycling (ecology); the phosphorus column gives what a system adds or draws. Heat follows power: each TW used heats the
Moon by 1.38 TW with power from the array (its receivers' losses included) and by 2.0–3.0 TW with fusion on the Moon,
and each TW of heat on the Moon is answered by 0.0087–0.0139% more dimming (population study). The product gives heat,
W/m² and dimming for every block.

**Fields** use ecology's crop model: the day fruit at the model's potential, clouds (0.89), losses of a quarter to a
third, a varied diet and the plant-rich diet's crop factor (people_and_land.md, item 11). Rain supports a share of the
model's growth in each latitude band (ecology), raised here to the water need the sunlight allows (Departures): 56–58%
within 15° of the equator, 32–34% at 15–30° and 17–18% at 30–45°. The atlas bands then hold 5.54–5.78 million km² of
rain-fed equivalent, 5.22–5.45 million km² of it in the bands the fill reaches (0–45°, the 30–45° band held to its wet
lowland, taken as the wet middle latitudes, and its own highlands). Watering lifts water from the seas to each zone's
mean height, pipes it the zone's mean distance through a 2 m main and presses it to 100–400 kPa at the emitters (main
and pressure are design guesses), at the 60–90% application efficiency of surface to drip irrigation (Brouwer et al.
1989), with one pump efficiency at each end, 0.85 to 0.7. A watered field grows what the sunlight lets it transpire, at
most the crop model's production at its latitude. Nitrogen, phosphorus, N2O and methane are ecology's lunar design per
person.

**Glass** is a Venlo house of 11.9 kg of glass, 11.0 kg of steel and 2.8 kg of aluminium a square metre, lasting 15
years (the EUPHOROS inventory, as given by Kan & Vieira), using 8–10 kWh/m² a year of electricity unlit (Raaphorst et
al. 2002; Montero et al. 2010). Its yield is the crop model's at the zone's latitude times the house's light factor
(transmission 0.66–0.73 raised to the crops' light response, 0.7–1.0; Marcelis et al. 2006) and the CO2 gain, held to
what the sunlight absorbed inside lets the crops transpire. A vented house holds CO2 at 480–700 ppm Earth-equivalent and
a closed one at 700–1,000 ppm (design guesses); the CO2 cuts the crops' water per unit of growth through the gain and by
closing their stomata 3–4% per 100 ppm (Nederhoff 1994), 1.0–1.21 times in vented houses and 1.18–1.49 times in closed
ones. Closed houses pump the absorbed sunlight out at a coefficient of performance of 4–6 (Kozai 2013), with the same
transmission at each end as the yield. As a check, the crop model's equatorial yield, 11,000–15,300 kcal/m² a year,
overlaps the top of the food energy of unlit Dutch glasshouse tomatoes, 9,900–12,400 kcal/m² a year, and stays below a
lit hybrid crop's 16,700 (*derived*; Raaphorst et al. 2002; de Gelder & Grootscholten 2011; Helmus-Schuddebeurs et al.
2019). Soilless culture discharges 0–10% of its nutrient solution's phosphorus a year (Beerling et al. 2014).

**Lamp-lit farms** grow the diet's food energy at a staple's edible yield per mole of light under lamps: potato
0.26–0.82 g dry per mol and wheat 0.13 g in NASA's Biomass Production Chamber (Wheeler et al. 2003; Wheeler 2020), wheat
0.32 g in Asseng et al.'s (2020) vertical farm (*derived*), at LED efficacies of 2.5–3.7 µmol/J (Kusuma et al. 2020)
with lighting 85–98% of the load (Kozai 2013; Asseng et al. 2020). Vegetables come at 7.8–44 g fresh per mol (NASA;
Carotti et al. 2021) and fruit at 4.2–13.5 g (NASA; Helmus-Schuddebeurs et al. 2019). Shelves stand 5–10 layers to a
floor (Graamans et al. 2018; Asseng et al. 2020), and the floors weigh 500 kg/m² (the summit tower study).

**Power-made protein and food.** Bacteria on electrolytic hydrogen give biomass of 50–80% protein at 30–87 kWh per kg
protein (*derived*). Precision fermentation of a fungal protein takes 10.5 kWh of electricity and 2.5 kg of glucose a kg
of protein (Järviö et al. 2021, Nature Food). Sugar made from CO2 and power returns 10–31% of the electricity as food
energy (García Martínez et al. 2021, Journal of CO2 Utilization; Dinger & Platt 2020), at glucose's 3.72 kcal a gram
(2,805 kJ a mole; NIST Chemistry WebBook); Hann et al.'s (2022) acetate route to yeast and algae returns 11–16%
(*derived*). The power-made diet takes half its protein from bacteria, the most their nucleic acids allow, and half from
fermentation on sugar made from power; electrochemical sugar supplies the 2,010–2,020 kcal a day that its protein (at
4.27–4.36 kcal a gram, USDA FoodData Central) and its lamp-grown lettuce and tomatoes (75 kcal) leave. Every protein
system grows the diet's losses, 1.3–1.5 times what is eaten, as the plant foods do.

**Aquaponics** raises tilapia in recirculating tanks at 3.1 kWh and 76 L of water per kg live (Bergman et al. 2020),
with feed conversion of 1.1–2.24 (Bergman et al. 2020; Jatobá et al. 2026), fillets 35–45% of the fish (Bergman et al.
2020; Fry et al. 2018), and feed of 29–35% protein and 0.71–1.41% phosphorus (Syandri et al. 2018, 2021). The fish and
plants keep 20–72% of the feed's phosphorus (Cerozi & Fitzsimmons 2017; Jatobá et al. 2026). The feed's crops take the
feed's own nitrogen at the chain's 40–60% efficiency, with the fields' or the glass's N2O. The feed's phosphorus,
1.0–6.0 Mt a billion, cycles through the fields and the sludge; the mixes count the fields' losses for the diet alone.

**Sea farms** grow seaweed at 16–30 t dry a hectare a year, China's mean over all farmed species to elite kelp (Xiao et
al. 2017; Li et al. 2016), with 10–19% protein (Angell et al. 2016; Xiao et al. 2017), and mussels at 7–60 t live a
hectare, a zebra-mussel farm in a near-fresh lagoon to a Danish fjord's mussel lines (Schernewski et al. 2012; Petersen
et al. 2014), 40% of it meat (Magnani et al. 2026). Their harvest draws its phosphorus against what resources' study
lets the waters' life spare, 30% of each water's dissolved input, a design guess of that study.

## Yields

| System | Food per area | Food grown per kWh |
|---|---|---|
| Rain-fed fields, at the equator | 2,490–7,090 residents a km² | – |
| Vented glass, fog desert | 1,010–4,070 residents a km² | – |
| Lamp-lit potato-like staple | 15,500–200,000 residents a km² of floor | 7.4–39 kcal |
| Lamp-lit wheat | 34,800–353,000 residents a km² of floor | 3.7–16 kcal |
| Microbial protein | – | 11–33 g of protein |
| Sugar from power | – | 86–267 kcal |
| Seaweed | 16–30 t dry a hectare a year | – |
| Mussels | 7–60 t live a hectare a year | – |

## Fit to the lunar cycle

Fields need no night energy. Evergreen crops store carbon by day and idle on it at a quarter of their daytime upkeep
once their own balance turns negative: 10.9 g C/m² carries an idling crop through the night, 48.6 g an Earth-like one
(canopy plant model). The day fruit is set at sunrise and ripe at sunset, so the harvest comes once a cycle at dusk and
the towns hold a cycle's food (lunar-cycle ecology register). Light stays under 1 µmol m⁻² s⁻¹ for 254 h of the cycle at
the equator, 237 h at 30° and 120 h at 60°, where twilight lights most of the night (canopy plant model).

Sunlit glass runs as the fields beneath it. Lamps through the dark hours (170 h in the fog desert, 240–245 h on the
highlands and in the rain belt) add 29–34% to a house's light over the cycle and save 37,800–250,000 km² of glass a
billion for 2.7–16 TW (*screen*): each square metre of glass that lamps replace costs 46–95 W around the clock.

Vertical farms keep their own clock, and a resident's lamps draw 1,970–21,600 kWh through each 354-hour night, 3–33
times a European's night use (650 kWh, population study), from the array, fusion or storage. Microbial reactors take
14–47 kWh a resident a night. Aquaponic fish are fed on a 24-hour routine from stored feed. Seaweed idles through the
15-day night on its stores: Arctic kelp survived seven months of continuous darkness (Henley & Dunton 1997), and
Svalbard's sugar kelp keeps its photosynthetic capacity through four months of polar night on its laminaran store (Li et
al. 2020). Shellfish feed on plankton that idles with the kelp.

## Mixes for 4, 6 and 8 billion

Half the residents live aloft (2, 3 and 4 billion, within the author's 2–4 billion) and the rest in highland towns at
the lunar design's 25–100 m² a person. Fields fill the bands where most land is wet (0–45°), wettest first, outside the
kept share and the towns, in proportion to each band's highland and lowland dry land; what the rain cannot grow comes
from watering those fields, and beyond that from watered fields in the fog desert. Each case is given at the low-cost
and high-cost ends of every input (*screen*). Rows give the kept-30% case unless they name the half. The residents'
other energy heats the Moon by 6–18, 9–27 and 12–36 TW (efficient to high-income use with industry in orbit, population
study).

A case fits when every zone it uses, grazing included, lies within that zone's dry land outside the kept share, its sea
farms within the seas shallower than 200 m, and their phosphorus within what the waters' life can spare. Where the
high-cost end overflows, a cell gives the feasible range beside the full one: it runs from the low-cost end to the case
where the first zone fills, every input moved the same fraction of the way towards its high-cost value.

**Rain-fed first.** The plant-rich diet on fields: rain-fed first, wettest land first; the rest by watering those fields
from the seas, then watered fields in the fog desert; ecology's grazing on the fog desert.

| | 4 billion | 6 billion | 8 billion |
|---|---|---|---|
| Fits its zones, kept 30% | to 84% of the way (fog desert full) | to 56% of the way (fog desert full) | to 42% of the way (fog desert full) |
| Fits its zones, kept half | to 60% of the way (fog desert full) | to 40% of the way (fog desert full) | to 30% of the way (fog desert full) |
| Land used, million km² | 1.04–3.12; feasible 1.04–2.79 | 1.56–5.24; feasible 1.56–3.29 | 2.08–8.16; feasible 2.08–3.8 |
| &nbsp;&nbsp;rain belt | 0.719–2.15; feasible 0.719–1.92 | 1.08–3.79; feasible 1.08–2.28 | 1.45–5.93; feasible 1.45–2.64 |
| &nbsp;&nbsp;wet middle latitudes | 0 | 0 | 0–0.249; feasible 0 |
| &nbsp;&nbsp;highlands, towns included | 0.322–0.969; feasible 0.322–0.868 | 0.48–1.45; feasible 0.48–1.01 | 0.637–1.98; feasible 0.637–1.16 |
| Grazing on the fog desert, million km² | 0–4; feasible 0–3.37 | 0–6; feasible 0–3.37 | 0–8; feasible 0–3.37 |
| Dry land left wild outside a kept 30% | 94% to 57%; feasible 94% to 62% | 91% to 31%; feasible 91% to 59% | 87% to 1%; feasible 87% to 56% |
| Dry land left wild outside a kept half | 91% to 37%; feasible 91% to 60% | 87% to over by 8%; feasible 87% to 56% | 82% to over by 31%; feasible 82% to 48% |
| Rain belt left wild outside a kept 30% | 88% to 64%; feasible 88% to 68% | 82% to 36%; feasible 82% to 62% | 76% to none; feasible 76% to 56% |
| Fog desert used outside a kept 30%, grazing included | 0%–119%; feasible 0%–100% | 0%–178%; feasible 0%–100% | 0%–237%; feasible 0%–100% |
| Rain-fed land left wild outside a kept 30% | 85% to 56%; feasible 85% to 61% | 78% to 34%; feasible 78% to 54% | 70% to 12%; feasible 70% to 47% |
| Power, TW | 0.025–0.0576; feasible 0.025–0.0525 | 0.0375–0.0864; feasible 0.0375–0.065 | 0.05–0.115; feasible 0.05–0.0775 |
| Added dimming, power from the array | 0.000299–0.0011%; feasible 0.000299–0.000944% | 0.000448–0.00165%; feasible 0.000448–0.00104% | 0.000597–0.0022%; feasible 0.000597–0.00116% |
| Added dimming, fusion | 0.000434–0.00242%; feasible 0.000434–0.00197% | 0.000651–0.00364%; feasible 0.000651–0.00194% | 0.000868–0.00485%; feasible 0.000868–0.00205% |
| Water from the seas, km³ a year | 0 | 0 | 0 |
| Phosphorus lost, Mt a year | 0.176–0.72; feasible 0.176–0.634 | 0.264–1.08; feasible 0.264–0.722 | 0.352–1.44; feasible 0.352–0.81 |
| N2O, ppb a year | 0.0225–0.0847; feasible 0.0225–0.0749 | 0.0338–0.127; feasible 0.0338–0.0862 | 0.0451–0.169; feasible 0.0451–0.0975 |
| Methane, Mt a year (steady ppm) | 1.2–12 (0.252–11.2); feasible 1.2–10.3 (0.252–6.24) | 1.8–18 (0.377–16.8); feasible 1.8–10.9 (0.377–4.05) | 2.4–24 (0.503–22.4); feasible 2.4–11.5 (0.503–3.58) |

**Greenhouse-heavy.** Half the plant foods in vented sunlit glasshouses on the highlands beside the towns, half on
rain-fed fields; the animal part as aquaponic fish fed from the fields.

| | 4 billion | 6 billion | 8 billion |
|---|---|---|---|
| Fits its zones, kept 30% | both ends | both ends | to 95% of the way (highlands full) |
| Fits its zones, kept half | both ends | to 90% of the way (highlands full) | to 63% of the way (highlands full) |
| Land used, million km² | 1.16–5.25 | 1.74–10.4 | 2.31–10.6 |
| &nbsp;&nbsp;rain belt | 0.584–3.26 | 0.893–5.93 | 1.22–5.85; feasible 1.22–5.93 |
| &nbsp;&nbsp;wet middle latitudes | 0 | 0–1.8 | 0–1.92; feasible 0–1.95 |
| &nbsp;&nbsp;highlands, towns included | 0.573–1.99 | 0.842–2.64 | 1.1–2.86; feasible 1.1–2.73 |
| &nbsp;&nbsp;of it glass | 0.329–1.23 | 0.494–1.85 | 0.659–2.46; feasible 0.659–2.34 |
| Dry land left wild outside a kept 30% | 93% to 68% | 89% to 37% | 86% to 35% |
| Dry land left wild outside a kept half | 90% to 47% | 85% to 35% | 80% to 34%; feasible 80% to 35% |
| Rain belt left wild outside a kept 30% | 90% to 45% | 85% to none | 80% to 1%; feasible 80% to none |
| Highlands used outside a kept 30% | 21%–73% | 31%–97% | 40%–105%; feasible 40%–100% |
| Highlands used outside a kept half | 29%–92% | 42%–110%; feasible 42%–100% | 54%–147%; feasible 54%–100% |
| Rain-fed land left wild outside a kept 30% | 88% to 47% | 83% to 20% | 77% to 25%; feasible 77% to 24% |
| Power, TW | 0.506–1.73 | 0.759–2.59 | 1.01–4.06; feasible 1.01–3.71 |
| Added dimming, power from the array | 0.00604–0.0331% | 0.00906–0.0496% | 0.0121–0.0775%; feasible 0.0121–0.0696% |
| Added dimming, fusion | 0.00878–0.0728% | 0.0132–0.109% | 0.0176–0.171%; feasible 0.0176–0.151% |
| Water from the seas, km³ a year | 38.9–57.7 | 58.3–86.5 | 77.8–3,750; feasible 77.8–2,760 |
| Phosphorus lost, Mt a year | 0.176–0.72 | 0.264–1.08 | 0.352–1.44; feasible 0.352–1.39 |
| Phosphorus discharged from glass, Mt a year | 0–0.18 | 0–0.27 | 0–0.36; feasible 0–0.335 |
| Phosphorus to fish sludge, Mt a year | 1.13–19.2 | 1.7–28.8 | 2.26–38.4; feasible 2.26–34.6 |
| N2O, ppb a year | 0.148–1.1 | 0.221–1.65 | 0.295–2.2; feasible 0.295–2.02 |
| Methane, Mt a year (steady ppm) | 1.2–12 (0.252–11.2) | 1.8–18 (0.377–16.8) | 2.4–24 (0.503–22.4); feasible 2.4–22.9 (0.503–18.3) |

**Fog-desert glass.** All plant foods and the fish feed in vented glasshouses in the fog desert, watered by rain on the
glass and pipes from the seas; the animal part as aquaponic fish; the rain belt left to the wild.

| | 4 billion | 6 billion | 8 billion |
|---|---|---|---|
| Fits its zones, kept 30% | to 57% of the way (fog desert full) | to 32% of the way (fog desert full) | to 16% of the way (fog desert full) |
| Fits its zones, kept half | to 36% of the way (fog desert full) | to 13% of the way (fog desert full) | neither end (fog desert short) |
| Land used, million km² | 1.28–6.77; feasible 1.28–3.51 | 1.92–10.2; feasible 1.92–3.52 | 2.56–13.5; feasible 2.56–3.52 |
| &nbsp;&nbsp;fog desert | 1.23–6.57; feasible 1.23–3.37 | 1.84–9.85; feasible 1.84–3.37 | 2.46–13.1; feasible 2.46–3.37 |
| &nbsp;&nbsp;highlands, towns included | 0.05–0.2; feasible 0.05–0.135 | 0.075–0.3; feasible 0.075–0.146 | 0.1–0.4; feasible 0.1–0.147 |
| &nbsp;&nbsp;of it glass | 1.23–6.57; feasible 1.23–3.37 | 1.84–9.85; feasible 1.84–3.37 | 2.46–13.1; feasible 2.46–3.37 |
| Dry land left wild outside a kept 30% | 92% to 59%; feasible 92% to 79% | 88% to 38%; feasible 88% to 79% | 84% to 17%; feasible 84% to 79% |
| Dry land left wild outside a kept half | 89% to 42%; feasible 89% to 79% | 84% to 13%; feasible 84% to 79% | 78% to over by 16% |
| Rain belt left wild outside a kept 30% | 100% | 100% | 100% |
| Fog desert used outside a kept 30%, grazing included | 36%–195%; feasible 36%–100% | 55%–292%; feasible 55%–100% | 73%–390%; feasible 73%–100% |
| Power, TW | 1.42–8.41; feasible 1.42–4.06 | 2.13–12.6; feasible 2.13–3.95 | 2.84–16.8; feasible 2.84–3.91 |
| Added dimming, power from the array | 0.017–0.161%; feasible 0.017–0.065% | 0.0255–0.241%; feasible 0.0255–0.0562% | 0.034–0.322%; feasible 0.034–0.0511% |
| Added dimming, fusion | 0.0247–0.354%; feasible 0.0247–0.122% | 0.037–0.531%; feasible 0.037–0.095% | 0.0494–0.708%; feasible 0.0494–0.0802% |
| Water from the seas, km³ a year | 983–4,580; feasible 983–2,490 | 1,470–6,880; feasible 1,470–2,570 | 1,970–9,170; feasible 1,970–2,630 |
| Phosphorus lost, Mt a year | 0.176–0.72; feasible 0.176–0.483 | 0.264–1.08; feasible 0.264–0.523 | 0.352–1.44; feasible 0.352–0.523 |
| Phosphorus discharged from glass, Mt a year | 0–2.76; feasible 0–0.861 | 0–4.14; feasible 0–0.501 | 0–5.51; feasible 0–0.256 |
| Phosphorus to fish sludge, Mt a year | 1.13–19.2; feasible 1.13–7.04 | 1.7–28.8; feasible 1.7–5.27 | 2.26–38.4; feasible 2.26–4.14 |
| N2O, ppb a year | 0.251–1.69; feasible 0.251–0.781 | 0.377–2.54; feasible 0.377–0.732 | 0.503–3.38; feasible 0.503–0.705 |
| Methane, Mt a year (steady ppm) | 1.2–12 (0.252–11.2); feasible 1.2–7.3 (0.252–2.72) | 1.8–18 (0.377–16.8); feasible 1.8–6.94 (0.377–1.93) | 2.4–24 (0.503–22.4); feasible 2.4–5.79 (0.503–1.38) |

**Power-heavy.** Vegetables and fruit in lamp-lit vertical farms, the animal part as microbial protein on power, the
rest of the plant foods on rain-fed fields.

| | 4 billion | 6 billion | 8 billion |
|---|---|---|---|
| Fits its zones, kept 30% | both ends | both ends | both ends |
| Fits its zones, kept half | both ends | both ends | both ends |
| Land used, million km² | 0.896–2.53 | 1.34–3.8 | 1.79–5.87 |
| &nbsp;&nbsp;rain belt | 0.614–1.72 | 0.923–2.6 | 1.23–4.27 |
| &nbsp;&nbsp;highlands, towns included | 0.282–0.813 | 0.421–1.19 | 0.558–1.6 |
| Dry land left wild outside a kept 30% | 95% to 85% | 92% to 77% | 89% to 64% |
| Dry land left wild outside a kept half | 92% to 78% | 89% to 61% | 85% to 37% |
| Rain belt left wild outside a kept 30% | 90% to 71% | 84% to 56% | 79% to 28% |
| Rain-fed land left wild outside a kept 30% | 87% to 65% | 81% to 48% | 75% to 30% |
| Power, TW | 0.539–3.39 | 0.808–5.09 | 1.08–6.78 |
| Added dimming, power from the array | 0.00644–0.0648% | 0.00966–0.0973% | 0.0129–0.13% |
| Added dimming, fusion | 0.00936–0.143% | 0.014–0.214% | 0.0187–0.286% |
| Night energy, times a European's night use | 0.0709–0.457 | 0.0709–0.457 | 0.0709–0.457 |
| Water from the seas, km³ a year | 0 | 0 | 0 |
| Phosphorus lost, Mt a year | 0.176–0.72 | 0.264–1.08 | 0.352–1.44 |
| Phosphorus in microbial biomass, Mt a year | 0.365–1.18 | 0.548–1.77 | 0.731–2.36 |
| N2O, ppb a year | 0.0259–0.0839 | 0.0388–0.126 | 0.0517–0.168 |
| Methane, Mt a year (steady ppm) | 1.2–12 (0.252–11.2) | 1.8–18 (0.377–16.8) | 2.4–24 (0.503–22.4) |

**Fully lamp-lit.** Every plant food in lamp-lit vertical farms on a potato-like staple, the animal part as microbial
protein: no farmland.

| | 4 billion | 6 billion | 8 billion |
|---|---|---|---|
| Fits its zones, kept 30% | both ends | both ends | both ends |
| Fits its zones, kept half | both ends | both ends | both ends |
| Land used, million km² | 0.05–0.2 | 0.075–0.3 | 0.1–0.4 |
| &nbsp;&nbsp;highlands, towns included | 0.05–0.2 | 0.075–0.3 | 0.1–0.4 |
| Dry land left wild outside a kept 30% | 100% to 99% | 100% to 98% | 99% to 98% |
| Dry land left wild outside a kept half | 100% to 98% | 99% to 97% | 99% to 97% |
| Rain belt left wild outside a kept 30% | 100% | 100% | 100% |
| Power, TW | 20.8–212 | 31.2–318 | 41.5–425 |
| Added dimming, power from the array | 0.248–4.06% | 0.372–6.09% | 0.496–8.11% |
| Added dimming, fusion | 0.361–8.93% | 0.541–13.4% | 0.721–17.9% |
| Night energy, times a European's night use | 2.83–28.9 | 2.83–28.9 | 2.83–28.9 |
| Water from the seas, km³ a year | 0 | 0 | 0 |
| Phosphorus lost, Mt a year | 0.176–0.72 | 0.264–1.08 | 0.352–1.44 |
| Phosphorus in microbial biomass, Mt a year | 0.365–1.18 | 0.548–1.77 | 0.731–2.36 |
| N2O, ppb a year | 0.0813–0.202 | 0.122–0.302 | 0.163–0.403 |
| Methane, Mt a year (steady ppm) | 1.2–12 (0.252–11.2) | 1.8–18 (0.377–16.8) | 2.4–24 (0.503–22.4) |

**Power-made.** No farmland: microbial protein for half the reference diet's protein and fungal protein fermented on
sugar made from power for the other half, sugar made by electrochemistry for the food energy, vegetables and fruit from
lamp-lit farms.

| | 4 billion | 6 billion | 8 billion |
|---|---|---|---|
| Fits its zones, kept 30% | both ends | both ends | both ends |
| Fits its zones, kept half | both ends | both ends | both ends |
| Land used, million km² | 0.05–0.2 | 0.075–0.3 | 0.1–0.4 |
| &nbsp;&nbsp;highlands, towns included | 0.05–0.2 | 0.075–0.3 | 0.1–0.4 |
| Dry land left wild outside a kept 30% | 100% to 99% | 100% to 98% | 99% to 98% |
| Dry land left wild outside a kept half | 100% to 98% | 99% to 97% | 99% to 97% |
| Rain belt left wild outside a kept 30% | 100% | 100% | 100% |
| Power, TW | 2.79–11.2 | 4.19–16.8 | 5.59–22.4 |
| Added dimming, power from the array | 0.0334–0.214% | 0.0501–0.321% | 0.0667–0.427% |
| Added dimming, fusion | 0.0485–0.471% | 0.0727–0.706% | 0.097–0.941% |
| Night energy, times a European's night use | 0.379–1.52 | 0.379–1.52 | 0.379–1.52 |
| Water from the seas, km³ a year | 0 | 0 | 0 |
| Phosphorus lost, Mt a year | 0.176–0.72 | 0.264–1.08 | 0.352–1.44 |
| Phosphorus in microbial biomass, Mt a year | 0.707–2.29 | 1.06–3.43 | 1.41–4.58 |
| N2O, ppb a year | 0.00663–0.0164 | 0.00995–0.0247 | 0.0133–0.0329 |
| Methane, Mt a year (steady ppm) | 1.2–12 (0.252–11.2) | 1.8–18 (0.377–16.8) | 2.4–24 (0.503–22.4) |

**Sea-supported.** The animal part from the fresh seas, half as macroalgae and half as filter-feeding shellfish; the
plant foods on rain-fed fields.

| | 4 billion | 6 billion | 8 billion |
|---|---|---|---|
| Fits its zones, kept 30% | neither end (the waters' phosphorus short) | neither end (the waters' phosphorus short) | neither end (the waters' phosphorus short) |
| Fits its zones, kept half | neither end (the waters' phosphorus short) | neither end (the waters' phosphorus short) | neither end (the waters' phosphorus short) |
| Land used, million km² | 0.971–2.74 | 1.46–4.2 | 1.94–6.61 |
| &nbsp;&nbsp;rain belt | 0.668–1.87 | 1–2.92 | 1.34–4.89 |
| &nbsp;&nbsp;highlands, towns included | 0.303–0.868 | 0.452–1.28 | 0.599–1.72 |
| Dry land left wild outside a kept 30% | 94% to 83% | 91% to 74% | 88% to 60% |
| Dry land left wild outside a kept half | 92% to 77% | 88% to 57% | 83% to 37% |
| Rain belt left wild outside a kept 30% | 89% to 69% | 83% to 51% | 77% to 18% |
| Rain-fed land left wild outside a kept 30% | 86% to 62% | 79% to 43% | 73% to 24% |
| Sea farms, thousand km² | 121–981 | 181–1,470 | 241–1,960 |
| &nbsp;&nbsp;phosphorus drawn against what the waters' life can spare | 15.4–97.4 times | 23.1–146 times | 30.8–195 times |
| &nbsp;&nbsp;share of the seas shallower than 50 m / 200 m | 36%–295% / 9%–72% | 54%–442% / 13%–107% | 73%–590% / 18%–143% |
| Power, TW | 0.025–0.0576 | 0.0375–0.0864 | 0.05–0.115 |
| Added dimming, power from the array | 0.000299–0.0011% | 0.000448–0.00165% | 0.000597–0.0022% |
| Added dimming, fusion | 0.000434–0.00242% | 0.000651–0.00364% | 0.000868–0.00485% |
| Water from the seas, km³ a year | 0 | 0 | 0 |
| Phosphorus lost, Mt a year | 0.176–0.72 | 0.264–1.08 | 0.352–1.44 |
| Phosphorus drawn from the seas, Mt a year | 0.865–2.05 | 1.3–3.08 | 1.73–4.1 |
| N2O, ppb a year | 0.0209–0.0735 | 0.0314–0.11 | 0.0419–0.147 |
| Methane, Mt a year (steady ppm) | 1.2–51.1 (0.252–47.8) | 1.8–76.7 (0.377–71.7) | 2.4–102 (0.503–95.6) |

## The author's wild share

On 11 October 2026 the author set the wild share: at least half of the land and of the dark sky wild is preferred, and
it can come down to 30% wild with another 20% hybrid, a compromise between keeping 30% and keeping half
([decisions](../../decisions.md#people-on-the-moon-and-beyond-11-october-2026)); this study applies each rule to every
zone. The hybrid land is heavily manicured megaparks and great gardens, food gardens among them, for leisure and
entertainment, with the small shops, restrooms and places to rest that so much land needs; residents on the ground and
aloft use them, and visitors too, and heavy automation keeps them. The mixes run under four land rules
(`basis.land_rules`; *screen*). The kept 30% and kept half close that share of each zone to fields, glass, towns and
grazing and keep it wild. Under the compromise (`wild_30_hybrid_20`) the parks take the hybrid fifth, so fields, glass,
towns and grazing keep to the same half as under a kept half. At the author's request the earlier reading, half wild
with another 20% hybrid (`wild_50_hybrid_20`), runs beside them for comparison with the same parks; it leaves them
30% of each zone. The food gardens' harvest stays outside these accounts.

| Mix, 4 / 6 / 8 billion | Kept half and the compromise | Half wild with 20% hybrid |
|---|---|---|
| Rain-fed first | to 60% / 40% / 30% of the way | to 36% / 24% / 18% of the way |
| Greenhouse-heavy | both ends / to 90% / to 63% | to 80% / 45% / 24% |
| Fog-desert glass | to 36% / 13% / neither end | to 8% / neither end / neither end |
| Power-heavy, fully lamp-lit, power-made | both ends | both ends |
| Sea-supported | neither end | neither end |

The rain alone feeds 6.5–19.3 billion under a kept half and the compromise, and 3.9–11.6 billion with half wild and
20% hybrid. At 8 billion the power-heavy mix's fields take 1.69–6.97 million km² under a kept half and the compromise,
with 0–235 km³ a year from the seas. With half wild and 20% hybrid they take 1.77–4.04 million km², and the high-cost
end grows 1.12 million km² of its crops' equivalent by watering those fields: 3,430 km³ a year from the seas, 0.56 TW of
lift and 7.3 TW of food power against 6.8.

## The author's blend

On 12 October 2026 the author chose a blend of the three mixes that fit under the wild share at 4–8 billion,
power-heavy, power-made and fully lamp-lit, left its weights to an optimum, and added small amounts of heritage meat
from the homesteads and hobby farms as a luxury
([decisions](../../decisions.md#automation-the-network-and-food-12-october-2026)). [blend.py](blend.py) builds each
blend as one case from the three mixes' rules (*screen*): power-heavy's share grows its fields around every resident's
towns, wettest land first, and the other two add their power, gases, phosphorus and buildings. A blend of one mix alone
repeats that mix's case.

**The rule.** Weights in steps of 5%, each mix at least a tenth of the diet, so all three stay in use and can grow when
another falters, and none above half, so that losing any one leaves at least half the food. Among those, the blend that
takes the least sunlight from wild life at 8 billion under the compromise, at the high-cost end and then over both ends:
the fields' area plus the added dimming under fusion, the planning assumption for local power, as a share of the whole
surface. Power-made uses no fields and little power and takes the cap; fully lamp-lit keeps the floor, since each tenth
of the diet under lamps takes 4.2–42 TW of power, 8.3–129 TW of heat with fusion; power-heavy carries the rest. The
limits set the weights: without them power-made alone takes least, and each ten points more of it would spare about
440,000 km² at the high-cost end. At the low-cost end lamps cost less than fields, and 10/50/40 would take 0.78% of the
Moon's sunlight against the blend's 1.9%.

| At 8 billion, the compromise | The blend, 40/50/10 | Least heat, 50/40/10 | Even thirds | Power-heavy alone |
|---|---|---|---|---|
| Fields, million km² | 0.68–1.86 | 0.85–2.33 | 0.56–1.55 | 1.69–6.97 |
| Heat with fusion, TW | 14.8–171 | 13.9–166 | 32.1–458 | 2.16–20.7 |
| Added dimming with fusion | 0.13–2.4% | 0.12–2.3% | 0.28–6.4% | 0.019–0.29% |
| Sunlight taken from wild life, share of the Moon's | 1.9–7.3% | 2.3–8.4% | 1.8–10% | 4.5–19% |

The blend's power is 7.4–56 TW, 10.2–77.5 TW of heat with power from the array. Its fields and towns leave 81–93% of the
dry land outside the kept half wild and the middle latitudes free of fields, so the homestead region, the middle-wet
zone's 1.39 million km² outside the kept half, is the homesteads' own. It needs no water from the seas, adds 0.044–0.124
ppb of N2O a year, against 0.052–0.17 for power-heavy alone, and its vertical farms take 1.9–22 Gt of building. At 4 and
6 billion its heat is 7.4–85 and 11–128 TW with fusion; under the comparison's half wild with 20% hybrid its fields take
0.68–2.17 million km² and it still fits.

**Heritage meat.** A basket of the author's examples, by weight 40% chicken, 30% heritage pork, 25% beef from dairy
herds (homestead dairy and vaca vieja) and 5% beef from beef herds (some wagyu), a design guess, counted as fat- and
bone-free meat. At Earth's global per-kilogram land use (Poore & Nemecek 2018) it takes 37.2 m² a kilogram. FAO GLEAM
gives beef 67.8 kg CO2-eq per kg of carcass from beef herds and 18.4 from dairy herds, 43% and 46% of it enteric
fermentation (Gerber et al. 2013; Opio et al. 2013); at a carcass's yield of 0.667–0.705 kg of boneless or retail beef
(USDA ERS 1992), the cattle emit 1.65–1.75 kg of methane a kilogram of beef from beef herds and 0.48–0.51 from dairy
herds, and the basket 0.20–0.21 kg. Pigs and poultry emit little.

| At 8 billion, a person a year | 0.5 kg | 1 kg | 2 kg |
|---|---|---|---|
| Land, thousand km² | 149 | 298 | 595 |
| Share of the homestead region | 11% | 21% | 43% |
| Share of the diet's animal protein, at 20–25% protein | 1.1–1.4% | 2.2–2.8% | 4.5–5.6% |
| Cattle methane, Mt a year | 0.81–0.86 | 1.6–1.7 | 3.2–3.4 |
| Steady methane, ppm | 0.17–0.80 | 0.34–1.6 | 0.68–3.2 |

At 8 billion each kilogram of beef a person a year emits 3.8–4.1 Mt of methane from dairy herds and 13–14 Mt from beef
herds. The heritage meat stands within the diet's animal part, so the blend's power is left as computed; its feed and
grazing are inside the land figure.

## What the physics allows a lamp-lit diet

The fully lamp-lit diet for 8 billion people uses 41.5–425 TW (*screen*), against 12–36 TW for the residents' other
energy. The array intercepts 1,900–2,400 PW and serves through the night, so the power is there; the heat is the bound.
With power from the array the Moon takes 1.38 W of heat for each watt used, and the population study answers each 72–115
TW of heat with 1% more dimming, so the diet needs 0.50–8.1% more dimming at 8 billion and 0.25–4.1% at 4 billion.
Fusion on the Moon releases 2.0–3.0 W of heat a watt and needs 0.72–17.9% at 8 billion. That dimming takes the same
share of sunlight from every wild plant, field and sky reef, and the night's lamps draw 2.8–29 times a European's night
use. At the low-cost end (the most efficient LEDs, the best potato light use, the smallest loss and variety allowances)
the whole diet costs 5.2 kW a resident; at the high-cost end 53 kW.

## Departures from the known designs

- **Crop water.** Ecology's crop water need, 4.79–10.8 mm a day at the equator, comes from Earth crops' water-use
  efficiency, which ecology takes from Tang et al. (2014), raised by the design air's CO2 to 1.26–2.83 g C per kg. The
  sunlight the ground absorbs under the design clouds evaporates 4.94–5.02 mm a day there (*screen*), so the crop
  model's yields need an efficiency of at least 2.7 g C per kg, near the top of that range. The study holds the need at
  4.79–5.02 mm and raises ecology's rain shares in the same proportion, 2.15 times at the high-cost end, which brings
  the rain-fed land from 2.57–5.78 to 5.54–5.78 million km² of equivalent on this screen's bands. At ecology's lowest
  efficiency the crop model's yields cannot be reached, watered fields grow at most 46% of them at the equator, and the
  rain feeds 4.2 billion people outside a kept 30% at that end.
- **The same energy limit by latitude and under glass.** The crop model's growth falls more slowly with latitude than
  the sunlight does, so watered fields grow 92–97% of the model's production at 15–30°, 87–92% at 30–45° and 73–78% at
  45–60° (*screen*). Under glass the sunlight absorbed inside caps transpiration at 2.0–3.4 mm a day by zone.
- **Band areas.** This screen counts each latitude band's dry land on the atlas grid (5.44 million km² within 15°, 5.58
  at 15–30°); ecology scaled the climate model's coarse land to the atlas's total (6.10 and 3.90). Both divide each
  band's production by the equator's. At ecology's water need the rain-fed equivalent comes to 2.57–5.78 million km²,
  against ecology's 2.55–5.73: the band areas add 1.5% and ecology's shares, given to two figures, take off 0.4–0.7%.
- **Zones.** The fog desert here is dry land at 30–60° below 4 km with under 0.5 mm of rain a day, 4.82 million km²
  (5.11 at all heights); ecology took the aerosol study's 6.2 million km². The upper grazing allowance therefore fits
  fewer people here than on ecology's area. The 30–45° band's fields hold to its wet lowland, taken as the wet middle
  latitudes (30–60°), and its own highlands.
- **Soilless N2O.** Glass and lamp farms emit at the measured soilless 0.72–1.19% of applied nitrogen, above the lunar
  design's 0.2–0.5% target for fields.
- **Microbial protein** is capped at half the reference diet's protein by nucleic acids; it carries its own ammonia in
  its electricity accounts, so the mixes' nitrogen fixation covers the protein grown on crops.
- **Farm yields.** The mixes use the crop model's potential; the farm-yield sensitivity sits beside the systems.

## Basis

| Zone | Dry land, million km² | Mean height above the sea, km | Mean distance to a sea or great lake, km | Rain, mm a day | Mean latitude |
|---|---|---|---|---|---|
| Highlands, ground 4 km and more above the sea | 3.90 | 5.46 | 352 | 4.49 | 24.5° |
| Rain belt, within 30° | 8.48 | 1.50 | 94.9 | 3.13 | 15.7° |
| Wet middle latitudes, 30–60°, rain 0.5 mm a day or more | 2.78 | 1.88 | 152 | 1.59 | 38.2° |
| Fog desert, 30–60°, less rain | 4.82 | 1.43 | 69.2 | 0.176 | 47.2° |
| Polar land, beyond 60° | 3.42 | 1.54 | 59.9 | 0.0555 | 70.2° |

The zones partition the 23.4 million km² of dry land on the 28% atlas, with the design climatology's rain (*screen*).
The seas and great lakes cover 10.2 million km², of which 333,000 km² is shallower than 50 m and 1.37 million km²
shallower than 200 m. The crop model's yield relative to the equator is 0.951 in the rain belt, 0.924 on the highlands
and 0.721 in the fog desert (canopy fruit model, near side).

The diet is ecology's plant-rich lunar design: 2,500 kcal a day eaten, 131–348 m² of crops for its plant foods and
141–401 m² with its animal part, and 0–1,000 m² of grazing. Its protein and water come from the EAT–Lancet reference
diet (Willett et al. 2019) with USDA FoodData Central composition: 94.2 g of protein a day, 24.3 g of it from animal
foods, and 300 g of vegetables and 200 g of fruit (*derived*; the food standing for each group is a design guess).

The pairs run from the low-cost end to the high-cost end: every uncertain input at its favourable value, then at its
unfavourable one. Design guesses, each a range: the composition-table food for each diet group, the CO2 held in glass,
the share of rain the gutters keep (80–95%), the emitter pressure and the water main.

## Bound inputs

| Input | Branch of origin | sha256 |
|---|---|---|
| `biosphere/ecology/results/people.json` | domain/ecology | `758a8dcf68c6eb58…` |
| `biosphere/ecology/results/first_screen.json` | domain/ecology | `d2a42b6e51f60250…` |
| `biosphere/ecology/people_and_land.md` | domain/ecology | `fbd83a2137a12eb5…` |
| `resources/results/phosphorus.json` | domain/resources | `4089899e2b70e8b7…` |
| `resources/phosphorus.py` | domain/resources | `48f3ef951dc6d781…` |
| `resources/phosphorus.md` | domain/resources | `afa8db4d5e8b698c…` |
| `resources/results/sea_phosphorus.json` | domain/resources | `47b600895919319d…` |

The runner reads them in the tree, which has held every former branch's inputs since the joint integration of 12
October 2026, and refuses bytes whose sha256 differs; the full hashes and Git blobs are in `run.py` and in the product's
`producer.bound_inputs`. Lifting the 1.32 kg of food a day that 4 billion residents aloft eat to 10 km takes 0.99 GW at
resources' ideal 16.1 kJ per kg lifted, 1.4 GW at its 70% efficiency, and the phosphorus rides in that food (*derived*).

## Reproduce and check

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 nice -n 10 python3 -m research.studies.food_systems.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 nice -n 10 python3 -m research.studies.food_systems.blend
python3 -m pytest -q research/studies/food_systems
```

The run reads three ignored grids (the atlas, drainage and design climatology), the committed products listed in
`producer.inputs` and the bound inputs above. The tests check that the product regenerates; the producer, input and
bound hashes; that every parameter names a source in `sources.json`; that the ecology basis reproduces the bound
product; and conservation and limits: heat and dimming follow power, water closes in every glasshouse, people and land
close in every mix, phosphorus balances, fields stay outside the kept land, fog nets need wind and fog, the Weibull
partial mean matches its integral, lamp power scales with light use, the fill spills to the fog desert when the wet
bands are full, every case's feasibility follows its zones, the land rules carry the author's wild share and keep a
kept half apart from the compromise, and the compromise's parkland gives the kept half's results in every case. The
regeneration test skips when the ignored grids or the bound inputs are absent. `checks.json` records the review's
findings, their fixes and the test results.

The blend reads the same basis and writes [results/food_blend.json](results/food_blend.json); its tests check the rule,
that each mix alone repeats the food study's case, that the blend lies between its parts, the heritage meat's arithmetic
and that the product regenerates.

## What remains

- C′ under the author's blend: the distributions study counts power-heavy's fields and heat; the blend's heat, 15–171 TW
  with fusion at 8 billion, and its fields replace them in the joint account. The heritage meat's homesteads, their
  grazing, feed and herds in detail.
- The crop model has no water yet (ecology's first open question). Its water-use efficiency decides between the rain-fed
  land here and ecology's lower end, and farm yields below potential add a quarter to a doubling of field land.
- Fog and wind paired snapshot by snapshot, from a CM1 ring rerun that writes them, would replace the fog-net bracket.
- The parks under the author's compromise: their water on dry land, in the dry highlands and in the fog desert and
  polar lowland, which hold 1.65 million km² of them; the energy of the automation that keeps them; and their night
  light within the hybrid fifth of the dark sky.
- The warming of N2O and methane in the deep column is uncomputed; N2O accumulates in every mix.
- The waters' harvest ceiling, 30% of each water's dissolved input, is a design guess of resources' study; particles
  that release phosphorus as in Earth's estuaries would raise it to 44–377 kt a year.
- Dew in the fog desert needs the night's surface energy balance, which no product gives yet.
- Methane under unfed seaweed and shellfish farms, energy per kg of fat made from CO2 and hydrogen, and an electricity
  figure per person for a fully lamp-lit staple diet were found in no source read.
- Not read in full, and so not used: Mortensen (1987) on CO2 response; Stanghellini's water per kilogram by greenhouse
  type; Davies & Paton (2005) on the seawater greenhouse; Sillman et al. (2019) on bacterial protein; Plappally &
  Lienhard (2012) on water's energy; Gultepe et al. (2007) and Tomaszkiewicz et al. (2015) on fog and dew; Scheschonk et
  al. (2019) on kelp in the polar night, read through Li et al. (2020).
