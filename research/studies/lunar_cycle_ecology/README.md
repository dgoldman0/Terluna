# Ecology of the lunar cycle

27 September 2026. A register of design ideas for the Open Moon's plants, animals and
ecosystems, and for the nutrient and mineral cycles that keep them going. The ideas were
gathered in the working sessions of 26–27 September 2026, some from the July 2026
bundle ([06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md](../../../archive/july_2026_knowledge_bundle/terluna_moon_terraforming_bundle/06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md)).
Each is set against the models and the literature, with a verdict, and the
[selection](#selection) recommends which to carry forward. The author decides; none
of this is in the [decisions register](../../decisions.md).

This is an ideas register with light calculations. It models no ecosystem.
`python -m research.studies.lunar_cycle_ecology.run` (about 20 seconds) writes
[results/lunar_cycle_ecology.json](results/lunar_cycle_ecology.json) (schema
`terluna.research.lunar-cycle-ecology/1`). It reads the biosphere's
[plant and fruit products](../../../biosphere/canopy/README.md), the geography
[atlas and drainage](../../../geography/README.md), the design air and the shared
constants. [sources.json](sources.json) records each source and how far it was
checked.

## How to read the register

Each idea carries one evidence label and one verdict.

- **Computed:** from this repository's models, with the product named.
- **Earth evidence:** measured on Earth and cited.
- **Design requirement:** what a model says an organism would have to do. Nobody has
  tested whether one can.
- **Idea:** reasoned but not yet computed or tested.

The verdicts are **adopt** (carry into the design and the models), **in places** (use
in particular landscapes or in farming), **test first** (promising, with a specific
test before adoption) and **set aside**.

## The setting

What the models say about the world these organisms would live in. Figures are for
the equator unless stated.

- **The cycle.** 708.7 hours from noon to noon, 354 of them with the Sun up. The tall
  sky keeps photosynthetic light above 1 µmol m⁻² s⁻¹ for 58 hours after sunset, so
  238 hours are dark at the equator, 217 at 30° and 34 at 60°
  ([plant.json](../../../biosphere/canopy/results/plant.json)).
- **Light.** Overhead sunlight carries 1,560 µmol m⁻² s⁻¹ of photosynthetic photons,
  against Earth's 2,250, with 37% of it diffuse and a UV index of 0.11
  ([surface light](../../../illumination/surface_light/README.md)). A deep canopy fixes
  18% more carbon over the cycle than the same canopy under Earth's sky
  ([canopy](../../../biosphere/canopy/README.md)).
- **Temperature.** Land air in the chosen climate runs 21–23.5 °C by day and 19–20 °C
  by night at every latitude on both sides; the coldest hour anywhere is 18.7 °C. The
  night is long, dark and warm.
- **Carbon.** The evergreen stand grows 229 g C per m² per cycle, 7.8 per 24 hours
  against Earth's 5.4 for the same stand, and needs a store of 46 g C per m² for the
  night, or 10 if it idles at a quarter of its daytime upkeep.
- **Darkness.** Earthlight gives the near side about 2 lux at midnight; far-side nights
  are dark ([illumination](../../../illumination/README.md)).
- **Air.** 1.2 atm at 1.43 kg/m³. Each m² of ground lies under 74,900 kg of air holding
  14,500 kg of O₂ and 46 kg of CO₂.
- **Water.** Seas cover 28% of the Moon, and rain-fed lakes a further 11.6%. Rivers
  deliver 279,000 m³/s to the seas, 0.32 m a year over the land
  ([geography](../../../geography/README.md)). The median point on land lies 476 km
  from a sea and 50 km from standing water.
- **Dusk moves.** The line between day and night crosses the ground westward at
  15.4 km/h at the equator, 13.3 at 30° and 7.7 at 60°.

## The register

### A. Light and the canopy

- **A1. Deep canopies under a diffuse sky.** *Computed · adopt.* The Moon's canopy
  fixes 18% more carbon than Earth's over the cycle (11% more at equal CO₂) because the
  diffuse sky feeds the shaded leaves. The gain grows with leaf area, from 9% at a leaf
  area index of 1 to 18% at 5–8, and leaves alone put the optimum near 6
  ([canopy](../../../biosphere/canopy/README.md)). Design canopies deep.
- **A2. Spectral tuning:** paler upper leaves, photosynthesis into the far red, clumped
  foliage. *Computed · set aside.* They add 3%, 6% and 3% on the Moon, against 8%, 11%
  and −2% on Earth. The diffuse sky already does most of what they do.
- **A3. C4 and "hyper-C4" plants.** *Computed · in places.* Land air stays at 19–26 °C,
  CO₂ is 48.6 Pa, and the canopy works mostly in dim diffuse light, where C4's extra
  two ATP per CO₂ cost most. With the leaf model's photorespiration (Bernacchi et al.
  2001 kinetics), the Moon's 48.6 Pa of CO₂ against Earth's 40.5 moves the temperature
  where a C4 leaf overtakes a C3 leaf 3.6 °C higher.
  C4's water economy can still earn it dry or hot sites.
- **A4. Extra photoprotection, reflective hairs, folding leaves** (July bundle §4).
  *Computed · set aside as general traits.* Peak light is 30% below Earth's and the UV
  index is 0.11. In tomato, the best-studied case of injury under continuous light,
  the trigger was timing, with light continuity and photo-oxidative pressure tested
  and excluded (Velez-Ramirez et al. 2017).
- **A5. Rapid chloroplast cycling, or a daily rest in the long day.** *Earth evidence
  and computed · set aside as a design goal.* Chloroplasts already switch protective
  quenching on and off within seconds to minutes. A rest that halves photosynthesis for
  eight hours in every 24 would cost about a quarter of the growth. Keeping water status
  in a two-week day is the one reason for a daily stomatal rhythm, and it needs the
  water model (see the next steps).

### B. Through the night

- **B1. Evergreen plants that store and idle, on a carbon trigger.** *Computed ·
  adopt.* Idling leaves, stems and roots at a quarter of their daytime upkeep cuts the
  store from 46 to 10 g C per m². Idling as soon as photosynthesis stops covering
  upkeep needs half the store of idling only in full darkness, 10 against 19
  ([plant.json](../../../biosphere/canopy/results/plant.json)).
- **B2. Using the long twilight.** *Computed · adopt.* Twilight is worth 14% of the
  Earth-like store and a quarter of the idling one, and it leaves only 34 dark hours at
  60°.
- **B3. Regrowing the canopy each cycle** (the July bundle's shedding before dusk).
  *Computed · set aside.* A new canopy every cycle costs 264 g C per m², built at dawn
  from reserves. Growth falls to 28 at the equator, 3 at 30° and none at 60°.
- **B4. Frost protection, antifreeze, insulation, heat sharing in plant mats** (July
  bundle §5). *Computed · set aside at low and middle latitudes.* The bundle assumed
  cold nights. The chosen climate's land nights stay at 18.7–20.1 °C everywhere the
  models reach, which makes the night warm enough for an active food web (E).
- **B5. Leaves that survive about 240 dark hours, and metabolism that slows and
  restarts on cue.** *Test first.* The physiology every evergreen strategy rests on. In
  Arabidopsis a leaf darkened on its own senesces, while darkening the whole plant holds
  senescence back (Weaver and Amasino 2001), and the lunar night darkens whole plants.
  The *kin10 kin11* mutants, which lack the SnRK1 energy sensor, mobilize starch poorly
  at night (Baena-González et al. 2007).

### C. Clocks and cues

- **C1. Earth's 24-hour clock under two weeks of light.** *Earth evidence · in
  places.* Cultivated tomato is injured under continuous light because its clock keeps
  scheduling a night that never comes and turns down the light-harvesting gene CAB-13
  on that schedule. A wild tomato's allele of that gene gives tolerance, and a daily
  temperature swing or alternating red and blue light prevents the injury
  (Velez-Ramirez et al. 2014, 2017). Farms can supply such cues and grow tolerant
  cultivars. Wild plants need clocks that loosen their hold on light harvesting, since
  the open Moon offers no 24-hour signal.
- **C2. A long, dusk-set night timer.** *Earth evidence and design requirement ·
  adopt.* Arabidopsis divides its leaf starch by the hours to the dawn its clock
  expects, and in 28-hour days its starch still runs out 24 hours after the last dawn
  (Graf et al. 2010; Scialdone et al. 2013). A plant on an Earth clock would budget the
  lunar night's store for one day. The design needs a timer sized for about 240 dark
  hours, or rationing on the plant's own carbon balance.
- **C3. The plant's own carbon balance as the idling cue.** *Computed · adopt.* Plants
  already switch between growth and starvation programs through the energy sensors
  SnRK1 and TOR (Baena-González et al. 2007; Xiong et al. 2013). Idling on the
  plant's own deficit halves the store (B1).
- **C4. Two clocks** (July bundle §2): daily rhythms in leaf chemistry inside a long
  environmental clock. *Idea · adopt with a condition.* It works if the daily clock
  does not command light harvesting (C1). People and farm animals keep 24-hour days.
- **C5. Animal life cycles timed by the Moon.** *Earth evidence · adopt.* Earth
  animals already keep lunar calendars:
  - The marine midge *Clunio marinus* emerges in millions at the low tides around new
    and full moon, times it with circalunar and circadian clocks, and lives as an adult
    for a few hours (Kaiser et al. 2016).
  - Samoan palolo worms cast off their spawning segments at the third quarter moon of
    October or November (Caspers 1984).
  - At least 32 coral species on the Great Barrier Reef spawn together a few nights after
    the late-spring full moons (Harrison et al. 1984).

  On the Open Moon the month is the day, and dawn and dusk are its sharpest cues.

### D. Where the long day's carbon goes

- **D1. Storage organs and wood as sinks.** *Computed · adopt.* Each cycle the
  equatorial stand grows 229 g C per m² and fills a store of up to 46 g C for the
  night, all while the Sun is up. A leaf that cannot export its sugar turns its own
  photosynthesis down, so the long day needs sinks sized for it: wood, roots, tubers
  and fruit.
- **D2. The day fruit.** *Design requirement · adopt as the fruit target.* A fruit whose
  program fits one sunlit half (165–199 degree-days above 10 °C, a third of today's
  watermelon's 560) is set at sunrise and ripe at sunset. It grows entirely on the
  daylight surplus and leaves the plant's night store as it was. A 7.5 kg fruit fills at
  up to 0.9 kg a day, 3.5 times today's watermelon. Giant pumpkins, bred for more
  phloem into the fruit, take in about 15 kg a day at their peak (Savage et al. 2015),
  and cucumbers grow from 5 to 30 cm in 10–15 days (Wiechers et al. 2011)
  ([fruit model](../../../biosphere/canopy/README.md#fruit-designed-for-the-lunar-day)).
  What it asks: the program compressed, fruit set on a dawn cue without pollination
  (as parthenocarpic cucumbers set), cells formed before sunrise, and sugar loading in
  the last days before dusk.
- **D3. How much goes to fruit.** *Computed · adopt 50–60% for farms.* With 50–60% of
  the plant's new tissue carbon in fruit, the equatorial stand yields 3.1–3.7 kg per m²
  every lunar cycle (3.5–4.2 when the plant idles at night), 51% more than the same
  stand on Earth with the same temperatures. Greenhouse cucumbers swing between 40% and
  90% of their new dry matter (Marcelis's thesis). Wild forests would put less into fruit.
- **D4. Leaves replaced at their own pace.** *Computed · adopt.* Leaves on Earth live
  from under a month to 24 years, and biome means run from 2.5 to 66 months (Wright et
  al. 2004; Reich et al. 1997). Replacing a canopy whose leaves live 6, 12, 24 or 36
  months takes 34, 17, 9 or 6 g C per m² per cycle.
  At a 50–60% fruit share the plant keeps 89–112 g C per cycle, and 66 at 70%. Leaves
  need not come and go each cycle; what the plant keeps pays for leaves, roots and wood.
- **D5. Fruit picked young.** *Earth evidence · in places.* Zucchini is picked a few
  days after flowering and cucumbers within about two weeks, so both fit inside one
  lunar day on today's programs.
- **D6. Today's watermelon.** *Computed · baseline.* It takes 45 days at the equator's
  temperatures and lives through a night, growing about a third of its mass in the
  dark. Full size then needs a store of 116 g C per m², two and a half times the
  plant's own.

### E. The night food web

- **E1. The dusk fruit fall feeds the night.** *Computed · adopt.* The day fruit ripens
  at sunset, so the fruit falls at dusk. A wild forest putting 30% of its new tissue
  into fruit drops 116 g of sugar per m² each cycle. Spread over the dark hours that is
  2.3 W per m², about what the plants themselves burn in the dark (2.4). Nights stay at
  19–20 °C, warm enough for insects, yeasts and fungi to keep working. On Earth,
  fermenting fruit draws fruit flies through its yeasts' volatiles, and the flies need
  the yeast to develop (Becher et al. 2012).
- **E2. No oxygen crisis at night.** *Computed · adopt, as a correction.* The July
  bundle (§9) warned of an oxygen-consuming decomposition pulse late in the night. The
  night's respiration plus the decay of a whole 60% fruit fall uses 3 × 10⁻⁵ of the O₂
  above each m², because the 1.2-atm air column holds 14,500 kg of it. The same carbon
  adds 1.5% to the CO₂ above that m² if the air stays put. What deserves modelling is
  CO₂ gathering under still night air near the ground.
- **E3. Glowing lures.** *Earth evidence and computed · test first.* Light already
  lures prey in the dark:
  - Glowworm larvae hang sticky threads beneath their light and catch mostly flies
    (Broadley and Stringer 2001; von Byern et al. 2019).
  - Artificial mushrooms glowing like the luminous fungus *Neonothopanus gardneri*
    caught 42 insects against 12 in dark ones over five nights (Oliveira et al. 2015).
  - Flying insects gather at lights because they turn their backs toward them (Fabian
    et al. 2024), and they come from about 10–25 m away (Truxa and Fiedler 2012; Degen
    et al. 2016).

  Plants have been engineered to glow with the fungal pathway, up to 6.5 × 10¹⁰ photons
  per minute per cm² of flower (Mitiouchkina et al. 2020). A square metre glowing that
  brightly at the fungal system's lowest measured efficiency (0.1%; Kaskova et al.
  2017), with ten times the photon energy allowed for making the luciferin, draws
  0.04 W: under 2% of the night upkeep of a square metre of plants. Spread over the
  ground one lure serves, the cost is negligible. Three cautions:
  - The engineered plants' glow follows their circadian clock and faded by the third
    and fourth day of darkness. A lunar-night lure needs a glow that lasts, or one that
    shines through the first nights after dusk, when the fruit falls and the migrants
    come.
  - Artificial light at night harms moths at most stages of their lives (Boyes et al.
    2021), so the lures' catch has to leave the insect populations standing.
  - Earthlight competes on the near side; far-side nights are dark.
- **E4. Carnivorous plants.** *Earth evidence · in places.* Carnivory pays in sunny,
  moist, nutrient-poor places (Givnish et al. 1984). Traps cost carbon, while their
  benefit levels off (Ellison and Adamec 2011).
  - **Uptake:** carnivorous plants take up 30–76% of a prey's nitrogen and 57–96% of
    its phosphorus, potassium and magnesium (Ellison and Adamec 2011).
  - **Share of their nitrogen:** sundews take about half their nitrogen from prey
    (Millett et al. 2003), and the pitcher plant *Nepenthes mirabilis* about 62%
    (Schulze et al. 1997).
  - **Weakness:** their photosynthesis per unit of nutrient is a fifth to a half of
    other plants' (Ellison 2006), so they lose on rich soil.

  The lunar day's surplus carbon makes traps and lures cheap. Wet ground, poor sands
  and the canopy's epiphytes are their places, and they are the forest end of the
  migrant cycle (H).
- **E5. Taking nutrients from animals without killing them.** *Earth evidence ·
  adopt.* *Nepenthes lowii* feeds tree shrews on the secretions of its lid and takes
  57–100% of its leaf nitrogen from their droppings (Clarke et al. 2009).
  *N. hemsleyana* houses Hardwicke's woolly bats and takes 34% from theirs, while
  catching few insects (Grafe et al. 2011). Roosts, lekking sites and feeding places that
  collect droppings take nutrients from larger migrants without a trap.
- **E6. Dusk followers.** *Computed · idea.* Dusk moves west at 15.4 km/h at the
  equator and 7.7 km/h at 60°. An animal flying west at that pace stays at dusk and meets
  the fruit fall as it happens; this is one form of the July bundle's migrators that
  follow the terminator.

### F. Two landscapes

- **F1. Plains of storage organs.** *Earth evidence · adopt.* On Earth, plants with
  bulbs, corms and tubers are most diverse in the Cape's winter-rainfall region
  (about 2,100 species; Procheş et al. 2006), where most sit out the dry summer dormant
  (Parsons 2000). They favour cooler, drier and more variable climates (Howard et al.
  2019), and large storage organs also mark arid, unpredictable places (Dafni et al.
  1981). The lunar night is the most regular unfavourable period there could be, and a
  storage organ carries a plant through it without wood. Perennial grasses and forbs
  seed at dusk, and grazers and seed-eaters follow.
- **F2. Forests of wood, deep canopies and day fruit.** *Computed and idea · adopt.*
  Wood holds the night's store: tens of g C per m², small next to what Earth's forests
  keep as non-structural carbohydrate (Martínez-Vilalta et al. 2016, as cited in the
  [plant model](../../../biosphere/canopy/README.md)). That store lets a forest spend on
  fruit. The forest keeps the carnivores and nutrient-catching plants of E4 and E5,
  fruit-eaters, decomposers and night insects.
- **F3. Wet margins.** *Computed and idea · in places.* Half the land lies within 50 km
  of standing water and 72% within 100 km. Lake shores and wetlands suit carnivorous
  plants, and lake insects carry nutrients onto the shore (G, H).

### G. Nutrients and minerals

- **G1. Nutrients flow one way.** *Earth evidence and computed · adopt as a design
  problem.* Without volcanism or uplift, weathering and runoff carry rock-derived
  nutrients downhill into lakes and seas, and nothing lifts them back. Soils lose their
  mineral phosphorus as they age:
  - In Hawaii, mineral phosphate fell from 82% of total phosphorus at a 300-year-old
    site to 1% at 20,000 years (Crews et al. 1995).
  - Growth at the oldest site, 4.1 million years, is limited by phosphorus (Vitousek
    and Farrington 1997).
  - Across six long soil sequences, forests decline as phosphorus runs short (Wardle et
    al. 2004).

  At the Moon's runoff and land temperature, the basalt weathering law of the
  [CO₂ study](../atmospheric_co2/README.md) dissolves about 50 g of soil per m² a year,
  and up to 500 from fresh glassy regolith. The top metre weathers through in
  3,000–30,000 years. Fresh regolith supplies phosphorus and potassium for millennia,
  and after that the Open Moon's soils age as Earth's old soils do, well within its
  billion-year horizon.
- **G2. Potassium is the scarce one.** *Earth evidence and computed · adopt.*
  - **The Moon's potassium:** Lunar Prospector puts the surface at 755 ppm on average,
    3% of Earth's upper continental crust at 23,240 ppm (Rudnick and Gao 2003). The
    Procellarum KREEP Terrane averages 2,004 ppm and the feldspathic highlands under 500
    (Prettyman et al. 2006, via Zhu et al. 2013).
  - **Soil samples:** Apollo and Luna soils hold 0.04–0.55% K₂O (McKay et al. 1991).
  - **Phosphorus:** it matches Earth's crust, at 0.05–0.5% P₂O₅ against 0.15%.
  - **Supply:** weathering releases 0.06–0.7 g of potassium per m² a year from mare or
    highland soil, and 0.2–2.3 g from soil rich in KREEP (the potassium-, rare-earth-
    and phosphorus-rich component of lunar rock). Earth's rivers, draining crust with
    thirty times the potassium, would carry 0.34–0.42 g off the same land at the Moon's
    runoff. Plants need as much potassium as on Earth, from rock that holds a thirtieth
    of it.
  - **Where the rich rock lies:** under water. The Near-side Sea floods 99% of Mare
    Imbrium and 60% of Oceanus Procellarum, the core of the KREEP terrane. That sea floor
    is the Moon's best source of potassium and phosphorus, and sea spray and marine life
    can bring some back (G5, H).
- **G3. Farms carry off far more than rivers.** *Computed · adopt.* At a 60% fruit
  share the equatorial stand's harvest is 46 kg of fruit per m² a year, at the model's
  clear-sky potential. With watermelon's composition (USDA FoodData Central 167765) it
  carries off, per m² a year:

  | Nutrient | Carried off by the harvest | Against Earth's rivers at the Moon's runoff |
  |---|---|---|
  | Nitrogen | 45 g | |
  | Potassium | 51 g | 120–150 times |
  | Phosphorus | 5.0 g | 600–1,600 times the dissolved phosphorus |
  | Magnesium | 4.6 g | |
  | Calcium | 3.2 g | |

  Those nutrients have to return to the fields from the settlements that eat the
  harvest, as food waste, sewage and manure. It is the largest nutrient loop on farmed
  land.
- **G4. Potassium from rock and the sea.** *Idea · adopt.* Crushed KREEP basalt or
  breccia (0.53–0.83% K₂O; Taylor et al. 1991) can be spread as rock dust, and potassium
  can be recovered from the seas, which collect it. Farms return what they harvest
  (G3).
- **G5. Dust and sea spray.** *Earth evidence · in places.* Dust and spray already
  carry nutrients inland on Earth:
  - **Dust:** African dust brings the Amazon 7–39 g of phosphorus per hectare a year,
    about what the forest loses to its rivers, 8–40 (Yu et al. 2015).
  - **Spray reach:** ocean aerosols can be traced 2,000 km inland (Meybeck 1994).
  - **Iodine:** most soil iodine comes from the sea through the air (Fuge and Johnson
    1986). Soils more than 50 km inland hold a geometric mean of 2.6 µg/g against 11.6
    near the coast, though soil iodine does not fall simply with distance (Johnson
    2003).

  The median point of lunar land lies 476 km from a sea: within spray's reach, but far
  enough that inland soils will be poor in iodine. People and animals will need iodine
  supplied until the spray and the migrants carry enough.
- **G6. Erosion and sediment.** *Computed · adopt as engineering.* Particulate matter
  makes up 95% of the phosphorus Earth's rivers carry naturally (Meybeck 1982). Erosion at
  Earth's rate would remove 0.06–0.17 g of phosphorus per m² of lunar land a year: 7–50
  times the dissolved loss, and more than any conveyor of animals could return (H4).
  Keep the land covered, trap sediment in lakes and deltas, and recover it; the
  biosphere seed already notes that nutrients buried in sediment might need recovery.

### H. The sea–forest migrant cycle

- **H1. The cycle.** *Idea · adopt as the candidate forest cycle.* Marine migrants grow
  at sea and fly to the forests for the dusk fruit fall and to mate. On land they leave
  nutrients as droppings, as carcasses and as the prey of glowing traps and pitcher
  plants, and the survivors fly back to breed at sea. The plants pay in sugar, which the
  long day makes in surplus. The migrants pay in the nitrogen and phosphorus the forest
  lacks, and since fruit carries little protein, they build their bodies at sea. Growth
  at sea and death on land give the largest net import, as with salmon.
- **H2. Earth has each piece.** *Earth evidence.*
  - **Seabirds:** colonies worldwide receive 591 Gg of nitrogen and 99 Gg of phosphorus
    a year in the birds' droppings (Otero et al. 2018). Near colonies on Spitsbergen,
    guano falls at a median of 0.3–0.4 g per m² a day (Zwolicki et al. 2013).
  - **Salmon:** trees and shrubs along salmon streams take 22–24% of their leaf nitrogen
    from the fish. Near spawning sites, trees reach 50 cm in diameter in about 86 years
    against 307 elsewhere (Helfield and Naiman 2001).
  - **Lake insects:** Lake Mývatn's midges emerge at 0.15–3.7 g of dry mass per m² of
    lake a year. In a high year they fertilise the first 50 m of shore with 1 kg of
    phosphorus per hectare (Dreyer et al. 2015).
  - **Bogong moths:** they fly about 1,000 km to Australia's alpine caves. By Green's
    estimate, about a billion die there each year, bringing 7.2 t of nitrogen and 0.97 t
    of phosphorus (Green 2011).
  - **High-flying migrants:** insects crossing southern Britain carry about 100 t of
    nitrogen and 10 t of phosphorus a year (Hu et al. 2016).
- **H3. Flight is cheap on the Moon.** *Computed · adopt.* For the same animal in the
  design air:
  - **Power and energy:** hovering takes 6% of the induced power it takes on Earth, and a
    kilometre of flight 17% of the energy at a fixed lift-to-drag ratio.
  - **Round trip:** at the median sea distance, 950 km there and back, a moth of Bogong
    size (0.33 g) burns a quarter of its mass in sugar. On Earth the same trip would
    cost one and a half times its mass.
  - **Speed:** flight is slower, at 0.38 of the Earth airspeed. The outward trip takes a
    moth three to five days and a 100 g bird one to one and a half.
  - **Keeping up with dusk:** dusk moves at 4.3 m/s at the equator. Small birds can keep
    pace with it (E6); moths only near the poles or with the wind.
- **H4. How large the conveyor must be.** *Computed · test first.*
  - **The target:** phosphorus returned as fast as Earth's rivers would carry it off in
    solution at the Moon's runoff, 0.003–0.008 g per m² of land a year.
  - **Bodies needed:** that takes 0.5–1.2 g of dry insect bodies left on each m² of land
    a year, or 0.13–0.32 g of vertebrate bodies whose bone makes them rich in phosphorus.
  - **What the seas must produce:** with seas on 28% of the Moon, they must send up
    1.2–3.0 g of dry insects per m² of sea a year, the top of Mývatn's range.
  - **Against Earth:** that is 5–12 times what all of Earth's seabird colonies deliver
    per m² of land, and 20–60 times southern Britain's insect migrants.
  - **Fuel:** the flight sugar is 0.05–0.7% of a wild forest's dusk fruit fall.
  - **Beyond its reach:** returning eroded phosphorus, or potassium at Earth's river
    rates, would take 21–89 g of insects per m² of sea. The conveyor suits dissolved
    phosphorus and nitrogen; potassium and eroded phosphorus need G4 and G6.
- **H5. Where the nutrients land.** *Earth evidence · a design point.* Animals deliver
  nutrients in patches:
  - guano deposition drops off within about 50–300 m of colonies (Zwolicki et al. 2013);
  - most midges land within 100–200 m of the shore (Dreyer et al. 2015);
  - salmon nitrogen reaches 25–100 m from the stream (Helfield and Naiman 2001).

  A lunar conveyor needs many landing places spread through the forests: roosts, lekking
  sites, fruiting groves, stands of pitcher plants and glowing traps, and the
  droppings-collecting plants of E5.
- **H6. Where carriers feed, and what they bring along.** *Earth evidence · prefer
  carriers that feed low.*
  - **Nitrogen and phosphorus hold roughly steady per gram of body:** insects are
    9–11% nitrogen (Wiesenborn 2011) and 0.6–1.3% phosphorus (Woods et al. 2004;
    Gratton et al. 2008). Bone makes fish 1–5% phosphorus, most of it in the skeleton
    (Sterner and George 2000; Hendrixson et al. 2007).
  - **Mercury builds up:** methylmercury rises about eightfold per step in the food
    chain (Lavoie et al. 2013).
  - **Selenium mostly doesn't:** across freshwater food webs its median magnification is
    1.01 (Pelletier et al. 2025).
  - **Carriers bring toxins with them:** pond sediments near Arctic seabird colonies
    hold 60 times the DDT and 25 times the mercury (Blais et al. 2005), and Bogong moths
    carry arsenic from lowland farms into their caves (Green 2008).

  Each step up a food chain leaves roughly a tenth as much biomass. Carriers that feed
  low in the sea's food web bring the most nitrogen and phosphorus with the fewest
  toxins. Vertebrate carriers add phosphorus through their bones.
- **H7. A life cycle one lunar month long.** *Idea.* Migrants leave the sea in the
  lunar afternoon, since the trip takes days, and arrive at dusk for the fruit fall.
  They feed and mate through the first nights while the lures glow, then fly back. Eggs
  laid at sea can hatch by dawn, when the sea's production restarts; the July bundle
  also proposes reproduction around dawn. A dusk-set timer (C2) and emergence timed by
  the Moon (C5) would govern it.

## Calculations

The figures quoted above, as the run computes them
([results/lunar_cycle_ecology.json](results/lunar_cycle_ecology.json)).

Flight, for the same animal in the design air against Earth's (1.43 against
1.23 kg/m³):

| Quantity | Moon ÷ Earth |
|---|---|
| Weight | 0.17 |
| Induced power to hover | 0.062 |
| Energy per km at a fixed lift-to-drag ratio | 0.17 |
| Minimum-power airspeed | 0.38 |

| Flier | Sugar per km, Moon / Earth | Airspeed on the Moon | One way to the median sea distance |
|---|---|---|---|
| 0.33 g moth (lift-to-drag 4, muscle efficiency 10%) | 0.086 / 0.52 mg | 1.1–1.9 m/s | 3–5 days |
| 100 g bird (lift-to-drag 10, muscle efficiency 18%) | 5.8 / 35 mg | 3.8–5.6 m/s | 1–1.5 days |

Distance from land, weighted by area, on the atlas grid (4 pixels per degree; the
seas are the four water bodies over 0.5% of the Moon; rain-fed lakes above sea level
would shorten the second row):

| To | Median | 90th percentile | Farthest | Land within 100 km | Within 500 km |
|---|---|---|---|---|---|
| The nearest sea | 476 km | 1,116 km | 1,735 km | 16% | 52% |
| Any standing water | 50 km | 280 km | 874 km | 72% | 96% |

Nutrients, per m² of land: the top metre's stock and what weathering releases
each year (basalt law, 1–10 times for glassy regolith):

| Soil | Potassium in the top metre | Phosphorus in the top metre | Potassium weathered a year | Phosphorus weathered a year |
|---|---|---|---|---|
| Mare (median of six soils) | 2.0 kg | 0.85 kg | 0.06–0.6 g | 0.03–0.3 g |
| Highland (Apollo 16) | 2.1 kg | 0.72 kg | 0.07–0.7 g | 0.02–0.2 g |
| KREEP-rich (Apollo 14) | 6.8 kg | 3.3 kg | 0.2–2.3 g | 0.1–1.1 g |
| Earth's upper crust | 34.9 kg | 0.98 kg | | |

Flows of phosphorus, per m² of land a year, for scale:

| Flow | Phosphorus |
|---|---|
| Earth's rivers at the Moon's runoff, in solution | 0.003–0.008 g |
| The same, as particles | 0.06–0.17 g |
| A 60% fruit harvest | 5.0 g |
| Earth's seabird colonies, spread over Earth's land | 0.0007 g |
| High-flying insects over southern Britain | 0.00014 g |
| Lake Mývatn's midges, first 50 m of shore, high year | 0.1 g |
| African dust on the Amazon | 0.0007–0.004 g |

The night and the lures:

| Quantity | Value |
|---|---|
| The plant's upkeep through the dark hours | 2.4 W/m² (0.70 idling) |
| Dusk fruit fall at a 30% fruit share, spread over the dark hours | 2.3 W/m² |
| Night respiration plus a whole 60% fruit fall, share of the O₂ above | 3 × 10⁻⁵ |
| The same, added to the CO₂ above if the air stays put | 1.5% |
| A 1 m² glowing surface at the fungal system's lowest efficiency | 0.04 W, under 2% of a m² of plants' night upkeep |

## Selection

A recommendation for the author, judged on five criteria:

1. It fits the environment the models compute.
2. Earth biology shows its parts exist.
3. It changes timing, stress tolerance or allocation before inventing a new metabolism,
   as the July bundle's engineering philosophy prefers (§18).
4. It solves a problem the models found: the night's carbon, the long day's sinks, the
   nutrients that flow one way.
5. Its risks can be bounded.

**Carry forward**, the core design:

1. **Evergreen plants that store and idle** on their own carbon balance, using the
   twilight, with their stores in wood and roots (B1, B2, C3, D1).
2. **Clocks for the lunar cycle:** a dusk-set night timer, daily clocks that leave
   light harvesting alone, and animal life cycles timed by dawn and dusk (C1, C2, C4,
   C5).
3. **The day fruit:** a program one sunlit half long, set at sunrise and ripe at
   sunset. On farms 50–60% of new tissue carbon goes to fruit, less in wild forests, and
   leaves are replaced at their own pace (D2–D4).
4. **Two landscapes:** plains of storage organs and forests of wood and fruit, with wet
   margins between (F1–F3).
5. **The dusk fruit fall and the warm night's food web.** The thing to watch is CO₂
   under still night air (E1, E2).
6. **Nutrients returned from settlements to farms,** the largest nutrient loop on
   farmed land (G3).
7. **Potassium from KREEP rock and the seas, erosion kept low and sediment recovered**
   (G2, G4, G6).
8. **The sea–forest migrant cycle** for nitrogen and phosphorus in forests near seas and
   lakes. Carriers feed low in the sea's food web, land in many places, and meet
   pitcher plants and droppings-collecting plants at the forest end (E4, E5, H1–H6).

**Test first:**

- leaves that survive about 240 dark hours, and metabolism that idles and restarts
  (B5);
- a large fruit's program compressed to one lunar day (D2);
- glowing lures that last through the first nights, with a catch the migrants can bear
  (E3);
- the conveyor's size: seas producing insects at the top of Mývatn's range, spread over
  hundreds of kilometres (H4).

**In places:**

- C4 plants on dry or hot ground (A3);
- Earth's 24-hour clocks, given daily cues, and cultivars that tolerate continuous
  light, on farms (C1);
- fruit picked young (D5);
- carnivorous plants on wet ground, on poor sands and in the canopy (E4);
- dusk followers (E6);
- dust and sea spray (G5).

**Set aside:**

- regrowing the canopy each cycle (B3);
- C4 or "hyper-C4" as the general strategy (A3);
- spectral tuning as a priority (A2);
- extra photoprotection, folding leaves and daily rests as design goals (A4, A5);
- frost protection at low and middle latitudes (B4);
- the night oxygen crisis (E2).

## Open questions and next work

- **The author's review.** The verdicts and the selection are recommendations; what
  enters the decisions register is the author's call.
- **Physiology to test:** B5, D2 and E3, and whether leaves use two weeks of continuous
  light (the canopy model's open question).
- **Models to build:**
  1. a land, lake and sea budget of nitrogen, phosphorus and potassium, over centuries
     to millions of years: weathering, leaching, erosion and sediment, the harvest and
     its return, dust, sea spray and the migrant conveyor;
  2. a migrant population tied to the fruit fall, the lures and the traps, to size the
     catch the carriers can bear;
  3. water and stomata in the plant model, which decides whether a daily rhythm pays
     (A5);
  4. CO₂ under still night air near the ground (E2).
- **Data:** weathering rates of lunar glass and soils; how much of their potassium and
  phosphorus plants can take up; erosion under lunar gravity and rain; iodine and
  selenium in lunar rock.
