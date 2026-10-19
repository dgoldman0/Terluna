# The aquatic programme: the waters' oxygen, methane, minerals and life

19 October 2026 · a plan, locked in on main at the author's word; no work has started, and it starts when the author
says ([register](../../decisions.md#the-design-climate-and-the-waters-19-october-2026)). It is the author's research
plan, drafted with ChatGPT and kept verbatim on the research drive (`session-tools/aquatic_plan_chatgpt_2026-10-19.md`),
with the changes agreed on 19 October ([below](#changes-from-the-draft)). The coupled study that draws the work together
will live in this folder; each model stays in its domain.

## The aim

A biologically rich, mainly fresh aquatic world at the design climate of 301 K: open water held nutrient-poor and rich
in species, productive habitats where nutrients gather, stable oxygen and methane, the dissolved minerals life needs,
reef and shellfish communities, and lasting harvests of biominerals. Living things take on much of the waters' mineral
and methane work, after much selective breeding: reef builders and other calcifiers take the calcium as calcite and, in
alkaline basins, the magnesium, their carbonate also serving as material; shellfish lock minerals in their shells; and
methanotrophs eat much of the methane.

Organisms start from Earth's, whose measured rates show what biology already reaches. Selection, breeding and
engineering shape them within the hard limits the register sets: light, carbon and energy supply, temperature, water,
transport and the 708.7-hour cycle ([register](../../decisions.md#life-and-people)). The work keeps four things apart:
what is ecologically plausible, what an organism has been shown to do, what a model predicts and what the design
requires.

## Where it starts

Checked against the repository on 19 October. These results are the reference cases; the work tests how far they hold,
what they mean for the waters' life and how they can be improved.

| Work | What it gives |
|---|---|
| [biosphere/ecology/lake.py](../../../biosphere/ecology/lake.py) | The water column's heat, mixing, oxygen and oxygen demand |
| [lake_oxygen](../lake_oxygen/README.md) | The lakes' oxygen in the three targets, under both forms of the deep mixing |
| [lake_kinds](../lake_kinds/README.md) | The lakes by kind at 301 K: the depth ladder, ventilation and phosphorus lines |
| [sea_oxygen](../sea_oxygen/README.md) | The seas' basins, stratification, renewal and oxygen |
| [lake_anoxia](../lake_anoxia/README.md) | Decay without oxygen, the floors' phosphorus release and the remedies |
| [lake_methane](../lake_methane/README.md) | Methane through the columns: production, three communities' oxidation, the anaerobic pathways, bubbles, the air |
| [sea_chemistry](../sea_chemistry/README.md) | Calcium, magnesium, alkalinity, salts and silica through a billion years, and the sinks that hold the seas |
| [sea phosphorus](../../../resources/sea_phosphorus.md) | Phosphorus delivered, kept by the lakes and passed to the seas |
| [food_systems](../food_systems/README.md) | Seaweed and shellfish farms and the waters' phosphorus ceiling |
| [resort_301k](../resort_301k/README.md) | Where the 301 K pieces meet: methane, the phosphorus room, the seas, dimming, watering |
| [the comparison](../climate_targets/comparison.md#the-lakes) | The lakes in the three target climates |

- **The climate.** 301 K, with the seas' median surface at 29.6 °C ([lake kinds](../lake_kinds/README.md)); the lock at
  295 K and the hot world near 307 K are the sensitivity cases.
- **The held fresh seas.** 0.142 g/kg at pH 8.32, with 0.84 mmol/kg of calcium, 1.68 meq/kg of alkalinity, 0.1 mmol/kg
  of silica and no magnesium, at three times calcite saturation. Four flows hold them: 150 Mt of calcite a year built by
  calcifiers; 95 Mt of magnesite made where the rivers enter the lakes, the managed route, or 182 Mt of magnesium clays
  at pH 9; 15 Mt of soda and sodium sulphate from works on the rivers; and 282 Mt of opal held by silica-building life.
  Left alone, the seas hold 0.42 mmol/kg of calcium and 0.78 of magnesium at year 10,000 and then turn alkaline
  ([the seas' make-up](../sea_chemistry/README.md)).
- **Methane at 301 K**, with Earth's methanotrophs at the methane share's centre: 0.57 Tg a year from nutrient-poor
  waters, 7.8 mesotrophic, 34 eutrophic and 57 as new reservoirs. Nutrient-poor waters hold the air at 1.01–1.75 ppm,
  under the held 1.87, with thin room; at the share's 0.8 they need 0.5–2.6 Tg a year of removal, and under Shih and
  Bouffard's form up to 0.77. Communities designed for sudden overturns bring mesotrophic waters to 4.9 Tg. In
  nutrient-poor water the shallow floors' bubbles carry most of the methane (the seas' water to 30 m, 2.2% of their
  10.6 million km², gives 0.25 of their 0.42 Tg a year); in productive water the giants' anoxic depths do, through
  oxyclines whose methanotrophs oxidise 88–99% of what reaches them. The anaerobic pathways are held to their electron
  acceptors' supply, and the ledger cuts every other source by biology's 90% while it takes the waters as traced
  ([the waters' methane](../lake_methane/README.md), [the resort at 301 K](../resort_301k/README.md)).
- **Oxygen.** In nutrient-poor water the deep lakes' deepest tenth keeps 5.7–6.3 mg/L by Osborn's rule and 3.5 under
  Shih and Bouffard's form; mesotrophic deep water runs low or anoxic under that form; the giants, 70% of the lake area,
  go anoxic at depth once mesotrophic, and the 3,500 m giant within thirteen years even nutrient-poor. The seas'
  area-weighted deepest tenth keeps 4.1 mg/L nutrient-poor
  ([the comparison](../climate_targets/comparison.md#the-lakes)).
- **Phosphorus.** The rain-fed lakes keep 96% of the land's phosphorus, and on the rivers alone the seas reach
  0.27–1.1 mg/m³, dust adding up to 18–28 at the top of its envelope. People's dissolved phosphorus has 57–71 g a person
  a year of room with 8 billion at the middle settling while the giants' floors stay oxic, and 5–10 g if their depths go
  anoxic ([the resort at 301 K](../resort_301k/README.md)). The waters' life can spare 21–56 kt a year, 1.4–3.7 kt of it
  in the seas, which holds sea farms to the animal part of the diet of 0.04–0.26 billion people
  ([food systems](../food_systems/README.md)).
- **The mixing.** The wind passes 0.3% (0.2–0.4%) of its energy to the stratified water. Osborn's rule against Shih and
  Bouffard's energetic regime is the largest physical unknown, moving the deep lakes' oxygen and methane two- to
  sevenfold.

## The questions

1. **Oxygen through the lunar cycle.** Which lakes and seas keep their oxygen from surface to floor, and where
   low-oxygen habitats persist: mixing across sills, the floors' demand and the deep water's renewal.
2. **Methane.** How much a mature aquatic biosphere makes, eats and releases, with production and oxidation counted
   apart along their physical paths and by community.
3. **Nutrient-poor water rich in life.** How the open water stays nutrient-poor while the waters hold high biodiversity:
   the habitats, recycling, delivery from the land and where production sits.
4. **The dissolved minerals.** The calcium, magnesium, alkalinity, silica and trace elements the waters hold for their
   life, with the salinity held.
5. **Living mineral sinks.** How far calcifiers, reef builders, shellfish, algae and microbes take the minerals out,
   with the carbonate chemistry, the habitat they build and the carbon.
6. **Harvests.** What pearls, nacre, carbonate stone and edible shellfish the waters can give lastingly, with
   production, finished yield, export and their effects counted apart.

Each is answered at 301 K, with the lock and the hot world as bounds where their results are cached or cheap.

## Who does what

The work divides by process, each part on its domain branch, with the coupled study in this folder.

- **Climate** (`domain/climate`): currents, stratification, the wind's mixing, overturning, density, gas exchange at the
  surface and transport between regions; the mixing form's check; the design case at 301 K.
- **Ecology** (`domain/ecology`): production, respiration, decay, the methanotrophs, the floors' communities and
  make-up, the calcifiers and the shellfish, and the organism evidence register.
- **Resources** (`domain/resources`): the dissolved and particulate balances of calcium, magnesium, carbonate,
  alkalinity, sodium, potassium, silica, phosphorus and the trace elements, with precipitation, burial, harvest,
  industry, dissolution and resupply.
- **Geography**, read on the branch that needs it: depths, hypsography, the basins' connections, rivers, groundwater,
  shores and sediment routing, taken from its products.
- **The coupled study**, in this folder at the next integration: the scenario matrix, the joint ledgers, the comparison
  and the recommendation.

Results cross as data products. The register and this plan change on main. Each branch is brought up to main whenever
main moves, by fast-forward while it has nothing of its own and by merging main once it has, so that the next
integration stays easy.

## The order

Agreed on 19 October, cheap work first.

1. **The lock-in** (main): the register's decisions and this plan. Done on 19 October.
2. **The seas' chemistry with magnesium kept** (resources). The sea chemistry with magnesium at 0.05–0.5 mmol/kg and
   calcium at 0.5–2 mmol/kg: salinity, pH, calcite and aragonite saturation, and what each sink must take (calcifiers,
   alkaline basins, the magnesite route and the clays), with the carbon each buries; the inventory each composition
   needs across the waters, the time to build it and the yearly throughput that holds it. Magnesium kept at
   0.05–0.5 mmol/kg raises the 0.14 g/kg to about 0.15–0.21 by the bicarbonate it brings, before the model's adjustment.
   The model runs in seconds.
3. **The calcifiers and the floors** (ecology).
   - *Habitat and rates.* Stonewort meadows in Earth's hard-water lakes lay down about 0.4–1.1 kg of carbonate a square
     metre in a growing season (Pełechaty et al. 2013; Pukacz et al. 2016). At that rate the seas' 150 Mt of calcite a
     year needs roughly 140,000–380,000 km² of meadow, against the 233,000 km² of the seas shallower than 30 m and the
     lakes' shallows. Microbial reefs (Bacalar's warm fresh lagoon; Salda and Alchichica, near pH 9, for magnesium) and
     shell beds join them.
   - *The night.* Photosynthetic calcifiers must carry 354 hours of dark at about 30 °C, which no Earth stonewort or
     microbial mat is adapted to; selection and engineering meet here.
   - *The floors and the methane share.* The share of the floors' decay that becomes methane, now a bracket of 0.3–0.8
     around 0.5, is tied to what the floors are made of: carbonate marl, the reactive iron of weathered basalt for
     iron-reducing microbes and iron-coupled methane oxidation, and the organic load. At 0.3 the floors of lakes to
     100 m stay below their bubble point and their methane falls ten to a thousand times; at 0.8 it rises two to seven
     times. The methane study reruns from its saved columns.
   - *Methanotroph guilds by habitat.* Aerobic candidates: *Methylocystis*, *Methylosinus*, *Methylobacter*,
     *Methylomonas* and *Methylosoma*, which the draft finds dominant through Lake Rotsee's methane-rich overturns;
     *Methylococcus capsulatus* Bath for warm, methane-rich water and treatment, its optimum above the seas' warmth; and
     *Methylocapsa gorgona*, which grows on the air's methane, for the soils. Without oxygen: *Candidatus*
     Methylomirabilis on nitrite and *Candidatus* Methanoperedens on nitrate, iron and manganese. For each: its
     temperatures, pH, oxygen, methane affinity, measured rates, growth and persistence, its nutrients and the copper
     its methane monooxygenase needs, its grazers, and its response to overturns and to the lunar cycle; whole
     communities beside single lineages; the N₂O of the anaerobic paths counted.
   - *The evidence register.* The measured values the models need, each in a sources file with its context and limits,
     bounded to those values; the draft's references, among them those dated 2026, checked as they are entered.
4. **The mixing form and the design case** (climate), each at the author's go-ahead with its machine time.
   - A lake turbulence model in columns (GOTM or similar), set to the Moon's gravity and rotation, the 1.2 atm air and
     the 301 K forcing, for one deep crater lake and one sea basin (Procellarum with its connected basins; Ingenii if
     cheap), checked against Earth's lake observations: the mixed layer through the cycle, mixing below it, overturns,
     the deep water's age and the oxygen's recovery, with the 0.2–0.4% wind share kept. A three-dimensional basin model
     (MITgcm or a lake and ocean package) follows only where the columns leave the two forms apart and the difference
     changes a conclusion.
   - The design case at 301 K ([the plan's item 1](../../plan.md#next-work-in-order)), then the glow pair from it and
     the lake studies on its weather.
5. **Harvests** (resources and ecology), once 2 and 3 stand: wild beds, selected cultivation and bulk mineral producers,
   each with its living biomass and area; its food, phosphorus, nitrogen and oxygen demand; gross deposition; finished
   yield by grade; rotation and recovery; recycling; the share in lasting use; and dissolution and return.
   - Phosphorus limits the shellfish: roughly a kilogram or more of shell grows for each gram of phosphorus in the
     flesh (a rough figure), so the food study's ceiling bounds them and their farms sit where nutrients gather.
   - Pearl mussels such as *Sinohyriopsis cumingii* lay nacre best at 28–30 °C by the draft's source, the seas' warmth;
     their larvae need host fish. Earth's freshwater pearl farms grow on fertilised water, so the Moon's keep their
     nutrients contained.
   - Whether the springs' water lays tufa and travertine, Earth's classic carbonate building stone, is checked against
     the springs study, which could give stone without taking it from living reefs.
   - The draft's ornamental targets, 2,450 and 24,500 t a year, are a hundred-thousandth and a ten-thousandth of the
     245 Mt of calcite and magnesite; at that scale the question is quality and habitat.
6. **The coupled study** (this folder, at the next integration): the methane, oxygen, nutrient, carbonate and mineral
   ledgers together; where living sinks meet the design's needs and where management stays; the comparison with the
   lock and the hot world; and proposed register entries. The design climate stays unless the coupled evidence gives a
   stated reason.

## Habitats

Five settings recur, each with its own chemistry and production, its methane and oxygen counted where they arise.

- **Open water**: little nutrient, modest plankton, efficient recycling and food webs on spread-out production.
- **Shelves and shallow reefs**: stonewort meadows, microbial reefs, algae, microbial films, shellfish and grazers,
  within the shelf's own oxygen and nutrient budget.
- **River mouths and deltas**: the fastest nutrient turnover and the most production, for intensive food webs and
  farms, with their sediment.
- **Deep lakes and basins**: surface and deep communities apart, methane eaten at the redox boundaries, and specialists
  where the oxygen stays low.
- **Alkaline and soda lakes**: the closed basins' own biome, with carbonate-forming microbes and the magnesium reefs.

The communities follow the register's condition of high biodiversity, which comes from variety of habitat and time
([register](../../decisions.md#the-seas-through-time-13-october-2026)). The world's map of production and trophic state
gives each lake and sea its own state, with its uncertainty. Each community is judged by its net methane and its whole
oxygen demand, since dense plant beds and shellfish can raise the floors' organic load: phosphorus held by oxic floors
and calcite and released when the floors turn reducing, uptake, grazing and remineralisation, the shellfish's filtering
and deposits, returns by migrating animals, and export to the sediment are counted together.

## Checks

- Closure for water, oxygen, methane carbon, dissolved inorganic carbon, calcium, magnesium, alkalinity, phosphorus,
  nitrogen, silica and the sediments' minerals.
- Methane oxidation uses its oxygen or other electron acceptor; bubbles leave the dissolved store; floors release
  phosphorus when their oxygen goes; dissolved minerals enter density and saturation; mineral-building removes the right
  ions and alkalinity; calcification, dissolution and lasting harvest are kept apart, with no uncounted sinks; nutrients
  moved between regions are conserved; the air's methane carries every source and sink; and settled states are told
  from transients.
- Calcification in bicarbonate water releases CO₂, so carbonate made and carbon taken from the air are separate
  quantities. Carbonate kept as stone keeps its carbon, and carbonate burned for lime returns it.
- Biology over many lunar cycles: growth, survival, reproduction, population stability, mineral-building and effects on
  the community. A short favourable measurement is a start; lasting performance needs the record over many cycles.
- Earth's reference systems, each with how it differs from the Moon: Lake Rotsee (methanotroph succession through
  overturn), Lake Chaohu (methane oxidation in a warm, shallow, eutrophic lake), Lake Lugano (nitrite-dependent methane
  oxidation in basins that mix differently), Lake Salda and Lake Alchichica (hydromagnesite and aragonite
  microbialites), Bacalar (warm fresh microbialites), the stonewort lakes and the freshwater pearl farms.
- The existing tests of code, input fingerprints, caches and conservation carry on and grow with the work.

## What the programme answers

1. Which lakes and seas keep their oxygen through repeated lunar cycles, and where low-oxygen habitats persist.
2. Whether living methane-eaters keep the waters' methane within the air's target, under unfavourable mixing and
   bubbling too.
3. Which nutrient conditions give the most biodiversity and lasting production with the least methane.
4. The calcium and magnesium the chosen organisms need, and the management that holds them.
5. Which calcifiers and shellfish best build reefs, pearls, lasting ornament and the bulk mineral sinks.
6. What production and harvest each habitat can sustain, and how much material stays in lasting use.
7. Which uncertainties move the answers most, and what measurement or simulation would settle them.

Recommendations keep apart supported Earth analogues, physically plausible extrapolations, selection and engineering
targets, and untested proposals. Where a result rests on a favourable assumption it carries the central estimate, a
less favourable bound and the assumption that sets the gap.

## Sources to read first

From the draft, each to be entered with its measurements before any calculation uses it, those dated 2026 checked
first:

- Mayr et al. (2020), [growth and rapid succession of methanotrophs limit methane release during lake
  overturn](https://doi.org/10.1038/s42003-020-0838-z)
- Mayr et al. (2020), [the mixing regime selects the methanotrophs' apparent oxidation
  kinetics](https://doi.org/10.5194/bg-17-4247-2020)
- [Aerobic methane oxidation in the water of a shallow eutrophic
  lake](https://www.sciencedirect.com/science/article/pii/S0301479726009023) (2026)
- [A widespread soil bacterium that oxidises atmospheric methane](https://doi.org/10.1073/pnas.1817812116) (2019)
- [Water column dynamics control nitrite-dependent anaerobic methane
  oxidation](https://doi.org/10.1038/s41396-023-01382-4) (2023)
- [Methane emissions from U.S. lakes dominated by ebullition](https://doi.org/10.1029/2026GL124792) (2026)
- Pukacz et al. (2016), [depth and the month in charophyte biomass and their calcium
  carbonate](https://doi.org/10.1007/s11356-016-7420-8)
- Pełechaty et al. (2013), [Chara vegetation in the precipitation of lacustrine calcium
  carbonate](https://doi.org/10.1111/sed.12020)
- [Charophytes as drivers of carbonate sedimentation, carbon sequestration and aquatic
  health](https://www.frontiersin.org/journals/freshwater-science/articles/10.3389/ffwsc.2025.1520590/full) (2025)
- Gunes et al. (2024), [hydromagnesite microbialites in Lake Salda](https://doi.org/10.1111/gbi.12619)
- [Carbonate–microbe interactions in Lake Alchichica's microbialites](https://www.nature.com/articles/ismej201381)
  (2013)
- [Temperature, pearl deposition and matrix protein genes in *Hyriopsis
  cumingii*](https://doi.org/10.1016/j.aaf.2021.06.006) (2023)
- [Algal symbionts of the freshwater sponge *Ephydatia muelleri*](https://doi.org/10.1007/s13199-023-00934-8) (2023)
- [Seasonal O₂, DIC and CH₄ fluxes in sediments with *Vallisneria
  spiralis*](https://doi.org/10.1016/j.aquabot.2011.01.003) (2011)
- [The nutrient–methane nexus mediated by submerged macrophytes](https://doi.org/10.1016/j.watres.2024.122316) (2024)

## Changes from the draft

- Its figures were checked against the repository and hold. Its first stage, the audit of main and the branches, was
  done the same day, every branch level with main.
- Its Earth-origin organisms and breeding sit within the register's engineered organisms: Earth's organisms are the
  starting stock and the evidence of reachable rates, and engineering stays available where Earth's range ends.
- With 301 K the design climate, its scenario matrix runs at 301 K, with the lock and the hot world where cached or
  cheap.
- The mixing starts with turbulence columns at the Moon's values; a basin model follows only where they leave the two
  forms apart.
- The food study's phosphorus ceiling on sea farms bounds the shellfish.
- Its ornamental targets are its own figures, small beside the bulk flux.
- The literature is bounded to the values the models need.
- Added: magnesium as the harder half, with its reefs in alkaline basins; the floors' make-up as the methane lever; the
  night as the photosynthetic calcifiers' hard limit; and the methane study's three communities and anaerobic pathways
  as the base the methane work extends.
