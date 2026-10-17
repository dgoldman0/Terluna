# Ecology of the Open Moon

The engineered organisms of the Open Moon and the communities they form: biomes and food webs, soils, the cycles
of nitrogen and the trace gases, the light and darkness life lives by, the waters, fire, dispersal and containment.
The work opened on 9 October 2026 on the branch `domain/ecology`, inside the biosphere domain, which owns biology
and ecology together. It starts from three things:
- what the project already holds ([state.md](state.md));
- the literature across ecology, plant and animal physiology, microbiology, photobiology and atmospheric chemistry
  ([literature.md](literature.md), with its [sources](sources.json));
- a first screen of the committed products ([screen.py](screen.py), [results/first_screen.json](results/first_screen.json)).

Its open questions, in order, are in [questions.md](questions.md).

The author's standing directions frame it.
- Organisms are engineered, with Earth species' traits as reference points and design variables.
- Surface materials become developed, weathered and biologically worked soils.
- The Moon should hold a mix of climates and biomes, with pleasantly warm days at the hottest.
- The ring fleet keeps its light off the night zones.
- The aerial biomes come almost first: the skies are a large set of biomes of their own and connect the others.
- Imports of carbon and other elements are part of the plan, and the design problem is closing each cycle with
  local return paths.

## The 9 October follow-ups

**The [coupled aerophyte study](../../research/studies/aerophytes/README.md) sets how large aerophytes grow, how
old they get, where they find water and how they live through storms** (*screen* values, product schema
`terluna.research.aerophytes/1`).
- **Size.** At 10 km, with joints, seams, storm gusts, spare pressure for trim, wet fire-barrier skins, gas cells and
  tendon renewal paid, a single round body lives from about 50 m to 700 m across with cellulose-class tendons and to
  1 km with 300 MPa fibre, and to 1.5 and 3 km with those fibres in canopy-like light. A colony of 60–150 m modules
  has no structural size limit.
- **Age.** A sky reef grows to 2 km across in 36–544 years and a round giant to 1 km in 73–1,186 years.
- **Water.** Mean rain meets the reference water need equatorward of 35–38°, and to 51° with cool leaves and a
  CO₂-concentrating trait. Sorbent skins take 0.2–2.4 kg/m² a day of vapour from the 64% air.
- **Storms and sailing.** Storms are escaped by route and hour, lightning burns small wet gas cells one at a time, and
  weighted tethered wings sail aerophytes 30–70° of latitude in a lunar day.
- **Nursery.** Megaforest crowns can raise juveniles to 50–60 m in one to seven years.
- **Open:** assembled wet materials, paths through the latitude and solar-hour winds, and recruitment.

**The [aerial food-web evaluation](../../research/studies/aerial_ecology/README.md) treats aeroplankton and large
photosynthetic aerophytes as producers in their own right**, with reproduction against removal, filter drag, wet mass,
optical depth, nutrients and harvest accounted separately. Its central case, aerophytes over 1% of the Moon at
1,000 g C/m² a year, gives 184.6 million food-energy equivalents before edibility and processing losses, so the sky's
wild harvest is a small food source and is credited with none.

**The [UV-window evaluation](../../research/studies/uv_biology/README.md) keeps the present shield,** with
biological adaptation and local UV as the working direction and a narrow UV-B window open as a candidate. The current
film also suppresses most short UV-A. The surface doses are computed; the window's chemistry and exobase effects are
open.

**The population figures below give what each billion people takes.** How many people the Moon's settings hold, on
the ground and aloft, comes from provisioning's comparison in the [shared population
cases](../../shared/scenarios/population.json). Provisioning owns the inhabited-volume and ways-of-living studies and
resources the phosphorus-return model.

## A taxonomic system, 9 October 2026

The author asked on 9 October to rename the aerophytes and to begin a taxonomic system for the biosphere.
[taxonomy.md](taxonomy.md) proposes one:
- **Descent.** A designed organism keeps the Earth lineage of its chassis, the lineage whose genome organizes the cell.
- **Degree of design.** A line within an Earth species, a designed species, a designed genus or a new lineage sets the
  form of its name.
- **Founding.** Every founding lineage is registered with a type genome and a preserved sample.
- **Descendants.** They are named as clades defined from the founders, which keep their meaning over the
  billion-year horizon.
- **Function.** Life form, setting, guild, night strategy and clock are classed apart from descent.

The first register, written by [taxa.py](taxa.py) to [results/taxa.json](results/taxa.json), holds 42 designed groups
with 11 subgroups, 31 design modules and 14 communities drawn from the project's studies, the sky reefs' tenants,
domatia, pools and distributed control among them. Each comes with its Earth reference organisms, placement, status and
sources. The author chose the everyday names on 10 October: aerophytes, with sky reefs for the colony giants. On 11
October the author named the giants *Welkinia*; the smaller kinds take names suited to their size and how they live, and
the class's scientific name and six choices about the system wait for the author.

## The living sky, 10 October 2026

The author asked how aerophytes join the aerial and global ecosystem: whether the small ones are eaten in the food
chain, and how much life the larger ones carry before they sink. The [sky-ecology
study](../../research/studies/sky_ecology/README.md) answers in one coupled screen (product schema
`terluna.research.sky-ecology/1`), checked by an independent reviewer in three passes and a confirmation:
- **Food web.** Over 1% of the Moon aerophytes build 162–413 Mt C of biomass a year while they grow and 134–397 Mt
  once they hold their size, and give their tenants 38 Mt C of sugar. On host food alone the tenants that eat the host
  reach 8–41% of the production their carbon allows, short of phosphorus, and the sky feeds 1.8–4.8 million
  albatross-scale grazing flyers. Small aerophytes are costly food: 5–8% of what they cost is edible tissue.
- **Tenants.** A lone body's underside gets 19–30% of the light on its top, and a 100-m round body can carry 196 t of
  tenants within its spare lift. A load is sudden when it lands within a few minutes, against a buoyancy period of 25
  minutes: left alone a body survives a step of 41–50% of its trim equivalent and a gradual load of 81–100%, venting
  takes 71% of it at once and dropping water all of it. Sky reefs can grow domatia for their tenants.
- **Phosphorus.** Aerophyte tissue holds 0.7–8% of Redfield's phosphorus. With pools that catch litter as a forest
  floor does, a mature sky needs 3.5–61 kt of new phosphorus a year at 1% cover and a growing one 13–165 kt; with
  nothing to catch it, 62–424 kt and 71–528 kt. Roosting sea-feeding flyers close a giant's phosphorus and can return
  the land's dissolved loss, and a fallen 2-km reef lays 0.5–5.6 g P/m² on the ground it drapes, the land's whale fall.
- **Demography.** Replacing the dead costs round giants up to two-thirds of their surplus, and fermented giants at the
  high hazard more than all of it. At a 200-year mean adult life 0.3–32% of 1-km round giants reach full size at the
  reference light, and 7–66% of 2-km sky reefs.

**Guards and partners.** At the author's direction the [guards and symbioses
analysis](../../research/studies/sky_ecology/symbioses.md) gives a reef's defence to several species sharing the work:
small patrols near green tissue and nests, hunters in sheltered margins, flying defenders around arrivals, recycling
communities in wet pockets, and visitors that feed elsewhere and roost on the reef. Response time sets the nests:
arrival within five minutes needs nests no more than 5.1 m apart at 0.02 m/s, or 25.5 m at 0.1 m/s, after a guessed
30-second alarm. Ten grams of live guards per square metre of reef take a quarter of the host's food for tenants, and a
hundred grams more than all of it.

**Living sails.** The author let a sky reef be partly animal-like, with a distributed control like a nervous system
and no brain. The [living-sail note](../../research/studies/aerophytes/sail_biology.md) builds the sails on that: a
lasting support system carries replaceable sails, each with its own sensing and contractile control under slower
coordination from the colony. A signal at the 0.3 m/s measured in a siphonophore's swimming-bell epithelium takes 4.6
hours to run 5 km of cable, so each sail steers itself, and a water column 5 km tall holds 8 MPa, so the cable carries
pressure-isolated service nodes. The earlier sailing speeds assume a wing that holds its angle; one specified living
sail, tested for stability and for failure, decides how far a reef steers.

Open: a population model on the routes, measured tissue phosphorus, the seas' food web, great flyers' energetics at
the lunar scale, and that one living sail.

## Living light, 11 October 2026

The author asked to rely more on bioluminescence and to let lighting follow the Moon's twilit nights. The [living-light
study](../../research/studies/living_light/README.md) gives what glowing organisms can do, what people need outdoors in
each setting, and how far people's light reaches the sky:
- **The glow.** Today's brightest glowing plants mark a path in the darkest hours and give a hundredth of the light
  needed to walk by. A designed plant with a beetle-class luciferin cycle, its glow facing down, lights a path at
  walking level for 2–10% of its night store; the fungal pathway as measured costs a third of the store to thirty
  times it, so its quantum yield is the gap.
- **The natural night.** Twilight and Earthlight keep the darkest moment above 1 lux, enough to walk by, on 48% of the
  dry land and above 8 lux for reading on 21%. Poleward of 48.6° the night keeps walking light all night; the far
  side's low latitudes are the dark belt, below 1 lux for 165 hours a cycle at the summit.
- **People's light in the sky.** The tall, dense air carries a town's skyglow farther than Earth's air does, yet the
  bright natural sky leaves ordinary towns' skies unpolluted on two-thirds of the dry land. The land under a
  light-polluted sky falls from the earlier 14–49% of the dry land to 0.9–40% at 4–8 billion lunar residents, and
  living light taking over 30–70% of outdoor lighting brings the low end to 0.5–0.6%.
- **The dark belt** holds Earth-like darkness on 13% of the dry land and is where luminous organisms signal; near
  towns, direct light matters more to night life than skyglow.

## Machines in the living world, 12 October 2026

The author set the wild under heavy automation, with the machines part of the whole world's technobiosphere, and chose
every way of joining the two that was offered
([register](../../research/decisions.md#machines-in-the-living-world-12-october-2026); [machines.md](machines.md)):
- **The machines' work.** Sensing everywhere, constant evaluation and hands-on work where the evaluation calls for it;
  automatons that behave as members of animal groups or as their predators; a shared signal vocabulary that sits mostly
  outside human senses; reporting traits and DNA watermarks; machines in ecological roles; field machines that compost
  when lost; and handling, killing included, under conservation's ethical standards and people's constant oversight.
- **What may live wild.** General constraints for an open world, with no containment, no hard limits and regions kept
  distinct loosely: climate fit, harmless Moon-wide, within the air's sinks, closed cycles, ready for the machines and
  robust to evolution.
- **N2O and methane.** Biology first, with technologies making up the difference, and no ultraviolet window in the
  shield. N2O comes mostly from the wild's soils and waters, 0.3–0.4 ppb a year, and reaches Earth's present level in
  about a thousand years. Nitrogen fixation feeds the nitrogen it comes from, and breaking it down with lamps in open
  air would take hundreds of terawatts. In the deep column both gases warm the Moon strongly: Earth's present levels
  force 4.7 W/m² all-sky and the screen's central steady state 7.9–8.3 W/m², 4.9–8.3 K before it is answered
  ([questions](questions.md), item 3). The author chose to hold them at 1.9 ppm of methane and 0.1 ppm of N2O, the
  balanced middle of three pairs weighed: 3.58 W/m², answered by 1.69% added dimming, with biology cutting most of
  each source and removal from open air taking the rest.

## Fresh shores and soda lakes, 15 October 2026

The author asked how lush the shores of fresh seas could be with lunar gravity's help, and then whether the load that
gathers where plants draw groundwater could make salt lakes. The water's side is geography's
[shore screen](../../geography/README.md#shores-by-the-fresh-seas) and the sea chemistry's
[closed lakes](../../research/studies/sea_chemistry/README.md#closed-lakes-soda-lakes-and-salt-seas), and Jura Island's
numbers are in [its study](../../research/studies/jura_island/README.md); these are the living world's side. All are
screens.
- **Fresh seas make green shores under desert skies.** The held seas, 0.14 g/kg, feed the water table under every shore
  with fresh water: no salt wedge rises under the shore when plants draw on it, and the spray waters what grows by the
  sea. Earth's analogues are the shores of its great lakes and the tidal reaches of its big rivers above the salt. The
  monthly tide, 3.7 m typical on the Near-side Sea, floods and drains a band of every shore once a month, a freshwater
  tidal marsh on the month's beat (finding 18).
- **Lunar gravity helps the roots and slows the ground.** Soil lifts a few millimetres a day to roots from about four
  times Earth's height above the water table, 12 m against 3 m in a sandy loam, and a 50 m tree holds its water column
  up with 0.08 MPa against Earth's 0.49. The sea's water moves sideways through the same rock six times more slowly.
  On steep shores the belts reach 15–25 m above the sea, against 6–10 m on Earth, and so run two to three times as
  wide; on flat shores the ground's supply holds them to 0.6–2.5 km, about Earth's width.
- **Small low islands can be wooded throughout.** The sea lies under the whole of them. Up to about a kilometre across
  in ordinary ground, or a few kilometres in sand, groundwater alone keeps them green from shore to shore. Where the
  draw ends, in the middle, the sea's load settles as a hardpan of limestone, magnesium minerals and opal, cementing
  about 18 cm of soil a millennium under a metre a year's draw. Salt builds slowly, since sodium and potassium are 2%
  of the load. The Okavango Delta's islands, whose water is fresher but carries more sodium, show the end state left
  alone: wooded rims round bare, salty centres (McCarthy and Ellery 1994; *recalled*). Rain, the tide or flushing keep
  the middles green.
- **High islands catch less cloud than Earth's.** Rising air cools by 1.6 K a kilometre at lunar gravity against
  Earth's 9.8, so for the same humidity cloud bases stand six times higher, and a 2.7 km island lifts air as a 450 m
  hill would on Earth. A high island's green comes mostly from the sea beneath it. Jura Island's box will show whether
  its arc catches cloud. *Derived; to test.*
- **Closed lakes become soda lakes, a biome of the dry belts.** The Moon's water carries calcium, magnesium and a
  little sodium with carbonate and almost no chloride, and rain from fresh seas brings none, so a lake with no outlet
  loses its calcium, magnesium and silica to its floor and keeps sodium carbonate. One 10 m deep reaches 1 g/kg in a
  few thousand years and the ocean's salinity, as soda at a pH near 10, in about 100,000, and the soda works could
  make one in years. On Earth soda lakes are among the most productive waters, with a specialist life: Arthrospira
  feeding lesser flamingos at Bogoria and Natron, brine shrimp and alkali flies feeding millions of migrating birds at
  Mono Lake, a tilapia at pH 10 in Magadi and the pearl mullet in Lake Van (*Earth evidence, recalled*). The author
  takes them as their own local ecosystems, a biome that may be fairly common on the Moon
  ([register](../../research/decisions.md#soda-lakes-15-october-2026)).
- **Saltwater life needs its chlorine brought.** Ocean-like seas take 36 Mt of salt per km³, 55% of it chloride that
  lunar rock cannot give: 2,600 Gt for Moscoviense, the smallest separate sea, and 12,600 for Smythii–Marginis. The
  author leans to aquariums, where a million cubic metres takes 36,000 t (open, resources' register).
- **The climate runs' land hides all of it.** Every run's land is bare placeholder ground at one dryness, with no water
  table and no plants, so no grid shows a green shore. The cover has to come from the plants and the groundwater
  ([questions](questions.md), items 6, 9 and 10).

## What the project already says

Each finding below follows from results the repository already holds, read together, with the literature that
bears on it. Numbers marked *screen* come from [results/first_screen.json](results/first_screen.json); the rest name
their source. The evidence ranges from solid physics to Earth analogues carried to the Moon, and each finding says
which.

### The light under the shield

1. **No vertebrate on the Open Moon makes vitamin D from sunlight.** The titania stack passes nothing below 320 nm.
   With the Sun overhead the ground gets no UV-B photons, against Earth's 6.8 µmol m⁻² s⁻¹, and the band that makes
   previtamin D carries 2.5×10⁻⁷ of Earth's light (*screen*). The previtamin D3 action spectrum lies below 315 nm
   (MacLaughlin et al. 1982; CIE 174:2006), and the revised in-vivo spectrum lies 5 nm shorter still (Young et al.
   2021). Earth's food webs make their vitamin D with UV-B somewhere in the chain. Here every vertebrate needs another
   route: engineered biosynthesis that runs without UV-B, calcium handling that does without the vitamin, or vitamin
   D supplied by people. Cats and dogs show that mammals living on dietary vitamin D are viable; basking reptiles are
   the hardest case. People need supplements or UV-B lamps. *Solid for the light and the chemistry.*
2. **Bees and birds lose most of their ultraviolet.** A honeybee's UV receptor (peak 344 nm, Peitsch et al. 1992)
   catches 6.9% of Earth's light with the Sun overhead and 6.4% at 30°, while its blue and green receptors catch 49%
   and 66%. Its UV signal against green is therefore a tenth of Earth's. A bird's UV-sensitive cone (370 nm) catches
   22% (*screen*). What UV remains lies at 350–400 nm. Flower UV guides and the UV skylight insects use fade.
   Engineered flowers and pollinators can signal at 350–400 nm or in the violet and visible; UV patterning in flowers
   follows the UV a habitat receives (Koski & Ashman 2015). Bees read the polarized skylight at moonlight levels
   (Dacke et al. 2003), so polarization navigation may survive on the light that is left. *Solid for the light; the
   behaviour is to test.*
3. **Sunlight barely disinfects.** Inactivation spectra falling a decade per 30–50 nm give the ground 1–6% of
   Earth's disinfecting dose (*screen*). Solar disinfection of water needs more than 108 kJ m⁻² at 295–385 nm
   (Ubomba-Jaswa et al. 2009). Pathogens therefore last longer in lakes, on wet leaves and in the air. CM1's land
   nights are 97–99% humid at 19–21 °C with fog and dew, leaf wetness of a hundred hours and more, and Earth's leaf
   fungi infect within 6–24 hours of wetness. Disease resistance becomes a design requirement of the evergreen
   strategy. *Solid in direction; the size depends on each pathogen's action spectrum.*
4. **Plants lose their UV-B cue, though only partly.** UVR8, the UV-B receptor, responds up to about 335–350 nm (Rai
   et al. 2020), where the stack passes 0.2–5%. Crops under UV-blocking greenhouse films grow as well or better, carry
   15–35% less flavonoid and anthocyanin, and see less pest pressure (Katsoulas et al. 2020). UV-B also primes the
   jasmonate defences against herbivores. The Moon's plants would grow well, colour paler and defend less. *Earth
   evidence; the size of the remaining UVR8 signal is uncomputed.*
5. **Nearside nights are brighter than any nocturnal threshold Earth's life knows.** Under a high Earth, Earthlight
   keeps the nearside at 2.3–2.9 lux at its darkest (decisions register), 10–60 times typical full moonlight on
   Earth, 0.05–0.3 lux (Kyba et al. 2017), and about a hundred times the 0.01–0.03 lux at which light at night
   suppresses melatonin in fishes and rodents (Grubisic et al. 2019). Toward the limbs it falls to 0.07–0.44 lux
   (joint synthesis, corrected on main on 9 October), still above those thresholds. Earth-calibrated nocturnal life
   fits the far side's low latitudes, about 0.001 lux at midnight, which is also where the summit metropolis stands. The fleet decision keeps the night
   dark; nearside organisms still need raised circadian thresholds. The night's requirement can be set by zone.
   Earth's lunar-cued behaviour, such as coral spawning in the dark hour after sunset (Lin et al. 2021), has to be
   rebuilt on the Moon's own month, where the 29.53-day cycle is the day itself. *Solid for the light.*

### An air without its cleanser

6. **Methane lasts centuries, and a wet Moon makes plenty of it.** Under the stack the air holds no OH (the joint
   synthesis), so methane leaves only into soils. Soils taking it up at Earth's rate per area over the 62% of the
   Moon that is dry land remove 1.1–4.8 Tg a year for each ppm in the air, a lifetime of 370–1,640 years (*screen*).
   Lakes on 10.3% of the Moon emitting at Earth's 2–30 g m⁻² a year, with 0.5–2 million km² of wet margins at
   10–50 g m⁻², give 13–218 Tg a year. Methane then settles near 20 ppm at central values (59 Tg a year), and
   anywhere from 3 to 200 ppm across the ranges, approached over several centuries (*screen*). Deep crater lakes
   (median 352 m) emit less per area than Earth's shallow ponds. Bubbles leave sediment more easily at a sixth of
   Earth's water pressure, though, and skip oxidation in the water column. Methane oxidizers in soils and waters are
   therefore a requirement from the start. Wetland extent, lake depth and anoxic volume become variables of the air's
   composition (Bastviken et al. 2011). In this deep column 21.6–29 ppm force 5.5–5.9 W/m² all-sky, 3.4–5.9 K
   ([questions](questions.md), item 3).
   *Solid in direction; the rates are Earth analogues, good to a factor of three to ten.*
7. **N2O has no sink at all.** N2O is destroyed by light below about 240 nm and by excited oxygen atoms made by
   ultraviolet, and the stack removes both. Soils take up at most 2% of Earth's N2O sources (Schlesinger 2013).
   Earth-like soils and seas would emit 1.0–1.3 Tg N a year, raising N2O by 0.3–0.4 ppb a year: 0.2 ppm over the
   500-year construction and 3–4 ppm in 10,000 years. The light passing the swarm's gaps gives a lifetime of 45,000
   years to 4 million years (*screen*). The biomes' nitrogen cycle has to end in N2, so complete denitrification,
   with the N2O reductase NosZ, is a design requirement. *Solid chemistry; the sources are Earth analogues.*
8. **Ethylene may decide whether cereals work.** Crop stands give off ethylene; closed chambers reach 40–120 ppb
   (Wheeler et al. 1996), and 50 ppb held continuously costs wheat 36% and rice 63% of their yield (Klassen & Bugbee
   2002). Without OH and ozone its sinks are soils and rain. The same holds for isoprene from trees and dimethyl
   sulfide from seas: seawater at 1–10 nM of dimethyl sulfide would hold the air near 2–17 ppb at equilibrium, around
   the human odour threshold (digest arithmetic, [state.md](state.md)). The aerosol study's sulfate and organic
   particles assumed oxidation that the joint synthesis later found absent. *Uncomputed; the hazard is Earth
   evidence.*
9. **"No oxidant" holds for the bare column.** Soils emitting NO at Earth's natural rate would put 15 to several
   thousand times lightning's NO into the air. NO2 still splits at a third of Earth's rate, because 350–420 nm
   passes, which can make some ozone, and ozone with alkenes or nitrous acid in light can make some OH. A living,
   peopled column needs its own chemistry run, with biological and urban emissions and their sinks (plan item 3).
   *The ratio is solid; the chemistry is open.*

### A world of long reach

10. **Pollen and spores stay aloft for weeks.** Particles fall at a sixth of Earth's speed. A 25 µm pollen grain
    falls at 3.1 mm/s (Earth 19), so a 5 km daytime mixed layer holds it for 19 days and CM1's 12 km noon layer for 45,
    against under a day on Earth. In a light 3 m/s wind 19 days carries it about 4,800 km. A 10 µm spore stays
    100–250 days, unless rain removes it first, and rain falls mostly in the tropics. Plumed seeds fall at 0.37 of
    Earth's speed and reach 2.7 times as far (*screen*). Founding biomes spread into one another; isolation distances
    for engineered plants fail; pathogens and allergens travel Moon-wide. Transgenic pollen already reaches 21 km on
    Earth (Watrud et al. 2004). Synthetic auxotrophy resists both mutation and gene transfer (Mandell et al. 2015),
    but it needs a supplied amino acid, which suits crops and bioreactors. For wild engineered organisms no
    demonstrated firewall exists, so what may live wild is an ethical and governance decision. *Solid physics.*
11. **Flight is cheap for animals too.** Hovering takes about a sixteenth of Earth's power per unit mass, and a
    kilometre on a wing a sixth of the energy (the ecology register and the sky-fleet study). Large flying grazers and
    long-range pollinators are natural on this world. Insects gain the aerodynamics, and diffusion through their
    tracheae at 1.2 atm runs at 0.83 of Earth's at the same oxygen partial pressure, which works against giant
    insects. *Solid.*

### Nutrients and soils

12. **Molybdenum may limit nitrogen fixation before nitrogen does.** Biology must fix almost all of the biomes' new
    nitrogen, since lightning gives at most 0.2% (the joint synthesis, confirmed from Schumann & Huntrieser 2007).
    Nitrogenase needs molybdenum, or vanadium or iron in its alternative forms. Molybdenum limits free-living fixation
    in weathered soils (Barron et al. 2009), and the bulk silicate Moon holds about 0.019 ppm (Sossi et al. 2024).
    Molybdenum and vanadium belong on the import list beside carbon and potassium, and the alternative nitrogenases
    belong in the engineered fixers. Fixation itself is affordable: Earth's cost of 6–12 g C per g N fixed is 0.2–2%
    of the stand's growth over the long day (digest arithmetic). *Literature cross-link; lunar molybdenum data are
    sparse.*
13. **Soil is the slow variable, and fresh basalt runs short of nitrogen first and phosphorus later.** Natural soil
    forms at 0.017–0.036 mm a year (Montgomery 2007), about 8,000–18,000 years for 30 cm, and lasting humus builds at
    0.7% of plant production (Schlesinger 1990). On Hawaii's basalt, young soils are short of nitrogen and old ones of
    phosphorus, and after 4 million years the forests live on phosphorus in dust from 6,000 km away (Chadwick et al.
    1999). The Moon has no continental dust, so its long-term phosphorus supply is a designed return. Pioneers on
    fresh basalt are nitrogen fixers and chemolithotrophs (Kelly et al. 2014). Biological crusts fix nearly half of
    Earth's land nitrogen (Elbert et al. 2012), and a moss mat held Hekla's succession for 170–700 years
    (Vilmundardóttir et al. 2018), so the order of introduction decides which state a region reaches. Microbes
    weather basalt strongly at any gravity (Olsson-Francis et al. 2012; Cockell et al. 2020). Soils are built by
    organic loading and fine material. *Earth evidence.*
14. **Long daylight breaks litter down without the soil.** Light breaks plant litter down through lignin's absorption
    in the UV and visible (Austin & Ballaré 2010). The stack removes the UV-B share, while 354 hours of daylight
    lengthen the exposure. Dry biomes, 42% of the land under 0.5 mm of rain a day, could lose litter carbon to the air
    without it entering the soil. *Direction from Earth evidence.*

### Climate, biomes and fire

15. **The biome mix comes from water and light.** The GCM's coldest 3-day mean anywhere is 14 °C, and the highest
    land is about 10 K cooler than the lowlands. No frost, alpine, boreal or seasonal biome exists on natural ground.
    The variety the author asked for on 2026-09-25 ("definitely more of a mix of climates and biomes") lies in:
    - rain: wet tropics at 4.8 mm a day, a fog desert at 30–60°, rainless poles;
    - light: the Earthlit nearside, the dark far-side equator, and poles that never fall below about
      1 µmol m⁻² s⁻¹ of light poleward of about 65°;
    - lakes, tides and soils.

    Cold places are the built ones, the winter parks at 13.5–15 km. The only year-like signals are the Earth's orbit,
    ±3.4% in sunlight, and the tides' beat, 347–412 days. *Solid on the models' climate.*
16. **Fire: ignition is rare and spread is uncertain.** Natural ignition is about a thousandth of Earth's
    (0.0016 ground strikes per km² a year from the joint synthesis). At 17.5% oxygen, tests at Earth's gravity put
    fire in the suppressed band, greatly suppressed below 18.5% and off below 16% (Belcher et al. 2010). Lunar gravity
    lowers materials' oxygen limit by 2–6 points (Ferkul & Olson 2011, in the port fire study), so natural fuels may
    burn here. No one has measured them at 17.5%, 1.2 atm and 0.16 g. Earth's fire-kept grassy biomes, over 70% of
    its burnt area (Bond & Zaloumis 2016), close over into forest where climate allows, unless people graze or burn
    them. *Ignition is solid; spread is open.*

### Waters

17. **The first seas are fresh, mixed each month, and eutrophic at first.** Weathering alone takes 1,100–87,000 years
    to bring the seas to 1 g/kg (the resources screen), so unless the delivered water carries salts the living seas
    begin as fresh to hard waters. Ice delivered fresh would leave them near pH 5.6 and undersaturated in carbonate,
    so calcifiers cannot build until weathering supplies alkalinity. New reservoirs run eutrophic for a decade
    (Paterson et al. 2019). The 354-hour heating and cooling give monthly stratification and overturn, which keeps
    waters oxygenated. Deep crater lakes, with a median depth of 352 m, risk anoxia, and productive shallows run short
    of oxygen late in the night at the night winds' gas exchange. The living-water guesses take Earth seawater at
    38.4 g/kg. *Solid if the water arrives fresh; its composition is open.*
18. **The monthly tide makes flats exposed for days.** Ranges of 2.5–3.7 m on shores sloping 1:50–1:100 expose about
    20,000–30,000 km² for days each month, the most productive kind of coast on Earth and a methane and N2O source.
    Low water drifts through the lunar day on the tidal lines' beat of about a year (digest arithmetic). *Rough.*

### Gravity and the long cycle

19. **Plants sense lunar gravity badly.** The threshold for plant gravity sensing lies at 0.1–0.3 g (Kiss 2014), and
    simulated Moon gravity disturbed root meristems more than simulated weightlessness (Manzano et al. 2018). Root
    growth regulation is an engineering target. Mice at lunar gravity keep their soleus muscle mass (2023).
    *Earth-based simulation.*
20. **Retuning crop clocks to the 708.7-hour cycle is probably the highest-value change.** Continuous-light injury in
    tomato comes from circadian asynchrony (Velez-Ramirez et al. 2017), and one locus gives tolerance with up to 20%
    more yield (Velez-Ramirez et al. 2014). Plants run their starch out about 24 hours after the last dawn whatever
    the real dawn (Graf et al. 2010). The canopy model's idling strategies already assume clocks set to the lunar
    day. *Earth evidence.*

### Sky, land and people

The author asked on 9 October for the aerial biomes to come almost first, and set a working scale of 20–30 billion
affluent people across the Solar System over the near horizons. The study is [people_and_land.md](people_and_land.md),
with its screen [people.py](people.py); numbers here are from it.

21. **The sky is the Moon's largest habitat.** Air above freezing fills 9.2×10⁸ km³, 24 km deep over the ground on
    average, and each square metre carries 7.1 Earth columns of air. It holds 130–240 mm of water vapour, 5–10 times
    Earth's, turned over in 98–158 days where it rains, and cosmic-ray ionization stays at 3–16% of Earth's sea level
    up to 90 km, so the sparse life of the high air answers to cold, dry air and scarce nutrients.
22. **Lift makes room for large flyers, and the land feeds them.** Earth's largest flyers carried over reach 60–66 m
    across and 44–55 t, and an animal hovers on 6.3% of Earth's power. Phosphorus falling from the air feeds at most
    four ten-thousandths of the land's growth, so the sky's food webs eat what the land and seas send up: pollen,
    spores, nectar, insects and the sea–forest migrants that carry phosphorus inland.
23. **Rain limits crops.** Earth's water-use efficiency, raised a fifth by the design air's CO2, lets rain support the
    equivalent of 2.5–5.7 million km² at the crop model's yields, 83% of it within 30° of the equator.
24. **The vegetated Moon's own ethylene settles at 94–233 ppb.** With soils its only sink, at the uptake estimated for
    Earth's soils, it passes the 50 ppb that costs wheat 36% of its yield, whatever the population; soils ten times
    faster hold it at 9–23 ppb. The plants' emission and the soils' uptake are design variables of the whole biosphere.
25. **People's gases build up with no sink.** At Earth's agriculture and waste, each billion people adds 6–28 ppm of
    methane at steady state and 0.23–0.33 ppm of N2O every thousand years; the lunar design's sealed waste and absent
    ruminants bring that to 0.06–2.8 ppm of methane a billion.
26. **The lunar design cuts what each billion people takes from the land tenfold.** Each billion living by Earth's
    practice grazes 1.2 times the fog desert, crops 6–36% of the rain-fed land and lights a third to three-fifths of
    the dry land. Each billion living by the lunar design's plant-rich diet, compact cities and shielded lamps crops
    2.5–16% of the rain-fed land and lights 3–6% of the dry land, and 0.035–0.084% more dimming answers its heat. How
    many people the Moon's settings hold, on the ground and aloft, comes from provisioning's comparison in the shared
    population cases.
27. **People's light reaches the sky before their crops fill the land.** Lit as affluent countries light themselves,
    each billion people lights more land than any diet crops, and at the metropolis's busy hour a flyer passes a point
    in the low sky every four to six minutes. Dark corridors for the migrants, quiet low traffic at dusk and shielded light keep
    the sky's biomes working.

## Statements elsewhere that need correcting

The evaluation ([state.md](state.md#statements-elsewhere-that-need-correcting)) lists stale or conflicting
statements in other folders. The joint integration of 12 October 2026 (7c230dd) corrected the register's inputs, the
forests' leaf area and air, plant twilight and the waters' fines, and the rest wait. The largest are:
- the ecology register's numbers, which predated the corrected design case (rerun on the corrected products in the
  joint integration, 7c230dd);
- the aerosol study's oxidation chemistry, which predates the no-oxidant finding;
- the forests' leaf area of 22, against the canopy model's optimum near 6 (set to 6 in the joint integration, 7c230dd);
- the joint synthesis's flat 2.3 lux of Earthlight at every nearside place, against 0.48 lux at Procellarum by Russell
  (fixed on main on 9 October, 78020e0).

## Files

| File | Holds |
|---|---|
| [screen.py](screen.py) | The first screen: ultraviolet by action spectrum, settling at lunar gravity, methane and N2O with soils as the only sink |
| [results/first_screen.json](results/first_screen.json) | Its product (schema `terluna.biosphere.ecology-first-screen/1`) |
| [state.md](state.md) | What the project already holds for ecology, what others assume of it, and statements to correct |
| [literature.md](literature.md) | The literature by topic, with what it means for the Open Moon |
| [questions.md](questions.md) | The open questions, in order |
| [sources.json](sources.json) | Every source cited, its use and how far it was read |
| [people_and_land.md](people_and_land.md) | Sky, land and people on the Open Moon: the sky's biomes first, the land, nutrients, gases and heat, and what binds first |
| [people.py](people.py) | Its screen, writing [results/people.json](results/people.json); sources in [people_sources.json](people_sources.json) |
| [taxonomy.md](taxonomy.md) | The taxonomic system: principles, the author's choices, the functional classification, the first register and names for the aerophytes |
| [taxa.py](taxa.py) | The register as code, writing [results/taxa.json](results/taxa.json) (schema `terluna.biosphere.ecology-taxa/1`); sources in [taxonomy_sources.json](taxonomy_sources.json) |

`python -m biosphere.ecology.screen`, `python -m biosphere.ecology.people` and `python -m biosphere.ecology.taxa`
rewrite the products; `python -m pytest biosphere/tests/test_ecology_screen.py biosphere/tests/test_ecology_people.py
biosphere/tests/test_taxa.py` checks them.
