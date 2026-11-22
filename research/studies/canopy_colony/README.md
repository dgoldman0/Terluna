# A canopy colony for the seas' oxygen

22 November 2026, on `domain/ecology`. The seas' deep water runs short of oxygen, and the physical supplies have been counted and found small: the winds' mixing, the storms' gusts and the rivers. What is left is the living side. The author's brief of 22 November, "Oxygen capsules and a floating canopy colony for Terluna's seas" (kept on the research drive under `conversation-exports/author-brief-oxygen-capsules-2026-11-22/`), proposes a floating colony whose members share the work as a siphonophore's do: a lit canopy held 20–50 m down by gas floats, a trunk that pumps surface water down, feeders that catch sinking matter, and capsules of calcium peroxide under a carbonate shell that carry oxygen below. This study refines that design into a system whose parts are each sized against the sea, asks what the colony gets back for the work, carries the design to the fresh seas, and lists what it would give the people of the Moon.

It is a study of design and budgets. Nothing here designs an organism: each part is described by the job it must do, the amounts it must carry and the physical limits on carrying them, with Earth's organisms as reference points for rates that life already reaches. The sea is a one-dimensional column ([column.py](column.py)) on the South Pole–Aitken Sea's own floor, mixed as three model years of that sea's circulation run under the Moon's own mixing have it ([the south sea study](../south_sea_ocean/README.md), whose reading of that run and whose correction for the Moon's sinking were recorded on `domain/hydrosphere` the same day). [run.py](run.py) writes [results/canopy_colony.json](results/canopy_colony.json) (schema `terluna.research.canopy-colony/1`); [seas.py](seas.py) measures the seven seas of the design surface ([results/seas.json](results/seas.json)) and [currents.py](currents.py) reads the south sea runs' flow by depth ([results/currents.json](results/currents.json)); [sources.json](sources.json) lists the outside sources and how far each was checked.

## What it finds

**Left alone, the salt sea loses its oxygen in a band that opens at 300–500 m and spreads downward** *(computed on the column)*. Under the Moon's own mixing the air renews the water down to about 250–330 m, and below about 400 m the water lives on the oxygen it holds. Sinking matter decays higher on the Moon than Earth's curve places it, because a particle sinks 3.6 times slower, and about 15 times shallower if the warm water also doubles the decay for each 10 K. Starting from a saturated sea (6.36 mg/L), the layer at 270–600 m, the water cool enough for the eels' eggs, holds:

| Oxygen at 270–600 m after 100, 300 and 1,000 years, mg/L | 12 g C/m² a year | 34 g C/m² a year | 64 g C/m² a year |
|---|---|---|---|
| Earth's curve, as the south sea study first had it | 4.6, 2.1, 0.6 | 1.3, none, none | none |
| The Moon's slower sinking | 5.4, 4.3, 1.8 | 3.8, 0.8, 0.2 | 1.8, none, none |
| Slower sinking and the warm water's faster decay | 6.0, 5.6, 4.7 | 5.4, 4.2, 1.8 | 4.5, 2.3, 0.3 |

The exports are the carbon sinking below 100 m: 34 g is the sea oxygen study's nutrient-poor state, 12 g Earth's poorest gyres and 64 g the middling state. In the central case, the Moon's sinking at 34 g, the lowest level still holds 2.3 mg/L after a century, at 330 m; by year 300 the water from 240 to 640 m is under 2 mg/L and that at 400–480 m holds none, and by year 1,000 everything from 240 m to 1,100 m is under 2 mg/L. The upper 220 m keep 3.5–5.5 mg/L throughout. With the fitted mixing above 300 m halved the band opens at 200 m within a century, and with it doubled the band still opens, a century or two later, so the finding holds across that range while its timing does not.

**The demand that matters is small and sits just under the water the air renews** *(computed)*. Of the demand below 100 m, taken over the sea's own floor:

| Share of the demand below 100 m | 100–250 m | 250–400 m | 400–850 m | Below 850 m | Of it on the floor |
|---|---|---|---|---|---|
| Earth's curve | 55% | 16% | 16% | 14% | 10% |
| The Moon's slower sinking | 80% | 8% | 7% | 5% | 4% |
| Slower sinking and faster decay | 93% | 3% | 2% | 1.5% | 1% |

Four fifths of the matter decays where the sea ventilates itself. One fifth, 0.6 mol of oxygen demand per m² a year at 34 g, passes 250 m into water that mixing reaches slowly or not at all, and that fifth empties the still layer over centuries because nothing replaces what it takes. So the breakdown of matter is handled where the matter lies by making sure of three things: that as little as possible passes 250 m, that the water from 250 m to the bottom of the still layer is flushed, and that the little reaching the water below has its oxygen sent after it.

**Three parts hold the sea, each at its own depth** *(computed: the steady state of the column)*. The sizes below are the least that keep every level below 100 m at 2 mg/L or more for good, and in the second trunk column also keep 4 mg/L in the 270–600 m layer, with the feeders taking half of the matter that reaches 250 m, the trunk letting its water out between 250 and 900 m, and the capsules releasing below 900 m:

| Per m² of sea a year | Feeders catch (as oxygen demand) | Trunk: water sent down, for 2 mg/L everywhere | For 4 mg/L at 270–600 m too | Capsules release below 900 m |
|---|---|---|---|---|
| The Moon's sinking, 34 g C (central) | 0.29 mol | 1.0 m | 1.8 m | 0.063 mol |
| 12 g C | 0.10 mol | 0.26 m | 0.50 m | 0.022 mol |
| 64 g C | 0.54 mol | 3.1 m | 3.6 m | 0.118 mol |
| Earth's curve, 34 g C | 0.65 mol | 2.8 m | 4.4 m | 0.177 mol |
| Slower sinking and faster decay, 34 g C | 0.09 mol | 0.24 m | 0.50 m | 0.019 mol |

In the central case that is a catch of a tenth of the export (3.4 g of carbon per m² a year), 5 litres of surface water per m² a day down the trunks, and 0.06 mol of capsule oxygen. The steady sea then holds 3.3 mg/L at its lowest, at 330 m, 4.0 mg/L across 270–600 m and 5.1 mg/L from 850 m to the floor. A trunk that stops at 500 m needs less water (1.5 m a year) and leaves twice the work to the capsules (0.12 mol); with no feeders the trunk must send 3.7 m and the capsules 0.13 mol; with three quarters caught, 0.9 m and 0.03 mol. The sizes follow the export nearly in proportion, the capsules exactly and the trunk within a quarter, which matters for the loop below.

**Each part is needed, on its own time scale** *(computed: the central system run for a thousand years from a saturated sea)*:

| After 1,000 years | 270–600 m | Lowest level | Under 2 mg/L |
|---|---|---|---|
| The whole system | 4.1 mg/L | 3.4 mg/L at 330 m | nowhere |
| Without the feeders | 2.0 | 1.3 at 400 m | 240–440 m |
| Without the trunk | 0.4 | none at 400 m | 360–770 m |
| Without the capsules | 4.0 | 3.4 at 330 m | nowhere yet |
| Nothing | 0.2 | none | 240–1,100 m |

The trunk decides the first centuries: without it the feeders alone only delay the band. The feeders halve what the trunk and capsules must carry. The capsules are for the long run. The water below 900 m holds enough oxygen for two to three thousand years of the little that reaches it, and then fails; the capsules' whole job is to keep that stock from running down.

### The feeders

The feeders hang along the trunks at 150–250 m, in water the air renews, so their own breathing costs the still layer nothing. What they catch they digest there and send up the trunk as dissolved food, and the colony respires it in the lit layer, where the air and the light pay for the oxygen. A litre of body fluid carrying 10 g of carbon moves as much demand as 4,000 litres of saturated water moves oxygen, so a third of a litre per m² a year carries the whole catch.

How much feeding surface half the catch needs follows from how long a particle stays in the tier *(computed from the south sea runs' flow; the sinking speeds are recalled)*. A particle sinking at speed w spends H/w in a tier H thick while the water carries it past the colony at the relative speed u, so feeders of efficiency e whose surface facing the flow is F m² for each m² of sea take the share 1 − exp(−eFu/w). The yearly-mean flow of the oxygen run moves at a median 5.3 cm/s at the canopy's depth and 1.2 cm/s at 220 m, so the water at the feeders passes a colony drifting with its canopy at a median 4.3 cm/s (1.0 cm/s at the slowest tenth of the sea). Matter that sinks at 20–100 m a day on Earth sinks at 6–28 m a day here. Half the catch then needs eF of 0.1–0.5% of the sea's surface in the median flow and up to 2% in the slow tenth: with feeders that take a tenth to a half of the particles in the water passing them, a facing surface of 0.2–5% of the sea, about 1% in the middle of the ranges. As a picture, that is a hundred trunks under each km², each fringed a metre wide over its 100 m of tier. The Moon's slow sinking helps twice: it keeps the matter high, and it gives the same net 3.6 times longer to meet each particle.

The brief models the feeders on deep-water corals and notes that the known reef builders are cold-water animals. A closer reference exists: corals without algae live at 200–1,000 m in the central Red Sea in water above 20 °C that is poor in both food and oxygen, and they do it with low respiration, slow calcification and little tissue (Roder and colleagues 2013). That is the temper the feeders need, since the catch is a thin income. At 3.4 g of carbon per m² a year, and a tissue that respires two to seven times its own carbon a year at 28 °C *(recalled rates, carried up from cold water)*, the catch keeps 0.5–1.7 g of living carbon under each m². The trunks and nets must therefore be mostly structure that does not respire, a cuticle or gel with a thin living layer, as a larvacean's house or a siphonophore's stem is. Colonies of that build reach the needed lengths: a siphonophore's stem has been filmed at about 46 m and Praya dubia is reported at 40–50 m.

### The trunk

The trunk's water carries the air's own oxygen, so it runs through the two-week night. Its difficulty is heat. Surface water is 2.5 K warmer than the water at 250–900 m, and the 270–600 m layer lies at the eels' limit of 28 °C, so water piped down uncooled would warm the water it enters by 0.7 K a century. The trunk must give its water's warmth to the water rising around it on the way down, as a counter-current exchanger does. A thin-walled channel whose water arrives a small lag warmer than its surroundings can carry Nu π κ times that lag over the temperature gradient, whatever its bore *(laminar flow in a tube; computed)*:

| Lag at the foot | One channel carries | Channels per km² for 1.8 m a year | Heat leaked downward | Warms the foot, if nothing carried it off |
|---|---|---|---|---|
| 0.3 K | 12 m³ a day | 410 | 0.07 W/m² | 0.08 K a century |
| 0.1 K | 4 m³ a day | 1,240 | 0.02 W/m² | 0.03 K a century |
| 0.03 K | 1.2 m³ a day | 4,140 | 0.007 W/m² | 0.008 K a century |

The floor's own heat is 0.011 W/m², so the slowest channels leak less than the Moon does. A hundred trunks per km² then carry four to forty channels each, with the water moving at 1–7 cm a second in a bore of 5 cm.

Pushing the water down costs almost nothing, and in part of the sea nothing at all *(computed from the oxygen run's yearly mean)*. Brought to the foot's temperature, the surface water is lighter than the foot's water by 0.2 kg/m³ at the median of the sea, since the rivers freshen the surface, and by 0.9 kg/m³ at the freshest twentieth. Over a trunk reaching to 575 m on average that is a head of 200–860 Pa, a third of a pascal to a pascal and a half for each metre of trunk, which cilia along the wall can hold; the work is 0.4–1.6 kJ per m² a year, 0.3–1.1% of the energy in the catch. Over 17% of the sea, where evaporation has made the surface the saltier, the cooled surface water is the denser and sinks through the trunk by itself: Stommel's salt fountain, run downward (Stommel, Arons and Blanchard 1956 describe the upward one).

Of the 0.37 mol of oxygen that the central trunk carries down per m² a year, the water gives up 0.13 mol at the foot and rises with the rest. The rising water is the other half of the trunk's work: it carries back up what decayed in its reach.

### The capsules

The capsules serve only the water below the trunk, so their job is the smallest in mass and the longest in time: 0.06 mol of oxygen per m² a year in the central case. With the brief's placeholder that 70% of capsules release where they are needed, the canopy makes capsules carrying 0.09 mol, which takes 0.18 mol of calcium and yields 13 g of calcium peroxide per m² a year.

The brief limits the capsules by alkalinity: each mol of oxygen leaves two of lime in the deep water, which returns only as fast as the deep water does. This study reads the lime's fate differently *(chemical reasoning, not computed and not tested)*. The oxygen a capsule releases is respired on the spot and yields carbon dioxide, half of what the capsule's own lime can take up; the rest is taken from the deep water's dissolved carbon, of which the sea holds some 4,000 mol under each m² *(seawater's usual 2 mmol/kg, recalled)* against a need of a tenth of a mol a year; and the product, calcium carbonate, stays solid, since at the Moon's pressures (34 atmospheres at 2.1 km) and at 27 °C the deep water stays saturated with it. The net of the whole round, canopy and deep together, is calcium and carbon dioxide leaving the water as carbonate mud, the same as a calcifier's shell sinking to the floor. So the capsules are a form of the calcifiers the seas already count on for their calcium ([the seas' make-up](../sea_chemistry/README.md): 150 Mt of calcite a year; the register's row on reef builders and shellfish), with oxygen riding on the way down. If that reading holds, what limits them is calcium supply and the canopy's rate of making them:

- **Calcium.** Shared out by river flow, the held seas' calcite gives this sea's rivers 0.13 mol of calcium per m² a year, enough for 0.046 mol of delivered oxygen. The central 0.063 mol takes 1.4 times that. The difference comes from the sea's dissolved calcium, 25,000 mol under each m², at 0.2% of the stock a millennium. In a salt sea the capsules are therefore not held to the rivers' rate for any time that matters here, and they do compete with the reef builders for the rivers' share.
- **The canopy.** At reefs' rates of laying down carbonate, 10–100 mol per m² of reef a year *(the brief's range, recalled)*, 0.18 mol needs canopy over 0.2–1.8% of the sea. Peroxide mineral has never been made by an organism, as the brief says, so this is a best case.
- **Electrons.** Each mol of delivered oxygen takes four electrons of the light's, the cost of fixing one carbon: 1.1 g of carbon forgone per m² of sea a year. A lit canopy in a nutrient-poor sea has light to spare, and symbiotic algae already pass surplus electrons to oxygen and make hydrogen peroxide when light outruns their carbon fixing (the Mehler reaction; corals release peroxide to the water).

What opens a capsule where it is needed is unproven. The seafloor ideation of 16 November has microbes open it; a shell of carbonate also thins fastest in the acid film around decaying matter, which would put the oxygen where the decay is.

## What the colony gets back, and why the loop holds

The work is real, and a colony that only served the sea would not last. Per m² of sea a year in the central case, the colony receives *(computed, with the recalled rates marked)*:

- **Food.** The catch is 3.4 g of carbon, 137 kJ, an income fed by the same sinking matter that makes the demand.
- **Nutrients.** The catch holds 43 mmol of nitrogen and 2.7 mmol of phosphorus. Open water held nutrient-poor ([the register's trophic design](../../decisions.md#the-design-climate-and-the-waters-19-october-2026)) is water in which a lit canopy starves for exactly these. The feeders pay the canopy in nitrogen and phosphorus and the canopy pays the feeders in sugar and oxygen: a coral's bargain between animal and alga, stretched over 250 m of trunk.
- **Light's carbon.** A canopy over 0.2–1.8% of the sea, handing on 300–1,000 g of carbon per m² of canopy a year *(reefs' production, recalled; the two-week night is not counted)*, brings 0.5–18 g of carbon per m² of sea, of the catch's order or above it.
- **A use for surplus light.** Making peroxide gives the canopy somewhere to put electrons it cannot use, and the acid that forming calcium peroxide releases at the surface turns bicarbonate into the carbon dioxide the algae take up, as calcification does for a coral's.
- **Breathable depth.** The trunk's water passes the colony's own deep members before the sea gets the rest.

Against that income, pumping costs 0.3–1.1% of the catch's energy and the capsules' electrons a carbon's worth of light for each mol delivered. The colony's upkeep is the large term, and it is the same thing as the service: the caught matter respired in the lit layer is the demand removed from the deep.

Four relations make the loop hold itself *(reasoning on the computed sizes; no population is modelled)*:

1. **The income follows the job.** More export means more sinking matter at the tier, a larger catch, and a colony that can keep more feeders and channels; the sizes the sea needs rise with the export in nearly the same proportion (table above). A colony whose mass follows its food sizes itself to the sea's demand, up and down.
2. **The colony cannot enrich the sea beyond its rivers.** It returns nutrients that would be lost, and it adds none. With the feeders and the trunk, 2.2% of the export's nutrients still sink below the trunk's reach each year, 18 mg of phosphorus per m²; rivers bringing 0.68 m of water a year over this sea replace that at 27 mg of available phosphorus per m³. Each mg/m³ in the rivers holds 1.3 g of carbon export. Without the colony 21% sinks below 250 m, 174 mg of phosphorus per m² a year, which no river replaces: at 10 mg/m³ the top 250 m hold 2.5 g of phosphorus per m², fourteen years of that loss, so the lit layer would strip itself of nutrients into the still layer within decades, and the sea would be both poorer at the top and emptier of oxygen below.
3. **Oxygen and nutrients travel the same loop in opposite directions.** The trunk's descending water carries oxygen and its rising water carries the regenerated nutrients of its reach, 8.6% of the export's, back to the light. Pumping alone is nearly oxygen-neutral over time, as the brief notes, because the nutrients it lifts make new matter that sinks again; the feeders break that circle by taking half the matter before it reaches the still layer, so the two together gain.
4. **A colony that dies costs little.** Its living carbon is 0.5–1.7 g per m², a twentieth of a year's export, so a die-off after a storm or a long night is a small pulse of demand.

The colony's catch is also a lever on how rich the sea is. Catching three quarters halves the loss below the trunk's reach, so the same rivers hold twice the export; the author's open question of how the salt sea's waters are kept rich enough for its shores then has one answer in the colony's density. The register's row of 22 November takes this way for all the seas: the open water is made richer by recycling what it receives, and feeding it from the deltas is set aside.

## The fresh seas and the lakes

The same three parts apply in fresh water, with four differences that change their sizes *(computed where a number is given)*:

- **Fresh water holds more oxygen**, 7.66 mg/L at the fresh seas' median surface against 6.36 in the salt sea, so each litre the trunk sends carries a fifth more.
- **Calcium is scarce.** The held fresh sea carries 0.84 mmol/kg, 300–1,400 mol under each m² against the salt sea's 25,000, so capsules there live on what the rivers bring and take it from the reef builders and shellfish that the register counts on for the same calcium.
- **No salt helps the trunk.** Nothing sinks by itself, so every trunk is pumped, against a smaller head since no fresh layer rides on top.
- **Anoxic fresh water makes methane and frees phosphorus**, where salt water makes hydrogen sulfide. Keeping the mud's surface oxic keeps its phosphorus locked, which is what holds the open water nutrient-poor, and stops the methane that the air's budget counts.

How large the job is in a fresh sea rests on a question no study has answered: how fast the sea's cold end renews its deep water. A fresh sea has no salt to hold its deep water down, so the coldest water its surface makes, at its high-latitude end on its coldest nights, sinks to the floor and lifts the rest. [The sea oxygen study](../sea_oxygen/README.md)'s columns bracket that renewal between years and centuries by the form of the mixing, and call for the seas' own circulation to decide. This study gives the two ends:

- **If the cold end renews the deep within years to decades**, the deep water needs no help while the sea stays nutrient-poor, and the colony's only work is the feeders', which lower the demand every level sees.
- **If a fresh sea is as still as the salt sea**, its need is the table below: the same column on each sea's own floor, with the salt sea's mixing borrowed and fresh water's saturation, at 34 g of carbon and the Moon's sinking.

| Fresh sea, as still as the salt sea | Mean depth | Floor shallower than 250 m | Unaided, 270–600 m after 1,000 years | Trunk for 2 mg/L everywhere | Capsules below 900 m | Rivers' calcium | Capsules' share of it |
|---|---|---|---|---|---|---|---|
| Near-side Sea, main part (Uamoana), 23° S–71° N | 725 m | 17% | 0.4 mg/L | 0.73 m a year | 0.022 mol/m² a year | 0.12 mol/m² a year | 51% |
| Near-side Sea, southern part (Bultengis), 32° S–9° N | 379 m | 39% | 0.2 | 0.51 | 0.003 | 0.04 | 21% |
| Humboldtianum (Yuldustengis), 48–66° N | 1,193 m | 9% | 0.4 | 0.72 | 0.043 | 0.14 | 88% |
| Moscoviense (Hanakosa), 12° S–37° N | 1,706 m | 8% | 0.4 | 0.70 | 0.048 | 0.38 | 36% |
| Orientale (Bahari Masariki), 32–6° S | 1,574 m | 13% | 0.4 | 0.68 | 0.042 | 0.38 | 32% |
| Smythii–Marginis (Dagat na Kaseh), 9° S–27° N | 1,411 m | 11% | 0.4 | 0.69 | 0.044 | 0.32 | 40% |

So at the still end a fresh sea asks half to two thirds of the salt sea's trunk and at most three quarters of its capsules, and the capsules fit inside the rivers' calcium everywhere, tightly in Humboldtianum. With the faster decay the trunk falls to 0.10–0.16 m a year and the capsules to a third. Which seas lie near which end can be guessed from their shapes: Humboldtianum lies wholly poleward of 45° and a quarter of the main Near-side Sea does, so both have cold ends; Moscoviense, Orientale, Smythii–Marginis and the southern Near-side Sea lie within 37° of the equator and have none to speak of, and the sea oxygen study found the equatorial basins' deep water the poorest. The rivers do not help: the rivers' oxygen study (`research/studies/river_oxygen` on `domain/hydrosphere`) finds they float. The fresh seas' need is an upper bound until one of them has a circulation run.

The Moon's slow sinking bears on the lakes and on the sea oxygen study's own demand as well. The lake model puts most of the deep demand on the floor, 78% of it in the south sea's basin, where the Moon's sinking puts 4%. In a crater lake 352 m deep, 14% of the matter passing 100 m would reach the floor by the same curve. The lake studies' mechanical remedies ([the deep lakes' anoxia](../lake_anoxia/README.md)) are sized on the floor's whole demand, 0.3–1.3 t of oxygen per km² a day; a colony in a deep lake would be a short one, feeders and a trunk to the floor, and the sizes there wait on the lake model taking the same correction. That correction is not made here.

## What it gives the people

- **The eels' water.** The layer cool enough for the eggs, from about 400 m down, kept at 4 mg/L across the whole sea, where the south sea study sized 20–50 MW of injected oxygen for one spawning ground 10 km across *(computed there and here)*.
- **A sea without a poisoned deep.** An anoxic salt deep gathers hydrogen sulfide and an anoxic fresh deep methane, and either is a hazard to shores and to the air if stirred upward. The oxygenated sea holds neither *(reasoning; the brief's Black Sea comparison is recalled)*.
- **Clear water kept clear, and fertility by choice.** Oxic mud keeps its phosphorus, which the nutrient-poor design rests on, and the colony's catch sets how much of the rivers' nutrients the lit layer keeps *(computed in the column's terms)*.
- **Carbonate.** The spent capsules lay down 18 g of calcium carbonate per m² a year, 43 Mt a year over the salt sea, as mud and sand; the register already counts the calcifiers' carbonate as material.
- **An oxidant to harvest.** The canopy's capsules are 13 g of calcium peroxide per m² a year, 31 Mt over the salt sea. On Earth, rice seed coated with calcium peroxide is sown straight into flooded soil, the coating releasing the oxygen the seedling needs (Westcott and Mikkelsen 1983 and later work); a share of the capsules would serve the Moon's wet fields the same way, and bleaching and disinfection besides *(the uses beyond rice are recalled)*.
- **Floating reefs.** Canopy over 0.2–1.8% of the salt sea is 4,000–43,000 km² of lit living structure in open water, 20–50 m down and so under any keel: habitat where the open sea has none *(a proposal)*.
- **A sounding line in every colony.** Lines that live wild carry reporting traits by the register's rule; a colony that spans the top 900 m could report the oxygen and warmth along its length *(a proposal)*.

## The brief's open questions

| The brief asks | Answer |
|---|---|
| The mixed layer's depth, and how much decays inside it | The air renews the water to about 250–330 m under the Moon's own mixing; 80% of the demand below 100 m falls between 100 and 250 m with the Moon's sinking (55% by Earth's curve). |
| Whether the demand slows sinking for a sixth of the gravity and speeds decay for warm water | It did not. The south sea study now carries both, and they move the case from "under 2 mg/L from 104 m" to the band at 300–500 m found here. |
| Which demand form the export scaling uses | The ocean form with Earth's curve. |
| The flux crossing 850 m | 16% of the export by Earth's curve, 6% with the Moon's sinking, 2% with faster decay: 0.45, 0.16 and 0.05 mol of oxygen per m² a year at 34 g, less over the sea's own floor. |
| The phosphorus stock planned, and the rivers' input | Not set by any study. The sea chemistry study gives 0.9–116 mg/m³ in a salted sea while the lakes trap the land's phosphorus. This study gives the relation instead: with the colony, each mg/m³ of available phosphorus in the rivers holds 1.3 g of carbon export. |
| The rivers' calcium and alkalinity | 0.13 mol of calcium per m² of this sea a year, the held seas' 150 Mt of calcite shared out by river flow; alkalinity arrives in excess of calcium (the seas' make-up). |
| How likely the deep is to be stirred upward, and the paused run's result | The run's three years show the still layer and leave the oxygen between wide bounds (the south sea study). The deep pool below 1.2 km is stirred within itself and cut off from above; the densest surface water stays 0.13–0.21 kg/m³ lighter than it. |
| Each fresh sea's depth, hardness and mixing | Depths and floors in the table above; hardness as the held fresh sea, 0.84 mmol/kg of calcium; mixing unknown, bracketed between the two ends. |

What the refinement changes in the brief: the trunk's flow falls from 30 to 5 litres per m² a day, because most of the 100–500 m demand is met by the sea itself; the trunk becomes a heat exchanger; the capsules' limit moves from alkalinity to calcium and the canopy's rate; the capsules' range moves from the whole deep to the water below the trunk; the feeders' model moves from cold-water corals to the Red Sea's warm deep corals; and the sizes are steady-state requirements where the brief's were first-order.

## How far Earth's organisms go toward each part

No organism does what the colony is asked to do, and nothing in this study shows that one can. What can be shown is how far living things already go toward each part, with a measured figure beside what the design asks *(search summaries of the sources named, 22 November; none read in full; the design's figures are this study's)*:

| Part | What the design asks | The nearest thing measured | The gap |
|---|---|---|---|
| Catching half of the sinking matter | The tier's water cleared of particles at 4–19% a day, since a particle takes 4–18 days to sink through 100 m here | In the North Atlantic's twilight zone the animals fragment and ingest about half of the fast-sinking particles (Giering and colleagues 2014). Giant larvaceans at their densest in Monterey Bay would clear every particle from their depths in 13 days, 8% a day, each filtering up to about 80 litres an hour (Katija and colleagues 2017). Flux feeders account for a substantial part of the flux lost in the shallow twilight zone (Stukel and colleagues 2019) | None in rate: midwater animals together already take about half. The gap is that one kind of colony does it, and holds its catch instead of repackaging it into faster-sinking pellets |
| Carrying the catch up 250 m | A third of a litre of fluid per m² of sea a year, at any speed | Giant kelp moves its sugars through sieve tubes at 0.35–0.72 m an hour (Schmitz and Srivastava 1979), and kelps as a group at 0.05–0.78 m an hour, over fronds tens of metres long | 250 m takes two to four weeks at kelp's speed, which the budget allows. The distance is five to ten times a kelp's |
| Pumping surface water down | 1.2–12 m³ a day in each channel, against a third of a pascal to a pascal and a half per metre of trunk | A sponge pumps 10–30 times its own volume a minute; one of its pumping chambers raises about 5.5 Pa and a whole sponge up to about 24 Pa; pumping costs 1–22% of its respiration (Riisgård, Larsen and colleagues) | Small. A litre of such tissue moves a channel's flow, and one chamber's pressure covers 4–16 m of trunk, so pumps every few metres suffice. The upkeep of that tissue has to fit the 0.5–1.7 g of living carbon the catch keeps |
| Giving up heat on the way down | Thin-walled channels in the rising water | Counter-current exchangers in the limbs and gills of warm-bodied fish and whales *(recalled)* | Not looked up |
| Living in warm, poor deep water | 27–28 °C, little food, 3–5 mg/L of oxygen | Corals without algae at 200–1,000 m in the Red Sea at 21.5 °C, on low respiration and little tissue (Roder and colleagues 2013) | Six kelvin warmer |
| A colony of this length | A trunk reaching 575 m on average and 900 m at its foot | A siphonophore's stem filmed at about 46 m; kelp of similar length | Ten to twenty times longer than any colonial animal known: the largest gap in form |
| The capsule | Calcium peroxide grown under a carbonate shell | Hydrogen peroxide made in the light by corals' algae and released by corals (Suggett and colleagues 2008; Armoza-Zvuloni and Shaked 2014). No organism is known to lay down calcium peroxide | The chemistry itself, as the brief says |

So the rates are within what animals reach; the form is not. The catch and the pumping are things midwater animals and sponges already do at the needed pace. What has no precedent is their union in one floating body half a kilometre long that keeps what it catches, and the peroxide mineral. Of the three parts, the order of confidence is therefore the trunk's pumping, then the feeders, then the capsules, and the sea's budget is kindest to that order: the capsules matter only after two to three thousand years, and on the supply the seas actually receive ([the seas' phosphorus budget](../sea_nutrients/README.md)) the colony's first use is the feeders' recycling.

## How it is set up

- **The column.** Thirty-two levels, the circulation run's, from 5 m thick at the surface to 500 m at depth. Each level's volume and the open area at its edges follow the sea's floor (the south sea study's stored basin), so matter sinking onto a slope is spent at the slope's depth. The surface level is held saturated. Levels mix at the sea-wide mean of the closure's diffusivity in model year 53 of the oxygen run, the weaker of two neighbours limiting their exchange.
- **The fit.** Above 300 m the mean diffusivity is multiplied by the one share, 0.125, with which the column follows that run's oxygen spent at 93–269 m through its three years (root-mean-square error 0.15 in the logarithm). A sea-wide mean of a mixing that comes in bursts and patches overstates what the mean water feels. Below 300 m nothing is fitted, and the column follows the run there within a quarter.
- **The demand.** The export leaves 100 m and falls as Martin's curve with the depth below 100 m stretched by 3.6, or 14.9 with the decay's factor of 2 for each 10 K ([run.py](run.py) repeats the south sea study's `sinking_stretch` until the branches join). A mol of oxygen is counted for each mol of carbon, as that study counts.
- **The parts.** The feeders remove a share of the flux at 250 m. The trunk lets surface water out evenly by thickness over its range, and each level below the surface passes its water to the one above. The capsules release below the trunk's foot in proportion to the demand left there. The trunk's size is the least flow at which the column's steady state meets the stated floor; `None` would mark a case that 60 m a year cannot hold, and none occurs.
- **The seas.** [seas.py](seas.py) measures each sea of the design surface at its own level on the redrawn floor at 16 pixels per degree, with the river water the design routing sends into it. The fresh seas' columns take the salt sea's levels and mixing on their own floors.
- **The currents.** [currents.py](currents.py) reads the yearly-mean flow of the 50-year run's last year and the oxygen run's year 53, without changing either run, and compares the top 35 m's water, cooled to the foot's temperature, with the foot's.
- **The input copy.** `inputs/south_sea_lunar_mixing_1deg.json` is a byte copy of the south sea study's product on `domain/hydrosphere` (commit `18916584`); the product records its hash, and a test compares it with the source once both are in one checkout.

## Assumptions and limits

- **One column, one run, three years.** The mixing is one circulation model's on the atlas's basin under the earlier fresh-sea climate, read over three model years; the design surface's sea is 300 m deeper on average. Halving or doubling the fitted share moves the unaided band's opening between one and three centuries.
- **A column has no sideways motion.** The colonies drift about 1,600 km a year with the water at their canopies; where the flow gathers or strands them, how sinking matter is focused, and how the cold ends and river arms differ from the mean are outside it. The trunk's upward return is taken with an upwind rule on levels 60–150 m thick, which smooths the foot's range by itself.
- **The sinking is Stokes's.** Aggregates large enough for their drag to grow with speed gain less from gravity, so 3.6 is the upper end for them; particle sizes on the Moon are untested. The decay's factor is a usual value, not a measurement; Earth's warm waters are reported to spend their sinking matter shallower than its cold ones (Marsay and colleagues 2015), which is the direction assumed.
- **The catch is a target.** No organism is shown to take half of a sinking flux over a sea; Earth's midwater animals do intercept a large part of theirs *(recalled)*. The feeding surface rests on yearly-mean currents, which are low values for the moving water, and on a range of efficiencies that is a guess.
- **The capsule is undemonstrated**, as in the brief: no organism makes calcium peroxide, the shell's life in warm water is unknown, and the 70% delivery is the brief's placeholder. The lime's fate as carbonate is this study's reasoning; if it stayed dissolved, the brief's alkalinity ceiling of 0.1–0.2 mol would still cover the central 0.06.
- **The export is given, not computed.** In a closed sea it depends on the nutrients. [The seas' phosphorus budget](../sea_nutrients/README.md), made the same day, finds that rivers and dust sustain 0.01–1.8 g of carbon per m² a year without a colony, so the 34 g of the central case stands for a sea fertilised to the lakes' nutrient-poor class by people or by the water it was filled with. On rivers and dust the seas keep their deep oxygen unaided for more than three thousand years, the colony's sizes fall with the export, and its first use becomes recycling the little phosphorus the sea receives.
- **The colony's biology is not modelled**: its growth, reproduction, grazers, disease and behaviour through the two-week night are open, and the rates marked as recalled were not looked up.
- **The fresh seas borrow the salt sea's mixing**, and the lakes are not computed.
- Nothing here is validated against an observation.

## The checks

[test_run.py](test_run.py): the Moon's sinking is gravity's ratio times the viscosities'; a column follows its sea's floor; the demand is what the flux loses, and the feeders' share is taken exactly at their tier with the whole kept; the steady column carries each level's demand up through the water above, and the stepped column reaches that state and never spends more than the water held; the trunk replaces its foot's water at the demand over the flushing rate and lifts the water above; the fit recovers the share a made-up run was mixed with; the channels' flow follows the wall's cooling; the feeding surface takes the share asked; the capsules' calcium is counted against rivers and stock; a sea is measured by depth and latitude; and the water past a canopy is the difference of two layers. On the products: the inputs are bound by hash; the column follows the paused run in the still layer, outside its fit; slower sinking keeps more oxygen and more export less; every sized system holds its sea and needs less of trunk and capsules the more is caught; the central system needs every part; the fresh seas stand on their own floors; and the product regenerates. The seas and the currents regenerate with `TERLUNA_CANOPY_COLONY_REGENERATE=1`, which reads the research drive.

## Running it

```sh
ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.canopy_colony.seas        # about 10 seconds
nice -n 10 /media/projectspace/terluna-research/venvs/ocean/bin/python -m research.studies.canopy_colony.currents   # seconds
OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.canopy_colony.run                              # about 5 seconds
```
