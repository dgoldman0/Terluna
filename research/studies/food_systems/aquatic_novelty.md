# Aquatic foods: from the waters' habitats to the table

4 November 2026, revised 7 November 2026 · `domain/food-systems`

This part of the food study designs aquatic foods for the Open Moon: fish, crustaceans and shellfish bred from Earth stock for the Moon's waters, followed from their habitats through food webs, founder lines and breeding to the dishes they could become. It sets out three hypotheses about what the lunar setting could give the food, the trials that would test each, and the culinary and safety work that turns a harvest into a product. Values carry their evidence: *published* (Earth literature, with its source and how far it was read), *product* (a project result, with its study and revision), *computed* (arithmetic on stated inputs, with who did it), *design guess* and *idea*.

## Fish in the Moon's diet

Animal foods are a luxury in the Moon's diet. The diet comes from a blend of three mixes, power-heavy, power-made and fully lamp-lit, which grow their protein as microbes and fungi ([register, 12 October](../../decisions.md#automation-the-network-and-food-12-october-2026)). Meat from animals comes to 1–2 kg a person a year, and wild culls and fishing can supply part of it where the parasite risk is kept low ([register, 12 October](../../decisions.md#the-domain-synthesiss-decisions-12-october-2026)). Wild-caught fish and shellfish therefore enter the diet at that small scale, and whether farmed fish count inside the 1–2 kg is a question the register leaves open. The food is made on the Moon, and the beams from the ring fleet are the Moon's bulk power supply ([register, 5–6 November](../../decisions.md#the-moons-energy-supply-the-loads-placement-and-the-fault-across-the-ring-56-november-2026)), so hatcheries, closed grow-out, chilling and freezing draw on the beams and add to their heat on the Moon.

The waters set the scale. Resources' sea-phosphorus study finds that the waters' life can spare 21–56 kt of phosphorus a year in all, to seaweed, fish and sorption, 1.4–3.7 kt of it from the common-level waters, so sea farms could supply the animal part of the diet for only 0.04–0.26 billion people (*product*; the production accounts' [finding 8](README.md#findings)). Fed fish release whatever phosphorus in their feed they do not retain: fed eel at 1 kg cooked a person a year for 8 billion people would release 167 kt a year, three to eight times what the waters' life can spare (*computed* in the review of 6 November on the eel calculator's reference values of 5 November; [register, 6 November](../../decisions.md#the-food-systems-branchs-evaluation-wild-lines-from-earth-stock-the-eels-seawater-spawning-and-organisms-crossing-between-worlds-6-november-2026)). The study therefore grows fed fish in closed systems whose effluent carries its phosphorus back to the land. Those systems carry the costs the production accounts find for aquaponic fish, four to six times the cropland of the animal food they replace and 0.8–1.9% of their feed's nitrogen emitted as N2O, which has no sink in the design air (findings 7 and 9). Open-water harvests take the fish the waters' own food webs grow, within the same phosphorus. The eel chapter (`eel_habitat_hatchery.md`, with `eel_budget.py`) works out per-person cases for 4–8 billion people against the heritage meat's 1–2 kg.

## Terms

A product is *qualified* for raw or for cooked service when its defined line, source, harvest phase, cut and process carry the hazard and quality evidence that [culinary evaluation](culinary_evaluation.md#three-qualification-states) sets out. Qualification belongs to a defined product and can be suspended.

*Gaean* is the world-crafting project's term for worlds with nitrogen–oxygen air and water as the principal condensable, at broadly Earthlike surface pressures and temperatures with liquid surface water and an active water cycle, whose gravity, size, rotation and starlight may differ from Earth's ([world-crafting README](https://github.com/dgoldman0/world-crafting#readme)). The Open Moon is a constructed Gaean world.

Working names such as "lunar-cycle fatfish", "alkaline sweetfish" and "butterperch" name the food a hypothesis aims at. Every line keeps its species name, and a fish other than eel prepared in the kabayaki manner is sold under its own.

## Read the study

| Part | What it holds |
|---|---|
| [Habitats and the environmental basis](aquatic_biomes.md) | The waters' temperature, oxygen, light and chemistry from the project's products, on the turned Moon's design climate; the habitats each hypothesis needs; the accounts that bound every harvest |
| [Lunar-cycle fisheries](lunar_cycle_fisheries.md) | H1: whether the 29.5-day light cycle can store energy in edible muscle and give harvests a rich and a lean phase |
| [Alkaline foods](alkaline_foods.md) | H2: which water chemistries change taste compounds, with the Earth evidence for and against |
| [Carbonate shores, shellfish and the butterperch](reef_shellfish_cuisine.md) | H3: carbonate habitats and their food webs, jade perch, shellfish and pearls, and kabayaki-style fish |
| [Founders, breeding and trials](breeding_and_trials.md) | The gates a founder passes, the trials that separate the three mechanisms under the lunar cycle, and selection |
| [Culinary evaluation and qualification](culinary_evaluation.md) | Product targets, fair sensory comparisons, raw and cooked qualification, and raw preparation across the waters |
| [`eel_habitat_hatchery.md`](eel_habitat_hatchery.md) | The eel chapter, with [`eel_budget.py`](eel_budget.py) and its product [`results/eel_budget.json`](results/eel_budget.json) |
| [`sturgeon_habitat_foods.md`](sturgeon_habitat_foods.md) | The sturgeon chapter |
| [aquatic_basis.json](aquatic_basis.json), [aquatic_checks.json](aquatic_checks.json) | The identity of every source the study read, and the checks run on it |

## Three hypotheses and what would count as progress

| Hypothesis | Mechanism to test | Food property sought | Evidence available | Result that weakens it |
|---|---|---|---|---|
| **H1: lunar-cycle fatfish** | The 708.7-hour light cycle changes the reserves of primary producers and the prey available, and some consumers store that energy in edible tissue | Fine muscle fat, a favourable melting profile, rich and lean harvest phases that recur | Earth plankton cycle their lipid daily; fish differ in where they store fat and in how it follows their diet | Energy accumulates outside edible tissue; food-web lags smooth the cycle away; the cyclic food web gives no useful muscle response |
| **H2: alkaline sweetfish** | Water ions and acid–base conditions change compatible solutes, nitrogen handling and muscle metabolism | Taste-active amino acids and nucleotides that people can taste, balanced sweetness and umami, firm flesh | Saline and alkaline rearing changes measured flavour chemistry in bass and prawns, with conflicting results between studies | Only stress, liver change or slower growth; calculated umami rises with no blinded sensory benefit; salinity alone explains the effect, which removes its carbonate-specific reading |
| **H3: carbonate-community pale fish** | Structured carbonate habitat supports particular prey and feeding; line and diet set flesh composition | Pale, sliceable or flaky flesh, localized fat and a clean flavour, with shellfish products beside | Jade perch's muscle lipid follows dietary lipid; living carbonate structures grow in warm fresh water on Earth | The habitat gives no sustainable distinctive diet or culinary benefit; prey demand exceeds production; off-flavour, biodeposition or oxygen cost outweighs the benefit |

Matched rations and matched diets separate the mechanisms. If a difference disappears when the diet is matched, the effect runs through food supply: that supports the food-web route and weakens a direct effect of photoperiod or habitat structure. Algal starch in place of lipid weakens the route through lipid-rich prey, while consumers can still make fat from other dietary energy. An H2 effect from salinity alone remains useful food production and loses its carbonate-specific reading.

H1 carries the most distinctively lunar forcing and the longest chain of untested links, from photons to muscle fat. H2 is the most measurable experiment, and no general rule makes more alkaline water taste better. H3 is the most practical way to build a local cuisine, and its products could also be made on Earth. The study looks for novelty in the coupled ecology, the timing of harvests and the regional food tradition, with every line descended from Earth species and the biochemistry made of familiar compounds.

## The first portfolio

These are research priorities for contained trials. The chapters hold each species' evidence and life-cycle conditions. Founders come from captive-bred stock with complete life cycles in culture, and threatened species enter only that way.

| Priority track | First tractable candidates | Wider ecological comparison | First culinary comparison |
|---|---|---|---|
| Cyclic edible fat | Jade perch (*Scortum barcoo*) and a pedigreed Nile tilapia (*Oreochromis niloticus*) reference | Freshwater clupeids such as *Limnothrissa miodon*, within their thermal and recruitment limits; plankton and prey before any predator fishery | Matched cuts and handling at several phases of the cycle: a fatty raw candidate once qualified, a fully cooked nigiri topping, plain grilled fish |
| Moderate alkaline chemistry | Nile tilapia and Mozambique tilapia (*O. mossambicus*, Vulnerable on the IUCN Red List) lines; giant freshwater prawn as an invertebrate comparison | Soda-lake *Alcolapia grahami* (Vulnerable) as a comparator of mechanism, from captive stock | Plain steamed fish and prawn first; raw fish only after qualification; flavoured recipes after the plain comparison |
| Carbonate-associated fish and shellfish | Jade perch or tilapia references; freshwater prawns and crayfish in compatible water; pearl mussels as habitat and material organisms | Barramundi, spawning in a held seawater lake or a hatchery; grazers, prey and host fish chosen together | Pale fish as a raw candidate; skin-on roast or kabayaki-style fillet; cooked crustaceans; bivalves where food use is justified |
| Eel and cooked rich flesh | Japanese eel (*Anguilla japonica*, Endangered since 2014) from captive-bred founders, with a held seawater lake as its leading spawning route (`eel_habitat_hatchery.md`); freshwater catfish as a preparation comparator | The migration between fresh sea and seawater lake, recruitment and the line's welfare | Eel for unagi; other species named and prepared kabayaki-style |

The threat statuses and their sources are in [alkaline foods](alkaline_foods.md) (A10, A11) and [carbonate shores](reef_shellfish_cuisine.md) (R15). The tracks share a reference line, so that a habitat effect can be told from a species difference. Marine herring, tuna and salmon serve as analytical and culinary comparators, in keeping with the leaning to keep saltwater life in aquariums ([register, 15 October](../../decisions.md#salt-water-and-the-fresh-seas-15-october-2026)); the held seawater lakes serve fish that live in the fresh seas and spawn in salt water ([register, 6 November](../../decisions.md#the-design-climate-at-the-target-on-the-turned-moon-its-comfort-and-the-seawater-lakes-6-november-2026)). Large predators wait for a trophic budget.

## Wild lines and the rules for what lives wild

Every line that lives wild meets the rules for what lives wild of 12 October, the fish and shellfish bred from Earth stock among them: selective breeding comes first, and a minimal engineering step then gives each wild line the shared signals, reporting traits and DNA watermark ([register, 6 November](../../decisions.md#the-food-systems-branchs-evaluation-wild-lines-from-earth-stock-the-eels-seawater-spawning-and-organisms-crossing-between-worlds-6-november-2026); the rules in the [register, 12 October](../../decisions.md#machines-in-the-living-world-12-october-2026) and [machines.md](../../../biosphere/ecology/machines.md#what-may-live-wild)). For an aquatic line the rules ask six things:

- **Climate fit:** it does best in the waters where it is placed and has no advantage everywhere.
- **Harmless Moon-wide:** it has no trait harmful at the Moon's scale and carries no pathogen of other organisms or of people.
- **Within the air's sinks:** the N2O and methane from its metabolism and its waters stay within what the air's sinks remove.
- **Closed cycles:** it returns what it takes, phosphorus above all, so a fishery's harvest and a spawning run that carries phosphorus downhill each have a return path.
- **Ready for the machines:** it carries the shared signals, reporting traits and DNA watermark.
- **Robust to evolution:** no safety feature rests on a genetic circuit that evolution can switch off.

The study applies the rules to every line that could live wild, farm lines that could escape included. Organisms on the Moon are engineered, with Earth species' traits as reference points and design variables ([register, Life and people](../../decisions.md#life-and-people)); the trials measure those traits, and breeding, or engineering where breeding stops short, changes them where the design needs it.

## Organisms crossing between worlds

Organisms will cross between Earth and the Moon over time, as they cross between continents. The design accepts the crossing, plans for it over long terms and mitigates its risks, the parasite risk of raw fish from the Moon's warm waters among them, and the DNA watermark tells the planned lines from those that crossed ([register, 6 November](../../decisions.md#the-food-systems-branchs-evaluation-wild-lines-from-earth-stock-the-eels-seawater-spawning-and-organisms-crossing-between-worlds-6-november-2026)). Each founder's pathogens and parasites, the hosts in each food web and the surveillance of raw products are therefore conditions that change as organisms arrive, and each line's effects are assessed in both directions, on the Moon and on Earth.

## What the study rests on

The habitats come from the project's products ([aquatic_biomes.md](aquatic_biomes.md)). On the forcings those products ran on, the seas' surface lies near 29.6 °C in the median and their deep water at 25–30 °C, and no fresh water falls below about 15 °C, the coolest being the deep highland lakes near 21 °C. The founders are therefore warm-water animals, and a life stage that needs cold water needs engineered water or a line selected for warmth. The light under the water sets H1's clock: just below the surface the night stays dark for 255 hours of each lunar day at the equator, 208 at 45° and 123 at 60°, and poleward of 65° the twilight never lets it darken that far. A food web's dark interval therefore depends on its latitude and depth, and fish keep feeding through it on zooplankton, benthos and detritus.

The Earth evidence supports the first link of each hypothesis, and the later links are what the trials test. Marine plankton store lipid by day and spend it at night; jade perch's muscle fat follows its diet; saline and alkaline rearing changes flavour chemistry in bass and prawns, while two prawn studies disagree on umami and a calculated umami index awaits a blinded tasting. Earth fish already vary in edible fat and in its response to feed, so choosing founders and breeding within them is the first route, with fertility, health and ecological performance selected beside growth and flavour.

Food, pearls and nacre, mineral sequestration and habitat can share one shellfish community, and each claims its own share of the biomass and minerals. A harvested mussel leaves the filter, and a shell is a lasting mineral sink only where it is kept or buried; filtration also moves organic matter into the sediment.

## Work packages and gates

| Stage | Deliverable | Gate to the next stage | State |
|---|---|---|---|
| FS0: basis and evidence | Versioned habitat basis, source-by-claim synthesis, candidate and exclusion lists | Every proposal has a mechanism, a counterexample and a stated evidence level | Synthesis written; each source read to the depth its register records |
| FS1: founder dossiers | Named populations and strains, reproduction and hatchery needs, tolerance of temperature, oxygen and ions by life stage, threat status and pathogen records | A complete life cycle and a sustainable supply of captive-bred founders, with no reliance on uncounted wild recruitment | Designed; strain acquisition and full-source extraction to come |
| FS2: environmental and food-web screen | Traces by site, time and depth; ledgers of primary production, prey, oxygen and nutrients; a chemistry matrix | Reproduction and feeding fit the whole exposure; diets close without exhausting prey or oxygen | Interfaces set; habitats read from existing products, with no coupled simulation |
| FS3: comparative pilot | Replicated line-by-environment tests with standard handling | Effect sizes and variability justify a powered trial over several generations; welfare limits met | Designed |
| FS4: breeding | Pedigree and genomic selection, reference lines, offspring tested across environments | Repeatable inherited gain at acceptable cost to health, diversity and ecology | Proposed; no heritability or generation count assumed |
| FS5: product qualification | Separate raw and cooked sensory and process evidence for each cut, fishery and harvest phase | A qualified process, repeatable quality and surveillance of the source | Proposed |
| FS6: managed ecosystems | Evidence from contained communities, a recruitment and harvest budget, an assessment of exchange between worlds, and the rules for what lives wild met | Populations and nutrients that persist at the proposed scale | Later |

The next desk work is three paired founder dossiers: jade perch against tilapia for edible fat, tilapia against prawn for moderate alkaline chemistry, and the chosen fish with a complete shellfish and prey community for carbonate habitats, beside the eel and sturgeon chapters. Each dossier is built to stay useful if its lunar advantage fails.

## Ideas from the ideation session of 4–5 November

The branch opened in an ideation session on 4–5 November. Its ideas shaped this study and are kept here as ideas; what stands is in main's register, cited beside each. The session's requests for documents became the chapters above.

| Idea (ideation, 4–5 November) | What stands in the register | Where the study takes it up |
|---|---|---|
| Work from habitats through food webs, founders and breeding to harvested tissue and cuisine, beginning with lunar-cycle fat, alkaline taste chemistry and carbonate-community foods | No row of its own; the habitats start from the waters' trophic design, nutrient-poor open water with productive habitats where nutrients gather ([register, 19 October](../../decisions.md#the-design-climate-and-the-waters-19-october-2026)) | The three hypotheses and this chapter |
| Start with Earth species and existing lines, breed intensively with pedigrees and genomics, and keep genetic engineering for needs that breeding leaves | Organisms are engineered, with Earth traits as design variables (Life and people); every wild line is bred first and then given the shared signals, reporting traits and DNA watermark by a minimal engineering step (6 November) | [Founders, breeding and trials](breeding_and_trials.md) |
| Cooked seafood as a full product category beside raw: unagi from eel, kabayaki-style fish named by species, cooked shellfish; for qualified wild catch, safety comparable to Earth's regulated sushi | Fishing can supply part of the heritage meat where the parasite risk is kept low (12 October); the parasite risk of raw fish from the warm waters is mitigated (6 November) | [Culinary evaluation](culinary_evaluation.md) |
| Earth and the Moon as connected worlds, with introductions recurring | Organisms will cross between the worlds, and the design plans for the crossing and mitigates it (6 November) | Founder, surveillance and qualification work throughout |
| Raw qualification reaching across the wider waters beyond the most controlled settings; eel grow-out in silty freshwater waterways | No row of its own; the parasite condition (12 October) and the crossing between worlds (6 November) bear on raw qualification | [Culinary evaluation](culinary_evaluation.md#raw-preparation-across-the-waters); `eel_habitat_hatchery.md` |
| Eels selected for a different spawning place, and in time for reproduction without a saline hatchery | Held seawater lakes for the eel's spawning lessen the need for an all-freshwater cycle, which stays a longer-term engineered option: a leaning, conditional on the lakes' warmth (6 November); the seawater lakes' study (6 November) | `eel_habitat_hatchery.md` |
| A complete freshwater eel life cycle whose populations persist in the Moon's Gaean waters without controlled reproduction | The seawater-lake migration leads, and the all-freshwater cycle is the longer-term engineered option (6 November) | `eel_habitat_hatchery.md` |
| A sturgeon assessment of species, habitats, beluga and its derivatives, with warm-water spawning as one part | The evaluation of 6 November records that no fresh water falls cold enough for any sturgeon to spawn, while adults could live deep in the coldest highland lakes | `sturgeon_habitat_foods.md` |

## Ideas from the review of 6 November

The review of 6 November set the older food work beside this study and drew seven ideas, none started, each a study of its own ([register, 6 November](../../decisions.md#the-food-systems-branchs-evaluation-wild-lines-from-earth-stock-the-eels-seawater-spawning-and-organisms-crossing-between-worlds-6-november-2026)): soda lakes as phosphate traps, with a spirulina food web; fish runs that carry the lakes' phosphorus back uphill; the waters' phosphorus as the measure of how much wild fish the heritage meat can hold; the lunar night as the sturgeon's winter in the highest lakes; tropical eels for the warm seas; captive lines of threatened eels and sturgeon as insurance populations; and a urea-excreting food fish for soda water.

## Interfaces with the domains

- **Biosphere and ecology:** organism traits, food webs, reproduction, community effects, the rules for what lives wild and the organisms that cross between worlds.
- **Climate, hydrosphere and illumination:** temperature, mixing, underwater light, currents and extremes, read from their products.
- **Resources:** phosphorus, essential ions, sediment, calcium and magnesium, with the food taken out, the feed brought in and the waste returned in one account.
- **Provisioning:** hatcheries, closed grow-out, feed, chilling, freezing and cooking, water treatment, the beams' power and heat, logistics and the return of waste to the land.
- **Human food culture:** regional dishes and harvest practice, grown from what the trials find.

This folder holds the synthesis and the trial designs; models of the habitats, nutrients and organisms stay in the domains that own them.
