# Eels on the Open Moon: the migration, the hatchery chain, fed eel and unagi

7 November 2026 · `domain/food-systems` · evidence synthesis, scenario accounts and research design

Eels can enter the Open Moon's waters through a migration of the Moon's own. They grow for years in the large fresh seas and their rivers, run down a one-way channel to spawn in a held seawater lake, and their glass eels climb back into the fresh sea. The hydrosphere's seawater-lakes study designed the lake and its channel, and recommends hollow 547 beside the Near-side Sea, at 37.5° N, 43.6° E on the turned Moon: a lake of 789 km², 40 m deep on average and 102 m at its deepest, fed by 25–29 m³/s of the sea's own water falling down a pool fishway, with an eel ramp beside the fishway for the climbing young (`research/studies/seawater_lakes/README.md` on `domain/hydrosphere`, f6fe4026). The register leans to this route so far as the lakes can be held and the eggs develop in their warmth, and keeps an all-freshwater life cycle as a longer-term engineered option (register, 6 November, the eel's seawater spawning).

Two conditions decide whether the migration can run. Japanese eel eggs hatch between 19 and 28 °C, best at 25–28 °C, and none hatched at 31 °C [E16]; the seawater-lakes screen set every lake at the design years' sea surface of 29.5 °C without modelling the lake, so the temperature of the lake's spawning layer is the first open check. The second is maturation: eels in captivity are matured with hormones, and the Japanese eel spawns at 150–200 m in the open ocean, a pressure that would take 0.9–1.2 km of water at lunar gravity, while the deepest water of hollow 547 presses like 17 m of Earth's sea *(computed, section 3)*. The line that runs to the lake has to mature and spawn there without hormones and without the ocean's depth, by selection first and a minimal engineering step where selection stalls.

Fed eel need a different setting from the wild run. Eaten at 1 kg of cooked eel a person a year by 8 billion people, fed eel would leave 167 kt of phosphorus a year that the fish do not retain, three to eight times the 21–56 kt a year the waters' life can spare in all, and even 0.1 kg a person leaves 16.7 kt *(computed, section 7)*. Fed eel are therefore kept in closed systems, whose sludge carries that phosphorus back to the fields, while the wild migration lives on the waters' own food web and meets the rules for what lives wild (register, 12 October and 6 November; section 6).

The Japanese eel (*Anguilla japonica*) is the reference founder. Its life cycle has been closed in captivity, its stages have been studied one by one [E04–E08], and unagi gives it a clear culinary destination [E12]. It has been listed as Endangered since 2014, a listing kept in 2020 [E17], so founders come from traceable captive-bred families and wild-caught glass eels stay out of the programme. Other *Anguilla* species enter with their own stage evidence, the tropical eels first among them (section 2).

## Scope and evidence

The chapter rests on Earth experiments cited by the IDs in the source table, on the project's products, and on scenario accounts: [eel_budget.py](eel_budget.py) reads [eel_parameters.json](eel_parameters.json) and writes [results/eel_budget.json](results/eel_budget.json), [eel_sources.json](eel_sources.json) records how far each source was read and [eel_checks.json](eel_checks.json) records the checks. Every biological rate, feed conversion and processing yield in the accounts is a scenario assumption chosen for the arithmetic. The per-person cases of section 7 read the older food study's product ([results/food_systems.json](results/food_systems.json)) for its factors: N2O at 0.8–1.9% of feed nitrogen, a recirculating farm's 3.09 kWh and 76 L of water per kg of live fish, the crop yields that price feed land, and the 21.1–56.2 kt a year of phosphorus the waters' life can spare. That product rests on the September climatology, which the turned design years at the target have replaced (ground 301.3 K, drylands on 59.5% of the land; register, 6 and 7 November), so its crop land and the waters' spare phosphorus await its refresh, while the N2O share, electricity and water come from Earth measurements and carry over. The branch's ideation session of 4–5 November and the review of 6 November contributed ideas, which section 12 lists as ideas; the decisions that stand are those of main's register. Jade perch, its fat storage and the raw preparation of fish across the aquatic habitats are treated in the [carbonate-food chapter](reef_shellfish_cuisine.md) and the [culinary evaluation](culinary_evaluation.md).

## 1. Four routes and how they rank

| Route | What it supplies | What must be achieved | Standing |
|---|---|---|---|
| R0. The hatchery chain | Captive-bred glass eels for nurseries, closed grow-out and breeding trials | Transfer the established saline hatchery methods; deliver healthy glass eels to fresh water; keep reference families | The start of every route |
| R1. The seawater-lake migration | A population that grows in the fresh seas and spawns in a held seawater lake | Maturation and spawning in the lake without hormones or ocean depth; eggs within their thermal range; larval food in the lake; glass eels finding the inflow and climbing the ramp | The leading route, the register's leaning conditional on the lake's warmth (6 November) |
| R2. Lower early-life salinity | Larval stages that need less salt | Each stage's performance curve, then inherited gains with fertility intact | A development line that serves R0 and R1: less salt in hatcheries and a wider brackish nursery |
| R3. The all-freshwater life cycle | Recruitment in the fresh waters themselves | Fertilization, embryos, larval food, metamorphosis and fertile descendants in the fresh seas' calcium-bicarbonate water | The longer-term engineered option |

R1 builds on what the Earth evidence already permits. Its eggs develop in seawater, where the Japanese eel's eggs and early embryos do best near 34–35 psu [E06]; its older larvae drift into a brackish reach of about 13–17 psu where the fishway's water spreads into the lake, since half-strength seawater gave 18.2% larval survival over two months against 8.2% in full strength [E05] and unfed larvae did best at 13–16 psu [E07]; and its glass eels, about 0.2 g each, find the inflow by its fresh water and climb the ramp. The seawater-lakes study puts the larval stage at 120–270 days, four to nine lunar days, over which a cohort spreads 4.6–21.6 km, so a lake about 20 km across or more keeps the cohort's core in open water. R1 therefore keeps the stage that has never been bred in fresh water in the medium where it works, and asks breeding for one change of behaviour: maturing and spawning in a nearby lake.

R2 and R1 proceed side by side. A line whose larvae tolerate 13 psu needs less salt in every hatchery and a wider brackish nursery in the lake, and changing where adults spawn does not wait on a freshwater embryo. R3 asks for every stage in fresh water, the eggs included, and no tested treatment has come near it: larvae in fresh water and in 10% seawater all died within six days in the short unfed comparison [E05], and the egg study's lowest salinity was 24 psu [E06]. Earth species' traits are design variables on the Moon (register, 27 September), so R3 stays open as an engineered option for the long term. Its cost goes beyond the breeding: an eel that recruits in any fresh water loses the limit the salt lake sets on where it breeds (section 6).

## 2. The warmth condition

Okamura and colleagues incubated artificially fertilized Japanese eel eggs at 19, 22, 25, 28 and 31 °C. Eggs hatched 24 hours after fertilization at 28 °C and 58 hours at 19 °C, none hatched at 31 °C, and 25–28 °C gave hatching of 76–86%, survival of 61–86% and fewer deformities than cooler water [E16]. No treatment fell between 28 and 31 °C, so the failure point lies somewhere in that interval. In the open ocean the eggs have been caught near the top of the thermocline, most at 150 m, under an upper 100 m at 26–30 °C (Aoyama et al. 2014, through the register's row of 6 November), so they develop below the warmest water.

For R1 the condition falls on the seawater lake. The seawater-lakes screen took every lake and its inflow at 29.5 °C, the design years' sea surface, between the eggs' best range and their failure. A lake's deep water settles near the coldest water its surface makes: the seas' deep water settles at 25–30 °C at 301 K, because the Moon's axis lies 1.5° from the ecliptic and gives no polar winter, while the lunar night cools the seas' surface by about a kelvin (`research/studies/cold_waters/README.md` on `domain/hydrosphere`, 1ac9357d). A seawater lake 40 m deep on average, with the fishway's fresh water spreading over it from one end, stratifies by salt as well as by heat, and its spawning layer could sit warmer or cooler than the open seas' deep water. The lake's temperature by depth through the lunar day is therefore the open check, for the hydrosphere's lake model with salinity added; this chapter runs no lake model.

If the spawning layer runs warm, two answers exist. The first is a line selected for eggs that hatch 1–2 K warmer (register, 6 November), tested in the same five-temperature design with treatments at 29, 30 and 31 °C added. The second is a tropical eel. The Indonesian shortfin eel (*A. bicolor*, Near Threatened since 2014, kept in 2020 [E17]), the giant mottled eel (*A. marmorata*) and the Celebes eel (*A. celebesensis*) live in warm rivers, and the Celebes eel spawns in Tomini Bay and the Celebes Sea, close to its rivers (Arai 2014, through the register), a short migration nearer the lunar pattern than the Japanese eel's thousands of kilometres. Their egg temperatures, their captive breeding and the status of the last two are the next sources to read.

Growth sets a second warmth question. Juvenile Japanese eels grew faster at 20 and 25 °C than at 15 °C [E02], the oxygen trials ran at 23.4–25 °C [E03] and the resting respiration was measured at 25 °C [E15], while the fresh seas' basins average 25–32.5 °C and their surface's median is 29.6 °C (the older sea-oxygen forcing, read as context). Growth and feeding at the seas' own temperatures, and the use of cooler refuges in deep water or shaded rivers, belong in the first grow-out trials.

## 3. Maturation without hormones and the missing pressure cue

In nature the Japanese eel matures during its ocean migration and spawns in the nights 2–4 days before the new moon at 150–200 m along the West Mariana Ridge (Tsukamoto et al. 2011, through the seawater-lakes study's register). In captivity it is matured by hormone treatment and spawns at surface pressure (recalled, as in the seawater-lakes study), which shows that pressure is a cue the spawning itself can do without. The cue still matters for maturation, and its size on the Moon is large: 150–200 m of seawater on Earth exerts 1.51–2.01 MPa, which at lunar gravity (1.6242 m/s²) in seawater of 1,022 kg/m³ takes 0.91–1.21 km of water, while hollow 547's greatest depth of 102 m exerts 0.17 MPa, the pressure 17 m of Earth's sea exerts *(computed: ρgh, with 1,025 kg/m³ and 9.80665 m/s² for Earth)*. The lunar night offers the darkness cue that the new-moon nights give in the ocean.

R1 therefore needs a line that matures in the fresh sea and the river, runs to the channel at the lake's own cues and spawns at a few tens of metres. Selection among captive families comes first: record for each female whether she matures unassisted, her egg quality and her fertilization, and keep assisted and unassisted spawnings apart in every record. Where selection stalls, a minimal engineering step makes the line mature in the lake (register, 6 November). Organisms on the Moon are engineered, and the founder's Earth cues are reference points to keep or replace (register, 27 September).

## 4. What the Earth evidence already permits

- **The captive life cycle is closed.** Nomura and colleagues describe the earlier closing of the Japanese eel's life cycle in captivity and measure genetic variation in size and timing at metamorphosis, with larvae reared in seawater [E04]. The reported heritability concerns metamorphic traits; it carries no estimate for freshwater reproduction.
- **Juvenile production has been repeated across facilities.** A 2026 study reports 560, 2,269 and 1,347 juveniles at three facilities that used the transferred tank technique [E08], counts from single campaigns.
- **Some larval stages do better in less salt.** Half-strength seawater gave 18.2% survival after two months against 8.2% in full strength [E05], and unfed larvae starting five days after hatch did best near 13–16 psu, with individual elements measured [E07, with its correction]. Larvae in fresh water and in 10% seawater all died within six days in the short unfed comparison [E05].
- **Eggs and early embryos need seawater strength.** Exposure across 24–42 psu changed development and favoured about 34–35 psu [E06]; hatch counts alone missed later deformities.

These results support a saline egg stage followed by brackish larval rearing, which is the nursery R1 builds into the lake. Low gravity changes egg suspension, the spread of feed particles, tank circulation and gas exchange, which each need testing before Earth tank geometry or stocking density carries over.

| Stage | Initial management | Evidence required to advance a lower-salinity or lake-spawning line |
|---|---|---|
| Broodstock and maturation | Separate reference families and reproductive conditioning | Repeated fertile gametes; cue response and assistance recorded; ancestry and reproductive contribution retained |
| Fertilization and embryos | Saline reproductive water at the eggs' temperature | Fertilization, buoyancy and suspension, normal development and hatch, judged beyond hatch counts |
| Yolk-sac larvae | Stage-specific water and handling | Normal jaw and organ development, survival and readiness to feed |
| Feeding leptocephali | Specialized larval feeding and controlled water | Ingestion, growth, sustained survival, deformities and the microbial and water-quality load through the whole stage |
| Metamorphosis and glass eel | Completion and condition measured | Functional juveniles, timing distribution and successful transfer to fresh water |
| Nursery and grow-out | Freshwater feed training, grading and habitat trials | Survival, feed conversion, health, tissue quality and no costly correlated responses |
| Reproductive descendants | Retained siblings returned to broodstock trials | Fertile offspring in the same target water across repeated generations, with failures counted |

Selection runs on split families: siblings divided across salinity, ion and temperature treatments in replicated tanks, an unselected reference line kept, and selected descendants tested in a common environment, so that family-by-environment effects are estimated before a family's good performance is read as heritable. Genomic parentage and breeding values preserve diversity and reveal trade-offs once phenotypes exist. A lower salinity is tried only after a stage performs at the preceding one, and maturation, development, larval duration, feed use, survival after transfer and cooked quality are tracked together, since an early gain that costs fertility advances no route.

## 5. Habitats for growth in fresh water

The fresh seas hold about 0.14 g/kg of dissolved salts at pH 8.3, dominated by calcium and bicarbonate, with their magnesium under study at 0.05–0.5 mmol/kg (register, 19 October, as the seawater-lakes study reads it); the same magnesium sets the seawater lake's upkeep. Larval and juvenile trials use the actual ions, alkalinity and pH, which a conductivity or salinity figure alone does not fix. Under the design air of 1.2 atm with 21.2 kPa of oxygen, water at saturation holds 7.66 mg/L of oxygen in the fresh seas at 29.6 °C, 8.94 mg/L at 21 °C and 6.31 mg/L in seawater of 35 g/kg at 29.6 °C *(computed with `biosphere/ecology/lake.py`)*. At the held seas' pH of 8.3 and 29.6 °C, 13.5% of ammonia nitrogen is un-ionized, against 0.8% at pH 7 and 89% at a soda lake's pH 10 *(computed, Emerson's relation; the share's rise with pH and temperature [E14])*, so ammonia limits fed eel more in the seas than in neutral water.

In an eight-eel pond experiment, Japanese eels used rock refuges and also burrowed into sand and mud, and filling the refuge gaps ended their use [E01]; juveniles were active at night [E02]. Silty waterways with a mosaic of refuges and feeding grounds suit eels, and sheltered channels, vegetated margins and lake edges join them. Daylight lasts 354 hours, so the eels need shelter through it, and their feeding is measured across the whole 709-hour cycle, since the long night's prey has yet to be measured.

| Habitat | Eel role | Conditions to establish before assigning production |
|---|---|---|
| Silty channels and connected backwaters | Grow-out and sheltered feeding territories | Oxygen through the full light and temperature cycle, accessible prey, sediment oxygen demand, connectivity and harvest access |
| Gravel and rock margins with crevices | Refuges within rivers, ponds and lake shores | Flow through the gaps, usable gaps, natural renewal, density and competition for prey |
| Vegetated lake margins and fresh-sea shores | The wild run's feeding grounds | Littoral production, predation effects and local oxygen minima, which deep-basin averages do not set |
| Recirculating systems | Fed grow-out, nurseries, family comparisons and quarantine | Filtration, oxygen transfer, CO2 removal, ammonia and nitrite control, reliability and sludge return (section 7) |
| Cooler upland and shaded water | Refuge from the seas' warmth and a temperature comparator | Growth, timing and the route back to the seas and the channel |
| The seawater lake and its brackish reach | Spawning, eggs, larvae and glass eels (R1) | The spawning layer's temperature, larval food, and the ramp's use by glass eels |
| Persistently hypoxic deep water or soda water | None | Direct life-stage evidence and a viable food and oxygen balance |

## 6. The wild line against the rules for what lives wild

Every line that lives wild meets the rules for what lives wild, the eels bred from Earth stock among them: selective breeding comes first, and a minimal engineering step then gives each wild line the shared signals, reporting traits and DNA watermark (register, 12 October and 6 November; [machines](../../../biosphere/ecology/machines.md)). The R1 migration is the eel line meant to live wild. Fed eel in closed systems are livestock, and the rules reach them through escapes.

| Rule | The R1 migration | What it needs |
|---|---|---|
| Fits the climate where it is placed, with no advantage everywhere | Breeds only where a seawater lake is held, and grows in connected fresh waters at temperatures beyond its founder's tested range | Growth and egg ranges matched to the lakes it is placed in; R3 would remove the salt lake's limit and needs another in its place |
| Harmless at the Moon's scale | A predator of invertebrates and small fish, and a host of eel parasites such as the swimbladder nematode carried by recruiting American eels [E11] | Founder screening, surveillance and predation trials in bounded waters before release |
| N2O, methane and ethylene within the air's sinks | Wild eels eat the food web's own production and bring no feed nitrogen; fed eel emit N2O at 0.8–1.9% of their feed's nitrogen, and the design air has no N2O sink | For the wild line, nothing beyond the food web's own cycling; fed eel stay in closed systems |
| Returns what it takes, phosphorus above all | The run carries 4 kg of phosphorus per tonne of spawners from the fresh sea down into the terminal lake; the returning glass eels carry back 0.11% of it at the culture survival of 72% from glass eel to spawner, 0.8% at 10% survival and 8% at 1% *(computed, product)* | A silver-eel harvest at the channel's head, which takes the run's phosphorus to the land; the phosphorus that reaches the lake feeds its nursery and leaves only through harvest |
| Carries the shared signals, reporting traits and a DNA watermark | Absent from Earth stock | The minimal engineering step after selective breeding |
| Its safety rests on nothing evolution can switch off | The eggs' need for salt holds breeding to the held lakes, a physiological limit | R2 and R3 select against that need, so their lines need limits that rest on more than the salt |

The run moves phosphorus in one direction. Spawners grow on the fresh seas' food web and, as anguillid eels spawn once and die (recalled), leave their bodies in a lake that stands at least a metre below the sea and has no outflow. Their young weigh 0.2 g against the parents' hundreds of grams, so the young carry back 8% or less of their parents' phosphorus for any survival of 1% or more from glass eel to spawner, and less still for spawners larger than the 250 g taken here. The phosphorus that reaches the lake stays in the waters' account and fertilises the nursery, which the waters' trophic design allows as one of the productive habitats where nutrients gather (register, 19 October). The silver-eel fishery is the return path, at the channel's head where the run gathers, and it draws on the same spare phosphorus as every wild catch (section 7). Fish that feed in the lakes and run up their inflows carry phosphorus the other way, which section 12 lists as an idea.

Organisms will cross between Earth and the Moon over time, and the design plans for the crossing and mitigates it (register, 6 November). The swimbladder nematode in newly recruited American eels [E11] shows how live eels moved between waters carry their parasites, so founders are screened, batches are traced and the waters are watched for arrivals. The DNA watermark tells the planned lines from eels that crossed.

## 7. Fed eel in the diet

Heritage meat is 1–2 kg a person a year, and wild culls and fishing can supply part of it if parasite risk is kept low (register, 12 October). The older food study counts it as boneless meat and prices it at 298,000 km² of land and 1.6–1.7 Mt of cattle methane a year at 1 kg a person for 8 billion people (`README.md`, heritage meat). The product's per-person cases set eel against that allowance, scaling the reference cohort's accounts to a population eating a set mass of cooked eel a year; whether fish count inside the allowance is open (section 12). A kilogram of cooked eel comes from 1.18–1.43 kg of prepared raw cut across the three processing examples, so the comparison with the boneless basket is close to like for like.

| People, cooked eel a person a year | Live harvest | Dry feed | Feed phosphorus | Phosphorus not retained | Times the waters' spare 21.1–56.2 kt | N2O, kt N a year | Closed systems' mean power | Feed crop land | Open-pond methane |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 billion, 0.1 kg | 0.77 Mt | 1.14 Mt | 11.4 kt | 8.4 kt | 0.15–0.40 | 0.7–1.7 | 0.27 GW | 370–850 km² | 0.002–0.48 Mt |
| 8 billion, 0.1 kg | 1.54 Mt | 2.29 Mt | 22.9 kt | 16.7 kt | 0.30–0.79 | 1.4–3.3 | 0.54 GW | 730–1,690 km² | 0.005–0.97 Mt |
| 4 billion, 1 kg | 7.69 Mt | 11.4 Mt | 114 kt | 83.5 kt | 1.5–4.0 | 7.0–16.7 | 2.7 GW | 3,660–8,470 km² | 0.025–4.8 Mt |
| 8 billion, 1 kg (or 4 billion, 2 kg) | 15.4 Mt | 22.9 Mt | 229 kt | 167 kt | 3.0–7.9 | 14.0–33.3 | 5.4 GW | 7,320–16,900 km² | 0.05–9.7 Mt |
| 8 billion, 2 kg | 30.8 Mt | 45.7 Mt | 457 kt | 334 kt | 5.9–15.8 | 28.1–66.7 | 10.9 GW | 14,600–33,900 km² | 0.10–19.4 Mt |

The table uses the reference processing example; the efficient and higher-resource examples leave 100 and 325 kt of phosphorus unretained at 8 billion and 1 kg, against the reference's 167 *(product)*. The waters' spare phosphorus is the older study's reading of resources' sea-phosphorus study, 30% of each water's dissolved input as a design guess for the largest harvest that leaves the food webs whole. Released to open water, the phosphorus that fed eel leave behind would exceed that whole allowance at any portion above about 0.13–0.34 kg a person at 8 billion, and at 0.1 kg it would take 30–79% of it, leaving little for any wild catch. Fed eel are therefore kept in closed systems: recirculating tanks whose solids and effluent water are captured and treated, with the sludge returned to the fields, as the older study returns aquaponic sludge. The feed's 229 kt of phosphorus at 8 billion and 1 kg comes from the land each year; 61.5 kt of it leaves in the harvested fish and the rest returns through the sludge, so with full capture the land's net draw is the phosphorus in the harvest, which then passes through processing and people to the project's recycling.

Against heritage meat the closed systems are light on land and power. The feed's crops take 7,300–16,900 km² at 8 billion and 1 kg, 2.5–5.7% of the heritage basket's 298,000 km², because that basket is largely grazed; the older study's finding that aquaponic fish take four to six times the cropland of the animal part they replace compares them with ecology's crop-fed animal part, a different baseline. The electricity, 5.4 GW on average at the older study's recirculating farm with its water held near 30 °C, is small against the 28–152 TW that stay on the Moon (register, 5–6 November), and the beams supply it as they supply the rest of the food, which is made on the Moon. N2O accumulates, since the design air has no sink for it: the fish emit 14–33 kt of N2O nitrogen a year at 8 billion and 1 kg, 0.005–0.011 ppb a year and 5–11 ppb over a thousand years, beside the fields' 0.045–0.17 ppb a year for the whole diet at 8 billion (the older study's finding 9); the feed crops' own N2O comes on top. A closed system's methane depends on how its sludge is treated, which the older study does not price; digesting the sludge and burning the captured gas is the treatment to price first. Open ponds show the range the closed systems avoid: Earth's carp ponds emit 3.2 kg of methane per tonne of fish and its semi-intensive bream ponds 630 kg (Yuan et al. 2019's paired fluxes and yields, through the older study's register, each flux held for a year), which at 8 billion and 1 kg spans 0.05–9.7 Mt a year against the heritage cattle's 1.6–1.7 Mt.

The wild run supplies eel on the waters' own terms. A silver-eel catch at the channel draws on the same 21–56 kt of spare phosphorus as every wild fish, and whole wild fish at 0.4% phosphorus could reach at most 5.3–14 Mt a year on it, 0.66–1.75 kg a person at 8 billion (computed in the review of 6 November, on a design-guess harvest share; section 12). Whether fish count inside the 1–2 kg of heritage meat, and how the wild share divides among eels, other fish and the sea-feeding flyers, are open.

## 8. Juvenile supply, feed and nutrients per tonne

For live harvest mass H, mean harvest mass m, grow-out survival s_g and nursery survival s_n, the harvest needs N_glass = (H/m)/(s_g s_n) glass eels; with s_l the survival from normal hatchling to glass eel, h the normal hatch per fertilized egg and f the fertilization success, it needs N_eggs = N_glass/(s_l h f) eggs. The reference takes 250 g harvest fish, 10 g at nursery exit, 0.2 g glass eels, 90% grow-out and 80% nursery survival, 1% survival from hatchling to glass eel, 60% normal hatch and 70% fertilization, all accounting inputs that await measurement.

| Stage | Expected count per 1,000 kg live harvest | Live biomass |
|---|---:|---:|
| Eggs before fertilization | 1,322,751 | |
| Fertilized eggs | 925,926 | |
| Normal hatchlings | 555,556 | |
| Glass eels entering the nursery | 5,556 | 1.11 kg |
| Juveniles entering grow-out | 4,444 | 44.44 kg |
| Harvested eels | 4,000 | 1,000 kg |

The larval stage is the bottleneck. With the later stages held fixed, survival from hatchling to glass eel of 0.1%, 1% and 5% takes 13.23 million, 1.32 million and 264,550 eggs, while the feed after the glass-eel stage stays the same for the same delivered juveniles. For throughput, 1,000 glass eels per successful batch, 80% of batches succeeding and 270 occupied larval days plus 30 for turnaround give 974 glass eels per batch position a year, so 1,000 kg of harvest a year needs 5.70 positions, six before redundancy; output, cycle time and success have not been shown together, and failed batches still consume feed, water and time.

Larval and grow-out diets are separate problems, and a palatable juvenile pellet does nothing for leptocephalus feeding; larval feed, ingestion, waste, water quality and survival are recorded per batch. For grow-out, the economic dry feed conversion ratio is dry feed supplied over net stage gain, which already includes lost feed and mortality. A 2025 trial at 24.4–24.7 °C compared paste and extruded diets under recirculating and flow-through conditions; pellets did well, but the diets differed in composition and paste intake was uncertain [E09]. Nursery FCR 1.2 and grow-out FCR 1.5 give 52 kg of nursery feed and 1,433 kg of grow-out feed, 1,485 kg of dry feed per tonne harvested. These sit inside the older study's tilapia range of 1.1–2.24; the eel feed's 1% phosphorus lies inside that study's 0.71–1.41%, and its 48% protein above that study's tilapia feed, as a carnivore's feed does.

| Element, kg per 1,000 kg live harvest | Feed input | Glass-eel input | Retained in harvested whole fish | Not retained |
|---|---:|---:|---:|---:|
| Nitrogen | 114.07 | 0.03 | 28.80 | 85.31 |
| Phosphorus | 14.85 | 0.004 | 4.00 | 10.86 |

The balance assumes feed of 48% crude protein and 1% phosphorus, whole fish of 18% protein and 0.4% phosphorus at entry and exit, and 6.25 kg of protein per kg of nitrogen. The unretained part is dissolved excretion, faeces, uneaten feed and mortalities, which treatment divides among captured solids, recovered nutrients, water and gaseous nitrogen; nitrification changes nitrogen's form and leaves it in the system. In open habitats, measured prey consumption and prey renewal replace the feed account, and a prey fish eaten by an eel is counted once.

## 9. Oxygen, water and salt in hatcheries and closed systems

Oxygen capacity follows the fish held and the feed supplied on the most demanding day. The standing-stock screen keeps three demands apart: fish respiration, O_fish = Bq(24)/10⁶; nitrification, 4.57 kg of oxygen per kg of ammonia nitrogen oxidized, with 7.14 kg of CaCO3-equivalent alkalinity consumed [E10]; and other biological demand per kg of feed. Its respiration rate of 200 mg O2/kg/h carries no temperature: the eel's resting rates were measured at 25 °C [E15], and a young-eel trial at 23.4–25 °C found poorer growth and feeding below about 40% saturation with comparable performance from 60% to above 100% [E03], so the rates for fed fish at 29–32 °C remain to be measured.

| Quantity, one tonne of fish held | Reference assumption or result |
|---|---:|
| Dry feed supplied | 10 kg/day |
| Fish respiration (temperature unassigned) | 200 mg O2/kg/h, 4.80 kg O2/day |
| Ammonia nitrogen at 0.03 kg N/kg feed, all nitrified | 0.30 kg N/day, 1.37 kg O2/day |
| Other biological demand at 0.20 kg O2/kg feed | 2.00 kg O2/day |
| Total demand | 8.17 kg O2/day |
| Capacity with an assumed factor of 1.5 | 12.26 kg O2/day, 0.511 kg/h |
| Alkalinity consumed by nitrification | 2.14 kg CaCO3-equivalent/day |

Respiration of 100, 200 and 300 mg/kg/h paired with 5, 10 and 15 kg of feed a day and 0.1, 0.2 and 0.3 kg of other demand per kg of feed gives capacities of 5.38, 12.26 and 20.63 kg O2/day, paired examples for one tonne held. Carried in flowing water from an inlet of 7 mg/L to an outlet floor of 5 mg/L, the reference capacity takes about 255 m³/h. That inlet is 91% of saturation in the fresh seas at 29.6 °C and above saturation in seawater of 35 g/kg, which holds 6.31 mg/L at saturation, so a seawater hatchery either adds oxygen or runs a lower inlet and a narrower drop. A connected 100 m³ holds 0.20 kg of oxygen above the floor, 23.5 minutes at the capacity demand if perfectly mixed, so oxygen supply and circulation need backup on the measured response time. Lower gravity leaves the biochemical demand unchanged and alters bubble transfer, circulation, settling and degassing, which are tested at the actual pressure and gravity before oxygen is converted to power.

| Salt concentration, 100 m³ inventory losing 1 m³/day as liquid with no recovery | Salt in inventory | Make-up replacing liquid losses |
|---|---:|---:|
| 35 g/kg | 3,500 kg | 35 kg/day |
| 17.5 g/kg | 1,750 kg | 17.5 kg/day |
| 13 g/kg | 1,300 kg | 13 kg/day |
| 0.14 g/kg | 14 kg | 0.14 kg/day |

The salt table is a physical comparison for a hatchery unit; the 0.14 g/kg row assigns no egg or larval success. Evaporation removes water and keeps the salts, recovery lowers make-up at a cost in treatment and energy, and a stage-separated hatchery with a small saline egg unit and brackish larval units captures much of the saving short of R3. Salt, water, feed, oxygen, energy and labour are compared per healthy glass eel delivered and per kg of cooked fish, with failed batches counted. The held lake carries salt on a larger scale: 1.15 Gt to fill hollow 547 and an upkeep from nothing to about 3 kt a year per km² of a 20 m lake depending on the feed water's magnesium (seawater-lakes study).

## 10. Cooked products, yield and breeding value

The culinary target is eel for cooked unagi, with separate trials of shirayaki, steamed-then-grilled kabayaki and direct grilling; the Japanese agriculture ministry describes both the steamed Kanto sequence and the Kansai grilling without steaming [E12]. The traits to select and measure are skin and fillet integrity, tenderness after steaming or grilling, fat retention, the effort of boning, flesh thickness, lipid distribution, aroma, resistance to oxidation, juiciness and reheating quality. Fat that renders away raises feed cost without adding to the valued serving. Tasting compares matched ages and sizes, blinded, sauced and unsauced, so that the glaze hides neither lineage nor habitat. Actual eel keeps the unagi name; catfish or jade-perch fillets cooked the same way are labelled by species as kabayaki-style.

Live weight, prepared skin-on deboned cut, cooked fish before sauce, rendered fat and other co-products are weighed separately. The reference tonne yields 650 kg of prepared cut and 350 kg of other material, and cooking leaves 520 kg of fish after 130 kg of loss: 5,200 portions of 100 g, about 130 g from each 250 g eel, with sauce, rice and the rest of the meal on top.

| Processing example | Grow-out dry FCR | Prepared cut of live | Cooked fish of cut | Cooked fish per tonne live | Dry feed after glass eel per kg cooked |
|---|---:|---:|---:|---:|---:|
| Efficient | 1.2 | 75% | 85% | 637.5 kg | 1.88 kg |
| Reference | 1.5 | 65% | 80% | 520 kg | 2.86 kg |
| Higher-resource | 2.0 | 55% | 70% | 385 kg | 5.10 kg |

Nursery FCR stays 1.2 in all three, and the pairs illustrate resource consequences until weighed processing trials replace them. Cooking widens the eel's culinary uses, while toxins and chemical contaminants from the source water still need source controls, since freezing and cooking address particular hazards [E13]; the thermal process is validated for the actual cut, and handling after cooking is controlled. The programme's eel products are cooked, and raw products of any fish follow the [culinary evaluation](culinary_evaluation.md).

## 11. The development sequence

1. **Establish the hatchery chain (R0).** Acquire traceable captive-bred founder families, reproduce saline hatchery success and deliver healthy glass eels to measured freshwater nursery and closed grow-out trials, recording failed batches and the full input boundary.
2. **Select for maturation in the lake (R1).** Rank families on unassisted maturation, egg quality and spawning at the lake's depth and darkness, beside the stage-by-stage salinity work of R2; add the minimal engineering step where selection stalls.
3. **Check the lake's warmth.** Obtain the seawater lake's temperature by depth, and test eggs at 28–31 °C; turn to a warmer-egg line or a tropical eel if the spawning layer runs warm.
4. **Run a bounded migration.** In one held lake and its fresh sea, follow cohorts through spawning, larval feeding in the lake, the ramp and the return, with no supplementary stocking; count locally born daughters reaching reproduction per breeding female, which must reach at least one at the tested density, and keep parentage records so that stocking or arrivals cannot pass as recruitment.
5. **Weigh the plate.** Pair the closed grow-out groups with standardized culinary trials, replace the assumed yields and estimate feed, oxygen and phosphorus per kg served.
6. **Keep R3 as an engineered option.** Hold the all-freshwater life cycle for the long term, with the limit on spread it would need settled before any line recruits in open fresh water.

## 12. Open checks and ideas

Open checks: the seawater lake's temperature by depth through the lunar day; egg hatching between 28 and 31 °C; the tropical eels' egg temperatures, captive breeding and status; maturation without hormones; growth and respiration at 25–32.5 °C; larval food in the lake and the eel ramp at lunar gravity; a closed system's sludge methane; the feed crops' N2O and the land of feed made from power; and whether fish count inside the heritage-meat allowance.

Ideas from the review of 6 November, none started: fish that feed in the lakes and run up their inflows to spawn carry phosphorus back toward the land, so lake-to-river runs return what the eel's run takes downhill; the waters' spare phosphorus as the measure of how much wild fish heritage meat can include; and captive lines of threatened eels, the Critically Endangered European eel and the Endangered Japanese eel [E17], kept in the held lakes as insurance populations. Ideas from the branch's ideation session of 4–5 November: a richness in unagi that follows the lunar cycle, through fat stored across the long day and night (register, July: life stores resources across the long day and night), tested across full cycles; and regional unagi traditions from local prey, husbandry and harvest timing.

## Sources and access

Primary studies and official accounts were read between 4 and 7 November 2026; the [source register](eel_sources.json) marks which were read in full and which at abstract level, and records the references read through other registers. No source concerns lunar breeding, lunar carrying capacity or the scenario yields.

| ID | Source and use |
|---|---|
| E01 | [Rock refuges and Japanese eel](https://doi.org/10.3825/ece.22.83), 2019; small pond telemetry experiment, abstract |
| E02 | [Juvenile growth and activity by temperature](https://doi.org/10.1111/j.1095-8649.2003.00211.x), 2003; abstract |
| E03 | [Dissolved oxygen and young-eel growth](https://doi.org/10.2331/suisan.52.597), 1986; abstract, with exposure conditions |
| E04 | [Captive-bred metamorphic-trait genetics](https://doi.org/10.1371/journal.pone.0201784), 2018; relevant full-text methods and results |
| E05 | [Low-salinity leptocephalus rearing](https://www.sciencedirect.com/science/article/pii/S0044848609007509), 2009; indexed publisher and FAO abstracts |
| E06 | [Salinity and early developmental deformities](https://www.sciencedirect.com/science/article/pii/S0044848609003366), 2009; indexed publisher abstract |
| E07 | [Elemental requirements of unfed larvae](https://doi.org/10.1155/are/6874048), 2026; relevant full-text results and [correction](https://doi.org/10.1155/are/9780765) |
| E08 | [Juvenile production repeated at three facilities](https://doi.org/10.2331/suisan.25-00038), 2026; final publisher abstract |
| E09 | [Feed types and culture systems](https://doi.org/10.3390/ani15162420), 2025; full-text design, outcomes and limitations |
| E10 | EPA, [Nutrient Control Design Manual](https://www.epa.gov/sites/default/files/2019-02/documents/nutrient-control-design-manual.pdf), 2010; nitrification stoichiometry |
| E11 | [Parasites in newly recruited American eels](https://academic.oup.com/icesjms/article/73/1/14/2458698), 2016; indexed primary abstract |
| E12 | MAFF, [Unagi no Kabayaki](https://www.maff.go.jp/e/policies/market/k_ryouri/search_menu/6616/index.html); official preparation description |
| E13 | FDA, [Fish and Fishery Products Hazards and Controls](https://www.fda.gov/food/seafood-guidance-documents-regulatory-information/fish-and-fishery-products-hazards-and-controls); parasite, harvest-area and natural-toxin guidance |
| E14 | UF/IFAS, [Ammonia in Aquatic Systems](https://edis.ifas.ufl.edu/publication/FA031/pdf), revised 2022; ammonia, pH, temperature and nitrification passages |
| E15 | [Oxygen consumption of pond-cultured eel](https://www.jstage.jst.go.jp/article/jji1950/7/2-4/7_2-4_49/_article), 1958; resting rate by size at 25 °C, publisher abstract |
| E16 | Okamura et al., [water temperature and the early development of Japanese eel](https://doi.org/10.1111/j.1444-2906.2007.01461.x), Fisheries Science 73, 2007; publisher abstract |
| E17 | Kaifu et al., [the IUCN Red List and threatened freshwater eels in Japan](https://doi.org/10.1371/journal.pone.0331511), PLOS ONE, 2025; full text, read for the Red List history |

## Reproducing the accounts

From the repository root:

```sh
python -m research.studies.food_systems.eel_budget
python -m pytest research/studies/food_systems/test_eel_budget.py -q
```

The product (schema `terluna.research.eel-budget/2`) carries its assumptions, the sha256 of its producer, parameters and source register, the sha256 of the older food study's product it reads, the shared year-length constant, and the older study's factors as read. The tests check the arithmetic, the mass balances, the per-person cases by hand, the boundaries and the product's reproduction; the biological assumptions await trials.
