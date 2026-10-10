# Guards and reciprocal partnerships on sky reefs

A sky reef could support several kinds of defenders alongside partners that retain nutrients, maintain surfaces and import food. Their value depends on which threats they prevent and which resources they return. A kilogram of occupants is a lifting load; it becomes a beneficial partnership only when its contribution offsets its food, water, shade, attachment and maintenance costs.

The most promising arrangement is a mosaic: small patrols near vulnerable tissue, larger hunters in sheltered spaces, mobile defenders around nesting sites, recycling communities in protected wet pockets, and externally feeding visitors at peripheral roosts. These are candidate functions, not selected organisms. Earth experiments establish some constituent relationships; they do not establish an assembled aerial community.

The earlier [tenant screen](tenants.md) supplies geometry and conditional budgets. The calculations in [symbioses.py](symbioses.py) keep collective live mass, individual mass, response distance and food allocation separate. Their [scenarios](symbioses_scenarios.json) are explicit illustrative choices. The [product](results/symbioses.json) records units, assumptions and provenance; the [source register](symbioses_sources.json) distinguishes accessible papers from abstracts and indexed passages.

## What the guard mass means

The original 1–10 g/m² entry is a design guess per square metre of reef footprint. It has no wet or dry qualification and was never a demonstrated defense requirement or upper limit. The new scenarios explicitly use **live mass**, and examine 1, 10 and 100 g/m² without claiming that any value is sustainable.

At 10 g/m², allocating all the mass to one size would give ten 1-g individuals per square metre, one 1-kg individual per 100 m², or one 100-kg individual per hectare. Mixtures share the same allowance. A circular 2-km reef would contain about 3.14, 31.4 or 314 tonnes at the three scenario densities. These are unit conversions, not population predictions.

The curved exterior is about 3.63 times the footprint in the earlier geometry. Uniformly spreading 10 g per footprint m² over that exterior gives about 2.76 g per surface m². Actual defenders would concentrate near nests, growing tissue, wounds and reproductive structures. Brood, resting animals and reproductive individuals also weigh something; total biomass cannot all be treated as active patrols.

Individual size, collective density and spatial deployment therefore need independent choices. Larger guards have fewer bodies at a fixed mass allowance. They may travel farther or deter larger intruders, but neither improvement follows from mass alone.

## Threats and complementary defenders

| Threat or loss | Candidate partner and service | What must work |
|---|---|---|
| Eggs, small grazers and organisms entering folds | Small social patrols remove or discourage them near sensitive host tissue | Reliable detection, access to crevices, recruitment and food that supports reproduction |
| Larger crawling grazers and concealed larvae | Resident hunters occupy junctions and sheltered margins | Access to prey without puncturing the host, and coexistence with small patrols |
| Flying herbivores and egg-layers | Flying predators hunt arrivals; social flyers could recruit around disturbances | Interception under the reef's wind, affordable flight and discrimination among visitors |
| Large browsers or animals stripping tissue | Swarming deterrents and territorial nesting flyers make feeding costly | Vulnerable attacker anatomy, timely response and an attacker willing to withdraw |
| Fouling of collectors or exposed exchange surfaces | Surface-maintenance animals remove obstructing material | Selective cleaning and delivery of waste to a catchment rather than overboard |
| Nutrient export in emerging insects | Resident pool predators retain part of the prey production | Their own residence, oxygen demand, diet and eventual nutrient return |
| Microbial disease and tissue decay | Candidate microbial competitors or selective sanitation | A demonstrated disease benefit; mere cleaning or microbial presence does not establish it |

**Small defenders can affect very large browsers.** Goheen and Palmer's ant removals and feeding trials linked ant abundance to elephant damage on acacias. Elephants avoided ant-occupied branches when alternative forage was available. Attacks on sensitive tissue were a proposed explanation. This supports deterrence despite a large size mismatch, with clear limits: it does not establish effectiveness against an insulated mouth, a rapid bite, or a determined hungry sky grazer. [Goheen and Palmer 2010](https://doi.org/10.1016/j.cub.2010.08.015).

**Larger defenders need not replace smaller ones.** In coral experiments, small guard-crabs defended against small predatory snails; intermediate and large crabs protected against different sea stars. Species as well as size affected performance. A mixture can cover gaps left by one kind, but complete protection by an assembled mixture was not demonstrated for every predator. [McKeon and Moore 2014](https://doi.org/10.7717/peerj.574).

For a reef, ant-like patrols and flying hunters are stronger starting candidates than a single universal guard. Larger territorial nesters are a further possibility: defending a nest can incidentally protect the host. That benefit would be local and seasonal unless another mechanism maintains it. Feeding elsewhere could support their larger bodies, while leaving landing loads, roost damage and nutrient delivery to be budgeted separately.

## Response distance and the first damaging bite

A first bite can occur before defenders arrive. A useful requirement is

`detection delay + travel time + deterrence time < time to unacceptable damage`.

The right-hand side depends on the threatened function. Removing renewable green tissue may be tolerable; opening a gas compartment, destroying a reproductive organ or severing a connection has a different consequence. No current number establishes a universal safe bite size or attack duration.

For an illustrative square grid of nests spaced `s` metres apart, the farthest point in a cell lies `s / sqrt(2)` from its nest. The screen calculates `delay + path factor × distance / speed`. It assumes a nest is staffed and connected to the target. The path factor represents indirect travel, not successful detection or combat.

| Nest spacing | At 0.02 m/s | At 0.1 m/s | At 2 m/s |
|---|---:|---:|---:|
| 5 m | 4.92 min | 1.38 min | 32.7 s |
| 20 m | 18.18 min | 4.04 min | 40.6 s |
| 100 m | 88.89 min | 18.18 min | 83.0 s |

All entries use a guessed 30-second activation delay and path factor of 1.5. Speeds and nest spacings are independent scenario inputs, not assigned Earth traits for future organisms. For an aerial route, wind-relative and ground-relative speed would need distinction; for a crawler, folds, gaps and wet surfaces could lengthen or block the route. The table stops at arrival and excludes deterrence time.

At these assumptions, arrival within five minutes requires square-grid nest spacing no greater than about
5.1 m at 0.02 m/s, or 25.5 m at 0.1 m/s. A ten-second deadline is impossible with the assumed 30-second activation,
even at zero travel distance. These are geometric requirements conditional on occupied nests, not occupancy estimates.

Earth response observations need their own interpretation. An ant-plant study measured local responses to disturbed
nest chambers in seconds, with response probability also varying sharply between settings. A social-wasp experiment
measured successful prey capture in minutes inside a small screenhouse. Neither supplies a fitted kilometre-scale
reef response constant. [Youngsteadt et al. 2024](https://doi.org/10.1002/ecy.4449),
[Southon et al. 2019](https://doi.org/10.1098/rspb.2019.1676).

This makes distributed small nests attractive even on a giant colony. One central nest does not serve a kilometre of skin quickly. Mobile defenders could connect neighbouring modules, while local patrols protect organs that cannot wait. Passive tissue resistance, wound closure and compartmentalisation remain necessary candidates wherever damage can precede recruitment. None of these guard response times is the colony's pressure-trim or landing response time.

## Food quality and the shared reward budget

The reference host allocates 0.1 kg C per footprint m² per year to tenants. This is one shared allowance. It cannot be assigned in full to ants, flying defenders, cleaners and pool fauna independently. Food from intercepted visitors or foraging elsewhere must be recorded at the reef boundary; predation on residents transfers food already in the local web.

For each guild the screen calculates annual replacement production from live biomass, dry fraction, tissue composition and turnover. Ingested carbon is production divided by assimilation and population production efficiency. The latter includes the respiratory spending implicit in that efficiency; adding a second complete maintenance allowance would double-count it. These efficiencies are hypothetical and particularly unsuitable as an unexamined substitute for the energy expenditure of a large warm-bodied flyer.

The following sensitivity fixes dry mass at 30% of live mass, carbon at 50% of dry mass, four annual biomass
turnovers, carbon assimilation at 80% and population production efficiency at 30%. Phosphorus is 0.795% of dry
mass with 65% assimilation, midpoints of the earlier study's insect-based ranges. These are illustrative tissue
and feeding assumptions, not a fitted diet or a selection of organism sizes.

| Live guard mass per footprint m² | Food C required per m² per year | Share of entire host tenant-food C budget | Dietary P required per m² per year |
|---|---:|---:|---:|
| 1 g | 2.5 g | 2.5% | 0.0147 g |
| 10 g | 25 g | 25% | 0.1468 g |
| 100 g | 250 g | 250% | 1.4677 g |

At these assumptions, 100 g live/m² exceeds the entire reference tenant-food carbon allowance. At one annual
turnover it would require 62.5% instead. Holding production efficiency fixed while changing turnover also changes
the implied respiratory expenditure, so this comparison cannot establish that a slowly reproducing animal has
low maintenance. Both parameters require independent calibration for each guild.

The mixed example divides 10 g live/m² among a small patrol, mobile defenders and larger nesters. It allocates
24% of the one host-food budget, uses 14.25%, and stipulates external food containing another 7 g C/m²/year. Carbon
requirements close to numerical precision, but the small patrol still lacks 15.1 mg of assimilable P/m²/year for its
prescribed replacement. These mass shares and supplies are scenarios, not a demonstrated available food web.

External food carries paired carbon and phosphorus at the chosen composition; only the consumed fraction's
phosphorus receives credit. Host sugar fills the remaining carbon requirement within each guild's allocation.
Unused supplies and assimilation losses remain explicit. Additional feeding solely to obtain phosphorus, with
surplus carbon disposed of, is not resolved. The candidate allocations replace portions of the tenant community;
they are not free extra residents added on top of the host's existing mass and food account.

Phosphorus demand is kept separately. Sugar can supply energy but cannot, by itself, supply the phosphorus for new animal tissue. Food bodies, prey, microbial foods or another nutrient-bearing diet are necessary. Guano delivered to a reef is not automatically edible or assimilable by every guard: decomposers, plants or prey may mediate its transfer.

The calculation is a steady replacement screen. It does not pay an unrecorded initial colonisation cost, guarantee a daily ration during darkness, resolve all essential nutrients or establish a stable age structure. The food quality of real ant rewards is consequential: nutrient supply limited food-body production in a Malaysian ant-plant experiment; colony size correlated with host food production, and more workers were associated with less leaf damage. [Heil et al. 2001](https://doi.org/10.1007/s004420000534).

A larger reward budget might be possible for a mature reef, but it competes with expansion, repair, replacement gas, sailing and extra housing. Protection should be judged by avoided damage and host performance at its cost. A dense guard colony could defend well and still cost its host more than the damage it prevents. An ant-exclusion
study reported that reducing some partners improved host growth or reproduction despite their protective service;
rare catastrophic damage may change the benefit over longer periods. This cost warning is based here on its
accessible abstract, not a transferred numerical cost. [Stanton and Palmer 2011](https://doi.org/10.1890/10-1239.1).

## Partnerships that keep the reef supplied

| Partnership | Potential host benefit | Partner's return and principal cost |
|---|---|---|
| Guards in chambers beside absorptive tissue | Reduced damage; recovery of nutrients from nest debris | Shelter and nutrient-bearing food; chamber mass and host investment |
| Hunters using fringes and crevices | Removal of some grazers; fertilisation from externally caught prey | Hunting sites and refuge; possible predation on other useful residents |
| Decomposers and pool predators | Conversion and retention of captured litter and animal waste | Organic substrate and water; respiration, oxygen demand and eventual export |
| Nitrogen-fixing films | Conversion of atmospheric N₂ into reactive nitrogen | Hydrated habitat, energy, phosphorus and trace elements; possible shading |
| Externally feeding roosting flyers | Import of nutrients from land or sea | Roosts and reproductive sites; concentrated loads, fouling and wet mass |
| Surface-maintenance partners | Keeping specialised collectors or exchange surfaces functional | Food and shelter; cleaning can remove useful microbes or export nutrients |
| Gas-management microbes | Candidate removal of oxygen entering a gas compartment | Hydrogen or another energy supply and controlled habitat conditions |
| Bud carriers and nursery associates | Dispersal and establishment of the next generation | A reward or habitat; reliable delivery of both host and partner propagules |

Ant-occupied *Dischidia* leaves provide a precedent for chambers near absorptive host tissues. Isotope estimates attributed part of the leaves' carbon and nitrogen to ant-associated respiration and debris. Those particular fractions do not become reef uptake rates, and recapturing respired carbon recycles carbon already fixed. [Treseder et al. 1995](https://doi.org/10.1038/375137a0).

Predators can fertilise hosts: a bromeliad-associated spider supplied nitrogen and improved plant growth in experiments. Whether that fertilisation imports nutrients depends on where the prey fed. The focal terrestrial bromeliad does not itself establish the proposed reef tank habitat. [Romero et al. 2006](https://doi.org/10.1890/0012-9658(2006)87%5B803:BSIHPN%5D2.0.CO;2).

Pool predators can also change what leaves. In a bromeliad food-web experiment, detritivores increased host nitrogen uptake when predators prevented nutrient-bearing prey from emerging and departing. This is a reason to examine residence time and emergence, not to assign every predator a fixed nutrient-retention benefit. [Ngai and Srivastava 2006](https://doi.org/10.1126/science.1132598).

Nitrogen fixation and nutrient recycling do different jobs. Epiphytic fixation on pelagic *Sargassum* varies substantially, with little or no dark fixation in the examined observations. It supports the possibility of diazotrophic associates, while providing no transferable sky-reef rate. Fixation adds reactive nitrogen from atmospheric N₂; it adds no phosphorus. [Johnson et al. 2023](https://doi.org/10.1371/journal.pone.0289485).

Sea-feeding flyers could import phosphorus across the reef boundary. The relevant factors are external feeding, nutrient intake, excretion at the roost, retention and subsequent availability, with export subtracted. The seabird literature supports concentrated deposition at colonies, including harmful enrichment and disturbance. Its global phosphorus figure is a modelled deposition at breeding colonies, not a guarantee for each flyer or reef. [Otero et al. 2018](https://doi.org/10.1038/s41467-017-02446-8).

## Water and maintenance can create conflicts

The proposed 1–5 L/m² of pools share the host's free-water store. Water assigned to ecological habitat cannot simultaneously be counted as freely disposable ballast without counting the organisms and nutrients lost with it. Persistent pool communities require a residual volume, water quality and night oxygen supply. Humus and epiphytes may retain moisture but also add wet mass and their own transpiration; the sign of their water benefit remains open.

Cleaning should preserve useful functions rather than sterilise every surface. Ants maintained specialised capture surfaces in a pitcher-plant experiment, restoring performance after contamination; some removed material was deposited outside the pitcher. A reef would benefit from directing such material into recycling pockets where possible. This precedent establishes surface maintenance, not antimicrobial protection or repair of a pressure-bearing membrane. [Thornham et al. 2012](https://doi.org/10.1111/j.1365-2435.2011.01937.x).

The [storm screen](../aerophytes/storms.md) already includes a hypothetical oxygen-scavenging lining, costing hydrogen through `2 H₂ + O₂ → 2 H₂O`. Its 4–18% additional hydrogen is relative to the model's permeation loss. Stoichiometry does not establish uptake kinetics, containment, night function or protection during a tear. The proposed guanine-like diffusion barrier is a separate host-material mechanism; a microbial film has not been shown to replace it.

Host dormancy also leaves other organisms with their own demands. Active guards, decomposers, pool fauna and gas-management microbes may continue consuming carbon or oxygen. An annual surplus cannot establish survival through the worst dark interval. Defensive coverage could change with light, temperature and activity, even when the year's food budget closes.

## Coexistence and inheritance

Specialised host modules are division of labour within the colony. Animals, microbes and epiphytes living on them are other organisms. Their residence may be mutualistic, neutral, parasitic or dependent on circumstances. The label guard should describe demonstrated or proposed behavior against a specified threat, not confer a benefit on every predator.

Several species may coexist by using different chambers, prey, body-size ranges or activity periods. They may also eat or displace one another. Spatially separated rewards could keep guards away from reproductive organs and useful visitors; experiments on extrafloral nectar show that reward placement can alter ant occupation of flowers and reproductive outcomes. [Villamil et al. 2019](https://doi.org/10.1111/1365-2745.13135).

Different partners may be useful at different host ages. Long-term acacia observations and demographic modelling found that sequential associations could give a greater predicted lifetime benefit than any single ant species. This supports succession across a reef's life, not peaceful simultaneous occupation by all those ants. [Palmer et al. 2010](https://doi.org/10.1073/pnas.1006872107).

Tolerance is also not proof of mutualism. Work on competing acacia ants distinguished host-monopolising
mutualists from an exploiting species whose unrelated workers could share a host. A reef-wide supercolony or
peaceful mixed-ant patrol needs its own mechanism. [Kautz et al. 2012](https://doi.org/10.1371/journal.pone.0037691).

The threat regime can change the exchange. Long exclusion of large herbivores altered host rewards and ant associations in an African savanna experiment. A mature reef on a quiet route may invest differently from a rapidly growing colony in a grazer-rich region. [Palmer et al. 2008](https://doi.org/10.1126/science.1151579).

A juvenile leaving a canopy nursery can carry associates, but transport is not complete inheritance. Workers without a reproductive individual or appropriate brood may not found a lasting colony. Microbial inocula need suitable habitat; aquatic residents need water through departure. Budding and module fragmentation must either transfer a reproductively sufficient partnership or encounter a reliable recruitment route.

Clustered attacks deserve separate tests. Observations of several predator types feeding simultaneously on a
guarded coral suggest possible overload as other food disappears, but do not supply a measured saturation
threshold. Many kilograms of defenders can still be poorly placed for concurrent damage at several modules.
[Kayal et al. 2018](https://doi.org/10.12688/f1000research.13118.2).

## What would distinguish workable communities

1. **Specify threats and damage.** Separate renewable grazing from dangerous punctures, attachment damage and reproductive loss. Measure attack frequency, damage before withdrawal and vulnerable organs.
2. **Measure each partner's contribution.** Compare host performance with and without a guild, and mixtures against their components. Record costs, predation among partners and effects on useful visitors.
3. **Close local response.** Measure detection probability, activation, travel, recruitment and deterrence under wet surfaces, wind and darkness. Nest spacing and guard density should follow those requirements.
4. **Close the shared food and nutrient budget.** Include brood, non-foraging individuals, replacement, external food, all guilds' reward shares, nutrient availability and export. Test the limiting interval as well as the annual total.
5. **Protect habitat through disturbances.** Couple pool survival, oxygen, wet mass and nutrient retention to rain, drought, ballast release, module loss and subsequent recolonisation.
6. **Follow a complete host life cycle.** Test establishment, partner turnover, reproduction and dispersal before treating a useful short-term association as a persistent reef community.

The present screen identifies requirements and tradeoffs. It leaves protection effectiveness, community persistence, empirical organism traits and the best mixture of species open.
