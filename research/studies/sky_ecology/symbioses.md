# Guards and reciprocal partnerships on sky reefs

**A mosaic of partners is the most promising defence for a sky reef.** Small patrols near vulnerable tissue, larger
hunters in sheltered spaces, mobile defenders around nesting sites, recycling communities in protected wet pockets and
externally feeding visitors at peripheral roosts each answer a different threat or loss. The author's direction lets the
guards be several species working together
([decisions](../../decisions.md#sky-reefs-partners-and-control-10-october-2026)). Ants deter elephants from occupied
acacia branches when other forage is at hand (Goheen & Palmer 2010), and coral guard-crabs of different sizes and
species protect against different predators (McKeon & Moore 2014). These are candidate functions for organisms still to
be chosen, and Earth experiments support some of the relationships they combine.

**A partner earns its place by the damage it prevents and the resources it returns, net of its own costs.** A kilogram
of occupants is a lifting load; it becomes a beneficial partnership when its contribution offsets its food, water,
shade, attachment and maintenance costs. In an ant-exclusion study, removing some partners improved host growth or
reproduction despite their protection (Stanton & Palmer 2011). Sea-feeding flyers can bring phosphorus in across the
reef boundary, pool predators can keep nutrients from leaving, and nitrogen fixers add nitrogen and no phosphorus.

**Guard mass, body size and placement are separate choices.** At 10 g/m² of live mass the allowance holds ten 1-g
guards per square metre, one 1-kg guard per 100 m² or one 100-kg guard per hectare, and real defenders would gather
near nests, growing tissue, wounds and reproductive structures.

**Guards must nest close to what they protect.** With a 30-second activation delay and a path factor of 1.5, arrival
within five minutes needs nests no more than about 5.1 m apart at 0.02 m/s, or 25.5 m at 0.1 m/s, and no spacing
meets a ten-second deadline. A giant colony therefore needs distributed small nests, mobile defenders between modules
and passive tissue defences wherever damage can come before the guards.

**All guards share the host's one tenant-food budget, and their phosphorus must come from food other than the host's
sugar.** At four turnovers a year, 10 g of live guards per m² need 25 g C/m² of food a year, a quarter of the
0.1 kg C/m² the host gives all its tenants, and 100 g need 250% of it. In the mixed example carbon closes and the
small patrol still lacks 15.1 mg of assimilable P/m² a year, because the sugar that fills its carbon carries no
phosphorus.

**Partners change with the reef's age and its threats.** Sequential ant partners gave acacias a greater predicted
lifetime benefit than any single species (Palmer et al. 2010), and long exclusion of large herbivores changed the hosts'
rewards and their ants (Palmer et al. 2008). A juvenile inherits a partnership by carrying partners able to reproduce or
by meeting a reliable recruitment route.

The [tenant screen](tenants.md) supplies the geometry and the conditional budgets. [symbioses.py](symbioses.py) keeps
collective live mass, individual mass, response distance and food allocation separate. Its
[scenarios](symbioses_scenarios.json) take the host's tenant allowance, the tissue carbon and phosphorus and the raft
geometry from the earlier screens; every other input, from guard densities, sizes, speeds and nest spacings to the
activation delay, path factor, dry fraction, turnover, efficiencies and the mixed guilds' shares and supplies, is a
design guess chosen to show how the requirements scale, labelled here once. The [product](results/symbioses.json)
records units, assumptions and provenance, and the [source register](symbioses_sources.json) records for each paper
whether its full text, its abstract or indexed passages were read. Numbers are *screen* values from symbioses.py unless
marked *derived*.

## What the guard mass means

The 1–10 g/m² guard entry in the [tenants note](tenants.md) is a design guess of collective mass per square metre of
reef footprint, wet or dry unspecified. The scenarios here use **live mass** at 1, 10 and 100 g/m².

At 10 g/m², all the mass in one size gives ten 1-g individuals per square metre, one 1-kg individual per 100 m² or one
100-kg individual per hectare; mixtures share the same allowance. A circular 2-km reef holds about 3.14, 31.4 or
314 tonnes at the three densities (*derived*). These figures convert mass into bodies and tonnes. How much guard mass
defence needs, and how much a reef can sustain, are open.

In the tenant screen's raft geometry the curved exterior is about 3.63 times the footprint, so 10 g per footprint m²
spread evenly gives about 2.76 g per m² of surface (*derived*). Defenders would concentrate near nests, growing tissue,
wounds and reproductive structures. Brood, resting animals and reproductive individuals carry part of the mass, so only
part of the total biomass patrols.

Individual size, collective density and placement therefore need independent choices. At a fixed mass allowance larger
guards are fewer. They may travel farther or deter larger intruders, and each of those gains needs a mechanism beyond
mass.

## Threats and complementary defenders

| Threat or loss | Candidate partner and service | What must work |
|---|---|---|
| Eggs, small grazers and organisms entering folds | Small social patrols remove or discourage them near sensitive host tissue | Reliable detection, access to crevices, recruitment and food that supports reproduction |
| Larger crawling grazers and concealed larvae | Resident hunters occupy junctions and sheltered margins | Access to prey without puncturing the host, and coexistence with small patrols |
| Flying herbivores and egg-layers | Flying predators hunt arrivals; social flyers could recruit around disturbances | Interception under the reef's wind, affordable flight and discrimination among visitors |
| Large browsers or animals stripping tissue | Swarming deterrents and territorial nesting flyers make feeding costly | Vulnerable attacker anatomy, timely response and an attacker willing to withdraw |
| Fouling of collectors or exposed exchange surfaces | Surface-maintenance animals remove obstructing material | Selective cleaning, with waste delivered to a catchment on the reef |
| Nutrient export in emerging insects | Resident pool predators retain part of the prey production | Their own residence, oxygen demand, diet and eventual nutrient return |
| Microbial disease and tissue decay | Candidate microbial competitors or selective sanitation | A measured disease benefit, beyond cleaning or the mere presence of microbes |

**Small defenders can deter very large browsers.** Goheen and Palmer's ant removals and feeding trials linked ant
abundance to elephant damage on acacias, and elephants avoided ant-occupied branches when alternative forage was
available; attacks on sensitive tissue were a proposed explanation. How such guards would fare against an insulated
mouth, a rapid bite or a determined, hungry sky grazer remains to be tested.
[Goheen and Palmer 2010](https://doi.org/10.1016/j.cub.2010.08.015).

**Defenders of different sizes cover different threats.** In coral experiments small guard-crabs defended against small
predatory snails, and intermediate and large crabs against different sea stars; species as well as size affected
performance. A mixture can cover the gaps one kind leaves; protection by an assembled mixture against every predator
remains to be shown. [McKeon and Moore 2014](https://doi.org/10.7717/peerj.574).

For a reef, ant-like patrols and flying hunters together make a stronger start than a single universal guard. Larger
territorial nesters are a further possibility: defending a nest can incidentally protect the host, locally and
seasonally unless another mechanism keeps the benefit going. Feeding elsewhere could support their larger bodies; their
landing loads, roost damage and nutrient delivery each need their own budget.

## Response distance and the first damaging bite

A first bite can come before defenders arrive. A useful requirement is

`detection delay + travel time + deterrence time < time to unacceptable damage`.

The right-hand side depends on the threatened function. Removing renewable green tissue may be tolerable; opening a gas
compartment, destroying a reproductive organ or severing a connection has a different consequence. Safe bite sizes and
attack durations differ from function to function, and each remains to be set.

On a square grid of nests spaced `s` metres apart, the farthest point in a cell lies `s / sqrt(2)` from its nest. The
screen calculates `delay + path factor × distance / speed` for a staffed nest connected to the target, with the
path factor standing for indirect travel.

| Nest spacing | At 0.02 m/s | At 0.1 m/s | At 2 m/s |
|---|---:|---:|---:|
| 5 m | 4.92 min | 1.38 min | 32.7 s |
| 20 m | 18.18 min | 4.04 min | 40.6 s |
| 100 m | 88.89 min | 18.18 min | 83.0 s |

Every entry uses a 30-second activation delay and a path factor of 1.5. Speeds and nest spacings are independent
scenario inputs, set apart from body size and from any Earth organism's traits. On an aerial route wind-relative and
ground-relative speeds differ; for a crawler, folds, gaps and wet surfaces could lengthen or block the route. The times
end at arrival, and successful detection, combat and deterrence time come on top.

Under these assumptions, arrival within five minutes needs square-grid nest spacing no greater than about 5.1 m at
0.02 m/s, or 25.5 m at 0.1 m/s. A ten-second deadline is impossible after the 30-second activation, even at zero
travel distance. These are geometric requirements for occupied nests; how many nests a reef keeps occupied is a separate
question.

Earth response times are local measurements. An ant-plant study measured responses to disturbed nest chambers in
seconds, with response probability varying sharply between settings, and a social-wasp experiment measured successful
prey capture in minutes inside a small screenhouse. A kilometre-scale reef needs its own fitted response constant.
[Youngsteadt et al. 2024](https://doi.org/10.1002/ecy.4449),
[Southon et al. 2019](https://doi.org/10.1098/rspb.2019.1676).

Distributed small nests therefore suit even a giant colony, since one central nest serves a kilometre of skin too
slowly. Mobile defenders could connect neighbouring modules while local patrols protect organs that cannot wait.
Wherever damage can come before recruitment, passive tissue resistance, wound closure and compartmentalisation stay
necessary candidates. These guard response times are separate from the colony's pressure-trim and landing response
times.

## Food quality and the shared reward budget

The reference host gives its tenants 0.1 kg C per footprint m² per year. That is one allowance, shared by ants, flying
defenders, cleaners and pool fauna together. Food from intercepted visitors or foraging elsewhere is recorded at the
reef boundary; predation on residents transfers food already in the local web.

For each guild the screen calculates annual replacement production from live biomass, dry fraction, tissue composition
and turnover. Ingested carbon is production divided by assimilation and population production efficiency, and that
efficiency already includes the respiration it implies, so a second complete maintenance allowance would count it twice.
A large warm-bodied flyer's energy expenditure needs its own estimate in place of these efficiencies.

The sensitivity below fixes dry mass at 30% of live mass, carbon at 50% of dry mass, four annual biomass turnovers,
carbon assimilation at 80% and population production efficiency at 30%. Phosphorus is 0.795% of dry mass with 65%
assimilation, midpoints of the earlier study's insect-based ranges. The guards' diet and sizes are left open.

| Live guard mass per footprint m² | Food C required per m² per year | Share of entire host tenant-food C budget | Dietary P required per m² per year |
|---|---:|---:|---:|
| 1 g | 2.5 g | 2.5% | 0.0147 g |
| 10 g | 25 g | 25% | 0.1468 g |
| 100 g | 250 g | 250% | 1.4677 g |

At 100 g live/m² the guards need more than the host's entire tenant-food carbon allowance; at one annual turnover they
would need 62.5% of it. Holding production efficiency fixed while changing turnover also changes the implied
respiration, so the comparison leaves a slowly reproducing animal's maintenance open. Turnover and efficiency each need
independent calibration for each guild.

The mixed example divides 10 g live/m² among a small patrol, mobile defenders and larger nesters. It allocates 24% of
the one host-food budget, uses 14.25%, and stipulates external food containing another 7 g C/m²/year. Carbon
requirements close to numerical precision, and the small patrol still lacks 15.1 mg of assimilable P/m²/year for its
prescribed replacement. Whether the sky supplies that external food is open.

External food carries carbon and phosphorus paired at the chosen composition, and only the phosphorus in the consumed
fraction receives credit. Host sugar fills the remaining carbon requirement within each guild's allocation. Unused
supplies and assimilation losses stay explicit. Feeding solely to obtain phosphorus, and disposing of the surplus
carbon, lie outside the screen. The guild allocations replace portions of the tenant community within the host's
existing mass and food account.

Phosphorus demand is kept separately. Sugar supplies energy and carries no phosphorus, so new animal tissue needs food
bodies, prey, microbial foods or another nutrient-bearing diet. Guano delivered to a reef can pass to guards through
decomposers, plants or prey; how much of it each guard can eat and assimilate directly needs measuring.

The calculation is a steady replacement screen. An initial colonisation cost, a daily ration through darkness, the other
essential nutrients and a stable age structure lie outside it. The food quality of real ant rewards matters: in a
Malaysian ant-plant experiment nutrient supply limited food-body production, colony size correlated with host food
production, and more workers were associated with less leaf damage.
[Heil et al. 2001](https://doi.org/10.1007/s004420000534).

A mature reef might afford a larger reward budget, and that budget competes with expansion, repair, replacement gas,
sailing and extra housing. Protection should be judged by the damage it avoids and by host performance at its cost: a
dense guard colony could defend well and still cost its host more than the damage it prevents. An ant-exclusion study
reported that reducing some partners improved host growth or reproduction despite their protective service, and rare
catastrophic damage may change the benefit over longer periods. This warning is qualitative and rests on the study's
accessible abstract. [Stanton and Palmer 2011](https://doi.org/10.1890/10-1239.1).

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

Ant-occupied *Dischidia* leaves provide a precedent for chambers near absorptive host tissues. Isotope estimates
attributed part of the leaves' carbon and nitrogen to ant-associated respiration and debris. Those fractions belong to
that epiphyte, so reef uptake rates need their own estimate, and recapturing respired carbon recycles carbon already
fixed. [Treseder et al. 1995](https://doi.org/10.1038/375137a0).

Predators can fertilise hosts: a bromeliad-associated spider supplied nitrogen and improved plant growth in experiments.
Whether that fertilisation imports nutrients depends on where the prey fed. The focal bromeliad was terrestrial, and the
proposed reef tank habitat remains to be tested.
[Romero et al. 2006](https://doi.org/10.1890/0012-9658(2006)87%5B803:BSIHPN%5D2.0.CO;2).

Pool predators can also change what leaves. In a bromeliad food-web experiment, detritivores increased host nitrogen
uptake when predators prevented nutrient-bearing prey from emerging and departing. The retention benefit depends on
residence time and emergence, so each predator's is examined on its own.
[Ngai and Srivastava 2006](https://doi.org/10.1126/science.1132598).

Nitrogen fixation and nutrient recycling do different jobs. Epiphytic fixation on pelagic *Sargassum* varies
substantially, with little or no dark fixation in the examined observations. It supports the possibility of diazotrophic
associates; a sky-reef rate needs its own measurement. Fixation adds reactive nitrogen from atmospheric N₂ and adds no
phosphorus. [Johnson et al. 2023](https://doi.org/10.1371/journal.pone.0289485).

Sea-feeding flyers could import phosphorus across the reef boundary. The import depends on external feeding, nutrient
intake, excretion at the roost, retention and subsequent availability, with export subtracted. The seabird literature
supports concentrated deposition at colonies, including harmful enrichment and disturbance. Its global phosphorus figure
is modelled deposition at breeding colonies; each flyer's and each reef's share follows from the factors above.
[Otero et al. 2018](https://doi.org/10.1038/s41467-017-02446-8).

## Water and maintenance can create conflicts

The proposed 1–5 L/m² of pools share the host's free-water store. Water kept as habitat serves as disposable ballast
only at the cost of the organisms and nutrients dropped with it. Persistent pool communities require a residual volume,
water quality and night oxygen supply. Humus and epiphytes may retain moisture and also add wet mass and their own
transpiration; the sign of their water benefit is open.

Cleaning should be selective and keep useful functions working. Ants maintained specialised capture surfaces in a
pitcher-plant experiment, restoring performance after contamination; some removed material was deposited outside the
pitcher. A reef would benefit from directing such material into recycling pockets where possible. This precedent covers
surface maintenance; antimicrobial protection and repair of a pressure-bearing membrane are untested.
[Thornham et al. 2012](https://doi.org/10.1111/j.1365-2435.2011.01937.x).

The [storm screen](../aerophytes/storms.md) already includes a hypothetical oxygen-scavenging lining, costing hydrogen
through `2 H₂ + O₂ → 2 H₂O`: 4–18% additional hydrogen, relative to the model's permeation loss. The stoichiometry
sets that cost; uptake kinetics, containment, night function and protection during a tear are open. The proposed
guanine-like diffusion barrier is a separate host-material mechanism, and a microbial film as its replacement is
untested.

A dormant host still carries other organisms with their own demands: active guards, decomposers, pool fauna and
gas-management microbes may continue consuming carbon or oxygen. Survival through the worst dark interval needs its own
budget beside the annual surplus. Defensive coverage could change with light, temperature and activity, even when the
year's food budget closes.

## Coexistence and inheritance

Specialised host modules are division of labour within the colony. Animals, microbes and epiphytes living on them are
other organisms, and their residence may be mutualistic, neutral, parasitic or dependent on circumstances. The label
guard describes demonstrated or proposed behaviour against a specified threat, and a predator earns it through that
behaviour.

Several species may coexist by using different chambers, prey, body-size ranges or activity periods. They may also eat
or displace one another. Spatially separated rewards could keep guards away from reproductive organs and useful
visitors; experiments on extrafloral nectar show that reward placement can alter ant occupation of flowers and
reproductive outcomes. [Villamil et al. 2019](https://doi.org/10.1111/1365-2745.13135).

Different partners may be useful at different host ages. Long-term acacia observations and demographic modelling found
that sequential associations could give a greater predicted lifetime benefit than any single ant species. This supports
succession across a reef's life, with partners following one another.
[Palmer et al. 2010](https://doi.org/10.1073/pnas.1006872107).

Tolerance and mutualism are distinct. Work on competing acacia ants distinguished host-monopolising mutualists from an
exploiting species whose unrelated workers could share a host. A reef-wide supercolony or a peaceful mixed-ant patrol
needs its own mechanism. [Kautz et al. 2012](https://doi.org/10.1371/journal.pone.0037691).

The threat regime can change the exchange. Long exclusion of large herbivores altered host rewards and ant associations
in an African savanna experiment. A mature reef on a quiet route may invest differently from a rapidly growing colony in
a grazer-rich region. [Palmer et al. 2008](https://doi.org/10.1126/science.1151579).

A juvenile leaving a canopy nursery can carry associates, and inheritance needs more than transport. Workers without a
reproductive individual or appropriate brood may found no lasting colony. Microbial inocula need suitable habitat;
aquatic residents need water through departure. Budding and module fragmentation must either transfer a reproductively
sufficient partnership or encounter a reliable recruitment route.

Clustered attacks deserve separate tests. Observations of several predator types feeding simultaneously on a guarded
coral suggest possible overload as other food disappears; a saturation threshold is unmeasured. Many kilograms of
defenders can still be poorly placed for concurrent damage at several modules.
[Kayal et al. 2018](https://doi.org/10.12688/f1000research.13118.2).

## What would distinguish workable communities

1. **Specify threats and damage.** Separate renewable grazing from dangerous punctures, attachment damage and
   reproductive loss. Measure attack frequency, damage before withdrawal and vulnerable organs.
2. **Measure each partner's contribution.** Compare host performance with and without a guild, and mixtures against
   their components. Record costs, predation among partners and effects on useful visitors.
3. **Close local response.** Measure detection probability, activation, travel, recruitment and deterrence under wet
   surfaces, wind and darkness. Nest spacing and guard density should follow those requirements.
4. **Close the shared food and nutrient budget.** Include brood, non-foraging individuals, replacement, external food,
   all guilds' reward shares, nutrient availability and export. Test the limiting interval as well as the annual total.
5. **Protect habitat through disturbances.** Couple pool survival, oxygen, wet mass and nutrient retention to rain,
   drought, ballast release, module loss and subsequent recolonisation.
6. **Follow a complete host life cycle.** Test establishment, partner turnover, reproduction and dispersal before
   treating a useful short-term association as a persistent reef community.

The screen sets out requirements and trade-offs. Protection effectiveness, community persistence, the organisms'
measured traits and the best mixture of species are open.
