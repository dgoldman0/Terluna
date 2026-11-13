# Engineering feasibility of the enclosed off-Earth programme

*Assessment date: 13 November 2026. Construction start: 2041. Evidence review, inherited requirements and bounded engineering calculations; no integrated design qualification.*

## Assessment

A staged industrial start around 2041 is a credible research and development objective. The architecture draws on established physical processes: solar electricity, chemical propulsion, excavation, thermal and electrochemical processing, pressure vessels, regenerative life support and industrial recycling. Fifteen years permits substantial development, including several generations of equipment and early operations in space. It does not establish that the particular production rates, lifetimes, autonomous maintenance and deployment sequence in the financial model will be available together.

The strongest qualification is at the beginning. The model's first fully occupied year has 12,500 residents in year 8, or 2049. Before that opening, the construction account already requires an average **37.4 million tonnes of asteroid-source feed per year for eight years**. A dated plan for installing and supplying that industry remains open. The same first occupied year requires **547 MW of firm lunar electricity** and **314 MW of source processing**, supported by much larger installed solar capacity. These are substantial industrial requirements even though the population is small.

At the other end of the schedule, approximately 6.47 million residents in year 200 are supported by an economy with extensive remote production and services. Reaching that scale need not require a new physical principle. It does require a sustained manufacturing, transport and maintenance industry whose expansion follows demonstrated performance. Continued technological development after 2041 is available to that industry. The calculations below preserve the current coefficients to expose their demands; they do not assume that 2026 equipment remains unchanged for two centuries.

The present assessment therefore supports **continued development of the staged architecture**, with a conditional 2041 commitment to production construction. It does **not yet establish the 2049 occupancy date or the published financial returns as engineering-qualified outcomes**. The most consequential unresolved conditions are the initial source/depot deployment, delivered power at actual sites, pressure-boundary construction and repair, reliable transport, and long-duration human habitation. Recycling and refurbishment improve the physical position materially, but their recovery and service-life assumptions need their own demonstrations.

## Boundary, method and evidence

This study covers the entire planned off-Earth economic zone served by **U**, through the enclosed buildout and continued operation at its population cap. U retains its established meaning as a general economic unit, potentially a currency or basket of U-related currencies. The engineering boundary includes lunar communities, mines and factories, near-Earth asteroid facilities, freight vehicles, depots, orbital infrastructure and the facilities supporting commercial services. Earth supplies equipment and expertise across the boundary. Open Moon, the habitat array and outer-system development are outside it. The orbital power stations considered here belong to the enclosed programme; they do not introduce the later habitat array.

The calculation reads the staged, demand-sized cohort product at provisioning revision `a8af8761`, including improved recovery and refurbishment. It preserves all 81 primary cases, six poor-source cases and six faster-rollout cases. Their existing capacity results are **81/81, 6/6 and 4/6 passing**, respectively. Both faster-rollout diversified cases encounter the resident-workforce limit. Those checks establish consistency with the producer's assumed capacities. They do not verify the hardware, source trajectories, operating reliability or construction sequence.

[The input manifest](feasibility_inputs.json) pins the physical product, engineering coefficients, constants, enclosure geometry and recycling assumptions by SHA-256. [The source register](feasibility_sources.json) records primary literature and programme reports, the sections consulted and the limits of each claim. It also identifies the resource-distribution and preservation studies on the existing `domain/resources` branch at revision `1d652941`; those records inform this assessment without replacing their producer or declaring sources released for extraction. The report is intended for engineering and economic integration: its recommendations are proposals, while inherited architecture and protection decisions retain their existing status.

Three kinds of uncertainty need different treatment. Physical laws constrain quantities such as pressure force and heat rejection. Engineering assumptions describe recoveries, availability, component life and operating procedures that can be tested. Programme uncertainty includes whether facilities, trained staff, source access and transport capacity arrive in the necessary order. The nine-portfolio ranges below are scenario ranges, **not probability intervals**. This study assigns no numerical probability that the entire programme succeeds; treating many correlated, unqualified subsystems as independent success percentages would produce unsupported precision.

## What must be built, and when

The following example is the strong-market propellant/logistics portfolio with 1% annual real market growth and combined recovery/refurbishment. It is useful because transport and source requirements remain visible. It is not a forecast that this portfolio will dominate. GW denotes a billion watts; Mt denotes a million tonnes. All flows in this table are annual operating requirements, not cumulative stocks.

| Requirement | Year 8 | Year 32 | Year 100 | Year 200 | Year 220 |
|---|---:|---:|---:|---:|---:|
| Residents | 12,500 | 15,763 | 222,240 | 6,468,289 | 6,468,289 |
| Firm lunar electricity, GW | 0.547 | 0.889 | 14.68 | 292.15 | 301.08 |
| Source processing electricity, GW | 0.314 | 1.044 | 19.62 | 245.38 | 249.06 |
| Asteroid-source feed, Mt/year | 3.59 | 11.55 | 204.81 | 2,681.43 | 2,722.03 |
| Cargo from sources, Mt/year | 0.119 | 0.382 | 6.76 | 88.35 | 89.65 |
| External mass delivered to low Earth orbit, Mt/year | 0.00152 | 0.0140 | 0.325 | 3.525 | 3.675 |
| Fresh lunar feed, Mt/year | 0.0346 | 0.174 | 6.82 | 2.00 | 2.15 |

Population is only one demand driver. Factories, exports, replacement equipment, transport propellant and remote power account for much of the material and electricity. Thus the services portfolios retain substantial industrial requirements, while differing industrial export strategies add different loads.

Year 200 also hides some construction peaks. In this reference case, fresh lunar feed reaches **62.33 Mt/year in year 168**, versus 2.00 Mt/year when construction has tapered and recovered material supplies much of renewal. External deliveries peak at **3.883 Mt/year in year 181**; source feed peaks at **2.725 Gt/year in year 190**. Ground electricity continues to increase slightly through year 220 as the economic assumptions continue operating beyond the population cap. A plant sized only to the year-200 lunar mining rate would miss the buildout requirement.

Across each portfolio's nine primary market/growth combinations, terminal requirements are:

| Portfolio | Lunar firm electricity, GW | Source processing, GW | Source feed, Gt/year | Mass to low Earth orbit, Mt/year |
|---|---:|---:|---:|---:|
| Raw-resource prices | 287–292 | 214–214 | 2.33–2.33 | 3.45–3.46 |
| Commodity refined materials | 287–292 | 215–325 | 2.34–2.52 | 3.45–4.00 |
| Specialty products | 288–349 | 217–591 | 2.35–5.22 | 3.53–13.79 |
| Premium lunar materials | 288–389 | 217–523 | 2.36–5.78 | 3.55–14.47 |
| Propellant and logistics | 287–292 | 214–245 | 2.33–2.68 | 3.45–3.52 |
| Sports, entertainment and other services | 287–292 | 214–214 | 2.33–2.33 | 3.45–3.46 |
| Computing and digital exports | 287–292 | 215–367 | 2.33–2.50 | 3.45–4.26 |
| Diversified production and services | 288–352 | 217–624 | 2.36–5.24 | 3.54–12.76 |
| Earlier manufacturing mix | 287–325 | 216–485 | 2.34–3.81 | 3.49–8.32 |

These loads include the producer's shared population programme. Similar values for raw resources and services reflect that shared infrastructure and the current export coefficients; they do not demonstrate equivalent profitability or detailed facility needs. In particular, an assumed expenditure per visitor, venue or computing facility is less engineering detail than a process flow sheet, equipment list and availability requirement.

The six poorer-source cases retain their modeled capacity passes by providing additional equipment. Their year-200 source-feed requirements rise to **6.69–15.04 Gt/year** across the services, propellant and diversified comparisons. Passing that capacity calculation therefore does not mean poor feed is inexpensive or equally easy to deploy. The resulting industrial scale, power and source-access requirements still need to be supplied.

### The first productive facilities

The eight-year seed contains **298.94 Mt of source feed** and **7.70 Mt of lunar raw material**. Dividing the source total by eight gives the 37.37 Mt/year construction average, over ten times the first occupied year's operating demand. At the inherited source-module rate of 25 tonnes/hour, 65% availability and the existing 25% capacity margin, that average corresponds to approximately **328 source mining-module equivalents**. The first occupied year's feed corresponds to about 31.5 modules. Neither number is a dated deployment schedule: some modules must be built, transported, powered and commissioned before they can supply later work.

The seed procurement is charged once. Its physical account totals 2.628 trillion constant USD2041, including a once-only common programme allowance of 155.8 billion USD2041. The broader financial initial-funding result includes other financing and reserve provisions; it is a different quantity. This study adds no second founding charge and does not treat the smaller operating-year requirement as the cost of establishing the system.

An executable initial sequence must identify, for each delivery, the departure depot and its fuel, available launch and assembly capacity, the destination's power and handling capability, commissioned production, crew accommodation, spares and the next departure window. Imported startup equipment and consumables can break dependencies between local manufacture and the facilities needed to manufacture it. Their masses and schedules must be explicit. The existing priced initial Earth-route fuel stock covers a specified first fleet fraction; the producer itself lists source/depot deployment, additional filling and commissioning losses, passenger access and preproduction accommodation as unresolved requirements.

A useful first commercial stage can be smaller than the full residential opening. Paid research, communications, remotely operated experiments and limited material products can accompany industrial demonstrations. Larger entertainment audiences on Earth do not require proportionate visitor transport, which is an important engineering advantage. Nevertheless, live sports leagues require trained participants, suitable venues and reliable life support; hospitality requires additional occupied facilities and passenger services. Each commercial opening needs its own minimum viable facility and operations plan. No change to the adopted population schedule is made here.

## Power delivery, reserves and heat

The model already includes a substantial constellation-capacity allowance. With ten stations and two assumed to serve the load, transmitter conversion, capture and receiver efficiencies of 0.8 each, and an end-of-life factor of 0.85, installed solar DC capacity is

\[
P_{DC}=P_{ground}\frac{10/2}{0.8\,0.8\,0.8\,0.85}=11.49P_{ground}.
\]

The reference opening consequently implies **6.29 GW of installed DC generation** for the lunar beam supply, and the year-200 case implies **3.36 TW**. At the assumed system-level 100 W/kg these represent approximately **62,900 tonnes** and **33.6 Mt** of power hardware. Source processing has its own solar supply: approximately 0.661 GW DC at opening and 516 GW at year 200, before separately checking depot systems. These ratios are already present in the equipment account; applying another fivefold visibility penalty would duplicate an allowance.

The fixed station geometry also needs examination as output grows. At the assumed 32% photovoltaic conversion and normal illumination near 1 AU, the beam supply implies approximately **14.4 km² of ideal collector area at opening** and **7,707 km² at year 200**, before layouts and deployment clearances. If ten stations retain one 1.2 km aperture each and two deliver the load, each transmitting station emits approximately **428 MW of radio-frequency power at opening** and **228 GW at year 200**. Averaged over the aperture, that rises from approximately 378 W/m² to **202 kW/m²**. The scalar watts-per-kilogram coefficient does not qualify that concentrated power architecture.

At the terminal output, 80% radio-frequency conversion produces about **57.1 GW of conversion heat per active station**. The same ideal 350 K radiator assumptions used below imply about **93 km² of emitting area per station** just for that heat. Converter radiators are already intended to fit within the inherited complete-system mass allowance; this screen checks that obligation rather than adding a second hardware charge. Photovoltaic waste heat is rejected separately by the collectors. More stations, multiple transmitters, larger apertures or another distribution arrangement are credible avenues for later development; their geometry, mass and operating costs need explicit accounting. This is a reason to redesign the terminal power network as it grows, rather than assume an unchanged station geometry is established by a capacity ledger. The adopted baseline remains unchanged in this study.

Power availability still needs a constellation and ground-grid solution. [Caltech's SSPD-1 report](https://www.caltech.edu/about/news/space-solar-power-project-ends-first-in-space-mission-with-successes-and-lessons) provides evidence of wireless transfer in space and transmission toward Earth, with subsequent degradation and thermal/electrical lessons. It does not qualify kilometre apertures, hundreds of megawatts per service link, continuous coverage or the model's mass and lifetime. A complete design must address aperture deployment, phase control, scan angle, sidelobes, receiver geometry, handover, outages and repair.

### Receiver size is a real constraint

For an ideal uniformly illuminated circular aperture focused at the receiver, encircled power is

\[
E(x)=1-J_0(x)^2-J_1(x)^2,\qquad x=\frac{\pi D r}{\lambda R}.
\]

At 2.45 GHz, a 1,200 m aperture and 10,000 km path, 80% capture requires a receiver radius of **915 m**, or **2.63 km²**. At 12,000 km the required area is 3.78 km². This is an ideal focusing calculation. The conventional far-field distance is about 23,500 km, so using an unfocused far-field antenna claim at the shorter distance would be inappropriate; the calculation assumes the required focusing phase profile. Real aperture illumination, scan and pointing must be evaluated separately.

The inherited 200 W/m² is **delivered electric flux**, after receiver conversion. At that flux, 2.63 km² supplies about **526 MW**. Sizing a 5 MW receiver solely from the same flux gives only 0.0250 km² and approximately **1.9% ideal capture**. For a 50 MW site the corresponding capture is about 17%. Small sites therefore cannot simultaneously retain the assumed flux-sized receiver area and 80% capture on that link. A large receiver serving several communities, larger receivers operating at lower mean flux, a changed transmitter/link design, or other local supplies are options requiring design and costing. The result does not show that every beam or the entire constellation fails.

A circular orbit at the adopted altitude has a period of about **31.7 hours**. A central lunar eclipse under a point-Sun approximation lasts about **1.50 hours** for one station. Neither quantity demonstrates simultaneous coverage by two stations or the longest correlated outage. Solar geometry, Earth eclipses, maintenance and the actual sites require an ephemeris-based power schedule.

### Reserves need a service requirement

The inherited design already provides a **72-hour critical reserve at 3.5 kW per resident**, including 0.5 kW for food preservation. Its nominal inventory is 105 GWh per 250,000 residents. The original sizing uses an 80% usable fraction, 90% discharge efficiency and a further 20% reserve allowance. Retaining all three definitions leaves 60% of nominal inventory for dispatch while preserving that allowance, and reproduces the 72-hour target exactly.

The same energy would supply the much larger **firm design load for 5.8 hours at opening or 5.6 hours at year 200**. If essential services instead require 25% of firm demand, autonomy falls to about 23 and 22 hours; at a 10% fraction it is about 58 and 56 hours. These sensitivities do not show an error in the existing 72-hour reserve. They identify the need to verify that the specified critical load covers the services actually required during an outage. The ratios assume charged inventory and adequate discharge power. Ageing, protected wiring, fuel-cell capacity and repair access still need equipment-level treatment.

Normal beamed supply does not require storing an entire lunar night. A prolonged common failure of that supply is a different operating condition. Essential services must be defined and supplied through the recovery time of a credible failure, including a failed receiver or damaged grid. Fission, local solar, stored reactants and geographically separated feeds can contribute, but their individual failures and dependencies must be considered. The [DOE's January 2026 programme description](https://www.energy.gov/ne/articles/5-things-you-need-know-about-fission-surface-power-systems) discusses a 40 kW lunar demonstration objective; it is evidence of development, not a qualified gigawatt substitute. Fusion is unnecessary to the assessed baseline.

### Heat rejection remains industrial infrastructure

Year-200 source waste heat in the reference case is **195.3 GW**. With emissivity 0.9, the existing 25% margin and an unobstructed cold sky,

\[
A=\frac{1.25Q}{\epsilon\sigma T^4}.
\]

| Radiator temperature | Ideal emitting area | Area if effective rejection is 70% of ideal |
|---|---:|---:|
| 300 K | 591 km² | 844 km² |
| 350 K | 319 km² | 455 km² |
| 400 K | 187 km² | 267 km² |
| 500 K | 76.6 km² | 109 km² |

The 70% column is an engineering sensitivity, not a measured degradation rate. At 350 K and 5 kg/m², the ideal surface allowance is approximately **1.59 Mt**. Emitting area is not necessarily land footprint or one-sided panel plan area. Coolant loops, pumps, deployment, micrometeoroid isolation, shadows, view factors and surface contamination need design. Raising the temperature requires an appropriate process or a heat pump; a low-temperature load cannot reject heat at 500 K for free.

Source power also depends on the selected orbit. The fixed 70% solar-flux assumption corresponds to about 1.20 AU under inverse-square illumination. At 1.5 AU the collector area for the same electrical service is 1.58 times that reference; at 2 AU it is 2.8 times. A representative flux is not a complete eccentric-orbit power supply. Source selection, production duty and cargo windows must be solved together.

## Sources, processing and manufacture

The newer resource work is important to this assessment. It integrates catalogue entries with missing sizes and compositions, spectroscopy and a debiased missing-population model. Its wider near-Earth orbital selection (semimajor axis at most 1.8 AU, eccentricity at most 0.5 and inclination at most 30°) contains a conditional **25–73 Gt of water equivalent**, with a median of 39 Gt, before recovery and restrictions. The broader near-Earth population contains **504–2,931 Gt**, with a median of 961 Gt. These overlapping pools cannot be added. Their existence supports a much larger material opportunity than the previously characterized candidates alone; it does not place the material in Earth orbit or qualify repeat cargo routes.

The resource study's lunar comparison also finds that bulk availability is large relative to the enclosed programme. Its prospective 1,000 km², five-metre layer contains approximately 7.84–9.21 Gt under the stated density assumptions. That is a regional inventory inference, not a surveyed mine: composition measurements near the surface do not verify grade, excavation properties or continuity throughout five metres. Accessible geology, preservation exclusions, haul distance and product chemistry determine which fraction can be used and at what cost.

For the reference year-200 water output of **187.7 Mt/year**, a recovered water yield of 7% requires approximately 2.68 Gt/year of feed. At 3% it requires 6.26 Gt/year; at 12%, 1.56 Gt/year. This holds water output fixed and exposes a feed sensitivity. It does not rerun the coupled energy, nitrogen, transport and replacement account. Bound hydroxyl must be released from minerals; water-equivalent abundance cannot automatically be treated as pumpable ice. In the current coefficients, 0.5 kWh/kg feed alone becomes about 7.1 kWh per kg of recovered water at 7% yield, before downstream processing.

The inherited source train handles 25 tonnes/hour with 65% availability and 250 tonnes of equipment. Year-200 reference feed implies approximately **23,530 train equivalents**, including the existing margin. A lower 40% availability would require 1.625 times as many equivalents at fixed feed. Combining that availability with a reduction in recovered yield from 7% to 3% increases the equivalent count by approximately 3.79 times at fixed water output. These are useful equipment sensitivities, not a cost forecast; shared power, logistics and staffing would need a fresh coupled calculation.

The near-term demonstrations support development without yet supporting those rates. [NASA's reported molten-regolith experiment](https://ntrs.nasa.gov/citations/20260006167) used 25 kg of simulant, averaging 0.07 kg/hour of oxygen during nine hours of sustained electrolysis under vacuum. [NASA's CaRD prototype](https://www.nasa.gov/centers-and-facilities/johnson/sunlight-extracts-oxygen-from-regolith-using-solar-chemistry/) integrates solar concentration, a reactor and controls, with carbon monoxide production confirmed. [Blue Origin's September 2025 announcement](https://www.blueorigin.com/news/blue-alchemist-hits-major-milestone-toward-permanent-sustainable-lunar-infrastructure) reports a critical design review and a subsequent demonstration objective. These are distinct achievements and evidence states. They do not establish continuous production of structural alloys or precision equipment.

The manufacturing qualification must follow products through separation, impurity removal, alloy additions, casting or forming, joining, inspection and service. Generic iron, aluminium or silicon abundance does not deliver a pressure-qualified laminate, cryogenic tank, bearing, seal, semiconductor or optical assembly. The model retains imported components, including the nonlocal share of source tugs; this is a viable division of production if the import supply, failure rates and replacement schedule are secured. Full industrial autarky is unnecessary, while unpriced dependence on critical imports would overstate resilience.

## Freight, depots and traffic safety

Near-Earth asteroids are heliocentric bodies. The model's 2.5 km/s one-way source velocity increment and three-year round trip are route assumptions. An orbital-element filter is insufficient to establish repeated departures, adequate delivered mass, affordable return propellant or a safe approach corridor. At least one source-to-depot route needs a fully dated mass, energy and navigation solution before its output supports the opening schedule. Several independently accessible sources would reduce interruption risk; counting correlated candidates as independent redundancy would not.

The assumed source payload is **10,000 tonnes per departure**. Reference year-200 cargo of 88.35 Mt requires approximately **8,835 departures per year**, or 24.2 per day. Three-year cycles, 80% availability and the existing margin imply about **41,414 vehicle equivalents**, close to the producer's surviving stock of 41,467. A 20-year mean life permits only several such cycles, so tanker inspection and refurbishment affect both procurement and availability. The opening stock of 592 tugs includes construction requirements; it must not be confused with the first occupied year's 11.9 departures or its much smaller steady throughput requirement.

Little's law gives **44–133 Mt of outbound cargo in transit** for the terminal flow and illustrative one-way travel times of 0.5–1.5 years. This is cargo only; vehicles and propellant are additional. The calculation identifies physical inventory and interruption exposure. It does not add another purchase cost to cargo already financed in the cash-flow model. A quarter-year depot stock also does not by itself cover a source outage longer than a year; the network needs inventories matched to disruptions and mission windows.

[NASA's LOXSAT preparation report](https://www.nasa.gov/directorates/stmd/tech-demo-missions-program/cryogenic-fluid-management-cfm/nasa-industry-prepare-cryogenic-fuel-technology-demo/) describes a nine-month liquid-oxygen mission intended to test eleven fluid-management technologies, including boiloff control, transfer, pressure and gauging. It is planning evidence in the reviewed source. Liquid-hydrogen storage, high-throughput transfer, contamination control and repeated fault recovery remain separate qualification requirements for this architecture. Installed refrigeration power does not prove that plumbing, tank insulation and docking seals achieve the assumed losses.

External delivery to low Earth orbit at year 200 is equivalent to approximately **35,248 launches/year carrying 100 tonnes each** in the reference case. At 200 launches per pad per year this is 176 pad-equivalents before outages, downrange constraints, refurbishment and environmental limits. These are century-scale requirements, not demands for 2041. Conversely, the opening year's approximately 15 annual operating launches exclude the much larger seed deployment. Neither number is a price quotation or a demonstrated launch-market capacity. The inherited 200 USD2026/kg is an economic scenario input.

Passenger access requires its own schedule and reliability case. Cargo throughput does not establish a way to transport 12,500 residents, rotating workers or visitors with adequate life support, radiation exposure management, emergency return capacity and medical care. Lunar landers also need to demonstrate the modeled 30 cycles/year and ten-year mean life under dusty conditions. Cargo and passenger vehicles can share facilities where engineered to do so, but their qualification evidence is different.

### High cadence changes the consequence of small failure probabilities

For independent cargo-loss probability \(p\) and \(N\) departures, the chance of at least one loss is \(1-(1-p)^N\). At 8,835 departures/year it is approximately 58.7% for \(p=10^{-4}\), 8.46% for \(10^{-5}\), and 0.88% for \(10^{-6}\). These are illustrations of cargo loss, **not Earth-impact probabilities or approved risk limits**. An illustrative 1% annual loss limit would require a per-departure probability below about 1.14 × 10⁻⁶. Establishing that level from zero failures alone at 95% confidence would require about 2.63 million independent trials. A few successful demonstrations cannot carry the full reliability claim.

The practical response is fault containment, safe failure geometry, conservative operating limits and evidence about failure mechanisms, combined with operational experience. The standing six traffic requirements remain: passive failure, Earth-clear catcher geometry, remote braking, packet-energy limits, tracking/interception/corridors, and verified capture before planetary intercept. A 10,000-tonne payload has 31.3 TJ of kinetic energy at 2.5 km/s and 125 TJ at 5 km/s. Those illustrative energies are not arrival-speed predictions or an adopted packet-energy cap. They show why payload size needs an explicit safety justification and why a navigation or software common cause cannot be represented solely by independent component reliability.

## Habitation and biological services

### Pressure structures and sites

The city comparison uses independent 60 m diameter pressure halls in larger cavities, with 500 m isolation lengths and hemispherical ends. At one atmosphere, the cylindrical hoop load is **3.04 MN per metre**, and the end force is **286 MN**. An assumed 600 MPa effective allowable gives a pressure-only membrane equivalent of 5.07 mm, before the existing detail factor and separate liners, supports and connections. This is a force balance, not a manufactured wall specification. Biaxial material behaviour, joints, penetrations, sustained loading, inspection, fire and repair must be established at the relevant scale.

The accepted small-town arrangement remains one contiguous dome. For a representative 300 m span and rise equal to 20% of span, pressure uplift is approximately **7.16 GN**. A three-metre shield at the assumed 1,600 kg/m³ density offsets only **8.93%** of that uplift. Below about **9.04 kPa** internal pressure difference, the shield's weight exceeds pressure support in the whole-roof balance. Roof and foundation systems must therefore support both normal inflation and an unpressurized condition. Regolith shielding does not remove the anchorage problem.

City cavities require geophysical surveys, access surveys, rock-mechanics investigations and construction trials. Roof cover provides shielding; it does not certify an airtight or stable pressure boundary. The large assumed cavity dimensions and connected usable lengths are site requirements, not verified properties of a selected location. Excavated or surface alternatives can preserve industrial development if a proposed cavity fails qualification, but their changed construction, shielding and access costs must be calculated. Internal refuges and isolated utility equipment can improve fault response while preserving the accepted contiguous town atmosphere.

Radiation, dust and fire need integrated designs. Shield thickness by itself is not a calculated lifetime dose: composition, penetrations, secondary radiation, shelter use, occupational exposure and the adopted population-health criteria matter. [NASA's lunar electrodynamic dust-shield demonstration](https://www.nasa.gov/image-article/nasas-dust-shield-successfully-repels-lunar-regolith-on-moon/) supports active cleaning of glass and radiator surfaces. It does not establish mine-equipment life or dust exclusion from all habitats, suits, seals and bearings. Clean/dirty interfaces, maintainable filtration, fire detection, smoke control and escape/refuge capacity belong in the enclosure design before occupancy.

### Recovery is a set of coupled loops

The current improvements raise domestic and industrial water recovery from 98% to 99.5%, crop-water recovery from 99.5% to 99.7%, and fish-purge recovery from 95% to 98%. For the same collected stream, moving from 98% to 99.5% reduces makeup by a factor of four. [The ISS's demonstrated 98% recovery](https://www.nasa.gov/missions/station/nasa-achieves-water-recovery-milestone-on-international-space-station/) is a relevant precedent within its collection and treatment boundary. It is not evidence that every household, industrial, agricultural and landscape flow is captured at that rate.

Qualification requires flow meters and inventories around each loop, water-quality monitoring, contaminant removal, membrane and filter replacement, brine and solids handling, and failure bypass arrangements. Carbon, oxygen, nitrogen, phosphorus and other nutrient accounts must connect people, crops, food processing, waste treatment and stores. The inherited 5 kg/person/day landscape/agriculture allowance remains in the primary model until its possible overlap with explicit crop and fish losses is resolved. Removing it wholesale would manufacture an engineering saving before the streams have been reconciled.

NASA's [2022 life-support baseline, Table 3-25](https://ntrs.nasa.gov/api/citations/20210024855/downloads/BVAD_2.15.22-final.pdf), gives 0.90 kg oxygen/day for its reference long-mission crewmember. Used only as a scale analogue, that is approximately **4,109 tonnes/year at opening** and **2.13 Mt/year at the cap**. These are gross metabolic flows, not required oxygen imports: plants and physicochemical processes can regenerate oxygen. The calculation demonstrates why an atmospheric leakage allowance cannot stand in for a complete metabolic and food-system balance. Population composition and activity require a separate design basis.

Industrial recovery is equally important. The current usable fractions are approximately 95.1% for collected/accepted goods and 98.0% for structural material under their assumed collection and quality factors. Recovered mass must also meet a usable specification. Alloy mixing, coatings, embedded electronics, refractory wear, chemical contaminants and repeated recycling can constrain reuse. Chemical rocket exhaust remains expended in the current architecture; a high local recycling rate cannot reclaim a plume after it has left the recoverable system.

### Long habitation needs direct evidence

The reviewed [NASA-hosted partial-gravity assessment](https://ntrs.nasa.gov/citations/20210019591) identifies the lack of direct long-duration human evidence and reasons for concern from indirect evidence. Its proposed health interpretation is not a demonstrated universal gravity threshold. It does make clear that experience with selected adults on short missions cannot establish lifelong safety for a general population.

During the lead period, habitat trials can establish many operational conditions: hygiene, trace contaminants, nutrition, food production, work patterns, emergency medicine and maintainability. Lunar or other relevant partial-gravity studies must address the remaining biological questions. Fifteen years cannot supply multigenerational human evidence. Residency duration, countermeasures, medical monitoring and return capability therefore remain design conditions, with any artificial-gravity medical facility requiring separate engineering and cost treatment. A whole habitat array is not assumed as a solution within this boundary.

## Automation, maintenance and service exports

Autonomous excavation is only part of an autonomous industry. Ore handling, jams, tool changes, abnormal product, leak repair, sensor replacement, cleaning and recovery after communications loss determine practical staffing. The source model assumes 65% mining availability; lunar processing assumes 80%. Those assumptions must be supported by maintainable modules, spares, local diagnostics and a realistic division between Earth-based operators and resident intervention. More capable software in fifteen years is a reasonable development opportunity, while a wage or asset-value staffing proxy remains an unverified productivity assumption.

The reference resident-workforce ratio peaks at **0.358 in year 28**, where one means the producer's allowed workforce has been fully used. At the cap it is about 0.083. These are model ratios, not fractions of the entire population. The faster diversified cases exceed that modeled allowance in year 24. Early staff and commissioning can therefore constrain a growth path even when the terminal economy appears to have ample labour. Task-level staffing, shift coverage, training and medical contingency are the next evidence needed.

Refurbishment is a consequential and testable improvement. Collection of 80% and repair acceptance of 60% give a 48% eligible-unit return. Rebuilds use new parts equal to 20% of the unit's mass, and the current producer restores the original lifetime. For a stationary, fully used population of eligible units, with the same acceptance rule at each retirement, the new material purchased per service-year relative to all-new equipment is

\[
F=\frac{1-q+qs}{1-q+q\alpha},
\]

where \(q=0.48\), \(s=0.20\), and \(\alpha\) is restored life divided by new-unit life.

| Restored service life | New material per service-year | Material saving before repair operations |
|---|---:|---:|
| 100% | 61.6% of all-new reference | 38.4% |
| 75% | 70.0% | 30.0% |
| 50% | 81.1% | 18.9% |
| 25% | 96.3% | 3.8% |

This sensitivity concerns eligible equipment. It is not a reduction of the whole operating budget. Workshops, inspection, consumables, energy, transport and downtime are additional. Restoring more than 20% of new life is the material-only break-even under these assumptions; economic break-even depends on those additional costs. A fifteen-year programme cannot observe every proposed 20–50-year equipment life. Accelerated tests must reproduce the relevant failure mechanisms, and operating designs need condition monitoring and life limits rather than indefinite automatic renewal credits.

Sports, entertainment, research and other services can earn revenue from worldwide customers while using relatively little material per additional customer. That advantage is compatible with substantial engineering requirements for the initial venues, production systems, communications and staff. Multiple sports and visiting competitors increase scheduling, facilities and passenger-service needs without making every Earth viewer a traveller. Computing adds continuous electrical and heat-rejection demand, radiation-sensitive hardware, replacement imports, communications capacity and latency constraints. Its competitive value must be tested against contemporaneous Earth facilities; vacuum does not eliminate cooling equipment. Premium lunar products similarly need demonstrated product advantages and repeatable quality, not merely lunar provenance.

## Preservation and source qualification

Preservation affects engineering before extraction begins. The resource study's protocol register distinguishes adopted decisions, proposals used as conservative study conditions and open questions. This assessment retains that distinction. Surveys, landings, exhaust, drilling and access construction can disturb a record before a mine opens, so the information and physical-protection work must precede those activities.

| Requirement carried into engineering | Evidence and provision needed |
|---|---|
| Pre-disturbance survey, representative samples and records | Survey coverage tied to the work area; documented sampling and characterization before irreversible alteration; protection of scientifically important settings |
| Physical, functional and heritage preservation | Site-specific assessment; records do not substitute automatically for a retained place or an ongoing scientific function |
| Native-condition samples and an unopened fraction | Appropriate vacuum, inert handling or temperature control; storage equipment and chain of custody; no invented universal retained percentage |
| Multiple locations and institutions | Independent physical stores and records, continuity arrangements, access and sustained maintenance |
| Scientific provenance and access | Traceable observations, calibration, uncertainty, versions, scientific review and usable archives; file hashes alone do not provide this programme |
| Future resource access and scarce-element recovery | Extraction plans that preserve future options and account for nutrients, residues, recovered equipment and material left in place |
| Dust, exhaust and polar records | Contamination and transport assessment for the actual activities; no commodity credit for polar ice before required scientific treatment |
| Scientific holds and traffic safety | Bennu, Ryugu, Itokawa and Eros retained consistently across source pools; all six cargo requirements applied to actual routes |

Buffers, released fractions, exact archive capacities and operating obligations remain site-specific. An undisturbed archive cannot preserve every function of a destroyed site. No listed candidate becomes authorized through this feasibility calculation, and a generic preservation allowance cannot substitute for a surveyed, staffed and funded programme. The later Open Moon atmosphere-loss budget is not imported as a numerical permission to contaminate the present environment.

## Development during 2026–2041

The lead time affects different technologies differently. Applying the same assumed improvement rate to every subsystem would obscure the distinction between improving an existing process and demonstrating a new operating system. The following judgments concern development paths; they are not numerical forecasts or formal technology-readiness ratings.

| Capability | Judgment with fifteen years of development | What still determines feasibility |
|---|---|---|
| Excavation, thermal processing, separations and routine fabrication | Credible candidates for productive modules by 2041, using established processes adapted through repeated trials | Actual feed properties, useful product quality, vacuum/dust performance, tool life and energy per saleable kilogram |
| Water treatment, industrial recovery and repairable equipment | Credible opportunities for material improvements over the inherited coefficients | Collection boundaries, contaminant accumulation, rejected products, repair labour and restored service life |
| Remote operations, inspection and robotic maintenance | Substantial improvement is a reasonable planning opportunity; productivity needs measurement task by task | Recovery from unplanned faults, communications loss, common software defects and the resident work that cannot be delegated remotely |
| Cryogenic depots, reusable freight and large power links | The interval permits an integrated demonstration programme and redesign; qualification of the assumed scale is still a critical condition | Delivered mass and power, operating losses, repeated transfers, trajectory windows, availability and maintainability |
| Large occupied pressure enclosures | Structural loads are calculable and modular construction is plausible; the selected structure and site remain major engineering work | Qualified joints and materials, foundation/cavity behaviour, depressurized support, penetrations, inspection and emergency response |
| Continuous food production and long habitation | Ground and lunar trials can resolve many operating questions before opening | Coupled biological stability, complete nutrient/gas accounts and human partial-gravity effects; lifelong evidence cannot be compressed into a short programme |
| The year-200 industrial fleet | Primarily a later manufacturing and operations challenge if the first systems perform; subsequent technology can reduce its burden | Demand, finance, useful source access, workforce and the ability to reproduce and maintain qualified capacity |

The following overlapping programme is a proposed way to use the fifteen-year allowance. Its dates are targets. Each stage must produce evidence before the dependent commitment; elapsed time by itself does not satisfy a gate.

| Period | Main work | Evidence required before the next dependent commitment |
|---|---|---|
| 2026–2029 | Define interface requirements, choose candidate sites and sources, begin preservation surveys and sample/record infrastructure; reconcile food/air/water and industrial ledgers | Traceable site and source shortlist, protected-area constraints, dated candidate routes, critical-import list and testable system requirements |
| 2028–2033 | Continuous ground processing with representative mixtures; integrated cryogenic, dust, robotics and repair rigs; pressure joints and full structural segments | Measured yield, product quality, energy, consumables, downtime and interventions across expected variation; pressure and fire failure modes; verified repair inspection methods |
| 2030–2036 | Lunar excavation/manufacture and enclosure demonstrators; orbital power and depot trials; source reconnaissance and production pilot | Operation through environmental cycles, source recovery under actual conditions, safely contained residues, reliable transfer and fault recovery; preservation records completed before disturbance |
| 2033–2038 | First complete source-to-depot cargo cycle, repeat departures, integrated power delivery and sustained habitat services | Closed mass/energy accounts over the mission cycle; dated route and contingency supply; delivered electrical availability at real sites; maintainable life-support and food systems |
| 2037–2041 | Qualify production modules and suppliers, demonstrate repair and spares logistics, complete the eight-year deployment plan and occupancy criteria | Measured capacity and availability sufficient for each initial production step; pressure-boundary acceptance; actual accommodation/passenger plan; no dependency on an uncommissioned supplier |

A three-year route cycle must begin early enough to complete, diagnose problems and repeat relevant operations before the production commitment. Ground rigs should run in parallel with flight preparation where their dependencies permit it. Source and habitat demonstrations need enough duration to expose wear, contamination and maintenance burdens. Proposed testing should include repeated continuous campaigns and at least a complete seasonal/operating cycle where meaningful, with durations justified by each failure mechanism; a single generic hours threshold would not qualify every subsystem.

The inherited common programme allowance equals **100 billion USD2026**, or about **6.67 billion USD2026/year** if spread evenly across fifteen years. That arithmetic shows the scale of an existing allowance. It does not establish an adequate research budget, and it must not be charged again as new expenditure. A work breakdown must allocate laboratories, flight demonstrations, survey campaigns, archive facilities, qualified suppliers, staff and contingency, then identify which costs are already included in seed procurement and which are genuinely additional. Moving expenditure into the pre-2041 development period also requires moving its financing dates in the next financial calculation; this assessment makes no silent change to that cash-flow timing.

By 2041, the programme needs qualified modules and an executable first production sequence. It does not need 23,530 source trains, 41,000 tugs or tens of thousands of annual Earth launches already operating. Subsequent expansion should be released in increments supported by measured source performance, maintenance burden, infrastructure availability and market demand. Standard interfaces, replaceable processing trains and reusable construction equipment make later technological improvements useful without requiring each site to repeat the original development programme.

## Consequences for the next integrated model

The engineering findings warrant additional scenarios rather than an arbitrary percentage reduction in every return. The current financial results should remain conditional on the physical coefficients until the following quantities can be propagated through their actual dependencies:

| Engineering uncertainty | Change to propagate | Principal economic effect |
|---|---|---|
| Initial deployment and source qualification | Dated commissioning, initial fuel/stocks, temporary imports and accommodation | Initial funding, first sales, financing duration and working inventory |
| Real sites and power links | Receiver footprints, distribution, scanning losses, independent reserves and outage duration | Capital, land use, firm production and service interruption |
| Grade, recovery and routes | Joint source assay, energy use, orbital power duty, flight time and delivered mass | Feed, power, fleet size, inventory and supply concentration |
| Product quality and qualified local manufacture | Usable alloy/product yields and critical component imports | Export volume, margins and replacement exposure |
| Repair and maintenance | Restored life, rejection rate, interventions and actual spare consumption | Replacement purchases, labour, availability and reserve stocks |
| Human habitation and commercial facilities | Qualified occupancy, passenger access, venues, health services and food closure | Growth timing, services capacity and operating obligations |

An efficient next calculation would compare the existing target coefficients with measured-pilot coefficients when available, plus explicit engineering stresses: lower source availability, poorer recovered yield, shorter repaired life, longer route cycles, larger power receivers and delayed initial commissioning. Improvements should also be admitted when demonstrated, including lower process energy, better recovery and cheaper maintainable hardware. Their savings must flow through the coupled producer; the report's single-system ratios cannot be applied directly to the whole cash-flow account.

The useful engineering prospect is a programme that learns while producing: imported first equipment establishes qualified local capacity, recovery reduces recurring freight, repaired modules extend the value of earlier investment, and worldwide services help pay for further development. The remaining work is to demonstrate that sequence at the first operational scale and make later increments depend on the evidence it produces. That makes the fifteen-year lead period a concrete engineering resource while keeping the long-run opportunity visible.

## Reproduction and review

Run from the provisioning checkout:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/studies/enclosed_feasibility/feasibility.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q research/studies/enclosed_feasibility/test_feasibility.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python visualization/enclosed_feasibility/plot.py
```

The [machine-readable product](results/feasibility.json) retains all 93 cases' checkpoints, operating peaks and inherited admission state, together with beam, storage, thermal, cargo, pressure and repair sensitivities. It records producer and input hashes. Generated PNG/PDF figures and rendering provenance go to ignored `research/runs/enclosed_feasibility/`; the renderer reads the study product and introduces no physical model.

[The check record](checks.md) distinguishes numerical controls, source reading and substantive review. This is a bounded feasibility study. It includes no new detailed trajectory integration, finite-element enclosure qualification, radiation-transport calculation, complete ecological simulation, hardware endurance trial or source-release decision.
