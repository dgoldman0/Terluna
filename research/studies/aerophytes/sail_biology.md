# Living sails and distributed control on sky reefs

A sky reef can be conceived as one photosynthetic colonial organism with excitable tissues, contractile organs and specialised growth zones. The author's direction on 10 October 2026 permits that combination without requiring the organism to fit an Earth plant or animal body plan. A central brain is optional; local sensing, coordination and memory still need physical tissues, energy and developmental instructions.

The strongest working proposal is a **persistent, renewably maintained support system carrying replaceable sail organs**. Each sail has local sensing and control, while the colony sends slower instructions about deployment and preferred direction. A substantial part of the sail and cable could be nonliving extracellular material maintained by a much smaller living component. The alternatives below retain seasonal sails, extensive retraction and independent partners where their costs or advantages differ.

This is biological ideation with explicit requirements. The [sailing screen](sailing.md) calculates a steady force balance after prescribing wing area and lift coefficients; it does not establish a controllable wing. The separate [requirements calculation](sail_biology.py) and [product](results/sail_biology.json) examine scale, signal delay, transport and replacement under stated scenarios. They do not simulate a nervous system, flight stability, a viable life cycle or a complete resource budget. Source access and limits are recorded in [sail_biology_sources.json](sail_biology_sources.json).

## What the sail actually does

The upper reef and lower sail occupy different currents. Their relative motion provides airflow past the sail; a shaped surface at an appropriate angle generates a sideways force, transmitted through the tether. The sail's aerodynamic lift is approximately horizontal. The reef's buoyancy carries the suspended mass. Equal currents remove the energy source for this particular mechanism.

A local steering tail could set the main sail's angle to the apparent wind. Its action needs an aerodynamic moment opposed by the main surface and its suspension, with a stable resulting orientation. Adjustable bridles, changing membrane curvature and a movable trailing flap are alternatives or supplements. They are different mechanisms to compare, not interchangeable decorations. Pitch, roll, yaw, cable swing and interactions between neighbouring sails must all remain bounded.

A flexible trunk chiefly transmits tension. Turning the reef cannot be assumed to turn a freely swivelling sail kilometres below it. The engineered StratoSail concept places its rudder actuator at the wing; its author draft describes a ribbed flexible skin, a boom and local controls. Its limited model tests and simulated routes provide a mechanical precedent, not a biological design at reef scale. [Aaron, Heun and Nock](https://gaerospace.com/projects/ULBDStratoSail/pdfs_docs/COSPAR2000StratoSail.pdf).

The earlier 0.4–0.8 m/s across the shear therefore remains a conditional result. A real living surface could have less lift, more drag, a smaller safe angle range and considerable time unavailable for steering. The six sails in the conversation sketch are a visual arrangement, not a selected number or a revised aerodynamic calculation.

## An anatomy with distinct jobs

| Part | Proposed organisation | Main obligation |
|---|---|---|
| Attachment region in the reef | Load paths spread through several float modules, with separate sealed service connections | Carry pull without peeling a single module or opening a gas chamber; tolerate loss of one attachment |
| Main suspension | Parallel structural fibre bundles with a protective sheath and smaller living service strands | Carry sustained weight and gust loads; permit staged fibre replacement under load |
| Branching hub | Short branches to several separated sails, with strain sensing at their roots | Limit unequal loads, collisions and twisting; keep a damaged branch from dragging down the others |
| Leading rib and membrane bays | Reinforced edges, short supported panels and replaceable extracellular film | Maintain an airfoil shape when wet, folded and repeatedly loaded |
| Steering organ | Local opposing contractile tissues acting on a tail, flap or adjustable bridle | Supply the required force, torque, stroke, speed and fatigue life |
| Local control tissue | Mechanosensory cells, an excitable network and a small energy reserve near each actuator | Stabilise the sail and reduce load without awaiting a distant command |
| Renewal zones | Protected growth tissue at panel margins and structural junctions | Replace specific parts while preserving load paths and service continuity |
| Lower dense tissue | Reserves or structural thickening placed where its weight helps orientation | Supply useful ballast without uncontrolled changes in centre of mass |

A membrane need not be densely cellular everywhere. Thin repair tissue could secrete and maintain a tougher extracellular sheet, while veins carry tension and service small patches. This follows the existing [carbon allocation](biology.md): inert material incurs construction and renewal costs, while active tissue additionally respires. A hollow or pressurised rib is another candidate, but maintaining its pressure adds seals, fluid and leakage costs. It cannot be counted as weightless stiffness.

The main cable should have distinguishable structural and service functions. A dead fibre bundle cannot conduct living electrical signals or transport nutrients by itself. Conversely, a thick wet nerve or vascular trunk should not be assigned the strength and dry mass of an ideal fibre. Protective layers, living sleeves and reserve stores belong in the suspended mass and renewal accounts.

## Control without a central brain

Earth plants provide precedents for long-distance electrical and calcium signalling, signal-triggered movement and resource integration between connected modules. Those processes do not establish plant neurons or a guidance system. Colonial animals provide a complementary precedent for specialised motor organs and signalling among them. A proposed reef may use excitable epithelia, nerve-like cords, local ganglia-like knots, chemical signals and mechanical feedback in different places.

Wound-triggered glutamate/calcium signalling in plants and induced defence among connected clover ramets support the idea of differentiated responses across a connected colony. The clover network's source–sink relations also caution against assuming that every module communicates equally with every other module. [Toyota et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30213912/), [Gómez et al. 2008](https://doi.org/10.1111/j.1469-8137.2008.02542.x).

The siphonophore *Nanomia bijuga* is especially useful as an organisational analogue: its swimming organs combine nerve rings, colony-connected tracts, excitable epithelia and muscle, and some local responses persist after the main tracts are interrupted. *Hydra* imaging also distinguishes functional neuronal networks associated with different behaviours. Neither supplies a reef controller, but both support examining more organisation than a single undifferentiated alarm network. [Norekian and Meech 2020](https://doi.org/10.1242/jeb.233494), [Dupre and Yuste 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5423359/).

Three levels of control suit the physical scale:

1. **Local protection.** Strain, hinge movement and apparent airflow affect the nearby actuator or a passive release. An overloaded panel reduces its exposed area or changes angle. Passive geometry should reduce load during a signal failure, with the actual response measured rather than assumed.
2. **Sail coordination.** Neighbouring organs share tension and angle information. They change tack gradually, avoid sweeping through one another and report sustained trouble. A heavily loaded branch overrides a request for more steering.
3. **Colony behaviour.** Slower signals express a preferred side, a deployment limit or a need for water. Changes in light, hydration, temperature, reserves and repeated load history alter that preference. This layer need not specify every hinge angle.

This architecture avoids a single decision point whose damage disables the whole rig. It introduces coordination problems: neighbouring controllers can oppose one another, delayed corrections can drive oscillations, and repair can reconnect a signal path with the wrong gain or timing. Refractory intervals, dead bands and limits on the rate of adjustment are candidate biological functions; their values must come from stability analysis and measured tissue response.

A useful local controller senses **apparent wind**, not a map of global weather. Hairs or flexible sensory projections could measure bending; paired pressure-sensitive patches could distinguish sides of the flow; strain cells could report tendon load; hydration and chemical signals could report resource state. Venus-flytrap experiments show that mechanical reception and threshold-dependent electrical/calcium propagation can be distinguished, with the MSL10 channel affecting weak-stimulus sensitivity. That supports examining sensory thresholds, not assigning an airflow sensitivity to the reef. [Suda et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485109/). Their sensitivity, fouling and mechanical interaction with the airfoil remain requirements. A local light compass could provide a reference direction, but the relationship between light, time and position changes along a drifting route.

The colony could test small sustained changes and retain a chemical or neural-like memory of whether hydration or light income improved. That is an exploratory control proposal. It does not establish learning, consciousness, precise geographic navigation, advance storm prediction or knowledge of where the nearest coast lies. Nutrient-importing visitors might provide cues, but their presence is not a reliable coordinate system.

## Signal distance and actuator effort

The following are one-way transit times from `distance / speed`, not response times of a proposed organ. The 5-km path is a round-distance illustration; the inherited cable for a 5-km vertical drop is 5.56 km long and would take about 11% longer to traverse.

| Assumed propagation speed | One-metre path | Five-kilometre path | Status |
|---|---:|---:|---|
| 0.005 m/s | 200 s | 11.6 days | Order of the small-organ electrical/calcium rates reported in *Mimosa* |
| 0.3 m/s | 3.33 s | 4.63 hours | Measured order for *Nanomia* nectophore epithelium, not its stem nerve |
| 1 m/s | 1 s | 83.3 minutes | Design sensitivity |
| 10 m/s | 0.1 s | 8.33 minutes | Design sensitivity |
| 100 m/s | 0.01 s | 50 s | Design sensitivity |

*Mimosa* measurements were tissue- and stimulus-specific, about 4–6 mm/s in the reported rachilla experiments. The *Nanomia* epithelial measurement was approximately 0.3 m/s at 11.5–15 °C. Extending either rate over kilometres demonstrates a distance penalty, not that the corresponding tissues can survive or propagate a signal over that distance. [Hagihara et al. 2022](https://doi.org/10.1038/s41467-022-34106-x), [Norekian and Meech 2020](https://doi.org/10.1242/jeb.233494).

These times exclude sensing, synapses or relays, muscle response and the motion of the sail itself. Faster dedicated cords could improve communication, but no chosen conduction speed demonstrates that such tissue can be maintained along a wet, loaded cable. Splitting the route into relays protects signal integrity; it does not eliminate total propagation time.

Local control also needs a suitably small physical organ. A hundred-metre flap cannot be treated as a millimetre-scale sensitive-plant leaf enlarged without changes in inertia, stiffness or fluid demand. A distributed row of shorter flaps may shorten control paths, at the cost of more living tissue and coordination. Passive structural damping and automatic unloading remain valuable even with faster nerves.

There are several actuator candidates. Contractile tissue can make repeated reversible changes but needs oxygen, fuel and a way to hold position economically. Osmotic or hydraulic cells can set tension and curvature, with water availability and transport limiting response. Elastic structures with a latch can switch quickly after slow charging, while reset costs and fatigue limit repetition. Differential growth suits persistent trim and seasonal remodelling; it cannot be the sole actuator for every gust. Flytrap snap-buckling supplies a small-scale precedent for releasing stored elastic energy; it does not establish safe snap-through of a large sail or equally rapid resetting. [Forterre et al. 2005](https://doi.org/10.1038/nature03185). Combining these mechanisms is more credible than assigning rapid, strong, low-cost motion to one generic tissue.

Steering effort depends on **hinge torque and travel**, not simply on the entire tether force. A small control surface can influence a larger surface through the airflow, but that does not guarantee a small actuator at every angle. The requirements product illustrates how aerodynamic moment, moment arm and muscle attachment determine force and work. Hinge imbalance, friction, elastic preload, stall, dynamic overshoot and sustained muscle metabolism can all increase the cost. No actuator mass or operating efficiency has been selected.

For scale, a deliberately hypothetical 100-m² panel with a 0.5-m force arm and 0.25-m tendon arm, moved through 30° in ten seconds under constant resisting torque, gives:

| Assumed resultant loading | Hinge torque | Tendon force | Mechanical work per motion |
|---|---:|---:|---:|
| Saved mean sail loads, 0.22–0.58 N/m² | 11–29 N m | 44–116 N | 5.8–15.1 J |
| Selected stress case, 143 N/m² | 7,150 N m | 28,600 N | 3,744 J |

This is a torque requirement example, not a simulated flap or a demonstrated muscle. Uniform load, moment arms, area, angular travel and timing are choices; holding, inertia and biological inefficiency are excluded. The stress case makes clear why ordinary trim and overload protection need different mechanisms.

## Behaviour through a change in weather

| Situation | Candidate response | What must be established |
|---|---|---|
| Suitable shear and adequate reserves | Expose some panels, stabilise their angle, then increase area within load limits | Positive benefit after service costs; stable wing and tether dynamics |
| Desired change of side | Reduce force, change local trim, pass through a low-force state, then rebuild force | Safe transition without entanglement, shock loading or loss of height |
| Weak or reversing shear | Reassess angle and reduce unnecessary area; accept drift when steering is ineffective | Reliable flow sensing; no assumption that all winds permit the requested heading |
| Rising or oscillating load | Local feathering or folding, followed by branch isolation if needed | Unloading actually lowers force; reaction is fast enough under load |
| Forecast or sustained severe weather | Long-term shortening where affordable, protective stowage, or avoidance by route | Adequate warning and a demonstrated retrieval mechanism |
| Control or service failure | A passive low-load state, redundant service route, or sacrificial loss of a small part | The failed part does not remain a violently flapping drag surface |
| Long darkness | Retain essential reflexes and tissue maintenance; reduce discretionary collecting and steering | Local fuel and oxygen balance throughout the longest relevant dark interval |

Feathering means reducing aerodynamic force by orientation. Furling means reducing exposed membrane area. Retraction means bringing the suspended assembly upward. These operations have different mechanics and must not share an uncosted instantaneous response. Even turning edge-on leaves ribs, lines and folded material exposed to drag.

Releasing a sail or cable relieves the reef of weight as well as pull. That can trigger an unwanted ascent, recoil and a change in the load distribution among remaining modules. A breakaway design needs a corresponding buoyancy and water response; losing a heavy sail is not automatically a benign event. Releasing a branch may also export its nutrient stock and create a hazard below.

## Supplying an organ kilometres below its parent

The lower sail needs water, carbon, minerals, oxygen and repair materials. Local reserves decouple short disturbances from supply interruptions. Photosynthetic patches or symbionts could contribute carbon when the sail receives useful light, but the upper reef can shade it and sail orientation is chosen for force. More photosynthetic tissue also adds wet mass and maintenance. Absorbing light below the reef cannot be counted again as light already used above.

A continuous water-filled conduit across a 5-km vertical separation develops roughly **8 MPa of hydrostatic pressure difference** under lunar gravity. A downward supply can release gravitational potential energy; it still needs containment and flow regulation. Returning water upward requires work. In a closed circulation the gravitational gains and losses can offset in the ideal energy balance, while local pressure, friction and leakage remain. A hanging organism does not evade the physical costs of long-distance fluid transport.

A candidate arrangement is a series of pressure-isolated reservoirs and service nodes, with valves or pumps between them. True pressure breaks require chambers and controlled exchange; drawing valves along one uninterrupted filled tube does not reset its hydrostatic head. More stages trade smaller pressure differences for added walls, fluid mass, failure points and metabolic work. Separate downfeed and return paths would need resistance and flow budgets. Local collection and use of water reduce transfer demand, but cannot deliver nutrients back to the upper reef without a return path.

Experiments on connected *Alternanthera* ramets under heterogeneous water availability demonstrate resource integration and unequal allocation benefits over an eight-week study. They support connected resource sharing as a biological concept, not any assumed flow capacity or pressure tolerance over kilometres. [You et al. 2016](https://doi.org/10.1038/srep29767).

This makes a mostly inert cable with spaced service nodes and narrow connecting living strands more attractive than a uniformly fleshy five-kilometre stalk. Electrical biological signalling still needs a continuous hydrated excitable pathway between nodes; inert gaps require a different communication mechanism. Motor and control organs also need gas exchange while retaining water, so their protective sheaths cannot be credited as both impermeable against drying and freely oxygenating active tissue. A fully separate lower organism could feed itself more independently, but then reliable attachment, communication, reproduction and reciprocal benefit become a symbiosis problem. Neither architecture can be chosen from tensile strength alone.

## Growing a complete rig

A juvenile can begin as a drifting buoyant unit with short sensory or collecting appendages. It need not possess the mature reef's steering ability. A short cable samples little vertical shear, so the adult's speed cannot be assigned to it merely by preserving the ratio of sail area to reef area. If juvenile survival requires steering, that becomes an additional nursery or juvenile-organ requirement.

An initial construction sequence worth developing is:

1. **Establish flotation and reserves.** Grow enough productive and buoyant modules to carry the next appendage increment while meeting repair, tenant and dark-period obligations.
2. **Build an attachment nursery.** Develop several load paths and protected growth tissue on the underside. Keep control and service connections separate from gas containment.
3. **Deploy a small pilot organ.** Grow a short supported pendant with a modest folded sail and local sensors. Its early benefit may be sensing or collecting; useful steering is a later possibility to calculate.
4. **Extend in affordable increments.** Lengthen the structural bundle alongside its service strands, strengthening the upper cable as the suspended load increases. Each intermediate geometry needs its own load and carbon allowance.
5. **Bud blades near the lower hub.** Establish the leading rib, control organ, supply connection and membrane supports before exposing a large new area. Open sectors gradually as their tissues mature.
6. **Rotate renewal.** Grow replacements while neighbouring blades remain useful, recover accessible nutrients from retiring tissue, and release only the residual parts.
7. **Provision reproduction or fission.** Each daughter needs working growth zones, stores and enough independent flotation. Sharing one cable between daughters that separate is not a completed inheritance scheme.

There are several ways to lengthen the cable, and none follows automatically from calling it a stem:

| Growth mechanism | Advantage | Difficult part |
|---|---|---|
| Distal growth near the sail | No large store of folded cable at the reef | Supply reaches an increasingly remote bud; newly formed material must mature before carrying large tension |
| Growth at many intermediate nodes | Growth and repair can occur near damage | Soft expanding regions interrupt a highly stressed tension member unless bypass fibres carry the load |
| Proximal production and controlled payout | Construction occurs near the main carbon source | Alternating grips or another load-transfer organ must let new material enter service without losing support or circulation |
| Deployment of pre-grown folds | Opening can be faster than new biosynthesis | The stored cable was already paid for and carried; folds require space, protection and controlled release |
| Serial budding of service modules | Repeated units can contain reserves, valves and repair tissue | Junctions, sealing and coordinated load transfer become principal failure sites |

The working cable hypothesis is an overlapping fibre bundle with protected service nodes. New strands take load before old strands are retired; the load-bearing bundle therefore persists while its constituents change. This still needs an account of how fibres attach, mature and transfer tension. It is not evidence that soft growth tissue can safely elongate an already loaded kilometre-scale cable.

In *Nanomia*, spatially restricted proliferation and candidate stem-cell populations support ordered budding and specialisation, while cautioning against assuming regeneration everywhere in a mature colony. Retained renewal sites are therefore a specific proposed organ. A separate colonial tunicate experiment demonstrated whole-body regeneration from small vascular fragments; that offers an alternative biological precedent, not a growth rate for an aerial cable. [Siebert et al. 2015](https://doi.org/10.1186/s13227-015-0018-2), [Rinkevich et al. 2007](https://doi.org/10.1371/journal.pbio.0050071).

## Permanent structure and replaceable parts

Permanence can apply to a location, a function, or the same physical material. A reef can retain a steering hub for decades while repeatedly replacing its cells and fibres. It can also keep the costly cable while replacing only its exposed membranes. These distinctions matter more than choosing a single lifespan for the entire organ.

| Life history | Best reason to consider it | Main cost or risk |
|---|---|---|
| Persistent cable with renewable, foldable blades | Keeps access to useful shear and avoids repeatedly building the support system | Cable and service tissue remain exposed even when blades are furled |
| Continuously exposed, continuously renewed rig | Available steering and sensing with few deployment transitions | Chronic drag, wetting, predation and storm exposure; no protected interval |
| Seasonal blades on a persistent cable | Reduces active surface and maintenance during a long unfavourable interval | Regrowth and nutrient recovery take time; reduced steering may prolong the bad conditions |
| Fully retractable rig | Changes working depth and brings distal parts close for repair | Hauling work, storage, bending damage, concentrated load and a substantial retrieval organ |
| Whole rig shed and regrown | Can isolate a severe failure without a retrieval system | Repeated cable construction, nutrient export and a potentially long period without steering |
| Separate living sail partner | Distal metabolism and specialised motor biology could be more autonomous | Recruitment, compatibility, reciprocal benefit and reliable communication must evolve or be designed together |
| Detachable reproductive appendage | Planned retirement might also disperse a provisioned offspring | Independent flotation, reserves and development must be supplied; once detached it cannot steer its parent |

The preferred starting variant is the first: retain attachment and renewal tissue, renew the cable through overlapping load paths, and replace membrane bays more readily than sensory and motor organs. It preserves expensive depth access while allowing ordinary damage to remain local. This preference could change if persistent cable drag, supply costs or frequent entanglement dominate.

Seasonal shedding is more plausible when a long predictable interval offers little steering benefit. A simple rule to shed every night is poorly supported: wind can power sailing in darkness, and navigation may help the organism avoid still worse conditions. Fast release should be reserved for failure states, with planned retirement allowing more time to reclaim resources.

For the inherited 5.56-km cable, selected constant retrieval speeds imply:

| Retrieval speed | Time to gather the paid-out line |
|---|---:|
| 0.01 m/s | 6.43 days |
| 0.1 m/s | 15.4 hours |
| 1 m/s | 1.54 hours |

These are arithmetic scenarios, not proposed biological performance. At fixed reef altitude, raising only the inherited 62–214-tonne bottom weight through the 5-km vertical drop requires **0.50–1.72 GJ**. Raising the cable, drag, friction, losses, storage and the hauling organ all add costs. Folding a sail leaves its weight in place; shortening the whole cable is a separate operation.

Even the engineered StratoSail draft omitted powered reel-up from its then-current design to save mass. A biological hauling organ, long-term shortening through tissue remodelling, and passive payout from stored folds are distinct development choices. Local furling remains the candidate fast response; no whole-rig retrieval rate has been demonstrated.

## Construction costs and phosphorus recovery

The representative 2-km circular reef has 3.14 million m² of footprint. A sail at the assumed 10% share has 314,159 m² of area, equivalent in area to a 560-m square. In the parent product's 5-km-drop cases, the prescribed sail force is about 70–182 kN. Its mean force per sail area is only 0.22–0.58 N/m², while the selected tether and bottom weight together range from about 97 to 328 tonnes across the two drag and fibre-strength choices. These are ideal support cases, not reef-specific anatomy.

The bottom weight alone ranges from 62 to 214 tonnes. A finished sail averaging an assumed 0.2–0.5 kg/m² would weigh 63–157 tonnes. That overlap makes detailed sizing worthwhile. It does not prove that wet membrane, ribs, muscles, conduits and fittings fit the budget. Sail tissue can replace some required ballast only when its location and load path provide the needed weight and balance; a distributed membrane is not mechanically identical to the model's point weight. A change in shape, drag, total mass or centre of mass requires the hanging geometry and host lift account to be recalculated. Across this sail, 0.1 kg/m² of additional retained water adds 31.4 tonnes; 1 kg/m² adds 314 tonnes. A collecting fringe therefore needs drainage and a wet-load limit alongside its capture function.

The prior sizing also combines simplified cable drag with an idealised hanging shape rather than solving a cable in distributed wind and gravity together. Attachments to the reef, membrane stiffness and actuator torque are outside its result. The separate [structure screen](structure.md) gives a 143 Pa design gust at 5 km: multiplying that by the full sail area and an assumed force coefficient of one gives about 45 MN if the entire surface meets that pressure. This is a load warning, not a simulated sail response or an estimate that every gust loads every panel equally.

The number of sails has a cost. Dividing a fixed ideal tensile load among N identical circular cables of the same material and length keeps their summed cross-sectional area constant but increases their summed diameter, and hence projected cylindrical drag area, by sqrt(N). Real cables also have their own weight and unequal loads. A few shared trunks with shorter terminal branches are worth comparing with many independent long lines; neither arrangement is automatically optimal for a deformable reef. A single trunk concentrates failure, while many trunks increase drag and entanglement paths.

Each component needs its own dry mass, living fraction, water content, tissue phosphorus fraction and renewal schedule. Structural polymer can be relatively nutrient-poor while motors, growing buds and service tissue are more demanding. In a steady inventory, planned phosphorus export from component i is

`P_loss_i = dry_mass_i × P_fraction_i × (1 − recovered_fraction_i) / renewal_interval_i`.

The recovered fraction here means phosphorus retained inside the colony boundary and available for reuse. Moving it from a blade into a hub that is shed immediately afterward is not recovery at that boundary. Emergency loss needs a separate occurrence rate and recovery fraction; planned and emergency replacement must not both charge the same event. Carbon construction, nutrient stock and recurring export are distinct quantities.

The requirements product varies dry sail loading and turnover with a single explicit, hypothetical tissue composition. It does not yet resolve the separate living and structural tissues or emergency losses. At 0.1 kg dry matter/m² of sail and 50% carbon by dry mass, the 314,159-m² surface contains 15.7 tonnes of carbon. With an assumed phosphorus fraction of 0.001 and 80% recovery before the residual tissue leaves the reef, replacement gives:

| Replacement schedule | Carbon incorporated in new tissue | Phosphorus exported in unrecovered material |
|---|---:|---:|
| 10% annually | 1.57 t C/year | 0.63 kg P/year |
| Complete annual replacement | 15.7 t C/year | 6.28 kg P/year |
| Complete replacement every lunar month | 194 t C/year | 77.7 kg P/year |

These are material flows, not total synthesis or maintenance costs, external carbon demand, affordable renewal rates or measured recovery. The product also records 0.02 and 0.5 kg dry matter/m² sensitivities. Emergency shedding with little recovery would export more phosphorus per kilogram lost; a replacement that remains inside the reef's food web is not an external loss. Resorbed phosphorus can help build the next blade but does not remove its carbon construction, water, transport and maturation costs.

Petioles of *Oxalis pes-caprae* provide a mechanical release precedent: preferential breakage occurred at a basal weak zone, and replacement leaves came from surviving meristems. That supports designing both a break plane and surviving renewal tissue; it does not show an actively commanded joint. Nutrient resorption varies with species and stress, so emergency detachment should not be credited with perfect nutrient recovery. [Shtein et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6408339/), [Urbina et al. 2021](https://doi.org/10.1002/ece3.7734), [Xu et al. 2020](https://doi.org/10.1093/jpe/rtaa053).

Budget allocation must close at the whole reef. Extra sail tissue, innervation, collectors and repair draw from the same surplus used for tenants, gas replacement, reproduction and other structures. A sail that imports phosphorus but consumes too much carbon, water or structural reserve can still reduce host persistence. The [symbioses analysis](../sky_ecology/symbioses.md) supplies the same distinction between a useful service and net benefit.

## Additional functions and their conflicts

A useful sail can do more than steer, but some combinations work best in different tissues or at different times.

| Additional function | Candidate location or mode | Benefit to examine | Tradeoff that must remain in the account |
|---|---|---|---|
| Flow and strain sensing | Leading edges, suspension roots and a small exposed pilot surface | Sensing already serves control and can inform deployment | A point sensor may miss a gust elsewhere; fouling changes calibration |
| Cloud-water capture | Replaceable fringes outside the clean airfoil and small local stores | Supplies distal tissues and may reduce downfeed | Droplets add weight and roughness; transporting excess upward costs work |
| Particle or prey collection | Dedicated retractable combs on separate edges | Net imports require externally sourced material; recaptured reef litter retains existing nutrients | Capture rate, digestibility and phosphorus content are unknown; drag and fouling can erase the gain |
| Photosynthesis | Selected illuminated ribs or separate small blades | Supports local metabolism in favourable light | Reef shade, sail angle and dark intervals; living wet mass and nutrient demand |
| Heat exchange | Service nodes with controllable exposed area | May help a locally warm control organ | Long-range cooling of the parent needs circulation; unwanted night cooling and water loss |
| Local defence | Guard chambers near permanent hubs and protected tendons | Protects vulnerable joints and maintenance tissue | Occupants consume food, add mass and may foul moving parts |
| Reproductive dispersal | A specifically provisioned bud near a renewal zone | Shares construction infrastructure with organ growth | An expendable sheet is not automatically a viable propagule |
| Stable signalling surface | Pigment or orientation changes visible to neighbouring organs or partners | Could supplement slower internal communication | Occlusion, darkness and required sensory machinery; no assumed long-range signalling rate |

Cactus experiments demonstrate fog-droplet collection and directional movement over spines into absorptive tissue. They support patterned collectors on short scales, not extracting liquid water from unsaturated air or lifting it kilometres to the parent. [Ju et al. 2012](https://doi.org/10.1038/ncomms2253). In the reef, lower collectors should first support lower organs; export to the upper colony remains a separate transport calculation.

A reasonable division would keep the main aerodynamic membrane relatively smooth, put collecting hairs on replaceable fringes, and protect muscles and growing tissue near permanent junctions. Capturing water and prey could be favoured when the sail is feathered; high-force steering might favour clean surfaces and reduced collecting area. Mode changes must preserve the combined mass and nutrient ledger.

Large roosting animals and persistent pools belong on stable reef modules rather than on frequently moving or sacrificial blades. Small guards, decomposers and surface cleaners could occupy protected hubs, provided their maintenance does not block hinges or add uncontrolled wet mass. Cleaning that throws nutrient-bearing material into the air defeats retention unless that export is paid for. A hanging sail also offers predators an approach to sensitive tissue; defensive benefit cannot be assumed from its size.

## One organism or a partnership

For the first design, the float modules, cable, control knots and sail organs can be differentiated parts of one genotype. Their local autonomy is physiological, and their development can share inherited rules. A tissue containing genuine neurons would be proposed neural tissue, while a calcium-signalling epithelium would remain a different mechanism. There is no need to settle an Earth kingdom label before comparing their costs.

Photosynthetic symbionts and external guards can still be different species. A separate sail organism is a richer alternative: the reef offers attachment and carbohydrate while its partner offers a motor organ, local sensing and perhaps imported resources. Reliable inheritance or recruitment, signal interpretation, load sharing, recognition and conflict resolution become additional requirements. A partner's incentive to detach and reproduce may conflict with the host's need for steering during a storm.

Fusion of independently founded reef modules could improve survival and construction, but compatible circulation, signalling and reproductive interests cannot be inferred from contact. Even a single-genotype colony needs controls on pathological proliferation and local resource monopolies. Protected growth zones should be regulated by whole-colony supply and load capacity.

## Forms worth carrying forward

Three contrasting ecotypes would expose useful differences without choosing one universal reef:

- **Persistent open-air sailor.** A small number of long support trunks, sparse living service tissue and repeatedly renewed low-mass panels. It invests heavily in local stability and retains an ability to steer through darkness.
- **Cloud-edge collector.** More collecting fringes, local water storage and short-range nutrient processing, with fewer assumptions about returning water to the upper reef. Its additional wet load and drag may require less exposed sail or a smaller depth separation.
- **Seasonal or migratory sailor.** Expands blade area ahead of a sustained useful wind regime and resorbs it during a long low-benefit interval. It retains control and growth hubs, with reduced steering while the surface regrows. Forecasting can be an inherited seasonal programme, but its success must be tested against actual routes.

A low coastal form using a sea drogue is a different mechanical ecology: the second medium is water, so it may obtain useful force without kilometre-scale air shear. It adds wave loads, immersion, salt and terrain clearance rather than solving the high-altitude form's problems. Its separate [sailing comparison](sailing.md) should remain separate until an anatomy is specified.

The promising visual anatomy is therefore a reef with durable suspension roots and distinct organs at their ends: some narrow clean steering blades, some collecting fringes, folded panels resting beside new buds, and older panels being retired. The number, dimensions and tissue mix should follow the next mechanical comparison rather than the appearance of the sketch.

## Requirements for the next model

The next useful calculation is a single shaped sail with its suspension, tail, finite actuator and local delay in changing airflow. It should report angle error, loads, oscillation growth, control effort and behaviour after a loss of power. Matching a commanded sign of lift is insufficient. Vary pivot and centre-of-mass positions, tail area and arm, material damping, membrane deformation and the mode used to reduce area.

Then connect several organs through compliant cables and a deformable buoyant reef. Include line drag on the actual hanging geometry, unequal gusts, wave propagation along the cable, mutual wakes and release of a loaded branch. An acceptable local controller can still make a coupled system oscillate or twist itself together.

A biological budget follows alongside that mechanics: measured or explicit hypothetical dry and wet masses, turnover, actuator metabolism, fluid pressures and flow, construction, phosphorus recovery and night reserves. Trace every stage from juvenile through cable extension, mature renewal, failure and reproduction. Only after those constraints can coexist should active trajectories be compared with passive drift in the latitude- and hour-resolved wind product.

The proposal can be falsified usefully: sustained service mass may exceed the lift allocation; safe furling may be slower than damaging loads; local actuators may require more food than the route improvement supplies; cable growth may fail before the useful shear is reached; or regenerative downtime may carry the reef into conditions where it cannot recover. Those outcomes would favour shorter appendages, smaller sails, another life history, or primarily passive reefs.

## Reproducing the calculations and checking the evidence

Run `python -m research.studies.aerophytes.sail_biology` from the repository root. The default output goes to ignored `research/runs/aerophytes/sail_biology.json`; `--output` selects another destination. The curated product names its schema, producer and input hashes, scenario inputs, outputs and reading limits; the runner contains the formulas. It reads the existing aerophyte product without changing its sailing speeds or claiming a new integrated phenotype.

[sail_biology_checks.json](sail_biology_checks.json) records numerical and repository checks. Primary-source reading, biological validation and the dynamic experiment described above remain separate evidence states.
