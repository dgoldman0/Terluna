# Living sails and distributed control on sky reefs

**The strongest working proposal is a persistent, renewably maintained support system carrying replaceable sail
organs.** Each sail senses and controls itself locally, while the colony sends slower instructions about deployment and
preferred direction. Much of the sail and cable could be nonliving extracellular material maintained by a much smaller
living component. Seasonal sails, extensive retraction and independent partners stay as alternatives where their costs
or advantages differ. The author's direction of 10 October 2026 lets a sky reef be partly animal-like, coordinated by a
distributed, nerve-like system with no brain and free to depart from any Earth plant or animal body plan, and asks for
the sails' biology to be developed ([decisions](../../decisions.md#sky-reefs-partners-and-control-10-october-2026)).

**Control has to be local, because signals cross kilometres slowly.** A 5-km path takes 4.63 hours at the 0.3 m/s
measured in the epithelium of *Nanomia*'s swimming organs, 11.6 days at the 0.005 m/s order of *Mimosa*, and 50 s
even at 100 m/s. Each sail protects itself, neighbouring sails coordinate, and the colony expresses slow preferences.

**Steering effort is hinge torque and travel, and overload protection needs a different mechanism from trim.** Moving a
100-m² panel through 30° in ten seconds takes 11–29 N m of hinge torque and 5.8–15.1 J of work at the saved mean sail
loads, and 7,150 N m and 3,744 J at a 143 N/m² stress case.

**A water column across the 5-km drop holds about 8 MPa of hydrostatic head, which favours a mostly inert cable with
spaced service nodes.** Narrow living strands between pressure-isolated nodes carry signals and supplies, a better fit
to the drop than a uniformly fleshy five-kilometre stalk.

**Local furling is the candidate fast response; retraction is slow and costly.** Gathering the inherited 5.56-km cable
takes 6.43 days at 0.01 m/s and 1.54 hours at 1 m/s, and raising its 62–214-tonne bottom weight alone takes
0.50–1.72 GJ. The preferred life history keeps the cable and renews foldable blades on it.

**A finished sail's weight falls within the range of the bottom weight it hangs with.** At 0.2–0.5 kg/m² the 314,159-m²
sail weighs 63–157 tonnes (*derived*), against a bottom weight of 62–214 tonnes, so detailed sizing is worthwhile.
Renewing the sail once a year at 0.1 kg dry matter/m² builds 15.7 t C of new tissue and, with 80% of its phosphorus
recovered, exports 6.28 kg P each year.

This note is biological design ideation with explicit requirements. A sky reef is treated as one photosynthetic colonial
organism with excitable tissues, contractile organs and specialised growth zones; its local sensing, coordination and
memory need physical tissues, energy and developmental instructions. The [sailing screen](sailing.md) calculates a
steady force balance from a prescribed wing area and lift coefficients. The [requirements calculation](sail_biology.py)
and its [product](results/sail_biology.json) examine scale, signal delay, transport and replacement under stated
scenarios; numbers are *screen* values from them unless marked *derived*, and the sailing loads and inherited cable come
from the aerophyte product they read. The [scenario inputs](sail_biology_scenarios.json) and the finished sail's mass
per area are design guesses, labelled here once: the signal speeds beyond the two Earth analogues, the panel and its
arms, the retrieval speeds, the sail loadings and tissue composition, the retained water and the cable counts. A
simulated nervous system, flight stability, a viable life cycle and a complete resource budget are open.
[sail_biology_sources.json](sail_biology_sources.json) records what was read of each source and its limits.

## What the sail does

The upper reef and lower sail occupy different currents. Their relative motion provides airflow past the sail; a shaped
surface at an appropriate angle generates a sideways force, transmitted through the tether. The sail's aerodynamic lift
is approximately horizontal, and the reef's buoyancy carries the suspended mass. Equal currents remove the energy source
for this mechanism.

A local steering tail could set the main sail's angle to the apparent wind. Its action needs an aerodynamic moment
opposed by the main surface and its suspension, with a stable resulting orientation. Adjustable bridles, changing
membrane curvature and a movable trailing flap are alternatives or supplements, each a distinct mechanism to compare.
Pitch, roll, yaw, cable swing and interactions between neighbouring sails must all stay bounded.

A flexible trunk chiefly transmits tension; turning the reef need not turn a freely swivelling sail kilometres below it.
The engineered StratoSail concept places its rudder actuator at the wing, and its author draft describes a ribbed
flexible skin, a boom and local controls. Its limited model tests and simulated routes make it a mechanical precedent
for the living design.
[Aaron, Heun and Nock](https://gaerospace.com/projects/ULBDStratoSail/pdfs_docs/COSPAR2000StratoSail.pdf).

The sailing screen's 0.4–0.8 m/s across the shear holds for a surface that can select and hold its angle. A real living
surface could have less lift, more drag, a smaller safe angle range and considerable time unavailable for steering. The
sketch Codex drew in the conversation has six sails for illustration; the number of sails is open, and the aerodynamics
remain those of the sailing screen's single wing.

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

A membrane can be sparsely cellular in places. Thin repair tissue could secrete and maintain a tougher extracellular
sheet, while veins carry tension and service small patches. This follows the existing [carbon allocation](biology.md):
inert material incurs construction and renewal costs, and active tissue also respires. A hollow or pressurised rib is
another candidate; holding its pressure costs seals, fluid and the loss through those seals, so its stiffness comes with
that mass and upkeep.

The main cable should have distinct structural and service parts. A dead fibre bundle carries load, and living strands
carry the electrical signals and nutrients. A thick wet nerve or vascular trunk is sized with its own strength and dry
mass, which differ from an ideal fibre's. Protective layers, living sleeves and reserve stores belong in the suspended
mass and renewal accounts.

## Control without a central brain

Earth plants provide precedents for long-distance electrical and calcium signalling, signal-triggered movement and
resource integration between connected modules, all without neurons or a guidance system. Colonial animals provide a
complementary precedent for specialised motor organs and signalling among them. A proposed reef may use excitable
epithelia, nerve-like cords, local ganglia-like knots, chemical signals and mechanical feedback in different places.

Wound-triggered glutamate/calcium signalling in plants and induced defence among connected clover ramets support the
idea of differentiated responses across a connected colony. The clover network's source–sink relations also mean that
modules may communicate unequally. [Toyota et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30213912/),
[Gómez et al. 2008](https://doi.org/10.1111/j.1469-8137.2008.02542.x).

The siphonophore *Nanomia bijuga* is especially useful as an organisational analogue: its swimming organs combine nerve
rings, colony-connected tracts, excitable epithelia and muscle, and some local responses persist after the main tracts
are interrupted. *Hydra* imaging also distinguishes functional neuronal networks associated with different behaviours.
Both support examining more organisation than a single undifferentiated alarm network, and the reef's controller remains
to be designed. [Norekian and Meech 2020](https://doi.org/10.1242/jeb.233494),
[Dupre and Yuste 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5423359/).

Three levels of control suit the physical scale:

1. **Local protection.** Strain, hinge movement and apparent airflow affect the nearby actuator or a passive release. An
   overloaded panel reduces its exposed area or changes angle. Passive geometry should reduce load during a signal
   failure, and that response needs measuring.
2. **Sail coordination.** Neighbouring organs share tension and angle information. They change tack gradually, avoid
   sweeping through one another and report sustained trouble. A heavily loaded branch overrides a request for more
   steering.
3. **Colony behaviour.** Slower signals express a preferred side, a deployment limit or a need for water. Changes in
   light, hydration, temperature, reserves and repeated load history alter that preference. Hinge angles stay with the
   sails.

This architecture avoids a single decision point whose damage disables the whole rig. It introduces coordination
problems: neighbouring controllers can oppose one another, delayed corrections can drive oscillations, and repair can
reconnect a signal path with the wrong gain or timing. Refractory intervals, dead bands and limits on the rate of
adjustment are candidate biological functions; their values must come from stability analysis and measured tissue
response.

A useful local controller senses the **apparent wind** at its sail. Hairs or flexible sensory projections could measure
bending; paired pressure-sensitive patches could distinguish sides of the flow; strain cells could report tendon load;
hydration and chemical signals could report resource state. Venus-flytrap experiments distinguish mechanical reception
from threshold-dependent electrical/calcium propagation, with the MSL10 channel affecting weak-stimulus sensitivity,
which supports examining sensory thresholds. [Suda et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485109/).
The sensors' airflow sensitivity, fouling and mechanical interaction with the airfoil remain requirements to set. A
local light compass could provide a reference direction, though the relationship between light, time and position
changes along a drifting route.

The colony could test small sustained changes and retain a chemical or neural-like memory of whether hydration or light
income improved. That is an exploratory control proposal; learning, consciousness, precise geographic navigation,
advance storm prediction and knowledge of where the nearest coast lies are outside it. Nutrient-importing visitors might
provide cues, though their presence is an unreliable coordinate system.

## Signal distance and actuator effort

The table gives one-way transit times from `distance / speed`. The 5-km path is a round-distance illustration; the
inherited cable for a 5-km vertical drop is 5.56 km long and would take about 11% longer to traverse (*derived*).

| Assumed propagation speed | One-metre path | Five-kilometre path | Status |
|---|---:|---:|---|
| 0.005 m/s | 200 s | 11.6 days | Order of the small-organ electrical/calcium rates reported in *Mimosa* |
| 0.3 m/s | 3.33 s | 4.63 hours | Measured order for *Nanomia* nectophore epithelium; the stem nerve is a separate pathway |
| 1 m/s | 1 s | 83.3 minutes | Design sensitivity |
| 10 m/s | 0.1 s | 8.33 minutes | Design sensitivity |
| 100 m/s | 0.01 s | 50 s | Design sensitivity |

*Mimosa* measurements were tissue- and stimulus-specific, about 4–6 mm/s in the reported rachilla experiments. The
*Nanomia* epithelial measurement was approximately 0.3 m/s at 11.5–15 °C. Extending either rate over kilometres shows
the distance penalty. [Hagihara et al. 2022](https://doi.org/10.1038/s41467-022-34106-x),
[Norekian and Meech 2020](https://doi.org/10.1242/jeb.233494).

An organ's response adds sensing, synapses or relays, muscle response and the motion of the sail itself to these times.
Faster dedicated cords could improve communication; whether any of these tissues can survive, conduct and be maintained
along a wet, loaded kilometre-scale cable is open. Splitting the route into relays protects signal integrity, and the
total propagation time remains.

Local control also needs a suitably small physical organ. Enlarging a millimetre-scale sensitive-plant leaf into a
hundred-metre flap changes its inertia, stiffness and fluid demand. A distributed row of shorter flaps may shorten
control paths, at the cost of more living tissue and coordination. Passive structural damping and automatic unloading
stay valuable even with faster nerves.

There are several actuator candidates. Contractile tissue can make repeated reversible changes but needs oxygen, fuel
and a way to hold position economically. Osmotic or hydraulic cells can set tension and curvature, with water
availability and transport limiting response. Elastic structures with a latch can switch quickly after slow charging,
while reset costs and fatigue limit repetition. Differential growth suits persistent trim and seasonal remodelling, and
gusts need the faster mechanisms. Flytrap snap-buckling supplies a small-scale precedent for releasing stored elastic
energy; safe snap-through of a large sail, and resetting as fast, remain to be shown.
[Forterre et al. 2005](https://doi.org/10.1038/nature03185). Combining these mechanisms is more credible than
assigning rapid, strong, low-cost motion to one generic tissue.

Steering effort depends on **hinge torque and travel** more than on the whole tether force. A small control surface can
steer a larger surface through the airflow, and at some angles its actuator may still need to be large. The requirements
product illustrates how aerodynamic moment, moment arm and muscle attachment determine force and work. Hinge imbalance,
friction, elastic preload, stall, dynamic overshoot and sustained muscle metabolism can all increase the cost. Actuator
mass and operating efficiency are open.

For scale, a 100-m² panel with a 0.5-m force arm and 0.25-m tendon arm, moved through 30° in ten seconds under constant
resisting torque, gives:

| Assumed resultant loading | Hinge torque | Tendon force | Mechanical work per motion |
|---|---:|---:|---:|
| Saved mean sail loads, 0.22–0.58 N/m² | 11–29 N m | 44–116 N | 5.8–15.1 J |
| Selected stress case, 143 N/m² | 7,150 N m | 28,600 N | 3,744 J |

This is a torque requirement for a uniform load with chosen moment arms, area, angular travel and timing; holding,
inertia and biological inefficiency come on top. The stress case makes clear why ordinary trim and overload protection
need different mechanisms.

## Behaviour through a change in weather

| Situation | Candidate response | What must be established |
|---|---|---|
| Suitable shear and adequate reserves | Expose some panels, stabilise their angle, then increase area within load limits | Positive benefit after service costs; stable wing and tether dynamics |
| Desired change of side | Reduce force, change local trim, pass through a low-force state, then rebuild force | Safe transition without entanglement, shock loading or loss of height |
| Weak or reversing shear | Reassess angle and reduce unnecessary area; accept drift when steering is ineffective | Reliable flow sensing; some winds may not permit the requested heading |
| Rising or oscillating load | Local feathering or folding, followed by branch isolation if needed | Unloading actually lowers force; reaction is fast enough under load |
| Forecast or sustained severe weather | Long-term shortening where affordable, protective stowage, or avoidance by route | Adequate warning and a demonstrated retrieval mechanism |
| Control or service failure | A passive low-load state, redundant service route, or sacrificial loss of a small part | The failed part stops flapping violently as a drag surface |
| Long darkness | Retain essential reflexes and tissue maintenance; reduce discretionary collecting and steering | Local fuel and oxygen balance throughout the longest relevant dark interval |

Feathering means reducing aerodynamic force by orientation, furling means reducing exposed membrane area, and retraction
means bringing the suspended assembly upward. Each has its own mechanics, response time and cost. Even turned edge-on,
ribs, lines and folded material stay exposed to drag.

Releasing a sail or cable relieves the reef of weight as well as pull, which can trigger an unwanted ascent, recoil and
a change in the load distribution among the remaining modules. A breakaway design therefore needs a corresponding
buoyancy and water response, since losing a heavy sail can itself do harm. Releasing a branch may also export its
nutrient stock and create a hazard below.

## Supplying an organ kilometres below its parent

The lower sail needs water, carbon, minerals, oxygen and repair materials. Local reserves carry it through short
interruptions of supply. Photosynthetic patches or symbionts could contribute carbon when the sail receives useful
light, though the upper reef can shade it and the sail's orientation is chosen for force. More photosynthetic tissue
also adds wet mass and maintenance. The sail receives only the light the reef above lets through, and the light account
counts it once.

A continuous water-filled conduit across a 5-km vertical separation develops roughly **8 MPa of hydrostatic pressure
difference** under lunar gravity. A downward supply can release gravitational potential energy and still needs
containment and flow regulation; returning water upward requires work. In a closed circulation the gravitational gains
and losses can offset in the ideal energy balance, while local pressure, friction and loss through seals remain. A
hanging organism pays the full physical costs of long-distance fluid transport.

A candidate arrangement is a series of pressure-isolated reservoirs and service nodes, with valves or pumps between
them. True pressure breaks require chambers and controlled exchange; valves along one uninterrupted filled tube leave
its hydrostatic head unchanged. More stages trade smaller pressure differences for added walls, fluid mass, failure
points and metabolic work. Separate downfeed and return paths would need resistance and flow budgets. Local collection
and use of water reduce transfer demand, and nutrients reach the upper reef only through a return path.

Experiments on connected *Alternanthera* ramets under heterogeneous water availability demonstrated resource integration
and unequal allocation benefits over an eight-week study. They support connected resource sharing as a biological
concept; flow capacity and pressure tolerance over kilometres are open.
[You et al. 2016](https://doi.org/10.1038/srep29767).

These costs make a mostly inert cable with spaced service nodes and narrow connecting living strands more attractive
than a uniformly fleshy five-kilometre stalk. Electrical biological signalling still needs a continuous hydrated
excitable pathway between nodes, and inert gaps need a different communication mechanism. Motor and control organs also
need gas exchange while retaining water, so each protective sheath trades protection from drying against oxygen supply
to the active tissue inside. A fully separate lower organism could feed itself more independently, and then reliable
attachment, communication, reproduction and reciprocal benefit become a symbiosis problem. The choice between these
architectures rests on more than tensile strength.

## Growing a complete rig

A juvenile can begin as a drifting buoyant unit with short sensory or collecting appendages, and it can do without the
mature reef's steering ability. A short cable samples little vertical shear, so a juvenile with the adult's ratio of
sail area to reef area sails slower than the adult. If juvenile survival requires steering, that becomes an additional
nursery or juvenile-organ requirement.

An initial construction sequence worth developing is:

1. **Establish flotation and reserves.** Grow enough productive and buoyant modules to carry the next appendage
   increment while meeting repair, tenant and dark-period obligations.
2. **Build an attachment nursery.** Develop several load paths and protected growth tissue on the underside. Keep
   control and service connections separate from gas containment.
3. **Deploy a small pilot organ.** Grow a short supported pendant with a modest folded sail and local sensors. Its early
   benefit may be sensing or collecting; useful steering is a later possibility to calculate.
4. **Extend in affordable increments.** Lengthen the structural bundle alongside its service strands, strengthening the
   upper cable as the suspended load increases. Each intermediate geometry needs its own load and carbon allowance.
5. **Bud blades near the lower hub.** Establish the leading rib, control organ, supply connection and membrane supports
   before exposing a large new area. Open sectors gradually as their tissues mature.
6. **Rotate renewal.** Grow replacements while neighbouring blades remain useful, recover accessible nutrients from
   retiring tissue, and release only the residual parts.
7. **Provision reproduction or fission.** Each daughter needs working growth zones, stores and enough independent
   flotation. A cable shared between daughters that separate leaves the inheritance scheme incomplete.

There are several ways to lengthen the cable, and calling it a stem selects none of them:

| Growth mechanism | Advantage | Difficult part |
|---|---|---|
| Distal growth near the sail | No large store of folded cable at the reef | Supply reaches an increasingly remote bud; newly formed material must mature before carrying large tension |
| Growth at many intermediate nodes | Growth and repair can occur near damage | Soft expanding regions interrupt a highly stressed tension member unless bypass fibres carry the load |
| Proximal production and controlled payout | Construction occurs near the main carbon source | Alternating grips or another load-transfer organ must let new material enter service without losing support or circulation |
| Deployment of pre-grown folds | Opening can be faster than new biosynthesis | The stored cable was already paid for and carried; folds require space, protection and controlled release |
| Serial budding of service modules | Repeated units can contain reserves, valves and repair tissue | Junctions, sealing and coordinated load transfer become principal failure sites |

The working cable hypothesis is an overlapping fibre bundle with protected service nodes. New strands take load before
old strands are retired, so the load-bearing bundle persists while its constituents change. How fibres attach, mature
and transfer tension, and whether soft growth tissue can safely elongate an already loaded kilometre-scale cable, are
open.

In *Nanomia*, spatially restricted proliferation and candidate stem-cell populations support ordered budding and
specialisation, and regeneration throughout a mature colony would need evidence of its own. Retained renewal sites are
therefore a specific proposed organ. A separate colonial tunicate experiment demonstrated whole-body regeneration from
small vascular fragments, an alternative biological precedent; a growth rate for an aerial cable remains to be found.
[Siebert et al. 2015](https://doi.org/10.1186/s13227-015-0018-2),
[Rinkevich et al. 2007](https://doi.org/10.1371/journal.pbio.0050071).

## Permanent structure and replaceable parts

Permanence can apply to a location, a function or the same physical material. A reef can retain a steering hub for
decades while repeatedly replacing its cells and fibres. It can also keep the costly cable while replacing only its
exposed membranes. These distinctions matter more than choosing a single lifespan for the entire organ.

| Life history | Best reason to consider it | Main cost or risk |
|---|---|---|
| Persistent cable with renewable, foldable blades | Keeps access to useful shear and avoids repeatedly building the support system | Cable and service tissue remain exposed even when blades are furled |
| Continuously exposed, continuously renewed rig | Available steering and sensing with few deployment transitions | Chronic drag, wetting, predation and storm exposure; no protected interval |
| Seasonal blades on a persistent cable | Reduces active surface and maintenance during a long unfavourable interval | Regrowth and nutrient recovery take time; reduced steering may prolong the bad conditions |
| Fully retractable rig | Changes working depth and brings distal parts close for repair | Hauling work, storage, bending damage, concentrated load and a substantial retrieval organ |
| Whole rig shed and regrown | Can isolate a severe failure without a retrieval system | Repeated cable construction, nutrient export and a potentially long period without steering |
| Separate living sail partner | Distal metabolism and specialised motor biology could be more autonomous | Recruitment, compatibility, reciprocal benefit and reliable communication must evolve or be designed together |
| Detachable reproductive appendage | Planned retirement might also disperse a provisioned offspring | Independent flotation, reserves and development must be supplied; once detached it cannot steer its parent |

The preferred starting variant is the first: retain attachment and renewal tissue, renew the cable through overlapping
load paths, and replace membrane bays more readily than sensory and motor organs. It preserves expensive depth access
while allowing ordinary damage to remain local. This preference could change if persistent cable drag, supply costs or
frequent entanglement dominate.

Seasonal shedding is more plausible when a long predictable interval offers little steering benefit. A simple rule to
shed every night is poorly supported: wind can power sailing in darkness, and navigation may help the organism avoid
still worse conditions. Fast release belongs to failure states, with planned retirement allowing more time to reclaim
resources.

For the inherited 5.56-km cable, selected constant retrieval speeds imply:

| Retrieval speed | Time to gather the paid-out line |
|---|---:|
| 0.01 m/s | 6.43 days |
| 0.1 m/s | 15.4 hours |
| 1 m/s | 1.54 hours |

These are arithmetic scenarios. At fixed reef altitude, raising only the inherited 62–214-tonne bottom weight through
the 5-km vertical drop requires **0.50–1.72 GJ**. Raising the cable, drag, friction, losses, storage and the hauling
organ all add costs. Folding a sail leaves its weight in place; shortening the whole cable is a separate operation.

The engineered StratoSail draft left powered reel-up out of its then-current design to save mass. A biological hauling
organ, long-term shortening through tissue remodelling and passive payout from stored folds are distinct development
choices. Local furling remains the candidate fast response; a whole-rig retrieval rate for a living rig is open.

## Construction costs and phosphorus recovery

The representative 2-km circular reef has 3.14 million m² of footprint. A sail at the assumed 10% share has 314,159 m²
of area, equivalent in area to a 560-m square (*derived*). In the parent product's 5-km-drop cases, the prescribed sail
force is about 70–182 kN. Its mean force per sail area is only 0.22–0.58 N/m², while the selected tether and bottom
weight together range from about 97 to 328 tonnes across the two drag and fibre-strength choices. These are ideal
support cases that leave the reef's own anatomy out.

The bottom weight alone ranges from 62 to 214 tonnes. A finished sail averaging 0.2–0.5 kg/m² would weigh
63–157 tonnes (*derived*). That overlap makes detailed sizing worthwhile; whether wet membrane, ribs, muscles, conduits
and fittings fit the budget is open. Sail tissue can replace some required ballast only where its location and load path
provide the needed weight and balance, and a distributed membrane acts differently from the model's point weight. A
change in shape, drag, total mass or centre of mass requires the hanging geometry and host lift account to be
recalculated. Across this sail, 0.1 kg/m² of additional retained water adds 31.4 tonnes; 1 kg/m² adds 314 tonnes. A
collecting fringe therefore needs drainage and a wet-load limit alongside its capture function.

The prior sizing combines simplified cable drag with an idealised hanging shape; the cable in distributed wind and
gravity together is unsolved, and attachments to the reef, membrane stiffness and actuator torque are outside its
result. The separate [structure screen](structure.md) gives a 143 Pa design gust at 5 km: multiplying that by the
full sail area and an assumed force coefficient of one gives about 45 MN if the entire surface meets that pressure
(*derived*). This is a load warning built on that whole-surface case; the sail's simulated response is open.

The number of sails has a cost. Dividing a fixed ideal tensile load among N identical circular cables of the same
material and length keeps their summed cross-sectional area constant but increases their summed diameter, and hence
projected cylindrical drag area, by sqrt(N). Real cables also have their own weight and unequal loads. A few shared
trunks with shorter terminal branches are worth comparing with many independent long lines: a single trunk concentrates
failure, while many trunks increase drag and entanglement paths. The best arrangement for a deformable reef is open.

Each component needs its own dry mass, living fraction, water content, tissue phosphorus fraction and renewal schedule.
Structural polymer can be relatively nutrient-poor while motors, growing buds and service tissue are more demanding. In
a steady inventory, planned phosphorus export from component i is

`P_loss_i = dry_mass_i × P_fraction_i × (1 − recovered_fraction_i) / renewal_interval_i`.

The recovered fraction here means phosphorus retained inside the colony boundary and available for reuse; phosphorus
moved from a blade into a hub that is shed immediately afterward leaves the boundary with the hub. Emergency loss needs
a separate occurrence rate and recovery fraction, and each event is charged once, to planned or to emergency
replacement. Carbon construction, nutrient stock and recurring export are distinct quantities.

The requirements product varies dry sail loading and turnover at a single explicit tissue composition; it treats living
and structural tissue together and leaves emergency losses out. At 0.1 kg dry matter/m² of sail and 50% carbon by dry
mass, the 314,159-m² surface contains 15.7 tonnes of carbon. With a phosphorus fraction of 0.001 and 80% recovery
before the residual tissue leaves the reef, replacement gives:

| Replacement schedule | Carbon incorporated in new tissue | Phosphorus exported in unrecovered material |
|---|---:|---:|
| 10% annually | 1.57 t C/year | 0.63 kg P/year |
| Complete annual replacement | 15.7 t C/year | 6.28 kg P/year |
| Complete replacement every lunar month | 194 t C/year | 77.7 kg P/year |

These are material flows; total synthesis and maintenance costs, external carbon demand, affordable renewal rates and
measured recovery are open. The product also records 0.02 and 0.5 kg dry matter/m² sensitivities. Emergency shedding
with little recovery would export more phosphorus per kilogram lost, while old tissue that stays inside the reef's food
web keeps its phosphorus within the colony. Resorbed phosphorus can help build the next blade, which still costs its
carbon construction, water, transport and maturation.

Petioles of *Oxalis pes-caprae* provide a mechanical release precedent: preferential breakage occurred at a basal weak
zone, and replacement leaves came from surviving meristems. That supports designing both a break plane and surviving
renewal tissue; the release there is passive, and an actively commanded joint remains to be shown. Nutrient resorption
varies with species and stress, so emergency detachment is credited with incomplete nutrient recovery.
[Shtein et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6408339/),
[Urbina et al. 2021](https://doi.org/10.1002/ece3.7734), [Xu et al. 2020](https://doi.org/10.1093/jpe/rtaa053).

Budget allocation must close at the whole reef. Extra sail tissue, innervation, collectors and repair draw from the same
surplus used for tenants, gas replacement, reproduction and other structures. A sail that imports phosphorus but
consumes too much carbon, water or structural reserve can still reduce host persistence. The
[symbioses analysis](../sky_ecology/symbioses.md) draws the same distinction between a useful service and net benefit.

## Additional functions and their conflicts

A useful sail can do more than steer, though some combinations work best in different tissues or at different times.

| Additional function | Candidate location or mode | Benefit to examine | Tradeoff that must remain in the account |
|---|---|---|---|
| Flow and strain sensing | Leading edges, suspension roots and a small exposed pilot surface | Sensing already serves control and can inform deployment | A point sensor may miss a gust elsewhere; fouling changes calibration |
| Cloud-water capture | Replaceable fringes outside the clean airfoil and small local stores | Supplies distal tissues and may reduce downfeed | Droplets add weight and roughness; transporting excess upward costs work |
| Particle or prey collection | Dedicated retractable combs on separate edges | Net imports require externally sourced material; recaptured reef litter retains existing nutrients | Capture rate, digestibility and phosphorus content are unknown; drag and fouling can erase the gain |
| Photosynthesis | Selected illuminated ribs or separate small blades | Supports local metabolism in favourable light | Reef shade, sail angle and dark intervals; living wet mass and nutrient demand |
| Heat exchange | Service nodes with controllable exposed area | May help a locally warm control organ | Long-range cooling of the parent needs circulation; unwanted night cooling and water loss |
| Local defence | Guard chambers near permanent hubs and protected tendons | Protects vulnerable joints and maintenance tissue | Occupants consume food, add mass and may foul moving parts |
| Reproductive dispersal | A specifically provisioned bud near a renewal zone | Shares construction infrastructure with organ growth | A viable propagule needs provisioning beyond an expendable sheet |
| Stable signalling surface | Pigment or orientation changes visible to neighbouring organs or partners | Could supplement slower internal communication | Occlusion, darkness and required sensory machinery; the long-range signalling rate is open |

Cactus experiments demonstrate fog-droplet collection and directional movement over spines into absorptive tissue. They
support patterned collectors on short scales; liquid water from unsaturated air, and lifting it kilometres to the
parent, lie beyond that evidence. [Ju et al. 2012](https://doi.org/10.1038/ncomms2253). In the reef, lower collectors
should first support lower organs; export to the upper colony remains a separate transport calculation.

A reasonable division would keep the main aerodynamic membrane relatively smooth, put collecting hairs on replaceable
fringes, and protect muscles and growing tissue near permanent junctions. Capturing water and prey could be favoured
when the sail is feathered; high-force steering might favour clean surfaces and reduced collecting area. Mode changes
must preserve the combined mass and nutrient ledger.

Large roosting animals and persistent pools belong on stable reef modules, away from frequently moving or sacrificial
blades. Small guards, decomposers and surface cleaners could occupy protected hubs, provided their upkeep keeps hinges
clear and adds no uncontrolled wet mass. Cleaning that throws nutrient-bearing material into the air defeats retention
unless that export is paid for. A hanging sail also offers predators an approach to sensitive tissue, and its defence
depends on its guards and tissues whatever its size.

## One organism or a partnership

For the first design, the float modules, cable, control knots and sail organs can be differentiated parts of one
genotype. Their local autonomy is physiological, and their development can share inherited rules. A tissue containing
genuine neurons would be proposed neural tissue, while a calcium-signalling epithelium would remain a different
mechanism. Their costs can be compared without settling an Earth kingdom label, as the author's direction allows
([decisions](../../decisions.md#sky-reefs-partners-and-control-10-october-2026)).

Photosynthetic symbionts and external guards can still be different species. A separate sail organism is a richer
alternative: the reef offers attachment and carbohydrate while its partner offers a motor organ, local sensing and
perhaps imported resources. Reliable inheritance or recruitment, signal interpretation, load sharing, recognition and
conflict resolution become additional requirements. A partner's incentive to detach and reproduce may conflict with the
host's need for steering during a storm.

Fusion of independently founded reef modules could improve survival and construction; compatible circulation, signalling
and reproductive interests need mechanisms beyond contact. Even a single-genotype colony needs controls on pathological
proliferation and local resource monopolies. Protected growth zones should be regulated by whole-colony supply and load
capacity.

## Forms worth carrying forward

Three contrasting ecotypes would expose useful differences without choosing one universal reef:

- **Persistent open-air sailor.** A small number of long support trunks, sparse living service tissue and repeatedly
  renewed low-mass panels. It invests heavily in local stability and retains an ability to steer through darkness.
- **Cloud-edge collector.** More collecting fringes, local water storage and short-range nutrient processing, with fewer
  assumptions about returning water to the upper reef. Its additional wet load and drag may require less exposed sail or
  a smaller depth separation.
- **Seasonal or migratory sailor.** Expands blade area ahead of a sustained useful wind regime and resorbs it during a
  long low-benefit interval. It retains control and growth hubs, with reduced steering while the surface regrows.
  Forecasting can be an inherited seasonal programme, but its success must be tested against actual routes.

A low coastal form using a sea drogue is a different mechanical ecology: the second medium is water, so it may obtain
useful force without kilometre-scale air shear. It brings problems of its own: wave loads, immersion, salt and terrain
clearance. Its [sailing comparison](sailing.md) stays separate until an anatomy is specified.

The promising visual anatomy is therefore a reef with durable suspension roots and distinct organs at their ends: some
narrow clean steering blades, some collecting fringes, folded panels resting beside new buds, and older panels being
retired. The next mechanical comparison sets their number, dimensions and tissue mix.

## Requirements for the next model

The next useful calculation is a single shaped sail with its suspension, tail, finite actuator and local delay in
changing airflow. It should report angle error, loads, oscillation growth, control effort and behaviour after a loss of
power. Matching a commanded sign of lift is only the first test. Vary pivot and centre-of-mass positions, tail area and
arm, material damping, membrane deformation and the mode used to reduce area.

Then connect several organs through compliant cables and a deformable buoyant reef. Include line drag on the actual
hanging geometry, unequal gusts, wave propagation along the cable, mutual wakes and release of a loaded branch. An
acceptable local controller can still make a coupled system oscillate or twist itself together.

A biological budget follows alongside that mechanics: measured or explicit hypothetical dry and wet masses, turnover,
actuator metabolism, fluid pressures and flow, construction, phosphorus recovery and night reserves. Trace every stage
from juvenile through cable extension, mature renewal, failure and reproduction. Once those constraints can coexist,
active trajectories can be compared with passive drift in the latitude- and hour-resolved wind product.

The proposal can be falsified usefully: sustained service mass may exceed the lift allocation; safe furling may be
slower than damaging loads; local actuators may require more food than the route improvement supplies; cable growth may
fail before the useful shear is reached; or regenerative downtime may carry the reef into conditions where it cannot
recover. Those outcomes would favour shorter appendages, smaller sails, another life history, or primarily passive
reefs.

## Reproducing the calculations and checking the evidence

Run `python -m research.studies.aerophytes.sail_biology` from the repository root. The default output goes to ignored
`research/runs/aerophytes/sail_biology.json`; `--output` selects another destination. The curated product names its
schema, producer and input hashes, scenario inputs, outputs and reading limits; the runner contains the formulas. It
reads the existing aerophyte product and leaves its sailing speeds and integrated phenotype unchanged.

[sail_biology_checks.json](sail_biology_checks.json) records numerical and repository checks. Primary-source reading,
biological validation and the dynamic experiment described above are separate evidence states.
