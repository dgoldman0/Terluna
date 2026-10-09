# Resources

Where every material of the Open Moon comes from and how its cycle closes. That covers:
- the imports from the Solar-System-wide operation (water, nitrogen, oxygen) and the elements shipped in as needed;
- what lunar rock gives;
- the inventories of air, water, carbon and nutrients;
- the flows: delivery, losses to space, the shield's upkeep, construction;
- the sinks that take material out of circulation;
- the return paths that bring it back.

The domain opened on 9 October 2026 on the branch `domain/resources`. It starts from:
- what the project already holds ([state.md](state.md));
- the literature across lunar geochemistry, weathering, planetary volatiles, space resources, the carbon and nitrogen
  cycles, industrial ecology and thin films ([literature.md](literature.md), with its [sources](sources.json));
- a first screen of the committed products ([screen.py](screen.py), [results/first_screen.json](results/first_screen.json)).

Its open questions, in order, are in [questions.md](questions.md).

The 9 October cross-domain follow-up is [Returning phosphorus through the aerial biosphere](phosphorus.md):
stock versus throughput, separate natural and food-return efficiencies, a conservative six-reservoir model,
concentrate-versus-water transport, and material accounts from the common population cases. No population allocation
is selected. The earlier industry scenario allocations below are historical illustrations and are superseded for
new joint planning by [shared/scenarios/population.json](../shared/scenarios/population.json).

The author's decisions frame it:
- Water, nitrogen and oxygen come from the Solar-System-wide operation; oxygen made on the Moon is at most a
  by-product of metal and glass production.
- Imports of carbon, potassium and other elements are part of the plan. The design problem is the steady state: what
  carries each element downhill into the seas and sediments, and which local processes return it. That is the rock
  cycle's missing return limb.
- The protection is a formation of replaceable units supplied over the long term.
- Construction takes up to 500 years, and operation at least 10⁹.

## What the project already says

Each finding follows from results the repository already holds, read together with the literature. Numbers marked
*screen* come from [results/first_screen.json](results/first_screen.json). Losses carry their cycle time, the time the
loss and the resupply that matches it take to replace the inventory.

### The imports

1. **Water dominates the import, and the nitrogen comes bound in ordinary icy bodies.**
   - **Water against air.** The water to deliver is 4.4–23 times the air's mass: seas of 1.07×10¹⁹ kg, rain-fed lakes
     of 3.1×10¹⁸ kg, and 0–5.9×10¹⁹ kg in the crust's pores, depending on how deep it saturates (*screen*). The air is
     3.14×10¹⁸ kg: 2.53×10¹⁸ kg of N2 and 6.13×10¹⁷ kg of O2.
   - **The stream.** The 19 September architecture sized its matter stream for the nitrogen, 2.7×10⁸ kg/s over 300
     years. Air and water over 500 years need 1.1–4.8×10⁹ kg/s (*screen*).
   - **Where the nitrogen lies.** Molecular nitrogen on the scale needed sits only in Titan's air, about 3.5 times the
     need, and Pluto's Sputnik Planitia, 0.33–3.4 times (literature). The conservation study's principle 5 leaves
     Titan untouched, and Sputnik Planitia is a unique world. Ordinary icy bodies hold their nitrogen bound, in
     ammonia, ammonium salts and organics:
     - comet 67P holds N2/CO of 0.006 (Rubin et al. 2015), with its nitrogen in ammonium salts (Altwegg et al. 2020);
     - Bennu holds 0.23–0.25 wt% N (Glavin et al. 2025), so the inner Solar System's material of its kind would need
       42–45% of the main belt's mass;
     - the Kuiper belt holds about 2% of Earth's mass (Pitjeva & Pitjev 2018), roughly 7×10²² kg of it outside the
       large dwarf planets. At a few tenths of a percent to a few percent nitrogen, its small bodies hold 80–800 times
       the air's need (*derived*).
   - **The scale.** Supplying the air means processing 10²⁰–10²¹ kg of small bodies, the mass of about 150–1,600
     bodies 100 km across, so the ore stream runs 30–300 times the nitrogen stream (*derived*). The bodies mined for
     the water may bring much of the nitrogen with them.
   - **What follows.** The conservation principles hold as written: the nitrogen comes from ordinary icy bodies of the
     middle and outer Solar System, as the 19 September source principle and principle 5 have it. Its costs are grade
     and chemistry: dilute, bound nitrogen, cracked to N2 at the source so that the ammonia's hydrogen stays there
     (finding 2). On 9 October the author asked for the extraction to be planned in detail, feasible and sustainable
     and protecting the major celestial bodies ([register](../research/decisions.md#atmosphere)).
   - *The sources' contents are measured; the Kuiper belt's count is a rough estimate from cited values. The first
     version of this finding weighed the need against the molecular reservoirs and the main belt alone, and set the
     nitrogen against the conservation principles; counting the Kuiper belt's small bodies removed that.*
2. **Every import carries an oxygen account.**
   - **Oxygen shipped as water.** It leaves 7.7×10¹⁶ kg of hydrogen. Released, that would make the air 26% hydrogen;
     burned, it would take back every molecule of the air's oxygen, since two hydrogen molecules go with each O2
     (*screen*). The sky fleet's hydrogen is a billionth of it.
   - **Nitrogen as ammonia.** It would carry hydrogen needing 7.1 times the air's oxygen to become water.
   - **Carbon.** Carbon as organics or CO also draws oxygen (literature register).
   - **What follows.** The hydrogen stays where the volatiles are processed, and imports arrive redox-balanced: water
     split at the source and oxygen shipped as O2, or nitrogen cracked to N2 there. The form of the imported oxygen and
     the use of its hydrogen is an open question of the conservation study.
   - *Solid stoichiometry.*
3. **The fall into the Moon's own gravity heats it.**
   - **Energy per kilogram.** Mass lowered to the ground releases 1.2–1.7 MJ/kg from the exobase (1.77–2.56 lunar
     radii in the loss response's titania cases), 2.6 MJ/kg from the ring radius and 2.7 MJ/kg from a far cislunar
     terminal (*screen*).
   - **Over the build.** Dissipated on the Moon, the 500-year build's 1.7–7.6×10¹⁹ kg comes to 35–349 W/m² averaged
     over the Moon, 12% to 119% of the design sunlight (*screen*).
   - **What remote braking does not cover.** The traffic rules brake packets far away, but the fall through the Moon's
     own well remains.
   - **What it needs.** Deliveries need energy recovery on the way down (tethers, elevators, regenerative lowering), a
     longer build, or both. They also need a rule, like S7's for the protection systems, for how much of their energy
     may reach the escape region.
   - *Solid physics; the delivery scheme is open.*

### The ground

4. **Rock takes the air's oxygen far faster than space does.**
   - **The sink.** Weathering oxidizes the ferrous iron of the rock it dissolves. At the corrected runoff the land
     weathers 0.9–9.3 Gt of rock a year. At 5–16 wt% FeO (highland 5%, mare 15%; Lunar Sourcebook) that takes 5–166
     Mt of O2 a year (*screen*).
   - **Cycle time.** 4–120 million years for the air's oxygen. That is 120–3,800 times what escape removes (1.4 kg/s
     of all air in the warmest titania case, a cycle time of 70 billion years).
   - **Over the operation.** Erosion keeps supplying fresh rock: 6–200 times the air's oxygen over 10⁹ years
     (*screen*). The regolith alone holds Fe(II), metal and sulfide for 0.3–3.9% of the inventory (literature register).
   - **The only internal source** is burial of organic carbon, one O2 for each carbon atom buried.
   - **What follows.** Oxygen's long-term resupply is set by rock. The industrial architecture's step "condition
     surface materials and establish the final oxygen regime" needs this number.
   - *Solid in order of magnitude; the FeO contents are measured averages.*
5. **The carbon import is larger than the biosphere's.**
   - **The weathering sink fits field measurements.** The CO2 study's 0.36–3.6 Gt CO2 a year is 0.13–1.3 t per
     hectare a year over the land. Field trials range from 0.1 (Swiss vineyards) to 2.6 t/ha/yr (US Corn Belt), and
     rates fall with time as t^−0.61 (White & Brantley 2003).
   - **The highlands bind as much as the maria.** Highland anorthite dissolves about 2,800 times faster than albite
     (Palandri & Kharaka 2004), and highland soil holds 15.7% CaO, so it binds as much CO2 per tonne as mare soil,
     179–194 kg.
   - **Capacity beyond the biosphere.** A metre of developed soil, fully carbonated, binds 3.5 times the air's CO2.
     The seas and lakes at Earth's dissolved inorganic carbon hold 0.6–0.7 times it (*screen*). Both come on top of
     the biosphere's 1.5–1.6×10¹⁵ kg.
   - *Solid chemistry; the rates are Earth field rates.*
6. **Burial in the seas can rival weathering.** Lakes and reservoirs bury 19–48 g C m⁻² a year (Mendonça et al. 2017).
   At those rates the lunar seas would bury 0.7–1.9 Gt CO2 a year and leave the matching oxygen in the air
   (literature register).
   - **The return.** Burning carbonate back takes 4.08 GJ per tonne of CO2 from calcite and 2.57 GJ from magnesite,
     with the best kilns at 3.2–4.2 GJ per tonne of lime (JRC BREF 2013).
   - **Steering.** Weathering products steered toward magnesium carbonate return more cheaply.
   - **Size.** The return is an Earth-scale lime industry, 27–610 GW depending on the sink. It sits on the Moon,
     since it releases its CO2 into the lunar air.
   - *Earth rates.*
7. **The sea floods most of the Moon's ore.** At 28% water the atlas puts these maria 99–100% under water:
   - Imbrium, Nubium, Serenitatis, Crisium, Humorum and Nectaris;
   - Smythii, Orientale, Ingenii, Moscoviense, Cognitum and Humboldtianum.

   Procellarum is 60% flooded, Fecunditatis 65%, Frigoris 73%, Marginis 83% and Tranquillitatis 37% (*screen*). Dry
   land is mostly feldspathic highland, the Moon's aluminium and calcium ground.
   - **What lies under the water.** The high-titanium basalts, whose ilmenite gives titania and iron, and most of the
     KREEP terrane, rich in potassium, phosphorus, rare earths and thorium.
   - **What follows.** Mining them comes before flooding, in the same window as the conservation programme's
     sampling of the future seafloor. The summit's local ore is aluminium, so its 98 Mt of steel comes from the dry
     mare remnants or from space.
   - *Solid on the flooding.*
8. **The seas start fresh.**
   - **How long.** Weathering alone brings the seas and lakes to 1 g/kg in 1,100–11,000 years counting everything it
     dissolves (cations with their bicarbonate and silica), or 9,000–89,000 years by the cations alone. It takes
     0.3–3 million years to reach Earth's salinity (*screen*).
   - **What follows.** Any salt the seas need for millennia comes with the delivered water, which makes the water's
     composition a design lever. Until then the plankton are freshwater, the spray is fresh and the density is near
     1,000 kg/m³.
   - **Assumptions to revisit.** The wave studies take 1,025 kg/m³ and the water optics Earth seawater at 38.4 g/kg.
   - *Solid for weathering alone.*

### The shield's material

9. **The shield's film is the largest steady material flow in the system.**
   - **Stock and upkeep.** The ring fleet holds 3.8 Gt of titania and 34 Gt of silica, and the film plant remakes
     2.3–4.6 Gt a year, 73–147 t/s (*screen*).
   - **Against Earth's industry.** That is 0.10–0.20 Gt of titania a year, 10–20 times the world's titania pigment
     capacity, and silica on the scale of the world's cement (literature register).
   - **Fresh feed.** Whatever the plant does not recover is fresh feed: 73–147 kg/s at 99.9% recovery and 730–1,470
     kg/s at 99% (*screen*). The latter exceeds the low end of the air the shield saves, 300 kg/s lost unshielded
     (requirement S6). Over 10⁹ years it is 0.7–15 times the mass of the atmosphere (*screen*).
   - **What follows.** The recovery fraction is the shield's central resource number, and nobody has set it yet. The
     film made each year covers 26–53 lunar cross-sections, which gives the protection work a direct check on the
     film's implied life (literature register).
   - **Damage.** Silica survives proton doses a hundred times those that crystallize titania, and a silica cap protects
     the titania (Corso et al. 2024).
   - *Solid arithmetic; the recovery fraction is unknown.*
10. **The shield is best built while the Moon is airless.**
    - **Before the air.** Lifting the 46 Gt fleet from the airless Moon to the ring radius takes 2.7 MJ/kg, about
      40 GW over a century (*screen*). The mass drivers of the space-settlement studies (2.4 km/s, 3.8–6.1 MJ/kg;
      NASA SP-413) need vacuum at the ground.
    - **After the air.** A launch from the surface crosses 75 t of air per m². The air above 100 km still weighs about
      as much per square metre as Earth's whole atmosphere, and orbits below the exobase, 1,300–2,600 km up, lie
      within the air.
    - **Sequence.** The protection report already places the shield's deployment (years 150–350) before the air
      builds (250–480).
    - **Later supply.** Fresh feed comes from recycling or from small bodies processed in orbit. Reducing ilmenite
      leaves a titania-rich residue (Sargeant et al. 2020), and evaporating lunar simulant in vacuum deposits
      transparent silica films (Freundlich et al. 2005).
    - *Solid; the sequence is a proposal.*
11. **A lunar space elevator would cross the ring stack.** The rings share a node line at right angles to the Sun in
    the Moon's orbit plane. An elevator lies along the Earth–Moon line, which meets the node line at the quarter
    phases, so every ring would cross it twice a month. Rings with phased gaps or no elevator: the choice belongs to
    the shield and transport work. A lunar elevator of today's materials lifts about 584 t a year (Pearson 2005, in
    the provisioning register). *Solid geometry.*

### Nutrients and critical elements

12. **Potassium is shipped in, phosphorus is local but slow, and molybdenum is scarce.**
    - **Soils.** Lunar soils hold 0.04–0.55 wt% K2O and 0.05–0.51 wt% P2O5 (Lunar Sourcebook).
    - **Depletion.** The Moon is depleted 75% in potassium and 99% in thallium and cadmium (Taylor & Wieczorek 2014).
    - **Molybdenum.** The bulk silicate Moon holds about 0.019 ppm (Sossi et al. 2024, in the ecology register), and
      nitrogen fixation needs it.
    - **What follows.** Every micronutrient needs an element-by-element check against lunar rock, and molybdenum and
      vanadium go on the import list. *Measured, with sparse lunar molybdenum data.*
13. **The array's hardware wants elements the Moon lacks.** Germanium for multijunction collectors, indium for
    transparent electrodes, copper and rare earths for magnets and motors are scarce in lunar rock (state.md). The trim
    devices are already required to be "made of common elements" (array_industry.md). *A check to make,
    element by element.*

### Industry and the near horizons

The author asked on 9 October how the array connects to industry and to the continuing flows of the
Solar-System-wide operation, over the first thousand years with the build and the several thousand after. The study
is [industry.md](industry.md), with its screen [industry.py](industry.py); numbers here are from it.

14. **Concentrated industry flies on platforms of its own.** The bundle's one sail loading leaves an annulus tile
    2.86 g/m² beyond its film for structure, trim, store and wiring, 5.9 Gt over the fleet, and photon keeping reaches
    a quarter of the sail force. Hardware spread over every tile adds film-plant feed with its mass and glow with its
    absorption: 1.66–2.91 W/m² on the Moon for each percent of the fleet's sunlight it absorbs.
15. **Where heat is released ranks the places across a factor of a million.** A terawatt used on the Moon adds
    0.026 W/m², released at the ring radius 5.0×10⁻⁵ W/m², and at the Sun–Earth L1/L2 hubs 0.6–1.6×10⁻⁸ W/m². The
    Moon's count applies to energy brought from outside its own sunlight, wind and rivers. With its collectors' heat,
    1,000 TW of computing near the fleet puts 0.15–0.22 W/m² on the Moon; at the hubs, 5–7×10⁻⁵ W/m².
16. **Computing at the integrated ledger's scale is a gigatonne industry.** On today's servers, 12.8 kg per kW,
    1,000 TW is 15–24 Gt of hardware renewed at 2.1–3.2 Gt a year, the film plant's size. Control and interactive
    services fit beside the fleet; bulk computing fits the hubs.
17. **The film plant's tiles travel on their own sails.** Renewing 2.3–4.6 Gt a year means 3,850–7,730 tile exchanges
    a day; from one plant 6–12% of the fleet is in transit, and eight plants across the stack's tilts bring that to
    0.7–1.5%. Electric tugs would release 69–1,070 kg/s of exhaust inside the protected region, which S7 governs.
18. **On the near horizons the air needs resupply of its CO2 alone.** Escape (cycle times of 70 billion years and
    more), rock (0.008–0.24% of the oxygen by year 10,000) and burial (0.004–0.05% of the nitrogen) leave the air as
    built, while weathering and burial turn its CO2 over in 400–1,860 years. Kilns of 38–380 GW return the weathered
    carbon; the buried carbon needs a return path or 6,400–16,200 kg/s imported.
19. **The array's habitats set the largest stream after the build.** The Stanford torus shields each person with
    990 t (Johnson & Holbrow 1977), so 0.1–1 billion residents take 99–990 Gt and 2–14 billion 2,000–14,000 Gt. As
    water the whole is 0.003–0.1% of what the build delivers.
20. **After the build the operation shrinks 2,400–380,000 times.** It falls from 1.07–4.83×10⁹ kg/s, braked at the
    hubs at 20–340 PW, to 1.3×10⁴–4.5×10⁵ kg/s for the air's top-up, the film's feed, the habitats and the Moon's
    copper, unless the crust takes its water late (1.6×10⁷–1.9×10⁹ kg/s for 1,000–10,000 years) or the carbon is
    imported. The Moon's 1–6 billion people would hold 43–2,010 Gt of stocks built from lunar rock and need copper
    near Earth's identified resources. Over 10⁹ years rock takes 6.5–206 times the air's oxygen, so the long frame
    turns on local return paths.

## Phosphorus and the common population comparison

The [return study](phosphorus.md) keeps land, sea, sediment and people in the aerial phosphorus account. Its assumed
1%-coverage biosphere turns over 9.228 Mt P/year, holds 4.153–41.527 Mt P and exports 0.4614 Mt P/year in a 5% harvest.
At 99.9% natural retention and 99% harvest-P return, the low deposition/capture case still requires another
10,726 t P/year returned from below. Atmospheric deposition redistributes the Moon's existing phosphorus; it is
not an import. Finite stocks, sediment recovery and irreversible losses determine long-term persistence.

Five shared population cases explicitly separate lunar surface and aerial residents. Their generic lunar material
proxies span 258–5,025 Gt; the orbital comparator assumes 990 t shield and 15 t structure per resident. These are
accounting sensitivities. Aerial structures, gas, ecological occupancy, food and practical capacity need their own
requirements; no historical Earth-heavy allocation or orbital-heavy allocation is adopted.

## Statements elsewhere that need correcting

The evaluation ([state.md](state.md#statements-elsewhere-that-need-correcting)) lists stale or conflicting statements
in other folders. They belong to main and are unchanged on this branch. The largest:
- the CO2 study's runoff and temperature, which run its weathering sink about a quarter high;
- constant-gravity air masses about 10% low in several documents;
- the conservation principle's claim that Saturn's rings alone could meet the need, which falls short of the water
  once lakes and crust are counted;
- a unit error in the joint ledger, which gives the film plant's power in TW under a W label.

## Files

| File | Holds |
|---|---|
| [screen.py](screen.py) | The first screen: inventories, delivery heat, the oxygen's hydrogen, oxygen against rock, the shield's film, flooded ore, carbon beyond the biosphere, fresh seas |
| [results/first_screen.json](results/first_screen.json) | Its product (schema `terluna.resources.first-screen/1`) |
| [state.md](state.md) | What the project already holds for resources, what others assume of it, and statements to correct |
| [literature.md](literature.md) | The literature by topic, with what it means for the Open Moon |
| [questions.md](questions.md) | The open questions, in order |
| [sources.json](sources.json) | Every source cited, its use and how far it was read |
| [industry.md](industry.md) | The array in industry and the flows of the near horizons |
| [industry.py](industry.py) | Its screen: places, heat by place, the film plant's traffic, computing's mass, the flows after the build, people's stocks |
| [results/industry.json](results/industry.json) | Its product (schema `terluna.resources.industry/1`) |
| [phosphorus.md](phosphorus.md) | Nutrient return, finite reservoirs, transport and shared population material accounts |
| [phosphorus.py](phosphorus.py) | Conservative phosphorus stock/flow and return calculation |
| [results/phosphorus.json](results/phosphorus.json) | Its product (schema `terluna.resources.phosphorus/1`) |
| [phosphorus_sources.json](phosphorus_sources.json) | Primary scholarly sources and actual access level |
| [industry_sources.json](industry_sources.json) | The industry study's sources and how far each was read |

`python -m resources.screen`, `python -m resources.industry` and `python -m resources.phosphorus` rewrite the products; `python -m pytest resources/tests`
checks them.
