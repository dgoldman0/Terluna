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

1. **Water dominates the import, and nitrogen is the hard one.**
   - **Water against air.** The water to deliver is 4.4–23 times the air's mass: seas of 1.07×10¹⁹ kg, rain-fed lakes
     of 3.1×10¹⁸ kg, and 0–5.9×10¹⁹ kg in the crust's pores, depending on how deep it saturates (*screen*). The air is
     3.14×10¹⁸ kg: 2.53×10¹⁸ kg of N2 and 6.13×10¹⁷ kg of O2.
   - **The stream.** The 19 September architecture sized its matter stream for the nitrogen, 2.7×10⁸ kg/s over 300
     years. Air and water over 500 years need 1.1–4.8×10⁹ kg/s (*screen*).
   - **Sources.** C-type bodies and comets are rich in water and carbon. Their nitrogen is scarce and bound:
     - comet 67P holds N2/CO of 0.006 (Rubin et al. 2015), with its nitrogen in ammonium salts (Altwegg et al. 2020);
     - Bennu holds 0.23–0.25 wt% N (Glavin et al. 2025), so its kind of material would need 42–45% of the main belt's
       mass to supply the N2;
     - only Titan's air, about 3.5 times the need, and Pluto's Sputnik Planitia, 0.33–3.4 times, hold N2 on that
       scale (literature).
   - **The tension.** The conservation study's principles leave Titan's atmosphere untouched, and Sputnik Planitia is a
     unique world. The nitrogen's source is therefore the domain's first question.
   - *Solid arithmetic; the sources' contents are measured.*
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

`python -m resources.screen` rewrites the product; `python -m pytest resources/tests` checks it.
