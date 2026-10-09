# Provisioning

How the Open Moon turns resources and energy into what people and the array need. That covers:
- energy, its storage through the night, and where its heat goes;
- food;
- building and making things;
- computing;
- transport;
- water, sanitation and the services cities run on;
- maintenance and repair.

The word reaches past industry to everything a settled world provides itself with. The domain opened on 9 October
2026 on the branch `domain/provisioning`. It starts from:
- what the project already holds ([state.md](state.md));
- the literature on human needs and their energy, urban metabolism, long-duration storage, power from space, waste
  heat, food systems, materials, computing, transport, maintenance, life support and water
  ([literature.md](literature.md), with its [sources](sources.json));
- a first screen of the committed products ([screen.py](screen.py), [results/first_screen.json](results/first_screen.json)).

Its open questions, in order, are in [questions.md](questions.md).

**Population follow-up, 9 October 2026:** the original Earth-heavy population placements are historical
illustrations, not capacity findings. The [inhabited-volume study](../research/studies/inhabited_volume/README.md)
now evaluates lunar surface and aerial habitation together using
[shared comparison cases](../shared/scenarios/population.json). It quantifies residential mass, buoyancy, projected
area, holding-power sensitivity, domestic water recovery, food and heat. No allocation or maximum capacity is
adopted, and runoff stress indices are not physical population ceilings.

The author's decisions frame it:
- **The Moon's cycle:** it keeps its 29.53-day cycle, and life stores resources across the night.
- **Power:** fusion supplies most local power, supplemented by regional solar and wind (a planning assumption).
- **The summit:** a commercial port of about 1.5 million people, in a beautiful, sustainable metropolis of about 100
  million within 32 km. Regional agriculture and foraging feed it.
- **Sky ships:** the Moon is a sky-ship world where lift is cheap, and the flyers stack by height.
- **Orbit:** energy-intensive industry and computing belong in orbit, and the planned branch `research/array-industry`
  takes the array's heat, computing and industry.
- **The night:** the ring fleet keeps its light off the night zones, redirecting it and using it where it can, for
  computing for instance.

## What the project already says

Each finding follows from results the repository already holds, read together with the literature. Numbers marked
*screen* come from [results/first_screen.json](results/first_screen.json).

### Heat and energy

1. **The ring fleet's own infrared heats the Moon more than all human use could.**
   - **The flow.** The fleet's films absorb 3–4.5% of the sunlight they pass, 65–83 PW in all, warm to about 175 K
     and radiate from both faces. The tiles stay radial-facing, so the same face looks at the Moon all orbit, and the
     films emit a little more there than sunward: one layer sends 56% of its heat inward. A ring hidden behind another
     absorbs that ring's glow and passes it on through both its faces, so the stack, about two layers deep in the
     Moon's view, sends between half and 56% inward. A tile 19,355–21,392 km out sends the Moon 0.66–0.81% of what it
     emits inward. The Moon receives 5.7–9.9 W/m² averaged over its surface (*screen*).
   - **In proportion.** That is 1.9–3.4% of the design sunlight, the heat of 215–376 TW used on the Moon.
   - **The warming.** If it acts like sunlight, the GCM's 1.9 K per 1% of sunlight puts the design case 3.7–6.4 K
     warmer (*screen*). The design case of 294.9 K leaves it out. The 5% dimmer, chosen to bring the Moon from 303 to
     294 K, removed 15.4 W/m² of sunlight at the top of the air.
   - **The decision.** On 9 October the author made the investigated climate the planned one, so the glow is answered
     with added solar protection, by redirecting it away from the Moon, or both
     ([register](../research/decisions.md#solar-shield-and-habitat-array)).
   - **The means.**
     - Films that absorb less cut the glow however the stack lies; most of what they absorb is the ultraviolet the
       titania stops.
     - A heat mirror on Moon-facing faces, clear to sunlight and reflective in the thermal infrared, sends a layer's
       glow outward, and the stack weakens it. By a two-layer estimate (*screen*), a mirror of emissivity 0.03–0.1 on
       every tile cuts the Moon's share 1.5–2 times and traps the hidden layer between two weak emitters at
       280–380 K. On the layer nearest the Moon alone it cuts the share 2–5 times, and 3–9 times with strong emitters
       on the outer faces. The films absorb only 3–4.5% of sunlight, so a mirror that absorbs a few percent more adds
       as much heat as it turns away; a sparse metal mesh or metasurface a few micrometres in pitch is the first
       candidate.
     - Added dimming answers what remains.
   - **Computing and industry against the glow: poor.** The Moon receives 0.33–0.45% of the heat the films absorb,
     so keeping a watt of glow off it means taking 221–303 W out of the films and rejecting it where the Moon does
     not see it (*screen*). Halving the glow that way would take 32,600–41,600 TW of computing; the world's data
     centres draw 0.047 TW (IEA 2025), and 100 TW of computing in orbit would cut the glow by 0.12–0.15%. The films'
     heat holds little work: Carnot across every film gives 0.25 TW, and the best measured thermoradiative diodes over
     every film 27 TW ([array_industry.md](../research/studies/solar_shield_array/array_industry.md)), and those
     narrow-gap devices would also absorb the sunlight the films pass. Titania made photovoltaic could carry part of
     its absorbed ultraviolet to outward radiators as current, but the contacts add absorption to films that absorb
     3–4.5% now, and without users the current is only heat moved elsewhere. Thermodynamic and optical computers
     raise the computing each watt buys; the heat still leaves wherever their radiators send it.
   - **Computing and industry on the light: significant, limited by demand.** The night half intercepts about
     850 PW the Moon never needs (joint synthesis), and its attitude keeps that light off the night. Each terawatt of
     computing takes 2,450–3,670 km² of collectors and 830 km² of radiators, so 1,000 TW, about 20,000 times today's
     data centres, takes 0.12–0.16% of the fleet's area (*screen*). With its collectors' heat, 2.9–4.4 W for each watt
     of computing, it puts 0.15–0.22 W/m² on the Moon released near the fleet, 5–7×10⁻⁵ W/m² from the Sun–Earth
     hubs (*screen*; the resources domain's industry screen), and less from radiators facing outward or edge-on to it. Energy is no limit there; demand, hardware mass
     and latency are, and computing or industry used for nothing is only waste heat.
   - **How industry changes the glow.** Industry and collection change the Moon's heat in three ways, and every
     design for the array's industry has to be weighed by them:
     - *What it puts on the films.* Each extra percent of sunlight the films absorb across the fleet adds
       1.66–2.91 W/m² of glow (*screen*), so photovoltaics, contacts, coatings and trim hardware on the tiles count
       against the planned climate.
     - *How it dims.* The window rings carry the climate stack round their whole orbits and intercept 45–57 times the
       sunlight the Moon receives, so a dimmer that absorbs on them gives 15–26% of its dimming back as glow where
       the light crosses one such layer (*screen*); a dimmer that reflects gives none back. That weighs on the dimmer
       as semitransparent photovoltaics and on the added dimming the glow needs. Collectors that take light bound for
       the Moon dim it as well, and can feed industry if most of their heat leaves from faces that do not see the
       Moon.
     - *Where its heat goes.* A terawatt released in orbit puts 5×10⁻⁵ W/m² on the Moon and a terawatt used on the
       Moon 0.026 W/m², so 100 TW of industry on the Moon would add 2.6 W/m², near half the glow's low end. Heavy
       industry belongs in orbit with outward radiators.
   - **Where it belongs.** With the shield's work on its films, sized from the protection domain's optical tables
     and the ring stack's radiative transfer. It surfaced here in the heat ledger.
   - *A first-order estimate: the stack's layering and how efficiently infrared absorbed aloft warms the ground decide
     its size. The range across those is 1.8–9.9 K.*
2. **Each terawatt used on the Moon warms it by 0.017–0.026 K, so a kelvin takes 38–61 TW** (*screen*).
   - **The metropolis.** Its 200 GW is 0.005 W/m² over the Moon and 62 W/m² over its 32 km circle, near the hottest
     city cells on Earth at 100–577 W/m² (Allen et al. 2011).
   - **Decent living.** A decent life takes about 15 GJ of final energy per person a year (Millward-Hopkins et al.
     2020), so a terawatt covers about two billion people at that level. Needs stop rising near 60 GJ (Vogel et al.
     2021).
   - **What follows.** The heat budget binds industry and computing, which belong in orbit, long before it binds
     living standards.
   - *Solid.*
3. **The night is the energy problem.**
   - **The size.** Carrying the metropolis through a 354-hour night takes 71 TWh, 7.9 times all of Earth's
     pumped-hydro storage, and 52 km³ of water through Korolev's 3.74 km head (*screen*). An all-solar supply needs
     about 6,100 km² of level panels for the mean alone, about twice the city, before storage.
   - **Storage.** Grid storage classes stop near 160 hours (US DOE 2023), and hydrogen returns about 40% (Headley &
     Schoenung 2020).
   - **A grid round the Moon.** Every night-side point lies within 2,729 km of the terminator. High-voltage direct
     current loses 2–3.5% per 1,000 km plus about 1% per converter (US EIA 2018), so such a grid loses 9–19% over
     4,000–5,500 km, better than any storage round trip.
   - **Fusion.** The register's planning assumption sidesteps the night if fusion holds.
   - **Microwaves from the fleet's night half.** They would keep the night dark. The public exposure limit of
     10 W/m² (ICNIRP 2020), though, puts 200 GW on 4,000–7,000 km² of receivers closed to people and sky ships.
   - *Solid arithmetic.*
4. **Rivers are the native source; tides are not.**
   - **Tides.** Tidal power per area is a 320th of Earth's for the same range: a sixth of the gravity, one cycle a
     month. Impounding the whole nearside sea would give 25–50 GW in theory (*screen*).
   - **Rivers.** They carry 229,410 m³/s from land averaging 2,345 m above the sea, at most 0.87 TW gross (*screen*).
   - **Waves.** They bring 87–731 W per metre of coast (the wave review).
   - **Fusion fuel.** The seas and lakes hold 4.8×10¹⁴ kg of deuterium, enough for 5 TW over 10⁹ years (*screen*),
     and the operation can deliver more with any water. Deuterium–tritium fusion also needs lithium to breed its
     tritium, which the industrial architecture counts among the source industry's materials.
   - *Solid.*

### Food

5. **Food takes about 60–80 m² a person at the crop model's potential and 130–340 m² in practice.**
   - **The potential.** The day fruit gives 31–42 kcal per m² a day at the far side's equator, 60–80 m² a person at
     2,500 kcal a day (*screen*).
   - **In practice.** Allowing for the GCM's clouds (89% of clear sky), losses of 30–50% and a varied diet at 1.5–2.5
     times the land of the energy-only fruit gives 131–341 m² a person (*screen*).
   - **For the metropolis.** It needs 13,000–34,000 km² of cropland: 11–28% of the dry land within 300 km of the
     summit and 4–10% within 500 km. Within 100 km the dry land, 18,700 km², is less than the need. Lakes, Korolev
     and Hertzsprung among them, cover more of the ground within 300 km than dry land does (*screen*).
   - **Against Earth.** Earth's crops feed about 6 people per hectare after feed, biofuels and other uses, about
     1,700 m² a person (Cassidy et al. 2013). NASA's growth chambers gave potato 18.4 and wheat 11.3 g of edible
     mass per m² a day (Wheeler et al. 2008). Biosphere 2 grew 81% of eight people's diet on 0.2 ha, about 250 m² a
     person (Allen & Nelson 1999).
   - **The timing.** The day fruit ripens at sunset, once a month, so the city holds about a month of food.
   - **What follows.** The register's "regional agriculture and foraging" fits the land, given water and engineered
     crops.
   - *A planning range on a clear-sky model with no water or nutrient limits.*
6. **Protein can run on the night's power, and food must not make methane.**
   - **Protein.** Microbial protein grown on solar electricity yields over ten times the protein per area of any crop
     (Leger et al. 2021), and it can run through the night on stored or firm power.
   - **Animal products.** They use 83% of farmland for 18% of calories (Poore & Nemecek 2018).
   - **Methane.** In an air without OH, methane lives for centuries (the ecology work), so ruminants, flooded paddies,
     landfills and open anaerobic treatment would raise it unless designed out or captured.
   - *Solid in direction.*

### Building, making and moving

7. **The air ends surface launch.**
   - **The column.** A launch from the surface crosses 75 t of air per m², 7.3 times Earth's column. A mass driver's
     2.4 km/s meets about 5 MPa of dynamic pressure near the ground.
   - **Orbits.** They lie within the air below the exobase, 1,300–2,600 km up.
   - **Freight to orbit.** It goes from altitude, by tether or by elevator, and a lunar elevator of today's materials
     lifts about 584 t a year (Pearson et al. 2005).
   - **What follows.** The surface and orbital economies trade little bulk mass: orbital industry draws on space
     resources, and the surface makes its own.
   - **The film plant.** Lifting its 73–147 t/s from the surface would carry 0.2–0.4 TW in the payload's orbital
     energy alone. By S7's diagnostic each kg/s of loss budget allows about 9.4 MW deposited in the protected
     region, 9.4–940 MW across the budgets.
   - **Exhaust.** Rocket black carbon warms about 500 times as much per kilogram as soot from the ground (Ryan et al.
     2022).
   - *Solid physics.*
8. **The far side's ground gives aluminium and glass.**
   - **The ore.** Dry land is mostly feldspathic highland, since the 28% sea floods most of the iron and titanium
     basalts (the resources work). Local building at the summit runs on aluminium, glass, ceramics and calcium
     silicates, and its steel comes from distant dry mare or from space.
   - **Stocks.** A decent life needs about 42 t of materials in use per person, against 335 t in developed countries
     (Vélez-Henao & Pauliuk 2023), so the metropolis holds 4–34 Gt.
   - **Embodied energy.** Steel 24 MJ/kg, cement 4.6 and concrete 0.95 (Hammond & Jones 2008).
   - *Solid on the ground; the stocks are Earth figures.*
9. **At a sixth of gravity, wings beat buoyancy on energy, and trains brake slowly.**
   - **Air freight.** Wings and rotors cost about a sixth as much energy per tonne-km, while airships stay at Earth's.
     The sky-fleet study already found winged liners matching airships per passenger-km. Its kilometre-class
     freighters at 12–38 Wh per tonne-km sit near rail (0.24 MJ, or 67 Wh, per tonne-km) and far below air freight's
     16.4 MJ (Krammer & Schäfer 2025).
   - **Ground transit.** Standing riders take a sixth of Earth's 1–1.5 m/s², 0.17–0.25 m/s². Stopping from 30 m/s then
     takes 1.8–2.7 km against 0.3–0.45 km on Earth, and accelerating at 1 m/s² tilts the apparent vertical by 32°
     (*screen*). Metro and rail need seated or held riders, tilting cars and wide station spacing, or linear motors
     with seated travel.
   - *Solid physics.*
10. **Computing in orbit is 0.13 s away.** That is the round trip to the ring radius (*screen*), so control that
    closes in milliseconds stays on the Moon: sky-boat separation, grid protection, robots. Data centres used 415 TWh
    in 2024, about 47 GW (IEA 2025), so a terawatt of orbital computing is about 21 times today's. *Solid.*

### Water, sanitation and work

11. **Water treatment carries the whole job of disinfection.**
    - **Sunlight.** Under the shield it disinfects at 1–6% of Earth's rate (the ecology work).
    - **Energy.** Water supply takes about 4% of the world's electricity, and wastewater holds 5–10 times the energy
      its treatment needs (IEA 2016).
    - **Nutrients.** Recovering nitrogen from urine takes 1.3–13 kWh per kg against about 9.4 for the best ammonia
      plants (Larsen et al. 2021; IEA 2021). Urine and faeces could supply 22% of the world's phosphorus demand
      (Mihelcic et al. 2011).
    - **The metropolis.** It passes about 44 kt of phosphorus and 0.44 Mt of nitrogen a year through food and sewage
      (digest arithmetic, [state.md](state.md)), so its sewage is one of the Moon's main phosphorus return paths.
    - *Solid.*
12. **Work and upkeep have measured sizes.**
    - **Working hours.** A decent-living basket takes an 18-hour paid working week per person, and a "good life" 46
      hours (McElroy & O'Neill 2025). Biosphere 2's crew worked 66 hours a week.
    - **Maintenance.** Routine maintenance runs 2–4% of replacement value a year (NRC 1990).
    - **Who works.** Who performs the maintenance of a settled world belongs to the human paper. These numbers give
      it its scale. *Earth figures.*

### People across the horizons

The author's working scale of 9 October is 20–30 billion people across the Solar System by the frontier of thousands
to tens of thousands of years, affluent by modern standards, mostly on Earth, the Moon and in the array. The study is
[population.md](population.md), with its screen [population.py](population.py); numbers here are from it.

13. **20–30 billion takes growth at the pre-industrial world's pace.** From the UN's 10.2 billion in 2100, 0.075–0.12%
    a year reaches it by the year 3000 and 0.023–0.037% by 5000; the world grew 0.055% a year from year 0 to 1750. From
    10 million in 2526, the Moon reaches a billion in 2756–2987 at 1–2% a year.
14. **Affluence is 2–6 kW of primary energy a person, centred on Europe's 3.7.** Needs stop rising at 1.9 kW of final
    energy (Vogel et al. 2021), 2.9 kW of primary energy, which gives the metropolis brief's 2 kW a basis.
15. **The selected heat tolerance constrains uncompensated local energy use.** With industry in orbit a
    person at Europe's energy releases 2.8 kW on the Moon: 0.1 K of warming holds 1.4–2.2 billion people and 0.5 K
    6.8–11 billion, while 1% more dimming answers the heat of 26–41 billion. Moving industry to orbit cuts a person's
    lunar heat by a quarter. Cities at 40,000 people per km² release 60–240 W/m², so a city's heat is answered in the
    city.
16. **The runoff stress index crosses a historical threshold near 4 billion.** Its 7,240 km³ a year reaches the stress line at 4.3 billion and
    scarcity at 7.2 billion; while the seas and lakes stay fresh, water drawn from them and returned supplies beyond it.
    Food land is ample: the dry land feeds 8.4–22 billion on Earth's 12% share of land in crops.
17. **Earth binds on food within planetary boundaries.** A transformed food system feeds 10.2 billion within four
    boundaries (Gerten et al. 2020), so 16–20 billion on Earth need part of their food grown off the land, for
    instance microbial protein on power the array can send. Their heat on thermal plants adds 0.12–0.14 W/m².
18. **The array's light is no limit; its habitats are its mass, and its commute is costly.** The whole population's
    energy, 40–260 TW, is 0.005–0.03% of the night half's 851 PW. Habitats at the Stanford torus's 990 t of shield a
    person take 248–990 Gt for a billion residents. The ring radius is 14.1 hours away on a minimum-energy transfer,
    and a weekly rocket commute costs 18–46 kW a person, so the array's workers rotate and its tinkerers live there.
19. **The array serves Earth with computing and power, and prototypes Earth's own arrays.** A 10.5 km transmitter
    reaches a rectenna of the 1978 reference size on Earth from the Moon's distance, and receiving land limits it at
    20,400 km² per TW; computing serves Earth over 1.28 s of light time. The film plant, photon-kept formations,
    collectors and links are the industry that would build Earth's own arrays.
20. **Affluence lasts the horizons as a plateau.** At 2% a year energy use would reach all the light the fleet
    intercepts in about 580 years, so energy per person holds near the needs plateau and growth near zero.

### Computing in the array

The author called computing the valuable resource that turns the array's otherwise wasted energy into use, and asked
what it is for. The study is [computing.md](computing.md), with its screen [computing.py](computing.py); numbers
here are from it, with "today" the H100's 1.4×10¹² FLOP per joule.

21. **The light outruns any computing the array can build, and mass binds first.** The night half's 851 PW powers
    1,000 TW of computing on 0.4–0.6% of it, while a terawatt needs 7–34 Mt in orbit for radiators, collectors and
    processors. So 1,000 TW weighs 15–74% of the fleet and renews 0.7–2.6 Gt of processors a year, the film plant's
    scale.
22. **Heat placement binds next.** With its collectors each watt computed releases 2.9–4.4 W, and from the ring radius
    460–1,100 TW warms the Moon 0.1 K. The Earth–Moon L1 and L2 points cut that 9.5 times, 100,000 km 25 times and the
    Earth–Sun hubs 5,600 times. Radiators should run hot: at the irreversible floor, operations per square metre grow
    as T³.
23. **A terawatt today is 1.4×10²⁴ FLOP/s, 94 times the world's AI computing.** CMOS has about 200-fold left, and the
    irreversible floor at 330 K lies 4,700–470,000-fold beyond today; reversible logic goes below it by running slower.
24. **Control, management and today's science take megawatts to terawatts.** The fleet's control is 20 MW. A
    weather and climate twin of the whole Moon, run as a forecast service, takes 6–12 MW at 1 km, 6–12 GW at 100 m and
    9–15 TW at 10 m; at 1 m it takes 38–44 TW in real time and 38,000–44,000 TW as a forecast service, the one
    management task that fills the array. Following one gram of soil cell by cell takes 350 TW. A far-side radio
    correlator for 100,000 dipoles takes 1 MW.
25. **AI and digital minds span the array, and the far-out tasks pass it.** AI services for 20–30 billion people take
    14–2,100 TW today and 0.07–10 TW at the CMOS limit, and a digital population the size of ours 14–21,000 TW. An
    ancestor simulation (Bostrom 2003) runs in 8 days to 23 years on 1,000 TW. A trillion emulated minds at a
    thousand times human speed, brains emulated molecule by molecule and a Matrioshka brain pass the array.
26. **Latency divides the work.** The 0.13 s round trip to the ring radius fits inside the 0.21 s people leave between
    turns in conversation, so the Moon's interactive services and forecasts can run in the array, with millisecond
    control on the Moon and Earth's 2.56 s round trip suiting batch work.

## Statements elsewhere that need correcting

The evaluation ([state.md](state.md#statements-elsewhere-that-need-correcting)) lists stale or conflicting statements
in other folders. They belong to main and are unchanged on this branch. The largest:
- the metropolis brief quotes the canopy's gross photosynthesis as if it were food;
- its 2 kW per person has no stated basis;
- the status record's crown berths and port power are stale;
- the joint ledger gives the film plant's power in terawatts under a label in watts.

## Files

| File | Holds |
|---|---|
| [screen.py](screen.py) | The first screen: the heat ledger with the fleet's infrared, its heat mirror, computing and industry against the glow, food land, the night, tides and rivers, transit, fusion fuel, latency |
| [population.md](population.md) | Population and living across the horizons: growth, the affluent band, what binds on the Moon and Earth, the array's habitats, commuting, Earth service |
| [population.py](population.py) | Its screen, writing [results/population.json](results/population.json); sources in [population_sources.json](population_sources.json) |
| [computing.md](computing.md) | Computing in the array, supply against demand, from the fleet's control to the far-out tasks |
| [computing.py](computing.py) | Its screen, writing [results/computing.json](results/computing.json); sources in [computing_sources.json](computing_sources.json) |
| [results/first_screen.json](results/first_screen.json) | Its product (schema `terluna.provisioning.first-screen/1`) |
| [state.md](state.md) | What the project already holds for provisioning, what others assume of it, and statements to correct |
| [literature.md](literature.md) | The literature by topic, with what it means for the Open Moon |
| [questions.md](questions.md) | The open questions, in order |
| [sources.json](sources.json) | Every source cited, its use and how far it was read |

`python -m provisioning.screen` rewrites the product; `python -m pytest provisioning/tests` checks it.
