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
   - **The flow.** The fleet's tiles absorb 65–83 PW of sunlight and stay facing the Moon. The face turned to the
     Moon carries half of their emission over the two halves of each orbit. A tile 19,355–21,392 km out sends the
     Moon 0.66–0.81% of what that face emits. The Moon receives 5.7–8.9 W/m² averaged over its surface, or 2.8–4.4 W/m²
     counting only one layer where tiles hide one another, since the stack shows the Moon about two layers (*screen*).
   - **In proportion.** That is 1–3% of the design sunlight, the heat of 106–336 TW used on the Moon.
   - **The warming.** If it acts like sunlight, the GCM's 1.9 K per 1% of sunlight puts the design case 1.8–5.7 K
     warmer (*screen*). The design case of 294.9 K leaves it out. The 5% dimmer, chosen to bring the Moon from 303 to
     294 K, removed 15.4 W/m² of sunlight at the top of the air.
   - **The levers.** A low-emissivity coating on the Moon-facing face, transparent in the visible as low-E window
     coatings are, would send most of the tiles' heat outward. More dimming is the other lever.
   - **Where it belongs.** With the shield and climate work. It surfaced here in the heat ledger.
   - *A first-order estimate: the stack's layering and how efficiently infrared absorbed aloft warms the ground decide
     its size. The range across those is 0.9–8.8 K.*
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
| [screen.py](screen.py) | The first screen: the heat ledger with the fleet's infrared, food land, the night, tides and rivers, transit, fusion fuel, latency |
| [results/first_screen.json](results/first_screen.json) | Its product (schema `terluna.provisioning.first-screen/1`) |
| [state.md](state.md) | What the project already holds for provisioning, what others assume of it, and statements to correct |
| [literature.md](literature.md) | The literature by topic, with what it means for the Open Moon |
| [questions.md](questions.md) | The open questions, in order |
| [sources.json](sources.json) | Every source cited, its use and how far it was read |

`python -m provisioning.screen` rewrites the product; `python -m pytest provisioning/tests` checks it.
