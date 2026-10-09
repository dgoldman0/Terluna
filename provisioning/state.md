# Where provisioning starts: what the project already holds

An evaluation of the repository at `main` `996ea6c` (9 October 2026) for provisioning. Paths are from the repository
root. The full working digest, with line references, is kept with the domain's run records on the research drive
(`research/runs/provisioning/repo_digest.md`).

## Decisions that bind provisioning

From [research/decisions.md](../research/decisions.md) unless named:
- **The cycle and the air.** The Moon keeps its 29.53-day cycle, and life and culture store resources across the
  night. The atmosphere is open, at 1.2 atm, with no dome.
- **Lifetimes.** Construction within 500 years, and operation for at least 10⁹.
- **Volatiles and water.** Water, nitrogen and oxygen come from the Solar-System-wide operation. Standing water
  covers 28%, with rain-fed lakes on 10.3% more.
- **Protection and its power.**
  - The protection architecture keeps an industrial hub and Earth–Sun L1/L2 hubs for freight, power and industry.
  - The shield system gives reasonable net-positive electricity for habitat. The 23.835 TW and 100 TW cases are
    comparison references.
  - The ring fleet is the lead design. Thermodynamic computing is kept for the array's control and services.
  - One later branch, `research/array-industry`, takes heat reuse, computing and industry with the trim and store
    power.
  - The fleet keeps its light off the night zones and puts the turned light to use where it can.
- **Traffic.** The six traffic-safety rules govern all bulk freight to the Moon.
- **The summit.**
  - The port is commercial, with about 1.5 million people and no homes or industry.
  - The metropolis holds about 100 million people at 40,000 per km² within 32 km.
  - Sky boats, gliders, a metro, rail and sky ferries carry it.
  - Fusion supplies most of the Moon's local power, supplemented by regional solar and wind (a planning assumption).
  - Regional agriculture and foraging feed the metropolis (a planning assumption; "the land per person has not been
    computed").
- **The flyers.** The findings of 2026-09-28 on the sky ships and flyers: the rule of similarity, rigid hydrogen ships
  to about 3 km, wings and rotors cheap, sky boats, ferries, the height stack, winged liners.
- **Climate.** Laptop-only compute, and the climate work paused with comfort carried as a range. The investigated climate is the
  planned climate (author, 2026-10-09): what the shield and the array add to the Moon's heat is answered so that the
  Moon keeps it.

The author's statement that this is a sky-ship world where lift is cheap had no row in the register. This branch adds
it to the register's section on sky ships and flyers, in the author's words.

## What is computed

| Result | Numbers | Evidence | Where |
|---|---|---|---|
| Metropolis demand | 200 GW: 100 million at 2 kW, 80 W/m² over its land, 62 W/m² over the 32 km circle | Planning figure; the 2 kW has no stated basis | `habitation/summit_metropolis`, the joint ledger |
| Port power | 411 MW at 100 kWh/m² a year; level panels 32.8 W/m² clear-sky mean; Korolev pumped storage 1.35 kWh/m³ over a 3.74 km head | Screening | `research/studies/summit_tower` |
| Wind and waves | Wind 16 W/m² at the summit's ground, 195 W/m² at the crown; screened rotors 287–927 MW; waves 87–731 W/m | Screening and SWAN runs | `summit_tower`, the infrastructure review, `climate/waves` |
| Power from the array | The aperture intercepts about 240 PW. Collectors take 2,450–3,670 km² per TW. Microwave links of 100 MW over 100 km run at 51% from bus to bus. Power delivered to habitat and industry is not yet modelled | Arithmetic and bounded calculation | `solar_shield_array`: integrated comparison, electromagnetic, array industry |
| Heat | A terawatt used on the Moon adds 0.026 W/m²; released in orbit, 0.19% of it reaches the Moon; the tiles absorb 65–83 PW | Closed-form screen | `research/studies/solar_shield_array/results/array_heat.json` |
| Food | Canopy 12.5 g C m⁻² per 24 h gross; whole plant 228–266 g C per cycle; day fruit 1.9–4.9 kg m⁻² per cycle | Clear-sky models | `biosphere/canopy` |
| Transport | Sky boat 0.74 kWh per 10 km; ferries 11.4 Wh per passenger-km; liners 12–27; winged liners 14; freighters 12–38 Wh per tonne-km; high platforms 115–750 kW | First-order sizing | `research/studies/sky_fleet` |
| Computing | About 1.4×10²⁴ operations a second per TW; 830 km² of radiator per TW at 330 K | Screen | `array_industry.md` |
| Materials | Fleet 46 Gt with film renewed at 2.3–4.6 Gt a year; the tower's chosen form 98 Mt of steel; regional magnets 1–4 Gt | Proposed designs | the shield, tower and magnet studies |
| Civilization scale | About 501 PW (K ≈ 1.17), three separate power ledgers, a bootstrap from fission to fusion and near-Sun solar | Recovered design history | `engineering/reference/industrial_architecture` |
| Drifting days and routes | Median drifting day on the equator 21.8 d at 10 km and 9.3 d at 40 km; about 34 d near the ground, with a long evening; the hour held within 1 m/s in 18–31% of 3-day means in the lowest 5 km at 40–75°; passive routes at 2.5 km share one track within 180 days; shear 4–6 m/s over 15 km | GCM 3-day means at T21; passive routes at fixed heights | `climate/gcm/zonal_winds.md` |
| Ways of living | Eighteen settings with their costs per billion; towns 0.6–10% of the dry land across the shared cases; ten billion lunar residents crop 35–224% of the rain-fed land outside a kept 30%; heat 0.08–1.8% more dimming; five billion aloft in about 870,000 towns with 48 Gt of lifting gas | Screens on the atlas, climate and light products | `research/studies/ways_of_living` |
| Daily life | Light states per setting through the cycle; home lighting 11–45 kWh a person a night; drifting days from the design winds; Earth 2.4–2.7 s and 2.1–5.0 days away | Clear-sky light screen and closed-form arithmetic | `research/studies/daily_life` |

No product gives the Moon's population, its demand by sector, its power system as a whole, power delivered from
orbit, its food system, water and sanitation, waste, its industrial base, or the economy and work.

## What other lines assume of provisioning

- **The metropolis and port:** fusion with solar, pumped storage, a 2 kW person, regional farms, interiors kept
  comfortable, snowmaking, and "which materials can be made on the Moon".
- **The flyers:** hydrogen as lift gas and fuel cells, batteries at 235 Wh/kg, 1.1 million shared sky boats with a
  traffic system, and yards for kilometre-long ships.
- **The shield and array:**
  - the film plant and spare tiles, collectors and power links, and computing;
  - trim devices of common elements;
  - habitats with their life support;
  - freight for the film plant, and the L1/L2 hubs.
- **The magnets:** 9–15 GW of refrigeration every night.
- **Conservation:** cryogenic sample stores, the undersea dome, and "services that inhabitants can build, repair and
  run themselves".
- **The ecology register:** salt works, dredging, settlements' wastes to land, and harvests from the sea.
- **The CO2 study:** a lime-kiln return of CO2.
- **The manuscripts:**
  - the human seed asks who performs maintenance;
  - the engineering seed asks for heat placement and buffered arrival;
  - the locked plan asks for a common system ledger.

## Open items already listed

- **Engineering:** source-to-use pathways, three power ledgers, the habitat power network to the delivered 100 TW
  case, and lightning protection.
- **Habitation:** the port's power, hydrogen or helium for the ships, the night, power delivered to the surface,
  transport and launch between the surface and orbit, and what launches release into the upper air.
- **The array industry:** trim and store hardware, the film plant, power and heat ledgers by place, computing, and
  which industries go to orbit.
- **The summit:** "A power plan for the summit city" and the land per person for food.

## Statements elsewhere that need correcting

These belong to main and stay unchanged on this branch.
1. **Gross photosynthesis quoted as food.** The metropolis brief uses the canopy's gross 12.5 g C per m² a day. Food
   comes from the whole plant's 7.7–9.0 g C, or the day fruit's 3.1–3.7 kg per m² a cycle.
2. **The 2 kW person** in the metropolis brief has no stated basis. Its "10–20 million people" from local sunlight
   does not state the share of land in panels.
3. **"80 W/m² over the city"** in the joint synthesis counts the city's land only; over the whole 32 km circle the
   figure is 62 W/m².
4. **Stale figures in the status record and studies:**
   - the status record's crown berths (27 berths with 772 m arms and 8 with 424 m) and port power (about 400 MW
     mixing September and corrected figures);
   - the infrastructure review's hydrogen at the crown (about 2,100 t against 2,300–4,100) and night water (105
     against 107.9 million m³).
5. **Storms.** The summit study says ships dock "so they never descend through the storms", yet the flight band lies
   inside the strongest storms. The metropolis brief's "the backbone carries on" needs the storm rules the hydrogen
   ferries carry.
6. **A unit error.** The joint ledger's `film_plant_w` holds terawatts (fixed on main on 9 October, 78020e0).
7. **Lunar oxygen.** The engineering seed and the locked plan-01 outline still name lunar oxygen routes, which D4
   replaced.
8. **Collectors on the tiles.** The integrated comparison puts them "on its tiles or flying free", against the
   one-sail-loading rule. Spread over a ten-thousandth of each tile they are harmless.
