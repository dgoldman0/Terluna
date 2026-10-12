# Inhabited volume and shared population comparisons

**9 October 2026.** The Moon's surface and aerial settlements require a joint capacity investigation. The
provisioning study's earlier 2–5 billion lunar residents and 16–20 billion on Earth were illustrative placements;
they did not establish capacity or justify preferring Earth. The author reopened distribution and required the
inhabited third dimension to be evaluated. This study supplies a first reproducible screen across habitation,
provisioning, engineering and ecology, using the existing air and service products.

The author also directs that much floating habitation should roam, while managing collision and storm risks. Continuous fixed-location holding is therefore a comparison, not a universal energy demand.

The result supports substantial aerial habitation as an engineering question. It does not establish a maximum
population. Buoyant residential mass, projected area and holding power matter well before the geometric volume of
warm air is filled. Powered position-holding of large round envelopes is particularly demanding; drifting,
streamlined, moored and structurally supported districts deserve separate designs.

## What is established and what remains hypothetical

[run.py](run.py) consumes schema-checked [sky-ship](../sky_ships/results/sky_ships.json) and
[historical provisioning](../../../provisioning/results/population.json) products, records their hashes and the shared
constants, and writes [results/inhabited_volume.json](results/inhabited_volume.json). The density, pressure and lifting
allowances come specifically from the sky-ship product. Similar ecology profiles use different approximations and
must not be silently substituted. This calculation uses tabulated levels without interpolation.

Archimedean lift, spherical geometry, drag scaling, gravitational potential and mass bookkeeping are physical
relations. District masses, envelope geometry, allocation of lift, water recovery, spacing and population placements
are imposed sensitivities. No structure, full weather record, collision-safe traffic system, settlement network,
complete life-support system or inhabited ecosystem has been validated.

The older population runner and JSON remain reproducible historical illustrations. Their demographic paths, heat
factors and service estimates remain useful inputs. The interpretation of their scenarios and water thresholds is
corrected in [population.md](../../../provisioning/population.md). Full historical source-admission records remain in
[population_sources.json](../../../provisioning/population_sources.json) and [sky_ships/sources.json](../sky_ships/sources.json);
this follow-up did not perform new full-text literature admission or rerun the climate.

## One population ledger

Every new case reads [shared/scenarios/population.json](../../../shared/scenarios/population.json).
Numbers below are billions of people. They are comparisons chosen to expose requirements; none is an adopted
allocation, forecast or certified capacity. The elsewhere remainder has no demonstrated capacity credited here.

| Case | Earth | Lunar surface | Lunar aerial | Orbital | Elsewhere | Total |
|---|---:|---:|---:|---:|---:|---:|
| total20 | 10 | 3 | 3 | 2 | 2 | 20 |
| total25 | 10 | 5 | 5 | 3 | 2 | 25 |
| total30 | 10 | 7.5 | 7.5 | 3 | 2 | 30 |
| total25_earth_heavier | 14 | 3 | 3 | 3 | 2 | 25 |
| total25_moon_heavier | 8 | 6 | 6 | 3 | 2 | 25 |

A second sensitivity holds each case's total lunar population fixed and moves 0%, 25%, 50% or 75% aloft. That separates
the engineering consequences of habitation style from an arbitrary change in total people. Food demand and baseline
human heat follow the whole lunar population; they do not disappear when a person moves upward.

## Turning lift into supported floor area

For a spherical gas volume of radius 500 m, its volume is 524 million m³ and its projected area 0.785 km².
Gross supported mass is `(ambient density − gas density) × volume`. Gravity cancels from this supported-*mass*
relation, although it changes forces and structural requirements. The tabulated lifting allowance retains the
sky-ship calculation's 98%-by-volume lifting species plus 2% ambient-air convention; gas-mixture and constituent masses are recorded separately. Tables use its rounded densities.

| Height above model sea level | Pressure | Temperature | Hydrogen gross lift | Helium gross lift |
|---|---:|---:|---:|---:|
| 0 km | 1.181 atm | 20.7 °C | 669,700 t | 620,500 t |
| 10 km | 0.988 atm | 14.2 °C | 574,900 t | 532,000 t |
| 20 km | 0.813 atm | 5.1 °C | 488,500 t | 452,400 t |
| 40 km | 0.542 atm | −13.9 °C | 350,300 t | 324,100 t |

The air at 0 km is a reference state, not a proposal to place a 1-km-diameter envelope through the ground. A real
500-m-radius sphere samples a kilometre of vertical atmospheric variation; this screen uses its centre's density.
The warmer lower atmosphere has the larger lift allowance. Pressure alone does not establish respiratory suitability:
oxygen partial pressure, weather, lightning, long-term low-gravity health and rescue also matter.

This screen allocates gross supported mass as 50% residential district payload, 25% envelope and flight structure,
10% control and flight services, and **15% unused load reserve**. Those shares are assumptions, not demonstrated
mass fractions. Reserve is a capacity margin, never counted as manufactured material. A fully gas-filled envelope would be
positively buoyant at the nominal load; ideal neutral trim instead puts ambient air in a ballonet occupying 15%
of the volume and lift gas in 85%. Expelling that air and adding lift gas accesses the reserve. The full-gas
inventory below is therefore a capacity inventory, with spare gas stored externally. Operational lifting-gas mass
at nominal trim is 85% of it. Pressure management, gas supply and compartment sizing are not designed here.
The district includes its own
building, people, contents, water and local services; the overhead is the buoyant carrier. No soil-rich landscape or
large onboard farm is silently included.

The metropolis brief provides 70 m² per resident for home, work and shared services. Three hypothetical wet payload
budgets expose sensitivity to construction and inventories:

| Per resident | Light | Central | Heavy |
|---|---:|---:|---:|
| Building, floors and partitions | 10.5 t | 28 t | 59.5 t |
| People and furnishings | 2.5 t | 4 t | 5.5 t |
| Water, food and other wet stocks | 3 t | 8 t | 15 t |
| District utilities and local stocks | 4 t | 10 t | 20 t |
| Total supported district | 20 t | 50 t | 100 t |
| Building mass per floor area | 150 kg/m² | 400 kg/m² | 850 kg/m² |

These are budgets that a future design must meet, not measured lightweight city construction. The Earth-derived
42–335 t/person economy-wide stock indicators in the earlier study are a different accounting boundary. They include
infrastructure outside homes and overlap these district inventories, so they cannot be added wholesale to them.

At 10 km, the reference envelope carries 287,456 t of district payload. For the shared total25 case's five billion
aerial residents:

| District mass/person | Residents/envelope | Envelope equivalents | Sum projected area / lunar area | Imposed single-layer grid area / lunar area |
|---|---:|---:|---:|---:|
| 20 t | 14,373 | 347,880 | 0.72% | 8.25% |
| 50 t | 5,749 | 869,699 | 1.80% | 20.63% |
| 100 t | 2,875 | 1,739,398 | 3.60% | 41.27% |

The last column places centres three envelope diameters apart on a square grid. That spacing is a diagnostic, not a
safe separation standard or proposed ecological occupation limit. The results make traffic, storm avoidance and
protected sky corridors concrete research requirements. The JSON includes rounded-up fleet counts for exact service
capacity and continuous equivalents for comparison.

The central five-billion-person case carries 250 Gt of districts, permits 175 Gt of carrier structure and flight
hardware, leaves 75 Gt of load reserve, and has a full-capacity lifting-gas-mixture inventory of 48.27 Gt (41.03 Gt at ideal nominal trim). The full mixture contains about 37.30 Gt hydrogen and 10.97 Gt entrained air; at nominal trim the hydrogen share is 31.71 Gt. These are approximate constituent masses inferred from rounded product densities. Its 350,000 km² of floor area does not come
from simply assigning people to 920 million km³ of warm atmosphere. It comes from the mass that the envelopes could
support if their assumed mass shares can be achieved. Gas production, replenishment and compartment-failure design
remain substantial requirements.

At fixed mass fractions, doubling radius reduces summed projected area by half and the envelope count by eight;
lifting-gas stock stays the same for the population. This algebra does not validate larger structures. The existing
sky-ship structural study already shows that simple geometric scaling eventually fails.

## Roaming habitation and the cost of powered control

For an unstreamlined sphere, this sensitivity takes `Cd = 0.47`, air density 1.204 kg/m³ and propulsive efficiency
70%. Drag is `½ ρ Cd A v²`; required propulsive input is `drag × v / efficiency`. Here **v is speed relative to the air**.
The chosen coefficient is an idealized sphere sensitivity, not measured city aerodynamics; appendages, wakes,
propulsors, deformation, turbulence and induced-power losses require a real design.

| Relative air speed | Drag per envelope | Input per envelope | All five billion aerial residents |
|---|---:|---:|---:|
| 0 m/s | 0 | 0 in this horizontal model | 0 |
| 5 m/s | 5.56 MN | 39.7 MW | 34.5 TW |
| 10 m/s | 22.2 MN | 317 MW | 276 TW |
| 20 m/s | 88.9 MN | 2.54 GW | 2,209 TW |

For the fixed-position comparison, a sustained 5 m/s relative wind would already exceed the inherited 27.8 TW baseline for all ten billion lunar
residents in total25. These rows assume every envelope continuously experiences the stated speed; they are neither
climatological means nor a prediction of actual fleet consumption. A real fleet account needs the mean of `v³`,
not the cube of mean wind, regional wind histories and operating choices. Added operating energy would also enter
the heat ledger, with its generation losses; it is not silently included in the existing per-person baseline.

Co-moving drift provides little relative wind and little horizontal holding demand, though navigation, vertical
control, docking and avoidance remain. The roaming sensitivity below applies powered manoeuvres to a stated fraction of time, with co-moving drift during the remainder. Navigation, collision avoidance, storm escape, rendezvous and occasional repositioning need both power and reserves.

| Active relative speed | Powered time | Mean control power, five billion residents | Per person |
|---|---:|---:|---:|
| 5 m/s | 1% | 0.345 TW | 69 W |
| 5 m/s | 5% | 1.73 TW | 345 W |
| 5 m/s | 10% | 3.45 TW | 690 W |
| 10 m/s | 1% | 2.76 TW | 552 W |
| 10 m/s | 5% | 13.8 TW | 2,761 W |
| 10 m/s | 10% | 27.6 TW | 5,522 W |

These are operating sensitivities, not weather-derived safe duty cycles. The product also reports energy for one-, six- and 24-hour powered events. A six-hour event at 5 m/s costs 238 MWh per reference envelope; at 10 m/s it costs 1,905 MWh. Firm power, reserve storage, propulsor sizing and the mass to provide them must close inside the assumed hardware/wet-stock budgets. A low average duty fraction cannot substitute for the ability to manoeuvre when collisions or storms require it. Shared drifting currents may reduce relative traffic speeds, while shear, storm outflows and converging routes can increase them. Regional forecasts, separation, alternative routes and safe meeting/evacuation sites remain to be designed.

A tether can react the horizontal load without paying `drag × wind speed`
as continuous electrical stationkeeping; it needs its own load path, anchorage, mass and fatigue calculation.
A tether from the ground does not support a negatively buoyant city. Tower-supported floors require compression
structures and foundations. Streamlined buoyant districts may reduce drag, at the cost of different structural and
floor layouts. The spherical held case therefore cannot determine all aerial architectures' energy requirements.

## Water, food, heat and ground comparisons

For five billion aerial residents, domestic use at 50 or 150 L/person/day circulates 91.3 or 273.9 km³/year.
Recovery of 90% leaves 9.13 or 27.39 km³/year of makeup; 99% leaves 0.913 or 2.739. These recovery rates refer to the
whole specified domestic loop and its losses; they are assumptions, not demonstrated citywide performance.

From sea level, ideal gravitational lift is 4.49 kWh/m³ to 10 km and 8.92 to 20 km, computed by integrating lunar
`GM/r²`. At 70% pump efficiency, supplying only the makeup takes 0.67–20.0 GW at 10 km or 1.33–39.8 GW at 20 km over
these cases. These are energy lower-bound components: treatment, pipe losses, pump/pipe/carrier mass, leakage,
reserves, initial inventories, routes and emergency supply remain. No recovery credit for descending water is taken.

Agricultural evapotranspiration is outside that domestic loop. It may be much larger and depends on diet, crop
physiology, climate and water recovery at cultivation sites. Clouds cannot be counted as a freely replenished
source without their collection efficiency and effect on precipitation. The Moon's runoff divided by population is
a stress indicator, **not a hard population ceiling**. Reservoirs, fresh standing water, recycling, redistribution
and later desalination require regional budgets, and rainfall cannot be simultaneously credited to competing users.

The inherited crop requirement is 131–341 m²/person, conditional on the assumed engineered crops and diets and still
without a complete water–nutrient–chemistry model. Both aerial and ground residents are included:

| Shared case | Lunar residents | Crop area | Baseline local human heat | Added dimming screen | Runoff/person/year |
|---|---:|---:|---:|---:|---:|
| total20 | 6 billion | 0.786–2.046 million km² | 16.67 TW | 0.145–0.232% | 1,207 m³ |
| total25 | 10 billion | 1.310–3.410 million km² | 27.78 TW | 0.241–0.386% | 724 m³ |
| total30 | 15 billion | 1.965–5.115 million km² | 41.67 TW | 0.362–0.579% | 483 m³ |
| total25_earth_heavier | 6 billion | 0.786–2.046 million km² | 16.67 TW | 0.145–0.232% | 1,207 m³ |
| total25_moon_heavier | 12 billion | 1.572–4.092 million km² | 33.34 TW | 0.289–0.463% | 603 m³ |

Heat uses the inherited 2.778 kW/person local thermal-supply case with industry in orbit. Dimming uses its broad
climate-response conversion, in percentage points of design sunlight. These values exclude the film's infrared,
city-specific holding power and other additional operating loads. They do not show that actual cooling or dimming
hardware has met the requirement.

At the metropolis's 40,000 residents/km², ten billion people housed on the ground would occupy 250,000 km² of
residential districts, 1.07% of the earlier product's 23.389 million km² dry land. With five billion aloft, surface
districts occupy 125,000 km². Those footprints exclude agriculture, conservation buffers and regional services.
The opportunity is therefore architectural and ecological diversity with different engineering demands; mere lack
of ground floor space never established the old low lunar allocation.

## Sharing the sky with a productive biosphere

The aerial biome is a potential food source and a habitat deserving its own ecological account. The ecology and
resources follow-ups examine plankton, aerophytes (photosynthetic organisms held aloft by their own lifting gas), nutrient recycling and harvest. This population
screen credits **zero guaranteed food** from wild aerial harvest until productivity and nutrient closure are
established. The conversation's 1% projected-cover, 1,000 g C/m²/year and 5%-harvest case gives about 185 million gross
person-energy equivalents; an assumed edible fraction 0.5 and processing retention 0.8 reduce that to about 74 million.
Those are sensitivity arithmetic, not a complete diet, safe harvest rate or demonstrated annual supply.

Multiple aerial levels share the same incident sunlight. Their projected areas cannot be multiplied into independent
solar budgets. Summed envelope area is neither a mapped shadow nor an ecological occupancy measure: overlap,
transparency, Sun angle, diffuse scattering and migration all matter. Habitats, crops and wild aerophytes
must be assessed together for shading, cloud changes, phosphorus return, artificial light, storm and migration corridors.
Neither a settlement altitude band nor a permitted sky-cover fraction is adopted here.

Roaming is the working direction for much floating habitation. The next investigation should compare several actual district designs: a drifting buoyant settlement with protected
wet stocks; a streamlined or moored buoyant settlement with a regional wind history; and elevated districts supported
by towers or terrain. Each needs load paths, full wet mass, gas losses, energy through the lunar night, accessibility,
traffic and rescue, food and water routes, nutrient return and ecological occupancy. Those designs can progressively
constrain population distributions while retaining the shared comparison ledger.

## Reproduction and checks

```sh
python -m research.studies.inhabited_volume.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest research/studies/inhabited_volume/test_run.py
make check
```

The focused tests check product schemas and hashes, shared constants, population and mass conservation, independent
mass-per-area calculation, drag force/power scaling, zero-drift and perfect-recovery limits, altitude effects and a
numerical gravity integral independent of the closed form. This is arithmetic validation; ecological viability and
engineering certification remain open. The [check record](checks.json) records 24 passing focused tests (15 new study tests), zero layer violations, calendar 13 pass/1 skip, and passing restored-provenance and ensemble-integrity checks. Required `make -k check` returned 2: its Python collection ran before optical-input restoration, and the captured run ended during the immersion sky bake. The separate remaining targets passed afterward; no full repository-suite pass is claimed. The unchanged first-screen provenance test also requires the absent atlas NPZ.
