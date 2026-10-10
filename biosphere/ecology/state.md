# Where ecology starts: what the project already holds

An evaluation of the repository at `main` `996ea6c` (9 October 2026) for ecology. Paths are from the repository root.
The full working digest, with line references, is kept with the domain's run records on the research drive
(`research/runs/ecology/repo_digest.md`).

## Decisions that bind ecology

From [research/decisions.md](../../research/decisions.md) unless named:
- The Moon keeps its 29.53-day cycle, and life adapts to it, storing resources across the night.
- An open atmosphere with no dome: 1.2 atm, O2 at Earth's 21 kPa as 17.5% by mole, and the shield's titania stack
  removing the ultraviolet below about 340 nm.
- Standing water covers 28% of the Moon, and the corrected climate's drier land is accepted: rain-fed lakes on 10.3%
  of the Moon, and 42% of the land under 0.5 mm of rain a day.
- The light rules:
  - practical dusk lasts 82 hours at the equator;
  - nearside nights keep 2.3–2.9 lux of Earthlight;
  - only the far side's equatorial seas grow dark at midnight.
- The ring fleet avoids lighting the Moon's night zones. The requirement's value "comes from the biosphere and
  human work".
- Climate work is paused, and comfort is carried as the range between CM1 and the GCM.
- Life and people: aerial and gliding life is a major feature; high-altitude life is sparse; the seas are living,
  with design guesses for their contents.
- Regional agriculture and foraging feed the metropolis; the land per person had not been computed.
- Conservation keeps two programmes (D1–D11 in [research/studies/conservation](../../research/studies/conservation/README.md)).
- The biosphere paper owns biology and ecology together ([AGENTS.md](../../AGENTS.md)).
- On 9 October 2026 the author asked for giant aerophytes that grow over centuries as far as the physics allows, for
  their water to be gathered from the air, and for a new name and a taxonomic system for the biosphere
  ([taxonomy.md](taxonomy.md)).

The author's standing directions on engineered organisms, developed soils, a mix of climates and imported elements
had no row in the register. This branch adds the first two to the register's section "Life and people", in the
author's words.

## What is computed

| Result | Numbers | Evidence | Where |
|---|---|---|---|
| Canopy photosynthesis | 12.5 g C m⁻² per 24 h at the equator, 18% above Earth's from 35% fewer photons; the diffuse sky adds 32% | Clear sky, hypothetical traits | [canopy/](../canopy/README.md), `canopy/results/canopy.json` |
| Whole plant through the cycle | 228–266 g C m⁻² per cycle; night store 47 g C, or 10 when idling; 238 dark hours at the equator and 34 at 60° | Clear sky, GCM temperatures, carbon scenarios | `canopy/results/plant.json` |
| Day fruit | 1.9–4.9 kg m⁻² per cycle at 30–80% fruit share; 50% above the same stand on Earth | Design requirement | `canopy/results/fruit.json` |
| Carbon reserve theorem | Minimum reserve 2.5–15.8 reference days | Proved under its model | [long_night.py](../long_night.py), [findings](../../research/findings.md) |
| The ecology register | Selections across the clocks, storage, fruit, landscapes, flight and the rock cycle's return limb | Ideas register with screening arithmetic | [research/studies/lunar_cycle_ecology](../../research/studies/lunar_cycle_ecology/README.md) |
| Forest mechanics | 300 m trees hold to 23.5 m/s at the roots; a 49-tree patch cascades at 22 m/s | Static screens, LAI 22 | `megaforest_wind`, `forest_patch` studies |
| Waters | Euphotic depth 68 m in the open sea, 9–12 m at productive coasts; design guesses for chlorophyll, CDOM and fines | Computed optics on Earth seawater | [living_water/](../living_water/), `sea_appearance` study |
| Tides | 3.7 m typical monthly range on the nearside sea, 2.1 m in South Pole–Aitken, 0.55 m in Smythii | Computed | [geography](../../geography/README.md#the-monthly-tide) |
| Nitrogen | Lightning gives 0.0043–0.0079 flashes per km² a year and at most 0.2% of the biomes' need | Computed with Earth yields | [research/studies/joint_synthesis](../../research/studies/joint_synthesis/README.md) |
| The air's chemistry | No OH, ozone or HO2 near the ground; lightning's NO held at 1–10 ppb by assumed deposition | Global-mean bare column | `joint_synthesis/results/air_chemistry.json` |
| CO2 | Plants need about 28–30 Pa for diverse ecosystems; the design's 49 Pa is comfortable; weathering draws it down in centuries to millennia | Literature synthesis | [research/studies/atmospheric_co2](../../research/studies/atmospheric_co2/README.md) |
| Light | 89,600 lux with the Sun overhead; UV index 0.11 against 14.8; twilight light per Sun depression; Earthlight by coast; the fleet's glow | Clear-sky solver | [illumination/surface_light](../../illumination/surface_light/README.md), `sea_appearance` |
| Zones | Seas, lakes, wet land, fog desert and polar dry land with their weather (the aerosol study's regions) | Structured estimate | [research/studies/open_moon_aerosol](../../research/studies/open_moon_aerosol/README.md) |
| Taxonomic register | 36 designed groups with 11 subgroups, 28 design modules and 12 communities on a 25-node Earth backbone; six candidate scientific names for the aerophytes | Design register with source bindings; cycle and size classes computed | [taxonomy.md](taxonomy.md), `results/taxa.json` |
| Giant aerophytes | Round bodies 50–700 m (cellulose-class tendons) and 1 km (300 MPa fibre) at 10 km, 1.5–3 km in canopy-like light; colonies of 60–150 m modules without a size limit; 2-km sky reefs in 36–544 years, 1-km round giants in 73–1,186; rain closes equatorward of 35–38° (51° with a CO₂-concentrating trait) | Screen values; design guesses stated as ranges; three independent numerical reviews | [aerophytes](../../research/studies/aerophytes/README.md), `results/aerophytes.json` |

No biome map, soil model, water model of plants, trace-gas budget of the biosphere, animal physiology or aquatic
food web exists yet.

## What other lines assume of ecology

- **The air's chemistry** assumes soils and plants remove H2, CH4, CO, organics and NO, and holds N2O fixed at
  330 ppb.
- **The joint synthesis** says biological fixation supplies the nitrogen.
- **The fleet decision** waits for the night requirement's value.
- **The aerosol study** takes dimethyl sulfide, isoprene, iodine, spores and crusts from the ecology register's
  landscapes.
- **The climate models** carry a placeholder land surface: albedo 0.2, moisture availability 0.5, 10 cm roughness.
- **The sea optics** take the living seas' contents as design guesses.
- **The CO2 study's box model** needs the biomes' carbon stocks.
- **The metropolis** needs its food land.
- **The wind devices** read flier speeds from the ecology register.
- **Requirement O3**, the shield's ultraviolet cut-off, is still open, and no biological criterion stands beside it.

## Open items already listed

- The biosphere's status lists:
  - leaves surviving about 240 dark hours, and the day fruit's untested programme;
  - nutrient and energy budgets, and the rock cycle's return limb;
  - biological fixation;
  - the light dark-night organisms need.
- The plan lists:
  - item 1, the night requirement;
  - item 3, the air's chemistry with the biomes' and cities' emissions;
  - item 4, nitrogen and the dark night.
- The ecology register's next work is the author's review, physiology tests, an element budget and CO2 under still
  night air.

## Statements elsewhere that need correcting

These belong to main and stay unchanged on this branch.
1. **The ecology register's numbers predate the corrected design case.** Its product's input hashes no longer match
   the plant, fruit and drainage products: lakes 11.6% against 10.3%, rivers 279,000 against 229,410 m³/s, the land
   at 294.0 K against 294.9 K, and several plant figures. Its "coldest hour anywhere is 18.7 °C" meets the GCM's
   coldest 3-day mean of 14.3 °C.
2. **The aerosol study's oxidation chemistry.** Its sulfate and secondary organic particles assumed dimethyl sulfide
   and isoprene oxidizing, which the joint synthesis's bare column, three days later, found no oxidant for. Its fire
   ignition (0.003–0.006 ground strikes per km² a year) is two to four times the Moon-wide rate.
3. **The forests' leaf area.** The megaforest and forest-patch studies use LAI 22 in air at 288 K. The canopy model
   puts the leaves' optimum near 6, and the climate is 294.9 K.
4. **The waters carry raw regolith.** The living-water fines are raw Apollo soils, against the developed-soils
   direction. The water optics take Earth seawater at 38.4 g/kg while the seas' composition is open, and the CDOM
   ranges assume Earth's ultraviolet.
5. **The joint synthesis's nearside Earthlight.** Its night step adds a flat 2.3 lux to every nearside place, so
   Procellarum by Russell reads 2.4 lux where its own light calendar gives 0.48 at the darkest, and the Smythii
   headland 0.001 against 0.085. Fixed on main on 9 October (78020e0): each place now takes the Earth's height and
   phase over it, 0.44 lux by Russell and 0.07 at the Smythii headland.
6. **The CO2 study's inputs.** It uses 26 °C and 349 mm of runoff against 21.75 °C and about 265 mm now, and the
   build-up's 72% land includes the lakes.
7. **Plant twilight** still comes from the exponential-column sky atlas, which the illumination status asks to
   replace with the solved-column sky.
8. **Fog gardens.** The locked plan's "highland fog gardens" and the root seed's clouded highlands meet a highland box
   that is dry and sunny.
9. **"High-altitude life is sparse"** rests partly on strong ultraviolet aloft, which the shield removes; radiation
   aloft now means cosmic rays and the light through gaps.
