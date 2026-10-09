# The Open Moon after the joint integration

On 8 October 2026 the five lines of work came together on one tree: main with its corrected climate, waves, tides
and optics; the sea-appearance study; the solar shield and habitat array; atmospheric electricity; and the
infrastructure of the summit port. The [integration record](../../integration/README.md) gives the merges. This
study joins the work itself. Each line had modelled its part of the Moon against assumptions about the others. Here
those assumptions meet the other lines' final answers, and eight joint calculations measure where the pieces touch.
The calculations ran in about an hour of machine time in all. Every number below comes from a product named beside
it, and each product states its evidence and reading rule.

The author asked for this as the second half of the integration: "taking the bits and pieces and seeing where they
fit together and inform on the whole and what comes next."

## In brief

- **The ring fleet lights the Moon's night.** The lead shield keeps its tiles facing the Moon round their whole
  orbits, so the night half of the fleet shows the Moon the face the Sun lights. At the 0.1% diffuse reflectance that
  requirement O5 sets for Sun-facing surfaces, the glow alone puts 25–66 lux on the night's ground, against 0.001 lux
  at a far-side equatorial midnight and about 5 lux under the full Earth over the nearside. The Moon would have no dark night
  anywhere. Tiles tilted 1.5–6° in the keeping envelope would add mirror glints averaging hundreds to over a
  thousand lux. Rolling the night half's tiles by 80°, which the shield's clearance study allows, keeps every glint
  off the Moon and brings the glow to about 2 lux at 0.1% scatter and 0.02 lux at 0.001%. The author's decision is
  that the fleet avoids lighting the night zones.
- **The air under the shield has no oxidant.** The titania stack removes the ultraviolet that makes OH, ozone and
  HO2, so the design column holds none of them. Gases released into the air (hydrogen, methane, carbon monoxide,
  organic compounds and lightning's NO) leave only through soils, rain and escape. Lightning's NO accumulates in the
  model without bound; with deposition, NO near 1–10 ppb is the bracket. It leaves the shield's loss screen
  unchanged: the upper air is too cold for NO to radiate.
- **Lightning over the whole Moon** averages 0.004–0.008 flashes per km² a year, about 0.01 a second, a third of the
  electrified box's rate and about 1/350 of Earth's per area. The summit port's cells flash at 0.011–0.083, up to four
  times the box. Lightning fixes 0.007–1.6 kg of nitrogen per km² a year, at most 0.2% of what productive ecosystems
  gain, so biological fixation supplies the biosphere's nitrogen.
- **The port's crown starts its own lightning.** In the box's storm climate the storm's potential at a 24 km crown's
  height passes the upward leader's inception potential about a tenth of the time, in some 100–120 separate episodes
  a year at any one place, and most of the time while storms are near. Potentials reach tens to hundreds of MV.
- **Haze barely touches the light.** The Open Moon's air already scatters half the noon sunlight (molecular optical
  depth 0.71 at 550 nm). The aerosol study's ordinary hazes add 0.002–0.035, keep visibility at 145–240 km and move
  practical dusk by under an hour of its 82; only its loaded cases cut visibility to 9–48 km.
- **The regional magnets reach the south pole's heritage.** Turned any way, the four 500 km loops' 0.5 mT lines take
  in four or five landing and impact sites and 13–14 science targets near the south pole; the 1,000 km design's
  narrower 62 km lines take in two to five sites.
- **The sky fleet's hydrogen stays a small part of the air's.** A fleet losing 1–10% of its lift gas a year releases
  480–7,100 t, six to ninety times what the air loses to space, and each megatonne held in the air raises its
  hydrogen by 0.005 ppm.

## Where the pieces meet

### 1. The night with the ring fleet

[illumination/fleet_light](../../../illumination/fleet_light/model.py) builds the lead fleet from the shield's own
products: 1,859 rings over the 8,625 km stack, each ring's radius by height and its forced eccentricity (0.11–0.16,
apolune toward the Sun, so the night half rides closer), the shared node line at right angles to the Sun, and the
window and annulus films' visible optics (2.1% and 1.3% visible reflectance, 97.8% and 98.7% transmission). Each ring
is a continuous strip of tiles 10 km wide. Sunlight reaches a night-side tile through the layers in front of it and the
Moon's penumbra. Results are in [night_light.json](../../../illumination/fleet_light/results/night_light.json) and,
against the night's own light from the sea-appearance light calendar, [results/night_sky.json](results/night_sky.json).

The night half of the fleet presents about 1.1 billion km² of sunlit tile to the Moon and intercepts about 8×10¹⁹
lumens, some 0.85 EW of sunlight, roughly a hundred times what the Moon itself absorbs.

| Light at local midnight, equator | Far side | Nearside, under the Earth |
|---|---|---|
| The night's own (twilight's last light, Earthlight) | 0.001 lux | 5.5 lux |
| Fleet, diffuse scatter at 0.1% (O5's value) | 66 lux | 66 lux |
| Fleet, diffuse at 0.001% | 0.66 lux | 0.66 lux |
| Fleet, night half rolled 80°, diffuse at 0.1% | 2.2 lux | 2.2 lux |

The glow falls from 66 lux at midnight to 25 lux near the terminator, and from the equator to 16 lux at 75°
latitude. At 0.1% scatter the far-side equator's 119 hours below a tenth of a lux and 190 hours below practical dusk
all disappear. Holding the fleet at a tenth of the night's own midnight light takes a Moon-facing scatter of about
10⁻⁹ at the far-side equator, 10⁻⁷ at 30°, 10⁻⁵ under the nearside's Earth and 10⁻⁴ at 60°. The night's own light
takes the Earth's height and phase over each place, in the mean geometry without libration: on the central meridian
the full Earth stands overhead at the equator's midnight, and those nights stay above practical dusk throughout,
4.1–4.9 lux at their darkest.

Mirror reflection depends on attitude. A tile facing exactly the Moon's centre returns light at the distance from
the Moon it came in, so it misses; only tiles in the penumbra land theirs, 100–150 lux 60–80° from midnight, inside
the evening twilight. Tilts turn the reflected ray through twice their angle:

| Night-side tilt, random direction | Share of mirrored light landing | Mean over the night hemisphere | Midnight band |
|---|---|---|---|
| 1.5° | 0.9% | 570 lux | 0–74 lux near midnight, 1,200–1,300 lux at 20–40° from it |
| 3° | 1.7% | 1,100 lux | about 5,900 lux |
| 6° | 2.2% | 1,400 lux | about 3,700 lux |

Each tile's glint is an image of the Sun about 170 km across, at 7–12 lux. Only tiles within a few thousand
kilometres of the Moon's shadow can land theirs: 3.5% of the night-side tiles for 1.5° tilts, 14% for 6°.

**The design options.** The author's decision ([register](../../decisions.md#solar-shield-and-habitat-array), 8 October)
is to avoid lighting the night zones: the light is redirected away from the Moon and used where it can be, for
computing for instance. Lighting over cities might be acceptable if fairly fine spotlights are possible; a tile's
glint, an image of the Sun about 170 km across, covers a city and its whole surroundings.
- *Steering the glints away* works with a rule for the tiles near the shadow, tilting only away from the Moon there,
  or with the roll below, which sends every reflection far from it.
- *Rolling the night half.* The shield's attitude study found rolls of 30°, 60° and 80° outside the service arc
  within the bundle's clearance, where edge-on turns collide over 82% of each orbit. Rolled 80°, the Moon-facing
  faces meet the Sun and the Moon obliquely: the glow falls thirtyfold, to about 2 lux at 0.1% scatter, and no
  mirror light lands. The cost is the shield's own: the push on the eccentricity falls to 0.48, turning a tile takes
  0.2–2 MN·m, and the momentum store grows to 2–6×10¹⁰ N·m·s.
- *Conditionally reflective tiles.* The trim films the shield already sizes (liquid-crystal films of the kind IKAROS
  flew, electrochromic films switched 2–20 times an orbit) can hold a night state of least reflection. That lowers
  the mirrored light; the diffuse scatter is a property of the films' surfaces and their wear, so it stays.
- *Using the light.* The night half's 0.85 EW never reaches the Moon. The light a rolled or switched tile turns away
  can feed collectors for computing, which the planned research/array-industry branch takes up.

Together, a roll near 80° and a Moon-facing scatter of 10⁻⁵ give about 0.02 lux at a far-side midnight, below a
tenth of a lux but twenty times the natural darkness there. A natural far-side night needs scatter near 10⁻⁷ even
when rolled. O5 covers Earth's night side; the Moon's own night needs a requirement of its own. Its value belongs to the biosphere and human work (the lunar-cycle ecology's dark nights,
circadian life, far-side astronomy).

The model stands the Moon's orbit plane in for its equator (they differ by up to 6.7°), takes tiles spaced evenly
along each ring (spaced evenly in time they would thin by up to a fifth on the night half), treats the films as
Lambertian scatterers and leaves out the Moon's atmosphere at the limb, which dims the penumbral tiles' reflections.

### 2. The air under the shield

The middle-atmosphere model's design column (1.2 atm, titania stack, 288 K) holds ozone at 10⁻¹¹ DU, OH near 10⁻²⁷
and HO2 near 10⁻²¹ as mole fractions ([air_chemistry.json](results/air_chemistry.json)). On Earth OH, made from
ozone's photolysis, is the air's cleanser. Under the titania stack the air keeps hydrogen at its 0.53 ppm and would
keep any methane, carbon monoxide or organic gas that biology, fires or cities release, until soils, rain or escape
take it. The column's own surface ultraviolet is an index of 0.12.

[air_chemistry.py](air_chemistry.py) gives the column lightning's NO, spread through the troposphere (whose top lies
at 86 km) from the electrified box's per-joule yield carried over the Moon: 2.2×10⁸ molecules per cm² a second, a
sixth of Earth's mean. With no oxidant and no sink for NO in the model, the NO grows until the solver stops
(139 ppm) and ozone and the ultraviolet stay as they were. The model lacks NO's slow deposition to soils and plants
and the lunar night's chemistry. Deposition at 0.001–0.01 cm/s would hold NO near 1–10 ppb at the central yield;
NO + NO + O2 alone settles it near 0.7 ppb. That NO, its slow conversion and its uptake by plants are how
lightning's nitrogen would reach the land.

[upper_air_no.py](upper_air_no.py) holds NO at the loss response's 0.3 Pa base at 1, 10 and 100 ppb in place of the
case's 10⁻¹⁶ and reruns the shield's thermal-column sweep for the design shield. The exobase temperature (145–210 K)
and the molecular loss stay the same to seven figures: at these temperatures NO's 5.3 µm band barely radiates, and
CO2 at 400 ppm cools the column ([upper_air_no.json](results/upper_air_no.json)). The shield's loss screen holds.

### 3. Lightning over the whole Moon

[atmosphere/electricity/occurrence.py](../../../atmosphere/electricity/occurrence.py) takes the box's 530 flashes in
two lunar days (0.022 per km² a year) per unit of the GCM's convective rain at the box's own two T21 cells, and
carries that rate over the Moon by the corrected run's ten settled years of convective rain
([climate/gcm/convection.py](../../../climate/gcm/convection.py); 2.0 mm a day over the Moon, 1.8 on land, 2.3 over
the seas). Three scalings (linear in the 3-day mean, its 1.5 power, its square) and a sea factor (sea storms at the
land rate or a tenth of it) bracket the translation.

| | Linear, seas at the land rate | Bracket |
|---|---|---|
| Moon mean | 0.0079 per km² a year | 0.0043–0.0079 |
| The whole Moon | 0.009 flashes a second, 0.0019 ground strikes | 0.005–0.009 |
| Summit port | 0.046 | 0.011–0.083 |
| Polar plain, 75° N | 0.0001 | 0–0.0001 |

The box sits on a convective spot of the GCM (6 mm of rain a day), so the Moon averages a third of its rate. The
summit's cells rain more and flash more. Earth's lightning averages about 2.7 flashes per km² a year.

Nitrogen ([ledgers.json](results/ledgers.json)): the box's three routes to fixed nitrogen, carried over the Moon, give
0.026 (per flash), 1.6 (per joule) and 0.007 (per metre of channel) kg N per km² a year, or 1, 62 and 0.25 kt over
the Moon. Productive ecosystems on Earth gain of order 1,000–5,000 kg N per km² a year to replace their losses, the
bracket taken here; the canopy model's equatorial stand fixes 12.5 g of carbon per m² a day. Lightning supplies at most 0.2% of that nitrogen. The
biomes' nitrogen comes from biological fixation, with the delivered nitrogen of the resource operation behind it.

### 4. The crown in the storm

[atmosphere/electricity/tower.py](../../../atmosphere/electricity/tower.py) solves the electrified CM1's field
equations again for each of the 905 saved charge fields of the canonical run and its two ten-minute windows (its
numpy solver matches the direct solve of the same difference equations), and reads the potential at a tower's top
height over every column. A grounded tower meets the air at its top with that potential. The thresholds are the
lunar storms' adopted leader crossing (0.225 MV positive, 0.4 MV negative) and Rizk's potential for stable upward
leaders from tall structures, 1.0–1.2 MV at the crown's air density
([tower_exposure.json](../../../atmosphere/electricity/results/tower_exposure.json)).

| Crown height above the box's flat ground | Share of time above Rizk's potential | Above 0.4 MV | Episodes a year at one place | 99.9th percentile |
|---|---|---|---|---|
| 24 km (the tower) | 8% (49–66% in the storm windows) | 12% | about 105–120 | 89 MV |
| 35.6 km (the crown above sea level) | 10% (60–75%) | 13% | about 115–125 | 151 MV |

The 3-hourly outputs undercount episodes; in the ten-minute windows a place passes the threshold anew every few hours.
The box is a flat lowland at 6 km cells, so the summit's terrain, the tower's own charge and the leaders' growth stay
outside this screen. On the occurrence map the summit's storms come 0.5–3.7 times as often as the box's.

### 5. Haze in the sky

[illumination/aerosol](../../../illumination/aerosol/haze.py) turns the aerosol study's modes into optics at 550 nm
with Mie theory ([mie.py](../../../illumination/aerosol/mie.py), checked against the Rayleigh limit, large spheres and
Bohren and Huffman's example), kappa-Koehler growth at each region's humidity and OPAC's refractive indices
([haze.json](../../../illumination/aerosol/results/haze.json)).

| Kind of region | Haze optical depth, day (clean / central / loaded) | Visibility, central | Noon direct share, clear → central |
|---|---|---|---|
| Seas | 0.006 / 0.022 / 0.033 | 146 km | 0.49 → 0.49 |
| Wet land | 0.012 / 0.030 / 0.59 | 198 km | 0.49 → 0.48 |
| Fog desert | 0.012 / 0.035 / 2.2 | 189 km | 0.49 → 0.48 |
| Polar dry land | 0.002 / 0.005 / 0.045 | 177 km | 0.49 → 0.49 |

The clear air alone passes 49% of the noon sunlight directly and 28% as skylight, and sees 282 km. Practical dusk comes
81.8 hours after an equatorial sunset in clear air and at 80.7–81.8 hours with any of the hazes. The light calendar's
clear sky describes ordinary days.

### 6. Ledgers

[ledgers.py](ledgers.py) gathers the joint quantities. The air holds 3.1×10¹⁸ kg.

| The air | Rate | Cycle time |
|---|---|---|
| Loss, warm titania case, cycle mean, no magnets | 1.4 kg/s | 70 billion years |
| Loss, with the September magnets | 0.044 kg/s | 2.2 trillion years |
| Loss, cooler titania cases | 0.61–0.80 kg/s | 120–160 billion years |
| Hydrogen to space | 0.0024–0.0027 kg/s (77–86 t a year) | the air's 114 Mt of hydrogen in 1.4 million years |

The sky fleet holds 48,000–71,000 t of hydrogen in its liners and ferries. Losing 1–10% of it a year releases
480–7,100 t; each megatonne held in the air adds 0.005 ppm to its 0.53. With no OH under the shield, hydrogen leaves
through soils that take it up, as Earth's do within about two years, or through escape.

Energy runs across six orders: the Moon absorbs about 8 PW of sunlight; the summit metropolis would use 200 GW, a
forty-thousandth of it but 80 W/m² over the city; the regional magnets' refrigeration takes 9–15 GW; the film plant
0.15–2.9 TW; the tiles absorb 65–83 PW in orbit; and the fleet's night half intercepts about 850 PW. A terawatt used
on the Moon adds 0.026 W/m² to its mean heat.

### 7. The regional magnets on the map

[magnets_map.py](magnets_map.py) turns the shield's screened loop designs, whose longitudes the shield left free,
through 0–85° and lays them on the 28% atlas with the heritage register and the science targets
([magnets_map.json](results/magnets_map.json)). The 0.5 mT line, which the guidance for implanted medical devices
keeps people outside, lies 206 km from a 500 km loop's cable and 62 km from a 1,000 km loop's.

| Design | Cable | Over water | Heritage sites inside the 0.5 mT line | Science targets inside |
|---|---|---|---|---|
| Four loops of 500 km | 12,566 km | 27–35% | 4–5 at every turn (LCROSS at Cabeus, the Lunar Prospector impact with Shoemaker's ashes, Chandrayaan-3, Luna 25, IM-1, Hakuto-R) | 13–14 |
| Four loops of 1,000 km | 25,133 km | 28–33% | 2–5, fewest at 15°, 70° and 75° | 2–6 |

The south pole gathers the landings and the cold traps, so the smaller loops reach them at every turn. Which turn
to choose weighs these sites, the seas the cables cross and the people who would live near the routes; that
judgement is people's, with a written evaluation for each site.

### 8. Places through a lunar day

[places.py](places.py) gathers the joint numbers at seven places the branches study
([places.json](results/places.json)).

| Place | Night's own light at midnight | Fleet glow (0.1%; rolled 80°) | Flashes per km² a year | Visibility | Nearest sea's monthly tide | Nearest magnet cable (500 / 1,000 km design) |
|---|---|---|---|---|---|---|
| Summit port, 5.4° N, 158.6° W | 0.001 lux | 66; 2.2 lux | 0.046 (0.011–0.083) | 198 km | Mare Ingenii, 2.1 m, 880 km off | 1,400 / 760 km |
| Procellarum by Russell, 32.9° N | 0.44 lux (the Earth 12° up) | 56; 1.9 | 0.0023 | 146 km | 3.7 m | 1,240 / 138 km |
| Eastern Smythii headland | 0.07 lux (the Earth on the horizon) | 65; 2.2 | 0.0067 | 146 km | 0.55 m | 1,560 / 577 km |
| Southern Mare Nubium | 4.3 lux (3.1 at the darkest) | 53; 1.8 | 0.0034 | 146 km | 3.7 m | 1,140 / 303 km |
| Mare Ingenii coast, 34° S | 0.2 lux | 45; 1.5 | 0.0056 | 146 km | 2.1 m | 865 / 580 km |
| Highland box, 44.7° S | 1.4 lux | 35; 1.2 | 0.0019 | 189 km | 2.1 m, 226 km off | 480 / 612 km |
| Polar plain, 75° N | 480 lux (twilight all night) | 16; 0.5 | 0.0001 | 177 km | 3.7 m, 211 km off | 230 / 66 km |

Each place's own night takes the Earth's height and phase over it. Under a high Earth the nearside's nights keep
2.3–5 lux at their darkest; toward the limbs the Earth stands low, and its light falls to 0.44 lux by Russell and
0.07 lux at the Smythii headland, where the libration lifts the Earth above the horizon and lowers it again (the
sea-appearance calendar's dated month gives 0.48 and 0.085 lux at their darkest). The far side and the south's coasts
hold the dark nights the fleet would take away. The climate at each place stays the range the paused climate work carries: CM1 and the GCM disagree on
the warmth and humidity of the air over land (42 against 387 strictly comfortable hours a lunar day on the rings'
land).

## What each line assumed about the others

| Assumption | Made by | Its source's answer now | State |
|---|---|---|---|
| Sunlight is the titania stack's, 5% dimmer, everywhere | Climate, illumination, the light calendar, biology | The shield's window gate stays conditional on the dimmer's form (absorbed or reflected, undecided), and the fleet's window ends at 1,830 km with no twilight band, so much of twilight's light passes the annulus film (98.7% visible, no dimmer) | Changed for the twilight (about 6% brighter); open for daylight |
| The night is dark | Sea appearance, lunar-cycle ecology, the night at the summit | The fleet's night half glows at 25–66 lux at O5's scatter | Changed |
| The flight band avoids storms; the crown is a structure in the air | Infrastructure | Storms flash in the band; the crown launches upward leaders about a tenth of the time | Changed |
| The magnetosphere is the September 1.5×10²¹ A·m² | Electricity (stage 3's boundary) | About 3×10¹⁹ A·m² may hold the warm case under 1 kg/s; muon ionization uses zero cutoff and holds | Holds for stages 1–2 |
| The box's land surface | Electricity | Its surface product names the September drainage grid; the box is all land with no standing water, so the corrected lakes do not enter it | Holds |
| The port's programme | Sky fleet, metropolis brief | It follows the square lattice, whose corrected base widened (long-haul 61% of the travel floor); the chosen form's programme is still to reconcile | Open |
| Hydrogen lift gas is harmless to the air | Infrastructure | Releases stay a few thousandths of a ppm per megatonne; sinks are soils and escape | Holds |
| NO cools the upper air as the middle atmosphere gives it | Shield | Lightning's NO at up to 100 ppb changes nothing at 145–210 K | Holds |
| Twilight light for plants | Biology | The canopy model already takes the sky atlas's twilight; Earth's 0.242 albedo lowers Earthlight by a third, which plants barely use | Holds |
| Comfort is a range | Every line | The climate pause stands | Holds |

## Superseded values

| Quantity | Was | Is | Where |
|---|---|---|---|
| Twilight | five hours | practical dusk 82 h after an equatorial sunset (2.98 lux) | research/decisions.md, Light and time |
| Earth's visible albedo | 0.367 | 0.242 with a measured phase curve | shared/constants.json |
| Lightning in the box | 504 flashes (465 in the old thunder) | 530 in the corrected run | climate/crm, atmospheric_electricity |
| Ring fleet tiles | 17–22 million | about 28 million (46 Gt) | solar_shield_array |
| Shield loss, warm case | 8.4 kg/s at cycle 21 (pre-infrared) | 1.4 kg/s cycle mean, 1.5 at cycle 21 | atmosphere/loss_response |
| Design gust at the summit | 27.3 m/s | 32.4 m/s | summit_tower |
| Chosen form's steel | 90 Mt | 98 Mt | summit_tower/form.py |
| Pressure hulls | 280–350 m | 260–330 m | sky_ships |
| The magnetic comparison's design point | current | the 26 September design point, kept in solar_shield_array/historical | magnetic_architecture.md |

## Against the research threshold for the core

[plan.md](../../plan.md) asks the research for several coherent environmental possibilities, a credible account of
light, functional biological requirements for a varied biome portfolio, quantified representative spatial
opportunities and compatible construction and renewal cases.
- *Environments:* the corrected climate with its range, two seas' waves and every sea's tide, storms with their
  lightning, regional aerosol and the places above. The climate's land disagreement stays open.
- *Light:* the dated light calendar, Earthlight measured, haze now bounded, and the fleet's light a new term that
  sets a requirement on the shield.
- *Biology:* the canopy, fruit and forest mechanics and the lunar-cycle ecology, now with lightning's nitrogen
  (small), the air's lack of an oxidant and the night's light as conditions. Complete life cycles and nutrient
  cycles stay open.
- *Spatial opportunities:* the summit port and metropolis, the coasts, the polar plains, now with storms, tides and
  the magnets' corridors.
- *Construction and renewal:* the ring fleet's inventory, renewal and heat, the summit port's concept sizing; the
  surface-to-orbit link and the fleet's night-side attitude law are the new gaps.

## What comes next

Ranked by what each would settle for the whole, with machine time where it can be given:
1. **The fleet's night.** The author has decided that the fleet avoids lighting the night zones. The requirement's
   value comes from the biosphere and human work; then a night-side attitude law flown in the ring dynamics (rolls
   with the keeping law; hours of flights, as the shield's twelve kept orbits took 55 CPU minutes), the films'
   scatter, and collectors for the light the night half turns away. This belongs in the planned
   research/array-industry branch with the trim and store hardware.
2. **The crown as a lightning conductor.** The tower's induced charge in these fields, its upward leaders and their
   charge, a lightning protection concept for 100-plus triggered flashes a year, and the hydrogen berths beside it.
   A terrain-aware electrified storm needs CM1's field solver over terrain, new work in climate/crm that touches the
   paused climate programme.
3. **The air's chemistry under the shield.** NO's deposition and the lunar night's chemistry in the middle
   atmosphere, then the biosphere's and cities' emissions (methane, organics, carbon monoxide) in an air with no OH.
   About an hour of runs once the chemistry is added.
4. **Nitrogen for the biomes.** Biological fixation as a design requirement of the engineered biosphere.
5. **The chosen form's programme**, so the crown's berths and traffic follow the building the author chose.
6. **The magnets' siting**, weighed with the heritage register's written evaluations.
7. **Atmospheric electricity's stage 3** (the global circuit and the upper boundary), which awaits the author's go-ahead.
8. **Sea appearance's open items** and the climate discrepancy, which stays paused.

## Decisions for the author

- Decided on 8 October: the fleet avoids lighting the night zones, its light redirected away from the Moon and used
  where it can be; lighting over cities might be acceptable if fairly fine spotlights are possible (section 1;
  [register](../../decisions.md#solar-shield-and-habitat-array)). The requirement's value, from the biosphere and
  human work, sets the means: the night-side roll, the films' scatter or both.
- The chosen form's programme, and with it the crown's berths (the 750 m class needs nine on the square lattice's
  programme).
- The regional magnets' design and turn (section 7).
- Whether the rigid ships' adopted limit of about 3 km becomes 2.5–3 km (its supports now sit at 2.6–2.9 km).
- Stage 3 of atmospheric electricity, and the planned array-industry branch.

## How the numbers were made

```sh
climate/gcm/.venv/bin/python -m climate.gcm.convection A28_dim5_moon:20-29          # 2 s
python -m atmosphere.electricity.occurrence                                          # 1 s
OPENBLAS_NUM_THREADS=1 python -m research.studies.joint_synthesis.air_chemistry --processes 3   # 20 min a case
python -m illumination.fleet_light.model                                             # 1 min
python -m research.studies.joint_synthesis.night_sky                                 # 1 s
OPENBLAS_NUM_THREADS=1 python -m atmosphere.electricity.tower                        # 5 min, reads 52 GB of saved fields
python -m illumination.aerosol.haze                                                  # 1 s
python -m research.studies.joint_synthesis.ledgers                                   # 1 s
python -m research.studies.joint_synthesis.magnets_map                               # 9 s
OPENBLAS_NUM_THREADS=1 python -m research.studies.joint_synthesis.upper_air_no       # 8 s
python -m research.studies.joint_synthesis.places                                    # 1 s
```

The air-chemistry run's second case (the per-joule yield's upper end) was stopped once the first showed NO with no
sink; its outcome would repeat the first's, and the product lists it under not_run. The first case was solved
before the runner gained its --only option and its not_run record, which change no physics; the product records
the runner as it now stands. The outputs of the electrified box live on the research drive under
climate/crm/runs.

Sources beyond the products: Bohren and Huffman (1983), Absorption and Scattering of Light by Small Particles,
appendix A (the Mie series); Petters and Kreidenweis (2007), Atmos. Chem. Phys. 7, 1961 (kappa-Koehler growth);
Hess, Koepke and Schult (1998), Bull. Amer. Meteor. Soc. 79, 831 (OPAC refractive indices); Bodhaine et al. (1999),
J. Atmos. Oceanic Technol. 16, 1854 (Rayleigh optical depth); Rizk (1990), IEEE Trans. Power Delivery 5, 1983
(upward leader inception). These were cited from the standard literature and not re-read for this study.
