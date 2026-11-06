# The Moon's orientation under the Open Moon's water

The Open Moon's seas and lakes hold 1.5–1.7×10¹⁹ kg of water, 0.02% of the Moon's mass, about the share Earth's ocean is of Earth. With the water the porous crust takes up over a million years, the total reaches 7×10¹⁹ kg, 0.1%. That is far too little to change the Moon's orbit, its tidal lock or the spin rate the lock sets. The Moon's own departure from a sphere is smaller still, though, in the quantity that fixes which way the Moon faces. The difference between its largest and smallest moments of inertia is 2.5×10⁻⁴ of MR², and the seas' water is 1.5×10⁻⁴ of its mass. That water lies unevenly: the Near-side Sea in the north of the Earth-facing side, the deep South Pole–Aitken sea in the south of the far side, roughly opposite each other through the Moon. That arrangement turns the Moon about 8° as the seas fill, and it moves the point beneath Earth by 10°.

The turn is a settling into a new stable orientation, with a lock 38% stiffer than today's. It threatens neither the orbit nor the Moon's stability. Its effect on the world comes through the grid: places shift up to 8° in latitude and Earth stands about 10° from where it stands today in every near-side sky. Every product drawn on today's latitudes and today's direction to Earth therefore describes a Moon that the full seas will not keep. The project's current work remains the world's starting condition, and it holds for as long as the filling leaves the Moon near today's orientation, a span that the order of filling controls. A second set of full simulations belongs on the turned Moon.

**Evidence state.** A calculation on committed products (the atlas's seas, the drainage product's rain-fed lakes, the groundwater uptake) and on the crustal-loading study's geoid changes. That study is committed on `domain/hydrosphere` as c2c8624a (4 November 2026) and reaches main at the next integration ([geography README, "The crust's response to the water"](../../../geography/README.md) once merged). The calculation passes the tests listed under "Checks", and a reviewer re-derived it independently with separate code (below). Its method is now a tested module on `domain/geology` (geology/orientation.py, commit 45c3d817), whose tests reproduce this report's elastic state and places, and that lane computes the long-term turn from the present Moon's thermal state; its 23 tests pass, and no independent re-derivation of the long-term response has been made. This report's own scripts are on the research drive with their hashes.

## The water's mass

| Part of the inventory | Mass |
|---|---|
| Seas at the atlas's 28% (one level on the dry Moon's geoid) | 1.07×10¹⁹ kg, a 283 m layer over the whole Moon |
| Extra water to hold 28% once the water's gravity and the crust's sag are counted | +0.09×10¹⁹ kg within years; +0.14–0.26×10¹⁹ kg once the mantle has flowed (0.22 at Te 40 km) |
| Rain-fed lakes | 0.31×10¹⁹ kg, and about 2% more as their floors sag |
| **Standing water** | **1.47×10¹⁹ kg at first; 1.52–1.66×10¹⁹ kg relaxed** |
| Water held in the crust's pores, central uptake | 0.39×10¹⁹ kg at 1,000 years, 1.95 at 10,000, 4.1 at 100,000, 5.1 at a million |
| Air at 1.2 atm | 0.28×10¹⁹ kg (a 75 m water equivalent) |
| **Everything after a million years** | **about 7×10¹⁹ kg, 0.095% of the Moon's 7.346×10²² kg** |

Earth's ocean is 0.022% of Earth's mass. The seas' extra water comes from the crustal-loading study: the water's own gravity raises the sea surface over deep water, and the crust sags under the load. The pores' water comes from the [groundwater uptake](../../../geography/README.md) at its central values.

## Why the water turns the Moon

A synchronously rotating Moon settles with its axis of least moment of inertia (A) pointing toward Earth and its axis of greatest moment (C) along the orbit's normal. From lunar laser ranging's β = (C−A)/B = 6.310×10⁻⁴ and γ = (B−A)/C = 2.277×10⁻⁴ (as tabulated by Murphy 2013) and C/MR² = 0.393112 at R = 1738 km (Williams et al. 2014), C − A is 5.50×10³¹ kg m² and B − A is 1.99×10³¹ kg m². These moments reproduce GRAIL's J2 and C22 (GL0420A) to 0.01%; read the other way, GRAIL's J2 and C22 with β give C/MR² = 0.393131. Water placed off those axes adds products of inertia that tilt them. The atlas's seas give I_xz = −7.35×10³⁰ kg m² (x toward 0° N, 0° E; z north), 13% of C − A, because the South Pole–Aitken sea's 4.8×10¹⁸ kg and the Near-side Sea's 3.4×10¹⁸ kg both lie in the quadrants of positive xz. The rain-fed lakes (3.1×10¹⁸ kg) add I_xy = −1.5×10³⁰ kg m², which moves the face eastward, and oppose the seas' I_xz slightly. The water settling under its own gravity adds about 9% to the seas' and lakes' products of inertia. The Moon's principal axes turn by the angle that diagonalises the new inertia tensor, and the Moon follows them. It has turned this way before: the South Pole–Aitken basin's own mass is thought to have reoriented it by about 15° (Keane & Matsuyama 2014), and the Procellarum region's heat to have moved its pole by about 6° (Siegler et al. 2016).

The degree-2 part of the crustal-loading study's geoid change N carries the whole redistribution: the water as it settles under its own gravity, the lakes, and the crust's response through the load Love numbers. Its inertia change is ΔI = −(2MR²/3)Q, Q being the degree-2 quadratic form of N/R. With seas held at 28%:

| Case | Pole moves | toward | Point under Earth moves to | B − A | C − A | Axis to the ecliptic | Free libration in longitude |
|---|---|---|---|---|---|---|---|
| Today | — | — | 0°, 0° | — | — | 1.54° | 1,045 days by the rigid formula; 1,056 observed |
| Rigid Moon, the water's own gravity | 8.14° | 167° W | 8.1° N, 6.2° E | +38.7% | +3.7% | 1.62° | 888 days |
| Elastic, within years | 8.05° | 167° W | 8.0° N, 6.1° E | +38.1% | +3.6% | 1.61° | 890 days |
| Relaxed, Te 150 km | 4.93° | 166° W | 4.9° N, 3.9° E | +20.7% | +0.8% | 1.56° | 952 days |
| Relaxed, Te 75 km | 3.83° | 166° W | 3.8° N, 3.0° E | +15.3% | +0.3% | 1.55° | 973 days |
| Relaxed, Te 40 km | 2.74° | 165° W | 2.7° N, 2.2° E | +10.4% | 0.0% | 1.54° | 995 days |
| Relaxed, Te 20 km | 1.55° | 165° W | 1.5° N, 1.3° E | +5.6% | −0.1% | 1.54° | 1,017 days |

Holding only today's atlas volume, without the extra water that keeps 28% cover, the elastic turn is 7.6°. The South Pole–Aitken sea does most of the turning. On its own, as computed directly from the atlas, it moves the pole 6.5°; the Near-side Sea alone moves it 2.2°. The overall share changes the result little, because the South Pole–Aitken basin fills first: the atlas's seas alone, without the lakes, turn the pole 6.5° at 20% cover, 8.4° at 28% and 10.1° at 35%. On the same basis, holding the South Pole–Aitken sea 1 km below the common level gives 5.6° (seas 26.8%), 2 km gives 3.5° (25.1%), and leaving it empty still leaves about 2°.

**Elastic and relaxed.** The crust's elastic response alters the turn by about 1%, so within years of the water's arrival the turn is the water's own: about 8°. The present Moon then yields slowly. Published present-day temperature profiles with olivine's laboratory flow laws give a mantle viscosity falling from about 10²⁸ Pa s at 100 km depth to 10¹⁷ at 1,000 km, so a degree-2 load meets a rigid layer that acts as an elastic lid about 730 km thick over a thousand years and 300 km over ten million (geology/README.md on `domain/geology`, commit 45c3d817, 4 November 2026). In the central case the pole eases to 7.62° at 1,000 years, 7.24° at 100,000 and 6.59° at 10 million. The stiff bound keeps 7.62° at 10 million years, and the far side's colder mantle under the South Pole–Aitken sea keeps 6.95°. The thin-shell states above (1.5–4.9°) are limits the present Moon approaches only if its shallow mantle is as warm as the warmest published profile, which brings the turn to 1.9° within about a thousand years; the measured heat flow leans toward the cooler profiles. The Moon's retained non-hydrostatic figure does not decide the question: a shell that carries a frozen-in figure can still let new loads sink in part, because the measured figure is already the net of whatever compensation has occurred. Keane & Matsuyama (2014) model the fossil figure itself as partly relaxed under 25–50 km lithospheres. A GRAIL-era elastic thickness under the near-side maria remains open.

**The hydrostatic part.** The hydrostatic response to today's rotation and tide is 4.5% of J2 and 12% of C22: J2 = (5/6)k_f q and C22 = (1/4)k_f q, with q = 7.59×10⁻⁶ and the fluid Love number k_f = 1.436 from the Darwin–Radau relation (Keane & Matsuyama 2014 use 1.44). That part re-forms about the new axes only where the interior flows, so it leaves the within-years turn unchanged and enlarges the long-term turns by about 6%, to 1.6–5.2° (2.9° at Te 40 km). Re-formed under an otherwise elastic response, it would bring the pole's move to 8.5° and the face's to 10.8°.

## The orientation through time

The water's gravity acts as soon as the water arrives, so the turn tracks the filling. As the inertia tensor changes, the Moon follows its new axes through its free librations, whose periods are 1,056 days in longitude, about 81 years for the free precession of the pole and 74.6 years for the free wobble (Rambaux & Williams 2011). A filling that takes many of the longer periods lets the Moon follow smoothly; a filling within decades sets it wobbling by degrees, which then damps over a time not computed here.

**The order of filling sets the path.** With the lakes in place and the South Pole–Aitken sea left for last, the Near-side Sea alone turns the pole 2.5° and every sea but the South Pole–Aitken sea 2.2°. The point beneath Earth then moves 5–5.5°, mostly in longitude because of the lakes. Without the lakes the turn is about 2° in each. The Moon stays near today's orientation until the South Pole–Aitken sea fills, which makes the order of delivery a design lever over how long today's grid holds.

**After the seas are in.** The crust's pores fill at low latitudes first, because the cold polar crust thaws last. That pore water adds little turn, but it raises C − A by up to 43%, which slightly stiffens the pole against the seas' pull and tilts the spin axis further from the ecliptic:

| Seas, lakes and pore water | Elastic: pole moves | Elastic: axis to the ecliptic | Relaxed Te 40 km: pole moves | Relaxed Te 40 km: axis to the ecliptic | Relaxed Te 150 km: pole moves | Relaxed Te 150 km: axis to the ecliptic |
|---|---|---|---|---|---|---|
| Seas and lakes only | 8.05° | 1.61° | 2.74° | 1.54° | 4.93° | 1.56° |
| Pores at 1,000 years | 7.46° | 1.80° | 2.69° | 1.59° | 4.73° | 1.66° |
| Pores at 10,000 years | 7.58° | 2.24° | 2.96° | 1.70° | 5.04° | 1.89° |
| Pores at 100,000 years | 7.56° | 2.51° | 3.09° | 1.77° | 5.16° | 2.03° |
| Pores at 1 million years | 8.05° | 2.30° | 3.17° | 1.72° | 5.38° | 1.92° |

In every case the pole moves toward 165–173° W and the point beneath Earth moves north and slightly east.

**The orientation for the second set of simulations.** The elastic state with the seas and lakes in place is the orientation set for the second set of full simulations (register, 4 November). In today's coordinates it places the north pole at 81.95° N, 167.14° W, the south pole at 81.95° S, 12.86° E, and the point beneath Earth at 8.00° N, 6.09° E; its rotation matrix and a tested latitude–longitude conversion are in `geology/results/moon_orientation.json` and `geology/orientation.py` on `domain/geology` (commit 45c3d817). It is the state people live with: in the central case the pole eases from it by 0.3° in a century, 0.4° in a thousand years and about 1° in a million, and in the stiff case by under 0.1° in a thousand years. Only a shallow mantle as warm as the warmest published profile would settle it near 1.9° within a millennium. With the current studies at today's orientation, the two sets bracket the range the Moon can take. The lakes move the face about 3° in longitude, and a climate on the turned Moon will route its lakes differently, so the orientation is recomputed once with that climate's lakes before the simulations' products are treated as settled.

## Where places move

Positions in the turned Moon's own latitude and longitude, from today's IAU centres, the conservation study's site register and the Jura box's centre:

| Place | Today | Elastic, seas and lakes | Relaxed Te 40 km, mature pores | Relaxed Te 150 km, mature pores |
|---|---|---|---|---|
| Jura Island's centre (the CM1 box) | 45.0° N, 36.1° W | 39.4° N, 37.0° W | 42.9° N, 36.0° W | 41.3° N, 36.3° W |
| Montes Jura | 47.5° N, 36.1° W | 41.9° N, 36.5° W | 45.4° N, 35.9° W | 43.8° N, 36.0° W |
| Mare Imbrium | 34.7° N, 14.9° W | 27.5° N, 18.9° W | 31.9° N, 16.2° W | 29.9° N, 17.3° W |
| Mare Serenitatis | 27.3° N, 18.4° E | 19.3° N, 11.9° E | 24.1° N, 16.0° E | 21.9° N, 14.2° E |
| Mare Tranquillitatis | 8.3° N, 30.8° E | 0.7° N, 24.5° E | 5.3° N, 28.5° E | 3.2° N, 26.8° E |
| Apollo 11 | 0.7° N, 23.5° E | 7.2° S, 17.4° E | 2.5° S, 21.2° E | 4.6° S, 19.6° E |
| Copernicus | 9.6° N, 20.1° W | 2.8° N, 25.8° W | 7.0° N, 22.1° W | 5.1° N, 23.6° W |
| Tycho | 43.3° S, 11.2° W | 50.5° S, 21.1° W | 46.2° S, 14.8° W | 48.2° S, 17.5° W |
| Korolev | 4.2° S, 157.4° W | 3.7° N, 163.6° W | 1.1° S, 159.7° W | 1.1° N, 161.3° W |
| Hertzsprung | 1.4° N, 128.7° W | 7.7° N, 134.4° W | 3.9° N, 130.8° W | 5.6° N, 132.4° W |
| Chang'e 4 | 45.4° S, 177.6° E | 37.6° S, 173.2° E | 42.4° S, 176.2° E | 40.2° S, 174.9° E |
| South Pole–Aitken sea's centre | 53° S, 169° W | 44.9° S, 174.9° W | 49.8° S, 171.1° W | 47.6° S, 172.7° W |
| Shackleton | 89.7° S | 81.8° S | 86.7° S | 84.5° S |
| Peary | 88.6° N | 80.6° N | 85.5° N | 83.3° N |

Near-side places move toward the south, about 8° on the central meridian and less toward the limbs. Far-side places move north. Today's polar craters end near 82° in the elastic state, where the noon Sun stands 6.5–10° high through the year. That shrinks their permanent shadow, although deep, fresh craters can keep shadowed floors at such latitudes, and new shadowed ground forms about the new poles in other terrain. Shadow makes a cold trap only in vacuum: under the design air a shadowed floor warms toward the air's temperature, and the turned design years hold the air poleward of 70° at 301.6 K, so the turned poles' shadow matters for the build's airless years and for sampling the ice. The polar cold-trap work therefore needs the turned poles and the higher tilt as inputs, and sampling the polar ice where it lies has one more reason to finish before the seas fill.

**Jura's new climate, read from the design climate.** The design run (`A28R_301k_design_l15_snow_s934`, years 10–19) has almost no temperature gradient between 39° and 45° N. Moving there therefore gives Jura more sunlight but no warmer air, and noticeably more rain on land:

| Around the Moon, design climate | Land at 45° N | Land at 39.4° N | Sea at 45° N | Sea at 39.4° N |
|---|---|---|---|---|
| Air at 2 m | 26.7 °C | 26.5 °C | 29.5 °C | 29.3 °C |
| Sunlight absorbed at the ground | 76 W/m² | 82 W/m² | 89 W/m² | 99 W/m² |
| Rain | 261 mm/yr | 400 mm/yr | 224 mm/yr | 223 mm/yr |
| Rain less evaporation | 0 | +34 mm/yr | −791 mm/yr | −867 mm/yr |
| Cloud cover | 32% | 35% | 32% | 31% |
| Range over the lunar day | 5.8 K | 5.6 K | 2.3 K | 2.4 K |
| Range of the lunar-day means over the year | 0.4 K | 0.4 K | 0.2 K | 0.3 K |

This reads the existing climate at Jura's new latitude around the whole Moon. On a turned Moon the seas and land also move relative to the air's circulation, so the rain in particular needs the second set of simulations. The model's T21 cells are about 170 km across, and Jura falls on the edge between a land cell and a sea cell.

## Seasons and the spin axis

For a synchronous rotator in Cassini state 2, the spin axis's angle ε from the orbit's normal satisfies tan ε = sin i / (cos i − x), with x = (3/2)(n/|Ω̇|)(M_E/(M_E + M_M))((C − A)/C) cos ε, i = 5.145° and an 18.6-year nodal period. It is the orbit average of the potential (3/4)n²(C − A)sin²ε that Wisdom (2006) uses, and from LLR's β it gives today's tilt to the ecliptic as 1.541° (1.550° with the orbit's eccentricity factor), against the measured 1.542°, with nothing fitted. The tilt follows C − A. The seas alone would lower it to 1.47°; the rain-fed lakes raise it again, to 1.61° with both. While the low latitudes' pores fill, the tilt rises to 1.8–2.5° if the crust stays elastic, and to 1.6–2.0° if it relaxes.

The tilt sets the seasons. Over the eclipse year (346.6 days), the lunar-day mean sunlight on flat ground at latitude φ varies by ±(π/2) tan φ sin ε relative to its mean, checked against direct integration within 2%. At 45° that is ±4.2% today, the same as moving 2.4° in latitude and back. It becomes ±4.4% once the seas are in, and at most ±6.9% at the elastic 2.5°, like moving 3.9° and back. Earth's changing distance from the Sun varies sunlight everywhere by ±3.4% over the year, and Earth's own seasons at 45° by about ±50%. The design climate's lunar-day means vary by 0.2–0.4 K over its year, so the Moon stays nearly seasonless at low and middle latitudes. At the poles the tilt sets the midsummer Sun's height: from 1.54° to 2.2–2.5°, flat polar ground receives 40–60% more midsummer sunlight. The pole's move outweighs the tilt there, as the places above show.

## Stability

The turn happens inside the Moon: the body moves beneath a spin axis that stays fixed in space, and the orbit does not change. The orbit can change only through the delivery's net momentum. All the water and air together, about 7×10¹⁹ kg, arriving at the Moon's escape speed of 2.38 km/s and all pushing one way along the orbit, would change the Moon's 1.02 km/s orbital speed by 2.3 m/s and the month by about four hours. A controlled delivery can hold that near zero; where the deliveries put their energy remains open in [resources/questions.md](../../../resources/questions.md), item 5. Earth's tides feed the Moon's orbit about 0.12 TW, moving it outward 3.8 cm a year. Tidal friction in the lunar seas, which nothing has yet computed, would take about ten billion years at a gigawatt to change the orbit's energy by 1%. The Moon's extra 0.1% of mass raises Earth's tides by 0.1%.

The lock stiffens as B − A grows by 38%, and the free libration in longitude shortens from 2.9 to 2.4 years. For an orbit inclined 5.145°, Cassini states 1 and 4 exist only when the ratio of the spin axis's precession constant to the node's regression rate exceeds (sin^(2/3) i + cos^(2/3) i)^(3/2) = 1.31. Today the ratio is 0.23, and the water raises it to at most 0.33, when the pores have filled for 100,000 years under an elastic crust. The Moon therefore stays in Cassini state 2, whose tilt changes smoothly with the moments.

The water opens no route back to volcanism. The Moon's large eruptions ended about a billion years ago from sources hundreds of kilometres down. A 10 MPa load raises melting points by about 1 K, and on Earth it is unloading, as ice caps melt, that raises volcanism. The water reaches at most about 20 km, where the pores close and the rock is at most about 600 K, well below the 900 K or more at which even wet rock begins to melt. Water circulating through the hot, shattered upper crust will make warm springs, which the groundwater work already notes. Moonquakes are the open hazard. Apollo's seismometers recorded 28 shallow moonquakes, with body-wave magnitudes estimated at 3.6–5.8 (Oberst 1987, as reviewed by Nunn et al. 2020), and the Moon carries young thrust faults (Watters et al. 2019). Kilometres of water, and pore pressure rising in the crust as it soaks in, could trigger slip on faults already near failure, as filling reservoirs does on Earth, mostly during the centuries to millennia of filling and soaking. The seas press 1.6 MPa per kilometre of depth, about 10 MPa on the deepest South Pole–Aitken floor. Nothing here has examined this question.

## What sits on today's grid

The turn leaves the atlas's geometry, the seas' volumes, the crustal loading, the shield and the air as they are. The shield's frame is the orbit and the Sun, and the spin axis keeps its place in space. The permanent tidal and rotational bulge, already part of the atlas's equipotential, turns with the face; at Earth's 12.9 m tide on the Moon, a 10° turn moves sea level against the ground by at most a few metres. Everything computed from sunlight by latitude, or from Earth's direction, sits on today's grid:
- the climate model's design cases and the products read from them, among them the climate targets, heat levers, the resort at 301 K, and the CM1 rings and boxes at fixed latitudes;
- geography's drainage and rain-fed lakes, the groundwater uptake's thaw times by latitude band and the water table;
- the coast water, the lakes' and seas' oxygen and methane studies, light limitation, food systems by band, and the biome map on `domain/ecology`;
- the tides' forcing pattern, and the views of Earth in sea appearance, living light and cloud twilight;
- the conservation study's statements about the face seen from Earth and the polar sites, and the polar cold-trap work still to come.

Those products remain the world's starting condition. A second set of full simulations on the turned Moon, followed by a plan for carrying the developed world across to it (biomes, lakes, coasts, places and their names), is the work this finding sets.

## Checks

Recorded in `verify_output.txt` beside the scripts:
- The grid's cell solid angles sum to the sphere within 10⁻¹². A 1° cap's inertia matches the point-mass value within 2.5×10⁻⁴. Degree 1 is recovered exactly, and a load's inertia comes out the same whether taken from its mass or from its geoid.
- A, B and C from LLR's β, γ and C/MR² reproduce GRAIL's J2 and C22 within 0.01%.
- The Cassini relation reproduces today's 1.542° within 0.0013° with nothing fitted, and the Darwin–Radau relation returns k_f = 3/2 for a homogeneous body.
- The seasonal formula matches direct integration of the lunar day's sunlight within 2% at 30°, 45° and 60°.
- The atlas's seas computed directly (I_xz = −7.35×10³⁰ kg m², pole 8.4°) agree with the study's rigid geoid change (−7.44×10³⁰ kg m², 8.1° with the lakes). The two grids share their row order.
- An independent reviewer, working from the brief and its own code (`review/r1_lakes.py`, `review/r2_rotation.py`, findings in `review/review.md`), reproduced the moments, the seas' and lakes' masses and products of inertia (the lakes' 3.1086×10¹⁸ kg from LOLA's 16 px/deg heights), every turn in the tables above, the Cassini relation and threshold, and Jura's new latitude. It found one error, since corrected: the Darwin–Radau relation had been written with 4/(1 + k_f) where 5/(1 + k_f) belongs, giving k_f = 0.949 for 1.436, which understated the hydrostatic part. It also showed that the 1.61° tilt with the seas depends on the lakes, that the threshold ratio peaks at 0.33, and that the Moon's retained figure leaves the long-term turn open. The fixed Darwin–Radau relation now has its own test (a homogeneous body gives k_f = 3/2).

Not done:
- a dynamical simulation of the turn and its wobble;
- the iteration of the orientation with lakes from a turned climate;
- full reading of the cited papers beyond the values used.

## What it leaves out

- The thin-shell relaxed cases use the crustal-loading study's lithospheres over a fluid mantle. The geology lane's response from today's thermal state replaces them for the long term; its range rests on the shallow mantle's temperature, grain size and water, through flow laws measured on Earth's olivine.
- The near side's warmer mantle (Park et al. 2025) and other lateral differences are absent, as are the fluid core's own inertia, under 0.1% of the Moon's, and the 0.02° between the topography's and the gravity field's frames.
- The Sun's torque on the Moon's figure, about 0.6% of Earth's, is left out of the Cassini relation.
- The pore water's load uses the uptake's latitude bands on the atlas's sea, lake and land columns. Its Moon-wide total exceeds the uptake's own by 11% at 1,000 years, 2% at 10,000 and under 1% later.
- The lakes, the pores and the seas' extra water keep their present-climate patterns, and the pore water takes the seas' load response.
- The free modes' damping time, which sets how long a fast filling's wobble lasts, is not computed.
- No moonquake estimate is made.

## Reproducing it

The scripts are on the research drive in `/media/projectspace/terluna-research/session-tools/rotation_check/`. `verify.py` holds the tests and every table above except Jura's climate, and its output is `verify_output.txt`; it supersedes the exploratory `rotation.py`, `check.py`, `levers.py`, `spa.py` and `jura.py`, which reached the same numbers with a fitted small-angle form of the Cassini relation. `jura_climate.py` reads the design run and needs `climate/gcm/.venv/bin/python` for netCDF4. The reviewer's scripts and findings are in `review/`. Run everything from main's checkout under `nice -n 10`, one thread and `ulimit -v 6000000`; each script takes under a minute.

| Script | SHA-256 (first 16) |
|---|---|
| verify.py | 9cf0eb38b33e6857 |
| jura_climate.py | 6c875a3fb7b920fd |
| rotation.py, check.py, levers.py, spa.py, jura.py | a431a5bafec9dbfc, 00705f3883e0209e, 6c19f2b92661390c, 0eeccb8f695981a3, bfa66b647115fb6a |
| review/r1_lakes.py, review/r2_rotation.py | 0ee635242d67c0c0, 0db93ee2a835943d |

| Input | SHA-256 (first 16) |
|---|---|
| geography/products/atlas_28pct_4ppd.npz | dc283827f349fdc1 |
| geography/products/drainage_28pct_16ppd.npz | 98842fc2a648f889 |
| crustal_loading_nmax359.npz (research drive, crustal-loading/) | 4b408e82d708c7af |
| geography/results/crustal_loading.json (domain/hydrosphere, c2c8624a) | e1b1b2fca964fb98 |
| geography/results/groundwater_uptake.json | 2e0cbd10dc37b906 |
| Design run `A28R_301k_design_l15_snow_s934`, years 10–19 | research drive, gcm-runs/ |

## References

- Keane, J. T. & Matsuyama, I. (2014). Evidence for lunar true polar wander and a past low-eccentricity, synchronous lunar orbit. *Geophys. Res. Lett.* 41, 6610–6619.
- Konopliv, A. S. et al. (2013). The JPL lunar gravity field to spherical harmonic degree 660 from the GRAIL Primary Mission. *J. Geophys. Res. Planets* 118, 1415–1434 (GL0420A's family of fields).
- Murphy, T. W. (2013). Lunar laser ranging: the millimeter challenge. *Rep. Prog. Phys.* 76, 076901 (arXiv:1309.6294).
- Nunn, C. et al. (2020). Lunar seismology: a data and instrumentation review. *Space Sci. Rev.* 216, 89.
- Oberst, J. (1987). Unusually high stress drops associated with shallow moonquakes. *J. Geophys. Res.* 92, 1397–1405.
- Park, R. S. et al. (2025). Thermal asymmetry in the Moon's mantle inferred from monthly tidal response. *Nature*, doi:10.1038/s41586-025-08949-5.
- Rambaux, N. & Williams, J. G. (2011). The Moon's physical librations and determination of their free modes. *Celest. Mech. Dyn. Astron.* 109, 85–100.
- Siegler, M. A. et al. (2016). Lunar true polar wander inferred from polar hydrogen. *Nature* 531, 480–484.
- Watters, T. R. et al. (2019). Shallow seismic activity and young thrust faults on the Moon. *Nature Geosci.* 12, 411–417.
- Williams, J. G. et al. (2014). Lunar interior properties from the GRAIL mission. *J. Geophys. Res. Planets* 119, 1546–1578.
- Wisdom, J. (2006). Dynamics of the lunar spin axis. *Astron. J.* 131, 1864–1871.
- Briaud, A. et al. (2023), Maxwell et al. (2023), Solomon & Head (1980) and Turcotte et al. (1981), through the crustal-loading study's sources file.
