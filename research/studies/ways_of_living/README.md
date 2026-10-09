# Ways of living on and around the Open Moon, and what each holds and costs

**9 October 2026.** Steps three and four of the population and lifestyle work the author approved: the ways of living,
built from the project's designs, and what each holds and costs. Each setting is placed on the 28%
[atlas](../../../geography/README.md#atlas-at-the-selected-water-share) with the drainage estimate's lakes and rivers,
and its physical situation is read from the products: the air, its warmth through the lunar day, rain and water,
the light calendar, hazards and what the third dimension offers there. Its costs per billion residents follow, and
the [shared population cases](../../../shared/scenarios/population.json) are set against them through mixes of
settings. Choosing a distribution is the author's, in step six; this study gives the costs that choice weighs.

[run.py](run.py) writes [results/ways_of_living.json](results/ways_of_living.json) (schema
`terluna.research.ways-of-living/1`) in about ten seconds. Numbers marked *screen* are computed here; *derived* marks
arithmetic on cited values. Other numbers come from the product named beside them.

## Findings

1. **Every shared case fits the Moon's land; the towns themselves take 0.6–10% of the dry land.** Across the five
   cases, the aerial shares of 0–75% and four mixes of surface settings, settlement covers 0.6–10.1% of the 23.4
   million km² of dry land at the ecology screen's densities (*screen*). The homes differ greatly in size. Each
   billion people takes 3.6–5.1% of the land along the seas and great lakes (the range runs from 30% kept
   wild to half), 14–20% of the highlands above 4 km and 25–35% of the narrow limb band (*screen*).

2. **Food is what grows with the Moon's people, and its answered cost is water lifted from the fresh seas.** Ten
   billion lunar residents on the lunar design's plant-rich diet need 1.4–3.9 million km² of crops at the crop
   model's yields, 35–224% of the rain-fed land outside a kept 30% (ecology's per-billion figures; *derived*). On the
   rain alone, wettest land first, the low end takes a third of the rain belt's land outside the kept share; the high
   end fills the wet land to 45° and leaves 2.2 million km² of crop-model yield to irrigate: 8,700 km³ a year, 1.2
   times the rivers' flow, lifted from the seas for 1.5 TW (*screen*). Irrigated to the crop model's yield from the
   start, crops and towns take 20–59% of the rain belt's land outside the kept 30% for 0.18–2.0 TW of pumping
   (*screen*). The seas take 1,100–89,000 years to reach a gram of salt per kilogram (ecology branch, from the
   resources screen), so the water stays fresh through the near horizons.

3. **People's heat costs 0.08–1.8% more dimming across the cases.** At 2.0–6.0 kW of primary energy a person with
   industry in orbit, 6–15 billion lunar residents release 9–68 TW on the Moon (*screen*). With the sky towns'
   operating power, the irrigation lift and the hydrogen remade after losses of 1–10% a year, the Moon takes 9–128 TW,
   answered by 0.08–1.8% more dimming at 0.009–0.014% a terawatt (*screen*). In total25 at its own split, 17–72 TW
   and 0.15–1.0% (*screen*).

4. **The dark night fills before the land does.** The lunar design's shielded lamps put 791–1,420 m² of land a
   person under a light-polluted sky (ecology), so 6–15 billion lunar residents, those aloft counted with those on
   the ground, light 20–91% of the dry land (*derived*). At the high end, 10 billion or more light more land than lies
   outside a kept half. Each billion living in the far side's twilit-night places lights 45–81% of their home
   (*derived*).

5. **Five billion people aloft live in about 870,000 towns 4–6 km apart; their mass, gas and power are answered costs
   and storms set their spacing.** In the inhabited-volume reference (500 m envelopes at 10 km, 50 t a resident), five
   billion aloft need 425 Gt of districts and carriers, a 75 Gt unused reserve and 48 Gt of lifting gas holding 37 Gt
   of hydrogen, and shade 1.8% of the Moon (*screen*). Spread over the sky clear of terrain outside a kept 30–50%,
   towns stand 4.7–5.5 km apart and 15–21 lie under each median storm, 24 km of rain; kept within the rain belt,
   3.3–3.9 km and 30–42 (*screen*). At 11 billion aloft (total30 at 75%) the rain belt's spacing falls to 2.2–2.6 km,
   inside the inhabited-volume three-diameter grid of 3 km, so towns would stack in layers (*screen*).

6. **The lifting gas's upkeep is the largest cost of living aloft.** Each percent a year of hydrogen lost takes
   0.43–0.44 TW per billion people aloft to remake by electrolysis and keeps 0.69 ppm of hydrogen in the air while
   soils take it up over about two years (*screen*; the joint ledger's sink). At the joint ledger's 1–10% a year, five
   billion aloft need 2.1–22 TW and keep 3.5–35 ppm of hydrogen in an air that holds 0.53 ppm now (*screen*). The gas
   barrier's tightness sets this cost.

7. **A town drifting at 10 km meets the Sun every 15–23 days and keeps one time of day for at most 10–20 days.** In
   the [zonal-winds diagnostic](../../../climate/gcm/zonal_winds.py)'s 3-day means, the eastward wind at 10 km carries
   a drifting town into the Sun's path: its solar day runs 15–23 days at the median by release latitude and 9–38
   days at the 10th and 90th percentiles, against the ground's 29.5 (the diagnostic's routes). Drifting towns gather
   over the equator, 1.4–1.9 times their share by area between 10° S and 10° N, and toward the north pole, 1.5–2.8
   times poleward of 60° N, while their share falls to half at 30–50° (the diagnostic's occupancy).

8. **Streamlined hulls keep a place for 0.09 kW a resident and a time of day for 0.04–0.5 kW.** The sky-ship study's
   3 km hull has a twenty-third of the 500 m sphere's drag for the same volume (*derived* from its cruise power).
   Holding a place against the 3-day-mean wind at 10 km costs 0.089 kW a resident, against 6.9 kW for a sphere held
   in 5 m/s; following the Sun takes 0.18–0.50 kW a resident at the best hour between 47° S and 47° N, where the wind
   blows against the motion, and 0.04–0.19 kW poleward of 60° (*screen*). Gusts and storms add to both.

9. **Water is answered everywhere for little energy.** The middle latitudes' homes generate 2.9–5.2 km³ of runoff a
   year each and the high latitudes 0.39 km³, with standing water 30–41 km away at the median (*screen*). Lifting a
   billion residents' 18–55 km³ of domestic water a year from the fresh seas takes 3.2–9.5 GW, and the makeup after
   90–99% recovery 0.03–0.95 GW (*screen*). Sky towns gather their own: an envelope's top catches 137 litres a
   resident for each millimetre of rain, 17–684 litres a day by latitude band at the ground's rain rate against
   0.5–15 litres of makeup (*screen*). A town at 10 km floats above the day's cloud base near 7 km, so its rain is
   lower than the ground's; rain by height is not yet computed.

10. **The settings differ most in their light.** Practical dusk ends 82 hours after sunset at the equator, 98 at 30°
    and 115–123 hours in the middle-latitude homes, and poleward of 48.6° it lasts the whole night (*screen*,
    reproducing the register). The far side's equator, where the metropolis stands, has the darkest night on the
    Moon, 0.0014 lux, with 118 hours below a tenth of a lux; under a high Earth the near side's nights keep about
    3 lux at their darkest ([main's places](#inputs-to-bind-at-integration)); at 63°, the high-latitude median,
    twilight keeps 45 lux all night (*screen*).

11. **Mixes of every kind carry each case at its own split; the surface-only end of total30 fills the rain belt.** At
    the cases' own 50% aloft, every surface mix uses at most 5–44% of any home's land outside a kept half, and every
    aerial mix fits the sky (*screen*). With all 15 billion of total30 on the surface, the rain-belt mix puts farming
    country on 88% of its home outside the kept half, and the high end of the diet range irrigated fills 99.5% of
    the rain belt's land outside a kept 30% (*screen*).

12. **Earth, the array and elsewhere carry their own costs.** Up to 10.2 billion people eat within four planetary
    boundaries on Earth (Gerten et al. 2020), so in the Earth-heavier case's 14 billion, 27% of Earth's food grows off
    the land (*derived*). The array's 2–3 billion need 1,980–2,970 Gt of shield at the Stanford torus's 990 t a
    person, 43–64 times the ring fleet's mass and a quarter to a third of a year of the build's matter stream, or a
    quarter of that in a denser design (*derived*, the provisioning population screen). Elsewhere's 2 billion take
    4–12 TW, with no capacity credited.

The cases' lunar populations of 6–15 billion arrive between 2849 and 3261 on the provisioning settlement path,
1–2% a year from 10 million in 2526 (*derived*).

## The ways of living

Each setting's home is a share of the atlas's dry land, assigned in this order: the metropolis's circle; land within
one atlas cell (7.6 km) of the seas and great lakes (1,000 km² and more); land beside the small lakes and the large
rivers (100 m³/s and more); ground 4 km and more above the sea; the rest of the rain belt within 30°; the middle
latitudes to 48.6° by side (nearside, limb, far side); and the high latitudes. Farming country and canopy towns share
the rain belt's interior, and highland towns and the districts held up by terrain share the highlands. Warmth is the
design climate's land air at the home's latitudes moved to its mean height; the GCM and CM1 disagree on the warmth
and humidity of the air over land, and both put the most comfortable land poleward of 30° and on high ground
([register](../../decisions.md#research-practice)). The light is the solved clear sky's.

### The summit metropolis and its port

A beautiful, sustainable metropolis of about 100 million round the Moon's largest sky port on the far side's summit.
Districts gather round crater lakes, and the city flows in under the tower between its legs. Sky boats serve as cars,
gliders fly throughout, and a backbone of metro, rail and sky ferries carries it. The port's rings and disks hold
parks to 15 km above the summit and the Moon's only winter, on the disks at 13.5–15 km
([metropolis brief](../../../habitation/summit_metropolis/README.md)).

- **Where.** 2,464 km² of dry land within 32 km of the tower, 7.7–11.1 km above the sea, 40,600 people per km² over
  it (*screen*). The port holds 1.5 million people on a regular basis, with no homes.
- **Air and warmth.** 0.947 atm at the tower's foot, oxygen like Earth's at 1,950 m, 12.1 °C (summit tower study).
  Over the city's mean ground, 9.5 km above the sea, the land air runs 14.8–19.5 °C through the lunar day, with a
  mean of 14.9–16.7 °C between the register's 1 °C per km and the free-air profile (*screen*).
- **Rain and water.** 13.1 mm a day, the rainiest ground on the Moon; lakes 15 km away at the median, Korolev's
  shore 84 km off and 3.7 km below.
- **Light.** Practical dusk ends 82 hours after sunset and 190 hours lie beyond it, 118 of them below a tenth of a
  lux; at its darkest the sky gives 0.0014 lux, and the Earth never rises (*screen*).
- **Hazards.** 0.046 flashes per km² a year, the most of the seven places (main's places); the crown starts its own
  lightning in some 100–120 episodes a year (joint synthesis).
- **The third dimension.** The tower to 35.6 km, parks for everyone to 3 km above the summit; sky boats to 1.5 km,
  ferries to 2.5 km, gliders on the day's thermals to about 10 km.
- **What limits it and what it costs.** Its planned scale on the summit's crown and upper slopes. Its heat over the
  city is 60–181 W/m² across the band, answered in the city. Its lamps light 79,000–142,000 km² of the Moon's darkest
  sky, a circle 160–210 km in radius, and its food takes 13,900–39,300 km² of crops at the model's yields (*derived*).

### Coastal and lake towns

Harbour towns on the seas and the great lakes: the near-side sea's long coasts, the far side's South Pole–Aitken sea,
islands such as Caucasus and Jura, and the shores of highland lakes such as Korolev and Hertzsprung. The related
coastal communities of the Apennine, Taurus and Le Monnier coasts belong here, beside the undersea heritage sites
([conservation](../conservation/README.md#9-candidate-community-settings)). The seas are gentle, and the monthly tide
lays out flats for days.

- **Where.** 3.95 million km² within an atlas cell of 361,000 km of sea coast and 293,000 km of great-lake shore,
  4–54° of latitude (10th–90th percentile); half on the far side (*screen*).
- **Air and warmth.** 116 kPa and 20.3 kPa of oxygen at the mean 2.3 km above the sea; 19.6–24.8 °C (*screen*).
- **Rain and water.** 3.3 mm a day; the fresh seas and lakes beside every town.
- **Light.** At the median latitude, 20°, practical dusk ends 88 hours after sunset and 178 hours lie beyond it. A
  fifth of this land lies under a high Earth, with about 3 lux at the darkest; half is on the far side.
- **Hazards.** A typical monthly tide of 2.8 m on the sea coasts (0.2–3.7 m by sea), which moves the waterline
  140–280 m on beaches of 1:50 to 1:100 (*derived*); significant waves of 1.4 m on average and 4.9 m at the most in
  the measured cycle on the near-side sea, loading structures a sixth as hard as Earth's waves of the same height
  (wave review, infrastructure review).
- **The third dimension.** Sky ferries and harbours for regional ships; from a 100 m headland the sea is in view for
  18.6 km.
- **What limits it and what it costs.** The shore itself. A billion people at 10,000 per km² cover 3.6–5.1% of it
  outside the kept share, and the whole shore built over would hold 20–28 billion (*screen*). Busy-hour flyers come to
  33 a km², and NASA's lane spacing holds a 23% share of trips by sky boat (*screen*).

### Wetland-edge towns

Towns at the wet margins: the shores of the many small crater lakes and the floodplains of the large rivers, among
reed belts and wetland mosaics, with water travel and boardwalk streets.

- **Where.** 3.58 million km² beside 498,000 km of small-lake shore and the large rivers, 3–45° (*screen*).
- **Air and warmth.** 115 kPa, 20.1 kPa of oxygen; 19.5–24.6 °C.
- **Rain and water.** 4.1 mm a day; 2,214 km³ of runoff a year.
- **Light.** At 16°, practical dusk ends 86 hours after sunset and 183 hours lie beyond it; half the land is on the far
  side.
- **Hazards.** Rain of 17–87 mm an hour beneath flashing storm cores (infrastructure review); wetland gases, which
  last centuries in an air without OH (ecology).
- **The third dimension.** Boardwalks and raised districts over the wet ground; personal wings over the reed beds.
- **What limits it and what it costs.** The wet margins, which the wetlands themselves need: 4.0–5.6% of them per
  billion (*screen*).

### Canopy-connected towns at forest edges

Towns at the edges of forests 100–500 m tall, whose routes rise through the forest's layers on their own supports
beside the trees, with gliders launched from the crowns and the dim lower forest kept for quiet trails
([biosphere outline](../../../ensemble/planning/Open_Moon_Ensemble_Outline_plan-01.md)).

- **Where, air, warmth, water and light.** The rain belt's interior: 4.79 million km² shared with farming country,
  5–27°, 1.5 km above the sea at the mean, 117 kPa, 20.1–25.3 °C, 2.6 mm of rain a day; practical dusk ends 87 hours
  after sunset at 19°. Over a third lies under a high Earth and a third on the far side (*screen*).
- **Hazards.** Wind on the giants: a 300 m tree's roots give way near 23.5 m/s of steady wind, or 40 m/s with a crown
  that streamlines (megaforest checkpoint, conditional). How fire spreads at 17.5% oxygen and 0.16 g is open
  (ecology).
- **The third dimension.** Walkways and stations through the canopy layers; a micro glider launched from a 300 m
  crown glides 4.5 km at 15 to 1 (*derived*).
- **What limits it and what it costs.** The forests' edges, on land the farms and crops also need. At 10,000 per km²
  a billion people take 6.0–8.3% of their share of the rain belt's interior (*screen*).

### Rain-belt farming country

Villages and market towns among the farms of the rain belt: plains of storage organs and orchards of day fruit set at
sunrise and ripe at sunset, so harvest, travel and the markets gather at dusk
([lunar-cycle ecology](../lunar_cycle_ecology/README.md)).

- **Where, air, warmth, water and light.** As canopy towns, on the same land.
- **The third dimension.** Cargo drones and sky vans between farm and market; harvest freighters to the cities.
- **What limits it and what it costs.** The rain belt's land, which the crops and the towns share. At affluent
  countries' 143–357 m² a person (Angel et al. 2010, via ecology) a billion people take 12–30% of their share of the
  interior outside a kept half (*screen*), and the crops themselves take the land of finding 2.

### Highland and ridge towns

Towns on ground 4 km and more above the sea, in cooler air: terraces on crater walls, ridges in the fog-garden country
of the middle latitudes, and the rainy upland rise round the summit.

- **Where.** 1.99 million km² shared with the terrain-held districts, 8–61°, 4.2–6.6 km above the sea; 86% on the far
  side (*screen*).
- **Air and warmth.** 110 kPa, 19.2 kPa of oxygen at the mean 5.2 km; 17.8–23.1 °C, a lunar-day mean of 18.5–19.7 °C
  (*screen*). The 3-D box over the highlands found comfortable hours rising with height, from 134 a lunar day below
  1 km to 240–243 at 4–6 km, with winds up the slopes by day and cold air draining into the valleys at night
  ([climate/crm](../../../climate/crm/README.md)).
- **Rain and water.** 2.7 mm a day; standing water 38 km away at the median.
- **Light.** At 34°, practical dusk ends 105 hours after sunset and 145 hours lie beyond it, never below 0.2 lux.
- **What limits it and what it costs.** The highland ground, a small home: a billion people at 10,000 per km² take
  14–20% of the towns' half of it outside the kept share (*screen*).

### Districts held up by towers or terrain

Districts carried along crater walls and ridges or held up on towers, where the same steel reaches two and a half to
three times as high as on Earth and sky boats reach the upper floors
([metropolis brief](../../../habitation/summit_metropolis/README.md)).

- **Where.** On the highlands with the ridge towns, and in any dense city.
- **What limits it and what it costs.** Structure. The metropolis brief's tower district, towers on a fifth of the
  ground at 30 storeys and 70 m² of floor a person, houses 85,700 people per km² (*derived*), and a billion take
  1.7–2.4% of the districts' half of the highlands outside the kept share (*screen*). Its buildings carry about 28 t
  a person (the inhabited-volume central budget); held 24 km up as on the summit tower, the frame adds 175 t of steel
  a person (*derived* from the tower form's first-order sizing). Its heat over the district is 129–388 W/m², and its
  flyers would fill NASA's lanes at a 2.7% share of trips by sky boat (*screen*).

### Nearside towns under the fixed Earth

Inland towns of the comfortable middle latitudes on the near side. The Earth stands in the same place in the sky all
month, near full at midnight, and where it stands high it keeps the night near the light at which Earth's civil
twilight ends.

- **Where.** 1.76 million km², 32–46°, 1.3 km above the sea at the mean (*screen*).
- **Air and warmth.** 118 kPa, 20.5 kPa of oxygen; 19.4–26.1 °C.
- **Rain and water.** 1.2 mm a day, 5.2 km³ of runoff a year; standing water 30 km away at the median.
- **Light.** At 39°, practical dusk ends 115 hours after sunset and twilight alone leaves 124 dark hours, from 152 at
  32° to 65 at 46°. On 59% of this land the Earth stands at least 30° up and holds the night at about 2–3 lux at its
  darkest; on the rest it stands lower and gives about 0.4 lux (main's places, the sea-appearance calendar;
  *screen*).
- **What limits it and what it costs.** The land: 8.1–11% of it per billion. A billion's domestic water is 3.5–11
  times the local runoff, drawn from the seas and lakes for 3.2–9.5 GW (*screen*).

### Limb towns where the Earth rises and sets

Towns along the limb, where the libration lifts the Earth above the horizon and lowers it again once a month.

- **Where.** 0.57 million km², 32–47° (*screen*); the eastern Smythii headland of the sea-appearance study lies on the
  same band at the equator.
- **Light.** At 42°, practical dusk ends 123 hours after sunset and 109 hours lie beyond it; twilight gives 0.9 lux at
  the darkest, to which the Earth near the horizon adds a few hundredths (main's places, at Smythii).
- **What limits it and what it costs.** The narrow band, the smallest home: 25–35% of it per billion (*screen*).

### The far side's twilit-night places

The far side's middle latitudes, with no Earth in the sky and twilit nights: dark skies and radio quiet for astronomy,
with towns lit by shielded, warm lamps, and the high platforms near 70 km above.

- **Where.** 1.75 million km², 32–47°, 2.2 km above the sea at the mean (*screen*).
- **Air and warmth.** 116 kPa, 20.3 kPa of oxygen; 19.4–25.6 °C.
- **Rain and water.** 0.7 mm a day, 2.9 km³ of runoff a year; standing water 41 km away at the median.
- **Light.** At 39°, practical dusk ends 116 hours after sunset, 122 hours lie beyond it and the darkest is 0.54 lux of
  twilight (*screen*).
- **What limits it and what it costs.** The dark night itself. A billion people's lamps light 45–81% of this home at
  the lunar design's lighting, while their towns at 143–357 m² a person cover 12–41% of it outside the kept share
  (*derived*, *screen*).

### High-latitude towns

Towns poleward of 48.6°, where a clear night never falls below the end of civil twilight and the Sun's light wraps
round the Moon all night; dry, with the polar highs such as Nobile and Peary as bases for the polar sampling
programme.

- **Where.** 5.00 million km², 51–77°, 1.7 km above the sea at the mean; 42% on the near side (*screen*).
- **Air and warmth.** 117 kPa, 20.4 kPa of oxygen; 19.0–23.1 °C.
- **Rain and water.** 0.1 mm a day and 0.39 km³ of runoff a year; standing water 33 km away at the median.
- **Light.** Practical dusk lasts the whole night; at 63°, day vision holds for 123 hours after sunset and twilight
  gives 45 lux at the darkest, 479 lux at 75° (*screen*).
- **Hazards.** If the regional magnets are built, their 0.5 mT lines, which the guidance for implanted medical
  devices keeps people outside, take in 53% of this land with the 500 km loops and 13% with the 1,000 km loops, at
  the turns that keep most heritage sites clear (*screen*, from the magnets map).
- **What limits it and what it costs.** Within the shared cases, the magnets' corridors if they are built. A billion
  people take 2.9–4.0% of the land outside the kept share; their domestic water, 47–141 times the local runoff, comes
  from the seas and lakes (*screen*).

### The undersea communities

Enclosed domed communities on the submerged heritage sites, beginning with Tranquility Base under 467 m of water,
reached by an access tower with a harbour and landing on the Sea of Tranquility (conservation D6 and D11).

- **Where.** Thirteen candidate sites, 378–4,027 m under the sea, at 7.3–66 atm outside (*screen*, from the heritage
  register).
- **What limits it and what it costs.** One dome for each site. A dome of 100 m radius with homes on half its floors
  at 70 m² each houses about 3,700 people, 48,600 across the thirteen (*derived*). At Tranquility's design depth of
  700 m the shell is 0.49 m of steel (conservation, section 5.6).

### Roaming sky towns

Buoyant towns drifting with the winds about 10 km up, moving under power only to keep clear of storms and one
another, to meet and to change course ([inhabited volume](../inhabited_volume/README.md)).

- **Air.** 100 kPa and 17.5 kPa of oxygen, 14.2 °C, changing by tenths of a degree through the lunar day as at the
  summit tower's heights (metropolis brief); winds of 2.6 m/s at the median.
- **Light.** The Sun stays up 12 hours longer at each end of the day than on the ground (*screen*). A drifting town's
  solar day runs 15–23 days at the median (finding 7).
- **Rain and water.** 137 litres a resident for each millimetre of rain on the envelope's top: 684 litres a day in the
  equatorial band, 17 poleward of 48.6°, at the ground's rain rate (*screen*). The climatology holds rain at the
  ground only. The freezing level stands near 25 km, where snow forms and melts on its way down, and storm tops reach
  55–80 km ([climate/gcm](../../../climate/gcm/README.md), [climate/crm](../../../climate/crm/README.md)); how much
  rain reaches 10 km is not yet computed.
- **Hazards.** Storms, which flash 0.7–1.5 hours after their cores form; a town moving at 5 m/s clears half a median
  storm in 40 minutes (*screen*). Ground above 8 km, 0.5% of the Moon and almost all of it on the far side's highland
  rise round the summit, stands in the way (*screen*). Hydrogen aloft needs cells kept from air and ignition.
- **Rescue.** A neighbour lies 2–3 minutes away by sky boat at five billion aloft; a canopy lowers each resident from
  10 km in 42 minutes (sky-fleet study; *screen*).
- **What limits it and what it costs.** The sky's share and the spacing in storms. A billion aloft are 173,900 towns
  with 85 Gt of districts and carriers, 15 Gt of reserve and 9.65 Gt of gas (7.46 Gt of hydrogen), shading 0.36% of
  the Moon; control at 5 m/s for 1–10% of the time takes 0.07–0.69 TW, trips down and back by sky boat for a fifth of
  residents a day 0.02 TW, and each percent a year of hydrogen lost 0.43 TW to remake (*screen*).

### Moored or streamlined sky towns

Buoyant towns that keep a place: moored by a tether over land or water, or streamlined hulls that hold a place, or
follow the Sun to keep a time of day, under power.

- **Moored.** A tether holds a sphere for no power. At the 50-year wind at 10 km raised by the summit tower's gust
  factor, 12.2 m/s, the sphere pulls 33 MN, which 3,230 t of steel cable 10 km long holds, 1.1% of the district's mass
  (*derived*). A tether from the ground into the storm layer is a tall grounded conductor, of the kind the
  infrastructure review finds starting lightning whenever a charged storm comes near.
- **Streamlined.** Finding 8: 0.089 kW a resident to hold a place, 0.04–0.5 kW to follow the Sun.
- **What limits it and what it costs.** As for roaming towns, with the tether's 0.56 Gt per billion or the hull's
  holding power in place of the control.

### High platforms near 70 km

Long-endurance platforms near 70 km for far-side astronomy and nearside earthshine photometry, with research crews.
A 100 m hull lifts 3.5 t and a 200 m hull 28 t, holding station on 115–750 kW in winds of 25–29 m/s (sky-fleet study);
their crews are counted outside the cases' residents.

### The array's spinning habitats

Spinning habitats on natural orbits beside the ring fleet. Workers commute in rotations, and the tinkerers who work at
the array's industry most of the time live there. One wheel gives Earth's gravity at 56 m radius at 4 rpm and the
Moon's at 9.3 m.

- **What limits it and what it costs.** Mass. A billion residents need 248–990 Gt of shield, 21 times the fleet's
  mass and 42 days of the build's matter stream at the full Stanford-torus shield, 20,000 km² of farm and 3 TW
  (provisioning population, after Johnson & Holbrow 1977). Commuting by rocket each week takes 18–46 kW a commuter;
  a rotation each lunar cycle 4.2–11 kW (provisioning population).

### Earth

Earth holds 8–14 billion in the cases. Each billion on a plant-rich diet farms 1.3–1.4 million km², 2.7–2.9% of
today's farmland, and up to 10.2 billion eat within four planetary boundaries (provisioning population, after Gerten
et al. 2020). Thermal plants warm Earth by 0.003–0.009 K a billion across the band.

### Elsewhere in the Solar System

Two billion in every case, taking 4–12 TW across the band, with no capacity credited.

## What every lunar resident costs

Per billion people on the Moon or aloft (*screen* unless named):

| Cost | Per billion | Basis |
|---|---|---|
| Energy | 2.0–6.0 TW | the provisioning band, 2.0, 3.68 and 5.99 kW a person |
| Heat on the Moon, industry in orbit | 1.5–4.5 TW, 0.040–0.12 W/m² | placement factor 0.755 |
| Added dimming that answers it | 0.013–0.063% | 0.009–0.014% a terawatt |
| Crops at the crop model's yields, plant-rich diet | 139,000–393,000 km² | ecology |
| Share of the rain-fed land | 2.5–16% of it; 3.5–22% outside a kept 30%; 4.9–31% outside half | ecology's 2.51–5.65 million km² |
| Grazing on the fog desert | 0–16% of it | ecology |
| Land under a light-polluted sky | 0.79–1.42 million km², 3.4–6.1% of the dry land | ecology |
| Domestic water | 18–55 km³ a year; makeup 0.18–5.5 km³ after 90–99% recovery | 50–150 litres a day |
| Crop water | 250–1,560 km³ a year; 0.04–0.27 TW if all of it is lifted from the seas | ecology's 4.85–10.9 mm a day |
| Energy through each 354-hour night | 350–1,060 TWh | provisioning population |
| Rotating floor at Earth's gravity | 420–3,330 km², 0.2–1.7 Gt | 1–8 hours a day at 10 m² an occupant |

Wild harvest from the sky's biomes is credited with no food (inhabited volume, ecology). A lunar centrifuge giving
Earth's gravity in resultant with the Moon's spins at 4 rpm on 55 m of radius with its floor tilted 80.5°, or at 2 rpm
on 220 m (*screen*). How much time at Earth's gravity a lifetime at a sixth of it needs is unknown: a centrifuge ridden
20 minutes every other day on a 16-day Shuttle mission reduced cardiovascular deconditioning, and in a lunar-mission
simulation found the Moon's gravity too weak to prevent it (Clément et al. 2015).

## The shared cases

Each case's lunar population is divided between the surface and the sky at its own split (50%) and at the
sensitivity's 0, 25 and 75%. Four surface mixes divide the surface's people among the settings, after the
metropolis's 100 million and the undersea communities' 48,600: **spread** in proportion to each home's land,
**waterfront** four times as dense on the shores and wet margins, **rain belt** four times as dense on farms, forest
edges and wet margins, and **cool ground** four times as dense on the highlands and the middle and high latitudes.
Three aerial mixes divide the sky's: all **roaming**, **half moored**, and all **streamlined** hulls holding their
places. Each mix is a rule for placing people, and every one carries every case at its own split.

At each case's own split, across the four surface mixes and three aerial mixes (*screen*):

| Case | Lunar people (surface + aloft) | Energy, all people | Heat on the Moon at 3.68 kW | Added dimming at 3.68 kW | Towns' share of dry land | Rain belt in towns and rain-fed crops, kept 30% | Lit sky, share of dry land | Towns aloft and spacing over the clear sky, kept 30% | Mass and gas aloft |
|---|---|---|---|---|---|---|---|---|---|
| total20 | 3 + 3 | 40–120 TW | 18–33 TW | 0.16–0.45% | 1.3–2.0% | 20% to full | 20–36% | 522,000, 7.1 km | 255 Gt, 29 Gt |
| total25 | 5 + 5 | 50–150 TW | 30–55 TW | 0.26–0.76% | 2.1–3.3% | 33% to full | 34–61% | 870,000, 5.5 km | 425 Gt, 48 Gt |
| total30 | 7.5 + 7.5 | 60–180 TW | 45–83 TW | 0.39–1.2% | 3.2–5.0% | 49% to full | 51–91% | 1,305,000, 4.5 km | 638 Gt, 72 Gt |
| total25_earth_heavier | 3 + 3 | 50–150 TW | 18–33 TW | 0.16–0.45% | 1.3–2.0% | 20% to full | 20–36% | 522,000, 7.1 km | 255 Gt, 29 Gt |
| total25_moon_heavier | 6 + 6 | 50–150 TW | 36–66 TW | 0.31–0.92% | 2.6–4.0% | 39% to full | 41–73% | 1,044,000, 5.0 km | 510 Gt, 58 Gt |

The heat counts the residents, the sky towns' control or holding, trips down and back, the irrigation lift of
finding 2 and the hydrogen remade after losses of 1–10% a year. "To full" means the high end of the diet range leaves
crops for irrigation once all the wet land is used; irrigated to the crop model's yield, crops and towns take
12–35% (total20), 20–59% (total25), 29–88% (total30) and 23–70% (moon heavier) of the rain belt's land outside the
kept 30%, for 0.11–3.2 TW of pumping (*screen*).

The aerial share moves the costs this way in total25 (*screen*):

| Aloft | Towns aloft | Spacing, clear sky, kept 30–50% | Towns under a median storm | Lifting gas | Hydrogen upkeep at 1–10% a year | Added dimming at 3.68 kW |
|---|---|---|---|---|---|---|
| 0% | none | – | – | – | – | 0.24–0.41% |
| 25% | 435,000 | 6.6–7.8 km | 7–10 | 24 Gt | 1.1–11 TW, 1.7–17 ppm | 0.25–0.59% |
| 50% | 870,000 | 4.7–5.5 km | 15–21 | 48 Gt | 2.1–22 TW, 3.5–35 ppm | 0.26–0.76% |
| 75% | 1,305,000 | 3.8–4.5 km | 22–31 | 72 Gt | 3.2–33 TW, 5.2–52 ppm | 0.27–0.94% |

Water drawn by region, in total25 at its own split: the domestic water of the equatorial band's people is 0.2–1.5%
of its runoff and the tropical band's 1–8%, while the middle latitudes draw 9–58% of their runoff and the high
latitudes 21–174 times theirs, the difference coming from the seas and lakes (*screen*). The full tables, by case,
share, mix and band, are in the product.

## The whole built from the parts

- **Land by band.** Towns and rain-fed crops share the rain belt; the towns' own footprint is small beside the crops',
  and the high end of the diet range needs irrigation in every case. The largest share any home's towns take outside
  a kept half is 9–29% across the mixes in total25 at its own split and 26–88% with all of total30 on the surface
  (*screen*).
- **Sky towns by band.** Spread evenly over the clear sky, five billion aloft put about 220,000 towns in each of the
  four latitude bands. Drifting with the 10 km winds they gather: 38% over the equatorial band, 22% the
  tropical, 13% the middle latitudes and 27% the high latitudes, so the equatorial towns stand 3.8–4.5 km apart and
  the middle latitudes' 6.5–7.7 km (*screen*).
- **Water by band.** Every band's domestic water closes on runoff or the fresh seas; at the ground's rain rate, sky
  towns' rain catchment exceeds their makeup in every band, narrowly poleward of 48.6° at 150 litres a day and 90%
  recovery (*screen*).

## Where this departs from the known designs

- **The metropolis** keeps the brief's 100 million at about 40,000 per km² within 32 km. The atlas gives 2,464 km² of
  dry land in that circle once the lakes are left open, so the density over dry land is 40,600 per km². Its warmth
  over the city's mean ground, 14.9–16.7 °C, is this study's estimate from the GCM's land air; the brief's 12 °C at the
  tower's foot comes from the site product and stands. Its heat over the city uses the provisioning band with industry
  in orbit (60–181 W/m²), where the brief took 2 kW a person released whole (80 W/m²). Its lit sky and its crops use
  ecology's per-person figures for the lunar design.
- **The port** keeps the summit tower study's 1.5 million regular people, counted outside the residents.
- **The flyers.** Sky boats, air taxis and their busy-hour numbers keep the sky-fleet study's, and the flyers per
  person at the busy hour are carried from the metropolis's Tokyo-like split to towns of other densities. The trip
  down from a sky town and back is the sky boat's climb from 10 km to the mean dry land at its drive efficiency,
  2.3 kWh a passenger. The streamlined town is the sky-ship study's 3 km hull, sized there for the flight band at
  40 km with 400 kg passengers; here it flies at 10 km with the hull's drag coefficient over propulsive efficiency
  (0.0357 on the volume's two-thirds power) and carries 50 t a resident under the inhabited-volume lift shares, 3,900
  residents a hull. The mooring tether takes the summit tower study's gust factor of 1.4 on the 50-year 3-day-mean
  wind at 10 km and the tower form's cable allowable.
- **The undersea dome** is the conservation study's 100 m dome; the residents on half its floors are this study's
  estimate.
- **The tower district** is the metropolis brief's: towers on a fifth of the ground averaging 30 storeys.

## Assumptions

Stated once; the product's `assumptions` block holds them all.

- **Densities.** The ecology screen's lunar design: 25 m² a person in the metropolis and 100 m² in a compact town;
  affluent countries' 143–357 m² in farming country and the far side's places; 11.7 m² in the tower district.
- **Homes.** The shore is land within one atlas cell (7.6 km) of a sea or lake of 1,000 km² or more; the wet margins
  lie beside smaller lakes and rivers of 100 m³/s or more; the highlands start 4 km above the sea; the rain belt lies
  within 30°.
- **Kept shares.** 30% and half, the Kunming–Montreal and Half-Earth levels in ecology's account, applied to every
  home and to the sky. How much of the Moon is kept wild and dark is the author's to set.
- **Mixes.** A mix places people four times as densely in its favoured settings, after the metropolis's planned 100
  million and the undersea domes.
- **Food.** Crops grow on the rain alone in the bands where most land is wet (0–45°), wettest first; beyond that they
  are irrigated at the crop model's yield, with water lifted from the seas to the mean land height at 70% pump
  efficiency.
- **Sky towns.** The inhabited-volume reference envelope at 10 km, its lift shares and its 50 t district; control at
  5 m/s for 1–10% of the time; a fifth of residents travel down and back each day by sky boat; the joint ledger's
  1–10% a year of lift gas lost, with released hydrogen held about two years by the soils; 2 km of clearance over
  terrain; the median storm's 24 km of rain for the storm count.
- **Lit sky.** Every lunar resident, aloft as on the ground, lights ecology's 791–1,420 m² of land under a
  light-polluted sky.
- **Health.** Facilities at Earth's gravity for 1–8 hours a day at 10 m² an occupant.

## Inputs to bind at integration

Two inputs live on other branches and are read at their pinned commits through Git, with their hashes checked:

- **Ecology's per-billion costs:** `biosphere/ecology/results/people.json` with
  [people_and_land.md](https://github.com/dgoldman0/Terluna/blob/11392b48c50b9af909cf8570f92e0bb50987efbe/biosphere/ecology/people_and_land.md),
  branch `domain/ecology` at `11392b4`.
- **Main's corrected Earthlight per place:** `research/studies/joint_synthesis/results/places.json`, branch `main` at
  `78020e0`; this branch's copy predates the correction.

The [zonal-winds diagnostic](../../../climate/gcm/zonal_winds.py) was built beside this study on this branch and is
read from `climate/results/gcm/zonal_winds_A28_dim5_moon.json`; the product records its hash, and a regenerated
diagnostic needs this study rerun. Without it, run.py falls back to the design run's stored global-mean wind by layer
and marks the drift numbers as a placeholder.

## Reproduction and checks

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.ways_of_living.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q research/studies/ways_of_living
```

The run reads the ignored grid products (the atlas and drainage grids and the design climatology) and the two bound
inputs; the tests that regenerate the product skip cleanly without them. The 19 tests check the product's
regeneration, the producer, input and bound hashes and the shared constants; that people across the settings sum to
every case, share and mix; the units per billion; that the homes partition the dry land; the twilight against the
register; the food fill's conservation; the hull's drag against the sphere's and the Sun-following power's cube of the
cosine of latitude; the sky towns' spacing, mass and gas scaling; the hydrogen upkeep; and the grid helpers.

What this leaves for later work: actual town and hull designs with their structure, gas cells and wet mass; the sky
towns' storm routing and traffic rules in the diagnostic's winds; the rain at their height; the lighting that meets a
dark-night requirement once its value is set; crop water in the plant model, which sets the rain-fed land within a
factor of two (ecology); and a lifetime's health at a sixth of Earth's gravity.
