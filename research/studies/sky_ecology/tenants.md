# Tenants: who lives on an aerophyte, and what it can carry

**The top is the host's income of light, and its tenants live on the flanks, the undersides, the fringes and in its
pools, where the light is far brighter than shade on Earth.** The dense air below an aerophyte at 10 km sends up 14–27%
of the daylight over ground of albedo 0.05–0.2 and 53% over cloud, so a lone body's underside gets 19–30% of the light
on a level surface at its top: 8–13 mol of photosynthetic photons per m² a day as a 24-hour mean at the equator, a half
to four-fifths of the 15.7–16.5 mol the ground gets at 60° (the surface-light product). A sky reef's undersides get
11–22%, and the rings where its modules touch 5–8%: the deep shade of the reef, where nests and domatia belong. Numbers
are *screen* values from [tenants.py](tenants.py) (`tenants` in the [product](results/sky_ecology.json)) unless marked
*derived*.

## The light on each surface

**How the light is estimated.** From the illumination domain's clear-sky [surface-light
product](../../../illumination/surface_light/README.md):
- **The downward light at 10 km** is the product's ground spectra, direct and diffuse, over the same ground, carried up
  through the air below by two-stream adding: a layer of reflectance R over ground of albedo A passes (1 − R)/(1 − AR)
  of the light arriving on it to the ground. That puts 9–10% more daylight on a level surface at 10 km than on ground
  of albedo 0.05–0.2, and 6% more over cloud. The direct beam and the diffuse light are scaled alike, which holds the
  direct share at 10 km at the ground's; giving the beam its own gain through the layer moves the zone shares by at most
  0.003 (the independent review).
- **The air below 10 km** holds 19% of the column's mass. Its reflectance comes from the column's computed reflectance
  for light from below: inverted through the conservative two-stream law R = (3/4)τ/(1 + (3/4)τ) to an effective depth
  that folds the column's absorption in, and scaled by that share. The same law on the column's Rayleigh depth matches
  the product's reflectance to a median 0.0035 at 400–560 nm and within 0.01 in 79% of those bins; where the column
  absorbs, the product lies below the law, by up to 0.068 (at 503.5 nm), and nowhere in 400–700 nm above it by more than
  0.0044. The inversion keeps that absorption.
- **The upward light at 10 km** is the downward light times the reflectance of that layer over a Lambertian ground of
  albedo 0.05–0.2, or 0.5 for cloud below (the product's albedo grid; green canopies are dark in this band and pale
  developed soils bright).
- **Diffuse light** from the sky and from below is taken as isotropic. The Sun crosses the zenith over the equatorial
  lunar day, sampled at the midpoints of 91 equal steps. Inside a raft, ray tests against the six nearest and twelve
  next modules block the direct beam and part of each surface's view.

| Zone of a body (share of the sphere) | Lone body, ground 0.05–0.2 | Over cloud | Raft module, ground 0.05–0.2 | Over cloud |
|---|---:|---:|---:|---:|
| Top cap, within 60° of the zenith | 84–85% | 89% | 79% | 79% |
| Upper flank, 60–90° | 54–59% | 70% | 27% | 27% |
| Lower flank, 90–120° (a raft's contact ring) | 32–41% | 58% | 5–8% | 15% |
| Under cap, 120–180° | 19–30% | 54% | 11–22% | 43% |

Shares are of the daytime light on a level surface at the top: 980–1,018 µmol/m²/s over the lunar day at 10 km (1,112
over cloud), 42–44 mol/m² a day as a 24-hour mean, where the product puts 38.7–40.7 mol/m² a day on the ground at the
equator. Beneath a sky reef the bright air lights a ceiling a kilometre or two across.

## Real estate

**Each projected square metre of a round body carries four of surface, one in each zone; a raft carries 3.6 per square
metre of raft, 0.91 in each, with pockets between its modules.**
- **The top cap** holds the host's green tissue. Tenants there shade it.
- **The flanks** take low Sun and the bright sky; on giants they can hold canopy soils.
- **The undersides** hold shade phototrophs, collectors and nests.
- **Fringes** of sorbent and fibre hang below, 0.5–2 m² per projected m² (a design guess; the water note's square metre
  of sorbent sits inside it).
- **Pools.** The reference organism's free water store, 5 kg/m² with a 1 kg floor, can lie as 1–5 litres of pool per m².
  A cloud forest's tank bromeliads hold up to 50,000 litres a hectare, 5 litres per m² of ground (Sugden & Robins 1979,
  as cited by later bromeliad studies): an aerophyte carries the same water as a bromeliad forest, per m².
- **A raft's pockets.** Above the ring where three modules touch, the pocket between them holds 0.40 r m³ per m² of raft
  (two triangles a cell of (√3 − π/3) r³ each): 12 m³/m² for 60-m modules and 30 m³/m² for 150-m ones. Each opens below
  through the raft's gaps, 9.3% of its area, so a tank membrane across it turns the pocket into a pool; the reef's water
  store sits there, as a bromeliad keeps its water in its leaf axils.

## What a body carries

**The pressure-trim and water figures below are static mass equivalents, and the next section gives their timing.**
The pressure allowance is shared with water management, and the starting trim state matters. Additional tenant masses
also need renewed checks of damage tolerance, attachment, respiration and water.

**The 100-m round body's spare-lift calculation allows 196 t of additional tenants, with a pressure-trim
mass equivalent of 1.9 t.** These are separate conditional limits. The aerophyte product supplies supported mass
and gas fraction; lift left within the gate is (0.7 − gas fraction) × gross lift (2/3)DΔρ. Extra residents require gas:
one kilogram of hydrogen lifts 13.4 kg at 10 km, made for 3.8 kg C by photolysis at reference light (9.5 in
canopy-like light) or 17 kg C fermented at 2.1 mol/mol.

The pressure allowance enters the trim rule Δp = ΔM p/(ρ V/A): 300 Pa in round bodies and 100 Pa in reef modules.
It is shared with the day's water swing; the usable direction and magnitude depend on the starting state. The
additional 4 kg/m² of dumpable water is also the damage reserve. Releasing it takes seconds to minutes (below), and
using it for added loads spends that reserve.

| Body at 10 km | Tenants within the reserve | Of which 10–30% of a century's surplus buys the gas | Pressure-trim mass equivalent | Trim plus water mass equivalent |
|---|---:|---:|---:|---:|
| Round, 50 m | 2.3 kg/m², 4.5 t | all | 0.24 t | 8.1 t |
| Round, 100 m | 25 kg/m², 196 t | 60–196 t | 1.9 t | 33 t |
| Round, 200 m | 68 kg/m², 2,130 t | 220–2,130 t | 15 t | 141 t |
| Round, 500 m | 179 kg/m², 35,100 t | 900–12,100 t | 236 t | 1,020 t |
| Round, 1 km, 300 MPa fibre | 427 kg/m², 335,000 t | 2,700–36,700 t | 1,890 t | 5,030 t |
| Round, 2 km, 300 MPa in canopy-like light | 777 kg/m², 2.4 Mt | 55,000–294,000 t | 15,100 t | 27,700 t |
| Reef module, 60 m | 3.4 kg/m², 10.5 t | all | 0.14 t | 12.6 t |
| Reef module, 100 m | 16 kg/m², 140 t | 43–140 t | 0.63 t | 35 t |
| Reef module, 150 m | 30 kg/m², 583 t | 33–443 t | 2.1 t | 80 t |

The century budget is a design guess: a body that spends 10% of a century's surplus on fermented gas for its tenants
and one that spends 30% on photolytic gas. Lift binds small bodies and gas carbon binds giants: filling a 1-km giant's
reserve with tenants would take 274–1,228 years of its surplus at the reference light. A mature body, holding its size,
pays for tenants from what its surplus leaves after replacing its dead at the adult hazard: a third to nine-tenths of a
1-km giant's surplus with photolytic gas, and nothing for 500-m and 1-km giants and 150-m modules on fermented gas at
the high hazard ([food web](food_web.md)).

**Larger bodies have greater static load allowances.** Dividing those allowances by individual mass gives
mass-equivalent counts; flyers arriving together are a step load, timed in the next section. For the 1.7-t
albatross-scale case, the 100-m round body has roughly one flyer equivalent by trim and 19 including water. The
corresponding counts are 9 and 82 at 200 m, 137 and 594 at 500 m, and 1,100 and 2,900 at 1 km. The 60-m reef module's
trim equivalent is below one such flyer.

## How sudden is sudden

**A load is sudden when it arrives within a few minutes: an aerophyte at fixed volume bobs on the air's density
gradient with a period of 25 minutes, whatever its size.** The author asked how sudden the 0.63–2.1 t that a reef
module's trim absorbs may be. A superpressured body keeps its volume, so a load makes it sink until the air it displaces
weighs as much as the body and the load: Δz = H_ρ ln(1 + ΔM/M), with M the mass of the air it displaces and H_ρ the
air's density scale height at 10 km, 63 km (the sky-ship product's levels). Displaced by δz it feels a restoring force
M g δz / H_ρ, against its own mass and the air it carries along, half its displaced mass for a sphere (C_a = 0.5), so it
swings with a period 2π √((1 + C_a) H_ρ / g) = 25 minutes (*screen*). A 2-km raft carries a broad disc of air with it,
an added mass of at most (8/3)ρR³ (a solid disc's; the raft's 9.3% gaps and its thickness lower it), and swings in up
to 67–102 minutes. Quadratic drag damps the swing, with C_d 0.2–0.5 for a sphere (a design range) and 1.2 for a raft
moving broadside; the deepest point of a step comes from integrating the float with that drag.

| At 10 km, a load equal to the trim mass equivalent | Load | Sinks, left alone | A step overshoots by | Nine-tenths of the swing gone in | The gas warms in | Venting ballonet air | Dropping the water |
|---|---:|---:|---:|---:|---:|---:|---:|
| Round, 50 m, 300 Pa spare | 0.24 t | 190 m | 42–66% | 1.5–3.8 h | 0.2–1.0 h | 1.2–12 min | 1–35 s |
| Round, 100 m | 1.9 t | 190 m | 61–80% | 3.0–7.5 h | 0.4–2.0 h | 2.5–25 min | 2–70 s |
| Round, 200 m | 15 t | 190 m | 76–89% | 6–15 h | 0.8–4.0 h | 5–50 min | 4–141 s |
| Round, 500 m | 236 t | 190 m | 89–95% | 15–38 h | 2.0–10 h | 12–124 min | 11–352 s |
| Round, 1 km, 300 MPa fibre | 1,890 t | 190 m | 94–97% | 30–75 h | 4.0–20 h | 25–249 min | 22–703 s |
| Reef module alone, 60 m, 100 Pa spare | 0.14 t | 63 m | 72–87% | 5–12 h | 0.2–1.2 h | 0.8–8 min | 0.4–13 s |
| Reef module alone, 100 m | 0.63 t | 63 m | 81–92% | 8–20 h | 0.4–2.0 h | 1.3–13 min | 0.7–21 s |
| Reef module alone, 150 m | 2.1 t | 63 m | 87–94% | 12–31 h | 0.6–3.0 h | 2–20 min | 1–32 s |

**Left alone, a body takes a load by sinking, and the sink costs more superpressure than the trim allows until its gas
warms.** The trim equivalent is the load whose share of the body's mass equals the spare superpressure's share of the
air's pressure, so every body sinks the same 190 m, or 63 m for a 100-Pa module, to take it. On the way down the ambient
pressure rises by ρ g Δz, 1.24 times the spare, because in air that cools upward density falls more slowly with height
than pressure (H_ρ of 63 km against a pressure scale height of 51 km). As the gas warms to the lower air the cost falls
back: to the trim rule's own amount with the lapse the density profile implies, and to 1.06 times it with the sky-ship
temperatures. The gas, heated at fixed volume, warms over 0.2–20 hours, larger bodies more slowly (an exchange through
the envelope of 2–10 W/m²/K, a design range). A body that does nothing thus keeps its spare superpressure through a load
spread over half an hour up to 81% of its trim equivalent, and through one spread over hours up to 94–100%.

**A step needs 1.6–2.0 times the margin of a gradual load.** A load that lands within a few minutes is a step: the body
falls through its new float height and overshoots by 42–97% of the sink, most in giants, whose drag is small against
their mass, and the swing dies away over 1.5 hours to 3 days. The deepest point comes about 13 minutes in, before the
gas warms, so a body that does nothing keeps its superpressure through a step only up to 41–50% of its trim equivalent.

**Dropping water answers a step; venting answers most of it.** Taking the load at the same height means shedding its
mass: venting that much ballonet air through valves at the spare superpressure, or dropping water. Venting within
minutes expands the gas adiabatically, so it costs 1.4 times the trim rule's superpressure and takes at most 71% of the
trim equivalent at once, the rest as the gas warms again. With valves of 0.001–0.01% of the projected area (a design
guess), venting the whole equivalent takes minutes for reef modules and small round bodies and up to four hours for a
1-km giant, whose valves would have to grow with its column. Drains of the same share under 0.1–1 m of head (a design
guess) drop the trim equivalent's water in seconds to minutes and cost no superpressure, spending the damage reserve.

**In a raft the ties share the load.** A module's 0.63–2.1 t spread through a 2-km raft sinks it 17–39 cm and takes
under 1% of its superpressure; a module carries its own trim equivalent alone only where its ties are slack.

## Domatia

**A reef can grow modules for its tenants, as a siphonophore specialises its zooids: a 60-m domatium carries 11 kg/m²,
three times a green module's 3.4, and the colony feeds it from two to six green modules.** A domatium leaves out the
green top (it keeps a living layer of 0.1 kg/m² that maintains its skin, fixes 0.2 kg C/m² a year and makes its gas by
photolysis, a design guess), carries the round bodies' 300 Pa of spare pressure and a full water store, and takes its
upkeep through the ties from the green modules around it, paid from 10–30% of their surplus (a design guess).

| Module | Tenants within the reserve | Trim alone, plus water mass equivalents | Upkeep the colony pays | Green modules a domatium needs |
|---|---:|---:|---:|---:|
| 60 m domatium | 11.0 kg/m², 34 t (green: 3.4, 10.5 t) | 0.41 t, 12.9 t | 0.43–0.53 kg C/m²/yr | 1.6–6.0 |
| 100 m domatium | 21.9 kg/m², 189 t (green: 16.1, 140 t) | 1.9 t, 36.5 t | 0.79–0.89 | 4.2–14 |
| 150 m domatium | 33.0 kg/m², 644 t (green: 29.9, 583 t) | 6.4 t, 84 t | 1.32–1.42 | 20–66 |

The upkeep's low end gives the guards no food and the high end 0.1 kg C/m² a year of food bodies. Small domatia are
cheap and carry most; 100–150 m landing domatia have mass equivalents of one to four albatross-scale flyers by trim and
21–49 including water. Their place at the reef's edge is proposed; how flyers land there remains to be worked out.

## Guilds and mutualisms

| Guild | Where | Earth precedent | Light | Weight | Damage | Nutrients |
|---|---|---|---|---|---|---|
| Guards | domatia, flanks | Macaranga triloba puts about 5% of its daily above-ground production, 9% of construction costs, into food bodies for its ants (Heil et al. 1997); ants give Dischidia major 29% of its nitrogen and 39% of its leaf carbon (Treseder et al. 1995) | none | 1–10 g/m² | keep grazers and egg-layers off the green top | debris and respiration in the domatia |
| Roosting and nesting flyers | giants' flanks, landing domatia | seabird colonies receive 99 Gg of P a year (Otero et al. 2018); rat-free islands carry 760 times the seabirds and 251 times the nitrogen, and their reefs 48% more fish (Graham et al. 2018); guano falls at 0.32–0.37 g/m² a day within 50–320 m of Spitsbergen colonies (Zwolicki et al. 2013) | droppings foul the top unless roosts are on the flanks | a 1.7-t flyer at a time on a 100-m body | claws and nests wear the skin; loads arrive at once | an albatross-scale flyer excretes 78 kg of P a year, 30–60% of it at its roost |
| Pool food webs | crevice tanks and catchments | tank bromeliads hold up to 5 L per m² of a cloud forest (Sugden & Robins 1979); a jumping spider supplies 18% of its bromeliad's nitrogen and lengthens its leaves 15% (Romero et al. 2006) | none | the host's own 1–5 kg/m² of water | larvae and pathogens in still water | decomposers mineralise litter and frass and the host takes up the phosphate: the recycling the sky needs |
| Canopy soils | giants' flanks | rain-forest trees root into the soils on their own branches (Nadkarni 1981); epiphytes weigh 1.3–6.0 kg dry per m² of central branch (Freiberg & Freiberg 2000) | the flanks' own | 2.6–12 kg dry per projected m², more when wet: within the spare lift of round bodies from 150 m up | weight and water on the skin | hold litter and water; the host roots into them |
| Drifting-weed community | fringes, undersides | Sargassum carries epiphytes, fungi, more than 100 invertebrate and over 100 fish species and four sea turtles (Coston-Clements et al. 1991); its nitrogen-fixing epiphytes may be the weed's largest source of reactive nitrogen (Johnson et al. 2023, after Phlips et al. 1986); 1,205 species raft (Thiel & Gutow 2005) | 19–30% (lone) and 11–22% (reef) of the top's | 0.05–0.5 kg/m² | epiphytes shade what they cover | nitrogen fixers feed the host |
| Grazers | the green top | tropical canopies lose 12–19% of their leaf production (Metcalfe et al. 2014); forests and shrublands 7.9% of above-ground production (Cebrian 2004) | remove green tissue | | cost 3.8–6.0% of a reef module's surplus | frass carries leaf P before resorption: caught in pools it returns, falling it is lost |
| Collectors | fringes, tethers | orb-weavers eat the pollen their webs catch (Eggs & Sanders 2013) | none | 1–10 g/m² | none | turn airborne pollen, spores and insects into frass on the host |
| Bud dispersers | visitors | the nursery partnership of the aerophyte growth note | none | | none | carry buds to the crowns and are paid in food |

Weights other than the flyers', the water and the epiphytes' are design guesses. The guard entry is collective mass per
footprint area, live or dry unspecified; individual size, the guard mass that defence requires and the most a reef can
carry are open. The [guards and symbioses analysis](symbioses.md) works with explicit live-mass scenarios,
complementary guilds, local response and shared food budgets. The reference organism allocates 0.1 kg C per m² a year
to all tenants, counted once across guilds. This is 12–18% of a 60-m reef module's net production (its biomass and the
sugar together), two and a half to three and a half times Macaranga's share.

## The canopy nursery inoculates the young

**A juvenile leaves its megaforest crown with lift to spare for the canopy's life.** The 50-m, 60-m and 80-m round
bodies have 2.3, 6.9 and 16 kg/m² of lift within their reserve, against a founder community of 0.01–0.1 kg dry per m²,
0.04–0.4 kg wet (a design guess), gathered over its one to seven years in the crown. Epiphyte faunas can double a
rainforest canopy's invertebrate biomass (Ellwood & Foster 2004), and 1,205 species are known to cross the sea on
floating objects (Thiel & Gutow 2005): the crown's guards, decomposers, pool-dwellers and epiphytes ride out with the
young as rafted life rides between coasts.

## Contact and disease

**Passive aerophytes at one height share their air and gather onto shared tracks, so a height band is one contact
group.** At 1% cover a 100-m body passes within 50 m of a neighbour 3.3–33 times a day at 0.1–1 m/s of relative drift,
and a 2-km reef 0.11–1.1 times (density × speed × swept width; the speeds and the 50-m margin are design guesses). The
zonal-wind product releases its routes 10° of latitude apart (300 km) and 30° of hour apart (910 km at the equator). At
10 km their median nearest neighbour is 125 km away after 30 days and 41 km after half a year; from then on it swings
between 0.2 and 92 km from one sample to the next as tracks merge and part, and 66–87% of routes share a 171-km model
cell with another. At 20 km the median settles at 20–35 km after the first year, with 91–94% sharing a cell. A 10-µm
spore falls at 0.55 mm/s in the lunar air (the ecology first screen): it stays aloft about 105 days in a 5-km mixed
layer of surface air, and one shed at 10 km takes about 210 days to reach the ground (*derived*), so it travels the
whole track with the bodies. The 1983 die-off of the sea urchin Diadema antillarum reached 3.5 million km² of the
Caribbean in 13 months, its agent carried by surface currents (Lessios 1988).

**Heights are separate contact groups.** At 10 km the winds gather bodies most strongly over the north pole (3.4 times a
uniform spread) and toward the equator (1.3–1.8 times); at 20 km over both poles (8.6–12 times). Populations at the two
heights rarely share air, and sailing, 30–70° of latitude in a lunar day (the aerophytes' sailing note), takes a body
off a shared track. Resistance that varies between heights and lineages keeps a pathogen in one band.

## Checks

[test_tenants.py](test_tenants.py) checks:
- the layer reflectance's and the underside factor's limits;
- the two-stream law against the product at every wavelength: never more than 0.005 below the product in 400–700 nm, a
  median difference under 0.005 at 400–560 nm, and the recorded check fields;
- the height factor against two-stream adding, and the light at 10 km above the ground's by under 15%;
- the zone areas, and a lone sphere's sky view factors against (1 + cos β)/2;
- that neighbours only take light away;
- the pocket volume against a Monte Carlo integral;
- the capacities against the aerophyte product's rows, and the domatia against green modules;
- the vertical response: the static sink against the density scale height, the superpressure it takes, the period
  with added mass, the drag damping of a sphere, and the order of venting and dropping times;
- the spores' fall from 10 km, the encounter rate and the landing counts.
