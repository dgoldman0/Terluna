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
  of the light arriving on it to the ground. That puts 6–10% more daylight on a level surface at 10 km than on the
  ground, most over dark ground. The direct beam and the diffuse light are scaled alike, which holds the direct share
  at 10 km at the ground's; the beam's own gain through the layer would raise it somewhat.
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
  lunar day, sampled at the midpoints of 91 equal steps. Inside a raft, ray tests against the six nearest and twelve next
  modules block the direct beam and part of each surface's view.

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
- **Pools.** The reference organism's free water store, 5 kg/m² with a 1 kg floor, can lie as 1–5 litres of pool
  per m². A cloud forest's tank bromeliads hold up to 50,000 litres a hectare, 5 litres per m² of ground (Sugden & Robins
  1979, as cited by later bromeliad studies): an aerophyte carries the same water as a bromeliad forest, per m².
- **A raft's pockets.** Above the ring where three modules touch, the pocket between them holds 0.40 r m³ per m² of raft
  (two triangles a cell of (√3 − π/3) r³ each): 12 m³/m² for 60-m modules and 30 m³/m² for 150-m ones. Each opens below
  through the raft's gaps, 9.3% of its area, so a tank membrane across it turns the pocket into a pool; the reef's water
  store sits there, as a bromeliad keeps its water in its leaf axils.

## What a body carries

**A 100-m body carries 196 t of tenants within its 30% reserve, as the 150–200 t estimate had it, and at most 1.9 t can
land at once.** The aerophyte product gives each passing body's supported mass and gas fraction; the lift left within
the gate is (0.7 − gas fraction) × the gross lift (2/3)DΔρ. Tenants weigh like anything else, so a resident load needs
gas: one kilogram of hydrogen lifts 13.4 kg at 10 km, made for 3.8 kg C by photolysis at the reference light (9.5 in
canopy-like light) or 17 kg C fermented at 2.1 mol/mol. A sudden load is what the hull's spare superpressure absorbs by
the trim rule, Δp = ΔM p/(ρ V/A), 300 Pa in round bodies and 100 Pa in reef modules; past that the body drops free
water, 4 kg/m² above its floor, and must later drink it back. Both are upper bounds. The spare pressure is the allowance
the aerophyte water note gives the day's water swing, so all of it is free only at the low point of the swing; the
dumpable water is the reserve whose size sets how many gas cells a body needs against damage (the aerophyte README), and
dropping it for a flock spends that reserve. Departures run the trim the other way.

| Body at 10 km | Tenants within the reserve | Of which 10–30% of a century's surplus buys the gas | At once by trim, at most | At once dropping water, at most |
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
has the whole of its surplus beyond replacing its dead for such uses.

**Heavy flocks belong on the big ones.** By trim alone a 100-m body takes at most one albatross-scale flyer (1.7 t) at a
time; dropping water it takes 19, two of the Argentavis scale (15 t) and no Quetzalcoatlus-scale flyer (44 t). A 200-m
body takes 9 albatross-scale flyers by trim and 82 dropping water; a 500-m body 137 and 594; a 1-km giant 1,100 and
2,900, or 43 and 114 of the Quetzalcoatlus scale. A reef's 60-m module takes no great flyer by trim at all.

## Domatia

**A reef can grow modules for its tenants, as a siphonophore specialises its zooids: a 60-m domatium carries 11 kg/m²,
three times a green module's 3.4, and the colony feeds it from two to six green modules.** A domatium leaves out the green
top (it keeps a living layer of 0.1 kg/m² that maintains its skin, fixes 0.2 kg C/m² a year and makes its gas by
photolysis, a design guess), carries the round bodies' 300 Pa of spare pressure and a full water store, and takes its
upkeep through the ties from the green modules around it, paid from 10–30% of their surplus (a design guess).

| Module | Tenants within the reserve | At once by trim, dropping water (at most) | Upkeep the colony pays | Green modules a domatium needs |
|---|---:|---:|---:|---:|
| 60 m domatium | 11.0 kg/m², 34 t (green: 3.4, 10.5 t) | 0.41 t, 12.9 t | 0.43–0.53 kg C/m²/yr | 1.6–6.0 |
| 100 m domatium | 21.9 kg/m², 189 t (green: 16.1, 140 t) | 1.9 t, 36.5 t | 0.79–0.89 | 4.2–14 |
| 150 m domatium | 33.0 kg/m², 644 t (green: 29.9, 583 t) | 6.4 t, 84 t | 1.32–1.42 | 20–66 |

The upkeep's low end gives the guards no food and the high end 0.1 kg C/m² a year of food bodies. Small domatia are
cheap and carry most; landing domatia of 100–150 m take at most one to four albatross-scale flyers by trim and 21–49
dropping water, and belong at the reef's edge.

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

Weights other than the flyers', the water and the epiphytes' are design guesses. The guards' food is the reference
organism's 0.1 kg C per m² a year of food to tenants, 12–18% of a 60-m reef module's net production (its biomass and
the sugar together), two and a half to three and a half times Macaranga's share.

## The canopy nursery inoculates the young

**A juvenile leaves its megaforest crown with lift to spare for the canopy's life.** The 50-m, 60-m and 80-m round bodies
have 2.3, 6.9 and 16 kg/m² of lift within their reserve, against a founder community of 0.01–0.1 kg dry per m²,
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
layer of surface air, and one shed at 10 km takes about 210 days to reach the ground (*derived*), so it travels the whole
track with the bodies. The 1983 die-off of the sea urchin Diadema antillarum reached 3.5 million km² of the Caribbean in
13 months, its agent carried by surface currents (Lessios 1988).

**Heights are separate contact groups.** At 10 km the winds gather bodies most strongly over the north pole (3.4 times a
uniform spread) and toward the equator (1.3–1.8 times); at 20 km over both poles (8.6–12 times). Populations at the two
heights rarely share air, and sailing, 30–70° of latitude in a lunar day (the aerophytes' sailing note), takes a body off
a shared track. Resistance that varies between heights and lineages keeps a pathogen in one band.

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
- the spores' fall from 10 km, the encounter rate and the landing counts.
