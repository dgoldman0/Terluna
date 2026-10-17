# The ring fleet's glow and its answers

The ring fleet's films absorb a few percent of the sunlight they pass and radiate that heat from both faces, and
the faces that look at the Moon send it part of the heat. The planned climate is the investigated design case near
295 K, which leaves this glow out. The author's decisions frame the answer
([register](../../decisions.md#solar-shield-and-habitat-array)):
- the glow is answered with added solar protection, by sending it away from the Moon, or both;
- added dimming is the easy part, since more of the light can be taken and turned into power, and the glow itself
  is the problem to work out;
- the fleet keeps its light off the Moon's night zones, turning the night half's light away and using it where it
  can.

[glow.py](glow.py) follows the fleet's heat onto the Moon from the shield study's own products and screens the
answers. It writes [results/glow.json](results/glow.json) in about two CPU minutes. Numbers marked *screen* are
computed in it. *Estimate* marks a result that rests on a stated approximation, and *design guess* a value chosen
for a coating, an absorber or the Moon's light. Sources, with how far each was read, are in
[glow_sources.json](glow_sources.json).

## What the screen finds

1. **As built, the stack puts 14.1 W/m² on the Moon, 1.4–2.5 times the first screen's glow** (*screen*).
   - **The glow.** The films absorb 67.7 PW, at the low end of array_heat's 65–83 PW, and their own emission
     gives the Moon 9.4 W/m²: 9.7 W/m² over the day hemisphere and 9.0 over the night hemisphere.
   - **Why it is higher.** The layout's radii are apolune radii, at the sunward crossing, and each ring's night
     half rides at perilune, 15,400–16,000 km out. The Moon's view, weighted by the glow, is 0.0092 against the
     first screen's 0.0066–0.0081. The films are partly transparent in the thermal infrared, so along its normal
     the Moon sees into the stack: 0.55–0.62 of a stack's heat leaves toward it at every depth.
   - **What the films send back.** Silica reflects strongly in its reststrahlen bands, 43% at 8.5–9.6 µm and
     30–32% at 20–25 µm. Weighted by the Moon's outgoing infrared, one layer reflects 19–21% of it at normal
     incidence and a deep stack 22–25%. A tile that faces the Moon images the Moon back onto itself. The tiles are
     held perpendicular to their velocity with a 1.5° shingle, 1.5° off the Moon's centre at the sunward crossing and
     up to 11° near the nodes, where the forced eccentricity tilts them, so 34% of what they reflect lands. That
     returns 4.0 W/m² of the Moon's own infrared and 0.71 W/m² of its reflected sunlight.
   - **In proportion.** The 14.1 W/m² is 4.8% of the design sunlight. The Moon absorbs all of the infrared and 72%
     of the sunlight, and the GCM answers the infrared at 0.71 of its watts, so dimming answers it at
     4.8% of the design sunlight (6.6% at an efficacy of 1).
2. **Heat mirrors on the Moon-facing faces raise the Moon's heat** (*screen* on *design guesses*).
   - On every tile, a face that reflects the thermal infrared cuts the glow to 5.2–8.9 W/m² at the candidates'
     favourable ends, and 2.6 W/m² on the nearest layer alone. It reflects the Moon's own infrared as well, and
     returns 11–19 W/m² of it.
   - Every candidate raises the heat on the Moon, by 6.7–21 W/m², and the dimming that answers it rises from
     4.8% to 7.0–11.8%.
   - The conducting oxides and silver stacks also absorb 3–10% of the sunlight, doubling or tripling the films'
     heat, and reflect the solar near infrared that the window has to pass.
   - Behind a mirror on every tile the hidden layer runs at 227–485 K at the sunward crossing.
3. **Rolling outside the service arc halves the heat, and the attitude study keeps it inside the clearance**
   (*screen*; the attitude study's laws along the kept orbits).
   - Rolled 80° outside the service arc on both halves, the attitude study's law, the tiles meet the Sun and the
     Moon obliquely there. On the layout re-flown for its push (ring_layout_rolled), the glow falls to 3.6 W/m²
     (6.7 by day, 0.48 by night), the returned infrared to 1.75 W/m² and the returned sunlight to 0.61 W/m²: 6.0
     W/m² in all. The forced eccentricity falls to 0.053–0.094 and the night halves ride at 17,060–17,720 km.
     Scaling the as-built eccentricity by the law's drive had given 5.7 W/m²: the re-flown stack sits 350–500 km
     further in, and its edges keep 0.58–0.62 of their eccentricity against the drive's 0.48.
   - Rolled on the night half alone, the fleet leaves 8.7 W/m² with its own eccentricity: the smaller eccentricity
     turns the day half's tiles closer to the Moon's centre, and more of the Moon's own infrared comes back.
   - Both laws keep the stack's 222 m between shingles. Rolling both halves cuts the drive on the eccentricity to
     0.48 and takes a momentum store of 1.0–2.3×10¹⁰ N·m·s a tile, 5–12 times the chosen one. The night half alone
     cuts it to 0.55 and takes 2.7–9.0×10¹⁰ N·m·s, 14–48 times, since its roll torque does not cancel over an orbit.
     The design's centred filter leaves both stores about as they are.
   - Tilting the tiles 3° or 6° toward the Sun through the nodes, where the roll passes through radial-facing,
     cuts the night's edge light by a third to a half and the momentum store to 1.2–1.7×10¹⁰ N·m·s, with the push
     unchanged. It breaks the clearance: with the roll ramping, neighbouring shingles meet over 3–7% of the orbit.
   - Tilts within the clearance's pitch envelope, −3° to +6°, change the glow by 0.1–0.4% and move the returns
     about: less infrared comes back (3.6 W/m²) and more of the day side's light reaches the night.
4. **What the films absorb is the titania's ultraviolet absorption, and absorbing less costs light pressure**
   (*screen*; the options *design guesses*).
   - Between 242 and 380 nm the titania absorbs 2.74 of the annulus film's 3.03 points of sunlight and 4.00 of the
     window stack's 4.50.
   - That band is what the window has to stop: behind the titania stack no light below about 340 nm reaches the
     air, which keeps it free of ozone and the upper air cold, as the loss chain's titania cases assume. The
     annulus holds O1's UV transmission below 175 nm, a band that carries 0.1 W/m².
   - A wide-gap absorber in the titania's place on the annulus, with its edge at 220–250 nm, cuts the glow to
     3.7–3.9 W/m² if the annulus may pass 250–400 nm along its paths above the limb. It takes two-thirds of the
     annulus film's light pressure (2R + A from 0.146 to 0.049), and one sail loading then needs about 5.5% of the
     sunlight reflected to restore it.
   - An ultraviolet mirror ahead of the titania keeps the band stopped by reflecting it: 5.7–5.8 W/m² with mirrors
     on the outer faces and 2.2 W/m² on both. Each face weighs 11 g/m² on the annulus and 6.8 g/m² on the window,
     80% of an annulus tile.
   - Thinner titania (0.05 µm) cuts the glow by 9% and gives up a level of ultraviolet margin.
5. **The night** (*screen*).
   - The glow brings the night no light. As built it brings the night hemisphere 9.0 W/m² of heat, 93% of the
     day's, and 13.5 W/m² with the returns.
   - The films' specular reflection of the Moon's day side lands on the night as well: 1.6 lux over the night
     hemisphere, 1.0 lux at midnight and 2.6 lux within 10° of the terminator, a thousand times a far-side
     midnight's own light. The joint synthesis counted diffuse scatter and glints of sunlight; this is a third source.
   - Rolling the night half keeps its reflections off the Moon. On the layout's eccentricity 0.08 lux is left, all
     within 20° of the terminator, from day-half tiles near the nodes. With the smaller eccentricity a roll's push
     leaves, those tiles pass radial-facing at the nodes closer to the Moon's centre and light the night's edge: on
     the re-flown layout 0.8 lux over the night hemisphere, 3.5 lux within 10° of the terminator, 1.0 lux at 10–20°
     and 0.12 at 20–30°, and none within 50° of midnight.
   - **That light falls in the Moon's own twilight**, which the tall sky keeps bright long after sunset: the solved
     sky gives 1,200–4,600 lux with the Sun 0–10° below the horizon, 180–1,200 lux at 10–20° and 24–180 lux at
     20–30° (illumination/sky). The returns add 0.3–0.6% of it at each band's darker edge, and the night decision's
     deep night stays dark. The tilt through the nodes that would cut it breaks the clearance (finding 3), and the
     night does not need it.
6. **The recommendation is to roll both halves outside the service arc and dim by extraction for the rest**
   (*screen*).
   - On the layout re-flown for its push, it leaves a glow of 3.6 W/m² (6.7 by day, 0.48 by night) and returns of
     2.4 W/m², 6.0 W/m² in all. Added dimming of 2.0% of the design sunlight answers it at the GCM's infrared
     efficacy (2.7% at an efficacy of 1) and gives back 45–67 TW of electricity.
   - It moves the attitude outside the service arc and the dimmer. The films, tile masses, sail loading,
     service-arc attitude and the protection's spectra stay as set. It takes a momentum store 5–12 times the chosen
     one, and the layout re-flown for half the push, 350–500 km further in, which keeps the design's 1 km between
     strips inside the photon roll's steering reach.
   - The annulus absorber (3.8 W/m² with this roll on the re-flown layout) and the ultraviolet mirror on the outer
     faces (3.3 W/m²) would cut it further, each with the costs above. The author set both aside on 13 October, the
     plants' water design making up the rain the glow takes
     ([register](../../decisions.md#the-gas-pair-13-october-2026)).

## The films: what they absorb and why

The screen evaluates both films with the protection domain's layers, optical constants and thin-film recursion,
taken to oblique incidence with s and p averaged. At normal incidence it reproduces annulus_film.json's
absorptance and transmittance and array_heat's face emissivities exactly. Absorbed sunlight, in points of the solar
constant, at normal incidence on the outer face (*screen*):

| Band | Annulus film, 0.1 µm titania on 4 µm silica | Window stack, 1 µm on 10 µm |
|---|---:|---:|
| Below 200 nm | 0.007 | 0.007 |
| 200–242 nm | 0.109 | 0.109 |
| 242–320 nm | 1.554 | 1.579 |
| 320–340 nm | 0.773 | 0.914 |
| 340–360 nm | 0.336 | 1.020 |
| 360–380 nm | 0.077 | 0.488 |
| 380–400 nm | 0.016 | 0.101 |
| Visible | 0.005 | 0.034 |
| Beyond 2.5 µm | 0.110 | 0.161 |
| **All** | **3.03** | **4.50** |

- **The solar infrared** beyond 2.5 µm carries 45.4 W/m² and adds 0.3% to array_heat's absorbed sunlight, which
  applied the 0.1–2,500 nm absorptance to the whole solar constant (*estimate*: a 5,772 K blackbody's shape).
- **Lit on the Moon-facing face**, as the night half is, the films absorb 3.21% and 4.77%. Up to about 70° of
  incidence the absorptance stays within 7% of these values. At 80° it falls to 2.7% and 4.0%, where the films
  reflect 47%.
- **What the protection needs**, as the protection domain states it:
  - O1, the UV transmission below 175 nm, within the loss budget's allowance: the annulus film passes 1.0×10⁻⁵ of
    that band, X-rays included, and the window stack 6.7×10⁻⁷. The swarm's gaps at its standard level, 2×10⁻⁴,
    set what reaches the air.
  - O3, the cut-off: the titania stack, with no ozone and a surface UV index of 0.1, or an ozone-forming edge near
    200 nm, which the author has not chosen. The climate design case and the loss chain are computed behind the
    titania stack; behind the 200-nm edge the oxygen atoms alone add 0.6–6 kg/s.
  - E2 and O4: the window passes the design sunlight, 5% less than the titania stack at every wavelength.
  - The annulus is a light film that blocks the ultraviolet and passes visible light, under the author's zoned
    aperture. Its light crosses the tall air above the limb.

So the window's absorption between 242 and 380 nm is its protection function. On the annulus the loss budget
counts O1's band, and the decision's wording covers the ultraviolet as a whole.

| Option (*design guess*) | Annulus absorbs | Window absorbs | Mass added, annulus / window | Annulus 2R + A | Glow |
|---|---:|---:|---:|---:|---:|
| As built | 3.03% | 4.50% | 0 | 0.146 | 9.37 W/m² |
| Annulus titania 0.05 µm | 2.73% | 4.50% | −0.2 g/m² | 0.162 | 8.56 W/m² |
| Annulus absorber, 220 nm edge | 0.14% | 4.50% | +1.1 g/m² | 0.048 | 3.71 W/m² |
| Annulus absorber, 250 nm edge | 0.25% | 4.50% | +1.1 g/m² | 0.049 | 3.94 W/m² |
| Ultraviolet mirror on the outer faces, 95% | 0.39% | 2.10% | +8.2 / +5.1 g/m² | 0.224 | 5.81 W/m² |
| Ultraviolet mirror on the outer faces, 98% | 0.30% | 2.03% | +11.0 / +6.8 g/m² | 0.226 | 5.71 W/m² |
| Ultraviolet mirror on both faces, 98% | 0.30% | 2.03% | +22.1 / +13.7 g/m² | 0.226 | 2.19 W/m² |

- **Thinner titania** passes 2.6×10⁻⁶ of the band between 2.5 and 175 nm. It meets the 10⁻⁵ level and misses the
  10⁻⁶ level that 0.1 µm meets. It reflects more, so its returned sunlight rises (0.85 W/m²), and its pressure rises
  by 11%: one sail loading would then need 1.7 g/m² more on each annulus tile.
- **A wide-gap absorber** (hafnia's edge near 220 nm, silicon nitride's or zirconia's near 250 nm) stops less of
  the sunlight. The screen takes it as a step absorber that holds the film's present opacity below 175 nm, which
  0.1–0.2 µm of a dense oxide would have to show. As built, the annulus film passes 58% of 320–400 nm and 1.4% of
  242–320 nm; with the absorber it would pass nearly all of 250–400 nm along its slant paths through the tall air
  above the limb, and the 220-nm edge part of O2's Herzberg band as well. The annulus tiles' sail loading falls to a
  third of the window rings', and nested rings of different loading are forced to eccentricities that cross.
  Restoring it takes 5.4–5.5% of the sunlight reflected by the annulus film, for instance by a near-infrared
  reflector band, with mass and rejected light of its own.
- **An ultraviolet mirror** of hafnia and silica reflects 245–390 nm on the annulus and 245–340 nm on the window,
  where the titania already passes nothing, so the window keeps its spectrum. Measured mirrors reach about 99% near
  250 nm and up to 98% across 245–400 nm. The mass comes from quarter-wave stacks: seven pairs a stack, three stacks
  on the annulus and two on the window, 2.0 and 1.2 µm thick. At 98% on the outer faces the added reflection keeps
  the two loadings matched within 0.1 g/m²; at 95% the annulus tiles need 1.3 g/m² more. The fleet grows by about
  60% in mass. A mirror on the outer face works only by day: lit on its Moon-facing face, the night half still
  absorbs in the titania. The concerns are stress in 2 µm of oxide on a 4 µm film, ultraviolet-induced absorption in
  hafnia, thermal cycling every orbit, and the reflected ultraviolet, which O5 and S7 have to trace.

## The stack: N layers and heat mirrors

At each point of the sky, the layers over it form a stack of parallel, specular plates, partly transparent in the
thermal infrared. Each plate has its faces' reflectance and its transmittance in 63 bands and 9 directions, from the
films' optics. The plates share the absorbed sunlight there evenly. The model solves every plate's temperature from
its energy balance and gives the radiance along the normal toward the Moon and the stack's normal reflectance. Its
energy closes to 10⁻¹³. With two grey opaque layers it reproduces the provisioning screen's two-layer estimate
exactly: every share, cut and hottest layer (*screen*).

Stacks of the as-built films at the sky's mean absorbed sunlight, 56 W/m² (annulus) and 81 W/m² (window):

| Layers | Toward the Moon, annulus / window | Inward over the hemisphere | Hottest layer | Moon's infrared reflected |
|---:|---:|---:|---:|---:|
| 1 | 0.554 / 0.617 | 0.533 / 0.589 | 187 / 189 K | 0.188 / 0.210 |
| 2 | 0.562 / 0.587 | 0.530 / 0.556 | 175 / 188 K | 0.215 / 0.220 |
| 3 | 0.565 / 0.574 | 0.526 / 0.541 | 177 / 198 K | 0.227 / 0.222 |
| 8 | 0.568 / 0.560 | 0.514 / 0.518 | 187 / 229 K | 0.242 / 0.224 |
| 32 | 0.575 / 0.567 | 0.504 / 0.505 | 235 / 296 K | 0.246 / 0.226 |

- **Over the hemisphere** a deep stack sends half its heat inward, as the first screen took it. Along the
  Moon's direction the films are clearer, so the Moon sees the stack's warmer interior and receives 0.55–0.62 at
  every depth.
- **Temperatures.** At the sunward crossing a single tile runs at 175 K (annulus) and 177 K (window), and two
  overlapping tiles at 191–207 K. The deep stacks near the nodes share the same sunlight among many layers. In 32
  of them the hottest annulus layer reaches 235 K, and the window stack, more opaque in the infrared, about 296 K.
- **Layers in view.** The screen counts the sky's layers as along-track tiles over a tile's angular length
  times rings over a strip's angular width: 2.6 on average in the Moon's view, against the first screen's 2.03,
  because the night half is closer.

Heat mirrors and emitters on every tile, at each candidate's favourable (low) and other (high) end (*design
guesses*, anchored to measured coatings):

| Coating | Emissivity | Sunlight absorbed | Taken from the beam | Mass | Glow | Infrared back | Total | Hidden layer |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Mesh, a few µm pitch | 0.10–0.30 | 0.5–2% | 5–20% | 0.02–0.2 g/m² | 5.41–7.90 | 13.3–17.1 | 23.9–24.4 | 239–296 K |
| Metasurface | 0.15–0.40 | 0.5–2% | 5–20% | 0.02–0.2 g/m² | 5.93–8.60 | 11.4–16.2 | 23.1–23.5 | 227–273 K |
| Conducting oxide | 0.11–0.21 | 3–10% | 5–15% | 0.5–1.5 g/m² | 8.91–18.7 | 15.0–16.9 | 26.8–35.0 | 304–333 K |
| Silver stack | 0.02–0.10 | 4–10% | 15–30% | 0.1–0.3 g/m² | 7.17–12.7 | 17.1–18.6 | 28.0–33.0 | 377–485 K |
| Strong outer emitter alone | 0.85 | 0.2% | 0 | 2 g/m² | 8.80 | 4.16 | 13.67 | 178–195 K |
| Mesh and outer emitter | 0.10, 0.85 | 0.7% | 5% | 2.0 g/m² | 5.17 | 17.11 | 23.65 | 269–293 K |
| Mesh on the nearest layer only, emitter | 0.10, 0.85 | 0.7% | 5% | 2.0 g/m² | 2.62 | 17.11 | 20.82 | 218–254 K |

Glow, infrared back and total are in W/m² over the Moon. The hidden layer is the hotter of two overlapping tiles at
the sunward crossing.

- **Why they fail.** A heat mirror facing the Moon reflects the Moon's own infrared along with the glow, and the
  tile's attitude images it back. The first screen's cuts (1.5–2 times on every tile, 3–9 times on the nearest
  layer) hold for the glow, and the returned infrared outweighs them for every candidate.
- **The nearest layer.** A mirror on the layer nearest the Moon alone is a bound. Every tile passes through every
  depth of the stack round its orbit, so a fixed coating cannot sit on the nearest layer only.
- **Light taken from the beam.** A coating that takes light out of the window's beam dims the Moon by that share,
  unextracted. The high ends take more than the heat calls for (the negative "to add" values below): they are
  dimmers that waste the light they take.
- **The candidates themselves.** The measured anchors are Fan and Bachner's TiO2/Ag/TiO2 heat mirror (84% passed at
  0.5 µm, 98–99% reflected at 10 µm), Woo and colleagues' nano-mesh (over 80% passed at 550 nm, emissivity under
  0.2 at 7–14 µm) and ITO films of emittance 0.11–0.21. Silver and conducting oxides reflect beyond about 1 µm,
  which takes the solar near infrared out of the window's light.

## Attitude

The dynamics hold each tile perpendicular to its velocity with the 1.5° shingle. On the forced eccentricity
(0.11–0.16) that puts the normal 1.5° off the Moon's centre at the sunward crossing, up to 6.1° within the service
arc and 10.9° near the nodes. From the tiles the Moon is a disc 4.7–6.5° in radius, and a tile's reflection of it
lands on it only where the normal lies within about that angle of its centre (*screen*).

A tilt by an angle cuts a face's view of the Moon by its cosine: 0.5% at 6°, 13% at 30°, 50% at 60° and 83% at 80°.
A rolled tile also takes less sunlight, with the cosine of its incidence, so a rolled night half's glow falls nearly
as the cosine squared. The bundle's clearance allows pitches of −3° to 0° and 3° to 6° at every phase. Outside
the service arc it allows rolls of up to 65–90°, one side on the day half and the other on the night half, passing
through radial-facing at the nodes.

| Attitude | Glow | Infrared back | Sunlight back | Total | Night lux | Its law in the attitude study |
|---|---:|---:|---:|---:|---:|---|
| As built | 9.37 | 4.04 | 0.71 | 14.12 | 1.64 | |
| Pitch +6° | 9.36 | 3.59 | 0.75 | 13.70 | 7.95 | |
| Pitch −3° | 9.34 | 3.59 | 0.75 | 13.68 | 7.93 | |
| Night roll 30° | 8.15 | 1.63 | 0.67 | 10.45 | 0.08 | not run |
| Night roll 60° | 5.73 | 1.63 | 0.67 | 8.03 | 0.08 | not run |
| Night roll 80° | 4.93 | 1.63 | 0.67 | 7.24 | 0.08 | push 0.55, store 9.0×10¹⁰ N·m·s |
| Night roll 80°, its own eccentricity | 4.72 | 2.83 | 1.12 | 8.67 | 0.78 | perilune 17,180–18,060 km |
| Rolled outside service 30° | 7.86 | 1.41 | 0.54 | 9.81 | 0.08 | push 0.86, store 5.8×10¹⁰ N·m·s |
| Rolled outside service 60° | 4.87 | 1.38 | 0.53 | 6.78 | 0.08 | push 0.57, store 6.2×10¹⁰ N·m·s |
| Rolled outside service 80° | 3.76 | 1.37 | 0.52 | 5.66 | 0.08 | push 0.48, store 2.3×10¹⁰ N·m·s |
| Rolled outside service 80°, its own eccentricity | 3.46 | 1.69 | 0.58 | 5.73 | 0.88 | perilune 17,440–18,410 km |
| Rolled 80°, tilted 3° through the nodes, its own eccentricity | 3.48 | 1.69 | 0.57 | 5.74 | 0.51 | store 1.7×10¹⁰ N·m·s; gap 0 m |
| Rolled 80°, tilted 6° through the nodes, its own eccentricity | 3.51 | 1.69 | 0.57 | 5.77 | 0.51 | store 1.2×10¹⁰ N·m·s; gap 0.05 m |
| Rolled outside service 80°, on the re-flown layout | 3.60 | 1.75 | 0.61 | 5.95 | 0.80 | perilune 17,060–17,720 km |

The last column is the attitude study's figure for each law along the kept orbits at 62.7 g/m²: the drive on the
eccentricity against the as-built tiles', and the median momentum store at 25% trim against the chosen
1.9×10⁹ N·m·s; at 5–26 g/m² the stores are 1.0–1.6×10¹⁰ (both halves) and 2.7–5.1×10¹⁰ N·m·s (night half). Turning
a tile against gravity gradient and the orbital rate takes 0.2–2 MN·m at 62.7 g/m². The rows with their own
eccentricity scale every ring's forced eccentricity by the law's drive about the sunward crossing's radius, the
attitude study's first-order rule. The re-flown layout (ring_layout_rolled) flies the calibration with the law's push
and lays the stack out again, so its rings carry their own radii and eccentricities. Both roll laws keep the 222 m
between shingles at every phase; the tilts through the nodes close it, the smallest gap in the last column.

## The answers screened

All in W/m² averaged over the Moon. Dimming is the added share of the design sunlight that answers the total, with
the infrared at the GCM's efficacy of 0.71 (below). The coating's share of it is light the coating takes out of the
window's beam; "to add" is the rest, which extraction can supply, with the electricity it gives at 20–30%
(*screen*).

| Answer | Glow | Infrared back | Sunlight back | Total | Dimming | To add | Electricity |
|---|---:|---:|---:|---:|---:|---:|---:|
| As built `as_built` | 9.37 | 4.04 | 0.71 | 14.12 | 4.75% | 4.75% | 106–159 TW |
| Annulus titania 0.05 µm `thin_titania` | 8.56 | 4.08 | 0.85 | 13.49 | 4.54% | 4.54% | 101–152 TW |
| Annulus absorber, 220 nm edge `wide_gap_220` | 3.71 | 4.04 | 0.63 | 8.38 | 2.82% | 2.82% | 63–94 TW |
| Annulus absorber, 250 nm edge `wide_gap_250` | 3.94 | 4.04 | 0.63 | 8.61 | 2.90% | 2.90% | 64–97 TW |
| UV mirror, outer faces, 95% `uv_mirror_outer_95` | 5.81 | 4.04 | 0.71 | 10.55 | 3.55% | 3.55% | 79–119 TW |
| UV mirror, outer faces, 98% `uv_mirror_outer_98` | 5.71 | 4.04 | 0.71 | 10.45 | 3.52% | 3.52% | 78–117 TW |
| UV mirror, both faces, 98% `uv_mirror_both_98` | 2.19 | 4.04 | 0.71 | 6.94 | 2.34% | 2.34% | 52–78 TW |
| Mesh heat mirror, low end `mesh_low` | 5.41 | 17.11 | 1.41 | 23.93 | 8.06% | 3.06% | 68–102 TW |
| Mesh heat mirror, high end `mesh_high` | 7.90 | 13.31 | 3.13 | 24.35 | 8.21% | −11.79% | none |
| Metasurface, low end `metasurface_low` | 5.93 | 16.16 | 1.41 | 23.50 | 7.91% | 2.91% | 65–97 TW |
| Metasurface, high end `metasurface_high` | 8.60 | 11.41 | 3.13 | 23.14 | 7.80% | −12.20% | none |
| Conducting oxide, low end `tco_low` | 8.91 | 16.92 | 1.00 | 26.83 | 9.03% | 4.03% | 90–135 TW |
| Conducting oxide, high end `tco_high` | 18.66 | 15.02 | 1.35 | 35.03 | 11.79% | −3.21% | none |
| Silver stack, low end `silver_low` | 7.17 | 18.64 | 2.24 | 28.04 | 9.45% | −5.55% | none |
| Silver stack, high end `silver_high` | 12.70 | 17.11 | 3.14 | 32.95 | 11.10% | −18.90% | none |
| Strong outer emitter `emitter` | 8.80 | 4.16 | 0.71 | 13.67 | 4.60% | 4.60% | 102–154 TW |
| Mesh and outer emitter `mesh_low_emitter` | 5.17 | 17.11 | 1.37 | 23.65 | 7.97% | 2.97% | 66–99 TW |
| Mesh on the nearest layer, emitter `mesh_low_nearest_emitter` | 2.62 | 17.11 | 1.09 | 20.82 | 7.01% | 2.01% | 45–67 TW |
| Pitch +6° `pitch_plus_6` | 9.36 | 3.59 | 0.75 | 13.70 | 4.61% | 4.61% | 103–154 TW |
| Pitch −3° `pitch_minus_3` | 9.34 | 3.59 | 0.75 | 13.68 | 4.61% | 4.61% | 103–154 TW |
| Night roll 30° `night_roll_30` | 8.15 | 1.63 | 0.67 | 10.45 | 3.52% | 3.52% | 78–117 TW |
| Night roll 60° `night_roll_60` | 5.73 | 1.63 | 0.67 | 8.03 | 2.71% | 2.71% | 60–90 TW |
| Night roll 80° `night_roll_80` | 4.93 | 1.63 | 0.67 | 7.24 | 2.44% | 2.44% | 54–81 TW |
| Rolled outside service 30° `rolled_outside_service_30` | 7.86 | 1.41 | 0.54 | 9.81 | 3.30% | 3.30% | 74–110 TW |
| Rolled outside service 60° `rolled_outside_service_60` | 4.87 | 1.38 | 0.53 | 6.78 | 2.28% | 2.28% | 51–76 TW |
| Rolled outside service 80° `rolled_outside_service_80` | 3.76 | 1.37 | 0.52 | 5.66 | 1.91% | 1.91% | 42–64 TW |
| Night roll, thin titania `night_roll_80_thin_titania` | 4.53 | 1.65 | 0.80 | 6.97 | 2.35% | 2.35% | 52–78 TW |
| Night roll, absorber 220 nm `night_roll_80_wide_gap_220` | 1.99 | 1.63 | 0.59 | 4.21 | 1.42% | 1.42% | 32–47 TW |
| Night roll, absorber 250 nm `night_roll_80_wide_gap_250` | 2.10 | 1.63 | 0.59 | 4.33 | 1.46% | 1.46% | 32–49 TW |
| Night roll, absorber, window mirror `night_roll_80_wide_gap_250_window_mirror` | 1.15 | 1.63 | 0.59 | 3.37 | 1.14% | 1.14% | 25–38 TW |
| Night roll, UV mirror outer 98% `night_roll_80_uv_mirror_outer_98` | 1.28 | 1.63 | 0.67 | 3.58 | 1.21% | 1.21% | 27–40 TW |
| Night roll, emitter `night_roll_80_emitter` | 4.62 | 1.69 | 0.67 | 6.98 | 2.35% | 2.35% | 52–79 TW |
| Night roll, mesh and emitter `night_roll_80_mesh_low_emitter` | 2.83 | 6.95 | 1.30 | 11.07 | 3.73% | −1.27% | none |
| Rolled outside service, absorber 220 nm `rolled_outside_service_80_wide_gap_220` | 1.49 | 1.37 | 0.46 | 3.32 | 1.12% | 1.12% | 25–37 TW |
| Rolled outside service, absorber 250 nm `rolled_outside_service_80_wide_gap_250` | 1.58 | 1.37 | 0.46 | 3.41 | 1.15% | 1.15% | 26–38 TW |
| Rolled outside service, mesh, emitter `rolled_outside_service_80_mesh_low_emitter` | 1.98 | 5.87 | 1.02 | 8.88 | 2.99% | −2.01% | none |
| Rolled outside service, UV mirror outer 98% `rolled_outside_service_80_uv_mirror_outer_98` | 1.00 | 1.37 | 0.52 | 2.90 | 0.98% | 0.98% | 22–33 TW |
| Night roll 80°, its own eccentricity `night_roll_80_own_eccentricity` | 4.72 | 2.83 | 1.12 | 8.67 | 2.92% | 2.92% | 65–98 TW |
| Rolled outside service 80°, its own eccentricity `rolled_outside_service_80_own_eccentricity` | 3.46 | 1.69 | 0.58 | 5.73 | 1.93% | 1.93% | 43–64 TW |
| Rolled 80°, tilted 3° through the nodes `rolled_outside_service_80_node_3_own_eccentricity` | 3.48 | 1.69 | 0.57 | 5.74 | 1.93% | 1.93% | 43–64 TW |
| Rolled 80°, tilted 6° through the nodes `rolled_outside_service_80_node_6_own_eccentricity` | 3.51 | 1.69 | 0.57 | 5.77 | 1.94% | 1.94% | 43–65 TW |
| Rolled outside service 80°, re-flown layout `rolled_outside_service_80_rolled_layout` | 3.60 | 1.75 | 0.61 | 5.95 | 2.01% | 2.01% | 45–67 TW |
| Re-flown layout, absorber 250 nm `rolled_outside_service_80_wide_gap_250_rolled_layout` | 1.52 | 1.75 | 0.53 | 3.80 | 1.28% | 1.28% | 28–43 TW |
| Re-flown layout, UV mirror outer 98% `rolled_outside_service_80_uv_mirror_outer_98_rolled_layout` | 0.96 | 1.75 | 0.61 | 3.32 | 1.12% | 1.12% | 25–37 TW |

## By lunar day and night

W/m² over each hemisphere; night lux is the films' specular reflection of the day side, averaged over the night
hemisphere (*screen*):

| Answer | Glow by day | Glow by night | Infrared back by day | Infrared back by night | Night, all | Night lux |
|---|---:|---:|---:|---:|---:|---:|
| As built `as_built` | 9.71 | 9.04 | 3.69 | 4.38 | 13.46 | 1.64 |
| Annulus absorber, 250 nm edge `wide_gap_250` | 4.19 | 3.69 | 3.69 | 4.38 | 8.11 | 2.02 |
| UV mirror, both faces, 98% `uv_mirror_both_98` | 2.37 | 2.01 | 3.69 | 4.38 | 6.44 | 1.64 |
| Mesh heat mirror, low end `mesh_low` | 5.79 | 5.03 | 15.61 | 18.61 | 23.72 | 5.88 |
| Pitch +6° `pitch_plus_6` | 9.70 | 9.02 | 3.65 | 3.53 | 12.73 | 7.95 |
| Night roll 30° `night_roll_30` | 9.50 | 6.80 | 3.26 | 0.01 | 6.81 | 0.08 |
| Night roll 80° `night_roll_80` | 9.01 | 0.86 | 3.26 | 0.01 | 0.87 | 0.08 |
| Rolled outside service 80° `rolled_outside_service_80` | 7.05 | 0.47 | 2.74 | 0.01 | 0.48 | 0.08 |
| Night roll 80°, its own eccentricity `night_roll_80_own_eccentricity` | 8.54 | 0.91 | 5.63 | 0.03 | 0.95 | 0.78 |
| Rolled outside service 80°, its own eccentricity `rolled_outside_service_80_own_eccentricity` | 6.45 | 0.47 | 3.35 | 0.04 | 0.53 | 0.88 |
| Rolled 80°, tilted 3° through the nodes `rolled_outside_service_80_node_3_own_eccentricity` | 6.47 | 0.49 | 3.33 | 0.04 | 0.54 | 0.51 |
| Rolled outside service 80°, re-flown layout `rolled_outside_service_80_rolled_layout` | 6.72 | 0.48 | 3.46 | 0.04 | 0.53 | 0.80 |
| Night roll, absorber 250 nm `night_roll_80_wide_gap_250` | 3.91 | 0.30 | 3.26 | 0.01 | 0.31 | 0.11 |
| Night roll, UV mirror outer 98% `night_roll_80_uv_mirror_outer_98` | 2.25 | 0.31 | 3.26 | 0.01 | 0.32 | 0.08 |

- **As built** the night hemisphere takes 91% as much of the fleet's heat as the day hemisphere. The night half
  absorbs less sunlight (26 of 68 PW), and rides closer to the Moon.
- **The returned day side** reaches the night as 1.0 lux at midnight, rising to 2.6 lux within 10° of the
  terminator. The far side's own midnight is 0.001 lux, and the joint synthesis's diffuse glow at O5's 0.1% is
  25–66 lux.
- **Rolled 80°**, the night half absorbs 4.2 PW, and its own glow comes to 0.12 W/m² averaged over the Moon. The
  night's 0.86 W/m² comes mostly from day-half tiles near the terminator, and its light from their reflections
  within 20° of it.
- **The night decision** keeps the fleet's light off the night zones, and the roll it calls for halves the heat the
  stack puts on the Moon without adding mass. What the roll leaves lands within 30° of the terminator, where the
  tall sky's own twilight is 180–1,500 times brighter. Pitch tilts send five times as much of the day side's light
  onto the night, and heat mirrors 2.1–9.4 times as much.

## What remains: dimming and the power it gives back

- **The conversion.** The Moon absorbs all of the infrared the stack sends it and 72% of the sunlight (planetary
  albedo 0.280), so each W/m² of infrared takes 1/0.72 = 1.39 times as much top-of-air sunlight away at an
  efficacy of 1. The infrared heats the upper air, and in the GCM's held glow branch
  ([gases, glow and dimming](../gases_glow_dimming/README.md), 13 October) the dimming that holds the design
  temperature is 0.71 of its watts matched as absorbed sunlight, so each W/m² of infrared takes 0.99 times as much
  top-of-air sunlight away, and the domain synthesis's convention, its W/m² over the design sunlight's 293.3 W/m²,
  holds for the glow to about 1%. At the design sunlight, 1% of dimming answers 2.97 W/m² of the infrared. On the
  corrected radiation of 17 October the pair gave about 0.59; its rerun on the final radiation waits on the
  target climate.
- **The electricity.** 1% of the design sunlight over the Moon's disc carries 111 TW of light, which collectors at
  20–30% turn into 22–33 TW. That is the synthesis's 23–35 TW per 1% of the titania stack's light, taken on the
  design sunlight.
- **The dimmer's form** decides its own costs (*screen*):
  - Absorbing in the window rings' films, it absorbs round their whole orbits. As built it gives back 25% of its
    dimming as glow, 14% with the night roll and 10% with both halves rolled.
  - Reflecting from the window tiles, it adds 2ΔR to their light pressure: the recommended answer's 2.0% more
    dimming raises the window's coefficient from 0.276 to about 0.32, and one sail loading moves with it.
  - Collectors that take light only while they cross the window, and reject their heat from faces that do not see
    the Moon, keep both small.

## The recommendation

**Roll both halves 80° outside the service arc and answer the rest by extraction.**

- **What it leaves** (*screen*), on the layout re-flown for its push: the glow 3.6 W/m² (6.7 by day, 0.48 by
  night), the returned infrared 1.75 W/m² and the returned sunlight 0.61 W/m², 6.0 W/m² in all. The deep night
  stays dark. The night's edge takes 0.1–3.5 lux within 30° of the terminator from day-half tiles passing
  radial-facing at the nodes, 0.3–0.6% of the twilight the tall sky already gives there.
- **The dimming** is 2.0% of the design sunlight at the GCM's efficacy, about 1% below the synthesis's convention
  (2.7% at an efficacy of 1). Taken by extraction it gives back 45–67 TW of electricity.
- **Its costs.** No film, mass or sail-loading change. The attitude outside the service arc carries the roll. Along
  the kept orbits the attitude study finds it inside the clearance, the drive on the eccentricity cut to 0.48, and a
  momentum store of 1.0–2.3×10¹⁰ N·m·s a tile, 5–12 times the chosen one: with rotors 400 m in radius at 1,000 m/s,
  0.7–1.6 Gt across the fleet against its 46 Gt. Re-flown for the smaller push, the stack's radii move 350–500 km in
  and its forced eccentricity falls to 0.053–0.094; the design's 1 km between strips stays inside the photon roll's
  steering reach. The night half rolled alone cuts the drive less (0.55), takes a store 14–48 times the chosen one and
  leaves 8.7 W/m².
- **What it needs proven:**
  - the rolled law's own out-of-plane forces flown with the layout (its strip turn and the keeping's plane drift,
    0.03 km a day of crossing height), and the steering through the full stack;
  - the films' infrared reflectance, measured on the coated film;
  - the films' flatness, since slopes of a few degrees spread the specular image;
  - the glow and the returns as an infrared forcing in the GCM's final radiation, whose run waits on the target
    climate (the held branch of 13 October set the efficacy at 0.71);
  - a form of the added dimming that neither absorbs in the window films nor reflects from the window tiles, and
    where its light goes (O5, S7).
- **The knobs it moves least:** the films and their protection, the tile masses and the one sail loading, the
  service arc's attitude, and the climate's design sunlight, which the dimming restores. It moves the attitude
  outside the service arc, which the night decision moves on the night half anyway, the momentum store, and the
  dimmer, from 5% to about 6.9% of the titania stack's light on the 5% dimmer's design case.
- **Set aside by the author** on 13 October, with the plants' water design making up the rain the glow takes:
  - the annulus absorber, which would bring the total to 3.8 W/m² and the dimming to 1.3%. It would need the
    atmosphere domain's check of 250–400 nm above the limb, the absorber's measured edge and opacity, and the
    annulus's light pressure restored;
  - the ultraviolet mirror on the outer faces, which would bring the total to 3.3 W/m² and the dimming to 1.1%, at
    about 60% more fleet mass.

## Assumptions

- **Films** (sourced): the comparison's films, 0.1 µm of titania on 4 µm of silica on the annulus rings and the
  window stack on the window rings, both with the protection domain's anti-reflection and matching layers.
- **Optics** (*estimate*): smooth specular films, s and p averaged, with the protection domain's constants. The
  solar infrared beyond 2.5 µm follows a 5,772 K blackbody's shape scaled to close the solar constant.
- **The stack** (*estimate*): at each point of the sky, parallel infinite plates sharing the local absorbed
  sunlight evenly, with no conduction between them. Thermal bands run 2–125 µm, and the values at 125 µm hold
  beyond.
- **Rings and tiles** (sourced): every ring of the layout on the flown matched rings' apolune radius and forced
  eccentricity by height. The inventory's 28.2 million tiles are shared evenly and spaced evenly in time, and held
  perpendicular to the velocity with the 1.5° shingle. Keeping tilts are left out.
- **Window rings** (sourced): the middle 26% of the stack by height, array_heat's share.
- **Sunlight on each tile** (*estimate*): attenuated by the films the projected counts put in front of it, by
  fleet_light's rule, with a linear penumbra.
- **Layers in the Moon's view** (*estimate*): along-track tiles over a tile's angular length times rings over a
  strip's angular width. A rolled stack opens to 1 + (n − 1)|cos roll| layers (*design guess*).
- **The Moon's outgoing infrared** (*estimate*): the atmosphere domain's correlated-k model on the clear
  titania-stack columns at 288 and 298 K, interpolated to the design case's 294.9 K, gives 227 W/m² clear. Its
  spectrum is scaled to the GCM design case's absorbed sunlight, (1 − 0.280) × 293.3 = 211.2 W/m², taken as
  uniform and Lambertian (*design guess*).
- **The Moon's reflected sunlight** (*design guess*): a Lambertian day side with the planetary albedo 0.28, flat in
  wavelength.
- **Infrared efficacy 0.71** (*GCM run*): the held glow branch of the gases, glow and dimming study, 13 October; about
  0.59 on the corrected radiation of 17 October. The answers at an efficacy of 1 are in the product beside them.
- **Wide-gap absorber, ultraviolet mirror, heat mirrors and emitters** (*design guesses*): as stated in the tables,
  anchored to the measured coatings in glow_sources.json.
- **Rolls** (sourced): attitude_schemes' law, zero in the 25° service arc and at the nodes, with 30° cosine ramps.
  The night roll keeps the night half only.

## What the records say that the screen changes

- **Provisioning finding 1 and the register's row of 9 October** give the Moon 5.7–9.9 W/m² of glow. The glow is
  9.4 W/m², and the films return 4.0 W/m² of the Moon's own infrared and 0.71 W/m² of its reflected sunlight.
- **The same records** have heat mirrors on the Moon-facing faces cutting the Moon's share 1.5–2 times, or 3–9 times
  on the nearest layer. Those cuts hold for the glow. A mirror facing the Moon images the Moon onto itself, and every
  candidate raises the Moon's heat by 6.7–21 W/m².
- **The provisioning screen** has a deep stack sending half its heat toward the Moon. Along the Moon's direction the
  stack sends it 0.55–0.62.
- **The provisioning screen** places the tiles 19,355–21,392 km out. Those are the layout's apolune radii, and the
  night half rides at 15,375–16,021 km.
- **The domain synthesis** answers heat with its W/m² over 293.3 W/m². For infrared at an efficacy of 1 it takes
  1.39 times that; at the GCM's 0.71, 0.99 times, so the convention holds for the glow to about 1%.
- **The joint synthesis's night light** counts diffuse scatter and glints of sunlight. The films' specular
  reflection of the Moon's day side adds 1.6 lux over the night hemisphere as built.
- **illumination/fleet_light** builds the rings with the layout's radius as the semi-major axis, which puts their
  sunward crossings 11–16% further out than the layout's.
- **The ring fleet's inventory** counts 28.2 million tiles on rings of the apolune radius. The eccentric rings are
  11% shorter, so the tiles overlap 1.51 times along the track at apolune against the designed 1.18. A count matched
  to the designed overlap would hold about 11% fewer tiles, with the glow in proportion.

## Limits

- **The returns rest on flat films held exactly as the dynamics hold them.** Keeping tilts of about 2°, the films'
  own slopes and their wrinkles move the specular image. Moved past the Moon's 5° they take the returns away, and
  nearer its centre they add to them.
- **The infrared spectrum** is the clear-sky columns', scaled. Clouds, at 23.5% cover, emit more evenly across the
  spectrum and take a little from the window band where silica reflects most.
- **The stack's plates** share their sunlight evenly. Near the nodes, sunlight at grazing incidence crosses
  neighbouring stacks, and the plane-parallel picture holds only roughly there.
- **The returns** are traced on every eighth ring and every second sample, 41,940 tiles; the glow uses every tile.
- **The window's share** is array_heat's 26% of the stack. Covering the window through the planes' common shift
  could take up to about 40%, which would raise the films' mean absorptance by about 6%.
- **The rolled stacks' opening** and the night-only roll's push and store are design guesses. The coatings, the
  absorber and the mirror are design guesses with measured anchors.
- **No climate model is run here.** The dimming takes the GCM's efficacy from the gases, glow and dimming study, whose
  runs carry the glow as infrared added at the top of the air.

## Next work

1. The rolled law's own out-of-plane forces. The layout re-flown for the law's push (ring_layout_rolled, 13 October)
   scales the radial-facing tile's sail force by the drive; the roll's out-of-plane pushes, which turn the strips
   by about 0.001° a day and drift the crossing heights by 0.03 km a day along the kept orbits, are next in the
   flights, with the keeping that holds them. The tilt through the nodes is set aside: it breaks the clearance, and
   the light it would take off the night's edge is under 1% of the twilight there.
2. In the atmosphere domain: the Moon's outgoing infrared with clouds. The GCM carried the glow and the returns as
   infrared at the top of its air on 13 October, night side included, and set the efficacy at 0.71
   ([gases, glow and dimming](../gases_glow_dimming/README.md)); the rerun on its final radiation waits on the target
   climate.
3. Only if the author returns to the absorber or the mirror: whether the annulus may pass 250–400 nm along its paths
   above the limb, a wide-gap absorber's measured edge and opacity below 175 nm at 1–2 g/m², and the ultraviolet
   mirror's measured broadband reflectance and stress on 4 µm of silica.
4. The films' measured infrared reflectance and their flatness at the scale of a tile.
5. In illumination: the films' specular reflection of the day side added to the night-light model.
6. The added dimming's form: collectors that take light only in the window, where they ride, and where their heat
   and the light they turn away go.
7. The tile inventory against the eccentric rings' length.

## Running

```sh
OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.solar_shield_array.glow    # about two CPU minutes
nice -n 10 python -m pytest research/studies/solar_shield_array/test_glow.py -q         # about three
```

The screen needs the protection domain's optical inputs (`python protection/fetch_inputs.py`), the WHI spectrum
(`python -m atmosphere.middle_atmosphere.fetch_inputs --download`), the middle atmosphere's committed columns and
the correlated-k tables (`atmosphere/radiative_convective`, built once). It runs no dynamics. Peak memory is about
550 MB on one core.
