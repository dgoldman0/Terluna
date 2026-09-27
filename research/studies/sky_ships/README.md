# How large a sky ship the air can carry

The first step toward the distribution of the Open Moon's flyers, from micro gliders to the largest ships: how large a
ship this atmosphere can reasonably sustain. The study answers with light calculations:
- the air ships fly in, by height, with its winds and storms;
- the rule of similarity that lunar gravity sets for every flyer;
- how a buoyant ship's structure grows with its size, on the Moon and on Earth, calibrated against the weight
  statements of four rigid airships;
- the largest envelope a fabric holds, for non-rigid ships;
- how a wing's structure grows with its span;
- what else grows with a ship: people, hydrogen, power, mooring loads, turning, berths at the summit port.

```
OPENBLAS_NUM_THREADS=1 python -m research.studies.sky_ships.run   # results/sky_ships.json, about a second
```

The models are in [engineering/flight](../../../engineering/flight/README.md), the runner is [run.py](run.py), the
numbers are in [results/sky_ships.json](results/sky_ships.json) and the sources in [sources.json](sources.json). The air
and winds over the whole Moon come from a new climate product,
[climate/gcm/global_winds.py](../../../climate/gcm/global_winds.py) (the design run A28_dim5, years 15–24); the storms
and gusts from the cloud-resolving equatorial ring ([climate/crm](../../../climate/crm/README.md)); the port's traffic
from the [summit tower study](../summit_tower/README.md#the-central-port).

**Evidence.** First-order sizing of the kind used to compare concepts.
- **What is tested.** The models reproduce Lamb's apparent-mass coefficients, Woodward's gust moment for the
  Shenandoah and the elliptic wing's bending. Fitted to the weight statements of four rigid airships, the hull model
  gives each one's empty weight within 7%; fitted to ten transports' wing groups, the wing model gives nine of them
  within 12%. The rule of similarity holds in both models exactly.
- **What is carried over.** The loads follow Earth's airship and aircraft practice: Woodward's gust envelope, a
  deflated gas cell, the FAA's envelope rules, Pratt's gust formula. Lunar gravity enters through the loads that
  depend on weight.
- **What is assumed.** The modern materials, the mass carried per passenger, the ships' speeds and the static bending
  coefficient (0.016 of weight times length with a cell deflated, from the lift one cell of sixteen leaves). The
  ships are sizing cases; designing them is later work.

## The answer

- **The Moon's sky is Earth's, six times larger.** Lunar gravity is a sixth of Earth's. From 9 km up, the Moon's air
  has at each height the density of Earth's at a sixth of the height above 9 km, so the flight band (35–45 km) holds
  the air of Earth's 4.3–5.8 km. Convection stretches the same way: the ring's storms are 30 km wide at the median,
  with updrafts of 12–14 m/s. A flyer scaled up six times in every length, flying at the same speed through air of the
  same density, carries the same stresses from its own weight and from the same gusts: the same ship, six times
  longer and 220 times heavier.
- **Rigid ships of up to about 3 km are reasonable.** Calibrated against the weight statements of four rigid airships
  (it gives each one's empty weight within 7%), the model puts the best size for the Hindenburg's own 1930s
  technology at about 270 m on Earth, beside the great rigid airships' 239–245 m. In the flight band the same
  technology is at its best near 750 m and keeps the Hindenburg's share of lift for payload, fuel and crew (45%) out
  to 3.0 km, where a ship's turning circle is as wide as a median storm. With a carbon-fibre frame and modern fabrics
  a ship keeps 81–86% of its lift for them from 245 m to 3 km, at its best near 870 m.
- **The air's hard limit is far larger.** The structure would take the whole lift at 7 km with the Hindenburg's
  materials and 30 km with modern ones; on Earth, at 1.2 km and 4.9 km.
- **Non-rigid ships are held by their fabric.** With today's strongest hull laminates, the largest envelope that
  flies at 25–36 m/s in the flight band is 280–350 m long, twice Earth's 140–160 m.
- **Winged ships.** A wing carries its weight with the An-225's share at a span of about 710 m on the Moon, eight
  times the An-225's: about 250,000 t at the An-225's wing loading, cruising at 195 m/s in the band.
- **Within that, a ship's size is set by what it carries.** A rigid ship in the band carries about 2,000 people at
  500 m, 7,000 at 750 m, 17,000 at 1 km and 57,000 at 1.5 km (at 400 kg each with their cabins and services), and
  holds 5, 15, 37 and 123 times the Hindenburg's hydrogen. The central port's long-haul traffic, 32,000–53,000
  passengers an hour, fills 16–26 berths of 500 m ships, 5–8 of 750 m ships, or one of 1.5 km ships.

## The air ships fly in

The air over the whole Moon, from the design run (area-weighted over ten model years; heights above sea level):

| Height | Air (kg/m³): mean, 1st–99th percentile | Pressure | Temperature | Lift of hydrogen / helium (kg per m³, 98% pure) | Earth's height with this air |
|---|---|---|---|---|---|
| 0 km | 1.406 (1.382–1.433) | 1.18 atm | 20 °C | 1.28 / 1.19 | 1.5 km below sea level |
| 10 km | 1.205 (1.192–1.216) | 0.99 atm | 14 °C | 1.10 / 1.02 | 0.2 km |
| 20 km | 1.024 (1.013–1.034) | 0.81 atm | 5 °C | 0.93 / 0.86 | 1.8 km |
| 25 km | 0.942 (0.931–0.952) | 0.73 atm | 1 °C | 0.86 / 0.80 | 2.7 km |
| 30 km | 0.866 (0.856–0.872) | 0.67 atm | −4 °C | 0.79 / 0.73 | 3.5 km |
| 35 km | 0.796 (0.787–0.802) | 0.60 atm | −8 °C | 0.73 / 0.67 | 4.3 km |
| 40 km | 0.733 (0.722–0.738) | 0.54 atm | −14 °C | 0.67 / 0.62 | 5.0 km |
| 45 km | 0.674 (0.663–0.679) | 0.49 atm | −19 °C | 0.61 / 0.57 | 5.8 km |
| 55 km | 0.570 (0.560–0.574) | 0.39 atm | −30 °C | 0.52 / 0.48 | 7.3 km |
| 70 km | 0.438 (0.429–0.442) | 0.28 atm | −48 °C | 0.40 / 0.37 | 9.5 km |

- **Lift.** A cubic metre of hydrogen lifts 1.28 kg at sea level and 0.61–0.73 kg in the flight band (35–45 km), where
  the air is Earth's at 4.3–5.8 km. A ship that must keep its lift to the top of the band carries 1.6 m³ of gas for
  every kilogram it lifts.
- **Steady air.** At any height the density varies by about 2% across the whole Moon, over the lunar day and over the
  years (1st to 99th percentile): a ship's lift barely changes along a route.

The wind, from the design run's 3-day means over the whole Moon and the cloud-resolving equatorial ring's 3-hourly
snapshots, which hold its storms:

| Height | 3-day means: median / 99th percentile / highest / 50-year | 3-hourly, with storms: median / 99th percentile / highest | Highest updraft |
|---|---|---|---|
| 1 km (ring), sea level (GCM) | 3.8 / 7.9 / 11.5 / 11.8 m/s | 2.4 / 8.4 / 14.1 m/s | 1.8 m/s |
| 10 km | 2.6 / 5.6 / 9.1 / 9.5 m/s | 1.7 / 6.1 / 12.7 m/s | 7.6 m/s |
| 20 km | 2.9 / 6.3 / 9.6 / 10.1 m/s | 2.6 / 6.4 / 13.4 m/s | 10.2 m/s |
| 30 km | 3.5 / 7.3 / 10.8 / 11.0 m/s | 3.8 / 7.2 / 14.2 m/s | 12.1 m/s |
| 40 km | 4.7 / 8.7 / 11.8 / 12.1 m/s | 5.0 / 8.9 / 20.6 m/s | 14.3 m/s |
| 55 km | 6.5 / 11.6 / 16.2 / 16.3 m/s | 6.7 / 10.0 / 16.4 m/s | 13.0 m/s |
| 70 km | 7.8 / 12.7 / 15.9 / 16.1 m/s | 7.9 / 12.4 / 17.9 m/s | 2.6 m/s |

- **Light winds.** In the flight band the wind's median is 4–5 m/s and its 99th percentile 8–10 m/s; the highest
  anywhere in the band in ten model years of 3-day means is 13.4 m/s, and the ring's highest 3-hourly wind at 40 km,
  in a storm, is 20.6 m/s. A ship cruising at 30 m/s makes headway against all of them.
- **Quiet air.** In the band the vertical wind stays under 0.6 m/s 99% of the time. Storm cloud reaches it in 0.4% of
  the ring's land columns and 0.2% of its sea columns, in the afternoon and evening, with updrafts up to 14.3 m/s
  inside.
- **Design winds for a moored ship.** Raised as the summit tower study raises its winds (1.4 for a 3-second gust,
  1.2 for decades), the band's highest modelled wind gives a design wind of 34.6 m/s and its 99th percentile a
  service wind of 13.6 m/s.

## Six times larger: the rule of similarity

Lunar gravity is 1/6.04 of Earth's, and three things scale by that factor together.

**The air.** The scale height is six times Earth's, so from 9 km up the Moon's air has, at each height, the density
of Earth's air at a sixth of the height above 9 km:

| Earth's height | 0 km | 1 km | 2 km | 3 km | 4 km | 5 km | 6 km | 8 km | 10 km |
|---|---|---|---|---|---|---|---|---|---|
| Air (kg/m³) | 1.225 | 1.112 | 1.006 | 0.909 | 0.819 | 0.736 | 0.660 | 0.525 | 0.413 |
| The Moon's height with that air | 9.0 km | 15.0 km | 21.0 km | 27.1 km | 33.3 km | 39.7 km | 46.3 km | 59.8 km | 73.4 km |
| 9 km + 6.04 × Earth's height | 9.0 km | 15.0 km | 21.1 km | 27.1 km | 33.2 km | 39.2 km | 45.2 km | 57.3 km | 69.4 km |

The temperature falls about 1.04 °C per km through the flight band, Earth's 6.5 °C per km stretched six times.

**The weather.** Convection at lunar gravity is Earth's stretched six times in size and time with the same speeds
(the premise of the [gravity pair](../../../climate/crm/README.md#the-gravity-pair)): the ring's storms rain over
30 km at the median and reach 28 km, with updrafts of 12–14 m/s.

**The flyer.** Scale a flyer up 6.04 times in every length and fly it at the same speed through air of the same
density. Its mass grows 220 times; its weight and its lift both grow 36 times. The stresses its weight puts in it stay
the same (they grow with g L), and so do the stresses from a gust of the same speed (they depend on ρ V U) and its
Froude number (V² / g L). It is the same flyer, six times longer. Its cover, cells and minimum gauges, which do not
scale with load, become a smaller share of it. [engineering/flight](../../../engineering/flight/README.md) holds this
rule in both of its models, and its tests check it to nine places.

So the Open Moon's air carries flyers six times the size of Earth's largest, with the same engineering. In air of
the same density, at the same speed:

| Earth's flyer | Its size | Its mass | On the Moon, six times larger |
|---|---|---|---|
| LZ 129 Hindenburg (flew 1936) | 245 m long | 207 t gross lift | 1.48 km, 45,500 t |
| CargoLifter CL160 (designed about 2000) | 260 m long | 160 t cargo | 1.57 km, 35,200 t |
| SkyCat 1000, a hybrid (designed) | 307 m long | 1,000 t payload | 1.85 km, 220,000 t |
| Antonov An-225 (flew 1988) | 88.4 m span | 640 t take-off | 530 m, 141,000 t |
| Stratolaunch Roc (flew 2019) | 117 m span | 590 t take-off | 710 m, 130,000 t |
| Lockheed CL-1201 (studied about 1970) | 340 m span | 5,375 t gross | 2.05 km, 1.18 million t |
| Mil Mi-26 (flew 1977) | 32 m rotor | 56 t take-off | 190 m, 12,300 t |

The flight band's air is 0.60 of Earth's at sea level, so a buoyant ship there lifts 0.59 of these masses.

Earth's own studies found that structure did not stop aircraft growing: NASA's span-loader comparison (Toll 1980)
found the payload share "almost no variation up to gross weights of at least 17.79 MN", and Kroo's study of very
large transports (1996) that "basic aerodynamics and structure do not limit the size of aircraft that can be operated
economically", with 600–800 passengers a reasonable bound for other reasons.

For one flyer of a given size, the design air against Earth's at sea level (the ecology register's figures, H3):

| | Moon ÷ Earth |
|---|---|
| Weight | 0.166 |
| Hover power | 0.063 |
| Energy per km at a lift-to-drag ratio | 0.166 |
| Speed at a lift coefficient | 0.38 |

## Buoyant ships by size

**The calibration.** Goodyear's 1975 study for NASA (CR-137692) lists the weights of four rigid airships by group, from
the German weight reports and the Macon's final weight statement. The model sizes each group from its loads and fits
one factor per group to all four:

| Group | What sizes it in the model | Factor |
|---|---|---|
| Longitudinals | Woodward's gust moment and the static moment with a cell deflated, at midship | 0.98 × ideal |
| Frames, wiring, gangways | the gas pressure the netting brings to the rings, and the hull's bending | 7.5 × ideal rings + 0.56 × ideal longitudinals |
| Empennage | the fins' gust load | 14.7 × ideal spar caps |
| Outer cover (hull and fins) | area | 0.25 kg/m² |
| Gas cells, nettings, valves | area, with the cells' bulkheads | 0.20 kg/m² |
| Power plant, fuel system, water recovery | power at top speed | 6.9 kg/kW |
| Fixed equipment | gross lift | 5.4% |

| Ship | Gross lift | Longitudinals: statement / model | Frames | Empty weight | Empty share |
|---|---|---|---|---|---|
| LZ 129 Hindenburg (hydrogen) | 206.9 t | 12.8 / 11.8 t | 37.3 / 37.3 t | 112.9 / 110.0 t | 55% / 53% |
| ZRS-5 Macon (helium) | 183.0 t | 10.4 / 11.0 t | 32.4 / 32.2 t | 107.3 / 100.0 t | 59% / 55% |
| LZ 127 Graf Zeppelin (hydrogen) | 118.1 t | 9.0 / 9.2 t | 18.7 / 17.2 t | 62.1 / 63.0 t | 53% / 55% |
| ZR-3 Los Angeles (helium) | 69.5 t | 5.1 / 5.0 t | 9.1 / 9.7 t | 40.5 / 39.2 t | 58% / 57% |

- **Weight sizes most of a rigid airship.** The frames follow the gas pressure: fitted to the bending alone they miss
  the four ships by up to 30%, to the gas pressure by up to 16%, to both by up to 8%. In the Hindenburg 17% of the lift
  went to structure sized by weight, 9% to structure sized by gusts, and 9% to cover and cells.
- **So lunar gravity cuts the largest part.** At 0.16 g the weight-driven structure of a ship of the same size is a
  sixth, and it only overtakes the gust structure at about 1 km.

**By size.** The share of the gross lift left for payload, fuel and crew, for hulls of fineness 6 with hydrogen at
36 m/s top speed ("none" where the structure outweighs the lift):

| Length | 100 m | 245 m | 500 m | 1 km | 1.5 km | 2 km | 3 km | 5 km | 10 km | 20 km |
|---|---|---|---|---|---|---|---|---|---|---|
| Earth, 1930s technology, sea level | 28% | 48% | 40% | 9% | none | none | none | none | none | none |
| Moon, 1930s technology, air of Earth's sea-level density (9 km) | 34% | 62% | 70% | 69% | 65% | 60% | 49% | 26% | none | none |
| Moon, 1930s technology, flight band | 16% | 53% | 64% | 64% | 61% | 56% | 45% | 22% | none | none |
| Moon, modern, flight band | 68% | 81% | 85% | 86% | 85% | 84% | 82% | 76% | 60% | 30% |
| Earth, modern, sea level | 73% | 81% | 80% | 72% | 63% | 54% | 36% | none | none | none |

| | Weight-driven structure overtakes gust structure | Best length (useful share) | Back to the Hindenburg's 45% | Structure takes all the lift |
|---|---|---|---|---|
| Earth, 1930s | 120 m | 270 m (48%) | 380 m | 1.2 km |
| Moon, 1930s, in air of Earth's sea-level density | 750 m | 650 m (70%) | 3.3 km | 7.2 km |
| Moon, 1930s, flight band | 1.0 km | 750 m (65%) | 3.0 km | 7.0 km |
| Moon, modern, flight band | 1.0 km | 870 m (86%) | 14.9 km | 30 km |
| Earth, modern, sea level | 120 m | 300 m (81%) | 2.5 km | 4.9 km |

- **The ships of the 1930s were built at their best size.** In the model the Hindenburg's technology does best near
  270 m on Earth; the Hindenburg was 245 m and the Akron and Macon 239 m.
- **The flight band.** The band's air lifts 0.61 kg per m³ of hydrogen (cells full at 45 km), 59% of Earth's at sea
  level, and its design gust (14.3 m/s) is a third above Earth's, so small ships do worse there than on Earth. From
  about 200 m up the Moon's weak gravity outweighs both.
- **Modern materials.** The modern case keeps the fitted factors and changes the materials: a carbon-fibre truss
  (300 MPa at 1.57 g/cm³ against the duralumin girders' 140 MPa at 2.79), the cover and film gas cells Goodyear
  proposed in 1975 (114 and 78 g/m² with nettings), and an electric drive with fuel cells at an assumed 2 kg per kW.

## Pressure hulls

A non-rigid ship holds its shape by pressure. Its envelope must stay taut under the gust's bending (FAA-P-8110-2,
§4.43) and hold the pressure at its crown with a fabric four times as strong as the limit load (§4.43(b)). The
largest hull of fineness 4 each fabric holds, with hydrogen, in the flight band's design gust (14.3 m/s) and Earth's
(10.7 m/s):

| Fabric | Top speed | Flight band | Moon at sea level | Earth at sea level |
|---|---|---|---|---|
| Today's strongest hull laminates (Vectran, Zylon), about 1,000 N/cm | 25 m/s | 350 m, 1.3 million m³, 870 t | 210 m, 260,000 m³, 340 t | 160 m, 125,000 m³, 140 t |
| | 30 m/s | 310 m, 0.92 million m³, 620 t | 180 m, 180,000 m³, 230 t | 150 m, 105,000 m³, 120 t |
| | 36 m/s | 280 m, 0.62 million m³, 410 t | 160 m, 110,000 m³, 140 t | 140 m, 84,000 m³, 90 t |
| A laminate twice as strong | 30 m/s | 540 m, 4.6 million m³, 3,100 t | 330 m, 1.0 million m³, 1,300 t | 240 m, 390,000 m³, 440 t |

- **The dynamic pressure sets it.** The envelope's tension is its pressure times its radius, and the pressure must
  exceed the dynamic pressure, which gravity leaves alone. The band's thinner air lowers the dynamic pressure at a
  given speed, and the gas's weight adds a sixth of Earth's head at the crown: non-rigid ships in the band reach about
  twice Earth's length.
- **Against Earth's ships.** The Zeppelin NT (75 m, 8,425 m³) and the Airlander 10 (98 m, 38,000 m³) sit well inside
  Earth's column. By this rule the CargoLifter CL160's envelope, 550,000 m³ and 65 m across, needs 2.2–2.5 times the
  strength of today's laminates at 90–125 km/h; she was designed with a keel carrying her loads on curtains inside
  the envelope, and never built.

## Winged ships

The wing model is fitted to the wing groups of ten transports from the DC-9 to the C-5A (Roskam, *Airplane Design*
part V): spar caps for the bending of an elliptic lift distribution at an ultimate load factor of 3.75, times 0.90,
plus 31.5 kg per m² of wing. It gives the ten within 12%, except the A300 B2 (24% light).

- **Wings grow eight times.** At the An-225's wing loading (6,935 Pa) and a manoeuvre load factor of 2.5, a wing
  carries its weight with the An-225's share (16% in the model) at a span of about 710 m on the Moon, eight times
  the An-225's 88.4 m, both in the flight band and in air of Earth's cruise density. The ship then has a mass of about
  250,000 t and cruises at 195 m/s in the band. The factor is more than six because the areal part of the wing is a
  sixth of the share at the same loading.
- **Lighter and slower.** At a quarter of that wing loading, a ship of the same span has a quarter of the mass
  (62,000 t), cruises at half the speed (97 m/s) and gives 18% of its mass to its wing.
- **Gusts.** A gust of 20 m/s (Earth's 50 ft/s design gust at the band's density, and the ring's 14.3 m/s updraft
  with the gust factor) adds less than the manoeuvre factor at this loading.
- **Hybrids.** A hull's own lift at a few degrees' pitch carries, at 30 m/s in the band, 3.7 times its buoyant lift
  at 40 m across, 1.8 times at 83 m and 0.9 times at 167 m; at Earth's sea level, 0.56, 0.27 and 0.13 times. A hybrid
  ship on the Moon can carry as much on its hull's lift as on its gas.
- **What wings need.** Winged ships need runways, water or lift fans to take off and land, and cannot hold still at a
  berth; the summit port's docks serve buoyant ships and hybrids that can hover.

## The atmosphere's own scale

The storms are the air's largest structures a ship meets. On the ring, over its second lunar day:

| | Median | Tenth largest | Largest |
|---|---|---|---|
| Width of a storm's rain | 30 km | 96 km | 246 km |
| Top of its cloud | 28 km | 68 km | 82 km |
| Its strongest updraft | 1.3 m/s | 5.2 m/s | 12.4 m/s |
| Its life | 6 hours | 12 hours | 36 hours |

- **Storms drift slowly**, between 2.7 m/s westward and 3.1 m/s eastward, and about three rain somewhere on the
  10,900 km ring at any time. The largest feed anvils about 1,000 km wide between 20 and 75 km.
- **Ships up to about 3 km are small against them.** A rigid airship turns on a radius of about 4.6 of its lengths
  (the Los Angeles, NACA Report 333), and the radius does not depend on speed. A 1.5 km ship turns on about 7 km and
  reverses course in 12 minutes at 30 m/s; a median storm is 30 km wide and moves at a tenth of the ship's speed. At
  3 km a ship's turning circle, 28 km across, is as wide as a median storm.
- **The worst gust is measured in ship lengths.** In Calligeros and McDavitt's theory a ship's bending peaks when a
  gust builds up over half its length, and changes little between a quarter and three-quarters of it; Woodward's
  envelope takes that worst case. The ring resolves updrafts only over its 6-km columns, so the sharper gusts inside
  them are unknown, and every ship here is sized for the worst case. The design gust for ships in the band is
  14.3 m/s, the ring's strongest updraft at 40 km, a third above the Guggenheim Airship Institute's 35 ft/s
  (10.7 m/s).

## What grows with a ship

Rigid hydrogen ships in the flight band, with modern materials, cruising at 30 m/s on trips of up to 48 hours (half
the Moon's circumference is 50 hours at that speed), with fuel cells and their hydrogen fuel and tanks; 250–600 kg
carried per passenger with a cabin, services and a share of the crew:

| Length | Across | Gas | Gross lift | Useful | Passengers (400 kg; 600–250 kg) | Hydrogen aboard | Cruise power | Energy per passenger-km | Moored in the design wind: nose-on / side-on | Turning radius, half-turn | Long-haul berths at the port |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 245 m | 41 m | 194,000 m³ | 119 t | 81% | 210 (140–340) | 9 t | 1.2 MW | 52 Wh | 0.04 / 2.1 MN | 1.1 km, 2 min | 150–250 |
| 500 m | 83 m | 1.65 million m³ | 1,010 t | 85% | 2,020 (1,350–3,200) | 76 t | 4.9 MW | 23 Wh | 0.15 / 8.8 MN | 2.3 km, 4 min | 16–26 |
| 750 m | 125 m | 5.6 million m³ | 3,420 t | 86% | 7,010 (4,690–11,100) | 257 t | 11 MW | 15 Wh | 0.34 / 20 MN | 3.4 km, 6 min | 5–8 |
| 1 km | 167 m | 13 million m³ | 8,100 t | 86% | 16,800 (11,200–26,600) | 609 t | 20 MW | 11 Wh | 0.61 / 35 MN | 4.6 km, 8 min | 2–3 |
| 1.5 km | 250 m | 44.5 million m³ | 27,300 t | 85% | 56,600 (37,900–89,900) | 2,060 t | 44 MW | 7 Wh | 1.4 / 79 MN | 6.9 km, 12 min | 1 |
| 2 km | 333 m | 106 million m³ | 64,800 t | 84% | 133,000 | 4,870 t | 79 MW | 5.5 Wh | 2.4 / 140 MN | 9.2 km, 16 min | under 1 |
| 3 km | 500 m | 356 million m³ | 219,000 t | 82% | 436,000 | 16,400 t | 177 MW | 3.8 Wh | 5.5 / 320 MN | 13.8 km, 24 min | under 1 |

- **People.** With the Hindenburg's own materials the numbers are a quarter to 30% lower from 750 m up (5,200 at
  750 m, 12,400 at 1 km, 40,000 at 1.5 km).
- **Hydrogen.** A 500 m ship holds 76 t of hydrogen, 4.6 times the Hindenburg's; a 1.5 km ship 2,060 t, 123 times.
  Helium would need about twice the mass and lift about 7% less.
- **Mooring.** Held at the nose, a ship turns into the wind, taking 0.15–1.4 MN from 500 m to 1.5 km in the band's
  design wind of 34.6 m/s; caught side-on before it turns, 9–79 MN. Those are the loads the crown's docking arms
  answer for.
- **Berths.** The port's long-haul traffic, 32,000–53,000 passengers an hour, with each call unloading and loading a
  full ship in two hours.
- **Energy.** Drag alone at 30 m/s; hotel loads of 500 W per passenger are in the fuel.

## What this leaves out, and what comes next

- **The distribution of flyers.** The next step: the classes from personal wings, canopies and micro gliders through
  air taxis, regional and long-haul ships and freighters to aerial platforms; for each its size, speed, energy,
  height band and capacity from these models; how many of each the Moon's trips and freight need; and how they share
  the air by height. It starts from the author's decisions for the summit metropolis
  ([research/decisions.md](../../decisions.md)): about 100 million people within about 32 km of the tower, small sky
  boats serving much as cars do, gliders throughout with glider and parkour zones, and large sky ferries on a
  backbone with metro and rail.
- **The long-haul class and the crown's port.** The chosen class sets the berths, the docking arms and their loads in
  the flight band.
- **Hydrogen.** The hydrogen aboard grows with a ship's volume. A hydrogen ship needs cells kept apart from the air and
  from ignition, and a helium ship needs a supply of helium. Fire at 0.16 g is worked out for buildings in the
  [port fire study](../port_fire/README.md); ships are the next case.
- **The lunar day.** Sunlight warms a ship's gas for 354 hours and the night cools it for as long; the lift that
  swings with it, and the ballast and water to balance it, are to be worked out.
- **Materials and non-structural parts.** The fit covers four ships of duralumin girders and cotton fabrics, built
  1924–1936. The modern case carries its factors to a carbon-fibre truss, laminates and electric drive, an assumption
  to check against a modern rigid design such as LTA Research's Pathfinder 1 once its weights are published.
- **Winged and hybrid ships.** Takeoff and landing (runways, water or lift fans), sweep and flutter, and the hull's own
  lift in hybrids, which at 0.16 g carries several times more of a ship than on Earth.
- **Building and keeping them.** Sheds or open-air yards for ships a kilometre long, and the handling of ships at
  berths in the design wind.
