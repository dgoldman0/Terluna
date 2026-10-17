# Building at a sixth of Earth's gravity

**17 October 2026.** A design note for the Open Moon's architecture: guards and railings, stairs, handholds, ramps
and halls for moving, at the Moon's gravity of 1.624 m/s², 0.166 of Earth's. It carries the author's rule for
low-gravity architecture, set in their ChatGPT work of 15 October and confirmed on 17 October
([register](../../research/decisions.md)): the architecture's guards and circulation show the Moon's gravity and the
added spring in a step, and do not look terrestrial. It builds on the
[daily-life study](../../research/studies/daily_life/README.md#5-homes-and-mobility-across-bodies-and-generations)'s
falls, sixfold jump, braking distances and access design, and on the 2.2 m railings of the summit port's
[ring 0 drawings](../../immersion/world/summit/ring0/draw.py).

[movement.py](movement.py) derives every figure from the shared gravity and writes
[results/low_gravity.json](results/low_gravity.json) (schema `terluna.habitation.low-gravity/1`) at once. The Earth
codes and measurements it reads are in [sources.json](sources.json). Each figure below is *derived*, closed-form
mechanics on the shared gravity and the cited values, or a *design choice*, marked as such.

## What a sixth of the weight changes

- **The same push goes six times as high.** A hop leaves the ground at the speed it would on Earth and rises 6.04
  times as high, staying aloft 6.04 times as long: a hop that rises 10 cm on Earth rises 60 cm and lasts 1.7 s. The
  countermovement jump of fit young men, 0.35 m on Earth (Markovic et al. 2004), rises 2.13 m and hangs 3.2 s; with
  an arm swing, 0.52 m becomes 3.15 m (*derived*).
- **Steady running stays low.** Runners in simulated low gravity take off more slowly, their vertical speed scaling
  as the square root of gravity, and bounce less high than on Earth (Polet et al. 2018), so each flight lasts about
  2.5 times as long (*derived* from their scaling). At lunar gravity skipping, the gait the Apollo crews chose, costs
  barely more than running, and both cost about what walking does on Earth (Pavei et al. 2015). The spring in a step
  is the hop people choose, which a light body invites at every edge and stair.
- **Falls land like a sixth of the drop on Earth.** A fall of 3 m lands like 0.5 m on Earth, 10 m like 1.7 m and 30 m
  like 5 m, and lasts six times as long as an Earth fall landing at the same speed (*derived*, without drag; long falls
  end at 17.9 m/s, like 16 m on Earth, in the daily-life study).
- **Grip on the floor stops and turns a sixth as hard.** From the lunar walk–run transition, 1.42 m/s (De Witt et al.
  2014), a person stops in 1.0 m against Earth's 0.17 m, and the tightest turn has a radius of 2.1 m against 0.34 m;
  from 3 m/s, 4.6 m and 9.2 m (*derived*, with the daily-life study's friction of 0.6).
- **A hand holds six times its grip.** The weak grip of 16 kg that a quarter of people have by 80 (Dodds et al. 2014)
  holds a 97 kg body on the Moon, and two hands of the median woman's 31 kg stop a 75 kg adult's fall of 1.2 m over
  0.3 m of give, where on Earth they cannot hold that adult's weight (*derived*, a hand's hold taken as its grip
  strength).
- **Impacts keep Earth's energy.** A 75 kg body meeting a rail at 2 m/s brings 150 J on either world, while a standing
  person's push, limited by the grip of their feet, falls from 441 N to 73 N (*derived*).

## Guards and railings

**The rule** (*derived*). A guard of Earth's height keeps Earth's margin over the centre of mass of a person who
stands, leans or stumbles at it. To keep that margin at the top of a hop that rises j on Earth, the lunar guard rises
by the extra height the same hop gains: H = H_Earth + (g_E/g_M − 1) j = 1.1 m + 5.04 j. Steady running needs no more
than Earth's guard, since it bounces less high than on Earth.

| Hop, as it rises on Earth | On the Moon | Guard keeping Earth's margin, from the UK's 1.1 m; from the IBC's 1.067 m |
|---|---:|---:|
| None: standing, leaning, stumbling, running | 0 | 1.10; 1.07 m |
| 5 cm | 0.30 m | 1.35; 1.32 m |
| 10 cm | 0.60 m | 1.60; 1.57 m |
| 20 cm | 1.21 m | 2.11; 2.07 m |
| 22 cm, the hop the summit drawings' guard covers | 1.32 m | 2.20; 2.17 m |
| 35 cm, a fit young man's full countermovement jump | 2.13 m | 2.88; 2.85 m |
| 52 cm, the same jump with an arm swing | 3.15 m | 3.73; 3.70 m |

The summit drawings' 2.2 m covers a hop that rises 0.22 m on Earth, about three-fifths of a fit young man's full
jump. A full jump clears it, as a vault clears Earth's 1.1 m. In the hop it covers, the hands rise 1.32 m, which puts
a 2.2 m top within almost anyone's reach, so the guard's form completes what its height begins: its top gives no grip.

**Where guards go** (*design choice*). Earth's codes leave unguarded the drops of 380 mm (UK, buildings other than
dwellings), 600 mm (UK dwellings) and 762 mm (IBC 1015.2). A lunar fall lands at those speeds from 2.3, 3.6 and 4.6 m.

| Drop below the edge | Lands like a fall on Earth of | Guard |
|---|---|---|
| Up to 2.3 m | Up to 0.38 m | None; a kerb, a change of paving, planting or water marks the edge |
| 2.3–4.6 m | 0.38–0.76 m, which Earth's codes guard in some places and not others | Earth's own 1.1 m |
| Above 4.6 m | Above 0.76 m | The lunar guard, 2.2 m |
| Decks aloft and the port's rims | Towards the 17.9 m/s end of a long fall | The lunar guard, with canopies at the rims as the daily-life study sets out |

Along stairs the same classes apply, measured from the line of the nosings, since people bound on stairs.

**The lunar guard** (*design choice*):
- continuous, with no opening a 100 mm sphere passes, as on Earth, since bodies are Earth's size;
- with no rail, ledge or pattern for a foot or a hand, since at a sixth of the weight every foothold is a ladder
  for a child;
- with a handrail along its inner face at Earth's 0.9–1.0 m, where the arm reaches;
- with a top that turns inward and down over the walk, smooth and without an edge or a shelf to press on;
- clear or finely open at eye level, so that the view and the light of the whole sky pass;
- built to Earth's code loads in full, since impacts keep Earth's energy.

## Stairs

**What the gravity does** (*derived*):
- Climbing a 150 mm riser takes 0.24 J/kg, what a 25 mm riser takes on Earth; Earth's steepest general-access riser,
  170 mm, takes 1.67 J/kg there.
- Stepping down 150 mm falls for 0.43 s and lands at 0.70 m/s, like stepping off 25 mm on Earth. The riser that lands
  like Earth's 170 mm is 1.03 m, so effort and landing leave the riser to the leg's reach, as on Earth.
- People bound: clearing four 150 mm risers in one bound takes a take-off of 1.4 m/s, the push of a 10 cm hop on
  Earth.
- A tumble down a flight ends no faster than a fall of the flight's rise, so the same end speed allows 6.04 times
  Earth's rise: the UK's general-access flight of 12 risers of 170 mm (2.04 m) becomes 12.3 m, and the IBC's 3.66 m
  becomes 22.1 m.
- A stop from 2 m/s takes 2.05 m of landing.

**Two stairs** (*design choice*):

| | Everyday stair | Bounding stair | Earth's general-access stair (UK) |
|---|---:|---:|---:|
| Rise | 150 mm | 300 mm | 150–170 mm |
| Going | 400 mm | 600 mm | 250–400 mm |
| Pitch | 20.6° | 26.6° | 21–34° |
| 2R + G | 700 mm | 1,200 mm | 550–700 mm |
| Work a step | 0.24 J/kg | 0.49 J/kg | 1.47–1.67 J/kg |
| A step down lands like Earth's | 25 mm | 50 mm | 150–170 mm |

- **The everyday stair** keeps each step within Earth's general-access range for children, older people and anyone
  carrying things, and its deep treads take a bound's landing on the whole foot.
- **The bounding stair** climbs a hop a step, each with the push of a 5 cm hop on Earth (take-off 0.99 m/s). It runs
  beside the everyday stair in parks, plazas and halls; with seat-height steps it becomes a terrace to sit on.
- **Flights** rise up to 12 m between landings, about the UK flight's lunar equivalent.
- **Landings** are at least 3 m deep in the direction of travel, and 3 m of level floor leads to the head of each
  descending flight.
- **Headroom** over stairs and through doors is 2.5 m: Earth's 2.0–2.03 m and the 0.5 m that a hop rising 10 cm on
  Earth gains at lunar gravity.
- **Handrails** run on both sides of every flight and through its landings, at Earth's 0.9–1.0 m and 32–51 mm across
  (Approved Document K; ADA 505).

## Handholds

On the Moon a handrail catches a fall as well as steadying a step (*derived*, for a 75 kg adult and 0.3 m of give):

| Grip (Dodds et al. 2014) | One hand holds, on the Moon | Fall stopped by one hand | By two hands | By two hands on Earth |
|---|---:|---:|---:|---|
| Women's weak grip, 16 kg | 97 kg | 0.09 m | 0.47 m | Cannot hold the weight |
| Men's weak grip, 27 kg | 163 kg | 0.35 m | 1.00 m | Cannot hold the weight |
| Women's peak median, 31 kg | 187 kg | 0.45 m | 1.20 m | Cannot hold the weight |
| Men's peak median, 51 kg | 308 kg | 0.93 m | 2.16 m | 0.11 m |

Rails run along every stair, ramp, guard and lane, with grab bars at doors, edges and landings and, in the halls for
moving, grips at several heights to catch, turn and launch from (*design choice*). Every handhold is a climbing hold
for a child, so the guards' faces carry none.

## Ramps and slopes

- **Climbing costs a sixth** (*derived*). A 1:12 ramp climbs with the force a 1:73 ramp takes on Earth, and a 30°
  slope climbs like Earth's 1:12.
- **Descending keeps Earth's proportions** (*derived*). Grip falls with gravity as the slope's pull does: at 1:12 the
  slope takes 14% of the grip on either world, and a wheel that rolls a whole run unbraked needs the same level run-out
  at its foot on both, 0.138 m for each metre of a 1:12 run. A landing of length L stops such a runaway from any run
  rising up to 0.6 L, whatever its slope. Stopping on the ramp itself from walking speed takes 1.2 m, against 0.2 m on
  Earth.
- **Grip sets the limit on a steep slope** (*derived*). At 18° the slope takes 55% of the grip and a stop from walking
  speed takes 2.4 m; at 25°, 78% and 5.1 m; at 30°, the slope that climbs like Earth's 1:12, 96% and 32 m.
- **Ramps** (*design choice*) keep Earth's gradients, 1:12 for short runs and 1:20 for long ones (Approved Document K,
  diagram 2.1), and rise up to 1.8 m a run between 3 m landings. Earth's ADA runs rise 0.76 m between 1.525 m landings,
  which stop a runaway from 0.91 m. At walking speed a ramp turns on a radius of 2.1 m, so ramps run straight or turn on
  wide landings, and steep paths in the parks take steps or high-grip surfaces.

## Halls for moving

The summit metropolis designs in glider and parkour zones, so that the rest of its architecture does not invite
either ([register](../../research/decisions.md#the-summit-port-and-metropolis)), and the scene the author accepted on
15 October had a separate hall for moving. A hall's dimensions follow from the jump (*derived*), and the sizes chosen
are a *design choice*:
- **Height.** The arm-swing jump rises 3.15 m: a 1.9 m adult's head reaches 5.05 m and raised arms about a metre more,
  so the halls stand 6.5 m clear for jumps from the floor, and higher below platforms.
- **Length.** A running hop at 3 m/s with a full jump's take-off flies 9.7 m in 3.2 s, against 1.6 m on Earth.
- **Turns.** Grip alone turns a runner at 3 m/s on a radius of 9.2 m. Banked floors tighten the turn: 17° of bank for a
  5 m radius at 3 m/s, 31° for 3 m, and 20° for walking speed round 1 m.
- **Landing.** A jump lands at its take-off speed, as on Earth, so sports floors serve.

## The visual criteria

The author's criteria: the guards and circulation show the Moon's gravity and the added spring in a step, and do not
look terrestrial. The scene accepted on 15 October carried continuous tall guards, inhabited spans with visible
support, human-scale spaces and a separate hall for moving, and the author found it suited the semi-arid parks,
gardens and their amenities, its setting left provisional (register). The figures above give these forms (*design
choice*):
- **Tall, continuous guards** of 2.2 m along every long drop, clear or finely open, their tops curving inward. Earth's
  waist-high rail appears only along drops of 2.3–4.6 m.
- **Open edges.** Level changes under 2.3 m carry no guard: terraces, stepped gardens, plinths and quays that people
  step and hop down.
- **Long, low stairs.** Deep treads in flights rising up to 12 m between deep landings, with bounding stairs of
  knee-high steps beside them.
- **High openings.** Door heads and stair headroom of 2.5 m.
- **Rounded and banked circulation.** Corners on walking routes rounded to 2 m or more, and banked lanes where people
  run.
- **Handholds everywhere**, shaped to be seen and used.
- **Stepping stones for a hop.** A hop at 1.4 m/s along and 1.0 m/s up, the push of a 5 cm hop on Earth, flies 1.7 m
  (0.3 m on Earth), so stones 1.5–1.8 m apart across the parks' rills and mirror pools make a lunar crossing.
- **Halls for moving**, tall and long, apart from calm rooms at human scale.

## The summit drawings

The ring 0 drawings carry the lunar guard already: 2.2 m railings along both edges of the roof park, over the port's
hollow and over the lanes' voids. Two elements of their section fall in its class and stand lower: the inner gallery's
balcony over the hollow at 1.2 m, and the galleries' balustrades along the 9 m nave at 1.1 m
([draw.py](../../immersion/world/summit/ring0/draw.py), lines 413 and 416). Raising them to 2.2 m belongs to the
drawings' next revision.

## Evidence and boundaries

The figures are closed-form mechanics at the shared gravity with the cited codes and measurements: falls without drag
(below 30 m drag lowers a lunar impact speed by under 7%), hops with Earth's take-off speeds, steady running after
Polet et al. (2018), floor grip of 0.6, and a hand's hold taken as its grip strength over 0.3 m of give. The
measurements of people in low gravity come from short simulations (harness, suspension and parabolic flight), and the
longest human stay on the Moon is 75 hours. Earth's codes carry Earth's experience of falls; a lunar fall lasts six
times as long, which may let a tumbling body turn further before it lands, and testing that is open. The figures use
adults' jumps and grips; children's are the next to add. The heights covered, the thresholds, and the stairs,
landings, ramps and halls are design choices for human-factors tests in low gravity. Where gravity is Earth's, in the
array's wheels, Earth's codes apply.

## Reproduction and checks

```sh
python -m habitation.low_gravity.movement                  # results/low_gravity.json, at once
OPENBLAS_NUM_THREADS=1 python -m pytest -q habitation/tests/test_low_gravity.py
```

The tests check that the product regenerates and binds its code and the shared constants; that hops rise, and falls
and stops scale, with the gravity ratio; that the short falls and the braking match the daily-life study's product;
that the guard rule starts at Earth's guard and gives the summit drawings' 2.2 m; the guard classes against the codes'
drops; the stairs' work, landing speeds, bounds, flight limits and landings; the handholds; that a ramp's runaway
run-out does not depend on gravity; and the banking and the halls. On 17 October: 13 passed.
