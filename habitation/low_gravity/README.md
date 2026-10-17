# Building at a sixth of Earth's gravity

**17 October 2026.** A design note for the Open Moon's architecture: guards and railings, stairs, handholds, ramps
and halls for moving, at the Moon's gravity of 1.624 m/s², 0.166 of Earth's. It carries the author's rule for
low-gravity architecture, set in their ChatGPT work of 15 October and confirmed on 17 October
([register](../../research/decisions.md#the-chatgpt-work-of-15-october-confirmed-17-october-2026)):
the architecture's guards and circulation show the Moon's gravity and the added spring in a step, in forms distinct
from Earth's. It builds on the
[daily-life study](../../research/studies/daily_life/README.md#5-homes-and-mobility-across-bodies-and-generations)'s
falls, sixfold jump, braking distances and access design, and on the 2.2 m railings of the summit port's
[ring 0 drawings](../../immersion/world/summit/ring0/draw.py).

[movement.py](movement.py) derives every figure from the shared gravity and writes
[results/low_gravity.json](results/low_gravity.json) (schema `terluna.habitation.low-gravity/1`) at once. It holds
the Earth codes and measurements it uses, each beside its source in [sources.json](sources.json), and a test checks
them against that file. Each figure below is *derived*, closed-form mechanics on the shared gravity and the cited
values, or a *design choice*, marked as such. The design choices await the author's review.

## What a sixth of the weight changes

- **The same take-off speed goes six times as high.** A hop that leaves the ground at its Earth speed rises 6.04
  times as high after take-off and stays aloft 6.04 times as long: a hop that rises 10 cm on Earth rises 60 cm and
  lasts 1.7 s. The countermovement jump of fit young men rises 0.35 m after take-off on Earth (Markovic et al. 2004),
  2.13 m here, and hangs 3.2 s. With an arm swing they reached 0.52 m from standing by a waist belt, 0.11 m of it the
  push-off's extension, which gravity leaves alone, so the jump reaches 2.60 m from standing here (*derived*).
- **The same effort leaves faster.** The legs lift a sixth of the weight, so the same leg work over a push-off of
  0.35–0.45 m (assumed) gives the countermovement jump 3.6–3.8 m/s in place of 2.6 m/s and a lunar rise of 3.9–4.4 m.
  That is an upper bound, since a muscle's force falls with its speed. Taking Earth's take-off speeds, as the figures
  here do, sets every jump low (*derived*).
- **Steady running stays low.** Runners in simulated low gravity take off more slowly, their vertical speed scaling
  as the square root of gravity, and bounce less high than on Earth (Polet et al. 2018), so each flight lasts about
  2.5 times as long (*derived* from their scaling). At lunar gravity skipping, the gait the Apollo crews chose, costs
  barely more than running, and both cost about what walking does on Earth (Pavei et al. 2015). The spring in a step
  is the hop people choose, which a light body invites at every edge and stair.
- **Falls land like a sixth of the drop on Earth.** With the daily-life study's drag, which ends a long fall at
  17.9 m/s, like 16 m on Earth, a fall of 3 m lands like 0.49 m on Earth, 10 m like 1.6 m and 30 m like 4.3 m, and
  lasts six times as long as an Earth fall landing at the same speed (*derived*).
- **Grip on the floor stops and turns a sixth as hard.** From the lunar walk–run transition, 1.42 m/s (De Witt et al.
  2014), a person stops in 1.0 m against Earth's 0.17 m, and the tightest turn has a radius of 2.1 m against 0.34 m;
  from 3 m/s, 4.6 m and 9.2 m (*derived*, at the daily-life study's friction of 0.6). The design figures below use
  half that grip and keep the other half in reserve (*design choice*): a stop from 1.42 m/s then takes 2.1 m and a
  level turn a radius of 4.1 m.
- **A hand holds six times its grip.** The weak grip of 16 kg that a quarter of people have by 80 (Dodds et al. 2014)
  holds a 97 kg body on the Moon, and two hands of the median woman's 31 kg stop a 75 kg adult's fall of 1.2 m over
  0.3 m of give, where on Earth they cannot hold that adult's weight (*derived*, a hand's hold taken as its grip
  strength).
- **Impacts keep Earth's energy, and lifting costs a sixth.** A 75 kg body meeting a guard at 2 m/s brings 150 J on
  either world, while carrying its centre of mass over the guard's top costs a sixth of Earth's for the same lift. A
  standing person's push, limited by the grip of their feet, falls from 441 N to 73 N (*derived*).

## Guards and railings

**Three rules** (*derived*):
- **Standing and leaning.** Tipping over a guard turns on where the centre of mass stands against the guard's top,
  which gravity leaves alone, so Earth's 1.1 m keeps Earth's margin for a person standing or leaning at it.
- **Hops.** A guard of H = H_Earth + (g_E/g_M − 1) j = 1.1 m + 5.04 j keeps the margin Earth's guard has at the top of
  the same hop, one that rises j after take-off on Earth.
- **Stumbling and running.** Carrying the centre of mass from its standing height c over the top costs g (H − c) per
  kilogram, while the body arrives with Earth's kinetic energy. The guard H = c + (g_E/g_M)(1.1 m − c) costs the same
  lift as Earth's 1.1 m, so that the same impact speed carries a body over it on both worlds if the same share of the
  impact's energy goes into the lift.

| Hop, rising after take-off on Earth | On the Moon | Guard keeping Earth's margin at its top, from the UK's 1.1 m; from the IBC's 1.067 m |
|---|---:|---:|
| None: standing and leaning | 0 | 1.10; 1.07 m |
| 5 cm | 0.30 m | 1.35; 1.32 m |
| 10 cm | 0.60 m | 1.60; 1.57 m |
| 20 cm | 1.21 m | 2.11; 2.07 m |
| 21.8 cm, the hop the summit drawings' 2.2 m keeps Earth's margin at | 1.32 m | 2.20; 2.17 m |
| 35.3 cm, a fit young man's countermovement jump | 2.13 m | 2.88; 2.85 m |
| 41.3 cm, the same jump with an arm swing | 2.49 m | 3.18; 3.15 m |

The energy rule, with the centre of mass at 0.55 of stature (assumed):

| Stature | Centre of mass | Lift over Earth's 1.1 m | Impact speed whose energy pays it, Earth | Guard costing that lift on the Moon | Impact speed paying the lift over 2.2 m, Moon |
|---|---:|---:|---:|---:|---:|
| 1.91 m | 1.05 m | 0.49 J/kg | 0.99 m/s | 1.35 m | 1.93 m/s |
| 1.75 m | 0.96 m | 1.35 J/kg | 1.64 m/s | 1.79 m | 2.00 m/s |
| 1.62 m | 0.89 m | 2.05 J/kg | 2.02 m/s | 2.15 m | 2.06 m/s |
| 1.40 m | 0.77 m | 3.24 J/kg | 2.54 m/s | 2.76 m | 2.16 m/s |
| 1.20 m | 0.66 m | 4.31 J/kg | 2.94 m/s | 3.32 m | 2.24 m/s |
| 1.00 m | 0.55 m | 5.39 J/kg | 3.28 m/s | 3.87 m | 2.32 m/s |

Both worlds' real thresholds lie above these speeds, since an impact spends most of its energy elsewhere; the rule
keeps them in Earth's proportion. The summit drawings' 2.2 m keeps Earth's margin at the top of a hop that rises
0.218 m on Earth (take-off 2.07 m/s), three-fifths of a fit young man's countermovement jump: with the hop rising from
standing height, an adult whose centre of mass stands at 0.96 m has it 8 cm above the top there, on both worlds. By
the energy rule it costs at least Earth's lift for every adult whose centre of mass stands at 0.88 m or higher, a
stature of 1.60 m or more; for children it costs less than Earth's 1.1 m does, and keeping Earth's lift for children
of 1.4–1.0 m takes 2.8–3.9 m. A full jump clears it, as a vault clears Earth's 1.1 m. In the hop it is sized for the
hands rise 1.32 m with the body, which brings its top within almost anyone's reach, so its top gives no grip.

**Guard bands** (*design choice, for the author's review*). Earth's codes leave unguarded the drops of 380 mm (UK,
buildings other than dwellings), 600 mm (UK dwellings) and 762 mm (IBC 1015.2). A lunar fall lands at those speeds
from 2.3, 3.7 and 4.7 m.

| Drop below the edge | Lands like a fall on Earth of | Guard |
|---|---|---|
| Up to 2.3 m | Up to 0.38 m | None; a kerb, a change of paving, planting or water marks the edge |
| 2.3–4.7 m | 0.38–0.76 m, which Earth's codes guard in some places and not others | One of two options: Earth's 1.1 m, which keeps Earth's margin for standing and leaning and costs stumbling and running a sixth of Earth's lift, or 2.2 m on every drop over 2.3 m |
| Above 4.7 m | Above 0.76 m | The lunar guard, 2.2 m |
| Above 4.7 m, where young children play | Above 0.76 m | Guards of 2.8–3.9 m, or play kept back from the edge |
| Decks aloft and the port's rims | Towards the 17.9 m/s end of a long fall | The lunar guard, with canopies at the rims as the daily-life study sets out |

Along stairs the same bands apply, measured from the line of the nosings, since people bound on stairs.

**Footholds** (*derived*; their use is a design choice). A surface a foot can use, at height s, raises a hop's paths:
those of a hop launched at speed v from a height y₀ reach no higher than y₀ + v²/2g − gx²/2v² at a distance x, so the
foothold gives a hop no higher crossing of the guard's plane than a hop from the floor at the guard once it stands at
least d = v√(2s/g) from the guard. For the hop the 2.2 m guard is sized for, that is 0.89 m for a 150 mm kerb,
1.54 m for a 450 mm bench seat or planter's rim and 2.30 m for a free-standing handrail at 1.0 m. Closer than that,
the guard's height counts from the foothold.

**The lunar guard** (*design choice*):
- continuous, with no opening a 100 mm sphere passes, as on Earth, since bodies are Earth's size;
- with no rail, ledge or pattern on its face for a foot or a hand, since at a sixth of the weight every foothold is a
  ladder for a child, and its height counted from any foothold within the reach above;
- its handhold a grip slot recessed into its face at Earth's 0.9–1.0 m, whose lower lip slopes away and down, so that
  it takes fingers and gives a toe nothing to stand on; a free-standing handrail stands at least 2.3 m back;
- with a top that turns inward and down over the walk, smooth and without an edge or a shelf to press on;
- clear or finely open at eye level, so that the view and the light of the whole sky pass;
- built to Earth's code loads in full, since impacts keep Earth's energy.

## Stairs

**What the gravity does** (*derived*):
- Climbing a 150 mm riser takes 0.24 J/kg, what a 25 mm riser takes on Earth; Earth's steepest general-access riser,
  170 mm, takes 1.67 J/kg there.
- Stepping down 150 mm falls for 0.43 s and lands at 0.70 m/s, like stepping off 25 mm on Earth. The riser that lands
  like Earth's 170 mm is 1.03 m, so effort and landing leave the riser to the leg's reach, as on Earth.
- People bound: clearing four 150 mm risers in one bound takes a take-off of 1.4 m/s, the speed of a 10 cm hop on
  Earth.
- A tumble down a flight ends no faster than a fall of the flight's rise, so the same end speed allows 6.04 times
  Earth's rise: the UK's general-access flight of 12 risers of 170 mm (2.04 m) becomes 12.3 m, and the IBC's 3.66 m
  becomes 22.1 m.

**Two stairs** (*design choice*):

| | Everyday stair | Bounding stair | Earth's general-access stair (UK) |
|---|---:|---:|---:|
| Rise | 150 mm | 300 mm | 150–170 mm |
| Going | 400 mm | 600 mm | 250–400 mm |
| Pitch | 20.6° | 26.6° | 21–34° |
| 2R + G | 700 mm | 1,200 mm | 550–700 mm |
| Work a step | 0.24 J/kg | 0.49 J/kg | 1.47–1.67 J/kg |
| A step down lands like Earth's | 25 mm | 50 mm | 150–170 mm |
| Landings, deep | 3.0 m | 4.5 m | At least the flight's width |

- **The everyday stair** keeps each step within Earth's general-access range for children, older people and anyone
  carrying things, and its deep treads take a bound's landing on the whole foot.
- **The bounding stair** climbs a hop a step, each with the take-off of a 5 cm hop on Earth (0.99 m/s). It runs beside
  the everyday stair in parks, plazas and halls; with seat-height steps it becomes a terrace to sit on.
- **Flights** rise up to 12 m between landings, about the UK flight's lunar equivalent.
- **Landings** stop a person at half the grip: the everyday stair's 3 m take a stop from the walk–run transition, 2.1 m
  and a foot, and the bounding stair's 4.5 m a stop from 2 m/s, 4.1 m and a foot. 3 m of level floor leads to the head
  of each descending flight.
- **Headroom** over stairs and through doors is 2.5 m: Earth's 2.0–2.03 m and the 0.5 m that a hop rising 10 cm on
  Earth gains at lunar gravity.
- **Handrails** run on both sides of every flight and through its landings, at Earth's 0.9–1.0 m and 32–51 mm across
  (Approved Document K; ADA 505). On a side with a guard, the guard's grip slot is the handrail.

## Handholds

On the Moon a handrail catches a fall as well as steadying a step (*derived*, for a 75 kg adult and 0.3 m of give):

| Grip (Dodds et al. 2014) | One hand holds, on the Moon | Fall stopped by one hand | By two hands | By two hands on Earth |
|---|---:|---:|---:|---|
| Women's weak grip, 16 kg | 97 kg | 0.09 m | 0.47 m | Cannot hold the weight |
| Men's weak grip, 27 kg | 163 kg | 0.35 m | 1.00 m | Cannot hold the weight |
| Women's peak median, 31 kg | 187 kg | 0.45 m | 1.20 m | Cannot hold the weight |
| Men's peak median, 51 kg | 308 kg | 0.93 m | 2.16 m | 0.11 m |

Rails run along every stair, ramp and lane, with grab bars at doors and landings and, in the halls for moving, grips at
several heights to catch, turn and launch from (*design choice*). Every handhold is a climbing hold for a child, so
along a guard the handhold is the grip slot in its face, and free-standing rails keep the footholds' reach from it.

## Ramps and slopes

- **Climbing costs a sixth** (*derived*). A 1:12 ramp climbs with the force a 1:73 ramp takes on Earth, and a 30°
  slope climbs like Earth's 1:12.
- **Descending keeps Earth's proportions** (*derived*). Grip falls with gravity as the slope's pull does: at 1:12 the
  slope takes 14% of the grip on either world, and a wheel that rolls a whole run unbraked needs the same level run-out
  at its foot on both, 0.138 m for each metre of a 1:12 run at the full grip. A landing of length L stops such a
  runaway from any run rising up to μL, whatever its slope: 0.6 L at the full grip and 0.3 L at half of it. Stopping on
  the ramp itself from walking speed takes 1.2 m at the full grip, against 0.2 m on Earth.
- **Steep slopes ask for grip both ways** (*derived*). Standing or walking steadily on a slope takes tan θ/μ of the
  grip, uphill or down and at any gravity: 55% at 18.4°, 78% at 25° and 96% at 30°, the slope that climbs like Earth's
  1:12. Going down, a stop from walking speed at the full grip takes 2.4, 5.1 and 32 m.
- **Ramps** (*design choice*) keep Earth's gradients, 1:12 to 1:20, and rise up to 0.9 m a run between 3 m landings,
  which stop a runaway at half the grip. The runs depart from Approved Document K, whose diagram 2.1 allows 1:12 only
  for flights of up to 2 m of going and 1:20 up to 10 m, 0.5 m of rise at most: here a run is 10.8 m long at 1:12,
  13.5 m at 1:15 and 18.0 m at 1:20. The climb costs a sixth of Earth's effort, which Earth's short flights between
  rests answer, and a landing stops a runaway from the same rise on either world. The ADA's 1:12 runs rise 0.76 m
  between 1.525 m landings, which stop a runaway from 0.91 m at the full grip.
- **Turns and steep paths** (*design choice*). At half the grip a level turn at walking speed takes a radius of 4.1 m,
  so ramps run straight or turn on wide landings, and paths steeper than 16.7°, where steady walking takes more than
  half the grip, take steps or high-grip surfaces.

## Halls for moving

The register designs glider and parkour zones into the summit metropolis
([register](../../research/decisions.md#the-summit-port-and-metropolis)), and its
[brief](../summit_metropolis/README.md) keeps the rest of its architecture from inviting either; the scene the author
accepted on 15 October had a separate hall for moving. A hall's dimensions follow from the jump (*derived*), and the
sizes chosen are a *design choice*:
- **Height** (provisional). The arm-swing jump at Earth's take-off speed rises 2.60 m from standing: a 1.9 m adult's
  head reaches 4.50 m and raised arms about a metre more, 5.5 m, so the halls stand 6.0 m clear for jumps from the
  floor, and higher below platforms. The same leg work would carry the jump 4.4–4.9 m from standing and ask 7.3–7.8 m,
  so tests of jumping at lunar gravity set the height.
- **Length.** A running hop at 3 m/s with the countermovement jump's take-off speed flies 9.7 m in 3.2 s, against
  1.6 m on Earth; at the same leg work's 3.6–3.8 m/s, 13–14 m.
- **Turns.** At half the grip a level turn at 3 m/s takes a radius of 18.5 m. Banked floors tighten it: 31° of bank
  for a 5 m radius at 3 m/s, 45° for 3 m, and 34° for walking speed round 1 m.
- **Landing.** A jump lands at its take-off speed, as on Earth, so sports floors serve.

## The visual criteria

The author's criteria: the guards and circulation show the Moon's gravity and the added spring in a step, in forms
distinct from Earth's. The scene accepted on 15 October carried continuous tall guards, inhabited spans with visible
support, human-scale spaces and a separate hall for moving, and the author found it suited the semi-arid parks,
gardens and their amenities, its setting left provisional
([register](../../research/decisions.md#the-chatgpt-work-of-15-october-confirmed-17-october-2026)). The figures above
give these forms (*design choice*):
- **Tall, continuous guards** of 2.2 m along every drop over 4.7 m, clear or finely open, their faces plain but for a
  grip slot and their tops curving inward; along drops of 2.3–4.7 m, Earth's waist-high rail or the same 2.2 m.
- **Open edges.** Level changes under 2.3 m carry no guard: terraces, stepped gardens, plinths and quays that people
  step and hop down.
- **Long, low stairs.** Deep treads in flights rising up to 12 m between deep landings, with bounding stairs of
  knee-high steps beside them.
- **High openings.** Door heads and stair headroom of 2.5 m.
- **Rounded and banked circulation.** Corners on walking routes rounded to 4.1 m or more, and banked lanes where people
  run.
- **Handholds everywhere**, shaped to be seen and used: rails along the routes, grip slots in the guards.
- **Stepping stones for a hop.** A hop at 1.4 m/s along and 1.0 m/s up, the take-off of a 5 cm hop on Earth, flies
  1.7 m (0.3 m on Earth), so stones 1.5–1.8 m apart across the parks' rills and mirror pools make a lunar crossing.
- **Halls for moving**, tall and long, apart from calm rooms at human scale.

## The summit drawings

The ring 0 drawings carry the lunar guard already: 2.2 m railings along both edges of the roof park; the lanes' voids
are railed and bridged, at a height the drawings leave open. Two elements of their section stand lower, beside long
drops: the inner gallery's balcony over the hollow at 1.2 m, and the galleries' balustrades along the 9 m nave at
1.1 m ([draw.py](../../immersion/world/summit/ring0/draw.py), lines 413 and 416). Under the bands above both rise to
2.2 m in the drawings' next revision.

## Evidence and boundaries

The figures are closed-form mechanics with the cited codes and measurements, each assumption stated beside the figures
it bounds. The measurements of people in low gravity come from short simulations (harness, suspension and parabolic
flight), and the longest human stay on the Moon is 75 hours. Earth's codes carry Earth's experience of falls; a lunar
fall lasts six times as long, which may let a tumbling body turn further before it lands, and testing that is open.
The jumps and grips are adults'. The guard bands, the footholds' reach, the stairs, landings, ramps and halls are
design choices for the author's review and for human-factors tests in low gravity. Where gravity is Earth's, in the
array's wheels, Earth's codes apply.

## Reproduction and checks

```sh
python -m habitation.low_gravity.movement                  # results/low_gravity.json, at once
OPENBLAS_NUM_THREADS=1 python -m pytest -q habitation/tests/test_low_gravity.py
```

The tests check that the product regenerates and binds its code and the shared constants; that the cited values match
sources.json; that hops rise and last six times as long; the belt jump's unscaled extension; the same leg work's bound
by its energy; falls with drag against a stepped integration and against the daily-life study; the hop rule's equal
margins and the energy rule's equal lift, each from its inputs; the guard bands against the codes' drops; the
footholds' reach against a search over launch angles; the stairs, landings and stops at half the grip; the handholds;
that a ramp's runaway run-out does not depend on gravity, and the runs' departure from Approved Document K; the banks'
force balance at half the grip; and the halls. On 17 October: 18 passed.
