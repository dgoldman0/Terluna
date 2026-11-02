# Dance on the Open Moon

**1 November 2026.** The aim is flowing, vibrant dance that the Moon alone allows: at a sixth of Earth's gravity, in
air as dense as Earth's, with a floor underfoot and a down direction to fall toward. This folder holds the dancer's
body, the physics of dancing there, the movement vocabulary and grammar that compose phrases, and the motion engine
that realises them; the garment solver sits in [garment/](garment/) and the renderer in
[visualization/dance](../../visualization/dance/). The work starts from two bundles of earlier sessions, assessed in
[intake.md](intake.md); the author's directions are in the decisions register, section "Dance and the movement arts,
1 November 2026"; the [plan](plan.md) and the [open questions](questions.md) carry what comes next.

## What the Moon does to a dancer

[envelope.py](envelope.py) works out, for one body on both worlds and against microgravity, what a dancer can do at
the shared lunar gravity (1.6242 m/s^2) in the design run's air (1.40 kg/m^3 at the sea-level datum, 1.17 at the summit
metropolis). The figures below are its product, [results/envelope.json](results/envelope.json); each rests on the
assumptions its docstring states beside it.

**Gravity's clock slows; the muscles' clock does not.** Every gravity-driven time stretches by sqrt(6.04) = 2.46: a
standing body tips from 2 to 15 degrees in 2.0 s (0.8 s on Earth), a flight that rises 1.5 m lasts 2.7 s, a leg swings
like a slower pendulum. Muscle-driven motion keeps Earth's tempo: an arm whips, the feet patter and a push leaves the
floor as fast as on Earth. Lunar dance can layer quick articulation inside slow arcs, and its phrasing has two clocks.

**Weight moves slowly from foot to foot.** Standing on both feet, the dancer can carry her centre of mass across onto
one only as fast as the floor pushes it, from the outer edge of the far foot and then of the near one, and gravity sets
that pace: a 10 cm shift in place takes at least 1.7 s on the Moon against 0.7 s on Earth (an inverted pendulum
x'' = (g/h)(x - p) with the push at the support's edges). Every change from two feet to one without a step carries this
cost, so the grammar times it as its own step.

**The body rises fast and sinks slowly.** Lifting is cheap, so the legs, pushing a sixth of the weight, leave the floor
at 3.0 m/s where they leave at 2.4 m/s on Earth (Samozino's force-velocity model with a profile measured in women
athletes, whose jumps rise 0.29 m): the dancer's highest jump rises 2.9 m and lasts 3.8 s. Sinking has a ceiling: the
centre of mass can never accelerate downward faster than g, so a 0.35 m grand plié takes at least 0.93 s to descend and
stop, where the rise back up can take a third of that. A spiral to the floor lowers the centre of mass 0.22 m into a
lunge and 0.25 m more into a kneel; at half of gravity's pull these take 1.0 and 1.1 s, against 0.4 and 0.45 s on
Earth.

**A fall asks the same grip on both worlds.** A dancer falling forward as one piece slows, but the ratio of the floor's
sideways push to its upward push is set by the angle and the body's shape alone. Her own moment of inertia keeps it
below a point mass's tan(angle): 0.49 at 35 degrees and 0.50 near 40, a dance floor's whole grip; past 55 degrees
the floor would carry no weight at all. The long fall therefore catches before 50 degrees, and slowly: from 2 to
35 degrees takes 2.7 s, where Earth's takes 1.1.

**Flights become sentences.** A jump at half the dancer's height (1.4 m, 2.7 s) holds three or four shapes of
0.7-0.8 s each, where Earth's 0.5 s holds one. A layout somersault over that flight needs 22 N m s of angular momentum,
a fifth of Earth's 121 for its best flight: the somersault turns slowly enough to read as a shape turning. With angular
momentum fixed in flight, drawing the arms from second to first speeds a twist 2.6 times, a change of tempo inside one
jump that the eye can follow.

**Rotations in the air become figures.** Posed in the vocabulary's shapes ([body/pose.py](body/pose.py)), the body
model shows how much a take-off can turn and what each figure takes. Leaning 0.25 m past the foot, a lunar take-off to
half height gives a somersault or a cartwheel about 37 N m s, and the floor's grip on its push gives about 22 about the
vertical; Earth's best standing take-off gives 74 and 44, but for a flight five times shorter. Five figures follow
that no standing jump on Earth allows:
- **Corkscrew.** A layout somersault with the arms gathered needs 18.5 N m s, and a tilt of 6 degrees out of its
  plane, made by lowering one arm, turns its spin into one full twist at 0.38 revolutions a second.
- **Star wheel.** A star spread wide turns over sideways once on 24 N m s, steady because it turns about its largest
  axis.
- **Flat spin.** Laid out flat and face down, the same spread spins once about the vertical on 16 N m s from a 2.1 m
  rise, within the grip's 25 (from 1.4 m it would need 19.6 of 22). It comes from the butterfly twist of Chinese
  martial arts, over in half a second on Earth.
- **Hanging salto.** Tucked, the body turns 4.0 times faster than laid out long, so a somersault that tucks, opens
  upside down and tucks again hangs inverted for 1.06 s, drifting through 66 degrees.
- **Sleeve flight.** Sleeves thrown up at the take-off are left standing up to 1.1 m above the hands as the body falls
  back, and they come down for 5 s after the landing; Earth's lift 0.14 m.

Height and distance open more, which wait on a future study of levels: a 2.9 m rise clears a standing partner or
reaches a 2 m ledge, and a drop from 3 m lands at 3.1 m/s, where Earth's lands at 7.7.

**Grip is weak, so momentum is precious.** With a dance shoe on vinyl (friction 0.5) the floor gives 0.8 m/s^2 of
horizontal acceleration where Earth gives 4.9: building or spending 2 m/s takes 2.5 m, and the tightest level turn at
that speed has a radius of 4.9 m. Bounding strides and landings steer instead: a stride's push and landing press the
floor hard enough to change the speed by about 0.8 m/s, and a landing at 2.5 m/s absorbed over 0.4 m presses with 5.8
lunar body weights and can turn 1.5 m/s of travel. Lunar travel bounds, curves in wide orbits and changes direction
through its landings.

**Turns start slowly and keep going.** The floor's grip gives a turn 8.9 N m of torque against Earth's 54, so two
revolutions a second need about 0.8 s of preparation; once turning on demi-pointe, pivot friction stops the spin only
after about 10 revolutions (under 2 on Earth). A slow whirl is another matter, since the revolutions friction allows
grow with the square of the spin: at 0.45 revolutions a second it stops an unpushed whirl within 1.3 revolutions with the
arms open in second (half a revolution with them gathered), so a whirl is kept going by a small push of the free foot
each revolution, as Mevlevi turning does. A rigid strip of skirt would flare to 45 degrees at 18 rpm (45 rpm on
Earth), but for a light skirt the air decides instead: in the garment solver a chiffon circle skirt turned steadily at
5-20 rpm winds up around the body, its hem lagging 18-95 degrees, and opens only 10-13 degrees, against 51-60 with the
air removed; satin, 3.4 times heavier, opens 16-23. A skirt for whirling needs weight in its fabric or hem.

**The air is a partner to the cloth and nearly absent to the body.** Air outweighs silk chiffon (26 g/m^2) 19 to 1 at
1 m/s on the Moon, against 3 to 1 on Earth; chiffon falls face-on at 0.23 m/s and takes 6.6 s to settle 1.5 m. In the garment
solver, water sleeves flung up as the arms rise to fifth hover level above the head while Earth's fall back to hang
from the hands, and in the realised leap the skirt opens into a bowl through the strides, spreads into a flat canopy
in the flight and blooms upward on the landing. The same air barely touches the body's flight: even a 2 m^2
spread of sleeves or cloak lowers a 1.4 m jump by about 3% and shortens it by about 2%. Ordinary garments shape the image,
and the body's arc stays its own.

**Partners weigh a sixth.** A 55 kg partner weighs 89 N, within a one-arm overhead hold, and a two-arm push of 350 N over
0.6 m throws them 1.75 m up for 2.9 s, where on Earth the same push cannot lift them.

**What microgravity lacks.** Without weight there is no grip, no landing, no fall and no direction for cloth to hang:
every push sets an endless drift. The Moon keeps the rhythm of weight and floor and gives it time.

## The vocabulary and the grammar

[vocabulary/poses.json](vocabulary/poses.json) holds named shapes by body region (legs, arms, torso, head, feet) after
ballet, modern, ballroom, acrobatics and water-sleeve teaching, combined into whole shapes the way a dancer combines
them. [vocabulary/units.json](vocabulary/units.json) holds 35 units in eight families (travel; rise, fall and balance;
turns; flight; floor; upper body; fabric; partner). Each names its source tradition, its entry and exit states, its
phases with physics-derived durations and key shapes, its Laban effort, what it does to the garment, and what the Moon
changes, in figures. Eighteen are Earth units adapted; seventeen are native to partial gravity in air: the orbit run,
skim, flight phrase, slow salto, unfolding turn, corkscrew, star wheel, flat spin, hanging salto, sleeve flight, carve,
long fall, long penché, whirl and canopy for a soloist, and a one-arm lift and a toss for partners, which wait on a
future study of partnering with a two-body model.

[grammar.py](grammar.py) places units into phrases. A unit enters only from a state it fits; otherwise the grammar
inserts the linking step a dancer would use, sized by the physics (a weight change, strides that build or spend speed,
a braking step for a spin, a wait for the sleeves to come within reach, a rise from the floor), and takes a unit on the
other side when only the supporting foot is wrong. Where units meet, it keeps one stance to a foot: a carve replaces
the landing of the jump before it, and a run that starts on the foot already down uses that stance as its first push
(adding a stride when the run must change speed). It never stands one foot through two stances, which at running speed
would carry the body a metre and more over a planted foot. Choreographic devices act on motifs (mirror, augmentation,
diminution, accumulation, retrograde of order), and a generator draws phrases from family transition weights for
comparison. Each phrase gets a count score and an evaluation of what a choreographer looks at: length, variety, effort
contrast, levels, time aloft, the lunar-native share, travel, curvature, stops, linking time and where the climax
falls, with the revolutions turned in the air about each axis. [results/phrases.json](results/phrases.json) holds
four: *Under the Canopy* (33.4 s, composed by hand: a port de bras, a canopy of sleeves, a weight shift and a long fall
under them, an orbit run, a suspended grand jeté, a carve, a run that spends the speed, a whirl and a spiral to a
kneel), *Suspension Study* (a motif, mirrored and augmented), *Air Study* (27.4 s, composed by hand from the aerial
figures, 56% of it in the air: a run into a star wheel, a hanging salto and a corkscrew on the spot, a sleeve flight,
the highest jump laid out flat and spinning, and a spiral down to a kneel) and the best of 40 generated phrases.

## The body, the motion engine, the garments and the renderer

[body/model.py](body/model.py) describes the dancer once and writes [body/body.json](body/body.json) and the MuJoCo
model [body/dancer.xml](body/dancer.xml): 21 segments and 48 degrees of freedom, ball-jointed hips and shoulders for
legs and arms that range over the whole sphere, a three-part spine, shoulder girdles, forearm rotation, wrists,
subtalar joints and toes for demi-pointe, sole points for contact, collision capsules for cloth, and joint strength
for the torque bounds. [motion_io.py](motion_io.py) fixes the motion file every consumer reads: segment world poses
through time.

[motion/](motion/) realises a phrase with the dancer's full rigid-body dynamics. [pinmodel.py](motion/pinmodel.py)
builds the body in Pinocchio from body.json, depth-first so that each subtree's velocities form one block, and maps
coordinates to and from MuJoCo; the two programs agree on positions and velocities to rounding error, on kinetic energy
to 1e-7 and on the joint accelerations the same torques produce to 1e-8. [realise.py](motion/realise.py) turns the
score into a multiple-shooting problem in Crocoddyl at 40 nodes a second: 6D contacts for feet flat or on demi-pointe
with a wrench-cone cost (friction, centre of pressure and yaw moment inside each foot's limits), free fall in flight so
that angular momentum is conserved and the centre of mass follows its arc, an impulse at each landing, joint torques
boxed by the body's strength, ranges as barriers, the shapes as posture targets, the travel as a centre-of-mass
reference and turns as the pelvis's spin. A turning foot pivots on a point under the ball that cannot slip but can
swivel, with its metatarsal heads pressing the floor without friction so that the foot cannot roll and the floor's
push can move across the ball. The floor's turning friction acts on that swivel as a moment (the envelope's spin
friction at full lunar weight, a Coulomb law smoothed below 1 rad/s), so a whirl slows unless the free foot pushes. A
kneeling knee holds the floor at a point on the front of the shank. The footholds are planned for the whole phrase
before any solving: a foot landing from the air goes under the centre of mass at mid-stance, a step or a kneel goes
where its shape puts it measured from the foot already down (the centre of mass moving over with it), and a turn's
preparation sets the turning foot under the centre of mass with the other behind it in fourth. Inverse kinematics
places the contacts and the centre of mass for the starting state and the warm start. Each node's model takes about
12 MB, so a whole phrase is solved in overlapping 3.5 s windows built one at a time (about 1.8 GB), keeping only each
node's state, torques, gains and contact forces. [replay.py](motion/replay.py) runs the optimised torques, with the
solver's feedback gains, through MuJoCo's soft contacts and fluid forces as an independent check.

What it has done so far. A plié and a suspended sauté (4.2 s) solve in 150 iterations with no dynamic gap: the take-off
push ramps from 120 to 800 N over 0.3 s, the body rises 0.45 m, the feet use friction up to 0.47 against the 0.5
available, and the joints at most 44% of their strength. The run, leap and carve at the centre of *Under the Canopy*
(6.5 s, 260 nodes) solve in 250 iterations with no dynamic gap: three bounding strides on a 4 m circle, a last stride
that sinks the centre of mass, a take-off whose push reaches the 900 N cap, a flight of 2.4 s that lifts the centre of
mass about 1.2 m and carries it 4.6 m through the split, and a landing at 6.0 s. The feet use a median friction of
0.19, reaching the 0.5 available at the hardest moments, and the joints at most 87% of their strength. A barrier keeps
every sole point that bears no weight above the floor; the renderer's review found the first solution's free toes
27 cm through the floor, and the barrier leaves a toe dipping at most 1.1 cm at the take-off for a few frames. The MuJoCo replay holds the plan for the first 0.4 s
of the sauté and then drifts: the optimiser's rigid contacts and MuJoCo's soft ones part, and the feedback gains alone
do not close the difference; a tracking controller or model-predictive control is the replay's next step. The garment
solver ([garment/](garment/)) and the renderer ([visualization/dance](../../visualization/dance/)) dressed and
rendered the leap. That leap was solved before the grammar's junction rules: its jeté landed, then carved on the same
foot, and its window ended 0.4 s into the carve, on one foot with the centre of mass travelling past it, which the
render shows as a fall.

The whole of *Under the Canopy* (33.4 s and a 1 s held finish, 1,375 nodes) is now planned end to end, and each new
piece has been solved on its own from the plan's state, for 80 iterations at most (light testing). Every piece reaches
zero dynamic gap. The fall (3.8 s) converges in 49 iterations, rolling onto the ball of the foot as its centre of mass
passes over it; friction reaches the floor's 0.5 only at the end of the weight shift, where the test begins, and the
joints use at most 32% of their strength. The jeté's push, flight and carve into the run (4.3 s) reaches the cap short
of convergence, friction passing the 0.5 its cone allows at two nodes (0.62 at the carve's landing, 0.52 in the run)
and the joints at most 76% of their strength. The whirl's preparation, pivot and
first push (2.8 s) spins at 0.43-0.45 revolutions a second against the plan's 0.45, the joints at most 26% of their
strength. The spiral, lunge, kneel and finish (4.2 s) turns at 0.46-0.5 revolutions a second and comes to rest
kneeling. These two also stop at the cap, with friction peaking at the floor's 0.5. Floor points dip at most 1-3 cm in these
short solves. Replayed through MuJoCo with the solver's gains, the fall holds its whole 3.8 s, the centre of mass within
1.3 cm of the plan and the pelvis within 3 degrees; the whirl holds about 1 s, then drifts once the pivot starts, since
MuJoCo's soft sphere contacts model the swivel and its friction differently, and falls at 2.7 s. The full solve, its
replay and its render wait for machine time. The engine turns the body in the air only about the vertical; the
somersaults, cartwheels and flat spins of the aerial units need orientation targets for the flight, still to come.

## Evidence

The envelope is closed-form and one-dimensional mechanics for one representative body with stated parameters; the
vocabulary's durations and limits are screens from it. The body's segment parameters follow de Leva (1996) scaled by
stature, from a transcription of his Table 4 checked two ways (sources.json); the printed page is still to be read.
Joint placements are modelling choices; the ranges and strengths a source could check agree with dancer and
young-women data (Almasi et al. 2025; Harbo et al. 2012), and the rest (hip flexion and abduction ranges, spine ranges,
trunk and hip-rotation torques, shoe-floor friction) are open in sources.json. The aerial figures take the body model's inertia in each shape (the computed tuck and pike shares, 0.25 and 0.50 of the layout's, sit a little below Yeadon's 0.30 and 0.55 that the older rotation figures assume) and an assumed lean of 0.25 m and spin lever of 0.30 m for what a take-off can turn; they are screens of angular momentum, not solved motions. Two sections are realised with whole-body optimal control, satisfying the rigid-body equations of motion with rigid contacts, friction cones and torque bounds; the sauté's replay through MuJoCo's softer contacts holds for under a second and the fall piece's for its whole 3.8 s, so robustness to another contact model is shown only for the fall. The rest of the first phrase is tested only in short solves capped at 80 iterations, three of them short of convergence and with friction costs a little past the floor's limit at single nodes. The pivot's turning friction assumes the full body weight on the ball (it overstates the friction while the other foot pushes) and leaves the frictionless metatarsal supports without any. No motion has been tested with dancers; the artistic worth of any phrase is for the author and dancers to judge.

## Reproduction and checks

```sh
python -m humanities.dance.body.model       # body.json and dancer.xml
python -m humanities.dance.envelope         # results/envelope.json
python -m humanities.dance.grammar          # results/phrases.json and the scores
OPENBLAS_NUM_THREADS=1 python -m pytest -q humanities/tests
# the dance environment: the whole first phrase in windows, ending held still for 1 s (about 1.8 GB at a time)
OMP_NUM_THREADS=1 python -m humanities.dance.motion.realise canopy --windows --finish 1 --tag under_the_canopy
# or one piece of it: --span 21.7:24.5 (seconds of the phrase) solves the whirl's preparation and first push
```

The tests check that the products regenerate and bind their code, constants and air; Samozino's take-off against its
own equation and the brackets around it; the body's mass, centre of mass and pose ordering of its inertias; traction
and toppling against their closed forms; drag flights against the vacuum and terminal speeds against weight; fabric
falls and the skirt's flare against force balances; partnering and the leap's stage; weight shifts against their
gravity scaling and the fall's grip against the floor's; that no phrase stands one foot through two stances; the
aerial figures against their mechanics (the corkscrew's tilt, the hanging salto's programme turning exactly once, each
lunar figure within its take-off and each standing figure beyond Earth's); the posed body against the standing one and
the Air Study's turns in the air; and the
body model's MuJoCo load where MuJoCo is installed (the dance environment,
`/media/projectspace/terluna-research/venvs/dance`). There the engine's tests also check the whole phrase's contact
schedule (turns on a pivot, two planted feet never asked to turn past what the hips allow, a push each revolution,
a still kneel at the end), its footholds, the pivot friction's direction and derivatives, and that the posed body
places every joint where the engine's Pinocchio model does.
