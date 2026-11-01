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

**The body rises fast and sinks slowly.** Lifting is cheap, so the legs, pushing a sixth of the weight, leave the floor
at 3.0 m/s where they leave at 2.4 m/s on Earth (Samozino's force-velocity model with a profile measured in women
athletes, whose jumps rise 0.29 m): the dancer's highest jump rises 2.9 m and lasts 3.8 s. Sinking has a ceiling: the
centre of mass can never accelerate downward faster than g, so a 0.35 m grand plié takes at least 0.93 s to descend and
stop, where the rise back up can take a third of that.

**Flights become sentences.** A jump at half the dancer's height (1.4 m, 2.7 s) holds three or four shapes of
0.7-0.8 s each, where Earth's 0.5 s holds one. A layout somersault over that flight needs 22 N m s of angular momentum,
a fifth of Earth's 121 for its best flight: the somersault turns slowly enough to read as a shape turning. With angular
momentum fixed in flight, drawing the arms from second to first speeds a twist 2.6 times, a change of tempo inside one
jump that the eye can follow.

**Grip is weak, so momentum is precious.** With a dance shoe on vinyl (friction 0.5) the floor gives 0.8 m/s^2 of
horizontal acceleration where Earth gives 4.9: building or spending 2 m/s takes 2.5 m, and the tightest level turn at
that speed has a radius of 4.9 m. Bounding strides and landings steer instead: a stride's push and landing press the
floor hard enough to change the speed by about 0.8 m/s, and a landing at 2.5 m/s absorbed over 0.4 m presses with 5.8
lunar body weights and can turn 1.5 m/s of travel. Lunar travel bounds, curves in wide orbits and changes direction
through its landings.

**Turns start slowly and keep going.** The floor's grip gives a turn 8.9 N m of torque against Earth's 54, so two
revolutions a second need about 0.8 s of preparation; once turning on demi-pointe, pivot friction stops the spin only
after about 10 revolutions (under 2 on Earth). A rigid strip of skirt would flare to 45 degrees at 18 rpm (45 rpm on
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
them. [vocabulary/units.json](vocabulary/units.json) holds 29 units in eight families (travel; rise, fall and balance;
turns; flight; floor; upper body; fabric; partner). Each names its source tradition, its entry and exit states, its
phases with physics-derived durations and key shapes, its Laban effort, what it does to the garment, and what the Moon
changes, in figures. Seventeen are Earth units adapted; twelve are native to partial gravity in air: the orbit run,
skim, flight phrase, slow salto, unfolding turn, carve, long fall, long penché, whirl and canopy for a soloist, and a
one-arm lift and a toss for partners, which wait on a two-body model.

[grammar.py](grammar.py) places units into phrases. A unit enters only from a state it fits; otherwise the grammar
inserts the linking step a dancer would use, sized by the physics (a weight change, strides that build or spend speed,
a braking step for a spin, a wait for the sleeves to come within reach, a rise from the floor), and takes a unit on the
other side when only the supporting foot is wrong. Choreographic devices act on motifs (mirror, augmentation,
diminution, accumulation, retrograde of order), and a generator draws phrases from family transition weights for
comparison. Each phrase gets a count score and an evaluation of what a choreographer looks at: length, variety, effort
contrast, levels, time aloft, the lunar-native share, travel, curvature, stops, linking time and where the climax
falls. [results/phrases.json](results/phrases.json) holds three: *Under the Canopy* (29.6 s, composed by hand: a
port de bras, a canopy of sleeves, a long fall under it, an orbit run, a suspended grand jeté, a carve, a whirl and a
spiral to the floor), *Suspension Study* (a motif, mirrored and augmented) and the best of 40 generated phrases.

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
score into a multiple-shooting problem in Crocoddyl at 40 nodes a second: 6D contacts for the feet on the floor with a
wrench-cone cost (friction, centre of pressure and yaw moment inside each foot's limits), free fall in flight so that
angular momentum is conserved and the centre of mass follows its arc, an impulse at each landing, joint torques boxed
by the body's strength, ranges as barriers, the shapes as posture targets, the travel as a centre-of-mass reference and
turns as the pelvis's yaw. Inverse kinematics lays the soles flat on their footsteps for the starting state and the
warm start. [replay.py](motion/replay.py) runs the optimised torques, with the solver's feedback gains, through
MuJoCo's soft contacts and fluid forces as an independent check.

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
solver ([garment/](garment/)) and the renderer ([visualization/dance](../../visualization/dance/)) are under way.

## Evidence

The envelope is closed-form and one-dimensional mechanics for one representative body with stated parameters; the
vocabulary's durations and limits are screens from it. The body's segment parameters follow de Leva (1996) scaled by
stature, from a transcription of his Table 4 checked two ways (sources.json); the printed page is still to be read.
Joint placements are modelling choices; the ranges and strengths a source could check agree with dancer and
young-women data (Almasi et al. 2025; Harbo et al. 2012), and the rest (hip flexion and abduction ranges, spine ranges,
trunk and hip-rotation torques, shoe-floor friction) are open in sources.json. Two sections are realised with whole-body optimal control, satisfying the rigid-body equations of motion with rigid contacts, friction cones and torque bounds; their replay through MuJoCo's softer contacts holds for under a second, so their robustness to another contact model is unproved. No motion has been tested with dancers; the artistic worth of any phrase is for the author and dancers to judge.

## Reproduction and checks

```sh
python -m humanities.dance.body.model       # body.json and dancer.xml
python -m humanities.dance.envelope         # results/envelope.json
python -m humanities.dance.grammar          # results/phrases.json and the scores
OPENBLAS_NUM_THREADS=1 python -m pytest -q humanities/tests
```

The tests check that the products regenerate and bind their code, constants and air; Samozino's take-off against its
own equation and the brackets around it; the body's mass, centre of mass and pose ordering of its inertias; traction
and toppling against their closed forms; drag flights against the vacuum and terminal speeds against weight; fabric
falls and the skirt's flare against force balances; partnering and the leap's stage; and the body model's MuJoCo load
where MuJoCo is installed (the dance environment, `/media/projectspace/terluna-research/venvs/dance`).
