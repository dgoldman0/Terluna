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
five: *Under the Canopy* (33.4 s, composed by hand: a port de bras, a canopy of sleeves, a weight shift and a long fall
under them, an orbit run, a suspended grand jeté, a carve, a run that spends the speed, a whirl and a spiral to a
kneel), *Suspension Study* (a motif, mirrored and augmented), *Air Study* (27.9 s, composed by hand from the aerial
figures, more than half of it in the air: a run into a star wheel, a hanging salto and a corkscrew on the spot, a sleeve flight,
the highest jump laid out flat and spinning, and a spiral down to a kneel), *Under the Canopy, with the hanging salto*
(35.7 s, the longer phrase of 2 November) and the best of 40 generated phrases.

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

The whole of *Under the Canopy* (33.4 s and a 1 s held finish, 1,375 nodes) was solved on 2 November in 32 overlapping windows of 3.5 s, 12 minutes on one core, with no dynamic gap anywhere, and it passes the checks that gate a render. Friction reaches the floor's 0.5 at the hardest moments (the fall's catch, the jeté's take-off and the carve, the whirl's preparation and the spiral) and never passes it; every contact presses on the floor (at least 1 N); the joints use at most 96% of their strength, in the carve; no part of the body but a sole point goes through the floor, and the soles dip at most 2.3 cm. The whirl turns at 0.43-0.47 revolutions a second against the plan's 0.45 and the spiral at 0.48. Soft limits are passed briefly: the left ankle bends 7.9 degrees past its 25 degrees of dorsiflexion for a quarter of a second in the jeté's take-off lunge (a dancer's demi-plié reaches 30-35), and the trunk passes its ranges by 1-3 degrees. Twelve of the 32 windows stop at the 200-iteration cap, most of them in the whirl and the spiral. The MuJoCo replay was not run again, since the replay needs a tracking controller first (plan, stage 2).

**Figures in the air.** The engine now turns the body over in flight. An aerial figure's flight carries a planned attitude for the pelvis and a planned rate of turn: in the air the body keeps its angular momentum L and turns at L / I, with I its moment of inertia about the turning axis in the shape it holds, so the planned angle grows with the integral of dt / I through the figure's shapes, and the take-off has to give the momentum the whole turn needs. A twist is spread over the middle of the flight about the body's long axis, and a flat spin lays the body over face down and spins it about the vertical. Each figure was solved as a short phrase of its own, with the steps the grammar adds to reach it and a 1.5 s held finish:
- **Hanging salto** (rise 1.6 m): zero dynamic gap at the 250-iteration cap. The take-off gives 16.8 N m s against the plan's 13.9; the body turns 0.97 of a backward revolution, the pelvis within 6 degrees of the planned attitude (median) and 16 degrees at the landing; friction at most 0.37, joints at most 83% of their strength.
- **Sleeve flight** (rise 1.4 m, in place): zero dynamic gap at the cap; friction reaches 0.5 and one joint its strength bound. It lands pitched far forward, the trunk nearly level. Its first solve, travelling at 0.6 m/s, put the left hand 71 cm through the floor after the landing, because the engine kept only the soles and knees above the floor; the palms are now held 3 cm above it, and in the final solve they come within 1 cm of it as the dancer lands deep in the lunge.
- **Star wheel** (rise 1.4 m from a 0.5 m/s run): zero dynamic gap, but the take-off gives 17.6 of the 24 N m s the cartwheel needs, so the body turns 0.82 of a revolution and lands 41 degrees short of upright. It was not rendered. Placing the take-off foot 0.25 m beside the body's lean, as the envelope assumes, did not change the result, and neither did a take-off stance 0.3 s longer (18 N m s, 0.84 of a revolution, in a test outside the repository), so the take-off's time is not what limits it; one joint already works at its strength bound.
- **Corkscrew and flat spin:** the solver stops after 29 and 41 iterations with dynamic gaps of 0.47 and 0.57, their take-offs giving almost none of the planned momentum (2 of 18 N m s; 1 of 15). Neither was rendered. Their analytic derivatives match finite differences at every node. A staged solve of the corkscrew (the somersault alone first, then the twist warm-started from it) reaches no dynamic gap but turns only a quarter twist, because its shapes keep both arms level: the twist the unit describes comes from lowering one arm, which the vocabulary's shapes do not yet hold.
A first version of the momentum target used Crocoddyl's centroidal-momentum residual, whose gradient in this setup disagreed with finite differences (0.19 against 0); the rate of turn of the pelvis took its place. Solve windows never cut an aerial figure: a window that reaches into one lengthens to hold it to its landing.

**Two junctions corrected.** The grammar drops a unit's first stance when the dancer already stands on that foot (a run that starts on the foot already down uses that stance as its first push). It had treated the end of a weight shift, with both feet still on the floor, as such a stance, and so took away the sleeve flight's push in the *Air Study*; a jump after a shift now keeps its own take-off. The closing step that stops the travel before a unit starting on both feet landed, out of a run's last flight, on the foot that had just pushed, 0.75 m ahead of the centre of mass and beyond the leg's reach; the first two solves of the longer phrase failed from that step on (every later window stopped after 11 iterations with a dynamic gap). It now lands on the foot the run comes down on. The *Air Study*'s flat spin still takes off within the sleeve flight's 0.38 s landing, a rebound from 2.1 m/s down to 2.7 m/s up that asks more push than the foot's 900 N bound allows; the study needs a step there before it is solved.

**The longer phrase.** *Under the Canopy, with the hanging salto* (35.7 s and a 1 s held finish; [grammar.py](grammar.py), `canopy_and_salto`) keeps the first phrase up to its whirl and ends on the salto: after the whirl the grammar brakes the turn with a pas de bourrée and closes to both feet, and the salto lands still. The salto was first set between the carve and the whirl, where the solve failed in the whirl that followed it in each of three window settings (4.5 and 5 s long; every window after the salto's landing stopped after 11 iterations with a dynamic gap), while the whirl entered from the run solves; why a whirl cannot start from the salto's landing is still open. Windows now lengthen to hold a figure from its take-off to its landing, so the rest of the phrase keeps the 3.5 s windows that solved *Under the Canopy*. Solved so (1,469 nodes, 34 windows, 14 minutes), the phrase has no dynamic gap and passes the checks: every contact presses on the floor, the joints stay within their strength, the palms stay at least 64 cm above the floor and the soles dip at most 4.1 cm, and friction passes the floor's 0.5 only at two nodes of the whirl, at 0.503 and 0.51. The left ankle's 7.9 degrees past its range in the jeté's take-off recurs. The whirl turns at 0.33-0.48 revolutions a second, and the salto rises to 2.45 m, hangs upside down at 34.2 s and comes upright by 35.7 s.

## Fabrics and cuts on the Moon

The first study of attire (1 November) asks how fabrics, and the shapes cut from them, behave at lunar gravity in the design run's air, as the ground for designing garments for particular dances. It has three parts: a library of thirteen fabrics ([garment/materials.py](garment/materials.py), each value's source and standing in [garment/sources.json](garment/sources.json)), a closed-form screen of every fabric on both worlds ([garment/screen.py](garment/screen.py), [results/garments.json](results/garments.json)), and trial runs in the garment solver ([garment/trials.py](garment/trials.py)) of what the screen can only estimate: Cusick's drape test, circle skirts in steady turns and a panel trailing in travel. The library spans 15 to 200 g/m^2: nine silks (organza, chiffon, georgette, habotai, crepe de chine, a soft satin, taffeta, a heavy satin and velvet), a plain wool, a plain cotton, a nylon tulle and a viscose jersey. Weights come from measurements or suppliers. Bending stiffness is measured for chiffon, georgette, the wool and the cotton; for the soft satin it is inferred from a measured drape coefficient through the solver's drape test; it is a stated bracket for crepe de chine and the heavy satin, and a modelling choice for the other six, which two literature passes could not find measured in open sources (question 5). Any result that rests on their stiffness carries that choice.

**Cloth keeps the dancer's two clocks.** Left to itself (hanging, swinging, falling, settling), a garment moves on the Moon as the same garment moves on Earth, slowed by the gravity clock of 2.46. Written in the time unit sqrt(L/g), the equation of a swinging panel keeps only the fabric's mass against the air's, so the motion keeps its shape on both worlds, and only the Moon's 14% denser air changes it: a chiffon panel 85 cm long let go from 30 degrees first reaches the vertical after 4.5 s on the Moon and 1.7 s on Earth. Whatever the dancer drives at the speed of muscle (a turn, a run, an arm's sweep) meets the centrifugal pull and the air at Earth's speeds while the weight is a sixth: against the weight, the pull counts 6.0 times more than on Earth and the air's push 6.9 times more.

**Weight against the air.** Falling and trailing turn on a fabric's weight against the air, so a lunar fabric falls and trails as one 6.9 times lighter does on Earth. Lunar wool (200 g/m^2) falls face-on at 0.63 m/s and, as a straight panel trailing in travel at 2 m/s, hangs 18 degrees below the horizontal; Earth's chiffon falls at 0.60 m/s and trails at 17. Lunar chiffon falls at 0.23 m/s, takes 6.7 s to settle from 1.5 m and trails 6.5 degrees below the horizontal at 2 m/s, so it streams nearly level behind any brisk movement. Its inertia is as much the air's as its own. The distance over which a moving fabric hands its momentum to the air, 2m/(rho C_N), is 3 cm for chiffon, 11 cm for the heavy satin and 24 cm for the wool on either world, and the air a 20 cm fold carries with it weighs 8.4 times the chiffon and 1.1 times the wool. Light silk moves with its air, as smoke does, and makes the air's motion around the dancer visible. In the solver, a 1.2 m panel trailing from a moving rod holds the angles the straight-panel estimate gives to within about a degree; on the Moon the satin and wool panels flutter at their far edge by 3-5 cm at 2-3 m/s, and Earth's by 2 cm at most.

**Stiffness against weight: the drape test.** At a sixth of the weight Peirce's bending length, (B/mg)^(1/3), is 1.82 times longer, so a fabric hangs on the Moon as one with six times its bending rigidity hangs on Earth. In Cusick's test (a disc 30 cm across clamped on an 18 cm support) the solver puts the measured silks' drape coefficients at 0.30-0.32 on Earth, with 14 folds round the hem, and at 0.76-0.80 on the Moon, with 8 or 9. Lunar chiffon stands out from its support further than the measured wool and cotton do on Earth (0.59-0.60), and those two barely drape on the Moon (0.94-0.95, with one to three folds). The softest drape on the Moon belongs to the soft satin, heavier than chiffon and less stiff for its weight: 0.64 there against the 0.21 it was matched to on Earth (Lee and colleagues measured 0.219 for a silk satin of 69 g/m^2). With their chosen stiffness, organza, taffeta and tulle stay almost flat on the Moon (0.96-0.99). On full skirts the same lengthening shows as broader folds: a hanging chiffon circle skirt has 9 folds round its hem on the Moon and 14 on Earth.

**Turning skirts.** In a steady turn the weight pulls a skirt down, the centrifugal pull draws it out and the air drags back the faces of its folds. Without the air a full circle skirt of any fabric opens to 75 degrees at 27 rpm on the Moon, the lunar whirl's rate; on Earth it opens 11 degrees and needs 60 rpm to reach 71. The air's share is set by an index of the folds' drag against the centrifugal pull, which holds neither gravity nor the rate of turn: 23 for chiffon and tulle, 7 for the heavy satin, 3.0-3.7 for velvet, jersey and wool. In the solver's still air at 27 rpm on the Moon, chiffon and tulle full circles open only 13-15 degrees and wind round the body, their hems about 125 degrees behind the waist; the heavy satin opens 30 degrees, and the wool and velvet fan out to 46-47 degrees in a spiral, each hem about 100 degrees behind. A half circle meets less of that drag, because at its 30-degree cone it is already smooth: the wool half circle (index 1.7) opens to 27 degrees at 15 rpm and at 27 rpm turns as a smooth, even cone of 31 degrees with no lag at all. The heavy satin half circle (index 3.7) lags 74 degrees, and with a 50 g/m weighted hem (index 1.5) it too turns as a clean cone of 34 degrees, its hem within 16 degrees of the waist. Weighted hems open light skirts as well (a chiffon full circle with 100-200 g/m opens to 66-70 degrees), but the light cloth between hem and waist then shears far beyond what a woven silk allows (99th-percentile strain of 7-40%, even with tethers holding the hem's distance to the waist), because the solver's shear does not lock as a weave jams; those shapes are not trustworthy. On Earth at 27 rpm the full circles hang at 7 degrees (chiffon) and 15 (wool).

**What this means for attire.** Light silks (15-35 g/m^2) are the air's fabric on the Moon: they settle slowly, stream nearly level and trace the air, which suits canopies, water sleeves, veils, floats and overlays, the effects the vocabulary's fabric units use. They cannot hold a turning silhouette by themselves; in a turn they wind round the body. A flowing line in travel and falls like Earth's chiffon asks for fabrics about seven times heavier: wool, velvet and jersey (160-200 g/m^2) fall and trail on the Moon as fabrics of 23-29 g/m^2 do on Earth, and their drag lengths, six to eight times chiffon's, let them keep their momentum through a turn or a stop. A turning silhouette comes from weight against fold drag: a half circle of a heavier fabric, or a heavy satin with a modest weighted hem, turns as a clean bell at the lunar whirl, and a heavy full circle fans out in a spiral. Every fabric drapes as a stiffer one does on Earth, so crisp fabrics hold sculpted shapes with little support, and the softest drape comes from weight with little stiffness. Layers can divide these roles, a heavier underskirt or half-circle shell carrying the turning shape while a light overlay traces the air, which is the combination of light and weighted fabrics the dance's direction asks for. These are readings of the study's runs; the paragraph below gives what the first dressed runs showed.

**Garments on the dancing body (2 November).** The study's garments were dressed on the solved motions and rendered. Water sleeves need tethers to their cuffs, as the skirts have to the waistband: untethered, a cast that lifts sleeves lying on the floor stretched them by 100% within a frame (99th percentile of edge strain), and in the final kneel of *Under the Canopy* a sleeve caught between forearm and shin stretched by 130%. With each vertex tethered to its cuff at 1.05 times its distance in the pattern, the sleeves stay within 10% through the salto and the sleeve flight and within 6% through the longer phrase (whether tethers free a sleeve caught in the kneel is untested: the longer phrase ends on the salto). In the longer phrase the satin half circle stays within 7% throughout but for one moment at 17%, as the run starts out of the fall. The heavy satin half circle with its 50 g/m hem passes the grand jeté's split within 3.5% and turns in the whirl within 6.7%; the chiffon full circle stays within 6% through the salto and 3.4% in the whirl. Holding cloth pinched between two parts of the body with the body was tried and left out: it caught the satin's hem on the leading shin in the jeté and tore the skirt by 57% there and 290% in the whirl. Two limits show in the renders. A circle skirt falls over the torso of an inverted body, so the hanging salto hides its upper body while it hangs, and a skirt for inversions needs another cut (shorter, divided, or fitted at the hip). A deep landing lunge traps a full skirt between the legs: the sleeve flight's chiffon stretches past 40% after its landing, so that clip ends at the landing.

**The cloth against the rendered body (2 November).** The first dressed renders showed a thigh through the skirt and a waistband sunk into the waist. The cloth collided with body.json's capsules, which lie inside the rendered skin by up to 5 cm at the thighs, and the waistband was a ring about the abdomen capsule's axis while the rendered waist sits 5 cm forward of it, so half the ring lay up to 3.2 cm inside the waist and its back stood 8 cm off the body. Measured as each cloth vertex's signed distance to the skinned dancer at 10 frames a second, 3.1% of the whirl's chiffon lay more than 3 mm inside the skin in every frame, up to 55 mm deep. Dressed runs now collide with a dress form measured from the rendered dancer (visualization/dance/README.md) and place the waistband on its waist outline. In the whirl the share inside falls to 0.06% (at most 0.25% in a frame, 13 mm deep at the crease of the hip), and the band sits 4-8 mm off the skin all round; in the salto the chiffon's share falls from 1.8% to 0.03% and the sleeves' from 1.1% to 0.4%, while in the crouch's deep fold a third of the band lies up to 25 mm inside the folded belly. In the whirl the chiffon now flares higher than it did on the capsules; contact with the fuller hips, which turn with the body, is the likely cause and is untested. The earlier renders' half-open mouth and creased throat came from the renderer's weighting of the jaw and neck, corrected the same day (visualization/dance/README.md). The form is a stand-in until the body model and the rendered mesh are fitted as one body (question 1).

**What the study does not yet show.** The air in the solver is still: a turning skirt drags its folds through air that does not turn with it, so the drag on folds is an upper bound, and the wound shapes may be partly its product. The solver carries no added mass of the air, so its transients run too fast for light silk; the steady shapes, drapes and trailing angles reported here do not depend on it, and the flutter does. Its shear does not lock, which made the weighted light skirts untrustworthy. The skirts' 2 cm mesh is coarser than chiffon's folds, and in air the flare of light skirts changes with it (the earlier validation found 21, 11 and 8 degrees at 3, 2 and 1.5 cm for chiffon at 15 rpm on the Moon); the drape tests use 5 mm. Each turn starts from a hanging skirt spun up over 2 s; whether a skirt flicked open holds a different shape was not tried. Bending stiffness is chosen for six fabrics. No result has been compared with filmed or measured garment motion. The runs, their shapes and the figures are in /media/projectspace/terluna-research/humanities/dance-runs/garment_study_2026-11-01 (trials.json, with the first weighted-hem runs kept under "turn_untethered"), the literature passes in humanities/agents/fabrics on the same drive.

## Evidence

The envelope is closed-form and one-dimensional mechanics for one representative body with stated parameters; the
vocabulary's durations and limits are screens from it. The body's segment parameters follow de Leva (1996) scaled by
stature, from a transcription of his Table 4 checked two ways (sources.json); the printed page is still to be read.
Joint placements are modelling choices; the ranges and strengths a source could check agree with dancer and
young-women data (Almasi et al. 2025; Harbo et al. 2012), and the rest (hip flexion and abduction ranges, spine ranges,
trunk and hip-rotation torques, shoe-floor friction) are open in sources.json. The aerial figures take the body model's inertia in each shape (the computed tuck and pike shares, 0.25 and 0.50 of the layout's, sit a little below Yeadon's 0.30 and 0.55 that the older rotation figures assume) and an assumed lean of 0.25 m and spin lever of 0.30 m for what a take-off can turn; they are screens of angular momentum, not solved motions. The first phrase, the longer phrase and two aerial figures are realised with whole-body optimal control, satisfying the rigid-body equations of motion with rigid contacts, friction cones and torque bounds; many of their solve windows stop at the iteration cap without a dynamic gap, so they are feasible motions, not proven optimal ones. The sauté's replay through MuJoCo's softer contacts holds for under a second and the fall piece's for its whole 3.8 s, so robustness to another contact model is shown only for the fall. The planned attitudes of the aerial figures follow the shapes' inertia, and the figures that do not yet solve (star wheel, corkscrew, flat spin) remain screens of angular momentum. The pivot's turning friction assumes the full body weight on the ball (it overstates the friction while the other foot pushes) and leaves the frictionless metatarsal supports without any. No motion has been tested with dancers; the artistic worth of any phrase is for the author and dancers to judge.

## Reproduction and checks

```sh
python -m humanities.dance.body.model       # body.json and dancer.xml
python -m humanities.dance.envelope         # results/envelope.json
python -m humanities.dance.grammar          # results/phrases.json and the scores
python -m humanities.dance.garment.screen   # results/garments.json, the fabric and cut screen
OPENBLAS_NUM_THREADS=1 python -m pytest -q humanities/tests
# the fabric and cut trials on the GPU (about 80 minutes for all of them), then their figures
python -m humanities.dance.garment.trials OUTDIR --only drape,satin,turn,stream,hem
python -m humanities.dance.garment.figures OUTDIR
# the dance environment: the whole first phrase in windows, ending held still for 1 s (about 1.8 GB at a time)
OMP_NUM_THREADS=1 python -m humanities.dance.motion.realise canopy --windows --finish 1 --tag under_the_canopy
# or one piece of it: --span 21.7:24.5 (seconds of the phrase) solves the whirl's preparation and first push
# the longer phrase, ending on the salto (about 14 minutes)
OMP_NUM_THREADS=1 python -m humanities.dance.motion.realise canopy_salto --windows --finish 1 --tag canopy_and_salto
```

The round of 2 November ran from humanities/dance-runs/round_2026-11-02 on the research drive: run_figure.py solves each aerial figure as a short phrase and checks its flight against its plan, run_cloth.py dresses a motion in the study's garments, check_motion.py applies the whole-body checks, and render_queue.py renders the clips; the renders, logs and the base commit are there too.

The tests check that the products regenerate and bind their code, constants and air; Samozino's take-off against its
own equation and the brackets around it; the body's mass, centre of mass and pose ordering of its inertias; traction
and toppling against their closed forms; drag flights against the vacuum and terminal speeds against weight; fabric
falls and the skirt's flare against force balances; partnering and the leap's stage; weight shifts against their
gravity scaling and the fall's grip against the floor's; that no phrase stands one foot through two stances; the
aerial figures against their mechanics (the corkscrew's tilt, the hanging salto's programme turning exactly once, each
lunar figure within its take-off and each standing figure beyond Earth's); the posed body against the standing one and
the Air Study's turns in the air; a jump after a weight shift keeping its own take-off; the fabric study's sewn part-circle cones against a cone's bending energy, open
weaves, weighted hems and tethers in the solver, the screen's closed forms (the falls, the trailing angle's balance,
the swing's gravity scaling, the wind-up index), its string against the validation's and its product's regeneration;
the water sleeves' tethers to their own cuffs;
and the body model's MuJoCo load where MuJoCo is installed (the dance environment,
`/media/projectspace/terluna-research/venvs/dance`). There the engine's tests also check the whole phrase's contact
schedule (turns on a pivot, two planted feet never asked to turn past what the hips allow, a push each revolution,
a still kneel at the end), its footholds, the pivot friction's direction and derivatives, that the posed body
places every joint where the engine's Pinocchio model does, and that the aerial figures' plans turn the body once on
the momentum their shapes need (the salto's backward revolution, the corkscrew's twist, the flat spin laid over about
the vertical) while phrases without figures are planned as before.
