# The petal gown: a design brief for the first full costume

Status: a proposal of 2 November 2026, with its open choices at the end. Nothing in it is built or accepted. The numbers are closed-form estimates from the garment study's own screen, applied to the petals and sleeves by /media/projectspace/terluna-research/humanities/dance-runs/costume_design/petal_screen.py (output beside it); they are screens, which the solver trials listed at the end will test.

## What the phrase asks of a costume

*Under the Canopy, with the hanging salto* gives the costume a role in every unit (units.json, "garment"). The sleeves carry the first half: they follow the port de bras with a lag of a second or more, open into the canopy overhead and settle over the dancer for 6-8 s, stay above her as she tips into the long fall, and stream behind her in the orbit run. The skirt carries the second: it opens with the jeté's split and blooms on the descent, swings wide round the carve, makes the whirl's figure with its flare, and floats open round the inverted body in the salto.

The first dressed renders showed what a costume for this phrase must avoid. A circle skirt falls over an inverted torso and hides it while the salto hangs; a full skirt is trapped between the legs in a deep landing lunge; chiffon on a 2 cm mesh renders as crumpled tissue; the sleeves render too thick; and on the dress form a skirt dropped from a flat start dresses bunched at the hips (README, "The cloth against the rendered body"). The costume is to be elegant, and the finer meshes that elegance needs are accepted as its cost.

## The design

The gown is one long line standing still and separates into petals when it moves.

- **Bodice.** Fitted, in heavy silk satin, with a softly wrapped neckline. It moves with the body and is rendered on it, as the leotard is now, so it costs no simulation.
- **Petals.** Eight gores cut on the bias from heavy silk satin (charmeuse weight, 90 g/m2), hanging from the high hip, 0.93 m above the floor in standing, to 0.13 m above it: 0.80 m long, 16.5 cm wide at the hip (half again its 11 cm share of the 0.88 m hip, so neighbours overlap by about 2.8 cm) and 40 cm wide at the hem. Below the hip they hang free of one another. The bias lets the satin follow the hip and fall in fluid folds, and it lengthens under load, which tethers to the hip line keep within 5%. A light weighted hem tape (10-50 g/m, chosen by trial) sets how far the petals rise in the whirl. An optional second ring of eight petals in silk georgette, set between the satin ones, would move on a slower clock than the satin and leave a sheer layer lingering as the satin settles.
- **Sleeves.** Sheer sleeves fitted to the wrist, rendered on the arm, ending in the existing georgette water sleeves (1.5 m long, 0.50 m round, tethered to their cuffs), which make the canopy and stream in the run.
- **Feet.** Soft slippers; the hem leaves the feet and ankles in view, so the foot work in the motion plan shows.
- **Colour (to choose).** Deep indigo shading to pale at the petals' hems, with silver-grey georgette; or ivory and champagne. Deep indigo would stand out against the warm cyclorama; ivory would catch the stage light.

The language's cultures (questions.md, question 8) enter only as vague influences, each taken on its merits: the wrapped bodice recalls the wrapped cloth of Austronesian dress, and the indigo the indigo-dyed cloth of Hakka dress. Beyond these notes the gown is its own design.

Two other elegant directions were weighed. A column gown with a wing cape has the most classical line, but its canopy repeats the sleeves' figure and the gown falls over the torso in the salto. A layered formal gown is the richest look, but its long full skirt brings back both of the first renders' failures.

## What the screen predicts

For a petal 0.80 m long hanging from a hip radius of 0.14 m, against the phrase's run (1.6 m/s) and whirl (0.45 turns a second). The whirl's range runs from fold drag (low) to skin friction alone (high).

| Fabric | Bending length, Moon (cm) | 0.5 m fall, Moon / Earth (s) | Trail below horizontal at 1.6 m/s, Moon / Earth | Swing to vertical, Moon / Earth (s) | Whirl flare, Moon: no hem / 50 g/m hem | Whirl flare, Earth |
|---|---|---|---|---|---|---|
| silk satin, heavy (90 g/m2) | 2.7 | 1.37 / 0.53 | 15° / 38° | 2.2 / 0.85 | 42-74° / 77° | 5-13° |
| silk crepe de chine (65 g/m2) | 2.4 | 1.55 / 0.60 | 13° / 33° | 2.5 / 0.97 | 55-74° / 77° | 3-13° |
| viscose jersey (160 g/m2) | 2.4 | 1.13 / 0.44 | 20° / 49° | 1.8 / 0.68 | 46-75° / 76° | 8-13° |
| silk georgette (34 g/m2) | 2.2 | 2.03 / 0.78 | 9° / 24° | 3.4 / 1.3 | 43-73° / 77° | 0-13° |

Each figure the phrase gives the skirt comes out larger and slower on the Moon. In the run the satin petals stream nearly level behind the dancer, 15° below the horizontal where on Earth they would hang at 38°. In the jeté the dancer descends from a 1.2 m rise in 1.2 s and lands at about 2 m/s, while a lifted satin petal falls at no more than 0.42 m/s, so the petals are left above her and bloom on the way down, as the unit asks. Let go, a satin petal takes 2.2 s to swing back to vertical against 0.85 s on Earth.

In the whirl, the screen holds the study's satin half circle at the 30° cone its own hoop allows, and a full circle, free to rise to 75°, gathers its spare width into folds whose drag outweighs the outward pull 6.6 times (the screen's wind-up index), enough to wind it round the body. Petals have no spare width to fold: their 3.2 m of hem is less than the 5.7 m round a hem flared to 0.9 m. They should therefore sit near the skin-friction figure, 74° from vertical without a hem weight and 77°, nearly level, with 50 g/m; the fold-drag figure (42°) bounds them from below. On Earth the same turn lifts them 3-13°. The whirl's rate is the choreographer's control of that flare, since the outward pull grows as its square.

In the salto, the tucks turn the body at about 250° a second, and at a petal's hem 0.5-1.0 m from the axis the turn's outward pull is 6-12 times lunar gravity, so the petals should fling outward and open round the body in the tucks; in the hang at 89° a second the pull is 0.7-1.5 times gravity, and the petals should float. Separated petals leave gaps through which the inverted torso shows. These last expectations rest on the ratio of pulls alone; the solver decides them.

The bending lengths set the mesh. On the Moon the satin's is 2.7 cm and the georgette's 2.2 cm, so folds are broad and a 1.5 cm mesh resolves the satin's, with 1.2 cm for the georgette. Organza (5.5 cm) would hold crisp shapes and is kept out of the petals for that reason.

## What it needs from the solver and renderer

- An open gore (patterns.gore, added 2 November): each petal is a ring sector of 0.294 rad between radii 0.561 and 1.361 m, with its warp at 45° to its centre line for the bias.
- Shear that jams. In the first design round (2 November) bias-cut satin petals overstretched in the whirl, 230% at the 99th percentile, because along a bias-cut petal only the weave's shear holds its length and the solver's shear keeps yielding where a woven fabric jams at large shear angles; on the straight grain the warp carries the length. Bias petals wait for shear that locks, and the round's petals are cut on the straight grain.
- The petals pinned along the dress form's hip outline (pelvis and thighs), alternate petals 4 mm apart so the overlaps start clear of one another, simulated together as one merged piece so that self-contact keeps them apart, and tethered to the hip line at 1.05 times their pattern distance.
- Dressing from a hanging start: each petal begins hanging straight from its pin line rather than dropped from a flat disc. The dress form's bunching must still be understood first, since its cause is open.
- Rendering: one level of subdivision smoothing on the cloth (cheap at 480p), satin and georgette looks in render.FABRIC_LOOK, sheer sleeves that read as silk, and the bodice and fitted sleeves as materials on the body.

## Trials and acceptance

Short clips first, each a few minutes on the GPU: the gown hanging in standing (hem height, overlaps, bias lengthening, bias against straight grain); the whirl (flare against the screen, at 0, 20 and 50 g/m of hem weight); the jeté (petals parting round the split, no petal trapped); the salto (petals opening in the tucks, the torso in view in the hang); the run and carve (trailing). Acceptance: 99th-percentile strain under 10% (12% if the bias's own lengthening needs it), under 0.5% of cloth vertices inside the rendered skin in any frame, the hem within 3 cm of its design height at the motion's start, no petal caught between the legs, and 1080p stills of each unit's moment judged by eye.

## Compute

About 28,000 cloth vertices: 9,300 in the satin petals at 1.5 cm, 14,500 in georgette petals at 1.2 cm if they are used, and 4,400 in the sleeves (about 14,000 without the georgette petals). The skirts of the first renders had 5,100-8,500 and simulated about 0.2 s of motion per second of compute, so the gown should take 10-15 minutes of GPU time for the 36 s phrase. The trials take about 30-45 minutes, a few full runs 30-40 minutes, and a 480p render 35-45 minutes: about 1.5-2 hours of GPU time in all, under 2 GB of system memory, alongside CPU work.

## Open choices

1. The petals' fabric: heavy silk satin (lustrous, swings), crepe de chine (matte, lighter), or viscose jersey (a heavy fluid knit, modern).
2. Whether to interleave sheer georgette petals.
3. How high the petals should rise in the whirl: nearly level (weighted hems), lower (lighter hems or a slower whirl), or by trial.
4. Colour: indigo shading to pale with silver-grey georgette, ivory and champagne, or another palette.

## Design rounds (2 November)

The first round dressed three gowns, each with layered water sleeves (a georgette tube inside a longer chiffon one), in one cloth each on a stand-in body (body.json's trunk and round waistband with the dress form's fitted limbs, since the full form still dresses skirts short), and rendered each standing and mid-whirl: a sheer chiffon full circle over a weighted satin half circle (layered), three sheer layers graded in length (which reads as tiers), and the petal gown on the straight grain. The graded layers were preferred, with two criticisms: nothing made the pieces one outfit, and the sleeves were plain tubes held at the wrists.

The second round makes the tiered gown one outfit. The bodice matches the tiers in matte silk with a sheer yoke over the shoulders, a satin sash with two tails marks the waist and hides the skirts' bands, and the sleeves dress the arm and echo the skirt: full sheer sleeves from the armhole to the wrist, gathered at both ends, then two flounces at the wrist graded in length like the skirt's tiers, then the long chiffon extension the canopy needs. The palette is provisional (champagne, ivory and gold).

A variant with stiff organza caps after the terno's butterfly sleeves was built and held back. Pinned along the top of the shoulder, the caps were swept by the arms in the port de bras and stretched over 200%; attached to the upper arm, as a sleeve cap is, a cap that stands above the shoulder swings into the head when the arm is raised to the side, which the dressing pose, the port de bras to fifth and the canopy's cast all do. A standing cap of that kind restricts a raised arm, so for this phrase a terno-like element needs another form: a lower crest at the neck clear of the arms, soft gathered caps that fold as the arms rise, or butterfly sleeves kept for a dance with the arms low.

The third round fixed the waist and gave the sleeves their ends. Round 2's skirts and sash hung from a round band about the old capsule axis, 5-6 cm behind the rendered waist, with the layers' tops held 1.8 cm apart by the solver's contact distance: the sash stood 1.2 cm off the belly, 5 cm off the sides and 13 cm off the small of the back, and from the side the waist read as a shelf with the skirts hung from its edge. Round 3 dresses on the measured dress form (with slim shoulders), whose bunching came from the skirts' starting layout and is cured (dance README, "The cloth against the rendered body"): the skirts' bands lie on the waist outline 4, 7 and 10 mm off the skin, their tops lie together in one contact group, and the sash is fitted 16 mm off the skin.

The round-2 wrist flounces, cut from three-quarter circles whose hems were two to three times as long as their rings, read as ruffles, and the plain extension ended in a blunt open tube. In their place three sheer tiers hang at each wrist, graded in length as the skirt is (1.40, 0.75 and 0.35 m; the longest, the canopy's, in champagne georgette inside, the others in ivory chiffon), each widening by a fifth toward its hem and each hem bound in 1 cm of satin, gold in the stills to echo the sash. Of the two variants, hems cut on a slant (12% longer at the back of the arm) fall to soft points that trail in the canopy and the whirl; level hems end square and read as the ends of tubes. Tiers lying one inside another must differ in circumference by more than about 7 cm at a 1.5 cm mesh, or the solver stretches the outer ones round the inner. The solver also resists compression as stiffly as stretch, so cloth gathered into a seam cannot fold as gathers do: the first centimetre below each tier's seam stays stretched (to 43-46% at the 99th percentile in the outer tier), as the bishop sleeves' gathered ends are (9-24%). Neither shows at the stills' framing; a fabric model soft in compression would cure both and needs the garment study's checks rerun. The longest tier pools on the floor while the arms are down, and the bodice is unchanged. The stills are in dance-runs/costume_design/stills3 (sheet_round3.png), rendered on the NVIDIA GPU at 256 samples, which leaves the sheer layers smooth.
