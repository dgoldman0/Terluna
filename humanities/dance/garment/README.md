# Garments: cloth with real aerodynamics for dance at lunar gravity

Ordinary flowing garments, a full circle skirt and a pair of Chinese water sleeves, simulated on any motion file
(`humanities/dance/motion_io.py`) for the Moon (1.6242 m/s², air 1.40 kg/m³) and an Earth control (9.80665 m/s²,
1.225 kg/m³). Nothing else changes between the two worlds. At a sixth of Earth's gravity the air's force on cloth
counts about six times more against the cloth's weight, so the aerodynamics carry the result. Measured validation
numbers and the plain list of shortcomings are in
`/media/projectspace/terluna-research/humanities/agents/garment/report.md`.

## Running

    python -m humanities.dance.garment.drive skirt moon MOTION.npz OUT.npz [--fabric silk_chiffon] [--h 0.02]
    python -m humanities.dance.garment.drive sleeves earth MOTION.npz OUT.npz [--fabric silk_georgette]
    python -m humanities.dance.garment.validate OUTDIR [--only terminal,cantilever,spin,convergence,spin_convergence]
    python -m humanities.dance.garment.screen          # the fabric and cut study's screen: results/garments.json
    python -m humanities.dance.garment.trials OUTDIR [--only drape,satin,turn,stream] [--fast]
    python -m humanities.dance.garment.figures OUTDIR
    python -m pytest humanities/tests/test_garment.py humanities/tests/test_garment_study.py

Use the dance venv from the worktree root. A sequence writes `faces` (F, 3) int32, `verts` (N, V, 3) float32 at the
motion's fps in the motion's world frame, `time`, `uv` (pattern coordinates in metres, for texturing), `pinned` and
`meta`. `meta` is a JSON string carrying the garment's dimensions, the fabric, gravity, air density, the solver
settings, the dressing sequence, timing and peak GPU memory, and the evidence statement (schema
`terluna.humanities.dance-garment/1`).

## What is in here

| File | Content |
|---|---|
| `cloth.py` | The solver: XPBD in small steps on the GPU (float32, one CUDA graph per substep), strain constraints in the fabric's grain, isometric bending, air, contact |
| `patterns.py` | Flat pattern pieces as triangle meshes: rectangles, tubes, the circle skirt's ring, and ring sectors sewn at a seam that turns the pattern (skirts cut from part of a circle or from more than one) |
| `garments.py` | The full circle skirt, circle skirts of any cut up to a full circle with or without a weighted hem, and the water sleeves: patterns, grain, attachment, starting placement |
| `body.py` | Capsules and attachment frames from `body.json`, posed per frame and interpolated per substep; the dressing sequence; the steady-turn pose |
| `materials.py`, `sources.json` | Thirteen fabrics (nine silks, a wool, a cotton, a nylon tulle, a viscose jersey) and every outside number with its basis and standing; each fabric's bending is marked measured, inferred (from a measured drape), bracket or chosen |
| `drive.py` | Dress a garment on a motion and write its sequence |
| `validate.py` | Terminal speed, Peirce bending length, steady-turn flare and convergence, with plots |
| `screen.py` | The fabric and cut study's closed-form screen of every fabric on the Moon and on Earth: `results/garments.json` |
| `trials.py` | The study's trial runs: Cusick's drape test, circle skirts in steady turns (weighted hems tethered to the waistband), a panel trailing in travel; findings in the dance README, "Fabrics and cuts on the Moon" |
| `costume.md` | The first full costume's design brief: the petal gown, its screen against the phrase and the solver work it needs (a proposal of 2 November) |
| `figures.py` | The study's figures from the screen and the trials |

## The physics, briefly

**Integration.** The solver is extended position-based dynamics in its small-steps form: 40 substeps per 60 Hz frame
(dt = 1/2400 s), with one constraint pass per substep, or two for stretch. Constraints run in parallel (Jacobi) with
mass splitting, so a stiff constraint converges like an average and a soft one keeps its true stiffness. The solver
updates velocity from the projections' increments (v + D/dt) rather than from position differences. A position
difference in float32 carries an error of ulp/dt, about 1 mm/s near 1–10 m, which is larger than lunar gravity's
gain per substep (0.7 mm/s); computed that way, a falling lunar sheet never moved.

**Fabric.** Each triangle carries warp and weft stretch, |F d| − 1, and the shear between the threads; the grain lies
in the flat pattern, so a circle skirt hangs on the straight grain at centre front and on the bias between. Bending
is the isometric quadratic energy written per vertex with the cotangent Laplacian of the flat pattern and Voronoi
(mixed) areas: E = Σ B/(2A_i) |(M x)_i|². At a boundary vertex the pattern's own in-plane Laplacian is carried
along by the local deformation gradient, so every flat placement of the piece is free of bending, edges included.
Cylinder tests reproduce (B/2)∫κ² dA within 2% on ring meshes and 6% on narrow grids, with the deficit at
boundaries. There is no fabric damping; the air does the damping.

**Air.** Each triangle feels pressure along its normal, −½ρC_N|v_n|v_n A (C_N = 1.17), and skin friction along its
plane, −½ρC_T|v_t|v_t A (C_T = 0.01, both faces), relative to still air. Both act semi-implicitly. Layers lying in
contact with parallel normals (folds, a flattened sleeve) share one face to the air. An open weave or net takes a
share of the normal pressure (the fabric's normal_factor: 0.6 for tulle, from porous-plate drag).
Mass sewn to the cloth, such as a weighted hem, adds to its vertices' inertia and weight but not to the air's push. Optional long-range tethers (Kim, Chentanez & Müller 2012) hold chosen vertices within a set distance of a pinned anchor; the weighted-hem trials tether every vertex to the waistband at 105% of its distance in the pattern, because one stretch pass cannot carry a heavy hem's pull up light cloth. The distance is measured the short way round the sewn piece, across its seam where that is shorter; until 2 November it was found by turning the band's pattern by the whole cut, which for cuts between a half and a full circle laid the turned band over the band itself and tethered cloth to the far side of the waist, tearing a three-quarter circle at its seam from the first second (half and full circles were unaffected). Shear does not lock as a weave jams, so a heavy hem can still shear light cloth further than a woven fabric allows.

**Contact.** Contact covers the body's 23 capsules with Coulomb friction relative to the moving surface, the floor,
and self-contact between non-adjacent vertices through a sorted spatial hash rebuilt every 4 substeps; the contact
distance is 0.75 of the mean edge. Contacts are inelastic. Penetration made within a substep is removed at once;
older penetration, such as cloth caught where the mannequin's forearm capsule overlaps its hip by 3.8 cm, is removed
at 0.5 m/s, in position only. For dressed renders the driver's `--form` replaces the capsules with a dress form measured from the rendered dancer
(visualization/dance/dress_form.py): 78 tapered pieces whose radius runs linearly along each, reaching the skin the
renderer shows, and a waist outline on which a skirt's waistband is placed, its points spaced by arc length as the
pattern's inner edge is. The study's runs and checks keep the capsules. Layers gathered into one seam may be given a
contact group (`contact_skip`, a group number per vertex): vertices of one group do not repel one another, so the tops
of skirts gathered into a waistband, or of tiers at a cuff, lie a few millimetres apart as gathered layers do, while
below the group the self-contact keeps the layers apart as before.

**Several fabrics in one cloth.** A Cloth may take a list of fabrics and `tri_material`, one index per triangle (cloth.merged_tri_material builds it for merged pieces): masses, air forces and stretch compliances are per triangle, bending and friction per vertex (the area-weighted mean where fabrics meet). With one fabric it computes exactly as before. Layers so simulated collide through self-contact, at one contact distance apart (0.75 of the mean edge), and slide on one another without friction; a light layer in a heavy layer's lee takes half its air force while stacked (tests: humanities/tests/test_garment_layers.py). garments.combine merges garments built on one body into one cloth, turning each pattern so its own grain is kept and carrying sewn seams over; garments.petals hangs open gores (patterns.gore) from an elliptical line at the high hip, cut on the bias or the straight grain, with tethers to their top edges. For a whole outfit, garments.arm_sleeves hangs full sleeves from the armhole to the wrist, their fullness gathered into both rings (a bishop sleeve); cuff_flounces sets small circle-cut flounces round the wrist, which graded in length make tiers; sleeve_tier hangs one tier of a tiered sleeve end at each wrist, a tube gathered into a ring at the wrist and widening by a stated fraction toward its hem (a cone cut from a narrow ring sector), its hem level or cut on a slant (patterns.sector takes a hem function), with the hem's last strip of triangles listed for a bound edge in another fabric; tiers that lie one inside another need each to be rounder than the one inside it by more than 2π times the self-contact distance (about 7 cm at a 1.5 cm mesh), or the solver stretches the outer ones round the inner; butterfly_shoulders gathers crisp panels into an arch over the shoulder, after the Philippine terno's standing sleeves, starting each folded as its gathers hold it (laid flat and squeezed, the solver fans a gathered panel out in its own plane and overstretches it), the folds opening as the panel rises; and sash holds a band round the waist outside the skirts' bands, on a dress form's waist outline when there is one, with two tails free at the back. The design runs that use them, with their stills, are in dance-runs/costume_design on the drive (costume.md).

**Garments.** The skirt is a full circle 0.85 m long on a waistband at the natural waist (2 cm above the omphalion)
that follows the pelvis and abdomen equally. Skirts of other cuts (a half circle, three quarters) are ring sectors sewn
at centre back with the straight grain at centre front; each starts as the cone it makes when laid smooth, and its
bending reproduces that cone's analytic energy within 1-2%. On a dress form's waist outline, which is not a circle, no
cone fits the pattern exactly: the skirt starts out along the outline's normal and turns within 6 cm into its pattern's
own cone about the outline's centre, and any starting point inside the body is moved out of it (garments.clear_of_body).
Laid out instead from the abdomen capsule's axis, the start carried a median strain of 17-20% and the dressed skirt
bunched (dance README, "The cloth against the rendered body"). The water sleeves are tubes of silk georgette 1.5 m long and 0.50 m
round, sewn to each forearm's cuff at the wrist. Each garment starts flat (the skirt as its own pattern, a horizontal
disc) or straight (the sleeves) in the first frame's pose with the arms raised 80° to the side. It settles for 6 s
with numerical damping of 2/s; the arms then come down over 1.5 s and hold for 1.5 s with damping of 0.5/s; the
motion runs with no numerical damping.
