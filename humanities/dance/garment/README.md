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
    python -m pytest humanities/tests/test_garment.py

Use the dance venv from the worktree root. A sequence writes `faces` (F, 3) int32, `verts` (N, V, 3) float32 at the
motion's fps in the motion's world frame, `time`, `uv` (pattern coordinates in metres, for texturing), `pinned` and
`meta`. `meta` is a JSON string carrying the garment's dimensions, the fabric, gravity, air density, the solver
settings, the dressing sequence, timing and peak GPU memory, and the evidence statement (schema
`terluna.humanities.dance-garment/1`).

## What is in here

| File | Content |
|---|---|
| `cloth.py` | The solver: XPBD in small steps on the GPU (float32, one CUDA graph per substep), strain constraints in the fabric's grain, isometric bending, air, contact |
| `patterns.py` | Flat pattern pieces as triangle meshes: rectangles, tubes, the circle skirt's ring |
| `garments.py` | The full circle skirt and the water sleeves: patterns, grain, attachment, starting placement |
| `body.py` | Capsules and attachment frames from `body.json`, posed per frame and interpolated per substep; the dressing sequence; the steady-turn pose |
| `materials.py`, `sources.json` | Four silks and every outside number with its basis and standing |
| `drive.py` | Dress a garment on a motion and write its sequence |
| `validate.py` | Terminal speed, Peirce bending length, steady-turn flare and convergence, with plots |

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
contact with parallel normals (folds, a flattened sleeve) share one face to the air.

**Contact.** Contact covers the body's 23 capsules with Coulomb friction relative to the moving surface, the floor,
and self-contact between non-adjacent vertices through a sorted spatial hash rebuilt every 4 substeps; the contact
distance is 0.75 of the mean edge. Contacts are inelastic. Penetration made within a substep is removed at once;
older penetration, such as cloth caught where the mannequin's forearm capsule overlaps its hip by 3.8 cm, is removed
at 0.5 m/s, in position only.

**Garments.** The skirt is a full circle 0.85 m long on a waistband at the natural waist (2 cm above the omphalion)
that follows the pelvis and abdomen equally. The water sleeves are tubes of silk georgette 1.5 m long and 0.50 m
round, sewn to each forearm's cuff at the wrist. Each garment starts flat (the skirt as its own pattern, a horizontal
disc) or straight (the sleeves) in the first frame's pose with the arms raised 80° to the side. It settles for 6 s
with numerical damping of 2/s; the arms then come down over 1.5 s and hold for 1.5 s with damping of 0.5/s; the
motion runs with no numerical damping.
