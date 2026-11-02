# The dance programme

Work in order; each stage ends in something the author can watch or read, with its checks.

1. **Foundation (1 November).** Intake of the two bundles; the body model; the physics envelope; the vocabulary and
   grammar with three example phrases; the motion engine (whole-body optimal control in Crocoddyl, replayed in
   MuJoCo), the garment solver (cloth with per-triangle aerodynamics, an Earth control and analytic validation) and the
   renderer (the CC0 female base mesh skinned to the body, a dance stage, MP4). The run, leap and carve at the centre
   of *Under the Canopy* were realised, dressed and rendered on 1 November. The stage closes with the whole phrase,
   which comes next, ahead of the replay controller: its grammar junctions and its engine (pivots with the floor's
   turning friction, kneeling, footholds planned for the whole phrase, solving in windows) were written and tested in
   short solves on 1 November. The full solve, a review of it, its dressed render and its replay lead the next round of
   heavy computation, every render gated on the whole-body checks (friction, centre of pressure, positive contact
   force, joint ranges, torques within strength).
2. **The body's fidelity.** First, one body for physics and picture. The rendered mesh's hips and knees sit 3 cm
   lower, its shank 3 cm shorter, hands 2.4 cm longer and spine joints 2-5 cm further forward than the physics body's
   (the renderer's landmark check). Either the physics body's joint centres and segment inertias are fitted to the
   mesh, with de Leva's table as the check (recommended), or the mesh is reshaped to the physics body; the choice is
   open in [questions.md](questions.md). Then the open literature values (hip and spine ranges, floor friction, silk
   shear) and De Leva's table read on the printed page. Then strength that falls with joint speed, a compliant sole and heel pad, the garments'
   loads fed back into the body (two-way coupling) and swing-cone limits at the hips and shoulders. The MuJoCo replay
   needs a tracking controller or model-predictive control: the optimiser's own feedback gains hold a plan for about
   0.4 s against MuJoCo's soft contacts.
3. **Garments and the picture.** Light and weighted fabrics both, since different dances suit different combinations: light silk for sleeves and canopies, heavier or hem-weighted skirts for turning, studied as a matrix of fabric against unit (whirl, leap, run, canopy). The first study of fabrics and cuts (1 November; [README](README.md), "Fabrics and cuts on the Moon") built a library of twelve fabrics, a closed-form screen of each on both worlds, half-circle cuts and weighted hems in the solver, and trial runs of Cusick's drape test, skirts in steady turns and panels trailing in travel. Next for the garments: the air's added mass and its push at shallow incidence in the solver, then air that turns with a turning skirt (each changes how light silk moves, which the study could only estimate); shear that locks as a woven fabric jams, so that a heavy hem needs no tethers; measured bending for the silks whose stiffness is a modelling choice (question 5); more cuts (a double circle, capes, trains, tiers, panels at the wrists) and layered garments (a light overlay on a weighted underskirt); then attire designed for particular dances, fabric against unit. The renderer's work goes to the meshes' physics (cloth folds at chiffon's scale, its added mass and the air it carries, the body mesh bending where the physics body bends) and to seeing the dance clearly, before fine detail such as face, hair and skin that does not change how the dance reads.
4. **The vocabulary's depth.** The aerial figures the long flights open were added on 1 November (corkscrew, star wheel, flat spin, hanging salto, sleeve flight, with the *Air Study* phrase); the engine needs orientation targets in flight to realise their somersaults, cartwheels and flat spins, and the figures that need height (a jump to a ledge, a drop from it, a leap over a standing partner) wait on a future study of levels. The figures are judged for beauty from their renders. Then more units per family, floor work, and a notation export (Laban-style symbols beside the count score); the generator's transition weights learned from composed phrases rather than set by hand. Partnering with a two-body model (lifts, tosses, catches, counterbalances) is a future study, once the solo work is much further along.
5. **A dance.** A first complete piece of two to three minutes, built from motifs and their development, with music
   timed to the lunar pulse, in the garments the author chooses, rendered at final quality.
6. **Beyond the solo.** Partnering, groups and canon in a hall's real dimensions (the habitation note's halls for
   moving), the other movement-art families (contact and momentum, suspended, aerodynamic, martial), and human-factors
   tests with dancers in the analogues that exist.
