# The dance programme

Work in order; each stage ends in something the author can watch or read, with its checks.

1. **Foundation (1 November).** Intake of the two bundles; the body model; the physics envelope; the vocabulary and
   grammar with three example phrases; then, under way, the motion engine (whole-body optimal control in Crocoddyl,
   replayed in MuJoCo), the garment solver (cloth with per-triangle aerodynamics, an Earth control and analytic
   validation) and the renderer (the CC0 female base mesh skinned to the body, a dance stage, MP4). The stage closes
   when one phrase of *Under the Canopy* is realised with dynamics, dressed and rendered, every render gated on the
   whole-body checks (friction, centre of pressure, positive contact force, joint ranges, torques within strength).
2. **The body's fidelity.** De Leva's female table read on the printed page (the transcription in sources.json
   corrected the first model's head, thorax and axes on 1 November); joint centres fitted to the rendered mesh, whose
   hips and knees sit 3 cm lower, shank 3 cm shorter, hands 2.4 cm longer and spine joints 2-5 cm further forward
   than the model's (the renderer's landmark check); strength that falls with joint speed; a compliant sole and heel
   pad; the garments' loads fed back into the body (two-way coupling); swing-cone limits at the hips and shoulders.
   The MuJoCo replay needs a tracking controller or model-predictive control: the optimiser's own feedback gains hold
   a plan for about 0.4 s against MuJoCo's soft contacts.
3. **The vocabulary's depth.** More units per family, partnering with a two-body model (lifts, tosses, catches,
   counterbalances), floor work, and a notation export (Laban-style symbols beside the count score); the generator's
   transition weights learned from composed phrases rather than set by hand.
4. **A dance.** A first complete piece of two to three minutes, built from motifs and their development, with music
   timed to the lunar pulse, in the garments the author chooses, rendered at final quality.
5. **Beyond the solo.** Groups and canon in a hall's real dimensions (the habitation note's halls for moving), the other
   movement-art families (contact and momentum, suspended, aerodynamic, martial), and human-factors tests with
   dancers in the analogues that exist.
