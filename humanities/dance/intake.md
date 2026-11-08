# The dance bundles of 1 November 2026: intake and assessment

Two archives from the author's online sessions start this domain. Their bytes live on the research drive, in
`/media/projectspace/terluna-research/humanities/intake/2026-11-01/` with a `SHA256SUMS` file; the copies there match
the downloads byte for byte. They are the author's working material from outside this repository and enter it as
assessed starting points. Their code, renders and prose stay on the drive; this file records what they hold, what
they established and what this domain carries forward.

| Archive | Size | SHA-256 | Made in |
|---|---:|---|---|
| `terluna_v0.5.1_handoff.zip` | 56.5 MB | `431502b6a92d32178a7ef28dfa8973f00f991a001f147a799ae77fadb58ef518` | a Claude session (dance motion v0.3 to v0.5.1) |
| `Terluna_Dance_Research_Full_Corrected_v2_2026-11-01.zip` | 24.1 MB | `1a4583c646e17ef0ac7f9c72c4ac26b91824ecbe1d5d4d13286f36082267d11c` | a ChatGPT session (research notes, v0.1 and v0.2 simulations) |

The second archive carries the morning's first bundle inside it (`ORIGINAL_ARCHIVES/Terluna_Dance_Research_Full_Bundle_2026-11-01.zip`,
sha256 `24dbf620…3f44`, byte-identical to the copy in Downloads). A third archive in Downloads,
`Terluna_Effort_and_Originality_Evidence.zip` (sha256 `5d01f869…a1e3`), audits the author's contributions across
projects and holds copies of the dance notes and the earlier sports reports; it was listed and left aside.

## The ChatGPT research bundle

**What it holds.** Ten expanded research chapters (six movement-art families, existing practices, physics and anatomy,
martial arts and acrobatics, costume and studio, one flowing dance, a movement library and transition grammar,
worked numbers, notation, falsification tests), a bibliography of 24 entries, a claims register, a chronology of the
conversation with the author's messages verbatim, and two simulations: v0.1, a kinematic wireframe walking a 3 m
curve at 0.78 m/s, and v0.2, a 3.3 s stick-figure scissor leap with a prescribed ballistic centre of mass.

**What it established.** The author's directions in the conversation, which this domain records in the decisions
register (section "Dance and the movement arts, 1 November 2026"): review existing practice first and keep the
Terlunan language work out of the choice of movement sources; the six families as the survey's options; one flowing
dance first, in large ordinary outfits with no unusual attachments; basic moves from the chosen styles, proving what
works at low gravity, which keeps its down direction; walking is not dancing; videos and code, with web interfaces
only on request. Its framing of the problem holds: weight falls while mass and momentum stay, the floor's grip sets
what travel can do, landings and braking belong to every movement, and a convincing animation proves nothing about
force. Its arithmetic checks out (jump tables, friction and turn radii, braking distances, pendulum periods, skirt
accelerations; the Claude session re-ran it). Its ten primitives from ballroom travel, supported whirling and water
sleeves, with entry and exit states and a first transition table, are the seed of this domain's vocabulary.

**What failed.** Both simulations. The Claude session's audit (in the other archive, `history/review_2026-11-01/`)
measured them from the moving body: the v0.1 walk needs a median friction coefficient of 0.60 and exceeds 0.70 in 46%
of frames, its planted feet hover 14 mm above the floor, its centre of pressure falls 16 cm outside the supporting
foot, and its Earth comparison reuses identical kinematics; the v0.2 leap's feet are set in world space, reaching 63-74
degrees of dorsiflexion, its trailing foot detaches from the shin during recovery, and its required centre of pressure
lies up to 29 cm outside the support. The prose is defensive throughout, as the author has found Codex and ChatGPT
writing to be, and its reconstructed assistant turns are summaries, labelled as such.

## The Claude handoff (v0.5.1)

**What it holds.** A 50-degree-of-freedom MuJoCo body (1.70 m, 69.8 kg); a curved-walk planner with a linear
inverted-pendulum centre of mass; whole-body inverse kinematics (mink); a full-phrase planner (walk, three lope strides,
split leap, landing, run-out, close); a free-fall correction of every flight's angular momentum; a torque-driven
task-space QP controller that simulates the phrase in MuJoCo; MuJoCo flex cloth with per-triangle air drag for a circle
skirt and water sleeves; renders with force overlays; and the audit of the ChatGPT simulations with a cloth
feasibility spike.

**What it established.** Physics that carries forward, each finding recomputed here where this domain relies on it:
- balance, more than friction, limits lunar walking: the floor's 113 N makes the centre of pressure about six times
  more sensitive to the limbs' angular momentum, so comfortable walking slows to about 0.3 m/s at a 1.4 s step, and
  faster travel lopes;
- balance gains must follow the lunar pendulum frequency sqrt(g/z), about 1.3 per second; Earth-stiff gains saturate
  the foot;
- flights made by inverse kinematics break angular-momentum conservation (28 degrees of trunk pitch in the leap) and
  must be corrected before tracking;
- the controller's notion of contact must match the physics, and timing must follow touch-down events;
- fabric is air-dominated at 1/6 g: MuJoCo's built-in air barely acts on cloth (a falling sheet 40 times too fast),
  per-triangle pressure drag fixes it, a skirt blooms upward on any descent faster than about 0.5 m/s, and slow turns of
  10-15 rpm flare a full skirt to about 40 degrees.

**What failed.** The finish after the leap: the landing twists the trunk about 30 degrees in yaw, kills the travel, and
the run-out and close stumble; seven fixes were tried, and the full run falls at 18.45 s. The leap succeeds only with
the shipped settings and floating-point details. The cloth runs only on the kinematic walk. The figure is a capsule
mannequin and the dress a coarse mesh without self-collision: in the author's words, little better than stick figures.
The phrase is a hand-set stance table with no vocabulary or grammar beyond it, by the session's request at the time.

## What this domain carries forward and what it leaves

Carried forward, as ideas and checked numbers rather than code: the author's directions; the ten primitives, extended
into a vocabulary of 29 units with lunar-native units added; the transition-state idea, grown into a grammar with
physics-sized linking steps and choreographic devices; the physics findings above; per-triangle aerodynamic drag for
cloth; the audit's rule that every render is gated on whole-body checks; and the Claude handoff's recommendation to
replace hand tuning with optimisation over the dynamic segments, which the motion engine follows with whole-body
optimal control (Crocoddyl) and an independent replay in MuJoCo.

Left behind: the stick figures and the capsule mannequin (the renderer skins a sculpted CC0 body instead), the
hand-tuned stance table and controller, the IK-first phrase planner, the coarse one-way cloth on kinematic motion, and
the HTML viewer. The body is rebuilt from de Leva's segment data for a 1.64 m, 55 kg female dancer whose proportions
match the rendered mesh, where the handoff's 1.70 m, 69.8 kg body had capsule geometry of its own.

Both archives used an air density of 1.45 kg/m^3 for the 1.2 atm scenario; this domain reads the design run's air,
1.40 kg/m^3 at the sea-level datum and 1.17 kg/m^3 at the summit metropolis, and the shared gravity GM/R^2 = 1.6242
m/s^2 where they used 1.62.
