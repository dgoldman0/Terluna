# Joint departure and next-service investigation

Continue the reviewed `c47c96040fc386b1d980c82a036ec91830f47390` checkpoint on
`research/solar-shield-habitat-array`. The author waived the preliminary remote
head/concurrency check. Main stays unchanged. This study treats the moving
array as an optical EUV/ionizing-radiation filter and climate dimmer; primary
magnetic protection is outside this experiment.

## Budget and candidate families

Initial allocation: **1,800 CPU seconds**, one numerical process and thread,
**2 GiB address space**, **512 MiB new raw output**. Required repository checks
have a separate five-minute allowance before each commit. Inspection/code
preparation are separate from simulation time; numerical partial attempts are
charged. Review usage before each stage. Do not silently expand the search.
The historical full tighter regression alone took about 1,960 seconds; use
an explicitly cross-checked faster force implementation and checkpointed
regressions rather than rerunning every historical producer.

Preserve all 361 physical 10 km squares, 9.89 km clear apertures, 50 g/m²
optical allocation, four layers spanning 36 km, 100 m surface floor, 50 m
replay gate, 0.001 m/s² translational cap and 0.1 degree/s attitude cap.
Retain the four-lunar-radius optical target, six-hour moving 20 km receiver
service, stored transmission products and the 13.7950037% ideal redirected
spectral allocation. Coverage is geometric evidence conditional on that
unchanged spectral model; it does not validate angular thin-film response.

Three meaningfully distinct trajectory families, at most twelve inexpensive
full-population proposals and three coupled candidates:

1. Earlier velocity shaping, beginning after the six-hour assigned service,
   with the previous arrival packing as comparator. Optimize early escape and
   the next service endpoint together; do not treat changed departure states
   as free initial conditions.
2. Changed service packing and assignment: vary depth/shear/relative velocity,
   preserve ray service, and test staggered approach times. Include the motion
   and energy of every change, and continue into the next six-hour service.
3. Low-dimensional off-service steering: common tilts and depth/row-group
   independent normals. Independent normals require generalized finite-square
   shadow, photon force, moments and surface-contact tests before any reliance.
   Document restrictions and retain failed proposals, including power failures.

Reserve a tighter independent replay for the best complete candidate or a
limiting executed obstruction. A reduced trajectory or imposed kinematic path
can propose controls and costs; it cannot pass dynamics acceptance. Publish
this checkpoint before numerical search and results afterward.

## Baselines and accounting gates

Read `cycling.md`, `packing.md`, `closure.md`, `return_repair.md`,
`collection.md`, `electromagnetic.md`, their producers and result products.
Verify point-tile cycling separately from finite-formation closure. Reproduce
the hardware-loaded passage and original return witness (325/347, hour
13.893535). Preserve the strongest previous repair (151/186, hour 15.390986,
160 installed-power overloads) as an additional comparator.

The 82.13% reduction compares optical-mass prefixes. Subsequent comparisons
use the **2.26330 billion kg** loaded formation and **51.859 TJ** twelve-hour
account. Attitude torque and energy remain chargeable. No finite-formation
return, next service, handover or global operating budget was accepted at the
starting checkpoint.

Require all-member dynamics with actual-date DE440 forces, finite Sun,
mutual shadowing, finite-square clearance, adaptive close-approach/transition
checks, beam exclusions, attitude rates, torque, individual actuator ratings
and available electrical power. Include numerical disagreement in clearance.
Changed hardware must be carried from its first installation and replayed.
Track initial service, departure and re-entry coverage separately; unassigned
or uncovered times need replacement tiles and receive no fleet credit.

The supply candidate is local Sun-tracking generation and a short buffer,
using the explicit engineering assumptions in `electromagnetic.md`. Price
collectors, PMAD, buffers, conductors and heat rejection; propagate their mass
and photon momentum before acceptance. Count first-intercept light once,
separating redirected optical power, gross generation, bus electricity,
storage state/losses, propulsion/attitude energy and heat. No external beam
capture or transmitted-band collection is free. Preserve unknown accounts.

Report local executed results first, terminal state/attitude, actual next
service and its endpoint's support for repetition. One return is distinct
from recurrence. Keep **238.35 TW** held-array comparison and **23.835 TW**
operating target. **17.29 million tiles** remains a capacity relaxation. No
failed-prefix extrapolation is an achieved fleet demand. Pause before a
substantially larger fleet campaign.

The early required `make check` passed the layer check and reported 718 Python
passes, 55 skips and 13 existing failures in 29.03 s. Missing pinned spectral
inputs, generated climate configuration and process lookup remain the same
baseline blockers. Later Makefile targets were not reached.
