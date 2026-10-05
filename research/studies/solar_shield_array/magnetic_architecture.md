# Primary magnetic architecture comparison

## Scope checkpoint — 2026-10-05

Continue from `2208ec2f03b18bb9ce0b3226bd63bae42c074712` on the existing
research branch. The author authorized a bounded comparison of larger current
paths and upstream protection after the local electromagnetic study. The
four regional lunar stations are a reference architecture, not a selected
optimum. This comparison preserves the optical four-lunar-radius footprint,
the 23.835-TW propulsion comparison and all failed return histories.

The existing decisions exclude a superconducting planetary ring. Regional
closed circuits and separate upstream structures remain candidates. A
2,000-km-radius loop is a scaling benchmark; placing it around the Moon is
outside the retained architecture constraint.

The atmosphere's stored loss calculations distinguish the optical footprint
from the magnetic requirement. Their weaker-field cases and 1–100 kg/s loss
budgets must be carried alongside the September 1.5e21 A m² reference and a
deliberately stricter four-radius, 2/20/100-nPa vacuum-pressure screen. A
pressure contour establishes neither ion retention nor radiation dose.

## Computation budget and methods

The environment supplies eight CPU equivalents and 8 GiB memory. New numerical
work is limited to **900 seconds aggregate wall time, one numerical process,
one thread, 2 GiB address space and 64 MiB new raw products**. Repository checks
have a separate five-minute cap per gate. Begin with finite circular-loop
fields, field/current scaling, conductor and support accounts, refrigeration,
replacement and simple plasma length/time scales. Retain unsuccessful cases.

Compare the four reference stations with larger, geometrically distinct
regional circuits, reduced moments supported conditionally by stored loss
screens, and separate upstream loops. Resolve the actual extended field where
a point dipole would be inaccurate. Include self-field/bundle size, mutual
loads and quench energy. Propulsion for a held upstream source uses the existing
DE440 force machinery and stated mass/power feedback; an optical shade's
sunward position is not an equilibrium for a heavy magnet. Solar-wind momentum,
absorbed sunlight and collection hardware belong in that account.

Upstream wake width, changing wind direction, refilling and plasma-current
sustainment are bounds or sensitivities. They are not a simulated protected
atmosphere. The actual 361-tile saved geometry supplies field/torque witnesses;
no new trajectories, plasma campaign, recurrence, dose or global coverage are
accepted by this study. Findings will specify the smallest discriminating
follow-up before a substantially larger simulation campaign.
