# Electromagnetic infrastructure feasibility

The bounded calculation favors **independent Sun-tracking generation and short
local buffers first**, with **microwave transfers as a conditional supplement**.
Finite tile coils can produce useful forces, but the tested single-winding
group has a narrow benefit after torque compensation and refrigeration.
A common power/control/protection field system has not earned selection.
Existing lunar protection magnets are potentially important external torque
sources and disturbances at this group's actual altitude. No new trajectory
or protected plasma volume is accepted.

## Checkpoint and bounded work

Continue the existing research branch from
`5dc49e11e9ba1ab782c3904793ce052f30b9e4b4`. The author explicitly waives the
preliminary remote-head/concurrent-work check. The clean local checkout is at
that commit. Main and all prior failed cases remain unchanged.

The decision is whether sharing power distribution, relative control and
solar-wind protection equipment earns its mass and operating complexity.
Compare the present electric-propulsion reference, tile-distributed equipment
and a hub-assisted group. Larger collectors and protection hubs are candidates.
Retain the 23.835 TW comparison target, the failed return and the distinction
between optical interception, generation, available bus power and delivery.

The environment exposes eight CPU equivalents and 8 GiB cgroup memory. New
numerical work is capped at **1,200 s wall time, one process, one numerical
thread, 2 GiB address space and 128 MiB new raw output**. Required repository
checks have a separate five-minute allowance per pre-commit gate. Charge
initial data inspection at 10 s and record subsequent numerical runs,
including failed attempts. No coupled plasma or new fleet propagation belongs
to this pass; pause for review before a substantially larger campaign.

Read the research/domain entry points, `return_repair.md`, `closure.md`,
`packing.md`, `collection.md`, their propagation and load implementations,
the protection architecture, and the September field/structure/plasma account.
Use source-bound saved hardware-loaded states and load histories, with the
strongest return correction's tighter replay. Analyze a central 3 by 3 tile
group (zero-based indices 160–162, 179–181, 198–200); retain all 361 members
when checking power availability, geometry and external field interactions.
Also inspect the original and migrated encounter pairs to avoid treating a
central group as an encounter bound.

1. Reconstruct time-dependent loads, actual separations and directions,
   individual power shortfalls and optimistic pooled-bus limits. Screen local
   collectors and batteries against microwave/laser links and near-field
   induction. Include receiving area, pointing, losses, heat, mass, bus limits
   and the source of energy. Existing propulsion capacity is not generation.
2. Integrate finite square-loop fields and forces at representative measured
   configurations. Check convergence and action/reaction, compare distant
   dipoles only where valid, include torque, inductance, support, conductor,
   cryogenics, fault energy, and unwanted interactions. Vary separation,
   force and group size. Internal forces receive no net-fleet momentum credit.
3. Screen solar-wind pressure, kinetic scales and particle rigidity. Compare
   tile coils and a larger hub with the existing lunar protection architecture.
   Pressure balance and single-particle bending are diagnostics; they establish
   neither a magnetopause nor atmosphere/dose protection. Keep EUV and energetic
   radiation requirements separate. Include external wind/photon momentum.
4. Rank candidates using these calculations and primary research. State
   demonstrated elements, engineering extrapolations and exploratory choices.
   Quantify replacement and propellant sensitivities. Specify the smallest
   coupled follow-up if supported, without claiming a changed trajectory until
   changed optical forces, added mass and field/plasma forces are propagated.

Compact, source-bound results and reproducible code will accompany this report.
Raw analysis products go in `research/runs/solar_shield_array/electromagnetic/`.
No full-cycle, global-coverage or accepted return result is implied by a
feasibility screen.

The early checkpoint's required `make check` reports 695 passes, 55 skips and
the same 13 baseline failures in 32.81 s. Missing spectral inputs, generated
climate configuration and process lookup account for those failures; later
Makefile targets are not reached. No numerical implementation has changed.

The early checkpoint was published as `bfb1c4e46c1e3db968d415149a65856021571176`
through the connected GitHub interface after command-line authentication
failed. The local branch was synchronized to that identical published tree.

## Measured group and reference

The nine selected central tiles carry **56.4146 million kg**. Their centre
separations range from **14.71–37.95 km initially**, **5.58–75.19 km at hour
12**, and **1.95–113.35 km at the stop**. Centre distances below 10 km are
possible because squares occupy different planes. They are not surface
clearances. The whole population lies **14,888–15,215 km from lunar centre**;
the old 78,000 km held-screen distance would be an inappropriate field input
for these actual trajectories.

The group is representative of the central packing, not the highest-load
members: its existing actuator allocations pass, whereas **160 of 361 tiles
remain overloaded**. All-member accounting retains those failures and the
**34.232-million-kg** additive actuator/power-equipment requirement for the
old allocations plus 25% margin. Generation and receiving equipment cannot
by themselves increase those actuator ratings.

| Measured interval | Nine-tile energy | Nine-tile simultaneous peak | 361-tile simultaneous peak |
|---|---:|---:|---:|
| Service, 0–6 h | 0.06077 TJ | 4.989 MW | 0.217 GW |
| Departure, 6–12 h | 1.22093 TJ | 128.719 MW | 5.179 GW |
| Failed corrected return, 12–15.391 h | 0.66832 TJ | 185.622 MW | 26.939 GW |

The actual group peaks at **1.295 kN translational force per tile** and
**2.332 MN m attitude torque** during the separate stages. The population's
measured maximum correction force is **5.153 kN**. The return is almost
entirely translational expenditure; departure is attitude expenditure.
Commands, per-member masses, gravity/radiation moments, finite-Sun collection
and the stopped trajectory come from the prior source-bound products.

The previous middle conversion assumption still produces **2.605 GW** during
return against **8.832 GW** mean propulsion. Its **76.006 TJ net deficit** is
not a battery capacity or peak-power rating. Time-resolved deficits reach
**24.321 GW**. Integrating only positive shortages gives **82.253 TJ** with
an ideal shared bus and **88.060 TJ** with isolated local buses. A chronological
storage calculation allowing recharge during surplus, 95% efficiency each
way and 80% usable depth needs **106.193 TJ nominal capacity**, or **147.49
million kg at 200 Wh/kg**, before routing and thermal equipment. The more
conservative no-return-recharge screen gives 150.32 million kg and another
8.36 million kg of peak-discharge radiator. The 150/400 Wh/kg cases bracket
that conservative battery mass at 200.42/75.16 million kg; 400 Wh/kg is an
exploratory pack assumption. All still assume the unbuilt optical collectors.

The 15.391-hour optical resource and propulsion account in `return_repair.md`
remain the comparison: 22.096 TW intercepted, 3.048 TW redirected optical,
and 2.882 GW propulsion for 361 tiles. None of those optical watts is credited
directly as a supplied electrical bus. Simultaneously feathered peers cannot
remove their collective deficit by trading electricity with one another.

## Generation, links and thermal equipment

The engineering cases use 30% solar conversion, 95% PMAD, 300 W/kg gross
array system mass, 1,000 W/kg converter equipment, 2 kg/m² receiving surface,
and 5 kg/m² transmitting aperture/radiator. These are explicit extrapolations,
with component precedent and access depth in [electromagnetic_sources.json](electromagnetic_sources.json).
They are not demonstrated MW spacecraft systems. NASA's component survey
supports the technology classes, not a particular whole-vehicle specific mass.

An independent collector is sized to each tile's peak plus 25%, uses fresh
sunlight with no assumed capture of outgoing shield beams, and follows the
Sun during feathering. Its array radiates its own PV heat from two emissive
faces; separate PMAD radiators avoid counting the same PV heat twice. At
0.9 solar absorptivity and 0.85 emissivity per face the imposed equilibrium
is about 303 K. Illumination, eclipse, deployment and flexible support remain
design constraints. One-minute buffers are limited mainly by the assumed
1 kW/kg discharge capability, rather than by stored energy.

| Additional power equipment for the same nine tiles | Local generation | Microwave hub | Laser hub |
|---|---:|---:|---:|
| Generation/PMAD/radiator mass | 1.284 million kg | 2.213 million kg | 5.267 million kg |
| Link equipment including its heat rejection | — | 3.468 million kg | 6.421 million kg |
| Common 60 s buffer and on-tile conductor allowance | 0.349 million kg | 0.349 million kg | 0.349 million kg |
| **Total additional power mass** | **1.633 million kg** | **6.030 million kg** | **12.037 million kg** |
| New solar collecting area | 0.682 km² | 1.176 km² | 2.799 km² |

The local case is **2.90% of the existing group mass**. It adds no credited
propellant saving. A full-population extrapolation needs the actual higher
outer-tile loads, not nine-tile multiplication. The hub totals exclude hub
propulsion, protection coils, redundant routes and structural mounting. No
existing 20% fixed-mass or 300 W/kg propulsion allowance is silently reclaimed:
an engineered bill of materials could later identify genuine overlap. Cable
and charging losses fit within the imposed 25% source margin but still consume
energy; cable cooling, insulation and complete HV gear are additional work.

The local arrays occupy about **0.076%** of the nine physical tile areas.
If their full projected area replaces useful transmitted sunlight during
service, it is an optical loss that must be included, not a free off-band
collector. Their absorbed photon forces total about **2.79 N**, before
reflected light and thermal anisotropy. Small sustained forces matter near
a 100 m encounter boundary. The change in sail acceleration from added mass
can matter more. These effects have not been propagated.

At 100 W/kg array specific power the local generation component rises to
3.142 million kg; at an exploratory 1,000 W/kg it falls to 0.634 million kg.
This keeps uncertain engineering assumptions visible. For every 100 MW of
waste heat, a one-sided 350 K radiator needs **0.131 km² and 0.653 million kg**.
At 300/400 K the masses become 1.210/0.383 million kg. Space-view factors,
solar heating, pumping and temperature compatibility would revise them.

### A useful 100 MW link at 100 km

| Quantity | 5.8 GHz microwave | 1064 nm laser |
|---|---:|---:|
| Transmitter electrical input | 196.57 MW | 467.84 MW |
| Radiated / captured beam | 137.60 / 123.84 MW | 233.92 / 210.53 MW |
| Useful receiving bus | 100 MW | 100 MW |
| DC-to-DC efficiency | **50.87%** | **21.38%** |
| Transmitting / receiving diameter | 100 / 635 m | 1 / 370 m |
| Transmitter / receiver waste heat | 58.97 / 23.84 MW | 233.92 / 110.53 MW |
| Missed beam power | 13.76 MW | 23.39 MW |
| Additional link mass, excluding generation | 1.510 million kg | 3.032 million kg |

The microwave efficiency assumes 70% source conversion and 85% rectification;
the laser assumes 50% at each conversion. Both use 90% capture and 95% receiving
PMAD. Dickinson and Brown's primary microwave experiment established about
54% at roughly 496 W over 1.7 m. DARPA's 2025 optical experiment delivered
over 800 W across 8.6 km; its greater-than-20% optical-to-electric result was
measured at shorter distances and excludes laser wall-plug conversion.
Neither demonstrates the installed mass or scale used here.

Receiver dimensions are primarily thermal/flux choices: Gaussian peak flux
is limited to 1 kW/m² microwave and 5 kW/m² laser. The corresponding diffraction
spot radii are only 49.36 m and 0.102 m at 100 km. A deliberately broadened
laser beam therefore avoids a need for microradian pointing. Allowing pointing
error of 0.1 beam radius gives **296 and 173 microradians**, respectively.
Measured pair tracking rates reach **97 microradians/s** at the sampled dates.
These require active ephemeris tracking, receiver attitude control, beam
acquisition and fast interruption on receiver loss. They are not full pointing
error budgets. The 1 m laser aperture emits nearly 298 MW/m²: laser-module
combining, coating absorption and damage limits remain explicit unresolved
hardware constraints. A larger transmitting optic is available at added mass.

The microwave aperture is in the Fresnel regime at 100 km (Fraunhofer distance
387 km). The calculation uses a focused paraxial Gaussian diffraction bound,
not far-field Friis extrapolation. Link cases also cover 15 and 1,000 km;
microwave receiving area begins to grow with diffraction at the longer range.
At fixed aperture it eventually scales as range squared.

Seven measured configurations and five hub placements were checked against
all 361 physical squares. Normal-offset hubs at 100 km have nine clear
geometric paths initially but none at the sampled return dates; typical paths
cross two to five intervening films. **These are potential obstructions, not
proven RF or optical link failures.** The stored film has 98.1595% normal-
incidence transmission at 1064 nm. A NIST-based bare 10 µm silica microwave
screen gives about 3 parts per million reflection and 12 parts per million
first-order absorption at its reported loss-tangent upper bound. Complete
coatings, grazing incidence, supports, conductors, receivers, phase distortion
and beam heating need their own measurements. Five 99%-transmitting crossings
retain 95.1%; five 99.9% crossings retain 99.5%. The conservative opaque
geometry cannot be used to reject a sufficiently transparent link.

Some heat-rejection mass can genuinely be shared with receiving surfaces;
the table conservatively keeps extra radiators. Even deleting the microwave
receiver radiator allowance does not close its mass gap to local generation.
Transmitters, collectors and receivers experience photon recoil: the 100 MW
microwave case exerts about 0.459 N at emission and 0.413 N at absorption.
Opposite endpoint reactions and escaping photons must enter a coupled model.

A 10 km internal route also needs a material account. Two aluminium conductors
carrying 100 MW with 1% resistance loss weigh **3.046 million kg at 10 kV**,
**30,456 kg at 100 kV**, or **305 kg at 1 MV**, before insulation/HV hardware.
The study uses the 100 kV conductor sensitivity. Electrical topology and
distributed generation can shorten these routes considerably.

## Extended coils, control and shared circuits

Each modeled winding follows the physical **10 km square perimeter**. Source
segments use the analytic finite-wire Biot–Savart field; target force, torque
and mutual inductance are integrated around the entire square. Dipoles are
used only as distant scaling checks. The closest investigated witness needs
up to 512 quadrature nodes per edge; its refinement is 2.34e-6 relative.
Action/reaction and total angular momentum residuals are checked separately.

A 100 kA-turn winding has approximately **240,000 kg conductor**, **200,000 kg
cryostat allowance**, **0.446 GJ stored energy**, **12.14 kN circular-equivalent
hoop tension**, and **0.628 MW refrigeration input**. The total screen is
**446,703 kg per winding**, before deployment, square-corner bending, flexible
supports, quench containment and current-converter mass. The ideal self-support
floor is only 486 kg and must not be mistaken for a realizable square frame.
The assumed J=100 A/mm², density 6,000 kg/m³ and working specific strength
1 MJ/kg follow the existing protection scenarios. SPARC supports high-current
superconductor precedent, not a 40 km cold winding.

| Actual hour-12 pair | Centre separation | Force at 100 kA-turns on each | Torque on receiving tile |
|---|---:|---:|---:|
| 180 / 181 | 14.62 km | **1,364 N** | **0.632 MN m** |
| 180 / 179 | 28.52 km | 53.79 N | 0.0145 MN m |
| 180 / 160 | 75.15 km | 0.961 N | 0.00244 MN m |
| Original encounter, 325 / 347 | 16.15 km | 789.66 N | 4.077 MN m |
| Migrated encounter, 151 / 186 | 14.81 km | 775.44 N | 5.739 MN m |

Force direction is fixed by the moment orientations and extended geometry.
Equal/opposite force acts on the partner, with associated torques. The encounter
pairs can produce large unwanted moments even at modest forces. At the final
migrated encounter the same currents produce about **21.79 kN and 45.78 MN m**:
a fixed-current controller becomes dangerous as spacing changes.

Energizing all nine coils identically generates net internal forces up to
33.6 kN per tile at the final snapshot; the other 352 identically energized
coils add disturbances up to approximately 89 kN. This explicitly includes
other-member interactions. The sums of internal force and total torque vanish
within numerical accuracy. A pair force is not independently commandable in
a populated array. Frequency-separated AC interactions have a recent
three-unit ground demonstration on linear air tracks ([Kamat et al.](https://arxiv.org/abs/2601.05408)), but
introduce AC conductor loss, voltage, harmonic and control requirements here.

For two coaxial 10 km square coils, the equal currents and per-coil installed
mass needed for **100 N** illustrate separation scaling:

| Separation | Current per coil | Mass per coil |
|---|---:|---:|
| 10 km | 22.1 kA-turn | 0.259 million kg |
| 15 km | 38.9 kA-turn | 0.300 million kg |
| 30 km | 126.7 kA-turn | 0.511 million kg |
| 100 km | 1.302 MA-turn | 3.413 million kg |
| 300 km | 11.63 MA-turn | 34.69 million kg |

Far away, F scales as I₁I₂A₁A₂/d⁴. Equal currents scale as sqrt(F)d² at fixed
coil area; conductor mass follows current, and self-energy/support scale as
current squared. The dipole force overestimates the actual coaxial force by
about 80% at 15 km; by 100 km its error is about 1.7%. This explains why a
small-dipole model is unsuitable for the close pairs.

### Small current-allocation tests

Eight bounded least-squares starts plus the zero-current choice were tested
at three actual group states. Each tile has one parallel winding; all other
coils are off. The model includes **54 force/torque requirements and nine
current controls**, residual electric thrusters, and cooling of all nine
installed coils. The objective's force-error norm differs from the final
propulsion account; the zero-current choice prevents an inferior fit from
being reported as a benefit. These are local static searches, not global
optima or executed maneuvers.

| Time and condition | Existing propulsion | Best found propulsion + coil cooling, 100 kA cap |
|---|---:|---:|
| Hour 7.833, peak departure load | 128.719 MW | **134.374 MW**, coils off |
| Hour 9, lower departure load | 5.485 MW | 9.852 MW |
| Hour 13, peak return load | 185.622 MW | **166.884 MW** |

The return snapshot improves about **10.1%** in this restricted vacuum account;
raising the cap to 300 kA does not improve the best found result. That higher
cap substantially increases conductor and support mass. Nine 100 kA coils
add **4.020 million kg**, plus their electricity supply, and require 5.655 MW
continuous nominal refrigeration. Heat-leak/COP sensitivities span **0.063–3.142
MW per coil**, so thermal engineering can change the conclusion. The departure
result is unfavorable even before mass feedback.

The existing return controls require a nonzero resultant force. Integrating
its magnitude gives a necessary external-impulse fraction of **3.57% for all
361 tiles** or **26.20% for the isolated nine-tile group**, relative to summed
individual translation impulse. The complements are not achievable savings:
torque, current limits and geometry impose further constraints. Internal
magnetic forces conserve fleet momentum. Residual thrusters, solar photons,
solar-wind momentum, an anchored external field, or a hub's changed orbit must
account for every external exchange.

### Induction and a magnetic hub

The actual 180/181 hour-12 geometry has mutual inductance **0.489 mH**, versus
about **89.13 mH self-inductance**, giving k=0.00548. Ideal resonant transfer
of 100 MW at 100 Hz requires roughly **18.0 kA RMS ampere-turns**, **1.01 MV
single-turn equivalent voltage**, and **18.23 GVAR per coil**. Keeping cold
AC losses below 10 kW per winding requires **Q greater than 1.82 million**;
that loss alone costs another 2 MW of refrigeration per winding at COP 0.005.
The same current-density warm-copper comparison dissipates about **1.21 GW**.
More turns trade terminal current for voltage without removing these burdens.
Kurs et al.'s metre-scale resonance experiment supplies no such AC-loss or
voltage capability. Induction remains useful to examine at docking/service
distances, with purpose-built small couplers. The kilometer control winding
does not presently earn a second role as the main power link.

The 361 reference coils together store **0.161 TJ**, just **0.212%** of the
net return deficit. Storing that deficit magnetically in the same geometry
would require **2.174 MA-turns per tile**, with much greater forces, material
and coupled fault energy. A static field itself distributes no continuous
real power. Changing currents requires converters, energy recovery and
protection; total stored energy uses the full mutual-inductance matrix.

A **100 km square, 1 MA-turn hub**, located 100 km along the tile normal,
produces only **295 N** on a 100 kA tile and takes the opposite reaction.
Its screened winding equipment weighs **26.55 million kg**, before the tile
coils or power hub. This is 47% of the nine existing tiles' mass. Constant
position would require externally supplied reaction force. Hub sharing can
amortize equipment over a larger group, while longer ranges reduce coupling.
For a hub expanded with group radius while tile size/current remain fixed,
maintaining force requires hub current proportional to radius squared and
conductor mass roughly proportional to radius cubed. A larger group is not
an automatic mass saving. Nine, 25 and 361 distributed coils cost about
4.02, 11.17 and 161.26 million kg, with 36, 300 and 64,980 pair interactions.

## Solar wind, radiation and existing protection

The quiet plasma case uses 5 protons/cm³, 400 km/s, 10 eV electrons and 5 nT
IMF. It gives **1.338 nPa**, **102 km ion inertial length**, **2.38 km electron
inertial length**, **10.5 m Debye length**, and **835 km bulk-proton gyroradius
in the IMF**. Compressed and severe cases use 20/100 cm⁻³ at 800 km/s,
giving 21.4/107 nPa. The inherited 2/20/100 nPa pressure comparisons are also
retained. This is an environment envelope, not a weather prediction.

Actual extended square fields give these vacuum pressure-equivalence contours:

| Winding and external pressure | Axial distance | Equatorial distance | Absorption to perfect-reversal momentum scale |
|---|---:|---:|---:|
| 10 km tile, 100 kA; 2 nPa | 29.89 km | 24.66 km | **3.82–7.64 N** |
| Same; 20 nPa | 19.92 km | 17.16 km | 18.51–37.02 N |
| Same; 100 nPa | 14.78 km | 13.47 km | 56.98–113.96 N |
| 100 km hub, 1 MA; 2 nPa | 298.90 km | 246.60 km | **382–764 N** |

The pressure is B²/(2μ₀), without an assumed Chapman–Ferraro enhancement.
The momentum column assumes a circular obstacle of the equatorial radius and
can only be used if plasma actually develops that obstacle. Axial/equatorial
crossings do not specify the surface, cusps, wake or protected volume.
Overlapping tile fields cannot be added as independent protective disks.

Bamford et al.'s observations, theory and plasma-wind-tunnel experiments show
why the electron scale matters: magnetized electrons can create an electric
barrier that reflects ions even when their gyroradii exceed the obstacle.
Thus the 25 km tile scale is neither rejected by its 102 km ion scale nor
validated by vacuum magnetic pressure. Electron-scale currents, reconnection,
field orientation, density, current closure and the altered plasma distribution
must be resolved. No unmeasured plasma-inflation factor or wind-energy harvest
is credited. Solar-wind thrust is external momentum, generally downwind with
limited steering, and is much smaller than the tested kN corrections at quiet
conditions. A perfect 1 kN reversal at 2 nPa alone requires a **282 km radius**
interaction area; actual plasma coupling can be weaker.

The local region is very different from the four-lunar-radius atmospheric
target (**6,950 km**). A 246.6 km hub envelope covers only 0.126% of that disk
even if perfectly projected. Neither the tile nor hub calculation establishes
a protected downstream lunar atmosphere. Wake refilling and field topology
over the actual approximately 15,000 km source-to-Moon distance remain
uncomputed. EUV filtering remains the optical aperture's job.

For energetic protons, one radian of bending needs approximately **0.145 T m
at 1 MeV, 1.483 T m at 100 MeV and 5.657 T m at 1 GeV**. A specific straight
path 1 km above the 100 kA tile opening integrates to 0.077 T m signed field
and 0.217 T m absolute field; the corresponding scaled hub path gives
0.769/2.175 T m. Field reversals and trajectories matter. These path integrals
are not dose reductions, omnidirectional cutoffs or protection against neutral
photons/neutrons. The existing atmospheric column and local habitat material
shielding remain separate radiation resources; published particle-transport
work requires secondary particles and structures to be included.

### Reuse of the lunar architecture

The four anchored regional magnets have combined moment **1.5e21 A m²**.
Their construction, cryogenic experience, power conversion, servicing and
replacement model are useful inputs. Anchoring supplies a real reaction path
into the Moon. They cannot be reassigned to a moving group without changing
their atmospheric-protection function and paying for mobile structure.

At the **actual 14.9–15.2 thousand km radius**, their central-dipole orientation
envelope is about **43–91 nT**. A 100 kA tile can therefore experience up to
about **0.91 MN m torque**, while its gradient-force scale is only **0.18 N**.
This creates a plausible additional attitude-control resource and a material
disturbance. It does not supply the kN translation correction. The exact four
displaced sources, lunar orientation and plasma distortion must replace this
envelope before using it in dynamics.

A separate 50 nT uniform body-y orientation sensitivity was included in the
small current-allocation tests. It exchanges angular momentum with its external
source. Its best 100 kA return result is **169.193 MW**, versus 166.884 MW in
vacuum; the tested peak-departure case again prefers zero current before paying
cooling. This is one field orientation, not the actual lunar field trace or
an exclusion of other magnetic attitude arrangements. It establishes that
the existing protection field cannot be omitted from a subsequent coil test.

The September distributed-field account gives a worst-direction scalar
pressure radius of **5,852 km at 100 nPa**, inside the current four-radius
target. Its old three-radius protection comparison cannot automatically be
carried forward to the current target during that pressure case. Both the
field requirement and atmospheric ion loss need renewed evaluation. No
previous magnetic hardware mass is reduced by the present local study.

## Decision, shared hardware and recurring cost

1. **Develop independent collection with short buffers first.** It supplies
   feathered tiles without bulk storage of the whole return deficit and has
   the lowest screened added mass. It still needs collector placement,
   shadows, full bus/thermal design and mass/force replay. Retain the installed
   per-tile actuator limits separately.
2. **Keep microwave assistance as the leading transfer option.** It is useful
   when a differently phased collector or existing hub provides demonstrable
   simultaneous surplus, or avoids particularly expensive local collection.
   Receiver and thermal equipment can share surfaces; tracking, metrology,
   PMAD, service equipment and short buffers can serve all functions. The
   source must pay the link's losses. Laser transfer remains a compact-
   transmitter alternative whose conversion heat dominates these cases.
3. **Retain selective magnetic trim as a bounded research candidate.** Force
   magnitude is sufficient nearby. The static return benefit earns a small
   discriminating coupled test, but the departure, mass, background-field and
   cooling results do not justify a whole-array installation. A field already
   required for local plasma work could share DC windings, cryostats and power
   electronics with slow control; required moment orientation and current
   cannot be independently assigned to both functions. Three-axis windings
   or compensation hardware have not been included in the one-winding mass.
4. **Defer a new magnetic/power mega-hub and kilometer inductive distribution.**
   Their forces, mass, AC losses and reaction obligations do not beat the
   reference group. A protection hub merits a separate plasma objective if
   local particle exposure later justifies it. Regional lunar protection and
   high-energy radiation remain separate requirements.

RF/laser transmitters and receivers are separate equipment from DC coils.
Photovoltaics, antennas, beam optics and low-frequency superconducting windings
operate under different frequency, thermal and pointing constraints. Static
field energy cannot replace converters or a timed source of electricity.

No propellant saving is credited to the power options. The coil snapshot's
lower residual thrust is an instantaneous possibility, with no integrated
fuel saving accepted. Its 4.020-million-kg winding inventory would require
about **0.00637 kg/s gross replacement and 0.637 MW manufacturing power**
at a 20-year life and 100 MJ/kg. At 99.9% recovery fresh feed is a thousandth
of gross throughput; manufacturing and deployment still recur. Local power
hardware's 1.633-million-kg inventory is about 0.00259 kg/s and 0.259 MW on
the same assumptions. Five- and 100-year lifetimes scale gross rates by
four and one fifth. Battery life also depends on equivalent full cycles;
none of these service lives is established for this hardware.

Nine nominal cold coils consume **5.655 MW continuously**. Repeating only the
measured group return once per two days would average **3.868 MW propulsion**:
even eliminating that entire partial return could not pay their continuous
nominal cooling on that invented cadence. Better thermal isolation or useful
control in other stages would be necessary. At the 17.29-million-tile capacity
relaxation, one nominal cold winding per tile would consume **10.86 TW**,
about 46% of the **23.835 TW** propulsion comparison target, before control
or generation losses. That inventory is still not an operational fleet.

The target is unchanged. The earlier **44.24 TW** thought experiment already
repeats only a failed prefix at an assumed two-day cadence. Extra generation
can pay electrical demand, but it cannot remove exhausted propellant or
validate that cadence. No 50/100 TW relaxation is selected, and no new global
energy/material saving follows from this feasibility account.

## Smallest discriminating follow-up, not executed

First perform **one nine-tile collector/short-buffer replay and its comparator**
from the original epoch through the measured stop, using the saved 352-member
environment as a declared external boundary. Add the actual collector and
receiver/heat-rejection geometry, per-tile mass/inertia, source visibility,
photon momentum and bus states. Test a compact collector placement that
preserves the assigned service rays. Charge all actuator overloads and check
clearance against every prescribed outside tile. This diagnostic can reject
a local supply design; changing the outside tiles' illumination or forces
requires a subsequent full-population replay before trajectory acceptance.

For magnetic trim, the smallest useful extension keeps that same nine-member
group, exact square-loop forces/torques and individual current/voltage/rate
limits, includes the real anchored-field orientation, mutual magnetic energy,
quench/fault shutdown and continuous cooling, and carries all residual external
forces. Start from the original epoch if coils are carried throughout; do not
add their mass for free at hour 12. Compare net electrical energy, propellant,
clearance, service coverage and remaining state debt with an equal-hardware
electric comparator. The three tested departure/return states must be included.
No plasma-protection benefit should be required to make its control budget pass.

Only a surviving hardware/control candidate would justify a kinetic plasma
patch: one extended coil or the nine-coil field, quiet and compressed winds,
both IMF orientations, electron-scale boundary resolution, and measured ion
transmission and momentum deposition into the coils. A purely MHD bubble is
insufficient for the tile scale. No such plasma or new fleet integration was
executed. Review these findings before that larger work.

## Reproduction and executed checks

The calculation code is [engineering/electromagnetic.py](../../../engineering/electromagnetic.py),
with [electromagnetic_run.py](electromagnetic_run.py) and
[electromagnetic_audit.py](electromagnetic_audit.py). The compact products
are [results/electromagnetic.json](results/electromagnetic.json) and
[results/electromagnetic_audit.json](results/electromagnetic_audit.json).
They pin sources, constants, parents, original trajectories and load files;
raw group histories retain force, torque, optical power, conditional generation
and actual demand. Missing or changed prior raw files stop reproduction.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m research.studies.solar_shield_array.electromagnetic_run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m research.studies.solar_shield_array.electromagnetic_audit
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest engineering/test_electromagnetic.py research/studies/solar_shield_array/test_electromagnetic_study.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 make check
```

Two producer evaluations and four successive audits, including the initial
10-second inspection charge, used **79.68 s numerical wall time**. The latest
reproducible chain is about 49.27 s including that inspection allowance;
the audit product's `total_numerical_charge_s` describes that chain, while
the check record includes all earlier refinement work. Peak RSS was **1,028
MiB** under 2 GiB address space; raw products occupy **10.66 MiB**. All runs
were sequential and single-threaded, within the 1,200 s / 128 MiB limits.
All original failed trajectories and unsuccessful current-allocation starts
remain explicitly represented. The budget did not launch a new trajectory
or plasma simulation.

Physics checks cover the analytic square-centre field, distant dipole limit,
force as the gradient of mutual inductance, reciprocal inductance, linear
and angular momentum, electrical/thermal partitions, storage chronology,
input identities and preservation of failed-return status. Numerical refinement
is evidence about these calculations; hardware, plasma and trajectory
validation remain distinct. Final repository gate results are recorded in
`electromagnetic_checks.json`.

The final `make check` reports **705 passes, 55 skips and the same 13 baseline
failures in 28.65 s**; later Makefile targets are not reached. The earlier
targeted dynamics/engineering/study suite passed **160 tests in 10.75 s**.
After staging, the layer check covers **544 files with zero violations** and
the staged whitespace check passes. Main and all prior failure products are
unchanged.
