# Primary magnetic architecture comparison

**Larger regional current paths earn further development.** At the same
1.5e21 A m² net moment, four 500-km-radius regional circuits reduce the
screened hardware mass by 59 times; four 1,000-km circuits reduce it by 237
times. A weaker primary field could reduce both mass and tile torque further,
conditional on atmospheric loss and storm response. Upstream sources retain a
large recurring propulsion cost, while their apparent mass advantage depends
strongly on the protected wake. A plasma-current source remains exploratory.

The executed calculation is [magnetic_architecture_run.py](magnetic_architecture_run.py),
with [results](results/magnetic_architecture.json),
[primary sources and assumptions](magnetic_architecture_sources.json), and
[verification record](magnetic_architecture_checks.json). This is an architecture
screen; the four-station layout is reopened without selecting replacement
hardware or accepting a new magnetosphere.

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

## Requirements that actually set the size

The four-radius optical footprint blocks sunlight that would dissociate and
ionize the extended atmosphere. Magnetic protection has a different task:
reduce direct wind access and the escape of atmospheric ions. It need not
automatically reproduce either the optical shadow or the September moment.
The retained design range is 1–100 kg/s total atmospheric loss. No loss-budget
change is selected here.

The stored `protection_architecture/results/design_point.json` already contains
a moment sweep at the standard film-leakage level and four-radius optical
shadow. Reading that product gives these total losses, including its thermal
loss term, in kg/s:

| Upper-air treatment | Solar spectrum | No primary magnets | September magnets | 3.162e19 A m² |
|---|---|---:|---:|---:|
| Collisional | Quiet | 0.731 | 0.000041 | 0.000155 |
| Collisional | Solar maximum | 0.876 | 0.00268 | 0.00503 |
| LTE | Quiet | 1.045 | 0.00534 | 0.0134 |
| LTE | Solar maximum | 1.636 | 0.173 | 0.281 |
| All near-infrared energy heats upper air | Quiet | 3.657 | 0.989 | 3.594 |
| Same warm treatment | Solar maximum | 28.620 | 16.249 | 28.557 |

These are inherited screening results, with quiet **1.338 nPa wind pressure
even in the solar-maximum spectrum rows**. The ion-fate model doubles the
boundary field to represent plasma currents and estimates draining,
convection, recombination and open-line escape. It omits self-consistent
plasma loading, cusp precipitation and weak-field enhancement of escape.
The new study does not rerun or validate those assumptions. Four cool cases
support testing a moment about 47 times smaller; the two warm cases retain
large losses. At 10–100 kg/s, accepting some wind-driven loss is also a
legitimate existing comparison, with replenishment and heating consequences.

The weak field's central-dipole stand-off is 2.75 lunar radii at 1.338 nPa with
that doubled-boundary-field assumption, falling to **1.75 radii at 20 nPa and
1.34 at 100 nPa**. Those radii fall inside the cool cases' 1.91–2.37-radius
exobases. A weak steady current cannot be promoted to a storm design from the
quiet result. Ramping current for storms trades installed conductor, support,
charging time, field energy and outages against average refrigeration.

For an independent engineering comparison below, circuits also meet a specified
**minimum vacuum field on a sphere of four lunar radii** at 2/20/100 nPa.
This deliberately stronger screen assumes no boundary-field amplification.
A factor-two amplification would halve required current, substantially changing
the budget; only a plasma calculation can assign it to a particular geometry.
Pressure equivalence and a geometric sphere do not establish retention.

## Finite paths, geometry and supporting equipment

Circular filaments use the exact elliptic-integral field [1]. Tests compare it
with independent Biot–Savart quadrature, axial limits, rotations, vacuum
Maxwell identities and mutual-inductance force derivatives. No point dipole
replaces a large nearby loop.

The historical comparator retains four parallel 100-km loops centered 200 km
off the polar axis. New four-circuit layouts consist of separate spherical
caps: each path lies on the reference lunar sphere, and its local normal is
tilted so its circuit closes on that regional cap. Their transverse moments
cancel, leaving the specified net axial moment. The cap centers have five
degrees more angular distance from the pole than their own angular radius;
northern and southern pairs use perpendicular meridians. Every tested cap
remains disjoint. Two polar-cap circuits are an additional comparison.
Their currents close regionally; no circuit encloses the entire lunar globe.
Actual relief, heritage areas, settlements and sea crossings are unassigned.

The equipment account includes conductor at 100 A/mm² engineering current
density and 6,000 kg/m³ composite density, ideal tensile support at 1 MJ/kg
working specific strength, a cold envelope, refrigeration, a 350 K radiator,
PMAD and an unshadowed 300 W/kg solar source. A 10 T bundle self-field limit
sets the envelope where it exceeds the conductor packing requirement. The
warm shell allowance is 10 kg/m². Heat leak is 0.05 W/m² and COP 0.005.
These are extrapolated engineering scenarios. The SPARC experiment supplies
metre-scale 20.1 T, 40.5 kA and joint data, rather than regional qualification [2].

At fixed bundle cross-section, the circular approximation is

\[
L=\mu_0a[\ln(8a/b)-2],\qquad
U=\tfrac12 L I^2,\qquad
T=\frac{\mu_0 I^2}{4\pi}[\ln(8a/b)-1].
\]

Here I denotes ampere-turns. The hoop tension is the fixed-current energy
derivative divided by 2π. Support mass is 2πaT divided by working specific
strength. The imported model retains its slightly different internal-inductance
convention and 1-km bundle. Its 2.14e14-kg ideal conductor/support total stays
unchanged. The new common-assumption comparator weighs **2.567e14 kg**.
Its explicit 10 T envelope is about 239 m in radius, rather than 1 km.

At fixed net moment, conductor mass scales approximately as **M/a**, and
support mass and stored energy as **M² ln(8a/b)/a³**. Regional tilt and mutual
coupling alter the coefficients. Larger current paths are a real material
advantage, also motivated by primary resource studies [3]. They increase route
length and civil-engineering reach. No global optimum is claimed at the tested
upper size of 1,000 km.

## Regional comparison

All rows in this table have **1.5e21 A m² net moment**. Mass includes the
equipment listed above, before foundations, terrain works, bending/gravity
supports, fault containment and firm night supply.

| Circuits and radius | Screened mass, kg | Refrigeration | Minimum vacuum field at 4 radii | Total self plus mutual energy |
|---|---:|---:|---:|---:|
| Four historical, 100 km | 2.567e14 | 67.7 GW | 337 nT | 2.185e20 J |
| Four regional, 250 km | 2.456e13 | 27.9 GW | 350 nT | 2.131e19 J |
| Four regional, 500 km | **4.337e12** | **14.8 GW** | 375 nT | 3.585e18 J |
| Four regional, 750 km | 1.775e12 | 10.8 GW | 411 nT | 1.365e18 J |
| Four regional, 1,000 km | **1.084e12** | **9.25 GW** | 454 nT | 7.774e17 J |
| Two polar caps, 500 km | 6.815e12 | 13.6 GW | 345 nT | 5.877e18 J |
| Two polar caps, 1,000 km | 1.139e12 | 6.97 GW | 376 nT | 8.868e17 J |

A 500-km regional circuit carries 514 MA-turns, with a 10.3-m bundle radius;
a 1,000-km circuit carries 156 MA-turns, with a 3.12-m bundle. The four circuits'
aggregate route lengths are 12,566 and 25,133 km respectively. Additional
construction is substantial even though the conductor/support comparison
improves sharply.

The same **four-radius vacuum pressure requirement**, rather than the same
moment, produces this second comparison:

| Layout | Mass at 2 nPa, kg | Mass at 20 nPa, kg | Mass at 100 nPa, kg |
|---|---:|---:|---:|
| Historical four × 100 km | 1.413e13 | 1.205e14 | 5.351e14 |
| Regional four × 500 km | **2.357e11** | **1.707e12** | **7.400e12** |
| Regional four × 1,000 km | **6.124e10** | **3.349e11** | **1.285e12** |
| Two polar caps × 1,000 km | 7.320e10 | 4.634e11 | 1.910e12 |

At 100 nPa, the required net moments are 2.231e21, 2.005e21 and 1.655e21 A m²
for the first three rows. Geometry changes the field distribution, so matching
moment alone would hide part of the comparison.

Mutual forces are carried into the reaction account. The largest pair force
at reference moment falls from **2.15e12 N** in the historical layout to
**4.35e10 N** at four × 500 km and **1.28e10 N** at four × 1,000 km.
The largest individual neighbor's line load is about 12% and 32% of each
large circuit's self line load. Those distributed loads, torques and gravity
require foundations and non-tensile structure; the table's ideal support term
does not pay for them. The Moon receives the external reaction. Pair forces
and total angular momentum balance numerically.

Four circuits do not automatically provide one-circuit-out capacity. Keeping
the nominal minimum four-radius field after the worst sampled single outage
requires remaining-current factors of **1.68** at 500 km and **2.02** at
1,000 km. Building all circuits for both 100 nPa and that outage raises
screened installed mass to **1.940e13** and **4.492e12 kg**, respectively.
That retains the size advantage while making redundancy explicit. Coupled
quench transients and field collapse have not been simulated.

At the conditional **3.162e19 A m²** moment, four × 500 km cost **1.088e10 kg
and 0.541 GW refrigeration**; four × 1,000 km cost **6.056e9 kg and 0.706 GW**.
The larger route's envelope reverses the refrigeration ranking at low current.
Two 1,000-km polar caps give 4.668e9 kg and 0.402 GW. These cases earn a loss
and storm test; none receives a retained-atmosphere credit in a fleet budget.

## Power, storage and renewal

Ten-kilometer superconducting turn sections and 50 kA terminal current give
32.3 MW joint heat across four 500-km circuits at 1 nΩ per joint. Together with
the cold envelope this produces the 14.8 GW electrical refrigeration load.
Varying splice spacing from 1 to 100 km and resistance from 0.5 to 2 nΩ gives
**8.64–137 GW** with the other thermal parameters fixed. These splice lengths
and lifetimes are unqualified; metre-scale experimental joint performance
cannot establish regional production yield.

Across the tested 2.4–20 T envelope, heat-leak and COP combinations, the
500-km layout spans **3.89e12–4.56e12 kg and 3.65–186 GW refrigeration**.
Reducing working specific strength from 1 to 0.1 MJ/kg increases its support
term tenfold. This worsens absolute construction but preserves the strong
radius scaling. Current-density sensitivities of 30 and 100 A/mm² are retained
in the product. Actual conductor qualification must couple field, temperature,
strain and current density.

The source in the main mass table is continuously illuminated by assumption.
Supplying 14.8 GW through half a lunar solar cycle requires **18.85 PJ** at the
bus. A 200 Wh/kg, 80%-depth, 95%-discharge battery would add **3.44e10 kg**,
before its separate radiator and lifetime account. A regional power grid or
firm generation may be preferable.

The primary magnets could also carry a limited energy reserve. Holding the
minimum current fixed and storing the night energy above it requires an
initial current about **0.276% higher** for the reference 500-km configuration,
in a one-pass 95%-delivery calculation. At the reduced moment it needs **15.8%
extra current**. These figures expose a possible shared function; peak support,
extra refrigeration, field changes, converters and fault energy must be sized
together before assigning any storage mass saving. RF/laser power distribution
still needs its own transmitters and receivers.

With 100-year replacement, 99.9% material recovery and 100 MJ/kg manufacturing:

| Layout at reference moment | Gross replacement | Fresh material | Manufacturing power |
|---|---:|---:|---:|
| Historical four × 100 km | 81,343 kg/s | 81.3 kg/s | 8.13 TW |
| Regional four × 500 km | 1,374 kg/s | 1.37 kg/s | 137 GW |
| Regional four × 1,000 km | 343 kg/s | 0.343 kg/s | 34.3 GW |

The lower-moment 500-km case gives 3.45 kg/s gross, 0.00345 kg/s fresh and
0.345 GW manufacturing. Lifetimes of 20 and 1,000 years multiply those rates
by five and one tenth. Full loss and recharge of the field on replacement is
recorded separately; its reference 500-km average is 1.14 GW at 100 years.
Material recovery does not remove manufacturing work. These are parameter
sensitivities rather than demonstrated service lives or a closed industrial
network. The optical fleet's 23.835-TW propulsion target remains unchanged.

## Upstream sources: field size and holding cost

The upstream literature motivates protecting a planet in an artificial
magnetotail [4]. A source must produce the necessary cross-section and sustain
it at the atmosphere. The new calculations therefore keep a local field
envelope separate from downstream coverage.

For **2 nPa**, sizing the minimum finite-loop vacuum field over a sphere of
radius 6,950 km around the source gives:

| Separate loop radius | Ampere-turns | Screened mass before holding, kg |
|---|---:|---:|
| 100 km | 7.57e9 | 2.76e13 |
| 500 km | 3.01e8 | 4.12e11 |
| 1,000 km | 7.40e7 | 7.85e10 |
| 2,000 km | 1.72e7 | 1.96e10 |
| 4,000 km | 3.64e6 | 6.73e9 |

The largest loop's minimum field occurs away from its equatorial witness.
An initial equatorial-only calculation underestimated its current; the final
calculation samples the entire polar interval and retains the first pass in
ignored run products. The finite source geometry matters at this size.

The 2,000-km benchmark with the full lunar **1.5e21 A m²** instead weighs
3.56e11 kg under the new equipment assumptions. Its small mass relative to
the historical ground layout is real as a winding comparison; its position,
support and external reaction remain additional requirements.

Prescribed Sun-following loops were evaluated over **30 days of DE440**,
starting at the existing study epoch. Sixty-four points around each uniformly
loaded extended loop contribute gravitational and required kinematic
acceleration and attitude torque. Both a
Sun-facing loop and a transverse loop were tested at 15,000, 78,000, 150,000
and 1,500,000 km. This is a holding-force history, with no trajectory search.

The conservative mass closure includes the winding, cryogenics, its power
source and PMAD, 1 kW/kg thrusters, a 600 K thruster-waste radiator, seven days
of upper-bound propellant supply and 25% installed power margin. It charges
solar-wind momentum up to perfect reversal and absorbed sunlight on the
winding envelope and new collectors. Ideal rim couples supply the rigid-body
and gravity-gradient attitude demand. A uniform 5 nT background adds its
maximum moment-cross-field torque as a separate conservative sensitivity;
its direction and shielding by the plasma are unresolved. For the 2,000-km
quiet-source loop the background-torque envelope is 1.08e12 N m. Wind pressure at 2 nPa on the ideal
6,950-km obstacle transfers **0.607 MN** at perfect reversal; at 100 nPa it
would transfer 30.3 MN. The reaction is external to the fleet and generally
points toward the Moon for an upstream source.

Collectors are assumed unshadowed. Eclipse backup and a detailed distribution
of generator, radiator and propellant mass are additional design inputs. The
reported time means use trapezoidal integration over the saved uniform grid.

At **15,000 km**, every tested loop fails this scalar mass/power closure at
both 300 and 1,000 W/kg collector specific power. Peak holding/control demand
asks for at least **749 W/kg of installed generation per total carried mass**
in the tested set, before the generator's own hardware, thrusters, heat
rejection and fuel are paid. This rejects the tested continuously held
configurations; it does not rule out a future free-orbit relay choreography.

At **78,000 km**, perfect downstream preservation of the ideal source envelope
gives these conditional Sun-facing results with 300 W/kg collection:

| Loop radius | Total carried mass | Mean total bus power | Propellant |
|---|---:|---:|---:|
| 500 km | 6.57e11 kg | 19.9 TW | 30,950 kg/s |
| 2,000 km | 3.17e10 kg | **0.998 TW** | **1,551 kg/s** |
| 4,000 km | 1.11e10 kg | **0.365 TW** | **566 kg/s** |

Propulsion dominates those totals. The wind-only atmospheric losses under
consideration are much smaller than this consumed propellant in the cold
cases. Even with nearly closed construction-material recycling, exhaust
continues to leave the industrial ledger. No propellant saving is credited.
The 1,500,000-km tests follow the Moon; they are not a freely parked Sun–Earth
L1 station. A fixed L1 source must separately accommodate the Moon's lateral
motion and wake direction.

### Wake sensitivity

No validated lunar wake law is available in this study. A transparent
kinematic erosion and pointing sensitivity uses

\[
R_{\rm source}=R_{\rm target}+v_\perp D/v_{\rm sw}+D\tan\delta.
\]

The sampled transverse speeds are 0/30/60/100 km/s, and direction errors
0/5/10 degrees, at 400 km/s bulk flow. The nonzero cases are exploratory
transport margins, not a prediction of magnetized wake refilling. Magnetic
guidance and reconnection require an actual plasma solution. The quiet
electron-temperature sound-speed scale is about 31 km/s; adding ion thermal
pressure changes it, motivating a range rather than a single refill speed.

For 60 km/s and five degrees, a source at 15,000 km needs a **10,512-km**
starting radius; at 78,000 km it needs **25,474 km**. The corresponding
78,000-km budgets grow to **696 TW and 1.08 million kg/s propellant** for
the 2,000-km loop, or **111 TW and 172,505 kg/s** for the 4,000-km loop.
These deliberately retained failing cases identify the decisive unknown.

In the far field, current scales as sqrt(pressure) times obstacle-radius
cubed divided by loop-radius squared. Conductor mass scales as obstacle-radius
cubed/loop-radius; support roughly as **pressure × obstacle-radius to the
sixth / loop-radius cubed**, apart from logarithms. Modest wake-width margins
can dominate the apparent source-mass saving. Collections of small tile
cavities cannot be added as independent protective disks.

## Plasma-current alternative

Primary research proposes replacing some solid conductor with a driven plasma
current and distributed guide stations [5]. Mini-magnetosphere observations
and laboratory work support local deflection [6]; fluid-only inflation can
misrepresent momentum transfer to the inner source [7]. No empirical
inflation factor is assigned here.

A diagnostic 10,000-km-radius plasma current needs **4.77 MA** to produce
1.5e21 A m². With a 100-km minor radius its thin-path magnetic energy is
**6.71e14 J** and its hoop-tension scale **13.0 MN**. At a 100-keV electron
drift energy the relativistic speed is 0.548c, the current inventory is
1.14e25 electrons and neutralizing protons weigh only 0.019 kg. That small
particle inventory leaves guide-field forces, pressure, return currents and
stability to be supplied by the architecture.

If the magnetic energy decays in one day, replenishing it takes **7.77 GW at
perfect efficiency**, or 77.7 GW at 10% drive efficiency. A 1,000-second decay
time raises those figures to 0.671 and 6.71 TW. The lifetime is an unsolved
input; particle replenishment, field-energy decay and current-drive efficiency
are separate quantities. The carrier energy inventory is also retained in
the product. This option earns a theoretical equilibrium/current-closure
screen, with no installed-mass estimate or atmosphere-protection credit yet.

## Interaction with the actual optical tiles

The original 361-tile centres and attitudes were sampled at service start,
six hours, twelve hours and the failed return's final state. The regional
fields were evaluated directly, with four initial lunar-longitude scenarios
and a spin axis set to ecliptic north. This omits true lunar pole/libration
and plasma currents; it is an orientation sensitivity using actual tile
geometry, rather than a lunar magnetic ephemeris.

For a 100 kA-turn tile coil, maximum sampled torque is **0.425 MN m** with
the historical primary sources, **0.436 MN m** with four × 500 km and
**0.458 MN m** with four × 1,000 km, all at the same net moment. Corresponding
translation forces stay below 0.092 N. These sampled torques sit within the
earlier orientation envelope reaching 0.91 MN m. Larger regional current
paths principally save equipment; they leave a substantial background field.
At 3.162e19 A m² the sampled torque falls to **9.0–9.6 kN m**, about 47 times
smaller. The energized tile coil's moment causes this interaction.

Replacing ground sources changes the reaction source and field topology.
Both must enter any future magnetic-control replay. Internal coil interactions
exchange fleet momentum; ground-anchored fields exchange it with the Moon;
wind and photons supply external momentum. Neither a field-pressure contour
nor a modified hardware inventory improves the saved trajectories without a
new dynamics solution. Their collision and installed-power failures remain.

## Decision and smallest useful follow-up

1. **Advance regional circuits in the 500–1,000 km range.** Their material
   advantage survives matched pressure requirements and the tested material
   sensitivities. Compare four regional paths with two polar-cap paths using
   real terrain, mutual loading and failure boundaries. Keep both the reference
   moment and a current range sized to the atmospheric requirement.
2. **Test the weaker-field atmospheric option before fixing installed moment.**
   It can remove most tile torque and markedly reduce construction and cooling.
   Storm compression, warm upper air, cusp access and open-line escape determine
   whether steady weak operation, reserve current or the stronger field is needed.
3. **Keep upstream protection conditional on a sustained wake and economical
   trajectories.** Its ideal large-loop winding is lighter, but continuously
   holding it near the moving array fails the tested hardware budget. Farther
   configurations spend considerable propellant even with perfect wake survival.
4. **Retain plasma currents as a separate exploratory source.** The first gate
   is a current-carrying equilibrium with quantified guide-station reactions,
   losses and stability. Particle mass alone cannot rank it against solid coils.

Useful hardware sharing now appears in regional power, refrigeration,
servicing and a possible field-energy reserve. RF/laser transfer remains
separate equipment. The preferred primary-source candidates need not be
carried by the moving optical tiles. EUV filtering, solar-particle dose,
galactic cosmic rays and neutral radiation retain their separate requirements.

The smallest **plasma** discrimination is one cold, standard-leakage lunar
atmosphere with the four × 500-km finite field, comparing 3.162e19 with
1.5e21 A m² under quiet wind and two IMF orientations. It must include the
stored ion-production/charge-exchange source and count ions returned to the
air versus escaping, with energy deposition. A cavity-size result alone
cannot settle the loss comparison [8]. A subsequent compressed-wind case
tests the weak field's identified failure. A mesh/particle/time-step cost
estimate and short benchmark must precede the hybrid-plasma run; this study
does not authorize or launch that larger campaign.

Before that campaign, a modest structural study can reject a regional layout
using actual route relief, distributed mutual/self loads, gravity, cold
supports, joints and a quench-energy matrix. These forces are now provided.
For upstream concepts the first plasma question is width/ion transmission at
the Moon over the 37.5–195-second wind transit, including wind-direction
changes; holding work should resume only if that result supports the narrow
envelope used in the favorable cases.

## Executed verification and evidence

Seven regional geometries, two primary moments, six pressure/radius sizing
cases per geometry, one-circuit-out reserves, 30 upstream winding cases,
48 wake sensitivities, 48 monthly held-path/configuration cases with two
collector specific powers, saved-tile field witnesses, material/thermal/joint
sensitivities and a parametric plasma-current inventory were executed.
Small numerical values do not imply measured engineering precision.

The angular minimum changes by less than 0.1% between the 61×120 and
121×240 regional sphere grids. Mutual-inductance quadrature changes by less
than 8e-15 relative, force and angular-momentum balance by less than 4e-14
relative. The monthly holding peak changes by at most 0.131% between six-hour
and three-hour samples. A near-zero torque normalization in the first audit
was corrected to use the force lever-arm scale. Both that correction and the
largest upstream loop's angular-field correction are recorded in the checks.
The 13 new tests also check rigid-body attitude inertia, source hashes, constants, failed-case retention
and the absence of trajectory/magnetosphere acceptance.

Detailed check counts, timing and baseline repository failures are in
[magnetic_architecture_checks.json](magnetic_architecture_checks.json).
No coupled plasma solution, full trajectory return, new global coverage,
radiation transport, geological route or qualified support/cryogenic system
was executed. Source reading depth is recorded alongside each primary source.

## Primary sources

1. Simpson et al., [finite circular-loop fields](https://ntrs.nasa.gov/api/citations/20010038494/downloads/20010038494.pdf), 2001.
2. Hartwig et al., [SPARC toroidal field model-coil program](https://arxiv.org/abs/2308.12301), 2023.
3. DuPont and Murphy, [physical/resource constraints on a Martian magnetic shield](https://arxiv.org/pdf/2006.05546), 2021. Its atmospheric and radiation assertions are not adopted by this calculation.
4. Green et al., [upstream Mars magnetotail proposal](https://www.hou.usra.edu/meetings/V2050/pdf/8250.pdf), 2017 workshop abstract.
5. Bamford et al., [artificial magnetosphere architectures](https://arxiv.org/pdf/2111.06887), 2022 publication / 2021 preprint, especially the plasma-current discussion.
6. Bamford et al., [lunar mini-magnetospheres](https://arxiv.org/html/1207.2076v1), 2012.
7. Khazanov et al., [MHD and kinetic plasma-sail studies](https://ntrs.nasa.gov/citations/20030066145), 2003 primary research abstract.
8. Gunell et al., [magnetization and atmospheric escape](https://doi.org/10.1051/0004-6361/201832934), 2018.
