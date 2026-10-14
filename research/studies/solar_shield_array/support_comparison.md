# Array materials and support architectures

14 October 2026 · `research/array-industry` · baseline commit `e1196f5`

**The original roughly 200-fold mass discrepancy belongs to the tested
straight-edge structure and hardware assumptions. It is not a bound on all
known materials or all sail architectures.** This broader comparison finds
that material substitutions improve that structure by a few times. Changing
the support geometry produces much larger component savings, but none of the
tested cases establishes a complete tile within the current mass, optical and
clearance requirements.

The [result](results/support_comparison.json) contains twelve whole-tile
material substitutions, two additional load-envelope cases, thirty-two
support cases and ten spinning-membrane diagnostics. The
[material catalogue](support_materials.json) distinguishes measured properties,
assumed engineering reductions, scalar proxies and experimental fibers. The
[architecture source register](support_sources.json) distinguishes flight
demonstrations, ground experiments and proposed designs. No prior flight,
glow or tile product is overwritten.

## Material substitution in the same tile

These cases retain the original 10 km square, 100 m frame depth, 5 m edge-bow
allowance, load envelope and electrical assumptions. Each re-sizes the
momentum store and its supply against the new assembled inertia. This is the
prior **whole-tile accounting model**, with its original unresolved local
buckling, equipment attachments and flexible-control limitations.

| Frame material scenario | Modulus used, GPa | Density, kg/m³ | Annulus whole-tile mass, kg/m² |
|---|---:|---:|---:|
| Original generic composite assumption | 100 | 1,600 | 2.684 |
| 6061-T6 aluminum | 68.3 | 2,700 | 6.603 |
| AS4/8552 UD carbon/epoxy | 128 | 1,580 | 2.074 |
| M55J/#2500 UD carbon/epoxy proxy | 274 | 1,666 | 1.031 |
| S-200-FH beryllium | 290 | 1,850 | 1.081 |
| Grade F silicon carbide sensitivity | 430 | 3,200 | 1.263 |

The inherited annulus allowance is **0.013765 kg/m²**. The lightest tested
substitution remains about **75 times** that allowance. The window material
cases are also retained in the product; its bare optical film already weighs
26.921 g/m² against a complete-tile allowance of 26 g/m².

The original frame is largely stiffness-controlled. The selected materials
provide at most 2.63 times its assumed modulus per unit density. That is a
different limitation from tensile strength. The M55J case uses the composite's
274 GPa flexural modulus as a scalar frame proxy and an inferred composite
density. Its bare fiber's 540 GPa modulus is not assigned to an isotropic
compression frame. AS4 uses the measured 128 GPa longitudinal compression
modulus. [Toray M55J data](https://www.toraycma.com/wp-content/uploads/M55J-Data-Sheet.pdf),
[Hexcel 8552 data](https://www.hexcel.com/wp-content/uploads/2025/12/HexPly_8552_eu_DataSheet1.pdf).

The SiC stress cap is a declared sensitivity: a flexural strength is not a
qualified tensile allowable. Ceramic flaws, composite hoop plies, shell
buckling, joints and long-term damage can reduce installed performance. None
of these catalogued values qualifies fifteen years of service near 175 K.
[3M SiC data](https://multimedia.3m.com/mws/media/993481O/3m-silicon-carbide-materials-properties.pdf).

## Resolve the load before attributing savings to materials

The earlier transverse-load screen used the norm of angular acceleration and
angular speed. For a flat sheet, the actual normal rigid acceleration has
coefficients

\[
c_x=-\alpha_y+\omega_z\omega_x-3n^2u_zu_x,\qquad
c_y=\alpha_x+\omega_z\omega_y-3n^2u_zu_y.
\]

Here `u` is the lunar radial unit vector in the commanded frame. The square
maximum is `L (|cx|+|cy|)/2`, maximized over time. The last terms subtract
lunar differential gravity. Pure rotation about the sheet normal contributes
in-plane stress, rather than the earlier norm's transverse centrifugal term.

On the ideal circular 80° roll command, this resolves the 10 km envelope from
1.891×10⁻⁴ to **7.913×10⁻⁵ m/s²**. Holding the material fixed, the old whole-tile
ledger falls to **1.459 kg/m²**. Combining the resolved load and the M55J
proxy gives **0.566 kg/m²**, still about **41 times** the allocation. These
are separately recorded model and material changes.

The calculation retains full unrelieved solar pressure. Uniform illumination
can partly accelerate a freely flying sail instead of bending it; the residual
requires the actual mass distribution and force balance. The factor of two
on static loading is a selected screen, not a general bound on resonant
forcing. The new boom cases additionally use the exact spatial maximum of the linear rigid-body/gravity
normal-acceleration field over a finite enclosing box, evaluated at sampled
command times and including the normal tidal term. In-plane dynamic loads and full flexible equilibrium remain
open.

## Cable edges with finite supports

A tension-only net cannot hold its own aperture open in static free flight.
The new model uses four inward-curved cable edges, four compression booms and
a finite square hub. For prescribed isotropic membrane tension `T`, circular
edge radius `R` needs cable force `T R`. All booms, hub members, braces, posts
and a 20% joint allowance carry mass. The member screen includes compression,
local Euler buckling, global cantilever buckling, web shear compliance, axial
shortening, finite transverse stiffness and iterated own-mass loading. Tube
shell buckling and the assembled hub's restraint are still unvalidated.

At the **same 9.316 N/m prestress**, the original straight perimeter weighs
149.95 million kg. With the same generic material and 500 m edge retreat,
the curved-edge support weighs **5.70 million kg**. Using the M55J proxy for
compression and Kevlar yarn for tension brings that support to **2.59 million
kg**. These are support-only, in-plane calculations. The new and old member
models have different completeness, so the ratio is not a controlled proof
of a universally achievable architectural gain.

Real tensile fibers are compared in their appropriate role: M55J fiber,
S-2 glass, Kevlar 49 and Spectra S-1000, with the generic reference retained.
At the same M55J-boom geometry, all tested cable choices give support masses
between **2.57 and 2.61 million kg**. The compression structure and its minimum
member sections dominate this case; a stronger cable alone changes little.
The measured CNT and graphene fibers are included only as experimental
tension sensitivities, with their millimetre-scale test gauges and assumed
stress reductions exposed. They supply no compression-frame qualification.
[CNT experiment](https://www.nature.com/articles/s41467-026-75150-1),
[graphene-fiber experiment](https://www.nature.com/articles/s41467-024-44692-7).

### Coverage and commanded-load subtotal

The scallops remove real film area. The model charges the actual retained
film, substrate and trim mass and integrates their area moments. A central
square inside the loaded membrane contains a uniformly scaled copy of the
old 8,900×9,890 m active rectangle. Both old pitches scale by the same factor;
the inventory required by this construction scales inversely with its square. This is a
sufficient planar construction, not a new ray trace or orbit solution.

The following cases use the resolved command load, M55J compression members,
Kevlar tension members and a scaled sag/displacement allowance. They charge
the full prescribed membrane normal resultant through the four boom tips
and the supports' own normal inertial loads. They are **sheet-plus-support
subtotals**. The denominator is a reference gross square metre: mass multiplied
by the inventory required by this construction relative to the old 10 km square. It permits
comparison with the 13.765 g/m² allocation while charging lost useful area.

| Side | Edge retreat | Boom depth | Subtotal, g/reference m² | Inventory multiplier | Nominal support clearance |
|---|---:|---:|---:|---:|---:|
| 10 km | 500 m | 50 m | 66.48 | 1.208 | 148.5 m |
| 10 km | 500 m | 100 m | 32.73 | 1.208 | 98.5 m |
| 10 km | 1,000 m | 100 m | 32.10 | 1.528 | 76.0 m |
| 3 km | 150 m | 30 m | 23.15 | 13.417 | 29.5 m |
| 1 km | 50 m | 5 m | 19.48 | 120.756 | 14.8 m |

The inherited clearance gate is **150 m**, unchanged. The 10 km, 50 m-deep
case already misses it by about 1.5 m before rotors, sag, deformation reserves
and swept-component checks. Counting only the distance between chord
centrelines would incorrectly suggest a pass; finite tube radii consume
another four metres across the pair of surfaces. Deepening the structure saves mass but consumes
clearance. Smaller units require a new formation; a scaled 15 m gate would be
a new precision requirement, not a pass of the existing gate.

The 1 km case needs about **3.41 billion units** in this area-normalized
inventory. If every unit retains the old 10 kW housekeeping and 5,000 kg
avionics assumptions, those alone imply **34.1 TW and 17.0 Gt**. At the same
15-year planned life and 1% yearly fault rate, they require about **261 million
exchanges per year**. These are consequences of the specified per-unit
overheads and lifetimes, not predictions for a redesigned bus or production
system. Independent cells have independent control; dividing a mathematical
membrane into cells does not create free supports or a shared controller.

Even the smallest subtotal exceeds the complete-tile allowance, before adding
power, rotors, avionics and assembly reserve. Equipment reactions also matter:
in the 10 km, 500 m-retreat, 100 m-depth case, imposing 100 N and 1 MN·m at
each boom tip raises the subtotal from 32.73 to **48.19 g/reference m²**;
1,000 N and 10 MN·m raise it to **251.50 g/reference m²**. These are explicit
load sensitivities, with no claim that the omitted hardware has been sized.

The membrane sag calculation assumes a fixed normal boundary. A curved cable
does not supply that boundary automatically. Its normal deflection, loaded
cut pattern, attachments, substrate/optical-film compatibility and complete
spacecraft dynamics still need solving. The support-obscuration output is an
upper bound from enclosing strips; opaque hardware reduces transmitted light
and is not automatically an unfiltered UV opening.

## Spinning membranes

The diagnostic disk has 5 km radius, full film area and independently
distributed rim ballast. Its edge prestress has an explicit inertial reaction:
`Tedge = Mrim * axial_spin² / (2π)`. The calculation includes centrifugal
membrane tension, concurrent static sag and stress, and the additional normal
loads and angular momentum associated with precessing the spin axis. Actual
axial spin is the added spin plus the command frame's normal component.

| Case | Added spin, rad/s | Film-plus-ballast subtotal, g/reference m² | Substrate stress / assumed cap |
|---|---:|---:|---:|
| Fixed-normal, solar-pressure benchmark | 0.010 | 16.23 | 0.192 |
| Actual rolling command | 0.001 | 3,355.65 | 0.495 |
| Actual rolling command | 0.003 | 1,070.51 | 1.374 |
| Actual rolling command | 0.010 | 325.52 | 4.564 |

The quiet benchmark is not the commanded lunar orbit. At 0.001 rad/s, the
rolling case requires about **566 MN·m** to redirect its spin angular momentum;
the full rigid non-gravitational torque diagnostic is about 559 MN·m. Increasing
spin lowers required ballast but raises membrane stress. The tested faster
cases exceed the assumed substrate stress cap. Stronger or thicker substrates
would need a new mass/load iteration and optical assessment. Even assigning
every incident solar photon twice its momentum at the full 5 km lever arm
gives only **3.57 MN·m**. The 0.001 rad/s rolling case requires about 157 times
that generous direct-photon torque bound. It needs counter-momentum hardware,
assisted external torque or a substantially changed normal history; their
complete hardware masses have not been supplied.

These numbers already charge **1.585 times** the original count. A 5 km-radius
disk leaves holes on the old rectangular grid; that grid needs a radius of at
least 6.295 km, or correspondingly reduced pitches. These are planar coverage
conditions. The homogenized trim distribution, planar rim under rolling loads,
spin-flexible dynamics, ballast attachments, controller, power and collision
geometry remain unresolved. Counter-rotation can reduce net angular momentum,
but each surface still needs support, control and clearance.

Spin support itself has flight evidence: IKAROS deployed a roughly 14 m square
membrane centrifugally and demonstrated a reflectivity-control experiment.
Its ordinary attitude operations used thrusters. That supports the physical
principle, not the present scale or roll schedule.
[JAXA mission](https://www.isas.jaxa.jp/en/missions/spacecraft/past/ikaros.html),
[JAXA attitude experiment](https://www.jaxa.jp/press/2010/07/20100723_ikaros_e.html).

## Architecture assessment and next comparison

| Architecture | What this work establishes | Main unresolved requirement |
|---|---|---|
| Straight rigid perimeter | Material substitution and load resolution leave a large mass excess | A different span, shape or load path |
| Curved cable edge with booms/hub | Much lighter conditional support components; coverage and depth costs quantified | Flexible boundary, equipment reactions and complete control/power ledger |
| Smaller independent sails | Lower support span and explicit inventory/overhead penalties | A new safe formation and much smaller per-unit bus/service costs |
| Spinning filled membrane | Finite ballast, stress and precession trade, with planar coverage penalty | Rolling shape and controller closure, or a different normal history |
| Striped net with booms | Existing ground-test architecture has finite external reactions | Full equivalent optical/support design, not a massless net |
| Heliogyro blades | Credible solar-propulsion research direction | Its swept disk is not instantaneous filled-area shielding |

NASA's ACS3 deployed an approximately 80 m² sail; Solar Cruiser demonstrated
a ground deployment of one quadrant of its approximately 1,652 m² design.
These are relevant steps in boom deployment, far from a 10 km optical sheet.
The L'Garde striped-net test retained four booms. CABLESSail studies actuated
cables on booms to change shape and manage momentum; that deformation would
also require a fresh filtering and coverage calculation here.
[ACS3 deployment](https://www.nasa.gov/blogs/smallsatellites/2024/08/29/nasa-composite-booms-deploy-mission-sets-sail-in-space/),
[Solar Cruiser test](https://www.nasa.gov/centers-and-facilities/marshall/nasa-solar-sail-technology-passes-crucial-deployment-test/),
[striped-net experiment](https://ntrs.nasa.gov/citations/20080015847),
[CABLESSail preprint](https://arxiv.org/abs/2512.13493v2).

The next useful engineering step is a coupled flexible cable/boom tile and
coverage model, with the commanded normal history allowed to change. Compare
10 km and intermediate spans with real bus placements, photon control,
periodic momentum unloading, deformation-dependent optics and a replacement
route. A second track can test a much quieter spin-normal history against the
glow requirement. The source review does not identify a missing extraordinary
material that alone rescues the current straight-edge layout. It identifies
support geometry, normal-motion requirements and the mass/clearance trade as
the more consequential design choices. No architecture is adopted here.

## Reproduction and evidence

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.support_comparison
python -m pytest engineering/test_array_support.py engineering/test_array_spin.py research/studies/solar_shield_array/test_support_comparison.py
make check
```

The recorded producer run used one numerical thread, about 2.4 seconds and a
770 MiB process high-water memory reading. It records all producer and input
hashes and the shared constants used. The square normal-load peak changes by
0.0253% from 2,881 to 5,761 samples; cases use 2,881 samples. Numerical tests
include direct cross-product and spherical-gravity comparisons, member load
and inertia checks, retained-area quadrature refinement, actual-spin limits,
static stress/sag identities, coverage normalization and baseline preservation.
**All 81 focused tests passed.** The final `make check` passed with **2,134
Python tests passed, 16 skipped**, and **191 JavaScript tests passed**. The
layer check covered 969 files with zero violations; provenance/protection
verification and ensemble integrity passed.
[support_checks.json](support_checks.json) records the exact snapshot and commands.
Physical qualification, a full scholarly-source review, flexible dynamic
simulation and a new fleet propagation have not been performed.
