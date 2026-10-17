# Array materials and support architectures

14 October 2026 · `research/array-industry` · baseline commit `e1196f5`

**The roughly 200-fold excess of the straight-edge tile comes from its structure and hardware assumptions, and a change
of support geometry moves it far more than a change of material.** Material substitutions lighten that structure a few
times. Curved cable edges on booms and a hub, smaller units and spinning membranes save much more component mass, at
costs in coverage, clearance and inventory, and none of the tested cases closes a complete tile within the current mass,
optical and clearance requirements.

The [result](results/support_comparison.json) holds twelve whole-tile material substitutions, two further load
envelopes, thirty-two support cases and ten spinning-membrane diagnostics. The
[material catalogue](support_materials.json) keeps measured properties, assumed engineering reductions, scalar
proxies and experimental fibres apart, and the [architecture source register](support_sources.json) keeps flight
demonstrations, ground experiments and proposed designs apart. The flight, glow and tile products are unchanged.

## Material substitution in the same tile

These cases keep the original 10 km square, the 100 m frame depth, the 5 m allowance of edge bow, the load envelope and
the electrical assumptions, and each resizes the momentum store and its supply against the new assembled inertia. This
is the prior **whole-tile accounting model**, with its local buckling, equipment attachments and flexible control
still open.

| Frame material scenario | Modulus used, GPa | Density, kg/m³ | Annulus whole-tile mass, kg/m² |
|---|---:|---:|---:|
| Original generic composite assumption | 100 | 1,600 | 2.684 |
| 6061-T6 aluminium | 68.3 | 2,700 | 6.603 |
| AS4/8552 UD carbon/epoxy | 128 | 1,580 | 2.074 |
| M55J/#2500 UD carbon/epoxy proxy | 274 | 1,666 | 1.031 |
| S-200-FH beryllium | 290 | 1,850 | 1.081 |
| Grade F silicon carbide sensitivity | 430 | 3,200 | 1.263 |

The annulus allowance is **0.013765 kg/m²**, and the lightest substitution stays about **75 times** above it. The
product also holds the window's material cases; the window's bare optical stack alone weighs 26.921 g/m² against its
complete-tile allowance of 26 g/m².

Stiffness sets the original frame, and the selected materials give at most 2.63 times its assumed modulus per unit
density, a limit of a different kind from tensile strength. The M55J case takes the composite's 274 GPa flexural
modulus as a scalar proxy for the frame, with an inferred composite density; the bare fibre's 540 GPa belongs to
tension members. AS4 takes its measured longitudinal compression modulus of 128 GPa
([Toray M55J data](https://www.toraycma.com/wp-content/uploads/M55J-Data-Sheet.pdf),
[Hexcel 8552 data](https://www.hexcel.com/wp-content/uploads/2025/12/HexPly_8552_eu_DataSheet1.pdf)).

The silicon carbide's stress cap is a declared sensitivity drawn from a flexural strength. Ceramic flaws, composite
hoop plies, shell buckling, joints and long-term damage lower installed performance, and fifteen years of service near
175 K lies beyond every catalogued value
([3M SiC data](https://multimedia.3m.com/mws/media/993481O/3m-silicon-carbide-materials-properties.pdf)).

## The load, resolved

The tile study's transverse-load screen took the norms of the angular acceleration and speed. On a flat sheet the
normal rigid acceleration is linear in position, with coefficients

\[
c_x=-\alpha_y+\omega_z\omega_x-3n^2u_zu_x,\qquad
c_y=\alpha_x+\omega_z\omega_y-3n^2u_zu_y,
\]

where `u` is the lunar radial unit vector in the commanded frame and the last terms subtract the Moon's differential
gravity. Its largest value on the square is `L (|cx|+|cy|)/2`, maximized over time. Pure rotation about the sheet's
normal stresses it in its own plane, with no transverse term.

On the ideal circular 80° roll command this brings the 10 km envelope from 1.891×10⁻⁴ to **7.913×10⁻⁵ m/s²**. With the
material unchanged the whole-tile ledger falls to **1.459 kg/m²**, and with the M55J proxy as well to **0.566 kg/m²**,
still about **41 times** the allocation; the product records the two changes separately.

The calculation keeps the full solar pressure on the sheet, although uniform light partly accelerates a free sail
instead of bending it (the part that bends it needs the actual mass distribution and force balance), and the factor of
two on static loading is a chosen screen for a single step. The boom cases also take the exact largest normal
acceleration of the linear rigid-body and gravity field over a finite enclosing box, at the sampled command times and
with the normal tidal term. In-plane dynamic loads and the full flexible equilibrium are open.

## Cable edges with finite supports

A net of tension members alone holds no aperture open in free flight, so the model reacts four inward-curved cable
edges against four compression booms and a finite square hub. Under a prescribed isotropic membrane tension `T`, an
edge of radius `R` carries a cable force `T R`. Booms, hub members, braces, posts and a 20% allowance for joints all
carry mass, and the members are sized for compression, local Euler buckling, global cantilever buckling, web shear
compliance, axial shortening, finite transverse stiffness and their own iterated mass. Tube shell buckling and the
assembled hub's restraint are open.

At the **same 9.316 N/m of prestress** the straight perimeter weighs 149.95 million kg; the curved-edge support with
the same generic material and a 500 m edge retreat weighs **5.70 million kg**, and with the M55J proxy in compression
and Kevlar yarn in tension **2.59 million kg**. These are in-plane support masses from member models of different
completeness, so their ratio indicates the gain of the geometry, which one member model would measure.

Real tensile fibres take the tension members: M55J fibre, S-2 glass, Kevlar 49 and Spectra S-1000, with the generic
reference beside them. On the same M55J booms every cable gives a support of **2.57–2.61 million kg**: the compression
structure and its minimum member sections set this case, and a stronger cable changes little. The measured carbon
nanotube and graphene fibres enter as experimental tension sensitivities, with their millimetre-scale test gauges and
assumed stress reductions stated; they bear on tension members only
([carbon nanotube experiment](https://www.nature.com/articles/s41467-026-75150-1),
[graphene-fibre experiment](https://www.nature.com/articles/s41467-024-44692-7)).

### Coverage and commanded-load subtotal

The scallops remove film area. The model charges the film, substrate and trim mass actually retained and integrates
their area moments. A central square inside the loaded membrane holds a uniformly scaled copy of the old
8,900 × 9,890 m active rectangle, both pitches scale by the same factor, and the inventory this construction needs
scales inversely with that factor's square: a sufficient planar construction, ahead of a new ray trace and orbit
solution.

The cases below take the resolved command load, M55J compression members, Kevlar tension members and a sag and
displacement allowance scaled with the side. They carry the full prescribed normal resultant of the membrane through
the four boom tips, with the supports' own normal inertial loads, and give **sheet-plus-support subtotals** per
reference gross square metre: the mass times the inventory this construction needs, over the area of the old 10 km
squares. That compares directly with the 13.765 g/m² allocation and charges the lost useful area.

| Side | Edge retreat | Boom depth | Subtotal, g/reference m² | Inventory multiplier | Nominal support clearance |
|---|---:|---:|---:|---:|---:|
| 10 km | 500 m | 50 m | 66.48 | 1.208 | 148.5 m |
| 10 km | 500 m | 100 m | 32.73 | 1.208 | 98.5 m |
| 10 km | 1,000 m | 100 m | 32.10 | 1.528 | 76.0 m |
| 3 km | 150 m | 30 m | 23.15 | 13.417 | 29.5 m |
| 1 km | 50 m | 5 m | 19.48 | 120.756 | 14.8 m |

The clearance gate stays at **150 m**. The 10 km case 50 m deep misses it by about 1.5 m before rotors, sag, deformation
reserves and swept components; measured between chord centre lines it would pass, and the finite tube radii take
another four metres across the two surfaces. Deeper structure saves mass and spends clearance. Smaller units need a
formation of their own, whose scaled gate (15 m for the 1 km units) is a new precision requirement.

The 1 km case needs about **3.41 billion units** in this area-normalized inventory. With the tile's 10 kW of
housekeeping and 5,000 kg of avionics on every unit, those alone come to **34.1 TW and 17.0 Gt**, and the same 15-year
life and 1% yearly faults ask for about **261 million exchanges a year**: the consequence of those per-unit overheads
and lifetimes, which a bus and production system designed for small units would change. Independent cells have their
own supports and control; dividing one membrane into cells on paper gives it neither.

The smallest subtotal is above the complete-tile allowance before power, rotors, avionics and assembly reserve.
Equipment reactions weigh as well: in the 10 km case with a 500 m retreat and 100 m depth, 100 N and 1 MN·m at each
boom tip raise the subtotal from 32.73 to **48.19 g/reference m²**, and 1,000 N and 10 MN·m to
**251.50 g/reference m²**, load sensitivities for hardware still to size.

The sag calculation takes the membrane's normal boundary as fixed, which a curved cable supplies only once its normal
deflection, loaded cut pattern, attachments, the substrate's compatibility with the optical film and the whole
spacecraft's dynamics are solved. The support's obscuration is an upper bound from enclosing strips; opaque hardware
blocks the light it covers, the ultraviolet included.

## Spinning membranes

The diagnostic disc is 5 km in radius, with full film area and rim ballast distributed on its own. Its edge prestress
has an explicit inertial reaction: `Tedge = Mrim * axial_spin² / (2π)`. The calculation includes the centrifugal
membrane tension, the concurrent static sag and stress, and the normal loads and angular momentum of precessing the spin
axis; the axial spin is the added spin plus the command frame's normal component.

| Case | Added spin, rad/s | Film-plus-ballast subtotal, g/reference m² | Substrate stress / assumed cap |
|---|---:|---:|---:|
| Fixed-normal, solar-pressure benchmark | 0.010 | 16.23 | 0.192 |
| Actual rolling command | 0.001 | 3,355.65 | 0.495 |
| Actual rolling command | 0.003 | 1,070.51 | 1.374 |
| Actual rolling command | 0.010 | 325.52 | 4.564 |

The fixed-normal benchmark carries the sunlight's pressure alone; the commanded lunar orbit adds the roll. At
0.001 rad/s the rolling disc needs about **566 MN·m** to redirect its spin angular momentum, against about 559 MN·m from
the full rigid non-gravitational torque diagnostic. Faster spin lowers the ballast and raises the membrane stress, and
the faster cases tested exceed the assumed substrate stress cap; stronger or thicker substrates would take a new
iteration of mass and load and a new optical assessment. Every incident solar photon given twice its momentum at the
full 5 km lever arm makes **3.57 MN·m**, so the 0.001 rad/s rolling case needs about 157 times that bound:
counter-momentum hardware, external torque or a much quieter normal history, whose hardware masses are still to size.

These numbers already charge **1.585 times** the original count: a disc 5 km in radius leaves holes in the old
rectangular grid, which needs a radius of at least 6.295 km or correspondingly smaller pitches (planar coverage
conditions). The homogenized trim, the planar rim under rolling loads, spin-flexible dynamics, ballast attachments, the
controller, power and the collision geometry are open. Counter-rotation lowers the net angular momentum, and each
surface still needs its own support, control and clearance.

Spin support has flight evidence. IKAROS deployed a membrane about 14 m square by spinning and controlled its attitude
in an experiment with liquid-crystal reflectivity, while its ordinary attitude operations used thrusters: the physical
principle at 14 m, on a schedule far quieter than the roll
([JAXA mission](https://www.isas.jaxa.jp/en/missions/spacecraft/past/ikaros.html),
[JAXA attitude experiment](https://www.jaxa.jp/press/2010/07/20100723_ikaros_e.html)).

## Architecture assessment and next comparison

| Architecture | What this work finds | Main open requirement |
|---|---|---|
| Straight rigid perimeter | Material substitution and load resolution leave a large mass excess | A different span, shape or load path |
| Curved cable edge with booms/hub | Much lighter conditional support components; coverage and depth costs quantified | Flexible boundary, equipment reactions and complete control/power ledger |
| Smaller independent sails | Lower support span, with explicit inventory and overhead penalties | A new safe formation and much smaller per-unit bus and service costs |
| Spinning filled membrane | Finite ballast, stress and precession trade, with a planar coverage penalty | Rolling shape and controller closure, or a different normal history |
| Striped net with booms | Existing ground-test architecture with finite external reactions | A full equivalent optical and support design |
| Heliogyro blades | Credible solar-propulsion research direction | Its blades shield only the part of the swept disc they cover at each instant |

NASA's ACS3 deployed a sail of about 80 m², and Solar Cruiser deployed one quadrant of its design of about 1,652 m² on
the ground: steps in boom deployment, far from a 10 km optical sheet. The L'Garde striped-net test kept four booms.
CABLESSail actuates cables along the booms to change the sail's shape and manage momentum; here that deformation would
need its own filtering and coverage calculation
([ACS3 deployment](https://www.nasa.gov/blogs/smallsatellites/2024/08/29/nasa-composite-booms-deploy-mission-sets-sail-in-space/),
[Solar Cruiser test](https://www.nasa.gov/centers-and-facilities/marshall/nasa-solar-sail-technology-passes-crucial-deployment-test/),
[striped-net experiment](https://ntrs.nasa.gov/citations/20080015847),
[CABLESSail preprint](https://arxiv.org/abs/2512.13493v2)).

The next engineering step is a coupled flexible model of a cable-and-boom tile and its coverage, with the commanded
normal history free to change, comparing 10 km and intermediate spans with real bus placements, photon control,
periodic unloading of momentum, optics that follow the deformation and a replacement route. A second track can test a
much quieter history of the spin normal against the glow requirement. The source review finds the support geometry,
the motion of the normal and the trade of mass against clearance more consequential than any single material for the
straight-edge layout. The choice of architecture is the author's.

## Reproduction and evidence

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.support_comparison
python -m pytest engineering/test_array_support.py engineering/test_array_spin.py research/studies/solar_shield_array/test_support_comparison.py
make check
```

The producer runs on one numerical thread in about 2.3 seconds with a 91 MiB high-water mark; it was regenerated on
14 October 2026 after the source check corrected the wording of its registers (the first run, inside a larger process,
read 2.4 seconds and 770 MiB). It records every producer and input hash and the shared constants used. The square's
normal-load peak changes by 0.0253% from 2,881 to 5,761 samples, and the cases use 2,881. The tests cover direct
cross-product and spherical-gravity comparisons, member loads and inertia, refinement of the retained-area quadrature,
the limits of actual spin, static stress and sag identities, coverage normalization and the preserved baseline; **all
81 focused tests passed**. The final `make check` passed with **2,134 Python tests passed and 16 skipped**, and **191
JavaScript tests passed**; the layer check covered 969 files with no violations, and provenance and protection
verification and ensemble integrity passed. [support_checks.json](support_checks.json) records the exact snapshot and
commands of that check. Physical qualification, a full reading of the sources, a flexible dynamic simulation and a new
fleet propagation are later steps.
