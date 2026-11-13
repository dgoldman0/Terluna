# Mining evidence and depletion sensitivities for the enclosed build

This note supports the material account for an enclosed lunar population of 6.468 million by programme year 200, together with the industry and transport that sustain it. Its purpose is to distinguish observed composition from the assumptions needed to estimate extraction effort as favorable sources are consumed. The sources support regional feedstock comparisons and conditional cost sensitivities; they do not demonstrate sufficient permitted, economically recoverable reserves for the prescribed build. [The source register](sources.json) records the inspected passages and the identity of the existing catalogue product.

Near-Earth sources here are suitable small bodies in heliocentric orbits, connected through selected flight opportunities to depots serving Earth and the Moon. Material already positioned in Earth orbit would require a separately costed capture and transfer programme. This account supplies the enclosed build and its current industrial support; the Open Moon, solar-shield array and outer Solar System programme have separate material demands.

## What the composition evidence establishes

The lunar samples show why material allocation must precede an extraction-cost estimate. Apollo 11 soil 10084 contains, by oxide mass, approximately 42% SiO₂, 15.81% FeO, 13.55% Al₂O₃, 11.96% CaO, 7.88% MgO and 7.54% TiO₂. The sampled mare material therefore offers a different processing mixture from the unusually aluminous Apollo 16 soil 67461, whose analyses include approximately 28–30% Al₂O₃, 16.5–17.6% CaO, 4.2–4.35% FeO and 0.29–0.50% TiO₂. These examples motivate separate regional feedstocks; neither defines a mine-wide grade distribution. [NASA sample 10084](https://www.lpi.usra.edu/lunar/samples/atlas/compendium/10084.pdf), [NASA sample 67461](https://www.lpi.usra.edu/lunar/samples/atlas/compendium/67461.pdf).

The inherited industrial allocation consequently has a useful physical basis: nearby lunar material supplies shielding, aggregate, glass-ceramics and ordinary structural products, while separation and fabrication determine which specified metals can be obtained. Oxide mass must first be converted to contained elemental mass and then reduced by process and quality losses. Specialist inputs need their own assessment: the reported 4.6 ppm copper in one analysis of soil 67461, for example, supplies no basis for assuming that local copper production is economical. Oxygen released during processing enters a material balance and acquires economic value only when a specified use displaces another supply.

Bennu samples contain 0.23–0.25 wt% nitrogen and 4.5–4.7 wt% carbon. This provides a chemical precedent for investigating volatile-rich asteroid material. The demonstrated laboratory extraction and hydrated mineralogy leave industrial recovery, energy demand and the variability among candidate bodies to be determined. The study's 0.24 wt% nitrogen assumption can therefore remain a reference composition, with Bennu retained as a scientific comparison rather than a designated source. [Glavin et al.](https://ntrs.nasa.gov/api/citations/20250001355/downloads/GlavinAbundantSTI.pdf).

Polar ice requires a different evidential treatment. Li and colleagues observed patchy surface ice, with signatures in approximately 3.5% of cold traps under their detection procedure; some spectral mixtures imply approximately 30 wt% ice at the optical surface. Their observations concern the upper few millimetres and do not determine excavatable thickness or continuity. The inherited source architecture therefore gains no committed polar-water inventory from these detections. A future polar extraction case would require its own survey, conservation assessment and processing account. [Li et al.](https://www.pnas.org/doi/full/10.1073/pnas.1802345115).

## Reference yields and proposed adverse sensitivities

The following sensitivities are proposed engineering tests. Their numerical changes are not measured depletion trends or predictions of conditions in year 200. Lunar yields describe usable product per unit feed; asteroid yields are explicitly the product of assumed grade and recovery. The baseline values come from the existing industrial input file identified in [sources.json](sources.json).

| Feed and product | Inherited reference | Proposed adverse case | Feed required per fixed unit of product, relative to reference |
|---|---:|---:|---:|
| Lunar bulk product | 0.90 kg/kg feed | 0.80 kg/kg feed | 1.125 times |
| Lunar metal product | 0.10 kg/kg feed | Half the reference yield: 0.05 kg/kg feed | 2 times |
| Asteroid water-equivalent product | Grade 0.10 × recovery 0.70 = 0.070 kg/kg feed | Grade 0.05 × recovery 0.50 = 0.025 kg/kg feed | 2.8 times |
| Asteroid nitrogen product | Grade 0.0024 × recovery 0.60 = 0.00144 kg/kg feed | Grade 0.0012 × recovery 0.50 = 0.00060 kg/kg feed | 2.4 times |
| Asteroid metal product | Grade 0.08 × recovery 0.60 = 0.048 kg/kg feed | Grade 0.04 × recovery 0.50 = 0.020 kg/kg feed | 2.4 times |

These multipliers apply to the relevant feed requirement. They are not multipliers on every operating expense or on the whole project's capital bill. When water, nitrogen and metal are co-products of a common source feed, the binding requirement determines throughput:

\[
F_t=\max_j\left(\frac{D_{j,t}}{g_j\eta_j q_j}\right),
\]

where \(D_{j,t}\) is the net requirement for product \(j\), \(g_j\) its grade, \(\eta_j\) its recovery and \(q_j\) the fraction meeting the required specification. The displayed sensitivities use the inherited usable-product boundary without an additional quality penalty, equivalent to \(q_j=1\). Qualification losses require a separate assumption. Unused co-products remain in controlled inventories or residues; their theoretical presence creates neither an automatic sale nor an additional feed requirement.

Annual virgin demand must include new construction, replacement losses and unrecovered consumption, while crediting usable returned material once. The inherited 95% lunar end-to-end material recovery reduces virgin replacement demand; the processing account still handles the full returning throughput. As grade falls, extraction, handling, process energy and supporting machinery must be recalculated. More distant or deeper lunar sources also change hauling and site development. Asteroid changes in departure opportunity, velocity requirement and journey duration change fleet and propellant needs, so the source-to-depot transport loop must be solved again. The NHATS documentation illustrates these date- and trajectory-dependent accessibility variables; its Earth-based human mission calculations are not direct substitutes for the model's cargo route. [NASA/JPL NHATS documentation](https://ssd-api.jpl.nasa.gov/doc/nhats.html).

## Required survey area and catalogue inventory

A transparent lunar footprint calculation uses the explicit assumptions of **1,500 kg/m³ bulk density, 5 m extraction depth and 25% accessible area**. For cumulative raw feed \(F\), the directly excavated area is

\[
A_{\rm excavated}=\frac{F}{1{,}500\times5},
\qquad
A_{\rm prospective}=\frac{F}{1{,}500\times5\times0.25}.
\]

The second quantity is the gross prospective area within which sufficient accessible material would have to be found. Per \(10^{12}\) kg of raw feed, the equations give approximately **133.3 km² excavated and 533.3 km² prospective area**. These are geometric requirements under assumed depth, density and access, without a geological survey. They exclude additional land occupied by processing, storage, transport and protection measures. A geology-informed calculation would vary depth, density and access by site.

The existing JPL catalogue calculation supplies two overlapping comparison inventories. It assigns spherical volumes and a common 1,500 kg/m³ density to bodies with reported positive diameters. The broad dark subset has albedo at most 0.1. The narrow subset also requires semimajor axis 0.8–1.3 AU, eccentricity at most 0.3 and inclination at most 15 degrees, then excludes the named science examples Bennu, Ryugu, Itokawa and Eros.

| Existing catalogue screen | Bodies | Mass at assumed density | Feed allowance if an illustrative 10% were extractable |
|---|---:|---:|---:|
| Narrow orbital-geometry and science-name screen | 17 | \(3.1188\times10^{12}\) kg | \(3.1188\times10^{11}\) kg |
| Broad dark subset with reported diameters | 498 | \(8.5112\times10^{15}\) kg | \(8.5112\times10^{14}\) kg |

These are screening upper inventories for the respective listed sets under the assumed density; the narrow set is contained within the broad set. They are neither additive nor statistical upper bounds on all near-Earth material. Missing diameters remain unknown, and the catalogue is incomplete. Low albedo and favorable orbital elements establish neither the required composition nor a specified transport route. The 10% allowance is a sensitivity assumption about extractable feed, distinct from chemical grade and recovery; it carries no mining approval or empirical reserve classification. The catalogue product and its underlying byte identity are recorded in [sources.json](sources.json).

A depletion comparison can therefore report cumulative required feed against each illustrative allowance and identify the point at which an assumed supply tier would need replacement. It cannot identify an actual exhaustion year without qualified sites and inventories. The resulting economic question is how much extra surveying, access work, processing and transport would maintain the population schedule when an earlier source becomes unavailable or less productive. Cost results belong to the coupled model; these observations and sensitivities supply its evidence boundary.
