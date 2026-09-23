# Finite forest-patch wind and structural interaction

23 September 2026. Parent checkpoint: `b5b8cc50914633b6e22587b6d3743108b10ad6f4`.

## Scope and interpretation

This checkpoint replaces the horizontally uniform canopy assumption with an
explicit group of trees, a three-dimensional conditional mean-flow calculation,
and an interacting beam/foundation network. Each tree supplies the geometry used
for its aerodynamic loading and its mechanical response. Damage changes the
geometry, and the airflow and support network are solved again.

The implemented experiment is a finite patch in a specified steady wind. Wind
climatology, turbulent gust production, nonlinear large deflection, impact
collisions, growth, reproduction and evolution remain separate work. The five
canopy shapes are imposed design comparisons. No spontaneous dome, evolutionary
optimum, storm probability or biologically viable megaforest is inferred here.

`immersion/` is outside the write scope. Existing lower-air and beam routines are
reused from `atmosphere/lower_air.py` and `biosphere/megaforest_mechanics.py`.
The existing A1-profile adapter remains available as a read-only thermodynamic
input. These runs use the standalone atmospheric scenario; they do not execute
or change the other team's renderer or atmospheric-profile generator.

## Files and reproduction

`climate/forest_patch_flow.py` owns the mean airflow. `biosphere/forest_patch.py`
owns individual-tree geometry, budgets and drag deposition.
`biosphere/forest_patch_mechanics.py` owns the interacting elastic structure and
shared soil accounting. `biosphere/forest_damage.py` owns damage state transitions.
`research/run_forest_patch.py` runs the comparisons and writes the records.
`tests/test_forest_patch.py` contains analytic, conservation and behavioral tests.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 python research/run_forest_patch.py --quick
```

The delivered numerical snapshot uses `--quick`, a 48 x 32 x 16 reference grid,
with separately recorded refinement to 96 x 64 x 32. Omitting `--quick` uses a
72 x 48 x 24 reference grid and regenerates the full study. Outputs default to
`research/runs/forest_patch/`. The output guard rejects `immersion/` and paths
resolving inside it. `--a1-profile FILE` accepts an existing A1 export and fails
on a missing or malformed input. It never silently substitutes another profile.

The full run produces individual-tree maps, flow fields, per-wave damage states,
configurations, source hashes, residuals and sensitivity results. The repository
retains compact tables and a checkpoint record; raw numerical fields accompany
the downloadable study bundle. No source or simulation has been placed in a ZIP
inside the repository.

## Shared geometry and budget controls

The reference patch contains 49 trees on a 7 x 7 grid, 105 m apart. Reference
tree height is 300 m, stem diameter H/35, crown radius 0.18 H, and leaf-area index
22 under the inherited hypothetical trait model. The reference canopy extends
about 738 m across including crowns, approximately 2.46 reference canopy heights.
This is a small patch with strong edge influence, not an established deep-forest
interior. The soil case uses 2 kPa cohesion and pore-pressure ratio 0.8.

The shape families are flat, domed, ramped, irregular, and a central-gap patch.
The dome is a rounded mound on a square footprint: its unscaled local height
factor is 0.6 + 0.4 (1 - max(|a|,|b|)^2), where a and b span -1 to 1 across the
planted positions. It is not an assumed hemispherical living structure. The ramp
rises from 0.6 to 1 of its unscaled height along one planted axis. The irregular
case uses a deterministic height variation. The gap removes nine central trees.

A scalar height multiplier then matches the **initial aboveground standing mass,
including the chosen retained water**, to the flat reference. Stem and crown
mass follow the inherited allometries. Leaf area is separately reported; it is
not held equal and is not a productivity calculation. Root construction cost,
new graft material, growth time, hydraulics and carbon maintenance are unpriced.
Consequently graft-on/graft-off comparisons bracket mechanical hypotheses and
do not constitute equal-cost engineering comparisons. Contact springs likewise
have assumed properties. The reference density is supplied by `LowerAir` at the
surface; compressible and buoyant flow are outside this checkpoint.

## Three-dimensional airflow

Let H0 and U0 be the shared reference height and imposed upstream wind speed.
Lengths and velocities are nondimensionalized by H0 and U0. The steady equations
are reached through a pseudo-time iteration of

    du/dt = -d(u)/dx - grad(pi) + nu_t laplacian(u) - b |u| u / 2
    div(u) = 0

The streamwise advecting speed in the first term is fixed at U0. This is an
**Oseen transport approximation**, not a full nonlinear Reynolds-averaged or
large-eddy turbulence calculation. The eddy viscosity is specified as
nu_t/(U0 H0) = 0.02 in the reference. Varying it is part of the sensitivity study.
There is no turbulence-production equation, stable/unstable boundary-layer
closure, surface energy balance or Coriolis term.

The solver uses staggered face velocities and cell pressures. A discrete cosine
transform solves the Neumann pressure-Poisson equation, and the matching discrete
gradient removes velocity divergence. Streamwise transport is second-order
upwind by default. The first-order option reports its leading artificial
viscosity. Convergence is measured by projected acceleration, independently of
the near-machine-zero divergence produced by the projection.

The finite channel is 12 H0 long, 8 H0 wide and 3 H0 high. Upstream and downstream
normal velocity is fixed. Ground, sides and top are impermeable and free slip.
The model therefore resolves flow around and above the patch within a confined
wind-tunnel-like domain. It does not generate a real atmospheric surface layer.
Separate doubled-length, doubled-width and doubled-top tests retain the physical
cell spacing. These checks assess enclosure sensitivity for the selected patch.

The distributed drag coefficient b = Cd a H0 comes from the same individual
crowns and stems used by the structural model. Quadrature integrates each layer's
projected area; normalized horizontal deposition preserves its Cd A even when a
narrow stem is unresolved. The force removed from the airflow is partitioned
back among the trees using those exact contributions. Thus the sum of tree loads
matches the volume drag. The spatial kernel remains a coarse-graining assumption,
and area conservation alone does not establish spatial accuracy.

The inherited coefficients are explicitly treated as porous-drag hypotheses.
A bulk whole-crown coefficient and a local volumetric drag coefficient need
calibration to be interchangeable. This issue is especially important when a
crown spans few cells. Fixed coefficients are used in this checkpoint; passing a
nonzero Vogel exponent is rejected. A future deforming-crown closure must change
the drag field and the structural exposure consistently.

## Collective mechanics

Each tree is represented by the existing tapered beam in two horizontal bending
planes, with linearized geometric stiffness from self-weight. Base translation
is fixed; base rotation now has a finite foundation stiffness. The foundation's
selected rotation at its soil-capacity threshold is 0.025 radians. This is a
constitutive assumption requiring root-soil calibration.

Crown interaction is an optional unilateral spring at a common height within the
two crowns' vertical overlap. It carries compression after the geometric gap
closes and never pulls crowns together. Initially overlapping crown envelopes
have zero preload. Springs have a prescribed 100 kN/m stiffness and 2.5 MN
breaking force. Equal and opposite forces are applied to the two beams at the
same height. Contacts can lower one tree's demand while increasing its neighbor's.
They do not add an external anchor to the stand.

The optional root-graft case adds finite rotational coupling between nearby
foundations, with a finite breaking torque. It is a reduced mechanical proxy for
a proposed connected root system. The range and stiffness do not establish that
such grafts occur naturally, that the roots reach one another, or that their
construction cost is negligible. Internal graft torques sum to zero. The model
reports stem/root-tissue demand and the torque delivered to each soil support
separately.

Where root plates overlap, a horizontal soil raster divides each shared cell
among its claimants. A tree's soil-capacity contribution is reduced by its
allocated-area fraction. Total allocated area equals the union of represented
root plates. This prevents duplicated soil capacity, while leaving the real
three-dimensional stress field and directional moment arms unresolved. Failed
root zones remain reserved during the event rather than instantly granting new
anchorage to surviving neighbors.

Stem stress, root/soil capacity, branch onset and link overload are recorded
separately. The previous small-displacement limits remain: displacement/H <= 0.1
and rotation <= 0.2 radians. The damage interpretation stops when a solution
crosses those limits or loses foundation/self-weight stability. Vertical wind
forces are retained in the records but not applied to axial stability in this
horizontal-bending checkpoint. Branch contact loads are aggregated; a resolved
branch network remains open.

## Damage feedback

At a specified upstream wind, overloaded components are identified together.
One damage wave simultaneously breaks overloaded links, removes failed trees,
and reduces damaged crowns. A crown-shedding case retains 35% of crown mass and horizontal footprint
area on its first event and removes the remainder on a subsequent event. Its
radius contracts with the square root of that fraction. Standing mass and debris
mass remain in an explicit ledger.

After a change to tree or crown geometry, airflow is solved again. After a link
change, the support network is solved again. This continues at the **same imposed
wind** until stable, empty, outside the model domain, or at an explicitly reported
iteration limit. Waves are algorithmic equilibrations and have no assigned real
duration. Simultaneous damage ordering is an assumption to test against finer
load stepping and alternative event orderings.

Removed trees contribute no further standing drag or support. Their debris mass
is retained in the accounting but debris blockage, falling-tree impacts,
root-plate excavation, regrowth and recovery are omitted. Elastic deflection
between damage events does not deform the airflow geometry. These limitations
matter before turning a cascade into a quantitative storm-loss prediction.

## Results and numerical checks

The completed final-source run includes **13 geometry/direction cases, 39
interaction comparisons, 1,884 individual-tree map records and 14 flow-grid,
domain or mixing sensitivities**. Intact comparisons use an imposed upstream
15 m/s; the damage comparison uses 22 m/s. All capacities and material traits
remain hypothetical.

### Load distribution and canopy shape

In the flat patch, mean root-capacity utilization falls from 0.526 in the
windward row to 0.334 in the last row when only aerodynamic shelter is included.
With crown contact, those row means become 0.498 and 0.345. Thus contact unloads
the leading trees while transferring demand to more sheltered neighbors.
Utilization equals calculated demand divided by the prescribed capacity; 1 is
the onset criterion. It is not a probability of failure.

Peak root utilization at 15 m/s is:

| Shape/orientation | Shelter only | Crown contact | Contact + hypothetical grafts |
|---|---:|---:|---:|
| Flat | 0.540 | 0.513 | 0.506 |
| Domed | 0.599 | 0.599 | 0.585 |
| Ramp, short edge windward | 0.617 | 0.617 | 0.603 |
| Ramp, wind across slope | 0.627 | 0.627 | 0.612 |
| Ramp, tall edge windward | 0.503 | 0.495 | 0.489 |
| Irregular, reference direction | 0.583 | 0.583 | 0.554 |
| Central gap | 0.528 | 0.519 | 0.509 |

The equal-mass dome has a maximum height near 392 m, while the flat reference is
300 m. Its redistribution of height and mass changes both exposure and capacity.
These cases provide no automatic advantage for a rounded envelope. They also
show a strong directional dependence for an asymmetric stand. The apparent
benefit of presenting the tall edge to a steady inflow applies to this model;
coherent edge gusts, acclimated edge trees and realistic atmospheric shear could
change that comparison. No evolutionary shape or preferred universal forest
architecture is established.

### Damage propagation

For the flat 49-tree patch at the same imposed 22 m/s, a single evaluation removes
12 trees. Re-solving after each damage wave removes 12, 12, 14, 7 and 4 trees,
respectively, ending with all 49 removed. The finer-grid repetition changes the
sequence to 7, 7, 7, 10, 12 and 6; it also reaches all 49. Exact onset and wave
membership are resolution-sensitive, while a cascade persists in these two runs.

Crown-contact and contact-plus-graft variants also lose all 49 in this selected
stress test. Their wave counts are five and four. Connection-mediated load
redistribution therefore does not guarantee survival under this particular
forcing. The model resolves successive exposure/support changes, rather than
holding each survivor's original sheltered loading fixed. The largest damage
sequence displacement-domain utilization is below 0.87; these runs did not cross
the stated small-displacement boundary before applying failure events.

The all-tree-loss result is conditional on the imposed wind, weak soil, small
patch, simultaneous event rules and removal of fallen obstacles. It is neither
a lunar storm prediction nor a forecast that a real megaforest would collapse.
The waves have no assigned physical duration.

### Numerical checks

All **62 available tests pass: 37 new patch tests and 25 earlier megaforest tests**.
The complete study was rerun from the final source after review. Source hashes
match the run record. The intact summary and per-tree CSV files, configurations,
and arrays in all five saved flow fields reproduce the development run exactly.
New hash metadata and schema guards make that a numerical-result comparison,
not a claim of byte identity for every JSON artifact.

Across intact geometries, maximum discrete divergence is 9.30e-16, relative
summed-tree/bulk-force discrepancy is 9.85e-16, and maximum structural equilibrium
residual is 5.16e-11. These are numerical consistency checks. The flow solver also
records its independently evaluated steady momentum ledger and acceleration
residual, so a pressure projection alone cannot certify convergence.

Refinement from 48 x 32 x 16 to 96 x 64 x 32 changes peak root utilization by
about -3.03% for the flat patch and -4.82% for the dome. Doubling the channel width
or height changes these peaks by roughly 1.0-1.14%; doubling length changes them
by about 0.03%. Halving or doubling the prescribed mixing coefficient changes
these peaks by less than 0.8%. These selected sensitivities do not constitute a
continuum limit or a complete physical uncertainty interval.

Review caught and corrected an integration bug that could free shared soil when
a tree failed. A regression test now exercises that behavior through the full
structural constructor. A separate guard binds every force extraction to the
exact drag-field hash, rejecting stale airflow after a geometry change.

The historical whole-repository checker and unrelated environmental regressions
were not run: this local worktree contains the scoped earlier checkpoint and new
files. No full-tree validation is claimed. Real A1 export integration, empirical
canopy/edge validation and full manuscript source admission remain deferred.


## Literature and evidence status

The implementation is independently written using the repository's existing
beam/air primitives and standard finite-difference/projection operations. The
following primary research informs the questions and validation targets; no
species parameters were fitted from these sources in this checkpoint.

Schelhaas et al. (2007), *Introducing tree interactions in wind damage simulation*,
Ecological Modelling 207, 197-209, DOI 10.1016/j.ecolmodel.2007.04.025. The author-posted
text distinguishes shelter, crown support, individual loading and progressive
wind damage. Abstract and relevant methods passages were inspected; complete
source admission for a manuscript remains outstanding.
https://www.researchgate.net/publication/51997058_Introducing_tree_interactions_in_wind_damage_simulation

Dupont and Brunet (2008), *Impact of forest edge shape on tree stability: a
large-eddy simulation study*, Forestry 81, 299-315, DOI 10.1093/forestry/cpn006.
The primary abstract documents edge treatments and downstream gust effects.
Only the abstract was accessible in this pass. Our mean-flow model does not
reproduce that study's resolved gust structures or establish its gust conclusions.
https://academic.oup.com/forestry/article-abstract/81/3/299/655616

Majumdar et al., *The drag length is key to quantifying tree canopy drag*,
arXiv:2411.01570. The primary abstract motivates explicit local-drag conventions
and grid-scale scrutiny. Its coefficients and calibration are not transplanted.
Abstract inspected; attempted full HTML access failed.
https://arxiv.org/abs/2411.01570

## Next discriminating work

First establish reference-quality spatial flow for selected edges using a
validated nonlinear mean-flow/turbulence treatment, a genuine upstream boundary
layer, and enough canopy heights to distinguish edge and interior behavior.
Use measured Earth canopy/edge cases as validation controls before extrapolation.
The current solver's conservation and refinement checks are numerical evidence,
not empirical calibration.

Develop geometrically nonlinear beams and foundation constitutive laws alongside
that airflow work. Add wind-dependent crown reconfiguration through the same
conserved drag/geometry interface, then compare coherent gust loading, contact
impacts and failure order. Tests of larger forest extents and multiple climate
forcing regimes can establish whether the small-patch patterns persist.

Evolutionary or engineered morphology comparisons then need light competition,
root and graft construction costs, carbon reserves, hydraulics and regeneration.
At that point canopy shape can become a model outcome rather than a prescribed
candidate. The atmospheric and biological task folders retain ownership;
immersion integration remains with the other team.
