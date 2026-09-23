# Interacting lunar forest patches

## Checkpoint and scope

This working checkpoint adds explicit tree positions, spatial mean flow, coupled
structural support and quasi-static damage updates to the earlier tree screen.
Its starting repository commit is `b5b8cc50914633b6e22587b6d3743108b10ad6f4`.
All new work belongs to `climate/`, `biosphere/`, `research/` and `tests/`.
`immersion/`, manuscript seeds and previous numerical implementations are unchanged.

This is an uncalibrated neutral-flow/mechanical requirements model. It supplies
conditional comparisons and identifies the next missing physics. No lunar wind
climatology, evolved canopy shape, biological feasibility or safety certification
is inferred. The previous individual-tree wind thresholds remain their original,
separate experiments.

## Components

- `climate/forest_flow.py`: three-dimensional projected Oseen flow with quadratic
  porous drag and explicitly prescribed neutral momentum mixing.
- `biosphere/forest_patch.py`: positioned individuals, six candidate layouts,
  matched aboveground mass proxies and conservative tree/cell drag quadrature.
- `biosphere/forest_mechanics.py`: tapered beams with self-weight geometric
  stiffness, compliant roots, compression-only crown contact, optional bonded
  crown/root connections, strength thresholds and shared soil-area accounting.
- `research/run_forest_patch.py`: layout comparisons, reference-wind envelopes,
  spatial/numerical sensitivities, a longer patch, damage sequences, an edge-root
  adaptation and a solitary-tree control.
- `tests/test_forest_patch.py`: focused analytic, conservation, interaction,
  domain, repeatability and scope checks.

## Spatial flow and load exchange

Let v = u/Uref, with prescribed ambient unit vector e. The steady equation is

    (e . grad) v = -grad(pi) + ell laplacian(v)
                  - (a_eff/2) |v| v + lambda(x,y) (e-v),
    div(v) = 0.

Here ell = nu_t/Uref is a prescribed mixing length in metres. The default is
12 m, corresponding to nu_t = 240 m2/s at Uref = 20 m/s. This coefficient is a
scenario input, with 6- and 24-m comparisons. Oseen advection uses the prescribed
ambient vector; nonlinear resolved self-advection and turbulent transport are
not solved. Consequently there are no computed intermittent gusts, canopy sweeps,
unsteady separation statistics, storm frequencies or gust-return periods.

A Fourier projection enforces incompressibility. Reflection gives impermeable
free-slip ground and lid. Horizontal periodicity is accompanied by a lateral
relaxation reservoir. Default dimensions are 2400 x 1920 x 1200 m with 40 x 32 x
20 physical cells; reflection doubles the internal vertical grid. Neither the
ground shear layer nor a genuine open atmospheric inflow boundary is represented.
The stated reference wind is the reservoir's target ambient wind, rather than a
claim that every upstream location attains exactly that speed.

The sparse drag map comes from each tree's tapered stem and vertically distributed
crown. The same Cd-weighted frontal areas extract momentum from the fluid and
supply loads to the individual structures. Foliage LAI enters inherited mass
accounting; it is not independently multiplied into this aerodynamic area. The
mass density supplied by `atmosphere/lower_air.py` enters the drag, while the
continuity/transport approximation remains a shallow, constant-reference-density
one. Atmospheric stratification and buoyancy production are absent.

Crown and stem drag are deposited into regularized horizontal footprints. Force
conservation is checked for all three Cartesian components. Structural horizontal
loads are applied at the silhouette quadrature heights. Their resulting force and
bending moment are preserved in the beam load map. A separate exact fluid-grid
angular-momentum exchange is not claimed: grid-cell and silhouette heights differ,
and crown torsion and off-centre vertical aerodynamic moments are omitted.
Vertical aerodynamic load is reported relative to gravitational load; a 5% guard
limits use of the horizontal mechanics approximation.

Reconfiguration can change the drag map consistently with speed. The reported
layout envelopes and damage ramps use fixed drag coefficients; their Uref-squared
force scaling is rejected when a reconfiguring crown is supplied. Deformed
positions remain frozen for aerodynamic loading until an explicit damage update.

## Structural connections and roots

The inherited Euler--Bernoulli beam matrices retain self-weight geometric
stiffness. Each root can rotate against a foundation spring. Its selected
stiffness equals its capacity divided by 0.04 rad. Capacity and compliance are
therefore linked assumptions in this checkpoint, not separately calibrated
properties. An individually unstable unengaged tree stops the analysis; possible
stabilization through preloaded connections is outside the implemented state.

Crown contacts engage in compression only. Optional bonded crown links are
explicit proposed connections. Contact displacement is solved through a convex
piecewise-quadratic elastic energy with active-set/Newton iteration and line
search. Optional root links exchange equal-and-opposite rotational couples.
Connections have separate finite-strength thresholds, and broken connections are
removed before the remaining stand is re-equilibrated.

The default crown-contact stiffness is 20,000 N/m and strength 400,000 N. The
reported clonal-root comparison uses 2e9 Nm/rad with a 1e8-Nm capacity and 120-m
connection range. These are hypothetical test values. No claim is made that
natural trees have these links, properties or costs. Connection and root-tissue
construction costs remain absent from the mass matching.

A 5-m horizontal soil grid partitions overlapping root footprints, so a ground
cell is allocated once across competing root plates. Each idealized soil capacity
is reduced by its allocation fraction; the root-tissue bottleneck remains separate.
This is an area-sharing proxy, not a three-dimensional soil constitutive model.
Depth-dependent competition, directional anchorage and collective soil failure
remain open. Initial ownership is frozen during damage: lost trees do not
immediately donate soil capacity to survivors.

## Layouts and fair comparisons

The main patches have 7 x 7 planting positions at 90-m spacing, covering 396,900 m2.
The flat reference has 300-m trees. A central gap removes three planting positions.
Candidate shapes are flat, tapered margins, dome, directional ramp, irregular
heights and a gap. Ramp winds are reversed as a separate case. The ordinary patch
is only about 2.1 reference tree heights across, so it represents a grove with
strong boundary influence. The additional 25 x 5 patch has 125 trees and about
7.5 reference heights of planted fetch, with lateral bypass still available.
Neither case establishes a large, horizontally homogeneous forest interior.

At a given H, the baseline stem diameter is H/35, crown radius 0.18 H, root radius
0.1 H and plate depth 4 sqrt(H/300) m. Both a common maximum-height scale and a
matched aboveground tissue/substrate mass proxy are compared. Matching includes
stem/crown wood, leaf and epiphyte mass under the selected densities; it excludes
water, belowground tissue and link construction. It is not an equal-carbon,
equal-productivity or equal-total-resources optimization.

Matching the dome to the flat patch's aboveground mass raises its peak to about
341.3 m and leaves its corner trees about 187.7 m tall. The tapered case peaks
near 402.0 m. These are experimental geometries, not inferred lunar tree heights.

## Results at the same 20-m/s ambient reference

All values below are conditional demand/capacity ratios. A value of 1 indicates
onset in the chosen anchorage model; it is not a calibrated safety factor.

| Layout/support | Maximum root demand/capacity |
|---|---:|
| Flat, independent structures | 0.642 |
| Flat, crown contacts | 0.632 |
| Flat, selected clonal-root links | 0.630 |
| Tapered, matched aboveground mass | 0.701 |
| Dome, matched aboveground mass | 0.726 |
| Dome, crown contacts | 0.726 |
| Dome, selected clonal-root links | 0.663 |
| Ramp, wind toward rising canopy | 0.719 |
| Same ramp, reversed wind | 0.620 |

For the flat patch, the windward-row mean is about 0.628 and the last-row mean
about 0.422. The same compliant-root tree in a solitary porous-flow control is
about 0.772. That control separates neighbor shelter from the change of foundation
model relative to the previous uniform-wind calculation.

The highest anchorage demand in the matched dome occurs at the short windward
corners. The result follows the coupled allometry: shorter margins also receive
thinner stems and smaller root plates. Thus a lowered margin alone is an
incomplete adaptation. This reduced neutral model does not establish that flat
forests outperform domes under realistic gusts or other allocations of tissues.

A separate dome experiment keeps perimeter root plates at least 30 m in radius
and 4 m deep while retaining exactly the same aboveground geometry and mass.
Maximum anchorage utilization falls from 0.726 to 0.599, and the most loaded
member shifts inward. The increased root/soil demand is recorded. This is a
location-specific engineering hypothesis, not a free improvement or an evolved
optimum. Its high-wind nominal failure lies beyond the displacement boundary.

The nominal first-damage speeds for selected intact stands range around
23--25 m/s. Several cases hit the small-displacement boundary first. The output
therefore records damage speed and domain-crossing speed separately. In particular,
the flat independent case crosses its displacement guard near 24.61 m/s before
its nominal 24.96-m/s damage root. Those numbers cannot certify high-wind survival.

## Damage and propagation

After a root or stem limit is reached, the selected removal rule deletes its
standing aerodynamic/structural contribution. A crown event removes distal crown
area once, with a prescribed 50% planform retention, changed mass/span, rebuilt
drag and changed contacts. Falling-body impacts, arrested falls, deadwood wind
obstruction and subsequent regeneration are absent. These limitations matter for
real forest damage.

In the reference-soil dome, two windward corners are removed at a 24-m/s ramp
sample. The changed airflow/stand is re-solved and is then below the selected
failure limits at that same wind. The subsequent 25-m/s sample reaches the
mechanical domain limit. Both flat reference-soil runs hit their domain boundary
before a resolved damage sequence.

The weaker-soil experiment uses 5000-Pa cohesion and a 0.6 pore-pressure ratio.
With 1-m/s wind increments, the 23-m/s sample triggers losses in three rounds:
5 trees, then 7, then 7. Wind remains fixed while geometry and shelter change.
The next round reaches the displacement guard, with 30 standing trees remaining.
These 30 are an unresolved remaining state, not an estimate of final survivors.

With 0.25-m/s wind increments, a propagating front starts at 22.5 m/s and removes
40 trees across ten rounds before reaching the same kind of domain guard. This
large difference in stopping-point counts means total damage is not numerically
resolved by this checkpoint. Earlier damage at a lower sampled wind allows the
smaller-step sequence to continue farther before the deformation guard intervenes.
Both sequences demonstrate a conditional propagation mechanism; neither yields
an ultimate mortality fraction or collapse probability.

## Numerical checks and sensitivity

The focused suite has 32 new tests plus the 25 unchanged earlier wind tests.
Controls cover a homogeneous analytic drag/reservoir equilibrium, zero drag,
projection/divergence, wind-reversal symmetry, conservative force transfer,
Hermite load force/moment identities, shared soil, unilateral contact, root-link
load transfer, broken links, domain stops and the immersion output guard.

The 13 layout/direction flows produce 17 structural comparisons and 827 individual
records. Their largest normalized steady momentum residual is about 1.27e-5;
fluid-to-tree force exchange closes to about 1.61e-15 and global horizontal
structural force balance to about 8.96e-12. These check numerical equations, not
material or ecological validity.

For the flat patch, refining 40 x 32 x 20 to 60 x 48 x 30 changes maximum root
utilization by about 1.81%. The coarser-grid result is also about 1.78% above the
base value; a monotonic convergence rate is not established. Expanding horizontal
dimensions by 50% at unchanged cell size changes the maximum by about 3.94%.
Raising the lid from 1200 to 1800 m changes it by about -0.045%.

The largest departure from target wind inside the strongly forced reservoir is
about 13.75% of reference speed for the base box and 8.39% for the expanded box.
These are boundary-condition diagnostics, not a uniform uncertainty interval.
Further domain/open-boundary work is required. Half/double the prescribed mixing
length changes maximum anchorage utilization by about -1.52%/+2.39% in the tested
flat case; this narrow sensitivity does not validate a turbulence closure.

A fresh repeat of the base flat flow reproduces the saved velocity and force
arrays exactly in the preparation environment. All 57 focused tests passed again
on the final numerical sources. The other scenarios were not all rerun twice.

Full-repository environmental/protection regressions, empirical forest benchmarks,
measured material calibration and full-source manuscript admission remain unrun.
The earlier wind modules are unchanged and their focused tests pass.

## Reproduction

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test*forest*.py' -v
OPENBLAS_NUM_THREADS=1 python research/run_forest_patch.py
```

The runner writes ordinary numerical files to `research/runs/forest_patch` by
default. `--out` chooses another directory outside `immersion/`. Optional
`--stage cases`, `--stage sensitivity`, `--stage damage` and `--stage extras`
separate the work without changing equations. All stages are needed for the full
checkpoint. Configurations, source hashes, per-tree tables, force fields and
stopping reasons are retained. A raw NPZ field is research output, not a renderer
asset or a change to the immersion team's work.

## Next discriminating work

First retain the spatial interaction layer while replacing the small-displacement
beam/contact geometry with a nonlinear formulation. Add event-localized wind
continuation, collision/entanglement and fallen-body obstruction before interpreting
total damage. Establish convergence for more than one geometry and replace the
periodic reservoir with a tested open-boundary treatment. Add resolved/parameterized
gust forcing and compare Earth-condition cases with admitted forest experiments.

Then evaluate local allocation of roots, stem taper and crown porosity under
hydraulic, carbon, resource and regeneration constraints. A developmental or
evolutionary simulation must include reproduction, competition and recruitment;
these shape comparisons alone do not predict a dome or any other evolved canopy.

## Source access and attribution

The inherited tree mechanics and lower-air implementation are reused by import.
No numerical trait values were newly calibrated from external literature.

Schelhaas et al. (2007), *Introducing tree interactions in wind damage simulation*,
Ecological Modelling 207, 197--209, DOI 10.1016/j.ecolmodel.2007.04.025. Publisher
abstract and available introduction/excerpts were inspected. They support the
research distinction between individual trees, shelter and support. This code is
not a reimplementation or validation of ForGEM-W. Full-source admission is pending.

Dupont and Brunet (2008), *Impact of forest edge shape on tree stability: a
large-eddy simulation study*, Forestry 81, 299--315, DOI 10.1093/forestry/cpn006.
Publisher abstract metadata were inspected. Its turbulent edge results motivate a
future gust-resolving comparison; they do not calibrate this neutral mean-flow
screen. The direct full-page attempt failed. Full-source admission is pending.
