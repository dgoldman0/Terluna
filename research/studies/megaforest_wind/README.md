# Megaforest wind: first computational checkpoint

## Scope and status

This checkpoint implements lower-air inputs, a static tree-mechanics screen and a
conditional one-dimensional mean-canopy flow calculation. It compares selected
architectural traits and produces inverse wind-load thresholds. Wind climate,
gust dynamics, evolution, biological viability and external support engineering
remain separate work. All material, root-soil and architectural parameters are
hypotheses. The tests establish numerical properties of the stated models.

Work is confined to atmosphere, climate, biosphere, research and tests. The
existing immersion atmospheric-profile format is consumed through a read-only
adapter. No immersion source, renderer, simulation, asset or generated output is
changed. The immutable ensemble seeds and historical numerical baselines remain
unchanged.

## Reproduce

From the repository root, using the existing NumPy/SciPy dependencies:

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest atmosphere/tests/test_lower_air.py climate/tests/test_canopy_flow.py biosphere/tests/test_megaforest_mechanics.py research/studies/megaforest_wind
OPENBLAS_NUM_THREADS=1 python -m research.studies.megaforest_wind.run
```

The default destination is the already ignored `research/runs/megaforest_wind/`.
It contains all 480 envelope cases, 3,360 imposed steady-load samples, 10 canopy
boundary-value solutions, 36 beam-refinement cases, full configurations and
source hashes. Compact reference cases and a check record are retained in
`research/studies/megaforest_wind/results/`. Each case records the wind reference height.
The run is deterministic for a fixed numerical environment. A missing input is an
error, and output paths inside `immersion/` are rejected.

An existing A1 atmospheric-profile export can supply the thermodynamics:

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.megaforest_wind.run \
  --a1-profile /path/to/existing-shared-profile.json \
  --out research/runs/megaforest_wind_a1
```

This command reads the supplied export; it neither rebuilds it nor runs the
immersion team's generators. The current check record covers synthetic adapter
fixtures. An end-to-end run against a real A1 export remains pending because the
referenced generated JSON is absent from the inspected Git tree.

## Ownership and interfaces

`atmosphere/lower_air.py` provides `LowerAir.sample(z)` and
`A1Profile.from_file(path).sample(z)`. Both return height, pressure, temperature,
density, gravity and vapor mixing ratio in SI units. Callers can provide another
object with the same sampling interface. Viscosity and acoustic properties have
not been added in this checkpoint.

`climate/canopy_flow.py` provides a selected mean-flow profile and its momentum
ledger. `biosphere/megaforest_mechanics.py` owns plant geometry, masses, distributed
wind loading, elastic response, prestress and idealized failure capacities.
`research/studies/megaforest_wind/run.py` coordinates cases and writes provenance. It
contains no independent copy of either physics model. Engineering interventions
and manuscript changes are deferred until the physical requirements are better
constrained.

## Lower atmosphere

The standalone atmosphere has surface pressure 121,590 Pa, temperature 288 K,
and a dry N2/O2 mixture with oxygen mole fraction 0.175. Radius, GM and molecular
constants follow the existing `atmosphere/thermal_column.py` lineage. These are
scenario inputs. The baseline is isothermal; a prescribed lapse or inversion and
a fixed vapor mixing ratio are available through `AirConfig`.

The model integrates

    d(log p)/dz = -g(z)/(R_mix T(z)),
    g(z) = g0 [R/(R+z)]^2,
    rho = p/(R_mix T).

The spherical isothermal control is analytic. Nonzero prescribed lapse rates use
16-point Gauss quadrature. The default lower-air domain ends at 5 km. Temperature,
height and mixture checks reject unsupported inputs.

The A1 adapter preserves its documented log-pressure and linear temperature /
vapor interpolation. It uses A1's Rd=287.05 and Rv=461.5 J/(kg K), rather than
silently replacing that closure with the standalone N2/O2 mixture. Its domain,
source hash and planetary parameters remain explicit. The adapter borrows the
existing schema and thermodynamic conventions; it performs no optical transport,
cloud simulation or climate inference.

## Conditional canopy flow

For a canopy of height Hc, the selected momentum balance is

    d[rho K dU/dz]/dz = (rho/2) Cd a_front U |U|,
    U(0)=0, U(Hc)=Uref, K=alpha Uref Hc.

The positive dimensionless mixing coefficient alpha is imposed. No prognostic
turbulence closure, pressure-gradient forcing, buoyant convection or gust
penetration is solved. The zero-wind solution is the zero-velocity limit.
The baseline alpha is 0.04. Area density has a smooth crown distribution above
0.4 Hc, with LAI 22 and an explicit frontal projection factor 0.3. The drag
coefficient uses the displayed one-half convention.

The numerical solve uses normalized speed and stress, retains variable air
density, and checks that the integrated drag balances the difference in boundary
momentum flux. An imposed logarithmic extension supplies wind above the canopy;
its displacement height and roughness are independent assumptions. A sheltered
tree shares the canopy top; an emergent rises above a stand with Hc=0.65 Htree.
The reported Uref is the speed at that stand's canopy top, so emergent and
sheltered reference speeds must be compared with their heights attached.

Tree crown drag and the stand's bulk area-density closure are presently coupled
one way. Changing one tree's traits does not automatically redesign the stand.
The selected mixing/drag parameters cannot establish lunar wind statistics or
validate a particular degree of shelter.

## Static plant mechanics

A solid, circular, linearly tapered trunk is represented by Euler--Bernoulli beam
elements. The baseline uses 60 elements, a tip/base diameter ratio of 0.12,
Young's modulus 8 GPa, tensile/bending and compressive thresholds of 20 MPa, and
wood density 650 kg/m3. These are illustrative inputs, with no species calibration.

Each segment carries the drag of the exposed stem and a selected fraction of an
elliptical crown silhouette. Leaf area contributes to leaf mass; it is not
substituted for the crown's frontal drag area. Crown woody mass, leaf mass,
epiphytes and retained water enter the distributed gravitational loading. At the
assumed water density, a millimeter of stored water contributes 1 kg per square
meter of crown footprint. Canopy hydrology and water loss are not solved.

Crown reconfiguration multiplies its drag by

    response(U) = max(floor, max(1, U/U0)^V).

The rigid control has V=0. The reconfiguring case uses V=-0.8, U0=5 m/s and a 0.2
floor. Only crown drag is reduced. This imposed response represents aerodynamic
streamlining without solving deforming individual leaves or branches. The floor
prevents an unlimited extrapolation of the reduced power law at high wind.

Element matrices include geometric stiffness from the mass above each segment.
The smallest load multiplier for buckling is computed through a generalized
eigenproblem; cases already unstable under their own weight are rejected before
wind thresholds are calculated. For stable cases, the small-displacement system
is `(K_elastic - K_geometric) displacement = wind_load`. Stress evaluation includes
axial compression and bending; the latter contains both direct wind moment and
the additional moment of weight acting through lateral displacement.

The reported displacement-domain boundary is reached at a lateral displacement
of 0.1 tree height or a rotation of 0.2 rad. These are explicit screening bounds,
not universal biological or engineering tolerances. A calculated failure beyond
this boundary is marked as an extrapolation requiring a nonlinear solve.

### Anchorage and crown damage

The root model is an idealized requirement screen. A selected root-soil plate
has a weight moment and a mobilized shear moment. Shear strength is prescribed as
`c + sigma_effective tan(phi)`, with explicit pore-pressure, area-mobilization and
lever-arm coefficients. The root tissue provides another capacity based on its
selected tensile area, strength and lever arm. The smaller of the tissue and
soil capacities governs this single load path. Their capacities are not summed.

These expressions do not replace root architecture, soil constitutive behavior,
foundation rotation, suction, cyclic weakening or actual uprooting measurements.
The assumed plate-weight lever arm also requires calibration. Root-soil parameter
uncertainty is a primary result limitation, especially under lunar gravity.

Crown damage uses a representative supporting branch. The crown mass and drag
are divided among a stated branch count; gravitational and wind bending moments
are combined as perpendicular components. This is an onset proxy. It does not
resolve a branch network, contact, torsion, shedding, regrowth or redistributed
loads after failure. Later thresholds retain the original intact tree and must
be read as counterfactual load requirements after the first damage threshold.

## Experiments and what they show

Five heights (100--500 m) are combined with eight trait scenarios, two soil
settings, dry / 30-mm-water crowns and three prescribed exposures. Baseline trunk
diameter scales as H/35, crown radius as 0.18 H, root-plate radius as 0.1 H, and
plate depth as 4 sqrt(H/300) m. A separate family keeps trunk diameter at 8 m to
expose the consequences of height growth without proportional thickening.

The eight scenarios are rigid reference, reconfiguring crown, wider root plate,
thicker stem, more porous crown, smaller branches, a combined design, and fixed
8-m diameter. These are matched design sensitivities, not an evolutionary
simulation or equal-carbon/equal-productivity optimization. Changes in trunk size
also change mass, branch size and root-tissue area under the declared allometry.
Reducing frontal area at fixed LAI leaves its biological attainability unresolved.

For the 300-m reference in uniform steady wind and the selected reference soil,
the root threshold is about 23.55 m/s. The reconfiguring case reaches about
40.42 m/s, the more porous crown about 30.73 m/s, and a 25% thicker trunk about
23.56 m/s. All four first thresholds occur inside the stated displacement range.
The nearly unchanged threshold after trunk thickening reflects the governing
anchorage bottleneck under these assumptions. These values are conditional model
outputs, not lunar storm predictions or validated wind tolerances.

Widening the root plate shifts the nominal first limit to the trunk, but that
threshold is beyond the small-displacement boundary. The combined design also
reaches that boundary first. Reporting either as a validated high-wind survivor
would overstate this model. The 500-m fixed-8-m-diameter wet case already crosses
the self-weight buckling boundary, whereas proportional thickening carries a
large and explicitly recorded biomass cost.

The standalone air density changes from approximately 1.45786 to 1.44373 kg/m3
between 0 and 500 m, about 0.97%. Under this particular lower-air scenario,
architecture and anchorage uncertainty are much larger effects than the density
change through the forest. This comparison gives no information about the wind
shear or gust distribution across the same height.

## Numerical checks and remaining work

The new suite has 25 passing tests. It checks spherical hydrostatics and its
pressure derivative, composition and humidity, A1 schema/interpolation fixtures,
analytic cantilever loads, Euler buckling, canopy momentum balance, zero-drag
flow, drag-response limits, root load-path bookkeeping, wet-mass effects,
threshold roots, right censoring, domain flags, grid sensitivity and the
immersion output guard.

Among the 10 canopy cases, the largest relative momentum residual is about
1.41e-10. Across the selected 36 refinement runs, the largest threshold difference
between 60 and 240 elements is about 0.474%. Those comparisons cover selected
geometries; they do not bound all model or parameter uncertainty. The preparation
worktree contains the scoped additions, so the historical whole-repository /
protection check and the existing environmental regression suite were not run.
The new work does not mark those checks PASS.

Of the 480 envelope cases, 6 are rejected for self-weight buckling, 18 reach a
damage threshold at zero wind, and 209 have a first nominal failure beyond the
displacement domain. Read the machine record for the complete classification;
none of these counts represents a prevalence or probability in a lunar forest.

The next mechanical step is a geometrically nonlinear beam with rotational
root-foundation compliance, followed by a branch network and explicit changes
in exposed crown area. That should precede stronger claims about the high-wind
combined designs. Root-soil calibration and sensitivity need to develop alongside
it. Independent climate work can then supply time-dependent pressure, temperature,
wind and stability fields; canopy turbulence and gust response will connect those
fields to dynamic tree loading. Carbon, hydraulics, long-night survival and
reproduction remain constraints on which mechanically attractive traits can
actually belong to a viable forest.
