# Atmospheric electricity on the Open Moon

The author adopted this modeling direction on 2026-10-02. It joins cloud
microphysics in [climate/](../../../climate/crm/README.md), the air column and
ionization in [atmosphere/](../../../atmosphere/README.md), and the transmitted
light and particles set by [protection/](../../../protection/README.md).
This is an agreed research plan. No graupel comparison, electrical model or new
simulation result is supplied by this study yet.

The immediate question is whether the corrected runs increased or decreased
graupel abundance, changed its sizes, or changed the volume and duration of
mixed-phase conditions. The committed summaries do not settle any of those
signs; an unmeasured change is not a neutral result. The wider question is
whether the deeper atmosphere and longer particle residence can offset weaker
local charging and slower charge separation.

Atmospheric electricity covers the wider subject. The global atmospheric
electric circuit describes the large-scale current system linking electrical
generators, the surface and conducting upper air. Lightning and transient
luminous events (TLEs, including sprites, halos, jets and elves) are parts of
the investigation. Their rates, heights and appearance remain open.

## Framework and geometry

| Task | Framework | Role and boundary |
|---|---|---|
| Graupel, particle sizes and mixed-phase conditions | Existing CM1 output, then matched 3-D CM1 boxes | Keep the current lunar CM1 and Morrison microphysics as the baseline. Diagnose microphysics before choosing an electrical implementation. |
| Locations, lunar phases and environmental forcing | Existing CM1 rings and GCM products | Select representative environments and supply large-scale forcing. Retain uncertainty where the GCM and CM1 disagree. |
| Charge separation, electric fields and lightning | Evaluate WRF-ELEC as a candidate 3-D storm framework | Published configurations include charging, a field solver and discharge treatment. Lunar adaptation and comparison with CM1 are future work. |
| Conductivity, the global circuit and TLEs | Ionization/conductivity columns, followed by suitable circuit and discharge-field models | Use storm electrical output and the shielded atmosphere as inputs; select and validate implementations at that stage. |
| Independent background climate | ROCKE-3D, if the separate climate investigation proceeds | Could constrain forcing and favorable environments. It does not directly answer the present graupel or charging question. |

Three-dimensional boxes are the primary follow-up geometry. Finite storm and
charge structures, flow around terrain, and separation in both horizontal
directions matter. The rings remain useful for environmental sampling; their
missing horizontal direction cannot establish a finite storm's electric field.
Ring integrals must state their assumed transverse width or be reported per
unit width, never silently converted into planetary storm volumes.

Begin with a convective equatorial environment. Retain the corrected highland
box as a dry comparison, and broaden to other environments when the diagnostics
justify it. The current highland run's small rainfall makes it a poor sole
test of active graupel charging. Box size and lateral electrical boundaries
will need sensitivity checks: periodic copies of a storm can affect its field.

## 1. Recover the microphysical comparison

Inventory the actual CM1 control files, raw fields, statistics and restarts
before specifying a rerun. Full output lives outside Git; access to the
committed summaries alone does not establish availability of every field.
[cm1_run.py](../../../climate/crm/cm1_run.py) selects Morrison double-moment
microphysics with graupel (`ptype=5`, `ihail=0`). CM1's `output_q` supports mass
and number fields for double-moment schemes; `output_fallvel` supports Morrison
fall speeds. Confirm each saved case's variables, units and switches [S1].

The first comparison should report:

- Graupel mass and number, vertical profiles and characteristic particle sizes,
  using the scheme's actual distribution and density assumptions. Bulk moments
  imply a size distribution; they do not track individual particles.
- Mixed-phase volume by temperature, distinguishing liquid-plus-ice overlap
  from the subset containing graupel, ice crystals and supercooled liquid.
  State and vary the condensate thresholds.
- Mixed-phase volume integrated over time, its fraction of the sampled domain,
  and occurrence through lunar phase. Report storm persistence separately from
  particle residence, which requires budgets, tracers or sufficiently frequent
  histories.
- Updrafts, differential particle fall speeds, liquid-water supply, riming and
  ice-collection rates where available. Distinguish local rates from totals
  integrated over a storm's volume and duration.

For mass, integrate mixing ratio with air density and physical cell volume.
For number, use the recorded number variable's units. Account for stretched
levels and terrain. Compare common sampled heights, phases and post-spin-up
windows, reporting any coverage mismatch.

Most current full-field output is three-hourly. It can support broad abundance
and overlap statistics; short-lived overlap and residence are undersampled.
If needed, rerun selected restart windows with process budgets accumulated
inside the model and more frequent saved fields. Establish an adequate output
cadence by checking convergence. Electrical evolution must eventually be
integrated at suitable model time steps, independently of how often files are
written.

**First deliverable:** an inventory and a before/after table with sampling
limits, identifying which comparisons are controlled and which are descriptive.
Missing variables or an unavailable baseline remain explicit gaps.

## 2. Run matched 3-D boxes where the inventory leaves a gap

Record the correction under test for each pair; several different changes have
been discussed as corrections:

| Comparison | Existing anchor | What a controlled pair can establish |
|---|---|---|
| Corrected GCM Sun, tilt and radiation feeding CM1 | Older `A28_dim5` forcing versus `A28_dim5_moon`; the corrected flat rings | Response to the specified forcing change. Existing rings also differ in some setup choices, and the old equatorial boxes retain earlier forcing. Inventory those differences before attribution. |
| Terrain forcing at each column's own height | `box_highland` and `box_highland_own_height` | Response to the terrain forcing correction within the documented pair. Keep its dry environment distinct from the equatorial storm case. |
| Lunar particle fall speeds | `box_0e_small` and `box_0e_small_earth_fall` | Sensitivity to particle fall speeds under the same older forcing. This is not a complete Earth-versus-Moon climate comparison. |

Hold site, surface, forcing, grid, microphysics and initialization procedure
fixed except for the treatment being tested. Use compatible restart states
for short controlled tests; give changed forcing enough adjustment time before
interpreting an equilibrium response. Match the lunar phase, quantify internal
storm variability, and retain the executable, input and forcing hashes.

Check selected active periods at finer horizontal and mixed-phase vertical
resolution. The current 6-km grid can supply an initial experiment; its adequacy
for graupel and charging statistics has to be demonstrated. Measure elapsed
time and memory on a short pilot before extending a case. Runs must fit the
standing laptop limit (8 cores, 31 GB shared with other work, sessions of hours
with rests; no rented compute).

Slower settling alone does not determine airborne graupel abundance. Residence,
growth, evaporation and transport also change. Larger mixed-phase volume can
increase the opportunity for collisions while local collision and separation
rates fall. Measure both effects. Keep the fall-speed control distinct from a
future controlled gravity comparison of the full atmospheric column.

**Second deliverable:** resolved changes in mass, number, inferred sizes and
mixed-phase volume/time, with variability and resolution limits. A charging
proxy may rank cases once its assumptions are stated; it cannot supply a flash
rate.

## 3. Evaluate an electrified storm model

Keep CM1 as the microphysical baseline while evaluating WRF-ELEC [S2, S3].
Sun et al. describe a configuration with NSSL double-moment microphysics,
collision-based charging, a Poisson electric-field solver and parameterized
discharges [S3]. This supplies a concrete candidate for the missing electrical
physics. Changing microphysics also changes the storm, so establish a matched
meteorological comparison before attributing differences to electricity.

Pin the candidate implementation and reproduce an Earth benchmark first.
Then audit gravity, sedimentation and collision velocities, gas composition,
number density, charge-transfer laws and their laboratory ranges. Express
breakdown conditions in the applicable gas-density regime rather than copying
an Earth height profile. Check charge conservation through transport and
microphysical transfers, field-solver boundaries, conductive leakage and charge
removal by discharges. Compare electrical behavior across plausible charging
parameterizations.

A separate lunar WRF-ELEC box and transferring electrical physics into CM1 are
implementation options to assess after that audit. Neither is an installed
capability or a simple post-processing switch for the current snapshots.
Record charge structure, field strength, currents and charge-moment changes
alongside any parameterized flashes; the wider electric system needs those
outputs.

**Third deliverable:** an applicability and resource assessment with benchmark
results, followed by a validated pilot if feasible. A framework name alone
does not establish lunar lightning frequency.

## 4. Extend to the conducting upper atmosphere and global circuit

Start with day/night ionization and conductivity columns using the modeled
neutral atmosphere, composition and the adopted protection scenario. Account
for transmitted solar radiation, particle access, background ionization,
attachment, recombination and aerosol/cloud effects as applicable. Existing
exospheric ion-production bounds do not provide a troposphere-to-ionosphere
conductivity profile. The appropriate conducting upper boundary must emerge
from that assessment.

Combine representative storm currents with column resistance and surface
conductivity to investigate a regional and then global circuit. Scaling a box
to the Moon requires the occurrence and area of each storm regime; a few rings
do not determine that population. Carry the GCM–CM1 disagreement into the range
of possible global generator strengths.

Use storm electrical histories to test upper-atmosphere responses. Sprites and
halos motivate a quasi-electrostatic response calculation [S4]; elves need an
electromagnetic-pulse treatment, while jets require their own discharge
propagation treatment. Detailed streamer physics has additional scale and
gas-density requirements [S5]. Choose those solvers when the source histories
and conductivity are available. The neutral CM1 domain reaching 150 km does
not provide those electrical processes, and familiar Earth event heights are
not lunar predictions.

**Fourth deliverable:** conditional circuit and TLE feasibility bounds, with
the source storms and conductivity assumptions attached. Event rates and
rendered appearances remain later outputs.

## Work order, records and evidence

The next task is the saved-output inventory and comparison design in step 1.
This commit records the plan; analysis, instrumented runs and electrical
implementation remain future work. The broader question of the GCM–CM1 surface
climate discrepancy retains its separate decision and uncertainty range.

CM1 diagnostics and any model extensions belong in their domain folders.
This study holds the cross-domain comparison and interpretation. Future
products follow [shared/README.md](../../../shared/README.md), including schema,
producer and input/model hashes, units, sampling, evidence state and limits.
Report what the calculations actually resolve before using them in the
manuscripts or immersive world.

Sources S1–S5 and their reading scope are recorded in [sources.json](sources.json).
They support framework selection and identify follow-up reading; they do not
constitute full scholarly source admission or validation for lunar conditions.
