# SWAN gravity audit

The first wave calculation uses SWAN 41.51 with Komen wind input and
whitecapping, Wu wind drag, and the discrete interaction approximation (DIA)
for four-wave interactions. These paths use the configured gravitational
acceleration. Source review also found a fixed-frequency wind-input cutoff
that changes an Earth–Moon comparison, so the study retains stock and modified
builds as separate parameterization choices.

## Source and license

The [official download page](https://swanmodel.sourceforge.io/download/download.htm)
identifies the release. The audited
[source archive](https://swanmodel.sourceforge.io/download/zip/swan4151.tar.gz)
has SHA-256
`325c229f3dde0b812db4ec9ae6fa1cd4086d2109562c389576590ce8f94b26bb`.
Source headers credit Delft University of Technology and specify GNU GPL
version 3 or later. Preserve those terms and attribution when distributing
modified SWAN code or executables. The original archive and the recorded
patch distinguish upstream code from Terluna's changes.

The [user manual](https://swanmodel.sourceforge.io/download/zip/swanuse.pdf)
and [technical documentation](https://swanmodel.sourceforge.io/download/zip/swantech.pdf)
describe the model formulations. File and line references below refer to the
archive's original `.ftn` and `.ftn90` files, before preprocessing or patching.

## Gravity and the selected physics

| Component | Source evidence | Consequence |
|---|---|---|
| Configured gravity and water density | `swanpre1.ftn:2179–2180` | `SET GRAV` and `SET RHO` reach model parameters. |
| Wave dispersion | `KSCIP1`, `swanser.ftn:694–717` | Deep, shallow and intermediate-water expressions use `GRAV`; regime boundaries are dimensionless. |
| Komen wind input | `SWIND3`, `swancom3.ftn:2068–2087` | Growth uses phase speed, friction velocity, frequency and the air/water density ratio. |
| Wu drag | `WINDP1`, `swancom3.ftn:1014–1024` | Drag is constant below 7.5 m/s and depends on wind speed above it. This remains an Earth-calibrated prescription. |
| Komen whitecapping | `SWCAP`, `swancom2.ftn:2485–2491,2611–2615` | Dissipation uses steepness, mean frequency and relative wave number. |
| DIA quadruplets | `swancom4.ftn:243` | The interaction coefficient includes `1/GRAV**4`. |

The stock initialization fixes air density to **1.28 kg/m³**
(`swanmain.ftn:1740–1743`). After input parsing, the density ratio is reset
using that value and the selected water density (`6624–6627`). The first strip
experiments hold air density at 1.28 kg/m³ and water density at 1025 kg/m³.
These are controlled coupling assumptions. `SET RHO` sets water density;
air density requires the separate source change described below.

## The AGROW sensitivity

AGROW supplies the initial linear wind source from an initially calm surface.
`SWIND0` reduces this source above **1 Hz** by a factor proportional to
frequency cubed in the denominator (`swancom3.ftn:1800–1807`). The source
comments describe the 1 Hz knee as largely arbitrary, chosen to leave
customary calculations ending at 1 Hz unchanged (`1714–1723`). A fixed knee
breaks gravity similarity when frequency bands are rescaled.

The modified build changes only the local frequency entering this reduction,
from `f` to `f * g_reference / g`. The knee becomes `g/g_reference` Hz.
Reference gravity comes from the project's shared constants. At that
reference gravity, the formula preserves stock behavior exactly. Stock and
modified results measure the consequence of this numerical choice;
dimensional similarity supplies its rationale, while the physical lunar
wind-input cutoff remains unresolved.

The runs use explicit nonstationary `INIT ZERO`, which bypasses another
dimensional default: SWAN's usual initial spectrum forces a 2 s period when
its initial significant wave height is below 0.05 m
(`swanmain.ftn:9411–9418`). **The chosen patch leaves initialization unchanged.**
The wind source creates the waves during the run.

## Disabled-breaking defect and the supported configuration

Initial runs with `OFF BREAKING` failed repeatability: identical executables
and inputs produced different wave histories. A repeated Earth control
differed by about 2.5% in significant wave height at some sampled locations.
Those outputs failed the repeatability criterion and were excluded.

Source inspection identified an uninitialized directional-partition value,
`KTETA`. `SINTGRL` declares it as an output (`swancom1.ftn:5706`) and assigns it
through `BRKPAR` only when depth breaking is enabled (`5944–5949`). With
breaking disabled, the routine still calls `FRABRE` after waves appear
(`6001`). `FRABRE` divides by `sqrt(KTETA)` (`swancom2.ftn:1744`). Its resulting
breaking fraction then controls the spectrum limiter
(`swancom1.ftn:7990`). A debug build trapped at this exact undefined read.
An isolated diagnostic with quadruplets disabled, which also disables that
limiter, produced identical repeated profiles.

The production configuration uses the supported command
`BREAKING CONSTANT 1.0 0.73`, retaining default DIA option 2. It initializes
`KTETA` through the normal path. At the deep-water conditions used here,
depth breaking should contribute zero; output checks verify the breaker
fraction throughout the runs. The breaker index remains Earth calibrated;
lunar surf-zone calibration remains open. This configuration uses the
separately recorded AGROW source change.

With that command, the diagnostic repeated Earth profiles matched exactly,
as did stock and AGROW-modified Earth profiles. The corresponding sixfold
gravity-scaling pair's maximum differences were 0.303% in significant height
and 0.326% in mean period over the retained finite-amplitude samples. Stock
differences were 5.44% and 3.91%, respectively. Its largest discrete peak-period
difference was one spectral bin: the roughly 10% frequency spacing allows
the selected maximum to move between adjacent, nearly equal peaks. Peak-bin
differences therefore need their spectral resolution stated alongside them.

## Numerical checks and remaining scope

For a gravity-similarity pair with fixed wind, scale length, depth and time by
the inverse gravity ratio, and frequency by the gravity ratio. Compare
`g Hs / U²` and `g Tp / U` at matching dimensionless positions and times.
This changes basin size as well as gravity. A separate comparison with the
same physical fetch and forcing duration answers the question about an
unchanged basin.

The first cases use Cartesian one-dimensional strips, uniform winds, flat
depth and the first-order backward-space/backward-time propagation scheme
(BSBT). Time and space refinement test its numerical diffusion. Frequency
range and resolution require separate checks: the first wider-band test
retains roughly 10% spacing between frequency bins to test band truncation.
Incoming boundary spectra are set to zero. Currents, bottom friction, depth
variation, vegetation, mud and ice are omitted. Depth breaking is enabled
with its contribution checked as described above.

The tested p90-wind lunar endpoint, at 100 km after 48 hours, changes by less
than 0.3% in significant height and 0.2% in mean period under the separate
time-step, spatial, directional and spectral-band changes. Early growth is
less settled: halving the time step changes significant height by as much
as 33.2% and mean period by 17.0% over the sampled profiles. The height
difference is largest at one hour; the mean-period difference at two hours.
Early growth and the other wind cases require their own refinement checks.

The spectral output marks very small variance densities with an exception
value: `SWCMSP`, `swanout2.ftn:2388–2403`, applies a floor of
`1e-12 m²/(rad/s)`, equivalent to `2π × 1e-12 m²/Hz`. Edge diagnostics can
bound the omitted resolved variance. Significant wave height also includes
SWAN's diagnostic high-frequency tail. Integrating the printed resolved
spectrum alone therefore gives a different value.

Several options need further work before expansion: ST6's FAN and ECMWF drag
paths set Earth gravity directly (`SdsBabanin.ftn90:2305,2466–2471`), and the
finite-depth exact XNL coupling does likewise
(`mod_xnl4v5.ftn90:8540,8563`). The selected configuration bypasses them.
Spherical geometry defaults to Earth size; a geographic lunar model must
change its radius. Gravity-wave dispersion also omits surface tension.

Numerical similarity, convergence and Earth benchmarks establish distinct
checks on the implementation. Terluna sea-state predictions still require
assessment of air–water coupling, wind persistence and direction, real basin
geometry and the empirical growth and dissipation laws at lunar gravity.

## Basin and slope follow-up

The pilot keeps the AGROW change and adds a separate air-density change in
`swanmain.ftn`, anchored to the pinned original file hash. The assignment to
`PWIND(16)` uses the mean moist surface density recovered from the equatorial
CM1 water columns near Smythii–Marginis. SWAN's existing post-parse reset of
`PWIND(9)` then supplies that air/water ratio to Komen growth. The coupled build
occupies a separate directory. The Wu drag law and empirical AGROW source
coefficient retain their existing values; their calibration at lunar gravity
remains open.

The radius is an input parameter. SWAN 41.51's
`swanpre1.ftn:827–828` reads `REARTH` from `COORDINATES SPHERICAL` and derives
the length of a degree. `swancom5.ftn:204–205` uses that length in propagation
metrics, and `1918` uses `REARTH` for geographic turning. The pilot supplies
`MOON_RADIUS` from shared constants. An executable control compares the same
small physical grid in Cartesian metres and lunar degrees, with a 1% tolerance
on height and mean period, covering the geographic options used by this pilot.

The chosen constant-index breaking path assigns `BRCOEF = PSURF(2)` in
`BRKPAR` (`swancom2.ftn:3120–3122`). Its limiting height is set by depth and
that dimensionless index. Dispersion and group speed use configured gravity;
the inherited default index 0.73 remains empirical. The slope experiments
retain depth breaking; quadruplets, whitecapping, bottom friction and triads
are disabled. A small-amplitude, narrow-band executable control checks
`Hs² Cg` conservation to 2% before breaking; the gravity-wave dispersion
relation supplies the independent group-speed calculation.

The [boundary-input documentation](https://swanmodel.sourceforge.io/online_doc/swanuse/node27.html)
allows spectral boundary files. The pilot transfers the final `SPEC1D`
variance, mean direction and spread from its offshore node into stationary
slope cases, retaining frequencies and density values. SWAN reconstructs the
directional distributions from their moments; transferring multiple directional
peaks requires the full two-dimensional spectrum. The separate
[numerical audit](results/pilot_checks.json) checks final timestamps, spectral
edges, stationary stopping and grid differences; it distinguishes Hs in a
region with `Qb ≥ 0.01` from individual-wave height. Wave groups and run-up
require a subsequent individual-wave calculation and a SWASH source audit.

An initial geographic smoke run exposed a data-reading issue: roundoff in
interpolation of exception-marked bottom nodes excluded some adjacent water
nodes. Those diagnostics were rejected. The accepted inputs use finite
negative depths on dry ground, and output the computational grid directly.
The reader checks every wet node's reported depth against the atlas input,
with explicit allowance for output precision, and rejects incomplete wet-node
output.
