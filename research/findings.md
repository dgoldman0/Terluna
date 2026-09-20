# First atmosphere, climate and long-night requirement calculations

**Scope:** executable research development against repository commit `79a60a1a51b4552ac9ae8e8bb295b567d9254995`, 20 September 2026. These are bounded theoretical and numerical results. No new empirical observations, complete aeronomy model, calibrated lunar climate or biological viability result is reported. The five papers and accepted seeds remain unchanged.

## Results in brief

This pass implements a solved spherical molecular thermal column, an explicit band-heating interface, a conservative periodic spatial thermal screen, and carbon-reserve/oxygen requirements models. All four are new code, separate from the pinned historical baselines.

The largest physical gap is now more specific: the inherited solar table begins at 202 nm, while the inherited filter calculation leaves about 24.8–120.18 nm outside its evaluated material-response treatments. A recovered historical solar band reference places at least 37.5% of its 0.1–118 nm energy entirely within that optical-response gap. Closing the full spectrum–chemistry–escape loop requires actual response data and species processes in that range.

Numerical counts are 232 molecular-column cases, 54 climate cases, 180 carbon/capacity cases, 20 oxygen cases, 18 climate-to-carbon examples, and 27 assumed band-response combinations. Every molecular BVP converged numerically; 73 cases carry a kinetic-escape warning. All cases retain their molecular-model and uncalibrated-boundary flags. The checks cover mathematical and numerical consistency, with full environmental compatibility still open.

## 1. Spectral inputs and the energy interface

The two inspected archive files have these actual ranges:

| Input/treatment | Range |
|---|---|
| TSIS1_HSRS_stride100.csv | 202–2729.9 nm |
| TiO2_Siefke.yml optical indices | approximately 120.181–125122.762 nm |
| Inherited independent-atom X-ray treatment | applied at energies >=50 eV, corresponding to wavelengths <=24.797 nm |

The X-ray limit is the inherited model's application boundary. Raw scattering-factor tables may contain entries at lower energies; their existence alone does not establish reliable coating response across the uncovered interval. The audit makes no claim that the material is transparent there. It establishes that the current input/model combination cannot justify a rejection fraction there.

Ribas et al., arXiv:astro-ph/0412253v1, Table 4, Sun column, gives integrated 1-AU fluxes of 0.15, 0.70, 2.05, 1.00 and 0.74 erg s^-1 cm^-2 in 1–20, 20–100, 100–360, 360–920 and 920–1180 Angstrom bands. The sum converts to 0.00464 W/m² over 0.1–118 nm. Section 3.5 describes a 1993-based mid-cycle composite with complementary and inferred bands. This is a historical benchmark, not a measurement of today's Sun, and excludes longer-wavelength FUV including Lyman-alpha.

The 36–92 and 92–118 nm bins lie wholly within the response gap: `(1.00+0.74)/4.64 = 0.375`. Counting the whole partly intersecting 10–36 nm bin gives the upper bound `3.79/4.64 = 0.81681`. No within-bin spectral distribution is assumed. These bounds concern incident energy in this reference only, not absorbed heat or escape.

For supplied bands the interface evaluates

`q_heat = sum [F_band (r_abs/R)^2 / 4 * transmission * absorptance * heating_fraction]`.

All quantities are energy fluxes normalized to the Moon's surface area. Absorbed power is partitioned into sensible heat and a recorded other-energy term. Missing parameters, gaps, overlaps, invalid fractions and altered pinned optical inputs raise errors. No unresolved band is set to zero. The 27 response combinations use explicit hypothetical transmissions, absorption radii and heat fractions; they are sensitivity calculations, not candidate-filter performance or self-consistent absorption-altitude solutions.

## 2. Solved molecular thermal column

### Boundary conditions and equations

The lower profile retains a prescribed dry-adiabatic/isothermal construction. Surface pressure is 121590 Pa, surface temperature 288 K, and the reference mixture is 82.5% N2/17.5% O2 by number. Base pressures 0.001, 0.1 and 10 Pa and base temperatures 150, 180 and 210 K are swept. Those base temperatures are inputs; surface habitability and the lower climate have not been solved by this upper-column model.

Above the base, the unknowns are radius and temperature versus log pressure, exobase pressure, molecular escape rates and the lower energy flux. Hydrostatic balance and steady conservation give

`dp/dr = -rho GM/r²`,

`E(r) = (r/R)² [-kappa dT/dr + sum_i rho_i v_i (c_p,i T - GM/r)]`,

`E(s) = E_base + q_heat C(s)`.

Here `s = ln(p_base/p)/ln(p_base/p_exo)` and `C(s)` is an explicitly prescribed cumulative deposition shape. Low, middle and high deposition shapes are compared. These are distributions in normalized log pressure, not optical absorption calculations. The fluxes `j_i=(r/R)² rho_i v_i` are constant in this model. Momentum acceleration and kinetic energy of the bulk flow are omitted; Mach and escape-parameter flags monitor important breakdowns.

The upper boundary sets `Kn=1` using the mixture pressure scale height and the inherited collision proxy. Each molecular flux equals a Jeans flux. The energy at infinity carried by escaping diatomic molecules is

`e_infinity,i = R_i T [(lambda_i+2)/(lambda_i+1) + 1]`,

where the final 1 is the two rotational degrees of freedom contribution and `lambda_i=GM/(r R_i T)`. The boundary condition is `E_exo=sum_i j_i e_infinity,i`. Gravitational potential and advected enthalpy are retained consistently, so the remaining conductive boundary flux supplies the difference. An independent Maxwellian integral test checks the translational escape-energy expression.

Conductivity is approximated as `kappa = 9.37e-5 T W/(m K)` with multipliers 0.5–2. The coefficient comes from the nitrogen proxy used by Tucker et al. in Pluto modelling, converted from cgs; applying it to this mixture and temperature range is an assumption. The inherited `5.5e-19 m²` collision cross-section is also a proxy, not a measured N2/O2 mixed collision model.

### Conditional numerical results

For base temperature 180 K, base pressure 0.1 Pa, middle deposition and the nominal conductivity:

| Deposited sensible heat, W/m² lunar surface | Solved exobase K | Exobase radius / R | Molecular loss kg/s | N2 Jeans parameter |
|---:|---:|---:|---:|---:|
| 0 | 180.00 | 2.198 | 0.0000071 | 24.03 |
| 0.000001 | 194.39 | 2.284 | 0.000083 | 21.41 |
| 0.000003 | 223.65 | 2.485 | 0.00470 | 17.11 |
| 0.000010 | 316.57 | 3.623 | 14.44 | 8.29 — kinetic warning |
| 0.000030 | 345.91 | 5.759 | 302.23 | 4.77 — kinetic warning |

These are steady solutions of the stated molecular equations, not predictions of the Sun/filter combination or total atmospheric losses. Low-gravity kinetic escape, atomic enrichment, radiative cooling and Earth tides could change them. A Jeans parameter below 15 triggers further kinetic review; that threshold is an indicative warning, not a validated sharp transition. Radius above one-tenth of a nominal lunar Hill scale triggers an additional conservative tide review flag. No stable lifetime is inferred.

At the same heat input of 3e-6 W/m², halving conductivity changes the solved exobase to approximately 268.10 K, while doubling it gives 201.64 K. Lower-atmospheric conditions and transport properties therefore remain material uncertainties. The solver removes the old fixed exobase-temperature floor; it does not remove dependence on physical boundary conditions.

The maximum absolute terminal energy residual in the sweep is about 9e-13 W/m². Boundary and collocation residuals meet the configured tolerance. A temperature-gradient reconstruction independently recovers the conductive flux. Grid checks compare 81 and 161 initial nodes. These checks do not establish physical stability or justify extending a hydrodynamic solution into a collisionless regime.

## 3. Carbon storage: a conditional necessary-and-sufficient result

Let `f(t)=usable_assimilation(t)-prescribed_demand(t)` be periodic with period `P`. A storage reservoir obeys `dS/dt=f(t)` within `0<=S<=C`; surplus may be discarded at `C`. There are no rate limits, state-dependent demand, growth or mortality in this mathematical problem.

A sustainable periodic operation exists exactly when

`integral_0^P f(t) dt >= 0`,

and

`C >= max_{circular intervals I, length(I)<=P} [-integral_I f(t) dt]_+`.

**Necessity.** A negative cycle balance accumulates an unbounded deficit over successive cycles. On any interval of net deficit, the reservoir must supply at least that interval's deficit; initial storage is at most `C`.

**Sufficiency.** Consider the maximum accumulated deficit ending at each phase, starting from any preceding phase. Removing whole earlier periods cannot increase this deficit when the cycle balance is nonnegative, so at most one period needs to be searched. Set the stored quantity equal to `C` minus that running deficit. It evolves with the prescribed net flow except at the upper boundary, where surplus is discarded. The capacity bound keeps the lower boundary nonnegative, and periodic forcing gives a periodic storage trajectory. This constructs a feasible initial phase/state.

The implementation computes all circular intervals in O(N) time using a sliding maximum of cumulative sums. Tests compare it with exhaustive interval enumeration and check that a reservoir just above the calculated limit meets demand while one below it fails. Floating-point feasibility comparisons include a documented numerical tolerance.

### Illustrative lunar-cycle demands

The numerical examples use fixed structural biomass and carbon units equal to one day of reference maintenance. Usable assimilation follows a half-sine daylight curve and averages three reference-maintenance units per illuminated day. The period is 29.53059 days. Temperature, a Q10 response, and nighttime demand suppression are hypothetical traits, not fits to lunar plants.

At 288 K throughout the cycle:

| Nighttime demand / reference day demand | Minimum reserve, reference-maintenance days |
|---:|---:|
| 1.00 | 15.77 |
| 0.50 | 8.38 |
| 0.25 | 4.69 |
| 0.10 | 2.48 |

At night demand 0.25, darkness alone costs 3.69 units. The full requirement is 4.69 because demand also exceeds production near the illuminated cycle's dawn and dusk. A positive monthly balance alone therefore misses a storage constraint even in this simple model.

These figures are requirements to compare with measured usable reserves and consumption rates. The model does not assign traits to species, prove survival, model injury, or treat spilled carbon as demonstrated growth. Nutrient supply, reproduction and community interactions remain separate requirements.

## 4. Aquatic oxygen with finite storage

The mixed box obeys `dO/dt=P(t)-R+k(O_sat-O)`. An imposed dissolved-oxygen ceiling represents rapid additional outgassing; at the lower boundary the unmet aerobic demand is recorded. The exact constant-forcing step integrates exchange and boundary contact times, preserving the oxygen ledger.

The ceiling prevents a bright interval from accumulating arbitrarily high dissolved oxygen to pay for later darkness. The periodic mathematical solution is tested separately for an oxygen threshold and whether all prescribed demand is supplied.

For the illustrative inputs `O_sat=8 g/m³`, ceiling `9.6 g/m³`, respiration `1 g/m³/day`, and daylight-mean production `3 g/m³/day`, an exchange coefficient of 0.05/day leaves unmet demand during part of the cycle. At 0.1/day the minimum is about 0.47 g/m³; at 0.2/day it is about 3.30 g/m³. The selected 2 g/m³ screening threshold is a scenario input, not a biological safety standard. Concentrations in g/m³ are numerically equal to mg/L.

Actual gas solubility, pressure, temperature, salinity, stratification, wind exchange and community demand need calibration before ecological inference. The box has no sediment or explicit nutrient dynamics. The full 20-case output exposes both outgassing and unmet respiration.

## 5. Spatial periodic climate screen

The finite-volume model uses equal-area cells in sine latitude and longitude:

`C_i dtheta_i/dt = Q_i(t) - A - B theta_i + D (L theta)_i`,

where `theta=T-273.15 K`. The symmetric Laplacian has zero row/column sums, so intercell transport conserves global heat. Positive capacities and `B>0` make the symmetrized operator negative definite. The period map is therefore a contraction with a unique periodic solution. Each constant-insolation interval is propagated by matrix exponentials; the periodic initial state is solved directly. This is a mathematical stability result for the specified linear system, not nonlinear climate stability.

The mean is fixed by `T_mean=273.15+[S*tau*(1-albedo)/4-A]/B`. The finite-grid solar quadrature is normalized to the exact spherical mean `S/4`. With `S=1361 W/m²`, transmission 0.90, albedo 0.30, `A=200 W/m²` and `B=2 W/m²/K`, the mean is 280.33 K by construction. Raising A by 20 W/m² lowers it by 10 K. This dependence must remain visible; an attractive mean cannot certify an atmosphere's greenhouse behavior.

At these radiative values, transport D=1, participating atmospheric mass fraction 0.1 and a 10-m water layer:

| Synthetic geography | Min K | Max K | Grid-time fraction 273.15–313.15 K |
|---|---:|---:|---:|
| Dry | 261.18 | 306.32 | 0.661 |
| Half-area concentrated water | 260.66 | 307.72 | 0.719 |
| Half-area alternating water strips | 267.88 | 299.42 | 0.793 |

The two wet cases have equal total water area and depth but different adjacency. They show a thermal-distribution sensitivity under the imposed transport law. Every case in this table reaches freezing somewhere. Latent heat and ice-albedo feedbacks are absent, so these are flagged thermal screens, not liquid-water or habitability results. Geographic water masks are synthetic; no measured terrain was used.

The maximum discrete energy residual is below 8e-11 W/m², and periodic residuals below 2e-12 K. Timestep refinement passes the chosen 0.01 K minimum-temperature check. A separate 6x24 to 12x48 comparison of the concentrated map changes coarse-cell averaged temperatures by up to about 2.04 K (RMS 0.72 K). Spatial uncertainty remains; detailed regional interpretation would be premature.

Eighteen examples pass climate temperature/incident-light traces into the carbon model using one fixed irradiance reference. Latitude-dependent irradiance is retained. The response is still a hypothetical linear assimilation and Q10 demand model; no cold injury, light saturation or photosynthetic action spectrum has been introduced.

## 6. Protection contamination and renewal boundary

The inherited candidate shield exhaust has kinetic power about 1.253e14 W. Multiplying that power by separately imposed intercepted and thermalized fractions and dividing by lunar surface area gives the new plume-heat sensitivity table. At both fractions combined equal to 1e-6, the deposited power is about 3.3e-6 W/m², comparable to the low-power molecular-column cases. This is an energy-accounting sensitivity, not a computed plume intersection, heating profile or demonstrated contamination limit.

The previous mass/energy/renewal budgets remain unchanged. Full orbit/plume transport, source-to-use logistics, manufacturing-network growth and total cumulative consumables were not newly solved in this pass.

## 7. Evidence status and next work

The completed results establish equations, explicit assumptions, tested numerical behavior and conditional thresholds. The unresolved high-value work is now:

1. Close the resolved solar/material-response and species absorption/chemistry/cooling inputs, retaining coverage checks. The historical band benchmark alone cannot provide those quantities.
2. Replace the prescribed lower atmosphere and linear OLR with calibrated radiative-convective calculations, then add water/ice feedback and actual geography.
3. Populate biological reserve, respiration, productivity and oxygen-exchange ranges from fully read sources; add growth, nutrients and life-cycle conditions. Experimental lunar-gravity viability remains a separate evidence need.
4. Test atomic/diffusive escape and an appropriate kinetic upper boundary, and quantify actual plume coupling before integrated retention claims.

Input gaps hold a complete coupled-atmosphere claim open; they do not prevent using these modules for explicit sensitivity and requirement studies. The root's expressive details should continue to inherit only the environmental and biological scope actually supported.

## Reproduction and provenance

Run `python -m unittest discover -s tests -v`, followed by `python research/run_environment_screens.py`. Supply the original protection ZIP with `--protection-archive` to inspect its data. The run records, source-access status and base revision are in [environment_checks.json](environment_checks.json) and [environment_sources.json](environment_sources.json).

The original feasibility script and protection verifier were rerun directly from their mounted archives. The repository-wide `research/check.py` was not executed on a complete checkout in this session: the connector is read-only and container GitHub DNS failed. That limitation is recorded separately from the 39 new tests. No paper, accepted seed, locked planning file, pinned baseline implementation or historical dataset is modified by this patch.
