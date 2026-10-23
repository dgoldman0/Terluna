# The improved crystal-growth control in the lunar snow zone

The improved McSnow plate control produces **3.14 mm crystals after a 54.18-hour lunar descent** when the existing 3 mm/h water allowance is shared among the prescribed incoming crystal flux. Reducing that flux tenfold and hundredfold gives **9.51 mm and 29.08 mm**, respectively. These trajectories start from identical 4 µm spheres. With continuously saturated air and no supply limit, the same seed instead reaches **147.47 mm in 18.07 hours**. The long lunar residence permits considerable growth, while the amount of water available per crystal determines whether that potential can be supplied.

This calculation follows the author's request to explore lunar behavior with the improved control despite the remaining Earth discrepancies. It transfers the explicit removal of the extra oblate mean-ventilation term, retaining IGF2+ directional growth, deposition density, geometry and native drag. The [preceding Earth diagnosis](benchmark_growth_trace.md) remains the reference for the intervention and its empirical errors. The new [product](results/mcsnow_lunar_growth.json) records 68 executions, including unchanged-source controls, independent gravity and pressure changes, falling trajectories, water allocation, seed sensitivity and numerical refinement. These are conditional single-crystal predictions at the freezing level; aggregation and a distribution of incoming habits are outside this experiment.

## Conditions and native-code reuse

The atmosphere is the existing [301 K reference column](scale_lunar_methods.md), using the same climate input and pressure profile as the SCALE lunar comparison. From −15 °C to freezing it spans **16.2805 km**, between heights **45.5256 and 29.2451 km**. Pressure rises from **46.677 to 65.925 kPa**. Gravity is the shared lunar value, **1.6242188766 m/s²**, held constant through the column. Temperature and nominal water saturation are prescribed, with no vertical air motion.

The [runner](mcsnow_lunar_growth.py) reuses the pinned 2024 McSnow source at revision `1e428d7e901cbc5b91f2040e575f736ec5963a74`; the original `mo_depo.f90` SHA-256 is `a0904ce5e1061b559d9fe3d9c5c14329803dc823c9179fe43a3048f87f402232`. It builds separate native and altered copies for Earth and Moon. The [new caller](mcsnow_lunar_driver.f90) invokes native `sp_depo`, geometry and `vterm_snow`. Lunar gravity replaces the source gravity parameter; the Earth control retains the source's 9.81 m/s² for exact benchmark parity. Native air transport and gas constants remain in use, so this experiment does not apply a new lunar-mixture transport correction.

The climate input is `research/studies/climate_targets/results/climate_targets.json`, SHA-256 `80a85def114b039cef9c9acd97b18d3bd47413f51f6693d1cc6f5e850ebdbea8`. Assembly verifies its identity and the continuous profile against all 30 saved cell-centre pressures in the preceding lunar product. That product supplies the **10,232.4266 crystals/m²/s** incoming number flux and **3 kg/m²/h** deposition allowance. Replacing its incoming distribution with equal spherical seeds is a declared sensitivity, rather than a reproduction of its population.

## Pressure and gravity act in opposite directions at a fixed age

The thirty-minute controls start from the Earth benchmark's 2 µm radius sphere at exact native water saturation. The local lunar pressures are 47.344 kPa at −14.4 °C and 58.483 kPa at −5.3 °C. Holding each temperature and pressure fixed isolates gravity's effect on fall speed and ventilation.

| Gravity and pressure | Plate width at −14.4 °C | Column length at −5.3 °C |
|---|---:|---:|
| Earth, 101 kPa | 5.643 mm | 1.838 mm |
| Earth, local lunar pressure | 9.030 mm | 2.261 mm |
| Moon, 101 kPa | 4.009 mm | 1.040 mm |
| Moon, local lunar pressure | 6.455 mm | 1.268 mm |

Lower pressure increases uptake, while weaker gravity reduces settling and ventilation. Their combined effect on the plate is a 14.4% increase over the Earth benchmark at thirty minutes. The column becomes about 31% shorter. The source change concerns oblate mean ventilation; native and altered column histories agree in this control. Lower gravity therefore does not increase every crystal's size at a given age. Its large additional opportunity for growth appears through the duration of the fall.

## Descent separates growth potential from water supplied

The matched Earth comparison compresses the same temperature and pressure history into a **2.696 km** cloud, scaling depth by the lunar-to-Earth gravity ratio. This preserves the imposed hydrostatic profile in temperature coordinates. It is a controlled Earth-gravity counterpart to the lunar cloud, distinct from the earlier generic Earth winter column.

For the supplied cases, let the incoming number flux be `F`, fall speed `v`, column depth `H`, and imposed deposition allowance per area and time `Q`. A steady collision-free population has local number density `N=F/v`. Distributing the allowance uniformly through depth gives

\[
\dot m_{\rm allowed}=\frac{Q}{H}\frac{v}{F},\qquad
\Delta m_{\rm allowed}(s)=\frac{Q}{F}\frac{s}{H},
\]

where `s` is distance fallen. The caller limits positive native uptake to this rate **before** the native habit, volume and mass updates. It therefore conserves the supplied deposition account while preserving the native response to admitted growth. The allowance bounds total deposited mass per crystal at **0.0814404 mg** for the original number flux. It is an external deposition budget, not a solved vapor reservoir; supersaturation, cooling and air motion do not respond dynamically.

| Improved control; 4 µm spherical seed | Final maximum dimension | Final mass | Time to freezing | Final fall speed |
|---|---:|---:|---:|---:|
| Earth counterpart, unlimited supply | 14.306 mm | 1.594 mg | 1.688 h | 0.551 m/s |
| Moon, unlimited supply | 147.471 mm | 231.865 mg | 18.071 h | 0.267 m/s |
| Earth counterpart, 3 mm/h allowance | 3.346 mm | 0.07915 mg | 2.588 h | 0.389 m/s |
| Moon, 3 mm/h allowance | 3.137 mm | 0.08127 mg | 54.177 h | 0.130 m/s |

At the original number flux, the unlimited lunar path would deposit **8,541 mm/h** of water equivalent, approximately **2,847 times** the prescribed allowance. It also requires **6.73 MW/m²** of latent-heat removal to maintain the imposed temperature. Its 14.7 cm crystal represents sustained isolated growth potential under those conditions. The supplied lunar path deposits **2.994 mm/h** and releases **2.36 kW/m²** of latent heat across the whole column. This latter energy demand is also imposed on the environment; no cooling mechanism or cloud lifetime is demonstrated here.

The recorded latent heat comes from admitted deposition. McSnow's separate surface-melting heat accumulator does not evolve air or ice temperature in this caller, and the descent ends at freezing. Extending the calculation below that level would require a consistent melting and energy treatment.

The native unlimited falling cases reach the deliberate **1 m numerical stop** before leaving either cloud: at 0.477 h on Earth and 0.826 h on the Moon. Their retained partial trajectories diagnose amplification and do not predict metre-sized physical crystals. The altered cases reach the freezing level. With the finite allowance, native and altered lunar dimensions are nearly equal, 3.1373 and 3.1367 mm, because the water budget dominates their uptake difference.

Crossed controls separate depth from gravity while deliberately prescribing nonhydrostatic air. Keeping Earth gravity but using the deep cloud gives 3.367 mm in 15.627 h; keeping lunar gravity in the shallow cloud gives 3.129 mm in 8.970 h. Depth strongly changes residence time, while the fixed total allowance per incoming crystal largely fixes final mass in these supplied cases.

## Incoming number and seed history determine the size

The same 3 mm/h allowance can support larger individual crystals when fewer enter the column. The following lunar sensitivities retain the 4 µm sphere and the improved growth control:

| Incoming number flux | Final maximum dimension | Final mass | Time to freezing | Allowance deposited |
|---|---:|---:|---:|---:|
| 10,232.43 crystals/m²/s | 3.137 mm | 0.08127 mg | 54.177 h | 99.79% |
| 1,023.24 crystals/m²/s | 9.506 mm | 0.80804 mg | 31.789 h | 99.22% |
| 102.32 crystals/m²/s | 29.084 mm | 7.94064 mg | 22.510 h | 97.50% |

These fluxes are specified inputs, not inferred lunar nucleation rates or percentiles of one snow population. Larger crystals settle faster, reducing their residence time while remaining able to consume most of the larger allocation per crystal.

Seed history also changes geometry at comparable final mass. Starting instead from a compact **0.2 mm sphere**, with initial mass 0.0038495 mg, gives a supplied lunar crystal **0.744 mm** wide and **0.403 mm** thick after **11.458 h**, with mass **0.08363 mg**. The 4 µm seed finishes at **3.137 mm** width and **0.0284 mm** thickness, with mass 0.08127 mg. Their similar masses therefore conceal very different shapes and fall speeds. With unlimited supply the larger seed reaches 4.610 mm and 2.893 mg in 6.902 h. These controls expose sensitivity to initialization and native shape history, rather than selecting a physical incoming habit.

The unlimited small-seed lunar crystal is only **0.0352 mm thick** at its 147.47 mm span. Mechanical survival of that very thin predicted geometry has not been tested. Aggregation, fragmentation, liquid collection, turbulent transport, a distribution of seeds and a time-dependent cloud would change the population problem. The existing Earth mass-and-axis errors remain relevant to these conditional dimensions.

## Numerical evidence and reproduction

The matrix contains 32 fixed-air runs, two Earth parity runs, 28 principal/crossed/seed descents, four source-flux sensitivities and two no-growth settling controls. All 68 executions return successfully and pass raw-output validation: 32 reach their requested duration, 30 reach freezing, and six stop at the declared diameter censor. The latter comprise one fixed-air and two falling native configurations, each retained at both refinement levels. No endpoint is filled beyond a censored trajectory.

Both Earth parity controls reproduce all compared physical fields exactly at common saved times from the preceding native and altered experiments. No-growth, constant-air settling agrees with the independent exit time `H/v` on both worlds and preserves mass exactly. The falling trajectories' maximum normalized deposition-balance residual is 2.22 × 10⁻¹⁶. The independent integrated allowance checks and instantaneous source-gate checks pass.

Fixed-air timesteps are 1/32 and 1/128 s. Falling controls compare maximum steps of 2 and 0.5 s together with mass-increment limits of 1% and 0.5%; output and exit boundaries further shorten steps. Across all 32 refinement pairs, the largest measured endpoint change is **0.333%**, including native cases and censored endpoints. For the improved control it is **0.054%**. These are numerical refinement measurements, separate from empirical accuracy and uncertainty in the prescribed environment.

Two near-freezing final steps initially failed a relative-only comparison between requested rate and the rate reconstructed by subtracting rounded masses. Their differences were below half a mass-storage unit divided by the step length. The retained validation now includes that explicit floating-point resolution, with a regression test that still rejects changes exceeding it; no simulation output was altered to pass the check.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -B -m atmosphere.snow.mcsnow_lunar_growth run
python -B -m atmosphere.snow.mcsnow_lunar_growth assemble
```

The earlier pinned McSnow build and benchmark traces are prerequisites; the runner verifies them and does not substitute a fresh source revision. All new source copies, executables, namelists, raw CSVs and logs live under `research/runs/snow_scale/mcsnow_lunar_growth/`, through the existing link to the project-space drive. The product records source, compiler, tables, builds, run inputs and raw-record identities. The original producers and results remain unchanged.

The [check ledger](results/mcsnow_lunar_growth_checks.json) records 35 passing focused tests and a passing repository-wide `make check`: 3,677 Python tests passed, 20 skipped, all 191 JavaScript tests passed, and 1,314 Python files checked with zero layer violations. Provenance and ensemble integrity also pass. Raw reassembly, source reading and substantive report reading are recorded separately from those software checks. The present result establishes a controlled lunar extrapolation of the improved component. A finite-water population with justified incoming numbers and habits remains the next scientific extension before treating these dimensions as a forecast of a lunar snowstorm.
