# SCALE snow: known limitations and discriminating controls

The small absolute flakes in the prescribed SCALE columns cannot yet be treated as a physical lunar limit. Five new full-column controls show that an experimentally more open contact geometry increases lunar dimensions substantially. A separate comparison with an Earth single-crystal growth measurement exposes a larger problem: at 30 minutes, the observed fitted dimension is 4.17 mm, while native SCALE produces 27.66 mm from the fitted thin initial crystal and 0.96 mm from a sphere with the same initial mass and span. The growth treatment's dependence on initial shape must therefore be resolved before a geometry change can become a predictive repair.

These results preserve the preceding finding that the matched full-column model produces larger lunar flakes than Earth flakes. They add evidence about why its absolute dimensions are unreliable, while keeping the observed single-crystal benchmark separate from aggregate and storm validation.

## What the literature and active code identify

SCALE's authors explicitly discuss limitations of the filling-in aggregation assumption and the Chen–Lamb deposition-density treatment. They propose contact separation as a way to extend aggregate span, and identify remaining questions in habit growth, capacitance, collection efficiency and transport coefficients. Their reported Earth evaluation supports the framework while leaving these processes open to improvement. [Shima et al. (2020), sections 7.2 and 9.3](https://gmd.copernicus.org/articles/13/4107/2020/).

The [active-source audit](scale_size_sources.json) identifies how those assumptions enter our executable. The following distinctions determine which controls can test them.

| Process or input | Active implementation and consequence | Evidence needed for a repair |
|---|---|---|
| Aggregate geometry | The outcome routine recomputes a porous spheroid from combined mass and an assumed arrangement. Its arguments contain neither gravity nor collision energy or contact strength. It can produce compact aggregates, but does not impose a millimetre diameter ceiling. | Measured component arrangement, mass, span, area and fall speed across aggregate sizes; contact mechanics for extrapolation. |
| Individual-crystal growth | Temperature selects an intrinsic growth ratio; deposition density and spheroid capacitance determine mass and axis growth. Terminal speed enters ventilation. | A comparison of mass and dimensions with measured crystal growth at the same temperature, pressure, humidity and age. |
| Injected population | Our driver supplies 0.1–0.4 mm spheres with the inherited mass–size law and a fixed seed number flux. Their initial habit is an adapter choice. | Observed Earth seed histories first, then a defensible lunar nucleation and growth history. |
| Water and weather | The prescribed driver holds air at water saturation, shares a local deposition allowance among particles, and omits atmospheric motion and liquid particles. | Coupled moisture and motion, measured Earth forcing, and independent checks of riming and turbulence. |
| Other native approximations | Ice collection efficiency remains 0.1, diffusion and heat transport use fixed reference coefficients, and this release supplies no ice-fracture law. | Separate process controls. Missing breakup cannot be the active mechanism imposing the present small size. |

Early-flake observations support a range of contact arrangements. Kajikawa found a primary separation-ratio peak near 0.2 for same-shape planar pairs and a peak near 0.65 for different-shape pairs. Large flakes were excluded. Consequently, using 0.5 everywhere below is an experimental perturbation, not a fit to a universal measured geometry. [Kajikawa (1985)](https://www.cambridge.org/core/journals/annals-of-glaciology/article/structure-and-falling-motion-of-early-snowflakes/D3CA015917FB1A33E820F3C43C27BAB6).

Airflow also affects crystal habit. Experiments including short reduced-gravity intervals found sensitivity of faceted versus skeletal growth to the surrounding flow. This gives a reason to examine ventilation when transferring a growth law to lunar gravity; the abstract provides no lunar aggregate-size limit. Full text remains unretrieved. [Alena, Hallett and Saunders (1990)](https://www.sciencedirect.com/science/article/pii/002202489090157G).

The [source register](scale_limitations_sources.json) records reading extent, retrieved hashes, the original laboratory paper and its corrigendum, and failed retrievals. The earlier [size diagnosis](scale_size_results.md) remains the reference for controlled gravity, depth, seed-flux and water-supply changes.

## Offset planar contacts

The [experimental adapter](scale_separation.py) changes only collisions in which both incoming particles are oblate. Let their equatorial radii be a1 and a2, with full planar spans D1 = 2a1 and D2 = 2a2. For a prescribed separation ratio s, the new equatorial radius is

```text
a_new = max(a1, a2, (1+s)(a1+a2)/2).
```

This is the span construction suggested in Shima et al.'s aggregation discussion. We extend the native density interpolation to that radius: combined mass and the sum of constituent volumes set the lower-volume endpoint, while the sum of the incoming polar radii sets the other endpoint. The resulting density is then calculated from the conserved mass and final spheroid volume. That interpolation extension is our sensitivity closure; no contact-energy model or empirical calibration is implied.

Two equal 1 mm planar particles therefore give a 1.5 mm planar span at s = 0.5. The zero-separation wrapper calls the preserved native routine directly. Columnar and mixed incoming pairs also retain the native outcome. For multiple physical encounters in a sampled pair, the offset branch applies successive binary outcomes; counters report how often that occurs so this numerical distinction remains visible.

All column controls retain the original 132.59-hour duration, with arrivals gathered over its latter half. Eight isolated horizontal columns receive the same water and seed supply. Lunar controls use the 16.28 km column; the Earth controls shorten its depth in inverse proportion to gravity while preserving temperature and pressure versus fractional height. Thus these are matched model experiments, with the same forcing limits described in [the column methods](scale_lunar_methods.md).

| Case | Step, s | Mass D50, mm | Mass D99, mm | Snow water flux, mm/h |
|---|---:|---:|---:|---:|
| Moon, native | 60 | 1.590 | 3.114 | 2.762 |
| Moon, offset | 60 | 2.687 | 7.746 | 2.930 |
| Moon, offset, finer step | 15 | 2.742 | 7.816 | 2.930 |
| Earth, native | 15 | 0.742 | 1.025 | 1.822 |
| Earth, offset | 15 | 0.746 | 1.428 | 1.869 |
| Earth, native, finer step | 5 | 0.746 | 1.033 | 1.832 |
| Earth, offset, finer step | 5 | 0.750 | 1.440 | 1.880 |

Dimensions are maximum spheroid dimensions weighted by arriving ice mass; D99 is a percentile, not the largest flake. At the matched 60 s step, offset contacts raise the lunar median by 69% and D99 by 149%. In the finer offset run, the mass-median apparent density is 326 kg/m³ and the polar/equatorial axis ratio is 0.249, compared with 470 kg/m³ and 1.044 in the native baseline. Mean arriving particle mass decreases from 0.582 to 0.485 mg as the particles occupy a larger span. Both cases conserve their supplied water.

The Earth median changes by only 0.46% at 5 s, while its D99 increases by 39.4%. About 21% of arriving Earth ice mass is in aggregates, compared with about 90% in the offset lunar case. Changing how crystals attach therefore affects the lunar median much more strongly; the Earth control's mostly unaggregated population continues to depend on its single-crystal growth and source assumptions.

The lunar 60→15 s refinement passes the declared screen: bulk changes are at most 2.06%, and the maximum difference between exact mass-weighted cumulative distributions is 1.25 percentage points. The Earth 15→5 s comparison changes bulk statistics by less than 0.87%, but its distribution difference is 6.42 percentage points, exceeding the five-point criterion. This Earth distribution remains numerically unresolved.

The zero-separation run reproduces all arrival bytes, water and collision ledgers, and 81 native history arrays exactly. No collision omissions were recorded in the new runs; the largest relative water residual is 1.02 × 10⁻¹². The coarse lunar offset run has two sampled pairs containing multiple events; all finer runs have zero, so the refined result does not depend on the new sequential treatment of those pairs. These counters count sampled numerical operations, not multiplicity-weighted physical encounters. Particle-count, vertical-resolution and additional random-seed refinements of the new geometry remain unperformed.

## Single-crystal growth against an Earth measurement

The [growth runner](scale_crystal_growth.py) uses the −14.4 °C row of Table 2 in Takahashi et al. The original measurements follow freely suspended crystals at water saturation and approximately 1010 hPa. Their empirical fit gives a 4.175 mm dendrite at 30 minutes. The table's axis lengths are full dimensions in millimetres and its mass is in grams; the comparison converts these explicitly to native radii and kilograms. [Takahashi et al. (1991)](https://www.jstage.jst.go.jp/article/jmsj1965/69/1/69_1_15/_article/-char/en).

The simulation begins at the three-minute fitted mass and dimensions, avoiding an assumed nucleation size. It then advances 27 minutes at constant temperature, pressure and water-vapour pressure. Particle positions are held fixed as in a suspension experiment, while native terminal speed still controls ventilation. Eight identical isolated crystals have no collisions, liquid population or new seed injection. Their numerical water allowance greatly exceeds actual demand, so the column's shared supply cap cannot explain a discrepancy here.

The observed-shape initialization maps the measured mass and full axes to a porous spheroid. The spherical control preserves the same initial mass and maximum dimension while changing thickness and density. An additional lunar-gravity case preserves the laboratory air and initial crystal; it is a gravity sensitivity, not an observed lunar growth curve.

| Crystal at age 30 min | Step, s | Maximum dimension, mm | Mass, mg | Thickness, µm |
|---|---:|---:|---:|---:|
| Observed Earth fit | — | 4.175 | 0.0729 | 25.80 |
| Earth, fitted initial axes | 10 | 27.167 | 1.3470 | 12.59 |
| Earth, fitted initial axes | 5 | 27.440 | 1.3717 | 12.57 |
| Earth, fitted initial axes | 1 | 27.662 | 1.3919 | 12.55 |
| Earth, same-mass/same-span sphere | 5 | 0.956 | 0.0571 | 497.39 |
| Moon gravity, fitted initial axes | 5 | 13.907 | 0.3995 | 14.24 |

The refined Earth run overshoots the fitted width by a factor of 6.63 and mass by a factor of 19.1. Its thickness reaches only 49% of the fitted value. The native axis rule explicitly gives d ln(c)/d ln(a) = Γ*, with the intrinsic temperature-dependent ratio multiplied by a ventilation correction. Over the simulated interval, the integrated ratio is 0.112, compared with 0.509 calculated from the observed fits. This identifies excessive preference for lateral growth in this configured representation. Capacitance and ventilation then change with the growing span, coupling the shape discrepancy back into vapour uptake.

The sphere starts with the same 0.398 mm width and 0.000384 mg mass, but it is initially much thicker and has a lower apparent density. The same growth laws retain a comparatively blunt shape and yield 0.956 mm. Initial shape thus changes the result by a factor of about 29 at the matched 5 s step. This is a controlled sensitivity, not evidence that real crystals would retain the same dependence or that the fitted spheroid exactly represents a branched crystal.

At fixed elapsed time and the same laboratory air, lunar gravity gives 13.91 mm from the thin initial crystal versus Earth's 27.44 mm at the matched step. That is the model's ventilation response in this isolated experiment. The full column also changes fall time and collision exposure, so its lunar/Earth ordering need not match this fixed-age comparison. Because the Earth growth curve fails, the ratio here cannot be adopted as a lunar growth prediction.

The fitted planar thickness was inferred from crystal mass and basal area. It is therefore not an independent third measurement, and the fitted curves have no pointwise measurement uncertainties. This one-temperature comparison can locate a growth discrepancy and test sensitivity to initialization; it does not establish validity across habits, temperatures or aggregate populations. The [corrigendum](https://www.jstage.jst.go.jp/article/jmsj1965/69/2/69_2_251/_article/-char/en) changes Figure 3, not the fitted coefficients in Table 2.

The 10→5→1 s sequence produces endpoint width changes of 1.00% and 0.81%, and mass changes of 1.83% and 1.47%. Across the saved growth trajectories, the final refinement changes width by at most 1.18% and mass by 2.29%. These are measured sensitivities, far below the discrepancy with the fitted experiment; they do not establish an asymptotic convergence rate.

Native histories maintain 258.75 K and 101000 Pa. The imposed Murphy–Koop water saturation reads 100.0475% on SCALE's saturation convention, an explicitly retained difference of 0.0475 percentage point. No case uses more than 6.4 × 10⁻⁷ of its numerical deposition allowance. Ice growth and deposited mass agree to relative error below 6 × 10⁻¹⁵. The forcing hook's raw water residual equals the known initial crystal mass, because that hook originally assumed an empty column; the verifier adds that initial inventory explicitly. Residual subtraction involving the much larger vapour inventory loses precision, giving a maximum absolute corrected residual below 8.8 × 10⁻¹⁵ kg/m². Particle masses at matching timestamps independently agree with the ledger.

## Consequence for the lunar prediction

The present absolute sizes are controlled by particle representation and imposed growth histories as well as available water. Offset contacts demonstrate one contribution to compact aggregates. The single-crystal comparison shows that adopting that change alone would leave a substantial empirical failure upstream of aggregation: the same native growth treatment can give either submillimetre blunt crystals or oversized thin crystals, depending on the supplied initial shape.

The next repair target is the coupled habit, deposition-density, capacitance and ventilation treatment. Its terms should be diagnosed against the measured mass and axis histories, then alternatives tested separately across the dendritic and columnar laboratory cases. Starting-shape uncertainty belongs in that comparison. Agreement at one fitted temperature would require checks at other temperatures before altering the production column. Contact geometry and sticking can then be assessed against observed aggregate populations under measured Earth forcing, with the resulting uncertainties carried into lunar runs.

This work identifies a specific empirical growth failure and a separate geometry sensitivity. A unique corrected growth closure, observed Earth aggregate distribution, mechanical survival law, and calibrated lunar snow/optical forecast remain open. The earlier numerical controls retain their value as solver checks; they did not establish realistic Earth crystal or aggregate sizes.

## Reproduction and verification

The [geometry product](results/scale_separation.json), [growth product](results/scale_crystal_growth.json) and [check record](results/scale_limitations_checks.json) preserve numerical outputs, source/build identity, inputs and verification. The baseline native source is unchanged; the geometry experiment lives in two separately identified builds under ignored `research/runs/snow_scale/`. At most two numerical or build workers ran concurrently, with one thread each.

```sh
python -B -m atmosphere.snow.scale_separation build --separation 0
python -B -m atmosphere.snow.scale_separation build --separation 0.5
python -B -m atmosphere.snow.scale_separation run --workers 2
# Run separately, in an environment containing NumPy and netCDF4:
python -B -m atmosphere.snow.scale_crystal_growth run
```

The growth runner requires the original publisher PDF at the hash and ignored cache location recorded in the source register. It does not silently substitute another source. The independent verifier regenerates all 50 physics input files for the ten new executions byte-for-byte and checks their recorded executable identities. Column arrivals agree with native accumulated precipitation within 2.0 × 10⁻¹⁰ relative error. The focused suite passes all 60 tests, including compiled conservation checks for 640 oblate pair outcomes, native parity for other encounter types, and the laboratory unit conversion. Layer and provenance results are recorded separately in the check record. The subsequent pre-commit `make check` also passes for this experimental work: 3251 Python tests passed, 20 skipped, 191 JavaScript tests passed, 1253 layer files checked with no violations, and provenance and ensemble integrity checks passed. These numerical and integrity checks do not supply the remaining empirical validation.
