# Full lunar-column comparison: SCALE and the custom solver

The subsequent [gravity/depth and source diagnosis](scale_size_results.md) finds that this SCALE configuration produces larger flakes under lunar conditions than in its matched Earth counterpart. Its small absolute sizes therefore require examining the shared population forcing and particle laws; the comparison below establishes no lunar-specific size ceiling.

The complete comparison produces millimetre snow with SCALE under the deep, continuously supplied conditions that produced decimetre snow in the custom model. **All ten SCALE executions and four additional custom controls complete the same 132.589-hour scenario.** SCALE's native particle, timestep and vertical refinements pass the declared practical screens. This is new evidence about the full prescribed column: the large-flake result depends strongly on the representation of particles and their interactions. It does not follow inevitably from lunar gravity, cloud depth and the water allowance.

The [methods](scale_lunar_methods.md), [runner](scale_lunar_comparison.py), [versioned product](results/scale_lunar_comparison.json) and [verification record](results/scale_lunar_checks.json) retain the exact inputs, source changes, executable identities and complete comparisons. Each SCALE execution contains eight horizontally isolated columns; the additional custom controls contain two realizations at each of two statistical volumes. The old event-column results remain preserved in their [original report](event_column_results.md).

## The full column changes the conclusion drawn from the component trials

Both frameworks receive the original 16.2805 km cloud from −15 °C to freezing, its pressure profile, 10232.43 injected crystals/m²/s, the same injected mass and maximum-dimension distribution, and the 3 mm/h local deposition allowance. Arrival statistics cover the final 66.2946 hours. The principal SCALE comparison maintains the original temperature and saturation, so its smaller particles cannot be attributed to allowing the atmosphere to warm or simply ending the experiment early.

| Complete-column calculation | Snowfall (mm/h water equivalent) | Mass-weighted median dimension (mm) | Mass-weighted 99th percentile (mm) |
|---|---:|---:|---:|
| Custom, ballistic/no chipping; 0.0004 m² | 2.716 ± 0.007 | 140.1 ± 2.3 | 225.8 ± 2.2 |
| Custom, ballistic/no chipping; 0.0016 m² | 2.743 ± 0.001 | 181.0 ± 1.7 | 317.9 ± 20.3 |
| SCALE, native ice efficiency 0.1 | 2.762 ± 0.001 | 1.593 ± 0.003 | 3.101 ± 0.023 |
| SCALE, original temperature-dependent efficiency | 2.181 ± 0.002 | 3.556 ± 0.011 | 7.467 ± 0.184 |
| SCALE, finite vapour reservoir and maintained temperature | 2.802 ± 0.002 | 1.263 ± 0.003 | 2.520 ± 0.031 |

Entries are means of column/realization quantiles ± one standard error, rather than quantiles of a pooled population. These errors describe numerical sampling under each model's assumptions. The custom areas specify statistical volumes, not walls around physically confined flakes. Native weighted-particle refinement is a different approximation from increasing that individual-particle volume.

The custom controls use ballistic settling and disable chipping, as SCALE does here. Their medians still rise from 14.0 to 18.1 cm as area, and therefore each layer's statistical volume, increases fourfold. These two changes therefore do not explain the contrast with SCALE. Both implementations carry similar water fluxes in their prescribed default cases: approximately 2.7–2.8 mm/h. The difference is chiefly in how that water is distributed among surviving particles and in the dimensions assigned to their mass.

With its native collection efficiency, SCALE's six prescribed-column configurations give mean medians of 1.580–1.593 mm and mean 99th percentiles of 3.101–3.136 mm. Their largest particle anywhere during a run is 4.93–5.56 mm. These are calculated distributions and sample maxima, not a physical maximum size for lunar snow.

## Refinement and the stronger sticking law

| Refinement | Maximum mass-CDF change (percentage points) | Flux change | Median change | 99th-percentile change | Screen |
|---|---:|---:|---:|---:|---|
| particle sampling, 2 to 4 packets/minute | 1.002 | +0.414% | +0.814% | -1.031% | Pass |
| particle sampling, 4 to 8 packets/minute | 0.436 | +0.271% | -0.124% | +0.467% | Pass |
| 60 to 30 s at fixed 8 packets/minute | 0.390 | -0.102% | +0.136% | +0.667% | Pass |
| 30 to 60 layers at fixed column depth and injection sampling | 0.663 | -0.534% | -0.217% | +0.307% | Pass |
| sticking sensitivity, 60 to 15 s at fixed 4 packets/minute | 0.789 | -0.384% | -0.153% | -0.729% | Pass |

The screen permits at most 5% changes in flux and the 50th, 90th and 99th percentiles, and at most five percentage points in the mass CDF. The retained 90th-percentile comparisons also pass. Changing the native case's seed changes the median by −0.484%, the 99th percentile by +0.316%, and the CDF by 0.631 percentage points. These checks support stability over the tested resolutions and ensemble; they supply no uniform error bound for arbitrarily rare extremes.

SCALE normally assigns ice collection efficiency 0.1, whereas the custom model uses `exp[0.09(T−273.15)]`, approximately 0.26–1 through this cloud. A separate build changes only that factor inside SCALE's existing ice collision kernel. It raises the mean median to about 3.56 mm and the 99th percentile to about 7.47 mm. Thus stronger sticking has a substantial effect, but does not recover the custom model's dominant decimetre population in these tests.

The initial 60 s stronger-sticking run clips 238 collision pairs, omitting 3,255,435 physical collision events, or 0.00866% of its attempted count. That small fraction alone cannot establish a small tail error. Reducing the step to 15 s, while holding four injected superparticles per physical minute, produces **zero omitted collisions** and passes the distribution screen. Its largest particle anywhere is 15.61 mm, compared with 21.73 mm in the coarse run; sample maxima remain more sensitive than the bulk quantiles. All native-efficiency cases also report zero omissions.

The retained morphology matters independently of the sticking factor. In the native four-packet case, the pooled mass-median apparent density is about 470 kg/m³. In the coarse stronger-sticking case it is about 312 kg/m³. Assigning that latter population's masses to the old Brown–Francis dimension law would give a mass-median dimension of 18.1 mm, compared with SCALE's actual 3.56 mm. The old geometry therefore magnifies the dimension of a given mass substantially. It still would not reproduce the custom column's approximately 140–181 mm median: the evolving particle masses and number distribution also differ. Under matched stronger sticking, SCALE retains roughly 619,000 arriving particles/m²/h, versus about 8,700–12,600 in the additional custom controls.

## Vapour and latent heat are separate effects

Allowing a vapour reservoir while maintaining temperature gives a mean median of 1.263 mm and a 99th percentile of 2.520 mm. Its initial cloud contains 48.341 kg/m² of vapour. The bounded source supplies another 354.734 kg/m² during the complete run; 356.642 kg/m² deposits on ice, leaving 46.433 kg/m² of vapour. These quantities distinguish a reservoir and source ceiling from the original instantaneous deposition ceiling. The smaller flakes persist under both choices.

The separate thermal control retains native latent heating and melting, with no external cooling or air dynamics. Its initially cold column finishes between 273.065 and 273.223 K, close to freezing throughout. All precipitation arriving during the gathering interval is liquid; the mean arriving-drop median is approximately 0.578 mm, with only 0.0633 mm/h of precipitation. This is a diagnostic response of an uncooled, motionless atmosphere. Maintaining a cold cloud for five days requires a heat-removal mechanism, and this experiment omits radiation and dynamical cooling. It does not predict that the actual lunar winter zone disappears.

The ledger accounts for initial vapour, injected ice, external vapour or prescribed deposition, suspended water and precipitation. The largest relative residual across the ten executions is 1.93 × 10⁻¹¹. An independent comparison against native SCALE precipitation histories differs from the arrival records by at most 7.31 × 10⁻⁸ kg/m², or 2.99 × 10⁻¹⁰ relative. Every prescribed-deposition total remains within the original 397.768 kg/m² allowance, plus separately accounted injected seed ice.

## What this resolves, and what remains distinct

SCALE now supplies an executed full-column alternative with stable millimetre-scale statistics for this imposed scenario. The original large-flake result is not robust to replacing the old particle treatment. Matching settling and chipping choices does not remove the discrepancy; matching sticking increases SCALE's sizes, but leaves a large difference. Atmospheric vapour and temperature feedback can further change the answer, yet neither feedback is needed to obtain the principal contrast.

This comparison changes the numerical representation together with particle geometry, projected area, drag, capacitance and growth laws. It therefore does not uniquely identify the cause of the original solver's failed volume limit. That earlier sparse-population and water-allocation diagnosis remains valid for its own equations. The remaining distinction lies in how the two systems aggregate particles, grow their mass, evolve their geometry and represent their population; matching those laws is a separate diagnostic from the comparison completed here.

The particle forms in SCALE are approximate porous spheroids, not resolved crystal branches. No contact-fracture calculation was added. A matched published Earth snow-population reproduction, justified lunar moisture and cooling histories, giant-aggregate mechanical survival and predictive storm optics remain open. The reported sizes support a working model comparison, not an adopted empirical forecast. Earlier snowstorm optical calculations retain their explicit dependence on the unvalidated population used to produce them.

## Implementation verification

`make check` passes: 3,226 Python tests with 20 skips, 191 JavaScript tests, layer enforcement, research provenance/numerical checks and the ensemble validator. The new focused controls cover the original forcing identity, conservative arrival weighting, phase separation and the isolated sticking modification. Independent native-history checks verify all ten runs, and the final adapter exactly reproduces all five native input files for each run using the original numerical environment. The verification record includes an additional layer scan of untracked sources, source hashes, retained setup failures and limits of the checks. No commit or empirical validation is implied by these results.
