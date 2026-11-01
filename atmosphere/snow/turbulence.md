# Motions and growth environments in the aerial winter zone

The archived cloud simulations distinguish relatively quiet falling snow from the more active air where snow
coexists with supercooled cloud water. Horizontal winds remain moderate in both settings, but the mixed-phase
cells usually rise faster than the roughly 0.6 m/s settling speed of large aggregates. This gives a reason to
separate a snowy community's surroundings from the cloud producing its snow when estimating flake residence,
growth and mechanical loading.

[turbulence.py](turbulence.py) reads 690 existing CM1 snapshots and writes
[results/turbulence.json](results/turbulence.json). These are diagnostics of kilometre-scale simulations, without
empirical lunar validation. The cases use the older **A28_dim5 / A28_dim5_moon forcing**, whereas the separate
flake-scale examples use the **301 K resort target's winter column**. The two calculations establish distinct
pieces of evidence; they do not form a coupled prediction of resort weather.

| Archived case | Grid and sampling | Model-time interval |
|---|---|---|
| `box_0e_small` — Morrison microphysics | 32 × 32 columns; 237 outputs every 3 hours | Days 29.5–59.0 |
| `box_0e_elec_corrected` — corrected electrified box | 64 × 64 columns; 237 outputs every 3 hours | Days 29.5–59.0 |
| `box_0e_elec_corrected_storms_second` — selected active window | 64 × 64 columns; 216 outputs every 10 minutes | Days 40.5069–42.0 |

All three have 6.01 km horizontal spacing and approximately 2 km layers in the sampled cold air. The diagnostic
selects cells at 20–60 km whose actual temperature is −20 to 0 °C. Its main snow threshold is
`qs ≥ 10⁻⁵ kg/kg` (0.01 g of snow per kg of dry air); mixed-phase cells additionally have
`qc ≥ 10⁻⁵ kg/kg` of cloud liquid. Temperature follows CM1's Exner relation. Vapour pressure comes from the
stored pressure and vapour mixing ratio, and relative humidity over ice uses the Murphy–Koop saturation formula
already implemented in [flakes.py](flakes.py).

The following percentages describe **equal-weight sampled cell-times**. The counts in parentheses are each
row's denominators; neighbouring cells and successive times are correlated. They are neither independent
observations nor probabilities of exposure along a flake's trajectory. The selected storm window also cannot
establish storm occurrence rates.

| Conditioning | Morrison: fraction with \|w\| > 0.6 m/s (n) | Electrified: same fraction (n) | Active window: same fraction (n) |
|---|---:|---:|---:|
| All cold air | 0.62% (2,314,058) | 0.34% (9,259,059) | 2.02% (7,656,283) |
| Snow ≥ 0.01 g/kg | 11.8% (83,631) | 13.5% (70,885) | 19.8% (313,499) |
| Snow and ice relative humidity ≥ 1 | 37.8% (9,717) | 44.0% (10,654) | 58.5% (44,611) |
| Snow and cloud liquid ≥ 0.01 g/kg | **88.2% (2,565)** | **89.3% (1,657)** | **94.8% (9,367)** |

The liquid-bearing snow cells all meet the ice-saturation criterion in this sample. Their median vertical winds
are upward at **3.24, 3.01 and 4.36 m/s**, respectively. They account for only **3.07%, 2.34% and 2.99%** of
snow-containing cell-times. Restricting the main snow mask to 0.1 g/kg raises the fraction with |w| > 0.6 m/s to
22.0%, 22.1% and 31.0%; the product also retains a 0.001 g/kg threshold. More concentrated snow therefore
occupies more vertically active air than the general cold-air sample.

For the two full-cycle boxes, horizontal winds in snow-containing cells have medians of **3.53 and 3.90 m/s**
and 99th percentiles of **6.13 and 6.65 m/s**. Median vertical shear of horizontal wind is
**2.76 × 10⁻⁴ and 3.27 × 10⁻⁴ s⁻¹**. A drifting flake shares the air's mean horizontal translation; these
wind speeds therefore do not specify its aerodynamic slip. The resolved vertical motion can nevertheless
change where a flake travels and how long it remains in a growth environment. Its actual trajectory requires
the evolving three-dimensional flow and a particle model.

Moisture further separates snowfall depth from growth depth. **75.8% and 64.5%** of snow-containing cells in
the full-cycle boxes have ice relative humidity below **95%**; the active window gives 69.8%. Snow in such air
can lose mass by sublimation. For every sampled horizontal column containing cold snow, the diagnostic also
finds the largest contiguous depth of each environment, keeping zero when that environment is absent:

| Maximum contiguous depth, km: median / 90th / 99th percentile | Morrison: 15,205 column-times | Electrified: 21,453 column-times |
|---|---:|---:|
| Snow ≥ 0.01 g/kg | 10 / 18 / 20 | 6 / 12 / 18 |
| Ice relative humidity ≥ 1 | 0 / 4 / 16 | 0 / 4 / 16 |
| Ice relative humidity ≥ 0.95 | 0 / 8 / 18 | 2 / 6 / 18 |
| Cloud liquid ≥ 0.01 g/kg | 0 / 0 / 12 | 0 / 0 / 8 |

Environment masks in this table need not themselves contain snow: they describe the cold air available in a
snowy column. Layer boundaries define the depths, so the 2 km grid cannot resolve thinner moist or dry layers.
The temperature interval is also wider than the aggregation model's −15 to 0 °C column. Even with that wider
interval, deep snow-bearing layers occur much more often than equally deep contiguous saturated layers.
The ideal model's persistent, liquid-saturated growth path remains a conditional experiment; these snapshots
do not measure its duration or the frequency with which flakes complete it.

Flake-scale velocity differences require a separate calculation because the saved winds do not resolve them.
The sensitivity cases prescribe turbulent energy dissipation ε from **10⁻⁵ to 10⁻² m²/s³**. Terrestrial Arctic
mixed-phase observations associate changes in aggregation and riming with ε moving from below 10⁻⁴ to above
10⁻³ m²/s³; the outer decades are our sensitivity choices, without assigned lunar probabilities.
[Chellini and Kneifel (2024)](https://doi.org/10.1029/2023GL106599).

For an approximately isotropic local cascade, the dimensional increment is `δu = (εD)^(1/3)` at separations
above the Kolmogorov length `η = (ν³/ε)^(1/4)`, matched to `δu = D√(ε/ν)` below it. The prefactor is fixed
to one, so these are characteristic scales rather than calibrated RMS values or extreme gusts. The inertial
expression also requires D well below the integral scale.
[Kolmogorov's similarity hypotheses](https://doi.org/10.1098/rspa.1991.0075).

The reference air is −10 °C at 40.10 km in the target column, with density 0.690 kg/m³ and kinematic viscosity
2.39 × 10⁻⁵ m²/s. The [corrected dendrite-area geometry](habits.md) gives a 60 cm aggregate a settling speed
of 0.614 m/s in that air. It remains an extrapolated geometry, not a measured giant snowflake.

| Prescribed ε, m²/s³ | δu across 15 cm, m/s | δu across 60 cm, m/s | Adverse aligned-slip pressure increment for 60 cm, Pa |
|---:|---:|---:|---:|
| 10⁻⁵ | 0.0114 | 0.0182 | 0.00782 |
| 10⁻⁴ | 0.0247 | 0.0391 | 0.0171 |
| 10⁻³ | 0.0531 | 0.0843 | 0.0382 |
| 10⁻² | 0.114 | 0.182 | 0.0884 |

The last column changes the settling slip from V to V + δu, giving
`Δq = ρVδu + ½ρ(δu)²`. The product separately reports the fluctuation-only pressure `½ρ(δu)²`; using that
quantity alone would omit the cross term around existing settling flow. Neither pressure specifies tensile
stress on an ice contact. Drag distribution, porosity, flutter, collisions and contact strength must determine
breakage, so the calculation supplies no survival threshold or maximum flake size. Turbulence can also increase
collisions and growth: Earth snowfall observations find aggregation and collisional fragmentation together.
[Kötsche et al. (2026)](https://acp.copernicus.org/articles/26/3277/2026/).

The missing dissipation field reflects the model configuration as well as its saved output. These CM1 cases use
`cm1setup=2`, `ipbl=2`, `sgsmodel=0`, `horizturb=0`, and disable TKE and diffusivity output. Inspection of
CM1 r22's `turb.F` identifies `ipbl=2` as the Bryan–Rotunno diagnostic vertical Smagorinsky closure.
`restart_write.F` stores LES TKE only when its LES/TKE flags are active, and stores MYNN/MYJ variables only for
`ipbl=4,5,6`. The selected cases therefore have no hidden prognostic cloud TKE field to recover from restarts.
They do include vertical mixing. Enabling `output_tke` alone would not establish the unresolved dissipation;
a subsequent calculation would need a justified turbulence closure or sufficiently resolved velocity statistics.

The product records producer hashes, named shared constants, the target-column input hash, and each case's
control-file hashes. For every snapshot it records the exact selected levels, field order and a SHA-256 hash
of the selected raw float32 fields, allowing the archived bytes to be checked without redistributing them.
[turbulence_sources.json](turbulence_sources.json) records source identity, reading extent, the CM1 source audit
and the checked 2026 corrigendum. Reproduction requires the existing external CM1 archives linked at
`climate/crm/runs`; it starts no climate simulation:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m atmosphere.snow.turbulence
OPENBLAS_NUM_THREADS=1 python -m pytest atmosphere/tests/test_snow_turbulence.py -q
```

The nine tests check conditioning, invariance of deformation under uniform translation, the viscous and inertial
scalings, geometry consistency with the legacy fall-law control, selected-field reading, contiguous-layer handling,
empty samples and product provenance. They test the calculation and its bookkeeping; the giant-aggregate
mechanics and the inferred experience of snowfall remain separate physical questions.

**Provenance, 25 October 2026.** [results/turbulence.json](results/turbulence.json) was built with the rain solver as it
stood before 25 October, whose Chapman–Enskog vapour diffusivity stood 17.6% below Massman's (1998) measured value. It
keeps its values; its record names `atmosphere/rain/drops.py` at that solver's hash, kept byte for byte as
[atmosphere/rain/baseline/drops_2026-10-24.py](../rain/baseline/drops_2026-10-24.py) and listed in
`research/provenance.json`, and its identity checks read that copy ([shared/baselines.py](../../shared/baselines.py)).
