# Terrestrial benchmarks for the snow calculation

The **Welss et al. (2024) deposition/sedimentation control** has been compared with independent characteristics; [welss_column.md](welss_column.md) records equations, results and controls. A separate [declared-source author aggregation rerun](welss_aggregation.md) and [component parity checks](welss_aggregation_physics.md) now extend that comparison. An independent coupled aggregation solver and the archived executable identity remain outstanding.

| Comparison | Evidence supplied |
|---|---|
| Python against pinned Fortran routines | Translation parity for tested inputs |
| Hydrodynamics against archived terrestrial outputs | Reproduction of author calculations under matched conditions |
| Local column under Earth gravity | Internal environmental control |
| Matched published population experiment | Conditional agreement between population solvers |
| Matched laboratory or field measurements | Empirical constraints within measured conditions |

An Earth-gravity swap does not reproduce published forcing. Lunar [geometry](geometry_literature.md) and [strength](strength_literature.md) constraints remain separate.

## Targets and recovered inputs

[Welss et al. (2024)](https://doi.org/10.1029/2023MS003805), Figures 17/19, supplies growth-and-sedimentation controls; Figure 21 adds aggregation. In the [author dataset](https://doi.org/10.5281/zenodo.10077724), `McSnow_paper_plotscripts.zip` contains `f19.ncl` and `f21.ncl`. These select `ob_hp0.nc`/`pr_hp0.nc` and `ob_agg0.nc`/`pr_agg0.nc` from `McSnow_paper_plot_data.zip`.

These inputs follow the [oblate](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/run/1d_obl) and [prolate](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/run/1d_pro) run scripts at [McSnow revision 1e428d7e901cbc5b91f2040e575f736ec5963a74](https://gitlab.dkrz.de/mcsnow/mcsnow/-/tree/1e428d7e901cbc5b91f2040e575f736ec5963a74). Run identity is qualified below.

| Setting | Oblate/prolate-labelled diagnostic cases |
|---|---|
| Computational domain | Both 0–5,000 m, 250 cells; area 5 m²; 20 m layers |
| Nucleation interval | 2,700–3,000 m / 4,500–5,000 m |
| Number production | 0.01 superparticles m⁻³ s⁻¹ × multiplicity 1,650 = **16.5 particles m⁻³ s⁻¹** |
| Five-second source update | 75 / 125 superparticles; floor operation introduces no rounding for these inputs |
| Integrated number source | 4,950 / 8,250 particles m⁻² s⁻¹ |
| Seed distribution | Gamma shape 2, mean mass **4.8 × 10⁻¹³ kg**; uniform nucleation height |
| Atmosphere | 273.15 K surface; 0.0062 K/m lapse rate; 101,325 Pa; g = 9.81 m/s² |
| Source humidity | `init_atmo_ideal`: 5% over ice above 600 / 1,000 m; water saturation below, followed by `get_atmo` interpolation |
| Process settings | `habit0`, `rt0`, `vt3`, `ckt0`; `agg0` then `agg2`; no boundary particle source |
| Duration | 36,000 s; five-second source/collision step; 600 s profile outputs |

[`nucl_simple`](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/src/mo_nucleation.f90) gives the source normalization; it is volumetric, despite an ambiguous script comment. Its surface number flux agrees within approximately 1.5% with the archive.

## Closures and comparison quantities

Match [`unrimed_geometry`/`m2d_fillin`](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/src/mo_mass2diam.f90) for geometry, [`boehmVel`](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/src/mo_velocity.f90) for velocity, and [`vapourDiffusion`, `ventilation_coef_i` and `capacitance`](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/src/mo_depo.f90) for deposition. Large-particle source relations are `m = 0.0443770 D^2.1` and `A = 0.131488 D^1.88` in SI, with the source's solid-ice transition and area floor. Capacitance distinguishes monomers and aggregates. Tested compact seeds select `C/D=1/π` through a positive rounding residual in void volume; nonpositive void selects `0.5`. Aggregation therefore requires monomer/aggregate identity even without habit prediction.

For aggregation, match [`hydrodynamic_kernel`](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/src/mo_coll_utils.f90), [`Ecbohm_spheroid` and `stick_effi_tempOnly`](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/src/mo_colleffi.f90) option 1. Use the source Mitchell sticking intervals; exclude local chipping, the lunar water-supply cap and its boundary seeds.

Compare `Fm` [kg m⁻² s⁻¹], `Fn` [m⁻² s⁻¹], `Md` [kg m⁻³], `Nd` [m⁻³], and `Fm/Md` [m/s]. **Archived `dep_tot` is kg/s per cell**: `sp_depo` sums multiplicity-weighted mass increments, and the [writer](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/src/mo_output.f90) applies no volume normalization. Divide by 100 m³. Summing `dep_tot` and dividing by 5 m² matches surface mass flux within 0.3–2.4%; the plot labels omit the volume conversion.

`Fm_rimed` is a residual after subtracting separately rounded category outputs; [preprocessing](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/postprocessing/NCL/hprofiles.ncl) masks nonpositive residuals before averaging. Its small positive values do not establish riming. `rt0` disables that process even when names contain `lwc2`. Preserve missing-value masks and compare total profiles first.

## Recovered identity and remaining limits

The [data archive](https://zenodo.org/api/records/10077724/files/McSnow_paper_plot_data.zip/content) contains `atmo_lwc.dat`, which `fig15.py` in the [script archive](https://zenodo.org/api/records/10077724/files/McSnow_paper_plotscripts.zip/content) reads for the background-atmosphere figure. Above 1.1 km, its ice supersaturation is 0.04999864–0.05000544; water-relative humidity falls to 78.0% at 5 km. This archived diagnostic supports the [public source's](https://gitlab.dkrz.de/mcsnow/mcsnow/-/blob/1e428d7e901cbc5b91f2040e575f736ec5963a74/src/mo_atmo.f90) **5% over ice**, despite the paper's water-reference wording. It does not identify each population run's exact atmosphere.

All four target `expdir` attributes contain `_vent`. Together with the publication and enabled source terms, this supports **base Hall–Pruppacher ventilation enabled** for the declared comparison. Table 1 distinguishes habit-dependent ventilation; diagnostic `phi=1` removes those additional shape terms. An older source comment does not establish disabled ventilation here.

The profiles were processed on 28 February 2023; public history starts at `1afcc5b51b7d172a9b67deee5879e5f6605a264f` on 26 March 2024. Atmospheric interpolation subsequently changed. Neither public version is identified as their producer. Public preprocessing averages four final snapshots, spanning 30 minutes at the stated cadence, while the paper specifies five hours. The archived profiles omit time axes, averaging windows and source hashes. Their finite endpoint bins exclude the examined full `opt=0` smoothing pass, and all 150 number-flux bins fit a 31-snapshot three-significant-digit mean. This supports direct raw-profile comparison and is consistent with five hours at ten-minute cadence; multiples of 31 also fit, so exact times remain unidentified. These are run-provenance limits, distinct from physical calibration.

Inspection covered publication sections 4.2.1 and 5.1–5.3, archive versions, figure scripts, four profiles, the atmosphere diagnostic, relevant source routines and public history. Extracted data-member lengths/CRC32 were verified and SHA256 values recorded below. The **full data-archive checksum remains unverified**. The complete script archive matched its published MD5, `a678b10d131f6a4059e7802a500a706b`. External bytes remain outside Git under `/tmp/lunar_snow_earth_population` and `/tmp/lunar_snow_bohm_sources`. The [raw comparison](results/welss_archive_comparison.json) gives full-profile differences of 0.55–1.77% and nucleation-layer number-density error of 7.62%. A separate, uncalibrated `C=D/2` compact-seed sensitivity reduces the latter to 3.14%, leaving substantial local differences. The subsequent author-source aggregation comparison is recorded separately in [welss_aggregation.md](welss_aggregation.md). Independent coupled aggregation and empirical validation remain unperformed.

```text
member          bytes   SHA256
ob_hp0.nc       40036   119871a2d269a90321a7fba21658b6f217aa2f8641dbda37f6cd9ae6ec41a304
pr_hp0.nc       40040   95d10776317a48eae6fa6936da46f4b3bcbed4618bf0e55529572ddb22821cca
ob_agg0.nc      40036   3766649cc9106aba999a78c603637dedf03b75e21f8bf21366fb652081bb3515
pr_agg0.nc      40040   81dc3826458c627811c464216c2657a340df59e612301087f14239c95d1444aa
atmo_lwc.dat   193193   6dbeca21e722b8d7c11d094446865235e04578e155417a125c7c428b344882cd
```
