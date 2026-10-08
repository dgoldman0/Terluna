# Atmospheric electricity branch review for joint integration

Review date: 2026-10-08. Reviewed branch `research/atmospheric-electricity-plan`, exact head `24d236b03c760171de8bbd98282c8902cd72b158`; merge base with current `main` is `88bb4f5e967804f6c0da01ba0ad2c89d95453516`. Main was `2f2e0a1`; comparison branches were infrastructure `625f60e`, solar shield `eeaecbc`, sea appearance `46681f4`. This is an analysis artifact only. No branch was merged, rebased, committed, or changed; no simulation or test suite was run for this review.

## Assessment

This branch contributes a substantial new domain capability and two cross-domain studies. Stages 1 and 2 of atmospheric electricity are complete at their declared screening/model scope; the global circuit, conducting upper boundary and transient luminous events remain stage 3 proposals. The cap is `d507275`; stage 2 synthesis is `dd38e10`; `24d236b` only dates the status record. It is suitable to retain in the unified research tree, with historical products preserved, explicit current-result pointers, provenance repairs, and current indexes reconciled. Its physical claims should retain their one-site, coarse-grid and empirical-transfer limits.

The branch's incremental diff contains 78 files, 152,761 insertions and 15 deletions. Most of the volume is diagnostic JSON. Its source/code does not directly overlap the other branches' changed code, but the shared research indexes and CM1 README do overlap. Scientific integration has more work than textual conflict resolution.

Authoritative current narratives: `research/studies/atmospheric_electricity/README.md` (“What stage 2 found”, “Where the study stands”); `climate/crm/README.md` (“The main run”, “Storm lives”, “Threads and reproducibility”); `atmosphere/electricity/README.md`; `research/studies/open_moon_aerosol/README.md`.

## What it supplies

| Component | New contribution | Evidence boundary |
| --- | --- | --- |
| `atmosphere/electricity/charging.py` | Collision-integrated graupel/ice/snow charging under four laboratory laws | Laboratory extrapolation to selected particle populations; low-speed measurements do not constrain all NSSL hail collisions |
| `muons.py`, `conductivity.py` | Deep-column cosmic-ray/muon ion production; ion balance; humidity-grown aerosol and conductivity by height | MCEq/CRAC/PDG and Earth checks; surface radioactivity omitted; assumed aerosol aloft; no conducting ionosphere |
| `climate/crm/cm1_elec.py`, `fortran/*`, `cm1_run.py` | Pinned WRF-ELEC/NSSL coupling, charge transport, Poisson field solve, branched/cylindrical discharges, leakage, optional screening/corona, ground leader crossing; restart/fine/sample machinery | Earth-calibrated discharge rules; periodic flat domain; terrain under field solver unsupported; CM1 multithread race unresolved |
| `mixed_phase_analysis.py`, `elec_analysis.py`, `storm_lives.py` | Particle/charging diagnostics, flash statistics and charge structure, ten-minute storm tracking | Output cadence and core threshold matter; rerun windows are fresh weather realizations |
| Atmospheric-electricity study | Stage-1 bound; corrected two-lunar-day simulation synthesis; thunder and nitrogen estimates | One equatorial lowland site at 6-km grid, with lunar outcomes calculated rather than measured |
| Open-Moon aerosol study | Regional day/night particles, fog conductivity and CCN built from ecology choices, GCM geography and CM1 weather | Ground-only structured Earth-analogue estimates; not a transported aerosol or global cloud model |

The implementation preserves a useful distinction between storm microphysics, electric discharge, acoustic propagation and nitrogen-yield inference. Keeping these separable avoids treating one benchmark as validation of all four.

## Current quantitative state

The canonical stage-2 run is `box_0e_elec_corrected`, analyzed in `climate/results/crm/elec_box_0e_elec_corrected.json`, with ordinary weather in `box_box_0e_elec_corrected.json`. It restarted from day zero after the ventilation repair (`1fa7f41`) and ran on October 7–8 for 27.2 wall hours on six threads with one executable throughout. It includes the author-adopted lunar rules from its start.

| Quantity | Current result | Qualification |
| --- | --- | --- |
| Lightning frequency | 530 flashes in two lunar days, on 12 model days; 0.022 flashes/km²/year | One 148,000-km² equatorial lowland box, none at night; cannot be multiplied into a global rate without occurrence weights |
| Ground strikes | 106 negative strikes, 20%; no positive strike | Saunders–Peck law plus Earth-calibrated inception and adopted 1-kV/m leader rule; not a universal fraction |
| Cloud flashes | Median 35.8-km start, 108 C, 86 GJ | 6-km horizontal grid, median channel height 29.8–41.8 km |
| Ground strikes | Median 149 C and 66 GJ; maximum 1,025 C, 325 GJ | Crossing model omits continuing current, cloud supply during descent, leader end emergence and current cutoff |
| Charging | Hail supplies 87%; main negative charge at 40–42 km, 90–96% on snow | Saunders–Peck normal tripole; Takahashi inverts polarity with similar flash rate |
| Storm lives | 14 of 139 tracked cores flash, first flash 0.7–1.5 h after core formation; flashing lives 5–12.2 h | Two ten-minute windows; their 68+125 flashes are not additional samples to add to main run's 530 |
| Leakage | Nighttime residual charge under 35 C; no terminal net box charge | Replaces the original nonconducting-ion artifact, not a global circuit solution |
| Near-ground field | Median maximum 4.9 kV/m within six hours of flashes; maximum 61.5 kV/m; 0.95 kV/m median where no particles reach lowest level | Transient box fields, not fair-weather climatology |
| Thunder | Cloud 60–68 dBA below; strikes 93–97 dBA nearby; 23–27 Hz cloud and about 53 Hz strike peaks | Acoustic share ×0.1–10 bracket, layered mean air, hard ground, no turbulence/soft ground/cold pools; 40–50-km daytime / 100–150-km quiet range is an analogue correction, not measured performance |
| Fixed nitrogen | Central routes give 0.02, 0.07 or 4.6 kg N/km²/year, about 0.2%, 0.7% or 47% of Earth's lightning per area | Per-channel-length, per-flash and per-joule routes diverge >100-fold; their full uncertainty intervals extend further; oxygen scaling is unmeasured |

The corrected NSSL box is 1°C warmer at 2 m than its Morrison counterpart, 0.8°C higher in dewpoint, rains 55% more by day and has about half the night fog. Neither resolves the older GCM/CM1 warmth and humidity disagreement. Existing climate pause and comfort range must stand.

Stage 1 remains useful as the published one-dimensional bound (`stage1.json`), with cloud conductivity around 10⁻¹⁶ S/m and roughly day-long charge relaxation. Its headline “slow graupel/law controls how often” is superseded for the NSSL storms by hail's dominance: stage 2 finds law choice mainly controls polarity. Retain the stage-1 output and its assumptions, not an undifferentiated final verdict.

At the ground, aerosol study central cases yield area-mean conductivity about 2.2–2.7 × 10⁻¹⁵ S/m and 55–67-minute relaxation, with fog conductivity 2–4% of clear air. Central CCN at 0.3% supersaturation average 110/cm³ by day and 53/cm³ by night. Values over the model's individual regions/cases vary much more; Earth checks lean high by up to two to three times. A surface aerosol number does not establish storm-level conductivity or observed visibility.

## Decisions actually adopted, results retained and models displaced

Preserve these author decisions as decisions, separately from calculated results:

- October 2: three-stage electricity program, storms already run before electrified CM1; ROCKE-3D remains a separate possible climate check (`75eb66c`, `a21a615`, decisions register).
- October 4: lift the 180-kV/m breakdown cap for lunar runs, preserving density scaling and 50-kV/m high-altitude floor (`62b4c1e`). Keep radiation's 140-µm ice cap to match Morrison, and log clipped ice (`0b9507d`). The cap does not change microphysics' own ice sizes.
- October 4: stop the 128-km-wide 2-km fine box after three days without a flash and collapsed deep convection (`0fc26fa`). This is neither evidence against lunar lightning nor a converged grid test.
- October 5: use WRF-ELEC's temperature-based ground rule (`da3253f`), leakage as lunar default (`4499e35`), and the 1-kV/m leader crossing as lunar default (`ac55dc2`). Point discharge remains an investigated optional setting, not an adopted default; the corrected run has `corona_v_m=0`, `screen=0`.
- October 5–6: hold regional aerosol for stage 3 (`184f320`), retaining the storm conductivity column's assumed aerosol and continuum attachment for now; put transition-regime attachment in stage 3's rebuild (`48b9af3`). Do not silently feed new aerosol into historical storm runs during integration.
- October 8: cap and await joint integration; stage 3 needs a later go-ahead. An integration commit would not itself authorize a new simulation campaign.

Displaced interpretations should remain explicitly historical:

- Zero ground flashes under a 5-km/matching-charge condition was a rule artifact; Earth benchmarks also initially had ground detection off. Do not let those counts override corrected ground-strike products.
- Kilocoulombs of persistent ions and roughly 20-kV/m domainwide surface field were artifacts of omitted conduction. Original point-discharge comparisons included that artifact; they do not establish the current fair-weather field.
- The original run ultimately counts 504 flashes in its tracked summary; old thunder analysis used 465. Current thunder is the corrected run's 530 with ground channels extended to ground. These are different result stages, not contradictory totals to average.
- Stage-1 200–500-fold charging-law contrast came from Morrison graupel. NSSL hail falls 4.8–7.7 m/s and strikes ice/snow at 4.4–7.6 m/s, supplies 82–88% in law windows, and moves the decisive laboratory question to −5 to −16°C and 0.02–0.2 g/m³ cloud water. The 1.2–1.8-m/s experiments bear on the smaller graupel share.
- WRF-ELEC's grid-point NOx yield is about fiftyfold low on the Earth benchmark and is rejected for the study's nitrogen estimate. Its stored diagnostic may remain, but is not a consumable nitrogen budget.
- Differences around 10–20% between threaded realizations do not establish a setting's effect; the exact source of the CM1 race is unpinned. Preserve that limitation rather than “repairing” it by smoothing the products.

## How the five states inform each other

| Other state | What electricity contributes to it | What it contributes to electricity / unresolved intersection |
| --- | --- | --- |
| Current main | Storm locations, charge, lightning, acoustic environment and aerosol/CCN provide missing physical conditions for climate, optics, ecology and habitation | Main's wave/tide fields replace generic marine roughness where sea spray is next developed; cloud-twilight tools can consume actual ice sizes. Main's existing climate uncertainty remains, and wave work is not an aerosol production/transport model |
| Infrastructure `625f60e` | Median flash starts at 35.8 km and fields 150–240 kV/m intersect the 35–45-km flight band and 35.6-km port crown. Charge/energy, hail, icing, 12–36-m/s storm updrafts and long thunder reach define design work | A 24-km grounded steel tower, hydrogen berths, long conductors and moving craft invalidate assumptions of unperturbed storm air. The cap review identifies this need but adopts no protective design. Need terrain-aware conductor-field/current studies, warning and closure rules; low average strike rate is not a safety claim |
| Solar shield `eeaecbc` | Stage 3 needs plasma conductivity, magnetospheric current closure, night-side persistence, storm current injection; these share the shield's plasma uncertainties | Shield branch changes upper-air heating/transmission, loss accounting, operating-radius/layout and magnetic alternatives. The electricity cap's four regional magnets / 1.5×10²¹ A·m² / ~10-radius standoff is an inherited scenario, not a fresh globally adopted field model. Circuit upper boundary requires agreed scenario-tagged UV/plasma inputs; neither branch solves it |
| Sea appearance `46681f4` | Ice larger than radiation cap and NSSL/Morrison fog spread constrain cloud light, dusk and night visibility; aerosol modes offer future extinction inputs | Sea appearance explicitly keeps microphysics' own sizes, with 140 µm a check. This resolves the choice for rendering existing cloud output, not the climate feedback caused by cap-thickened radiation. Its water optical choices and seas' chemistry cannot be inferred from aerosol salt/DMS analogues |

Particularly valuable shared questions:

1. **Storms, traffic and the port.** Can warning/route closure protect moving craft while the tower remains a conductor into storm charge? Does the tower trigger leaders without a natural flash? Charge fields exist, but the current field solver is flat/periodic and does not support terrain beneath it. This is new research, not a merge conflict to resolve by a constant.
2. **Microphysics, aerosol, optical climate.** NSSL halves fog compared with the Morrison rings supplying aerosol and some cloud scenes. A common scenario registry must identify scheme, site and output window; swapping products wholesale would erase a valuable model spread. Later aerosol transport and cloud activation could change both electrical leakage and visibility.
3. **Shield, plasma and circuit.** Above 149 km/22 hPa, there is no electron/conductivity model closing the circuit. Shield-filtered UV, cosmic rays, Lyman-α and charged particles enter together, with night-side and polar currents. The shield branch's new solar-history cases should be available as inputs, not treated as already solved ionization profiles.
4. **Lightning and ecology.** Nitrogen yield routes span over two orders of magnitude, at one site. Do not close a global nitrogen inventory or imply terrestrial-level biological supply; later include occurrence maps, NOx chemistry, deposition and losses.
5. **Waves and aerosol.** Main's resolved coastal waves/cycle now give a path to replace an analogue sea-spray estimate. Gravity-adjusted bubble production, breaking/whitecaps, sea chemistry and vertical mixing still require a model; the sea-appearance water-colour classes do not supply those quantities.
6. **Radiation beyond visible light.** Runaway-threshold fields imply a gamma-ray/glow question in the flight band; no dose calculation exists. Elves need peak currents, which the branched scheme does not output. Neither can be replaced by flash energy alone.

## Merge and provenance hazards

### Textual and semantic conflicts

Both sides changed `climate/crm/README.md`, `research/README.md`, `research/decisions.md` and `research/status.json`. Shield also changes `atmosphere/README.md`. Against their actual common ancestor, infrastructure and electricity overlap in `climate/gcm/README.md` and the three research-index paths; the pairwise conflicts are the GCM README, research README and status, while decisions auto-combines. Comparing only branch-authored changes against main misses inherited edits. These are merge-by-topic files: retain all domain/study entries and author decisions, not one complete branch's version.

Current cap prose is newer than several other statements in the same branch. `research/status.json`'s electricity `unresolved` still says the full two-day comparison/windows remain, despite their completion, and retains superseded benchmark/ground-rule language. Its condition string records chronological states such as “running” as well as “done”. `research/plan.md` still describes the original box as running from October 4; the long decisions row ends at “main run … start held” on October 6; stage-2 item 2 likewise retains the original run's start-time wording. Reconcile current-state pointers with the cap while preserving dated history. This is authorized cleanup for a future integration, not evidence that stage 3 ran or the general climate pause ended.

### Constants

The branch does not edit shared constants, but introduces physical literals already available or later centralized on main: `conductivity.py` defines Boltzmann, elementary charge and rounded universal gas constant locally; `muons.py` defines speed of light in cm/s and `RD=287.04` (CM1 closure); aerosol `inputs.py` uses literal 29.53 days. Vacuum permittivity is also local and would need a deliberate shared value/derivation if consolidated. The prospective unified code should use `shared/constants.json` and preserve the distinction between universal values and CM1's solver closures. Do not change historical scenario durations or byte-pinned third-party formulas casually. Any producer file edit changes hashes even if results are numerically identical: reconcile products honestly through reproduction or explicit historical linkage, not a header-only hash replacement.

### Products

Read-only producer audit performed for six current product families:

- Conductivity, stage 1, nitrogen, thunder and aerosol producer hashes match the files recorded in their headers at this branch head.
- `elec_box_0e_elec_corrected.json` records producer hash `6f158af5afeaa644`; that is exactly `climate/crm/elec_analysis.py` at `94082be`. Current code is `7eef6eaf777860d2`, changed in `496ab7f` to permit reusing an already-read snapshot. This is traceable historical provenance, not proof the result is wrong. Retain source commit or regenerate the diagnostic in an isolated output directory after confirming inputs; do not claim the historical product was made by current bytes.
- `stage1.json` lacks top-level evidence/reading-rule/input hashes; its generator likewise omits them. The narrative supplies assumptions but consumers should not have to infer the product contract.
- Conductivity records producer and evidence but no source-input hashes (mixed-phase profile, muon products and tables).
- Aerosol records climatology, atlas, drainage and conductivity hashes, but not the ring JSONs its regional weather reads or the imported conductivity/physics helpers. Nitrogen/thunder similarly need their raw-flash/profile/source-product dependencies pinned, beyond naming cases or outputs.
- Several study readers directly parse JSON without checking schemas. Plan a small dependency-manifest/schema-validation reconciliation rather than treating a same-name file as interchangeable across branches.

The atlas hash `78cfbf2f469eb339`, drainage `18494e3d3656ead7` and muon upper-profile CSV `e7079c911c9e70c4` were verified byte-identical across electricity, current main, shield and sea heads. Their age alone does not require a rerun. The solar branch's changed upper-air loss models do not automatically mean this particular source CSV changed. A rerun decision should follow actual consumed values and declared scenario, not every shared-file timestamp.

### Reproducibility outside Git

The electricity worktree contains ignored `atmosphere/electricity/inputs/`; all three input files were read and matched the manifest's exact hashes/sizes: Takahashi 6,070 bytes, CRII ZIP 712,048 bytes, PDG muon table 18,086 bytes. Restore through the pinned fetcher or copy verified bytes into the future worktree; never substitute inputs or commit outside tables by accident.

`climate/crm/runs` is an ignored symlink to `/media/projectspace/terluna-research/crm-runs`. Current canonical, original, fine, storm-window and Earth-supercell directories exist with `case.json` and `progress.json`; raw grids, logs, restarts, namelists and executable attribution are outside Git. `climate/gcm/products/climatology_A28_dim5_moon.npz` in this worktree is an ignored symlink to main's product. A fresh Git worktree does not recreate either. Avoid copying the entire raw archive; inventory exact required files and mount/link as read-only inputs where practical.

`TERLUNA_CM1_HOME` defaults to `/media/projectspace/terluna-research/cm1`. Pinned CM1 and WRF-ELEC sources, build records and executable files are present there; the MCEq venv is `/media/projectspace/terluna-research/venvs/mceq`. The branch's documented `climate/gcm/.venv/bin/python` does not exist locally. Select an existing interpreter by dependencies and record it rather than copying that command blindly.

The main corrected run and its windows link to the shared `cm1/build/moon_omp_elec/cm1.exe`; original/fine cases now point to a hash-named preserved executable. `cm1_run build moon_omp_elec` overwrites that shared build and can change what existing cases run. Any future build must have an isolated `TERLUNA_CM1_HOME` and unique build/run destination, or an explicitly preserved executable. Do not use setup/run commands against these archive paths merely to test a merge.

Literature and smoke histories are under `/media/projectspace/terluna-research/atmospheric-electricity/`, including aerosol notes, charging, ionization, humidity, leader, nitrogen and benchmark notes; `smoke/race_bisect.py`, restart and determinism checks. Source manifests contain 57 electricity and 26 aerosol entries. The electricity manifest's top-level “none imported/installed” and its old S2 overview entry are inherited planning text alongside later implementation entries; reconcile their time scope. No literature was re-read externally for this review, and none of these internal source records admits evidence to manuscripts automatically.

## Checks: historical record versus this review

Historical final log actually inspected: `/media/projectspace/terluna-research/atmospheric-electricity/make_check_2026-10-08_status_date.log`, explicitly rooted in this worktree and ending `exit 0`. It records:

- layer check: 314 files, zero violations;
- Python: 486 passed, 29 skipped, 26.59 seconds;
- JavaScript: 128 immersion tests plus 48 column/light-transport/viewer tests, all passed;
- provenance: `PASS_WITH_BLOCKED_INPUTS`; protection verifier not executed because `TSIS1_HSRS_stride100.csv`, `TiO2_Siefke.yml`, `o.nff`, `si.nff`, `ti.nff`, `SiO2_Franta.yml` were absent or changed;
- ensemble validator: PASS, 14 snapshot files and five paper workspaces, integrity only, no scientific/manuscript clearance.

Cap and storm-lives logs were also inspected and show the same pass/skip counts. These are historical logs, not a new `make check` result and not full empirical validation. Historical documented solver checks include numerical Poisson comparisons, synthetic positive/negative ground flashes, periodic boundary cases, screening conservation, Earth benchmark comparisons and single-thread/restart reproducibility. Multi-thread bit reproducibility specifically fails in some cases; Earth charging remains faster than published scheme runs.

Actual checks performed during this review: Git history/diff/file-overlap reads; JSON structure inspection; the listed producer/source-input hash comparisons; exact three external-input integrity checks; archived case/progress/executable-path presence; final log inspection. No pytest, model run, output regeneration, downloaded literature, manuscript review or PDF review was performed.

## Concrete verification and integration sequence to propose

1. Pin source heads and retain cap commits; preserve imported byte-pinned science and the five-paper UUID/planning files. The electricity code can merge largely intact; resolve shared docs by subject, making the corrected run current and marking stage 3 pending. Incorporate other branch findings as cited questions/inputs rather than adopting their recommendations by implication.
2. Record exact external input and archive availability in the unified workspace; keep current archive paths read-only for analysis. Check `gfortran`, NumPy/SciPy/xarray as required, verified CM1/WRF source paths and Takahashi/CRII tables before counting skipped tests as coverage.
3. Run the focused tests once, serialized with other sessions and with BLAS constrained. Candidate command from the unified checkout (choose a verified interpreter):

   ```sh
   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -ra \
     atmosphere/tests/test_electricity_charging.py \
     atmosphere/tests/test_electricity_conductivity.py \
     atmosphere/tests/test_electricity_laws.py \
     atmosphere/tests/test_electricity_thunder.py \
     climate/tests/test_cm1_elec.py \
     climate/tests/test_cm1_elec_field.py \
     climate/tests/test_cm1_elec_branched.py \
     climate/tests/test_cm1_elec_screen.py \
     climate/tests/test_cm1_run.py \
     climate/tests/test_elec_analysis.py \
     climate/tests/test_mixed_phase_analysis.py \
     climate/tests/test_storm_lives.py \
     research/studies/atmospheric_electricity \
     research/studies/open_moon_aerosol
   ```

   Fortran tests compile small harnesses in temporary directories; some skip without gfortran/CM1/WRF source. Conductivity/table and aerosol integration tests also skip if ignored inputs are missing. Report these separately. Do not start a full CM1 storm, MCEq cascade or global climate run to satisfy this check.
4. Reconcile constants and product headers/dependency pinning as a separate reviewable change. If a diagnostic must be rerun, redirect its output into an isolated comparison area; several current entry points otherwise overwrite tracked results. Compare numerical values before deciding to replace historical products. Fast recomputation can cover aerosol/nitrogen/stage 1 once inputs are verified; thunder takes about ten minutes, MCEq 4/16 minutes, and full storms are a later campaign, not routine merge validation.
5. Run `make check` on the unified tree before the eventual integration commit, with actual skip/block states reported. Run `python research/check.py --require-inputs` only when the exact protection inputs have been restored; a permissive provenance PASS_WITH_BLOCKED_INPUTS is not a successful protection-model reproduction. Preserve historical logs separately from unified checks.
6. At the end of integration, require no unintended source-worktree changes, no live archive executable overwritten, no stale “running/start held/full analysis pending” current statuses, no stage-3/global-circuit claim, and no unqualified global lightning/nitrogen/safety claim. The next research program can then be chosen from stage-3 current diagnostics, conductor/craft exposure, finer/wider lunar storms, aerosol transport and plasma closure, without hiding the remaining climate uncertainty.
