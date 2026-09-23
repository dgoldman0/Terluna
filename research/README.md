# Terluna research

Shared computational and theoretical support for the five-paper [Open Moon ensemble](../ensemble/). **Constructing and Sustaining an Open Moon is the next full manuscript.** Companion analyses develop alongside it; the core's opening landscape does not define the limits of the research.

Read [the new findings and equations](findings.md) and [numerical results](results/environment_screens/). The current computational pass adds a solved molecular thermal limit, a periodic spatial climate screen, and carbon/oxygen requirements models. It also identifies the unfilled spectrum/material/chemistry interface.

| Topic | New/current calculation | Principal open condition |
|---|---|---|
| [Atmosphere](../atmosphere/) | Solved conduction/advection/Jeans column and band ledger | Actual spectrum-to-chemistry-to-escape closure |
| [Climate](../climate/) | Conservative periodic latitude–longitude thermal screen | Calibrated radiation, moisture/ice, dynamics and terrain |
| [Biosphere](../biosphere/) | Carbon reserve theorem/model and bounded oxygen box | Measured traits, complete life cycles and ecological interactions |
| [Protection](../protection/) | Original component model plus spectral coverage audit | EUV response, clean operation, particle transport and lifetime resources |
| [Engineering](../engineering/) | Original transport/renewal accounts plus plume-heat sensitivity | Complete industrial network and safe source-to-use routes |
| [Geography](../geography/), [illumination](../illumination/), [habitation](../habitation/) | Data leads, angular diagnostic and design concepts | Actual terrain, visual appearance and practical inhabited capacity |


Read [status.json](status.json) for condition and next-task detail, [plan.md](plan.md) for the core-first work order, and [archive_status.md](archive_status.md) for unresolved recovery. [provenance.json](provenance.json) lists every original archive member and its disposition, including hashes for retained files and omitted material. Model code and CSV tables are verbatim; the two JSON reference files explicitly select original fields with unchanged values.

## New models

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 python research/run_environment_screens.py
```

Supply `--protection-archive /path/to/Lunar_Protection_Model.zip` for direct inspection of the original optical inputs, or restore `protection/sources/` using the existing helper. New runs use their own ordinary result paths. Source-access and scientific-validation status are recorded in [environment_sources.json](environment_sources.json) and [environment_checks.json](environment_checks.json).

## Historical baseline reproduction

```sh
python -m pip install -r research/requirements.txt
python research/check.py
# Restore the exact five protection inputs from the original project ZIP:
python protection/fetch_inputs.py --archive /path/to/Lunar_Protection_Model.zip
python research/check.py --require-inputs
```

Upstream retrieval is also available through `python protection/fetch_inputs.py --download`; the exact expected hash and size must match. Direct network retrieval was unavailable in the preparation container and was not newly verified. With missing protection data the checker reports `BLOCKED` explicitly; it never substitutes synthetic inputs.

The [checked snapshot](checks.json) records an actual run using the original archived inputs. Future checks write to ignored `research/runs/` by default. Numerical reproduction and file integrity do not grant physical validation, full-source admission, or manuscript clearance.

## Organization

The September feasibility implementation is a historical multi-domain script. Its single intact copy lives in [baselines/feasibility](baselines/feasibility/); topic folders hold its original compact reference tables. Protection has its own intact implementation. Large generated grids, rendered figures/PDFs, duplicate baseline copies and external spectral bytes are omitted from the curated commit; their source records and reproduction/restoration paths remain explicit.

Ordinary file names and Git history manage ongoing changes. The earlier dated planning snapshot and accepted seeds stay intact. New calculations should state assumptions, track conservation/residuals, test convergence and identify which conclusion their outputs can change.

## Megaforest wind checkpoint

The [wind methods and findings](megaforest_wind.md) add lower-air inputs, a
read-only A1 profile adapter, conditional canopy momentum flow, and static
plant/root/branch load envelopes. The [compact check record](results/megaforest_wind/checkpoint.json)
and [reference cases](results/megaforest_wind/reference_cases.csv) accompany
480 explicit scenarios. Material and soil traits remain hypothetical; wind
climate, gust dynamics and evolution are unmodelled. The methods identify
small-displacement limits and the next nonlinear-mechanics work.

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_megaforest_wind.py' -v
OPENBLAS_NUM_THREADS=1 python research/run_megaforest_wind.py
```

Full generated results go to ignored `research/runs/megaforest_wind/`. Source
ownership remains in `atmosphere/`, `climate/` and `biosphere/`; the runner stays
in `research/`. `immersion/` is unchanged and accepts no output from this runner.

## Spatial interacting-forest checkpoint

[Methods and findings](forest_patch.md) extend the earlier tree screen to explicitly
positioned forest patches, conservative spatial airflow/drag, root and crown
interactions, and prescribed damage sequences. The [runner](run_forest_patch.py)
compares edge shapes and matched aboveground-mass cases, with separate domain
limits. See [compact checkpoint](results/forest_patch/checkpoint.json).

This neutral steady-flow screen has no lunar wind climatology, resolved gusts or
evolved canopy prediction. Damage fronts occur in selected cases, while final
loss counts remain unresolved at the small-displacement boundary.
