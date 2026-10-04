# Terluna research

The hub for the project's shared research: what exists in each domain, the
cross-domain studies, the findings, the plan, the condition register and the
provenance of imported material. It supports the five-paper
[Open Moon ensemble](../ensemble/). **Constructing and Sustaining an Open Moon is
the next full manuscript**; the core's opening landscape does not define the
limits of the research.

Models live in their domain folders. This folder holds what spans domains:

| File | Holds |
|---|---|
| [findings.md](findings.md) | Equations, bounds, source limits and numerical examples from the environment screens |
| [plan.md](plan.md) | The core-first work order |
| [status.json](status.json) | Condition, open questions and next task for each topic |
| [decisions.md](decisions.md) | The author's decisions, where each is recorded, and whether it stands |
| [provenance.json](provenance.json), [archive_status.md](archive_status.md) | Every original archive member and its disposition; unresolved recovery |
| [check.py](check.py), [checks.json](checks.json) | Byte integrity of the pinned imports and reproduction of the historical baseline |
| [baselines/feasibility](baselines/feasibility/) | The September multi-domain feasibility script and its full report, kept intact |
| [studies/](studies/) | Cross-domain studies, each with its runner, results and write-up |

## Domains

| Topic | Current calculation | Principal open condition |
|---|---|---|
| [Atmosphere](../atmosphere/) | Solved conduction/advection/Jeans column, band ledger, lower air, column soundings, line-by-line radiative–convective column with shield-filtered sunlight, correlated-k thermal scheme, global-mean radiative–photochemical equilibrium (ozone, stratosphere, surface UV) per shield with solved balance temperatures, non-LTE CO2 cooling bounds, exobase with idealised-filter and titania-film heating and atomic-oxygen escape, line-by-line radiation benchmarks for a GCM, fluorinated trace-gas forcing from measured bands, the exobase with less CO2, the loss response to protection (screening) | 3-D middle-atmosphere transport and its effect on the exobase; light bypassing the shield film through aperture gaps |
| [Climate](../climate/) | Conservative periodic latitude–longitude thermal screen; canopy and forest-patch flow; clear-sky line-by-line budget from the atmosphere domain; single-column month-long day and night over land and sea; GCM boundary files, experiment plan and a restartable ExoPlaSim runner (design case near 299.6 K with PlaSim's cloud scheme corrected for lunar gravity; water and cloud sensitivities); [the seas' waves](../climate/waves/): two lunar cycles of SWAN waves over the nearside sea and Smythii–Marginis from the GCM's surface stress (mean Hs 1.39 and 0.88 m), coastal nests to 118 m terrain, and SWASH breaking and run-up | Climate programme paused on 2026-09-30 (decisions register): CM1 and the GCM disagree on warmth and humidity over land, so comfort is carried as a range; cloud amount, dynamics and terrain. Waves: further lunar cycles and the South Pole–Aitken sea, the monthly tide's water level at the coasts, a finer basin grid, sediment on regolith shores |
| [Biosphere](../biosphere/) | Carbon reserve theorem/model, bounded oxygen box, tree and forest-patch mechanics, clear-sky canopy photosynthesis by wavelength under the Moon's and Earth's skies, the whole plant's carbon through the lunar night with twilight, fruit designed for the lunar day on that plant (program, set time, size and share) | Measured traits, complete life cycles and ecological interactions; whether leaves survive weeks of darkness and use two weeks of continuous light; whether a fruit's program can be compressed to one lunar day |
| [Illumination](../illumination/) | Spectral clear-sky solver, A1–A3 light-transport references, site geometry and earthlight, bright stars, surface spectra, a shielded spherical sky on the solved columns and regional evening cloud radiance with evolving views | Regional profiles, three-dimensional cloud structure, aerosols, refraction, polarization and the earthlight spectrum |
| [Protection](../protection/) | Original component model and its full September design report, plus spectral coverage audit; the stored film's transmission from hard X-rays to the far ultraviolet | Aperture leakage (gaps, pinholes, edges), clean operation, particle transport and lifetime resources |
| [Engineering](../engineering/) | Original transport/renewal accounts plus plume-heat sensitivity | Complete industrial network and safe source-to-use routes |
| [Geography](../geography/) | LOLA topography above the GRAIL geoid; hydrostatic water storage (level curves, basin joins); the atlas at the selected 28% water share, named from the IAU gazetteer; a first estimate of rivers and rain-fed lakes from the climate run's runoff; explicit cross-shore profiles for the wave studies; [every sea's monthly tide](../geography/README.md#the-monthly-tide) (3.7 m typical range over the nearside sea) with a dynamic check | How the water inventory divides between seas and rain-fed lakes (lakes hold about a third of the seas' volume in the first estimate); groundwater; Mare Fecunditatis's strait, which sets its fortnightly tides |
| [Habitation](../habitation/) | Design concepts; a first optical-comfort screen with surface luminance, eye illumination, shadow contrast and Earth controls | Practical inhabited capacity; visual comfort in resolved scenes and human validation |

Rendering that displays these results is in [visualization](../visualization/);
the explorable experience, which reads them as baked products, is in
[immersion](../immersion/). Constants and scenarios shared by every lane are in
[shared](../shared/).

## Studies

| Study | What it couples | Boundary |
|---|---|---|
| [environment_screens](studies/environment_screens/) | Thermal column, band ledger, climate screen, carbon and oxygen requirements, protection plume heat | Conditional consequences of selected equations; the spectrum/material/chemistry interface is unfilled |
| [megaforest_wind](studies/megaforest_wind/README.md) | Lower air (atmosphere), canopy flow (climate), tree mechanics (biosphere): 480 static load envelopes | Material and soil traits hypothetical; wind climate, gusts and evolution unmodelled |
| [conservation](studies/conservation/README.md) | Geography (the atlas), atmosphere (design pressure), climate (albedo, polar temperature), engineering (enclosures, sourcing) and a written heritage evaluation: the two preservation programmes, gates, register, scientific zones and community settings | Heritage evaluation is a first pass awaiting the author's review; numbers are hydrostatic geometry and screening estimates |
| [forest_patch](studies/forest_patch/README.md) | Three-dimensional patch airflow, tree load sharing, compliant foundations, shared soil and damage | Wind climatology, gusts, nonlinear failure and evolved morphology open |
| [atmospheric_co2](studies/atmospheric_co2/README.md) | Biology (plant CO2 needs, a C3 leaf model), weathering (a basalt law) and atmosphere (CO2 forcing, exobase): the CO2 level plants need, drawdown by biosphere build-up and weathering, and the return a settled Moon needs | Literature synthesis with Earth rates; no box model, ocean uptake or lunar soils |
| [protection_architecture](studies/protection_architecture/README.md) | Atmosphere (the loss response), climate (dimming), protection (hardware), engineering (supply) and the decisions register: the requirements for the optical shield and charged-particle protection, and what each loss budget asks of them | Screening loss model with the exosphere's losses; the loss budget, the ultraviolet cut-off and the protected radius are the author's to set; hardware, holding and supply are still to design |
| [lunar_cycle_ecology](studies/lunar_cycle_ecology/README.md) | Biosphere (plant and fruit carbon), geography (seas, lakes, runoff), atmosphere (air column, flight) and illumination (darkness): a register of design ideas for plants, animals, the night food web and the nutrient and mineral cycles, with light calculations and a recommended selection | Ideas and screening arithmetic; no ecosystem, population or nutrient-cycle model; the selection is the author's to make |
| [optical_comfort](studies/optical_comfort/README.md) | Illumination, surfaces and human vision: baseline screens plus 630 gaze cases and 24 finite matte scenes using a consistent shielded spherical sky, with an Earth control and local cloud and water-slope inputs | Conditional clear molecular atmosphere and imposed scene geometry; cloud radiance, measured surface reflectance, rough-water reflection and human-factors evidence remain open |
| [cloud_twilight](studies/cloud_twilight/README.md) | Second-cycle regional cloud survey, spectral cloud radiance and colour, an 8–80 h fixed-observer history and sheltered viewing deep into evening | One sampled cycle, explicit cloud extrusion and particle-optics scenarios; three-dimensional regional structure, terrain and cloud Earthlight remain further inputs |
| [sea_appearance](studies/sea_appearance/README.md) | Climate (the waves), geography (the tide and the seas), illumination (sky, clouds, earthlight) and biosphere (living water): how the seas look by day, through the long evening, under earthlight and in the far side's twilit night | Plan adopted 2026-10-04; lighting calendar of six coasts and the Earth over the seas done, giving dawn and dusk their photometric definition; the sea's slopes through the month, the reflection model and the waters' colours done; clear sky, uniform column, earthlight with sunlight's colour |

The forest-patch folder also holds a [reviewed checkpoint](studies/forest_patch/reviewed/README.md)
from a separate lineage (57 tests, 90-m spacing, 20-m/s reservoir flow), whose
source is kept at commit `363b161` and on `checkpoint/forest-patch-review-57-tests`.
The two lineages differ in their assumptions. In both, the dome was an imposed
candidate: neither simulates an evolved canopy shape or establishes an optimum.

Run everything from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest
OPENBLAS_NUM_THREADS=1 python -m research.studies.environment_screens.run   # results in studies/environment_screens/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.megaforest_wind.run       # results in research/runs/megaforest_wind
OPENBLAS_NUM_THREADS=1 python -m research.studies.forest_patch.run --quick  # results in research/runs/forest_patch
python -m research.studies.atmospheric_co2.run                             # results in studies/atmospheric_co2/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.protection_architecture.run # results in studies/protection_architecture/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.lunar_cycle_ecology.run     # results in studies/lunar_cycle_ecology/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.run
OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.directional  # results in studies/optical_comfort/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.cloud_twilight.run      # results in studies/cloud_twilight/results
```

For the environment screens, supply `--protection-archive /path/to/Lunar_Protection_Model.zip`
to inspect the original optical inputs directly. Source-access and validation
status are in [sources.json](studies/environment_screens/sources.json) and
[checks.json](studies/environment_screens/checks.json).

## Historical baseline reproduction

```sh
python -m pip install -r research/requirements.txt
python research/check.py
# Restore the exact five protection inputs from the original project ZIP:
python protection/fetch_inputs.py --archive /path/to/Lunar_Protection_Model.zip
python research/check.py --require-inputs
```

Upstream retrieval is also available through `python protection/fetch_inputs.py --download`;
the exact expected hash and size must match. With missing protection data the
checker reports `BLOCKED` explicitly; it never substitutes synthetic inputs.

The [checked snapshot](checks.json) records an actual run using the original
archived inputs. Future checks write to ignored `research/runs/`. Numerical
reproduction and file integrity do not grant physical validation, full-source
admission, or manuscript clearance.

The September feasibility implementation is a historical multi-domain script;
its single intact copy lives in [baselines/feasibility](baselines/feasibility/),
and topic folders hold its original compact reference tables. Protection has its
own intact implementation. Model code and CSV tables are verbatim; the two JSON
reference files select original fields with unchanged values. New calculations
should state assumptions, track conservation and residuals, test convergence and
say which conclusion their outputs can change.
