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
| [provenance.json](provenance.json), [archive_status.md](archive_status.md) | Every original archive member and its disposition; unresolved recovery |
| [check.py](check.py), [checks.json](checks.json) | Byte integrity of the pinned imports and reproduction of the historical baseline |
| [baselines/feasibility](baselines/feasibility/) | The September multi-domain feasibility script, kept intact |
| [studies/](studies/) | Cross-domain studies, each with its runner, results and write-up |

## Domains

| Topic | Current calculation | Principal open condition |
|---|---|---|
| [Atmosphere](../atmosphere/) | Solved conduction/advection/Jeans column, band ledger, lower air, column soundings, line-by-line radiative–convective column with shield-filtered sunlight, correlated-k thermal scheme, global-mean radiative–photochemical equilibrium (ozone, stratosphere, surface UV) per shield with solved balance temperatures, non-LTE CO2 cooling bounds, exobase with idealised-filter leakage heating and atomic-oxygen escape, line-by-line radiation benchmarks for a GCM | 3-D middle-atmosphere transport and its effect on the exobase; the titania stack's EUV response; measured fluorinated-gas cross-sections |
| [Climate](../climate/) | Conservative periodic latitude–longitude thermal screen; canopy and forest-patch flow; clear-sky line-by-line budget from the atmosphere domain; single-column month-long day and night over land and sea; GCM boundary files and experiment plan (no GCM run) | Clouds, moisture/ice, dynamics and terrain |
| [Biosphere](../biosphere/) | Carbon reserve theorem/model, bounded oxygen box, tree and forest-patch mechanics | Measured traits, complete life cycles and ecological interactions |
| [Illumination](../illumination/) | Spectral clear-sky solver, A1–A3 light-transport references, site sky geometry and earthlight, bright stars | Date-accurate ephemeris, earthlight spectrum, profiles tied to the solved column |
| [Protection](../protection/) | Original component model plus spectral coverage audit | EUV response, clean operation, particle transport and lifetime resources |
| [Engineering](../engineering/) | Original transport/renewal accounts plus plume-heat sensitivity | Complete industrial network and safe source-to-use routes |
| [Geography](../geography/) | LOLA topography above the GRAIL geoid; hydrostatic water storage (level curves, basin joins) | Where water stands once rain, runoff and evaporation act |
| [Habitation](../habitation/) | Design concepts | Practical inhabited capacity |

Rendering that displays these results is in [visualization](../visualization/);
the explorable experience, which reads them as baked products, is in
[immersion](../immersion/). Constants and scenarios shared by every lane are in
[shared](../shared/).

## Studies

| Study | What it couples | Boundary |
|---|---|---|
| [environment_screens](studies/environment_screens/) | Thermal column, band ledger, climate screen, carbon and oxygen requirements, protection plume heat | Conditional consequences of selected equations; the spectrum/material/chemistry interface is unfilled |
| [megaforest_wind](studies/megaforest_wind/README.md) | Lower air (atmosphere), canopy flow (climate), tree mechanics (biosphere): 480 static load envelopes | Material and soil traits hypothetical; wind climate, gusts and evolution unmodelled |
| [forest_patch](studies/forest_patch/README.md) | Three-dimensional patch airflow, tree load sharing, compliant foundations, shared soil and damage | Wind climatology, gusts, nonlinear failure and evolved morphology open |

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
