# Biosphere

The [accepted companion seed](../ensemble/papers/74a946e2-8f98/seed-02.md) owns biological viability, human biological requirements, ecology and biogeochemical persistence across the whole Moon.

## Executable requirements models

[long_night.py](long_night.py) now provides two bounded calculations:

- Periodic carbon reserves: cycle-integrated surplus and the maximum cumulative deficit determine the minimum reservoir capacity for a prescribed production/demand cycle. The necessary-and-sufficient condition is proved under its model assumptions in [findings](../research/findings.md), with an O(N) implementation tested against exhaustive searches.
- Mixed aquatic oxygen: production, respiration, exchange and a finite concentration ceiling are balanced exactly for each interval. Outgassing and unmet aerobic demand are reported explicitly. A mathematical periodic state can still fail the selected concentration or demand criteria.

Trait rates, suppression, Q10, exchange coefficients, dissolved-oxygen ceiling and thresholds are hypothetical inputs. Units and every parameter are recorded. The models identify requirements; their outputs are not survival, health, reproduction, nutrient-cycle or ecosystem evidence. Excess carbon is unallocated surplus, not proven growth. The temperature/irradiance interface demonstrates how a climate trace alters these requirements; tissue injury and photosynthetic saturation remain unmodelled.

[Results](../research/results/environment_screens/) include 180 carbon cases, 20 oxygen cases and 18 climate-to-carbon demonstrations. Species-level calibration and long-night experiments are the next evidence step.

## Remaining biological work

The wider portfolio retains soils, aquatic communities, detrital/subsurface habitats, varied plant architectures, aerial exchange and human developmental requirements. Nutrient compartments, ecological interactions, plant hydraulics, structural support and complete life cycles still need separate models and empirical tests. Megaforests remain one candidate within this scope.

Run `python -m unittest discover -s tests -v` and `python research/run_environment_screens.py` from the root. Provisional research sources and their actual access status are in [environment_sources.json](../research/environment_sources.json); older sources and archives remain linked through [archive_status.md](../research/archive_status.md).
