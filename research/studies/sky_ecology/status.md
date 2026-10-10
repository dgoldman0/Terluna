# Sky ecology status

## Follow-up on guards and symbioses

The interrupted study is preserved in commit `980d203`. The current [analysis](symbioses.md) and separate
[requirements product](results/symbioses.json) add live-mass scenarios, individual sizes and active fractions,
local response geometry, shared carbon allocations, paired dietary carbon and phosphorus, and nutrient imports.
Seventeen selected primary studies support the mechanisms and their limits; access details are in
[symbioses_sources.json](symbioses_sources.json). No organism or population is validated.

The saved sky ecology product has been regenerated with the current source binding. Its phosphorus-limited
consumer ledger now reports actual limited production and an explicit unallocated-assimilate pool, preserving
carbon accounting without inventing the fate of that carbon. Static trim estimates are labelled as mass
equivalents with no landing or recovery time. Default study runs now write to ignored research/runs/sky_ecology/.
The research and biosphere indexes link the study.

The response and food scenarios are conditional requirements. Defensive efficacy, coexistence, partner
inheritance, net host benefit, complete diets, gas-control biology, wet habitat persistence and added-load damage
tolerance remain open. Current validation is recorded in [symbioses_checks.json](symbioses_checks.json);
[checkpoint_checks.json](checkpoint_checks.json) preserves the earlier failing checkpoint unchanged.

## Partial checkpoint on 10 October 2026

The author requested a commit of the interrupted sky ecology work before deeper analysis of guards and general sky reef symbioses. This checkpoint preserves the existing model code, source register and saved product, with a status notice added to the README. It does not declare the study complete.

The study contains coupled accounting for host production, trophic transfers, tenant habitat and loading, phosphorus recycling, and transport along prescribed wind routes. Its species roles, tissue traits, retention efficiencies and population equivalents remain hypotheses or conditional requirements.

## Known unfinished work at the checkpoint

- **Source binding.** The saved product records sources.json SHA-256 `5e2e2423359f511014137d3a287f96fb99636869b7f789cef82cdcd9bab313c0`, while the current file is `39f59e82b8c88def4842215832a19c32ea756d6c736c13c15cec53a86489e868`. The source register changed after regeneration. The old code and product are deliberately retained in this checkpoint; an exact-regeneration check cannot pass until the product is regenerated and reviewed.
- **Consumer carbon.** food_web.trophic_chain limits production by phosphorus after calculating the carbon-only consumer ledger. The carbon excluded from production has no assigned fate in the returned phosphorus-limited ledger. Carbon-only closure does not demonstrate closure of the phosphorus-limited food web.
- **Trim without a timescale.** Tables labelled sudden or at once convert pressure allowance into equivalent mass changes. They contain no arrival duration, airflow or valve rate, landing impulse, vertical response or recovery time. The 0.63 and 2.1 tonne allowances for 100 and 150 metre reef modules are not demonstrated instantaneous landing limits.
- **Added resident loads.** Spare-lift and gas-carbon arithmetic does not rerun gas-cell damage tolerance, attachment mechanics, water use or resident respiration at every proposed extra load. These masses are not food-supported populations.
- **Guard assumptions.** The 1 to 10 grams per square metre entry is a guessed collective guard mass per projected reef area. It does not specify individual size or wet versus dry mass, and does not establish defense against any threat. There is no selected guard organism or requirement that one species do every job.
- **Community persistence.** No population model, complete coupled nutrient cycle, night survival, pool oxygen history or weather trajectory establishes a self-sustaining community. Strong litter retention is assumed rather than derived from catchment geometry and flow.
- **Integration.** The research and biosphere indexes have not yet incorporated this study. The original checks.json records earlier focused tests and their limits; current checkpoint checks will be recorded separately there.

## Planned work at the checkpoint

The checkpoint called for a multi-species analysis of guards and the broader symbioses: complementary threats and services, host and partner costs, residence versus reciprocal benefit, body size versus total biomass, local response times, shared carbon and nutrient budgets, nutrient import and recycling, recruitment, night survival and conflicts between partners. Earth precedents must remain separate from proposed sky reef organisms.

## Checkpoint validation

[checkpoint_checks.json](checkpoint_checks.json) records the required repository checks: 1,645 Python tests passed, two failed on the known stale source binding, and 14 were skipped. All 191 JavaScript tests, the layer check, historical provenance and ensemble workspace integrity passed. The partial checkpoint intentionally retains those two failures; it is not a clean scientific release.
