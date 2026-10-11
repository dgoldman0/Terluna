# Sky ecology status

**The guards and symbioses analysis extends the study from the checkpoint at commit `980d203`.** It answers the author's
direction of 10 October 2026 for a deeper analysis of the sky reefs' guards and wider symbioses, with the guards free to
be several species ([decisions](../../decisions.md#sky-reefs-partners-and-control-10-october-2026)). The
[analysis](symbioses.md) and its separate [requirements product](results/symbioses.json) add live-mass scenarios,
individual sizes and active fractions, local response geometry, shared carbon allocations, paired dietary carbon and
phosphorus, and nutrient imports. Seventeen selected primary studies support the mechanisms and set their limits;
[symbioses_sources.json](symbioses_sources.json) records how each was accessed. The organisms and populations are
unvalidated proposals.

**The saved sky ecology product is regenerated with the current source binding, and its accounting is repaired.** Its
phosphorus-limited consumer ledger reports the production the phosphorus limit allows and an explicit
unallocated-assimilate pool, which keeps the carbon accounted for while its fate stays open. Static trim estimates are
labelled as mass equivalents, and the tenants note gives their timing. Default study runs write to ignored
research/runs/sky_ecology/, and the research and biosphere indexes link the study.

**The independent review's second and third passes found problems in the committed code, and both sets are fixed.** Each
edge of the phosphorus loop takes the harvest share that bounds the need, which puts the sky's need at 1% cover at
62–424 kt a year mature and 71–528 kt growing with nothing to catch the litter, and at 3.5–61 kt and 13–165 kt with
forest-floor pools. Replacing the dead costs every body the hazard times its construction carbon: at the high hazard
500-m and 1-km round giants and 150-m reef modules on fermented gas cannot pay it from their surplus, the mature figures
leave them out, and the hazard is set against the growth note's ages and the storm note's lightning. The sudden-load
tables carry their timing: a 25-minute buoyancy period, a step that overshoots its 190-m sink by 42–97%, venting that
takes 71% of the trim equivalent at once, and dropping water that takes all of it in seconds to minutes.
[checks.json](checks.json) lists every finding and fix.

**Defensive efficacy, coexistence, partner inheritance, net host benefit, complete diets, gas-control biology, wet
habitat persistence and added-load damage tolerance are open.** The response and food scenarios are conditional
requirements. [symbioses_checks.json](symbioses_checks.json) records the current validation, and
[checkpoint_checks.json](checkpoint_checks.json) preserves the earlier failing checkpoint unchanged.

## Partial checkpoint on 10 October 2026

The author requested a commit of the interrupted sky ecology work before the deeper analysis of guards and general sky
reef symbioses. The checkpoint preserved the model code, source register and saved product as they stood, with a status
notice added to the README, and recorded the study as partial. It held coupled accounting for host production, trophic
transfers, tenant habitat and loading, phosphorus recycling, and transport along prescribed wind routes; its species
roles, tissue traits, retention efficiencies and population equivalents were hypotheses or conditional requirements.

[checkpoint_checks.json](checkpoint_checks.json) records the required repository checks: 1,645 Python tests passed, two
failed on the known stale source binding, and 14 were skipped. All 191 JavaScript tests, the layer check, historical
provenance and ensemble workspace integrity passed. The checkpoint kept those two failures deliberately, as a partial
record of work in progress.

The unfinished work the checkpoint listed now stands as follows:

- **Source binding.** The saved product recorded sources.json SHA-256
  `5e2e2423359f511014137d3a287f96fb99636869b7f789cef82cdcd9bab313c0` while the file had become
  `39f59e82b8c88def4842215832a19c32ea756d6c736c13c15cec53a86489e868`, because the source register changed after
  regeneration; the exact-regeneration check could pass only after the product was regenerated and reviewed. The product
  is now regenerated with the current binding.
- **Consumer carbon.** food_web.trophic_chain limited production by phosphorus after calculating the carbon-only
  consumer ledger, and the carbon excluded from production had no assigned fate in the returned ledger, so carbon closed
  only in the carbon-only ledger. The unallocated-assimilate pool now carries that carbon.
- **Trim without a timescale.** Tables labelled sudden or at once converted the pressure allowance into equivalent mass
  changes, with no arrival duration, airflow or valve rate, landing impulse, vertical response or recovery time; the
  0.63 and 2.1 tonne allowances of 100 and 150 metre reef modules were static, untested as instantaneous landing limits.
  The tenants note now gives the timing.
- **Added resident loads.** The spare-lift and gas-carbon arithmetic leaves gas-cell damage tolerance, attachment
  mechanics, water use and resident respiration to be rerun at each proposed extra load, and the masses are lift
  allowances with no food supply behind them. This remains open.
- **Guard assumptions.** The 1 to 10 grams per square metre entry was a guessed collective guard mass per projected reef
  area, with individual size, wet or dry mass and the threats it answers unspecified, no guard organism selected and no
  single species required to do every job. The symbioses analysis now works in live mass, individual sizes and several
  guilds.
- **Community persistence.** A self-sustaining community needs a population model, a complete coupled nutrient cycle,
  night survival, pool oxygen histories and weather trajectories, and the strong litter retention is assumed, awaiting
  derivation from catchment geometry and flow. This remains open.
- **Integration.** The research and biosphere indexes now link the study. The original checks.json kept the earlier
  focused tests and their limits, and the checkpoint's own checks were recorded separately in checkpoint_checks.json.

The checkpoint also planned the multi-species analysis of guards and broader symbioses that [symbioses.md](symbioses.md)
now holds: complementary threats and services, host and partner costs, residence against reciprocal benefit, body size
against total biomass, local response times, shared carbon and nutrient budgets, nutrient import and recycling,
recruitment, night survival and conflicts between partners, with Earth precedents kept separate from the proposed sky
reef organisms.
