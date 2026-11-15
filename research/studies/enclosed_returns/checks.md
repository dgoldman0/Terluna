# Verification record

This conditional research appraisal is written for the project author and researchers reviewing the industrial and economic assumptions. It distinguishes the prescribed architectural target, numerical consequences of stated scenarios, and engineering or institutional requirements that remain unqualified. The study and its current synthesis were read as complete prose, including the interpretation of tables and the connection between financial and source-supply conditions.

## Admitted inputs

The provisioning source is commit `3e04591148223b079ecd255416555b8d4a4d4f7e` on `domain/provisioning`. The resources source is commit `38cd446f361f19ff732e8f511ee87104dd51a083` on `domain/resources`. Exact compressed product hashes, schemas and source paths are recorded in `finance_inputs.json`; the resources producer additionally verifies decompressed identities. The consumer verifies those products and the shared customer orderbook before calculating returns.

The main physical product contains 81 cases. Three additional products supply two cases each for poorer lunar grade, poorer source grade and deterministic component lives. All 87 cases pass their annual rated-equipment and workforce screens. Those numerical results do not qualify mines, transport trajectories, installation or complete deployment.

## Numerical and repository checks

The final financial product was generated with:

```sh
PYTHONDONTWRITEBYTECODE=1 python -m research.studies.enclosed_returns.finance --output research/studies/enclosed_returns/results/returns.json
```

It contains 783 checkpoint accounts: 87 physical cases, three financing modes and years 100, 200 and 220. The saved-product test reproduces the complete result exactly from the admitted inputs. Twenty-eight study tests cover customer-budget competition, the residential schedule, independent commercial demand, capital and renewal timing, seed payments, service-quote overlap, cash/debt identities, fixed debt ceilings, refusal of loss-funding credit, resource conditions, financing-neutral social value, forecast retention and an unanticipated receipt shock.

The main checkout's affected-check selector was directed at the humanities worktree, using the selector from `/home/kir/Documents/Projects/terluna/shared/check_changed.py` with its `ROOT` set to this worktree. Its selected commands were:

```sh
make check-layers PYTHON=/home/kir/miniconda3/bin/python
python -m pytest atmosphere/tests/test_snow_crystal_transport_replay.py humanities/tests/test_cultural_buildout.py::test_buildout_regenerates_and_preserves_bound_annual_products research/studies/enclosed_returns/tests/test_finance.py research/studies/enclosed_returns/tests/test_markets.py research/studies/enclosed_returns/tests/test_saved.py
```

**Result: 74 tests passed; the layer check inspected 1,672 files and found zero violations.** The initial invocation stopped because this branch's older layer checker tried to read deleted files still present in the unstaged Git index. Staging the intended deletions resolved that file-inventory issue; the selected checks then passed. The selector's uncovered paths are the withdrawn implementations and the figure renderer; withdrawal was inspected, and the renderer was checked separately below.

The owning provisioning work reported 92 passing repository-selected tests, one optional large GCM-extraction skip, and passing layer/provenance checks. The owning resources work reported 64 passing tests, one unrelated optional phosphorus-grid skip, exact reproduction of its 87-case assessment and zero layer violations. These are separate check scopes, not a count of unique tests across the repository. Full repository verification was not requested or run for this bounded correction.

An independent accounting review compared the financial consumer with every annual row of all 81 base physical cases. Capital plus operating costs reconcile to physical costs, commercial operations, residual service allowances, hotel payroll and the explicit capital allowance, with a maximum floating-point residual below USD2041 0.01. That residual measures arithmetic consistency, not forecast accuracy. Identified resource provisions and resident payroll are counted once financially and removed once where required by the socioeconomic account.

## Figures and editorial review

The final renderer ran successfully:

```sh
python visualization/enclosed-returns/render.py
```

It generated PNG, SVG and a plotting manifest in ignored `research/runs/enclosed_returns/figures/`. Every plotted series was compared exactly with the final financial product, and the manifest's input hash matched that product. The PNG was inspected at its delivered size for labels, units, legends, year markers and the source/deployment qualifications. Titles were shortened and the 35% capital allowance made explicit during that visual review. Generated images remain outside Git.

The results and methods were read together. That review clarified the potential audience ceilings, the difference between annual operating cash and construction funding, the narrow mature margin of the diversified case, and the dependence of the earlier manufacturing mix on its higher product price. It also brought the source-water/fleet cost chain and cumulative feed requirements into the integrated study. The current humanities, research and visualization introductions now point to the replacement account. Obsolete generated figures were removed; the withdrawn code and figures remain reproducible from their historical commits.

## Evidence limits

`finance_sources.json` records the selected primary-source appraisal sections consulted. This was not a complete scholarly review of either publication or of the broader economic literature. Customer spending, component lives, U settlement and credit rules remain explicit scenario assumptions. The capital allowance is treated as spent expenditure and adds no productive capacity. It does not establish a complete installation quote.

Cash, reserves and custody earn zero real return. The payout rule uses the stated future programme, while the separate receipt shock is first revealed in year 40. Post-gap completion accounts require the reported outside funding and do not estimate lender default recovery. Source access, grades, actual delivery routes, transit inventories, passenger transport and initial deployment still require qualification. Financial or socioeconomic net present value does not remove those requirements.
