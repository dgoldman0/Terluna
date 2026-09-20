# September feasibility baseline

`model.py` is a verbatim import from the 4 September 2026 project package. `reference.json` selects its original baseline, constants, inventory, tidal diagnostic and numerical checks; values are unchanged and duplicate table arrays are omitted. The original script combines atmospheric, thermal, transport, growth and maintenance calculations. One canonical copy is retained here; the topic folders contain its original CSV tables and links.

## Run

From the repository root, after installing `research/requirements.txt`:

```sh
python research/baselines/feasibility/model.py --out research/runs/feasibility
```

This regenerates the tables, aggregate results and four original plots outside the retained reference data. The original plotting code and its scientific assumptions are unchanged.

## Condition

- The temperature profile is prescribed; its baseline uses a 180 K middle atmosphere, a 250 K asymptotic upper target, 0.1 Pa heating boundary and 17.5% oxygen by number at 1.2 atm.
- The molecular profile holds composition fixed. The separate diffusive calculation imposes atomic fractions; neither predicts photochemistry or an energy-balanced outflow.
- The collision cross section is an illustrative neutral proxy. Very extended profiles require tidal and kinetic treatment beyond this calculation.
- The thermal response is a linear, spatially unresolved sensitivity model. Its heat capacity and outgoing-radiation slope are inputs.
- Hohmann routes are Sun-only circular/coplanar references. Propulsion omits dry mass, return fleets and complete source/destination operations.
- Growth solves the rate needed to deliver an inventory under imposed capacity assumptions. It does not forecast technological progress.
- Inventory/loss ratios are diagnostics. Energy-limited tables are separate ceilings and must not be added to Jeans escape.

The mathematical consistency tests cover grid convergence, an independent mass integral, an isothermal limit and zero-distance transfer. Their passage does not validate a habitable Moon. See [research status](../../status.json) and [current check record](../../checks.json).

The source archive's report, PDF builder and generated graphics are indexed in [archive manifests](../../provenance.json). They are outside this selected import; the code can regenerate its plots. External scholarly sources still require complete reading and claim-level admission.
