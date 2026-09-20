# Engineering and industry

**Current code:** the original cross-domain [feasibility model](../research/baselines/feasibility/model.py). **Reference data:** [logistics](reference/logistics.csv), [routes](reference/transfer_routes.csv), [scalar growth](reference/growth.csv) and [maintenance](reference/maintenance.csv). These are verbatim September outputs, not new calculations or demonstrated hardware.

The route model is Sun-only and circular/coplanar. Electric propulsion includes exhaust and propellant but omits dry mass, return fleets and complete extraction/capture systems. Growth is a required capacity curve with an imposed cap. Maintenance tables impose depletion timescales rather than predict escape. Oxygen process yields and energy benchmarks retain their separate process assumptions.

## Next work

- Add water, nutrients, surface sinks and infrastructure to a shared material ledger.
- Compare a few complete source-to-use pathways, including containment, braking, reuse and failure geometry.
- Develop local power/fuel/heat accounts. Fusion is a conditional technology scenario; no reactor has been designed or validated here.
- Replace scalar growth with a small network of manufacturing bottlenecks, replacements and critical imports.
- Coordinate atmosphere buildup, major heat-intensive operations, surface conditioning, biological establishment and service renewal.

The [protection](../protection/) calculations supply component hardware and operating budgets, not a complete atmospheric-retention solution. Sagan K≈1.17 remains an illustrative civilization-allocation scenario pending consistent whole-system energy accounting. Planetary traffic safety is a design requirement; no certified trajectory network has been established.

Run the shared baseline with `python research/baselines/feasibility/model.py --out research/runs/feasibility`. The [core-first plan](../research/plan.md) keeps industry one part of the whole-world argument.
