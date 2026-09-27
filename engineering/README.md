# Engineering and industry

**Current code:** the original cross-domain [feasibility model](../research/baselines/feasibility/model.py), with its full [report](../research/baselines/feasibility/report.md), the start of the [industrial network model](network/README.md): its vocabulary and a supply ledger for the protection hardware, [tower models](towers/README.md) that size very tall lattice frames for gravity and wind and screen wind devices built into them, used by the [summit tower study](../research/studies/summit_tower/README.md), and [fire models](fire/README.md) that carry Earth's fire-engineering correlations to lunar gravity and the Open Moon's air, used by the [port fire study](../research/studies/port_fire/README.md). **Reference data:** [logistics](reference/logistics.csv), [routes](reference/transfer_routes.csv), [scalar growth](reference/growth.csv) and [maintenance](reference/maintenance.csv). These are verbatim September outputs, not new calculations or demonstrated hardware.

The route model is Sun-only and circular/coplanar. Electric propulsion includes exhaust and propellant but omits dry mass, return fleets and complete extraction/capture systems. Growth is a required capacity curve with an imposed cap. Maintenance tables impose depletion timescales rather than predict escape. Oxygen process yields and energy benchmarks retain their separate process assumptions.

## Sourcing and conservation

- **Where the air and water come from.** Water, nitrogen and oxygen all come from the Solar-System-wide resource operation ([conservation study](../research/studies/conservation/README.md), decision D4). The baseline's local oxygen production at 24.3 kWh/kg is superseded.
- **Why oxygen is shipped.** Making the design's oxygen from lunar rock would mean processing 0.5–1 million km³, two to three times the Moon's whole regolith, and would leave about 1–2.5×10¹⁸ kg of reduced residue.
- **Shipped as water.** The oxygen for the 17.5% design comes to 6.2×10¹⁷ kg of water, about 6% of the 28% seas. The 7×10¹⁶ kg of hydrogen left over serves as propellant or as feedstock for ammonia and biomass.
- **Source bodies.** Sources follow the study's principle 5: common icy bodies, with Saturn's rings, Titan's atmosphere and the ocean worlds left untouched.

## Industrial architecture from 19 September

A reconstruction of the industrial design worked out with the author on
19 September 2026 is in
[reference/industrial_architecture/](reference/industrial_architecture/README.md),
written on 2026-09-26 from the conversation record; no file from 19 September
survives. It keeps three power ledgers apart
(civilization-wide power, the Open Moon's allocation, and where heat is
released), sequences power from fission bootstrap through fusion to near-Sun
collection with power beaming, and builds industry at the source bodies. Bulk
material moves as a stream of standardized packets, about 2.7×10⁸ kg/s with
some 5×10¹⁶ kg in transit, and a cislunar terminal layer brakes, captures,
inspects, repacks, buffers and meters it; the Moon itself is never the brake.
The author confirmed its six traffic-safety rules as decisions on 2026-09-26
([research/decisions.md](../research/decisions.md)). The document ends with a
specification for the network model listed below.

## Next work

- Add water, nutrients, surface sinks and infrastructure to a shared material ledger.
- Compare a few complete source-to-use pathways, including containment, braking, reuse and failure geometry.
- Develop local power/fuel/heat accounts. Fusion is a conditional technology scenario; no reactor has been designed or validated here.
- Replace scalar growth with a small network of manufacturing bottlenecks, replacements and critical imports.
- Coordinate atmosphere buildup, major heat-intensive operations, surface conditioning, biological establishment and service renewal.
- Heritage works from the conservation study: the undersea Tranquility community (pressure hull, anchoring against ground uplift, access tower), enclosed land sites and the sealed reference area.
- Replacements for the Moon's scientific functions: traceable calibration satellites, a shielded deep-space radio array and ranging stations on high ground.
- Return of eroded sediment and nutrients uphill, 0.7–7 Gt a year at 0.1–1.1 GW ideal power.

The [protection](../protection/) calculations supply component hardware and operating budgets, not a complete atmospheric-retention solution. Sagan K≈1.17 remains an illustrative civilization-allocation scenario pending consistent whole-system energy accounting. Planetary traffic safety is a design requirement; no certified trajectory network has been established.

Run the shared baseline with `python research/baselines/feasibility/model.py --out research/runs/feasibility`. The [core-first plan](../research/plan.md) keeps industry one part of the whole-world argument.
