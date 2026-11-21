# Resource requirements of bounded interventions

The [diagnostic adapter](bottlenecks.py) applies the existing finite-source and lunar-footprint calculation to the commercial-rollout and source-tug-life comparison. It reads a pinned physical product and preserves the existing catalogue, extraction fractions, assumed availability dates, lunar density and worked depth. It introduces no population-based feed multiplier, recovery credit, new mineral inventory or cost adjustment.

Each case retains its exact experiment name and physical feasibility status. An unconverged case without annual requirements receives no substitute feed trajectory. A converged case that fails workforce or equipment capacity still has a calculable requirement, which is labelled as the demand of that unresolved physical plan rather than actual production. Seed feed is included once in the eight average preproduction years, followed by the producer's dated annual feed.

The source-tug-life variant can reduce imported replacement mass and the associated chemical freight and source-feed requirements. Its feed changes come from the physical rerun, not a proportional reduction in this resource account. A catalogue mass allowance remains conditional on composition, extraction, permission and actual route availability, regardless of the variant's financial results.

From the existing resources checkout:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m resources.enclosed_mining.bottlenecks
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest resources/enclosed_mining/tests -q
```

The separate `bottlenecks_inputs.json` pins the physical product and manifest, the existing resource method and its reference inputs. The compressed diagnostic output is `results/bottlenecks.json.gz`; the accepted baseline assessment and its inputs remain unchanged.

## Computed comparison

The diagnostic output assesses all 28 dated requirement paths. Seven retain workforce failures from the physical producer, so their material totals describe unresolved plans. All four longer-tug-life cases pass the physical screen. Their cumulative source feed through year 200 falls by 1.46–1.47 percent relative to the corresponding baseline, including seed feed once; the services case falls from 1.6459×10¹⁴ kg to 1.6217×10¹⁴ kg, and diversified exports from 2.3101×10¹⁴ kg to 2.2764×10¹⁴ kg.

None of the 28 requirement paths fits the original narrow catalogue's entire-mass ceiling or the wider geometry screen's 10 percent allowance through year 200. The broad 10 percent catalogue comparison remains only a hypothetical upper stock allowance and does not qualify the assumed near-Earth delivery route. Longer component life changes replacement and freight requirements, but this bounded improvement does not resolve source qualification.
