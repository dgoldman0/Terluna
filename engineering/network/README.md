# Industrial network model

The start of the integrated source-to-use model the engineering domain lists as
next work. It holds two pieces so far.

[model.py](model.py) fixes the vocabulary: node kinds (source extraction and
power, factories, launch, relay depots, cislunar capture, lunar orbital depots,
construction, protection service, storage and release, recycling), commodities,
link properties and node constraints. `Network.validate()` checks a network
built from them, and it enforces the author's traffic-safety rule that passive
failures miss planets. The vocabulary follows section 10 of the reconstructed
19 September architecture
([reference/industrial_architecture](../reference/industrial_architecture/README.md)),
which supplies names and relationships only.

[ledger.py](ledger.py) reads the protection module catalogue
([protection/modules/catalogue.json](../../protection/modules/catalogue.json)) as
demand and writes
[results/protection_supply_ledger.json](results/protection_supply_ledger.json)
(schema `terluna.engineering.protection-supply-ledger/1`): gross replacement,
fresh material after recovery, irreversible propellant and power, by commodity,
for replacement every 20, 100 or 1,000 years at 99% or 99.9% recovery, with
totals over 10⁹ years.

```sh
python -m protection.modules.catalogue       # the demand
python -m engineering.network.ledger          # the ledger
python -m pytest engineering
```

## What the ledger shows for the September reference

Replacing the film every century takes 1,547 kg/s, all the screen's dry hardware
1,883 kg/s, and the magnets 95,000–951,000 kg/s at the report's planning
allowance; at 99.9% recovery the fresh material is a thousandth of that. The
propellant is the large irreversible flow: 2.8×10⁵ kg/s, 8.8×10²¹ kg over 10⁹
years, 3.5 times the highest loss of air the loss response finds with no shield
at all. That is requirement S6 failing
([research/studies/protection_architecture](../../research/studies/protection_architecture/requirements.md)).
These figures reproduce the September report's own and inherit all its
assumptions.

## Next

Source, transport and capture nodes with their capacities; the delivery routes
that carry the air, water and replacement units; and the demand of a redesigned
protection system once the holding trade has chosen one.
