# Cultural production and lunar export income

Humanities economics studies how people could earn income through lunar practices,
performances and cultural work, and how that income reaches the settlement. It
uses the domain's physical research as evidence about possible activities and
keeps commercial assumptions separate from measured demand. This is research
support for the human paper, with no change to the five manuscript roles.

| Study | Contribution | Evidence state |
|---|---|---|
| [Professional lunar sports](sports.md) | Resident clubs, a season, audience acquisition, rights, digital adaptations, production, trade and project cash | Reproducible conditional accounts; no selected sport, qualified venue, observed lunar demand or financed league |

The [world-economy study](../../research/studies/world_economy/README.md) owns the
cross-domain planetary accounts. These studies identify cultural activity within
those accounts; their turnover is never added wholesale to existing world product.
The earlier [venue finance study](../../research/studies/world_economy/sports_finance.md)
remains a separate, unchanged scenario. Physical sport research remains in
[sport](../sport/README.md), and art and motion research in [dance](../dance/README.md).

All modeled financial amounts use **U**, the project's hypothetical real
accounting unit. Published market comparisons retain their original US dollars
and observation years. There is no assumed exchange rate between them.

Run from the repository root:

```sh
python3 -m humanities.economics.sports
python3 -m pytest -q humanities/tests/test_cultural_economics.py
```

Normal runs write to ignored `research/runs/cultural_economics/`. The reviewed
product is [results/sports.json](results/sports.json); updating it requires an
explicit `--output humanities/economics/results/sports.json`. Its versioned schema,
producer hashes, input identities, units and reading rule accompany the results.
The source register records targeted reading and unresolved retrievals;
[checks.json](checks.json) records verification separately from commercial evidence.
