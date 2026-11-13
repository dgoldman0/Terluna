# U monetary-pressure verification

13 November 2026. The study is a conditional numerical experiment for the planned off-Earth economic zone. Verification establishes arithmetic, input identity and reproducibility; it does not establish observed money demand, a qualified supply response, reserve adequacy, an adopted currency constitution or a floating exchange-rate forecast.

The project change selector from the current main checkout was loaded with its ROOT set to the isolated humanities worktree, through `/tmp/check_u_pressure.py`. Its `--plan` selected the layer check and all three monetary-system test files. This uses the current repository selection procedure while preserving the humanities branch's existing Makefile.

Executed checks:

- `python /tmp/check_u_pressure.py`: **53 tests passed in 3.53 seconds**. Its selected pytest command covered `test_monetary_system.py`, `test_pressure.py` and `test_trade_transition.py`; the layer check inspected **1,657 files with zero violations**. The initial focused pressure run passed 12 tests before the final reserve/dealer refinements; the 53-test run covers those refinements.
- The new product test regenerated `results/pressure.json` exactly, reproduced all 27 underlying financial-path summaries, and verified 108 records with 73 financially continuing paths and 35 failed-path exclusions.
- Numerical controls check local/outside stock conservation, neutral issuance, persistent monetary shortfalls, outside-dealer redistribution, monetary retirement, exhausted domestic balances, single allocation of reserve cash, invariance to common unit scaling, import-quote direction, and compounding of annual changes into endpoint indices.
- Byte comparisons against the parent commit confirmed unchanged strategy, earlier-price, institutional-reference and trade-transition products. Their respective SHA-256 values are `685a27c5b6820adf4b1e188e19906a1c747f2a34ef8b0da8f97c8326b62ed4dc`, `6a67dd4da7e1e2f3026cecf0610ae0230b41a5c6c32221a717913b8eba4e011e`, `447ac1d96d70b60d1385b1a82c1394344ff1f8eea3be88de93d3cb498646c10d` and `34eabb915ab5bc6247e7df23ef932d12e391a992f409ffaf76b20103bb4a0a9f`.

Two bounded read-only reviews examined monetary transmission and the financial inputs. They identified the distinction between acquisition not supplied by issuance and genuinely unfilled exchange orders; the requirement to fund outside dealers' opening inventory; the private external-asset assumption in incomplete monetary retirement; and endpoint payouts versus the prior study's decade-average flows. The implementation and write-up now retain those distinctions. The reversal also reports automatic adjustment of domestic transaction balances alongside retirement of returned outside balances.

The finished study was read for its definitions, assumptions, equations, comparative tables and interpretation. Numeric tables were compared with the saved product. Primary literature received targeted reading recorded in `pressure_sources.json`, rather than a full literature review. No figures or PDF were produced. The full repository suite was not run because the change selector bounded the affected work to the monetary study and structural checks. Monetary inputs, supplier currencies, long-term saving, constituent-specific clearing and feedback into export demand and investment remain research requirements.

To reproduce the product from the repository root:

```sh
python -m research.studies.monetary_system.pressure --output research/runs/monetary_system/pressure.json
python -m pytest research/studies/monetary_system/tests
```
