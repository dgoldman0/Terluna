# Atmospheric protection

`model.py`, `verify.py` and the compact result tables are verbatim imports from the September 2026 protection package. `results/results.json` retains selected original constants, optical summaries, reference/storage cases and input hashes for the verifier; duplicate table arrays and unrelated summaries are omitted, with values unchanged. They calculate candidate optics, aperture geometry, lunar-phase holding forces, propulsion/power/storage feedback, magnetic structures, particle-rigidity diagnostics and renewal budgets.

The protected atmosphere and transport comparison are inherited from the [September feasibility baseline](../research/baselines/feasibility/), with fixed inventory constants in the original code. Their use is not an independent atmospheric validation.

## Inputs and reproduction

Five external optical/solar input files are pinned in [inputs.json](inputs.json). They are deliberately outside this source commit; data attribution, retrieval URLs, sizes and exact checksums are retained. The original files are available in the recovered project ZIP in the current session. Restore them from that package or fetch the recorded upstream bytes:

```sh
python protection/fetch_inputs.py --archive /path/to/Lunar_Protection_Model.zip
# Alternatively, when network access is available:
python protection/fetch_inputs.py --download
python protection/verify.py
```

The verifier imports `model.py`, which immediately loads the TiO2 table; restoring inputs is required even for that original verifier. A SHA-256 mismatch stops retrieval. Upstream changes require an explicit reviewed update; hashes are never silently rewritten.

After installing `research/requirements.txt`, run `python protection/model.py` to regenerate the optical optimization and all reference outputs. It writes into `protection/results/`; preserve the selected reference snapshot with Git or run in a copied working directory. The fixed optimization seed aids reproducibility but floating-point/optimizer versions can change results. `python research/check.py` performs the safe checks in a temporary copy.

## Condition and remaining questions

- Optical layers use measured constituent constants plus effective-medium and normal-incidence approximations. This is not a measured complete coating, irradiated lifetime test, or gap-free aperture.
- The holding-force calculation is an instantaneous circular-phase sweep, not full ephemeris propagation or a passive lunar orbit.
- Specific electrical power, material properties, fuel buffer, exhaust cant and storage capability are assumed design parameters. Clean plume operation remains unestablished.
- Magnetic moment, pressure-balance and rigidity estimates are component diagnostics. Global plasma performance, reconnection, trapped particles and human radiation dose remain unmodelled.
- Spectrum, thermospheric chemistry, lower climate, water loss and species escape remain to be coupled. The original model does not establish that its filter produces a 250 K exobase.
- Hardware replacement may be much smaller than construction throughput while accumulated propellant/resource use remains substantial over geological time.

The preserved results include successful and failed propulsion closures. Large sampled optical/phase grids and rendered report files are omitted from this curated import and recorded in [provenance](../research/provenance.json). Their numerical tables can be regenerated from the original code once inputs are restored. The older [shield geometry table](reference/legacy_shield_geometry.csv) remains explicitly separate from this later design.

[Checks](../research/checks.json) report precisely which implementation tests were run, without granting environmental or engineering validation.
