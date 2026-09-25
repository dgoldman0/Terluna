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

## Spectral transmission product

[spectra/transmission.py](spectra/transmission.py) writes
[spectra/shield_transmission.json](spectra/shield_transmission.json) (schema
`terluna.protection.shield-transmission/1`). It holds T(λ) from 200 nm to 5 µm
for the stored titania–silica design and for idealised edge filters. The design
is evaluated with `model.py`'s own thin-film code; the edge filters are
scenarios, not designs. The atmosphere domain reads this product to filter the
sunlight in its climate and photochemistry calculations. Regenerate it with
`python -m protection.spectra.transmission` after restoring the inputs.

## Short-wave transmission (2026-09-25)

[spectra/short_wave.py](spectra/short_wave.py) writes
[spectra/stack_short_wave.json](spectra/stack_short_wave.json) (schema
`terluna.protection.stack-short-wave/1`): the stored design's T(λ) from 0.1 to
210 nm, with the design's own layers and thin-film code. The optical constants
are published:
- CXRO atomic scattering factors (Henke, Gullikson and Davis 1993) for both
  oxides below 24.8 nm;
- fused silica from 24.8 nm (Franta et al. 2016, pinned in
  [inputs.json](inputs.json));
- titania from 120 nm (the design's Siefke data).

Between 24.8 and 120 nm no measured titania set was found, so the CXRO
estimate stands in; the 10 µm of silica alone is opaque there (below 10⁻¹⁴⁶),
so this does not affect the result. The rebuilt layer sequence reproduces
`model.optical_stack` to 10⁻¹², and the X-ray absorption reproduces the design's
[xray_absorption.csv](results/xray_absorption.csv).

The film passes hard X-rays only:

| Wavelength | Transmission |
|---|---|
| 0.1 nm | 96% |
| 0.5 nm | 9% |
| 1 nm | 0.5% |
| 1.5 nm | 2×10⁻⁷ |
| Anything beyond 2.5 nm | 4×10⁻⁸ at most |
| 5–200 nm | below 10⁻²⁰ |

The atmosphere domain's escape calculation reads this. The film lets at most
10⁻⁹ W/m² of heat into the upper air for the quiet Sun and 10⁻⁷ W/m² at solar
maximum, so the exobase stays within 1.5 K of its base. How much extreme
ultraviolet reaches the Moon is therefore set by light that bypasses the film:
gaps, pinholes and edges of the aperture, and off-normal incidence. None of
these is modelled; they are design parameters. Regenerate it with
`python -m protection.spectra.short_wave` after restoring the inputs.

## Condition and remaining questions

- Optical layers use measured constituent constants plus effective-medium and normal-incidence approximations. This is not a measured complete coating, irradiated lifetime test, or gap-free aperture.
- The holding-force calculation is an instantaneous circular-phase sweep, not full ephemeris propagation or a passive lunar orbit.
- Specific electrical power, material properties, fuel buffer, exhaust cant and storage capability are assumed design parameters. Clean plume operation remains unestablished.
- Magnetic moment, pressure-balance and rigidity estimates are component diagnostics. Global plasma performance, reconnection, trapped particles and human radiation dose remain unmodelled.
- Thermospheric chemistry, lower climate, water loss and species escape remain to be coupled. The original model does not establish that its filter produces a 250 K exobase. The film's own short-wave transmission is now computed (above), and the atmosphere domain finds that it holds the exobase at 142–193 K. Light bypassing the film through the aperture is not yet specified.
- Hardware replacement may be much smaller than construction throughput while accumulated propellant/resource use remains substantial over geological time.

The preserved results include successful and failed propulsion closures. Large sampled optical/phase grids and rendered report files are omitted from this curated import and recorded in [provenance](../research/provenance.json). Their numerical tables can be regenerated from the original code once inputs are restored. The older [shield geometry table](reference/legacy_shield_geometry.csv) remains explicitly separate from this later design.

[Checks](../research/checks.json) report precisely which implementation tests were run, without granting environmental or engineering validation.
