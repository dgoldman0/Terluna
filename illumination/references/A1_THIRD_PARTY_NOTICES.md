# A1 third-party notices

The 48-value solar/ozone tables, selected CIE 1931 values and XYZ-to-linear-sRGB matrix in `data/a1/spectral-inputs.json` were transcribed from Eric Bruneton's precomputed atmospheric scattering reference implementation. Source provenance is recorded in that JSON file and A1_METHODS.md. The originating source is distributed under BSD-3-Clause. Its copyright, conditions and disclaimer are retained verbatim in `data/a1/BRUNETON_LICENSE.txt` (upstream LICENSE blob `bfad612c567c7b3a2772f614b609fecb16385243`).

A1's transport code was independently implemented. PBRT is cited for algorithms; its code is outside this distribution. NumPy, SciPy, Numba and Playwright are runtime dependencies listed separately, with their own upstream licenses. Existing baseline Three.js and asset notices remain unchanged.
