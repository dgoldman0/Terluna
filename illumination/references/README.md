# Sky and cloud light-transport references (A1–A3)

Reference radiative-transfer calculations for the Open Moon atmosphere, and the
measurements that compare the immersion's real-time sky and cloud lighting against
them. This work was developed inside the shoreline experience in September 2026
(checkpoints A1–A3) and moved here because it is illumination science: it computes
how light moves through the air and clouds, independent of any one renderer.

| Step | What it computes | Main result |
|---|---|---|
| A1 ([methods](A1_METHODS.md), [record](A1_VALIDATION.json)) | A shared thermodynamic and optical column (the atmosphere domain's [column model](../../atmosphere/column/) continued hydrostatically to 600 km) exported as 9 profiles; three independent optical-column integrators; a 48-band single-scattering spectral sky; a GPU read-back of the experience's cloud field and a Monte Carlo cloud multiple-scattering reference | Column integrators agree to better than 1e-6; the shared column's vertical air path is about 11.6% longer than the exponential proxy used by the [clear-sky solver](../sky/) |
| A2 ([methods](A2_METHODS.md), [record](A2_VALIDATION.json)) | Successive-orders molecular multiple scattering through the A1 profile (48 bands, up to 56 orders in the blue); cloud voxel-grid and crop convergence | On seven sampled rays, multiple scattering raises luminance 1.9–5.4× over first order (a per-ray ratio, not a global multiplier) |
| A3 ([methods](A3_METHODS.md), [record](A3_VALIDATION.json)) | Coupled gas-plus-cloud backward Monte Carlo transport (about 50 million histories) | The experience's real-time cloud-lighting expressions come out about 39% low on the dense-cloud benchmark rays; A2's molecular field is about 8% low at 400 nm |

The limitation register that motivated the program is [CLOUD_ACCURACY.md](CLOUD_ACCURACY.md);
its plan is [CLOUD_NEXT_ACCURACY_MILESTONE.md](CLOUD_NEXT_ACCURACY_MILESTONE.md).
The offline scene path tracer built on A3 (B1) is scientific rendering and lives in
[visualization/reference-renderer](../../visualization/reference-renderer/); the
interactive lab pages that display these results are in
[visualization/labs](../../visualization/labs/).

## Layout

- `tools/a1`, `tools/a2`, `tools/a3`: Python reference programs (numba) and the
  A1 profile export (`tools/a1/export_profiles.cjs`, Node).
- `src/profile-optics.js`: independent JavaScript path integration of the shared
  profile. `src/frozen-cloud-field.js` and `tests/*_frozen_probe.js`: probes that
  read the experience's cloud field back from the GPU for comparison.
- `data/a1/spectral-inputs.json`: the 48 Bruneton solar/ozone/CIE bins
  ([notice](data/a1/BRUNETON_LICENSE.txt)). Generated outputs go to
  `data/a*/generated/` (ignored).
- `validation/a1..a3`: compact records of what was run.

## Run

```sh
python -m pip install -r requirements.txt       # numpy, scipy, numba, playwright
node tools/a1/export_profiles.cjs                 # writes data/a1/generated/profiles.json
node --test tests/*.test.cjs                      # profile, optics and probe unit tests
```

The capture tools (`tools/a1/capture_field.py`, `tools/a2/capture_cloud_convergence.py`)
measure a built copy of the shoreline experience in headless Chromium.

`validate_a1.py`, `validate_a2.py` and `validate_a3.py` re-verify a checkpoint
against the exact files it was published with, including a pinned build of the
experience. They document how each `A*_VALIDATION.json` was produced. After the
repository reorganization their baseline-hash gates report the changed files; to
reproduce a checkpoint exactly, check out the publication commit named in its record.
