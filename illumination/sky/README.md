# Clear-sky solver

A scalar, spherical, spectral multiple-scattering solver for sky radiance and
surface irradiance under Earth and Open Moon optical profiles. It produces the
clear-sky atlases that the [month-of-light viewer](../../visualization/month-of-light/)
displays and that the immersion bakes into its sky.

- Method, profiles and limits: [METHODS.md](METHODS.md); the original algorithm
  is in [validation/legacy_v1/METHODS_v1.md](validation/legacy_v1/METHODS_v1.md).
- Sources and the Bruneton data notice: [SOURCES.md](SOURCES.md), [LICENSES.txt](LICENSES.txt).
- Check records: [validation/](validation/).

## Run

numba is required. From this folder:

```sh
python -m pip install -r requirements.txt
mkdir -p work
NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 python calculate_full_cycle.py   # writes data/*_atlas.npz
python -m unittest test_solver test_atlas -v
python check_noon.py        # six independent noon Monte Carlo spot checks
```

From the repository root, `python -m pytest illumination` runs the same tests
when numba is installed; without numba they are not collected and the pytest
header says so. `test_atlas.py` skips until the atlases have been generated.

| File | Role |
|---|---|
| `solver.py` | Optical profiles, ray paths, successive-orders transport |
| `render_data.py` | Builds the atlas product: panoramas, irradiance, spherical harmonics, cloud-altitude light |
| `calculate_full_cycle.py` | Runs the solve and atlas build for all three profiles |
| `reference_mc.py`, `check_noon.py`, `validate_spots.py` | Independent Monte Carlo spot references |
| `global_energy.py` | Whole-sphere energy audit |
