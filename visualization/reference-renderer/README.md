# Reference scene renderer (B1)

An offline spectral path tracer that renders full-frame light studies of an
authored 48-km coastal heightfield under clear and cloudy Open Moon skies. It
uses the coupled gas/cloud transport from
[illumination/references](../../illumination/references/) (A3) unchanged, adds a
finite solar disk, Lambertian terrain and a Fresnel sea, and writes coarse
radiometric images (160 × 96 frames, sky hemispheres and irradiance maps) with
their sampling errors. Method and limits: [B1_METHODS.md](B1_METHODS.md); record of
what was run: [B1_VALIDATION.json](B1_VALIDATION.json).

This is scientific rendering: its images are radiance estimates with stated
uncertainty, not the experience. The experience lives in [immersion](../../immersion/).

```sh
python -m pip install -r requirements.txt
python tools/b1/run_checkpoint.py        # full checkpoint sequence (long; numba)
python tools/b1/test_b1.py               # boundary and transport checks
python build_b1_painter.py               # Open_Moon_B1_Scene_Painter.html
```

Inputs from the reference program keep their recorded relative names
(`data/a1/...`, `tools/a3/...`); `tools/b1/scene_transport.py` resolves them to
`illumination/references/`. `validate_b1.py` re-verifies the published checkpoint,
including a pinned build of the experience; after the repository reorganization
that gate reports the changed files.
