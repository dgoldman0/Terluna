# Open Moon — Shoreline immersion prototype

This directory contains the **normal editable source tree** for the current shoreline immersion checkpoint. The source is committed as ordinary files; it is not stored behind an archive wrapper.

## Status

This is a WIP preservation checkpoint, not a visual-quality milestone. The browser render has known serious weaknesses:

- distant terrain can read as floating/disconnected ribbons rather than coherent shoreline relief;
- water reflections are coarse and visibly aliased;
- vegetation, rocks, ground materials, and placement remain placeholder-like;
- lighting/exposure can become muddy enough to obscure the structure of the scene;
- the integrated render is not yet calibrated to a lifelike visual standard.

The next milestone should make one shoreline view convincing at noon, then verify the same view under low Sun and twilight before expanding environmental complexity.

## Layout

- `src/` — editable runtime code
- `data/atmosphere.json` — packed inherited full-cycle optical data used by this checkpoint
- `index.template.html` — editable HTML shell
- `Open_Moon_Shoreline.html` — generated first-open build for this checkpoint
- `scenario.json` — authored local scene/weather inputs
- `build.py` / `pack_atmosphere.py` — build/data tools
- `VALIDATION.json`, `SOURCES.md`, `LICENSES.txt` — evidence boundaries and dependency records

The preceding full-cycle viewer remains unchanged.

The generated build retained here has SHA-256 `8cecea048a67040c532f1660b4b071dd0fe1db1927dfc063aba862e45ca3b218`. It loads pinned Three.js r180 on first opening. The current visual weaknesses above remain part of this checkpoint and should not be interpreted as properties of the proposed Open Moon environment.
