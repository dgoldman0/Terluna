# Renderer lab (retired)

A Three.js r186 WebGPU/TSL port of the shoreline viewer, kept as an engine
experiment until it showed value over the WebGL path. It never did, and it froze
when its Sun slider was dragged (each input rebuilt the environment lighting and
rendered a full frame). It was retired on 2026-09-23 when the explorable world moved
to Unreal Engine 5 (github.com/dgoldman0/terluna-game).

It no longer builds from here: its imports point at the immersion engine by
relative path. The last commit where it ran in place is `7f997cf`
(`immersion/experiences/renderer-lab/`). `ROADMAP.md` is its evaluation plan.
