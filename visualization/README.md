# Visualization

Scientific and engineering rendering: tools whose job is to show what a model
computed, faithfully and with its evidence boundary attached. This is distinct
from the [immersion](../immersion/) lane, which builds a believable world to
experience and labels what is computed, informed by research, or artistic.

| Tool | What it shows | Model it displays |
|---|---|---|
| [atmospheric-columns](atmospheric-columns/) | Preset column soundings, Open Moon beside Earth: temperature, humidity, condensate and cloud diagnostics | [atmosphere/column](../atmosphere/column/) column set |
| [month-of-light](month-of-light/) | A fixed-viewpoint 360° view of clear-sky light across the synodic month, Earth vs Open Moon | [illumination/sky](../illumination/sky/) clear-sky atlases |
| [labs](labs/) | Interactive pages comparing the A1–A3 light-transport references with real-time approximations | [illumination/references](../illumination/references/) |
| [reference-renderer](reference-renderer/) | Offline spectral path-traced light studies of an authored coast (B1) | A3 transport from [illumination/references](../illumination/references/) |

Rules for this lane:

- Display model outputs; do not become the source of truth for them. Physics
  belongs to its domain folder, and a display approximation says it is one.
- Nothing else depends on visualization code.
- Keep each tool's build reproducible from the domain products it reads.
