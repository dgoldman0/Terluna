# Visualization

Scientific and engineering rendering: tools whose job is to show what a model
computed, faithfully and with its evidence boundary attached. This is distinct
from the [immersion](../immersion/) lane, which builds a believable world to
experience and labels what is computed, informed by research, or artistic.

| Tool | What it shows | Model it displays |
|---|---|---|
| [month-of-light](month-of-light/) | A fixed-viewpoint 360° view of clear-sky light across the synodic month, Earth vs Open Moon | [illumination/sky](../illumination/sky/) clear-sky atlases |

Rules for this lane:

- Display model outputs; do not become the source of truth for them. Physics
  belongs to its domain folder, and a display approximation says it is one.
- Nothing else depends on visualization code.
- Keep each tool's build reproducible from the domain products it reads.
