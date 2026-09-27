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
| [atlas](atlas/) | Map sheets of the Open Moon at its selected water share: seas, islands, heritage sites and scientific targets; a 3D globe of how it looks from space and from Earth, with its air, clouds and a guessed plant cover | [geography atlas](../geography/README.md#atlas-at-the-selected-water-share), the [conservation study](../research/studies/conservation/), the [engine atmosphere and sky atlas](../illumination/sky/) and the [GCM climatology](../climate/gcm/) |
| [equator-weather](equator-weather/) | The Open Moon's equator through the month-long day from a cloud-resolving model: weather by local time over land and sea, a storm system in section, where and when it rained, the circulation that follows the Sun, winds by height and the flight band | [climate/crm](../climate/crm/) equatorial ring products |
| [summit-port](summit-port/) | The central port on the summit as an architectural model: leg towers, wings with berths and gardens, gates and the crown, from the whole tower down to human scale, and the programme by zone | The [summit tower study](../research/studies/summit_tower/)'s port (frame, bands, zones, programme); the architecture on it is a proposal |
| [reference-renderer](reference-renderer/) | Offline spectral path-traced light studies of an authored coast (B1) | A3 transport from [illumination/references](../illumination/references/) |

Rules for this lane:

- Display model outputs; do not become the source of truth for them. Physics
  belongs to its domain folder, and a display approximation says it is one.
- Nothing else depends on visualization code.
- Keep each tool's build reproducible from the domain products it reads.
