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
| [reference-renderer](reference-renderer/) | Offline spectral path-traced light studies of an authored coast (B1) | A3 transport from [illumination/references](../illumination/references/) |
| [waves](waves/) | The wave studies' figures: Earth and lunar gravity over a strip, the Smythii–Marginis basin and its lunar cycle, coastal nests and shore exposure, breaking and run-up, the nearside sea's waves and coasts, a month at four shores, and every sea's monthly tide | [climate/waves](../climate/waves/) products and the [geography tide product](../geography/README.md#the-monthly-tide) |
| [optical-comfort](optical-comfort/) | Surface luminance, eye illumination and shaded light against Earth; the directional scenes' gaze, finite scenes, workplane light and tread–riser contrast | The [optical comfort study](../research/studies/optical_comfort/) |
| [cloud-twilight](cloud-twilight/) | Evening clouds: their illumination and occurrence, the regional survey's maps, computed cloud views in colour and brightness, a cloud system through 72 hours and a distant cloud in the deep evening, with a local gallery | The [cloud twilight study](../research/studies/cloud_twilight/) |
| [sea-appearance](sea-appearance/) | How the seas look: the light through a month at six coasts, from the Sun's long twilight and the Earth, where the Earth stands over every sea, the sea's slopes at both gravities and through the month, the waters' colours, and panoramas of sky and sea by regime at four coasts | The [sea appearance study](../research/studies/sea_appearance/) |

Rules for this lane:

- Display model outputs; do not become the source of truth for them. Physics
  belongs to its domain folder, and a display approximation says it is one.
- Nothing else depends on visualization code.
- Keep each tool's build reproducible from the domain products it reads.
