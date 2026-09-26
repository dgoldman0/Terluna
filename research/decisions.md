# Decisions

The author's decisions about the Open Moon (Terluna), with where each is recorded
and whether it stands. A later decision replaces an earlier one, and the
replaced entry names what replaced it. Modelling choices, reference values and
results appear here only where the author adopted them.

Sources:

- **July**: the recovered canon of the July 2026 knowledge bundle
  ([archive/july_2026_knowledge_bundle/](../archive/july_2026_knowledge_bundle/README.md);
  its file `01_RECOVERED_CANON_AND_CORRECTIONS.md` and `data/decision_log.csv`).
- **Sep 4**: the user requirements of the September feasibility report
  ([baselines/feasibility/report.md](baselines/feasibility/report.md), section 1).
- **Sep 9**: the September protection report
  ([protection/report.md](../protection/report.md), section 3).
- The conservation study keeps its own decisions, D1–D11, in
  [its README](studies/conservation/README.md#4-decisions); they are not repeated here.

## The world

| Decision | Source | Status |
|---|---|---|
| The world is Earth's own Moon, transformed in place, and the project is called Terluna. | July | Stands. The paper ensemble's working title is *Constructing and Sustaining an Open Moon*. |
| The Moon keeps its orbit, gravity, synchronous rotation and 29.53-day solar cycle. Life and culture adapt to the long day and night, storing resources across them. | July | Stands |
| Construction within 500 years of 2026, and operation for at least 10⁹ years. | Sep 4 | Stands |
| An open atmosphere with no pressure dome over the Moon. | Sep 4 | Stands. Local enclosed places, such as the undersea Tranquility community (conservation D6), are separate. |

## Light and time

| Decision | Source | Status |
|---|---|---|
| Practical dawn and dusk each last about five hours, inside the 29.53-day cycle. This replaced an earlier 30-hour figure. At the equator five hours is about 2.54° of solar elevation. | July | Stands; settled on 2026-09-26 against the ensemble plan's nominal 12.5-hour sunset, a design case. It still needs a photometric definition for the actual atmosphere. |

## Atmosphere

| Decision | Source | Status |
|---|---|---|
| Surface pressure 1.2 atm (121,590 Pa). | Sep 4 | Stands; the design pressure of the climate runs. The July bundle's 80 kPa was its own reference value. |
| An exobase near 250 K is preferred, with 260 K a candidate limit. | Sep 4 | Stands. The solved upper air behind the titania film gives 142–193 K before the sky's Lyman-alpha glow is counted ([middle atmosphere](../atmosphere/middle_atmosphere/README.md)). |
| Water, nitrogen and oxygen come from the Solar-System-wide resource operation. | Conservation D4 | Stands; replaced the September baseline's oxygen from lunar rock. |
| Atmospheric layers take descriptive names, never L1–L5, which stay reserved for Lagrange points. | July | Stands |
| "Gravity wave" means an atmospheric buoyancy wave; write "atmospheric gravity wave" on first use. | July | Stands |
| A band at 35–45 km was placed as the top of the main troposphere, provisionally, pending a climate model. | July | Settled on 2026-09-26: the band stays as a named flight band, and its physics comes from the models. The 1-D column puts the tropopause near 86 km (288 K surface); in the 5% climate the band holds air at −10 to −20 °C and 0.58–0.48 atm ([climate/gcm](../climate/gcm/README.md)). |

## Water and climate

| Decision | Source | Status |
|---|---|---|
| Standing water covers 28% of the surface. | [shared/scenarios/water.json](../shared/scenarios/water.json); conservation D10 | Stands; selected on 2026-09-25 within the author's 25–35% range. |
| The shield passes 5% less sunlight than the titania stack at every wavelength, and the Moon settles near 294 K. | [climate/gcm](../climate/gcm/README.md) | Stands; chosen on 2026-09-26. |

## Protection

| Decision | Source | Status |
|---|---|---|
| The protection architecture: a solar-filter complex with an industrial hub; magnetic protection; protection that grows around occupied destinations and Earth–Moon traffic corridors, secondary to protecting the lunar atmosphere and surface; Earth–Sun L1/L2 hubs for freight, power and industry. | Sep 9, carrying the original 2025 architecture | Stands |
| No superconducting planetary ring. | Sep 4 | Stands |
| A surface radiation dose of at most 0.027 mSv/day. | Sep 4 | Stands |

## Life and people

| Decision | Source | Status |
|---|---|---|
| Aerial and gliding life is a major feature of the world. | July | Stands |
| High-altitude life is sparse: microbes or engineered radiation-hard films. | July | Stands |
| Long-endurance high-altitude platforms serve farside astronomy and nearside earthshine photometry. | July | Stands |
