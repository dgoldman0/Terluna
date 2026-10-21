# Crystal forms, aggregate geometry and survival

The large particles in the existing calculation are **aggregates of crystals**.
Their 15–60 cm maximum dimensions describe the span of those aggregates. The
calculation does not grow a single, sixfold symmetric crystal to that diameter.
This distinction matters both for appearance and for survival: separated facets
and branches scatter light, while contacts between component crystals transmit
uneven aerodynamic loads.

This note records a targeted review for the 21 October snow analysis. Its reader
is assumed to know the existing [column calculation](README.md). The equations
below define sensitivity calculations and identify their limits; they establish
neither a lunar snow climatology nor empirical validation. The full source
identities and actual reading extent are in [habits_sources.json](habits_sources.json).

## Where different crystal forms can grow

[Libbrecht's growth study](https://arxiv.org/pdf/1211.5555) identifies plate-like
growth near −2 °C, needles and slender columns near −5 °C, hollow columns near
−7 °C, and thin plates near −15 °C. Greater supersaturation generally permits
more elaborate branching. Applying those temperature examples to the present
301 K resort column gives:

| Temperature | Height above sea level | Pressure | Supersaturation over ice at liquid saturation | Growth context |
|---|---:|---:|---:|---|
| −15 °C | 45.53 km | 467 hPa | 15.74% | Thin plates and, with sufficient vapour, dendritic forms |
| −7 °C | 36.84 km | 563 hPa | 7.04% | Hollow columns |
| −5 °C | 34.67 km | 589 hPa | 4.98% | Slender columns and needles |
| −2 °C | 31.42 km | 630 hPa | 1.97% | Thin plates |

These are growth examples at particular temperatures, not boundaries separating
uniform layers of one habit. A falling particle experiences a sequence of
temperatures and humidities, and an aggregate retains component crystals grown
under earlier conditions. Clouds that contain supercooled droplets may also rime
the particles. A deep cloud therefore permits mixed structures and substantial
growth histories; its geometric depth alone supplies no fraction of dendrites,
needles, plates or rimed particles. The current column imposes liquid saturation
but has no prognostic droplet population or riming calculation.

`habit_ladder(column)` reproduces the table from the column, while
`humidity_regime` distinguishes sublimation from potential vapour growth. At
−15 °C, for example, 90% humidity relative to liquid water still permits ice
growth; at 80%, ice sublimates. Neither function assigns a frequency to a habit.
The temperature examples also come from terrestrial laboratory conditions;
pressure and growth history remain relevant when carrying them to this column.

## Two definition corrections and the remaining geometry uncertainty

The historical calculation uses `m = 0.0185 D^1.9` with metres and kilograms and
labels D as maximum dimension. [Hogan et al. (2012), equation 4](https://www.met.reading.ac.uk/~swrhgnrj/publications/oblate_ice.pdf)
show that Brown and Francis used the mean of two orthogonal image dimensions.
Their conversion to maximum dimension gives `m = 0.0121 Dmax^1.9`. At fixed mass
this increases the inferred maximum dimension by about 25%. Changing the law
inside the population model also changes fall speeds, collision rates and
residence times, so its final size distribution requires a fresh calculation.

The original area citation also conflates two definitions.
[Seifert et al. (2026), sections 3.2–3.3](https://arxiv.org/html/2601.10608v1)
define `q` relative to a circumscribing **ellipse**, while the fall law here needs
`Ar` relative to a circle of diameter `Dmax`. For a horizontal aspect ratio `phi`,
the conversion is `Ar = q phi`. That preprint is by Seifert and colleagues; its
aggregation algorithm builds on Leinonen's earlier work. For 1,024 monomers, its
separate fitted means give these representative geometries:

| Component crystals | Circle area ratio `q phi` |
|---|---:|
| Needles | 0.2451 |
| Plates | 0.3128 |
| Dendrites | 0.2602 |

`GEOMETRIES` preserves the historical closure as `legacy_bf95`, provides the
mass-convention correction alone as `bf95_dmax`, and combines that corrected
mass law with each representative area in three additional sensitivity cases.
The source fits cover synthetic unrimed aggregates with 2–2,048 monomers.
Multiplying separate fitted means does not reproduce their joint distribution;
holding the resulting ratio constant at every particle size is a further model
assumption. These cases consequently describe a change in area, without claiming
to simulate the full population of that named habit.

The shared mass exponent still extrapolates far beyond ordinary observed
snowflakes. Under the corrected law, a 15 cm aggregate has 0.329 g of ice and a
60 cm aggregate has 4.584 g. These are conditional mass–dimension conversions,
not measurements or evidence of mechanical survival. Depositional growth,
riming, compaction and the contact geometry can change the mass–area relation.
The existing 0.25 capacitance-to-diameter ratio is retained as an aggregate
diffusion closure supported by [Westbrook et al. (2008)](https://journals.ametsoc.org/view/journals/atsc/65/1/2007jas2315.1.xml).

## What a mechanical test must compare

A flake travelling with a steady horizontal wind does not experience that wind's
whole speed as slip. Its downward settling already produces aerodynamic loading,
however, and uneven loading can break contacts even before ambient turbulence
is considered. In [Grzegorczyk et al.'s experiments](https://acp.copernicus.org/articles/23/13505/2023/),
artificial 10 mm aggregates released without a collision shed component crystals
or branches. Their preparation and release differ from lunar snowfall, so the
experiment establishes a failure mechanism without supplying a transferable
strength for giant flakes.

`differential_load` evaluates a stated reference load, `F = Cd rho du² A / 2`,
where `du` is an imposed difference in air velocity across an aggregate. For
`rho = 0.7 kg/m³`, `Cd = 1` and the retained `Ar = 0.35`, the calculation gives:

| Aggregate span | `du = 0.01 m/s` | `du = 0.1 m/s` | `du = 0.3 m/s` |
|---|---:|---:|---:|
| 15 cm | 0.216 µN | 21.6 µN | 0.195 mN |
| 60 cm | 3.46 µN | 0.346 mN | 3.12 mN |

These forces describe a diagnostic convention: the velocity-difference dynamic
pressure acts over the whole projected area. They are **not calculated internal
bond stresses**. A full drag perturbation about the existing settling velocity
can contain terms proportional to settling speed times `du`, as well as `du²`.
Orientation, permeability and the distribution of contacts determine how that
load becomes bending or tension. `critical_velocity_difference` therefore
requires the caller to supply a force capacity; it has no default snow strength.

Residence time can change that capacity as neighbouring crystals bond together.
[Bahaloo et al. (2022)](https://ltu.diva-portal.org/smash/get/diva2%3A1628847/FULLTEXT01.pdf)
find that ice-contact forces depend on temperature, contact load, duration and
particle size. Their loaded contacts do not specify the contact histories or
connectivity of suspended flakes. Sublimating air can also erode contacts. The
long lunar transit consequently gives more time for growth, contact evolution
and destructive encounters; choosing a single breakup strength from solid ice
or deposited snow would hide the main unknown.

The next discriminating calculation needs a joint history of ice mass,
projected area, component habits, contact formation and sublimation, coupled to
the air motions experienced by falling particles. The present distribution
and turbulence sensitivity calculations can identify which sizes demand which
survival conditions. They cannot determine how often an actual lunar aggregate
meets those conditions.

## Checks

Eight focused tests passed for the size-convention conversion, legacy closure,
ellipse-to-circle area conversion, regression range, humidity distinction,
temperature-to-height mapping, load inversion and source identifiers. The table
values were evaluated with the current column and the new functions. This note
was read for consistency with those calculations; no laboratory validation,
visual inspection of published geometry figures or complete source-paper review
was performed.
