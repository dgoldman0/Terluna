# UV transmission of the swarm

How much sunlight below 175 nm the optical shield lets through once it is a
formation of replaceable cells. The intact titania–silica film passes at most
1.6×10⁻⁴³ of it ([spectra/stack_short_wave.json](../spectra/stack_short_wave.json),
10–175 nm), so the UV transmission is the share of the aperture left uncovered,
averaged over time and area. [model.py](model.py) adds four terms: seams,
micrometeoroid holes, missing cells and manufacturing defects. It does this for
three levels of design choices; the method is in its docstring.

```sh
python -m protection.transmission.model   # writes results/uv_transmission.json in a few seconds
python -m pytest protection/tests/test_transmission.py
```

The product is [results/uv_transmission.json](results/uv_transmission.json)
(schema `terluna.protection.uv-transmission/1`). The protection design study
turns it into a loss rate, a cycle time and a protected radius
([research/studies/protection_architecture](../../research/studies/protection_architecture/README.md)).

## Results

| Level | Choices | UV transmission |
|---|---|---|
| Tight | neighbours within 1 m, 50 m overlap; cells patched or replaced every 10 years; 1% of cells fail a year, each covered within a day | 2.9–4.0×10⁻⁵ |
| Standard | within 3 m, 50 m overlap; every 20 years; 1% a year, covered within a week | 1.9–2.2×10⁻⁴ |
| Relaxed | within 10 m, 20 m overlap; every 100 years; 3% a year, covered within a month; defects 10⁻⁵ | 2.5–2.6×10⁻³ |

The ranges span micrometeoroid holes 1–4 times the particle diameter, and each
level includes manufacturing defects of 10⁻⁶ (10⁻⁵ at the relaxed level) as a
specification.

**Seams vanish once neighbours overlap.** Cells that overlap by several times
their relative position error leave nothing open; a 50 m overlap costs 1% more
film. Butt-jointed cells with a 3 m error would leave 2.4×10⁻⁴ open.

**Micrometeoroid holes are slow.** The Grün et al. (1985) flux holes
1.5×10⁻⁷ of the film a year if each hole is the particle's own size, mostly
through particles near 60 µm. Patching or replacing cells every 10–20 years
keeps holes at 0.7–23×10⁻⁶; left 1,000 years they would reach 7×10⁻⁵–1.2×10⁻³.

**Missing cells set the level.** A failed cell opens 10⁻⁶ of the aperture until
it is covered, so the average is the failure rate times the time to cover.
Covering within a day needs spare cells or overlays kept in the formation. The
light through a missing cell spreads over the Sun's blur: at the Moon it lights
a patch 368 km in radius at 2.4×10⁻⁴ of full sunlight, so a single gap makes
no concentrated leak.

## Limits

The design levels are assumptions to test. No formation-control design, failure
record or hypervelocity test of the oxide film exists yet. The hole size comes
from aluminium and Teflon targets, known here from their abstracts; a brittle
oxide film may crack beyond the hole. The Grün flux is the interplanetary
average at 1 AU, without meteor showers or focusing by Earth and the Moon.
Contamination and darkening of the film, exhaust plumes crossing the aperture,
and light scattered by the film's edges and framing are not included.

Sources: Grün, Zook, Fechtig and Giese (1985), *Icarus* 62, 244–272, in the form
given by [SPENVIS](https://www.spenvis.oma.be/help/background/metdeb/metdeb.html);
Hörz and colleagues' penetration experiments in
[aluminium](https://www.researchgate.net/publication/223543633_Dimensionality_scaled_penetration_experiments_Aluminum_targets_and_glass_projectiles_50_mm_to_32_MM_in_diameter)
and [Teflon](https://www.sciencedirect.com/science/article/abs/pii/0734743X9599867Q) targets.
