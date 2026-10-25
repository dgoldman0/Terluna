# A breakup law for snowflakes: what Earth's measurements support

A 45-minute search on 25 October 2026 found no law, measured or fitted on Earth, for how a whole snowflake breaks while
it falls, either in collisions with other flakes and crystals or under drag and turbulence. The measurements that exist
concern flakes struck by dense ice or blown against a snow surface, and they count fragments chipped off or record the
impact speed at which flakes come apart there. The load analogy therefore remains the estimate of where lunar flakes
break ([survival](turbulence_survival.md)), with the limits set out at the end of this note. No breakup patch was added
to McSnow, and the aggregation runs of 24 October stand as they were.

This note serves the snow work's next steps. It supplements the 21 October review of contact strength and fragmentation
([strength literature](strength_literature.md)), which reached the same conclusion for decimetre aggregates from the
contact and sintering side. The sources, the quantities taken from each and how each was obtained are recorded in
[breakup_sources.json](breakup_sources.json).

## The closest Earth measurements

| Study | What was measured | Quantities | What it constrains |
|---|---|---|---|
| Sato, Kosugi, Mochizuki and Nemoto (2006), cold wind tunnel and field | Artificial aggregates of dendrites, 7 and 20 mm, at −7 °C, falling onto soft and hard snow surfaces in wind; natural flakes up to about 10 mm at −0.1 to −1.4 °C in the field | Normal impact 0.8–1.0 m/s; horizontal impact speed equal to the wind at 1 m. No shear fracture below about 2 m/s, complete decomposition into crystals above about 5 m/s, partial fracture between. The compression of flakes by the normal impact itself was seen and not scored | The impact speed at which loose aggregates fail against a rigid surface |
| Grzegorczyk and co-authors (2023), laboratory | Aggregates of dendrites, 10 ± 2 mm, at −15 °C, struck by 4 mm graupel | Collision kinetic energy 10⁻⁷–10⁻⁵ J; 50 to more than 400 fragments, size mode 0.4 mm; central hits broke the flake completely. Fit to the Phillips form: A = 2 × 10⁷ m⁻², C = 5.9 × 10⁸ J⁻¹, γ = 0.78, as printed on the article page | Fragments from dense-ice impacts, and the energy at which a centimetre flake is destroyed by one |
| Gautam and co-authors (2024), field chamber | Natural snow striking ice spheres, northern Sweden | About 5–40 fragments per collision for dendritic snow of 3–12 mm, 1–10 for non-dendritic snow of 1–3 mm (abstract only) | Fragment number from snow–graupel contacts |
| Phillips, Yano and Khain (2017), theory fitted to earlier experiments | Collisions among graupel, hail, snow and crystals | Fragment number N = αA[1 − exp(−(C·K/αA)^γ)] for collision energy K (formulation only) | Chipping of branches; tested in this project on 20 October, where it changed lunar flake sizes by 10–25% without capping them ([README](README.md)) |
| Lo and Passarelli (1982), aircraft spiral descents | Height evolution of snow size spectra in winter storms | Spectra pass through deposition, aggregation and breakup stages; breakup limits the number of large particles and holds the spectral slope at a limiting value (abstract only) | That breakup acts in real snow, with no rate or strength law |

Recent work points the same way. Kötsche and co-authors (2026) infer ice-collisional fragmentation inside turbulent
layers of Colorado orographic snow from a second radar mode of small anisotropic ice, using a dissipation-rate threshold
of 10⁻³ m²/s³; they measure no change in the large aggregates. Von Terzi and co-authors (2026, preprint) need
collisional fragmentation in McSnow to reproduce radar signatures in the dendritic growth zone, again as a source of
small ice. BinMod1D (2026) can run snow aggregation with collisional breakup, with fragment distributions borrowed from
rain or chosen by the user and no snow calibration. Guala and Hong (2015) propose that clustering at a critical Stokes
number limits how large aggregates grow in turbulence, a growth limit that involves no breakup. The turbulent breakup of
colloidal and marine aggregates has a developed theory, in which an aggregate fails where the flow's stress exceeds a
strength that falls with size; its strengths come from those materials and none has been measured for snow.

## Collisions aloft are gentler than any measured breakage

In free fall a flake meets crystals and other flakes at the difference of their fall speeds. On the Moon a 0.2 mm
particle falls at 0.12 m/s and a 5 cm flake at 0.64 m/s; on Earth the same sizes fall at 0.54 and 1.50 m/s
([README](README.md)). The largest closing speed from settling is therefore about 0.5 m/s on the Moon and 1 m/s on
Earth, below the 2 m/s at which Sato's aggregates began to fracture against a snow surface. The energies are small as
well: the 3 mm, 0.081 mg crystal of the lunar column at 3 mm/h striking a 10 cm flake carries about 10⁻⁸ J, and a 2 cm
and a 5 cm lunar flake, of 11 and 62 mg by Brown and Francis's mass law, close at 0.034 m/s with about 5 × 10⁻⁹ J.
Turbulent velocity differences of (εD)^⅓ raise flake-on-flake energies to about 3 × 10⁻⁸ J at a dissipation rate of 10⁻²
m²/s³, still below 10⁻⁷ J, the lowest energy of Grzegorczyk's graupel impacts, which already broke centimetre flakes
into tens of fragments. These are order-of-magnitude calculations from fall speeds and masses. No experiment has
measured two loose aggregates colliding, so the comparison bounds the collisional route by analogy with denser
impactors.

## What remains: the load analogy and its limits

Without a breakup law, the size at which lunar flakes fail is estimated by comparing the load on their junctions with
that on Earth's largest observed flakes ([survival](turbulence_survival.md)): 2 cm Earth flakes correspond to lunar
flakes of about 5 cm, and the 4–6 cm CASP II aggregates to about 8–15 cm. The estimate holds under four assumptions,
each untested:

1. The weakest junction fails at the same load on both worlds, whatever the flake's size; contact strength measured
   on millimetre ice spheres does not reach decimetre networks.
2. Earth's largest observed flakes were limited by breakup. If time or water limited them instead, the lunar spans are
   conservative.
3. Lunar turbulence lies within a tenth to ten times Earth's. The archived cloud simulations save no dissipation rate.
4. Neither sintering during one to four days aloft nor loosening by repeated sublimation changes the junctions; both
   are documented in the laboratory and neither is quantified for flakes.

A breakup law would need measurements of whole aggregates failing under sustained aerodynamic load or in flake-on-flake
collisions, reported as the largest surviving remnant against load or energy. The strength review's proposal of a
contact-network model that reports the junction strength each geometry needs remains the route to one
([strength literature](strength_literature.md)).

## Search record

The search ran from 14:40 to 15:20 on 25 October with web search, publisher pages and one full PDF. The American
Meteorological Society's pages refused access (HTTP 403), as they did on 24 October, so Lo and Passarelli (1982) and
Gautam and co-authors (2024) are represented by abstracts and search-engine summaries. Kraus's (1968) MIT thesis on
snowflake aggregation was located and not read, and the MIT repository refused Lo's thesis on snow growth for rate
limiting. Sato and co-authors (2006) was read in full from the Montana State University snow-science archive. No papers
were added to the repository.