# The salt for an ocean-like south sea

On 15 November 2026 the author asked to test the viability and the cost and benefit of making the South Pole–Aitken Sea an Earth-like ocean. The Moon's rock gives its water almost no chloride ([the sea chemistry](../research/studies/sea_chemistry/README.md#closed-lakes-soda-lakes-and-salt-seas)), so the sea's salt is an import, and this note asks whether the feed the project already plans to process carries it. The sea's water balance, rivers and mixing are the companion study on `domain/hydrosphere` (`research/studies/south_sea_ocean` there). Numbers marked *screen* come from [sea_salt.py](sea_salt.py) and its product [results/sea_salt.json](results/sea_salt.json) (schema `terluna.resources.sea-salt/1`). The grades and their evidence labels are in [sea_salt_sources.json](sea_salt_sources.json), from a bounded lookup of the same day whose full record, with each value's table or page, is on the research drive (`session-tools/agents/south_sea_salt/`).

## What it finds

**The sea needs 95,650 Gt of chloride and 53,300 Gt of sodium** *(screen)*. Its 4.84 million km³ at seawater's 35.165 g/kg hold 174,000 Gt of salt, 1.3% of the mass of the standing seas and lakes. Half strength takes half of each, and the companion study adds what the crust draws down with the water it takes from the floor, a fifth more within 10,000 years at geography's central case. The five other large seas would take 1,410 Gt of chloride (Moscoviense) to 67,700 Gt (Near-side).

**Lunar ground can give the sodium at great cost and cannot give the chloride.** Highland soil holds 0.34–0.37% sodium (Lunar Sourcebook, Table 7.15), bound in feldspar, so the sea's sodium is the sodium of 1.4–1.6×10¹⁹ kg of soil, three times the sea's own mass of water. Mare basalts and anorthosites hold 3–16 ppm of chlorine (Gargano et al. 2020), about a tenth of what Earth's ocean-ridge basalts hold, and soils are inferred at 7–70 ppm; at that grade the sea's chloride is the chlorine of 2–19% of the whole Moon *(screen)*.

**Rock like the CI chondrites carries seven times more sodium than chlorine, and sea water half as much, so chloride decides and the sodium comes with it.** The compilations put CI chlorine at 717 ± 110 ppm (Lodders and Fegley 2023), with an independent activation study of four falls at 510–710 ppm and the Ryugu sample at 776 ± 21 ppm. One 2017 study by another method on milligram samples found 115 ppm; the compilers reject it for stated analytical reasons, and the lookup read their case and not the 2017 paper, so it is carried here as a disputed low grade. At 717 ppm the sea's chloride is the chlorine of 1.3×10²⁰ kg of such rock (1.2–1.6×10²⁰ kg across 800–600 ppm, 8.3×10²⁰ kg at the disputed grade), which holds 13 times the sea's sodium *(screen)*. That rock is 13–133% of the 10²⁰–10²¹ kg of small bodies finding 1 has processed for the air, and 0.2% of the Kuiper belt's small bodies.

**What rides in the planned imports runs from a sixth of the sea's chloride to eight times it, by the kind of body mined** *(screen; every atom of chlorine in the feed counted as kept)*:

| Import, if it came wholly from | Chlorine carried | Share of the sea's chloride |
|---|---|---|
| The air's nitrogen, from rock like CI chondrites | 7.3×10¹⁷ kg (6.1–10×10¹⁷) | 7.6 times (6.3–10.6); 1.2 times at the disputed grade |
| The air's nitrogen, from ice like comet 67P's | 8.9×10¹⁶ kg, gas chlorine over gas nitrogen; 1.3–2.8×10¹⁶ kg over the gas and dust's nitrogen | 93%; 13–29% |
| The standing water, from rock like CI chondrites | 6.0×10¹⁶ kg (5.0–6.7×10¹⁶) | 62% (52–70%); 10% at the disputed grade |
| The standing water, from ice like comet 67P's | 3.9×10¹⁵ kg (0.6–17×10¹⁵) | 4% (0.6–17%); 12% at the Sun's proportion of chlorine to oxygen |

The two rows of an import are alternatives, and the imports add. Ice like the comet's throughout gives 17% of the sea's chloride on the measured floors and about the whole of it on the gas ratio; rock like CI chondrites throughout gives eight times it. The comet figures are lower bounds: only the gas of one comet has been measured, ammonium chloride led the salts seen in one of its dust grains (Altwegg et al. 2020), and no published budget counts a nucleus's gas, salts and dust together.

**So the chloride for the south sea lies within the scale of processing already planned, and whether it arrives without mining more rests on the chlorine of the bodies chosen, which nobody has assayed.** Sourcing the air's nitrogen from carbonaceous rock would bring the sea's salt with it several times over; sourcing it from comet-like ice may bring most of it or a sixth. In the second case the shortfall is the chlorine of about 10²⁰ kg more of carbonaceous rock.

**The near-Earth asteroids of the source-access study hold salt for aquariums and one small lake** *(screen; upper bounds, all feed taken as carbonaceous rock)*. The vacuum period's 164–291 Gt of feed through year 200 hold the chloride of 4–11 km³ of sea water at the CM to CI grades and of about 1 km³ at the ordinary chondrites' 99 ppm. The wider orbital pool at its median, 1,780 Gt, holds 42–64 km³ of sea water's chloride, enough for the smaller of the seawater lakes sited on 6 November (30–160 km³ each), and the whole near-Earth population 2,300–3,500 km³, 0.05–0.07% of the south sea.

**The salt is one to five years of the delivery stream.** 174,000 Gt is 0.2–1.0% of the 1.7–7.6×10¹⁹ kg of air and water the first screen delivers over 500 years *(screen)*.

## What it leaves open

- **The chlorine of the bodies to be mined.** A bulk figure for a comet nucleus or a Kuiper-belt body does not exist in what was found; the disputed CI grade wants the 2017 paper and any reply read.
- **Keeping the chlorine.** About a third or more of a carbonaceous chondrite's chlorine leaches in water and the rest leaves on heating, and in comet-like feed it comes off with the ammonia as hydrogen chloride (the lookup, heading 6). No study sizing its recovery was found, and none is designed here.
- **The salt's other ions.** Sulphate, magnesium, calcium and potassium are 23,800 Gt together; lunar rock holds all four, and their extraction is not sized.
- **The delivered water's own chloride** ([question 17](questions.md)): on the comet's gas value it is 4% of this sea's need, which is also what would salt the fresh seas to 0.28 g/kg if it were not taken out.

## How it is set up

The need is each sea's volume in geography's atlas times seawater's density and make-up. Each import's chlorine is its mass in the first screen's inventories times a ratio of chlorine to water or to nitrogen: the CI chondrites' group means, with all their hydrogen counted as water, or comet 67P's gas. Shares divide by the south sea's chloride. The asteroid figures read the source-access study's feed and pool masses and multiply by a chondrite group's chlorine. Nothing is modelled, and Earth's values serve only as the seawater recipe.

## The checks

[tests/test_sea_salt.py](tests/test_sea_salt.py): the product binds its inputs and regenerates; seawater's ions stay within its salinity; the need matches the sea chemistry's 174,000 Gt; each import's chlorine is its grade times its mass; the rock for the chloride brings more than ten times the sodium; the asteroid figures scale with their mass and stay under a thousandth of the need; and every value in the sources file names its evidence.
