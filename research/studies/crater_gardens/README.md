# Crater garden settlements: the first catalogue count

This study counts published crater bowls against the [adopted garden-town design](../../../habitation/crater_gardens/README.md): a built ring 500–1,000 m across, with a lake, arcaded terraces and gardens, holding about 5,000 residents in the representative 750 m town. The further 200–500 m comparison examines smaller versions and possible sites for vacuum-era enclosures. The census measures available terrain under stated exclusions; the floor-space calculation measures architectural population.

The subsequent [vacuum-to-open study](../vacuum_crater_buildout/README.md) uses this census for the adopted 2.5% eventual and 0.1% initial buildout fractions across all 200–1,000 m craters. Its initial 500–1,000-person communities occupy one contiguous dome and retain the future-water and climate requirements from the beginning. The earlier abundance scenario below remains the comparison made by this first census.

## Results

The catalogue contains **3,443,077 detections in the adopted 500–1,000 m band**. Of these, **2,177,704 centres fall outside the coarse standing-water masks**, before depth and quality exclusions. The reference depth, source quality, water and regional climate criteria retain **372,128 candidates**. The initial 20,000–100,000-settlement scenario would use about **5–27% of this candidate set** after detailed terrain, spacing, water-supply and engineering assessment.

| Diameter band | Catalogue detections | Candidates under the selected depth family | Selected centres adjacent to water |
|---|---:|---:|---:|
| 200–300 m | 2,550,759 | 249,677 | 44,713 |
| 300–500 m | 2,892,722 | 280,486 | 49,718 |
| 500–1,000 m | 3,443,077 | 372,128 | 70,986 |

The 500–1,000 m row uses the 40–70 m absolute depth window. The smaller rows use the proportional relief comparison described below. They add **530,163 candidates**, with **902,291 selected records** in the combined exported list. About 18–19% of selected centres have a water cell among their eight neighbours; those records deserve early shoreline and flood inspection. The counts retain separate crater detections, and the space for outer gardens and routes still has to be allocated among neighbouring bowls.

| Diameter band and its selected depth family | Rain ≥ 0.5 mm/day | Rain ≥ 1 mm/day | Rain ≥ 2 mm/day |
|---|---:|---:|---:|
| 200–300 m; proportional relief | 340,521 | 249,677 | 141,926 |
| 300–500 m; proportional relief | 380,646 | 280,486 | 161,126 |
| 500–1,000 m; 40–70 m depth | 476,453 | 372,128 | 236,368 |

For 500–1,000 m craters, widening the depth allowance to 30–100 m gives **715,470 candidates** at the reference rain threshold. Keeping every positive measured depth gives **1,074,886** under the same source-quality, water and climate exclusions. These depth choices change the retained count by nearly a factor of three, while the catalogue’s 3.44 million craters in the size band sit within the initial 2–4 million abundance estimate. Every selected region passes the broad mean-ground-temperature window; the regional climate exclusion here comes from rainfall.

Evaluating the adopted floor-space formula at each selected crater’s diameter and summing its residents gives the following architectural population scenarios. These sums assume that all selected records become occupied settlements with the specified coverage, floors and allowance per resident; local siting and construction capacity remain to be assessed. The smaller bands extend the floor-space assumptions below the adopted 500 m lower bound and require their own architectural designs.

| Diameter band | Spacious layout, million residents | Representative layout, million residents | Fuller layout, million residents |
|---|---:|---:|---:|
| 200–300 m | 55.5 | 138.8 | 238.0 |
| 300–500 m | 149.1 | 372.7 | 639.0 |
| 500–1,000 m | 686.7 | 1,716.8 | 2,943.0 |

The full 10,699,651-record count took **64.4 seconds on one CPU worker**, peaking at **123.4 MiB of resident memory**. The corrected-climate drainage and input preparation took **37.4 seconds**, peaking at **1.84 GiB**. [results/first_count.json](results/first_count.json) records every stage, sensitivity, source and hash. The selected records are in `/media/projectspace/terluna-research/crater-settlements/candidates.csv.gz`, with their identifier, present-day coordinates, diameter, depth, source flags and regional values.

## Catalogue and input identity

The input is **LUC-GRAS200 version 1.0.0**, the official April 2026 release of 10,699,651 crater detections at diameters of 200 m and above, published by Zuo, Gao, Liu and colleagues ([article](https://doi.org/10.1029/2025JE009551), [dataset](https://doi.org/10.12350/CLPDS.GRAS.CE2.LUC-GRAS200.vA), [download page](https://moon.bao.ac.cn/Moon/LUC-GRAS200.jsp)). Its CSV supplies location, diameter, rim-to-floor depth, detection confidence and resolution flags. The study streams the complete archive, retains the requested size range and checks the archive's recorded SHA-256, ZIP integrity and published row count.

[inputs.json](inputs.json) records source bytes, hashes, the dataset identifier and the prepared regional context. External bytes and the selected crater list remain on the research drive. The [source register](sources.json) records exactly which text and data were examined: the publisher's abstract and availability statement, the official README and CSV fields. The full article methods have yet to be read. The README describes three additional scores or labels that are absent from the downloaded twelve-column CSV; this calculation uses the fields actually supplied. The source's secondary-crater flag is retained in the selected list because formation as a secondary impact alone supplies no exclusion for habitation.

Wang and colleagues’ 2021 catalogue begins at 1 km. LUC-GRAS200 supplies the smaller craters required by this design. The newer LU6M371TGT catalogue is another available input, recorded in the source register; this count uses one identified catalogue throughout.

## Selection and comparison

| Quantity | Reference selection | Comparisons retained |
|---|---|---|
| Diameter | 500–1,000 m | Separate 200–300 m and 300–500 m bands |
| Published rim-to-floor depth | 40–70 m, a proxy for the 750 m example's relief | 30–100 m; every positive measured depth; relief scaled with diameter |
| Detection confidence | High or medium | High alone; every reported class |
| Resolution flag | 0, the source's unaffected class | All classes |
| Standing water | Centre outside the common-level and rain-fed water masks | Selected centres adjacent to a water cell are reported |
| Regional mean ground temperature | 283.15–308.15 K | Counts before the temperature exclusion |
| Regional mean precipitation | At least 1 mm per 24 hours | At least 0.5 and 2 mm per 24 hours |

The 40–70 m depth window uses the representative town's adopted rim-to-lake relief as a first terrain proxy. The published depth reaches the crater floor, so a lake level and any bed shaping will change the relief seen from its shore. The smaller bands' selected lists use the same proportions as that town: depth/diameter from 40/750 to 70/750. These scaled values define a comparison of similar bowl shapes; the absolute depth windows and all measured depths remain available for considering other forms. The size bands share no records: each lower bound is included, each upper bound is excluded, and the 1,000 m cap is included.

The rain and ground-temperature cuts are analysis choices for a broadly watered, warm garden setting. The rain threshold ranks regional moisture; the lake and garden water balances remain open. The current climate's mean ground temperature is a surface quantity, and human thermal comfort requires the air's temperature, humidity, wind and radiant environment. The model discrepancy in near-ground air remains part of the project's climate uncertainty. Irrigation can extend a settlement beyond the rain cut, so each rain threshold's count is retained separately.

## Water and climate

The regional context takes the corrected design climate `A28R_301k_design_l15_snow_s938_turned_refit_precession`, model years 10–49. Its producer carries the fields back to today's lunar coordinates, where the crater catalogue and terrain are defined. [prepare.py](prepare.py) runs geography's existing drainage calculation with this exact climate product, retaining its output on the research drive, and combines those rain-fed lakes with the atlas's common-level water.

Common-level water comes from the atlas at 4 pixels per degree, about 7.6 km at the equator. Rain-fed lakes use 16 pixels per degree, about 1.9 km. The count samples categorical water membership at each crater centre and interpolates the regional climate bilinearly on its Gaussian grid. A selected centre beside any water cell is flagged for later terrain inspection. Both water resolutions exceed the town's diameter: the first count supplies regional exclusions, while detailed shoreline, flood, groundwater and water-supply assessments need finer terrain and local hydrology. A water-flagged centre can lie in a wholly submerged crater or in a lake whose dry rim could carry the pictured terraces. Resolving that distinction requires the rim profile and water level; the report retains the counts before water exclusion for this reason. The atlas retains the dry Moon's reference geoid pending the planned shoreline revision for the water's gravity and crustal response.

## Population and construction phase

Population is built floor divided by the adopted allowance for each permanent resident. The lake occupies one-quarter of the crater's plan area; buildings cover 30%, 35% or 45% of the rest and average two, three or four storeys. The paired floor allowances are 100, 70 and 70 m² per resident, including homes, work and services. This reproduces the adopted population table. Values below 500 m extend those same architectural assumptions for comparison.

Vacuum-era enclosures control their own internal air and water. Their terrain inventory therefore starts before the planned Open Moon's water and rainfall cuts. For the proportional-relief comparison, the 200–500 m bands contain **1,480,315 candidates after source-quality cuts**, of which **995,649 lie outside the projected waters**. The **530,163** that also pass the future rainfall screen describe a narrower setting for gardens after the open atmosphere arrives. For other bowl forms, the smaller bands contain **4,584,818 records with a positive depth and passing source quality**, before relief, future water and future climate exclusions.

| Smaller-crater band | Proportional relief and source quality | Also outside projected waters | Also passing future garden climate |
|---|---:|---:|---:|
| 200–300 m | 696,058 | 479,983 | 249,677 |
| 300–500 m | 784,257 | 515,666 | 280,486 |
| Total | 1,480,315 | 995,649 | 530,163 |

The future water classification identifies places whose rim and projected water level must be compared before deciding their long-term use. A suitable profile could support a town on the lake shore; other profiles would require water exclusion or a plan for inundation. A feasible enclosure still needs its pressure load, shielding, anchors and failure containment assessed on the chosen terrain.

The supplied discussion of vacuum-era enclosures proposes a 300 m dome and compares pressure, reinforcement and shielding. It is preserved with its source identity on the research drive. Those engineering proposals require their own source and structural assessment. This census supplies crater size, relief and regional setting. A roof’s feasible span and its pressure-vessel design belong to the enclosure assessment; the Open Moon retains its established atmosphere design.

## Reproduction and checks

```sh
python -m research.studies.crater_gardens.fetch_inputs --download
python -m research.studies.crater_gardens.prepare --source-root /path/to/current/main
python -m research.studies.crater_gardens.count
python -m pytest research/studies/crater_gardens
```

The default work directory is `/media/projectspace/terluna-research/crater-settlements`; each command accepts `--work-dir` for another persistent location. The prepared context is pinned by its byte count and SHA-256. Re-preparing it creates a new input snapshot whose identity must be reviewed and recorded in inputs.json before it is accepted by the counter.

The tests use independent annular-area sums, the diameter-squared population scaling, coordinate seam and pole controls, an analytically interpolated climate field, unsorted synthetic source records, exclusion boundaries, missing measurements and a changed-input failure. The aggregate product binds the counter, manifest, design and context to their exact bytes. Software checks, source reading and site inspection remain distinct records of evidence.

[checks.json](checks.json) records the verification: all seven study tests passed, and the final `make check` passed with 5,219 Python tests, 31 Python skips and 191 JavaScript tests, together with the layer, calendar, provenance and ensemble-workspace checks. The full check took 15.3 minutes (0.254 machine-hours), peaking at 2.54 GiB of resident memory on one numerical thread. The original article's full methods and the local terrain and enclosure assessments remain the further evidence work described above.
