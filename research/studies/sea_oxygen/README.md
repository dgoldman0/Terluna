# The held fresh seas' oxygen and mixing in their basins at 301 K

18 October 2026. The author set the resort's ground at 301 K on 18 October (the register on domain/climate) and asked
for this study after [the deep lakes' anoxia](../lake_anoxia/README.md), to run the seas' own basins, which the lake
studies left. On the lake model as corrected that day, [the waters' oxygen](../lake_oxygen/README.md) finds the wind's
energy below the night's mixing turning its 352 m crater lakes over every second to fourth lunar month by Osborn's rule
and about once a year by Shih and Bouffard's energetic regime, their deep water keeping its oxygen in nutrient-poor
water and, in richer water, resting on the mixing's form and the warmth. The seas are the common-level waters, 28% of
the surface, held fresh at 0.14 g/kg ([the seas' make-up](../sea_chemistry/README.md)); they are far larger and deeper
than the lakes and span latitudes, on a Moon whose equator lies 1.5° from the ecliptic, so that the lunar night of about
354 hours is their only regular cooling. [seas.py](seas.py) labels the seas' basins on LOLA and writes
[results/seas.json](results/seas.json) (schema `terluna.research.sea-oxygen-seas/1`); [extract.py](extract.py) reads the
GCM's weather at their water cells; [run.py](run.py), with [setting.py](setting.py), [columns.py](columns.py) and
[analysis.py](analysis.py), writes [results/sea_oxygen.json](results/sea_oxygen.json) (schema
`terluna.research.sea-oxygen/1`); [sources.json](sources.json) records the sources beyond those it shares with the
waters' oxygen study.

**The model.** The columns run on the lake model as the independent review of lake.py, lake_oxygen and lake_anoxia on
18 October left it (biosphere/ecology/lake.py, sha256 9a34ba69…), with the surface exchange paired at the GCM's lowest
level, the anaerobic products in the density, the night's respiration running down with a labile pool, and the deep
mixing from the wind, the author's central treatment, with its bracket. The full set of columns turns the dissipation
below the mixed layer into a diffusivity by Osborn's rule, as lake.py does (the two Lake-number rows aside). Two further
columns, at Procellarum's and Ingenii's typical cells, take Shih et al.'s (2005) form for energetic turbulence instead,
as Bouffard & Boegman (2013) carry it to lakes, which mixes the deep water less; the deep-oxygen findings below give the
range between the two forms where both ran. The findings are runs of one lake model forced by one GCM at T21, checked by
their settle tests and budgets; no sea like these has been measured. Two columns of Ingenii run on the model before the
review are kept as the before case, in `research/runs/sea_oxygen/before_cases.json`.

## What it finds

**The seas are deep, their shelves narrow, and a quarter of them lie in crater basins.** At the common level of
−1,654.2 m they cover 10.62 million km², 28.0% of the surface on LOLA's 16 pixel/degree grid, and hold 10.7 million km³
at a mean depth of 1,009 m. Half their area lies over water deeper than 719 m; 0.75% of it is shallower than 10 m, 2.2%
than 30 m and 7.2% than 100 m. They form 5,708 bodies. The four seas over 100,000 km² hold 76% of the area and 83% of
the volume; the 5,704 crater seas under 100,000 km² hold the rest, their deepest points at an area-weighted median of
1,709 m and their median latitude 36.5°. A third of the seas' area, and nearly half their volume, lies poleward of 45°:

| | 0–15° | 15–30° | 30–45° | 45–60° | 60–75° | 75–90° |
|---|---|---|---|---|---|---|
| Share of the seas' area | 18.7% | 23.1% | 22.8% | 22.0% | 10.2% | 3.3% |
| Share of their volume | 11.3% | 18.0% | 23.7% | 27.5% | 15.6% | 3.8% |

**Each large sea is several deep basins sharing only their top hundred metres or so.** The Near-side Sea, 5.19 million
km² from 32° S to 71° N, divides below 119–181 m into four basins: the northern one of Oceanus Procellarum, Mare
Imbrium, Serenitatis and Frigoris; the south-western one of Mare Nubium and Humorum; Mare Crisium, the deepest at
3,891 m; and Mare Fecunditatis. Inside the northern basin Nectaris joins at 411 m and Imbrium meets Serenitatis at
759 m, and Nubium meets Humorum at 391 m. The South Pole–Aitken Sea, 2.31 million km² from 22° S to 87.5° S, holds the
deepest water, 6,795 m in Ingenii's basin at 70° S, which Schrödinger's polar basin joins at 76 m and Drygalski's at
164 m. Humboldtianum (48–66° N) and Smythii (on the equator) each centre on one deep basin, 4,388 and 4,087 m deep, with
smaller ones beside it, Marginis's beside Smythii's. Below its sill a basin's deep water is renewed only by its own
surface or by denser water spilling over the sill. The basins run:

| Basin | Area | Mean depth | Deepest | Latitudes | Sills to the others of its sea | Stands for |
|---|---|---|---|---|---|---|
| Procellarum, Imbrium, Serenitatis, Frigoris | 3,146,800 km² | 726 m | 3,237 m | 23° S–71° N | 119–128 m | 29.8% |
| Ingenii, South Pole–Aitken | 2,068,300 km² | 2,180 m | 6,795 m | 88–23° S | 76 and 164 m | 19.7% |
| Nubium, with Humorum and Cognitum | 891,200 km² | 478 m | 3,217 m | 32° S–7° N | 119–236 m | 8.7% |
| Humboldtianum | 253,700 km² | 1,175 m | 4,388 m | 48–66° N | | 2.6% |
| Crisium | 242,200 km² | 1,641 m | 3,891 m | 8–28° N | 119–181 m | 2.4% |
| Fecunditatis | 228,400 km² | 504 m | 2,501 m | 12° S–10° N | 119–181 m | 2.2% |
| Smythii | 168,500 km² | 1,734 m | 4,087 m | 9° S–10° N | | 2.5% |
| Schrödinger, South Pole–Aitken | 106,600 km² | 1,569 m | 3,610 m | 87–70° S | 76 m | 1.9% |
| Western Procellarum's shallows, 7° N 302° E | 22,100 km² | 177 m | 475 m | 4–10° N | 119–236 m | 5.8% |
| Drygalski, South Pole–Aitken | 16,900 km² | 1,922 m | 3,611 m | 83–77° S | 76 and 164 m | 0.2% |
| Rozhdestvenskiy, by the north pole | 14,100 km² | 615 m | 2,819 m | 83–87° N | | 0.1% |
| Jules Verne, the crater seas | 10,200 km² | 578 m | 2,120 m | 37–33° S | | 24.1% |

The rule takes every basin of at least 100,000 km², which brings in the largest seas and the deepest point; the basin
nearest each pole among those of at least 10,000 km²; one crater sea, the body of 10,000–100,000 km² within 10° of the
crater seas' median latitude whose deepest point lies nearest their median; and one shallow sea, the largest region of a
large sea whose deepest point lies above 500 m. Each stands for its own area, the other regions of its sea nearest it in
mean depth and, for Jules Verne, every crater sea not run, so that together they stand for all the seas. The large seas'
other regions, 0.90 million km², are mostly shallow: 0.42 million km² of them lie in regions whose deepest point is
above 500 m, at a mean depth of 131 m.

**At 301 K the held seas hold 7.66 mg/L at saturation, 1.5% less than at the resort's 300 K.** Under the Moon's air,
17.46% oxygen at 1.2 atm, the seas' median surface at 301 K, 29.6 °C on climate_targets' basis (0.520 of the way from
the lock to the hot world), holds 7.658 mg/L in the held seas' 0.14 g/kg, 21% more than seawater at the same warmth
would hold. Beside lake_oxygen's three targets:

| | 301 K | The resort, 300 K | The hot world | The lock |
|---|---|---|---|---|
| Seas' surface, median | 29.6 °C | 28.7 °C | 34.6 °C | 24.2 °C |
| Oxygen at saturation | 7.66 mg/L | 7.78 mg/L | 7.06 mg/L | 8.43 mg/L |
| Band medians, 0–15° to 60–90° | 7.49–8.00 mg/L | 7.59–8.16 mg/L | 6.96–7.21 mg/L | 8.15–9.05 mg/L |
| At the seas' 99th-percentile warmth | 7.31 mg/L (32.4 °C) | 7.41 mg/L (31.6 °C) | 6.82 mg/L (36.8 °C) | 7.92 mg/L (27.7 °C) |
| Every water cell's 3-day means at its own pressure: median; 1st–99th percentile | 7.58; 7.05–8.18 mg/L | 7.69; 7.15–8.34 mg/L | 6.96; 6.53–7.25 mg/L | 8.25; 7.62–9.19 mg/L |
| Earth's seawater (35 g/kg) at the seas' median warmth | 6.32 mg/L | 6.41 mg/L | 5.86 mg/L | 6.90 mg/L |

The bands' medians run from 31.0 °C within 15° of the equator to 27.1 °C poleward of 60°. The sea cells alone, at their
mean 1.197 atm, have a median of 7.63 mg/L across their 3-day means.

**Each large basin has a cold end, 0.5–2.9 K colder in its coldest spell than its typical water.** In the GCM's own
water at 301 K (its 50 m slab's surface, shifted by the share of the way as the forcing is), the run basins' cells
average 25.4–32.5 °C, from Rozhdestvenskiy by the north pole to Smythii on the equator. A cell's coldest 3-day mean lies
only 0.4–1.1 K below its mean, since the lunar night cools the slab little. Each basin's typical and coldest cells:

| Basin | Typical cell: place; mean, coldest 3-day mean | Coldest cell: place; mean, coldest 3-day mean | Basin's area within 0.5 K of its coldest spell |
|---|---|---|---|
| Procellarum | 36° N 6° E; 30.1, 29.4 °C | 69° N 338° E, Frigoris's northern edge; 27.1, 26.5 °C | 0.2% |
| Ingenii | 53° S 214° E; 28.6, 27.9 °C | 86° S 225° E; 26.4, 25.9 °C | 10% |
| Nubium | 14° S 338° E; 31.2, 30.1 °C | 31° S 349° E; 29.8, 28.8 °C | 12% |
| Humboldtianum | 58° N 73° E; 28.3, 27.7 °C | 64° N 107° E; 27.2, 26.6 °C | 18% |
| Crisium | 19° N 62° E; 31.5, 30.7 °C | 25° N 51° E; 30.9, 30.0 °C | 40% |
| Fecunditatis | 3° N 62° E; 31.5, 30.6 °C | 8° S 51° E; 30.9, 29.8 °C | 35% |
| Smythii | 3° S 84° E; 32.0, 30.9 °C | 8° S 84° E; 31.5, 30.4 °C | 33% |
| Schrödinger | 80° S 124° E; 26.5, 26.1 °C | 75° S 118° E; 26.1, 25.6 °C | 48% |

Western Procellarum's shallows (3° N; 31.3, 30.3 °C), Drygalski (80° S; 26.0, 25.6 °C), Rozhdestvenskiy (86° N; 25.5,
25.1 °C) and Jules Verne (36° S; 29.4, 28.4 °C) run one column each, their cells within 0.5 K of one another. A basin's
densest water forms at its coldest end; in Procellarum and Ingenii, the two basins spanning the most latitude, that
end's coldest spells run 2.9 and 2.0 K below the coldest their typical water reaches.

**Under Osborn's rule a column's deep water turns over within years, and its oxygen follows the trophic state.** In the
saved decade every typical column mixes to its floor 2–32 times, 20–330 times a century, as often as in the coarse
passes' last half to within 3% except at Rozhdestvenskiy, 9% less often, and Crisium and Smythii, a quarter to a third
less often at the analysis step. Counted by the lunar cycles they touch, as the waters' oxygen study counts, the
decade's overturns fill 2–50 of its 120 cycles, up to 1.75 times as many, since an overturn often spans the turn of a
cycle. The deepest tenth is 0.3–3.8 years from its last contact with the surface. The deep water lies within 0.02 K of
the coldest the cell's surface reaches (0.07 K in the western shallows), so that a cold night's convection runs to the
floor, in Procellarum, Nubium, Humboldtianum and the western shallows on more than one night in ten; below the mixed
layer the diffusivity the wind's dissipation gives over so weak a stratification reaches the cap of 0.1 m²/s across
33–69% of the interfaces. The deepest tenth's oxygen in the saved decade (mg/L), with the overturns a century in the
saved decade and in the coarse passes' last half and the deepest tenth's age, at the wind's central share and, in
brackets, at the ends of its 0.2–0.4% bracket (whose columns pass the record 30 times from the central one, so that
their coarse last half is 150 years):

| Basin | Stands for | Overturns a century: saved decade; coarse passes | Deepest tenth's age | Oligotrophic | Mesotrophic | Eutrophic | New reservoir |
|---|---|---|---|---|---|---|---|
| Procellarum | 29.8% | 143 (102–163); 140 (101–160) | 0.6 (0.4–0.7) y | 5.57 (4.78–5.92) | 4.05 (2.93–4.58) | 2.74 (1.45–3.35) | 0.87 (0–1.73) |
| Ingenii | 19.7% | 41 (41–41); 40 (40–40) | 3.8 (2.7–4.5) y | 2.49 (1.54–3.13) | 0.01 (0–0.73) | 0 (0–0) | 0 (0–0) |
| Nubium | 8.7% | 184 (143–204); 180 (150–200) | 0.4 (0.4–0.6) y | 5.93 (5.50–6.17) | 4.64 (4.01–5.00) | 3.54 (2.82–3.95) | 1.96 (0.95–2.53) |
| Humboldtianum | 2.6% | 163 (123–174); 160 (120–180) | 0.8 (0.7–1.0) y | 5.87 (5.44–6.09) | 4.39 (3.82–4.71) | 2.98 (2.31–3.36) | 1.18 (0.27–1.70) |
| Crisium | 2.4% | 51 (31–72); 70 (70–70) | 1.7 (1.7–3.8) y | 2.45 (1.75–3.17) | 0 (0–0.38) | 0 (0–0) | 0 (0–0) |
| Fecunditatis | 2.2% | 143 (143–184); 140 (140–170) | 0.4 (0.4–1.0) y | 5.22 (4.37–5.59) | 3.32 (2.00–3.91) | 1.71 (0.23–2.43) | <0.01 (0–0.39) |
| Smythii | 2.5% | 41 (31–41); 59 (44–60) | 1.7 (1.7–2.6) y | 0.55 (0–1.22) | 0 (0–0) | 0 (0–0) | 0 (0–0) |
| Schrödinger | 1.9% | 72 (61–72); 70 (65–70) | 1.9 (1.9–2.5) y | 1.63 (0.68–2.35) | 0 (0–0) | 0 (0–0) | 0 (0–0) |
| Western shallows | 5.8% | 327 (225–511); 330 (230–500) | 0.3 (0.2–0.4) y | 6.30 (5.79–6.53) | 4.96 (4.15–5.36) | 3.70 (2.73–4.19) | 2.10 (0.76–2.80) |
| Drygalski | 0.2% | 20 (10–20); 20 (11–20) | 3.7 (2.4–4.8) y | 2.47 (1.18–3.32) | 0 (0–0.84) | 0 (0–0) | 0 (0–0) |
| Rozhdestvenskiy | 0.1% | 82 (72–112); 90 (70–110) | 1.5 (1.1–2.5) y | 4.37 (3.35–4.90) | 1.73 (0.28–2.55) | <0.01 (0–0.40) | 0 (0–0) |
| Jules Verne, the crater seas | 24.1% | 92 (72–143); 90 (80–150) | 0.7 (0.7–1.6) y | 3.04 (0.96–3.65) | <0.01 (0–0.41) | 0 (0–0) | 0 (0–0) |

Oligotrophic water keeps 1.6–6.3 mg/L in the deepest tenth of every basin but Smythii (0.55). Mesotrophic water keeps
1.7–5.0 mg/L in Procellarum, Nubium, Fecunditatis, the western shallows, Humboldtianum and Rozhdestvenskiy and next to
none in Ingenii, Crisium, Smythii, Schrödinger, Drygalski and the crater seas; eutrophic water keeps 1.7–3.7 mg/L only
in the Near-side Sea's shallower basins and Humboldtianum. By the share of the seas each basin stands for, the deepest
tenth averages 4.1 mg/L in oligotrophic water, with 4% of the seas' area under 2 mg/L; 2.1 mg/L in mesotrophic water,
with half the area under 0.5; and 1.5 mg/L in eutrophic water. At the wind's low share these are 3.1, 1.6 and 0.9 mg/L,
with half the seas' area under 2 mg/L even in oligotrophic water, and at its high share 4.6, 2.6 and 1.7. Two of the
bracket's columns had not settled after their ten passes at the analysis step. Crisium's at 0.4% was still losing
oxygen, its mesotrophic deepest tenth falling from 1.28 to 0.39 mg/L over the passes, so its 0.38 is an upper value;
Fecunditatis's at 0.2% was still gaining, from 1.61 to 1.97 mg/L over its last eight passes while its bottom water's age
fell from 1.5 to 1.2 years, so its 2.00 is a lower value. The coarse passes, in 4-step means, leave more oxygen in the
deep water than the analysis step: their last half gives the mesotrophic deepest tenth 0.67 mg/L in Ingenii's typical
column against the saved decade's 0.01, 0.34 in Crisium's against 0, 0.47 in Drygalski's against 0 and 1.35 in Ingenii's
coldest against 0.58. The findings read the analysis step.

**Under the other form of the deep mixing the two largest basins' deep water turns over rarely and holds little or no
oxygen.** With Shih et al.'s form for energetic turbulence in place of Osborn's rule, the diffusivity reaches the cap at
0.1–0.4% of the interfaces below the mixed layer against 63–64% under Osborn's, and the nights' mixing stops short of
the floor (its 99th percentile at 168 m in Procellarum and 445 m in Ingenii). Ingenii's column settles at the analysis
step: it turns over once in its saved decade, touching two lunar cycles, and 21 times a century in the coarse passes,
against 4 times in the decade and 40 a century under Osborn's rule, and its deepest tenth, 4.9 years old against 3.8, is
anoxic in every trophic state. Procellarum's does not settle. In its ten passes at the analysis step, 100 years, it
turns over only in the last: through the 90 years before, no overturn reached its floor, its bottom water aged to
39 years and its mesotrophic deepest tenth was anoxic throughout. Its saved decade is the decade after that overturn,
and turns over again 3.5 years in; it holds the column just renewed. Its oligotrophic deepest tenth, falling about
0.17 mg/L a year between the overturns (from 1.11 to 0.77 mg/L in the decade's first three years, from 2.26 to 1.56 in
its last five), averages 1.54 mg/L, an upper value, and its deepest tenth's 4.7 years are a lower value; its 20
overturns a century, in the saved decade as in the coarse passes, stand against one pass of the ten at the analysis step
with an overturn. A longer run at the analysis step would settle it. Between the two forms the deepest tenth's oxygen in
the saved decade runs (mg/L):

| Basin | Oligotrophic | Mesotrophic | Eutrophic | New reservoir |
|---|---|---|---|---|
| Procellarum | 1.54 or less to 5.57 | 0.01–4.05 | 0–2.74 | 0–0.87 |
| Ingenii | 0–2.49 | 0–0.01 | 0 | 0 |

With the second form at these two basins, half the seas' area, and Osborn's rule elsewhere, the seas' deepest tenth
averages at most 2.5 mg/L in oligotrophic water, 54% of the area under 2 mg/L; 0.9 mg/L in mesotrophic water, 81% of the
area under 0.5; and 0.6 mg/L in eutrophic water, against Osborn's 4.1, 2.1 and 1.5. Procellarum's renewed decade makes
these means, and its own oligotrophic 1.54 mg/L, upper values; the shares under 2 and 0.5 mg/L stand. The other basins
ran under Osborn's rule alone.

**At their cold ends the basins make colder deep water and keep less oxygen in it.** Each basin's coldest column makes
deep water 0.6–2.7 K colder than its typical column, colder than the typical column's surface ever gets, and turns over
20–153 times a century in its saved decade (20–160 in the coarse passes). Its deepest tenth in the saved decade, under
Osborn's rule (mg/L):

| Basin | Overturns a century: saved decade; coarse passes | Deepest tenth's age | Deep water | Oligotrophic | Mesotrophic | Eutrophic | New reservoir |
|---|---|---|---|---|---|---|---|
| Procellarum | 102; 100 | 1.2 y | 25.74 °C (typical 28.43) | 2.71 | 0 | 0 | 0 |
| Ingenii | 20; 20 | 4.2 y | 26.31 °C (typical 27.02) | 3.04 | 0.58 | 0 | 0 |
| Nubium | 153; 160 | 0.5 y | 27.04 °C (typical 28.70) | 5.25 | 3.24 | 1.55 | <0.01 |
| Humboldtianum | 61; 60 | 1.4 y | 26.09 °C (typical 26.92) | 3.32 | 0.26 | 0 | 0 |
| Crisium | 61; 68 | 1.6 y | 27.70 °C (typical 28.93) | 1.99 | 0 | 0 | 0 |
| Fecunditatis | 133; 130 | 0.9 y | 27.92 °C (typical 28.54) | 5.40 | 3.67 | 2.25 | 0.19 |
| Smythii | 41; 60 | 1.5 y | 27.35 °C (typical 28.14) | 0.36 | 0 | 0 | 0 |
| Schrödinger | 82; 80 | 1.4 y | 25.34 °C (typical 25.90) | 2.63 | 0 | 0 | 0 |

With each basin's coldest column standing for it where it has one (70% of the seas' area), the seas' deepest tenth
averages 3.3 mg/L in oligotrophic water, 5% of the area under 2 mg/L, and 0.8 mg/L in mesotrophic water, 83% of the area
under 2 mg/L and 64% under 0.5. Ingenii's coldest column alone among the central columns does not settle at the analysis
step: its mesotrophic deepest tenth was still falling over the ten passes, from 0.95 to 0.60 mg/L, 3–4% a pass at the
end, so its 0.58 mg/L and the coldest columns' seas-wide 0.8 mg/L are upper values.

**An overturn renews the deep water's age more than its oxygen.** Procellarum's coldest column turns over ten times in
its saved decade and its deepest tenth is 1.2 years old, yet its mesotrophic deep water is anoxic. The whole column
mixes for a night or two at a time, the column's mean oxygen, mostly that of its anoxic depths, comes to the surface,
and the air restores little before the column stratifies again: where a basin's deep water is anoxic, its mesotrophic
column's surface after an overturn falls to 0.001–0.8 mg/L. A column turns its whole basin over at once; a sea's
overturn is local, at its cold end and on its slopes.

**Without the wind's mixing the deep water is cut off.** With Yeates & Imberger's diffusivity from the Lake number, or
the first round's fixed dissipation of 10⁻¹⁰ W/kg under a cap of 10⁻⁴ m²/s, in its place, Procellarum's and Ingenii's
typical columns stop turning over: none of the four turns over in its 100 years at the analysis step or in its saved
decade (in the coarse passes Procellarum's under the fixed dissipation overturned 8 times a century, the others not at
all). Their bottom water ages to 99–394 years by the end of those passes, their deepest tenth is 79–404 years old in the
saved decade and still ageing, and it is anoxic in every trophic state. Leaving Ingenii's respiration products out of
the density, or giving Procellarum oligotrophic water's clarity, moves the deep oxygen by up to 1.5 mg/L.

**A basin's renewal rests on its cold end, which in Procellarum is small.** No typical column's surface reaches its
basin's coldest deep water, so with exchange the deep water takes the cold end's temperature and only the cold end
renews it. The share of a basin whose coldest 3-day means lie within 0.5, 1 and 2 K of its coldest cell's is 0.2%, 12%
and 31% in Procellarum, 10%, 25% and 50% in Ingenii and 12%, 29% and 98% in Nubium. Spread over those shares, the
coldest column's renewal, ten overturns in Procellarum's saved decade, once in 0.98 years, becomes once in
445, 7.8 and 3.1 years (455, 8.0 and 3.2 at the coarse passes' rate), and Ingenii's, once in 4.9 years, becomes once in
48, 20 and 9.8. Taking the deepest tenth's reference demand over the added interval, oligotrophic deep water keeps
0–2.3 mg/L in Procellarum and 0.3–2.7 in Ingenii, mesotrophic water is anoxic in Procellarum and keeps at most 0.21 mg/L
in Ingenii, and Nubium, its cold end broad, keeps 1.8–3.2 mg/L in mesotrophic water. The reference demand is the central
columns', fixed by supply at any warmth; the warm kinetics would raise it 5.3–6.5 times at the coldest columns' deep
water (25.3–27.9 °C), and with it the screen's losses.

**Procellarum's cold north can fill the Near-side Sea's other basins.** Its coldest column's water at the sills,
25.7 °C at 119–128 m, is 2.2–3.2 K colder than the deep water of Nubium, Crisium, Fecunditatis and the western
shallows, and so could spill over and sink into them; the cold ends of Nubium, Crisium and Fecunditatis make water up to
1.9 K colder than their neighbours' typical deep water as well, so the Near-side Sea's basins share their densest water
over their sills, Procellarum's sinking beneath the rest. In the South Pole–Aitken Sea the dense water runs from the
small polar basins into the large one: Schrödinger's and Drygalski's water at their 76 and 164 m sills is 1.4–1.7 K
colder than Ingenii's typical deep water, while Ingenii's coldest water there, 26.3 °C, is 0.4–0.7 K warmer than theirs.

**The shelves and bays keep their nights oxygenated; the bays' floors rest on the wind's mixing.** Under the wind's
share lake_oxygen's 3 m productive pan falls in the median night to 3.7–6.2 mg/L at the surface in eutrophic water at
every cell, and the 10 m pan to 4.0–6.7, with no hypoxic night at the surface and the eutrophic floors above 2 mg/L; the
lakes' own 3 m pans at 301 K fall to 3.3–3.9. Under the first round's fixed dissipation, the bracket for the wind's
share concentrated in thin water, the surfaces barely change (3.7–5.4 mg/L) while the floors go hypoxic: 13–44% of the
time in the 3 m pans and 10–17% in the 10 m pans in eutrophic water, the 10 m floors falling to 0.5–0.7 mg/L in the
median night and staying under 2 mg/L for 28–54 hours of it (148–150 hours in new reservoirs).

**The wind's stress outweighs the stratification, and the Moon's rotation shapes what it does.** By the lunar month's
thermocline the Lake number has medians of 0.003–0.14 and lies under 1 on 95–100% of days in every central column; by
the day's warm layer its medians are 0.0015–0.05. From each typical column's own stratification taken as two layers, the
day's mixed layer, a median 3–6 m, over the deepest tenth, 0.6–3.9 K colder, the first internal mode runs at a median
0.03–0.085 m/s, and at the basins' median latitudes the internal Rossby radius is 6–170 km, 8–200 km at the
90th-percentile speed; with the nights' 20 m mixing and a 3 K step, the review's estimate, it is 40–130 km at 60–15° and
about 240 km at the equator. The large basins and the polar ones are 10–90 radii across (Procellarum 74, Ingenii 90,
Nubium 15, Humboldtianum 39, Crisium 10, Schrödinger 41), so the wind drives Ekman transport across the wind, upwelling
along the shores within a radius of them, and Kelvin waves along the shores and Poincaré waves across the basins, and no
tilt grows across a whole basin with its length. The equatorial basins, Fecunditatis and Smythii 2.5–3 radii across and
the western shallows about one, and the crater seas, Jules Verne 4, come nearer a lake's tilt. A column holds none of
these.

**The model before the review cut the deep water off.** On lake.py before the review, with its fixed dissipation and
cap and the surface exchange paired at 10 m, Ingenii's typical and coldest columns last turned over in the 51st and
40th of 120 passes, and their deepest tenth aged to 306 and 336 years and lay anoxic in every trophic state, 1.1 and
1.5 K colder than under the corrected model's wind.

**What the seas need.** Their deep water's oxygen rests first on how fast the water below the mixed layer mixes, which a
column can only bracket: under Osborn's rule every basin's deep water is renewed within years and oligotrophic seas keep
oxygen in their depths almost everywhere; under the published form for energetic turbulence the two largest basins' deep
water turns over ten to twenty times a century, Procellarum's once in its century at the analysis step, and holds little
or no oxygen; under the first round's fixed dissipation or the Lake number's diffusivity it is cut off for centuries.
The first need is the seas' own circulation: a three-dimensional model of the basins on their LOLA bathymetry under the
GCM's forcing, with the cold ends' dense water, the sills, the slopes and the Moon's rotation, decides where between
these the deep water lies. The second is low nutrient loading. In every column run for them, at either end of the wind's
bracket, under either form and at their coldest cells, mesotrophic water leaves the deepest tenth of Ingenii, Crisium,
Smythii, Schrödinger, Drygalski and the crater seas, together half the seas' area, under 1 mg/L; only oligotrophic water
keeps the depths above 2 mg/L over most of the seas' area under Osborn's rule (96% of it), a share that falls to 46%
with the second form at the two largest basins, as only oligotrophic lakes fit the air's held methane ([the deep lakes'
anoxia](../lake_anoxia/README.md)). The third is the cold ends' cooling: the Near-side Sea's densest water forms at
Frigoris's northern edge, the 0.2% of Procellarum's area within 0.5 K of its coldest, and the South Pole–Aitken Sea's in
Schrödinger and Drygalski, so heat put into the water there would weaken the renewal each sea draws on. The shelves and
bays need least: their surfaces' median night stays at 3.7 mg/L or more in eutrophic water in every case run, and only
their floors, under the weaker mixing, go hypoxic, the 10 m bays' floors a tenth to a sixth of the time in eutrophic
water. What the anoxic depths make, the lake methane study's tracer takes from the saved columns.

## How it is set up

**The seas.** [seas.py](seas.py) takes every cell of LOLA's 16 pixel/degree grid whose height above the GRAIL geoid lies
below the atlas's common level, −1,654.2 m, as sea, as geography's drainage takes its outlets, and labels the cells'
connected bodies with geography's own labeller, joined across the 0/360 seam and the poles. The four seas over
100,000 km² are divided into basins at 250 m, the depth at which the Near-side Sea falls into its four main parts; each
basin takes the shallower water nearest it through the sea, growing cell by cell. A sill is the depth of the deepest
water path between two deepest points, found by bisection to 2 m. Each basin carries its hypsography (the share of its
surface area deeper than each fraction of its deepest point, as lake_oxygen's lakes.py gives it), the shares of its area
shallower than 2–200 m, and the T21 cells that hold its water.

**The forcing.** [extract.py](extract.py) reads the GCM's 3-day means, with lake_oxygen's readers, at every T21 cell
holding a run basin's water that PlaSim takes as water at the seas' level: 486 cells, in the design case (the lock) and
the hot world's test, years 10–19, and in the design case's two years written every 6.4 hours, which set the wind's eddy
factor; and at every water cell, the surface temperature and pressure lake_oxygen keeps for its waters. The resort at
301 K is the lake kinds study's [Setting301](../lake_kinds/setting.py): the lock's weather moved by 0.520 of the way to
the hot world's test, the share climate_targets' interpolation gives 301 K on the ground's temperature between its
design case (295.02 K) and its test (306.52 K). [setting.py](setting.py) keeps everything else of it (its lapse rate,
0.84 K/km, carries air down to a lake's level and has nothing to carry at a sea cell) and reads this study's cells, from
which it recomputes the Sun's phase and the eddy factor as lake_oxygen does; lake_oxygen's windows(), resort_windows()
and hourly() run on them unchanged. A sea cell sits at the seas' level, so its air, sky and pressure are the cell's own,
as lake_oxygen reads its water cells. Each basin takes its typical cell, whose mean surface temperature at 301 K lies
nearest the area-weighted median of its cells, and its coldest, whose coldest 3-day mean is the lowest among the cells
holding at least 1,000 km² of it; where the two coldest 3-day means lie within 0.5 K the basin runs one column.

**The columns.** [columns.py](columns.py) runs the lake model,
[biosphere/ecology/lake.py](../../../biosphere/ecology/lake.py), on each basin's own hypsography to its deepest point,
with lake_oxygen's layers (0.5 m to 20 m, then growing to 25 m) and its physics() and scenarios(), so that the seas'
settings are the lakes' own: the clarity and solutes of mesotrophic water, as lake_oxygen's deep lakes take them, with
aerobic respiration's CO2 and the debt's anaerobic CO2 and methane in the density; below the mixed layer the wind's
share of its energy flux at 10 m, 0.3% central (Wüest, Piepke & Van Senden 2000), spread over the water below, under a
cap of 0.1 m²/s; the surface exchange paired at the GCM's lowest level, as PlaSim pairs it; and 24 oxygen columns, four
trophic states by two transfer velocities by two kinetics with the euphotic water's labile pool decaying over 4 days,
and the pool's 2- and 7-day brackets. The held seas' 0.14 g/kg sets the density's salinity term, and the steps are
7,150 s. lake_oxygen settles its deep lakes' coarse passes by its settle test and then settles them again at the
analysis step; this study keeps a fixed number of recorded coarse passes and then runs the same settle test at the
analysis step. A central column starts warm and saturated, 1 K above the air's mean, as lake_oxygen's deep lakes start,
and passes over the ten-year record in 4-step means 80 times, 800 years of the record repeated, each pass recorded and
the column saved after each so a stopped run resumes; its long state is read from the last 40. It then passes over the
record at the analysis step until lake_oxygen's settle test passes there, the deep water's temperature, the coupled
column's deepest-tenth oxygen, the bottom water's age and its solute each steady from one pass to the next, or until
lake_oxygen's ten passes; a column cut off from the surface keeps ageing and does not pass. The record is then run once
more and saved with the column every day in the lake kinds study's arrays, with the labile pool, for the lake methane
study's tracer; the decade saved is the one the findings read. Each typical column also runs at the wind's share's
bracket, 0.2% and 0.4%, and the coldest columns at the central 0.3%, as the author set the brackets. At Procellarum and
Ingenii, Yeates & Imberger's diffusivity from the Lake number and the first round's fixed dissipation (10⁻¹⁰ W/kg under
a cap of 10⁻⁴ m²/s) run as sensitivity, with Ingenii's respiration products left out of the density and Procellarum's
physics that of oligotrophic water. These variants start from their basin's central typical column and pass the record
30 times before the same settle test. The shelves and bays are lake_oxygen's 3 m productive pan, with eutrophic physics,
at every column's cell, and a 10 m pan at each typical cell; each spins up over the record's first year in hour steps,
as lake_oxygen's shallows do, before its decade is saved.

**The two forms of the deep mixing.** lake.py turns the dissipation ε below the mixed layer into a diffusivity by
Osborn's rule, K = 0.2 ε/N², at every buoyancy Reynolds number Re_b = ε/(νN²). In the seas' nearly neutral deep water
Re_b runs to 10³–10⁵, where Shih et al. (2005) found the mixing's efficiency falling, K = 2ν Re_b^½ above Re_b = 100,
Osborn's rule from 7 to 100 and molecular diffusion alone below 7, the regimes Bouffard & Boegman (2013) carry to lakes.
Above Re_b = 100 the second form gives 0.1 Re_b^½ times less: 3 times at 10³ and 32 at 10⁵. lake_oxygen's
[mixing_forms.py](../lake_oxygen/mixing_forms.py) swaps lake.py's diffusivity() for those regimes inside the process
that runs a column (ν = 10⁻⁶ m²/s, the cap and the unstable rule unchanged), and Procellarum's and Ingenii's typical
columns run under it from the warm start, as the central columns run, with the same passes at the analysis step; the
record of those passes, whether each overturned and the coupled column's deep oxygen and bottom water's age at its end,
says for every column that did not settle which way its saved decade leans. In thin water the wind's share is spread
over little mass, so its dissipation, and the floor's ventilation, grows as the water below the mixed layer thins; the
3 m and 10 m pans at Procellarum's and Ingenii's typical cells and the 3 m pans of the western shallows and Drygalski
run again under the first round's fixed dissipation as the bracket. The plan is 90 columns.

**Bounds on the circulation.** A column carries no horizontal circulation. A sea spanning latitudes makes its densest
water at its coldest end, and that water sinks and fills its depths from below, as a filling box fills
(Baines & Turner 1969). The coldest column sets the floor to the deep water's temperature and the most its coldest
weather could renew the whole basin; the typical column sets the deep water's temperature from above and the renewal the
basin's ordinary weather gives alone. With exchange the deep water takes the cold end's temperature and only the cold
end, the share of the basin whose coldest 3-day means lie within 0.5, 1 or 2 K of the coldest cell's, renews it; the
screen stretches the coldest column's interval between renewals, from its saved decade's overturns as its oxygen is
read, with the coarse passes' rate beside it, by that share, and takes the deepest tenth's reference demand over the
added interval (half of it on average). The reference demand is the central columns', fixed by supply at any warmth; the
warm kinetics would raise it from the rates' 6 °C by ecosystem respiration's 0.62 eV, 5–7 times at the deep water's
25–29 °C. Whether connected basins share deep water is read from each basin's coldest water at its sill, the coldest the
coldest column holds at the sill's depth in its decade, against its neighbour's typical deep water.

**What a column can and cannot tell about a basin.** A one-dimensional column is a horizontally uniform basin under one
cell's weather. It can tell how deep that weather's nights and winds mix the water, how much of the floor they reach,
how fast the floor and the water use oxygen against what the surface takes up, the temperature the deep water settles at
when the whole basin feels one cell's coldest spell, and how often such weather would turn the whole column over. Run at
a basin's typical and coldest cells it bounds the deep water's temperature and the renewal its weather could give. It
cannot tell anything that moves water sideways: the cold end's dense water running down the slopes and along the floor;
the wind's stress, which the review found tilting the stratification to the surface even in a 581 km² crater lake (its
Lake number under 1 on 93% of days and its Wedderburn number under 1 on 99%), and which in basins many internal Rossby
radii across drives Ekman transport, upwelling along the shores and Kelvin and Poincaré waves; those waves' breaking on
the slopes, where deep lakes mix much of their deep water (Wüest & Lorke 2003); and the water spilling over sills
between basins. Below the nights' mixing, the column's exchange is its diffusivity, the wind's dissipation over the
water's stratification, and in the seas' nearly neutral deep water that reaches the cap of 0.1 m²/s across 33–69% of the
interfaces below the mixed layer under Osborn's rule (0.1–0.4% under the second form); there the cap sets the exchange.
A column's deep oxygen, its overturns and its deep water's age are therefore what the vertical physics alone gives each
cell's weather, and the seas' own deep water lies wherever the sideways circulation takes it between the bounds.

## Assumptions and limits

- **One column for a basin.** A column has no horizontal circulation, no sills inside it and no inflows (what that
  leaves a column able to tell is set out above); each basin's deep parts below its own inner sills (Imbrium and
  Serenitatis below 759 m, Nectaris below 411 m, Humorum below 391 m) are lumped into one hypsography. The screen of
  exchange scales the coldest column's renewal by its cold end's share and takes the deepest tenth's reference demand,
  which the warm kinetics would raise 5–7 times at the deep water's 25–29 °C, so its losses are the central columns' and
  lower ones than warm kinetics would give.
- **One record, repeated.** Each column passes the GCM's ten years 80 times, or 30 from its basin's central column. The
  deep water settles to the coldest its surface gets in those ten years, and its renewal turns on that cold recurring;
  a longer record holds colder spells.
- **The deep mixing.** The wind's share of its energy flux reaching the stratified water is Earth lakes' 0.3%, with its
  bracket of 0.2–0.4%, spread evenly over all the water below the mixed layer, which in a basin kilometres deep spreads
  it thin and in a pan's few metres concentrates it. Osborn's rule turns it into a diffusivity in every column but the
  two Shih–Bouffard rows; at the buoyancy Reynolds numbers of the seas' deep water the published form for energetic
  turbulence mixes less, so of the two forms Osborn's renews the deep water the faster. The seas' tides (near
  their equilibrium at periods of 14–32 days, with monthly ranges of 2–4 m in the largest seas; geography's tides) and
  the internal waves they and the wind drive are not in the column. Where the water below is nearly unstratified, the
  diffusivity reaches the cap of 0.1 m²/s, which then sets the exchange there.
- **Solutes.** Aerobic respiration's CO2 and the debt's anaerobic CO2 and methane, half its oxygen equivalents ending
  as methane, weigh on the density at Lake Kivu's coefficients, as in lake_oxygen; iron reduction's dissolved load,
  calcite's bicarbonate, the seas' own salts, the rivers' inflow and the groundwater through the floors are left out.
- **The weather at the seas.** 3-day means at T21 made hourly; the resort mixes two runs' physics. The GCM's surface
  pressure at the cells chosen is 1.163–1.207 atm against its sea-level pressure of 1.199–1.200 atm there: T21's
  truncated orography raises the model's surface beside high land, and saturation follows the cell's own pressure, as
  lake_oxygen takes it, up to 3% below what 1.2 atm would hold.
- **Earth's rates.** The trophic states' production and demand, the transfer velocities, the GLM coefficients, the
  wind's share and the labile pool's 2–7 days are Earth lakes' values or design values; the sediment's demand acts on
  every square metre of floor at every depth, as in a lake. Biology is undesigned; the 2 and 0.5 mg/L lines are the
  usual hypoxic and anoxic marks.
- **Pressure on the water.** The lake model's density is Kell's at one atmosphere, so the deep water's compressibility
  and the shift of its temperature of greatest density at 50–110 atm, the pressure 3–7 km down at lunar gravity, are
  left out; the deep water lies more than 20 K above that temperature.

## The checks

- **The construction.** At 300 K the lake kinds study's Setting301 gives climate_targets' own share, ladder and lapse
  rate; this extract holds lake_oxygen's waters (every water cell's surface temperature and pressure in both runs) and
  its Sun's phase identically; and lake_oxygen's saturation figures, the resort's 3-day-mean median of 7.688 mg/L and
  the seas' 7.751 among them, come out of this extract unchanged.
- **The seas.** The labelling gives the atlas's 28% share on the finer grid (27.998%) and its largest bodies' areas
  within 10%; the run basins stand for all 10.62 million km².
- **The runs.** Every column keeps its heat to 4.6 × 10⁻⁹ of the column's heat (1.5 × 10⁻⁸ W/m² over its surface)
  and its oxygen to 2.1 × 10⁻⁷ of each oxygen column's content, lake_oxygen's measures and within its tests' limits
  (10⁻⁸, 10⁻⁶ W/m² and 10⁻⁶). columns.py saved the heat error against the decade's net heat, which near equilibrium is
  almost nothing; run.py's heat_budget() rebuilds the net heat from the saved decade (its first day less the heat of
  the record's first step, run again from the decade's end) and so the error against the column's own heat. No sea
  freezes. lake_oxygen's 50 m slab, forced as the columns are, follows the GCM's own water at the typical cells of the
  four largest basins and the two polar ones in both runs: paired at the GCM's lowest level, as the columns take it,
  its daily temperature correlates at 0.94–0.98 with the GCM's and within 0.2 K of its mean with PlaSim's
  free-convection constant, and more closely than with the 10 m pairing at every cell (0.77–0.98); carried to lunar
  gravity the free convection leaves the slab 0.0–1.0 K warmer than the GCM's water. The Shih–Bouffard columns ran
  mixing_forms.py as its owner left it (sha256 9ba0a0bb…), and the index records the fingerprint of the code they
  ran.
- **The passes at the analysis step.** 44 of the 52 deep columns settle there. Of the eight that do not, the four with
  the Lake number's diffusivity or the fixed dissipation are cut off and ageing, anoxic throughout; Procellarum's under
  the second form turns over in the last of its ten passes only; Crisium's at the wind's high share and Ingenii's
  coldest are still losing oxygen and Fecunditatis's at the low share still gaining. The product's analysis_step names
  each and the way its decade's figures lean, and the seas-wide views name the unsettled columns they take.
- [test_run.py](test_run.py) checks the products' binding to their code and inputs, the seas against the atlas, the
  basins' hypsography, shelves, sills and the areas they stand for, the T21 grid, the extract's binding and a sea cell
  read at its own level, the 300 K reproduction, the 301 K share, saturation between the targets, the cells' rule, the
  budgets in lake_oxygen's measures, the spin-up and the settle test, the labels on the unsettled columns, the physics
  against lake_oxygen's and the form of each column's deep mixing, the other form's weaker mixing, the pans' bracket,
  the coldest column's densest water and the screen on the saved decade's overturns, the Rossby radius against its
  formula and the large basins' many radii, oxygen against the trophic state, the slab check and the index's binding of
  every saved column and of the other form's code.

## Running it

```sh
python -m research.studies.sea_oxygen.seas                    # about a minute: LOLA, the geoid and the atlas
/home/kir/Documents/Projects/terluna/climate/gcm/.venv/bin/python -m research.studies.sea_oxygen.extract  # netCDF4
python -m research.studies.sea_oxygen.run                     # the product from the cached columns; runs no column
python -m research.studies.sea_oxygen.run --columns           # the columns, on two cores (hours; resumable)
python -m pytest research/studies/sea_oxygen
```

The extract writes `research/runs/sea_oxygen/forcing.npz` with `forcing.json` (ignored, on the research drive),
recording the sha256 of every file it read. The columns are cached in `research/runs/sea_oxygen/cases`, each spin-up
saved after every pass, and `research/runs/sea_oxygen/columns.json` indexes them with their sha256s for the studies that
read them, the lake methane study's tracer first; each case's `diffusivity` names the form of its deep mixing. The cache
key holds this study's column code and lake.py whole, the source of what the columns use from lake_oxygen, the lake
kinds study and mixing_forms.py (so that an edit elsewhere in those files leaves the columns standing), and each input's
sha256 (so that an input written again, even with the same figures, re-keys them). The 90 columns take about four and a
half hours on two cores.
