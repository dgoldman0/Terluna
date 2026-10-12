# Living light on the Open Moon

**11 October 2026.** The author's direction of that day: "Lighting does not have to work the same way as on Earth for
many areas because many areas get twilight at darkest. We can also rely more heavily on bioluminsecence" (the
[register](../../decisions.md), People on the Moon and beyond, 11 October 2026; the second sentence also under
Welkinia, strange taxa and living light). This study asks three things of it: what living
light can give and what it costs; what people need outdoors at night in each setting, against that setting's own night;
and how much land people's light puts under a light-polluted sky once twilight, Earthlight and living light are counted.

[run.py](run.py) writes [results/living_light.json](results/living_light.json) (schema
`terluna.research.living-light/1`) in about 20 seconds. Numbers marked *screen* are computed there; *derived* marks
arithmetic on cited values. Every source is in [sources.json](sources.json) with what was read. Design guesses are
ranges, labelled once where they enter: the share of a glow that reaches the path in section 1, the lighting guesses in
section 3.

## Findings

1. **Today's brightest glowing plants mark a path in the far side's darkest hours and give a hundredth of the light
   needed to walk by.** The brightest measured plants, tobacco flowers carrying the enhanced fungal pathway, emit
   5.0×10¹³ photons a second per m² (Zheng et al. 2023): 0.0026 cd/m², or 0.0072 cd/m² to dark-adapted eyes. A closed
   canopy glowing that brightly would light the path beneath at 0.008 lux (*screen*), against the 1 lux at which walkers
   of every age see a 10-mm trip hazard (Fotios & Uttley 2018) and the 0.05 lux of a summer full moon (Kyba, Mohar &
   Posch 2017). The improved pathway of 2024, ten to a hundred times brighter (Shakhova et al. 2024), would give
   0.018–0.18 lux (*derived*). A waymark at 0.0026 cd/m² stands out from the ground wherever the darkest light stays
   under 0.023 lux: over the summit metropolis and on 22–30% of the other low-latitude homes' ground, and nowhere under
   the Earth or in the twilit middle and high latitudes (*screen*). Bacteria already glow more brightly over whole seas:
   a milky sea seen from orbit shines at 0.04 cd/m² or more over 15,400–100,000 km² (Miller et al. 2005, 2021).

2. **The chemistry and the glow's direction set the cost: firefly chemistry facing down lights a path through the
   summit's dark hours for 2–10% of a plant's night store, and the fungal pathway as measured for a third of the store
   to thirty times it.** A closed glowing canopy that lights its path at 1 lux through the 165 hours a cycle the
   summit's night stays below that level takes 0.78–0.93 g C/m² with firefly chemistry, quantum yield 0.41–0.48 (Ando et
   al. 2008; Niwa et al. 2010), when 90% of its light reaches the path. That is 1.6–9.6% of the canopy model's night
   store: 44–49 g C/m² for an Earth-like plant and 9.6–11 g C/m² for one that idles through a quarter of the cycle. Glow
   made inside green leaves, a tenth of it reaching the path, costs nine times as much: 7.0–8.3 g C/m², 14–87% of the
   store. Bacterial chemistry facing down takes 1.1–2.0 g C/m². The fungal pathway takes 15–275 g C/m² facing down,
   because its quantum yield is 0.00096–0.017 (Teranishi 2018) (*screen*). Per photon made, the fungal cycle spends
   960–17,800 ATP equivalents, firefly chemistry 60–74 and bacterial chemistry 39–70 (*screen*).

3. **At the physical limits light is cheap; quantum yield, luciferin supply and light lost in the leaf set how brightly
   a plant glows.** A tenth of the night store, spent through a whole 354-hour night, lights the ground under a glowing
   canopy at 0.48–2.9 lux with firefly chemistry facing down, 0.69–3.6 lux with the best beetle luciferase (quantum
   yield 0.61; Niwa et al. 2010) and 13–64 lux at the thermodynamic floor, one photon per luciferin with nothing lost
   (*screen*). Glow that loses nine-tenths of its light in the leaf gives a ninth as much. One micromole of firefly
   luciferase per m², 60 mg of protein, running at its measured 0.2 photons a second a molecule (Klein et al. 2023),
   would give 24 lux under a closed canopy with all its light reaching the path, 2.4–22 lux at the 0.1–0.9 share
   (*derived*). A designed plant with a beetle-class luciferin cycle, its glow facing down, can light a path at walking
   level for 2–10% of its night store; the fungal pathway's yield is the gap.

4. **The natural night already serves walking on half the Moon's dry land.** Twilight and Earthlight keep the night's
   darkest moment above 1 lux on 48% of the dry land, above a summer full moon's 0.05 lux on 80%, above the 1.6 lux that
   lets a face be recognised at 4 m on 40% and above 8 lux for reading on 21% (*screen*). Poleward of 48.6° twilight
   keeps walking light and a face's light at 4 m all night, and reading light from about 54°. Under an Earth 50° and
   more up the darkest moment keeps 1.5–4.2 lux, enough for walking, and with the Earth 60° and more up 2.4–4.2 lux,
   enough to know a face at 4 m. The far side's low latitudes are the dark belt: over the summit metropolis the night
   stays below 1 lux for 165 hours a cycle and below 0.05 lux for 102.

5. **Lamps remain for faces, reading and work in every setting, and for walking in the dark belt where living light is
   not grown; they need add only what the natural night lacks.** Over the summit that is 1 lux for walking, 1.6–5.4 lux
   at gathering places so faces are known at 4–10 m (Caminada & van Bommel 1980, as Fotios, Yang & Cheal 2015 report it)
   and 8–100 lux of task light (*screen*). Under an Earth 50° and more up walkways need no lamps and reading needs
   3.8–6.5 lux more; poleward of 48.6° walkways and faces at 4 m need none. Twilight-aware lamps burn through the hours
   below each level: for walking, 165 hours a cycle over the summit, 84–104 on average in the other low-latitude homes,
   16 in the nearside towns and none in the high latitudes (*screen*). Living light's colours suit the eye at these
   levels: at 1 lux on the ground the fungal pathway's green counts 1.6 times its photopic value in mesopic vision,
   bacterial blue-green 2.4 times and firefly yellow-green 1.1 times (CIE 191:2010's MES2 as Goodman 2013 sets it out;
   *screen*).

6. **The Moon's tall, dense air carries a town's skyglow farther than Earth's.** Above each square metre stands 7.8
   times Earth's air, an optical depth of 0.76 at 550 nm against 0.097. In single scattering a shielded town's zenith
   skyglow is 1.2–1.5 times Earth's at 1 km, 2.3–2.9 times at 5 km, 4.9–9.2 times at 20 km and 12–53 times at 100 km.
   Between 10 and 50 km it falls with distance to the power −2.0 on the Moon, against −2.6 in Earth's clear air and −3.0
   in its hazy air (*screen*; Walker 1977 measured −2.5 and Garstang 1986 modelled −2.3).

7. **Bright natural skies remove the polluted sky round ordinary towns on two-thirds of the dry land.** Falchi et al.
   (2016) count a sky as light-polluted where artificial light adds more than 8% to the natural zenith. At its darkest
   the Moon's natural sky is 1.1–7.8 times Earth's on the far side within 15° of the equator, 6–106 times at 15–30°,
   85–205 times on the limb at low latitudes, 100–3,200 times on the far side and the limb at 30–48.6° and 230–5,300
   times wherever the Earth is up; poleward of 48.6° twilight makes it 3,100 to a million times (*screen*). Towns of up
   to a million people light their sky above Falchi's level on 33% of the dry land with shielded lamps, against 84% when
   they light it as affluent Earth does (*screen*).

8. **The land under a light-polluted sky falls from ecology's 14–49% of the dry land to 0.9–40% at 4–8 billion lunar
   residents.** With shielded lamps a person lights 23–1,730 m² in the low-latitude homes and 0–36 m² in the Earthlit
   and twilit ones, against ecology's Earth-sky 791–1,420 m² (*screen*). Spread with the land, 4, 6 and 8 billion light
   1.0–17%, 1.2–24% and 1.3–31% of the 23.4 million km² of dry land; denser on the waterfronts or in the rain belt, up
   to 23%, 33% and 40%; on the cool ground of the highlands and higher latitudes, 0.9–11%, 1.0–15% and 1.1–19%
   (*screen*). Over the Moon that is 31–1,340 m² a resident. Town size sets most of the range: in towns of 10,000 the
   spread case gives 1.0–9.1% across 4–8 billion, in cities of a million 2.4–31%. The summit metropolis lights a circle
   230–470 km in radius, 0.7–3.0% of the dry land, when it covers Earth's area a person; when each of its people lights
   the sky as a resident of a town of 10,000 does, 120–180 km and 0.2–0.4%, and the totals become 0.4–38%. Living light
   that takes over 30–70% of outdoor lighting brings the low ends to 0.5–0.6% (*screen*).

9. **The long twilight shields most of a town's night.** A town of 100,000 on the far side's equator lights its sky
   above Falchi's level for 76–108 of the night's 354 hours, and over the night its polluted area averages 9–12% of its
   darkest-hour extent; the metropolis at the summit's latitude for 208–292 hours, averaging 19–26% (*screen*).

10. **The dark belt holds Earth-like darkness and the luminous organisms' niche, and near towns direct light matters
    more than skyglow.** The night falls below the 0.01 lux at which light suppresses melatonin in fishes and rodents
    (Grubisic et al. 2019) on 13% of the dry land, 3.1 million km²: for 66 hours a cycle at the summit and 9–15 hours on
    average over the other low-latitude homes (*screen*). Skyglow alone brings melatonin-level light to the ground
    within 53–229 km of the metropolis's centre and just past the edge of a city of a million in hazy air (*screen*).
    Lit hedgerows held 47% fewer caterpillars (Boyes et al. 2021), and at about 5 lux an unshielded LED caught 5.4 times
    a dark control's insects where a shielded one caught 1.3 times (Dietenberger et al. 2024). Glow-worm mate attraction
    may fall at 0.02–0.1 lux (Van den Broeck et al. 2021), so luminous signalling belongs to the dark belt and to forest
    floors, which keep 4–5×10⁻⁵ lux in the dark belt and 0.013–0.19 lux under the Earth (*screen*). The ring fleet's
    rolled night half, about 0.02 lux at a far-side midnight, would light the dark belt at 2.5–18 times its natural
    darkest (*screen*).

## 1. What living light can give

### Light from organisms and engineered plants

A glowing surface emits light into the half-space above it; its luminance is its photopic output over π. What it gives a
path depends on how much of the sky it fills: a closed canopy overhead gives the path its full exitance, hedges 1 m tall
on both sides of a 2-m path 29% of it (configuration factors, *derived*). Emission spectra are Gaussian stand-ins with
the cited peaks (fungal 525 nm, firefly 560, bacterial 490) and assumed widths; the 132 luminous fungi known glow at
520–530 nm (Perry, Desjardin & Stevani 2025).

| Source | Photons s⁻¹ m⁻² | Luminance, cd/m² (dark-adapted) | Path under a closed canopy | Brighter needed to walk by | 354-hour night, g C/m² |
|---|---:|---:|---:|---:|---:|
| *N. benthamiana* flowers, fungal pathway, 2020 (Mitiouchkina et al. 2020) | 1.1×10¹³ | 0.00056 (0.0016) | 0.0018 lux | 570 times | 0.057–1.0 |
| *N. tabacum* flowers, 2020 | 5.3×10¹² | 0.00028 (0.00076) | 0.0009 lux | 1,200 times | 0.028–0.51 |
| Leaves, 2020 | 9.3×10¹¹ | 0.00005 (0.00013) | 0.00015 lux | 6,600 times | 0.005–0.09 |
| Enhanced flowers, 2023 (Zheng et al. 2023) | 5.0×10¹³ | 0.0026 (0.0072) | 0.008 lux | 120 times | 0.27–4.8 |
| Improved pathway, 2024, 10–100 times 2020 (*derived*) | 1.1×10¹⁴–1.1×10¹⁵ | 0.0056–0.056 (0.016–0.16) | 0.018–0.18 lux | 6–57 times | 0.57–104 |
| Milky sea, luminous bacteria (Miller et al. 2005) | 1.4×10¹⁵ | 0.040 (0.22) | 0.12 lux | 8 times | 0.31–0.53 |

Values are *screen* except where marked. The night's cost takes 90% of the glow to reach the path, at the fungal
pathway's quantum yields and the milky sea's at bacterial ones; the share that reaches the path is set out under what it
costs. The 2020 plants' glow followed their clock and kept on through four days of darkness, its rhythm switched off on
the third and fourth (Mitiouchkina et al. 2020); the enhanced plants' glow fell slowly in the dark and rose 2.9 times
when they were fed sucrose (Zheng et al. 2023). Luminous fungi's own glow, glow-worms' and firefly lanterns' steady
output have no measured luminance in physical units in the sources reached.

Other organisms give light in flashes or from single cells (*screen*, from the per-cell values):
- **Dinoflagellates.** A *Pyrodinium bahamense* cell flashes at 2.9×10⁸ photons a second (Letendre et al. 2025): a flash
  at 1 cd/m² needs 2.7×10⁸ cells under each square metre flashing together, at 0.01 cd/m² 2.7×10⁶. *Noctiluca
  scintillans*, 9×10¹⁰ photons a second (Letendre et al. 2024), needs 870,000 cells for 1 cd/m².
- **Bacteria.** At 10³–2×10⁴ photons a second a cell (Pérez et al. 2010; Dunlap 1985), a milky sea's glow takes
  7×10¹⁰–1.4×10¹² cells under each square metre.
- **Fireflies.** A *Photinus pyralis* flash holds about 10⁸–10¹¹ photons (Silver 2026, as press reports give it). Over a
  lantern of 0.7–3.5 mm², the size of *Photinus palaciosi*'s (López-Palafox et al. 2020), one flash holds 3–12,500
  seconds of the 2020 flower's glow from each square millimetre (*derived*).
- **Other luciferins.** Beetle luciferases and the marine imidazopyrazinone systems come nearest firefly chemistry's
  yield: NanoLuc 0.45 (Russo, Civili & Winssinger 2024, citing Machleidt et al. 2015), *Cypridina* 0.28 and aequorin
  about 0.17 (Niwa et al. 2023), *Renilla* 0.053 (Loening, Dragulescu-Andrasi & Gambhir 2010); the red-glowing railroad
  worm's luciferase gives 0.15 (Niwa et al. 2023).

### What it costs

Each luciferin cycle's turnover is costed from its reactions in ATP equivalents, at P/O ratios of 2.3–2.5 for NAD(P)H
(Hinkle 2005), and divided by the quantum yield (*screen*):
- **Fungal caffeic acid cycle** (Kotlobay et al. 2018): hispidin synthase adenylates caffeic acid and adds two
  malonyl-CoA, made at acetyl-CoA carboxylase from the cytosol's acetyl-CoA, which ATP-citrate lyase supplies; the
  hydroxylase spends NADPH, 30% of it uncoupled (Tong et al. 2020); the luciferase releases CO₂ and light; the hydrolase
  returns caffeic acid and pyruvate. 16.3–17.1 ATP equivalents a turnover; 960–17,800 a photon. Yu et al. (2026) count
  the same cycle as one ATP, one NAD(P)H and two malonyl-CoA.
- **Firefly** (Marques & Esteves da Silva 2009; Gomi & Kajiyama 2001): ATP to AMP, and oxyluciferin back to luciferin
  with one cysteine, 26.8–28.5 ATP equivalents by the sulfate and serine pathway with its nitrogen from glutamine
  synthetase and GOGAT (*derived*). 28.8–30.5 a turnover; 60–74 a photon. Whether fireflies recycle their luciferin is
  unproven (Oba et al. 2013); made new from two cysteines each time it costs about twice as much.
- **Bacterial** (Brodl et al. 2018): reduced flavin from NAD(P)H, and the aldehyde regenerated with ATP to AMP and
  NADPH. 6.6–7.0 a turnover; 39–70 a photon at quantum yields of 0.10–0.17 (Shimomura, Johnson & Kohama 1972; Kaku et
  al. 2021).
- **The floor**: 0.45–0.52 mol of carbon per mole of photons, the photon's energy paid at FAO's 15.7 kJ per gram of
  sugar.

The costs leave out making the enzymes and credit nothing for the two-carbon fragment released when luciferin is
recycled. Each then scales with the share of the glow that reaches the path, a design guess of 0.1–0.9: up to 0.9 for
glowing cells on leaves' undersides beneath a reflecting layer, 0.1 for glow made inside green leaf tissue, which
absorbs much of it and sends the rest both ways.

A closed glowing canopy lighting its path at each need's level, through the hours the summit's night stays below it, its
glow facing down (*screen*):

| Need on the path | Hours at the summit | Fungal pathway | Firefly chemistry | Bacterial chemistry | Floor |
|---|---:|---:|---:|---:|---:|
| Orientation, 0.05 lux | 102 | 0.47–8.5 g C/m² | 0.024–0.029 | 0.035–0.061 | 0.001 |
| Colour, 0.32 lux | 141 | 4.2–75 | 0.21–0.25 | 0.31–0.54 | 0.010 |
| Walking, 1 lux | 165 | 15–275 | 0.78–0.93 | 1.1–2.0 | 0.035 |
| A face at 4 m, 1.6 lux | 177 | 26–472 | 1.3–1.6 | 2.0–3.4 | 0.060 |
| A face at 10 m, 5.4 lux | 204 | 101–1,840 | 5.2–6.2 | 7.7–13 | 0.24 |
| Reading, 8 lux | 213 | 156–2,830 | 8.0–9.5 | 12–20 | 0.36 |
| Walking through the whole night | 354 | 33–590 | 1.7–2.0 | 2.5–4.2 | 0.075 |

The columns take 90% of the light to reach the path; at a tenth every chemistry costs nine times as much, the floor
unchanged. The canopy model's night store is 44–49 g C/m² for an Earth-like plant and 9.6–11 g C/m² for one that idles
through a quarter of the cycle, and its growth 203–256 g C/m² a cycle. At the brightest measured glow a night costs
0.27–4.8 g C/m² facing down, 0.008–0.15 W through the night; the register's lasting glow (M09) puts a glowing square
metre at about 0.04 W (*derived*).

Glowing raised a glow-worm larva's CO₂ output by a quarter, from 3.49 to 4.38 µl an hour (Willis, White & Merritt 2011),
and flashing raised a firefly's metabolic rate 37% over basal (Woods et al. 2007, as López-Palafox et al. 2020 report
it).

## 2. What people need outdoors at night

The needs are Earth's standards and experiments, levels of light on the ground:

| Need | Level | Basis |
|---|---:|---|
| Orientation: moving about open ground by eye | 0.05 lux | A summer full moon at temperate latitudes, 0.05–0.1 lux (Kyba, Mohar & Posch 2017) |
| Telling colours | 0.32 lux | Paper samples sorted by colour above 0.32 lux; at 0.1–0.01 lux only reds and oranges named (Pokorny et al. 2006) |
| Walking safely | 1 lux | Trip hazards seen by walkers of all ages under any lamp (Fotios & Uttley 2018); the lowest pedestrian class asks 2 lux on average, 0.4 at least (CIE 115:2010, in Fotios & Gibbons 2018) |
| A face at 4 m and at 10 m | 1.6 and 5.4 lux | Semi-cylindrical 0.8 and 2.7 lux (Caminada & van Bommel 1980, in Fotios, Yang & Cheal 2015), doubled to the ground under a uniform sky (*derived*) |
| Reading ordinary print | 8 lux | Full reading speed at 2 cd/m² (Peñaloza et al. 2025) on paper of reflectance 0.8 (*derived*) |
| Close work outdoors | 100 lux | Reading addresses outdoors (EN 12464-2); 500 lux for indoor reading and writing (EN 12464-1) |

Night vision sets the rest. The rods take over after 5–8 minutes in the dark and reach their threshold, about 10⁻⁵
cd/m², after about 40 (Kalloniatis & Luu 2007). Below 0.005 cd/m², the foot of CIE 191:2010's mesopic range (Crumey
2014) and 0.087 lux on 18% grey, vision is the rods' alone and colourless; between 0.005 and 5 cd/m² rods and cones
share it, and green and blue-green light count for more than their photopic value (finding 5). People feel about as safe
as by day at about 30 lux (Boyce et al. 2000). The face levels assume a uniform sky, which lights a face with half the
ground's light; under a high Earth most light comes from above, a face gets less, and the verdicts for faces lean
generous there.

The natural night is the calendar's evaluator on its solved clear sky in the mean geometry, the same calculation as
provisioning's daily-life study, which it reproduces at all 22 of that study's places and latitude rows, and
main's per-place Earthlight at all seven of main's places (*screen*; tests). The light is read on level ground,
without libration, terrain, cloud or haze; ordinary hazes move practical dusk by under an hour (illumination/aerosol).
Hours below each need's level a cycle (*screen*):

| Place | Darkest | Orientation | Walking | A face at 4 m | Reading |
|---|---:|---:|---:|---:|---:|
| Summit metropolis, far side 5.4° N | 0.0014 lux | 102 h | 165 h | 177 h | 213 h |
| Far side 15° | 0.0057 lux | 87 h | 158 h | 170 h | 208 h |
| Far side 30° | 0.090 lux | 0 | 123 h | 140 h | 187 h |
| Ingenii coast, far side 34° S | 0.20 lux | 0 | 106 h | 126 h | 178 h |
| Eastern Smythii headland, the limb | 0.061 lux | 0 | 164 h | 177 h | 213 h |
| Procellarum by Russell, the Earth 12° up | 0.42 lux | 0 | 99 h | 121 h | 179 h |
| Highland box, far side 44.7° S | 1.44 lux | 0 | 0 | 31 h | 136 h |
| Southern Mare Nubium, the Earth 61° up | 3.1 lux | 0 | 0 | 0 | 186 h |
| Tranquillity's harbour, the Earth 66° up | 2.4 lux | 0 | 0 | 0 | 204 h |
| Polar plain, 75° N | 486 lux | 0 | 0 | 0 | 0 |

Over the ways-of-living homes, read at five latitudes in each band and by where the Earth stands (*screen*; hours are
means over each home's land). Within each home the latitude bands and the Earth's classes are taken as independent, and
poleward of 60° the class with the Earth 30° and more up is read at the Earth's highest elevation there. The classes are
weighted by the cosine of the Earth's elevation, which at low latitudes gives the 30° readings more weight than their
share of each latitude row; weighted by that share instead, the land never below the walking level rises from 48% to
49–51% and the land never below the level for faces at 4 m from 40% to 42–44%:

| Home | Dry land | Darkest | Never below 1 lux | Hours below 1 lux | Never below 8 lux | Hours below 8 lux |
|---|---:|---:|---:|---:|---:|---:|
| Summit metropolis | 2,464 km² | 0.0014 lux | 0 | 165 | 0 | 213 |
| Coastal and lake towns | 3.95 million km² | 0.0011–726 lux | 32% | 92 | 11% | 169 |
| Wetland-edge towns | 3.58 million km² | 0.0011–726 lux | 26% | 102 | 7% | 182 |
| Canopy towns and farming country | 4.79 million km² | 0.0011–3.7 lux | 20% | 104 | 0 | 201 |
| Highland towns and elevated districts | 1.99 million km² | 0.0011–726 lux | 36% | 84 | 17% | 148 |
| Nearside towns | 1.76 million km² | 0.38–5.2 lux | 79% | 16 | 0 | 152 |
| Limb towns | 0.57 million km² | 0.20–2.2 lux | 40% | 56 | 0 | 159 |
| Far side's twilit-night places | 1.75 million km² | 0.14–2.1 lux | 40% | 58 | 0 | 159 |
| High-latitude towns | 5.00 million km² | 4.2–726 lux | 100% | 0 | 80% | 17 |

Natural light never falls below the walking level on 48% of the dry land: all the high latitudes, 79% of the nearside
towns' land, 40% of the limb and far-side middle latitudes and 20–36% of the low-latitude homes. Under the Earth the
darkest moment depends on its height and on the place's longitude: 0.6–3.4 lux with the Earth 30° up at low and middle
latitudes, 1.5–3.9 lux with it 50° up and 2.4–4.2 lux with it 60° and more up. Poleward of 48.6° the darkest light is
4.2–6.9 lux at 50.8° and 10–13 lux at 55.6° (*screen*).

**Where each kind of light serves.**
- **Natural twilight and Earthlight alone**: walking, colour and faces at 4 m poleward of 48.6°, with reading from about
  54°; orientation everywhere outside the far side's low latitudes; walking under an Earth 50° and more up, and faces at
  4 m with it 60° and more up.
- **Living light**: waymarks for orientation wherever the night grows darker than the waymark; at today's measured
  brightness that is the summit and 22–30% of the other low-latitude homes' ground, and at the 2024 pathway's
  hundredfold also 56–63% of the low-latitude homes' ground and 40% of the middle latitudes' far side and limb
  (*screen*). Walking light from glowing canopies needs a beetle-class luciferin cycle facing down: 2–10% of the night
  store for the summit's 165 hours.
- **Lamps**: faces at gathering places, reading and work everywhere; walking wherever living light is not grown in the
  dark belt. They add what the natural night lacks at its darkest: over the summit 1 lux for walking, 1.6 and 5.4 lux
  for faces, 8 lux for reading and 100 for close work; under an Earth 50° and more up nothing for walking, up to 0.12
  lux for faces at 4 m, 1.2–3.9 lux for faces at 10 m and 3.8–6.5 lux for reading; poleward of 48.6° nothing for walking
  or faces at 4 m (*screen*).

## 3. People's light reaching the sky

### The skyglow model

A town is a uniformly lit disc whose fully shielded lamps send light upward only by reflection from the ground, a
Lambertian source (Garstang 1986's G term; full cut-off fixtures emit no uplight, Kinzey et al. 2017). The light
scatters once towards each observer's zenith, on a sphere, through the design run's global-mean air at the solved
column's cross-section, and through Earth's 1976 standard air. Each world is read clear and hazy: the Moon with the
design run's central night haze, an aerosol depth of 0.006 (illumination/aerosol), and Earth with the haze of Falchi et
al.'s atlas (clarity K = 1, 0.33 magnitudes of vertical extinction: an aerosol depth of 0.21, *derived*; its 1.5-km
scale height, asymmetry of 0.7 and single-scattering albedo of 0.95 are assumed). Multiple scattering, which adds light
far from towns in a column of depth 0.76, is left out, so the Moon's distances lean short. On Earth the model's curve
falls with distance to the power −2.6 between 10 and 50 km in clear air, and in hazy air −2.9 between 10 and 35 km and
−3.0 between 10 and 50; Walker (1977) measured −2.5, and Garstang (1986) found −2.3 for his standard haze between 10 and
35 km.

Each town is calibrated on Earth to ecology's Earth-sky area a person (791–1,420 m² with shielded lamps, a tenth of
Earth's 7,910–14,200 m², from Falchi et al. 2016 over the World Bank's land and people), as an isolated town of 10,000,
100,000 or a million people at 10,000 a km², and the metropolis's 100 million at 40,000 a km². Earth practice is the
same calibration at Earth's areas, with the same ground-reflected (Lambertian) pattern. On the Moon a town's polluted
zone is where its skyglow exceeds the same share of each place's natural sky at the night's darkest moment. Holding the
area a person fixed sets each size's upward light a person: towns of 100,000 and a million emit 0.6–1.9 times a town of
10,000's light a person, and the metropolis 10–230 times (*screen*), so the metropolis is also read with its people
lighting the sky as a small town's do.

**The natural night sky.** Twilight's zenith is 0.4–0.8 of the night sky's mean luminance: the solved sky gives 0.82 at
6° of depression and 0.65 at 30°, and is unsampled deeper; if the share keeps falling in the dark belt's deepest night,
its polluted zones grow. The Earth's scattered light takes the solved sky's ratio at its height. Starlight, zodiacal and
galactic light and airglow, which the light products omit, add 57–174 µcd/m²: at most Earth's natural 174 µcd/m² (Falchi
et al. 2016), at least its extraterrestrial 52% (Duriscoe & Hung 2016) dimmed by the Moon's deeper air, with no airglow,
whose lunar counterpart is unknown (*derived*).

**Design guesses** (ranges, *guess*): twilight-aware lamps add only what the natural light lacks of the 1-lux walking
level; windows, vehicles and task lights that no outdoor rule governs keep 10–50% of the upward light (street lights
give 13% of a city's light after midnight on Earth; Kyba et al. 2021); living light takes over 30–70% of the outdoor
lighting with 1–10% of the shielded lamps' upward light for the same service, its glow facing down from canopies and
hedges.

### Per person and over the dry land

At the night's darkest moment, m² of land under a light-polluted sky a person, shielded lamps (*screen*):

| Home | Dry land | Per person | Towns of 10,000 | Cities of a million |
|---|---:|---:|---:|---:|
| Coastal and lake towns | 3.95 million km² | 34–1,410 m² | 34–310 | 180–1,410 |
| Wetland-edge towns | 3.58 million km² | 43–1,730 | 43–390 | 220–1,730 |
| Canopy towns and farming country | 4.79 million km² | 23–1,090 | 23–220 | 130–1,090 |
| Highland towns and elevated districts | 1.99 million km² | 28–1,260 | 28–260 | 155–1,260 |
| Far side's twilit-night places | 1.75 million km² | 0–36 | 0 | 0–36 |
| Limb towns | 0.57 million km² | 0–6 | 0 | 0–6 |
| Nearside and high-latitude towns | 6.76 million km² | 0 | 0 | 0 |
| Summit metropolis, at Earth's area a person | 2,464 km² | 1,600–6,930 | | |
| Summit metropolis, lit like a town of 10,000 | | 480–990 | | |
| Residents aloft, over the whole Moon | | 13–540 | | |

Twilight-aware lamps barely change these darkest-hour areas: where the night is darkest the lamps run fully, and where
twilight would dim them the natural sky already hides their glow. Their gain is in the hours lamps burn and the light on
the ground. Living light lowers the low-latitude homes to 4–1,530 m² a person (*screen*).

Share of the dry land under a light-polluted sky, by provisioning's mixes of where people live (*screen*):

| Lunar residents | Spread with the land | Waterfronts | Rain belt | Cool ground | Highlands | Spread, half aloft | Spread, living light |
|---|---:|---:|---:|---:|---:|---:|---:|
| 4 billion | 1.0–17% | 1.2–23% | 1.1–21% | 0.9–11% | 1.0–18% | 0.95–14% | 0.5–15% |
| 6 billion | 1.2–24% | 1.4–33% | 1.3–31% | 1.0–15% | 1.2–26% | 1.1–20% | 0.5–21% |
| 8 billion | 1.3–31% | 1.7–40% | 1.6–38% | 1.1–19% | 1.4–31% | 1.2–26% | 0.5–28% |

The 4–8 billion lunar residents are the author's working range (the register, People on the Moon and beyond). The
mixes are the ways-of-living study's: people in proportion to each home's land, or four times as dense on the
waterfronts, in the rain belt or on the cool ground of the highlands and higher latitudes, with the metropolis at its
100 million; half aloft puts 2–4 billion in sky towns over the whole Moon (the register). The
author's direction of 11 October puts land residents mostly in the drier highlands, the summit metropolis round its
tower above all (the same register); four times as dense on the highland homes, a mix derived here, that gives 1.0–18%,
1.2–26% and 1.4–31%, close to the spread case, since 86% of the highland homes' land lies where the Earth never rises
and 43% within 30° of the equator (ways of living). The low ends are towns of 10,000 under the brightest natural floor;
the high ends cities of a million, an airglow-free sky and Earth's hazy air as the reference. Overlapping zones are
summed and capped at each home's land. The metropolis's circle counts as dry land whole, while standing water lies 15–24
km from the summit; the summit stands 9.45 km above the sea, under less air than the global-mean column read here (ways
of living). With the metropolis lit like a town of 10,000 every shielded share falls by 0.5 points at its low end and
2.5 at its high end. Lit as affluent Earth lights itself, the same people would light 10–83% of the dry land (*screen*).

## 4. The dark-night requirement and night life

**The dark belt is the Moon's dark-night reserve.** Natural light falls below 0.01 lux on 13% of the dry land, 3.1
million km², in the far side's low latitudes: averaged over the low-latitude homes for 9–15 hours a cycle below 0.01 lux
and 19–26 below 0.03, and over the summit for 66 and 92 hours, with the darkest light 0.0011–0.008 lux (*screen*).
Elsewhere twilight or Earthlight keeps every night above the 0.01–0.03 lux at which light suppresses melatonin in fishes
and rodents (Grubisic et al. 2019), and under an Earth 30° and more up above the 0.05–0.3 lux of Earth's full moon
(Kyba, Mohar & Posch 2017); organisms there need raised thresholds, as ecology's finding 5 has it. The dark belt is also
where the metropolis stands.

**The dark-night requirement.** The fleet decision keeps the ring fleet's light off the night zones and leaves the
requirement's value to the biosphere and human work (the register). Rolled near 80° with a Moon-facing scatter
near 10⁻⁵, the fleet's night half gives about 0.02 lux at a far-side midnight (the same register, from the joint
synthesis): 2.5–18 times the dark belt's natural darkest, and above the melatonin levels. Held within Falchi's 8% of the
dark belt's natural darkest, the fleet's light there stays under 0.00009–0.0006 lux; under the Earth 8% of the natural
darkest is 0.03–0.4 lux (*derived*).

**Outside the metropolis skyglow adds little to the night animals meet.** At the polluted boundary the artificial sky
adds 8% of the natural sky's light: 0.00006–0.0014 lux on the dark belt's ground, below the melatonin levels, and up to
0.04 lux where twilight keeps the natural night above them (*derived* from the natural zenith). Reaching 0.01–0.03 lux
by skyglow alone takes a sky 91–550 times Earth's polluted threshold: towns of 100,000 never reach it, a city of a
million within 6.3 km of its centre in hazy air, just past its edge, and the metropolis within 53–229 km of its centre
(*screen*). The sky near a town is brightest low towards it, so ground light taken from the zenith understates these
radii. Direct light is what animals meet: lit hedgerows at 1.4–16 lux held 47% fewer caterpillars and lit grass margins
33% fewer (Boyes et al. 2021); blackbirds at 0.3 lux developed their reproductive system up to a month early (Dominoni,
Quetting & Partecke 2013); LED light cut slow-flying bats' activity even at 3.6 lux (Stone, Jones & Harris 2012);
effects appear below 1 lux across taxa (Sanders et al. 2021); night visits to flowers fell 62% under street lights (Knop
et al. 2017). Flowering responses to night light saturate at under about 3.5–28 lux (Bennie et al. 2016, after Whitman
et al. 1998), so the near side's Earthlit nights, 0.4–5 lux at their darkest, lie within the range that steers Earth
plants' flowering.

**Insects drawn to light.** Insects turn their backs to a light and circle it, mostly within half a metre (Fabian et al.
2024), from within about 23 m of a street light (Degen et al. 2016) and often under 10 m of a weak one (Truxa & Fiedler
2012). Ultraviolet and blue draw most: at equal radiant flux UV took 84% of moths against blue, green and white, and
with UV replaced by red, blue took 63% while green, white and red did not differ (Brehm et al. 2021); metal-halide and
mercury lamps draw three and five times the moths of high-pressure sodium (Boyes et al. 2021, review), and amber and red
light draw relatively few insects (Owens et al. 2020). Shielding mattered more than brightness: an LED at 30.6 lux
caught 7.7 times a dark control, dimmed to 5.6 lux 5.4 times, and a shielded LED at 5.2 lux 1.3 times (Dietenberger et
al. 2024). Moth numbers fell 14% over three years under white, green and red light alike (van Grunsven et al. 2020). Dim
living light still draws some insects: traps lit like a glowing mushroom caught 42 insects against 12 in dark traps over
five nights (Oliveira et al. 2015), while the ghost fungus's steady glow drew no more than dark controls (Weinstein et
al. 2016). White light's blue also brightens the skyglow the dark-adapted eye sees: replacing sodium lamps with white
ones raises it 2.5–5 times (Falchi et al. 2011). Living light for people does best facing down, green or amber; living
light for lures can use the draw.

**Luminous organisms in the night food web.** Glows signal where nights are dark: glow-worm females start courting below
1.3 lux and are never active above 15 (Dreisig 1975, as Owens & Lewis 2018 report it), and mate attraction may fall at
0.02–0.1 lux (Van den Broeck et al. 2021); red light spares fireflies best and green, yellow and amber disturb them most
(Owens et al. 2022). Glows lure prey: glow-worm larvae catch mostly flies (Broadley & Stringer 2001), and the register's
glowing lures (G28) and lasting glow (M09) draw insects to the dusk fruit fall and its night food web (C07). Glows warn
and betray: glowing millipede models were attacked 22% of the time against 49% for dark ones (Marek et al. 2011); a
fish's strikes succeed over 75% of the time on prey first lit by a dinoflagellate's flash (Mensinger & Case 1992);
copepods reject dinoflagellates induced to glow brighter (Prevett et al. 2019); zooplankton that eat luminous bacteria
glow and are eaten eight times as fast (Zarubin et al. 2012, through Tanet et al. 2020). Off California 76% of the water
column's animals can glow (Martini & Haddock 2017), and the register's dark-night luminous plankton (G39) carry that
night in the seas' plankton community (C12). On the Open Moon these webs work where the night is dark enough to see a
glow: the dark belt's ground, forest floors, which keep 4–5×10⁻⁵ lux in the dark belt and 0.013–0.19 lux under the
Earth, and the seas off the far side's low latitudes (*screen*).

## Files

| File | Holds |
|---|---|
| [run.py](run.py) | The screen: natural light, living light and its cost, needs, skyglow and night life |
| [results/living_light.json](results/living_light.json) | Its product (schema `terluna.research.living-light/1`) |
| [test_run.py](test_run.py) | Bindings, exact regeneration, photon and energy budgets, limits and the skyglow model |
| [sources.json](sources.json) | Every source, its use and what was read |
| [checks.json](checks.json) | The checks run and the independent review's fixes |

The calendar's evaluator (`illumination/calendar`) is called live. Read in-tree with their Git blob and SHA-256
checked, now that the branches are merged: main's per-place Earthlight
(`research/studies/joint_synthesis/results/places.json`) and provisioning's daily-life and ways-of-living products.

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.living_light.run \
    --output research/studies/living_light/results/living_light.json
OPENBLAS_NUM_THREADS=1 python -m pytest research/studies/living_light
```
