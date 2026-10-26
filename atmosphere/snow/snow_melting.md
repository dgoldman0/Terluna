# How far below the freezing level lunar snow falls, and where its water and heat go

Below the freezing level of the design-climate storms, lunar snowflakes first sublimate, then melt, and then fall as
rain for 15–25 km, giving much of their water back to the air on the way. The air under falling snow is 93–97% humid
at the freezing level and warms about 1 K for each kilometre below it, so melting begins 0.3–0.55 km down, where a
flake's ice-bulb temperature reaches 0 °C. Flakes of 2 cm finish melting about 1.6 km below the freezing level, and
flakes of 10–20 cm 3.7–5.2 km below it, three to four times the 1.1–1.4 km the flake model of 20 October gives at 90%
humidity. The
drier tenth of the columns lowers the large flakes' melting by 0.4–0.7 km, and the crystal model's Earth error band
moves it by −7% to +9%. The large flakes melt to drops of 9–15 mm, which hold together under lunar gravity, where the
stability screen stands at 16.6 mm, and they bring 81–93% of their water to the sea or the ground; the 1–3 mm drops of
flakes up to 2 cm evaporate before they reach any ground. No snow reaches the ground: the largest flakes of the
populations finish melting at least 15 km above the sea, and the highest summit stands 11.6 km up.

For the snow the growth runs deliver, light snow of 0.3 mm/h melts within 2.3–2.4 km and sends 11–13% of its water to
the surface. Snow whose flakes reach about 10 cm, the heavy and moderate cases, melts within 3.4–3.8 km and sends
55–60%, and light snow forming near −30 °C sends 47–51%. Sublimation takes 2–8% of the snow over land and 10–23% over
the Near-side Sea before and while it melts. At 2.83 MJ/kg against 0.33 for melting, it draws about as much heat from
the air over the sea as the melting does in the heavier snow, and more than twice as much in light snow. Melting
draws 71–90 W/m² for each mm/h of snowfall in a layer 2–3.5 km deep below the freezing level, and the drops'
evaporation 230–550 W/m², spread over the rain's whole fall, with half of it within 5–9 km of the surface in the
heavier snow. Of its water, the snow returns 40–53% to the air on the way down in the larger-flake cases and 87–89%
in light snow, against the 61% CM1's own storms return of their falling water.

On Earth, with Earth's gravity and air, the same code begins melting within 3–5% of the depths Matsuo and Sasyo (1981)
computed and within 0.17 K of the temperatures Heymsfield and co-authors (2021) fitted to observations, and its
ordinary flakes finish where observed melting layers end. The lunar results are runnable calculations through archived
model air; no measurement checks them at lunar gravity or at these sizes. The flakes come from McSnow runs whose lunar
builds diffused vapour at 0.651 of McSnow's own fit, through a patch found on 25 October; how far that moves them is
measured separately.

## How a flake is followed

A flake leaves the freezing level with the mass McSnow gave it there, on McSnow's law for unrimed aggregates (after
Mitchell 1996: m = 0.0028 D^2.1 g and a projected area of 0.2285 D^1.88 cm² for a span D in cm), and falls at its
terminal speed by Böhm's law ([boehm.py](boehm.py)), which reproduces McSnow's own speeds for the arriving flakes
within 0.3%. It falls through still air, a column taken from the design-climate storm boxes, and passes through three
stages.

While it is dry ice, the flake's surface sits at the temperature where the heat the air conducts to it balances the
latent heat of the vapour it exchanges, its ice-bulb temperature. The air below the freezing level is warmer than 0 °C
but subsaturated, so the flake sublimates there, losing mass and cooling itself and the air, and its surface stays
below freezing. Melting begins only where the ice-bulb temperature reaches 0 °C (Heymsfield and co-authors 2021). That
level lies a little below the one where the wet-bulb temperature passes 0 °C, because the latent heat of sublimation
exceeds that of evaporation by an eighth.

Once melting begins the surface holds at 0 °C, and the model is Mitra and co-authors' (1990), built on their
wind-tunnel flakes, as Frick, Seifert and Wernli (2013) give it equation by equation. The heat the air conducts to the
flake, less the latent heat of the vapour it exchanges, melts ice. The flake is an oblate spheroid whose axis ratio
rises from 0.3 to 1 with its liquid fraction l, its density moves from the dry flake's to water's as l grows, and its
capacitance carries Mitra's factor 0.8 + 0.2 l. Its fall speed moves from the dry flake's to the drop's by Mitra's
measured function 0.246 l + 0.754 l⁷, which keeps the flake slow until about seven-tenths of it has melted. Its
ventilation is Hall and Pruppacher's (1976), at a Reynolds number taken between the dry flake's and the drop's in
proportion to l. Meltwater gathers first in a flake's hollows and then covers it (Leinonen and von Lerber 2018), so the
vapour the melting flake exchanges leaves exposed ice, as sublimation, in proportion to the ice's share of its mass,
and leaves the meltwater, as evaporation, in the rest. That split moves mass between the sublimated, melted and
evaporated records; the flake's ice and water, and the heat the air gives up, are the same for any split.

When the last ice has melted, the drop falls by the rain solver's physics ([drops.py](../rain/drops.py)): drag on its
equilibrium shape, its ventilated evaporation and its own temperature, to the ground or until it has evaporated. A drop
larger than its stability screen breaks into eight drops of half its diameter, again until they hold together. The
heights here are above each box's surface, the sea for the sea box and the box's own ground for the land box, whose
land the geography atlas puts 0.5–1.0 km above sea level.

Water vapour's diffusivity is Massman's (1998) measured fit, 0.2178 cm²/s at 0 °C and 1 atm, which the rain solver
takes since 25 October in place of its kinetic theory; that theory stood 17.6% lower, which weakened sublimation and
evaporation and started melting too high. The crystal model's Earth error band ([earth error band](earth_error_band.md))
enters as its multipliers on capacitance, 0.836–1.109, which scale the flake's heat and vapour exchange, and on terminal
speed, 0.555–1.159, which scale the dry flake's fall; the band's deep edge pairs the smaller capacitance with the larger
speed, and its shallow edge the reverse.

The sources, the quantities taken from each and how each was read are in [melting_sources.json](melting_sources.json).

## The air the snow falls through

The two boxes carry snow through their freezing level in 83,627 columns over the land box's second lunar day and 43,876
over the sea box's third, a third and a fifth of the columns they save. Weighted by that snow, the freezing level stands
25.0–26.8 km above the land box's ground (10th to 90th percentile; median 25.9 km) and 21.7–22.7 km above the sea
(median 22.3 km), 3–7 km below the 29.2 km of the climate targets' ladder that the calculation of 20 October took. Only
10–14% of these columns hold cloud at their freezing level, so most snow falls out of its cloud into clear air a few per
cent short of saturation:

| Depth below the freezing level | 0 | 1 km | 2 km | 4 km | 8 km | 12 km | 20 km |
|---|---|---|---|---|---|---|---|
| Land: temperature, relative humidity | 0 °C, 97% | +0.9 °C, 96% | +1.8 °C, 94% | +4.0 °C, 89% | +8.5 °C, 81% | +12.9 °C, 76% | +21.8 °C, 82% |
| Sea: the same | 0 °C, 94% | +1.0 °C, 92% | +2.1 °C, 91% | +4.4 °C, 86% | +9.2 °C, 82% | +13.8 °C, 87% | +24.4 °C, 70% |
| Land: humidity, 10th to 90th percentile | 90–100% | 89–100% | 87–99% | 81–97% | 73–91% | 68–88% | |
| Sea: the same | 88–100% | 87–99% | 85–98% | 81–93% | 75–90% | 72–95% | |

The first two rows are the snow-weighted medians the composite columns take. Most of the snow falls in the afternoon:
over land 60% between noon and sunset and 23% in the morning, over the sea 57% between noon and sunset, 37% before
midnight and almost none in the morning. The boxes' air already carries the cooling of CM1's own melting and
evaporation, at the rates of its own snow and rain.

## How far down each flake melts

In the composite columns melting begins 0.34 km below the freezing level over land, at +0.3 °C, and 0.55 km below it
over the sea, at +0.5 °C, for flakes of every size. Each flake then takes longer to melt the larger it is:

| Span at the freezing level | Its drop | Melted by, land | Melted by, sea | Water reaching the surface: land; sea |
|---|---:|---:|---:|---:|
| 4.7 mm | 1.0 mm | 0.79 km, +0.7 °C | 0.87 km, +0.9 °C | none; none |
| 9.2 mm | 1.7 mm | 1.10 km, +1.0 °C | 1.17 km, +1.2 °C | none; none |
| 1.8 cm | 2.6 mm | 1.54 km, +1.4 °C | 1.59 km, +1.6 °C | none; none |
| 4.4 cm | 4.9 mm | 2.46 km, +2.3 °C | 2.45 km, +2.6 °C | 42%; 46% |
| 8.5 cm | 7.8 mm | 3.42 km, +3.4 °C | 3.34 km, +3.6 °C | 74%; 77% |
| 10.6 cm | 9.1 mm | 3.81 km, +3.8 °C | 3.71 km, +4.1 °C | 81%; 83% |
| 13.2 cm | 10.6 mm | 4.23 km, +4.2 °C | 4.10 km, +4.5 °C | 85%; 87% |
| 16.5 cm | 12.4 mm | 4.68 km, +4.8 °C | 4.53 km, +5.1 °C | 89%; 91% |
| 20.6 cm | 14.5 mm | 5.17 km, +5.3 °C | 4.98 km, +5.6 °C | 92%; 93% |
| 25.7 cm | 17.0 mm, broken into 8.5 mm drops | 5.70 km, +5.9 °C | 5.48 km, +6.2 °C | 79%; 82% |

A flake's mass grows as its span to the power 2.1 while its fall speed barely changes, 0.6–0.7 m/s for every span
above a few centimetres, and its capacitance and ventilation grow more slowly than its mass. Each doubling of span
therefore deepens its melting, by about 0.4 km near a centimetre and by 0.9–1.4 km above 4 cm, and the largest flakes
melt 3–4 times deeper than flakes of 2 cm.

Humidity moves where melting begins, and the large flakes' depths with it. In the moister tenth of the columns, near
saturation at the freezing level, melting begins there, and the 10.6 cm flake finishes 3.1–3.2 km down; in the drier
tenth, at 88–90%, it begins 1.1 km down, and the same flake finishes 4.1–4.5 km down and the 20.6 cm flake 5.4–5.9 km
down. Drawn by their snow, 120 actual columns from each box put the melting of 11 cm flakes 3.3–4.3 km below the
freezing level (10th to 90th percentile) and of 19 cm flakes 4.3–5.4 km below it, with melting beginning anywhere from
the freezing level itself to 0.9–1.1 km below it. By time of day the air is moistest after sunset over land, where
melting begins 0.07 km down, and driest after midnight in both boxes, where it begins 0.67–0.74 km down; through all
four quarters the 10.6 cm flake finishes between 3.4 and 4.0 km down. The crystal model's error band moves the same
flake's depth to 3.5–3.6 km at its shallow edge and 4.0–4.2 km at its deep edge, −7% to +9%.

The drops of flakes up to about 2 cm evaporate on the way down: in the sampled columns those of 5 mm flakes vanish
1.3–11 km below the freezing level and those of 2 cm flakes 12–24 km below it, so that 2 cm flakes bring a mean 1–3% of
their water to the surface. Meltwater drops stay below the lunar stability screen, 16.6 mm at the melting level ([rain
breakup](../rain/README.md#breakup-at-lunar-gravity)), up to flakes of about 24 cm. Earth's screen, 6.6 mm, would break
the drop of every flake above about 7 cm. The largest flakes of the populations, 30–34 cm, melt to drops of 19–21 mm,
which break.

## The snow of the growth runs: its bands, its water and its heat

The flakes McSnow delivers to the freezing level ([growth and aggregation](mcsnow_growth_aggregation.md)), weighted by
their mass flux, melt and evaporate in the composite columns as follows (land; sea):

| Snow | Half the ice gone | 99% of the ice gone | Sublimated | Evaporated | Reaching the surface |
|---|---:|---:|---:|---:|---:|
| Light, 0.3 mm/h | 0.91–0.96; 0.96–1.01 km | 2.33–2.40; 2.32–2.38 km | 7–8%; 21–23% | 80–81%; 66% | 11–13%; 11–13% |
| Light, a tenth of the crystals | 0.94; 1.00 km | 2.46; 2.43 km | 7%; 21% | 81%; 67% | 12%; 12% |
| Light, formed near −30 °C | 1.34–1.40; 1.37–1.43 km | 3.57–4.12; 3.47–3.98 km | 3–4%; 12–13% | 45–49%; 38–40% | 48–51%; 47–50% |
| Moderate, 1 mm/h | 1.51; 1.53 km | 3.84; 3.72 km | 2%; 10% | 38%; 32% | 59%; 58% |
| Heavy, 3 mm/h, a hundredth of the crystals | 1.48–1.49; 1.50–1.51 km | 3.60–3.75; 3.50–3.64 km | 2%; 10% | 40%; 34% | 57–58%; 56% |
| Held at saturation, a hundredth of the crystals | 1.47; 1.49 km | 3.52; 3.43 km | 2%; 10% | 41%; 35% | 56%; 55% |

Ranges give the two superparticle numbers of each converged case; depths are below the freezing level. Sublimation is
the ice lost to vapour before and while the flakes melt, and evaporation the meltwater lost from melting flakes and
drops. These bands are where the sky's winter weather lies. Over the sea, with the freezing level 22.3 km up, dry snow
falls above 21.75 km; wet snow and sleet fill the band below it, down to 19.9 km in light snow and 18.3–18.8 km in the
heavier cases; rain falls beneath. Over land the same bands stand about 3.6 km higher above the box's ground, which
itself stands 0.5–1.0 km above the sea.

The snow's latent heat, in W/m² for each mm/h of snowfall arriving at the freezing level:

| Snow | Dry layer: sublimation | Melting layer: melting; sublimation; evaporation | Rain: evaporation | All |
|---|---:|---:|---:|---:|
| Light, land | 40–45 | 85; 13–15; 0–3 | 551–554 | 690 |
| Light, sea | 133–147 | 71–73; 30–33; 11–16 | 436–442 | 690 |
| Light from −30 °C, land | 24–26 | 89; 1–2; 0–2 | 308–333 | 430 |
| Light from −30 °C, sea | 79–85 | 80–81; 16–17; 6–7 | 251–269 | 430 |
| Heavy, a hundredth of the crystals, land | 20 | 90; −1; −11 | 283–285 | 380 |
| Heavy, a hundredth of the crystals, sea | 66 | 82; 13; −3 to −4 | 234–235 | 390 |

The dry layer reaches from the freezing level to where melting begins, the melting layer from there to where 99% of
the ice has gone, and the rain from there to the surface. Negative values are vapour condensing on melting flakes where
the air holds more vapour than saturation at the flakes' 0 °C surface. Over the sea the sublimation, at 2.83 MJ/kg,
takes about as much heat as the melting in heavy snow and more than twice as much in light snow. The water goes back
to the air at heights set by the drops' size: light snow over the sea releases a quarter of its water as vapour within
2 km of the freezing level, by sublimation, and half by 7–8 km down; heavy snow releases a quarter by 5 km down, half
by 14 km down and three-quarters by 19 km down, 3 km above the sea.

## Rain at each ground height

No snow reaches any ground in these storms. Rain does, and the higher the ground the less of the fall its drops must
survive. Over the sea box's column, by height above the sea (free air above the sea standing in for the air above
higher ground):

| Ground height above the sea | 0 | 1 km | 2 km | 3 km | 5 km | 8 km | 11.6 km |
|---|---:|---:|---:|---:|---:|---:|---:|
| Light, 0.3 mm/h | 13% | 16% | 20% | 23% | 30% | 39% | 48% |
| Light, formed near −30 °C | 47–50% | 51–54% | 54–58% | 58–61% | 64–67% | 70–73% | 75–77% |
| Heavy, a hundredth of the crystals | 56% | 60% | 64% | 67% | 73% | 78% | 82% |

Over the land box's column, by height above its own ground, light snow delivers 11–13% there and 36–42% at 11.6 km,
and heavy snow 57–58% and 82–83%. The geography atlas puts half the land below 1.9 km above the sea and nine-tenths
below 5.1 km, with the highest summit at 11.6 km, so most ground receives 15–30% of light snow's water and 60–73% of
heavy snow's.

## Against the flakes of 20 October and CM1's own snow

In the straight column of 20 October (the climate targets' ladder, 0.92 K/km from a freezing level at 29.2 km, at 90%
humidity), the flake model of that day, which let no ice sublimate before it melted, finished 2 mm to 2 cm flakes
1.0–1.3 km below the freezing level; with the rain solver's corrected diffusivity it finishes them 1.1–1.4 km down, and
3.1–3.4 km down at 70%. This model, in the same column at 90%, begins melting 0.9 km down, after the slow flakes have
spent half an hour or more sublimating: the 2 mm flake vanishes before it can melt, the 5 mm flake melts by 1.0 km with
97% of its mass already sublimated, and the 2 cm flake melts by 1.9 km after losing 44%. Flakes of 4.5 cm finish 2.8 km
down, and of 8–19 cm 3.6–5.3 km down, losing 10–18%. At 70% humidity every flake up to 4.5 cm sublimates away before
melting, and 11–19 cm flakes melt 4.8–6.1 km down after losing 75–89% of their mass. Sublimation before melting
therefore matters far more on the Moon than on Earth. On Earth a flake crosses the layer above its melting in about 90 s
at 90% humidity; a lunar flake takes 10–15 minutes in the storm boxes' air and half an hour in that column.

CM1's own snow and graupel, in the same boxes, melt with half their cooling within 1.6–1.8 km of the freezing level and
nine-tenths within 3.3–3.9 km: close to the heavier snow's 90% (2.6–2.9 km) and 99% (3.5–4.1 km) here, and deeper than
light snow's 1.7 km. CM1's rain evaporates with half its cooling within 6.6–8.0 km of the freezing level and nine-tenths
within 18.6–18.7 km. Its falling water, rain from its clouds and melted snow and graupel, loses 60.4–60.9% on the way
down ([climate/crm](../../climate/crm/README.md)), where the heavier snow here loses 40–45% and light snow 87–89%. The
difference is the drops' size: CM1's rain falls in exponential spectra whose breakup limit, carried to lunar gravity,
settles near a mean diameter of 1.7 mm ([rain breakup](../rain/README.md#breakup-at-lunar-gravity)), while the large
flakes melt to drops of 6–15 mm that keep most of their water through a fall of 20 km.

## Above the freezing level: the storm paths

The growth-and-aggregation runs already carry the vapour exchange above the freezing level
([growth and aggregation](mcsnow_growth_aggregation.md)). The reference columns are held at water saturation, so their
snow sublimates only below the freezing level. Along the archived storm's 354 paths the air is often subsaturated over
ice, and each path column's ice can be followed every 30 minutes: its change plus the snow it delivers at 0 °C is its
net growth, and the sums of the positive and negative parts bound from below the ice it grew and the ice it lost to
sublimation, since growth and sublimation in one column within one interval cancel.

| Storm paths, crystal concentration | Paths delivering snow | Ice grown that sublimates before 0 °C, at least | Arrives at 0 °C | Leaves as single crystals at the 1 m censor | Aloft at the end |
|---|---:|---:|---:|---:|---:|
| Archived humidity, a hundredth | 28–34 | 12–14% | 6–7% | 79–81% | 0.1% |
| Archived humidity, a tenth | 63–71 | 14–16% | 37–42% | 44–46% | 0.4–0.6% |
| Archived humidity, reference | 77 | 14% | 82% | 2% | 1.0% |
| Held at water saturation, a hundredth to reference | 276–300 | none | 47–88% | 6–46% | 6–12% |

Ranges give the base runs and those with four times the superparticles. In the archived humidity at least an eighth of
the ice the storm grows sublimates before it reaches the freezing level, and of the ice that is neither censored nor
still aloft, 15% sublimates at the reference concentration, 24–30% at a tenth and 65–69% at a hundredth. Where
crystals are few, most of the ice that does not run away into metre-long single crystals sublimates before it can
reach the freezing level. Each kilogram that sublimates there takes back from the air the 2.83 MJ its growth
released.

## The Earth check

With Earth's gravity and air, the same code melts McSnow's aggregates in a straight winter sounding with its freezing
level at 1.5 km and 850 hPa. Where melting begins tests the flake's heat and vapour exchange before it melts; where it
ends tests the melting model with them.

| Relative humidity below the freezing level | 90% | 80% | 70% | 50% |
|---|---:|---:|---:|---:|
| Melting begins, metres below the freezing level at 6 K/km: this model | 125 | 258 | 402 | 735 |
| The same, computed by Matsuo and Sasyo (1981) | 120 | 250 | 390 | 700 |
| Air temperature there: this model | 0.75 °C | 1.55 °C | 2.41 °C | 4.41 °C |
| Where the ice-bulb temperature is 0 °C, from Heymsfield and co-authors' (2021) fits at 850 hPa | 0.86 °C | 1.72 °C | 2.57 °C | 4.29 °C |

The model begins melting within 3–5% of Matsuo and Sasyo's depths and within 0.17 K of the temperatures Heymsfield and
co-authors fitted to aircraft and surface observations.

At the standard 6.5 K/km, flakes of 2 mm to 2 cm, which melt to drops of 0.6–2.8 mm, finish melting 140–580 m below the
freezing level in saturated air, at +0.9 to +3.8 °C, and 190–640 m below it at 95% humidity, at +1.2 to +4.1 °C. The
centimetre flake, a 1.75 mm drop, finishes at +2.5 to +3.2 °C between saturation and 90%. Flakes with 2.2 times the
mass at each span, Heymsfield and co-authors' mass law, finish 0.2–0.5 K warmer. Observed melting layers are 300–350 m
deep at the median over four radar sites (Wolfensberger and co-authors 2016, as Heymsfield and co-authors cite them),
melting completes 1–2 K above the level where it begins and in moist aircraft spirals by +3 °C, and the largest
aggregates survive to about +5 °C (Willis and Heymsfield 1989, as cited). The model's 5 mm flakes, 1.1 mm drops,
finish 250–350 m down and 1.6–2.0 K above where they began, inside those observations. Its centimetre flakes finish
390–490 m down and 2.4–2.5 K above, slightly past the observed spread, and its 2 cm flakes at +3.8 to +4.5 °C, toward
the warmest observed aggregates. In air at 70% humidity its smallest flakes sublimate away before they melt, as Matsuo
and Sasyo found for small light flakes in dry air.

## Limits

Each flake falls alone through air that its melting and evaporation do not change. A real storm's snow cools and
moistens the layers it melts and evaporates in, which raises the humidity, lowers the freezing level and shortens the
dry layer; the boxes' air carries CM1's own melting and evaporation of its own snow, at its own rates, and no more.
The air is still: the boxes' downdraughts, which lengthen a flake's time in each layer of descending air, and their
updraughts are left out. Below the freezing level flakes do not aggregate, although wet flakes stick well, and they
melt to one drop each; drops fall alone, so collisions between drops, and their coalescence and breakup, which the rain
work's encounter diagnostic shows would reshape the spectrum ([rain](../rain/README.md)), are left out. A drop above
its screen breaks into equal drops, a stand-in for a breakup spectrum.

The flake populations at the freezing level are those of the growth-and-aggregation runs of 25 October. Their lunar
McSnow builds took the freezing level's pressure, 659 hPa, as the reference of McSnow's vapour diffusivity, which is 1
atm, so they diffused vapour at 0.651 of McSnow's fit; where diffusion limits growth, their crystals grew at 0.68–0.79
of its rate. The flakes' sizes, and with them the depths and shares here, carry whatever that changes until the
corrected runs replace them.

Mitra's model carries its laboratory flakes to the lunar sizes. Those flakes were centimetre-sized; lunar flakes of
10–20 cm fall at Reynolds numbers of several thousand, beyond the range of the ventilation law, and their melting
geometry and fall speed are the laboratory forms extrapolated. The flakes' masses come from McSnow's single law for
unrimed aggregates, so rimed or denser flakes, which melt deeper, are covered only by the Earth check's denser variant.
On Earth the model reproduces where melting begins and how warm the air is where ordinary flakes finish melting, which
tests its heat and vapour exchange and its onset; no measurement checks it at lunar gravity, in lunar air or at these
sizes.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting_environment   # the boxes' snowy columns, a minute
OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting earth          # seconds
OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting ladder         # seconds
OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting paths          # a minute
OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting moon           # several minutes
OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.melting assemble       # results/snow_melting.json
```

The environment step reads the two boxes' saved fields through `climate/crm/runs` on the project-space drive and keeps
the snowy columns, with their hashes, in `research/runs/snow_agent_3/environment`; the other steps write their full
records to `research/runs/snow_agent_3/melting`, from which the compact [product](results/snow_melting.json) is
assembled. The McSnow populations and storm paths are read from `research/runs/snow_agent_2/growagg`
([growth and aggregation](mcsnow_growth_aggregation.md)), each run recorded by the hash of its `run.json`.
