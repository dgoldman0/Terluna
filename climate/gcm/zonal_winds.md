# The days of a drifting town

The author has adopted roaming sky towns: "Allow much floating habitation to roam, with collision and storm risks
managed" ([register](../../research/decisions.md#life-and-people)). A town that drifts with the wind lives a solar day
set by the wind it rides, and the same holds for floating life. This note reads that day from the design climate's
winds, by height and latitude. The numbers come from [zonal_winds.py](zonal_winds.py) and its product
[zonal_winds_A28_dim5_moon.json](../results/gcm/zonal_winds_A28_dim5_moon.json), built from years 20–29 of the design
run `A28_dim5_moon`.

The subsolar point moves west over the Moon at v☉ = 2π(R + z) cos φ / 29.53 days: 4.30 m/s on the equator 10 km up,
2.15 m/s at 60° and 0.75 m/s at 80°. A body drifting east at u passes through the solar hours at u + v☉ and lives a
day of 29.53 v☉ / |u + v☉| Earth days. At rest it lives the 29.53-day lunar day; an eastward wind shortens the day; a
westward wind slower than v☉ lengthens it; u = −v☉ holds the hour; a faster westward wind turns the Sun back, so that
it rises in the west. Hour angles below are degrees east of the subsolar point: 0° is noon, −90° sunrise and +90°
sunset, and on a 24-hour clock of the lunar day +60° reads 16:00.

## Findings

**Above about 8 km a drifting town lives a day shorter than the lunar day, and the higher it floats the shorter the
day.** The dense, tall air superrotates. Its mean eastward wind is 0.85 m/s at 7.5 km, 3.5 m/s at 20 km, 10.7 m/s at
40 km and 21.9 m/s at 75 km, and from 20 km up it blows eastward in at least 98% of 3-day means. On the equator the
median drifting day is 21.8 days at 10 km, 17.1 at 20 km, 12.7 at 30 km, 9.3 at 40 km and 4.8 at 75 km. At 60° it is
17.3, 10.4, 6.5, 4.3 and 2.5 days: the wind keeps its speed toward the poles while the Sun's own motion slows with the
cosine of latitude. In the 30–45 km band, around the 35–45 km flight band, the day runs 7.9–12.6 days within 20° of
the equator (10th to 90th percentile) and 4.5–8.2 days at 40–60°.

**Holding the hour belongs to the lowest 5 km poleward of about 40°, from the late afternoon into the night.** There
the return flow toward the afternoon runs westward at the Sun's own speed: at 2.5 km and 64° the mean Sun-relative
speed is 0.0–0.2 m/s from +65° to +105°. Pooled over the lowest 5 km, 3-day means hold the hour within 1 m/s 18% of
the time at 40–60° and 31% at 60–75°, and 49% and 66% at the best hour (+85° and +65°). At 2.5 km the hour is held in
more than 30% of 3-day means at 25–75° from +15° through the evening and night to −135°, and in 85% at 69° N and +55°.
Within 20° of the equator the lowest 5 km hold it 5% of the time, 18% just after sunset. From 7.5 km up only the polar
caps and 69–75° hold it that often, and from 20 km up nowhere equatorward of 76° does.

**The Sun turns back for a drifting town only near the ground at high latitude, mostly in the evening.** In the lowest
5 km the local wind turns the Sun back in 9% of 3-day means at 40–60°, 17% at 60–75° and 38% poleward of 75°, from
+45° to midnight. Along a route this comes as an episode around sunset: the Sun sets, rises again in the west and sets
once more. In its first year 6% of routes at 2.5 km and 43% at 5 km see it at least once, and at 5 km 5–12% of sunsets
are followed by it. Every route at every height passed through all hours the Sun's ordinary way within its first year;
none lived a whole day backward.

**Near the ground a drifting town lives a long evening and a short morning.** On the equator at 2.5 km the air parts
15° before dawn and meets 64° past noon, and the Sun-relative speed falls from 6.9 m/s at noon to 2.2–2.3 m/s from
+95° to +125°. A town at 2.5 km, which the inflow soon brings near the equator, takes 37 days to pass once through all
hours (median; 30–43 days from the 10th to the 90th percentile), and over ten years its day averages 34 days: 12 for
the evening and first half of the night, 10 for the afternoon, 7 after midnight and 5 for the morning, against 7.4
each at rest. The Sun is up 44% of its time, and it gathers 0.80 of the sunlight fixed places at its latitudes get.

**The hour a town lingers in moves from the evening to the morning as it rises.** In its first year a town lingers
most between +60° and +120° at 2.5 km, +30° and +90° at 5 km, noon and +60° at 7.5 km, around noon at 10 km, and −90°
and −30° at 15–20 km; from 30 km up it passes all hours at an even pace. At 7.5 km the Sun is up 55% of its time and
it gathers 1.13 times the light of fixed places at its latitudes; at 10–20 km, 1.05–1.14 times. A town chooses which
part of its day stretches by choosing its height.

**The steady flow of the Sun's frame stops only at high latitude, and the 3-day weather moves towns on within weeks.**
The mean Sun-frame flow stops near the ground in the evening at 65–73° (+77° to +117°; weak spiral sinks in the north
and spiral sources in the south) and near the poles at 7.5–20 km (77–85°), mostly in the morning (−34° to −77°; one at
+63°); above 20 km it never stops. Its sinks gather parcels over e-folding times of 21–159 days. In the 3-day means no
route stays within a 30° span of hours longer than 45 days. The median route's longest such stay in a year is 10 days
at 2.5–5 km, 5 days at 7.5–20 km and under a day from 40 km up away from the poles, against 2.5 days at rest.

**Drifting towns gather onto shared routes.** At a fixed height the horizontal flow converges where the air rises
(near the ground) or sinks (aloft), and a body that holds its height stays behind as the air leaves. Routes released
300 km apart at 2.5 km come within 12 km of a neighbour (median) after 30 days, and after 180 days every route at that
height runs within one 171 km model cell of another, on a single track at the model's resolution; at 5 and 7.5 km
98–99% do within a year and all within three. Near the ground the equatorward inflow, 0.5–0.9 m/s between 10° and 70°,
brings 98% of the routes at 2.5 km within 20° of the equator in 60 days, and in years 5–10 they spend 60% of their
time there, where 34% of the Moon's area lies. Aloft a poleward drift of 0.03–0.08 m/s carries routes at 40–60 km into
the polar caps within a few years: in years 5–10 they spend 97–100% of their time poleward of 70°, where 6% of the
area lies, and at 40 km the median distance to the nearest neighbour falls from 310 km to 80 km in a year and 1.4 km
in five. Routes spread most evenly at 30 km, where the meridional drift vanishes; after ten years their nearest
neighbour is still 117 km away (median). Collision risk for roaming towns concentrates wherever the wind gathers
routes: within weeks near the ground and within years at 40–60 km.

**Near the ground the wind a town meets is set mostly by the hour and the place; aloft it is set by the steady
eastward flow.** At 2.5 km the mean by hour at each place explains 85% of the eastward wind's variation about its
zonal mean, and the mean by hour in the Sun's frame alone 50%. From 10 km up these explain 19–39% and 7–14%, and the
rest is weather, but the variation there is 1.0–1.6 m/s rms about mean winds of 1.4–22 m/s.

**The shear to sail on is about 4 m/s over 15 km in the lowest 20 km and 6 m/s in and below the flight band.** The
median difference between the winds 15 km apart is 3.9–4.9 m/s for a lower end between the ground and 15 km, 5.3 m/s
from 20 km and 6.3–6.4 m/s from 30–40 km, with 90th percentiles of 6.1–8.4 m/s; below 55 km it is 1.6–2.7 m/s over 5
km and 0.8–1.8 m/s over 2.5 km. From 20 km up nearly all of it is the eastward wind's growth with height, 5.1–6.2 m/s
over 15 km on average, and it changes by about 10% through the solar day. Nearer the ground the shear turns with the
hour: from 2.5 km it is 2.1–2.6 m/s in the late morning and 5.4–7.0 m/s in the late afternoon (+55° to +95°), the
hours when the wind there also holds the hour. StratoSail hangs a wing on a tether up to 15 km below a stratospheric
balloon and biases the balloon's drift by a few m/s, drawing on relative winds of 11–26 m/s between 20 and 35 km in
the cases its designers cite, with air about ten times denser at the wing than at the balloon (Aaron, Heun & Nock
2002, Advances in Space Research 30, 1227). The Moon's shear in the flight band is a quarter to a half of those, and
its air is 1.3 times denser 15 km down; the air itself is 0.73–0.94 kg/m³ between 25 and 40 km ([global
winds](../results/gcm/global_winds_A28_dim5_moon.json)), as dense as Earth's air 3–5 km up and ten times StratoSail's
0.088 kg/m³ at 20 km. A wing hung 15 km below a town at 40 km meets 19–33 Pa of dynamic pressure at the median and
90th-percentile shear, against StratoSail's 5–30 Pa, so a wing of a given size draws forces as large as StratoSail's
or larger (*derived*). For a town at 40 km on the equator the mean wind 15 km lower is
4.2 m/s against 9.4 m/s at the town: drifting at either speed it would live 15.0 or 9.4 days, and a sail between them
chooses within that range.

## By height and latitude

Days are those of a body drifting with each 3-day mean wind (10th, 50th and 90th percentiles, Earth days); "longer
than 29.5 d" is the share of those days longer than the lunar day, the Sun's direction either way. "Held" is
|u + v☉| under 1 m/s, given for the whole band and for its best hour. Shear is the difference between the winds 15 km
apart, with its lower end in the band. Poleward of 75° the Sun moves at under 0.75 m/s, so a 1 m/s margin there admits
most light winds; the days describe those caps better than the held share does.

| Height | Latitude | Day: 10th, 50th, 90th percentile | Longer than 29.5 d | Held within 1 m/s; at its best hour | Sun turned back | Shear over 15 km: median, 90th |
|---|---|---|---|---|---|---|
| 0–5 km | 0–20° | 17.1, 27.9, 75.1 d | 44% | 5%; 18% at +95° | 1% | 3.7, 6.5 m/s |
| 0–5 km | 20–40° | 17.1, 28.8, 92.2 d | 47% | 8%; 31% at +75° | 3% | 4.2, 7.1 m/s |
| 0–5 km | 40–60° | 15.6, 29.6, 140.7 d | 49% | 18%; 49% at +85° | 9% | 4.7, 7.7 m/s |
| 0–5 km | 60–75° | 12.3, 26.7, 146.7 d | 45% | 31%; 66% at +65° | 17% | 4.3, 7.2 m/s |
| 0–5 km | 75–90° | 5, 13.9, 74.9 d | 24% | 41%; 60% at +25° | 38% | 3.9, 6.9 m/s |
| 5–15 km | 0–20° | 16, 21.6, 38 d | 19% | 1% | 0 | 3.5, 5.8 m/s |
| 5–15 km | 20–40° | 15, 20.9, 39.3 d | 19% | 2%; 7% at +65° | 0% | 4.3, 6.7 m/s |
| 5–15 km | 40–60° | 12.5, 18.8, 42.3 d | 18% | 4%; 13% at +65° | 2% | 4.8, 7.2 m/s |
| 5–15 km | 60–75° | 9, 16.5, 62.1 d | 23% | 15%; 20% at +55° | 7% | 4.1, 6.6 m/s |
| 5–15 km | 75–90° | 4.4, 11.6, 62.7 d | 21% | 35%; 45% at +5° | 27% | 2.9, 5.3 m/s |
| 15–30 km | 0–20° | 11.9, 15.3, 21.9 d | 2% | 0% | 0 | 5.0, 7.2 m/s |
| 15–30 km | 20–40° | 10, 13.2, 19.2 d | 1% | 0% | 0 | 5.6, 7.7 m/s |
| 15–30 km | 40–60° | 7.5, 10.5, 16.2 d | 0% | 0% | 0 | 6.3, 8.4 m/s |
| 15–30 km | 60–75° | 5.6, 8.4, 16.5 d | 3% | 2%; 6% at −65° | 1% | 6.0, 8.6 m/s |
| 15–30 km | 75–90° | 3.5, 7.4, 35 d | 12% | 22%; 33% at −65° | 15% | 3.5, 6.4 m/s |
| 30–45 km | 0–20° | 7.9, 9.7, 12.6 d | 0% | 0% | 0 | 5.9, 7.5 m/s |
| 30–45 km | 20–40° | 6.4, 8.1, 10.5 d | 0% | 0% | 0 | 6.2, 7.7 m/s |
| 30–45 km | 40–60° | 4.5, 6, 8.2 d | 0% | 0% | 0 | 6.7, 8.4 m/s |
| 30–45 km | 60–75° | 3, 4.2, 6.3 d | 0% | 0% | 0 | 7.1, 9.1 m/s |
| 30–45 km | 75–90° | 2.1, 3.7, 10.3 d | 3% | 7%; 12% at −75° | 4% | 4.8, 7.6 m/s |
| 45–60 km | 0–20° | 5.8, 6.7, 8 d | 0% | 0% | 0 | 5.5, 7.1 m/s |
| 45–60 km | 20–40° | 4.8, 5.6, 6.7 d | 0% | 0% | 0 | 4.2, 6.2 m/s |
| 45–60 km | 40–60° | 3.2, 4.1, 5.1 d | 0% | 0% | 0 | 3.6, 6.0 m/s |
| 45–60 km | 60–75° | 2, 2.6, 3.4 d | 0% | 0% | 0 | 4.6, 6.9 m/s |
| 45–60 km | 75–90° | 1.3, 2, 3.6 d | 0% | 2% | 1% | 3.7, 6.3 m/s |
| 60–90 km | 0–20° | 4.8, 5.3, 5.9 d | 0% | 0% | 0 | 3.1, 5.9 m/s |
| 60–90 km | 20–40° | 4.4, 4.8, 5.4 d | 0% | 0% | 0 | 2.1, 4.3 m/s |
| 60–90 km | 40–60° | 2.9, 3.6, 4.3 d | 0% | 0% | 0 | 1.8, 3.5 m/s |
| 60–90 km | 60–75° | 1.6, 2.2, 2.6 d | 0% | 0% | 0 | 2.0, 4.0 m/s |
| 60–90 km | 75–90° | 1, 1.4, 2.3 d | 0% | 1% | 1% | 2.2, 4.2 m/s |

The routes, by height. First passes and the mean day come from the ten-year route set; the longest stays and the
western sunrises from the ten one-year releases. Light is the route's mean sunlight over the mean that fixed places at
the latitudes it passed would get.

| Height | First pass through all hours: median (10th–90th) | Mean day over ten years | Longest stay within 30° of hour in a year: median, longest | Sun up | Light against fixed places | Routes with a western sunrise in the first year | Nearest neighbour after 1 / 10 years | Where routes spend years 5–10 |
|---|---|---|---|---|---|---|---|---|
| 2.5 km | 37 d (30–43) | 33.9 d | 10, 20 d | 44% | 0.80 | 6% | 0 / 0 km | 60% within 20° of the equator |
| 5 km | 28 d (22–41) | 30.9 d | 10, 45 d | 48% | 0.82 | 43% | 0 / 0 km | one track in the north, 83% of its time at 20–70° |
| 7.5 km | 22 d (18–32) | 29.4 d | 5, 20 d | 55% | 1.13 | 7% | 0 / 0 km | 79% within 20° of the equator |
| 10 km | 21 d (14–26) | 21.3 d | 5, 20 d | 54% | 1.09 | 10% | 32 / 52 km | 46% within 20° of the equator, 11% poleward of 70° |
| 15 km | 16 d (11–21) | 13.1 d | 5, 15 d | 55% | 1.14 | 9% | 94 / 35 km | 81% poleward of 50° |
| 20 km | 13 d (9–18) | 10.9 d | 5, 15 d | 52% | 1.05 | 28% | 20 / 27 km | 94% poleward of 50° |
| 30 km | 9 d (5–13) | 10.2 d | 1, 7 d | 50% | 0.99 | 3% | 157 / 117 km | 40% within 20° of the equator, 6% poleward of 70° |
| 40 km | 6 d (3–9) | 3.9 d | 0, 7 d | 51% | 1.04 | 1% | 80 / 7 km | 97% poleward of 70° |
| 50 km | 5 d (2–7) | 3.1 d | 0, 7 d | 51% | 1.03 | 1% | 46 / 0 km | all poleward of 70° |
| 60 km | 4 d (2–6) | 2.6 d | 0, 10 d | 50% | 1.03 | 3% | 45 / 0 km | all poleward of 70° |

The western sunrises at 20 km and from 40 km up come from routes near the poles, where the Sun moves at under 0.75 m/s
and stays within a few degrees of the horizon.

## For floating life

**A floater can hold the noon Sun on the wind only at high latitude, where that Sun stands low.** From 2.5 km up,
within 30° of noon, the wind holds the hour within 1 m/s in more than 20% of 3-day means only poleward of 52°, and in
more than 40% only poleward of about 69°, where the noon Sun stands at most about 21° high. Elsewhere a floater
drifting with the wind passes noon at the pace of the table above.

**Drifting gains a floater most light at 7.5–15 km and costs it light in the lowest 5 km.** At 7.5–15 km a floater has
the Sun up 54–55% of its time and 1.09–1.14 times the light of fixed places at its latitudes, because it lingers in
daylight: the afternoon at 7.5 km, noon at 10 km and the morning at 15 km. In the lowest 5 km it has the Sun up 44–48%
of its time and 0.80–0.82 times the light, because it lingers through the evening. A floater that changes height
between these layers chooses between them; one at 30 km and above meets the light of a fixed place, 0.99–1.04 times.

## Reading the product

[zonal_winds_A28_dim5_moon.json](../results/gcm/zonal_winds_A28_dim5_moon.json) (`terluna.climate.gcm-zonal-winds/1`)
holds the statistics by height (0–90 km above the model's sea level, every 2.5 km), latitude (the model's 32 rows) and
hour angle (36 bins of 10°). The eastward wind by height and latitude is `ground_frame.u_mean_m_s[height][lat]`, with
percentiles beside it; the Sun-frame means are `sun_frame`, the held and reversed shares `holding`, the drifting day's
percentiles and classes `solar_day`, the shear `shear`, the band table above `bands`, and the routes, their gathering
and the steady flow's stopping points `trajectories`. Its reading rule names every array. The `.npz` beside it holds
the Sun-frame percentiles and spread, the shear by solar hour and the routes' positions.

## What the numbers rest on

The winds are ExoPlaSim's 3-day means at T21, about 170 km cells: storms, gusts and turbulence are averaged out, and
each mean spans 36° of the Sun's path, so the days and holds are those of the large-scale flow. The model's ground is
its cells' mean elevation; at 2.5 km routes spend about 40% of their time over ground higher than their height, where
they ride the lowest model layer, about 0.9 km above the ground. The ten layers stand near 0.9, 4, 10, 18, 29, 42, 56,
74 and 98 km above the ground, so a shear over less than their spacing is the linear share of the shear between two
layers. Routes are horizontal at fixed heights and carried by the 3-day mean wind alone, as for a body that holds its
height and does not steer; their joining onto one track holds at the model's resolution, within which the model has no
structure to keep them apart.
