# The power system: beamed power from the ring fleet

On 12 October 2026 the author chose beamed power from the ring fleet for food and the sky towns, on the design in the
[domain synthesis](../research/studies/domain_synthesis/README.md#power-from-the-array-the-chosen-distribution)
([register](../research/decisions.md#the-domain-synthesiss-decisions-12-october-2026)):
- 2.45 GHz beams straight to the sky towns' crowns;
- closed receivers at the ground's loads, with planting beneath the mesh, joined by a grid;
- several beams for every receiver from different parts of the ring, moved before a sky reef crosses a column;
- retrodirective beams, with exposure above 10 W/m² kept to closed ground, decks and columns.

This note screens that design for provisioning question 3 ([questions](questions.md)). The screen,
[power.py](power.py), writes [results/power.json](results/power.json) in about ten seconds from the synthesis's
[ledgers](../research/studies/domain_synthesis/results/ledgers.json), the shield study's ring layout, integrated
ledger, array heat and links, the distributions study's C′ with the food blend, the food blend, the design winds' sky
towns, the cloud-resolving rings' storms, the ecology work's sky bands, the inhabited volume, sky ecology, the summit
tower and the atmosphere domain's moon_1.2atm column. Numbers marked *screen* come from that product; *estimate* marks a
result resting on a stated approximation, and *design guess* a value chosen for the design. Sources, with how far each
was read, are in [power_sources.json](power_sources.json).

Pairs are the low-cost and high-cost ends of C′ with the blend, as the synthesis gives them. The intake is the beam
power arriving at receivers and crowns, 30–145 TW with close computing, which the synthesis counts as heat on the Moon;
the delivered power is what reaches the loads, 22–105 TW.

## Findings

1. **The deep air takes 5.9% of a beam overhead and 6.3–7.8% of the beams as they are used, so each watt delivered
   leaves 1.47–1.49 W of heat on the Moon** (*screen*).
   - **The model.** ITU-R P.676-13's line-by-line model, coded here, reproduces ITU-R Study Group 3's own Earth values
     to 10⁻¹⁴ at fifteen frequencies from 1 to 350 GHz, and its slant-path example at 28 GHz through the mean annual
     global reference atmosphere, 0.4708 dB, to 10⁻¹².
   - **The column.** The moon_1.2atm column at 295 K holds 1.22 atm at the surface, 17.5% oxygen at Earth's partial
     pressure, 206 kg/m² of water vapour and its tropopause at 85 km. A 2.45 GHz beam straight down it loses 0.266 dB,
     5.9%, 7.4 times Earth's 0.036 dB. At 30° it loses 11.0%, at 20° 14.9% and at 10° 23%. A town at 10 km loses 4.4%
     overhead and 8.2% at 30°, and the summit at 11.6 km 4.2% overhead.
   - **The beams as used.** Over the lunar cycle at C′'s loads, every beam from the highest tile in view loses 6.3% and
     beams from any tile above 30° 7.7–7.8%: the towns 5.1–6.3% and the ground 6.8–8.3%.
   - **The link.** The DC-to-DC efficiency falls from 0.509 to 0.470–0.477. The air absorbs 2.1–2.6 TW at the low-cost
     end and 9.7–12.1 TW at the high-cost end, so the array's heat with close computing rises 6.7–8.4%, to 32.5–33.0 TW
     and 155–157 TW, and the dimming that answers it by 0.02–0.17% of the sunlight.
   - **At 5.8 GHz** a beam loses 7.2% overhead and the ground's beams 9.9% as used.
   - **Against the synthesis's estimate.** It put the loss at 5–8%, scaling Earth's oxygen by the pressure squared and
     the scale height, and the computed loss lies inside that range. Oxygen at Earth's partial pressure absorbs 17%
     less than oxygen at Earth's mole fraction would (0.266 against 0.317 dB overhead).
2. **Rain barely touches 2.45 GHz, and storms stop 5.8 GHz** (*screen*; the slow fall an *estimate*).
   - **The model.** ITU-R P.838-3 reproduces its own Table 5 at 2–6 GHz and the validation file at 14.25 and 29 GHz.
   - **Slow rain.** At a sixth of the gravity raindrops fall 2.9 times slower (the CM1 runs' microphysics), so a
     lunar millimetre an hour holds 2.9 times the water of Earth's and absorbs as 2.9 mm/h would there. Through 10 km
     at 50–100 mm/h a 2.45 GHz beam loses 0.26–0.54 dB; Earth's fall gives the synthesis's 0.1–0.2 dB. At 50 mm/h a
     5.8 GHz beam loses 96% through 10 km, against 44% with Earth's fall.
   - **The storms.** Rain above 1 mm/h covers 1.4% of the equatorial ring and 0.18–0.22% of the 45° rings; the 80°
     rings, under the GCM's sinking air, have none. A ground beam at its operating elevation, 61° on average over C′'s
     loads, passes through a storm 2.4% of the time at the equator and 0.33% at 45°, a town's at 71° 1.8% and 0.24%.
     At 30° they cross storms 1.7–1.8 times as often.
   - **The losses at 2.45 GHz.** The median storm (3.7 mm/h at its peak) costs 0.5–1.1%. A storm at the 90th percentile
     (18 mm/h, 90 km wide) costs 3.6–6.1% at the operating elevations and 6.6–10% at 30°. The largest storm on the
     ring (69 mm/h, 264 km wide) costs 14–23%, and 24–36% at 30°. At 5.8 GHz the 90th-percentile storm takes 68–97%
     and the largest all of it.
   - **The margin.** Three beams each rated for half the load carry 1.5 loads, 1.76 dB of headroom; only the largest
     storm, met at 30°, outruns it, and by 4%.
3. **The stack's tilted rings keep tiles above every horizon, and beams stand above 30° at all times to 60° of
   latitude** (*screen*).
   - **At the equator** 13.3 million of the fleet's 28.2 million tiles stand above the horizon and 8.6 million above
     30°; the highest averages 89° and never falls below 81°.
   - **At 45°** 4.7–5.7 million stand above 30°, and the highest averages 58–63°, never below 39–43°.
   - **Poleward of 60°** times come with no tile above 30°: 3% of the time at 62.5°, 26% at 72.5°, 76% at 82.5° and
     all the time at the pole, where 8.2–9.3 million tiles still stand above the horizon, the highest 20–25° up.
   - **The correction.** The synthesis and the network brief put the ring below the horizon within about 5° of its
     plane's poles. That holds for the middle ring alone. The held matched stack's rings are tilted up to 24° from the
     Moon's orbit plane, so their sunward and anti-sunward stretches climb 24° from it, and that plane stands 6.7° from
     the lunar pole.
4. **A tile stays above 30° for about 20 hours, and a receiver hands a beam over every 2–7 hours** (*screen*).
   - **Passes.** The retrograde rings turn round the Moon in 2.8–3.2 days, 2.7 days against its turning surface. A
     tile passing overhead is in view 30 hours and above 30° for 20. Passes above 30° last 20 hours at the equator, 15
     at 45° and 8 at 72.5°.
   - **Handovers.** With three beams each riding its tile through a pass, a receiver hands over a beam every 6.7 hours
     at the equator, 5 at 45° and 3 at 70°. Tiles follow one another along a ring every 16–19 seconds, so a beam can
     also step from tile to tile and stay near the top of the sky.
5. **The intake takes 47–224 TW of electricity at the tiles, which the dimmer's light covers at the low-cost end, and
   3,500–25,000 transmitting tiles, 0.013–0.09% of the fleet** (*screen* on *design guesses*).
   - **Power.** The radio power leaving the tiles is 33–157 TW, the intake over the share the air leaves, and 50–235
     TW rated. At the transmitters' 70% it takes 47–224 TW of electricity. The light the planned 5% dimmer takes gives
     117–176 TW (the shield's integrated ledger), so the high-cost end also draws on collectors in light the Moon never
     receives, which the night half intercepts at 851 PW (joint synthesis).
   - **A tile's collector.** At the solar constant a collector gives 0.27–0.41 GW of electricity and 0.19–0.29 GW of
     radio power per km². A whole tile, its 97.8 km² of clear face a sandwich of sun-pointing collector and phased
     array, beams 18.6–28 GW; the dimmer's 5% share of a tile's light beams 0.9–1.4 GW.
   - **Apertures** by the 1978 reference's tau of 2.04 (US DOE & NASA 1978; Rodenbeck et al. 2021), at a range of
     17,600–20,500 km: a town's 1 km crown needs a transmitter 5.0–5.8 km on a side (2.1–2.4 km at 5.8 GHz), and a
     3 km ground patch 1.5–1.7 km. A whole tile's aperture focuses to 500–580 m and holds the town beams at a tau of
     3.5–4.1, where almost nothing spills; shaped flat over a crown, a town's beams keep the synthesis's 15–82 W/m².
   - **Beams.** The 696,000 towns and 10,700–48,300 ground patches of 3 km, each 2.07 GW at 230 W/m², take three beams
     each: 2.1–2.2 million beams, a town's carrying 4–22 MW and a ground patch's 0.69 GW.
   - **Tiles.** At half their rating in use (*design guess*), whole-tile sandwiches number 3,540–25,300, 0.013–0.09%
     of the fleet, each carrying 88–600 beams, with 0.35–2.5 million km² of collector and aperture installed and
     0.12–0.82 million km² of collector in use. With collectors over a tenth of each transmitting tile the fleet needs
     35,000–253,000 of them, 9–60 beams each.
   - **Their heat.** The collectors' and transmitters' losses, 105–830 TW released in orbit, put 0.005–0.04 W/m² on
     the Moon, 0.1–0.4% of the fleet's glow.
   - **The high latitudes.** Spread evenly through the stack, the whole-tile transmitters above 30° over a receiver at
     57.5° fall to 7–48 at the worst moment (67–475 with tenth-tile collectors), and to none poleward of 60°.
     Transmitters placed on the most tilted rings serve those homes best.
6. **The receivers at the loads take land equal to 36–59% of C′'s settlements, and the high-latitude homes take their
   own low beams, which lose less than lines from 60°** (*screen* on *design guesses*).
   - **Where the loads are.** C′ with the blend places the ground's 4 billion as the metropolis's 0.1 billion, the dry
     highland cities' 1.58 billion at 30–45°, 0.61 at 45–60° and 0.41 at 60–90°, and the scattered communities' 0.62
     billion within 30° and 0.68 billion on the high-latitude lowland homes poleward of 48.6° (the middle of each pair).
     A fifth of the ground's load lies poleward of 60°.
   - **Land.** At 230 W/m² the ground's 22–100 TW take 96,000–436,000 km², 36–59% of the 264,000–736,000 km² that C′
     settles. One 3 km patch stands for every 15–25 km² of settled land, 3.9–5.0 km apart. A patch's edge, blurred by
     a whole tile's spot, falls below the public 10 W/m² within 0.18–0.20 km (*estimate*), so the closed ground grows
     by a quarter, to 120,000–561,000 km², 45–76% of the settled land.
   - **The high latitudes.** Every band to 60° sees a tile above 30° at all times. From there a line to a home at
     62.5–87.5° runs 150–910 km and loses 2.3–5.2% with its two converters (US EIA 2018), on top of the 8.5% the
     57.5° beams lose: 10.7–13.3%. The homes' own beams at 20–45°, from receivers 4 km up, lose 8.8–11.5% (9.9–12.9% at
     sea level), so the polar homes take their own beams and the grid backs them.
   - **Backing up a receiver.** Lines rated for its 2.07 GW to two neighbours 4–5 km away lose 2%, almost all in the
     converters.
7. **Sky reefs cross the beams' columns every day, so beams step aside as a matter of routine** (*estimate*).
   - **The columns.** A beam crosses the reefs' layer at 10–20 km along a strip: 1 km by 4.4 km for a town's beam at
     its operating 71°, 3 by 8.5 km for a ground patch's at 61°, and 18–20 km long at 30°.
   - **Meetings.** At 1% cover, half of it reefs a kilometre or two across drifting 0.1–1 m/s (the synthesis's rule), a
     town's beam meets 0.08–2.5 reefs a day and a ground beam 0.13–4.6. Left in place, a beam would have a reef in its
     column 3–7% of the time for a town and 8–22% for the ground, 9–23% and 17–43% at 30°. A crossing lasts 1–26
     hours at the operating elevations and up to 47 hours at 30°.
   - **Switching.** A receiver therefore moves one of its beams 0.2–14 times a day at 1% cover and twice as often at
     2%, to tiles whose columns miss the reefs the network's radar follows. Clearing a reef 5–10 km above a town takes
     10–31° of offset, near the synthesis's 10–20°.
   - **What a reef takes.** Caught by a late switch, a reef receives 5–27 W/m² from a town's beam (8–41 at its rating)
     or 77 W/m² from a ground beam (115 at its rating), 1–18% of the 650–730 W/m² of overhead sun in its layer, in a
     band that wet tissue absorbs.
8. **Single tiles fail every day without touching supply** (*design guess*).
   - At one whole-tile outage a decade per transmitting tile, 1–7 transmitting tiles fail a day (10–70 with
     tenth-tile collectors). Each receiver loses a beam 0.3 times a year, and its other two carry the load at their
     rating.
   - Two of a receiver's three beams fail within the 10 minutes a fresh beam takes once in about 900,000
     receiver-years: 0.8 times a year across the Moon, each time leaving one receiver at half its load for minutes.
9. **A fault across the whole ring, the design's one common point of failure, leaves 2.5–7.5 TW of critical loads, and
   a week of them takes 425–1,270 TWh** (*design guess*).
   - **The critical loads.** A tenth of the residents' power, 0.8–2.4 TW, for air, water and sewage pumps,
     refrigeration of food and medicine, safety lighting, communications and health care; a fifth of the towns'
     control, 0.06–0.55 TW, the towns drifting with the wind and steering only to keep clear; all the water lift; a
     tenth of the food's living cultures, 0.29–1.2 TW; a twentieth of the lamp farms, dark with their air and water
     moving, 0.22–2.2 TW; and a quarter of close computing, 1.2 TW, the millisecond control kept on the Moon. Together
     they are 7–11% of the 22–105 TW delivered. The towns' hydrogen remaking and their trips are shed.
   - **Days.** Seven, Earth's practice for safety-related emergency power (US NRC Regulatory Guide 1.137); each day
     takes 61–181 TWh.
   - **The towns** carry their own 0.46–1.75 TW for the week on 4.6–17.7 Mt of the 29.8 Gt of hydrogen that lifts
     them, burned in fuel cells at half its heating value: 0.015–0.06% of their hydrogen and 0.1–0.4% of their reserve
     lift. A week unremade loses another 0.02–0.2% of it.
   - **The ground's** 2.1–5.8 TW for the week, 348–972 TWh, take 21–58 Mt of stored hydrogen, made with 870–2,430 TWh
     of electricity at a 40% round trip (Headley & Schoenung 2020); or 258–720 km³ of water through Korolev's 3.74 km
     head, 1.8–5.1 m of its level; with batteries of 8–50 Mt carrying the first hour.
   - **Food.** Production stops while the cultures are kept alive, and the week's food comes from stores, which the
     field crops' harvest once a lunar day already requires (README, finding 5).
10. **The summit metropolis takes 0.25–1.06 TW on 1,500–6,300 km² of receivers outside the city, with a week's
    critical store of 5–15 TWh** (*screen*; its reserves on *design guesses*).
    - **Its load.** Its 100 million residents at C′'s use per person take 0.10–0.30 TW, which holds the brief's 200 GW;
      close computing 0.06 TW, and its share of the blend's food power 0.09–0.70 TW.
    - **Receivers.** Its 0.34–1.46 TW of intake at 230 W/m² take 166–705 patches, a band 7–23 km wide round the 32 km
      city, beyond its streets and its flyers' airspace.
    - **Beams.** At 5.4° N the highest tile averages 89° and never falls below 78°. From the summit's ground at 11.6 km
      the beams lose 5.2%, so each watt leaves 1.45 W of heat.
    - **Reserves.** The critical 31–87 GW for a week take 5.2–14.6 TWh: 3.8–10.8 km³ through Korolev's head or
      0.3–0.9 Mt of hydrogen. The night needs no store, since the night side's tiles beam through it, and the beams
      make no light; the first screen's 71 TWh night store was for local power alone.
11. **What changes in the synthesis's figures.**
    - **Absorption:** the 5–8% estimate holds, at 6.3–7.8% for the beams as used. The heat per watt delivered becomes
      1.47–1.49 W against the 1.38 W the synthesis takes from Earth service, and the array's heat with close computing
      32.5–33.0 TW and 155–157 TW against 30.4–145 TW.
    - **Supply:** the synthesis sets the dimmer's 117–176 TW of electricity against the Moon's 30–145 TW of intake.
      That intake takes 47–224 TW of electricity at the tiles, so the dimmer's light covers the low-cost end and the
      night half's light the rest of the high-cost end.
    - **The polar caps:** no cap is blind; poleward of 60° the beams come from 20–45° up.
    - **Rain:** with the Moon's slow fall, 2.45 GHz through 10 km at 50–100 mm/h loses 0.26–0.54 dB and 5.8 GHz at
      50 mm/h 96%, which strengthens the choice of 2.45 GHz for towns and ground alike.

## The deep air's absorption

**The model** is ITU-R P.676-13 Annex 1, coded in [power.py](power.py) with its Tables 1 and 2 in full:
- specific attenuation, equations 1–9: 44 oxygen lines with line mixing, the non-resonant Debye spectrum and
  pressure-induced nitrogen, 35 water-vapour lines including the 1,780 GHz pseudo-line that carries the wet continuum;
- slant paths, equations 13–19b: layers 10 cm thick at the receiver, thickening by e^(1/100), to the column's top at
  1 Pa, about 577 km, with Snell's law in polar coordinates on the refractivity of ITU-R P.453-14 and the Moon's
  radius.

**The Earth check.** The model reproduces the validation file's gamma_o and gamma_w at 1,013.25 hPa, 288.15 K and
7.5 g/m³ to 4×10⁻¹⁵ at fifteen frequencies, and its 28 GHz path at 30° through 922 layers of the P.835-6 mean annual
global atmosphere, 0.47081 dB, to 1.4×10⁻¹². Through the same atmosphere Earth's zenith takes 0.036 dB at 2.45 GHz and
0.040 dB at 5.8 GHz (*screen*).

**The Moon's column** is the atmosphere domain's moon_1.2atm scenario, rebuilt here with its own code: 1.2 atm of dry
air with Earth's 21.2 kPa of oxygen, 400 ppm of CO₂, a saturated moist adiabat at 295 K with Manabe–Wetherald humidity,
and a 200 K stratosphere. At 290 and 300 K it reproduces the product's water columns (136 and 306 kg/m²), tropopauses
(75 and 96 km) and surface pressures. Within 2.2% of pressure and 1.3 K of temperature it matches the inhabited-volume
air at 10–40 km; at the ground the GCM's air there holds 1.18 atm against the column's 1.22. P.676's oxygen line
strengths and Debye term take oxygen's mole fraction over Earth's 0.2095, 0.833 here, and its widths the total pressure.

Share of a beam lost, by elevation and receiver height (*screen*):

| Elevation | 2.45 GHz, sea level | from 4 km | town at 10 km | summit, 11.6 km | 5.8 GHz, sea level | town at 10 km |
|---|---:|---:|---:|---:|---:|---:|
| 90° | 5.9% (0.266 dB) | 5.2% | 4.4% | 4.2% | 7.2% (0.324 dB) | 5.0% |
| 60° | 6.8% | 6.0% | 5.0% | 4.8% | 8.2% | 5.7% |
| 45° | 8.2% | 7.2% | 6.0% | 5.7% | 9.9% | 6.9% |
| 30° | 11.0% | 9.8% | 8.2% | 7.8% | 13.3% | 9.3% |
| 20° | 14.9% | 13.3% | 11.1% | 10.6% | 18.0% | 12.7% |
| 10° | 23% | 21% | 17% | 17% | 28% | 20% |
| 5° | 31% | 28% | 24% | 23% | 38% | 27% |

Over a 2.45 GHz beam's path from sea level the water vapour takes 3% of the absorption; at 5.8 GHz it takes 15%.
From 290 to 300 K the zenith loss at 2.45 GHz falls by 5% of itself: oxygen's absorption weakens with warmth more than
the water's grows.

The beams as used, by latitude over the lunar cycle, beams drawn from any tile above 30° (or the highest tile where
none stands above it) (*screen*):

| Latitude | Highest tile, mean / least | 2.45 GHz, sea level | from 4 km | town | 5.8 GHz, sea level | town |
|---|---|---:|---:|---:|---:|---:|
| 2.5° | 89° / 81° | 7.4% | 6.6% | 5.5% | 9.0% | 6.3% |
| 17.5° | 85° / 66° | 7.7% | 6.8% | 5.7% | 9.3% | 6.5% |
| 32.5° | 73° / 52° | 8.3% | 7.4% | 6.1% | 10.0% | 7.0% |
| 47.5° | 58° / 39° | 9.1% | 8.1% | 6.7% | 11.0% | 7.7% |
| 57.5° | 48° / 31° | 9.6% | 8.5% | 7.1% | 11.7% | 8.1% |
| 67.5° | 38° / 24° | 10.3% | 9.2% | 7.7% | 12.5% | 8.7% |
| 77.5° | 30° / 20° | 11.6% | 10.3% | 8.6% | 14.0% | 9.8% |
| 87.5° | 24° / 20° | 12.9% | 11.5% | 9.6% | 15.6% | 11.0% |

The towns follow the design winds' gathering at 10 km (ways of living); the ground's loads follow its people, with the
dry highland cities' receivers at 4 km and the metropolis on its summit (*design guess*). Weighted so, the beams lose
5.1–6.3% at the towns and 6.8–8.3% on the ground, from the highest tiles to any usable tile.

## Rain

P.838-3's coefficients at 2.45 GHz are k = 1.35×10⁻⁴ and alpha = 1.055 for circular polarization, and at 5.8 GHz
4.80×10⁻⁴ and 1.602. Specific attenuation with the Moon's fall, the rate scaled by 2.94 (*estimate*):

| Rain rate | 1 mm/h | 5 | 10 | 25 | 50 | 100 |
|---|---:|---:|---:|---:|---:|---:|
| 2.45 GHz, dB/km | 0.0004 | 0.0023 | 0.0048 | 0.013 | 0.026 | 0.054 |
| 5.8 GHz, dB/km | 0.0027 | 0.036 | 0.11 | 0.47 | 1.4 | 4.3 |

The cloud-resolving rings (CM1, two-dimensional, 3-hourly snapshots over a lunar day) give the storms:

| Ring | Rain above 0.1 mm/h | 99th percentile | Wettest column | Storms in a snapshot | Rain above 1 mm/h covers |
|---|---:|---:|---:|---:|---:|
| Equator | 5.3–5.7% | 3.7–4.5 mm/h | 69 mm/h | 6.6, 24 km wide at the median | 1.4% |
| 45° N and S | 0.2–0.9% | below 0.1 mm/h | 8–23 mm/h | 0.6 | 0.18–0.22% |
| 80° N and S, with the GCM's sinking air | none | none | 0.006 mm/h | none | none |

The ecology work's sky bands put cloud over the equatorial ring 3.4% of the time at 10 km by day, the most at any
height a beam crosses, and 1.0–1.4% by night, in line with the storms' 1.4% of the ring and a slant path's 1.8–2.4%.

A storm's rain is taken at its peak surface rate from the ground to its top or the 0 °C level at 26.6 km, as wide as
the storm; the rain aloft is heavier, since most of the slow rain evaporates on the way down (*estimate*). Losses at
2.45 GHz and 5.8 GHz (*screen*):

| Storm | Ground at 61° | Ground at 30° | Town at 71° | Town at 30° |
|---|---|---|---|---|
| 99th-percentile rate, 4.5 mm/h | 1.2% / 16% | 1.3% / 18% | 0.6% / 8.6% | 1.2% / 16% |
| Median storm, 3.7 mm/h | 1.0% / 12% | 1.1% / 13% | 0.5% / 6.3% | 0.9% / 12% |
| 90th-percentile storm, 18 mm/h, 90 km | 6.1% / 86% | 10% / 97% | 3.6% / 68% | 6.6% / 88% |
| Largest storm, 69 mm/h, 264 km | 23% / all | 36% / all | 14% / all | 24% / all |

## Tiles and beams

The held matched stack of 1,859 rings, sampled here with 93 rings at 1° of argument, each ring circular at its layout
radius (19,360–21,390 km) and tilt (−24° to +24°), its nodes at right angles to the Sun and its tiles moving
retrograde. The Moon turns under it once a synodic month while the ring plane's pole circles the lunar pole at 6.7°;
the screen samples both. Tiles in view, in millions of the fleet's 28.2 million (*screen*):

| Latitude | Above the horizon | Above 30°, mean / least | Above 45° | Time with none above 30° | Pass above 30° | Handover |
|---|---:|---:|---:|---:|---:|---:|
| 2.5° | 13.3 | 8.6 / 8.4 | 6.3 | 0 | 20 h | 6.7 h |
| 22.5° | 13.2 | 8.0 / 7.3 | 5.2 | 0 | 19 h | 6.3 h |
| 42.5° | 13.0 | 5.7 / 3.8 | 1.9 | 0 | 15 h | 5.1 h |
| 52.5° | 12.7 | 3.5 / 0.97 | 0.63 | 0 | 13 h | 4.5 h |
| 57.5° | 12.5 | 2.4 / 0.05 | 0.29 | 0 | 12 h | 4.1 h |
| 62.5° | 12.2 | 1.45 / 0 | 0.11 | 3% | 11 h | 3.7 h |
| 72.5° | 11.1 | 0.39 / 0 | 0 | 26% | 8 h | 2.7 h |
| 82.5° | 9.7 | 0.02 / 0 | 0 | 76% | 5 h | 1.8 h |
| 87.5° | 9.3 | 0 | 0 | all | | |

The beams' sizing (*screen* on the *design guesses* below):

| | Low-cost end | High-cost end |
|---|---:|---:|
| Intake at towns and ground | 30.4 TW | 145 TW |
| Radio power leaving the tiles | 33 TW | 157 TW |
| Rated, three half-load beams a receiver | 50 TW | 235 TW |
| Electricity at the tiles | 47 TW | 224 TW |
| The dimmer's light at 30% / 20% | 176 TW | 117 TW |
| Ground patches of 3 km | 10,700 | 48,300 |
| Beams | 2.12 million | 2.23 million |
| Whole-tile transmitters at half their rating | 3,540 | 25,300 |
| Their share of the fleet | 0.013% | 0.09% |
| Beams per transmitter | 600 | 88 |
| Collector and aperture installed | 0.35 million km² | 2.5 million km² |
| Transmitters with collectors over a tenth of the tile | 35,400 | 253,000 |

Each beam may use a transmitter's whole aperture, so the beams a tile carries are set by its power and its
beamforming channels. The Moon lies within 5° of a tile's nadir, so a radial-facing aperture steers its beams with
almost no loss of area.

## The grid

The receivers sit at the loads, the ground's intake spread as its people. Poleward of the last band that always sees a
tile above 30° the loads could take power by line from 57.5° (*screen*):

| Home | Line from 57.5° | Line and converters | Through the 57.5° beams | Own beams from 4 km | at sea level |
|---|---:|---:|---:|---:|---:|
| 62.5° | 150 km | 2.3–2.5% | 10.7–10.9% | 8.8% | 9.9% |
| 72.5° | 450 km | 2.9–3.6% | 11.2–11.8% | 9.6% | 10.8% |
| 87.5° | 910 km | 3.8–5.2% | 12.0–13.3% | 11.5% | 12.9% |

So the high homes keep their own receivers, with low beams, and the grid ties them to the band equatorward for backup.
Within the settled land a receiver's two neighbours stand 4–5 km away, and lines rated for its 2.07 GW lose 2% through
their converters. Across a storm the grid moves power round the 24–264 km a storm spans.

## Reserves

| Event | What it does | What holds the load |
|---|---|---|
| A storm over a receiver | All its beams cross the same rain: 0.5–1.1% in the median storm, 4–10% at the 90th percentile, 14–36% in the largest (2.45 GHz) | The beams' 1.5-load rating, 1.76 dB; the grid beyond it |
| A sky reef in a column | Each beam meets reefs 0.08–4.6 times a day at 1% cover | The beam moves to a tile whose column misses the reef, make-before-break |
| A tile lost | 1–7 transmitting tiles a day; each receiver loses a beam 0.3 times a year | Its other two beams at their rating; a fresh beam within 10 minutes |
| A fault across the ring | All beams stop | Critical loads on stored energy for a week |

The critical loads for a fault across the ring (*design guess*), TW delivered:

| Load | Share kept | Low-cost end | High-cost end |
|---|---:|---:|---:|
| Residents, ground and aloft | 10% | 0.80 | 2.39 |
| The sky towns' control | 20% | 0.055 | 0.55 |
| Water lift | all | 0.001 | 0.028 |
| The power-made food's living cultures | 10% | 0.29 | 1.17 |
| The lamp-lit farms | 5% | 0.22 | 2.23 |
| Close computing | 25% | 1.16 | 1.16 |
| **Total** | | **2.53** | **7.53** |
| A week | | 425 TWh | 1,266 TWh |

| Store | Low-cost end | High-cost end |
|---|---:|---:|
| The towns' lifting hydrogen in fuel cells | 4.6 Mt, 0.015% of 29.8 Gt | 17.7 Mt, 0.059% |
| The ground's stored hydrogen | 21 Mt | 58 Mt |
| or pumped water through Korolev's 3.74 km head | 258 km³ | 720 km³ |
| Batteries for the first hour | 2.5 TWh, 8 Mt | 7.5 TWh, 50 Mt |

## The summit metropolis

| | Low-cost end | High-cost end |
|---|---:|---:|
| Residents, 100 million | 0.10 TW | 0.30 TW |
| Close computing | 0.06 TW | 0.06 TW |
| Its share of the blend's food power | 0.09 TW | 0.70 TW |
| Delivered | 0.25 TW | 1.06 TW |
| Intake | 0.34 TW | 1.46 TW |
| Receivers at 230 W/m² | 1,500 km², 166 patches | 6,300 km², 705 patches |
| A band round the 32 km city | 7 km wide | 23 km wide |
| Critical loads | 31 GW | 87 GW |
| A week's store | 5.2 TWh: 3.8 km³ through Korolev or 0.31 Mt of hydrogen | 14.6 TWh: 10.8 km³ or 0.88 Mt |

The summit's highland is the rainiest ground on the Moon, 9–17 mm a day in the GCM's cells round it
([metropolis brief](../habitation/summit_metropolis/README.md)), and its storms reach 22 km. Its beams come in from
near overhead, where the rain path is shortest, and its receivers sit on the closed slopes beyond the city and its
flyers' airspace, the lakes within 25 km taking 24% of the ground.

## Assumptions

**Design guesses**, each labelled in the product:
- beams planned from 30° up;
- three beams a receiver from different parts of the ring, each rated for half the load, switched make-before-break;
- ground receivers as flat-filled patches 3 km on a side at the 1978 reference's 230 W/m²;
- transmitters as sandwiches of sun-pointing collector and phased array over a tile's clear face, at half their rating
  in use; the variant with collectors over a tenth of each transmitting tile;
- one whole-tile outage a decade per transmitting tile, and 10 minutes to put a fresh beam on a receiver;
- the critical shares in the table above, a week of storage after Earth's 7-day practice, fuel cells at half
  hydrogen's heating value, batteries of 0.15–0.3 kWh/kg bridging one hour;
- the ground's loads where its people are, the dry highland cities' receivers at 4 km, the middle of each C′ pair for
  the composition, and reefs 5–10 km above a town for the switch's offset.

**Estimates:** the Moon's slow rain scaled to Earth's rate holding the same water, with Earth's drop spectra; storms as
uniform cells at their peak surface rate up to their tops or the 0 °C level; circular rings, leaving out the layout's
forced eccentricities of 0.11–0.16, which swing the distances by that share about each orbit.

**Published values:** HVDC losses of 2–3.5% per 1,000 km and 1% a converter (US EIA 2018); the public reference level
of 10 W/m² (ICNIRP 2020); hydrogen storage's 40% round trip (Headley & Schoenung 2020); a 7-day fuel supply for
safety-related emergency power (US NRC Regulatory Guide 1.137); hydrogen's 33.3 kWh/kg.

## Limits

- P.676's line parameters are Earth's air's. The screen changes the oxygen's share; the Moon's air is otherwise Earth's
  dry air. The recommendation is built for Earth's air up to its surface pressure, and the Moon's lowest kilometres
  reach 1.22 atm.
- The column is clear. The clouds' liquid water is left out; below 10 GHz it absorbs in proportion to the water held,
  as rain does. The heavier rain aloft, the melting layer and the GCM's regional climates are left out too.
- The storms come from two-dimensional rings, which organise convection into lines more readily than real storms, and
  from surface rates; the rain's drop sizes at lunar gravity are Earth's in P.838-3.
- The ring is sampled, circular, with the Sun in its plane; the year's 5° swing of the Sun about it and the stack's
  eccentricities are left out.
- The beams are counted. Placing them on tiles through time, against the tiles' ratings, the reefs moving on the design
  winds, the storms and the faults, is next work.
- The loads' siting by latitude follows people, and the food plants' own siting is open.
- The reefs' encounter rule is the synthesis's, the reefs spread evenly through 10–20 km.

## Next work

1. **A beam assignment model through the lunar cycle:** tiles against receivers with their ratings and visibility,
   reefs moving on the design winds, storms from the rings and faults drawn at their rates. Minutes on one core.
2. **The transmitting tiles:** sandwich or separate collector, size, beamforming channels and beams per tile,
   transmitters weighted to the tilted rings for the high latitudes, sidelobes and the spill round each receiver
   against 10 W/m² at the closed ground's edge.
3. **The receivers' land:** 36–59% of C′'s settled land, 45–76% with the ground closed round each patch, which the
   food plants' and lamp farms' roofs, the lakes and closed slopes could carry, sited with the cities' look and their
   flyers' airspace.
4. **Reefs and beams:** what a late switch does to a reef's tissue, and the switching rate under the author's
   preference for a somewhat higher reef density.
5. **The fault across the ring:** its causes and the time a restart takes, which set the days of storage; the stores'
   siting (hydrogen caverns or tanks, Korolev's head, the towns' fuel cells).
6. **Rain aloft:** the CM1 storms' liquid water by height, for fades from the products themselves.

## How to run it

```
python -m provisioning.power                          # writes results/power.json in about ten seconds
python -m pytest provisioning/tests/test_power.py
```

The product (schema `terluna.provisioning.power/1`) binds its producer, the atmosphere column code it imports and its
22 inputs by sha256, with the evidence statement, the reading rule and the assumptions. The tests regenerate it,
check the inputs' hashes, check the absorption and rain models against the ITU's own Earth values, and check the sums
this note states.
