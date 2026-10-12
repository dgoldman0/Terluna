# Growth over centuries

**A giant's age is the time its own light takes to make its gas: round giants of 1–2 km are 70–1,200 years old, and
colonies reach 2 km in 36–544 years.** Gas is 37–92% of the carbon an aerophyte spends building itself, and a round
body's gas and structure per square metre grow with its size. A round body's radius therefore grows almost linearly
and slows as its tendons grow; it stops where renewing them takes its whole surplus. A flat colony adds modules that
each earn their own keep, so its area grows exponentially, doubling every 3.6–55 years. Both live on the timescale of
great trees, and the round giant is the slower and older of the two.

Round bodies' ages start from a 40-m juvenile and colonies' from one founding 60-m module; all are *screen* values
from [run.py](run.py) (`ages`). The growth share is the part of the surplus a body spends on growth (1 or 0.5, a
design guess); the rest goes to reproduction and repair.

| Diameter | Round body | Colony of 60-m modules |
|---|---:|---:|
| 40–60 m | 1–7 years in a canopy nursery; 7 years on its own light | the founding module, 2–8 years in a nursery |
| 200 m | 17–304 years | 12–182 years |
| 1 km | 73–1,186 years (300 MPa fibres, or cellulose fibres in canopy-like light) | 29–435 years |
| 2 km | 225–702 years (300 MPa fibres in canopy-like light) | 36–544 years |

The low ends take canopy-like light (GPP 5 kg C/m²/yr), the four-mole acetate ceiling or photolytic hydrogen, and
all surplus into growth; the high ends the reference light (GPP 2), the measured 2.1 mol H₂ per mol glucose and half
the surplus. With photolytic gas at the reference light and all surplus to growth, a 300 MPa round body reaches 1 km
in 226 years and a colony 2 km in 83. The [structure](structure.md) note gives the fibres.

## What sets the pace

**Making a kilogram of hydrogen costs 3.8–9.5 kg of carbon by photolysis and 8.9–21.9 kg by fermentation.** Area and
time given to photolytic hydrogen fix no carbon, so 1 kg H₂ at 1% of 200 W/m² displaces GPP × LHV/(ηI) of carbon:
3.8 kg C at the reference GPP of 2 kg C/m²/yr and 9.5 kg at 5 (*derived*). Fermentation commits 22.3–54.8 kg of
glucose per kg H₂ across 4.0–1.63 mol/mol ([biology](biology.md)), 8.9–21.9 kg C.

**A colony's every new square metre costs 10.7–34.8 kg C and earns 0.89 kg C a year at the reference light.** The
60-m module carries 1.83 kg H₂ per square metre of raft, 5.6 kg/m² of dry structure and 1.1 kg/m² of living tissue.
At the reference light its surplus after upkeep, renewal and gas replacement is 0.89 kg C/m²/yr, so the raft doubles
every 8.4 years with photolytic gas and 27.3 years with fermented gas; in canopy-like light it doubles every 3.6–6.3
years. A 2-km raft is about 1,000 modules of 60 m.

**A round body's construction per square metre grows with size until its surplus runs out.** For the
cellulose-class body it rises from 9.4 kg C/m² at 40 m to 52 kg C/m² at 500 m while its surplus falls from 0.99 to
0.58 kg C/m²/yr; growth stops at about 850 m at the reference light and 1.9 km in canopy-like light. With 300 MPa
fibres it stops at 1.4 km and passes 3 km in canopy-like light. Integrating dE/dR against the surplus over the
growing collecting area gives the ages; for constant costs per area the integral reduces to t = (2C/S) ln(R/R₀),
which the tests check. An independent integration on a 400,000-step grid reproduced every age and stop radius to
0.1 year.

## The canopy nursery

**Megaforest crowns can raise juveniles to 50–60 m, the smallest round bodies that float with their reserve and wet
skins, in one to seven years at the measured 2.1 mol H₂ per mol glucose, and from nine months at the four-mole
ceiling; beyond that the gas outgrows any host, and the adult grows on its own light.** The host passes sugar through
a phloem connection and the juvenile ferments it, while the juvenile's own photolysis runs in full sunlight at the
crown top.

- **The host's production.** The canopy model's equatorial stand fixes 405.7 g C/m² gross and grows 228.4 g C/m²
  per lunar cycle, a carbon-use efficiency of 0.56: 5.0 kg C/m² gross and 2.8 kg C/m² net a year at 12.37 cycles a
  year (*derived*). At 60° it is 3.5 and 1.7 kg C/m². The model has no water limit; rain supports 26–58% of the crop
  model's growth within 15° of the equator ([people_and_land](../../../biosphere/ecology/people_and_land.md),
  finding 10), so a rain-limited host grows 0.73–1.6 kg C/m² net.
- **The crowns.** In the megaforest study's allometry a tree's crown radius is 0.18 of its height: 0.92 ha for a
  300-m tree and 2.5 ha for a 500-m tree, growing 26 and 72 t C a year at the clear-sky equatorial rate.
- **The share.** The project's day fruit already puts 30–80% of a stand's growth into fruit (the fruit model's
  shares), so a host passing 10–50% of its growth to juveniles is within the project's own designs (a design guess).
- **Fill times.** A 50-m juvenile holds 3.4 t of hydrogen at its 10-km trim and loses 38 kg a year through its skin
  while it fills. On a 300–500-m host passing 10–50% of clear-sky equatorial growth and fermenting at 2.1 mol/mol,
  with its own photolysis on half its area, it fills in 1.3–5.4 years; on rain-limited hosts, 2.0–6.6 years; on its
  own light alone, 7.2 years. A 40-m juvenile (2.1 t) fills in 0.9–4.7 years but floats with its reserve only after
  growing past about 45 m; a 60-m juvenile (5.0 t) takes 1.8–6.0 years. A colony's founding 60-m module (5.7 t) takes
  2.0–6.8 years.
- **The 200-m body.** Its trimmed gas, 81 t of hydrogen, is 1,380 t of sugar carbon at 2.1 mol/mol: 49
  hectare-years of clear-sky forest growth, or 188 rain-limited (26–243 across yields and growth). Its full envelope,
  343 t, would be 109–1,025 hectare-years. One 500-m host passing 10–50% of its growth would take 38–190 years to
  make the trimmed gas, 150–740 years if rain-limited.

**Host wind loads cap nursery juveniles at about 60 m on 300–500-m trees.** A juvenile at the crown top adds its drag
times the tree's height to the crown's overturning moment. A 40-m juvenile on a 300-m tree raises that moment by 24%
and lowers the reference tree's first failure (root, 23.5 m/s in uniform wind) to 21.1 m/s; a 60-m juvenile raises it
55%, to 18.9 m/s. On a 500-m tree the same juveniles add 9% and 20%. A 100-m tree carries juveniles of about 10–20 m
across: a 20-m juvenile adds 55% to its moment and a 40-m one triples it. Juveniles nested inside the crown at
0.6–0.9 of the height meet less wind.

**Release at dawn rides the morning into the sky.** Over land the CM1 ring's air at 10 m blows at 0.3 m/s at dawn, and
the mixed layer deepens from 0.2 km at dawn to 9.9 km just before noon, while the cloud base rises from 1.1 to 7.3 km.
A juvenile freed in the calm dawn climbs with the deepening mixed layer and reaches 10 km in the first week of its
first day, with the whole sunlit half-cycle ahead.

**The host's partnership has two clean forms.** As a mutualism, the aerophytes carry the host species' seeds or
spores across the Moon and return phosphorus as litter wherever they range. As one species with two generations, the
nursery stage is the aerophyte's own rooted or epiphytic generation that buds aerophytes, as a cnidarian polyp buds
medusae, and the question of a partner's gain disappears. Mistletoes show the route for water and minerals: xylem-
tapping hemiparasites draw water and inorganic nutrients through the haustorium and fix their own carbon, and the
leaf phosphorus of *Loranthus europaeus* tracks that of its oak host (Glatzel & Geils 2009). Dodder shows the route for
sugar: its phloem joins the host's and passes sucrose and amino acids without selection (Birschwilks et al. 2006).

## Phosphorus, up and back down

The contents, losses and routes below are design guesses, labelled once here, with their ranges.

**An aerophyte holds 1.2–4.5 g of phosphorus per square metre and must import 0.2–0.7 g a year.** Active tissue (1.1 kg
dry/m²) at 1–3 g P/kg holds 1.1–3.3 g; the mostly inert structure (2–4 kg/m²) at 0.05–0.3 g P/kg holds 0.1–1.2 g.
Turning over half the active tissue a year and resorbing 65% of its phosphorus, the mean resorption of terrestrial
leaves (Vergutz et al. 2012), loses 0.19–0.58 g/m²; renewing a tenth of the structure loses 0.02–0.09 g/m². A round
giant's structure is heavier, 60 kg/m² at 1 km with 300 MPa fibre and 167 kg/m² with cellulose-class fibre, mostly
tendon: at the same content it holds 3–50 g P/m² more, renewed once a century (*derived*).

**Routes up.**
- **Through the host, once.** A 50-m juvenile leaves the nursery with 2.4–8.8 kg of phosphorus, 0.7–5.5% of a 500-m
  host's yearly uptake of 160–320 kg (*derived*: 72 t C of growth a year at 1–2 g P per kg of dry growth).
- **Rain.** At 0.005–0.05 mg P per litre, the design climate's 2.24 mm/day brings 0.004–0.04 g/m² a year and the
  equatorial 5.5 mm/day 0.01–0.1 g: 1–50% of the need.
- **Airborne particles.** Collectors sweeping 1–10 µg/m³ of pollen, spores and dust holding 0.1–0.5% P at 0.1–1 m/s
  and 10–50% efficiency gather 0.0003–0.8 g/m² a year; sailing supplies the relative flow.
- **Roosting animals.** One 1-kg flyer per 100–1,000 m² passing 0.1–1 g P a day would bring 0.04–3.7 g/m² a year.
  Aerophytes as roosts for the sea–forest migrants close the cycle best, and are the most beautiful answer: the
  migrants carry the seas' and forests' phosphorus up to the aerophytes, as seabirds enrich their islands.

**Routes down.** Unresorbed litter (0.2–0.6 g/m² a year) falls wherever the aerophyte drifts; rain washes some from its
surfaces; offspring carry their stock away; and a dead aerophyte falls with all of it. A 2-km colony holds 3.8–14 t of
phosphorus, so its fall feeds a patch of forest or sea for years, as a whale fall feeds the deep sea floor.

## Lifespan, replacement and recruitment

**Reproduction is cheap for a giant.** A 50-m juvenile costs 8.7–32.8 kg C/m² over its 1,960 m², 17–64 t C, and a
colony's founding module 10.7–34.8 kg C/m² over its 3,118 m² cell, 33–109 t C. A 1-km round body in canopy-like light
earns 2,080–2,660 t C a year (2.6–3.4 kg C/m² over 0.79 km²), enough for a juvenile every two to eleven days. A
population replaces adult deaths when each adult releases its annual death rate divided by the juveniles' survival to
maturity: at a hazard of 0.001–0.005 a year and survival of 0.1–0.5, that is 0.002–0.05 juveniles per adult a year.

**Lightning sets the hazard of a round body, and colonies lose modules instead.** A 200-m body meets 0.02–0.13 strikes
in 300 years over Moon-mean storms and 0.05–1.4 over the summit's convective cells; a 2-km body 0.9–2.1 and 2.3–23
([storms](storms.md)). If every strike spread fire through a shared skin, a 200-m body would live 300 years with
probability 0.87–0.98 at the Moon-mean flash density and down to 0.24 over the summit's cells, and a 2-km body would
rarely reach 300 years over convective ground. Each body carries the skin water that keeps a burning cell's fire from
its neighbours ([structure](structure.md)), so a strike costs one cell. A 2-km colony of separate modules meets
0.3–0.7 strikes a century over Moon-mean storms and up to 7.5 over the summit's cells, loses the modules each strike
reaches and regrows them; its age is limited by fragmentation and drift into storms, which the present products
cannot rate.

## Sources and checks

[growth.py](growth.py) holds the pure functions; [growth_sources.json](growth_sources.json) records Vergutz et al.
(2012), Glatzel & Geils (2009) and Birschwilks et al. (2006) with what was read. Inputs are the canopy carbon cycle
([plant.json](../../../biosphere/canopy/results/plant.json)), the fruit model's shares
([fruit.json](../../../biosphere/canopy/results/fruit.json)), the megaforest reference cases
([reference_cases.csv](../megaforest_wind/results/reference_cases.csv)) and the CM1 equatorial ring.
[test_growth.py](test_growth.py) checks the hydrogen costs, the closed-form round age, the stop radius, exponential
doubling, the crown silhouette against a numerical integral, the wind-moment ratio, hydrogen conservation in the
nursery fill with the juvenile's gas loss, and the phosphorus arithmetic. The host's wind moment leaves out the
stem's drag, about a fifth of the tree's own moment, so the stated thresholds err on the cautious side.
