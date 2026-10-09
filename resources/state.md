# Where resources start: what the project already holds

An evaluation of the repository at `main` `996ea6c` (9 October 2026) for resources. Paths are from the repository
root. The full working digest, with line references, is kept with the domain's run records on the research drive
(`research/runs/resources/repo_digest.md`).

## Decisions that bind resources

From [research/decisions.md](../research/decisions.md) unless named:
- **Duration:** construction within 500 years of 2026, and operation for at least 10⁹ years.
- **The air:** 1.2 atm with no dome, 82.5% N2 and 17.5% O2 by mole.
- **The volatiles:** water, nitrogen and oxygen come from the Solar-System-wide operation. Conservation D4 adds that
  lunar oxygen arrives only as a by-product of local metal and glass production; shipped as water, the oxygen leaves
  its hydrogen over. Polar ice is sampled and stored (D5).
- **The water:** standing water covers 28% of the surface (D10: 1.07×10¹⁹ kg), and the corrected climate's rain-fed
  lakes cover 10.3% more.
- **Protection:**
  - It is a formation of replaceable units, supplied over the long term.
  - Candidates are weighed on hardware, fresh feedstock, recycled throughput and delivery.
  - Requirement S7 keeps the protection systems' outflows out of the escape region.
- **Loss budgets:** 1–100 kg/s, each stated with its cycle time: about 100 billion years at 1 kg/s and 1 billion at
  100 kg/s.
- **Traffic:** the six rules for packets of 10⁶–10⁸ kg, captured and braked remotely, with Earth never behind the
  lunar catcher.
- **The ring fleet's form:** annulus tiles of about 13 g/m², carrying 0.1 µm of titania on 4 µm of silica.
- **Elsewhere:**
  - fusion supplies most local power (a planning assumption);
  - the sky ferries and liners are hydrogen ships;
  - the planned branch `research/array-industry` takes the film plant and the array's industry.
- **Conservation principle 5** (proposed): "Draw on common materials and protect unique ones". Water and gases come
  from ordinary icy bodies; Saturn's rings, Titan's atmosphere and the ocean worlds stay untouched.

The author's direction that imports are part of the plan, and that the problem is closing each cycle locally, had no
row in the register. This branch adds it to the register's section "The world", in the author's words.

## What is computed

| Result | Numbers | Evidence | Where |
|---|---|---|---|
| Air inventory | 3.14×10¹⁸ kg; column 78,700 kg/m², 7.6 Earth columns; scale height 51 km | Runnable baseline | `research/baselines/feasibility/reference.json` |
| Water | Seas 1.07×10¹⁹ kg; rain-fed lakes 3.1×10¹⁸ kg; the crust's pores 5.0×10¹⁸–5.9×10¹⁹ kg if saturated to 1–20 km | Hydrostatic storage and a first runoff estimate | `geography/results/atlas.json`, `drainage.json`, `groundwater.json` |
| CO2 | 400 ppm, 1.9×10¹⁵ kg (the study's constant-gravity 1.73×10¹⁵); an Earth-like biosphere holds 1.5–1.6×10¹⁵ kg | Literature synthesis | `research/studies/atmospheric_co2` |
| Weathering | 0.36–3.6 Gt CO2 a year net in the design-like case; the air from 49 to 28 Pa in 210–2,080 years unreturned; calcination 128 GW per Gt a year in theory | Literature synthesis with small calculations | same |
| Soil elements | Top metre (kg/m²) of mare, highland and KREEP-rich soils: Ca, Mg, Na, K, P; weathering exports per m² | Screening on Apollo averages | `research/studies/lunar_cycle_ecology` |
| Losses to space | 0.044–1.4 kg/s over the solar cycle, cycle times of 70 billion to 2 trillion years; hydrogen 0.0024–0.0027 kg/s | Runnable screening | `research/studies/joint_synthesis/results/ledgers.json`; the loss response |
| Shield film | 28.1–28.2 million tiles, 46 Gt; film remade at 2.3–4.6 Gt a year; the film plant 0.15–2.9 TW near 2,000 K | Arithmetic on a proposed design | `research/studies/solar_shield_array/results/array_heat.json`, `array_industry.md` |
| Build flows | The air over 500 years 1.99×10⁸ kg/s; nitrogen over 300 years 2.67×10⁸ kg/s; nitrogen transport 33–235 PW | Historical model | `research/baselines/feasibility`, `engineering/reference` |
| Erosion | 10–100 m per million years of land lowering, 0.71–7.1 Gt of rock a year; lifting it 3 km takes 0.11–1.1 GW | Screening on an assumed rate | `research/studies/conservation/results/gates.json` |
| Tonnages | The tower's chosen form 98 Mt of steel plus 58 Mt of floors; regional magnets 1–4 Gt; momentum stores 0.13–1.7 Gt; lift gas 48,000–71,000 t of hydrogen | Proposed designs | `summit_tower`, `solar_shield_array`, `sky_fleet` studies |
| Flooding | Mare shares under the 28% sea | Runnable | `geography/results/atlas.json` |

No product gives the Moon's rocks by element and region, the seas' composition, or a carbon, oxygen or element budget
through time.

## What other lines assume of resources

- **The shield:**
  - feedstock for 2.3–4.6 Gt of film a year, and for the initial 46 Gt;
  - the unrecovered fraction;
  - the fate of failed tiles;
  - freight to the film plant, which rides the surface-to-orbit link if the feedstock comes from the Moon.
- **The tower:** 98 Mt of steel; "which materials can be made on the Moon" is open.
- **The sky fleet:** 48,000–71,000 t of hydrogen, or a helium supply.
- **The magnets:** 1–4 Gt of conductor and structure.
- **The air:** resupply of 1.4–44 kt a year.
- **Carbon:** the biosphere's 1.5–1.6×10¹⁵ kg of CO2, "imported with the biosphere".
- **The biomes:** nitrogen from biological fixation, with the delivered N2 behind it.
- **The seas:** their composition and density, which the waves (1,025 kg/m³), the tides' self-attraction and the
  aerosol's spray all assume.
- **Developed soils:** clays, carbonates and some salt.
- **Fusion fuel.**
- **The habitat array:** its materials, so far unbudgeted.

## Open items already listed

- The CO2 study: Earth rates, ocean carbonate chemistry, how weathering slows, and the box model with a managed
  return.
- The biosphere: the rock cycle's return limb, at 0.3–3 times the world's salt output, and a land, lake and sea budget
  of every rock-derived element.
- Geography: how fast and how deep the crust takes up water; erosion and sediment.
- Engineering: source-to-use pathways, the dependency network, a shared material ledger, and the return of eroded
  sediment.
- Conservation: the form of imported oxygen and the use of its hydrogen.
- Habitation and the infrastructure review: a launch rule in the spirit of S7, where failed tiles go, and freight for
  the film plant.

## Statements elsewhere that need correcting

These belong to main and stay unchanged on this branch.
1. **Constant-gravity air masses** about 10% low sit beside the spherical 3.14×10¹⁸ kg:
   - geography's 2.8×10¹⁸ kg;
   - the middle-atmosphere README;
   - the ecology register;
   - the CO2 study's 35.6 Gt per Pa;
   - the infrastructure review's 105 Mt of hydrogen;
   - D4's water and hydrogen, which with 6.13×10¹⁷ kg of O2 become 6.9×10¹⁷ and 7.7×10¹⁶ kg.
2. **The CO2 study's runoff and temperature.** It uses 349 mm a year from run A and 26 °C, against about 265 mm and
   21.75 °C now. Its weathering sink runs about a quarter high on runoff alone.
3. **The ecology register's runoff** is 279,000 m³/s with lakes on 11.6%, against 229,410 m³/s and 10.3%.
4. **A unit error.** The joint ledger's `film_plant_w` holds terawatts under a label in watts (fixed on main on
   9 October, 78020e0).
5. **Two figures for film replacement,** 2.2–4.5 and 2.3–4.6 Gt a year.
6. **Erosion.** The ecology register says erosion "has not been estimated", while conservation and engineering use
   0.7–7 Gt a year from an assumed rate.
7. **Lunar oxygen.** The engineering seed and the locked plan-01 outline still name lunar oxygen routes; D4 replaced
   them.
8. **Carbonate chemistry.** The ecology register says calcium carbonate settling "draws CO2". Precipitating carbonate
   releases CO2, as the CO2 study's "returns half" has it.
9. **Principle 5.** It says Saturn's rings or Titan alone could meet the need. The rings' 1.5×10¹⁹ kg of ice falls
   short of the seas, lakes and any saturated crust (1.4–7.3×10¹⁹ kg).
10. **The fleet in the research README** is given as 27–33 million tiles and 45–53 Gt, against 28.1–28.2 million and
    46 Gt now.
11. **Water density.** The atlas and groundwater use 1,000 kg/m³, the waves and tides 1,025.
