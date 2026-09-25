# Conservation of the Moon through its transformation

The Open Moon will be transformed into a living world, and this study shapes how that transformation keeps faith with what exists today. It sets out:

- the ethical basis;
- two preservation programmes, one scientific and one for heritage;
- replacements for the scientific functions the Moon loses;
- the order in which the transformation removes each record;
- a first heritage register;
- scientific sampling zones;
- candidate settings for communities.

Its geography comes from the [atlas](../../../geography/README.md#atlas-at-the-selected-water-share) at the selected water share, 28% ([shared/scenarios/water.json](../../../shared/scenarios/water.json)).

## Status and evidence

This is the working framework of 25 September 2026. It records the author's decisions, the proposals awaiting decision and the open questions.

The heritage register has two layers:

- **The written evaluation** ([sites.json](sites.json)): what each site marks, its values, how much its meaning depends on the exact spot, what survives, who it matters to, and a proposed treatment. This is human judgment. The current text is a first pass by Claude awaiting the author's review, and each designation remains open to the people each site matters to.
- **The computed situation** ([run.py](run.py)): each site's height, depth, pressure, distance to land or water, and the water body it lies in, from the atlas product.

The numbers are of two kinds:

- hydrostatic storage geometry, from the atlas;
- screening estimates from [screens.py](screens.py), each with its source and assumptions in its docstring. Domain models supersede these as they are built.

The sources in [sources.json](sources.json) were read as abstracts, summaries or primary web pages. Full-text reading and claim locators come before any manuscript citation.

## Reproduce

From the repository root:

```sh
python -m geography.fetch_inputs --download        # LOLA, GRAIL and the IAU gazetteer, hash-checked
OPENBLAS_NUM_THREADS=1 python -m geography.atlas   # geography/results/atlas.json and the grid product
OPENBLAS_NUM_THREADS=1 python -m research.studies.conservation.run
OPENBLAS_NUM_THREADS=1 python -m pytest geography/tests/test_atlas.py research/studies/conservation
python visualization/atlas/render.py               # map sheets and the 3D globe, kept out of Git
```

The runner writes these files to [results/](results/):

- [heritage_register.json](results/heritage_register.json) and a [CSV summary](results/heritage_register.csv);
- [science_targets.csv](results/science_targets.csv);
- [gates.json](results/gates.json), holding the screening numbers;
- [community_settings.json](results/community_settings.json);
- [checkpoint.json](results/checkpoint.json), holding the product and file hashes.

## 1. Starting point and principles

The transformation is given; conservation shapes how it is carried out. The work rests on ethics and sustainability, with law as supporting argument. Stewardship and decisions belong to people, to communities, and to nations understood as peoples.

The many nations of Earth are the progenitors of the Moon's nations. Those nations form as the Moon's own indigenous nations, centred on being of the Moon.

**Principles (proposed):**

1. **Record before changing.** Each irreversible step waits until the record it destroys has been taken.
2. **Keep some of the old Moon intact.** Samples stay sealed and unopened, and one or two small places stay sealed.
3. **Respect what places mean to people**, as those people understand it.
4. **Prefer low-upkeep measures that fail slowly.** More than one organization keeps each site, the archive and every life-support function.
5. **Draw on common materials and protect unique ones.**
   - Water and gases come from ordinary icy bodies.
   - Saturn's rings (about 1.5×10¹⁹ kg of ice), Titan's atmosphere (about 9×10¹⁸ kg, mostly nitrogen) and the ocean worlds stay untouched, although the rings or Titan alone could meet the need.
   - The project weighs its effects on Earth as part of its own costs.
6. **Leave room for the future.** Choices are reversible, records are open, and some land stays unassigned and some unmanaged.
7. **Decide with the consent of those affected**, in proportion to how much they are affected. *A rewording is pending: the author asked for one that gives those affected and those with relevant expertise even standing.*

## 2. Three kinds of value, handled three ways

| | Scientific preservation | Heritage preservation | Replacing lost functions |
|---|---|---|---|
| What is valued | Information about nature | What places and things mean to people | Services the Moon provides |
| Where the value lies | In material and records, which can leave the place | In places, their contents, their setting and continued care | In a function that can move elsewhere |
| How it is kept | Survey, excavate, sample and measure; store in original condition; research through the transformation | Designate sites; keep them in place where the place carries the meaning; enclose, move within the landscape, mark, replicate, document | Systems at least as useful, built elsewhere before the loss |
| When | Before the stage that destroys the record, with research continuing throughout | Designated early; cared for indefinitely | Ready before the loss |
| Who decides | Open research review; an archive open to anyone who can use it | The people and communities to whom it matters, including any community living there | Those who use the function |
| How it fails | A missed deadline makes the loss permanent | When care stops, decline is slow and visible by design | A gap in service |

**The test.** Take a complete sample and record of a thing, then ask what they carry:

- if the sample and record carry everything important, the thing belongs to scientific preservation;
- if its meaning stays with the place, it is heritage;
- if it provides a function, the function gets a replacement.

**Overlaps.** One thing can carry several kinds of value, and each programme handles its own part. At Tranquility Base:

- heritage governs the place;
- science takes minimal, documented samples where the heritage process agrees;
- the laser reflector's ranging moves to new stations, while the reflector stays as a heritage object.

## 3. What changes, what survives, what can be kept

### 3.1 Changes the transformation makes everywhere

These changes are the purpose of the project. The figures come from [gates.json](results/gates.json).

- **Vacuum becomes atmosphere.** Outside sealed enclosures, the natural exosphere, the black daytime sky, space weathering, and the solar wind and micrometeoroids reaching the ground all end.
- **Surface chemistry.**
  - Metallic and nanophase iron oxidise.
  - Glass and minerals hydrate and weather, and soluble coatings dissolve at first wetting.
  - The solar-wind rims of grains (60–200 nm) dissolve in 1.5–5 years at laboratory rates, or 10² to 10⁵ times slower at field rates.
- **Temperature.**
  - The month-long day–night swing is damped.
  - The polar cold traps warm to the global temperature (the GCM puts the poles near 295 K).
  - Warming reaches 18–56 m into the ground per century.
- **Small landforms.** Wind, rain, soil and vegetation erase or soften footprints, tracks, fresh crater rims and bright rays.
- **Lowlands flood.** At 28% water cover, the Near-side Sea and the South Pole–Aitken Sea form, and every depression below −1,654 m holds water.
- **The crust fills with water.** The upper crust is about 12% porous and becomes saturated to kilometres deep.
- **Life spreads everywhere.** The lifeless surface and its organic baseline from before life both end.
- **Quiet conditions end.**
  - Seismic quiet ends.
  - Radio quiet ends as an ionosphere forms.
  - Low lunar orbit ends, because the air reaches thousands of kilometres up.
- **The Moon seen from Earth.** With the GCM's planetary albedo of 0.287, full Moon grows about 2.0 times brighter and quarter Moon 3.2–4.4 times. Clouds appear.
- **Over millions of years.**
  - Land lowering of 10–100 m per million years carries 0.7–7 Gt of rock a year into the seas.
  - Lifting it back up 3 km would take 0.1–1.1 GW, a small fraction of the roughly 5×10¹⁷ W available at K ≈ 1.17.

### 3.2 What survives on its own

- **Large landforms and their names.** Tycho and Copernicus remain as landforms while their rays fade.
- **Deep geology**, below the reach of water and heat.
- **The orbit, month, phases, eclipses and librations**, so calendars keep their timing.
- **The broad pattern of the Moon's face.**
  - Sea covers 37% of the near-side hemisphere, against 31% for today's maria.
  - The great round maria fill completely: Imbrium, Serenitatis, Crisium, Nectaris, Nubium and Humorum.
  - Procellarum, Frigoris and Fecunditatis are mostly sea.
  - Tranquillitatis is about one-third sea, and Vaporum and Insularum stay land.
  - Basalt weathers to darker, iron-rich soils, so the maria that stay land probably still read darker than the highlands. This is a hypothesis for the illumination and biosphere work.
- **The 1651 names for the maria.** Most become literally true. The Sea of Tranquility becomes true in its southwest arm, where Tranquility Base lies.

### 3.3 What can be kept, and how

| Kept through | What | How |
|---|---|---|
| Scientific preservation | Samples of every geological unit; polar ice; ancient soils sealed between lava flows; rock from Earth; implanted gases; volcanic glass with its coatings | Excavation, coring and drilling; sealed storage in original condition (vacuum, dry inert gas, cryogenic); part of each category left unopened; copies in more than one place; open access |
| | Measurements that need the old conditions | Seismic, heat flow, electrical and magnetic sounding, gravity; exosphere, dust and plasma; surface photometry |
| | Complete surface records | Imagery, laser altimetry, radar and spectra of the whole surface; micrometre-scale records of heritage sites |
| | The transformation itself | Monitoring networks, sediment cores from the new seas, and the crust's response to the water load |
| Heritage preservation | Designated places | Undersea communities, enclosed land sites and markers (section 5) |
| | Surfaces and objects | Cut out and kept sealed at the site or nearby |
| | Access and memory | Towers, markers and replicas; full documentation of every site |
| | A sample of the pre-transformation surface | A sealed vacuum reference area (proposed, section 5.6) |
| Replacement | Calibration standard, radio quiet, ranging targets, low-orbit science | Section 7 |

### 3.4 What survives as records

- The vacuum Moon as a whole, and its natural sky.
- Space weathering and cold trapping as ongoing processes.
- Seismic and radio quiet.
- The Moon's present appearance from Earth.

## 4. Decisions

- **D1. The transformation is given (decided).** Conservation concerns how it is done.
- **D2. Two programmes (decided).** Scientific preservation and heritage preservation are separate programmes, each with its own criteria and justification.
- **D3. Replacements (decided).** Each scientific function the transformation removes gets a replacement at least as useful.
- **D4. Where the air comes from (decided).** Water, nitrogen and oxygen all come from the Solar-System-wide resource operation. Lunar oxygen arrives only as a by-product of local metal and glass production.
  - This corrects the September feasibility baseline, which ships the nitrogen and produces the oxygen locally.
  - Oxygen for the 17.5% design shipped as water comes to 6.2×10¹⁷ kg of water, about 6% of the 28% seas, leaving 7×10¹⁶ kg of hydrogen.
  - Making the same oxygen from rock would mean processing 0.5–1 million km³, two to three times the whole regolith.
- **D5. Polar ice (decided).** Polar ice is excavated, analysed and stored as samples before the poles warm. Its value is scientific: as water it would supply less than a millionth of the need.
- **D6. Tranquility Base (decided).** Tranquility Base becomes an undersea domed community.
  - It is built early as a surface dome that protects the site.
  - It becomes submerged when the sea fills.
  - Key surface material and hardware are cut out and preserved inside.
- **D7. Laser ranging (decided).** Ranging is a function, handled as a replacement.
- **D8. Who decides (decided).** Stewardship and decisions rest with people, communities and nations as peoples. Law serves as supporting argument.
- **D9. The Moon's nations (decided).** Earth's nations are the progenitors of the Moon's nations, which form as the Moon's own indigenous nations. The early missions' sites belong to the shared heritage of all of them.
- **D10. Water share (selected at the author's request).** Standing water covers 28% of the surface. Sea level is −1,654 m above the geoid, the equivalent global layer is 282 m and the water mass is 1.07×10¹⁹ kg. The basis is in [water.json](../../../shared/scenarios/water.json).
- **D11. Submerged heritage sites (decided).** Heritage sites that the seas cover become submerged communities. A community on a nearby dry coast can be a related, second community.

## 5. Heritage preservation

### 5.1 How the register is made

Heritage is a matter of meaning, so the register starts from written evaluation. For each site [sites.json](sites.json) records:

- what it marks;
- its historic, scientific, social, memorial and aesthetic values, after the Burra Charter;
- how much its meaning depends on the exact spot;
- what physically survives;
- who it matters to;
- a proposed treatment, with its reason;
- the status of that evaluation;
- whom to consult;
- open questions.

The runner then places each site on the atlas.

Where the evaluation departs from the author's standing rule, the entry says so. Luna 2, for example, left only scattered impact debris, so its evaluation proposes a marker. Such cases stay open for the author.

### 5.2 Categories

1. Sites of human activity on the Moon.
2. Memorials and human remains.
3. Natural landmarks, and the Moon's face as a cultural landscape.
4. Sealed reference areas of the pre-transformation surface.
5. Sites of the transformation and of lunar life as they arise.

The early missions open the register; most lunar heritage will be made later, and the Moon's nations will designate their own.

### 5.3 Criteria

- **Significance:** how widely and how deeply people hold the site meaningful.
- **Dependence on place:** whether the meaning is tied to the exact spot or travels with the objects.
- **Integrity:** how much of the original material and setting survives.
- **Wishes:** those of the people the site matters to, including any community living there, those who hold the Moon sacred, and families of those commemorated.
- **Upkeep:** how much continuing care the site needs, and how it fails when care stops.
- **Equal standing** for all nations (D9).

### 5.4 Methods, most preferred first where the place carries the meaning

1. Keep the site in place, lived in or visited: an undersea community, or an enclosed land site.
2. Keep it in place under a sealed enclosure, holding dry gas, or vacuum for a reference area.
3. Move it within the same landscape.
4. Cut out key surfaces and objects, keep them sealed, and make replicas for public access.
5. Document it fully and mark it, or let it change, where the people it matters to choose that.

### 5.5 The register at 28% water

The full register is in [heritage_register.json](results/heritage_register.json).

| Site | Situation | Treatment (evaluation status) |
|---|---|---|
| Tranquility Base (Apollo 11) | 467 m under, 8.7 atm; land 35 km | Undersea community (decided) |
| Apollo 15, Fallen Astronaut memorial | 453 m under; Apennine coast 8 km | Undersea community, memorial kept in place (proposed) |
| Apollo 17 | 1,037 m under, 17.8 atm; Taurus coast 6 km | Undersea community (proposed) |
| Lunokhod 1 / Lunokhod 2 | 865 / 1,141 m under | Undersea community / community or relic chamber (proposed) |
| Chang'e 3 / Chang'e 5 | 1,211 / 747 m under | Undersea communities (proposed) |
| Chang'e 4 / Chang'e 6 | 4,027 / 3,312 m under, 66 / 54 atm | Deep undersea communities; depth weighed against place (open) |
| Surveyor 1, Luna 16, Luna 24, Blue Ghost 1 | 378 to 2,233 m under | Community or relic chamber; scale open |
| Luna 2 | 480 m under | Marker over the impact point (open) |
| Apollo 12 and Surveyor 3 | 50 m above the sea, water 11 km | Coastal site built to become undersea if the shore rises |
| Apollo 14 / Apollo 16 | 436 / 1,416 m above | Land sites with enclosed cores; reference-area candidates |
| Luna 9, SLIM, Chandrayaan-3, IM-1 | 394 to 4,480 m above | Land sites with enclosed cores |
| Shoemaker's ashes | under the Shoemaker crater lake | Marked scattering site; impact zone left out of excavation; family and Diné consulted |
| LCROSS, Beresheet, Hakuto-R, Luna 25, Resilience | impact and crash sites | Documented and marked; relics and the Beresheet archive recovered where found |

### 5.6 Undersea Tranquility community

These figures are from [gates.json](results/gates.json).

- **Loads.**
  - The dome pushes outward at +1 atm in vacuum and is about neutral under the new air.
  - At 467 m of water the net inward load is 7.7 atm.
  - The design meets the final state, and anchoring meets the first.
- **Shell.** A dome of 100 m radius needs about 0.40 m of steel or 1.95 m of lunar concrete, sized conservatively for buckling.
- **Ground.**
  - Saturated ground pushes up on a 1-atm floor with the weight of 178 m of basalt, which calls for rock anchors or a sealed cut-off wall with drainage.
  - The preserved surface sits in sealed trays, isolated from groundwater.
- **Design depth.** About 700 m, covering the shoreline's uncertainty. At that depth the net load is 11.4 atm and the shell 0.49 m of steel.
- **Access tower.** It stresses its structure about as much as a 77 m tower on Earth. A harbour and landing mark the spot on the Sea of Tranquility.
- **Diving.** 62.4 m of lunar water adds 1 atm, so the site sits at the pressure of a 70 m dive on Earth.

### 5.7 Sealed reference area (proposed)

- **What it is.** A vacuum enclosure a few hundred metres across on high, dry ground. It holds a young crater, lowland and highland regolith, and a sample of polar ice kept at native cold.
- **Candidates.**
  - Apollo 16, with North Ray and South Ray craters, 1,416 m above the sea;
  - Apollo 14, with Cone crater, 436 m above the sea.

  Apollo samples already come from all three craters.
- **Engineering.**
  - The air alone loads the enclosure with 74.9 t/m², the weight of 28 m of rock, which needs about 0.24 m of steel at 150 m radius.
  - Holding polar ice at 40 K takes about 42 MW per km².
- **Upkeep.** Its failure mode, collapse, is the most severe of the methods, so it receives the most dependable care.

## 6. Scientific preservation

### 6.1 Zones and targets

[science_targets.csv](results/science_targets.csv) places 31 targets.

- **Zone A, the future seafloor.** This is 28% of the Moon, everything below −1,654 m, sampled and drilled before flooding. It includes the lava stacks of the near-side maria and their ancient soils sealed between flows, the young basalts around Mons Rümker, and the swirl at Reiner Gamma. It also includes the volcanic glass of Sulpicius Gallus and Aristarchus, the South Pole–Aitken floor at Von Kármán, Apollo, Poincaré and Schrödinger, and the dark floor of Tsiolkovskiy, about 110 m below sea level, which becomes a crater lake.
- **Zone B, the polar cold traps.** These are sampled before the first gas arrives. Most become crater lakes.
  - Near the pole, the 4 and 16 px/deg cylindrical grids disagree by hundreds of metres for small craters such as de Gerlache, so the polar programme adopts LOLA's polar stereographic grids.
  - Modelling shows lander exhaust already reaching the cold traps, so this zone is urgent now.
- **Zone C, land.** These are sampled before the ground is wetted: Ina, the young craters at Apollo 14 and 16, and the highlands generally.

### 6.2 Programme, most urgent first

1. Clean excavation and coring of the cold traps at both poles. Cores are stored at native temperature (a sunshaded vault in space holds about 40 K passively), and some traps are set aside from landing and mining until sampled.
2. Survey and sampling of every geological unit, before the dusty stage and before water.
3. Deep drilling through the lowland lava stacks, before groundwater reaches them.
4. A search for rock from Earth, Mars and Venus in the regolith, and for gases implanted from Earth's atmosphere, before weathering and before life arrives.
5. Geophysics under the old conditions: seismic, heat flow, electrical and magnetic sounding, gravity, exosphere, dust and plasma.
6. Complete surface records, and the Moon's photometry as seen from Earth.
7. Documentation of heritage sites. Where the heritage process agrees, this includes studies of long-exposed hardware and of the Earth microbes left at landing sites.
8. Research through the transformation: monitoring networks, sediment cores from the new seas, and the crust's response to loading.

**Storage.**

- Samples stay in original condition, in sealed vacuum or dry inert gas, and cryogenic for ice.
- Part of each category stays unopened for future methods.
- Copies are held in more than one place.
- Provenance is recorded and access is open.

### 6.3 Gates

Each irreversible step waits until the record it destroys has been taken.

| Stage | What ends | What must be complete before it |
|---|---|---|
| Now: landers and prospecting | Clean polar volatiles | First clean polar cores; untouched cold traps set aside |
| First gas, days into bulk delivery | The exosphere; solar wind at the ground; nitrogen frosting on 40 K ground above 6 Pa, oxygen above 0.2 Pa | Polar coring; exosphere, dust and plasma baselines |
| About 0.1 atm, year 42 of a linear 500-year build-up | Footprints, tracks and surface fabric. Bare regolith moves in 10 m winds of 14 m/s, falling to 8 m/s at 0.3 atm and 4 m/s at full pressure (Earth sand: about 7 m/s) | Heritage enclosures; the surface sampling grid |
| Warm poles | Cold traps; cold ground | Deep polar cores |
| Humidity, rain and groundwater | Coatings, metal and grain rims; the dry-rock seismic character; the crust to kilometres deep | Deep drilling; the search for rock from Earth; a sealed regolith archive; geophysical records |
| Seas and clouds | Lowland sites; the old reflectors; the calibration standard | Site decisions; ranging stations tied in; a calibration replacement |
| Settlement and life | Radio quiet; the organic baseline from before life | A radio replacement; sterile sealed samples |

**Erosion under rain**, relative to Earth at the same discharge and slope:

| Measure | Lunar value relative to Earth |
|---|---|
| Raindrop impact energy per kilogram of rain | 0.07–0.14 |
| Shields number | 1.8 times |
| Bedload far above threshold | about 1.0 |
| Settling speed of fine sediment | 0.17 |
| Ease of suspension | 3.3 times |
| Stable height of cohesive slopes | 6 times |

Fresh regolith behaves like fresh volcanic ash, which set world sediment-yield records for years after the 1991 Pinatubo eruption. Soil crusts and vegetation early in the wet stage shorten that period.

## 7. Replacing lost functions

| Function | Why it ends | Replacement |
|---|---|---|
| Calibration standard for Earth-observing satellites (reflectance stable to about 10⁻⁸ per year) | Clouds, water, vegetation | Traceable reference satellites and engineered targets. Existing efforts: the Landolt artificial star (NASA with NIST), planned for launch from 2028; ESA's TRUTHS, suspended at the end of 2025 |
| Far-side radio quiet | A new ionosphere; settlement noise | A radio array in deep space behind an engineered shield |
| Laser-ranging targets | Submerged sites; about 17 m of zenith delay in the new air, and clouds | New stations on high ground tied to the old reflectors; radio transponders; local weather instruments |
| Low lunar orbit | Air to thousands of kilometres | Gravity mapping and close imaging completed beforehand; high orbits and surface gravimetry afterwards |
| Seismic quiet, dry-rock physics, the natural exosphere | Ocean, wind, water | Recorded beforehand (section 6) |

## 8. Stewardship and decisions

- **Sites.** Any community living with a site decides together with the communities to whom the site matters. These include the Moon's nations as they form, those who hold the Moon sacred, and, for remains and memorials, families and those who knew the people commemorated.
- **The archive.** An archive organization with open review, and access for anyone who can use it.
- **Continuity.** Organizations that train their own staff, with more than one organization keeping each site, the archive and every life-support function. Services that inhabitants can build, repair and run themselves keep the world free of chokepoints.
- **Earth.** Everyone who sees the Moon shares in its change. Changes to moonlight on Earth, and to cultural practices tied to the Moon, call for study and consultation.
  - Moonlight cues that depend on brightness would shift: Arctic zooplankton migration, the activity of nocturnal mammals, and lunar clocks in marine worms.
  - Cues that depend on timing would hold, such as coral spawning in the dark interval after sunset.

## 9. Candidate community settings

[community_settings.json](results/community_settings.json) lists settings from geography alone. People choose community sites, with climate, ecology and the wishes of those who would live there. The settings are:

- **Undersea heritage communities**, from the register.
- **Related coastal communities** on the Apennine coast by Hadley (8 km), the Taurus massifs by Taurus–Littrow (6 km) and the Le Monnier coast (15 km).
- **Islands**, including Caucasus Island (68,000 km², summit 4.2 km) and Jura Island (58,000 km²), which encloses the flooded bay of Sinus Iridum.
- **Near-side coasts** with Earth in the sky.
- **Far-side waters** of the South Pole–Aitken Sea, Moscoviense and Orientale.
- **Dry polar highs** such as Nobile and Peary, bases for the polar sampling programme.

## 10. Open questions

- **The water inventory and shoreline for 28% cover.** Groundwater uptake by the porous crust and crustal loading decide how much water to deliver and where the shore falls.
- **Which closed depressions hold lakes.** This follows from each catchment's rainfall and evaporation in a climate run at 28%.
- **The form of imported oxygen**, and the use of its hydrogen.
- **Confirmation of the register**, including the deep far-side communities, the Apollo 12 shoreline design, and consultation on Shoemaker's ashes.
- **The sealed reference area**, at Apollo 16 or Apollo 14.
- **The archive's storage locations**, including a vault in space.
- **Replacement systems** and their timing.
- **Effects on Earth**: moonlight, nocturnal ecology, astronomy, and sighting of the crescent Moon.
- **Consultation with people on Earth** about a change everyone will see.
- **The rewording of principle 7.**
