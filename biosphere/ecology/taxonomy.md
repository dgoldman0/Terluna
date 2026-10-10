# A taxonomic system for the Open Moon's biosphere

On 9 October 2026 the author asked to rename the aerophytes and to begin a taxonomic system: "let's also change the
name. Indeed, let's start working out a taxonomic system for this world's biosphere"
([decisions](../../research/decisions.md)). This document proposes the system, sets out the choices it leaves to the
author, registers the biosphere's designed groups as the project has worked them out, and offers names for the
aerophytes. [taxa.py](taxa.py) holds the register as code and writes [results/taxa.json](results/taxa.json) (schema
`terluna.biosphere.ecology-taxa/1`); [taxonomy_sources.json](taxonomy_sources.json) records what was read and how
far.

**Everything in the register is design.** No organism, group or trait named here has been built, tested or
released. Each entry carries the project's own verdict on it, taken from the files it cites.

## The system in brief

- **Descent.** A designed organism keeps the Earth lineage of its chassis: the lineage whose genome organizes the
  cell and is inherited. Genes, pathways and symbionts from other lineages are recorded as contributions. Until a
  chassis is chosen, the register places a group at the smallest clade of an Earth backbone that holds every
  candidate.
- **Degree of design.** A design is a line within an Earth species, a designed species in an Earth genus, a designed
  genus in an Earth family, or a new lineage that no Earth genus holds. Each degree takes a naming form Earth already
  uses for engineered organisms, cultivars or new taxa.
- **The founding register.** Every lineage released on the Moon is registered once, with a permanent identifier, its
  designed genome deposited as the type, a preserved sample, its purpose, the organisms it was built from, its release
  and every later version.
- **Descendants.** Over the billion-year horizon the founders' descendants are named by descent. A clade is an
  ancestor and all its descendants, so "the clade originating in founding lineage X" keeps its meaning however far
  they diverge. Descendant species receive binomials when they become distinct.
- **Function.** Life form, setting, guild, night strategy and clock are classed apart from descent, after Raunkiær,
  by how an organism lives through the 354-hour night and where. The aerial forms come first: aeroplankton,
  drifters, sailors, giants, flyers, gliders, films and the tenants of aerophytes.
- **Names.** Scientific names are treated as Latin and drawn from any language, the Moon's own names included, with
  the etymology stated, no homonym of any Earth genus and nothing derogatory to any people. Common names belong to the
  lunar peoples and are recorded as they arise.
- **Identifiers.** Every entry has a permanent identifier, and names can change around it. The aerophytes kept G01
  when they were renamed.

The first register holds 36 designed groups with 11 subgroups, 28 design modules and 12 communities, hung from a
backbone of 25 Earth clades. Six candidate scientific names for the aerophytes wait for the author's choice.

## What Earth's systems name, and what each suits here

| System | What it names | What fixes a name | What it suits on the Open Moon |
|---|---|---|---|
| ICN, Madrid Code 2025 | Algae, fungi and plants, kingdom to form | A preserved specimen or illustration; metabolically inactive cultures of algae and fungi (Art. 8.1, 8.4) | Plant, algal and fungal chassis; the Latin grammar of names |
| ICZN, 4th edition 1999 | Animals: family, genus and species groups | A name-bearing type; ZooBank registration for electronic works (Art. 8.5) | Animal chassis; names from any language (Art. 11.3) |
| ICNP, 2022 revision | Prokaryotes, phylum to subspecies, with domain and kingdom added by 2024 | A type strain living in two collections in different countries (Rule 30) | Bacterial chassis; strains below the Rules (Rule 5d) |
| SeqCode 1.2.0, 2026 | Prokaryotes described from sequences | A genome sequence in INSDC (Rule 18a) and registration (Rule 26) | The model for genome types and registration |
| ICNCP, 9th edition 2016 | Cultivars, Groups and grexes | Characters that stay distinct, uniform and stable (Art. 2.3) | Lines on farms; names in living languages |
| PhyloCode 6, 2020 | Clades, without ranks | A phylogenetic definition with specifiers, registered (Art. 8.1, 9) | The founders' descendants over a billion years |
| Draft BioCode 2011 | All organisms | Types and registration (Art. 5.2) | The precedent for one lunar code |
| Raunkiær 1934 | Plant life forms | Where the buds that survive the unfavourable season sit | Forms classed by how they carry the night |
| Sieburth et al. 1978 | Plankton size classes | Decades of size | Aeroplankton size classes |
| WMO Cloud Atlas 2017 | Clouds | Genera, species and varieties in Latin | The sky's other Latin classification |

**No Earth code can validly publish a name for an organism that does not yet exist.** The ICZN excludes names
proposed for hypothetical concepts (Art. 1.3.1), which its glossary defines as concepts that contained no animal then
known to exist in nature. The ICN, ICNP and SeqCode tie each name to a type that exists: a preserved specimen, a
living type strain or a deposited genome. The names in this register are therefore proposals in the codes' forms,
ready for valid publication once the organisms and their types exist.

**The codes agree on their foundations and differ in detail.** The ICN's six principles are independence, types,
priority, a single correct name, names treated as Latin whatever their derivation, and retroactivity. The ICZN starts
on 1 January 1758 with Linnaeus and Clerck (Art. 3); the ICN starts seed plants and ferns with Linnaeus's Species
plantarum of 1 May 1753 (Art. 13.1). Since 2012 a new name under the ICN needs a
Latin or English description or diagnosis (Art. 39.2). The ICN declares itself independent of zoological and
prokaryotic nomenclature, while the ICNP and the SeqCode avoid names already regulated by the others (ICNP Principle 2; SeqCode
Principle 2 and Rule 9b). Cyanobacteria fall under the ICN and the ICNP at once (ICN Preamble 8), and their phylum
name, Cyanobacteriota, rests on a type genus validated under both (Oren et al. 2022).

**Earth names engineered organisms below the species, each with an identifier.** The chemically synthesized
*Mycoplasma mycoides* JCVI-syn1.0 genome, booted in a *M. capricolum* cell, made new *M. mycoides* cells, and its
sequence carries watermarks (Gibson et al. 2010): the species followed the genome. The minimal JCVI-syn3.0, 531
kilobase pairs and 473 genes, is still a strain of *M. mycoides* (Hutchison et al. 2016). Sc2.0 redesigns the
*Saccharomyces cerevisiae* genome with design software that enforces version control of every edit (Richardson et
al. 2017); its chromosomes are named synII to synXVI (Goold et al. 2025), all 16 now exist in functional synthetic
versions (Archer et al. 2026), and the organism is still *S. cerevisiae*. *Escherichia coli* with a 61-codon synthetic
genome is a variant of *E. coli* (Fredens et al. 2019). A crop's engineered insertion is a transformation event:
*Zea mays* event MON810 carries the unique identifier MON-ØØ81Ø-6 and the trade name YieldGard (Biosafety
Clearing-House). The identifier has nine characters in three parts, an applicant code, an event code and a check
digit, and a stack of events receives a new one (Commission Regulation (EC) No 65/2004, following the OECD). An
assemblage of genetically modified plants may form a cultivar (ICNCP Art. 2.19), though such lines are usually sold
under trademarks.

**A designed lineage's genome is its natural type, and a register is its natural place of publication.** The SeqCode
types prokaryote species by a genome deposited in INSDC and publishes a name only when it is registered. The
International Committee on Systematics of Prokaryotes rejected sequences as types in 2020, and the Madrid congress
rejected DNA-sequence types under the ICN and sent the question to a committee reporting in 2029 (Turland 2025). Registration already fixes names elsewhere: ZooBank for
animals in electronic works, nomenclatural repositories under the ICN (Art. 42.1), the PhyloCode's registration
database (Art. 8.1), and the draft BioCode's registration as the condition of establishing any name (Art. 5.2).

**Phylogenetic nomenclature names clades by definition, independent of rank, and the founders are known ancestors.**
Under the PhyloCode a clade is an ancestor and all its descendants (Art. 2.1). A name is established with a
definition in English or Latin (Art. 9.3): the smallest clade containing A and B, the largest containing A but not Z,
the clade characterized by an apomorphy as inherited by A, or the clade originating in a directly specified ancestor
A (Note 9.5.1). Total clades take the prefix Pan- (Art. 10.3). Hybridization and symbiogenesis can make clades
overlap, and the definitions allow it (Note 2.1.3). Species names stay with the rank-based codes (Art. 21.1). On the
Open Moon every founder is a registered ancestor, which suits the directly-specified-ancestor definition exactly.

**Life-form systems class organisms by how they meet the unfavourable season.** Raunkiær classed plants by where the
buds that survive it sit: phanerophytes carry them on stems in the air, chamaephytes near the ground,
hemicryptophytes at the soil surface, cryptophytes buried in the ground or in water (most bulbous and tuberous
plants), and therophytes complete their lives in one favourable season and wait as seed. He divided phanerophytes by
height: nano under 2 m, micro 2–8 m, meso 8–30 m, mega over 30 m (Raunkiaer 1934). Plankton is classed by decades of
size: pico 0.2–2 µm, nano 2–20 µm, micro 20–200 µm (Sieburth et al. 1978). Aeroecology studies the organisms that
depend on the aerosphere and move through it by passive and active displacement (Kunz et al. 2008), and current
work calls the wind-carried organisms, living nematodes, mites and thrips among them, aeroplankton (Ptatscheck et al.
2018). On the Open Moon the unfavourable season is the 354-hour night, warm and frost-free at low and middle
latitudes, so the forms here are classed by how an organism carries it.

**The sky's other inhabitants already have Latin names.** The WMO's International Cloud Atlas classes clouds in
genera, species and varieties with Latin names, "similar to the systems used in the classification of plants or
animals". The candidate names below avoid the cloud genera.

**The codes welcome names from every language and now refuse offensive ones.** A zoological name "may be a word in or
derived from Latin, Greek or any other language (even one with no alphabet)" (ICZN Art. 11.3), and the Code's
examples include opossum from Algonquian, *Abudefduf* from Arabic, nakpo from Tibetan and canguru from Kokoimudji. The
ICN takes generic names and epithets "from any source whatever" (Art. 20.1, 23.2); the SeqCode takes roots from "any
language in use or extinct" (Rule 26); the ICNP prefers Latin or Greek where an equivalent exists (Rule 6). Since
1 January 2026 a new name under the ICN that is derogatory to a group of people can be rejected (Art. 51.2);
epithets on the root caffr- are now spelled afr- (Art. 61.6); and the ICZN's Code of Ethics asks that no author propose a name likely to
give offence on any grounds (Appendix A.4). The Moon has its own Latin root stock: the IAU names of its seas, bays,
lakes and marshes and their stated origins, such as Mare Nubium, "Sea of Clouds", Sinus Roris, "Bay of Dew", Mare
Insularum, "Sea of Islands", and Lacus Temporis, "Lake of Time". In the atlas the Near-side Sea floods 99% of Mare
Imbrium, the "Sea of Showers" ([ecology register](../../research/studies/lunar_cycle_ecology/README.md), H10).

## Principles

**P1. A designed organism keeps the Earth lineage of its chassis.** The chassis is the lineage whose genome organizes
the cell and is inherited, as JCVI-syn1.0 is *M. mycoides* by its genome. Genes, pathways and symbionts from other
lineages are contributions, recorded with their donors' Earth names: an aerophyte built on a brown alga, carrying a
green alga's hydrogenase and bacterial symbionts, is a brown alga with two contributions. Where the chassis is open, the
register places a group at the smallest backbone clade holding every candidate, and the placement narrows as the
design firms up.

**P2. How far a design departs sets its degree, and each degree has a name form Earth already uses.**

| Degree | What it is | Name form | Earth precedent | Anchor |
|---|---|---|---|---|
| line | Timing, tolerance or allocation changed within one Earth species | Earth binomial, founding identifier, optional line name: *Solanum lycopersicum* OM-L-0007 'Long Dusk' | JCVI-syn3.0, Sc2.0, event MON810, a GM cultivar | the species or above |
| designed species | Traits beyond any species of an Earth genus | A new binomial in that genus | A new species | a genus or above |
| designed genus | A design no Earth genus holds, within an Earth family | A new genus and species in the family | A new genus | a family or above |
| new lineage | A body plan or metabolism no Earth genus holds | A clade name defined from its founders, with genera and species inside it | A PhyloCode clade | a family or above |

The July 2026 bundle prefers to "modify stress tolerance and timing before inventing new metabolisms" (§18), and most
of the register's groups would found designed species or lines. The aerophytes are a new lineage.

**P3. Every founding lineage is registered with a type genome and a preserved sample.** The designed genome is
deposited as a sequence, as the SeqCode types prokaryotes, and a cryopreserved sample is kept in two collections in
different places, as the ICNP asks of type strains and the ICN accepts for inactive cultures. Registration fixes the
name's date and priority. The register also records the lineage's purpose, its reference organisms, its release and
every later version; Sc2.0's software already tracks each edit, and JCVI-syn1.0's genome carries a watermark. The
July bundle asks to "record every released genome and location" (§18).

**P4. Descendants are named by descent from the founders.** A clade defined as originating in a registered founder
keeps its meaning through a billion years of divergence, whatever ranks its descendants come to deserve. A descendant species receives a binomial when it becomes distinct, typified by a sampled genome
and specimen, and the founder's registered genome stays the reference for every clade defined from it. Gene flow
between lineages makes partly overlapping clades, which the PhyloCode's definitions accommodate.

**P5. Function is classed apart from descent.** Two lineages can share a form, and one lineage can take several
forms in its life: an aerophyte begins as a drifting propagule and may end as a giant. As Raunkiær did, the register
classes each group by its adult.

**P6. Scientific names are Latin in grammar and open in origin.** They follow the strictest Earth rule wherever the
codes differ, so a name stays publishable under any of them:
- treated as Latin, in the 26-letter alphabet, diacritics dropped (ICZN Art. 11.2; ICN Art. 60.7);
- species names binary, with no tautonyms and epithets of 2–30 characters (ICN Art. 23.1, 23.2, 23.4);
- roots from any language, the Moon's own names included, with the etymology stated (ICZN Art. 11.3; SeqCode
  Rule 26; ICNP Rule 6; ICN Rec. 60I.1);
- no homonym of any genus under any Earth code, checked against IRMNG and the GBIF backbone (SeqCode Principle 2);
- nothing derogatory to any people, and a word from a living language adopted with its speakers (ICN Art. 51.2,
  Rec. 51A.1; ICZN Appendix A.4);
- no standard prefix stamped on Earth names to mark them lunar, which the ICZN treats as formulae outside
  nomenclature (Art. 1.3.7).

**P7. Common names belong to the peoples who live with the organisms.** The Moon's nations descend from all of
Earth's nations and form as their own lunar nations, and the names they give are recorded as they arise: the name,
its language and community, its date and its meaning. Earth already keeps such names beside Latin ones, as the
ICNCP's cultivar epithets in living languages show (Art. 21.11). The author chose the everyday names aerophytes and
sky reefs on 10 October 2026; the project's other English names stay marked as working names until the peoples name
them.

**P8. Identifiers never change.** Every group, module, community and founding lineage keeps its identifier through
renaming, regrouping and new placements, as a transformation event keeps its unique identifier. The data products and
the immersive world key on identifiers.

## Choices for the author

**O1. Which rules govern new designed taxa?** (a) The Earth code of each chassis; (b) one lunar code for all designed
taxa, in forms every Earth code accepts; (c) the Earth codes alone, with every design kept as a line within an Earth
species. **Option (b) is recommended.** The Earth codes publish names only for organisms with types; chimeric designs
and groups the codes share, as they share the cyanobacteria, need one rulebook; the draft BioCode shows the form; and
taking the strictest Earth rule leaves each name publishable later under its chassis's code.

**O2. Ranks.** (a) Linnaean ranks for founders and rankless clade names for descendants; (b) rankless clade names
throughout, with binomials only for species. **Option (a) is recommended:** binomials need genera, and the descendant
clades need names that keep their meaning through deep time.

**O3. The founding identifier.** (a) A register number with a check digit, such as OM-L-0001-7; (b) three parts after
the OECD and EU format, a designer code, a lineage code and a check digit. Both carry a check digit; (b) also records
who designed the lineage.

**O4. Roots for scientific names.** (a) Latin and Greek, as the ICNP advises; (b) any language in Latin grammar, as the
ICZN, ICN and SeqCode allow; (c) a lunar root stock for founders, from the IAU names and the Moon's own phenomena, and
any language for descendants named by the lunar peoples.

**O5. Where founding ends and descent begins.** (a) Every release is a founding lineage, a re-engineered release a new
one derived from the old, and change after release is descent; (b) founding ends with the 500-year build. **Option (a)
is recommended**, since releases can continue after the build and each needs its type.

**O6. Giants and lesser aerophytes.** (a) One lineage whose oldest members reach the giant form; (b) two sister
lineages, the giants designed for growth that never stops. The coupled study gives the giants two forms, sky reefs
(colonies of wet modules on a raft) and round bodies in many compartments, and either can sit under (a) or (b). The
register holds the giants as a subgroup of the aerophytes (G03) under either answer.

**O7. The aerophytes' scientific name**, from the candidates below. The author chose the everyday names on 10 October
2026: aerophytes, and sky reefs for the colony giants.

## The functional classification

### Realms and settings

| Realm | Setting | What it is |
|---|---|---|
| sky | low sky | 0–10 km above sea level: the day's mixed layer, cloud base 4–7 km |
| sky | storm layer | 10–22 km: storms to their median top |
| sky | cool air | 22–35 km: the air freezes near 26 km |
| sky | flight band | 35–45 km: the long-haul flight band |
| sky | high air | 45–90 km |
| land | canopy | Crowns, branches and epiphyte perches of forests |
| land | ground | The ground surface, its litter, soils and plains |
| water | seas | Open seas, productive coasts, river mouths and tidal flats |
| water | lakes and wetlands | Rain-fed lakes, their shores and wet margins |
| water | rivers | Rivers from the lakes to the seas |

The sky's bands are those of [people and land](people_and_land.md) (table 1). Each group lists every setting it uses
and names one realm as its home.

### Life forms

| Class | Form | What it is |
|---|---|---|
| aerial | aeroplankton | Organisms the air carries, too small to hold a course; resident aeroplankton complete their lives aloft, transient aeroplankton are dispersing stages |
| aerial | drifter | Buoyant or ballooning organisms that move with the wind and control mainly their height |
| aerial | sailor | Buoyant organisms that steer across the wind with a wing or drogue hung in another layer of air or in the sea, as *Physalia* sails the sea surface (Iosilevskii & Weihs 2009; [aerophyte study](../../research/studies/aerophytes/sailing.md)) |
| aerial | giant | Buoyant organisms that keep growing for centuries, like the oldest trees |
| aerial | flyer | Animals that fly under their own power |
| aerial | glider | Animals that launch from canopies and cliffs and glide |
| aerial | film | Sparse microbial or algal films on particles, droplets or platforms in the high air |
| aerial | tenant | Organisms living on or within aerophytes and other aerial hosts, in their retained-water pockets |
| plant | tree | Woody plants whose living crowns stand in the air through the night (phanerophytes) |
| plant | vine | Climbing or trailing plants |
| plant | herb | Non-woody plants whose shoots stand through the night |
| plant | geophyte | Plants that carry the night in bulbs, corms or tubers below ground (cryptophytes) |
| plant | epiphyte | Plants rooted in the canopy |
| plant | crust | Crusts and mats of cyanobacteria, mosses and lichens |
| ground | fungal | Yeasts and mycelial fungi |
| ground | ground animal | Animals that walk, climb or burrow |
| microbial | microbe | Free-living microbes of soils, sediments and waters |
| aquatic | plankton | Organisms the water carries |
| aquatic | nekton | Animals that swim against currents |
| aquatic | benthos | Organisms of the beds of seas and lakes |
| aquatic | macroalga | Seaweeds |

The tenant form extends the requested list of aerial forms, for the small communities the aerophyte study places in its
organisms' retained-water pockets.

### Size

**At lunar gravity particles fall about six times slower, so an aeroplankter's size class sets how long it can stay
aloft.** The first screen's particles, placed in Sieburth's classes, settle as follows. The times are dry settling out
of a mixed layer of the stated depth; rain, grazing and death remove organisms much sooner.

| Class | Size | Example | Falls on the Moon (Earth) | Stays in a 5 km layer on the Moon | Stays in Earth's 1.5 km layer |
|---|---|---|---|---|---|
| pico | 0.2–2 µm | 1 µm bacterium | 0.0062 mm/s (0.039) | 9,364 days | 447 days |
| nano | 2–20 µm | 5 µm spore | 0.14 mm/s (0.86) | 414 days | 20 days |
| nano | 2–20 µm | 10 µm spore | 0.55 mm/s (3.4) | 105 days | 5.1 days |
| micro | 20–200 µm | 25 µm pollen grain | 3.1 mm/s (19) | 19 days | 0.9 days |

Trees take Raunkiær's height classes: nano under 2 m, micro 2–8 m, meso 8–30 m and mega over 30 m. The megaforest
studies' trees of 100–500 m are megaphanerophytes.

### Guilds, night strategies and clocks

- **Guilds:** producer, grazer, fruit-eater, seed-eater, predator, collector, consumer (diet open), decomposer,
  pollinator, fixer, soil builder, air keeper (removes methane, N2O, hydrogen, carbon monoxide or ethylene, which
  the air under the shield cannot), carrier (brings elements back from the seas to the land), receiver (takes up
  what carriers bring), lure. Adl et al. (2019) give protists a guide to trophic functional guilds.
- **Night strategies:** store, idle, dormant, twilight, night-active, terminator-following (travelling to change the
  light: with dusk or the Sun, east on the wind to meet the next dawn sooner, or across the wind toward the twilit
  poles), and shedding, which the register sets aside.
- **Clocks:** dusk-set, dawn-set, carbon-cued, two clocks (daily rhythms inside the lunar clock), and an Earth clock
  kept on daily cues on farms.
- **Carriers** name the elements they carry: phosphorus, nitrogen, sulfur, iodine, potassium and trace elements, by
  biology or through the air ([ecology register](../../research/studies/lunar_cycle_ecology/README.md), G3 and H).

"Undetermined" marks what the cited files leave open.

### The speed of dusk

Dusk and the Sun cross the Moon westward at 2π(R + h) cos(latitude) per synodic month of 708.73 hours. Dusk followers
fly it at the ground. Aloft the mean wind blows east, so the aerophyte study's aerophytes ride it to shorten their
nights, linger in the evening flow near the ground, or sail north and south across the shear
([aerophyte study](../../research/studies/aerophytes/README.md)).

| Latitude | Dusk at the ground | Sun at 10 km |
|---|---|---|
| 0° | 15.4 km/h, 4.28 m/s | 4.30 m/s |
| 30° | 13.3 km/h, 3.71 m/s | 3.73 m/s |
| 60° | 7.70 km/h, 2.14 m/s | 2.15 m/s |
| 80° | 2.67 km/h, 0.74 m/s | 0.75 m/s |

## The founding register

The founding register holds one record per released lineage:
- its identifier, name, degree and group;
- its chassis, with the Earth name and authority, and its contributions, with their donors' Earth names;
- the modules it carries and the purpose it was designed for;
- its type genome (a sequence accession) and its preserved sample (two collections in different places);
- its genome's watermark;
- its release: date, place, scale and containment;
- its versions, each edit with a link to the lineage it derives from;
- its status, and the common names given to it (name, language, community, date and meaning).

No lineage has been designed or released, so the register's list of lineages is empty and its fields are fixed.

## The first register

The register draws on the [ecology register](../../research/studies/lunar_cycle_ecology/README.md) (its codes A1–H10
are cited), this branch's [ecology findings](README.md), [literature](literature.md) and
[people and land](people_and_land.md), the [aerial](../../research/studies/aerial_ecology/README.md),
[aerophyte](../../research/studies/aerophytes/README.md),
[megaforest](../../research/studies/megaforest_wind/README.md),
[forest-patch](../../research/studies/forest_patch/README.md) and
[sea-appearance](../../research/studies/sea_appearance/README.md) studies, the
[living-water](../living_water/waters.json) guesses and the [decisions register](../../research/decisions.md).
Statuses: carried forward, test first, hypothesis and set aside, with "in places" where the selection limits a design
to particular landscapes or to farms. The placement column gives the backbone anchor, the degree of design and the
life form; subgroups are marked ↳.

### The sky

| ID | Group | Placement: anchor · degree · form | Earth reference points | Guilds · settings | Night · clock | Status | Source |
|---|---|---|---|---|---|---|---|
| G01 | aerophytes | Eukaryota · lineage · drifter (or sailor) | *Physalia physalis*; *Nereocystis luetkeana*; holopelagic *Sargassum*; siphonophores; teleost swim bladders; mistletoes; *Chlamydomonas reinhardtii*; *Clostridium* | producer, receiver · low sky, storm layer, canopy | twilight, terminator following, store, idle · undetermined | hypothesis | [aerophyte study](../../research/studies/aerophytes/README.md) What a giant looks like; [aerophyte biology](../../research/studies/aerophytes/biology.md); [aerophyte sailing](../../research/studies/aerophytes/sailing.md) Riding the eastward wind; [aerophyte gas biology](../../research/studies/aerophytes/gas_biology.md) Siphonophores; [aerial study](../../research/studies/aerial_ecology/README.md); [decisions](../../research/decisions.md) Cross-domain ecology follow-up |
| G02 | ↳ lesser aerophytes | Eukaryota · lineage · drifter (or sailor) | as its parent | producer · low sky, canopy | twilight, terminator following, store, idle · undetermined | hypothesis | [aerophyte study](../../research/studies/aerophytes/README.md) The reference organism; [aerophyte growth](../../research/studies/aerophytes/growth.md) The canopy nursery |
| G03 | ↳ giant aerophytes | Eukaryota · lineage · giant | *Pinus longaeva*; trees of 403 species; siphonophores | producer, receiver · low sky, storm layer | twilight, terminator following, store, idle · undetermined | hypothesis | [decisions](../../research/decisions.md); [aerophyte study](../../research/studies/aerophytes/README.md) What a giant looks like; [aerophyte gas biology](../../research/studies/aerophytes/gas_biology.md) |
| G04 | resident aeroplankton phototrophs | Life · species · aeroplankton, pico/nano | cloud cyanobacteria and algae | producer · low sky, storm layer | undetermined · undetermined | hypothesis | [aerial study](../../research/studies/aerial_ecology/README.md) What the evidence supports; [aerial sources](../../research/studies/aerial_ecology/sources.json); [decisions](../../research/decisions.md) Cross-domain ecology follow-up |
| G05 | ↳ aerial cyanobacteria | Cyanobacteriota · species · aeroplankton, pico/nano | as its parent | producer · low sky, storm layer | undetermined · undetermined | hypothesis | [aerial study](../../research/studies/aerial_ecology/README.md) |
| G06 | ↳ aerial microalgae | Eukaryota · species · aeroplankton, nano | as its parent | producer · low sky, storm layer | undetermined · undetermined | hypothesis | [aerial study](../../research/studies/aerial_ecology/README.md) |
| G07 | resident aeroplankton heterotrophs | Bacteria · species · aeroplankton, pico | cloud bacteria; fog heterotrophs | decomposer · low sky, storm layer | undetermined · undetermined | hypothesis | [aerial study](../../research/studies/aerial_ecology/README.md) What the evidence supports; [aerial sources](../../research/studies/aerial_ecology/sources.json); [people and land](people_and_land.md) item 3 |
| G08 | aerial collectors | Araneae · species · tenant | orb-weaving spiders; ballooning spiders | collector · low sky, canopy | undetermined · undetermined | hypothesis | [aerial study](../../research/studies/aerial_ecology/README.md); [aerial sources](../../research/studies/aerial_ecology/sources.json) |
| G09 | large aerial filter feeders | Animalia · lineage · flyer | none named | collector · low sky, storm layer | undetermined · undetermined | hypothesis | [aerial study](../../research/studies/aerial_ecology/README.md) |
| G10 | high-altitude films | Life · species · film | *Deinococcus radiodurans*; engineered algal films | producer · flight band, high air | dormant · undetermined | carried forward | [decisions](../../research/decisions.md) Life and people; [people and land](people_and_land.md) item 2; [July bundle](../../archive/july_2026_knowledge_bundle/terluna_moon_terraforming_bundle/06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md) section 12 |
| G11 | great flyers | Aves · genus · flyer | *Diomedea exulans*; *Argentavis magnificens*; *Quetzalcoatlus northropi* | consumer · low sky, storm layer | undetermined · undetermined | hypothesis | [people and land](people_and_land.md) item 5, item 22, Table 2; [decisions](../../research/decisions.md) Life and people; [July bundle](../../archive/july_2026_knowledge_bundle/terluna_moon_terraforming_bundle/06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md) section 11 |
| G12 | gliders | Animalia · species · glider | arboreal vertebrate gliders; arboreal ants | consumer · canopy, low sky | undetermined · undetermined | carried forward | [decisions](../../research/decisions.md) Life and people; [July bundle](../../archive/july_2026_knowledge_bundle/terluna_moon_terraforming_bundle/06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md) section 11 |
| G13 | sea-forest migrant insects | Insecta · species · flyer | *Agrotis infusa*; Chironomidae; *Clunio marinus*; high-flying insect migrants | carrier · seas, low sky, canopy | night active · dusk set, dawn set | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) H1-H7 |
| G14 | seabird-like colonial fliers | Aves · species · flyer | seabird colonies; Spitsbergen colonies | carrier · seas, low sky, ground, canopy | undetermined · undetermined | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) H1, H2, H5 |
| G15 | dusk followers | Animalia · species · flyer | small birds; moths | fruit eater · low sky, canopy | terminator following · dusk set | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) E6, H3 |

### The land

| ID | Group | Placement: anchor · degree · form | Earth reference points | Guilds · settings | Night · clock | Status | Source |
|---|---|---|---|---|---|---|---|
| G16 | evergreen idlers | Tracheophyta · species · tree | *Arabidopsis thaliana*; evergreen trees under polar winters | producer · canopy, ground | store, idle, twilight · carbon cued, dusk set, two clocks | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) B1, B2, C3, D1, F2; [biosphere](../README.md) Canopy photosynthesis; [literature](literature.md) The long day and night |
| G17 | day-fruit plants | Angiospermae · species · vine or tree | Cucurbitaceae | producer · ground, canopy | store, idle · dawn set, carbon cued | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) D2-D4; [biosphere](../README.md) Canopy photosynthesis |
| G18 | ↳ day-fruit crops | Cucurbitaceae · line · vine | *Citrullus lanatus*; *Cucurbita maxima*; *Cucumis sativus* | producer · ground | store, idle · dawn set, carbon cued | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) D3; [biosphere](../README.md) Canopy photosynthesis |
| G19 | ↳ day-fruit trees | Angiospermae · species · tree | as its parent | producer · canopy | store, idle, twilight · dawn set, dusk set | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) E1, F2 |
| G20 | farm crops on Earth programs | Angiospermae · line · vine or herb | see its subgroups | producer · ground | undetermined · earth clock with cues | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) C1, D5 |
| G21 | ↳ fruit picked young | Cucurbitaceae · line · vine | *Cucurbita pepo*; *Cucumis sativus* | producer · ground | undetermined · undetermined | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) D5 |
| G22 | ↳ continuous-light tomatoes | *Solanum lycopersicum* · line · herb | *Solanum lycopersicum* | producer · ground | undetermined · earth clock with cues | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) C1 |
| G23 | storage-organ plants | Angiospermae · species · geophyte | Cape geophytes | producer · ground | store, idle · carbon cued | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) F1 |
| G24 | dusk-seeding grasses and forbs | Angiospermae · species · herb | perennial grasses and forbs | producer · ground | store, idle · dusk set | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) F1 |
| G25 | C4 plants of dry ground | Angiospermae · species · herb | C4 plants | producer · ground | store, idle · undetermined | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) A3 |
| G26 | trap plants | Caryophyllales · species · epiphyte (or herb) | *Nepenthes mirabilis*; *Drosera* | producer, predator, receiver · canopy, lakes and wetlands, ground | store, idle · carbon cued | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) E4 |
| G27 | droppings-collecting plants | *Nepenthes* · species · epiphyte | *Nepenthes lowii*; *Nepenthes hemsleyana* | producer, receiver · canopy | store, idle · carbon cued | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) E5, H5 |
| G28 | glowing lures | Eukaryota · species · herb or epiphyte or fungal | *Neonothopanus gardneri*; *Nicotiana tabacum*; glowworm larvae | lure, receiver · canopy, ground | night active · dusk set | test first | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) E3 |
| G29 | megaforest emergents | Tracheophyta · species · tree, mega | *Sequoia sempervirens* | producer · canopy, ground | store, idle · carbon cued | hypothesis | [megaforest](../../research/studies/megaforest_wind/README.md) Static plant mechanics; [forest patch](../../research/studies/forest_patch/README.md); [biosphere](../README.md) Remaining biological work; [state](state.md) |
| G30 | basalt pioneers | Bacteria · species · crust | Icelandic basalt pioneers; cyanobacteria; biological soil crusts | fixer, soil builder, producer · ground | idle, dormant · undetermined | hypothesis | [ecology](README.md) finding 13; [literature](literature.md) Ecosystems on fresh basalt |
| G31 | crust mosses | Bryophyta · species · crust | Hekla moss mats | producer, soil builder · ground | idle · undetermined | hypothesis | [ecology](README.md) finding 13; [literature](literature.md) Ecosystems on fresh basalt |
| G32 | night decomposers of the fruit fall | Fungi · species · fungal | fruit yeasts | decomposer · ground | night active · undetermined | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) E1, F2 |
| G33 | night insects of the fruit fall | Insecta · species · flyer | fruit flies | fruit eater, decomposer · ground, canopy | night active · dusk set | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) E1, F2 |
| G34 | fruit-eaters | Animalia · species · ground animal or flyer | none named | fruit eater · canopy, ground | night active · undetermined | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) F2 |
| G35 | grazers and seed-eaters of the plains | Animalia · species · ground animal | none named | grazer, seed eater · ground | undetermined · undetermined | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) F1 |
| G36 | pollinators of the violet sky | Insecta · species · flyer | *Apis mellifera*; insects that orient by polarized moonlight | pollinator · canopy, ground, low sky | undetermined · undetermined | hypothesis | [ecology](README.md) finding 2; [literature](literature.md) Light; [July bundle](../../archive/july_2026_knowledge_bundle/terluna_moon_terraforming_bundle/06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md) section 11 |
| G37 | lake insects | Diptera · species · flyer | Chironomidae | carrier · lakes and wetlands, ground | undetermined · undetermined | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) F3, H2 |

### The waters

| ID | Group | Placement: anchor · degree · form | Earth reference points | Guilds · settings | Night · clock | Status | Source |
|---|---|---|---|---|---|---|---|
| G38 | sea phytoplankton | Life · species · plankton, pico/nano/micro | Earth's ocean phytoplankton | producer · seas | undetermined · undetermined | carried forward | [decisions](../../research/decisions.md) Life and people; [living water](../living_water/waters.json); [sea appearance](../../research/studies/sea_appearance/README.md) The waters' colour; [ecology](README.md) finding 17 |
| G39 | dark-night luminous plankton | Life · species · plankton | luminous marine organisms | undetermined · seas | night active · undetermined | hypothesis | [decisions](../../research/decisions.md) Life and people; [sea appearance](../../research/studies/sea_appearance/README.md) Decisions |
| G40 | sulfur- and iodine-releasing seaweeds | Eukaryota · species · macroalga | Earth's sea plankton and seaweeds | producer, carrier · seas | undetermined · undetermined | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) H8; [ecology](README.md) finding 8 |
| G41 | sea grazers | Eukaryota · species · plankton | none named | grazer · seas | night active · undetermined | hypothesis | [living water](../living_water/waters.json) |
| G42 | run fish | Salmonidae · species · nekton | salmon | carrier · seas, rivers | undetermined · undetermined | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) H1, H2, H5 |

### Microbes of soils and waters

| ID | Group | Placement: anchor · degree · form | Earth reference points | Guilds · settings | Night · clock | Status | Source |
|---|---|---|---|---|---|---|---|
| G43 | keepers of the air | Bacteria · species · microbe | see its subgroups | air keeper · ground, lakes and wetlands, seas | undetermined · undetermined | hypothesis | [ecology](README.md) findings 6-9; [state](state.md) |
| G44 | ↳ methane oxidizers | Bacteria · species · microbe | Earth's soils | air keeper · ground, lakes and wetlands, seas | undetermined · undetermined | hypothesis | [ecology](README.md) finding 6; [literature](literature.md) Trace gases without OH |
| G45 | ↳ complete denitrifiers | Bacteria · species · microbe | as its parent | air keeper · ground, lakes and wetlands, seas | undetermined · undetermined | hypothesis | [ecology](README.md) finding 7 |
| G46 | ↳ trace-gas consumers of soils | Bacteria · species · microbe | soil bacteria of 51 phyla; soils | air keeper · ground | undetermined · undetermined | hypothesis | [literature](literature.md) Trace gases without OH; [people and land](people_and_land.md) item 16 |
| G47 | nitrogen fixers | Bacteria · species · microbe | free-living fixers | fixer · ground, seas | undetermined · undetermined | hypothesis | [ecology](README.md) finding 12; [literature](literature.md) Nitrogen |

### Design modules

| ID | Module | What it does | Carried by | Status | Source |
|---|---|---|---|---|---|
| M01 | carbon-balance idling | Leaves, stems and roots idle at about a quarter of their daytime upkeep as soon as photosynthesis stops covering upkeep, cued by the energy sensors SnRK1 and TOR; this halves the night store against idling only in full darkness (10 against 19 g C per m2) | G16, G17, G23, G26 | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) B1, C3; [biosphere](../README.md) Canopy photosynthesis |
| M02 | dusk-set night timer | A timer sized for about 240 dark hours, started at dusk, that rations the night's store; a plant on an Earth clock would budget the lunar night's store for one day | G16, G17 | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) C2; [ecology](README.md) finding 20 |
| M03 | two clocks | Daily rhythms in leaf chemistry inside the long lunar clock, with the daily clock kept off light harvesting | G16 | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) C1, C4; [July bundle](../../archive/july_2026_knowledge_bundle/terluna_moon_terraforming_bundle/06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md) section 2 |
| M04 | long-dark leaf survival | Leaves that survive about 240 dark hours without senescing, with metabolism that slows and restarts on cue | G16 | test first | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) B5; [biosphere](../README.md) Canopy photosynthesis |
| M05 | twilight use | Photosynthesis through the long twilight, worth 14% of an Earth-like night store and a quarter of the idling one; it leaves 238 dark hours at the equator and 34 at 60 degrees | G16 | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) B2 |
| M06 | day-fruit program | A fruit program one sunlit half long: set at sunrise on a dawn cue without pollination, cells formed before sunrise, sugar loaded in the last days before dusk, ripe at sunset | G17 | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) D2-D4; [biosphere](../README.md) Canopy photosynthesis |
| M07 | lunar calendar of animal lives | Emergence, spawning and migration timed by dawn and dusk, the month's sharpest cues: migrants leave the sea in the lunar afternoon, arrive at dusk and lay eggs that hatch by dawn | G13 | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) C5, H7 |
| M08 | continuous-light tolerance on farms | Cultivars that tolerate continuous light, and daily cues such as a temperature swing or alternating red and blue light | G22 | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) C1; [ecology](README.md) finding 20 |
| M09 | lasting glow | Bioluminescence from the fungal pathway that lasts through the first nights after dusk; a glowing square metre draws about 0.04 W, under 2% of the night upkeep of a square metre of plants | G28 | test first | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) E3 |
| M10 | violet and visible signals | Flowers that signal at 350-400 nm or in the violet and visible, and pollinators tuned to them: a honeybee's UV receptor catches 7% of Earth's light under the shield | G36 | hypothesis | [ecology](README.md) finding 2; [literature](literature.md) Light |
| M11 | raised night-light thresholds | Circadian and melatonin thresholds raised for nearside nights of 2-3 lux of Earthlight, set zone by zone; far-side low latitudes keep Earth-calibrated thresholds | G13, G33, G34, G36 | hypothesis | [ecology](README.md) finding 5; [literature](literature.md) Light |
| M12 | vitamin D without UV-B | A route to vitamin D for every vertebrate without UV-B, which the shield removes: biosynthesis that runs without the light, calcium handling without the vitamin, or a dietary supply | G11, G14, G42 | hypothesis | [ecology](README.md) finding 1; [literature](literature.md) Light |
| M13 | complete denitrification | Denitrification that ends in N2 through the N2O reductase NosZ, since N2O has no sink under the shield | G45 | hypothesis | [ecology](README.md) finding 7 |
| M14 | alternative nitrogenases | Vanadium and iron-only nitrogenases beside the molybdenum enzyme, since molybdenum may limit fixation before nitrogen does | G30, G47 | hypothesis | [ecology](README.md) finding 12; [literature](literature.md) Nitrogen |
| M15 | wet-leaf disease resistance | Resistance to leaf pathogens through a hundred hours and more of leaf wetness, with sunlight disinfecting at 1-6% of Earth's rate | G16 | hypothesis | [ecology](README.md) finding 3 |
| M16 | root regulation at lunar gravity | Root growth regulated by design, since plants sense gravity poorly below 0.1-0.3 g | G16, G18 | hypothesis | [ecology](README.md) finding 19 |
| M17 | ethylene balance | Low ethylene emission by plants and fast ethylene uptake by soils: the vegetated Moon's own ethylene settles at 94-233 ppb with Earth's soil uptake, above the 50 ppb that costs wheat 36% of its yield | G16, G46 | hypothesis | [ecology](README.md) finding 24; [people and land](people_and_land.md) item 16 |
| M18 | gas-barrier envelope | A mostly inert, replaceable hydrogen barrier grown and repaired by living tissue, behind a tougher outer surface | G01 | hypothesis | [aerophyte study](../../research/studies/aerophytes/README.md) The reference organism; [aerophyte gas biology](../../research/studies/aerophytes/gas_biology.md) |
| M19 | hydrogen organ or symbiont | Hydrogen made by a dedicated organ or symbiont, from light or by dark fermentation, with its light, substrates and oxygen management counted | G01 | hypothesis | [aerophyte biology](../../research/studies/aerophytes/biology.md); [aerophyte gas biology](../../research/studies/aerophytes/gas_biology.md) The swim bladder's gas gland |
| M20 | water from the air | Water drawn from rain, cloud and the humid air through sorbent skins and fringes, in preference to descending to drink | G01 | hypothesis | [decisions](../../research/decisions.md); [aerophyte study](../../research/studies/aerophytes/README.md) The limits, in brief |
| M21 | attached buds | Daughters grown as attached buds, fed carbon, nutrients and gas by the parent until they can carry themselves | G01 | hypothesis | [aerophyte biology](../../research/studies/aerophytes/biology.md) |
| M22 | spectral tuning | Paler upper leaves, photosynthesis into the far red and clumped foliage; they add 3%, 6% and 3% on the Moon, where the diffuse sky already does most of what they do | — | set aside | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) A2 |
| M23 | extra photoprotection and daily rests | Reflective hairs, folding leaves and a daily rest in the long day: peak light is 30% below Earth's and the UV index 0.11, and a rest halving photosynthesis for eight hours in 24 costs about a quarter of the growth | — | set aside | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) A4, A5 |
| M24 | canopy regrowth each cycle | Shedding the canopy before dusk and regrowing it at dawn costs 264 g C per m2 a cycle and leaves little or no growth | — | set aside | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) B3 |
| M25 | frost protection | Antifreeze, insulation and heat sharing in plant mats: land nights stay near 19-20 degrees C at low and middle latitudes | — | set aside | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) B4 |
| M26 | hyper-C4 as the general strategy | C4 everywhere: the canopy works mostly in dim diffuse light, where C4's extra two ATP per CO2 cost most | — | set aside | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) A3 |
| M27 | tethered sails | A wing hung on a tether in another layer of air, or a drogue in the sea, that sails a buoyant body across the wind | G01 | hypothesis | [aerophyte sailing](../../research/studies/aerophytes/sailing.md) Sailing between two layers |
| M28 | compartmented gas | Gas held in many small wet compartments or modules, so that a tear or a lightning strike empties or burns one at a time | G01 | hypothesis | [aerophyte study](../../research/studies/aerophytes/README.md) The limits, in brief |

### Communities

| ID | Community | Members | Status | Source |
|---|---|---|---|---|
| C01 | the sky archipelago | aerophytes (G01), lesser aerophytes (G02), giant aerophytes (G03), resident aeroplankton phototrophs (G04), resident aeroplankton heterotrophs (G07), aerial collectors (G08) | hypothesis | [aerophyte study](../../research/studies/aerophytes/README.md) What a giant looks like; [aerial study](../../research/studies/aerial_ecology/README.md) The cycle and the community |
| C02 | cloud-water aeroplankton | resident aeroplankton phototrophs (G04), aerial cyanobacteria (G05), aerial microalgae (G06), resident aeroplankton heterotrophs (G07) | hypothesis | [aerial study](../../research/studies/aerial_ecology/README.md) What the evidence supports; [people and land](people_and_land.md) item 3 |
| C03 | the films of the high air | high-altitude films (G10) | carried forward | [decisions](../../research/decisions.md) Life and people |
| C04 | forests of wood and fruit | evergreen idlers (G16), day-fruit trees (G19), trap plants (G26), droppings-collecting plants (G27), glowing lures (G28), night decomposers of the fruit fall (G32), night insects of the fruit fall (G33), fruit-eaters (G34) | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) F2 |
| C05 | plains of storage organs | storage-organ plants (G23), dusk-seeding grasses and forbs (G24), grazers and seed-eaters of the plains (G35) | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) F1 |
| C06 | wet margins | trap plants (G26), lake insects (G37) | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) F3 |
| C07 | the dusk fruit fall and the night food web | day-fruit trees (G19), night decomposers of the fruit fall (G32), night insects of the fruit fall (G33), fruit-eaters (G34), glowing lures (G28), sea-forest migrant insects (G13), dusk followers (G15) | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) E1, E2 |
| C08 | the sea-forest return path | sea-forest migrant insects (G13), seabird-like colonial fliers (G14), run fish (G42), trap plants (G26), droppings-collecting plants (G27), glowing lures (G28), sulfur- and iodine-releasing seaweeds (G40) | carried forward | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) G1-G4, H1-H10 |
| C09 | farms of the day fruit | day-fruit crops (G18), fruit picked young (G21), continuous-light tomatoes (G22) | carried forward, in places | [ecology register](../../research/studies/lunar_cycle_ecology/README.md) C1, D3, D5; [people and land](people_and_land.md) item 11 |
| C10 | megaforests | megaforest emergents (G29), lesser aerophytes (G02) | hypothesis | [megaforest](../../research/studies/megaforest_wind/README.md); [forest patch](../../research/studies/forest_patch/README.md); [biosphere](../README.md); [aerophyte growth](../../research/studies/aerophytes/growth.md) The canopy nursery |
| C11 | pioneer crusts | basalt pioneers (G30), crust mosses (G31) | hypothesis | [ecology](README.md) finding 13 |
| C12 | the living seas | sea phytoplankton (G38), dark-night luminous plankton (G39), sulfur- and iodine-releasing seaweeds (G40), sea grazers (G41), run fish (G42), sea-forest migrant insects (G13) | carried forward | [decisions](../../research/decisions.md) Life and people; [living water](../living_water/waters.json); [sea appearance](../../research/studies/sea_appearance/README.md) Decisions |

The backbone ([taxa.py](taxa.py), `BACKBONE`) has 25 nodes from Life down to *Solanum lycopersicum*: the domains
Bacteria (with Cyanobacteriota) and Eukaryota; Plantae with Chlorophyta, Bryophyta and Tracheophyta, the rankless
Angiospermae and, within it, Caryophyllales, Nepenthaceae, *Nepenthes*, Cucurbitaceae and *Solanum*; Fungi; and
Animalia with Arthropoda (Insecta, Diptera, Arachnida, Araneae) and Chordata (Aves, Salmonidae). Ranked nodes follow
the GBIF backbone; the rankless ones are cited.

## Names for the aerophytes

The author chose the everyday names on 10 October 2026: aerophytes for the class (G01), and sky reefs for the colony
giants, the rafts of gas modules (G03); the round giants have no common name yet. Until then the working name was
Sagan & Salpeter's "floaters", their word in 1976 for buoyant organisms of Jupiter's air. "Aerophyte" coincides with
an older botanical term for epiphytes; the project uses it for its own class. The scientific names stay open: N1–N6
remain candidates for that layer, which the author said cloud roots could serve.

The system distinguishes three stages of the Open Moon's aerophytes: the young, drifting propagules; the lesser
aerophytes of the coupled study's reference organism, 40–150 m across (G02); and the giants that keep growing for
centuries, as sky reefs or round bodies (G03). Each candidate names the aerophytes and their giants, and some name the
young or the formation they drift in. Every scientific form was checked against IRMNG and the GBIF backbone on
9 October 2026, and neither lists any of them as a genus.

| ID | The aerophytes | The giants | Scientific form | Roots |
|---|---|---|---|---|
| N1 | welkins | elder welkins | *Welkinia* | Old English wolcen, cloud, sky |
| N2 | holms | great holms | *Aeronesia* | Old English holm, Old Norse holmr, islet; Greek aēr, nēsos |
| N3 | meghas | meghadutas | *Meghaphyta*, *Meghaduta* | Sanskrit megha, cloud; Kālidāsa's Meghadūta |
| N4 | kunpeng (the young: kun) | peng | *Kunpengia* | Classical Chinese kūn and péng, Zhuangzi |
| N5 | serenes | serenissimas | *Sereniphyta*, *Serenissima* | Latin serenus; Mare Serenitatis |
| N6 | mawingu (one: a wingu) | great mawingu | *Mawingua* | Swahili wingu, plural mawingu, cloud |

### N1. Welkins and elder welkins

- **Meaning and roots.** English *welkin*, the vault of the sky in poetry, comes through Middle English *welkne*
  (weather, the heavens, earlier a cloud) from Old English *wolcen*: "a cloud; the clouds; the heavens; the sky"
  (Bosworth–Toller). Old English also has *wolcen-faru*, "the cloud-host, the moving clouds", free to name the
  formation they drift in.
- **Sound.** Two soft syllables, WEL-kin, an old word that ends on "kin".
- **Image.** The dome of the sky in old poetry, and living pieces of it drifting over the seas.
- **Logic.** The word meant cloud before it meant sky; these organisms are the sky's living clouds.
- **Place.** The English working name of the aerophytes (G01); *Welkinia* the clade name of the new lineage; elder
  welkins the giants (G03); the wolcenfaru their formation.
- **Checks.** The word is rare in modern English and carries no biological meaning.

### N2. Holms and great holms

- **Meaning and roots.** Old English *holm* means "a mound, hill, rising ground; wave, ocean, water, sea; land rising
  from the water, an island in a river" (Bosworth–Toller), and Old Norse *holmr* an islet. *Aeronesia* joins Greek
  ἀήρ (aēr), air, and νῆσος (nēsos), island, as Polynesia joins "many" and "island".
- **Sound.** One round syllable, "holm"; *Aeronesia* flows, air-oh-NEE-zha.
- **Image.** A drifting archipelago of green islets over the sea, each holding its own small pools and residents.
- **Logic.** The coupled study's concept is a light host carrying small communities in retained water: an island. Its
  colonial giants are rafts of modules, an archipelago of holms.
- **Place.** The English working name of the aerophytes (G01); *Aeronesia* the clade; great holms the giants (G03), a
  colony's modules its holms; the sky archipelago (C01) their formation.
- **Checks.** *Holmia*, the obvious Latin form, is preoccupied by a trilobite genus and others, so the clade takes the
  Greek form. "Holm" also names the holm oak; context separates them.

### N3. Meghas and meghadutas

- **Meaning and roots.** Sanskrit मेघ *megha* is "'sprinkler', a cloud", from the root *mih*, and also "a mass, a
  multitude"; *meghadūta* is "cloud-messenger", the name of a celebrated poem by Kālidāsa (Monier-Williams 1899). The
  Meghadūta "describes the complaint of an exiled lover, and the message he sends to his wife by a cloud"
  (Encyclopaedia Britannica 1911).
- **Sound.** MAY-gha, with a breathed g; me-gha-DOO-ta flows in four syllables. To English ears "megha" also carries
  "mega".
- **Image.** Kālidāsa's monsoon cloud carrying a message across the land; giants carrying their communities across the
  Moon.
- **Logic.** A cloud that bears water and messages; a root that means one who sprinkles, for organisms that gather and
  hold water.
- **Place.** Meghas the working name of the aerophytes (G01), with *Meghaphyta* the clade; meghadutas the giants (G03),
  with *Meghaduta* free for a genus if the giants become a sister lineage.
- **Checks.** The macron of meghadūta drops under the codes' alphabet rule. The poem is secular.

### N4. Kun and peng

- **Meaning and roots.** Zhuangzi's first chapter opens: "In the northern darkness there is a fish, its name is Kun. ...
  It changes into a bird, its name is Peng. The back of the Peng, no one knows how many thousand li it measures; when it
  rises and flies, its wings are like clouds hanging from the sky" (北冥有魚，其名爲鯤 ... 化而爲鳥，其名爲鵬 ...
  其翼若垂天之雲). The Peng rises ninety thousand li on the whirlwind, and one old reading has it rest only after six
  months. The same chapter says that "if
  the wind is not piled up deep, it has no strength to bear great wings" (風之積也不厚，則其負大翼也无力). In the
  Erya, 鯤 *kun* means fish roe, and Guo Pu glosses it as the young of all fish (Erya zhushu); the commentary on Yang
  Shen's Yiyu tuzan says that Zhuangzi made the smallest thing into the greatest. Translations here are this work's.
- **Sound.** Two short, bell-like syllables, koon and pung.
- **Image.** The Peng's back thousands of li across, rising on deep air and, in that reading, flying half a year
  before it rests.
- **Logic.** A speck-sized propagule becomes a giant that lives for centuries; the Moon's dense, tall air is Zhuangzi's
  wind piled deep enough to bear great wings. The coupled study raises its young in megaforest crowns before they rise.
- **Place.** Kun the young (propagules and juveniles), peng the giants (G03), kunpeng the lineage (G01), with
  *Kunpengia* its clade name.
- **Checks.** *Kunpengia* is free. *Peng* (Lu & Li, 2023) and *Pengia* (Sharpe, 1883; Geyer & Corbacho, 2015) are animal
  genera, so the giants' scientific form needs another word, and the pterosaur genus *Kunpengopterus* (Wang et al.,
  2010) already uses the root. Kunpeng is a familiar compound in Chinese and a trade name; "peng" is also English slang.

### N5. Serenes and serenissimas

- **Meaning and roots.** Latin *serenus* is "clear, fair, bright", said of the sky: *caelo sereno* (Lewis & Short). Mare
  Serenitatis, the "Sea of Serenity", is one of the seas the Near-side Sea joins. English keeps *serene* as a poetic noun,
  "the serene of heaven", and has a homonym *serene* or *serein* for a fine rain falling from a cloudless sky after
  sunset (Wiktionary). *Serenissima* is the Latin superlative, "the most serene", once the title of Venice.
- **Sound.** se-REEN, smooth and long; se-re-NIS-si-ma for the giants, five open syllables.
- **Image.** Calm daylight over the Sea of Serenity, and a fine rain falling from clear air after sunset.
- **Logic.** The Moon's own sea name; organisms that keep to clear daylight and drink from the air at dusk.
- **Place.** Serenes the working name of the aerophytes (G01), with *Sereniphyta* the clade; serenissimas the giants
  (G03), with *Serenissima* free for a genus.
- **Checks.** *Serenitas* (Wells, 2009) is preoccupied, and *Serenia* sounds too close to Sirenia, the sea cows.

### N6. Mawingu and great mawingu

- **Meaning and roots.** Swahili *wingu*, plural *mawingu*, means cloud (Wiktionary).
- **Sound.** ma-WEEN-goo, three open syllables stressed on the second; the singular *wingu*.
- **Image.** Tall clouds over East Africa's plains and lakes.
- **Logic.** A word from one of Africa's most widely spoken languages, for organisms that cross every nation's sky.
- **Place.** Mawingu the working name of the aerophytes (G01), a wingu one of them, with *Mawingua* the clade; great
  mawingu the giants (G03).
- **Checks.** Swahili is a living language, so the name is adopted with its speakers (P6).

The register carries the author's common name, aerophytes, for G01, and the identifier stays the same when the
scientific name is chosen.

## Files and checks

| File | Holds |
|---|---|
| [taxa.py](taxa.py) | The system's vocabularies, the backbone, the register, the name candidates, the computed cycle and size classes, and `check()` |
| [results/taxa.json](results/taxa.json) | Its product (schema `terluna.biosphere.ecology-taxa/1`): the files it reads bound by hash, the files it cites listed by path |
| [taxonomy_sources.json](taxonomy_sources.json) | What was read for the system and the names, and how far |
| [../tests/test_taxa.py](../tests/test_taxa.py) | The checks |

`python -m biosphere.ecology.taxa` rewrites the product; it refuses to write a register that fails `check()`.
`python -m pytest biosphere/tests/test_taxa.py` checks the schema, unique names and identifiers, that every entry cites
an existing project file, the ranks and parentage of the backbone and the groups, the vocabularies, the cycle against
the shared constants, the size classes against the first screen, the bindings of producer and inputs, that the product
matches what the producer writes, and that this document lists every entry. The product hash-binds the two files
taxa.py reads, the first screen and the sources register, and lists the project files its entries cite by path in
`producer.cited_files`, so routine edits to those records leave it current; renaming or removing a cited file calls
for an update.
