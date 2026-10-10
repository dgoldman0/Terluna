"""A taxonomic system for the Open Moon's biosphere, and a first register of its designed groups.

    python -m biosphere.ecology.taxa
    # -> biosphere/ecology/results/taxa.json

The author asked on 9 October 2026 to rename the aerophytes and to begin a taxonomic system for the biosphere
(research/decisions.md). This module holds the system's vocabularies, an Earth backbone the designed groups hang
from, the register of designed groups, design modules and communities drawn from the project's studies, the common
names the author chose for the aerophytes on 10 October 2026 and the candidates for their scientific names.
taxonomy.md explains it; taxonomy_sources.json records what was read.

Everything registered here is a design. Each entry carries the project's own verdict on it, the Earth lineages it
draws on as reference points, its proposed placement by descent (an anchor in the backbone and a degree of design)
and by function (life form, settings, guilds, night strategy, clock), and the project files it comes from. The
product hash-binds the files this module reads (the sources register and the first screen) and lists the other cited
files by path, so routine edits to the records it cites leave it current. check() validates the register's
vocabularies, identifiers, ranks and parentage, and is run before the product is written.

Two blocks are computed: the lunar cycle and the speed of dusk by latitude (shared constants), and the residence
times of the aeroplankton size classes at lunar gravity (read from this branch's first screen).
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

from shared.constants import MOON_RADIUS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.biosphere.ecology-taxa/1'
OUT = HERE / 'results' / 'taxa.json'
SOURCES = HERE / 'taxonomy_sources.json'
SCREEN = HERE / 'results' / 'first_screen.json'
SCREEN_SCHEMA = 'terluna.biosphere.ecology-first-screen/1'
VERSION = '0.1, first proposal of 9 October 2026'

# Project files the entries cite (paths from the repository root).
LCE = 'research/studies/lunar_cycle_ecology/README.md'
BIO = 'biosphere/README.md'
ECO = 'biosphere/ecology/README.md'
LIT = 'biosphere/ecology/literature.md'
STATE = 'biosphere/ecology/state.md'
PEOPLE = 'biosphere/ecology/people_and_land.md'
AERIAL = 'research/studies/aerial_ecology/README.md'
AERIAL_SRC = 'research/studies/aerial_ecology/sources.json'
AEROPHYTES = 'research/studies/aerophytes/README.md'
AEROPHYTES_BIO = 'research/studies/aerophytes/biology.md'
AEROPHYTES_SAIL = 'research/studies/aerophytes/sailing.md'
AEROPHYTES_GAS = 'research/studies/aerophytes/gas_biology.md'
AEROPHYTES_GROW = 'research/studies/aerophytes/growth.md'
MEGA = 'research/studies/megaforest_wind/README.md'
PATCH = 'research/studies/forest_patch/README.md'
SEA = 'research/studies/sea_appearance/README.md'
WATERS = 'biosphere/living_water/waters.json'
DEC = 'research/decisions.md'
JULY = 'archive/july_2026_knowledge_bundle/terluna_moon_terraforming_bundle/06_BIOSPHERE_AND_EVOLUTIONARY_DESIGN.md'

EVIDENCE = ('A proposed taxonomic system and a first register of design groups drawn from the project\'s studies. '
            'Every entry is a design: no organism named here has been built, tested or released. Each status repeats '
            'the verdict of the project files the entry cites. The cycle block is arithmetic on shared constants; the '
            'aeroplankton residence times are dry settling at lunar gravity from the first screen.')
READING_RULE = ('Entries are keyed by permanent identifiers (G groups, M design modules, C communities, N name '
                'candidates); working names can change and identifiers never do. A group\'s anchor is the smallest '
                'backbone clade holding all its candidate chassis; its degree is the rank its founding lineages would '
                'take (system.degrees). form, settings, guilds, night and clock use the vocabularies in '
                'system.functional; "undetermined" means the cited sources leave it open. A subgroup inherits its '
                'parent\'s modules. status is the project\'s verdict (carried_forward, test_first, hypothesis, '
                'set_aside) with scope throughout or in_places, and basis says where the verdict comes from. '
                'An experience that shows an entry labels it as informed by a design. producer.inputs hash-binds the '
                'two files the producer reads; producer.cited_files lists the project files the entries cite, by path. '
                'cycle gives hours, km/h and m/s; size classes give micrometres and days.')

# ---- Vocabularies ----------------------------------------------------------------------------------------------
STATUSES = {
    'carried_forward': 'Selected for the design: the lunar-cycle ecology register\'s carry-forward and in-places lists, '
                       'or a standing entry in the decisions register',
    'test_first': 'Selected with a named test to pass before adoption',
    'hypothesis': 'Proposed, required or opened for research by a project study, without a selection yet',
    'set_aside': 'Recorded and not pursued',
}
STATUS_ORDER = ('set_aside', 'hypothesis', 'test_first', 'carried_forward')
SCOPES = ('throughout', 'in_places')
BASES = {
    'selected': 'The selection of the lunar-cycle ecology register (research/studies/lunar_cycle_ecology)',
    'decision': 'A standing entry in research/decisions.md',
    'direction': 'A research direction the author has opened',
    'requirement': 'A project screen finds the biosphere needs it; how it is met is open',
    'idea': 'Proposed or used as a reference case in a project study',
}
RANKS = ('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
DEGREES = {
    'line': dict(
        meaning='A founding line within an Earth species: timing, tolerance or allocation changed, the species kept',
        name_form='Earth binomial + founding identifier + optional line name in quotes, e.g. Solanum lycopersicum '
                  'OM-L-0007 \'Long Dusk\'',
        precedents='Mycoplasma mycoides JCVI-syn3.0; Saccharomyces cerevisiae Sc2.0; Zea mays event MON810 '
                   '(MON-ØØ81Ø-6); genetically modified cultivars (ICNCP Art. 2.19)',
        anchor_at_or_above='species'),
    'species': dict(
        meaning='A designed species within an Earth genus, with traits beyond any Earth species of it',
        name_form='a new binomial in the Earth genus', precedents='a new species under the chassis\'s code',
        anchor_at_or_above='genus'),
    'genus': dict(
        meaning='A designed genus within an Earth family or larger clade',
        name_form='a new genus and species placed in the Earth family', precedents='a new genus under the chassis\'s code',
        anchor_at_or_above='family'),
    'lineage': dict(
        meaning='A new lineage: a body plan or metabolism no Earth genus holds, placed within its chassis\'s Earth clade',
        name_form='a clade name defined from its founding lineages, with genera and species inside it',
        precedents='a clade under the PhyloCode', anchor_at_or_above='family'),
}
REALMS = {'sky': 'The air, from the ground to 90 km', 'land': 'The ground, its soils and its canopies',
          'water': 'Seas, lakes, wetlands and rivers'}
SETTINGS = {
    'low_sky': dict(realm='sky', definition='0-10 km above sea level: the day\'s mixed layer, cloud base 4-7 km'),
    'storm_layer': dict(realm='sky', definition='10-22 km: storms to their median top'),
    'cool_air': dict(realm='sky', definition='22-35 km: the air freezes near 26 km'),
    'flight_band': dict(realm='sky', definition='35-45 km: the long-haul flight band'),
    'high_air': dict(realm='sky', definition='45-90 km'),
    'canopy': dict(realm='land', definition='Crowns, branches and epiphyte perches of forests'),
    'ground': dict(realm='land', definition='The ground surface, its litter, soils and plains'),
    'seas': dict(realm='water', definition='Open seas, productive coasts, river mouths and tidal flats'),
    'lakes_and_wetlands': dict(realm='water', definition='Rain-fed lakes, their shores and wet margins'),
    'rivers': dict(realm='water', definition='Rivers from the lakes to the seas'),
}
FORM_CLASSES = {
    'aerial': 'Organisms that live substantially in the air (the aerosphere of Kunz et al. 2008)',
    'plant': 'Rooted or attached photosynthetic forms, classed after Raunkiaer by how they carry the night',
    'ground': 'Fungi and animals of the land',
    'microbial': 'Free-living microbes of soils and waters',
    'aquatic': 'Forms of the waters',
}
FORMS = {
    'aeroplankton': dict(cls='aerial', definition='Organisms the air carries, too small to hold a course; resident '
                         'aeroplankton complete their lives aloft, transient aeroplankton are dispersing stages'),
    'drifter': dict(cls='aerial', definition='Buoyant or ballooning organisms that move with the wind and control '
                    'mainly their height'),
    'sailor': dict(cls='aerial', definition='Buoyant organisms that steer across the wind with a wing or drogue hung '
                   'in another layer of air or in the sea, as Physalia sails the sea surface'),
    'giant': dict(cls='aerial', definition='Buoyant organisms that keep growing for centuries, like the oldest trees'),
    'flyer': dict(cls='aerial', definition='Animals that fly under their own power'),
    'glider': dict(cls='aerial', definition='Animals that launch from canopies and cliffs and glide without powered '
                   'flight'),
    'film': dict(cls='aerial', definition='Sparse microbial or algal films on particles, droplets or platforms in '
                 'the high air'),
    'tenant': dict(cls='aerial', definition='Organisms living on or within aerophytes and other aerial hosts, in their '
                   'retained-water pockets'),
    'tree': dict(cls='plant', definition='Woody plants whose living crowns stand in the air through the night '
                 '(Raunkiaer\'s phanerophytes)'),
    'vine': dict(cls='plant', definition='Climbing or trailing plants'),
    'herb': dict(cls='plant', definition='Non-woody plants whose shoots stand through the night'),
    'geophyte': dict(cls='plant', definition='Plants that carry the night in bulbs, corms or tubers below ground '
                     '(Raunkiaer\'s cryptophytes)'),
    'epiphyte': dict(cls='plant', definition='Plants living on other plants, rooted in the canopy'),
    'crust': dict(cls='plant', definition='Crusts and mats of cyanobacteria, mosses and lichens on the ground'),
    'fungal': dict(cls='ground', definition='Yeasts and mycelial fungi'),
    'ground_animal': dict(cls='ground', definition='Animals that walk, climb or burrow'),
    'microbe': dict(cls='microbial', definition='Free-living microbes of soils, sediments and waters'),
    'plankton': dict(cls='aquatic', definition='Organisms the water carries'),
    'nekton': dict(cls='aquatic', definition='Animals that swim against currents'),
    'benthos': dict(cls='aquatic', definition='Organisms of the beds of seas and lakes'),
    'macroalga': dict(cls='aquatic', definition='Seaweeds'),
}
GUILDS = {
    'producer': 'Fixes carbon from light (or from chemical energy)',
    'grazer': 'Eats living plant or algal tissue',
    'fruit_eater': 'Eats fruit', 'seed_eater': 'Eats seeds',
    'predator': 'Catches and eats animals (carnivorous plants included)',
    'collector': 'Gathers airborne particles and small organisms on webs or filters',
    'consumer': 'Eats other organisms; diet not yet specified',
    'decomposer': 'Breaks down dead organic matter',
    'pollinator': 'Carries pollen between flowers',
    'fixer': 'Fixes atmospheric nitrogen',
    'soil_builder': 'Weathers rock and builds soil',
    'air_keeper': 'Removes trace gases the air cannot remove: methane, N2O, hydrogen, carbon monoxide, ethylene',
    'carrier': 'Carries elements back from the seas to the land on a return path',
    'receiver': 'Takes up the nutrients carriers deliver, at the forest end of a return path',
    'lure': 'Draws animals by light or scent',
    'undetermined': 'The cited sources leave the role open',
}
NIGHT = {
    'store': 'Carries the night on stored carbon',
    'idle': 'Slows its upkeep through the night, on a cue',
    'dormant': 'Passes the night or longer in a dormant state: seed, spore, cyst or dried film',
    'twilight': 'Uses the long twilight\'s light',
    'night_active': 'Active through the night on stores, the fruit fall or prey',
    'terminator_following': 'Travels to change its light: with dusk or the Sun, east on the wind to meet the next '
                            'dawn sooner, or across the wind toward the twilit poles',
    'shedding': 'Sheds its canopy before dusk and regrows it',
    'undetermined': 'The cited sources leave the strategy open',
}
CLOCKS = {
    'dusk_set': 'A timer started at dusk',
    'dawn_set': 'A program set at dawn',
    'carbon_cued': 'The organism\'s own carbon balance as the cue',
    'two_clocks': 'Daily rhythms inside the lunar environmental clock',
    'earth_clock_with_cues': 'An Earth 24-hour clock kept on daily cues, on farms',
    'undetermined': 'The cited sources leave the clock open',
}
ELEMENTS = ('P', 'N', 'S', 'I', 'K', 'trace elements')
TREE_HEIGHT_CLASSES = (  # Raunkiaer 1934: phanerophytes by height, m
    dict(name='nano', range_m=[0.0, 2.0]), dict(name='micro', range_m=[2.0, 8.0]),
    dict(name='meso', range_m=[8.0, 30.0]), dict(name='mega', range_m=[30.0, None]))
AEROPLANKTON_CLASSES = (  # Sieburth et al. 1978 by decades; the three smallest as reproduced in Brotas et al. 2022
    dict(name='pico', range_um=[0.2, 2.0]), dict(name='nano', range_um=[2.0, 20.0]),
    dict(name='micro', range_um=[20.0, 200.0]))
SIZE_NAMES = {c['name'] for c in TREE_HEIGHT_CLASSES} | {c['name'] for c in AEROPLANKTON_CLASSES}

# ---- The system --------------------------------------------------------------------------------------------------
LAYERS = (
    dict(name='descent', rule='Each designed lineage keeps the Earth lineage of its chassis; the register places it '
         'at an anchor in the backbone', sources=['gibson_2010', 'hutchison_2016', 'richardson_2017']),
    dict(name='design', rule='How far a design departs sets its degree (line, species, genus, lineage) and its name '
         'form', sources=['icncp_2016', 'eu_65_2004', 'bch_mon810', 'icnp_2022']),
    dict(name='founding', rule='Every founding lineage is registered with a permanent identifier, a type genome and '
         'a preserved sample, its purpose, reference organisms, release and versions', sources=['seqcode_1_2',
         'biocode_2011', 'icn_madrid_2025', 'richardson_2017']),
    dict(name='descendants', rule='Clades after release are named by phylogenetic definitions with the founding '
         'lineages as specifiers; descendant species get binomials when they become distinct',
         sources=['phylocode_6']),
    dict(name='function', rule='Life form, settings, guilds, night strategy and clock, classed separately from '
         'descent', sources=['raunkiaer_1934', 'sieburth_1978', 'kunz_2008', 'adl_2019']),
    dict(name='names', rule='Scientific names in Latin grammar from any language, with stated etymology and no Earth '
         'homonym; common names recorded as the lunar peoples give them', sources=['icn_madrid_2025', 'iczn_1999',
         'seqcode_1_2', 'icnp_2022', 'icncp_2016']),
)
PRINCIPLES = (
    dict(id='P1', text='A designed organism keeps the Earth lineage of its chassis, the lineage whose genome organizes '
         'the cell and is inherited; genes, pathways and symbionts from other lineages are recorded as contributions',
         sources=['gibson_2010', 'hutchison_2016', 'richardson_2017', 'fredens_2019']),
    dict(id='P2', text='How far a design departs sets its degree, and each degree takes a naming form Earth already '
         'uses: a line within an Earth species, a designed species, a designed genus, or a new lineage',
         sources=['icncp_2016', 'eu_65_2004', 'bch_mon810', 'hutchison_2016']),
    dict(id='P3', text='Every founding lineage is registered with a permanent identifier, a type genome deposited as '
         'a sequence and a preserved sample, with its purpose, reference organisms, release and versions',
         sources=['seqcode_1_2', 'icnp_2022', 'icn_madrid_2025', 'biocode_2011', 'richardson_2017']),
    dict(id='P4', text='Descendants are named by descent from the founders: a clade is an ancestor and all its '
         'descendants, so a clade defined from a founding lineage keeps its meaning however far its descendants '
         'diverge', sources=['phylocode_6']),
    dict(id='P5', text='Function is classified apart from descent, by how an organism lives through the 354-hour '
         'night and where it lives; one lineage can take several forms in its life, classed by the adult',
         sources=['raunkiaer_1934']),
    dict(id='P6', text='Scientific names are treated as Latin, spelled in the 26-letter alphabet, drawn from any '
         'language including the Moon\'s own names, with the etymology stated, no homonym of any Earth genus and '
         'nothing derogatory to any group of people', sources=['icn_madrid_2025', 'iczn_1999', 'seqcode_1_2',
         'icnp_2022', 'iau_gazetteer']),
    dict(id='P7', text='Common names belong to the peoples who live with the organisms and are recorded as they arise, '
         'with language, community, date and meaning', sources=['icncp_2016', 'iczn_1999']),
    dict(id='P8', text='Identifiers never change; names may', sources=['eu_65_2004', 'seqcode_1_2']),
)
CHOICES = (
    dict(id='O1', question='Which rules govern new designed taxa (degrees species, genus and lineage)?',
         options={'a': 'The Earth code of each chassis (ICN, ICZN, ICNP or SeqCode)',
                  'b': 'One lunar code for all designed taxa, in forms every Earth code accepts',
                  'c': 'The Earth codes alone, keeping every design as a line within an Earth species'},
         recommendation='b'),
    dict(id='O2', question='Ranks for designed taxa and their descendants',
         options={'a': 'Linnaean ranks at founding, rankless clade names for descendants',
                  'b': 'Rankless clade names throughout, with binomials only for species'},
         recommendation='a'),
    dict(id='O3', question='The founding identifier',
         options={'a': 'A register number with a check digit, e.g. OM-L-0001-7',
                  'b': 'Three components after the OECD and EU format: designer code, lineage code, check digit'},
         recommendation=None),
    dict(id='O4', question='Roots for scientific names',
         options={'a': 'Latin and Greek, as the ICNP advises',
                  'b': 'Any language in Latin grammar, as the ICZN, ICN and SeqCode allow',
                  'c': 'A lunar root stock (the IAU names and the Moon\'s phenomena) for founders, any language for '
                       'descendants'},
         recommendation=None),
    dict(id='O5', question='Where founding ends and descent begins',
         options={'a': 'Every release is a founding lineage, a re-engineered release a new one derived from the old; '
                       'change after release is descent',
                  'b': 'Founding ends with the 500-year build; later releases are descendants'},
         recommendation='a'),
    dict(id='O6', question='Giants and lesser aerophytes',
         options={'a': 'One lineage, the giants a form its oldest members reach',
                  'b': 'Two sister lineages, the giants designed for indeterminate growth'},
         recommendation=None),
    dict(id='O7', question='The aerophytes\' scientific name',
         options={c: c for c in ('N1', 'N2', 'N3', 'N4', 'N5', 'N6')}, recommendation=None,
         note='The author chose the common names on 10 October 2026: aerophytes, and sky reefs for the colony giants. '
              'N1-N6 remain candidates for the scientific names, a layer the author said cloud roots could serve.'),
)
NAMING_RULES = (
    dict(rule='Names are treated as Latin, whatever their origin', sources=['icn_madrid_2025:Principle V',
         'icnp_2022:Principle 3', 'seqcode_1_2:Principle 3', 'biocode_2011:Principle VI']),
    dict(rule='Only the 26 letters of the Latin alphabet; diacritics are dropped (ä becomes ae, é becomes e)',
         sources=['iczn_1999:Art. 11.2', 'icn_madrid_2025:Art. 60.7', 'seqcode_1_2:Rule 8']),
    dict(rule='A species name is a binary combination: a genus name and one epithet', sources=[
         'icn_madrid_2025:Art. 23.1', 'iczn_1999:Art. 5.1', 'icnp_2022:Rule 30']),
    dict(rule='No tautonyms: the epithet does not repeat the genus name', sources=['icn_madrid_2025:Art. 23.4']),
    dict(rule='Epithets of 2 to 30 characters', sources=['icn_madrid_2025:Art. 23.2', 'icncp_2016:Art. 21.13']),
    dict(rule='Genus names and epithets may come from any language, or be arbitrary words usable as words',
         sources=['iczn_1999:Art. 11.3', 'icn_madrid_2025:Art. 20.1, 23.2', 'seqcode_1_2:Rule 26']),
    dict(rule='The etymology of every new name is stated', sources=['icnp_2022:Rule 6', 'seqcode_1_2:Rule 26',
         'icn_madrid_2025:Rec. 60I.1']),
    dict(rule='No homonym of any genus under any Earth code, checked against IRMNG and the GBIF backbone',
         sources=['seqcode_1_2:Principle 2, Rule 9b', 'icnp_2022:Principle 2', 'iczn_1999:Rec. 1A', 'irmng',
                  'gbif_backbone']),
    dict(rule='No name derogatory to a group of people; names from a living language are adopted with its speakers',
         sources=['icn_madrid_2025:Art. 51.2, Rec. 51A.1', 'iczn_1999:Appendix A.4']),
    dict(rule='No standard prefix or suffix stamped on Earth names to mark lunar membership', sources=[
         'iczn_1999:Art. 1.3.7']),
    dict(rule='Registration fixes a name\'s date and priority', sources=['seqcode_1_2:Rule 26', 'iczn_1999:Art. 8.5',
         'icn_madrid_2025:Art. 42.1', 'biocode_2011:Art. 5.2', 'phylocode_6:Art. 8.1']),
    dict(rule='Common names and line names are words of living languages and may be recorded in any of them',
         sources=['icncp_2016:Art. 21.11']),
)
LUNAR_ROOTS = (  # IAU lunar water features with their stated origins (iau_gazetteer)
    ('Mare Nubium', 'Sea of Clouds'), ('Mare Vaporum', 'Sea of Vapors'), ('Mare Imbrium', 'Sea of Showers'),
    ('Oceanus Procellarum', 'Ocean of Storms'), ('Sinus Roris', 'Bay of Dew'), ('Sinus Iridum', 'Bay of Rainbows'),
    ('Mare Serenitatis', 'Sea of Serenity'), ('Mare Tranquillitatis', 'Sea of Tranquility'),
    ('Mare Insularum', 'Sea of Islands'), ('Mare Humorum', 'Sea of Moisture'), ('Mare Undarum', 'Sea of Waves'),
    ('Mare Spumans', 'Foaming Sea'), ('Mare Fecunditatis', 'Sea of Fecundity'), ('Mare Nectaris', 'Sea of Nectar'),
    ('Lacus Somniorum', 'Lake of Dreams'), ('Lacus Temporis', 'Lake of Time'),
    ('Lacus Perseverantiae', 'Lake of Perseverance'), ('Lacus Spei', 'Lake of Hope'), ('Palus Somni', 'Marsh of Sleep'),
    ('Sinus Aestuum', 'Seething Bay'),
)
FOUNDING_FIELDS = (
    'identifier', 'name', 'degree', 'group', 'chassis (Earth name with authority)', 'contributions (donor Earth names)',
    'modules', 'purpose', 'type genome (sequence accession)', 'preserved sample (two collections in different places)',
    'genome watermark', 'release (date, place, scale, containment)', 'versions (each edit, derived-from link)',
    'status', 'common names (name, language, community, date, meaning)',
)

# ---- The Earth backbone ------------------------------------------------------------------------------------------
# Ranked nodes as the GBIF backbone places them (gbif_backbone), intermediate ranks omitted; rankless nodes cited.
BACKBONE = (
    dict(name='Life', rank=None, parent=None, code=None, sources=['phylocode_6'],
         note='The clade of all organisms, assuming one origin (PhyloCode Note 2.1.1)'),
    dict(name='Bacteria', rank='domain', parent='Life', code='ICNP', sources=['woese_1990', 'goker_oren_2024']),
    dict(name='Cyanobacteriota', rank='phylum', parent='Bacteria', code='ICNP and ICN',
         sources=['oren_2022_cyanobacteriota', 'icn_madrid_2025'],
         note='Governed by both codes (ICN Pre. 8); the GBIF backbone lists it as Cyanobacteria'),
    dict(name='Eukaryota', rank='domain', parent='Life', code=None, sources=['woese_1990', 'adl_2019'],
         note='Woese et al. named it Eucarya; eukaryote names above the codes\' ranks are not regulated'),
    dict(name='Plantae', rank='kingdom', parent='Eukaryota', code='ICN', sources=['gbif_backbone']),
    dict(name='Chlorophyta', rank='phylum', parent='Plantae', code='ICN', sources=['gbif_backbone']),
    dict(name='Bryophyta', rank='phylum', parent='Plantae', code='ICN', sources=['gbif_backbone']),
    dict(name='Tracheophyta', rank='phylum', parent='Plantae', code='ICN', sources=['gbif_backbone', 'cantino_2007']),
    dict(name='Angiospermae', rank=None, parent='Tracheophyta', code=None, sources=['cantino_2007'],
         note='The flowering plants, a phylogenetically defined clade name'),
    dict(name='Caryophyllales', rank='order', parent='Angiospermae', code='ICN', sources=['gbif_backbone']),
    dict(name='Nepenthaceae', rank='family', parent='Caryophyllales', code='ICN', sources=['gbif_backbone']),
    dict(name='Nepenthes', rank='genus', parent='Nepenthaceae', code='ICN', sources=['gbif_backbone']),
    dict(name='Cucurbitaceae', rank='family', parent='Angiospermae', code='ICN', sources=['gbif_backbone']),
    dict(name='Solanum', rank='genus', parent='Angiospermae', code='ICN', sources=['gbif_backbone']),
    dict(name='Solanum lycopersicum', rank='species', parent='Solanum', code='ICN', sources=['gbif_backbone']),
    dict(name='Fungi', rank='kingdom', parent='Eukaryota', code='ICN', sources=['gbif_backbone']),
    dict(name='Animalia', rank='kingdom', parent='Eukaryota', code='ICZN', sources=['gbif_backbone', 'iczn_1999']),
    dict(name='Arthropoda', rank='phylum', parent='Animalia', code='ICZN', sources=['gbif_backbone']),
    dict(name='Insecta', rank='class', parent='Arthropoda', code='ICZN', sources=['gbif_backbone']),
    dict(name='Diptera', rank='order', parent='Insecta', code='ICZN', sources=['gbif_backbone']),
    dict(name='Arachnida', rank='class', parent='Arthropoda', code='ICZN', sources=['gbif_backbone']),
    dict(name='Araneae', rank='order', parent='Arachnida', code='ICZN', sources=['gbif_backbone']),
    dict(name='Chordata', rank='phylum', parent='Animalia', code='ICZN', sources=['gbif_backbone']),
    dict(name='Aves', rank='class', parent='Chordata', code='ICZN', sources=['gbif_backbone']),
    dict(name='Salmonidae', rank='family', parent='Chordata', code='ICZN', sources=['gbif_backbone'],
         note='The GBIF backbone gives ray-finned fishes no class'),
)


def ref(name, role, lit=None):
    """An Earth lineage drawn on as a reference point, with the literature key where this work read it."""
    return dict(name=name, role=role, lit=list(lit or []))


def src(path, at=''):
    return dict(path=path, at=at)


def entry(id, kind, working_name, summary, status, basis, sources, scope='throughout', tests=(), notes='', **kw):
    return dict(id=id, kind=kind, working_name=working_name, summary=summary, status=status, scope=scope, basis=basis,
                tests=list(tests), notes=notes, sources=list(sources), **kw)


def group(id, working_name, summary, *, anchor, degree, chassis, references, form, settings, realm, guilds, night,
          clock, status, basis, sources, parent=None, form_candidates=(), size=(), modules=(), carries=(), scope='throughout',
          tests=(), notes='', name_candidates=()):
    return entry(id, 'subgroup' if parent else 'group', working_name, summary, status, basis, sources, scope=scope,
                 tests=tests, notes=notes, parent=parent, anchor=anchor, degree=degree, chassis=list(chassis),
                 references=list(references), form=form, form_candidates=list(form_candidates), size=list(size),
                 settings=list(settings), realm=realm, guilds=list(guilds), night=list(night), clock=list(clock),
                 modules=list(modules), carries=list(carries), name_candidates=list(name_candidates))


def module(id, working_name, summary, *, status, basis, sources, references=(), clock=(), night=(), scope='throughout',
           tests=(), notes=''):
    return entry(id, 'module', working_name, summary, status, basis, sources, scope=scope, tests=tests, notes=notes,
                 references=list(references), clock=list(clock), night=list(night))


def community(id, working_name, summary, *, members, settings, realm, status, basis, sources, described=(),
              scope='throughout', tests=(), notes=''):
    return entry(id, 'community', working_name, summary, status, basis, sources, scope=scope, tests=tests, notes=notes,
                 members=list(members), described_members=list(described), settings=list(settings), realm=realm)


# ---- Design modules: designed traits that groups carry ----------------------------------------------------------
MODULES = (
    module('M01', 'carbon-balance idling', 'Leaves, stems and roots idle at about a quarter of their daytime upkeep '
           'as soon as photosynthesis stops covering upkeep, cued by the energy sensors SnRK1 and TOR; this halves '
           'the night store against idling only in full darkness (10 against 19 g C per m2)',
           status='carried_forward', basis='selected', clock=['carbon_cued'], night=['idle'],
           references=[ref('Arabidopsis thaliana', 'SnRK1 mutants mobilize starch poorly at night; TOR and SnRK1 '
                           'switch growth and starvation programs')],
           sources=[src(LCE, 'B1, C3; Selection 1'), src(BIO, 'Canopy photosynthesis')]),
    module('M02', 'dusk-set night timer', 'A timer sized for about 240 dark hours, started at dusk, that rations the '
           'night\'s store; a plant on an Earth clock would budget the lunar night\'s store for one day',
           status='carried_forward', basis='selected', clock=['dusk_set'], night=['store'],
           references=[ref('Arabidopsis thaliana', 'divides its leaf starch by the hours to the dawn its clock expects')],
           sources=[src(LCE, 'C2; Selection 2'), src(ECO, 'finding 20')]),
    module('M03', 'two clocks', 'Daily rhythms in leaf chemistry inside the long lunar clock, with the daily clock '
           'kept off light harvesting', status='carried_forward', basis='selected', clock=['two_clocks'],
           references=[ref('Solanum lycopersicum', 'continuous-light injury comes from the clock turning down the '
                           'light-harvesting gene CAB-13 on schedule')],
           sources=[src(LCE, 'C1, C4; Selection 2'), src(JULY, 'section 2')]),
    module('M04', 'long-dark leaf survival', 'Leaves that survive about 240 dark hours without senescing, with '
           'metabolism that slows and restarts on cue', status='test_first', basis='selected', night=['idle'],
           references=[ref('Arabidopsis thaliana', 'a leaf darkened alone senesces; darkening the whole plant holds '
                           'senescence back')],
           tests=['leaves that survive about 240 dark hours, and metabolism that idles and restarts (B5)'],
           sources=[src(LCE, 'B5; Test first'), src(BIO, 'Canopy photosynthesis')]),
    module('M05', 'twilight use', 'Photosynthesis through the long twilight, worth 14% of an Earth-like night store '
           'and a quarter of the idling one; it leaves 238 dark hours at the equator and 34 at 60 degrees',
           status='carried_forward', basis='selected', night=['twilight'], sources=[src(LCE, 'B2; Selection 1')]),
    module('M06', 'day-fruit program', 'A fruit program one sunlit half long: set at sunrise on a dawn cue without '
           'pollination, cells formed before sunrise, sugar loaded in the last days before dusk, ripe at sunset',
           status='carried_forward', basis='selected', clock=['dawn_set'],
           references=[ref('Citrullus lanatus', 'today\'s watermelon takes 45 days and lives through a night'),
                       ref('Cucurbita maxima', 'giant pumpkins take in about 15 kg a day at their peak'),
                       ref('Cucumis sativus', 'parthenocarpic cucumbers set fruit without pollination')],
           tests=['a large fruit\'s program compressed to one lunar day (D2)'],
           sources=[src(LCE, 'D2-D4; Selection 3; Test first'), src(BIO, 'Canopy photosynthesis')]),
    module('M07', 'lunar calendar of animal lives', 'Emergence, spawning and migration timed by dawn and dusk, the '
           'month\'s sharpest cues: migrants leave the sea in the lunar afternoon, arrive at dusk and lay eggs that '
           'hatch by dawn', status='carried_forward', basis='selected', clock=['dusk_set', 'dawn_set'],
           references=[ref('Clunio marinus', 'emerges at the low tides around new and full moon on circalunar and '
                           'circadian clocks'),
                       ref('Samoan palolo worm', 'casts its spawning segments at the third quarter moon'),
                       ref('Great Barrier Reef corals', 'at least 32 species spawn together after late-spring full '
                           'moons')],
           sources=[src(LCE, 'C5, H7; Selection 2')]),
    module('M08', 'continuous-light tolerance on farms', 'Cultivars that tolerate continuous light, and daily cues '
           'such as a temperature swing or alternating red and blue light', status='carried_forward', basis='selected',
           scope='in_places', clock=['earth_clock_with_cues'],
           references=[ref('Solanum lycopersicum', 'a wild tomato\'s CAB-13 allele gives tolerance with up to 20% '
                           'more yield')],
           sources=[src(LCE, 'C1; In places'), src(ECO, 'finding 20')]),
    module('M09', 'lasting glow', 'Bioluminescence from the fungal pathway that lasts through the first nights after '
           'dusk; a glowing square metre draws about 0.04 W, under 2% of the night upkeep of a square metre of '
           'plants', status='test_first', basis='selected', night=['night_active'], clock=['dusk_set'],
           references=[ref('Neonothopanus gardneri', 'luminous fungus; glowing artificial mushrooms caught 42 insects '
                           'against 12 in dark ones'),
                       ref('Nicotiana tabacum', 'engineered to glow with the fungal pathway; the glow followed its '
                           'clock and faded by the third and fourth day of darkness'),
                       ref('glowworm larvae', 'hang sticky threads beneath their light and catch mostly flies')],
           tests=['glowing lures that last through the first nights, with a catch the migrants can bear (E3)'],
           sources=[src(LCE, 'E3; Test first')]),
    module('M10', 'violet and visible signals', 'Flowers that signal at 350-400 nm or in the violet and visible, and '
           'pollinators tuned to them: a honeybee\'s UV receptor catches 7% of Earth\'s light under the shield',
           status='hypothesis', basis='requirement',
           references=[ref('Apis mellifera', 'UV receptor peaking near 344 nm'),
                       ref('flowering plants', 'floral UV patterning follows the UV a habitat receives')],
           sources=[src(ECO, 'finding 2'), src(LIT, 'Light')]),
    module('M11', 'raised night-light thresholds', 'Circadian and melatonin thresholds raised for nearside nights of '
           '2-3 lux of Earthlight, set zone by zone; far-side low latitudes keep Earth-calibrated thresholds',
           status='hypothesis', basis='requirement',
           references=[ref('fishes and rodents', 'melatonin suppressed at 0.01-0.03 lux of light at night')],
           sources=[src(ECO, 'finding 5'), src(LIT, 'Light')]),
    module('M12', 'vitamin D without UV-B', 'A route to vitamin D for every vertebrate without UV-B, which the shield '
           'removes: biosynthesis that runs without the light, calcium handling without the vitamin, or a dietary '
           'supply', status='hypothesis', basis='requirement',
           references=[ref('cats and dogs', 'mammals that live on dietary vitamin D')],
           sources=[src(ECO, 'finding 1'), src(LIT, 'Light')]),
    module('M13', 'complete denitrification', 'Denitrification that ends in N2 through the N2O reductase NosZ, since '
           'N2O has no sink under the shield', status='hypothesis', basis='requirement',
           sources=[src(ECO, 'finding 7')]),
    module('M14', 'alternative nitrogenases', 'Vanadium and iron-only nitrogenases beside the molybdenum enzyme, since '
           'molybdenum may limit fixation before nitrogen does', status='hypothesis', basis='requirement',
           sources=[src(ECO, 'finding 12'), src(LIT, 'Nitrogen')]),
    module('M15', 'wet-leaf disease resistance', 'Resistance to leaf pathogens through a hundred hours and more of '
           'leaf wetness, with sunlight disinfecting at 1-6% of Earth\'s rate', status='hypothesis',
           basis='requirement', sources=[src(ECO, 'finding 3')]),
    module('M16', 'root regulation at lunar gravity', 'Root growth regulated by design, since plants sense gravity '
           'poorly below 0.1-0.3 g', status='hypothesis', basis='requirement', sources=[src(ECO, 'finding 19')]),
    module('M17', 'ethylene balance', 'Low ethylene emission by plants and fast ethylene uptake by soils: the '
           'vegetated Moon\'s own ethylene settles at 94-233 ppb with Earth\'s soil uptake, above the 50 ppb that '
           'costs wheat 36% of its yield', status='hypothesis', basis='requirement',
           sources=[src(ECO, 'finding 24'), src(PEOPLE, 'item 16')]),
    module('M18', 'gas-barrier envelope', 'A mostly inert, replaceable hydrogen barrier grown and repaired by living '
           'tissue, behind a tougher outer surface', status='hypothesis', basis='idea',
           references=[ref('goldbeater\'s skin', 'processed intestinal collagen of cattle used as hydrogen gas-cell '
                           'barriers in historical airships'),
                       ref('teleost swim bladders', 'guanine plates make the bladder wall a hundred times tighter'),
                       ref('Nereocystis luetkeana', 'the wall of the bull kelp\'s float makes gas as it grows',
                           ['liggan_martone_2020'])],
           tests=['one assembled wet barrier: hydrogen, oxygen and nitrogen permeation, sustained load and creep, seams, '
                  'tear arrest and healing, with and without guanine-like plates'],
           sources=[src(AEROPHYTES, 'The reference organism; What would decide it next'),
                    src(AEROPHYTES_GAS, 'Guanine plates in the bladder wall; Goldbeater\'s skin')]),
    module('M19', 'hydrogen organ or symbiont', 'Hydrogen made by a dedicated organ or symbiont, from light or by dark '
           'fermentation, with its light, substrates and oxygen management counted', status='hypothesis',
           basis='idea',
           references=[ref('Chlamydomonas reinhardtii', 'sulfur-deprived cultures enter a hydrogen-producing state'),
                       ref('Clostridium', 'dark fermentation yields about 2.1 mol H2 per mol glucose'),
                       ref('Physalia physalis', 'secretes its float gas from a dedicated gas gland'),
                       ref('teleost swim bladders', 'a gas gland and rete mirabile secrete gas against high pressure')],
           sources=[src(AEROPHYTES_BIO, 'Hydrogen from light and from sugar'),
                    src(AEROPHYTES_GAS, 'The swim bladder\'s gas gland; The Portuguese man-of-war\'s gas gland')]),
    module('M20', 'water from the air', 'Water drawn from rain, cloud and the humid air through sorbent skins and '
           'fringes, in preference to descending to drink', status='hypothesis', basis='direction',
           sources=[src(DEC, 'Sky giants and the biosphere\'s names'), src(AEROPHYTES, 'The limits, in brief')]),
    module('M21', 'attached buds', 'Daughters grown as attached buds, fed carbon, nutrients and gas by the parent '
           'until they can carry themselves', status='hypothesis', basis='idea',
           sources=[src(AEROPHYTES_BIO, 'Buds, offspring and replacement')]),
    module('M27', 'tethered sails', 'A wing hung on a tether in another layer of air, or a drogue in the sea, that '
           'sails a buoyant body across the wind', status='hypothesis', basis='idea',
           references=[ref('Physalia physalis', 'sails up to 50-55 degrees from downwind with its tentacles as the sea '
                           'anchor', ['iosilevskii_weihs_2009'])],
           sources=[src(AEROPHYTES_SAIL, 'Sailing between two layers; Sailors between air and water')]),
    module('M28', 'compartmented gas', 'Gas held in many small wet compartments or modules, so that a tear or a '
           'lightning strike empties or burns one at a time', status='hypothesis', basis='idea',
           sources=[src(AEROPHYTES, 'The limits, in brief')]),
    module('M22', 'spectral tuning', 'Paler upper leaves, photosynthesis into the far red and clumped foliage; they add '
           '3%, 6% and 3% on the Moon, where the diffuse sky already does most of what they do', status='set_aside',
           basis='selected', sources=[src(LCE, 'A2; Set aside')]),
    module('M23', 'extra photoprotection and daily rests', 'Reflective hairs, folding leaves and a daily rest in the '
           'long day: peak light is 30% below Earth\'s and the UV index 0.11, and a rest halving photosynthesis for '
           'eight hours in 24 costs about a quarter of the growth', status='set_aside', basis='selected',
           sources=[src(LCE, 'A4, A5; Set aside')]),
    module('M24', 'canopy regrowth each cycle', 'Shedding the canopy before dusk and regrowing it at dawn costs '
           '264 g C per m2 a cycle and leaves little or no growth', status='set_aside', basis='selected',
           night=['shedding'], sources=[src(LCE, 'B3; Set aside')]),
    module('M25', 'frost protection', 'Antifreeze, insulation and heat sharing in plant mats: land nights stay near '
           '19-20 degrees C at low and middle latitudes', status='set_aside', basis='selected',
           sources=[src(LCE, 'B4; Set aside')]),
    module('M26', 'hyper-C4 as the general strategy', 'C4 everywhere: the canopy works mostly in dim diffuse light, '
           'where C4\'s extra two ATP per CO2 cost most', status='set_aside', basis='selected',
           sources=[src(LCE, 'A3; Set aside')]),
)

MODULES = tuple(sorted(MODULES, key=lambda m: m['id']))

# ---- Designed groups -------------------------------------------------------------------------------------------
G = []
# The sky
G.append(group('G01', 'aerophytes', 'Large photosynthetic aerial organisms with gas-filled, lobed bodies: a light, '
    'mostly inert and replaceable envelope, distributed photosynthetic tissue and small protected pockets of retained '
    'water that hold their own small communities', anchor='Eukaryota', degree='lineage',
    chassis=['brown algae (Phaeophyceae), whose kelps inflate their own floats',
             'green algae (Chlorophyta), for photosynthetic hydrogen'],
    references=[ref('Physalia physalis', 'a colony of polyps with a gas-filled, sail-like float and trailing '
                    'tentacles, which sails', ['iosilevskii_weihs_2009']),
                ref('Nereocystis luetkeana', 'bull kelp: a gas-filled float holds it upright to compete for light',
                    ['liggan_martone_2018', 'liggan_martone_2020']),
                ref('holopelagic Sargassum', 'brown algae that live afloat and drift for months', ['machado_2024']),
                ref('siphonophores', 'colonies whose zooids specialise as float, feeders, swimmers and reproductive '
                    'members, budded from growth zones'),
                ref('teleost swim bladders', 'secrete gas into a bladder walled with guanine plates'),
                ref('mistletoes', 'draw water and nutrients from a host\'s sap, as the young draw on megaforest crowns'),
                ref('Chlamydomonas reinhardtii', 'photosynthetic hydrogen'),
                ref('Clostridium', 'dark-fermentation hydrogen')],
    form='drifter', form_candidates=['sailor'], settings=['low_sky', 'storm_layer', 'canopy'], realm='sky',
    guilds=['producer', 'receiver'], night=['twilight', 'terminator_following', 'store', 'idle'],
    clock=['undetermined'], modules=['M18', 'M19', 'M20', 'M21', 'M27', 'M28'], status='hypothesis',
    basis='direction',
    tests=['one assembled wet gas barrier', 'a trajectory of water, carbon, gas and trim along a sailing path through '
           'storms, twilight and humid layers'],
    notes='The author chose the common name aerophytes on 10 October 2026; the scientific name waits for the author '
          'among N1-N6. Until then the working name was Sagan & Salpeter\'s "floaters" (1976), their word for buoyant '
          'organisms of Jupiter\'s air. Aloft the mean wind blows east, so aerophytes ride it to shorten their nights, '
          'linger in the evening flow near the ground or sail north and south across the shear. Migrating flyers roost '
          'on them and bring phosphorus up.',
    name_candidates=['N1', 'N2', 'N3', 'N4', 'N5', 'N6'],
    sources=[src(AEROPHYTES, 'What a giant looks like; The limits, in brief; Twilight and the Sun'),
             src(AEROPHYTES_BIO, 'all sections'),
             src(AEROPHYTES_SAIL, 'Riding the eastward wind; Sailing between two layers'),
             src(AEROPHYTES_GAS, 'Siphonophores; The swim bladder\'s gas gland'),
             src(AERIAL, 'Large photosynthetic organisms and their wet mass'),
             src(DEC, 'Cross-domain ecology follow-up; Aerophyte light-history correction; Sky giants')]))
G.append(group('G02', 'lesser aerophytes', 'Aerophytes of the coupled study\'s reference organism: spheres of '
    '40-150 m at about 10 km, the size the giants grow from; the young can grow to 40-60 m in megaforest crowns, '
    'fed by host sap',
    parent='G01', anchor='Eukaryota', degree='lineage', chassis=['as the aerophytes'], references=[], form='drifter',
    form_candidates=['sailor'], settings=['low_sky', 'canopy'], realm='sky', guilds=['producer'],
    night=['twilight', 'terminator_following', 'store', 'idle'], clock=['undetermined'], status='hypothesis',
    basis='idea', sources=[src(AEROPHYTES, 'The reference organism'), src(AEROPHYTES_GROW, 'The canopy nursery')]))
G.append(group('G03', 'giant aerophytes', 'Aerophytes that keep growing for centuries, like the oldest trees: sky '
    'reefs, colonies of wet modules 60-150 m across on a raft one or two kilometres wide, which regrow what storms '
    'take, or round bodies of hundreds of metres to a few kilometres in eight to over a hundred compartments', parent='G01', anchor='Eukaryota',
    degree='lineage',
    chassis=['as the aerophytes'],
    references=[ref('Pinus longaeva', 'the oldest known living eukaryote, 4,770 years old in 2005', ['flanary_2005']),
                ref('trees of 403 species', 'mass growth rate rises continuously with size as leaf area outpaces '
                    'the fall in productivity per leaf area', ['stephenson_2014']),
                ref('siphonophores', 'the colonial plan: specialised modules added at growth zones')],
    form='giant', settings=['low_sky', 'storm_layer'], realm='sky', guilds=['producer', 'receiver'],
    night=['twilight', 'terminator_following', 'store', 'idle'], clock=['undetermined'], status='hypothesis',
    basis='direction', tests=['fittings, seams, gusts and damage, then water and storms'],
    notes='A giant\'s age is the time its own light takes to make its gas: decades to more than a millennium from a '
          '40-m juvenile in the coupled study. Lightning and tears take one module or compartment at a time. The '
          'colony form carries the common name sky reefs, the author\'s choice of 10 October 2026; the round form has '
          'no common name yet.',
    sources=[src(DEC, 'Sky giants and the biosphere\'s names'),
             src(AEROPHYTES, 'What a giant looks like; The limits, in brief'),
             src(AEROPHYTES_GAS, 'Siphonophores: colonies of specialised parts')]))
G.append(group('G04', 'resident aeroplankton phototrophs', 'Cyanobacteria and algae that complete their lives aloft, '
    'in cloud water and in retained water', anchor='Life', degree='species',
    chassis=['cyanobacteria (Cyanobacteriota)', 'eukaryotic microalgae'],
    references=[ref('cloud cyanobacteria and algae', 'photosynthetic algae and cyanobacteria recovered from clouds and '
                    'rain')],
    form='aeroplankton', size=['pico', 'nano'], settings=['low_sky', 'storm_layer'], realm='sky', guilds=['producer'],
    night=['undetermined'], clock=['undetermined'], status='hypothesis', basis='direction',
    tests=['reproduction outrunning removal: at 0.01 a day of added mortality a 3.6-day doubling needs favourable '
           'conditions for 22.5% of the cycle if rain removes cells in 30 days'],
    sources=[src(AERIAL, 'What the evidence supports; Reproduction must outrun removal'), src(AERIAL_SRC),
             src(DEC, 'Cross-domain ecology follow-up')]))
G.append(group('G05', 'aerial cyanobacteria', 'Cyanobacteria of the resident aeroplankton', parent='G04',
    anchor='Cyanobacteriota', degree='species', chassis=['cyanobacteria'], references=[], form='aeroplankton',
    size=['pico', 'nano'], settings=['low_sky', 'storm_layer'], realm='sky', guilds=['producer'],
    night=['undetermined'], clock=['undetermined'], status='hypothesis', basis='direction', sources=[src(AERIAL)]))
G.append(group('G06', 'aerial microalgae', 'Eukaryotic algae of the resident aeroplankton', parent='G04',
    anchor='Eukaryota', degree='species', chassis=['eukaryotic microalgae'], references=[], form='aeroplankton',
    size=['nano'], settings=['low_sky', 'storm_layer'], realm='sky', guilds=['producer'], night=['undetermined'],
    clock=['undetermined'], status='hypothesis', basis='direction', sources=[src(AERIAL)]))
G.append(group('G07', 'resident aeroplankton heterotrophs', 'Bacteria that grow in cloud and fog water on organic '
    'carbon other organisms make', anchor='Bacteria', degree='species', chassis=['cloud-water bacteria'],
    references=[ref('cloud bacteria', 'grow in supercooled cloud droplets and keep transcribing in cloud water'),
                ref('fog heterotrophs', 'C1-consuming heterotrophs growing with fog droplets')],
    form='aeroplankton', size=['pico'], settings=['low_sky', 'storm_layer'], realm='sky', guilds=['decomposer'],
    night=['undetermined'], clock=['undetermined'], status='hypothesis', basis='direction',
    sources=[src(AERIAL, 'What the evidence supports'), src(AERIAL_SRC), src(PEOPLE, 'item 3')]))
G.append(group('G08', 'aerial collectors', 'Small animals that gather airborne pollen, spores and small organisms on '
    'webs or filters, on substrates that take the loads: aerophytes and canopies', anchor='Araneae', degree='species',
    chassis=['spiders (Araneae)'],
    references=[ref('orb-weaving spiders', 'eat the pollen their webs catch'),
                ref('ballooning spiders', 'rise on the charge of their silk in an electric field',
                    ['morley_gorham_2020'])],
    form='tenant', settings=['low_sky', 'canopy'], realm='sky', guilds=['collector'], night=['undetermined'],
    clock=['undetermined'], status='hypothesis', basis='idea',
    sources=[src(AERIAL, 'Filter-feeding has a speed and drag budget'), src(AERIAL_SRC)]))
G.append(group('G09', 'large aerial filter feeders', 'Large mobile animals that filter food from the air; at 10 m/s a '
    'filter breaks even at 0.9-45 mg of dry food per m3 across drag coefficients of 0.02-1', anchor='Animalia',
    degree='lineage', chassis=['undetermined'], references=[], form='flyer', settings=['low_sky', 'storm_layer'],
    realm='sky', guilds=['collector'], night=['undetermined'], clock=['undetermined'], status='hypothesis',
    basis='idea', notes='A demanding hypothesis in the aerial study, which names no Earth model and infers no '
    'whale-sized animal.', sources=[src(AERIAL, 'Filter-feeding has a speed and drag budget')]))
G.append(group('G10', 'high-altitude films', 'Sparse microbial and algal films on particles, droplets and platforms '
    'in the flight band and the high air', anchor='Life', degree='species',
    chassis=['desiccation-tolerant bacteria', 'engineered algae'],
    references=[ref('Deinococcus radiodurans', 'its radiation resistance follows from its tolerance of desiccation',
                    ['mattimore_battista_1996']),
                ref('engineered algal films', 'films attached to research platforms (July 2026 bundle)')],
    form='film', settings=['flight_band', 'high_air'], realm='sky', guilds=['producer'], night=['dormant'],
    clock=['undetermined'], status='carried_forward', basis='decision',
    notes='Cosmic-ray ionization is 7% of Earth\'s sea level in the flight band and 16% at 45-90 km, so the '
          'sparseness answers to cold, dry air and scarce nutrients.',
    sources=[src(DEC, 'Life and people'), src(PEOPLE, 'item 2'), src(JULY, 'section 12')]))
G.append(group('G11', 'great flyers', 'Large soaring animals of the dense lower sky; a flyer six times Earth\'s in '
    'every length carries the same stresses, which carries Earth\'s largest flyers to 60-66 m across and 44-55 t',
    anchor='Aves', degree='genus', chassis=['birds (Aves)'],
    references=[ref('Diomedea exulans', 'wandering albatross, 3.0-3.1 m and 7.8-9.4 kg'),
                ref('Argentavis magnificens', 'Miocene, 7 m and 70-72 kg'),
                ref('Quetzalcoatlus northropi', 'Cretaceous pterosaur, 10-11 m and 200-250 kg')],
    form='flyer', settings=['low_sky', 'storm_layer'], realm='sky', guilds=['consumer'], night=['undetermined'],
    clock=['undetermined'], modules=['M12'], status='hypothesis', basis='idea',
    notes='Extinct flyers serve as reference points; a chassis is a living lineage.',
    sources=[src(PEOPLE, 'item 5, item 22, Table 2'), src(DEC, 'Life and people; Sky ships and flyers'),
             src(JULY, 'section 11')]))
G.append(group('G12', 'gliders', 'Animals that launch from canopies and cliffs and glide, at a sixth of Earth\'s '
    'weight in the dense air', anchor='Animalia', degree='species', chassis=['arboreal vertebrates', 'insects'],
    references=[ref('arboreal vertebrate gliders', 'more than 30 independent lineages', ['dudley_2007']),
                ref('arboreal ants', 'directed aerial descent without wings', ['dudley_2007'])],
    form='glider', settings=['canopy', 'low_sky'], realm='sky', guilds=['consumer'], night=['undetermined'],
    clock=['undetermined'], status='carried_forward', basis='decision',
    sources=[src(DEC, 'Life and people'), src(JULY, 'section 11')]))
G.append(group('G13', 'sea-forest migrant insects', 'Insects that grow at sea and fly to the forests for the dusk '
    'fruit fall and to mate, leaving phosphorus, nitrogen and trace elements in droppings, carcasses and the prey of '
    'lures and traps, then fly back to breed at sea', anchor='Insecta', degree='species',
    chassis=['moths (Lepidoptera)', 'midges (Diptera)'],
    references=[ref('Agrotis infusa', 'Bogong moths fly about 1,000 km to alpine caves and bring nitrogen and '
                    'phosphorus where they die'),
                ref('Chironomidae', 'Lake Myvatn midges fertilise the first 50 m of shore'),
                ref('Clunio marinus', 'emerges in millions on lunar and tidal timing'),
                ref('high-flying insect migrants', 'carry 100 t of nitrogen and 10 t of phosphorus a year over '
                    'southern Britain')],
    form='flyer', settings=['seas', 'low_sky', 'canopy'], realm='sky', guilds=['carrier'], night=['night_active'],
    clock=['dusk_set', 'dawn_set'], modules=['M07', 'M11'], carries=['P', 'N', 'trace elements', 'K'],
    status='carried_forward', basis='selected',
    tests=['the conveyor\'s size: seas producing insects at the top of Myvatn\'s range, spread over hundreds of '
           'kilometres (H4)'],
    notes='Returning the dissolved phosphorus takes 0.5-1.2 g of dry insect bodies on each m2 of land a year; a '
          'Bogong-sized moth flies the 950 km round trip to the median sea distance on a quarter of its mass in sugar.',
    sources=[src(LCE, 'H1-H7; Selection 6; Test first')]))
G.append(group('G14', 'seabird-like colonial fliers', 'Vertebrate fliers that feed at sea and nest in colonies on '
    'land, whose droppings and bones carry nitrogen and phosphorus inland', anchor='Aves', degree='species',
    chassis=['seabirds (Aves)'],
    references=[ref('seabird colonies', 'receive 591 Gg of nitrogen and 99 Gg of phosphorus a year in droppings'),
                ref('Spitsbergen colonies', 'guano at a median of 0.3-0.4 g per m2 a day, dropping off within '
                    '50-300 m')],
    form='flyer', settings=['seas', 'low_sky', 'ground', 'canopy'], realm='sky', guilds=['carrier'],
    night=['undetermined'], clock=['undetermined'], modules=['M12'], carries=['N', 'P'], status='carried_forward',
    basis='selected', sources=[src(LCE, 'H1, H2, H5; Selection 6')]))
G.append(group('G15', 'dusk followers', 'Animals that fly west at dusk\'s own pace, staying at dusk to meet the fruit '
    'fall as it happens', anchor='Animalia', degree='species', chassis=['small birds (Aves)', 'moths (Lepidoptera)'],
    references=[ref('small birds', 'can keep pace with dusk at the equator'),
                ref('moths', 'keep pace only near the poles or with the wind')],
    form='flyer', settings=['low_sky', 'canopy'], realm='sky', guilds=['fruit_eater'],
    night=['terminator_following'], clock=['dusk_set'], status='carried_forward', scope='in_places', basis='selected',
    notes='Dusk crosses the equator at 15.4 km/h (cycle block).', sources=[src(LCE, 'E6, H3; In places')]))
# The land
G.append(group('G16', 'evergreen idlers', 'Evergreen trees and shrubs that keep their leaves through the 354-hour '
    'night, store its carbon in wood and roots, idle on their own carbon balance and use the long twilight',
    anchor='Tracheophyta', degree='species', chassis=['seed plants'],
    references=[ref('Arabidopsis thaliana', 'the starch timer, the energy sensors SnRK1 and TOR, and whole-plant '
                    'darkening'),
                ref('evergreen trees under polar winters', 'shedding a canopy costs more carbon than evergreen '
                    'respiration through polar winters')],
    form='tree', settings=['canopy', 'ground'], realm='land', guilds=['producer'], night=['store', 'idle', 'twilight'],
    clock=['carbon_cued', 'dusk_set', 'two_clocks'], modules=['M01', 'M02', 'M03', 'M04', 'M05', 'M15', 'M16', 'M17'],
    status='carried_forward', basis='selected',
    tests=['leaves that survive about 240 dark hours, and metabolism that idles and restarts (B5)'],
    notes='Idling at a quarter of daytime upkeep cuts the night store from 46 to 10 g C per m2.',
    sources=[src(LCE, 'B1, B2, C3, D1, F2; Selection 1'), src(BIO, 'Canopy photosynthesis'), src(LIT, 'The long day and '
             'night')]))
G.append(group('G17', 'day-fruit plants', 'Plants carrying the day-fruit program: fruit set at sunrise and ripe at '
    'sunset, grown on the daylight surplus without touching the night\'s store', anchor='Angiospermae',
    degree='species', chassis=['flowering plants'],
    references=[ref('Cucurbitaceae', 'watermelon, giant pumpkin and cucumber as the fruit model\'s references')],
    form=None, form_candidates=['vine', 'tree'], settings=['ground', 'canopy'], realm='land', guilds=['producer'],
    night=['store', 'idle'], clock=['dawn_set', 'carbon_cued'], modules=['M06', 'M01', 'M02'],
    status='carried_forward', basis='selected', tests=['a large fruit\'s program compressed to one lunar day (D2)'],
    sources=[src(LCE, 'D2-D4; Selection 3'), src(BIO, 'Canopy photosynthesis')]))
G.append(group('G18', 'day-fruit crops', 'Farm plants putting 50-60% of their new tissue carbon into day fruit: '
    '3.1-3.7 kg per m2 every lunar cycle at the equator, 51% more than the same stand on Earth', parent='G17',
    anchor='Cucurbitaceae', degree='line', chassis=['cucurbits'],
    references=[ref('Citrullus lanatus', 'watermelon'), ref('Cucurbita maxima', 'giant pumpkin'),
                ref('Cucumis sativus', 'cucumber')],
    form='vine', settings=['ground'], realm='land', guilds=['producer'], night=['store', 'idle'],
    clock=['dawn_set', 'carbon_cued'], modules=['M16'], status='carried_forward', scope='in_places', basis='selected',
    sources=[src(LCE, 'D3; Selection 3'), src(BIO, 'Canopy photosynthesis')]))
G.append(group('G19', 'day-fruit trees', 'Forest trees that put about 30% of their new tissue into day fruit, which '
    'falls at dusk and feeds the night: 116 g of sugar per m2 each cycle, about what the plants burn in the dark',
    parent='G17', anchor='Angiospermae', degree='species', chassis=['flowering trees'], references=[], form='tree',
    settings=['canopy'], realm='land', guilds=['producer'], night=['store', 'idle', 'twilight'],
    clock=['dawn_set', 'dusk_set'], status='carried_forward', basis='selected', sources=[src(LCE, 'E1, F2')]))
G.append(group('G20', 'farm crops on Earth programs', 'Earth crops kept on today\'s programs on farms: fruit picked '
    'young within one lunar day, and cultivars that tolerate continuous light given daily cues',
    anchor='Angiospermae', degree='line', chassis=['crop species'], references=[], form=None,
    form_candidates=['vine', 'herb'], settings=['ground'], realm='land', guilds=['producer'], night=['undetermined'],
    clock=['earth_clock_with_cues'], status='carried_forward', scope='in_places', basis='selected',
    sources=[src(LCE, 'C1, D5; In places')]))
G.append(group('G21', 'fruit picked young', 'Crops whose fruit is picked a few days to two weeks after flowering, '
    'inside one lunar day on today\'s programs', parent='G20', anchor='Cucurbitaceae', degree='line',
    chassis=['cucurbits'], references=[ref('Cucurbita pepo', 'zucchini, picked a few days after flowering'),
                                       ref('Cucumis sativus', 'cucumber, picked within about two weeks')],
    form='vine', settings=['ground'], realm='land', guilds=['producer'], night=['undetermined'],
    clock=['undetermined'], status='carried_forward', scope='in_places', basis='selected', sources=[src(LCE, 'D5')]))
G.append(group('G22', 'continuous-light tomatoes', 'Tomatoes carrying a wild tomato\'s tolerance of continuous light, '
    'grown with daily cues', parent='G20', anchor='Solanum lycopersicum', degree='line', chassis=['tomato'],
    references=[ref('Solanum lycopersicum', 'a wild tomato\'s allele of CAB-13 gives tolerance')],
    form='herb', settings=['ground'], realm='land', guilds=['producer'], night=['undetermined'],
    clock=['earth_clock_with_cues'], modules=['M08'], status='carried_forward', scope='in_places', basis='selected',
    sources=[src(LCE, 'C1')]))
G.append(group('G23', 'storage-organ plants', 'Plants of the plains that carry the night in bulbs, corms and tubers, '
    'without wood', anchor='Angiospermae', degree='species', chassis=['bulbous and tuberous flowering plants'],
    references=[ref('Cape geophytes', 'about 2,100 species of bulbs, corms and tubers, most dormant through the dry '
                    'summer')],
    form='geophyte', settings=['ground'], realm='land', guilds=['producer'], night=['store', 'idle'],
    clock=['carbon_cued'], modules=['M01'], status='carried_forward', basis='selected',
    sources=[src(LCE, 'F1; Selection 4')]))
G.append(group('G24', 'dusk-seeding grasses and forbs', 'Perennial grasses and forbs of the plains that set seed at '
    'dusk, followed by grazers and seed-eaters', anchor='Angiospermae', degree='species', chassis=['grasses', 'forbs'],
    references=[ref('perennial grasses and forbs', 'the plains\' herbs')], form='herb', settings=['ground'],
    realm='land', guilds=['producer'], night=['store', 'idle'], clock=['dusk_set'], status='carried_forward',
    basis='selected', sources=[src(LCE, 'F1; Selection 4')]))
G.append(group('G25', 'C4 plants of dry ground', 'C4 plants on dry or hot ground, where their water economy earns '
    'its place', anchor='Angiospermae', degree='species', chassis=['C4 grasses and other C4 lineages'],
    references=[ref('C4 plants', 'the Moon\'s 48.6 Pa of CO2 moves the temperature where a C4 leaf overtakes a C3 '
                    'leaf 3.6 degrees C higher')],
    form='herb', settings=['ground'], realm='land', guilds=['producer'], night=['store', 'idle'],
    clock=['undetermined'], status='carried_forward', scope='in_places', basis='selected',
    sources=[src(LCE, 'A3; In places')]))
G.append(group('G26', 'trap plants', 'Carnivorous plants of wet ground, poor sands and the canopy, waiting at the '
    'forest end of the migrant cycle', anchor='Caryophyllales', degree='species',
    chassis=['pitcher plants (Nepenthes)', 'sundews (Drosera)'],
    references=[ref('Nepenthes mirabilis', 'takes about 62% of its nitrogen from prey'),
                ref('Drosera', 'sundews take about half their nitrogen from prey')],
    form='epiphyte', form_candidates=['herb'], settings=['canopy', 'lakes_and_wetlands', 'ground'], realm='land',
    guilds=['producer', 'predator', 'receiver'], night=['store', 'idle'], clock=['carbon_cued'], modules=['M01'],
    status='carried_forward', scope='in_places', basis='selected', sources=[src(LCE, 'E4; Selection 6; In places')]))
G.append(group('G27', 'droppings-collecting plants', 'Plants that take nitrogen and phosphorus from the droppings of '
    'visiting animals without trapping them, at roosts, lekking sites and feeding places', anchor='Nepenthes',
    degree='species', chassis=['pitcher plants (Nepenthes)'],
    references=[ref('Nepenthes lowii', 'feeds tree shrews and takes 57-100% of its leaf nitrogen from their '
                    'droppings'),
                ref('Nepenthes hemsleyana', 'houses woolly bats and takes 34% from theirs')],
    form='epiphyte', settings=['canopy'], realm='land', guilds=['producer', 'receiver'], night=['store', 'idle'],
    clock=['carbon_cued'], status='carried_forward', basis='selected', sources=[src(LCE, 'E5, H5; Selection 6')]))
G.append(group('G28', 'glowing lures', 'Plants or fungi that glow through the first nights after dusk and draw '
    'insects, to the fruit fall or into traps', anchor='Eukaryota', degree='species',
    chassis=['flowering plants carrying the fungal pathway', 'luminous fungi'],
    references=[ref('Neonothopanus gardneri', 'a luminous fungus; glowing artificial mushrooms caught more insects '
                    'than dark ones'),
                ref('Nicotiana tabacum', 'engineered to glow with the fungal pathway'),
                ref('glowworm larvae', 'catch flies on sticky threads beneath their light')], form=None,
    form_candidates=['herb', 'epiphyte', 'fungal'], settings=['canopy', 'ground'], realm='land',
    guilds=['lure', 'receiver'], night=['night_active'], clock=['dusk_set'], modules=['M09'], status='test_first',
    basis='selected', tests=['glowing lures that last through the first nights, with a catch the migrants can bear '
                             '(E3)'], sources=[src(LCE, 'E3; Test first')]))
G.append(group('G29', 'megaforest emergents', 'Trees 100-500 m tall, screened so far for static wind loads in single '
    'trees and a 49-tree patch; root grafts are an untested hypothesis', anchor='Tracheophyta', degree='species',
    chassis=['seed plants'],
    references=[ref('Sequoia sempervirens', 'coast redwood, the tallest known tree, with a hydraulic maximum of '
                    '122-130 m on Earth', ['koch_2004'])],
    form='tree', size=['mega'], settings=['canopy', 'ground'], realm='land', guilds=['producer'],
    night=['store', 'idle'], clock=['carbon_cued'], status='hypothesis', basis='idea',
    notes='The forest studies use a leaf area index of 22 in air at 288 K; the canopy model puts the optimum near 6.',
    sources=[src(MEGA, 'Static plant mechanics'), src(PATCH, 'Shared geometry and budget controls'),
             src(BIO, 'Remaining biological work'), src(STATE, 'Statements elsewhere that need correcting, item 3')]))
G.append(group('G30', 'basalt pioneers', 'Nitrogen fixers, chemolithotrophs and cyanobacteria that first colonise '
    'fresh basalt and weather it', anchor='Bacteria', degree='species', chassis=['bacteria and cyanobacteria'],
    references=[ref('Icelandic basalt pioneers', 'the first colonists were nitrogen fixers and chemolithotrophs'),
                ref('cyanobacteria', 'release basalt\'s Ca, Mg, Si and K over five times faster than abiotic '
                    'leaching'),
                ref('biological soil crusts', 'fix about 49 Tg of nitrogen a year, nearly half of Earth\'s land '
                    'fixation')],
    form='crust', settings=['ground'], realm='land', guilds=['fixer', 'soil_builder', 'producer'],
    night=['idle', 'dormant'], clock=['undetermined'], modules=['M14'], status='hypothesis', basis='idea',
    sources=[src(ECO, 'finding 13'), src(LIT, 'Ecosystems on fresh basalt')]))
G.append(group('G31', 'crust mosses', 'Moss mats that hold fresh ground and build soil', anchor='Bryophyta',
    degree='species', chassis=['mosses'],
    references=[ref('Hekla moss mats', 'a moss mat held Hekla\'s succession for 170-700 years')],
    form='crust', settings=['ground'], realm='land', guilds=['producer', 'soil_builder'], night=['idle'],
    clock=['undetermined'], status='hypothesis', basis='idea',
    notes='The order of introduction decides which state a region reaches.',
    sources=[src(ECO, 'finding 13'), src(LIT, 'Ecosystems on fresh basalt')]))
G.append(group('G32', 'night decomposers of the fruit fall', 'Yeasts and fungi that ferment and decompose the dusk '
    'fruit fall through the warm night', anchor='Fungi', degree='species', chassis=['yeasts', 'mycelial fungi'],
    references=[ref('fruit yeasts', 'their volatiles draw fruit flies, which need the yeast to develop')],
    form='fungal', settings=['ground'], realm='land', guilds=['decomposer'], night=['night_active'],
    clock=['undetermined'], status='carried_forward', basis='selected', sources=[src(LCE, 'E1, F2; Selection 5')]))
G.append(group('G33', 'night insects of the fruit fall', 'Insects drawn to the fermenting fruit fall in the first '
    'nights after dusk', anchor='Insecta', degree='species', chassis=['flies (Diptera)'],
    references=[ref('fruit flies', 'drawn to fermenting fruit by its yeasts\' volatiles')],
    form='flyer', settings=['ground', 'canopy'], realm='land', guilds=['fruit_eater', 'decomposer'],
    night=['night_active'], clock=['dusk_set'], modules=['M11'], status='carried_forward', basis='selected',
    sources=[src(LCE, 'E1, F2; Selection 5')]))
G.append(group('G34', 'fruit-eaters', 'Animals that eat the dusk fruit fall', anchor='Animalia', degree='species',
    chassis=['undetermined'], references=[], form=None, form_candidates=['ground_animal', 'flyer'],
    settings=['canopy', 'ground'], realm='land', guilds=['fruit_eater'], night=['night_active'],
    clock=['undetermined'], modules=['M11'], status='carried_forward', basis='selected', sources=[src(LCE, 'F2')]))
G.append(group('G35', 'grazers and seed-eaters of the plains', 'Animals that follow the plains\' dusk seeding',
    anchor='Animalia', degree='species', chassis=['undetermined'], references=[], form='ground_animal',
    settings=['ground'], realm='land', guilds=['grazer', 'seed_eater'], night=['undetermined'],
    clock=['undetermined'], status='carried_forward', basis='selected', sources=[src(LCE, 'F1')]))
G.append(group('G36', 'pollinators of the violet sky', 'Pollinators whose vision and flowers signal at 350-400 nm or '
    'in the violet and visible, serving broad landscapes', anchor='Insecta', degree='species', chassis=['bees'],
    references=[ref('Apis mellifera', 'UV receptor near 344 nm catches 7% of Earth\'s light under the shield'),
                ref('insects that orient by polarized moonlight', 'navigation on the light left at night')],
    form='flyer', settings=['canopy', 'ground', 'low_sky'], realm='land', guilds=['pollinator'],
    night=['undetermined'], clock=['undetermined'], modules=['M10', 'M11'], status='hypothesis', basis='requirement',
    sources=[src(ECO, 'finding 2'), src(LIT, 'Light'), src(JULY, 'section 11')]))
G.append(group('G37', 'lake insects', 'Insects whose larvae grow in the lakes and whose adults carry nutrients onto '
    'the shore', anchor='Diptera', degree='species', chassis=['midges (Chironomidae)'],
    references=[ref('Chironomidae', 'Lake Myvatn midges put 1 kg of phosphorus per hectare on the first 50 m of shore '
                    'in a high year')],
    form='flyer', settings=['lakes_and_wetlands', 'ground'], realm='land', guilds=['carrier'], night=['undetermined'],
    clock=['undetermined'], carries=['P', 'N'], status='carried_forward', scope='in_places', basis='selected',
    sources=[src(LCE, 'F3, H2; In places')]))
# The waters
G.append(group('G38', 'sea phytoplankton', 'The phytoplankton of the living seas by kind of water: the nutrient-poor '
    'open sea, productive coasts and river mouths, greenest at dusk after the 15-day day', anchor='Life',
    degree='species', chassis=['cyanobacteria', 'eukaryotic algae'],
    references=[ref('Earth\'s ocean phytoplankton', 'chlorophyll 0.02-0.3 mg/m3 in open oceans and 1-10 mg/m3 on '
                    'shelves and upwelling coasts, the design guesses\' basis')],
    form='plankton', size=['pico', 'nano', 'micro'], settings=['seas'], realm='water', guilds=['producer'],
    night=['undetermined'], clock=['undetermined'], status='carried_forward', basis='decision',
    notes='Design guesses: 0.1, 2 and 5 mg/m3 of chlorophyll, twice as much at dusk as at dawn. Calcifiers cannot '
          'build until weathering supplies alkalinity if the delivered water is fresh.',
    sources=[src(DEC, 'Life and people'), src(WATERS), src(SEA, 'The waters\' colour'), src(ECO, 'finding 17')]))
G.append(group('G39', 'dark-night luminous plankton', 'Plankton that glow in the seas where the night grows dark',
    anchor='Life', degree='species', chassis=['undetermined'],
    references=[ref('luminous marine organisms', 'bioluminescence evolved many times in the sea, from bacteria to '
                    'fish', ['haddock_2010'])],
    form='plankton', settings=['seas'], realm='water', guilds=['undetermined'], night=['night_active'],
    clock=['undetermined'], status='hypothesis', basis='decision',
    notes='The decision admits bioluminescence where the lunar-cycle ecology supports it; no value is computed yet.',
    sources=[src(DEC, 'Life and people'), src(SEA, 'Decisions')]))
G.append(group('G40', 'sulfur- and iodine-releasing seaweeds', 'Seaweeds and plankton that release dimethyl sulfide '
    'and iodine, returning sulfur and iodine to the land through the air', anchor='Eukaryota', degree='species',
    chassis=['brown and red seaweeds', 'DMS-producing plankton'],
    references=[ref('Earth\'s sea plankton and seaweeds', 'the ocean releases 28 (18-34) Tg of sulfur a year as '
                    'dimethyl sulfide, and iodine as methyl iodide and inorganic iodine')],
    form='macroalga', settings=['seas'], realm='water', guilds=['producer', 'carrier'], night=['undetermined'],
    clock=['undetermined'], carries=['S', 'I'], status='carried_forward', basis='selected',
    notes='Seawater at 1-10 nM of dimethyl sulfide would hold the air near 2-17 ppb, around the human odour threshold.',
    sources=[src(LCE, 'H8; Selection 6'), src(ECO, 'finding 8')]))
G.append(group('G41', 'sea grazers', 'Zooplankton that graze the phytoplankton through the 15-day night, leaving '
    'about half the biomass at dawn', anchor='Eukaryota', degree='species', chassis=['zooplankton'], references=[],
    form='plankton', settings=['seas'], realm='water', guilds=['grazer'], night=['night_active'],
    clock=['undetermined'], status='hypothesis', basis='idea', sources=[src(WATERS)]))
G.append(group('G42', 'run fish', 'Fish that grow at sea and run up the rivers to spawn and die, carrying the sea\'s '
    'nitrogen and phosphorus inland', anchor='Salmonidae', degree='species', chassis=['salmon (Salmonidae)'],
    references=[ref('salmon', 'trees along salmon streams take 22-24% of their leaf nitrogen from the fish')],
    form='nekton', settings=['seas', 'rivers'], realm='water', guilds=['carrier'], night=['undetermined'],
    clock=['undetermined'], modules=['M12'], carries=['N', 'P'], status='carried_forward', basis='selected',
    sources=[src(LCE, 'H1, H2, H5; Selection 6')]))
# Microbes of soils and waters
G.append(group('G43', 'keepers of the air', 'Microbes of soils and waters that remove the gases the air cannot: '
    'under the shield the air holds no OH, so methane, N2O, hydrogen, carbon monoxide and ethylene leave through '
    'soils, waters and rain', anchor='Bacteria', degree='species', chassis=['soil and water bacteria'],
    references=[], form='microbe', settings=['ground', 'lakes_and_wetlands', 'seas'], realm='land',
    guilds=['air_keeper'], night=['undetermined'], clock=['undetermined'], status='hypothesis', basis='requirement',
    sources=[src(ECO, 'findings 6-9'), src(STATE, 'What other lines assume of ecology')]))
G.append(group('G44', 'methane oxidizers', 'Methane oxidizers in soils and waters, a requirement from the start: '
    'with soils the only sink, methane settles near 20 ppm at central rates', parent='G43', anchor='Bacteria',
    degree='species', chassis=['methanotrophic bacteria'],
    references=[ref('Earth\'s soils', 'take up 30 (11-49) Tg of methane a year')], form='microbe',
    settings=['ground', 'lakes_and_wetlands', 'seas'], realm='land', guilds=['air_keeper'], night=['undetermined'],
    clock=['undetermined'], status='hypothesis', basis='requirement',
    sources=[src(ECO, 'finding 6'), src(LIT, 'Trace gases without OH')]))
G.append(group('G45', 'complete denitrifiers', 'Denitrifiers whose pathway ends in N2, since N2O has no sink under the '
    'shield', parent='G43', anchor='Bacteria', degree='species', chassis=['denitrifying bacteria'], references=[],
    form='microbe', settings=['ground', 'lakes_and_wetlands', 'seas'], realm='land', guilds=['air_keeper'],
    night=['undetermined'], clock=['undetermined'], modules=['M13'], status='hypothesis', basis='requirement',
    sources=[src(ECO, 'finding 7')]))
G.append(group('G46', 'trace-gas consumers of soils', 'Soil microbes that take up hydrogen, carbon monoxide and '
    'ethylene', parent='G43', anchor='Bacteria', degree='species', chassis=['soil bacteria'],
    references=[ref('soil bacteria of 51 phyla', 'carry high-affinity hydrogenases and remove about three-quarters of '
                    'Earth\'s atmospheric hydrogen'),
                ref('soils', 'take 10-25% of the global carbon monoxide flux')],
    form='microbe', settings=['ground'], realm='land', guilds=['air_keeper'], night=['undetermined'],
    clock=['undetermined'], modules=['M17'], status='hypothesis', basis='requirement',
    sources=[src(LIT, 'Trace gases without OH'), src(PEOPLE, 'item 16')]))
G.append(group('G47', 'nitrogen fixers', 'Free-living and symbiotic fixers that supply almost all the biomes\' new '
    'nitrogen, since lightning gives at most 0.2%', anchor='Bacteria', degree='species',
    chassis=['diazotrophic bacteria', 'cyanobacteria'],
    references=[ref('free-living fixers', 'molybdenum limits fixation in weathered tropical soils; fixers grow down to '
                    '5 mbar of N2')],
    form='microbe', settings=['ground', 'seas'], realm='land', guilds=['fixer'], night=['undetermined'],
    clock=['undetermined'], modules=['M14'], status='hypothesis', basis='requirement',
    sources=[src(ECO, 'finding 12'), src(LIT, 'Nitrogen')]))
GROUPS = tuple(G)

# ---- Communities -------------------------------------------------------------------------------------------------
COMMUNITIES = (
    community('C01', 'the sky archipelago', 'Aerophytes and the small communities in their retained-water pockets, '
              'drifting through the low sky, with the migrating flyers that roost on them',
              members=['G01', 'G02', 'G03', 'G04', 'G07', 'G08'],
              described=['microbes, algae, small grazers and decomposers in retained-water pockets',
                         'hydrogen-oxidising microbes lining the gas cells',
                         'migrating flyers that roost on aerophytes and bring phosphorus up'],
              settings=['low_sky', 'storm_layer'], realm='sky', status='hypothesis', basis='direction',
              sources=[src(AEROPHYTES, 'What a giant looks like; The limits, in brief'),
                       src(AERIAL, 'The cycle and the community')]),
    community('C02', 'cloud-water aeroplankton', 'Resident phototrophs and heterotrophs of cloud and fog water, with '
              'the transient pollen, spores and small animals the air carries', members=['G04', 'G05', 'G06', 'G07'],
              settings=['low_sky', 'storm_layer'], realm='sky', status='hypothesis', basis='direction',
              sources=[src(AERIAL, 'What the evidence supports'), src(PEOPLE, 'item 3')]),
    community('C03', 'the films of the high air', 'Sparse films in the flight band and the high air',
              members=['G10'], settings=['flight_band', 'high_air'], realm='sky', status='carried_forward',
              basis='decision', sources=[src(DEC, 'Life and people')]),
    community('C04', 'forests of wood and fruit', 'Evergreen forests whose wood holds the night\'s store and whose '
              'fruit feeds the night, with their trap plants, lures, fruit-eaters, decomposers and night insects',
              members=['G16', 'G19', 'G26', 'G27', 'G28', 'G32', 'G33', 'G34'], settings=['canopy', 'ground'],
              realm='land', status='carried_forward', basis='selected', sources=[src(LCE, 'F2; Selection 4')]),
    community('C05', 'plains of storage organs', 'Plains of bulbs, corms and tubers, with grasses and forbs seeding '
              'at dusk and the grazers and seed-eaters that follow', members=['G23', 'G24', 'G35'],
              settings=['ground'], realm='land', status='carried_forward', basis='selected',
              sources=[src(LCE, 'F1; Selection 4')]),
    community('C06', 'wet margins', 'Lake shores and wetlands, where carnivorous plants and lake insects meet',
              members=['G26', 'G37'], settings=['lakes_and_wetlands', 'ground'], realm='land',
              status='carried_forward', scope='in_places', basis='selected', sources=[src(LCE, 'F3; Selection 4')]),
    community('C07', 'the dusk fruit fall and the night food web', 'The dusk fruit fall feeding yeasts, fungi, '
              'insects, fruit-eaters and arriving migrants through the warm night',
              members=['G19', 'G32', 'G33', 'G34', 'G28', 'G13', 'G15'], settings=['canopy', 'ground'], realm='land',
              status='carried_forward', basis='selected', tests=['CO2 gathering under still night air near the ground'],
              sources=[src(LCE, 'E1, E2; Selection 5')]),
    community('C08', 'the sea-forest return path', 'Carriers grown at sea that land in the forests, and the plants '
              'that wait for them: biology\'s share of the rock cycle\'s missing return limb',
              members=['G13', 'G14', 'G42', 'G26', 'G27', 'G28', 'G40'], settings=['seas', 'canopy', 'rivers'],
              realm='land', status='carried_forward', basis='selected',
              tests=['the conveyor\'s size (H4)', 'a land, lake and sea budget of every rock-derived element'],
              sources=[src(LCE, 'G1-G4, H1-H10; Selection 6')]),
    community('C09', 'farms of the day fruit', 'Farms of day-fruit crops and Earth crops on today\'s programs',
              members=['G18', 'G21', 'G22'], settings=['ground'], realm='land', status='carried_forward',
              scope='in_places', basis='selected', sources=[src(LCE, 'C1, D3, D5'), src(PEOPLE, 'item 11')]),
    community('C10', 'megaforests', 'Forests of 100-500 m trees, whose crowns can nurse young aerophytes',
              members=['G29', 'G02'], described=['young aerophytes of 40-60 m fed by host sap in the crowns'],
              settings=['canopy', 'ground'], realm='land', status='hypothesis', basis='idea',
              sources=[src(MEGA), src(PATCH), src(BIO), src(AEROPHYTES_GROW, 'The canopy nursery')]),
    community('C11', 'pioneer crusts', 'Crusts and mats that colonise fresh basalt and start its soils',
              members=['G30', 'G31'], settings=['ground'], realm='land', status='hypothesis', basis='idea',
              sources=[src(ECO, 'finding 13')]),
    community('C12', 'the living seas', 'Phytoplankton by kind of water, their grazers, the seaweeds that feed the '
              'air its sulfur and iodine, luminous plankton of dark nights and the carriers grown at sea',
              members=['G38', 'G39', 'G40', 'G41', 'G42', 'G13'], settings=['seas'], realm='water',
              status='carried_forward', basis='decision',
              sources=[src(DEC, 'Life and people'), src(WATERS), src(SEA, 'Decisions')]),
)


# ---- Candidate names for the aerophytes ---------------------------------------------------------------------------
def name(id, group_name, giants, roots, meaning, sound, image, logic, place, checks, sources, young=None,
         formation=None):
    return dict(id=id, group=group_name, giants=giants, young=young, formation=formation, roots=roots,
                meaning=meaning, sound=sound, image=image, logic=logic, place=place, checks=checks, sources=sources)


def root(form, language, meaning, source):
    return dict(form=form, language=language, meaning=meaning, source=source)


NAME_CANDIDATES = (
    name('N1', dict(english='welkin', plural='welkins', scientific='Welkinia'),
         dict(english='elder welkin', plural='elder welkins', scientific=None),
         [root('welkin', 'English', 'the vault of the sky; poetic', 'wiktionary'),
          root('welkne', 'Middle English', 'weather; the heavens; earlier, cloud', 'wiktionary'),
          root('wolcen, wolcn', 'Old English', 'a cloud; the clouds; the heavens; the sky', 'bosworth_toller'),
          root('wolcen-faru', 'Old English', 'the cloud-host; the moving clouds', 'bosworth_toller')],
         'The clouds, and then the sky they fill',
         'Two soft syllables, WEL-kin, an old word that ends on "kin"',
         'The dome of the sky in old poetry; living pieces of that sky drifting over the seas',
         'The word meant cloud before it meant sky, and these organisms are the sky\'s living clouds',
         'English working name for the aerophytes (G01); Welkinia the clade name of the new lineage; elder welkins the '
         'giant form (G03); the wolcenfaru, the cloud-host, available for their drifting formation',
         dict(homonyms='Welkinia: no genus in IRMNG or the GBIF backbone (9 October 2026)',
              notes='Rare in modern English, so it arrives without a biological meaning'),
         ['bosworth_toller', 'wiktionary', 'irmng', 'gbif_backbone'],
         formation=dict(english='the wolcenfaru', meaning='the cloud-host, the moving clouds')),
    name('N2', dict(english='holm', plural='holms', scientific='Aeronesia'),
         dict(english='great holm', plural='great holms', scientific=None),
         [root('holm', 'Old English', 'a mound or hill; wave, ocean, sea; land rising from the water, an island in a '
               'river', 'bosworth_toller'),
          root('holmr', 'Old Norse', 'islet', 'wiktionary'),
          root('aēr + nēsos (ἀήρ + νῆσος)', 'Greek', 'air + island, formed like Polynesia', 'lsj_1940')],
         'Islets, here islands of life rising in the air',
         'One round syllable, "holm"; Aeronesia flows, air-oh-NEE-zha',
         'A drifting archipelago of green islets over the sea, each holding its own small pools and residents',
         'The coupled study\'s concept is a light host carrying small retained-water communities: an island',
         'English working name for the aerophytes (G01); Aeronesia the clade name; great holms the giants (G03); the '
         'sky archipelago (C01) their formation',
         dict(homonyms='Aeronesia: no genus in IRMNG or GBIF; Holmia, the obvious Latin form, is preoccupied by a '
                       'trilobite genus and others, so the clade takes the Greek form',
              notes='"Holm" also names the holm oak in English; context separates them'),
         ['bosworth_toller', 'wiktionary', 'lsj_1940', 'irmng', 'gbif_backbone']),
    name('N3', dict(english='megha', plural='meghas', scientific='Meghaphyta'),
         dict(english='meghaduta', plural='meghadutas', scientific='Meghaduta'),
         [root('megha (मेघ)', 'Sanskrit', '"sprinkler", a cloud; also a mass, a multitude', 'monier_williams_1899'),
          root('meghadūta (मेघदूत)', 'Sanskrit', 'cloud-messenger: Kālidāsa\'s poem of an exiled lover who sends his '
               'wife a message by a cloud', 'eb1911_kalidasa')],
         'Clouds; the giants as cloud messengers',
         'MAY-gha, with a breathed g; me-gha-DOO-ta flows in four syllables; to English ears "megha" also carries '
         '"mega"',
         'Kālidāsa\'s monsoon cloud carrying a lover\'s message across the land; giants carrying their communities '
         'across the Moon',
         'A cloud that bears water and messages; the root means one who sprinkles, which suits organisms that gather '
         'and hold water',
         'Meghas the working name for the aerophytes (G01) with Meghaphyta as the clade; meghadutas the giants (G03), '
         'with Meghaduta available as a genus if the giants become a sister lineage',
         dict(homonyms='Megha, Meghaphyta and Meghaduta: no genus in IRMNG or GBIF',
              notes='The macron of meghadūta drops under the codes\' alphabet rule; the poem is secular'),
         ['monier_williams_1899', 'eb1911_kalidasa', 'irmng', 'gbif_backbone']),
    name('N4', dict(english='kunpeng', plural='kunpeng', scientific='Kunpengia'),
         dict(english='peng', plural='peng', scientific=None),
         [root('kūn (鯤)', 'Classical Chinese', 'the fish of the northern darkness in Zhuangzi; in the Erya, fish roe, '
               'the young of all fish', 'erya_zhushu'),
          root('péng (鵬)', 'Classical Chinese', 'the bird Kun becomes, whose wings are like clouds hanging from the '
               'sky', 'zhuangzi_1')],
         'The smallest thing becoming the greatest',
         'Two short, bell-like syllables, koon and pung',
         'Zhuangzi\'s Peng: its back thousands of li across, rising ninety thousand li on deep-piled wind and, in one old '
         'reading, resting only after six months',
         'From a speck-sized propagule to a centuries-old giant; Zhuangzi says wind piled too thin cannot bear great '
         'wings, and the Moon\'s dense, tall air is wind piled deep',
         'Kun the young (propagules and juveniles), peng the giants (G03), kunpeng the lineage (G01) with Kunpengia '
         'as its clade name',
         dict(homonyms='Kunpengia: no genus in IRMNG or GBIF. Peng (Lu & Li, 2023) and Pengia (Sharpe, 1883; Geyer & '
                       'Corbacho, 2015) are animal genera, so the giants\' scientific form needs another word; '
                       'Kunpengopterus (Wang et al., 2010), a pterosaur, already uses the root',
              notes='Kunpeng is a familiar compound in Chinese and a trade name; "peng" is also English slang'),
         ['zhuangzi_1', 'erya_zhushu', 'irmng', 'gbif_backbone'],
         young=dict(english='kun', plural='kun')),
    name('N5', dict(english='serene', plural='serenes', scientific='Sereniphyta'),
         dict(english='serenissima', plural='serenissimas', scientific='Serenissima'),
         [root('serenus', 'Latin', 'clear, fair, bright, of the sky: caelo sereno', 'lewis_short_1879'),
          root('Mare Serenitatis', 'IAU lunar name', 'Sea of Serenity, one of the seas the Near-side Sea joins',
               'iau_gazetteer'),
          root('serene', 'English', 'poetic noun, "the serene of heaven"; and serene or serein, a fine rain from a '
               'cloudless sky after sunset', 'wiktionary')],
         'The clear, calm sky; the giants as the most serene',
         'se-REEN, smooth and long; se-re-NIS-si-ma for the giants, five open syllables',
         'Calm daylight sky over the Sea of Serenity, and a fine rain falling from clear air after sunset',
         'The Moon\'s own sea name; organisms that keep to clear daylight and drink from the air at dusk',
         'Serenes the working name for the aerophytes (G01) with Sereniphyta as the clade; serenissimas the giants '
         '(G03), with Serenissima available as a genus',
         dict(homonyms='Sereniphyta and Serenissima: no genus in IRMNG or GBIF; Serenitas (Wells, 2009) is '
                       'preoccupied; Serenia sounds too close to Sirenia, the sea cows',
              notes='Serenissima was the old title of Venice'),
         ['lewis_short_1879', 'iau_gazetteer', 'wiktionary', 'irmng', 'gbif_backbone']),
    name('N6', dict(english='wingu', plural='mawingu', scientific='Mawingua'),
         dict(english='great wingu', plural='great mawingu', scientific=None),
         [root('wingu, plural mawingu', 'Swahili', 'cloud, clouds', 'wiktionary')],
         'Clouds',
         'ma-WEEN-goo, three open syllables stressed on the second; the singular wingu',
         'Tall clouds over East Africa\'s plains and lakes',
         'A word from one of Africa\'s most widely spoken languages for organisms that cross every nation\'s sky',
         'Mawingu the working name for the aerophytes (G01), a wingu one of them, with Mawingua as the clade; great '
         'mawingu the giants (G03)',
         dict(homonyms='Mawingua: no genus in IRMNG or GBIF',
              notes='A living language: adopted with its speakers, as principle P6 asks'),
         ['wiktionary', 'irmng', 'gbif_backbone']),
)


# ---- Computed blocks -----------------------------------------------------------------------------------------------
def cycle():
    """The lunar cycle, and the speed at which dusk and the Sun cross the Moon by latitude."""
    hours = SYNODIC_MONTH_DAYS * 24.0
    seconds = hours * 3600.0

    def speed(lat_deg, height_m=0.0):
        return 2.0 * math.pi * (MOON_RADIUS + height_m) * math.cos(math.radians(lat_deg)) / seconds

    lats = (0, 30, 60, 80)
    return dict(
        synodic_hours=round(hours, 2), half_cycle_hours=round(hours / 2.0, 2),
        dusk_ground_speed_km_h={str(lat): round(speed(lat) * 3.6, 2) for lat in lats},
        dusk_ground_speed_m_s={str(lat): round(speed(lat), 3) for lat in lats},
        sun_following_speed_m_s_at_10_km={str(lat): round(speed(lat, 10e3), 3) for lat in lats},
        note='Westward speed that keeps local solar time: 2 pi (R + h) cos(latitude) / synodic month. Dusk followers '
             'fly it at the ground; an aerophyte holding the Sun would have to make it good aloft, where the mean wind '
             'blows east.')


def size_classes(screen):
    """Aeroplankton size classes with the first screen's particles placed in them; Raunkiaer's tree heights."""
    particles = screen['settling']['particles']
    rows = []
    for c in AEROPLANKTON_CLASSES:
        lo, hi = c['range_um']
        examples = {}
        for key, p in particles.items():
            size = float(re.search(r'_(\d+(?:\.\d+)?)um$', key).group(1))
            if lo <= size < hi:
                examples[key] = dict(diameter_um=size, fall_m_s_moon=p['fall_m_s']['moon'],
                                     fall_m_s_earth=p['fall_m_s']['earth'],
                                     residence_days_moon_5_km=p['residence_days']['moon_wet_land_day'],
                                     residence_days_moon_12_km=p['residence_days']['moon_cm1_noon'],
                                     residence_days_earth_1_5_km=p['residence_days']['earth_day'])
        rows.append(dict(name=c['name'], range_um=c['range_um'], examples=examples))
    return dict(aeroplankton=rows, trees=[dict(c) for c in TREE_HEIGHT_CLASSES],
                note='Residence is dry settling out of a mixed layer of the stated depth; rain, grazing and death '
                     'remove organisms much sooner. Classes above micro are not used yet.')


# ---- Product ---------------------------------------------------------------------------------------------------
def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def read_inputs():
    """The files this module reads to compute or check the product: hash-bound in producer.inputs."""
    return sorted(str(p.relative_to(ROOT)) for p in (SOURCES, SCREEN))


def cited_files():
    """Every other project file an entry cites, recorded by path in producer.cited_files without a hash."""
    paths = {s['path'] for e in (*MODULES, *GROUPS, *COMMUNITIES) for s in e['sources']}
    return sorted(paths - set(read_inputs()))


def load_screen():
    screen = json.loads(SCREEN.read_text(encoding='utf-8'))
    if screen.get('schema') != SCREEN_SCHEMA:
        raise ValueError(f'{SCREEN}: schema {screen.get("schema")!r}, expected {SCREEN_SCHEMA!r}')
    return screen


def build():
    screen = load_screen()
    return dict(
        schema=SCHEMA,
        producer=dict(domain='biosphere', files={'biosphere/ecology/taxa.py': digest(__file__)},
                      inputs={p: digest(ROOT / p) for p in read_inputs()},
                      cited_files=cited_files(),
                      constants=constants_used(['biosphere/ecology/taxa.py'])),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        system=dict(
            version=VERSION, layers=list(LAYERS), principles=list(PRINCIPLES), choices=list(CHOICES),
            statuses=STATUSES, status_order=list(STATUS_ORDER), scopes=list(SCOPES), bases=BASES, ranks=list(RANKS),
            degrees=DEGREES,
            founding_register=dict(fields=list(FOUNDING_FIELDS), lineages=[],
                                   note='No lineage has been designed or released; the register holds the fields.'),
            naming=dict(rules=list(NAMING_RULES),
                        lunar_roots=[dict(name=n, origin=o) for n, o in LUNAR_ROOTS]),
            functional=dict(realms=REALMS, settings=SETTINGS, form_classes=FORM_CLASSES, forms=FORMS,
                            guilds=GUILDS, night=NIGHT, clocks=CLOCKS, elements=list(ELEMENTS),
                            size_classes=size_classes(screen)),
            cycle=cycle()),
        backbone=[dict(b) for b in BACKBONE],
        modules=list(MODULES), groups=list(GROUPS), communities=list(COMMUNITIES),
        aerophyte_names=dict(
            group='G01', common_name='aerophytes', common_names_chosen='2026-10-10',
            giant_common_names=dict(group='G03', colony_form='sky reefs', round_form=None),
            scientific_name=None, candidates=list(NAME_CANDIDATES),
            note='The author chose the common names on 10 October 2026: aerophytes for the class and sky reefs for the '
                 'colony giants; the round giants have no common name yet. The scientific name waits for the author, '
                 'and N1-N6 remain candidates for it.'))


# ---- Validation ----------------------------------------------------------------------------------------------------
REQUIRED = dict(
    module=('id', 'kind', 'working_name', 'summary', 'status', 'scope', 'basis', 'tests', 'notes', 'sources',
            'references', 'clock', 'night'),
    group=('id', 'kind', 'working_name', 'summary', 'status', 'scope', 'basis', 'tests', 'notes', 'sources', 'parent',
           'anchor', 'degree', 'chassis', 'references', 'form', 'form_candidates', 'size', 'settings', 'realm',
           'guilds', 'night', 'clock', 'modules', 'carries', 'name_candidates'),
    community=('id', 'kind', 'working_name', 'summary', 'status', 'scope', 'basis', 'tests', 'notes', 'sources',
               'members', 'described_members', 'settings', 'realm'),
)
ID_PATTERN = dict(module=r'M\d\d', group=r'G\d\d', subgroup=r'G\d\d', community=r'C\d\d')


def source_keys():
    return {s['key'] for s in json.loads(SOURCES.read_text(encoding='utf-8'))['sources']}


def backbone_ancestors(nodes, name):
    """Names from a backbone node up to the root, the node included; raises on a cycle or a missing parent."""
    chain, seen = [], set()
    while name is not None:
        if name in seen:
            raise ValueError(f'backbone cycle at {name}')
        if name not in nodes:
            raise ValueError(f'backbone node {name} missing')
        seen.add(name)
        chain.append(name)
        name = nodes[name]['parent']
    return chain


def check(product) -> list[str]:
    """Problems with the register, as messages; an empty list means it is consistent."""
    errors = []
    err = errors.append
    keys = source_keys()
    system = product['system']
    fn = system['functional']
    for top in ('schema', 'producer', 'evidence', 'reading_rule', 'system', 'backbone', 'modules', 'groups',
                'communities', 'aerophyte_names'):
        if top not in product:
            err(f'missing top-level key {top}')
    if product.get('schema') != SCHEMA:
        err('schema mismatch')
    producer = product.get('producer', {})
    if sorted(producer.get('inputs', {})) != read_inputs():
        err(f'producer.inputs should hold exactly the files this module reads: {read_inputs()}')
    for path in producer.get('cited_files', []):
        if Path(path).is_absolute() or '..' in Path(path).parts or not (ROOT / path).is_file():
            err(f'producer.cited_files: {path} is not a file in the repository')

    # Vocabularies hang together.
    for name_, form in fn['forms'].items():
        if form['cls'] not in fn['form_classes']:
            err(f'form {name_}: unknown class {form["cls"]}')
    for name_, s in fn['settings'].items():
        if s['realm'] not in fn['realms']:
            err(f'setting {name_}: unknown realm {s["realm"]}')

    # The backbone: unique names, one root, parents exist, ranks descend.
    nodes = {}
    for b in product['backbone']:
        if b['name'] in nodes:
            err(f'backbone name {b["name"]} repeated')
        nodes[b['name']] = b
        if b['rank'] is not None and b['rank'] not in RANKS:
            err(f'backbone {b["name"]}: unknown rank {b["rank"]}')
        for k in b['sources']:
            if k not in keys:
                err(f'backbone {b["name"]}: source key {k} not in taxonomy_sources.json')
    roots = [b for b in nodes.values() if b['parent'] is None]
    if [r['name'] for r in roots] != ['Life']:
        err(f'backbone roots {[r["name"] for r in roots]}, expected only Life')
    for b in nodes.values():
        try:
            chain = backbone_ancestors(nodes, b['name'])
        except ValueError as e:
            err(str(e))
            continue
        if b['rank'] is not None:
            level = RANKS.index(b['rank'])
            for anc in chain[1:]:
                r = nodes[anc]['rank']
                if r is not None and RANKS.index(r) >= level:
                    err(f'backbone {b["name"]} ({b["rank"]}) sits under {anc} ({r})')

    # Entries: kinds, identifiers, names, vocabularies, sources.
    entries = [*product['modules'], *product['groups'], *product['communities']]
    ids, names_seen = set(), {}
    for e in entries:
        kind = e.get('kind')
        req = REQUIRED['group' if kind in ('group', 'subgroup') else kind] if kind in ('module', 'group', 'subgroup',
                                                                                       'community') else None
        if req is None:
            err(f'{e.get("id")}: unknown kind {kind}')
            continue
        for k in req:
            if k not in e:
                err(f'{e["id"]}: missing field {k}')
        if not re.fullmatch(ID_PATTERN[kind], e['id']):
            err(f'{e["id"]}: identifier does not match {ID_PATTERN[kind]}')
        if e['id'] in ids:
            err(f'identifier {e["id"]} repeated')
        ids.add(e['id'])
        low = e['working_name'].strip().lower()
        if not low:
            err(f'{e["id"]}: empty working name')
        if low in names_seen:
            err(f'working name "{e["working_name"]}" used by {names_seen[low]} and {e["id"]}')
        names_seen[low] = e['id']
        if not e['summary'].strip():
            err(f'{e["id"]}: empty summary')
        if e['status'] not in STATUSES:
            err(f'{e["id"]}: unknown status {e["status"]}')
        if e['scope'] not in SCOPES:
            err(f'{e["id"]}: unknown scope {e["scope"]}')
        if e['basis'] not in BASES:
            err(f'{e["id"]}: unknown basis {e["basis"]}')
        if e['status'] == 'test_first' and not e['tests']:
            err(f'{e["id"]}: test_first without a named test')
        if not e['sources']:
            err(f'{e["id"]}: cites no project file')
        for s in e['sources']:
            p = Path(s['path'])
            if p.is_absolute() or '..' in p.parts:
                err(f'{e["id"]}: source path {s["path"]} must be relative to the repository root')
            elif not (ROOT / p).is_file():
                err(f'{e["id"]}: cited file {s["path"]} does not exist')
            elif s['path'] not in product['producer']['inputs'] and \
                    s['path'] not in product['producer'].get('cited_files', []):
                err(f'{e["id"]}: cited file {s["path"]} is in neither producer.inputs nor producer.cited_files')
        for r in e.get('references', []):
            if not r.get('name') or not r.get('role'):
                err(f'{e["id"]}: reference without name or role')
            for k in r.get('lit', []):
                if k not in keys:
                    err(f'{e["id"]}: literature key {k} not in taxonomy_sources.json')
        for k in e.get('night', []):
            if k not in fn['night']:
                err(f'{e["id"]}: unknown night strategy {k}')
        for k in e.get('clock', []):
            if k not in fn['clocks']:
                err(f'{e["id"]}: unknown clock {k}')

    modules = {m['id']: m for m in product['modules']}
    groups = {g['id']: g for g in product['groups']}
    candidates = {c['id'] for c in product['aerophyte_names']['candidates']}

    # Groups: placement by descent and by function, and parentage.
    for g in groups.values():
        gid = g['id']
        if g['anchor'] not in nodes:
            err(f'{gid}: anchor {g["anchor"]} not in the backbone')
            continue
        if g['degree'] not in DEGREES:
            err(f'{gid}: unknown degree {g["degree"]}')
        else:
            floor = DEGREES[g['degree']]['anchor_at_or_above']
            rank = nodes[g['anchor']]['rank']
            if rank is not None and RANKS.index(rank) > RANKS.index(floor):
                err(f'{gid}: degree {g["degree"]} needs an anchor at {floor} or above, got {g["anchor"]} ({rank})')
        if g['form'] is None:
            if not g['form_candidates']:
                err(f'{gid}: no form and no form candidates')
        elif g['form'] not in fn['forms']:
            err(f'{gid}: unknown form {g["form"]}')
        for f in g['form_candidates']:
            if f not in fn['forms']:
                err(f'{gid}: unknown form candidate {f}')
        for s in g['size']:
            if s not in SIZE_NAMES:
                err(f'{gid}: unknown size class {s}')
        if not g['settings']:
            err(f'{gid}: no settings')
        for s in g['settings']:
            if s not in fn['settings']:
                err(f'{gid}: unknown setting {s}')
        if g['realm'] not in fn['realms']:
            err(f'{gid}: unknown realm {g["realm"]}')
        elif g['realm'] not in {fn['settings'][s]['realm'] for s in g['settings'] if s in fn['settings']}:
            err(f'{gid}: realm {g["realm"]} is not the realm of any of its settings')
        if not g['guilds']:
            err(f'{gid}: no guild')
        for k in g['guilds']:
            if k not in fn['guilds']:
                err(f'{gid}: unknown guild {k}')
        for m in g['modules']:
            if m not in modules:
                err(f'{gid}: unknown module {m}')
            elif modules[m]['status'] == 'set_aside':
                err(f'{gid}: carries the set-aside module {m}')
        for k in g['carries']:
            if k not in ELEMENTS:
                err(f'{gid}: unknown element {k}')
        if g['carries'] and 'carrier' not in g['guilds']:
            err(f'{gid}: carries elements without the carrier guild')
        for c in g['name_candidates']:
            if c not in candidates:
                err(f'{gid}: unknown name candidate {c}')
        if g['kind'] == 'group' and g['parent'] is not None:
            err(f'{gid}: a group has no parent')
        if g['kind'] == 'subgroup':
            parent = groups.get(g['parent'])
            if parent is None or parent['kind'] != 'group':
                err(f'{gid}: parent {g["parent"]} is not a group')
                continue
            if parent['anchor'] in nodes and parent['anchor'] not in backbone_ancestors(nodes, g['anchor']):
                err(f'{gid}: anchor {g["anchor"]} lies outside its parent\'s anchor {parent["anchor"]}')
            if STATUS_ORDER.index(g['status']) > STATUS_ORDER.index(parent['status']):
                err(f'{gid}: status {g["status"]} above its parent\'s {parent["status"]}')
            if parent['scope'] == 'in_places' and g['scope'] != 'in_places':
                err(f'{gid}: a subgroup of an in-places group is in places')
            if not g['parent'] < gid:
                err(f'{gid}: numbered before its parent {g["parent"]}')
    used = {m for g in groups.values() for m in g['modules']}
    for m in modules.values():
        if m['status'] != 'set_aside' and m['id'] not in used:
            err(f'{m["id"]}: no group carries this module')

    # Communities.
    for c in product['communities']:
        if not c['members']:
            err(f'{c["id"]}: no members')
        for m in c['members']:
            if m not in groups:
                err(f'{c["id"]}: unknown member {m}')
        for s in c['settings']:
            if s not in fn['settings']:
                err(f'{c["id"]}: unknown setting {s}')
        if c['realm'] not in fn['realms']:
            err(f'{c["id"]}: unknown realm {c["realm"]}')

    # The aerophytes' names: the author chose the common names; the scientific name waits for the author.
    an = product['aerophyte_names']
    if an['scientific_name'] is not None:
        err('a scientific name is marked chosen for the aerophytes; the author chooses')
    if groups.get(an['group'], {}).get('working_name') != an['common_name']:
        err(f'the aerophytes group {an["group"]} does not carry the common name {an["common_name"]!r}')
    if set(groups.get(an['group'], {}).get('name_candidates', [])) != candidates:
        err('the aerophytes group does not list every candidate')
    taken = {}
    for c in an['candidates']:
        if not re.fullmatch(r'N\d', c['id']):
            err(f'{c["id"]}: identifier does not match N followed by a digit')
        for k in ('group', 'giants', 'roots', 'meaning', 'sound', 'image', 'logic', 'place', 'checks', 'sources'):
            if not c.get(k):
                err(f'{c["id"]}: missing {k}')
        forms = {f.lower() for part in ('group', 'giants', 'young') if c.get(part)
                 for f in (c[part].get('english'), c[part].get('plural'), c[part].get('scientific')) if f}
        for f in forms:
            if f in taken:
                err(f'{c["id"]}: name {f} also offered by {taken[f]}')
            taken[f] = c['id']
        for r in c['roots']:
            if r['source'] not in keys:
                err(f'{c["id"]}: root source {r["source"]} not in taxonomy_sources.json')
        for k in c['sources']:
            if k not in keys:
                err(f'{c["id"]}: source key {k} not in taxonomy_sources.json')
        sci = c['group'].get('scientific')
        if sci and not re.fullmatch(r'[A-Z][a-z]+', sci):
            err(f'{c["id"]}: scientific name {sci} is not one capitalised word of the Latin alphabet')

    # System references to literature.
    for item in (*system['layers'], *system['principles']):
        for k in item['sources']:
            if k not in keys:
                err(f'system: source key {k} not in taxonomy_sources.json')
    for rule in system['naming']['rules']:
        for k in rule['sources']:
            if k.split(':')[0] not in keys:
                err(f'naming rule: source key {k} not in taxonomy_sources.json')
    return errors


def main(argv=None) -> int:
    product = build()
    problems = check(product)
    if problems:
        print('\n'.join(problems), file=sys.stderr)
        return 1
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
    counts = {k: len(product[k]) for k in ('backbone', 'modules', 'groups', 'communities')}
    print(json.dumps(dict(counts=counts, cycle=product['system']['cycle'],
                          candidates=[c['group']['plural'] for c in NAME_CANDIDATES]), indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
