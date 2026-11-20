"""The Moon's shared sounds: the framework's checks and repairs, each parent's rule sheet, and the lexicons, the
names table and the draft record standing as passed."""
import json

import pytest

from research.studies.atlas_round import district_names as N
from research.studies.atlas_round import lunar_sounds as L


def test_the_framework_tells_a_word_that_fits_from_one_that_does_not():
    assert L.syllables('kostengis') == (['kos', 'te', 'ngis'], [])             # ng is one consonant and opens a syllable
    assert L.syllables('kaainga') == (['kaa', 'i', 'nga'], [])                  # a doubled vowel is one nucleus
    assert L.syllables("ma'ruf") == (["ma'", 'ruf'], [])                        # the glottal stop is a consonant and may close one
    assert L.syllables("tau'am")[1] == [] and L.syllables('haud')[0] == ['haud']           # au is one nucleus
    assert L.faults('kratero') == ['two consonants open a syllable']
    assert L.faults('bahr') == ['two consonants close a syllable']
    assert L.faults('zahav') == ['outside the inventory: z']
    assert L.faults('shol') == ['two consonants open a syllable']               # sh is two sounds here, not one
    assert L.said('karabogas') == 'ka.ra.BO.gas' and L.said('aral') == 'A.ral'  # the syllable before the last
    assert L.said('haukaainga') == 'hau.KAA.i.nga'                              # a long vowel takes the stress
    assert L.said("H'Naharo Keraha") == 'ha na.HA.ro ke.RA.ha'                  # the article said as Hebrew says its own
    assert L.faults("H'Naharo Keraha") == []                                    # and kept apart from its noun
    assert L.said('Dagat na Kaseh') == 'DA.gat na KA.seh'                       # the linker as written


def test_the_repairs_open_what_the_framework_cannot_hold():
    assert L.repair('kratero') == 'karatero'                                    # a copy of the next vowel
    assert L.repair('kratero', 'e') == 'keratero'                               # or the parent's own vowel
    assert L.repair('pwani') == 'puwani' and L.repair('kyo') == 'kiyo'          # a glide takes its own vowel
    assert L.repair('bult') == 'bulut' and L.repair('els') == 'eles'            # a copy of the vowel before
    assert L.repair('akkala') == 'akala'                                        # a doubled consonant is one
    assert L.repair('impya') == 'impiya' and L.repair('vintra') == 'vintara'    # three consonants running
    assert L.repair('kubwa') == 'kubwa' and L.repair('adrar') == 'adrar'        # a coda meeting an onset stands
    assert L.repair('wadi') == 'wadi'                                           # no vowel is added at the end


def test_parts_fuse_as_the_reference_fuses_compounds():
    assert L.fuse(['bult', 'tengis']) == 'bultengis'                            # the same consonant meeting itself
    assert L.fuse(['bult', 'kol']) == 'bulutkol'                                # a closing cluster opened before a third
    assert L.fuse(['hanak', 'kosa']) == 'hanakosa'
    assert L.fuse(['ak', 'kala']) == 'akala'
    assert L.fuse(['tang', 'ngar']) == 'tangar'                                 # ng meets ng as one consonant
    assert L.fuse(['mata', 'awa']) == 'mataawa'                                 # equal vowels stay a doubled vowel
    assert L.voiced('kawa') == 'gawa' and L.voiced('tani') == 'dani' and L.voiced('sima') == 'sima'


def test_each_parent_keeps_what_the_framework_holds():
    # Kazakh and Kyrgyz: ng, final consonants and the rounding harmony stay; the sibilants fall together; j is y
    assert L.kazakh('жұлдыз') == 'yuldus' and L.kazakh('құдық') == 'kuduk'
    assert L.kazakh('таң') == 'tang' and L.kazakh('жаңа') == 'yanga' and L.kazakh('теңіз') == 'tengis'
    assert L.kazakh('қызыл') == 'kisil' and L.kazakh('шөл') == 'sol'
    assert L.kazakh('жайлау') == 'yailau' and L.kazakh('ой') == 'oi' and L.kazakh('тау') == 'tau'
    assert L.kazakh('ұйқы') == 'uyku' and L.kazakh('уақыт') == 'wakit' and L.kazakh('суық') == 'suuk'
    assert L.kazakh('үсті') == 'ustu' and L.kazakh('асты') == 'asti'            # the harmony follows the vowel before
    # Mongolian: long vowels stay doubled, the old final -i returns
    assert L.mongolian('говь') == 'gobi' and L.mongolian('нуур') == 'nuur' and L.passed('mn:элс') == 'eles'
    # Quechua: q is k and leaves the lowered vowels it gave
    assert L.quechua('qucha') == 'kosa' and L.quechua('urqu') == 'orko' and L.quechua('sunqu') == 'sonko'
    assert L.quechua('llaqta') == 'lakta' and L.quechua('ñawi') == 'nawi' and L.quechua('rimaq') == 'rimak'
    assert L.quechua('wayra') == 'waira' and L.quechua('mayu') == 'mayu' and L.quechua("p'isqu") == 'pesko'
    # Arabic and Tamazight: the throat stops stay as the glottal stop; j is g; closing clusters open
    assert L.arabic('jazira') == 'gasira' and L.arabic('ʿamiq') == "'amik" and L.arabic('dafiʾ') == "dafi'"
    assert L.passed('ar:bahr') == 'bahar' and L.passed('ar:shams') == 'samas' and L.passed('ar:hayya') == 'haya'
    assert L.arabic('layl') == 'lail' and L.arabic('tawil') == 'tawil' and L.passed('ar:tamazirt') == 'tamasirit'
    # Swahili: the syllabic nasal gets its vowel back, by its class
    assert L.swahili('m°to') == 'muto' and L.swahili('mvua') == 'imvua' and L.swahili('nchi') == 'insi'
    assert L.swahili('nyota') == 'inyota' and L.passed('sw:pwani') == 'puwani' and L.passed('sw:mpya') == 'impiya'
    assert L.swahili('jangwa') == 'yangwa' and L.passed('sw:kubwa') == 'kubwa' and L.passed('sw:mwezi') == 'muwesi'
    # Tagalog is the framework's own; Hebrew keeps f and v and opens a cluster with its own e
    assert L.tagalog("lo'ok") == "lo'ok" and L.hebrew('shemesh') == 'semes' and L.hebrew('kokhav') == 'kohav'
    assert L.hebrew('katzeh') == 'kaseh' and L.passed('he:stav') == 'setav' and L.passed('he:zman') == 'seman'
    # Japanese in the spelling of its own sound system
    assert L.japanese('tsuchi') == 'tuti' and L.japanese('shima') == 'sima' and L.japanese('fuka') == 'huka'
    assert L.japanese('mizu') == 'misu' and L.japanese('ō') == 'oo' and L.japanese('kyo') == 'kiyo'
    # the author's Esperanto-based tongue: its letters and sibilants, and its clusters opened as Hebrew opens them
    assert L.esperanto('ĥadaŝa') == 'hadasa' and L.esperanto('jamo') == 'yamo' and L.esperanto('lajla') == 'laila'
    assert L.passed('eo:kratero') == 'keratero' and L.passed('eo:profunda') == 'perofunda' and L.esperanto('ĝemela') == 'gemela'
    assert L.terlunan('whakapono') == 'fakapono' and L.terlunan("feti'a") == "feti'a"
    assert L.passed('la:kaguya') == 'kaguya' and L.passed('la:zāhiya') == 'sahiya'


@pytest.mark.parametrize('tongue', sorted(N.TONGUES))
def test_every_word_of_every_lexicon_stands_as_passed(tongue):
    lex = L.lexicon(tongue)
    assert len(lex) >= len(N.HEADS) + len(N.MODIFIERS)
    for (role, key), e in lex.items():
        assert e['faults'] == [], (tongue, key, e['passed'], e['faults'])      # in the framework, and its parent form passed
        assert L.repair(e['passed'], L.repair_vowel(e['parent'])) == e['passed'], (tongue, key)
    assert L.collisions(tongue) == []                                          # no two different parent words fall together
    for head in N.HEADS:
        for mod in N.MODIFIERS:
            assert L.faults(N.compose(tongue, head, mod)) == [], (tongue, head, mod)   # nor does a joined name break it


def test_the_tables_hold_the_words_as_the_moons_tongues_say_them():
    assert N.TONGUES['dala']['heads']['sea'] == 'tengis' and N.TONGUES['dala']['parents'][('head', 'sea')] == 'kk:теңіз'
    assert N.TONGUES['dala']['mods']['cloud'] == 'bulut' and N.TONGUES['dala']['fuse'][('modifier', 'cloud')] == 'bult'
    assert N.TONGUES['runa']['heads']['lake'] == 'kosa' and N.TONGUES['runa']['parents'][('head', 'lake')] == 'qu:qucha'
    assert N.TONGUES['sahran']['heads']['sea'] == 'Bahar' and N.TONGUES['sahran']['heads']['desert'] == "'Irik"
    assert N.TONGUES['pwani']['heads']['river'] == 'Muto' and N.TONGUES['pwani']['mods']['rain'] == 'Imvua'
    assert N.TONGUES['tagalhebrew']['heads']['bay'] == "Lo'ok" and N.TONGUES['tagalhebrew']['mods']['sun'] == 'Semes'
    assert N.TONGUES['satoyama']['heads']['island'] == 'sima' and N.TONGUES['satoyama']['mods']['great'] == 'oo'
    assert N.TONGUES['olivolingvo']['heads']['sea'] == 'yamo' and N.TONGUES['olivolingvo']['heads']['bowl'] == 'keratero'
    assert N.CALQUE_WORDS['terlunan']['faith'] == 'fakapono' and N.CALQUE_WORDS['tagalhebrew']['edge'] == 'Kaseh'
    moved = [(t, key) for t in N.TONGUES for (role, key), e in L.lexicon(t).items()
             if e['passed'] != e['parent'].split(':', 1)[1].lower()]
    assert len([m for m in moved if m[0] == 'terlunan']) == 1                  # Terlunan's words are the framework's own


def test_the_pins_are_the_names_their_words_make():
    table = N.load_table()
    got = {p['anchor']: p for p in L.pins(table)}
    assert set(L.PINS) | set(L.PINS_BY_PARTS) | set(L.PINS_STAND) == set(table)   # every pin has its keys
    for anchor, p in got.items():
        assert p['stands'], (anchor, p['name'], p['from_the_words'])            # the table and the lexicons agree
        assert p['faults'] == [], (anchor, p['name'], p['faults'])
    names = {p['name'] for p in got.values()}
    assert {'Uamoana', 'Bultengis', "H'Yamo Melaha", 'Yuldustengis', 'Hanakosa', 'Bahari Masariki', 'Dagat na Kaseh'} <= names
    assert {'Hukagawa', 'Muto Kiyani', 'Wairamayu', 'Koilurmayu', 'Yarikgol', 'Rereawa', 'Anuanuafanga', 'Jura', 'Ingenii'} <= names


def test_the_draft_record_stands_as_passed():
    if not N.RECORD.exists():
        pytest.skip('run python -m research.studies.atlas_round.districts first')
    record = json.loads(N.RECORD.read_text())
    units = record['regions'] + record['districts'] + record['coastal_seas']
    for u in units:
        if u['name'] not in ('Jura', 'Ingenii'):                               # the author's own stand as written
            assert L.faults(u['name']) == [], (u['code'], u['name'])
    found, made = L.record_names(record)
    assert found > 600 and made == found                                       # each generated name is what its words make
    seas = {r['anchor']: r['name'] for r in record['regions'] if r['kind'] == 'sea'}
    assert seas['sea:southern near-side seas'] == 'Bultengis' and seas['sea:Mare Moscoviense'] == 'Hanakosa'


def test_the_record_is_built():
    record = L.build(log=lambda *a: None)
    assert record['schema'] == L.SCHEMA and set(record['tongues']) == set(N.TONGUES)
    assert sum(len(t['words']) for t in record['tongues'].values()) == 789
    assert all(t['faults'] == 0 and t['fallen_together'] == [] and t['joined_names_at_fault'] == 0 for t in record['tongues'].values())
    assert len(record['pins']) == 31 and all(p['stands'] for p in record['pins'])
    assert set(record['rule_sheets']) == set(L.READERS)
