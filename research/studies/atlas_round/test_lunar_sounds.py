"""The Moon's shared sounds: the framework's checks and repairs, each parent's rule sheet, and the lexicons and pins
passed through."""
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
    assert L.said("H'Naharo Keraha") == "h'na.HA.ro ke.RA.ha" and L.faults("H'Naharo Keraha") == []   # the article apart
    assert L.said('Dagat na Kaseh') == 'DA.gat na KA.seh'                       # the linker as written


def test_the_repairs_open_what_the_framework_cannot_hold():
    assert L.repair('kratero') == 'karatero'                                    # a copy of the next vowel
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
    # Tagalog is the framework's own; Hebrew keeps f and v
    assert L.tagalog("lo'ok") == "lo'ok" and L.hebrew('shemesh') == 'semes' and L.hebrew('kokhav') == 'kohav'
    assert L.hebrew('katzeh') == 'kaseh' and L.passed('he:stav') == 'satav'
    # Japanese in the spelling of its own sound system
    assert L.japanese('tsuchi') == 'tuti' and L.japanese('shima') == 'sima' and L.japanese('fuka') == 'huka'
    assert L.japanese('mizu') == 'misu' and L.japanese('ō') == 'oo' and L.japanese('kyo') == 'kiyo'
    # the author's Esperanto-based tongue; Terlunan itself
    assert L.esperanto('ĥadaŝa') == 'hadasa' and L.esperanto('jamo') == 'yamo' and L.esperanto('lajla') == 'laila'
    assert L.passed('eo:kratero') == 'karatero' and L.esperanto('ĝemela') == 'gemela'
    assert L.terlunan('whakapono') == 'fakapono' and L.terlunan("feti'a") == "feti'a"
    assert L.passed('la:kaguya') == 'kaguya' and L.passed('la:zāhiya') == 'sahiya'


@pytest.mark.parametrize('tongue', sorted(N.TONGUES))
def test_every_word_of_every_lexicon_passes(tongue):
    lex = L.lexicon(tongue)
    assert len(lex) >= len(N.HEADS) + len(N.MODIFIERS)
    for (role, key), e in lex.items():
        assert e['faults'] == [], (tongue, key, e['passed'], e['faults'])      # it stands in the framework
        assert L.repair(e['passed']) == e['passed'], (tongue, key)             # and the repairs leave it alone
        assert L.faults(L.repair(e['read'])) == []
    assert L.collisions(tongue) == []                                          # no two words the draft kept apart fall together


def test_terlunan_is_the_frameworks_own():
    lex = L.lexicon('terlunan')
    moved = {key: e['passed'] for (role, key), e in lex.items() if e['passed'] != e['draft'].lower()}
    assert moved == {'faith': 'fakapono'}                                      # one coined word still in its Māori spelling


def test_names_are_joined_as_before_from_the_passed_words():
    assert L.compose('dala', 'lake', 'black') == 'Karakol'                     # nothing to change
    assert L.compose('dala', 'sea', 'cloud') == 'Bultengis'                    # t meets t; z is s
    assert L.compose('dala', 'lake', 'cloud') == 'Bulutkol'                    # the closing cluster opened
    assert L.compose('runa', 'lake', 'high') == 'Hanakosa'
    assert L.compose('satoyama', 'river', 'black') == 'Kurogawa'               # the voicing the framework can show
    assert L.compose('satoyama', 'island', 'black') == 'Kurosima'              # s has no voiced partner in it
    assert L.compose('tagalhebrew', 'lake', 'light') == 'Lawang Or'
    assert L.compose('tagalhebrew', 'bay', 'black') == "Lo'ok na Sahor"
    assert L.compose('sahran', 'sea', 'black') == 'Bahar Aswad'                # each word repaired on its own
    assert L.compose('sahran', 'desert', 'red') == "'Irik Ahmar"
    assert L.compose('pwani', 'river', 'green') == 'Muto Kiyani'
    assert L.compose('olivolingvo', 'river', 'cold') == "H'Naharo Keraha"      # the article as the author writes it
    assert L.compose('terlunan', 'bay', 'salt') == 'Mitifanga'
    for tongue in N.TONGUES:
        for head in N.HEADS:
            for mod in N.MODIFIERS:
                assert L.faults(L.compose(tongue, head, mod)) == [], (tongue, head, mod)


def test_the_pins_are_found_again_from_their_keys_and_passed():
    table = N.load_table()
    for anchor, (tongue, head, mod) in L.PINS.items():
        assert L.draft_compose(tongue, head, mod) == table[anchor]['name'], anchor      # the keys make today's spelling
    got = {p['name']: p for p in L.pins(table)}
    assert all(p['passed'] for p in got.values())                               # every pin has a rule
    assert all(p['faults'] == [] for p in got.values())
    assert got['Bultengiz']['passed'] == 'Bultengis' and got['Juldyztengiz']['passed'] == 'Yuldustengis'
    assert got['Hanaqucha']['passed'] == 'Hanakosa' and got['Bahari Mashariki']['passed'] == 'Bahari Masariki'
    assert got['Dagat na Katzeh']['passed'] == 'Dagat na Kaseh' and got["H'Jamo Melaĥa"]['passed'] == "H'Yamo Melaha"
    assert got['Mto Kijani']['passed'] == 'Muto Kiyani' and got['Fukagawa']['passed'] == 'Hukagawa'
    for name in ('Uamoana', 'Paramayu', 'Rereawa', 'Anuanuafanga', 'Jura', 'Ingenii', 'Asif Dahab', 'Niigawa'):
        assert got[name]['passed'] == name and not got[name]['changed'], name   # the author's own stand; these already fit
    assert sum(p['changed'] for p in got.values()) == 15


def test_the_record_is_built_and_changes_no_name():
    before = (N.TONGUES['dala']['heads']['sea'], N.load_table()['sea:southern near-side seas']['name'])
    record = L.build(log=lambda *a: None)
    assert record['schema'] == L.SCHEMA and set(record['tongues']) == set(N.TONGUES)
    assert sum(len(t['words']) for t in record['tongues'].values()) == 789
    assert sum(t['changed'] for t in record['tongues'].values()) == 324
    assert all(t['faults'] == 0 and t['fallen_together'] == [] for t in record['tongues'].values())
    assert len(record['pins']) == 31 and set(record['rule_sheets']) == set(L.READERS)
    assert (N.TONGUES['dala']['heads']['sea'], N.load_table()['sea:southern near-side seas']['name']) == before
