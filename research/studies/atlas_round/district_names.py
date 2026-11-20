"""Working names for the atlas's regions and districts in the Moon's own tongues (17 November 2026, first draft).

The author asked for names with no Latin, English or mainland European roots, and for regions mapped to their own
languages, loosely as landscapes on the Moon resemble those of peoples on Earth: Levantine, North African and Arabian
drylands, lush forest lands, Pacific-island seas. The tongues below are lunar languages, each descended in the setting
from Earth's families as Terlunan is, the author's main conlang (terlunan-language-v0.6.md: Polynesian-rooted with
Hakka, Sumerian and Akkadian). Two come from the author's own work: Terlunan, and Proto-Olivolingvo (an
Esperanto-based tongue with Hebrew loans, the article h' and the plural -im, from the author's poem). The Tagalog-Hebrew
tongue follows the author's idea of such a hybrid. The rest are built here from everyday landscape words of their
source languages, without sacred or protected terms. Each tongue's words are listed with their meaning; Terlunan's
words outside its v0.6 dictionary are marked COINED, proposals for the author.

A name is a head chosen by a place's centre (crater, lake, river, plain, coast...) and a modifier chosen by what sets it
apart (its moisture zone, rock, height, lakes, falls, the near side's Earthlight...), joined as each tongue joins them.
The choice is deterministic: the same draft gives the same names. Names are working labels, for the people who come to
live there to keep or replace.

Every word is written below as its parent language spells it, and passes when this module loads through the Moon's
shared sounds (lunar_sounds.py; register, 20 November 2026): a name keeps its parent language's sounds to a degree and
comes out in the Moon's own phonetic framework. After loading, the tables hold the words as the Moon's tongues say
them (tengis, kosa, bahar), each tongue's `parents` holds the parent forms and its `fuse` the forms names are joined
from, and every name is made from the passed words.
"""
from __future__ import annotations

import hashlib

from research.studies.atlas_round import lunar_sounds as LS

HEADS = ('bowl', 'basin', 'lake', 'lakes', 'river', 'valley', 'plain', 'hills', 'mountains', 'ridge', 'cliff', 'island',
         'bay', 'coast', 'strait', 'sea', 'forest', 'desert', 'field', 'land', 'falls', 'spring')
MODIFIERS = ('black', 'white', 'red', 'green', 'gold', 'blue', 'great', 'small', 'high', 'deep', 'long', 'old', 'new',
             'dry', 'wet', 'warm', 'cold', 'light', 'dawn', 'night', 'star', 'cloud', 'mist', 'rain', 'wind', 'sun',
             'earth', 'stone', 'sand', 'salt', 'water', 'leaf', 'flower', 'thousand', 'twin')

TONGUES = {
    'terlunan': dict(
        label='Terlunan', roots='the author\'s conlang v0.6 (Tahitian, Māori, Hawaiian, Hakka, Sumerian, Akkadian)',
        order='modifier-head, fused, as Terlunan compounds (dungu + tap = dungutap)',
        heads=dict(bowl='poka', basin='galpoka', lake='loko', lakes='siyenloko', river='awa', valley='awaawa',
                   plain='mania', hills='puke', mountains='kur', ridge='hiwi', cliff='pari', island='motu', bay='fanga',
                   coast='tahatai', strait='koawa', sea='moana', forest='lim', desert='puehu', field='matie',
                   land='fenua', falls='rere', spring='puna'),
        mods=dict(black='pooiri', white="'uo'uo", red="'ute'ute", green='liyukset', gold='kim', blue='ninamu', great='gal',
                  small='iti', high='teitei', deep='hohonu', long='roa', old='ul', new='gibil', dry='maroke',
                  wet='maaku', warm='mehana', cold='maatao', light='nuru', dawn='ata', night='musu', star="feti'a",
                  cloud='dungu', mist='kohu', rain='ua', wind="mata'i", sun='raa', earth='honua', stone='abnu',
                  sand='one', salt='miti', water='tama', leaf='rau', flower='tiare', thousand='siyen', twin='liyong'),
        coined=dict(loko='lake (Hawaiian loko; roto is taken by inside)', siyenloko='thousand lakes (siyen + loko)',
                    awa='river (Māori awa)', awaawa='valley (Māori awaawa)', mania='plain (Māori mania)',
                    hiwi='ridge (Māori hiwi)', tahatai='shore, coast (taha side + tai tide)',
                    koawa='channel, strait (Māori kōawa)', matie='grass, grassland (Tahitian matie)',
                    rere='waterfall (Māori rere)', puna='spring (Māori puna)', teitei='high (Tahitian teitei)',
                    hohonu='deep (Tahitian hohonu)', maroke='dry (Māori maroke)', maaku='wet (Māori mākū)',
                    kohu='mist (Māori kohu)', **{"'ute'ute": "red (Tahitian 'ute'ute, beside 'ere'ere black and 'uo'uo white)"},
                    galpoka='great crater, basin (gal + poka)'),
        join='fused'),
    'olivolingvo': dict(
        label='Proto-Olivolingvo', roots='the author\'s Esperanto-based tongue with Hebrew loans (h\' the, -im plural)',
        order='h\' + noun + adjective, as in the author\'s poem (h\'naharo keraha, the frozen river)',
        heads=dict(bowl='kratero', basin='emeko', lake='lago', lakes='lagim', river='naharo', valley='valo',
                   plain='ebeno', hills='montetim', mountains='harim', ridge='kresto', cliff='klifo', island='insulo',
                   bay='golfo', coast='bordo', strait='markolo', sea='jamo', forest='jaaro', desert='midbaro',
                   field='kampo', land='erco', falls='akvofalo', spring='eino'),
        mods=dict(black='shaxora', white='lavana', red='aduma', green='jaroka', gold='zahava', blue='kaĥola',
                  great='gadola', small='ketana', high='alta', deep='profunda', long='longa', old='atika',
                  new='ĥadaŝa', dry='seka', wet='malseka', warm='varma', cold='keraha', light='luma', dawn='ŝaĥara',
                  night='lajla', star='kokava', cloud='anana', mist='nebula', rain='geŝema', wind='ruaĥa',
                  sun='ŝemeŝa', earth='tera', stone='ŝtona', sand='sabla', salt='melaĥa', water='majima',
                  leaf='folia', flower='flora', thousand='mila', twin='ĝemela'),
        join='olivo'),
    'sahran': dict(
        label='Sahran (Arabic- and Amazigh-rooted)', roots='Arabic and Tamazight landscape words, simplified sounds',
        order='head + modifier, as Arabic names run (Wadi Nur); Tamazight heads where they fit (Adrar, Asif, Tala)',
        heads=dict(bowl='Hufra', basin='Hawd', lake='Tamda', lakes='Tamdawin', river='Asif',
                   valley='Wadi', plain='Hamada', hills='Tilal', mountains='Adrar', ridge='Ighil', cliff='Jurf',
                   island='Jazira', bay='Halij', coast='Sahil', strait='Madik', sea='Bahr', forest='Gaba',
                   desert='Erg', field='Marj', land='Tamazirt', falls='Shallal', spring='Tala'),
        mods=dict(black='Aswad', white='Abyad', red='Ahmar', green='Ahdar', gold='Dahab', blue='Azrak', great='Kabir',
                  small='Sagir', high='Afella', deep='Amik', long='Tawil', old='Kadim', new='Jadid', dry='Yabis',
                  wet='Ratib', warm='Dafi', cold='Barid', light='Nur', dawn='Fajr', night='Layl', star='Najm',
                  cloud='Sahab', mist='Dabab', rain='Matar', wind='Rih', sun='Shams', earth='Ard', stone='Hajar',
                  sand='Raml', salt='Milh', water='Aman', leaf='Waraka', flower='Zahr', thousand='Alf', twin='Tawam'),
        join='spaced'),
    'tagalhebrew': dict(
        label='Tagalog-Hebrew', roots='Tagalog heads and Hebrew modifiers, the author\'s hybrid idea',
        order='Tagalog head, the linker -ng after a vowel or na after a consonant, Hebrew modifier (Lawang Or)',
        heads=dict(bowl='Uka', basin='Lambak', lake='Lawa', lakes='Lawalawa', river='Ilog', valley='Lambak',
                   plain='Kapatagan', hills='Burol', mountains='Bundok', ridge='Gulod', cliff='Bangin', island='Pulo',
                   bay='Look', coast='Baybay', strait='Kipot', sea='Dagat', forest='Gubat', desert='Ilang',
                   field='Parang', land='Lupa', falls='Talon', spring='Bukal'),
        mods=dict(black='Shahor', white='Lavan', red='Adom', green='Yarok', gold='Zahav', blue='Kahol', great='Gadol',
                  small='Katan', high='Gavoah', deep='Amok', long='Aroh', old='Atik', new='Hadash', dry='Yavesh',
                  wet='Ratuv', warm='Ham', cold='Kar', light='Or', dawn='Shahar', night='Laila', star='Kohav',
                  cloud='Anan', mist='Arafel', rain='Geshem', wind='Ruah', sun='Shemesh', earth='Erets',
                  stone='Even', sand='Hol', salt='Melah', water='Mayim', leaf='Ale', flower='Perah', thousand='Elef',
                  twin='Teom'),
        join='tagalog'),
    'pwani': dict(
        label='Pwani (Swahili-rooted)', roots='Swahili landscape words, without the noun-class agreement',
        order='head + modifier (Siwa Inyota)',
        heads=dict(bowl='Shimo', basin='Bonde', lake='Ziwa', lakes='Maziwa', river='Mto', valley='Bonde',
                   plain='Tambarare', hills='Vilima', mountains='Milima', ridge='Kingo', cliff='Genge', island='Kisiwa',
                   bay='Guba', coast='Pwani', strait='Mlango', sea='Bahari', forest='Msitu', desert='Jangwa',
                   field='Mbuga', land='Nchi', falls='Maporomoko', spring='Chemchemi'),
        mods=dict(black='Eusi', white='Eupe', red='Ekundu', green='Kijani', gold='Dhahabu', blue='Samawati',
                  great='Kubwa', small='Dogo', high='Juu', deep='Kina', long='Refu', old='Kale', new='Mpya', dry='Kavu',
                  wet='Rutuba', warm='Joto', cold='Baridi', light='Nuru', dawn='Alfajiri', night='Usiku', star='Nyota',
                  cloud='Wingu', mist='Ukungu', rain='Mvua', wind='Upepo', sun='Jua', earth='Dunia', stone='Jiwe',
                  sand='Mchanga', salt='Chumvi', water='Maji', leaf='Jani', flower='Ua', thousand='Elfu', twin='Pacha'),
        join='spaced'),
    'dala': dict(
        label='Dala (Turkic- and Mongolic-rooted)', roots='Kazakh, Kyrgyz and Mongolian landscape words',
        order='modifier + head, fused (Karakol, black lake)',
        heads=dict(bowl='oy', basin='ova', lake='kol', lakes='kollar', river='gol', valley='say', plain='dala',
                   hills='tobe', mountains='tau', ridge='kyr', cliff='jar', island='aral', bay='koltyk', coast='jaga',
                   strait='bogaz', sea='tengiz', forest='orman', desert='gobi', field='jaylau', land='jer',
                   falls='sarkyrama', spring='bulak'),
        mods=dict(black='kara', white='ak', red='kyzyl', green='jasyl', gold='altyn', blue='kok', great='ulug',
                  small='kishi', high='biik', deep='teren', long='uzun', old='eski', new='jana', dry='kurgak',
                  wet='ylgal', warm='jyly', cold='suyk', light='jaryk', dawn='tan', night='tun', star='juldyz',
                  cloud='bult', mist='tuman', rain='jangbyr', wind='jel', sun='kun', earth='jer', stone='tas',
                  sand='kum', salt='tuz', water='su', leaf='japyrak', flower='gul', thousand='myng', twin='egiz'),
        join='fused'),
    'runa': dict(
        label='Runa (Quechua-rooted)', roots='Quechua landscape words, spelled without ejective marks',
        order='modifier + head, fused (Yanakosa, black lake)',
        heads=dict(bowl='hutku', basin='kancha', lake='qucha', lakes='quchakuna', river='mayu', valley='wayqu',
                   plain='pampa', hills='pata', mountains='urqu', ridge='qasa', cliff='qaqa', island='wata',
                   bay='qucha', coast='chimpa', strait='chaka', sea='mamaqucha', forest='sacha', desert='purun',
                   field='chakra', land='llaqta', falls='paqcha', spring='pukyu'),
        mods=dict(black='yana', white='yuraq', red='puka', green='qumir', gold='quri', blue='anqas', great='hatun',
                  small='huchuy', high='hanaq', deep='ukhu', long='suni', old='nawpa', new='musuq', dry='chaki',
                  wet='huqu', warm='quni', cold='chiri', light='kanchay', dawn='paqarin', night='tuta',
                  star='quyllur', cloud='puyu', mist='pukutay', rain='para', wind='wayra', sun='inti', earth='allpa',
                  stone='rumi', sand='tiyu', salt='kachi', water='yaku', leaf='rapi', flower='tika',
                  thousand='waranqa', twin='iskay'),
        join='fused'),
    'satoyama': dict(
        label='Satogo (Japanese-rooted)', roots='Japanese landscape words, with sequential voicing in compounds',
        order='modifier + head, fused, the head voiced (Kurogawa, black river)',
        heads=dict(bowl='kubo', basin='bonchi', lake='ko', lakes='ikeno', river='kawa', valley='tani', plain='hara',
                   hills='oka', mountains='yama', ridge='one', cliff='gake', island='shima', bay='ura', coast='hama',
                   strait='seto', sea='umi', forest='mori', desert='suna', field='no', land='kuni', falls='taki',
                   spring='izumi'),
        mods=dict(black='kuro', white='shira', red='aka', green='midori', gold='kogane', blue='ao', great='o',
                  small='ko', high='taka', deep='fuka', long='naga', old='furu', new='nii', dry='kara', wet='nure',
                  warm='nuku', cold='samu', light='hikari', dawn='asa', night='yo', star='hoshi', cloud='kumo',
                  mist='kiri', rain='ame', wind='kaze', sun='hi', earth='tsuchi', stone='ishi', sand='suna',
                  salt='shio', water='mizu', leaf='ha', flower='hana', thousand='chi', twin='futa'),
        join='rendaku'),
}

WHERE = {'terlunan': 'shores and islands, and the seas', 'olivolingvo': 'the polar lands and the semi-arid uplands',
         'sahran': 'the deserts', 'tagalhebrew': 'the lowland rainforests', 'pwani': 'the humid woodlands and savannas',
         'dala': 'the steppe plains', 'runa': 'the high mountains and their cloud forests',
         'satoyama': 'the humid hill country'}

RENDAKU = LS.RENDAKU                 # the sequential voicing the Moon's sounds can show: k, t and h to g, d and b


def _cap(word):
    """Capital on the first letter, past a leading glottal stop ('uo'uo -> 'Uo'uo)."""
    i = 1 if word[:1] == "'" else 0
    return word[:i] + word[i:i + 1].upper() + word[i + 1:]


def compose(tongue, head, mod):
    """A name from a head and a modifier key in a tongue; the modifier may be the tongue's word for an old name's
    meaning. A fused name is joined from the words' fusing forms and repaired as one word, as the Moon's sounds have
    it (lunar_sounds.py): identical consonants meeting at the join become one, and a cluster the join makes is opened."""
    t = TONGUES[tongue]
    role = 'modifier' if mod in t['mods'] else 'meaning'
    h = t['heads'][head]
    m = t['mods'][mod] if role == 'modifier' else CALQUE_WORDS[tongue][mod]
    join = t['join']
    if join in ('fused', 'rendaku'):
        hf, mf = t['fuse'][('head', head)], t['fuse'][(role, mod)]
        if join == 'rendaku' and head not in ('lake', 'lakes', 'basin', 'desert'):
            hf = LS.voiced(hf)
        return _cap(LS.repair(LS.fuse([mf, hf])))
    if join == 'olivo':
        return f"H'{_cap(h)} {_cap(m)}"
    if join == 'tagalog':
        return f"{h}ng {m}" if h[-1] in 'aeiou' else f"{h} na {m}"
    return f'{h} {m}'


def _hash(code, k):
    return int(hashlib.sha256(code.encode()).hexdigest()[:8], 16) % max(k, 1)


def tongue_for(region):
    """The tongue of a land region, from its character, loosely as Earth's peoples sit in like landscapes:
    - the polar lands: Proto-Olivolingvo, whose poem is of winter;
    - the dry lands, by what dominates: deserts Sahran; semi-arid broken uplands of a Levantine kind Proto-Olivolingvo;
      semi-arid plains, the steppe, Dala;
    - the lush lands: high mountains with their cloud forests Runa; shores and islands Terlunan, the seafaring tongue;
      humid woodland and savanna Pwani; hill country Satogo; lowland rainforest Tagalog-Hebrew;
    - an island group out in international waters, outside the polar lands: Terlunan, the islands' tongue."""
    z = region['zones']
    lat = abs(region['at']['turned'][0])
    desert = z.get('hyper-arid', 0) + z.get('arid', 0)
    semi = z.get('semi-arid', 0)
    relief = region['relief_m']
    if lat > 62:
        return 'olivolingvo'
    if region.get('island_group'):                          # an island group in international waters: the islands' tongue
        return 'terlunan'
    if desert + semi >= 0.6:
        if desert >= semi:
            return 'sahran'
        return 'olivolingvo' if relief >= 3000 else 'dala'
    if relief >= 3500:
        return 'runa'
    if region['shore_share'] >= 0.6:
        return 'terlunan'
    if z.get('humid', 0) + z.get('dry sub-humid', 0) > z.get('wet', 0):
        return 'pwani'
    return 'satoyama' if relief >= 2900 else 'tagalhebrew'


def district_head(d):
    """The head of a land district's name from its centre and character."""
    if d['lakes'] >= 8 and d['lake_share'] >= 0.05:
        return 'lakes'
    if d['lake_share'] >= 0.2 or (d['centre'] == 'lake' and d['lake_share'] >= 0.05):
        return 'lake'
    if d['falls_gw'] and d['falls_gw'][0] >= 10:
        return 'falls'
    if d['centre'] == 'island':
        return 'island'
    if d['centre'] == 'crater':
        return 'bowl' if (d['centre_km'] or 0) < 150 else 'basin'
    if d['largest_flow_m3s'] >= 5000:
        return 'valley' if d['relief_m'] > 2000 else 'river'
    if d['relief_m'] > 3500:
        return 'mountains'
    if d['relief_m'] > 2000:
        return 'hills' if _hash(d['code'], 2) else 'ridge'
    zone = d['zone']
    if zone in ('wet', 'humid', 'dry sub-humid'):
        return 'forest'
    if zone == 'semi-arid':
        return 'field'
    if zone in ('arid', 'hyper-arid'):
        return 'desert'
    return 'plain'


ZONE_MODS = {'wet': ('rain', 'cloud', 'mist'), 'humid': ('green', 'leaf', 'rain'),
             'dry sub-humid': ('gold', 'flower', 'leaf'), 'semi-arid': ('wind', 'gold', 'sun'),
             'arid': ('red', 'sand', 'stone'), 'hyper-arid': ('white', 'salt', 'sun')}
FALLBACK = ('dawn', 'night', 'light', 'long', 'old', 'twin', 'blue', 'new', 'warm', 'water', 'star', 'great', 'small',
            'high', 'deep', 'cold', 'dry', 'wet', 'earth', 'stone', 'gold', 'green', 'red', 'white', 'black', 'thousand')


def district_mods(d, near_side):
    out = []
    if d['heritage']:
        out.append('star')
    if d['cut']:
        out.append('new')
    if d['falls_gw'] and d['falls_gw'][0] >= 10:
        out.append('mist')
    if d['lakes'] >= 8:
        out.append('thousand')
    if d['basalt_share'] >= 0.5:
        out.append('black')
    zm = list(ZONE_MODS.get(d['zone'], ()))
    k = _hash(d['code'], len(zm))
    out += zm[k:] + zm[:k]
    if d['mean_height_m'] > 3500:
        out.append('high')
    if d['mean_height_m'] < -300:
        out.append('deep')
    if near_side:
        out.append('earth')
    out.append('great' if d['area_km2'] > 100_000 else 'small' if d['area_km2'] < 10_000 else 'old')
    return out + [m for m in FALLBACK if m not in out]


def coastal_head(c):
    return 'strait' if c.get('strait') else 'bay' if c.get('enclosed') else 'coast'


def coastal_mods(c):
    out = ['salt'] if c['salt'] else []
    out += ['warm', 'green'] if c['mean_depth_m'] < 200 else ['blue', 'deep'] if c['mean_depth_m'] > 1500 else ['light']
    if c['fed_by']:
        out.append('water')
    return out + [m for m in FALLBACK if m not in out]


def unique_name(tongue, heads, mods, used):
    for head in heads:
        for mod in mods:
            name = compose(tongue, head, mod)
            if name not in used:
                used.add(name)
                return name, head, mod
    for mod2 in ('great', 'small', 'old', 'new', 'twin', 'long'):
        for mod in mods:
            name = compose(tongue, heads[0], mod) + ' ' + compose(tongue, 'land', mod2).split()[-1]
            if name not in used:
                used.add(name)
                return name, heads[0], mod
    raise ValueError('no name left')


SEA_NAMES = {
    'near-side sea': ('Honuamoana', 'Earth Sea: the sea under Earth\'s gaze'),
    'southern near-side seas': ('Itimoana', 'little seas: the near-side sea\'s southern basins'),
    'South Pole–Aitken Sea': ('Mitimoana', 'salt sea'),
    'Mare Humboldtianum': ('Fetiamoana', 'star sea, in the north towards the limb'),
    'Mare Moscoviense': ('Teiteimoana', 'high sea, 3.4 km above the common level'),
    'Mare Orientale': ('Pi\'omoana', 'arch sea, in the rings of its basin'),
    'Smythii–Marginis Sea': ('Tahamoana', 'side sea, at the eastern limb'),
}


def name_all(record, body_names, near_side_of, reserved=()):
    """Tongues for the regions and names for every region, district and coastal sea in the draft record, in place,
    none of them taking a reserved name."""
    import collections
    used = set(reserved)
    seen = collections.defaultdict(collections.Counter)      # heads and modifiers already taken, by scope

    def fresh(keys, scope):
        """The keys in their order of fit, those not yet taken in the scope first, so that neighbours differ."""
        return sorted(keys, key=lambda k: seen[scope][k])

    def take(scope, head, mod):
        seen[(scope, 'head')][head] += 1
        seen[(scope, 'mod')][mod] += 1

    regions = {r['code']: r for r in record['regions']}
    for r in record['regions']:
        if r['kind'] == 'sea':
            continue
        r['tongue'] = tongue_for(r)
    land_regions = [r for r in record['regions'] if r['kind'] == 'land']
    for r in land_regions:
        z = max(r['zones'], key=r['zones'].get)
        k = _hash(r['code'], 3)
        if r['relief_m'] >= 3500:
            fits = ['mountains', 'land', 'valley', 'ridge']
        elif r['relief_m'] >= 2900:
            fits = ['hills', 'land', 'valley', 'ridge']
        elif z in ('arid', 'hyper-arid'):
            fits = ['desert', 'land', 'plain', 'hills']
        elif z == 'semi-arid':
            fits = ['field', 'land', 'plain', 'river']
        else:
            fits = ['forest', 'land', 'plain', 'river']
        if r.get('lake_share', 0) >= 0.12:
            fits.insert(1, 'lakes')
        heads = fits[k:] + fits[:k]
        zm = list(ZONE_MODS.get(z, ()))
        mods = (['cold', 'night'] if abs(r['at']['turned'][0]) > 62 else []) + \
            (['black'] if r['basalt_share'] > 0.3 else []) + zm[k:] + zm[:k] + \
            (['earth'] if near_side_of(r) else []) + [m for m in FALLBACK]
        scope = ('regions', r['tongue'])                    # a tongue's regions differ in head and modifier
        r['name'], head, mod = unique_name(r['tongue'], fresh(heads, (scope, 'head')), fresh(mods, (scope, 'mod')), used)
        take(scope, head, mod)
    for r in record['regions']:
        if r['kind'] == 'sea':
            r['tongue'] = 'terlunan'
            r['name'], r['name_meaning'] = SEA_NAMES.get(body_names[r['code']], (r['code'], ''))
            used.add(r['name'])
    for d in record['districts']:
        tongue = regions[d['region']]['tongue']
        d['tongue'] = tongue
        head = district_head(d)
        heads = [head] + [h for h in ('land', 'plain', 'hills', 'field') if h != head]
        scope = ('districts', d['region'])                  # a region's districts keep their landform, vary the rest
        d['name'], d['name_head'], d['name_mod'] = unique_name(
            tongue, heads, fresh(district_mods(d, near_side_of(d)), (scope, 'mod')), used)
        take(scope, d['name_head'], d['name_mod'])
    land_by_code = {d['code']: d for d in record['districts']}
    for c in record['coastal_seas']:
        facing = [land_by_code[x] for x in c['faces'] if x in land_by_code]
        tongue = facing[0]['tongue'] if facing else 'terlunan'
        c['tongue'] = tongue
        head = coastal_head(c)
        scope = ('coasts', c['region'])
        c['name'], c['name_head'], c['name_mod'] = unique_name(tongue, [head, 'coast', 'bay', 'sea'],
                                                               fresh(coastal_mods(c), (scope, 'mod')), used)
        take(scope, c['name_head'], c['name_mod'])
    return record


# ---------------------------------------------------------------- names kept by anchor, and the old names' meanings
#
# Names live apart from the draw. Each unit has an anchor that survives a redraw: the sea's body, the named feature it
# carries (districts.feature_units), the crater at a district's centre, or a centre's place; a region takes the anchor
# of its largest district. place_names.json pins a name to an anchor, so a hand-set name stays through any redraw whose
# anchor still exists, and `python -m research.studies.atlas_round.district_names` renames the record in seconds
# without drawing it again. A unit that carries a named mare, ocean, bay, lake or marsh keeps that old name's meaning
# (a calque): its own head in its tongue, with the meaning as the modifier. Everything else takes the generated name.

from pathlib import Path as _Path

PLACE_NAMES = _Path(__file__).resolve().parent / 'place_names.json'
RECORD = _Path(__file__).resolve().parent / 'results' / 'districts.json'

# The IAU names' meanings: (gloss, modifier key). Names after people, probes and Earth's places carry none.
CALQUES = {
    'Oceanus Procellarum': ('Ocean of Storms', 'storm'), 'Mare Imbrium': ('Sea of Showers', 'rain'),
    'Mare Frigoris': ('Sea of Cold', 'cold'), 'Mare Tranquillitatis': ('Sea of Tranquillity', 'calm'),
    'Mare Serenitatis': ('Sea of Serenity', 'clear'), 'Mare Fecunditatis': ('Sea of Fertility', 'plenty'),
    'Mare Nubium': ('Sea of Clouds', 'cloud'), 'Mare Humorum': ('Sea of Moisture', 'wet'),
    'Mare Nectaris': ('Sea of Nectar', 'nectar'), 'Mare Crisium': ('Sea of Crises', 'trial'),
    'Mare Vaporum': ('Sea of Vapours', 'mist'), 'Mare Insularum': ('Sea of Islands', 'islands'),
    'Mare Cognitum': ('Sea That Has Become Known', 'known'), 'Mare Australe': ('Southern Sea', 'south'),
    'Mare Orientale': ('Eastern Sea', 'east'), 'Mare Marginis': ('Sea of the Edge', 'edge'),
    'Mare Undarum': ('Sea of Waves', 'wave'), 'Mare Spumans': ('Foaming Sea', 'foam'),
    'Mare Anguis': ('Serpent Sea', 'serpent'), 'Sinus Roris': ('Bay of Dew', 'dew'),
    'Sinus Medii': ('Central Bay', 'middle'), 'Sinus Aestuum': ('Bay of Billows', 'wave'),
    'Sinus Amoris': ('Bay of Love', 'love'), 'Sinus Asperitatis': ('Bay of Roughness', 'rough'),
    'Sinus Concordiae': ('Bay of Harmony', 'harmony'), 'Sinus Fidei': ('Bay of Faith', 'faith'),
    'Sinus Successus': ('Bay of Success', 'success'), 'Lacus Somniorum': ('Lake of Dreams', 'dream'),
    'Lacus Mortis': ('Lake of Death', 'death'), 'Lacus Temporis': ('Lake of Time', 'time'),
    'Lacus Doloris': ('Lake of Sorrow', 'sorrow'), 'Lacus Veris': ('Lake of Spring', 'springtime'),
    'Lacus Autumni': ('Lake of Autumn', 'autumn'), 'Lacus Aestatis': ('Lake of Summer', 'summer'),
    'Lacus Hiemalis': ('Wintry Lake', 'winter'), 'Lacus Spei': ('Lake of Hope', 'hope'),
    'Lacus Lenitatis': ('Lake of Tenderness', 'gentle'), 'Lacus Bonitatis': ('Lake of Goodness', 'good'),
    'Lacus Perseverantiae': ('Lake of Perseverance', 'endure'), 'Lacus Solitudinis': ('Lake of Solitude', 'alone'),
    'Lacus Oblivionis': ('Lake of Forgetfulness', 'forget'), 'Lacus Timoris': ('Lake of Fear', 'fear'),
    'Lacus Felicitatis': ('Lake of Happiness', 'joy'), 'Lacus Gaudii': ('Lake of Joy', 'joy'),
    'Lacus Luxuriae': ('Lake of Luxury', 'plenty'), 'Lacus Excellentiae': ('Lake of Excellence', 'excellent'),
    'Lacus Tenebrarum': ('Lake of Darkness', 'dark'), 'Palus Somni': ('Marsh of Sleep', 'sleep'),
    'Palus Putredinis': ('Marsh of Decay', 'decay'), 'Palus Epidemiarum': ('Marsh of Epidemics', 'plague'),
}

# The meanings' words in each tongue, beside the tongues' own modifiers (which already hold rain, cold, cloud, wet and
# mist). Terlunan takes the v0.6 dictionary's words where it has them (here, mung, maarama, sung, salalu, mong,
# monamona, damku, taime, oto, maru, wareware, pooiri, taha, motu); the rest are COINED in its manner from its donors.
CALQUE_WORDS = {
    'terlunan': dict(storm='vero', calm='marino', clear='maarama', plenty='momona', nectar='monamona',
                     trial='tamata', islands='motu', known="'ite", south='tonga', east='hikina', edge='taha',
                     wave="'are", foam='huka', serpent='nakahi', dew='tomairangi', middle='sung', love='here',
                     rough='taratara', harmony='aumaru', faith='whakapono', success='angitu', dream='mung',
                     death='mate', time='taime', sorrow='oto', springtime='koanga', autumn='ngahuru',
                     summer='raumati', winter='hotoke', hope='mong', gentle='maru', good='damku',
                     endure='manawanui', alone='mokemoke', forget='wareware', fear='mataku', joy='koa',
                     excellent='hiranga', dark='pooiri', sleep='salalu', decay='pirau', plague='mate'),
    'olivolingvo': dict(storm='ŝtorma', calm='trankvila', clear='serena', plenty='fekunda', nectar='nektara',
                        trial='kriza', islands='insula', known='konata', south='suda', east='orienta', edge='randa',
                        wave='onda', foam='ŝauma', serpent='serpenta', dew='rosa', middle='meza', love='ama',
                        rough='malglata', harmony='konkorda', faith='fida', success='sukcesa', dream='sonĝa',
                        death='morta', time='tempa', sorrow='dolora', springtime='printempa', autumn='aŭtuna',
                        summer='somera', winter='vintra', hope='espera', gentle='milda', good='bona',
                        endure='persista', alone='sola', forget='forgesa', fear='timema', joy='ĝoja',
                        excellent='bonega', dark='malluma', sleep='dorma', decay='putra', plague='pesta'),
    'sahran': dict(storm='Asifa', calm='Hadi', clear='Safi', plenty='Khisb', nectar='Shahd', trial='Mihna',
                   islands='Juzur', known='Maruf', south='Janub', east='Sharq', edge='Hafa', wave='Mawj',
                   foam='Zabad', serpent='Hayya', dew='Nada', middle='Wasat', love='Hubb', rough='Khashin',
                   harmony='Wiam', faith='Iman', success='Najah', dream='Ahlam', death='Mawt', time='Zaman',
                   sorrow='Huzn', springtime='Rabi', autumn='Kharif', summer='Sayf', winter='Shita', hope='Amal',
                   gentle='Latif', good='Khayr', endure='Sabr', alone='Uzla', forget='Nisyan', fear='Khawf',
                   joy='Farah', excellent='Fadl', dark='Zalam', sleep='Nawm', decay='Afan', plague='Waba'),
    'tagalhebrew': dict(storm='Sufa', calm='Shaket', clear='Tsalul', plenty='Shefa', nectar='Tsuf', trial='Mashber',
                        islands='Iyim', known='Yadua', south='Darom', east='Kedem', edge='Katzeh', wave='Gal',
                        foam='Ketsef', serpent='Nachash', dew='Tal', middle='Tikhon', love='Ahava', rough='Gas',
                        harmony='Harmonia', faith='Emuna', success='Hatslaha', dream='Halom', death='Mavet',
                        time='Zman', sorrow='Tsaar', springtime='Aviv', autumn='Stav', summer='Kayits',
                        winter='Horef', hope='Tikva', gentle='Adin', good='Tov', endure='Hatmada', alone='Bdidut',
                        forget='Shikhecha', fear='Pahad', joy='Simha', excellent='Metsuyan', dark='Hoshekh',
                        sleep='Shena', decay='Rikavon', plague='Dever'),
    'pwani': dict(storm='Dhoruba', calm='Shwari', clear='Angavu', plenty='Wingi', nectar='Asali', trial='Shida',
                  islands='Visiwa', known='Julikana', south='Kusini', east='Mashariki', edge='Ukingo', wave='Mawimbi',
                  foam='Povu', serpent='Nyoka', dew='Umande', middle='Kati', love='Upendo', rough='Kparuparu',
                  harmony='Maelewano', faith='Imani', success='Fanaka', dream='Ndoto', death='Mauti', time='Wakati',
                  sorrow='Huzuni', springtime='Masika', autumn='Vuli', summer='Kiangazi', winter='Kipupwe',
                  hope='Tumaini', gentle='Pole', good='Wema', endure='Subira', alone='Upweke', forget='Sahau',
                  fear='Hofu', joy='Furaha', excellent='Bora', dark='Giza', sleep='Usingizi', decay='Uozo',
                  plague='Tauni'),
    'dala': dict(storm='boran', calm='tynysh', clear='ashyk', plenty='mol', nectar='bal', trial='synak',
                 islands='araldar', known='belgili', south='ontustik', east='shygys', edge='shet', wave='tolkyn',
                 foam='kobik', serpent='jylan', dew='shyk', middle='orta', love='mahabbat', rough='kedir',
                 harmony='yntymak', faith='senim', success='jenis', dream='tus', death='olim', time='uakyt',
                 sorrow='kaygy', springtime='koktem', autumn='kuz', summer='jaz', winter='kys', hope='umit',
                 gentle='jumsak', good='jaksy', endure='tozim', alone='jalgyz', forget='umyt', fear='korky',
                 joy='kuanysh', excellent='asyl', dark='karangy', sleep='uiky', decay='shirik', plague='oba'),
    'runa': dict(storm='wayrapara', calm='thak', clear='chuya', plenty='achka', nectar='mishki', trial='sasa',
                 islands='watakuna', known='riqsisqa', south='urin', east='inti', edge='kuchu', foam='pusuqu',
                 serpent='machaqway', dew='sulla', middle='chawpi', love='munay', harmony='huklla', faith='iñiy',
                 success='atipay', dream='musquy', death='wanuy', time='mita', sorrow='llaki', springtime='pawqar',
                 autumn='aymuray', summer='ruphay', winter='chiraw', hope='suyay', gentle='llampu', good='alli',
                 endure='kallpa', alone='sapan', forget='qunqay', fear='manchay', joy='kusi', excellent='sumaq',
                 dark='tutayaq', sleep='puñuy', decay='ismu', plague='unquy'),
    'satoyama': dict(storm='arashi', calm='nagi', clear='hare', plenty='yutaka', nectar='mitsu', trial='kiki',
                     islands='shima', known='shiru', south='minami', east='higashi', edge='fuchi', wave='nami',
                     foam='awa', serpent='hebi', dew='tsuyu', middle='naka', love='ai', rough='ara', harmony='wa',
                     faith='shin', success='sachi', dream='yume', time='toki', sorrow='kanashi', springtime='haru',
                     autumn='aki', summer='natsu', winter='fuyu', hope='nozomi', gentle='yasa', good='yoi',
                     endure='shinobu', alone='hitori', forget='wasure', fear='osore', joy='yoro', excellent='masa',
                     dark='yami', sleep='nemuri', decay='kuchi', plague='eki'),
}


# ---------------------------------------------------------------- the words passed through the Moon's shared sounds
#
# The tables above write each word as its parent language spells it in Latin letters. PARENTS gives a fuller parent
# form where that spelling had lost a sound the Moon's framework can keep (Kazakh ng and its vowels, the Arabic throat
# stops, the class of a Swahili syllabic nasal, Tagalog's glottal stop, Quechua ñ) or where the first draft's word was
# not the parent's: a code for the parent language's rule sheet (lunar_sounds.READERS), a colon and the word in the
# parent's own letters. Dala's words are all given so, since the Latin spelling could not tell its vowels apart.

PARENTS = {
    'dala': dict(
        bowl='kk:ой', basin='la:owa', lake='kk:көл', lakes='kk:көлдөр', river='mn:гол', valley='kk:сай', plain='kk:дала', hills='kk:төбе',
        mountains='kk:тау', ridge='kk:қыр', cliff='kk:жар', island='kk:арал', bay='kk:қолтық', coast='kk:жаға', strait='kk:боғаз',
        sea='kk:теңіз', forest='kk:орман', desert='mn:говь', field='kk:жайлау', land='kk:жер', falls='kk:сарқырама', spring='kk:бұлақ',
        black='kk:қара', white='kk:ақ', red='kk:қызыл', green='kk:жасыл', gold='kk:алтын', blue='kk:көк', great='la:ulug', small='kk:кіші',
        high='kk:биік', deep='kk:терең', long='kk:ұзын', old='kk:ескі', new='kk:жаңа', dry='kk:құрғақ', wet='kk:ылғал', warm='kk:жылы',
        cold='kk:суық', light='kk:жарық', dawn='kk:таң', night='kk:түн', star='kk:жұлдыз', cloud='kk:бұлт', mist='kk:тұман',
        rain='kk:жаңбыр', wind='kk:жел', sun='kk:күн', earth='kk:жер', stone='kk:тас', sand='kk:құм', salt='kk:тұз', water='kk:су',
        leaf='kk:жапырақ', flower='kk:гүл', thousand='kk:мың', twin='kk:егіз',
        storm='kk:боран', calm='kk:тыныш', clear='kk:ашық', plenty='kk:мол', nectar='kk:бал', trial='kk:сынақ', islands='kk:аралдар',
        known='kk:белгілі', south='kk:оңтүстік', east='kk:шығыс', edge='kk:шет', wave='kk:толқын', foam='kk:көбік', serpent='kk:жылан',
        dew='kk:шық', middle='kk:орта', love='kk:махаббат', rough='kk:кедір', harmony='kk:ынтымақ', faith='kk:сенім', success='kk:жеңіс',
        dream='kk:арман', death='kk:өлім', time='kk:уақыт', sorrow='kk:қайғы', springtime='kk:көктем', autumn='kk:күз', summer='kk:жаз',
        winter='kk:қыс', hope='kk:үміт', gentle='kk:жұмсақ', good='kk:жақсы', endure='kk:төзім', alone='kk:жалғыз', forget='kk:унут',
        fear='kk:қорқы', joy='kk:қуаныш', excellent='kk:асыл', dark='kk:қараңғы', sleep='kk:ұйқы', decay='kk:шірік', plague='kk:оба'),
    'runa': dict(old='qu:ñawpa', death='qu:wañuy', faith='qu:iñiy', sleep='qu:puñuy'),
    'sahran': dict(desert='ar:ʿirq', deep='ar:ʿamiq', warm='ar:dafiʾ', twin='ar:tawʾam', storm='ar:ʿasifa', calm='ar:hadiʾ',
                   known='ar:maʿruf', harmony='ar:wiʾam', springtime='ar:rabiʿ', winter='ar:shitaʾ', alone='ar:ʿuzla', decay='ar:ʿafan',
                   plague='ar:wabaʾ', dark='ar:ʿatma'),
    'pwani': dict(river='sw:m°to', strait='sw:m°lango', forest='sw:m°situ', sand='sw:m°changa', rough='sw:paruparu'),
    'tagalhebrew': dict(bay="tl:lo'ok", sorrow="he:tsa'ar"),
    'satoyama': dict(great='ja:ō'),
    'olivolingvo': {},
    'terlunan': {},
}
# What a parent form puts right beyond the sounds: the draft's word was another language's, a slip, or would fall
# together with another word of the same tongue once passed.
PARENT_NOTES = {
    ('dala', 'lakes'): "the draft's kollar is the Uzbek plural; көлдөр is the Kyrgyz",
    ('dala', 'dream'): "the draft's tus (түс) would fall together with salt (тұз); арман is Kazakh for a dream or longing",
    ('dala', 'forget'): "the draft's umyt (ұмыт) would fall together with hope (үміт); унут is the Kyrgyz",
    ('sahran', 'dark'): "the draft's Zalam would read as the Arabic for peace once z is s; ʿatma is darkness",
    ('pwani', 'rough'): "the draft's Kparuparu was a slip for paruparu",
    ('tagalhebrew', 'sorrow'): "the throat stop between the two vowels of tsaʿar keeps them apart",
}
# The rule sheet a tongue's own table passes through where PARENTS gives no fuller form; the Tagalog-Hebrew tongue
# reads its heads as Tagalog and its modifiers as Hebrew.
_READER = dict(terlunan='te', olivolingvo='eo', sahran='ar', pwani='sw', dala='kk', runa='qu', satoyama='ja')


def parent_spec(tongue, key, word, role):
    """The parent form a word passes from: PARENTS' where it gives one, else the table's own spelling read by the
    tongue's rule sheet."""
    spec = PARENTS.get(tongue, {}).get(key)
    if spec:
        return spec
    if tongue == 'tagalhebrew':
        return ('tl:' if role == 'head' else 'he:') + word
    return f'{_READER[tongue]}:{word}'


def _pass_through():
    """Every word passed from its parent form into the Moon's shared sounds, in place: the tables then hold the words
    as the Moon's tongues say them (capitalised where the table's were), each tongue's `parents` the parent forms and
    its `fuse` the forms names are joined from."""
    for tongue, t in TONGUES.items():
        t['parents'], t['fuse'] = {}, {}
        for role, table in (('head', t['heads']), ('modifier', t['mods']), ('meaning', CALQUE_WORDS.setdefault(tongue, {}))):
            for key, word in list(table.items()):
                spec = parent_spec(tongue, key, word, role)
                alone = LS.passed(spec)
                t['parents'][(role, key)] = spec
                t['fuse'][(role, key)] = LS.fuse(LS.read(spec))
                table[key] = _cap(alone) if word[:1].isupper() else alone


_pass_through()


def load_table(path=PLACE_NAMES):
    """The names table: {anchor: {name, tongue, meaning, note}}; empty where there is none yet."""
    import json
    if not _Path(path).exists():
        return {}
    return json.loads(_Path(path).read_text()).get('names', {})


def pinned_features(path=PLACE_NAMES):
    """The gazetteer features the names table pins by name beyond the maria, oceans, bays, lakes and marshes."""
    return sorted(a.split(':', 1)[1] for a in load_table(path) if a.startswith('feature:'))


def anchor_of(u, units=None):
    """A unit's anchor: what it is tied to that a redraw keeps."""
    if u.get('kind') == 'sea':
        return f"sea:{u['body']}"
    if u.get('kind') == 'land' and 'districts' in u and not isinstance(u['districts'], list):
        mine = [d for d in (units or {}).get('districts', []) if d['region'] == u['code']]
        if not mine:
            return f"region:{u['at']['today'][0]:.1f},{u['at']['today'][1]:.1f}"
        return 'region:' + max(mine, key=lambda d: d['area_km2'])['anchor']
    if u.get('feature'):
        return f"feature:{u['feature']['name']}"
    if u.get('centre_name'):
        return f"crater:{u['centre_name']}"
    if u.get('kind') == 'land':
        la, lo = u.get('centre_at') or u['at']['today']
        return f"{u.get('centre', 'ground')}:{la:.1f},{lo:.1f}"
    la, lo = u['at']['today']
    return f"coast:{la:.1f},{lo:.1f}"


def calque(u, used):
    """The name that keeps the meaning of the old name the unit carries, in its tongue on its own head, if the tongue
    has the word; with the old name and its gloss."""
    f = u.get('feature')
    if not f or f['name'] not in CALQUES:
        return None
    gloss, key = CALQUES[f['name']]
    t = TONGUES[u['tongue']]
    if key not in CALQUE_WORDS.get(u['tongue'], {}) and key not in t['mods']:
        return None
    head = u.get('name_head') or 'land'
    name = compose(u['tongue'], head, key)
    if name in used:
        return None
    return name, f"{key} {head}, after {f['name']}, the {gloss}"


def apply(units, table=None, log=print, previous=None):
    """Names for every region, district and coastal sea, in place: the generated names first; then a unit whose
    anchor held a generated name in the previous record keeps it, so a redraw does not reshuffle names; then each unit
    that carries an old name with a meaning takes its calque; then the names table's pins by anchor. Descriptions then
    read names for codes. Records each name's source and meaning. previous is {anchor: (name, meaning, source)}."""
    import re
    table = load_table() if table is None else table
    previous = previous or {}
    body = {r['code']: r.get('body') for r in units['regions'] if r['kind'] == 'sea'}
    near_side = lambda x: abs(x['at']['turned'][1]) < 90.0                        # noqa: E731
    name_all(units, body, near_side)                       # the tongues and heads the calques take
    pins = {v['name'] for v in table.values()}             # every pinned name, the seas' and regions' included
    calques, missing, taken = {}, [], set(pins)
    for x in units['districts'] + units['coastal_seas']:
        if x.get('anchor') in table or not x.get('feature') or x['feature']['name'] not in CALQUES:
            continue
        got = calque(x, taken)
        if got:
            calques[x['code']] = got
            taken.add(got[0])
        else:
            missing.append((x['code'], x['tongue'], CALQUES[x['feature']['name']][1]))
    every = units['regions'] + units['districts'] + units['coastal_seas']
    kept = {}                                               # the previous record's generated names, by anchor
    for x in every:
        old = previous.get(x.get('anchor'))
        if old and old[2] == 'generated' and old[3] == x.get('tongue') and x['code'] not in calques \
                and x.get('anchor') not in table and old[0] not in taken and old[0] not in kept.values():
            kept[x['code']] = old[0]
    name_all(units, body, near_side, reserved=taken | set(kept.values()))   # new names keep clear of all those
    for x in every:
        x['name_source'] = 'generated'
        x.setdefault('name_meaning', '')
        if x['code'] in kept:
            old = previous[x['anchor']]
            x['name'] = kept[x['code']]
            x['name_meaning'] = old[1] or x['name_meaning']
            if old[4]:
                x['name_head'], x['name_mod'] = old[4], old[5]
        if x['code'] in calques:
            x['name'], x['name_meaning'] = calques[x['code']]
            x['name_source'] = 'calque'
    for x in every:
        pin = table.get(x.get('anchor'))
        if pin:
            x['name'] = pin['name']
            x['name_meaning'] = pin.get('meaning', '')
            x['name_source'] = 'pinned'
            if pin.get('tongue'):
                x['tongue'] = pin['tongue']
    names = {x['code']: x['name'] for x in every}
    for x in every:
        x['description'] = re.sub(r'\b([LS]\d+(?:\.\d+)?)\b', lambda m: names.get(m.group(1), m.group(1)),
                                  x.get('description_codes', x['description']))
    if missing:
        log('calques left out, the word missing in the tongue or the name taken: '
            + ', '.join(f'{c} ({t}: {k})' for c, t, k in missing))
    return units


def disambiguate(units):
    """Anchors made unique, in place: where several units share one (a large crater's floor split into districts, a
    fill-in point's ground split by the size rules), each takes its own label point, to a degree, else to a tenth,
    else its order among them by code."""
    from collections import defaultdict
    groups = defaultdict(list)
    for u in units:
        groups[u['anchor']].append(u)
    for base, us in groups.items():
        if len(us) < 2:
            continue
        us.sort(key=lambda u: u['code'])
        for digits in (0, 1):
            keys = [f"{base}@{u['at']['today'][0]:.{digits}f},{u['at']['today'][1]:.{digits}f}" for u in us]
            if len(set(keys)) == len(keys):
                break
        else:
            keys = [f'{base}#{k}' for k in range(len(us))]
        for u, k in zip(us, keys):
            u['anchor'] = k


def anchor_all(units):
    """Every unit's anchor, unique: the districts' and coastal seas' first, then the regions' from their largest
    districts."""
    for x in units['districts'] + units['coastal_seas']:
        x['anchor'] = anchor_of(x)
    disambiguate(units['districts'] + units['coastal_seas'])
    for x in units['regions']:
        x['anchor'] = anchor_of(x, units)
    disambiguate(units['regions'])


def previous_names(record):
    """{anchor: (name, meaning, source, tongue, head, modifier)} for every unit of a districts record, its anchors
    found again by today's rules, so that an older record's units match today's like with like."""
    import copy
    units = copy.deepcopy(dict(regions=record.get('regions', []), districts=record.get('districts', []),
                               coastal_seas=record.get('coastal_seas', [])))
    if any('anchor' not in x for x in units['districts'] + units['coastal_seas']):
        return {}
    anchor_all(units)
    return {x['anchor']: (x['name'], x.get('name_meaning', ''), x.get('name_source', 'generated'), x.get('tongue'),
                          x.get('name_head'), x.get('name_mod'))
            for x in units['regions'] + units['districts'] + units['coastal_seas']}


def main(argv=None) -> int:
    """Rename the draft in place from the names table and the rules, without drawing it again."""
    import argparse
    import json
    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument('--record', type=_Path, default=RECORD)
    p.add_argument('--previous', type=_Path, help='a record whose generated names to keep by anchor (default: the record '
                                                 'itself, so a rename moves only what the table changes)')
    p.add_argument('--fresh', action='store_true', help='generate every name afresh')
    args = p.parse_args(argv)
    record = json.loads(args.record.read_text())
    prev = None if args.fresh else previous_names(json.loads(args.previous.read_text()) if args.previous else record)
    units = dict(regions=record['regions'], districts=record['districts'], coastal_seas=record['coastal_seas'])
    if any('anchor' not in x for x in units['districts'] + units['coastal_seas']):
        raise SystemExit(f'{args.record} predates the anchors; draw it again with research.studies.atlas_round.districts')
    apply(units, previous=prev)
    record['names'] = dict(table=str(PLACE_NAMES.relative_to(PLACE_NAMES.parents[3])),
                           sources={s: sum(1 for x in units['regions'] + units['districts'] + units['coastal_seas']
                                           if x['name_source'] == s) for s in ('pinned', 'calque', 'generated')})
    args.record.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    print(f"renamed {args.record}: {record['names']['sources']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
