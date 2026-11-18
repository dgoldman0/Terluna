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
"""
from __future__ import annotations

import hashlib

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
        order='head + modifier (Ziwa Nyota)',
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
        order='modifier + head, fused (Yanaqucha, black lake)',
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

RENDAKU = (('sh', 'j'), ('ch', 'j'), ('ts', 'z'), ('k', 'g'), ('t', 'd'), ('s', 'z'), ('h', 'b'))


def _cap(word):
    """Capital on the first letter, past a leading glottal stop ('uo'uo -> 'Uo'uo)."""
    i = 1 if word[:1] == "'" else 0
    return word[:i] + word[i:i + 1].upper() + word[i + 1:]


def compose(tongue, head, mod):
    """A name from a head and a modifier key in a tongue."""
    t = TONGUES[tongue]
    h = t['heads'][head]
    m = t['mods'][mod]
    join = t['join']
    if join == 'fused':
        if m[-1:] == h[:1] and m[-1:] not in 'aeiou':
            h = h[1:]                                   # identical consonants coalesce at the boundary
        if m[-1:] == h[:1] and m[-1:] in 'aeiou' and m[-2:-1] == m[-1:]:
            h = h[1:]                                   # equal vowels cap at two letters
        return _cap(m + h)
    if join == 'rendaku':
        if head not in ('lake', 'lakes', 'basin', 'desert'):
            for a, b in RENDAKU:
                if h.startswith(a):
                    h = b + h[len(a):]
                    break
        return _cap(m + h)
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
      humid woodland and savanna Pwani; hill country Satogo; lowland rainforest Tagalog-Hebrew."""
    z = region['zones']
    lat = abs(region['at']['turned'][0])
    desert = z.get('hyper-arid', 0) + z.get('arid', 0)
    semi = z.get('semi-arid', 0)
    relief = region['relief_m']
    if lat > 62:
        return 'olivolingvo'
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


def name_all(record, body_names, near_side_of):
    """Tongues for the regions and names for every region, district and coastal sea in the draft record, in place."""
    used = set()
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
        r['name'], _, _ = unique_name(r['tongue'], heads, mods, used)
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
        d['name'], d['name_head'], d['name_mod'] = unique_name(tongue, heads, district_mods(d, near_side_of(d)), used)
    land_by_code = {d['code']: d for d in record['districts']}
    for c in record['coastal_seas']:
        facing = [land_by_code[x] for x in c['faces'] if x in land_by_code]
        tongue = facing[0]['tongue'] if facing else 'terlunan'
        c['tongue'] = tongue
        head = coastal_head(c)
        c['name'], c['name_head'], c['name_mod'] = unique_name(tongue, [head, 'coast', 'bay', 'sea'], coastal_mods(c),
                                                               used)
    return record
