"""The Moon's shared sounds: every name passes through them and keeps what it can of its parent language's.

The framework is the sound section of the Terlunan reference (v0.6, section 2), taken as the sounds all the Moon's
tongues share: the vowels a e i o u, a long vowel written double; the consonants p b t d k g m n ng s h l r w y and the
glottal stop, with f and v admitted in established words; no sh, z, emphatic or affricate series; syllables (C)V(C)
with ai, au, oi, ou as one nucleus, so a consonant may close a syllable and meet the next one's opening consonant but
no syllable opens or closes on two; no vowel added at the end of a borrowed word; identical consonants meeting where
two parts are fused become one; stress on the syllable before the last unless a long vowel takes it.

Each parent language has a rule sheet below (a reader): a sound the framework has is kept, a sound it lacks goes to its
nearest neighbour inside it, and where there is a choice the reader takes the one that keeps most of the parent. Every
rule is a choice that can be turned; lunar_sounds.md gives each with its reason. The rules are a first set and the
parent words are written from general knowledge of the languages: no dictionary or speaker has checked them.

The lexicons are district_names.py's (TONGUES and CALQUE_WORDS). A word passes from the draft's spelling, or from the
fuller parent form in PARENTS where the draft's spelling had lost a sound the framework can keep or was not a word of
the parent. This module changes no name: district_names.py and place_names.json stand as they are until the passed
forms are adopted.

    python -m research.studies.atlas_round.lunar_sounds      -> results/lunar_sounds.json and a printed summary
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from pathlib import Path

from research.studies.atlas_round import district_names as DN

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'lunar_sounds.json'
SCHEMA = 'terluna.research.atlas-round-lunar-sounds/1'

VOWELS = 'aeiou'
CONSONANTS = ('ng', 'p', 'b', 't', 'd', 'k', 'g', 'm', 'n', 's', 'h', 'l', 'r', 'w', 'y', "'", 'f', 'v')
DIPHTHONGS = ('ai', 'au', 'oi', 'ou')


# ------------------------------------------------------------------------------------------------ the framework

def segments(w):
    """A word's sounds in order: ng is one consonant."""
    out, i = [], 0
    while i < len(w):
        if w[i:i + 2] == 'ng':
            out.append('ng'); i += 2
        else:
            out.append(w[i]); i += 1
    return out


def _v(s):
    return s in VOWELS


def repair(w):
    """The framework's repairs on one word: a doubled consonant becomes one; a cluster that opens the word gets a
    vowel inside it (i before y and u before w, as the reference's liyong and kuwo, else a copy of the next vowel); a
    cluster that closes it is opened with a copy of the vowel before it; a run of three consonants is opened after
    its second, before a glide with the glide's own vowel."""
    s = segments(w)
    s = [x for i, x in enumerate(s) if not (i and x == s[i - 1] and not _v(x))]
    if len(s) > 1 and not _v(s[0]) and not _v(s[1]):
        nxt = next((x for x in s[2:] if _v(x)), 'a')
        s.insert(1, {'y': 'i', 'w': 'u'}.get(s[1], nxt))
    if len(s) > 2 and not _v(s[-1]) and not _v(s[-2]):
        prev = next((x for x in reversed(s[:-2]) if _v(x)), 'a')
        s.insert(len(s) - 1, prev)
    i = 0
    while i + 2 < len(s):
        if not _v(s[i]) and not _v(s[i + 1]) and not _v(s[i + 2]):
            nxt = next((x for x in s[i + 3:] if _v(x)), 'a')
            s.insert(i + 2, {'y': 'i', 'w': 'u'}.get(s[i + 2], nxt))
        i += 1
    return ''.join(s)


def fuse(parts):
    """Parts fused into one word as the reference fuses compounds: identical consonants meeting at the join become
    one; a part that ends on two consonants is opened before a part that begins with a third (bult + kol, bulutkol);
    equal vowels stay as a doubled vowel of at most two letters."""
    out = ''
    for p in parts:
        if not p:
            continue
        if out:
            a, b = segments(out), segments(p)
            if a[-1] == b[0] and not _v(a[-1]):
                b = b[1:]
            if len(a) > 1 and b and not _v(a[-1]) and not _v(a[-2]) and not _v(b[0]):
                prev = next((x for x in reversed(a[:-2]) if _v(x)), 'a')
                a.insert(len(a) - 1, prev)
            out, p = ''.join(a), ''.join(b)
        out += p
    return re.sub(r'([aeiou])\1\1+', r'\1\1', out)


def syllables(w):
    """Syllables (C)V(C), a doubled vowel or ai, au, oi, ou being one nucleus; with the faults against the framework:
    a sound outside it, no vowel, two consonants opening or closing a syllable."""
    s = segments(w)
    faults = [f'outside the inventory: {x}' for x in s if not _v(x) and x not in CONSONANTS]
    nuclei, i = [], 0
    while i < len(s):
        if _v(s[i]):
            if i + 1 < len(s) and _v(s[i + 1]) and (s[i] == s[i + 1] or s[i] + s[i + 1] in DIPHTHONGS):
                nuclei.append((i, i + 2)); i += 2
            else:
                nuclei.append((i, i + 1)); i += 1
        else:
            i += 1
    if not nuclei:
        return [w], faults + ['no vowel']
    syl, start = [], 0
    for k, (a, b) in enumerate(nuclei):
        last = k == len(nuclei) - 1
        between = (len(s) if last else nuclei[k + 1][0]) - b
        take = between if last else max(0, between - 1)
        if a - start > 1:
            faults.append('two consonants open a syllable')
        if take > 1:
            faults.append('two consonants close a syllable')
        syl.append(''.join(s[start:b + take]))
        start = b + take
    return syl, faults


ARTICLE = "h'"                       # the author's Esperanto-based tongue's article, left as the author writes it
PARTICLES = ('na',)                  # the Tagalog linker between a head that ends on a consonant and its modifier


def _bare(word):
    """A word without the article h': no word of the framework opens on h and a glottal stop, so that opening is
    always the article before its noun."""
    w = word.lower()
    return w[len(ARTICLE):] if w.startswith(ARTICLE) and len(w) > len(ARTICLE) else w


def said(w):
    """A word's syllables with the stressed one in capitals: the one before the last, or the rightmost long vowel's.
    The article h' and the linker na are shown as written."""
    out = []
    for word in w.split():
        if word.lower() in PARTICLES:
            out.append(word.lower())
            continue
        bare = _bare(word)
        lead = ARTICLE if bare != word.lower() else ''
        syl, _ = syllables(bare)
        k = max(0, len(syl) - 2)
        longs = [i for i, x in enumerate(syl) if re.search(r'([aeiou])\1', x)]
        if longs:
            k = longs[-1]
        out.append(lead + '.'.join(x.upper() if i == k else x for i, x in enumerate(syl)))
    return ' '.join(out)


def faults(w):
    """A name's faults against the framework, word by word, the article h' apart."""
    return [f for word in w.split() for f in syllables(_bare(word))[1]]


# ------------------------------------------------------------------------------------------------ the rule sheets

_KAZAKH = {'а': 'a', 'ә': 'e', 'б': 'b', 'в': 'w', 'г': 'g', 'ғ': 'g', 'д': 'd', 'е': 'e', 'ж': 'y', 'з': 's', 'и': 'i', 'й': 'Y',
           'к': 'k', 'қ': 'k', 'л': 'l', 'м': 'm', 'н': 'n', 'ң': 'ng', 'о': 'o', 'ө': 'o', 'п': 'p', 'р': 'r', 'с': 's', 'т': 't',
           'у': 'U', 'ұ': 'u', 'ү': 'u', 'ф': 'p', 'х': 'h', 'һ': 'h', 'ц': 's', 'ч': 's', 'ш': 's', 'щ': 's', 'ы': 'I', 'і': 'I',
           'э': 'e', 'ю': 'yu', 'я': 'ya', 'ь': '', 'ъ': ''}
_MONGOLIAN = {'а': 'a', 'б': 'b', 'в': 'b', 'г': 'g', 'д': 'd', 'е': 'ye', 'ё': 'yo', 'ж': 'y', 'з': 's', 'и': 'i', 'й': 'Y',
              'к': 'k', 'л': 'l', 'м': 'm', 'н': 'n', 'о': 'o', 'ө': 'o', 'п': 'p', 'р': 'r', 'с': 's', 'т': 't', 'у': 'u', 'ү': 'u',
              'ф': 'p', 'х': 'h', 'ц': 's', 'ч': 's', 'ш': 's', 'ы': 'i', 'э': 'e', 'ю': 'yu', 'я': 'ya', 'ь': 'i', 'ъ': ''}


def _glides(tokens):
    """Kazakh's short high vowels and glides settled: a short high vowel is u after o or u and i otherwise, the
    rounding harmony of speech (құдық, kuduk); й after a or o closes a diphthong and is y elsewhere; у is a vowel
    after a consonant, closes a diphthong after a or o, and is w before a vowel at the start of a word."""
    out, last = [], ''
    for k, ch in enumerate(tokens):
        nxt = tokens[k + 1] if k + 1 < len(tokens) else ''
        if ch == 'I':
            ch = 'u' if last and last in 'ou' else 'i'
            last = ch
        elif ch == 'Y':
            ch = 'i' if out and out[-1] in ('a', 'o') else 'y'
        elif ch == 'U':
            if out and out[-1] in ('a', 'o'):
                ch = 'u'
            elif out and _v(out[-1]):
                ch = 'w'
            elif not out and nxt and nxt[0] in 'aeiouI':
                ch = 'w'
            else:
                ch = 'u'; last = 'u'
        elif ch and ch[-1] in VOWELS:
            last = ch[-1]
        out.append(ch)
    return ''.join(out)


def kazakh(word):
    """Kazakh and Kyrgyz in their own letters. Kept: ng, final consonants, the rounding harmony. Changed: the uvular
    stop and fricative to k and g; sh and z to s; the j of жылан to y, the sound it answers to across the Turkic
    languages; front rounded ö and ü to o and u."""
    return _glides([_KAZAKH.get(ch, ch) for ch in word.lower()])


def mongolian(word):
    """Mongolian in its own letters. Kept: long vowels, written double. Changed: its sibilants and affricates to s; a
    final soft sign is the old -i (говь, gobi)."""
    w = ''.join(_MONGOLIAN.get(ch, ch) for ch in word.lower())
    w = re.sub(r'([ao])Y', r'\1i', w)
    return w.replace('Y', 'y')


def quechua(word):
    """Quechua in its standard spelling. Kept: the lowering of i and u to e and o beside the uvular q, which Quechua
    speaks and does not write, so that qucha is kosa when q becomes k. Changed: q to k; ch and sh to s; ll and ñ to l
    and n; the aspirated and ejective stops to the plain ones."""
    w = word.lower().replace('ñ', 'N').replace('chh', 'C').replace("ch'", 'C').replace('ch', 'C').replace('sh', 's').replace('ll', 'L')
    w = re.sub(r"([ptkq])[h']", r'\1', w)
    low = {'i': 'e', 'u': 'o'}
    w = re.sub(r'[iu](?=[^aeiou]?q)', lambda m: low[m.group(0)], w)
    w = re.sub(r'(?<=q)[iu]', lambda m: low[m.group(0)], w)
    w = w.replace('q', 'k').replace('C', 's').replace('L', 'l').replace('N', 'n')
    w = re.sub(r'([aeo])y(?![aeiou])', r'\1i', w)
    return re.sub(r'aw(?![aeiou])', 'au', w)


def arabic(word):
    """Arabic and Tamazight in Latin letters, with ʿ and ʾ where the parent has them. Kept: the two throat stops as
    the framework's glottal stop; f. Changed: the emphatics to their plain partners and q to k, as the reference does
    with its Akkadian words; kh and the deep h to h; gh to g; sh, z and th to s, s and t; j to g, as Cairo says it;
    long vowels to short."""
    w = word.lower()
    for a, b in (('ʿ', "'"), ('ʾ', "'"), ('ā', 'a'), ('ī', 'i'), ('ū', 'u'), ('ḥ', 'h'), ('ṣ', 's'), ('ḍ', 'd'), ('ṭ', 't'), ('ẓ', 's'),
                 ('kh', 'h'), ('gh', 'g'), ('sh', 's'), ('th', 't'), ('dh', 'd'), ('j', 'g'), ('z', 's'), ('q', 'k')):
        w = w.replace(a, b)
    w = re.sub(r'([^aeiou])\1', r'\1', w)
    w = re.sub(r'aw(?![aeiou])', 'au', w)
    return re.sub(r'ay(?![aeiou])', 'ai', w)


def swahili(word):
    """Swahili in its standard spelling; m° marks the syllabic m of the tree-and-river class. Kept: its five vowels,
    which are the framework's; f and v; ng as the framework's ng. Changed: sh, ch and z to s; j to y; a syllabic nasal
    gets back the vowel it lost, mu- in the tree-and-river class (mto, muto) and i- before the nasal elsewhere, as
    its sister languages have it (mvua, imvua); a consonant and a glide that open a word get a vowel between."""
    w = word.lower().replace("ng'", 'ng')
    for a, b in (('sh', 's'), ('ch', 's'), ('th', 't'), ('dh', 'd'), ('gh', 'g'), ('kh', 'h'), ('z', 's'), ('j', 'y')):
        w = w.replace(a, b)
    if w.startswith('m°'):
        return 'mu' + w[2:]
    w = w.replace('°', '')
    s = segments(w)
    if len(s) > 1 and s[0] in ('m', 'n', 'ng') and not _v(s[1]) and s[1] != 'w':
        w = 'i' + w
    return w


def tagalog(word):
    """Tagalog in its standard spelling, with its glottal stop written where it stands between vowels. Tagalog's
    sounds are the framework's own."""
    return word.lower()


def hebrew(word):
    """Hebrew in Latin letters. Kept: f and v. Changed: sh, ts and z to s; kh and the deep h to h; q to k."""
    w = word.lower()
    for a, b in (('ʿ', "'"), ('ʾ', "'"), ('tz', 's'), ('ts', 's'), ('sh', 's'), ('kh', 'h'), ('ch', 'h'), ('ḥ', 'h'), ('z', 's'), ('q', 'k')):
        w = w.replace(a, b)
    return w


def japanese(word):
    """Japanese in the Hepburn spelling. Kept: its five vowels and long vowels. Changed to the spelling that follows
    its own sound system: shi, chi, tsu, fu are si, ti, tu, hu; z and j are s; a consonant and y before a vowel get
    an i between (kyo, kiyo)."""
    w = word.lower()
    for a, b in (('ō', 'oo'), ('ū', 'uu'), ('shi', 'si'), ('chi', 'ti'), ('tsu', 'tu'), ('fu', 'hu'), ('ji', 'si'),
                 ('sh', 'siy'), ('ch', 'tiy'), ('j', 'siy'), ('ts', 't'), ('z', 's'), ('f', 'h')):
        w = w.replace(a, b)
    return re.sub(r'([^aeiouy])y(?=[aeiou])', r'\1iy', w)


def esperanto(word):
    """The author's Esperanto-based tongue in its own spelling. Kept: f and v. Changed: j, its y sound, is written y;
    the hatted ĥ and ŝ to h and s; ĝ to g; ĉ, ĵ, c and z to s."""
    w = word.lower()
    for a, b in (('ŭ', 'u'), ('ĥ', 'h'), ('ŝ', 's'), ('ĝ', 'g'), ('ĉ', 's'), ('ĵ', 's'), ('sh', 's'), ('x', 'h'), ('c', 's'), ('z', 's'), ('j', 'y')):
        w = w.replace(a, b)
    return re.sub(r'([ao])y(?![aeiou])', r'\1i', w)


def terlunan(word):
    """Terlunan's words are the framework's own; a coined word still in its Māori spelling takes f for wh, as the
    reference's fanga does."""
    return word.lower().replace('wh', 'f')


def plain(word):
    """A name from elsewhere in Latin letters with its sounds marked (zh, sh, ch, ts, z, j, q, long vowels)."""
    w = word.lower()
    for a, b in (('ā', 'a'), ('ī', 'i'), ('ū', 'u'), ('zh', 's'), ('sh', 's'), ('ch', 's'), ('ts', 's'), ('kh', 'h'), ('q', 'k'),
                 ('z', 's'), ('j', 'y'), ('x', 'h')):
        w = w.replace(a, b)
    return w


READERS = dict(kk=kazakh, mn=mongolian, qu=quechua, ar=arabic, sw=swahili, tl=tagalog, he=hebrew, ja=japanese, eo=esperanto,
               te=terlunan, la=plain)
PARENT_OF = dict(kk='Kazakh or Kyrgyz', mn='Mongolian', qu='Quechua', ar='Arabic or Tamazight', sw='Swahili', tl='Tagalog', he='Hebrew',
                 ja='Japanese', eo='Esperanto with Hebrew loans', te='Terlunan', la='another language')
# The reader each tongue's draft spellings pass through where no fuller parent form is given.
DRAFT_READER = dict(terlunan='te', olivolingvo='eo', sahran='ar', pwani='sw', dala='kk', runa='qu', satoyama='ja')


def read(spec):
    """'kk:көл' or 'kk:қара+mn:говь': each part read by its parent's rule sheet, unrepaired, for fusing."""
    out = []
    for part in spec.split('+'):
        code, word = part.split(':', 1)
        out.append(READERS[code](word))
    return out


def passed(spec):
    """A parent spec passed through whole: read, fused and repaired."""
    return repair(fuse(read(spec)))


# ------------------------------------------------------------------------------------------------ the parent forms
#
# A fuller parent form for a lexicon word: where the draft's spelling had lost a sound the framework can keep (Kazakh
# ng and its vowels, the Arabic throat stops, the Swahili syllabic nasal's class, Tagalog's glottal stop, Quechua ñ),
# or where the draft's word was not the parent's. Words not listed pass from the draft's spelling.

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
# The Tagalog-Hebrew tongue reads its heads as Tagalog and its modifiers as Hebrew.


def parent_spec(tongue, key, word, role):
    """The parent form a lexicon word passes from: PARENTS' where it gives one, else the draft's own spelling read by
    the tongue's rule sheet."""
    spec = PARENTS.get(tongue, {}).get(key)
    if spec:
        return spec
    if tongue == 'tagalhebrew':
        return ('tl:' if role == 'head' else 'he:') + word
    return f'{DRAFT_READER[tongue]}:{word}'


def lexicon(tongue):
    """A tongue's words passed through: {key: dict(role, draft, parent, read, passed, said, faults)}. read is the form
    names are fused from; passed is the word standing alone."""
    t = DN.TONGUES[tongue]
    out = {}
    for role, table in (('head', t['heads']), ('modifier', t['mods']), ('meaning', DN.CALQUE_WORDS.get(tongue, {}))):
        for key, word in table.items():
            if key in out and role == 'meaning':
                continue                                    # the tongue's own modifier already holds the meaning's word
            spec = parent_spec(tongue, key, word, role)
            got = read(spec)
            alone = repair(fuse(got))
            entry = dict(role=role, draft=word, parent=spec, read=fuse(got), passed=alone, said=said(alone), faults=faults(alone))
            note = PARENT_NOTES.get((tongue, key))
            if note:
                entry['note'] = note
            out[(role, key)] = entry
    return out


_LEX = {}


def _lex(tongue):
    if tongue not in _LEX:
        _LEX[tongue] = lexicon(tongue)
    return _LEX[tongue]


def word(tongue, key, role=None):
    """A tongue's word for a key as names are fused from it: its head, its modifier, or its word for an old name's
    meaning."""
    lex = _lex(tongue)
    for r in ((role,) if role else ('modifier', 'meaning', 'head')):
        if (r, key) in lex:
            return lex[(r, key)]['read']
    raise KeyError((tongue, key))


RENDAKU = (('k', 'g'), ('t', 'd'), ('h', 'b'))      # the voicing the framework can show; s has no voiced partner in it


def compose(tongue, head, mod):
    """A name from a head and a modifier key, as district_names.compose joins them, from the passed words."""
    join = DN.TONGUES[tongue]['join']
    h, m = word(tongue, head, 'head'), word(tongue, mod)
    cap = DN._cap
    if join == 'fused':
        return cap(repair(fuse([m, h])))
    if join == 'rendaku':
        if head not in ('lake', 'lakes', 'basin', 'desert'):
            for a, b in RENDAKU:
                if h.startswith(a):
                    h = b + h[len(a):]
                    break
        return cap(repair(fuse([m, h])))
    h, m = repair(h), repair(m)
    if join == 'olivo':
        return f"H'{cap(h)} {cap(m)}"
    if join == 'tagalog':
        return f'{cap(h)}ng {cap(m)}' if h[-1] in VOWELS else f'{cap(h)} na {cap(m)}'
    return f'{cap(h)} {cap(m)}'


# The names table's pins as a head and a modifier (or a meaning's word) of their tongue; district_names.compose makes
# each pin's present spelling from the same keys (a test holds this), and compose above its passed form.
PINS = {
    'sea:near-side sea': ('terlunan', 'sea', 'rain'), 'sea:southern near-side seas': ('dala', 'sea', 'cloud'),
    'sea:South Pole–Aitken Sea': ('olivolingvo', 'sea', 'salt'), 'sea:Mare Humboldtianum': ('dala', 'sea', 'star'),
    'sea:Mare Moscoviense': ('runa', 'lake', 'high'), 'sea:Mare Orientale': ('pwani', 'sea', 'east'),
    'sea:Smythii–Marginis Sea': ('tagalhebrew', 'sea', 'edge'),
    'river:-22.9,173.3': ('runa', 'river', 'great'), 'river:47.5,49.7': ('terlunan', 'river', 'blue'),
    'river:2.3,177.7': ('runa', 'river', 'rain'), 'river:21.3,105.0': ('satoyama', 'river', 'deep'),
    'river:15.0,287.7': ('pwani', 'river', 'green'), 'river:-13.4,253.1': ('tagalhebrew', 'river', 'high'),
    'river:-23.0,192.3': ('runa', 'river', 'mist'), 'river:25.8,103.0': ('satoyama', 'river', 'new'),
    'river:3.6,302.2': ('dala', 'river', 'light'), 'river:-45.1,230.0': ('pwani', 'river', 'stone'),
    'river:57.1,111.8': ('sahran', 'river', 'gold'), 'river:16.8,13.5': ('terlunan', 'river', 'flower'),
    'river:18.1,184.8': ('runa', 'river', 'wind'), 'river:3.3,156.4': ('runa', 'river', 'cloud'),
    'river:8.8,146.9': ('runa', 'river', 'star'), 'river:1.4,302.5': ('pwani', 'river', 'leaf'),
    'river:-14.1,311.5': ('pwani', 'river', 'twin'), 'river:-50.2,131.6': ('olivolingvo', 'river', 'long'),
    'river:12.3,346.0': ('terlunan', 'river', 'black'), 'river:47.4,186.0': ('olivolingvo', 'river', 'old'),
}
# Pins that are not a head and a modifier: fused from named parts, or the author's own, which stand.
PINS_BY_PARTS = {'river:19.9,184.6': 'te:rere+te:awa'}
PINS_STAND = ('feature:Sinus Iridum', 'feature:Montes Jura', 'feature:Mare Ingenii')


def draft_compose(tongue, head, mod):
    """district_names.compose with the word for an old name's meaning allowed as the modifier, as its calque does."""
    t = DN.TONGUES[tongue]
    if mod in t['mods']:
        return DN.compose(tongue, head, mod)
    saved = t['mods']
    t['mods'] = dict(saved, **{mod: DN.CALQUE_WORDS[tongue][mod]})
    try:
        return DN.compose(tongue, head, mod)
    finally:
        t['mods'] = saved


def pins(table=None):
    """Every pinned name with its passed form: [dict(anchor, name, tongue, passed, said, changed, how)]."""
    table = DN.load_table() if table is None else table
    out = []
    for anchor, pin in table.items():
        name = pin['name']
        if anchor in PINS_STAND:
            new, how = name, "the author's own name, which stands"
        elif anchor in PINS_BY_PARTS:
            new, how = DN._cap(passed(PINS_BY_PARTS[anchor])), PINS_BY_PARTS[anchor]
        elif anchor in PINS:
            tongue, head, mod = PINS[anchor]
            new, how = compose(tongue, head, mod), f'{tongue}: {mod} + {head}'
        else:
            new, how = None, 'no rule: named by hand since this table was made'
        out.append(dict(anchor=anchor, name=name, tongue=pin.get('tongue'), passed=new, said=said(new) if new else None,
                        changed=bool(new) and new != name, how=how,
                        faults=faults(new) if new and anchor not in PINS_STAND else []))
    return out


def collisions(tongue):
    """Pairs of a tongue's words that the draft kept apart and the pass makes one, among its heads and among its
    modifiers and meanings' words."""
    lex = _lex(tongue)
    out = []
    for roles in (('head',), ('modifier', 'meaning')):
        seen = {}
        for (role, key), e in lex.items():
            if role not in roles:
                continue
            other = seen.setdefault(e['passed'], (key, e['draft']))
            if other[0] != key and other[1].lower() != e['draft'].lower():
                out.append((other[0], key, e['passed']))
    return out


def record_names(record):
    """How many of a districts record's generated names the pass changes: (generated names found again from their
    head and modifier, of those changed)."""
    found = changed = 0
    for u in record.get('regions', []) + record.get('districts', []) + record.get('coastal_seas', []):
        if u.get('name_source') != 'generated' or not u.get('name_head') or not u.get('name_mod'):
            continue
        if DN.compose(u['tongue'], u['name_head'], u['name_mod']) != u['name']:
            continue                                        # a name with a second word added to tell it from another
        found += 1
        changed += compose(u['tongue'], u['name_head'], u['name_mod']) != u['name']
    return found, changed


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(log=print):
    t0 = time.time()
    tongues = {}
    for tongue, t in DN.TONGUES.items():
        lex = lexicon(tongue)
        words = [dict(key=key, **e) for (role, key), e in lex.items()]
        tongues[tongue] = dict(label=t['label'], words=words, changed=sum(e['passed'] != e['draft'].lower() for e in lex.values()),
                               faults=sum(bool(e['faults']) for e in lex.values()), fallen_together=collisions(tongue),
                               from_fuller_parent_forms=sum(1 for (role, key) in lex if key in PARENTS.get(tongue, {})))
        log(f"{t['label']:36s} {len(words):3d} words, {tongues[tongue]['changed']:3d} change, {tongues[tongue]['faults']} against the "
            f"framework, {len(tongues[tongue]['fallen_together'])} fallen together")
    pinned = pins()
    record = json.loads(DN.RECORD.read_text()) if DN.RECORD.exists() else {}
    found, changed = record_names(record)
    log(f"pins: {len(pinned)}, {sum(p['changed'] for p in pinned)} change; the draft record's generated names: {changed} of {found} change")
    return dict(
        schema=SCHEMA, producer='research.studies.atlas_round.lunar_sounds',
        inputs={str(p.relative_to(ROOT)): _sha(p) for p in (HERE / 'district_names.py', DN.PLACE_NAMES, DN.RECORD) if p.exists()},
        evidence=("A rule sheet and its results, nothing in it validated: the framework is the sound section of the author's Terlunan "
                  "reference (v0.6, section 2, outside the repository); each parent language's rules are choices made on 20 November 2026; "
                  "the parent words are written from general knowledge of Kazakh, Kyrgyz, Mongolian, Quechua, Arabic, Tamazight, Swahili, "
                  "Tagalog, Hebrew, Japanese and Esperanto, and no dictionary or speaker has checked them. No name is changed by this record: "
                  "district_names.py and place_names.json stand as they are."),
        reading_rule=("tongues gives, for every word of district_names.py's lexicons, its key, its role (head, modifier, or the word for an "
                      "old name's meaning), the draft's spelling, the parent form it passes from (a code for the parent's rule sheet, a colon "
                      "and the word in the parent's letters or the draft's), the form names are fused from (read), the word standing alone "
                      "(passed), its syllables with the stressed one in capitals (said) and any fault against the framework; changed counts "
                      "the words whose spelling moves, fallen_together lists pairs the draft kept apart and the pass makes one. pins gives "
                      "every name in place_names.json with its passed form. record gives how many of the draft record's generated names "
                      "the pass would change."),
        framework=dict(vowels=list(VOWELS), consonants=list(CONSONANTS), diphthongs=list(DIPHTHONGS), admitted_in_established_words=['f', 'v'],
                       syllable='(C)V(C)', stress='the syllable before the last, or the rightmost long vowel'),
        rule_sheets={code: (fn.__doc__ or '').strip().replace('\n    ', ' ') for code, fn in READERS.items()},
        tongues=tongues, pins=pinned,
        record=dict(generated_names_found=found, changed=changed),
        seconds=round(time.time() - t0, 2))


def main(argv=None) -> int:
    record = build()
    OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
