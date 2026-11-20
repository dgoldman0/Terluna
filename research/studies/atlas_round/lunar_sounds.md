# The Moon's shared sounds: every name passes through them

A place name on the Open Moon keeps its parent language's sounds to a degree and comes out in the Moon's own phonetic framework (register, 20 November 2026). [lunar_sounds.py](lunar_sounds.py) is that passage: one framework, a rule sheet for each parent language, and every word of the eight tongues' lexicons passed through. The record is [results/lunar_sounds.json](results/lunar_sounds.json).

Of the lexicons' 789 words 324 change their spelling. Every passed word stands in the framework, and no two words that a tongue's draft kept apart fall together. Of the 31 names pinned in the names table 15 would change, among them six of the seven seas; of the 658 generated names in the page's draft record 401 would change. **No name is changed yet**: [district_names.py](district_names.py) and [place_names.json](place_names.json) stand as they are until the passed forms are adopted.

*A rule sheet and its results, nothing in it validated. The rules are choices of one day, the parent words are written from general knowledge of the languages, and no dictionary or speaker has checked them.*

## The framework

The framework is the sound section of the author's Terlunan reference (v0.6, section 2, kept outside the repository), taken as the sounds all the Moon's tongues share:

- five vowels, *a e i o u*, a long vowel written double;
- the consonants *p b t d k g m n ng s h l r w y* and the glottal stop, with *f* and *v* admitted in established words; no *sh*, *z*, emphatic or affricate series;
- syllables of the shape (C)V(C), with *ai, au, oi, ou* as one nucleus: a consonant may close a syllable and meet the next one's opening consonant, but no syllable opens or closes on two;
- no vowel added at the end of a borrowed word;
- where two parts of a name are fused and the same consonant meets itself, one remains;
- stress on the syllable before the last, unless a long vowel takes it.

What the framework cannot hold is repaired as the reference repairs its own loans. A cluster that opens a word gets a vowel inside it, *i* before *y* and *u* before *w* (the reference's *liyong* and *kuwo*) and otherwise a copy of the next vowel. A cluster that closes a word is opened with a copy of the vowel before it. A doubled consonant becomes one.

## What each parent keeps, and what changes

A sound the framework has is kept. A sound it lacks goes to its nearest neighbour inside the framework, and where there is a choice the rule takes the one that keeps most of the parent.

| Tongue | Parent | Kept | Changed | Draft to passed |
|---|---|---|---|---|
| Dala | Kazakh, Kyrgyz, Mongolian | *ng*, which the draft's spelling had often lost; final consonants; the rounding harmony of speech (a short high vowel is *u* after *o* or *u*); Mongolian's long vowels | *q*, *gh* to *k*, *g*; *sh*, *z* to *s*; *j* to *y*; *ö*, *ü* to *o*, *u* | tengiz, *tengis*; juldyz, *yuldus*; jana, *yanga*; kyzyl, *kisil*; bult, *bulut*; uiky, *uyku* |
| Runa | Quechua | The lowering of *i* and *u* to *e* and *o* beside the uvular *q*, which Quechua speaks and does not write | *q* to *k*; *ch*, *sh* to *s*; *ll*, *ñ* to *l*, *n*; the aspirated and ejective stops to the plain ones | qucha, *kosa*; urqu, *orko*; llaqta, *lakta*; quyllur, *koilur*; wayra, *waira* |
| Sahran | Arabic, Tamazight | The two throat stops, as the framework's glottal stop; *f* | The emphatics to their plain partners and *q* to *k*; *kh* to *h*, *gh* to *g*; *sh*, *z* to *s*; *j* to *g*; long vowels to short; closing clusters opened | Bahr, *bahar*; Jazira, *gasira*; Shams, *samas*; Amik, *'amik*; Najm, *nagam*; Tamazirt, *tamasirit* |
| Pwani | Swahili | Its five vowels, the framework's own; *f*, *v*; *ng* | *sh*, *ch*, *z* to *s*; *j* to *y*; a syllabic nasal gets back the vowel it lost; a consonant and a glide that open a word get a vowel between | Mto, *muto*; Nchi, *insi*; Mvua, *imvua*; Kijani, *kiyani*; Pwani, *puwani*; Ziwa, *siwa* |
| Tagalog-Hebrew | Tagalog heads, Hebrew modifiers | Tagalog whole, its sounds being the framework's, with its glottal stop written between vowels; Hebrew's *f*, *v* | Hebrew *sh*, *ts*, *z* to *s*; *kh* to *h* | Look, *lo'ok*; Shemesh, *semes*; Erets, *eres*; Katzeh, *kaseh*; Zahav, *sahav* |
| Satogo | Japanese | Its five vowels and long vowels | The spelling that follows its own sound system: *shi*, *chi*, *tsu*, *fu* are *si*, *ti*, *tu*, *hu*; *z*, *j* to *s* | shima, *sima*; tsuchi, *tuti*; mizu, *misu*; fuka, *huka*; o, *oo* |
| Proto-Olivolingvo | Esperanto with Hebrew loans | *f*, *v*; the article *h'* as the author writes it | *j*, its *y* sound, written *y*; *ĥ*, *ŝ* to *h*, *s*; *ĝ* to *g*; *c*, *z* to *s*; opening clusters repaired | jamo, *yamo*; ĥadaŝa, *hadasa*; kratero, *karatero*; profunda, *porofunda*; erco, *erso* |
| Terlunan | itself | Everything: its words are the framework's own | One coined word still in its Māori spelling | whakapono, *fakapono* |

## The choices, and what else was possible

Each rule is a choice that can be turned. These change most how the names look.

**The sibilants fall together in *s*.** The reference does the same with its Sumerian, Akkadian and Hakka words. It makes the names *s*-heavy, and it is why two parent words were changed for others (below): once *z* is *s*, the draft's Arabic word for darkness would have read as the word for peace.

**The sound written *j* goes three ways, by parent.** In the Kazakh-rooted words it becomes *y*, the sound Kazakh *j* answers to across the Turkic languages (*yilan*, *yuldus*). In the Arabic-rooted words it becomes *g*, as Cairo says it (*gasira*, *gadid*). In the Swahili-rooted words it becomes *y*. The other candidate everywhere is *s*, which the reference uses for an affricate and which would pile a fourth sound onto *s*.

**Long vowels stay where the parent's own spelling shows them** (Mongolian, Japanese, Terlunan) and go where it does not (Arabic, Hebrew), as the reference drops the length of its Akkadian words.

**The throat stops of Arabic stay** as the glottal stop the framework has: *'amik*, *ma'ruf*, *rabi'*. Dropping them is the other choice.

**Closing clusters open with a copy vowel**: *bahar*, *nagam*, *ramal*, *bulut*, *eles*. Several living forms agree (Maltese *baħar*, Kyrgyz *bulut*).

**Quechua keeps its lowered vowels.** Quechua writes *qucha* and says the vowel beside the uvular as *o*. With *q* gone to *k* the pass writes what is said: *kosa*, *orko*, *hanak*.

**Swahili's syllabic nasal gets its vowel back by its class**: *mu-* in the tree-and-river class (*muto*, *musitu*) and *i-* before the nasal elsewhere (*imvua*, *insi*, *inyota*), as its sister languages have it.

**Japanese takes the spelling of its own sound system**: *si*, *ti*, *tu*, *hu*. Sequential voicing shows where the framework can show it (*Kurogawa*) and not on *s* (*Kurosima*).

**The author's own tongue is passed like the rest**, as a proposal. Its opening clusters are repaired (*karatero*, *kilifo*, *porofunda*, *sotona*), which changes its look more than any other rule. Its article *h'* is left as the author writes it; how it is said inside the framework is the author's to settle. One homograph follows: its word for dawn, *ŝaĥara*, passes to *sahara*.

Six parent forms put right more than sounds, and each is noted in the record: the draft's Dala plural for lakes was Uzbek (*koldor*, the Kyrgyz, replaces it); its words for dream and forget would have fallen together with salt and hope (*arman* and *unut* replace them); the Sahran word for darkness is *'atma*; the Pwani word for rough had a slip; and a throat stop keeps the two vowels of the Hebrew word for sorrow apart.

## What it does to the lexicons and the pins

| Tongue | Words | Change | From a fuller parent form |
|---|---|---|---|
| Terlunan | 99 | 1 | 0 |
| Proto-Olivolingvo | 99 | 35 | 0 |
| Sahran | 99 | 57 | 14 |
| Tagalog-Hebrew | 99 | 28 | 2 |
| Pwani | 99 | 39 | 5 |
| Dala | 99 | 66 | 99 |
| Runa | 97 | 68 | 4 |
| Satogo | 98 | 30 | 1 |

A word passes from the draft's spelling, or from a fuller parent form where that spelling had lost a sound the framework can keep. Dala's words are all given in Kazakh, Kyrgyz and Mongolian letters, because the draft's Latin spelling could not tell its vowels apart.

The names pinned in the names table, found again from their head and modifier and joined from the passed words:

| Pinned name | Passed | Said |
|---|---|---|
| Bultengiz | Bultengis | bul.TE.ngis |
| H'Jamo Melaĥa | H'Yamo Melaha | h'YA.mo me.LA.ha |
| Juldyztengiz | Yuldustengis | yul.dus.TE.ngis |
| Hanaqucha | Hanakosa | ha.na.KO.sa |
| Bahari Mashariki | Bahari Masariki | ba.HA.ri ma.sa.RI.ki |
| Dagat na Katzeh | Dagat na Kaseh | DA.gat na KA.seh |
| Fukagawa | Hukagawa | hu.ka.GA.wa |
| Mto Kijani, Mto Jiwe, Mto Jani, Mto Pacha | Muto Kiyani, Muto Yiwe, Muto Yani, Muto Pasa | MU.to ki.YA.ni |
| Pukutaymayu | Pukutaimayu | pu.ku.tai.MA.yu |
| Jarykgol | Yarikgol | ya.RIK.gol |
| Wayramayu | Wairamayu | wai.ra.MA.yu |
| Quyllurmayu | Koilurmayu | koi.lur.MA.yu |

Sixteen stay: Uamoana; the rivers Hatunmayu, Ninamuawa, Paramayu, Ilog na Gavoah, Niigawa, Rereawa, Asif Dahab, Tiareawa, Puyumayu, Pooiriawa, H'Naharo Longa and H'Naharo Atika; and the author's own Anuanuafanga, Jura and Ingenii, which stand whatever the pass would say.

## What is not done

- **No name is changed.** Adopting the pass means putting the passed words into district_names.py's lexicons and the passed names into place_names.json, renaming the draft record (seconds), filling the regions' texts and the waterways again, and building and publishing the page.
- **The naming trial's two regions** are on the research drive with their own parent words, made before this module and by the same rules.
- **The lexicons themselves are the draft's.** The pass changes how a word sounds, not which word a tongue has; the reworked way of making names chooses the words.

## Limits

- The framework is read from the Terlunan reference, which sets it for Terlunan. Taking it as the sounds every lunar tongue shares is this study's reading of the register's row.
- The parent words and the rules for eleven languages are written from general knowledge. A check against dictionaries, and of the Kazakh rounding harmony and the Quechua lowering against descriptions of those languages, has not been made.
- A passed word can look like a word of another language: Dala's *sol* (desert, thirst), the author's tongue's *sahara* (dawn), the Hebrew *kesef* from the word for foam. The module checks only that two words of one tongue do not fall together.
- Stress moves with the pass: the framework stresses the syllable before the last, where Kazakh stresses the last and Arabic follows the weight of its syllables.

## Running it

```sh
python -m research.studies.atlas_round.lunar_sounds      # under a second: results/lunar_sounds.json and a summary by tongue
```

[test_lunar_sounds.py](test_lunar_sounds.py) holds the framework's checks and repairs on words that fit and words that do not, each parent's rules on words worked by hand, every word of every lexicon and every joined name against the framework, the pins found again from their keys, and that building the record changes no name.
