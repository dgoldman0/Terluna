# The Moon's shared sounds: every name passes through them

A place name on the Open Moon keeps its parent language's sounds to a degree and comes out in the Moon's own phonetic framework (register, 20 November 2026). [lunar_sounds.py](lunar_sounds.py) is that passage: one framework and a rule sheet for each parent language. [district_names.py](district_names.py) writes every word of the eight tongues as its parent language spells it and passes them through when it loads, so every name is made from the passed words. The record is [results/lunar_sounds.json](results/lunar_sounds.json).

Of the lexicons' 789 words 324 change their spelling. Every passed word and every name joined from them stands in the framework, and no two different parent words of one tongue fall together. The passed forms were adopted the same day: 15 of the 31 names pinned in [place_names.json](place_names.json) are respelled, six of the seven seas among them, and 446 of the 755 names in the page's draft record.

*A rule sheet and its results, nothing in it validated. The rules are choices of one day, the parent words are written from general knowledge of the languages, and no dictionary or speaker has checked them.*

## The framework

The framework is the sound section of the author's Terlunan reference (v0.6, section 2, kept outside the repository), taken as the sounds all the Moon's tongues share:

- five vowels, *a e i o u*, a long vowel written double;
- the consonants *p b t d k g m n ng s h l r w y* and the glottal stop, with *f* and *v* admitted in established words; no *sh*, *z*, emphatic or affricate series;
- syllables of the shape (C)V(C), with *ai, au, oi, ou* as one nucleus: a consonant may close a syllable and meet the next one's opening consonant, but no syllable opens or closes on two;
- no vowel added at the end of a borrowed word;
- where two parts of a name are fused and the same consonant meets itself, one remains;
- stress on the syllable before the last, unless a long vowel takes it.

What the framework cannot hold is repaired as the reference repairs its own loans. A cluster that opens a word gets a vowel inside it, *i* before *y* and *u* before *w* (the reference's *liyong* and *kuwo*) and otherwise a copy of the next vowel. A cluster that closes a word is opened with a copy of the vowel before it. A doubled consonant becomes one. Where a parent language has its own vowel for opening a cluster, that vowel is used: Hebrew words open a cluster with *e*.

## What each parent keeps, and what changes

A sound the framework has is kept. A sound it lacks goes to its nearest neighbour inside the framework, and where there is a choice the rule takes the one that keeps most of the parent.

| Tongue | Parent | Kept | Changed | As the parent spells it, and passed |
|---|---|---|---|---|
| Dala | Kazakh, Kyrgyz, Mongolian | *ng*; final consonants; the rounding harmony of speech (a short high vowel is *u* after *o* or *u*); Mongolian's long vowels | *q*, *gh* to *k*, *g*; *sh*, *z* to *s*; *j* to *y*; *ö*, *ü* to *o*, *u* | теңіз, *tengis*; жұлдыз, *yuldus*; жаңа, *yanga*; қызыл, *kisil*; бұлт, *bulut*; ұйқы, *uyku* |
| Runa | Quechua | The lowering of *i* and *u* to *e* and *o* beside the uvular *q*, which Quechua speaks and does not write | *q* to *k*; *ch*, *sh* to *s*; *ll*, *ñ* to *l*, *n*; the aspirated and ejective stops to the plain ones | qucha, *kosa*; urqu, *orko*; llaqta, *lakta*; quyllur, *koilur*; wayra, *waira* |
| Sahran | Arabic, Tamazight | The two throat stops, as the framework's glottal stop; *f* | The emphatics to their plain partners and *q* to *k*; *kh* to *h*, *gh* to *g*; *sh*, *z* to *s*; *j* to *g*; long vowels to short; closing clusters opened | baḥr, *bahar*; jazīra, *gasira*; shams, *samas*; ʿamīq, *'amik*; najm, *nagam*; tamazirt, *tamasirit* |
| Pwani | Swahili | Its five vowels, the framework's own; *f*, *v*; *ng* | *sh*, *ch*, *z* to *s*; *j* to *y*; a syllabic nasal gets back the vowel it lost; a consonant and a glide that open a word get a vowel between | mto, *muto*; nchi, *insi*; mvua, *imvua*; kijani, *kiyani*; pwani, *puwani*; ziwa, *siwa* |
| Tagalog-Hebrew | Tagalog heads, Hebrew modifiers | Tagalog whole, its sounds being the framework's, with its glottal stop written between vowels; Hebrew's *f*, *v* | Hebrew *sh*, *ts*, *z* to *s*; *kh* to *h*; a Hebrew cluster opened with *e* | look, *lo'ok*; shemesh, *semes*; erets, *eres*; katzeh, *kaseh*; zman, *seman* |
| Satogo | Japanese | Its five vowels and long vowels | The spelling that follows its own sound system: *shi*, *chi*, *tsu*, *fu* are *si*, *ti*, *tu*, *hu*; *z*, *j* to *s* | shima, *sima*; tsuchi, *tuti*; mizu, *misu*; fuka, *huka*; ō, *oo* |
| Proto-Olivolingvo | Esperanto with Hebrew loans | *f*, *v*; the article *h'* as the author writes it | *j*, its *y* sound, written *y*; *ĥ*, *ŝ* to *h*, *s*; *ĝ* to *g*; *c*, *z* to *s*; opening clusters opened by the framework's own repair for now | jamo, *yamo*; ĥadaŝa, *hadasa*; kratero, *karatero*; profunda, *porofunda*; erco, *erso* |
| Terlunan | itself | Everything: its words are the framework's own | One coined word still in its Māori spelling | whakapono, *fakapono* |
| No lexicon yet: tried on one region's names | Louisiana French, with the Creoles of Louisiana and Haiti | *f*, *v*; the *h* that Louisiana still sounds; a consonant closing a syllable before the next one's opening consonant | The hushing sounds and *z* to *s*; *gn* to *n* and *y*; the rounded front vowels unrounded; a nasal vowel to its vowel with *n*; a final stop loses the *l* or *r* after it | bleu, *bele*; lune, *lin*; lagniappe, *lanyap*; grand lac, *garanlak*; ouaouaron, *wawaron*; sable, *sab* |

| No lexicon yet: tried on one region's names | Somali | Long vowels, written double; the glottal stop | *c*, the deep throat sound, to the glottal stop; *x*, *kh* to *h*; *q* to *k*; *sh* to *s*; *dh* to *d*; *j* to *y* | xeeb, *heeb*; cidla, *'idla*; qabow, *kabow*; gacan, *ga'an* |
| No lexicon yet: tried on one region's names | Uyghur | *ng*; final consonants | *q*, *gh* to *k*, *g*; *x* to *h*; *ch*, *sh*, *zh*, *z* to *s*; *j* to *y*; *ö*, *ü*, *é* to *o*, *u*, *e* | qumluq, *kumluk*; köl, *kol*; quruq tagh, *kuruktag* |
| No lexicon yet: tried on one region's names | Setswana | *ng*; its open syllables | *kg* to *k* and *g*, its throat sound, to *h*; the breathed stops to plain ones; *tl* to *t*; *ts*, *š* to *s* | Kgalagadi, *kalahadi*; thaba, *taba*; kgolo, *kolo* |
| No lexicon yet: tried on one region's names | English of the American West | *f*, *v*, *h*, the *r* the West sounds; the glides *ai*, *au*, *oi*, *ou* | *th* to *t* or *d*; the hushing sounds, *z*, *ch*, *j* to *s*; the short vowels to the nearest of the five | rimrock, *rimrak*; lonesome, *lounsem*; Dyson, *daisen* |
| No lexicon yet: tried on one region's names | Spanish of New Mexico | Its five vowels | *ll* to *y*; *j* and the soft *g* to *h*; *z*, the soft *c* and *ch* to *s*; *v* to *b*; the silent *h* dropped | jornada, *hornada*; hoyo, *oyo*; grande, *garande*; verde, *berde* |

The Louisiana French and English sheets take each word as it is said, written in the phonetic alphabet's letters, because their spellings carry letters that are not sounded. These six sheets were made as regions were named on the map of parent cultures (register, 20 November): Louisiana French for the second naming trial, the other five for the first batch of regions after it. Tamasheq and Hassaniya Arabic pass by the Arabic and Tamazight sheet. No lexicon uses the six yet.

## The choices, and what else was possible

Each rule is a choice that can be turned. These change most how the names look.

**The sibilants fall together in *s*.** The reference does the same with its Sumerian, Akkadian and Hakka words. It makes the names *s*-heavy, and it is why two parent words were changed for others (below): once *z* is *s*, the first draft's Arabic word for darkness would have read as the word for peace.

**The sound written *j* goes three ways, by parent.** In the Kazakh-rooted words it becomes *y*, the sound Kazakh *j* answers to across the Turkic languages (*yilan*, *yuldus*). In the Arabic-rooted words it becomes *g*, as Cairo says it (*gasira*, *gadid*). In the Swahili-rooted words it becomes *y*. The other candidate everywhere is *s*, which the reference uses for an affricate and which would pile a fourth sound onto *s*.

**Long vowels stay where the parent's own spelling shows them** (Mongolian, Japanese, Terlunan) and go where it does not (Arabic, Hebrew), as the reference drops the length of its Akkadian words.

**The throat stops of Arabic stay** as the glottal stop the framework has: *'amik*, *ma'ruf*, *rabi'*. Dropping them is the other choice.

**Closing clusters open with a copy vowel**: *bahar*, *nagam*, *ramal*, *bulut*, *eles*. Several living forms agree (Maltese *baħar*, Kyrgyz *bulut*).

**Quechua keeps its lowered vowels.** Quechua writes *qucha* and says the vowel beside the uvular as *o*. With *q* gone to *k* the pass writes what is said: *kosa*, *orko*, *hanak*.

**Swahili's syllabic nasal gets its vowel back by its class**: *mu-* in the tree-and-river class (*muto*, *musitu*) and *i-* before the nasal elsewhere (*imvua*, *insi*, *inyota*), as its sister languages have it.

**Japanese takes the spelling of its own sound system**: *si*, *ti*, *tu*, *hu*. Sequential voicing shows where the framework can show it (*Kurogawa*) and not on *s* (*Kurosima*).

**Louisiana French unrounds its front vowels as the Creoles do, and writes a nasal vowel with *n*.** The Creoles of Louisiana and Haiti say *lalin* for *la lune* and *ble* for *bleu*, so the pass gives *lin* and *bele*; keeping *u* for the French *u* is the other choice. A nasal vowel keeps its nasal as *n* (*garan*, *wawaron*); dropping it is the other choice. A final stop loses the *l* or *r* after it, as Louisiana speech does (*sab* for *sable*), which spares the framework a closing cluster.

**Proto-Olivolingvo takes the letter and sibilant changes, and the rest waits for the local languages** (the author's decisions, 20 November). It is the proposed common language of Earth, as Terlunan is the Moon's, and it has no second parent of its own: on the Moon a word of it is shaped by the local language that takes it up, so how its clusters open and how its article *h'* is said depend on that language's other parent. Until the local languages are set out its clusters take the framework's own repair (*karatero*, *kilifo*, *porofunda*, *sotona*), and its article is written as the author writes it with no saying given. A first reading that took Hebrew for that other parent everywhere was withdrawn the same day. One homograph follows from the sibilants: its word for dawn, *ŝaĥara*, passes to *sahara*. Hebrew's own words, in the Tagalog-Hebrew tongue, open a cluster with Hebrew's *e* (*seman*, *setav*, *bedidut*).

Six parent forms put right more than sounds, and each is noted in the record: the first draft's Dala plural for lakes was Uzbek (*koldor*, the Kyrgyz, replaces it); its words for dream and forget would have fallen together with salt and hope (*arman* and *unut* replace them); the Sahran word for darkness is *'atma*; the Pwani word for rough had a slip; and a throat stop keeps the two vowels of the Hebrew word for sorrow apart.

## What it did to the lexicons and the names

| Tongue | Words | Changed by the pass | Given in a fuller parent form |
|---|---|---|---|
| Terlunan | 99 | 1 | 0 |
| Proto-Olivolingvo | 99 | 35 | 0 |
| Sahran | 99 | 57 | 14 |
| Tagalog-Hebrew | 99 | 28 | 2 |
| Pwani | 99 | 39 | 5 |
| Dala | 99 | 66 | 99 |
| Runa | 97 | 68 | 4 |
| Satogo | 98 | 30 | 1 |

A word passes from its table's spelling, or from a fuller parent form where that spelling had lost a sound the framework can keep. Dala's words are all given in Kazakh, Kyrgyz and Mongolian letters, because the Latin spelling could not tell its vowels apart.

The names pinned in the names table, each found again from its head and modifier and made from the passed words:

| Was | Is | Said |
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

Sixteen stand as they were: Uamoana; the rivers Hatunmayu, Ninamuawa, Paramayu, Ilog na Gavoah, Niigawa, Rereawa, Asif Dahab, Tiareawa, Puyumayu, Pooiriawa, H'Naharo Longa and H'Naharo Atika; and the author's own Anuanuafanga, Jura and Ingenii, which stand as written whatever the pass would say.

In the draft record 424 generated names, 16 calques and the six seas changed. Each was found again from its head and modifier with the lexicons as they were and made again with the lexicons as they are, so no name was drawn afresh and no two names fell together. The regions' fact sheets and the waterways carry the new spellings with every figure as it was.

## What is left

- **The scheme of tongues is to be set out again** (register, 20 November). A region's names are to fit the distribution of the related parent cultures of Earth; most names have a specific parent language; within the broader regions local languages take those parent cultures and shape the Moon's dialects and languages; Terlunan is the Moon's common language and Proto-Olivolingvo the proposed common language of Earth; and the local languages are more diverse than the eight tongues here, of which one, Proto-Olivolingvo, now stands as a region's tongue where it should be a contributing layer.
- **Rule sheets for the other languages of the map of parent cultures** are made as their regions are named; six are here so far. The map and both naming trials are working products on the research drive.
- **The redraw of 20 November and the naming trial's two regions** are on the research drive; the redraw's names are the generator's and its record keeps their earlier spellings.
- **The lexicons themselves are the first draft's.** The pass changes how a word sounds, not which word a tongue has; the reworked way of making names chooses the words.
- **Older records keep the spellings of their date**: the register's earlier rows and the lake rule's two records.

## Limits

- The framework is read from the Terlunan reference, which sets it for Terlunan. Taking it as the sounds every lunar tongue shares is this study's reading of the register's row.
- The parent words and the rules for seventeen languages are written from general knowledge. The words of one batch of names in Somali, Setswana, Uyghur and Tamasheq were looked up in web dictionaries on 20 November: most were found, and four names rest on words that were not. The rule sheets are sound rules and do not depend on that check. A check against dictionaries, and of the Kazakh rounding harmony and the Quechua lowering against descriptions of those languages, has not been made.
- A passed word can look like a word of another language: Dala's *sol* (desert, thirst), the author's tongue's *sahara* (dawn), the Hebrew *kesef* from the word for foam. The tests check only that two words of one tongue do not fall together.
- Stress moves with the pass: the framework stresses the syllable before the last, where Kazakh stresses the last and Arabic follows the weight of its syllables.

## Running it

```sh
python -m research.studies.atlas_round.lunar_sounds      # under a second: results/lunar_sounds.json and a summary by tongue
```

[test_lunar_sounds.py](test_lunar_sounds.py) holds the framework's checks and repairs on words that fit and words that do not, each parent's rules on words worked by hand, every word of every lexicon and every joined name against the framework, the names table against the names its words make, and the draft record's names against theirs.
