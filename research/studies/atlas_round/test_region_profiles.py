"""Tests of the regions' fact sheets and the texts written from them (region_profiles.py)."""

import json

import pytest

from research.studies.atlas_round import district_names as DN
from research.studies.atlas_round import districts as D
from research.studies.atlas_round import region_profiles as RP


def test_display_forms():
    nb = lambda x: x.replace(' ', '\u00a0', 1) if any(x.endswith(k) for k in ('km²', 'm³/s', ' m', '°C', 'mm', 'GW', ' km')) else x
    assert RP.display('area', 677288) == nb('677,000 km²')
    assert RP.display('area', 4259) == nb('4,300 km²')
    assert RP.display('lakes', 651) == nb('650')
    assert RP.display('lakes', 17) == nb('17')
    assert RP.display('lake_share', 0.27) == nb('27%')
    assert RP.display('lake_share', 0.004) == nb('under 1%')
    assert RP.display('earth_share', 0.998) == nb('99%')
    assert RP.display('river', 7095) == nb('7,100 m³/s')
    assert RP.display('river', 688) == nb('690 m³/s')
    assert RP.display('rain', 1639.77) == nb('1,640 mm')
    assert RP.display('air', 26.65) == nb('27 °C')
    assert RP.display('swing', 2.63) == nb('2.6 °C')
    assert RP.display('earth', -1.2) == nb('−1°')
    assert RP.display('highest', 5312) == nb('5,300 m')
    assert RP.display('centre', [8.385, -0.072]) == nb('8° N, 0°')
    assert RP.display('centre', [-56.2, 186.0]) == nb('56° S, 174° W')


def test_gloss_recovers_each_tongues_words():
    for tongue, t in DN.TONGUES.items():
        head, mod = next(iter(t['heads'])), next(iter(t['mods']))
        name = DN.compose(tongue, head, mod)
        g = RP.gloss(dict(name=name, tongue=tongue, name_meaning=''))
        assert g.startswith(f'{mod} {head}'.title()), (tongue, name, g)
        assert t['heads'][head].lower() in g and t['mods'][mod].lower() in g
    assert RP.gloss(dict(name='Uamoana', tongue='terlunan', name_meaning='Rain Sea')) == 'Rain Sea'


def test_fill_follows_renames_and_stale_lists_moved_facts():
    entry = dict(text='Some {lakes} lakes drain to [Uamoana], under {cloud} cloud.',
                 places={'Uamoana': 'sea:near-side sea'},
                 rests_on={'main_zone': 'wet', 'lake_share': 0.27, 'rain': 1600.0, 'side': 'near'})
    figures = dict(lakes='650', cloud='32%')
    names = {'sea:near-side sea': 'Uamoana'}
    assert RP.fill(entry, figures, names) == 'Some 650 lakes drain to Uamoana, under 32% cloud.'
    assert RP.fill(entry, figures, {'sea:near-side sea': 'Ua'}) == 'Some 650 lakes drain to Ua, under 32% cloud.'
    facts = dict(main_zone='wet', lake_share=0.30, rain=1700.0, side='near')
    assert RP.stale(entry, facts, names) == []
    moved = dict(main_zone='humid', lake_share=0.40, rain=2000.0, side='near')
    reasons = RP.stale(entry, moved, names)
    assert [r.split(':')[0] for r in reasons] == ['main_zone', 'lake_share', 'rain']
    assert RP.stale(entry, facts, {}) == ['[Uamoana]: sea:near-side sea gone']
    with pytest.raises(KeyError):
        RP.fill(dict(entry, text='{nonsense}'), figures, names)
    with pytest.raises(KeyError):
        RP.fill(dict(entry, text='to [Elsewhere]'), figures, names)


@pytest.mark.skipif(not RP.OUT.exists() or not D.OUT.exists(), reason='the fact sheets have not been built')
def test_fact_sheets_cover_the_regions_and_the_texts_resolve():
    product = json.loads(RP.OUT.read_text())
    record = json.loads(D.OUT.read_text())
    land = [r for r in record['regions'] if r['kind'] == 'land']
    assert product['districts_product'] == record['product']['sha256'], 'the fact sheets predate the districts draw'
    assert [p['anchor'] for p in product['regions']] == [r['anchor'] for r in land]
    for p in product['regions']:
        for key in ('area', 'districts', 'air', 'rain', 'earth', 'centre'):
            assert p['figures'].get(key), (p['name'], key)
        assert p['gloss'], p['name']
        assert p['facts']['side'] in ('near', 'far', 'limb')
    texts = RP.load_texts()
    anchors = {p['anchor'] for p in product['regions']}
    assert set(texts) <= anchors, sorted(set(texts) - anchors)
    names = RP.names_by_anchor(record)
    by_anchor = {p['anchor']: p for p in product['regions']}
    for anchor, entry in texts.items():
        p = by_anchor[anchor]
        RP.fill(entry, p['figures'], names)
        assert RP.stale(entry, p['facts'], names) == [], (p['name'], RP.stale(entry, p['facts'], names))
