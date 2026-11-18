"""The draft regions and districts against their rules, and the names' composition in the Moon's tongues."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import district_names as N
from research.studies.atlas_round import districts as D


def test_every_tongue_names_every_head_and_modifier():
    for key, t in N.TONGUES.items():
        assert set(t['heads']) == set(N.HEADS), key
        assert set(t['mods']) == set(N.MODIFIERS), key
        assert key in N.WHERE


def test_each_tongue_joins_as_its_grammar_does():
    assert N.compose('terlunan', 'bay', 'salt') == 'Mitifanga'                  # modifier-head, fused
    assert N.compose('terlunan', 'lake', 'white') == "'Uo'uoloko"               # capital past the glottal stop
    assert N.compose('olivolingvo', 'river', 'cold') == "H'Naharo Keraha"      # h' + noun + adjective, as the poem
    assert N.compose('tagalhebrew', 'lake', 'light') == 'Lawang Or'            # -ng after a vowel
    assert N.compose('tagalhebrew', 'mountains', 'white') == 'Bundok na Lavan'  # na after a consonant
    assert N.compose('satoyama', 'river', 'black') == 'Kurogawa'               # sequential voicing
    assert N.compose('sahran', 'valley', 'light') == 'Wadi Nur'
    assert N.compose('dala', 'lake', 'black') == 'Karakol'


def test_landscapes_take_their_tongues():
    region = lambda **k: dict(zones=k.pop('zones'), at=dict(turned=[k.pop('lat', 0.0), 0.0]),
                              relief_m=k.pop('relief', 2000), shore_share=k.pop('shore', 0.2))
    assert N.tongue_for(region(zones={'arid': 0.8}, lat=70)) == 'olivolingvo'
    assert N.tongue_for(region(zones={'arid': 0.7, 'semi-arid': 0.2})) == 'sahran'
    assert N.tongue_for(region(zones={'semi-arid': 0.7}, relief=1500)) == 'dala'
    assert N.tongue_for(region(zones={'semi-arid': 0.7}, relief=3200)) == 'olivolingvo'
    assert N.tongue_for(region(zones={'wet': 0.9}, relief=3800)) == 'runa'
    assert N.tongue_for(region(zones={'wet': 0.9}, shore=0.7)) == 'terlunan'
    assert N.tongue_for(region(zones={'humid': 0.7, 'wet': 0.3})) == 'pwani'
    assert N.tongue_for(region(zones={'wet': 0.9}, relief=3000)) == 'satoyama'
    assert N.tongue_for(region(zones={'wet': 0.9}, relief=2000)) == 'tagalhebrew'


@pytest.fixture(scope='module')
def draft():
    if not D.OUT.exists() or not D.PRODUCT.exists():
        pytest.skip('run python -m research.studies.atlas_round.districts first')
    record = json.loads(D.OUT.read_text())
    if record['product']['sha256'] != D.sha256(D.PRODUCT):
        pytest.fail('the districts product no longer matches its record')
    with np.load(D.PRODUCT, allow_pickle=False) as z:
        return record, z['district'].astype(np.int64), z['region'].astype(np.int64)


def test_districts_nest_in_regions_and_cover_the_land(draft):
    record, district, region = draft
    land = (region > 0) & (region < 100)
    assert np.all(district[land] > 0)                                       # every land cell is in a district
    pairs = np.unique(np.stack([district[district > 0], region[district > 0]]), axis=1)
    assert np.unique(pairs[0]).size == pairs.shape[1]                       # each district lies in one region


def test_districts_keep_to_their_sizes(draft):
    record, _, _ = draft
    sizes = np.array([d['area_km2'] for d in record['districts']])
    assert np.median(sizes) > 10_000
    assert (sizes > 1.05 * D.MAX_KM2).sum() == 0
    assert (sizes < D.MIN_KM2).sum() <= 10                                  # small islands only
    assert 30 <= len(record['regions']) <= 60


def test_every_unit_has_a_unique_name_and_every_coast_its_land(draft):
    record, _, _ = draft
    units = record['regions'] + record['districts'] + record['coastal_seas']
    names = [u['name'] for u in units]
    assert len(set(names)) == len(names)
    for u in record['districts'] + record['coastal_seas']:
        if u['name_source'] == 'generated':
            assert N.compose(u['tongue'], u['name_head'], u['name_mod']) in u['name'], u['name']   # built from its tongue
    assert all(c['faces'] for c in record['coastal_seas'])


def test_names_keep_to_their_sources(draft):
    record, _, _ = draft
    table = N.load_table()
    units = record['regions'] + record['districts'] + record['coastal_seas']
    assert {u['name_source'] for u in units} <= {'generated', 'calque', 'pinned'}
    for u in units:
        if u['anchor'] in table:                                            # a pin wins over the rules
            assert u['name'] == table[u['anchor']]['name'] and u['name_source'] == 'pinned'
        elif u['name_source'] == 'calque':                                  # the old name's meaning, on its own head
            gloss, key = N.CALQUES[u['feature']['name']]
            assert u['name_meaning'].endswith(f"after {u['feature']['name']}, the {gloss}")
    carried = [u['feature']['name'] for u in record['districts'] + record['coastal_seas'] if u.get('feature')]
    assert len(carried) == len(set(carried))                                # each feature names one unit


def test_a_calque_takes_its_tongues_word_on_the_units_head():
    u = dict(tongue='pwani', name_head='coast', feature=dict(name='Mare Orientale'))
    name, meaning = N.calque(u, set())
    assert name == 'Pwani Mashariki' and meaning == 'east coast, after Mare Orientale, the Eastern Sea'
    assert N.calque(dict(u, feature=dict(name='Mare Moscoviense')), set()) is None   # named after a place: no meaning
    assert N.calque(u, {'Pwani Mashariki'}) is None                                  # a name already taken


def test_anchors_survive_a_redraw():
    district = dict(kind='land', feature=None, centre_name='Plato', centre='crater', at=dict(today=[51.6, -9.4]))
    assert N.anchor_of(district) == 'crater:Plato'
    assert N.anchor_of(dict(district, feature=dict(name='Montes Jura'))) == 'feature:Montes Jura'
    assert N.anchor_of(dict(kind='sea', body='near-side sea')) == 'sea:near-side sea'
    units = dict(districts=[dict(region='L01', area_km2=5, anchor='crater:A'), dict(region='L01', area_km2=9, anchor='crater:B')])
    assert N.anchor_of(dict(kind='land', code='L01', districts=2, at=dict(today=[0, 0])), units) == 'region:crater:B'


def test_islands_take_the_territory_whose_waters_they_lie_in(monkeypatch):
    # an equatorial strip at 8 px/deg, cells of about 14 km2: a mainland, an island 30 km off it, a second island 38 km
    # beyond the first (more than 60 km from the mainland), two islands far out whose waters touch, and a lone rock
    monkeypatch.setattr(D, 'REGION_MIN_KM2', 5000.0)
    rows, cols = 64, 240
    lat = 4.0 - (np.arange(rows) + 0.5) / 8.0
    lon = (np.arange(cols) + 0.5) / 8.0
    area = (D.MOON_KM * np.radians(1 / 8.0)) ** 2 * np.cos(np.radians(lat))
    land = np.zeros((rows, cols), bool)
    land[10:50, 10:40] = True                       # the mainland, about 17,000 km2
    land[28:31, 48:51] = True                       # 8 cells (30 km) off its shore
    land[28:31, 61:64] = True                       # 10 cells beyond the first island
    land[10:26, 140:156] = True                     # two islands of about 3,700 km2, out of reach,
    land[30:46, 166:182] = True                     # 10 cells apart
    land[40, 220] = True                            # a rock alone
    territory, comps, joined, through, groups, rocks = D.island_territories(land, area, lat, lon, log=lambda *a: None)
    main = territory[30, 20]
    assert territory[29, 49] == main and territory[29, 62] == main            # the chain joins with its first link
    assert joined.index(comps[29, 49]) < joined.index(comps[29, 62])
    assert len(groups) == 1 and territory[18, 148] == territory[38, 174] == groups[0]
    assert rocks[40, 220] and territory[40, 220] == -1 and not rocks[30, 20]


def test_a_rename_keeps_every_name(draft):
    import copy
    record, _, _ = draft
    units = copy.deepcopy(dict(regions=record['regions'], districts=record['districts'],
                               coastal_seas=record['coastal_seas']))
    N.apply(units, previous=N.previous_names(record), log=lambda *a: None)
    for key in ('regions', 'districts', 'coastal_seas'):
        assert [u['name'] for u in units[key]] == [u['name'] for u in record[key]], key
