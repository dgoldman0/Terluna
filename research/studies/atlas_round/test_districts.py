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
        assert N.compose(u['tongue'], u['name_head'], u['name_mod']) in u['name'], u['name']   # built from its tongue
    assert all(c['faces'] for c in record['coastal_seas'])
