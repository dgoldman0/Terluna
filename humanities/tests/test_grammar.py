"""The vocabulary is well formed and the grammar places units with physics-derived durations, links them where their
states disagree, and applies the choreographic devices as stated."""
import hashlib
import json
import math

import pytest

from humanities.dance import grammar as g
from shared.constants import MOON_SURFACE_GRAVITY as GM

PRODUCT = json.loads(g.OUT.read_text())


def test_product_regenerates_and_binds_its_code_and_vocabulary():
    assert PRODUCT == json.loads(json.dumps(g.results()))
    assert PRODUCT['schema'] == g.SCHEMA
    root = g.HERE.parents[1]
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest()[:16] == h, path


def test_every_unit_is_well_formed_and_uses_known_shapes():
    regions = g.POSES['regions']
    for uid, u in g.UNITS.items():
        assert u['family'] in ('travel', 'level', 'turn', 'flight', 'floor', 'upper', 'fabric', 'partner'), uid
        assert u['lunar'].split(':')[0] in ('adapted', 'native', 'native in degree'), uid
        assert set(u['effort']) == {'weight', 'time', 'space', 'flow'}, uid
        if not u.get('realisable', True):
            continue
        p = g.params_for(uid, None)
        for ph in g.expand_phases(u, p, 'R'):
            assert ph['duration'] > 0, (uid, ph)
            shapes = ph.get('shapes', [])
            for s in shapes if isinstance(shapes, list) else []:
                for region in ('legs', 'arms', 'torso', 'head', 'feet'):
                    if region in s:
                        assert s[region] in regions[region], (uid, region, s[region])


def test_durations_follow_the_physics():
    jete = g.place('grand_jete', 'R', {'rise_m': 1.5, 'speed_m_s': 1.5}, g.State(support='L', speed=1.5))
    flight = [p for p in jete.phases if p['kind'] == 'flight'][0]
    assert flight['duration'] == pytest.approx(2 * math.sqrt(2 * 1.5 / GM), rel=1e-9)
    # The body can never sink faster than g: a 0.35 m descent at full gravity takes 2 sqrt(d/g).
    assert g.descend_s(0.35, 1.0) == pytest.approx(2 * math.sqrt(0.35 / GM))
    # A fall from 2 to 45 degrees takes the envelope's toppling time.
    env = json.loads((g.HERE / 'results' / 'envelope.json').read_text())
    assert g.topple_s(2, 45) == pytest.approx(env['balance']['from_2_to_45_deg_s']['moon'], abs=0.01)


def test_the_grammar_links_what_does_not_fit():
    ph = g.Phrase('test')
    ph.add('canopy', 'R')                       # sleeves go up and are gathered again inside the unit
    ph.add('sleeve_cast', 'R')                  # cast: sleeves aloft
    ph.add('canopy', 'R')                       # needs gathered sleeves: a gather is inserted, waiting on the fabric
    links = [p for p in ph.placed if p.connector]
    assert any(p.unit == 'sleeve_recover' for p in links)
    ph2 = g.Phrase('speed')
    ph2.add('grand_jete', 'R', speed_m_s=1.8)   # from standing still: strides build the approach speed first
    assert ph2.placed[0].connector and ph2.placed[0].unit == 'orbit_run'
    build = ph2.placed[0]
    assert build.exit.speed == pytest.approx(1.8)
    # Each stride's stance carries the weight and reverses the vertical speed: the floor's grip over all strides
    # covers the speed built.
    vz = g.takeoff(build.params['stride_rise_m'])
    per_stride = sum(g.MU * (GM * p['duration'] + 2 * vz) for p in build.phases if p['kind'] == 'stance')
    assert per_stride >= 1.8 - 1e-9


def test_devices():
    m = g.motif(('skim', 'R', {'hops': 4}), ('flight_phrase', 'R', {'rise_m': 1.0}))
    assert [s for _, s, _ in g.mirror(m)] == ['L', 'L']
    aug = g.augment(m, 1.5)
    assert aug[1][2]['rise_m'] == pytest.approx(1.5)
    assert [u for u, _, _ in g.accumulate(m)] == ['skim', 'skim', 'flight_phrase']
    assert [u for u, _, _ in g.retrograde(m)] == ['flight_phrase', 'skim']


def test_the_example_phrases_are_lunar_and_flowing():
    first = PRODUCT['phrases'][0]['evaluation']
    assert first['lunar_native_share'] > 0.5 and first['flights_over_1_s'] >= 1 and first['stops'] == 0
    assert 'air' in first['levels'] and 'falling' in first['levels']
    gen_a, _ = g.generate(5)
    gen_b, _ = g.generate(5)
    assert g.evaluate_phrase(gen_a) == g.evaluate_phrase(gen_b)          # reproducible from its seed
