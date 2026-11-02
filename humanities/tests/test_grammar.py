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


def test_weight_shifts_and_falls_follow_gravity():
    # A shift onto one foot is an inverted pendulum pushed from the feet's edges: its time scales as 1/sqrt(g).
    moon = g.shift_s(0.10)
    g.G = 9.80665
    try:
        earth = g.shift_s(0.10)
        grip_earth = g.fall_grip(35)
    finally:
        g.G = GM
    assert 1.5 < moon < 2.0 and moon / earth == pytest.approx(math.sqrt(9.80665 / GM), rel=1e-6)
    # Falling as one piece asks the same grip on both worlds, under a point mass's tan(angle), and the vocabulary
    # stops falls where the dance floor's grip still holds.
    assert g.fall_grip(35) == pytest.approx(grip_earth, rel=1e-9)
    assert g.fall_grip(35) < math.tan(math.radians(35)) and g.fall_grip(35) <= g.MU
    assert g.fall_grip(g.UNITS['long_fall']['params']['to_deg']['max']) <= g.MU + 0.01


def test_junctions_never_stand_one_foot_through_two_stances():
    for ph in (g.first_phrase(), g.suspension_study()):
        seq = [(pl, p) for pl in ph.placed for p in pl.phases]
        for (pa, a), (pb, b) in zip(seq, seq[1:]):
            if pa is not pb and a['kind'] == b['kind'] == 'stance':
                assert a['support'] != b['support'] or a['support'] not in ('L', 'R'), (ph.name, pa.unit, pb.unit)
    first = g.first_phrase()
    units = [pl.unit for pl in first.placed]
    # Standing on both feet, the weight shifts in place before the fall; the jete's landing is the carve.
    assert units[units.index('long_fall') - 1] == 'weight_shift'
    jete = first.placed[units.index('grand_jete')]
    assert jete.dropped_landing and jete.exit.support == 'none' and jete.phases[-1]['kind'] == 'flight'
    # The whirl starts from rest and turns on the spot; the spiral ends kneeling, still.
    whirl = first.placed[units.index('whirl')]
    assert whirl.entry.speed <= 0.15 and whirl.exit.speed == 0.0 and whirl.phases[0]['support'] == 'both'
    assert first.placed[-1].exit.support == 'floor' and first.placed[-1].phases[-1]['support'] == 'kneel'


def test_the_air_study_turns_in_the_air_and_brakes_on_the_floor():
    air = g.air_study()
    e = g.evaluate_phrase(air)
    assert e['flights_over_1_s'] >= 5 and e['share_aloft'] > 0.5 and e['climax_at'] > 0.5
    assert e['aerial_revolutions'] == {'cartwheel': 1.0, 'somersault': 2.0, 'twist': 1.0, 'turn': 1.0}
    # A closing step that stops the travel takes as long as the floor's grip needs.
    for pl in air.placed:
        if pl.connector and pl.note.startswith('close to both feet, stopping'):
            assert pl.duration >= pl.entry.speed / (g.MU * GM) - 1e-9
