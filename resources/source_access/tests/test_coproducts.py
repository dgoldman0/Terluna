"""Coproduct units, conserved elements, evidence boundaries and input identities."""
import copy
import json

import pytest

from resources.source_access import coproducts as m
from resources.source_access import protection


@pytest.fixture(scope='module')
def inputs():
    return m.load()


def test_one_ppm_means_one_tonne_per_million_tonnes_and_conserves_element():
    r = {'retention': .7, 'separation': .8, 'finishing': .95}
    out = m.element_balance(1e9, 1, .1, r)
    assert out['contained_kg'] == 1000
    assert out['refined_kg'] == pytest.approx(53.2)
    assert sum(v for k, v in out.items() if k != 'contained_kg') == pytest.approx(1000)
    assert out['untreated_kg'] == 900
    assert out['separation_residue_kg'] == pytest.approx(14)


def test_absent_plant_and_stage_failure_deliver_zero_without_erasing_inventory():
    r = {'retention': .7, 'separation': .8, 'finishing': .95}
    for share in (0, 1):
        for stage in r:
            result = m.element_balance(10, 100, share, {**r, stage: 0})
            assert result['refined_kg'] == 0
            assert result['contained_kg'] == pytest.approx(.001)
            assert sum(v for k, v in result.items() if k != 'contained_kg') == pytest.approx(.001)


def test_invalid_recovery_grade_or_mass_is_rejected():
    r = {'retention': 1, 'separation': 1, 'finishing': 1}
    for feed, grade, share, recovery in [(-1, 1, 1, r), (1, 1e6+1, 1, r),
        (1, 1, 1.1, r), (1, 1, 1, {**r, 'retention': float('nan')}),
        (1, 1, 1, {'recovery': .5}), (float('inf'), 1, 1, r)]:
        with pytest.raises(ValueError):
            m.element_balance(feed, grade, share, recovery)


def test_transcription_units_groups_and_missing_measurements(inputs):
    _, assays, _ = inputs
    p, g = assays['profiles'], assays['group_definitions']
    assert len(p['ci']['elements']) == 40
    assert len(p['kreep15386']['elements']) == 40
    assert m.grade(p['ci'], 'PGM6', g) == pytest.approx(3.225)
    assert m.grade(p['kreep15386'], 'REY', g) == pytest.approx(902.45)
    assert m.grade(p['kreep15386'], 'NdPr', g) == pytest.approx(155.6)
    assert set(g['NdPr']) < set(g['REY'])
    assert not set(g['PGM6']) & set(g['REY'])
    assert m.grade(p['kreep15386'], 'PGM6', g) is None
    assert m.grade(p['hoba'], 'Pt', g) == 30.78  # ng/g to ppm, not 30,780 ppm
    assert m.grade(p['filomena'], 'Au', g) == .610
    assert p['kreep15386']['elements']['Ga']['reported_sigma_ppm'] is None


def test_transfer_and_recovery_are_separate_common_scenarios(inputs):
    c, assays, _ = inputs
    p, g = assays['profiles']['ci'], assays['group_definitions']
    sums = m.envelope(1e9, p, 'PGM6', .01, c, g, [.5, 1, 2])
    assert sums[0] == pytest.approx(1e9*.01*3.225e-6*.5*.12)
    assert sums[1] == pytest.approx(1e9*.01*3.225e-6*.532)
    assert sums[2] == pytest.approx(1e9*.01*3.225e-6*2*.893475)
    members = [m.envelope(1e9, p, e, .01, c, g, [.5, 1, 2]) for e in g['PGM6']]
    assert [sum(a[i] for a in members) for i in range(3)] == pytest.approx(sums)


def test_selected_feed_and_recovery_scale_monotonically(inputs):
    c, assays, _ = inputs
    p, g = assays['profiles']['ci'], assays['group_definitions']
    a = m.envelope(1e9, p, 'Ni', .01, c, g, [.5, 1, 2])
    b = m.envelope(1e9, p, 'Ni', .1, c, g, [.5, 1, 2])
    assert b == pytest.approx([v*10 for v in a])
    assert sorted(a) == a
    changed = copy.deepcopy(c)
    changed['recovery_scenarios']['middle']['finishing'] = 0
    assert m.envelope(1e9, p, 'Ni', .01, changed, g, [.5, 1, 2])[1] == 0


def test_marginal_break_even_uses_joint_feed_once():
    r = {'retention': .7, 'separation': .8, 'finishing': .95}
    # A whole marginal circuit costing $0.01/kg feed, supported by PGM alone.
    assert m.break_even(3.225, r, .01) == pytest.approx(.01/(3.225e-6*.532))
    assert m.break_even(0, r, .01) is None


def test_saved_product_reproduces_exactly_and_pins_all_consumed_inputs():
    result = m.read(m.HERE/'results/coproducts.json.gz')
    assert result == m.run()
    assert result['producer']['sha256'] == m.sha(m.__file__)
    for name, expected in result['dependencies'].items():
        assert m.sha(m.HERE/name) == expected
    assert result['additional_saleable_output_kg'] is None
    assert result['additional_export_revenue_USD2041'] is None
    assert result['extraction_authorized'] is False


def test_annual_products_use_fresh_feed_once_and_keep_case_feasibility(inputs):
    c, _, flows = inputs
    result = m.read(m.HERE/'results/coproducts.json.gz')
    lookup = {(r['name'], r['year']): r for r in result['annual_cases']}
    assert len(lookup) == len(result['annual_cases']) == 243
    for case in flows['cases']:
        for row in case['checkpoints']:
            out = lookup[case['name'], row['year']]
            assert out['physical_feasible'] == row['physical_feasible']
            assert out['admitted_completed_path'] == case['admitted_completed_path']
            expected = row['source_feed_kg_year']*.01*3.225e-6*.532
            assert out['ci']['refined_kg_year_low_middle_high']['PGM6'][1] == pytest.approx(expected)
    reference = lookup[c['reference_case'], 200]
    assert reference['source_feed_kg_year'] == pytest.approx(2681426291609.669)
    assert reference['ci']['refined_kg_year_low_middle_high']['PGM6'][1]/1000 == pytest.approx(46.005230885147085)


def test_patch_uses_square_kilometres_and_the_original_density_product(inputs):
    c, _, _ = inputs
    result = m.read(m.HERE/'results/coproducts.json.gz')
    lunar = m.read(m.HERE/'results/lunar_potential.json.gz')
    row = next(r for r in lunar['bulk_inventory'] if r['depth_m'] == 5)
    out = result['hypothetical_kreep_patch']
    assert out['area_km2'] == 100
    assert out['bulk_mass_kg_p05_p50_p95'] == pytest.approx([x*100e6 for x in row['column_kg_m2_p05_p50_p95']])
    assert out['contained_REY_kg_low_middle_high'][1] == pytest.approx(out['bulk_mass_kg_p05_p50_p95'][1]*902.45e-6)
    assert out['refined_REY_kg_low_middle_high'][1] == pytest.approx(out['contained_REY_kg_low_middle_high'][1]*.532)


def test_named_metal_candidates_remain_outside_current_orbits_and_unqualified():
    result = m.read(m.HERE/'results/coproducts.json.gz')
    assert len(result['named_candidates']) == 2
    for target in result['named_candidates']:
        assert set(target['failed_wider_orbit_axes']) == {'a', 'e'}
        assert target['in_wider_orbital_screen'] is False
        assert target['release_status'] == 'unresolved'
        assert target['qualified_stock_kg'] is None
        assert target['precious_metal_assay_available'] is False
    for name in ('Bennu', 'Ryugu', 'Itokawa', 'Eros'):
        release = protection.evaluate_release({'name': name, 'region': 'near_earth'})
        assert release['release_status'] == 'blocked'
        assert release['extraction_authorized'] is False
