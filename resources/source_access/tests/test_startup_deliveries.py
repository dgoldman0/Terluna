"""Eight elapsed years, once-only seed inventory, and conservative arrival timing."""
import math

import pytest

from resources.source_access import startup_deliveries as m


RECOVERY = {'retention': 1., 'separation': 1., 'finishing': 1.}


def test_fractional_delivery_delay_and_conservation_at_year_eight():
    result = m.deliveries([1.]*8, [1.]*8, 1e6, RECOVERY, .8, 1.5, .01)
    assert result['refined_total_kg'] == 8
    assert result['received_on_Earth_by_year_8_kg'] == pytest.approx(5.148)
    assert result['awaiting_or_in_transit_at_year_8_kg'] == pytest.approx(1.2)
    assert result['retained_off_Earth_kg'] == pytest.approx(1.6)
    assert result['lost_on_Earth_return_by_year_8_kg'] == pytest.approx(.052)
    assert sum(result[k] for k in ['retained_off_Earth_kg', 'awaiting_or_in_transit_at_year_8_kg', 'lost_on_Earth_return_by_year_8_kg', 'received_on_Earth_by_year_8_kg']) == pytest.approx(8)
    assert sum(result['annual_Earth_receipts_kg']) == pytest.approx(5.148)
    assert result['annual_Earth_receipts_kg'][:2] == pytest.approx([0, .396])


def test_absent_refinery_no_export_and_arrival_after_cutoff():
    assert m.deliveries([1.]*8, [0.]*8, 1e6, RECOVERY, .8, 1.5, .01)['received_on_Earth_by_year_8_kg'] == 0
    assert m.deliveries([1.]*8, [1.]*8, 1e6, RECOVERY, 0, 1.5, .01)['retained_off_Earth_kg'] == 8
    late = m.deliveries([1.]*8, [1.]*8, 1e6, RECOVERY, .8, 8, .01)
    assert late['received_on_Earth_by_year_8_kg'] == 0
    assert late['awaiting_or_in_transit_at_year_8_kg'] == pytest.approx(6.4)
    with pytest.raises(ValueError):
        m.deliveries([1.]*9, [1.]*9, 1e6, RECOVERY, .8, 1.5, .01)
    with pytest.raises(ValueError):
        m.deliveries([1.]*8, [1.]*8, 1e6, RECOVERY, .8, -1, .01)


def test_feed_allocation_counts_late_procurement_once():
    record = {'base_kg': 36., 'late_years_6_7_kg': 4., 'total_kg': 40.}
    assert m.allocate_feed(record, list(range(1, 9))) == [1, 2, 3, 4, 5, 6, 9, 10]
    assert sum(m.allocate_feed(record, [1]*8)) == 40
    with pytest.raises(ValueError):
        m.allocate_feed({**record, 'total_kg': 44}, [1]*8)
    with pytest.raises(ValueError):
        m.allocate_feed(record, [0]*8)


def test_unchanged_opening_and_exact_seed_totals():
    config, feeds, _, _ = m.load()
    assert len(feeds['cases']) == 81
    assert {r['opening_population'] for r in feeds['cases']} == {12500}
    case = next(r for r in feeds['cases'] if r['name'] == config['reference_case'])
    assert case['source_feed']['total_kg'] == 298938653315.0383
    assert case['lunar_raw']['total_kg'] == 7702686198.0891485
    for row in feeds['cases']:
        for key in ('source_feed', 'lunar_raw'):
            assert math.fsum(m.allocate_feed(row[key], config['feed_weights'])) == pytest.approx(row[key]['total_kg'], rel=1e-13)


def test_capacity_scaling_does_not_change_the_lunar_source_selection():
    result = m.c.read(m.HERE/'results/startup_deliveries.json.gz')
    ref = result['reference_by_final_source_share']
    small = ref['0.01']['source']['PGM6']['middle']['received_on_Earth_by_year_8_kg']
    medium = ref['0.1']['source']['PGM6']['middle']['received_on_Earth_by_year_8_kg']
    large = ref['1']['source']['PGM6']['middle']['received_on_Earth_by_year_8_kg']
    assert medium == pytest.approx(10*small)
    assert large == pytest.approx(10*medium)
    assert ref['0.01']['lunar'] == ref['0.1']['lunar'] == ref['1']['lunar']
    assert medium/1000 == pytest.approx(10.718485377576856)
    assert ref['0.1']['source']['PGM6']['middle']['annual_Earth_receipts_kg'][:3] == [0, 0, 0]


def test_group_elements_and_all_arrival_balances_close():
    result = m.c.read(m.HERE/'results/startup_deliveries.json.gz')
    for r in result['reference_by_final_source_share'].values():
        for scenario in m.c.SCENARIOS:
            group = r['source']['PGM6'][scenario]['received_on_Earth_by_year_8_kg']
            assert group == pytest.approx(sum(r['source'][e][scenario]['received_on_Earth_by_year_8_kg'] for e in ('Pt', 'Pd', 'Rh', 'Ru', 'Ir', 'Os')))
        for origin in ('source', 'lunar'):
            for metric in r[origin].values():
                for b in metric.values():
                    assert b['refined_total_kg'] == pytest.approx(sum(b[k] for k in ('retained_off_Earth_kg', 'awaiting_or_in_transit_at_year_8_kg', 'lost_on_Earth_return_by_year_8_kg', 'received_on_Earth_by_year_8_kg')))


def test_saved_result_and_dependencies_reproduce():
    result = m.c.read(m.HERE/'results/startup_deliveries.json.gz')
    assert result == m.run()
    assert len(result['case_summaries']) == 243
    assert result['producer']['sha256'] == m.c.sha(m.__file__)
    for name, digest in result['dependencies'].items():
        assert m.c.sha(m.HERE/name) == digest
    assert result['incremental_revenue_USD2041'] is None
