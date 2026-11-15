"""Material recovery and refurbishment: conservation, timing and what the measures add to the industry's work."""
from copy import deepcopy
import json
from pathlib import Path
import pytest
from provisioning.vacuum_industry import industry, model, recovery

HERE = Path(__file__).resolve().parents[1]
OPERATING = ('operating_external_replacement_USD2041', 'operating_external_logistics_materials_USD2041',
             'operating_external_services_USD2041', 'operating_local_replacement_USD2041',
             'operating_local_payroll_USD2041', 'operating_resource_provision_USD2041')


@pytest.fixture(scope='module')
def data():
    return json.loads((HERE / 'engineering.json').read_text())


@pytest.fixture(scope='module')
def spec():
    return json.loads((HERE / 'recovery.json').read_text())


def small_book():
    """Thirty years of a small settlement that grows 3.3% a year and sells some goods, propellant and computing."""
    rows = []
    for year in range(8, 38):
        people = 12500. * 1.033**(year - 8)
        rows.append(dict(year=year, population=people, new_residential_people=people * 1.033**2 * .033,
                         goods=2e5 * (year - 7), propellant_sales=5e5 * (year - 7), compute=2e6,
                         service_jobs=150., visitor_beds=40., visitor_nights_year=40 * 365.25 * .7))
    return dict(name='small', annual=rows)


@pytest.fixture(scope='module')
def solved(data, spec):
    book = small_book()
    return {measures: industry.build(data, spec, book, measures) for measures in ('none', 'combined')}


def test_the_stated_streams_reproduce_the_model_makeup_water_and_goods_are_counted_once(data, spec):
    account = recovery.water_account(spec, spec['normalization_population'])
    assert account['operational_makeup_kg'] == pytest.approx(data['demand']['yearly_water'])
    for v in account['streams'].values():
        assert v['recovered_kg'] + v['makeup_kg'] == pytest.approx(v['gross_kg'])
        assert v['loss_kg'] == v['makeup_kg']
    closed = recovery.configure(data, spec, True, False)
    g = spec['goods']
    gross = g['gross_kg_person_year'] * spec['normalization_population']
    assert closed['demand']['goods_virgin'] + closed['settlement']['operating_import_kg_year']['everyday_carbon_feed'] == pytest.approx(gross * (1 - g['collection'] * g['quality']))


def repair(spec, retirement_year):
    cats = {'asteroid_tugs': dict(destination='source', owner='industry', price=100.)}
    q = deepcopy(spec)
    q['refurbishment']['eligible_assets'] = ['asteroid_tugs']
    rows = [dict(year=y) for y in range(8, 16)]
    flows = {y: model.blank_flow() for y in range(6, 16)}
    events = [dict(asset='asteroid_tugs', component='active', mass_kg=100., replacement_kg=100., order_year=10, commission_year=12)]
    for y in (10, 11):
        flows[y]['imports']['source'] = 50.
        flows[y]['external_replacement'] = 5000.
        flows[y]['industry_replacement_import_kg'] = 50.
    retired = {retirement_year: {('asteroid_tugs', 'active'): 100.}}
    return flows, events, recovery.refurbish(cats, q, rows, flows, events, retired)


def test_a_machine_retired_in_the_year_of_an_order_cannot_fill_it(spec):
    flows, events, audit = repair(spec, 10)
    assert events[0]['refurbished_mass_kg'] == 0
    assert audit[10]['source']['pending_inspection_kg'] == 80
    assert audit[11]['source']['qualified_inventory_kg'] == 48
    assert audit[12]['source']['refurbished_commissioned_kg'] == 0
    assert sum(a['source']['imported_mass_avoided_kg'] for a in audit.values()) == 0


def test_accepted_machines_repair_parts_and_rejects_close_in_mass_and_money(spec):
    flows, events, audit = repair(spec, 9)
    assert events[0]['refurbished_mass_kg'] == 48
    assert events[0]['new_imported_mass_kg'] == pytest.approx(61.6)
    assert audit[9]['source']['collected_kg'] == 80
    assert audit[10]['source']['quality_reject_kg'] == 32
    assert audit[10]['source']['discarded_old_parts_kg'] == pytest.approx(9.6)
    assert audit[12]['source']['refurbished_commissioned_kg'] == 48
    assert audit[11]['source']['repair_WIP_kg'] == pytest.approx(48)
    assert audit[12]['source']['repair_WIP_kg'] == pytest.approx(0)
    assert sum(flows[y]['external_replacement'] for y in (10, 11)) == pytest.approx(6160)
    assert sum(a['source']['import_purchase_avoided_USD2026'] for a in audit.values()) == pytest.approx(3840)
    # What is left uncollected, what inspection rejects and the old parts the repair displaces are three losses.
    lost = sum(a['source']['uncollected_kg'] + a['source']['quality_reject_kg'] + a['source']['discarded_old_parts_kg'] for a in audit.values())
    assert lost == pytest.approx(61.6)
    for y, a in audit.items():
        assert a['source']['mass_balance_residual_kg'] == pytest.approx(0, abs=1e-12)
        assert all(a[k]['withdrawn_for_repair_kg'] == 0 for k in ('surface', 'depot', 'leo'))
        assert all(flows[y]['imports'][k] == 0 for k in ('surface', 'depot', 'leo'))
    assert sum(sum(f['loop_consume_kg'].values()) for f in flows.values()) > 0


def test_with_nothing_accepted_no_import_and_no_purchase_is_removed(spec):
    rejecting = deepcopy(spec)
    rejecting['refurbishment']['quality'] = 0
    flows, events, audit = repair(rejecting, 9)
    assert events[0]['refurbished_mass_kg'] == 0
    assert sum(flows[y]['external_replacement'] for y in (10, 11)) == 10000
    assert sum(a['source']['quality_reject_kg'] for a in audit.values()) == 80


def test_without_the_residual_allowance_every_stated_stream_keeps_its_losses(data, spec):
    e = recovery.Recovery(spec, residual_water_allowance=False)
    a = recovery.water_account(e.spec, 250000, improved=True)
    b = recovery.water_account(spec, 250000, improved=True)
    assert a['streams'] == b['streams']
    assert a['operational_makeup_kg'] > 0
    assert b['operational_makeup_kg'] - a['operational_makeup_kg'] == pytest.approx(spec['landscape_agriculture_allowance_kg_year'])
    assert e.configure(data)['demand']['yearly_water'] == pytest.approx(a['operational_makeup_kg'])


def test_the_measures_add_their_equipment_to_a_copy_of_the_coefficients(data, spec):
    before = json.dumps(data, sort_keys=True)
    configured = recovery.Recovery(spec).configure(data)
    assert json.dumps(data, sort_keys=True) == before
    assert 'loop_water_equipment' not in model.catalogue(data)
    assert {'loop_water_equipment', 'loop_sort_surface', 'loop_repair_depot'} <= set(model.catalogue(configured))


def test_recovery_fractions_outside_the_unit_interval_are_refused(spec):
    wrong = deepcopy(spec)
    wrong['goods']['collection'] = 1.2
    with pytest.raises(ValueError):
        recovery.Recovery(wrong)
    late = deepcopy(spec)
    late['refurbishment']['inspection_delay_years'] = 0
    with pytest.raises(ValueError):
        recovery.Recovery(late)


def test_no_measures_is_the_plain_model(data, solved):
    physical, loops = solved['none']
    plain = industry.price_opening_route_fuel(deepcopy(industry.compact(small_book(), model.solve(data, small_book()['annual'], max_iterations=50))), data)
    assert physical == plain
    assert loops['variant'] == 'none'


@pytest.mark.parametrize('measures', ['none', 'combined'])
def test_every_year_closes_its_water_and_material_balances(solved, measures):
    physical, loops = solved[measures]
    assert [r['year'] for r in loops['annual']] == list(range(38))
    for r in loops['annual']:
        assert r['chemical_exhaust_recovery_kg'] == 0
        if r['phase'] != 'operation':
            continue
        w = r['water']
        assert w['new_stock_and_construction_kg'] >= 0 and abs(w['rounding_residual_kg']) <= .001
        for stream in w['streams'].values():
            assert stream['gross_kg'] == pytest.approx(stream['recovered_kg'] + stream['loss_kg'])
        for v in r['materials']['structural'].values():
            assert v['retired_kg'] == pytest.approx(v['recovered_kg'] + v['collection_loss_kg'] + v['quality_reject_kg'], abs=.001)
            assert v['opening_inventory_kg'] + v['recovered_kg'] == pytest.approx(v['used_kg'] + v['closing_inventory_kg'], abs=.001)
            assert v['mass_balance_residual_kg'] == pytest.approx(0, abs=.001)
        for v in r['materials']['refurbishment'].values():
            assert v['mass_balance_residual_kg'] == pytest.approx(0, abs=.001)
            assert all(v[field] >= -.001 for field in ('pending_inspection_kg', 'qualified_inventory_kg', 'repair_WIP_kg'))
            assert v['collected_kg'] + v['uncollected_kg'] == pytest.approx(v['retired_kg'])
            assert v['accepted_kg'] + v['quality_reject_kg'] == pytest.approx(v['inspection_input_kg'])


@pytest.mark.parametrize('measures', ['none', 'combined'])
def test_the_founding_years_and_the_feed_are_counted_once(data, solved, measures):
    physical, loops = solved[measures]
    source = (recovery.Recovery(json.loads((HERE / 'recovery.json').read_text()), *recovery.OPTIONS[measures]).configure(data))['industrial']['source']
    for r in loops['annual']:
        for key, grade, share in [('source_water_equivalent_kg', 'water_grade', 'water_recovery'), ('source_nitrogen_kg', 'nitrogen_grade', 'nitrogen_recovery'),
                                  ('source_virgin_metal_kg', 'metal_grade', 'metal_recovery')]:
            assert r[key] <= r['source_feed_kg'] * source[grade] * source[share] * (1 + 1e-12) + .001
    seed, audit, founding = physical['seed'], loops['seed_audit'], loops['annual'][:8]
    assert sum(r['source_feed_kg'] for r in founding) == pytest.approx(seed['source_feed_kg'])
    assert sum(r['lunar_raw_kg'] for r in founding) == pytest.approx(seed['lunar_raw_kg'])
    for new, old in [('source_water_equivalent_kg', 'source_water_kg_year'), ('source_nitrogen_kg', 'source_nitrogen_kg_year'), ('source_virgin_metal_kg', 'source_virgin_metal_kg_year')]:
        assert sum(r[new] for r in founding) == pytest.approx(8 * audit['physical_annual_average'][old] + audit['late_preproduction_physical_total'][old])
    assert loops['cumulative_source_feed_kg'] == pytest.approx(sum(r['source_feed_kg'] for r in loops['annual']))


def test_the_measures_work_is_priced_and_lowers_the_feed_and_the_imports(solved):
    plain, combined = solved['none'][0], solved['combined'][0]
    for r in combined['annual']:
        cash = r['cash']
        assert cash['operating_total_USD2041'] == pytest.approx(sum(cash[k] for k in OPERATING), abs=.01)
        assert abs(cash['accounting_residual_USD2041']) < .02
        if r['year'] > 20:
            assert r['direct_refurbishment_jobs'] > 0 and r['direct_refurbishment_payroll_USD2041'] > 0
            assert sum(r['physical']['loop_extra_process_kWh_year'].values()) > 0
            assert sum(r['physical']['loop_repair_consumables_kg_year'].values()) > 0
        if r['maximum_capacity_use'] > 1 + 1e-6:
            assert not r['feasible'] and r['actual_sales_capacity_fraction'] == 0
    assert solved['combined'][1]['cumulative_source_feed_kg'] < solved['none'][1]['cumulative_source_feed_kg']
    late = lambda case: sum(r['cash']['operating_external_replacement_USD2041'] for r in case['annual'] if r['year'] > 25)
    assert late(combined) < late(plain)
