"""The ledger's array form against a plain cohort-by-cohort account of the same rules: equal to the last bit.

The plain account below follows each component and each of its cohorts one at a time, in the order the array form
adds them, so the two agree exactly, not merely within rounding."""
from copy import deepcopy
import json
import math
from pathlib import Path
import pytest
from provisioning.vacuum_industry import model as m, recovery

HERE = Path(__file__).resolve().parents[1]


def plain_ledger(data, rows, targets, seed, price_index=None):
    """The ledger kept component by component and cohort by cohort."""
    cats = m.catalogue(data)
    start = rows[0]['year']
    end = rows[-1]['year']
    lead = 2
    flows = {y: m.blank_flow() for y in range(start - lead, end + 1)}
    stocks = {y: {} for y in range(start, end + 1)}
    events = []
    working = m.local_working_parts(data)
    shape = data.get('lifetime_shape', 3.0)
    for name, a in cats.items():
        for part, share, life in [('structure', a['structural_share'], a['structural_life']), ('active', 1 - a['structural_share'], a['life'])]:
            if share <= 0:
                continue
            phases = max(1, min(8, math.ceil(life)))
            denominator = sum(m.survival(j + 1, life, shape) for j in range(phases))
            if denominator <= 0:
                raise ValueError('Seed deployment ages exceed every component lifetime')
            cohorts = [(start - j - 1, seed.get(name, 0.0) * share / denominator) for j in range(phases)]
            curve = [m.survival(age, life, shape) for age in range(end - start + phases + 2)]
            for row in rows:
                y = row['year']
                expired = alive = 0.0
                surviving = []
                for birth, mass in cohorts:
                    age = y - birth
                    fraction = curve[age]
                    if y > start:
                        expired += mass * (curve[max(0, age - 1)] - fraction)
                    if fraction > 1e-12:
                        surviving.append((birth, mass))
                        alive += mass * fraction
                cohorts = surviving
                required = targets[y].get(name, 0.0) * share
                added = max(0.0, required - alive)
                if added > 1e-06:
                    cohorts.append((y, added))
                    order = y - lead
                    replacement = min(added, expired)
                    growth = added - replacement
                    able = 0.0 if part != 'active' or name in working['excluded_assets'] else m.at(working['share'], order)
                    made = able * m.at(working['staffed'], order)
                    labour = m.at(working['labour_fraction_of_import_price'], order)
                    if not (0 <= made <= 1 and labour >= 0):
                        raise ValueError('The locally made share of working parts lies between 0 and 1, and its labour is not negative')
                    learned = price_index.get(order, 1.0) if price_index and part == 'active' else 1.0
                    for workyear in (order, order + 1):
                        f = flows[workyear]
                        amount = added / 2
                        growth_work = growth / 2
                        replace_work = replacement / 2
                        unit = a['price'] * learned if part == 'active' or a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
                        if part == 'active':
                            bought = 1 - made
                            f['imports'][a['destination']] += amount * bought
                            f['external_capital'] += growth_work * bought * unit
                            f['external_replacement'] += replace_work * bought * unit
                            f['capital_imports'][a['destination']] += growth_work * bought
                            if a['owner'] == 'industry':
                                f['industry_replacement_import_kg'] += replace_work * bought
                            if a['owner'] in ('service', 'hospitality'):
                                f[a['owner'] + '_external_capital'] += growth_work * bought * unit
                                f[a['owner'] + '_external_replacement'] += replace_work * bought * unit
                            if name == 'return_capsule_buses':
                                f['capsule_purchases_kg'] += amount
                            f['local'][a['destination']] += amount * made
                            f['capital_local'][a['destination']] += growth_work * made
                            f['local_working_capital'] += growth_work * made * unit * labour
                            f['local_working_replacement'] += replace_work * made * unit * labour
                            f['local_working_labour_if_staffed'] += amount * able * unit * labour
                        elif a['material'] == 'bulk':
                            f['bulk'] += amount
                            f['capital_bulk'] += growth_work
                        else:
                            f['local'][a['destination']] += amount
                            f['capital_local'][a['destination']] += growth_work
                        if part == 'structure':
                            f['local_capital'] += growth_work * unit
                            f['local_replacement'] += replace_work * unit
                        if a['owner'] == 'service':
                            f['service_capital'] += growth_work * unit
                            f['service_replacement'] += replace_work * unit
                        if a['owner'] == 'hospitality':
                            f['hospitality_capital'] += growth_work * unit
                            f['hospitality_replacement'] += replace_work * unit
                        if name in ('commercial_civil_bulk', 'service_civil_bulk'):
                            f['new_floor'] += growth_work / data['economy']['exports']['factory_bulk_kg_m2']
                    events.append(dict(asset=name, component=part, order_year=order, commission_year=y, mass_kg=added, replacement_kg=replacement,
                                       locally_made_share=made, local_labour_fraction=labour, price_index=learned))
                if expired:
                    f = flows[y]
                    f['retired_mass'] += expired
                    if part == 'structure':
                        if a['material'] == 'bulk':
                            f['bulk_recycled'] += expired * data['industrial']['lunar']['recovery']
                        else:
                            f['recycled'][a['destination']] += expired * (data['industrial']['lunar']['recovery'] if a['destination'] == 'surface' else data['economy']['production']['source_recovery'])
                available = alive + added
                old = stocks[y].get(name, math.inf)
                stocks[y][name] = min(old, available / share)
    previous_beds = rows[0].get('visitor_beds', 0.0)
    for row in rows:
        beds = row.get('visitor_beds', 0.0)
        added = max(0.0, beds - previous_beds)
        previous_beds = beds
        for y in (row['year'] - 2, row['year'] - 1):
            flows[y]['new_visitor_beds'] += added / 2
    return flows, stocks, events


def plain_retirements(data, rows, seed, events):
    """The mass of each component that expires in each year, cohort by cohort."""
    cats = m.catalogue(data)
    start = rows[0]['year']
    end = rows[-1]['year']
    shape = data.get('lifetime_shape', 3.)
    out = {y: {} for y in range(start, end + 1)}
    births = {}
    for n, a in cats.items():
        for part, share, life in [('active', 1 - a['structural_share'], a['life']), ('structure', a['structural_share'], a['structural_life'])]:
            if not share:
                continue
            phases = max(1, min(8, math.ceil(life)))
            den = sum(m.survival(j + 1, life, shape) for j in range(phases))
            births[n, part] = [(start - j - 1, seed.get(n, 0.) * share / den) for j in range(phases)]
    for ev in events:
        births[ev['asset'], ev['component']].append((ev['commission_year'], ev['mass_kg']))
    for (n, part), cohorts in births.items():
        life = cats[n]['life' if part == 'active' else 'structural_life']
        curve = [m.survival(a, life, shape) for a in range(end + 2)]
        values = [0.] * (end + 1)
        for birth, mass in cohorts:
            for y in range(max(start + 1, birth + 1), end + 1):
                age = y - birth
                values[y] += mass * (curve[age - 1] - curve[age])
        for y in range(start, end + 1):
            out[y][n, part] = values[y]
    return out


def same(a, b, path=''):
    """Whether two records are equal value for value, type for type and in the same order; raises at the first difference."""
    assert type(a) is type(b), path
    if isinstance(a, dict):
        assert list(a) == list(b), path
        for k in a:
            same(a[k], b[k], '%s/%s' % (path, k))
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            same(x, y, '%s[%d]' % (path, i))
    else:
        assert a == b, '%s: %r and %r' % (path, a, b)
    return True


@pytest.fixture(scope='module')
def data():
    """The model's coefficients with the recovery measures' equipment, and working parts made in part from year 20."""
    plain = json.loads((HERE / 'engineering.json').read_text())
    spec = json.loads((HERE / 'recovery.json').read_text())
    configured = recovery.configure(plain, spec, True, True)
    configured['economy']['production']['local_working_parts'] = dict(
        share={'by_year': [[20, 0.0], [60, 0.6]]}, staffed={'by_year': [[8, 1.0], [40, 0.5], [70, 0.8]]},
        labour_fraction_of_import_price=0.25, excluded_assets=['source_computer'])
    return configured


def growing(data, last=110):
    """Rows and targets that grow, pause, shrink for a while and grow again: idle surplus, retirements unreplaced and
    cohorts that wear out within the rows."""
    rows = [dict(year=y, population=12500. * 1.03 ** (y - 8), visitor_beds=40. + (y - 8)) for y in range(8, last + 1)]
    names = list(m.catalogue(data))
    targets = {}
    for row in rows:
        y = row['year']
        level = 1.03 ** (y - 8) if y < 50 else 1.03 ** 42 * (0.98 ** (y - 50) if y < 70 else 0.98 ** 20 * 1.02 ** (y - 70))
        targets[y] = {n: 1e4 * (1 + i % 7) * level for i, n in enumerate(names) if (i % 5 or y >= 30)}
    seed = {n: 1e4 * (1 + i % 7) * 0.9 for i, n in enumerate(names) if i % 3}
    return rows, targets, seed


@pytest.mark.parametrize('shape', [3.0, None])
@pytest.mark.parametrize('learning', [False, True])
def test_the_array_ledger_keeps_the_cohort_by_cohort_account_exactly(data, shape, learning):
    data = dict(data, lifetime_shape=shape)
    rows, targets, seed = growing(data)
    index = {o: 0.97 ** max(0, o - 8) for o in range(6, 110)} if learning else None
    plain = plain_ledger(data, rows, targets, seed, index)
    array = m.cohort_ledger(data, rows, targets, seed, price_index=index)
    assert len(plain[2]) > 1000 and any(e['replacement_kg'] < e['mass_kg'] for e in plain[2])
    assert same(plain, array)
    assert same(plain_retirements(data, rows, seed, plain[2]), recovery.retirements(data, rows, seed, array[2]))


def test_the_orders_worked_in_a_year_are_those_placed_in_it_and_the_year_before_in_the_ledgers_order():
    events = [dict(asset=a, order_year=y) for a, y in [('a', 6), ('a', 7), ('a', 9), ('b', 7), ('b', 8)]]
    worked = m.worked_by_year(events)
    assert worked[7] == [events[0], events[1], events[3]]
    assert worked[8] == [events[1], events[3], events[4]]
    assert worked[10] == [events[2]] and 11 not in worked


def test_the_array_ledger_raises_where_the_cohort_by_cohort_account_does(data):
    rows, targets, seed = growing(data, 40)
    wrong = deepcopy(data)
    wrong['economy']['production']['local_working_parts']['share'] = 1.4
    with pytest.raises(ValueError, match='locally made share'):
        plain_ledger(wrong, rows, targets, seed)
    with pytest.raises(ValueError, match='locally made share'):
        m.cohort_ledger(wrong, rows, targets, seed)
    with pytest.raises(ValueError, match='Consecutive'):
        m.cohort_ledger(data, rows[:3] + rows[4:], targets, seed)
