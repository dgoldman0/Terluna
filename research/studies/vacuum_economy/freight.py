"""The tugs to the near-Earth sources, chosen by what the industry costs with them, period by period, in each scenario.

usage: python3 -m research.studies.vacuum_economy.freight [--scenario NAME ...] [--jobs N]      (in the repository root)

The tugs mine the nearest water of the source survey first, so the leg they fly lengthens as the water near the depot
is drawn (provisioning/vacuum_industry, with resources/source_access/reach.py). For each scenario of demand the
constant-growth case's sales are decided once with the founding design's chemical tug, and the industry that fills
those orders is then sized with each kind of tug throughout:

- `chemical`, the founding design's tug on hydrogen and oxygen from the sources' water, its tanks following the leg;
- `electric_water`, an electric tug on the sources' water at the performance a laboratory thruster has measured;
- `electric_xenon`, an electric tug on xenon bought on Earth, at the performance of a thruster built for flight;
- `electric_water_target`, an electric tug on water at the 12 November comparison's development target.

In each reporting period the kind chosen is the one whose industry costs least over the period's years, capital and
operating together, undiscounted, at the prices of each year with learning; the founding deployment counts in the
first period. A kind stands in a period only where it flies every year of the period and, on xenon, needs no more
than the world's yearly output. A kind that cannot fly some year is sized again flying until the period before, the
cheapest kind that flies throughout taking over from then, and stands only in the periods before. The scenario's
freight, `chosen`, is chosen among the first three; `chosen_with_target` shows the choice with the development target
among them. Sales decided with the chosen tugs may differ from those the comparison holds fixed.
"""
from __future__ import annotations
import argparse
import json
import math
import multiprocessing
from pathlib import Path

from provisioning.vacuum_industry import industry, model
from . import inputs, run
from .inputs import HERE, ROOT, digest

SCHEMA = 'terluna.research.vacuum-economy-freight/1'
OUTPUT = HERE / 'results/freight.json'
CANDIDATES = ('chemical', 'electric_water', 'electric_xenon')
TARGET = 'electric_water_target'
SCENARIOS = ('low', 'central', 'high')
# The industry's costs by category, the founding deployment apart, as the industry model gives them each year.
CATEGORIES = ('capital_external_USD2041', 'capital_local_USD2041', 'capital_resource_provision_USD2041', 'operating_external_replacement_USD2041',
              'operating_external_logistics_materials_USD2041', 'operating_external_services_USD2041', 'operating_local_replacement_USD2041',
              'operating_local_payroll_USD2041', 'operating_resource_provision_USD2041')


def dependencies():
    return inputs.dependencies((HERE / 'freight.py', HERE / 'run.py'))


def schedule(kinds, periods):
    """A schedule by year from one kind for each period, merging periods that keep the same kind."""
    points = []
    for (first, _), kind in zip(periods, kinds):
        if not points or points[-1][1] != kind:
            points.append([first, kind])
    return points[0][1] if len(points) == 1 else {'by_year': points}


def size(case, book, kind):
    """The industry for a book of orders with the given tugs: (its years, the year it cannot fly, if any)."""
    s = inputs.sources()
    data = inputs.engineering(dict(case, freight='chemical'))
    data['industrial']['transport']['source_propulsion'] = kind
    try:
        physical, _ = industry.build(data, s['recovery'], book, 'combined')
    except model.Unreachable as exc:
        return None, exc.year, str(exc)
    return physical, None, None


def account(physical, periods, supply):
    """Each period's industry cost, whether every year of it is feasible, and the largest yearly xenon from Earth."""
    out = []
    for first, last in periods:
        rows = [r for r in physical['annual'] if first <= r['year'] <= last]
        cost = math.fsum(r['cash']['capital_total_USD2041'] + r['cash']['operating_total_USD2041'] for r in rows)
        if first <= physical['annual'][0]['year'] - 1:
            cost += physical['seed']['cost']['total_USD2041']
        xenon = max((r['physical']['source_tugs']['xenon_from_Earth_kg_year'] for r in rows), default=0.0)
        legs = [r['physical']['source_tugs']['leg_m_s'] for r in rows]
        parts = {k: math.fsum(r['cash'][k] for r in rows) for k in CATEGORIES}
        out.append(dict(first_year=first, last_year=last, industry_cost_USD2041=cost, by_category_USD2041=parts,
                        source_feed_kg=math.fsum(r['physical']['source_feed_kg_year'] for r in rows), feasible=all(r['feasible'] for r in rows),
                        xenon_from_Earth_largest_kg_year=xenon, within_world_xenon_supply=xenon <= supply,
                        leg_m_s=[min(legs), max(legs)] if legs else None))
    return out


def choose(scenario):
    s = inputs.sources()
    periods = [tuple(p) for p in s['periods']]
    case = dict(name='%s:constant_growth' % scenario, scenario=scenario, ramp='constant_growth', recovery='combined')
    # The orders are decided once, with the founding freight, which flies every leg; each kind of tug then fills them.
    print(scenario, 'deciding the sales with the founding freight', flush=True)
    book = run.decide(dict(case, freight='founding'))['book']
    supply = inputs.sources()['engineering']['transport']['source_propulsion']['options']['electric_xenon']['world_supply_kg_year']
    kinds = (*CANDIDATES, TARGET)
    results = {}
    for kind in kinds:
        print(scenario, 'sizing with', kind, flush=True)
        physical, failed, message = size(case, book, kind)
        results[kind] = dict(physical=physical, failed_year=failed, failure=message)
    whole = sorted((k for k in kinds if results[k]['physical'] is not None),
                   key=lambda k: math.fsum(p['industry_cost_USD2041'] for p in account(results[k]['physical'], periods, supply)))
    accounts = {k: account(results[k]['physical'], periods, supply) for k in whole}
    # A kind that cannot fly some year stands only in the periods before that year's. Its costs there come from sizing it
    # until then and handing over to the cheapest kind that flies throughout; where no kind flies throughout, from the
    # orders up to the end of the period before, whose last years buy no equipment for the growth after.
    for kind in kinds:
        failed = results[kind]['failed_year']
        if failed is None:
            continue
        cut = next(i for i, (first, last) in enumerate(periods) if first <= failed <= last)
        if cut == 0:
            continue
        if whole:
            successor = whole[0]
            print(scenario, 'sizing with', kind, 'until year', periods[cut][0] - 1, 'then', successor, flush=True)
            physical, again, message = size(case, book, schedule([kind] * cut + [successor] * (len(periods) - cut), periods))
            results[kind]['handover'] = dict(to=successor, from_year=periods[cut][0])
        else:
            # Orders that end earlier draw less, so the cut moves back a period at a time until the kind flies them.
            physical = None
            while physical is None and cut > 1:
                last = periods[cut - 1][1]
                print(scenario, 'sizing with', kind, 'for the orders to year', last, flush=True)
                physical, again, message = size(case, dict(book, annual=[r for r in book['annual'] if r['year'] <= last]), kind)
                results[kind]['handover'] = dict(to=None, orders_to_year=last)
                if physical is None:
                    cut = min(cut - 1, next(i for i, (first, end) in enumerate(periods) if first <= again <= end))
        if physical is None:
            results[kind]['failure'] += '; sized again: ' + message
            continue
        unknown = [dict(first_year=f, last_year=l, industry_cost_USD2041=None, by_category_USD2041=None, source_feed_kg=None, feasible=False,
                        xenon_from_Earth_largest_kg_year=None, within_world_xenon_supply=None, leg_m_s=None) for f, l in periods[cut:]]
        accounts[kind] = [dict(p, available=i < cut) for i, p in enumerate(account(physical, periods[:cut], supply) + unknown)]

    def pick(options):
        """The cheapest kind standing in each period; where none stands, the kind before it, recorded as unreachable."""
        chosen, unreachable = [], []
        for i in range(len(periods)):
            usable = [k for k in options if k in accounts and accounts[k][i].get('available', True) and accounts[k][i]['feasible']
                      and accounts[k][i]['within_world_xenon_supply']]
            if not usable:
                unreachable.append(dict(first_year=periods[i][0], last_year=periods[i][1]))
                chosen.append(chosen[-1] if chosen else options[0])
                continue
            chosen.append(min(usable, key=lambda k: accounts[k][i]['industry_cost_USD2041']))
        return dict(by_period=[dict(first_year=p[0], last_year=p[1], kind=k) for p, k in zip(periods, chosen)], schedule=schedule(chosen, periods),
                    periods_no_kind_reaches=unreachable)

    return scenario, dict(
        comparison={k: dict(periods=accounts.get(k), failed_year=results[k]['failed_year'], failure=results[k]['failure'],
                            handover=results[k].get('handover'),
                            total_industry_cost_USD2041=math.fsum(p['industry_cost_USD2041'] for p in accounts[k]) if k in whole else None)
                    for k in kinds},
        chosen=pick(CANDIDATES), chosen_with_target=pick((*CANDIDATES, TARGET)))


def product(scenarios, jobs=1):
    if jobs > 1:
        with multiprocessing.Pool(jobs) as pool:
            found = dict(pool.map(choose, scenarios, chunksize=1))
    else:
        found = dict(choose(s) for s in scenarios)
    return dict(schema=SCHEMA, producer=dict(path=str(Path(__file__).resolve().relative_to(ROOT)), sha256=digest(__file__)),
                dependencies=dependencies(), candidates=list(CANDIDATES), development_target=TARGET,
                evidence='A comparison of conditional sizings: the same orders filled with each kind of tug, the tugs mining the '
                         'nearest water of the source survey first. Costs follow the stated coefficients; no tug of these kinds '
                         'has been built for bulk cargo, and the velocity changes are figures of merit, not trajectories.',
                reading_rule='For each scenario, the constant-growth orders as decided with chemical tugs, sized with each kind of '
                             'tug. Each period gives the industry\'s capital and operating costs over its years, undiscounted, the '
                             'founding deployment in the first, whether every year was feasible, the largest yearly xenon bought '
                             'on Earth against the world\'s output, and the range of the leg flown. A kind that could not fly some '
                             'year was sized again with the cheapest kind that flies throughout from the period of that year, and '
                             'stands only before it. The schedule chosen takes the cheapest kind standing in each period.',
                scenarios=found)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--scenario', action='append', choices=SCENARIOS, help='a scenario of demand; repeat for several (default: all)')
    parser.add_argument('--jobs', type=int, default=1, help='scenarios chosen at once (default: one)')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = product(args.scenario or SCENARIOS, args.jobs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1, sort_keys=True, allow_nan=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
