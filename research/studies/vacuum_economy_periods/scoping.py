"""The measures of the scoping study, from the period model: the base case, the account in dollars alone and the
programme's ramp year by year; what the funders put in and get back, with the internal rate of return of their flows;
and how far each uncertain input moves the residents the settlement carries by 2341.

usage (repository root):
  python3 -m research.studies.vacuum_economy_periods.scoping          results/scoping.json

The base case is central demand and central technology with the residents a result, under the default policy of U
(the liberal policy with U's rate floored at 75% of par and the funders' credit line behind the reserve). The
sensitivities read the saved product (results/periods.json) where it holds the case and run the model where it does
not; each moves one input, all else as in the base case. The technology groups move one group of the stack at a time
to its slow and its fast path (technology.GROUPS), the rest central.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

from . import economy, run, technology

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PRODUCT = HERE / 'results/periods.json'
OUT = HERE / 'results/scoping.json'
SCHEMA = 'terluna.research.vacuum-economy-periods.scoping/1'
YEARS = range(0, run.HORIZON + 1)

# Groups of payments to Earth and of Earth's holdings of U, for reading the year-by-year series.
PAYMENTS = {'partners_and_services': ('lines',), 'residents_purchases': ('pay',), 'supplies': ('supplies',),
            'equipment_plant_power_tugs': ('equipment', 'plant', 'power', 'fleet'),
            'fares_teleoperation_and_other': ('fares', 'teleoperation', 'zone_goods', 'dollars_held')}
HOLDINGS = {'working_balances': ('u_balances',), 'bonds': ('u_bonds',), 'investors_stakes': ('u_staked',),
            'holders_for_use': ('u_holders_U', 'u_holders_staked'), 'funders': ('u_funders_U', 'u_funders_staked')}

# One input moved at a time from the base case: the case in the saved product, or the options to run.
SENSITIVITIES = (
    ('Technology', 'slow', 'central:slow:steady', 'fast', 'central:fast:steady'),
    ("Earth's demand", 'low', 'low:central:steady', 'high', 'high:central:steady'),
    ("Share of residents' pay spent on Earth", '20%', (('share_spent_on_Earth', 0.2),), '5%', (('share_spent_on_Earth', 0.05),)),
    ("Floor on U's rate", 'par', 'central:central:steady:U_liberal:floor_1', '50%', 'central:central:steady:U_liberal:floor_0.5'),
    ('Money', 'dollars alone', 'central:central:steady:dollar', None, None),
    ('Policy of U', 'conservative', 'central:central:steady:U_conservative:floor_0.75', None, None),
    ("Earth's purchases when paid in U", '10% fewer', 'central:central:steady:U_liberal:trade_effect_-0.1', None, None),
    ("Funders' returns paid in U", 'none', 'central:central:steady:U_liberal:funders_paid_in_U_0', None, None),
    ("Holders' spending of their U in the zone", 'none', 'central:central:steady:U_liberal:spending_share_0', '60%',
     'central:central:steady:U_liberal:spending_share_0.6'),
    ('Months of spending the holders for use keep in U', 'one', 'central:central:steady:U_liberal:holders_months_1', 'six',
     'central:central:steady:U_liberal:holders_months_6'),
    ('Bonds', 'one-year', 'central:central:steady:U_liberal:bond_term_1', None, None),
    ("Investors' premium for a new currency", '1 point', 'central:central:steady:U_liberal:investor_premium_0.01', None, None),
    ('Metals in the reserve', 'none', 'central:central:steady:U_liberal:metal_share_0', None, None),
    ("The funders' backstop", 'none', 'central:central:steady:U_liberal:no_backstop', None, None),
    ("Water set aside on the sources", 'none', 'central:central:steady:whole_pool', None, None),
    ('Computing contracted at the opening', '1 GW', (('computing_W_2049', 1e9),), '20 GW', (('computing_W_2049', 2e10),)),
    ("Visits' response to their price", '0.5', (('visits_price_elasticity', 0.5),), '2', (('visits_price_elasticity', 2.0),)),
)

# Each group of the technology stack moved alone to its slow and its fast path, the rest of the stack central.
GROUP_LABELS = {'launch': 'Launch from Earth and the Moon', 'tugs_and_sources': 'Tugs and water from the sources', 'power': 'Power',
                'industry': 'Local industry, equipment and output per job', 'food_and_life_support': 'Food and life support',
                'computing': 'Computation per watt'}


def irr(flows):
    """The rate at which the flows' present value is zero: searched between 0 and 100% a year when the flows sum to a
    gain and between -20% and 0 when they sum to a loss, since deficits near the horizon would otherwise make the present
    value change sign again at rates far below zero; None when it does not change sign over the range searched."""
    npv = lambda r: sum(f / (1 + r) ** t for t, f in enumerate(flows))
    lo, hi = (0.0, 1.0) if npv(0.0) >= 0 else (-0.2, 0.0)
    if npv(lo) * npv(hi) > 0:
        return None
    for _ in range(200):
        mid = (lo + hi) / 2
        if npv(lo) * npv(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def series(rows):
    """The year-by-year measures the study reads: residents, the account, the funders' flows, U and its reserve, the
    payments to Earth and Earth's holdings of U by group, the receipts by line and the beamed stations."""
    group = lambda r, keys, src: sum(src.get(k, 0.0) for k in keys)
    out = dict(
        residents=[r['residents'] for r in rows], receipts=[r['receipts'] for r in rows], paid_to_Earth=[r['paid_to_Earth'] for r in rows],
        net=[r['net'] for r in rows], funders_flow=[r.get('funders_flow', r['net']) for r in rows],
        workforce=[r['workforce'] for r in rows], counted_jobs=[sum(r['jobs'].values()) for r in rows],
        firm_GW=[r['firm_GW'] for r in rows],
        power_by_option_GW={o: [r['power_by_option_kW'].get(o, 0.0) / 1e6 for r in rows] for o in sorted({o for r in rows for o in r['power_by_option_kW']})},
        payments={g: [group(r, keys, r['to_Earth']) for r in rows] for g, keys in PAYMENTS.items()},
        receipts_by_line={k: [r['lines'][k]['receipts'] if k in r['lines'] and r['lines'][k]['served'] else 0.0 for r in rows]
                          for k in ('computing', 'visits', 'culture', 'institutions', 'propellant', 'materials', 'metals')})
    cum, total = [], 0.0
    for f in out['funders_flow']:
        total += f
        cum.append(total)
    out['funders_cumulative'] = cum
    if 'u_rate' in rows[0]:
        out.update(u_rate=[r['u_rate'] for r in rows], u_reserve=[r['u_reserve'] for r in rows],
                   u_reserve_minimum=[r['u_reserve_minimum'] for r in rows], u_run_need=[r['u_run_need'] for r in rows],
                   u_reserve_metals=[r['u_reserve_metals'] for r in rows], u_earth_pay=[r['u_earth_pay'] for r in rows],
                   u_output=[r['u_output'] for r in rows],
                   holdings_USD={g: [group(r, keys, r) * r['u_rate'] for r in rows] for g, keys in HOLDINGS.items()})
    return out


def measures(rows):
    s = economy.summary(rows)
    flows = [r.get('funders_flow', r['net']) for r in rows]
    year = lambda y: None if y is None else 2041 + y
    return dict(plan_reached=year(s['plan_reached_year']), repaid=year(s['repaid_year']), stays_repaid_from=year(s['payback_year']),
                residents_2241=rows[200]['residents'], residents_2341=rows[run.HORIZON]['residents'],
                funding_needed_USD2041=s['funding_needed_USD2041'], opening_funding_USD2041=s['opening_funding_USD2041'],
                NPV3_to_2241_USD2041=s['NPV_USD2041']['0.03']['200'], NPV3_to_2341_USD2041=s['NPV_USD2041']['0.03'][str(run.HORIZON)],
                undiscounted_to_2341_USD2041=s['NPV_USD2041']['0.0'][str(run.HORIZON)],
                funders_IRR_to_2241=irr(flows[:201]), funders_IRR_to_2341=irr(flows),
                u_rate_lowest=s.get('u_rate_lowest'), u_years_at_floor=s.get('u_years_at_floor'),
                funders_lent_after_repayment_USD2041=-sum(f for r, f in zip(rows, flows) if r.get('u_repaid') and f < 0),
                u_test_failure_years=s['u_test_failure_years'], exhausted=year(s['exhausted_year']))


def product():
    saved = json.loads(PRODUCT.read_text())
    if saved['producer']['sha256'] != digest(HERE / 'run.py') or saved['dependencies'] != run.dependencies():
        raise RuntimeError('results/periods.json is stale; run `python3 -m research.studies.vacuum_economy_periods.run --write` first')
    cases = {c['name']: c for c in saved['cases']}
    paths = {}
    for name, scenario, tcase, ramp, opts in (('base', 'central', 'central', 'steady', ()),
                                              ('dollar', 'central', 'central', 'steady', (('u_policy', 'dollar'),)),
                                              ('ramp', 'central', 'central', 'constant_growth', ()),
                                              ('ramp_dollar', 'central', 'central', 'constant_growth', (('u_policy', 'dollar'),))):
        rows, rates = run.case_rows(scenario, tcase, ramp, opts)
        paths[name] = dict(scenario=scenario, technology=tcase, ramp=ramp, options=dict(opts), growth_rates=rates,
                           measures=measures(rows), series=series(rows))
    # What the funders earn depends on how fast the settlements grow, since growth is paid for from the account with
    # Earth: the programme's ramp held at the plan, the base case's growth to the plan held there, and the base case.
    to_plan = paths['base']['growth_rates']['to_plan']
    returns = {}
    for name, label, kw in (('programme_pace_held', "the programme's ramp of 3.3% a year to the plan, then held at it",
                             dict(ramp='constant_growth', past_plan='hold')),
                            ('fastest_to_plan_held', "the base case's growth to the plan, then held at it",
                             dict(schedule=run.steady_schedule(to_plan, 0.0))),
                            ('base', 'the base case, growing past the plan', None)):
        m = paths['base']['measures'] if kw is None else measures(economy.Run('central', 'central', **kw).simulate())
        returns[name] = dict(label=label, **{k: m[k] for k in ('plan_reached', 'repaid', 'residents_2241', 'residents_2341',
                                                                'funding_needed_USD2041', 'funders_IRR_to_2241', 'funders_IRR_to_2341')})
    base = paths['base']['measures']['residents_2341']
    swings = []
    for label, low_name, low, high_name, high in SENSITIVITIES:
        ends = []
        for name, case in ((low_name, low), (high_name, high)):
            if case is None:
                continue
            if isinstance(case, str):
                residents = cases[case]['summary']['residents'][str(run.HORIZON)]
                source = case
            else:
                rates, sched = run.solve('central', 'central', 'steady', tuple(sorted(case)))
                residents = sched[run.HORIZON]
                source = dict(case)
            ends.append(dict(setting=name, residents_2341=residents, change=residents - base, source=source))
        swings.append(dict(input=label, ends=ends))
    swings.sort(key=lambda s: -max(abs(e['change']) for e in s['ends']))
    groups = []
    for group in technology.GROUPS:
        ends = []
        for path in ('slow', 'fast'):
            tcase = 'central|%s=%s' % (group, path)
            rates, sched = run.solve('central', tcase, 'steady', ())
            m = measures(economy.Run('central', tcase, schedule=sched).simulate())
            ends.append(dict(setting=path, residents_2341=sched[run.HORIZON], change=sched[run.HORIZON] - base, technology=tcase,
                             growth_rates=rates, plan_reached=m['plan_reached'], repaid=m['repaid'],
                             funding_needed_USD2041=m['funding_needed_USD2041'], NPV3_to_2341_USD2041=m['NPV3_to_2341_USD2041'],
                             funders_IRR_to_2341=m['funders_IRR_to_2341']))
        groups.append(dict(group=group, label=GROUP_LABELS[group], members=list(technology.GROUPS[group]), ends=ends))
    groups.sort(key=lambda s: -max(abs(e['change']) for e in s['ends']))
    grid = {}
    for policy in ('', ':dollar'):
        for scenario in run.SCENARIOS:
            for tcase in ('static', 'slow', 'central', 'fast'):
                c = cases['%s:%s:steady%s' % (scenario, tcase, policy)]
                s = c['summary']
                grid['%s:%s%s' % (scenario, tcase, policy or ':U')] = dict(
                    residents_2241=s['residents'].get('200'), residents_2341=s['residents'].get(str(run.HORIZON)),
                    plan_reached=None if s['plan_reached_year'] is None else 2041 + s['plan_reached_year'],
                    funding_needed_USD2041=s['funding_needed_USD2041'], NPV3_to_2341_USD2041=s['NPV_USD2041']['0.03'][str(run.HORIZON)])
    return dict(
        schema=SCHEMA, producer=dict(path=str((HERE / 'scoping.py').relative_to(ROOT)), sha256=digest(HERE / 'scoping.py')),
        dependencies=dict(run.dependencies(), **{str(PRODUCT.relative_to(ROOT)): digest(PRODUCT)}),
        evidence=('Measures read from the period model and its saved product on the model\'s stated inputs; the internal rate of '
                  'return is that of the funders\' yearly flows, in constant dollars of 2041, so a real rate.'),
        units='Residents; constant dollars of 2041; years of the calendar; gigawatts; rates as fractions a year.',
        reading_rule=('paths holds the base case (central demand and technology, residents a result, the default policy of U), the '
                      'same in dollars alone, and the programme\'s ramp under each, year by year from 2041 (index 0) to 2341. '
                      'sensitivities moves one input at a time from the base case and gives the residents by 2341 and the change '
                      'from the base case, largest first. technology_groups moves each group of the technology stack alone to '
                      'its slow and its fast path, the rest central, with the same measure, largest first, and beside it the year '
                      'the plan is reached, the funding needed and the funders\' present value and rate of return. grid gives '
                      'residents as a result for each demand and technology, under the default policy of U (:U) and in dollars '
                      'alone (:dollar). returns_by_path gives the funders\' measures for three growth paths on central demand '
                      'and technology: the programme\'s ramp held at the plan, the base case\'s growth to the plan held there, '
                      'and the base case growing past it.'),
        paths=paths, returns_by_path=returns, sensitivities=swings, technology_groups=groups, grid=grid)


def main():
    OUT.write_text(json.dumps(product(), indent=1) + '\n')
    print(OUT)


if __name__ == '__main__':
    main()
