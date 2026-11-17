"""The cases of the vacuum economy period by period, the growth of residents as a result, the rough time to the
exhaustion of the sites and the sources' water, and the saved product.

usage (repository root):
  python3 -m research.studies.vacuum_economy_periods.run                      a table of every case
  python3 -m research.studies.vacuum_economy_periods.run --case NAME          one case, year by year
  python3 -m research.studies.vacuum_economy_periods.run --exhaustion         the rough time to exhaustion
  python3 -m research.studies.vacuum_economy_periods.run --write              results/periods.json

A case is a scenario of demand (low, central, high), a technology case (static, slow, central, fast) and a schedule of
residents to the horizon, 2341: one of the programme's ramps (constant_growth, front_loaded), carried on past the planned
6.47 million in 2241 at its final rate, or residents as a result (steady: the highest single growth rate from the opening
to the plan, then the highest single rate past it; fastest: the coupled study's rule, the highest steady rate period by
period). Residents as a result grow from the opening's 12,500 in year 8 as fast as the settlement can while its funders
still earn the programme's 3% (the net present value of the account with Earth at 3% to the horizon is not negative) and,
once the account is repaid, no later year is in deficit with Earth, with the settlement's own work staffed and its
sources' water not exhausted. In every schedule the sites bound the residents: the caves' cities and the eventual set of
crater sites in towns. Every case reports the year it reaches the plan, which is a milestone, not an end.

The rough time to exhaustion carries each path's growth on, technology held at 2241, until the water open to mining runs
out or to 3541, and reports when the sites fill and when the water runs out: with the survey's whole conditional pool and
with the project's set-asides for future use, protection and uncertainty.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

from research.studies.vacuum_economy import programme as definition
from . import economy, engineering, technology
from .designs import inputs, v
from .economy import OPENING

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.vacuum-economy-periods/2'
SCENARIOS = ('low', 'central', 'high')
RAMPS = ('constant_growth', 'front_loaded', 'steady')
HORIZON = economy.horizon()
PERIODS = ((9, 12), (13, 18), (19, 25), (26, 37), (38, 50), (51, 65), (66, 80), (81, 100), (101, 125), (126, 150), (151, 175), (176, 200),
           (201, 225), (226, 250), (251, 275), (276, 300))
# Sensitivities around the central case, each changing one option.
SENSITIVITIES = {
    'pay_all': dict(pay_counted='all'),
    'no_births': dict(births=0.0),
    'migrants_pay_fares': dict(fares_paid_by='migrants'),
    'towns_only': dict(allocation='towns'),
    'twenty_caves': dict(caves=20),
    'reach_p05': dict(reach='p05'),
    'reach_p95': dict(reach='p95'),
    'pay_on_Earth_20pct': dict(share_spent_on_Earth=0.2),
    'held_at_plan': dict(past_plan='hold'),
    'towns_of_1000': dict(town_residents=1000),
    'whole_pool': dict(water_set_aside=0.0),
}
TOLERANCE = 1e-6


def planned():
    return definition.load()['population']['planned']


def schedule(rates, last=None):
    """Residents by year 0 to `last` (the horizon) from period rates: (first year, last year, rate) triples after the
    opening's 12,500, with no growth outside them. The sites bound the residents in the run."""
    p = definition.load()['population']
    last = HORIZON if last is None else last
    n, out = p['opening'], [0.0] * OPENING + [p['opening']]
    for y in range(OPENING + 1, last + 1):
        n = n * (1 + next((rate for a, b, rate in rates if a <= y <= b), 0.0))
        out.append(n)
    return out


def steady_schedule(to_plan, past_plan=0.0, last=None):
    """Residents by year 0 to `last`: growth at `to_plan` a year from the opening until the planned 6.47 million is
    reached, then at `past_plan` a year."""
    p = definition.load()['population']
    last = HORIZON if last is None else last
    n, out, reached = p['opening'], [0.0] * OPENING + [p['opening']], False
    for y in range(OPENING + 1, last + 1):
        if reached:
            n = n * (1 + past_plan)
        else:
            n = min(n * (1 + to_plan), p['planned'])
            reached = n >= p['planned'] * (1 - 1e-12)
        out.append(n)
    return out


def sustained(scenario, tcase, sched, options):
    """Whether a schedule is carried: the funders earn the programme's 3% (net present value at 3% to the horizon not
    negative); the account with Earth is repaid and no ten years after that are together in deficit, so a year's
    reinvestment, such as a change of power supply that the factories cannot yet make in full, is covered within the
    decade and a lasting deficit is not; the settlements' own work is staffed; the sources' water is not exhausted.
    Under a policy of U the funders' flows carry Earth's purchases of U and the reserve, and every one of the five tests
    of a safe and healthy policy holds in every year. Returns the test and the rows."""
    rows = economy.Run(scenario, tcase, schedule=sched, **dict(options)).simulate()
    s = economy.summary(rows)
    ok = (s['NPV_USD2041']['0.03'][str(s['last_year'])] >= 0 and s['repaid_year'] is not None and s['deficit_decades_after_repayment'] == 0
          and s['hands_short_years'] == 0 and s['exhausted_year'] is None and s['u_test_failure_years'] == 0)
    return ok, rows


def highest(test):
    """The highest rate a year, to a hundredth of a percentage point, for which `test(rate)` holds, searched by bisection
    between 0 and 30%; None if even no growth fails."""
    lo, hi = 0.0, 0.3
    if not test(lo):
        return None
    for _ in range(14):
        mid = (lo + hi) / 2
        if test(mid):
            lo = mid
        else:
            hi = mid
    return lo


def solve(scenario, tcase, rule='steady', options=()):
    """Residents as a result. 'steady': the highest single rate from the opening to the plan, with residents held at the
    plan meanwhile, then the highest single rate past it. 'fastest': period by period, the highest steady rate that is
    sustained to the horizon with no growth after the period. Returns the rates and the schedule."""
    ok = lambda sched: sustained(scenario, tcase, sched, options)[0]
    if rule == 'steady':
        to_plan = highest(lambda r: ok(steady_schedule(r, 0.0)))
        if to_plan is None:
            return dict(rule='steady', to_plan=0.0, past_plan=0.0, sustained=False), steady_schedule(0.0)
        reached = next((y for y, n in enumerate(steady_schedule(to_plan)) if n >= planned() * (1 - 1e-9)), None)
        past = highest(lambda r: ok(steady_schedule(to_plan, r))) if reached is not None else None
        rates = dict(rule='steady', to_plan=to_plan, past_plan=past or 0.0, plan_year=reached, sustained=True)
        return rates, steady_schedule(to_plan, past or 0.0)
    rates = []
    for a, b in PERIODS:
        rate = highest(lambda r: ok(schedule(rates + [(a, b, r)])))
        rates.append((a, b, rate or 0.0))
    return rates, schedule(rates)


def case_rows(scenario, tcase, ramp, options=()):
    if ramp in ('steady', 'fastest'):
        rates, sched = solve(scenario, tcase, ramp, options)
        return economy.Run(scenario, tcase, schedule=sched, **dict(options)).simulate(), rates
    return economy.Run(scenario, tcase, ramp, **dict(options)).simulate(), None


# The rough time to exhaustion: central demand under each improving technology, on the programme's pace carried on and
# on growth as a result; with the survey's whole pool at its median, and with the set-asides at the median and at the
# 5th percentile; towns of 750, and of 1,000 on the programme's pace.
READINGS = (('whole_pool_p50', dict(reach='p50', water_set_aside=0.0)),
            ('set_aside_p50', dict(reach='p50', water_set_aside=0.9)),
            ('set_aside_p05', dict(reach='p05', water_set_aside=0.9)))


def exhaustion_horizon():
    return int(v(inputs()['growth']['exhaustion_horizon_year']))


def exhaustion(scenario, tcase, path, options=()):
    """One path carried on to the exhaustion horizon: the year the plan is reached, the year from which the account is in
    deficit for good (every ten years from it to the end together in deficit), the year the sites are full, the year the
    water open to mining runs out or lies beyond every tug (or None by 3541), the residents and the water drawn then. `path` is 'ramp' (the
    programme's constant growth, carried on at its final rate whatever the account) or 'result' (steady growth as a
    result to the horizon, where its own test puts the account at its limit, and held at that size after it)."""
    last = exhaustion_horizon()
    opts = dict(options)
    if path == 'ramp':
        run = economy.Run(scenario, tcase, 'constant_growth', last_year=last, stop_at_exhaustion=True, **opts)
        rates = None
    else:
        rates, sched = solve(scenario, tcase, 'steady', tuple(sorted(opts.items())))
        run = economy.Run(scenario, tcase, schedule=sched + [sched[-1]] * (last - HORIZON), last_year=last, stop_at_exhaustion=True, **opts)
    rows = run.simulate()
    s = economy.summary(rows)
    end = rows[-1]
    full = s['sites_full_year']
    nets = [r['net'] for r in rows]
    # The account is in deficit for good from the first year after which every ten-year window to the end of the run is
    # in deficit; a decade in deficit that the account recovers from, as at a change of power supply, does not count.
    lasting = None
    if s['repaid_year'] is not None:
        for y in range(len(rows) - economy.DECADE, s['repaid_year'], -1):
            if sum(nets[y:y + economy.DECADE]) >= 0:
                break
            lasting = y
    because = None
    pool = dict(economy.defaults(), **opts)
    if s['exhausted_year'] is not None:
        drawn_out = engineering.leg_for(end['contained_drawn_t'] * 1000.0, pool['reach'], pool['water_set_aside']) is None
        because = 'all drawn' if drawn_out else 'the rest beyond every tug'
    return dict(scenario=scenario, technology=tcase, path=path, options=opts, growth_rates=rates,
                usable_water_kg=engineering.usable_water_kg(pool['reach'], pool['water_set_aside']),
                capacity_residents=run.capacity(), plan_reached_year=s['plan_reached_year'], lasting_deficit_from_year=lasting,
                sites_full_year=full, residents_when_sites_full=None if full is None else rows[full]['residents'],
                water_exhausted_year=s['exhausted_year'], water_out_because=because, residents_at_end=end['residents'], last_year=end['year'],
                contained_drawn_kg_at_end=end['contained_drawn_t'] * 1000.0, leg_m_s_at_end=end['leg_m_s'])


def exhaustion_cases():
    out = []
    for tcase in ('slow', 'central', 'fast'):
        for path in ('ramp', 'result'):
            for name, opts in READINGS:
                out.append(('central:%s:%s:%s' % (tcase, path, name), 'central', tcase, path, tuple(sorted(opts.items()))))
        out.append(('central:%s:ramp:whole_pool_p50:towns_of_1000' % tcase, 'central', tcase, 'ramp',
                    tuple(sorted(dict(READINGS[0][1], town_residents=1000).items()))))
    return out


SERIES = ('residents', 'receipts', 'paid_to_Earth', 'net', 'cumulative', 'firm_GW', 'power_chosen', 'local_share', 'imported_t',
          'launches_year', 'pads', 'tug', 'leg_m_s', 'water_delivered_t', 'contained_drawn_t', 'visitors', 'visitor_share_of_world',
          'workforce', 'teleoperated_jobs', 'lines_served_share', 'housed_town', 'housed_city', 'housed_beyond_initial_quota', 'housing_short',
          'fare', 'landed_USD_kg', 'kWh_USD')
AT_YEARS = (8, 25, 50, 100, 150, 200, 250, 300)
# The ledger of U, in cases with a policy of U.
U_SERIES = ('trade_multiple', 'u_settlement_share', 'u_real_yield', 'u_output', 'u_balances', 'u_investment', 'u_investment_target',
            'u_funders_U', 'u_holdings', 'u_reserve', 'u_reserve_floor', 'u_reserve_own', 'u_from_balances', 'u_from_investors',
            'u_funders_converted', 'u_yield_income', 'u_spent_in_zone', 'u_paid_to_funders', 'funders_flow', 'funders_cumulative', 'u_growth_20')


def record(name, scenario, tcase, ramp, options, rows, rates):
    s = economy.summary(rows)
    keys = SERIES + (U_SERIES if 'funders_flow' in rows[0] else ())
    at = {}
    for h in AT_YEARS:
        r = rows[h]
        at[str(h)] = dict({k: r[k] for k in keys}, to_Earth=r['to_Earth'], jobs=r['jobs'], power_by_option_kW=r['power_by_option_kW'],
                          lines={k: dict(receipts=x['receipts'], served=x['served']) for k, x in r['lines'].items()},
                          production_kg=r['production_kg'], **({'u_tests': r['u_tests']} if 'u_tests' in r else {}))
    series = {k: [rows[y][k] for y in range(0, HORIZON + 1, 5)] for k in keys}
    return dict(name=name, scenario=scenario, technology=tcase, ramp=ramp, options=dict(options), growth_rates=rates,
                summary=s, at_year=at, every_fifth_year=series)


def cases():
    out = []
    for scenario in SCENARIOS:
        for tcase in technology.CASES:
            for ramp in RAMPS:
                out.append(('%s:%s:%s' % (scenario, tcase, ramp), scenario, tcase, ramp, ()))
    out.append(('central:central:fastest', 'central', 'central', 'fastest', ()))
    for key, opts in SENSITIVITIES.items():
        out.append(('central:central:constant_growth:%s' % key, 'central', 'central', 'constant_growth', tuple(sorted(opts.items()))))
    out.append(('central:central:steady:whole_pool', 'central', 'central', 'steady', (('water_set_aside', 0.0),)))
    # U beside the account with Earth: the two policies on the programme's ramp and with growth as a result, the liberal
    # one under each technology and on low demand, and its least certain holders' behaviour moved one at a time.
    for policy in ('conservative', 'liberal'):
        for ramp in ('constant_growth', 'steady'):
            out.append(('central:central:%s:U_%s' % (ramp, policy), 'central', 'central', ramp, (('u_policy', policy),)))
    for scenario, tcase in (('central', 'slow'), ('central', 'fast'), ('low', 'central'), ('low', 'fast')):
        out.append(('%s:%s:steady:U_liberal' % (scenario, tcase), scenario, tcase, 'steady', (('u_policy', 'liberal'),)))
    for key, value in (('u_trade_effect', 0.0), ('u_trade_effect', 0.3), ('u_spending_share', 0.0), ('u_spending_share', 0.6)):
        out.append(('central:central:steady:U_liberal:%s_%g' % (key[2:], value), 'central', 'central', 'steady',
                    tuple(sorted((('u_policy', 'liberal'), (key, value))))))
    return out


# The most liberal policy that stays safe and healthy: the liberal policy's yield offered to Earth's investors and the
# share of the funders' returns paid in U, searched over a grid, each with growth as a result under the five tests;
# the one that carries the most residents by 2341, the higher yield and share where they tie.
SEARCH_YIELDS = (0.0, 0.01, 0.02, 0.03, 0.04)
SEARCH_FUNDERS_IN_U = (0.0, 0.1, 0.25, 0.5)


def policy_search(scenario='central', tcase='central'):
    grid = []
    for real_yield in SEARCH_YIELDS:
        for in_U in SEARCH_FUNDERS_IN_U:
            opts = tuple(sorted((('u_policy', 'liberal'), ('u_real_yield', real_yield), ('u_funders_paid_in_U', in_U))))
            rates, sched = solve(scenario, tcase, 'steady', opts)
            rows = economy.Run(scenario, tcase, schedule=sched, **dict(opts)).simulate()
            s = economy.summary(rows)
            grid.append(dict(real_yield=real_yield, funders_paid_in_U=in_U, sustained=rates['sustained'], growth_rates=rates,
                             residents_2241=rows[200]['residents'], residents_2341=rows[HORIZON]['residents'],
                             plan_reached_year=s['plan_reached_year'], funding_needed_USD2041=s['funding_needed_USD2041'],
                             NPV3_USD2041=s['NPV_USD2041']['0.03'][str(HORIZON)], u_test_failures=s['u_test_failures'], options=opts))
    carried = [g for g in grid if g['sustained']]
    best = max(carried, key=lambda g: (round(g['residents_2341']), g['real_yield'], g['funders_paid_in_U'])) if carried else None
    return grid, best


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dependencies():
    from research.studies.vacuum_economy import demand
    from . import designs
    files = [HERE / f for f in ('run.py', 'economy.py', 'designs.py', 'engineering.py', 'money.py', 'sales.py', 'technology.py', 'period_technology.json', 'period_inputs.json')]
    files += [designs.TOWNS, designs.CITIES, engineering.REACH, ROOT / 'research/studies/vacuum_economy/scenarios.json', Path(demand.__file__),
              ROOT / 'research/studies/vacuum_economy/programme.json', Path(definition.__file__)]
    return {str(Path(p).resolve().relative_to(ROOT)): digest(p) for p in files}


def product():
    records = [record(name, sc, tc, ramp, opts, *case_rows(sc, tc, ramp, opts)) for name, sc, tc, ramp, opts in cases()]
    grid, best = policy_search()
    if best is not None:
        records.append(record('central:central:steady:U_most_liberal_safe', 'central', 'central', 'steady', best['options'],
                              *case_rows('central', 'central', 'steady', best['options'])))
    projections = [dict(name=name, **exhaustion(sc, tc, path, opts)) for name, sc, tc, path, opts in exhaustion_cases()]
    return dict(
        schema=SCHEMA, producer=dict(path=str((HERE / 'run.py').relative_to(ROOT)), sha256=digest(HERE / 'run.py')),
        dependencies=dependencies(),
        evidence=('A model of the vacuum buildout\'s economy built on stated paths of technology (period_technology.json) and stated inputs '
                  '(period_inputs.json), each with its basis; the settlements\' quantities are their design studies\'. The paths beyond their '
                  'sourced values, the rates of work, the plant masses, the bounds on visits and depot sales, the set-aside shares of the '
                  'sources\' water and the open design choices are the study\'s own. Results hold on those inputs and show what they imply; '
                  'no plant, market, price, site or source is qualified. The time to exhaustion is a rough projection with technology held '
                  'at its 2241 values.'),
        units='Years of the programme (zero is 2041); residents; kilograms and tonnes as named; watts and gigawatts; constant dollars of 2041.',
        reading_rule=('Each case is one scenario of demand, one technology case and one schedule of residents to 2341 (year 300). Money is '
                      'the account with Earth: receipts from Earth\'s customers against what the settlement pays Earth. Funding needed is the '
                      'deepest the cumulative account goes; the opening funding is how deep it goes by year 8; the net present values discount '
                      'each year\'s net to year zero at the rate named, to years 100, 200 and 300; repaid_year is the first year after the '
                      'opening in which the cumulative account is not negative and payback_year the first from which it stays so to 2341. A '
                      'solved schedule gives its growth rates. at_year holds chosen years in full and every_fifth_year a series of years 0, '
                      '5, ..., 300. A case\'s years after its water runs out (exhausted_year) are not physical. Each exhaustion projection '
                      'carries one path on to 3541 at most and names the years the plan is reached, the sites are full and the water open '
                      'to mining runs out. In a case with a policy of U (options u_policy), the funders\' measures come from their own '
                      'flows, which carry Earth\'s purchases of U net of conversions and the dollar reserve, the ledger\'s series are '
                      'given beside the others, and u_test_failures counts the years failing each of the five tests; the policy search '
                      'lists every yield and share of the funders\' returns paid in U tried, and the most liberal policy that stays safe '
                      'and healthy is the case central:central:steady:U_most_liberal_safe.'),
        cases=records, exhaustion=projections, u_policy_search=[{k: x for k, x in g.items() if k != 'options'} for g in grid])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--case', help='scenario:technology:ramp, year by year')
    parser.add_argument('--exhaustion', action='store_true', help='the rough time to exhaustion')
    parser.add_argument('--write', action='store_true', help='write results/periods.json')
    args = parser.parse_args()
    if args.write:
        out = HERE / 'results/periods.json'
        out.parent.mkdir(exist_ok=True)
        out.write_text(json.dumps(product(), indent=1) + '\n')
        print(out)
        return
    if args.exhaustion:
        print('%-50s %6s %6s %9s %9s %9s %10s %8s' % ('projection', 'plan', 'deficit', 'sites', 'water out', 'res. end', 'drawn Gt', 'of Gt'))
        for name, sc, tc, path, opts in exhaustion_cases():
            e = exhaustion(sc, tc, path, opts)
            year = lambda y: '-' if y is None else str(2041 + y)
            print('%-50s %6s %6s %9s %9s %8.1fM %10.1f %8.0f' % (name, year(e['plan_reached_year']), year(e['lasting_deficit_from_year']), year(e['sites_full_year']),
                                                              year(e['water_exhausted_year']), e['residents_at_end'] / 1e6, e['contained_drawn_kg_at_end'] / 1e12,
                                                              e['usable_water_kg'] / 1e12))
        return
    if args.case:
        scenario, tcase, ramp = args.case.split(':')[:3]
        rows, rates = case_rows(scenario, tcase, ramp)
        if rates:
            print('growth rates:', rates)
        print('%5s %11s %9s %9s %9s %7s %13s %7s %8s %7s %6s' % ('year', 'residents', 'receipts', 'to Earth', 'cumul.', 'GW', 'power', 'local', 'tele', 'pads', 'leg'))
        for r in rows:
            if r['year'] % 10 == 0 or r['year'] in (8, HORIZON):
                print('%5d %11.0f %9.3f %9.3f %9.2f %7.1f %13s %7.2f %8.0f %7.0f %6s' % (
                    r['year'], r['residents'], r['receipts'] / 1e12, r['paid_to_Earth'] / 1e12, r['cumulative'] / 1e12, r['firm_GW'],
                    r['power_chosen'] or '-', r['local_share'], r['teleoperated_jobs'], r['pads'], r['leg_m_s']))
        return
    print('%-44s %9s %9s %7s %9s %8s %6s %11s %11s' % ('case', 'fund 2241', 'fund 2341', 'open T', 'NPV3 2341', 'repaid', 'plan', 'res 2241', 'res 2341'))
    year = lambda y: '-' if y is None else str(2041 + y)
    for name, sc, tc, ramp, opts in cases():
        rows, rates = case_rows(sc, tc, ramp, opts)
        s = economy.summary(rows)
        f = s['funding_needed_by_year_USD2041']
        print('%-44s %9.2f %9.2f %7.2f %9.2f %8s %6s %11.0f %11.0f' % (name, f['200'] / 1e12, f[str(HORIZON)] / 1e12, s['opening_funding_USD2041'] / 1e12,
                                                                    s['NPV_USD2041']['0.03'][str(HORIZON)] / 1e12, year(s['repaid_year']), year(s['plan_reached_year']),
                                                                    rows[200]['residents'], rows[HORIZON]['residents']))


if __name__ == '__main__':
    main()
