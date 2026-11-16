"""The cases of the vacuum economy period by period, the growth of residents as a result, and the saved product.

usage (repository root):
  python3 -m research.studies.vacuum_economy_periods.run                      a table of every case
  python3 -m research.studies.vacuum_economy_periods.run --case NAME          one case, year by year
  python3 -m research.studies.vacuum_economy_periods.run --write              results/periods.json

A case is a scenario of demand (low, central, high), a technology case (static, slow, central, fast) and a schedule of
residents: one of the programme's ramps (constant_growth, front_loaded) or residents as a result (steady: the highest
steady growth rate; fastest: the coupled study's rule, the highest steady rate period by period). Residents as a result
grow from the opening's 12,500 in year 8 as fast as the settlement can while its funders still earn the programme's 3%
(the net present value of the account with Earth at 3% to year 220 is not negative), its own work staffed and its
sources' water not exhausted; they are held at the planned 6.47 million once reached. Residents cost the account more
than they bring, since the receipts follow Earth's demand, so this bounds their growth.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from research.studies.vacuum_economy import programme as definition
from . import economy, technology
from .economy import OPENING

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.vacuum-economy-periods/1'
SCENARIOS = ('low', 'central', 'high')
RAMPS = ('constant_growth', 'front_loaded', 'steady')
PERIODS = ((9, 12), (13, 18), (19, 25), (26, 37), (38, 50), (51, 65), (66, 80), (81, 100), (101, 125), (126, 150), (151, 175), (176, 200))
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
}
TOLERANCE = 1e-6


def planned():
    return definition.load()['population']['planned']


def schedule(rates, cap=True):
    """Residents by year 0-220 from period rates: (first year, last year, rate) triples after the opening's 12,500."""
    p = definition.load()['population']
    n, out = p['opening'], [0.0] * OPENING + [p['opening']]
    for y in range(OPENING + 1, 221):
        r = next((rate for a, b, rate in rates if a <= y <= b), 0.0)
        n = n * (1 + r)
        if cap:
            n = min(n, p['planned'])
        out.append(n)
    return out


def sustained(scenario, tcase, sched, options):
    """Whether a schedule is carried: the funders earn the programme's 3% (net present value at 3% to year 220 not
    negative), the settlements' own work is staffed and the sources' water is not exhausted. Returns the test and rows."""
    rows = economy.Run(scenario, tcase, schedule=sched, **dict(options)).simulate()
    s = economy.summary(rows)
    ok = s['NPV_USD2041']['0.03']['220'] >= 0 and s['hands_short_years'] == 0 and s['exhausted_year'] is None
    return ok, rows


def solve(scenario, tcase, rule='steady', options=()):
    """Residents as a result. 'steady': the highest single growth rate from the opening to the plan. 'fastest': period
    by period, the highest steady rate that is sustained to year 220 with no growth after the period. Each rate is found
    by bisection to a hundredth of a percentage point a year."""
    def best(prefix, first, last):
        lo, hi = 0.0, 0.3
        if not sustained(scenario, tcase, schedule(prefix + [(first, last, lo)]), options)[0]:
            return 0.0
        for _ in range(14):
            mid = (lo + hi) / 2
            if sustained(scenario, tcase, schedule(prefix + [(first, last, mid)]), options)[0]:
                lo = mid
            else:
                hi = mid
        return lo
    if rule == 'steady':
        rate = best([], OPENING + 1, 200)
        return [(OPENING + 1, 200, rate)]
    rates = []
    for a, b in PERIODS:
        rates.append((a, b, best(rates, a, b)))
        if schedule(rates)[b] >= planned() * (1 - 1e-9):
            break
    return rates


def case_rows(scenario, tcase, ramp, options=()):
    if ramp in ('steady', 'fastest'):
        rates = solve(scenario, tcase, ramp, options)
        return economy.Run(scenario, tcase, schedule=schedule(rates), **dict(options)).simulate(), rates
    return economy.Run(scenario, tcase, ramp, **dict(options)).simulate(), None


SERIES = ('residents', 'receipts', 'paid_to_Earth', 'net', 'cumulative', 'firm_GW', 'power_chosen', 'local_share', 'imported_t',
          'launches_year', 'pads', 'tug', 'leg_m_s', 'water_delivered_t', 'contained_drawn_t', 'visitors', 'visitor_share_of_world',
          'workforce', 'teleoperated_jobs', 'lines_served_share', 'housed_town', 'housed_city', 'housing_short', 'fare', 'landed_USD_kg', 'kWh_USD')


def record(name, scenario, tcase, ramp, options, rows, rates):
    s = economy.summary(rows)
    at = {}
    for h in (8, 25, 50, 100, 150, 200, 220):
        r = rows[h]
        at[str(h)] = dict({k: r[k] for k in SERIES}, to_Earth=r['to_Earth'], jobs=r['jobs'], power_by_option_kW=r['power_by_option_kW'],
                          lines={k: dict(receipts=x['receipts'], served=x['served']) for k, x in r['lines'].items()},
                          production_kg=r['production_kg'])
    series = {k: [rows[y][k] for y in range(0, 221, 5)] for k in SERIES}
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
    return out


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dependencies():
    from research.studies.vacuum_economy import demand
    from . import designs, engineering
    files = [HERE / f for f in ('run.py', 'economy.py', 'designs.py', 'engineering.py', 'sales.py', 'technology.py', 'period_technology.json', 'period_inputs.json')]
    files += [designs.TOWNS, designs.CITIES, engineering.REACH, ROOT / 'research/studies/vacuum_economy/scenarios.json', Path(demand.__file__),
              ROOT / 'research/studies/vacuum_economy/programme.json', Path(definition.__file__)]
    return {str(Path(p).resolve().relative_to(ROOT)): digest(p) for p in files}


def product():
    records = [record(name, sc, tc, ramp, opts, *case_rows(sc, tc, ramp, opts)) for name, sc, tc, ramp, opts in cases()]
    return dict(
        schema=SCHEMA, producer=dict(path=str((HERE / 'run.py').relative_to(ROOT)), sha256=digest(HERE / 'run.py')),
        dependencies=dependencies(),
        evidence=('A model of the vacuum buildout\'s economy built on stated paths of technology (period_technology.json) and stated inputs '
                  '(period_inputs.json), each with its basis; the settlements\' quantities are their design studies\'. The paths beyond their '
                  'sourced values, the rates of work, the plant masses, the bounds on visits and depot sales and the five open design '
                  'choices are the study\'s own. Results hold on those inputs and show what they imply; no plant, market or price is qualified.'),
        units='Years of the programme (zero is 2041); residents; kilograms and tonnes as named; watts and gigawatts; constant dollars of 2041.',
        reading_rule=('Each case is one scenario of demand, one technology case and one schedule of residents. Money is the account with '
                      'Earth: receipts from Earth\'s customers against what the settlement pays Earth. Funding needed is the deepest the '
                      'cumulative account goes; the opening funding is how deep it goes by year 8; the net present values discount each '
                      'year\'s net to year zero at the rate named; the payback year is the first after the opening in which the cumulative '
                      'account is not negative. A solved ramp gives its growth rates by period. at_year holds chosen years in full and '
                      'every_fifth_year a series of years 0, 5, ..., 220.'),
        cases=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--case', help='scenario:technology:ramp, year by year')
    parser.add_argument('--write', action='store_true', help='write results/periods.json')
    args = parser.parse_args()
    if args.write:
        out = HERE / 'results/periods.json'
        out.parent.mkdir(exist_ok=True)
        out.write_text(json.dumps(product(), indent=1) + '\n')
        print(out)
        return
    if args.case:
        scenario, tcase, ramp = args.case.split(':')[:3]
        rows, rates = case_rows(scenario, tcase, ramp)
        if rates:
            print('growth rates by period:', [(a, b, round(r * 100, 2)) for a, b, r in rates])
        print('%5s %11s %9s %9s %9s %7s %13s %7s %8s %7s %6s' % ('year', 'residents', 'receipts', 'to Earth', 'cumul.', 'GW', 'power', 'local', 'tele', 'pads', 'leg'))
        for r in rows:
            if r['year'] % 10 == 0 or r['year'] in (8, 220):
                print('%5d %11.0f %9.3f %9.3f %9.2f %7.1f %13s %7.2f %8.0f %7.0f %6s' % (
                    r['year'], r['residents'], r['receipts'] / 1e12, r['paid_to_Earth'] / 1e12, r['cumulative'] / 1e12, r['firm_GW'],
                    r['power_chosen'] or '-', r['local_share'], r['teleoperated_jobs'], r['pads'], r['leg_m_s']))
        return
    print('%-34s %8s %8s %8s %7s %9s %9s' % ('case', 'fund T', 'open T', 'NPV3 T', 'payback', 'res 100', 'res 200'))
    for name, sc, tc, ramp, opts in cases():
        rows, rates = case_rows(sc, tc, ramp, opts)
        s = economy.summary(rows)
        print('%-34s %8.2f %8.2f %8.2f %7s %9.0f %9.0f' % (name, s['funding_needed_USD2041'] / 1e12, s['opening_funding_USD2041'] / 1e12,
                                                         s['NPV_USD2041']['0.03']['200'] / 1e12, s['payback_year'], rows[100]['residents'], rows[200]['residents']))


if __name__ == '__main__':
    main()
